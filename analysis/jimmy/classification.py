#!/usr/bin/env python3

import os
import re
import sys
import pysam
from cyvcf2 import VCF
from tqdm import tqdm
    
def check_snps(vcf_file_obj, chromosome, left_bound, right_bound):
    """
    Generate all SNPs on given chromosome and VCF file within the
    left_bound and right_bound using cyvcf2

    SNPs must have all calls with GQ >= 30 and no heterozygous calls.
    Please use the filter script to preprocess parental VCF prior to
    phase change detection.

    vcf_file_obj: cyvcf2 VCF file object
    chromosome: string: chromosome name
    left_bound/right_bound: integer: 0-based index of reference sequence

    Return list of cyvcf Variant objects
    """

    # 1 is added to record.reference_start and the following parameter because vcf is 1 indexed
    # in order to keep code consistent
    region = '{c}:{l}-{r}'.format(c=chromosome, l=left_bound+1, r=right_bound+1)

    # list comp with if statement to only include SNPs
    with SilentVCF():
        records = [rec for rec in vcf_file_obj(region)]
    return records

class SilentVCF:
    def __enter__(self):
        self._original_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stderr.close()
        sys.stderr = self._original_stderr

def cigar(record):
    """
    Build the query segment using the cigar tuple given by the bam record so that it aligns with
    indexes of SNPs in the VCF and the reference sequence

    record: bam record from pysam alignment file

    Returns a dna sequence string
    """
    cigar_tuples = record.cigartuples

    # initialize segment for building
    segment = []

    # index to keep track of where we are in the query_segment
    query_segment = record.query_sequence
    index = 0

    # no alignment, cigar_tuple is equal to None
    if cigar_tuples is None:
        return segment

    for cigar_tuple in cigar_tuples:

        # 0 is a match, add it onto the segment
        if cigar_tuple[0] == 0:
            segment.append(query_segment[index:index+cigar_tuple[1]])
            index += cigar_tuple[1]

        # 1 is an insertion to query segment
        # skip it because SNPs are aligned to reference and do not exist in this region
        elif cigar_tuple[0] == 1:
            index += cigar_tuple[1]

        # 2 is an deletion to query segment, add a gap to realign it to reference
        # don't add index to not moving further through the query segment
        elif cigar_tuple[0] == 2:
            segment.append('-' * cigar_tuple[1])

        # 4 is soft clipping, query sequence includes this but is not aligned to reference
        elif cigar_tuple[0] == 4:
            index += cigar_tuple[1]

        # 5 = hard clipping, record.query_sequence does not have the portion that is
        # hard clipping so skip it and don't add anything
        elif cigar_tuple[0] == 5:
            continue

        else:
            raise ValueError(
                'No condition for tuple {} of {}'.format(
                    cigar_tuples.index(cigar_tuple), cigar_tuples))

    return ''.join(segment)

def downstream_phase_detection(variants, segment, record):
    """
    Takes a segment and a list in the region of the segment and generates a list of
    strings to represent in order the parent of each variant on the segment

    snps: list of cyvcf2 variants in the region of the segment given by the function check_snps(),
        can include variants that are outside of the region when there is soft clipping
    segment: string sequence that is the segment built by the function cigar()
    record: bam record from pysam alignment file

    Returns a list of strings ('1', '2', or 'N')
    """

    detection_results = []

    cigar_tuples = record.cigartuples

    # check for no alignment
    if not cigar_tuples:
        return detection_results

    for variant in variants:
        # Using SNP.start and record.reference_start since they are both 0 based
        # SNP.start grabs vcf positions in 0 index while vcfs are 1 indexed
        idx = variant.start - record.reference_start

        parent1 = variant.gt_bases[0].split('/')[0]
        parent2 = variant.gt_bases[1].split('/')[0]

        if parent1 == parent2: # ignore uninformative SNPs
            continue
        
        # ignore if snp is before sequence
        if idx < 0:
            continue
        
        if variant.is_indel:
            # check if indel is outside of segment
            if idx + max(len(parent1), len(parent2)) > len(segment):
                continue

            parent1_match = segment[idx:idx + len(parent1)] == parent1
            parent2_match = segment[idx:idx + len(parent2)] == parent2

            if len(parent1) > len(parent2):
                if parent1_match:
                    detection_results.append(('1', idx))
                elif parent2_match:
                    detection_results.append(('2', idx))
                else:
                    detection_results.append(('N', idx))
            else:
                if parent2_match:
                    detection_results.append(('2', idx))
                elif parent1_match:
                    detection_results.append(('1', idx))
                else:
                    detection_results.append(('N', idx))
        else:
            if idx >= len(segment):
                break

            if segment[idx] == parent1:
                detection_results.append(('1', idx))

            elif segment[idx] == parent2:
                detection_results.append(('2', idx))

            else:
                detection_results.append(('N', idx))
    
    return detection_results

class Pair():
    '''
    bam read class
    '''
    def __init__(self, record1, record2, vcf_filepath):

        # put pairs in order
        if record1.reference_start < record2.reference_start:
            self.rec_1 = record1
            self.rec_2 = record2
        else:
            self.rec_1 = record2
            self.rec_2 = record1
        
        self.vcf_filepath = vcf_filepath

    def call(self, vcf=None):
        '''
        Takes in optional vcf arguement to take in a VCF object for use in multiprocessing/threading
        - Do not use the same VCF object across threads/processes
        '''
        self.segment_1 = cigar(self.rec_1)
        self.segment_2 = cigar(self.rec_2)

        if not vcf:
            vcf = VCF(self.vcf_filepath)

        self.variants_1 = check_snps(vcf, self.rec_1.reference_name, 
                                self.rec_1.reference_start,
                                self.rec_1.reference_start + self.rec_1.query_alignment_length)
        self.variants_2 = check_snps(vcf, self.rec_2.reference_name, 
                                self.rec_2.reference_start,
                                self.rec_2.reference_start + self.rec_2.query_alignment_length)

        self.detection_1 = downstream_phase_detection(self.variants_1, self.segment_1, self.rec_1)
        self.detection_2 = downstream_phase_detection(self.variants_2, self.segment_2, self.rec_2)

        # simplification of results
        prev_strand = None
        prev_location = None
        begin_idx = None
        self.simplify = []
        for variant in self.detection_1 + self.detection_2:
            strand = variant[0]
            location = variant[1]
            # first variant
            if begin_idx == None:
                begin_idx = location
                prev_strand = strand
                prev_location = location
            # strand change
            if strand != prev_strand:
                self.simplify.append((prev_strand, begin_idx, location))
                begin_idx = location

            prev_strand = strand
            prev_location = location
        
        # last variant
        self.simplify.append((prev_strand, begin_idx, location))
        
        # classification
        strands = [tupl[0] for tupl in self.simplify]
        ranges = [(tupl[1], tupl[2]) for tupl in self.simplify]

        cross_over = r'(12)|(21)'
        gene_conversion = r'(12+1)|(21+2)'

        if 'N' in strands:
            self.classify = 'no_match'
        elif re.match(cross_over, ''.join(strands)):
            self.classify = 'cross_over'
        elif re.match(gene_conversion, ''.join(strands)):
            self.classify = 'gene_conversion'
        else:
            self.classify = 'no_classification'



def pairs_creation(bam_filepath, vcf_filepath):

    bam = pysam.AlignmentFile(bam_filepath, 'r')

    pairs = []
    prev_rec = None
    for rec in bam:
        # first record
        if not prev_rec:
            prev_rec = rec
        elif rec.query_name == prev_rec.query_name:
            pairs.append(Pair(rec, prev_rec, vcf_filepath))
            prev_rec = None

    return pairs

if __name__ == '__main__':
    pairs = pairs_creation('readcomb/recomb_diagnosis.sam', 'tests/parental_filtered.vcf.gz')
    
    for pair in pairs[:10]:
        pair.call()
        print(pair.simplify)
        print(pair.classify)

    
