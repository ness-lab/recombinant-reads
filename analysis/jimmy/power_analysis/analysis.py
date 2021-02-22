#!/usr/bin/env python3

import numpy.random as r
import re
import argparse
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def arg_parser():
    parser = argparse.ArgumentParser(
        description='generate sam reads for power analysis of readcomb',
        usage='python3 analysis.py [options]'
    )

    parser.add_argument('-f', '--fasta', required=True,
                        type=str, help='vcf2fasta output fasta file with 2 sequences each representing a haplotype')
    
    parser.add_argument('-s', '--start', required=False, default=1,
                        type=int, help='1-based index start of read generation (default is 1)')

    parser.add_argument('-e', '--end', required=False, default=0,
                        type=int, help='1-based index end of read generation (default is 0/all of chromosome)')

    parser.add_argument('-p', '--percent', required=False, default=0.004,
                        type=int, help='probability of phase change per base (default 0.004)')

    args = parser.parse_args()

    return {
        'fasta': args.fasta,
        'start': args.start,
        'end': args.end,
        'percent': args.percent
    }

def main():

    args = arg_parser()
    
    # parse results from vcf2fasta
    fasta_file = list(SeqIO.parse(args['fasta'], 'fasta'))

    # if end, change to all of chromosome
    if args['end'] == 0:
        args['end'] = len(fasta_file[0].seq)

    events = r.poisson((args['end'] - args['start']) * args['percent'])
    phase_changes = r.choice((args['end']- args['start']), size=events)

    # regex to identify correct fasta reference_name
    chrom = re.search(r'chromosome_[0-9]{1,2}', str(fasta_file[0].id)).group()

    haplotype0 = str(fasta_file[0].seq)
    haplotype1 = str(fasta_file[1].seq)

    current = 0
    idx = args['start']
    recomb_seq = []

    for base in zip(haplotype0, haplotype1):

        # switch haplotype if idx in list phase_changes    
        current = (current + 1) % 2 if idx in phase_changes else current
        recomb_seq.append(base[current])

        if idx >= args['end']:
            break
        else:
            idx += 1

    recomb_record = SeqRecord(
        Seq(''.join(recomb_seq)),
        id=chrom,
        name=chrom,
        description=chrom,
        dbxrefs=fasta_file[0].dbxrefs,
        features=fasta_file[0].features,
        annotations=fasta_file[0].annotations,
        letter_annotations=fasta_file[0].letter_annotations
    )

    SeqIO.write(recomb_record, 'recomb_fasta', 'fasta')

    # art_illumina command
    subprocess.run(['art_illumina',
                    '-f', '400',
                    '-l', '250',
                    '-ss', 'Msv1',
                    '-i', 'recomb_fasta',
                    '-m', '600',
                    '-s', '1',
                    '-o', 'generated_reads',
                    '-M', '-p', '-sam', '-na'])

    
if __name__ == '__main__':
    main()