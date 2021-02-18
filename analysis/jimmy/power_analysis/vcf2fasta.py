#!/bin/env python2.7

import vcf, sys, re, argparse
from Bio import SeqIO
from collections import OrderedDict
import ness_vcf


def hap2dip_converter(bases):
	converter = {'AA': 'A', 'TT':'T', 'GG':'G', 'CC':'C', 'AT':'W', 'GC':'S', 'AG':'R', 'CT':'Y', 'GT': 'K', 'AC': 'M', 'AN':'A','CN':'C','GN':'G','TN':'T','TA':'W', 'CG':'S', 'GA':'R', 'TC':'Y', 'TG': 'K', 'CA': 'M', 'NA':'A','NC':'C','NG':'G','NT':'T', 'NN':'N'}
	hap_base = converter[bases]
	return hap_base


def rc(sequence):
	complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'R':'Y', 'Y':'R', 'K':'M', 'M':'K' ,'W':'W', 'S':'S', '?': '?', 'X':'X', '.': '.', '-':'-', 'B':'V', 'D':'H', 'H':'D', 'V':'B'}
	return ''.join([complements[base] for base in sequence][::-1])


def get_genotype_calls(record, samples,recordAttribute_filters, recordINFO_filters, callFormat_filters, callAttribute_filters, site_lenient=True,call_lenient=True, consensus=False, sample_ploidy=2, vcf_ploidy=2, fill_w_reference=False, missing = 'N',  parse_indels = False):
	invariant = False
	variant = False
	chromosome = record.CHROM
	#now if it were a new gene everything is empty otherwise it just carries on
	position = record.POS
	#make allele dict
	if not record.ALT[0]:
		allele_dict = {0:record.REF}
		invariant = True
		#invariant
	elif record.ALT[0]: #variant
		variant = True
		allele_dict = {0:record.REF}
		idx=1
		for alt in record.ALT: 
			allele_dict[idx] = str(record.ALT[idx-1])
			idx+=1
	#now make genotype calls
	site_genotypes = OrderedDict()
	if ness_vcf.test_site(record,recordAttribute_filters, recordINFO_filters, lenient = site_lenient):
		# if consensus == True:
		# 	consensus_allele = consensus_position(record, position, min_GQ=min_GQ, samples=samples)
		for s in samples:
			good_genotype=ness_vcf.test_genotype(record.genotype(s),callFormat_filters, callAttribute_filters, lenient=call_lenient)
			gt = record.genotype(s)['GT']
			if vcf_ploidy == 2:
				if sample_ploidy==2:
					if good_genotype: site_genotypes[s] = (allele_dict[int(gt.split("/")[0])], allele_dict[int(gt.split("/")[1])])
					else:#uncallable
						if fill_w_reference: site_genotypes[s] = (record.REF, record.REF) 
						else: site_genotypes[s] = (missing, missing)
				elif sample_ploidy ==1:
					if good_genotype:
						site_genotypes[s] = (allele_dict[int(gt.split("/")[0])],)
					else: 
						if fill_w_reference: site_genotypes[s] = (record.REF,) 
						else: site_genotypes[s] = (missing,)
				else:
					sys.stderr.write("we don't deal with your weird ploidy")
					sys.exit()
			elif vcf_ploidy==1:
				if good_genotype: site_genotypes[s] = (allele_dict[int(gt.split("/")[0])],)
				else:#uncallable
					if fill_w_reference: site_genotypes[s] = (record.REF,) 
					else: site_genotypes[s] = (missing,)
	else: #this site is low qual - what do we do?
		for s in samples:
			if sample_ploidy==2:
				if fill_w_reference: site_genotypes[s] = (record.REF, record.REF) 
				else: site_genotypes[s] = (missing, missing)
			elif sample_ploidy==1:
				if fill_w_reference: site_genotypes[s] = (record.REF,) 
				else: site_genotypes[s] = (missing,)
	return site_genotypes


def vcf2fasta(vcf_file, chromosome, start, end, reference_fasta, recordAttribute_filters, recordINFO_filters, callFormat_filters, callAttribute_filters, site_lenient=True,call_lenient=True, samples='all', consensus=False, sample_ploidy=2, vcf_ploidy=2, fill_w_reference=False, missing = 'N', parse_indels = False):
	# things to think about
	# if parse_indels and sample_ploidy ==1 and vcf_ploidy==1 and filter_heterozygote_haploids ==False <- this is an impossible set of situations 
	vr = vcf.Reader(filename=vcf_file)
	current_chromosome = chromosome
	skipable = []
	first_position=True 
	ref_dict = SeqIO.to_dict(SeqIO.parse(open(reference_fasta, 'r'), 'fasta'))
	output_seqs = OrderedDict()
	if samples == 'all': samples = vr.samples
	for s in samples:
		output_seqs[s] = list(ref_dict[chromosome].seq[start-1:end])
	#output_seqs['reference'] = str(ref_dict[chromosome].seq[start-1:end])
	if consensus: consensus_seq = list(ref_dict[chromosome].seq[start-1:end])
	previous_position = start-1
	for record in vr.fetch(chromosome, start-1, end):
		position = record.POS
		if position == previous_position +1 and not record.ALT[0] and fill_w_reference:
			pass #if its invariant and fill_w_reference - there is nothing to do so skip the rest
		else: #site could be variant or invariant but only invariant if fill_w_ref  is True
			site_genotypes = get_genotype_calls(record, samples, recordAttribute_filters=recordAttribute_filters, recordINFO_filters=recordINFO_filters, callFormat_filters=callFormat_filters, callAttribute_filters=callAttribute_filters, site_lenient=True,call_lenient=True, consensus=consensus, sample_ploidy=sample_ploidy, vcf_ploidy=vcf_ploidy, fill_w_reference=fill_w_reference, missing = missing,  parse_indels = parse_indels)
			if first_position and position != start: #fill in the 5' end
				if not fill_w_reference:
					#swap in missing symbol for all the missing bases
					for s in output_seqs:
						for base in range(0,position-start):
							output_seqs[s][base] = missing
				else: pass #this is assuming the reference is all parsed into output seqs
			first_position = False
			#### Missing reference data #################################################################################################
			if position > previous_position + 1:
				#fill in with N's, ref or if it was deleted leave it alone
				for s in samples:
					if fill_w_reference: pass
					else:
						for i in range((previous_position-start)+1, (previous_position-start)+(position-previous_position)):
							output_seqs[s][position-start] = missing
				# now if we set previous position to position-1 it will carry on normally
				previous_position = position-1
			#### Normal Position #################################################################################################
			if position == previous_position + 1 and not record.is_indel: #normal_position
				# if record.ALT and fill_w_reference: #
				# 	pass
				if not record.ALT[0] and fill_w_reference: #
					continue
				else:
					for s in samples:
						#if the site is deleted ("-") then skip it
						if sample_ploidy == 1:
							output_seqs[s][position-start] = site_genotypes[s][0]
						elif sample_ploidy == 2:
							#this currently is for making just one diploid sequence ie with ambiguity calls
							output_seqs[s][position-start] = hap2dip_converter("".join(site_genotypes[s]))
							# pseudocode
							# if sample == homozygous:
							# 	output_seqs[s][position-start] = site_genotypes[s][0]
							# elif sample == heterozygous:
							# 	if using 2 seqs:
							# 		output_seqs[s][0][position-start] = site_genotypes[s][0]
							# 		output_seqs[s][1][position-start] = site_genotypes[s][1]
							# 	elif using 1 diploid seq:
							# 		output_seqs[s][position-start] = hap2dip(site_genotypes[s])
			####   INDEL   ##################################################################################################################
			elif position == previous_position or record.is_indel:
				if parse_indels == True:
					#we are in the dreaded INDEL!!!
					if not record.ALT[0]: #invariant????
						pass #forget it and move on!!!
					elif sample_ploidy == 1 and record.num_het == 0:
						#if there are heterozygotes skip site
						#now there should only be homozygotes
						### insertion ###
						if len(record.REF) < max(len(i) for i in record.ALT):
							if len(record.REF) > 1:
								sys.stderr.write("Not equipped to deal with these complex indels" + \
												 "\t".join([str(x) for x in [record.CHROM, record.POS, record.REF, ",".join([str(i) for i in record.ALT])]]) + "\n") 
							else:
								#alter the base at position to match the insertion
								#get len of alleles, pad the non-inserted individuals with a string of ----'s after the allele
								alleles = [i[0] for i in set(site_genotypes.values())]
								idict={}
								for a in alleles:
									idict[a] = str(a) + "-"*(max(len(i) for i in alleles)-len(a))
								for s in samples:
									output_seqs[s][position-start] = idict[site_genotypes[s][0]]
						### deletion ###
						if len(record.REF) > max(len(i) for i in record.ALT):
							if max(len(i) for i in record.ALT) > 1:
								sys.stderr.write("Not equipped to deal with these complex indels" + "\t".join([str(x) for x in [record.CHROM, record.POS, record.REF, ",".join([str(i) for i in record.ALT])]]) + "\n") 
							else:
								# the number of deleted bases needs to be removed from the following positions
								# The current position is left and len(ALT) - len(REF) deletions are made to the following bases
								for s in samples:
									del_length = len(record.REF) - len(site_genotypes[s][0])
									for i in range (1,1+del_length):
										idx = (position-start) + i
										output_seqs[s][idx] = "-"
					elif sample_ploidy ==2:
						sys.stderr.write("I haven't written the code for this")
						sys.exit()
						# pseudocode
						# if its a heterozygote indel we need two sequences to represent it.
						# otherwise a flag needs to go up maybe lower case for hemizygous bases.
				else:pass
			
			####   something is weird ie position isn't indel isn't missing and isn't normal...VCF is broken? ################################
			else:
				pass
				#print "ORPHAN:"
		previous_position = position
	
	return output_seqs



parser = argparse.ArgumentParser(description="Export FASTA format sequences from a VCF", usage="vcf2fasta.py [options] my_vcf.vcf")
parser.add_argument("-v", "--vcf_file",
					required = True,
					type=str,
					help="This is the vcf input file. It must tabix indexed [required]")
parser.add_argument("-r", "--reference",
					required = True,
					type=str,
					help="This is the genome reference file used to make the vcf. Only the chromosome overlapping the region of interest is required [required]")
parser.add_argument("-i", "--regions",
					required = True,
					type=str,
					nargs="+",
					help="This is the regions of the vcf to generate a FASTA for. In SamTools format chromosome:start-end [required]")
parser.add_argument("-s", "--samples",
					required = False,
					default='all',
					type=str,
					nargs="+",
					help="This is the samples of the vcf to generate a FASTA for [all]")
parser.add_argument("--consensus",
					action="store_true",
					required=False,
					help="If invoked a consensus will be generated [False]")
parser.add_argument("--concatenate",
					action="store_true",
					required=False,
					help="If invoked this will patch all the sequences derived from the regions end 2 end - useful for CDS [False]")
parser.add_argument("--reverse_complement",
					action="store_true",
					required=False,
					help="If invoked this will reverse complement the extracted sequence useful for CDS on - strand [False]")
parser.add_argument("--sample_ploidy",
					default=1,
					type=int,
					help="The ploidy of the sample [1]")
parser.add_argument("--vcf_ploidy",
					default=1,
					type=int,
					help="The ploidy assumed when the vcf was called [1]")
parser.add_argument("--fill_w_reference",
					action="store_true",
					required=False,
					help="If invoked missing and low quality sites will be assumed to match the reference [False]")
parser.add_argument("--filter_heterozygote_haploids",
					action="store_false",
					required=False,
					help="If True apparently heterozygous genotypes in haploids will be assumed to be errors and marked as missing or low quality [True]")
parser.add_argument("--filter_heterozygote_haploid_sites",
					action="store_false",
					required=False,
					help="If True if any individual is apparently heterozygous the whole site in all individuals is assumed to be an error and marked as missing or low quality [True]")
parser.add_argument("--missing",
					required = False,
					type=str,
					default = "N",
					help="This character will be used to denote missing or low quality sites [N]")
parser.add_argument("-g", "--min_GQ",
					dest="min_GQ",
					default=None,
					type=int,
					help="This is the minimum GQ score that a variant site can have [None]")
parser.add_argument("-rgq", "--min_RGQ",
					dest="min_RGQ",
					default=None,
					type=int,
					help="This is the minimum RGQ score that an invariant site can have - specific to GenotypeGVCF [None]")
parser.add_argument("--min_QUAL",
					dest="min_QUAL",
					default=None,
					type=int,
					help="This is the minimum QUAL score that a site can have to be considered [None]")
parser.add_argument("--min_MQ",
					dest="min_MQ",
					default=None,
					type=int,
					help="This is the minimum MQ score that a site can have to be considered [None]")
parser.add_argument("--min_DP",
					dest="min_DP",
					default=None,
					type=int,
					help="This is the minimum depth (DP) that a site can have to be considered [None]")
parser.add_argument("--parse_indels",
					action="store_true",
					required=False,
					help="If invoked indels will be called and inserted into sequences. Note that this will mean the annotations of genes may be altered [False]")
parser.add_argument("--site_lenient",
					action="store_false",
					required=False,
					help="If invoked any missing filter for the site will kill the process [True]")
parser.add_argument("--call_lenient",
					action="store_false",
					required=False,
					help="If invoked any missing filter for a genotype call will kill the process [True]")




#############
args = parser.parse_args()
vcf_file = args.vcf_file
reference = args.reference
regions = args.regions
samples = args.samples
consensus = args.consensus
concatenate = args.concatenate
reverse_complement = args.reverse_complement
sample_ploidy = args.sample_ploidy
vcf_ploidy = args.vcf_ploidy
fill_w_reference = args.fill_w_reference
missing = args.missing
parse_indels = args.parse_indels
### These all need to be converted to filters
#recordAttribute_filters,  
filter_heterozygote_haploid_sites = args.filter_heterozygote_haploid_sites
min_QUAL = args.min_QUAL
#recordINFO_filters
min_MQ = args.min_MQ
#callFormat_filters
min_GQ = args.min_GQ
min_DP = args.min_DP
min_RGQ = args.min_RGQ
#callAttribute_filters
filter_heterozygote_haploids = args.filter_heterozygote_haploids
call_lenient=args.call_lenient
site_lenient=args.site_lenient



#construct filters
recordAttribute_filters={}
if min_QUAL != None: recordAttribute_filters['QUAL']=">%i" %(min_QUAL)
if filter_heterozygote_haploid_sites: recordAttribute_filters['num_hets']="==0"

recordINFO_filters={}
if min_MQ != None: recordINFO_filters['MQ']=">%i" %(min_MQ) 

callFormat_filters={}
if min_GQ != None: callFormat_filters['GQ']=">%i" %(min_GQ) 
if min_RGQ != None: callFormat_filters['RGQ']=">%i" %(min_RGQ)
if min_DP != None: callFormat_filters['DP']=">%i" %(min_DP) 

callAttribute_filters={}
if filter_heterozygote_haploids: callAttribute_filters['is_het']="==False"




# common samples
# samples = 'all'
# samples = ['CC3060',  'CC3064',  'CC3065',  'CC3068',  'CC3069',  'CC3071',  'CC3072',  'CC3076',  'CC3078',  'CC3086'] #quebec mt+
# samples = ['CC3059', 'CC3061', 'CC3062', 'CC3063', 'CC3073', 'CC3075', 'CC3079', 'CC3082', 'CC3083', 'CC3084'] #quebec mt-
# samples = ['CC124', 'CC1373', 'CC1690', 'CC1952', 'CC2342', 'CC2343', 'CC2344', 'CC2931', 'CC2932', 'CC2935', 'CC2936', 'CC2937', 'CC2938'] #all wt species wide
# samples = ['CC124','CC1952', 'CC2938','CC2342', 'CC2935','CC2931' ] #wild mt-
# samples = ['CC2343','CC2936', 'CC1373','CC2344', 'CC2937','CC1690', 'CC2932'] #wild mt+
# reference_fasta = '/home2/data/genomes/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fasta'

output_seqs = OrderedDict()
for region in regions:
	chromosome, start, end = region.split(":")[0],  int(re.split(r"[-:]",region)[-2]), int(re.split(r"[-:]",region)[-1])
	seqs = vcf2fasta(vcf_file, chromosome, start, end, reference, \
					recordAttribute_filters, recordINFO_filters, callFormat_filters, callAttribute_filters, \
					site_lenient=site_lenient,call_lenient=call_lenient, \
					samples=samples, consensus=consensus, sample_ploidy=sample_ploidy, vcf_ploidy=vcf_ploidy, \
					fill_w_reference=fill_w_reference, missing = missing, parse_indels = parse_indels)
	if concatenate:
		for s in samples:
			if s not in output_seqs: output_seqs[s]=seqs[s]
			else: output_seqs[s] += seqs[s]
	else:
		for s in seqs:
			if reverse_complement: sys.stdout.write(">%s|%s\n%s\n" %(s,region,rc("".join(seqs[s]))))
			else: sys.stdout.write(">%s|%s\n%s\n" %(s,region,"".join(seqs[s])))



# for s in seqs:
# 	print s,len("".join(seqs[s]))


if concatenate: #### REGION IS WRONG
	coords = []
	for region in regions:
		chromosome, start, end = re.split("[-:]",region)[0], int(re.split("[-:]",region)[1]), int(re.split("[-:]",region)[2]) 
		coords+=[start,end]
	region="%s:%i-%i" %(chromosome,min(coords),max(coords))
	for s in output_seqs:
		if reverse_complement: sys.stdout.write(">%s|%s\n%s\n" %(s,region,rc("".join(output_seqs[s]))))
		else: sys.stdout.write(">%s|%s\n%s\n" %(s,region,"".join(output_seqs[s])))




"""
eg deletion (quebec_wt.vcf.gz)
19083 A [None]
19083 ACATCTACGAGCC [A]
1/1, None, None, None, None, None, None, None, None, 1/1, None, 0/1, 1/1, None, 0/0, 0/0, 0/0, 0/0, 0/0, 0/0, 


eg insertion (quebec_wt.vcf.gz)
get_genotype_calls(insertion, samples,min_GQ=0,min_DP=0,min_QUAL=0,min_MQ=0, consensus=False, sample_ploidy=1, vcf_ploidy=2, fill_w_reference=False, missing = 'N', filter_heterozygote_haploids =True, parse_indels = True)

19097 A [None]
19097 A [AGGGTATGTCGGAACGCAAGCGATC, AAAAGGCAA]
1/2, None, 2/2, 2/2, 2/2, 2/2, 2/2, 2/2, 2/2, 1/1, None, 2/2, 2/2, 2/2, 0/0, 2/2, 0/0, 0/0, 0/0, 2/2, 

>>> alleles = set(site_genotypes.values())
>>> idict={}
>>> for a in alleles:
...     idict[a] = str(a) + "-"*(max(len(i) for i in alleles)-len(a))
... 
>>> 
>>> for i in idict:
...     print idict[i], i
... 
A------------------------ A
AAAAGGCAA---------------- AAAAGGCAA
N------------------------ N
AGGGTATGTCGGAACGCAAGCGATC AGGGTATGTCGGAACGCAAGCGATC

# 1234567
# 0123456
# 4567890
# ABCDEFG
# start = 4
# pp =6
# position =9
# missing
# from pp 6 
# 7,8 are missing
# pp-start  to (pp-start)+(pos-pp)
# 6-4 to (6-4)+(9-6)
# 2:5
# deletion
position = 6 
7 and 8 are deleted
so REF = 678 ALT = 6
for i in range (1,1+del_length)
	idx = (position-start) + i
"""