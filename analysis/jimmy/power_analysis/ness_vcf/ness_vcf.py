import vcf, sys, re
from collections import Counter

from SFS import SFS


def test_site(pyVCF_RecordObject, recordAttribute_filters, recordInfo_filters, lenient=True, verbose = False):
	"""
	this will filter based on whole site categories
	(1)Attributes of the record that can be accessed like this "record.ATTRIBUTE"
		e.g. QUAL, POS, FILTER, REF, ALT,  num_het, num_called, is_indel, is_sv etc
	(2)INFO field filters accessed like this record.INFO[field[
		e.g. MQ, DP, AF, AN etc
	#############
	THE challenge with whole SITE filters is that they won't necessarily all be numbers:
		e.g. we may want a site with QUAL > 30 and MQ > 20 BUT with ReF
	We need the filter dict to include the name of the thing plus the test
	then use eval to do the test:
		eg. f = 'QUAL': '> 30.0'
			if not eval('pyVCF_RecordObject.%s %s' %(f, SiteProperty_filter_dict[f]))
			#which is thq equivalent of eval(record.QUAL > 30)
	# recordAttribute_filters={"num_called": ">8", "num_hets": "==0"}
	# recordINFO_filters={"MQ": ">30"}
	# callFormat_filters={"GQ": ">30", "RGQ": ">30", "DP":">5"}
	# callAttribute_filters={}
	# SFSs = SFS_from_genotypes(vr, region, recordAttribute_filters=recordAttribute_filters, recordINFO_filters=recordINFO_filters, callFormat_filters=callFormat_filters, callAttribute_filters=callAttribute_filters,  ploidy=1, vcf_ploidy=1)
	########
	""" 
	for f in recordAttribute_filters:
		try:
			if not eval('pyVCF_RecordObject.%s %s' %(f, recordAttribute_filters[f])): #this means that along the way if any of the filters fail its over for this call
				return False
		except AttributeError:
			if verbose: sys.stderr.write("%s does not exist in record at site %s:%i \n" %(f, pyVCF_RecordObject.CHROM, pyVCF_RecordObject.POS))
			if lenient: continue
			else: return False
	for f in recordInfo_filters:
		try:
			if not eval("pyVCF_RecordObject.INFO['%s'] %s" %(f, recordInfo_filters[f])): #this means that along the way if any of the filters fail its over for this call
				return False
		except KeyError:
			if verbose: sys.stderr.write("%s does not exist in the INFO field of site %s:%i \n" %(f, pyVCF_RecordObject.CHROM, pyVCF_RecordObject.POS))
			if lenient: continue
			else: return False
	#if a site got here without returning false its a keeper
	return True

			
def test_genotype(pyVCF_CallObject, callFormat_filters, callAttribute_filters, lenient=True, verbose =False):
	"""
	we want this function to take a PyVCF call object and test whether it passes all the defined filters
	if that field doesn't exist in this call object either return fail or go onto the other filters
	filters = {'min_GQ': min_GQ, 'min_DP': min_DP, 'min_RGQ': min_RGQ, 'min_purity': min_purity, 'called': called, ' heterozygote': heterozygote}
	# ???how do we translate the filters you want to a series of tests ???
		so things in the genotype field will be stored in a dictionary
			eg call['GT'] = 0, call['RGQ'] = 99 #need to be sure these are all minimum values
			these things can be accessed using the filter name to get the value
		other filters that are tests can be accessed from pyVCFs built in convenience functions
			eg. call.is_snp, call.called - NEED TO BE CAREFUL to ensure these are all True/FALSE
			to access these we can use the filter name eg 'called'  or 'is_het' 
			along with the 'getattr(c, filter)' to return the value
	perhaps the best way to set this up is to pass two dicts:
		(1) GT call filters ie anything that might show up in the record field
		(2) pyVCF attributes 
	"""
	#####################
	if not pyVCF_CallObject.called:
		return False
	for f in callFormat_filters:
		try:
			if not eval("pyVCF_CallObject['%s'] %s" %(f, callFormat_filters[f])): #this means that along the way if any of the filters fail its over for this call
			   return False
		except (AttributeError, TypeError): #this catches cases where the field doesn't exist in that call - 
			if verbose: sys.stderr.write("%s does not exist in genotype %s at site %s:%i \n" %(f, pyVCF_CallObject.sample,pyVCF_CallObject.site.CHROM, pyVCF_CallObject.site.POS))
			if lenient : continue
			else: return False
	for f in callAttribute_filters:
		try:
			if not eval('pyVCF_CallObject.%s %s' %(f, callAttribute_filters[f])): #this means that along the way if any of the filters fail its over for this call
				return False
		except (AttributeError, TypeError): #this catches cases where the field doesn't exist in that call - 
			if verbose: sys.stderr.write("%s is not an attribute of genotype %s at site %s:%i \n" %(f, pyVCF_CallObject.sample,pyVCF_CallObject.site.CHROM, pyVCF_CallObject.site.POS))
			if lenient : continue
			else: return False  
	# if you've reached here without returning False the site must have passed
	return True

def genotypes_to_MAF(record, recordAttribute_filters=None, recordINFO_filters=None, callFormat_filters=None, callAttribute_filters=None, samples='all', site_lenient=True, call_lenient=True, ploidy=2, vcf_ploidy=2, minHQ_GTs =2 ):
	##########################################################################################
	""" takes record and returns minor allele freq and num_called which exceed quality threshold """
	##########################################################################################
	if record.is_monomorphic and test_site(record,recordAttribute_filters, recordINFO_filters, lenient = site_lenient):
		if samples == 'all' or samples == ['all']:
			GTs = [i['GT'] for i in record.samples if test_genotype(i,callFormat_filters, callAttribute_filters, lenient=call_lenient)]
		else:
			GTs = [record.genotype(i)['GT'] for i in samples if test_genotype(i,callFormat_filters, callAttribute_filters, lenient=call_lenient)]
		if len(GTs) > minHQ_GTs:
			return (0.0, len(GTs)*ploidy)
		else: return None
	elif test_site(record,recordAttribute_filters, recordINFO_filters, lenient = site_lenient) and not record.is_monomorphic:
		if samples == 'all' or samples == ['all']:
			GTs = [i['GT'] for i in record.samples if test_genotype(i,callFormat_filters, callAttribute_filters, lenient=call_lenient)]
		else:
			GTs = [record.genotype(i)['GT'] for i in samples if test_genotype(i,callFormat_filters, callAttribute_filters, lenient=call_lenient)]
		# now make some calculations
		try:
			if len(GTs) < minHQ_GTs: return None
			elif ploidy == 2: 
				MAF = Counter(re.findall('[0-9]', str(GTs))).most_common(2)[-1][-1]/(2.0*len(GTs))
				return (MAF, 2*len(GTs))
			elif ploidy == 1 and vcf_ploidy==2:
				MAF = Counter(re.findall('[0-9]', str(GTs))).most_common(2)[-1][-1]/(2.0*len(GTs))
				return (MAF, len(GTs))
			elif ploidy==1 and vcf_ploidy==1:
				MAF = Counter(re.findall('[0-9]', str(GTs))).most_common(2)[-1][-1]/(1.0*len(GTs))
				return (MAF, len(GTs))
			else: return None
		except ZeroDivisionError:
			return None
	else:
		return None


def genotypes_to_allele_counts(record, recordAttribute_filters=None, recordINFO_filters=None, callFormat_filters=None, callAttribute_filters=None, samples='all', site_lenient=True, call_lenient=True, ploidy=2, vcf_ploidy=2, minHQ_GTs =2 ):
	##########################################################################################
	""" takes record and returns the top two allele freqs and the alleles thtat pass quality thresholds """
	##########################################################################################
	#this is monomorphic sites that pass filter
	if record.is_monomorphic and test_site(record,recordAttribute_filters, recordINFO_filters, lenient = site_lenient):
		if samples == 'all' or samples == ['all']:
			GTs = [i['GT'] for i in record.samples if test_genotype(i,callFormat_filters, callAttribute_filters, lenient=call_lenient)]
		else:
			GTs = [record.genotype(i)['GT'] for i in samples if test_genotype(record.genotype(i), callFormat_filters, callAttribute_filters, lenient=call_lenient)]
		if len(GTs) > minHQ_GTs:
			return ([str(record.REF), len(GTs)*ploidy],  )
		else: return None
	elif test_site(record,recordAttribute_filters, recordINFO_filters, lenient = site_lenient) and not record.is_monomorphic:
		if samples == 'all' or samples == ['all']:
			GTs = [i['GT'] for i in record.samples if test_genotype(i,callFormat_filters, callAttribute_filters, lenient=call_lenient)]
		else:
			GTs = [record.genotype(i)['GT'] for i in samples if test_genotype(record.genotype(i),callFormat_filters, callAttribute_filters, lenient=call_lenient)]
		# now make some calculations
		try:
			allele_counts = Counter(re.findall('[0-9]', str(GTs))).most_common()
			alleles = [record.REF]+[str(i) for i in record.ALT]
			if len(GTs) < minHQ_GTs:
				return None
			ploidy_correction = vcf_ploidy/float(ploidy) #if you have a diploid VCF and a haploid organism this corrects the real number of alleles you saw.
			return [(alleles[int(i[0])], int(i[-1])/ploidy_correction)  for i in allele_counts ] # this divides the number of alleles by two because each haploid  
		except ZeroDivisionError:
			return None
	else:
		return None


def SFS_from_genotypes(vcf_reader, region, recordAttribute_filters=None, recordINFO_filters=None, callFormat_filters=None, callAttribute_filters=None, samples='all', site_lenient=True, call_lenient=True, ploidy=2, vcf_ploidy=2,minHQ_GTs=2):
	"""this goes to each site in a VCF and calculates the alternate allele frequency then summarises those calculations in SFS objects
	It returns a dictionary of SFSs - where the key is the number of alleles and the value is the SFS of sites with that many alleles"""
	chromosome = str(region.split(":")[0])
	start_coord = int(region.split(":")[1].split("-")[0])
	end_coord = int(region.split(":")[1].split("-")[1])
	AFs =[genotypes_to_MAF(record, recordAttribute_filters=recordAttribute_filters, recordINFO_filters=recordINFO_filters, callFormat_filters=callFormat_filters, callAttribute_filters=callAttribute_filters, samples=samples, site_lenient=site_lenient, call_lenient=call_lenient, ploidy=ploidy, vcf_ploidy=vcf_ploidy,minHQ_GTs=minHQ_GTs) \
			for record in vcf_reader.fetch(chromosome, start_coord, end_coord) \
			if genotypes_to_MAF(record, recordAttribute_filters=recordAttribute_filters, recordINFO_filters=recordINFO_filters, callFormat_filters=callFormat_filters, callAttribute_filters=callAttribute_filters, samples=samples, site_lenient=site_lenient, call_lenient=call_lenient, ploidy=ploidy, vcf_ploidy=vcf_ploidy,minHQ_GTs=minHQ_GTs) !=None]
	SFSs={}
	if samples == ['all'] or samples == 'all':
		for i in range(minHQ_GTs,ploidy*len(vcf_reader.samples)+1):
			SFSs[i] = SFS([0]*(i+1))
	else:
		for i in range(minHQ_GTs, ploidy*len(samples)+1):
			SFSs[i] = SFS([0]*(i+1))
	for a in AFs:
		num_alleles = a[1]
		freq = a[0]
		SFSs[num_alleles].add(freq, num_alleles)
	return SFSs #dictionary of SFSs for multiple allele depths


def AF1(record, min_called=2, ploidy=1, max_hets=None):
	""" for a record this will return an AF1 even if its invariant"""
	if max_hets ==None: max_hets = record.samples
	if record.num_called >= min_called and record.is_snp:
		return (record.INFO['AF'][0], record.num_called*ploidy)
	elif record.num_called >= min_called:
		return (0.0, record.num_called*ploidy)
	else:
		return None

def AF2pi(AF, num_alleles):
	"""pi for a single AF"""
	if AF!=None:
		pi = (float(num_alleles)/(num_alleles - 1.0)) * (2.0 * AF * (1-AF))
		return pi
	else:
		return None


def consensus_sequence(vcf_file, chromosome, position, callFormat_filters={}, callAttribute_filters={}, purity=0.5, indels=False, samples='all'):
	#################################################################
	""" takes a vcf_file and position and returns the majority allele """
	#################################################################
	vcf_reader = vcf.Reader(filename=vcf_file)
	vcf_record = None
	for vcf_record in vcf_reader.fetch(chromosome, position-1, position+10):
		if vcf_record.POS == position:
			break
		else: vcf_record = None
	if vcf_record ==None or vcf_record.POS != position:
		consensus_allele = 'N'
		return consensus_allele
	if samples == 'all': samples = vcf_record.samples
	if vcf_record.is_indel or vcf_record.is_snp:
		#there are some decisions to make
		GTs = [i['GT'] for i in samples if test_genotype(i, callFormat_filters, callAttribute_filters)]
		alleles = re.findall('[0-9]',"".join(GTs)) #so this is a bunch of numbers 
		majority = Counter(re.findall('[0-9]',"".join(GTs))).most_common(1)[0] #this will be tuple (allele[str], count[int]) most common number
		#print vcf_record.POS, Counter(re.findall('[0-9]',"".join(GTs)))
		if majority[1]/float(len(alleles)) > purity: #if its sufficiently common
			maj_allele = majority[0]
			if maj_allele == '0':
				consensus_allele = str(vcf_record.REF)
				if indels == False:
					consensus_allele = str(vcf_record.REF)[0]
			else:
				consensus_allele=str(vcf_record.ALT[int(maj_allele)-1])
				if indels == False and len(consensus_allele) > 1:
					consensus_allele = str(vcf_record.REF)[0]
		else:
			consensus_allele=str(vcf_record.REF)
			if indels == False:
				consensus_allele = str(vcf_record.REF)[0]
	else:
		consensus_allele=str(vcf_record.REF)
		if indels == False:
			consensus_allele = str(vcf_record.REF)[0]
	return consensus_allele 
	#it must be the minor allele unless there are two genotypes (ie 50:50) or many alt alleles


def pi_from_AF(vcf_reader, region, min_called=2, ploidy=1):
	"""  calculates theta pi for a region of vcf assuming AF is defined in the INFO field"""
	chromosome = str(region.split(":")[0])
	start_coord = int(region.split(":")[1].split("-")[0])
	end_coord = int(region.split(":")[1].split("-")[1])
	num_alleles = ploidy*len(vcf_reader.samples)
	pis =  [AF2pi(AF1(i, min_called)[0],num_alleles) for i in vcf_reader.fetch(chromosome, start_coord, end_coord) if AF1(i, min_called)!=None]
#	return [region, len([i for i in pis if i>0]), len(pis), sum(pis)/len(pis)]
	if len(pis) > 0:
		return [region, len([i for i in pis if i>0.0]), len(pis), sum(pis)/len(pis)]
	else: 
		return [region, 'n/a', 'n/a', 'n/a']





