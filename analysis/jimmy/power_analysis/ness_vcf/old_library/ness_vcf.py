import vcf, sys, math, re
from collections import defaultdict, Counter
from copy import deepcopy


def consensus_sequence(vcf_file, chromosome, position, purity=0.5, min_quality=0, indels=False, samples='all'):
	#################################################################
	""" takes a vcf_file and position and returns the majority allele """
	#################################################################
	vcf_reader = vcf.Reader(filename=vcf_file)
	for vcf_record in vcf_reader.fetch(chromosome, position-1, position+10):
		if vcf_record.POS == position:
			break
	if vcf_record.POS != position:
		consensus_allele = 'N'
		return consensus_allele
	if indel(vcf_record) or vcf_record.is_snp:
		#there are some decisions to make
		GTs = [i['GT'] for i in vcf_record.samples if i['GQ']>min_quality ]
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

def genotypes_to_AAF(record, min_called=2, ploidy=1, min_GQ=0, min_DP=0, min_MQ=0, max_hets=None, samples='all', vcf_ploidy=2):
	##########################################################################################
	""" takes record and returns alt allele freq and num_called which exceed quality threshold """
	##########################################################################################
	if record.num_called >= min_called and record.is_snp and record.INFO['MQ']>=min_MQ:
		if max_hets== None: max_hets = len(record.samples)
		if samples == 'all' or samples == ['all']:
			GTs = [i['GT']  for i in record.samples if (i['GT'] and i['GQ']>min_GQ and i['DP']>min_DP)]
		else:
			GTs = [record.genotype(i)['GT'] for i in samples if (record.genotype(i)['GT'] and record.genotype(i)['GQ']>min_GQ and record.genotype(i)['DP']>min_DP)]
		if len(GTs) >= min_called:
			if ploidy == 2:
				#this sums non-reference alleles
				x = lambda x: 1 if x>0 else 0
				AAF = sum([x(int(i.split("/")[0]))+x(int(i.split("/")[1])) for i in GTs])/(len(GTs)*float(ploidy))
				return (AAF, 2*len(GTs))
			elif ploidy==1 and vcf_ploidy==1:
				x = lambda x: 1 if x>0 else 0
				AAF = sum([x(int(i)) for i in GTs])/(len(GTs)*float(ploidy))
				return (AAF, len(GTs))
			elif ploidy==1 and vcf_ploidy==2 and record.num_het<=max_hets:
				x = lambda x: 1 if x>0 else 0
				AAF = sum([x(int(i.split("/")[0]))+x(int(i.split("/")[1])) for i in GTs])/(len(GTs)*2.0)
				return (AAF, len(GTs))
			else: return None
		else:
			return None
	#if the site is invariant we only use depth filters
	elif record.num_called >= min_called and record.INFO['MQ']>=min_MQ and not indel(record):
		if samples == 'all' or samples == ['all']: callable = len([i for i in record.samples if (i and i['DP']>min_DP)])
		else: callable = len([record.genotype(i)['DP'] for i in samples if (record.genotype(i)['DP']>min_DP)])
		if callable >= min_called:
			return (0.0, callable*ploidy)
		else: return None
	else: 
		return None

def genotypes_to_MAF(record, min_called=2, ploidy=1, min_GQ=0, min_DP=0, min_MQ=0, max_hets=None, samples='all',vcf_ploidy=2):
	##########################################################################################
	""" takes record and returns minor allele freq and num_called which exceed quality threshold """
	##########################################################################################
	if record.num_called >= min_called and record.is_snp and record.INFO['MQ']>=min_MQ:
		if max_hets ==None: max_hets = len(record.samples)
		if samples == 'all' or samples == ['all']:
			GTs = [i['GT']  for i in record.samples if (i['GT'] and i['GQ']>min_GQ and i['DP']>min_DP)]
		else:
			GTs = [record.genotype(i)['GT'] for i in samples if (record.genotype(i)['GT'] and record.genotype(i)['GQ']>min_GQ and record.genotype(i)['DP']>min_DP)]
		if len(GTs) >= min_called:
			if ploidy == 2: 
				MAF = Counter(re.findall('[0-9]', str(GTs))).most_common(2)[-1][-1]/(2.0*len(GTs))
				return (MAF, 2*len(GTs))
			elif ploidy == 1 and vcf_ploidy==2 and record.num_het<=max_hets:
				MAF = Counter(re.findall('[0-9]', str(GTs))).most_common(2)[-1][-1]/(2.0*len(GTs))
				return (MAF, len(GTs))
			elif ploidy==1 and vcf_ploidy==1:
				MAF = Counter(re.findall('[0-9]', str(GTs))).most_common(2)[-1][-1]/(1.0*len(GTs))
				return (MAF, len(GTs))
			else: return None
		else:
			return None
	#if the site is invariant we only use depth filters
	elif 'MQ' in record.INFO and record.num_called >= min_called and record.INFO['MQ']>=min_MQ and not indel(record):
		if samples == 'all' or samples == ['all']: callable = len([i for i in record.samples if (i and i['DP']>min_DP)])
		else: callable = len([record.genotype(i)['DP'] for i in samples if (record.genotype(i)['DP']>min_DP)])
		if callable >= min_called:
			return (0.0, callable*ploidy)
	elif record.num_called >= min_called and not indel(record):
		if samples == 'all' or samples == ['all']: callable = len([i for i in record.samples if (i and i['DP']>min_DP)])
		else: callable = len([record.genotype(i)['DP'] for i in samples if (record.genotype(i)['DP']>min_DP)])
		if callable >= min_called:
			return (0.0, callable*ploidy)
	else: 
		return None

def indel(record):
	"""take a vcf record and return TRUE if the record is an indel"""
	for i in record.ALT:
		if i:
			if len(i) > 1:
				return True
		else:
			return False
	if len(record.REF) >1:
		return True
	else:
		return False

def one_site_pi(record, min_called=2, ploidy=1, min_GQ=0, min_DP=0, min_MQ=0, max_hets=0, samples='all', vcf_ploidy=2 ):
	"""Calculate theta Pi for one site"""
	#AAF = genotypes_to_AAF(record, min_called, ploidy,min_GQ,min_DP,min_MQ, max_hets, samples, vcf_ploidy)
	MAF = genotypes_to_MAF(record, min_called, ploidy,min_GQ,min_DP,min_MQ, max_hets, samples,vcf_ploidy)
	if MAF!=None:
		num_alleles =MAF[1]
		if num_alleles >= min_called:
			return (float(num_alleles)/(num_alleles - 1.0)) * (2.0 * MAF[0] * (1-MAF[0]))
		else: return None
	else: return None

def pi_from_genotypes(vcf_reader, region, min_called=2, ploidy=1, min_GQ=0, min_DP=0,min_MQ=0, max_hets=0, samples='all', vcf_ploidy=2 ):
	# expects samtools style region and a tabixed vcf "chromosome:start-end"
	# Note the tabixed vcf has to end in .gz rather than .bgz for some buggy reason
	chromosome = str(region.split(":")[0])
	start_coord = int(region.split(":")[1].split("-")[0])
	end_coord = int(region.split(":")[1].split("-")[1])
	pis = [one_site_pi(record, min_called,ploidy, min_GQ, min_DP,min_MQ, max_hets, samples, vcf_ploidy) for record in vcf_reader.fetch(chromosome, start_coord, end_coord) if one_site_pi(record, min_called,ploidy, min_GQ, min_DP,min_MQ, max_hets, samples, vcf_ploidy)!=None]
	if len(pis) > 0:
		return [region, len([i for i in pis if i>0.0]), len(pis), sum(pis)/len(pis)]
	else: 
		return [region, 'NA', 'NA', 'NA']

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

def pi_from_AF(vcf_reader, region, min_called=2, ploidy=1):
	"""  calculates theta pi for a region of vcf """
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

def SFS_from_AF(vcf_reader, region, min_called=2, ploidy=1, max_hets=None):
	chromosome = str(region.split(":")[0])
	start_coord = int(region.split(":")[1].split("-")[0])
	end_coord = int(region.split(":")[1].split("-")[1])
	AFs = [(AF1(record, min_called, ploidy)[0], ploidy*record.num_called) \
			for record in vcf_reader.fetch(chromosome, start_coord, end_coord) \
			if AF1(record,min_called, ploidy)!=None]
	num_alleles = ploidy*len(vcf_reader.samples)
	#initialize empty SFSs
	SFSs={}
	for i in range(min_called, ploidy*len(vcf_reader.samples)+1):
		SFSs[i] = SFS([0]*(i+1)) #zeros 
	for a in AFs:
		num_alleles = a[1]
		freq = a[0]
		SFSs[num_alleles].add(freq, num_alleles)
	return SFSs #dictionary of SFSs for multiple allele depths

def SFS_from_genotypes(vcf_reader, region, min_called=2, ploidy=1, min_GQ=0, min_DP=0, min_MQ=0, max_hets=None, samples='all', vcf_ploidy=2):
	"""this goes to each site in a VCF and calculates the alternate allele frequency then summarises those calculations in SFS objects
	It returns a dictionary of SFSs - where the key is the number of alleles and the value is the SFS of sites with that many alleles"""
	chromosome = str(region.split(":")[0])
	start_coord = int(region.split(":")[1].split("-")[0])
	end_coord = int(region.split(":")[1].split("-")[1])
	AFs = [genotypes_to_AAF(record, min_called, ploidy, min_GQ, min_DP, min_MQ, max_hets, samples, vcf_ploidy) \
	for record in vcf_reader.fetch(chromosome, start_coord, end_coord) \
	if genotypes_to_AAF(record, min_called, ploidy, min_GQ, min_DP, min_MQ, max_hets, samples, vcf_ploidy) !=None]
	SFSs={}
	if samples == ['all'] or samples == 'all':
		for i in range(min_called, ploidy*len(vcf_reader.samples)+1):
			SFSs[i] = SFS([0]*(i+1))
	else:
		for i in range(min_called, ploidy*len(samples)+1):
			SFSs[i] = SFS([0]*(i+1))
	for a in AFs:
		num_alleles = a[1]
		freq = a[0]
		SFSs[num_alleles].add(freq, num_alleles)
	return SFSs #dictionary of SFSs for multiple allele depths

def taj_D(n, k, S):
	# n = num_alleles_sampled
	# k = k ~ pi -avg pairwise diffs
	# S = number of segregating sites
	if float(k)==float(S):
		D=0.0
	elif n < 2:
		return None
	else:
		a1 = sum([1.0/i for i in range(1, n)])
		a2 = sum([1.0/(i**2) for i in range(1, n)])
		b1 = float(n+1)/(3*(n-1))
		b2 = float(2 * ( (n**2) + n + 3 )) / (9*n*(n-1))
		c1 = b1 - 1.0/a1
		c2 = b2 - float(n+2)/(a1 * n) + float(a2)/(a1 ** 2)
		e1 = float(c1) / a1
		e2 = float(c2) / ( (a1**2) + a2 )
		if ((e1 * S )+ ((e2 * S) * (S - 1))) == 0.0:
			return None	
		else: 
			D = (float(k - (float(S)/a1))/math.sqrt((e1 * S )+ ((e2 * S) * (S - 1) )))
	return D

def pi_variance(pi, n, L):
	"""
	This calculates the sampling variance of nucleotide diversity
	following M Nei 1987 "Molecular Evolutionary Genetics" equations 10.7 and 10.9
	pi = nucleotide diversity
	n = number of alleles ie number of individuals * ploidy
	L = length of sequence
	"""
	var_pi = (pi*(n+1)/(3*L*(n-1))) + (pi**2)*(2*(n**2+n+3)/(9*n*(n-1)))
	return var_pi

class SFS:
	def __init__(self, l):
		self.sfs = l #l is a list 
		self.alleles = False
		self.folded = False
	def allele_count(self):
		if self.alleles: return self.alleles
		if self.folded: self.alleles = 2* (len(self.sfs)-1)
		else: self.alleles = len(self.sfs)-1
		return self.alleles
	def invariant(self): return self.sfs[0] + self.sfs[-1]
	def sites(self): return sum(self.sfs)
	def variant(self): return self.sites() - self.invariant()
	def taj_D(self):
		if self.sites() > 0:
			D = taj_D(self.alleles, self.theta_pi()*self.sites(), self.variant())
		else:
			D = None
		return D
	def theta_pi(self):
		if self.sites() > 0:
			pi = sum([(self.sfs[bin] * (2.0 * self.allele_count()/ (self.allele_count()-1.0)) * (float(bin)/self.allele_count())*(1-(float(bin)/self.allele_count()))) for bin in range(len(self.sfs))])
			pi = pi/self.sites()
		else:
			pi= None
		return pi
	def theta_w(self):
		if self.sites() > 0:
			if self.folded ==True:
				a1 = sum([1.0/i for i in range(1,2*self.alleles)])
			else:
				a1 = sum([1.0/i for i in range(1,self.alleles)])
			w = (self.variant()/a1) / float(self.sites())
		else:
			w = None
		return w
 	def add(self, AF, alleles, sites=1):
		if alleles != self.allele_count(): 
			#print "what the shit! You can't add that its the wrong number of alleles"
			pass
		bin = int(round(AF*alleles))
		self.sfs[bin] +=1*sites
		return self
	def fold(self):
		for i in range(len(self.sfs)/2): self.sfs[i], self.sfs[-i-1] = self.sfs[i]+self.sfs[-i-1], 0
		self.folded = True
		return self
		# i = 1
		# while i < (len(self.sfs)+1)/2.0:
		# 	self.sfs[i-1] = self.sfs[i-1]+self.sfs[-1*i]
		# 	self.sfs[-1*i] = 0
		# 	i +=1
	def summarise(self, header=True):
		"""Prints a few summary stats for and SFS object, but returns NOTHING"""
		if header==True:
			print "\t".join(['alleles', 'folded', 'sites', 'invariant', 'variant', 'theta_pi', 'theta_w', 'taj_D', 'sfs'])
		sys.stdout.write("\t".join( [str(i) for i in [self.alleles, self.folded, \
		self.sites(), self.invariant(), self.variant(), \
		self.theta_pi(), self.theta_w(), self.taj_D(), self.sfs]]) + "\n")
