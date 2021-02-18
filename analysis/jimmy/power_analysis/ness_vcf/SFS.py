import sys, math


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
		###
		def td(n, k, S):
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
		####
		if self.sites() > 0:
			D = td(self.alleles, self.theta_pi()*self.sites(), self.variant())
		else:
			D = None
		return D
	def theta_pi(self):
		"""
		f = bin = frequency count ie 1:singleton, 2:doubleton etc #
		x = self.sfs[bin] #number of sites with count bin
		n = self.allele_count() # how many samples you have
		p = float(bin)/self.allele_count()
		pi = sum([(x * (2n/(n-1)) * p*q) for bin in range(len(self.sfs))])
		"""
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
			print("what the shit! You can't add that its the wrong number of alleles")
			pass
		bin = int(round(AF*alleles))
		self.sfs[bin] +=1*sites
		return self
	def fold(self):
		for i in range(int(len(self.sfs)/2)): self.sfs[i], self.sfs[-i-1] = self.sfs[i]+self.sfs[-i-1], 0
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
			print("\t".join(['alleles', 'folded', 'sites', 'invariant', 'variant', 'theta_pi', 'theta_w', 'taj_D', 'sfs']))
		sys.stdout.write("\t".join( [str(i) for i in [self.alleles, self.folded, \
		self.sites(), self.invariant(), self.variant(), \
		self.theta_pi(), self.theta_w(), self.taj_D(), self.sfs]]) + "\n")