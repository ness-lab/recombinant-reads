import gzip

class VCF_ReferenceFileWalker:
	"returns an iterator where each iteration is a group of lines \
	 from a number of sorted files that have an 'index' that matches \
	 the refernce file"
	def __init__(self, ref_file, files, index_position=1, delimiter='\t'): #initialize the instance
		self.ref_file = ref_file #this is the file that we are trying to get matching lines for (ie no group of lines returned can have an empty ref_line)
		self.files = files #this is the group of files that we are matching to the reference file
		self.index_position = index_position # the position of the column that we are using to match lines up
		self.delimiter = '\t' #the delimiter between columns
	def __iter__(self):  #this section basically creates the iterator class but doesn't really doing the looping. It simply sets up the object to make it ready to loop
		if self.ref_file[-2:] in ["gz", "GZ"]: self.ref = gzip.open(self.ref_file, 'rb') #reference file handle
		else: self.ref = open(self.ref_file, 'r') #reference file handle		
		self.fs = [(lambda f: gzip.open(f, 'rb') if f[-2:] in ["gz", "GZ"] else open(f, 'r'))(f) for f in self.files]
		self.ref_line = self.ref.readline() # a reference line
		while self.ref_line[0]=="#": self.ref_line = self.ref.readline() #advance through the ### header lines of the VCF
		self.lines = []# inititiate in  __init__ 
		for i in range(len(self.fs)):    # this advances all the other files through their headerlines
			f_line = self.fs[i].readline()
			while f_line[0] == "#":
				f_line = self.fs[i].readline()
			self.lines.append(f_line) #finish with lines holding all the non-reference first lines
		self.good_lines = ["." for i in range(len(self.fs))] #initialize good_lines variable like this ["." , "." , "."    ]
		self.current_ch = [self.ref_line.split(self.delimiter)[0]] #tells us the chromosome we are currently on
		self.previous_index = 0 # this is the index of the previous iteration (ie the previous genomic position) initialized to zero.
		return self
	def next(self): #this is the real meat of the iterator - this is where the looping comes in and it keeps going until the reference file is finished reading
		self.ref_line = self.ref.readline() #advances the reference file one line for the next iteration
		if self.ref.closed: #this is how you stop the iteration when the reference file is done being read
			raise StopIteration
		if not self.ref_line:  #not sure why I have this twice but one of these two must work
			self.ref.close()
			for f in self.fs: f.close()
			raise StopIteration
		ref_index = int(self.ref_line.split(self.delimiter)[self.index_position]) #defines the index we are on (the genome position of the reference)
		ref_ch = self.ref_line.split(self.delimiter)[0] #defines the chromosome the reference is on (which may or may not match )
		# this section will advance the other file lines to the point where their chromosome matches the chromosome of the reference
		# if they advance so far that their chromosome has not occureed in the reference
		# ie they have no information for the current reference chromosome they will stop there and wait until the reference catches up
		# if however the reference is missing the information so that it is missing a chromosome - this may still be a problem
		# if the sort order is different  - you're fucked
		if ref_ch != self.current_ch[-1]: #we have a new reference ch
			self.current_ch.append(ref_ch) #this keeps a list of used reference chromosomes
			for i in range(len(self.fs)): #so bring all the other files up to the new reference
				f = self.fs[i]
				f_ch = self.lines[i].split(self.delimiter)[0]
				while f_ch != ref_ch:
					if f_ch not in self.current_ch:break #this means the files chromosome hasn't been discovered in the reference yet so stop advancing the file
					else:
						self.lines[i] = f.readline() #if f is still on the old chromosome keep going forward
						f_ch = self.lines[i].split(self.delimiter)[0] #this is the new chromosome for the next line in the while loop
		#ref_line_alleles = (self.ref_line.split(self.delimiter)[3], self.ref_line.split(self.delimiter)[4])
		if ref_index == self.previous_index:  #REFERENCE INDEL
			# if we advance each f_index - if it is equal to the ref-index than its\
			# there is information - take it and use it as usual
			# if it is not the same there is no information meaning its either missing \
			# or its moved to the next base because there is no indel in that individual or its a deletion
			# -either way insert the good_lines base because it will either match properly or be empty
			#INSERTION - grab the invariant individuals' IQ from previous position (ie its still in good_lines)
			#DELETION - grab all the IQs from good_lines
			for i in range(len(self.fs)):#the index i refers to the file
				f = self.fs[i]
				f_index = int(self.lines[i].split(self.delimiter)[self.index_position])
				f_ch = self.lines[i].split(self.delimiter)[0]
				if f_index == ref_index and f_ch == ref_ch: #there is information
					self.good_lines[i] = self.lines[i]
					self.lines[i]=self.fs[i].readline()
				elif f_index > ref_index and f_ch == ref_ch: #this means that there is no indel info in the vcf file f
					self.good_lines[i] = '.'
					self.lines[i]=self.fs[i].readline() 
					#there is no information
				else:
					self.good_lines[i] = '.'
			self.previous_index = ref_index
			return [self.ref_line] + self.good_lines
		#should this be an else statement?
		for i in range(len(self.fs)):
		#the index i refers to the file 
			#maybe I should check that the next line exists in the file - YES you moron you should
			f = self.fs[i]
			if self.lines[i]:
				# print (self.ref_line.strip())
				# print (self.lines[i])
				f_index = int(self.lines[i].split(self.delimiter)[self.index_position])
				f_ch = self.lines[i].split(self.delimiter)[0]
				if f_ch == ref_ch:
					#so what is happening is that once entered into this loop file f flips over to the next chromosome and then goes to the end of the vcf without ever being > ref_index
					while f_index < ref_index and f_ch == ref_ch: #and len(self.lines[i].split(self.delimiter)) > 2:
						self.lines[i] = f.readline()
						if self.lines[i]:#if you're at the end of a non-reference file this will not be true
							f_index = int(self.lines[i].split(self.delimiter)[self.index_position])
							f_ch = self.lines[i].split(self.delimiter)[0]
						else: #this means the file self.fs[i] is finished so just return an empty good line and move on
							self.good_lines[i] = '.'
							break
					# so this loop has ended either because (1) f_index >= ref_index OR (2) the f_index is less but the f_ch != ref_ch
					# in the case of (1): business as usual
					# in the case of (2): hold line as is by putting "." into good_lines and stop advancing line in the file f
					if ref_ch == f_ch:
						if f_index == ref_index:
							self.good_lines[i] = self.lines[i]
							self.lines[i]=f.readline()
						elif f_index < ref_index:
							self.good_lines[i] = '.'
							self.lines[i]=f.readline()
						elif f_index > ref_index:
							self.good_lines[i] = '.'
						else:	
							print(self.ref_line, ref_index, self.lines[i])
							sys.exit()
					else:
						self.good_lines[i] = '.'
				else:
					self.good_lines[i] = '.'
			else:
				self.good_lines[i] = '.'
		self.previous_index = ref_index
		return [self.ref_line] + self.good_lines
