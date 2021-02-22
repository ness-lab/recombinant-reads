
# January 28th

Webpages
```
Samtools: http://biobits.org/samtools_primer.html
Pysam: https://pysam.readthedocs.io/en/latest/api.html
Pandas: https://pandas.pydata.org/pandas-docs/stable/getting_started/10min.html
```

Git configuration
```bash
git config --list
```

Location of breakpoint files
```
/recombinant-reads/data/references/liu-data-reference/all_break_points.txt
/recombinant-reads/data/references/liu-data-reference/liu_sample_lookup.tsv
```

Opening a tsv file using pandas
```python
df = pandas.open_csv('filepath', sep='\t')
```

Iterate through dataframe
```python
for index, row in df_breakpoints.iterrows():
    #do stuff
```

Select value by row
```python
df.iloc[row_index][column_name]
```

Open alignment file with pysam
```python
pysam.AlignmentFile(filepath, 'rb')
```

# February 10th

[tmux cheat sheet](https://tmuxcheatsheet.com/)

[cyvcf2 documentation](https://brentp.github.io/cyvcf2/)

# February 18th

How to run jupyter notebook locally by tunneling into the hpcnode1 server:

Sign into the server and run jupyter notebook, make note of the port and token
```bash
ssh liujiyu@hpcnode1.utm.utoronto.ca
tmux attach
jupyter notebook --no-browser
```

Then locally:
```bash
ssh -N -f -L localhost:[local_port]:localhost:[port] liujiyu@hpcnode1.utm.utoronto.ca
```

Note, if you are using a public wifi like UTM wifi, set your local port
to something uncommon like 23232 so that your port isn't also being used by someone else

Then go to the webpage:
```
localhost:8888/tree
```

Copy the token from the server to authenticate yourself

[vcf 4.2 documentation](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

[paper for example data we are using](https://www.nature.com/articles/s41559-017-0372-7)

For the example data, we only care about variants in cross H. The file is:
```python
Wvcf_files_path + 'H.vcf.gz'
```

# February 24th

When writing a bam file, if you want to transfer over the the header from another pysam alignment file, make sure to also include mode 'h'

```python
f = pysam.AlignmentFile('recomb_mock.bam', 'wh', template=recomb_file)
```

Keeps disconnecting because the school is not letting me use a port so I'm saving this chunk here just incase it hasn't saved on the server. Just give me a port bruh it's not a big deal you have like 40,000 of them

```python
'''
Given a bam file object, return the percentage of recombination in the bam files
'''

def recomb_percentage(bam_file_obj):
    
    recomb = 0
    total = 0
    
    for record in bam_file_obj:
        SNPs = check_SNPs(vcf_file_path + 'H.vcf.gz', record.reference_name, 
                          record.reference_alignment_start,
                          record.reference_alignment_start + record.query_alignment_length)
        print(SNPs)
```

[samtools view manual](https://www.htslib.org/doc/samtools-view.html)

[SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf)

[pysam documentation](https://pysam.readthedocs.io/en/latest/api.html)

Issues encountered:
- The VCF has two or more alternate alleles, not sure what that means
- For variants, just assuming that one parent is reference and other is alternate
- For records in bam file some don't have reference_end attribute, not sure what that means
- Not matching to either reference or alternate allele

# February 25th

REF is the reference string, not the actual lineages, use gt_bases to get the bases of the lineages.

For gaps and stuff look at cigar string so we can line up the query and the SNPs.

Write function to check how many SNPs each breakpoint has that's based on the mock breakpoint function

# March 5th

Did some fixing of the code and some SNP analysis. The results are in SNP_10k.txt which is where the SNP frequency per 10,000 bases was found; SNP_breakpoint.txt is the SNP frequency for the breakpoints that were provided. This file is alot shorter than all_breakpoints.txt because this only looked at SNPs on chromosome 1.

Analyzing the number of SNPs, it's quite low (< 2) for a lot of the breakpoints. Also for making the mock recombination, if the left bound and right bound are not in positions with lots of reads, the file may not be as long as we would like (200 per file, 2000 in total for 10 files). So I've pushed the bounds out of the breakpoints a bit which increases the reads but not by much.

The algorithm now uses gt_bases instead of REF and ALT which was wrong, but the recomb detection is far from perfect. All the recombinations should be at the start of the file but they are spread through the whole file and there are more mis-matches than recombinations. So the simplified algorithm still needs lots of work.

# March 10th

To find samples, run VCF(fname).samples

When the gt_bases is [T/T, G/G], it means that one allele is T and the other allele is T for the first sample and so on.
You can check uniqueness of the elements using a set.

record.query_alignment_qualities to get the array of qualities

record.cigartuples to get the array of tuples of alignment

Use SeqIO and print out 50 from reference and 50 from sample from a suspect sequence, and see if they are aligned

# March 17th

[UTM wiki](https://wiki.utm.utoronto.ca/display/WF/Wiki+FAQ)

# March 31st

Ask Ahmed where is the bam document with the reference sequence.

- Fasta file where that is, use SeqIO to parse the reference sequence

Working on cigar tuples

# April 14th

Created a diagnosis bam file with just sequences with no matches

There seems to be some problems with hard clipping and soft clipping that have to be fixed

VCF SNPs seems to be 0-based indexing while bam seems to be 1-based indexed, have to be careful when interacting the both of them

- This seems to be more complicated than I thought, something weird with the indexing will have to look into this later

record below refers to bam records in a bam file alignment object:

- Soft-clipping: record.query_sequence includes the part that soft clips with the sequence (the [(4, x), ...] section)

- Hard-clipping record.query_sequence does not include the part that hard clips with the sequence

# April 23

Doing some exploratory work on other breakpoints. There definitely aren't as many SNPs in the other breakpoints so it's a lot harder to find recombinations

Wrote code so that the diagnosis prints no matches in a bam format

Need to try and integrate diagnosis and recomb detection so that you can print all recombinations

Paralogous genes: VCF gt_bases filter out heterozygous calls because we're working with haploid organisms that shouldn't have heterozygous calls

Bug with SNPs: SNP index is the REF index but when you have an insertion the alignment of the ALT is then off because of that so the SNP checking needs to compensate for it

# May 5

VCFs are 1-indexed while bams are 0-indexed. I have decided to change the check_snp() function to add a 1 to the left-bound and right-bound input in order to try and keep the code consistent although extra has to be taken if working more directly with SNPs (using snp.start instead of snp.POS).

Bam reads that are mate pairs have the same read id that you can grab, most likely going to use a dictionary to store these values (cached reads are the most efficient)

# May 12th

Mate pairs are working for the human readable recombination and it has detected a lot more phase changes across the mate pairs. It also allows us to use all of our sequences that only have 1 SNP.

The current algorithm has to run 3 checks to write
- Check if 1st pair has a PC/NM if 2nd pair doesn't exist
    - Missing case here where 2nd pair has no snps but that is created afterwards
- Check if 2nd pair has a PC/NM
    - 2nd check can be removed? It can be removed if I always build snp_str
- Check if the two mate pairs combined create a phase change

Want to work on whether we can optimize these 3 checks

Because of new introduction of sequences with only 1 SNP, there are now sequences with deletions (cigar = 2), so need to work on that and checking if the SNPs line up or not.

Pairwise 2: Python-based aligner
from Bio import pairwise2

# May 13

Old human readable recombination detection code:

```python
def human_read_recomb(bam_file_obj, mode='no_match', output_filename='recomb_diagnosis'):
    '''
    Given a bam file object, return a human-readable file with aligned reference, quality string, 
        bam sequence, and SNPs    

    If mode='no_match', create a file with only no_match sequences.
    If mode='all', create a file with all sequences
    If mode='phase_change', create a file with only sequences with phase changes
    
    Makes assumptions on the following file locations:
        VCF: vcf_files_path + 'parental_filtered.vcf.gz'
        Reference: reference_files_path + 'chlamy.5.3.w_organelles_mtMinus.fasta'
    '''
    
    f_obj = open(output_filename + '.txt', 'w')

    f_obj.write('Key: \nSequence: Start - End \nReference sequence \nPhred Scale Quality \nQuery alignment sequence \n')
    f_obj.write('1: CC2935, 2: CC2936, N: Does not match SNP \n\n')
    
    # get reference segment
    seq_obj = SeqIO.parse(reference_files_path + 'chlamy.5.3.w_organelles_mtMinus.fasta', 'fasta')

    # grab chromosome 1
    for seq in seq_obj:
        chrom_1 = seq
        break    
        
    #counters
    no_match_counter = 0
    phase_change_counter = 0
    phase_change_no_match_counter = 0
    seq_with_snps_counter = 0
    all_seq_counter = 0
    
    for record in bam_file_obj:
        
        all_seq_counter += 1
                    
        snps = check_snps(vcf_files_path + 'parental_filtered.vcf.gz', record.reference_name, 
                          record.reference_start,
                          record.reference_start + record.query_alignment_length)
        
        if len(snps) > 1:
            
            seq_with_snps_counter += 1
            
            # tuple-checking
            cigar_tuples = record.cigartuples
            
            # initialize segment for building
            segment = ''
            
            # reference sequence
            ref = chrom_1[record.reference_start:record.reference_start + record.query_alignment_length] 
            
            # 0 is match and sequences are length of 150, so it's a full match
            if cigar_tuples == [(0, 150)]:
                segment = record.query_sequence
            else:
                # index to keep track of where we are in the query_segment
                query_segment = record.query_sequence
                index = 0
                
                for cigar_tuple in cigar_tuples:
                    # 4 = soft clipping, the record.query_sequence has the portion that is soft clipping
                    # so we need to skip it with index
                    if cigar_tuple[0] == 4:
                        index += cigar_tuple[1]
                        
                    # 5 = hard clipping, record.query_sequence does not have the portion that is
                    # hard clipping so we don't skip it and we don't add anything
                    elif cigar_tuple[0] == 5:
                        continue
                    
                    # if it is a match(0), then just add it onto the segment 
                    elif cigar_tuple[0] == 0:
                        segment += query_segment[index:index+cigar_tuple[1]]
                        index += cigar_tuple[1]
                    
                    # 1 is an insertion, we will add gaps to the reference
                    elif cigar_tuple[0] == 1:
                        segment += query_segment[index:index+cigar_tuple[1]]
                        index += cigar_tuple[1]
                        
                        ref = ref[:index] + '-' * cigar_tuple[1] + ref[index:]
                        
                    else:
                        print('oops forgot to consider this: ' + str(cigar_tuple))
                        print(cigar_tuples)
            
            snp_lst = [' '] * record.query_alignment_length
            
            for snp in snps:
                # Using SNP.start and record.reference_start since they are both 0 based
                # SNP.start grabs vcf positions in 0 index while vcfs are 1 indexed
                
                start = snp.start - record.reference_start
                
                # extra calculations to realign start if there is an insertion                    
                current_tuple = 0
                current_base = 0

                while current_base < start and current_tuple < len(cigar_tuples):
                    if cigar_tuples[current_tuple][0] == 1:
                        # shift the start over by the amount of insertion to compensate for it
                        start += cigar_tuples[current_tuple][1]

                    current_base += cigar_tuples[current_tuple][1]
                    current_tuple += 1
                            
                        
                
                # indexing for VCF seems to be a bit weird and will sometimes be -1
                if start < 0:
                    raise Exception('VCF indexing is off. Check SNP at {}'.format(snp))
                
                strand1 = snp.gt_bases[0][0]
                strand2 = snp.gt_bases[1][0]
                
                if start >= len(segment):
                    break
                
                if segment[start] == strand1:
                    snp_lst[start] = '1'
                    
                elif segment[start] == strand2:
                    snp_lst[start] = '2'
                    
                else:
                    snp_lst[start] = 'N'
                    
            snp_str = ''.join(snp_lst)
            
            # qualities string
            qualities_str = ''
            for quality in record.query_alignment_qualities:
                qualities_str += chr(quality + 33)
            
            #phase_change_counter update
            if '1' in snp_str and '2' in snp_str:
                phase_change_counter += 1
                
            #no_match_counter update
            if 'N' in snp_str:
                no_match_counter += 1
                
            #phase_change_no_match_counter update
            if '1' in snp_str and '2' in snp_str and 'N' in snp_str:
                phase_change_no_match_counter += 1            
            
            
            #there is not a no_match and the user only wants no_matches
            if 'N' not in snp_str and mode == 'no_match':
                continue
            
            #there is not a phase change and the user only wants phase changes
            elif not ('1' in snp_str and '2' in snp_str) and mode == 'phase_change':
                continue
                
            else:
                f_obj.write('Sequence: {start} - {end} \n'.format(start=record.reference_start, 
                                                               end=record.reference_start+record.query_alignment_length))

                f_obj.write(str(cigar_tuples) + '\n')

                f_obj.write(str(ref.seq) + '\n' + qualities_str + '\n' + segment + '\n' + snp_str + '\n \n')
    
    # print counters
    print('Sequences with phase changes: ' + str(phase_change_counter))
    print('Sequences with a no match: ' + str(no_match_counter))
    print('Sequences with a phase change and a no match: ' + str(phase_change_no_match_counter))
    print('Sequences with more than 1 SNP: ' + str(seq_with_snps_counter))
    print('No. of sequences: ' + str(all_seq_counter))
```

# May 17th

Ran into situation where record.cigartuples = None. The documentation says this is None when the alignment is not present. Made a check for it in both the phase_change function and the human_cigar function

# May 21st

The [MAF](http://genome.ucsc.edu/FAQ/FAQformat#format5) format is quite similar to the custom one I built for human readable, so it would be good to implement my human readable in that form

# May 25th

I've kind of come to realize that the cigar rebuilding is really useful for displaying a human readable recomb but it doesn't have to be as complicated for bams. Going to rewrite it so that we don't have to recompensate for insertions and deletions in phase_change

Bash command I've been using to test inside recomb folder

```bash
python3 phase_change_filter.py -f ../tests/mock_sams/mock_1.sam -v ../tests/parental_filtered.vcf.gz -o recomb_diagnosis
```

VCF Files
- AC = 2, AF = AC / AN, AN = 4
- For info, reference is 0 and alternate is 1
- GT is whether it's REF or ALT, there's 0/0 or 1/1 because we're pretending that the haploid specimen is a diploid
- Set GQ to 99, max quality
- AD the number of each, just set to 100, 0 or 50, 50 doesn't matter
- DP is the sum of the two numbers at AD
- PL phred scale likelihood - just don't do this
- Format defines the format of the info that comes after

# June 4th
- Test for soft clipping and hard clipping
- Insertions and deletions
- Test for 250 length reads
- Unaligned reads: about 1% just don't match and are encoded in the bams
     - It is in bit 4 that shows it is unaligned

# June 17th
- Fixed problem with deletions where it continued indexing while adding gaps
- Apparently for soft clipping the reference start is where the matching starts
- I have chosen to make mock_sequence 1-based sequencing because VCFs are 1-based indexed so it makes it easier to work with
- Still need to do all the argument stuff with human readable phase change
- Fixed bugs with recombinant reads

# June 18th
- Change function names for human readable version for organization
- Change sequences to mock_sequence functions
- Use setattr() for the mock_create

# June 23rd
- Change paired reads to read pairs to be more clear on what the variable means
- Change reference_files_path so it's not global 
- Make reference required

# July 2nd
- Add more verbosity to stuff that takes a long time like cachepairs
- It seems that there is name repeating, check if it happens over chromosomes
- Add argument so that chromosome is an arguement
- Add this check to cache pairs
```python
if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
    continue
```

# July 9th
- Maybe recreating the VCF everytime we run check_snps is slow? Find some way to pass around a generator
    - Watch out if reads might get eaten (they're not eaten on initial testing)
- Also create something to check if interval exists, don't know how right now

# July 16th
- Work on implementing pypy
- To install pip with pypy run these commands
```bash
wget https://bootstrap.pypa.io/get-pip.py
pypy3 get-pip.py
```
- Pypy does not work well with programs that use C, unable to pip instally either pysam or cyvcf2

# July 21st
- Working on multiprocessing vs multithreading, threading is good with programs with heavy I/O while processing is good with programs with heavy CPU load
- Was hoping to be able to use module `concurrent.futures` but pysam objects can't be pickled so we may have to run threads then write to separate files
- [SO page](https://stackoverflow.com/questions/18114285/what-are-the-differences-between-the-threading-and-multiprocessing-modules) about multi-threading vs multi-processing
- Did some performance analysis using cProfile and it seems that 85% of our performance load is cyvcf2 grabbing reads
- Used this [webpage](https://julien.danjou.info/guide-to-python-profiling-cprofile-concrete-case-carbonara/#:~:text=Profiling%20a%20Python%20program%20is,area%20might%20be%20worth%20optimizing.) to find out what we needed

# July 26th
- Using a mutable object and an index to pass through each thread in order to manipulate it without having to pass information back and forth: [SO page](https://stackoverflow.com/questions/6893968/how-to-get-the-return-value-from-a-thread-in-python)
    - Kind of confusing having a list of dictionaries of tuples of pysam records, gotta remember to write this in documentation
- Using itertools to split the dictionary of pairs in a kind of confusing but concise way: [SO page](https://stackoverflow.com/questions/12988351/split-a-dictionary-in-half)

# July 28th
- Consider what happens when `split_index` is 0 `if len(pairs) < arg["threads"]`
- Update documentation

# August 6th
- Fix bug on line 354, writes pair 0 twice instead of 0 and 1
- Remove sam extension at the end when outputting a file
- Add exception to cache reads if a read shows up 3 times

# August 11th
- Check for is_proper_pair when the sequence is_paired
- If there is a phase_change in one of the read but it's not a proper_pair should we include the one read
- Filter out snps where parent1 == parent2

# August 18th
- Insertions causing an offset by 1? Look into it to see what's actually happening (Ahmed is looking at this)
- Note on alignment: if there is an insertion and new mutation from no where, it is likely a misalignment issue
    - Not really applicable to the code, but worthy of note nonetheless
- Quality of reads filtering, shouldn't be in main script
- Argparse: action='store_true' or action='store_false' for whether or not we include improper pairs or not
- Update output so they no longer give percentages because they cause some annoying divide by zero errors

# August 20th
- Insertions weren't an issue, just some problems with the testing code by Ahmed
- Currently when you filter out non-proper reads, there are still some unpaired reads, gotta investigate this

# August 26th - Lab Meeting
- Filtering VCF pre-readcomb to catch weird SNPs
- It seems a lot of trouble comes from quality of alignment
- Hard filter on not looking at SNPs near indel?
- Unpaired reads could be signs of some bigger issue in the gene structure
- Create subprocess of vcftools that filters vcfs depending on the arguements
- scikit allele packaging data for compact storage of recombinantions

# August 27th
- Look into multiprocessing over multithreading
- Indel polymorphism detection in addition to SNPs

# September 2nd
- You might have one parent be AT and the other be ATAT and thus they might be both be correct
- For simple insertions and deletions in which either the REF or one of the ALT alleles would otherwise be null/empty, the REF and ALT Strings must include the base before the event (which must be reflected in the POS field)
    - I believe this is why we're getting errors with VCF indexing

# September 8th
- Weird anomaly where one of the bases is ./.
- sometimes snps are lining up with deletions
    - Nevermind that's just pre filtering

# September 10th
- Readable, add parent1 and parent2 sequences and align
- Insertions in readable use index and soft cliipping seems to be misalign it
    - Look at weird_read.sam

# September 14th
- Looked into weird_read.sam, found that we weren't compensating for soft clipping correctly in phase_change_readable
    - We need a separate ref_id to keep track of where we are on the reference because soft cliped parts isn't included when we retrieve reference by `chrom_1[record.reference_start:record.reference_start + record.query_alignment_length]` where before I assumed that it was
    - Also, soft clipped and hard clipped parts aren't considered when we're doing extra calculations for recompensating for insertions during phase change detection
- Surprisingly, adding a ref_id makes the cigar parsing sequence more intuitive because it more closely matches how BAM files consider insertions and deletions through whether each sequence is iterated through shown on the sam specification table
- Still got a no match in this sequence inside of an insertion so gotta look into that
- Not sure about the parent1 parent2 thing, it can either be based off of the reference or the segment I don't know which one will be more useful

# September 15th
- Do we need to readjust for deletions after we adjusted the cigar algorithm?
- Rebuilt cigar_human so that it is closer to sam specifications
- Do this after phase detection

# September 23rd
- Some weird stuff happening with parent 1 and parent 2 alignment, gotta look into it'

# October 7th
- Looking into multiprocessing, could use some shared memory but the buffer system is kind of weird and Python recommends not to share memory space
- We could create classes for bam objects in order to pass them around
- There could be a very clean implementation if we don't consider pairs
    - Gotta have Ahmed look into if pairs are actually giving us good results

# October 8th
- Bam files, organize them in a way where reads are right after each other
    - Samtools sort by name?

# November 3rd
- Logging support
- Documentation
- VCF suppression
- Smarter iteration counter

# November 5th
- Improve phase change algorithm
- tqdm bars for each process

# November 24th
- SNP/phase change density
    - Usable variants on x, phase change reads on y
- Recombination classification

# November 26th
- Make histogram of number of variants in each bam sequence
    - Full bam / filtered reads comparison
- Downstream script
    - Chlamydomonas median NCO-GC length: 73, average: 40-50, use 100 bp where the distribution tapers off
    - Create read objects that proces on initiation and can call to give information

# December 2nd
- Histogram of number of variants in each bam PAIR

# December 8th
- Change name check in scheduler of filter.py to use query_name instead of reference_name; reference_name refers to the chromosome
- Change order of parents in phase_change so that they differ when the length of indels of parent 2 is longer
- Change the variable snp to variant in phase_change to be more accurate in terminology

# December 9th
- [Regex tester](https://regex101.com/)

# December 10th
- Change strand to haplpotype/parent
- Fix issue with indexing pairs
- Extend ranges to the end of the reads
- Use length of list instead of regex for CO/NCO identification

# December 15th
- Masking - Ahmed is asking Rob
- 2121/1212 - complex classification
- You can't have crossovers next to each other - crossover interference
- SDSA - gene conversion
- [Super intense phase change link](https://www.nature.com/articles/nrm2008)

# December 18th
- ambiguous/unambiguous crossovers
- different masking sizes - parameter
- clear crossovers, gene conversions, ambiguous, complex
- actual clear crossover: happens 70 to 100 bp away from ends

# December 22nd
- Instead of counting, figure out the base number that the start and end are at
- midpoint method for self.simplify

# January 4th
- Update usage of filter.py, currently still says phase_change_Filter.py
    - Probably gonna be obsolete after we make it into bash program
- Do stuff noted in December 22nd
- Look into setuptools to see if that's what we want to use for the command line program

# January 11th
- [setuptools](https://packaging.python.org/tutorials/packaging-projects/)

# January 18th
- [setuptools entry points](https://setuptools.readthedocs.io/en/latest/userguide/entry_point.html?highlight=scripts)
- [setuptools packages](https://setuptools.readthedocs.io/en/latest/userguide/package_discovery.html)
- [setuptools config](https://setuptools.readthedocs.io/en/latest/userguide/declarative_config.html)
- [old setup script](https://docs.python.org/3/distutils/setupscript.html)
- [distutil examples](https://docs.python.org/3/distutils/examples.html)
- Readcomb pacakge is working
    - readcomb-filter --bam [bam] --vcf [vcf]
    - readcomb-vcfprep --bam [bam] --vcf [vcf]
    - in python: `import readcomb.classification` as rc or `from readcomb import classification`

# January 21st
- Create the build releases: `python3 setup.py sdist bdist_wheel`
- Upload to pypi: `twine upload --repository-url https://upload.pypi.org/legacy/ dist/*`
- Updating to pypi, create the new build then `twine upload dist/*`
- Look at changing description-content-type to markdown - done

# January 24th
- Ask Ahmed about bamprep, is it done? Also does it need argument for samtools binary?
- Currently pairs in classification cannot be passed in processes, maybe write pickle and unpickle functions?

# January 25th
- Change call to classify and then classify to call
- Package and unpackage script for classification
- Look for hardcoded stuff
- There are sequences with no cigartuples, just skip those - already had this done nice
cd ..
# January 28th
- try/except block in classification for from readcomb.filter, have warning and then import from local directory
- power analysis: 1 - false negatives

# February 1st
- [Comparison of read simulation]https://www.nature.com/articles/nrg.2016.57)
- [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)

# February 4th
- ART needs g++ and gsl `sudo apt-get install libgsl-dev`

# February 8th
- Random chance: numpy.random.choice [link](https://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice)
- Dunno about to_dict, could just iterate through the one at a time
- ART Illumina ask about:
    - fold of read coverage -f --fcov
    - insertion/deletion rate
    - cigarM: seems to be the stuff in sam cigar numbers 7 and 8
    - quality profile builder - seems useful
- apt-get install art-nextgen-simulation-tools
- raw fastq come from sequencing company, then bam/sams are sequences aligned to the reference sequence
- use default miseq 250bp profile and use sequencing profile later
- stick with cigarM
- watch out for matepair vs paired, matepair is something where the space between two pairs that's >5000 bp and not what we're using
- len is the length of the sequence, mflen is the mean size of both sequences in the pair plus the average space in the middle
- f/fcov is just coverage, start with 20, final should be 400
- stick to MSv1 for builtin reads

# February 18th
- Include art_illumina command into analysis.py and include argparse to take arguements
- Draft out a bash script

# February 21st
- vcf2fasta command `python2 vcf2fasta.py -v ../data/filtered_full.vcf.gz -r ../data/chlamy.5.3.w_organelles_mtMinus.fasta -i chromosome_1:1-1000000 > vcf2fasta.fasta`
- art_illumina command `art_illumina -f 400 -M -p -sam -l 250 -ss MSv1 -i recomb_fasta -m 600 -s 1 -o generated_reads -na`