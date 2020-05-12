
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
- Check if the two mate pairs combined create a phase change

Want to work on whether we can optimize these 3 checks

Because of new introduction of sequences with only 1 SNP, there are now sequences with deletions (cigar = 2), so need to work on that and checking if the SNPs line up or not.