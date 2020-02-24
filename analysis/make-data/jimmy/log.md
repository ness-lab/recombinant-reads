
## Jan 28th

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

## February 10th

[tmux cheat sheet](https://tmuxcheatsheet.com/)

[cyvcf2 documentation](https://brentp.github.io/cyvcf2/)

## February 18th

How to run jupyter notebook locally by tunneling into the hpcnode1 server:

Sign into the server and run jupyter notebook, make note of the port and token
```bash
ssh liujiyu@hpcnode1.utm.utoronto.ca
tmux attach
jupyter notebook
```

Then locally:
```bash
ssh -N -f -L localhost:[local port]:localhost:[port] liujiyu@hpcnode1.utm.utoronto.ca
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

## February 24th

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
