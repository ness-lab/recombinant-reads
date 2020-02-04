
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