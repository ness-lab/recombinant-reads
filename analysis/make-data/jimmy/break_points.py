#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import pysam

df_breakpoints = pd.read_csv('../../../data/references/liu-data-reference/all_break_points.txt', sep='\t')
df_lookup = pd.read_csv('../../../data/references/liu-data-reference/liu_sample_lookup.tsv', sep='\t')
bam_files_path = '../../../data/liu-data/bam/'
f = open('out.txt', 'w')

counter = 0
for index, row in df_breakpoints.iterrows():
    
    if counter > 0:
        break;
    counter += 1
    
    individual = row[2]
    chromosome = row[3]
    left_bound = int(row[4])
    right_bound = int(row[5])
    
    run = df_lookup[df_lookup['individual'] == individual].iloc[0]['Run'] + '.dup.bam'
    
    alignment_file = pysam.AlignmentFile(bam_files_path + run, 'rb')
    f.write(alignment_file.fetch(chromosome, left_bound, right_bound))
    


# In[ ]:




