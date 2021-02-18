#!/usr/bin/env python3

import numpy.random as r
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# vcf2fasta command
# python2 vcf2fasta.py -v ../data/filtered_full.vcf.gz -r ../data/chlamy.5.3.w_organelles_mtMinus.fasta -i chromosome_1:1-1000000 > vcf2fasta.fasta

PHASE_CHANGE_PROBABILITY = 0.004
BASE_COUNT = 1000000
BASE_START = 100000
FASTA_FILEPATH = 'analysis/jimmy/power_analysis/vcf2fasta.fasta'
CHROMOSOME = 'chromosome_1'

events = r.poisson(BASE_COUNT * PHASE_CHANGE_PROBABILITY)
phase_changes = r.choice(BASE_COUNT, size=events)

# parse results from vcf2fasta
fasta_file = list(SeqIO.parse(FASTA_FILEPATH, 'fasta'))
haplotype0 = str(fasta_file[0].seq)
haplotype1 = str(fasta_file[1].seq)

current = 0
idx = 0
recomb_seq = []

for base in zip(haplotype0, haplotype1):

    # switch haplotype if phase_change    
    current = (current + 1) % 2 if idx in phase_changes else current
    recomb_seq.append(base[current])

    if idx >= BASE_COUNT:
        break
    else:
        idx += 1

ex_fasta = list(SeqIO.parse(FASTA_FILEPATH, 'fasta'))[0]

chrom = re.search(r'chromosome_[0-9]{1,2}', str(ex_fasta.id)).group()

recomb_record = SeqRecord(
    Seq(''.join(recomb_seq)),
    id=chrom,
    name=chrom,
    description=chrom,
    dbxrefs=ex_fasta.dbxrefs,
    features=ex_fasta.features,
    annotations=ex_fasta.annotations,
    letter_annotations=ex_fasta.letter_annotations
)

SeqIO.write(recomb_record, 'analysis/jimmy/power_analysis/recomb_fasta', 'fasta')

# art_illumina command
# art_illumina -f 400 -M -p -sam -l 250 -ss MSv1 -i recomb_fasta -m 600 -s 1 -o out -na

    

