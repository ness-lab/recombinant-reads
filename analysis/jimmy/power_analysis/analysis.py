#!/usr/bin/env python3

import numpy.random as r
import numpy as np
import re
import os
import argparse
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def arg_parser():
    parser = argparse.ArgumentParser(
        description='generate sam reads for power analysis of readcomb',
        usage='python3 analysis.py [options]'
    )

    parser.add_argument('-f', '--fasta', required=True,
                        type=str, help='vcf2fasta output fasta file with 2 sequences each representing a haplotype')
    
    parser.add_argument('-s', '--start', required=False, default=1,
                        type=int, help='1-based index start of read generation (default is 1)')

    parser.add_argument('-e', '--end', required=False, default=0,
                        type=int, help='1-based index end of read generation (default is 0/all of vcf2fasta fasta output)')

    parser.add_argument('-p', '--percent', required=False, default=9.15 * (10 ** -8),
                        type=float, help='probability of phase change per base (default 9.15e-8)')

    parser.add_argument('-c', '--coverage', required=False, default=1,
                        type=int, help='read coverage of art_illumina (default 1)')

    parser.add_argument('-i', '--insert', required=False, default=100,
                        type=int, help='average insert size between paired reads (default 100)')

    parser.add_argument('-o', '--output', required=False, default='generated_reads',
                        type=str, help='prefix of art_illumina output (default generated_reads)')

    args = parser.parse_args()

    return {
        'fasta': args.fasta,
        'start': args.start,
        'end': args.end,
        'percent': args.percent,
        'coverage': args.coverage,
        'insert': args.insert,
        'output': args.output,
    }

def main():

    args = arg_parser()
    
    # parse results from vcf2fasta
    fasta_file = list(SeqIO.parse(args['fasta'], 'fasta'))

    # if end, change to all of chromosome
    if args['end'] == 0:
        args['end'] = len(fasta_file[0].seq)

    events = r.poisson((args['end'] - args['start']) * args['percent'])
    phase_changes = list(r.choice((args['end']- args['start']), size=events))

    if os.path.exists('output/phase_change_log.txt'):
        with open('output/phase_change_log.txt', 'a') as f:
            f.write('\n' + str(sorted(phase_changes)))

    # regex to identify correct fasta reference_name
    chrom = re.search(r'chromosome_[0-9]{1,2}', str(fasta_file[0].id)).group()

    haplotype0 = str(fasta_file[0].seq)
    haplotype1 = str(fasta_file[1].seq)

    current = 0
    idx = args['start']
    recomb_seq = []

    for base in zip(haplotype0, haplotype1):

        # switch haplotype if idx in list phase_changes    
        current = (current + 1) % 2 if idx in phase_changes else current
        recomb_seq.append(base[current])

        if idx >= args['end']:
            break
        else:
            idx += 1

    recomb_record = SeqRecord(
        Seq(''.join(recomb_seq)),
        id=chrom,
        name=chrom,
        description=chrom,
        dbxrefs=fasta_file[0].dbxrefs,
        features=fasta_file[0].features,
        annotations=fasta_file[0].annotations,
        letter_annotations=fasta_file[0].letter_annotations
    )

    SeqIO.write(recomb_record, '.recomb_fasta', 'fasta')

    # art_illumina command
    subprocess.run(['art_illumina',
                    '-f', str(args['coverage']),
                    '-l', '250',
                    '-ss', 'MSv1',
                    '-i', '.recomb_fasta',
                    '-m', str(500 + args['insert']),
                    '-s', '1',
                    '-o', args['output'],
                    '-M', '-p', '-sam', '-na', '-q'])

    # delete temp recomb_fasta file
    if os.path.exists('.recomb_fasta'):
        os.remove('.recomb_fasta')
    
if __name__ == '__main__':
    main()