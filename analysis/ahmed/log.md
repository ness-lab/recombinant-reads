
## 21/1/2020

today: configuring the shared repo! 

first - needed to set `git config core.sharedRepository' to world

then - change permissions of `analysis` (and `data`, though
that's not being committed) 

finally - head into `.git` and run `chmod a+rwX` to fix
permissions

## 4/3/2020

today - adapting Jimmy's `break_points.ipynb` notebook into a script to demo
argparse functionality

## 23/4/2020

sam -> sorted bam

```bash
samtools view -S -b recomb_mock.sam > recomb_mock.bam
samtools sort -T ../../../data/references/chlamy.5.3.w_organelles_mtMinus.fasta. -o recomb_mock.sorted.bam recomb_mock.bam
samtools index recomb_mock.sorted.bam # creates a .bam.bai file
```

## 28/4/2020

today: 

- test `parse_phase_changes.py` 
- append counters to it

here goes:

```bash
time python3.5 analysis/ahmed/phase_change_filter.py \
--filename analysis/ahmed/recomb_mock.sorted.bam \
--vcf data/liu-data/vcf/parental.vcf.gz \
--out analysis/ahmed/test_out.sam
```
