
## 21/1/2020

today: configuring the shared repo! 

first - needed to set `git config core.sharedRepository` to world

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

## 12/5/2020

today: 

- SNP density tests for reads in `recomb_mock.bam` and maybe others as well

could try with mate pairs and without

brute force, without mate pairs:

```python
import pysam
from cyvcf2 import VCF
from tqdm import tqdm

fname = 'recomb_mock.sorted.bam'
vcf = '../../data/liu-data/vcf/parental_filtered.vcf.gz'

bam_reader = pysam.AlignmentFile(fname)

def check_snps(f_name, chromosome, left_bound, right_bound):
    vcf_in = VCF(f_name)
    region = '{c}:{l}-{r}'.format(c=chromosome, l=left_bound+1, r=right_bound+1)
    records = [rec for rec in vcf_in(region)]
    return len(records)

counts = {}
for r in tqdm(bam_reader):
    count = check_snps(vcf, 'chromosome_1', r.reference_start+1, r.reference_end+1)
    if count in counts:
        counts[count] += 1
    else:
        counts[count] = 1
# {0: 182, 1: 56, 2: 9, 3: 85, 4: 62, 5: 85, 6: 53, 7: 5, 8: 70, 9: 86, 10: 17}

# trying this with mate pairs:
def cache_pairs(bam_file_obj):
    
    cache = {}
    
    paired = 0
    total = 0
    
    for record in tqdm(bam_file_obj):
        name = record.query_name
        total += 1
        
        if name not in cache:
            cache[name] = [record,None]
        else:
            cache[name][1] = record
            paired += 1
            
    print('Paired: {paired}, Unpaired: {unpaired}'.format(paired=paired, unpaired=total-paired))    
    return cache

bam_reader = pysam.AlignmentFile(fname)
cache = cache_pairs(bam_reader)
# Paired: 241, Unpaired: 469

counts_paired = {}
for key, val in cache.items():
    read1, read2 = val
    if not read2: # unpaired
        continue
    combined_count = sum([check_snps(vcf, 'chromosome_1', r.reference_start+1, r.reference_end+1) for r in val])
    if combined_count in counts_paired:
        counts_paired[combined_count] += 1
    else:
        counts_paired[combined_count] = 1

# {0: 43, 1: 12, 2: 1, 3: 8, 4: 8, 5: 14, 6: 20, 7: 16, 8: 12, 9: 7, 10: 14, 11: 16, 12: 15, 13: 15, 14: 28, 15: 12}

# how many reads/read pairs have 3 or more SNPs? 
sum([v for k, v in counts.items() if k > 2])
# 463
sum([v for k, v in counts_paired.items() if k > 2])
# 185

# trying this with one of the original bams
# SRR5243250.dup.bam has 25433936 (!) reads -

full_bam = '../../data/liu-data/bam/SRR5243250.dup.bam'
full_bam_reader = pysam.AlignmentFile(full_bam)

full_counts = {}
full_counts_paired = {}

cache = cache_pairs(full_bam)
# 25433936it [02:08, 198357.15it/s]
# Paired: 13406587, Unpaired: 12027349
```

shoot - the parental vcfs only contain chromosome 1 - need to create full VCFs

this was a holdover from when I made this VCF in `rcmb-seq-schemas`, where I just
made chromosome 1 to save time (and yet it still took 4 hours, lol)

in the meantime, let's try the above again but for just chromosome 1

```
from tqdm import tqdm
import pysam
from cyvcf2 import VCF

vcf_fname = 'data/liu-data/vcf/parental_filtered.vcf.gz'
bam_fname = 'data/liu-data/bam/SRR5243250.dup.bam'

def count_snps(f_name, chromosome, left_bound, right_bound):
    vcf_in = VCF(f_name)
    region = '{c}:{l}-{r}'.format(c=chromosome, l=left_bound+1, r=right_bound+1)
    records = [rec for rec in vcf_in(region)]
    return len(records)


def cache_pairs_chr1(bam_file_obj):
    cache = {}
    paired = 0
    total = 0
    bonked_reads = 0
    
    for record in tqdm(bam_file_obj):
        if record.reference_id == -1:
            bonked_reads += 1
            continue
        if record.reference_name != 'chromosome_1':
            continue
        name = record.query_name
        total += 1
        
        if name not in cache:
            cache[name] = [record,None]
        else:
            cache[name][1] = record
            paired += 1
                
    print('Paired: {paired}, Unpaired: {unpaired}'.format(paired=paired, unpaired=total-paired))
    print('Bonked reads: {}'.format(bonked_reads))
    return cache

bam_file = pysam.AlignmentFile(bam_fname)
cache = cache_pairs_chr1(bam_file)
# 25433936it [01:11, 356735.88it/s]
# Paired: 676920, Unpaired: 840929
# Bonked reads: 512350 -> had reference_id = -1 -> what does this mean? 

counts_paired = {}

for key, val in tqdm(cache.items()):
    read1, read2 = val
    if not read2: # unpaired
        continue
    combined_count = sum([count_snps(vcf, 'chromosome_1', r.reference_start+1, r.reference_end+1) for r in val])
    if combined_count in counts_paired:
        counts_paired[combined_count] += 1
    else:
        counts_paired[combined_count] = 1

```

also - I need to make a full, genomic VCF between just 2935 and 2936. this is probably
best done a few chromosomes at a time

I've already got the requisite files ready in the `rcmb-seq-schemas` folder, so I'm just going
to run the commands from 24/07/2019 in that project's `tetrad-draws` log in there, combine the files,
and then symlink the final file here

```bash
# THESE PATHS ARE IN genomewide_recombination/rcmb-seq-schemas SO RUN THIS THERE

time for i in {2..9}; do
    time java -jar bin/GenomeAnalysisTK.jar \
    -R data/references/chlamy_reference.fasta \
    -T HaplotypeCaller \
    -I data/tetrad-draws/bams/CC2935.RG.bam \
    -I data/tetrad-draws/bams/CC2936.RG.bam \
    -L chromosome_${i} \
    -o data/tetrad-draws/vcfs/chromosome_${i}.vcf;
done # took 1.25 days...

time for i in {10..17}; do
    time java -jar bin/GenomeAnalysisTK.jar \
    -R data/references/chlamy_reference.fasta \
    -T HaplotypeCaller \
    -I data/tetrad-draws/bams/CC2935.RG.bam \
    -I data/tetrad-draws/bams/CC2936.RG.bam \
    -L chromosome_${i} \
    -o data/tetrad-draws/vcfs/chromosome_${i}.vcf;
done
```

## 2/7/2020

shoot - need to add mtDNA and cpDNA

```bash
echo "mtDNA" > organelles.intervals
echo "cpDNA" >> organelles.intervals
time java -jar bin/GenomeAnalysisTK.jar \
-R data/references/chlamy_reference.fasta \
-T HaplotypeCaller \
-I data/tetrad-draws/bams/CC2935.RG.bam \
-I data/tetrad-draws/bams/CC2936.RG.bam \
-L organelles.intervals \
-o data/tetrad-draws/vcfs/organelles.vcf
```

## 9/7/2020

add scaffolds and organelles to `parental_full`

```bash
touch scaffolds.intervals
for i in {18..54}; do
    echo "scaffold_${i}" >> scaffolds.intervals
done
echo "mtMinus" >> scaffolds.intervals
    
time java -jar bin/GenomeAnalysisTK.jar \
-R data/references/chlamy_reference.fasta \
-T HaplotypeCaller \
-I data/tetrad-draws/bams/CC2935.RG.bam \
-I data/tetrad-draws/bams/CC2936.RG.bam \
-L scaffolds.intervals \
-o data/tetrad-draws/vcfs/scaffolds.vcf
```

## 18/7/2020

install pypy on server:

```bash
cd ~/apps
wget https://bitbucket.org/pypy/pypy/downloads/pypy3.6-v7.3.1-linux64.tar.bz2
tar xjf pypy3.6-v7.3.1-linux64.tar.bz2
mv -v pypy3.5-v7.3.1-linux64 pypy3.6

# to project folder
mkdir -p bin/
cd bin/
ln -sv ~/apps/pypy3.6/bin/pypy3 .
```

## 21/7/2020

today: test whether a bam can be written out of order
and then put 'back in order' via `samtools sort`

```python
>>> import pysam
>>> from tqdm import tqdm
>>> fname = 'data/liu-data/bam/CC2935.RG.bam'
>>> chrs = ['chromosome_3', 'chromosome_2', 'chromosome_1']
>>> bam_in = pysam.AlignmentFile(fname, 'r')
Warning: The index file is older than the data file: data/liu-data/bam/CC2935.RG.bai
>>> bam_out = pysam.AlignmentFile('test.sam', 'wh', template=bam_in)
>>> for chrom in chrs:
  2     reader = bam_in.fetch(chrom)
  3     for record in tqdm(reader):
  4         bam_out.write(record)
10096504it [00:29, 341206.27it/s]
10015743it [00:34, 287836.71it/s]
9491147it [00:33, 286559.67it/s]
>>> bam_out_test = pysam.AlignmentFile('test_presorted.sam', 'wh', template=bam_in)
>>> for chrom in chrs[::-1]: # in correct order
  2     reader = bam_in.fetch(chrom)
  3     for record in tqdm(reader):
  4         bam_out_test.write(record)
9491147it [00:27, 349873.16it/s]
10015743it [00:29, 333985.53it/s]
10096504it [00:29, 337003.14it/s]
```

diffing the entire files is a little impractical, so let's compare
the first 100 lines:

```bash
diff <(head -n 1000 test.sam) <(head -n 1000 test_presorted.sam) | wc -l
# 1844
```

sorting and then trying again:

```bash
samtools sort -O sam -o test_sorted.sam test.sam # took 3m44s for 3 chromosomes

diff <(head -n 1000 test_sorted.sam) <(head -n 1000 test_presorted.sam) | wc -l
# 4

diff <(head -n 1000 test_sorted.sam) <(head -n 1000 test_presorted.sam)
# 1c1
# < @HD   VN:1.4  SO:coordinate
# ---
# > @HD   VN:1.4  GO:none SO:coordinate
```

looks good! 

now for the string lookup test. here's a standard run:

```bash
# use python to create a chr1 only version of SRR5243250.dup.bam first

time python3.5 readcomb/phase_change_filter.py \
--bam chr1.sam \
--vcf tests/data/parental_full.vcf.gz \
--out test.out
```

didn't let it get to completion, but it was projected to take 20 min
at ~650 it/s

now to refactor with the string lookup

moment of truth:

```bash
time python3.5 readcomb/phase_change_filter.py \
--bam chr1.sam \
--vcf tests/data/parental_full.vcf.gz \
--out test.out
```

## 6/8/2020

(it's 1 am, lol)

quick test with log on a test bam:

```bash
time python3.5 readcomb/phase_change_filter.py \
--bam chr1.sam \
--vcf tests/data/parental_full.vcf.gz \
--threads 4 \
--chrom chromosome_1 \
--mode phase_change \
--log test.log \
--out test.out.sam
```

## 10/8/2020

need to remake these files since there was a bug that was writing the first read in a
pair twice

with and without threading:

```bash
time python3.5 readcomb/phase_change_filter.py \
--bam chr1.sam \
--vcf tests/data/parental_full.vcf.gz \
--threads 4 \
--chrom chromosome_1 \
--mode phase_change \
--log test.log \
--out test.out.threaded.sam # done in 7 min

time samtools sort -@4 -O sam -o test.out.threaded.sorted.sam test.out.threaded.sam

time python3.5 readcomb/phase_change_filter.py \
--bam chr1.sam \
--vcf tests/data/parental_full.vcf.gz \
--chrom chromosome_1 \
--mode phase_change \
--log test.log \
--out test.out.nothread.sam

time samtools sort -@4 -O sam -o test.out.nothread.sorted.sam test.out.nothread.sam

```

## 11/8/2020

so we're letting lots of crappy SNPs through! updating the filters in `check_snps`:

```python
time python3.5 readcomb/phase_change_filter.py \
--bam chr1.bam --vcf tests/data/parental_full.vcf.gz \
--chrom chromosome_1 --mode phase_change --threads 4 \
--log test.log --out test.out.snps.sam
```

once more after fixing another issue with SNPs where parent1 == parent2:

```python
time python3.5 readcomb/phase_change_filter.py \
--bam chr1.bam --vcf tests/data/parental_full.vcf.gz \
--chrom chromosome_1 --mode phase_change --threads 4 \
--log test.log --out test.out.diffs.sam
```

so it looks like insertions/deletions, whether in the parents or
the recombinant, are making it hard to clearly see phase changes - need
to dig into this more

also, it seems that the bams aren't lining up properly with the depth values
I see in IGV here - might need to download the full bams instead of just working
with this filtered chromosome 1 bam? and also check how much concordance we expect
in these values to begin with

## 17/8/2020

alright - many hours spent staring at IGV later, found a few bugs and a few fixes, detailed
in notion rn

re depth - these are expected to differ, since variant callers will filter out reads for whatever
reason and thus have lower depth values than bam coverage

now that the bug re proper pairs has been fixed:

```bash
time python3.5 readcomb/phase_change_filter.py \
--bam chr1.bam --vcf tests/data/parental_full.vcf.gz \
--chrom chromosome_1 --mode phase_change --threads 4 \
--log test.log --out test.out.pairs.sam
```

maybe need to regen the VCF with more specific parameters? (eg heterozygosity)

just making chr1, back in `rcmb-seq-schemas`:

```bash
time java -jar bin/GenomeAnalysisTK.jar \
-R data/references/chlamy_reference.fasta \
-T HaplotypeCaller \
-I data/tetrad-draws/bams/CC2935.RG.bam \
-I data/tetrad-draws/bams/CC2936.RG.bam \
-L chromosome_1 \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/tetrad-draws/vcfs/chromosome_1_redone.vcf;
```

seems this yields more variants (221k vs 212k) than the original, although
that misaligned area near 103k is not treated any differently 

still, should remake `parental_full` in light of this, and maybe delete some old files

## 21/8/2020

for now - symlinking the new chr1 vcf and creating a new output bam from it

we've also fixed some bugs in the phase change script that should lead to different output

```bash
time python3.5 readcomb/phase_change_filter.py \
--bam chr1.bam --vcf tests/data/chromosome_1_redone.vcf.gz \
--chrom chromosome_1 --mode phase_change --threads 4 \
--log test.log --out test.out.redone.sam
```

## 24/8/2020

need to get stats for lab meeting presentation after pulling latest changes

```bash
time python3.5 readcomb/phase_change_filter.py \
--bam chr1.bam --vcf tests/data/chromosome_1_redone.vcf.gz \
--chrom chromosome_1 --mode phase_change --threads 4 \
--log test.new.log --out test.out.proper.sam
```

and on a full bam:

```bash
time python3.5 readcomb/phase_change_filter.py \
--bam ../recombinant-reads/data/liu-data/bam/SRR5243250.dup.bam \
--vcf tests/data/parental_full.vcf.gz \
--mode phase_change --threads 3 \
--log test.new.log --out test.out.full.sam
``` 

took 2:37:07 with 3 threads - not bad

```
Done.
0 phase changes reads from 0 total unpaired
106236 phase changes reads across mate pairs from 19944774 total read pairs
166808 reads had no-match variants.
10971995 reads did not have enough SNPs (> 0) to call
time taken: 2:37:07
```

using this fxn:

```python
def check_rec(read_list):
    print(read_list[0].reference_start, read_list[0].reference_end, read_list[1].reference_start, read_list[1].reference_end)
    print(check_snps(vcf_in, 'chromosome_1', read_list[0].reference_start, read_list[0].reference_end))
    print(phase_detection(check_snps(vcf_in, 'chromosome_1', read_list[0].reference_start, 
            read_list[0].reference_end), cigar(read_list[0]), read_list[0]))
    print(check_snps(vcf_in, 'chromosome_1', read_list[1].reference_start, read_list[1].reference_end))
    print(phase_detection(check_snps(vcf_in, 'chromosome_1', read_list[1].reference_start, 
            read_list[1].reference_end), cigar(read_list[1]), read_list[1]))
```

## 28/8/2020

trying this with a prefiltered VCF while removing filters from the script itself:

```bash
./bin/bcftools filter -i 'MIN(GQ)>29 && TYPE="snp"' \
tests/data/chromosome_1_redone.vcf.gz > \
tests/data/chromosome_1_redone_filtered.vcf

bgzip tests/data/chromosome_1_redone_filtered.vcf
tabix -p vcf tests/data/chromosome_1_redone_filtered.vcf.gz
```

wait - can't filter out sites where the genotypes differ between parents
using just bcftools

might have to change readcomb logic to account for this, but in the interim
will just use `cyvcf2` to create the new VCF instead 

```python
from cyvcf2 import VCF
from cyvcf2 import Writer
from tqdm import tqdm
vcf_in = VCF('tests/data/chromosome_1_redone.vcf.gz')
vcf_out = Writer('tests/data/chromosome_1_filtered.vcf', vcf_in)
vcf_out.write_header()
for rec in tqdm(vcf_in):
    if rec.is_snp and len(rec.ALT) > 0 and rec.num_het == 0:
        if all(rec.gt_quals >= 30) and rec.gt_bases[0] != rec.gt_bases[1]:
            vcf_out.write_record(rec)
```

and then after bgzipping and tabixing, let's give this a go:

```bash
time python3.5 readcomb/phase_change_filter.py \
--bam chr1.bam --vcf tests/data/chromosome_1_filtered.vcf.gz \
--threads 3 --chrom chromosome_1 --mode phase_change \
--log test.updated.log --out test.out.vcf.sam
```

so that was WAY faster - making a VCF preprocessing script - here's a test run:

```bash
time python3.5 readcomb/vcfprep.py \
--vcf tests/data/chromosome_1_redone.vcf.gz \
--snps_only --no_hets --min_GQ 30 \
--out test.vcf
```

looks good! 

## 2/9/2020

testing out indel only functionality:

```bash
time python3.5 readcomb/vcfprep.py \
--vcf tests/data/chromosome_1_redone.vcf.gz \
--indels_only --no_hets --min_GQ 30 \
--out indel.test.vcf
```

## 10/9/2020

trying out jimmy's indel script, first with readable:

```bash
time python3.5 readcomb/phase_change_readable.py --bam chr1.bam \
--vcf tests/data/parental_full.vcf.gz \
--chrom chromosome_1 \
--reference ../recombinant-reads/data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--out readable.out.txt
```

## 1/11/2020

pushing jimmy's code to the server and giving it a test run

before I can do that, have to make an updated pre-filtered VCF that
contains the entire genome and not just chromosome 1:


```python
from cyvcf2 import VCF
from cyvcf2 import Writer
from tqdm import tqdm
vcf_in = VCF('tests/data/parental_full.vcf.gz')
vcf_out = Writer('tests/data/parental_full_filtered.vcf', vcf_in)
vcf_out.write_header()
for rec in tqdm(vcf_in):
    if rec.is_snp and len(rec.ALT) > 0 and rec.num_het == 0:
        if all(rec.gt_quals >= 30) and rec.gt_bases[0] != rec.gt_bases[1]:
            vcf_out.write_record(rec)
```

back to bash:

```bash
# bgzip and tabix
bgzip tests/data/parental_full_filtered.vcf
tabix -p vcf tests/data/parental_full_filtered.vcf.gz

# naive run without setting processes
# needed to upgrade pysam for this
time python3.5 readcomb/phase_change_process.py \
--bam ../recombinant-reads/data/liu-data/bam/test_sort.sam \
--vcf tests/data/parental_full_filtered.vcf.gz \
--out test_out/processing_test.sam
```

this seems to have trouble dealing with some sequences
that don't have SNPs - raising warnings like this:

```
Process 0 at 1000 iterations
Process 0 at 2000 iterations
no intervals found for b'tests/data/parental_full_filtered.vcf.gz' at scaffold_42:3012-3162
no intervals found for b'tests/data/parental_full_filtered.vcf.gz' at scaffold_42:3207-3357
no intervals found for b'tests/data/parental_full_filtered.vcf.gz' at scaffold_27:43694-43844
no intervals found for b'tests/data/parental_full_filtered.vcf.gz' at scaffold_27:43884-44034
```

seems to have taken four hours:

```bash
21206
1501845 phase change reads pairs from total 9972387 read pairs
535855 reads had no-match variants
5315828 reads did not have enough SNPs (> 0) to call
time taken: 4:00:02

real    240m2.459s
user    257m34.914s
sys     20m36.379s
```

trying this with 4 processes instead:

```bash
time python3.5 readcomb/phase_change_process.py \
--bam ../recombinant-reads/data/liu-data/bam/test_sort.sam \
--vcf tests/data/parental_full_filtered.vcf.gz \
--processes 4 \
--out test_out/processing_test_4.sam
```

seems time taken is reduced to 1/4 the original:

```
1501845 phase change reads pairs from total 9972387 read pairs
535855 reads had no-match variants
5315828 reads did not have enough SNPs (> 0) to call
time taken: 1:02:35

real    62m35.874s
user    265m5.745s
sys     21m30.657s
```

## 03/11/2020

16 processes:

```bash
time python3.5 readcomb/phase_change_process.py \
--bam ../recombinant-reads/data/liu-data/bam/test_sort.sam \
--vcf tests/data/parental_full_filtered.vcf.gz \
--processes 16 \
--out test_out/processing_test_16.sam
```

holy shit this was fast! 

```
1501845 phase change reads pairs from total 9972387 read pairs
535855 reads had no-match variants
5315828 reads did not have enough SNPs (> 0) to call
time taken: 0:16:01

real    16m1.603s
user    263m33.411s
sys     12m23.357s
```

need to check whether all these output files are identical, but so
far this is looking great

## 4/11/2020

trying out this hidden prints solution on stack overflow
to suppress cyvcf2 output:

```python
# https://stackoverflow.com/questions/8391411/how-to-block-calls-to-print
import os
import sys

class Silent:
    def __enter__(self):
        self._original_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stderr.close()
        sys.stderr = self._original_stderr

# cyvcf2 code -
# https://github.com/brentp/cyvcf2/blob/master/cyvcf2/cyvcf2.pyx
```

the key is to suppress stderr and not stdout, based on
the actual cyvcf2 code (which calls on `sys.stderr.write`)

## 5/11/2020

running with 32 processes

```bash
time python3.5 readcomb/phase_change_process.py \
--bam ../recombinant-reads/data/liu-data/bam/test_sort.sam \
--vcf tests/data/parental_full_filtered.vcf.gz --processes 32 \
--out test_out/processing_test_32.sam
```

this seems to have been slower somehow:

```
Creating processes
Processes created
25433936it [21:09, 20041.22it/s]
Waiting for processes to complete
1501845 phase change reads pairs from total 9972387 read pairs
535855 reads had no-match variants
5315828 reads did not have enough SNPs (> 0) to call
time taken: 0:21:09

real    21m10.665s
user    324m4.088s
sys     34m18.964s
```

## 8/11/2020

checking whether at least this output file can
be sorted correctly:


```bash
time samtools sort -n test_out/processing_test_32.sam > 
test_out/processing_test_32.sorted.sam
```

still getting the same error that seems to be common
to the earlier files as well (didn't take notes during
the meeting the other day but this *was* happening then) -
it seems the last file isn't written in full:

```
[W::sam_read1] Parse error at line 3003750
samtools sort: truncated file. Aborting

real    0m29.317s
user    0m23.690s
sys     0m2.185s
```

this could be an issue with the process ending
before the last file is written - added a `time.sleep`
bandaid the other day but seems that didn't make a difference

going to try to debug this by running on a smaller input
file (`test_sort_20k.sam`). let's try a naive run without
changing anything else:

```bash
time python3.5 readcomb/phase_change_process.py \
--bam ../recombinant-reads/data/liu-data/bam/test_sort_20k.sam \
--vcf tests/data/parental_full_filtered.vcf.gz \
--processes 16 \
--out test_out/test_short_16.sam
```

that was done in 2 seconds - nice! 

sorting:

```bash
time samtools sort -n test_out/test_short_16.sam > \
test_out/test_short_16.sorted.sam
```

same issue:

```
[E::sam_parse1] SEQ and QUAL are of different length
[W::sam_read1] Parse error at line 2369
samtools sort: truncated file. Aborting

real    0m0.010s
user    0m0.007s
sys     0m0.002s
```

increasing sleep time to 5 seconds and rerunning:

```bash
time python3.5 readcomb/phase_change_process.py \
--bam ../recombinant-reads/data/liu-data/bam/test_sort_20k.sam \
--vcf tests/data/parental_full_filtered.vcf.gz \
--processes 16 \
--out test_out/test_short_16.sam
```

same error... ten seconds?? 

no dice still - and it also doesn't appear from the eye
test that increasing the time allows for 'more' of the record to be
written

wait - does this still happen if we just use one process?
removing all calls to `time.sleep` for this:

```bash
time python3.5 readcomb/phase_change_process.py \
--bam ../recombinant-reads/data/liu-data/bam/test_sort_20k.sam \
--vcf tests/data/parental_full_filtered.vcf.gz \
--processes 1 --out test_out/test_short_1.sam
```

WELL THEN

```
[E::sam_parse1] SEQ and QUAL are of different length
[W::sam_read1] Parse error at line 2369
samtools sort: truncated file. Aborting
```

so this happens even with one process... back to the drawing board

checking if `SilentVCF` slows the script down, it seems that
using one process on the shorter file completes in 13.526 seconds
with SilentVCF enabled and 12.704 seconds without - which is
pretty tiny. should test this with a larger file as well
to see if the issue compounds, though - in the meantime I
could reroute stderr using `2>` to make sure the screen doesn't
explode with cyvcf2 messages

## 17/11/2020

writing a bam preprocessing script - this should:

1. sort file and remove 'undesirable' reads
2. make bam (if not already bam)
3. create `.bai` file
4. use `idxstats` to get expected number of total reads?

before we do that, let's make sure `idxstats` is giving us
the numbers we expect:

```python
import pysam
from tqdm import tqdm
fname = '../recombinant-reads/data/liu-data/bam/SRR5243250.dup.bam'
reader = pysam.AlignmentFile(fname)
counter = dict.fromkeys(['chromosome_' + str(i) for i in range(1, 18)], 0)

for record in tqdm(reader):
    if record.reference_name in counter:
        if record.is_proper_pair:
            counter[record.reference_name] += 1
```

all the counter values seem a bit lower than the idxstats output - 

```
>>> counter
{'chromosome_6': 1311925, 'chromosome_1': 1208415, 'chromosome_3': 1325325, 'chromosome_16': 1151404, 
'chromosome_14': 670269, 'chromosome_2': 1350132, 'chromosome_9': 1138041, 'chromosome_15': 198757, 
'chromosome_7': 910801, 'chromosome_4': 611120, 'chromosome_13': 823079, 'chromosome_11': 490728, 
'chromosome_5': 411561, 'chromosome_17': 1038973, 'chromosome_12': 1401302, 'chromosome_10': 1018450, 'chromosome_8': 701364}

      1 chromosome_1    8033585 1501550 16299
      2 chromosome_2    9223677 1662704 16812
      3 chromosome_3    9219486 1587010 13894
      4 chromosome_4    4091191 755475  9601
      5 chromosome_5    3500558 528501  8320
      6 chromosome_6    9023763 1554211 13822
      7 chromosome_7    6421821 1121690 9996
      8 chromosome_8    5033832 863367  11153
      9 chromosome_9    7956127 1354110 12934
     10 chromosome_10   6576019 1195224 8688
     11 chromosome_11   3826814 613003  9919
     12 chromosome_12   9730733 1648508 13002
     13 chromosome_13   5206065 974336  9449
     14 chromosome_14   4157777 789405  8409
     15 chromosome_15   1922860 270740  5586
     16 chromosome_16   7783580 1358496 11352
     17 chromosome_17   7188315 1217401 10535
```

what happens if we do this without the proper pair check? 

```python
reader = pysam.AlignmentFile(fname)
counter = dict.fromkeys(['chromosome_' + str(i) for i in range(1, 18)], 0)
for record in tqdm(reader):
    if record.reference_name in counter: # this will count unpaired as well...
        counter[record.reference_name] += 1
```

this adds up - each of the values here are paired + unpaired (2nd and 3rd
column) from the `idxstats` output

```python
>>> counter
{'chromosome_6': 1568033, 'chromosome_1': 1517849, 'chromosome_3': 1600904,
 'chromosome_16': 1369848, 'chromosome_14': 797814, 'chromosome_2': 1679516, 
'chromosome_9': 1367044, 'chromosome_15': 276326, 'chromosome_7': 1131686, 
'chromosome_4': 765076, 'chromosome_13': 983785, 'chromosome_11': 622922, 
'chromosome_5': 536821, 'chromosome_17': 1227936, 'chromosome_12': 1661510, 
'chromosome_10': 1203912, 'chromosome_8': 874520}
```

finally, doing this with just paired reads: 

```python
reader = pysam.AlignmentFile(fname)
counter = dict.fromkeys(['chromosome_' + str(i) for i in range(1, 18)], 0)
for record in tqdm(reader):
    if record.reference_name in counter:
        if record.is_paired:
            counter[record.reference_name] += 1
```

wait - this gives the same dict

I'm forgetting now, but wasn't every single read
designated 'paired' here even if the mate was
removed after the fact? in either case, those
reads will be filtered out - so let's instead try
this on an already filtered bam:

```python
fname = 'CC2935_chr1.RG.bam' # only contains chr1
reader = pysam.AlignmentFile(fname)
counter = 0
for record in tqdm(reader):
    if record.reference_name == 'chromosome_1':
        counter += 1
```

this looks great:

```python
>>> for record in tqdm(reader):
  2     if record.reference_name == 'chromosome_1':
  3         counter += 1
7737472it [00:21, 351957.65it/s]
>>> counter
7737472
```

and here's the first bit of the `idxstats` output:

```
      1 chromosome_1    8033585 7737472 0
      2 chromosome_2    9223677 0       0
      3 chromosome_3    9219486 0       0
      4 chromosome_4    4091191 0       0
      5 chromosome_5    3500558 0       0
```

let's get to work on this bam preprocessing script
then - will call it `bamprep.py`

## 23/11/2020

nts - `>>> {l[0]: int(l[2]) for l in [line.split('\t') for line in out.decode('utf-8'
  1 ).split('\n') if len(line) > 1]}`

alright - time to run the command over a varying number
of processes and save the output of `time` to create a
comparison table with

```bash
time python3.5 readcomb/timetest.py
```

need to rerun with larger file - but seems that
user value being pulled in is incorrect - regex
might be screwing up a little here - do another
check in the console tmo

## 2/12/2020

testing out the bam preprocessing script:

```bash
ln -sv ../recombinant-reads/data/liu-data/bam/test_sort_20k.sam .
samtools view test_sort_20k.sam | wc -l
# 19940

samtools view -f 0x2 test_sort_20k.sam | wc -l
# 15837

samtools view -f 0x2 -F 0x100 test_sort_20k.sam | wc -l
# still 15837

samtools view -f 0x2 -F 0x100 -F 0x800 test_sort_20k.sam | wc-l
# 15619

# trying out bam filtering script
time python3.5 readcomb/bamprep.py \
--fname test_sort_20k.sam \
--samtools bin/samtools \
--threads 4 \
--outdir test_out

samtools view test_out/test_sort_20k.sorted.bam | wc -l
# 15619 - and looks correctly sorted in less!
```

trying with alt options just to make sure they work:

```bash
time python3.5 readcomb/bamprep.py \
--fname test_sort_20k.sam \
--samtools bin/samtools \
--threads 4

time python3.5 readcomb/bamprep.py \
--fname test_sort_20k.sam \
--samtools bin/samtools \
--outdir .

time python3.5 readcomb/bamprep.py \
--fname test_sort_20k.sam \
--samtools bin/samtools \
--index_csi
```

all looks good! need to figure out whether the progress 
bar side of things should be made optional though

## 8/3/2021

getting data for ART error rate profiling

initially found this:

```bash
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra31/SRZ/013767/SRR13767223/Juvenile3_R2.fq.gz
```

but this is from a metagenomics paper

fortunately found this S. pombe 2 x 250 dataset - only one of its kind I could find!

the SRA Run Selector is pretty useful to quickly filter for read length (AvgSpotLen)

```bash
cd data/novaseq
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra25/SRR/011509/SRR11785254
ln -sv ~/apps/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump .
time fastq-dump -I --split-files SRR11785254 # creates fwd and rev files
# fastq conversion took 138 min
time gzip SRR11785254_1.fastq # 103 min
time gzip SRR11785254_2.fastq
```
this is a pretty large dataset too, at 15 gigs! 

associated preprint: https://www.biorxiv.org/content/10.1101/2020.03.03.972455v3.full.pdf+html

aligned to ENSEMBL reference genome ASM294v2 - let's see if we can get that downloading too

two main tasks here:

1. create error profile using art profiler
2. get insert sizes after alignment

genome obtained from NCBI: 

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_genomic.fna.gz
```

it's crazy small too - just 4 Mb gzipped, 11 g-unzipped

```
NC_003424.3 5579133 # chr1
NC_003423.3 4539804 # chr2
NC_003421.2 2452883 # chr3
NC_001326.1 19431 # mtDNA
```

## 9/3/2021

today - aligning

using commands from the salt MA project:

```bash
mv -v GCF_000002945.1_ASM294v2_genomic.fna pombe_ref.fna
time bwa mem -t 12 pombe_ref.fna \
SRR11785254_1.fastq.gz \
SRR11785254_2.fastq.gz | \
samtools sort -@4 -T pombe.sorting.tmp -O bam -o pombe.sorted.bam
```

## 11/3/2021

need to fix error with bbmap:

```bash
time ~/.conda/env/work/bin/bbrename.sh in=SRR11785254_1.fastq.gz out=SRR11785254_1_fixed.fastq prefix=SRR11785254
time ~/.conda/env/work/bin/bbrename.sh in=SRR11785254_2.fastq.gz out=SRR11785254_2_fixed.fastq prefix=SRR11785254
```

## 16/3/2021

finally getting around to wrapping this up

first, index the bam:

```bash
samtools index pombe.sorted.bam
```

next, look at insert sizes. looks like a LOT of the reads in the alignment
are soft clipped, so let's look at the 'max' size (e.g. if 250M)

```python
import pysam
from tqdm import tqdm
chrname = 'NC_003424.3'

reader = pysam.AlignmentFile('pombe.sorted.bam', 'rb')
reads = []
for record in tqdm(reader):
    if record.cigarstring == '250M':
        reads.append(record)
    if len(reads) == 10:
        break
```

this yielded one pair, with the id `SRR11785254_79432193`:

```python
r = reads[1:3]
(r[1].reference_start + r[1].query_alignment_length) - r[0].reference_start
# 789
```

the next one I could find (`SRR1178254_28114860`) was 589 bases instead though, so
there's a lot more variation than in the Liu dataset - I think 600 is a fine proxy for ART
in that case

## 20/3/2021

today: running Jimmy's power analysis code, found
[here](https://github.com/ness-lab/recombinant-reads/blob/dev/analysis/jimmy/power_analysis/analysis.ipynb)

going to run this in ptpython though:

```python
import subprocess
import pysam
import numpy.random as r
import re
import ast
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

VCF2FASTA_FILEPATH = 'output/vcf2fasta.fasta'
ANALYSIS_SCRIPT = 'analysis.py'

# run analysis.py to create fasta sequence with phase change
# and generate reads with art_illumina

COVERAGE = 1
PERCENTAGE = 0.0004
BASE_COUNT = 1000000
ANALYSIS_OUTPUT = 'output/unfiltered'
FILTERED_OUTPUT = 'output/filtered'
window_counts = open('output/window_counts.tsv', 'w')

iteration = 1

for i in range(5):

    subprocess.run(map(str, ['python3', 'analysis.py',
                    '-f', VCF2FASTA_FILEPATH,
                    '-e', BASE_COUNT,
                    '-c', COVERAGE,
                    '-p', PERCENTAGE,
                    '-o', ANALYSIS_OUTPUT + str(iteration)]))

    # read last line of output/phase_change_log.txt 
    # that analysis.py writes to get list of phase changes

    phase_changes = subprocess.check_output(
    ['tail', '-1', 'output/phase_change_log.txt']).decode('utf-8').strip()

    phase_changes = ast.literal_eval(phase_changes)

    # run readcomb on output of art_illumina
    subprocess.run(['readcomb-filter',
                    '-b', ANALYSIS_OUTPUT + str(iteration) + '.sam',
                    '-v', '../../../data/liu-data/vcf/filtered_full.vcf.gz',
                    '-p', '4',
                    '-o', FILTERED_OUTPUT + str(iteration) + '.sam'])
    ###

    # sort and index the readcomb-filtered and unfiltered sam files

    with open(ANALYSIS_OUTPUT + str(iteration) + 'sorted.bam', 'w') as f:
        subprocess.call(['samtools', 'sort', ANALYSIS_OUTPUT + str(iteration) + '.sam'], stdout=f)

    subprocess.run(['samtools', 'index', ANALYSIS_OUTPUT + str(iteration) + 'sorted.bam'])

    with open(FILTERED_OUTPUT + str(iteration) + 'sorted.bam', 'w') as f:
        subprocess.call(['samtools', 'sort', FILTERED_OUTPUT + str(iteration) + '.sam'], stdout=f)

    subprocess.run(['samtools', 'index', FILTERED_OUTPUT + str(iteration) + 'sorted.bam'])

    # read in indexed files and fetch reads inside windows
    #### 
    unfiltered = pysam.AlignmentFile(ANALYSIS_OUTPUT + str(iteration) + 'sorted.bam', 'r')
    filtered = pysam.AlignmentFile(FILTERED_OUTPUT + str(iteration) + 'sorted.bam', 'r')

    for i in range(0, BASE_COUNT, 1000):
        filtered_count = len(list(filtered.fetch('chromosome_1', i, i+2000)))
        unfiltered_count = len(list(unfiltered.fetch('chromosome_1', i, i+2000)))

        window_counts.write('\t'.join(map(str, [
            iteration,
            i,
            i + 2000,
            filtered_count,
            unfiltered_count,
            list(filter(lambda x: x>i and x<i+2000, phase_changes))
        ])) + '\n')
        
    iteration += 1
```

some issues -

- needed to replace `base_count` with `BASE_COUNT` in the final loop
- needed to fix `filtered_full.vcf.gz` path
- needed to add `.sam` as output to the first run of `readcomb-filter`
- seems like the expected behaviour if no bai files is still a little off
    - was getting this error:

```bash
Bai file not found, continuing without progress bars
Creating processes
Processes created
Process Counter-1:
Traceback (most recent call last):
  File "/home/hasans11/.conda/env/work/lib/python3.6/multiprocessing/process.py", line 258, in _bootstrap
    self.run()
  File "/home/hasans11/.conda/env/work/bin/filter.py", line 343, in run
    self.progress.close()
AttributeError: 'Counter' object has no attribute 'progress'
```

process still ran correctly, to be clear, but need to look at this bug at some pointo

next up - generating the 400 fastas:

```python
args = {}
args['end'] = 1100000
args['start'] = 1000000
args['percent'] = 0.00004
args['fasta'] = 'output/vcf2fasta.fasta'

iteration = 1

# parse results from vcf2fasta
fasta_file = list(SeqIO.parse(args['fasta'], 'fasta'))

# tsv for phase changes
phase_change_tsv = open('output400/phase_changes.tsv', 'w')

# regex to identify correct fasta reference_name
chrom = re.search(r'chromosome_[0-9]{1,2}', str(fasta_file[0].id)).group()

fastas = []

for i in tqdm(range(400)):

    events = r.poisson((args['end'] - args['start']) * args['percent'])
    phase_changes = sorted(list(r.choice((args['end']- args['start']), size=events)))
    
    phase_change_tsv.write('sequence' + str(iteration) + '\t' + str(phase_changes) + '\n')

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

    recomb_record = SeqRecord(Seq(''.join(recomb_seq)), description=chrom + '_' + str(iteration),
        name='sequence_' + str(iteration), id=chrom, dbxrefs=fasta_file[0].dbxrefs, 
        features=fasta_file[0].features, annotations=fasta_file[0].annotations, 
        letter_annotations=fasta_file[0].letter_annotations
    )

    fastas.append(recomb_record)
    iteration += 1

SeqIO.write(fastas, 'output400/400_generated_sequences.fasta', 'fasta')
phase_change_tsv.close()
```

and finally, the actual generation of sams and such:

```python
from classification import *

# 400 coverage with 400 sequences

VCF_FILEPATH = '../../../data/liu-data/vcf/filtered_full.vcf.gz'
window_counts_tsv = open('output400/window_counts.tsv', 'w')
SEQUENCE400 = 'output400/400_generated_sequences.fasta'

COVERAGE = 5
INSERT = 100
ART_OUTPUT = 'output400/unfiltered'
FILTERED_OUTPUT = 'output400/filtered'

iteration = 1

for i in range(1):
    
    # run art_illumina on 400 sequences with coverage 400
    
    subprocess.run(['art_illumina',
                    '-f', str(COVERAGE),
                    '-l', '250',
                    '-ss', 'MSv1',
                    '-i', SEQUENCE400,
                    '-m', str(500 + INSERT),
                    '-s', '1',
                    '-o', ART_OUTPUT + str(iteration),
                    '-M', '-p', '-sam', '-na', '-q'])

    # run readcomb on 400 sequence output generated by art_illumina

    subprocess.run(['readcomb-filter',
                    '-b', ART_OUTPUT + str(iteration) + '.sam',
                    '-v', VCF_FILEPATH,
                    '-p', '4',
                    '-o', FILTERED_OUTPUT + str(iteration) + '.sam'])
    

    # read in sam files and 
    
    unfiltered_pairs = pairs_creation(ART_OUTPUT + str(iteration) + '.sam', VCF_FILEPATH)
    filtered_pairs = pairs_creation(FILTERED_OUTPUT + str(iteration) + '.sam', VCF_FILEPATH)
    
    # call get_midpoint on all of them, this will take quite a long time so use tqdm
    
    #[[filtered, unfiltered], ...]
for i in range(1):
    window_counts = [[0,0] for i in range(100)]
    
    for pair in tqdm(filtered_pairs):
        midpoint = pair.get_midpoint()
        window = int(midpoint // 1000)
        window_counts[window][0] += 1
    
    for pair in tqdm(unfiltered_pairs):
        midpoint = pair.get_midpoint()
        window = int(midpoint // 1000)
        window_counts[window][1] += 1
        
    for i in range(len(window_counts)):
        window_counts_tsv.write('\t'.join(map(str, [
            iteration,
            i * 1000,
            (i + 1) * 1000,
            window_counts[i][0],
            window_counts[i][1]
        ])) + '\n'
            
        )
        
    iteration += 1
```

## 22/3/2021

so this was breaking specifically because I forgot Jimmy had sent updated +
fixed code! 

removing all the 400 coverage files, reinstalling readcomb, and then
running the analysis script

```bash
pip install --upgrade readcomb
python 400coverage.py
```

had to symlink `art_illumina` into `~/bin` instead of adding the
ART folder to my PATH

uh oh - need to regenerate the fasta as well - might have been preemptive to delete
that folder! 

```python
from Bio import SeqIO
import numpy.random as r
from tqdm import tqdm
import re
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

args = {}
args['end'] = 200000
args['start'] = 100000
args['percent'] = 0.00004
args['fasta'] = 'data/vcf2fasta.fasta'

iteration = 1

# parse results from vcf2fasta
fasta_file = list(SeqIO.parse(args['fasta'], 'fasta'))

# tsv for phase changes
phase_change_tsv = open('output/phase_changes.tsv', 'w')

# regex to identify correct fasta reference_name
chrom = re.search(r'chromosome_[0-9]{1,2}', str(fasta_file[0].id)).group()

fastas = []

for i in tqdm(range(400)):

    events = r.poisson((args['end'] - args['start']) * args['percent'])
    phase_changes = sorted(list(r.choice((args['end']- args['start']), size=events)))
    
    phase_change_tsv.write('sequence' + str(iteration) + '\t' + str(phase_changes) + '\n')

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

    recomb_record = SeqRecord(Seq(''.join(recomb_seq)), description=chrom + '_' + str(iteration),
        name='sequence_' + str(iteration), id=chrom + '_' + str(iteration), dbxrefs=fasta_file[0].dbxrefs, 
        features=fasta_file[0].features, annotations=fasta_file[0].annotations, 
        letter_annotations=fasta_file[0].letter_annotations
    )

    fastas.append(recomb_record)
    iteration += 1

SeqIO.write(fastas, 'output/400_generated_sequences.fasta', 'fasta')
phase_change_tsv.close()
```

## 26/3/2021

first pass: seems the reads in the filtered sam are almost 20x the 'phase change rate':

```
wc -l *sam
       495 filtered1.sam
       497 filtered2.sam
       471 filtered3.sam
       493 filtered4.sam
       473 filtered5.sam
    637015 unfiltered1.sam
    637053 unfiltered2.sam
    638019 unfiltered3.sam
    637089 unfiltered4.sam
    637197 unfiltered5.sam

```

which is roughly 0.00077 - likely cause certain phase changes have multiple reads spanning them

## 31/3/2021

running this again - the previous stuff worked, but perhaps if I set the id in the generation
of `400_sequences.fasta` above to include the iteration count I can trace reads back to their
'original' sequence - plus I have to add `SQ:chromosome_1 LN:8033585` to the top of each sam anyways

split `400coverage.py` into that and `400coverage-2.py` since the steps post sam generation
require for the header to be modified

set coverage to 100 and it still took 55 minutes - woof! now for the fun part

update: never mind - the reference names are still `chromosome_1_1` etc - which means
readcomb is not able to search for phase changes in the VCF at all - going to have
to stick with the original files then

have to repeat them with 100 coverage though - was run with 5 earlier! 

## 1/4/2021

fixing the sam headers - can't just manually do this in vim this time,
so will have to pull from the earlier files

```bash
head -n 2 output-old/unfiltered1.sam > corrected_header
tail -n +3 output/unfiltered1.sam | cat corrected_header - > output/unfiltered1.corrected.sam

# looks good - repeating
for i in {2..5}; do
    tail -n +3 output/unfiltered${i}.sam | cat corrected_header - > output/unfiltered${i}.corrected.sam
done

rm -v output/unfiltered?.sam

for i in {1..5}; do
    mv -v output/unfiltered${i}.corrected.sam output/unfiltered${i}.sam
done

# and now we're ready to go:
time python 400coverage-2.py
```

## 3/4/2021

looks good - took 71 minutes to complete, which isn't bad

things to do:

- how many phase change reads did we actually end up getting?
- use `window_counts.tsv` to make an 'expected' landscape and then match up with filtered reads 

first, let's figure some things out re: art

```
# in art examples folder
mkdir -p tests
ls testSeq.fa # 2 seqs, 7207 length of first, 3057 second

# sim 4 in example script
# PE of length 100, mean fragment size 500 and SD 10
# coverage of 1
../art_illumina -ss HS20 -i testSeq.fa -o tests/paired_end_sep -l 100 -f 1 -p -m 500 -s 10 -sp -sam

# coverage 2
../art_illumina -ss HS20 -i testSeq.fa -o tests/paired_end_sep_2 -l 100 -f 2 -p -m 500 -s 10 -sp -sam
```

looks like art spans the entire sequence each time - here's coverage 1:

```python
>>> first = [r for r in pysam.AlignmentFile('paired_end_sep.sam', 'r')]
>>> second = [r for r in pysam.AlignmentFile('paired_end_sep_2.sam', 'r')]
>>> len(first)
102

>>> len(second)
204

>>> first[0].reference_name
'seq1'

>>> set([r.reference_name for r in first])
{'seq2', 'seq1'}

>>> seq2_first = [r for r in first if r.reference_name == 'seq2']
>>> sorted([r.reference_start for r in seq2_first])
'''
[419, 554, 739, 744, 836, 918, 961, 1136, 1142, 1155, 1314, 1335, 1345, 1538,
1619, 1710, 1722, 1760, 1976, 2009, 2115, 2188, 2266, 2281, 2351, 2389, 2600,
2662, 2701, 2743]
'''

>>> seq2_second = [r for r in second if r.reference_name == 'seq2']
>>> sorted([r.reference_start for r in seq2_second])
'''
[5, 50, 70, 126, 244, 248, 256, 387, 416, 426, 426, 442, 474, 531, 560, 618,
650, 659, 666, 683, 799, 808, 821, 833, 903, 949, 1009, 1036, 1090, 1100, 1108,
1165, 1201, 1281, 1300, 1355, 1404, 1407, 1490, 1522, 1526, 1548, 1672, 1691,
1788, 1797, 1933, 1950, 2045, 2068, 2120, 2188, 2273, 2351, 2446, 2519, 2553,
2589, 2676, 2942]
'''
```

31 reads with coverage 1 and fragment size 500, 61 with coverage 2

what happens if we reduce fragment size to 250? 

```bash
>>> short = [r for r in pysam.AlignmentFile('paired_end_sep_short.sam', 'r')]
>>> seq2_short = [r for r in short if r.reference_name == 'seq2']
>>> sorted([r.reference_start for r in seq2_short])
'''
[10, 175, 236, 357, 398, 404, 429, 503, 541, 550, 586, 708, 807, 851, 970, 994,
1174, 1185, 1192, 1316, 1331, 1353, 1810, 1868, 1977, 2015, 2150, 2290, 2681,
2851]
'''
```

## 4/4/2021

need to generate a new fasta with separate IDs for its names and then run art
with coverage 1 to understand how many reads are generated per sample

used the fasta code above and with `id` set to `chrom + '_' + str(iteration)` into
`output-test`

art time:

```bash
art_illumina -f 1 -l 250 -ss MSv1 -i output-test/400_generated_sequences.fasta \
-m 600 -s 1 -o output-test/unfiltered -M -p -sam -na -q
```

400 reads per sample - makes sense given the original sequences are 100 kb long!
(100000/250 = 400)

## 6/4/2021

now - figuring out what we need to set coverage to to get an accurate readcomb estimation

some numbers:

- sequence is 100000 bases long - might bump this up to first 200 kb of chromosome 1 actually
- read length is 250 - that's 400 expected reads per 'coverage', or 800 for new length
- we have 400 individuals in the dataset

our actual sequencing will likely be at roughly 100x, which means that on average,
a given locus will have 100 reads from 100 different individuals spanning it

if `-f` is set to 0.25, we get 200 reads from each indiv -> 80000 total reads

I think this makes sense? let's just give it a go for now

this time, the fasta should have iterations as part of the sequence name (since we
can edit these later) and we also need to add in a 200 kb snippet of chromosome 1
proper for the SAM header (the reads themselves can be removed later on) 

prepping fasta in a new `output` folder after deleting the contents of the original
(yes they took a while to run but those files are useless now)

```bash
# in power_analysis
mkdir output

```

generating the fasta -

```python
from Bio import SeqIO
import numpy.random as r
from tqdm import tqdm
import re
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

args = {}
args['end'] = 200000
args['start'] = 0
args['percent'] = 0.00004
args['fasta'] = 'data/vcf2fasta.fasta'

iteration = 1

# parse results from vcf2fasta
fasta_file = list(SeqIO.parse(args['fasta'], 'fasta'))

# tsv for phase changes
phase_change_tsv = open('output/phase_changes.tsv', 'w')

# regex to identify correct fasta reference_name
chrom = re.search(r'chromosome_[0-9]{1,2}', str(fasta_file[0].id)).group()

fastas = []
fastas.append(
    SeqRecord(Seq(str(fasta_file[0].seq))[args['start']:args['end']], 
        id='chromosome_1', description='chromosome_1')
)

for i in tqdm(range(400)):

    events = r.poisson((args['end'] - args['start']) * args['percent'])
    phase_changes = sorted(list(r.choice((args['end']- args['start']), size=events)))
    
    phase_change_tsv.write('sequence' + str(iteration) + '\t' + str(phase_changes) + '\n')

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

    recomb_record = SeqRecord(Seq(''.join(recomb_seq)[:-1]), description=chrom + '_' + str(iteration),
        name='sequence_' + str(iteration), id=chrom + '_' + str(iteration), dbxrefs=fasta_file[0].dbxrefs, 
        features=fasta_file[0].features, annotations=fasta_file[0].annotations, 
        letter_annotations=fasta_file[0].letter_annotations
    )

    fastas.append(recomb_record)
    iteration += 1

SeqIO.write(fastas, 'output/400_generated_sequences.fasta', 'fasta')
phase_change_tsv.close()
```

needed to add `[:-1]` to the final record construction cause the code was catching one
extra base

will also need to remove the 'vanilla' chromosome 1 reads when fixing the reference names

looks good - and now for `400coverage.py`, fixing the sam reference names,
and then finally `400coverage-2.py`, all with coverage set to 0.25

`fix_names.py`:

```python
# fixing the reference names and removing chr1
import sys
import pysam
from tqdm import tqdm

fname = sys.argv[-2]
outname = sys.argv[-1]

reader = pysam.AlignmentFile(fname, 'r')
writer = pysam.AlignmentFile(outname, 'wh', template=reader)

for record in tqdm(reader):
    if record.reference_name == 'chromosome_1':
        continue
    else:
        record.reference_name = 'chromosome_1'
        record.next_reference_name = 'chromosome_1'
        writer.write(record)

writer.close()
```

followed by:

```bash
for i in {2..5}; do
    python fix_names.py output/unfiltered${i}.sam output/unfiltered${i}.corrected.sam;
done

# 2 to 5 since I tested this on the first one already

for i in {1..5}; do 
    mv -v output/unfiltered${i}.sam output/unfiltered${i}.old.sam
    mv -v output/unfiltered${i}.corrected.sam output/unfiltered${i}.sam
done

time python 400coverage-2.py
```

all done - at first pass, lots more phase change reads than I expected, about ~800
in each filtered file - e.g. 1 percent of all reads 

things to check tomorrow:

1. how many phase changes for each sequence were found? how many missed? how many false calls? 
2. did the coverage calculation 'work'? (is average depth roughly 100?)
3. how many phase changes in a given indiv showed up in more than one read?
    - e.g. how close is this to reality - ideally, no indiv should have more than one read at a locus

## 7/4/2021

doing the above analysis in a jupyter notebook cause plotting might be helpful here


