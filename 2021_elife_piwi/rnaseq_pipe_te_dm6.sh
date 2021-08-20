#!/bin/bash

## see raw scripts: https://raw.githubusercontent.com/mhf27/hannon_roo_fabry2021/main/rnaseq_ver8_slurm_te_first_dm6_sec.sh

## RNAseq analysis

## global variables
te_index='~/data/genome/dm6/STAR_index/dm6_transposon'
te_gtf='~/data/genome/dm6/dm6_transposon/dm6_transposon.gtf'
dm6_index='~/data/genome/dm6/STAR_index/dm6'
dm6_gene_gtf='~/data/genome/dm6/annotation_and_repeats/Drosophila_melanogaster.BDGP6.92.gtf'

# input: RNA-Seq*fq.gz
echo 'Start' ${1}
filename=$(basename ${1%%.fq.gz})

## trimming reads first, last base
zcat $1 | fastx_trimmer -Q33 -f 2 -l 49 -i - -o $filename.trimmed

#2 mismatches per read, add --alignIntronMax 1 to suppress spliced reads, add --alignEndsType EndToEnd to account for difficult to map indels (makes it more relatable to bowtie)
STAR --genomeDir $te_index --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1000 --outMultimapperOrder Random --outSAMmultNmax 1 --alignIntronMax 100 --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.04 --outReadsUnmapped Fastx --readFilesIn $filename.trimmed --runThreadN 6 --outFileNamePrefix $filename.te.
echo 'STAR dm6 done' $filename;
samtools index $filename.te.Aligned.sortedByCoord.out.bam 
echo 'Index bai done' $filename
STAR --genomeDir $dm6_index --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outMultimapperOrder Random --outSAMmultNmax 1 --alignIntronMax 100 --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.04 --readFilesIn $filename.te.Unmapped.out.mate1 --runThreadN 6 --outFileNamePrefix $filename.dm6.
samtools index $filename.dm6.Aligned.sortedByCoord.out.bam
echo 'Index bai done' $filename

# 16 = antisense, 0 = sense transcripts mapping to TE genome
# Use sorted bam file as input
# Generates list with read number for each TE

#Split te reads into sense (0) and antisense (16)
samtools view -f 0x10 -b $filename.te.Aligned.sortedByCoord.out.bam > $filename.te.0.bam;
samtools index $filename.te.0.bam;
samtools idxstats $filename.te.0.bam | cut -f 1,3 > $filename.te.0.chrom_reads.txt;

samtools view -F 0x10 -b $filename.te.Aligned.sortedByCoord.out.bam > $filename.te.16.bam;
samtools index $filename.te.16.bam;
samtools idxstats $filename.te.16.bam | cut -f 1,3 > $filename.te.16.chrom_reads.txt;

#Split dm6 reads into sense (0) and antisense (16)
samtools view -f 0x10 -b $filename.dm6.Aligned.sortedByCoord.out.bam > $filename.dm6.0.bam;
samtools index $filename.dm6.0.bam;
samtools idxstats $filename.dm6.0.bam | cut -f 1,3 > $filename.dm6.0.chrom_reads.txt;

samtools view -F 0x10 -b $filename.dm6.Aligned.sortedByCoord.out.bam > $filename.dm6.16.bam;
samtools index $filename.dm6.16.bam;
samtools idxstats $filename.dm6.16.bam | cut -f 1,3 > $filename.dm6.16.chrom_reads.txt;

# Normalisation with scaling factor!!!
echo 'Start Normalization'

#Calculate total reads
var1=$(samtools view -c -F 260 $filename.dm6.Aligned.sortedByCoord.out.bam)
var2=$(samtools view -c -F 260 $filename.te.Aligned.sortedByCoord.out.bam)
var3=$((var1+var2))

#calculate reads for condition (te.0) / total reads
var4=$(samtools view -c -F 260 $filename.te.0.bam)
var5=$(bc <<<"scale=10; $var4 / $var3")
echo $var5 > $filename.te.0.scaling.factor.txt
bamCoverage -b $filename.te.0.bam -bs 10 -p 4 --scaleFactor $var5 --normalizeUsing CPM -o $filename.te.0.cpm.norm.bw

#calculate reads for condition (te.16) / total reads
var4=$(samtools view -c -F 260 $filename.te.16.bam)
var5=$(bc <<<"scale=10; $var4 / $var3")
echo $var5 > $filename.te.16.scaling.factor.txt
bamCoverage -b $filename.te.16.bam -bs 10 -p 4 --scaleFactor $var5 --normalizeUsing CPM -o $filename.te.16.cpm.norm.bw

#calculate reads for condition (dm6.0) / total reads
var4=$(samtools view -c -F 260 $filename.dm6.0.bam)
var5=$(bc <<<"scale=10; $var4 / $var3")
echo $var5 > $filename.dm6.0.scaling.factor.txt
bamCoverage -b $filename.dm6.0.bam -bs 10 -p 4 --scaleFactor $var5 --normalizeUsing CPM -o $filename.dm6.0.cpm.norm.bw

#calculate reads for condition (dm6.16) / total reads
var4=$(samtools view -c -F 260 $filename.dm6.16.bam)
var5=$(bc <<<"scale=10; $var4 / $var3")
echo $var5 > $filename.dm6.16.scaling.factor.txt
bamCoverage -b $filename.dm6.16.bam -bs 10 -p 4 --scaleFactor $var5 --normalizeUsing CPM -o $filename.dm6.16.cpm.norm.bw

#htseq on te and dm6
# fix gene GTF, add chr to chr name
awk '{if($1!~"#"){$0="chr"$0} print $0}' $dm6_gene_gtf > dm6.ensembl.gtf
htseq-count -s reverse -f bam -i gene_name $filename.te.Aligned.sortedByCoord.out.bam $te_gtf > $filename.te.count.htseq
htseq-count -s reverse -f bam -i gene_name $filename.dm6.Aligned.sortedByCoord.out.bam dm6.ensembl.gtf > $filename.gene.count.htseq
# merge te and gene counts
cat <(head -n -5 $filename.te.count.htseq) <(head -n -5 $filename.gene.count.htseq) > $filename.te_gene.count.htseq

echo "Done: $filename"
