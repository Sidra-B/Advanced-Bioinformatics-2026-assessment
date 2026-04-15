#!/bin/bash

# Define inputs as Positional Parameters ($1 and $2)
R1=$1
R2=$2

#Project Paths
REF=~/ngs_workflow/dnaseq/data/reference/hg19.fa
RESULTS=~/ngs_workflow/dnaseq/results
DATA=~/ngs_workflow/dnaseq/data
ANNOVAR_DIR=~/annovar

echo "--- STARTING PIPELINE ---"

# Step 2.2: QC and Trimming
fastp -i $R1 -I $R2 -o ${R1}_trimmed.fq.gz -O ${R2}_trimmed.fq.gz

# Step 2.3: Alignment & Filtering
bwa mem -t 4 -R '@RG\tID:WES01\tSM:WES01\tPL:ILLUMINA' $REF ${R1}_trimmed.fq.gz ${R2}_trimmed.fq.gz > $DATA/aligned_data/NGS0001.sam
samtools sort -o $DATA/aligned_data/NGS0001_sorted.bam $DATA/aligned_data/NGS0001.sam
picard MarkDuplicates I=$DATA/aligned_data/NGS0001_sorted.bam O=$DATA/aligned_data/NGS0001_marked.bam M=$RESULTS/marked_dup_metrics.txt
samtools view -F 1796 -q 20 -b -o $DATA/aligned_data/NGS0001_final.bam $DATA/aligned_data/NGS0001_marked.bam
samtools index $DATA/aligned_data/NGS0001_final.bam

# Step 2.4: Variant Calling (FreeBayes)
freebayes -f $REF $DATA/aligned_data/NGS0001_final.bam > $RESULTS/NGS0001_variants.vcf

# Step 2.5: Annotation (ANNOVAR)
cd $ANNOVAR_DIR
./table_annovar.pl $RESULTS/NGS0001_variants.vcf humandb/ -buildver hg19 -out $RESULTS/NGS0001_annovar -remove -protocol refGene -operation g -nastring . -csvout

echo "--- PIPELINE COMPLETE ---"
