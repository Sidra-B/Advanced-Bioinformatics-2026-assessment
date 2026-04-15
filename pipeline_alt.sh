#!/bin/bash

# Define inputs as Positional Parameters ($1 and $2)
R1=$1
R2=$2

# Project Paths
REF=~/ngs_workflow/dnaseq/data/reference/hg19.fa
RESULTS=~/ngs_workflow/dnaseq/results
DATA=~/ngs_workflow/dnaseq/data
ANNOVAR_DIR=~/annovar

echo "--- STARTING ALTERNATIVE PIPELINE (BCFTOOLS) ---"

# Step 2.2: Trimming (fastp)
echo "Step 1: Trimming..."
fastp -i $R1 -I $R2 -o ${R1}_trimmed.fq.gz -O ${R2}_trimmed.fq.gz

# Step 2.3: Alignment & Filtering (BWA & Samtools)
echo "Step 2: Alignment..."
bwa mem -t 4 -R '@RG\tID:WES01\tSM:WES01\tPL:ILLUMINA' $REF ${R1}_trimmed.fq.gz ${R2}_trimmed.fq.gz > $DATA/aligned_data/NGS0001.sam

echo "Step 3: Sorting and Marking Duplicates..."
samtools sort -o $DATA/aligned_data/NGS0001_sorted.bam $DATA/aligned_data/NGS0001.sam

picard MarkDuplicates I=$DATA/aligned_data/NGS0001_sorted.bam O=$DATA/aligned_data/NGS0001_marked.bam M=$RESULTS/marked_dup_metrics.txt

echo "Step 4: Filtering high quality reads..."
samtools view -F 1796 -q 20 -b -o $DATA/aligned_data/NGS0001_final.bam $DATA/aligned_data/NGS0001_marked.bam
samtools index $DATA/aligned_data/NGS0001_final.bam

# Step 2.6: ALTERNATIVE Variant Calling (Replacing FreeBayes with BCFTools)
echo "Step 5: Calling variants with BCFTools..."
# mpileup generates likelihoods, call performs the actual calling
bcftools mpileup -f $REF $DATA/aligned_data/NGS0001_final.bam | \
    bcftools call -mv -Ob -o $RESULTS/NGS0001_alternative_bcftools.bcf

# Convert the binary BCF file to VCF for annotation
bcftools view $RESULTS/NGS0001_alternative_bcftools.bcf > $RESULTS/NGS0001_variants.vcf

# Step 2.5: Annotation (ANNOVAR)
echo "Step 6: Annotating variants..."
cd $ANNOVAR_DIR
./table_annovar.pl $RESULTS/NGS0001_variants.vcf humandb/ -buildver hg19 -out $RESULTS/NGS0001_annovar -remove -protocol refGene -operation g -nastring . -csvout

echo "--- ALTERNATIVE PIPELINE COMPLETE ---"
