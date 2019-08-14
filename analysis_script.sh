#!/bin/bash
#PBS -l select=1:ncpus=10:mem=40gb
#PBS -l walltime=48:00:00
#PBS -N analysis_combined_scripts
#PBS -m abe
#PBS -M guyf0601@gmail.com
#PBS -q bep

cd $PBS_O_WORKDIR

### step 1: cut adapters and filter
## run fastqc
module load java/8.0_161 Trimmomatic/0.38

java -jar /software/Trimmomatic/0.38/trimmomatic-0.38.jar PE -trimlog logfile \
-phred33 -threads 4 \
$PBS_O_WORKDIR/NexteraA_1.fastq $PBS_O_WORKDIR/NexteraA_2.fastq \
$PBS_O_WORKDIR/NexteraA_1.trimmed.fastq $PBS_O_WORKDIR/NexteraA_1.unpaired.fasq \
$PBS_O_WORKDIR/NexteraA_2.trimmed.fastq $PBS_O_WORKDIR/NexteraA_2.unpaired.fasq \
ILLUMINACLIP:/software/Trimmomatic/0.38/adapters/NexteraPE-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar /software/Trimmomatic/0.38/trimmomatic-0.38.jar PE -trimlog logfile \
-phred33 -threads 4 \
$PBS_O_WORKDIR/NexteraB_1.fastq $PBS_O_WORKDIR/NexteraB_2.fastq \
$PBS_O_WORKDIR/NexteraB_1.trimmed.fastq $PBS_O_WORKDIR/NexteraB_1.unpaired.fasq \
$PBS_O_WORKDIR/NexteraB_2.trimmed.fastq $PBS_O_WORKDIR/NexteraB_2.unpaired.fasq \
ILLUMINACLIP:/software/Trimmomatic/0.38/adapters/NexteraPE-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
## run fastqc

### step:2 map to sam file, then sort to bam file
module load bwa/0.7.17 samtools/1.8

## use wget to download the GRCh38 data file
## use bwa to create its index
#bwa index $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \

# for NexteraA
bwa mem -t 4 $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
$PBS_O_WORKDIR/NexteraA_1.trimmed.fastq \
$PBS_O_WORKDIR/NexteraA_2.trimmed.fastq > $PBS_O_WORKDIR/NexteraA.mapped.sam

samtools view -S -b NexteraA.mapped.sam > NexteraA.mapped.bam
samtools sort -o NexteraA.mapped.sorted.bam NexteraA.mapped.bam
samtools index NexteraA.mapped.sorted.bam

# for NexteraB
bwa mem -t 4 $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
$PBS_O_WORKDIR/NexteraB_1.trimmed.fastq \
$PBS_O_WORKDIR/NexteraB_2.trimmed.fastq > $PBS_O_WORKDIR/NexteraB.mapped.sam

## ncpus is required as 3.98, so must specify ncpus>=4
samtools view -S -b NexteraB.mapped.sam > NexteraB.mapped.bam
samtools sort -o NexteraB.mapped.sorted.bam NexteraB.mapped.bam
samtools index NexteraB.mapped.sorted.bam


### step 3: mark and remove duplicates
module load Picard/2.18.9

java -jar /software/Picard/2.18.9/picard.jar MarkDuplicates \
  I=NexteraA.mapped.sorted.bam \
  O=NexteraA.marked_duplicates.bam \
  M=NexteraA.marked_dup_metrics.txt \
  REMOVE_DUPLICATES=true

java -jar /software/Picard/2.18.9/picard.jar MarkDuplicates \
  I=NexteraB.mapped.sorted.bam \
  O=NexteraB.marked_duplicates.bam \
  M=NexteraB.marked_dup_metrics.txt \
  REMOVE_DUPLICATES=true

### step 4; base quality recalibration

## wget .vcf file from database and its coresponding index file

## create .dict file
#java -jar /software/Picard/2.18.9/picard.jar CreateSequenceDictionary \
#      REFERENCE= $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
#      OUTPUT= $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.dict

## create .fai file
#samtools faidx $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa

## add read groups in marked_duplicates bam file
java -jar /software/Picard/2.18.9/picard.jar AddOrReplaceReadGroups \
       I=NexteraA.marked_duplicates.bam \
       O=NexteraA.marked_duplicates.add.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=20

java -jar /software/Picard/2.18.9/picard.jar AddOrReplaceReadGroups \
       I=NexteraB.marked_duplicates.bam \
       O=NexteraB.marked_duplicates.add.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=20

## run BaseRecalibrator & generate recal_data.table file
java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar BaseRecalibrator \
        -I $PBS_O_WORKDIR/NexteraA.marked_duplicates.add.bam \
        -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
        --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
        -O $PBS_O_WORKDIR/NexteraA.recal_data.table

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar BaseRecalibrator \
        -I $PBS_O_WORKDIR/NexteraA.marked_duplicates.add.bam \
        -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
        --known-sites $PBS_O_WORKDIR/Homo_sapiens_assembly38.dbsnp138.vcf \
        -O $PBS_O_WORKDIR/NexteraA.recal_data.table

## apply recalibration
java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar ApplyBQSR \
        -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
        -I $PBS_O_WORKDIR/NexteraA.marked_duplicates.add.bam \
        --bqsr-recal-file $PBS_O_WORKDIR/NexteraA.recal_data.table \
        -O $PBS_O_WORKDIR/NexteraA.marked_duplicates.add.apply.bam

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar ApplyBQSR \
        -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
        -I $PBS_O_WORKDIR/NexteraB.marked_duplicates.add.bam \
        --bqsr-recal-file $PBS_O_WORKDIR/NexteraA.recal_data.table \
        -O $PBS_O_WORKDIR/NexteraB.marked_duplicates.add.apply.bam


## step 5, generate .vcf .gvcf file
## be careful when specify the -R reference file,
## different test/dataset requires different reference .fa file 
module load java/8.0_161 GenomeAnalysisTK/4.1.2.0

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar HaplotypeCaller \
	-R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
	-I $PBS_O_WORKDIR/NexteraA.marked_duplicates.add.apply.bam \
	-O NexteraA.output.g.vcf.gz \
	-ERC GVCF

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar HaplotypeCaller \
  -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
  -I $PBS_O_WORKDIR/NexteraA.marked_duplicates.add.apply.bam \
  -O NexteraA.output.g.vcf.gz \
  -ERC GVCF

## generate .vcf file + .bam file
java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar HaplotypeCaller \
  -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
  -I $PBS_O_WORKDIR/NexteraA.marked_duplicates.add.apply.bam \
  -O NexteraA.output.vcf.gz \
  -bamout NexteraA.bamout.bam

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar HaplotypeCaller \
  -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
  -I $PBS_O_WORKDIR/NexteraB.marked_duplicates.add.apply.bam \
  -O NexteraB.output.vcf.gz \
  -bamout NexteraB.bamout.bam
