#!/bin/bash
#PBS -l select=1:ncpus=10:mem=40gb
#PBS -l walltime=48:00:00
#PBS -N analysis_combined_scripts
#PBS -m abe
#PBS -M youremailaddress
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
java -jar /software/Picard/2.18.9/picard.jar CreateSequenceDictionary \
      REFERENCE= $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
      OUTPUT= $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.dict

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
## use the tool HaplotypeCaller
## the memory required for generate .gvcf better set higher
## in this case, use mem = 20g (run time error when mem = 10g)
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

## step 6: follow the GATK best practise, generate database
### follow the guide in GATK, set color by tag 'HC', however, the visualization
### did not make sense unless the INDEL been defined

mkdir database

## here run into issue...
## how to specify the intervals
java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar \
  --java-options "-Xmx4g" GenomicsDBImport \
    -V $PBS_O_WORKDIR/NexteraA.output.g.vcf.gz \
    #-V $PBS_O_WORKDIR/NexteraB.output.g.vcf.gz \
    --genomicsdb-workspace-path database/NexteraA_database \
    -L 20

java  -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar \
  --java-options "-Xmx4g -Xms4g" GenotypeGVCFs \
  -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
  -V $PBS_O_WORKDIR/database/NexteraA_database \
  -O output.NexteraA.vcf.gz

### maybe skip this step, use the previous generated .vcf file to run VariantRecalibrator


## step 7: VariantRecalibrator
### notice the difference in syntax
### depending on the operation environment, sometimes manualy load modules is required


cd $PBS_O_WORKDIR

module load java/8.0_161 GenomeAnalysisTK/4.1.2.0
module load R/3.4.3 ggplot2/r3.4.3_2.2.1

java  -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar VariantRecalibrator \
        -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
	-V $PBS_O_WORKDIR/NexteraA.output.vcf.gz \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
	-resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode BOTH \
	--output $PBS_O_WORKDIR/NexteraA.output.all.recal \
	--tranches-file $PBS_O_WORKDIR/NexteraA.output.all.tranches \
	--rscript-file $PBS_O_WORKDIR/NexteraA.output.all.plots.R


## step 8: applyVQSR

#!/bin/bash
#PBS -l select=1:ncpus=4:mem=10gb
#PBS -l walltime=48:00:00
#PBS -N ApplyVQSR
#PBS -m abe
#PBS -M ....@gmail.com
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161 GenomeAnalysisTK/4.1.2.0

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar ApplyVQSR \
	-R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
	-V $PBS_O_WORKDIR/NexteraA.output.vcf.gz \
	-O $PBS_O_WORKDIR/NexteraA.output.vqsr.vcf.gz \
	-ts-filter-level 99.0 \
	--tranches-file $PBS_O_WORKDIR/NexteraA.output.all.tranches \
	--recal-file $PBS_O_WORKDIR/NexteraA.output.all.recal \
	-mode BOTH


java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar ApplyVQSR \
	-R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
	-V $PBS_O_WORKDIR/NexteraB.output.vcf.gz \
	-O $PBS_O_WORKDIR/NexteraB.output.vqsr.vcf.gz \
	-ts-filter-level 99.0 \
	--tranches-file $PBS_O_WORKDIR/NexteraB.output.all.tranches \
	--recal-file $PBS_O_WORKDIR/NexteraB.output.all.recal \
	-mode BOTH



## step 9: CalculateGenotypePosteriors
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
wget https:////storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx

module load java/8.0_161 GenomeAnalysisTK/4.1.2.0

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar CalculateGenotypePosteriors \
      -V NexteraA.output.vqsr.vcf.gz \
      -O NexteraA.output.1000G_PPs.vcf.gz \
      -supporting 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar CalculateGenotypePosteriors \
      -V NexteraB.output.vqsr.vcf.gz \
      -O NexteraB.output.1000G_PPs.vcf.gz \
      -supporting 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf


## step 10: VariantAnnotator to generate new .vcf file (note its different from the other selection)
cd $PBS_O_WORKDIR

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar VariantAnnotator \
     -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
     -I $PBS_O_WORKDIR/NexteraA.marked_duplicates.add.apply.bam \
     -V $PBS_O_WORKDIR/NexteraA.output.1000G_PPs.filter.vcf.gz \
     -O NexteraA.Refinement.vcf \
     -A Coverage \
     --dbsnp $PBS_O_WORKDIR/Homo_sapiens_assembly38.dbsnp138.vcf


java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar VariantAnnotator \
     -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
     -I $PBS_O_WORKDIR/NexteraB.marked_duplicates.add.apply.bam \
     -V $PBS_O_WORKDIR/NexteraB.output.1000G_PPs.filter.vcf.gz \
     -O NexteraB.Refinement.vcf \
     -A Coverage \
     --dbsnp $PBS_O_WORKDIR/Homo_sapiens_assembly38.dbsnp138.vcf


## step 11: SelectVariants
cd $PBS_O_WORKDIR

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants\
  -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
  -V $PBS_O_WORKDIR/NexteraA.Refinement.vcf \
  -select-type SNP \
  -O $PBS_O_WORKDIR/NexteraA.RawSNP.vcf


java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
  -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
  -V $PBS_O_WORKDIR/NexteraA.Refinement.vcf \
  -select-type INDEL \
  -O $PBS_O_WORKDIR/NexteraA.RawINDEL.vcf


java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants\
  -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
  -V $PBS_O_WORKDIR/NexteraB.Refinement.vcf \
  -select-type SNP \
  -O $PBS_O_WORKDIR/NexteraB.RawSNP.vcf


java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
  -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
  -V $PBS_O_WORKDIR/NexteraB.Refinement.vcf \
  -select-type INDEL \
  -O $PBS_O_WORKDIR/NexteraB.RawINDEL.vcf



# step 12: ANNOVAR
cd $PBS_O_WORKDIR

module load ANNOVAR/2018Apr16

table_annovar.pl $PBS_O_WORKDIR/NexteraA.RawSNP.vcf \
/home/BEP/2019/WGS_Exome/Reference/ANNOVAR/humandb_2018 \
-buildver hg38 \
-out $PBS_O_WORKDIR/NexteraA.RawSNP.annovar \
-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,exac03,exac03nontcga,exac03nonpsych,kaviar_20150923,avsnp150,dbnsfp35a,cosmic70,clinvar_20180603,nci60,hrcr1,mcap,revel,1000g2015aug_eur,intervar_20180118,regsnpintron,dbscsnv11 \
-operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . --thread 2 -vcfinput


table_annovar.pl $PBS_O_WORKDIR/NexteraA.RawINDEL.vcf \
/home/BEP/2019/WGS_Exome/Reference/ANNOVAR/humandb_2018 \
-buildver hg38 \
-out $PBS_O_WORKDIR/NexteraA.RawINDEL.annovar \
-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,exac03,exac03nontcga,exac03nonpsych,kaviar_20150923,avsnp150,dbnsfp35a,cosmic70,clinvar_20180603,nci60,hrcr1,mcap,revel,1000g2015aug_eur,regsnpintron,dbscsnv11 \
-operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . --thread 2 -vcfinput
