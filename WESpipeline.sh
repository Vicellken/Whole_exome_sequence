###FastQC
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=40gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N FastQC
#PBS -q bep

cd $PBS_O_WORK/DIR
mkdir FastQC

#Load the module required
module load java/8.0_161
module load FastQC/0.11.2

#Run FastQC in a for loop for multiple Fastq file
for file in *.fastq.bz2
do
    echo "This script is about to run another script."
    export file
    sh ./fastqc.sh
    echo "This script has just run another script."
done

###fastqc.sh
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=40gb
#PBS -l walltime=1:00:00
#PBS -m abe
#PBS -N fastqc
#PBS -q bep

cd $PBS_O_WORKDIR


fastqc --outdir FastQC --noextract $file

###Trimmomatic
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=40gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N trimmomatic
#PBS -q bep

module load java/8.0_161
module load Trimmomatic

cd $PBS_O_WORKDIR
#use trimmomatic to do quality filtering
java -Xms4g -Xmx4g -jar /software/Trimmomatic/0.38/trimmomatic-0.38.jar PE -threads 2 -phred33 \
NexteraA_1.fastq.bz2 NexteraA_2.fastq.bz2 NexteraA_1_trimmed.fastq NexteraA_1_unpaired.fastq \
NexteraA_2_trimmed.fastq NexteraA_2_unpaired.fastq \
ILLUMINACLIP:/software/Trimmomatic/0.38/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:5:11 MINLEN:50 \
&> trimmomatic.log

###Alignment
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=40gb
#PBS -l walltime=10:00:00
#PBS -m abe
#PBS -N alignment
#PBS -q bep

cd $PBS_O_WORKDIR

module load bwa/0.7.17
#echo [MSG] Building Index (FYI)
#bwa index reference/GRCh38.primary_assembly.genome.fa
bwa mem -t 2 reference/GRCh38.primary_assembly.genome.fa NexteraA_1_trimmed.fastq NexteraA_2_trimmed.fastq > NexteraA.mapped.sam

module load samtools/1.8
samtools view -S -b NexteraA.mapped.sam > NexteraA.mapped.bam
samtools sort -o NexteraA.mapped.sorted.bam NexteraA.mapped.bam
samtools index NexteraA.mapped.sorted.bam
samtools rmdup NexteraA.mapped.sorted.bam NexteraA.rmdup.mapped.sorted.bam

###Mark duplicates
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l ncpus=2
#PBS -l mem=40gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N picard
#PBS -q bep

#load modules for picard
module load java/8.0_161
module load Picard/2.18.9

cd $PBS_O_WORKDIR

#mark duplicates (picard and samtools rmdup can do the same thing)
java -jar /software/Picard/2.18.9/picard.jar MarkDuplicates \
    I=NexteraA.mapped.sorted.bam \
    REMOVE_DUPLICATES=true \
    O=NexteraA.marked_duplicates.bam \
    M=NexteraA.marked_dup_metrics.txt
#gatk-package-4.beta.x-spark.jar is the jar for running Spark tools on a Spark cluster,
#while gatk-package-4.beta.x-local.jar is the jar that is used for everything else (including running Spark tools "locally", ie on a regular server or cluster).

###create index file for reference
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l ncpus=2
#PBS -l mem=40gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N prepare_reference
#PBS -q bep

cd $PBS_O_WORKDIR

module load samtools/1.8

#creat fai.file index
samtools faidx GRCh38.primary_assembly.genome.fa

#load modules for picard
module load java/8.0_161
module load Picard/2.18.9

#create dict. file
java -jar /software/Picard/2.18.9/picard.jar CreateSequenceDictionary \
      R=GRCh38.primary_assembly.genome.fa \
      O=GRCh38.primary_assembly.genome.dict

###Base Recalibration
#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=40gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N baserecalibrator
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161
module load Picard/2.18.9

java -jar /software/Picard/2.18.9/picard.jar AddOrReplaceReadGroups \
       I=NexteraA.marked_duplicates.bam \
       O=NexteraA.marked_duplicates_R.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=20

module load GenomeAnalysisTK/4.1.2.0
cd $PBS_O_WORKDIR

#GATK base recalibrator after marking duplicates
java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar BaseRecalibrator \
      -I NexteraA.marked_duplicates_R.bam \
      --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
      -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
      -O NexteraA.recal_data.table

###applyBQSR
#Specifically, it recalibrates the base qualities of the input reads based on the recalibration table produced by the BaseRecalibrator tool, and outputs a recalibrated BAM or CRAM file.
#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=40gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N applyBQSR
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0

#GATK apply BQSR(base quality score recalibration) after base recalibration
#the PrintReads tool from GATK3.8 has become ApplyBQSR in GATK4
java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar ApplyBQSR \
     -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
     -I NexteraA.marked_duplicates_R.bam \
     --bqsr-recal-file NexteraA.recal_data.table \
     -O NexteraA.marked_duplicates_R_BQSR.bam

###HaplotypeCaller
#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=40gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N haplotypecaller
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar HaplotypeCaller \
     -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
     -I NexteraA.marked_duplicates_R_BQSR.bam \
     -O NexteraA.raw.snps.indels.vcf
#Primary assembly refers to the collection of (i) assembled chromosomes, (ii) unlocalized and (iii) unplaced sequences. It represents a non-redundant haploid genome.
#(i) Assembled chromosomes for hg38 are chromosomes 1–22 (chr1–chr22), X (chrX), Y (chrY) and Mitochondrial (chrM).
#(ii) Unlocalized sequence are on a specific chromosome but with unknown order or orientation. Identify by _random suffix.
#(iii) Unplaced sequence are on an unknown chromosome. Identify by chrU_ prefix.
#Primary assembly contains all toplevel sequence regions excluding haplotypes and patches. This file is best used for performing sequence similarity searches where patch and haplotype sequences would confuse analysis.


###VariantRecalibrator
#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=40gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N VariantRecalibrator
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0
module load R/3.4.3
module load ggplot2/r3.4.3_2.2.1

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar VariantRecalibrator \
     -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
     -V NexteraA.raw.snps.indels.vcf \
     --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
     --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz \
     --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
     --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf \
     -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
     -mode BOTH \
     #Use either SNP for recalibrating only SNPs (emitting indels untouched in the output VCF) or INDEL for indels (emitting SNPs untouched in the output VCF). \
     #There is also a BOTH option for recalibrating both SNPs and indels simultaneously, but this is meant for testing purposes only and should not be used in actual analyses.
     -O NexteraA.recal \
     --tranches-file NexteraA.tranches \
     --rscript-file NexteraA.plots.R

###ApplyVQSR
#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=40gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N ApplyVQSR
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar ApplyVQSR \
     -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
     -V NexteraA.raw.snps.indels.vcf \
     -O NexteraA.ApplyVQSR.vcf.gz \
     --truth-sensitivity-filter-level 99.0 \
     --tranches-file NexteraA.tranches \
     --recal-file NexteraA.recal \
     -mode BOTH
     #The --mode argument is an enumerated type (Mode), which can have one of the following values

###CalculateGenotypePosteriors
#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=40gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N CalculateGenotypePosteriors
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0

#Genotype Refinement workflow
#Step 1: Derive posterior probabilities of genotypes
java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar CalculateGenotypePosteriors \
      -V NexteraA.ApplyVQSR.vcf.gz \
      -O NexteraA.1000G_PPs.vcf.gz \
      -supporting 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf

###VariantFiltration
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=20gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N VariantFiltration
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0

#Genotype Refinement workflow
#Step 2: Filter low quality genotypes
#Hard Filter
java -Xmx30g -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar VariantFiltration \
     -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
     --variant $PBS_O_WORKDIR/NexteraA.1000G_PPs.vcf.gz \
     -O $PBS_O_WORKDIR/NexteraA.1000G_PPs.filtered.vcf.gz \
     -G_filter "GQ < 20.0" -G_filterName "lowGQ"

###VariantAnnotator
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=20gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N VariantAnnotator
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0

#Genotype Refinement workflow
#Step 3: Annotate possible de novo mutations
java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar VariantAnnotator \
     -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
     -I $PBS_O_WORKDIR/NexteraA.marked_duplicates_R_BQSR.bam \
     -V $PBS_O_WORKDIR/NexteraA.1000G_PPs.filtered.vcf.gz \
     -O NexteraA.Refinement.vcf \
     -A Coverage \
     --dbsnp $PBS_O_WORKDIR/Homo_sapiens_assembly38.dbsnp138.vcf

###SelectVariants
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=20gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N SelectVariants
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0

# SNP
java -Xmx30g -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
     -R /home/BEP/2019/WGS_Exome/Reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
     --variant $PBS_O_WORKDIR/NexteraA.Refinement.vcf \
     -O $PBS_O_WORKDIR/NexteraA.RawSNP.vcf \
     --select-type-to-include SNP

#INDEL
java -Xmx30g -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
     -R /home/BEP/2019/WGS_Exome/Reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
     --variant $PBS_O_WORKDIR/NexteraA.Refinement.vcf \
     -O $PBS_O_WORKDIR/NexteraA.RawINDEL.vcf \
     --select-type-to-include INDEL

###ANNOVARSNP
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=20gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N ANNOVARSNP
#PBS -q bep

cd $PBS_O_WORKDIR

module load ANNOVAR/2018Apr16

table_annovar.pl $PBS_O_WORKDIR/NexteraA.RawSNP.vcf \
/home/BEP/2019/WGS_Exome/Reference/ANNOVAR/humandb_2018 \
-buildver hg38 \
-out $PBS_O_WORKDIR/NexteraA.RawSNP.annovar \
-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,exac03,exac03nontcga,exac03nonpsych,kaviar_20150923,avsnp150,dbnsfp35a,cosmic70,clinvar_20180603,nci60,hrcr1,mcap,revel,1000g2015aug_eur,intervar_20180118,regsnpintron,dbscsnv11 \
-operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . --thread 2 -vcfinput

###ANNOVARINDEL
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=20gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N ANNOVARSNP
#PBS -q bep

table_annovar.pl $PBS_O_WORKDIR/NexteraA.RawINDEL.vcf \
    /home/BEP/2019/WGS_Exome/Reference/ANNOVAR/humandb_2018 \
    -buildver hg38 \
    -out $PBS_O_WORKDIR/NexteraA.RawINDEL.annovar \
-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,exac03,exac03nontcga,exac03nonpsych,kaviar_20150923,avsnp150,dbnsfp35a,cosmic70,clinvar_20180603,nci60,hrcr1,mcap,revel,1000g2015aug_eur,regsnpintron,dbscsnv11 \
-operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . --thread 2 -vcfinput

###VariantEval
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=20gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N VariantEval
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0

#creat index for vcf file
java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar IndexFeatureFile \
     -F NexteraA.RawSNP.annovar.hg38_multianno.vcf

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar VariantEval \
     -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
     -O NexteraA.eval.grp \
     --eval NexteraA.RawSNP.annovar.hg38_multianno.vcf

###CNNScoreVariants
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=20gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N CNNScoreVariants
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true \
-Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 \
-jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar CNNScoreVariants \
-R /home/xinfeng/Exome/Project_1and2/primary_seq/reference/GRCh38.primary_assembly.genome.fa \
--disable-avx-check -V NexteraA.raw.snps.indels.vcf \
-O NexteraA.annotated.vcf


rm -r /home/xinfeng/.conda/envs/gatk

###FilterVariantTranches
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=20gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N FilterVariantTranches
#PBS -q bep

module load java/8.0_161
module load GenomeAnalysisTK/4.1.2.0

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar FilterVariantTranches \
   -V NexteraA.annotated.vcf.gz \
   --resource hapmap_3.3.hg38.vcf.gz \
   --resource Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   --info-key CNN_1D \
   --tranche 99.9 --tranche 99.0 --tranche 95 \
   #--tranche -t
   #The levels of truth sensitivity at which to slice the data. (in percents, i.e. 99.9 for 99.9 percent and 1.0 for 1 percent)
   -O filtered.vcf

###DepthOfCoverage
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=10gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -N DepthOfCoverage
#PBS -q bep

cd $PBS_O_WORKDIR

module load java/8.0_161 GenomeAnalysisTK/3.8.1.0

java -jar /software/GenomeAnalysisTK/3.8.1.0/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
-o DepthOfCoverage \
-I NexteraA.marked_duplicates_R_BQSR.bam \
-L hglft_genome_3f517_e53cd0.bed \
--omitDepthOutputAtEachBase


###Qualimap
#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=20gb
#PBS -l walltime=120:00:00
#PBS -m abe
#PBS -q bep
#PBS -N AlignStat_qualimap

cd $PBS_O_WORKDIR

module load R/3.4.3
module load qualimap/2.2.2_26-08-18

qualimap bamqc \
-bam $PBS_O_WORKDIR/NexteraA.marked_duplicates_R_BQSR.bam \
-nt 12 \
-outdir $PBS_O_WORKDIR/Qualimap \
--feature-file $PBS_O_WORKDIR/hglft_genome_3f517_e53cd0.bed \
--collect-overlap-pairs \
--skip-duplicated \
--paint-chromosome-limits
