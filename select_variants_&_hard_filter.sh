#!/bin/bash
#PBS -l select=1:ncpus=4:mem=10gb
#PBS -l walltime=48:00:00
#PBS -N ApplyVQSR
#PBS -m abe
#PBS -M ....@gmail.com
#PBS -q bep

cd $PBS_O_WORKDIR

# generate select variants

module load java/8.0_161 GenomeAnalysisTK/4.1.2.0

## notice the difference in syntax for different GATK version
java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants\
  -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
  -V $PBS_O_WORKDIR/NexteraA.output.vcf.gz \
  -select-type SNP \
  -O $PBS_O_WORKDIR/RawSNP_A.vcf


java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
  -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
  -V $PBS_O_WORKDIR/NexteraA.output.vcf.gz \
  -select-type INDEL \
  -O $PBS_O_WORKDIR/RawINDEL_A.vcf


  java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
    -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
    -V $PBS_O_WORKDIR/NexteraB.output.vcf.gz \
    -select-type SNP \
    -O $PBS_O_WORKDIR/RawSNP_B.vcf


  java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
    -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
    -V $PBS_O_WORKDIR/NexteraB.output.vcf.gz \
    -select-type INDEL \
    -O $PBS_O_WORKDIR/RawINDEL_B.vcf


### VariantFiltration in GATK4 (this is hard filtering step)
module load java/8.0_161 GenomeAnalysisTK/4.1.2.0

java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar VariantFiltration \
  -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
  -V $PBS_O_WORKDIR/RawSNP_A.vcf \
  -O $PBS_O_WORKDIR/RawSNP_Filtered_A.vcf \
  --filterExpression "QD < 2.0" --filterName "LowQD" \
  --filterExpression "FS > 60.0" --filterName "StrandBias" \
  --filterExpression "MQ < 40.0" --filterName "LowMQ" \
  --filterExpression "MQRankSum < -12.5" --filterName "LowMQRankSum" \
  --filterExpression "ReadPosRankSum < -8.0" --filterName "LowReadPosRankSum"


  java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar VariantFiltration \
    -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
    -V $PBS_O_WORKDIR/RawINDEL_A.vcf \
    -O $PBS_O_WORKDIR/RawINDEL_Filtered_A.vcf \
    --filterExpression "QD < 2.0" --filterName "LowQD" \
    --filterExpression "FS > 200.0" --filterName "StrandBias" \
    --filterExpression "ReadPosRankSum < -20.0" --filterName "LowReadPosRankSum"


    java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar VariantFiltration \
      -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
      -V $PBS_O_WORKDIR/RawSNP_B.vcf \
      -O $PBS_O_WORKDIR/RawSNP_Filtered_B.vcf \
      --filterExpression "QD < 2.0" --filterName "LowQD" \
      --filterExpression "FS > 60.0" --filterName "StrandBias" \
      --filterExpression "MQ < 40.0" --filterName "LowMQ" \
      --filterExpression "MQRankSum < -12.5" --filterName "LowMQRankSum" \
      --filterExpression "ReadPosRankSum < -8.0" --filterName "LowReadPosRankSum"


      java -jar /software/GenomeAnalysisTK/4.1.2.0/gatk-package-4.1.2.0-local.jar VariantFiltration \
        -R $PBS_O_WORKDIR/reference/GRCh38.primary_assembly.genome.fa \
        -V $PBS_O_WORKDIR/RawINDEL_B.vcf \
        -O $PBS_O_WORKDIR/RawINDEL_Filtered_B.vcf \
        --filterExpression "QD < 2.0" --filterName "LowQD" \
        --filterExpression "FS > 200.0" --filterName "StrandBias" \
        --filterExpression "ReadPosRankSum < -20.0" --filterName "LowReadPosRankSum"
