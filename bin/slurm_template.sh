#!/bin/bash -l
#SBATCH -A b2011097
#SBATCH -p node -n 1
#SBATCH -t 26:00:00
#SBATCH -J ccchr9q_w2

#export BASE=$HOME
export BASE=/home/danieln

module load bioinfo-tools
module load samtools/1.1
module load bcftools/1.1


samtools_bin=`which samtools`
export SAMTOOLS=$samtools_bin
export SAMTOOLSBIN=`dirname $samtools_bin`
export BCFTOOLS=bcftools
#export VCFUTILS=$SAMTOOLSBIN/vcfutils.pl

#export SAMTOOLS=$BASE/src/samtools-0.1.18/samtools
#export BCFTOOLS=$BASE/src/samtools-0.1.18/bcftools/bcftools
#export VCFUTILS=$BASE/src/samtools-0.1.18/bcftools/vcfutils.pl

# FastQC wants java
module load java/sun_jdk1.6.0_18

# use local Mosaik build (it's the latest version anyway)
module load bioinfo-tools

module load mosaik-aligner/2.2.30
mosaik_aligner=`which MosaikAligner`
export MOSAIKBIN=`dirname $mosaik_aligner`
export REFERENCE=/home/danieln/glob/private/etiologica_ref/hs37d5.fa.gz

export mismatches=14
export sw_bandwidth=33
export clustersize=35
export JUMP=yes
export mjump=15
export mhp=100

export MOSAIK_TMP=$SNIC_TMP
export TMP=$SNIC_TMP
export TEMP=$SNIC_TMP

# export ANN_PATH=../src/networkFile
# -annpe $ANN_PATH/2.1.26.pe.100.0065.ann -annse $ANN_PATH/2.1.26.se.100.005.ann
export MOSAIK_DAT_ON_SCRATCH="yes"
#MOSIK_DAT_ON_SCRATCH=yes

export BINDIR=$BASE/sandbox/etiologica/bin

export ANNOVARBIN=$BASE/src/annovar
export ANNOVAR_DISPENSABLE=$BASE/src/annovar/example/dispensable.all
export ANNOVAR_1KG_MAF=0.02
export ANNOVAR_PP2_BENIGN=0.85
export AVDBDIR=$BASE/src/annovar/humandb
export LOCAL_CLIN_DB="hg19_100clinical.real.avdb"
export LOCAL_DANES_DB="hg19_200danes.avdb"
export DB_SNP_VERSION="snp135NonFlagged"



# UPPMAX avdbdir is not up to date
#export AVDBDIR=/bubo/nobackup/uppnex/annotations/annovar/humandb/

export FASTQC=$BASE/src/FastQC/fastqc

#cd /proj/b2011097/private/chr9q

# first pass..
#for file in *fastq ; 
#do 
# pigz $file
#done

$BINDIR/run_secondary_analysis.sh
$BINDIR/run_tertiary_analysis.sh
