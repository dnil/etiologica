#!/bin/bash -l
#SBATCH -A b2011097
#SBATCH -p node -n 1
#SBATCH -t 26:00:00
#SBATCH -J ccchr9q_w2

# FastQC wants java
module load java/sun_jdk1.6.0_18

# use local Mosaik build (it's the latest version anyway)
module load bioinfo-tools
module load mosaik-aligner/1.1.0021
mosaik_aligner=`which MosaikAligner`
export MOSAIKBIN=`dirname $mosaik_aligner`

export REFERENCE=/proj/b2011146/private/reference/human_g1k_v37.fasta.gz

export mismatches=14
export sw_bandwidth=33
export clustersize=35
export JUMP=yes
export mjump=15
export mhp=100

export MOSAIK_TMP=$SNIC_TMP
export TMP=$SNIC_TMP
export TEMP=$SNIC_TMP

export MOSAIK_DAT_ON_SCRATCH="yes"
#MOSIK_DAT_ON_SCRATCH=yes

export BINDIR=$HOME/sandbox/etiologica/bin

export ANNOVARBIN=$HOME/src/annovar
export ANNOVAR_DISPENSABLE=$HOME/src/annovar/example/dispensable.all
export ANNOVAR_1KG_MAF=0.02
export ANNOVAR_PP2_BENIGN=0.85
export AVDBDIR=$HOME/src/annovar/humandb
export LOCAL_CLIN_DB="hg19_100clinical.real.avdb"
export LOCAL_DANES_DB="hg19_200danes.avdb"
export DB_SNP_VERSION="snp135NonFlagged"

# UPPMAX avdbdir is not up to date
#export AVDBDIR=/bubo/nobackup/uppnex/annotations/annovar/humandb/

export SAMTOOLS=$HOME/src/samtools-0.1.18/samtools
export BCFTOOLS=$HOME/src/samtools-0.1.18/bcftools/bcftools
export VCFUTILS=$HOME/src/samtools-0.1.18/bcftools/vcfutils.pl
export FASTQC=$HOME/src/FastQC/fastqc

#cd /proj/b2011097/private/chr9q

# first pass..
#for file in *fastq ; 
#do 
#    cat $file | gzip -c > ${file}.gz
#done

$BINDIR/run_secondary_analysis.sh
