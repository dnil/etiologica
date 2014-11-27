#!/bin/bash -l
#SBATCH -A cust002
#SBATCH -n 6
#SBATCH -t 24:00:00
#SBATCH -J trim_clin_exomes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.nilsson@scilifelab.se

export BASE=$HOME
export BINDIR=$BASE/sandbox/etiologica/bin

PIPELINE=etiologica
PIPELINEFUNK=$BINDIR/pipelinefunk.sh

. $PIPELINEFUNK
#log "Running $PIPELINE." "main"
#log "BINDIR: $BINDIR" "main"

for pat in `ls *fastq.gz |cut -f4 -d_ |sort |uniq` ; 
do 
    echo $pat; 
    mkdir $pat ; 
    mv *$pat*fastq.gz $pat/ ;
done


for pat in `ls */*fastq.gz |cut -f4 -d_ |sort |uniq` ; 
do
    cd $pat
    for fwd in *_1.fastq.gz; 
    do
	if workLockOk $fwd
	then     
	    rev=${fwd%%_1.fastq.gz}_2.fastq.gz
	    
	    fwd_trimmed=${fwd%%.fastq.gz}.trimmed.fastq.gz
	    rev_trimmed=${rev%%.fastq.gz}.trimmed.fastq.gz
  
	    fwd_unpaired=${fwd%%.fastq.gz}.unpaired.trimmed.fastq.gz
	    rev_unpaired=${rev%%.fastq.gz}.unpaired.trimmed.fastq.gz

	    java -jar ~/src/Trimmomatic-0-3.32/trimmomatic-0.32.jar PE $fwd $rev $fwd_trimmed $fwd_unpaired $rev_trimmed $rev_unpaired ILLUMINACLIP:$HOME/src/Trimmomatic-0-3.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	fi
    done
    cd ..
done
