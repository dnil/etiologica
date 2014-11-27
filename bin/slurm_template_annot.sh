#!/bin/bash -l
#SBATCH -A cust002
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J annotate_clin_exomes

#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.nilsson@scilifelab.se

export BASE=$HOME
export BINDIR=$BASE/sandbox/etiologica/bin

export ANNOVARBIN=$BASE/src/annovar

/home/daniel.nilsson/sandbox/etiologica/bin/convert_and_annotate.sh
