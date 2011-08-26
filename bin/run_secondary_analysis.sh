#!/bin/bash

# (c)2011 Daniel Nilsson, daniel.nilsson@scilifelab.se, daniel.nilsson@karolinska.se, daniel.nilsson@ki.se
# Copyright 2011, held by the author.
# Released under the Perl Artistic License.
#
# Documentation available in this file. Use e.g. perldoc to view, or pod2man to create a manual.

# . etiologica_site_config.sh

### BEGIN USER SPECIFIC CONFIG

BINDIR=/home/daniel/sandbox/etiologica/bin

PIPELINEFUNK=$BINDIR/pipelinefunk.sh
TMP=$SCRATCH

REFERENCE=human_g1k_v37.fasta.gz

MOSAIKBIN=/home/daniel/src/mosaik-aligner/bin
SAMTOOLS=/home/daniel/src/samtools-0.1.12a/samtools
BCFTOOLS=../src/samtools-0.1.12a/bcftools/bcftools
VCFUTILS=../src/samtools-0.1.12a/bcftools/vcfutils.pl

### END CONFIG

### BEGIN MAIN DOCUMENTATION 

: << 'POD_INIT'

=head1 NAME

run_secondary_analysis.sh

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@scilifelab.se, daniel.nilsson@ki.se, daniel.nilsson@karolinska.se

=head1 SYNOPSIS

 mkdir myproject
 ln -s patientid.fastq myproject/
 run_secondary_analysis.sh

=head1 DESCRIPTION

Process a patient fastq files to generate a list of variants.

=head1 COPYRIGHT AND LICENSE

Copyright Daniel Nilsson, 2011. 

The package is released under the Perl Artistic License.

=head1 BUGS

Surely!

=cut

POD_INIT

### END MAIN DOCUMENTATION

. $PIPELINEFUNK

# align

export MOSAIK_TMP = $TMP

referencedat=${REFERENCE%%.fasta.gz}.dat
mjump=15
mhp=100
reference_jump=${reference_dat%%.dat}.$mjump

if needsUpdate $referencedat $REFERENCE $MOSAIKBIN
then
    $MOSAIKBIN/MosaikBuild -fr $REFERENCE -oa $reference_dat
    $MOSAIKBIN/MosaikJump -ia $reference_dat -hs $mjump -out $reference_jump -mhp $mhp
fi

for dir in patient group metadata
do 
    if [ ! -d "$dir" ]
    then
	mkdir $dir 
    fi
done

if [  "" != "`ls -1 | grep fastq.gz`"  ]
then
    for patientfastq in *fastq.gz
    if [ "$patientfastq" != "$REFERENCE" ]	
    then 
	mv $patientfastq patient/
    fi
fi

for patient_fastq_gz in patient/*fastq.gz
do 
    patient_dat=${patient_fastq_gz%%fastq.gz}dat
    if needsUpdate $patient_dat $patient_fastq_gz 
    then
	$MOSAIKBIN/MosaikBuild -q $patient_fastq_gz -out $patient_dat -st illumina	
	registerFile $patient_dat temp
    fi

    patient_dat=${patient_dat%%dat}mosaik.dat
    if needsUpdate $patient_aln_dat $patient_dat $reference_jump $MOSAIKBIN 
    then
	$MOSAIKBIN/MosaikAligner -in $patient_dat -ia $reference_dat -out $patient_aln_dat -m unique -hs $mjump -bw 29 -j $reference_jump -mhp $mhp -mm 12 -act 35 -p 10
	registerFile $patient_aln_dat temp
    fi

    patient_sorted=${patient_aln_dat%%dat}sorted.dat
    if needsUpdate ${patient_sorted} $patient_aln_dat $MOSAIKBIN
    then
	$MOSAIKBIN/MosaikSort -in $patient_aln_dat -out $patient_sorted
	registerFile $patient_sorted temp
    fi

# ../src/mosaik-aligner/bin/MosaikDupSnoop 
# ../src/mosaik-aligner/bin/MosaikSort -dup 

    patient_bam=${patient_sorted%%sorted.dat}bam
    if needsUpdate ${patient_bam} ${patient_sorted}
    then
	$MOSAIKBIN/MosaikText -in $patient_sorted $patient_sorted -bam $patient_bam
	registerFile $patient_bam result
    fi

    patient_bcf=${patient_fastq_gz%%fastq.gz}var.raw.bcf
    $SAMTOOLS mpileup -ugf $REFERENCE $patient_bam | $BCFTOOLS view -bvcg - > $patient_bcf

    patient_vcf=${patient%%raw.bcf}flt.vcf
    if needsUpdate $patient_vcf $patient_bcf $BCFTOOLS $VCFUTILS
    then
	$BCFTOOLS view $patient_bcf | $VCFUTILS varFilter -D200 > $patient_vcf
	registerFile $patient_vcf result
    fi

done