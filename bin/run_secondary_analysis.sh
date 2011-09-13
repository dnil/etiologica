#!/bin/bash
# (c)2011 Daniel Nilsson, daniel.nilsson@scilifelab.se, daniel.nilsson@karolinska.se, daniel.nilsson@ki.se
# Copyright 2011, held by the author.
# Released under the Perl Artistic License.
#
# Documentation available in this file. Use e.g. perldoc to view, or pod2man to create a manual.

# . etiologica_site_config.sh
#
# TODO: decide on parameter config structures 
#

### BEGIN USER SPECIFIC CONFIG


BINDIR=/home/daniel/sandbox/etiologica/bin

PIPELINE=etiologica
PIPELINEFUNK=$BINDIR/pipelinefunk.sh

REFERENCE=human_g1k_v37.fasta.gz

MOSAIKBIN=/home/daniel/src/mosaik-aligner-read-only/bin
ANNOVARBIN=/home/daniel/src/annovar
SAMTOOLS=/home/daniel/src/samtools-0.1.17/samtools
BCFTOOLS=/home/daniel/src/samtools-0.1.17/bcftools/bcftools
VCFUTILS=/home/daniel/src/samtools-0.1.17/bcftools/vcfutils.pl
FASTQC=/home/daniel/src/FastQC/fastqc
AVDBDIR=/home/daniel/src/annovar/humandb/

# MOSAIK settings
mismatches=14
sw_bandwidth=33
clustersize=35
JUMP=yes
mjump=15
mhp=100

# VCFTOOLS
var_filter_settings="-D1000 -d6"

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
 cd myproject
 run_secondary_analysis.sh

=head1 DESCRIPTION

Process patient fastq.gz files to generate a list of variants.
Uses FastQC, MOSAIK, SAMTOOLS and ANNOVAR.

=head1 MULTIPLE WORKER INSTANCES

Multiple worker safe using fs lock / bash noclobber. Feel free to start multiple work instances 
in the same directory, and they should coexist on a first come-first served basis. One worker per
node may be appropriate - also see the MOSAIC_CORES variable. Work blocks are large, roughly a complete pipleine for lane in each exclusive block.

=cut

POD_INIT

### END INITIAL DOCUMENTATION -- more to follow inline where appropriate.

. $PIPELINEFUNK

: <<POD_ENV

=head1 ENVIRONMENT VARIABLES

The pipeline adheres to environment exported variables.

=over 4

=item C<SCRATCH> I<path>

Designated temp area. Please make it local and sufficiently big for the sake of MOSAIK.
Can be overruled by directly setting TMP or MOSAIK_TMP, with the latter having the highest precedence rank.

=item C<MOSAIK_CORES> I<int>

Number of cores to use for MOSAIK. 
Defaults to the number of available cores (see NPROC in the pipelinefunc library).

=back

=cut

POD_ENV

if [ -z "$SCRATCH" ]
then
    TMP=$SCRATCH
fi 

if [ -z "$TMP" ] 
then 
    export TMP=/tmp
fi

if [ -z "$MOSAIK_TMP" ]
then 
    export MOSAIK_TMP=$TMP
fi 

if [ -z "$MOSAIK_CORES" ]
then
    if [ -z "$NPROC" ]
    then 
	MOSAIK_CORES=10
    else 
	MOSAIK_CORES=$NPROC
    fi
fi

if [ -z "$JUMP" ]
then
    JUMP="yes"
fi

for dir in patient group metadata
do 
    if [ ! -d "$dir" ]
    then
	mkdir $dir 
    fi
done

reference_dat=${REFERENCE%%.fasta.gz}.dat

if needsUpdate $reference_dat $REFERENCE $MOSAIKBIN
then
    runme="$MOSAIKBIN/MosaikBuild -fr $REFERENCE -oa $reference_dat"
    vanillaRun "$runme" "$reference_dat" "temp" "MosaikBuild"
fi

if [ "$JUMP" == "yes" ]
then
    reference_jump=${reference_dat%%.dat}.$mjump
    if needsUpdate $reference_jump $reference_dat $MOSAIKBIN/MosaikJump
    then
	runme="$MOSAIKBIN/MosaikJump -ia $reference_dat -hs $mjump -out $reference_jump -mhp $mhp"
	vanillaRun "$runme" "$reference_jump" "temp" "MosaikJump"
    fi
fi

if [ "" != "`ls -1 | grep fastq.gz`" ]
then
    for patientfastq in *fastq.gz
    do
	if [ "$patientfastq" != "$REFERENCE" ]	
	then 
	    # ok with links to fastq.gz?
	    mv $patientfastq patient/
	fi
    done
fi

: <<'POD_FASTQC'

=head1 QUALITY CHECK

Runs FastQC for each fastq.gz file.

POD_FASTQC

for patientfastq in patient/*fastq.gz
do    
    patientfastqc_zip=${patientfastq%%.fastq.gz}_fastqc.zip

    if needsUpdate $patientfastqc_zip $patientfastq $FASTQC 
    then
	# ensure that the zip archive exists before trying to lock it (first run).
	touch $patientfastqc_zip

	if workLockOk $patientfastqc_zip
	then
	    runme="$FASTQC $patientfastq"
	    vanillaRun "$runme" "$patientfastqc_zip" "result" "FastQC"

	    releaseLock $patientfastqc_zip
	fi
    fi
done

: << 'POD_MATE'

=head1 MATEPAIRS

The pipeline assumes mated reads are found as _1.fastq.gz and _2.fastq.gz archives.
If a mate lib is missing for any, it is assumed all reads in the project are unpaired.

=cut

POD_MATE

# single reads or mates?
if [ -z "$MATEPAIRS" ] 
then 
    
    MATEPAIRS=1
    for patient_fastq_gz in patient/*_1.fastq.gz 
    do
	possible_mate=${patient_fastq_gz%%_1.fastq.gz}_2.fastq.gz
	if [ ! -s $possible_mate ] 
	then
	    # possible mate according to convention not found.
	    echo "Did not find any ($possible) mate for $patient_fasta_gz. Assuming all reads are unpaired."

	    MATEPAIRS=0
	fi
    done

    for patient_fastq_gz in patient/*_2.fastq.gz 
    do
	possible_mate=${patient_fastq_gz%%_2.fastq.gz}_1.fastq.gz
	if [ ! -s $possible_mate ] 
	then
	    # possible mate according to convention not found.
	    echo "Did not find any ($possible) mate for $patient_fasta_gz. Assuming all reads are unpaired."

	    MATEPAIRS=0
	fi
    done

fi

if [ "$MATEPAIRS" == "0" ]
then
    patient_fastq_list=( patient/*fastq.gz )
else
    patient_fastq_list=( patient/*_1.fastq.gz )
fi

for patient_fastq_gz in ${patient_fastq_list[@]}
do 
    # 
    # File locking for multiple script instances concurrently in same directory.
    # One instance per patient fastq.gz is a rather appropriate work block right now - if 
    # multilane samples start popping in this recommendation may change.
    #

    if workLockOk $patient_fastq_gz
    then
	patient_dat=${patient_fastq_gz%%fastq.gz}dat

	if [ "$MATEPAIRS" == "0" ]
	then
	    if needsUpdate $patient_dat $patient_fastq_gz
	    then
		runme="$MOSAIKBIN/MosaikBuild -q $patient_fastq_gz -out $patient_dat -st illumina"
		vanillaRun "$runme" "$patient_dat" "temp" "MosaikBuild"
	    fi
	else
	    patient_2_fastq_gz=${patient_fastq_gz%%_1.fastq.gz}_2.fastq.gz
	    if needsUpdate $patient_dat $patient_fastq_gz $patient_2_fastq_gz
	    then
		runme="$MOSAIKBIN/MosaikBuild -q $patient_fastq_gz -q2 $patient_2_fastq_gz -out $patient_dat -st illumina"
		vanillaRun "$runme" "$patient_dat" "temp" "MosaikBuild"
	    fi
	fi

	MOSAIK_ALIGN_MODE=all
	if [ "$MATEPAIRS" == "0" ]
	then
	    MOSAIK_ALIGN_MODE=unique
	    dupstring=""
	fi
       
	patient_aln_dat=${patient_dat%%dat}mosaik.dat
	if needsUpdate $patient_aln_dat $patient_dat $reference_jump $MOSAIKBIN/MosaikAligner
	then
	    if [ "$JUMP" == "yes" ]
	    then
		runme="$MOSAIKBIN/MosaikAligner -in $patient_dat -ia $reference_dat -out $patient_aln_dat -m $MOSAIK_ALIGN_MODE -hs $mjump -bw $sw_bandwidth -j $reference_jump -mhp $mhp -mm $mismatches -act $clustersize -p $MOSAIK_CORES"
	    else
		runme="$MOSAIKBIN/MosaikAligner -in $patient_dat -ia $reference_dat -out $patient_aln_dat -m $MOSAIK_ALIGN_MODE -bw $sw_bandwidth -mm $mismatches -act $clustersize -p $MOSAIK_CORES"
	    fi
	    vanillaRun "$runme" "$patient_aln_dat" "temp" "MosaikAligner"
	fi

	: <<POD_MOSAIKDUP

=head1 DEAL WITH DUPLICATES

MosaikDupSnoop is enabled to deal with duplicate reads in paired mode.
For unpaired sequences, MOSAIK-Aligner is run in C<unique> mode and 
duplicate checking is turned off.

=cut

POD_MOSAIKDUP

	patient_lib_dupdata_dir=${patient_aln_dat%%.mosaik.dat}_DupData
	if [ "$MATEPAIRS" != 0 ] 
	then
	    if needsUpdate ${patient_lib_dupdata_dir}/.db $patient_aln_dat $MOSAIKBIN/MosaikDupSnoop
	    then
		if [ ! -d $patient_lib_dupdata_dir ]
		then
		    mkdir $patient_lib_dupdata_dir
		fi
		
		runme="$MOSAIKBIN/MosaikDupSnoop -in $patient_aln_dat -od $patient_lib_dupdata_dir"
		vanillaRun "$runme" "$patient_lib_dupdata_dir/.db" "temp" "MosaikDupSnoop"
	    fi
	fi

	patient_sorted=${patient_aln_dat%%dat}sorted.dat
	if needsUpdate ${patient_sorted} $patient_aln_dat $MOSAIKBIN/MosaikSort
	then
	    if [ "$MATEPAIRS" == 0 ]
	    then 
		dupstring=""
	    else
		dupstring="-dup $patient_lib_dupdata_dir"
	    fi
	    
	    runme="$MOSAIKBIN/MosaikSort -in $patient_aln_dat -out $patient_sorted $dupstring"
	    vanillaRun "$runme" "$patient_sorted" "temp" "MosaikSort"
	fi

       # TODO: mosaik merge, when several libs for one patient are available.

	patient_bam=${patient_sorted%%sorted.dat}bam
	if needsUpdate ${patient_bam} ${patient_sorted} $MOSAIKBIN/MosaikText
	then
	    runme="$MOSAIKBIN/MosaikText -in $patient_sorted -bam $patient_bam"
	    vanillaRun "$runme" "$patient_bam" "result" "MosaikText -bam"
	fi

	patient_bcf=${patient_fastq_gz%%fastq.gz}var.raw.bcf
	if needsUpdate $patient_bcf $patient_bam $SAMTOOLS $BCFTOOLS
	then
	    runme="$SAMTOOLS mpileup -ugf $REFERENCE $patient_bam | $BCFTOOLS view -bvcg - > $patient_bcf"
	    vanillaRun "$runme" "$patient_bcf" "temp" "samtools mpileup - bcftools view"
	fi

	patient_vcf=${patient_bcf%%raw.bcf}flt.vcf
	if needsUpdate $patient_vcf $patient_bcf $BCFTOOLS $VCFUTILS
	then
	    if [ -z "$var_filter_settings" ]
	    then
		var_filter_settings="-D1000 -d6"
	    fi
	    runme="$BCFTOOLS view $patient_bcf | $VCFUTILS varFilter $var_filter_settings > $patient_vcf"
	    # -D200? 
	    vanillaRun "$runme" "$patient_vcf" "result" "bcftools view |vcfutils varFilter"
	fi
	
	patient_q20_vcf=${patient_vcf%%vcf}q20.vcf
	if needsUpdate $patient_q20_vcf $patient_vcf
	then
	    # generic baq-filter
	    awk '($6>=20) { print; }' < $patient_vcf > $patient_q20_vcf	    
	    checkExitStatus "Filter VCF to q20." "$patient_q20_vcf"
	    registerFile $patient_q20_vcf result
	fi

	patient_q20_avlist=${patient_q20_vcf%%vcf}avlist
	if needsUpdate $patient_q20_avlist $patient_q20_vcf $ANNOVARBIN/convert2annovar.pl
	then
	    runme="$ANNOVARBIN/convert2annovar.pl -format vcf4 $patient_q20_vcf > $patient_q20_avlist"
	    vanillaRun "$runme" "$patient_q20_avlist" "result" "convert2annovar q20"
	fi

	patient_avlist=${patient_vcf%%vcf}avlist
	if needsUpdate $patient_avlist $patient_vcf $ANNOVARBIN/convert2annovar.pl
	then
	    runme="$ANNOVARBIN/convert2annovar.pl -format vcf4 $patient_vcf > $patient_avlist"
	    vanillaRun "$runme" "$patient_avlist" "result" "convert2annovar"
	fi

	patient_exonic_variant=${patient_avlist}.exonic_variant_function
	if needsUpdate $patient_exonic_variant $patient_avlist $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --buildver hg19 $patient_avlist $AVDBDIR"
	    vanillaRun "$runme" "$patient_exonic_variant" "result" "ANNOVAR --geneanno"
	    
#	    $ANNOVARBIN/annotate_variation.pl --regionanno --buildver hg19 --dbtype dgv $patient_avlist $AVDBDIR
#	    export PATH=$PATH:"$ANNOVARBIN"
#	    $ANNOVARBIN/summarize_annovar.pl --buildver hg19 --verdbsnp 132 --outfile ${patient_fastq_gz%%.fastq.gz} $patient_avlist $AVDBDIR 
#	    registerFile other_annovar_files temp
	    
	fi

	releaseLock $patient_fastq_gz
    fi
done

: << 'POD_END'

=head1 DEPENDENCIES

Uses FastQC, MOSAIK, SAMTOOLS and ANNOVAR.

=head1 COPYRIGHT AND LICENSE

Copyright Daniel Nilsson, 2011. 

The package is released under the Perl Artistic License.

=head1 BUGS

Surely!

=cut

POD_END
    