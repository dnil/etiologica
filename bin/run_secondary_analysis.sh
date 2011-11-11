#!/bin/bash
# (c)2011 Daniel Nilsson, daniel.nilsson@scilifelab.se, daniel.nilsson@karolinska.se, daniel.nilsson@ki.se
# Copyright 2011, held by the author.
# Released under the Perl Artistic License.
#
# Documentation available in this file. Use e.g. perldoc to view, or pod2man to create a manual.
#
# TODO: decide on parameter config structures 
#

my_path=$_

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

### ENVIRONMENT VARIABLES

: <<POD_ENV

=head1 ENVIRONMENT VARIABLES

The pipeline adheres to environment exported variables.

=over 4

=item C<SNIC_TMP> I<path>

Designated temp area. Please make it local and sufficiently big for the sake of MOSAIK.
Can be overruled by directly setting TMP or MOSAIK_TMP, with the latter having the highest precedence rank.

=item C<BINDIR> I<path>

The path to the etiologia pipeline bin dir. Defaults to the (absolute) dir of the main script.

=item C<REFERENCE> I<genome.fasta.gz>

Reference genome sequence. Defaults to human_g1k_v27.fasta.gz. 

=item C<MOSAIKBIN> I<path>

Path to the Mosaik-aligner binary files. 

=item C<MOSAIK_CORES> I<int>

Number of cores to use for MOSAIK. 
Defaults to the number of available cores (see NPROC in the pipelinefunc library).

=item C<ANNOVARBIN> I<path>

Path to the annovar main scripts.

=back

=cut

POD_ENV

if [ -z "$BINDIR" ]
then 
    BINDIR=`dirname $my_path`
fi

PIPELINE=etiologica
PIPELINEFUNK=$BINDIR/pipelinefunk.sh

. $PIPELINEFUNK
log "Running $PIPELINE." "main"
log "BINDIR: $BINDIR" "main"

if [ -z "$REFERENCE" ]
then 
    REFERENCE=human_g1k_v37.fasta.gz
fi
log "REFERENCE: $REFERENCE" "main"

if [ -z "$SAMTOOLS" ]
then 
    SAMTOOLS=/home/daniel/src/samtools-0.1.17/samtools
fi
log "SAMTOOLS: $SAMTOOLS" "main"

if [ -z "$BCFTOOLS" ]
then 
    BCFTOOLS=/home/daniel/src/samtools-0.1.17/bcftools/bcftools
fi
log "BCFTOOLS: $BCFTOOLS" "main"

if [ -z "$FASTQC" ]
then 
    FASTQC=/home/daniel/src/FastQC/fastqc
fi
log "FASTQC: $FASTQC" "main"

if [ -z "$ANNOVARBIN" ]
then
    ANNOVARBIN=/home/daniel/src/annovar
fi
log "ANNOVARBIN: $ANNOVARBIN" "main"

if [ -z "$AVDBDIR" ]
then
    AVDBDIR=/home/daniel/src/annovar/humandb/
fi
log "ANNOVAR AVDBDIR: $AVDBDIR" "main"

if [ -z "$ANNOVAR_DISPENSABLE" ]
then
    ANNOVAR_DISPENSABLE=/home/daniel/src/annovar/example/dispensable.all
fi
log "ANNOVAR AVDBDIR: $AVDBDIR" "main"


if [ -z "$ANNOVAR_1KG_MAF" ]
then
    ANNOVAR_1KG_MAF=0.1
fi
log "ANNOVAR 1KG MAF: $ANNOVAR_1KG_MAF" "main"

if [ -z "$ANNOVAR_PP2_BENIGN" ]
then
    ANNOVAR_PP2_BENIGN=0.85
fi
log "ANNOVAR PP2 BENIGN: $ANNOVAR_PP2_BENIGN" "main"

if [ -z "$TMP" ] 
then
    if [ -z "$SNIC_TMP" ]
    then
	export TMP=/tmp
    else
	export TMP=$SNIC_TMP
    fi
fi
log "TMP: $TMP" "main"

# MOSAIK settings

if [ -z "$MOSAIKBIN" ]
then
    mosaik_aligner=`which MosaikAligner`
    if [ ! -z "$mosaik_aligner" ]
    then
	MOSAIKBIN=`dirname $mosaik_aligner`
    else
	MOSAIKBIN=/home/daniel/src/mosaik-aligner-read-only/bin
    fi
fi
log "MOSAIKBIN: $MOSAIKBIN" "main"

if [ -z "$MOSAIK_DAT_ON_SCRATCH" ]
then
    MOSAIK_DAT_ON_SCRATCH="no"
fi
log "MOSAIK_DAT_ON_SCRATCH: $MOSAIK_DAT_ON_SCRATCH" "main"

if [ -z "$MOSAIK_mismatches" ]
then 
    MOSAIK_mismatches=14
fi
log "MOSAIK - mismatches: $MOSAIK_mismatches" "main"

if [ -z "$MOSAIK_sw_bandwidth" ]
then
    MOSAIK_sw_bandwidth=33
fi
log "MOSAIK - SW bandwidth: $MOSAIK_sw_bandwidth" "main"

if [ -z "$MOSAIK_clustersize" ]
then
    MOSAIK_clustersize=35
fi
log "MOSAIK - clustersize: $MOSAIK_clustersize" "main"

if [ -z "$MOSAIK_JUMP" ]
then
    MOSAIK_JUMP="yes"
fi
log "MOSAIK - JUMP: $MOSAIK_JUMP" "main"

if [ -z "$MOSAIK_mjump" ]
then 
    MOSAIK_mjump=15
fi
log "MOSAIK - mjump: $MOSAIK_mjump" "main"

if [ -z $MOSAIK_mhp ]
then
    MOSAIK_mhp=100
fi
log "MOSAIK - mhp: $MOSAIK_mhp" "main"

if [ -z "$MOSAIK_TMP" ]
then 
    export MOSAIK_TMP=$TMP
fi 
log "MOSAIK - MOSAIK_TMP: $MOSAIK_TMP" "main"

if [ -z "$MOSAIK_CORES" ]
then
    if [ -z "$NPROC" ]
    then 
	MOSAIK_CORES=10
    else 
	MOSAIK_CORES=$NPROC
    fi
fi
log "MOSAIK - MOSAIK_CORES: $MOSAIK_CORES" "main"

# VCFUTILS (from samtools/bcftools)
if [ -z "$VCFUTILS" ]
then
    VCFUTILS=/home/daniel/src/samtools-0.1.17/bcftools/vcfutils.pl
fi
log "VCFUTILS: $VCFUTILS" "main"


if [ -z "$VCFUTILS_min_call_cov" ]
then
    VCFUTILS_min_call_cov=6
fi
log "VCFUTILS - min_call_cov $VCFUTILS_min_call_cov" "main"

if [ -z "$VCFUTILS_max_call_cov" ] 
then
    VCFUTILS_max_call_cov=400
fi
log "VCFUTILS - min_call_cov $VCFUTILS_max_call_cov" "main"

VCFUTILS_var_filter_settings="-D${VCFUTILS_max_call_cov} -d${VCFUTILS_min_call_cov}"

#log() all settings!
cd $BINDIR; 
pipeline_git_status=`git show --pretty=format:"%h %ci %s"; git status --porcelain`
cd -

log "pipeline git status: $pipeline_git_status" "main"

### end environment variable processing

for dir in patient group metadata
do 
    if [ ! -d "$dir" ]
    then
	mkdir $dir 
    fi
done

reference_dat=${REFERENCE%%.fasta.gz}.dat
reference_jump=${reference_dat%%.dat}.$MOSAIK_mjump

if needsUpdate $reference_dat $REFERENCE $MOSAIKBIN
then
    runme="$MOSAIKBIN/MosaikBuild -fr $REFERENCE -oa $reference_dat"
    vanillaRun "$runme" "$reference_dat" "common" "MosaikBuild"

    if [ "$JUMP" == "yes" ]
    then
	runme="$MOSAIKBIN/MosaikJump -ia $reference_dat -hs $MOSAIK_mjump -out $reference_jump -mhp $MOSAIK_mhp"
	vanillaRun "$runme" "$reference_jump" "common" "MosaikJump"
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
	    log "Locked $patient_fastqc_zip for FastQC." "main"
	    runme="$FASTQC $patientfastq"
	    vanillaRun "$runme" "$patientfastqc_zip" "result" "FastQC"

	    releaseLock $patientfastqc_zip
	    log "Workblock exit: released $patient_fastqc_zip lock for FastQC." "main"
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
	    echo "Did not find any ($possible_mate) mate for $patient_fasta_gz. Assuming all reads are unpaired."

	    MATEPAIRS=0
	fi
    done

    for patient_fastq_gz in patient/*_2.fastq.gz 
    do
	possible_mate=${patient_fastq_gz%%_2.fastq.gz}_1.fastq.gz
	if [ ! -s $possible_mate ] 
	then
	    # possible mate according to convention not found.
	    echo "Did not find any ($possible_mate) mate for $patient_fasta_gz. Assuming all reads are unpaired."

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
	log "Locked $patient_fastq_gz for alignment, variant calling and annotation." "main"

#	patient_trim_fastq_gz=${patient_fastq_gz%%fastq.gz}trim.fastq.gz
#	if needsUpdate $patient_trim_fastq_gz $patient_fastq_gz $BINDIR/trim.pl
#	then
#	    runme="zcat $patient_fastq_gz | $BINDIR/trim.pl|$PGZIP -c > $patient_trim_fastq_gz"
#	    vanillaRun "$runme" "$patient_trim_fastq_gz" "temp" "Trim read ends"
#	fi
	# trim

	patient_dat=${patient_fastq_gz%%fastq.gz}dat

	# If mosaik_dat_on_scratch is enabled, we only really want to redo all the dat files in case the final product mosaik bam is outdated. Yay, let us reinvent make!
	run_align="yes"
	if [ $MOSAIK_DAT_ON_SCRATCH == "yes" ]
	then
	    mkdir -p $MOSAIK_TMP/`dirname $patient_dat`
	    patient_dat=$MOSAIK_TMP/$patient_dat

	    final_patient_bam=${patient_dat%%dat}mosaik.bam
	    if needsUpdate ${final_patient_bam} ${patient_fastq_gz} $MOSAIKBIN/MosaikBuild $MOSAIKBIN/MosaikAligner $MOSAIKBIN/MosaikText $reference_jump $MOSAIKBIN/MosaikDupSnoop $MOSAIKBIN/MosaikSort
	    then
		run_align="yes"
	    else
		run_align="no"
	    fi
	fi

	if [ $run_align == "yes" ] 
	then
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
		    runme="$MOSAIKBIN/MosaikAligner -in $patient_dat -ia $reference_dat -out $patient_aln_dat -m $MOSAIK_ALIGN_MODE -hs $MOSAIK_mjump -bw $MOSAIK_sw_bandwidth -j $reference_jump -mhp $MOSAIK_mhp -mm $MOSAIK_mismatches -act $MOSAIK_clustersize -p $MOSAIK_CORES"
		else
		    runme="$MOSAIKBIN/MosaikAligner -in $patient_dat -ia $reference_dat -out $patient_aln_dat -m $MOSAIK_ALIGN_MODE -bw $MOSAIK_sw_bandwidth -mm $MOSAIK_mismatches -act $MOSAIK_clustersize -p $MOSAIK_CORES"
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

	    if [ $MOSAIK_DAT_ON_SCRATCH == "yes" ] 
	    then
		patient_bam=${patient_bam##${MOSAIK_TMP}/}
	    fi

	    if needsUpdate ${patient_bam} ${patient_sorted} $MOSAIKBIN/MosaikText
	    then
		runme="$MOSAIKBIN/MosaikText -in $patient_sorted -bam $patient_bam"
		vanillaRun "$runme" "$patient_bam" "result" "MosaikText -bam"
	    fi
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
	    runme="$BCFTOOLS view $patient_bcf | $VCFUTILS varFilter $VCFUTILS_var_filter_settings > $patient_vcf"
	    vanillaRun "$runme" "$patient_vcf" "result" "bcftools view |vcfutils varFilter"
	fi
	
	patient_pass_vcf=${patient_vcf%%vcf}pass.vcf
	if needsUpdate $patient_pass_vcf $patient_vcf
	then
	    # generic baq-filter : NB need to avoid $6 being evaluated already att passing... PASS only?
	    # note: mpileup does not give PASS, only UnifiedGenotyper.. Go GATK.
	    runme="awk '(\$6>=20) { print; }' < $patient_vcf > $patient_pass_vcf"
#	    runme="awk '(\$7==\"PASS\") { print; }' < $patient_vcf > $patient_pass_vcf"
	    vanillaRun "$runme" "$patient_pass_vcf" "result" "Filter VCF to pass."
	fi
	
	patient_avlist=${patient_vcf%%vcf}avlist
	if needsUpdate $patient_avlist $patient_vcf $ANNOVARBIN/convert2annovar.pl
	then
	    runme="$ANNOVARBIN/convert2annovar.pl -allallele -format vcf4 $patient_vcf > $patient_avlist"
	    vanillaRun "$runme" "$patient_avlist" "result" "convert2annovar"
	fi

	patient_exonic_variant=${patient_avlist}.exonic_variant_function
	if needsUpdate $patient_exonic_variant $patient_avlist $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --buildver hg19 $patient_avlist $AVDBDIR"
	    vanillaRun "$runme" "$patient_exonic_variant" "result" "ANNOVAR --geneanno"
	fi

	patient_pass_avlist=${patient_pass_vcf%%vcf}avlist
	if needsUpdate $patient_pass_avlist $patient_pass_vcf $ANNOVARBIN/convert2annovar.pl
	then
	    runme="$ANNOVARBIN/convert2annovar.pl -allallele -format vcf4 $patient_pass_vcf > $patient_pass_avlist"
	    vanillaRun "$runme" "$patient_pass_avlist" "result" "convert2annovar pass"
	fi

	patient_exonic_variant=${patient_pass_avlist}.exonic_variant_function
	if needsUpdate $patient_exonic_variant $patient_pass_avlist $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --buildver hg19 $patient_pass_avlist $AVDBDIR"
	    vanillaRun "$runme" "$patient_exonic_variant" "result" "ANNOVAR --geneanno"
	fi

	patient_maf_filter=${patient_pass_avlist}.hg19_ALL.sites.2010_11_filtered
	if needsUpdate $patient_maf_filter $patient_pass_avlist $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --filter --dbtype 1000g2010nov_all --maf_threshold $ANNOVAR_1KG_MAF --buildver hg19 $patient_pass_avlist $AVDBDIR"
	    vanillaRun "$runme" "$patient_maf_filter" "result" "ANNOVAR --filter 1000g maf $ANNOVAR_1KG_MAF"
	fi

	patient_maf_exonic_variant=${patient_maf_filter}.exonic_variant_function
	if needsUpdate $patient_maf_exonic_variant $patient_maf_filter $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --buildver hg19 $patient_maf_filter $AVDBDIR"
	    vanillaRun "$runme" "$patient_maf_exonic_variant" "result" "ANNOVAR 1000g MAF $ANNOVAR_1KG_MAF --geneanno"
	fi

	patient_pp2_filter=${patient_maf_filter}.hg19_ljb_pp2_filtered
	if needsUpdate $patient_pp2_filter $patient_maf_filter $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --filter --dbtype ljb_pp2 --reverse --score_threshold $ANNOVAR_PP2_BENIGN --buildver hg19 $patient_maf_filter $AVDBDIR"
	    vanillaRun "$runme" "$patient_pp2_filter" "result" "ANNOVAR --filter PolyPhen2 benign level $ANNOVAR_PP2_BENIGN"
	fi

	patient_pp2_exonic_variant=${patient_pp2_filter}.exonic_variant_function
	if needsUpdate $patient_pp2_exonic_variant $patient_pp2_filter $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --buildver hg19 $patient_pp2_filter $AVDBDIR"
	    vanillaRun "$runme" "$patient_pp2_exonic_variant" "result" "ANNOVAR PolyPhen2 benign $ANNOVAR_PP2_BENIGN --geneanno"
	fi
	
	patient_pp2_variant_function=${patient_pp2_filter}.variant_function
	patient_splicing=${patient_pp2_variant_function}.splicing
	if needsUpdate $patient_splicing $patient_pp2_variant_function
	then
	    runme="awk '(\$1 == \"exonic;splicing\" || \$1 == \"splicing\") { print }' $patient_pp2_variant_function > $patient_splicing"
	    vanillaRun "$runme" "$patient_splicing" "result" "Extract splicing variants."
	fi	    

        # gather exonic;splicing from variant_function
	patient_all_variant_function=${patient_pass_avlist}.variant_function
	patient_all_splicing=${patient_all_variant_function}.splicing
	if needsUpdate $patient_all_splicing $patient_all_variant_function
	then
	    runme="awk '(\$1 == \"exonic;splicing\" || \$1 == \"splicing\") { print }' $patient_all_variant_function > $patient_all_splicing"
	    vanillaRun "$runme" "$patient_all_splicing" "result" "Extract passed splicing variants."
	fi

	patient_not_syn=${patient_pp2_exonic_variant}.not_syn
	if needsUpdate $patient_not_syn $patient_pp2_exonic_variant
	then
	    runme="cut -f2- $patient_pp2_exonic_variant | awk '(\$1 != \"synonymous\") { print }' > $patient_not_syn"
	    vanillaRun "$runme" "$patient_not_syn" "result" "Ignore synonymous variants for now."
	fi

	# hom | 2xhet (remember splicing!)	
	patient_recessive_genelist=${patient_pp2_variant_function}.recessive_model
	if needsUpdate $patient_recessive_genelist $patient_not_syn $patient_splicing
	then
	   runme="cat ${patient_not_syn} $patient_splicing |perl -e 'while(<STDIN>) { chomp; @r=split(/\t/); @transcripts=split(/[:;(]+/,\$r[1]); \$gene=\$transcripts[0]; \$gene_row=\$gene.\"\t\".join(\"\t\",@r).\"\n\"; \$nvars{\$gene}++; \$vars{\$gene}.=\$gene_row; if((\$r[7] eq \"hom\") || (\$r[7] eq \"het\" && \$nvars{\$gene} >1) ) { \$modelok{\$gene}=1; } } foreach \$gene (keys %modelok) {print \$vars{\$gene}};' > $patient_recessive_genelist" 
	   vanillaRun "$runme" "$patient_recessive_genelist" "result" "Produce recessive model gene list"
	fi

	patient_indisp=${patient_recessive_genelist}.indisp
	if needsUpdate $patient_indisp $patient_recessive_genelist
	then
	  runme="grep -v -w -f $ANNOVAR_DISPENSABLE $patient_recessive_genelist > $patient_indisp"
	  vanillaRun "$runme" "$patient_indisp" "result" "Indispensable genes on recessive gene list."
	fi

# sift
#	patient_sift_filter=${patient_maf_filter}.hg19_ALL.sites.2010_11_filtered
#	if needsUpdate $patient_sift_filter $patient_avlistOA $ANNOVARBIN/annotate_variation.pl
#	then
#	    runme="$ANNOVARBIN/annotate_variation.pl --filter --dbtype avsift --buildver hg19 $patient_maf_filter $AVDBDIR"
#	    vanillaRun "$runme" "$patient_maf_filter" "result" "ANNOVAR --filter 1000g maf $ANNOVAR_1KG_MAF"
#	fi
# structural variants
#	    $ANNOVARBIN/annotate_variation.pl --regionanno --buildver hg19 --dbtype dgv $patient_avlist $AVDBDIR
# to csv, many sources --- slow!
#	    export PATH=$PATH:"$ANNOVARBIN"
#	    $ANNOVARBIN/summarize_annovar.pl --buildver hg19 --verdbsnp 132 --outfile ${patient_fastq_gz%%.fastq.gz} $patient_avlist $AVDBDIR 
#	    registerFile other_annovar_files temp
#
        # cnv analysis
	
	releaseLock $patient_fastq_gz
	log "Workblock exit: released lock on $patient_fastq_gz." "main"
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
    