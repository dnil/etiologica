#!/bin/bash
# (c)2011 Daniel Nilsson, daniel.nilsson@scilifelab.se, daniel.nilsson@karolinska.se, daniel.nilsson@ki.se
# Copyright 2011, held by the author.
# Released under the Perl Artistic License.
#
# Documentation available in this file. Use e.g. perldoc to view, or pod2man to create a manual.

my_path=$_

### BEGIN MAIN DOCUMENTATION

: << 'POD_INIT'

=head1 NAME

run_tertiary_analysis.sh

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@scilifelab.se, daniel.nilsson@ki.se, daniel.nilsson@karolinska.se

=head1 SYNOPSIS

 mkdir myproject
 ln -s patientid.fastq myproject/
 cd myproject
 run_secondary_analysis.sh
#  iro iro attanda
 run_tertiary_analysis.sh

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

=item C<BINDIR> I<path>

The path to the etiologia pipeline bin dir. Defaults to the (absolute) dir of the main script.

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

if [ -z "$VCFUTILS" ]
then
    VCFUTILS=/home/daniel/src/samtools-0.1.18/bcftools/vcfutils.pl
fi
log "VCFUTILS: $VCFUTILS" "main"

if [ -z "$ANNOVARBIN" ]
then
    ANNOVARBIN=/home/daniel/src/annovar
fi
log "ANNOVARBIN: $ANNOVARBIN" "main"

if [ -z "$ANNOVAR_SUMMARIZE" ]
then 
    ANNOVAR_SUMMARIZE="$ANNOVARBIN/summarize_annovar_custom.pl"
fi

if [ -z "$AVDBDIR" ]
then
    AVDBDIR=/home/daniel/src/annovar/humandb/
fi
log "ANNOVAR AVDBDIR: $AVDBDIR" "main"

if [ -z "$MAX_MAF" ]
then
    MAX_MAF=0.02
fi

if [ -z "$ANNOVAR_1KG_MAF" ]
then
    ANNOVAR_1KG_MAF=$MAX_MAF
fi
log "ANNOVAR 1KG_MAF: $ANNOVAR_1KG_MAF" "main"

if [ -z "$LOCAL_MAF" ]
then
    LOCAL_MAF=$MAX_MAF
fi
log "LOCAL DB MAF: $LOCAL_MAF" "main"

if [ -z "$ANNOVAR_DISPENSABLE" ]
then
    ANNOVAR_DISPENSABLE=/home/daniel/src/annovar/example/dispensable.all
fi
log "ANNOVAR AVDBDIR: $AVDBDIR" "main"

LOCAL_DB_FILTER="yes"
if [ -z "$LOCAL_CLIN_DB" ] 
then 
    LOCAL_DB_FILTER="no"
fi 

if [ -z "$LOCAL_DANES_DB" ] 
then
    LOCAL_DB_FILTER="no"
fi

if [ -z "$DB_SNP_VERSION" ]
then
    DB_SNP_VERSION="snp137NonFlagged"
fi

if [ -z "$DB_SNP_VERSION_ANNOTATE" ]
then
    DB_SNP_VERSION_ANNOTATE="137"
fi

if [ -z "$DB_ESP_VERSION" ]
then
    DB_ESP_VERSION="6500"
fi

if [ -z "$KGENOMES_VERSION" ]
then     
    KGENOMES_VERSION="hg19_ALL.sites.2012_04"
    KGENOMES_DB="1000g2012apr"
fi

if [ -z "$ANNOVAR_SPLICING_THRESHOLD" ]
then
    ANNOVAR_SPLICING_THRESHOLD=5
fi

if [ -z "$PASS_QUAL" ]
then
    PASS_QUAL=25
fi

if [ -z "$BCFTOOLS" ]
then 
    BCFTOOLS=/home/daniel/src/samtools-0.1.18/bcftools/bcftools
fi
log "BCFTOOLS: $BCFTOOLS" "main"

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

# VCFTOOLS?
#log() all settings!
cd $BINDIR; 
pipeline_git_status=`git show --pretty=format:"%h %ci %s"; git status --porcelain`
cd -

log "pipeline git status: $pipeline_git_status" "main"

### end environment variable processing

patient_vcf_list=( patient/*flt.vcf )

for patient_vcf in ${patient_vcf_list[@]}
do 
    # 
    # File locking for multiple script instances concurrently in same directory.
    # One instance per patient fastq.gz is a rather appropriate work block right now - if 
    # multilane samples start popping in this recommendation may change.
    #

    if workLockOk $patient_vcf
    then
	log "Locked $patient_vcf for annotation." "main"

	# generic baq-filter : NB need to avoid $6 being evaluated already att passing... PASS only?
	# note: mpileup does not give PASS, only UnifiedGenotyper.. Go GATK.

	patient_left_vcf=${patient_vcf%%.vcf.gz}.leftAlign.vcf
	patient_pass_vcf=${patient_left_vcf%%vcf}pass.vcf
	if needsUpdate $patient_pass_vcf $patient_left_vcf
	then
            runme="$BCFTOOLS filter -O z -o $patient_pass_vcf -s LOWQUAL -i'%QUAL>$PASS_QUAL' $patient_left_vcf"
	    vanillaRun "$runme" "$patient_pass_vcf" "result" "Filter VCF to pass."
	fi
	
	patient_pass_vcf_nogz=${patient_pass_vcf%%.gz}
	patient_avlist=${patient_pass_vcf%%vcf}avlist
	if needsUpdate $patient_avlist $patient_pass_vcf $ANNOVARBIN/convert2annovar.pl
	then
	    
	    runme="gunzip -c $patient_pass_vcf > ${patient_pass_vcf_nogz}"	    
	    vanillaRun "$runme" "$patient_pass_vcf_nogz" "result" "gunzip vcf"

	    runme="$ANNOVARBIN/convert2annovar.pl -format vcf4 $patient_pass_vcf_nogz > $patient_avlist"
	    vanillaRun "$runme" "$patient_avlist" "result" "convert2annovar vcf"
	fi

	patient_annovar_summarize_base=${patient_avlist%%avlist}sum
	patient_annovar_summarize_csv=${patient_avlist%%avlist}sum.genome_summary.csv
	if needsUpdate $patient_annovar_summarize_csv $patient_avlist $ANNOVAR_SUMMARIZE
	then
	    runme="$TABLE_ANNOVAR --buildver hg19 -vcfdbfile $LOCAL_CLIN_DB -protocol refGene,phastConsElements46way,genomicSuperDups,popfreq_max,exac02,esp6500si_all,1000g2012apr_all,vcf,snp138,cosmic,caddgt10,ljb2_all,clinvar_20131105 -operation g,r,r,f,f,f,f,f,f,f,f,f,f -remove -otherinfo -csvout $AVDBDIR"
	    vanillaRun "$runme" "$patient_annovar_summarize_csv" "result" "ANNOVAR SUMMARIZE"
	fi
	# --localdb $LOCAL_CLIN_DB (vcf..), ExAC, cosmic, ... 

	# inheritance models

	releaseLock $patient_vcf
	log "Workblock exit: released lock on $patient_vcf." "main"
    fi
done

# cnv analysis...bam file to variant vcf in 2ndary, filtering in 3rtiary?

: << 'POD_END'

=head1 DEPENDENCIES

Uses ANNOVAR.

=head1 COPYRIGHT AND LICENSE

Copyright Daniel Nilsson, 2011. 

The package is released under the Perl Artistic License.

=head1 BUGS

Surely!

=cut

POD_END


# # sift
# #	patient_sift_filter=${patient_maf_filter}.hg19_ALL.sites.2010_11_filtered
# #	if needsUpdate $patient_sift_filter $patient_avlistOA $ANNOVARBIN/annotate_variation.pl
# #	then
# #	    runme="$ANNOVARBIN/annotate_variation.pl --filter --dbtype avsift --buildver hg19 $patient_maf_filter $AVDBDIR"
# #	    vanillaRun "$runme" "$patient_maf_filter" "result" "ANNOVAR --filter 1000g maf $ANNOVAR_1KG_MAF"
# #	fi
# # structural variants
# #	    $ANNOVARBIN/annotate_variation.pl --regionanno --buildver hg19 --dbtype dgv $patient_avlist $AVDBDIR
# # to csv, many sources --- slow!
# #	    export PATH=$PATH:"$ANNOVARBIN"
# #	    $ANNOVARBIN/summarize_annovar.pl --buildver hg19 --verdbsnp 132 --outfile ${patient_fastq_gz%%.fastq.gz} $patient_avlist $AVDBDIR 
# #	    registerFile other_annovar_files temp
# #
#         # cnv analysis
    