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
    VCFUTILS=/home/daniel/src/samtools-0.1.17/bcftools/vcfutils.pl
fi
log "VCFUTILS: $VCFUTILS" "main"

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

if [ -z "$ANNOVAR_1KG_MAFR" ]
then
    ANNOVAR_1KG_MAF=0.02
fi
log "ANNOVAR 1KG_MAF: $ANNOVAR_1KG_MAF" "main"

if [ -z "$ANNOVAR_DISPENSABLE" ]
then
    ANNOVAR_DISPENSABLE=/home/daniel/src/annovar/example/dispensable.all
fi
log "ANNOVAR AVDBDIR: $AVDBDIR" "main"

if [ -z "$ANNOVAR_PP2_BENIGN" ]
then
    ANNOVAR_PP2_BENIGN=0.85
fi
log "ANNOVAR PP2 BENIGN: $ANNOVAR_PP2_BENIGN" "main"

if [ -z "$ANNOVAR_SPLICING_THRESHOLD" ]
then
    ANNOVAR_SPLICING_THRESHOLD=5
fi

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

	patient_pass_vcf=${patient_vcf%%vcf}pass.vcf
	if needsUpdate $patient_pass_vcf $patient_vcf
	then
	    # generic baq-filter : NB need to avoid $6 being evaluated already att passing... PASS only?
	    # note: mpileup does not give PASS, only UnifiedGenotyper.. Go GATK.
	    runme="awk '(\$6>=30) { print; }' < $patient_vcf > $patient_pass_vcf"
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
	    runme="$ANNOVARBIN/annotate_variation.pl --splicing_threshold $ANNOVAR_SPLICING_THRESHOLD --buildver hg19 $patient_avlist $AVDBDIR"
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
	    runme="$ANNOVARBIN/annotate_variation.pl --splicing_threshold $ANNOVAR_SPLICING_THRESHOLD --buildver hg19 $patient_pass_avlist $AVDBDIR"
	    vanillaRun "$runme" "$patient_exonic_variant" "result" "ANNOVAR --geneanno"
	fi

        # gather exonic;spilcing from variant_function
	patient_all_variant_function=${patient_pass_avlist}.variant_function
	patient_all_splicing=${patient_all_variant_function}.splicing
	if needsUpdate $patient_all_splicing $patient_all_variant_function
	then
	    runme="awk '(\$1 == \"exonic;splicing\" || \$1 == \"splicing\") { print }' $patient_all_variant_function > $patient_all_splicing"
	    vanillaRun "$runme" "$patient_all_splicing" "result" "Extract passed splicing variants."
	fi

	patient_rare_filter=${patient_pass_avlist}.rare
	patient_maf_filter=${patient_pass_avlist}.hg19_ALL.sites.2010_11_filtered
	if needsUpdate $patient_maf_filter $patient_pass_avlist $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --filter --dbtype 1000g2010nov_all --maf_threshold $ANNOVAR_1KG_MAF --buildver hg19 $patient_pass_avlist $AVDBDIR"
	    vanillaRun "$runme" "$patient_maf_filter" "result" "ANNOVAR --filter 1000g maf $ANNOVAR_1KG_MAF"
	fi

	patient_100clinical_filter_orig=${patient_maf_filter}.hg19_generic_filtered
	patient_100clinical_filter=${patient_maf_filter}.hg19_100clinical_filtered
	if needsUpdate $patient_100clinical_filter $patient_maf_filter $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --filter --dbtype generic --genericdbfile hg19_100clinical.avdb --buildver hg19 $patient_maf_filter $AVDBDIR"
	    vanillaRun "$runme" "$patient_100clinical_filter" "result" "ANNOVAR --filter 100clinical"
	    mv $patient_100clinical_filter_orig $patient_100clinical_filter
	fi

	patient_danes_filter_orig=${patient_100clinical_filter}.hg19_generic_filtered
	patient_danes_filter=${patient_100clinical_filter}.hg19_200danes_filtered
	if needsUpdate $patient_danes_filter $patient_clinical_filter $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --filter --dbtype generic --genericdbfile hg19_200danes.avdb --buildver hg19 $patient_100clinical_filter $AVDBDIR"
	    vanillaRun "$runme" "$patient_danes_filter" "result" "ANNOVAR --filter 200danes"
	    mv $patient_danes_filter_orig $patient_danes_filter
	fi
	
	if needsUpdate $patient_rare_filter $patient_danes_filter
	then
	    cp $patient_danes_filter $patient_rare_filter
	fi

	patient_rare_exonic_variant=${patient_rare_filter}.exonic_variant_function
	if needsUpdate $patient_rare_exonic_variant $patient_rare_filter $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --geneanno --splicing_threshold $ANNOVAR_SPLICING_THRESHOLD --buildver hg19 $patient_rare_filter $AVDBDIR"
	    vanillaRun "$runme" "$patient_rare_exonic_variant" "result" "ANNOVAR 1000g MAF $ANNOVAR_1KG_MAF and 200danes --geneanno"
	fi

	patient_uncommon_not_syn=${patient_rare_exonic_variant}.not_syn
	if needsUpdate $patient_uncommon_not_syn $patient_rare_exonic_variant
	then
	    runme="cut -f2- $patient_rare_exonic_variant | awk '(\$1 != \"synonymous\") { print }' > $patient_uncommon_not_syn"
	    vanillaRun "$runme" "$patient_uncommon_not_syn" "result" "Ignore synonymous variants without any other severity prediction."
	fi

	patient_pp2_filter=${patient_rare_filter}.hg19_ljb_pp2_filtered
	if needsUpdate $patient_pp2_filter $patient_rare_filter $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --filter --dbtype ljb_pp2 --reverse --score_threshold $ANNOVAR_PP2_BENIGN --buildver hg19 $patient_rare_filter $AVDBDIR"
	    vanillaRun "$runme" "$patient_pp2_filter" "result" "ANNOVAR --filter PolyPhen2 benign level $ANNOVAR_PP2_BENIGN"
	fi

	patient_pp2_exonic_variant=${patient_pp2_filter}.exonic_variant_function
	if needsUpdate $patient_pp2_exonic_variant $patient_pp2_filter $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --splicing_threshold $ANNOVAR_SPLICING_THRESHOLD --buildver hg19 $patient_pp2_filter $AVDBDIR"
	    vanillaRun "$runme" "$patient_pp2_exonic_variant" "result" "ANNOVAR PolyPhen2 benign $ANNOVAR_PP2_BENIGN --geneanno"
	fi
	
	patient_pp2_variant_function=${patient_pp2_filter}.variant_function
	patient_pp2_splicing=${patient_pp2_variant_function}.splicing
	if needsUpdate $patient_pp2_splicing $patient_pp2_variant_function
	then
	    runme="awk '(\$1 == \"exonic;splicing\" || \$1 == \"splicing\") { print }' $patient_pp2_variant_function > $patient_pp2_splicing"
	    vanillaRun "$runme" "$patient_pp2_splicing" "result" "Extract splicing variants."
	fi
        
	patient_not_syn=${patient_pp2_exonic_variant}.not_syn
	if needsUpdate $patient_not_syn $patient_pp2_exonic_variant
	then
	    runme="cut -f2- $patient_pp2_exonic_variant | awk '(\$1 != \"synonymous\") { print }' > $patient_not_syn"
	    vanillaRun "$runme" "$patient_not_syn" "result" "Ignore synonymous variants for now."
	fi

	patient_pp2_splicing_indisp=${patient_pp2_splicing}.indisp
	if needsUpdate $patient_pp2_splicing_indisp $patient_pp2_splicing $ANNOVAR_DISPENSABLE
	then
	  runme="grep -v -w -f $ANNOVAR_DISPENSABLE $patient_pp2_splicing > $patient_pp2_splicing_indisp"
	  vanillaRun "$runme" "$patient_pp2_splicing_indisp" "result" "Filter splicing variants in dispensable genes."
	fi

	patient_not_syn_indisp=${patient_not_syn}.indisp
	if needsUpdate $patient_not_syn_indisp $patient_not_syn $ANNOVAR_DISPENSABLE
	then
	  runme="grep -v -w -f $ANNOVAR_DISPENSABLE $patient_not_syn > $patient_not_syn_indisp"
	  vanillaRun "$runme" "$patient_not_syn_indisp" "result" "Filter variants in dispensable genes."
	fi
	
	# hom | 2xhet (remember splicing!)	
	patient_recessive_variant_list=${patient_pp2_variant_function}.indisp.recessive_model
	if needsUpdate $patient_recessive_variant_list $patient_not_syn_indisp $patient_pp2_splicing_indisp
	then
	   runme="cat ${patient_not_syn_indisp} $patient_pp2_splicing_indisp |perl -e 'while(<STDIN>) { chomp; @r=split(/\t+/); @transcripts=split(/[:;(]+/,\$r[1]); \$gene=\$transcripts[0]; \$gene=~s/\)//; \$gene_row=\$gene.\"\t\".join(\"\t\",@r).\"\n\"; \$nvars{\$gene}++; \$vars{\$gene}.=\$gene_row; if((\$r[7] eq \"hom\") || (\$r[7] eq \"het\" && \$nvars{\$gene} >1) ) { \$modelok{\$gene}=1; } } foreach \$gene (keys %modelok) {print \$vars{\$gene}};' > $patient_recessive_variant_list" 
	   vanillaRun "$runme" "$patient_recessive_variant_list" "result" "Produce recessive model variant list"
	fi
	
	# cut recessive model output to avlist
	patient_recessive_avlist=${patient_recessive_variant_list}.avlist
	if needsUpdate $patient_recessive_avlist $patient_recessive_variant_list
	then
	    runme="perl -ne 'chomp; @r=split(/\t+/); print join(\"\t\",@r[3..10], @r[0..2]),\"\n\";' $patient_recessive_variant_list > $patient_recessive_avlist"
	    vanillaRun "$runme" "$patient_recessive_avlist" "temp" "Produce recessive model avlist"
	fi

	# RECESSIVE: filter against dbSNP - will retain for most apps; regarded mostly as an rsID annotation step.
	patient_recessive_dbSNP_filter=${patient_recessive_avlist}.hg19_snp132_filtered
	if needsUpdate $patient_recessive_dbSNP_filter $patient_recessive_avlist $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --filter --dbtype snp132 --buildver hg19 $patient_recessive_avlist $AVDBDIR"
	    vanillaRun "$runme" "$patient_recessive_dbSNP_filter" "result" "ANNOVAR --filter dbSNP132"
	fi

	# DOMINANT: filter against dbSNP - will retain for most apps; regarded mostly as an rsID annotation step.
	patient_dbSNP_filter=${patient_pp2_filter}.hg19_snp132_filtered
	if needsUpdate $patient_dbSNP_filter $patient_pp2_filter $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --filter --dbtype snp132 --buildver hg19 $patient_pp2_filter $AVDBDIR"
	    vanillaRun "$runme" "$patient_dbSNP_filter" "result" "ANNOVAR --filter dbSNP132"
	fi
	
	# annotate... 
	patient_dbSNP_exonic_variant=${patient_dbSNP_filter}.exonic_variant_function
	if needsUpdate $patient_dbSNP_exonic_variant $patient_dbSNP_filter $ANNOVARBIN/annotate_variation.pl
	then
	    runme="$ANNOVARBIN/annotate_variation.pl --splicing_threshold $ANNOVAR_SPLICING_THRESHOLD --buildver hg19 $patient_dbSNP_filter $AVDBDIR"
	    vanillaRun "$runme" "$patient_dbSNP_exonic_variant" "result" "ANNOVAR PolyPhen2 & snp132-filtered --geneanno"
	fi
	
	patient_dbSNP_variant_function=${patient_dbSNP_filter}.variant_function
	patient_dbSNP_splicing=${patient_dbSNP_variant_function}.splicing
	if needsUpdate $patient_dbSNP_splicing $patient_dbSNP_variant_function
	then
	    runme="awk '(\$1 == \"exonic;splicing\" || \$1 == \"splicing\") { print }' $patient_dbSNP_variant_function > $patient_dbSNP_splicing"
	    vanillaRun "$runme" "$patient_dbSNP_splicing" "result" "Extract splicing variants."
	fi
        
	patient_not_syn=${patient_dbSNP_exonic_variant}.not_syn
	if needsUpdate $patient_not_syn $patient_dbSNP_exonic_variant
	then
	    runme="cut -f2- $patient_dbSNP_exonic_variant | awk '(\$1 != \"synonymous\") { print }' > $patient_not_syn"
	    vanillaRun "$runme" "$patient_not_syn" "result" "Ignore synonymous variants for now."
	fi

	patient_dbSNP_splicing_indisp=${patient_dbSNP_splicing}.indisp
	if needsUpdate $patient_dbSNP_splicing_indisp $patient_dbSNP_splicing $ANNOVAR_DISPENSABLE
	then
	  runme="grep -v -w -f $ANNOVAR_DISPENSABLE $patient_dbSNP_splicing > $patient_dbSNP_splicing_indisp"
	  vanillaRun "$runme" "$patient_dbSNP_splicing_indisp" "result" "Filter splicing variants in dispensable genes."
	fi

	patient_not_syn_indisp=${patient_not_syn}.indisp
	if needsUpdate $patient_not_syn_indisp $patient_not_syn $ANNOVAR_DISPENSABLE
	then
	  runme="grep -v -w -f $ANNOVAR_DISPENSABLE $patient_not_syn > $patient_not_syn_indisp"
	  vanillaRun "$runme" "$patient_not_syn_indisp" "result" "Filter variants in dispensable genes."
	fi

       # cnv analysis...
	
	# variants in common, genes in common, both for dominant and recessive

	releaseLock $patient_vcf
	log "Workblock exit: released lock on $patient_vcf." "main"
    fi
done

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



# 	patient_maf_filter=${patient_pass_avlist}.hg19_ALL.sites.2010_11_filtered
# 	if needsUpdate $patient_maf_filter $patient_pass_avlist $ANNOVARBIN/annotate_variation.pl
# 	then
# 	    runme="$ANNOVARBIN/annotate_variation.pl --filter --dbtype 1000g2010nov_all --maf_threshold $ANNOVAR_1KG_MAF --buildver hg19 $patient_pass_avlist $AVDBDIR"
# 	    vanillaRun "$runme" "$patient_maf_filter" "result" "ANNOVAR --filter 1000g maf $ANNOVAR_1KG_MAF"
# 	fi

# 	patient_maf_exonic_variant=${patient_maf_filter}.exonic_variant_function
# 	if needsUpdate $patient_maf_exonic_variant $patient_maf_filter $ANNOVARBIN/annotate_variation.pl
# 	then
# 	    runme="$ANNOVARBIN/annotate_variation.pl --buildver hg19 $patient_maf_filter $AVDBDIR"
# 	    vanillaRun "$runme" "$patient_maf_exonic_variant" "result" "ANNOVAR 1000g MAF $ANNOVAR_1KG_MAF --geneanno"
# 	fi

# 	patient_pp2_filter=${patient_maf_filter}.hg19_ljb_pp2_filtered
# 	if needsUpdate $patient_pp2_filter $patient_maf_filter $ANNOVARBIN/annotate_variation.pl
# 	then
# 	    runme="$ANNOVARBIN/annotate_variation.pl --filter --dbtype ljb_pp2 --reverse --score_threshold $ANNOVAR_PP2_BENIGN --buildver hg19 $patient_maf_filter $AVDBDIR"
# 	    vanillaRun "$runme" "$patient_pp2_filter" "result" "ANNOVAR --filter PolyPhen2 benign level $ANNOVAR_PP2_BENIGN"
# 	fi

# 	patient_pp2_exonic_variant=${patient_pp2_filter}.exonic_variant_function
# 	if needsUpdate $patient_pp2_exonic_variant $patient_pp2_filter $ANNOVARBIN/annotate_variation.pl
# 	then
# 	    runme="$ANNOVARBIN/annotate_variation.pl --buildver hg19 $patient_pp2_filter $AVDBDIR"
# 	    vanillaRun "$runme" "$patient_pp2_exonic_variant" "result" "ANNOVAR PolyPhen2 benign $ANNOVAR_PP2_BENIGN --geneanno"
# 	fi
	
# 	patient_pp2_variant_function=${patient_pp2_filter}.variant_function
# 	patient_splicing=${patient_pp2_variant_function}.splicing
# 	if needsUpdate $patient_splicing $patient_pp2_variant_function
# 	then
# 	    runme="awk '(\$1 == \"exonic;splicing\" || \$1 == \"splicing\") { print }' $patient_pp2_variant_function > $patient_splicing"
# 	    vanillaRun "$runme" "$patient_splicing" "result" "Extract splicing variants."
# 	fi	    

#         # gather exonic;splicing from variant_function
# 	patient_all_variant_function=${patient_pass_avlist}.variant_function
# 	patient_all_splicing=${patient_all_variant_function}.splicing
# 	if needsUpdate $patient_all_splicing $patient_all_variant_function
# 	then
# 	    runme="awk '(\$1 == \"exonic;splicing\" || \$1 == \"splicing\") { print }' $patient_all_variant_function > $patient_all_splicing"
# 	    vanillaRun "$runme" "$patient_all_splicing" "result" "Extract passed splicing variants."
# 	fi

# 	patient_not_syn=${patient_pp2_exonic_variant}.not_syn
# 	if needsUpdate $patient_not_syn $patient_pp2_exonic_variant
# 	then
# 	    runme="cut -f2- $patient_pp2_exonic_variant | awk '(\$1 != \"synonymous\") { print }' > $patient_not_syn"
# 	    vanillaRun "$runme" "$patient_not_syn" "result" "Ignore synonymous variants for now."
# 	fi

# 	# hom | 2xhet (remember splicing!)	
# 	patient_recessive_genelist=${patient_pp2_variant_function}.recessive_model
# 	if needsUpdate $patient_recessive_genelist $patient_not_syn $patient_splicing
# 	then
# 	   runme="cat ${patient_not_syn} $patient_splicing |perl -e 'while(<STDIN>) { chomp; @r=split(/\t/); @transcripts=split(/[:;(]+/,\$r[1]); \$gene=\$transcripts[0]; \$gene_row=\$gene.\"\t\".join(\"\t\",@r).\"\n\"; \$nvars{\$gene}++; \$vars{\$gene}.=\$gene_row; if((\$r[7] eq \"hom\") || (\$r[7] eq \"het\" && \$nvars{\$gene} >1) ) { \$modelok{\$gene}=1; } } foreach \$gene (keys %modelok) {print \$vars{\$gene}};' > $patient_recessive_genelist" 
# 	   vanillaRun "$runme" "$patient_recessive_genelist" "result" "Produce recessive model gene list"
# 	fi

# 	patient_indisp=${patient_recessive_genelist}.indisp
# 	if needsUpdate $patient_indisp $patient_recessive_genelist
# 	then
# 	  runme="grep -v -w -f $ANNOVAR_DISPENSABLE $patient_recessive_genelist > $patient_indisp"
# 	  vanillaRun "$runme" "$patient_indisp" "result" "Indispensable genes on recessive gene list."
# 	fi

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
    