#!/bin/bash

# (c)2011 Daniel Nilsson, daniel.nilsson@scilifelab.se, daniel.nilsson@karolinska.se, daniel.nilsson@ki.se
# Copyright 2011, held by the author.
# Released under the Perl Artistic License.
#
# Documentation available in this file. Use e.g. perldoc to view, or pod2man to create a manual.
#
# library QC:
#

: <<'POD_INIT'

=head1 NAME

capture_qc.sh

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@scilifelab.se, daniel.nilsson@ki.se, daniel.nilsson@karolinska.se

=head1 SYNOPSIS

cd project
./capture_qc.sh

=head1 DESCRIPTION

Process patient secondary_analysis dir to generate a QC report.

Target capture enrichment status can be reported if desired. 
Please set the CAPTURE_TARGET_BED environment variable to a space separated list of bed target range files that you want to investigate.

=head1 ENVIRONMENT VARIABLES

The pipeline adheres to environment exported variables.

=over 4

=item C<BINDIR> I<path>

The path to the etiologia pipeline bin dir. Defaults to the (absolute) dir of the main script.

POD_INIT

if [ -z "$BINDIR" ]
then 
    BINDIR=`dirname $my_path`
fi

PIPELINE=etiologica_qc
PIPELINEFUNK=$BINDIR/pipelinefunk.sh

. $PIPELINEFUNK
log "Running $PIPELINE." "main"
log "BINDIR: $BINDIR" "main"

if [ -z "$REFERENCE" ]
then 
    REFERENCE=human_g1k_v37.fasta.gz
fi
log "REFERENCE: $REFERENCE" "main"

if [ -z "$GATK_REFERENCE" ]
then
    ~/glob/private/GATK/hg19/human_g1k_v37.fasta
fi 
log "GATK_REFERENCE - a GATK style reference: $GATK_REFERENCE" "main"

if [ -z "$GATK_JAR" ]
then    
    GATK_JAR=/bubo/sw/apps/bioinfo/GATK/1.4.5/GenomeAnalysisTK.jar
fi
log "GATK_JAR - path to GATK jar file: $GATK_JAR" "main"

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

if [ -z "$BEDTOOLS_BIN" ]
then 
    export BEDTOOLS_BIN=~/src/BEDTools/bin
fi
log "BEDTOOLS_BIN: $TMP"

cd $BINDIR; 
pipeline_git_status=`git show --pretty=format:"%h %ci %s"; git status --porcelain`
cd -

log "pipeline git status: $pipeline_git_status" "main"

CAPTURE_TARGET_BED="~/cocaclingen/private/SureSelect_All_Exon_50mb_with_annotation_hg19_bed_num.bed ~/cocaclingen/private/cocaclingen_transcripts.exons.nochr.bed SureSelect_All_Exon_50mb_with_annotation_hg19_bed_num.cocapanel.bed"

if [ -z "$CAPTURE_TARGET_BED" ] 
then
    log "No capture target given. Assuming WGS."
fi

for dir in patient metadata/capture_qc
do
    if [ ! -d "$dir" ]
	mkdir -p $dir
    fi
done 

for fastq_file in *_1.fastq.gz ; 
do 
    fastqc_data=${fastq_file%%.fastq.gz}_fastqc/fastqc_data.txt
    input_readpairs=`grep Total\ Sequences ${fastqc_data}`
    echo $fastq_file has $input_readpairs readpairs.
do

# reference genome BEDtools genome file
for target_bed_file in $CAPTURE_TARGET_BED
do
# target BED file
    target_bed_file_summary=${target_bed_file}.summary

    if needsUpdate $target_bed_file_summary $target_bed_file

	runme="date > $target_bed_file_summary; echo $target_bed_file >> $target_bed_file_summary; awk 'BEGIN {sum=0; count=0;} {sum=sum+\$3-\$2; count=count+1} END {print sum,count,sum/count;}' $target_bed_file >> $target_bed_file_summary"
	# plot? come on, this is a one-off-ish thing.

	vanillaRun $runme
    fi
done


for target_bed_file in $CAPTURE_TARGET_BED
do 

    for alignment in *mosaik.bam ; 
    do 
    # work lock!
	if workLockOk $alignment
	then
	    log "Locked $alignment for capture_qc." "main"

	    alignment_target_summary=${alignment%%.mosaik.bam}.${CAPTURE_TARGET_BED%%.bed}.CLW.summary
	    alignment_target_bed=${alignment%%.mosaik.bam}.${CAPTURE_TARGET_BED%%.bed}.CLW.bed

	    if needsUpdate ${alignment_target_summary} ${alignment} or needsUpdate ${alignment_target_bed} ${alignment}
	    then
		runme="java -jar $GATK_JAR  -T CallableLoci -I $alignment -summary ${alignment_target_summary} -o ${alignment_target_bed} -R $GATK_REFERENCE" 
		vanillaRun $runme
	    fi
	fi
    done
done


java -jar /bubo/sw/apps/bioinfo/GATK/1.4.5/GenomeAnalysisTK.jar -T CallableLoci -I Sample_1_11_1.mosaik.bam -summary Sample_1_11_1.CLW.summary -o Sample_1_11_1.CLW.bed -R ~/glob/private/GATK/hg19/human_g1k_v37.fasta 
sed 's/ /\t/g' Sample_1_11_1.CLW.bed > Sample_1_11_1.CLW.tab.bed
 ~/src/BEDTools-Version-2.13.1/bin/intersectBed -a Sample_1_11_1.CLW.tab.bed -b ../../SureSelect_All_Exon_50mb_with_annotation_hg19_bed_num.bed > Sample_1_11_1.CLW.sureselect_ae_50.intersect.bed
 awk 'BEGIN {lc=0; noc=0; c=0} ($4=="CALLABLE"){c=c+$3-$2} ($4=="LOW_COVERAGE") {lc=lc+$3-$2} ($4=="NO_COVERAGE") {noc=noc+$3-$2} END { print "callable: ",c, "low coverage: ",lc," no_coverage: ",noc; }' Sample_1_11_1.CLW.sureselect_ae_50.intersect.bed



# input bam


# input bed regions
head -2 
tail -n +3 OID33040.hg19.3.bed > OID33040.hg19.3.actual.bed

# total reads
#zcat fastq.gz |grep -c "^@"
grep -c "^@" 1_110617_Ac00hfabxx_CO666_index5_1.fastq

# total reads mapped

bamToBed -i 3_110617_Ac00hfabxx_4-sort-dup-gatkrecal-realign.bam > 3_110617_Ac00hfabxx_4-sort-dup-gatkrecal-realign.bed
wc -l 3_110617_Ac00hfabxx_1-sort-dup-gatkrecal-realign.bed

# reads mapped to region?
outer_limits_region=../OID33040.hg19.3.outer.bed
intersectBed -abam 1_110617_Ac00hfabxx_1-sort-dup-gatkrecal-realign.bam -b ../OID33040.hg19.3.outer.bed > 1_110617_Ac00hfabxx_1-sort-dup-gatkrecal-realign.OID33040.outer.bed
wc -l 1_110617_Ac00hfabxx_1-sort-dup-gatkrecal-realign.OID33040.outer.bed

# >=400x
genomeCoverageBed -i 1_110617_Ac00hfabxx_4-sort-dup-gatkrecal-realign.outer.bed -g ../hg19.genome -max 200 > 1_110617_Ac00hfabxx_4-sort-dup-gatkrecal-realign.outer.max200.genomecov
# >=6x
genomeCoverageBed -i 1_110617_Ac00hfabxx_4-sort-dup-gatkrecal-realign.outer.bed -g ../hg19.genome -max 6 > 1_110617_Ac00hfabxx_4-sort-dup-gatkrecal-realign.outer.max6.genomecov
# >=1x? (also gives #reads in region)
coverageBed -a 1_110617_Ac00hfabxx_1-sort-dup-gatkrecal-realign.OID33040.outer.bed -b ../OID33040.hg19.3.outer.bed


# fraction target bases covered >=1x
awk 'BEGIN {sum=0; count=0; } {sum=sum+$10; count=count+1;} END {print sum,count,sum/count;}' 6_110825_AD035WACXX_927_ind_12_1.sureselect.bedcov

# fraction target bases covered >=6x
intersectBed  --abam

#bed out flag available to modern intersectBed --abam
bamToBed

genomeCoverageBed -i 6_110825_AD035WACXX_927_ind_12_1.sureselect.bed -g human.hg19.GLdot.nochr.genome -max 6 > 6_110825_AD035WACXX_927_ind_12_1.mosaik.hg19.sureselect.genomecov.max6
for file in *max6 ; do echo $file ; awk '($2 ==6) {print }' $file |grep genome ;done

# total number of bases callable

# total number of variants called

# confidently

# breakdown on indel/ss/ns syn UTR intergenic intronic

# pp2 FDR<0.1



