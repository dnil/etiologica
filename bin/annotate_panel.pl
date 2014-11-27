#!/bin/bash
#
# Daniel Nilsson, 2013-08-08
#
#
# Given a current directory with vcf files of variants and a list of gene symbols, filter variants to return only those that annotate as belonging to the genes on the list.
#
# Usage: annotate_panel.pl [gene_symbol_list_file]
#

MYLIST=$1

if [ -z "$MYLIST" ]
then 
    if [ -e "NMD_and_CMP" ]
    then
	MYLIST="NMD_and_CMP"
	echo "No gene list given! Using $MYLIST."
    elif [ -e "ciliopat_skel" ]
    then
	MYLIST="ciliopat_skel"
	echo "No gene list given! Using $MYLIST."
    else
	echo "No gene list given!"
	exit 1;
    fi
fi

#perl -ne 'chomp; print "\"exonic\",\"",$_,"\",\n";' $MYLIST > ${MYLIST}.exonic_csv

for file in *vcf.gz ; do gunzip -c $file > ${file%%.gz} ; done
for file in *vcf ; do ~/src/annovar/convert2annovar.pl -format vcf4 -allallele $file > ${file%%vcf}avlist ; done
for file in *avlist; do ~/src/annovar/summarize_annovar_custom.pl --buildver hg19 --ver1000g 1000g2012apr --verdbsnp 137NonFlagged --veresp 6500si --alltranscript -remove --outfile $file.var.sum $file ~/src/annovar/humandb/; done

for file in *exome_summary.csv; 
do 
    csv_out="${file%%exome_summary.csv}exome_summary.${MYLIST}.csv" ; 
    grep Func,Gene $file > $csv_out;
    grep -wf ${MYLIST}.exonic_csv $file >> $csv_out; 
done
