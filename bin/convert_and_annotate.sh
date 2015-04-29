#!/bin/bash
#
# Daniel Nilsson, 2013-08-08
#
# Given a current directory with vcf files of variants, use annovar to filter variants
#
# Usage: convert_and_annotate.sh
#

for file in *vcf.gz ; do gunzip -c $file > ${file%%.gz} ; done
for file in *vcf ; do ~/src/annovar/convert2annovar.pl -format vcf4 $file > ${file%%vcf}avlist ; done

for file in *avlist; do ~/src/annovar/table_annovar.pl --buildver hg19 $file ~/src/annovar/humandb/ -vcfdbfile $LOCAL_CLIN_DB -protocol refGene,phastConsElements46way,genomicSuperDups,popfreq_max,exac02,esp6500siv2_all,1000g2014oct_all,vcf,snp138,ljb2_all,clinvar_20140929 -operation g,r,r,f,f,f,f,f,f,f,f -outfile $file.sum -remove -otherinfo -csvout; done

