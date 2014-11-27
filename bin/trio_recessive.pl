#!/usr/bin/perl -w

# runme="cat ${patient_not_syn_indisp} $patient_pp2_splicing_indisp |perl -e 'while(<STDIN>) { chomp; @r=split(/\t+/); @transcripts=split(/[:;(]+/,\$r[1]); \$gene=\$transcripts[0]; \$gene=~s/\)//; \$gene_row=\$gene.\"\t\".join(\"\t\",@r).\"\n\"; \$nvars{\$gene}++; \$vars{\$gene}.=\$gene_row; if((\$r[7] eq \"hom\") || (\$r[7] eq \"het\" && \$nvars{\$gene} >1) ) { \$modelok{\$gene}=1; } } foreach \$gene (keys %modelok) {print \$vars{\$gene}};' > $patient_recessive_variant_list" 


#cat ${patient_not_syn} $patient_splicing |

my $mom =$ARGV[0];
my $dad =$ARGV[1];
my $kid = $ARGV[2];

open(MOM, $mom);
while (<MOM>) {
    chomp; @r=split(/\t/); 
    
    @transcripts=split(/[:;(]+/,$r[1]); 
    my $gene=$transcripts[0]; 
    $gene=~s/\)//;

    my $gene_row=$gene."\t".join("\t",@r)."\tMOM\n";
    
    my $chr=$r[2];
    my $start=$r[3];
    my $obs=$r[6];
    my $chr_start_obs = join("_",$chr,$start,$obs);

    if($r[7] eq "het") {
	$modelok_mom{$chr_start_obs}=1;
	$vars_mom{$chr_start_obs}=$gene_row;
    }
}
   
open(DAD, $dad);
while (<DAD>) {
    chomp; @r=split(/\t/); 
    
    @transcripts=split(/[:;(]+/,$r[1]); 
    my $gene=$transcripts[0]; 
    $gene=~s/\)//;

    my $gene_row=$gene."\t".join("\t",@r)."\tDAD\n";
    
    my $chr=$r[2];
    my $start=$r[3];
    my $obs=$r[6];
    my $chr_start_obs = join("_",$chr,$start,$obs);

    if($r[7] eq "het") {
	$modelok_dad{$chr_start_obs}=1; 
	$vars_dad{$chr_start_obs}=$gene_row;
    }
}

open(KID, $kid);
while (<KID>) {
    chomp; @r=split(/\t/); 
    
    @transcripts=split(/[:;(]+/,$r[1]); 
    my $gene=$transcripts[0]; 
    $gene=~s/\)//;

    my $gene_row=$gene."\t".join("\t",@r)."\tKID\n";
    
    my $chr=$r[2];
    my $start=$r[3];
    my $obs=$r[6];
    my $chr_start_obs = join("_",$chr,$start,$obs);
    
    if($r[7] eq "hom") {
	if( defined ( $modelok_dad{$chr_start_obs} ) && $modelok_dad{$chr_start_obs} &&
	    defined ( $modelok_mom{$chr_start_obs} ) && $modelok_mom{$chr_start_obs} ) {
	    
	    print $gene_row;
	    print $vars_mom{$chr_start_obs};
	    print $vars_dad{$chr_start_obs};
	}
    }
}
