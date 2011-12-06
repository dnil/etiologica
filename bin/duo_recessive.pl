#!/usr/bin/perl -w

my $mom =$ARGV[0];
my $kid = $ARGV[1];

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
	if( defined ( $modelok_mom{$chr_start_obs} ) && $modelok_mom{$chr_start_obs} ) {
	    print $gene_row;
	    print $vars_mom{$chr_start_obs};
	}
    }
}
