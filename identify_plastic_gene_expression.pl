#!/usr/bin/perl
#use strict;
use warnings;

my @pop_pair=('EF_EA','FR_EG','SD_SP');
my $index=1; #select for population pair
my $pair= $pop_pair[$index]; my @pop_name=split("_", $pop_pair[$index]);my $pop_c= $pop_name[0];my $pop_w= $pop_name[1];

my $EF=2; my $EA=10; my $FR=18; my $EG=26; my $SD=34; my $SP=42;
my $focal_c=$FR; my $focal_w=$focal_c+8;
my $focal=$focal_c; #select for cold or warm population

my $file = '';
my @rRNA = ('FBgn0013686','FBgn0013688','FBgn0085802',	'FBgn0085753',	'FBgn0267496',	'FBgn0267497',	'FBgn0267498',	'FBgn0267499','FBgn0267500','FBgn0267501','FBgn0267502','FBgn0267503','FBgn0267504','FBgn0267505','FBgn0267506','FBgn0267507','FBgn0267508','FBgn0267509','FBgn0267510','FBgn0267511','FBgn0267512','FBgn0000556');
my %gene_rRNA = map { $_ => 1 } @rRNA;
my $outfile="plas_direction_consistent7Of8_exp_".$pop_w;
open(OUT, ">$outfile");

my %RPM_warm;
for ($f = $focal; $f < $focal+8; $f++){
  #print "$f\n";
$files="warm_RSEM_combined_count.txt";
my $sum_count=0; my $sum_count_length_tpm=0;
  open I, "</mnt/sas0/AD/yhuang349/expression_counts/$files" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
if($count =~ m/FBgn/){
	my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
$sum_count+=$b[$f];
}}
}

 open I, "</mnt/sas0/AD/yhuang349/expression_counts/$files" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
if($count =~ m/FBgn/){
	my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
my $rpm = $b[$f]*1000000/$sum_count;
my $ID=$b[0]."\t".$b[1];
  $RPM_warm{$ID}=$rpm;
}
}}

$files="cold_RSEM_combined_count.txt";
$sum_count=0; $sum_count_length_tpm=0;
  open I, "</mnt/sas0/AD/yhuang349/expression_counts/$files" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
if($count =~ m/FBgn/){
	my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
$sum_count+=$b[$f];
}}
}

 open I, "</mnt/sas0/AD/yhuang349/expression_counts/$files" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
if($count =~ m/FBgn/){
	my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
my $rpm = $b[$f]*1000000/$sum_count;
my $ID=$b[0]."\t".$b[1];
if (exists $RPM_warm{$ID}){
  if(($RPM_warm{$ID}>0)&&($rpm>0)){
  my $dif_ratio=$rpm/$RPM_warm{$ID};

  if (exists $plas_ratio{$ID}){
  $plas_ratio{$ID}=$plas_ratio{$ID}."\t".$dif_ratio;
  }else{
  $plas_ratio{$ID}=$dif_ratio;
  }
}}}
}}}

my $tot=0; my $plas=0;
foreach my $key (keys %plas_ratio) {
my @ratio=split("\t", $plas_ratio{$key});
if (scalar(@ratio)==8){
$tot++;
my $i=0;my $sum_r=0;
foreach my $r (@ratio){
$sum_r+=$r;
if ($r>1){
#if ($r>0){
$i++;
}}
my $avg_r=$sum_r/8;

if(($i<=1)&&($avg_r<1)){
$plas++;
print OUT "$key\t","$i\t","$plas_ratio{$key}\t","$avg_r\n";
}elsif(($i>=7)&&($avg_r>1)){
$plas++;
print OUT "$key\t","$i\t","$plas_ratio{$key}\t","$avg_r\n";
}
}
}
my $prop=$plas/$tot;
print "$plas\t","$tot\t","$prop\n";
