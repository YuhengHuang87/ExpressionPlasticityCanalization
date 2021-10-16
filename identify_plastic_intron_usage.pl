#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);

my @pop_pair=('EF_EA','FR_EG','SD_SP');
my $index=2; #select for the population pair
my $pair= $pop_pair[$index]; my @pop_name=split("_", $pop_pair[$index]);my $pop_c= $pop_name[0];my $pop_w= $pop_name[1];

my $outfile1="Plas_intron_cold_warm_direction7of8_".$pop_c;
open(OUT1, ">$outfile1");
my $outfile="All_intron_cold_warm_direction_".$pop_c;
open(OUT, ">$outfile");
my $h=0; my $j=0; my $k=0;
my $size1=8; my $size2=8;my $EA=1; my $EF=9; my $EG=17; my $FR=25; my $SD=33; my $SP=41; my $rep=49;

my $SampleStart1=$SD;#select for the focal population
my $SampleEnd1=$SampleStart1+$size1-1;

my %total_cold;
open(FILE1,'<', 'total_Intron_count_gene_reads_cold')||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
	$total_cold{$b[0]}=$count;
	}
close FILE1;
my %intron_perce_cold;
open(FILE1,'<', 'intron_count_cold.txt')||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @line = split("\r", $count);
for my $line (@line){
my @a=split(" ", $line);
		my $id = $a[0];
if($id ne "location"){
	my @clu = split(":", $id);
my $intron=$clu[0].":".$clu[1].":".$clu[2];
#if (exists $hash1{$intron}){
my $clu_name=$clu[3];
if (exists $total_cold{$clu_name}){
my @sum_type = split("\t", $total_cold{$clu_name});
for (my $i = $SampleStart1; $i <= $SampleEnd1; $i++){
if ($sum_type[$i]>0){
my $perct=$a[$i]/$sum_type[$i];
if (exists $intron_perce_cold{$intron}){
$intron_perce_cold{$intron}=$intron_perce_cold{$intron}."\t".$perct;
}else{
$intron_perce_cold{$intron}=$perct;
}}}
}}}}
close FILE1;

my %total_warm;my %intron_perce_warm;
open(FILE1,'<', 'total_Intron_count_gene_reads_warm')||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
	$total_warm{$b[0]}=$count;
	}
close FILE1;
open(FILE1,'<', 'intron_count_warm.txt')||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @line = split("\r", $count);
for my $line (@line){
my @a=split(" ", $line);
		my $id = $a[0];
if($id ne "location"){
	my @clu = split(":", $id);
my $intron=$clu[0].":".$clu[1].":".$clu[2];
my $clu_name=$clu[3];
if (exists $total_warm{$clu_name}){
my @sum_type = split("\t", $total_warm{$clu_name});
for (my $i = $SampleStart1; $i <= $SampleEnd1; $i++){
if ($sum_type[$i]>0){
my $perct=$a[$i]/$sum_type[$i];
if (exists $intron_perce_warm{$intron}){
$intron_perce_warm{$intron}=$intron_perce_warm{$intron}."\t".$perct;
}else{
$intron_perce_warm{$intron}=$perct;
}}}
}}}}
close FILE1;

my $plas=0;my $tot=0;
foreach my $key (keys %intron_perce_cold) {
my @intron_freq_cold=split("\t",$intron_perce_cold{$key});
my $len_cold=scalar @intron_freq_cold;
if (exists $intron_perce_warm{$key}){
my @intron_freq_warm=split("\t",$intron_perce_warm{$key});
my $len_warm=scalar @intron_freq_warm;

if (($len_cold == 8)&&($len_warm == 8)){
	$tot++;
my $i=0; my $alldif=$key; my $count_dir=0; my $sum_dif=0;

for ($i = 0; $i < $size1; $i++){
	my $dif = $intron_freq_cold[$i]-$intron_freq_warm[$i];
	$sum_dif+=$dif;
	$alldif=$alldif."\t".$dif;
	if ($dif > 0){
$count_dir++;
}}
my $avg_dif=($sum_dif)/8;
print OUT "$count_dir\t","$alldif\t","$avg_dif\n";

if(($count_dir<=1)&&($avg_dif<0)){
$plas++;
print OUT1 "$count_dir\t","$alldif\t","$avg_dif\n";
}elsif(($count_dir>=7)&&($avg_dif>0)){
$plas++;
print OUT1 "$count_dir\t","$alldif\t","$avg_dif\n";
}
}}
}

my $prop=$plas/$tot;
print "$plas\t","$tot\t","$prop\n";

sub average {
my @array = @_; # save the array passed to this function
my $sum; # create a variable to hold the sum of the array's values
foreach (@array) { $sum += $_; } # add each element of the array
# to the sum
return $sum/@array; # divide sum by the number of elements in the
# array to find the mean
}

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
