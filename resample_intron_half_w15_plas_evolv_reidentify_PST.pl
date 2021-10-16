#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
use List::Util 'shuffle';

my @pop_pair=('EF_EA','FR_EG','SD_SP');
my $index=2; #specify population pair
my @pop_name=split("_", $pop_pair[$index]); my $pop_c= $pop_name[0]; my $pop_w= $pop_name[1];
my $pair= $pop_pair[$index];
my $EA=1; my $EF=9; my $EG=17; my $FR=25; my $SD=33; my $SP=41;
my $focal_c=$SD; my $focal_w=$SP; #specify the cold and warm populations

my @Cold_15=($focal_c..$focal_c+7); my @Warm_15=($focal_w..$focal_w+7);
my %plas; my $j=0;
my $infile="Plas_intron_cold_warm_direction7of8_".$pop_w;
print "$infile\n";
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my $id = $b[1];
$j+=1;
	$plas{$id} = $b[10];
}

my $i=0;
my $pass_top=0;
my %permute;
	while ($i<70){
my @Warm_15_shuffled = shuffle(@Warm_15);
my @plas_warm_15=(@Warm_15_shuffled[0..3]);
my @evol_warm_15=(@Warm_15_shuffled[4..7]);
    my @sorted_plas_warm_15 = sort { $a <=> $b } @plas_warm_15;
    my @sorted_evol_warm_15 = sort { $a <=> $b } @evol_warm_15;

my @permute_set=(@Cold_15[0..3],@sorted_plas_warm_15,@sorted_evol_warm_15);

my $combine = join( ' ', @permute_set );
if (exists $permute{$combine}){
}else{
  $permute{$combine}=$i;
  $i++;

	my $k=0;
	my %qst; my $num_gene=0; my %exp_median;
my $infile="/mnt/sas0/AD/yhuang349/intron_analysis/cold_warm_comb_count_total.txt";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my @clu=split(":", $b[0]);
		my $id = $clu[0].":".$clu[1].":".$clu[2];
		if ($b[$focal_c+49]*$b[$focal_c+1+49]*$b[$focal_c+2+49]*$b[$focal_c+3+49]*$b[$focal_c+4+49]*$b[$focal_c+5+49]*$b[$focal_c+6+49]*$b[$focal_c+7+49]*$b[$permute_set[4]+49]*$b[$permute_set[5]+49]*$b[$permute_set[6]+49]*$b[$permute_set[7]+49]*$b[$permute_set[8]+49]*$b[$permute_set[9]+49]*$b[$permute_set[10]+49]*$b[$permute_set[11]+49]*$b[$permute_set[4]+148]*$b[$permute_set[5]+148]*$b[$permute_set[6]+148]*$b[$permute_set[7]+148] > 0){
	my @C15 = ($b[$focal_c]/$b[$focal_c+49],$b[$focal_c+1]/$b[$focal_c+1+49],$b[$focal_c+2]/$b[$focal_c+2+49],$b[$focal_c+3]/$b[$focal_c+3+49],$b[$focal_c+4]/$b[$focal_c+4+49],$b[$focal_c+5]/$b[$focal_c+5+49],$b[$focal_c+6]/$b[$focal_c+6+49],$b[$focal_c+7]/$b[$focal_c+7+49]);
	my @P_W15 = ($b[$permute_set[4]]/$b[$permute_set[4]+49],$b[$permute_set[5]]/$b[$permute_set[5]+49],$b[$permute_set[6]]/$b[$permute_set[6]+49],$b[$permute_set[7]]/$b[$permute_set[7]+49]);
  my @E_W15 = ($b[$permute_set[8]]/$b[$permute_set[8]+49],$b[$permute_set[9]]/$b[$permute_set[9]+49],$b[$permute_set[10]]/$b[$permute_set[10]+49],$b[$permute_set[11]]/$b[$permute_set[11]+49]);
  my @P_W25 = ($b[$permute_set[4]+99]/$b[$permute_set[4]+148],$b[$permute_set[5]+99]/$b[$permute_set[5]+148],$b[$permute_set[6]+99]/$b[$permute_set[6]+148],$b[$permute_set[7]+99]/$b[$permute_set[7]+148]);

my $exp_P_W15=median(@P_W15);my $exp_P_W25=median(@P_W25);
my $exp_C15=median(@C15);my $exp_E_W15=median(@E_W15);

my @sorted_p1 = sort { $a <=> $b } @C15;my @sorted_p2 = sort { $a <=> $b } @E_W15;
my @pop1 = @sorted_p1[1..6];
my @pop2 = @sorted_p2;
my @pop_z = (@pop1,@pop2);
my $z1 = average(@pop1); #here is just the average transcript level for population1.
my $z2 = average(@pop2);
my $z = average(@pop_z);
my $size1=scalar(@pop1);
my $size2=scalar(@pop2);
# var within and between and Qst
my $sum_var=0;
my $sum_var1 = 0;
for ($k = 0; $k < @pop1; $k++){
$sum_var1+=($pop1[$k]-$z1)**2;
$sum_var+=($pop1[$k]-$z)**2;
}
my $sum_var2 = 0;
for ($k = 0; $k < @pop2; $k++){
$sum_var2+=($pop2[$k]-$z2)**2;
$sum_var+=($pop2[$k]-$z)**2;
}
my $var_within=($sum_var1/$size1+$sum_var2/$size2)/2;
my $var_within1=$sum_var1/$size1;
my $var_within2=$sum_var2/$size2;

my $var_betw = (($z1-$z)**2+($z2-$z)**2)/(2-1);#from Lande 1992. Evolution
my $var_total=$var_betw+2*$var_within;

my $var_total2=$sum_var/($size1+$size2);

if($var_total>0){
my $Qst = $var_betw/$var_total; #same as Qst_Wright: ($var_total2-$var_within)/$var_total2;

$num_gene+=1;
$qst{$b[0]}=$Qst;
$exp_median{$b[0]}=$exp_C15."\t".$exp_E_W15."\t".$exp_P_W15."\t".$exp_P_W25;
}}}

my $outfile1="permute_cold_warm_intron_".$pair."/".$pair."_intron_QST_quantile_top5_resample_W15_".$i;
open(OUT1, ">$outfile1");
my $outfile2="permute_cold_warm_intron_".$pair."/".$pair."_intron_QST_quantile_nontop5_resample_W15_".$i;
open(OUT2, ">$outfile2");

my $top=int($num_gene*0.05);
my $rank=0; my %cluster; my $outlier_num=0;my %cluster_bg;
	foreach my $key (sort { $qst{$b} <=> $qst{$a} } keys %qst) {
my $quantile=$rank/$num_gene;
my @exon = split(":", $key); my $name=$exon[0].":".$exon[1].":".$exon[2];
if($rank<=$top){
if(exists $cluster{$exon[3]}){
}else{
	if (exists $plas{$name}){
  print OUT1 "$key\t","$qst{$key}\t","$quantile\t","$exp_median{$key}\t","plas\n";
}else{
	print OUT1 "$key\t","$qst{$key}\t","$quantile\t","$exp_median{$key}\t","nonplas\n";
}
$cluster{$exon[3]}=$qst{$key};
$outlier_num++;
}

}else{
if((exists $cluster_bg{$exon[3]})||(exists $cluster{$exon[3]})){

}else{
	if (exists $plas{$name}){
  print OUT2 "$key\t","$qst{$key}\t","$quantile\t","$exp_median{$key}\t","plas\n";
}else{
	print OUT2 "$key\t","$qst{$key}\t","$quantile\t","$exp_median{$key}\t","nonplas\n";
}
$cluster_bg{$exon[3]}=$qst{$key};
}
        }
$rank+=1;
}
}}

my $outfile=$pair."_resample_W15_concordant_intron";
open(OUT, ">$outfile");

my $extrem=0;my $tot_pass_plas=0;
	for (my $i=1; $i <= 70; $i++) {
my $plas_num_top=0;my $cons_top=0; my $anti_top=0;
my $file="permute_cold_warm_intron_".$pair."/".$pair."_intron_QST_quantile_top5_resample_W15_".$i;
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my $id = $b[0];
if ($b[7] eq "plas"){
	$plas_num_top++;
	my $exp_c_15=$b[3];
	my $E_w_15=$b[4];
	my $P_w_15=$b[5];
		my $exp_w_25=$b[6];

my $plas_dif=$P_w_15-$exp_w_25;
my $evol_dif=$exp_c_15-$E_w_15;

if ($evol_dif*$plas_dif>0){
$cons_top++;
}else{
$anti_top++;
}
}}
close FILE1;

my $plas_num_nontop=0;my $cons_nontop=0; my $anti_nontop=0;

$file="permute_cold_warm_intron_".$pair."/".$pair."_intron_QST_quantile_nontop5_resample_W15_".$i;
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my $id = $b[0];
if ($b[7] eq "plas"){
	$plas_num_nontop++;
	my $exp_c_15=$b[3];
	my $E_w_15=$b[4];
	my $P_w_15=$b[5];
		my $exp_w_25=$b[6];

	my $plas_dif=$P_w_15-$exp_w_25;
	my $evol_dif=$exp_c_15-$E_w_15;

if ($evol_dif*$plas_dif>0){
$cons_nontop++;
}else{
$anti_nontop++;
}
}
}
close FILE1;

	$tot_pass_plas++;
my $ratio_top=$cons_top/($cons_top+$anti_top);
my $ratio_nontop=$cons_nontop/($cons_nontop+$anti_nontop);
my $dif_r_top_nontop=$ratio_top-$ratio_nontop;
print OUT "$i\t","$plas_num_top\t","$ratio_top\t","$plas_num_nontop\t","$ratio_nontop\t","$dif_r_top_nontop\n";
if ($dif_r_top_nontop >0){
	$extrem++;
}
}

print "$pair\t","$extrem\t","$tot_pass_plas\n";



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
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
