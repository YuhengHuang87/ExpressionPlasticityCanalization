#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);

my @pop_pair=('EF_EA','FR_EG','SD_SP');
my $index=2; #specify population pair
my $pair= $pop_pair[$index];my @pop_name=split("_", $pop_pair[$index]);my $pop_c= $pop_name[0];my $pop_w= $pop_name[1];
my $EF=2; my $EA=10; my $FR=18; my $EG=26; my $SD=34; my $SP=42;
my $focal_c=$SD; my $focal_w=$focal_c+8; #indicate the cold and warm population columns
my $plas_pop=$pop_w; #specify warm or cold population for showing plasticity
my $condition_warm_c=$focal_c+50; my $condition_warm_w=$focal_w+50;

my @Cold=($focal_c..$focal_c+7); my @Warm=($focal_w..$focal_w+7);

my $j=0; my $k=0;
my @r; my $line=0;
my $infile="cold_warm_comb_RPM_rRNA_removed_median_count_200";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my $id = $b[0];
	my @p1 = ($b[$focal_c],$b[$focal_c+1],$b[$focal_c+2],$b[$focal_c+3],$b[$focal_c+4],$b[$focal_c+5],$b[$focal_c+6],$b[$focal_c+7]);
	my @p2 = ($b[$focal_w],$b[$focal_w+1],$b[$focal_w+2],$b[$focal_w+3],$b[$focal_w+4],$b[$focal_w+5],$b[$focal_w+6],$b[$focal_w+7]);

my @sorted_p1 = sort { $a <=> $b } @p1;
my @sorted_p2 = sort { $a <=> $b } @p2;

	my @pop1 = @sorted_p1[1..6];#exclude the extremely high and low values
    my @pop2 = @sorted_p2[1..6];
		my $z1 = average(@pop1);
		my $z2 = average(@pop2);
		my $cold_warm_ratio=$z1/$z2;
		push @r, $cold_warm_ratio;
		$line++;
}
my $ratio = median(@r);

my %qst; my $num_gene=0; my %type;my %type_neutra;
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my $id = $b[0];
    my @p1 = ($b[$focal_c],$b[$focal_c+1],$b[$focal_c+2],$b[$focal_c+3],$b[$focal_c+4],$b[$focal_c+5],$b[$focal_c+6],$b[$focal_c+7]);
    my @p2 = ($b[$focal_w],$b[$focal_w+1],$b[$focal_w+2],$b[$focal_w+3],$b[$focal_w+4],$b[$focal_w+5],$b[$focal_w+6],$b[$focal_w+7]);
    my @p3 = ($b[$focal_w+50],$b[$focal_w+1+50],$b[$focal_w+2+50],$b[$focal_w+3+50],$b[$focal_w+4+50],$b[$focal_w+5+50],$b[$focal_w+6+50],$b[$focal_w+7+50]);
    my $Exp_C15=median(@p1);my $Exp_W15=median(@p2);my $Exp_W25=median(@p3);
    my $plas_dif=$Exp_W15-$Exp_W25;
    my $evol_dif=$Exp_C15-$Exp_W15;

    if ($evol_dif*$plas_dif>0){
    $type{$id}="concordant";
		$type_neutra{$id}="NA";
    }else{
			if (abs($evol_dif)>2*abs($plas_dif)){
			$type{$id}="reversing";
			$type_neutra{$id}="NA";
		}else{
			$type{$id}="neutralizing";
		if(abs($evol_dif)<(abs($plas_dif))){
			$type_neutra{$id}="undercompensate";
		}else{
			$type_neutra{$id}="overcompensate";
		}
  }}

my @sorted_p1 = sort { $a <=> $b } @p1;
my @sorted_p2 = sort { $a <=> $b } @p2;

	my @pop1 = @sorted_p1[1..6];
    my @pop2 = @sorted_p2[1..6];

foreach my $x (@pop2) { $x = $x * $ratio; }
	my @pop_z = (@pop1,@pop2);

#mean phenotype for each pop and the grant mean
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
my $cold_warm_ratio=$z1/$z2;
my $evolve_prop=$z1/($z1+$z2);
$num_gene+=1;
$qst{$b[0]}=$Qst;
}}

my %plas;
$infile="plas_direction_consistent7Of8_exp_".$plas_pop;
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my $id = $b[0];
	$plas{$id} = $b[10];
}

#my $outfile1=$pair."_QST_quantile_top5_type_plas_".$plas_pop;
my $outfile1=$pair."_QST_quantile_top5";
open(OUT1, ">$outfile1");
my $outfile2=$pair."_QST_quantile_nontop5";
#my $outfile2=$pair."_QST_quantile_nontop5_type_plas_".$plas_pop;
open(OUT2, ">$outfile2");

my $top=int($num_gene*0.05);
my $rank=0; my $check_o=0;my $total_o=0;my $check_no=0;my $total_no=0;my $outlier_num=0;my $nonoutlier_num=0;
foreach my $key (sort { $qst{$b} <=> $qst{$a} } keys %qst) {
my $quantile=$rank/$num_gene;
if($rank<=$top){
	#if (exists $plas{$key}){
		$total_o++;
		if ($type{$key} eq "concordant"){
			$check_o++;
		}
            print OUT1 "$key\t","$qst{$key}\t","$quantile\t","Outlier\t","$type{$key}\t","$type_neutra{$key}\n";

#}

}else{
	$nonoutlier_num++;
	#if (exists $plas{$key}){
		$total_no++;
		if ($type{$key} eq "concordant"){
			$check_no++;
		}
print OUT2 "$key\t","$qst{$key}\t","$quantile\t","Nonout\t","$type{$key}\t","$type_neutra{$key}\n";
#}
}
$rank+=1;
}
print "$plas_pop\t","$top\t","$nonoutlier_num\t","$check_o\t","$total_o\t","$check_no\t","$total_no\n";



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
