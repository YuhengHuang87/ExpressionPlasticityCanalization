#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
my @t=('Outbred','Inbred');
#need to specify the type of evolutionary changes, choosing one from line 20-23 and one from line 32-35; or not choose any, which means using all genes
my @pop_pair=('EF_EA','FR_EG','SD_SP');
my @index=('2_6','10_14','18_22');
for (my $pair_index = 0; $pair_index <= 2; $pair_index++){
my $pair= $pop_pair[$pair_index];my @pop_name=split("_", $pop_pair[$pair_index]); my $pop_c= $pop_name[0];my $pop_w= $pop_name[1];
my @star_end_index=split("_", $index[$pair_index]);
my $infile;

my %novelist;my %nonlist;
$infile=$pair."_QST_quantile_top5";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
	#if (($b[4] eq "concordant")||($b[4] eq "reversing")){
	#if ($b[4] eq "neutralizing"){
		#if ($b[5] eq "undercompensate"){
			#if ($b[5] eq "overcompensate"){
$novelist{$b[0]}=$count;
#}
}
$infile=$pair."_QST_quantile_nontop5";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
	#if (($b[4] eq "concordant")||($b[4] eq "reversing")){
	#if ($b[4] eq "neutralizing"){
		#if ($b[5] eq "undercompensate"){
			#if ($b[5] eq "overcompensate"){
		$nonlist{$b[0]}=$count;
#}
}

my $size1=4; my $size2=4;my $EF=2; my $EA=6; my $FR=10; my $EG=14; my $SD=18; my $SP=22;
my $OUTStart1=$star_end_index[0];my $OUTEnd1=$OUTStart1+$size1-1;
my $OUTStart2=$star_end_index[1];my $OUTEnd2=$OUTStart2+$size2-1;

my $INStart1=$star_end_index[0]+24;my $INEnd1=$INStart1+$size1-1;
my $INStart2=$star_end_index[1]+24;my $INEnd2=$INStart2+$size2-1;

my $V_cold_outlier=0; my $V_warm_outlier=0;my $V_cold_nonoutlier=0;;my $V_warm_nonoutlier=0;
my $k=0; my $InOut_cold=0; my $InOut_warm=0; my $total=0;

my $cold_ratio;my $warm_ratio; my @cold_var_OUT; my @warm_var_OUT;my @cold_var_IN; my @warm_var_IN;
$infile="OUT_IN_RPM_rRNA_removed";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my $id = $b[0];
	my @pop1_OUT = @b[$OUTStart1..$OUTEnd1];my @pop2_OUT = @b[$OUTStart2..$OUTEnd2];
	my @pop1_IN = @b[$INStart1..$INEnd1];my @pop2_IN = @b[$INStart2..$INEnd2];

	my $z1_OUT = average(@pop1_OUT); my $z2_OUT = average(@pop2_OUT);
	my $z1_IN = average(@pop1_IN); my $z2_IN = average(@pop2_IN);
my $sum_var1_OUT = 0;
for ($k = 0; $k < @pop1_OUT; $k++){
$sum_var1_OUT+=($pop1_OUT[$k]-$z1_OUT)**2;
}
my $sum_var1_IN = 0;
for ($k = 0; $k < @pop1_IN; $k++){
$sum_var1_IN+=($pop1_IN[$k]-$z1_IN)**2;
}
my $sum_var2_OUT = 0;
for ($k = 0; $k < @pop2_OUT; $k++){
$sum_var2_OUT+=($pop2_OUT[$k]-$z2_OUT)**2;
}
my $sum_var2_IN = 0;
for ($k = 0; $k < @pop2_IN; $k++){
$sum_var2_IN+=($pop2_IN[$k]-$z2_IN)**2;
}
my $var1_OUT=$sum_var1_OUT/$size1;my $var1_IN=$sum_var1_IN/$size1;
my $var2_OUT=$sum_var2_OUT/$size2;my $var2_IN=$sum_var2_IN/$size2;

if(($var1_OUT>0)&&($var2_OUT>0)&&($var1_IN>0)&&($var2_IN>0)){

	if (exists $novelist{$id}){
		$total++;
		push @cold_var_OUT,$var1_OUT ; push @warm_var_OUT, $var2_OUT;push @cold_var_IN, $var1_IN; push @warm_var_IN, $var2_IN;
	my $cold_InTot_ratio=$var1_IN/$var1_OUT;
	my $warm_InTot_ratio=$var2_IN/$var2_OUT;
	if ($cold_InTot_ratio>1){
	$InOut_cold++;
}
if ($warm_InTot_ratio>1){
$InOut_warm++;
}
	if ($cold_InTot_ratio>$warm_InTot_ratio){
	$V_cold_outlier++;
	}else{
	$V_warm_outlier++;
	}
	}
	if (exists $nonlist{$id}){
		$total++;
		push @cold_var_OUT,$var1_OUT ; push @warm_var_OUT, $var2_OUT;push @cold_var_IN, $var1_IN; push @warm_var_IN, $var2_IN;
		my $cold_InTot_ratio=$var1_IN/$var1_OUT;
		my $warm_InTot_ratio=$var2_IN/$var2_OUT;
		if ($cold_InTot_ratio>1){
		$InOut_cold++;
	}
	if ($warm_InTot_ratio>1){
	$InOut_warm++;
	}
	if ($cold_InTot_ratio>$warm_InTot_ratio){
	$V_cold_nonoutlier++;
	}else{
	$V_warm_nonoutlier++;
	}
	}
}
}
my $OutIn_cold= $total-$InOut_cold;my $OutIn_warm= $total-$InOut_warm;
my $mean_Vout_cold=median(@cold_var_OUT);my $mean_Vout_warm=median(@warm_var_OUT);
my $mean_Vin_cold=median(@cold_var_IN);my $mean_Vin_warm=median(@warm_var_IN);

print "$pair\t","$total\t","$InOut_cold\t","$OutIn_cold\t","$InOut_warm\t","$OutIn_warm\t","$V_cold_outlier\t","$V_warm_outlier\t","$V_cold_nonoutlier\t","$V_warm_nonoutlier\n";
print "$pair\t","$mean_Vin_warm\t","$mean_Vin_cold\t","$mean_Vout_warm\t","$mean_Vout_cold\n";
}



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
