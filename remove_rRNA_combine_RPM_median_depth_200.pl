#!/usr/bin/perl
#use strict;
use warnings;

my $file = '';
my @rRNA = ('FBgn0013686','FBgn0013688','FBgn0085802',	'FBgn0085753',	'FBgn0267496',	'FBgn0267497',	'FBgn0267498',	'FBgn0267499','FBgn0267500','FBgn0267501','FBgn0267502','FBgn0267503','FBgn0267504','FBgn0267505','FBgn0267506','FBgn0267507','FBgn0267508','FBgn0267509','FBgn0267510','FBgn0267511','FBgn0267512','FBgn0000556');
my %gene_rRNA = map { $_ => 1 } @rRNA;

my $outfile="cold_RPM_rRNA_removed_median_count_200";
open(OUT, ">$outfile");
my $outfile1="cold_median_count_all";
open(OUT1, ">$outfile1");

my %RPM;
for ($f = 2; $f < 50; $f++){
  #print "$f\n";
$files="expression_count_cold_condition";
my $sum_count=0; my $sum_count_length_tpm=0;
open(I,"<$files")||die"$!";
	while(my $count = <I>){
  	chomp($count);
if($count =~ m/FBgn/){
	my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
$sum_count+=$b[$f];
}}
}
#print "$sum_count\n";

open(I,"<$files")||die"$!";
	while(my $count = <I>){
  	chomp($count);
if($count =~ m/FBgn/){
	my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
$rpm = $b[$f]*1000000/$sum_count;
my $ID=$b[0]."\t".$b[1];
if (exists $RPM{$ID}){
  $RPM{$ID}=$RPM{$ID}."\t".$rpm;
}else{
  $RPM{$ID}=$rpm;
}
}}}
}

my $line=0; my $screen=0;
open(I,"<$files")||die"$!";
 while(my $count = <I>){
   chomp($count);
if($count =~ m/FBgn/){
 my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
  my $ID=$b[0]."\t".$b[1];
  if (exists $RPM{$ID}){
    @a=@b[2..49];
#print "@a\n";
if (@a == 48){
  $line++;
  my $median_count=median(@a);
  my $avg_count=average(@a);
  print OUT1 "$ID\t","$median_count\n";
  if ($median_count>200){
    $screen++;
    print OUT "$ID\t","$RPM{$ID}\n";
  }
}}}}}
print "$line\t","$screen\n";

$outfile="warm_RPM_rRNA_removed";
open(OUT, ">$outfile");
$outfile1="warm_median_count_all";
open(OUT1, ">$outfile1");

my %RPM_warm;
for ($f = 2; $f < 50; $f++){
  #print "$f\n";
$files="expression_count_warm_condition";
my $sum_count=0; my $sum_count_length_tpm=0;
open(I,"<$files")||die"$!";
	while(my $count = <I>){
  	chomp($count);
if($count =~ m/FBgn/){
	my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
$sum_count+=$b[$f];
}}
}
#print "$sum_count\n";

open(I,"<$files")||die"$!";
	while(my $count = <I>){
  	chomp($count);
if($count =~ m/FBgn/){
	my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
$rpm = $b[$f]*1000000/$sum_count;
my $ID=$b[0]."\t".$b[1];
if (exists $RPM_warm{$ID}){
  $RPM_warm{$ID}=$RPM_warm{$ID}."\t".$rpm;
}else{
  $RPM_warm{$ID}=$rpm;
}
}}}
}

$line=0; $screen=0;
open(I,"<$files")||die"$!";
 while(my $count = <I>){
   chomp($count);
if($count =~ m/FBgn/){
 my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
  my $ID=$b[0]."\t".$b[1];
  if (exists $RPM_warm{$ID}){
    @a=@b[2..49];
if (@a == 48){
  $line++;
  my $median_count=median(@a);
  my $avg_count=average(@a);
  print OUT1 "$ID\t","$median_count\n";
    $screen++;
    print OUT "$ID\t","$RPM_warm{$ID}\n";
}}}}}
print "$line\t","$screen\n";

$outfile="OUT_IN_RPM_rRNA_removed";
open(OUT, ">$outfile");
$outfile1="OUT_IN_median_count_all";
open(OUT1, ">$outfile1");
my %RPM_OUT_IN;
for ($f = 2; $f < 50; $f++){
  #print "$f\n";
$files="Outbred_Inbred_expression_count";
my $sum_count=0; my $sum_count_length_tpm=0;
open(I,"<$files")||die"$!";
	while(my $count = <I>){
  	chomp($count);
if($count =~ m/FBgn/){
	my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
$sum_count+=$b[$f];
}}
}
#print "$sum_count\n";

open(I,"<$files")||die"$!";
	while(my $count = <I>){
  	chomp($count);
if($count =~ m/FBgn/){
	my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
$rpm = $b[$f]*1000000/$sum_count;
my $ID=$b[0]."\t".$b[1];
if (exists $RPM_OUT_IN{$ID}){
  $RPM_OUT_IN{$ID}=$RPM_OUT_IN{$ID}."\t".$rpm;
}else{
  $RPM_OUT_IN{$ID}=$rpm;
}
}}}
}

$line=0; $screen=0;
open(I,"<$files")||die"$!";
 while(my $count = <I>){
   chomp($count);
if($count =~ m/FBgn/){
 my @b = split("\t", $count);
if (exists $gene_rRNA{$b[0]}){
}else{
  my $ID=$b[0]."\t".$b[1];
  if (exists $RPM_OUT_IN{$ID}){
    @a=@b[2..49];
if (@a == 48){
  $line++;
  my $median_count=median(@a);
  my $avg_count=average(@a);
  print OUT1 "$ID\t","$median_count\n";
    $screen++;
    print OUT "$ID\t","$RPM_OUT_IN{$ID}\n";
}}}}}
print "$line\t","$screen\n";



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
