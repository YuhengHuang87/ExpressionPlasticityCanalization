#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);

my @condition=("cold","warm");
foreach my $condition (@condition){
my $outfile="total_Intron_count_gene_reads_".$condition;
open(OUT, ">$outfile");

my %hash1;my %hash2;my %hash4;my %hash3;
my $infile1="intron_count_".$condition.".txt";
for (my $j = 1; $j <= 48; $j++){
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @line = split("\r", $count);
	for my $line (@line){
		my @b = split(" ", $line);
		my $id = $b[0];
if($line =~ m/genome.bam/){
}else{
my @reg=split(":", $id);
my $trID=$reg[3]."\t".$j;
 if(exists $hash1{$trID}){
	$hash1{$trID}=$hash1{$trID}+$b[$j];
	$hash4{$trID}+=1;
 }else{
	$hash1{$trID}=$b[$j];
	$hash4{$trID}=1;
}
}
}}
close FILE1;
}

open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @line = split("\r", $count);
	for my $line (@line){
		my @b = split(" ", $line);
		my $id = $b[0];
if($line =~ m/genome.bam/){
}else{
my @reg=split(":", $id);
print OUT "$reg[3]\t";
my $ID;
for (my $j = 1; $j <= 48; $j++){
$ID = $reg[3]."\t".$j;
print OUT "$hash1{$ID}\t";
}
print OUT "$hash4{$ID}\n";
}}
}
close OUT;
}

my $outfile="total_Intron_count_gene_reads_outbred_inbred";
open(OUT, ">$outfile");
my %hash1;my %hash2;my %hash4;my %hash3;
my $infile1="intron_count_outbred_inbred";
for (my $j = 1; $j <= 48; $j++){
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);
		my $id = $b[0];
if($count =~ m/genome.bam/){
}else{
my @reg=split(":", $id);
my $trID=$reg[3]."\t".$j;
 if(exists $hash1{$trID}){
	$hash1{$trID}=$hash1{$trID}+$b[$j];
	$hash4{$trID}+=1;
 }else{
	$hash1{$trID}=$b[$j];
	$hash4{$trID}=1;
}
}
}
close FILE1;
}

open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);
		my $id = $b[0];
if($count =~ m/genome.bam/){
}else{
my @reg=split(":", $id);
print OUT "$reg[3]\t";
my $ID;
for (my $j = 1; $j <= 48; $j++){
$ID = $reg[3]."\t".$j;
print OUT "$hash1{$ID}\t";
}
print OUT "$hash4{$ID}\n";
}
}
