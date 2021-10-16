#!/usr/bin/perl
use strict;
use warnings;

my $outfile="cold_warm_comb_count_total.txt";
open(OUT, ">$outfile");

my $shared_line=0; my $cold_line=0; my $warm_line=0;

my %warm_count;
my $infile="warm_total_intron_count.txt";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my @c=split(":", $b[0]);
		my $id = $c[0]."\t".$c[1]."\t".$c[2];
    $warm_count{$id}=$count;
		$warm_line++;
  }

$infile="cold_total_intron_count.txt";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my @c=split(":", $b[0]);
		my $id = $c[0]."\t".$c[1]."\t".$c[2];
		$cold_line++;
    if (exists $warm_count{$id}){
      print OUT "$count\t","$warm_count{$id}\n";
			$shared_line++;
    }
}
close OUT;
print "$warm_line\t","$cold_line\t","$shared_line\n";

$shared_line=0; $cold_line=0; $warm_line=0;
$outfile="cold_warm_comb_RPM_rRNA_removed_median_count_200";
open(OUT, ">$outfile");
my %warm_exp;
$infile="warm_RPM_rRNA_removed";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my $id = $b[0]."\t".$b[1];
    $warm_exp{$id}=$count;
		$warm_line++;
  }

$infile="cold_RPM_rRNA_removed_median_count_200";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
    my $id = $b[0]."\t".$b[1];
		$cold_line++;
    if (exists $warm_exp{$id}){
      print OUT "$count\t","$warm_exp{$id}\n";
			$shared_line++;
    }
}
print "$warm_line\t","$cold_line\t","$shared_line\n";
