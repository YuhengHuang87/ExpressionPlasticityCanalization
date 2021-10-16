#!/usr/bin/perl
use strict;
use warnings;

my $derived_copy=5;
my $outfile="cold_total_intron_count.txt";
open(OUT, ">$outfile");
my %total_count_cold;
my $infile="total_Intron_count_gene_reads_cold";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
    my $id = $b[0];
    $total_count_cold{$id}=$count;
  }

$infile="intron_count_cold";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @line = split("\r", $count);
		for my $line (@line){
		my @a=split(" ", $line);
		my @b=split(":", $a[0]);
		my $id = $b[3];
    if (exists $total_count_cold{$id}){
			my @sum_type = split("\t", $total_count_cold{$id});
			my $avg_total=average(@sum_type[1..48]);

			my $f=0;
			foreach my $sum_t (@sum_type[1..48]){
			if ($sum_t>0){
			$f++;
			}}

			if ($f == 48){
			my $avg_intron = average(@a[1..48]);
			my $alter=$avg_total-$avg_intron;
			my @list_intron=($avg_intron,$alter);

			my @sorted_intron = sort { $a <=> $b } @list_intron;
			if ($sorted_intron[0]>=$derived_copy){

			foreach my $a (@a){
				print OUT "$a\t";
			}
      print OUT "$total_count_cold{$id}\n";
    }
}
}}}
close OUT;

$outfile="warm_total_intron_count.txt";
open(OUT, ">$outfile");
my %total_count_warm;
$infile="total_Intron_count_gene_reads_warm";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
    my $id = $b[0];
    $total_count_warm{$id}=$count;
  }

$infile="intron_count_warm";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @line = split("\r", $count);
		for my $line (@line){
		my @a=split(" ", $line);
		my @b=split(":", $a[0]);
		my $id = $b[3];
    if (exists $total_count_warm{$id}){
			foreach my $a (@a){
				print OUT "$a\t";
			}
			print OUT "$total_count_warm{$id}\n";
			}
			}}


			my $outfile="outbred_inbred_total_intron_count.txt";
			open(OUT, ">$outfile");

			my %total_count;
			my $infile="total_Intron_count_gene_reads_outbred_inbred";
			open(FILE1,"<$infile")||die"$!";
				while(my $count = <FILE1>){
					chomp($count);
					my @b = split("\t", $count);
			    my $id = $b[0];
			    $total_count{$id}=$count;
			  }

			$infile="intron_count_outbred_inbred";
			open(FILE1,"<$infile")||die"$!";
				while(my $count = <FILE1>){
					chomp($count);
					my @a=split(" ", $count);
					my @b=split(":", $a[0]);
					my $id = $b[3];
			    if (exists $total_count{$id}){
						foreach my $a (@a){
							print OUT "$a\t";
						}
			      print OUT "$total_count{$id}\n";
			    }
			}

sub average {
my @array = @_; # save the array passed to this function
my $sum; # create a variable to hold the sum of the array's values
foreach (@array) { $sum += $_; } # add each element of the array
# to the sum
return $sum/@array; # divide sum by the number of elements in the
# array to find the mean
}
