#!/bin/perl
use warnings;
use strict;

#This script takes the *.out files of syri and creates an output that can be read into R by putting the larger block info onto the smaller block info (i.e. which larger block is each smaller block in)

print "chr_1\tregion_start_1\tregion_end_1";
print "\tchr_2\tregion_start_2\tregion_end_2";
print "\tblock_start_1\tblock_end_1\tblock_start_2\tblock_end_2";
my $region_string;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($_ =~ "^#"){
    $region_string = "$a[1]\t$a[2]\t$a[3]\t$a[5]\t$a[6]\t$a[7]";
  }else{
    print "\n$region_string\t$a[0]\t$a[1]\t$a[2]\t$a[3]";
  }
}
