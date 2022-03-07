#!/bin/perl
use strict;
use warnings;
#usage
# cat file.gff | perl replace_contig_names.pl names_info.txt > file_fixed.gff
# cat genome.fa | perl replace_contig_names.pl names_info.txt > file_fixed_genome.fa
my $info_file = $ARGV[0];
my %hash;
open INFO, $info_file;
while(<INFO>){
  chomp;
  my @a = split(/,/,$_);
  $hash{$a[0]} = $a[1];
}
while(<STDIN>){
  chomp;
  my $line = $_;
  foreach my $i (sort keys %hash){
    $line =~ s/$i/$hash{$i}/g;
  }
  print "$line\n";
}
