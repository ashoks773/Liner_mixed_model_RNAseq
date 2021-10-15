#!/usr/bin/perl -w
  
use strict;
use Data::Dumper;

my $input = <$ARGV[0]>;
open (IN, $input) || die "Cannot opent he file1\n";
my $hash1 = {};
while (chomp(my $line = <IN>))
{
#ProCoding_Gene_counts_mitoremove_Riboremove_Filtered.txt (Third column is Gsymbol)                
my @arr = split("\t", $line);
                my $id = $arr[0];
		my $symbol = $arr[2];
		$hash1->{$id}="$symbol";
}
print Dumper $hash1;             

my $input1 = <$ARGV[1]>;
open (IN1, $input1) || die "Cannot opent he file2\n";
open (OUT1, ">$input1.Symbols");
while (chomp(my $line1 = <IN1>))
{               
        my @arr1 = split ("\t", $line1);
        my $eid = shift @arr1;
        if (exists ($hash1->{$eid}))
        {
		print OUT1 $hash1->{$eid}."\t$line1\n";
        }
}

