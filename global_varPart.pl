#!/usr/bin/perl -w
#--- Arrange all genes (16494) in the Master File. For global VariancePartition results
 
use strict;
use Data::Dumper;

my $input = <$ARGV[0]>;
open (IN, $input) || die "Cannot opent he file1\n";
my $hash1 = {};
while (chomp(my $line = <IN>))
{
#Open master file
	$hash1->{$line}=$line;
}
print Dumper $hash1;             

my $input1 = <$ARGV[1]>;
#--- Open Global Variance Partition result
open (IN1, $input1) || die "Cannot opent the file2\n";
open (OUT1, ">$input1.master.txt");
while (chomp(my $line1 = <IN1>))
{
#Gnames  Condition       Epithelial_cell Immunization    Isolator        Mouse   SAC     Sex     Tissue.type     Residuals
#Chtop   0.32232594      0.119877267     0.006597329     9.92E-11        0.005597541     0.483122256     9.62E-11        0.032515147     0.02996452               
        my @arr1 = split ("\t", $line1);
        if (exists ($hash1->{$arr1[0]}))
        {
		print OUT1 "$line1\n";
        }
#	else
#	{
#		print OUT1 "$hash1->{$arr1[0]}"."NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
#	}
}

