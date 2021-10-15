#!/usr/bin/perl -w
#- Get MBCluster Clustering Info
 
use strict;
use Data::Dumper;

my $input = <$ARGV[0]>;
open (IN, $input) || die "Cannot opent he file1\n";
my $hash1 = {};
while (chomp(my $line = <IN>))
{
# Open File which contains Cluster Information - (/Users/sharmaa4/Box/David_Casero_Lab/Laxmi_Yeruva/RNA_seq/Complex_ana_Final/Clustering_on_Sel_Genes/LI_Gname_cluter_FC.txt)
#Ganmes  clsTargetsALL.cluster   X1      X2      X3      X4      X5      X6      X7      X8
#Sox17   8       0.588721283418646       0.127391326243468       0.219391252682175       0.408554539449426       -0.105605401603504      -0.402337849757953      -0.406426620647327      -0.429688529784933
	my @arr = split("\t", $line);
	$hash1->{$arr[0]}="$arr[1]";
}
print Dumper $hash1;             

my $input1 = <$ARGV[1]>;
#--- Open Load Average TPM File (LI_TPMs.txt)
open (IN1, $input1) || die "Cannot opent he file2\n";
open (OUT1, ">LI_cluster_genes_TPM.txt");
while (chomp(my $line1 = <IN1>))
{               
        my @arr1 = split ("\t", $line1);
        if (exists ($hash1->{$arr1[0]}))
        {
		print OUT1 $hash1->{$arr1[0]}."\t$line1\n";
        }
	else
	{
		print OUT1 "C0\t$line1\n";
	}
}

