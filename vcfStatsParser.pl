#!/usr/bin/perl 
use Data::Dumper;
use warnings;
use strict;
#******************************************************************************
# vcfStatsParser.pl -- Transform vcf-stats dumper file into a results table.
# 
# http://vcftools.sourceforge.net/perl_module.html#vcf-stats
# It is assumed that you will be using the vcf-stats file.vcf.gz > dumpfile 
# command results.  These results are produced by perl Data::Dumper.     	
#
# Input : dump file to convert to table.
# Output: print tab delimited table where to stdout.
#         rows = sample names, columns = type 
#
# author: Mike Place  mplace@wisc.edu, Gasch Lab
# Date  : May 1, 2015
#
#******************************************************************************
sub usage {
    print STDERR qq (
	Usage: vcfStatsParser.pl <vcf-stats file>
	
	Input 'dump' file to convert to table.
	Input file generated via 'vcf-stats file.vcf.gz > dumpfile' cmd.
	
	Output print tab delimited table to stdout, 
	rows = sample names, columns = type

	Sample	count het_AA_count	het_RA_count	etc.....
	YJM1526	195	18	23	...........

	Try again!
);
}

# Check for command line inputs
if ( (scalar(@ARGV) < 1) ) {
	print STDERR "\n Required inputs missing.\n";
	usage();
	exit(0);
};

# Command Line Vars
my $inputFile = $ARGV[0];         # name of vcf-stats dump file

# if asked for help give it.
if( $inputFile =~ /^-h/ ) {
	 usage();
         exit(0);
}

$Data::Dumper::Purity = 1;      # set to fill self-referenced data structure
open (FILE, "<$inputFile")  || die "$! $@";

my $file;
{
	undef $/;
	$file = <FILE>;
}

my $VAR1;
eval $file;   # call to eval recreates the data structure associated with $VAR1

# get a list of sample names
my @sampleNames;
for my $k ( keys %{ $VAR1->{'samples'} } ) {
	push(@sampleNames, $k);
}

# find a complete list column headers
my %uName;
for my $k ( keys %{ $VAR1->{'samples'} } ) {
	my @r; 
	my @finalRow;
	for my $ks ( sort keys %{ $VAR1->{'samples'}->{$k} } ) {
	       if( $ks =~ /^indel$/) { next; }
               if( $ks =~ /^snp$/ )  { next; }
	       push(@r, $ks );
	       if( not exists $uName{$ks} ) {
			$uName{$ks} = 1;
		}else {
			$uName{$ks} += 1;
		}		
	}
}

# Get a sort list of column names
my @colNames;
for my $k ( sort keys (%uName) ) {
	push(@colNames, $k);
}

my @header = @colNames;
unshift( @header, 'sample' );
print join("\t",@header); 
print "\n";

for my $k ( keys %{ $VAR1->{'samples'} } ) {
	print "$k\t";
	my @row;
	foreach my $name (@colNames) {
	my $ct;
	if( not exists $VAR1->{'samples'}->{$k}{$name} ) {
	    $ct = 'NA';
	}else {
	    $ct = $VAR1->{'samples'}->{$k}{$name};
	}
	push( @row, $ct);
	}

	print join("\t", @row);
	print "\n";
}




