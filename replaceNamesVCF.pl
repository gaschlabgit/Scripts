#!/usr/bin/perl
use strict;
use warnings;
use 5.10.0;
#*******************************************************************************
#    replaceNamesVCF.pl  --  Replace the vcf file strain names and converts
#                            ref|NC_001133| names to roman numerals.
#	INPUT:  
#         vcf file  --standard vcf file
#         csv file  --comma separated file, one line per name combination
#		    ex:	  name,nametoReplace
#                             
#    OUTPUT:    
#		standard vcf file with strain names replaced
#
#    author: Mike Place  mplace@wisc.edu
#    Date:     1/22/2015
#    
#*******************************************************************************		          
sub usage {
		print STDERR
qq(
Usage:  		replaceNamesVCF.pl  <vcfFile> <namesFile>
	INPUT:  
         vcf file  --standard vcf file
         csv file  --comma separated file, one line per name combination
				   	 name,nametoReplace
                             
    OUTPUT:    
		standard vcf file with strain names replaced
    );
 }
#*******************************************************************************
# chromosome name conversion table
my %chr = ("ref|NC_001133|"=>"I","ref|NC_001134|"=>"II","ref|NC_001135|"=>"III","ref|NC_001136|"=>"IV",
"ref|NC_001137|"=>"V","ref|NC_001138|"=>"VI","ref|NC_001139|"=>"VII","ref|NC_001140|"=>"VIII",
"ref|NC_001141|"=>"IX","ref|NC_001142|"=>"X","ref|NC_001143|"=>"XI","ref|NC_001144|"=>"XII",
"ref|NC_001145|"=>"XIII","ref|NC_001146|"=>"XIV","ref|NC_001147|"=>"XV","ref|NC_001148|"=>"XVI",
"ref|NC_001224|"=>"Mito",
);

#*******************************************************************************
if ( (scalar(@ARGV) < 2) ) {
    print STDERR "\nRequired inputs missing!\n";
    usage();
    exit(0);
};

my $vcfFile      = $ARGV[0];                     # input vcf file
my $nameFile     = $ARGV[1];                     # strain list , optional
my %strainMap;									 # Strain Name Data structure

# open the strain names file and get name information
%strainMap =  getStrain($nameFile);
 
# open snpeff vcf file 
open my $fh, '<', $vcfFile or die "Could not open file  ' $vcfFile  ' :  $! ";

while( my $line = <$fh> ) {
		if ( $line =~ /^##/ ) {     
				print $line;
		}
		elsif ( $line =~ /^#CHROM/ ) {
			    chomp( $line );		# Get column names
                my @h = split('\t', $line);
			    #my @header   = splice( [split('\t', $line)],0,9);
				my @header = splice( @h,0,9 );
                my @c = split('\t', $line);
                my @colNames = splice( @c, 9);
                
                foreach( @header ){ print "$_\t"; }
                foreach ( @colNames ){
                print "$strainMap{$_}\t";          
                }
                say "";
        }
        else {  #replace ref|NC_001133| type names with roman numerals
         		#print $line;
                my @row = split('\t',$line);
                my $chrom = $chr{$row[0]};
                shift(@row);
                print "$chrom\t";
                print join("\t", @row );                
                
		}
	} 
      
#*******************************************************************************
# Functions
#*******************************************************************************
#  getStrain() -- accepts file of strain names as inputs.
#                 returns an hash of strain names
#*******************************************************************************
sub getStrain  {
	my $sFile = shift;
	my %strainlist;
	
if( $sFile  ) {
		open my $lh, '<', $sFile or die "Could not open file  ' $sFile ' : $!";
		while( <$lh> ) {
				if ($_ =~ /^$/) {         # skip blank lines
					next;
				}
				else { 
					chomp($_);
					my @n = split(',', $_);
					$strainlist{$n[1]} = $n[0];
				}
			}
	}
return %strainlist;
}















