#!/usr/bin/perl
use English;
use Getopt::Long;
use strict;
use warnings;
use 5.10.1;

#*******************************************************************************
#    parseVCF.pl  -- converts Snpeff vcf file into a more readable format, 
#    by removing some columns and optionally allowing user to subset
#    data by strain, chromosome (one at a time) , user supplied start/stop positions
#    
#   To get a list of names from the vcf file run the following from the command line:
#    
#    grep  "^#CHROM" snpeff.vcf | perl -F'\t' -nae 'my @d = splice(@F,9); foreach(@d){ print "$_\n"; };   ' > strainNames.txt
#    This strips everything but the strain name from the column name list.
#    
#   INPUTS:  <snpEff.vcf file>  <strain list file> <Chromosome> <start pos> <stop pos>
#   1)  Snpeff vcf file to parse.  Assumed to be the output of snpeff
#   2)  List of strains to select from vcf file.   One strain name per line.  *** OPTIONAL
#   3) Chromosome to use for filtering. Use I, II, III, IV, V ..... Mito   *** OPTIONAL
#   4) Start position of desired sequence window.  ***OPTIONAL but if used you must supply a stop position
#   5) Stop position of desired sequence window.    ***OPTIONAL but if used you must supply a start position
#
#    OUTPUT:    Tab delimited output file:
#
#    CHROM	POS	REF	ALT	LOC	UPGENE	DOWNGENE	TYPE	STRAIN1	STRAIN2......
#
#    author: Mike Place  mplace@wisc.edu
#    Date:     1/22/2015
#    
#*******************************************************************************		          
#  Chromosome & Strain Name Data structures
#*******************************************************************************
my %chr =  (
"I" => 1, "II" => 1, "III" =>1,	"IV"=>1, "V"=>1,"VI" =>1,"VII"	=>1,"VIII"=>1, "IX"=>1, "X"=>1, "XI"=>1, 	
"XII"	=>1, "XIII"=>1, "XIV"=>1,"XV"=>1, "XVI"=>1, "Mito"	=>1,
);

#******************************************************************************* 
my @strainlist;				# list of stain names to filter provided by user
my %strain_column;			# hash, key = strain name, value = column index 
my @colIndex;				# final columns to get for each row
my @nameIndex;                          # strain name column headers
my @totalNames; 			# holds all strain names found in header if you must print them all

GetOptions( );

#*******************************************************************************
sub usage {
		print STDERR
qq(
Usage: parseSplitVCF.pl   <snpEff.vcf file>  <strain list file> <Chromosome> <start pos> <stop pos>

   Snpeff vcf file to parse.  Assumed to be the output of snpeff.         REQUIRED
   List of strains to select from vcf file.   One strain name per line.   REQUIRED
   
   Chromosome to use for filtering. Use I, II, III, IV, V ..... Mito      OPTIONAL
   Start position of desired sequence window.   OPTIONAL 
   Stop position of desired sequence window.    OPTIONAL , use only w/ a start position

    OUTPUT:    Tab delimited output file:
    CHROM	POS	REF	ALT	LOC	UPGENE	DOWNGENE	TYPE	STRAIN1	STRAIN2......
    

);
 }
#*******************************************************************************

if ( (scalar(@ARGV) < 2) ) {
    print STDERR "\nRequired inputs missing!\n";
    usage();
    exit(0);
};

# Command Line Vars
my $file        = $ARGV[0];		# input snpeff annotated vcf file
my $strains     = $ARGV[1];		# strain list , optional
my $chrFilter; 					# chromosome to use
my $start;						# start position
my $stop;						# stop position
my $numArgs = @ARGV;
$| = 1;                     	# turn off buffering of STDOUT

# Filter results by Chromosome
if( $numArgs == 3 ) {
     $chrFilter = $ARGV[2];
     say "Filtering on chr: $chrFilter";
     checkChr($chrFilter); 
    }

# Filter results by start position, requires stop 
if( $numArgs > 3) {
    $chrFilter = $ARGV[2];
    checkChr($chrFilter); 
    $start = $ARGV[3];
    $stop = $ARGV[4];
    if( $start &&  $stop ){
		#say "start: $start     stop: $stop ";	
		}
		else
		{  
		say  " Must provide both start and stop positions";
		exit(0);
		}
 }
 
@strainlist = getStrain($strains);	# list of stain names to filter

# open snpeff vcf file 
open my $fh, '<', $file or die "Could not open file  ' $file  ' :  $! ";

# process snpeff file line by line
while ( my $line = <$fh> ) {
		if ( $line =~ /^##/ ) {     
				next;         						# skip header information 
		} 
		elsif ( $line =~ /^#CHROM/ ){
			    chomp( $line );		# Get column names
			    ####$line =~ s{//}{/}g;
			    my @header =  split ('\t', $line);
				my @c = split('\t', $line);
                my @col = splice( @c, 9);
				foreach(@col) { 
					state $index = 8;
					#my @name = split('\/', $_);		
					my $item = $_;
					push(@totalNames, $item) ;
					$strain_column{$item} = $index;			
					$index++;		
					}
              #print first part of header
              print "$header[0]\t$header[1]\t$header[3]\t$header[4]\tINFO\tFORMAT\t";

              # print strain names 
             if( @strainlist  ) {
				 foreach ( @strainlist ) { 
				           push(@colIndex, int( $strain_column{$_} ) );
				           push(@nameIndex, $_);
	  	            }
	  	            print  	join ("\t", @nameIndex), "\n";
		      }
		       else {
				    my @list;
				     foreach( @totalNames ) {
								push(@list, $_);
						 }
						   print  	join ("\t", @list), "\n";
				   }
		       
			
		}
		else {															#  here after parse each line for required information.
			 my @data = split( '\s', $line );
			 splice( @data, 2, 1);
			 
			 # Start Filtering data
			 if( $chrFilter ) {                     					# Print Start and Stop positions for a chromosome
                    if( $start ) {
						if(  $data[0]  eq $chrFilter  && $data[1] >= $start && $data[1] <= $stop) {
							    my @result = getAnnotation( $data[6] ); 
								print "$data[0]\t$data[1]\t$data[2]\t$data[3]\t";                
								print join (",", @result), "\t$data[7]\t";   
							 
								foreach( @colIndex) {
									#print "$data[$_]\t";
									my $gt = recodeAllele( $data[2],$data[3],$data[$_] );
								    print "$gt\t";
								}
								print "\n";
							}
			 	
			 	} 
				elsif( $data[0]  eq $chrFilter ) {                      # Print only the selected chromosome
						my @result = getAnnotation( $data[6] ); 
						print "$data[0]\t$data[1]\t$data[2]\t$data[3]\t";    
						print join (",", @result), "\t$data[7]\t"; 
						
						foreach( @colIndex) {
							#print "$data[$_]\t";
							my $gt = recodeAllele( $data[2],$data[3],$data[$_] );
							print "$gt\t";
							}
						print "\n";
						  
				}		
				
			  }else {													# print all rows for chosen strains
						
						if( $strains eq "all") {
							print join("\t", @data); 
							print "\n";
							}
							else { 
						    my @result = getAnnotation( $data[6] ); 
						    print "$data[0]\t$data[1]\t$data[2]\t$data[3]\t";   
						    print join (",", @result), "\t$data[7]\t"; 
						    
						    foreach( @colIndex) {
							my $gt = recodeAllele( $data[2],$data[3],$data[$_] );
							print "$gt\t";
							#print "$data[$_]\t";
							}
						print "\n";
					}
			  	}
			}
	}
	

#*******************************************************************************
# Functions
#*******************************************************************************
#  getStrain() -- accepts file of strain names as inputs. File has single strain name per line.
#                              returns an array of strain names
#*******************************************************************************
sub getStrain  {
	my $sFile = shift;
	my @strainlist;
	
if( $sFile  ) {
		open my $lh, '<', $sFile or die "Could not open file  ' $sFile ' : $!";
		#say "selecting for strains from file: $sFile";
		while( <$lh> ) {
				if ($_ =~ /^$/) { 
					next;
				}
				else { 
					chomp($_);
					push( @strainlist, $_);
				}
			}
	}
	else {
				$strains = "all";
		}
return @strainlist;
}


#*******************************************************************************
# Check if the chrom filter is a valid character, must be a roman numeral
#*******************************************************************************
sub checkChr { 
	my $ch = shift;
	
    if ( $ch ) {
        if( not exists( $chr{$ch} ) ) {
                print "Chromosome name invalid,  use I, II, III, IV.... Mito \n";
                exit(0);
        }
        else {
                #print "Filtering on  chr $chrFilter  \n";
        }
    }
}

#*******************************************************************************
# split snpeff ANN field and return required fields
#*******************************************************************************
sub getAnnotation {
  my $ann     = shift ;
  my @ln      = split(',', $ann);
  my $ctAnn   = @ln;
  my $result ;
  my $ct = 0;

		foreach(@ln) {  
			   my @field = split('\|', $_);
			   if (!$field[3]) {
				$field[3] = "NONE";
			   }
			   
			   $result  .=  $field[1].",".$field[2].",".$field[3];   # remove 7
						$ct++;
			   if ($ct < $ctAnn ) {
								$result .= ":";
						}
			}
		return $result;		
	}

#*******************************************************************************
# recode the genotype field for each strain
# GATK genotype field: 0 = REF allele 1 = ALT allele 2= 2nd ALT allele etc...
# args: REF and ALT alleles and strain genotype column
#*******************************************************************************
sub recodeAllele {
  my $ref     = shift;				# Reference allele
  my $alt	  = shift;				# alternate alleles, 1 or more 
  my @alt	  = split(',',$alt);	# 
  my $pheno   = shift;				# strain genotype data: 1/1:0,69:69:99:1943,183,0
  my @ln      = split(':', $pheno); 
  my $result;
  my $ct      = 1; 
  my %gt;							# key = allele (G,C,T,A) value = 0,1,2,3
  
  if( not exists $gt{0} ) {  $gt{0} = $ref; }
  
  foreach( @alt ){
				if ( not exists $gt{$ct} ) {
				    $gt{$ct} = $_;
					$ct++;
				}
  }
  
  if ($pheno =~ /^\.\/\./) {                    # if genotype column data is ./.:.:.:.:. (empty) return 
                return $pheno;
  }
  
  my @code = split( '/', $ln[0] );				# break up the 0/0 field
  my $ck = @code;
    
  my $newGT = $gt{$code[0]}."/".$gt{$code[1]};  # convert the field number to genotype 
  shift(@ln);						# remove the old genotype field
  unshift(@ln,$newGT);				# add the new genotype coding field
  
  $result = join(':', @ln );
  return $result;
  
	}












