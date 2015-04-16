#!/usr/bin/perl 

use warnings;
use strict;
use diagnostics;
use Cwd;

#  This script will parse the files in an directory and call bwa on each fastq file in the directory.  All
#  *.fastq files are stored in @INfiles.  A directory inside the input directory called bwa will store
#  all the bwa output, including the resulting *.bam file and a text file containing the std out from bwa.
#
# Modified Feb 2015 by Mike Place
#
######################################################################################################################
#  Input Arguments and Variables

unless (scalar @ARGV == 1) {
  print "\nUSAGE:  BwaMem2inputPE.pl <Reference.basename>\n";
  print "Assumes input directory is current directory which\n";
  print "fastq files and runs bwa mem on each pair of fastq files.\n";
  print "basename of the reference bwa index file required.\n";
  exit;
}

#my $INPUT_DIR = shift;
my $INPUT_DIR = getcwd();
my $REFERENCE = shift;
die "Undefined input directory!\n" unless defined $INPUT_DIR;

print "input dir:  $INPUT_DIR \n";

opendir my ($IN), $INPUT_DIR or die "Cannot open input directory '$INPUT_DIR'\n";

#  Collect all fastq files from the directory and store in @INfiles
my @INfiles = grep { /\.fastq$/} readdir $IN;
die "No input .fastq files found!\n" unless @INfiles > 0;

#  Create output direcotry to store all bwa output files
my $OUTPUT_DIR = "$INPUT_DIR/Bwa/";
die "Output directory '$OUTPUT_DIR' already exists.  To continue, rename existing directory!\n" if (-d $OUTPUT_DIR);

unless (mkdir $OUTPUT_DIR) {
	die "Could not make output directory '$OUTPUT_DIR'\n$!\n";
}

my %PAIRS;

######################################################################################################################
#  Body of Script

$|=1;									#Autoflush standard out
&main();								#Script entry point

######################################################################################################################
#  Main subroutine
sub main {
  
	bwa ();
	
}
######################################################################################################################
#  bwa2 subroutine
#  This subroutine loops through all fastq files in the directory and runs bwa mem.
#  This script requires the reference index files, the fastq file, and the output base name.  The output
#  base name is created from the name of the fastq file.  The output .bam file and std out from bwa is saved in
#  an output directory: InputDirectory/Bwa/

sub bwa {
		 # Find matching pairs and store them in a hash.
		foreach my $file (@INfiles) {
		  print $file."\n";
			my $name;
			if ($file =~ /(.*)_1\.fastq$/) {
			$name = $1;
			my $file2 = "$name"."_2".".fastq";
			push @{$PAIRS{$name}}, $file, $file2;
			}
		}
		
		open(my $ph, '>>', "$OUTPUT_DIR/picard.log") or die "Unable to write to picard.log : $!";
		open(my $sh, '>>', "$OUTPUT_DIR/samtools.log") or die "Unable to write to samtools.log : $!"; 		
		 
	  foreach my $strain (keys %PAIRS) {
		my $file1 = "$INPUT_DIR/"."$PAIRS{$strain}[0]";
		my $file2 = "$INPUT_DIR/"."$PAIRS{$strain}[1]";
		my $out = "$OUTPUT_DIR"."$strain";
		my $fh;
		open( $fh,'>>', "$OUTPUT_DIR/Bwa_run.log" ) or die "Unable to write to file $out.run.log : $!";
	    print $fh "Running bwa mem on '$file1'\n";
		close($fh);
		
		# Run bwa 
		my $btlog = `bwa mem -t 8 -M $REFERENCE $file1 $file2 1> $out.sam 2>>$OUTPUT_DIR/Bwa_run.log`;
		
		###    Clean the SAM file
		###    This soft-clips an alignment that hangs off the end of its reference sequence.
		###    This will print out all the errors that it ignores (MAPQ errors
		print $ph "\n\n I=out.sam = $out.sam     O= $out.cleaned.sam\n";
		my $cleanSam = `java -Xmx8g -jar /opt/bifxapps/picard-tools/CleanSam.jar I=$out.sam O=$out.cleaned.sam 2>&1`;
		print $ph "$cleanSam";
		
		# REMOVE original sam file  #rm $OUT.sam
		unlink("$out.sam");
		
		###    Add the RG header and sort the SAM file
		###    This will print out all the errors that it ignores (MAPQ errors)
		my $sortSam  = `java -Xmx8g -jar /opt/bifxapps/picard-tools/AddOrReplaceReadGroups.jar I=$out.cleaned.sam O=$out.final.sam SO=coordinate LB=$REFERENCE.fasta PL=ILLUMINA PU=unknown SM=$out VALIDATION_STRINGENCY=LENIENT 2>&1`;
		print $ph "$sortSam";
		
		# REMOVE cleaned sam file
		unlink("$out.cleaned.sam");
		
		###    Make the BAM file, sort and index it
		my $makeBam  = `samtools view -uS -t $REFERENCE.fsa.fai $out.final.sam | samtools sort - $out.sorted 2>&1`;
		print $sh "$makeBam";
		
		### Index bam file
		my $indexBam = `samtools index $out.sorted.bam 2>&1`;
		print $sh "$indexBam";
		
		unlink("$out.final.sam");
		
		
		
		print "Finished running bowtie2 on '$file1'\n\n";
	}
	  close($ph);
	  close($sh);
}

######################################################################################################################
