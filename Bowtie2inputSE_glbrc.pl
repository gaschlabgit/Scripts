#!/usr/bin/perl 

use warnings;
use strict;
use diagnostics;
use Cwd;

#  This script will parse the files in an directory and call bowtie2 on each fastq file in the directory.  All
#  *.fastq files are stored in @INfiles.  The script then loops through the array and calls my bowtie2 script
#  (runBowtie2_genetics_se.sh) on each file. A directory inside the input directory called Bowtie2 will store
#  all the bowtie2 output, including the resulting *.bam file and a text file containing the std out from bowtie2.
#
# Modified Feb 2015 by Mike Place
# This combines Katie's original bash (runBowtie2_glbrc_se.sh) and perl script (Bowtie2input_glbrc_se.pl) into one script.
#
######################################################################################################################
#  Input Arguments and Variables

unless (scalar @ARGV == 1) {
  print "\nUSAGE:  Bowtie2input.pl <Reference.basename>\n";
  print "Assumes input directory is current directory which\n";
  print "fastq files and runs runBowtie2 shell script on each file. Also input the\n";
  print "basename of the reference index files.\n";
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

#  Create output direcotry to store all bowtie2 output files
my $OUTPUT_DIR = "$INPUT_DIR/Bowtie2/";
die "Output directory '$OUTPUT_DIR' already exists.  To continue, rename existing directory!\n" if (-d $OUTPUT_DIR);

unless (mkdir $OUTPUT_DIR) {
	die "Could not make output directory '$OUTPUT_DIR'\n$!\n";
}



######################################################################################################################
#  Body of Script

$|=1;									#Autoflush standard out
&main();								#Script entry point


######################################################################################################################
#  Main subroutine

sub main {
	bowtie2 ();	
	
}

######################################################################################################################
#  bowtie2 subroutine
#  This subroutine loops through all fastq files in the directory and sends them to runBowtie2_genetics_se.sh shell
#  script.  This scripte requires the reference index files, the fastq file, and the output base name.  The output
#  base name is created from the name of the fastq file.  The output .bam file and std out from bowtie2 is saved in
#  an output directory: InputDirectory/Bowtie2/



sub bowtie2 {
	
	foreach my $file (@INfiles) {
		my $name;
		if ($file =~ /(.*)\.fastq$/) {
			$name = $1;
		} else {
			die "$file is not a *.fastq file.\n";
		}
		my $out = "$OUTPUT_DIR"."$name";
		open( my $fh,'>>', "$OUTPUT_DIR/Bowtie_run.log" ) or die "Unable to write to file $out.run.log : $!";
		print $fh "Running bowtie2 on '$file'\n";
		close($fh);
		print "Running bowtie2 on '$file'\n";
		
		# Run Bowtie2 
		my $btlog = `bowtie2 -p 8 --phred33 -N 1 -x $REFERENCE -U $file -S $out.sam 2>>$OUTPUT_DIR/Bowtie_run.log`;
		
		###    Clean the SAM file
		###    This soft-clips an alignment that hangs off the end of its reference sequence.
		###    This will print out all the errors that it ignores (MAPQ errors
		print"\n\n I=out.sam = $out.sam     O= $out.cleaned.sam\n";
		my $cleanSam = `java -Xmx8g -jar /opt/bifxapps/picard-tools/CleanSam.jar I=$out.sam O=$out.cleaned.sam`;
		
		# REMOVE original sam file  #rm $OUT.sam
		unlink("$out.sam");
		
		###    Add the RG header and sort the SAM file
		###    This will print out all the errors that it ignores (MAPQ errors)
		my $sortSam  = `java -Xmx8g -jar  /opt/bifxapps/picard-tools/AddOrReplaceReadGroups.jar I=$out.cleaned.sam O=$out.final.sam SO=coordinate LB=$REFERENCE.fasta PL=ILLUMINA PU=unknown SM=$out VALIDATION_STRINGENCY=LENIENT`;
		
		# REMOVE cleaned sam file
		unlink("$out.cleaned.sam");
		
		###    Make the BAM file, sort and index it
		my $makeBam  = `samtools view -uS -t $REFERENCE.fsa.fai $out.final.sam | samtools sort - $out.sorted`;
		
		### Index bam file
		my $indexBam = `samtools index $out.sorted.bam`;
		unlink("$out.final.sam");
		
		
		
		print "Finished running bowtie2 on '$file'\n\n";
	}
	
}

######################################################################################################################
