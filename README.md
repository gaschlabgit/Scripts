# Scripts 
# 

Bowtie and bwa -- combo of the Gasch lab perl and bash scripts
		 ( runBowtie2_glbrc_se.sh Bowtie2input_glbrc_se.pl)

	Bowtie2inputPE_glbrc.pl
	BwaMeminputSE_glbrc.pl
	Bowtie2inputSE_glbrc.pl
	BwaMeminputPE_glbrc.pl
	
	Bowtie2inputPE.pl
	Bowtie2inputSE.pl
	BwaMeminputPE.pl
	BwaMeminputSE.pl

UNIFIEDGENOTYPER.sh -- Run GATK unifiedgenotyper
	
	assumes high memory 32GB ram with 8 cores.
	
	Works on scarcity 6-10
	
	-- GATK raw variants vcf files as input, useful to rerun unifiedgenotyper without
	running entire GATK pipeline.

runHTSeqCount.py  
		
		Run HTSeq for all gff feature types.

rnaCountAlignments.py  

		 Systematically search and count all mapped reads in bam files
                 associated w/ all rRNA loci.


remove3primeAFastq.py 
	
	 Remove all consecutive 3 prime A's and quality scores from fastq
        
	  where pattern is something like:
        
	  nnnnnnnnnnnnnnnnnAAAGATCTAAAAAAAA 
	
		or
        
	  nnnnnnnnnnnnnnnAAAAAAAGATCGGAAGAGC

utilities to change chrom names:
	
	convertChromNameFasta.py

	convertChromName.py
	
	convertChromName_SamFile.py

convertGene.py 

	Create a strain specific gene sequence from a vcf file.

convertGene_Genetics.py 
	
	same as above but runs on deepthought
	
	yeast_Gene_name_to_ORF.py -- required by convertGene
 	

parseSplitVCF.pl 
	
	converts Snpeff vcf file into a more readable format, 
	
	by removing some columns and optionally allowing user to subset
 	
	data by strain, chromosome (one at a time) , user supplied start/stop positions

Utilities for fasta files, initially used to help kevin 
create a unique sequence file for y22-3 project.

	unionFasta.py

	uniqueFasta.py

	findFastaDups.py

	findDiffFasta.py
	
	sequenceLength.py 
	
		input a fasta file, will print sequence name and length to stdout

yeastmine.py  

	how to use python to query yeastmine, just an example.

getGeneNamesGff.py 
	
	Create a gff file from Sean's annotation files for Kevin.

runSnpSift.py  

	Run SnpSift on a list of files to extract variant specific files.

vcfStatsParser.pl 

	Convert vcf-stats dump file into a table.

getPrivateSnps.py 

	Run vcf-contrast on vcf files, to get private variant vcf file for sample.

calcHeterozygosityVCF.py 

	Calculate % Heterozygosity on a list of private vcf files.

coverageAnalysis.py

	Copy Number Analysis, simple count reads for 1kb non-overlaping bins.

fixFastaName.py 

	Create a new fasta file, renaming each chrom. 

coverageAnalysis.py 

	Calculate Copy Number Variants using a sorted bam file.
	Uses log2(window Count/median of count for chromosome)


matchLinesIn2files.py

	Match 2 lines in 2 files. 
	Easy to modify to alter matching criteria.

runCntrlFreec.py

	Automates the use of control-freec for CNV calling.

parseBlastn.py

	Use BioPython to parse Blastn xml results files.

extractGene.py
	
	Extract sequence plus window from genome assembly

blastnExtract.py

	Extract a gene plus & minus user defined window from a sequence assembly
	by parsing blastn results generated using -outfmt 6 


