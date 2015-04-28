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

UNIFIEDGENOTYPER.sh -- Run GATK unifiedgenotyper, assumes high memory 32GB ram with 8 cores.
	Works on scarcity 6-10
	-- GATK raw variants vcf files as input, useful to rerun unifiedgenotyper without
	running entire GATK pipeline.

runHTSeqCount.py -- Run HTSeq for all gff feature types.

rnaCountAlignments.py -- Systematically search and count all mapped reads in bam files
                         associated w/ all rRNA loci.


remove3primeAFastq.py -- Remove all consecutive 3 prime A's and quality scores from fastq
          where pattern is something like:
          nnnnnnnnnnnnnnnnnAAAGATCTAAAAAAAA 
		or
          nnnnnnnnnnnnnnnAAAAAAAGATCGGAAGAGC

utilities to change chrom names:
convertChromNameFasta.py
convertChromName.py
convertChromName_SamFile.py

convertGene.py -- Create a strain specific gene sequence from a vcf file.
convertGene_Genetics.py -- same as above but runs on deepthought
	yeast_Gene_name_to_ORF.py -- required by convertGene
 	

parseSplitVCF.pl -- converts Snpeff vcf file into a more readable format, 
	   by removing some columns and optionally allowing user to subset
 	   data by strain, chromosome (one at a time) , user supplied start/stop positions

Utilities to find unique sequences in fasta file, initiall used to help kevin 
create a unique sequence file for y22-3 project.
unionFasta.py
uniqueFasta.py
findFastaDups.py
findDiffFasta.py

yeastmine.py  -- how to use python to query yeastmine, just an example.

getGeneNamesGff.py --  Create a gff file from Sean's annotation files for Kevin.

runSnpSift.py  -- Run SnpSift on a list of files to extract variant specific files.
