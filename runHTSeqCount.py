#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
Program: runHTSeqCount.py 

Purpose: run HTseq-count on a list of bam files

Input  : none
 default gff:
  gff    = "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208_noFasta.gff"
   
Output: HTSeq results for each gff feature type.

*******************************************************************************
HTSeq

Given a file with aligned sequencing reads and a list of genomic features, 
count how many reads map to each feature.

example: 
    htseq-count -t CDS -i Parent samFile  gff 

    OUTPUT: htseq text file

http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
*******************************************************************************
@author: mplace
"""

import os
import re
import sys
import subprocess      
import argparse       
import logging

# list of gff feature types
feature = { "ARS", "binding_site", "CDS", "centromere", "chromosome", "external_transcribed_spacer_region",
            "five_prime_UTR_intron", "gene", "gene_cassette", "insertion", "internal_transcribed_spacer_region",
            "intron", "long_terminal_repeat", "LTR_retrotransposon", "ncRNA", "noncoding_exon", "nucleotide_match",
            "pseudogene", "region", "repeat_region", "rRNA", "snoRNA", "snRNA", "telomere", "transposable_element_gene",
            "tRNA" }     
 
def runHTSeq( files, ftype, gff, rvse=0 ):
    """
    Call HTSeq-count
    files = list of files to process
    ftype = sam or bam
    rvse  = 0 (no reverse) or 1 (reverse)
    """
    program  = "/opt/bifxapps/python/bin/htseq-count"
        
    parent = { "binding_site", "external_transcribed_spacer_region", "five_prime_UTR_intron",  "insertion", 
               "internal_transcribed_spacer_region", "intron", "noncoding_exon", "nucleotide_match", "region", 
               "repeat_region" }
               
    for i in files:
        if ftype == 'bam':
            samName  = re.sub(r"bam", "sam", i )
            samcmd  = [ 'samtools', 'view', '-o', samName, i ]
            output  = subprocess.Popen( samcmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
        else:
            samName = i
        
        for item in feature:
            
            if item in parent:
                term = 'Parent'
            else:
                term = 'ID'
                
            htseqOut    = i + "_" + item + "_HTseqOutput.txt"  
            if rvse == 1:
                cmd     = [ program, '-t', item, '-i', term, '-s', 'reverse' , samName, gff ]
            else:
                cmd     = [ program, '-t', item, '-i', term , samName, gff ]
            
            logging.info(' htseq commands used: ')
            logging.info(' %s' %(cmd))
            output  = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
            result1 = output[0].decode( 'utf-8' )
        
            with open( htseqOut, 'w' ) as out:
                out.write( result1 )      

def main():
    cmdparser = argparse.ArgumentParser(description="Run HTSeq for all gff feature types.",
                                        usage='%(prog)s [optional args: -g -i -r -s ]' ,prog='runHTSeqCount.py'  )
    cmdparser.add_argument('-g', '--gff',       action='store',      dest='GFF',       help='GFF file to use with HTSeq, full path to file.')                                    
    cmdparser.add_argument('-s', '--samFiles' , action='store_true', dest='SAM',       help='Use sam files instead of bam files.')
    cmdparser.add_argument('-r', '--reverse',   action='store_true', dest='REVERSE',   help='HTSeq reverse parameter, anti-sense library preps, optional.')
    cmdparser.add_argument('-i', '--info',      action='store_true', dest='INFO',      help='Print a more detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
        
    if cmdResults['INFO']:
        print("\n  runHTSeqCount.py")
        print("\n  Purpose: Run HTSeq-Count on every bam/sam file in directory.")
        print("\n  Input  : gff file to use, program assumes you are in the directory with bam/sam files.")
        print("\n  Usage  : Move to the directory with bam/sam files.\n")
        print("\trunHTSeqCount.py \n")
        print("  Optinal arguments:")
        print("\t-g , full path to gff file")
        print("\t-r , use for reads that map to the anti-sense strand")
        print("\t-s , use sam files instead of bam files")
        print("")
        print("\tSee Mike Place for any problems or suggestions.")
        sys.exit(1)
        
    if cmdResults['GFF']:
        gff = cmdResults['GFF']
    else:
        gff = "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208_noFasta.gff"
    
    cwd      = os.getcwd() 
    ftype    = ''            # sam or bam
    LOG_FILE =  'htseq.log'    
    
    if os.path.exists(LOG_FILE):
        os.unlink(LOG_FILE)
        logging.basicConfig(filename=LOG_FILE, level=logging.INFO,)
        logging.info(' DELETED previous log file')
    else:
        logging.basicConfig(filename=LOG_FILE, level=logging.INFO,)
    
    logging.info(' Running HTSeq-count')
    logging.info(' Directory %s' %(cwd))
    logging.info(' Using GFF file : %s' %(gff))
        
    if cmdResults['SAM']:
        ftype = 'sam'
        files = [ fn for fn in os.listdir(cwd) if fn.endswith(".sam") ]       # Get a list of all sam files in current directory
        logging.info(' Using SAM files')
    else:
        ftype = 'bam'
        files = [ fn for fn in os.listdir(cwd) if fn.endswith(".bam") ]    # Get a list of all bam files in current directory
        logging.info(' Using BAM files')

    if cmdResults['REVERSE']:
        logging.info(' using htseq-count with -s option, for anti-sense library.')
        rvse = 1
        runHTSeq( files, ftype, gff,  rvse )
        
    else:
        runHTSeq  ( files, ftype, gff )
        
if __name__ == "__main__":
    main()

