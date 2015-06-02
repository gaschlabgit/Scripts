#!/home/mplace/anaconda3/bin/python
"""
@Program: fixFastaName.py

@Purpose: Reformat fasta sequence headers to match bed file.
          Write a new fasta file using strain name in chromosome header.

@Input: fasta file downloaded from ncbi
        input fasta file header looks like:
            CP005027.1 Saccharomyces cerevisiae YJM1615 chromosome IX sequence

@Dependencies: Python 3
               BioPython      
               
@Output: fasta file
@author: Mike Place
@Date:   5/27/2015
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse	
import os
import re
import sys

chromDict = { "I" : "Chr01", "II" : "Chr02", "III" : "Chr03",
              "IV" : "Chr04", "V" : "Chr05", "VI" : "Chr06", 
              "VII" : "Chr07", "VIII" : "Chr08", "IX" : "Chr09",
              "X" : "Chr10", "XI" : "Chr11", "XII" : "Chr12",
              "XIII" : "Chr13", "XIV" : "Chr14", "XV" : "Chr15",
              "XVI" : "Chr16", "complete" : "mito", "2" : "plasmid2micron",
              "2micron" : "plasmid2micron"
}

def writeRefFasta( rec, geneName, chrom, left, right ):
    """
    Write Gene Reference sequence to file.
    """
    start = int(left) + 1
    left = str(start)
    file = geneName + ".reference.fasta"
    geneRec = SeqRecord(rec)
    geneRec.id = geneName + "|chr" + chrom + ":" + str(left) + "-" + str(right)
    geneRec.name = geneName
    geneRec.description = "Seq extracted by convertGene.py"
    SeqIO.write(geneRec, file, "fasta")
    return file 

def main():
    """
    Main, parse cmd line args and call writeRefFasta. 
    """
#******************************************************************************
# Command line args
#******************************************************************************    
    cmdparser = argparse.ArgumentParser(description="Create a new fasta file, renaming each chrom.",
                                        usage='%(prog)s -f <filelist> ' ,prog='extractGene.py'  )
    cmdparser.add_argument('-f', '--file',   action='store', dest="FILE",  help='Required: fasta input file (.fasta, .fsa, .fa )')
    cmdparser.add_argument('-d', '--detail', action='store_true', dest="DETAIL" )
             
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)

#******************************************************************************
# Check and define input variables
# Print each parameter used to screen to give the user a record.
#******************************************************************************    
    # Print detailed program information to screen
    if cmdResults['DETAIL']:
        print("\n\tfixFastaName.py program information:\n")
        print("\tPurpose: Given a list of Complete Genome Fasta files, rename each")
        print("\t         chromosome sequence and write a new file.\n")
        print("\tSequence name is expected to look like: \n")
        print("\tCP005630.1 Saccharomyces cerevisiae YJM1244 chromosome XV sequence\n")
        print("\tThe strain name YJM1244 is used to name the new fasta file.")
        print("\tThe chrom roman numeral is mapped to a Chr15 to name the sequence.")
        print("\n  Required parameters:")
        print("\t-f List of fasta files, a text file, 1 file name per line.")
        print("\n\toutput file:")
        print("\t 1) Strain_Name.fasta - complete fasta sequence\n\n")
        sys.exit(1)
    
    print("\n\tRunning fixFastaName.py ....")
    
    if cmdResults['FILE'] is not None:
        inFile = cmdResults['FILE']
        print("\tFile     : %s" %(inFile))  
    
    fastaFiles = []
        
    # read in the fasta file list
    with open( inFile, 'r') as f:
        for line in f:
            line = line.rstrip()
            fastaFiles.append(line)
    
    for i in fastaFiles:
        for seq in SeqIO.parse(str(i), "fasta"):
            item       = seq.description
            header     = item.split()
            strainName = header[3]
            seq.id = chromDict[header[5]]
            seq.description = chromDict[header[5]]
            #print(seq.format("fasta"))
            out = open( strainName + ".fasta", "a")
            SeqIO.write(seq, out, "fasta")
            out.close()

if __name__ == "__main__":
    main()

