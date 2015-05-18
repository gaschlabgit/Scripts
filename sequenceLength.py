#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
program: sequenceLength.py

purpose: print out seqence name and length of sequence
         for a fasta file.

usage: sequenceLength.py -f input.fasta

input: a fasta file 

output: text to stdout, 
        seqName,Length
        seq1,129383

By Mike Place

"""
from Bio import SeqIO
import argparse


def main():
    """
    Process command line arguments
    """
    cmdparser = argparse.ArgumentParser( description="Count sequence lengths in a fasta file",
                                        prog='sequenceLength.py')
    cmdparser.add_argument('-f', '--file', action='store', required='true',
                           dest='FILE', help='Required: fasta input file (.fasta, .fsa, .fa )')
    
    cmdResults = vars( cmdparser.parse_args() )
    
    if cmdResults['FILE'] is not None:
        inFile = cmdResults['FILE']
    
    for seqRec in SeqIO.parse(str(inFile), "fasta"):
        print( "%s,%i" %(seqRec.id, len(seqRec)))
    
    


if __name__ == "__main__":
    main()
    
