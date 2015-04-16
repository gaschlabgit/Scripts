#!/home/mplace/anaconda3/bin/python
"""
uniqueFasta.py
Print out unique sequences within a fasta file.
Uses  Biopython 
The sequence is used as a key to a dictionary.

"""
from Bio import SeqIO
import argparse	
import textwrap

    

def main():
    """
    Start Here
    """    
    # handle command line arguments
    cmdparser = argparse.ArgumentParser(description="Produce a unique fasta sequence file", prog='uniqueFasta.py'  )
    cmdparser.add_argument('-f1','--file1', action='store', required='true',
                           dest='FILE1', help='REQUIRED, fasta file (.fasta, .fa , .fsa)')        
   
    cmdResults = vars(cmdparser.parse_args())
    firSeq   = {}
    
    if cmdResults['FILE1'] is not None:
        f1 = cmdResults['FILE1']

    
    for seq1 in SeqIO.parse(f1,"fasta"):
        tseq = str(seq1.seq)

        if tseq in firSeq:
            firSeq[ tseq ] = firSeq[ tseq ] + "|" + seq1.id
        else:
            firSeq[ tseq ] = seq1.id
            
    #{ print (">%s\n%s" %(firSeq[k], textwrap.wrap(k,width=50) ) for k in firSeq.keys()}
    for k,v in firSeq.items():
        print (">%s" %(v))
        data = (textwrap.wrap(k,width=80) )
        [ print(i) for i in data ]
        
        
        


if __name__ == "__main__":
    main()