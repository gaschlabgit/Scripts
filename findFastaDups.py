#!/home/mplace/anaconda3/bin/python
"""
findFastaDups.py
Print out duplicate sequences within one fasta file.
Uses  Biopython 
The sequence is used as a key to a dictionary.
Duplicate keys are printed out.

"""
from Bio import SeqIO
import argparse	


#for liftover in SeqIO.parse("liftover.allpass.aa.2015_01_20.fasta", "fasta"):
#    print(liftover.id)
#    print(repr(liftover.seq))
#    print(len(liftover))
    

def main():
    """
    Start Here
    """    
    # handle command line arguments
    cmdparser = argparse.ArgumentParser(description="Compare fasta sequences in a file & produce a duplicate list", prog='findFastaDups.py'  )
    cmdparser.add_argument('-f1','--file1', action='store', required='true',
                           dest='FILE1', help='REQUIRED, fasta file (.fasta, .fa , .fsa)')        
   
    cmdResults = vars(cmdparser.parse_args())
    firSeq   = {}
    
    if cmdResults['FILE1'] is not None:
        f1 = cmdResults['FILE1']

    
    for seq1 in SeqIO.parse(f1,"fasta"):
        tseq = str(seq1.seq)

        if tseq in firSeq:
            print( "%s,%s" %( seq1.id, seq1.seq ) )
            print( "Dup:%s,%s" %( firSeq[ tseq ], tseq) )
        else:
            firSeq[ tseq ] = seq1.id
        


if __name__ == "__main__":
    main()