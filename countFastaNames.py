#!/home/mplace/anaconda3/bin/python
"""
Name: countFastaNames.py 
Purpose: Get Gene and Seq Count info by comparing the sequence names of 2 fasta files 
Input:  2 fasta files
Output: Print to screen 
Author: Mike Place
Date: 1/29/2015
"""
from Bio import SeqIO
import argparse
import re	
import textwrap
    

def main():
    """
    Start Here
    """    
    # handle command line arguments
    cmdparser = argparse.ArgumentParser(description="Compare 2 fasta files & produce a unique file", prog='compareFasta.py'  )
    cmdparser.add_argument('-f1','--file1', action='store', required='true',
                           dest='FILE1', help='REQUIRED, fasta file (.fasta, .fa , .fsa)')
    cmdparser.add_argument('-f2','--file2', action='store', required='true',
                           dest='FILE2', help='REQUIRED, fasta file (.fasta, .fa , .fsa)')
    cmdResults = vars(cmdparser.parse_args())

    firSeq   = {}
    firName  = {}
    secSeq   = {}    
    secName  = {}
    result   = {}
    
    if cmdResults['FILE1'] is not None:
        f1 = cmdResults['FILE1']
        f2 = cmdResults['FILE2']
        fname = f1.split('.')
        sname = f2.split('.')

       
    # parse the first fasta file, get unique set
    for seq1 in SeqIO.parse(f1,"fasta"): 
        tseq    = str(seq1.seq)
        name    = str(seq1.id)
        if name in firName:
            firName[name] = firName[name] + 1        
        else:
            firName[name] = 1
        
        if tseq in firSeq:      # if sequence is found, concatenate gene name
            firSeq[tseq] = firSeq[tseq] + "|" + seq1.id + "_" + fname[0]
        else:
            firSeq[tseq] = seq1.id + "_" + fname[0]
        
            
    # parse the second fasta file, get unique set
    for seq2 in SeqIO.parse(f2,"fasta"):
        nseq = str(seq2.seq)
        name = str(seq2.id)
        if name in secName:
            secName[name] = secName[name] + 1        
        else:
            secName[name] = 1
        
        if nseq in secSeq:
            secSeq[nseq] = secSeq[nseq] + "|" + seq2.id + "_" + sname[0]
        else:
            secSeq[nseq] = seq2.id + "_" + sname[0]
            
    #Search secSeq dict with firSeq keys
    for k,v in firSeq.items():
        if k in secSeq:
            result[k] = v + "|" + secSeq[k] 
        else:
            result[k] = v 
    
    #Now do reciprocal search second dict
    for k,v in secSeq.items():
        if k in firSeq:
            result[k] = v + "|" + firSeq[k]            
        else:
            result[k] = v
            
    bothNameCount = 0
    #match names
    for k,v in firName.items():
        if k in secName:
            bothNameCount += 1  

    firNameLen = len(firName)
    firSeqLen  = len(firSeq)
    secNameLen = len(secName)
    secSeqLen  = len(secSeq)
    print("\n%s num genes: %d Unique sequences: %d" %(f1, firNameLen, firSeqLen) )
    print("\n%s num genes: %d Unique sequences: %d" %(f2, secNameLen, secSeqLen) ) 
    
    print("Sequence count by group:")
    
    groups = { 'liftover' : 0, 'ygap' : 0, 'both' : 0 }    
    
    for k,v in result.items():
        #print( "%s" %(v))
        if re.search('\|',v):
            groups['both'] += 1
        elif v.endswith('ygap'):
            groups['ygap'] += 1
        else:
            groups['liftover'] += 1
    
    for k,v in groups.items():
        print( "%s : %s" %(k,v))
    
    resultLen = len(result)
    print("\nTotal number of genes found: %d" %(resultLen))  

    print("Gene names in common (exact match): %d" %(bothNameCount))      

if __name__ == "__main__":
    main()