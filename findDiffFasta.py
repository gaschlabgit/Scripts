#!/home/mplace/anaconda3/bin/python
"""
Name: findDiffFasta.py 
Purpose: compare the sequences of 2 fasta files, output sequences in the first
         file but not in the second
Input:  2 fasta files
Output: single fasta file, sequences in first but not the second 
Author: Mike Place
Date: 1/29/2015
"""
from Bio import SeqIO
import argparse	
import textwrap
    

def main():
    """
    Start Here
    """    
    # handle command line arguments
    cmdparser = argparse.ArgumentParser(description="Compare 2 fasta files & output seq in first not in second file", prog='compareFasta.py'  )
    cmdparser.add_argument('-f1','--file1', action='store', required='true',
                           dest='FILE1', help='REQUIRED, fasta file (.fasta, .fa , .fsa)')
    cmdparser.add_argument('-f2','--file2', action='store', required='true',
                           dest='FILE2', help='REQUIRED, fasta file (.fasta, .fa , .fsa)')
    cmdparser.add_argument('-o', '--output', action='store', required='false', dest='OUT', help='Diff.fasta')         
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
    
    if cmdResults['OUT'] is not None:
        outFile = cmdResults['OUT']
    else:
        outFile = "Unique.fasta"
       
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
            continue #result[k] = v + "|" + secSeq[k] 
        else:
            result[k] = v 
    
        
    # print output in fasta format     
    fileHandle = open(outFile, 'w')
    
    for k,v in result.items():
        fileHandle.write(">%s\n" %(v))
        #print (">%s" %(v))
        data = (textwrap.wrap(k,width=80) )
        #[ print(i) for i in data ]        
        for i in data:
            fileHandle.write("%s\n" %(i))
            

if __name__ == "__main__":
    main()