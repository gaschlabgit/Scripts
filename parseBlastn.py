#!/home/mplace/anaconda/bin/python
"""
Created on Fri Jul 24 15:06:07 2015

@Program: parseBlastn.py

@Purpose: Create a table from blastn results text file.

@Input:  text file, blastn results run on command line

@Output: table of results
@author: Mike Place
@Date:   7/24/2015
@Dependencies:  python 2.7

Bio.Blast.Record.HSP

Stores information about one hsp in an alignment hit.

Members:

        score BLAST score of hit. (float)
        bits Number of bits for that score. (float)
        expect Expect value. (float)
        num_alignments Number of alignments for same subject. (int)
        identities Number of identities (int) if using the XML parser. Tuple of numer of identities/total aligned (int, int) if using the (obsolete) plain text parser.
        positives Number of positives (int) if using the XML parser. Tuple of numer of positives/total aligned (int, int) if using the (obsolete) plain text parser.
        gaps Number of gaps (int) if using the XML parser. Tuple of numer of gaps/total aligned (int, int) if using the (obsolete) plain text parser.
        align_length Length of the alignment. (int)
        strand Tuple of (query, target) strand.
        frame Tuple of 1 or 2 frame shifts, depending on the flavor.
        query The query sequence.
        query_start The start residue for the query sequence. (1-based)
        query_end The end residue for the query sequence. (1-based)
        match The match sequence.
        sbjct The sbjct sequence.
        sbjct_start The start residue for the sbjct sequence. (1-based)
        sbjct_end The end residue for the sbjct sequence. (1-based)


"""
import argparse
import sys
from Bio.Blast import NCBIXML

def main():
    """
    Main
    """    
#******************************************************************************
# Command line args
#******************************************************************************    
    cmdparser = argparse.ArgumentParser(description="Create table from a list of blastn text files.",
                                        usage='%(prog)s -f <file> list of blastn file to process.' ,prog='parseBlastn.py'  )
    cmdparser.add_argument('-f', '--file', action='store', dest='FILE',   help='File listing blastn files to process, one per line.')
    cmdResults = vars(cmdparser.parse_args())
    
    # variables
    files = []    
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
        
    if cmdResults['FILE']:
        inFile = cmdResults['FILE']

    with open(inFile,'r') as data:
        for item in data:
            files.append(item.rstrip())
    
    for xmlFile in files:
        for record in NCBIXML.parse( open(xmlFile) ):
            if record.alignments:
                print "QUERY: %s" % record.query[:60]
                for align in record.alignments:
                    for hsp in align.hsps:
                        print "bits %s" % hsp.bits 
                        print hsp.expect
                        print hsp.score
                        print hsp.identities
                        print hsp.gaps
        
        
            
    
        

if __name__ == "__main__":
    main()