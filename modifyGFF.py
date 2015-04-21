#!/home/mplace/anaconda3/bin/python
"""
@Program: modifyGFF.py

@Purpose: Add a 3'UTR feature to the GFF file for every gene, starting at the base 
following the stop codon and ending at the base before the  start codon of the next gene.
          
@Input:  gene_only.gff

gawk ' ( $3 ~ /^gene$/ ) { print $0} ' saccharomyces_cerevisiae_R64-1-1_20110208_noFasta.gff  > R64-1-1_genes_only.gff

@Dependencies: Python 3

@Output:  
@author: Mike Place
@Date:   4/21/2015

Example:

~/scripts/modifyGFF.py original.gff

"""
from itertools import tee, islice, chain, zip_longest
import re
import sys
    
def peek(some_iterable):
    """
    return a pair of iterable objects
    """
    current, nexts = tee(some_iterable,2)         # provides 2 pointers to the list
    nexts = chain(islice(nexts, 1, None), [None]) # creates a list right shifted, with last item = None
    return zip_longest(current, nexts)
   
def main():
    """ 
    """
    # if no args print help
    if len(sys.argv) < 1:
        print("Input file required.")
        sys.exit(1)
    else:
        inFile = sys.argv[1]
    
    gene = []
    
    with open( inFile, 'r') as f:
            for line in f:
                line = line.rstrip()
                data = line.split("\t")
                if len(data) >= 9:
                    if re.match( r'^gene$', data[2] ):
                        gene.append(data)   
                
    for current, nxt in peek(gene):
        if nxt == None:
            break        
        if current[0] != nxt[0]:
            continue          
        if current[4] >= nxt[3] and current[4] <= nxt[4]:
            continue
        if current[4] > nxt[4]:
            continue
        nStart = int(current[4]) + 1
        nEnd   = int(nxt[3]) -1
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" %(current[0], current[1], "3'_next_gene", nStart, nEnd, 
                   current[5], current[6], current[7], current[8]))
        #print(current[0])
        #print( "next ", end=" ")
        #print(nxt[3])
        
  

    

if __name__ == "__main__":
    main()

