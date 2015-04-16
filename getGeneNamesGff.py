#!/home/mplace/anaconda3/bin/python
"""
@Program: getGeneNamesGff.py

@Purpose: Create a gff file from Sean's annotation files

@Input:  file with list of gene names, Sean's annotation
         only does one type at at time, only LIFTOVER or YGAP
         so run them separately and cat later.

@Output: GFF file 
@author: Mike Place
@Date:   2/19/2015

"""
import sys
import re

def main():
    """
    main() 
    """
    # if no args print help
    if len(sys.argv) == 1:
        print("file required")
        sys.exit(1)
    else:
        file    = sys.argv[1]   # list of gene names to parse
        annFile = sys.argv[2]   # annotation file, csv format
        names   = set()
        
    with open(file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('Homologue'):
                #
                #nl = line.replace("_ygap", "")              
                names.add(line)
            else:
                gname = re.split("_", line)
                names.add(gname[0])

    with open("homoCheck.txt",'w') as fo:
        for x in names:
            fo.write("%s\n" %(x))
    
    with open(annFile) as a:
        next(a)
        for line in a:
            line     = line.rstrip()
            first    = line.split(',')            
            gene = first[9]
            #first[9] = gene + "_liftover"
            first[9] = gene + "_ygap"            
            att = first[8].split(';')
            att[1] = gene + "_ygap"
            #att[1] = gene + "_liftover"
            first[8] = ';'.join(att)
            if gene in names:
                row = first[:10]
                print('\t'.join(row))              
             

if __name__ == "__main__":
    main()

