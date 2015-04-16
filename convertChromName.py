#!/home/mplace/anaconda3/bin/python
"""
@Program: convertChromName.py

@Purpose: Convert chromosome names in gff file.
          Change chrI  to ref|NC_
          
@Input:  Fastq file, unzipped

@Dependencies: Python 3

@author: Mike Place


Example:


"""
import re
import sys

chromDict = { "chrI" : "ref|NC_001133|", "chrII" : "ref|NC_001134|", "chrIII" : "ref|NC_001135|",
              "chrIV" : "ref|NC_001136|", "chrV" : "ref|NC_001137|", "chrVI" : "ref|NC_001138|", 
              "chrVII" : "ref|NC_001139|", "chrVIII" : "ref|NC_001140|", "chrIX" : "ref|NC_001141|",
              "chrX" : "ref|NC_001142|", "chrXI" : "ref|NC_001143|", "chrXII" : "ref|NC_001144|",
              "chrXIII" : "ref|NC_001145|", "chrXIV" : "ref|NC_001146|", "chrXV" : "ref|NC_001147|",
              "chrXVI" : "ref|NC_001148|", "chrMito" : "ref|NC_001224|", "chrmt" : "ref|NC_001224|"
} 

def main():
    """
    main() 
    """
  
    # if no args print help
    if len(sys.argv) == 1:
        print("\n\tInput file required.")
        print("\n")
        print("\tconvertChromName.py <gff file>")
        print("\tConverts chrI name in .gff file to ref|NC_001133| to match")
        print("\tthe SGD reference fasta file naming conventions\n")
        sys.exit(1)
    else:
        fastq = sys.argv[1]
        
    with open(fastq) as f:
        for line in f:
            if line.startswith('#'):
                line = line.rstrip()
                print(line),
                continue
            else:
                line  = line.rstrip()
                row   = line.split('\t')
                chrom = row[0]
                if chrom in chromDict:
                    row[0]= chromDict[chrom]
                    print( "\t".join(row) )
                else:
                    print( "\t".join(row) )

if __name__ == "__main__":
    main()

