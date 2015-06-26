#!/home/mplace/anaconda3/bin/python
"""
@Program: convertChromName.py

@Purpose: Convert chromosome names in fasta file.
          Change ref|NC_.... to chr1
          or chrI to ref|NC_.... (must modify code see below)
          
@Input:  Fasta file

@Dependencies: Python 3

@author: Mike Place


Example:


"""
import re
import sys

chromSGD = { "chrI" : "ref|NC_001133|", "chrII" : "ref|NC_001134|", "chrIII" : "ref|NC_001135|",
              "chrIV" : "ref|NC_001136|", "chrV" : "ref|NC_001137|", "chrVI" : "ref|NC_001138|", 
              "chrVII" : "ref|NC_001139|", "chrVIII" : "ref|NC_001140|", "chrIX" : "ref|NC_001141|",
              "chrX" : "ref|NC_001142|", "chrXI" : "ref|NC_001143|", "chrXII" : "ref|NC_001144|",
              "chrXIII" : "ref|NC_001145|", "chrXIV" : "ref|NC_001146|", "chrXV" : "ref|NC_001147|",
              "chrXVI" : "ref|NC_001148|", "chrMito" : "ref|NC_001224|", "chrmt" : "ref|NC_001224|"
} 

chromArabic = {  "ref|NC_001133|" : "chr1", "ref|NC_001134|" : "chr2", "ref|NC_001135|" : "chr3",
           "ref|NC_001136|" : "chr4", "ref|NC_001137|" : "chr5", "ref|NC_001138|" : "chr6", 
           "ref|NC_001139|" : "chr7", "ref|NC_001140|" : "chr8", "ref|NC_001141|" : "chr9",
           "ref|NC_001142|" : "chr10", "ref|NC_001143|" : "chr11", "ref|NC_001144|" : "chr12",
           "ref|NC_001145|" : "chr13", "ref|NC_001146|" : "chr14", "ref|NC_001147|" : "chr15",
           "ref|NC_001148|" : "chr16", "ref|NC_001224|" : "chrMito" } 
           
chromRoman = {  "ref|NC_001133|" : "chrI", "ref|NC_001134|" : "chrII", "ref|NC_001135|" : "chrIII",
           "ref|NC_001136|" : "chrIV", "ref|NC_001137|" : "chrV", "ref|NC_001138|" : "chrVI", 
           "ref|NC_001139|" : "chrVII", "ref|NC_001140|" : "chrVIII", "ref|NC_001141|" : "chrIX",
           "ref|NC_001142|" : "chrX", "ref|NC_001143|" : "chrXI", "ref|NC_001144|" : "chrXII",
           "ref|NC_001145|" : "chrXIII", "ref|NC_001146|" : "chrXIV", "ref|NC_001147|" : "chrXV",
           "ref|NC_001148|" : "chrXVI", "ref|NC_001224|" : "chrMito" } 
           
chromDict64 = { "ref|NC_001133|[R64]" : "chr1",  "ref|NC_001134|[R64]": "chr2", "ref|NC_001135|[R64]" : "chr3",
              "ref|NC_001136|[R64]" : "chr4", "ref|NC_001137|[R64]" : "chr5", "ref|NC_001138|[R64]" : "chr6", 
              "ref|NC_001139|[R64]" : "chr7", "ref|NC_001140|[R64]" : "chr8", "ref|NC_001141|[R64]" : "chr9",
              "ref|NC_001142|[R64]" : "chr10", "ref|NC_001143|[R64]" : "chr11", "ref|NC_001144|[R64]" : "chr12",
              "ref|NC_001145|[R64]" : "chr13", "ref|NC_001146|[R64]" : "chr14", "ref|NC_001147|[R64]" : "chr15",
              "ref|NC_001148|[R64]" : "chr16", "ref|NC_001224|[R64]" : "chr17" } 


def main():
    """
    main() 
    """
  
    # if no args print help
    if len(sys.argv) == 1:
        print("\n\tInput file required.")
        print("\n")
        print("\tconvertChromNameFasta.py <fasta file>")
        print("\tConverts ref|NC_001133| name in .fasta file to chr1 etc...")
        sys.exit(1)
    else:
        fasta = sys.argv[1]
        
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                line = line.rstrip()
                line = line.lstrip('>')
                row  = line.split()
                if row[0] in chromArabic:
                    print(">%s" %( chromArabic[row[0]] ) )
            else:
                line = line.rstrip()
                print(line),
               

if __name__ == "__main__":
    main()

