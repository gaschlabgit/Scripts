#!/home/mplace/anaconda3/bin/python
"""
@Program: convertChromName_SamFile.py

@Purpose: Convert chromosome names sam file.
          Change ref|NC_XXXXX to chrI, etc....
          
@Input:  Fastq file, unzipped

@Dependencies: Python 3

@author: Mike Place


Example:


"""
import re
import sys

chromDict = { "ref|NC_001133|" : "chrI",  "ref|NC_001134|": "chrII", "ref|NC_001135|" : "chrIII",
              "ref|NC_001136|" : "chrIV", "ref|NC_001137|" : "chrV", "ref|NC_001138|" : "chrVI", 
              "ref|NC_001139|" : "chrVII", "ref|NC_001140|" : "chrVIII", "ref|NC_001141|" : "chrIX",
              "ref|NC_001142|" : "chrX", "ref|NC_001143|" : "chrXI", "ref|NC_001144|" : "chrXII",
              "ref|NC_001145|" : "chrXIII", "ref|NC_001146|" : "chrXIV", "ref|NC_001147|" : "chrXV",
              "ref|NC_001148|" : "chrXVI", "ref|NC_001224|" : "chrMito", "ref|NC_001224|" : "chrmt"
} 


def main():
    """
    main() 
    """
  
    # if no args print help
    if len(sys.argv) == 1:
        print("\n\tInput file required.")
        print("\n")
        print("\tconvertChromName.py <sam file>")
        print("\tConverts ref|NC_001133| in sam file to match chrI , etc..")
        print("\n")
        sys.exit(1)
    else:
        fastq = sys.argv[1]        
        
    with open(fastq) as f:
        for line in f:
            if line.startswith('@'):
                line = line.rstrip()
                print(line),
                continue
            else:
                line  = line.rstrip()
                row   = line.split('\t')
                if row[2] in chromDict:
                    row[2]= chromDict[row[2]]
                    print( "\t".join(row) )

if __name__ == "__main__":
    main()

