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

chromDictRoman = { "ref|NC_001133|" : "chrI",  "ref|NC_001134|": "chrII", "ref|NC_001135|" : "chrIII",
              "ref|NC_001136|" : "chrIV", "ref|NC_001137|" : "chrV", "ref|NC_001138|" : "chrVI", 
              "ref|NC_001139|" : "chrVII", "ref|NC_001140|" : "chrVIII", "ref|NC_001141|" : "chrIX",
              "ref|NC_001142|" : "chrX", "ref|NC_001143|" : "chrXI", "ref|NC_001144|" : "chrXII",
              "ref|NC_001145|" : "chrXIII", "ref|NC_001146|" : "chrXIV", "ref|NC_001147|" : "chrXV",
              "ref|NC_001148|" : "chrXVI", "ref|NC_001224|" : "chrMito", "ref|NC_001224|" : "chrmt"
} 

chromDict64 = { "ref|NC_001133|[R64]" : "chr1",  "ref|NC_001134|[R64]": "chr2", "ref|NC_001135|[R64]" : "chr3",
              "ref|NC_001136|[R64]" : "chr4", "ref|NC_001137|[R64]" : "chr5", "ref|NC_001138|[R64]" : "chr6", 
              "ref|NC_001139|[R64]" : "chr7", "ref|NC_001140|[R64]" : "chr8", "ref|NC_001141|[R64]" : "chr9",
              "ref|NC_001142|[R64]" : "chr10", "ref|NC_001143|[R64]" : "chr11", "ref|NC_001144|[R64]" : "chr12",
              "ref|NC_001145|[R64]" : "chr13", "ref|NC_001146|[R64]" : "chr14", "ref|NC_001147|[R64]" : "chr15",
              "ref|NC_001148|[R64]" : "chr16", "ref|NC_001224|[R64]" : "chr17" } 

chromDict = { "ref|NC_001133|" : "chr1",  "ref|NC_001134|" : "chr2",  "ref|NC_001135|" : "chr3",
              "ref|NC_001136|" : "chr4",  "ref|NC_001137|" : "chr5",  "ref|NC_001138|" : "chr6", 
              "ref|NC_001139|" : "chr7",  "ref|NC_001140|" : "chr8",  "ref|NC_001141|" : "chr9",
              "ref|NC_001142|" : "chr10", "ref|NC_001143|" : "chr11", "ref|NC_001144|" : "chr12",
              "ref|NC_001145|" : "chr13", "ref|NC_001146|" : "chr14", "ref|NC_001147|" : "chr15",
              "ref|NC_001148|" : "chr16", "ref|NC_001224|" : "chr17" } 
              
chrom  = { "chr1" : "chrI",  "chr2": "chrII", "chr3" : "chrIII",
              "chr4" : "chrIV", "chr5" : "chrV", "chr6" : "chrVI", 
              "chr7" : "chrVII", "chr8" : "chrVIII", "chr9" : "chrIX",
              "chr10" : "chrX", "chr11" : "chrXI", "chr12" : "chrXII",
              "chr13" : "chrXIII", "chr14" : "chrXIV", "chr15" : "chrXV",
              "chr16" : "chrXVI", "chr17" : "chrMito" }


def main():
    """
    main() 
    """
  
    # if no args print help
    if len(sys.argv) == 1:
        print("\n\tInput file required.")
        print("\n")
        print("\tconvertChromName.py <sam file>")
        print("\tConverts ref|NC_001133| in sam file to match chr1 , etc..")
        print("\n")
        sys.exit(1)
    else:
        fastq = sys.argv[1]        
        
    with open(fastq) as f:
        for line in f:
            if line.startswith('@'):
                line = line.rstrip()
                header = line.split('\t')
                name = header[1].split(':')
                if name[1] in chrom:
                    name[1] = chrom[name[1]]
                    out = header[0] + "\t" + name[0] + ":" + name[1] + "\t" + header[2]
                    print(out)
                else:
                    print(line)
                continue
            else:
                line  = line.rstrip()
                row   = line.split('\t')
                if row[2] in chrom:
                    row[2]= chrom[row[2]]
                    print( "\t".join(row) )

if __name__ == "__main__":
    main()

