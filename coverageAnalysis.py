#!/home/mplace/anaconda3/bin/python
"""
Created May 29th 2015

@author: mplace

Program: coverageAnalysis.py 

Purpose: Use samtools to count reads in 1kb window, stepping 100bp.

Input  : bam file
 
Output: Table where columns are position, count, %GC

dependencies: samtools, must be in path

To visual with R:

library(ggplot2)

d <- read.table("chr1.copynumber", header=T)

data looks like:
 pos   ct       gc
1 1010  647 0.424419
2 2022 1523 0.354627
3 3001 1449 0.398337
4 4000  686 0.337660

qplot(d$pos,d$ct )



"""
import argparse 
import os
import sys
import re
import subprocess 
from collections import defaultdict     

chromSet = {"ref|NC_001133|", "ref|NC_001134|", "ref|NC_001135|", "ref|NC_001136|", 
             "ref|NC_001137|", "ref|NC_001138|", "ref|NC_001139|", "ref|NC_001140|",
             "ref|NC_001141|", "ref|NC_001142|", "ref|NC_001143|", "ref|NC_001144|",
             "ref|NC_001145|", "ref|NC_001146|", "ref|NC_001147|", "ref|NC_001148|", 
             "ref|NC_001224|", "ref|NC_001224|"
} 

chromLen = { "ref|NC_001133|" : 230218, "ref|NC_001134|" : 813184, "ref|NC_001135|" : 316620,
             "ref|NC_001136|" : 1531933, "ref|NC_001137|" : 576874, "ref|NC_001138|" : 270161,
             "ref|NC_001139|" : 1090940, "ref|NC_001140|" : 562643, "ref|NC_001141|" : 439888,
             "ref|NC_001142|" : 745751, "ref|NC_001143|" : 666816, "ref|NC_001144|" : 1078177,
             "ref|NC_001145|" : 924431, "ref|NC_001146|" : 784333, "ref|NC_001147|" : 1091291,
             "ref|NC_001148|" : 948066, "ref|NC_001224|" : 85779 
}

# key = chromosome  each value is a list of reads counts by window
chromDict = defaultdict(list)

def main():
    """
    Read in bam file, parse and count reads for each 1kb non-overlaping bins.
    """
    cmdparser = argparse.ArgumentParser(description="Calculate Copy Number Variants using a sorted bam file.",
                                        usage='%(prog)s <file.bam>' ,prog='coverageAnalysis.py'  )
    cmdparser.add_argument('-b', '--bam',  action='store',      dest='BAM',  help='bam file')         
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='Program information')
    cmdResults = vars(cmdparser.parse_args())
    
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)    
    
    if cmdResults['INFO']:
        print("\n coverageAnalysis.py")
        print("\n Purpose: Copy Number Analysis.")
        print("\t  parse and count reads for each 1kb non-overlaping bins.")
        print(" Required Parameters: -b sampleName.bam")
        print("\n Input : Bam file")
        print("\n\tAssumes Chrom names are like: ref|NC_001133|\n")
        print(" To run:  /home/GLBRCORG/mplace/scripts/coverageAnalysis.py -b samplename.bam\n")
        print(" Output: Tab delimited table w/ columns WindowStartPosition WindowCount %GC ")
        print("         Start of new chromosome marked by chromosome name. ")
        print(" Example:")        
        print("\tref|NC_001133|" )
        print("\t1010    647     0.424419")
        print("\t2022    1523    0.354627")
        print("\n")  
        print(" To visual with R:")
        print(" Edit file: add header and split by chromosome")
        print("\n library(ggplot2)")
        print(" d <- read.table( inputfile, header=T)")
        print(" qplot(d$pos,d$ct )" )
        print("")

        sys.exit(1)
    
    if cmdResults['BAM']:
        infile = cmdResults['BAM']   
        samName = re.sub(r"bam", "sam", infile )        
        
        if ( not os.path.exists(samName) ):
            cmd  = [ 'samtools', 'view', '-o', samName, infile ]
            output  = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
            
        first     = 0
        chrName   = ""
        startPos  = 0
        endPos    = 1000
        ctr       = 0
        gcContent = 0
        totLen    = 0
        with open(samName, 'r') as reads:
            for i in reads:
                line = i.split()
                if first == 0:
                    chrName   = line[2]
                    first     = 1
                    ctr       += 1
                    gcContent += len(re.findall("G|C",line[9]))
                    totLen    += len(line[9])
                    print("%s" %(chrName) )
                elif chrName != line[2]:
                    chrName   = line[2]
                    startPos  = 0
                    endPos    = 1000
                    ctr       = 1
                    gcContent += len(re.findall("G|C",line[9]))
                    totLen    += len(line[9])
                    print("%s" %(chrName) )
                elif ( int(line[3]) < endPos and int(line[3]) >= startPos ):
                    ctr += 1
                    gcContent += len(re.findall("G|C",line[9]))
                    totLen    += len(line[9])
                else:
                    ctr += 1
                    gcContent += len(re.findall("G|C",line[9]))
                    totLen    += len(line[9])    
                    percentage = gcContent/totLen
                    print( "%s\t%d\t%f" %(line[3], ctr, percentage))                                
                    chromDict[chrName].append(ctr)
                    startPos += 1000
                    endPos   += 1000
                    ctr       = 1
                    gcContent = 0
                    totLen    = 0
          
if __name__ == "__main__":
    main()

