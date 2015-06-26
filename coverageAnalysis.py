#!/home/mplace/anaconda3/bin/python
"""
Program: coverageAnalysis.py 

Purpose: Use samtools to count reads in a user specified window size.
         Default window size is 500 bp.

Input  : bam file
 
Output : Table where columns are chrom, start_position, end_position, count, count/median_count, log2(count/median_count)

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

Created June 26th 2015

@author: mplace

"""
import argparse 
import math
import os
import sys
import re
import subprocess 
from collections import defaultdict  
from statistics import median

# key = chromosome  each value is a dict of lists, the keys are 'pos', 'cnt'  
chrom    = defaultdict(dict)


def main():
    """
    Read in bam file, parse and count reads for each non-overlaping bin.
    """
    cmdparser = argparse.ArgumentParser(description="Calculate Copy Number Variants using a sorted bam file.",
                                        usage='%(prog)s <file.bam>' ,prog='coverageAnalysis.py'  )
    cmdparser.add_argument('-b', '--bam',  action='store',      dest='BAM',  help='bam file', metavar='') 
    cmdparser.add_argument('-w', '--win',  action='store',      dest='WIN',  help='Window size for counting reads', metavar='')        
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='Program information')
    cmdResults = vars(cmdparser.parse_args())
    
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)    
    
    if cmdResults['INFO']:
        print("\n coverageAnalysis.py")
        print("\n Purpose: Copy Number Analysis.")
        print("\t  parse and count reads for non-overlaping bins.")
        print(" Required Parameters: -b sampleName.bam")
        print(" Optional Parameters: -w window size")
        print("\n Input : Bam file")
        print(" To run:  /home/GLBRCORG/mplace/scripts/coverageAnalysis.py -b samplename.bam \n")
        print(" Output: Tab delimited table w/ columns Chrom Start End Count Count/median log2(Count/median)  ")
        print(" Example:")        
        print("\tref|NC_001133|    1010    1510    647     0.424419  -1.2364")
        print("\n")  
        print(" To visual with R:")
        print(" Edit file: add header and split by chromosome")
        print("\n library(ggplot2)")
        print(" d <- read.table( inputfile, header=T)")
        print(" qplot(d$pos,d$ct )" )
        print("")

        sys.exit(1)
    
    if cmdResults['WIN']:
        window = int(cmdResults['WIN'])
    else:
        window = 500          # default window size
    
    if cmdResults['BAM']:
        infile = cmdResults['BAM']   
        samName = re.sub(r"bam", "sam", infile )   
        # create sam file for parsing        
        if ( not os.path.exists(samName) ):
            cmd  = [ 'samtools', 'view', '-o', samName, infile ]
            output  = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
            
        first     = 0
        chrName   = ""
        startPos  = 0
        endPos    = window
        ctr       = 0
        with open(samName, 'r') as reads:
            for i in reads:
                line = i.split()
                # find first chromosome
                if first == 0:
                    chrName    = line[2]
                    first      = 1
                    ctr       += 1
                    chrom[chrName]['cnt'] = []
                    chrom[chrName]['pos'] = []
                # when chromosome changes
                elif chrName != line[2]:
                    chrName   = line[2]
                    startPos  = 0
                    endPos    = window
                    ctr       = 1
                    chrom[chrName]['cnt'] = []
                    chrom[chrName]['pos'] = []    
                # count if within current window
                elif ( int(line[3]) < endPos and int(line[3]) >= startPos ):
                    ctr += 1
                # change to next window
                else:
                    ctr += 1
                    windowPos = int(line[3]) - window
                    
                    chrom[chrName]['cnt'].append(ctr)
                    chrom[chrName]['pos'].append(windowPos)
                                        
                    startPos += window
                    endPos   += window
                    ctr       = 1
               
        # sometimes mapped reads have a chromosome named "*", get rid of those
        remove = [ key for key in chrom if key == '*']
        for key in remove:
            del chrom[key]
        
        chromList = chrom.keys()
            
        # print file header
        print("%s\t%s\t%s\t%s\t%s\t%s" %("chr", "Startpos","EndPos", "cnt", "cnt/median", "log2(cnt/median)" ))
                  
        # loop through all chromosomes in dictionary
        for item in chromList:
            cntMedian = median ( chrom[item]['cnt'] )

            for position, count in zip(chrom[item]['pos'], chrom[item]['cnt']):
                end = int(position) + window
                if count != 0:
                    ratio = count/cntMedian
                else:
                    ratio = 0
                
                if ratio == 0:
                    logRatio = "NA"
                else:
                    logRatio = math.log2(ratio)
                # write results to stdout
                print("%s\t%s\t%s\t%s\t%s\t%s" %(item, int(position), end, count, ratio, logRatio  ) )
                
            
            
        
          
if __name__ == "__main__":
    main()

