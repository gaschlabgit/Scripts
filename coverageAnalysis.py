#!/home/mplace/anaconda3/bin/python
"""
Program: coverageAnalysis.py 

Purpose: Use samtools to count reads in a user specified window size.
         Default window size is 500 bp.

Input  : bam file
 
Output : Table where columns are chrom, start_position, end_position, count, count/median_count, log2(count/median_count)

dependencies: samtools, must be in path

To visual with R:

read in file, for example separate out a chromosome.
data looks like:

chr	Startpos	EndPos	cnt	cnt/median	log2(cnt/median)
chr9	0	500	10583	2.3591172536781095	1.2382471256973449
chr9	500	1000	7440	1.6584930896121266	0.7298730009091214

library(ggplot2)
library(scales)

qplot(chr8$Startpos,chr8$log2.cnt.median) + geom_point() 
+ ggtitle("Y22-3 Chr8 CNVs 500 bp windows") 
+ labs(x="Base Pair Position", y="log2(count/median(count))") 
+ scale_x_continuous( breaks = seq(0,max(chr8$Startpos), by = 100000 ), labels = comma )


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
from statistics import stdev

# key = chromosome  each value is a dict of lists, the keys are 'pos', 'cnt'  
chrom    = defaultdict(dict)

# used to convert chrom names if required
chromDict = {  "ref|NC_001133|" : "chr1", "ref|NC_001134|" : "chr2", "ref|NC_001135|" : "chr3",
           "ref|NC_001136|" : "chr4", "ref|NC_001137|" : "chr5", "ref|NC_001138|" : "chr6", 
           "ref|NC_001139|" : "chr7", "ref|NC_001140|" : "chr8", "ref|NC_001141|" : "chr9",
           "ref|NC_001142|" : "chr10", "ref|NC_001143|" : "chr11", "ref|NC_001144|" : "chr12",
           "ref|NC_001145|" : "chr13", "ref|NC_001146|" : "chr14", "ref|NC_001147|" : "chr15",
           "ref|NC_001148|" : "chr16", "ref|NC_001224|" : "chr17" , "ref|NC_001133|[R64]" : "chr1",
           "ref|NC_001134|[R64]" : "chr2", "ref|NC_001135|[R64]" : "chr3",
           "ref|NC_001136|[R64]" : "chr4", "ref|NC_001137|[R64]" : "chr5", "ref|NC_001138|[R64]" : "chr6", 
           "ref|NC_001139|[R64]" : "chr7", "ref|NC_001140|[R64]" : "chr8", "ref|NC_001141|[R64]" : "chr9",
           "ref|NC_001142|[R64]" : "chr10", "ref|NC_001143|[R64]" : "chr11", "ref|NC_001144|[R64]" : "chr12",
           "ref|NC_001145|[R64]" : "chr13", "ref|NC_001146|[R64]" : "chr14", "ref|NC_001147|[R64]" : "chr15",
           "ref|NC_001148|[R64]" : "chr16", "ref|NC_001224|[R64]" : "chr17" } 

# use to print ordered result to file
chromSet = [ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
             "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
             "chr16", "chr17" ]
             
def main():
    """
    Read in bam file, parse and count reads for each non-overlaping bins.
    For each bin calculate 
        ratio = (count/median count for chromosome)
        log2 of ratio
        print to stdout: chr   Startpos   EndPos   cnt   cnt/median  log2(cnt/median)
        
        Write Copy number calls to file.
        Proposed CNVs are filtered using:  
            
            abs(logR) > (1 + 2*chromStDev)   
    
    """
    cmdparser = argparse.ArgumentParser(description="Calculate Copy Number Variants using a sorted bam file.",
                                        usage='%(prog)s <file.bam>' ,prog='coverageAnalysis.py'  )
    cmdparser.add_argument('-b', '--bam',    action='store',      dest='BAM',  help='bam file', metavar='') 
    cmdparser.add_argument('-w', '--win',    action='store',      dest='WIN',  help='Window size for counting reads', metavar='') 
    #cmdparser.add_argument('-c', '--cutoff', action='store',      dest='CUT',  help='CNV Cutoff Score', metavar='' )       
    cmdparser.add_argument('-i', '--info',   action='store_true', dest='INFO', help='Program information')
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
        print(" Optional Parameters:")
        print(" \t-w window size")
        #print(" \t-c cutoff score [Default = .6]")
        print("\n Input : Bam file")
        print(" To run:  /home/mplace/scripts/coverageAnalysis.py -b samplename.bam \n")
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
    # handle bin/window size argument
    if cmdResults['WIN']:
        window = int(cmdResults['WIN'])
    else:
        window = 500          # default window size
        
    # handle CNV cutoff score 
    #if cmdResults['CUT']:
    #    score = float(cmdResults['CUT'])
    #else:
    #    score = 0.60
    
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
                    chrom[chrName]['cnt']      = []
                    chrom[chrName]['pos']      = []
                    chrom[chrName]['logRatio'] = []
                # when chromosome changes
                elif chrName != line[2]:
                    chrName   = line[2]
                    startPos  = 0
                    endPos    = window
                    ctr       = 1
                    chrom[chrName]['cnt'] = []
                    chrom[chrName]['pos'] = [] 
                    chrom[chrName]['logRatio'] = []
                # count if within current window
                elif ( int(line[3]) <= endPos and int(line[3]) >= startPos ):
                    ctr += 1
                # change to next window
                else:
                    # need to check the position 
                    ctr += 1
                    
                    chrom[chrName]['cnt'].append(ctr)
                    chrom[chrName]['pos'].append(startPos)
                    chrom[chrName]['logRatio'].append(0.0)
                                        
                    startPos += window
                    endPos   += window
                    ctr       = 0
                    while int(line[3]) >= endPos:
                        chrom[chrName]['cnt'].append(ctr)
                        chrom[chrName]['pos'].append(startPos + window)
                        chrom[chrName]['logRatio'].append(0.0)
                        startPos += window
                        endPos   += window
                        ctr       = 1
                        
               
        # sometimes mapped reads have a chromosome named "*", get rid of those
        remove = [ key for key in chrom if key == '*']
        for key in remove:
            del chrom[key]
        
        # Get list of chromsomes
        chromList = chrom.keys()
            
        # print file header
        print("%s\t%s\t%s\t%s\t%s\t%s" %("chr", "Startpos","EndPos", "cnt", "cnt/median", "log2(cnt/median)" ))
        
        # open file to write genomic regions that pass cutoff
        cutOff    = re.sub(r"bam", "calls", infile )
        cutOffout = open(cutOff, 'w' )
        cutOffout.write("chrom\tstart\tend\tcount\tlog2Ratio\tchromStdDev\n")
        result    = defaultdict(list)
        
                  
        # loop through all chromosomes in dictionary
        for item in chromList:
            cntMedian  = median( chrom[item]['cnt'] )
            logList    = []          # used to calculate standard deviation
            logPos     = 0           # used to index chrom[item]['logRatio']

            for position, count, in zip(chrom[item]['pos'], chrom[item]['cnt']):
                end = int(position) + window
                if count != 0:
                    ratio = count/cntMedian
                else:
                    ratio = 0
                
                if ratio == 0:
                    logRatio = math.log2(1)    # can't take the log of 0 so use 1 instead
                else:
                    logRatio = math.log2(ratio)
                
                chrom[item]['logRatio'][logPos] = logRatio
                logPos += 1      
                logList.append(logRatio)
                # write results to stdout
                print("%s\t%s\t%s\t%s\t%s\t%s" %(item, int(position), end, count, ratio, logRatio  ) )
            chromStDev = stdev(logList)
            
            for position, count, logR in zip(chrom[item]['pos'], chrom[item]['cnt'], chrom[item]['logRatio']):
                if ( abs(logR) > (1 + 2*chromStDev) ): 
                    end = int(position) + window
                    answer = str(position) + "\t" + str(end) + "\t" + str(count) +"\t" + str(logR) + "\t" + str(chromStDev)
                    result[chromDict[item]].append(answer) 
        # write final calls to file in chrom order       
        for i in chromSet:
            for x in result[i]:
                cutOffout.write("%s\t%s\n" %(i,x))
            
         
        cutOffout.close()
        
          
if __name__ == "__main__":
    main()

