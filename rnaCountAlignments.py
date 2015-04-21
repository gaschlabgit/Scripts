#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
Created on Mon Apr 20 09:49:58 2015

@author: mplace

Program: rnaCountAlignments.py 

Purpose: Use samtools to count reads for all rRNA loci using positions found in 
         SGD GFF file.

Input  : list of bam files , one per line in a text file
 
Output: Table where columns are the rRNA loci and rows are the bam file name

The rRNA positions were taken from the SGD R64-1-1 GFF file using:
  gawk '{ if ($3 ~/rRNA/) print $1,$3,$4,$5, substr($9, 1,11)  }' saccharomyces_cerevisiae_R64-1-1_20110208_noFasta.gff
  
chrXII rRNA 451575 451785 ID=ETS2-1
chrXII rRNA 451575 458432 ID=RDN37-1
chrXII rRNA 451786 455181 ID=RDN25-1
chrXII rRNA 455182 455413 ID=ITS2-1
chrXII rRNA 455414 455571 ID=RDN58-1
chrXII rRNA 455572 455932 ID=ITS1-1
chrXII rRNA 455933 457732 ID=RDN18-1
chrXII rRNA 457733 458432 ID=ETS1-1
chrXII rRNA 458433 459675 ID=NTS2-1
chrXII rRNA 459676 459796 ID=RDN5-1
chrXII rRNA 459797 460711 ID=NTS1-2
chrXII rRNA 460712 460922 ID=ETS2-2
chrXII rRNA 460712 467569 ID=RDN37-2
chrXII rRNA 460923 464318 ID=RDN25-2
chrXII rRNA 464319 464550 ID=ITS2-2
chrXII rRNA 464551 464708 ID=RDN58-2
chrXII rRNA 464709 465069 ID=ITS1-2
chrXII rRNA 465070 466869 ID=RDN18-2
chrXII rRNA 466870 467569 ID=ETS1-2
chrXII rRNA 467570 468812 ID=NTS2-2
chrXII rRNA 468813 468931 ID=RDN5-2
chrXII rRNA 472465 472583 ID=RDN5-3
chrXII rRNA 482045 482163 ID=RDN5-4
chrXII rRNA 485697 485815 ID=RDN5-5
chrXII rRNA 489349 489469 ID=RDN5-6

chrMito rRNA 6546 8194 ID=15S_rRNA
chrMito rRNA 58009 62447 ID=21S_rRNA

chrXII = ref|NC_001144|
chrmt  = ref|NC_001224|

Dependencies:  requires samtools, expected to be in your $PATH 

samtools basic command to query reads aligned to a region:
    samtools view -c myFile.bam chr:startPos-endPos
    samtools view -c Gasch136_ATCACG_R1_MERGED.sort.bam "ref|NC_001144|:455933-457732"

***WARNING*** samtools requires that the bam file be accompanied by the index file (.bai)

"""
import argparse 
import os
import sys
import re
import subprocess 
from collections import defaultdict     


# rRNA gene names w/ start & end positions
rnaPos = { "ETS2-1" : "451575-451785", "RDN37-1" : "451575-458432", "RDN25-1" : "451786-455181",
       "ITS2-1" : "455182-455413", "RDN58-1" : "455414-455571", "ITS1-1" : "455572-455932",
       "RDN18-1" : "455933-457732", "ETS1-1" : "457733-458432", "NTS2-1" : "458433-459675",
       "RDN5-1" : "459676-459796", "NTS1-2" : "459797-460711", "ETS2-2" : "460712-460922",
       "RDN37-2" : "460712-467569", "RDN25-2" : "460923-464318", "ITS2-2" : "464319-464550",
       "RDN58-2" : "464551-464708", "ITS1-2" : "464709-465069", "RDN18-2" : "465070-466869",
       "ETS1-2" : "466870-467569", "NTS2-2" : "467570-468812", "RDN5-2" : "468813-468931",
       "RDN5-3" : "472465-472583", "RDN5-4" : "482045-482163", "RDN5-5" : "485697-485815",
       "RDN5-6" : "489349-489469", "15S_rRNA" : "6546-8194", "21S_rRNA" : "58009-62447"
       }

rnaNames = [ "15S_rRNA", "21S_rRNA", "ETS1-1", "ETS1-2", "ETS2-1", "ETS2-2", "ITS1-1",
             "ITS1-2", "ITS2-1", "ITS2-2", "NTS1-2", "NTS2-1", "NTS2-2", "RDN18-1", "RDN18-2",
             "RDN25-1", "RDN25-2", "RDN37-1", "RDN37-2", "RDN5-1", "RDN5-2", "RDN5-3", "RDN5-4",
             "RDN5-5", "RDN5-6", "RDN58-1", "RDN58-2"
           ] 
           
chr12  =  "ref|NC_001144|"          
chrM   =  "ref|NC_001224|"

def main():
    cmdparser = argparse.ArgumentParser(description="Systematically search and count all reads associated w/ all rRNA loci.",
                                        usage='%(prog)s <bamFileList.txt>' ,prog='rnaCountAlignments.py'  )
    cmdparser.add_argument('-f', '--file', action='store',     dest='FILE', help='File listing the bam files to process, one per line')         
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='List of rRNA positions used for search')
    cmdResults = vars(cmdparser.parse_args())
    
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)    
    
    table = defaultdict(list)    
    
    if cmdResults['INFO']:
        print("\n rnaCountAlignments.py")
        print("\n Purpose: Count all reads associated w/ all rRNA loci.")
        print("\n Input : A text file listing bam files, one per line.")
        print(" Generate input file by listing your files:  /bin/ls *.bam > input.txt\n")
        print(" Required Parameters: -f input.txt")
        print(" To run:  /home/GLBRCORG/mplace/scripts/rnaCountAlignments.py -f input.txt\n")
        print(" Output: Text table where column names = gene and row names = bam file name\n")
        print("\n")     
        sys.exit(1)
    
    if cmdResults['FILE']:
        infile = cmdResults['FILE']       
        if os.path.exists(infile):
            with open( infile, 'r') as f:
                for item in f:
                    item = item.rstrip()
                    print("Processing file: %s" %(item))
                    for gene in rnaNames:
                        positions = rnaPos[gene]
                        
                        if( re.match(r'21S_rRNA', gene) or re.match(r'15S_rRNA', gene) ):
                            region    = chrM + ":" + positions 
                        else:
                            region    =  chr12 + ":" + positions 
                            
                        cmd     = [ 'samtools', 'view', '-c', item, region ]
                        output  = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
                        result1 = output[0].decode( 'utf-8' )
                        result1 = result1.rstrip()
                        table[item].append(result1)
                    
                    totalCmd    = [ 'samtools', 'view', '-c', item ]
                    totalOutput = subprocess.Popen( totalCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
                    total       = totalOutput[0].decode( 'utf-8' )
                    total       = total.rstrip()
                    table[item].append(total)
                    
    
    fileHandle = open("rnaCounts.tab", 'w')
    #print header
    rnaNames.append("total")
    fileHandle.write("file\t%s\n" %("\t".join(rnaNames)))
   
    for key, value in table.items():
        row = "\t".join(table[key])
        fileHandle.write("%s\t%s\n" %(key,row) )



if __name__ == "__main__":
    main()

