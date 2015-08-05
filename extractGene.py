#!/home/mplace/anaconda3/bin/python
"""
@Program: extractGene.py

@Purpose: Extract a gene plus & minus user defined window from a sequence
          assembly.  
          Program assumes you are running in the cds file directory.

@Input:  1) User list of cds fasta files, assumes header for each read looks
            something like :
                >yjm1252.ORF1 yjm1252_Chr05 6747-7661
                
                fasta files name format =  sample.fa
                Assembly file are assumed to be in a different dir
                and named the same way, 
                
         2) Directory where sequence assemblies are located. Sequence Assemblies
            are assumed to have been indexed with samtools.
             example: samtools faidx sample.fa 
         
         3) Window around cds to extract from assembly, 1000 = 1000 base pairs 
            up & down stream of cds start and stop position.
         
@Dependencies: 
        Python 3, 
        BioPython, 
        samtools (for indexing assembly and extracting region)
      
@Output: DNA fasta files
@author: Mike Place
@Date:   7/30/2015

"""
from Bio import SeqIO
import argparse	
import os
import re
import subprocess
import sys

def extract( refFile, orf, loc, chrom, strain, assemblyDir  ):
    """
    Use samtools faidx reference.fasta chr:start-end to extract cds region
    from reference.  Then call writeFasta().
    example:
    
        samtools faidx yjm555.fa Chr07:1060052-1060879
    
    returns:
    
        >Chr07:1060052-1060879program = '/opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar'
        
    """
    #print( refFile, orf, strain, chrom, loc )
    reference = assemblyDir + "/" + refFile
    cmd    = [ 'samtools', 'faidx', reference, chrom + ":" + loc ]  
    output = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
    result = output[0].decode('utf-8')
    
    return result             # type = str
    

def processFiles( files, assemblyDir, window ):
    """
    Loop through the cds files, get list of regions to extract for each file
    then call extract().  The region to extract is expected to be in each
    fasta sequence header.
        files       = cds files to process
        assemblyDir = directory (different from cds file dir) contains sequence assemblies
        window      = window to extract on both sides of orf
    """
    for f in files:
        #outFile =re.sub(r"fa","gene.fa", f)
        #with open(outFile, 'w') as output:
        print(f)
        for seqrec in SeqIO.parse( f, "fasta" ):
            nameTmp   = seqrec.id.split('.')
            orfName   = nameTmp[1]
            orfInfo   = seqrec.description.split()
            orfInfo.pop(0)
            # set window positions for getting region
            orfTmp    = orfInfo[0].split('_')
            orfLocTmp = [int(x) for x in orfInfo[1].split("-")]
            # set left site, can't have a value less than zero
            # don't worry about having a right side greater than chrom length
            # samtools just returns the end of the chromosome.
            if orfLocTmp[0] - window < 0:
                orfLocTmp[0] = 0
            else:
                orfLocTmp[0] = orfLocTmp[0] - window
            orfLocTmp[1] += window    # add window to right side
                
            orfChrom  = orfTmp[1]
            strain    = orfTmp[0]
            orfLoc    = str(orfLocTmp[0]) + "-" + str(orfLocTmp[1])
            seq = extract( f , orfName, orfLoc, orfChrom, strain, assemblyDir )
            seqRegion = seq.replace('>', '>' + strain + " " + orfName + " ")
            outFile = orfName + "." + "fa"
            with open(outFile, 'a') as output:
                if not re.search(r'[GATC]{10,100}$', seqRegion):                   
                    seqRegion = seqRegion.rstrip() + " BAD\n"
                    
                output.write(seqRegion)

def main():
    """
    Main 
    """
    # Variables
    assemblies = []         # list of assembly files to extract sequence from
    cdsFiles   = []         # cds files to process, in curr working dir
    
    cmdparser = argparse.ArgumentParser(description="Extract gene/cds sequence from Yeast Genome assemblies",
                                        usage='%(prog)s -f <list_file> -d <directory> [-w size_bp]' ,prog='extractGene.py'  )                                  
    cmdparser.add_argument('-f', '--file', action='store',      dest='FILE', 
                           help='List of gene/cds fasta files to process.' , metavar='')
    cmdparser.add_argument('-d', '--dir' , action='store',      dest='DIR' , help='Directory path containing assembly fasta files.', metavar='')
    cmdparser.add_argument('-w', '--win' , action='store',      dest='WIN' , help='Window size, optional default= 1000bp', metavar='')
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='Detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
             
    if cmdResults['INFO']:
        print("\n  extractGene.py -f list.txt -w window [ -i ]")
        print("\n  Purpose: Extract a gene plus & minus user defined window from a sequence assembly.")
        print("\n  Input  : list file, full path to assemblies, window size (default 1000 bp)")
        print("\n  Output : fasta files for each orf.")    
        print("\n  Usage  : extractGene.py -f list.txt -d /home/GLBRCORG/user/assemblyDir -w 200")        
        print("  ")       
        print("\tTo see Python Docs for this program:")
        print("\n\tOpen python console and enter")
        print("\timport sys")
        print("\tsys.path.append('/full/path/to/script')")
        print("\timport extractGene")
        print("\thelp(extractGene)")
        print("\n\tSee Mike Place for any problems or suggestions.")
        sys.exit(1)
    
    if cmdResults['FILE']:
        inFile = cmdResults['FILE']
        
    if not os.path.exists(inFile):
        print("\n\t-f input file does not exist.\n")
        cmdparser.print_help()
        sys.exit(1)
    else:
        # open file containing a list of cds files to process.
        with open( inFile, 'r') as filelist:
            for cds in filelist:
                cds = cds.rstrip()
                cdsFiles.append(cds)
    
    if cmdResults['DIR']:
        assemblyDir = cmdResults['DIR']
    
    if not os.path.exists(assemblyDir):
        print("\n\t-d assembly directory does not exist.\n")
        cmdparser.print_help()
        sys.exit(1)
    else:
        # Here we are assuming the .fai index file for reference exists
        for f in os.listdir(assemblyDir):
            if re.search(r'\.fa$',f):
                assemblies.append(f)
    
    if cmdResults['WIN']:
        window = int(cmdResults['WIN'])
    else:
        window = 1000

    # start the real work                 
    processFiles(cdsFiles, assemblyDir, window)
    
    

if __name__ == "__main__":
    main()
    




