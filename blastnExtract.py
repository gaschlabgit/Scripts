#!/home/mplace/anaconda3/bin/python
"""
@Program: blastnExtract.py

@Purpose: Extract a gene plus & minus user defined window from a sequence
          assembly.  Program assumes you are running in the blast file directory.

@Input:  1) User list of blast tab files, assumes header for each read looks
            something like :
                qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
                
                Assembly file are assumed to be in a different dir
                and named the same way, 
                
         2) Directory where sequence assemblies are located. Sequence Assemblies
            are assumed to have been indexed with samtools.
             example: samtools faidx sample.fa 
         
         3) Window around blast coordinates to extract from assembly, 1000 = 1000 base pairs 
            up & down stream of cds start and stop position.
         
@Dependencies: 
        Python 3, 
        BioPython, 
        samtools (for indexing assembly and extracting region)
      
@Output: ORF specific fasta files, one for each orf with Sequence extracted
         in fasta format by samtools
@author: Mike Place
@Date:   8/05/2015

"""
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections   import OrderedDict
import argparse	
import os
import re
import subprocess
import sys

def extract( refFile, strainName, assemblyDir, window ): #,  assemblyDir  ):
    """
    Use samtools faidx reference.fasta chr:start-end to extract cds region
    from reference.  Then call writeFasta().
    example:
    
        samtools faidx yjm555.fa Chr07:1060052-1060879
    
    BioPython Seq and SeqIO are used to write output.
    
    Each ORF file created has all the sequences across all strains
    for that ORF.
        
    """
    reference = assemblyDir + strainName + ".fa"                               
    print(strainName)                                                      # hold final sequences for printing
    for k,v in refFile.items():
        orfResults = k + "-results.fa"
        out = open(orfResults, 'a')
        rvse = False
        # check for start position > end position (reverse strand)
        if v[2] > v[3]:
            tmp = v[3]
            v[3] = v[2]
            v[2] = tmp
            rvse =True                
        # check for start position less than zero    
        if v[2] - window < 0:
            v[2] = 0
        else:
            v[2] = v[2] - window
        v[3] += window                                                         # Set right side of window                   
            
        loc = str(v[2]) + "-" + str(v[3])                                      # create start-end for samtools extraction
        # Now run samtools to extract sequence from reference assembly
        cmd = [ 'samtools', 'faidx', reference, v[0] + ":" + loc ]
        output = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
        result = output[0].decode('utf-8')
        # format result for writing to file
        tmpSeq  = result.split()
        seqName = re.sub(r">", "", tmpSeq.pop(0) )
        orfSeq  = ''.join(tmpSeq)
        orfSeq  = Seq(orfSeq)
            
        if rvse:            
            orfSeq  = orfSeq.reverse_complement()
            rvse = False
            
        record =  SeqRecord(orfSeq, id = seqName, description = strainName) 
        SeqIO.write(record, out, "fasta")      


def processFiles( blastFiles, assemblyDir, window ):
    """
    Loop through the blast tab files, get list of regions to extract for each file
    then call extract().  The region to extract is defined as the "winner", i.e.
    the best hit for a given ORF based on percent ID and bitscore
        blastFiles  = list of blast tab files to process one at a time
        assemblyDir = directory (different from cds file dir) contains sequence assemblies
        window      = window to extract on both sides of orf
    """    
    for f in blastFiles:
        with open(f) as blastF:
    # orderedDict, key = orf (ORF1) value = list, [0] = chrom [1] = %id [2] = start , [3] = end [4] = bitscore
            orf = OrderedDict() 
            strainName = f.split('.')[0]                                       # get strain name
            for line in blastF:
                if line.startswith('qseqid'):                                  # skip header
                    continue
                else:
                    line  = line.rstrip()
                    tabs  = line.split()
                    names = tabs[0].split('.')                                 # split strain.ORF#
                    if names[1] in orf:
                        if float(tabs[11]) > orf[names[1]][4] :                # is current bitscore > old bitscore
                            orf[names[1]][0] = tabs[1]                          
                            orf[names[1]][1] = float(tabs[2])
                            orf[names[1]][2] = int(tabs[8])
                            orf[names[1]][3] = int(tabs[9])
                            orf[names[1]][4] = float(tabs[11])
                    else:                                                      # if orf not a key, add it
                        bestHit = []
                        bestHit.append(tabs[1])                                # chrom
                        bestHit.append(float(tabs[2]))                         # percent id
                        bestHit.append(int(tabs[8]))                               # start position
                        bestHit.append(int(tabs[9]))                                # end position
                        bestHit.append(float(tabs[11]))                        # bitsscore
                        orf[names[1]] = bestHit
            
            extract( orf, strainName, assemblyDir, window )
            #for k,v in orf.items():
            #    print(k,v)

def main():
    """
    Main 
    """
    # Variables
    assemblies = []         # list of assembly files to extract sequence from
    blastFiles   = []         # cds files to process, in curr working dir
    
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
        print("\n  extractBlastn.py -f list.txt -w window [ -i ]")
        print("\n  Purpose: Extract a gene plus & minus user defined window from a sequence assembly.")
        print("\n  Input  : file w/ list of blastn table files, 1 per strain, path to assemblies, window size (default 1000 bp)")
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
            for blast in filelist:
                blast = blast.rstrip()
                blastFiles.append(blast)
    
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
    processFiles(blastFiles, assemblyDir, window)
    
    

if __name__ == "__main__":
    main()
    




