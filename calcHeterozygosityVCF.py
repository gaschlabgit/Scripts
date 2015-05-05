#!/home/mplace/anaconda3/bin/python
"""
Program: calcHeterozygosityVCF.py

Purpose: Calculate % heterozygosity on a private snp or indel vcf file.
         Private means, the sample in question has an allele not shared
         by other samples.  

Input  : Text file listing sample names, one per line , assumes you are in the directory
with the vcf files you want to process. 

Heterozygosity: 
    for SNPs:  anything other than 0,X  where is X is any number.
              i.e. 1,21  or 9,43  are heterozygous
    for Indels:  anything other than for examples 0,0,0,X  or 0,0,X  etc...
              i.e.  1,0,45,3 or 2,0,9 are heterozygous    

It is also assumed that the sample of interest is the only one in the file.
This may be achieved using:
      cat List | while read i; do vcf-subset -c $i ${i}.private.vcf > individual_vcfs/$i-snps.vcf; done

example:
273614N
378604X
BC187
 
Output: A table with the % heterozygosity for private snps or indels.

Results can be cleaned up for easier processing by running:

*******************************************************************************
@author: mplace
"""
import os
import re
import sys
import subprocess      
import argparse       

# list of files created
outFileList = set()
variant = set()
 
def calHet( inFile, varType ):
    """
    Calculate heterozygosity 
    
    samples = list of sample names
    
    vcf = VCF file    
    """
    names = []
    print("Sample\tfracHet\thetCt\thomCt")                  # print header
    
    with open( inFile, 'r') as files:                       # open sample name file
        for i in files:
            i = i.rstrip()
            vcf = i + "." + varType + ".vcf"                            
            with open( vcf, 'r' ) as data:
                hom = 0.0                                   # count homozygous sites
                het = 0.0                                   # count heterozygous sites
                fractionHet = 0.0                           # fraction heterozygous
                
                for var in data:
                    if var.startswith("#"):                 # skip header
                        continue
                    else:   
                        var = var.rstrip()
                        line = var.split("\t")
                        stats = line[9].split(':')          # 
                        alleles = list( map( int, stats[1].split(',') ) )     # create list of allele counts
                        check = [ i for i in alleles if i > 0]                # put any counts > 0 into a list
                        if not check:                                         # if all allele counts == 0
                            continue    # all alleles are set to zero wtf? Result of a quality score that is low.
                        elif len(check) > 1:                                  # multiple allele counts , must be heterozygous
                            het += 1    # more than one allele 
                        elif len(check) == 1:                                 # only one allele has a count
                            hom += 1
                            #print("%s\t%s\t%s\t%s\t%s\t%s" %(i, line[0], line[1], stats[0], stats[1], check ) )  
                if hom == 0:
                    fractionHet = 100
                else:
                    fractionHet = het/(hom + het)                             # calculate fraction heterozygous
                print("%s\t%f\t%f\t%f" %(i, fractionHet, het,hom ))                
                
    files.close()  

def main():
    cmdparser = argparse.ArgumentParser(description="Calculate % Heterozygosity on a list of private vcf files.",
                                        usage='%(prog)s -f samplelist.txt [optional args: -h ]' ,prog='calcHeterozygosityVCF.py'  )                                  
    cmdparser.add_argument('-f', '--File', action='store', dest='FILE', help='List of samples names process, text file, 1 per line.')
    cmdparser.add_argument('-t', '--Type', action='store', dest='TYPE', help='Variation type: snp or indel, if blank default = snp.',
                           default='snp')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
             
    if cmdResults['FILE']:
        inFile   = cmdResults['FILE']
        cwd      = os.getcwd() 
    
    if not os.path.exists(inFile):
        print("\n\t-f input file does not exist.\n")
        cmdparser.print_help()
        sys.exit(1)
    
    if cmdResults['TYPE'] == 'indel':
        varType = "indel"
    else:
        varType = "snp"
        
    
    calHet(inFile, varType)
        
if __name__ == "__main__":
    main()

