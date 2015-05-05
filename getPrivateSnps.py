#!/home/mplace/anaconda3/bin/python
"""
Program: getPrivateSnps.py

Purpose: run vcf-contrast -n  on scarcity

Input  : text file listing sample names, one per line

example:
273614N
378604X
BC187
CBS7960
CLIB324
 
Output: A new vcf file for each sample containing the private snps or indels.

command line example:
vcf-contrast +YJM1573 -273614N,378604X,BC187,YJM1078,YJM1083 -n Total-pass-Snps-rmRefGen_anno.vcf 

This will select the private variants for YJM1573 (w/ + ) and exclude the other
samples (w/ -).  The exclusion list can be very long.  In this program the list 
of strain names is run through selected an individual sample and excluding all 
others.  A vcf file will be created for each sample name.

vcftools is assumed to be in your path.

This script automates the use of vcf-contrast from vcftools.

http://vcftools.sourceforge.net/perl_module.html#vcf-contrast

About: Finds differences amongst samples adding NOVELGT, NOVELAL and NOVELTY annotations to INFO field.
       Note that haploid genotypes are internally treated as homozygous diploid genotypes, therefore
       "0/1" and "1" are considered different genotypes.

Usage: vcf-contrast +<list> -<list> [OPTIONS] file.vcf.gz

Options:
   +<list>                             List of samples where unique variant is expected
   -<list>                             List of background samples
   -d, --min-DP <int>                  Minimum depth across all -<list> samples
   -f, --apply-filters                 Skip sites with FILTER column different from PASS or "."
   -n, --novel-sites                   Print only records with novel genotypes
   -h, -?, --help                      This help message.
   
Example:
   # Test if any of the samples A,B is different from all C,D,E
   vcf-contrast +A,B -C,D,E -m file.vcf.gz

Results can be cleaned up for easier processing by running:
    nohup cat SampleNameList | while read i; do vcf-subset -c $i sampleName.private.vcf > individual_vcfs/$i-snps.vcf; done &

*******************************************************************************
@author: mplace
"""
import os
import re
import sys
import subprocess      
import argparse       
import logging

# list of files created
outFileList = set()
variant = set()
 
def vcfContrast( samples, vcf ):
    """
    Call vcf-contrast
    Path is defined for scarcity:  
    
    samples = list of sample names
    
    vcf = VCF file    
    vcf-contrast +YJM1573 -273614N,378604X,BC187,YJM1078,YJM1083 -n Total-pass-Snps-rmRefGen_anno.vcf 
    """
    program  = "vcf-contrast"
    names    = []

    with open( samples, 'r') as files:
        for i in files:
            i = i.rstrip()
            names.append(i)
    files.close()
    
    for x in range(len(names)):
        item    = names[x]
        sName   = item
        outFile     = item + ".private.vcf"
        item    = "+" + str(item)
        back    = names[x+1:]
        front   = names[:x]
        control = ",".join((back + front) )
        control = "-" + control
        
        cmd = [ program, item, control, '-n', vcf  ]
        logging.info(' cmd: %s' %(cmd) )
        output  = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
        result1 = output[0].decode( 'utf-8' )
        outFileList.add(outFile)
        with open( outFile, 'w' ) as out:
            out.write( result1 )      


def main():
    cmdparser = argparse.ArgumentParser(description="Run vcf-contrast on vcf file given a list of samples.",
                                        usage='%(prog)s -f samplelist.txt -v vcfFile [optional args: -h -i ]' ,prog='getPrivateSnps.py'  )                                  
    cmdparser.add_argument('-f', '--File', action='store',      dest='FILE', help='List of samples names process, text file, 1 per line.')
    cmdparser.add_argument('-v', '--vcf',  action='store',      dest='VCF',  help='VCF file.' )
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='Print a more detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
             
    if cmdResults['INFO']:
        print("\n  getPrivateSnps.py -f samplelist.txt -v vcfFile")
        print("\n  Purpose: Run vcf-contrast on a list of samples.")
        print("\t   For more vcf-contrast information see:")
        print("\t   http://vcftools.sourceforge.net/perl_module.html#vcf-contrast")
        print("\n  Input  : sample list file, vcfFile")
        print("\n  Output : Individual vcf files like: strainName.private.vcf")
        print("\n  Usage  : getPrivateSnps.py -f list.txt -v vcfFile [-i]")
        print("\t   list.txt = text file with one vcf file name per line.\n")
        print("  Optinal arguments:")
        print("\t-i , --info  Print this message.\n")
        print("\n  This produces a file with all sample columns present.")
        print("  For each file you are probably only interested in the specific sample's column.")
        print("  To clean up these files run something like the following:")
        print("  cat SampleNameList | while read i; do vcf-subset -c $i sampleName.private.vcf > $i-snps.vcf; done")
        print("\n\tSee Mike Place for any problems or suggestions.")
        sys.exit(1)
    
    if cmdResults['FILE']:
        inFile   = cmdResults['FILE']
        cwd      = os.getcwd() 
        LOG_FILE =  'vcfContrast.log' 
    
    if not os.path.exists(inFile):
        print("\n\t-f input file does not exist.\n")
        cmdparser.print_help()
        sys.exit(1)
        
    if cmdResults['VCF']:
        vcf = cmdResults['VCF']
    
    if not os.path.exists(vcf):
        print("\n\t-v vcf file does not exist.\n")
        cmdparser.print_help()
        sys.exit(1)    
    
    if os.path.exists(LOG_FILE):
        os.unlink(LOG_FILE)
        logging.basicConfig(filename=LOG_FILE, level=logging.INFO,)
        logging.info(' DELETED previous log file')
    else:
        logging.basicConfig(filename=LOG_FILE, level=logging.INFO,)
    
    logging.info(' Running vcf-contrast')
    logging.info(' Directory %s' %(cwd))
    logging.info(' Sample file : %s' %(inFile))
    logging.info(' VCF file    : %s' %(vcf) )
    
    vcfContrast(inFile, vcf)
    
    
    logging.info("Files Created: %s" %(outFileList))
        
if __name__ == "__main__":
    main()

