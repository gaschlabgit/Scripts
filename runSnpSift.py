#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
Program:  runSnpSift.py

Purpose: run SnpSift on scarcity

Input  : annotated vcf files to process
 
     
Output: a new vcf file for each variant type produced by snpEff

See  http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf for definitions.

*******************************************************************************
SnpSift

SnpSift is a toolbox that allows you to filter and manipulate annotated files.

http://snpeff.sourceforge.net/SnpSift.html

*******************************************************************************
@author: mplace
"""

import os
import re
import sys
import subprocess      
import argparse       
import logging

# list of variant types
# to get a list of variants for your vcf file:
# run:   cat YOUR_annotated.vcf  | awk '{ n=split($8,array,"|"); print array[2];  }' | sort | uniq
# then edit file to remove blank lines

variant = { "disruptive_inframe_deletion", "disruptive_inframe_insertion", "disruptive_inframe_insertion&splice_region_variant",
"frameshift_variant", "frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&splice_region_variant&intron_variant",
"frameshift_variant&splice_region_variant", "frameshift_variant&start_lost",
"frameshift_variant&stop_gained", "frameshift_variant&stop_lost&splice_region_variant",
"inframe_deletion","inframe_insertion", "inframe_insertion&splice_region_variant",
"intergenic_region", "intragenic_variant", "intron_variant", "non_coding_exon_variant", "splice_donor_variant&intron_variant",
"splice_region_variant&intron_variant", "splice_region_variant&non_coding_exon_variant",
"start_lost&disruptive_inframe_deletion", "start_lost&inframe_deletion", "start_lost&inframe_insertion",
"stop_gained&disruptive_inframe_deletion", "stop_gained&disruptive_inframe_insertion",
"stop_gained&inframe_insertion", "stop_lost&disruptive_inframe_deletion&splice_region_variant",
"stop_lost&inframe_deletion&splice_region_variant", "initiator_codon_variant",
"intergenic_region", "intragenic_variant", "intron_variant", "missense_variant",
"missense_variant&splice_region_variant", "non_coding_exon_variant", "splice_acceptor_variant&intron_variant",
"splice_donor_variant&intron_variant", "splice_region_variant&intron_variant", "splice_region_variant&non_coding_exon_variant",
"splice_region_variant&splice_region_variant&intron_variant", "splice_region_variant&stop_retained_variant",
"splice_region_variant&synonymous_variant", "start_lost", "stop_gained", "stop_lost",
"stop_lost&splice_region_variant", "synonymous_variant"    }     

# list of files created
outFileList = set()
 
def runSnpSift( inFile, varFilt = variant ):
    """
    Call SnpSift
    Path is defined for scarcity:  /home/GLBRCORG/mplace/bin/snpEff/SnpSift.jar
    
    java -Xmx8g -jar /home/GLBRCORG/mplace/bin/snpEff/SnpSift.jar filter "ANN[0].EFFECT  = 'synonymous_variant'"
    Total-pass-Snps-rmRefGen_anno.vcf > SNP_synonymous_rmRefGen.vcf
    
    files = list of files to process
    """
    program  = "/home/GLBRCORG/mplace/bin/snpEff/SnpSift.jar"
    #program = "/home/mplace/bin/snpEff/SnpSift.jar"

    with open( inFile, 'r') as files:
        for i in files:
            logging.info(' Using file: %s' %(i) )
            i = i.rstrip()

            for var in varFilt:      
                name = re.sub(r"&", "_", var)
                item = "Ann[0].EFFECT = \'" + var + "\'"
                cmd     = [ 'java', '-Xmx4g', '-jar' , program, "filter", item, i   ]
           
                logging.info(' cmd: %s' %(cmd) )
                output  = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
                result1 = output[0].decode( 'utf-8' )
            
                outFile = re.sub(r"vcf", name , i)
                outFileList.add(outFile)
                with open( outFile, 'w' ) as out:
                    out.write( result1 )     



def main():
    cmdparser = argparse.ArgumentParser(description="Run SnpSift on a list of files.",
                                        usage='%(prog)s input.txt [optional args: -h -i ]' ,prog='runSnpSift.py'  )                                  
    cmdparser.add_argument('-f', '--File', action='store',      dest='FILE', help='List of files to process in a text file, one per line.')
    cmdparser.add_argument('-l', '--list', action='store',      dest='LIST', help='Variant types for filter, text file, one per line.' )
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='Print a more detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
             
    if cmdResults['INFO']:
        print("\n  runSnpSift.py")
        print("\n  Purpose: Run SnpSift on each file listed in file.")
        print("\t   For more SnpSift information see:")
        print("\t   http://snpeff.sourceforge.net/SnpSift.html")
        print("\n  Input  : Text filefiles.")
        print("\n  Usage  : runSnpSift.py -f input.txt ")
        print("\t   input.txt = text file with one vcf file name per line.\n")
        print("  Optinal arguments:")
        print("\t-l , --list  Variant type filter, one per line in text file.\n")
        print("\t-i , --info  Print this message.\n")
        print("\tList of variant types:")
        for i in variant:
            print("\t%s" %(i))
        print("")
        print("  To generate your own list of variant types run:")
        print("  cat YOUR_annotated.vcf  | awk \'{ n=split($8,array,\"|\"); print array[2];  }\' | sort | uniq\n" )
        print("  To use that list:\n\topen this program and replace the list in the dictionary  variant = { }")
        print("\n\tSee Mike Place for any problems or suggestions.")
        sys.exit(1)
    
    if cmdResults['FILE']:
        inFile   = cmdResults['FILE']
        cwd      = os.getcwd() 
        LOG_FILE =  'SnpSift.log'    
    
    if os.path.exists(LOG_FILE):
        os.unlink(LOG_FILE)
        logging.basicConfig(filename=LOG_FILE, level=logging.INFO,)
        logging.info(' DELETED previous log file')
    else:
        logging.basicConfig(filename=LOG_FILE, level=logging.INFO,)
    
    logging.info(' Running SnpSift')
    logging.info(' Directory %s' %(cwd))
    logging.info(' Input file: %s' %(inFile))
    
    if cmdResults['LIST']:
        logging.info(' Variant filter file: %s' %(cmdResults['LIST']) )
        logging.info(' Items to filter: ')
        varFilter = set()
        with open(cmdResults['LIST']) as item:
            for x in item:
                x = x.rstrip()
                varFilter.add(x)
        runSnpSift( inFile, varFilter)
    else:
        runSnpSift( inFile )
    
    logging.info("Files Created: %s" %(outFileList))
        
if __name__ == "__main__":
    main()

