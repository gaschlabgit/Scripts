#!/home/mplace/anaconda3/bin/python
"""
Program: runCntrFreec.py

Purpose: This script automates the use of control-freec on a list of bam files.
         Samples run are assumed to not have a control sample.

Control-FREEC is a tool for detection of copy-number changes and allelic imbalances
(including LOH) using deep-sequencing data developed by the Bioinformatics
Laboratory of Institut Curie (Paris).

 http://bioinfo-out.curie.fr/projects/freec/

 To cite please use:

 Boeva V, Zinovyev A, Bleakley K, Vert JP, Janoueix-Lerosey I, Delattre O, Barillot E. (2011)
 Control-free calling of copy number alterations in deep-sequencing data using GC-content normalization.
 Bioinformatics 2011; 27(2):268-9. PMID: 21081509.

 Boeva V, Popova T, Bleakley K, Chiche P, Cappo J, Schleiermacher G, Janoueix-Lerosey I,
 Delattre O, Barillot E. (2011)
 Control-FREEC: a tool for assessing copy number and allelic content using next generation
 sequencing data. Bioinformatics. 2011 Dec 6. [Epub ahead of print] PubMed PMID: 22155870. 

Input  : text file, GC_content file, ( a chromosome length file is also needed, 
         but has been hard coded for scarcity S288C )
    text file listing sample names, one per line: 
    
W303.sorted.Realigned.bam
YPS163_JH.sorted.Realigned.bam

A GC content file must be created prior to running.

make the GC_content profile:
	in config_GC.txt change the following as needed:
		ploidy 
		window 
	then run:
	
	~/bin/gccount -conf config_GC.txt

config_GC.txt: 
	
	[general]
	chrFiles   = /home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/individu
alChroms
	chrLenFile = /home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/individu
alChroms/S288C_R64-1-1.len 
	ploidy = 2 
	GCcontentProfile = GC_profile_2kb.cnp
	window = 2000
	
Output:

 _CNVs: file with coordinates of predicted copy number alterations. Information for each region:

    chromosome
    start position
    end position
    predicted copy number
    type of alteration
    genotype (if [BAF] option is used)
    * precentage of uncertainty of the predicted genotyp (max=100).
    *&** status (somatic or germline)
    *&** precentage of germline CNA/LOH 

 _ratio.txt: file with ratios and predicted copy number alterations for each window. Information for each window:

    chromosome
    start position
    ratio
    median ratio for the whole fragment resulted from segmentation
    predicted copy number
    * median(abs(BAF-0.5)) per window
    * estimated BAF
    * genotype
    * precentage of uncertainty of the predicted genotype, max=100. 

* _BAF.txt: file B-allele frequencies for each possibly heterozygous SNP position. Information for each window:

    chromosome
    start position
    BAF
    fitted frequency of A-allele
    fitted frequency of B-allele
    inferred frequency of A-allele
    inferred frequency of B-allele
    precentage of uncertainty of the predicted genotype, max=100. 

 _sample.cnp and _control.cnp files: files with raw copy number profiles. Information for each window:

    chromosome
    start position
    number of read starts 

 GC_profile.cnp: file with GC-content profile. Information for each window:

    chromosome
    start position
    GC-content
    percentage of ACGT-letter per window (1-poly(N)%)
    *** percentage of uniquely mappable positions per window 

**** _ratio.BedGraph: file with ratios in BedGraph format for visualization in the UCSC genome browser. The file contains tracks for normal copy number, gains and losses, and copy neutral LOH (*). Information for each window:

    chromosome
    start position
    end position
    ratio*ploidy 

*   if [BAF] option is used
**  if [control] option is used
*** if gemMappabilityFile is used
**** if [general] BedGraphOutput=TRUE


Calculate significance of Control-FREEC predictions with an R script

One can add p-value to the predicted CNVs by running assess_significance.R:

    cat assess_significance.R | R --slave --args < *_CNVs > < *_ratio.txt > 

Visualize Control-FREEC's output

    cat makeGraph.R | R --slave --args < ploidy > < *_ratio.txt >

#*******************************************************************************
@author: mplace
"""
import argparse 
import os
import glob
import logging
import subprocess      
import sys

class CntrlFreec( object ):
    """
    This script automates the use of control-freec on a list of bam files.
    Control-Freec allows CNV calling w/o a control sample. 
    
    Required:
    chromosome length file formatted as:
        
    1       ref|NC_001133|  230218
    2       ref|NC_001134|  813184
    
    GC content file formatted as:
    
    ref|NC_001133|  0       0.4015  2
    ref|NC_001133|  2000    0.3735  2
    ref|NC_001133|  4000    0.331   2

    config file:
    
    [general]

    BedGraphOutput         = TRUE
    chrLenFile             = /home/mplace/Projects/copynumber/S288C_R64-1-1.len 
    coefficientOfVariation = 0.01
    GCcontentProfile       = /home/mplace/Projects/copynumber/GC_profile-2kb.cnp
    maxThreads             = 8 
    samtools               = samtools
    #Yeast telomeres are ~300bp , centeromeres ~150bp, telecentromeric value sets
    #exclusion zone around these regions
    telocentromeric        = 5000
    window                 = 2000
    outputDir              = /home/mplace/Projects/copynumber 
    ploidy                 = 2 
    breakPointType         = 4

    [sample]

    mateFile               = W303.sorted.Realigned.bam
    inputFormat            = BAM
    mateOrientation        = 0
    
    [control]
   
    """
    
    def __init__( self, inFile, currDir, gcprofile, window="2000" ):
        """
        inFile    = text file with bam files to process
        currDir   = current working directory
        gcprofile = file with GC content profiles
        window    = length in bp of bins, defaults to 2kb
        """
        self.file      = inFile         # bam files  
        self.currDir   = currDir
        self.gcprofile = gcprofile      # user supplied gc profile
        self.window    = window
        self.bamList   = []
        
        with open(inFile, 'r') as data:
            for item in data:
                item = item.rstrip()
                self.bamList.append(item)
                
    def printList( self ):
        """
        Print bam file names
        """
        for name in self.bamList:
            print(name)
    

    def runCF( self ):
        """
        Call control-freec
        Path is defined for scarcity:  
    
        samples = list of sample names to process
        
        Command Line example of running Control-Freec
        
        ~/bin/freec/freec -conf generic_config.txt
        """
        for item in self.bamList:
            self.writeConfig(item)
            program = "/home/mplace/bin/freec/freec"
            config  = "CF_config.txt"        
            cmd     = [ program, '-conf', config ]
                
            output  = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
            result1 = output[0].decode( 'utf-8' )
            logging.info(result1)        
 
    def writeConfig( self, bam ):
        """
        Control Freec requires a config file for each job.
        This is a template which inserts the bam file name
        and the output directory name.
        
        chrLenFile       = Chromosome len file, hard coded.
        GCcontentProfile = GC content file, supplied by user.
        
        Yeast telomeres are ~300bp, exclude from CNV calling.
        telocentromeric  = 5000
        
        """
        with open("CF_config.txt", 'w') as config:
            config.write("[general]\n")
            config.write("BedGraphOutput         = TRUE\n")
            config.write("chrLenFile             = /home/mplace/Projects/copynumber/S288C_R64-1-1.len\n") 
            config.write("coefficientOfVariation = 0.01\n")
            config.write("GCcontentProfile       = /home/mplace/Projects/copynumber/%s \n" %(self.gcprofile))
            config.write("maxThreads             = 8 \n")
            config.write("samtools               = samtools\n")      # CHANGE HERE FOR SAMTOOLS ON SCARCITY
            #Yeast telomeres are ~300bp , centeromeres ~150bp, telecentromeric value sets
            #exclusion zone around these regions
            config.write("telocentromeric        = 5000\n")
            config.write("window                 = %s \n" %(self.window) )
            config.write("outputDir              = %s \n" %(self.currDir) )
            config.write("ploidy                 = 2 \n")
            config.write("breakPointType         = 4 \n")
            config.write("\n")         
            config.write("[sample]\n")
            config.write("mateFile               = %s \n" %(bam))
            config.write("inputFormat            = BAM \n")
            config.write("mateOrientation        = 0 \n")
            config.write("\n")
            config.write("[control] \n")

    def matchGene( self ):
        """
        Use gff file to match CNV calls to known gene locations.
        
        """
        
        for name in glob.glob(  '/*.bam_ratio.txt'): #os.path.abspath(self.currDir) + '/*.bam_ratio.txt'):
            print( name )
            item = name.split(".")
            print(item[0])
            
                   
        
def main():
    cmdparser = argparse.ArgumentParser(description="Automates the use of control-freec for CNV calling.",
                                        usage='%(prog)s -f samplelist.txt -g GC_profile.cpn [-i]' ,prog='runCntrlFreec.py'  )                                  
    cmdparser.add_argument('-f', '--File', action='store',      dest='FILE', help='List of bam files to process, 1 name per line.',
                            metavar='')
    cmdparser.add_argument('-g', '--GC'  , action='store',      dest='GC'  , help='GC content file, see control-freec docs.',
                           metavar='')    
    cmdparser.add_argument('-w', '--Win' , action='store',      dest='WIN' , help='Window size', metavar='')
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='Print a more detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
             
    if cmdResults['INFO']:
        print("\n  runCntrlFreec.py -f samplelist.txt -g GC_profile.cpn -w window [ -i ]")
        print("\n  Purpose: Run CNV calling with Control-Freec.")
        print("\t   For more Control-Freec information see:")
        print("\t   http://bioinfo-out.curie.fr/projects/freec/")
        print("\n  Input  : bam file list file, GC content file, window size (default 2kb)")
        print("\n  Output : 5 files for each sample.")
        print("    _CNVs           : coordinates of predicted copy number alterations." )
        print("    _ratio.txt      : ratios and predicted copy number alterations for each window.")
        print("    _sample.cnp     : raw copy number profiles.")
        print("    _ratio.BedGraph : ratios in BedGraph format for visualization.")
        print("    .png            : R plot of CNV by chromosome.")
        print("\n  Usage  : runCntrlFreec.py -f list.txt [-i]")        
        print("  ")       
        print("\tTo see Python Docs for this program:")
        print("\n\tOpen python console and enter")
        print("\timport sys")
        print("\tsys.path.append('/full/path/to/script')")
        print("\timport runCntrlFreec")
        print("\thelp(runCntrlFreec)")
        print("\n\tSee Mike Place for any problems or suggestions.")
        sys.exit(1)
    
    if cmdResults['FILE']:
        inFile   = cmdResults['FILE']
        cwd      = os.getcwd() 
        LOG_FILE =  'runCntrlFreec.log'
    
    if not os.path.exists(inFile):
        print("\n\t-f input file does not exist.\n")
        cmdparser.print_help()
        sys.exit(1)
    
    if cmdResults['GC']:
        gcfile = cmdResults['GC']
    else:
        print("\n\t-g GC content file is required.\n")
        cmdparser.print_help()
        sys.exit(1)
    
    if not os.path.exists(gcfile):
        print("\n\t-g GC content file does not exist.\n")
        cmdparser.print_help()
        sys.exit(1)       
    
    if os.path.exists(LOG_FILE):
        os.unlink(LOG_FILE)
        logging.basicConfig(filename=LOG_FILE, level=logging.INFO,)
        logging.info(' DELETED previous log file')
    else:
        logging.basicConfig(filename=LOG_FILE, level=logging.INFO,)
    
    if cmdResults['WIN']:
        window = cmdResults['WIN']
        job = CntrlFreec(inFile, cwd, gcfile, window)
    else:
        job = CntrlFreec(inFile, cwd, gcfile)        
     
    logging.info(" Bam file list: %s" %(inFile) )
    logging.info(" Working Directory: %s" %(cwd) )
    logging.info(" GC Content file: %s" %(gcfile) )
    
    job.runCF()    
    job.matchGene()

if __name__ == "__main__":
    main()

