#!/mnt/gv0/Groupspace/GaschLab/anaconda3/bin/python
"""
@Program: convertGene.py

@Purpose: Create a strain specific gene sequence by replacing snps
          found vcf into the reference gene sequence.

@Input:  Gene Name in yeast ORF format or as standard SGD gene name 
            example: YAL046C = AIM1
         List of strain names or all strains available
         VCF file provided by user
	   Length in bps up & down stream to include with gene sequence.
         User coordinates can be used as well instead of gene name
         WARNING:  
             Be aware that any DNA codon with a IUPAC value that is not the
             standard will be output as a nonstandard character.

@Dependencies: Python 3
               Must have vcftools install and have the perl library in the path.
               One way to ensure vcftools lib is in the path:
               At the command line run:
               export PERL5LIB=/FULLPATHTO/bin/vcftools_0.1.12b/lib/perl5/site_perl
                    This is set in os.environ (see below)               
               
               gffutils information at https://pythonhosted.org/gffutils/index.html
               yeast_Gene_name_to_ORF must be present (maps gene names to SGD ORF name)
               path to reference gff set in refDict (see below)
               path to reference directory set at pathToRef (see below)

               
@Installation:  Need to set the path to reference directory to find the gff
                databases created by gffutils. 

@Important Variables: 
        reference  -- reference genome name , R64-1-1
        vcfFile    -- vcf file name
        window     -- number of base pairs left and right of gene to include
        strainList -- list of strain names sequences to output
        outType    -- dna or aa
        geneInfo   -- Chrom, left & right positions and strand (+,-) for the gene, from gff file  
                                      

@Output:  DNA or AA fasta file, default is dna
@author: Mike Place
@Date:   2/5/2015

"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import gffutils
import argparse	
import os
import re
import sys

import yeast_Gene_name_to_ORF      # must be present, contains gene name mapping to ORF name

# vcftools requires vcf.pm in path 
os.environ['PERL5LIB'] = '/mnt/gv0/AD/mplace/bin/vcftools_0.1.12b/lib/perl5/site_perl'
pathToRef              = "/mnt/gv0/AD/mplace/reference"


# paths to reference gff files
refDict = { 'R64-1-1' : "/mnt/gv0/AD/mplace/reference/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208.gff",
            'R64-2-1' : "/mnt/gv0/AD/mplace/reference/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff"
          }
# paths to reference fasta files
fastaDict = { 'R64-1-1' :   "/mnt/gv0/AD/mplace/reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa",         
              'R64-2-1' :   "/mnt/gv0/AD/mplace/reference/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa"
    }
    
chromDict = { "I" : "ref|NC_001133|", "II" : "ref|NC_001134|", "III" : "ref|NC_001135|",
              "IV" : "ref|NC_001136|", "V" : "ref|NC_001137|", "VI" : "ref|NC_001138|", 
              "VII" : "ref|NC_001139|", "VIII" : "ref|NC_001140|", "IX" : "ref|NC_001141|",
              "X" : "ref|NC_001142|", "XI" : "ref|NC_001143|", "XII" : "ref|NC_001144|",
              "XIII" : "ref|NC_001145|", "XIV" : "ref|NC_001146|", "XV" : "ref|NC_001147|",
              "XVI" : "ref|NC_001148|", "Mito" : "ref|NC_001224|", "mt" : "ref|NC_001224|"
}
           
# Create a set for strain names
strainNames ={"273614N",	"378604X", "BC187", "CBS7960",	"CLIB324", "CLIB382", "DBVPG1106", 
"DBVPG1373", "DBVPG1788", "DBVPG1853", "DBVPG6044", "DBVPG6765", "DCM16",	"DCM6", "EthanolRed",
"FL100", "I14", "IL01", "K1", "K10",	"K11",	"K9",	"L1374", "L1528", "NC02", "NCY3290",
"NCY3455", "NCYC110", "NCYC361", "PE2", "PW5", "SB", "SK1", "T73", "UC5",	"UWOPS03_461_4",
"UWOPS05_217_3", "UWOPS05_227_2", "UWOPS83_787_3", "UWOPS87_2421", "W303", "WE372", "Y1",
"Y12",	 "Y2", "Y2189", "Y2209",	"Y3", "Y389", "Y55", "Y6", "Y7568", "YB210", "YB4081",
"YB4082", "YB908", "YJM1129", "YJM269", "YJM308", "YJM320", "YJM326", "YJM339",
"YJM421", "YJM428", "YJM440", "YJM451", "YJM454", "YJM653", "YJM975", "YJM978",
"YJM981", "Yllc17_E5", "YPS1000", "YPS1009", "YPS128", "YPS163", "YPS606", "YS2", "YS4", "YS9"
}

# assumes key  is in alphabetical order
iupacDNA = { 'AG':'R', 'CT':'Y', 'CG':'S', 'AT':'W', 'GT':'K' , 'AC':'M',
            'CGT':'B', 'AGT' : 'D', 'ACT':'H', 'ACG':'V', 'ACGT':'N'      
}

# yeast chromosome sizes, sizes taken from SGD (http://www.yeastgenome.org/strain/S288C/overview)
yeastSize = { "I" : 230218, "II" : 813184, "III" : 316620, "IV" : 1531933, "V" : 576874,
             "VI" : 270161, "VII" : 1090940, "VIII" : 562643, "IX" : 439888, "X" : 745751,
             "XI" : 666816, "XII" : 1078177, "XIII" : 924431, "XIV" : 784333, "XV" : 1091291,
             "XVI" : 948066, "Mito" : 85779, "mt" : 85779
}

def printStrainList():
    """
    Print list of all supported strain names
    """
    print("\n\tCurrent Strain List:")
    for s in strainNames:
        print("  %s" %(s) )
    print("\n\tSee Mike Place about adding strain names")

def searchGff( gene, gffFile, reference ):
    """
    Parse the reference gff file for specified gene ( in yeast ORF format, YAL046C)
    Return left & right positions, chromosome and strand as a list
    use gfftutils     
    """
    refdb = reference + '.db'
    geneInfo = []
    
    # if gff database doesn't exist create one
    if not os.path.exists(refdb):
        db = gffutils.create_db(gffFile, dbfn=refdb, force=True)

    try:
        tdb = gffutils.FeatureDB(refdb, keep_order=True)
        item = tdb[gene] 
        geneInfo.append( item.chrom.replace('chr','') )
        geneInfo.append( item.start -1 )
        geneInfo.append( item.end )
        geneInfo.append( item.strand )
        return geneInfo
    except gffutils.exceptions.FeatureNotFoundError:
        return geneInfo
    
def extractVCF( strainVCF, seqRec, chrom, left, right, window, strain ):
    """
    Extract the segment of VCF file and replace bases in reference sequence.
    
    Resulting in a strain specific gene sequence.   
    
    To get index of base to replace --
        ( SNP_POSITION - GENE_START ) + WINDOW == seq_index_for_base
    
    """
    # can't use negative left boundary
  
    if left - window < 0:
        lf = 0
    else: 
        lf = left - window
    rt = right + window
    
    base = ""                                                                  # hold the vcf snp base
    strainSeq  = seqRec.tomutable()                                             # copy reference sequence as parameters are passed by reference in python    
        
    
    with open(strainVCF) as fp:
        for line in fp:
            line = line.rstrip()
            row = line.split("\t")
            # find correct chrom and left & right bounds in VCF file
            if chrom == row[0]:             
                if int(row[1]) >= lf and int(row[1]) <= rt:
                    gt = row[9].split(":")
                    newPos = ( int(row[1]) - left -1 ) + window
                    baseList = row[4].split(",")                               # GT (genotypes can have mutliple values: A,T or G,T,C sometimes)
                    if gt[0] == "0/0":
                        pass                                                   # DO NOTHING, Ref Allele no change to make                        
                    elif gt[0] == "1/1":
                        base = baseList[0]
                        strainSeq[newPos] = base
                        #print("HOMOZYGOUS  ALT %s  base: %s  np: %d" %(gt[0], base, newPos))
                        #print(line)
                    elif gt[0] == "0/1":
                        base = ''.join( sorted( row[3]+baseList[0] ) )
                        code = iupacDNA[base]
                        strainSeq[newPos] = code
                        #print("HETEROZYGOUS %s    base: %s iupac: %s np: %d"  %(gt[0], base, code, newPos))
                        #print(line)
                    elif gt[0] == "0/2":
                        base = ''.join( sorted( row[3]+baseList[0] ) )
                        code = iupacDNA[base]
                        strainSeq[newPos] = code
                        #print("HETEROZYGOUS %s    base: %s iupac: %s np: %d"  %(gt[0], base, code, newPos))
                        #print(line)
                    elif gt[0] == "2/2":
                        base = baseList[1]
                        strainSeq[newPos] = base
                        #print("HOMOZYGOUS  ALT %s  base: %s  np: %d" %(gt[0], base, newPos))
                        #print(line)                                                                
                    else:
                        print("UNKNOWN  %s" %(gt[0]))

    newStrainSeq = strainSeq.toseq()
    
    return newStrainSeq            
 
def extractRefSeq( ref, chrom, left, right, window, strand ):
    """
    Given fasta file input, extract a sequence and create a SeqRecord object
    need: reference name, chrom, left & right positions, window
    Window not supported for proteins.
    """
    print("\n\tExtracting reference sequence from fasta .....")
    ncChrom = chromDict[chrom]  
    # check window boundaries, set defaults if bounds are exceeded
    if left - window < 0:
        leftWindow = 0
    else:
        leftWindow = left - window
    
    if right + window > yeastSize[chrom]:
        rightWindow = yeastSize[chrom]
    else:
        rightWindow = right + window      
        
    if ref in fastaDict:
        print("\tRef File: %s" %(fastaDict[ref]))
        for record in SeqIO.parse(fastaDict[ref], "fasta"):
            if ncChrom == record.id:
                geneSeq   = record.seq[left:right]
                if strand == '-':
                    geneSeq = geneSeq.reverse_complement()
                # check for start and stop codon on the extracted Gene Sequence
                check = checkSeq(geneSeq, chrom)
                if check != "":
                    print("\n\tInvalid Start or Stop codon in Gene sequence")
                    print("Sequence Found: ")
                    print(geneSeq)
                    print("\t%s" %(check))
                    sys.exit(1)
                
                # Gene checks out so now extract entire window
                geneSeq = record.seq[leftWindow:rightWindow]
                if strand == '-':
                    geneSeq.reverse_complement()    
                return geneSeq
    else:
        return None
                #print(geneSeq)  

def checkSeq( seq, chrom ):
    """
    Check for the presence of start (ATG) on the left end.
            Yeast Mitochondria start codes differ (ATG & ATA)
    Check for the presence of stop  (TAA,TAG,TGA) on the right end.
    """
    valid = ""
    start = seq[0:3]
    length = len(seq)
    end = seq[length-3:length]     
    
    if chrom == "Mito" or chrom == "mt":
        if not re.search(r"ATG|ATA", str(start)):
            valid = "Invalid Start codon: " + start + " " 
    else:
        if not re.search(r"ATG", str(start)):
            valid = "Invalid Start codon: " + start + " "        
         
    if not re.search(r"T(AA|GA|AG)", str(end) ):
        valid = valid + "\n\tInvalid Stop: " + end
        
    return valid


def seqTranslate( seq, chrom ):
    """
    Convert a DNA sequence to protein sequence
    Standard codons used for chr I-XVI
    Yeast Mitochondrial codons used for Mito
    """
    remainder = len(seq)%3                                                     # must be divisible by 3 to translate
    if remainder:
        newRecord = seq[remainder:]                                            # remove the remainder from beginning of sequence 
    else:
        newRecord = seq          
             
    if chrom == 'Mito' or chrom == 'mt':
        recordAA = newRecord.translate(table="Yeast Mitochondrial")
    else:
        recordAA = newRecord.translate()
    
    return recordAA

def writeRefFasta( rec, geneName, chrom, left, right ):
    """
    Write Gene Reference sequence to file.
    """
    start = int(left) + 1
    left = str(start)
    file = geneName + ".reference.fasta"
    geneRec = SeqRecord(rec)
    geneRec.id = geneName + "|chr" + chrom + ":" + str(left) + "-" + str(right)
    geneRec.name = geneName
    geneRec.description = "Seq extracted by convertGene.py"
    SeqIO.write(geneRec, file, "fasta")
    return file
    
def writeFasta( rec, geneName, chrom, left, right ):
    """
    Write gene specific results to fasta file
    """
    start = int(left) + 1
    left  = str(start)
    file  = "Sequence_results.fasta"
    geneRec    = SeqRecord(rec)
    geneRec.id = geneName + "|chr" + chrom + ":" + str(left) + "-" + str(right)
    geneRec.name        = geneName
    geneRec.description = "Seq extracted by convertGene.py"
    output = open(file, "a")
    SeqIO.write(geneRec, output, "fasta")


def verifyStrainName( strains ):
    """
    Check for valid strain name by searching strainNames set
    """
    for l in strains:
        if l not in strainNames:
            return l
    return None 
        
def verifyRef( userRef,  pathToRef ):
    """
    Check if user provided Reference exists in Reference Directory
    """
    file = pathToRef + "/refTable"

    refMatch = None        
    with open(file, 'r') as refFile:
        for line in refFile:
            line = line.rstrip()
            refMatch = re.match(line, userRef)
            if refMatch:
                return userRef
                
    if refMatch is None:
        return None      

def main():
    """
    main() 
    """
#******************************************************************************
# Command line args
#******************************************************************************    
    cmdparser = argparse.ArgumentParser(description="Create a strain specific gene sequence from a vcf file.",
                                        usage='%(prog)s -g <gene> -s <strain> -v <vcf file> [options]' ,prog='convertGene.py'  )
    cmdparser.add_argument('-g', '--gene',    action='store', dest='GENE',   help='Gene name: YAL046C or AIM1 supported.')
    cmdparser.add_argument('-s', '--strain',  action='store', dest='STR',    help='Strain file, comma separated on one line.')
    cmdparser.add_argument('-v', '--vcf',     action='store', dest='VCF',    help='VCF file name')
    cmdparser.add_argument('-r', '--ref',     action='store', dest='REF',    help='Genome Reference, default = R64-1-1 (R64-2-1 also available)', default='R64-1-1' )
    cmdparser.add_argument('-w', '--window',  action='store', dest="WIN",    help='Number Bases around gene to include in sequence. default = 0, Not supported for protein.', type=int, default=0)
    cmdparser.add_argument('-t', '--type',    action='store', dest="TYPE",   help='Output DNA or AA sequnece', default='dna' )
    cmdparser.add_argument('-p', '--print',   action='store_true', dest="PRINT",  help='Print Current list of supported strains.')
    cmdparser.add_argument('-d', '--detail',  action='store_true', dest="DETAIL", help='Print out detailed program information.')
             
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)

#******************************************************************************
# Check and define input variables
# Print each parameter used to screen to give the user a record.
#******************************************************************************    
    # Print strain list if requested and exit
    if cmdResults['PRINT']:
        printStrainList()
        sys.exit(1)
        
    # Print detailed program information to screen
    if cmdResults['DETAIL']:
        print("\tconvertGene.py program information:\n")
        print("  Purpose: Create a strain specific gene sequence by replacing snps")
        print("\tfound vcf into the reference gene sequence.")
        print("\n  Required parameters:")
        print("\t-g Gene Name in yeast ORF format or as standard SGD gene name")
        print("\texample: YAL046C = AIM1\n")
        print("\t-s List of strain names, this is a text file, with strain names comma delimited.")
        print("\texample line: K10,YJM653,YPS606,YS2 \n")
        print("\t-v VCF file provided by user, the strain names need to be cleaned up.")
        print("\t which means DBVPG1788 is valid, but /path/to/DBVPG1788/ is not.\n")
        print("\n  Optional parameters:")
        print("\t-r Genome Reference, default is R64-1-1, use -r R64-2-1 for latest reference\n")
        print("\t-w Window length in bps up & down stream to include with gene sequence.")
        print("\texample: -w 100 returns 100bp up and down stream of gene for DNA.\n")
        print("\tWindow not supported for proteins\n")
        print("\tWARNING: ")
        print("\tBe aware that any DNA codon with a IUPAC value that is not")
        print("\tstandard will be output as a nonstandard character.")
        print("\tThis will be caused by heterozygosity, a call of A/G will be output as R.\n")
        print("\t-t type of output, defaults to DNA, use aa if protein required.\n")
        print("\t-p y Print the currently supported strain list")
        print("\t\tSee Mike Place about adding more strains.\n")
        print("  Dependencies:")
        print("\tPython 3, VCFtools (v0.1.12b), gffutils")
        print("\tMust have vcftools installed and have the perl library in the path.\n")
        print("\tOne way to ensure vcftools lib is in the path:")
        print("\tAt the command line run:")
        print("\texport PERL5LIB=/FULLPATHTO/bin/vcftools_0.1.12b/lib/perl5/site_perl")
        print("\tThis is set in os.environ automatically in script.\n")
        print("\tgffutils information at https://pythonhosted.org/gffutils/index.html")
        print("\tyeast_Gene_name_to_ORF.py must be present (maps gene names to SGD ORF name) \n")
        print("\tExample run and screen output:\n")
        print("\t ~/scripts/convertGene.py -g YBL100W-C -s strains.input -v test.vcf\n")
        print("\tRunning convertGene.py ....")
        print("\tUsing the following parameters")
        print("\tReference    : R64-1-1")
        print("\tVCF File     : test.vcf")
        print("\tWindow size  : 0")
        print("\tStrain Names : DBVPG1788,K10,YJM653,YPS606,YS2,YS9")
        print("\tType         : dna")
        print("\tGene         : YBL100W-C")
        print("\tGFF information")
        print("\tGene Chrom   : II")
        print("\tGene Start   : 28427")
        print("\tGene End     : 28546")
        print("\tGene Strand  : +")
        print("\tExtracting reference sequence from fasta .....")
        print("\tRef File: /home/mplace/Projects/convertGene/reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa")
        print("")
        print("\tSometimes the following will print, it is not an error just output from vcftools.")
        print("\t\"Empty fields in the header line, the column 90 is empty, removing.\"")
        print("\n\tExpect 2 output files:")
        print("\t 1) YBL100W-C.reference.fasta - original gene reference sequence file")
        print("\t 2) Sequence_results.fasta    - contains all the strain specific sequences requested\n\n")

        sys.exit(1)
    
    # check if reference exists 
    if cmdResults['REF'] is not None:
        ref = cmdResults['REF'].rstrip()
        reference = verifyRef( ref, pathToRef )
        if reference is None:
            print("\n\tReference Name is invalid!\n")
            cmdparser.print_help()
            sys.exit(1)
    
    print("\n\tRunning convertGene.py ....")
    print("\tUsing the following parameters")  
    print("\tReference    : %s" %(reference))
    
    if cmdResults['VCF'] is not None:
        vcfFile = cmdResults['VCF']
        print("\tVCF File     : %s" %(vcfFile))    
    
    if cmdResults['WIN'] is not None:
        window = cmdResults['WIN']
        print("\tWindow size  : %d" %(window))
    
    # Strain file assumed to have only one line of comma separated values
    # no spaces at all
    if cmdResults['STR'] is not None:
        strainFile = cmdResults['STR']
        strainList = []
        with open(strainFile) as FILE:
            line = FILE.readline()
            line = line.rstrip()
            strainList = line.split(",")
        valid = verifyStrainName( strainList )
        if valid is not None:
            print("\n\tStrain name: %s Not in list" %(valid))
            print("\tDo you need to add new strain to list?")
            sys.exit(1)
        else:
            print("\tStrain Names : %s" %( ','.join(strainList) ))
    
    if cmdResults['TYPE'] is not None:
        outType = cmdResults['TYPE']
        outType = outType.lower()
        if outType == 'dna' or outType == 'aa':
            print("\tType         : %s" %(outType))
        else:
            print("\n\tInvalid output type!: %s" %(outType))
            sys.exit(1)
            
    # check if gene exists in GFF file
    if cmdResults['GENE'] is not None:
        gene = cmdResults['GENE']
        if not gene.startswith("Y"):                                           # if gene is not ORF format convert to it
            if gene in yeast_Gene_name_to_ORF.geneToOrf:
                orfgene = yeast_Gene_name_to_ORF.geneToOrf[gene]
                gene = orfgene
        geneInfo = searchGff( gene, refDict[reference], reference )            # geneInfo list has chrom, left, right positions and strand
        print("\tGene         : %s" %(gene)) 
            

    if len(geneInfo) != 0:
        print("\n\tGFF information")
        print("\tGene Chrom   : %s" %(geneInfo[0]))
        gffStart = int(geneInfo[1] + 1)
        print("\tGene Start   : %d" %(gffStart) )
        print("\tGene End     : %s" %(geneInfo[2]))
        print("\tGene Strand  : %s" %(geneInfo[3]))
    else:
        print("\n\tInvalid Gene Name!\n")
        print("\tGene Name must be SGD ORF")
        print("\texample: YAL046C")
        sys.exit(1)
        
    # EXTRACT REFERENCE GENOME SEQUENCE, pass in reference name, chrom, left & right positions 
    # if aa window is set to 0 
    if outType == 'aa':
        record = extractRefSeq( reference, geneInfo[0], geneInfo[1], geneInfo[2] , 0, geneInfo[3] )  # <class 'Bio.Seq.Seq'>
    else:
        record = extractRefSeq( reference, geneInfo[0], geneInfo[1], geneInfo[2] , window, geneInfo[3] )        
    # Write Original Gene ref sequence to fasta file for later use
    if outType == 'aa':
        newRecord = seqTranslate(record, geneInfo[0])
        geneRefSeqFile = writeRefFasta(newRecord, gene, geneInfo[0], geneInfo[1], geneInfo[2] )
    else:
        geneRefSeqFile = writeRefFasta(record, gene, geneInfo[0], geneInfo[1], geneInfo[2] )

    # CREATE TEMP vcf SUBSET files for each strain
    for strain in strainList:
        sName = strain + '.temp.vcf'
        os.system('vcf-subset -c %s %s > %s' %(strain, vcfFile, sName ))

        # create an strain specific gene sequence
        strainSeq = extractVCF( sName, record, geneInfo[0], geneInfo[1], geneInfo[2], window, strain )

        if os.path.exists(sName):
            os.remove(sName)
            
        # if amino acid sequence translate 
        if outType == 'aa':
            strainAA = seqTranslate(strainSeq, geneInfo[0])
            writeFasta( strainAA, strain, geneInfo[0], geneInfo[1], geneInfo[2] )
        else:
            writeFasta( strainSeq, strain, geneInfo[0], geneInfo[1], geneInfo[2] )
                

if __name__ == "__main__":
    main()

