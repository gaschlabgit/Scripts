#! /bin/bash 
#  Shell script for running GATK unifiedgenotyper 
#  Katie Clowers/ Maria Sardi
#  Altered to run only unifiedgenotyper and downstream vcftools 
#  modified by Mike Place
#  April 2015
#
#  Assumes a high memory machine, set to 32Gb ram and 8 threads for unified genotyper
#******************************************************************************


#  Set variables
#BAM=$1  #path to BAM file don't include .bam Test! Need to create script to add this sequentially
#OUT=$2 #path to ref fasta
#MBQ=$3 #path to SNPdb
#
# Reference is SGD standard R64-1-1 
REF='/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/gatkRef/S288C_reference_sequence_R64-1-1_20110203.fasta'
# FILE = file w/ a list of *.Recalibrated.bam files (one per line ) to process with unifiedgenotyper 
FILE=$1
# INPUT = string of all files to process
INPUT=""        

display_usage() {
	echo ""
	echo " Missing arguments"
	echo " 2 required arguments"
	echo "  Usage: UNIFIEDGENETYPER.sh <inputfile> <outputPrefix> [MBQ value optional]"
	echo "" 
}

echo "$#"
echo "$FILE"

# if user provides < 2 arguments
if [ $# -le 1 ]
then
	display_usage
	exit 1
fi

# if user wants help
if [[ ( $# == "--help") ||  $# == "-h" ]] 
then 
	display_usage
	exit 0
fi 

OUT=$2

# check if input file exists
if [ -f $FILE ]
then
	for line in $(cat $FILE)
		do 
			INPUT+="-I "$line" "
		done
else
	echo "Input file doesn't exist"
fi

echo "Files used: "
echo $INPUT

# Check if MBQ parameter is present, if not use default
if [ -z "${3}" ]; then
	MBQ='-mbq 25'
else
	MBQ="-mbq $3"
fi

GATK=/opt/bifxapps/gatk3/GenomeAnalysisTK.jar # Path to GATK
echo "MBQ $MBQ"

#  Java parameters
JAVAPARAM='-Xmx32g'

#  Parameters for GATK
#  -mbq 30 (default value 17) Minimum base quality required to consider a base for calling
#  -stand_call_conf 50 (default value 30) The minimum phred-scaled confidence threshold for genotype calling
#  -glm (default value SNP) Should we call indels (INDEL), SNPs (SNP), or both (BOTH)?
#  -out_mode (default value EMIT_VARIANTS_ONLY) Should we output confident genotypes (EMIT_ALL_CONFIDENT_SITES) or just the variants?
#  -hets (default value 0.001) The expected heterozygosity used to compute prior likelihoods for genotype calling

#  UnifiedGenotyper parameters
# SNPPARAM='-mbq 30'
SNPPARAM=$MQB

#  Ploidy parameter
PLOIDY= '-ploidy 2'

#  VCFToTable parameters
TABPARAM='-F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F DP -F GQ -F GT -AMD'

############################################################################################################################################

#  5. Do single sample variant calling using UnifiedGenotyper
#  Can also do multi-sample variant calling after multiple files have been recalibrated. Maybe write a script for this?!?!?
echo "Step 5: Unified Genotyper"
echo
java $JAVAPARAM -jar $GATK -R $REF -T UnifiedGenotyper -glm BOTH -o $OUT.RawVariants.vcf -I $INPUT $SNPPARAM $PLOIDY -nt 8
echo


############################################################################################################################################ 
#  6. Apply hard filters to a call set
#  Apply hard filters to a variant callset that is too small for VQSR or for which truth/training sets are not available

echo "Step 6: Apply filters to call set"
echo
#6.1 Extract SNPs from the call set
#  Creates a VCF file called $BAM.RawSNPs.vcf containing just the SNPs from the original file of raw variants
echo "6.1 Extract SNPs from call set"
echo
java $JAVAPARAM -jar $GATK -T SelectVariants -R $REF -V $OUT.RawVariants.vcf -selectType SNP -o $OUT.RawSNPs.vcf $PLOIDY
echo

#  6.2 Determine parameters for filtering SNPs
#  SNPs matching any of these conditions will be considered bad and filtered out (marked FILTER in the output vcf).  Specifies which
#  parameter was chiefly responsible for exclusion of the SNP using the culprit annotation.  SNPs that do not match any of these conditions
#  will be considered good and marked PASS in the output vcf.
#  QualByDepths (QD) 2.0- variant confidences divided by unfiltered depth
#  FisherStrand (FS) 60.0- Phred-caled p-value using fisher's exact test to detect strand bias... more bias->false positives
#  RMSMappingQuality (MQ) 40.0- root mean square of the mapping quality of the reads across all samples
#  HaplotypeScore 13.0- consistency of the site with two (and only two) segregating haplotypes (not applicable for calls made using
#   UnifiedGenotyper on non-diploid organisms.
#  MappingQualitRankSumTest (MQRankSum) 12.5- u-based z-approx from the Mann-Whitney Rank Sum Test for mapping qualities (reads with ref
#   bases vs those with alt allele).  Can not be calculated for sites without a mixture of reads showing both the ref and alt alleles ie
#   this will only be applied to heterozygous calls.
#  ReadPosRankSumTest (ReadPosRankSum) 8.0- u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of
#  	the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, indicative of error.  Can
#   not be calculated for sites without a mixture of reads showing both the reference and alternate alleles ie this will only be applied
#   to heterozygous calls.


#  6.3 Apply the filter to the SNP call set
#  Creates a VCF file called $BAM.FiltSNPs.vcf containing all original SNPs, but annotated with either PASS or FILTER depending on whether
#  or not they passed the filters.  For SNPs that fail, the variant annotation also includes the name of the filter.  That way, if you apply
#  several different filters, you can keep track of which filter(s) each SNP failed, and later you can retireive specific subsets of your
#  calls using the SelectVariants tool.

echo "6.3 Apply filter to SNP call set"
echo
java $JAVAPARAM -jar $GATK -T VariantFiltration -R $REF -V $OUT.RawSNPs.vcf -o $OUT.FiltSNPs.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0 || HaplotypeScore > 13.0" --filterName "AllFilters" $PLOIDY
#  Add if diploid:  || HaplotypeScore > 13.0

echo
echo "Filtered SNPs in *.FiltSNPs.vcf file"
echo

#  6.4 Extract the indels from the call set
#  Creates a VCF file called $BAM.RawIndels.vcf containg just the indels from original file of raw variants.

echo "6.4 Extract Indels from call set"
echo
java $JAVAPARAM -jar $GATK -T SelectVariants -R $REF -V $OUT.RawVariants.vcf -selectType INDEL -o $OUT.RawIndels.vcf $PLOIDY
echo


#  6.5 Determine parameters for filtering indels
#  Indels matching any of these coniditions will be considered bad and filtered out (marked FILTER in the output vcf).  Specifies which
#  parameter was chiefly responsible for exclusion of the indel using the culprit annotation.  Indels that do not match any of these
#  conditions will be considered good and marked PASS in the output vcf.
#  QualByDepth (QD) 2.0- variant confidence (from the QUAL field) divided by the unfiltered depth of non-ref samples
#  FisherStrand (FS) 200.0- Phred-caled p-value using fisher's exact test to detect strand bias... more bias->false positives
#  ReadPosRankSumTest (ReadPosRankSum) 20.0- 8.0- u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end
#   of the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, indicative of error.  Can
#   not be calculated for sites without a mixture of reads showing both the reference and alternate alleles ie this will only be applied
#   to heterozygous calls.


#  6.6 Apply the filter to the indel call set
#  Creates a VCF file called $BAM.FiltIndels.vcf containing all the original indels but they are now annotated with either PASS or FILTER
#  depending on whether or not they passed the filters.  For indels that failed the filter, the variant annotation also includes the name
#  of the filter.  That way, if you apply several different filteres, you can keep track of which filter(s) each indel failed, and later you
#  can retireve specific subsets of your calls using the SelectVariants tool echo "6.6 Apply filter to indel call set"

echo
java $JAVAPARAM -jar $GATK -T VariantFiltration -R $REF -V $OUT.RawIndels.vcf -o $OUT.FiltIndels.vcf \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "AllFilters" $PLOIDY
echo
echo "Filtered indels in *.FiltIndels.vcf file"
echo

############################################################################################################################################
# 7 Create an additional vcf containing variants that passed all filter

echo "7.1 Keep SNPs that passed all filters"
echo
vcftools --vcf $OUT.FiltSNPs.vcf --keep-filtered PASS --out $OUT.PASS_SNPs --recode
echo

echo
echo "7.2 Keep Indels that passed all filters"
vcftools --vcf $OUT.FiltIndels.vcf --keep-filtered PASS --out $OUT.PASS_Indels --recode
echo
echo

############################################################################################################################################

#  8. Combine SNP and indel vcfs?
#  I don't understand all of the parameters here...
#java $JAVAPARAM -jar $GATK  -T CombineVariants --variant $BAM.FiltSNPs.vcf --variant $BAM.FiltIndels.vcf -o $BAM.FiltVariants.vcf

############################################################################################################################################

