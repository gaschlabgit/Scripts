#!/home/mplace/anaconda/bin/python2.7
"""
PROGRAM: yeastmine.py

PURPOSE: Accept a file with a list of gene names in found in SGD.
         Run queries via Yeastmine retrieving a table of gene interactions.
         Will create 2 levels of interactions, the first gene's interactions
         and the interactions with a list of secondary actors (i.e. Trey's genes).
         
INPUT:  Plain text file with one gene name per line.
        azf1
        act1
        ....

OUTPUT: Produces a tab delimited text file for each gene.
        level1 table info for primary gene
        level2 table info for secondary genes
        
Yeastmine:  http://yeastmine.yeastgenome.org/yeastmine/begin.do

AUTHOR: Mike Place

DATE: 7/16/2015
"""
import argparse
import os
import re
import sys
from intermine.webservice import Service
service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")

class interact( object ):
    """
    Automate the querying of Yeast mine using a list of genes.
    file = input file 
    dir = current working directory
    geneList = List of genes for level1 query.
    """
    def __init__( self, inFile ):
        """
        Initialize interact object
        """
        self.file = inFile
        self.dir  = os.getcwd()
        self.geneList = []    # list of genes to query
        self.geneInteractions = [ 'ISU1', 'HOG1', 'GSH1', 'GRE3', 'IRA2', 'SAP190' ]
        
        with open( self.file ) as f:
            for line in f:
                ln = line.strip('\r\n\s')
                self.geneList.append(ln)
        
    def callQueryYM( self ):
        """
        Run through gene name list to query interactions from yeastmine.
        """
        count = {}
        
        for item in self.geneList:
            result = self.queryYM(item, level="level1" )
        
        for d in result:
            if re.findall('\\b'+"level"+'\\b',d[0]):
                continue
            else:
                if d[11] in self.geneInteractions:
                    print d[11]
                    if d[11] in count:
                        continue
                    else:
                        count[d[11]] = 1
                        result += self.queryYM( d[11], level="level2")
                #    result += resultlevel2
        header =  "level", "primaryIdentifier", "symbol", "secondaryIdentifier", "sgdAlias", "name",\
               "organism.shortName", "interactions.details.annotationType",\
               "interactions.details.phenotype", "interactions.details.role1", \
               "interactions.details.experimentType", "interactions.participant2.symbol",\
               "interactions.participant2.secondaryIdentifier", "interactions.details.experiment.name", \
               "interactions.details.relationshipType"
        headerList = list(header)
        
        for line in result:
            print line
      
    
    def queryYM( self, geneName, level ):
        """
        Use Yeastmine API to look for gene interactions.
        """
        result = []        
            
    # Get a new query on the table to query
        query = service.new_query("Gene")
        query.add_constraint("interactions.participant2", "Gene")    
           
        query.add_view(   "primaryIdentifier", "symbol", "secondaryIdentifier", "sgdAlias", "name",
                       "organism.shortName", "interactions.details.annotationType",
                       "interactions.details.phenotype", "interactions.details.role1",
                       "interactions.details.experimentType", "interactions.participant2.symbol",
                       "interactions.participant2.secondaryIdentifier",
                       "interactions.details.experiment.name",
                       "interactions.details.relationshipType"    )
    
        query.add_constraint("organism.shortName", "=", "S. cerevisiae", code = "B")
        query.add_constraint("Gene", "LOOKUP", geneName, code = "A")
    
        for row in query.rows():
            data = level, row["primaryIdentifier"], row["symbol"], row["secondaryIdentifier"], row["sgdAlias"], \
                row["name"], row["organism.shortName"], row["interactions.details.annotationType"], \
                row["interactions.details.phenotype"], row["interactions.details.role1"], \
                row["interactions.details.experimentType"], row["interactions.participant2.symbol"], \
                row["interactions.participant2.secondaryIdentifier"], \
                row["interactions.details.experiment.name"], row["interactions.details.relationshipType"]
            item = list(data)
            result.append(item)    
        
        return result
    
    def writeTable( self, level ):
        """
        Write Yeastmine gene interaction query results to text file.
        level = indicates the depth of search
        """
        pass
    
    def findInteractions( targetGenes, geneTable ):
        """
        targetGenes = secondary genes used to search list of genes 
                      default are Trey's genes.
        geneTable really a list of gene names returned from last
        query of Yeastmine.
        """
        pass
    
    def getGeneList( self ):
        """
        Return initial gene list
        """
        return self.geneList
    
    def getDir( self ):
        """
        Return current working directory
        """
        return self.dir
    
    def getInFile( self ):
        """
        Return input file name
        """
        return self.file
    

def main():
    """
    
    """
    # HANDLE COMMAND LINE ARGS
    cmdparser = argparse.ArgumentParser( description="Query SGD Yeastmine with a list of gene names.", prog='yeastmain.py' )
    cmdparser.add_argument( '-f', '--file', action='store', dest='FILE', help='REQUIRED, text file with list of gene names' , metavar='')
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='Print a more detailed description of program.')
    
    cmdResults = vars( cmdparser.parse_args() )
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
    
    if cmdResults['INFO']:
       print "\n  yeastmine.py "
       print "\n  Purpose: Accept a file with a list of gene names in found in SGD.\n\
         Run queries via Yeastmine retrieving a table of gene interactions.\n\
         Will create 2 levels of interactions, the first gene's interactions\n\
         and the interactions with a list of secondary actors (i.e. Trey's genes)."

       print "\n  Input  : gene name list file, one gene per line"
       print "\n  Output : Table (tab delimited) for each gene in gene list input file"
       print "       level indicates the relationship depth, level1 = primary interaction"
       print "       returned by SGD."
       print "       level2 genes that are related to level1 via Trey's genes"
       print " "
       print " Trey's genes are: ISU1, HOG1, GSH1, GRE3, IRA2, SAP190 "       
        
       print "\n  Usage  : yeastmine.py -f genelist.txt"
       print "  "       
       print "\tTo see Python Docs for this program:"
       print "\tOpen python (version 2.7) console and enter"
       print " import sys"
       print " sys.path.append('/path/to/yeastmine.py') "
       print " import yeastmine"
       print " help(yeastmine)\n"
       sys.exit(1)
   
   
   # check that input file exists
    if cmdResults['FILE'] is not None:
        inFile = cmdResults['FILE'] 
        if not os.path.exists(inFile):
            print "\n\t-f input file does not exist.\n"
            cmdparser.print_help()
            sys.exit(1)
    
    data = interact( inFile )
    
    data.callQueryYM(  )
            


    
   
    
    
        


if __name__ == "__main__":
    main()