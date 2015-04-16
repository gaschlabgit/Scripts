#!/home/mplace/anaconda/bin/python2.7
"""
PROGRAM: yeastmine.py
PURPOSE: Accept a file with a list of gene names in long SGD Name form (like YPR187W)
         and queries the database producing a user friendly table of required information.
         
INPUT: Plain text file with one gene name per line.
        YPR187W
        YPR202W

OUTPUT: csv text file

AUTHOR: Mike Place

DATE: 1/30/2015

"""

import argparse
from intermine.webservice import Service
service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")

def main():
    """
    
    """
    names  = []              # list of gene names
    result = []

    # HANDLE COMMAND LINE ARGS
    cmdparser = argparse.ArgumentParser( description="Query SGD Yeastmine with a list of gene names.", prog='yeastmain.py' )
    cmdparser.add_argument( '-f', '--file', action='store', required='true', 
                           dest='FILE', help='REQUIRED, text file with list of gene names' )
    
    cmdResults = vars( cmdparser.parse_args() )
    
    if cmdResults['FILE'] is not None:
        inFile = cmdResults['FILE']   
    
    with open( inFile ) as f:
        for line in f:
            ln = line.strip('\r\n\s')
            names.append(ln)
    
    query = service.new_query("Gene")
    query.add_view(  "secondaryIdentifier", "symbol", "name",
                   "organism.shortName", "proteins.symbol", "sgdAlias", "featureType", "description" )
    
    query.add_sort_order("Gene.secondaryIdentifier", "ASC")
    
    print "secondaryIdentifier\tsymbol\tname\torganism.shortName\tproteins.symbol\tsgdAlias\tfeatureType\tdescription"
    
    for row in query.rows():
        if row["secondaryIdentifier"] in names:
            data = row["secondaryIdentifier"], row["symbol"], row["name"], row["organism.shortName"],\
            row["proteins.symbol"], row["sgdAlias"], row["featureType"], row["description"]
            item = list(data)
            result.append(item)
            
    for d in result:
        for i in d:
            print "%s\t" %(i),
        print ""
    
    
    
        


if __name__ == "__main__":
    main()