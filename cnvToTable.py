#!/home/mplace/anaconda3/bin/python
"""
@Program: cnvToTable.py

@Purpose: create table for graphing control freec's _CNV P-value file.
         
@Input:  Fasta file

@Dependencies: Python 3

@author: Mike Place

"""
import sys
    
chromArabicR64 = {  "ref|NC_001133|[R64]" : "1", "ref|NC_001134|[R64]" : "2", "ref|NC_001135|[R64]" : "3",
           "ref|NC_001136|[R64]" : "4", "ref|NC_001137|[R64]" : "5", "ref|NC_001138|[R64]" : "6", 
           "ref|NC_001139|[R64]" : "7", "ref|NC_001140|[R64]" : "8", "ref|NC_001141|[R64]" : "9",
           "ref|NC_001142|[R64]" : "10", "ref|NC_001143|[R64]" : "11", "ref|NC_001144|[R64]" : "12",
           "ref|NC_001145|[R64]" : "13", "ref|NC_001146|[R64]" : "14", "ref|NC_001147|[R64]" : "15",
           "ref|NC_001148|[R64]" : "16", "ref|NC_001224|[R64]" : "17" }  

chromArabic = {  "ref|NC_001133|" : "1", "ref|NC_001134|" : "2", "ref|NC_001135|" : "3",
           "ref|NC_001136|" : "4", "ref|NC_001137|" : "5", "ref|NC_001138|" : "6", 
           "ref|NC_001139|" : "7", "ref|NC_001140|" : "8", "ref|NC_001141|" : "9",
           "ref|NC_001142|" : "10", "ref|NC_001143|" : "11", "ref|NC_001144|" : "12",
           "ref|NC_001145|" : "13", "ref|NC_001146|" : "14", "ref|NC_001147|" : "15",
           "ref|NC_001148|" : "16", "ref|NC_001224|" : "17" }  

 
# lengths for Y22-3 reference fasta           
chromLen = { "ref|NC_001133|[R64]" : 214974, "ref|NC_001134|[R64]" : 826451, "ref|NC_001135|[R64]" : 316888,
            "ref|NC_001136|[R64]" : 1556355,"ref|NC_001137|[R64]" : 573699,"ref|NC_001138|[R64]" : 305960,
            "ref|NC_001139|[R64]" : 1119015,"ref|NC_001140|[R64]" : 545823,"ref|NC_001141|[R64]" : 461146,
            "ref|NC_001142|[R64]" : 764844,"ref|NC_001143|[R64]" : 690853,"ref|NC_001144|[R64]" : 1041051,
            "ref|NC_001145|[R64]" : 908319,"ref|NC_001146|[R64]" : 784798,"ref|NC_001147|[R64]" : 1063911,
            "ref|NC_001148|[R64]" : 910614,"ref|NC_001224|[R64]" : 85779 }
            
# offset for each chromosome, i.e. distance to add to start & stop gene positions            
offset64 = { "ref|NC_001133|[R64]" : 0, "ref|NC_001134|[R64]" : 214974, "ref|NC_001135|[R64]" : 1041425,
            "ref|NC_001136|[R64]" :1358313 ,"ref|NC_001137|[R64]" : 2914668,"ref|NC_001138|[R64]" : 3488367,
            "ref|NC_001139|[R64]" :3794327 ,"ref|NC_001140|[R64]" : 4913342 ,"ref|NC_001141|[R64]" :  5459165,
            "ref|NC_001142|[R64]" :5920311 ,"ref|NC_001143|[R64]" : 6685155,"ref|NC_001144|[R64]" : 7376008,
            "ref|NC_001145|[R64]" :8417059 ,"ref|NC_001146|[R64]" : 9325378,"ref|NC_001147|[R64]" : 10110176,
            "ref|NC_001148|[R64]" :11174087 ,"ref|NC_001224|[R64]" :12084701  }

# offset for each chromosome, i.e. distance to add to start & stop gene positions            
offset = { "ref|NC_001133|" : 0, "ref|NC_001134|" : 214974, "ref|NC_001135|" : 1041425,
            "ref|NC_001136|" :1358313 ,"ref|NC_001137|" : 2914668,"ref|NC_001138|" : 3488367,
            "ref|NC_001139|" :3794327 ,"ref|NC_001140|" : 4913342 ,"ref|NC_001141|" :  5459165,
            "ref|NC_001142|" :5920311 ,"ref|NC_001143|" : 6685155,"ref|NC_001144|" : 7376008,
            "ref|NC_001145|" :8417059 ,"ref|NC_001146|" : 9325378,"ref|NC_001147|" : 10110176,
            "ref|NC_001148|" :11174087 ,"ref|NC_001224|" :12084701  }      


def main():
    """
    main() 
    """
  
    # if no args print help
    if len(sys.argv) == 1:
        print("\n\t_CNV input file required.")
        print("\n")
        print("\tcnvToTable.py <_CNV file>")
        print("\tCreate a table from Control Freec's _CNV file.")
        sys.exit(1)
    else:
        cnvFile = sys.argv[1]
        
    # print header
    print("chr\tstart\tend\tcopynum\tstatus\tWilcoxonRankSumTestPvalue\tKolmogorovSmirnovPvalue")
    with open(cnvFile) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#'):
                continue
            row  = line.split()
            
            if row[0] in offset64:
                row[1] = str(int( row[1]) + int( offset64[row[0]] ) )
                row[2] = str(int( row[2]) + int( offset64[row[0]] ))
                row[0] = chromArabicR64[row[0]]
                print( "\t".join(row) )

               

if __name__ == "__main__":
    main()

