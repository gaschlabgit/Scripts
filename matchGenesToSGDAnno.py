#!/home/mplace/anaconda3/bin/python
"""
Created on Wed Mar 11 12:29:53 2015

Program: match.py

Purpose: Match 1st column in SGD-summary file with the 7th column in 
         GFF-Gene-results.txt file. 

         example of files:

    SGD-summary:
         
         Q0045  S000007260 COX1 Cytochrome c OXidase 12884 ORF Q0045 S. cerevisiae
         ....
    GFF-Gene-results.txt:
        
        Y127-bam_CNVs   ref|NC_001134|[R64]     0       14000   0       loss    GTK30B00110     gene    6528    6560
        Y127-bam_CNVs   ref|NC_001134|[R64]     0       14000   0       loss    YBL113W-A       gene    3391    3873
        .....
        

Input: 2 text files delimited by white space.

Output: 2 text files, one with lines that match, the second with the lines
        from both files that didn't match.

@author: mplace
"""
import sys

def main():
    if len(sys.argv) == 1:
        print("\n\t2 Input files required.")
        print("\trun: match.py summary.file GFF-Gene-results.txt\n")
        print("\n")
        sys.exit(1)
    else:
         file1  = sys.argv[1]  # smaller file 1 column ex.  NC_001133_100000
         file2  = sys.argv[2]
         
    
    found = {}
    
    #with open("matched.txt", 'w') as match:
    matchOut = open("match.txt", 'w')
    notMatchedOut = open("not-matched.txt", 'w')
   
   # fileHandle.write("file\t%s\n" %("\t".join(rnaNames)))    
        
    # create a unique list of values in first file, gene name
    # something like: YGL258W
    with open(file1) as f:
        for line in f:                
            line = line.rstrip()
            data = line.split()
            found[data[0]] = data        
    
    with open(file2) as b:
        for line in b:
            line = line.rstrip()
            item = line.split()
                        
            if item[6] in found:
                tmp = found[item[6]]
                tmp[2:] = [' '.join(tmp[2:])]
                out = item[1:2] + item[4:5] + item[5:6] + tmp
                matchOut.write( "%s\n" %("\t".join(out))) 
            else:
                notMatchedOut.write(line)
                
    #for k,v in found.items():
    #    print(k, v)

            
                    
if __name__ == "__main__":
    main()                    