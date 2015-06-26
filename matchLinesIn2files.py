#!/home/mplace/anaconda3/bin/python
"""
Created on Wed Mar 11 12:29:53 2015

This is a quick and dirty script which compares  the first column of 2 files.
Prints out lines that match

@author: mplace
"""
import sys
import re

def main():
    if len(sys.argv) == 1:
        print("\n\tInput file required.")
        print("\trun: matchLineIn2files.py file1 file2")
        print("\n")
        sys.exit(1)
    else:
         file1  = sys.argv[1]  # smaller file 1 column ex.  NC_001133_100000
         file2  = sys.argv[2]
        
    found = {}
    
    with open(file1) as f:
        for line in f:
            line = line.rstrip()
            data = line.split()
            found[data[0]] = 1
    
    with open(file2) as b:
        for line in b:
            line = line.rstrip()
            item = line.split()
            if item[0] in found:
                #print(line)
                found[item[0]] = line
                
    for k,v in found.items():
        print(k, v)

            
                    
if __name__ == "__main__":
    main()                    