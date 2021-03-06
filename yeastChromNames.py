# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 08:58:26 2015

@author: mplace
"""

chromSGD = { "chrI" : "ref|NC_001133|", "chrII" : "ref|NC_001134|", "chrIII" : "ref|NC_001135|",
              "chrIV" : "ref|NC_001136|", "chrV" : "ref|NC_001137|", "chrVI" : "ref|NC_001138|", 
              "chrVII" : "ref|NC_001139|", "chrVIII" : "ref|NC_001140|", "chrIX" : "ref|NC_001141|",
              "chrX" : "ref|NC_001142|", "chrXI" : "ref|NC_001143|", "chrXII" : "ref|NC_001144|",
              "chrXIII" : "ref|NC_001145|", "chrXIV" : "ref|NC_001146|", "chrXV" : "ref|NC_001147|",
              "chrXVI" : "ref|NC_001148|", "chrMito" : "ref|NC_001224|", "chrmt" : "ref|NC_001224|"
} 

chromArabic = {  "ref|NC_001133|" : "chr1", "ref|NC_001134|" : "chr2", "ref|NC_001135|" : "chr3",
           "ref|NC_001136|" : "chr4", "ref|NC_001137|" : "chr5", "ref|NC_001138|" : "chr6", 
           "ref|NC_001139|" : "chr7", "ref|NC_001140|" : "chr8", "ref|NC_001141|" : "chr9",
           "ref|NC_001142|" : "chr10", "ref|NC_001143|" : "chr11", "ref|NC_001144|" : "chr12",
           "ref|NC_001145|" : "chr13", "ref|NC_001146|" : "chr14", "ref|NC_001147|" : "chr15",
           "ref|NC_001148|" : "chr16", "ref|NC_001224|" : "chr17" } 

chromArabic64 = {  "ref|NC_001133|[R64]" : "chr1", "ref|NC_001134|[R64]" : "chr2", "ref|NC_001135|[R64]" : "chr3",
           "ref|NC_001136|[R64]" : "chr4", "ref|NC_001137|[R64]" : "chr5", "ref|NC_001138|[R64]" : "chr6", 
           "ref|NC_001139|[R64]" : "chr7", "ref|NC_001140|[R64]" : "chr8", "ref|NC_001141|[R64]" : "chr9",
           "ref|NC_001142|[R64]" : "chr10", "ref|NC_001143|[R64]" : "chr11", "ref|NC_001144|[R64]" : "chr12",
           "ref|NC_001145|[R64]" : "chr13", "ref|NC_001146|[R64]" : "chr14", "ref|NC_001147|[R64]" : "chr15",
           "ref|NC_001148|[R64]" : "chr16", "ref|NC_001224|[R64]" : "chr17" }
           
chromRoman = {  "ref|NC_001133|" : "chrI", "ref|NC_001134|" : "chrII", "ref|NC_001135|" : "chrIII",
           "ref|NC_001136|" : "chrIV", "ref|NC_001137|" : "chrV", "ref|NC_001138|" : "chrVI", 
           "ref|NC_001139|" : "chrVII", "ref|NC_001140|" : "chrVIII", "ref|NC_001141|" : "chrIX",
           "ref|NC_001142|" : "chrX", "ref|NC_001143|" : "chrXI", "ref|NC_001144|" : "chrXII",
           "ref|NC_001145|" : "chrXIII", "ref|NC_001146|" : "chrXIV", "ref|NC_001147|" : "chrXV",
           "ref|NC_001148|" : "chrXVI", "ref|NC_001224|" : "chrMito" } 
           
chromRoman64 = {  "ref|NC_001133|[R64]" : "chrI", "ref|NC_001134|[R64]" : "chrII", "ref|NC_001135|[R64]" : "chrIII",
           "ref|NC_001136|[R64]" : "chrIV", "ref|NC_001137|[R64]" : "chrV", "ref|NC_001138|[R64]" : "chrVI", 
           "ref|NC_001139|[R64]" : "chrVII", "ref|NC_001140|[R64]" : "chrVIII", "ref|NC_001141|[R64]" : "chrIX",
           "ref|NC_001142|[R64]" : "chrX", "ref|NC_001143|[R64]" : "chrXI", "ref|NC_001144|[R64]" : "chrXII",
           "ref|NC_001145|[R64]" : "chrXIII", "ref|NC_001146|[R64]" : "chrXIV", "ref|NC_001147|[R64]" : "chrXV",
           "ref|NC_001148|[R64]" : "chrXVI", "ref|NC_001224|[R64]" : "chrXVII" }
           
chromNum =  { "1" : "ref|NC_001133|", "2" : "ref|NC_001134|", "3" : "ref|NC_001135|",
              "4" : "ref|NC_001136|", "5" : "ref|NC_001137|", "6" : "ref|NC_001138|", 
              "7" : "ref|NC_001139|", "8" : "ref|NC_001140|", "9" : "ref|NC_001141|",
              "10" : "ref|NC_001142|", "11" : "ref|NC_001143|", "12" : "ref|NC_001144|",
              "13" : "ref|NC_001145|", "14" : "ref|NC_001146|", "15" : "ref|NC_001147|",
              "16" : "ref|NC_001148|", "mito" : "ref|NC_001224|"}      

chromNum64 =  { "1" : "ref|NC_001133|[R64]", "2" : "ref|NC_001134|[R64]", "3" : "ref|NC_001135|[R64]",
              "4" : "ref|NC_001136|[R64]", "5" : "ref|NC_001137|[R64]", "6" : "ref|NC_001138|[R64]", 
              "7" : "ref|NC_001139|[R64]", "8" : "ref|NC_001140|[R64]", "9" : "ref|NC_001141|[R64]",
              "10" : "ref|NC_001142|[R64]", "11" : "ref|NC_001143|[R64]", "12" : "ref|NC_001144|[R64]",
              "13" : "ref|NC_001145|[R64]", "14" : "ref|NC_001146|[R64]", "15" : "ref|NC_001147|[R64]",
              "16" : "ref|NC_001148|[R64]", "mito" : "ref|NC_001224|[R64]"} 
