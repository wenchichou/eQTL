#!/usr/bin/python
# use Python-3.4
# python /home/unix/wcchou/gsapWenChi/gautvik/bin/gwas_bagging/GWAS_enrichment.arg.py ../bagging/bagSNP.chr1.print.txt ../../data/gwas.2.5m/all.pvalue/HEIGHT_hg19_Pvalue 1
# headers
import sys
import pandas as pd
import numpy as np
import re
import os
import gzip
from pandas import DataFrame

## Problem shooting -- p value need to match bags, generate p value file

#read in  the SNP bag file
# wcc arg3 = chromosome number
# wcc arg1
#bagSNP_file = open("C:\\Users\\User\\Desktop\\eQTL_project\\bagging\\finished.Bagging.NoMissing\\bagSNP.chr1.print.txt", 'r').read()
#sys.argv[1] = "/home/simon/eQTL_project/GTeX_enrichment/eQTL_finshed_bagging/bagSNP.chr1.print.txt"
inputeQTLLD_filePath = sys.argv[1]
#sys.argv[2] = "/home/simon/eQTL_project/GTeX_enrichment/Muscle_Skeletal_Analysis_chrpos_pval.txt"
inputGWAS_filePath = sys.argv[2]
#sys.argv[3] = '1'
chromosome = sys.argv[3]
outputDirName = sys.argv[4]

bagSNP_file = open(inputeQTLLD_filePath, 'r').read()
#GWAS_file = gzip.open(inputGWAS_filePath, 'r')
print ("working on chromosome " + chromosome)
print ("input eQTL LD file: " + inputeQTLLD_filePath)
print ("input GTeX file: " + inputGWAS_filePath)
# mkdir a output dir
#outputDirName = "/restricted/projectnb/chengcancer/Simon/project/T2D/bulk_RNA-seq/merge_170415_NB501164_170506_NB501164/2mapping/QT/outdir"

#This original code would occur race conditional error sometimes
if not os.path.exists(outputDirName):
            os.makedirs(outputDirName)

#while True:
#     try:
#          os.makedirs(outputDirName)
#          break
#     except OSError, e:
#        if e.errno != 17:
#            raise   
        # time.sleep might help here
#        pass


print ("output Folder is at " + outputDirName)

outputBagSNPFilePath=outputDirName + '/GTeXbag.chr' + sys.argv[3] + '_bagSNP'
outputBagPvalueFilePath=outputDirName + '/GTeXbag.chr' + sys.argv[3] + '_bagPvalue'
print ("output bagSNP file is at " + outputBagSNPFilePath)
print ("output bagPvalue file is at " + outputBagPvalueFilePath)

# split the strings of the bagSNP by  spaces and newlines
print ("Loading eQTL bag ...")
bagSNP = []
bagidx = []
bag_number = 1
for line in range(len(bagSNP_file.splitlines())):
    temp = (bagSNP_file.splitlines()[line])
    words = temp.split()
    idx = 0
    for idx in range(len(words)):
        bagSNP.append(words[idx])
        bagidx.append(bag_number)
    #print (bag_number)
    bag_number += 1

print ("There are ", bag_number, " LD bags in the input eQTL LD file.")
# For chr1,  bqg_number should be 93418


# Create a Pandas dataframe table
# this table have each row of bagSNPs and index as their bag numbers
bag_df = pd.DataFrame({"SNPs" : bagSNP, "bag index" : bagidx})
#wcc020417 bag_df.sort_values(by="SNPs")
#bag_df

#import GWAS file
# wcc arg2
#GWAS_file = open("C:\\Users\\User\\Desktop\\eQTL_project\\python_data\\HEIGHT_hg19_Pvalue", 'r')
# read the GWAS line by line


GWAS_SNP = []
GWAS_Pval = []
count = 1
gwasFileLineCount = 0
    
#with open(inputGWAS_filePath) as GWAS_file:
with gzip.open(inputGWAS_filePath, 'rt') as GWAS_file:
  for line in GWAS_file:
    line = line.strip()
    gwasFileLineCount += 1
    p = line.split()
# wcc 1=arg3=chr## Need to change the number if we bag GWAS with different chromosome !!!
    #if(p[0].split(":")[0] == '1'):      
    if (p[0].split(":")[0] == chromosome):
        GWAS_SNP.append(p[0].split(":")[1])
        GWAS_Pval.append(p[1])
        count += 1
print ("There ", gwasFileLineCount, " SNPs in the input GWAS file.")
print ("There ", count, " GTeX SNPs in selected chromosome.")
# for HEIGHT GWAS, there should be 2481248 lines, for chromosome 1, there should be 188142 lines
# already check that it doesn't contain chr 1X in the list p

# Create a Pandas dataframe table
# this table have each row of GWASSNPs and index as their p values

GWAS_df = pd.DataFrame({"SNPs" : GWAS_SNP, "P Value" : GWAS_Pval})
#wcc020417 GWAS_df = GWAS_df.sort_values(by="SNPs")

## create a dataframe for storing the mapped results
# Create a Pandas dataframe table
# this table have each row of SNPs and index as their bag index
GWASbag_df = pd.merge(bag_df, GWAS_df, on="SNPs") # 161068 rows for chr1 and HEIGHT
#GWASbag_df  = GWASbag_df[~GWASbag_df["P Value"].str.contains("NA")] #160951
#GWASbag_df  = GWASbag_df[~GWASbag_df["P Value"].str.contains("NaN")]
GWASbag_df.dropna()  ## Need to reset the bag index
#wcc020417 GWASbag_df = GWASbag_df.sort_values(by = "bag index")
GWASbag_df = GWASbag_df.reset_index(drop=True)  # reset the index of rows, not index of bags


# Create a dataframe table for storing (bag index - representative p value) information

GWAS_bag_amounts = len(GWASbag_df["bag index"].unique()) # 44708  bags for chr1 and HEIGHT
Mapped_GWAS_SNP_amounts = len(GWASbag_df["SNPs"])
Missed_GWAS_SNP_amounts = len(GWAS_df["SNPs"]) - Mapped_GWAS_SNP_amounts
bag_number   # original eQTL bag numbers


# Solve the jump number of bag index problem
newbagidx = []
idx = 1
newbagidx.append(idx)
for j in range (1, len(GWASbag_df)): # len(GWASbag_df) = 161068 for HEIGHT chr1    ## 1~161067
    if (GWASbag_df["bag index"] [j] != GWASbag_df["bag index"] [j - 1]):
        idx += 1
        newbagidx.append(idx)
    else:
        newbagidx.append(idx)

GWASbag_df = pd.DataFrame({"SNPs" : GWASbag_df["SNPs"], "P Value" : GWASbag_df["P Value"], "bag index" : newbagidx})

# for each bag index in merged GWASbag_df , find the smallest p value, assign to the new GWAS_bag_result_df dataframe
GWASbagidx = []
smallestPval = []
for idx in range(GWAS_bag_amounts):
    df = (GWASbag_df[GWASbag_df["bag index"] == idx + 1])
    smallestPval.append(df["P Value"].min())
    GWASbagidx.append(idx + 1)
GWAS_bag_result_df = pd.DataFrame({"bag index": GWASbagidx, "P Value": smallestPval})
GWAS_bag_result_df

#  space delimited output


#  space delimited output
# wcc arg4 pvalue outputfile
#GWAS_bag_result_df['P Value'].to_csv('C:\\Users\\User\\Desktop\\eQTL_project\\python_data\\%s_bagPval_Chr1.txt'%"HEIGHT_hg19"                                     , header=None, index=None, mode='w')
GWAS_bag_result_df['P Value'].to_csv(outputBagPvalueFilePath, header=None, index=None, sep=' ', mode='w')




# wcc arg5 SNP outputfile"
#outFile = open('C:\\Users\\User\\Desktop\\eQTL_project\\python_data\\%s_bagSNP_Chr1.txt'% "HEIGHT_hg19", 'w')
outFile = open(outputBagSNPFilePath, 'w')

for bagidx in range(GWAS_bag_amounts):
    test = GWASbag_df[GWASbag_df["bag index"] == (bagidx + 1)].sort_values(by = "P Value")
    a = test.groupby("bag index")["SNPs"].apply(tuple)
    for result in a:
     result = ' '.join(result)
     outFile.write(result + '\n')

GWAS_file.close()
outFile.close()


