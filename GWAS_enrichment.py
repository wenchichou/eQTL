
# headers
import pandas as pd
import numpy as np
from pandas import DataFrame

## Problem shooting -- p value need to match bags, generate p value file

#read in  the SNP bag file
bagSNP_file = open("C:\\Users\\User\\Desktop\\eQTL_project\\bagging\\finished.Bagging.NoMissing\\bagSNP.chr1.print.txt", 'r').read()

# split the strings of the bagSNP by  spaces and newlines
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
    print (bag_number)
    bag_number += 1
# For chr1,  bqg_number should be 93418


# Create a Pandas dataframe table
# this table have each row of bagSNPs and index as their bag numbers
bag_df = pd.DataFrame({"SNPs" : bagSNP, "bag index" : bagidx})
bag_df.sort_values(by="SNPs")
bag_df

#import GWAS file

GWAS_file = open("C:\\Users\\User\\Desktop\\eQTL_project\\python_data\\HEIGHT_hg19_Pvalue", 'r')
# read the GWAS line by line
for line in GWAS_file:
    lines = GWAS_file.readlines()

GWAS_SNP = []
GWAS_Pval = []
count = 1
for line in lines:
    p = line.split()
    if (p[0].split(":")[0] == '1'):      ## Need to change the number if we bag GWAS with different chromosome !!!
            GWAS_SNP.append(p[0].split(":")[1])
            GWAS_Pval.append(p[1])
            count += 1
# for HEIGHT GWAS, there should be 2481248 lines, for chromosome 1, there should be 188142 lines
# already check that it doesn't contain chr 1X in the list p

# Create a Pandas dataframe table
# this table have each row of GWASSNPs and index as their p values

GWAS_df = pd.DataFrame({"SNPs" : GWAS_SNP, "P Value" : GWAS_Pval})
GWAS_df = GWAS_df.sort_values(by="SNPs")

## create a dataframe for storing the mapped results
# Create a Pandas dataframe table
# this table have each row of SNPs and index as their bag index
GWASbag_df = pd.merge(bag_df, GWAS_df, on="SNPs") # 161068 rows for chr1 and HEIGHT
#GWASbag_df  = GWASbag_df[~GWASbag_df["P Value"].str.contains("NA")] #160951
#GWASbag_df  = GWASbag_df[~GWASbag_df["P Value"].str.contains("NaN")]
GWASbag_df.dropna()  ## Need to reset the bag index
GWASbag_df = GWASbag_df.sort_values(by = "bag index")
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
GWAS_bag_result_df['P Value'].to_csv('C:\\Users\\User\\Desktop\\eQTL_project\\python_data\\%s_bagPval_Chr1.txt'%"HEIGHT_hg19"
                                     , header=None, index=None, mode='w')


outFile = open('C:\\Users\\User\\Desktop\\eQTL_project\\python_data\\%s_bagSNP_Chr1.txt'% "HEIGHT_hg19", 'w')

for bagidx in range(GWAS_bag_amounts):
    test = GWASbag_df[GWASbag_df["bag index"] == (bagidx + 1)].sort_values(by = "P Value")
    a = test.groupby("bag index")["SNPs"].apply(tuple)
    for result in a:
     result = ' '.join(result)
     outFile.write(result + '\n')

GWAS_file.close()
outFile.close()


