
# headers
import pandas as pd
import numpy as np
from pandas import DataFrame

#read in  the SNP bag file
bagSNP_file = open("C:\\Users\\User\\Desktop\\eQTL_project\\python_data\\bagSNP.chr22.print.txt", 'r').read()

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
    if (p[0].split(":")[0] == '22'):
            GWAS_SNP.append(p[0].split(":")[1])
            GWAS_Pval.append(p[1])
            print (count)
            count += 1

# Create a Pandas dataframe table
# this table have each row of GWASSNPs and index as their p values

GWAS_df = pd.DataFrame({"SNPs" : GWAS_SNP, "P Value" : GWAS_Pval})
bag_df = bag_df.sort_values(by="SNPs")
GWAS_df = GWAS_df.sort_values(by="SNPs")

## create a dataframe for storing the mapped results
# Create a Pandas dataframe table
# this table have each row of SNPs and index as their bag index
GWASbag_df = pd.merge(bag_df, GWAS_df, on="SNPs")

# Create a dataframe table for storing (bag index - representative p value) information

GWAS_bag_amounts = len(GWASbag_df["bag index"].unique())
Mapped_GWAS_SNP_amounts = len(GWASbag_df["SNPs"])
Missed_GWAS_SNP_amounts = len(GWAS_df["SNPs"]) - Mapped_GWAS_SNP_amounts
bag_number   # original eQTL bag numbers
# for each bag index in merged GWASbag_df , find the smallest p value, assign to the new GWAS_bag_result_df dataframe
GWASbagidx = []
smallestPval = []
for idx in range(GWAS_bag_amounts):
    df = (GWASbag_df[GWASbag_df["bag index"] == idx + 1])
    smallestPval.append(df["P Value"].min())
    GWASbagidx.append(idx + 1)
GWAS_bag_result_df = pd.DataFrame({"bag index": GWASbagidx, "P Value": smallestPval})
GWAS_bag_result_df
GWAS_bag_result_df = GWAS_bag_result_df[pd.notnull(GWAS_bag_result_df['P Value'])]

#  space delimited output
GWAS_bag_result_df['P Value'].to_csv(r'C:\\Users\\User\\Desktop\\eQTL_project\\python_data\\HEIGHT_hg19_bagged_Pval.txt'
                                     , header=None, index=None, sep=' ', mode='a')

GWAS_file.close()
