
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
lines = GWAS_file.readlines()
GWAS_SNP = []
GWAS_Pval = []
count = 1
for line in lines:
    p = line.split()
    GWAS_SNP.append(p[0].split(":")[1])
    GWAS_Pval.append(p[1])
    print (count)
    count += 1

# Create a Pandas dataframe table
# this table have each row of GWASSNPs and index as their p values

GWAS_df = pd.DataFrame({"SNPs" : GWAS_SNP, "P Value" : GWAS_Pval})
bag_df.sort_values(by="SNPs")
GWAS_df

## create a dataframe for storing the mapped results
# Create a Pandas dataframe table
# this table have each row of SNPs and index as their bag index

GWAS_bag_SNP = []
GWAS_bag_idx = []
unmapped_GWAS_SNP = []

for i in range(GWAS_df.__len__()):
    for j in range(bag_df.__len__()):
        if GWAS_df["SNPs"][i] == bag_df["SNPs"][j]:
            GWAS_bag_SNP.append(GWAS_df["SNPs"][i])
            GWAS_bag_idx.append(bag_df["bag index"][j])
            break
        else:
            unmapped_GWAS_SNP.append(GWAS_df["SNPs"][i])
        print(j)

    print(i)

GWASbag_df = pd.DataFrame({"SNPs": GWAS_bag_SNP, "bag index": GWAS_bag_idx})
GWASbag_df.sort_values(by="bag index")
GWASbag_df

# Create a dataframe table for storing (bag index - representative p value) information

GWAS_bag_amounts = len(GWASbag_df["bag index"].unique())
Mapped_GWAS_SNP_amounts = len(GWASbag_df["SNPs"])
Missed_GWAS_SNP_amounts = len(GWAS_df["SNPs"]) - Mapped_GWAS_SNP_amounts
#merge GWASbag_df table with GWAS_df
GWASbag_df = pd.merge(GWAS_df, GWASbag_df, on="SNPs")
# for each bag index in merged GWASbag_df , find the smallest p value, assign to the new GWAS_bag_result_df dataframe
GWASbagidx = []
smallestPval = []
for idx in range(GWAS_bag_amounts):
    df = GWASbag_df["bag index"] == idx + 1
    smallestPval.append(df["P Value"].min())
    GWASbagidx.append(idx + 1)
GWAS_bag_result_df = pd.DataFrame({"bag index": GWASbagidx, "P Value": smallestPval})

# assign  smallest p value to each GWAS bag

bagSNP_file.close()
bagPval_file.close()
GWAS_file.close()
