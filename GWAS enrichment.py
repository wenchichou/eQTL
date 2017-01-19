
# Try to finish the first 3 steps by Thursday Taiwan Time.
##1.1  import eQTL bags in python panda dataframe
import pandas as pd
from pandas import DataFrame

##read in the bag and p value file in tab delimited format
bagSNP_file = pd.read_csv("bagSNP.chr22.print.txt.gz", sep='\t', header=None)
bagPval_file = pd.read_csv("bagPvalue.chr22.print.txt.gz", sep='\t', header=None)
bagSNP_file.head()
bagPval_file.head()

## create a panda dataframe table 
## this table have each row of bagSNPs and index as thier bag numbers
bagnum = sum(1 for row in bagSNP_file)
bag_df = pd.DataFrame(bagSNP_file, index = bagnum)

# Try to finish the last 3 steps by Friday Taiwan Time.
##2.1 import GWAS file
GWAS_file = pd.read_csv("LSBMD_CHRPOS_PVALUE_done.gz", sep='\t', header=None)
GWAS_df = pd.DataFrame(bagSNP_file)
GWAS_df.head()
##2.2 use each GWAS SNP to search through the SNPs in the bags.
## If it matches the position, then assign the SNP to the new GWAS bag table
## else drop it to a list to collect the unmapped GWAS SNPs
for i in range(bagnum):
  bag_df.loc[[i], :]




bagSNP_file.close()
bagPval_file.close()
GWAS_file.close()
