
# Try to finish the first 3 steps by Thursday Taiwan Time.
##1.1  import eQTL bags in python panda dataframe
import pandas as pd
bagSNP_file = open("bagSNP.chr22.print.txt.gz", "r")
bagPval_file = open("bagPvalue.chr22.print.txt.gz", "r")



##1.2  set the first column as each SNPs positions

##1.3  set the second column as the bag index




# Try to finish the last 3 steps by Friday Taiwan Time.
##2.1 import GWAS file
GWAS_file = open("LSBMD_CHRPOS_PVALUE_done.gz", "r")
##2.2 map each GWAS SNPs with the first column of the dataframe

##2.3 if mapped, put into each bag list, else drop it

bagSNP_file.close()
bagPval_file.close()
GWAS_file.close()
