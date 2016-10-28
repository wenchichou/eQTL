## GWAS preprocessing in cluster, then set GWAS P value cut off in R program. After it, bagging GWAS
## At the end, conduct hypergeometric test
## GWAS data should have chr:pos and p value, bagging results should be revised to chr:pos



setwd("C:\\Users\\User\\Desktop\\eQTL_project")

list.files() 

## 1. Load in bagging files, eQTL and GWAS dataset

chr1_Pvalue <- read.table("finished.Bagging.NoMissing\\bagPvalue.chr1.print.txt")
head(chr1_Pvalue)
dim(chr1_Pvalue)

# Read in the data
x <- scan("finished.Bagging.NoMissing\\bagSNP.chr1.print.txt", what="", sep="\n")
# Separate elements by one or more whitepace
y <- strsplit(x, "[[:space:]]+")
# Extract the first vector element and set it as the list element name
names(y) <- sapply(y, `[[`, 1)
chr1_bag_SNP <- y
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
#y <- lapply(y, `[`, -1)
#y <- lapply(y, function(x) x[-1]) # same as above


GWAS <- read.table("C:\\Users\\User\\Desktop\\eQTL_project\\openGWAS\\T2D_RELATED_TRAITS.merged.data.chr.pos.sorted.2014.06.11.BP1.txt.gz")
head(GWAS)
dim(GWAS)

GWAS_LSBMD <- read.table("C:\\Users\\User\\Desktop\\eQTL_project\\openGWAS\\LSBMD.ALLCHR.CHRPOS.gz")
head(GWAS_LSBMD)
dim(GWAS_LSBMD)


##2. Assign GWAS and eQTL depending on the bagging results





##3. Compare the size of GWAS and eQTL and conduct hypergeometric test






##4. Print out results
