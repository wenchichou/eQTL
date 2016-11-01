## GWAS preprocessing in cluster, then set GWAS P value cut off in R program. After it, bagging GWAS
## At the end, conduct hypergeometric test
## GWAS data should have chr:pos and p value, bagging results should be revised to chr:pos


#========================================
# Configure these parameters
#========================================
dirBagging <- "/home/unix/wcchou/gsapWenChi/gautvik/results/bagging/"
#GWASfilePath <- "/home/unix/wcchou/gsapWenChi/gautvik/data/gwas.2.5m/LSBMD_CHRPOS_PVALUE_done.gz"
GWASfilePath <- "/home/unix/wcchou/gsapWenChi/gautvik/data/gwas.2.5m/ALZ_CHRPOS_PVALUE_done.gz"
sigEQTLpvalueCutoff <- 1e-8 # 10e-8 is 0.0000001; 1e-8 is 0.00000001
GWASpvalueRangeSmall <- 1e-6
GWASpvalueRangeLarge <- 1e-1000
#========================================
## 1. Load in bagging files, eQTL and GWAS dataset

# get significant eQTL bags
# get non-significant eQTL bags
cat("Reading bagged eQTL SNPs ...\n")
sig.bagSNP.all <- list()
sig.bagPvalue.all <- NULL
nonSig.bagSNP.all <- list()
nonSig.bagPvalue.all <- NULL
for(CHR in seq(1,22,1)){
	bagSNP.CHR.path <- paste(dirBagging,"bagSNP.chr",CHR,".print.txt2",sep="")
	bagPvalue.CHR.path <- paste(dirBagging,"bagPvalue.chr",CHR,".print.txt",sep="")
	x <- scan(bagSNP.CHR.path, what="", sep="\n")
	y <- strsplit(x, "[[:space:]]+")
	bagSNP.CHR <- y
	# exclude non-significant eQTL bags
	bagPvalue.CHR <- read.table(bagPvalue.CHR.path)
	sig.index <- which(bagPvalue.CHR[,1] < sigEQTLpvalueCutoff) 
	sig.bagSNP.CHR <- bagSNP.CHR[sig.index]
	sig.bagPvalue.CHR <- bagPvalue.CHR[sig.index,]
	sig.bagSNP.all <- c(sig.bagSNP.all, sig.bagSNP.CHR)
	sig.bagPvalue.all <- c(sig.bagPvalue.all, sig.bagPvalue.CHR)
	nonSig.index <- which(bagPvalue.CHR[,1] >= sigEQTLpvalueCutoff) 
	nonSig.bagSNP.CHR <- bagSNP.CHR[nonSig.index]
	nonSig.bagPvalue.CHR <- bagPvalue.CHR[nonSig.index,]
	nonSig.bagSNP.all <- c(nonSig.bagSNP.all, nonSig.bagSNP.CHR)
	nonSig.bagPvalue.all <- c(nonSig.bagPvalue.all, nonSig.bagPvalue.CHR)
}
length(sig.bagSNP.all)
length(sig.bagPvalue.all)
length(nonSig.bagSNP.all)
length(nonSig.bagPvalue.all)

# list all bag number of each significant eQTL SNPs in all bags
length.sig.bagSNP.all <- unlist(do.call(rbind, lapply(sig.bagSNP.all, function(x) length(x)))[,1])
sig.bagSNP.all.expendedBag <- rep(c(1:length(sig.bagSNP.all)), length.sig.bagSNP.all)

# list all bag number of each non-significant eQTL SNPs in all bags
length.nonSig.bagSNP.all <- unlist(do.call(rbind, lapply(nonSig.bagSNP.all, function(x) length(x)))[,1])
nonSig.bagSNP.all.expendedBag <- rep(c(1:length(nonSig.bagSNP.all)), length.nonSig.bagSNP.all)



# get GWAS data ready
cat("Reading GWAS SNPs ...\n")
GWAS <- read.table(GWASfilePath, header=T)
colnames(GWAS) <- c("CHR.POS","pvalue")
head(GWAS)
dim(GWAS)
str(GWAS)

## use p-value ranges to select subset of GWAS SNPs
subsetGWAS <-subset(GWAS[complete.cases(GWAS),], pvalue < GWASpvalueRangeSmall & pvalue >= GWASpvalueRangeLarge) #windows system may need && rather than &
head(subsetGWAS)
dim(subsetGWAS)

# how many subset GWAS SNPs are not covered by eQTL SNPs
(numSubsetGWASnotINeQTL <- length(which(is.na(match(subsetGWAS[,1], c(unlist(sig.bagSNP.all), unlist(nonSig.bagSNP.all)))))))

## get number of bags belonging to subsetGWAS and significant eQTL
index.intersect_sigeQTL_subsetGWAS <- match(subsetGWAS$CHR.POS, unlist(sig.bagSNP.all), nomatch=0)
length(which(index.intersect_sigeQTL_subsetGWAS!=0))
bag.intersect_sigeQTL_subsetGWAS <- sig.bagSNP.all.expendedBag[index.intersect_sigeQTL_subsetGWAS]
length(bag.intersect_sigeQTL_subsetGWAS)
(numberOfbag.intersect_sigeQTL_subsetGWAS <- length(table(bag.intersect_sigeQTL_subsetGWAS)))

## get number of bags belonging to subsetGWAS and non-significant eQTL
index.intersect_nonSigeQTL_subsetGWAS <- match(subsetGWAS$CHR.POS, unlist(nonSig.bagSNP.all), nomatch=0)
length(which(index.intersect_nonSigeQTL_subsetGWAS!=0))
bag.intersect_nonSigeQTL_subsetGWAS <- nonSig.bagSNP.all.expendedBag[index.intersect_nonSigeQTL_subsetGWAS]
length(bag.intersect_nonSigeQTL_subsetGWAS)
(numberOfbag.intersect_nonSigeQTL_subsetGWAS <- length(table(bag.intersect_nonSigeQTL_subsetGWAS)))

##3. Compare the size of GWAS and eQTL and conduct hypergeometric test
(q <- numberOfbag.intersect_sigeQTL_subsetGWAS)
x=q
#m  
#the number of white balls in the urn.
#(eQTL SNPs with p-value <= 10-8 )
(m <- length(sig.bagSNP.all))
#n  
#the number of black balls in the urn.
#(eQTL SNPs with p-value > 10-8 )
(n <- length(nonSig.bagSNP.all))
m+n
#k  
#the number of balls drawn from the urn.
#(GWAS SNPs with p-value < 10-5 and exclude SNPs not in eQTL results)
(k <- numberOfbag.intersect_sigeQTL_subsetGWAS + numberOfbag.intersect_nonSigeQTL_subsetGWAS) 

## Print out results

enrichmentPvalue <- (1-phyper(x, m, n, k, lower.tail = T))
cat("Analyzed GWAS: ",gsub("^.*/","",GWASfilePath, perl=TRUE),"\n" )
cat("bagged eQTL significant p-value cutoff is ",sigEQTLpvalueCutoff,"\n")
cat("How many GWAS SNPs with p-values? ",nrow(GWAS[complete.cases(GWAS),]),"\n")
cat("GWAS p-value ranges used to get subset GWAS are from ",GWASpvalueRangeSmall,"to ",GWASpvalueRangeLarge,"\n")
cat("How many subset GWAS SNPs? ",nrow(subsetGWAS),"\n")
cat("How many subset GWAS SNPs not covered by eQTL bags? ",numSubsetGWASnotINeQTL,"\n")
cat("x =",x,"\n")
cat("m =",m,"\n")
cat("n =",n,"\n")
cat("k =",k,"\n")
cat("enrichment Pvalue = ",enrichmentPvalue,"\n")

