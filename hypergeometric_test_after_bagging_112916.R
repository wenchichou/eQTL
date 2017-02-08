#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
## GWAS preprocessing in cluster, then set GWAS P value cut off in R program. After it, bagging GWAS
## At the end, conduct hypergeometric test
## GWAS data should have chr:pos and p value, bagging results should be revised to chr:pos

#========================================
# Configure these parameters
#========================================
dirBagging <- "/home/unix/wcchou/gsapWenChi/gautvik/results/bagging/"
GWASfilePath <-"/home/unix/wcchou/gsapWenChi/gautvik/bin/eQTL/data/LSBMD_CHRPOS_PVALUE_done.gz"
sigGWASpvalueCutoff <- 1e-4 # 10e-8 is 0.0000001; 1e-8 is 0.00000001
EQTLpvalueRangeSmall <- 1e-4
EQTLpvalueRangeLarge <- 1e-1000

dirBagging <- args[1]
GWASfilePath <- args[2]  
sigGWASpvalueCutoff <- args[3]
sigGWASpvalueCutoff <- as.numeric(as.character(sigGWASpvalueCutoff))
EQTLpvalueRangeSmall <- args[4]
EQTLpvalueRangeSmall <- as.numeric(as.character(EQTLpvalueRangeSmall))
EQTLpvalueRangeLarge <- args[5]
EQTLpvalueRangeLarge <- as.numeric(as.character(EQTLpvalueRangeLarge))
outputFilePath <- args[6]
#========================================
## 1. Load in bagging files, eQTL and GWAS dataset

#1. Read all GWAS SNPs
cat("Reading all GWAS SNPs ...\n")
GWAS <- read.table(GWASfilePath, header=T)
colnames(GWAS) <- c("CHR.POS","pvalue")
head(GWAS)
dim(GWAS)

#2. Read all eQTL bags including all SNPs
cat("Reading bagged eQTL SNPs ...\n")
all.eQTL.bagSNP <- list()
all.eQTL.bagPvalue <- NULL
for(CHR in seq(1,22,1)){
	cat("loading chr",CHR,"\n")
	bagSNP.CHR.path <- paste(dirBagging,"bagSNP.chr",CHR,".print.txt",sep="")
	bagPvalue.CHR.path <- paste(dirBagging,"bagPvalue.chr",CHR,".print.txt",sep="")
	x <- scan(bagSNP.CHR.path, what="", sep="\n")
	y <- strsplit(x, "[[:space:]]+")
	bagSNP.CHR <- lapply(y, function(list) paste(CHR,":",list,sep=""))
	bagPvalue.CHR <- read.table(bagPvalue.CHR.path)
	all.eQTL.bagSNP <- c(all.eQTL.bagSNP, bagSNP.CHR)
	all.eQTL.bagPvalue <- c(all.eQTL.bagPvalue, unlist(bagPvalue.CHR))
}
rm(x, y, bagSNP.CHR, bagPvalue.CHR)
length(all.eQTL.bagSNP)
length(all.eQTL.bagPvalue)

#3. put GWAS SNPs into LD bags

# find GWAS SNPs overlapped to eQTL SNPs
# for loop all list of all.eQTL.bagSNP to find index.intersect_eQTL_GWAS 

#index.intersect_eQTL_GWAS <- match(GWAS$CHR.POS, unlist(all.eQTL.bagSNP), nomatch=0)
#length(which(index.intersect_eQTL_GWAS!=0))
#overlappedGWAS.SNPs <- GWAS[index.intersect_eQTL_GWAS,]


	res <- lapply(head(all.eQTL.bagSNP), function(list) match(GWAS[i,1], list))


for(i in 1:nrow(GWAS)){
	res <- lapply(all.eQTL.bagSNP, function(ch) match(GWAS$CHR.POS[i], ch))
	# which vectors contain a search term
	sapply(res, function(x) length(x) > 0)
}

all.eQTL.bagSNP
all.eQTL.bagPvalue



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


#3. output1 GWASbag.SNP matrix
#3. output2. GWASbag.smallestPvalue matrix

#4. get significant and non-significant GWAS bags
sig.GWAS.bag.all <- list()
sig.GWAS.bag.Pvalue.all <- NULL
nonSig.GWAS.bag.all <- list()
nonSig.GWAS.bag.Pvalue.all <- NULL

sig.index <- which(GWASbag.smallestPvalue <= sigGWASpvalueCutoff) 
sig.GWAS.bag.all <- GWASbag.SNP[sig.index]
sig.GWAS.bag.Pvalue.all <- GWASbag.smallestPvalue[sig.index]	

nonSig.index <- which(GWASbag.smallestPvalue > sigGWASpvalueCutoff) 
nonSig.GWAS.bag.all <- GWASbag.SNP[nonSig.index]
nonSig.GWAS.bag.Pvalue.all <- GWASbag.smallestPvalue[nonSig.index]	

#5. get eQTL bags according to p-value ranges 

#5. output subsetEQTL.bag
#5. output subsetEQTL.bag.Pvalue



















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
#str(GWAS)

## use p-value ranges to select subset of GWAS SNPs
subsetGWAS <-subset(GWAS[complete.cases(GWAS),], pvalue <= sigGWASpvalueCutoff) #windows system may need && rather than &
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
#wcc#cat("Analyzed GWAS: ",gsub("^.*/","",GWASfilePath, perl=TRUE),"\n" )
#wcc#cat("bagged eQTL significant p-value cutoff is ",sigGWASpvalueCutoff,"\n")
#wcc#cat("How many GWAS SNPs with p-values? ",nrow(GWAS[complete.cases(GWAS),]),"\n")
#wcc#cat("GWAS p-value ranges used to get subset GWAS are from ",EQTLpvalueRangeSmall,"to ",EQTLpvalueRangeLarge,"\n")
#wcc#cat("How many subset GWAS SNPs? ",nrow(subsetGWAS),"\n")
#wcc#cat("How many subset GWAS SNPs not covered by eQTL bags? ",numSubsetGWASnotINeQTL,"\n")
#wcc#cat("x =",x,"\n")
#wcc#cat("m =",m,"\n")
#wcc#cat("n =",n,"\n")
#wcc#cat("k =",k,"\n")
#wcc#cat("enrichment Pvalue = ",enrichmentPvalue,"\n")

#cat(gsub("^.*/","",GWASfilePath, perl=TRUE),"\t",sigGWASpvalueCutoff,"\t",nrow(GWAS[complete.cases(GWAS),]),"\t",EQTLpvalueRangeSmall,"\t",EQTLpvalueRangeLarge,"\t",nrow(subsetGWAS),"\t",numSubsetGWASnotINeQTL,"\t",x,"\t",m,"\t",n,"\t",k,"\t",enrichmentPvalue,"\n",sep="")

sink(outputFilePath)
cat(gsub("^.*/","",GWASfilePath, perl=TRUE),"\t",sigGWASpvalueCutoff,"\t",nrow(GWAS[complete.cases(GWAS),]),"\t",EQTLpvalueRangeSmall,"\t",EQTLpvalueRangeLarge,"\t",nrow(subsetGWAS),"\t",numSubsetGWASnotINeQTL,"\t",x,"\t",m,"\t",n,"\t",k,"\t",enrichmentPvalue,"\n",sep="")
sink()



