#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## This file perform enrichment analysis of significant GWAS with eQTL

#========================================
# Configure these parameters
#========================================
setwd("/restricted//projectnb//chengcancer//Simon//project//T2D//bulk_RNA-seq//merge_170415_NB501164_170506_NB501164//2mapping//QT//")
## take the whole eQTL bags and chr1 of HEIGHT GWAS bags as an example  
#eQTLfilePath <- "/restricted//projectnb//chengcancer//Simon//project//T2D//bulk_RNA-seq//merge_170415_NB501164_170506_NB501164//2mapping//QT//outdir//"
eQTLfilePath <- args[1]
#eQTLfilePath <- "/home/unix/wcchou/gsapWenChi/gautvik/results/bagging/" 
#GWASfilePath <- "/restricted//projectnb//chengcancer//Simon//project//T2D//bulk_RNA-seq//merge_170415_NB501164_170506_NB501164//2mapping//QT//eQTL/"
GWASfilePath <- args[2]

## set eQTL in different bins

#eQTLpvalueUpperBound <- 1e-2
eQTLpvalueUpperBound <- args[3]
eQTLpvalueUpperBound <- as.numeric(as.character(eQTLpvalueUpperBound))
#eQTLpvalueLowerBound <- 1e-100
eQTLpvalueLowerBound <- args[4]
eQTLpvalueLowerBound <- as.numeric(as.character(eQTLpvalueLowerBound))
sigGWASpvalueCutoff <- args[5]
sigGWASpvalueCutoff <- as.numeric(as.character(sigGWASpvalueCutoff))
outputFilePath <- args[6]
#========================================


## 1. read in whole eQTL bags data (append from chr 1 to chr 22)

# get significant eQTL bags
# get non-significant eQTL bags
cat("Reading bagged eQTL SNPs ...\n")
sig.eQTL_bagSNP.all <- list()
sig.eQTL_bagPvalue.all <- NULL

for(CHR in seq(1,22,1)){
  cat("Reading CHR",CHR,"...\n")
  eQTL_bagSNP.path <- paste(eQTLfilePath,"//GTeXbag.chr",CHR,"_bagSNP",sep="")
  eQTL_bagPvalue.path <- paste(eQTLfilePath,"//GTeXbag.chr",CHR,"_bagPvalue",sep="")
  x <- scan(eQTL_bagSNP.path, what="", sep="\n")
  y <- strsplit(x, "[[:space:]]+")
  for (bag in 1:length(y)){
    for (element in 1:length(y[[bag]])){
      y[[bag]][element] <- paste(CHR, ":", y[[bag]][element], sep='')
    }
  }
  eQTL_bagSNP <- y
  # exclude non-significant eQTL bags
  eQTL_bagPvalue <- read.table(eQTL_bagPvalue.path)
  sig.index <- which(eQTL_bagPvalue[,1] < eQTLpvalueUpperBound & eQTL_bagPvalue[,1] >= eQTLpvalueLowerBound) 
  sig.eQTL_bagSNP <- eQTL_bagSNP[sig.index]
  sig.eQTL_bagPvalue <- eQTL_bagPvalue[sig.index,]
  sig.eQTL_bagSNP.all <- c(sig.eQTL_bagSNP.all, sig.eQTL_bagSNP)
  sig.eQTL_bagPvalue.all <- c(sig.eQTL_bagPvalue.all, sig.eQTL_bagPvalue)
  
}
length(sig.eQTL_bagSNP.all)
length(sig.eQTL_bagPvalue.all)

# wcc020917LoopAllGWAS #for loop all GWSA
# wcc020917LoopAllGWAS for(GWASDIR in paste(grep("_Pvalue$",list.dirs(path = "/home/unix/wcchou/gsapWenChi/gautvik/results/bagging.gwas/", full.names = TRUE), value=T),"/",sep="")){
# wcc020917LoopAllGWAS 	cat("running ",GWASDIR,"\n")
# wcc020917LoopAllGWAS GWASfilePath <- GWASDIR

## 2. read in specific type of GWAS bags data (append from chr 1 to chr 22)
## HEIGHT GWAS as an example

# get significant GWAS bags
# get non-significant GWAS bags
cat("Reading bagged GWAS SNPs ...\n")
sig.GWAS_bagSNP.all <- list()
sig.GWAS_bagPvalue.all <- NULL
nonSig.GWAS_bagSNP.all <- list()
nonSig.GWAS_bagPvalue.all <- NULL
bagnum = 0
for(CHR in seq(1,22,1)){
  cat("Reading CHR",CHR,"...\n")
  GWAS_bagSNP.path <- paste(GWASfilePath,"//GTeXbag.chr",CHR,"_bagSNP",sep="")
  GWAS_bagPvalue.path <- paste(GWASfilePath,"//GTeXbag.chr",CHR,"_bagPvalue",sep="")
  x <- scan(GWAS_bagSNP.path, what="", sep="\n")
  y <- strsplit(x, "[[:space:]]+")
  for (bag in 1:length(y)){
    for (element in 1:length(y[[bag]])){
      y[[bag]][element] <- paste(CHR, ":", y[[bag]][element], sep='')
    }
  }
  GWAS_bagSNP <- y
  bagnum = bagnum + length(y)
  # exclude non-significant GWAS bags
  GWAS_bagPvalue <- read.table(GWAS_bagPvalue.path)
  sig.index <- which(GWAS_bagPvalue[,1] < sigGWASpvalueCutoff) 
  sig.GWAS_bagSNP <- GWAS_bagSNP[sig.index]
  sig.GWAS_bagPvalue <- GWAS_bagPvalue[sig.index,]
  sig.GWAS_bagSNP.all  <- c(sig.GWAS_bagSNP.all, sig.GWAS_bagSNP)
  sig.GWAS_bagPvalue.all <- c(sig.GWAS_bagPvalue.all, sig.GWAS_bagPvalue)
  nonSig.GWAS_index <- which(GWAS_bagPvalue[,1] >= sigGWASpvalueCutoff) 
  nonSig.GWAS_bagSNP <- GWAS_bagSNP[nonSig.GWAS_index]
  nonSig.GWAS_bagPvalue <- GWAS_bagPvalue[nonSig.GWAS_index,]
  nonSig.GWAS_bagSNP.all <- c(nonSig.GWAS_bagSNP.all, nonSig.GWAS_bagSNP)
  nonSig.GWAS_bagPvalue.all <- c(nonSig.GWAS_bagPvalue.all, nonSig.GWAS_bagPvalue)
}
length(sig.GWAS_bagSNP.all)
length(sig.GWAS_bagPvalue.all)
length(nonSig.GWAS_bagSNP.all)
length(nonSig.GWAS_bagPvalue.all)


# list number of SNPs in each significant GWAS bag
length.sig.GWAS_bagSNP.all <- unlist(do.call(rbind, lapply(sig.GWAS_bagSNP.all, function(x) length(x)))[,1])
# generate bag index for each SNP
sig.GWAS_bagSNP.all.expendedBag <- rep(c(1:length(sig.GWAS_bagSNP.all)), length.sig.GWAS_bagSNP.all)

# list number of SNPs in each significant eQTL bag
length.sig.eQTL_bagSNP.all <- unlist(do.call(rbind, lapply(sig.eQTL_bagSNP.all, function(x) length(x)))[,1])
# generate bag index for each SNP
sig.eQTL_bagSNP.all.expendedBag <- rep(c(1:length(sig.eQTL_bagSNP.all)), length.sig.eQTL_bagSNP.all)

# list number of SNPs in each non-significant GWAS bag
length.nonSig.GWAS_bagSNP.all <- unlist(do.call(rbind, lapply(nonSig.GWAS_bagSNP.all, function(x) length(x)))[,1])
# generate bag index for each SNP
nonSig.GWAS_bagSNP.all.expendedBag <- rep(c(1:length(nonSig.GWAS_bagSNP.all)), length.nonSig.GWAS_bagSNP.all)


sigGWAS <- unlist(sig.GWAS_bagSNP.all)


## make sure the order of eQTL and GWAS
length(sig.eQTL_bagSNP.all)
length(sig.eQTL_bagPvalue.all)
length(sig.GWAS_bagSNP.all)
length(sig.GWAS_bagPvalue.all)
length(nonSig.GWAS_bagSNP.all)
length(nonSig.GWAS_bagPvalue.all)

# get all first elements of the list of significant GWAS bag 
firstElement_sigGWASBag <- unlist(lapply(sig.GWAS_bagSNP.all, `[[`, 1))
plot.sigGWASbag2sigeQTLbag.index <- sig.eQTL_bagSNP.all.expendedBag[match(firstElement_sigGWASBag, unlist(sig.eQTL_bagSNP.all))]
plot.sigGWASbag2sigeQTLbag.pvalue <- plot.sigGWASbag2sigeQTLbag.index
plot.sigGWASbag2sigeQTLbag.pvalue[which(!is.na(plot.sigGWASbag2sigeQTLbag.index))] <- sig.eQTL_bagPvalue.all[(plot.sigGWASbag2sigeQTLbag.index[!is.na((plot.sigGWASbag2sigeQTLbag.index))])]

# get all first elements of the list of non-significant GWAS bag 
firstElement_nonSigGWASBag <- unlist(lapply(nonSig.GWAS_bagSNP.all, `[[`, 1))
plot.nonSigGWASbag2sigeQTLbag.index <- sig.eQTL_bagSNP.all.expendedBag[match(firstElement_nonSigGWASBag, unlist(sig.eQTL_bagSNP.all))]
plot.nonSigGWASbag2sigeQTLbag.pvalue <- plot.nonSigGWASbag2sigeQTLbag.index
plot.nonSigGWASbag2sigeQTLbag.pvalue[which(!is.na(plot.nonSigGWASbag2sigeQTLbag.index))] <- sig.eQTL_bagPvalue.all[(plot.nonSigGWASbag2sigeQTLbag.index[!is.na((plot.nonSigGWASbag2sigeQTLbag.index))])]

# merge two set of pvalues of eQTL bags 
plot.eQTL.GWAS.pvalues <- data.frame(eQTLbagPvalue=c(plot.sigGWASbag2sigeQTLbag.pvalue, plot.nonSigGWASbag2sigeQTLbag.pvalue), GWASbagPvalue=c(sig.GWAS_bagPvalue.all, nonSig.GWAS_bagPvalue.all))
plot.eQTL.GWAS.pvalues <- plot.eQTL.GWAS.pvalues[complete.cases(plot.eQTL.GWAS.pvalues),]
dim(plot.eQTL.GWAS.pvalues)
head(plot.eQTL.GWAS.pvalues)
# make a plot
pdf(file=paste(basename(eQTLfilePath), "_" , basename(GWASfilePath), "_", eQTLpvalueUpperBound, "_", eQTLpvalueLowerBound,  "scatter.pdf", sep=""))

plot(-log10(plot.eQTL.GWAS.pvalues$eQTLbagPvalue), -log10(plot.eQTL.GWAS.pvalues$GWASbagPvalue), asp=1, col=rgb(0,0,0,0.1), pch=19, xlim=c(0,20), ylim=c(0,20), xlab=paste("-log10", " ", basename(eQTLfilePath), "pvalue"), ylab=paste("-log10", " ", basename(GWASfilePath), "pvalue"), main=paste(basename(eQTLfilePath),"; ",nrow(plot.eQTL.GWAS.pvalues)," overlapped bags\n",length(sig.eQTL_bagPvalue.all)," and ",length(sig.GWAS_bagPvalue.all)+length(nonSig.GWAS_bagPvalue.all)," bags in input eQTL and GWAS",sep=""))
abline(h = 9, col = "blue")
abline(h = 5, col = "blue")
abline(h = 7, col = "blue")
abline(v = 9, col = "red")
abline(v = 5, col = "red")
abline(v = 7, col = "red")
dev.off()


## make sig eQTL list for dataframe
eQTL_Pval_list.intersect_sigeQTL_GWAS <- sig.eQTL_bagPvalue.all[unique(sig.eQTL_bagSNP.all.expendedBag[match(sigGWAS, unlist(sig.eQTL_bagSNP.all), nomatch=0)])]

## get number of bags belonging to sigGWAS and sig eQTL
index.intersect_sigeQTL_sigGWAS <- match(unlist(sig.eQTL_bagSNP.all), sigGWAS, nomatch=0)

length(which(index.intersect_sigeQTL_sigGWAS!=0))
bag.intersect_sigeQTL_sigGWAS <- sig.GWAS_bagSNP.all.expendedBag[index.intersect_sigeQTL_sigGWAS]
GWAS_Pval_list.intersect_sigeQTL_sigGWAS<- sig.GWAS_bagPvalue.all[unique(bag.intersect_sigeQTL_sigGWAS)]

length(bag.intersect_sigeQTL_sigGWAS)
(numberOfbag.intersect_sigeQTL_sigGWAS <- length(table(bag.intersect_sigeQTL_sigGWAS)))

## get number of bags belonging to non-sig.GWAS and sig.eQTL
index.intersect_nonSigGWAS_SigeQTL <- match(unlist(sig.eQTL_bagSNP.all), unlist(nonSig.GWAS_bagSNP.all),nomatch=0)
length(which(index.intersect_nonSigGWAS_SigeQTL!=0))
bag.intersect_nonSigGWAS_SigeQTL <- nonSig.GWAS_bagSNP.all.expendedBag[index.intersect_nonSigGWAS_SigeQTL]

GWAS_Pval_list.intersect_sigeQTL_nonsigGWAS<- nonSig.GWAS_bagPvalue.all[unique(bag.intersect_nonSigGWAS_SigeQTL)]

length(bag.intersect_nonSigGWAS_SigeQTL)
(numberOfbag.intersect_nonSigGWAS_SigeQTL <- length(table(bag.intersect_nonSigGWAS_SigeQTL)))

##3. Compare the size of GWAS and eQTL and conduct hypergeometric test
## P Value:  the probability of obtaining a result equal to or "more extreme" 
##          than what was actually observed, when the null hypothesis is true. 

(q <- numberOfbag.intersect_sigeQTL_sigGWAS)
x=q

#m  
#the number of white balls in the urn.
#(GWAS bags with p-value <= sigGWASpvalueCutoff )
(m <- length(sig.GWAS_bagSNP.all))
#n  
#the number of black balls in the urn.
#(GWAS bags with p-value > sigGWASpvalueCutoff )
(n <- length(nonSig.GWAS_bagSNP.all))
m+n

#k  
#the number of balls drawn from the urn.
# sig. eQTL bags overlapped with significant and non-significant GWAS bags
(k <- numberOfbag.intersect_sigeQTL_sigGWAS + numberOfbag.intersect_nonSigGWAS_SigeQTL)

## Print out results

## The meaning of lower.tail and 1 - phyper?

enrichmentPvalue <- (1-phyper(x, m, n, k, lower.tail = T))
cat("Analyzed eQTL: ", basename(eQTLfilePath), "\n")
cat("Analyzed GWAS: ", basename(GWASfilePath),"\n" )
cat("sig. GWAS bags p-value cutoff is ",sigGWASpvalueCutoff,"\n")
#wcc#cat("How many GWAS SNPs with p-values? ",nrow(GWAS[complete.cases(GWAS),]),"\n")
cat("eQTL bags P value range from ",eQTLpvalueLowerBound,"to ",eQTLpvalueUpperBound,"\n")
cat("How many sig. eQTL bags? ",length(sig.eQTL_bagSNP.all),"\n")
cat("How many sig. GWAS bags? ",length(sig.GWAS_bagSNP.all),"\n")
cat("How many non-sig. GWAS bags? ",length(nonSig.GWAS_bagSNP.all),"\n")

cat("x =",x,"\n")
cat("m =",m,"\n")
cat("n =",n,"\n")
cat("k =",k,"\n")
cat("enrichment Pvalue = ",enrichmentPvalue,"\n")




#cat(gsub("^.*/","",GWASfilePath, perl=TRUE),"\t",sigGWASpvalueCutoff,"\t",nrow(GWAS[complete.cases(GWAS),]),"\t",EQTLpvalueRangeSmall,"\t",EQTLpvalueRangeLarge,"\t",nrow(subsetGWAS),"\t",numSubsetGWASnotINeQTL,"\t",x,"\t",m,"\t",n,"\t",k,"\t",enrichmentPvalue,"\n",sep="")

#sink(outputFilePath)
#cat(gsub("^.*/","",GWASfilePath, perl=TRUE),"\t",sigGWASpvalueCutoff,"\t",nrow(GWAS[complete.cases(GWAS),]),"\t",EQTLpvalueRangeSmall,"\t",EQTLpvalueRangeLarge,"\t",nrow(subsetGWAS),"\t",numSubsetGWASnotINeQTL,"\t",x,"\t",m,"\t",n,"\t",k,"\t",enrichmentPvalue,"\n",sep="")
#sink()


## produce scatter plot of matched SNP's GWAS and eQTL P value distribution
## sig eQTL overlapped with both sig and non-sig GWAS



