setwd("C://Users//User//Desktop//eQTL_enrichment//")

list.files() 

testGWAS<-read.table(gzfile("C:\\Users\\User\\Desktop\\eQTL_project\\openGWAS\\T2D_RELATED_TRAITS.merged.data.chr.pos.sorted.2014.06.11.BP1.txt.gz"))
head(testGWAS)

testGWASnew <- subset(testGWAS, select=c(4, 5, 2, 3, 1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30))
head(testGWASnew)

testGWASnew$V4 <- sub("^", "chr", testGWASnew$V4 )
head(testGWASnew)



testGWASnew["V31"] <- NA

head(testGWASnew)

testGWASnew[1,31] = "startpos"

head(testGWASnew)

testGWASnew$V31 <- (as.numeric(testGWASnew$V5)-1)

head(testGWASnew)

testGWASnew <- subset(testGWAS, select=c(4, 31, 5, 2, 3, 1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30))

write.csv(file="C://Users//User//Desktop//eQTL_enrichment//data//BEDformat_test_T2D_merged")



eQTL <- read.table("C://Users//User//Desktop//eQTL_enrichment//data//all.cis-eQTL.chrPos.smallestPvalue.txt")

colnames(eQTL)=c("chrPos","pvalue")
head(eQTL)
dim(eQTL)
eQTL <- subset(eQTL, pvalue<0.05)

directory="./data"
(GWASFileList <- grep("lessThan10-5.chrPosition.pvalue", list.files(directory),value=T))

hyperPvalue<-NULL
for(i in 1:length(GWASFileList)){
cat(i,"\n")
#GWAS <- read.table("C://Users//User//Desktop//eQTL_enrichment//data//P_TG_lessThan10-5.chrPosition.pvalue")
GWAS <- read.table(paste("C://Users//User//Desktop//eQTL_enrichment//data//",GWASFileList[i],sep=""))

colnames(GWAS)=c("chrPos","pvalue")
head(GWAS)
dim(GWAS)
GWAS <- subset(GWAS, pvalue< 10e-8 && pvalue > 10e-10 )
GWAS$pvalue
##phyer (sig.eQTL & sig. GWAS, sig. eQTL, non-sig eQTL, sig. GWAS)
#x, q	
#vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.
#(GWAS SNPs with p-value < 10-5 and exclude SNPs not in eQTL results; with significant eQTL p-value)
sig_eQTL <- subset(eQTL, eQTL$pvalue < 10e-8)
head(sig_eQTL)
dim(sig_eQTL)
intersect_sigeQTL_GWAS  <- intersect(sig_eQTL$chrPos, GWAS$chrPos)

q <- length(intersect_sigeQTL_GWAS)
q
x=q
#m	
#the number of white balls in the urn.
#(eQTL SNPs with p-value <= 10-8 )
m <- length(which(eQTL$pvalue <= 10e-8))
m
#n	
#the number of black balls in the urn.
#(eQTL SNPs with p-value > 10-8 )
n <- length(which(eQTL$pvalue > 10e-8))
n
m+n
dim(eQTL)
#k	
#the number of balls drawn from the urn.
#(GWAS SNPs with p-value < 10-5 and exclude SNPs not in eQTL results)
inter_GWAS <- intersect(eQTL$chrPos, GWAS$chrPos)
inter_GWAS
k <- length(inter_GWAS)
k

hyperPvalue[i]<-1-phyper(x, m, n, k, lower.tail = T)
c(x,m,n,k)
}


names(hyperPvalue)=sub("_lessThan10-5.chrPosition.pvalue", "", GWASFileList)
sort(hyperPvalue)


write.csv(sort(hyperPvalue), file="C://Users//User//Desktop//eQTL_enrichment//data//result_GWAS_10-8-10.csv")



x=0:k

sum(dhyper(x, m, n, k))

plot(c(0:k),dhyper(x, m, n, k), xlab="number of  sig. GWAS & sig. eQTL", ylab="probability(frequency)", xlim=c(0,100))
lines(c(0:k),dhyper(x, m, n, k),col="red")
abline(h=0)
abline(v=2, col="green")

1-phyper(x, m, n, k, lower.tail = T)
phyper(x, m, n, k, lower.tail = F)
c(x,m,n,k)
