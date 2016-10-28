## GWAS preprocessing in cluster, then set GWAS P value cut off in R program. After it, bagging GWAS
## At the end, conduct hypergeometric test
## GWAS data should have chr:pos and p value, bagging results should be revised to chr:pos



setwd("C:\\Users\\User\\Desktop\\eQTL_project")

list.files() 

## 1. Load in bagging files, eQTL and GWAS dataset

#chr1_Pvalue <- read.table("finished.Bagging.NoMissing\\bagPvalue.chr1.print.txt")
#head(chr1_Pvalue)
#dim(chr1_Pvalue)

#WCC1028## Read in the data
#WCC1028#x <- scan("finished.Bagging.NoMissing\\bagSNP.chr1.print.txt", what="", sep="\n")
#WCC1028## Separate elements by one or more whitepace
#WCC1028#y <- strsplit(x, "[[:space:]]+")
#WCC1028## Extract the first vector element and set it as the list element name
#WCC1028#names(y) <- sapply(y, `[[`, 1)
#WCC1028#bag_SNP <- y
#WCC1028##names(y) <- sapply(y, function(x) x[[1]]) # same as above
#WCC1028## Remove the first vector element from each list element
#WCC1028##y <- lapply(y, `[`, -1)
#WCC1028##y <- lapply(y, function(x) x[-1]) # same as above

dirBagging <- "/home/unix/wcchou/gsapWenChi/gautvik/results/bagging/"
sig.bagSNP.all <- list()
sig.bagPvalue.all <- NULL
sigEQTLpvalueCutoff <- 1e-8 # 10e-8 is 0.0000001; 1e-8 is 0.00000001
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
}
length(sig.bagSNP.all)
length(sig.bagPvalue.all)


#GWAS <- read.table("C:\\Users\\User\\Desktop\\eQTL_project\\openGWAS\\T2D_RELATED_TRAITS.merged.data.chr.pos.sorted.2014.06.11.BP1.txt.gz")
#head(GWAS)
#dim(GWAS)

#GWAS_LSBMD <- read.table("C:\\Users\\User\\Desktop\\eQTL_project\\openGWAS\\LSBMD.ALLCHR.CHRPOS.gz")
GWAS <- read.table("/home/unix/wcchou/gsapWenChi/gautvik/data/gwas.2.5m/LSBMD_CHRPOS_PVALUE_done.gz", header=F)
colnames(GWAS) <- c("CHR.POS","pvalue")
head(GWAS)
dim(GWAS)

## 2.0 use p-value range to select subset of GWAS results
subsetGWAS <-subset(GWAS, pvalue < 1e-8 & pvalue > 1e-100) #windows system may need && rather than &
head(subsetGWAS)
dim(subsetGWAS)

##2. Assign GWAS and eQTL depending on the bagging results
# test with all GWAS of ALT
intersect_sigeQTL_GWAS  <- intersect(unlist(sig.bagSNP.all), subsetGWAS$CHR.POS)
intersect_sigeQTL_GWAS  <- intersect(unlist(sig.bagSNP.all), GWAS$CHR.POS)
length(intersect_sigeQTL_GWAS)

overlappedGWAS <- GWAS[match(intersect_sigeQTL_GWAS, GWAS$CHR.POS),]

overlappedGWAS.bag <- NULL
for(i in overlappedGWAS[,1]){
	#res <- lapply(sig.bagSNP.all, function(ch) grep(i, ch))
	res <- lapply(sig.bagSNP.all, function(ch) match(i, ch))
	#bagIndex <- which(sapply(res, function(x) length(x) > 0))
	bagIndex <- which(!is.na(res))
	overlappedGWAS.bag <- c(overlappedGWAS.bag, bagIndex)
}
length(overlappedGWAS.bag)
dim(overlappedGWAS)

overlappedGWAS <- (data.frame(overlappedGWAS, bag=overlappedGWAS.bag))

#find bag count >1 and use the samllest p-value
overlappedGWAS.bagged.multi <- list()
j=1
for(i in names(which(table(overlappedGWAS$bag) > 1))){
	i.oneBag <- overlappedGWAS[overlappedGWAS$bag == i,]
	i.oneBag <- i.oneBag[complete.cases(i.oneBag),]
	if(nrow(i.oneBag)>0){
		overlappedGWAS.bagged.multi[[j]] <- i.oneBag[which.min(i.oneBag$pvalue),];
		j=j+1;
	}
}
overlappedGWAS.bagged.multi <- do.call(rbind.data.frame, overlappedGWAS.bagged.multi)

# put the bag count ==1 back 
overlappedGWAS.bagged.uniq <- list()
j=1
#find bag count >1 and use the samllest p-value
for(i in names(which(table(overlappedGWAS$bag) == 1))){
	i.oneBag <- overlappedGWAS[overlappedGWAS$bag == i,]
	i.oneBag <- i.oneBag[complete.cases(i.oneBag),]
	if(nrow(i.oneBag)>0){
		overlappedGWAS.bagged.uniq[[j]] <- i.oneBag[which.min(i.oneBag$pvalue),];
		j=j+1;
	}
}
overlappedGWAS.bagged.uniq <- do.call(rbind.data.frame, overlappedGWAS.bagged.uniq)

# merge two data.frames
overlappedGWAS.bagged <-NULL
if(nrow(overlappedGWAS.bagged.multi)>0){
	overlappedGWAS.bagged <- rbind(overlappedGWAS.bagged.multi, overlappedGWAS.bagged.uniq)
}else{
	overlappedGWAS.bagged <- overlappedGWAS.bagged.uniq
}


##3. Compare the size of GWAS and eQTL and conduct hypergeometric test

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

hyperPvalue[i]<-1-phyper(x, m, n, k, lower.tail = T)
c(x,m,n,k)
}


names(hyperPvalue)=sub("_lessThan10-5.chrPosition.pvalue", "", GWASFileList)
sort(hyperPvalue)


write.csv(sort(hyperPvalue), file="C://Users//User//Desktop//eQTL_enrichment//data//result_GWAS_10-8-10.csv")





##4. Print out results
