putGWAS.SNPs.into.Bags <-function(inputSubsetGWAS.ChrPos, inputWholeGWAS.ChrPos, input.bagSNP.all){ 
	#overlappedGWAS <- GWAS[match(intersect_sigeQTL_GWAS, GWAS$CHR.POS),]
	overlappedGWAS <- GWAS[match(inputSubsetGWAS.ChrPos, inputWholeGWAS.ChrPos),]
	dim(overlappedGWAS)
	overlappedGWAS.bag <- NULL
	for(i in 1:nrow(overlappedGWAS)){
		cat(paste(round(1/nrow(overlappedGWAS)*100,2),"%\n",sep=""))
		#res <- lapply(sig.bagSNP.all, function(ch) grep(i, ch))
		res <- lapply(input.bagSNP.all, function(ch) match(overlappedGWAS[i,1], ch))
		#bagIndex <- which(sapply(res, function(x) length(x) > 0))
		bagIndex <- which(!is.na(res))
		overlappedGWAS.bag <- c(overlappedGWAS.bag, bagIndex)
	}
	length(overlappedGWAS.bag)
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
	dim(overlappedGWAS.bagged.multi)
	head(overlappedGWAS.bagged.multi)
	
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
	dim(overlappedGWAS.bagged.uniq)
	head(overlappedGWAS.bagged.uniq)
	
	# merge two data.frames
	overlappedGWAS.bagged <-NULL
	if(nrow(overlappedGWAS.bagged.multi)>0){
		overlappedGWAS.bagged <- rbind(overlappedGWAS.bagged.multi, overlappedGWAS.bagged.uniq)
	}else{
		overlappedGWAS.bagged <- overlappedGWAS.bagged.uniq
	}
	dim(overlappedGWAS.bagged)
	head(overlappedGWAS.bagged)
	return(overlappedGWAS.bagged)
}#END function

