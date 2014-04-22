rpm2timing3 <-
function ( bgfiles, reference=1,cores="max"){
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	
	
	fractions=c("G1b","G2","S1","S2","S3","S4")
	
	cat("checking file count\n")
	if(floor(length(bgfiles)/6) != ceiling(length(bgfiles)/6) ){ stop("# of files not divisible by 6. missing a file?") }
	
	#cat("checking file lines\n")
	#if(length(unique(unlist(lapply(bgfiles,filelines)))) != 1){stop("files have different numbers of lines")}
	
	cat("finding files for each cell line\n")
	cells<-unique(remove.suffix(remove.prefix(bgfiles[grep("G2",bgfiles)],"Seq"),"G2"))
	print(cells)
	numcells<-length(cells)
	rb=rainbow(numcells)
	cellfiles<-lapply(1:numcells,function(x) bgfiles[grep(cells[x],bgfiles)] )
	#print(cellfiles)
	
	ubgnames<-paste(cells,"_union.bg",sep="")
	outnames<-paste(cells,"_weightedscore.bg",sep="")
	finalnames<-paste(cells,"_iqr_qnorm.bg",sep="")
	#print(finalnames)
	
	cat("checking files for each cell line\n")
	numpercell=unique(unlist(lapply(cellfiles,length)))
	if(length(numpercell) != 1 | numpercell[1] != 6 ){print(cellfiles);stop("each cell line must match only 6 files")}
	checko<-lapply(1:numcells,function(x){
		lapply(1:6, function(y){
			if(grepl(fractions[y],cellfiles[[x]][y])==FALSE){
				stop(paste("cell files do not match expected fraction for",cells[x],cellfiles[[x]][y]))
			}
		} )
	})
	
	cat("making unionbg\n")
	#celldata<-lapply(1:numcells,function(x){
	celldata<-mclapply(1:numcells,function(x){
		cd<-read.delim(pipe(paste("bedtools unionbedg -filler NA -i",paste(cellfiles[[x]],collapse=" "))))
		colnames(cd)<-c("chrom","start","stop",fractions)
		goodrows<-which(rowSums(is.na(cd)) != 6)
		numrows<-nrow(cd)
		numgoodrows<-length(goodrows)
		cat(numrows-numgoodrows,"windows have no data in at least 1 fraction (",(numrows-numgoodrows)/numrows,"% )\n")
		cd[is.na(cd)]<-0
		wa<-(0.917*cd$G1b)+(0.750*cd$S1)+(0.583*cd$S2)+(0.417*cd$S3)+(0.250*cd$S4)+(0*cd$G2)
		cd<-data.frame("V1"=cd$chrom,"V2"=cd$start,"V3"=cd$stop,"V4"=wa,stringsAsFactors=FALSE)
		rm(wa)
		return(cd)
		print(head(cd))
	},mc.cores=cores)
	#})
	ref<-unlist(lapply(1:numcells,function(x) celldata[[x]][,4] ))
	#numrefs<-length(ref)
	
	cat("normalizing data\n")
	celldata<-mclapply(1:numcells,function(x){
		cd<-celldata[[x]]
		numpoints<-nrow(cd)
	#	if(numpoints < numrefs){
			curref<-sample(ref,numpoints)
			cd[,4][order(cd[,4])]<-curref[order(curref)]
	#	}
	#	if(numpoints > numrefs){
	#		curref<-sample(rep(ref,10),numpoints)
	#		cd[,4][order(cd[,4])]<-curref[order(curref)]
	#	}
	#	if(numpoints == numrefs){
	#		curref<-ref
	#		cd[,4][order(cd[,4])]<-curref[order(curref)]
	#	}
		cd[,4]<-( (cd[,4]-median(cd[,4],na.rm=TRUE) ) ) * 1.59  / IQR(cd[,4],na.rm=TRUE)
		write.tsv(cd,file=finalnames[x])
		return(cd)
	},mc.cores=cores)
}
