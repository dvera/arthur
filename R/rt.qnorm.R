rt.qnorm <- function(rtmatrix, cols = NA , excludestrings = NA , recenter=FALSE, setIQR=FALSE , iqr=1.59){

	options(scipen=9999)

	read.tsv <- function( tsv, ... ){
		read.table(tsv,stringsAsFactors=FALSE,sep="\t", ... )
	}

	write.tsv <- function( tsv, colnames=FALSE, rownames=FALSE, ... ){
		write.table(tsv,sep="\t",quote=FALSE,col.names=colnames,row.names=rownames, ... )
	}

	removeext <- function( filenames ){
		filenames<-as.character(filenames)
		for(i in 1:length(filenames)){
			namevector<-unlist(strsplit(filenames[i],"\\."))
			filenames[i]<-paste(namevector[1:(length(namevector)-1)],collapse=".")
		}
		filenames
	}
	'%ni%' <- Negate('%in%')

	filename<-basename(removeext(rtmatrix))

	cat("reading in matrix\n")
	mat<-read.tsv(rtmatrix,header=T)
	numcols<-ncol(mat)
	if(is.na(cols)){cols<-3:numcols}
	cat("processing matrix\n")
	matsums<-apply(mat[,cols],1,sum)
	goodrows<-which(matsums != 0)

	badcols<-unique(unlist(lapply(excludestrings,function(x) grep(x,colnames(mat)))))
	if(length(badcols)>1){cols<-cols[which(cols %ni% badcols)]}
	cat("excluding",colnames(mat)[badcols],"\n")
	cat("using",length(cols),"/",numcols,"columns\n")

	mat<-mat[goodrows,]
	smat<-mat[,cols]
	smat<-apply(smat,2,sort)
	vmat<-rowMeans(smat)
	vmat<-vmat[order(vmat)]
	nmat<-mat
	for(i in 3:numcols){
		nmat[order(nmat[,i]),i]<-vmat
	}

	if(recenter){
		nmat[,3:numcols]<-nmat[,3:numcols]-median(nmat[,3])

	}

	if(setIQR){
		nmat[,3:numcols]<-nmat[,3:numcols] * iqr / IQR(nmat[,3])
	}


	cat("saving normalized matrix\n")
	write.tsv(nmat,colnames=TRUE,file=paste0(basename(removeext(rtmatrix)),"_qnormToPlatformAvg.txt"))

	cat("calculating histogram limits\n")
	dl<-lapply(mat[,3:numcols], density)
	dlx<-lapply(dl,"[[",1)
	dly<-lapply(dl,"[[",2)
	ylims<-c(0,1.15*max(unlist(dly)))
	xrange<-ceiling(quantile(abs(unlist(dlx)),probs=0.95))
	xlims<-c(-xrange,xrange)
	if(min(unlist(dlx)) >= 0){xlims[1]=0}

	cat("plotting histograms\n")
	pdf(file=paste0(filename,"_histogram.pdf"))
	plot(0,type="n",ylim=ylims,xlim=xlims,xlab="RT score",ylab="density",main=paste("score histograms for",filename))
	for(i in 3:numcols){
		lines(density(mat[,i]),col=rgb(0,0,0,25,maxColorValue=255),lwd=2)
	}
	lines(density(vmat),lwd=3,col="red")
	legend("topright",legend=c("individual samples","average distribution"),col=c("grey50","red"),lwd=3)
	dev.off()

}
