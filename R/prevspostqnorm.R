prevspostqnorm <- function (ifiles, qfiles, pdfname = "prevspostqnorm.pdf" , samplenames = basename(removeext(ifiles)) , lspan=0.2, chrom="chr19" , chromrange = c(1,5000000 )){

	pdf(file=pdfname)
		

	for(l in 1:length(ifiles)) {

		# read in files
		i<-read.tsv(ifiles[l])
		q<-read.tsv(qfiles[l])
		
		# find rows within range
		ic<-which(i[,1]==chrom & i[,2] >= chromrange[1] & i[,2] <= chromrange[2])
		qc<-which(q[,1]==chrom & q[,2] >= chromrange[1] & q[,2] <= chromrange[2])

		# grab bases
		ib<-i[ic,2]
		qb<-q[qc,2]


		# grap data within range
		iu<-i[ic,4]
		qu<-q[qc,4]

		# smooth data
		il<-loess(iu~ib,span=lspan)$fitted
		ql<-loess(qu~qb,span=lspan)$fitted
		
		
		plot(qb,qu,type="l",col="lightblue",main=paste(samplenames[l],"RepliSeq before and after qnorm to pool"),xlab=paste(chrom,"coordinate (bp)"),ylab="repliseq score")
		lines(qb,ql,col="blue",lwd=3)
		lines(ib,iu,col="pink",lwd=3)
		lines(ib,il,col="red",lwd=3)
		legend("topleft",legend=c("post qnorm","post qnorm smoothed","pre  qnorm","pre  qnorm smoothed"),col=c("lightblue","blue","pink","red"),lwd=3)

		}

	dev.off()

}
