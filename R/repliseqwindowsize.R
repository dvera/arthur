repliseqwindowsize <-
function( early , late , scalar="auto" , genome="hg19" , stepsize=1000 , windowsizes = 1000 * 2 ^(0:10) , cores="max"){
	
	library(parallel)
	if(cores=="max"){cores=detectCores()-1}
	
	numwins<-length(windowsizes)
	rb<-rainbow(numwins)
	
	cat("windowing genome\n")
	windowfiles<-unlist(mclapply(1:numwins,function(x) bed.makewindows("hg19",genome=TRUE,windowsize=windowsizes[x],stepsize=stepsize,outname=paste(early,windowsizes[x],"windows",sep="")) , mc.cores=cores ))
	
	cat("counting windows\n")
	windowcounts<-unlist(mclapply(windowfiles,filelines,mc.cores=cores))
	print(windowcounts)
	
	# SORT BEDS!
	
# 	cat("calculating genome coverage\n")
# 	gcovs<-unlist(mclapply(c(early,late),bed.genomecov , genome=genome , scalar=scalar , makebigwig=FALSE , mc.cores=2 ))
# 	
# 	cat("calculating windowed coverage for early\n")
# 	ecovs<-unlist(mclapply(1:numwins,function(x) bg.window(gcovs[1], windowbed=windowfiles[x] , windowsize=windowsizes[x] , stepsize=stepsize , genome=genome , operation="mean") , mc.cores=cores ))
# 	lcovs<-unlist(mclapply(1:numwins,function(x) bg.window(gcovs[2], windowbed=windowfiles[x] , windowsize=windowsizes[x] , stepsize=stepsize , genome=genome , operation="mean") , mc.cores=cores ))
# 	
# 	cat("loading data\n")
# 	covs<-mclapply(ecovs,read.tsv,mc.cores=cores)
# 	
# 	covs<-mclapply(1:numwins,function(x) { 
# 		covs[[x]]$V5<-as.numeric(readLines(pipe(paste("cut -f 4",lcovs[x]))))
# 		covs[[x]]$V6<-log2(covs[[x]][,4]/covs[[x]][,5])
# 		covs[[x]]$V6[is.infinite(covs[[x]][,6])]<-NA
# 		covs[[x]]
# 	} , mc.cores=cores )
	
	cat("calculating windowed coverage for early\n")
	ecovs<-unlist(mclapply(1:numwins,function(x) bed.windowcov(early, windowbed=windowfiles[x] , windowsize=windowsizes[x] , stepsize=stepsize , genome=genome , scalar=scalar) , mc.cores=cores ))
	
	cat("calculating windowed coverage for late\n")
	lcovs<-unlist(mclapply(1:numwins,function(x) bed.windowcov(late, windowbed=windowfiles[x] , windowsize=windowsizes[x] , stepsize=stepsize , genome=genome , scalar=scalar) , mc.cores=cores ))
	
	ecovs<-gsub("bw","bg",ecovs)
	lcovs<-gsub("bw","bg",lcovs)
	
	cat("loading coverages\n")
	cat("loading data\n")
	covs<-mclapply(ecovs,read.tsv,mc.cores=cores)
	
	covs<-mclapply(1:numwins,function(x) { 
		covs[[x]]$V5<-as.numeric(readLines(pipe(paste("cut -f 4",lcovs[x]))))
		covs[[x]]$V6<-log2(covs[[x]][,4]/covs[[x]][,5])
		covs[[x]]$V6[is.infinite(covs[[x]][,6])]<-NA
		covs[[x]]
	} , mc.cores=cores )
	
	
	cat("counting and plotting zeroes\n")
	earlyzeroes<-(mclapply(1:numwins,function(x) which(covs[[x]][,4] == 0 ), mc.cores=cores ))
	numearlyzeroes<-unlist(lapply(earlyzeroes,length))
	
	latezeroes<-(mclapply(1:numwins,function(x) which(covs[[x]][,5] == 0 ), mc.cores=cores ))
	numlatezeroes<-unlist(lapply(latezeroes,length))
	
	bothzeroes<-(mclapply(1:numwins,function(x) which(covs[[x]][,4] == 0 & covs[[x]][,5] == 0) , mc.cores=cores ))
	numbothzeroes<-unlist(lapply(bothzeroes,length))
	
	pdf(file=paste(basename(removeext(early)),"_noscores-vs-windowsize.pdf",sep=""))
	plot(0,type="n",xlab="windowsize (bp)",ylab="% windows with no reads",xlim=c(0,128000) , ylim=c(0,50) )
	abline(v=windowsizes,col="grey70")
	points(windowsizes,100*numearlyzeroes/windowcounts,col="blue",lwd=3)
	lines(windowsizes,100*numearlyzeroes/windowcounts,col="blue",lwd=3)
	lines(windowsizes,100*numlatezeroes/windowcounts,col="red",lwd=3)
	points(windowsizes,100*numlatezeroes/windowcounts,col="red",lwd=3)
	lines(windowsizes,100*numbothzeroes/windowcounts,lwd=3)
	points(windowsizes,100*numbothzeroes/windowcounts,lwd=3)
	legend("topright",legend=c("early","late","both"),col=c("blue","red","black"),lwd=3)
	
	
	cat("calculating score distributions and plotting\n")
	equants<-unlist(mclapply(1:numwins,function(x) quantile(coverages[[x]][,4],probs=0.95,na.rm=TRUE) , mc.cores=cores ))
	equants<-1000*equants/windowsizes
	
	lquants<-unlist(mclapply(1:numwins,function(x) quantile(coverages[[x]][,5],probs=0.95,na.rm=TRUE) , mc.cores=cores ))
	lquants<-1000*lquants/windowsizes
	
	maxscore<-max(c(equants,lquants))
	
	earlydensities<-mclapply(1:numwins,function(x) density(1000*coverages[[x]][,4]/windowsizes[x],na.rm=TRUE,from=0,to=maxscore) , mc.cores=cores)
	latedensities<-mclapply(1:numwins,function(x) density(1000*coverages[[x]][,5]/windowsizes[x],na.rm=TRUE,from=0,to=maxscore) , mc.cores=cores)
	ratiodensities<-mclapply(1:numwins,function(x) density(coverages[[x]][,6],na.rm=TRUE,from=-8,to=8) , mc.cores=cores)
	
	
	plot(0,type="n",xlab="RPM",ylab="frequency of windows (kernel density estimate)",main="read density histograms for early fraction",xlim=c(0,quantile(unlist(lapply(earlydensities,"[","x")),probs=0.9)) , ylim=c(0,max(unlist(lapply(earlydensities,"[","y")))) )
	for(i in 1:numwins){
		lines(earlydensities[[i]][["x"]],earlydensities[[i]][["y"]],lwd=3,col=rb[i])
	}
	legend("topright",legend=paste(windowsizes,"bp windows") , col=rb , lwd=3, lty=1)
	
	
	plot(0,type="n",xlab="RPM",ylab="frequency of windows (kernel density estimate)",main="read density histograms for late fraction",xlim=c(0,quantile(unlist(lapply(earlydensities,"[","x")),probs=0.9)) , ylim=c(0,max(unlist(lapply(earlydensities,"[","y")))) )
	for(i in 1:numwins){
		lines(latedensities[[i]][["x"]],latedensities[[i]][["y"]],lwd=3,col=rb[i])
	}
	legend("topright",legend=paste(windowsizes,"bp windows") , col=rb , lwd=3, lty=1)
	
	
	plot(0,type="n",xlab="RT score , log2(early/late)",ylab="frequency of windows (kernel density estimate)",main="RT score histograms",xlim=c(-8,8) , ylim=c(0,max(unlist(lapply(ratiodensities,"[","y")))) )
	for(i in 1:numwins){
		lines(ratiodensities[[i]][["x"]],ratiodensities[[i]][["y"]],lwd=3,col=rb[i])
	}
	legend("topright",legend=paste(windowsizes,"bp windows") , col=rb , lwd=3, lty=1 , cex=0.75)
	
	latewherenoearly<-mclapply(1:numwins,function(x){
		coverages[[x]][,5][which(coverages[[x]][,4]==0)]
	},mc.cores=cores)
	
	earlywherenolate<-mclapply(1:numwins,function(x){
		coverages[[x]][,4][which(coverages[[x]][,5]==0)]
	},mc.cores=cores)
	dev.off()
	
	
}
