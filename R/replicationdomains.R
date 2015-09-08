replicationdomains <-
function(bgfile , lspan = 0 , mergewithin = 100000 , removeY = T ){
	
	
	ttrname<-paste0(basename(removeext(bgfile)),"_ttr.bed")
	rdname<-paste0(basename(removeext(bgfile)),"_rds.bed")
	

	if(lspan != 0){
		cat(basename(bgfile),": smoothing bedGraph\n")
		r<-bg.loess(r,lspan=lspan)
	}
	
	cat(basename(bgfile),": loading bedGraph\n")
	r<-read.tsv(bgfile)
	rl<-nrow(r)
	#gf<-read.tsv(genomefile)
	
	if(removeY==TRUE){
		r<-r[which(r$V1 != "chrY"),]
	}

	cat(basename(bgfile),": finding TTRs\n")
	
	#find sign-changes
	r$prod=c(r$V4[1:(rl-1)]*r$V4[2:rl],1)
	
	#calculate slope at each point
	slopes=diff(r$V4)/diff(r$V2)
	r$slope=c(slopes,slopes[rl-1])
	rm(slopes)

	#find slop inversions
	r$slopeprod=c(r$slope[1],r$slope[1:(rl-1)]*r$slope[2:rl])
	slopeinvs<-which(r$slopeprod<0)
	signchanges<-which(r$prod<0)

	cat(basename(bgfile),": finding left TTR boundaries\n")
	leftttrb<-unlist(lapply(1:length(signchanges), function(x) signchanges[x]-which(r$slopeprod[signchanges[x]:1]<0)[1]+1 ))
	cat(basename(bgfile),": finding right TTR boundaries\n")
	ritettrb<-unlist(lapply(1:length(signchanges), function(x) signchanges[x]+which(r$slopeprod[(signchanges[x]):(signchanges[length(signchanges)])]<0)[1]-1 ))
	
	#t<-data.frame(signchanges,c(NA,leftttrb),c(ritettrb,NA))
	t<-data.frame(signchanges,leftttrb,ritettrb)

	cat(basename(bgfile),": removing incomplete TTRs at chromosome ends\n")
	cbreaks1 <- which(  r[t[,2],2] > r[t[,1],2]  )
	if(length(cbreaks1)>0){ t[cbreaks1,2] <- NA }
	
	cbreaks2 <- which(  r[t[,3],2] < r[t[,1],2]  )
	if(length(cbreaks2)>0){ t[cbreaks2,3] <- NA }

	cat(basename(bgfile),": finding TTR directionality\n")
	t$strand <- NA
	t[which(r[signchanges,4]<0),4]<-"-"
	t[which(r[signchanges,4]>0),4]<-"+"

	cat(basename(bgfile),": saving TTRs to ",ttrname,"\n")
	ttr<-data.frame("V1"=r[t[,1],1],"V2"=r[t[,2],2],"V3"=r[t[,3],2],"V4"=1:nrow(t),"V5"=1,"V6"=t[,4],"V7"=r[t[,1],2],"V8"=r[t[,1],2]+1)
	ttr<-ttr[which(complete.cases(ttr)),]
	write.tsv(ttr,file=ttrname)

	cat(basename(bgfile),": saving RDs to ",rdname,"\n")
	rd<-data.frame("V1"=r[t[1:(nrow(t)-1),1],1],"V2"=r[t[1:(nrow(t)-1),1],2],"V3"=r[t[2:(nrow(t)),1],2], "V4"=1:(nrow(t)-1) , "V5"=1 , "V6"=t[2:(nrow(t)),4] , "V7"=r[t[1:(nrow(t)-1),3],2],"V8"=r[t[2:(nrow(t)),2],3])
	rd<-rd[which(rd$V3>rd$V2),]
	rd<-rd[which(complete.cases(rd)),]
	write.tsv(rd,file=rdname)
	return(c(ttrname,rdname)
	#chromends<-prerd[which(prerd$V3<prerd$V2),]
	#chrombegs<-chromends
	#chrombegs[,2]<-1
	#chromends[,3]<-gf[match(chromends$V1,gf$V1),2]
		#otherchrom=r[t[2:(nrow(t)),1],1]
	#rd<-rd[which(rd$V1==otherchrom),]
	#rm(otherchrom)

}
