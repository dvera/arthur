imatToIbed <- function( imat , binsize , prefix , minInteractions=1 , threads=getOption("threads",1L) ){

  require(Matrix)
  options(scipen=99999)

  chroms <- names(imat)

  outnames <- paste0(prefix,"_w",binsize,".ibed")

  triplet <- mclapply(1:length(chroms),function(x){
    s <- summary(imat[[x]])
    s <- s[which(s[,3]>=minInteractions),]
    s[,1]=s[,1]*binsize
    s[,2]=s[,2]*binsize
    data.frame(
      V1=rep(chroms[x],nrow(s)),
      V2=c(s[,1],s[,2]),
      V3=c(s[,1]+binsize,s[,2]+binsize),
      V4=paste0(rep(chroms[x],2),":",s[,2],"-",s[,2]+binsize,",",s[,3]),
      V5=1:(2*nrow(s)),
      V6="."
    )
  },mc.cores=threads,mc.preschedule=FALSE)

  triplet <- do.call(rbind,triplet)

  #triplet <- triplet[order(triplet$V1,triplet$V2,triplet$V3),]

  write.table(triplet,outnames,sep="\t",quote=F,row.names=F,col.names=F)
  outnames <- bedSort(outnames)
  bgz <- bgzip(outnames)
  tabix(bgz,"bed")
  return(bgz)
}
