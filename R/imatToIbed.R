imatToIbed <- function( imat , binsize , minInteractions=1 , threads=getOption("threads",1L) ){

  require(Matrix)
  options(scipen=99999)
  
  if(any(is.null(names(imat)))){ stop("matrices in supplied imat list must have names")}
  outnames <- paste0(names(imat),"_w",binsize,".ibed")
  numfiles <- length(imat)
  
  
triplet <- lapply(1:numfiles,function(f){
  chroms <- names(imat[[f]])
  numchroms <- length(chroms)
  triplets <- lapply(imat[[f]],summary)
  numbins <- unlist(lapply(triplets,nrow))
  chromcol <- unlist(lapply(1:numchroms,function(x) rep(chroms[x],numbins[x])))
  names(triplets) <- chroms
  triplet <- do.call(rbind,triplets)
  triplet <- cbind(chromcol,triplet,stringsAsFactors=F)
  triplet <- triplet[which(abs(triplet[,4])>=minInteractions),]
  col_chrom  = triplet[,1]
  col_lstart = triplet[,2]*binsize
  col_lstop  = triplet[,2]*binsize+binsize
  col_rstart = triplet[,3]*binsize
  col_rstop  = triplet[,3]*binsize+binsize
  col_count  = triplet[,4]
  #self <- which(col_lstart==col_rstart)
  triplet1 <- data.frame(
      V1=col_chrom,
      V2=col_lstart,
      V3=col_lstop,
      V4=paste0(col_chrom,":",col_rstart,"-",col_rstop,",",col_count),
      V5=1:length(col_chrom),
      stringsAsFactors=FALSE
    )
  triplet2 <- data.frame(
      V1=col_chrom,
      V2=col_rstart,
      V3=col_rstop,
      V4=paste0(col_chrom,":",col_lstart,"-",col_lstop,",",col_count),
      V5=length(col_chrom)+(1:length(col_chrom)),
      stringsAsFactors=FALSE
  )
  triplet <- rbind(triplet1,triplet2)
  triplet <- triplet[order(triplet$V1,triplet$V2),]
  write.table(triplet,outnames[f],sep="\t",quote=F,row.names=F,col.names=F)
  bgz <- bgzip(outnames[f])
  tabix(bgz,"bed")
  return(bgz)
})
  
return(triplet)
}
