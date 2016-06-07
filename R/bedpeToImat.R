bedpeToImat <- function(bedFile , binsize , prefix , threads=getOption("threads",1L) ){
  #  awk '{if($1==$4){print $1,int($2/100000),int($5/100000)}}' g1all.bedpe | sort -k1,1 -k2,2n -k3,3n | uniq -c | awk '{print $2,$3,$4,$1}' OFS='\t'  > g1all.triplet

  require(Matrix)

  if(missing(prefix)){ prefix <- basename(removeext(bedFile)) }
  bedpe  = tsvRead(bedFile)
  bedpe  = bedpe[which(bedpe[,1]==bedpe[,4]),]
  chroms = bedpe[,1]
  chromnames = unique(chroms)

  outnames <- paste0(prefix,"_",chromnames,"_w",binsize,".imat")

  # subset only coordinates
  bedpe  = bedpe[,c(2,5)]

  # flip columns so only one diagnal of matrix is filled
  flipme = which(bedpe[,1]>bedpe[,2])
  bedpe[flipme,] <- bedpe[flipme,c(2,1)]



  # split by chromosome
  bedpel = split.data.frame(bedpe,factor(chroms),drop=TRUE)
  rm(bedpe)

  # convert to coord tables to matrices
  bedpel = lapply(bedpel,data.matrix)

  dump <- mclapply(1:length(chromnames),function(i){
    s=seq(from=0,to=max(bedpel[[i]]),by=binsize)
    x1=findInterval(bedpel[[i]][,1],s,all.inside=T)
    x2=findInterval(bedpel[[i]][,2],s,all.inside=T)
    maxbin <- max(c(x1,x2))
    m=Matrix(0,nrow=maxbin,ncol=maxbin,sparse=TRUE)
    for( seq_along(x1) ){
      m[ x1[j] , x2[j] ] = m[ x1[j] , x2[j] ] + 1
    }

    write.table(as.matrix(m),file=outnames[i],sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)


  } , mc.cores=threads, mc.preschedule=FALSE)



}
