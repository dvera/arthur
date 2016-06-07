mat2bg <- function( matlist , binsize , chromnames , prefix , maxbins=20 , threads=getOption("threads",1l) ){

  options(scipen=99999)

  numfiles <- length(matlist)

  threads1 <- floor(threads/numfiles)
  if(threads1<1){threads1=1}
  threads2 <- floor(threads/threads1)
  if(threads2<1){threads2=1}

  bins <- mclapply(1:numfiles, function(i){

    # read in matrix
    mat <- matlist[[i]]
    # mat <- Matrix(read_tsv(matfiles[i],col_names=FALSE))
    #mat <- data.matrix(read_tsv(matfiles[i],col_names=FALSE))

    #trim matrix to make square
    if(ncol(mat) != nrow(mat)){
      cat("WARNING: MATRIX",matlist[i],"IS NOT SQUARE, TRIMMING\n")
      mindim <- min(dim(mat))
      mat <- mat[ 1:mindim , 1:mindim ]
    }


    numbins <- ncol(mat)

    # lower maxbins if necessary
    if( maxbins > (numbins-1)){ maxbins <- numbins-1 }
    distwidth <- nchar(maxbins * binsize)


    flankleft  = rep( c(0:maxbins), each=2 )
    flankright = rotvec( flankleft )
    bincount <- length(diag(mat))

    y <- data.matrix(as.data.frame(mclapply(1:maxbins, function(w){
      c( rep(NA, flankleft[w] ) , diag(mat[,w:numbins]) , rep( NA, flankright[w] ) )
    },mc.cores=threads2,mc.preschedule=FALSE)))

    y[y==0] <- NA
    colnames(y) <- paste0("b",seq_len(maxbins))


    bgs <- mclapply(1:maxbins, function(w){

      bg <- data.frame(
        chromnames[i] ,
        ((0:bincount)*binsize+1)[1:bincount] ,
        (1:bincount)*binsize ,
        y[,w]
      )
      #bg[seq(2,bincount,2),2:3] <- bg[seq(2,bincount,2),2:3] - (binsize/2)
      if(w %in% seq(2,maxbins,2)){
        bg[,2:3] <- bg[,2:3] + (binsize/2)
      }

      bg <- bg[complete.cases(bg),]
      colnames(bg) <- paste0("V",1:4)
      bg
    }, mc.cores=threads2, mc.preschedule=FALSE)

    return(bgs)

  }, mc.cores=threads1, mc.preschedule=FALSE)

  dump<-mclapply(1:maxbins, function(w){
    # if(length(numfiles)>1){
      bg<-do.call(rbind,lapply(bins,"[[",w))
    # } else{
      # bg=bins[[w]]
    # }
    bg<-bg[order(bg[,1],bg[,2]),]
    outname<-paste0(prefix,"_d",sprintf(paste0("%0",distwidth,"d"), binsize*(w-1)),".bg")
    tsvWrite(bg,file=outname)
    outname
  }, mc.cores=threads1, mc.preschedule=FALSE)

  return(unlist(dump))

}
