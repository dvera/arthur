repliSig <- function( bedGraphs, clusters=5, subsetMethod="sd", subsetParameter="auto", excludeFromFiltering=NULL, exludeFromClustering=NULL, pairing=NULL, merge=NULL, heatmapLims=c(-5,5), heatmapColors=c("red","black","green"), threads=getOption("threads",1L), ... ){
  require(travis)
  bgs <- tsvRead(bedGraphs,threads=threads)
  bgl <- data.matrix(as.data.frame(lapply(bgs,"[",4)))

  if( !is.null(pairing) && length(unique(unlist(lapply(pairing,length)))) != 1 ){
    stop("pairing must be a list of indices of each replicate")
  }
  if( !is.null(excludeFromFiltering)){
    pairing = lapply(1:length(pairing),function(p){
      pairing[[p]][which(pairing[[p]] %ni% excludeFromFiltering)]
    })
  }
  if( !is.null(pairing) && length(unique(unlist(lapply(pairing,length)))) != 1 ){
    stop("pairing must be a list of indices of each replicate")
  }

  if(subsetMethod="maxDiff"){

    if(!is.numeric(subsetParameter)){stop("subsetParameter must be a number specifying the minimum-consistent maximum difference in RT")}

    maxdiffs=as.data.frame(lapply(1:length(pairing),function(p){
      apply(bgl[,pairing[[p]]],1,max)-apply(bgl[,pairing[[p]]],1,min)
    }))

    selectRows <- which(apply(maxdiffs,1,min) >= subsetParameter)

    bgls <- bgl[selectRows]

  }

  k <- kmeans(bgls,clusters)$cluster
  ko <- order(k)
  kok <- k[ko]
  bglso <- bgls[ko]
  brks <- which(diff(kok)==1)+1


  heatmap.2(bglso,Rowv=F,trace='none',breaks=seq(heatmapLims[1],heatmapLims[2],(heatmapLims[2]-heatmapLims[1])/100) , col=colorRampPalette(heatmapColors),margins=c(15,5) )

}
