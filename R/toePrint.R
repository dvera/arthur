toePrint <- function( bgFilesList, removeSuffix=NULL, threads=getOption("threads",1L) )

require(travis)

numFiles <- lapply(bgFilesList,length)
numSets <- length(bgFilesList)
fileIndex <- rep(1:numSets,numFiles)
allFiles <- unlist(bgFilesList)
if(any(grepl(" ",allFiles))){stop("I hate spaces in file names, quiting.")}
setNames <- names(bgFilesList)
cat("counting file lines...\n")
allFileLines <- filelines(allFiles,threads=threads)

if(!is.list(bgFilesList) || numSets<2){
  stop("ERROR: bgFilesList should be a list of vectors of paths to bedGraph files")
}

if(length(unique(setNames))!=numSets){
  stop("must name lists of files uniquely")
}

# check if files have consistent rows
if(length(unique(allFileLines))>1){
  cat("files have inconsistent rows, unifying files\n")
  allFiles <- bgUnify(allFiles)
  newFileLines <- filelines(allFiles,threads=threads)
  bgFilesList <- lapply(1:numSets,function(x){ allFiles[which(fileIndex==x)] })
  print(data.frame(group=fileIndex,file=allFiles,linesBeforeUnifying=allFileLines,linesAfterUnifying=newFileLines,percentRemaining=100*newFileLines/allFileLines))
}

# define comparisons
i=1:length(bgFilesList)
eg=unique(t(apply(expand.grid(i,i),1,sort)))
eg=eg[which(eg[,1]!=eg[,2]),,drop=F]
nm=eg
nm[] <- setNames[nm]
numComps <- nrow(eg)

# file names
compNames <- paste0(nm[,1],"-vs-",nm[,2])
bedNames <- paste0(compNames,".bed")
bgNames <- paste0(compNames,"_phred.bg")
pdfNames <- paste0(compNames,".pdf")



#cat("reading ",numFiles," bedGraph files...\n")
bgs <- tsvRead(allFiles,threads=threads)
bgl <- as.data.frame(lapply(bgs,"[",4))
numRows <- nrow(bgl)
names(bgl) <- basename(removeext(allFiles))
names(bgs) <- setNames
sleepTimes <- sort(1:numComps)/100

wtf=lapply(1:numComps,function(compNum){


  Sys.sleep(sleepTimes[compNum])
  cat(compNames[compNum],": starting...\n")

  group1=which(fileIndex==eg[compNum,1])
  group2=which(fileIndex==eg[compNum,2])

  bgl1=bgl[,group1,drop=F]
  bgl2=bgl[,group2,drop=F]

  bga1=rowMeans(bgl1)
  bga2=rowMeans(bgl2)

  numSamples1=length(group1)
  numSamples2=length(group2)

  diag1=diag(matrix(1:(numSamples1^2),nrow=numSamples1))
  diag2=diag(matrix(1:(numSamples2^2),nrow=numSamples2))

  index1=1:numSamples1
  index2=(numSamples1+1):sum(numSamples1+numSamples2)

  numDists1=((numSamples1^2)-numSamples1)
  numDists2=((numSamples2^2)-numSamples2)
  numDists=numSamples1*numSamples2

  pb <- txtProgressBar(min = 0, max = numRows, style = 3)

  results <- do.call(rbind,mclapply(1:numRows,function(x){
#    results <- as.data.frame(t(as.data.frame(lapply(1:numRows,function(x){
    setTxtProgressBar(pb, x)
    allDistances <-  as.matrix(dist(bgl[x,c(group1,group2),drop=T]))
    withinDist1  <-  as.vector(allDistances[index1,index1])[-diag1]
    withinDist2  <-  as.vector(allDistances[index2,index2])[-diag2]
    acrossDist   <-  as.vector(allDistances[index1,index2])
    within_mean1 <- mean(withinDist1)
    within_mean2 <- mean(withinDist2)
    across_mean  <- sum(acrossDist)  / numDists
    #withinMeanMax <- max(c(within_mean1,within_mean2))
    withinMeanMean <- mean(c(within_mean1,within_mean2))
    localPvalue <- t.test(c(as.vector(withinDist1),as.vector(withinDist2)),as.vector(acrossDist),alternative="l")$p.value
    scorePvalue <- t.test(as.vector(bgl[x,group1]),as.vector(bgl[x,group2]),alternative="two.sided")$p.value
    deltaScore= mean(unlist(as.vector(bgl[x,group1,drop=T]))) - mean(unlist(as.vector(bgl[x,group2,drop=T])))
    
    
    res <- data.frame(
      rowNum         = x,
      within_mean1   = within_mean1,
      within_mean2   = within_mean2,
      acrossMeanDist = across_mean,
      #withinMeanMax  = withinMeanMax,
      withinMeanMean = withinMeanMean,
      localPvalue    = localPvalue,
      within = paste(c(withinDist1,withinDist2),collapse=","),
      across = paste(acrossDist,collapse=","),
      scorePvalue    = scorePvalue,
      deltaScore = deltaScore
      stringsAsFactors=F
    )
    return(res)
    
  },mc.cores=40,mc.preschedule=T))
  rownames(results)=NULL
  cat("\n")

  #results <- results[which(results$localPvalue <=0.05),]
  #results$localQvalue=p.adjust(results$localPvalue)
  withinMeanCdf <- ecdf(results$withinMeanMean)
  withinDists=as.numeric(unlist(lapply(results$within,strsplit,",")))
  withinCdf <- ecdf(withinDists)
  acrossDists=as.numeric(unlist(lapply(results$across,strsplit,",")))
  results$globalPvalue <- unlist(lapply(results$acrossMeanDist,withinMeanCdf))
  results$globalPvalue <- 1-results$globalPvalue
  results$globalQvalue <- p.adjust(results$globalPvalue,method="fdr")
  results$scoreQvalue  <- p.adjust(results$scorePvalue,method="fdr")
  #results$localQvalue  <- p.adjust(results$localPvalue)

  www=as.numeric(unlist(strsplit(paste(results$within,collapse=","),",")))


  rgbm <- matrix(rep(0,(3*numRows)),ncol=3)
  lowLocQ <- which(results$localQvalue  <= 0.05)
  lowGloQ <- which(results$globalQvalue <= 0.05)
  lowGloP <- which(results$globalPvalue <= 0.05)

  rgbm[lowLocQ,1]=0.7
  rgbm[lowGloP,3]=0.7
  colors=unlist(lapply(1:numRows,function(r){ do.call(rgb,as.list(rgbm[r,])) }))



  cat(compNames[compNum],": finished identifying regions, starting to save output...\n")

  pdf(file=pdfNames)
  rageHist(bgl[,c(group1,group2)])
  plot()


})
}
