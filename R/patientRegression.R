patientRegression <- function( bedGraphs, metadata, labelColumn=1, variableColumn=9, top=100, categoricalMetadata=FALSE, bedGraphIsPredictor=TRUE, patientNames=NULL, threads=getOption("threads",1L) ){

  # load data
  bgl = bgRead(bedGraphs,enforceEquality=TRUE)
  bg  = tsvRead(bedGraphs[1])
  vars = tsvRead(metadata,col_names=T)

  #if(removeBad){vars=vars[vars$Exclude2==0,]}

  # label scores to match metadata labels
  if(is.null(patientNames)){sampleNames=basename(removeext(bedGraphs))}
  colnames(bgl) <- patientNames

  # create matching sets of metadata and bedGraph scores
  bgToVars=match(colnames(bgl),vars[,labelColumn])
  bgToVars2=which(!is.na(bgToVars))
  vars <- vars[na.omit(bgToVars),]
  bgl  <- bgl[,bgToVars2]
  if(!identical(vars[,labelColumn],colnames(bgl))){stop("error parsing patients")}

  # convert categorical variables to factors
  if(categoricalMetadata){vars[,variableColumn]<-as.factor(vars[,variableColumn])}

  # perform regressions and extract pvals and slopes
  metaData <- vars[,variableColumn]

  if(bedGraphIsPredictor){
    glms <- lapply(1:nrow(bgl),function(i){
      bgscores <- bgl[i,]
      glm(metaData~bgscores, family = binomial(link = "logit"))
    })
  } else{
    glms <- lapply(1:nrow(bgl),function(i){
      bgscores <- bgl[i,]
      glm(bgscores~metaData, family = binomial(link = "logit"))
    })
  }

  pvals <- unlist(lapply(1:length(glms),function(i) summary(glms[[i]])$coefficients[2,4]))
  slopes <- unlist(lapply(1:length(glms),function(i) summary(glms[[i]])$coefficients[1,1]))
  signs <- sign(slopes)
  logp <- abs(-signs*log10(pvals))

  porder <- order(pvals,decreasing=F)
  varorder <- order(vars[,variableColumn])


  vars2 <- vars[varorder,]
  bgl2 <- bgl[porder,varorder]
  bg2 <- bg[porder,1:3]
  pvals2 <- pvals[porder]
  slopes2 <- slopes[porder]
  glms2 <- glms[porder]
  metaData2 <- metaData[varorder]

  pdf(file="patientRegression.pdf")
  hist(pvals,main="histogram of unadjusted p-values",xlim=c(0,1),breaks=500)
  abline(v=0.05,col="red")
  text(0.05,0,"0.05")
  hist(p.adjust(pvals,method="fdr"),main="histogram of fdr-adjusted p-values",xlim=c(0,1),breaks=500)
  #plot(density(slopes),main="histogram of model slopes")
  for(i in 1:top){
    bgscores2 <- bgl2[i,]
    plot(bgscores2,as.numeric(metaData2)-1,ylab=colnames(vars2)[variableColumn],xlab="RT",main=paste(paste(bg2[i,1:3],collapse="-"),", p=",pvals2[i]))
    curve(predict.glm(glms2[[i]],data.frame(bgscores=x),type="resp"),add=TRUE) # draws a curve based on prediction from logistic regression model
    #points(bgscores,fitted(glms2[[i]]),pch=20) # optional: you could skip this draws an invisible set of points of body size survival based on a 'fit' to glm model. pch= changes type of dots.

  }
  dev.off()


}
