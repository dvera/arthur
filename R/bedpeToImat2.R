bamToImat <- function(bedFile , binsize , minInteractions=1 , minQual=20 , threads=getOption("threads",1L) ){

  require(Matrix)

  numfiles <- length(bedFile)
  
  chroms <- sort(as.character(unlist(cmdRun(paste("cut -f 1",bedFile[1],"| uniq"),lines=T))))
  chroms <- chroms[which(chroms!="chrY")]

  outnames <- lapply(1:numfiles,function(x){
    paste0(basename(removeext(bedFile[x])),"_q",minQual,"_min",minInteractions,"_",chroms,"_w",binsize,".imat")
  })

  
  cmdStrings <- paste(
    "bedtools bamtobed -bedpe -i",bedFile," | awk '{",
      "if($1==$4 && $8 >= ",minQual,"){",
        "left=int($2/",binsize,");",
        "right=int($6/",binsize,");",
        "if(left<right){print $1,left,right}",
        "else{print $1,right,left}",
      "}' OFS='\t' | sort -T . -S 10G -k1,1 -k2,2n -k3,3n | uniq -c | awk '{",
        "if($1>=",minInteractions,"){",
          "print $2,$3,$4,$1",
        "}",
      "}' OFS='\t'"
  )

  res <- cmdRun(cmdString,threads=threads,tsv=T)

  res <- mclapply(1:numfiles,function(x){
    y=split(res[[x]],res[[x]][,1])
    y=lapply(y,"[",-1)
    y=lapply(y,data.matrix)
    matlist = lapply(seq_along(chroms),function(m){
      numbins <- max(y[[m]][,1:2])
      spMatrix(numbins,numbins,y[[m]][,1],y[[m]][,2],y[[m]][,3])
    })
    return(matlist)
  },mc.cores=threads)
  
  names(res) <- bedFile
  return(matlist)
}
