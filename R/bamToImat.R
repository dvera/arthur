bamToImat <- function(bedFile , binsize , minInteractions=1 , minQual=20 , threads=getOption("threads",1L) ){

  require(Matrix)
  options(scipen=99999)
  numfiles <- length(bedFile)
  
  cmdStrings <- paste(
    "bedtools bamtobed -bedpe -i",bedFile," | awk '{",
      "if($1==$4 && $8 >= ",minQual,"){",
        "left=int($2/",binsize,");",
        "right=int($6/",binsize,");",
        "if(left<right){print $1,left,right}",
        "else{print $1,right,left}",
      "}}' OFS='\t' | sort -T . -S 10G -k1,1 -k2,2n -k3,3n | uniq -c | awk '{",
        "if($1>=",minInteractions,"){",
          "print $2,$3,$4,$1",
        "}",
      "}' OFS='\t'"
  )

  res <- cmdRun(cmdStrings,threads=threads,tsv=T)
  
  res <- mclapply(1:numfiles,function(x){
    
    y=split(res[[x]],res[[x]][,1])
    yl=unlist(lapply(y,nrow))
    y=y[which(yl>100)]
    chroms <- names(y)
    y=lapply(y,"[",-1)
    y=lapply(y,data.matrix)
    matlist = lapply(seq_along(chroms),function(m){
      numbins <- max(y[[m]][,1:2])
      z=spMatrix(numbins,numbins,y[[m]][,1],y[[m]][,2],y[[m]][,3])
      return(z)
    })
    names(matlist) <- chroms
    return(matlist)
  },mc.cores=threads)
  
  
  
  names(res) <- bedFile
  return(matlist)
}
