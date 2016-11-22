bedpeToImat <- function(bedFile , binsize , prefix , minInteractions=1 , threads=getOption("threads",1L) ){

  require(Matrix)

  if(missing(prefix)){ prefix <- basename(removeext(bedFile)) }
  #tripletName <- paste0(basename(removeext(bedFile)),".triplet")


  chroms <- sort(as.character(unlist(cmdRun(paste("cut -f 1",bedFile,"| uniq"),lines=T))))
  chroms <- chroms[which(chroms!="chrY")]
  cmdString <- paste0(
    "awk '{if($1==\"",chroms,"\" && $4== \"",chroms,"\" ){print $1,int($2/",binsize,"),int($5/",binsize,")}}' ", bedFile," | awk '{if($2<$3){print $1,$2,$3} else{print $1,$3,$2}}' OFS='\t' | sort -T . -S 10G -k1,1 -k2,2n -k3,3n | uniq -c | awk '{if($1>=",minInteractions,"){print $3,$4,$1}}' OFS='\t'"
  )

  res <- cmdRun(cmdString,threads=threads,tsv=T)
  names(res)=chroms

  outnames <- paste0(prefix,"_",chroms,"_w",binsize,".imat")

  # convert to coord tables to matrices
  res <- lapply(res,data.matrix)

  matlist = mclapply(seq_along(chroms),function(m){
    numbins <- max(res[[m]][,1:2])
    spMatrix(numbins,numbins,res[[m]][,1],res[[m]][,2],res[[m]][,3])
  }, mc.cores=threads, mc.preschedule=FALSE)

  names(matlist) <- chroms
  return(matlist)
}
