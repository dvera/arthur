coords.2.matrix <-
function ( coordfile , genomefile , outdir , windowsize = 10000 , cores="max" ){
	
	library(parallel)
	library(gtools)
	if(cores=="max"){cores=detectCores()-1}
	
	cat("finding chromosomes\n")
	chroms <-  mixedsort (readLines ( pipe ( paste ( "cut -f 1" , coordfile , "| sort | uniq"))))
	chroms <- chroms[which(chroms != "chrM")]
	chroms <- chroms[which(chroms != "chrY")]
	numchroms <- length(chroms)
	chromsizes <- read.tsv ( genomefile )
	if(sum(chroms %in% chromsizes[,1]) != length(chroms)){stop("one or more chromosomes not found in genome file")}
	chromsizes<-chromsizes[match(chroms,chromsizes[,1]),]

	dir.create(outdir)
	outnames<-paste0(outdir,"/",chroms,"_win",windowsize,".mat")
		

	mclapply(1:numchroms,function(r) {
		cat("getting",chroms[r],"data\n")
		numchromwins<-floor(chromsizes[r,2]/windowsize)
		brks <- (0:numchromwins)*windowsize
		coords <- read.delim ( pipe ( paste0 ( "awk '($1==\"",chroms[r],"\" && $3==\"",chroms[r],"\")' ",coordfile," | cut -f 2,4" ) ) , stringsAsFactors=F , header = F )
		coords <- coords[which(coords[,1] <= brks[length(brks)] & coords[,2] <= brks[length(brks)] ),]
		a=coords[,1]
		b=coords[,2]
		
	
		cat("binning x-axis data\n")
		xbinind<-cut(a,brks,labels=F)
		
		cat("binning y-axis data\n")
		h<-mclapply(1:numchromwins,function(x){
			hist(b[which(xbinind==x)],breaks=brks,plot=F)
		},mc.cores=cores)
		
		cat("creating interaction matrix\n")
		binmat<-data.matrix(as.data.frame(lapply(h,"[[",2)))

		cat("saving binmat\n")
		write.tsv(binmat,file=outnames[r])
	} , mc.cores=cores , mc.preschedule=F )
}
