\name{coords.2.matrix}
\alias{coords.2.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
coords.2.matrix(coordfile, genomefile, outdir, windowsize = 10000, cores = "max")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordfile}{
%%     ~~Describe \code{coordfile} here~~
}
  \item{genomefile}{
%%     ~~Describe \code{genomefile} here~~
}
  \item{outdir}{
%%     ~~Describe \code{outdir} here~~
}
  \item{windowsize}{
%%     ~~Describe \code{windowsize} here~~
}
  \item{cores}{
%%     ~~Describe \code{cores} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (coordfile, genomefile, outdir, windowsize = 10000, 
    cores = "max") 
{
    library(parallel)
    library(gtools)
    if (cores == "max") {
        cores = detectCores() - 1
    }
    cat("finding chromosomes\n")
    chroms <- mixedsort(readLines(pipe(paste("cut -f 1", coordfile, 
        "| sort | uniq"))))
    chroms <- chroms[which(chroms != "chrM")]
    chroms <- chroms[which(chroms != "chrY")]
    numchroms <- length(chroms)
    chromsizes <- read.tsv(genomefile)
    if (sum(chroms \%in\% chromsizes[, 1]) != length(chroms)) {
        stop("one or more chromosomes not found in genome file")
    }
    chromsizes <- chromsizes[match(chroms, chromsizes[, 1]), 
        ]
    dir.create(outdir)
    outnames <- paste0(outdir, "/", chroms, "_win", windowsize, 
        ".mat")
    mclapply(1:numchroms, function(r) {
        cat("getting", chroms[r], "data\n")
        numchromwins <- floor(chromsizes[r, 2]/windowsize)
        brks <- (0:numchromwins) * windowsize
        coords <- read.delim(pipe(paste0("awk '($1==\"", chroms[r], 
            "\" && $3==\"", chroms[r], "\")' ", coordfile, " | cut -f 2,4")), 
            stringsAsFactors = F, header = F)
        coords <- coords[which(coords[, 1] <= brks[length(brks)] & 
            coords[, 2] <= brks[length(brks)]), ]
        a = coords[, 1]
        b = coords[, 2]
        cat("binning x-axis data\n")
        xbinind <- cut(a, brks, labels = F)
        cat("binning y-axis data\n")
        h <- mclapply(1:numchromwins, function(x) {
            hist(b[which(xbinind == x)], breaks = brks, plot = F)
        }, mc.cores = cores)
        cat("creating interaction matrix\n")
        binmat <- data.matrix(as.data.frame(lapply(h, "[[", 2)))
        cat("saving binmat\n")
        write.tsv(binmat, file = outnames[r])
    }, mc.cores = cores, mc.preschedule = F)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
