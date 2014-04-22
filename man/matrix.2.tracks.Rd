\name{matrix.2.tracks}
\alias{matrix.2.tracks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
matrix.2.tracks(matdir, genomefile, outname, windowsize = 10000, maxdist = windowsize * 40, cores = "max")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matdir}{
%%     ~~Describe \code{matdir} here~~
}
  \item{genomefile}{
%%     ~~Describe \code{genomefile} here~~
}
  \item{outname}{
%%     ~~Describe \code{outname} here~~
}
  \item{windowsize}{
%%     ~~Describe \code{windowsize} here~~
}
  \item{maxdist}{
%%     ~~Describe \code{maxdist} here~~
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
function (matdir, genomefile, outname, windowsize = 10000, maxdist = windowsize * 
    40, cores = "max") 
{
    library(parallel)
    library(gtools)
    if (cores == "max") {
        cores = detectCores() - 1
    }
    cat("finding chromosomes\n")
    matfiles <- mixedsort(files(paste0(matdir, "/*.mat")))
    chroms <- unique((remove.suffix(basename(matfiles), "_")))
    numchroms <- length(chroms)
    chromsizes <- read.tsv(genomefile)
    if (sum(chroms \%in\% chromsizes[, 1]) != length(chroms)) {
        stop("one or more chromosomes not found in genome file")
    }
    chromsizes <- chromsizes[match(chroms, chromsizes[, 1]), 
        ]
    numdistances = floor(maxdist/windowsize)
    d <- (1:numdistances) * windowsize
    distancenames <- formatC(d, width = nchar(maxdist), format = "d", 
        flag = "0")
    outnames <- lapply(1:numchroms, function(o) paste0(matdir, 
        "/", outname, "_", chroms[o], "_", distancenames, ".bg"))
    foutnames <- paste0(matdir, "/", outname, "_", distancenames, 
        ".bg")
    for (b in 1:numchroms) {
        mat <- read.tsv(matfiles[b])
        matcols <- ncol(mat)
        if (ceiling(matcols/2) != floor(matcols/2)) {
            mat <- mat[1:(matcols - 1), 1:(matcols - 1)]
        }
        matcols <- ncol(mat)
        l1 <- lapply(1:matcols, function(x) 1:x)
        l2 <- lapply(1:matcols, function(x) matcols + 1 - (x:1))
        rotmat <- matrix(NA, nrow = matcols, ncol = matcols)
        starts <- (matcols/2) - floor((0:(matcols - 1)/2))
        cat("rotating", chroms[b], "matrix\n")
        fillers <- mclapply(1:matcols, function(i) {
            v <- vector(length = i)
            for (j in 1:i) {
                v[j] <- mat[l1[[i]][j], l2[[i]][j]]
            }
            return(v)
        }, mc.cores = cores, mc.preschedule = F)
        for (i in 1:matcols) {
            rotmat[i, starts[i]:(starts[i] + i - 1)] <- fillers[[i]]
        }
        rotmat[is.na(rotmat)] <- 0
        oddstarts <- (0:(matcols - 1)) * windowsize + 1
        oddstops <- oddstarts + windowsize
        evenstarts <- oddstarts - windowsize/2
        evenstops <- evenstarts + windowsize
        md <- which(d > maxdist)[1] - 1
        for (i in ((1:(numdistances/2)) * 2) - 1) {
            df <- data.frame(chroms[b], evenstarts, evenstops, 
                rotmat[(matcols + 1 - i), ])
            write.tsv(df, file = outnames[[b]][i])
        }
        for (i in ((1:(numdistances/2)) * 2)) {
            df <- data.frame(chroms[b], oddstarts, oddstops, 
                rotmat[(matcols + 1 - i), ])
            write.tsv(df, file = outnames[[b]][i])
        }
    }
    for (b in 1:numdistances) {
        infiles <- paste(unlist(lapply(outnames, "[", b)), collapse = " ")
        system(paste("cat", infiles, "| awk '($2 > 0)' >", foutnames[b]))
        system(paste("bedGraphToBigWig", foutnames[b], genomefile, 
            gsub(".bg", ".bw", foutnames[b])))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
