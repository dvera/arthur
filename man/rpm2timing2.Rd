\name{rpm2timing2}
\alias{rpm2timing2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rpm2timing2(bgfiles, reference = 1, cores = "max", fractionhist = TRUE, scorehist = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bgfiles}{
%%     ~~Describe \code{bgfiles} here~~
}
  \item{reference}{
%%     ~~Describe \code{reference} here~~
}
  \item{cores}{
%%     ~~Describe \code{cores} here~~
}
  \item{fractionhist}{
%%     ~~Describe \code{fractionhist} here~~
}
  \item{scorehist}{
%%     ~~Describe \code{scorehist} here~~
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
function (bgfiles, reference = 1, cores = "max", fractionhist = TRUE, 
    scorehist = TRUE) 
{
    library(parallel)
    if (cores == "max") {
        cores = detectCores() - 1
    }
    fractions = c("G1b", "G2", "S1", "S2", "S3", "S4")
    cat("checking file count\n")
    if (floor(length(bgfiles)/6) != ceiling(length(bgfiles)/6)) {
        stop("# of files not divisible by 6. missing a file?")
    }
    cat("checking file lines\n")
    if (length(unique(unlist(lapply(bgfiles, filelines)))) != 
        1) {
        stop("files have different numbers of lines")
    }
    cat("finding files for each cell line\n")
    cells <- unique(remove.suffix(remove.prefix(bgfiles[grep("G2", 
        bgfiles)], "Seq"), "G2"))
    print(cells)
    numcells <- length(cells)
    rb = rainbow(numcells)
    cellfiles <- lapply(1:numcells, function(x) bgfiles[grep(cells[x], 
        bgfiles)])
    print(cellfiles)
    outnames <- paste(cells, "_weightedscore.bg", sep = "")
    finalnames <- paste(cells, "_iqr_qnorm.bg", sep = "")
    cat("checking files for each cell line\n")
    numpercell = unique(unlist(lapply(cellfiles, length)))
    if (length(numpercell) != 1 | numpercell[1] != 6) {
        print(cellfiles)
        stop("each cell line must match only 6 files")
    }
    checko <- lapply(1:numcells, function(x) {
        lapply(1:6, function(y) {
            if (grepl(fractions[y], cellfiles[[x]][y]) == FALSE) {
                stop(paste("cell files do not match expected fraction for", 
                  cells[x], cellfiles[[x]][y]))
            }
        })
    })
    coords <- read.delim(pipe(paste("cut -f 1,2,3", cellfiles[[1]][1])), 
        header = F)
    cat("getting scores\n")
    celldata <- mclapply(1:numcells, function(x) {
        cd <- as.data.frame(lapply(1:6, function(y) {
            as.numeric(readLines(pipe(paste("cut -f 4", cellfiles[[x]][y]))))
        }))
        colnames(cd) <- fractions
        return(cd)
    }, mc.cores = cores)
    celldata <- as.data.frame(mclapply(1:numcells, function(x) {
        cd <- celldata[[x]]
        wa <- (0.917 * cd$G1b) + (0.75 * cd$S1) + (0.583 * cd$S2) + 
            (0.417 * cd$S3) + (0.25 * cd$S4) + (0 * cd$G2)
        wa[which(wa == 0)] <- NA
        return(wa)
    }, mc.cores = cores), stringsAsFactors = FALSE)
    colnames(celldata) <- cells
    numwindows <- nrow(celldata)
    goodblocks <- which(complete.cases(celldata))
    cat("discarding", numwindows - length(goodblocks), "windows that have no data in at least 1 sample (", 
        (numwindows - length(goodblocks))/numwindows, "\% )\n")
    celldata <- celldata[goodblocks, ]
    coords <- coords[goodblocks, ]
    reference <- rowMeans(as.data.frame(lapply(celldata, sort)))
    celldata <- as.data.frame(mclapply(1:numcells, function(x) {
        celldata[, x][order(celldata[, x])] <- reference
        return(celldata[, x])
    }, mc.cores = cores))
    celldata <- ((celldata - median(celldata, na.rm = TRUE))) * 
        1.59/IQR(unlist(celldata), na.rm = TRUE)
    celldata <- mclapply(1:numcells, function(x) {
        write.tsv(cbind(coords, celldata[, x], stringsAsFactors = FALSE), 
            file = finalnames[x])
    }, mc.cores = cores)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
