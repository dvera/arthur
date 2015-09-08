\name{rpm2timing3}
\alias{rpm2timing3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rpm2timing3(bgfiles, reference = 1, cores = "max")
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
function (bgfiles, reference = 1, cores = "max") 
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
    cat("finding files for each cell line\n")
    cells <- unique(remove.suffix(remove.prefix(bgfiles[grep("G2", 
        bgfiles)], "Seq"), "G2"))
    print(cells)
    numcells <- length(cells)
    rb = rainbow(numcells)
    cellfiles <- lapply(1:numcells, function(x) bgfiles[grep(cells[x], 
        bgfiles)])
    ubgnames <- paste(cells, "_union.bg", sep = "")
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
    cat("making unionbg\n")
    celldata <- mclapply(1:numcells, function(x) {
        cd <- read.delim(pipe(paste("bedtools unionbedg -filler NA -i", 
            paste(cellfiles[[x]], collapse = " "))))
        colnames(cd) <- c("chrom", "start", "stop", fractions)
        goodrows <- which(rowSums(is.na(cd)) != 6)
        numrows <- nrow(cd)
        numgoodrows <- length(goodrows)
        cat(numrows - numgoodrows, "windows have no data in at least 1 fraction (", 
            (numrows - numgoodrows)/numrows, "\% )\n")
        cd[is.na(cd)] <- 0
        wa <- (0.917 * cd$G1b) + (0.75 * cd$S1) + (0.583 * cd$S2) + 
            (0.417 * cd$S3) + (0.25 * cd$S4) + (0 * cd$G2)
        cd <- data.frame(V1 = cd$chrom, V2 = cd$start, V3 = cd$stop, 
            V4 = wa, stringsAsFactors = FALSE)
        rm(wa)
        return(cd)
        print(head(cd))
    }, mc.cores = cores)
    ref <- unlist(lapply(1:numcells, function(x) celldata[[x]][, 
        4]))
    cat("normalizing data\n")
    celldata <- mclapply(1:numcells, function(x) {
        cd <- celldata[[x]]
        numpoints <- nrow(cd)
        curref <- sample(ref, numpoints)
        cd[, 4][order(cd[, 4])] <- curref[order(curref)]
        cd[, 4] <- ((cd[, 4] - median(cd[, 4], na.rm = TRUE))) * 
            1.59/IQR(cd[, 4], na.rm = TRUE)
        write.tsv(cd, file = finalnames[x])
        return(cd)
    }, mc.cores = cores)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
