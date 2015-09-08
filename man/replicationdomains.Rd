\name{replicationdomains}
\alias{replicationdomains}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
replicationdomains(bgfile, lspan = 0, mergewithin = 1e+05, removeY = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bgfile}{
%%     ~~Describe \code{bgfile} here~~
}
  \item{lspan}{
%%     ~~Describe \code{lspan} here~~
}
  \item{mergewithin}{
%%     ~~Describe \code{mergewithin} here~~
}
  \item{removeY}{
%%     ~~Describe \code{removeY} here~~
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
function (bgfile, lspan = 0, mergewithin = 1e+05, removeY = T) 
{
    if (lspan != 0) {
        r <- bg.loess(r, lspan = lspan)
    }
    r <- read.tsv(bgfile)
    if (removeY == TRUE) {
        ttr <- ttr[which(ttr$V1 != "chrY"), ]
        r <- r[which(r$V1 != "chrY"), ]
    }
    r$Vprod = c(1, r$V4[1:(rl - 1)] * r$V4[2:rl])
    r$diff <- c(1, diff(r$V4))
    r$slope = c(1, diff(r$V4)/diff(r$V2))
    r$slopeprod = c(1, r$slope[1:(rl - 1)] * r$slope[2:rl])
    spc <- which(r$slopeprod < 0)
    pc <- which(r$prod < 0)
    ttrb1 <- unlist(lapply(1:length(pc), function(x) pc[x] - 
        which(r[pc[x]:1, 8] < 0)[1] + 1))
    ttrb2 <- unlist(lapply(1:length(pc), function(x) pc[x] + 
        which(r[pc[x]:pc[length(pc)], 8] < 0)[1] - 1))
    chroms <- sort(unique(r$V1))
    numchroms <- length(chroms)
    chromlines <- lapply(1:numchroms, function(x) which(r[, 1] == 
        chroms[x]))
    t <- data.frame(pc, ttrb1, ttrb2)
    t[which(r[t[, 2], 2] > r[t[, 1], 2]), 2] <- NA
    t[which(r[t[, 3], 2] < r[t[, 1], 2]), 3] <- NA
    t$strand <- NA
    t[which(r[pc, 4] < 0), 4] <- "+"
    t[which(r[pc, 4] > 0), 4] <- "-"
    ttr <- data.frame(V1 = r[t[, 1], 1], V2 = r[t[, 2], 2], V3 = r[t[, 
        3], 3], V4 = 1:nrow(t), V5 = 1, V6 = t[, 4], V7 = r[t[, 
        1], 2], V8 = r[t[, 1], 3])
    ttr <- ttr[which(complete.cases(ttr)), ]
    rd <- data.frame(V1 = r[t[1:(nrow(t) - 1), 1], 1], V2 = r[t[1:(nrow(t) - 
        1), 1], 2], V3 = r[t[2:(nrow(t)), 1], 2], V4 = 1:(nrow(t) - 
        1), V5 = 1, V6 = t[2:(nrow(t)), 4], V7 = r[t[1:(nrow(t) - 
        1), 3], 2], V8 = r[t[2:(nrow(t)), 2], 3])
    otherchrom = r[t[2:(nrow(t)), 1], 1]
    rd <- rd[which(rd$V1 == otherchrom), ]
    rm(otherchrom)
    rd <- rd[which(complete.cases(rd)), ]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
