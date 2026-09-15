% copyright (C) 2001 Frisch/Maurer
%
\name{st.def.hblocks}
\alias{st.def.hblocks}
\title{st.def.hblocks()}
\description{
Building of haplo blocks either based on number of markers, based on a given
genetic distance in cM, or for dissecting an introgression library in
disjunct segments.
}
\usage{
st.def.hblocks(hap = 5, hap.unit = 1, ld.threshold = 0.8,
               ld.criterion = "none", tolerance = 0, hap.symbol = "b",
               out.filename = "blocks", auxfiles = TRUE, data.set = "default")
}
\arguments{
   \item{hap}{Defines either a number of markers (if \emph{hap.unit=1}) or a genetic distance
            in cM (if \emph{hap.unit=2}) that should be used for building haplo blocks}
   \item{hap.unit}{Set to 1 for building haplo blocks based on number of markers, set to 2 for
                  building haplo blocks based on markers lying in a
		  specified stretch of the chromosome (in cM). Set to 0
		  to build disjunct chromosome segments in an
		  introgression population}
   \item{ld.threshold}{Function argument.}
   \item{ld.criterion}{Function argument.}
   \item{tolerance}{Function argument.}
   \item{hap.symbol}{Function argument.}
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{character. Name of the file that is used for building haplo blocks}
}

%\details{}

\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{Pos}, \code{Name}, \code{Class}, and \code{Markers} describing the defined haplotype blocks when \code{auxfiles = TRUE} and the operation succeeds; otherwise \code{NULL}.
}

%\note{
%
%}
%
