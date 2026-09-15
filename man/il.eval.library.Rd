\name{il.eval.library}
\alias{il.eval.library}
\title{Construct or evaluate introgression-line libraries}

\description{
Construct or evaluate introgression-line libraries.
}

\usage{
il.eval.library(rp.allele = 1, dp.allele = 2, chrom = 0, lower = 0, upper = 0,
                data.set = "default")
}

\arguments{
   \item{rp.allele}{Allele-related parameter.}
   \item{dp.allele}{Allele-related parameter.}
   \item{chrom}{Chromosome identifier.}
   \item{lower}{Function argument used by this operation.}
   \item{upper}{Function argument used by this operation.}
   \item{data.set}{Name of the SelectionTools data set.}
}
\value{
A named numeric vector summarizing the evaluated introgression-line library. It contains \code{cov}, \code{dep}, \code{seglib}, \code{dpg}, \code{ndps}, and \code{ldps}; when target regions are supplied it additionally contains \code{tr.dpg}, \code{nt.rpg}, \code{nt.ndps}, and \code{nt.ldps}. Returns \code{NULL} on validation failure.
}
