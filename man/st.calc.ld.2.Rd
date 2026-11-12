\name{st.calc.ld.2}
\alias{st.calc.ld.2}
\title{Calculates linkage disequilibrium from marker data}

\description{
Calculates linkage disequilibrium from marker data.
}

\usage{
st.calc.ld.2(ld.measure = "r2", auxfiles = TRUE, keep.ldfile = FALSE,
             data.set = "default", ld.filename = NULL)
}

\arguments{
   \item{ld.measure}{Function argument used by this operation.}
   \item{auxfiles}{Logical. Controls availability of the LD result. With \code{keep.ldfile = FALSE}, the result table is returned directly from memory.}
   \item{keep.ldfile}{Logical. If \code{FALSE}, the LD result is returned to R directly from memory. If \code{TRUE} and \code{auxfiles = TRUE}, a persistent LD file is written and retained; \code{ld.filename} must then be supplied. The function returns \code{NULL} invisibly when \code{keep.ldfile = TRUE}.}
   \item{data.set}{Name of the SelectionTools data set.}
   \item{ld.filename}{Character string naming the persistent LD output file. Required only when \code{keep.ldfile = TRUE} and \code{auxfiles = TRUE}.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{Locus1}, \code{Locus2}, \code{Name1}, \code{Name2}, and \code{LD} when the LD result is requested in memory. If \code{keep.ldfile = TRUE}, the R return value is \code{NULL}; when \code{auxfiles = TRUE}, the LD file named by \code{ld.filename} is retained.
}
