\name{st.calc.ld.2}
\alias{st.calc.ld.2}
\title{Calculate linkage disequilibrium from marker data}

\description{
Calculate linkage disequilibrium from marker data.
}

\usage{
st.calc.ld.2(ld.measure = "r2", auxfiles = TRUE, keep.ldfile = FALSE,
             data.set = "default", ld.filename = NULL)
}

\arguments{
   \item{ld.measure}{Function argument used by this operation.}
   \item{auxfiles}{Logical. Controls generation of the LD result file. With \code{keep.ldfile = FALSE}, this file is created in the R session temporary directory, read back, and removed before the function returns.}
   \item{keep.ldfile}{Logical. If \code{FALSE}, the LD result is transferred through a session-temporary file and returned to R. If \code{TRUE} and \code{auxfiles = TRUE}, a persistent LD file is written and retained; \code{ld.filename} must then be supplied. The function returns \code{NULL} invisibly when \code{keep.ldfile = TRUE}.}
   \item{data.set}{Name of the SelectionTools data set.}
   \item{ld.filename}{Character string naming the persistent LD output file. Required only when \code{keep.ldfile = TRUE} and \code{auxfiles = TRUE}.}
}
\value{
Invisibly returns a \code{data.frame} with columns \code{Chrom}, \code{Locus1}, \code{Locus2}, \code{Name1}, \code{Name2}, and \code{LD} when the generated LD file is read back. If \code{keep.ldfile = TRUE}, the R return value is \code{NULL}; when \code{auxfiles = TRUE}, the LD file named by \code{ld.filename} is retained.
}
