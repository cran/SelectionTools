\name{masbc.input}
\alias{masbc.input}
\docType{data}
\title{Input Files for Marker-Assisted Backcrossing Examples}
\description{
A collection of linkage maps and parental marker-genotype files used by the
marker-assisted backcrossing examples. The files are stored as character
vectors in one R data object for the file-based \code{mab.*} interface.
}
\usage{
data(masbc.input)
}
\format{
A named list of 16 character vectors.  Each list element contains the exact
text lines of one input file.  The element names are the file names:
\code{ex-mab1.map}, \code{ex-mab2.map}, \code{ex-mab3.map},
\code{ex-mab4.map}, \code{ex-mabi1.mda} through \code{ex-mabi8.mda},
\code{ex-mabp1.mda}, \code{ex-mabp1_m.mda}, \code{ex-mabp2.mda}, and
\code{ex-mabp2_m.mda}.
}
\details{
The \code{*.map} elements contain linkage-map input files.  The \code{*.mda}
elements contain marker genotypes for donor, recurrent-parent, or individual
backcross plants.  Keeping the files as character vectors preserves the
original input-file representation while making the complete example data
available through the standard R \code{data()} mechanism.

The file-based marker-assisted backcrossing interface identifies these
input files by their list-element names.
}
\source{
Input files used by the marker-assisted backcrossing examples supplied with
SelectionTools.
}
\seealso{
\code{\link{mab.simulate}}, \code{\link{mab.compare}},
\code{\link{mab.load.data}}, \code{\link{mab.tabulate}},
\code{\link{mab.Input.files}}
}
\keyword{datasets}
