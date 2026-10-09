\name{masbc.input}
\alias{masbc.input}
\docType{data}
\title{Input Files for Marker-Assisted Backcrossing Examples}
\description{
A collection of linkage maps and parental marker-genotype files used by the
marker-assisted backcrossing examples.  The files are stored as an R data
object so that they can be loaded with \code{data()} and written to a
temporary directory when the file-based \code{mab.*} interface is used.
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

The marker-assisted backcrossing functions still use file names.  Therefore,
a vignette or analysis can load \code{masbc.input}, write the required list
elements to \code{tempdir()}, and set \code{st.input.dir} to that temporary
directory.
}
\examples{
data("masbc.input")

mab.input.dir <- file.path(tempdir(), "masbc-input")
dir.create(mab.input.dir, showWarnings = FALSE)

invisible(Map(
    function(x, file) writeLines(x, file.path(mab.input.dir, file)),
    masbc.input,
    names(masbc.input)
))

list.files(mab.input.dir)
unlink(mab.input.dir, recursive = TRUE)
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
