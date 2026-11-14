\name{mab.Input.files}
\alias{mab.Input.files}


\title{
Structure of input files for the optimize-mab routines
}

\description{
  Structure of the files describing linkage maps and crossing parents.
  Files \emph{must} be placed in the directory \code{st.input.dir}.
}

\details{
A linkage-map file contains chromosome number, position in cM, and marker
name. The old package documentation gave, for example:
\preformatted{
1       0      umc94
1       0002   bnl805
1       0067   umc76
1       0070   umc137
2       0      bnl845
2       0012   umc53
}

For a homozygous crossing parent, the first line is the genotype name and
each following line contains a marker name and one allele, for example:
\preformatted{
dh.2010.12345
s1e2855 ade
s1e3177 cyt
s1e2083 ade
s1e4552 gua
}

For a heterozygous crossing parent, the first line is again the genotype
name and each following line contains a marker name and two alleles, for
example:
\preformatted{
bc1i3
umc94   1  8
bnl805  1  8
umc76   8  8
umc137  8  8
}
}

\value{
No return value. This is a documentation page describing the input-file
formats used by the marker-assisted backcrossing routines.
}

\seealso{
  \code{\link{mab.simulate}}, \code{\link{mab.Examples}}
}
