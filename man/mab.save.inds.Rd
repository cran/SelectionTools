\name{mab.save.inds}
\alias{mab.save.inds}
\title{Helper used by the marker-assisted backcross simulation routines}

\description{
Writes one marker-data file for each individual in a population. Output file
names use the explicitly supplied \code{file.prefix}.
}

\usage{
mab.save.inds(pop, pname, file.prefix)
}

\arguments{
   \item{pop}{Population name or population specification.}
   \item{pname}{Prefix used for the individual names written inside the marker-data files.}
   \item{file.prefix}{Required prefix for the output file names. Each file receives a three-digit individual number and the extension \code{.mda}. Relative prefixes are resolved against \code{st.output.dir}; absolute paths are used directly.}
}
\value{
\code{NULL}. The function is called for its side effect of writing one marker-data file for each individual in the supplied population.
}
