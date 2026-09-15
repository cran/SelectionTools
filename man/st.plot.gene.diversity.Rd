% copyright (C) 2001 Frisch/Maurer
%
\name{st.plot.gene.diversity}
\alias{st.plot.gene.diversity}
\title{st.plot.gene.diversity()}
\description{
The function \emph{st.plot.gene.diversity} plots gene diversity along the
chromosomes and returns the marker-level values used for the plot.
}
\usage{
st.plot.gene.diversity(data.set ="default", plt.pdf = FALSE,
                       plt.fname = NULL, plt.width = 10, plt.height = 6,
                       plt.ptsize = 14)
}
\arguments{
   \item{data.set}{character. Name of the data set that is used by the function.}
   \item{plt.pdf}{logical. If \code{TRUE}, a PDF file is generated in the output directory.}
   \item{plt.fname}{Character string giving the base name of the PDF file. Required when \code{plt.pdf = TRUE}.}
   \item{plt.width}{Width of the generated PDF file.}
   \item{plt.height}{Height of the generated PDF file.}
   \item{plt.ptsize}{Point size of the generated PDF file.}
}

\details{
The gene diversity (expected heterozygosity) at a marker is calculated as
\deqn{1-\sum_i p_i^2,}
where \eqn{p_i} is the frequency of allele \eqn{i} at that marker. The values
used for plotting are obtained from \code{st.marker.data.statistics()}. The
average gene diversity for each chromosome is printed at the top of the plot.
}

\value{
Invisibly returns a \code{data.frame} with columns \code{chrom}, \code{pos}, \code{name}, \code{cpos}, and \code{gene.div}, containing the marker-level positions and gene-diversity values used for plotting.
}

%\note{
%
%}
%
