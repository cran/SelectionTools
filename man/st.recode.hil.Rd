% copyright (C) 2001 Frisch/Maurer
%
\name{st.recode.hil}
\alias{st.recode.hil}
\title{Recode Defined Haplotype Blocks as Marker Loci}
\description{
Recodes marker sequences within previously defined haplotype blocks and uses
the resulting haplotype variants as alleles of new haplotype-block loci.
}
\usage{
st.recode.hil(out.filename = "blocks", auxfiles = TRUE, data.set = "default")
}
\arguments{
   \item{out.filename}{Character string retained for compatibility and used as part of the name of an internal file in the R session temporary directory. The temporary file is removed before the function returns.}
   \item{auxfiles}{Logical. If \code{TRUE}, the compiled result is written to the R session temporary directory, read back by the R wrapper, and removed before the function returns.}
   \item{data.set}{Name of data set that is used for recoding.}
}

\details{
Haplotyping is a two-step procedure. First, \code{st.def.hblocks} defines the
marker boundaries of the haplotype blocks. Then \code{st.recode.hil} compares
the stored marker sequence on each homologue within each block and identifies
identical haplotype sequences.

Within a block, identical marker sequences receive the same haplotype allele
code. Haplotype alleles are ordered by observed frequency, with the most
frequent haplotype receiving the first code. A haplotype containing a missing
marker allele is represented as missing after recoding. If a block would
produce more than 254 distinct haplotype alleles, recoding stops with an error
rather than merging or truncating haplotypes.

The operation modifies \code{data.set}. After successful recoding, the
original marker-level representation is replaced by a new marker data set in
which the defined haplotype blocks are the loci and the observed haplotype
variants are their alleles. Make a copy of a marker data set before
haplotyping if the original marker-level representation is also required.

The returned haplotype description identifies each block and haplotype allele
and gives the sequence of original marker allele codes defining that
haplotype. The external haplotype description uses \code{-1} for missing
marker alleles.

Using \code{st.plot.ggt} after recoding provides a graphical overview of
haplotype-block sharing among individuals.
}

\value{
Invisibly returns a \code{data.frame} with columns \code{Block}, \code{AlleleNr}, and \code{AlleleDef} describing recoded haplotype-block alleles when \code{auxfiles = TRUE} and the operation succeeds; otherwise \code{NULL}.
}

\seealso{
\code{st.def.hblocks}, \code{st.calc.ld}, \code{st.plot.ggt},
\code{st.plot.gene.diversity}
}
