\name{write.vcf}
\alias{write.vcf}
\title{Write SelectionTools Marker Data for Beagle}

\description{
Writes a SelectionTools marker data set as a minimal GT-only VCF intended for
direct use as Beagle target genotype input.
}

\usage{
write.vcf(vcffile, data.set = "default")
}

\arguments{
  \item{vcffile}{
    Character string giving the output VCF file.  If the name ends in
    \code{.gz}, the native writer first creates a temporary plain VCF and base
    R \code{gzfile()} compresses it to the requested file.
  }
  \item{data.set}{
    Character string naming an existing SelectionTools marker data set.
  }
}

\details{
The output is deliberately small and Beagle-oriented. It contains VCF 4.3
header information and one FORMAT field, \code{GT}. A marker data set whose
stored phase state is known is written with the phased separator \code{|}; an
unphased data set is written with \code{/}. Missing diploid genotypes are
written consistently as \code{.|.} or \code{./.}, respectively. Thus VCF
round trips retain the SelectionTools data-set phase state.

SelectionTools allele codes are normalized independently at each marker to VCF
allele indices 0, 1, and so on.  REF and ALT are valid synthetic DNA allele
labels used only for interchange; they must not be interpreted as the original
biological nucleotide labels.  Marker IDs and individual names are retained.

SelectionTools genetic positions are in centiMorgans.  To make a VCF usable by
Beagle without a separate genetic-map file, positions are converted to
pseudo-base-pair positions with one centiMorgan equal to one megabase, which
is Beagle's documented no-map scale.  Tied map positions are made strictly
increasing within a chromosome.  The file records the scale in a
\code{SelectionToolsPositionScale} metadata line.

No zlib, htslib, gzip executable, or other additional system library is linked
by SelectionTools.  Compression of \code{.vcf.gz} output uses only the public
base-R \code{gzfile()} interface.
}

\value{
No return value.  The VCF file is written for its side effect.
}

\seealso{
\code{read.vcf}, \code{st.destroy.phase}, \code{st.phase}
}
