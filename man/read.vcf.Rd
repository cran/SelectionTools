\name{read.vcf}
\alias{read.vcf}
\title{Read a Beagle-Oriented VCF into SelectionTools}

\description{
Reads diploid genotype data from a VCF file directly into a SelectionTools
marker data set.  The function is intentionally focused on interchange with
Beagle rather than on preserving every VCF annotation.
}

\usage{
read.vcf(vcffile, data.set = "default")
}

\arguments{
  \item{vcffile}{
    Character string giving the input VCF file.  Files ending in \code{.gz}
    are decompressed with the documented base-R \code{gzfile()} connection to
    a temporary plain VCF before the native parser is called.
  }
  \item{data.set}{
    Character string naming the SelectionTools marker data set to create or
    replace.
  }
}

\details{
The native parser reads the VCF \code{GT} subfield and ignores other FORMAT,
INFO, QUAL, and FILTER information.  Diploid phased (\code{|}) and unphased
(\code{/}) GT fields are accepted.  The two allele indices are stored in their
input order; SelectionTools marker data do not maintain a separate VCF phase
flag.  If either allele of a diploid genotype is missing, the complete diploid
genotype is stored as missing.

VCF REF/ALT allele indices are normalized to the compact SelectionTools allele
codes used internally.  Marker and sample names are retained.  The VCF
chromosome labels are represented by consecutive SelectionTools chromosome
numbers in first-appearance order.

Files written by \code{write.vcf()} contain the metadata line
\code{SelectionToolsPositionScale=1000000}.  For ordinary Beagle-style VCF
input without this metadata, physical positions are interpreted with the
Beagle no-map convention of one centiMorgan per megabase.  Consequently this
reader is intended for Beagle interchange, not for reconstructing an arbitrary
external genetic map from physical coordinates.

No external compression library is linked by SelectionTools.  Gzip transport
is handled only by base R; plain VCF parsing is implemented in C.
}

\value{
No return value.  The requested SelectionTools marker data set is created or
replaced for its side effect.
}

\seealso{
\code{write.vcf}, \code{st.phase}, \code{st.switch.error}
}
