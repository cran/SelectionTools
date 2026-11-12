\name{read.vcf}
\alias{read.vcf}
\title{Reads a Beagle-Oriented VCF into SelectionTools}

\description{
Reads diploid genotype data from a VCF file directly into a SelectionTools
marker data set.  The function is intentionally focused on interchange with
Beagle rather than on preserving every VCF annotation.
}

\usage{
read.vcf(vcffile, data.set = "default", pedigree = TRUE,
         strictACTG = FALSE)
}

\arguments{
\item{pedigree}{Logical; controls writing or reading quoted \code{SelectionToolsPedigree} records in the VCF header.}
\item{strictACTG}{Logical, default \code{FALSE}. When enabled, the reader accepts only
A, C, G and T allele labels or their abbreviations Ade, Cyt, Gua and Thy,
and uses fixed internal codes A/Ade=1, C/Cyt=2, G/Gua=3 and T/Thy=4,
independent of REF/ALT order. Matching is case-insensitive. Other allele
labels are treated as missing.}

  \item{vcffile}{
    Character string giving the input VCF file.  Files ending in \code{.gz}
    are decompressed with the documented base-R \code{gzfile()} connection to
    a temporary plain VCF before the compiled parser is called.
  }
  \item{data.set}{
    Character string naming the SelectionTools marker data set to create or
    replace.
  }
}

\details{
Full simulated pedigrees are stored in quoted \code{SelectionToolsPedigree}
header records within the VCF, including compressed files. Sample names retain
stable individual identifiers and may contain long pedigree-based names.
Without these records, imported individuals are founders whose pedigree equals
their sample name. Incomplete, duplicate or inconsistent pedigree records cause
import to fail before replacing the destination. \code{pedigree=FALSE} omits
these records on output and ignores them on input. Fitted prediction state is
retained by the in-memory STvcf representation, rather than plain VCF.

The compiled parser reads the VCF \code{GT} subfield and ignores other FORMAT,
INFO, QUAL, and FILTER information. Diploid phased (\code{|}) and unphased
(\code{/}) GT fields are accepted. The two allele indices are stored in their
input order. A marker data set is recorded as phased only when genotype
separators are consistently \code{|}; mixed \code{|} and \code{/} input is
accepted but recorded as unphased. If either allele of a diploid genotype is
missing, the complete diploid genotype is stored as missing.

VCF REF/ALT allele indices are normalized to the compact SelectionTools allele
codes used internally.  Marker and sample names are retained.  The VCF
chromosome labels are represented by consecutive SelectionTools chromosome
numbers in first-appearance order.

With \code{strictACTG=TRUE}, labels in REF and ALT determine the fixed
nucleotide codes instead. The option retains a marker if at least one of
its REF/ALT labels is A, C, G, T, Ade, Cyt, Gua or Thy (case-insensitive).
It discards markers with no such label; an accepted label need not be
observed in a nonmissing genotype.
Indels, symbolic alleles, spanning-deletion stars, ambiguous bases and
numeric allele labels are not accepted as nucleotides. If either allele
of a diploid call is excluded, the entire call is missing. Genotype order,
retained marker identifiers, map positions and individual pedigrees are
preserved; chromosome codes are consecutive after filtering. If no marker
remains, import fails without replacing an existing destination.
The default \code{FALSE} preserves the original import behavior.

Fixed nucleotide coding aligns matching marker identifiers with different
REF/ALT orders. It does not match marker identifiers or positions, complement
opposite strands, or normalize indel representations.

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
