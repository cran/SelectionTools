\name{mab.simulate}
\alias{mab.simulate}

\title{Simulates a marker-assisted backcrossing program}

\description{
  Simulates one marker-assisted backcrossing program. The simulation result
  is stored in the explicitly supplied \code{result.file}. A relative file
  name is resolved against \code{st.data.dir}; an absolute path is used as
  supplied.
}

\usage{
mab.simulate(simulation.name, simulation.run, repetitions = 1000,
             recurrent.parent = "H", donor.parent = "H", linkage.map,
             target.loci = NULL, flanking.loci = NULL, recipient.loci = NULL,
             gen.type = NULL, population.size = NULL, sel.strategy = NULL,
             no.selected = NULL, no.preselected = NULL, reg.chr = NULL,
             reg.begin = NULL, reg.end = NULL, missing.allele = "9",
             result.file, success.factor = 1, recode.infiles = FALSE,
             recode.file = NULL)
}

\arguments{
   \item{simulation.name}{Name of the simulation stored in the returned result object. The argument is required and has no default value.}
   \item{simulation.run}{Name of the simulation run. Usually a single character. This can be
    used to run a certain simulation for several times. For example, in a
    test run with 1000 repetitions, the character \code{x} is used. And in a
    final run with 50000 repetitions the character \code{f} is used.  The
    argument is required and has no default value.}
   \item{repetitions}{Number of repetitions of the simulation. The default is
    1000. For final analyses, a larger number of repetitions may be used.}
   \item{recurrent.parent}{Description of the recurrent parent. If \code{H} is given, then a
    homozygous parent is assumed. Polymorphism between donor and recurrent
    parent is assumed for all markers on the linkage map.  If the name of a
    file is given, the genotype described by this file is used.  The file
    must be located in the directory \code{st.input.dir}.  Examples are
    given in the section 'Input Files'.  The argument is optional. Default
    value is \code{H}.  \emph{Note:} \code{recurrent.parent} and
    \code{donor.parent} need to be specified consistently. Either for both a
    \code{H}, or for both a file name.}
   \item{donor.parent}{Description of the donor parent. Further details appear in \code{recurrent.parent}.}
   \item{linkage.map}{Name of the file, in which the linkage map is described. Further details appear in the
    section 'Input Files' for a description of the file structure.The file
    must be located in the directory \code{st.input.dir}.  The argument is
    required and has no default value.}
   \item{target.loci}{Character vector with names of target loci. Selection for the donor
    allele at \code{target.loci} is indicated by a \code{t} in the
    definition of a selection strategy. The default is \code{NULL}.}
   \item{flanking.loci}{Character vector with names of loci, which are used to carry out
    pre selection for the recurrent parent allele before genome wide
    background selection is carried out. Selection for the recipient allele
    at \code{flanking.loci} is indicated by an \code{f} in the definition of
    a selection strategy. \code{flanking.loci} are used to reduce the length
    of the chromosome segment attached to the target genes. In general, they
    are flanking the target loci, but this in no requirement. For example,
    loci at a defined map distance from the target loci can be used. The
    default is \code{NULL}.}
   \item{recipient.loci}{Character vector, which contains the names of loci wich are used to carry
    out preselection for the recipient alleles. Selection for the recipient
    allele at \code{recipient.loci} is indicated by \code{r} in the
    definition of a selection strategy. \code{recipient.loci} are used to
    define chromosome regions of the recipient, which must be present in
    the backcrossing product. The default is \code{NULL}.}
   \item{gen.type}{A character vector consisting of a description of the generation
    types of the backcrossing program. The following generation types are
    implemented:

    \tabular{cl}{

      \code{f1} \tab
      F1 population
      \cr
      
      \code{bc} \tab
      Backcross population
      \cr

      \code{s} \tab
      Selfing population
      \cr
    
      \code{dh} \tab
      Doubled haploids
    }}
   \item{population.size}{Numeric vector containing the population sizes of the
    generations in the backcrossing program. The formal default is
    \code{NULL}.}
   \item{sel.strategy}{Character vector consisting of the selection strategies for the
  different generations of a backcrossing program.  The following
  selection strategies are implemented.  SI: selection index.

  \tabular{lll}{

    Strategy \tab
    Step \tab
    Description \cr

    \code{n}\tab
    1  \tab
    Random selection of individuals \cr

    \code{t}\tab
    1 \tab
    SI: Donor alleles at target loci \cr

    \code{tb}\tab
    1 \tab
    SI: Donor alleles at target loci \cr
    \tab
    2 \tab
    SI: Recurrent parent alleles genome wide \cr
   
    \code{tfb}\tab
    1 \tab
    SI: Donor alleles at target genes \cr
    \tab
    2 \tab
    SI: Recurrent parent alleles at flanking loci \cr
    \tab
    3 \tab
    SI: Recurrent parent alleles genome wide \cr

    \code{tf} \tab
    1 \tab
    SI: Donor alleles at target loci \cr
    \tab
    2 \tab
    SI: Recurrent parent alleles flanking loci \cr

    \code{tfs}\tab
    1 \tab
    SI: Donor alleles at target genes \cr
    \tab
    2 \tab
    SI: Recurrent parent alleles at flanking loci \cr
    \tab
    3 \tab
    SI: Low number of donor chromosome sgements \cr

    
    \code{tr} \tab
    1 \tab
    SI: Donor alleles at target loci \cr
    \tab
    2 \tab
    SI: Recurrent parent alleles at recipient loci \cr

    \code{trb} \tab
    1 \tab
    SI: Donor alleles at target loci \cr
    \tab
    2 \tab
    SI: Recurrent parent alleles at recipient loci \cr
    \tab
    3 \tab
    SI: Recurrent parent alleles genome wide \cr

    \code{trfb} \tab
    1 \tab
    SI: Donor alleles at target loci \cr
    \tab
    2 \tab
    SI: Recurrent parent alleles at recipient loci \cr
    \tab
     3 \tab
    SI: Recurrent parent alleles flanking loci \cr
    \tab
    4 \tab
    SI: Recurrent parent alleles genome wide \cr

    \code{tfrb} \tab
    1 \tab
    SI: Donor alleles at target loci \cr
    \tab
    2 \tab
    SI: Recurrent parent alleles at recipient loci \cr
    \tab
    3 \tab
    SI: Recurrent parent alleles genome wide \cr
    \tab
    4 \tab
    SI: Recurrent parent alleles flanking loci \cr
    
  }


  \emph{Selection strategies} \code{tb} and \code{tfb} are the
  two-stage and three-stage selection strategies described by Frisch et
  al. (1999).

  \emph{Counting marker data points:}

    \code{t}: 
    No markers are counted for the foreground selection for the target
    genes(s). (This selection step requires the same number of marker
    analyses for all selection stragegies.)

    \code{f}:
    Pre selection for markers flanking the target genes is assumed to be
    carries out with single marker assays. For each individual marker locus,
    a marker data point is counted. In advanced generations of a
    backcrossing program, only such markers are assayed, which were not
    already fixed for the recurrent parent allele in the non-recurrent
    parent (= the plant selected in the previous generation).

    \code{b}: 
    Background selection is assumed to be carried out with high
    throughput assays. For each plant one marker assay is counted for all
    loci on the linkage map (irrespective of their marker genotype).
    
    \code{s}:
    Selection for a low number of donor segments is assumed to be carried
    out with high throughput assays. For each plant one marker assay is
    counted for all loci on the linkage map. 

    \code{r}: 
    Pre selection for recipint alleles ia assumed to be carried out with
    high HT assays. For strategies \code{tr}, \code{trb}, \code{trfb}, one
    HT assay is counted for each plant carrying the target genes. For
    \code{tfrb}, the selection for the flanking markers is assumed to be
    carried out with single marker assays.}
   \item{no.selected}{Numerical vector specifying  the number of selected individuals in
    the backcross generations.  \emph{Note:} Specifying values larger than
    one makes only sense for some selection strategies.  The argument is
    optional, default value is a vector consisting of ones.}
   \item{no.preselected}{Numerical vector specifying  the minimum number of individuls that are
    preselected when selection for flanking marker (step "f") is carried out.
    \emph{Note:} Only used for selection strategies including
    preselection for flanking markers. The argument is
    optional. If not specified, all individuals are preselected that show recombination
    at the maximum number of flanking markers.}
   \item{reg.chr}{Definition of chromosome regions, which are evaluated separately for
    the recurrent parent genome.  \code{reg.chr} defines the chromosome
    number of the region to be evaluated. The default is \code{NULL}.}
   \item{reg.begin}{Definition of chromosome regions, which are evaluated separately for
    the recurrent parent genome.  \code{reg.begin} defines the begin of the
    region to be evaluated in cM distance from the telomere. The default is
    \code{NULL}. If \code{reg.chr} is specified, a corresponding value is
    required.}
   \item{reg.end}{End of the region to be evaluated, in cM distance from the
    telomere. If 0 is given, the entire chromosome is evaluated. The default
    is \code{NULL}. If \code{reg.chr} is specified, a corresponding value is
    required.}
   \item{missing.allele}{Character string that specifies missing marker data.
    The default is \code{"9"}.}
   \item{result.file}{Required name of the file in which the simulation result is stored. Relative names are resolved against \code{st.data.dir}; absolute paths are used directly. No file name is generated automatically.}
   \item{success.factor}{Numeric success probability used after selection.
    The default is 1. If it is below 1, the number of successful selected
    individuals is drawn from a binomial distribution and at least one
    individual is retained.}
   \item{recode.infiles}{Logical. If \code{TRUE}, recoded parental marker data are written to the explicitly supplied \code{recode.file}.}
   \item{recode.file}{File name for recoded parental marker data when \code{recode.infiles = TRUE}. It is required in that case. Relative names are resolved against \code{st.output.dir}; absolute paths are used directly.}
}

\value{
On successful completion, invisibly returns a named \code{list} containing
both the simulation definition and the simulation summaries. The parameter
part includes \code{simulation.name}, \code{simulation.run},
\code{repetitions}, \code{recurrent.parent}, \code{donor.parent},
\code{linkage.map}, \code{target.loci}, \code{flanking.loci},
\code{recipient.loci}, \code{gen.type}, \code{population.size},
\code{sel.strategy}, \code{no.selected}, \code{no.preselected},
\code{reg.chr}, \code{reg.begin}, \code{reg.end}, and
\code{success.factor}.

Depending on the requested analyses, the result part can contain:
\describe{
  \item{\code{TargetAlleles}}{Summary of the frequency of target alleles in the selected plants over
  the generations of the backcrossing program.}
  \item{\code{RPG.gw}}{Summary of recurrent-parent genome content of the
  selected plants over the whole genome.}
  \item{\code{RPG.reg1}, \code{RPG.reg2}, \dots}{Recurrent-parent genome
  content in chromosome regions defined by \code{reg.chr}, \code{reg.begin},
  and \code{reg.end}; the suffix identifies the region.}
  \item{\code{DonorSegments.gw}}{Number of donor chromosome segments in the
  selected individuals over the whole genome. Segments are counted per
  homologue, so the two homologues of each chromosome in a diploid individual
  are counted separately.}
  \item{\code{DonorSegments.reg1}, \code{DonorSegments.reg2}, \dots}{Number
  of donor chromosome segments within separately evaluated chromosome
  regions.}
  \item{\code{LinkageDrag1}, \code{LinkageDrag2}, \dots}{Length of the donor
  chromosome segment attached to each target locus. Target loci are evaluated
  separately and linkage drag is averaged over the two homologues.}
  \item{\code{SM}}{Summary of the number of marker data points attributed to
  single-marker assays. In the implemented accounting these arise from
  selection steps using flanking-marker assays.}
  \item{\code{HT}}{Summary of the number of marker data points attributed to
  high-throughput assays, including genome-wide background selection and
  selection for a low number of donor segments as specified by the selection
  strategy.}
}

Validation failures return \code{NULL} invisibly. The same list is saved in
the explicitly supplied \code{result.file}.
}

\seealso{
  \code{\link{mab.compare}}, \code{\link{mab.load.data}},
  \code{\link{mab.tabulate}}, \code{\link{mab.Input.files}},
  \code{\link{mab.Examples}}
}
