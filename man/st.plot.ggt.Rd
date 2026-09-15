% copyright (C) 2001 Frisch/Maurer
%
\name{st.plot.ggt}
\alias{st.plot.ggt}
\title{st.plot.ggt()}
\description{
  The function \emph{st.plot.ggt} plots graphical genotypes for all or a subset
  of individuals.
}
\usage{
st.plot.ggt(data.set = "default", ifilename = "", i.list = "", f.ind = 0,
            l.ind = 0, f.lr = 2, f.tb = 2, d.h = 0.001, d.v = 0.001, p.t = FALSE,
            p.s = FALSE, d.t = 50, d.map = 100, c.nme = "", i.nme = "",
            cex.chrom = 1, cex.ind = 1,
            color = c ("yellow","blue","red","green", "wheat",
                       "skyblue","tomato","palegreen",
                       "yellow4","darkblue",
                       "darkblue","darkgreen" ),
            z.min = 0, z.max = 12, plt.pdf = FALSE, plt.png = FALSE,
            plt.fname = NULL, plt.width = 10, plt.height = 10,
            plt.ptsize = 10)
}
\arguments{
   \item{data.set}{character. Uses the data set that was created by \emph{st.read.marker.data}.
                If no data set was specified, the "default" data set is used.}
   \item{ifilename}{Name of a file containing a list of individuals of which the graphical
                  genotypes should be plotted. The file is by default read from the input directory.}
   \item{i.list}{Function argument.}
   \item{f.ind}{Position of first individual used for plotting.}
   \item{l.ind}{Position of last individual used for plotting.}
   \item{f.lr}{Width between color bars and labels (individual names).}
   \item{f.tb}{Height between color bars and top and bottom labels (chromosome names and marker names).}
   \item{d.h}{Horizontal distance between chromosomes.}
   \item{d.v}{Vertical distance between individuals.}
   \item{p.t}{logical. If TRUE marker tickmarks are plotted.}
   \item{p.s}{logical. If TRUE, a cM scale is plotted.}
   \item{d.t}{Length of marker tickmarks.}
   \item{d.map}{Distance between marker labels.}
   \item{c.nme}{Character string specifying the names of chromosomes.}
   \item{i.nme}{Character string specifyin the names of individuals.}
   \item{cex.chrom}{Character expansion factor for chromosome names.}
   \item{cex.ind}{Character expansion factor for individual names.}
   \item{color}{Character vector specifying the colors used to represent genotype or allele codes.}
   \item{z.min}{Minimum genotype or allele code represented by the color scale.}
   \item{z.max}{Maximum genotype or allele code represented by the color scale.}
   \item{plt.pdf}{logical. If TRUE a pdf file is generated.}
   \item{plt.png}{logical. If TRUE, a PNG file is generated.}
   \item{plt.fname}{Character string giving the base name of the output file. Required when \code{plt.pdf = TRUE} or \code{plt.png = TRUE}.}
   \item{plt.width}{Width of generated pdf file.}
   \item{plt.height}{Height of generated pdf file.}
   \item{plt.ptsize}{Point size of generated pdf file.}
}
\value{
On a successful run, an invisible implementation-level \code{list} returned by the final timing/information helper; the meaningful result is the side effect of producing the graphical genotype plot. Validation failures may return \code{NULL}.
}

\details{The color coding for SNP data is by default A=yellow, C=blue, G=red, T=green. For
        haplo blocks the coding is dependent on the frequency of a particular haplo block.
        The most frequent haplo block is plotted in yellow, the second frequent in blue
        and so on.\cr
        Using the function \emph{st.plot.ggt} after the function \emph{st.recode.hbc} gives a good graphical
        overview, which haplo blocks are not inherited from the reference genotype (if reference
        individual = 1, the reference parent genome is plotted in yellow. The haplo
        blocks inherited from the donor parent are plotted in blue.\cr
        Using the function \emph{st.plot.ggt} after the function \emph{st.recode.hil} gives a good graphical
        overview, which haplo blocks are inherited from which parent. For plotting, the individuals are sorted
        by their haplotype frequency. Thus, the parents are usually plotted on top in yellow and blue (default values).
        The haplo blocks of the individuals are also plotted in yellow or blue, if they are completely identical
        to a parent's haplo block. If they are plotted in a different color, a cross-over did occur and thus,
        this haplo block is not identical one of the parent's haplo blocks.\cr
        If an external list of individuals is used, the ordering of individuals
        within that list is used for plotting. When plotting on the current
        graphics device, the function stops if that device is already in
        split-screen mode, so that an existing user-defined split-screen layout
        is not changed.}

%\note{
%
%}
%
