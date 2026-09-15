\name{st.plot.ggt.src}
\alias{st.plot.ggt.src}
\title{Create a graphical genotype plot from SelectionTools source data}

\description{
Create a graphical genotype plot from SelectionTools source data.
}

\usage{
st.plot.ggt.src(map, pop, f.lr, f.tb, d.h, d.v, p.t, p.s, d.t, d.map, c.nme,
                i.nme, color, z.min, z.max, cex.chrom, cex.ind, plt.png,
                plt.pdf, plt.name, plt.width, plt.height, plt.ptsize)
}

\arguments{
   \item{map}{Linkage map or map object.}
   \item{pop}{Population name or population specification.}
   \item{f.lr}{Function argument used by this operation.}
   \item{f.tb}{Function argument used by this operation.}
   \item{d.h}{Function argument used by this operation.}
   \item{d.v}{Function argument used by this operation.}
   \item{p.t}{Function argument used by this operation.}
   \item{p.s}{Function argument used by this operation.}
   \item{d.t}{Function argument used by this operation.}
   \item{d.map}{Map-related parameter.}
   \item{c.nme}{Function argument used by this operation.}
   \item{i.nme}{Function argument used by this operation.}
   \item{color}{Function argument used by this operation.}
   \item{z.min}{Function argument used by this operation.}
   \item{z.max}{Function argument used by this operation.}
   \item{cex.chrom}{Chromosome-related parameter.}
   \item{cex.ind}{Function argument used by this operation.}
   \item{plt.png}{Function argument used by this operation.}
   \item{plt.pdf}{Function argument used by this operation.}
   \item{plt.name}{Function argument used by this operation.}
   \item{plt.width}{Function argument used by this operation.}
   \item{plt.height}{Function argument used by this operation.}
   \item{plt.ptsize}{Function argument used by this operation.}
}
\value{
No useful return value is intended. The function is called for its side effect of producing a graphical genotype plot and, when requested, closing the graphics device.
}
