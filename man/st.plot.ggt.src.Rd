\name{st.plot.ggt.src}
\alias{st.plot.ggt.src}
\title{Plot Graphical Genotypes from Map and Genotype Objects}
\description{
Low-level plotting routine used by \code{st.plot.ggt}.  It draws a graphical
genotype representation from an already prepared linkage map and genotype
matrix on the currently active R graphics device.
}
\usage{
st.plot.ggt.src(map, pop, f.lr, f.tb, d.h, d.v, p.t, p.s, d.t, d.map,
                c.nme, i.nme, color, z.min, z.max, cex.chrom, cex.ind)
}
\arguments{
\item{map}{Data frame or matrix describing chromosome, position, marker name,
and marker class in the format used by \code{st.plot.ggt}.}
\item{pop}{Genotype matrix; rows are markers and columns are individuals.}
\item{f.lr}{Relative width of the left individual-label area.}
\item{f.tb}{Relative height of the top and bottom label areas.}
\item{d.h}{Horizontal distance between chromosome panels.}
\item{d.v}{Vertical distance between individual panels.}
\item{p.t}{Logical.  If \code{TRUE}, marker tick marks are plotted.}
\item{p.s}{Logical.  If \code{TRUE}, a map-position scale is plotted.}
\item{d.t}{Length of marker tick marks.}
\item{d.map}{Distance between map-position labels.}
\item{c.nme}{Character vector of chromosome labels.}
\item{i.nme}{Character vector of individual labels.}
\item{color}{Character vector of colors used for genotype or allele codes.}
\item{z.min}{Minimum genotype or allele code represented by the color scale.}
\item{z.max}{Maximum genotype or allele code represented by the color scale.}
\item{cex.chrom}{Character expansion factor for chromosome labels.}
\item{cex.ind}{Character expansion factor for individual labels.}
}
\details{
This is a low-level plotting function.  It does not open or close a graphics
device.  The caller controls file output with the standard R graphics-device
functions.  The function uses \code{split.screen()} internally, refuses to run
inside an already active split-screen layout, and restores graphics parameters
when it exits.  Missing genotype codes are represented as missing values and
are not painted.
}
\value{
No useful value is intended.  The function is called for its side effect of
drawing on the active graphics device.
}
