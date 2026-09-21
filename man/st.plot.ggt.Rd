% copyright (C) 2001 Frisch/Maurer
%
\name{st.plot.ggt}
\alias{st.plot.ggt}
\title{Plot Graphical Genotypes}
\description{
Plots graphical genotypes for all individuals or for a selected subset.  The
individuals can optionally be ordered by performance, and stored performance
values can be appended automatically to the individual labels.
}
\usage{
st.plot.ggt(data.set = "default", ifilename = "", i.list = "", f.ind = 0,
            l.ind = 0, f.lr = 2, f.tb = 2, d.h = 0.001, d.v = 0.001,
            p.t = FALSE, p.s = FALSE, d.t = 50, d.map = 100,
            c.nme = "", i.nme = "", cex.chrom = 1, cex.ind = 1,
            color = c("yellow", "blue", "red", "green", "wheat",
                      "skyblue", "tomato", "palegreen", "yellow4",
                      "darkblue", "darkblue", "darkgreen"),
            z.min = 0, z.max = 12,
            sort.performance = FALSE, show.performance = FALSE,
            performance.decreasing = TRUE, performance.format = "\%4.2f")
}
\arguments{
\item{data.set}{Character string giving the SelectionTools marker data set.}
\item{ifilename}{Character string giving a file with one individual identifier
per row.  The file is read from the input directory.  Ignored when
\code{i.list} is supplied.}
\item{i.list}{Vector of individual identifiers to plot.  The order is retained
unless \code{sort.performance = TRUE}.}
\item{f.ind}{Position of the first individual to plot.  Without performance
sorting this has the established meaning in the current data-set/list order.
With \code{sort.performance = TRUE}, it refers to the order after sorting.}
\item{l.ind}{Position of the last individual to plot.  Without performance
sorting this has the established meaning in the current data-set/list order.
With \code{sort.performance = TRUE}, it refers to the order after sorting.}
\item{f.lr}{Relative width of the left label area.}
\item{f.tb}{Relative height of the top and bottom label areas.}
\item{d.h}{Horizontal distance between chromosome panels.}
\item{d.v}{Vertical distance between individual panels.}
\item{p.t}{Logical.  If \code{TRUE}, marker tick marks are plotted.}
\item{p.s}{Logical.  If \code{TRUE}, a map-position scale is plotted.}
\item{d.t}{Length of marker tick marks.}
\item{d.map}{Distance between map-position labels.}
\item{c.nme}{Character vector with chromosome labels.  The default uses
consecutive chromosome labels.}
\item{i.nme}{Character vector with individual labels.  When supplied, it takes
precedence over automatic performance labels.  If performance sorting is
requested, these labels are reordered together with their individuals.}
\item{cex.chrom}{Character expansion factor for chromosome labels.}
\item{cex.ind}{Character expansion factor for individual labels.}
\item{color}{Character vector of colors used for genotype or allele codes.}
\item{z.min}{Minimum genotype or allele code represented by the color scale.}
\item{z.max}{Maximum genotype or allele code represented by the color scale.}
\item{sort.performance}{Logical.  If \code{TRUE}, order the eligible individuals
by the performance values stored in \code{data.set} before applying
\code{f.ind} and \code{l.ind}.}
\item{show.performance}{Logical.  If \code{TRUE} and \code{i.nme} is not
supplied, append the stored performance to each individual identifier.}
\item{performance.decreasing}{Logical.  Direction of performance sorting.
\code{TRUE} places the largest performance first.}
\item{performance.format}{A single \code{sprintf}-style format string used for
the numeric performance value in automatically generated labels.}
}
\details{
The plot is drawn on the currently active R graphics device.  File output is
therefore controlled by the caller, for example with \code{pdf()}, \code{png()},
or another graphics device followed by \code{dev.off()}.

If neither \code{sort.performance} nor \code{show.performance} is requested,
the established individual-selection and ordering behavior is unchanged.
When \code{sort.performance = TRUE}, sorting is restricted to the individuals
selected by \code{i.list} or \code{ifilename}; if neither is supplied, all
individuals in the data set are eligible.  The range \code{f.ind:l.ind} is then
applied after sorting.  With \code{sort.performance = FALSE}, \code{f.ind} and
\code{l.ind} retain their established meaning in the current data-set/list
order.

Performance values are obtained from the performance data stored in the
specified marker data set.  If performance sorting or automatic performance
labels are requested but performance data are unavailable, the function stops
with an error.  A supplied \code{i.nme} always takes precedence over automatic
performance labels.

For SNP data the default colors correspond to successive allele codes.  After
haplotype recoding, colors represent the resulting haplotype-allele codes.
Missing genotypes are left unpainted.  The function cannot be entered while the
current device is already in \code{split.screen()} mode, because it uses
split-screen panels internally.  Graphics parameters changed by the function
are restored on exit.
}
\value{
No useful value is intended.  The function is called for its side effect of
drawing the graphical genotype plot on the active graphics device.
}
