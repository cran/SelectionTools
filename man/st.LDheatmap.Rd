\name{st.LDheatmap}
\alias{st.LDheatmap}
\title{Plots a Linkage Disequilibrium Heat Map}
\description{
Plots a triangular heat map of pairwise linkage disequilibrium (LD)
values using only the standard R graphics system. Marker cells are equally
spaced in the LD matrix. If marker positions are supplied, their genetic or
physical spacing is displayed separately on a map line connected to the
matrix.
}
\usage{
st.LDheatmap(LD, map = NULL,
             distances = "genetical",
             LDmeasure = "r",
             title = "",
             col = grDevices::heat.colors(20),
             zlim = c(0, 1),
             add.map = TRUE,
             add.key = TRUE,
             depth = Inf,
             max.dist = Inf,
             write.ld = NULL,
             write.snp.id = FALSE,
             text = FALSE,
             digits = 2,
             geneMapLocation = 0.15,
             geneMapLabelX = NULL,
             geneMapLabelY = NULL,
             cex.title = 0.9,
             cex.map = 0.7,
             cex.key = 0.6,
             cex.ld = NULL,
             cex.snp = NULL)
}
\arguments{
  \item{LD}{
    A numeric square matrix containing pairwise LD values. The marker order
    is given by the row and column order of the matrix. Values from the upper
    triangle are plotted. If an upper-triangle entry is \code{NA} and the
    corresponding lower-triangle entry is available, the lower-triangle
    value is used. The diagonal is not plotted.
  }
  \item{map}{
    A numeric vector containing one map position for each marker, in the same
    order as the rows and columns of \code{LD}. Positions must be finite and
    in non-decreasing order; equal positions are allowed. \code{map} is
    required when \code{add.map = TRUE} or when \code{max.dist} is finite.
  }
  \item{distances}{
    Character string specifying the units represented by \code{map}.
    \code{"genetic"} and \code{"genetical"} specify genetic positions in
    centiMorgans (cM). \code{"physical"} specifies physical positions in
    base pairs. For physical positions, the total map length is reported in
    kilobases. The spelling \code{"genetical"} is accepted as an alias for
    \code{"genetic"}.
  }
  \item{LDmeasure}{
    Character string used to label the color key. Values \code{"r"},
    \code{"r2"}, and \code{"r^2"} produce an \eqn{R^2} label. Values
    \code{"Dp"}, \code{"Dprime"}, \code{"D'"}, and \code{"D"} produce
    a \eqn{D'} label. Other strings are used literally. This argument does
    not calculate or transform LD values.
  }
  \item{title}{
    A single character string giving the plot title. An empty string omits
    the title.
  }
  \item{col}{
    A vector of at least two valid R colors. The first supplied color is used
    for the highest LD values and the last supplied color for the lowest LD
    values. Thus \code{grDevices::heat.colors(20)} gives high LD in red and
    low LD in pale yellow.
  }
  \item{zlim}{
    Numeric vector of length two giving the lower and upper limits of the LD
    color scale. The default is \code{c(0, 1)}. Values outside this interval
    generate an error, apart from tiny floating-point deviations at the
    limits, which are clipped to the limits.
  }
  \item{add.map}{
    Logical. \code{TRUE} draws a line representing the supplied marker
    positions and connects each marker in the equally spaced LD matrix to its
    position on that line.
  }
  \item{add.key}{
    Logical. \code{TRUE} draws the color key.
  }
  \item{depth}{
    A non-negative integer or \code{Inf}. Only marker pairs separated by at
    most \code{depth} marker intervals are displayed. For example,
    \code{depth = 1} displays only adjacent marker pairs. The default
    \code{Inf} applies no restriction based on marker index.
  }
  \item{max.dist}{
    A non-negative number or \code{Inf}. Only marker pairs separated by at
    most \code{max.dist} map units are displayed. The units are the units of
    \code{map}: cM for a genetic map and base pairs for a physical map. A
    finite value requires \code{map}. The default \code{Inf} applies no
    restriction based on map distance.
  }
  \item{write.ld}{
    \code{NULL} or a function used to construct text labels for displayed LD
    values. The function is called separately for each displayed value and
    must return one printable value. Returning \code{""} or \code{NA}
    suppresses the label for that cell. For example,
    \code{function(x) if (x >= 0.8) sprintf("\%.2f", x) else ""} prints only
    LD values of at least 0.8. If supplied, \code{write.ld} takes precedence
    over \code{text} and \code{digits}.
  }
  \item{write.snp.id}{
    Logical. \code{TRUE} prints marker identifiers. Row names of
    \code{LD} are used if available, otherwise column names, otherwise the
    marker numbers \code{1, ..., n}. For densely spaced markers, labels can
    overlap; \code{cex.snp} can be supplied to control their size.
  }
  \item{text}{
    Logical convenience option. If \code{TRUE} and \code{write.ld} is
    \code{NULL}, all displayed LD values are printed using a fixed number of
    decimal places specified by \code{digits}.
  }
  \item{digits}{
    Non-negative integer giving the number of decimal places used when
    \code{text = TRUE} and \code{write.ld = NULL}.
  }
  \item{geneMapLocation}{
    Non-negative numeric value controlling the perpendicular distance of the
    map line from the matrix diagonal. Larger values move the map line farther
    from the matrix.
  }
  \item{geneMapLabelX}{
    Optional finite numeric x-coordinate for the text reporting the total map
    length. If \code{NULL}, a value of 0.52 is used.
  }
  \item{geneMapLabelY}{
    Optional finite numeric y-coordinate for the text reporting the total map
    length. If \code{NULL}, a value of 0.24 is used.
  }
  \item{cex.title}{
    Positive finite numeric character expansion factor for the plot title.
  }
  \item{cex.map}{
    Positive finite numeric character expansion factor for the map-length
    label.
  }
  \item{cex.key}{
    Positive finite numeric character expansion factor for the color-key
    title and tick labels.
  }
  \item{cex.ld}{
    \code{NULL} or a positive finite numeric character expansion factor for
    LD values printed in the cells. If \code{NULL}, a size is estimated from
    the cell size and the labels to be printed.
  }
  \item{cex.snp}{
    \code{NULL} or a positive finite numeric character expansion factor for
    marker identifiers. If \code{NULL}, a size is estimated from the number
    of markers and the identifier lengths.
  }
}
\details{
\code{st.LDheatmap} is a plotting function; it does not calculate LD. A
precomputed pairwise LD matrix is supplied in \code{LD}. In a typical
SelectionTools workflow the matrix and map can be obtained with
\code{st.LDplot.ld} and \code{st.LDplot.map}, respectively.

The heat-map cells represent marker order rather than map distance. The marker
centers are therefore equally spaced from zero to one along both matrix axes.
When \code{add.map = TRUE}, the actual values in \code{map} are rescaled to a
separate line parallel to the matrix diagonal, and each matrix position is
connected to its corresponding map position. This keeps the LD cells equally
sized while still displaying uneven marker spacing. Markers with identical map
positions are supported and are connected to the same position on the map
line.

Only the upper triangle of \code{LD} is displayed, and the diagonal is omitted.
An upper-triangle value is authoritative when it is present. If it is missing
and the corresponding lower-triangle value is non-missing, that lower-triangle
value is copied to the displayed upper triangle. Missing LD values remain blank.
The lower triangle is otherwise ignored.

The restrictions specified by \code{depth} and \code{max.dist} are applied
simultaneously. A cell is displayed only if it satisfies both restrictions.
The restrictions affect only which LD cells are drawn; they do not change the
marker order or the map line.

The color scale covers \code{zlim}. The order of \code{col} is from high LD to
low LD. This convention makes \code{heat.colors()} produce the customary red
for high LD and yellow for low LD. Missing or filtered cells are not painted.

If \code{write.ld} is supplied, it is evaluated once for each displayed LD
value. This permits selective labels without changing the plotted matrix. If
\code{cex.ld = NULL}, the function estimates a label size intended to fit into
one heat-map cell. Automatic sizing is necessarily approximate on very dense
plots; an explicit \code{cex.ld} can be supplied when exact control is needed.
The same consideration applies to marker identifiers and \code{cex.snp}.

The function uses standard R graphics and starts one new plot in the current
graphics layout. It does not modify global \code{par} settings. Consequently,
standard layouts such as \code{par(mfrow = c(2, 2))} can be used to place
several LD heat maps on one graphics device.
}
\value{
An invisible named \code{list} with components:

  \item{LDmatrix}{The upper-triangular LD matrix actually displayed after the plotting restrictions are applied.}
  \item{map}{The supplied map vector, or \code{NULL}.}
  \item{snp.position}{Plotting coordinates of markers in the LD matrix.}
  \item{map.x, map.y}{Plotting coordinates for the map line, or \code{NULL} when no map is added.}
  \item{distances}{The map-distance type used internally.}
  \item{depth}{The marker-depth restriction.}
  \item{max.dist}{The map-distance restriction.}
  \item{zlim}{The color-scale limits.}
  \item{col}{The color vector used for the heat map.}

}
\examples{
data("v-tropmaize-vcf", package = "SelectionTools")
marker.frame <- st.STvcf.to.dataframe(v.tropmaize.vcf)[1:12, ]
st.load.vcf.data(st.dataframe.to.STvcf(marker.frame), data.set = "maize.ld")
Z <- gs.get.Z("maize.ld", names = TRUE)
LD <- cor(Z, use = "pairwise.complete.obs")^2
map <- st.get.map(data.set = "maize.ld")
st.LDheatmap(LD, map$Pos, title = "Tropical maize LD",
             col = grDevices::heat.colors(20))
st.LDheatmap(LD, map$Pos, depth = 2, write.snp.id = TRUE)
}
\seealso{
\code{st.calc.ld}, \code{st.LDplot.ld}, \code{st.LDplot.map}
}
\keyword{hplot}
