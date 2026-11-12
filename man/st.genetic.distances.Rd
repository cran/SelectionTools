% copyright (C) 2001 Frisch/Maurer
%
\name{st.genetic.distances}
\alias{st.genetic.distances}
\alias{st.plot.genetic.distances}
\alias{st.plot.pcoa}
\title{Genetic Distances and Population Plots}
\description{
  The function \emph{st.genetic.distances} calculates genetic distances for each possible
  pair of individuals. The plotting helpers draw a genetic-distance heatmap
  and principal-coordinate score or variance plots using standard R graphics.
}
\usage{
st.genetic.distances(measure = "mrd", format = "l",
                     filename = NULL, auxfiles = FALSE,
                     data.set = "default")
st.plot.genetic.distances(x, labels = FALSE, cluster = TRUE,
                         col = grDevices::heat.colors(200),
                         hclust.method = "complete", ...)
st.plot.pcoa(x, which = c("scores", "variance", "cumulative"),
             axes = c(1L, 2L), groups = NULL, col = NULL,
             pch = NULL, legend = TRUE, cex = 0.8, ...)
}
\arguments{
   \item{measure}{character. Defines the measure to be used for calculation. Incorporated
                measures are modified Rogers distance ("mrd"), Rogers distance ("rd")
                and Euclidean distance ("euc").}
   \item{format}{character. Specifies the format of the output file. The value "l" specifies
               long format or to "m" for matrix format.}
   \item{filename}{Character string giving the base name of the genetic-distance file. Required when \code{auxfiles = TRUE}; the file is written with extension \file{.gdi} in the output directory.}
   \item{auxfiles}{Logical. If \code{TRUE}, a persistent output file with the genetic distances is generated and \code{filename} must be supplied.}
   \item{data.set}{character. Name of the data set that is used by the function.}
   \item{x}{For the distance heatmap, a \code{dist} object or a finite,
       nonnegative symmetric distance matrix with zero diagonal. For the
       principal-coordinate plots, a \code{stats::cmdscale} result calculated
       with \code{eig = TRUE}, containing \code{points} and \code{eig}.}
   \item{labels}{TRUE shows individual names in the distance heatmap;
       FALSE omits them. A character vector supplies one label per individual.}
   \item{cluster}{TRUE clusters the distance heatmap; FALSE retains input order.}
   \item{col}{Heatmap colours, group colours for coordinate scores, or the
       colour of variance-plot points. Group colours are recycled as needed.}
   \item{hclust.method}{Linkage method passed to \code{stats::hclust}.
       Clustering uses the supplied pairwise distances directly.}
   \item{which}{Coordinate scores, individual percentages of explained
       variance, or cumulative percentages of explained variance.}
   \item{axes}{Two different principal-coordinate indices for a score plot.}
   \item{groups}{NULL, or one nonmissing group per individual. If this
       vector has names and the coordinates have row names, groups are
       aligned to individuals by name. Otherwise their order must agree.}
   \item{pch}{Symbols for the groups in a score plot, recycled as needed,
       or the symbol for a variance plot.}
   \item{legend}{TRUE adds a group legend to a score plot when groups are supplied.}
   \item{cex}{Point size for principal-coordinate plots.}
   \item{...}{Named arguments passed to \code{stats::heatmap} or
       \code{graphics::plot}, such as margins, a title or axis labels.
       Data, distance scaling, clustering and group-aesthetic arguments
       must be supplied through the explicit arguments above.}
}

\details{
         The Euclidean distance is calculated as:\cr
         \deqn{D_E=\sqrt{\sum_{i=1}^m\sum_{j=1}^{n_j} (p_{ij}-q_{ij})^2}}
         
         where \eqn{p_{ij}} and \eqn{q_{ij}} are allele frequencies of the jth
         allele at the ith locus in the two individuals under consideration, \eqn{n_i}
         is the number of alleles at the ith locus, and \eqn{m} refers to the number of loci.
         The \eqn{D_E} ranges from zero to \eqn{\sqrt{2m}}, the limits being assumed when the
         two individuals have identical allele frequencies or are fixed for different alleles.\cr
         The \eqn{D_E} is appropriate if allelic informative marker data are available
         and the relationships between individuals are investigated in combination
         with multivariate methods that require dissimilarities possessing the Euclidean property (Reif 2005).\cr
         

         The Rogers distance (Rogers 1972) is a modification of \eqn{D_E} and is calculated as follows:\cr
         \deqn{D_R=\frac{1}{m}\sum_{i=1}^m \sqrt{\frac{1}{2}\sum_{j=1}^{n_j} (p_{ij}-q_{ij})^2}}
         
         The \eqn{D_R} is the average \eqn{D_E} across all loci standardized with the factor 1/2
         to restrict the values to the interval [0,1]. It is one only if two individuals
         are fixed for different alleles, but if one or both individuals are not fixed
         and they have no alleles in common, \eqn{D_R} is not equal to one.\cr
         Melchinger et al. (1991) derived theoretical results that \eqn{D_R} estimates
         between two homozygous inbreds are linearly related to the coancestry coefficient.
         Consequently, \eqn{D_R} is suitable for studying the relationship between
         the genetic dissimilarity of inbreds based on allelic informative marker
         data and the coefficient of coancestry.\cr

         
         The modified Rogers distance (Wright 1978) is calculated with the formula:\cr
         \deqn{D_W=\frac{1}{\sqrt{2m}}\sqrt{\sum_{i=1}^m\sum_{j=1}^{n_j} (p_{ij}-q_{ij})^2}}

         As an Euclidean distance with the values in [0,1] it can be used for the
         same applications as recommended for \eqn{D_E}. Like \eqn{D_R}, \eqn{D_W }
         is not equal to one in the case of multiple alleles, even if the two individuals
         have no allele in common.\cr
         \eqn{D_W} is especially suitable in studies based on allelic informative
         marker data for examining the prediction of heterosis with genetic dissimilarities
         or the establishment of heterotic groups. Furthermore, \eqn{D_W} can be used
         for the same applications as suggested for \eqn{D_E} (e.g., multivariate
         methods), owing to its Euclidean property.\cr
         

          The Rogers distance and the modified Rogers distance are
          standardized measures in the interval [0,1]. The standardization is
          achieved by the divisions by m and by 2 in the respective equations.
          This requests codominant markers! These measures can't be used with
          dominant markers. However the Euclidean distance is not standardized and
          can also be used with dominant marker data.\cr
          For homozygous genotypes Rogers distance equals the Nei-Li Distance
          (Nei and Li 1979), which is the same as 1 minus the Dice coefficient
          (Dice 1945), further the modified Rogers distance is the square root of the Rogers distance.

The distance heatmap uses the same ordering on both axes and never scales
rows or columns. With clustering enabled, hierarchical clustering is based
on the original pairwise genetic distances, not distances between rows of
the distance matrix. Its dendrogram can therefore differ from an unmodified
call to \code{heatmap(as.matrix(x))}.

The principal-coordinate helper draws already calculated coordinates and
does not repeat the distance calculation or ordination. Explained-variance
percentages divide each positive eigenvalue by the sum of positive
eigenvalues; negative eigenvalues are assigned zero for these display
percentages. For a non-Euclidean distance this convention does not measure
variation associated with negative axes. These helpers open no graphics
device and add no package dependencies.
}

\value{
\code{st.genetic.distances} invisibly returns either an object of class
\code{"dist"} when \code{format = "m"}, representing pairwise genetic
distances, or a long-format \code{data.frame} with columns \code{OTU1},
\code{OTU2}, and \code{Measure} when \code{format = "l"}.

\code{st.plot.genetic.distances} invisibly returns the \code{heatmap}
result, including row and column order indices. \code{st.plot.pcoa}
invisibly returns a list with the supplied \code{points}, the displayed
\code{variance.percent}, and the aligned \code{groups} when a score plot
is drawn.
}

\references{
           Dice, L.R. 1945. Measures of the amount of ecologic association
           between species. Ecology 26:197--302.

           Nei, M. and W.H. Li. 1979. Mathematical Models for studying genetic
           variation in terms of restriction endonucleases. Proc. Natl. Acad. Sci. USA. 76:5268--5371.
           
           Reif, J., Melchinger, A. and Frisch, M. 2005. Genetical and mathematical
           properties of similarity and dissimilarity coefficients applied in
           plant breeding and seed bank management. Crop Sci.45:1-7.

           Rogers, J.S. 1972. Measures of genetic similarity and genetic distance.
           VII. Univ. Tex. Publ. 7213:145--153.

           Wright, S. 1978. Evolution and the Genetics of Populations. Volume 4:
           Variability Within and Among Natural Populations. University of
           Chicago Press, Chicago, Illinois.
}


%\note{
%
%}
%
