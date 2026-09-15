% copyright (C) 2001 Frisch/Maurer
%
\name{st.genetic.distances}
\alias{st.genetic.distances}
\title{st.genetic.distances()}
\description{
  The function \emph{st.genetic.distances} calculates genetic distances for each possible
  pair of individuals.
}
\usage{
st.genetic.distances(measure = "mrd", format = "l",
                     filename = NULL, auxfiles = FALSE,
                     data.set = "default")
}
\arguments{
   \item{measure}{character. Defines the measure to be used for calculation. Incorporated
                measures are modified Rogers distance ("mrd"), Rogers distance ("rd")
                and Euclidean distance ("euc").}
   \item{format}{character. Specifies the format of the output file. Set to "l" for
               long format or to "m" for matrix format.}
   \item{filename}{Character string giving the base name of the genetic-distance file. Required when \code{auxfiles = TRUE}; the file is written with extension \file{.gdi} in the output directory.}
   \item{auxfiles}{Logical. If \code{TRUE}, a persistent output file with the genetic distances is generated and \code{filename} must be supplied.}
   \item{data.set}{character. Name of the data set that is used by the function.}
}

\details{
         The Euclidean distance is calculated as:\cr
         \deqn{D_E=\sqrt{\sum_{i=1}^m\sum_{j=1}^{n_j} (p_{ij}-q_{ij})^2}}
         
         where \eqn{p_{ij}} and \eqn{q_{ij}} are allele frequencies of the jth
         allele at the ith locus in the two individuals under consideration, \eqn{n_i}
         is the number of alleles at the ith locus, and \eqn{m} refers to the number of loci.
         The \eqn{D_E} ranges from zero to 2m, the limits being assumed when the
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
}

\value{
Invisibly returns either an object of class \code{"dist"} when \code{format = "m"}, representing pairwise genetic distances, or a long-format \code{data.frame} with columns \code{OTU1}, \code{OTU2}, and \code{Measure} when \code{format = "l"}.
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
