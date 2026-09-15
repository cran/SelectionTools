\name{gd.genetic.distance}

\alias{gd.genetic.distance}

\title{gd.genetic.distance()}

\description{
Calculates pairwise genetic distances between all pairs of populations
in the dataset.
}

\usage{
gd.genetic.distance(dta, measure="euc", par="dist", boot=0, no.pop.u=0,
                    pop.u=0, no.ind.u=0, ind.u=0, no.mar.u=0, mar.u=0)
}

\arguments{
   \item{dta}{Name of the dataset containing the observed marker
   genotypes. The format of the dataset is illustrated by the datasets
 gd.dummy.populations and gd.maize.populations.}
   \item{measure}{\code{"mrd"} is the Modified Rogers Distance according to
   Wright (1978, p. 78). \code{"rd"} is the Rogers Distance (Rogers
   (1972). \code{"euc"} is the Euclidean Distance with respect to the
   allele frequencies, sometimes this distance measure is referred to as
   Gowers Distance. If measure is ommitted the Euclidean distance is computed.}
   \item{par}{Parameters to be calculated:
   \tabular{ll}{
     \code{"dist"} \tab Genetic distance measure \cr
     \code{"sdev/jackm"} \tab Jackkniving over markers \cr
     \code{"sdev/bootm"} \tab Bootstrapping over markers \cr
     \code{"sdev/booti"} \tab Bootstrapping over individuals \cr
     \code{"sdev/bootmi"} \tab Bootstrapping over markers and individuals \cr
     }}
   \item{boot}{Number of bootstrap runs for estimating the standard
   deviation of the  estimates.}
   \item{no.pop.u}{Function argument.}
   \item{pop.u}{Function argument.}
   \item{no.ind.u}{Function argument.}
   \item{ind.u}{Function argument.}
   \item{no.mar.u}{Function argument.}
   \item{mar.u}{Function argument.}
}


\value{
A \code{data.frame} with columns \code{OTU1}, \code{OTU2}, and \code{Measure}, containing the requested pairwise genetic distances; an additional \code{SDev} column is included when resampling-based standard deviations are requested. Returns \code{NULL} invisibly on validation failure.
}

\note{

\emph{1.) } The Rogers distance and the Modified Rogers distance are
standardized measures in the interval [0,1]. The standardization is
achieved by the divisions by m and by 2 in the respective equations.
This requests codominant markers! These measures can't be used with
dominant markers. However the Euclidean distance is not standardized and
can also be used with dominant marker data (such as the AFLP data in the
dataset \code{gd.maize.aflp}).


  
\emph{2.) } For homozygous genotypes Rogers' Distance (Rogers 1972)
equals the Nei-Li Distance (Nei and Li 1979), which is the same as 1
minus the Dice coefficient (Dice 1945), further the Modified Rogers'
Distance is the square root of the Rogers' Distance. 


\emph{3.) }
For calculating pairwise genetic distances between individuals you can
use populations consisting of one single individual each. The column
names could be e.g., \code{ p1.1 p2.1 p3.1 p4.1}

}

\references{
Dice, L.R. 1945. Measures of the amount of ecologic
  association between species. Ecology 26:197--302.
  
Nei, M. and W.H. Li. 1979. Mathematical Models for
studying genetic variation in terms of restriction
endonucleases. Proc. Natl. Acad. Sci. USA. 76:5268--5371.  
  
Rogers, J.S. 1972. Measures of genetic similarity and
genetic distance. VII. Univ. Tex. Publ. 7213:145--153.

Wright, S. 1978. Evolution and the Genetics of Populations. Volume 4:
Variability Within and Among Natural Populations. University of
Chicago Press, Chicago, Illinois.
}
