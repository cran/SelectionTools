% copyright (C) 2001 Frisch/Maurer
%
\name{gs.plot.model.fit}
\alias{gs.plot.model.fit}
\title{Plot the Fit of a Genomic Prediction Model}
\description{
Plots observed performance values against genomic predictions for the
individuals used to fit the marker model.
}

\usage{
gs.plot.model.fit(training.set ="default", title = "")
}
\arguments{
   \item{training.set}{Character string naming the data set used for marker-effect estimation and model-fit assessment.}
   \item{title}{Character string added to the plot title.}
}

\details{
The function calls \code{gs.predict.genotypes} with \code{training.set} used
as both training and prediction set. The plotted \code{y} values are the
observed performance values stored in that data set, whereas \code{yhat} are
the direct genomic predictions obtained from the estimated marker effects.
In the simulation workflow, these \code{yhat} values represent estimated
genotypic values of the training population.

The plot includes the identity line, horizontal and vertical lines through
the respective means, and reports both the Pearson correlation \code{r} and
the Spearman rank correlation \code{rho} between \code{y} and \code{yhat}.
Marker effects must have been estimated before this function is called, and
performance data must be available in \code{training.set}.
}

\value{
No useful return value is intended. The function is called for its side effect of plotting model-fit information; invalid input returns \code{NULL} invisibly.
}

\seealso{
\code{gs.esteff.rr}, \code{gs.predict.genotypes}, \code{gs.plot.validation}
}
%
