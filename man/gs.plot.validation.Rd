% copyright (C) 2001 Frisch/Maurer
%
\name{gs.plot.validation}
\alias{gs.plot.validation}
\title{gs.plot.validation()}
\description{
Creates a graphical output wit two plots. The left plot shows the model fit as produced
with the function \emph{gs.plot.model.fit}. The right plot shows a scatter plot of the predicted
DGV and the real phenotypes for a set of individuals that was not used for estimation of marker effects.
}
\usage{
gs.plot.validation(training.set, prediction.set, title="")
}
\arguments{
   \item{training.set}{character. Name of the data set that was used for estimation
                       of marker effects.}
   \item{prediction.set}{character. Name of the data set that is used for prediction of DGVs
                       within this function.}
   \item{title}{character. Title of the plot.}
}
\value{
No useful return value is intended. The function is called for its side effect of plotting validation information; invalid input returns \code{NULL} invisibly.
}

\details{
Within this function, the function \emph{gs.estimate.gv} is called twice. For the model fit,
the estimation set is both used as estimation and prediction set in \emph{gs.estimate.gv}.
To validate the model, the function \emph{gs.estimate.gv} is called again with the
validation set used as prediction set.\cr
The validation set is a set of individuals that was not used for building the model,
but having phenotypes to be correlated to the predicted DGVs.
}


%\note{
%
%}
%
