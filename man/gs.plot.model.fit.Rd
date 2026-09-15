% copyright (C) 2001 Frisch/Maurer
%
\name{gs.plot.model.fit}
\alias{gs.plot.model.fit}
\title{gs.plot.model.fit()}
\description{
Direct genomic values are predicted for those individuals that were used for
building the model and estimating marker effects. A scatter plot is created showing
the true breeding values (y, phenotypes) and the predicted DGVs (yhat).
}

\usage{
gs.plot.model.fit(training.set ="default", title = "")
}
\arguments{
   \item{training.set}{character. Name of the data set to be used for prediction.}
   \item{title}{character. Title of the plot.}
}

%\details{
%
%}

\value{
No useful return value is intended. The function is called for its side effect of plotting model-fit information; invalid input returns \code{NULL} invisibly.
}
%\note{
%
%}
%
