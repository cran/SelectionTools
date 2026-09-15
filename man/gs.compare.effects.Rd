\name{gs.compare.effects}
\alias{gs.compare.effects}
\title{gs.compare.effects()}
\description{
Creates a scatter plot of effects of two different models to visualize differences in estimated effects.
}

\usage{
gs.compare.effects(data.set.a, data.set.b, label.a = data.set.a,
                   label.b = data.set.b)
}
\arguments{
   \item{data.set.a}{character. Name of the first data set containing marker effects.
                  These effects are plotted on the x axis.}
   \item{data.set.b}{character. Name of the second data set containing marker effects.
                  These effects are plotted on the y axis.}
   \item{label.a}{character. Label of the x axis.}
   \item{label.b}{character. Label of the y axis.}
}
\value{
No useful return value; the function is called for its side effect of comparing and plotting marker-effect results. Validation failures return \code{NULL} invisibly.
}

%\details{
%
%}

%\note{
%
%}
