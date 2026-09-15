\name{set.co.freq}
\alias{set.co.freq}
\title{Set crossover-frequency information}

\description{
Set crossover-frequency information.
}

\usage{
set.co.freq(cofreq=1)
}

\arguments{
   \item{cofreq}{Function argument used by this operation.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting crossover-frequency information.
}
