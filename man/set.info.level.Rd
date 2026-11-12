% file set.info.level.Rd
% copyright (C) 2001 Frisch/Maurer
%
\name{set.info.level}
\alias{set.info.level}
\title{set.info.level()}
\description{
  The function \emph{set.info.level()} sets the information level. The
  smaller the level number is the less program information is printed on
  the screen.
  \tabular{rl}{
    level\tab information printed\cr
    -2\tab error messages\cr
    -1\tab warning messages\cr
    0\tab default messages\cr
    1\tab verbouse mode
  }
}
\usage{
set.info.level(level)
}
\arguments{
   \item{level}{information level}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of setting the simulation information level.
}
\note{
  If you want to read only error messages then you have to set your
  information level to -2!
}
