\name{define.effects}
\alias{define.effects}
\title{define.effects()}

\description{
  Defines effects that can be used for selection
}

\usage{
define.effects(name, description)
}

\arguments{
   \item{name}{Name of the effect}
   \item{description}{Description of the effect}
}
\value{
An integer status returned by the final \code{unlink()} call that removes the temporary effect-definition file. The meaningful result is the side effect of defining/loading the requested effects.
}

\details{
  The syntax for the \code{description} string is
  \cr
  Mean, Classname Allele [Allele] Distribution Parameter [,\dots]
  \cr
  \tabular{rl}{
    Mean \tab The population mean \cr
    Classname \tab A class of Loci to which effects are assigned to \cr
    Allele \tab The allels to which the effects are assigned to \cr
    [Allele] \tab Optional. Is used to define dominance effects \cr
    [, \dots] \tab An arbitrary number of effects can be defined
    separated by commas \cr
    Distribution \tab normal 0/0.01 or uniform 0
  }
  }
