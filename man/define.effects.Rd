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
An integer compatibility status of zero. The main effect is defining and loading the requested effects directly in memory.
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

  For \code{normal mean/variance}, the first number is the mean and the
  second number is the variance of the normal distribution from which an
  effect is drawn independently for every matching locus. Thus
  \code{normal 1/0.05} uses mean 1 and variance 0.05 (standard deviation
  \eqn{\sqrt{0.05}}). The variance must be nonnegative.

  For \code{uniform value}, the supplied number is used as a fixed effect
  value at every matching locus; despite the historical keyword
  \code{uniform}, no random draw is made for this form.

  If one allele is supplied, the specification defines an allele effect.
  If two alleles are supplied before the distribution keyword, the
  specification defines the effect of that allele combination and can be
  used for dominance effects. Several effect specifications can be included
  in one description, separated by commas.
  }
