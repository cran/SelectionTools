\name{load.effmap}
\alias{load.effmap}
\title{load.effmap()}

\description{ Loads an effect file }

\usage{
load.effmap(name, file = NA)
}

\arguments{
   \item{name}{Name of the effect}
   \item{file}{File that contains the effect description}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of loading the requested effect map.
}

\details{ Effect files are text files can be generated either by hand with
    an editor or with the command \code{generate.effect.file()}.

  An effect file consists of a value, which is assigned to each individual
  of a population (the population mean) and a list of effects which are
  added to the population mean if the certain allelic combinations occur
  in an individual. The order of the effects is arbitrary, but the
  population mean must occur first in the map file.\cr
  Mean;\cr

Structure of an effect file: \cr
  Locus Allele [Locus Allele \dots] Effect;\cr
  [,\dots]
  \tabular{ll}{
    Mean \tab The population mean \cr
    Classname \tab A class of Loci to which effects are assigned to \cr
    Allele \tab The alleles to which the effects are assigned \cr
    [Allele] \tab Optional. Is used to define dominance effects epistatic effects\cr
    [,\dots] \tab An arbitrary number of effects can be defined seperated by commas\cr
    } 
  }

\note{Further details appear in \code{define.effects} for an example }
