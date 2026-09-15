\name{ssd.mating}
\alias{ssd.mating}
\title{ssd.mating()}

\description{
Generates single-seed-descent progeny by repeated selfing of each individual in a diploid parental population.
}

\usage{
ssd.mating(NamePg, NameP, NoPg=1, maxcycles=1)
}

\arguments{
\item{NamePg}{Character string naming the progeny population to create.}
   \item{NameP}{Character string naming the parental population.}
   \item{NoPg}{Positive number of independently generated descendants per parental individual.}
   \item{maxcycles}{Positive maximum number of successive selfing generations for each descendant.}
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of creating progeny using the SSD mating routine.
}

\details{
For every active individual in \code{NameP}, \code{NoPg} descendant lines are generated. The total progeny size is therefore the parental population size multiplied by \code{NoPg}. Descendants belonging to one parent occupy consecutive positions in the new population.

Each descendant is repeatedly selfed for at most \code{maxcycles} generations. Selfing stops earlier when the simulated chromosomes have become completely homozygous. Large values of \code{maxcycles} can therefore be used to approach fully homozygous recombinant inbred lines without requiring unnecessary generations after homozygosity is reached.

The routine requires a diploid simulation genome. The destination is rebuilt from chromosome data and does not retain stale cached genotypes or evaluated values. Allocation or meiosis failure removes the incomplete destination population.
}
