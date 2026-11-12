\name{generate.effect.file}
\alias{generate.effect.file}
\title{generate.effect.file()}

\description{ Generates a file with effect descriptions that can be used
  for selection.  Effects can be loaded.  Additive and dominance effects
  can be }

\usage{
generate.effect.file(fName, description)
}

\arguments{
   \item{fName}{Name of the file}
   \item{description}{Description of the effect}
}
\value{
An invisible \code{list} returned by the compiled routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of generating the effect file.
}

\note{Further details appear in \code{define.effects} for the syntax of the
  descripton and for an example }
