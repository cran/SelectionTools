\name{define.map}
\alias{define.map}
\title{define.map()}

\description{
Defines and installs a simulation linkage map from a compact textual map description.
}

\usage{
define.map(description)
}

\arguments{
   \item{description}{A string containing a description of the linkage map}
}
\value{
An invisible \code{list} returned by \code{generate.map.file()} or \code{linkage.map.load()}. Its \code{retval} component is the status code of the map-generation/loading step; the main effect is defining the simulation map.
}

\details{
The description is first converted to a temporary map file by \code{generate.map.file} and is then loaded through the simulation linkage-map loader. The operation therefore uses the same map validation and ordering rules as file-based simulation-map loading.

Description elements specify a locus name, class, and position rule. A position can identify a chromosome and map position directly or can request generated placement such as random, terminal, or nonterminal positions according to the existing description syntax. A multiplicity prefix can request several loci from one specification.

The function replaces the currently active simulation linkage-map definition only after its generated description has been accepted by the underlying map-generation and loading operations.
}
