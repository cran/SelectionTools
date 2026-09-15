\name{generate.map.file}
\alias{generate.map.file}
\title{generate.map.file()}

\description{
Generates a simulation linkage-map text file from a compact map description without itself installing the map.
}

\usage{
generate.map.file(fName, description)
}

\arguments{
   \item{fName}{The name of the map file that is generated}
   \item{description}{Description of the linkage map}
}
\value{
An invisible \code{list} whose \code{retval} component is the compiled map-generation status. Invalid map dimensions return \code{list(retval = -2L)} invisibly. The main effect is writing the generated map file.
}

\details{
The description syntax is the same syntax consumed by \code{define.map}. The function validates the file name and description, normalizes a multi-element description to one comma-separated specification, and asks the simulation backend to write the map file.

The generated file can subsequently be loaded as a simulation linkage map. Generating the file does not by itself replace the currently loaded map.
}
