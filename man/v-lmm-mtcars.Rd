\name{v-lmm-mtcars}
\alias{v-lmm-mtcars}
\alias{mtcars}
\docType{data}
\title{Motor Trend Car Measurements for Regression}
\description{Fuel consumption and vehicle characteristics for 32 cars,
used for ordinary regression and fixed-effect subset comparison.}
\usage{data("v-lmm-mtcars", package = "SelectionTools")}
\format{A \code{data.frame} named \code{mtcars} with 32 rows and 11 columns.
Every column is stored as a numeric double. Car names are row names.
\describe{
\item{\code{mpg}}{Fuel economy in miles per US gallon.}
\item{\code{cyl}}{Number of engine cylinders.}
\item{\code{disp}}{Engine displacement in cubic inches.}
\item{\code{hp}}{Gross horsepower.}
\item{\code{drat}}{Rear axle ratio.}
\item{\code{wt}}{Vehicle weight in thousands of pounds.}
\item{\code{qsec}}{Quarter-mile time in seconds.}
\item{\code{vs}}{Engine configuration: 0 for V-shaped and 1 for straight.}
\item{\code{am}}{Transmission: 0 for automatic and 1 for manual.}
\item{\code{gear}}{Number of forward gears.}
\item{\code{carb}}{Number of carburettors.}
}}
\details{The observations concern 1973--1974 car models and were originally
compiled from the US magazine \emph{Motor Trend}. The data contain no missing
values. Original values, row names, column names and row order are retained.
The coded variables \code{cyl}, \code{vs}, \code{am}, \code{gear} and
\code{carb} remain numeric; they are not automatically converted to factors.

The numeric values match the Dryad release \emph{Car parameters},
dedicated to the public domain under CC0 1.0.}
\source{Henderson H, Velleman PF (2016). Car parameters. Dryad dataset.
DOI: 10.15146/R32S37. File \file{CarParameters.tsv}, released 28 September 2016.

Henderson HV, Velleman PF (1981). Building Multiple Regression Models
Interactively. \emph{Biometrics} 37:391--411.}
\value{Loading the dataset creates the \code{mtcars} data frame.}
\keyword{datasets}
