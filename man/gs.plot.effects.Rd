% copyright (C) 2001 Frisch/Maurer
%
\name{gs.plot.effects}
\alias{gs.plot.effects}
\title{gs.plot.effects()}
\description{
Plotting function for the estimated marker effects.
}

\usage{
gs.plot.effects(data.set = "default", absolute = TRUE, alpha = 0.05,
                p.adjust.method = "none", pvals = NULL, xlim = NULL,
                ylim = NULL, plt = TRUE, col.sp = "blue", pch.sp = 19,
                col.sn = "blue", pch.sn = 19, col.ns = "blue", pch.ns = 1,
                ...)
}
\arguments{
   \item{data.set}{character. Name of the data set that is used by the function.}
   \item{absolute}{Function argument.}
   \item{alpha}{Function argument.}
   \item{p.adjust.method}{Function argument.}
   \item{pvals}{Function argument.}
   \item{xlim}{Function argument.}
   \item{ylim}{Function argument.}
   \item{plt}{Function argument.}
   \item{col.sp}{Function argument.}
   \item{pch.sp}{Function argument.}
   \item{col.sn}{Function argument.}
   \item{pch.sn}{Function argument.}
   \item{col.ns}{Function argument.}
   \item{pch.ns}{Function argument.}
   \item{...}{Additional arguments.}
}

%\details{
%
%}

\value{
Invisibly returns the \code{data.frame} used for plotting, with columns \code{chrom}, \code{pos}, \code{name}, \code{cpos}, and \code{effect}, and a \code{p} column when p-values are available. Returns \code{NULL} if the required effect or map data are unavailable.
}

%\note{
%
%}
%
