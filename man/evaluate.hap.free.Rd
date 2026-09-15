\name{evaluate.hap.free}
\alias{evaluate.hap.free}
\title{Release or reset temporary haplotype-evaluation state}

\description{
Release or reset temporary haplotype-evaluation state.
}

\usage{
evaluate.hap.free()
}
\value{
An invisible \code{list} returned by the underlying \code{.C} routine. The list is an implementation-level result containing the arguments returned by the compiled routine; the function is primarily called for its side effect of releasing haplotype evaluation state.
}
