\name{gs.make.A}
\alias{gs.make.A}
\title{Constructs a genomic relationship matrix from marker dosages}
\description{
Constructs the centred and scaled VanRaden relationship matrix from a
matrix with individuals in rows and markers in columns.
}
\usage{
gs.make.A(Z, coding = "0 1 2")
}
\arguments{
\item{Z}{Nonempty numeric marker matrix without missing or infinite values.}
\item{coding}{Marker coding: \code{"0 1 2"}, \code{"-1 0 1"},
or \code{"0 0.5 1"}. Values must belong to the specified coding.
\code{"detect"} requires exactly one compatible coding.}
}
\details{
After conversion to diploid allele dosage, frequencies are column means
divided by two. With \eqn{W_{ij}=Z_{ij}-2p_j} and
\eqn{c=2\sum_j p_j(1-p_j)}, the result is \eqn{A=WW^T/c}.
The scale must be positive. Monomorphic columns contribute zero.
The matrix can be singular; no inverse or diagonal adjustment is applied.
Neither the input nor an internal SelectionTools data set is modified.
}
\value{
A numeric square matrix in input individual order. If \code{Z} has row
names, they label both dimensions. Returns \code{NULL} with a diagnostic
for invalid input or a relationship-matrix construction failure.
}
