\name{.st_predict_method}
\alias{.st_predict_method}
\title{Dispatches compiled Prediction Methods}
\description{
Dispatches summaries, coefficient and covariance accessors, and printing
for random-effect prediction objects to the compiled implementation.
}
\usage{
.st_predict_method(object, action, ...)
}
\arguments{
\item{object}{An \code{st.predict} prediction object.}
\item{action}{A character string identifying the operation to perform.}
\item{...}{Named arguments required by the selected operation.}
}
\value{
A prediction summary table, a vector of predicted means or their
prediction-error covariance matrix, according to the selected operation.
Printing returns the supplied prediction object.
}
\keyword{internal}
