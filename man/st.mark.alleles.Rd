\name{st.mark.alleles}
\alias{st.mark.alleles}
\title{Identify markers by effect and replace selected marker alleles in marker data}

\description{
Identify markers by effect and replace selected marker alleles in marker data.
}

\usage{
st.mark.alleles(markers = NULL, data.set="default", effect.set=NULL,
                target.set, hap.lst = NULL, alpha=0.05,
                p.adjust.method="none", replace = "3")
}

\arguments{
   \item{markers}{Function argument used by this operation.}
   \item{data.set}{Name of the SelectionTools data set.}
   \item{effect.set}{Data set supplying marker effects.}
   \item{target.set}{Target data set.}
   \item{hap.lst}{Function argument used by this operation.}
   \item{alpha}{Method-specific tuning or significance parameter.}
   \item{p.adjust.method}{Function argument used by this operation.}
   \item{replace}{Logical or coding value controlling replacement.}
}
\value{
On successful completion, invisibly returns the character value used to restore the previous SelectionTools input directory. Validation failures return \code{NULL} invisibly. The main effect is the requested allele recoding and marker-data update.
}
