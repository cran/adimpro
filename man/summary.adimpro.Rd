\name{summary.adimpro}
\alias{summary.adimpro}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O functions }
\description{
  'summary' method for class '"adimpro"'.
}
\usage{
\method{summary}{adimpro}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of class \code{adimpro}, usually, a result of a call to
    \code{read.raw}, or \code{read.image}. }
  \item{\dots}{ further arguments passed to or from other methods. }
}
\details{
  The method tries to print information on the image, like image dimension,
  color space, value range, etc.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\examples{
\dontrun{demo(io)}
\dontrun{demo(manipulate)}
}
\keyword{ print }% at least one, from doc/KEYWORDS
