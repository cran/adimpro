\name{shrink.image}
\alias{shrink.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Image Processing }
\description{
  Shrink resolution of an image
}
\usage{
shrink.image(img, method = "median",
             xt = img$dim[1], yt = img$dim[2], ratio = TRUE,
             compress=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{ image data, an object of class "adimpro". }
  \item{method}{ method to be used to shrink the image. "median" (default),
    "mean", or "nearest". "median" is supposed to give best results. For a considerably faster result use "nearest". }
  \item{xt}{ target x-dimension }
  \item{yt}{ target y-dimension }
  \item{ratio}{ logical. preserve x-y ratio? default: \code{TRUE} }
  \item{compress}{ logical, determines if image data are stored in raw-format. }
}
\details{
  This function shrinkes the resolution of the image such that the x-y
  dimension of the resulting image is smaller than the original one.
  \code{xt} and \code{yt} give the target dimension
  of the image. If \code{ratio == TRUE} (default) the ratio between x- and
  y-dimension is preserved.

  \code{method} "nearest" selects a pixel nearest to the new coordinates.
  \code{method} "mean" defines the color of a pixel as the mean of 
  all pixel identified with the new coordinate.
  \code{method} "median" set the color of a pixel as color of an 
  L1-generalized median of all pixel identified with the new coordinate.
}
\value{
  shrinked image
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de} }
\examples{
\dontrun{demo(manipulate)}
}
\keyword{ manip }% at least one, from doc/KEYWORDS
