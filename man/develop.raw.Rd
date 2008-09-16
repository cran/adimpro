\name{develop.raw}
\alias{develop.raw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create a color image from RAW image data. }
\description{
  The function generates a color image from RAW image data.
}
\usage{
develop.raw(object, method = "BILINEAR", wb = c(1, 1, 1), maxrange= TRUE, compress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class \code{"adimpro"} containing RAW image data, e.g. created by function \code{\link{read.raw}}. Such an object is characterized by \code{object$type=="RAW"}. }
  \item{method}{ Method to fill missing color values. Currently implemented are \code{method="HALF"} (reduce image size by factor of 2),
  \code{method="FULL"} (use color from neighbor within (2x2) Bayer mask)),
  \code{method="BILINEAR"} (bilinear interpolation), \code{method="Median4"} and \code{method="Median16"}, the latter two being based on L1-MM over shifted Bayer masks. \code{method="Median16"} delivers much smoother results, but is considerably slower then the other methods. }
  \item{wb}{ numerical vector of length 3 containing multiplicative factors for the three color channels. }
  \item{maxrange}{If TRUE increase range of values to maximum.}
  \item{compress}{ logical, determines if image data are stored as \code{raw}.} 
}
\details{
  RAW image data usually contain only one color value at each pixel,
  with colors arranged in a so called Bayer mask. Converting RAW images into color images requires to fill the missing entries in the color
  channels. 
}
\value{
  object of class "adimpro" containing the image.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de} }
\note{ The function requires dcraw to be installed.}
\seealso{ \code{\link{read.image}} }
\examples{
\dontrun{demo(io)}
\dontrun{demo(raw)}
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
