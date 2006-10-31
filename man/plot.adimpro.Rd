\name{plot.adimpro}
\alias{plot.adimpro}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O functions }
\description{
  Visualize image data.
}
\usage{
\method{plot}{adimpro}(x, new=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{image object of class "adimpro" }
  \item{new}{ should new X11() be opened? default \code{FALSE}}
  \item{...}{ not used }
}
\details{
  This functions shows information on the image. This includes
  histograms of color values in each channel of the specified (\code{x$type}) color space, a thumbnail,
  some information on the image and, if x was produced by \code{awsimage},
  an image illustrating the local adaptation.
}
\value{
  nothing is returned.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{summary.adimpro}} }
\examples{
img.noise <- read.image(system.file("img/wias-noise.ppm",package="adimpro"))
plot(img.noise)
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ iplot }% at least one, from doc/KEYWORDS
