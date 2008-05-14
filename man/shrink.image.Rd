\name{shrink.image}
\alias{shrink.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Image Processing }
\description{
  Shrink resolution of an image
}
\usage{
shrink.image(img, method = "gap",
             xt = img$dim[1], yt = img$dim[2], ratio = TRUE,
             compress=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{ image data, an object of class "adimpro". }
  \item{method}{ method to be used to shrink the image. "gap" (default),
    "mean", or "nearest" }
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

  \code{method} "gap" selects a pixel nearest to the new coordinates.
  \code{method} "mean" defines the color of a pixel as the mean of 
  all pixel identified with the new coordinate.
  \code{method} "nearest" set the color of a pixel as color of an 
  representative pixel from the original image.
}
\value{
  shrinked image
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de} }
\examples{
img <- read.image(system.file("img/wias.ppm",package="adimpro"))
X11(height=3,width=9)
par(mfrow=c(1,3))
show.image(img,main="Original Image")
img.shrink1 <- shrink.image(img,method="gap")
img.shrink2 <- shrink.image(img,method="nearest",xt=20,yt=20)
show.image(img.shrink1,main="Shrinked Image (gap)")
show.image(img.shrink2,main="Shrinked Image (nearest)")
}
\keyword{ manip }% at least one, from doc/KEYWORDS
