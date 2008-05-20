\name{rotate.image}
\alias{rotate.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Image Processing }
\description{
  Rotate an image by 0, 90, 180 or 270 degrees.
}
\usage{
rotate.image(img, angle = 90, compress=NULL)
}
\arguments{
  \item{img}{ image data, that is an object of class "adimpro". }
  \item{angle}{ 0, 90, 180 or 270 degrees }
  \item{compress}{ store result as compressed image if \code{compress=TRUE}. \code{compress=NULL}
   uses the same format as the original object. }
}
\details{
  The function rotates the image \code{img} by 0, 90, 180 or 270
  degrees. Any other value for \code{degree} will cause an exception.
  The returned object contains an additional component \code{rotate} 
  determining the rotation.
}
\value{
  An object of class "adimpro" containing the rotated image.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\seealso{ \code{\link{show.image}} }
\examples{
if(Sys.getenv("ImageMagick")==""){
cat("Please install ImageMagick\n")
} else {
img <- read.image(system.file("img/wias.ppm",package="adimpro"))
X11(height=5,width=10)
par(mfrow=c(1,2))
show.image(img,main="Original Image")
img.rot <- rotate.image(img)
show.image(img.rot,main="Rotated Image")
}
}
\keyword{ manip }% at least one, from doc/KEYWORDS
