\name{show.image}
\alias{show.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O functions }
\description{
  Display an image on the screen.
}
\usage{
show.image(img, max.x = 1000, max.y = 1000, color.par = NULL, channel=NULL,
           alg = 1, new = FALSE,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{ image data, an object of class "adimpro".}
  \item{max.x}{ maximum value of pixels in x dimension to be displayed. }
  \item{max.y}{ maximum value of pixels in y dimension to be displayed. }
  \item{color.par}{ list of parameters for the gamma correction and some
    white balance. List entries:
    \tabular{rll}{
      ga \tab gamma value \tab default 2.4\cr
      bp  \tab break point \tab default 0.00304\cr
      alg \tab algorithm for gamma correction \tab default 1\cr
      red \tab factor to multiply with red channel. \tab default 1.0\cr
      blue \tab factor to multiply with blue channel. \tab default 1.0\cr
      brightness \tab factor to multiply with all channels \tab default 1.0}
  }  
  \item{channel}{ allows to select a color channel (1: red, 2: green, 3: blue) for display,
   for rgb-images only.}
  \item{alg}{ chooses computing algorithms which differ in speed and precision. }
  \item{new}{ should new X11() be opened? default \code{FALSE}}
  \item{\dots}{ additional arguments to \code{\link{image}} can be
    passed here. }
}
\details{
  This function displayes greyscale and color images
  on the screen. If the actual dimension of the
  image exceeds \code{max.x} or \code{max.y} the image is shrinked by
  displaying only part of the pixels (every second/third/... such that
  the resulting dimension is smaller than \code{max.x} or \code{max.y})

  If \code{color.par$red}, \code{color.par$blue} or \code{color.par$brightness} differ from 1.0 the
  corresponding channels are multiplied with the provided
  values. Saturated values are set to 1.

  If \code{img$gamma} is FALSE, perform gamma correction with
  \code{color.par$ga} and \code{color.par$bp}. \code{alg} chooses between three different
  computing algorithms (approximations) with increasing computation time
  and precision (\code{alg} is 1,2, or 3).

}
\value{
  Nothing is returned
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{read.image}}, \code{\link{write.image}} }
\examples{
img <- read.image(system.file("img/wias.ppm",package="adimpro"))
show.image(img)
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
