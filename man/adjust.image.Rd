\name{adjust.image}
\alias{adjust.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Image Processing }
\description{
  Adjust color channels and apply gamma correction (if applicable).
}
\usage{
adjust.image(img, color.par = NULL, alg = 1, compress= TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{ image object, class"adimpro". }
  \item{color.par}{ list of parameters for the gamma correction and some
    white balance. list entries:
    \tabular{rll}{
      ga \tab gamma value \tab default 2.4\cr
      bp  \tab break point \tab default 0.00304\cr
      alg \tab algorithm for gamma correction \tab default 1\cr
      red \tab factor for red channel. \tab default 1.0\cr
      blue \tab factor for blue channel. \tab default 1.0\cr
      brightness \tab factor for all channels \tab default 1.0}
  }  
  \item{alg}{ determines the approximation for the gamma correction. Select 1 for fastest 
  computation and 3 for maximum accuracy, or 2 for a compromise. }
  \item{compress}{ logical, determines if image data are stored in raw-format. }
}
\details{
  This function adjusts color channels and applies gamma correction (if applicable).

  If \code{color.par$red} or \code{color.par$blue} or \code{color.par$brightness} differ from 1.0 the
  corresponding channels are multiplied with the provided
  values. Saturated values are set to 1.

  If \code{img$gamma==FALSE}, perform gamma correction with
  \code{color.par$ga} and \code{color.par$bp}. \code{alg} chooses between three different
  computing algorithms (approximations) with increasing computation time
  and precision (\code{alg} is 1,2, or 3).
}
\value{
  Adjusted image object of class "adimpro".
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{show.image}}, \code{\link{write.image}} }
\examples{
img <- read.image(system.file("img/wias.ppm",package="adimpro"))
img2 <- adjust.image(img,color.par=list(red=1.5,blue=2))
show.image(img2)
}
\keyword{ manip }% at least one, from doc/KEYWORDS
