\name{imganiso2D}
\alias{imganiso2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ create an image that visualizes anisotropy}
\description{
  The function creates an object of class \code{adimpro} that
  visualizes anisotropy information  using the HSI color space for 
  main direction of anisotropy (H), \code{maximum eigenvalue^satexp} (S)
  and log(excentricity) (I).
}
\usage{
imganiso2D(x, satexp = 0.25, g=3, rho=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a field of 2D tensors (Dimension \code{c(3,n1,n2)}) or an object of class \code{"adimpro"}.  }
  \item{satexp}{ exponent for maximum eigenvalue in saturation channel. Determines the contrast in this channel.}
  \item{g}{ Bandwidth for anisotropic smoothing gradient estimates,
    preferably \eqn{g \ge 3}{g >= 3} for images with line type texture and small 
    \eqn{g \approx 1}{g approx 1} for improving edges between homogeneous regions. }
  \item{rho}{ Regularization parameter for anisotropic smoothing gradient estimates,
    preferably \eqn{\rho = 0}{rho = 0} for images with line type texture and large 
    \eqn{\rho \approx 3*variance}{rho approx 3*variance} for improving edges between homogeneous regions. }
}
\value{
an object (image) of class \code{adimpro}.
}
\references{
  Polzehl, J. and Tabelow, K. (2007). Adaptive smoothing of digital
  images, Journal of Statistical Software 19 (1). 
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\seealso{ \code{\link{awsaniso}},\code{\link{read.image}}, \code{\link{read.raw}}, \code{\link{make.image}}, \code{\link{show.image}}, \code{\link{clip.image}} }
\examples{
\dontrun{demo(awsimage)}
\dontrun{demo(manipulate)}
}
\keyword{ manip }% at least one, from doc/KEYWORDS
