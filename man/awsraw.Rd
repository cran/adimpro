\name{awsraw}
\alias{awsraw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothing and demosaicing of RAW images }
\description{
  The function integrates smoothing and demosaicing of RAW image data.
}
\usage{
awsraw(object, hmax = 4, aws = TRUE, wb = c(1, 1, 1), cspace = "Adobe", 
       ladjust = 1, maxrange=TRUE, lkern = "Triangle", graph = FALSE, 
       max.pixel = 400, compress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of class \code{adimpro} containing the RAW image data. 
                See \code{\link{read.raw}} for creating such objects.  }
  \item{hmax}{ maximal bandwidth to use in the smoothing algorithm.}
  \item{aws}{ use adaptive weights if \code{aws==TRUE}.}
  \item{wb}{ Vector containing factors for the three color chanels, allows to change the white balance. }
  \item{cspace}{ Color space of the result,  }
  \item{ladjust}{ Factor for the critical value \eqn{\lambda}{lambda}. Defaults to \code{1}, smaller values 
   increase sensitivity but may result in isolated noisy pixel. Values larger than \code{1} give smoother up to cartoon
   like results.}
  \item{maxrange}{If TRUE increase range of values to maximum.}
  \item{lkern}{ Specifies the location kernel. Defaults to
    "Triangle", other choices are "Quadratic", "Cubic" and "Uniform".
    The use of "Triangle" corresponds to the Epanechnicov kernel 
    nonparametric kernel regression. }
  \item{graph}{ (logical). If  \code{graph=TRUE} intermediate results are
    illustrated after each iteration step. Defaults to \code{FALSE}. }
  \item{max.pixel}{ Maximum dimension of images for display
    if \code{graph=TRUE}. If the true dimension is larger, the 
    images are downscaled for display. See also \code{\link{show.image}}.}
  \item{compress}{ logical, determines if image data are stored in raw-format. }
}
\details{
   Adaptive smoothing is performed on the original RAW data, restricting positive weights to
   pixel corresponding to the same color channel. Noise is assumed to have a variance 
   depending linearly on the mean. Weights are determined by weigthed distances between 
   color vectors. These color vectors are obtained by demosaicing that is applied to the smoothed
   RAW data after each iteration of the smoothing algorithm. The demosaicing algorithm is 
   a 3D generalized median, see \code{method="Median4"} in function \code{\link{develop.raw}}.
}
\value{
Object of class \code{"adimpro"}
  \item{img }{Contains the reconstructed image. }
  \item{ni }{Contains the sum of weights, i.e. \code{trace(W_i)}, in all grid points \code{i}.}
  \item{ni0 }{Contains the maximum sum of weights for an nonadaptive kernel estimate with the same bandwidth.}
  \item{hmax }{Bandwidth used in the last iteration.}
  \item{call }{The arguments of the function call.}
  \item{varcoef }{Estimated coefficients in the linear variance model for the
    color channels.}
}
\references{   Polzehl, J. and Tabelow, K. (2007). Adaptive smoothing of digital
  images, Journal of Statistical Software 19 (1). }
 \author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\seealso{ \code{\link{read.raw}},\code{\link{awsimage}}, \code{\link{make.image}}, \code{\link{show.image}}, \code{\link{clip.image}} }
\examples{
\dontrun{demo(raw)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
\keyword{ manip }% __ONLY ONE__ keyword per line
