\name{awsimage}
\alias{awsimage}
\alias{awspimage}
\alias{awsaniso}
\title{Propagation-Separation approach for smoothing of 2D images}
\description{This functions implement the Propagation-Separation approach 
(local constant and local polynomial model) for smoothing images.
Function \code{awsaniso} uses anisotropic location weights. This is done by evaluating local gradient estimates obtained from the actual estimated color values. 
}
\usage{
awsimage(object, hmax=4, aws=TRUE, varmodel=NULL, ladjust=1.25,
         mask=NULL, xind = NULL, yind = NULL,
         wghts=c(1,1,1,1), scorr=TRUE,
         lkern="Plateau", plateau=NULL, homogen=TRUE, earlystop=TRUE,
         demo=FALSE, graph=FALSE,
         max.pixel=4.e2, clip = FALSE, compress=TRUE)
awspimage(object, hmax=12, aws=TRUE, degree=1, varmodel = NULL,
          ladjust=1.0, xind = NULL, yind = NULL, 
          wghts=c(1,1,1,1), scorr= TRUE,
          lkern="Plateau", plateau=NULL, homogen=TRUE, earlystop=TRUE,
          demo=FALSE, graph=FALSE, 
          max.pixel= 4.e2, clip = FALSE, compress=TRUE)
awsaniso(object, hmax = 4, g = 3, rho = 0, aws = TRUE, varmodel = NULL,
          ladjust = 1, xind = NULL, yind = NULL, wghts = c(1, 1, 1, 1), 
          scorr = TRUE, lkern = "Triangle", demo = FALSE, graph = FALSE,
          satexp = 0.25, max.pixel = 400, clip = FALSE, compress = TRUE)
}
\arguments{
  \item{object}{Image object, class "adimpro", as from
    \code{read.image}, \code{read.raw}, or \code{make.image}.}
  \item{hmax}{ Maximum bandwidth to use in the iteration procedure.}
  \item{g}{ Bandwidth for anisotropic smoothing gradient estimates,
    preferably \eqn{g \ge 3}{g >= 3} for images with line type texture and small 
    \eqn{g \approx 1}{g approx 1} for improving edges between homogeneous regions (function \code{awsaniso} only).}
  \item{rho}{ Regularization parameter for anisotropic smoothing gradient estimates,
    preferably \eqn{\rho = 0}{rho = 0} for images with line type texture and large 
    \eqn{\rho \approx 3}{rho approx 3} for improving edges between homogeneous regions. (function \code{awsaniso} only)}
  \item{aws}{ (logical). If \code{TRUE} the propagation - separation
    (PS) approach from Polzehl and Spokoiny (2006) is used.
    \code{aws=FALSE} turns off the statistical penalty resulting in a
    nonadaptive kernel estimate using a kernel with bandwidth \code{hmax}.}
  \item{degree}{ Degree of the local polynomial model for
    \code{awspimage}. 0, 1, or 2 only.}
  \item{varmodel}{ \code{varmodel} specifies how variances are to be
    estimated. This can be a homogeneous variance estimate
    (\code{varmodel="None"}) assuming uncorrelated errors (both spatial
    and between channels). Alternatives are an adaptive homogeneous or linear 
    (function of the mean) variance estimate that depends on estimated
    correlations and on residuals from the last iteration step. 
    The default \code{varmodel=NULL} corresponds to
    \code{varmodel == "Linear"} if \code{img$gamma==FALSE} and
    \code{varmodel == "Constant"} otherwise.}
  \item{ladjust}{ adjustment factor for lambda (>=1). Default values for
    lambda are selected for Gaussian distributions and default settings of 
    parameters \code{lkern} and \code{plateau}. Skewed or heavy
    tailed distributions may require slightly larger values for lambda
    to meet the propagation condition. \code{ladjust} allows to increase
    lambda in such situations. }
  \item{mask}{ logical array of the same size as the image or
    \code{NULL} (default). Smoothing is restricted to the smallest rectangle
    including all pixel where \code{mask==TRUE} and restricts
    computations to these pixel. This need not be a
    connected area (Typical usage: smooth all bright regions)!
    Only used if \code{is.null(xind)} \code{&&}
    \code{is.null(yind)}. Inactive if 
    \code{mask==NULL}. Can only be used if \code{varmodel="None"}. }
  \item{xind, yind}{Restrict smoothing to rectangular area defined by pixel
    indices \code{xind,yind} in x- and y-direction. Full range
    if \code{NULL} (default).}
  \item{wghts}{ allows to weight the information from
    different (up to 4) color channels. The weights are used in the
    statistical penalty of the PS-procedure. Note that lambda-values are selected for \code{wghts==c(1,1,1,1)},
    please use parameter \code{ladjust} to set an appropriate value.}
  \item{scorr}{ (logical). Specifies whether  spatial correlation is to be
    estimated. Defaults to \code{TRUE}.  Is set to \code{FALSE} if
    \code{mask} is not \code{NULL}. }
  \item{lkern}{ Specifies the location kernel. Defaults to
    "Triangle", other choices are "Quadratic", "Cubic" and "Uniform".
    The use of "Triangle" corresponds to the Epanechnicov kernel 
    nonparametric kernel regression. }
  \item{plateau}{Extension of the plateau in the statistical kernel. Can take
    values from (0,1), defaults to \code{0.25}. }
  \item{homogen}{If TRUE the algorithm determines, in each design point i, a circle of maximum radius,
    such that the statistical penalty \code{s_{ij}} for all points j within the
    circle is less than the value specified in \code{plateau}. In subsequent 
    iteration steps the statistical penalty  for such points is set to zero. 
    This is only used if \code{plateau>0}. This results in more stable intermediate estimates and in a smoother reconstruction. \code{homogen=TRUE}
    leads to increased memory requirements.}
  \item{earlystop}{If TRUE the algorithm determines, in each design point i, a circle of minimal radius,
    such that the circle includes all point j with positive weights \code{w_{ij}}.
    if this radius is considerably smaller than the actual bandwidth then the
    estimate in point i is fixed. This should considarably reduce computing time 
    in case of large \code{hmax}.\code{earlystop=TRUE}
    slightly increases memory requirements. }
  \item{demo}{ (logical). If \code{demo=TRUE} the function pauses after each
    iteration. Defaults to \code{FALSE}.  }
  \item{graph}{ (logical). If  \code{graph=TRUE} intermediate results are
    illustrated after each iteration step. Defaults to \code{FALSE}. }
  \item{max.pixel}{ Maximum dimension of images for display
    if \code{graph=TRUE}. If the true dimension is larger, the 
    images are downscaled for display. See also \code{\link{show.image}}.}
  \item{satexp}{exponent used for scaling saturation in anisotropy visualization (function \code{awsaniso} only)}
  \item{clip}{ (logical). If \code{TRUE} a clipping region is selected, see
    \code{\link{clip.image}}, using the information contained in
    \code{xind} or \code{yind}. If both are \code{NULL} a clipping
    region can be defined by left mouse clicks. The image \code{object}
    is reduced to the clipping region before smoothing.}
  \item{compress}{ logical, determines if image data are stored in raw-format. }
}
\details{The function implements the Propagation-Separation (PS) approach to 
  nonparametric smoothing (formerly introduced as Adaptive Weights Smoothing) 
  for varying coefficient likelihood (awsimage) and local polynomial
  (awspimage) models for greyscale and color images. 

  The distribution of grey (color) values is considered to be
  Gaussian. Noise can be colored.

  The numerical complexity of the procedure is mainly determined by
  \code{hmax}. The number of iterations is \code{2*log(hmax)/log(1.25)}.
  Comlexity in each iteration step is \code{Const*hakt*n} with \code{hakt}
  being the actual bandwith in the iteration step and \code{n} the number of pixels. 
  \code{hmax} determines the maximal possible variance reduction.

  All other parameters of the approach only depend on the specified
  values for \code{skern/lkern} and are therefore set internally to
  meaningful default values.

  For a detailed description of the procedure see references below.

  The script used to control the values of parameter \code{lambda} is stored in
  directory inst/adjust.
}
\value{ Object of class \code{"adimpro"}
  \item{img }{Contains the reconstructed image. }
  \item{ni }{Contains the sum of weights, i.e. \code{trace(W_i)}, in all grid points \code{i}.}
  \item{ni0 }{Contains the maximum sum of weights for an nonadaptive kernel estimate with the same bandwidth.}
  \item{hmax }{Bandwidth used in the last iteration.}
  \item{call }{The arguments of the function call.}
  \item{varcoef }{Estimated coefficients in the variance model for the
    color channels, if \code{varmodel} is "Constant" or "Linear". }
  \item{wghts }{The weights used for the color channels.}
  \item{scorr}{Estimated spatial correlations for each channel, if \code{scorr=TRUE} }
  \item{chcorr}{Estimated correlations between
    color channels, if \code{scorr=TRUE} }
}
\references{
  Polzehl and Spokoiny (2006). Propagation-Separation Approach for Local
  Likelihood Estimation. Probability Theory and Related Fields. 3 (135)
  335 - 362.
  
  Polzehl and Spokoiny (2005). Structural adaptive smoothing
  adaptive smoothing by
  Propagation-Separation-methods. WIAS-Preprint No. 1068.
  
  Polzehl, J. and Tabelow, K. (2007). Adaptive smoothing of digital
  images, Journal of Statistical Software 19 (1). 
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\seealso{ \code{\link{read.image}}, \code{\link{read.raw}}, \code{\link{make.image}}, \code{\link{show.image}}, \code{\link{clip.image}} }
\examples{
\dontrun{demo(awsimage)}
}
\keyword{ manip }% at least one, from doc/KEYWORDS
\keyword{ smooth }% __ONLY ONE__ keyword per line
