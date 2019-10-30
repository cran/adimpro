\name{segment}
\alias{segment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Segmentation by AWS}
\description{
  The function allows to segment an image into two or three level sets.
}
\usage{
segment(object, level=0.5, delta = 0, thresh = 3, fov = NULL, channel = 0, 
        hmax = 4, aws = TRUE, varmodel = NULL, ladjust = 1.25, xind = NULL, 
        yind = NULL, wghts = c(0.299, 0.587, 0.114, 0), scorr = TRUE, 
        lkern = "Triangle", plateau = NULL, homogen = TRUE, 
        earlystop = TRUE, demo = FALSE, select = FALSE, sext = 1.4, 
        connected = FALSE, graph = FALSE, max.pixel = 400, compress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Image object, class "adimpro", as from
    \code{read.image}, \code{read.raw}, or \code{make.image}.}
  \item{level}{ center of gray/color-values of the second segment, will not be used if \code{select=TRUE}. May be specified such that either 
  \code{level-delta} and \code{level+delta} are within the interval (0,1) or such that they are within the interval (0,65535) (2 Byte integers).}
  \item{delta}{ half width of gray/color-values of the second segment, nay be increased if \code{select=TRUE}. May be specified such that either 
  \code{level-delta} and \code{level+delta} are within the interval (0,1) or such that they are within the interval (0,65535) (2 Byte integers).}
  \item{thresh}{ Critical value for final assignment to segment 1 or 3 , should be specified as a quantile of the standard Gaussian distribution.}
  \item{fov}{ size of field of view in pixel }
  \item{channel}{ specifies which information to use for segmentation. 
   0: use grey valued image obtained from color images, 1-3: use the specified color channel.  }
  \item{hmax}{ Maximum bandwidth to use in the iteration procedure.}
  \item{aws}{ (logical). If \code{TRUE} the propagation - separation
    (PS) approach from Polzehl and Spokoiny (2006) is used.
    \code{aws=FALSE} turns off the statistical penalty resulting in a
    nonadaptive kernel estimate using a kernel with bandwidth \code{hmax}.}
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
    lambda are selected for Gaussian distributions. Skewed or heavy
    tailed distributions may require slightly larger values for lambda
    to meet the propagation condition. \code{ladjust} allows to increase
    lambda in such situations. }
  \item{xind, yind}{Restrict smoothing to rectangular area defined by pixel
    indices \code{xind,yind} in x- and y-direction. Full range
    if \code{NULL} (default).}
  \item{wghts}{ allows to weight the information from
    different (up to 4) color channels. The weights are used in the
    statistical penalty of the PS-procedure. }
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
  \item{select}{ if TRUE a homogeneous rectangular region can be specified interactively. A value of \code{level} is the generated 
    as the median of values within the selected region.}
  \item{sext}{ if \code{select==TRUE} the value of delta is increased by 
   \code{sext} times the standard deviation (estimated by IQR) of the values in the selected region.}
  \item{connected}{ if TRUE the set of pixel within the same segment connected to the specified pixel is extracted. }
  \item{graph}{ (logical). If  \code{graph=TRUE} intermediate results are
    illustrated after each iteration step. Defaults to \code{FALSE}. }
  \item{max.pixel}{ Maximum dimension of images for display
    if \code{graph=TRUE}. If the true dimension is larger, the 
    images are downscaled for display. See also \code{\link{show.image}}.}
  \item{compress}{ logical, determines if image data are stored in raw-format. }
}
\details{
  The image is segmented into three parts by performing multiscale tests 
of the hypotheses H1
\code{value >= level - delta} and H2 \code{value <= level + delta}.
Pixel where the first hypotesis is rejected are classified as \code{-1} (segment 1)
while rejection of H2 results in classification  \code{1} (segment 3).
Pixel where neither H1 or H2 are rejected ar assigned to a value \code{0} (segment 2). Critical values for the tests are adjusted for smoothness at the different scales inspected in the iteration process using results from multiscale testing,
see e.g. Duembgen and Spokoiny (2001). Critical values also depend on the 
size of the region of interest specified in parameter \code{fov}.

Within segment 2 structural adaptive smoothing is performed while if a pair of pixel belongs to segment 1 or segment 3 the corresponding weight will be nonadaptive. 

If \code{connected==TRUE} pixel in segment 2 \code{0} are reassigned to a value \code{2} if they belong to a maximal connected subset of segment2 that contains the center of the specified homogeneous set.
}
\value{ Object of class \code{"adimpro"} with 
  \item{img }{containing a greyvalued image with 3 or 4 levels corresponding to the identified segments.}
and additional list elements
\item{hsegm}{ containing the maximal bandwidth used} 
\item{level}{ the value of parameter \code{level} used}
\item{delta}{ the value of parameter \code{delta} used}
\item{thresh}{ the value of parameter \code{thresh} used}
}
\note{ This function is still experimental and may be changes considerably in future.}
\references{ 
  Duembgen, L. and Spokoiny, V. (2001). Multiscale testing of qualitative hypoteses. Ann. Stat. 29, 124--152.

  Polzehl, J. and Spokoiny, V. (2006). Propagation-Separation Approach for Local
  Likelihood Estimation. Probability Theory and Related Fields. 3 (135)
  335 - 362.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\seealso{ \code{\link{read.image}}, \code{\link{read.raw}}, \code{\link{make.image}}, \code{\link{show.image}}, \code{\link{clip.image}} }
\examples{
\dontrun{demo(segment)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
