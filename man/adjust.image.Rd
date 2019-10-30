\name{adjust.image}
\alias{adjust.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Image Processing }
\description{ Color space transformations, change of white balance and exposure,
  gamma correction and histogram equalization.
}
\usage{
adjust.image(img, gammatype=NULL, cspace = NULL, whitep = NULL, 
                 temp = NULL, black=0, exposure = 1, kind = "Bradford", 
                 alg = 1, compress= TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{ image object, class"adimpro". }
  \item{gammatype}{ character, determines the type of gamma correction within the image.
     "ITU" stands for ITU-R BT.709-3 as e.g. used by \code{dcraw}. Alternatives recognized
     within the package are "None", "sRGB" and "CIE" (CIE L*). \code{NULL} keeps the actual setting.
     \code{gammatype="histogram"} forces histogram equalization based on the corresponding greyvalue image.}
  \item{cspace}{ defines the output color space, default "sRGB" (sRGB D65), alternatives are 
     "Adobe" (Adobe 1998 D65), "wGamut" (Wide Gamut  D65), "kodak" (Kodak ProPhoto D65)
     "xyz", "yuv", "yiq" and "hsi". \code{NULL} keeps the actual setting.}
  \item{whitep}{White point in \code{xyY} space. Can be given as one of (character)
  \code{c("A","B","C","E","D50","D55","D65","D75","F2","F7","F11")} or as a two element numeric
   vector of chromatic \code{xy} coordinates. \code{"D65"} corresponds to the default white point of \code{"sRGB"} and \code{"Adobe"} RGB-spaces. \code{NULL} keeps the actual setting.}
  \item{temp}{Color temperature. Can be used to specify chromatic \code{xy} coordinates of the
   whitepoint. Only used if \code{is.null(whitep)}. }
  \item{black}{Adjustment for black color. Color values with luminance <= black will be assigned to
   black in RGB. Adjustment ist done in xyY space.}
  \item{exposure}{Multiplicative factor for all color channels (in xyz or rgb spaces). Applied
  in linear color space, i.e. if the image is gamma corrected the gamma correction is reversed first. }
  \item{kind}{Algorithm for chromatic adaptation. Alternatives are \code{"Bradford"}, \code{"VonKries"} and \code{"XYZscaling"}}
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
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\seealso{ \code{\link{show.image}}, \code{\link{write.image}} }
\examples{
\dontrun{demo(color)}
}
\keyword{ manip }% at least one, from doc/KEYWORDS
