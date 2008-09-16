\name{plot.adimpro}
\alias{plot.adimpro}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O functions }
\description{
  Visualize image data.
}
\usage{
\method{plot}{adimpro}(x, new=FALSE, gammatype = NULL, cspace = NULL,
                 whitep = NULL, temp = NULL, black = 0, exposure = 1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{image object of class "adimpro" }
  \item{new}{ should new X11() be opened? default \code{FALSE}}
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
  \item{...}{ not used }
}
\details{
  This functions shows information on the image. This includes
  histograms of color values in each channel of the specified (\code{x$type}) color space, a thumbnail (in "sRGB"
  with gamma correction \code{gammatype} if specified and \code{gammatype="ITU"} elsewhere),
  some information on the image and, if x was produced by \code{awsimage},
  an image illustrating the local adaptation.
}
\value{
  nothing is returned.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\seealso{ \code{\link{summary.adimpro}}, \code{\link{adjust.image}} }
\examples{
\dontrun{demo(color)}
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ iplot }% at least one, from doc/KEYWORDS
