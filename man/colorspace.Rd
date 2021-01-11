\name{colorspace}
\alias{rgb2grey}
\alias{rgb2hsi}
\alias{rgb2yuv}
\alias{rgb2yiq}
\alias{rgb2xyz}
\alias{hsi2rgb}
\alias{yuv2rgb}
\alias{yiq2rgb}
\alias{xyz2rgb}
\title{ Color Space Conversion }
\description{
 Color space conversion functions for RGB, YUV, YIQ, XYZ, and HSI.
}
\usage{
rgb2grey(obj, compress=TRUE)
rgb2hsi(obj)
hsi2rgb(obj, cspace = "Adobe", compress=TRUE)
rgb2yuv(obj)
yuv2rgb(obj, cspace = "Adobe", compress=TRUE)
rgb2yiq(obj)
yiq2rgb(obj, cspace = "Adobe", compress=TRUE)
rgb2xyz(obj)
xyz2rgb(obj, cspace = "Adobe", black= 0, exposure=1, compress=TRUE)
}
\arguments{
  \item{obj}{an object of class "adimpro", \code{obj$type} should coincide with 
    the first three letters of the function name, that is \code{obj$img} containes 
    color values in the corresponding color space.}
  \item{cspace}{Target color space, one of \code{c("sRGB","Adobe","wGamut","kodak","xyz")}.
  For function \code{xyz2rgb} also \code{c("yuv","yiq")} are implemented.}
  \item{exposure}{Multiplicative factor for all color channels (in xyz or rgb spaces). Applied
  in linear color space, i.e. if the image is gamma corrected the gamma correction is reversed first. }
  \item{black}{Adjustment for black color. Color values with luminance <= black will be assigned to
   black in RGB.}
  \item{compress}{ logical, determines if image data are returned in raw-format. }
}
\value{
  an object of class "adimpro", with \code{value$type} specifying the color space
  (last three letters of the function name or 'greyscale' for \code{rgb2grey})
  and \code{value$img} containing the color values. 
 }
\details{
  The functions convert an image \code{obj} from one color space into
  another. 

  \code{rgb2grey} converts from RGB to GREYSCALE.

  \code{rgb2hsi} and \code{hsi2rgb} convert from RGB to HSI and vice
  versa.
  
  \code{rgb2yuv} and \code{yuv2rgb} convert from RGB to YUV and vice
  versa.
  
  \code{rgb2yiq} and \code{yiq2rgb} convert from RGB to YIQ and vice
  versa.

  \code{rgb2xyz} and \code{xyz2rgb} convert from RGB to CIE XYZ and vice
  versa.

  Conversion to \code{XYZ}, \code{YIQ}, \code{YUV} and \code{HSI}
  involves an inverse gamma correction if required.
}
\references{
  Gonzalez, R.C., and Woods, R.E. (2002). Digital Image
  Processing. Prentice Hall. 
  
  Polzehl, J., and Tabelow, K. (2007). Adaptive smoothing of digital
  images, Journal of Statistical Software 19 (1). 
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\examples{
\dontrun{demo(color)}
}
\keyword{color}
