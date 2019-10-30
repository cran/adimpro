\name{write.image}
\alias{write.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O Functions }
\description{
  Write an image file.
}
\usage{
write.image(img, file = "tmp.ppm", max.x=NULL, max.y=NULL, depth = NULL, gammatype="ITU", 
        whitep = NULL, temp = NULL, cspace = NULL, black = 0, exposure = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{ image data, an object of class "adimpro". }
  \item{file}{ file name, the extension determines the image file format. }
  \item{max.x}{ maximum value of pixels in x dimension. }
  \item{max.y}{ maximum value of pixels in y dimension. }
  \item{depth}{ color depth, either 8 or 16 (bit)! The value
       is reset to 8 if the specified image file format does not allow for 16 Bit encoding. }
  \item{gammatype}{ character, determines the type of gamma correction within the image.
     "ITU" stands for ITU-R BT.709-3 as e.g. used by \code{dcraw}. Alternatives recognized
     within the package are "None", "sRGB" and "CIE" (CIE L*). \code{NULL} keeps the actual setting.
     \code{gammatype="histogram"} forces histogram equalization based on the corresponding greyvalue image.}
  \item{cspace}{ defines the output color space, default is \code{img$cspace} can be set to "sRGB", "Adobe" (Adobe 1998), "wGamut" (Wide Gamut), "kodak" (Kodak ProPhoto)
     "xyz", "yuv", "yiq" and "hsi"}
  \item{whitep}{White point in \code{xyY} space. Can be given as one of (character)
  \code{c("A","B","C","E","D50","D55","D65","D75","F2","F7","F11")} or as a two element numeric
   vector of chromatic \code{xy} coordinates. \code{"D65"} corresponds to the default white point of \code{"sRGB"} and \code{"Adobe"} RGB-spaces. \code{NULL} keeps the actual setting.}
  \item{temp}{Color temperature. Can be used to specify chromatic \code{xy} coordinates of the
   whitepoint. Only used if \code{is.null(whitep)}. }
  \item{black}{Adjustment for black color. Color values with luminance <= black will be assigned to
   black in RGB. Adjustment ist done in xyY space.}
  \item{exposure}{Multiplicative factor for all color channels (in xyz or rgb spaces). Applied
  in linear color space, i.e. if the image is gamma corrected the gamma correction is reversed first. }
}
\details{
  This function writes the image data in \code{img} to the file
  \code{file}. Color depth \code{depth} is used for writing, but if
  image has a component "depth", this argument will be ignored. Note: Not
  all target formats support 16bit coding. The target format is
  determined from the file extension, and should be one of the many that
  ImageMagick supports.
  
  Note that write image by default  applies a gamma correction with \code{gammatype="ITU"}.
  This provides a good standard for display on a screen. For printing 
  \code{cspace="Adobe"} should be prefered. Images that are intended for further editing should preferably be saved in an image format that allows for
  16Bit depth (tiff, png) using either  \code{cspace="Adobe"},  \code{cspace="wGamut"} or \code{cspace="kodak"} and preferably no gamma
  correction to prevent from additional loss of information.
}
\value{
  Nothing is returned.}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\note{ ImageMagick has to be installed on the system to write "tif", "tiff",
  "png", "gif", "jpg" or "jpeg" files.}

\seealso{ \code{\link{read.image}}, \code{\link{adjust.image}} }
\examples{
\dontrun{demo(io)}
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
