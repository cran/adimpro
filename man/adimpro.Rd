\name{adimpro}
\alias{read.raw}
\alias{read.image}
\alias{make.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O Functions }
\description{
  Create image objects of class "adimpro" from arrays, RAW-format files
  and other image formats.
}
\usage{
read.raw(filename, type="PPM",
         wb="CAMERA",cspace="Adobe",interp="Bilinear",maxrange=TRUE,
         rm.ppm=TRUE, compress=TRUE)
read.image(filename, compress=TRUE)
make.image(x,compress=TRUE, gammatype="None", whitep = "D65", 
         cspace="Adobe", scale="Original",xmode="RGB")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ file name }
  \item{x}{ Array or matrix containing RGB or greyscale values in the
    range (0,1) or (0,65535). }
  \item{type}{ option settings for dcraw. default "PPM". \code{type="png"} allows to read greyvalue png images as RAW-data (used as internal solution to
  store RAW information) }
  \item{wb}{ white balance. default "CAMERA"}
  \item{cspace}{ defines the output color space, default "sRGB" (sRGB D65), alternatives are 
    "RAW" (Camera specific), "Adobe"
   (Adobe 1998 D65), "wGamut" (Wide Gamut  D65), "kodak" (Kodak ProPhoto D65)
     and "XYZ", see manpages of dcraw.}
  \item{interp}{ defines the interpolation method, default "Bilinear", Alternatives are
   "VNG", "AHD", "FourC" (Four color interpolation) and "Halfsize", see manpages of dcraw. 
   "VNG" seems to provide the smallest spatial correlations.}
  \item{maxrange}{If TRUE increase range of values to maximum.}
  \item{rm.ppm}{ remove intermediate tmp file? default \code{TRUE}}
  \item{gammatype}{ character, determines the type of gamma correction within the image.
     "ITU" stands for ITU-R BT.709-3 as e.g. used by \code{dcraw}. Alternatives recognized
     within the package are "None", "sRGB" and "CIE" (CIE L*). Please specify if you know
     that your image is not gamma corrected using ITU-R BT.709-3.}
  \item{whitep}{White point in \code{xyY} space. Can be given as one of (character)
  \code{c("A","B","C","E","D50","D55","D65","D75","F2","F7","F11")} or as a two element numeric
   vector of chromatic \code{xy} coordinates. \code{"D65"} corresponds to the default white point of \code{"sRGB"} and \code{"Adobe"} RGB-spaces. }
  \item{compress}{ logical, determines if image data are stored in raw-format. }
  \item{scale}{\code{"Original"} scales to \code{(0,max(img$img))} if 
  \code{min(img$img)<0}, otherwise keeps the original scale.
  \code{"Maxcontrast"} scales each channel to maximum contrast}
  \item{xmode}{\code{xmode} determines how to interpret the values in \code{x} if \code{length(dim(x))==3}.
    Implemented are \code{xmode="RGB"} (default) and \code{xmode="HSI"}}
}
\details{
  If ImageMagick is available on the system, \code{read.image} reads any of
  the following image file formats: c("tif", "tiff", "pgm", "ppm",
  "png", "pnm", "gif", "jpg", "jpeg") converts it into a temporary "pgm" or
  "ppm" file. This file is removed after reading the image.  
  If ImageMagick is not available  only "pgm", "ppm" and "pnm" formats
  can be processed.

  If dcraw is available on the system, \code{read.raw} reads many RAW
  formats. \code{type} sets options to dcraw: "PPM" sets "-4", "RAW" sets "-4
  -d", "HALFSIZE" sets "-h", "INFO" sets "-i -v". \code{wb} indicates,
  which white balance should be used: "NONE", "AUTO", "CAMERA".

  Functions \code{\link{read.raw}(file,type="RAW")} and \code{\link{read.image}(file)} provide identical results on png-mages. If the result is a color, greyvalued
  or RAW image depends on the contend of the comment associated with the png-image.

  \code{make.image} converts an appropriate 2 or 3 dimensional array to
  an image object of class "adimpro".
}
\value{
  object of class "adimpro" containing the image. The object has the following components:
  \item{img}{ array containing the color values in the color space specified by \code{value$type}.}
  \item{type}{ the color space.}
  \item{depth}{ color depth, here "16bit".}
  \item{dim}{ vector of length 2 containing the number of pixel in horizontal and vertival direction.}
  \item{file}{ the argument \code{file} identifying the image.}
  \item{cspace}{ the type of rgb space used, as specified by \code{cspace}.}
  \item{interp}{ interpolation applied by dcraw, as specified by \code{interp}.}
  \item{gamma}{ has a gamma correction been applied, here FALSE for \code{read.raw} and TRUE for \code{read.image} }
  \item{gammatype}{ type of gamma correction \code{read.image.}}
  \item{wb}{ type of white balance, as specified by \code{wb}.}
  \item{compressed}{ image data are stored as raw-vector (TRUE) or array of integers (FALSE).}
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de} }
\note{ The function \code{read.raw} requires dcraw to be installed.}
\seealso{ \code{\link{read.image}} }
\examples{
\dontrun{demo(io)}
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
