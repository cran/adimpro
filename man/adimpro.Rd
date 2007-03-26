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
         wb="NONE",cspace="sRGB",interp="Bilinear",
         rm.ppm=TRUE, compress=TRUE)
read.image(filename, compress=TRUE, convert.path="convert")
make.image(x,gamma=FALSE,compress=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ file name }
  \item{x}{ Array or matrix containing RGB or greyscale values in the
    range (0,1) or (0,65535). }
  \item{type}{ option settings for dcraw. default "PPM", alternatives are "RAW", "Halfsize" and "Info".}
  \item{wb}{ white balance. default "CAMERA", alternatives are "AUTO" and "NONE".}
  \item{cspace}{ defines the output color space, default "sRGB", alternatives are "RAW", "Adobe", "wGamut"
     and "XYZ", see manpages of dcraw.}
  \item{interp}{ defines the interpolation method, default "Bilinear", Alternatives are
   "VNG", "" and "FourC" (Four color interpolation), see manpages of dcraw. 
   "VNG" seems to provide the smallest spatial correlations.}
  \item{rm.ppm}{ remove intermediate tmp file? default \code{TRUE}}
  \item{gamma}{ logical, determines \code{value$gamma}}
  \item{compress}{ logical, determines if image data are stored in raw-format. }
  \item{convert.path}{ Windows system calls seem to change environment
    variable PATH, such that "convert" points to
    c:\\windows\\system32\\convert.exe instead of
    ImageMagick-x.x.x\\convert. Therefore here the full path must be
    provided. Backslash needs to be escaped: c:\\\\programs\\\\imagemagick-6.3.0\\\\convert }
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
  \item{gamma}{ has a gamma correction been applied, here FALSE.}
  \item{wb}{ type of white balance, as specified by \code{wb}.}.
  \item{compressed}{ image data are stored as raw-vector (TRUE) or array of integers (FALSE).}
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de}, Joerg Polzehl \email{polzehl@wias-berlin.de} }
\note{ The function requires dcraw to be installed.}
\seealso{ \code{\link{read.image}} }
\examples{
\dontrun{read.raw("test.crw")}
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
