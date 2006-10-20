\name{write.image}
\alias{write.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O Functions }
\description{
  Write an image file.
}
\usage{
write.image(img, file = "tmp.ppm", max.x=NULL, max.y=NULL, depth = NULL,color.par = NULL,
            alg = 1, convert.path = "convert")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{ image data, an object of class "adimpro". }
  \item{file}{ file name, the extension determines the image file format. }
  \item{max.x}{ maximum value of pixels in x dimension. }
  \item{max.y}{ maximum value of pixels in y dimension. }
  \item{depth}{ color depth, either 8 or 16 (bit)! The value
       is reset to 8 if the specified image file format does not allow for 16 Bit encoding. }
  \item{color.par}{ list of parameters for the gamma correction and some
    white balance. list entries:
    \tabular{rll}{
      ga \tab gamma value \tab default 2.4\cr
      bp  \tab break point \tab default 0.00304\cr
      alg \tab algorithm for gamma correction \tab default 1\cr
      red \tab factor to multiply with red channel. \tab default 1.0\cr
      blue \tab factor to multiply with blue channel. \tab default 1.0\cr
      brightness \tab factor to multiply with all channels \tab default 1.0}
  }  
  \item{alg}{ chooses computing algorithms which differ in speed and precision. }
  \item{convert.path}{ Windows system calls seem to change environment
    variable PATH, such that "convert" points to
    c:\\windows\\system32\\convert.exe instead of
    ImageMagick-x.x.x\\convert. Therefore here the full path must be
    provided. Backslash needs to be escaped: c:\\\\programs\\\\imagemagick-6.3.0\\\\convert }
}
\details{
  This function writes the image data in \code{img} to the file
  \code{file}. Color depth \code{depth} is used for writing, but if
  image has a component "depth", this argument will be ignored. Note: Not
  all target formats support 16bit coding. The target format is
  determined from the file extension, and should be one of the many that
  ImageMagick supports.

  If \code{color.par$red} or \code{color.par$blue} or \code{color.par$brightness} differ from 1.0 the
  corresponding channels are multiplied with the provided
  values. Saturated values are set to 1.

  If \code{img$gamma} is FALSE, perform gamma correction with
  \code{color.par$ga} and \code{color.par$bp}. \code{alg} chooses between three different
  computing algorithms (approximations) with increasing computation time
  and precision (\code{alg} is 1,2, or 3).

}
\value{
  Nothing is returned.}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\note{ ImageMagick has to be installed on the system to write "tif", "tiff",
  "png", "gif", "jpg" or "jpeg" files.
}

\seealso{ \code{\link{read.image}} }
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
