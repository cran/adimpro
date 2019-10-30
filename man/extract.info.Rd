\name{extract.info}
\alias{extract.info}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract EXIF information and additional characteristics from an 
object of class "adimpro". }
\description{
  Extract EXIF information and additional characteristics from an 
object of class "adimpro". 
}
\usage{
extract.info(object, what = "Bayer")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of class "adimpro" or a character string that
  was written as a comment by functions \code{\link{write.image}} or 
  \code{\link{write.raw}}. }
  \item{what}{ A character string specifying which information is to be extracted. Currently implemented are \code{"Bayer"} (Bayer mask), \code{"Daymulti"} (Daylight multiplyers), \code{"Cammulti"} (camera multipliers), \code{"Camera"} (camera model), \code{"Isize"} (Image size), \code{"Osize"} (image output size), \code{"File"} (source filename), \code{"Interpolation"}(Interpolation used to convert RAW to RGB), \code{"Gammatype"} (type of gamma correction), \code{"WhiteBalance"} (white balance), \code{"WhitePoint"} (white point), \code{"Type"} (image type), \code{"xind"} (horizontal subindex  used by function \code{\link{clip.image}}) and \code{"yind"} (vertical subindex  used by function \code{\link{clip.image}}).}
}
\details{
  This function is mainly used to access and manipulate information when reading or developing images. A more comprehensive overview is provided by function \code{\link{summary}}.
}
\value{
  A character string or integer vector depending on the argument \code{"what"}.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\note{ ImageMagick has to be installed on the system to write "tif", "tiff",
  "png", "gif", "jpg" or "jpeg" files.}

\seealso{ \code{\link{read.raw}}, \code{\link{read.image}}, \code{\link{develop.raw}}, \code{\link{summary}} }
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
