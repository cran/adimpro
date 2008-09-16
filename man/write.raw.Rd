\name{write.raw}
\alias{write.raw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write image RAW data as greyvalue png image }
\description{
  Image RAW data is saved as a 16-Bit greyvalue png image.
EXIF information contained in the original RAW image is stored as 
a comment.
}
\usage{
write.raw(img, filename = "tmp.png")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{ object of class \code{"adimpro"} containing 
   image RAW data (\code{img$type=="RAW"}.) }
  \item{filename}{ Name of the resulting png-image. If filename does
   not include an extension ".png" the extension ".png" is added.   }
}
\details{ EXIF information contained in the original RAW image as well 
as other available information in object \code{img} are added as a 
comment to the resulting png-image. This comment is evaluated when the
image is read by functions \code{\link{read.raw}} or \code{\link{read.image}} 
}
\value{
  Nothing is returned.}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\note{ ImageMagick has to be installed on the system to write 
  "png" files.}

\seealso{ \code{\link{read.raw}}, \code{\link{develop.raw}} }
\examples{
\dontrun{demo(raw)}
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
