\name{extract.ni}
\alias{extract.ni}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ extract information about effective size of neighborhoods}
\description{
  The function allows to extract information about the effective
size of neighborhoods used in each pixel from objects 
generated by \code{\link{awsimage}} or \code{\link{awspimage}}. 
The result is converted into a greyscale image. 
}
\usage{
extract.ni(object, gammatype = "ITU", compress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ object returned from \code{\link{awsimage}} or \code{\link{awspimage}}. }
  \item{gammatype}{character, determines the type of gamma correction within the image.
     "ITU" stands for ITU-R BT.709-3 as e.g. used by \code{dcraw}. Alternatives recognized
     within the package are "None", "sRGB" and "CIE" (CIE L*).  }
  \item{compress}{ logical, if \code{TRUE} the returned image will be compressed.}
}
\value{ an object of class "adimpro". 
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\seealso{ \code{\link{awsimage}}, \code{\link{awspimage}}, \code{\link{show.image}}, \code{\link{write.image}} }
\keyword{ manip }% at least one, from doc/KEYWORDS
