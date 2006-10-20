\name{clip.image}
\alias{clip.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create an image by clipping }
\description{
The function allows to define a clipping region by arguments \code{xind} and \code{yind}
or interactively by mouseclicks. A new image is created by cutting out 
the clipping region.
}
\usage{
clip.image(img, xind = NULL, yind = NULL, compress=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{ Object of class "adimpro" as created by
    \code{\link{read.image}}, \code{\link{read.raw}}
    or \code{\link{make.image}}. }
  \item{xind}{ \code{xind} defines the horizontal extension of the clipping region.  }
  \item{yind}{ \code{yind} defines the vertical extension of the clipping region. }
  \item{compress}{image data are stored as raw-vector (TRUE) or array of integers (FALSE). \code{compress=NULL} 
  keeps the format used in img. }
}
\details{
  If both \code{xind==NULL} and \code{yind==NULL}, \code{\link{show.image}} is called and the clipping region
  can be set by two left mouse clicks to opposite corners. 
}
\value{ An object of class "adimpro". }
\author{ Joerg Polzehl \email{polzehl@wias-berlin.de} }
\seealso{ \code{\link{read.image}}, \code{\link{read.raw}}, \code{\link{make.image}}, \code{\link{awsimage}}, \code{\link{awspimage}}}
\keyword{ manip }% at least one, from doc/KEYWORDS
