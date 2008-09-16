\name{adimpro.options}
\alias{adimpro.options}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Set parameters for graphical display}
\description{
  On systems capabale of X11 the function sets the X11-type
  (preferably "Xlib"). It also sets a default size for graphical displays
  opended by functions from the package. 
}
\usage{
adimpro.options(xsize = 1280, ysize = 1024, X11.type = "Xlib")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xsize}{ display width in pt }
  \item{ysize}{ display height in pt }
  \item{X11.type}{ sets X11.options(type=X11.type if the system was build with X11 capabilities. We strongly suggest to use \code{X11.type = "Xlib"} which turns out to be
up to 30 times faster than the default X11.options(type="cairo") of R. }
}
\details{
   The function assignes the specified values to as a list to the varaiable name
   ".adimpro". This variable is, if it exists, evaluated by several other functions.
}
\value{ returns \code{invisible(NULL)} }
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de} }
\keyword{ utilities }
