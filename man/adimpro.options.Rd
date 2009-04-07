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
adimpro.options(xsize = NULL, ysize = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xsize}{ display width in pt }
  \item{ysize}{ display height in pt }
}
\details{
   The function assignes the specified values to as a list to the variable name
   ".adimpro". This variable is, if it exists, evaluated by several other functions.

   On some systems the default \code{X11.options()$type} ''cairo'' leads to significant slower 
   image display. You may try to use \code{X11.options(type=''Xlib'')} instead. To automatically 
   choose this option set the system environment variable \code{R_X11type} (\code{setenv R_X11type Xlib} 
   or \code{export R_X11type=Xlib}) before loading the package.
}
\value{ returns \code{invisible(NULL)} }
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de} }
\keyword{ utilities }
