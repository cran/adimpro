\name{edges}
\alias{edges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Image Processing }
\description{
  Edge detection using Laplacian, Sobel, or Robert Cross filter.
}
\usage{
edges(img, type = "Laplacian", ltype=1, abs=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{ an object of class "adimpro". }
  \item{type}{ type of edges detection filter. "Laplacian" (default), "Sobel" , or "Robertcross". }
  \item{ltype}{ type of laplacian filter. 1,2,3, or 4 }
  \item{abs}{ take absolute values of results. This has only an effect for 
  \code{tyoe="Laplacian"} }
}
\details{
  This function applies the  Laplacian, Sobel, or Robert Cross filter to
  the input image \code{img}. The filter is applied to each color channel separately. \code{ltype} determines
  the different matrices for Laplacian filter used in the
  literature. \code{ltype == 1} will use:

  \code{conv <- matrix(c(-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,
    -1,-1,24,-1,-1,
    -1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1),5,5)}
  
  \code{ltype == 2} will use:

  \code{conv <- matrix(c(0,-1,0,-1,4,-1,0,-1,0), 3, 3)}

  \code{ltype == 3} will use:

  \code{conv <- matrix(c(-1,-1,-1,-1,8,-1,-1,-1,-1), 3, 3)}

  \code{ltype == 4} (default) will use:

  \code{conv <- matrix(c(1,-2,1,-2,4,-2,1,-2,1), 3, 3)}  

}
\value{
  Array containing the values for the edge detector in each pixel and 
  color channel.
  
}
\references{
  Gonzalez, R.C., and Woods, R.E. (2002). Digital Image
  Processing. Prentice Hall. 
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de} }
\examples{
\dontrun{demo(manipulate)}
\dontrun{demo(awspimage)}
}
\keyword{ manip }% at least one, from doc/KEYWORDS
