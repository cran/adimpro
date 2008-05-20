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
if(Sys.getenv("ImageMagick")==""){
cat("Please install ImageMagick\n")
} else {
img <- read.image(system.file("img/wias.ppm",package="adimpro"))
X11(height=5,width=10)
par(mfrow=c(2,4),mar=c(2,2,2,.2),mgp=c(2,1,0))
show.image(img,main="Original Image")
img.edge <- edges(img)
show.image(img.edge,main="Edges in red channel",channel=1)
show.image(img.edge,main="Edges in green channel",channel=2)
show.image(img.edge,main="Edges in blue channel",channel=3)
img.edge <- edges(img)
show.image(img.edge,main="Edges (Laplacian)")
img.edge <- edges(img,abs=TRUE)
show.image(img.edge,main="Edges (abs(Laplacian))")
img.edge <- edges(img,"Sobel")
show.image(img.edge,main="Edges (Sobel)")
img.edge <- edges(img,"Robertcross")
show.image(img.edge,main="Edges (Robertcross)")
}
}
\keyword{ manip }% at least one, from doc/KEYWORDS
