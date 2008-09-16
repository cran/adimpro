\name{show.image}
\alias{show.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O functions }
\description{
  Display an image on the screen.
}
\usage{
show.image(img, max.x = 1000, max.y = 1000, gammatype = "ITU", 
           whitep = NULL, temp = NULL, cspace = "sRGB", black=0, exposure = 1, 
           channel=NULL, new = FALSE,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{ image data, an object of class "adimpro".}
  \item{max.x}{ maximum value of pixels in x dimension to be displayed. }
  \item{max.y}{ maximum value of pixels in y dimension to be displayed. }
  \item{gammatype}{ character, determines the type of gamma correction within the image.
     "ITU" stands for ITU-R BT.709-3 as e.g. used by \code{dcraw}. Alternatives recognized
     within the package are "None", "sRGB" and "CIE" (CIE L*). 
     \code{gammatype="histogram"} forces histogram equalization based on the corresponding greyvalue image.}
  \item{cspace}{ defines the output color space, default "sRGB" (sRGB D65), alternatives are 
     "Adobe" (Adobe 1998 D65), "wGamut" (Wide Gamut  D65), "kodak" (Kodak ProPhoto D65)
     "xyz", "yuv", "yiq" and "hsi". \code{NULL} keeps the actual setting. If \code{color.space\%in\%c("hsi","yuv","yiq","xyz")} the individual channels are rescaled to provide maximum 
    contrast. Information from the three channels is coded as "red", "green" and "blue" providing a 
    miscolored image in these cases. Individual channels can be displayed as greyvalue images by
    specifying the \code{channel} argument. 
    \code{color.space="greyvalue"} provides a greyvalue image.}
  \item{whitep}{White point in \code{xyY} space. Can be given as one of (character)
  \code{c("A","B","C","E","D50","D55","D65","D75","F2","F7","F11")} or as a two element numeric
   vector of chromatic \code{xy} coordinates. \code{"D65"} corresponds to the default white point of \code{"sRGB"} and \code{"Adobe"} RGB-spaces. \code{NULL} keeps the actual setting.}
  \item{temp}{Color temperature. Can be used to specify chromatic \code{xy} coordinates of the
   whitepoint. Only used if \code{is.null(whitep)}. }
  \item{black}{Adjustment for black color. Color values with luminance <= black will be assigned to
   black in RGB. Adjustment ist done in xyY space.}
  \item{exposure}{Multiplicative factor for all color channels (in xyz or rgb spaces). Applied
  in linear color space, i.e. if the image is gamma corrected the gamma correction is reversed first. }
  \item{channel}{ allows to select a color channel (1: red, 2: green, 3: blue in case of "rgb") for display.}
  \item{new}{ should new X11() be opened? default \code{FALSE}}
  \item{\dots}{ additional arguments to \code{\link{image}} can be
    passed here. }
}
\details{
  This function displayes greyscale and color images
  on the screen. If the actual dimension of the
  image exceeds \code{max.x} or \code{max.y} the image is shrinked by
  displaying only part of the pixels (every second/third/... such that
  the resulting dimension is smaller than \code{max.x} or \code{max.y})
  
  

}
\value{
  Nothing is returned
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\seealso{ \code{\link{read.image}}, \code{\link{write.image}}, \code{\link{adjust.image}} }
\examples{
\dontrun{demo(io)}
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
