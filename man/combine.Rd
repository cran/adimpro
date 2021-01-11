\name{combine}
\alias{combine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pixelwise operations on a pair of images}
\description{
 The function allows to perform pixelwise operations, specified by a
supplied function, on a pair of images. 
}
\usage{
combine(img1, img2, fun = "+", rescale = TRUE, compress = TRUE, gammatype = "None",
        whitep = "D65", cspace = "Adobe", xmode = "RGB", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img1}{ image, object of class adimpro }
  \item{img2}{ image, object of class adimpro, need to have the same dimension as img1 }
  \item{fun}{ A function or primitive of two (or more) arguments
  specifying the operation. The first argument corresponds to 
  grey/color-values in img1, the second to img2. Auxiliary parameters can passed throug \code{"..."} }
  \item{rescale}{ logical: if TRUE the resulting image is rescaled
  to fit into the range of possible grey/color-values, if FALSE
  values outside the range are truncated. }
  \item{compress}{ logical, determines if image data are stored in raw-format. }
  \item{gammatype}{ character, determines the type of gamma correction within the image.
     "ITU" stands for ITU-R BT.709-3 as e.g. used by \code{dcraw}. Alternatives recognized
     within the package are "None", "sRGB" and "CIE" (CIE L*). Please specify if you know
     that your image is not gamma corrected using ITU-R BT.709-3.}
  \item{whitep}{White point in \code{xyY} space. Can be given as one of (character)
  \code{c("A","B","C","E","D50","D55","D65","D75","F2","F7","F11")} or as a two element numeric
   vector of chromatic \code{xy} coordinates. \code{"D65"} corresponds to the default white point of \code{"sRGB"} and \code{"Adobe"} RGB-spaces. }
  \item{cspace}{ defines the output color space, default "sRGB" (sRGB D65), alternatives are 
    "RAW" (Camera specific), "Adobe"
   (Adobe 1998 D65), "wGamut" (Wide Gamut  D65), "kodak" (Kodak ProPhoto D65)
     and "XYZ", see manpages of dcraw.}
  \item{xmode}{\code{xmode} determines how to interpret the values in \code{x} if \code{length(dim(x))==3}.
    Implemented are \code{xmode="RGB"} (default) and \code{xmode="HSI"}}
  \item{\dots}{ additional parameters for function \code{fun} }
}
\details{
There are two mayor applications for this function. First it allows to add noise to an image by first creating an image that contains the noise  and then adding this image using \code{fun="+"}. Second it offers a way to replace parts of an image, see examples.
}
\value{
  object of class "adimpro" containing the image. The object has the following components:
  \item{img}{ array containing the color values in the color space specified by \code{value$type}.}
  \item{type}{ the color space.}
  \item{depth}{ color depth, here "16bit".}
  \item{dim}{ vector of length 2 containing the number of pixel in horizontal and vertival direction.}
  \item{file}{ the argument \code{file} identifying the image.}
  \item{cspace}{ the type of rgb space used, as specified by \code{cspace}.}
  \item{interp}{ interpolation applied by dcraw, as specified by \code{interp}.}
  \item{gamma}{ has a gamma correction been applied, here FALSE for \code{read.raw} and TRUE for \code{read.image}}
  \item{gammatype}{ type of gamma correction \code{read.image.}}
  \item{wb}{ type of white balance, as specified by \code{wb}.}
  \item{compressed}{ image data are stored as raw-vector (TRUE) or array of integers (FALSE).}
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de} }
\seealso{ \code{\link{make.image}} }
\examples{
\dontrun{demo(combine)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
