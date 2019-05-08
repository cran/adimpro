\name{rimage}
\alias{rimage}
\alias{rimage.options}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Slighty enhanced image function.
}
\description{
The function builds upon function from package graphics, but allows to change some of
it's defaults. These changes can be made through \code{rimage.options}
for subsequent calls of \code{rimage} or by direct specification in \code{\dots}.
}
\usage{
rimage(x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, length.out = ncol(z)), z, ...)
rimage.options(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
  locations of grid lines at which the values in 'z' are
            measured.  These must be finite, non-missing and in
            (strictly) ascending order.  By default, equally spaced
            values from 0 to 1 are used. If 'x' is a 'list', its
          components 'x$x' and 'x$y' are used for 'x' and 'y',
          respectively. If the list has component 'z' this is used for
          'z'. (Same as for function image)

}
  \item{z}{
  a numeric or logical matrix containing the values to be
            plotted ('NA's are allowed).  Note that 'x' can be used
            instead of 'z' for convenience. (Same as for function image)
}
  \item{\dots}{
  The following arguments can be supplied to both \code{rimage} and \code{rimage.options}:
  \itemize{
  \item zquantiles - numeric(2): quantiles of image intensity values to be usd to determine
     ,zlim as \code{zlim <- quantile(z, zquantiles)}. not used if zlim is supplied directly.
     default \code{c(0,1)}.
  \item up - color for intensity values larger than \code{zlim[2]}, default ''white''.
  \item low - color for intensity values smaller than \code{zlim[1]}, default ''black''.
  \item NAcolor - color for intensity NA values, default 0 (transparent).
  \item col - color scheme for values in the range of \code{zlim}, default \code{grey(0:255/255)}.
  \item asp - aspect ratio, default TRUE.
  \item xlab - name for x-axis, default ''x''
  \item ylab - name for y-axis, default ''y''
  \item xaxt - axis type for x-axis, default ''s'', see \link{par}.
  \item yaxt - axis type for y-axis, default ''s'', see \link{par}.
  \item bty - type of box to draw, default ''o'', see \link{par}.
  \item swapx - swap x axis, default FALSE
  \item swapy - swap y axis, default FALSE
  }

    Additionally all parameters that can be passed to function \code{image} via
  \code{\dots} can be used with function \code{rimage}.
}
}
\details{
   This function exists just for convenience to be used if sequences of
   images are to be plotted using the same settings / style. Function
   \code{rimage.options} uses an hidden object \code{.rimage} within an
   evironment \code{.adimproOpts} in the space of package adimpro to store
   the options.

}
\value{
   Both functions return \code{invisible(NULL)}.
}
\author{Joerg Polzehl \email{polzehl@wias-berlin.de}}



\seealso{
 \code{\link{image}}, \code{\link{par}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")
\keyword{ environment }% __ONLY ONE__ keyword per line
