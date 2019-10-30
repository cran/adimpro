\name{mask.create}
\alias{mask.create}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a mask for use within function awsimage}
\description{ Select part of image according to greyscale (or color) value.}
\usage{
mask.create(img,
            range1 = c(0, 1), range2 = c(0, 1), range3 = c(0, 1),
            locate = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{ Image object of class "adimpro", usually the result
    returned from \code{read.image}, \code{read.raw}, or \code{make.image}.  }
  \item{range1}{ Range of gray values or values in the
    first color channel. \code{mask} is set to \code{FALSE} for all
    pixel with values outside this range. Defaults to \code{c(0,1)}.  }
  \item{range2}{ Range of  values in the second color
    channel. \code{mask} is set to \code{FALSE} for all pixel with
    values outside this range. Defaults to \code{c(0,1)}.  }
  \item{range3}{ Range of  values in the third color
    channel. \code{mask} is set to \code{FALSE} for all pixel with
    values outside this range. Defaults to \code{c(0,1)}.  }
  \item{locate}{ (logical). If \code{TRUE} (default), the image is displayed and
    two opposite corners of a rectangular region can be selected using the
    mouse (call of \code{locator(2)}). }
}
\details{ Function to create a mask. The returned object \code{mask} can
  be used to restrict computations in function \code{\link{awsimage}} to
  a region characterized by \code{mask==TRUE}.
}
\value{logical matrix of image dimension. 
}
\references{ Polzehl, J., and Tabelow, K. (2007). Adaptive smoothing of
  digital images, Journal of Statistical Software 19 (1). 
 }
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} and Joerg Polzehl \email{polzehl@wias-berlin.de}}
\seealso{ \code{\link{awsimage}}  }
\keyword{ manip }% at least one, from doc/KEYWORDS
