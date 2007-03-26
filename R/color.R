rgb2grey <- function(obj, compress=TRUE) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type!="rgb") {
    warning("Error: image type is not rgb\n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    dim(obj$img) <- c(prod(dm),3)
    obj$img <- as.integer(c(0.299, 0.587, 0.114) %*% t(obj$img))
    dim(obj$img) <- dm
    obj$type <- "greyscale"
  }
  invisible(if(compress) compress.image(obj) else obj)
}

rgb2yiq <- function(obj) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type!="rgb") {
    warning("Error: image type is not rgb \n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    obj$img <- obj$img / 65535
    dim(obj$img) <- c(prod(dm),3)
    
    conv <- c(0.299,  0.595716,  0.211456,
              0.587, -0.274453, -0.522591,
              0.114, -0.321263,  0.311135)
    dim(conv) <- c(3,3)
    
    obj$img <- obj$img %*% t(conv)
    dim(obj$img) <- c(dm,3)
    obj$type <- "yiq"
  }
  invisible(obj)
}

yiq2rgb <- function(obj, compress=TRUE) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type!="yiq") {
    warning("Error: image type is not yiq\n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    dim(obj$img) <- c(prod(dm),3)
    
    conv <- c(1,          1,          1,
              0.9562957, -0.2721221, -1.1069890,
              0.6210244, -0.6473806,  1.7046150)
    dim(conv) <- c(3,3)
    
    obj$img <- as.integer(65535 * pmax(0, pmin(1, obj$img %*% t(conv))))
    dim(obj$img) <- c(dm,3)
    obj$type <- "rgb"
  }
  invisible(if(compress) compress.image(obj) else obj)
}

rgb2yuv <- function(obj) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type!="rgb") {
    warning("Error: image type is not rgb \n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    obj$img <- obj$img / 65535
    dim(obj$img) <- c(prod(dm),3)
    
    conv <- c(0.299, -0.147,  0.615,
              0.587, -0.289, -0.515,
              0.114,  0.436, -0.100)
    dim(conv) <- c(3,3)
    
    obj$img <- obj$img %*% t(conv)
    dim(obj$img) <- c(dm,3)
    obj$type <- "yuv"
  }
  invisible(obj)
}

yuv2rgb <- function(obj, compress=TRUE) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type!="yuv") {
    warning("Error: image type is not yuv\n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    dim(obj$img) <- c(prod(dm),3)
    
    conv <- c( 1,         1,         1,
              -0.000039, -0.394610,  2.032000,
              1.139828, -0.580500, -0.000481)
    dim(conv) <- c(3,3)
    
    obj$img <- as.integer(65535 * pmax(0, pmin(1, obj$img %*% t(conv))))
    dim(obj$img) <- c(dm,3)
    obj$type <- "rgb"
  }
  invisible(if(compress) compress.image(obj) else obj)
}

hsi2rgb <- function (obj, compress=TRUE) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type!="hsi") {
    warning("Error: image type is not hsi\n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    proddm <- prod(dm)
    dim(obj$img) <- c(proddm, 3)
    
    obj$img[, 1] <- obj$img[, 1] * 2 * pi
    obj$img[obj$img[,1]<0,1] <- 0
    obj$img[obj$img[,1]>2*pi,1] <- 2*pi
    
    r <- rep(0, proddm)
    g <- rep(0, proddm)
    b <- rep(0, proddm)
    
    ind1<-(1:proddm)[(obj$img[,1] < 2 * pi/3)]
    ind2<-(1:proddm)[(obj$img[,1] >= 2 * pi/3) & (obj$img[,1] < 4 * pi/3)] # use & operator here to compare vectors!! not &&
    ind3<-(1:proddm)[(obj$img[,1] >= 4 * pi/3)]
    
    b[ind1] <- obj$img[ind1,3] * (1 - obj$img[ind1,2])
    r[ind1] <- obj$img[ind1,3] * (1 + obj$img[ind1,2] * cos(obj$img[ind1,1])/cos(pi/3 - obj$img[ind1,1]))
    g[ind1] <- 3 * obj$img[ind1,3] - (b[ind1] + r[ind1])
    
    r[ind2] <- obj$img[ind2,3] * (1 - obj$img[ind2,2])
    g[ind2] <- obj$img[ind2,3] * (1 + obj$img[ind2,2] * cos(obj$img[ind2,1] - 2 * pi/3)/cos(pi - obj$img[ind2,1]))
    b[ind2] <- 3 * obj$img[ind2,3] - (r[ind2] + g[ind2])
    
    g[ind3] <- obj$img[ind3,3] * (1 - obj$img[ind3,2])
    b[ind3] <- obj$img[ind3,3] * (1 + obj$img[ind3,2] * cos(obj$img[ind3,1] - 4 * pi/3)/cos(5 * pi/3 - obj$img[ind3,1]))
    r[ind3] <- 3 * obj$img[ind3,3] - (g[ind3] + b[ind3])
    
    obj$img <- as.integer(65535 * pmax(0, pmin(1, c(r, g, b))))
    dim(obj$img) <- c(dm, 3)
    obj$type <- "rgb"
  }
  invisible(if(compress) compress.image(obj) else obj)
}

rgb2hsi <- function (obj) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type!="rgb") {
    warning("Error: image type is not rgb\n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    obj$img <- obj$img / 65535
    dim(obj$img) <- c(prod(dm), 3)
    
    z <- obj$img[,1] - 0.5 * obj$img[,2] - 0.5 * obj$img[,3]
    n <- (obj$img[,1] - obj$img[,2])^2 + (obj$img[,1] - obj$img[,3]) * (obj$img[,2] - obj$img[,3])
    frac <- z[n != 0]/sqrt(n[n != 0])
    frac <- signif(frac, 5)
    
    h <- rep(0, prod(dm))
    h[n != 0] <- acos(frac) # h not defined for n==0
    h[obj$img[, 3] > obj$img[,2]] <- 2 * pi - h[obj$img[,3] > obj$img[,2]]
    h <- h/(2 * pi)
    
    i <- (obj$img[,1] + obj$img[,2] + obj$img[,3])/3
    
    s <- rep(0, prod(dm))
    s[i != 0] <- 1 - 1/i[i != 0] * pmin(obj$img[,1], obj$img[,2], obj$img[,3])[i != 0] # s not defined for i==0
    
    obj$img <- c(h, s, i)
    dim(obj$img) <- c(dm, 3)
    obj$type <- "hsi"
  }
  invisible(obj)
}

rgb2xyz <- function(obj) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type!="rgb") {
    warning("Error: image type is not rgb\n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    obj$img <- obj$img / 65535
    dim(obj$img) <- c(prod(dm),3)
    
    conv <- c(0.49,  0.17697, 0.0,
              0.31,  0.81240, 0.01,
              0.20,  0.01063, 0.99)
    dim(conv) <- c(3,3)
    conv <- conv/0.17697
    
    obj$img <- obj$img %*% t(conv)
    dim(obj$img) <- c(dm,3)
    obj$type <- "xyz"
  }
  invisible(obj)
}

xyz2rgb <- function(obj, compress=TRUE) {
  if(!check.adimpro(obj)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (obj$type!="xyz") {
    warning("Error: image type is not xyz\n")
  } else {
    if(obj$compressed) obj <- decompress.image(obj)
    dm <- obj$dim
    dim(obj$img) <- c(prod(dm),3)
    
    conv <- c( 0.41846571, -0.09116896,  0.0009208986,
              -0.15866078,  0.25243144, -0.0025498125,
              -0.08283493,  0.01570752,  0.1785989139)
    dim(conv) <- c(3,3)
    
    obj$img <- as.integer(65535 * pmax(0, pmin(1, obj$img %*% t(conv))))
    dim(obj$img) <- c(dm,3)
    obj$type <- "rgb"
  }
  invisible(if(compress) compress.image(obj) else obj)
}

gamma.correction <- function (img, ga = 2.4, bp = 0.00304,
                              nbins = 65536, alg = 1, log = FALSE) {
  sls <- 1/(ga/bp^(1/ga - 1) - ga * bp + bp)
  fs <- ga/bp^(1/ga - 1)    # slope divided by sls to easy computation
  c0 <- fs * bp^(1/ga) - bp # segment offset divided by sls to easy computation
  if (log) {
    cat("Gamma correction - Gamma:",ga,"\n")
    cat("             Break Point:",bp,"\n")
    cat("                   Slope:",sls,"\n")
    cat("   Slope matching factor:",sls * fs,"\n")
    cat("          Segment offset:",sls * c0,"\n")
  }
  
  img$img <- img$img / 65535
  di <- dim(img$img)
  
  breaks <- seq(0,1+1/(nbins-1),length=nbins+1)
  nimg <- (nbins-1)*img$img
  iimg <- as.integer(nimg)
  ind <- breaks>bp
  if (alg == 1) {
    midpoints <- (breaks[-1]+breaks[-nbins-1])/2
    gammaofmidpoints <- sls * midpoints
    gammaofmidpoints[ind] <- sls * (fs * midpoints[ind]^(1/ga) - c0)
    img$img <- gammaofmidpoints[iimg+1]
  } else if (alg == 2) {
    aimg <- nimg - iimg
    gammaofbreaks <- sls * breaks
    gammaofbreaks[ind] <- sls * (fs * breaks[ind]^(1/ga) - c0)
    img$img <- (1-aimg)*gammaofbreaks[iimg+1] + aimg*gammaofbreaks[iimg+2]
  } else {
    ind <- (1:length(img$img))[img$img > bp]
    img$img[ind] <- fs * img$img[ind]^(1/ga) - c0
    img$img <- sls * img$img
  }
  
  img$img <- as.integer(65535 * img$img)
if(any(img$img>65535)) cat("large values in gamma\n")
  dim(img$img) <- di
  if(ga!=1) img$gamma <- TRUE
  
  invisible(img)
}


white.balance <- function(img, red = 1.0, blue = 1.0, brightness =
                          1.0, log = FALSE) {
  if (log) cat("White balance factors - R:",brightness * red, "G:",brightness,"B:",brightness * blue,"\n")
  max.value <- as.integer(65535)
  if(img$type == "rgb") {
    img$img[,,1] <- as.integer(img$img[,,1] * brightness * red )
    img$img[,,2] <- as.integer(img$img[,,2] * brightness)
    img$img[,,3] <- as.integer(img$img[,,3] * brightness * blue)
    
    img$img[,,1][img$img[,,1]>max.value] <- max.value
    img$img[,,2][img$img[,,2]>max.value] <- max.value
    img$img[,,3][img$img[,,3]>max.value] <- max.value
  } else if (img$type == "greyscale") {
    img$img <- as.integer(img$img * brightness)
    img$img[img$img>max.value] <- max.value
  } else {
    warning("cannot proceed color space",img$type)
  }
  
  invisible(img)
}

adjust.image <- function(img, color.par = NULL, alg = 1, compress=TRUE) {
  
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if(img$compressed) img <- decompress.image(img)
  cp <- colorpar(color.par)
  
  # white balance
  if ( (cp[3] != 1.0) || (cp[4] != 1.0) || (cp[5] != 1.0)) {
    img <- white.balance(img, red = cp[3], blue = cp[4], brightness = cp[5])
  }
  
  # gamma correction
  if (!img$gamma) {
    img <- gamma.correction(img, ga = cp[1], bp = cp[2], alg = alg)
  }
  
  invisible(if(compress) compress.image(img) else img)
}

colorpar <- function(color.par = NULL) {
  ga <- if (is.null(color.par$ga)) 2.4 else color.par$ga
  bp <- if (is.null(color.par$bp)) 0.00304 else color.par$bp
  red <- if (is.null(color.par$red)) 1. else color.par$red
  blue <- if (is.null(color.par$blue)) 1. else color.par$blue
  brightness <- if (is.null(color.par$brightness)) 1. else color.par$brightness
  c(ga, bp, red, blue, brightness)
}
