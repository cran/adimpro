read.image <- function(filename, compress=TRUE, convert.path = "convert") {
  fileparts <- strsplit(filename,"\\.")[[1]]
  ext <- tolower(fileparts[length(fileparts)])
  name <- filename
  
  if (ext %in% c("tif","tiff","pgm","ppm","png","pnm","gif","jpg","jpeg")) {
    if (ext %in% c("tif","tiff","png","gif","jpg","jpeg")) {
      tmpfile <- paste(c(fileparts[-length(fileparts)],"ppm"),collapse=".")
      if (file.exists(filename)) {
        if (.Platform$OS.type == "windows") {
          a <- system(paste(convert.path,"-version"),FALSE)
          if (a >= 0) {
            system(paste(convert.path,"-compress None",filename,tmpfile),wait=TRUE)
            filename <- tmpfile            
          } else {
            stop(paste("Error: need ImageMagick (convert). install from www.imagemagick.org\n"))
          }
        } else {
          if (.Platform$OS.type != "unix") warning("never tested this OS. maybe we cannot proceed here.\n")
          a <- system("convert -version",TRUE,TRUE)
          if (length(grep("ImageMagick",a,ignore.case=TRUE)) > 0) {
            system(paste("convert -compress None",filename,tmpfile))
            filename <- tmpfile
          } else {
            stop(paste("Error: need ImageMagick (convert). install from www.imagemagick.org\n"))
          }          
        }
      } else {
        stop(paste("Error: file",filename,"does not exist!"))
      }
    }
    
    object <- list()
    
    if (ext == "pgm") {
      if (file.exists(filename)) {
        object$img <- read.pgm(filename)
      } else {
        stop("cannot find ",filename,"\n")
      }
      object$type <- "greyscale"
    } else {
      if (file.exists(filename)) {
        object$img <- read.ppm(filename)
      } else {
        stop("cannot find ",filename,"\n")
      }
      object$type <- "rgb"
      object$cspace <- "sRGB"
    }
    object$depth <- attr(object$img, "depth")
    attr(object$img, "depth") <- NULL
    object$dim <- dim(object$img)[1:2]
    
    if (ext %in% c("tif","tiff","png","gif","jpg","jpeg")) file.remove(tmpfile)
    
    # in case of length(dim(img))==3  test if image contains same information in all channels
    if(length(dim(object$img)) == 3)
      if(sum(abs(range(object$img[,,1]-object$img[,,2]))+abs(range(object$img[,,1]-object$img[,,3])))==0) {
        object$img <- object$img[,,1]
        object$type <- "greyscale"
      }
    if(compress) {
      dim(object$img) <- NULL 
      object$img <- writeBin(as.integer(object$img),raw(),2)
    }
    object$gamma <- TRUE
    object$wb <- "UNKNOWN"
    object$file <- name
    object$compressed <- compress
    
    class(object) <- "adimpro"
    invisible(object)
    
  } else {
    warning(paste("Error: cannot handle file extension",ext))
    invisible(NULL)
  }
}

write.image <- function(img, file="tmp.ppm", max.x = NULL, max.y =NULL, depth=NULL, color.par = NULL, alg = 1, convert.path = "convert") {
  
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  
  if(img$compressed) img <- decompress.image(img)

  if(!is.null(max.x)||!is.null(max.y)) img <- shrink.image(img, xt=max.x, yt=max.y, method="nearest",compress=FALSE)
  cp <- colorpar(color.par)
  
  fileparts <- strsplit(file,"\\.")[[1]]
  ext <- tolower(fileparts[length(fileparts)])
  
  # supported color depths
  if (is.null(depth)) {
    depth <- switch(img$depth,
                    "8bit" = 8,
                    "16bit" = 16,
                    8)
  }

  # image dimension
  dimg <- dim(img$img)

  # convert colorspace if neccessary
  img <- switch(img$type,
                "hsi" = hsi2rgb(img,compress=FALSE),
                "yuv" = yuv2rgb(img,compress=FALSE),
                "yiq" = yiq2rgb(img,compress=FALSE),
                "xyz" = xyz2rgb(img,compress=FALSE),
                img)
  
  # determine file name for PPM intermediate
  fileparts <- strsplit(file,"\\.")[[1]]
  ext <- tolower(fileparts[length(fileparts)])
  if (img$type == "greyscale") {
    tmpfile <- paste(c(fileparts[-length(fileparts)],"pgm"),collapse=".")
  } else if (img$type == "rgb") {
    tmpfile <- paste(c(fileparts[-length(fileparts)],"ppm"),collapse=".")
  }
  
  # white balance
  if ( (cp[3] != 1.0) || (cp[4] != 1.0) || (cp[5] != 1.0)) {
    img <- white.balance(img, red = cp[3], blue = cp[4], brightness = cp[5])
  }

  # gamma correction
  if (!img$gamma) {
    img <- gamma.correction(img, ga = cp[1], bp = cp[2], alg = alg)
  }

  # rotate image appropriately
  pimg <- switch(img$type,
                 "greyscale" = img$img[,dimg[2]:1],
                 "rgb" = img$img[,dimg[2]:1,],
                 NULL)
  
  # now write
  if(img$type != "unknown") {
    if (depth == 8) {
      pimg <- as.integer(pimg / 256)
      dim(pimg) <- dimg
      maximg <- 255
    } else {
      maximg <- 65535
    }
    ptype <- switch(img$type,
                    "greyscale" = "P5",
                    "rgb" = "P6")
    
    con <- file(tmpfile, "wb")
    writeChar(ptype,con,eos=NULL)
    writeBin(charToRaw("\n"),con)
    writeChar(paste(dimg[1],dimg[2]),con,eos=NULL)
    writeBin(charToRaw("\n"),con)
    writeChar(as.character(maximg),con,eos=NULL)
    writeBin(charToRaw("\n"),con)
    close(con)
    
    if (img$type == "rgb") pimg <- aperm(pimg,c(3,1,2))
    
    if (depth == 8) {
      con <- file(tmpfile, "ab")
      writeBin(as.vector(pimg), con, 1)
      close(con)      
    } else {
      con <- file(tmpfile, "ab")
      writeBin(as.vector(pimg), con, 2, endian="big")
      close(con)              
    }
    
    if (tmpfile != file) {
      if (.Platform$OS.type == "windows") {
        a <- system(paste(convert.path,"-version"),FALSE)
        if (a >= 0) {
          system(paste(convert.path,"-compress None",tmpfile,file),wait=TRUE)
          file.remove(tmpfile)
        } else {
          warning(paste("could not convert",
      tmpfile,"into",file,"need ImageMagick (convert). install from www.imagemagick.org\n",tmpfile,"is kept\n"))
        }
      } else {
        if (.Platform$OS.type != "unix") warning("never tested this OS. maybe we cannot proceed here.\n")
        a <- system("convert -version",TRUE,TRUE)
        if (length(grep("ImageMagick",a,ignore.case=TRUE)) > 0) {
          system(paste("convert -compress None",tmpfile,file))
          file.remove(tmpfile)
        } else {
          warning(paste("could not convert", tmpfile,"into",file,"need ImageMagick (convert). install from www.imagemagick.org\n",tmpfile,"is kept\n"))
        }          
      }
    }
  } else {
    stop("Error: unknown colorspace type! exiting!")
  }
  invisible(NULL)
}

show.image <- function (img, max.x = 1e+03, max.y =1e+03,
                        color.par = NULL, channel=NULL, alg = 1, new = FALSE, ...) {

  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }

  cp <- colorpar(color.par)
  
  # we need of course!
  di <- dim(img$img)
  
  # convert colorspace if neccessary
  img <- switch(img$type,
                "hsi" = hsi2rgb(img,compress=FALSE),
                "yuv" = yuv2rgb(img,compress=FALSE),
                "yiq" = yiq2rgb(img,compress=FALSE),
                "xyz" = xyz2rgb(img,compress=FALSE),
                img)
  
  if (!is.null(channel)&&img$type=="rgb") {
    if(channel %in% (1:3)){
      if(img$compressed){
        nimg <- length(img$img)%/%3
        img$img <- img$img[(channel-1)*nimg+(1:nimg)]
      } else img$img <- img$img[,,channel]
      img$type <- "greyscale"
    }
  }
  if (new) X11()
    
  # now plot according to image type attribute
  switch(img$type,
         "greyscale" = show.greyscale(img,max.x=max.x,max.y=max.y, ga = cp[1], bp = cp[2], 
           alg = alg, brightness = cp[5],...),
         "rgb" = show.rgb(img,max.x=max.x,max.y=max.y, ga = cp[1], bp = cp[2], alg = alg, 
           red = cp[3], blue = cp[4], brightness = cp[5],...),
         return())
  
}

show.rgb <- function(img, max.x=max.x, max.y=max.y,
                     ga = 2.4, bp = 0.00304, alg = 1, 
                     red = 1.0, blue = 1.0, brightness = 1.0,...) {

  # if max.pixel is exceeded, shrink the image to a displayable size 
  dimg0 <- img$dim
  img <- shrink.image(img,xt=max.x,yt=max.y,method="gap",compress=FALSE)

  # check colorspace
  if(any(is.na(img$img))) cat("NA's in img$img\n")
  img$img[img$img < 0] <- as.integer(0)
  img$img[is.na(img$img)] <- as.integer(0)
  minimg <- min(img$img)
  maximg <- max(img$img)
  
  # should never be executed!!!
  if (maximg > 65535) {
    warning("Found values smaller than 0 or larger than 65535. Please check!")
    img$img <- 65535*(img$img - minimg) / (maximg - minimg)
  }
  
  # end check colorspace
  dimg <- img$dim
  
  # define an image z of same size as img and fill with increasing numbers
  z <- matrix(1:prod(dimg),nrow = dimg[1],ncol = dimg[2])
    
  # white balance
  if ( (red != 1.0) || (blue != 1.0) || (brightness != 1.0)) {
    img <- white.balance(img, red = red, blue = blue, brightness = brightness)
  }
   # gamma correction
  if (!img$gamma) {
    img <- gamma.correction(img, ga = ga, bp = bp, alg = alg)
  }
    
  # define the clor map such that for every pixel the color is set
  # according to the value of z
  #    color <-
  #      rgb(img$img[,,1]/65535,img$img[,,2]/65535,img$img[,,3]/65535)
  # we do _not_ use build-in R function rgb() due to memory usage!
  # define 0 to 255 in hexadecimal notation
  hex <- c(0:9, LETTERS[1:6])
  hex <- paste(hex[(0:255)%/%16+1],hex[(0:255)%%16+1],sep="")
  color <- paste("\#",hex[img$img[,,1]%/%256+1],hex[img$img[,,2]%/%256+1],hex[img$img[,,3]%/%256+1],sep="")
    
  x <- seq(1,dimg0[1],length=dimg[1])
  y <- seq(1,dimg0[2],length=dimg[2])
  # display the image
  image(x, y, z, col = color, asp = 1, xlab="",ylab="", ...)
}

show.greyscale <- function(img,max.x=max.x,max.y=max.y,
                           ga = 2.4, bp = 0.00304, alg = 1, 
                           brightness = 1.0,...) {
  # if max.pixel is exceeded, shrink the image to a displayable size 
  dimg0 <- img$dim
  img <- shrink.image(img,xt=max.x,yt=max.y,method="gap",compress=FALSE)
  dimg <- img$dim
  
  # check colorspace
  img$img[img$img < 0] <- 0
  img$img[is.na(img$img)] <- 0
  minimg <- min(img$img)
  maximg <- max(img$img)
  
  # should never be executed!!!
  if (maximg > 65535) {
    warning("Found values smaller than 0 or larger than 65535. Please check!")
    img$img <- (img$img - minimg) / (maximg - minimg)
  }
  
  # end check colorspace
  # define an image z of same size as img and fill with increasing numbers
  z <- matrix(1:prod(dimg),nrow = dimg[1],ncol = dimg[2])

  # white balance
  if (brightness != 1.0) {
    img <- white.balance(img, red = 1, blue = 1, brightness = brightness)
  }

  # gamma correction
  if (!img$gamma) {
    img <- gamma.correction(img, ga = ga, bp = bp, alg = alg)
  }

  color <- grey(img$img/65535)
    
  x <- seq(1,dimg0[1],length=dimg[1])
  y <- seq(1,dimg0[2],length=dimg[2])

  # display the image
  image(x,y,z,col=color, asp=1, xlab="",ylab="",...)
}

read.ppm <- function(filename) {
  # read header
  con <- file(filename,"r")
  type <- readLines(con,n=1)
  count <- 0
  while (TRUE) {
    nl <- readLines(con,n=1)
    if ((regexpr("^ *\#", nl) != -1) || (regexpr("^ *$", nl) != -1)) {
      count <- count + 1
    } else {
      size <- nl
      break
    }
  }
  size <- as.integer(strsplit(size," ")[[1]])
  sizex <- as.integer(size[1])
  sizey <- as.integer(size[2])
  maxval <- as.integer(readLines(con,n=1))
  close(con)

  # read all
  if (type == "P3") {
    img <- as.integer(scan(filename,"int",sizex*sizey*3,skip=3+count,quiet=TRUE))
    if (maxval < 256) {
      img <- img * 256 # make "16bit"
      dd <- "8bit"
    } else {
      dd <- "16bit"      
    }
  } else if (type == "P6") {
    con <- file(filename,"rb")
    ttt <- readLines(con,n=3) # header again
    if (maxval > 255) {
      # endianess not clear. use quasi-standard
      img <- readBin(con,"int",n=sizex*sizey*3,2,signed=FALSE,endian="big") 
      dd <- "16bit"
    } else {
      img <- readBin(con,"int",n=sizex*sizey*3,1,signed=FALSE)
      img <- img * 256 # make "16bit"
      dd <- "8bit"
    }
    close(con)
  } else {
    stop("Cannot handle PPM type",type,"\n")
  }

  img <- as.integer(img)
  dim(img) <- c(3,sizex,sizey)
  img <- aperm(img,c(2,3,1))[,sizey:1,]

  attr(img,"depth") <- dd
  invisible(img)
}

read.pgm <- function(filename) {
  # read header
  con <- file(filename,"r")
  type <- readLines(con,n=1)
  count <- 0
  while (TRUE) {
    nl <- readLines(con,n=1)
    if ((regexpr("^ *\#", nl) != -1) || (regexpr("^ *$", nl) != -1)) {
      count <- count + 1
    } else {
      size <- nl
      break
    }
  }
  size <- as.integer(strsplit(size," ")[[1]])
  sizex <- as.integer(size[1])
  sizey <- as.integer(size[2])
  maxval <- as.integer(readLines(con,n=1))
  close(con)

  # read all
  if (type == "P2") {
    img <- as.integer(scan(filename,"int",sizex*sizey,skip=3+count,quiet=TRUE))
    if (maxval < 256) {
      img <- img * 256 # make "16bit"
      dd <- "8bit"
    } else {
      dd <- "16bit"
    }
  } else if (type == "P5") {
    con <- file(filename,"rb")
    ttt <- readLines(con,n=3) # header again
    if (maxval > 255) {
      # endianess not clear. use quasi-standard
      img <- readBin(con,"int",n=sizex*sizey*3,2,signed=FALSE,endian="big") 
      dd <- "16bit"
    } else {
      img <- readBin(con,"int",n=sizex*sizey*3,1,signed=FALSE)
      img <- img * 256 # make "16bit"
      dd <- "8bit"
    }
    close(con)
  } else {
    stop("Cannot handle PGM type",type,"\n")
  }
  
  img <- as.integer(img)
  dim(img) <- c(sizex,sizey)
  img <- img[,sizey:1]
  
  attr(img,"depth") <- dd
  invisible(img)
}

read.raw <- function (filename,type="PPM",wb="NONE",cspace="sRGB",interp="Bilinear",rm.ppm=TRUE,compress=TRUE) {
  opt1 <- switch(toupper(type),PPM="-4",RAW="-4 -d",HALFSIZE="-h",INFO="-i -v","-4")
  opt2 <- if (opt1 == "-i -v") NULL else switch(toupper(wb),NONE=NULL,AUTO="-a",CAMERA="-w")
  opt3 <- switch(cspace,RAW="-o 0",sRGB=NULL,Adobe="-o 2",wGamut="-o 3",XYZ="-o 5",NULL)
  opt4 <- switch(interp,Bilinear="-q 0",VNG="-q 2",AHD=NULL,FourC="-f","-q 2")
  #  VNG seems to provide minimal spatial correlation
  system(paste("dcraw", opt1, opt2, opt3, opt4, filename))

  object <- list()
  
  if(opt1 != "-i -v") {
    if(opt1 == "-4 -d") {
      filename <- paste(strsplit(filename,"\\.")[[1]][1],".pgm",sep="") 
      object$img <- read.pgm(filename)
      object$type <- "greyscale"
    } else {
      filename <- paste(strsplit(filename,"\\.")[[1]][1],".ppm",sep="") 
      if (file.exists(filename)) {
        object$img <- read.ppm(filename)
      } else {
        # ImageMagick under WINDOWS contains version of dcraw that
        # creates image.ppm instead of filename.ppm
        filename2 <- "image.ppm"
        if (file.exists(filename2)) {
          object$img <- read.ppm(filename2)
          filename <- filename2
        } else {
          stop("cannot find neither ",filename," nor ",filename2,"\n")
        }
      }
      object$type <- "rgb"
    }
    if(rm.ppm) file.remove(filename)
  }
  
  object$depth <- attr(object$img, "depth")
  attr(object$img, "depth") <- NULL
  object$dim <- dim(object$img)[1:2]
  object$file <- filename
  object$interpolation <- interp
  object$cspace <- cspace
  object$gamma <- FALSE
  object$wb <- wb
  if(compress) {
    dim(object$img) <- NULL 
    object$img <- writeBin(as.integer(object$img),raw(),2)
  }
  object$compressed <- compress
  class(object) <- "adimpro"
  invisible(object)
}

  
plot.adimpro <- function(x, new = FALSE, ...) {
  if(!check.adimpro(x)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if(x$compressed) x <- decompress.image(x)
  if ("ni" %in% names(x)) {
    aws <- TRUE
    ni <- make.image(65535* x$ni / x$ni0, compress=FALSE, gamma=FALSE)
  } else {
    aws <- FALSE
  }
  if (x$type == "rgb") {
    col1 <- "red"
    col2 <- "green"
    col3 <- "blue"
    tt <- c("Red channel","Green channel","Blue channel")
    xlim1 <- c(0,65535)
    xlim2 <- c(0,65535)
    xlim3 <- c(0,65535)
  } else if (x$type == "hsi") {
    col1 <- hsv(0:99/99,rep(1,100),rep(1,100))
    col2 <- hsv(rep(1,100),0:99/99,rep(1,100))
    col3 <- grey(0:99/99)
    tt <- c("Hue","Saturation","Intensity")
    xlim1 <- c(0,1)
    xlim2 <- c(0,1)
    xlim3 <- c(0,1)
  } else if (x$type == "yuv") {
    col1 <- grey(0:99/99)
    col2 <- "green1"
    col3 <- "green4"
    tt <- c("Intensity","U channel","V channel")
    xlim1 <- c(0,1)
    xlim2 <- c(-.436,.436)
    xlim3 <- c(-0.615,0.615)
  } else if (x$type == "yiq") {
    col1 <- grey(0:99/99)
    col2 <- "blue1"
    col3 <- "blue4"
    tt <- c("Intensity","I channel","Q channel")
    xlim1 <- c(0,1)
    xlim2 <- c(-.596,.596)
    xlim3 <- c(-.523,.523)
  } else if (x$type == "xyz") {
    col1 <- grey(0:99/99)
    col2 <- "blue1"
    col3 <- "blue4"
    tt <- c("X channel","Y channel","Z channel")
    xlim1 <- NULL
    xlim2 <- NULL
    xlim3 <- NULL
  } else if (x$type == "greyscale") {
    col <- grey(0:99/99)
    tt <- "Grey level"
    xlim <- c(0,65535)
  } else {
    stop("Really dont know what to do with this image type! Sorry!")
  }
  
  if (x$type == "greyscale") {
    if (aws) { 
      if (new) X11(width=5,height=5)
      oldpar <- par(mfrow = c(2,2),mar=c(3,3,3,1))
      on.exit(par(oldpar))
    } else {
      if (new) X11(width=7,height=3)
      oldpar <- par(mfrow = c(1,3),mar=c(3,3,3,1))
      on.exit(par(oldpar))
    }
    hist(x$img,100,col=col,border=col,xlim=xlim,main=tt)
    show.image(x,max.x=400,max.y=400)
    
    plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",frame.plot=FALSE)
    text(0,1,paste("Image file",x$file),pos=4)
    text(0,0.9,paste("Min",min(x$img)),pos=4)
    text(0,0.8,paste("Max",max(x$img)),pos=4)
    
    if (aws) show.image(ni,max.x=400,max.y=400)
    
  } else {
    if (new) X11(width=7,height=5)
    oldpar <- par(mfrow = c(2,3),mar=c(3,3,3,1))
    on.exit(par(oldpar))
    
    hist(x$img[,,1],100,col=col1,border=col1,xlim=xlim1,main=tt[1])
    hist(x$img[,,2],100,col=col2,border=col2,xlim=xlim2,main=tt[2])
    hist(x$img[,,3],100,col=col3,border=col3,xlim=xlim3,main=tt[3])
    show.image(x,max.x=400,max.y=400)
    
    plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",frame.plot=FALSE)
    text(0,1,paste("Image file",x$file),pos=4)
    text(0,0.9,paste("Min",min(x$img)),pos=4)
    text(0,0.8,paste("Max",max(x$img)),pos=4)
    
    if (aws) show.image(ni,max.x=400,max.y=400)
  }

}

summary.adimpro <- function(object, ...) {
  if(!check.adimpro(object)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  #  if(object$compressed) object <- decompress.image(object)
  cat("         Image file:", object$file,"\n")
  if(!is.null(object$xind))
    cat("horizontal clipping:", min(object$xind),":",max(object$xind),"\n")
  if(!is.null(object$yind))
    cat("  vertical clipping:", min(object$yind),":",max(object$yind),"\n")
  cat("    Image dimension:", object$dim,"\n")
  cat("        Color space:", object$type,"\n")
  cat("        Color depth:", object$depth,"\n")
  cat("   Gamma correction:", object$gamma,"\n")
  cat("      White balance:", object$wb,"\n")
  if(!object$compressed) cat("              Range:", as.integer(range(object$img)),"\n")
  if(object$compressed)  cat("   Compressed image\n")
  if (!is.null(object$hmax))
    cat("     max. bandwidth:", signif(object$hmax,2),"\n")
  if (!object$compressed && !is.null(object$ni))
    cat("mean rel adaptation:", signif(mean(object$ni)/object$ni0,2),"\n")
  if (!is.null(object$scorr))
    cat("spatial correlation:", object$scorr,"\n")
  if (!is.null(object$chcorr))
    cat("channel correlation:", object$chcorr,"\n")
  if (!is.null(object$varcoef))
    cat("est. variance param:", object$varcoef,"\n")
}

make.image <- function(x, gamma = FALSE, compress=TRUE){
  dimg <- dim(x)
  if(is.null(dimg) || !(length(dimg) %in% 2:3)) return(warning("x is not an array of appropriate dimensions."))
  if(min(x) < 0) x <- (x-min(x))/(max(x)-min(x))
  if( max(x) <= 1) x <- 65535 * x
  dim(x) <- NULL
  x <- if(compress) writeBin(as.integer(x),raw(),2) else array(as.integer(x),dimg)
  img <- list(img=x, type=switch(length(dimg)-1,"greyscale","rgb"),depth="16bit",
              dim=dimg[1:2], gamma=gamma, wb="UNKNOWN", file="artificial",compressed=compress)
  class(img) <- "adimpro"
  invisible(img) 
}

extract.ni <- function (object, gamma = FALSE, compress=TRUE) {
  if (!check.adimpro(object)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (object$compressed)
    object <- decompress.image(object)
  if ("ni" %in% names(object)) {
    aws <- TRUE
    ni <- make.image(65535 * object$ni/object$ni0, compress = compress, gamma = gamma)
  }
  else {
    stop("image was not processed by awsimage or awspimage")
  }
  invisible(ni)
}
