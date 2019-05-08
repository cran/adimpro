#
rimage <- function(x = seq(0, 1, length.out = nrow(z)),
                   y = seq(0, 1, length.out = ncol(z)),
                   z, ...){
  impars <- get(".rimage",envir=.adimproOpts)
  args <- list(...)
  nargs <- length(args)
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        if (is.null(dim(x)))
          stop("argument must be matrix-like")
        z <- x
        x <- 1:nrow(z)
        y <- 1:ncol(z)
      }
    }
  }
  zlim <- quantile(z, impars[["zquantiles"]], na.rm=TRUE)
  low <- impars[["low"]]
  up <- impars[["up"]]
  NAcolor <- impars[["NAcolor"]]
  col <- impars[["col"]]
  asp <- impars[["asp"]]
  xlab <- impars[["xlab"]]
  ylab <- impars[["ylab"]]
  xaxt <- impars[["xaxt"]]
  yaxt <- impars[["yaxt"]]
  bty <- impars[["bty"]]
  swapx <- impars[["swapx"]]
  swapy <- impars[["swapy"]]
  if(length(args)>0){
     nargs <- names(args)
     if("zquantiles" %in% nargs){
       zlim <- quantile(z, args[["zquantiles"]], na.rm=TRUE)
       args[["zquantiles"]] <- NULL
     }
     if("zlim" %in% nargs){
        zlim <- args[["zlim"]]
        args[["zlim"]] <- NULL
      }
     if("up" %in% nargs){
        up <- args[["up"]]
        args[["up"]] <- NULL
      }
     if("low" %in% nargs){
        low <- args[["low"]]
        args[["low"]] <- NULL
      }
     if("NAcolor" %in% nargs){
        NAcolor <- args[["NAcolor"]]
        args[["NAcolor"]] <- NULL
      }
     if("col" %in% nargs){
        col <- args[["col"]]
        args[["col"]] <- NULL
      }
     if("asp" %in% nargs){
         asp <- args[["asp"]]
         args[["asp"]] <- NULL
       }
     if("xlab" %in% nargs){
        xlab <- args[["xlab"]]
        args[["xlab"]] <- NULL
      }
     if("ylab" %in% nargs){
        ylab <- args[["ylab"]]
        args[["ylab"]] <- NULL
      }
     if("xaxt" %in% nargs){
         xaxt <- args[["xaxt"]]
         args[["xaxt"]] <- NULL
       }
     if("yaxt" %in% nargs) {
         yaxt <- args[["yaxt"]]
         args[["yaxt"]] <- NULL
       }
     if("bty" %in% nargs) {
         yaxt <- args[["bty"]]
         args[["bty"]] <- NULL
       }
  }
  eps <- diff(zlim)/length(col)
  if(zlim[1] > min(z,na.rm=TRUE)){
     zlim[1] <- zlim[1]-eps
     z[z<zlim[1]] <- zlim[1]
     col <- c(low,col)
  }
  if(zlim[2] < max(z,na.rm=TRUE)){
     zlim[2] <- zlim[2]+eps
     z[z>zlim[2]] <- zlim[2]
     col <- c(col,up)
  }
  if(swapx) {
     z <- z[length(x):1,]
  }
  if(swapy) {
     z <- z[,length(y):1]
  }
  args <- c(list(x=x, y=y, z=z, zlim=zlim, col=col, asp=asp,
    xlab=xlab, ylab=ylab, xaxt=xaxt, yaxt=yaxt, bty=bty),args)
  do.call(image,args)
  if(any(is.na(z))) image(x,y,is.na(z),col=c(NA,NAcolor),add=TRUE)
  invisible(NULL)
}
