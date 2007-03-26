awsimage <- function (object, hmax=4, aws=TRUE, varmodel=NULL,
                      ladjust=1.0 , mask = NULL, xind = NULL, 
                      yind = NULL, wghts=c(1,1,1,1), scorr=TRUE, lkern="Triangle", demo=FALSE, graph=FALSE, max.pixel=4.e2,clip=FALSE,compress=TRUE) {

  #
  #          Auxilary functions
  #
  IQRdiff <- function(data) IQR(diff(data))/1.908

  #####################################################################################
  ###    
  ###    function body
  ###    
  ###    first check arguments and initialize
  ###
  #####################################################################################    
  if(!check.adimpro(object)) {
    cat(" Consistency check for argument object failed (see warnings). object is returned.\"n")
    return(invisible(object)) 
  }
  object <- switch(object$type,
                   "hsi" = hsi2rgb(object),
                   "yuv" = yuv2rgb(object),
                   "yiq" = yiq2rgb(object),
                   "xyz" = xyz2rgb(object),
                   object)
  if(!is.null(object$call)) { 
    warning("argument object is result of awsimage or awspimage. You should know what you do.")
  }
  args <- match.call()
  if(!is.null(mask)) {
    varmodel <- "None"
    scorr <- FALSE
    clip <- FALSE
  }
  if(clip) object <- clip.image(object,xind,yind,compress=FALSE)
  if(object$compressed) object <- decompress.image(object)
  if(is.null(varmodel)) {
    varmodel <- if(object$gamma) "Constant" else "Linear"
  }
  bcf <- c(-.4724,1.1969,.22167,.48691,.13731,-1.0648)
# coefficients for  bias correction for spatial correlation
  estvar <- toupper(varmodel) %in% c("CONSTANT","LINEAR")
  hpre <- 2.
  if(estvar) {
    dlw<-(2*trunc(hpre)+1)
    nvarpar <- switch(varmodel,Constant=1,Linear=2,1)
  }
  #
  #   Check image type
  #
  dimg <- dimg0 <- dim(object$img)
  imgtype <- object$type 
  dv <- switch(imgtype,rgb=dimg[3],greyscale=1)
  if(length(wghts)<dv) wghts <- c(wghts,rep(1,dv-length(wghts))) else wghts <- pmin(.1,wghts) 
  wghts <- switch(imgtype, greyscale=1, rgb = wghts[1:dv])
  spcorr <- matrix(0,2,dv)
  h0 <- c(0,0)
  #
  #   specify the range of data needed if !is.null(mask) || !is.null(xind) || !is.null(yind)
  #
  if( valid.index(xind,dimg[1]) || valid.index(yind,dimg[2]))   mask <- NULL
  if(!is.null(mask)) {
    if(!is.logical(mask)||dim(mask)!=dimg[1:2]) {
      mask<-NULL
      warning("Argument mask changed to NULL")
    } else if(sum(mask)==0) mask<-NULL
  }
  if(!is.null(mask)){
    xmask <- range((1:dimg[1])[apply(mask,1,any)])
    ymask <- range((1:dimg[2])[apply(mask,2,any)])
    xind <- xmask[1]:xmask[2]
    yind <- ymask[1]:ymask[2]
    n1 <- length(xind)
    n2 <- length(yind)
    n <- n1*n2
  } else {
    if(is.null(xind)) xind <- 1:dimg[1]
    if(is.null(yind)) yind <- 1:dimg[2]
    n1 <- length(xind)
    n2 <- length(yind)
    n <- n1*n2
  }
  dimg <- switch(imgtype,greyscale=c(n1,n2),rgb=c(n1,n2,dv))
  #
  #     set approriate defaults
  #
  if (aws) qlambda <- .9 else qlambda <- 1
  ladjust <- max(1,ladjust)
  lseq <- pmax(1,c(1,1,1,1,1,1,1,1.6,1.5,1.4,1.3,1.2,1.1)/ladjust)
  lkern <- switch(lkern,
                  Triangle=2,
                  Quadratic=3,
                  Cubic=4,
                  Uniform=1,
                  2)
  if (is.null(hmax)) hmax <- 4
  wghts <- wghts/sum(wghts)
  dgf <- sum(wghts)^2/sum(wghts^2)
  if (qlambda<1) lambda <- ladjust*2*qchisq(qlambda,dgf)/dgf else lambda <- 1e50
  #
  #      in case of colored noise get the corresponding bandwidth (for Gaussian kernel)
  #
  sigma2 <- switch(imgtype,
                   "greyscale" = IQRdiff(object$img)^2,
                   "rgb" = apply(object$img,3,IQRdiff)^2,
                   IQRdiff(object$img)^2)
  cat("Estimated variance (assuming independence): ", signif(sigma2/65635^2,4),"\n")
  if (!estvar) {
    if (length(sigma2)==1) {
      #   homoskedastic Gaussian case
      if(dv>1) sigma2 <- rep(sigma2,1)
      lambda <- lambda*sigma2 
    } else if (length(sigma2)==dv) {
      wghts <- wghts[1:dv]/sigma2
    } 
  }
    # set the support of the statistical kernel to (0,1), set spmin and spmax
    lambda <- 2*lambda
    spmin <- 0
    spmax <- 1
  #     now set hinit and hincr if not provided
  if(aws) hinit <- 1 else {
    cat("No adaptation method specified. Calculate kernel estimate with bandwidth hmax.\n")
    hinit <- hmax
  }
  hincr <- sqrt(1.25)
  if (demo && !graph) graph <- TRUE
  if(graph){
    oldpar <- par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
    on.exit(par(oldpar))
    graphobj0 <- object[-(1:length(object))[names(object)=="img"]]
    graphobj0$dim <- c(n1,n2)
  }
  # now check which procedure is appropriate
  #
  #    Initialize  list for theta
  #
  hmax <- hmax
  bi <- rep(1,n)
  theta <- switch(imgtype,greyscale=object$img[xind,yind],rgb=object$img[xind,yind,])
  bi0 <- 1
  #
  #  if varmodel specified prepare for initial variance estimation
  #
  if(estvar){
    coef <- matrix(0,nvarpar,dv)
    coef[1,] <- sigma2
    vobj <- list(coef=coef,meanvar=sigma2)
    imgq995 <- switch(imgtype,greyscale=quantile(object$img[xind,yind],.995),
                      rgb=apply(object$img[xind,yind,],3,quantile,.995))
  }
  if(scorr){
    twohp1 <- 2*trunc(hpre)+1
    pretheta <- .Fortran("awsimg",
                         as.integer(switch(imgtype,
                                           greyscale=object$img[xind,yind],
                                           rgb=object$img[xind,yind,])),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         as.double(hpre),
                         as.double(1.e40),
                         as.integer(theta),
                         bi=as.double(bi),
                         as.double(bi0),
                         theta=integer(prod(dimg)),
                         as.integer(lkern),
                         as.integer(1),
                         as.double(spmin),		       
                         as.double(spmax),
                         double(twohp1*twohp1),# array for location weights
                         as.double(wghts),
                         double(dv),DUP=FALSE,
                         PACKAGE="adimpro")$theta
    dim(pretheta) <- dimg
    spchcorr <- .Fortran("estcorr",
                         as.double(switch(imgtype,
                                          greyscale=object$img[xind,yind],
                                          rgb=object$img[xind,yind,])-pretheta),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         scorr=double(2*dv),
                         chcorr=double(max(1,dv*(dv-1)/2)),
                         as.double(hpre),
                         DUP=FALSE,
                         PACKAGE="adimpro")[c("scorr","chcorr")]
    spcorr <- spchcorr$scorr
    srh <- sqrt(hpre) 
    spcorr <- matrix(pmin(.9,spcorr+
                          bcf[1]/srh+bcf[2]/hpre+
                          bcf[3]*spcorr/srh+bcf[4]*spcorr/hpre+
                          bcf[5]*spcorr^2/srh+bcf[6]*spcorr^2/hpre),2,dv)
    #  bias correction for spatial correlation
    chcorr <- spchcorr$chcorr
    for(i in 1:dv) cat("Estimated spatial correlation in channel ",i,":",signif(spcorr[,i],2),"\n")
    if(dv>1) cat("Estimated correlation between channels (1,2):",signif(chcorr[1],2),
        " (1,3):",signif(chcorr[2],2),
        " (2,3):",signif(chcorr[3],2),"\n")
  } else {
    spcorr <- matrix(0,2,dv)
    chcorr <- numeric(max(1,dv*(dv-1)/2))
  }
  #
  #         fix values of the image in inactiv pixel
  #
  if(!is.null(mask)) fix <- !mask[xind,yind]
  ###
  ###              gridded   2D
  ###
  steps <- as.integer(log(hmax/hinit)/log(hincr)+1)
  if (length(lseq)<steps) lseq <- c(lseq,rep(1,steps-length(lseq)))
  lseq <- lseq[1:steps]
  hakt0 <- hakt <- hinit
  if(estvar) { 
    if(aws) hakt <- hakt*hincr
    step <- 1
  } else step <- 0 
  lambda0 <- 1e50
  progress <- 0
  total <- (hincr^(2*ceiling(log(hmax/hinit)/log(hincr)))-1)/(hincr^2-1)
  if (total == 0) total <- hincr^(2*step) # for (hmax == hinit)
  #
  #   run single steps to display intermediate results
  #
  while (hakt<=hmax) {
    twohp1 <- 2*trunc(hakt)+1
    spcorr <- pmax(apply(spcorr,1,mean),0)
    if(any(spcorr)>0) {
      h0<-numeric(length(spcorr))
      for(i in 1:length(h0))
        h0[i]<-geth.gauss(spcorr[i])
      if(length(h0)<2) h0<-rep(h0[1],2)
      # cat("Corresponding bandwiths for specified correlation:",h0,"\n")
    }
    if (any(spcorr>0)) {
      lcorr <- Spatialvar.gauss(hakt0,h0)/
        Spatialvar.gauss(h0,1e-5)/Spatialvar.gauss(hakt0,1e-5)
      # Correction C(h0,hakt) for spatial correlation depends on h^{(k-1)} 
      lambda0 <-lambda0*lcorr
      if(varmodel=="None") 
        lambda0 <- lambda0*Varcor.gauss(h0)
    } 
    hakt0 <- hakt
    if(is.null(mask)){
      if(estvar){
        varcase <- 1
        zobj <- .Fortran("awsvimg",
                         as.integer(switch(imgtype,
                                           greyscale=object$img[xind,yind],
                                           rgb=object$img[xind,yind,])),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         as.double(vobj$coef),
                         as.integer(nvarpar),
                         as.double(vobj$meanvar),
                         as.double(chcorr),
                         hakt=as.double(hakt),
                         as.double(lambda0),
                         as.integer(theta),
                         bi=as.double(bi),
                         bi0=as.double(bi0),# just take a scalar here
                         theta=integer(prod(dimg)),
                         as.integer(lkern),
                         as.integer(1),
                         as.double(spmin),		       
                         as.double(spmax),
                         as.double(sqrt(wghts)),
                         double(twohp1*twohp1),# array for location weights
                         double(dv),DUP=FALSE,
                         PACKAGE="adimpro")[c("bi","bi0","theta","hakt")]
      } else {
        # all other cases
        varcase <- 3
        zobj <- .Fortran("awsimg",
                         as.integer(switch(imgtype,
                                           greyscale=object$img[xind,yind],
                                           rgb=object$img[xind,yind,])),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         hakt=as.double(hakt),
                         as.double(lambda0),
                         as.integer(theta),
                         bi=as.double(bi),
                         bi0=as.double(bi0),
                         theta=integer(prod(dimg)),
                         as.integer(lkern),
                         as.integer(1),
                         as.double(spmin),
                         as.double(spmax),
                         double(twohp1*twohp1),# array for location weights
                         as.double(wghts),
                         double(dv),DUP=FALSE,
                         PACKAGE="adimpro")[c("bi","bi0","theta","hakt")]
      }
    } else {
      # all other cases
      varcase <- 3
      zobj <- .Fortran("mawsimg",
                       as.integer(switch(imgtype,
                                         greyscale=object$img[xind,yind],
                                         rgb=object$img[xind,yind,])),
                       as.logical(fix),
		       as.logical(mask[xind,yind]),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(dv),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.integer(theta),
                       bi=as.double(bi),
                       bi0=as.double(bi0),
                       theta=integer(prod(dimg)),
                       as.integer(lkern),
		       as.integer(1),
                       as.double(spmin),
                       as.double(spmax),
                       double(twohp1*twohp1),# array for location weights
                       as.double(wghts),
                       double(dv),DUP=FALSE,
                       PACKAGE="adimpro")[c("bi","bi0","theta","hakt")]
    }
    theta <- zobj$theta
    bi <- zobj$bi
    bi0 <- zobj$bi0
    rm(zobj)
    gc()
    dim(bi) <- dimg[1:2]
    if (graph) {
      graphobj <- graphobj0
      class(graphobj) <- "adimpro"
      graphobj$img <- switch(imgtype,
                             greyscale=object$img[xind,yind],
                             rgb=object$img[xind,yind,])
      show.image(graphobj,max.x=max.pixel,xaxt="n",yaxt="n")
      title("Observed Image")
      graphobj$img <- array(as.integer(theta),dimg)
      show.image(graphobj,max.x=max.pixel,xaxt="n",yaxt="n")
      title(paste("Reconstruction  h=",signif(hakt,3)))
      graphobj$img <- matrix(as.integer(65534*bi/bi0),n1,n2)
      graphobj$type <- "greyscale"
      graphobj$gamma <- FALSE
      show.image(graphobj,max.x=max.pixel,xaxt="n",yaxt="n")
      title(paste("Adaptation (rel. weights):",signif(mean(bi)/bi0,3)))
      rm(graphobj)
      gc()
    }
    progress <- progress + hincr^(2*step)
    cat("Bandwidth",signif(hakt,3)," Progress",signif(progress/total,2)*100,"% \n")
    if (scorr) {
      #
      #   Estimate Correlations  (keep old estimates until hmax > hpre)
      #
      if(hakt > hpre){
      spchcorr <- .Fortran("estcorr",
                           as.double(switch(imgtype,
                                            greyscale=object$img[xind,yind],
                                            rgb=object$img[xind,yind,]) - theta),
                           as.integer(n1),
                           as.integer(n2),
                           as.integer(dv),
                           scorr=double(2*dv),
                           chcorr=double(max(1,dv*(dv-1)/2)),
                           DUP=FALSE,
                           PACKAGE="adimpro")[c("scorr","chcorr")]
    spcorr <- spchcorr$scorr
#    spcorr <- matrix(pmin(.9,0.8817*spcorr+0.231/hakt+6.018*spcorr/hakt^2+
#                            1.753*spcorr^2/hakt-10.622*spcorr^2/hakt^2),2,dv)
    srh <- sqrt(hakt) 
    spcorr <- matrix(pmin(.9,spcorr+
                         bcf[1]/srh+bcf[2]/hakt+
                         bcf[3]*spcorr/srh+bcf[4]*spcorr/hakt+
                         bcf[5]*spcorr^2/srh+bcf[6]*spcorr^2/hakt),2,dv)
      #  bias correction for spatial correlation
      chcorr <- spchcorr$chcorr
      for(i in 1:dv) cat("Estimated spatial correlation in channel ",i,":",signif(spcorr[,i],2),"\n")
      if(dv>1) cat("Estimated correlation between channels (1,2):",signif(chcorr[1],2),
                   " (1,3):",signif(chcorr[2],2),
                   " (2,3):",signif(chcorr[3],2),"\n")
      } else {
           spcorr <-   matrix(pmin(.9,spchcorr$scorr+
                         bcf[1]/srh+bcf[2]/hpre+
                         bcf[3]*spchcorr$scorr/srh+bcf[4]*spchcorr$scorr/hpre+
                         bcf[5]*spchcorr$scorr^2/srh+bcf[6]*spchcorr$scorr^2/hpre),2,dv)         
# spcorr <- matrix(pmin(.9,.06+0.9*spchcorr$scorr+1.108/hpre),2,dv)
      #  bias correction for spatial correlation
         chcorr <- spchcorr$chcorr
      }
    } else {
      spcorr <- matrix(0,2,dv)
      chcorr <- numeric(max(1,dv*(dv-1)/2))
    }
    if (estvar) {
      #
      #   Create new variance estimate
      #
      vobj <- .Fortran(switch(varmodel,Constant="esigmac",Linear="esigmal"),
                       as.integer(switch(imgtype,
                                         greyscale=object$img[xind,yind],
                                         rgb=object$img[xind,yind,])),
                       as.integer(n1*n2),
                       as.integer(dv),
                       as.integer(theta),
                       as.double(bi),
                       as.integer(imgq995),
                       coef=double(nvarpar*dv),
                       meanvar=double(dv),
                       DUP=FALSE,
                       PACKAGE="adimpro")[c("coef","meanvar")]
      dim(vobj$coef) <- c(nvarpar,dv)
      cat("Estimated mean variance",signif(vobj$meanvar/65635^2,3),"\n")
    }
    step <- step + 1
    if (demo) readline("Press return")
    hakt <- hakt*hincr
    lambda0 <- lambda*lseq[step]
  }
  ###                                                                       
  ###            end cases                                                  
  ###                                 .....................................
  if(graph) par(oldpar)
  if(imgtype=="rgb") {
    object$img[xind,yind,] <- theta
  } else if(imgtype=="greyscale") {
    object$img[xind,yind] <- theta
  }
  #  if(dv==1) dim(img) <- dim(img)[1:2]
  ni <- array(1,dimg0[1:2])
  ni[xind,yind] <- bi
  object$ni <- ni
  object$ni0 <- bi0
  object$hmax <- hakt/hincr
  object$call <- args
  if(estvar) {
    if(varmodel=="Linear") {
       vobj$coef[1,] <- vobj$coef[1,]/65535^2
       vobj$coef[2,] <- vobj$coef[2,]/65535
    } else vobj$coef <- vobj$coef/65535^2
    object$varcoef <- vobj$coef
    object$wghts <- vobj$wghts
  }
  if(scorr) {
    object$scorr <- spcorr
    object$chcorr <- chcorr
  }
  invisible(if(compress) compress.image(object) else object)
}


#########################################################################################
#
#    local polynomial version
#
#########################################################################################  
awspimage <- function(object, hmax=12, aws=TRUE, degree=1,
                      varmodel=NULL, ladjust=1.0, xind = NULL, yind = NULL, 
                      wghts=c(1,1,1,1), scorr=TRUE, lkern="Triangle",demo=FALSE,
                      graph=FALSE, max.pixel=4.e2, clip=FALSE, compress=TRUE){
  #
  #          Auxilary functions
  #
  IQRdiff <- function(img) IQR(diff(img))/1.908
  #
  #     Compute theta
  #
  gettheta <- function(ai,bi){
    n1 <- dim(ai)[1]
    n2 <- dim(ai)[2]
    n <- n1*n2
    dp1 <- dim(ai)[3]
    dp2 <- dim(bi)[3]
    dv <- dim(ai)[4]
    ind <- matrix(c(1, 2, 3, 4, 5, 6,
                    2, 4, 5, 7, 8, 9,
                    3, 5, 6, 8, 9,10,
                    4, 7, 8,11,12,13,
                    5, 8, 9,12,13,14,
                    6, 9,10,13,14,15),6,6)[1:dp1,1:dp1]
    theta <- ai
    for(i in 1:dv) {
      theta[,,,i] <- array(.Fortran("mpaws2",
                                    as.integer(n),
                                    as.integer(dp1),
                                    as.integer(dp2),
                                    as.double(ai[,,,i]),
                                    as.double(bi),
                                    theta=double(dp1*n),
                                    double(dp1*dp1),
                                    as.integer(ind),PACKAGE="adimpro")$theta,c(n1,n2,dp1))
    }
    theta
  }
  #####################################################################################
  ###    
  ###    function body
  ###    
  ###    first check arguments and initialize
  ###
  #####################################################################################    
  if(!check.adimpro(object)) {
    cat(" Consistency check for argument object failed (see warnings). object is returned.\"n")
    return(invisible(object)) 
  }
  object <- switch(object$type,
                   "hsi" = hsi2rgb(object),
                   "yuv" = yuv2rgb(object),
                   "yiq" = yiq2rgb(object),
                   "xyz" = xyz2rgb(object),
                   object)
  if(!is.null(object$call)) { 
    warning("argument object is result of awsimage or awspimage. You should know what you do.\"n")
  }
  args <- match.call()
  if(!(degree %in% c(1,2))) {
    warning("only polynomial degrees 1 and 2 implemented, original image is returned")
    return(img)
  }
  if(clip) object <- clip.image(object,xind,yind,compress=FALSE)
  if(object$compressed) object <- decompress.image(object)
  if(is.null(varmodel)) {
    varmodel <- if(object$gamma) "Constant" else "Linear"
  }
  estvar <- toupper(varmodel) %in% c("CONSTANT","LINEAR")
  hpre <- switch(degree,2.,3.)
  bcf <- switch(degree,c(-.4724,1.1969,.22167,.48691,.13731,-1.0648),
                             c(-.69249,2.0552,.05379,1.6063,.28891,-1.907))
# coefficients for  bias correction for spatial correlation
  if(estvar) {
    dlw<-(2*trunc(hpre)+1)
    nvarpar <- switch(varmodel,Constant=1,Linear=2,1)
  } else nvarpar <- 1
  #
  #   Check image type
  #
  dimg <- dimg0 <- dim(object$img)
  imgtype <- object$type 
  dv <- switch(imgtype,rgb=dimg[3],greyscale=1)
  if(length(wghts)<dv) wghts <- c(wghts,rep(1,dv-length(wghts))) else wghts <- pmin(.1,wghts) 
  wghts <- switch(imgtype, greyscale=1, rgb = wghts[1:dv])
  spcorr <- matrix(0,2,dv)
  h0 <- c(0,0)
  #
  #   specify region of interest
  #
  if(is.null(xind)) xind <- 1:dimg[1]
  if(is.null(yind)) yind <- 1:dimg[2]
  n1 <- length(xind)
  n2 <- length(yind)
  n <- n1*n2
  dimg <- c(n1,n2,dv)
  #
  #     set approriate defaults
  #
  dp1 <- switch(degree+1,1,3,6)
  dp2 <- switch(degree+1,1,6,15)
  if(aws) hinit <- 1 else {
    cat("No adaptation method specified. Calculate local polynomial estimate with bandwidth hmax.\n")
    hinit <- hmax
  }
  if (aws) qlambda <- switch(degree,.9,.995) else qlambda <- 1
  lkern <- switch(lkern,
                  Triangle=2,
                  Quadratic=3,
                  Cubic=4,
                  Uniform=1,
                  2)
  if (is.null(hmax)) hmax <- 12
  wghts <- wghts/max(wghts)
  #
  #    Initialize  list for theta
  #
  if (qlambda<1) lambda <- ladjust*2.e0*qchisq(qlambda,sum(wghts)^2/sum(wghts^2)*dp2) else lambda <- 1e50
  sigma2 <- switch(imgtype,
                   "greyscale" = IQRdiff(object$img)^2,
                   "rgb" = apply(object$img,3,IQRdiff)^2, IQRdiff(object$img)^2)
  cat("Estimated variance: ", signif(sigma2/65635^2,4),"\n")
  if (!estvar) {
    vobj <- list(coef= sigma2, meanvar=sigma2)
    spcorr <- matrix(0,2,dv)
  }
    #
    # set the support of the statistical kernel to (0,1), set spmin and spmax
    #
    lambda <- 2*lambda
    spmin <- 0
    spmax <- 1
  bi <- array(rep(1,n*dp2),c(dimg[1:2],dp2))
  theta <- array(0,c(dimg,dp1))
  theta[,,,1] <- switch(imgtype,greyscale=object$img[xind,yind],rgb=object$img[xind,yind,])
  bi0 <- 1
  ind <- matrix(c(1, 2, 3, 4, 5, 6,
                  2, 4, 5, 7, 8, 9,
                  3, 5, 6, 8, 9,10,
                  4, 7, 8,11,12,13,
                  5, 8, 9,12,13,14,
                  6, 9,10,13,14,15),6,6)[1:dp1,1:dp1]
  #
  #  if varmodel specified prepare for initial variance estimation
  #
  if(estvar){
    coef <- matrix(0,nvarpar,dv)
    coef[1,] <- sigma2
    vobj <- list(coef=coef,meanvar=sigma2)
    imgq995 <- switch(imgtype,greyscale=quantile(object$img[xind,yind],.995),
                      rgb=apply(object$img[xind,yind,],3,quantile,.995))
  }
  if(scorr){
    twohp1 <- 2*trunc(hpre)+1
    pretheta <- .Fortran("awsimg",
                         as.integer(switch(imgtype,
                                           greyscale=object$img[xind,yind],
                                           rgb=object$img[xind,yind,])),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         as.double(hpre),
                         as.double(1.e40),
                         as.integer(switch(imgtype,
                                           greyscale=object$img[xind,yind],
                                           rgb=object$img[xind,yind,])),
                         bi=as.double(bi),
                         as.double(bi0),
                         theta=integer(prod(dimg)),
                         as.integer(lkern),
                         as.integer(1),
                         as.double(spmin),		       
                         as.double(spmax),
                         double(twohp1*twohp1),# array for location weights
                         double(twohp1*twohp1),# array for general weights
                         as.double(wghts),
                         double(dv),DUP=FALSE,
                         PACKAGE="adimpro")[c("theta","bi")]
     prebi <- pretheta$bi
     pretheta <- pretheta$theta
    #
    # just initialize, value of sigma2 has no influence on estimates if hakt==hinit==1
    #
    #     gamma <- 0
    spchcorr <- .Fortran("estcorr",
                         as.double(switch(imgtype,
                                          greyscale=object$img[xind,yind],
                                          rgb=object$img[xind,yind,])-pretheta),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         scorr=double(2*dv),
                         chcorr=double(max(1,dv*(dv-1)/2)),
                         DUP=FALSE,
                         PACKAGE="adimpro")[c("scorr","chcorr")]
    spcorr <- spchcorr$scorr
    srh <- sqrt(hpre) 
    spcorr <- matrix(pmin(.9,spcorr+
                         bcf[1]/srh+bcf[2]/hpre+
                         bcf[3]*spcorr/srh+bcf[4]*spcorr/hpre+
                         bcf[5]*spcorr^2/srh+bcf[6]*spcorr^2/hpre),2,dv)
#  bias correction for spatial correlation
    chcorr <- spchcorr$chcorr
    for(i in 1:dv) cat("Estimated spatial correlation in channel ",i,":",signif(spcorr[,i],2),"\n")
    if(dv>1) cat("Estimated correlation between channels (1,2):",signif(chcorr[1],2),
        " (1,3):",signif(chcorr[2],2),
        " (2,3):",signif(chcorr[3],2),"\n")
  } else {
    spcorr <- matrix(0,2,dv)
    chcorr <- numeric(max(1,dv*(dv-1)/2))
  }
    if (estvar) {
      #
      #   Create new variance estimate
      #
      vobj <- .Fortran(switch(varmodel,Constant="epsigmac",Linear="epsigmal"),
                       as.integer(switch(imgtype,
                                         greyscale=object$img[xind,yind],
                                         rgb=object$img[xind,yind,])),
                       as.integer(n1*n2),
                       as.integer(dv),
                       as.integer(pretheta),
                       as.double(prebi),
                       as.integer(imgq995),
                       coef=double(nvarpar*dv),
                       meanvar=double(dv),
                       as.integer(dp1),
                       DUP=FALSE,
                       PACKAGE="adimpro")[c("coef","meanvar")]
      if(any(is.na(vobj$coef))) vobj <- oldvobj
      dim(vobj$coef) <- c(nvarpar,dv)
      cat("Estimated mean variance",signif(vobj$meanvar/65635^2,3),"\n")
    }
  rm(prebi)
  #
  #     now set hinit and hincr if not provided
  #
  hincr <- sqrt(1.25)
  if (demo && !graph) graph <- TRUE
  if(graph){
    oldpar <- par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
    on.exit(par(oldpar))
    graphobj0 <- object[-(1:length(object))[names(object)=="img"]]
    graphobj0$dim <- c(n1,n2)
  }
  hw <- degree+.1
  steps <- as.integer(log(hmax/hinit)/log(hincr)+1)
  hakt0 <- hakt <- hinit
  lambda0 <- 1e50
  progress <- 0
  step <- 0
  total <- (hincr^(2*ceiling(log(hmax/hinit)/log(hincr)))-1)/(hincr^2-1)
  #
  #   run single steps to display intermediate results
  #
  while (hakt<=hmax) {
    twohp1 <- 2*trunc(hakt)+1
    twohhwp1 <- 2*trunc(hakt+hw)+1
    spcorr <- pmax(apply(spcorr,1,mean),0)
    if(any(spcorr)>0) {
      h0<-numeric(length(spcorr))
      for(i in 1:length(h0))
        h0[i]<-geth.gauss(spcorr[i])
      if(length(h0)<2) h0<-rep(h0[1],2)
      # cat("Corresponding bandwiths for specified correlation:",h0,"\n")
    }
    if (any(spcorr>0)) {
      lcorr <- Spatialvar.gauss(hakt0,h0)/
        Spatialvar.gauss(h0,1e-5)/Spatialvar.gauss(hakt0,1e-5)
      # Correction C(h0,hakt) for spatial correlation depends on h^{(k-1)} 
      lambda0 <-lambda0*lcorr
      if(varmodel=="None") 
        lambda0 <- lambda0*Varcor.gauss(h0)
    } 
    zobj <- .Fortran("awspimg",
                     as.integer(switch(imgtype,
                                       greyscale=object$img[xind,yind],
                                       rgb=object$img[xind,yind,])),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(dv),
                     as.integer(degree),
                     as.double(hw),
                     as.double(vobj$coef),
                     as.integer(nvarpar),
                     as.double(vobj$meanvar),
                     hakt=as.double(hakt),
                     as.double(lambda0),
                     as.double(theta),
                     bi=as.double(bi),
                     bi0=double(1),
                     ai=double(n*dp1*dv),
                     as.integer(lkern),
                     as.integer(1),
                     as.double(spmin),
                     as.double(spmax),
                     double(twohp1*twohp1),# array for location weights
                     double(twohp1*twohp1),# array for general weights
                     double(twohhwp1*twohhwp1),# array for smoothed location weights
                     double(twohhwp1*twohhwp1),# array for smoothed general weights
                     as.double(wghts),
                     as.integer(ind),
                     PACKAGE="adimpro")[c("bi","bi0","ai","hakt")]
    dim(zobj$ai) <- c(dimg[1:2],dp1,dv)
    bi <- zobj$bi
    dim(bi)<-c(n1,n2,dp2)
    theta <- gettheta(zobj$ai,bi)
    dim(theta)<-c(dimg[1:2],dp1,dv)
    bi0 <- max(zobj$bi0,bi[,,1])
    rm(zobj)
    gc()
    if (graph) {
      graphobj <- graphobj0
      class(graphobj) <- "adimpro"
      graphobj$img <- switch(imgtype,
                             greyscale=object$img[xind,yind],
                             rgb=object$img[xind,yind,])
      show.image(graphobj,max.x=max.pixel,xaxt="n",yaxt="n")
      title("Observed Image")
      graphobj$img <- switch(imgtype,
                             greyscale=array(pmin(65535,pmax(0,as.integer(theta[,,1,]))),dimg[1:2]),
                             rgb=array(pmin(65535,pmax(0,as.integer(theta[,,1,]))),dimg))
      show.image(graphobj,max.x=max.pixel,xaxt="n",yaxt="n")
      title(paste("Reconstruction  h=",signif(hakt,3)))
      graphobj$img <- matrix(pmin(65535,pmax(0,as.integer(65534*bi[,,1]/bi0))),n1,n2)
      graphobj$type <- "greyscale"
      graphobj$gamma <- FALSE
      show.image(graphobj,max.x=max.pixel,xaxt="n",yaxt="n")
      title(paste("Adaptation (rel. weights):",signif(mean(bi[,,1])/bi0,3)))
      rm(graphobj)
      gc()
    }
    progress <- progress + hincr^(2*step)
    cat("Bandwidth",signif(hakt,3)," Progress",signif(progress/total,2)*100,"% \n")
    if (scorr) {
      #
      #   Estimate Correlations
      #
      if(hakt > hpre){
      spchcorr <- .Fortran("estcorr",
                           as.double(switch(imgtype,
                                            greyscale=object$img[xind,yind],
                                            rgb=object$img[xind,yind,])- theta[,,1,]),
                           as.integer(n1),
                           as.integer(n2),
                           as.integer(dv),
                           scorr=double(2*dv),
                           chcorr=double(max(1,dv*(dv-1)/2)),
                           DUP=FALSE,
                           PACKAGE="adimpro")[c("scorr","chcorr")]
    spcorr <- spchcorr$scorr
    srh <- sqrt(hakt) 
    spcorr <- matrix(pmin(.9,spcorr+
                         bcf[1]/srh+bcf[2]/hakt+
                         bcf[3]*spcorr/srh+bcf[4]*spcorr/hakt+
                         bcf[5]*spcorr^2/srh+bcf[6]*spcorr^2/hakt),2,dv)
      #  bias correction for spatial correlation
      chcorr <- spchcorr$chcorr
      for(i in 1:dv) cat("Estimated spatial correlation in channel ",i,":",signif(spcorr[,i],2),"\n")
      if(dv>1) cat("Estimated correlation between channels (1,2):",signif(chcorr[1],2),
                   " (1,3):",signif(chcorr[2],2),
                   " (2,3):",signif(chcorr[3],2),"\n")
      } else {
        spcorr <- matrix(pmin(.9,spcorr+
                         bcf[1]/srh+bcf[2]/hakt+
                         bcf[3]*spcorr/srh+bcf[4]*spcorr/hakt+
                         bcf[5]*spcorr^2/srh+bcf[6]*spcorr^2/hakt),2,dv)
        chcorr <- spchcorr$chcorr
      }
    } else {
      spcorr <- matrix(0,2,dv)
      chcorr <- numeric(max(1,dv*(dv-1)/2))
    }
    if (estvar && hakt > hpre) {
      #
      #   Create new variance estimate
      #
      oldvobj <- vobj
      vobj <- .Fortran(switch(varmodel,Constant="epsigmac",Linear="epsigmal"),
                       as.integer(switch(imgtype,
                                         greyscale=object$img[xind,yind],
                                         rgb=object$img[xind,yind,])),
                       as.integer(n1*n2),
                       as.integer(dv),
                       as.integer(theta[,,1,]),
                       as.double(bi),
                       as.integer(imgq995),
                       coef=double(nvarpar*dv),
                       meanvar=double(dv),
                       as.integer(dp1),
                       DUP=FALSE,
                       PACKAGE="adimpro")[c("coef","meanvar")]
      if(any(is.na(vobj$coef))) vobj <- oldvobj
      dim(vobj$coef) <- c(nvarpar,dv)
      cat("Estimated mean variance",signif(vobj$meanvar/65635^2,3),"\n")
    }
    if (demo) readline("Press return")
    step <- step + 1
    hakt <- hakt*hincr
    lambda0 <- lambda
  }
  ###                                                                       
  ###            end cases                                                  
  ###                                 
  ###   component var contains an estimate of Var(theta) 
  ###   
  if(graph) par(oldpar)
  if(imgtype=="rgb") {
    object$img[xind,yind,] <- as.integer(pmin(65535,pmax(0,theta[,,1,])))
  } else if(imgtype=="greyscale") {
    object$img[xind,yind] <- as.integer(pmin(65535,pmax(0,theta[,,1,1])))
  }
  ni <- array(1,dimg0[1:2])
  ni[xind,yind] <- bi[,,1]
  object$ni <- ni
  object$ni0 <- bi0
  object$hmax <- hakt/hincr
  object$call <- args
  if(estvar) {
    if(varmodel=="Linear") {
       vobj$coef[1,] <- vobj$coef[1,]/65535^2
       vobj$coef[2,] <- vobj$coef[2,]/65535
    } else vobj$coef <- vobj$coef/65535^2
    object$varcoef <- vobj$coef
    object$wghts <- vobj$wghts
  }
  if(scorr) {
    object$scorr <- spcorr
    object$chcorr <- chcorr
  }
  invisible(if(compress) compress.image(object) else object)
}








