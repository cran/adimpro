#########################################################################################
#
#    local polynomial version
#
#########################################################################################  


awspimage.old <- function(img, hmax=12, aws=TRUE, sagg=FALSE, degree=0, 
                      sigma2= NULL,  wghts=c(1,1,1,1), scorr=0,
                      lkern="Triangle", aggkern="Triangle", colorspace="YUV", demo=FALSE,
                      graph=FALSE, max.pixel=1.e3){
  #
  #          Auxilary functions
  #
  IQRdiff <- function(img) IQR(diff(img))/1.908
  #
  #     Quadratic distance for memory penalty
  #
  Pardist <- function(Bi0,dtheta,dv,tau){
  #  local polynomial bi  
    dp1 <- dim(dtheta)[3]
    dp2 <- dim(Bi0)[3]
    ind <- matrix(c(1, 2, 3, 4, 5, 6,
                    2, 4, 5, 7, 8, 9,
                    3, 5, 6, 8, 9,10,
                    4, 7, 8,11,12,13,
                    5, 8, 9,12,13,14,
                    6, 9,10,13,14,15),6,6)[1:dp1,1:dp1,drop=FALSE]
    dist <- 0
   # control in first channel only
    for(k in 1:dv)
    for(i in 1:dp1) for(j in 1:dp1) dist <-
      dist+dtheta[,,i,k]*Bi0[,,ind[i,j]]*dtheta[,,j,k]/tau[k]
    dist
  }
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
  #
  #  update theta
  #
  updtheta <- function(zobj,tobj,cpar,aggkern) {
    heta <- cpar$heta
    tau1 <- cpar$tau1
    tau2 <- cpar$tau2
    kstar <- cpar$kstar
    hakt <- zobj$hakt
    tau <- 2*(tau1+tau2*max(kstar-log(hakt),0))
    bi <- zobj$bi
    bi0 <- zobj$bi0
    dimg <- cpar$dimg
    dv <- cpar$dv
    thetanew <- gettheta(zobj$ai,bi)
    theta <- tobj$theta
    dim(theta)<-dim(thetanew)
    n<-prod(dim(tobj$fix))
    thetanew[array(tobj$fix,dim(thetanew))] <- theta[rep(tobj$fix,dp1*dv)]
    if (hakt>heta) {
      KLdistth<-Pardist(bi0,thetanew-theta,dv,tau)
      eta <- switch(aggkern,
                    "Uniform"=as.numeric(KLdistth>1),
                    "Triangle"=pmin(1,KLdistth),
         	    as.numeric(KLdistth>1))
      eta <- rep(eta,dp2)
    } else {
      eta <- rep(0,n*dp2)
    }
    if(any(tobj$fix)) eta[rep(tobj$fix,dp2)] <- 1
    bi <- (1-eta)*bi + eta * tobj$bi
    eta <- array(eta[1:n],dim(theta))
    theta <- (1-eta)*thetanew + eta * theta
    list(theta=theta,bi=bi,eta=eta[,,1,1],fix=(eta[,,1,1]==1))
  }
  #
  #    first check arguments and initialize
  #
  args <- match.call()
  #
  #   Check image type
  #
  dimg <- dim(img)
  if(is.null(attr(img,"type"))) {
    attr(img,"type") <- switch(as.character(length(dimg)),
                               "2" = "greyscale",
                               "3" = "RGB",
                               "unknown")
  }
  if (toupper(attr(img,"type"))=="GREYSCALE" && length(dimg) != 2) return(warning("incorrect grayscale image dimensions"))
  if (toupper(attr(img,"type"))=="RGB" && (length(dimg) != 3 || !(dimg[3] %in% 3:4))) return(warning("incorrect RGB image dimensions"))
  if (!(toupper(attr(img,"type")) %in% c("GREYSCALE","RGB"))) return(warning("incorrect image type"))
  if (toupper(attr(img,"type"))=="GREYSCALE") wghts <- 1 else wghts <- wghts[1:dimg[3]]
  dv <- switch(toupper(attr(img,"type")),
               "GREYSCALE"=1,
               "RGB"=dimg[3])
  nwghts <- switch(toupper(attr(img,"type")),
                   "GREYSCALE"=1,
                   "RGB"=max(dv,sum(wghts>0)))
  #
  #   convert image if necessary (assumes original image in RGB
  #
  if (toupper(attr(img,"type"))=="RGB") {
    img <- switch(colorspace,
                  "YUV"=rgb2yuv(img),
                  "YIQ"=rgb2yiq(img),
                  img)
    attr(img,"type")<-toupper(colorspace)
  }
  #
  #   set cut off point in K_{st}(x) = exp(-x) I_{x<spmax}
  #
  #
  #     set approriate defaults
  #
  dp1 <- switch(degree+1,1,3,6)
  dp2 <- switch(degree+1,1,6,15)
  if(sagg||aws) hinit <- 1 else {
      cat("No adaptation method specified. Calculate local polynomial estimate with bandwidth hmax.\n")
      hinit <- hmax
  }
  heta <- degree+3
  if (aws) qlambda <- switch(degree+1,.975,.85,.99) else qlambda <- 1
  lseq <- switch(degree+1,c(1.85,1.3,1.1,1.1),
                          c(1.85,1.3,1.1,1.1),c(1.85,1.3,1.1,1.1))
  lkern <- switch(lkern,
                  Triangle=2,
                  Quadratic=3,
                  Cubic=4,
                  Uniform=1,
                  2)
  if (qlambda>=1) {
    # thats stagewise aggregation with kernel specified by aggkern
    if (sagg) qtau <- switch(degree+1,.25,.15,.2) else qtau <- 1
    if (qtau==1) tau1 <- 1e50 else tau1 <- qchisq(qtau,dv)
    if (aggkern=="Triangle") tau1 <- 2.5*tau1
    tau2 <- tau1/2
  } else {
    if (sagg) qtau<-switch(degree+1,.95,.6,.92) else qtau <- 1
    if (qtau>=1) {
      tau1 <- 1e50 
      heta <- 1e10
    } else {
      tau1 <- qchisq(qtau,dv)
    }
    if (aggkern=="Triangle") tau1 <- 2.5*tau1
    tau2 <- tau1/2
  }
  cat("tau1",tau1,"heta",heta,"\n")
  if (is.null(hmax)) hmax <- 12
  cpar <- list(heta=heta,tau1=tau1,tau2=tau2,dimg=dimg,dv=dv,kstar=log(15*dp1))
  n <- length(img)
  wghts0 <- wghts/max(wghts)
  #
  #  in case of colored noise get the corresponding bandwidth (for Gaussian kernel)
  #
  #  specify which statistics are needed and transform data if necessary
  #
  #  now set hinit and hincr if not provided
  #
  hincr <- sqrt(1.25)
  if (demo && !graph) graph <- TRUE
  # now check which procedure is appropriate
  ##  this is the version on a grid
  n1 <- dimg[1]
  n2 <- dimg[2]
  n <- n1*n2
  ddim <- 2
  #
  #    Initialize  list for theta
  #
  if (qlambda<1) lambda <- qchisq(qlambda,sum(wghts0)^2/sum(wghts0^2)*dp1) else lambda <- 1e50
  if (scorr>0) {
    h0 <- geth0(scorr,,1)
    cat("Corrsponding bandwidth for specified correlation:",h0,"\n")
  }
  if (is.null(sigma2)) {
    sigma2 <- switch(attr(img,"type"),
                     "greyscale" = IQRdiff(img)^2,
                     "YIQ" = apply(img[,,],3,IQRdiff)^2, 
                     "YUV" = apply(img[,,],3,IQRdiff)^2,
                     "RGB" = apply(img[,,],3,IQRdiff)^2,
                     IQRdiff(img)^2)
    if(scorr>0) sigma2<-sigma2*Varfactor("Gaussian",h0,d=2)
    cat("Estimated variance: ", signif(sigma2,4),"\n")
  }
  if (length(sigma2)==1) {
    #   homoskedastic Gaussian case
    lambda <- lambda*sigma2*2 
    cpar$tau1 <- cpar$tau1*sigma2*2 
    cpar$tau2 <- cpar$tau2*sigma2*2 
    cat("value of lambda",lambda/sigma2,"\n")
  } else if (length(sigma2)!=n) {
    wghts <- wghts0/sigma2
    cpar$tau1 <- cpar$tau1*sigma2*2 
    cpar$tau2 <- cpar$tau2*sigma2*2 
    lambda <- lambda*2 
    cat("value of lambda",lambda,"\n")
  } else {
    #   heteroskedastic Gaussian case
    if (length(sigma2)!=n) stop("sigma2 does not have length 1 or same length as img")
    lambda <- lambda*2 
    cat("value of lambda",lambda,"\n")
    cpar$tau1 <- cpar$tau1*2 
    cpar$tau2 <- cpar$tau2*2 
    sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes 
  }
  spmin <- 0
  spmax <- 3
  tobj <- list(bi= rep(1,n*dp2), theta= array(img,c(dim(img),dp1)), fix=rep(FALSE,n))
  bi0old <- rep(1,n*dp2)
  ind <- matrix(c(1, 2, 3, 4, 5, 6,
                  2, 4, 5, 7, 8, 9,
                  3, 5, 6, 8, 9,10,
                  4, 7, 8,11,12,13,
                  5, 8, 9,12,13,14,
                  6, 9,10,13,14,15),6,6)[1:dp1,1:dp1]
  hw <- degree+.1
  ###
  ###              gridded   ( 2D )
  ###
  steps <- as.integer(log(hmax/hinit)/log(hincr)+1)
  if (length(lseq)<steps) lseq <- c(lseq,rep(1,steps-length(lseq)))
  lseq <- lseq[1:steps]
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
    if (length(sigma2)==n) {
      # heteroskedastic Gaussian case
      if (scorr>0) lambda0<-lambda0*Spatialvariance(lkern,"Gaussian",hakt0,h0,2)/Spatialvariance(lkern,"Gaussian",hakt0,1e-5,2)
      hakt0 <- hakt
      zobj <- .Fortran("awsphimg",as.double(img),
                       as.double(sigma2),
                       as.logical(tobj$fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(dv),
                       as.integer(degree),
		       as.double(hw),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta[1:(n*dp1*nwghts)]),
                       bi=as.double(tobj$bi),
                       bi0=double(n*dp2),
                       ai=double(n*dp1*dv),
                       as.integer(lkern),
                       as.double(spmin),
                       as.double(spmax),
                       double(twohp1*twohp1),# array for location weights
                       double(twohp1*twohp1),# array for general weights
                       double(twohhwp1*twohhwp1),# array for smoothed location weights
                       double(twohhwp1*twohhwp1),# array for smoothed general weights
                       as.double(wghts),
                       as.integer(nwghts),
                       as.integer(ind),
                       PACKAGE="adimpro")[c("bi","bi0","ai","hakt")]
    } else {
      # all other cases
      zobj <- .Fortran("awspimg",as.double(img),
                       as.logical(tobj$fix),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(dv),
                       as.integer(degree),
		       as.double(hw),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(tobj$theta[1:(n*dp1*nwghts)]),
                       bi=as.double(tobj$bi),
                       bi0=double(n*dp2),
                       ai=double(n*dp1*dv),
                       as.integer(lkern),
                       as.double(spmin),
                       as.double(spmax),
                       double(twohp1*twohp1),# array for location weights
                       double(twohp1*twohp1),# array for general weights
                       double(twohhwp1*twohhwp1),# array for smoothed location weights
                       double(twohhwp1*twohhwp1),# array for smoothed general weights
                       as.double(wghts),
                       as.integer(nwghts),
                       as.integer(ind),
                       PACKAGE="adimpro")[c("bi","bi0","ai","hakt")]
    }
    dim(zobj$ai) <- c(dimg[1:2],dp1,dv)
    if (hakt>n1/2) zobj$bi0 <- hincr^ddim*biold
    biold <- zobj$bi0
    dim(zobj$bi0)<-c(n1,n2,dp2)
    tobj <- updtheta(zobj,tobj,cpar,aggkern)
    rm(zobj)
    gc()
    dim(tobj$theta) <- c(dimg[1:2],dp1,dv)
    dim(tobj$bi) <- c(dimg[1:2],dp2)
    dim(tobj$eta) <- dimg[1:2]
    if (graph) {
      par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
      switch(toupper(attr(img,"type")),
             "GREYSCALE" = show.image(img,max.x=max.pixel),
             "RGB" = show.image(img,max.x=max.pixel),
             "YUV" = show.image(yuv2rgb(img),max.x=max.pixel),
             "YIQ" = show.image(yiq2rgb(img),max.x=max.pixel))
      title("Observed Image")
      switch(toupper(attr(img,"type")),
             "GREYSCALE" = show.image(tobj$theta[,,1,],max.x=max.pixel),
             "RGB" = show.image(tobj$theta[,,1,],max.x=max.pixel),
             "YUV" = show.image(yuv2rgb(tobj$theta[,,1,]),max.x=max.pixel),
             "YIQ" = show.image(yiq2rgb(tobj$theta[,,1,]),max.x=max.pixel))
      title(paste("Reconstruction  h=",signif(hakt,3)))
      show.image(tobj$bi[,,1],max.x=max.pixel)
      title(paste("Sum of weights: min=",signif(min(tobj$bi[,,1]),3)," mean=",signif(mean(tobj$bi[,,1]),3)," max=",signif(max(tobj$bi[,,1]),3)))
      show.image(tobj$eta,max.x=max.pixel)
      title(paste("eta   max=",signif(max(tobj$eta),3)))
    }
    if (demo) readline("Press return")
    progress <- progress + hincr^(2*step)
    step <- step + 1
    cat(signif(progress/total,2)*100,"% . ",sep="")
    hakt <- hakt*hincr
    lambda0 <- lambda*lseq[step]
   }
  ###                                                                       
  ###            end cases                                                  
  ###                                 
  ###   component var contains an estimate of Var(tobj$theta) if and aggkern="Uniform", or if qtau1=1 
  ###   
  z <- list(hat=switch(toupper(attr(img,"type")),
                       "GREYSCALE" = tobj$theta[,,1,],
                       "RGB" = tobj$theta[,,1,],
                       "YUV" = yuv2rgb(tobj$theta[,,1,]),
                       "YIQ" = yiq2rgb(tobj$theta[,,1,])),
            ni=tobj$bi[,,1],
	    sigma2=sigma2,
            hmax=hakt/hincr,
            call=args)
  class(z) <- "awspimg.gaussian"
  z
}
