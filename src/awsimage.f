CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant  aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mawsimg(y,fix,mask,n1,n2,dv,hakt,lambda,theta,bi,
     1       bi0,thnew,kern,skern,spmin,spmax,lw,wght,swjy)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew       \sum  Wi Y     (output)
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldistd,lkern
      real*8 kldistd,lkern
      integer n1,n2,dv,kern,skern,y(1),theta(1),thnew(1)
      logical aws,fix(1),mask(1)
      real*8 bi(1),bi0,lambda,spmax,spmin,wght(dv),hakt,lw(1)
      integer ih,ih1,i1,i2,j1,j2,k,n,
     1        iind,jind,jind2,jwind2,dlw,clw,jw1,jw2
      real*8 bii,sij,swj,swj0,swjy(dv),z1,z2,wj,hakt2,bii0,spf
      hakt2=hakt*hakt
C      spf=spmax/(spmax-spmin)
      spf=spmax/(spmax-spmin)
      ih=hakt
      dlw=2*ih+1
      clw=ih+1
      aws=lambda.lt.1d40
      n=n1*n2
      bii0=bi0
      swj0=0.d0
C   compute location weights first
      DO j2=1,dlw
         z2=clw-j2
         z2=z2*z2
         ih1=dsqrt(hakt2-z2)
         jind2=(j2-1)*dlw
         DO j1=clw-ih1,clw+ih1
C  first stochastic term
            jind=j1+jind2
            z1=clw-j1
            wj=lkern(kern,(z1*z1+z2)/hakt2)
            swj0=swj0+wj
            lw(jind)=wj
         END DO
      END DO
      bi0=swj0
      call rchkusr()
      DO i2=1,n2
         DO i1=1,n1
            iind=i1+(i2-1)*n1
            IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
            bii=bi(iind)/lambda
C   scaling of sij outside the loop
            swj=0.d0
            DO k=1,dv
               swjy(k)=0.d0
            END DO
            DO jw2=1,dlw
	       j2=jw2-clw+i2
	       if(j2.lt.1.or.j2.gt.n2) CYCLE
	       jind2=(j2-1)*n1
               jwind2=(jw2-1)*dlw
               z2=clw-jw2
               ih1=dsqrt(hakt2-z2*z2)
               DO jw1=clw-ih1,clw+ih1
		  j1=jw1-clw+i1
	          if(j1.lt.1.or.j1.gt.n1) CYCLE
		  jind=j1+jind2
		  if(.not.mask(jind)) CYCLE
		  wj=lw(jw1+jwind2)
                  IF (aws) THEN
                     sij=bii*kldistd(theta(iind),theta(jind),n,wght,dv)
                     IF (sij.gt.spmax) CYCLE
		     IF (skern.eq.1) THEN
C  skern == "Triangle"
                        wj=wj*(1.d0-sij)
		     ELSE
C  skern == "Exp"
		        IF (sij.gt.spmin) wj=wj*dexp(-spf*(sij-spmin))
		     ENDIF
                  END IF
                  swj=swj+wj
                  DO k=1,dv
                     swjy(k)=swjy(k)+wj*y(jind+(k-1)*n)
                  END DO
               END DO
            END DO
            DO k=1,dv
               thnew(iind+(k-1)*n)=swjy(k)/swj
            END DO
            bi(iind)=swj
            call rchkusr()
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant  aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsimg(y,n1,n2,dv,hakt,lambda,theta,bi,bi0,
     1       thnew,kern,skern,spmin,spmax,lw,wght,swjy)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew       \sum  Wi Y     (output)
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldistd,lkern
      real*8 kldistd,lkern
      integer n1,n2,dv,kern,skern,y(1),theta(1),thnew(1)
      logical aws
      real*8 bi(1),bi0,lambda,spmax,spmin,wght(dv),hakt,lw(1)
      integer ih,ih1,i1,i2,j1,j2,k,n,
     1        iind,jind,jind2,jwind2,dlw,clw,jw1,jw2
      real*8 bii,sij,swj,swj0,swjy(dv),z1,z2,wj,hakt2,bii0,spf
      hakt2=hakt*hakt
C      spf=spmax/(spmax-spmin)
      spf=spmax/(spmax-spmin)
      ih=hakt
      dlw=2*ih+1
      clw=ih+1
      aws=lambda.lt.1d40
      n=n1*n2
      bii0=bi0
C   compute location weights first
      swj0=0.d0
      DO j2=1,dlw
         z2=clw-j2
         z2=z2*z2
         ih1=dsqrt(hakt2-z2)
         jind2=(j2-1)*dlw
         DO j1=clw-ih1,clw+ih1
C  first stochastic term
            jind=j1+jind2
            z1=clw-j1
            wj=lkern(kern,(z1*z1+z2)/hakt2)
            swj0=swj0+wj
            lw(jind)=wj
         END DO
      END DO
      bi0=swj0
      call rchkusr()
      DO i2=1,n2
         DO i1=1,n1
            iind=i1+(i2-1)*n1
            bii=bi(iind)/lambda
C   scaling of sij outside the loop
            swj=0.d0
            DO k=1,dv
               swjy(k)=0.d0
            END DO
            DO jw2=1,dlw
	       j2=jw2-clw+i2
	       if(j2.lt.1.or.j2.gt.n2) CYCLE
	       jind2=(j2-1)*n1
               jwind2=(jw2-1)*dlw
               z2=clw-jw2
               ih1=dsqrt(hakt2-z2*z2)
               DO jw1=clw-ih1,clw+ih1
		  j1=jw1-clw+i1
	          if(j1.lt.1.or.j1.gt.n1) CYCLE
		  jind=j1+jind2
		  wj=lw(jw1+jwind2)
                  IF (aws) THEN
                     sij=bii*kldistd(theta(iind),theta(jind),n,wght,dv)
                     IF (sij.gt.spmax) CYCLE
		     IF (skern.eq.1) THEN
C  skern == "Triangle"
                        wj=wj*(1.d0-sij)
		     ELSE
C  skern == "Exp"
		        IF (sij.gt.spmin) wj=wj*dexp(-spf*(sij-spmin))
		     ENDIF
                  END IF
                  swj=swj+wj
                  DO k=1,dv
                     swjy(k)=swjy(k)+wj*y(jind+(k-1)*n)
                  END DO
               END DO
            END DO
            DO k=1,dv
               thnew(iind+(k-1)*n)=swjy(k)/swj
            END DO
            bi(iind)=swj
            call rchkusr()
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsvimg(y,n1,n2,dv,vcoef,nvpar,meanvar,chcorr,
     1                   hakt,lambda,theta,bi,bi0,thnew,kern,skern,
     2                   spmin,spmax,wghts,lw,swjy)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew       \sum  Wi Y     (output)
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldistgc,lkern
      real*8 kldistgc,lkern
      integer n1,n2,dv,kern,skern,nvpar,y(1),theta(1),thnew(1)
      logical aws
      real*8 bi(1),lambda,spmax,spmin,hakt,lw(1),wghts(dv),bi0,
     2       vcoef(nvpar,dv),chcorr(1),meanvar(dv)
      integer ih,ih1,i1,i2,j1,j2,ja1,je1,l,k,n,info,i2n1,kdv,
     1        iind,jind,jind2,jwind2,dlw,clw,jw1,jw2,m0,thi(4)
      real*8 bii,sij,swj,swjy(dv),z1,z2,wj,hakt2,spf,thij(4),
     1       s2i(16),si(4),swj0
C  s2i, s2ii temporay stor sigma^2_i and its inverse (nneded for KL-distance)
C  maximaum dv = 4
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      ih=hakt
      dlw=2*ih+1
      clw=ih+1
      aws=lambda.lt.1d40
      n=n1*n2
C   compute location weights first
      swj0=0.d0
      DO j2=1,dlw
         z2=clw-j2
         z2=z2*z2
         ih1=dsqrt(hakt2-z2)
         ja1=max0(1,clw-ih1)
         je1=min0(dlw,clw+ih1)
         jind2=(j2-1)*dlw
         DO j1=ja1,je1
C  first stochastic term
            jind=j1+jind2
            z1=clw-j1
            wj=lkern(kern,(z1*z1+z2)/hakt2)
            swj0=swj0+wj
            lw(jind)=wj
         END DO
      END DO
      bi0=swj0
      call rchkusr()
      DO i2=1,n2
         i2n1=(i2-1)*n1
         DO i1=1,n1
            iind=i1+i2n1
            bii=bi(iind)/lambda
C   scaling of sij outside the loop
            swj=0.d0
            DO k=1,dv
               swjy(k)=0.d0
               thi(k)=theta(iind+(k-1)*n)
               si(k) = vcoef(1,k)
               if(nvpar.gt.1) THEN 
                  si(k) = si(k) + vcoef(2,k) * thi(k)
               END IF
               si(k) = dsqrt(dmax1(si(k),0.1*meanvar(k)))
C set small variances to  0.1 * mean variance
            END DO
C  Now fill estimated Covariancematrix in pixel i
            m0=1
            DO k=1,dv
               kdv = (k-1)*dv
               DO l=1,k
                  s2i(l+kdv)=si(k)*si(l)/wghts(k)/wghts(l)
                  if(l.ne.k) THEN
                     s2i(l+kdv)=s2i(l+kdv)*chcorr(m0)
                     m0=m0+1
                  END IF
               END DO
            END DO
            call dpotrf("U",dv,s2i,dv,info)
            IF (info.ne.0) call intpr("non-definite matrix 1",21,i,1)
	    call dpotri("U",dv,s2i,dv,info)
            IF (info.ne.0) call intpr("non-definite matrix 2",21,i,1)
            IF(dv.gt.1) THEN
               DO k=2,dv
                  kdv = (k-1)*dv
                  DO l=1,k-1
                     s2i(k+(l-1)*dv)=s2i(l+kdv)
                  END DO
               END DO
            END IF
            DO jw2=1,dlw
	       j2=jw2-clw+i2
	       if(j2.lt.1.or.j2.gt.n2) CYCLE
	       jind2=(j2-1)*n1
               jwind2=(jw2-1)*dlw
               z2=clw-jw2
               ih1=dsqrt(hakt2-z2*z2)
               DO jw1=clw-ih1,clw+ih1
		  j1=jw1-clw+i1
	          if(j1.lt.1.or.j1.gt.n1) CYCLE
		  jind=j1+jind2
                  DO k=1,dv
                     thij(k)=thi(k)-theta(jind+(k-1)*n)
                  END DO
		  wj=lw(jw1+jwind2)
                  IF (aws) THEN
                     sij=bii*kldistgc(thij,s2i,dv)
                     IF (sij.gt.spmax) CYCLE
		     IF (skern.eq.1) THEN
C  skern == "Triangle"
                        wj=wj*(1.d0-sij)
		     ELSE
C  skern == "Exp"
		        IF (sij.gt.spmin) wj=wj*dexp(-spf*(sij-spmin))
		     ENDIF
                  END IF
                  swj=swj+wj
                  DO k=1,dv
                     swjy(k)=swjy(k)+wj*y(jind+(k-1)*n)
                  END DO
               END DO
            END DO
            DO k=1,dv
               thnew(iind+(k-1)*n)=swjy(k)/swj
            END DO
            bi(iind)=swj
            call rchkusr()
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute Location Kernel (Compact support only, based on x^2
C                                   ignores scaling)
C
C          Kern=1     Uniform
C          Kern=2     Epanechnicov
C          Kern=3     Biweight
C          Kern=4     Triweight
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function lkern(kern,xsq)
      implicit logical (a-z)
      integer kern
      real*8 xsq,z
      IF (xsq.ge.1) THEN
         lkern=0.d0
      ELSE IF (kern.eq.1) THEN
         lkern=1.d0
      ELSE IF (kern.eq.2) THEN
         lkern=1.d0-xsq
      ELSE IF (kern.eq.3) THEN
         z=1.d0-xsq
         lkern=z*z
      ELSE IF (kern.eq.4) THEN
         z=1.d0-xsq
         lkern=z*z*z
      ELSE
C        use Epanechnikov
         lkern=1.d0-xsq
      ENDIF
      RETURN 
      END   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute the Kullback-Leibler Distance
C
C          Gaussian   
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function kldistgc(thij,s2ii,dv)
      implicit logical (a-z)
      integer k,l,dv
      real*8 thij(dv),s2ii(dv,dv),z,thijk
      z= thij(1)*thij(1)*s2ii(1,1)
      IF (dv.gt.1) THEN
         DO k=2,dv
            thijk=thij(k)
            DO l=1,k-1
               z = z + 2.d0*thijk*thij(l)*s2ii(l,k)
            END DO
            z = z + thijk*thijk*s2ii(k,k)
         END DO
      END IF
      kldistgc=z
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute the Kullback-Leibler Distance
C
C          Gaussian, Diagonal covariance matrix
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function kldistd(thi,thj,n,wght,nwght)
      implicit logical (a-z)
      integer n,nwght,i,k,thi(1),thj(1)
      real*8 z,wght(nwght)
      kldistd=0.d0
      i=1
      DO k=1,nwght
         z=thi(i)-thj(i)
         kldistd=kldistd+z*z*wght(k)
         i=i+n
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
C
C    Estimate variance parameters
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine esigmac(y,n,dv,theta,bi,quant,varcoef,mvar)
      implicit logical (a-z)
      integer n,dv,y(n,dv),theta(n,dv),quant(dv)
      real*8 bi(n),varcoef(dv),mvar(dv)
      integer i,k
      real*8 z,bii,sumres,sumwght,wght,res
      DO k=1,dv
         sumres=0.d0
         sumwght=0.d0
         DO i=1,n
            bii=bi(i)
            if(bii.le.1.d0.or.y(i,k).ge.quant(k)) CYCLE
            wght=bii-1.d0
            res=(y(i,k)-theta(i,k))
            res=res*res*bii/wght
            sumres=sumres+res*wght
            sumwght=sumwght+wght
         END DO
         z=sumres/sumwght
         varcoef(k)=z
         mvar(k)=z
      END DO
      RETURN
      END
      subroutine esigmal(y,n,dv,theta,bi,quant,varcoef,mvar)
      implicit logical (a-z)
      integer n,dv,y(n,dv),theta(n,dv),quant(dv)
      real*8 bi(n),varcoef(2,dv),mvar(dv),res
      integer i,k
      real*8 z,bii,s0,s1,s2,t0,t1,d,wght,thi,mth
      DO k=1,dv
         s0=0.d0
         s1=0.d0
         s2=0.d0
         t0=0.d0
         t1=0.d0
         mth=0.d0
         DO i=1,n
            bii=bi(i)
            mth=mth+theta(i,k)
            if(bii.le.1.d0.or.y(i,k).ge.quant(k)) CYCLE
            wght=bii-1.d0
            thi=theta(i,k)
            res=(y(i,k)-thi)
            res=res*res*bii/wght
            s0=s0+wght
            z=wght*thi
            s1=s1+z
            s2=s2+z*thi
            t0=t0+wght*res
            t1=t1+z*res
         END DO
         d=s2*s0-s1*s1
         varcoef(1,k)=(s2*t0-s1*t1)/d
         varcoef(2,k)=(-s1*t0+s0*t1)/d
         mvar(k)=varcoef(1,k)+varcoef(2,k)*mth/n
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
C
C    Estimate correlations
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine estcorr(res,n1,n2,dv,scorr,chcorr)
      implicit logical (a-z)
      integer n1,n2,dv
      real*8 res(n1,n2,dv),scorr(2,dv),chcorr(1)
      integer i,j,k,n,m,l
      real*8 vres(4),z,z1,z2,resij
      n=n1*n2
      DO k=1,dv
         z1=0.d0
         z2=0.d0
         DO i=1,n1
            DO j=1,n2
               resij=res(i,j,k)
               z1=z1+resij
               z2=z2+resij*resij
            END DO
         END DO
         z2=z2/n
         z1=z1/n
         vres(k)=n/(n-1)*(z2-z1*z1)
         z=0.d0
         DO i=1,n1
            DO j=1,n2
               res(i,j,k)=res(i,j,k)-z1
            END DO
         END DO
         DO i=1,n1-1
            DO j=1,n2
               z=z+res(i,j,k)*res(i+1,j,k)
            END DO
         END DO
         scorr(1,k)=z/n2/(n1-1)/vres(k)
         z=0.d0
         DO i=1,n1
            DO j=1,n2-1
               z=z+res(i,j,k)*res(i,j+1,k)
            END DO
         END DO
         scorr(2,k)=z/n1/(n2-1)/vres(k)
      END DO
C   between channels
      chcorr(1)=0.d0
      IF(dv.eq.1) RETURN
      m=1
      DO  k=1,dv-1
         DO l=k+1,dv
            z=0.d0
            DO i=1,n1
               DO j=1,n2
                  z=z+res(i,j,k)*res(i,j,l)
               END DO
            END DO
            chcorr(m)=z/n/dsqrt(vres(l)*vres(k))
            m=m+1
         END DO
      END DO
      RETURN
      END
