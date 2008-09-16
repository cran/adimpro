CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant  aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mawsimg(y,fix,mask,n1,n2,dv,hakt,lambda,theta,bi,
     1       bi0,thnew,kern,spmin,lw,wght,swjy)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew       \sum  Wi Y     (output)
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldistd,lkern
      real*8 kldistd,lkern
      integer n1,n2,dv,kern,y(n1,n2,dv),theta(n1,n2,dv),
     1        thnew(n1,n2,dv)
      logical aws,fix(n1,n2),mask(n1,n2)
      real*8 bi(n1,n2),bi0,lambda,spmin,wght(dv),hakt,lw(1)
      integer ih,ih1,i1,i2,j1,j2,k,n,
     1        jind,jind2,jwind2,dlw,clw,jw1,jw2
      real*8 bii,sij,swj,swj0,swjy(dv),z1,z2,wj,hakt2,bii0,spf
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
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
         ih1=sqrt(hakt2-z2)
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
C            iind=i1+(i2-1)*n1
            IF (fix(i1,i2)) CYCLE
C    nothing to do, final estimate is already fixed by control 
            bii=bi(i1,i2)/lambda
C   scaling of sij outside the loop
            swj=0.d0
            DO k=1,dv
               swjy(k)=0.d0
            END DO
            DO jw2=1,dlw
               j2=jw2-clw+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=(jw2-1)*dlw
               z2=clw-jw2
               ih1=sqrt(hakt2-z2*z2)
               DO jw1=clw-ih1,clw+ih1
                  j1=jw1-clw+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  if(.not.mask(j1,j2)) CYCLE
                  wj=lw(jw1+jwind2)
                  IF (aws) THEN
              sij=bii*kldistd(theta(i1,i2,1),theta(j1,j2,1),n,wght,dv)
                     IF (sij.gt.1.d0) CYCLE
                        wj=wj*(1.d0-sij)
                  END IF
                  swj=swj+wj
                  DO k=1,dv
                     swjy(k)=swjy(k)+wj*y(j1,j2,k)
                  END DO
               END DO
            END DO
            DO k=1,dv
               thnew(i1,i2,k)=swjy(k)/swj
            END DO
            bi(i1,i2)=swj
            call rchkusr()
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute nonadaptive kernel estimate
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsimg(y,n1,n2,dv,hakt,thnew,bi,kern,lw,swjy)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   bi       \sum  Wi   (output)
C   thnew    non-adaptive estimates    (output)
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      external kldistd,lkern
      real*8 kldistd,lkern
      integer n1,n2,dv,kern,y(n1,n2,dv),thnew(n1,n2,dv)
      real*8 bi(n1,n2),hakt,lw(1)
      integer ih,ih1,i1,i2,j1,j2,k,n,
     1        jind,jind2,jwind2,dlw,clw,jw1,jw2
      real*8 swj,swj0,swjy(dv),z1,z2,wj,hakt2
      hakt2=hakt*hakt
      ih=hakt
      dlw=2*ih+1
      clw=ih+1
      n=n1*n2
      swj0=0.d0
      DO j2=1,dlw
         z2=clw-j2
         z2=z2*z2
         ih1=sqrt(hakt2-z2)
         jind2=(j2-1)*dlw
         DO j1=clw-ih1,clw+ih1
            jind=j1+jind2
            z1=clw-j1
            wj=lkern(kern,(z1*z1+z2)/hakt2)
            swj0=swj0+wj
            lw(jind)=wj
         END DO
      END DO
      call rchkusr()
      DO i2=1,n2
         DO i1=1,n1
            swj=0.d0
            DO k=1,dv
               swjy(k)=0.d0
            END DO
            DO jw2=1,dlw
               j2=jw2-clw+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=(jw2-1)*dlw
               z2=clw-jw2
               ih1=sqrt(hakt2-z2*z2)
               DO jw1=clw-ih1,clw+ih1
                  j1=jw1-clw+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  wj=lw(jw1+jwind2)
                  swj=swj+wj
                  DO k=1,dv
                     swjy(k)=swjy(k)+wj*y(j1,j2,k)
                  END DO
               END DO
            END DO
            DO k=1,dv
               thnew(i1,i2,k)=swjy(k)/swj
            END DO
            bi(i1,i2)=swj
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
      subroutine awsvimg(y,fix,n1,n2,dv,vcoef,nvpar,meanvar,chcorr,
     1                   hakt,hhom,lambda,theta,bi,bi0,thnew,kern,
     2                   spmin,wghts,lw,swjy,early,homogen)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew       \sum  Wi Y     (output)
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      external kldistgc,lkern
      real*8 kldistgc,lkern
      integer n1,n2,dv,kern,nvpar,y(n1,n2,dv),theta(n1,n2,dv),
     1        thnew(n1,n2,dv)
      logical aws,fix(n1,n2),early,homogen,fixi
      real*8 bi(n1,n2),lambda,spmin,hakt,lw(1),wghts(dv),bi0,
     2       vcoef(nvpar,dv),chcorr(1),meanvar(dv),hhom(n1,n2)
      integer ih,ih1,i1,i2,j1,j2,ja1,je1,l,k,info,kdv,
     1        jind,jind2,jwind2,dlw,clw,jw1,jw2,m0,thi(4)
      real*8 bii,sij,swj,swjy(dv),z1,z2,wj,hakt2,spf,thij(4),
     1       s2i(16),si(4),swj0,hhomi,hhommax,hfixmax,hnfix,hmax2
C  s2i, s2ii temporay stor sigma^2_i and its inverse (nneded for KL-distance)
C  maximaum dv = 4
      hakt2=hakt*hakt
      hnfix=max(2.d0,6.d0-hakt)
      spf=1.d0/(1.d0-spmin)
      ih=hakt
      dlw=2*ih+1
      clw=ih+1
      aws=lambda.lt.1d40
C      n=n1*n2
C   compute location weights first
      swj0=0.d0
      hhomi=1.d0
      fixi=.FALSE.
      hmax2=0.d0
      DO j2=1,dlw
         z2=clw-j2
         z2=z2*z2
         ih1=sqrt(hakt2-z2)
         ja1=max(1,clw-ih1)
         je1=min(dlw,clw+ih1)
         jind2=(j2-1)*dlw
         DO j1=ja1,je1
C  first location weight
            jind=j1+jind2
            z1=clw-j1
            wj=lkern(kern,(z1*z1+z2)/hakt2)
            if(wj.gt.0) hmax2=max(hmax2,z1*z1+z2)
            swj0=swj0+wj
            lw(jind)=wj
         END DO
      END DO
      bi0=swj0
      call rchkusr()
      DO i2=1,n2
         DO i1=1,n1
            if(early) fixi=fix(i1,i2)
            if(fixi) THEN
               DO k=1,dv
               thnew(i1,i2,k)=theta(i1,i2,k)
               END DO
               CYCLE
            END IF
            if(homogen) THEN
               hhomi=hhom(i1,i2)
               hhomi=hhomi*hhomi
            END IF
            hhommax=hmax2
            hfixmax=hhomi
            bii=bi(i1,i2)/lambda
C   scaling of sij outside the loop
            swj=0.d0
            DO k=1,dv
               swjy(k)=0.d0
               thi(k)=theta(i1,i2,k)
               si(k) = vcoef(1,k)
               if(nvpar.gt.1) THEN 
                  si(k) = si(k) + vcoef(2,k) * thi(k)
               END IF
               if(nvpar.gt.2) THEN 
                  si(k) = si(k) + vcoef(3,k) * thi(k) * thi(k)
               END IF
               si(k) = sqrt(max(si(k),0.1*meanvar(k)))
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
         IF (info.ne.0) call intpr("non-definite matrix 1",21,info,1)
            call dpotri("U",dv,s2i,dv,info)
         IF (info.ne.0) call intpr("non-definite matrix 2",21,info,1)
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
               jwind2=(jw2-1)*dlw
               z2=clw-jw2
               z2=z2*z2
               ih1=sqrt(hakt2-z2)
               DO jw1=clw-ih1,clw+ih1
                  j1=jw1-clw+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  z1=clw-jw1
                  z1=z1*z1+z2
                  DO k=1,dv
                     thij(k)=thi(k)-theta(j1,j2,k)
                  END DO
                  wj=lw(jw1+jwind2)
                  IF (aws.and.z1.ge.hhomi) THEN
                     sij=bii*kldistgc(thij,s2i,dv)
                     IF (sij.gt.1.d0) THEN
                        if(homogen) hhommax=min(hhommax,z1)
                        CYCLE
                     END IF
                     if(early) hfixmax=max(hfixmax,z1)
                     IF (sij.gt.spmin) THEN
                         wj=wj*(1.d0-spf*(sij-spmin))
                         if(homogen) hhommax=min(hhommax,z1)
                     END IF 
                  END IF
                  swj=swj+wj
                  DO k=1,dv
                     swjy(k)=swjy(k)+wj*y(j1,j2,k)
                  END DO
               END DO
            END DO
            DO k=1,dv
               thnew(i1,i2,k)=swjy(k)/swj
            END DO
            bi(i1,i2)=swj
            if(homogen) hhom(i1,i2)=sqrt(hhommax)
            IF(early.and.hakt-sqrt(hfixmax).ge.hnfix) THEN
               fix(i1,i2)=.TRUE.
            END IF
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
         lkern=1.d0-xsq
      ELSE IF (kern.eq.2) THEN
         z=1.d0-xsq
         lkern=z*z
      ELSE IF (kern.eq.3) THEN
         z=1.d0-xsq
         lkern=z*z*z
      ELSE IF (kern.eq.4) THEN
C   Plateau
         IF(xsq.le.0.5d0) THEN
            lkern=1.d0
         ELSE
            lkern=2.d0*(1.d0-xsq)
         END IF
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
         if(sumwght.gt.0.d0) THEN
            z=sumres/sumwght
         ELSE
            z=1d-2
         END IF
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
         IF(d.gt.0.d0) THEN
            varcoef(1,k)=(s2*t0-s1*t1)/d
            varcoef(2,k)=(-s1*t0+s0*t1)/d
         ELSE
            varcoef(1,k)=1d-2
            varcoef(2,k)=0d0
         END IF
         mvar(k)=varcoef(1,k)+varcoef(2,k)*mth/n
      END DO
      RETURN
      END
      subroutine esigmaq(y,n,dv,theta,bi,quant,varcoef,mvar)
      implicit logical (a-z)
      integer n,dv,y(n,dv),theta(n,dv),quant(dv)
      real*8 bi(n),varcoef(3,dv),mvar(dv),res,mat(3,3),imat(3,3)
      integer i,k,info
      real*8 z,bii,s0,s1,s2,s3,s4,t0,t1,t2,wght,thi,mth,mthn,tt(3)
      DO k=1,dv
         s0=0.d0
         s1=0.d0
         s2=0.d0
         s3=0.d0
         s4=0.d0
         t0=0.d0
         t1=0.d0
         t2=0.d0
         mth=0.d0
         DO i=1,n
            bii=bi(i)
            if(theta(i,k).le.0.025*65535) CYCLE
            if(theta(i,k).gt.0.975*65535) CYCLE
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
            s3=s3+z*thi*thi
            s4=s4+z*thi*thi*thi
            t0=t0+wght*res
            t1=t1+z*res
            t2=t2+z*thi*res
         END DO
         mat(1,1)=s0
         mat(1,2)=s1
         mat(1,3)=s2
         mat(2,2)=s2
         mat(2,3)=s3
         mat(3,3)=s4
         imat(1,1)=1.d0
         imat(2,2)=1.d0
         imat(3,3)=1.d0
         imat(1,2)=0.d0
         imat(1,3)=0.d0
         imat(2,1)=0.d0
         imat(2,3)=0.d0
         imat(3,1)=0.d0
         imat(3,2)=0.d0
C     now calculate theta as B_i^{-1} A_i
         call dposv("U",3,3,mat,3,imat,3,info)
C    if info>0 just keep the old estimate
         IF (info.gt.0) THEN
             call intpr("info",4,info,1)
             varcoef(1,k)=1d-2
             varcoef(2,k)=0d0
             varcoef(3,k)=0d0
             mvar(k)=1d-2
             CYCLE  
         END IF 
         tt(1)=t0
         tt(2)=t1
         tt(3)=t2
         varcoef(1,k)=imat(1,1)*t0+imat(1,2)*t1+imat(1,3)*t2
         varcoef(2,k)=imat(2,1)*t0+imat(2,2)*t1+imat(2,3)*t2
         varcoef(3,k)=imat(3,1)*t0+imat(3,2)*t1+imat(3,3)*t2
         varcoef(3,k)=max(0.d0,varcoef(3,k))
         mthn=mth/n
         mvar(k)=varcoef(1,k)+varcoef(2,k)*mthn+varcoef(3,k)*mthn*mthn
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
C  just to avoid problems with images without noise !!!
         DO i=1,n1
            DO j=1,n2
               res(i,j,k)=res(i,j,k)-z1
            END DO
         END DO
         z=0.d0
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
            chcorr(m)=z/n/sqrt(vres(l)*vres(k))
            m=m+1
         END DO
      END DO
      RETURN
      END
