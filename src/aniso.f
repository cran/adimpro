      subroutine sm2dtens(x,n1,n2,h,rho,xhat)
      implicit none
      integer n1,n2
      double precision x(3,n1,n2),h,rho,xhat(3,n1,n2)
      integer i1,i2,j1,j2,ja1,je1,ja2,je2,jx,jy
      double precision a(3),b(3),sw,wij,adist2,h2,d,sqrd
      external adist2
      h2=h*h
      DO i1=1,n1
         DO i2=1,n2
            a(1)=x(1,i1,i2)*(1.d0+rho)
            a(2)=x(2,i1,i2)
            a(3)=x(3,i1,i2)*(1.d0+rho)
            d=a(1)*a(3)-a(2)*a(2)
            if(d.le.1e-15) THEN
               a(1)=1.d0
               a(2)=0.d0
               a(3)=1.d0
            ENDIF
            sqrd=sqrt(d)
            a(1)=a(1)/sqrd
            a(2)=a(2)/sqrd
            a(3)=a(3)/sqrd
            b(1)=0.d0
            b(2)=0.d0
            b(3)=0.d0
            sw=0.d0
            call rangex(a,h,ja1,je1)
            DO j1=ja1,je1
               jx=i1+j1
               if(jx.lt.1.or.jx.gt.n1) CYCLE
               call rangey(a,j1,h,ja2,je2)
               DO j2=ja2,je2
                  jy=i2+j2
                  if(jy.lt.1.or.jy.gt.n2) CYCLE
                  wij=1.d0-adist2(a,j1,j2)/h2
                  wij=max(0.d0,wij)
                  sw=sw+wij
                  b(1)=b(1)+wij*x(1,jx,jy)
                  b(2)=b(2)+wij*x(2,jx,jy)
                  b(3)=b(3)+wij*x(3,jx,jy)
              END DO
            END DO
            xhat(1,i1,i2)=b(1)/sw
            xhat(2,i1,i2)=b(2)/sw
            xhat(3,i1,i2)=b(3)/sw
         END DO
      END DO
      return
      end
      subroutine rangex(a,h,ja,je)
      implicit none
      integer ja,je
      double precision a(3),h,z
      z=sqrt(a(3))*h
      ja=int(-z)
      je=int(z)
      return
      end
      subroutine rangey(a,ix,h,ja,je)
      implicit none
      integer ix,ja,je
      double precision a(3),h,z1,z2
      z1=-a(2)/a(3)*ix
      z2=sqrt(a(3)*h*h-ix*ix)/a(3)
      ja=int(z1-z2)
      if(z1-z2.gt.0.d0) ja=ja+1
      je=int(z1+z2)
      if(z1+z2.lt.0.d0) je=je-1
C this accounts for rounding effects !!!
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute anisotropic distance
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function adist2(a,x,y)
C
      implicit none
      integer x,y
      double precision a(3)
      adist2=a(1)*x*x+2.d0*a(2)*x*y+a(3)*y*y
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant  aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine aniawsim(y,n1,n2,dv,ani,hakt,lambda,theta,bi,
     1       thnew,kern,skern,spmin,spmax,wght,swjy)
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
      implicit none
      external kldistd,lkern,adist2
      double precision kldistd,lkern,adist2
      integer n1,n2,dv,kern,skern,y(n1,n2,dv),theta(n1,n2,dv),
     1        thnew(n1,n2,dv)
      logical aws
      double precision bi(n1,n2),lambda,spmax,spmin,wght(dv),hakt,
     1        ani(3,n1,n2)
      integer i1,i2,j1,j2,k,n,jx,jy,ja1,ja2,je1,je2
      double precision bii,sij,swj,swjy(dv),wj,hakt2,spf,a(3),d,sqrd
      hakt2=hakt*hakt
C      spf=spmax/(spmax-spmin)
      spf=spmax/(spmax-spmin)
      aws=lambda.lt.1d40
      n=n1*n2
C   compute location weights first
      DO i2=1,n2
         DO i1=1,n1
C Characterize anisotropic neighborhood
            a(1)=ani(1,i1,i2)
            a(2)=ani(2,i1,i2)
            a(3)=ani(3,i1,i2)
            d=a(1)*a(3)-a(2)*a(2)
            if(d.le.1e-15) THEN
               a(1)=1.d0
               a(2)=0.d0
               a(3)=1.d0
            ENDIF
            sqrd=sqrt(d)
            a(1)=a(1)/sqrd
            a(2)=a(2)/sqrd
            a(3)=a(3)/sqrd
C now a descibrs an ellipse with unit volume
            bii=bi(i1,i2)/lambda
C   scaling of sij outside the loop
            swj=0.d0
            DO k=1,dv
               swjy(k)=0.d0
            END DO
            call rangex(a,hakt,ja1,je1)
            DO j1=ja1,je1
               jx=i1+j1
               if(jx.lt.1.or.jx.gt.n1) CYCLE
               call rangey(a,j1,hakt,ja2,je2)
               DO j2=ja2,je2
                  jy=i2+j2
                  if(jy.lt.1.or.jy.gt.n2) CYCLE
                  wj=lkern(kern,adist2(a,j1,j2)/hakt2)
                  IF (aws) THEN
                     sij=bii*kldistd(theta(i1,i2,1),theta(jx,jy,1),n,
     1                               wght,dv)
                     IF (sij.gt.spmax) CYCLE
                     IF (skern.eq.1) THEN
C  skern == "Triangle"
                        IF (sij.gt.spmin) wj=wj*spf*(1.d0-sij)
                     ELSE
C  skern == "Exp"
                     IF (sij.gt.spmin) wj=wj*dexp(-spf*(sij-spmin))
                     ENDIF
                  END IF
                  swj=swj+wj
                  DO k=1,dv
                     swjy(k)=swjy(k)+wj*y(jx,jy,k)
                  END DO
               END DO
            END DO
            DO k=1,dv
               thnew(i1,i2,k)=int(swjy(k)/swj)
            END DO
            bi(i1,i2)=swj
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine aniawsv(y,n1,n2,dv,ani,vcoef,nvpar,meanvar,chcorr,
     1                   hakt,lambda,theta,bi,thnew,kern,skern,
     2                   spmin,spmax,wghts,swjy)
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
      implicit none
      external kldistgc,lkern,adist2
      double precision kldistgc,lkern,adist2
      integer n1,n2,dv,kern,skern,nvpar,y(n1,n2,dv),theta(n1,n2,dv),
     1        thnew(n1,n2,dv)
      logical aws
      double precision bi(n1,n2),lambda,spmax,spmin,hakt,wghts(dv),
     2       vcoef(nvpar,dv),chcorr(*),meanvar(dv),ani(3,n1,n2)
      integer i1,i2,j1,j2,ja1,je1,ja2,je2,l,k,n,info,kdv,jx,jy,i,
     1        m0,thi(4)
      double precision bii,sij,swj,swjy(dv),wj,hakt2,spf,thij(4),
     1       s2i(16),si(4),sqrd,d,a(3)
C  s2i, s2ii temporay stor sigma^2_i and its inverse (nneded for KL-distance)
C  maximaum dv = 4
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      aws=lambda.lt.1d40
      n=n1*n2
C   compute location weights first
      DO i2=1,n2
         DO i1=1,n1
            a(1)=ani(1,i1,i2)
            a(2)=ani(2,i1,i2)
            a(3)=ani(3,i1,i2)
            d=a(1)*a(3)-a(2)*a(2)
            if(d.le.1e-15) THEN
               a(1)=1.d0
               a(2)=0.d0
               a(3)=1.d0
            ENDIF
            sqrd=sqrt(d)
            a(1)=a(1)/sqrd
            a(2)=a(2)/sqrd
            a(3)=a(3)/sqrd
C now a descibrs an ellipse with unit volume
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
C            IF (info.ne.0) call intpr("non-definite matrix 1",21,i,1)
            call dpotri("U",dv,s2i,dv,info)
C            IF (info.ne.0) call intpr("non-definite matrix 2",21,i,1)
            IF(dv.gt.1) THEN
               DO k=2,dv
                  kdv = (k-1)*dv
                  DO l=1,k-1
                     s2i(k+(l-1)*dv)=s2i(l+kdv)
                  END DO
               END DO
            END IF
            call rangex(a,hakt,ja1,je1)
            DO j1=ja1,je1
               jx=i1+j1
               if(jx.lt.1.or.jx.gt.n1) CYCLE
               call rangey(a,j1,hakt,ja2,je2)
               DO j2=ja2,je2
                  jy=i2+j2
                  if(jy.lt.1.or.jy.gt.n2) CYCLE
                  wj=lkern(kern,adist2(a,j1,j2)/hakt2)
                  DO k=1,dv
                     thij(k)=thi(k)-theta(jx,jy,k)
                  END DO
                  IF (aws) THEN
                     sij=bii*kldistgc(thij,s2i,dv)
                     IF (sij.gt.spmax) CYCLE
                     IF (skern.eq.1) THEN
C  skern == "Triangle"
                        if (sij.gt.spmin) wj=wj*spf*(1.d0-sij)
                     ELSE
C  skern == "Exp"
                     IF (sij.gt.spmin) wj=wj*dexp(-spf*(sij-spmin))
                     ENDIF
                  END IF
                  swj=swj+wj
                  DO k=1,dv
                     swjy(k)=swjy(k)+wj*y(jx,jy,k)
                  END DO
               END DO
            END DO
            DO k=1,dv
               thnew(i1,i2,k)=int(swjy(k)/swj)
            END DO
            bi(i1,i2)=swj
         END DO
      END DO
      RETURN
      END
