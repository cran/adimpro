CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant  aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mawsimg0(y,fix,mask,n1,n2,dv,hakt,lambda,theta,bi,
     1       bi0,thnew,kern,spmin,lw,wght)
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
      implicit none
      external kldistd,lkern
      double precision kldistd,lkern
      integer n1,n2,dv,kern,y(dv,n1,n2),theta(dv,n1,n2),
     1        thnew(dv,n1,n2),fix(n1,n2),mask(n1,n2)
      logical aws
      double precision bi(n1,n2),bi0,lambda,spmin,wght(dv),hakt,lw(*)
      integer ih,ih1,i1,i2,j1,j2,k,n,
     1        jind,jind2,jwind2,dlw,clw,jw1,jw2
      double precision bii,sij,swj,swj0,swjy(3),z1,z2,wj,hakt2,bii0,spf
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih=int(hakt)
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
         ih1=int(sqrt(hakt2-z2))
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
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,dv,kern,y,theta,thnew,fix,mask,bi,bi0,lambda,
C$OMP& spmin,wght,hakt,lw)
C$OMP& FIRSTPRIVATE(hakt2,dlw,clw,ih,aws,n,bii0,spf)
C$OMP& PRIVATE(ih1,i1,i2,j1,j2,k,jind,jind2,jwind2,jw1,jw2,
C$OMP& bii,sij,swj,swj0,z1,z2,wj,swjy)
C$OMP DO SCHEDULE(GUIDED)
      DO i2=1,n2
         DO i1=1,n1
C            iind=i1+(i2-1)*n1
            IF (fix(i1,i2).ne.0) CYCLE
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
               ih1=int(sqrt(hakt2-z2*z2))
               DO jw1=clw-ih1,clw+ih1
                  j1=jw1-clw+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  if(mask(j1,j2).eq.0) CYCLE
                  wj=lw(jw1+jwind2)
                  IF (aws) THEN
              sij=bii*kldistd(theta(1,i1,i2),theta(1,j1,j2),1,wght,dv)
                     IF (sij.gt.1.d0) CYCLE
                        wj=wj*(1.d0-sij)
                  END IF
                  swj=swj+wj
                  DO k=1,dv
                     swjy(k)=swjy(k)+wj*y(k,j1,j2)
                  END DO
               END DO
            END DO
            DO k=1,dv
               thnew(k,i1,i2)=int(swjy(k)/swj)
            END DO
            bi(i1,i2)=swj
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(bi,thnew)
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
      implicit none
      external kldistd,lkern
      double precision kldistd,lkern
      integer n1,n2,dv,kern,y(n1,n2,dv),thnew(n1,n2,dv)
      double precision bi(n1,n2),hakt,lw(*)
      integer ih,ih1,i1,i2,j1,j2,k,n,
     1        jind,jind2,jwind2,dlw,clw,jw1,jw2
      double precision swj,swj0,swjy(dv),z1,z2,wj,hakt2
      hakt2=hakt*hakt
      ih=int(hakt)
      dlw=2*ih+1
      clw=ih+1
      n=n1*n2
      swj0=0.d0
      DO j2=1,dlw
         z2=clw-j2
         z2=z2*z2
         ih1=int(sqrt(hakt2-z2))
         jind2=(j2-1)*dlw
         DO j1=clw-ih1,clw+ih1
            jind=j1+jind2
            z1=clw-j1
            wj=lkern(kern,(z1*z1+z2)/hakt2)
            swj0=swj0+wj
            lw(jind)=wj
         END DO
      END DO
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
               ih1=int(sqrt(hakt2-z2*z2))
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
               thnew(i1,i2,k)=int(swjy(k)/swj)
            END DO
            bi(i1,i2)=swj
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute nonadaptive kernel estimate
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsimg0(y,n1,n2,dv,hakt,thnew,bi,kern,lw,swjy)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   bi       \sum  Wi   (output)
C   thnew    non-adaptive estimates    (output)
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
C$    use omp_lib
      implicit none
      external kldistd,lkern
      double precision kldistd,lkern
      integer n1,n2,dv,kern,y(dv,n1,n2),thnew(dv,n1,n2)
      double precision bi(n1,n2),hakt,lw(*),swjy(dv,*)
      integer ih,ih1,ii,i1,i2,j1,j2,k,n,thrednr,
     1        jind,jind2,jwind2,dlw,clw,jw1,jw2
      double precision swj,swj0,z1,z2,wj,hakt2
      hakt2=hakt*hakt
      ih=FLOOR(hakt)
      dlw=2*ih+1
      clw=ih+1
      n=n1*n2
      swj0=0.d0
      DO j2=1,dlw
         z2=clw-j2
         z2=z2*z2
         ih1=FLOOR(sqrt(hakt2-z2))
         jind2=(j2-1)*dlw
         DO j1=clw-ih1,clw+ih1
            jind=j1+jind2
            z1=clw-j1
            wj=lkern(kern,(z1*z1+z2)/hakt2)
            swj0=swj0+wj
            lw(jind)=wj
         END DO
      END DO
      thrednr=1
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,dv,kern,y,thnew,bi,hakt,lw,swj0,swjy,clw,dlw,
C$OMP&        ih,n,hakt2,z1,jind,jind2)
C$OMP& PRIVATE(ii,ih1,i1,i2,j1,j2,k,jwind2,jw1,
C$OMP&         jw2,swj,z2,wj,thrednr)
C$OMP DO SCHEDULE(GUIDED)
      DO ii=1,n1*n2
C!$         thrednr = omp_get_thread_num()+1
         i1=mod(ii,n1)
         if(i1.eq.0) i1=n1
         i2=(ii-i1)/n1+1
         swj=0.d0
         DO k=1,dv
            swjy(k,thrednr)=0.d0
         END DO
         DO jw2=1,dlw
            j2=jw2-clw+i2
            if(j2.lt.1.or.j2.gt.n2) CYCLE
            jwind2=(jw2-1)*dlw
            z2=clw-jw2
            z2=z2*z2
            ih1=FLOOR(sqrt(hakt2-z2))
            DO jw1=clw-ih1,clw+ih1
               j1=jw1-clw+i1
               if(j1.lt.1.or.j1.gt.n1) CYCLE
               wj=lw(jw1+jwind2)
               swj=swj+wj
               DO k=1,dv
                  swjy(k,thrednr)=swjy(k,thrednr)+wj*y(k,j1,j2)
               END DO
            END DO
         END DO
         DO k=1,dv
            thnew(k,i1,i2)=FLOOR(swjy(k,thrednr)/swj)
         END DO
         bi(i1,i2)=swj
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(bi,thnew)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsvimg0(y,fix,n1,n2,n,dv,vcoef,nvpar,meanvar,chcorr,
     1                   hakt,hhom,lambda,theta,bi,bi0,thnew,kern,
     2                   spmin,wghts,lw,early,homogen)
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
      implicit none
      external kldistgc,lkern
      double precision kldistgc,lkern
      integer n,n1,n2,dv,kern,nvpar,y(dv,n1,n2),theta(dv,n1,n2),
     1        thnew(dv,n1,n2),fix(n),early,homogen
      logical aws
      double precision bi(n),lambda,spmin,hakt,lw(*),wghts(dv),bi0,
     1       vcoef(nvpar,dv),chcorr(*),meanvar(dv),hhom(n),
     2       thi(3),thij(3),si(3),s2i(9),swjy(3)
      integer ih,ih1,ii,i1,i2,j1,j2,ja1,je1,l,k,info,kdv,
     1        jind,jind2,jwind2,dlw,clw,jw1,jw2,m0,fixi
      double precision bii,sij,swj,z1,z2,wj,hakt2,spf,
     1       swj0,hhomi,hhommax,hfixmax,hnfix,hmax2
C  s2i, s2ii temporay stor sigma^2_i and its inverse (nneded for KL-distance)
C  maximaum dv = 4
      hakt2=hakt*hakt
      hnfix=max(2.d0,6.d0-hakt)
      spf=1.d0/(1.d0-spmin)
      ih=int(hakt)
      dlw=2*ih+1
      clw=ih+1
      aws=lambda.lt.1d40
C      n=n1*n2
C   compute location weights first
      swj0=0.d0
      fixi=0
      hmax2=0.d0
      DO j2=1,dlw
         z2=clw-j2
         z2=z2*z2
         ih1=int(sqrt(hakt2-z2))
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
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n1,n2,dv,kern,nvpar,y,theta,thnew,fix,early,homogen,
C$OMP&        bi,lambda,spmin,hakt,lw,wghts,vcoef,chcorr,meanvar,
C$OMP&        hhom,bi0,dlw,clw,n)
C$OMP& FIRSTPRIVATE(aws,ih,hnfix,hakt2,spf,swj0,fixi,hmax2)
C$OMP& PRIVATE(ii,ih1,i1,i2,j1,j2,ja1,je1,l,k,info,kdv,jind,jind2,
C$OMP&         jwind2,jw1,jw2,m0,bii,sij,swj,z1,z2,wj,hhommax,hfixmax,
C$OMP&         hhomi,thi,thij,si,s2i,swjy)
C$OMP DO SCHEDULE(GUIDED)
      DO ii=1,n1*n2
         hhomi=1.d0
         i1=mod(ii,n1)
         if(i1.eq.0) i1=n1
         i2=(ii-i1)/n1+1
         if(early.ne.0) fixi=fix(ii)
         if(fixi.ne.0) THEN
            DO k=1,dv
               thnew(k,i1,i2)=theta(k,i1,i2)
            END DO
            CYCLE
         END IF
         if(homogen.ne.0) THEN
            hhomi=hhom(ii)
            hhomi=hhomi*hhomi
         END IF
         hhommax=hmax2
         hfixmax=hhomi
         bii=bi(ii)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         DO k=1,dv
            swjy(k)=0.d0
            thi(k)=theta(k,i1,i2)
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
         call dpotrf("U",dv,s2i(1),dv,info)
C         IF (info.ne.0) call intpr("non-definite matrix 1",21,info,1)
         call dpotri("U",dv,s2i(1),dv,info)
C         IF (info.ne.0) call intpr("non-definite matrix 2",21,info,1)
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
            ih1=int(sqrt(hakt2-z2))
            DO jw1=clw-ih1,clw+ih1
               j1=jw1-clw+i1
               if(j1.lt.1.or.j1.gt.n1) CYCLE
               z1=clw-jw1
               z1=z1*z1+z2
               DO k=1,dv
                  thij(k)=thi(k)-theta(k,j1,j2)
               END DO
               wj=lw(jw1+jwind2)
               IF (aws.and.z1.ge.hhomi) THEN
                  sij=bii*kldistgc(thij,s2i,dv)
                  IF (sij.gt.1.d0) THEN
                     if(homogen.ne.0) hhommax=min(hhommax,z1)
                     CYCLE
                  END IF
                  if(early.ne.0) hfixmax=max(hfixmax,z1)
                  IF (sij.gt.spmin) THEN
                     wj=wj*(1.d0-spf*(sij-spmin))
                     if(homogen.ne.0) hhommax=min(hhommax,z1)
                  END IF
               END IF
               swj=swj+wj
               DO k=1,dv
                  swjy(k)=swjy(k)+wj*y(k,j1,j2)
               END DO
            END DO
         END DO
         DO k=1,dv
            thnew(k,i1,i2)=int(swjy(k)/swj)
         END DO
         bi(ii)=swj
         if(homogen.ne.0) hhom(ii)=sqrt(hhommax)
         IF(early.ne.0.and.hakt-sqrt(hfixmax).ge.hnfix) THEN
            fix(ii)=1
         END IF
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(bi,thnew,fix,hhom)
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
      double precision function lkern(kern,xsq)
      implicit none
      integer kern
      double precision xsq,z
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
      double precision function kldistgc(thij,s2ii,dv)
      implicit none
      integer k,l,dv
      double precision thij(dv),s2ii(dv,dv),z,thijk
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
      double precision function kldistd(thi,thj,n,wght,nwght)
      implicit none
      integer n,nwght,i,k,thi(*),thj(*)
      double precision z,wght(nwght)
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
      implicit none
      integer n,dv,y(n,dv),theta(n,dv),quant(dv)
      double precision bi(n),varcoef(dv),mvar(dv)
      integer i,k
      double precision z,bii,sumres,sumwght,wght,res
      DO k=1,dv
         sumres=0.d0
         sumwght=0.d0
C$OMP PARALLEL DO
C$OMP& DEFAULT(SHARED)
C$OMP& PRIVATE(i,bii,res,wght)
C$OMP& REDUCTION(+:sumres,sumwght)
         DO i=1,n
            bii=bi(i)
            if(bii.le.1.d0.or.y(i,k).ge.quant(k)) CYCLE
            wght=bii-1.d0
            res=(y(i,k)-theta(i,k))
            res=res*res*bii/wght
            sumres=sumres+res*wght
            sumwght=sumwght+wght
         END DO
C$OMP END PARALLEL DO
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
      implicit none
      integer n,dv,y(n,dv),theta(n,dv),quant(dv)
      double precision bi(n),varcoef(2,dv),mvar(dv),res
      integer i,k
      double precision z,bii,s0,s1,s2,t0,t1,d,wght,thi,mth
      DO k=1,dv
         s0=0.d0
         s1=0.d0
         s2=0.d0
         t0=0.d0
         t1=0.d0
         mth=0.d0
C$OMP PARALLEL DO
C$OMP& DEFAULT(SHARED)
C$OMP& PRIVATE(i,bii,thi,res,wght,z)
C$OMP& REDUCTION(+:s0,s1,s2,t0,t1,mth)
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
C$OMP END PARALLEL DO
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
      implicit none
      integer n,dv,y(n,dv),theta(n,dv),quant(dv)
      double precision bi(n),varcoef(3,dv),mvar(dv),res,mat(3,3),
     1        imat(3,3)
      integer i,k,info
      double precision z,bii,s0,s1,s2,s3,s4,t0,t1,t2,wght,thi,mth,
     1        mthn,tt(3)
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
C$OMP PARALLEL DO
C$OMP& DEFAULT(SHARED)
C$OMP& PRIVATE(i,bii,thi,res,wght,z)
C$OMP& REDUCTION(+:s0,s1,s2,s3,s4,t0,t1,t2,mth)
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
C$OMP END PARALLEL DO
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
C             call intpr("info",4,info,1)
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
      subroutine estcorr(res,n1,n2,dv,scorr,chcorr)
      implicit none
      integer n1,n2,dv
      double precision res(n1,n2,dv),scorr(2,dv),chcorr(1)
      integer i,j,k,n,m,l
      double precision vres(4),z,z1,z2,resij
      n=n1*n2
      DO k=1,dv
         z1=0.d0
         z2=0.d0
C$OMP PARALLEL DO
C$OMP& DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,resij)
C$OMP& REDUCTION(+:z1,z2)
         DO j=1,n2
            DO i=1,n1
               resij=res(i,j,k)
               z1=z1+resij
               z2=z2+resij*resij
            END DO
         END DO
C$OMP END PARALLEL DO
         z2=z2/n
         z1=z1/n
         vres(k)=n/(n-1)*(z2-z1*z1)
C  just to avoid problems with images without noise !!!
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j)
C$OMP DO SCHEDULE(GUIDED)
         DO j=1,n2
            DO i=1,n1
               res(i,j,k)=res(i,j,k)-z1
            END DO
         END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(res)
         z=0.d0
C$OMP PARALLEL DO
C$OMP& DEFAULT(SHARED)
C$OMP& PRIVATE(i,j)
C$OMP& REDUCTION(+:z)
         DO j=1,n2
            DO i=1,n1-1
               z=z+res(i,j,k)*res(i+1,j,k)
            END DO
         END DO
C$OMP END PARALLEL DO
         scorr(1,k)=z/n2/(n1-1)/vres(k)
         z=0.d0
C$OMP PARALLEL DO
C$OMP& DEFAULT(SHARED)
C$OMP& PRIVATE(i,j)
C$OMP& REDUCTION(+:z)
         DO j=1,n2-1
            DO i=1,n1
               z=z+res(i,j,k)*res(i,j+1,k)
            END DO
         END DO
C$OMP END PARALLEL DO
         scorr(2,k)=z/n1/(n2-1)/vres(k)
      END DO
C   between channels
      chcorr(1)=0.d0
      IF(dv.eq.1) RETURN
      m=1
      DO  k=1,dv-1
         DO l=k+1,dv
            z=0.d0
C$OMP PARALLEL DO
C$OMP& DEFAULT(SHARED)
C$OMP& PRIVATE(i,j)
C$OMP& REDUCTION(+:z)
            DO j=1,n2
               DO i=1,n1
                  z=z+res(i,j,k)*res(i,j,l)
               END DO
            END DO
C$OMP END PARALLEL DO
            chcorr(m)=z/n/sqrt(vres(l)*vres(k))
            m=m+1
         END DO
      END DO
      RETURN
      END
