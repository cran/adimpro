CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant  aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awspimg(y,n1,n2,dv,degr,hw,vcoef,nv,mvar,hakt,lambda,
     1           theta,bi,bi0,ai,kern,spmin,lw,w,slw,sw,
     2           wght,ind)
C
C   y        observed values of regression function
C   fix      logical TRUE fro points where we have nothing to do
C   n1,n2    design dimensions
C   dv       1 (grayscale) or 3 (color) numer of channels
C   degr     degree of polynomials 0,1 or 2
C   hw       bandwidth used to smooth weights
C   hakt     actual bandwidth in aws
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       Matrix Bi dim(n1,n2,dp2)
C   bi0      Matrix Bi0 dim(n1,n2,dp2) (with location weights only)
C   ai       \sum  Wi Y     (output) dim(n1,n2,dp1,dv)
C   kern     specifies the location kernel
C   lw       array of location weights dim(dlw,dlw) dlw=2*ih+1
C   w        array of weights dim(dlw,dlw)
C   sw       array of "smoothed" weights dim(dls,dls) dls=2*(ih+ihw)+1
C   wght     weight for color channels
C
C   temporary arrays set for maximum degree 2
C
      implicit none
      external kldistp2,lkern
      double precision kldistp2,lkern
      integer n1,n2,dv,kern,nv,degr,ind(*),y(*)
      logical aws
      double precision theta(*),bi(*),bi0,ai(*),lambda,spmin,mvar(dv),
     1       wght(dv),hakt,lw(*),w(*),hw,sw(*),slw(*),vcoef(nv,dv)
      integer ih,ih1,i1,i2,j1,j2,k,n,km1dp1,
     1        iind,jind,jind2,jwind,jwind2,dlw,clw,jw1,jw2,
     2        dp1,dp2,dnp1,ihs,csw,dsw,kw,l,lk,dsw2
      double precision bii(15),sij,swj(15),swjy(18),z1,z2,wj,swj0,si,
     1  thi(18),hakt2,thij(18),zz(15),lwj,hs2,hs,z,cc,wjy,spf,wghti(4)
C   arrays with variable length are organized as
C   theta(n1,n2,dp1,dv)
C   bi(n1,n2,dp2)
C   arrays of fixed length correspond to degr=2
C   first set dimensions for arrays depending on degree
      aws=lambda.lt.1.d20
      spf=1.d0/(1.d0-spmin)
      if(degr.eq.0) THEN
         dp1=1
         dp2=1
      ELSE IF (degr.eq.1) THEN
         dp1=3
         dp2=6
      ELSE
         dp1=6
         dp2=15
      END IF
      dnp1=dv*dp1
      hakt2=hakt*hakt
      ih=int(hakt)
      dlw=2*ih+1
      clw=ih+1
      hs=hakt+hw
      hs2=hs*hs
      ihs=int(hs)
      dsw=2*ihs+1
      csw=ihs+1
      n=n1*n2
C   compute location weights first  sum in slw
      DO j2=1,dlw
         z2=j2-clw
         z2=z2*z2
         ih1=int(sqrt(hakt2-z2))
         jind2=(j2-1)*dlw
         DO j1=clw-ih1,clw+ih1
C  first stochastic term
            jind=j1+jind2
            z1=j1-clw
            lw(jind)=lkern(kern,(z1*z1+z2)/hakt2)
         END DO
      END DO
      cc=0.0d0
      call smwghts2(lw,hakt,hw,slw,dlw,dsw,cc)
      swj0=0.d0
      DO j2=1,dsw
         jind2=(j2-1)*dsw
         DO j1=1,dsw
            swj0=swj0+slw(j1+jind2)
         END DO
      END DO
      bi0=swj0
C  now stochastic term
      zz(1)=1.d0
      DO i2=1,n2
         DO i1=1,n1
            iind=i1+(i2-1)*n1
            DO k=1,dv
               thi(k)=theta(iind+(k-1)*n*dp1)
               si = vcoef(1,k)
               if(nv.gt.1) si = si + vcoef(2,k) * thi(k)
               si = max(si,0.1*mvar(k))
               wghti(k)=wght(k)/si
C               call dblepr("si",2,si,1)
            END DO
            DO k=1,dp2
               bii(k)=bi(iind+(k-1)*n)/lambda
            END DO
            DO k=1,dnp1
               thi(k)=theta(iind+(k-1)*n)
            END DO
C   scaling of sij outside the loop
            DO jw2=1,dlw
               jwind2=(jw2-1)*dlw
               DO jw1=1,dlw
                  w(jw1+jwind2)=0.d0
               END DO
            END DO
            DO jw2=1,dlw
               j2=jw2-clw+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1
               jwind2=(jw2-1)*dlw
               z2=jw2-clw
C  get directional differences that only depend on i2-j2
               IF(dp1.gt.1) THEN
                  zz(3)=z2
                  zz(6)=z2*z2
               END IF
               ih1=int(sqrt(hakt2-z2*z2))
               DO jw1=clw-ih1,clw+ih1
                  j1=jw1-clw+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  jwind=jw1+jwind2
                  wj=lw(jwind)
                  z1=jw1-clw
C  get rest of directional differences
                  IF(dp1.gt.1) THEN
                     zz(2)=z1
                     zz(4)=z1*z1
                     zz(5)=z1*z2
                  END IF
                  IF (aws) THEN
                     DO k=1,dnp1
                        thij(k)=theta(jind+(k-1)*n)
                     END DO
                     IF(dp1.gt.1) THEN
                        DO k=1,dv
                           kw=(k-1)*dp1
                           thij(1+kw)=thij(1+kw)-thij(2+kw)*z1-
     1                              thij(3+kw)*z2
                           IF (dp1.gt.3) THEN
                              thij(1+kw)=thij(1+kw)+thij(4+kw)*zz(4)+
     1                           thij(5+kw)*zz(5)+thij(6+kw)*zz(6)
                              thij(2+kw)=thij(2+kw)-thij(5+kw)*z2-
     1                                 2.d0*thij(4+kw)*z1
                              thij(3+kw)=thij(3+kw)-thij(5+kw)*z1-
     1                                 2.d0*thij(6+kw)*z2
                           END IF
                        END DO
                     END IF
C
C           get difference of thetas
C
                     DO k=1,dnp1
                        thij(k)=thi(k)-thij(k)
                     END DO
                   sij=kldistp2(dp1,thij,bii,wghti,dv,ind)
                     w(jwind)=0.d0
                     IF (sij.gt.1.d0) CYCLE
                        IF(sij.gt.spmin) THEN
                           w(jwind)=wj*spf*(1.d0-sij)
                        ELSE
                           w(jwind)=wj
                        END IF
                  ELSE
                     w(jwind)=wj
                  END IF
               END DO
            END DO
C
C      Smooth the weights
C
            call testwght(w,dlw,dp1,hw,z)
            z=max(.1d0,min(z,hw))
            cc=min(z-1.d0,1.d0/hakt2)
            call smwghts2(w,hakt,z,sw,dlw,dsw,cc)
            DO k=1,dp2
               swj(k)=0.d0
            END DO
            DO k=1,dnp1
               swjy(k)=0.d0
            END DO
            dsw2=dsw*dsw
            DO jw2=1,dsw
               j2=jw2-csw+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1
               jwind2=(jw2-1)*dsw
               z2=jw2-csw
               IF(dp1.gt.1) THEN
                  zz(3)=z2
                  zz(6)=z2*z2
               END IF
               IF(dp1.gt.3) THEN
                  zz(10)=z2*zz(6)
                  zz(15)=z2*zz(10)
               END IF
               ih1=int(sqrt(hs2-z2*z2))
               DO jw1=csw-ih1,csw+ih1
                  j1=jw1-csw+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  z1=jw1-csw
                  jwind=jw1+jwind2
                  jind=j1+jind2
                  lwj=slw(jwind)
                  wj=sw(jwind)
                  if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE
                  IF(dp1.gt.1) THEN
                     zz(2)=z1
                     zz(4)=z1*z1
                     zz(5)=z1*z2
                  END IF
                  IF(dp1.gt.3) THEN
                     zz(7)=z1*zz(4)
                     zz(8)=z1*zz(5)
                     zz(9)=z1*zz(6)
                     zz(11)=z1*zz(7)
                     zz(12)=z1*zz(8)
                     zz(13)=z1*zz(9)
                     zz(14)=z1*zz(10)
                  END IF
                  if(wj.le.0.d0) CYCLE
                  DO k=1,dp2
                     swj(k)=swj(k)+wj*zz(k)
                  END DO
                  DO k=1,dv
                     wjy=wj*y(jind+(k-1)*n)
                     km1dp1=(k-1)*dp1
                     DO l=1,dp1
                        lk=l+km1dp1
                        swjy(lk)=swjy(lk)+wjy*zz(l)
                     END DO
                  END DO
               END DO
            END DO
            DO k=1,dnp1
               ai(iind+(k-1)*n)=swjy(k)
            END DO
            DO k=1,dp2
               bi(iind+(k-1)*n)=swj(k)
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function kldistp2(dp1,thij,bii,wght,dv,ind)
C
C  search for maximum in w within bandwidth hw, result in sw
C
C     dp1  polynomial degree +1
C     thij parameter estimate in j dim(dp1*nwght) for basis in i
C     bii  XTX dim(dp2)
C     wght     weight for color channels
C     nwght    number of positive weights (<=dv)
C     ind   index matrix to access the correct elements in bii
C
      implicit none
      integer dp1,dv,ind(dp1,dp1)
      double precision thij(*),bii(*),wght(dv),thijl
      integer i,j,l,k
      double precision z,d
      z=0.d0
      DO i=1,dv
         d=0.d0
         j=(i-1)*dp1
         DO l=1,dp1
            thijl=thij(l+j)
            d=d+bii(ind(l,l))*thijl*thijl
            IF(l.eq.dp1) CYCLE
            DO k=l+1,dp1
               d=d+2.d0*bii(ind(k,l))*thijl*thij(k+j)
            END DO
         END DO
         z=z+d*wght(i)
      END DO
      kldistp2=z
      RETURN
      END
      subroutine smwghts2(w,hakt,hw,sw,dw,dsw,cc)
C
C  smooth w with epakern and bandwidth hw, result in sw
C
C     w  array of weights dim(dw,dw)   dw=2*ih+1
C     hakt aktual bandwidth in w
C     hw   bandwidth for smoothing of w
C     sw   array of smoothed weights dim(dsw,dsw)   dsw=2*(ihw+ih)+1
C     cc   dumping factor of weights
C
      implicit none
      integer dw,dsw,cw,csw,cdiff
      double precision w(dw,dw),sw(dsw,dsw),hw,hakt,cc
      integer i1,i2,id,jd,ja1,je1,ja2,je2,j1,j2,i10,i20
      double precision z,z0,z1,z2,hw2,zmax,hakt2,hsw,hsw2,ww
      cw=(dw+1)/2
      csw=(dsw+1)/2
      cdiff=csw-cw
      hsw=hw+hakt
      hsw2=hsw*hsw
      hakt2=hakt*hakt
      hw2=hw*hw
      zmax=0.d0
      DO i1=1,dsw
         DO i2=1,dsw
            sw(i1,i2)=0.d0
         END DO
      END DO
      IF(cc.le.0.d0) THEN
         DO j1=1,dw
            DO j2=1,dw
               sw(j1+cdiff,j2+cdiff)=w(j1,j2)
            END DO
         END DO
      ELSE
         DO i1=1,dsw
            z1=i1-csw
            i10=i1-cdiff
            ja1=max(i1-2*cdiff,1)
            je1=min(i1,dw)
            id=int(sqrt(hsw2-z1*z1))
            if(csw-id.lt.1) CYCLE
            DO i2=csw-id,csw+id
               i20=i2-cdiff
               z=0.d0
               z0=0.d0
               DO j1=ja1,je1
                  z1=(i10-j1)
                  z1=z1*z1
                  if(hw2-z1.lt.0.d0) CYCLE
                  jd=int(sqrt(hw2-z1))
                  ja2=max(i20-jd,1)
                  je2=min(i20+jd,dw)
                  DO j2=ja2,je2
                     z2=(i20-j2)
                     ww=(1.d0-(z1+z2*z2)/hw2)
                     if(ww.lt.1.d0) ww=cc*ww
                     z=z+ww*w(j1,j2)
                  END DO
               END DO
               sw(i1,i2)=z
               zmax=max(zmax,z)
            END DO
         END DO
         DO i1=1,dsw
            DO i2=1,dsw
               sw(i1,i2)=sw(i1,i2)/zmax
            END DO
         END DO
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      test regulatity of weight matrix
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine testwght(w,dlw,dp1,hw,z)
      integer dlw,dp1
      double precision w(dlw,dlw),hw,z
      integer clw,i,ip,im,cp1,cp2,cm1,cm2
      double precision zh,zv
      clw=(dlw+1)/2
      cp1=clw+1
      cm1=clw-1
      z=hw
      IF(clw.gt.2.and.dp1.eq.3) THEN
         cp2=clw+2
         cm2=clw-2
         zh=w(clw,cp1)*w(clw,cp2)+w(clw,cm1)*w(clw,cm2)
         zv=w(cp1,clw)*w(cp2,clw)+w(cm1,clw)*w(cm2,clw)
         IF(zh*zv.gt.0.125d0) THEN
            z=hw-2.d0
         ELSE
            DO i=1,clw-1
               ip=clw+i
               im=clw-i
               zh=zh+w(ip,cp1)*w(ip,cp2)+w(ip,cm1)*w(ip,cm2)+
     1               w(im,cp1)*w(im,cp2)+w(im,cm1)*w(im,cm2)
               zv=zv+w(cp1,ip)*w(cp2,ip)+w(cm1,ip)*w(cm2,ip)+
     1               w(cp1,im)*w(cp2,im)+w(cm1,im)*w(cm2,im)
               IF(zh*zv.gt.0.125d0) THEN
                  z=hw-2.d0
                  CYCLE
               END IF
            END DO
            IF(zh*zv.le.0.125d0) THEN
               zh=w(clw,cp1)+w(clw,cm1)
               zv=w(cp1,clw)+w(cm1,clw)
               DO i=1,clw-1
                  ip=clw+i
                  im=clw-i
                  zh=zh+w(ip,cp1)+w(ip,cm1)+w(im,cp1)+w(im,cm1)
                  zv=zv+w(cp1,ip)+w(cm1,ip)+w(im,cp1)+w(im,cm1)
                  IF(zh*zv.gt.0.5d0) THEN
                     z=hw-1.d0
                     CYCLE
                  END IF
               END DO
            END IF
         END IF
      ELSE IF(clw.gt.1.and.dp1.eq.2) THEN
         zh=w(clw,cp1)+w(clw,cm1)
         zv=w(cp1,clw)+w(cm1,clw)
         IF(zh*zv.gt.0.5d0) THEN
            z=hw-1.d0
         ELSE
            DO i=1,clw-1
               ip=clw+i
               im=clw-i
               zh=zh+w(ip,cp1)+w(ip,cm1)+w(im,cp1)+w(im,cm1)
               zv=zv+w(cp1,ip)+w(cm1,ip)+w(im,cp1)+w(im,cm1)
               IF(zh*zv.gt.0.5d0) THEN
                  z=hw-1.d0
                  CYCLE
               END IF
            END DO
         END IF
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (bivariate case)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mpaws2(n,dp1,dp2,ai,bi,theta,dmat,ind)
C
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     ai         \sum \Psi^T Wi^k Y
C     bi         \sum \Psi^T Wi^k \Psi
C     theta      new parameter estimate
C     dmat       working arrays
C     restricted to dp2<=20
      implicit none
      integer n,dp1,dp2
      double precision ai(n,dp1),bi(n,dp2),theta(n,dp1),dmat(dp1,dp1)
      integer i,j,k,info,ind(dp1,dp1)
      double precision aa(20)
      DO i=1,n
         DO k=1,dp1
            DO j=k,dp1
               dmat(k,j)=bi(i,ind(k,j))
            END DO
            aa(k)=ai(i,k)
         END DO
C     now calculate theta as B_i^{-1} A_i
         call dposv("U",dp1,1,dmat,dp1,aa,dp1,info)
C    if info>0 just keep the old estimate
         IF (info.gt.0) CYCLE
         DO j=1,dp1
            theta(i,j)=aa(j)
         END DO
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
      subroutine epsigmac(y,n,dv,theta,bi,quant,varcoef,mvar,dp1)
      implicit none
      integer n,dv,y(n,dv),theta(n,dv),quant(dv),dp1
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
            if(bii.le.dp1.or.y(i,k).ge.quant(k)) CYCLE
            res=(y(i,k)-theta(i,k))
            wght=bii-dp1
            res=res*res*bii/(bii-dp1)
            sumres=sumres+res*wght
            sumwght=sumwght+wght
         END DO
C$OMP END PARALLEL DO
         z=sumres/sumwght
         varcoef(k)=z
         mvar(k)=z
      END DO
      RETURN
      END
      subroutine epsigmal(y,n,dv,theta,bi,quant,varcoef,mvar,dp1)
      implicit none
      integer n,dv,y(n,dv),theta(n,dv),quant(dv),dp1
      double precision bi(n),varcoef(2,dv),
     1       mvar(dv),res
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
            thi=theta(i,k)
            mth=mth+thi
            if(bii.le.dp1.or.y(i,k).ge.quant(k)) CYCLE
            wght=bii-dp1
            res=(y(i,k)-thi)
            res=res*res*bii/(bii-dp1)
            s0=s0+wght
            z=wght*thi
            s1=s1+z
            s2=s2+z*thi
            t0=t0+wght*res
            t1=t1+z*res
         END DO
C$OMP END PARALLEL DO
         d=s2*s0-s1*s1
         varcoef(1,k)=(s2*t0-s1*t1)/d
         varcoef(2,k)=(-s1*t0+s0*t1)/d
         mvar(k)=varcoef(1,k)+varcoef(2,k)*mth/n
      END DO
      RETURN
      END
