      subroutine segment(y,fix,level,delta,n1,n2,hakt,lambda,
     1                   theta,vcoef,nvc,meanvar,bi,s2i,thnew,kern,
     2                   spmin,lw,pvalue,segm,thresh,ext,fov,varest)
      implicit logical (a-z)
      external kldistgc,lkern,fpchisq
      real*8 kldistgc,lkern,fpchisq
      integer n1,n2,y(n1,n2),theta(n1,n2),thnew(n1,n2),kern,
     1        segm(n1,n2),nvc
      logical fix(n1,n2)
      real*8 level,delta,hakt,lambda,bi(n1,n2),spmin,lw(1),ext,
     1       pvalue(n1,n2),thresh,vcoef(nvc),meanvar,s2i(n1,n2),
     1       varest(n1,n2),fov
      integer ih,dlw,clw,ih1,ja1,je1,jind2,jind,j1,j2,i1,i2,jw2,jw1,
     1        jwind2
      real*8 hakt2,spf,swj0,swj,swjy,z2,z1,wj,thi,thij,si,bii,sij,
     1        a,b,ti,swjw,swjw2,cofh
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih=hakt
      dlw=2*ih+1
      clw=ih+1
      aws=lambda.lt.1d40
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih=hakt
      dlw=2*ih+1
      clw=ih+1
      aws=lambda.lt.1d40
      swj0=0.d0
      DO j2=1,dlw
         z2=clw-j2
         z2=z2*z2
         ih1=sqrt(hakt2-z2)
         ja1=max(1,clw-ih1)
         je1=min(dlw,clw+ih1)
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
      call rchkusr()
C  fill field of pvalues
      a = level-delta
      b = level+delta
      DO i1=1,n1
         DO i2=1,n2
            if(fix(i1,i2)) CYCLE
            thi = theta(i1,i2)
            si = vcoef(1)
            if(nvc.gt.1) THEN 
               si = si + vcoef(2) * thi
            END IF
            if(nvc.gt.2) THEN 
               si = si + vcoef(3) * thi * thi
            END IF
            si = max(si,0.1*meanvar)
C set small variances to  0.1 * mean variance
C  Now fill estimated Covariancematrix in pixel i
            s2i(i1,i2)=1.d0/si
            cofh = 2.d0*log(2.d0*varest(i1,i2)/si*fov)
            cofh=sqrt(cofh)
        if(max(a-thi,thi-b)/sqrt(varest(i1,i2))-cofh.gt.thresh) THEN
               fix(i1,i2)=.TRUE.
               if(segm(i1,i2).eq.0) segm(i1,i2)=sign(1.d0,thi-level)
C wee need to assign a value to segment before we can fix the decision
            ELSE
               fix(i1,i2)=.FALSE.
               ti=max(0.d0,max(a-thi,thi-b))
               pvalue(i1,i2)=fpchisq(ti*s2i(i1,i2),1.d0,1,0)
            END IF
         END DO
      END DO
      DO i2=1,n2
         DO i1=1,n1
            if(fix(i1,i2)) THEN
               thnew(i1,i2)=theta(i1,i2)
               CYCLE
            END IF
            bii=bi(i1,i2)/lambda
C   scaling of sij outside the loop
            swj=0.d0
            swjw=0.d0
            swjw2=0.d0
            swjy=0.d0
            thi=theta(i1,i2)
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
                  thij=thi-theta(j1,j2)
                  wj=lw(jw1+jwind2)
                  IF (aws) THEN
                  IF(segm(i1,i2)*segm(j1,j2).gt.0) THEN
                      sij=max(pvalue(i1,i2),pvalue(j1,j2))*
     1                        bii*kldistgc(thij,s2i(i1,i2),1)
                  ELSE
                     sij=bii*kldistgc(thij,s2i(i1,i2),1)
                  END IF
                     IF (sij.gt.1.d0) THEN
                        CYCLE
                     END IF
                     IF (sij.gt.spmin) THEN
                         wj=wj*(1.d0-spf*(sij-spmin))
                     END IF 
                  END IF
                  wj=wj
                  swj=swj+wj
                  swjy=swjy+wj*y(j1,j2)
                  swjw=swjw+wj*s2i(i1,i2)
                  swjw2=swjw2+wj*wj*s2i(i1,i2)
               END DO
            END DO
            thi=swjy/swj
            thnew(i1,i2)=thi
            bi(i1,i2)=swj
            si=swjw2/swjw/swjw
            varest(i1,i2)=si
            cofh = 2.d0*log(2.d0*si*s2i(i1,i2)*fov)
            si=sqrt(si)
            cofh=sqrt(cofh)
            IF((thi-a)/si+cofh.lt.-ext) THEN
               segm(i1,i2)=-1
            ELSE IF ((thi-b)/si-cofh.gt.ext) THEN
               segm(i1,i2)=1
            ELSE
               segm(i1,i2)=0
            END IF
            call rchkusr()
         END DO
      END DO
      RETURN
      END
      subroutine connect(segm,n1,n2,i1,i2,ind1,ind2,checked)
      implicit logical (a-z)
      integer n1,n2,segm(n1,n2),i1,i2,ind1(1),ind2(1)
      logical final,checked(1)
      integer j1,j2,k,l1,l2,lind,lind0
C     first find pixel close to (i1,i2) with segm(j1,j2)=0
      if(segm(i1,i2).ne.0) THEN
         final=.FALSE.
         DO k=1,n1
            DO l1=-k,k
               DO l2=-k,k
                  if(max(abs(l1),abs(l2)).ne.k) CYCLE
                  j1=i1+l1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  j2=i2+l2
                  if(j2.lt.1.or.j2.gt.n2) CYCLE
                  if(segm(j1,j2).eq.0) THEN
                     final=.TRUE.
                     i1=j1
                     i2=j2
                  END IF
                  if(final) EXIT
               END DO
               if(final) EXIT
            END DO
            if(final) EXIT
         END DO
      END IF
      segm(i1,i2)=2
      ind1(1)=i1
      ind2(1)=i2
      lind=1
      lind0=1
      DO k=1,n1*n2
         checked(k)=.FALSE.
      END DO
      final=.FALSE.
      DO while(.not.final)
         DO k=1,lind0
            if(checked(k)) CYCLE
            DO l1=-1,1
               DO l2=-1,1
                  if(l1.eq.0.and.l2.eq.0) CYCLE
                  j1=ind1(k)+l1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  j2=ind2(k)+l2
                  if(j2.lt.1.or.j2.gt.n2) CYCLE
                  if(segm(j1,j2).eq.0) THEN
                     segm(j1,j2)=2
                     lind=lind+1
                     ind1(lind)=j1
                     ind2(lind)=j2
                  END IF
               END DO
            END DO 
         END DO
         if(lind.eq.lind0) THEN
            final=.TRUE.
         ELSE
            lind0=lind
         END IF
      END DO
      RETURN
      END
      subroutine connect1(segm,n1,n2,i1,i2,ind1,ind2,checked)
      implicit logical (a-z)
      integer n1,n2,segm(n1,n2),i1,i2,ind1(1),ind2(1)
      logical final,checked(1)
      integer j1,j2,k,l1,l2,lind,lind0,isegm
C     first find pixel close to (i1,i2) with segm(j1,j2)=0
      isegm=segm(i1,i2)
      segm(i1,i2)=2
      ind1(1)=i1
      ind2(1)=i2
      lind=1
      lind0=1
      DO k=1,n1*n2
         checked(k)=.FALSE.
      END DO
      final=.FALSE.
      DO while(.not.final)
         DO k=1,lind0
            if(checked(k)) CYCLE
            DO l1=-1,1
               DO l2=-1,1
                  if(l1.eq.0.and.l2.eq.0) CYCLE
                  j1=ind1(k)+l1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  j2=ind2(k)+l2
                  if(j2.lt.1.or.j2.gt.n2) CYCLE
                  if(segm(j1,j2).eq.isegm) THEN
                     segm(j1,j2)=2
                     lind=lind+1
                     ind1(lind)=j1
                     ind2(lind)=j2
                  END IF
               END DO
            END DO 
         END DO
         if(lind.eq.lind0) THEN
            final=.TRUE.
         ELSE
            lind0=lind
         END IF
      END DO
      RETURN
      END
