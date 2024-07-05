      subroutine smsensor(s,shat,th,n1,n2,nt1,nt2,bayer,vcoef,
     1                    meanvar,hakt,lambda,bi,kern,spmin,lw)
      implicit none
      integer n1,n2,nt1,nt2,s(n1,n2),shat(n1,n2),th(nt1,nt2,3),
     1        kern,bayer
      double precision vcoef(2,3),bi(n1,n2),hakt,spmin,lw(*),lambda,
     1       meanvar(3)
      logical aws
      external channel,lkern,kldisdem
      double precision lkern,kldisdem
      integer channel
      integer i1,i2,i1th,i2th,j1,j2,ch,dlw,clw,ih,ih1,ja1,je1,
     1        jind2,jind,k,j1th,j2th,jw2,jw1,jwind2
      double precision spf,hakt2,z1,z2,wj,sw,swy,thi(3),thij(3),s2i(3),
     1       bii,sij
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih=int(hakt)
      dlw=2*ih+1
      clw=ih+1
      aws=lambda.lt.1d40
      DO j2=1,dlw
         z2=clw-j2
         z2=z2*z2
         ih1=int(sqrt(hakt2-z2))
         ja1=max(1,clw-ih1)
         je1=min(dlw,clw+ih1)
         jind2=(j2-1)*dlw
         DO j1=ja1,je1
C  first location weights
            jind=j1+jind2
            z1=clw-j1
            wj=lkern(kern,(z1*z1+z2)/hakt2)
            lw(jind)=wj
         END DO
      END DO
      DO i1=1,n1
         i1th=min(n1-1,max(2,i1))-1
         DO i2=1,n2
            i2th=min(n2-1,max(2,i2))-1
            ch=channel(i1,i2,bayer)
            bii=bi(i1,i2)/lambda
            sw=0.d0
            swy=0.d0
            DO k=1,3
               thi(k)=th(i1th,i2th,k)
C  k in vcoef refers to variance parameters in the channels of the sensor data
               s2i(k) = 1.d0/max(vcoef(1,k) + vcoef(2,k) * thi(k),
     1                               .1*meanvar(k))
C  thats inverse standard deviation
            END DO
            DO jw2=1,dlw
               j2=jw2-clw+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               j2th=min(n2-1,max(2,j2))-1
               jwind2=(jw2-1)*dlw
               z2=clw-jw2
               z2=z2*z2
               ih1=int(sqrt(hakt2-z2))
               DO jw1=clw-ih1,clw+ih1
                  j1=jw1-clw+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  if(channel(j1,j2,bayer).ne.ch) CYCLE
C  we don't smooth sensor data from different channels
                  j1th=min(n1-1,max(2,j1))-1
                  DO k=1,3
                     thij(k)=thi(k)-th(j1th,j2th,k)
                  END DO
                  wj=lw(jw1+jwind2)
                  IF (aws) THEN
                     sij=bii*kldisdem(thij,s2i)
                     IF (sij.gt.1.d0) CYCLE
                     IF (sij.gt.spmin) wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
                  sw=sw+wj
                  swy=swy+wj*s(j1,j2)
               END DO
            END DO
            shat(i1,i2)=int(swy/sw)
            bi(i1,i2)=sw
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   search for homogeneous regions
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dhomogen(img,n1,n2,imghom,bayer)
      implicit none
      external channel
      integer i,j,n1,n2,bayer,ch,channel
      integer img(n1,n2),imghom(n1,n2)
      double precision thi,thi2,r1,r2,r3,r4,r5,r6,r7,r8,r9
      DO i=3,n1-2
         DO j=3,n2-2
            ch=channel(i,j,bayer)
            if(ch.eq.2) THEN
               r1=img(i,j)
               r2=img(i-1,j-1)
               r3=img(i+1,j-1)
               r4=img(i-1,j+1)
               r5=img(i+1,j+1)
               r6=img(i,j-2)
               r7=img(i,j+2)
               r8=img(i-2,j)
               r9=img(i+2,j)
               thi=(r1+r2+r3+r4+r5+r6+r7+r8+r9)/9.d0
               if(thi.gt.0) THEN
                  thi2=(r1*r1+r2*r2+r3*r3+r4*r4+r5*r5+
     1                  r6*r6+r7*r7+r8*r8+r9*r9)/9.d0
                  imghom(i,j)=int(thi2/thi-thi)
               ELSE
                  imghom(i,j)=0
               END IF
            ELSE
               r1=img(i,j)
               r2=img(i-2,j)
               r3=img(i+2,j)
               r4=img(i,j-2)
               r5=img(i,j+2)
               thi=(r1+r2+r3+r4+r5)/5.d0
               if(thi.gt.0) THEN
                  thi2=(r1*r1+r2*r2+r3*r3+r4*r4+r5*r5)/5.d0
                  imghom(i,j)=int(thi2/thi-thi)
               ELSE
                  imghom(i,j)=0
               END IF
            END IF
         END DO
      END DO
      RETURN
      END
      subroutine smsens0(s,shat,bi,n1,n2,bayer)
      implicit none
      integer n1,n2,s(n1,n2),shat(n1,n2),bayer
      double precision bi(n1,n2)
      external channel
      integer channel
      integer i1,i2,j1,j2,ch,ih1,jw2,jw1
      double precision hakt,hakt2,z2,sw,swy
      hakt=2.1d0
      hakt2=hakt*hakt
      DO i1=1,n1
         DO i2=1,n2
            ch=channel(i1,i2,bayer)
            sw=0.d0
            swy=0.d0
            DO jw2=-2,2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               z2=jw2*jw2
               ih1=int(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  if(channel(j1,j2,bayer).ne.ch) CYCLE
C  we don't smooth sensor data from different channel
                  sw=sw+1.d0
                  swy=swy+s(j1,j2)
               END DO
            END DO
            shat(i1,i2)=int(swy/sw)
            bi(i1,i2)=sw
         END DO
      END DO
      RETURN
      END
      subroutine senvar(s,n1,n2,shat,bi,bayer,vcoef,mvar,nothom)
      implicit none
      integer n1,n2,s(n1,n2),shat(n1,n2),bayer,nothom(n1,n2)
      double precision vcoef(2,3),bi(n1,n2),mvar(3)
      external channel
      integer channel,ch,i1,i2,n(3)
      double precision s0(3),s1(3),s2(3),t0(3),t1(3),ms(3),bii,wght,
     1       si,d,z,res
      DO ch=1,3
         s0(ch)=0.d0
         s1(ch)=0.d0
         s2(ch)=0.d0
         t0(ch)=0.d0
         t1(ch)=0.d0
         ms(ch)=0.d0
         n(ch)=0
      END DO
      DO i1=1,n1
         DO i2=1,n2
            if(nothom(i1,i2).ne.0) CYCLE
            ch=channel(i1,i2,bayer)
            ms(ch)=ms(ch)+s(i1,i2)
            n(ch)=n(ch)+1
            bii=bi(i1,i2)
            if(bii.le.1.d0) CYCLE
            si=shat(i1,i2)
            wght=bii-1.d0
            res=s(i1,i2)-si
            res=res*res*bii/wght
            s0(ch)=s0(ch)+wght
            z=wght*si
            s1(ch)=s1(ch)+z
            s2(ch)=s2(ch)+z*si
            t0(ch)=t0(ch)+wght*res
            t1(ch)=t1(ch)+z*res
         END DO
      END DO
      DO ch=1,3
         d=s2(ch)*s0(ch)-s1(ch)*s1(ch)
         IF(d.gt.0.d0) THEN
            vcoef(1,ch)=(s2(ch)*t0(ch)-s1(ch)*t1(ch))/d
            vcoef(2,ch)=(-s1(ch)*t0(ch)+s0(ch)*t1(ch))/d
         ELSE
            vcoef(1,ch)=1d-2
            vcoef(2,ch)=0d0
         END IF
         mvar(ch)=vcoef(1,ch)+vcoef(2,ch)*ms(ch)/n(ch)
      END DO
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
      double precision function kldisdem(thij,s2ii)
      implicit none
      double precision thij(3),s2ii(3)
      kldisdem = thij(1)*thij(1)*s2ii(1)+
     1           thij(2)*thij(2)*s2ii(2)+
     2           thij(3)*thij(3)*s2ii(3)
      RETURN
      END
