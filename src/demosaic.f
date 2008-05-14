      subroutine halfsize(sensor,theta,n1,n2,h1,h2,bayer)
C
C   h1 = (n1%/%2)-1
C   h2 = (n2%/%2)-1
C
      implicit logical (a-z)
      external channel
      integer n1,n2,h1,h2,sensor(n1,n2),theta(h1,h2,3),bayer
      integer j1,j2,i1,i2,channel,s(3,3),k1,k2,ch
      DO j1=1,h1
         DO j2=1,h2
            i1=2*j1
            i2=2*j2
            DO k1=1,3
               DO k2=1,3
                  s(k1,k2)=0
               END DO
            END DO
            ch = channel(i1,i2,bayer)
            s(ch,1)=s(ch,1)+sensor(i1,i2)
            ch = channel(i1+1,i2,bayer)
            s(ch,1)=s(ch,1)+sensor(i1+1,i2)
            ch = channel(i1+1,i2+1,bayer)
            s(ch,1)=s(ch,1)+sensor(i1+1,i2+1)
            ch = channel(i1,i2+1,bayer)
            s(ch,1)=s(ch,1)+sensor(i1,i2+1)
            ch = channel(i1-1,i2,bayer)
            s(ch,2)=s(ch,2)+sensor(i1-1,i2)
            ch = channel(i1-1,i2+1,bayer)
            s(ch,2)=s(ch,2)+sensor(i1-1,i2+1)
            ch = channel(i1+2,i2,bayer)
            s(ch,2)=s(ch,2)+sensor(i1+2,i2)
            ch = channel(i1+2,i2+1,bayer)
            s(ch,2)=s(ch,2)+sensor(i1+2,i2+1)
            ch = channel(i1,i2-1,bayer)
            s(ch,2)=s(ch,2)+sensor(i1,i2-1)
            ch = channel(i1+1,i2-1,bayer)
            s(ch,2)=s(ch,2)+sensor(i1+1,i2-1)
            ch = channel(i1,i2+2,bayer)
            s(ch,2)=s(ch,2)+sensor(i1,i2+2)
            ch = channel(i1+1,i2+2,bayer)
            s(ch,2)=s(ch,2)+sensor(i1+1,i2+2)
            ch = channel(i1-1,i2-1,bayer)
            s(ch,3)=s(ch,3)+sensor(i1-1,i2-1)
            ch = channel(i1-1,i2+2,bayer)
            s(ch,3)=s(ch,3)+sensor(i1-1,i2+2)
            ch = channel(i1+2,i2+2,bayer)
            s(ch,3)=s(ch,3)+sensor(i1+2,i2+2)
            ch = channel(i1+2,i2-1,bayer)
            s(ch,3)=s(ch,3)+sensor(i1+2,i2-1)
            theta(j1,j2,1)=(9*s(1,1)+3*s(1,2)+s(1,3))/16
            theta(j1,j2,2)=(18*s(2,1)+9*s(2,2)+4*s(2,3))/80
            theta(j1,j2,3)=(9*s(3,1)+3*s(3,2)+s(3,3))/16
         END DO
      END DO
      return
      end
      subroutine fullsize(sensor,theta,n1,n2,h1,h2,bayer)
      implicit logical (a-z)
      external channel
      integer h1,h2,n1,n2,sensor(n1,n2),theta(h1,h2,3),bayer
      integer i1,i2,channel,s(3,3),k1,k2,ch,j1,j2
      DO i1=3,n1-2
         DO i2=3,n2-2
            DO k1=1,3
               DO k2=1,3
                  s(k1,k2)=0
               END DO
            END DO
            j1=i1-2
            j2=i2-2
            ch = channel(i1,i2,bayer)
            s(ch,1)=s(ch,1)+sensor(i1,i2)
            ch = channel(i1+1,i2,bayer)
            s(ch,1)=s(ch,1)+sensor(i1+1,i2)
            ch = channel(i1+1,i2+1,bayer)
            s(ch,1)=s(ch,1)+sensor(i1+1,i2+1)
            ch = channel(i1,i2+1,bayer)
            s(ch,1)=s(ch,1)+sensor(i1,i2+1)
            ch = channel(i1-1,i2,bayer)
            s(ch,2)=s(ch,2)+sensor(i1-1,i2)
            ch = channel(i1-1,i2+1,bayer)
            s(ch,2)=s(ch,2)+sensor(i1-1,i2+1)
            ch = channel(i1+2,i2,bayer)
            s(ch,2)=s(ch,2)+sensor(i1+2,i2)
            ch = channel(i1+2,i2+1,bayer)
            s(ch,2)=s(ch,2)+sensor(i1+2,i2+1)
            ch = channel(i1,i2-1,bayer)
            s(ch,2)=s(ch,2)+sensor(i1,i2-1)
            ch = channel(i1+1,i2-1,bayer)
            s(ch,2)=s(ch,2)+sensor(i1+1,i2-1)
            ch = channel(i1,i2+2,bayer)
            s(ch,2)=s(ch,2)+sensor(i1,i2+2)
            ch = channel(i1+1,i2+2,bayer)
            s(ch,2)=s(ch,2)+sensor(i1+1,i2+2)
            ch = channel(i1-1,i2-1,bayer)
            s(ch,3)=s(ch,3)+sensor(i1-1,i2-1)
            ch = channel(i1-1,i2+2,bayer)
            s(ch,3)=s(ch,3)+sensor(i1-1,i2+2)
            ch = channel(i1+2,i2+2,bayer)
            s(ch,3)=s(ch,3)+sensor(i1+2,i2+2)
            ch = channel(i1+2,i2-1,bayer)
            s(ch,3)=s(ch,3)+sensor(i1+2,i2-1)
            theta(j1,j2,1)=(9*s(1,1)+3*s(1,2)+s(1,3))/16
            theta(j1,j2,2)=(18*s(2,1)+9*s(2,2)+4*s(2,3))/80
            theta(j1,j2,3)=(9*s(3,1)+3*s(3,2)+s(3,3))/16
         END DO
      END DO
      return
      end
      subroutine indemos4(sensor,theta,n1,n2,bayer,bi,bi3)
C
C   this is bilinear interpolation
C
      implicit logical (a-z)
      integer n1,n2,sensor(n1,n2),theta(n1,n2,3),bi(n1,n2),
     1        bi3(n1,n2,3),bayer
      integer i1,i2,icolor,channel,sn(8),bni(8),which
      external channel
      DO i1=1,n1
         DO i2=1,n2
            icolor=channel(i1,i2,bayer)
            call neighbor(sensor,bi,n1,n2,i1,i2,bayer,sn,bni,which)
            if(icolor.eq.1) THEN
               call inred4(sn,sensor(i1,i2),bni,bi(i1,i2),bi3(i1,i2,1),
     1                    bi3(i1,i2,2),bi3(i1,i2,3),theta(i1,i2,1),
     2                    theta(i1,i2,2),theta(i1,i2,3))
            ELSE IF(icolor.eq.2) THEN
            call ingreen4(sn,sensor(i1,i2),bni,bi(i1,i2),bi3(i1,i2,1),
     1                    bi3(i1,i2,2),bi3(i1,i2,3),theta(i1,i2,1),
     2                    theta(i1,i2,2),theta(i1,i2,3),which)
            ELSE
               call inblue4(sn,sensor(i1,i2),bni,bi(i1,i2),bi3(i1,i2,1),
     1                    bi3(i1,i2,2),bi3(i1,i2,3),theta(i1,i2,1),
     2                    theta(i1,i2,2),theta(i1,i2,3))
            ENDIF
         END DO
      END DO
      RETURN
      END
      subroutine neighbor(sensor,bisen,n1,n2,i1,i2,bayer,sn,bi,which)
C
C    copy sensor data from neighboring pixel clockwise into sn
C
      implicit logical (a-z)
      external channel
      logical i1a,i1e,i2a,i2e
      integer n1,n2,sensor(n1,n2),sn(8),i1,i2,j,which,bayer
      integer bisen(n1,n2),bi(8),channel
      i1a=i1.gt.1
      i1e=i1.lt.n1
      i2a=i2.gt.1
      i2e=i2.lt.n2
      DO j=1,8
         sn(j)=-65536
         bi(j)=-65536
      END DO
      which=channel(i1,i2+1,bayer)
      IF(i1a.and.i2e) THEN
          sn(1)=sensor(i1-1,i2+1)
          bi(1)=bisen(i1-1,i2+1)
      END IF
      IF(i2e) THEN
          sn(2)=sensor(i1,i2+1)
          bi(2)=bisen(i1,i2+1)
      END IF
      IF(i1e.and.i2e) THEN
          sn(3)=sensor(i1+1,i2+1)
          bi(3)=bisen(i1+1,i2+1)
      END IF
      IF(i1e) THEN
          sn(4)=sensor(i1+1,i2)
          bi(4)=bisen(i1+1,i2)
      END IF
      IF(i1e.and.i2a) THEN
          sn(5)=sensor(i1+1,i2-1)
          bi(5)=bisen(i1+1,i2-1)
      END IF
      IF(i2a) THEN
          sn(6)=sensor(i1,i2-1)
          bi(6)=bisen(i1,i2-1)
      END IF
      IF(i1a.and.i2a) THEN
          sn(7)=sensor(i1-1,i2-1)
          bi(7)=bisen(i1-1,i2-1)
      END IF
      IF(i1a) THEN
          sn(8)=sensor(i1-1,i2)
          bi(8)=bisen(i1-1,i2)
      END IF
      DO j=1,8
C   make edges of the image homogeneous by mirroring
         IF(sn(j).lt.0) THEN
              sn(j)=sn(mod(j-1+4,8)+1)
              bi(j)=bi(mod(j-1+4,8)+1)
         END IF
         IF(sn(1).lt.0) sn(1)=max(sn(3),sn(5),sn(7))
         IF(sn(3).lt.0) sn(3)=max(sn(1),sn(5),sn(7))
         IF(sn(5).lt.0) sn(5)=max(sn(1),sn(3),sn(7))
         IF(sn(7).lt.0) sn(7)=max(sn(1),sn(3),sn(5))
         IF(bi(1).lt.0) bi(1)=max(bi(3),bi(5),bi(7))
         IF(bi(3).lt.0) bi(3)=max(bi(1),bi(5),bi(7))
         IF(bi(5).lt.0) bi(5)=max(bi(1),bi(3),bi(7))
         IF(bi(7).lt.0) bi(7)=max(bi(1),bi(3),bi(5))
      END DO
      return
      end
      subroutine ingreen4(sn,sni,bi,bii,bir,big,bib,red,green,blue,
     1                   which)
C
C   demosaicing for green pixel 
C   sensori contains the observed green pixel
C   sn the sendor data from neighboring pixel (clockwise)
C   which contains information on wether 
C   sn(2) corresponds to a red (which=1) or blue (which=3) pixel
C
      implicit logical (a-z)
      integer sn(8),sni,red,green,blue,which,bi(8),bii,bib,bir,big
C   first check if we have homogeneity based on green
      green=sni
      big=bii
      if(which.eq.1) THEN
         red=0.5d0*(sn(2)+sn(6))
         blue=0.5d0*(sn(4)+sn(8))
         bir=0.5d0*(bi(2)+bi(6))
         bib=0.5d0*(bi(4)+bi(8))
      ELSE
         blue=0.5d0*(sn(2)+sn(6))
         red=0.5d0*(sn(4)+sn(8))
         bib=0.5d0*(bi(2)+bi(6))
         bir=0.5d0*(bi(4)+bi(8))
      END IF
      return
      end
      subroutine inred4(sn,sni,bi,bii,bir,big,bib,red,green,blue)
C
C   demosaicing for red pixel 
C   sni contains the observed red pixel
C   sn the sendor data from neighboring pixel (clockwise)
C   
      implicit logical (a-z)
      integer sn(8),sni,red,green,blue,bi(8),bii,bir,big,bib
      red=sni
      bir=bii
      blue=(sn(1)+sn(3)+sn(5)+sn(7))*.25d0
      green=(sn(2)+sn(4)+sn(6)+sn(8))*.25d0
      bib=(bi(1)+bi(3)+bi(5)+bi(7))*.25d0
      big=(bi(2)+bi(4)+bi(6)+bi(8))*.25d0
      return
      end
      subroutine inblue4(sn,sni,bi,bii,bir,big,bib,red,green,blue)
C
C   demosaicing for blue pixel 
C   sni contains the observed blue pixel
C   sn the sendor data from neighboring pixel (clockwise)
C   
      implicit logical (a-z)
      integer sn(8),sni,red,green,blue,bi(8),bii,
     1       bir,big,bib
      blue=sni
      bib=bii
      red=(sn(1)+sn(3)+sn(5)+sn(7))*.25d0
      green=(sn(2)+sn(4)+sn(6)+sn(8))*.25d0
      bir=(bi(1)+bi(3)+bi(5)+bi(7))*.25d0
      big=(bi(2)+bi(4)+bi(6)+bi(8))*.25d0
      return
      end

      subroutine wbalance(sensor,n1,n2,wb,bayer)
      implicit logical (a-z)
      external channel
      integer n1,n2,sensor(n1,n2),bayer,channel,z
      real*8 wb(3)
      integer i1,i2
      DO i1=1,n1
         DO i2=1,n2
            z=sensor(i1,i2)*wb(channel(i1,i2,bayer))
            if(z.gt.65535) z = 65535
            sensor(i1,i2)=z
         END DO
      END DO
      return
      end
C##########################################################################
C
C        identify is a pixel is r, g or b in a Bayer map
C
C##########################################################################
      integer function channel(i,j,bayer)
      implicit logical (a-z)
      integer i,j,k,l,bayer
         k=mod(i,2)
         l=mod(j,2)
         channel=1
         IF(bayer.eq.1) THEN
C  e.g. Canon Powershot S30 (Bayer RGGB)
            IF(k+l.ne.1) THEN 
               channel=2
            ELSE IF(k.eq.1) THEN
               channel=1
            ELSE
               channel=3
            ENDIF
C  (Bayer GRBG)
         ELSE IF(bayer.eq.2) THEN
            IF(k+l.eq.1) THEN 
               channel=2
            ELSE IF(k.ne.0) THEN
               channel=3
            ELSE
               channel=1
            ENDIF
         ELSE IF(bayer.eq.3) THEN
C  (Bayer BGGR)
            IF(k+l.ne.1) THEN 
               channel=2
            ELSE IF(k.eq.1) THEN
               channel=3
            ELSE
               channel=1
            ENDIF
         ELSE IF(bayer.eq.4) THEN
C  e.g. Lumix LX 2  (Bayer GBRG)
            IF(k+l.eq.1) THEN 
               channel=2
            ELSE IF(k.ne.0) THEN
               channel=1
            ELSE
               channel=3
            ENDIF
         END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  array(as.integer(pmax(0,pmin(zobj$theta %*% out.cam,65535))),c(dimg,3))
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cam2rgb(theta,n,outcam,thetanew)
      implicit logical (a-z)
      integer n,theta(n,3),thetanew(n,3) 
      real*8 outcam(3,3)
      integer i,j,k,iz
      real*8 z
      DO i=1,n
         DO j=1,3
            z=0.d0
            DO k=1,3
               z=z+theta(i,k)*outcam(k,j)
            END DO
            iz=z
            iz=max(min(iz,65535),0)
            thetanew(i,j)=iz
         END DO
      END DO
      RETURN
      END
