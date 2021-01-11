      subroutine halfsize(sensor,theta,n1,n2,h1,h2,bayer)
C
C   h1 = (n1%/%2)-1
C   h2 = (n2%/%2)-1
C
      implicit none
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
      subroutine median2(x,n,y,tol)
      implicit none
      integer n,x(2,n),y(2)
      double precision tol
      integer i,j
      double precision etaofy,di1,di2,rofy,dxy,z,delta,normy,
     1       y1,y2,r1,r2,t1,t2,t0,c1,c2
C  use mean as init
      y1=x(1,1)
      y2=x(2,1)
      DO i=2,n
         y1=y1+x(1,i)
         y2=y2+x(2,i)
      END DO
      y1=y1/n
      y2=y2/n
C  iterate until convergence
      rofy = 1.d10
      j=0
      DO while (rofy.gt.tol)
C  compute r(y) and check for y=x_k
         etaofy=0.d0
         r1=0.d0
         r2=0.d0
         t0=0.d0
         t1=0.d0
         t2=0.d0
         DO i=1,n
            di1=x(1,i)-y1
            di2=x(2,i)-y2
            dxy=sqrt(di1*di1+di2*di2)
            if(dxy.lt.1e-8) THEN
               etaofy=etaofy+1.d0
            ELSE
               r1=r1+di1/dxy
               r2=r2+di2/dxy
               t0=t0+1.d0/dxy
               t1=t1+x(1,i)/dxy
               t2=t2+x(2,i)/dxy
            END IF
         END DO
         rofy=sqrt(r1*r1+r2*r2)
         if(rofy.le.tol) EXIT
         t1=t1/t0
         t2=t2/t0
         etaofy=etaofy/rofy
         c1=max(0.d0,1.d0-etaofy)
         c2=min(1.d0,etaofy)
         z=c1*t1+c2*y1
         delta=abs(y1-z)
         normy=1.d0+abs(z)
         y1=z
         z=c1*t2+c2*y2
         delta=delta+abs(y2-z)
         normy=normy+abs(z)
         y2=z
         if(delta.lt.tol*normy) EXIT
         j=j+1
         if(j.gt.20) EXIT
      END DO
      y(1)=int(y1)
      y(2)=int(y2)
      RETURN
      END
      subroutine demmed4(sensor,theta,n1,n2,h1,h2,bayer)
      implicit none
      external channel
      integer n1,n2,h1,h2,sensor(n1,n2),theta(h1,h2,3),bayer
      integer i1,i2,channel,ch,ch1,z(2,8),y(2)
C  h1 == n1-2,  h2 == n2-2
      DO i1=2,n1-1
         DO i2=2,n2-1
            ch = channel(i1,i2,bayer)
            theta(i1-1,i2-1,ch)=sensor(i1,i2)
            if(ch.eq.2) THEN
C   green channel
               ch1 = channel(i1-1,i2,bayer)
               if(ch1.eq.1) THEN
C   red channel in i1-1,i2
                  z(1,1)=sensor(i1-1,i2)
                  z(2,1)=sensor(i1,i2-1)
                  z(1,2)=sensor(i1+1,i2)
                  z(2,2)=sensor(i1,i2-1)
                  z(1,3)=sensor(i1-1,i2)
                  z(2,3)=sensor(i1,i2+1)
                  z(1,4)=sensor(i1+1,i2)
                  z(2,4)=sensor(i1,i2+1)
               ELSE
                  z(1,1)=sensor(i1,i2-1)
                  z(2,1)=sensor(i1-1,i2)
                  z(1,2)=sensor(i1,i2-1)
                  z(2,2)=sensor(i1+1,i2)
                  z(1,3)=sensor(i1,i2+1)
                  z(2,3)=sensor(i1-1,i2)
                  z(1,4)=sensor(i1,i2+1)
                  z(2,4)=sensor(i1+1,i2)
                END IF
               call median2(z,4,y,1.d-3)
               theta(i1-1,i2-1,1)=y(1)
               theta(i1-1,i2-1,3)=y(2)
            ELSE
C  red or blue channel
               z(1,1)=sensor(i1-1,i2)
               z(2,1)=sensor(i1-1,i2-1)
               z(1,2)=sensor(i1-1,i2)
               z(2,2)=sensor(i1-1,i2+1)
               z(1,3)=sensor(i1+1,i2)
               z(2,3)=sensor(i1+1,i2-1)
               z(1,4)=sensor(i1+1,i2)
               z(2,4)=sensor(i1+1,i2+1)
               z(1,5)=sensor(i1,i2-1)
               z(2,5)=sensor(i1-1,i2-1)
               z(1,6)=sensor(i1,i2-1)
               z(2,6)=sensor(i1+1,i2-1)
               z(1,7)=sensor(i1,i2+1)
               z(2,7)=sensor(i1-1,i2+1)
               z(1,8)=sensor(i1,i2+1)
               z(2,8)=sensor(i1+1,i2+1)
               call median2(z,8,y,1.d-3)
               theta(i1-1,i2-1,2)=y(1)
               theta(i1-1,i2-1,4-ch)=y(2)
            END IF
         END DO
      END DO
      RETURN
      END
      subroutine median2b(x,n,y,tol)
      implicit none
      integer n,x(2,n),y(2)
      double precision tol
      integer i,j
      double precision etaofy,di1,di2,rofy,dxy,z,delta,normy,
     1       y1,y2,r1,r2,t1,t2,t0,c1,c2
C  use mean as init
      y1=y(1)
      y2=y(2)
C  iterate until convergence
      rofy = 1.d10
      j=0
      DO while (rofy.gt.tol)
C  compute r(y) and check for y=x_k
         etaofy=0.d0
         r1=0.d0
         r2=0.d0
         t0=0.d0
         t1=0.d0
         t2=0.d0
         DO i=1,n
            di1=x(1,i)-y1
            di2=x(2,i)-y2
            dxy=sqrt(di1*di1+di2*di2)
            if(dxy.lt.1e-8) THEN
               etaofy=etaofy+1.d0
            ELSE
               r1=r1+di1/dxy
               r2=r2+di2/dxy
               t0=t0+1.d0/dxy
               t1=t1+x(1,i)/dxy
               t2=t2+x(2,i)/dxy
            END IF
         END DO
         rofy=sqrt(r1*r1+r2*r2)
         if(rofy.le.tol) EXIT
         t1=t1/t0
         t2=t2/t0
         etaofy=etaofy/rofy
         c1=max(0.d0,1.d0-etaofy)
         c2=min(1.d0,etaofy)
         z=c1*t1+c2*y1
         delta=abs(y1-z)
         normy=1.d0+abs(z)
         y1=z
         z=c1*t2+c2*y2
         delta=delta+abs(y2-z)
         normy=normy+abs(z)
         y2=z
         if(delta.lt.tol*normy) EXIT
         j=j+1
         if(j.gt.5) EXIT
      END DO
      y(1)=int(y1)
      y(2)=int(y2)
      RETURN
      END
      subroutine demmed4b(sensor,theta,n1,n2,h1,h2,bayer)
      implicit none
      external channel
      integer n1,n2,h1,h2,sensor(n1,n2),theta(h1,h2,3),bayer
      integer i1,i2,channel,ch,ch1,z(2,8),y(2)
C  h1 == n1-2,  h2 == n2-2
      DO i1=2,n1-1
         DO i2=2,n2-1
            ch = channel(i1,i2,bayer)
            theta(i1-1,i2-1,ch)=sensor(i1,i2)
            if(ch.eq.2) THEN
C   green channel
               ch1 = channel(i1-1,i2,bayer)
               y(1) = theta(i1-1,i2-1,1)
               y(2) = theta(i1-1,i2-1,3)
               if(ch1.eq.1) THEN
C   red channel in i1-1,i2
                  z(1,1)=sensor(i1-1,i2)
                  z(2,1)=sensor(i1,i2-1)
                  z(1,2)=sensor(i1+1,i2)
                  z(2,2)=sensor(i1,i2-1)
                  z(1,3)=sensor(i1-1,i2)
                  z(2,3)=sensor(i1,i2+1)
                  z(1,4)=sensor(i1+1,i2)
                  z(2,4)=sensor(i1,i2+1)
               ELSE
                  z(1,1)=sensor(i1,i2-1)
                  z(2,1)=sensor(i1-1,i2)
                  z(1,2)=sensor(i1,i2-1)
                  z(2,2)=sensor(i1+1,i2)
                  z(1,3)=sensor(i1,i2+1)
                  z(2,3)=sensor(i1-1,i2)
                  z(1,4)=sensor(i1,i2+1)
                  z(2,4)=sensor(i1+1,i2)
                END IF
               call median2b(z,4,y,1.d-3)
               theta(i1-1,i2-1,1)=y(1)
               theta(i1-1,i2-1,3)=y(2)
            ELSE
C  red or blue channel
               z(1,1)=sensor(i1-1,i2)
               z(2,1)=sensor(i1-1,i2-1)
               z(1,2)=sensor(i1-1,i2)
               z(2,2)=sensor(i1-1,i2+1)
               z(1,3)=sensor(i1+1,i2)
               z(2,3)=sensor(i1+1,i2-1)
               z(1,4)=sensor(i1+1,i2)
               z(2,4)=sensor(i1+1,i2+1)
               z(1,5)=sensor(i1,i2-1)
               z(2,5)=sensor(i1-1,i2-1)
               z(1,6)=sensor(i1,i2-1)
               z(2,6)=sensor(i1+1,i2-1)
               z(1,7)=sensor(i1,i2+1)
               z(2,7)=sensor(i1-1,i2+1)
               z(1,8)=sensor(i1,i2+1)
               z(2,8)=sensor(i1+1,i2+1)
               y(1) = theta(i1-1,i2-1,2)
               y(2) = theta(i1-1,i2-1,4-ch)
               call median2b(z,8,y,1.d-3)
               theta(i1-1,i2-1,2)=y(1)
               theta(i1-1,i2-1,4-ch)=y(2)
            END IF
         END DO
      END DO
      RETURN
      END
      subroutine median16(x,n,th,tol)
      implicit none
      integer n,x(16,n),th(3)
      double precision y(16),tol
      integer i,j,k
      double precision etaofy,di(16),rofy,dxy,
     1       r(16),t(16),t0,c1,c2,z,delta,normy
C  use mean as init
      DO j=1,16
         y(j)=x(j,1)
      END DO
      DO i=1,n
         DO j=2,16
            y(j)=y(j)+x(j,i)
         END DO
      END DO
      DO j=1,16
         y(j)=y(j)/n
      END DO
C  iterate until convergence
      rofy = 1.d10
      k=0
      DO while (rofy.gt.tol)
C  compute r(y) and check for y=x_k
         etaofy=0.d0
         DO j=1,16
            r(j)=0.d0
            t(j)=0.d0
         END DO
         t0=0.d0
         DO i=1,n
            dxy=0.d0
            DO j=1,16
               z=x(j,i)-y(j)
               dxy=dxy+z*z
               di(j)=z
            END DO
            dxy=sqrt(dxy)
            if(dxy.lt.1e-8) THEN
               etaofy=etaofy+1.d0
            ELSE
               DO j=1,16
                  r(j)=r(j)+di(j)/dxy
                  t(j)=t(j)+x(j,i)/dxy
               END DO
               t0=t0+1.d0/dxy
            END IF
         END DO
         rofy=r(1)*r(1)
         DO j=2,16
            rofy=rofy+r(j)*r(j)
         END DO
         rofy=sqrt(rofy)
         if(rofy.le.tol) EXIT
         DO j=1,16
            t(j)=t(j)/t0
         END DO
         etaofy=etaofy/rofy
         c1=max(0.d0,1.d0-etaofy)
         c2=min(1.d0,etaofy)
         delta=0.d0
         normy=1.d0
         DO j=1,16
            z=c1*t(j)+c2*y(j)
            delta=delta+abs(y(j)-z)
            normy=normy+abs(z)
            y(j)=z
         END DO
         if(delta.lt.tol*normy) EXIT
         k=k+1
         if(k.gt.20) EXIT
      END DO
      th(2)=int((y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8))/8)
      th(1)=int((y(9)+y(10)+y(11)+y(12))/4)
      th(3)=int((y(13)+y(14)+y(15)+y(16))/4)
      RETURN
      END
      subroutine demmed16(sensor,theta,n1,n2,h1,h2,bayer)
      implicit none
      external channel
      integer n1,n2,h1,h2,sensor(n1,n2),theta(h1,h2,3),bayer
      integer i1,i2,channel,k1,k2,ch,j1,j2,z(16,16)
      integer ig,ir,ib,k,kk1,kk2,jj1,jj2,y(3)
C  h1 == n1-6,  h2 == n2-6
      DO i1=1,h1
         DO i2=1,h2
            k=1
            DO j1=0,3
               jj1=i1+j1
               DO j2=0,3
                  jj2=i2+j2
                  ig=1
                  ir=9
                  ib=13
                  DO k1=0,3
                     kk1=jj1+k1
                     DO k2=0,3
                        kk2=jj2+k2
                        ch = channel(kk1,kk2,bayer)
                        select case (ch)
                           case (1)
                           z(ir,k)=sensor(kk1,kk2)
                           ir=ir+1
                           case (2)
                           z(ig,k)=sensor(kk1,kk2)
                           ig=ig+1
                           case (3)
                           z(ib,k)=sensor(kk1,kk2)
                           ib=ib+1
                           CASE DEFAULT
C                           call intpr("wrong channel",13,channel,1)
                        END SELECT
                     END DO
                  END DO
                  k=k+1
               END DO
            END DO
            call median16(z,16,y,1.d-3)
            theta(i1,i2,1)=y(1)
            theta(i1,i2,2)=y(2)
            theta(i1,i2,3)=y(3)
         END DO
      END DO
      RETURN
      END
      subroutine fullsize(sensor,theta,n1,n2,h1,h2,bayer)
      implicit none
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
      implicit none
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
      implicit none
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
      implicit none
      integer sn(8),sni,red,green,blue,which,bi(8),bii,bib,bir,big
C   first check if we have homogeneity based on green
      green=sni
      big=bii
      if(which.eq.1) THEN
         red=int(0.5d0*(sn(2)+sn(6)))
         blue=int(0.5d0*(sn(4)+sn(8)))
         bir=int(0.5d0*(bi(2)+bi(6)))
         bib=int(0.5d0*(bi(4)+bi(8)))
      ELSE
         blue=int(0.5d0*(sn(2)+sn(6)))
         red=int(0.5d0*(sn(4)+sn(8)))
         bib=int(0.5d0*(bi(2)+bi(6)))
         bir=int(0.5d0*(bi(4)+bi(8)))
      END IF
      return
      end
      subroutine inred4(sn,sni,bi,bii,bir,big,bib,red,green,blue)
C
C   demosaicing for red pixel
C   sni contains the observed red pixel
C   sn the sendor data from neighboring pixel (clockwise)
C
      implicit none
      integer sn(8),sni,red,green,blue,bi(8),bii,bir,big,bib
      red=sni
      bir=bii
      blue=int((sn(1)+sn(3)+sn(5)+sn(7))*.25d0)
      green=int((sn(2)+sn(4)+sn(6)+sn(8))*.25d0)
      bib=int((bi(1)+bi(3)+bi(5)+bi(7))*.25d0)
      big=int((bi(2)+bi(4)+bi(6)+bi(8))*.25d0)
      return
      end
      subroutine inblue4(sn,sni,bi,bii,bir,big,bib,red,green,blue)
C
C   demosaicing for blue pixel
C   sni contains the observed blue pixel
C   sn the sendor data from neighboring pixel (clockwise)
C
      implicit none
      integer sn(8),sni,red,green,blue,bi(8),bii,
     1       bir,big,bib
      blue=sni
      bib=bii
      red=int((sn(1)+sn(3)+sn(5)+sn(7))*.25d0)
      green=int((sn(2)+sn(4)+sn(6)+sn(8))*.25d0)
      bir=int((bi(1)+bi(3)+bi(5)+bi(7))*.25d0)
      big=int((bi(2)+bi(4)+bi(6)+bi(8))*.25d0)
      return
      end

      subroutine wbalance(sensor,n1,n2,wb,bayer)
      implicit none
      external channel
      integer n1,n2,sensor(n1,n2),bayer,channel,z
      double precision wb(3)
      integer i1,i2
      DO i1=1,n1
         DO i2=1,n2
            z=int(sensor(i1,i2)*wb(channel(i1,i2,bayer)))
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
      implicit none
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
      implicit none
      integer n,theta(n,3),thetanew(n,3)
      double precision outcam(3,3)
      integer i,j,k
      double precision z
      DO i=1,n
         DO j=1,3
            z=0.d0
            DO k=1,3
               z=z+theta(i,k)*outcam(k,j)
            END DO
            thetanew(i,j)=max(min(int(z),65535),0)
         END DO
      END DO
      RETURN
      END
