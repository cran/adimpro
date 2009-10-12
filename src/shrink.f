      subroutine median3(x,n,y,tol)
      implicit logical (a-z)
      integer n
      real*8 x(3,n),y(3),tol
      integer i,j
      real*8 etaofy,di1,di2,di3,rofy,dxy,z,delta,normy,
     1       y1,y2,y3,r1,r2,r3,t1,t2,t3,t0,c1,c2
C  use mean as init
      y1=x(1,1)
      y2=x(2,1)
      y3=x(3,1)
      DO i=2,n
         y1=y1+x(1,i)
         y2=y2+x(2,i)
         y3=y3+x(3,i)
      END DO
      y1=y1/n
      y2=y2/n
      y3=y3/n
C  iterate until convergence
      rofy = 1.d10
      j=0
      DO while (rofy.gt.tol)
C  compute r(y) and check for y=x_k 
         etaofy=0.d0
         r1=0.d0
         r2=0.d0
         r3=0.d0
         t0=0.d0
         t1=0.d0
         t2=0.d0
         t3=0.d0
         DO i=1,n
            di1=x(1,i)-y1
            di2=x(2,i)-y2
            di3=x(3,i)-y3
            dxy=sqrt(di1*di1+di2*di2+di3*di3)
            if(dxy.lt.1e-8) THEN
               etaofy=etaofy+1.d0
            ELSE
               r1=r1+di1/dxy
               r2=r2+di2/dxy
               r3=r3+di3/dxy
               t0=t0+1.d0/dxy
               t1=t1+x(1,i)/dxy
               t2=t2+x(2,i)/dxy
               t3=t3+x(3,i)/dxy
            END IF
         END DO      
         rofy=sqrt(r1*r1+r2*r2+r3*r3)
         if(rofy.le.tol) EXIT
         t1=t1/t0
         t2=t2/t0
         t3=t3/t0
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
         z=c1*t3+c2*y3
         delta=delta+abs(y3-z)
         normy=normy+abs(z)
         y3=z
         if(delta.lt.tol*normy) EXIT
         j=j+1
         if(j.gt.20) EXIT
      END DO
      y(1)=y1
      y(2)=y2
      y(3)=y3
      RETURN
      END
      subroutine median1(x,n,y,tol)
      implicit logical (a-z)
      integer n
      real*8 x(n),y,tol
      integer i,j
      real*8 etaofy,di1,rofy,dxy,r1,t1,t0,c1,c2,y0
C  use mean as init
      y=0.d0
      DO i=1,n
         y=y+x(i)
      END DO
      y=y/n
C  iterate until convergence
      rofy = 1.d10
      j=0
      y0=y
      DO while (rofy.gt.tol)
C  compute r(y) and check for y=x_k 
         etaofy=0.d0
         r1=0.d0
         t0=0.d0
         t1=0.d0
         DO i=1,n
            di1=x(i)-y
            dxy=abs(di1)
            if(dxy.lt.1e-8) THEN
               etaofy=etaofy+1.d0
            ELSE
               r1=r1+di1/dxy
               t0=t0+1.d0/dxy
               t1=t1+x(i)/dxy
            END IF
         END DO      
         rofy=abs(r1)
         if(rofy.le.tol) EXIT
         t1=t1/t0
         etaofy=etaofy/rofy
         c1=max(0.d0,1.d0-etaofy)
         c2=min(1.d0,etaofy)
         y=c1*t1+c2*y
         if(abs(y0-y).lt.tol*max(1.d0,y)) EXIT
         y0=y
         j=j+1
         if(j.gt.20) EXIT
      END DO
      RETURN
      END
      subroutine shrinkg(x,nx1,nx2,y,ny1,ny2,tol,z,nz,method)
      implicit logical (a-z)
      integer nx1,ny1,nx2,ny2,nz,x(nx1,nx2),y(ny1,ny2)
      real*8 z(nz),tol
      integer iy1,iy2,ja1,ja2,je1,je2,jx1,jx2,k,l,method
      real*8 yy,d1,d2
C
C   x - original image
C   y - new image
C   method = 1 use nearest observed value
C   method = 2 use (weighted) median of corresponding x pixel
C   method = 3 use weighted mean of corresponding x pixel
C
      d1=nx1
      d1=d1/ny1
C      rd1=d1*d1/4.d0
C   d1  contains the factor for shrinkage in first dimension
      d2=nx2
      d2=d2/ny2
C      rd2=d2*d2/4.d0
C   d1  contains the factor for shrinkage in second dimension
      DO iy1=1,ny1
         ja1=max(1.d0,0.5d0+(iy1-1)*d1)
         je1=max(1.d0,0.5d0+iy1*d1)
         je1=min(je1,nx1)
         DO iy2=1,ny2
            ja2=max(1.d0,0.5+(iy2-1)*d2)
            je2=max(1.d0,0.5d0+iy2*d2)
            je2=min(je2,nx2)
            if(ja1.eq.je1.and.ja2.eq.je2) THEN
               y(iy1,iy2)=x(ja1,ja2)
               CYCLE
            END IF
            if(method.eq.1) THEN
               jx1=(ja1+je1)/2
               jx2=(ja2+je2)/2
               y(iy1,iy2)=x(jx1,jx2)
            ELSE
               k=1
               DO jx1=ja1,je1
                  DO jx2=ja2,je2
                     z(k)=x(jx1,jx2)
                     k=k+1
                  END DO
               END DO
               k=k-1
               if(method.eq.2) THEN
                  call median1(z,k,yy,tol)
                  y(iy1,iy2)=yy
               END IF
               if(method.eq.3) THEN
                  yy = z(1)
                  DO l=2,k
                     yy=yy+z(l)
                  END DO
                  y(iy1,iy2)=yy/k
               END IF
            END IF
         END DO
      END DO
      RETURN
      END
      subroutine shrinkc(x,nx1,nx2,y,ny1,ny2,tol,z,nz,method)
      implicit logical (a-z)
      integer nx1,ny1,nx2,ny2,nz,x(nx1,nx2,3),y(ny1,ny2,3)
      real*8 z(3,nz),tol
      integer iy1,iy2,ja1,ja2,je1,je2,jx1,jx2,k,l,method
      real*8 yy(3),d1,d2
C
C   x - original image
C   y - new image
C   method = 1 use nearest observed value
C   method = 2 use (weighted) median of corresponding x pixel
C   method = 3 use weighted mean of corresponding x pixel
C
      d1=nx1
      d1=d1/ny1
C      rd1=d1*d1/4.d0
C   d1  contains the factor for shrinkage in first dimension
      d2=nx2
      d2=d2/ny2
C      rd2=d2*d2/4.d0
C   d1  contains the factor for shrinkage in second dimension
      DO iy1=1,ny1
         ja1=max(1.d0,0.5d0+(iy1-1)*d1)
         je1=max(1.d0,0.5d0+iy1*d1)
         je1=min(je1,nx1)
         DO iy2=1,ny2
            ja2=max(1.d0,0.5d0+(iy2-1)*d2)
            je2=max(1.d0,0.5d0+iy2*d2)
            je2=min(je2,nx2)
            if(ja1.eq.je1.and.ja2.eq.je2) THEN
               y(iy1,iy2,1)=x(ja1,ja2,1)
               y(iy1,iy2,2)=x(ja1,ja2,2)
               y(iy1,iy2,3)=x(ja1,ja2,3)
               CYCLE
            END IF
            if(method.eq.1) THEN
               jx1=(ja1+je1)/2
               jx2=(ja2+je2)/2
               y(iy1,iy2,1)=x(jx1,jx2,1)
               y(iy1,iy2,2)=x(jx1,jx2,2)
               y(iy1,iy2,3)=x(jx1,jx2,3)
C  this should be the central x point in the pixel
            ELSE
               k=1
               DO jx1=ja1,je1
                  DO jx2=ja2,je2
                     z(1,k)=x(jx1,jx2,1)
                     z(2,k)=x(jx1,jx2,2)
                     z(3,k)=x(jx1,jx2,3)
                     k=k+1
                  END DO
               END DO
               k=k-1
               if(method.eq.2) THEN
                  call median3(z,k,yy,tol)
                  y(iy1,iy2,1)=yy(1)
                  y(iy1,iy2,2)=yy(2)
                  y(iy1,iy2,3)=yy(3)
               END IF
               if(method.eq.3) THEN
                  yy(1) = z(1,1)
                  yy(2) = z(2,1)
                  yy(3) = z(3,1)
                  DO l=2,k
                     yy(1)=yy(1)+z(1,l)
                     yy(2)=yy(2)+z(2,l)
                     yy(3)=yy(3)+z(3,l)
                  END DO
                  y(iy1,iy2,1)=yy(1)/k
                  y(iy1,iy2,2)=yy(2)/k
                  y(iy1,iy2,3)=yy(3)/k
               END IF
            END IF
         END DO
      END DO
      RETURN
      END
