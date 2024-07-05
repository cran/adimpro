      subroutine getvofh2(bw,kern,vol)
      implicit none
      integer kern
      double precision bw,vol,sofw2D
      external sofw2D
      vol=sofw2D(bw,kern)
      RETURN
      END
      double precision function sofw2D(bw,kern)
      implicit none
      integer kern
      double precision bw
      integer j1,j2,clw,ih1,ih
      double precision sw,sw2,h2,lkern,z1,z2,z
      external lkern
      h2=bw*bw
C
C   first calculate location weights
C
      ih=int(bw)
      clw=ih+1
      sw=0.d0
      sw2=0.d0
      DO j2=clw-ih,clw+ih
         z2=(clw-j2)
         z2=z2*z2
         ih1=int(sqrt(h2-z2))
         DO j1=clw-ih1,clw+ih1
            z1=clw-j1
            z=lkern(kern,(z1*z1+z2)/h2)
            sw=sw+z
            sw2=sw2+z*z
         END DO
      END DO
      sofw2D=sw*sw/sw2
      RETURN
      END
      subroutine geth2(x,y,kern,value,eps,bw)
      implicit none
      integer kern
      double precision x,y,value,eps,bw
      double precision fw1,fw2,fw3,z
      double precision sofw2D
      external sofw2D
      if(x.ge.y) RETURN
      fw1=sofw2D(x,kern)
      fw2=sofw2D(y,kern)
      DO WHILE(fw1.gt.value)
         x=x*x/y
         fw1=sofw2D(x,kern)
      END DO
      DO WHILE(fw2.le.value)
         y=y*y/x
         fw2=sofw2D(y,kern)
      END DO
      DO WHILE(min(fw2/value,value/fw1).gt.1.d0+eps)
         z=x+(value-fw1)/(fw2-fw1)*(y-x)
         fw3=sofw2D(z,kern)
         if(fw3.le.value) THEN
            x=z
            fw1=fw3
         ENDIF
         if(fw3.ge.value) THEN
            y=z
            fw2=fw3
         ENDIF
      END DO
      if(fw2/value.gt.value/fw1) THEN
          bw=x+(value-fw1)/(fw2-fw1)*(y-x)
      ELSE
          bw=y-(fw2-value)/(fw2-fw1)*(y-x)
      ENDIF
      RETURN
      END
      subroutine hequalg(x,n,y,yi)
      integer n,x(n),xi(65536),y(n),yi(65536)
      integer i
      double precision z, sz
      DO i=1,65536
         xi(i)=0
      END DO
      DO i=1,n
         j=x(i)+1
         xi(j)=xi(j)+1
      END DO
      sz=0.d0
      z=n
      z=65536/z
      DO i=1,65536
         sz=sz+xi(i)
         yi(i)=min(int(sz*z),65535)
      END DO
      DO i=1,n
         j=x(i)+1
         y(i)=yi(j)
      END DO
      RETURN
      END
      subroutine cumhist(x,n,yi)
      integer n,x(n),xi(65536),yi(65536)
      integer i
      double precision z, sz
      DO i=1,65536
         xi(i)=0
      END DO
      DO i=1,n
         j=x(i)+1
         xi(j)=xi(j)+1
      END DO
      sz=0.d0
      z=n
      z=65536/z
      DO i=1,65536
         sz=sz+xi(i)
         yi(i)=min(int(sz*z),65535)
      END DO
      RETURN
      END
      subroutine hequalc(x,n,y,yi)
      integer n,x(n,3),y(n,3),yi(65536)
      integer i,j,k
      DO i=1,n
         DO k=1,3
            j=x(i,k)+1
            y(i,k)=yi(j)
         END DO
      END DO
      RETURN
      END
      subroutine ihequal(x,n,y,yi)
      integer n,x(n),xi(65536),y(n),yi(65536)
      integer i,j,k
      k=1
      DO i=1,65536
         if(k.le.yi(i)) THEN
            DO j=k,yi(i)
               xi(j)=i-1
            END DO
            k=yi(i)+1
         END IF
      END DO
      IF(k.le.65536) THEN
         DO j=k,65536
            xi(j)=65535
         END DO
      END IF
      DO i=1,n
         y(i)=xi(x(i))
      END DO
      RETURN
      END
      subroutine ihequalc(x,n,y,yi)
      integer n,x(n,3),xi(65536),y(n,3),yi(65536)
      integer i,j,k
      k=1
      DO i=1,65536
         if(k.le.yi(i)) THEN
            DO j=k,yi(i)
               xi(j)=i-1
            END DO
            k=yi(i)+1
         END IF
      END DO
      IF(k.le.65536) THEN
         DO j=k,65536
            xi(j)=65535
         END DO
      END IF
      DO i=1,n
         DO k=1,3
            y(i,k)=xi(x(i,k))
         END DO
      END DO
      RETURN
      END
