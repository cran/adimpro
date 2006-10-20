      subroutine convolve(img, kernel, eimg, height, width, ksize)

      integer height, width, ksize, i, j
      real*8 img(width,height), kernel(ksize,ksize), eimg(width,height)
      real*8 tmp(5)

      if (ksize.eq.5) then
         do j=3,height-2
            do i=3,width-2
               tmp(1) = img(i-2,j-2) * kernel(1,1)
     1                + img(i-2,j-1) * kernel(1,2)
     2                + img(i-2,j) * kernel(1,3)
     3                + img(i-2,j+1) * kernel(1,4)
     4                + img(i-2,j+2) * kernel(1,5)
               tmp(2) = img(i-1,j-2) * kernel(2,1)
     1                + img(i-1,j-1) * kernel(2,2)
     2                + img(i-1,j) * kernel(2,3)
     3                + img(i-1,j+1) * kernel(2,4)
     4                + img(i-1,j+2) * kernel(2,5)
               tmp(3) = img(i,j-2) * kernel(3,1)
     1                + img(i,j-1) * kernel(3,2)
     2                + img(i,j) * kernel(3,3)
     3                + img(i,j+1) * kernel(3,4)
     4                + img(i,j+2) * kernel(3,5)
               tmp(4) = img(i+1,j-2) * kernel(4,1)
     1                + img(i+1,j-1) * kernel(4,2)
     2                + img(i+1,j) * kernel(4,3)
     3                + img(i+1,j+1) * kernel(4,4)
     4                + img(i+1,j+2) * kernel(4,5)
               tmp(5) = img(i+2,j-2) * kernel(5,1)
     1                + img(i+2,j-1) * kernel(5,2)
     2                + img(i+2,j) * kernel(5,3)
     3                + img(i+2,j+1) * kernel(5,4)
     4                + img(i+2,j+2) * kernel(5,5)
               eimg(i,j) = tmp(1) + tmp(2) + tmp(3) +tmp(4) + tmp(5)
            end do
         end do
      else if (ksize.eq.3) then 
         do j=2,height-1
            do i=2,width-1
               tmp(1) = img(i-1,j-1) * kernel(1,1)
     1                + img(i-1,j) * kernel(1,2)
     2                + img(i-1,j+1) * kernel(1,3)
               tmp(2) = img(i,j-1) * kernel(2,1)
     1                + img(i,j) * kernel(2,2)
     2                + img(i,j+1) * kernel(2,3)
               tmp(3) = img(i+1,j-1) * kernel(3,1)
     1                + img(i+1,j) * kernel(3,2)
     2                + img(i+1,j+1) * kernel(3,3)
               eimg(i,j) = tmp(1) + tmp(2) + tmp(3)
            end do
         end do

      else if (ksize.eq.2) then
         do j=1,height-1
            do i=1,width-1
               eimg(i,j) = img(i,j) * kernel(1,1) 
     1              + img(i,j+1) * kernel(1,2) 
     2              + img(i+1,j) * kernel(2,1) 
     3              + img(i+1,j+1) * kernel(2,2)
            end do
         end do
      end if


      return
      end


      subroutine shrnkrgb(img,nx,ny,dv,imgnew,nxnew,nynew,indx,indy,
     1                    method)
C
C   shrink an RGB image
C      
C   indx, indy  -  index vectors of length nxnew+1 and nynew+1
C
      implicit logical (a-z)
      integer nx,ny,dv,img(nx,ny,dv),nxnew,nynew,
     1        imgnew(nxnew,nynew,dv),indx(1),indy(1),method
      integer i,j,inew,jnew,k,nij,ibest,jbest
      real*8 z,znew,gap,zmean(4),dist,bestdist
C     
C     First generate index vectors
C 
      z = nx
      znew = nxnew
      gap = z/znew
      indx(1)=1
      DO inew=2,nxnew
         indx(inew)=(inew-1)*gap+1
      END DO
      indx(nxnew+1)=nx+1
      z = ny
      znew = nynew
      gap = z/znew
      indy(1)=1
      DO inew=2,nynew
         indy(inew)=(inew-1)*gap+1
      END DO
      indy(nynew+1)=ny+1
C
C     Now fill imgnew
C
      IF(method.eq.1) THEN
C
C       select representative (central) pixel
C
         DO inew=1,nxnew
            DO jnew=1,nynew
               i=(indx(inew)+indx(inew+1)-1)/2
               j=(indy(jnew)+indy(jnew+1)-1)/2
               DO k=1,dv               
                  imgnew(inew,jnew,k)=img(i,j,k)
               END DO
            END DO
         END DO
      END IF
      IF(method.eq.2) THEN
C
C       select pixel as the mean
C
         DO inew=1,nxnew
            DO jnew=1,nynew
               nij=0
               DO k=1,dv
                  zmean(k)=0.d0
               END DO
               DO i=indx(inew),indx(inew+1)-1
                  DO j=indy(jnew),indy(jnew+1)-1
                     nij=nij+1
                     DO k=1,dv
                        zmean(k)=zmean(k)+img(i,j,k)
                     END DO
                  END DO
               END DO
                DO k=1,dv
                  imgnew(inew,jnew,k)=zmean(k)/nij
               END DO
            END DO
         END DO
      END IF
      IF(method.eq.3) THEN
C
C       select pixel most similar to the mean
C
         DO inew=1,nxnew
            DO jnew=1,nynew
               nij=0
               DO k=1,dv
                  zmean(k)=0.d0
               END DO
               DO i=indx(inew),indx(inew+1)-1
                  DO j=indy(jnew),indy(jnew+1)-1
                     nij=nij+1
                     DO k=1,dv
                        zmean(k)=zmean(k)+img(i,j,k)
                     END DO
                  END DO
               END DO
               DO k=1,dv
                  zmean(k)=zmean(k)/nij
               END DO
               bestdist=1.d40
               jbest=1
               ibest=1
               DO  i=indx(inew),indx(inew+1)-1
                  DO j=indy(jnew),indy(jnew+1)-1
                     dist=0.d0
                     DO k=1,dv
                        dist=dist+dabs(img(i,j,k)-zmean(k))
                     END DO
                     IF(dist.lt.bestdist) THEN
                         ibest=i
                         jbest=j
                         bestdist=dist
                     END IF
                  END DO
               END DO
               DO k=1,dv
                  imgnew(inew,jnew,k)=img(ibest,jbest,k)
               END DO
            END DO
         END DO
      END IF
      RETURN
      END
      subroutine shrnkcsp(img,nx,ny,dv,imgnew,nxnew,nynew,indx,indy,
     1                    method)
C
C   shrink an RGB image
C      
C   indx, indy  -  index vectors of length nxnew+1 and nynew+1
C
      implicit logical (a-z)
      integer nx,ny,dv,nxnew,nynew,indx(1),indy(1),method
      real*8 img(nx,ny,dv),imgnew(nxnew,nynew,dv)
      integer i,j,inew,jnew,k,nij,ibest,jbest
      real*8 z,znew,gap,zmean(4),dist,bestdist
C     
C     First generate index vectors
C 
      z = nx
      znew = nxnew
      gap = z/znew
      indx(1)=1
      DO inew=2,nxnew
         indx(inew)=(inew-1)*gap+1
      END DO
      indx(nxnew+1)=nx+1
      z = ny
      znew = nynew
      gap = z/znew
      indy(1)=1
      DO inew=2,nynew
         indy(inew)=(inew-1)*gap+1
      END DO
      indy(nynew+1)=ny+1
C
C     Now fill imgnew
C
      IF(method.eq.1) THEN
C
C       select representative (central) pixel
C
         DO inew=1,nxnew
            DO jnew=1,nynew
               i=(indx(inew)+indx(inew+1)-1)/2
               j=(indy(jnew)+indy(jnew+1)-1)/2
               DO k=1,dv               
                  imgnew(inew,jnew,k)=img(i,j,k)
               END DO
            END DO
         END DO
      END IF
      IF(method.eq.2) THEN
C
C       select pixel as the mean
C
         DO inew=1,nxnew
            DO jnew=1,nynew
               nij=0
               DO k=1,dv
                  zmean(k)=0.d0
               END DO
               DO i=indx(inew),indx(inew+1)-1
                  DO j=indy(jnew),indy(jnew+1)-1
                     nij=nij+1
                     DO k=1,dv
                        zmean(k)=zmean(k)+img(i,j,k)
                     END DO
                  END DO
               END DO
                DO k=1,dv
                  imgnew(inew,jnew,k)=zmean(k)/nij
               END DO
            END DO
         END DO
      END IF
      IF(method.eq.3) THEN
C
C       select pixel most similar to the mean
C
         DO inew=1,nxnew
            DO jnew=1,nynew
               nij=0
               DO k=1,dv
                  zmean(k)=0.d0
               END DO
               DO i=indx(inew),indx(inew+1)-1
                  DO j=indy(jnew),indy(jnew+1)-1
                     nij=nij+1
                     DO k=1,dv
                        zmean(k)=zmean(k)+img(i,j,k)
                     END DO
                  END DO
               END DO
               DO k=1,dv
                  zmean(k)=zmean(k)/nij
               END DO
               jbest=1
               ibest=1
               bestdist=1.d40
               DO  i=indx(inew),indx(inew+1)-1
                  DO j=indy(jnew),indy(jnew+1)-1
                     dist=0.d0
                     DO k=1,dv
                        dist=dist+dabs(img(i,j,k)-zmean(k))
                     END DO
                     IF(dist.lt.bestdist) THEN
                         ibest=i
                         jbest=j
                         bestdist=dist
                     END IF
                  END DO
               END DO
               DO k=1,dv               
                  imgnew(inew,jnew,k)=img(ibest,jbest,k)
               END DO
            END DO
         END DO
      END IF
      RETURN
      END
      subroutine shrnkgr(img,nx,ny,imgnew,nxnew,nynew,indx,indy,
     1                    method)
C
C   shrink an RGB image
C      
C   indx, indy  -  index vectors of length nxnew+1 and nynew+1
C
      implicit logical (a-z)
      integer nx,ny,img(nx,ny),nxnew,nynew,
     1        imgnew(nxnew,nynew),indx(1),indy(1),method
      integer i,j,inew,jnew,nij,ibest,jbest
      real*8 z,znew,gap,zmean,dist,bestdist
C     
C     First generate index vectors
C 
      z = nx
      znew = nxnew
      gap = z/znew
      indx(1)=1
      DO inew=2,nxnew
         indx(inew)=(inew-1)*gap+1
      END DO
      indx(nxnew+1)=nx+1
      z = ny
      znew = nynew
      gap = z/znew
      indy(1)=1
      DO inew=2,nynew
         indy(inew)=(inew-1)*gap+1
      END DO
      indy(nynew+1)=ny+1
C
C     Now fill imgnew
C
      IF(method.eq.1) THEN
C
C       select representative (central) pixel
C
         DO inew=1,nxnew
            DO jnew=1,nynew
               i=(indx(inew)+indx(inew+1)-1)/2
               j=(indy(jnew)+indy(jnew+1)-1)/2
               imgnew(inew,jnew)=img(i,j)
            END DO
         END DO
      END IF
      IF(method.eq.2) THEN
C
C       select the mean
C
         DO inew=1,nxnew
            DO jnew=1,nynew
               nij=0
               zmean=0.d0
               DO i=indx(inew),indx(inew+1)-1
                  DO j=indy(jnew),indy(jnew+1)-1
                     nij=nij+1
                     zmean=zmean+img(i,j)
                  END DO
               END DO
               imgnew(inew,jnew)=zmean/nij
            END DO
         END DO
      END IF
      IF(method.eq.3) THEN
C
C       select the pixel closest to the mean
C
         DO inew=1,nxnew
            DO jnew=1,nynew
               nij=0
               zmean=0.d0
               DO i=indx(inew),indx(inew+1)-1
                  DO j=indy(jnew),indy(jnew+1)-1
                     nij=nij+1
                     zmean=zmean+img(i,j)
                  END DO
               END DO
               zmean=zmean/nij
               bestdist=1.d40
               jbest=1
               ibest=1
               DO  i=indx(inew),indx(inew+1)-1
                  DO j=indy(jnew),indy(jnew+1)-1
                     dist=dabs(img(i,j)-zmean)
                     IF(dist.lt.bestdist) THEN
                         ibest=i
                         jbest=j
                         bestdist=dist
                     END IF
                  END DO
               END DO
               imgnew(inew,jnew)=img(ibest,jbest)
            END DO
         END DO
      END IF
      RETURN
      END
