      subroutine convolve(img, kernel, eimg, height, width, ksize)

      integer height, width, ksize, i, j
      double precision img(width,height), kernel(ksize,ksize),
     1       eimg(width,height)
      double precision tmp1,tmp2,tmp3,tmp4,tmp5

      if (ksize.eq.5) then
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,tmp1,tmp2,tmp3,tmp4,tmp5)
C$OMP DO SCHEDULE(GUIDED)
         do j=3,height-2
            do i=3,width-2
               tmp1 = img(i-2,j-2) * kernel(1,1)
     1                + img(i-2,j-1) * kernel(1,2)
     2                + img(i-2,j) * kernel(1,3)
     3                + img(i-2,j+1) * kernel(1,4)
     4                + img(i-2,j+2) * kernel(1,5)
               tmp2 = img(i-1,j-2) * kernel(2,1)
     1                + img(i-1,j-1) * kernel(2,2)
     2                + img(i-1,j) * kernel(2,3)
     3                + img(i-1,j+1) * kernel(2,4)
     4                + img(i-1,j+2) * kernel(2,5)
               tmp3 = img(i,j-2) * kernel(3,1)
     1                + img(i,j-1) * kernel(3,2)
     2                + img(i,j) * kernel(3,3)
     3                + img(i,j+1) * kernel(3,4)
     4                + img(i,j+2) * kernel(3,5)
               tmp4 = img(i+1,j-2) * kernel(4,1)
     1                + img(i+1,j-1) * kernel(4,2)
     2                + img(i+1,j) * kernel(4,3)
     3                + img(i+1,j+1) * kernel(4,4)
     4                + img(i+1,j+2) * kernel(4,5)
               tmp5 = img(i+2,j-2) * kernel(5,1)
     1                + img(i+2,j-1) * kernel(5,2)
     2                + img(i+2,j) * kernel(5,3)
     3                + img(i+2,j+1) * kernel(5,4)
     4                + img(i+2,j+2) * kernel(5,5)
               eimg(i,j) = tmp1 + tmp2 + tmp3 +tmp4 + tmp5
            end do
         end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(eimg)
      else if (ksize.eq.3) then
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j,tmp1,tmp2,tmp3)
C$OMP DO SCHEDULE(GUIDED)
         do j=2,height-1
            do i=2,width-1
               tmp1 = img(i-1,j-1) * kernel(1,1)
     1                + img(i-1,j) * kernel(1,2)
     2                + img(i-1,j+1) * kernel(1,3)
               tmp2 = img(i,j-1) * kernel(2,1)
     1                + img(i,j) * kernel(2,2)
     2                + img(i,j+1) * kernel(2,3)
               tmp3 = img(i+1,j-1) * kernel(3,1)
     1                + img(i+1,j) * kernel(3,2)
     2                + img(i+1,j+1) * kernel(3,3)
               eimg(i,j) = tmp1 + tmp2 + tmp3
            end do
         end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(eimg)
      else if (ksize.eq.2) then
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(i,j)
C$OMP DO SCHEDULE(GUIDED)
         do j=1,height-1
            do i=1,width-1
               eimg(i,j) = img(i,j) * kernel(1,1)
     1              + img(i,j+1) * kernel(1,2)
     2              + img(i+1,j) * kernel(2,1)
     3              + img(i+1,j+1) * kernel(2,2)
            end do
         end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(eimg)
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
      implicit none
      integer nx,ny,dv,img(nx,ny,dv),nxnew,nynew,
     1        imgnew(nxnew,nynew,dv),indx(*),indy(*),method
      integer i,j,inew,jnew,k,nij,ibest,jbest
      double precision z,znew,gap,zmean(4),dist,bestdist
C
C     First generate index vectors
C
      z = nx
      znew = nxnew
      gap = z/znew
      indx(1)=1
      DO inew=2,nxnew
         indx(inew)=int((inew-1)*gap+1)
      END DO
      indx(nxnew+1)=nx+1
      z = ny
      znew = nynew
      gap = z/znew
      indy(1)=1
      DO inew=2,nynew
         indy(inew)=int((inew-1)*gap+1)
      END DO
      indy(nynew+1)=ny+1
C
C     Now fill imgnew
C
      IF(method.eq.1) THEN
C
C       select representative (central) pixel
C
         DO jnew=1,nynew
            DO inew=1,nxnew
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
         DO jnew=1,nynew
            DO inew=1,nxnew
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
                  imgnew(inew,jnew,k)=int(zmean(k)/nij)
               END DO
            END DO
         END DO
      END IF
      IF(method.eq.3) THEN
C
C       select pixel most similar to the mean
C
         DO jnew=1,nynew
            DO inew=1,nxnew
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
      implicit none
      integer nx,ny,dv,nxnew,nynew,indx(*),indy(*),method
      double precision img(nx,ny,dv),imgnew(nxnew,nynew,dv)
      integer i,j,inew,jnew,k,nij,ibest,jbest
      double precision z,znew,gap,zmean(4),dist,bestdist
C
C     First generate index vectors
C
      z = nx
      znew = nxnew
      gap = z/znew
      indx(1)=1
      DO inew=2,nxnew
         indx(inew)=int((inew-1)*gap+1)
      END DO
      indx(nxnew+1)=nx+1
      z = ny
      znew = nynew
      gap = z/znew
      indy(1)=1
      DO inew=2,nynew
         indy(inew)=int((inew-1)*gap+1)
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
      implicit none
      integer nx,ny,img(nx,ny),nxnew,nynew,
     1        imgnew(nxnew,nynew),indx(*),indy(*),method
      integer i,j,inew,jnew,nij,ibest,jbest
      double precision z,znew,gap,zmean,dist,bestdist
C
C     First generate index vectors
C
      z = nx
      znew = nxnew
      gap = z/znew
      indx(1)=1
      DO inew=2,nxnew
         indx(inew)=int((inew-1)*gap+1)
      END DO
      indx(nxnew+1)=nx+1
      z = ny
      znew = nynew
      gap = z/znew
      indy(1)=1
      DO inew=2,nynew
         indy(inew)=int((inew-1)*gap+1)
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
               imgnew(inew,jnew)=int(zmean/nij)
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
