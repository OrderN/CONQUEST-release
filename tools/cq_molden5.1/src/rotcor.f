      subroutine rotcod(b,coo)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat

      common/rotxyz/amat(3,3),bmat(3,3),cmat(3,3),rmat(3,3),todeg
      dimension temp(3,3),b(3,*),coo(3,*)

      do i=1,3
          do j=1,3
             temp(i,j) = amat(i,1)*bmat(1,j) +
     &                   amat(i,2)*bmat(2,j) +
     &                   amat(i,3)*bmat(3,j) 
          end do
      end do

      do i=1,3
          do j=1,3
             rmat(i,j) = temp(i,1)*cmat(1,j) +
     &                   temp(i,2)*cmat(2,j) +
     &                   temp(i,3)*cmat(3,j) 
          end do
      end do

      do i=1,iatoms
          do j=1,3
             b(j,i) = rmat(j,1)*coo(1,i) + rmat(j,2)*coo(2,i)
     &                + rmat(j,3)*coo(3,i)
          end do
      end do

      return
      end
