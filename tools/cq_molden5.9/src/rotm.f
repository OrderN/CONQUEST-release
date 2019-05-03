      subroutine rotm
      implicit double precision (a-h,o-z)
      common/rotxyz/amat(3,3),bmat(3,3),cmat(3,3),rmat(3,3),todeg

      do i=1,3
         do j=1,3
           rmat(i,j) = 0.0d0
         end do
      end do
      do i=1,3
         do j=1,3
            do k=1,3
               rmat(i,j)=rmat(i,j)+amat(i,k)*bmat(k,j)
            end do
         end do
      end do

      return
      end
