      subroutine rotden(p,occ)
      implicit double precision (a-h,o-z)
      common/rotxyz/amat(3,3),bmat(3,3),cmat(3,3),rmat(3,3),todeg
      dimension p(3,3),occ(3)

      do i=1,3
          do j=i,3
             sum = 0.d0
             do k=1,3
                 sum = sum + occ(k)*rmat(k,i)*rmat(k,j)
             end do
             p(j,i)=sum
             p(i,j)=sum
          end do
      end do

      return
      end
