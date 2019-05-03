      subroutine mdout
      implicit double precision (a-h,o-z)
      common /del1/dmolin(169),dofscl(169)
      common /del2/d11(3,3),d12(3,3),d21(3,3),d22(3,3),
     &             f11(3),f12(3),f21(3),f22(3),occ(3),
     &             p(3,3),g(3,3),ibasis,idim

      do i=1,3
         k1 = (idim - 8 + i) * idim
         k2 = (idim - 4 + i) * idim
         do j=1,3
            l1 =idim - 7 + j
            l2 =l1 + 4
            dofscl(k1+l1) = d11(i,j)
            dofscl(k1+l2) = d12(i,j)
            dofscl(k2+l1) = d21(i,j)
            dofscl(k2+l2) = d22(i,j)
         end do
      end do

      return
      end
