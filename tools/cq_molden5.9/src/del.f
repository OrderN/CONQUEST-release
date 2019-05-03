      double precision function del(var)
      implicit double precision (a-h,o-z)
      common /del1/dmolin(169),dofscl(169)
      common /del2/d11(3,3),d12(3,3),d21(3,3),d22(3,3),
     &             f11(3),f12(3),f21(3),f22(3),occ(3),
     &             p(3,3),g(3,3),ibasis,idim
      dimension var(2)

      call rota(var(1))
      call rotb(var(2))
      call rotm

      if (ibasis.eq.0) then

         call rotden(p,occ)
         del=0.0d0

         do i=1,3
            do j=1,3
               del=del+(p(i,j)-g(i,j))**2
            end do
         end do

      else

         call rotden(d11,f11)
         call rotden(d12,f12)
         call rotden(d21,f21)
         call rotden(d22,f22)
         call mdout
         del=0.0d0

         do i=1,idim
            do j=1,idim
                 iaddr=j+(i-1)*idim
                 del=del+(dofscl(iaddr)-dmolin(iaddr))**2
            end do
         end do

      endif

      return
      end
