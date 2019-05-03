      subroutine rotatg(nat,idebug)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /rotxyz/ amat(3,3),bmat(3,3),cmat(3,3),rmat(3,3),todeg
      common /orient/ o321g(9,9),o431g(9,9),o631g(9,9),f321g(9,9),
     &                f431g(9,9),f631g(9,9),s321g(13,13),s431g(13,13),
     &                s631g(13,13),cl321g(13,13),cl431g(13,13),
     &                cl631g(13,13)
      common /orihlp/ ori(numatm),iuser(numatm),oalpha(numatm),
     $                obeta(numatm),norien,ibal
      common /del1/   dmolin(169),dofscl(169)
      common /del2/   d11(3,3),d12(3,3),d21(3,3),d22(3,3),
     &                f11(3),f12(3),f21(3),f22(3),occ(3),
     &                p(3,3),g(3,3),ibasis,idim
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension angle(2)

c###### split valence part ###############
c     dmolin contains the atomic part off the molecular
c            density matrix
c     dofscl   contains the oriented atomic density matrix
c            on return
c     idim   dimension of atomic density matrix without
c            polarisation function (=9 for o and f,
c                                   =13 for s and cl )
c     d11,d12,d21,d22
c            contain the p*p p*p' p'*p and p'*p' parts
c            of the atomic density matrix
c     f11,f12,f21,f22
c            idem,but now only diagonal terms
c
c--- initialisation -----
c
      todeg = 45.0d0 / datan(1.0d0)
      delold = 1.0d5

      amat(1,3) = 0.0d0
      amat(2,3) = 0.0d0
      amat(3,1) = 0.0d0
      amat(3,2) = 0.0d0
      amat(3,3) = 1.0d0
      bmat(1,2) = 0.0d0
      bmat(2,1) = 0.0d0
      bmat(2,2) = 1.0d0
      bmat(2,3) = 0.0d0
      bmat(3,2) = 0.0d0

      if (ibasis.eq.0) then

c----- sto3g ----------

         if (nat.eq.8.or.nat.eq.16) then
            occ(1) = 2.0d0
            occ(2) = 1.0d0
            occ(3) = 1.0d0
         endif
         if (nat.eq.9.or.nat.eq.17) then
            occ(1) = 2.0d0
            occ(2) = 2.0d0
            occ(3) = 1.0d0
         endif

      else

c------- spit valence -------

         if (nat.eq.8) then

           if (ibasis.eq.1) then

              do i=1,9
                 do j=1,9
                    dofscl(j+(i-1)*idim) = o321g(i,j)
                 end do
              end do

           elseif (ibasis.eq.2) then

              do i=1,9
                 do j=1,9
                    dofscl(j+(i-1)*idim) = o431g(i,j)
                 end do
              end do

           elseif (ibasis.eq.3) then

              do i=1,9
                 do j=1,9
                    dofscl(j+(i-1)*idim) = o631g(i,j)
                 end do
              end do

           endif

         endif

         if (nat.eq.9) then
           if (ibasis.eq.1) then

              do i=1,9
                 do j=1,9
                    dofscl(j+(i-1)*idim) = f321g(i,j)
                 end do
              end do

           elseif (ibasis.eq.2) then

              do i=1,9
                 do j=1,9
                    dofscl(j+(i-1)*idim) = f431g(i,j)
                 end do
              end do

           elseif (ibasis.eq.3) then

              do i=1,9
                 do j=1,9
                    dofscl(j+(i-1)*idim) = f631g(i,j)
                 end do
              end do

           endif
         endif

         if (nat.eq.16) then

           if (ibasis.eq.1) then

              do i=1,13
                 do j=1,13
                    dofscl(j+(i-1)*idim) = s321g(i,j)
                 end do
              end do

           elseif (ibasis.eq.2) then

              do i=1,13
                 do j=1,13
                    dofscl(j+(i-1)*idim) = s431g(i,j)
                 end do
              end do

           elseif (ibasis.eq.3) then

              do i=1,13
                 do j=1,13
                    dofscl(j+(i-1)*idim) = s631g(i,j)
                 end do
              end do

           endif

         endif

         if (nat.eq.17) then

           if (ibasis.eq.1) then

              do i=1,13
                 do j=1,13
                    dofscl(j+(i-1)*idim) = cl321g(i,j)
                 end do
              end do

           elseif (ibasis.eq.2) then

              do i=1,13
                 do j=1,13
                    dofscl(j+(i-1)*idim) = cl431g(i,j)
                 end do
              end do

           elseif (ibasis.eq.3) then

              do i=1,13
                 do j=1,13
                    dofscl(j+(i-1)*idim) = cl631g(i,j)
                 end do
              end do

           endif

         endif

         idim2=idim*idim
         do i=1,3
            f11(i) = dofscl(idim2+(i-7)*(idim+1))
            f12(i) = dofscl(idim2+(i-3)*idim-7+i)
            f21(i) = f12(i)
            f22(i) = dofscl(idim2+(i-3)*(idim+1))
         end do
      endif

c---- end initialisation -----

      write(iun3,*)' '
      write(iun3,*)'p-part atomic density matrix before orientation'
      write(iun3,*)' '

      if (ibasis.ne.0) then

          do i=1,3
             k1=(idim-8+i)*idim
             k2=(idim-4+i)*idim
             write(iun3,130)(dofscl(k1+idim-7+j),j=1,3),
     &                   (dofscl(k1+idim-3+j),j=1,3)
          end do

          write(iun3,*)' '

          do i=1,3
             k1=(idim-8+i)*idim
             k2=(idim-4+i)*idim
             write(iun3,130)(dofscl(k2+idim-7+j),j=1,3),
     &                   (dofscl(k2+idim-3+j),j=1,3)
          end do

      else

          write(iun3,'(f9.4,a)')occ(1),'   0.0000   0.0000'
          write(iun3,'(a,f9.4,a)')'   0.0000',occ(2),'   0.0000'
          write(iun3,'(a,f9.4)')'   0.0000   0.0000',occ(3)

      endif

      write(iun3,*)' '

      if (iuser(ibal).eq.1) then
         angle(1) = oalpha(ibal)
         angle(2) = obeta(ibal)
         delta    = del(angle)
         goto 1400
      endif

c----- start of coarse search ----

      do ia=0,180,10

         angle(1) = dfloat(ia)

         do ib=0,180,10

            angle(2) = dfloat(ib)
            delta = del(angle)

            if (delta.lt.delold) then
                 delold = delta
                 alopt  = angle(1)
                 beopt  = angle(2)
            endif

         end do

      end do

      angle(1) = alopt
      angle(2) = beopt

c----- end of coarse search ------

c----- start of medium search ----

      alco = alopt
      beco = beopt

      do ia=1,21

         angle(1) = alco + dfloat(ia-1)-10.0d0

         do ib=1,21

            angle(2) = beco + dfloat(ib-1)-10.0d0
            delta    = del(angle)

            if (delta.lt.delold) then
              delold = delta
              alopt  = angle(1)
              beopt  = angle(2)
            endif

         end do

       end do

c----- end of medium search ------

c----- start of fine search ----

      alco = alopt
      beco = beopt

      do ia=1,21

         angle(1) = alco + dfloat(ia-1)*0.1d0-1.0d0

         do ib=1,21

            angle(2) = beco + dfloat(ib-1)*0.1d0-1.0d0
            delta    = del(angle)

            if (delta.lt.delold) then
              delold = delta
              alopt  = angle(1)
              beopt  = angle(2)
            endif

         end do

      end do

c----- end of fine search ------

      angle(1) = alopt
      angle(2) = beopt
      delta    = del(angle)

1400  continue

      write(iun3,*)' '
      write(iun3,*)'delta squared                  =',delta
      write(iun3,*)'alfa  optimised                =',angle(1)
      write(iun3,*)'beta  optimised                =',angle(2)
      write(iun3,*)' '
      write(iun3,*)'-------- rot-matrix ----------'
      write(iun3,*)' '

      do i=1,3
         write(iun3,120)(rmat(i,j),j=1,3)
      end do

      write(iun3,*)' '

      if (ibasis.eq.0) then

          write(iun3,*)'p-part molecular density matrix'
          write(iun3,*)' '

          do i=1,3
             write(iun3,120)(g(i,j),j=1,3)
          end do

          write(iun3,*)' '
          write(iun3,*)'p-part oriented atomic density matrix'
          write(iun3,*)' '

          do i=1,3
             write(iun3,120)(p(i,j),j=1,3)
          end do

      else

          write(iun3,*)'p-part molecular density matrix'
          write(iun3,*)' '

          do i=1,3
             k1=(idim-8+i)*idim
             k2=(idim-4+i)*idim
             write(iun3,130)(dmolin(k1+idim-7+j),j=1,3),
     &                      (dmolin(k1+idim-3+j),j=1,3)
          end do

          write(iun3,*)' '

          do i=1,3
             k1=(idim-8+i)*idim
             k2=(idim-4+i)*idim
             write(iun3,130)(dmolin(k2+idim-7+j),j=1,3),
     &                   (dmolin(k2+idim-3+j),j=1,3)
          end do

          write(iun3,*)' '
          write(iun3,*)'p-part oriented atomic density matrix'
          write(iun3,*)' '

          do i=1,3
             write(iun3,130)(d11(i,j),j=1,3),(d12(i,j),j=1,3)
          end do

          write(iun3,*)' '

          do i=1,3
             write(iun3,130)(d21(i,j),j=1,3),(d22(i,j),j=1,3)
          end do

          write(iun3,*)' '
      endif

120   format(10x,3f9.4)
130   format(10x,3f9.4,2x,3f9.4)

      return
      end
