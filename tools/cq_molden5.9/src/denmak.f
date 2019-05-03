      subroutine denmdd(ido,
     &                  vectrs,vectrb,focc,focb,p,paa,averag,q)
c THIS IS REALLY denmak
c     returns density matrix in p common /densty/
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      character bstrng(6)*8
      logical valenc,bonds,ovrlap,atomic,doori,orient,denok,dolap,
     &        doelf
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /orbhlp/ mxorb,iuhf,ispd
      common /atoms/  nbeg(68),natn(68),atmd(5440)
      common /bond/   nbas(numatm),ipol(numatm)
      common /orihlp/ ori(numatm),iuser(numatm),oalpha(numatm),
     $                obeta(numatm),norien,ibal
      common /del1/   dmolin(169),dofscl(169)
      common /del2/   d11(3,3),d12(3,3),d21(3,3),d22(3,3),
     &                f11(3),f12(3),f21(3),f22(3),occ(3),
     &                pr(3,3),ph(3,3),ibasis,idim
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap,doelf
      common /basopt/ extbas
      common /elem/elemnt(mxel)
      logical extbas
      character*2 elemnt
      dimension ioffst(4)
      dimension vectrs(*),vectrb(*),focc(*),focb(*),p(*),paa(*),
     &          averag(*),q(*)

      denok = .true.
      ido = 1
      if (ipsi.ne.0) then
         write(iun3,20)ipsi
20       format(//10x,'M.O. NUMBER ',i4,'  IS TO BE PLOTTED',/)
      endif
      if (ipsi.eq.0) call inferr('Calculating Density Matrix',0)

      ioffst(1) = 0
      ioffst(2) = 769
      ioffst(3) = 2608
      ioffst(4) = 3605

      bstrng(1) = 'sto3g   '
      bstrng(2) = '3-21g   '
      bstrng(3) = '4-31g   '
      bstrng(4) = '6-31g   '
      bstrng(5) = '4-31g*  '
      bstrng(6) = '6-31g*  '

c************ distribute electrons over mos *******************
*
      if (valenc) then
        if (iuhf.eq.1) call inferr('dont use VALENCE with uhf',1)
        if (iftyp.eq.1) call inferr('dont use VALENCE with MOPAC',1)
        izel = 0
        do j=1,natoms
          if (nat(j).lt.11.and.nat(j).gt.2)  izel = izel + 1
          if (nat(j).lt.19.and.nat(j).gt.10) izel = izel + 5
          if (nat(j).gt.18) 
     &      call inferr('Use VALENCE ONLY with atoms up to ARGON',1)
        end do
        do j=1,izel
          focc(j) = 0.0d0
        end do
      endif

      if (idebug.eq.1) then
         write(iun3,*)'----------------------------------------------'
         write(iun3,*)'Orbital Occupancies'
         write(iun3,*)' '

         if (iuhf.eq.1) then
           write(iun3,*)' '
           write(iun3,*)'Alpha set'
           write(iun3,*)' '
         endif

         do i=1,norbs
           write(iun3,'(i3,a,f9.4)')i,'  ',focc(i)
         end do

         if (iuhf.eq.1) then

           write(iun3,*)' '
           write(iun3,*)'Beta set'
           write(iun3,*)' '

           do i=1,norbs
               write(iun3,'(i3,a,f9.4)')i,'  ',focb(i)
           end do

         endif

         write(iun3,*)'----------------------------------------------'

      endif

      if (bonds .and. ipsi.ne.0) 
     &    call inferr('BONDS WITH A M.O. IS NOT PLOTTABLE',1)
      if (ipsi.ne.0) goto 2222
c
c   construct the DENSITY matrix
c
      ibonds = 0
      if (iftyp.eq.1.and.bonds) ibonds = 1

      do i=1,norbs
          do j=i,norbs

              sum = 0.d0

              do k=1,norbs
                  if (focc(k).gt.0.0d0) then
                     sum = sum + focc(k)*vectrs((k-1)*mxorb+i)
     &                                  *vectrs((k-1)*mxorb+j)
                  endif
                  if (iuhf.eq.1) then
                    if (ispd.eq.0) then
                       if (focb(k).gt.0.0d0) then
                          sum = sum + focb(k)*vectrb((k-1)*mxorb+i)
     &                                       *vectrb((k-1)*mxorb+j)
                       endif
                    else
                       if (focb(k).gt.0.0d0) then
                          sum = sum - focb(k)*vectrb((k-1)*mxorb+i)
     &                                       *vectrb((k-1)*mxorb+j)
                       endif
                    endif
                  endif
              end do

              p((i-1)*mxorb+j) = sum
              p((j-1)*mxorb+i) = sum
          end do

          if (ibonds.eq.1) p((i-1)*mxorb+i) = 
     &                     p((i-1)*mxorb+i) - averag(i)
      end do

      if ((bonds.or.ovrlap.or.atomic).and.iftyp.ne.1) then
c
c############# - substract atomic densities from molecular ####
c############# - or do overlap density                     ####
c############# - or do atomic density within molecule      ####
c
c    call new routines for external atomic density (Bjoern Pedersen)

         if (extbas) then
             call newdenmak
             if (denok) then
                ido = 1
             else
                ido = 0
             endif
             return
         endif

         if (iuhf.eq.1.and.bonds) then
             write(iun3,*)'BONDS with UHF wavefunction using '//
     +       'ROHF derived atomic densities'
         endif
         if (bonds.or.atomic)
     1   write(iun3,'(/,''ELECTRON DENSITIES DUE TO'',
     1  '' ATOMS ARE TO BE SUBTRACTED FROM DENSITY PLOT'',/)')
         if (ovrlap) write(iun3,*)'INTERATOMIC OVERLAP DENSITY USED'
         if (atomic) write(iun3,*)'ATOMIC PART LOOKED AT ONLY'
c
c IX     the start position for an atom within the molecular 
c        density matrix
c ATMD   a large array containing atomic densities.
c N      the number of basis functies per atom possibly minus 
c        polarisation functions
c ISTART the start position of an atom within ATMD 
c 
c  H  He   Li-F      Na-Cl
c  4  0    81        169
c 25  0    225       361
c-------------------------- -
c 21  0    144       192     iplus
c
c      ibasis=     (0,1,2,3 sto3g,321g,431g,631g)
c      ipolf=   0 or 1 ; are there polarisationfunctions
c
          write(iun3,*)' '
          call setbas(denok)

          if (denok) then
             ido = 1
          else
             ido = 0
          endif

          if (.not.denok) return

          do i=1,norbs
             do j=1,norbs
               q((j-1)*mxorb+i) = 0.0d0
             end do
          end do

          ix = 0

          do 700 k=1,natoms

            ibasis = nbas(k)
            ipolf  = ipol(k)

            if (nat(k).le.2) then
               iplus = 3
            elseif (nat(k).ge.3.and.nat(k).le.17) then
               iplus = 6
            else
               call inferr('BONDS: not standard basisset',1)
            endif

            n = natn(17*ibasis + nat(k))
            istart = ioffst(ibasis+1) + nbeg(17*ibasis+nat(k))
            write(iun3,*)elemnt(nat(k)),' ',bstrng(ibasis+1+2*ipolf)

            if (ovrlap) then
c
c============ do overlap densities only ========================
c
               do i=1,n+ipolf*iplus
                  do j=1,n+ipolf*iplus
                     p((ix+j-1)*mxorb+ix+i) = 0.0d0
                  end do
               end do
            else
c
c============ molecule minus atoms =============================
c
              orient = .false.
              do i=1,norien
                 if (ori(i).eq.k) then
                     orient = .true.
                     ibal   = i
                 endif
              end do
              if (orient.and.doori) then
              if (ibasis.eq.0) then
c
c------------ substract oriented sto3g o,f,s or cl ---------
c
                do i=1,3
                   do j=1,3
                      ph(i,j) = p((ix+n-3+j-1)*mxorb+(ix+n-3+i))
                   end do
                end do

                call rotatg(nat(k),idebug)

                do i=1,n-3
                 do j=1,n-3
                  p((ix+j-1)*mxorb + (ix+i)) = p((ix+j-1)*mxorb +(ix+i))
     &                  - atmd(istart+(j-1)+(i-1)*n)
                  q((ix+j-1)*mxorb + (ix+i)) = q((ix+j-1)*mxorb +(ix+i))
     &                  + atmd(istart+(j-1)+(i-1)*n)
                 end do
                end do

                do i=1,3
                 do j=1,3
                  p((ix+n-3+j-1)*mxorb + (ix+n-3+i)) = 
     &                p((ix+n-3+j-1)*mxorb + (ix+n-3+i)) - pr(i,j)
                  q((ix+n-3+j-1)*mxorb + (ix+n-3+i)) = pr(i,j)
                 end do
                end do

              else
c
c------------ substract oriented split valence o,f,s,cl -----
c
                idim=9
                if (nat(k).eq.16.or.nat(k).eq.17) idim = 13

                do i=1,idim
                 do j=1,idim
                  dmolin(j+(i-1)*idim) = p((ix+j-1)*mxorb + (ix+i))
                 end do
                end do

                call rotatg(nat(k),idebug)

                do i=1,idim
                 do j=1,idim
                  p((ix+j-1)*mxorb + (ix+i)) = p((ix+j-1)*mxorb +(ix+i))
     &                 - dofscl(j+(i-1)*idim)
                  q((ix+j-1)*mxorb + (ix+i)) = dofscl(j+(i-1)*idim)
                 end do
                end do

              endif
              else
c
c----- for all other atoms substract spherical density -----
c
                do i=1,n
                 do j=1,n
                  p((ix+j-1)*mxorb + (ix+i)) = p((ix+j-1)*mxorb +(ix+i))
     &                 - atmd(istart+(j-1)+(i-1)*n)
                  q((ix+j-1)*mxorb + (ix+i)) = 
     &                   atmd(istart+(j-1)+(i-1)*n)
                 end do
                end do

              endif
c
c============== clear out the overlap region of the density matrix
c               if ATOMIC is supplied
c
              if (atomic) then
               idavje = n + ipolf*iplus
               do i=ix+1,ix+idavje

                  if (k.ne.0) then
                     do j=1,ix
                        p((j-1)*mxorb + i) = 0.0d0
                     end do
                  endif

                  if (k.ne.natoms) then
                     do j=ix+idavje+1,norbs
                        p((j-1)*mxorb + i) = 0.0d0
                     end do
                  endif

               end do

              endif
              endif
                ix = ix + n + ipolf*iplus
700       continue

          if (bonds.and.idebug.eq.1) then
            write(iun3,*)' '
            write(iun3,*)'***** Atomic Density Matrix *****'
            write(iun3,*)' '
            call prev(q,norbs,norbs,mxorb)
            write(iun3,*)' '
          endif

      endif

c############### end bonds , overlap , atomic ##################

      if (idebug.eq.1) then
          write(iun3,'(''   DENSITY MATRIX USED BY MAP'')')
          call prev(p,norbs,norbs,mxorb)
      endif

      return

c~~~~~~~~~~~~~~~~~   PLOT PSI   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
2222  continue 
      if (ipsi.lt.0) then
         do i=1,norbs
             paa(i) = vectrb((-ipsi-1)*mxorb+i)
         end do
      else
         do i=1,norbs
             paa(i) = vectrs((ipsi-1)*mxorb+i)
         end do
      endif
      if (idebug.eq.1) then
        write(iun3,'(''  COEFFICIENTS OF M.O. TO BE PLOTTED'')')
        write(iun3,'(6f12.6)')(paa(i),i=1,norbs)
      endif

      return
      end

      logical function onden(idum)
      implicit double precision (a-h,o-z)
      logical valenc,bonds,ovrlap,atomic,doori,dolap,doelf
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap,doelf
      onden = .true.

      if (ipsi.ne.0) onden = .false.
      if (bonds.or.ovrlap.or.atomic.or.doelf) onden = .false.

      return
      end
