      subroutine addd(imo,ipsi,iao,
     &                p,paa,psi)
c THIS IS REALLY atmd
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /atopla/ xsym(numatm),ysym(numatm),zsym(numatm),
     &                isym(numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /orbhlp/ mxorb,iuhf,ispd
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      parameter (mxel=100)
      character*2 elemnt
      common /elem/elemnt(mxel)
      dimension p(*),paa(*),psi(*)

      ipreca = 0
      write(iun3,'(//''ELECTR. DENSITY/INTENSITY AT ATOMS '',
     &   ''LYING IN THE PLANE OF THE PLOT''//)')
      write(iun3,
     &   '(''  ATOM    X         Y         Z           VALUE'')')
      write(iun3,'(//)')
      do l=1,natoms
        if (isym(l).eq.1) then
           if (imo.eq.1) then
              call slater(xyz(1,l),xyz(2,l),xyz(3,l),psi)
           elseif (iao.eq.1.and.isgau.eq.0) then
              call adffun(xyz(1,l),xyz(2,l),xyz(3,l),psi)
           else
              call gaussian(xyz(1,l),xyz(2,l),xyz(3,l),psi,norbs,
     &                      ipreca,0,0)
           endif
           sum = 0.0d0
           if (ipsi.eq.0)  then
              do i=1,norbs
                 sum = sum - psi(i)*psi(i)*p((i-1)*mxorb+i)*0.5d0
                 do j=1,i
                    sum = sum + psi(i)*psi(j)*p((j-1)*mxorb+i)
                 end do
              end do
              sum = sum + sum
           else
              do i=1,norbs
                 sum = sum + paa(i)*psi(i)
              end do
           endif
          write(iun3,222)elemnt(nat(l)),(xyz(ig,l),ig=1,3),sum
        endif
      end do

222   format(4x,a2,3f10.5,2x,f13.5)

      return
      end
