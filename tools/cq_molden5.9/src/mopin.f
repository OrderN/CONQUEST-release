      subroutine mopdd(istat,ibin,impas,
     &                 vectrs,averag,p,focc,eiga,halfs,psi,nocc,ncols)

c THIS IS REALLY mopin

      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
c**************************************************************************
      logical dat3ln
      character*80 temps
      character*2 gstr
      character*137 lstr
      common /orbhlp/ mxorb,iuhf,ispd
      common /coord / xyz(3,numatm)
      common /moldat/ natoms,norbs,nelecs,nat(numatm)
      common /mopac/  nfirst(numatm),nlast(numatm),
     &                npq(numatm),pqn(numatm),emus(numatm),
     &                emup(numatm),emud(numatm),consts(numatm),
     &                constp(numatm),constd(numatm),npqref(54)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real eiga
      dimension halfs(*), psi(*), vectrs(*),averag(*),focc(*),p(*),
     &          eiga(*)
      dimension emu(numatm,3)
      equivalence (emu(1,1),emus(1)),
     &            (emu(1,2),emup(1)),
     &            (emu(1,3),emud(1))
c
c*********************** mopac file read **************************
c

      istat = 1

      do i=1,mxorb
         psi(i) = 0.d0
      end do

      natoms = 0

      if (ibin.eq.1) then

             read(iun2,end=1000,err=1000) natoms
             if (natoms.gt.numatm.or.natoms.le.0) then
                 natoms = 0
                 istat = 0
                 return
             endif
             rewind iun2
             read(iun2,end=1000,err=1000)
     &       natoms,norbs,nelecs,((xyz(i,j),j=1,natoms),i=1,3)
      else
             call nxtlin(lstr,jstat)
             if (jstat.eq.1.or.jstat.eq.2) goto 1000
             if (.not.dat3ln(lstr)) goto 1000
             call bckfil
             impas = 2
             if (linlen(lstr).gt.25) impas = 1
             if (impas.eq.1) then
                call nxtlin(lstr,jstat)
                if (jstat.eq.1.or.jstat.eq.2) goto 1000
                read(lstr,*,end=1000,err=1000)
     &          natoms,norbs,nelecs,((xyz(i,j),j=1,natoms),i=1,3)
             else
                il = linlen(lstr)
                if (il.eq.25) then
                   call nxtlin(lstr,jstat)
                   if (jstat.eq.1.or.jstat.eq.2) goto 1000
                   read(lstr,'(1x,3i8)',end=1000,err=1000)
     &                natoms,norbs,nelecs
                elseif (il.eq.0) then
                   goto 1000
                else
                   call nxtlin(lstr,jstat)
                   if (jstat.eq.1.or.jstat.eq.2) goto 1000
                   read(lstr,'(3i8)',end=1000,err=1000)
     &                natoms,norbs,nelecs
                endif
                do i=1,natoms
                   call nxtlin(lstr,jstat)
                   if (jstat.eq.1.or.jstat.eq.2) goto 1000
                   read(lstr,'(1x,3f16.8)',end=1000,err=1000)
     &                (xyz(j,i),j=1,3)
                end do
             endif
      endif

      if (norbs.gt.mxorb) then
              temps = 'Max. No. of Orbitals '//gstr(mxorb)//
     &                ' Exceeded'
              if (ibin.eq.0) call inferr(temps,1)
              istat = 0
              return
      endif

      if (ibin.eq.1) then
             read(iun2,end=1000,err=1000)
     &         (nlast(i),nfirst(i),i=1,natoms)
             read(iun2,end=1000,err=1000)
     &         ((emu(j,i),j=1,natoms),i=1,3),(nat(i),i=1,natoms)
      else
             if (impas.eq.1) then
                call nxtlin(lstr,jstat)
                if (jstat.eq.1.or.jstat.eq.2) goto 1000
                read(lstr,*,end=1000,err=1000)
     &            (nlast(i),nfirst(i),i=1,natoms)
                call nxtlin(lstr,jstat)
                if (jstat.eq.1.or.jstat.eq.2) goto 1000
                read(lstr,*,end=1000,err=1000)
     &         ((emu(j,i),j=1,natoms),i=1,3),(nat(i),i=1,natoms)
             else
                do i=1,natoms
                   call nxtlin(lstr,jstat)
                   if (jstat.eq.1.or.jstat.eq.2) goto 1000
                   read(lstr,'(1x,2i8)',err=1000)
     &               nfirst(i),nlast(i)
                end do
                do i=1,natoms
                   call nxtlin(lstr,jstat)
                   if (jstat.eq.1.or.jstat.eq.2) goto 1000
                   read(lstr,'(1x,3f16.8,i8)',err=1000)
     &            (emu(i,j),j=1,3),nat(i)
                end do
             endif
      endif


      linear = norbs*norbs
c
c  need to read in un-normalised eigenvectors and
c  inverse of overlap matrix
c
      if (ibin.eq.1) then
             read(iun2,end=1000,err=1000)
     &          ((p((i-1)*mxorb+j),j=1,norbs),i=1,norbs)
             read(iun2,end=1000,err=1000)(halfs(i),i=1,linear)
      else
             if (impas.eq.1) then
                call nxtlin(lstr,jstat)
                if (jstat.eq.1.or.jstat.eq.2) goto 1000
                read(lstr,*,end=1000,err=1000)
     &             ((p((i-1)*mxorb+j),j=1,norbs),i=1,norbs)
                call nxtlin(lstr,jstat)
                if (jstat.eq.1.or.jstat.eq.2) goto 1000
                read(lstr,*,end=1000,err=1000)(halfs(i),i=1,linear)
             else
                do i=1,norbs
                   do j=1,norbs
                       call nxtlin(lstr,jstat)
                       if (jstat.eq.1.or.jstat.eq.2) goto 1000
                       read(lstr,'(1x,f16.8)',end=1000,err=1000)
     &                    p((i-1)*mxorb+j)
                   end do
                end do
                do i=1,linear
                   call nxtlin(lstr,jstat)
                   if (jstat.eq.1.or.jstat.eq.2) goto 1000
                   read(lstr,'(1x,f16.8)',end=1000,err=1000) halfs(i)
                end do
             endif
      endif

      call averg(averag,p)
      call mulpxs(p,halfs,psi,vectrs)
      call setcst

      do i=1,norbs
         focc(i) = 0.0d0
      end do

      nelec = nelecs / 2

      do i=1,nelec
         focc(i) = 2.0d0
      end do

      nocc = nelec

      if (nelec*2.ne.nelecs) then
         focc(nelec+1) = 1.0d0
         nocc = nocc + 1
      endif

      ncols = norbs

      do i=1,ncols
         eiga(i) = 0.0e0
      end do

      return
1000  istat = 0
      call inferr('ERROR reading .gpt file',0)
c      call inferr('Unknown File Type !',1)
      return
      end

      subroutine averg(averag,p)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /orbhlp/ mxorb,iuhf,ispd
      common /moldat/ natoms,norbs,nelecs,nat(numatm)
      common /mopac/  nfirst(numatm),nlast(numatm),
     &                npq(numatm),pqn(numatm),emus(numatm),
     &                emup(numatm),emud(numatm),consts(numatm),
     &                constp(numatm),constd(numatm),npqref(54)
      dimension averag(*),p(*)

      nelec = nelecs/2
      jj = nelec + 1
      fract = 0.5d0*(nelecs-nelec*2)

      do i=1,norbs
        x = 0.d0
        do j=1,nelec
           x = x + p((j-1)*mxorb+i)**2
        end do
        x = x + p((jj-1)*mxorb+i)**2*fract
        averag(i) = x*2
      end do

c
c loop to calculate spherical-average atomic orbital occupancy
c
      do i=1,natoms

           il = nfirst(i)
           iu = nlast(i)
           ir = iu - il + 1

           goto (120,130,130,130,140,140,140,140,140),ir
140            x = 0.d0

               do j=1,5
                   ji = j + il + 3
                   x = x + averag(ji)
               end do

               x = x*0.2d0

               do j=1,5
                   ji = j + il + 3
                   averag(ji) = x
               end do

130            x = 0.d0

               do j=1,3
                   ji = j + il
                   x = x + averag(ji)
               end do

               x = x*0.333333d0

               do j=1,3
                   ji = j + il
                   averag(ji) = x
               end do

120       continue
      end do

      return
      end

      subroutine setcst
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /orbhlp/ mxorb,iuhf,ispd
      common /moldat/ natoms,norbs,nelecs,nat(numatm)
      common /mopac/  nfirst(numatm),nlast(numatm),
     &                npq(numatm),pqn(numatm),emus(numatm),
     &                emup(numatm),emud(numatm),consts(numatm),
     &                constp(numatm),constd(numatm),npqref(54)
      dimension fa(20)
      data npqref/1,0, 2,2,2,2,2,2,2,0, 
     & 2,2,2,2,2,2,2,0,
     & 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,0, 
     & 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5/

      do i=1,numatm
          emus(i) = emus(i)/0.529167d0
          emup(i) = emup(i)/0.529167d0
          emud(i) = emud(i)/0.529167d0
      end do

      do i=1,natoms
         npq(i) = npqref(nat(i))
      end do

      fa(2) = dsqrt(2.d0)

      do i=4,20,2
         fa(i) = fa(i-2)*dsqrt(i*i-i*1.d0)
      end do

c
c setting up of orbital constants for density map
c
      do i=1,natoms

          j = 2*npq(i)
          pqn(i) = npq(i)
          if (j.gt.6) pqn(i) = pqn(i) - 1.d0
          k = pqn(i)*2.01
c
c  0.282095  =  sqrt( 1/(4*pi))
c  0.48860   =  sqrt( 3/(4*pi))
c  1.092547  =  sqrt(15/(4*pi))
c
c   for slater-type orbitals radial normalisation constant 
c
c    2**(pqn+1/2)*sqrt(factorial((2*pqn)))
c
c   angular normalisation constants
c
c   s:   1,   p:  sqrt(3),   d: sqrt(15)
c
          if (fa(j).eq..0d0) then
             goto 1000
          else
             consts(i)=(2.d0*emus(i))**(npq(i)+0.5d0)/fa(j)*0.282095
             constp(i)=(2.d0*emup(i))**(npq(i)+0.5d0)/fa(j)*0.48860
          endif

          if (fa(k).eq..0d0) then
             goto 1000
          else
             constd(i)=(2.d0*emud(i))**(pqn(i)+0.5d0)/fa(k)*1.092547
          endif

      end do

      do i=1,natoms
          pqn(i) = pqn(i) - 1
          npq(i) = npq(i) - 1
      end do

      return

1000  print*,'error fa'
      return
      end

      subroutine mulpxs(p,halfs,psi,vectrs)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /orbhlp/ mxorb,iuhf,ispd
      common /moldat/ natoms,norbs,nelecs,nat(numatm)
      common /mopac/  nfirst(numatm),nlast(numatm),
     &                npq(numatm),pqn(numatm),emus(numatm),
     &                emup(numatm),emud(numatm),consts(numatm),
     &                constp(numatm),constd(numatm),npqref(54)
      dimension halfs(*),psi(*),vectrs(*),p(*)
c
c matrix multiply eigenvectors by inverse square-root of overlap matrix
c
      do i=1,norbs

           l = 0
           do j=1,norbs
               x = 0.d0
               do k=1,norbs
                   l = l + 1
                   x = x + p((i-1)*mxorb+k)*halfs(l)
               end do
               psi(j) = x
           end do

           do j=1,norbs
               p((i-1)*mxorb+j) = psi(j)
           end do

      end do

      do i=1,norbs
           do j=1,norbs
               vectrs((i-1)*mxorb+j) = p((i-1)*mxorb+j)
           end do
      end do

      return
      end

      subroutine mopxyd(iun,coo,ianz,iconn)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (mxcon=10)
      parameter (maxsym=108)
      parameter (mxel=100)
      character*2 elemnt,tocapf
      character*3 atom,atomt
      common /athlp/  iatoms, mxnat
      common /elem/   elemnt(mxel)
      character*137 line
      common /curlin/ line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension coo(3,*),ianz(*),iconn(mxcon+1,*)

      toang = 0.52917706d0
      natoms = iatoms
      nbnds = 0
      idochg = 0

      write(iun,'(a)') 'PM6 XYZ MMOK '
      write(iun,'(a)') '  '
      write(iun,'(a)') '  '

      do i=1,natoms
         atomt = '   '
         atomt(1:2) = tocapf(elemnt(ianz(i)))
         call leftj(atomt,atom)
         write(iun,'(a3,a1,f10.4,a3,f10.4,a3,f10.4,a3)') 
     &      atom,' ',coo(1,i)*toang,' 1 ',coo(2,i)*toang,' 1 ',
     &               coo(3,i)*toang,' 1 '
      end do

      return
      end
