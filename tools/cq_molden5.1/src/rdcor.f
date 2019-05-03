      subroutine rdcor(idebug,istat)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      character*137 line,str
      integer getlin
      logical gnreal
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /gauori/ nzm,nso,nio,nzo,ioropt,ifor,
     &                ixyz98,iopr,isymm,irc,imp2,icntp,itd
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /curlin/ line
      dimension r(3)

      ig98 = 0

      if (isymm.eq.0) goto 1000

      if (idebug.eq.1)
     &   call inferr('looking for Standard orientation',0)

      call search(line,'Standard orientation:',istat)
      if (istat.eq.0) then
          if (idebug.eq.1)
     &      call inferr('Standard orientation not found!',0)
          return
      else
          call readel(line,2)
          read(iun2,'(a)',err=20) line
          if (icdex(line,'Type').ne.0) ig98 = 1
          call readel(line,1)
      endif

      natoms = 0

      do while ( .true. )
        if (getlin(1).eq.1) then
         if (line(2:5).eq.'----') goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.2) goto 20
         natoms = itype
         if ( natoms .gt. numatm ) then
             call inferr('Exceeding Max Atoms!',0)
             goto 20
         endif
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.2) goto 20
         nat(natoms) = itype
         if (ig98.eq.1) then
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.2) goto 20
         endif
         if (gnreal(r,3,.false.)) then
             do j=1,3
                xyz(j,natoms) = r(j)
             end do
         else 
             goto 20
         endif
        endif
      end do


c
c===== Z-Matrix/ Input orientation
c
1000  if (idebug.eq.1) call inferr(
     &    'looking for Z-Matrix/Input orientation',0)

      call searchd(line,'Z-Matrix orientation:',
     &             'Input orientation:',istat)
      if (istat.eq.0) then
          if (idebug.eq.1) call inferr(
     &    'Z-Matrix/Input orientation not found!',0)
          return
      else
          call readel(line,2)
          read(iun2,'(a)',err=20) line
          if (icdex(line,'Type').ne.0) ig98 = 1
          call readel(line,1)
      endif

      natoms = 0
      do while ( .true. )
         if (getlin(1).eq.1) then
            if (line(2:4).eq.'---') goto 100
            natoms = natoms + 1
            ktype = nxtwrd(str,nstr,itype,rtype)
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.2) goto 20
            nat(natoms) = itype
            if (ig98.eq.1) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.ne.2) goto 20
            endif
            if (gnreal(r,3,.false.)) then
                do j=1,3
                   xyz(j,natoms) = r(j)
                end do
            else 
                goto 20
            endif
         endif
      end do

100   continue
      if (idebug.eq.1) then
         do i=1,natoms
            write(iun3,'(i3,3f12.6)') nat(i),(xyz(j,i),j=1,3)
         end do
         write(iun3,*) ' '
      endif
c
c get rid off dummy atoms
c
200   continue
      do i=1,natoms
         if (nat(i).lt.0) then
            do j=i+1,natoms
               nat(j-1) = nat(j)
               do k=1,3
                  xyz(k,j-1) = xyz(k,j)
               end do
            end do
            natoms = natoms - 1
            goto 200
         endif
      end do

      toang = 0.52917706d0
c
c convert to atomic units
c
      do i=1,natoms
         do j=1,3
            xyz(j,i) = xyz(j,i) / toang
         end do
      end do

      istat = 1
      return

20    call inferr('Error reading Standard orientation!',1)
      istat = 0
      return
      end
