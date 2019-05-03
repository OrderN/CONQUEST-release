      subroutine gaupod(istat,coo,ianz)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      parameter (numat1=20000)
      parameter (maxfat=1000)
      parameter (numatm=2000)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat) 
      common /athlp/ iatoms, mxnat
      common /gauori/ nzm,nso,nio,nzo,ioropt,ifor,
     &                ixyz98,iopr,isymm,irc,imp2,icntp,itd
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character lstr*137
      common /curlin/ lstr
      character*23 stmp
      character*25 sirc
      integer getlin
      logical gnreal
      dimension itemp(numat1),ctemp(3,numat1)
      dimension coo(3,*),ianz(*),r(3)

      toang = 0.52917706d0
      ig98 = 0

      stmp = 'Input orientation:'
      ns = 18
      if (ixyz98.gt.0.and.(ioropt.eq.2)) then
          stmp = 'Standard orientation:'
          ns = 22
      endif

      sirc = 'Item               Value'
      nis = 24
      if (irc.eq.2) then
         sirc = '  Delta-x Convergence'
         nis = 21
      endif

c     istat = 0     no Z-Matrix orientation found
c     istat = -1    no Cartesian Forces found
c     istat = ge 1  both Standard orientation and forces found

c
c     Gaussian forces are in the Z-matrix orientation (not in standard)
c
5     if ((nzm.eq.nzo.or.2*nzm.eq.nzo).and.nzm.gt.0.and.ixyz98.ne.2) 
     &   then
         call search(lstr,'Z-MATRIX (ANGSTROMS',istat)
      elseif (nzo.gt.0.and.ixyz98.ne.2) then
         call search(lstr,'Z-Matrix orientation:',istat)
      else
         call search(lstr,stmp(1:ns),istat)
      endif

      if (istat.eq.0 ) return
      if (index(lstr,'Z-MATRIX (ANGSTROMS').eq.0) then
         call haszm(.false.)
         backspace iun2
      else
         call readel(lstr,2)
         call convzmat(coo,ianz,iatoms,1,0,1)
      endif

      call searchd(lstr,'Z-Matrix orientation:',stmp(1:ns),istat)
      if (istat.eq.0) then
         return
      else
         call readel(lstr,2)
         read(iun2,'(a)',err=100) lstr
         if (icdex(lstr,'Type').ne.0) ig98 = 1
         call readel(lstr,1)
      endif

      nz = 0

      do while ( .true. )
         if (getlin(1).eq.1) then
            if (lstr(2:4).eq.'---') goto 10
            nz = nz + 1
            ktype = nxtwrd(lstr,nstr,itype,rtype)
            ktype = nxtwrd(lstr,nstr,itype,rtype)
            if (ktype.ne.2) goto 100
            itemp(nz) = itype
            if (ig98.eq.1) then
               ktype = nxtwrd(lstr,nstr,itype,rtype)
               if (ktype.ne.2) goto 100
            endif
            if (gnreal(r,3,.false.)) then
                do j=1,3
                   ctemp(j,nz) = r(j)
                end do
            else
                goto 100
            endif
         endif
      end do

10    continue

c get rid off dummy atoms

      iatoms = 0
      do 400 i=1,nz
         if(itemp(i))400,500,500
500      iatoms = iatoms + 1
         ianz(iatoms) = itemp(i)
         do j=1,3
            coo(j,iatoms) = ctemp(j,i)
         end do
400   continue 

      do i=1,iatoms
c         print*,'coo ',i,(coo(j,i),j=1,3)
         do j=1,3
             coo(j,i) = coo(j,i)/toang
         end do
      end do

c gaussian 98 doesnt do more than 50 atoms forces print

c      if (ixyz98.gt.0.and.iatoms.gt.50.and.iopr.eq.0) goto 100

      if (((nzm.eq.nzo.or.2*nzm.eq.nzo).and.nzm.gt.0)
     &   .or.nzo.gt.0) then
         call searcht(lstr,'Forces (Hartrees/Bohr)',
     &   'Z-MATRIX (ANGSTROMS',
     &   'Z-Matrix orientation:',istat)
      else
         call searchd(lstr,'Forces (Hartrees/Bohr)',
     &   stmp(1:ns),istat)
      endif

      if (istat.eq.0) goto 200
      if (icdex(lstr,'Forces (Hartrees/Bohr)').eq.0) then
           istat = -1
           backspace iun2
           return
      endif
      call readel(lstr,2)
      do i=1,iatoms
         read(iun2,'(23x,3(f15.9))',err=200)(fxyz(j,i),j=1,3)
c         print*,(fxyz(j,i),j=1,3)
      end do

      if (irc.ne.0) then

         if (((nzm.eq.nzo.or.2*nzm.eq.nzo).and.nzm.gt.0)
     &      .or.nzo.gt.0) then
            call searcht(lstr,sirc(1:nis),
     &      'Z-MATRIX (ANGSTROMS',
     &      'Z-Matrix orientation:',istat)
         else
            call searchd(lstr,sirc(1:ns),stmp(1:ns),istat)
         endif
         if (istat.eq.0) goto 300

         if (irc.eq.1) then

            if (icdex(lstr,'Item').eq.0) then
                 backspace iun2
                 goto 300
            endif

            icv1 = 0
            icv2 = 0
            icv3 = 0
            icv4 = 0

            call readel(lstr,1)
            if (icdex(lstr,'YES').ne.0) icv1 = 1
            call readel(lstr,1)
            if (icdex(lstr,'YES').ne.0) icv2 = 1
            call readel(lstr,1)
            if (icdex(lstr,'YES').ne.0) icv3 = 1
            call readel(lstr,1)
            if (icdex(lstr,'YES').ne.0) icv4 = 1
            icv = icv1 + icv2 + icv3 + icv4
            if (icv.ne.4) goto 5
            istat = 1

         elseif (irc.eq.2) then

            if (icdex(lstr,' NOT met').ne.0) goto 5
            istat = 1

         endif

      endif

      return

100   istat = 0
      return
200   istat =-1 
      return
300   istat = 1
      return
      end
