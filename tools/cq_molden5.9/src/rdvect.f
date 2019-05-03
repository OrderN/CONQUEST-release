      subroutine rdvecd(idebug,ig94,istats,
     &                  vectrs,vectrb,focc,eiga,eigb,ncols,ncolb)
      implicit double precision (a-h,o-z)
      parameter (mxcol=10)
      parameter (numatm=2000)
      character*137 line
      common /orbhlp/ mxorb,iuhf,ispd
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      logical natorb,nlmo,reo
      real eiga,eigb
      dimension vectrs(*),vectrb(*),focc(*),eiga(*),eigb(*)
      dimension vt(15,mxcol),irord(15)

      istats = 1
      reo = .false.

      irord(1)  = 15
      irord(2)  = 5
      irord(3)  = 1
      irord(4)  = 14
      irord(5)  = 13
      irord(6)  = 9
      irord(7)  = 4
      irord(8)  = 6
      irord(9)  = 2
      irord(10) = 12
      irord(11) = 10
      irord(12) = 3
      irord(13) = 11
      irord(14) = 8
      irord(15) = 7

      natorb = .false.
      nlmo = .false.
      npl = 5
      if (idebug.eq.1) 
     &      call inferr('Looking for Molecular Orbital Coefficients',0)

      call searcht(line,'Molecular Orbital Coefficients',
     &    'Natural Orbital Coefficients',
     &    'NLMOs in the AO basis:',istat)
      if (istat.eq.0) then
       call inferr('no MO COEFFICIENTS found on GAUSSIAN outputfile',0)
        goto 150
      endif
      if (index(line,'Natural').ne.0.or.index(line,'NATURAL').ne.0) then
          natorb = .true.
          iuhf = 0
      endif
      if (index(line,'NLMO').ne.0) then
          nlmo = .true.
          npl = 8
          iuhf = 0
          ido5d = 1
          ido7f = 1
          ido9g = 1
      endif
      imin = 0
      imax = 0
      npass = (ncols - 1) / npl + 1
      do ipass=1,npass
         imin = imax + 1
         imax = imax + npl
         imax = min0(ncols,imax)
         do i =1,10
            call nxtlin(line,jstat)
            call tocap(line,80)
            if (nlmo) then
               if (index(line,' AO ').ne.0) then
                  call redel(line,1)
                  goto 100
               endif
            else
               if (index(line,'EIGENVALUES').ne.0) goto 100
            endif
         end do
         if (imax.lt.ncols) goto 150
         return
100      continue
         if (natorb) then
            read(line,'(21x,5f10.5)',err=150) (focc(j),j=imin,imax)
            do j=imin,imax
                eiga(j) = 0.0d0
            end do
         else
            if (nlmo) then
               do j=imin,imax
                   eiga(j) = 0.0d0
               end do
            else
             read(line,'(21x,5f10.5)',err=150) (eiga(j),j=imin,imax)
            endif
         endif

         istrt = 0
         itel = 0
         nc = imax - imin + 1

         do irow=1,norbs
           call nxtlin(line,jstat)
           if (jstat.eq.2) goto 150
cnull           call rmnull(line)
           if (nlmo) then
              read(line,'(16x,8f8.4)',err=150)
     &         (vectrs((j-1)*mxorb+irow),j=imin,imax)
              if (index(line,'D6').ne.0.or.index(line,'d6').ne.0)
     &         ido5d = 0
              if (index(line,'F8').ne.0.or.index(line,'f8').ne.0)
     &         ido7f = 0
              if (index(line,'G10').ne.0.or.index(line,'g10').ne.0)
     &         ido9g = 0
           else
              read(line,'(21x,5f10.5)',err=150)
     &         (vectrs((j-1)*mxorb+irow),j=imin,imax)
              if (index(line,'D-2').ne.0.or.index(line,'d-2').ne.0)
     &         ido5d = 1
              if (index(line,'F-2').ne.0.or.index(line,'f-2').ne.0)
     &         ido7f = 1
              if (index(line,'G-2').ne.0.or.index(line,'g-2').ne.0)
     &         ido9g = 1
           endif

           if (istrt.eq.0.and.index(line,'ZZZZ').ne.0) then
               istrt = irow
               itel = 0
               reo = .true.
           endif

           if (istrt.ne.0) then
              itel = itel + 1
              do k=1,nc
                 vt(itel,k) = vectrs((imin+k-2)*mxorb+irow)
              end do
           endif

           if (itel.eq.15) then
              do l=1,15
                 do k=1,nc
                    vectrs((imin+k-2)*mxorb + (istrt+l-1)) = 
     &                   vt(irord(l),k)
                 end do
              end do
              istrt = 0
              itel = 0
           endif

         end do
      end do

      if (iuhf.eq.1) then
         call nxtlin(line,jstat)
         if (jstat.eq.2) goto 150
         imin = 0
         imax = 0
         do ipass=1,npass
            imin = imax + 1
            imax = imax + npl
            imax = min0(ncols,imax)
            do i =1,10
               call nxtlin(line,jstat)
               if (jstat.eq.2) goto 150
               call tocap(line,80)
               if (index(line,'EIGENVALUES').ne.0) goto 200
            end do
            if (imax.lt.ncolb) goto 150
            return
200         continue

            read(line,'(21x,5f10.5)',err=150) (eigb(j),j=imin,imax)

            istrt = 0
            itel = 0
            nc = imax - imin + 1

            do irow=1,norbs
cnull              call nxtlin(line,jstat)
cnull              if (jstat.eq.2) goto 150
cnull              call rmnull(line)
cnull              read(line,'(21x,5f10.5)',err=150)
              call nxtlin(line,jstat)
              read(line,'(21x,5f10.5)',err=150)
     &             (vectrb((j-1)*mxorb+irow),j=imin,imax)


              if (istrt.eq.0.and.index(line,'ZZZZ').ne.0) then
                  istrt = irow
                  itel = 0
                  reo = .true.
              endif

              if (istrt.ne.0) then
                 itel = itel + 1
                 do k=1,nc
                    vt(itel,k) = vectrb((imin+k-2)*mxorb+irow)
                 end do
              endif

              if (itel.eq.15) then
                 do l=1,15
                    do k=1,nc
                       vectrb((imin+k-2)*mxorb + (istrt+l-1)) = 
     &                      vt(irord(l),k)
                    end do
                 end do
                 istrt = 0
                 itel = 0
              endif

            end do
         end do
      endif

      if (reo) then
          print*,' '
          print*,'========================================'//
     &           '======================================'
          print*,'Changed order of G functions:'
          print*,' '
          print*,'zzzz,yzzz,yyzz,yyyz,yyyy,xzzz,xyzz,xyyz,',
     &           'xyyy,xxzz,xxyz,xxyy,xxxz,xxxy,xxxx ->'
          print*,' '
          print*,'xxxx,yyyy,zzzz,xxxy,xxxz,yyyx,yyyz,zzzx,',
     &           'zzzy,xxyy,xxzz,yyzz,xxyz,yyxz,zzxy'
          print*,' '
          print*,'========================================'//
     &           '======================================'
          print*,' '
      endif

      return

150   call inferr('error in reading vectors',1)
      if (ig94.eq.1) then
         call inferr('use IOP(6/7=3) in gaussian inputfile',1)
      else
         call inferr('use IOP(6/7=1) in gaussian inputfile',1)
      endif
      istats = 0
      return
      end

      subroutine rdvcd(idebug,istats,
     &           vectrs,vectrb,focc,focb,eiga,eigb,ncols,ncolb)
      implicit double precision (a-h,o-z)
      parameter (mxcol=10)
      parameter (numatm=2000)
      character*137 line,str
      common /curlin/ line
      common /orbhlp/ mxorb,iuhf,ispd
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      logical natorb,nlmo,reo
      real eiga,eigb
      integer getlin
      dimension vectrs(*),vectrb(*),focc(*),focb(*),eiga(*),eigb(*)
      dimension v(mxorb)

      istats = 1
      ist = 0

      call search(line,'%molecular orbital range',ist)
      if (getlin(0).ne.1) goto 100

      ktype = nxtwrd(str,nstr,itype,rtype)

      if (ktype.eq.2) then
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            norbs = itype
            ncols = norbs
c            print*,"norbs ",norbs
         else
            goto 100
         endif
      else
         goto 100
      endif
         
      call rewfil

      call search(line,'%molecular orbital energies',ist)

      irow = 1

      do while (.true.)
10       if (getlin(0).ne.1) goto 100
         do i=1,4
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.3) then
               eiga(irow) = rtype
               irow = irow + 1
               if (irow.gt.norbs) goto 20
            else
               goto 10
            endif
         end do
      end do

20    continue

      call search(line,'%molecular orbital occupations',ist)

      irow = 1

      do while (.true.)
30       if (getlin(0).ne.1) goto 100
         do i=1,4
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.3) then
               focc(irow) = rtype
               irow = irow + 1
               if (irow.gt.norbs) goto 40
            else
               goto 30
            endif
         end do
      end do

40    continue

      call search(line,'%molecular orbital vectors',ist)

      irow = 1
      icol = 1

      do while (.true.)
         if (getlin(0).ne.1) goto 100
         do while (.true.)
            ktype = nxtwrz(str,nstr,itype,rtype)

            if (ktype.eq.3.or.ktype.eq.4) then
               if (ktype.eq.3) then
                  v(irow) = rtype
                  irow = irow + 1
               else
                  do j=0,itype-1
                       v(irow+j) = 0.0d0
                  end do
                  irow = irow + itype
               endif
            else if (ktype.eq.0) then
               if (irow-1.eq.norbs) then
                  do j=1,irow-1
c                     print*,"v(",j,")=",v(j)
                     vectrs((icol-1)*mxorb+j) = v(j)
                  end do
                  irow = 1
                  icol = icol + 1
                  if (icol.gt.norbs) goto 51
               endif
               goto 50

            else
               goto 51
            endif

         end do
50       continue
      end do

51    idbl = nelecs/2

c      if (nelecs-idbl*2.eq.0) return

c search for beta information

      call rewfil

60    call search(line,'%molecular orbital range',ist)
      if (ist.eq.0) goto 100
      if (icdex(line,'beta').eq.0) goto 60
      iuhf = 1

      if (getlin(0).ne.1) goto 100

      ktype = nxtwrd(str,nstr,itype,rtype)

      if (ktype.eq.2) then
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            ncolb = itype
c            print*,"norbs ",norbs
         else
            goto 100
         endif
      else
         goto 100
      endif
         
      call rewfil

70    call search(line,'%molecular orbital energies',ist)
      if (ist.eq.0) goto 100
      if (icdex(line,'beta').eq.0) goto 70

      irow = 1

      do while (.true.)
80       if (getlin(0).ne.1) goto 100
         do i=1,4
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.3) then
               eigb(irow) = rtype
               irow = irow + 1
               if (irow.gt.norbs) goto 90
            else
               goto 80
            endif
         end do
      end do

90    continue

      call search(line,'%molecular orbital occupations',ist)

      irow = 1

      do while (.true.)
81       if (getlin(0).ne.1) goto 100
         do i=1,4
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.3) then
               focb(irow) = rtype
               irow = irow + 1
               if (irow.gt.norbs) goto 82
            else
               goto 81
            endif
         end do
      end do

82    continue

      call search(line,'%molecular orbital vectors',ist)

      irow = 1
      icol = 1

      do while (.true.)
         if (getlin(0).ne.1) goto 100
         do while (.true.)
            ktype = nxtwrz(str,nstr,itype,rtype)

            if (ktype.eq.3.or.ktype.eq.4) then
               if (ktype.eq.3) then
                  v(irow) = rtype
                  irow = irow + 1
               else
                  do j=0,itype-1
                       v(irow+j) = 0.0d0
                  end do
                  irow = irow + itype
               endif
            else if (ktype.eq.0) then
               if (irow-1.eq.norbs) then
                  do j=1,irow-1
c                     print*,"v(",j,")=",v(j)
                     vectrb((icol-1)*mxorb+j) = v(j)
                  end do
                  irow = 1
                  icol = icol + 1
                  if (icol.gt.norbs) return
               endif
               goto 99

            else
               goto 100
            endif

         end do
99       continue
      end do

      return

100   istats = 0
      return
      end
