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
            read(iun2,'(a)') line
            call tocap(line,80)
            if (nlmo) then
               if (index(line,' AO ').ne.0) then
                  call readel(line,1)
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
           read(iun2,'(a)',err=150) line
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
         read(iun2,'(a)',err=150) line
         imin = 0
         imax = 0
         do ipass=1,npass
            imin = imax + 1
            imax = imax + npl
            imax = min0(ncols,imax)
            do i =1,10
               read(iun2,'(a)',err=150) line
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
cnull              read(iun2,'(a)',err=150) line
cnull              call rmnull(line)
cnull              read(line,'(21x,5f10.5)',err=150)
              read(iun2,'(21x,5f10.5)',err=150)
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
