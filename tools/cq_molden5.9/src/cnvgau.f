      subroutine cnvgau
      parameter (MAXITER=1000)
      implicit double precision (a-h,o-z)
      logical g92,g03
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      common /convrg/ convg1(MAXITER),convg2(MAXITER),cnvmax,cnvmin,
     &jstrt1,jend1,jstrt2,jend2,icvav1,icvav2
      character lstr*137
      character tstr*137

      call rewfil
      icvav1 = 1
      icvav2 = 1
      g92 = .false.
      g03 = .false.

      call search(tstr,'03 program',istat)
      if (istat.eq.1) then
         g03 = .true.
      else
         call rewfil
      endif

c     SCF convergence first point

      call search(tstr,'nuclear repulsion energy',istat)
      if (istat.eq.0) goto 300

      call searchd(lstr,'iter  electronic-energy',
     &                  ' Cycle ',istat)
      if (istat.eq.0) goto 300
      if (lstr(1:6).eq.' Cycle'.or.lstr(1:6).eq.' CYCLE') then
         call bckfil
         g92 = .true.
      endif

c     defer reading nuclear repulsion energy until g92 is known
      if (g92) then
         read(tstr,'(37x,f15.10)',err=300) enuc
      else
         read(tstr,'(30x,f15.10)',err=300) enuc
      endif

      jstrt1 = 0

10    call nxtlin(lstr,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 300
      ilen = ifblen(lstr)
      if (index(lstr,'SCF Done').ne.0
     &    .or.index(lstr,'SCF DONE').ne.0) goto 20
      if (g92) then
         if (lstr(1:6).ne.' Cycle'.and.lstr(1:6).ne.' CYCLE')
     &   goto 10
         read(lstr,'(7x,i3)',err=300) jend1
         call search(lstr,' E=',istat)
         read(lstr,'(3x,f22.15)',end=300,err=300) convt
      else
         if (ilen.ne.4) goto 10
         read(lstr,'(1x,i3)',err=300) jend1
         call nxtlin(lstr,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 300
         read(lstr,'(6x,f22.15)',end=300,err=300) convt
      endif
      convg1(jend1) = convt
      if (jstrt1.eq.0) jstrt1 = jend1
      goto 10

20    continue
      if (.not.g03) then
         do j=jstrt1,jend1
            convg1(j) = convg1(j) + enuc
         end do
      endif

c     SCF convergence last point

      if (g92) then
         call search(lstr,' Cycle ',istat)
      else
         call search(lstr,'iter  electronic-energy',istat)
      endif
      if (istat.eq.0) goto 400
100   call nxtlin(lstr,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 400
      goto 100
110   call bckfil
      if (g92) then
         call scback(lstr,' Cycle ',istat)
      else
         call scback(lstr,'iter  electronic-energy',istat)
      endif
      if (istat.eq.0) goto 400
      call scback(lstr,'nuclear repulsion energy',istat)
      if (istat.eq.0) goto 400
      if (g92) then
         read(lstr,'(37x,f15.10)',err=300) enuc
      else
         read(lstr,'(30x,f15.10)',err=300) enuc
      endif

      if (g92) then
         call search(lstr,' Cycle ',istat)
         call bckfil
      else
         call search(lstr,'iter  electronic-energy',istat)
      endif
      if (istat.eq.0) goto 400
      jstrt2 = 0

120   call nxtlin(lstr,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 400

      ilen = ifblen(lstr)
      if (index(lstr,'SCF Done').ne.0
     &    .or.index(lstr,'SCF DONE').ne.0) goto 200
      if (g92) then
         if (lstr(1:6).ne.' Cycle'.and.lstr(1:6).ne.' CYCLE')
     &   goto 120
         read(lstr,'(7x,i3)',err=400) jend2
         call search(lstr,' E=',istat)
         read(lstr,'(3x,f22.15)',end=400,err=400) convt
      else
         if (ilen.ne.4) goto 120
         read(lstr,'(1x,i3)',err=400) jend2
         call nxtlin(lstr,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 400
         read(lstr,'(6x,f22.15)',end=400,err=400) convt
      endif
      convg2(jend2) = convt
      if (jstrt2.eq.0) jstrt2 = jend2
      goto 120

200   continue
      if (.not.g03) then
         do j=jstrt2,jend2
            convg2(j) = convg2(j) + enuc
         end do
      endif

      return

300   icvav1 = 0
      icvav2 = 0
      return

400   icvav2 = 0
      return

      end
