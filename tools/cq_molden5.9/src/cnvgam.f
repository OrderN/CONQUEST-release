      subroutine cnvgam
      parameter (MAXITER=1000)
      implicit double precision (a-h,o-z)
      logical direct,dft,casscf
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      common /convrg/ convg1(MAXITER),convg2(MAXITER),cnvmax,cnvmin,
     &jstrt1,jend1,jstrt2,jend2,icvav1,icvav2
      common /orbhlp/ mxorb,iuhf,ispd
      character lstr*137

cd     write(iun3,'(a)')'enter subroutine cnvgam'
      call rewfil
      icvav1 = 1
      icvav2 = 1
      jstrt1 = 0
      jstrt2 = 0
      jend1 = 0
      jend2 = 0
      iels  = 2
      iels2 = 10
      if (iuhf.eq.1) then
         iels  = 3
         iels2 = 12
      endif

      direct = .false.
      dft    = .false.
      casscf = .false.
      call search(lstr,'casscf options',istat)
      if (istat.eq.1) casscf = .true.
      call search(lstr,'* Direct-scf',istat)
      if (istat.eq.1) direct = .true.
      call search(lstr,'CCP1 DFT MODULE',istat)
      if (istat.eq.1) dft = .true.
      if (casscf) goto 60

c     SCF convergence first point

      call searchd(lstr,'  cycle          total     electronic',
     &                  '  cycle           total      electronic',istat)
      if (istat.eq.0) goto 300
      call redel(lstr,iels)
      jstrt1 = 0

10    if (direct) then
         call searchd
     &        (lstr,'  cycle          total     electronic',
     &              '  cycle           total      electronic',istat)
         if (istat.eq.0) goto 20
         call redel(lstr,iels)
      endif
      if (.not.direct.and.dft) then
         call redel(lstr,iels2)
      endif
      call nxtlin(lstr,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 300
      if (lstr(1:8).eq.'        ') goto 20
      read(lstr,'(i4,4x,f15.8)',err=300) jend1,convt
      convg1(jend1) = convt
      if (jstrt1.eq.0) jstrt1 = jend1
      goto 10

60    continue
      jstrt1 = 1
      jend1  = 0
70    call search(lstr,'cycle    time   tester               energy',
     &            istat)
      if (istat.ne.0) then
         call redel(lstr,2)
80       call nxtlin(lstr,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 300
         if (lstr(1:8).ne.'        ') then
            jend1=jend1+1
            read(lstr,'(24x,f20.7)',err=300)convg1(jend1)
            goto 80
         endif
         goto 70
      endif

20    continue

c     Possibly SCF convergence info when incomplete scf

      call search(lstr,'terminating due to incomplete scf',istat)
      if (istat.eq.0) goto 400

50    call bckfil
      call bckfil
      call nxtlin(lstr,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 400
      if (lstr(1:8).eq.'        ') goto 50
      call bckfil
      jend2 = 0

100   call nxtlin(lstr,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 400
      if (lstr(1:8).eq.'        ') goto 200
      read(lstr,'(i4,4x,f15.8)',err=400) jstrt2,convt
      convg2(jstrt2) = convt
      if (jend2.eq.0) jend2 = jstrt2
      call bckfil
      call bckfil
      goto 100
c
200   continue
cd     write(iun3,'(a)')'leave subroutine cnvgam'
      return

300   icvav1 = 0
      icvav2 = 0
cd     write(iun3,'(a)')'leave subroutine cnvgam'
      return

400   icvav2 = 0
cd     write(iun3,'(a)')'leave subroutine cnvgam'
      return

      end
