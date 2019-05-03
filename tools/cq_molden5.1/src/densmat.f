      subroutine densmad(idebug,
     &                   vectrs,vectrb,focc,focb,p,nocc)

c THIS IS REALLY densmat

c NEW style atomic density matrices to file (Bjoern Pedersen)

      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs, nat(numatm)
      common /orbhlp/ mxorb,iuhf,ispd
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*137 line
      logical gamess,gamus,gauss,statio
      dimension vectrs(*),vectrb(*),focc(*),focb(*),p(*)

c--------------read vectors off atom and create
c--------------atomic density matrix in data form
c============read norbs nelecs ==================

c---------is it gamess, gaussian or mopac ?-----------
cd     write(iun3,'(a)')'enter subroutine densmat'
      gamess = .false.
      gauss = .false.
      gamus = .false.

      rewind iun2

      call searchd(line,'g a m e s s','M.W.SCHMIDT',istat)
      if (istat.ne.0) then
         gamess = .true.
         if (index(line,'M.W.SCHMIDT').ne.0) gamus = .true.
      else
         rewind iun2
         call seardu(line,'Gaussian System',
     &        'part of the Gaussian',istat)
         if (istat.ne.0) then
            gauss = .true.
         endif
      endif


c     try gamess-US, then gaussian then gamess-UK(old routine)    

      rewind iun2
      istat=0
      if (gamus) then 
         call rdgamu(idebug,.false.,statio,irtype,.false.,istat)
      else if (gauss) then
         call rdgaus(idebug,.false.,statio,irtype,istat)
      else if (gamess)then
         call search(line,'total number of basis func',istat)
         if (index(line,'cartesian').ne.0) then
            read(line,'(45x,i5)') norbs
            read(iun2,'(a)') line
            read(iun2,'(a)') line
            read(line,'(45x,i5)') nelecs
         else
            read(line,'(38x,i5)') norbs
            read(iun2,'(a)') line
            read(line,'(38x,i5)') nelecs
         endif
      
c
c     read occupancies
c
         call search(line,'m.o.  irrep  orbital',istat)
         call readel(line,2)
         nocc = 0
         do i=1,norbs
            read(iun2,'(a)') line
            read(line,'(11x,f16.8,f20.7)') eig,focc(i)
            if (focc(i).gt.0.0d0) nocc = nocc + 1
         end do
c
c     read vectors
c
         if (idebug.eq.1) write(iun3,*)'read in occupancies'
         call searchd(line,'gvb natural orbital','eigenvectors'
     &                ,istat)
         if ( index(line,'gvb natural orbital').ne.0)then
           call readel(line,11)
         else
           call readel(line,9)
         endif
         call readvv(vectrs,norbs,nocc,.false.)

c     end of old input(gamess-UK)

      endif


      if (idebug.eq.1) then
          write(iun3,*)' norbs=',norbs,' nelecs=',nelecs
          write(iun3,*)' vectors'
          call prev(vectrs,norbs,norbs,mxorb)
      endif

c
c     construct density matrix
c
      call dmat(p,vectrs,vectrb,focc,focb)

      if (idebug.eq.1) then
          write(iun3,*)'p-matrix'
          call prev(p,norbs,norbs,mxorb)
      endif

      open(unit=99,file='basiinf.mbi',form='formatted',status='unknown')

      write (99,'(A4,A3)') 'ATOM',' XX'
      write (99,'(A8,A8)') 'CHARGE','NORBS'
      write (99,'(i4,i10)') nelecs,norbs
      do i=1,norbs
         write (99,'(6f15.10)') (p((j-1)*mxorb + i),j=1,norbs) 
      end do
      write (99,*) '*******************************************' 
      close(99)

cd     write(iun3,'(a)')'leave subroutine densmat'
      return
      end

c OLD style gamess-UK atomic density matrices, as data statements

      subroutine densmtd(idebug,
     &                   v,pop,p,nocc)

c THIS ISREALLY densmto

      implicit double precision (a-h,o-z), integer ( i-n)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*137 line
      common /orbhlp/ mxorb,iuhf,ispd
      dimension v(*),pop(*),p(*)

c-------------- read vectors off atom and create
c-------------- atomic density matrix in data form
c============ read norbs nelecs ==================

cd     write(iun3,'(a)')'enter subroutine densmto'
c      call search(line,'total number of basis func',istat)
      call search(line,'basis functions',istat)
      if (index(line,'cartesian').ne.0) then
         read(line,'(45x,i5)') norbs
         read(iun2,'(a)') line
         read(iun2,'(a)') line
         read(line,'(45x,i5)') nelecs
      else
         read(line,'(38x,i5)') norbs
         read(iun2,'(a)') line
         read(line,'(38x,i5)') nelecs
      endif

      if (idebug.eq.1) 
     &    write(iun3,*)' norbs=',norbs,' nelecs=',nelecs
c
c     read occupancies
c
      call search(line,'m.o.  irrep  orbital',istat)
      call readel(line,2)
      nocc = 0
      do i=1,norbs
         read(iun2,'(a)') line
         read(line,'(11x,f16.8,f20.7)') eig,pop(i)
         if (pop(i).gt.0.0d0) nocc = nocc + 1
      end do
c
c     read vectors
c
      if (idebug.eq.1) write(iun3,*)'read in occupancies'
      call searchd(line,'gvb natural orbital','eigenvectors'
     &             ,istat)
      if ( index(line,'gvb natural orbital').ne.0)then
        call readel(line,11)
      else
        call readel(line,9)
      endif
      call readvv(v,norbs,nocc,.false.)
      if (idebug.eq.1) write(iun3,*)'vectors'
c
c     construct density matrix
c
      if (idebug.eq.1) call prev(v,norbs,norbs,mxorb)
      call dmato(v,pop,norbs,p)
      if (idebug.eq.1) write(iun3,*)'p-matrix'
      if (idebug.eq.1) call prev(p,norbs,norbs,norbs)
c
c     print densitymatrix in data form
c
      n2    = norbs * norbs
      nrow  = n2 / 3
      nlast = n2 - nrow*3
      noff  = 0
      if (nlast.eq.0) noff = -1
      ntel  = 1
      write(iun3,'(a)')'      data string/'

      do k=1,nrow+noff
         write(iun3,'(a6,3(f15.10,a2,a1))')'     &',p(ntel),'d0',
     &   char(44),p(ntel+1),'d0',char(44),p(ntel+2),'d0',char(44)
         ntel = ntel + 3
      end do

      if (nlast.eq.0) then
         write(iun3,'(a6,3(f15.10,a2,a1))')'     &',p(ntel),'d0',
     &   char(44),p(ntel+1),'d0',char(44),p(ntel+2),'d0','/'
      endif

      if (nlast.eq.1) then
         write(iun3,'(a6,f15.10,a2,a1)')'     &',p(ntel),'d0','/'
      elseif (nlast.eq.2) then
         write(iun3,'(a6,2(f15.10,a2,a1))')'     &',p(ntel),'d0',
     &   char(44),p(ntel+1),'d0','/'
      endif

cd     write(iun3,'(a)')'leave subroutine densmto'
      return
      end
