      subroutine rdinfd(npts1,npts2,isubtr,istat,dens,denst)
      implicit double precision (a-h,o-z)

c THIS IS REALLY rdinfo

      parameter (numatm=2000)
      logical valenc,bonds,ovrlap,atomic,doori,ostep,fine,dolap
      logical oscal
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /coord / xyz(3,numatm)
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap
      common /plspec/ step,stepf,scale,scalf,cut,pmin,pmax,one,
     &                oscal,ostep,fine
      common /grdhlp/ mx3d,mx3d2
      character*80 mapfil
      character*80 grdfil
      common /maphlp/ mapfil,grdfil
      logical denok,opfil
      character lstr*137
      dimension txyz(3),dens(*),denst(*)

      small = 1.0d-4
      istat = 1

      m = linlen(grdfil)
      if (.not.opfil(21,grdfil(1:m),m,0,1,0)) then
         lstr = 'gridfile non-existent'
         goto 100
      endif
      lstr = 'trying '//grdfil(1:m)

      pmax   = -100000.d0
      pmin   =  100000.d0
c      ipsi = 0
c      bonds  = .false.
c      atomic = .false.
c      ovrlap = .false.

      read(21,end=100,err=100) katoms
      if (isubtr.eq.0) then
         natoms = katoms
      else
         if (katoms.ne.natoms) then
            lstr = 'number of read atoms incorrect'
            goto 100
         endif
      endif

      taccu = 0
      do i=1,natoms
         read(21,end=100,err=100)(txyz(j),j=1,3)
         if (isubtr.eq.0) then
            do j=1,3
               xyz(j,i) = txyz(j)
            end do
         else
            taccu = taccu + (vlen(txyz) - vlen(xyz(1,i)))**2
         endif
      end do
      if (isubtr.ne.0) then
         if (taccu.gt.small) then
            lstr = 'cartesian coordinates read incorrect'
            goto 100
         endif
      endif

      read(21,end=100,err=100) 
     &    px, py, pz, cx, cy, cz, r(1),r(2),npts1,npts2

      if (npts1.gt.mx3d.or.npts2.gt.mx3d) then
          lstr = 'dimension greater than maximum !'
          print*,'npts1 ',npts1,' npts2 ',npts2,' maxdim ',mx3d
          goto 100
      endif

      r(3) = r(1)

c      ipsi = 0
c      bonds  = .false.
c      atomic = .false.
c      ovrlap = .false.
      if (ipsi.ne.0.or.
     &       (ipsi.eq.0.and.(bonds.or.atomic.or.ovrlap))) then
         call messg(11)
      endif
c      call pareul
c      call proato
c      call denmak(denok)
c      call grdcal(dens,npts1,npts2,1,0)

      read(21,end=100,err=100) (denst(i),i=1,npts1*npts2)
      do i=1,npts1*npts2
         if (isubtr.eq.0) then
            dens(i) = denst(i)
         else
            dens(i) = dens(i) - denst(i)
         endif
         if (dens(i).gt.pmax) pmax = dens(i)
         if (dens(i).lt.pmin) pmin = dens(i)
      end do

      pmax = max(pmax,-pmin)

      step = 0.005d0

      close(21)

      lstr = 'found '//grdfil(1:m)
      call inferr(lstr,0)
      return

100   close(21)
      istat = 0
      natoms = 0
      call inferr(lstr,0)
      return
      end

      subroutine tstrd3
      implicit double precision (a-h,o-z)
      common /grdhlp/ mx3d,mx3d2
      character*80 mapfil
      character*80 grdfil
      common /maphlp/ mapfil,grdfil

      logical opfil

      small = 1.0d-4
      istat = 1

      m = linlen(grdfil)
      if (.not.opfil(21,grdfil(1:m),m,0,1,0)) return

      rewind(21)
      read(21,end=100,err=100) katoms

      read(21,end=100,err=100) (natt,i=1,katoms)

      read(21,end=100,err=100) adjust

      do i=1,katoms
         read(21,end=100,err=100)(txyz,j=1,3)
      end do

      read(21,end=100,err=100) pxt, pyt, pzt, cxt, cyt, czt, 
     &                 rt,rt,rt, npts1,npts2,npts3,iplt

      nptsmx = npts1
      if (npts2.gt.nptsmx) nptsmx = npts2
      if (npts3.gt.nptsmx) nptsmx = npts3

      if (nptsmx.gt.mx3d) call allgrd(nptsmx)

100   close(21)

      return
      end

      subroutine rd3ind(npts1,npts2,npts3,isubtr,adjus,istat,
     &                  denn,pmnn,denst)
      implicit double precision (a-h,o-z)

c THIS REALLY rd3inf

      parameter (numatm=2000)
      logical valenc,bonds,ovrlap,atomic,doori,ostep,fine,dolap
      logical oscal
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /coord / xyz(3,numatm)
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap
      common /plspec/ step,stepf,scale,scalf,cut,pmin,pmax,one,
     &                oscal,ostep,fine
      common /grdhlp/ mx3d,mx3d2
      logical oclose,osshad,ofill,ofrst
      common /spacom/ grey,iscol,oclose,osshad,ofill,ofrst
      character*80 mapfil
      character*80 grdfil
      common /maphlp/ mapfil,grdfil

      logical opfil
      character lstr*137
      dimension txyz(3),natt(numatm),rt(3),np(3)
      dimension denn(*),pmnn(*),denst(*)

      small = 1.0d-4
      istat = 1

      if ((ipsi.ne.0.or.bonds.or.atomic.or.ovrlap)
     &    .and.isubtr.eq.1) then
          call inferr('Only Normal Density Can be Subtracted!',0)
          istat = 0
          return
      endif

      m = linlen(grdfil)
      if (.not.opfil(21,grdfil(1:m),m,0,1,0)) then
         lstr = '3dgridfile non-existent'
         goto 100
      endif
      rewind(21)

      read(21,end=100,err=100) katoms
      if (isubtr.eq.1) then
         if (katoms.ne.natoms) then
            lstr = 'number of read atoms incorrect'
            goto 100
         endif
      else
         natoms = katoms
      endif

      read(21,end=100,err=100) (natt(i),i=1,natoms)

      read(21,end=100,err=100) adjust
      if (isubtr.eq.0) adjus = adjust
c      print*,'adjust=',adjust

      taccu = 0
      do i=1,natoms
         read(21,end=100,err=100)(txyz(j),j=1,3)
         if (isubtr.eq.1) then
            taccu = taccu + (vlen(txyz) - vlen(xyz(1,i)))**2
         else
            nat(i) = natt(i)
            do j=1,3
                xyz(j,i) = txyz(j)
            end do
         endif
      end do
      if (taccu.gt.small.and.isubtr.eq.1) then
         lstr = 'cartesian coordinates read incorrect'
         goto 100
      endif

      read(21,end=100,err=100) pxt, pyt, pzt, cxt, cyt, czt, 
     &                 rt(1),rt(2),rt(3), np(1),np(2),np(3),iplt

      write(*,'(a,3(f9.3,1x))') 'Center',pxt,pyt,pzt
      write(*,'(a,3(f9.3,1x))') 'Line  ',cxt,cyt,czt
      write(*,'(a,3(f9.3,1x))') 'edge  ',rt(1),rt(2),rt(3)

      if (np(1).gt.mx3d.or.np(2).gt.mx3d.or.np(3).gt.mx3d) then
          lstr = 'dimension greater than maximum !'
          print*,'npts1 ',np(1),' npts2 ',np(2),' npts3 ',np(3),
     &           ' maxdim ',mx3d
          goto 100
      endif


      if (isubtr.eq.1) then
         if (pxt.ne.px.or.pyt.ne.py.or.pzt.ne.pz
     &   .or.cxt.ne.cx.or.cyt.ne.cy.or.czt.ne.cz
     &   .or.rt(1).ne.r(1).or.rt(2).ne.r(2).or.rt(3).ne.r(3)
     &   .or.np(1).ne.npts1.or.np(2).ne.npts2.or.np(3).ne.npts3
     &   .or.iplt.ne.iplat) then
            lstr = 'NOT the same CENTER,LINE or EDGE used !'
            print*,'old px py pz ',px,py,pz,' cx cy cz ',cx,cy,cz
            print*,'    edx,edy,edz ',r(1),r(2),r(3),' iplat ',iplat
            print*,'    npts1,npts2,npts3 ',npts1,npts2,npts3
            print*,'new px py pz ',pxt,pyt,pzt,' cx cy cz ',cxt,cyt,czt
            print*,'    edx,edy,edz ',rt(1),rt(2),rt(3),' iplat ',iplat
            print*,'    npts1,npts2,npts3 ',np(1),np(2),np(3)
            goto 100
         endif
      else
         iplat = iplt
         px = pxt
         py = pyt
         pz = pzt
         cx = cxt
         cy = cyt
         cz = czt
         do i=1,3
            r(i) = rt(i)
         end do
         npts1 = np(1)
         npts2 = np(2)
         npts3 = np(3)

         call pareul
         call parrat
         call proato
         ipsi = 0
         ofrst = .false.
      endif

      do j=1,npts3
         pmnn(j) = 100000.d0
         read(21,end=100)(denst(i),i=1,npts1*npts2)
         do i=1,npts1*npts2
            if (isubtr.eq.1) then
               denn((j-1)*mx3d2 + i) = denn((j-1)*mx3d2 + i) 
     &                                    - denst(i)
            else
               denn((j-1)*mx3d2 + i) = denst(i)
            endif
         end do
      end do

      close(21)

      lstr = 'found '//grdfil(1:m)
      call inferr(lstr,0)
      return

100   close(21)
      natoms = 0
      istat = 0
      call inferr(lstr,0)
      return
      end

      subroutine rd3chd(npts1,npts2,npts3,igauss,impas,istat,fmap)
      implicit double precision (a-h,o-z)
      common /grdhlp/ mx3d,mx3d2
      character*80 mapfil
      character*80 grdfil
      common /maphlp/ mapfil,grdfil
      character lstr*137
      logical opfil,oampac
      dimension np(3),dtmp(6)
      dimension fmap(*)

      istat = 1
      itype = 0
      iform = 0
      if (igauss.eq.1) iform = 1
      iun = 21

      m = linlen(mapfil)
      if (.not.opfil(iun,mapfil(1:m),m,iform,1,0)) then
         lstr = 'mapfile non-existent'
         call inferr(lstr,0)
         call messg(8)
         istat = 0
         return
      endif

      lstr = 'read error'

      if (igauss.eq.1) then
         read(iun,'(a)',end=100) lstr
         read(iun,'(a)',end=100) lstr
         if (index(lstr,'Density').ne.0) itype = 0
         if (index(lstr,'MO coefficients').ne.0) itype = 1
         read(iun,'(a)') lstr
         oampac = (linlen(lstr).eq.44)
         if (oampac) then
            read(lstr,'(i5,3f13.6)') natoms,(rdum,i=1,3)
         else
            read(lstr,'(i5,3f12.6)') natoms,(rdum,i=1,3)
         endif
         if (natoms.lt.0) itype = 1
         natoms = iabs(natoms)

         if (oampac) then
            read(iun,'(i5,3f13.6)') np(1),(rdum,i=1,3)
            read(iun,'(i5,3f13.6)') np(2),(rdum,i=1,3)
            read(iun,'(i5,3f13.6)') np(3),(rdum,i=1,3)
         else
            read(iun,'(i5,3f12.6)') np(1),(rdum,i=1,3)
            read(iun,'(i5,3f12.6)') np(2),(rdum,i=1,3)
            read(iun,'(i5,3f12.6)') np(3),(rdum,i=1,3)
         endif

         do i=1,natoms
            read(iun,'(a)') lstr
         end do

         if (itype.eq.1) then
            read(iun,'(i5)') no
            if (no.gt.1) then
               lstr = 'ONLY single orbital cube files are supported !'
               goto 100
            endif
            nlines = (no + 1) / 10
            if (no+1 - nlines*10.gt.0) nlines = nlines + 1
            do i=1,nlines-1
               read(iun,'(a)') lstr
            end do
         endif

      else
         read(iun,end=100) katoms
         read(iun,end=100) (idum,i=1,katoms)
         read(iun,end=100) rdum
         do i=1,katoms
            read(iun,end=100) (rdum,j=1,3)
         end do
         read(iun,end=100) (rdum,i=1,9),np(1),np(2),np(3),idum
      endif

      if (np(1).ne.npts1.or.np(2).ne.npts2.or.np(3).ne.npts3) then
        print*,'npts1-3 ',npts1,npts2,npts3,' np1-3 ',np(1),np(2),np(3)
          lstr = 'dimensions grid dont match first grids !'
          call messg(10)
          close(iun)
          istat = 0
          return
      endif

      if (impas.eq.0) then

         call almgrd

      else

         lstr = 'read error'

         if (igauss.eq.1) then

            ij = 0
            do i=1,npts1
               do j=1,npts2
                  ij = ij + 1
                  kk = 0
                  do while (kk.lt.npts3)
                     n = 6
                     if (npts3-kk.lt.n) n = npts3-kk
                     read(iun,'(6E13.5)') (dtmp(k),k=1,n)
                     do k=1,n
                        fmap((npts3-(kk+k))*mx3d2 + ij) = dtmp(k)
                     end do
                     kk = kk + 6
                  end do
               end do
            end do

         else

            do j=1,npts3
               read(iun,end=100)(fmap((j-1)*mx3d2 + i),
     &              i=1,npts1*npts2)
            end do

         endif

      endif

      
      close(iun)

      return

100   istat = 0
      close(iun)
      call inferr(lstr,0)
      call messg(10)
      return
      end

