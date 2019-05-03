      subroutine spacdd(npts1,npts2,npts3,valcnt,idofil,
     &                  adjus,ipsprt,idisml,idvrml,mapit,
     &                  denn,pmnn,iedlog,fmap)
      implicit double precision (a-h,o-z)

c THIS IS REALLY spaced

      parameter (numatm=2000)
      parameter (mxel=100)
      logical euclid, yes, euctmp
      logical oscal,ostep,fine,lapdbl
      real xx, yy
      integer*2 ixx
      character str*100, esc
      character monstr*100
      character*2 gstr
      character*80 vfile
      logical oclose,osshad,ofill,ofrst
      common /spacom/ grey,iscol,oclose,osshad,ofill,ofrst
      common /fill /  ixx(1000),nixx,ihight,icol1
      common /plot/   iplot,iplwin,icolps
      common /cntval/ cntval, euclid, yes, icnt, iflag
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /grdhlp/ mx3d,mx3d2
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /plspec/ step,stepf,scale,scalf,cut,pmin,pmax,one,
     &                oscal,ostep,fine
      common /eul/    ca,cb,sa,sb,ccab,scab,ssab,scba
      common /atopla/ xsym(numatm),ysym(numatm),zsym(numatm),
     &                isym(numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /zproj/  zval(numatm),nindx(numatm),nconn(numatm)
      common /coord / xyz(3,numatm)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      logical valenc,bonds,ovrlap,atomic,doori,dolap,doelf
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap,doelf
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /lapcom/ lapbox,lapdbl
      common /orbhlp/ mxorb,iuhf,ispd
      common /vropt/ ivtwo,ihand,ivadd
      character*80 spcdfil
      common /vrmlhl/ spcdfil
      dimension origin(3),rt(3)
      dimension denn(*),pmnn(*),iedlog(*),fmap(*)

      esc = char(27)

      if (npts1.gt.mx3d.or.npts2.gt.mx3d.or.npts3.gt.mx3d) 
     &    stop 'this routine only handles 61 points'

      if (dolap.and.valcnt.eq.0.0d0) valcnt = -0.0001

      if (ofrst) print*,'Contour Value = ',valcnt

      osshad = .true.
      nixx = 0
      if (iplot.eq.6) call gethei(ihight)
      icol1 = 15

      euctmp = euclid
      euclid = .false.
      pxt = px
      pyt = py
      pzt = pz
      pmint = pmin
      pmaxt = pmax

c     This will stop the Ultrix compiler from optimising pmaxt away
c     We could use -O0 compiler flag but this keeps the makefile simple
      pmax = pmaxt

      if (idofil.eq.1) ofill = .true.

      rf3 = r(3)/dble(npts3-1)
      iflag = 0
      if (lapbox.eq.0) then
         l1 = 1
         l2 = npts3
      elseif (lapbox.eq.1) then
         l1 = 1
         l2 = npts3/2 + 1
      elseif (lapbox.eq.2) then
         l1 = npts3/2 + 1
         l2 = npts3
      endif
         
     
      if (ofrst) then 
         if (iplot.eq.6) call curs(1)
         do i=1,natoms
            nconn(i) = 0
            do j=1,natoms
              if (i.ne.j) then
               dmaxsq = (vdwr(nat(i)) + vdwr(nat(j)))**2
               dijsq = ((xyz(1,i)-xyz(1,j))*adjus)**2
     +                +((xyz(2,i)-xyz(2,j))*adjus)**2
     +                +((xyz(3,i)-xyz(3,j))*adjus)**2
               if (dijsq.lt.dmaxsq) then
                   nconn(i) = nconn(i) + 1
               endif
              endif
            end do
        end do

        call precal(npts1,npts2)

        do i=l1,l2
            px = pxt
            py = pyt
            pz = pzt
            call precz(r(3),npts3,i)
            call dolift(0.5d0*r(3) - dble(i-1)*rf3)
            if (iplot.eq.6.or.ipsprt.eq.1) then
               monstr = 'Progress Monitor [1-'//gstr(l2)//'] '//gstr(i)
               call inferr(monstr,1)
            endif
            call grdcal(denn( (i-1)*mx3d2 + 1),npts1,npts2,0,1)
            pmnn(i) = pmax
         end do
         if (iplot.eq.6) call inferr('Space is Done !!',0)
         ofrst = .false.
         if (iplot.eq.6) call curs(0)
      endif

      do i=1,natoms
         zval(i) = scba*(0.5d0 + xsym(i) / r(1)) +
     &             ssab*(0.5d0 + ysym(i) / r(2)) + cb*zsym(i) / r(3)
      end do

      call shsort(natoms,zval,nindx)

      do i=l1,l2
        height = 0.5d0 - dble(i-1)/dble(npts3-1)
        if (dabs(valcnt).le.pmnn(i)) then
         call selsol
         if (iplot.eq.6) then
            xx = 2.0
            yy = 0.0
            call xwin(xx,yy,99,str,nstr,idum1,idum2)
            call xwin(xx,yy,8,str,nstr,idum1,idum2)
         elseif (iplot.eq.4.and.icolps.eq.1) then
            write(iun4,*) 'poscontour setcol'
         endif
         grey = 1.0
         cntval = valcnt
         call cntour(denn( (i-1)*mx3d2 + 1),npts2,npts2,npts1,
     &               height,cntval,x,iedlog)

c------ select other line type 

         if (iplot.eq.0) then
             write(iun4,'(a)') '.nc 3'
             write(iun4,'(a)') '.ls 2 0.2'
         endif
         if (iplot.eq.1) write(iun4,'(a)') 'SP3;'
         if (iplot.eq.2) write(iun4,'(a,''*m7b'')')esc
         if (iplot.eq.3) then
             idum = 1
             call plotgh(idum,0.0d0,0.0d0)
             write(iun4,*) esc//'a'
         endif
         if (iplot.eq.4) then
             write(iun4,'(''s'')')
             if (icolps.eq.1) then
                write(iun4,*) 'negcontour setcol'
             else
                write(iun4,'(''0 setlinecap [4 7] 0 setdash'')')
             endif
             write(iun4,'(''n'')')
         endif
         if (iplot.eq.6) then
             xx = 1.0
             yy = 0.0
             call xwin(xx,yy,99,str,nstr,idum1,idum2)
             call xwin(xx,yy,9,str,nstr,idum1,idum2)
             idum = 1
             call plotgh(idum,0.0d0,0.0d0)
         endif

         if (osshad) then
            grey = 0.8
         else
            grey = 1.0
         endif
         cntval = -1.0d0*valcnt
         if (.not.(dolap.and..not.lapdbl))
     &   call cntour(denn( (i-1)*mx3d2 + 1),npts2,npts2,npts1,
     &               height,cntval,x,iedlog)

       endif
       if (ofill.and.idisml.eq.1) then
           call selsol
           call pl3dm(adjus,.true.,height)
       endif
      end do
 
      if (iplot.eq.4) then
        write(iun4,'(''s'')')
        write(iun4,'(''   3 setlinewidth'')')
        write(iun4,'(''n'')')
      endif
      if (iplot.eq.3) then
         idum=1
         call plotgh(idum,0.0d0,0.0d0)
      endif
      if (iplot.eq.6) then
         idum=1
         call plotgh(idum,0.0d0,0.0d0)
      endif
      call selsol
      if ((idisml.eq.1.or.(idisml.eq.0.and.iplot.ne.6)).and..not.ofill)
     &     call pl3dm(adjus,.false.,rdum)

      euclid = euctmp 
      cntval = 99.999d0
      px = pxt
      py = pyt
      pz = pzt
      pmin = pmint
      pmax = pmaxt
      call parstp

      ofill = .false.

      if (idvrml.eq.1) then
         vfile = spcdfil
         if (iplot.eq.6) call curs(1)
         if (ivtwo.eq.3) then
            idum = 1
         elseif (ivtwo.eq.2) then
            call inferr('Generating POVRAY file',0)
         elseif (ivtwo.eq.0.or.ivtwo.eq.1) then
            call inferr('Generating VRML',0)
         endif
         itrans = 0
         if (ipsi.eq.0.and..not.bonds.and..not.dolap.and..not.elpot
     &       .and..not.molpot.and..not.chpot.and.ispd.eq.0) itrans = 1
         if (ivtwo.eq.4) then
            do i=1,3
               rt(i) = 1.0
               origin(i) = 0.0d0
            end do
            ny = npts1
            nx = npts2
            nz = npts3
            call cvtcom
         else
            origin(1) = -0.5d0 * r(2) / r(1)
            origin(2) = -0.5d0
            origin(3) = -0.5d0 * r(3) / r(1)
            do i=1,3
               rt(i) = r(i)
            end do
            nx = npts2
            ny = npts1
            nz = npts3
         endif


         if (ivtwo.eq.3) then
            iwrt = 1
            call initog(iwrt)
            call oginsp(rt,adjus,natoms,nat,idum,icol,xsym,ysym,zsym,
     &                  vdwr,idum,idum,idum,idum,iwrt)
         elseif (ivtwo.eq.4) then
            idum = 1
         else

            iun = 25

            open(unit=iun,file=vfile,form='formatted',
     &        status='unknown',err=1000)
         endif

         if (lapbox.eq.0) then
            call mcubes(nx,ny,nz,denn,fmap,mapit,origin,valcnt,
     &                      -valcnt,rt,ipsi,dolap,lapdbl,itrans,iun)
         elseif (lapbox.eq.1) then
            nz = npts3/2 + 1
            call mcubes(nx,ny,nz,denn,fmap,mapit,origin,valcnt,
     &                      -valcnt,rt,ipsi,dolap,lapdbl,itrans,iun)
         elseif (lapbox.eq.2) then
            nz = npts3/2 + 1
            call mcubes(nx,ny,nz,denn( (nz-1)*mx3d2 + 1),
     &                  fmap( (nz-1)*mx3d2 + 1),mapit,origin,
     &                  valcnt,-valcnt,rt,ipsi,dolap,lapdbl,itrans,iun)
         endif

         if (ivtwo.eq.3) then
            call ogspst
         elseif (ivtwo.eq.4) then
            idum = 1
         else
            if (idisml.eq.1.or.(idisml.eq.0.and.iplot.ne.6))
     &        call plvmol(iun,npts2,adjus)
c
c --- Close out the object, and display it --- FPA 2/27/2000
c
            if (ivtwo .eq. 2) then
              write(iun,*) '}'
              write(iun,*) 'object {molecule}'
            end if
c
            close(iun)
            if (ivtwo.ne.4) call inferr('Wrote file: '//vfile,0)
         endif
         if (iplot.eq.6) call curs(0)
c
c write file points with xyz coordinates of each points of the 3d grid
c         if (ivtwo.eq.4) call wcubev(nx,ny,nz,denn,origin,rt)

         idvrml = 0
      endif

      mapit = 0
      return

1000  call inferr('Couldnt open file: '//vfile,0)
      mapit = 0
      return
      end

      subroutine spasrd(npts1,npts2,npts3,valcnt,
     &                  denn,pmnn,dens,iedlog)
      implicit double precision (a-h,o-z)

c THIS IS REALLY spasrf

      parameter (numatm=2000)
      logical euclid, yes, euctmp
      logical oscal,ostep,fine,lapdbl
      logical oclose,osshad,ofill,ofrst
      common /spacom/ grey,iscol,oclose,osshad,ofill,ofrst
      common /cntval/ cntval, euclid, yes, icnt, iflag
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /grdhlp/ mx3d,mx3d2
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /plspec/ step,stepf,scale,scalf,cut,pmin,pmax,one,
     &                oscal,ostep,fine
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)
      logical valenc,bonds,ovrlap,atomic,doori,dolap,doelf
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap,doelf
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /lapcom/ lapbox,lapdbl
      common /pntsta/ iproj,idrcol,nesp
      dimension vt1(3),vt2(3),rt(3)
      dimension denn(*),pmnn(*),dens(*),iedlog(*)

      if (ofrst)  return

      if (dolap.and.valcnt.eq.0.0d0) valcnt = -0.0001

      call pareul
      euctmp = euclid
      euclid = .false.
      pxt = px
      pyt = py
      pzt = pz
      pmint = pmin
      pmaxt = pmax
      do i=1,3
         vt1(i) = v1(i)
         vt2(i) = v2(i)
      end do
      cxt = cx
      cyt = cy
      czt = cz


c     This will stop the Ultrix compiler from optimising pmaxt away
c     We could use -O0 compiler flag but this keeps the makefile simple
      pmax = pmaxt

      rf3 = r(3)/dble(npts3-1)
      iproj = 0
      iflag = 1
      height = 0.0d0
      if (lapbox.eq.0) then
         l1 = 1
         l2 = npts3
      elseif (lapbox.eq.1) then
         l1 = 1
         l2 = npts3/2 + 1
      elseif (lapbox.eq.2) then
         l1 = npts3/2 + 1
         l2 = npts3
      endif
         
      call precal(npts1,npts2)

      do i=l1,l2
        px = pxt
        py = pyt
        pz = pzt
        call precz(r(3),npts3,i)
        call dolift(0.5d0*r(3) - dble(i-1)*rf3)

        if (dabs(valcnt).le.pmnn(i)) then
         idrcol = 11
         cntval = valcnt
         call cntour(denn( (i-1)*mx3d2 + 1),npts2,npts2,npts1,
     &               height,cntval,x,iedlog)

         cntval = -1.0d0*valcnt
         if (.not.(dolap.and..not.lapdbl)) then
            idrcol = 1
            call cntour(denn( (i-1)*mx3d2 + 1),npts2,npts2,npts1,
     &                  height,cntval,x,iedlog)
         endif

        endif
      end do
 
      rf1 = r(1)/dble(npts1-1)
      cx = vt1(1)
      cy = vt1(2)
      cz = vt1(3)
      v1(1) = cxt
      v1(2) = cyt
      v1(3) = czt
      do i=1,3
         rt(i) = r(i)
      end do
      r(1) = rt(3)
      r(3) = rt(1)

      call vsc1(v1,1.0d0,1.0d-4)

      do i=1,npts1
         px = pxt
         py = pyt
         pz = pzt
         ij = 0
         do k=1,npts3
            do j=1,npts2
               ij = ij + 1
               dens(ij) = 
     &           denn( (npts3-k)*mx3d2 + (j+npts2*((npts1-i+1)-1)))
            end do
         end do
         call precz(rt(1),npts1,i)
         call dolift(0.5d0*rt(1) - dble(i-1)*rf1)

         cntval = valcnt
         idrcol = 11
         call cntour(dens,npts2,npts2,npts3,height,cntval,x,iedlog)
         cntval = -1.0d0*valcnt
         if (.not.(dolap.and..not.lapdbl)) then
            idrcol = 1
            call cntour(dens,npts2,npts2,npts3,height,cntval,x,iedlog)
         endif
      end do

100   continue
      euclid = euctmp 
      cntval = 99.999d0
      px = pxt
      py = pyt
      pz = pzt
      pmin = pmint
      pmax = pmaxt
      cx = cxt
      cy = cyt
      cz = czt
      do i=1,3
         v1(i) = vt1(i)
         v2(i) = vt2(i)
      end do
      do i=1,3
         r(i) = rt(i)
      end do

      iproj = 1

      call docent
      call doscal
      call setxyv

      return
      end

      subroutine isoded(valc,nvalc,scincr,nespt,iwhere,
     &                  dens,iedlog)
      implicit double precision (a-h,o-z)

c THIS IS REALLY isoden

      parameter (lnbuck=10)
      parameter (max3d=61)
      parameter (mesp = (max3d*max3d*max3d+max3d)/4)
      parameter (mdum = (max3d*max3d*max3d+max3d) - mesp*4)
      parameter (mxvalc=10)
      logical euclid, yes, euctmp
      common /cntval/ cntval, euclid, yes, icnt, iflag
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      logical oscal,ostep,fine
      common /plspec/ step,stepf,scale,scalf,cut,pmin,pmax,one,
     &                oscal,ostep,fine
      common /pntsta/ iproj,idrcol,nesp
      common /spa3d/  connl(3,mesp), esp(mesp), dum(mdum)
      common /srfhlp/ edge,ctval(mxvalc),pxyz(3),mvalc,nspts,istyp
      common /plot/   iplot,iplwin,icolps
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension vt1(3),vt2(3),valc(*),rtmp(3)
      dimension dens(*),iedlog(*)

      if (nspts.gt.max3d) then
         call inferr('To many grid points, adjusting !',0)
         nspts = max3d
      endif

      call pareul
      euctmp = euclid
      euclid = .false.
      pxtt = px
      pytt = py
      pztt = pz
      px = pxyz(1)
      py = pxyz(2)
      pz = pxyz(3)
      pxt = px
      pyt = py
      pzt = pz
      pmint = pmin
      pmaxt = pmax
      do i=1,3
         rtmp(i) = r(i)
         r(i) = edge
         vt1(i) = v1(i)
         vt2(i) = v2(i)
      end do
      cxt = cx
      cyt = cy
      czt = cz
      iflag = 1

      rf3 = r(3)/(nspts-1)

      iproj = iwhere
      idrcol = 13
      nesp = 0
      height = 0.0d0
      if (iplot.eq.6) call curs(1)

      call precal(nspts,nspts)

      do i=1,nspts
         px = pxt
         py = pyt
         pz = pzt
         call precz(r(3),nspts,i)
         call dolift(0.5d0*r(3) - dble(i-1)*rf3)
         call grdcal(dens,nspts,nspts,0,1)
         do j=1,nvalc
            call cntour(dens,nspts,nspts,nspts,height,valc(j),x,iedlog)
         end do
      end do


      rf1 = r(1)/(nspts-1)
      cx = vt1(1)
      cy = vt1(2)
      cz = vt1(3)
      v1(1) = cxt
      v1(2) = cyt
      v1(3) = czt
      call vsc1(v1,1.0d0,1.0d-4)

      px = pxt
      py = pyt
      pz = pzt

      call precal(nspts,nspts)

      do i=1,nspts
         px = pxt
         py = pyt
         pz = pzt
         call precz(r(1),nspts,i)
         call dolift(0.5d0*r(1) - dble(i-1)*rf1)
         call grdcal(dens,nspts,nspts,0,1)
         do j=1,nvalc
            call cntour(dens,nspts,nspts,nspts,height,valc(j),x,iedlog)
         end do
      end do

      rf2 = r(2)/(nspts-1)
      cx = vt2(1)
      cy = vt2(2)
      cz = vt2(3)
      v1(1) = vt1(1)
      v1(2) = vt1(2)
      v1(3) = vt1(3)
      v2(1) = cxt
      v2(2) = cyt
      v2(3) = czt
      call vsc1(v2,1.0d0,1.0d-4)

      px = pxt
      py = pyt
      pz = pzt
      call precal(nspts,nspts)

      do i=1,nspts
         px = pxt
         py = pyt
         pz = pzt
         call precz(r(2),nspts,i)
         call dolift(0.5d0*r(2) - dble(i-1)*rf2)
         call grdcal(dens,nspts,nspts,0,1)
         do j=1,nvalc
            call cntour(dens,nspts,nspts,nspts,height,valc(j),x,iedlog)
         end do
      end do

      if (iplot.eq.6) call curs(0)
      iproj = 1

      euclid = euctmp
      px = pxtt
      py = pytt
      pz = pztt
      pmin = pmint
      pmax = pmaxt
      do i=1,3
         r(i) = rtmp(i)
         v1(i) = vt1(i)
         v2(i) = vt2(i)
      end do
      cx = cxt
      cy = cyt
      cz = czt

      nespt = nesp

      return
      end

      subroutine wrsrd(iun,nesp,iesp,
     &                 ianz,iatclr,iconn,coo)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxel=100)
      parameter (mxcon=10)
      parameter (max3d=61)
      parameter (mxesp=max3d*max3d*max3d+max3d)
      common /surf/ natorg,nosncd
      common /spa3d/ esp(mxesp)
      dimension icnn(mxcon)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*)

c write surface file

      toang = 0.52917706d0

      write(iun,'(a)') '[Molden Format]'
      write(iun,'(a)') '[SURFACE]'
      write(iun,'(i5,1x,i5,1x,i1)') nesp,natorg,iesp

      do i=1,nesp

         ibnds = 0

         do j=1,iconn(1,i)

            if (iconn(j+1,i).gt.0) then
               ibnds = ibnds + 1
               icnn(ibnds) = iconn(j+1,i)
            endif

         end do

         if (i.gt.natorg) then

            if (iesp.eq.1) then
              write(iun,'(i3,1x,3(f12.6),1x,f12.6,1x,6(i5,1x))')
     &        ianz(i),(coo(j,i)*toang,j=1,3),
     &        esp(i),(icnn(j),j=1,ibnds)
            else
              write(iun,'(i3,1x,3(f12.6),1x,i2,1x,6(i5,1x))')
     &        ianz(i),(coo(j,i)*toang,j=1,3),
     &        iatclr(i),(icnn(j),j=1,ibnds)
            endif

         else

            write(iun,'(i3,1x,3(f12.6),1x,6(i5,1x))')
     &        ianz(i),(coo(j,i)*toang,j=1,3),
     &        (icnn(j),j=1,ibnds)

         endif

      end do

      return
      end

      subroutine rdsrd(iun,istats,iesp,iaddprv,idebug,
     &                 ianz,iaton,iatclr,iresid,iconn,coo)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (max3d=61)
      parameter (mxesp=max3d*max3d*max3d+max3d)
      common /athlp/ iatoms, mxnat
      common /surf/  natorg,noscnd
      common /spa3d/ esp(mxesp)
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character*137 line,str
      common /curlin/ line
      logical gnreal
      integer getlin
      dimension r(4)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*),
     &          iresid(*),iaton(*)

c read surface file

      istats = 1
      iuntmp = iun2
      iun2 = iun

      toang = 0.52917706d0

      call rewfil
      call searchu(line,'[SURFACE]',istat)
      if (istat.eq.0) goto 110
      iesp = 0

      if (getlin(0).eq.1) then
          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.2) then
             natoms = itype
             if (natoms.gt.mxnat) then
                 natoms = mxnat
                 call inferr('Exceeding MaxNum of atoms !',0)
             endif
             if (idebug.eq.1) print*,'iatoms=',natoms
             if (iaddprv.eq.1) then
                if (natoms+iatoms.gt.mxnat) then
                    call inferr('Exceeding MaxNum of atoms !',0)
                    goto 100
                endif
                iscst = iatoms
                nscnd = natoms
                ialtyp = 0
             else
                iscst = 0
             endif
          else 
             goto 100
          endif
          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.2) then
             natot = itype
             if (idebug.eq.1) print*,'natorg=',natot
          else 
             goto 100
          endif
          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.2) then
             iesp = itype
             if (idebug.eq.1) print*,'iesp=',iesp
          else 
             iesp = 1
          endif
      endif

      do i=1+iscst,natoms+iscst

         if (getlin(0).eq.1) then

             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.2) then
                ianz(i) = itype
             else
                goto 100
             endif

             if (gnreal(r,3,.false.)) then
                do j=1,3
                   coo(j,i) = r(j) / toang
                end do
             else
                goto 100
             endif

             if (i.gt.natot+iscst) then
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (iesp.eq.1.and.ktype.eq.3) then
                    esp(i) = rtype
                    iatclr(i) = 12
                elseif (iesp.eq.0.and.ktype.eq.2) then
                    iatclr(i) = itype
                else
                   goto 100
                endif
             else
                iatclr(i) = 12
             endif

             nc = 0
             do j=1,mxcon
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.2) then
                   if (itype.le.mxnat) then
                      iconn(j+1,i) = itype + iscst
                      nc = nc + 1
                   endif
                endif
             end do
             iconn(1,i) = nc
            
             iresid(i) = -4
             iaton(i) = 1

         else
            goto 100
         endif

      end do

      if (iaddprv.eq.1) then
         iatoms = iatoms + natoms
      else
         iatoms = natoms
         natorg = natot
      endif

      iun2 = iuntmp

      call doscal

      return

100   call inferr('Error reading surface file !',0)
      if (idebug.eq.1) print*,'line=',line
110   iuntmp = iun2
      istats = 0

      return
      end

