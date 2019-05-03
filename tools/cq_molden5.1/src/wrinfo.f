      subroutine wrinfd(npts1,npts2,dens)
      implicit double precision (a-h,o-z)

c THIS IS REALLY wrinfo

      parameter (numatm=2000)
      logical valenc,bonds,ovrlap,atomic,doori,dolap
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /coord / xyz(3,numatm)
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap
c      logical denok
      dimension dens(*)

c      ipsi = 0
c      bonds  = .false.
c      atomic = .false.
c      ovrlap = .false.
c      call denmak(denok)
c      call grdcal(dens,npts1,npts2,1,0)

      if (ipsi.ne.0.or.
     &       (ipsi.eq.0.and.(bonds.or.atomic.or.ovrlap))) then
         call messg(11)
      endif

      open(unit=21,form='unformatted',file='gridfile',status='unknown')

      write(21) natoms
      do i=1,natoms
         write(21)(xyz(j,i),j=1,3)
      end do
      write(21) px, py, pz, cx, cy, cz, r(1),r(2),npts1,npts2
      write(21)(dens(i),i=1,npts1*npts2)

      close(21)

      if (ipsi.eq.0.and..not.(bonds.or.atomic.or.ovrlap)) then
         call inferr('Normal Density written',0)
      endif

      return
      end

      subroutine wr3ind(npts1,npts2,npts3,adjus,denn)
      implicit double precision (a-h,o-z)

c THIS IS REALLY wr3inf

      parameter (numatm=2000)
      logical valenc,bonds,ovrlap,atomic,doori,dolap
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /coord / xyz(3,numatm)
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap
      common /grdhlp/ mx3d,mx3d2
      dimension denn(*)

      if (ipsi.ne.0.or.bonds.or.atomic.or.ovrlap) then
          call inferr('Only Normal Density Can be Written!',0)
          return
      endif

      open(unit=21,form='unformatted',file='3dgridfile',
     &     status='unknown')

      write(21) natoms
      write(21) (nat(i),i=1,natoms)
      write(21) adjus
      do i=1,natoms
         write(21)(xyz(j,i),j=1,3)
      end do
      write(21) px, py, pz, cx, cy, cz, r(1),r(2),r(3),
     &          npts1,npts2,npts3,iplat
      do i=1,npts3
         write(21)(denn((i-1)*mx3d2 + j),j=1,npts1*npts2)
      end do

      close(21)
      call inferr('3D Density written',0)

      return
      end
