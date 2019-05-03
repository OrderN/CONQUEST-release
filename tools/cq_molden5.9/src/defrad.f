      subroutine defrad(doiso)
c     supplies default radius ( see keyword EDGE )
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxvalc=10)
      parameter (max3d=61)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)
      common /srfhlp/edge,ctval(mxvalc),pxyz(3),nvalc,nspts,istyp
      logical doiso
      dimension vec(3)

      if (doiso) then

         call cntvc2(pxyz,xyz,natoms)
         edge = 0.0d0
         do i=1,natoms
            dist = dist2(pxyz,xyz(1,i))
            if (dist.gt.edge) edge = dist
         end do
         edge = dsqrt(edge)*2.0d0

         edge = edge + 7.0d0

         nspts = edge / 0.3d0
         if (nspts.gt.max3d) nspts = max3d

      else

         call cntvec(vec,xyz,nat,natoms)
         r(1) = 0.0d0
         do i=1,natoms
            dist = (xyz(1,i)-px)**2 + 
     &             (xyz(2,i)-py)**2 + 
     &             (xyz(3,i)-pz)**2
            if (dist.gt.r(1)) r(1) = dist
         end do
         r(1) = dsqrt(r(1))*2.4d0
         r(2) = r(1)
         r(3) = r(1)

      endif

      return
      end
