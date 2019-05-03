      subroutine defpc
c     count atoms in the xy,xz,yz planes
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)

      px  = 0.0d0
      py  = 0.0d0
      pz  = 0.0d0

      cx  = 0.0d0
      cy  = 0.0d0
      cz  = 0.0d0

      nxy = 0
      nxz = 0
      nyz = 0

      do l=1,natoms
         xinpxy    = xyz(1,l)*.0d0 + xyz(2,l)*.0d0 + xyz(3,l)*1.d0
         xinpxy    = xinpxy * xinpxy
         xinpxz    = xyz(1,l)*.0d0 + xyz(2,l)*1.d0 + xyz(3,l)*.0d0
         xinpxz    = xinpxz * xinpxz
         xinpyz    = xyz(1,l)*1.d0 + xyz(2,l)*.0d0 + xyz(3,l)*.0d0
         xinpyz    = xinpyz * xinpyz
         if (xinpxy.lt.1.d-10) nxy = nxy + 1
         if (xinpxz.lt.1.d-10) nxz = nxz + 1
         if (xinpyz.lt.1.d-10) nyz = nyz + 1
      end do

      if (natoms.gt.2) then
         if (nxy.ge.3.or.nxz.ge.3.or.nyz.ge.3) then
            iplat = 0
            if (nxy.ge.nxz.and.nxy.ge.nyz) then
                cz = 1.0d0
            elseif (nxz.ge.nxy.and.nxz.ge.nyz) then
                cy = 1.0d0
            else
                cx = 1.0d0
            endif
         else
            call parpla(1,2,3,istat)
            if (istat.eq.1) then
               do i=4,natoms
                  call parpla(1,2,i,istat)
                  if (istat.eq.0) goto 100
               end do
            endif
         endif
100      continue
      else
         iplat = 0
         if (nxy.eq.2) then
            cz = 1.0d0
         elseif (nxz.eq.2) then
            cy = 1.0d0
         elseif (nyz.eq.2) then
            cz = 1.0d0
         else
            cx = xyz(1,2) - xyz(1,1)
            cy = xyz(2,2) - xyz(2,1)
            cz = xyz(3,2) - xyz(3,1)
         endif
      endif

      return
      end
