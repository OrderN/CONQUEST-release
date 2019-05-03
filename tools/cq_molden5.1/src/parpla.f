      subroutine parpla(iat1,iat2,iat3,istat)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /coord / xyz(3,numatm)
      dimension d1(3),d2(3),c12(3)

      istat = 0
      px = (xyz(1,iat1) + xyz(1,iat2) + xyz(1,iat3)) / 3.0d0
      py = (xyz(2,iat1) + xyz(2,iat2) + xyz(2,iat3)) / 3.0d0
      pz = (xyz(3,iat1) + xyz(3,iat2) + xyz(3,iat3)) / 3.0d0
      do i=1,3
        d1(i) = xyz(i,iat2) - xyz(i,iat1)
        d2(i) = xyz(i,iat3) - xyz(i,iat1)
      end do
      call impsc(d1,d2,cosb)
      cc = dabs(cosb) - 1.0d0
      if (dabs(cc).lt.1.0d-3) istat = 1
      call crprod(d1,d2,c12)
      cx = c12(1)
      cy = c12(2)
      cz = c12(3)
      iplat = iat1

      return
      end
