      subroutine dolift(rlift)
      implicit double precision (a-h,o-z)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      dimension d2(3)

      d2(1) = cx
      d2(2) = cy
      d2(3) = cz

      d2norm = dsqrt(d2(1)*d2(1)+d2(2)*d2(2)+d2(3)*d2(3))
      if (d2norm.lt.1.0d0-8) return

      do i=1,3
         d2(i) = d2(i)*rlift/d2norm
      end do

      px = px + d2(1)
      py = py + d2(2)
      pz = pz + d2(3)

      return
      end
