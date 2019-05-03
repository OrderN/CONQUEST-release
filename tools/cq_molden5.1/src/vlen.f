      double precision function vlen(a)
      implicit double precision (a-h,o-z)
      dimension a(3)

      tot = a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
      vlen = 0.0d0
      if (tot.gt.0.0d0) vlen = dsqrt(tot)

      return
      end
