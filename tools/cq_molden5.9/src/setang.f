      subroutine setang(px,py,pz,theta,phi)
      implicit double precision (a-h,o-z)
      common /eul/ ca,cb,sa,sb,ccab,scab,ssab,scba

      todeg = 45.0d0 / datan(1.0d0)
*
* SET UP COS AND SINE ANGLES 
*
      xy = px*px + py*py
      if (xy.gt.1.d-15) then
          r1 = dsqrt(xy + pz*pz)
          xy = dsqrt(xy)
          ca = px / xy
          cb = pz / r1
          sa = py / xy
          sb = xy / r1
      else
          ca = 1.0d0
          cb = 1.0d0
          sa = 0.0d0
          sb = 0.0d0
      endif
      theta = dasin(sa) * todeg
      phi   = dasin(sb) * todeg
      ccab  = ca * cb
      scab  = sa * cb
      ssab  = sa * sb
      scba  = sb * ca

      return
      end
