      subroutine parang(theta,phi)
      implicit double precision (a-h,o-z)
      common /eul/ ca,cb,sa,sb,ccab,scab,ssab,scba

      todeg = 45.0d0 / datan(1.0d0)
      sa    = dsin(theta / todeg)
      ca    = dcos(theta / todeg)
      sb    = dsin(phi / todeg)
      cb    = dcos(phi / todeg) 
      ccab  = ca * cb
      scab  = sa * cb
      ssab  = sa * sb
      scba  = sb * ca

      return
      end
