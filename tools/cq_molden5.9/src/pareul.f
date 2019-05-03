      subroutine pareul
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /coord / xyz(3,numatm)
      common /eulx/   ca,cb,sa,sb,cc,sc
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      dimension t(3)

      if (cx.eq.0.d0.and.cy.eq.0.d0) cz = 1.d0
c
c calculate direction cosines and sines that determine the 
c orientation of the plot plane with respect to the coordinate
c system
c
      xy = cx*cx + cy*cy
      if (xy.gt.1.d-15) then
          r1 = dsqrt(xy+cz*cz)
          xy = dsqrt(xy)
          ca = cx/xy
          cb = cz/r1
          sa = cy/xy
          sb = xy/r1
      else
          ca = 1.d0
          cb = 1.d0
          sa = 0.d0
          sb = 0.d0
      endif

      v1(1) = ca*cb
      v1(2) = sa*cb
      v1(3) = -sb

      v2(1) = -sa
      v2(2) = ca
      v2(3) = 0.0d0

      if (iplat.gt.0) then
         t(1) = xyz(1,iplat) - px
         t(2) = xyz(2,iplat) - py
         t(3) = xyz(3,iplat) - pz
         call vsc1(t,1.0d0,1.0d-4)
         call impsc(v1,t,cc)
         call impsc(v2,t,sc)
         do j=1,3
            v1(j) = t(j)
         end do
         t(1) = cx
         t(2) = cy
         t(3) = cz
         call crprod(t,v1,v2)
         call vsc1(v2,1.0d0,1.0d-4)
      else
         cc = 1.0d0
         sc = 0.0d0
      endif

      return
      end

      subroutine parrat
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /plrat/ rat(3)

c      print*,'edx,edy,edz, ',(r(i),i=1,3)

      rat(1) = 1.0d0
      rat(2) = r(2) / r(1)
      rat(3) = r(3) / r(1)

      if (rat(2).gt.1.0d0) then
         rat(1) = r(1) / r(2)
         rat(2) = 1.0d0
         rat(3) = r(3) / r(2)
      endif

      return
      end
