      subroutine rys(norder,x,root,weight)
      implicit double precision (a-h,o-z)
      dimension root(*),weight(*)
      common /ryscom/ et(3745),ew(3745),cut(2,8),ltw(5),istrt(5)
      common /hermit/ hr(36),hw(36)
      dimension d(6)
 
      if (norder.le.0.or.norder.gt.8) stop 'error in rys'
 
      choose  =  x / (4.0d0 + x)

      if (choose.ge.cut(2,norder)) then

c     hermite polynomial

         indx = (norder*(norder-1))/2
         do i = 1,norder
            indx = indx + 1
            root(i) = hr(indx)**2 / x
            weight(i) = hw(indx) / dsqrt(x)
         end do
         return
      endif

      if (choose.ge.cut(1,norder).or.norder.gt.5) then
         call ryspol(norder,x,root,weight)
         return
      endif

      m      = ltw(norder)*choose
      alfa   = ltw(norder)*choose-dfloat(m)
      beta   = 1.d0-alfa
      at     = alfa*(alfa*alfa-1.d0)/6.0d0
      ap     = beta*(beta*beta-1.d0)/6.0d0
      d(6)   = at*(alfa*alfa-4.d0)*0.05d0
      d(1)   = ap*(beta*beta-4.d0)*0.05d0
      d(2)   =              d(6) - 4.0d0*d(1)            +       ap
      d(3)   = beta - 4.0d0*d(6) + 6.0d0*d(1) +       at - 2.0d0*ap
      d(4)   = alfa + 6.0d0*d(6) - 4.0d0*d(1) - 2.0d0*at +       ap
      d(5)   =      - 4.0d0*d(6) +       d(1) +       at
      indx = istrt(norder)
      do i = 1,norder
         root(i) = evertt(et(indx),d,m)
         weight(i) = evertt(ew(indx),d,m)
         weight(i) = dsqrt(weight(i))
         indx = indx + ltw(norder) + 3
      end do

      return
      end
      double precision function evertt(t,d,m)
      implicit double precision (a-h,o-z)
c
c     everett interpolation
c
      dimension t(*),d(*)
 
      evertt = 0.0d0
      do i=1,6
         evertt = evertt + d(i)*t(m+i)
      end do
 
      return
      end
