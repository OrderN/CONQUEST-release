      subroutine fmtset
      implicit double precision (a-h,o-z)
      common /fm/    c(17),dflt(17),cutoff
      dimension ftmp(15)

      c(1) = dsqrt(4.0d0*datan(1.0d0))

      tmp = 0.5d0
      do i=2,17
         c(i) = c(i-1)*tmp
         tmp = tmp + 1.0d0
      end do

      tmp = 1.0d0
      dflt(1) = 1.0d0
      do i = 2,17
         tmp = tmp + 2.0d0
         dflt(i) = 1.0d0 / tmp
      end do

      tmp = 20.0d0
10    call fmtgen(ftmp,tmp,1,istatus)
      if (istatus.eq.0) then
         tmp = tmp - 1.0d0
         if (tmp-10.0d0.lt.0.0d0) then
            return
         else
            goto 10
         endif
      endif
      cutoff = tmp + 1.0d0

      return
      end
      subroutine fmtgen(f,t,m,istatus)
      implicit double precision (a-h,o-z)
      dimension f(m)
      common /fm/  c(17),dflt(17),cutoff

      tol = 10.0d-66
c      tol = 10.0d-39
c      VMS only tolerates this value
      istatus = 0

      if (dabs(t).le.0.0d0) then
         do i=1,m
            f(i) = dflt(i)
         end do
         return
      endif

      exx = 0.0d0
      if (dabs(t).ge.50.0d0) then
         f(m) = 0.5d0*c(m) / (t**(dfloat(m) - 0.5d0))
         goto 200
      endif

      exx = dexp(-t)
      if (dabs(t).lt.cutoff) then
         a = dfloat(m-1) + 0.5d0
         deel = 1.0d0 / a
         sig = deel
         do i=2,400
            a = a + 1.0d0
            deel = deel*t / a
            sig = sig + deel
            if (dabs(deel/sig)-tol.lt.0.0d0) then
                f(m) = 0.5d0*sig*exx
                goto 200
            endif
         end do
      endif

      b = dfloat(m-1) + 0.5d0
      a = dfloat(m-1) - 0.5d0
      ca = 0.5d0*dsqrt(4.0d0*datan(1.0d0))*(t**(1-m))/dsqrt(t)
      if (m-1.ne.0) then
         do i=1,m-1
            b = b - 1.0d0
            ca = ca*b
         end do
      endif
      wim = 0.5d0*exx / t
      sig = 0.0d0
      if (wim.eq.0.0d0) goto 100
      wp = wim / ca
      deel = 1.0d0
      sig = 1.0d0
      no = dint(t) + m - 1
      do i=2,no
         deel = deel*a / t
         sig = sig + deel
         if (dabs(deel*wp / sig) - tol.le.0.0d0) goto 100
         a = a - 1.0d0
      end do
      istatus = 1
      return

100   f(m) = ca - wim*sig

200   tmp = dfloat(2*m-3)
      if (m-1.eq.0) return

      do i=1,m-1
         f(m-i) = (2.0d0*t*f(m-i+1) + exx) / tmp
         tmp = tmp - 2.0d0
      end do

      return
      end
