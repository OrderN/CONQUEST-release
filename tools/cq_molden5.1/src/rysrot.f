      subroutine rysrot(ndeg,a,r,ap)
      implicit double precision (a-h,o-z)
      dimension a(*),r(*),ap(*)

      step = .01d0
      call fdpol(ndeg,nfd,a,ap)
      nroots = ndeg / 2
      iroot = 1
      r(1) = 1.0d0

10    x = r(iroot) - .0001d0
      sgn1 = fpol(ndeg,x,a)

20    x = x - step
      if (x.lt.0.0d0) stop 'x is neg. in rysrot'
      sgn2 = fpol(ndeg,x,a)
      if (sign(1.d0,sgn1).eq.sign(1.d0,sgn2)) goto 20

      r(iroot) = x + .5d0*step

      do j=1,30
         f = fpol(ndeg,r(iroot),a)
         fp = fpol(nfd,r(iroot),ap)
         if (dabs(fp).lt.1.d-10) goto 60
         extrap = f / fp
         if (dabs(extrap/r(iroot)).le.1.d-5) goto 60
         r(iroot) = r(iroot) - extrap
      end do

60    iroot = iroot + 1
      if (iroot.gt.nroots) goto 70

      r(iroot) = r(iroot-1)

      goto 10

70    nhalf = nroots / 2
      if (1.gt.nhalf) return
      do  i=1,nhalf
         x = r(i)
         r(i) = r(nroots+1-i)
         r(nroots+1-i) = x
      end do

      return
      end
