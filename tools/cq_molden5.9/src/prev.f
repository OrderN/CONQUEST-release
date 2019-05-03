      subroutine prev(v,m,n,ndim)
c
c     ----- print out matrices
c     vectors density etc.
c
      implicit double precision (a-h,o-z)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      dimension v(ndim,*)

      max  = 12
      imax = 0

100      imin = imax + 1
         imax = imax + max
         if (imax .gt. m) imax = m
         write (iun3,9008)
         write (iun3,8028) (i,i = imin,imax)
         write (iun3,9008)
         do j = 1,n
            write (iun3,8048) j,(v(j,i),i = imin,imax)
         end do
         if (imax .lt. m) go to 100

      return
 9008 format(/)
 8028 format(7x,12(3x,i3,3x))
 8048 format(i5,2x,12f9.5)

      end
