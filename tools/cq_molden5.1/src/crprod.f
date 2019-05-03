      subroutine crprod(a,b,c)
c
c       calculates cross product:   a x b = c
c                                   -   -   -
      implicit double precision (a-h,o-z)
      dimension a(3),b(3),c(3)
c
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)

      return
      end
      subroutine crpros(a,b,c)
c
c       calculates cross product:   a x b = c
c                                   -   -   -
      implicit double precision (a-h,o-z)
      real b
      dimension a(3),b(3),c(3)
c
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)

      return
      end

      subroutine arbper(a,b)
      implicit double precision (a-h,o-z)
      dimension a(3),b(3)

      a1 = dabs(a(1))
      a2 = dabs(a(2))
      a3 = dabs(a(3))

      do i=1,3
         b(i) = 0.0d0
      end do

      if (a1.lt.a2) then
         if (a1.lt.a3) then
            b(1) = 1.0d0
         else
            b(3) = 1.0d0
         endif
      else
         if (a2.lt.a3) then
            b(2) = 1.0d0
         else
            b(3) = 1.0d0
         endif
      endif

      return
      end

