      subroutine impsc(a,b,c)
      implicit double precision (a-h,o-z)
      dimension a(3),b(3)

      rimp = 0.0d0
    
      do i=1,3
         rimp = rimp + a(i)*b(i)
      end do

      al = vlen(a)
      bl = vlen(b)

      if (al.gt.0.0d0.and.bl.gt.0.0d0) then
         c = rimp/(vlen(a)*vlen(b))
      else
         c = 0.0d0
      endif

      return
      end

      subroutine timpsc(a,b,c)
      implicit double precision (a-h,o-z)
      dimension a(3),b(3)

      c = 0.0d0
    
      do i=1,3
         c = c + a(i)*b(i)
      end do

      return
      end
