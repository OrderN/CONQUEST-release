      subroutine vsc1(a,scale,tol)
      implicit double precision (a-h,o-z)
      dimension a(3)

      rlen = vlen(a)

      if (rlen.gt.tol) then
         do i=1,3
            a(i) = a(i)*scale/rlen
         end do
      endif

      return
      end

      subroutine vnrm(a)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(3)

      vl = vlen(a)

      if (vl.gt.0.0d0) then
         do i=1,3
            a(i) = a(i)/vl
         end do
      endif

      return
      end

      subroutine vscal(a,n,scale)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(n)


      do i=1,n
         a(i) = a(i)*scale
      end do

      return
      end

      subroutine vcpy(a,b,n)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(n),b(n)

      do i=1,n
         a(i) = b(i)
      end do

      return
      end

      subroutine v3cpy1(a,b,n,iel)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(3,n),b(n)

      do i=1,n
         a(iel,i) = b(i)
      end do

      return
      end

      subroutine v3cpy2(a,b,n,iel)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(n),b(3,n)

      do i=1,n
         a(i) = b(iel,i)
      end do

      return
      end

      subroutine vadd(a,b,c,n)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(n),b(n),c(n)

      do i=1,n
         c(i) = a(i) + b(i)
      end do

      return
      end

      subroutine v3add1(a,b,c,n,iel)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(3,n),b(n),c(n)

      do i=1,n
         c(i) = a(iel,i) + b(i)
      end do

      return
      end

      subroutine v3add2(a,b,c,iel)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(3,*),b(3),c(3)

      do i=1,3
         c(i) = a(i,iel) + b(i)
      end do

      return
      end

      subroutine vsub(a,b,c,n)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(n),b(n),c(n)

      do i=1,n
         c(i) = a(i) - b(i)
      end do

      return
      end

      subroutine vseti(ia,n,ival)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension ia(n)

      do i=1,n
         ia(i) = ival
      end do

      return
      end

      subroutine vsetr(a,n,rval)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(n)

      do i=1,n
         a(i) = rval
      end do

      return
      end

