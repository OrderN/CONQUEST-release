      subroutine shsort(n,arr,iarr)
      implicit double precision (a-h,o-z)
      parameter (aln2i = 1.0d0 / 0.69314718d0, tiny=1.d-6)
      dimension arr(*),iarr(*)
      if (n.eq.0) return
      lognb2 = int(dlog(dfloat(n))*aln2i+tiny)
      do i=1,n
         iarr(i) = i
      end do
      m = n
      do nn=1,lognb2
         m = m / 2
         k = n - m
         do j=1,k
            i = j
3           continue
            l = i + m
            if (arr(iarr(l)).lt.arr(iarr(i))) then
               it = iarr(i)
               iarr(i) = iarr(l)
               iarr(l) = it
               i = i - m
               if (i.ge.1) goto 3
            endif
         end do
      end do
       
      return
      end
