      subroutine vclr(a,incr,n)
      implicit double precision (a-h,o-z)
      dimension a(*)

      if(incr.eq.1) then
        do loop=1,n
          a(loop) = 0.0d0
        end do
      else
        loopi=1
        do loop=1,n
           a(loopi)=0.0d0
           loopi=loopi+incr
        end do
      endif

      return
      end
