      character*5 function tk4014(ix,iy)
      integer highx,highy,lowx,lowy,extrax,extray,extra
c
c     tektronix-4014
c
c     highx  = <bit-7> ... <bit-11> von ix
c     highy  = <bit-7> ... <bit-11> von iy
c     lowx   = <bit-2> ... <bit-6>  von ix
c     lowy   = <bit-2> ... <bit-6>  von iy
c     extrax = <bit-0> ... <bit-1>  von ix
c     extray = <bit-0> ... <bit-1>  von iy
c
      highx=mod(ix/128,32)
      highy=mod(iy/128,32)
      lowx=mod(ix/4,32)
      lowy=mod(iy/4,32)
      extrax=mod(ix,4)
      extray=mod(iy,4)
c
c     extra  <bit-0> <bit-1>  extrax
c            <bit-2> <bit-3>  extray
c
      extra=4*extray+extrax
c
c
c     tag bits <bit-6> <bit-5>
c                 0       1      highx/y
c                 1       0      lowx
c                 1       1      lowy/extra
c
c
      highx=highx+32
c     set the 5. bit (beginning with bit 0)
      highy=highy+32
      lowx=lowx+64
c     set the 6. bit (beginning with bit 0)
      lowy=lowy+32+64
c     set the 5.+6. bit (beginning with bit 0)
      extra=extra+32+64
c     set the 5.+6. bit (beginning with bit 0)
c
      tk4014=char(highy)//char(extra)//char(lowy)//
     &        char(highx)//char(lowx)
      return
      end
