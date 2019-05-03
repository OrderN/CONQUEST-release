      double precision function reada(a,istart,ll)
      implicit double precision (a-h,o-z)
      character*1 a(ll)
      logical multi,scient

      nine  = ichar('9')
      izero = ichar('0')
      minus = ichar('-')
      iplus = ichar('+')
      idot  = ichar('.')
      ie    = ichar('E')
      iee   = ichar('e')
      id    = ichar('D')
      idd   = ichar('d')
      idig  = 0
      k1    = 0
      k2    = 0
      k3    = 0
      one   = 1.d0
      x     = 1.d0
      do j=istart,ll
         n = ichar(a(j))
         if((n.le.nine.and.n.ge.izero).or. n.eq.minus.or.n.eq.idot)
     &   goto 7
      end do
      reada = 0.d0
      return

7     continue
c     before the dot
      do i=j,ll
         n = ichar(a(i))
         if (n.le.nine.and.n.ge.izero) then
            idig = idig + 1
            if (idig.gt.10) then
               ii = i
               goto 90
            endif
            k1 = k1*10 + n - izero
         elseif (n.eq.minus.and.i.eq.j) then
             one = -1.d0
         elseif (n.eq.idot) then
             goto 30
         else
             ii = i
             goto 90
         endif
      end do
30    continue
c     after the dot
      idig = 0
      do ii=i+1,ll
         n = ichar(a(ii))
         if (n.le.nine.and.n.ge.izero) then
            idig = idig + 1
            if (idig.gt.9) goto 90
            k2 = k2*10 + n - izero
            x = x /10
         elseif (n.eq.minus.and.ii.eq.i) then
            x = -x
         else
            goto 90
         endif
      end do
90    continue
      scient = .false.
      do jj=ii,ll
         n = ichar(a(jj))
         if (n.eq.ie.or.n.eq.iee.or.n.eq.id.or.n.eq.idd) goto 95
      end do
      goto 100
95    continue
         if (jj.lt.159) then
            scient = .true.
            n = ichar(a(jj+1))
            if (n.eq.minus.or.n.eq.iplus) then
                multi = .true.
                if (n.eq.minus) multi = .false.
                idig = 0
                do j=jj+2,ll
                   n = ichar(a(j))
                   if (n.le.nine.and.n.ge.izero) then
                       idig = idig + 1
                       if (idig.gt.9) goto 90
                       k3 = k3*10 + n - izero
                   else
                      goto 100
                   endif
                end do
            endif
         endif

100   continue
      reada = one * ( k1 + k2 * x)
      if (scient) then
          if (multi) then
              reada = reada * 10**k3
          else
              if (k3.gt.25) then
                 reada = 0.0d0
              else
                 reada = reada / 10**k3
              endif
          endif
      endif

      return
      end
