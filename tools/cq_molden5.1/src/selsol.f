      subroutine selsol
      implicit double precision (a-h,o-z)
      character str*100, esc
      real xx, yy
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /plot/   iplot,iplwin,icolps

c----- select solid line again 

      esc    = char(27)

      if (iplot.eq.0) write(iun4,'(a)') '.ls 1 0'

      if (iplot.eq.1) write(iun4,'(a)') 'SP1;'

      if (iplot.eq.2) write(iun4,'(a,''*m1b'')')esc

      if (iplot.eq.4) then
         write(iun4,'(''s'')')
         write(iun4,'(''[] 0 setdash'')')
         write(iun4,'(''   1 setlinewidth'')')
         write(iun4,'(''n'')')
      endif

      if (iplot.eq.3) then
         idum=1
         call plotgh(idum,0.0d0,0.0d0)
         write(iun4,*) esc//'`'
      endif

      if (iplot.eq.6) then
         xx = 15.0
         yy = 0.0
         call xwin(xx,yy,99,str,nstr,idum1,idum2)
         call xwin(xx,yy,8,str,nstr,idum1,idum2)
         idum = 1
         call plotgh(idum,0.0d0,0.0d0)
      endif

      return
      end
