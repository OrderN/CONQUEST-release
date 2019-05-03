      subroutine plbox
      implicit double precision (a-h,o-z)

      common /plot/   iplot,iplwin,icolps
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
*
*  DRAW A BOX AROUND THE PLOT
*
      if (iplot.eq.0) write(iun4,'(a)') '.nc 1'
      if (iplot.eq.1) write(iun4,'(a)') 'SP1;'

      call euler(0.d0,0.d0,0.d0,1)
      call euler(1.d0,0.d0,0.d0,2)
      call euler(1.d0,1.d0,0.d0,2)
      call euler(0.d0,1.d0,0.d0,2)
      call euler(0.d0,0.d0,0.d0,2)

      if (iplot.eq.4) write(iun4,'(''s'')')

      return
      end
