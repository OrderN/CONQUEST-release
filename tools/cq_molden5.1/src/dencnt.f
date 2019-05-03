      subroutine dencnd(npts1,npts2,fcnt,dens,iedlog)
      implicit double precision (a-h,o-z)

c THIS IS REALLY dencnt

      parameter (mxcntr=100)
      logical euclid, yes, oscal, ostep, fine
      logical valenc,bonds,ovrlap,atomic,doori,dolap
      real xx, yy
      character str*100, esc
      common /plspec/ step,stepf,scale,scalf,cut,pmin,pmax,one,
     &                oscal,ostep,fine
      common /cntval/ cntval, euclid, yes, icnt, iflag
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap
      common /plot/   iplot,iplwin,icolps
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension fcnt(*),dens(*),iedlog(*)

      esc = char(27)

      if (iplot.eq.6) then
         xx = 2.0
         yy = 0.0
         call xwin(xx,yy,99,str,nstr,idum1,idum2)
      elseif (iplot.eq.4.and.icolps.eq.1) then
         write(iun4,*) 'poscontour setcol'
      else
         write(iun3,111) step
         write(iun3,'('' CONTOUR     CONTOURVALUE''//)')
      endif
111   format(//20x,'INTERVAL BETWEEN CONTOURS IS',f8.5,/,20x,
     +          '   ELECTRONS PER CUBIC A.U.',/)

      icnt = 0

c----- Positive Contours 

      is = 1
      if (bonds.or.ovrlap.or.atomic) is=2

      do i=is,100
          cntval = (i-1)*step
          if (cntval.lt.pmin*cut) goto 107
          if (cntval.gt.pmax*cut) goto 107
          iflag = 1
          if (i.eq.2) iflag = 0
          if (icnt.lt.mxcntr) then
             icnt = icnt + 1
             fcnt(icnt) = cntval
             if (iplot.ne.6) write(iun3,'(4x,i3,15x,f13.5)')
     +                    icnt,cntval
          endif
          height = cntval*scale
          if (iplot.eq.4.and.i.eq.1) write(iun4,*) 'nullcontour {'
          call cntour(dens,npts2,npts2,npts1,height,cntval,x,iedlog)
          if (iplot.eq.4.and.i.eq.1) write(iun4,*) '} if'
107       continue
      end do

c------ select other line type 

      if (iplot.eq.0) then
          write(iun4,'(a)') '.nc 3'
          write(iun4,'(a)') '.ls 2 0.2'
      endif
      if (iplot.eq.1) write(iun4,'(a)') 'SP3;'
      if (iplot.eq.2) write(iun4,'(a,''*m7b'')')esc
      if (iplot.eq.3) then
          idum = 1
          call plotgh(idum,0.0d0,0.0d0)
          write(iun4,*) esc//'a'
      endif
      if (iplot.eq.4) then
          write(iun4,'(''s'')')
          if (icolps.eq.1) then
              write(iun4,*) 'negcontour setcol'
          else
              write(iun4,'(''0 setlinecap [4 7] 0 setdash'')')
          endif
          write(iun4,'(''n'')')
      endif
      if (iplot.eq.6) then
          xx = 1.0
          yy = 0.0
          call xwin(xx,yy,99,str,nstr,idum1,idum2)
          call xwin(xx,yy,9,str,nstr,idum1,idum2)
          idum = 1
          call plotgh(idum,0.0d0,0.0d0)
      endif

c------ Negative Contours 

      do i=1,20
          cntval = -step*i
          if (cntval.lt.pmin*cut) goto 112
          if (cntval.gt.pmax*cut) goto 112
          iflag = 1
          if (i.eq.1) iflag = 0
          if (icnt.lt.mxcntr) then
              icnt = icnt + 1
              fcnt(icnt) = cntval
              if (iplot.ne.6) 
     +        write(iun3,'(4x,i3,15x,f13.5)') icnt,cntval
          endif
          height = cntval*scale
          call cntour(dens,npts2,npts2,npts1,height,cntval,x,iedlog)
112       continue
      end do

      if (iplot.eq.4.and.icolps.eq.1) then
          write(iun4,'(''s'')')
          write(iun4,'(''0 setgray'')')
          write(iun4,'(''n'')')
      endif

      cntval=99.999

      return
      end
