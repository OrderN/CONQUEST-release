      subroutine eulstr(px,py,pz,rval)
      implicit double precision (a-h,o-z)
c     print string
      real xx, yy
      character*80 str
      logical euclid, yes
      common /cntval/ cntval, euclid, yes, icnt, iflag
      common /eul/    ca,cb,sa,sb,ccab,scab,ssab,scba
      common /plrat/  rat(3)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /plot / iplot,iplwin,icolps

      roff = 0.01d0
      if (rval.lt.10000.0d0.and.rval.ge.1.0d0) then
          write(str,'(f7.1)') rval
          nstr = 7
      else
          write(str,'(f7.4)') rval
          nstr = 7
      endif
      cy = ccab*(px-0.5d0)*rat(1) + scab*(py-0.5d0)*rat(2) 
     &     - sb*pz*rat(3)
      cx =  -sa*(px-0.5d0)*rat(1) +   ca*(py-0.5d0)*rat(2)

      if (euclid) then
         r2 = 1.0d0
      else
         r2 = dsqrt(2.0d0)
      endif
      cx = 0.5d0 + cx / r2
      cy = 0.5d0 + cy / r2

      if (iplot.eq.6) then
         if (rval.le.0.d0) then
            xx = 5.0
         else
            xx = 7.0
         endif
         call xwin(xx,yy,99,str,nstr,idum1,idum2)
         xx = cx - roff
         yy = cy
         call xwin(xx,yy,2,str,nstr,idum1,idum2)
         xx = cx + roff
         yy = cy
         call xwin(xx,yy,3,str,nstr,idum1,idum2)
         xx = cx 
         yy = cy - roff
         call xwin(xx,yy,2,str,nstr,idum1,idum2)
         xx = cx
         yy = cy + roff
         call xwin(xx,yy,3,str,nstr,idum1,idum2)
         xx = 5.0
         call xwin(xx,yy,99,str,nstr,idum1,idum2)
         xx = cx
         yy = cy 
         call xwin(xx,yy,4,str,nstr,idum1,idum2)
      endif

      if (iplot.eq.4) then
         write(iun4,'(''   0.5 setgray'')')
         write(iun4,'(''n '',i4,'' '',i4,'' 10 0 360 arc fill'')')
     +   int(cx*2000),int(cy*2000+125)
         write(iun4,'(''   0 setgray'')')
         write(iun4,'(i4,'' '',i4,'' m (   '',a,'') show'')')
     +   int(cx*2000),int(cy*2000+125),str(1:nstr)
      endif

      return
      end
