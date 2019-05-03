      subroutine plotgh(ipen,ex,ey)
      implicit double precision (a-h,o-z)
      character*64 carray
      character*5 tk4014
      character*13 str
      real xx, yy
      common /plot / iplot,iplwin,icolps
      common /tektro/ carray, ic
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
************************************************************************
*
* PLOTGH IS THE HIDE PLOTTING SUBROUTINE, AND SHOULD BE MODIFIED
*     t  TO SUIT THE LOCAL CONDITIONS.
* CONVENTIONS USED HERE ARE:
*   IPEN   = 1   PEN UP (NO LINE DRAWN, BUT PEN MOVES TO 
*                THE POSITION EX, EY)
*   IPEN   = 2   PEN DOWN (A LINE IS DRAWN FROM THE LAST POSITION TO 
*                THE POINT EX, EY)
*   EX, EY       'X' AND 'Y' COORDINATES. PEN WILL MOVE TO THIS POINT.
*
************************************************************************


      r2 = dsqrt(2.0d0)
      dx = 0.5d0 + ex/r2
      dy = 0.5d0 + ey/r2
      if (iplot.eq.6) then
         xx = dx
         yy = dy
         if (ipen.eq.1) then
             call xwin(xx,yy,2,str,nstr,inct,incp)
         else
             call xwin(xx,yy,3,str,nstr,inct,incp)
         endif
         ipen = 3 - ipen
         return
      endif
      dx = max(0.d0,min(1.d0,dx))
      dy = max(0.d0,min(1.d0,dy))
      if (ipen.eq.1) then
        if (iplot.eq.0) write(iun4,1)dx,dy
        if (iplot.eq.1) write(iun4,30)dx,dy
        if (iplot.eq.2) write(iun4,50)char(27),int(389*dx)
     +                             ,int(389*dy)
        if (iplot.eq.4) then
           write(iun4,'(''s'')')
           write(iun4,'(''n'')')
           write(iun4,90)int(dx*2000),int(dy*2000)+125
        endif
      else
        if (iplot.eq.0) write(iun4,2)dx,dy
        if (iplot.eq.1) write(iun4,40)dx,dy
        if (iplot.eq.2) write(iun4,60)char(27),int(389*dx)
     +                             ,int(389*dy)
        if (iplot.eq.4) then
          write(iun4,100)int(dx*2000),int(dy*2000)+125
        endif
 
      endif
      if (iplot.eq.3) then
         if (ipen.eq.1) then
c-- voer laaste line series uit
           if (carray(1:1).eq.char(29)) then
            carray(7+ic*5:7+ic*5) = char(13)
            write(iun4,'(a)') carray(1:7+ic*5)
           endif
           carray(1:1) = char(29)
           carray(2:6) = tk4014(int(dx*3100),int(dy*3100))
           ic = 0
         else
           if (ic.eq.9) then
             carray(7+ic*5:7+ic*5) = char(13)
             write(iun4,'(a)') carray(1:7+ic*5)
             carray(2:6)  = carray(2+ic*5:6+ic*5)
             carray(7:11) = tk4014(int(dx*3100),int(dy*3100))
             ic = 1
           else
             ic = ic + 1
             carray(2+ic*5:6+ic*5) = tk4014(int(dx*3100),
     &       int(dy*3100))
           endif
         endif
      endif
1     format('.m',f8.5,f8.5)  
2     format('.d',f8.5,f8.5)
30    format(' PU',f8.4,f8.4,';')  
40    format(' PD',f8.4,f8.4,';')
50    format(a,'*pa',i3,',',i3,'Z')
60    format(a,'*pb',i3,',',i3,'Z')
90    format(' ',i4,' ',i4,' m')
100    format(' ',i4,' ',i4,' l')
      ipen = 3 - ipen

      return
      end
