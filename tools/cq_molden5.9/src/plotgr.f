      subroutine plotgr(ipen,dy,dx)
      implicit double precision (a-h,o-z)
      parameter (mxfill = 1000)
      logical euclid, yes
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*64 carray
      character*5 tk4014
      character*13 str
      real xx, yy
      integer*2 ixx
      logical oclose,osshad,ofill,ofrst
      common /spacom/ grey,iscol,oclose,osshad,ofill,ofrst
      common /fill /  ixx(mxfill),nixx,ihight,icol1
      common /plot / iplot,iplwin,icolps
      common /cntval/ cntval, euclid, yes, icnt, iflag
      common /tektro/ carray, ic
c      integer istmnt
      character*2 gstr

      if (euclid) then
         r2 = 1.0d0
      else
         r2 = dsqrt(2.0d0)
      endif

      dx = 0.5d0 + dx / r2
      dy = 0.5d0 + dy / r2

      isgrey = 0
      if (grey.lt.1.0) isgrey = 1

      if (iplot.eq.6) then
         xx = dx
         yy = dy
         if (ipen.eq.1) then
             if (oclose.and.ofill) then
                if (nixx.gt.1) then
                   nixx = nixx + 1
                   ixx(2*nixx-1) = ixx(1)
                   ixx(2*nixx)   = ixx(2)
                   icol1 = 10
                   if (isgrey.eq.1) icol1  = 6
                   call drwpol(ixx,nixx,icol1,0,isgrey,1)
                endif
                nixx = 1
                ixx(1) = dx*ihight
                ixx(2) = (1.0-dy)*ihight
             endif
             call xwin(xx,yy,2,str,nstr,inct,incp)
             if (yes.and.euclid) then
                 str = gstr(icnt)
                 call xwin(xx,yy,4,str,2,inct,incp)
             endif 
         else
             if (oclose.and.ofill) then
                 if (nixx.le.mxfill-2) then
                    nixx = nixx + 1
                    ixx(2*nixx-1) = dx*ihight
                    ixx(2*nixx)   = (1.0-dy)*ihight
                 endif
             endif
             call xwin(xx,yy,3,str,nstr,inct,incp)
         endif
         return
      endif

      dx = max(0.d0,min(1.d0,dx))
      dy = max(0.d0,min(1.d0,dy))

c      istmnt = ipen
c      if (ipen.eq.2) then
c        if (iplot.eq.0) assign 20 to istmnt
c        if (iplot.eq.1) assign 40 to istmnt
c        if (iplot.eq.2) assign 60 to istmnt
c        if (iplot.eq.4) assign 100 to istmnt
c      else
c        if (iplot.eq.0) assign 10 to istmnt
c        if (iplot.eq.1) assign 30 to istmnt
c        if (iplot.eq.2) assign 50 to istmnt
c        if (iplot.eq.4) assign 90 to istmnt
c      endif

      if (iplot.eq.0) then
          if (ipen.eq.2) then
             write(iun4,'(a,2f8.5)') '.d',dx,dy
          else
             write(iun4,'(a,2f8.5)') '.m',dx,dy
          endif
      endif

      if (iplot.eq.1) then
          if (ipen.eq.2) then
             write(iun4,'(a,2f8.4,a)') 'PD',dx,dy,';'
          else
             write(iun4,'(a,2f8.4,a)') 'PU',dx,dy,';'
          endif
      endif

      if (iplot.eq.2) then
          if (ipen.eq.2) then
              write(iun4,'(a,a,i3,a,i3,a)') 
     &           char(27),'*pb',int(dx*389),',',int(dy*389),'Z'
          else
              write(iun4,'(a,a,i3,a,i3,a)') 
     &           char(27),'*pa',int(dx*389),',',int(dy*389),'Z'
          endif
      endif

      if (iplot.eq.4) then
        if (ipen.ne.2) then
          if (oclose.and.ofill) then
              write(iun4,*) 'closepath'
              write(iun4,*) 'gsave'
              if (icolps.eq.1) then
                 if (isgrey.eq.1) then
                    write(iun4,*) 'negfill setcol'
                 else
                    write(iun4,*) 'posfill setcol'
                 endif
              else
                 write(iun4,*)  grey,' setgray'
              endif
              write(iun4,*) 'fill'
              write(iun4,*) 'grestore'
              if (icolps.eq.0) then
                 write(iun4,*) '0 setgray'
              endif
          endif

          write(iun4,'(''s'')')

          if (euclid) then
             if (yes) then
                write(iun4,'(''tabel {'',i4,'' '',i4,'' m ('',
     &                     i2,'') show } if'')')
     &           int(dx*2000),int(dy*2000+125),icnt
             else
                if (cntval.ne.99.999)
     &           write(iun4,'(''contourvalues {'',i4,'' '',i4,'' m ('',
     &                     f10.6,'') show } if'')')
     &           int(dx*2000),int(dy*2000+125),cntval
             endif
          endif

          write(iun4,'(''n'')')

        endif

        if (ipen.eq.2) then
          write(iun4,'(a,i4,a,i4,a)') 
     &       ' ',int(dx*2000),' ',int(dy*2000)+125,' l'
        else
          write(iun4,'(a,i4,a,i4,a)') 
     &       ' ',int(dx*2000),' ',int(dy*2000)+125,' m'
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
c10    format('.m',f8.5,f8.5)  
c20    format('.d',f8.5,f8.5)
c30    format('PU',f8.4,f8.4,';')  
c40    format('PD',f8.4,f8.4,';')
c50    format(a,'*pa',i3,',',i3,'Z')
c60    format(a,'*pb',i3,',',i3,'Z')
c90    format(' ',i4,' ',i4,' m')
c100    format(' ',i4,' ',i4,' l')

      if (yes.and.euclid) then
        if (iplot.eq.0) then
           write(iun4,'(''.to 3'')')
           write(iun4,'(''.ch 1.2'')')
           write(iun4,'(''.pt '',i2)') icnt
           write(iun4,'(''.to 5'')')
           write(iun4,'(''.ch 1.5'')')
        endif
        if (iplot.eq.1) then
           write(iun4,'(''LB'',i2,a,'';'')')icnt,char(3)
        endif
      endif

      return
      end
