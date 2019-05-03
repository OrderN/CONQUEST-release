      subroutine pltab(fcnt)
      implicit double precision (a-h,o-z)
      logical  euclid, yes, oscal,ostep,fine
      real xx, yy
      character str*100, esc, etx, lf, gs, us
      character tk4014*5

      common /plot/   iplot,iplwin,icolps
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /cntval/ cntval, euclid, yes, icnt, iflag
      common /plspec/ step,stepf,scale,scalf,cut,pmin,pmax,one,
     &                oscal,ostep,fine
      dimension fcnt(*)

      etx    = char(3)
      lf     = char(10)
      esc    = char(27)
      gs     = char(29)
      us     = char(31)
                                                                                
      if (iplot.eq.0) then
        write(iun4,'(a)') '.m 1.25 1.0' 
        write(iun4,'(''.pt NUMBER OF CONTOURS ='',i4,a)') icnt
        write(iun4,'(a)') '.m 1.15 1.0' 
        write(iun4,'(''.pt CONTOUR INTERVAL ='',f13.6)') step
        write(iun4,'(a)') '.ch 1.2'
        write(iun4,'(a)') '.td 0.0' 
        write(iun4,'(a)') '.m 1.05 0.58' 
        write(iun4,'(''.pt  CONTOUR    VALUE'')')
      endif

      if (iplot.eq.1) then
         write(iun4,'(''PU46,18;LBNUMBER OF CONTOURS ='',i4,a,'';'')')
     +   icnt,etx
         write(iun4,'(''PU44,18;LBCONTOUR INTERVAL ='',f13.6,a,'';'')') 
     +   step,etx
         write(iun4,'(a)') 'SC0,1,0,1;DI;SL0;' 
         write(iun4,'(''PU1.05,0.58;LB CONTOUR    VALUE'',a,'';'')')
     +   etx
      endif

      if (iplot.eq.2) then
         write(iun4,'(a,''*pa408,222Z'')')esc
         write(iun4,'(a,''*lCONTOUR  VALUE'',a)')esc,lf
      endif

      if (iplot.eq.3) then
c         write(iun4,*)esc//'/0d'
         write(iun4,*)gs//tk4014(int(1.05*3100),int(0.58*3100))
     +              //us//'CONTOUR VALUE'
      endif

      if (iplot.eq.4) then
         write(iun4,'(''tabel {'')')
         write(iun4,'(''   1 setlinewidth'')')
         write(iun4,'(''2150 1285 m'')')
         write(iun4,'(''( CONTOUR   VALUE) show'')')
      endif

      if (iplot.eq.6) then
        xx= 15.0
        yy= 0.0
        call xwin(xx,yy,99,str,nstr,idum1,idum2)
        xx= 1.02
        yy= 0.85
        call xwin(xx,yy,4,'CONTOUR VALUE',13,idum1,idum2)
      endif

c-------------- plot contour value tabel 

        do i=1,icnt
         fpos = 0.55d0 - 0.02d0*dble(i)
         if (i.lt.26) then

           if (iplot.eq.0) then
              write(iun4,'(''.m 1.05 '',f8.5)') fpos
              write(iun4,'(''.pt '',i3,1x,f13.5)')i,fcnt(i)
           endif

           if (iplot.eq.1) then
              write(iun4,'(''PU1.05,'',f8.5,'';'')') fpos
              write(iun4,'(''LB'',i3,1x,f13.5,a,'';'')')i,fcnt(i),etx
           endif

           if (iplot.eq.2) then
              write(iun4,'(a,''*pa408,'',i3,''Z'')')esc,int(fpos*389)
              write(iun4,'(a,''*l'',i3,1x,f9.5,a)')esc,i,fcnt(i),lf
           endif

           if (iplot.eq.3) then
c              write(iun4,*)esc//'/0d'
              write(iun4,'(a,i3,1x,f9.5)')gs//tk4014(int(1.05*3100),
     +                 int(fpos*3100))//us,i,fcnt(i)
           endif

           if (iplot.eq.4) then
              write(iun4,'(''2150 '',i4, '' m'')')int(fpos*2000+125)
              write(iun4,'(''('',i3,10x,f13.5,'') show'')')i,
     +               fcnt(i)
           endif

           if (iplot.eq.6) then
             write(str,'(i3,1x,f9.5)') i,fcnt(i)
             if (fcnt(i) .ge. 0.0) then
               xx = 2.0
             else
               xx = 1.0
             endif
             yy = 0.0
             call xwin(xx,yy,99,str,nstr,idum1,idum2)
             xx = 1.02
             fpos = 0.80d0 - 0.03d0*dble(i)
             yy = fpos
             call xwin(xx,yy,4,str,13,idum1,idum2)
           endif

         endif
         if (i.eq.26) then

           if (iplot.eq.0) then
              write(iun4,'(''.m 1.05 '',f8.5)') fpos
              write(iun4,'(''.pt ............'')')
           endif

           if (iplot.eq.1) then
              write(iun4,'(''PU1.05,'',f8.5,'';'')') fpos
              write(iun4,'(''LB ............'',a,'';'')')etx
           endif

           if (iplot.eq.2) then
              write(iun4,'(a,''*pa408,'',i3,''Z'')')esc,int(fpos*389)
              write(iun4,'(a,''*l ............'',a)')esc,lf
           endif

           if (iplot.eq.3) then
c              write(iun4,*)esc//'/0d'
              write(iun4,*)gs//tk4014(int(1.05*3100),
     +                 int(fpos*3100))//us//'  ............'
           endif

           if (iplot.eq.4) then
              write(iun4,'(''2150 '',i4, '' m'')')int(fpos*2000+125)
              write(iun4,'(''(   ............) show'')')
           endif

           if (iplot.eq.6) then
             xx = 15.0
             yy = 0.0
             call xwin(xx,yy,99,str,nstr,idum1,idum2)
             xx = 1.02
             fpos = 0.80d0 - 0.03d0*dble(i)
             yy = fpos
             call xwin(xx,yy,4,'  ...........',13,idum1,idum2)
           endif

         endif
      end do

      if (iplot.eq.0) then
        write(iun4,'(a)') '.m 1.05 0.60' 
        write(iun4,'(a)') '.d 1.30 0.60' 
        write(iun4,'(a)') '.d 1.30 0.00' 
        write(iun4,'(a)') '.d 1.05 0.00' 
        write(iun4,'(a)') '.d 1.05 0.60' 
      endif

      if (iplot.eq.1) then
        write(iun4,'(''PU1.05,0.60;'')') 
        write(iun4,'(''PD1.30,0.60;'')') 
        write(iun4,'(''PD1.30,0.00;'')') 
        write(iun4,'(''PD1.05,0.00;'')') 
        write(iun4,'(''PD1.05,0.60;'')') 
      endif

      if (iplot.eq.2) then
      write(iun4,'(a,''*pa408,233,506,233,506,0,408,0,408,233Z'')')
     +        esc
      endif

      if (iplot.eq.4) then
        write(iun4,'(''n'')')
        write(iun4,'(''2100 1325 m'')')
        write(iun4,'(''2600 1325 l'')')
        write(iun4,'(''2600 125 l'')')
        write(iun4,'(''2100 125 l'')')
        write(iun4,'(''2100 1325 l'')')
        write(iun4,'(''s'')')
        write(iun4,'(''} if'')')
      endif

      return
      end
