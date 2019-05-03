      subroutine plhead(title)
      implicit double precision (a-h,o-z)
      real xx, yy
      character str*100, esc, etx, lf, gs, us
      character keywrd*320, keyori*320, title*80
      character tk4014*5

      common /plot/   iplot,iplwin,icolps
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /keywrd/ keywrd,keyori

      etx    = char(3)
      lf     = char(10)
      esc    = char(27)
      gs     = char(29)
      us     = char(31)

c------- Plot Header

      if (iplot.eq.0) then
         write(iun4,'(a)') '.nc 1'
         write(iun4,'(a)') '.to 3'
         write(iun4,'(a)') '.ch 2.0'
         write(iun4,'(a)') '.fo 1' 
         write(iun4,'(a)') '.td 270.0' 
         write(iun4,'(a)') '.m 1.47 1.0' 
         write(iun4,'(''.pt '',a)')'********* M O L D E N  **********'
         write(iun4,'(a)') '.m 1.42 1.0' 
         write(iun4,'(''.pt '',a)')title
         write(iun4,'(a)') '.m 1.37 1.0' 
         write(iun4,'(a)') '.ch 1.5'
         write(iun4,'(''.pt '',a)')keywrd(1:80)
         write(iun4,'(a)') '.m 1.32 1.0' 
         write(iun4,'(''.pt '',a)')keywrd(81:160)
      endif

      if (iplot.eq.1) then
        write(iun4,'(a)') 'SC0,36,0,18;SP1;DI0,-1;' 
        write(iun4,'(''PU51,18;LB********* M O L D E N  **********''
     +          ,a,'';'')')etx
        write(iun4,'(''SL0.6;PU50,18;LB'',a,a,'';'')')title,etx
        write(iun4,'(''PU49,18;LB'',a,a,'';'')')keywrd(1:80),etx
        write(iun4,'(''PU48,18;LB'',a,a,'';'')')keywrd(81:160),etx
      endif

      if (iplot.eq.2) then
c       text vertical no slanting
        write(iun4,'(a,''*m1nO'')')esc
c       position for text write
        write(iun4,'(a,''*pa408,379Z'')')esc
c       write text
        write(iun4,'(a,''*lM O L D E N '',a)') esc,lf
        write(iun4,'(a,''*l'',a)')esc,lf
c       text slanting
        write(iun4,'(a,''*mO'')')esc
        write(iun4,'(a,''*l'',a,a)')esc,title(1:15),lf
c       text horizontal
        write(iun4,'(a,''*m1nO'')')esc
      endif

      if (iplot.eq.3) then
c------ flush last instructions tek4014 ----------------
       idum=1
       call plotgh(idum,0.0d0,0.0d0)
c       write(iun4,*)esc//'/0d'
       write(iun4,*)gs//tk4014(int(1.05*3100),int(0.95*3100))
     +            //us//'M O L D E N'
      endif

      if (iplot.eq.6) then
        xx= 15.0
        yy= 0.0
        call xwin(xx,yy,99,str,nstr,idum1,idum2)
        xx= 1.02
        yy= 0.95
c        call xwin(xx,yy,4,'M O L D E N',11,idum1,idum2)
      endif

      return
      end
