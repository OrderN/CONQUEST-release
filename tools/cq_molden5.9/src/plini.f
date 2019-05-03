      subroutine plini(iplot,fine,icolps)
      implicit double precision (a-h,o-z)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character esc, ff
      character*64 carray
      character line*80
      logical fine
      common /tektro/ carray, icn

      esc = char(27)
      ff = char(12)
c
c--------- plot initialisation ----------------
c
      if(iplot.eq.0)then
         write(iun4,'(a)') '.vp 0 150 0 100'
         write(iun4,'(a)') '.wn 0 1.5 0 1'
      endif
      if(iplot.eq.1)then
        line(1:1)=esc
        line(2:3)='.('
        line(4:4)=esc
        line(5:13)='.I81;;17:'
        line(14:14)=esc
        line(15:20)='.N;19:'
        write(iun4,'(a)')line(1:20)
        write(iun4,'(a)')';IN;IP200,430,7400,7630;PA;' 
        write(iun4,'(a)') 'SC0,1,0,1;SL0;SR0.8,1.5;SP1;DI;'
      endif
      if(iplot.eq.2)then
c       turn graphics on
        write(iun4,'(a,''*dC'')')esc
c       turn alfa off
        write(iun4,'(a,''*dF'')')esc
c       clear graphics
        write(iun4,'(a,''*dA'')')esc
c       select solid line type
        write(iun4,'(a,''*m1b'')')esc
      endif

      if(iplot.eq.3)then
        icn=0
        carray(1:1)=' '
        write(iun4,*) esc//ff
        write(iun4,*) esc//'8'//esc//'`'
      endif

      if(iplot.eq.4.or.iplot.eq.6)then
       write(iun4,'(''%!PS-Adobe-2.0 EPSF-2.0'')')
       write(iun4,'(''%%Title: Molden'')')
       write(iun4,'(''%%For: Schaft'')')
       write(iun4,'(''%%Creator: Drs G Schaftenaar'')')
       write(iun4,'(''%%DocumentFonts: Courier'')')
       write(iun4,'(''%%Pages (atend)'')')
       write(iun4,'(''%%BoundingBox: 0 0 612 792'')')
       write(iun4,'(''%%EndComments'')')
       write(iun4,'(''%'')')
       write(iun4,'(''%###### User Preferences ############'')')
       write(iun4,'(''%'')')
       write(iun4,'(''%---- SIZE AND ORIENTATION OF THE PLOT ---'')')
       write(iun4,'(''%'')')
       write(iun4,'(''/size    {  0.24 } def'')')
       write(iun4,'(''%---- These number can be negative -------'')')
       write(iun4,'(''/originx {  39.0 } def'')')
       write(iun4,'(''/originy { 753.0 } def'')')
       write(iun4,'(''/angle   { -90.0 } def'')')
       write(iun4,'(''%For Portrait use'')')
       write(iun4,'(''%/originx { 40.0 } def'')')
       write(iun4,'(''%/originy { 240.0 } def'')')
       write(iun4,'(''%/angle   { 0.0 } def'')')
       write(iun4,'(''%and BoundingBox: 25 255 535 765'')')
       write(iun4,'(''%'')')
       write(iun4,'(''%---- COLORS FOR DENSITY CONTOURS  -------'')')
       write(iun4,'(''%'')')
       write(iun4,'(''/poscontour { 18 } def'')')
       write(iun4,'(''/negcontour { 19 } def'')')
       write(iun4,'(''%---- FILL COLORS OF DENSITY SPACE MODE  -'')')
       write(iun4,'(''/posfill { 16 } def'')')
       write(iun4,'(''/negfill { 17 } def'')')
       write(iun4,'(''%---- COLOR OF UNIT CELL DEF. BLACK -'')')
       write(iun4,'(''/cellcol { 0 setcol } def'')')
       write(iun4,'(''/labelcol { 0 setcol } def'')')
       write(iun4,'(''%/cellcol { 1 setcol } def'')')
       write(iun4,'(''%'')')
       write(iun4,'(''%---- COLORS HAVE A HUE AND SATURATION ---'')')
       write(iun4,'(''%'')')
       write(iun4,'(''/hues {[.0 .0 .17 .33 .66 .5 .12 .54 .0 '',
     + ''.83 .33 .1 .08 .15 .10 .0 .33 .1 .66 .0]} def'')')
       write(iun4,'(''/satus {[.0 .93 .95 1. .5 1. 1. .38 .0 '',
     + ''1. .9 .41 1. .6 .95 .0 .5 .7 1. 1.]} def'')')
       write(iun4,'(''%'')')
       write(iun4,'(''%---- SET BOND RENDERING:  ---------------'')')
       write(iun4,'(''%---- shadedrod, whiterod, blackrod  -----'')')
       write(iun4,'(''%'')')
       write(iun4,'(''/dorod { shadedrod } def'')')
       write(iun4,'(''%'')')
       write(iun4,'(''%---- Include Tabel & Logo, Fontsize -----'')')
       write(iun4,'(''/tabel {true} def'')')
       write(iun4,'(''/titleandlogo {true} def'')')
       write(iun4,'(''/nullcontour {false} def'')')
       write(iun4,'(''/contourvalues {false} def'')')
       write(iun4,'(''/stickwidth {8} def'')')
       write(iun4,'(''/fontwidth {27} def'')')
       write(iun4,'(''/fontheight {35} def'')')
       write(iun4,'(''/dobackground {false} def'')')
       write(iun4,'(''%'')')
       write(iun4,'(''%###### END User Preferences ########'')')
       write(iun4,'(''/hue {hues col get} def'')')
       write(iun4,'(''/satu {satus col get} def'')')
c       write(iun4,
c     +'(''/setcol {/col exch def hue satu 1.0 sethsbcolor} def'')')
       write(iun4,'(''/setcol {'')')
       write(iun4,'('' /col exch def col 0 eq'')')
       write(iun4,'('' {hue satu 0.0 sethsbcolor}'')')
       write(iun4,'('' {hue satu 1.0 sethsbcolor}'')')
       write(iun4,'(''ifelse} def'')')
       write(iun4,'(''/m { moveto } def'')')
       write(iun4,'(''/l { lineto } def'')')
       write(iun4,'(''/s { stroke } def'')')
       write(iun4,'(''/n { newpath } def'')')
       write(iun4,'(''/lc { setlinecap } def'')')
       write(iun4,'(''/offset { 0 0 moveto (A) false charpath'')')
       write(iun4,'(''flattenpath pathbbox'')')
       write(iun4,'(''pop pop pop -1 mul /xoff exch def } def'')')
       write(iun4,'(''/doatom'')')
       write(iun4,'(''{ gsave'')')
       write(iun4,'(''  rx ry translate'')')
       if (fine) then
           write(iun4,'(''  90 -1 1'')')
       else
           write(iun4,'(''  90 -5 1'')')
       endif
       write(iun4,'(''  { gsave'')')
       if (icolps.eq.1) then
          write(iun4,
     &'(''    dup cos hue exch satu exch sethsbcolor sin dup scale'')')
       else
          write(iun4,'(''    dup cos setgray sin dup scale'')')
       endif
       write(iun4,'(''    newpath'')')
       write(iun4,'(''    0 0 rad 0 360 arc'')')
       write(iun4,'(''    closepath fill grestore } for'')')
       write(iun4,'(''    grestore } def'')')
       write(iun4,'(''/shadedrod'')')
       write(iun4,'(''{ gsave'')')
       write(iun4,'(''  x1 y1 translate'')')
       write(iun4,'(''  x2 x1 neg add'')')
       write(iun4,'(''  y2 y1 neg add'')')
       write(iun4,'(''  {atan neg rotate} stopped not {'')')
       if (fine) then
          write(iun4,'(''  87 -3 0'')')
       else
          write(iun4,'(''  85 -5 0'')')
       endif
       write(iun4,'(''  {dup'')')
       write(iun4,'(''  gsave'')')
       write(iun4,'(''  newpath'')')
       if (icolps.eq.1) then
          write(iun4,'(''   cos 1.0 cosb 0.5 mul neg add mul'')')
          write(iun4,'(''   hue exch satu exch sethsbcolor'')')
       else
          write(iun4,
     &         '(''   cos 1.0 cosb 0.5 mul neg add mul setgray'')')
       endif
       write(iun4,'(''   sin 1.0 scale'')')
       write(iun4,'(''   1 cosb scale'')')
       write(iun4,'(''   0 0 hd 0 180 arcn'')')
       write(iun4,'(''   x2 x1 neg add dup mul'')')
       write(iun4,'(''   y2 y1 neg add dup mul'')')
       write(iun4,'(''   add sqrt'')')
       write(iun4,
     & '(''  0 cosb eq {/cosb 1.0 def} if 0 exch cosb div translate'')')
       write(iun4,'(''   0 0 hd 180 360 arc'')')
       write(iun4,'(''  closepath fill'')')
       write(iun4,'(''  grestore } for'')')
       write(iun4,'(''  } if'')')
       write(iun4,'(''  grestore } def'')')
       write(iun4,'(''/blackrod'')')
       write(iun4,'(''{ gsave'')')
       write(iun4,'(''  x1 y1 translate'')')
       write(iun4,'(''  x2 x1 neg add'')')
       write(iun4,'(''  y2 y1 neg add'')')
       write(iun4,'(''  {atan neg rotate} stopped not {'')')
       write(iun4,'(''  newpath'')')
       write(iun4,'(''   0 setgray'')')
       write(iun4,'(''   1 cosb scale'')')
       write(iun4,'(''   0 0 hd 0 180 arcn'')')
       write(iun4,'(''   x2 x1 neg add dup mul'')')
       write(iun4,'(''   y2 y1 neg add dup mul'')')
       write(iun4,'(''   add sqrt'')')
       write(iun4,
     & '(''  0 cosb eq {/cosb 1.0 def} if 0 exch cosb div translate'')')
       write(iun4,'(''   0 0 hd 180 360 arc'')')
       write(iun4,'(''  closepath stroke'')')
       write(iun4,'(''  } if'')')
       write(iun4,'(''  grestore } def'')')
       write(iun4,'(''/whiterod'')')
       write(iun4,'(''{ gsave'')')
       write(iun4,'(''  x1 y1 translate'')')
       write(iun4,'(''  x2 x1 neg add'')')
       write(iun4,'(''  y2 y1 neg add'')')
       write(iun4,'(''  {atan neg rotate} stopped not {'')')
       write(iun4,'(''  newpath'')')
       write(iun4,'(''   1 setgray'')')
       write(iun4,'(''   1 cosb scale'')')
       write(iun4,'(''   0 0 hd 0 180 arcn'')')
       write(iun4,'(''   x2 x1 neg add dup mul'')')
       write(iun4,'(''   y2 y1 neg add dup mul'')')
       write(iun4,'(''   add sqrt'')')
       write(iun4,
     & '(''  0 cosb eq {/cosb 1.0 def} if 0 exch cosb div translate'')')
       write(iun4,'(''   0 0 hd 180 360 arc'')')
       write(iun4,'(''  closepath stroke'')')
       write(iun4,'(''  } if'')')
       write(iun4,'(''  grestore } def'')')
       write(iun4,'(''/dobond'')')
       write(iun4,'(''{ gsave'')')
       write(iun4,'(''  0 setlinecap'')')
       write(iun4,'(''  newpath moveto lineto'')')
       write(iun4,'(''  90 -5 5'')')
       write(iun4,
     & '(''  { gsave dup cos hue exch satu exch sethsbcolor'')')
       write(iun4,'(''    sin 23 mul setlinewidth'')')
       write(iun4,'(''    stroke grestore } for'')')
       write(iun4,'(''  grestore } def'')')
       write(iun4,'(''/dobond2'')')
       write(iun4,'(''{ gsave'')')
       write(iun4,'(''  1 setlinecap'')')
       write(iun4,'(''  2 setlinejoin'')')
       write(iun4,'(''  newpath moveto lineto lineto'')')
       write(iun4,'(''  90 -5 5'')')
       write(iun4,
     & '(''  { gsave dup cos hue exch satu exch sethsbcolor'')')
       write(iun4,'(''    sin 23 mul setlinewidth'')')
       write(iun4,'(''    stroke grestore } for'')')
       write(iun4,'(''  grestore } def'')')
       write(iun4,'(''/doshadedstick'')')
       write(iun4,'(''{ gsave'')')
       write(iun4,'(''  stickwidth setlinewidth'')')
       write(iun4,'(''  1 setlinecap'')')
       write(iun4,'(''  2 setlinejoin'')')
       if (icolps.eq.1) then
          write(iun4,'(''   1.0 cosb 0.8 mul neg add '')')
          write(iun4,'(''   hue exch satu exch sethsbcolor'')')
       else
          write(iun4,
     &         '(''   1.0 cosb 0.8 mul neg add setgray'')')
       endif
       write(iun4,'(''  newpath'')')
       write(iun4,'(''  x1 y1 moveto'')')
       write(iun4,'(''  x2 y2 lineto'')')
       write(iun4,'(''  stroke'')')
       write(iun4,'(''  grestore } def'')')
       write(iun4,'(''/dostick'')')
       write(iun4,'(''{ gsave'')')
       write(iun4,'(''  stickwidth setlinewidth'')')
       write(iun4,'(''  1 setlinecap'')')
       write(iun4,'(''  2 setlinejoin'')')
       write(iun4,'(''  newpath'')')
       write(iun4,'(''  x1 y1 moveto'')')
       write(iun4,'(''  x2 y2 lineto'')')
       write(iun4,'(''  stroke'')')
       write(iun4,'(''  grestore } def'')')
       write(iun4,'(''/dobg'')')
       write(iun4,'(''{ gsave'')')
       write(iun4,'(''  setgray'')')
       write(iun4,'(''  newpath'')')
       write(iun4,'(''  0 0 m 2000 0 l 2000 2000 l 0 2000 l'')')
       write(iun4,'(''  closepath fill'')')
       write(iun4,'(''  grestore } def'')')
       write(iun4,'(''/Helvetica-Bold findfont '',
     + ''[ fontwidth   0   0  fontheight   0   0] makefont setfont'')')
       write(iun4,'(''originx originy translate'')')
       write(iun4,'(''angle rotate'')')
       write(iun4,'(''   3 setlinewidth'')')
       write(iun4,'(''2 setlinecap'')')
       write(iun4,'(''size size scale'')')
       write(iun4,'(''%%EndProlog'')')
       write(iun4,'(''%%Page: ? 1'')')
       write(iun4,'(''dobackground {0.5 dobg} if'')')
       write(iun4,'(''n'')')
      endif

      return
      end
