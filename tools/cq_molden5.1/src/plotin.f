      subroutine plotin
      implicit double precision (a-h,o-z)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /plot /  iplot,iplwin,icolps
      common /keywrd/ keywrd,keyori
      character*320 keywrd,keyori

c-------- determine graphics device ----------------------
c
c     Keyword        Description                        iplot
c
c     FIGURE         figure (gcg software)              0
c     HPGL           Hewlett Packard Graphics Language  1
c     HP2393         Hewlett Packard terminal           2
c     TEK4014        Tektronics 4014                    3
c     POSTSCRIPT     Postscript for laserwriters        4
c     SILLY,OPENGL   OpenGL Version                     5
c     XWINDOWS       Xwindow                            6
c
      if (index(keywrd,'FIG').ne.0) iplot = 0
      if (index(keywrd,'HPGL').ne.0) iplot = 1
      if (index(keywrd,'HP23').ne.0) iplot = 2
      if (index(keywrd,'TEK4').ne.0) iplot = 3
      if (index(keywrd,'POST').ne.0) iplot = 4
      if (index(keywrd,'SILLY').ne.0) iplot = 5
      if (index(keywrd,'OPENGL').ne.0) iplot = 5
      if (index(keywrd,'XWIN').ne.0) iplot = 6

      return
      end
