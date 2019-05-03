      subroutine eulerh(px,py,pz,cx,cy)
      implicit double precision (a-h,o-z)
************************************************************************
*
* EULER DOES A SIMPLE EULERIAN TRANSFORM ON THE POINT PX,PY,PZ TO
* PRODUCE A POINT FOR PLOTTING. IT HAS EXTRA CODE TO ALLOW THE USER
* TO SPECIFY WHICH TWO DIMENSIONS ARE TO BE PLOTTED.
* ON INPUT: 
*    THE FIRST CALL OF EULER (NOT EULERH) INITIALISES COMMON /EUL/
*    SET THE PLANE OF THE PLOT
*            PX,PY,PZ = POINT IN 3-D SPACE TO MOVE THE PEN TO.

************************************************************************
      common /eul/ ca,cb,sa,sb,ccab,scab,ssab,scba
      common /plrat/ rat(3)

      cy = ccab*(px-0.5d0)*rat(1) + scab*(py-0.5d0)*rat(2) - 
     &     sb*pz*rat(3)
      cx =  -sa*(px-0.5d0)*rat(1) +   ca*(py-0.5d0)*rat(2)

      return
      end
