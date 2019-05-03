      subroutine cntoud(a,mdim,imax,jmax,pz,value,r11,id,ib)
      implicit double precision (a-h,o-z)

c THIS IS REALLY cntour

      dimension  a(mdim,*)
      logical oclose,osshad,ofill,ofrst
      integer*2 ixx
      common /fill /  ixx(1000),nixx,ihight,icol1
      common /spacom/ grey,iscol,oclose,osshad,ofill,ofrst
      common /cntr/   r1,im,jm
      dimension id(mdim,*),ib(*),ic(2000,2)

************************************************************************
*
* GENERAL CONTOUR-DRAWING ROUTINE.
*
*  ON INPUT:
*      A      = 2-D ARRAY OF POINTS DEFINING HEIGHT IN MAP-SPACE.
*      MDIM   = SIZE OF FIRST DIMENSION OF A.
*      IMAX   = LENGTH OF ONE SIDE OF ARRAY WITHIN A.
*      JMAX   = LENGTH OF THE OTHER SIDE OF THE ARRAY WITHIN A.
*      PZ     = VALUE OF Z DIRECTION IF 3-D PLOTS ARE BEING DRAWN.
*      VALUE  = VALUE OF CONTOUR.
*  ON EXIT R11= LENGTH OF CONTOUR. THIS IS THE TOTAL LENGTH OF ALL LINES
*               MAKING UP THE CONTOUR, OPEN AND CLOSED.
************************************************************************
      im = imax
      r1 = 0.d0
      jm = jmax
C
C  SET UP STRING FOR LATER ANNOTATION OF CONTOURS.
C
      number = dabs(value)
C
C  THE FOLLOWING CODE SUBTRACTS THE VALUE OF THE VALUE BEING SOUGHT
C  FROM ALL ELEMENTS OF A, SO THAT THE ZERO VALUE WILL BE SOUGHT
C  HEREAFTER. 
C
      do i=1,imax
         do j=1,jmax
            a(i,j) = a(i,j) - value
         end do
      end do

      r11 = r1
      jmax1 = jmax-1
      imax1 = imax-1

      do i=1,imax
          do j=1,jmax1
             aa = a(i,j)
             ab = a(i,j+1)
C
C  IF A CONTOUR  PASSES CLOSE TO AN ELEMENT OF A, THE CORRESPONDING
C  ELEMENT OF LOGICAL ARRAY B IS SET.
C
             if ((aa.ge.0..and.ab.lt.0.).or.
     &           (aa.lt.0..and.ab.ge.0.)) then
                ib(imax*(j-1)+i) = 1
             else
                ib(imax*(j-1)+i) = 0
             endif

          end do
      end do

      K = 1
      do j=1,jmax,jmax1
          ab = a(1,j)
          do i=1,imax1
              aa = ab
              ab = a(i+1,j)
C
C  THE ARRAY C HAS THE SAME FUNCTION AS B, BUT ON THE BOUNDARY OF THE
C  AREA.
C
              if ((aa.ge.0..and.ab.lt.0.).or.
     &            (aa.lt.0..and.ab.ge.0.)) then
                  ic(i,k) = 1
              else
                  ic(i,k) = 0
              endif
          end do
          k = 2
      end do
C
C   DRAW UNCLOSED CONTOURS (IE STARTING ON A BOUNDARY).
C
      oclose = .false.
      do i=1,imax1
          if (ic(i,1).eq.1) call draw(i,1,1,a,ib,ic,id,pz,mdim)
          if (ic(i,2).eq.1) call draw(i+1,jmax,3,a,ib,ic,id,pz,mdim)
      end do
C
C  DRAW INTERIOR CONTOURS PASSING CLOSE TO BOUNDARY POINTS.
C
      do j=1,jmax1
        if (ib((j-1)*imax+1).eq.1) 
     &     call draw(1,j+1,4,a,ib,ic,id,pz,mdim)
        if (ib(j*imax).eq.1) 
     &     call draw(imax,j,2,a,ib,ic,id,pz,mdim)
      end do

      if (ofill) then
         idum = 1
         cdum1 = 0.0d0
         cdum2 = 0.0d0
         call plotgr(idum,cdum1,cdum2)
      endif
C
C  DRAW OTHER INTERIOR CONTOURS.
C
      oclose = .true.
      if (ofill) nixx = 0

      do i=2,imax1
         do j=1,jmax1
            if (ib((j-1)*imax+i).eq.1) 
     &         call draw(i,j,2,a,ib,ic,id,pz,mdim)
         end do
      end do

      if (ofill) then
         idum = 1
         cdum1 = 0.0d0
         cdum2 = 0.0d0
         call plotgr(idum,cdum1,cdum2)
      endif
      oclose = .false.
*
*  RESTORE "A" MATRIX TO IT'S ORIGINAL VALUE.
*
      do i=1,imax
          do j=1,jmax
              a(i,j) = a(i,j) + value
          end do
      end do

      return
      end
