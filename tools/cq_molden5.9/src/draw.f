      subroutine draw(iz,jz,jdir,a,ib,ic,id,pz,mdim)
      implicit double precision (a-h,o-z)
      dimension a(mdim,mdim), idirx(6), idiry(6)
      logical ai,aj,ak,lc(4), yes, euclid
***********************************************************************
*
*   DRAW DRAWS CONTOURS, STARTING AT POINT IZ,JZ.
*   JDIR = STARTING DIRECTION, IF 1 THEN IN +X DIRECTION
*                                 2 THEN IN +Y DIRECTION
*                                 3 THEN IN -X DIRECTION
*                                 4 THEN IN -Y DIRECTION
*    A   = ARRAY TO BE PLOTTED
*    B,C = WORK SPACES OF SIZE SUITABLE TO HOLD 5000 AND 4000 LOGICAL
*          ELEMENTS.
*    PZ  = VALUE OF THIRD DIMENSION. SET TO ZERO IF NOT WANTED.
*    MDIM= SIZE OF FIRST DIMENSION OF A.
***********************************************************************
      common /cntval/ cntval, euclid, yes, icnt, iflag
      common/cntr/r1,imax,jmax
      equivalence (lc(1),ilc)
      dimension id(mdim,*),ib(*),ic(2000,2)
      data idirx/1,0,-1,0,1,0/,idiry/0,1,0,-1,0,1/

      icount =-5
      idir = jdir
      ilc = 0
      ipen = 1
      i = iz
      j = jz
      iflab = 1
      if (iflag.eq.0) iflab = 0
    1 aa = a(i,j)
      ai = (aa.lt.0.)
      idx = idirx(idir)
      idy = idiry(idir)
      ab = a(i+idx,j+idy)
      factor= aa/(aa-ab)
      py = (imax-i-factor*idx)/(imax-1)
      px = (jmax-j-factor*idy)/(jmax-1)

      if (iflab.eq.0) then
         iflab = 1
         yes = .true.
         call euler(px,py,pz,ipen)
         yes = .false.
      endif

      if (ipen.eq.2.or.r1.gt.0.d0) 
     &      r1 = r1 + dsqrt((cx-px)**2 + (cy-py)**2)

      if (id(i,j).eq.1.and.id(i+idx,j+idy).eq.1)then
        call propnt(px,py,pz,ipen)
        ipen = 2
      else
        ipen = 1
      endif

      cy = py
      cx = px

      icount = icount + 1
      if (idir/2*2.ne.idir) goto 5
      if (idir.ne.2) goto 6
      number = j
      m = 1
      goto 7
    6 number= j-1
      m = imax
    7 k = (number-1)*imax +i
      if (ib(k).eq.0) return
      ib(k) = 0
      if (i.eq.m) return
      goto 3
    5 if (j.ne.1.and.j.ne.jmax) goto 3
      number = idir*j
      if (number.ne.3.and.number.ne.jmax) goto 3
      if (number.eq.3) then
          ic(i-1,1) = 0
      else
          ic(i,2) = 0
      endif
      return
    3 continue
      iddx = i + idirx(idir+1)
      iddy = j + idiry(idir+1)
      ac = a(iddx+idx,iddy+idy)
      aj = (ac.lt.0.)
      ad = a(iddx,iddy)
      ak = (ad.lt.0.)
      fa = 1
      if (aj.and.ak) goto 10
      if (.not.(aj.or.ak)) goto 11
      fa=0.
      if ((ai.or..not.ak).and.(.not.ai.or.ak)) goto 4
      fa= aa*ac-ab*ad
      if (fa.ge.0.) goto 4
   15 idir = idir + 1
      if (idir.eq.5) idir = 1
      goto 1
 4    i = iddx
      j = iddy
      if (fa.eq.0) goto 1
      i = i + idx
      j = j + idy
      idir = idir - 1
      if (idir.eq.0) idir = 4
      goto 1
   10 if (ai) goto 4
      goto 15
 11   if(ai) goto 15
      goto 4

      end
