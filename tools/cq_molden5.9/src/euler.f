      subroutine euler(px,py,pz,ipen)
      implicit double precision (a-h,o-z)
************************************************************************
*
* EULER DOES A SIMPLE EULERIAN TRANSFORM ON THE POINT PX,PY,PZ TO
* PRODUCE A POINT FOR PLOTTING. 
*    SET THE PLANE OF THE PLOT
*            PX,PY,PZ = POINT IN 3-D SPACE TO MOVE THE PEN TO.
*            IPEN     = 1 DO NOT DRAW A LINE, = 2 DRAW A LINE.
************************************************************************
      common /eul/ ca,cb,sa,sb,ccab,scab,ssab,scba
      common /plrat/ rat(3)

      cx = ccab*(px-0.5d0)*rat(1) + scab*(py-0.5d0)*rat(2) - 
     &     sb*pz*rat(3)
      cy =  -sa*(px-0.5d0)*rat(1) +   ca*(py-0.5d0)*rat(2)
      call plotgr(ipen,cx,cy)

      return
      end

      subroutine propnd(qx,qy,qz,ipen,
     &                  coo,ianz,iaton,iatclr,iresid,iconn)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (max3d=61)
      parameter (mesp = (max3d*max3d*max3d+max3d)/4)
      parameter (mdum = (max3d*max3d*max3d+max3d) - mesp*4)
      common /athlp/ iatoms, mxnat
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /spa3d/  connl(3,mesp), esp(mesp), dum(mdum)
      common /pntsta/ iproj,idrcol,nesp
      dimension coo(3,*),ianz(*),iaton(*),iatclr(*),iresid(*),
     &          iconn(mxcon+1,*)

      if (iproj.eq.1) then

          call euler(qx,qy,qz,ipen)

      elseif (iproj.eq.0) then

          if (iatoms.lt.mxnat) then
             iatoms = iatoms + 1

             x = (0.5d0-qx)*r(1)
             y = (0.5d0-qy)*r(2)

             coo(1,iatoms) = v1(1)*x + v2(1)*y + px
             coo(2,iatoms) = v1(2)*x + v2(2)*y + py
             coo(3,iatoms) = v1(3)*x + v2(3)*y + pz

             iaton(iatoms) = 1
             ianz(iatoms) = 100
             iatclr(iatoms) = idrcol
             iresid(iatoms) = -4
             if (ipen.eq.2) then
                iconn(2,iatoms) = iatoms - 1
                iconn(1,iatoms) = 1
                l = iconn(1,iatoms-1)
                if (l.lt.mxcon) then
                   iconn(l+2,iatoms-1) = iatoms
                   iconn(1,iatoms-1) = l + 1
                endif
             else
                iconn(1,iatoms) = 0
             endif
          endif

      elseif (iproj.eq.-1) then

          if (nesp.lt.mesp) then
             nesp = nesp + 1

             x = (0.5d0-qx)*r(1)
             y = (0.5d0-qy)*r(2)

             connl(1,nesp) = v1(1)*x + v2(1)*y + px
             connl(2,nesp) = v1(2)*x + v2(2)*y + py
             connl(3,nesp) = v1(3)*x + v2(3)*y + pz
          endif

      endif

      return
      end

