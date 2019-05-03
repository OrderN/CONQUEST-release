      subroutine pl3dm(adjus,ovis,height)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      character str*100
      real xx, yy
      logical ovis
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /atopla/ xsym(numatm),ysym(numatm),zsym(numatm),
     &                isym(numatm)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /coord / xyz(3,numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /plot/   iplot,iplwin,icolps
      common /zproj/  zval(numatm),nindx(numatm),nconn(numatm)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)

c     Plot Molecule in the 3d mode

      if (iplot.eq.1) write(iun4,'(a)') 'SP2;'

      if (iplot.eq.4) then
         if (.not.(icolps.eq.1.and.ovis)) 
     &      write(iun4,'(''0.5 setgray'')')
         write(iun4,'(''14 setlinewidth'')')
         write(iun4,'(''1 setlinecap'')')
      endif

      if (iplot.eq.6) then
         xx = 4.0
         call xwin(xx,yy,10,str,nstr,idum1,idum2)
         call xwin(xx,yy,8,str,nstr,idum1,idum2)
      endif

      do m=1,natoms
         i = m
         if (ovis) i = nindx(m)
         xa = 0.5d0 + xsym(i) / r(1)
         ya = 0.5d0 + ysym(i) / r(2)
         za = zsym(i) / r(3)

         do j=1,natoms
            dmaxsq = (vdwr(nat(i)) + vdwr(nat(j)))**2
            dijsq = ((xyz(1,i)-xyz(1,j))*adjus)**2
     +             +((xyz(2,i)-xyz(2,j))*adjus)**2
     +             +((xyz(3,i)-xyz(3,j))*adjus)**2
            if (dijsq.lt.dmaxsq) then
               if (iplot.eq.6) then
                  if (isym(i).eq.1.and.isym(j).eq.1) then
                     xx = 12.0
                     call unstip
                  else
                     xx = 11.0
                     call ststip
                  endif
                  yy = 0.0
                  call xwin(xx,yy,99,str,nstr,idum1,idum2)
               endif
               xb = 0.5d0 + xsym(j) / r(1)
               yb = 0.5d0 + ysym(j) / r(2)
               zb = zsym(j) / r(3)
               if (.not.ovis) then
                  call euler(xa,ya,za,1)
                  call euler(xb,yb,zb,2)
               else
                 if (zval(i).le.zval(j)) then
                  xt = (xb - xa)/2.0d0 + xa
                  yt = (yb - ya)/2.0d0 + ya
                  zt = (zb - za)/2.0d0 + za

                  isend = 0
                  if (nconn(i).eq.1) isend = 1

                  if (iplot.eq.6) then
                     xx = icol(nat(i))
                     call xwin(xx,yy,99,str,nstr,idum1,idum2)
                  elseif (iplot.eq.4.and.icolps.eq.1) then
                     icolt = icol(nat(i))
                  endif

                  if (za.le.height.and.zt.le.height) then
                     if (iplot.eq.4.and.icolps.eq.1) then
                        call psbond(xa,ya,za,xt,yt,zt,icolt,isend)
                     else
                        call euler(xa,ya,za,1)
                        call euler(xt,yt,zt,2)
                     endif
                  elseif (za.le.height) then
                     call snypnt(xa,ya,za,xt,yt,zt,height,sx,sy,sz)
                     if (iplot.eq.4.and.icolps.eq.1) then
                        call psbond(xa,ya,za,sx,sy,sz,icolt,isend)
                     else
                        call euler(xa,ya,za,1)
                        call euler(sx,sy,sz,2)
                     endif
                  elseif (zt.le.height) then
                     call snypnt(xa,ya,za,xt,yt,zt,height,sx,sy,sz)
                     if (iplot.eq.4.and.icolps.eq.1) then
                        call psbond(xt,yt,zt,sx,sy,sz,icolt,0)
                     else
                        call euler(xt,yt,zt,1)
                        call euler(sx,sy,sz,2)
                     endif
                  endif

                  isend = 0
                  if (nconn(j).eq.1) isend = 1

                  if (iplot.eq.6) then
                     xx = icol(nat(j))
                     call xwin(xx,yy,99,str,nstr,idum1,idum2)
                  elseif (iplot.eq.4.and.icolps.eq.1) then
                     icolt = icol(nat(j))
                  endif

                  if (zb.le.height.and.zt.le.height) then
                     if (iplot.eq.4.and.icolps.eq.1) then
                        call psbond(xb,yb,zb,xt,yt,zt,icolt,isend)
                     else
                        call euler(xb,yb,zb,1)
                        call euler(xt,yt,zt,2)
                     endif
                  elseif (zb.le.height) then
                     call snypnt(xb,yb,zb,xt,yt,zt,height,sx,sy,sz)
                     if (iplot.eq.4.and.icolps.eq.1) then
                        call psbond(xb,yb,zb,sx,sy,sz,icolt,isend)
                     else
                        call euler(xb,yb,zb,1)
                        call euler(sx,sy,sz,2)
                     endif
                  elseif (zt.le.height) then
                     call snypnt(xb,yb,zb,xt,yt,zt,height,sx,sy,sz)
                     if (iplot.eq.4.and.icolps.eq.1) then
                        call psbond(xt,yt,zt,sx,sy,sz,icolt,0)
                     else
                        call euler(xt,yt,zt,1)
                        call euler(sx,sy,sz,2)
                     endif
                  endif
                 endif
               endif
            endif 
         end do
      end do

      if (iplot.eq.6) then
         call unstip
         xx = 1.0
         call xwin(xx,yy,10,str,nstr,idum1,idum2)
         call xwin(xx,yy,8,str,nstr,idum1,idum2)
      endif

      if (iplot.eq.4) then
         if (.not.(icolps.eq.1.and.ovis)) write(iun4,'(''s'')')
         write(iun4,'(''0 setgray'')')
         write(iun4,'(''n'')')
      endif

      return
      end

      subroutine psbond(px,py,pz,qx,qy,qz,icol,icap)
      implicit double precision (a-h,o-z)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /plot / iplot,iplwin,icolps

      call pstran(px,py,pz,ix,iy)
      call pstran(qx,qy,qz,jx,jy)

      if (ix.eq.jx.and.iy.eq.jy) return

      write(iun4,*)'/col ',icol,' def'
      if (icap.eq.0) then
         write(iun4,*) ix,iy,jx,jy,' dobond'
      else
         write(iun4,*) ix,iy,jx,jy,ix,iy,' dobond2'
      endif

      return
      end

      subroutine pstran(px,py,pz,ix,iy)
      implicit double precision (a-h,o-z)
      logical euclid,yes
      common /cntval/ cntval, euclid, yes, icnt, iflag

      call eulerh(px,py,pz,cx,cy)

      if (euclid) then
         r2 = 1.0d0
      else
         r2 = dsqrt(2.0d0)
      endif
      cx = 0.5d0 + cx / r2
      cy = 0.5d0 + cy / r2

      cx = max(0.d0,min(1.d0,cx))
      cy = max(0.d0,min(1.d0,cy))

      ix = int(cx*2000)
      iy = int(cy*2000)+125

      return
      end
