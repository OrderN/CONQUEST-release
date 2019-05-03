      subroutine plmod(coo,qat,rzp,ixp,iyp,
     &                 iconn,ianz,iaton,iatclr,iresid,inat,lab,labhet,
     &                 xv,yv,
     &                 icalf,ncalf,icxp,icyp,
     &                 scal,scali)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      parameter (mxcon=10)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /athlp/ iatoms, mxnat
      common /cllab/ iclon,iclpnt(4)
      character*2 clstr,ccell
      common /clchr/ clstr(4)
      common /surf/  natorg,noscnd
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /gracom/ uscl,colscd,colscpd,ivdwpl

      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      logical dozme
      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme
      parameter (mxres=42)
      character*3 aminos
      common /amino/aminos(mxres)
      logical dsurf,scnd,dcell
      integer dolabs,fancy,shade,atcol,persp,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo

      real xx
      integer*2 ixx(24)
      dimension zvect(3)
      dimension coo(3,*),rzp(*),ixp(*),iyp(*),iconn(mxcon+1,*),
     &          ianz(*),iaton(*),iatclr(*),iresid(*),inat(*),qat(*)
      dimension icalf(6,*),icxp(*),icyp(*),lab(*),labhet(*)

      toang  = 0.52917706d0
      roddef = 0.7d0*vrad(1)
      zvect(1) = 0.0d0
      zvect(2) = 0.0d0
      zvect(3) = 1.0d0
c
c     colsc = color darker when object further away scaling
c
      colsc = colscd*uscl
c
c     icltan = color darker when a line under angle 
c              with light direction scaling factor
c
      icltan = 8

      call gethei(ihigh)

      if (iatoms.eq.0) return

      do j=1,iatoms
          call rott(coo(1,j),coo(2,j),coo(3,j),xc,yc,zc,1)
          rzp(j) = zc
          if (iaton(j).ge.1) then
             ixp(j) = (0.5d0 -(xc - xv)/scal)*ihigh
             iyp(j) = (0.5d0 -(yc - yv)/scal)*ihigh
          else
             ixp(j) = 0
             iyp(j) = 0
          endif
      end do

c      call shsort(iatoms,rzp,inat)
      call rqsrt(iatoms,rzp,inat)

      if (ipdbon.eq.1) then
         do k=1,ncalf
             j = icalf(1,k)
             call rott(coo(1,j),coo(2,j),coo(3,j),xc,yc,zc,1)
             icxp(k) = (0.5d0 -(xc - xv)/scal)*ihigh
             icyp(k) = (0.5d0 -(yc - yv)/scal)*ihigh
         end do
      endif

      inr = 0
      call labnr(inr)

      if (fancy.eq.1) then

c     color + spheres for atoms

         do i=1,iatoms
            k = inat(i)
            dsurf = .false.
            if (k.gt.natorg.and.natorg.ne.0) dsurf = .true.
            scnd = .false.
            if (backb.eq.1.and.iresid(k).le.0.and.iresid(k).ge.-3) 
     &          scnd = .true.
            if (iaton(k).ge.1) then
            xs = (rzp(k)/scali-1.0d0)/2.0d0
c            rfac = (1.0d0-colsc*xs*xs)
            rfac = (1.0d0-1.d0*xs*xs)
            if (rfac.lt.0.0d0) rfac = 0.d0
            itemp = 9*rfac
            if (.not.dsurf) then
               ia = 0
               if (iaton(k).ge.2) ia = iaton(k)
               if (ivdwpl.eq.1.and..not.scnd) then
                  isph = (vdwr(ianz(k))/(toang*scal))*ihigh*2
               else
                  isph = (vrad(ianz(k))/scal)*ihigh*2
               endif
            
               if (shade.eq.1) then
                  if (atcol.eq.1) then
                     i1 = 16+(icol(ianz(k))-1)*10
                  else
                     i1 = 16+(iatclr(k)-1)*10
                  endif
                  i2 = itemp
                  call plsph3(iyp(k),ixp(k),isph,ia,i1,i2)
               else
                  if (atcol.eq.1) then
                     kcol = icol(ianz(k))
                  else
                     kcol = iatclr(k)
                  endif
                  call setcol(kcol)
                  call plsph(iyp(k),ixp(k),isph,ia)
               endif
               
               call stlab(iresid(k),ilab,idres,iqon,lab,labhet)
               if (dolabs.eq.1.or.idres.gt.0) then
                   if (shade.eq.1) then
                      kcol = 56.0 + itemp
                   else
                      kcol = 15.0
                   endif
                   call setcol(kcol)
                   call pllab(ixp(k),iyp(k),ianz(k),k,
     &                        qat(k),inr,ilab,iresid(k),0)
               endif
            endif
            if (dsurf.and.iclon.eq.1) then
               dcell = .false.
               do ik=1,4
                  if (k.eq.iclpnt(ik)) then
                     dcell = .true.
                     ccell = clstr(ik)
                  endif
               end do
                
               if (dcell) then
                  if (shade.eq.1) then
                     kcol = 56 + itemp
                  else
                     kcol = 15
                  endif
                  call setcol(kcol)
                  call drwstr(ixp(k),iyp(k),ccell,2,k)
               endif
            endif
            xx = 3.0
            call cwidth(xx)
            call sollin

            if (iconn(1,k).gt.mxcon)
     &          write(iun3,*) 'plmol bigger than mxcon'
            do j=1,iconn(1,k)
              im = iconn(j+1,k)
              m = abs(im)
              if (iaton(m).ge.1.and.rzp(k).le.rzp(m)) then
                  if (dsurf) then
                     if (atcol.eq.1) then
                       if (iatclr(m).ne.iatclr(k)) then
                         call astick(zvect,colsc,icltan,m,k,ihigh,
     &                             shade,0,1)
                       else
                         call mstick(zvect,colsc,icltan,m,k,shade,0)
                       endif
                     else
                         call astick(zvect,colsc,icltan,
     &                               m,k,ihigh,shade,0,0)
                     endif
                  else
                     if (im.lt.0) then
                         call astick(zvect,colsc,icltan,
     &                               m,k,ihigh,shade,1,0)
                     else
                       if (ivdwpl.eq.0.or.scnd) then
                         call fstick(roddef,rfac,zvect,colsc,icltan,
     &                           m,k,ihigh,shade,atcol,scnd)
                       endif
                     endif
                  endif
              endif
            end do

            if (fyesno.eq.1) 
     &          call plfc(shade,ixx,k,ihigh,colsc,icltan,zvect)
            if (idipon.eq.1) call pldip(shade,ihigh,colsc,icltan,zvect,
     &                                  dipo,coo,ianz,xv,yv,scal,scali)

            endif
         end do
      else

        if (atcol.eq.1) then
c
c        multi colored sticks
c
         do i=1,iatoms
            k = inat(i)
            if (iaton(k).ge.1) then
            scnd = .false.
            if (backb.eq.1.and.iresid(k).le.0.and.iresid(k).ge.-3) 
     &          scnd = .true.
            dsurf = .false.
            if (k.gt.natorg.and.natorg.ne.0) dsurf = .true.
            ia = 0
            if (iatoms.eq.1) ia = 2
            if (iaton(k).ge.2) ia = iaton(k)
            if (ia.ne.0) call plsel(iyp(k),ixp(k),ia)
            call stlab(iresid(k),ilab,idres,iqon,lab,labhet)
            if (((dolabs.eq.1.or.idres.ge.1).and..not.dsurf).or.
     &          (dsurf.and.iclon.eq.1)) then
                if (shade.eq.1) then
                   xs = (rzp(k)/scali-1.0d0)/2.0d0
                   rfac = (1.0d0-1.d0*xs*xs)
                   if (rfac.lt.0.0d0) rfac = 0.d0
                   itemp = 9*rfac
                   kcol = 56 + itemp
                else
                   kcol = 15
                endif
                call setcol(kcol)
                if (dsurf.and.iclon.eq.1) then
                   dcell = .false.
                   do ik=1,4
                      if (k.eq.iclpnt(ik)) then
                         dcell = .true.
                         ccell = clstr(ik)
                      endif
                   end do
                   if (dcell)
     &                call drwstr(ixp(k),iyp(k),ccell,2,k)
                else
                   call pllab(ixp(k),iyp(k),ianz(k),k,
     &                        qat(k),inr,ilab,iresid(k),0)
                endif
            endif
            nc = 0
            do j=1,iconn(1,k)
              im = iconn(j+1,k)
              m = abs(im)
              if (iaton(m).ge.1.and.rzp(k).le.rzp(m)) then
                 if (scnd.and..not.dsurf) then
                    call mstick(zvect,colsc,icltan,m,k,shade,0)
                 else
                    id = 0
                    if (im.lt.0) id = 1
                    call astick(zvect,colsc,icltan,m,k,ihigh,shade,id,0)
                 endif
              endif
              if (iaton(m).ge.1.and.im.gt.0) nc = nc + 1
            end do
            if (nc.eq.0.and..not.dsurf) 
     &           call snglat(zvect,scal,colsc,icltan,k,ihigh,
     &                               shade,atcol,0,0,0)
            if (fyesno.eq.1) 
     &          call plfc(shade,ixx,k,ihigh,colsc,icltan,zvect)
            if (idipon.eq.1) call pldip(shade,ihigh,colsc,icltan,zvect,
     &                                  dipo,coo,ianz,xv,yv,scal,scali)
            endif
         end do
        else
c
c        mono colored sticks
c
         do i=1,iatoms
            k = inat(i)
            if (iaton(k).ge.1) then
            dsurf = .false.
            if (k.gt.natorg.and.natorg.ne.0) dsurf = .true.
            ia = 0
            if (iatoms.eq.1) ia = 2
            if (iaton(k).ge.2) ia = iaton(k)
            if (ia.ne.0) call plsel(iyp(k),ixp(k),ia)
            call stlab(iresid(k),ilab,idres,iqon,lab,labhet)
            if (((dolabs.eq.1.or.idres.ge.1).and..not.dsurf).or.
     &          (dsurf.and.iclon.eq.1)) then
                if (shade.eq.1) then
                   xs = (rzp(k)/scali-1.0d0)/2.0d0
                   rfac = (1.0d0-1.d0*xs*xs)
                   if (rfac.lt.0.0d0) rfac = 0.d0
                   itemp = 9*rfac
                   kcol = 56 + itemp
                else
                   kcol = 15
                endif
                call setcol(kcol)
                if (dsurf.and.iclon.eq.1) then
                   dcell = .false.
                   do ik=1,4
                      if (k.eq.iclpnt(ik)) then
                         dcell = .true.
                         ccell = clstr(ik)
                      endif
                   end do
                   if (dcell)
     &                call drwstr(ixp(k),iyp(k),ccell,2,k)
                else
                   call pllab(ixp(k),iyp(k),ianz(k),k,
     &                        qat(k),inr,ilab,iresid(k),0)
                endif
            endif
            nc = 0
            do j=1,iconn(1,k)
              im = iconn(j+1,k)
              m = abs(im)
              if (m.ne.0.and.iaton(m).ge.1.and.rzp(k).le.rzp(m)) then
                 if (im.lt.0) then
                    call astick(zvect,colsc,icltan,m,k,ihigh,shade,1,1)
                 else
                    if (iatclr(m).ne.iatclr(k)) then
                       call astick(zvect,colsc,icltan,m,k,ihigh,
     &                             shade,0,1)
                    else
                       call mstick(zvect,colsc,icltan,m,k,shade,0)
                    endif
                 endif
              endif
              if (iaton(m).ge.1.and.im.gt.0) nc = nc + 1
            end do
            if (nc.eq.0.and..not.dsurf) 
     &          call snglat(zvect,scal,colsc,icltan,k,ihigh,
     &                               shade,atcol,0,0,0)
            if (fyesno.eq.1) 
     &          call plfc(shade,ixx,k,ihigh,colsc,icltan,zvect)
            if (idipon.eq.1) call pldip(shade,ihigh,colsc,icltan,zvect,
     &                                  dipo,coo,ianz,xv,yv,scal,scali)
            endif
         end do
        endif

      endif

      xx = 1.0
      call cwidth(xx)
      call sollin

      call plalab(iqon)
      call reslab
      
      return
      end

      subroutine plfd(shade,ixx,k,ihigh,colsc,icltan,zvect,
     &                coo,rzp,ixp,iyp,
     &                xv,yv,scal,scali)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat)
      real xx
      integer*2 ixx(24)
      integer shade
      character str*13
      dimension tmp1(3),tmp2(3),zvect(3)
      dimension coo(3,*),rzp(*),ixp(*),iyp(*)

      do l=1,3
          tmp1(l) = fc(l,k) - coo(l,k)
      end do

      call rott(tmp1(1),tmp1(2),tmp1(3),tmp2(1),tmp2(2),tmp2(3),0)
      call impsc(tmp2,zvect,cot1)

      itemp = abs(icltan*cot1)
      xs    = (rzp(k)/scali-1.d0)/2.0d0

      rfac  = (1.0d0-colsc*xs*xs)
      if (rfac.lt.0.0d0) rfac = 0.d0
      if (rfac.gt.1.0d0) rfac = 1.d0

      itmp = 5*rfac

      if (shade.eq.1) then
         kcol = 16+(9-itemp)*rfac
      else
         kcol = 25
      endif

      call setcol(kcol)

      xx = (itmp+5)/2+0.5
      if (xx.lt.0.5) xx = 1.0

      call cwidth(xx)
      call sollin

      ixx(2) = ixp(k)
      ixx(1) = iyp(k)

      call rotts(fc(1,k),fc(2,k),fc(3,k),xc,yc,zc,1)

      ixx(4) = (0.5d0 - (xc-xv)/scal)*ihigh
      ixx(3) = (0.5d0 - (yc-yv)/scal)*ihigh

      call drawseg(ixx,1,0)

      return
      end

      subroutine pldip(shade,ihigh,colsc,icltan,zvect,dipo,
     &                 coo,ianz,xv,yv,scal,scali)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      common /athlp/ iatoms, mxnat
      real xx
      integer*2 ixx(4)
      integer shade
      character str*13
      dimension tmp1(3),tmp2(3),tmp3(3),zvect(3)
      dimension coo(3,*),ianz(*),dipo(*)

      call cntvec(tmp1,coo,ianz,iatoms)

      do l=1,3
          tmp2(l) = dipo(l) - tmp1(l)
      end do

      call rott(tmp2(1),tmp2(2),tmp2(3),tmp3(1),tmp3(2),tmp3(3),0)
      call impsc(tmp3,zvect,cot1)

      call rott(tmp1(1),tmp1(2),tmp1(3),xc,yc,zc,1)
      ixx(2) = (0.5d0 -(xc - xv)/scal)*ihigh
      ixx(1) = (0.5d0 -(yc - yv)/scal)*ihigh

      itemp = abs(icltan*cot1)
      xs    = (zc/scali-1.d0)/2.0d0

      rfac  = (1.0d0-colsc*xs*xs)
      if (rfac.lt.0.0d0) rfac = 0.d0
      if (rfac.gt.1.0d0) rfac = 1.d0

      itmp = 5*rfac

      if (shade.eq.1) then
         kcol = 16+(9-itemp)*rfac
      else
         kcol = 25
      endif

      call setcol(kcol)

      xx = (itmp+5)/2+0.5
      if (xx.lt.0.5) xx = 1.0

      call cwidth(xx)
      call sollin


      call rott(dipo(1),dipo(2),dipo(3),xc,yc,zc,1)

      ixx(4) = (0.5d0 - (xc-xv)/scal)*ihigh
      ixx(3) = (0.5d0 - (yc-yv)/scal)*ihigh

      call drawseg(ixx,1,0)

      return
      end

      subroutine msticd(zvect,colsc,icltan,m,k,shade,idash,
     &                  coo,rzp,ixp,iyp,iatclr,scali)
      implicit double precision (a-h,o-z)
      integer shade
      real xx
      integer*2 ixx(24)
      character str*13
      dimension tmp1(3),tmp2(3),zvect(3)
      dimension coo(3,*),rzp(*),ixp(*),iyp(*),iatclr(*)

      do l=1,3
          tmp1(l) = coo(l,m) - coo(l,k)
      end do

      call rott(tmp1(1),tmp1(2),tmp1(3),tmp2(1),tmp2(2),tmp2(3),0)
      call impsc(tmp2,zvect,cot1)

      itemp = icltan*cot1
      xs = (((rzp(k) + 0.5d0*(rzp(m) - rzp(k)))/scali)-1.d0)/2.0d0

      rfac = (1.0d0-colsc*xs*xs)
      if (rfac.lt.0.0d0) rfac = 0.d0

      itmp = 5*rfac
      if (shade.eq.1) then
         kcol = 16+(iatclr(k)-1)*10+(9-itemp)*rfac
      else
         kcol = iatclr(k)
      endif

      call setcol(kcol)

      xx = (itmp+5)/2+0.5
      if (xx.lt.0.5) xx = 1.0

      call cwidth(xx)
      call dash(idash)

      ixx(1) = iyp(k)
      ixx(2) = ixp(k)
      ixx(3) = iyp(m)
      ixx(4) = ixp(m)

      call drawseg(ixx,1,0)
      call dash(0)

      return
      end

      subroutine ssticd(coo,ixp,iyp,iaton,iconn,xy,yv,scal)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxseg=1000)
      common /athlp/ iatoms, mxnat
      real xx
      integer*2 ixx(mxseg*4)
      dimension coo(3,*),ixp(*),iyp(*),iaton(*),iconn(mxcon+1,*)

      xx = 1.0
      call cwidth(xx)
      call setcol(15)

      call gethei(ihigh)

      if (iatoms.eq.0) return

      do j=1,iatoms
          call rott(coo(1,j),coo(2,j),coo(3,j),xc,yc,zc,1)
          if (iaton(j).ge.1) then
             ixp(j) = (0.5d0 -(xc - xv)/scal)*ihigh
             iyp(j) = (0.5d0 -(yc - yv)/scal)*ihigh
          else
             ixp(j) = 0
             iyp(j) = 0
          endif
      end do

      nseg = 0

      do i=1,iatoms

         if (iaton(i).ge.1) then

            do j=1,iconn(1,i)
               m = iconn(1+j,i)
               nn = nseg*4

               ixx(nn+1) = iyp(i)
               ixx(nn+2) = ixp(i)
               ixx(nn+3) = iyp(m)
               ixx(nn+4) = ixp(m)

               nseg = nseg + 1

               if (nseg.eq.mxseg) then
                  call drwseg(ixx,mxseg,0)
                  nseg = 0
               endif

            end do

         endif

      end do

      if (nseg.gt.0) call drwseg(ixx,mxseg,0)

      return
      end

      subroutine asticd(zvect,colsc,icltan,m,k,ihigh,shade,idash,
     &                  imon,coo,rzp,ixp,iyp,ianz,iatclr,xv,yv,
     &                  scal,scali)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      integer shade
      real xx
      integer*2 ixx(24)
      dimension temp(6),tmp1(3),tmp2(3),zvect(3)
      dimension coo(3,*),rzp(*),ixp(*),iyp(*),ianz(*),iatclr(*)

      ik = ianz(k)
      im = ianz(m)

      do l=1,3
         temp(l) = (coo(l,m) - coo(l,k))/2.0d0 + coo(l,k)
      end do

      call rott(temp(1),temp(2),temp(3),xc,yc,zc,1)

      ixt = (0.5d0 -(xc - xv)/scal)*ihigh
      iyt = (0.5d0 -(yc - yv)/scal)*ihigh

      do l=1,3
         tmp1(l) = coo(l,m) - coo(l,k)
      end do

      call rott(tmp1(1),tmp1(2),tmp1(3),tmp2(1),tmp2(2),tmp2(3),0)
      call impsc(tmp2,zvect,cot1)

      itemp = icltan*cot1
      xs = (rzp(k)/scali-1.0d0)/2.0d0

      rfac = (1.0d0-colsc*xs*xs)
      if (rfac.lt.0.0d0) rfac = 0.d0
      itmp = 5*rfac

      if (imon.eq.1) then
         icolt = iatclr(k)
      else
         if (ik.eq.6.and.icol(ik).eq.14) then
            icolt = 10
         else
            icolt = icol(ik)
         endif
      endif

      if (shade.eq.1) then
         kcol = 16+(icolt-1)*10+(9-itemp)*rfac
      else
         kcol = icolt
      endif

      call setcol(kcol)

      xx = (itmp+5)/2+0.5
      if (xx.lt.0.5) xx = 1.0

      call cwidth(xx)
      call dash(idash)

      ixx(1) = iyp(k)
      ixx(2) = ixp(k)
      ixx(3) = iyt
      ixx(4) = ixt

      call drawseg(ixx,1,0)

      xs = (rzp(m)/scali-1.0d0)/2.0d0
      rfac = (1.0d0-colsc*xs*xs)
      if (rfac.lt.0.0d0) rfac = 0.d0
      itmp = 5*rfac

      if (imon.eq.1) then
         icolt = iatclr(m)
      else
         if (im.eq.6.and.icol(im).eq.14) then
            icolt = 10
         else
            icolt = icol(im)
         endif
      endif

      if (shade.eq.1) then
         kcol = 16+(icolt-1)*10+(9-itemp)*rfac
      else
         kcol = icolt
      endif

      call setcol(kcol)

      xx = (itmp+5)/2+0.5
      if (xx.lt.0.5) xx = 1.0

      call cwidth(xx)
      call dash(idash)

      ixx(1) = iyp(m)
      ixx(2) = ixp(m)
      ixx(3) = iyt
      ixx(4) = ixt

      call drawseg(ixx,1,0)
      call dash(0)

      if (idash.eq.1) call pldst(k,m,ixt,iyt,0)

      return
      end

      subroutine snglad(zvect,scalp,colsc,icltan,k,ihigh,shade,
     &                  atcol,ipersp,ipost,icolps,
     &                  coo,rzp,ixp,iyp,ianz,iatclr,
     &                  xv,yv,zv,c0,scal,scali)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /rdwr/  iun1,iun2,iun3,iun4,iun5

      integer shade,atcol
      real xx
      integer*2 ixx(4)
      dimension tmp1(3),tmp2(3),zvect(3)
      dimension cv(3,6),rzt(6),ixt(6),iyt(6),inat(6)
      dimension coo(3,*),rzp(*),ixp(*),iyp(*),ianz(*),iatclr(*)

      ik = ianz(k)
      rsiz = 0.4d0

      do i=1,6
         do j=1,3
            cv(j,i) = 0.0d0
         end do
      end do
      do i=1,3
         cv(i,i) = rsiz
         cv(i,i+3) = -rsiz
      end do
      

      do i=1,6

         do j=1,3
            tmp1(j) = coo(j,k) + cv(j,i)
         end do

         call rott(tmp1(1),tmp1(2),tmp1(3),xc,yc,zc,1)

         rzt(i) = zc
         if (ipersp.eq.1) then
            rp = c0/(zv - zc + c0)
            ixt(i) = (0.5d0 -rp*(xc - xv)/scalp)*ihigh
            iyt(i) = (0.5d0 -rp*(yc - yv)/scalp)*ihigh
         else
            if (ipost.eq.1) then
               ixt(i) = (0.5d0 +(xc - xv)/scalp)*ihigh
               iyt(i) = (0.5d0 -(yc - yv)/scalp)*ihigh
            else
               ixt(i) = (0.5d0 -(xc - xv)/scal)*ihigh
               iyt(i) = (0.5d0 -(yc - yv)/scal)*ihigh
            endif
         endif
      end do
      call shsort(6,rzt,inat)

      if (ipost.eq.1) then
         call rott(coo(1,k),coo(2,k),coo(3,k),xc,yc,zc,1)
         ixpp = (0.5d0 +(xc - xv)/scalp)*ihigh
         iypp = (0.5d0 -(yc - yv)/scalp)*ihigh
      endif

      do i=1,6
         m = inat(i)

         do l=1,3
            tmp1(l) = cv(l,m)
         end do

         call rott(tmp1(1),tmp1(2),tmp1(3),tmp2(1),tmp2(2),tmp2(3),0)
         call impsc(tmp2,zvect,cot1)

         cot1 = dabs(cot1)
         itemp = icltan*cot1

         if (ipersp.eq.1) then
            xs = (zv - rzp(k))/scali
         else
            xs = (rzt(m)/scali-1.0d0)/2.0d0
         endif

         rfac = (1.0d0-colsc*xs*xs)
         if (rfac.lt.0.0d0) rfac = 0.d0

         if (atcol.eq.1) then
            if (ik.eq.6.and.icol(ik).eq.14) then
               icolt = 10
            else
               if (ik.le.0) then
                  icolt = 1
               else
                  icolt = icol(ik)
               endif
            endif
         else
            icolt = iatclr(k)
         endif

c postscript
         if (ipost.eq.1) then

            call plbnd(iun4,iypp,ixpp,iyt(m),ixt(m),1,1,icolps,
     &       1,icolt,cot1,0,atcol,shade,.false.,.false.,.false.)

         else
c xwindows
            if (shade.eq.1) then
               kk = (9-itemp)*rfac
               kcol = 16+(icolt-1)*10+kk
            else
               kcol = icolt
            endif

            call setcol(kcol)

            itmp = 5*rfac
            xx = (itmp+5)/2+0.5
            if (xx.lt.0.5) xx = 1.0

            call cwidth(xx)
            call dash(0)

            ixx(1) = iyp(k)
            ixx(2) = ixp(k)
            ixx(3) = iyt(m)
            ixx(4) = ixt(m)

            call drawseg(ixx,1,0)
   
         endif
      
      end do

      return
      end

      subroutine fsticd(roddef,rfac,zvect,colsc,icltan,m,k,ihigh,
     &                  shade,atcol,scnd,
     &                  coo,rzp,ianz,iatclr,
     &                  xv,yv,scal,scali)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /helpar/ helrod, ihtype
      logical scnd,twostk
      integer shade,atcol
      real cosa, sina
      dimension tmp1(3),tmp2(3),tmp3(3),zvect(3)
      dimension rtmp1(3),rtmp2(3),rtmp3(3)
      dimension coo(3,*),rzp(*),ianz(*),iatclr(*)

      tol = 1.0d-5

      ik = ianz(k)
      im = ianz(m)
      v1 = vrad(ik)
      v2 = vrad(im)
      twostk = .false.
      if (atcol.eq.1) then
         if (ik.ne.im) twostk = .true.
      else
         if (iatclr(k).ne.iatclr(m)) twostk = .true.
      endif

      if (ik.eq.100.and.im.eq.100) then
          rodrad = helrod
      elseif (ik.ne.1.and.im.ne.1) then
c          rodrad = roddef*2
          rodrad = 0.46d0*v1
          if (v2.lt.v1) rodrad = 0.46d0*v2
      else
          rodrad = roddef
      endif
      rodrd2 = rodrad*rodrad

      do l=1,3
          tmp1(l) = coo(l,m) - coo(l,k)
      end do

      call rott(tmp1(1),tmp1(2),tmp1(3),tmp2(1),tmp2(2),tmp2(3),0)
      call impsc(tmp2,zvect,cot1)

      do l=1,3
          tmp1(l) = coo(l,k) - coo(l,m)
      end do

      if (v1.lt.rodrad.and.v2.lt.rodrad) then

          do l=1,3
             tmp2(l) = coo(l,k)
          end do

          do l=1,3
             tmp3(l) = coo(l,m)
          end do 

      else

          sc1 = dsqrt(dabs(v1*v1-rodrd2))

          call vsc1(tmp1,sc1,tol)

          do l=1,3
             tmp2(l) = coo(l,k) - tmp1(l)
          end do

          sc2 = dsqrt(dabs(v2*v2-rodrd2))

          call vsc1(tmp1,sc2,tol)

          do l=1,3
             tmp3(l) = coo(l,m) + tmp1(l)
          end do 

      endif

      do l=1,3
          if (twostk) then
             tmp1(l) = (tmp3(l) - tmp2(l))/2.0d0 + tmp2(l)
          else
             tmp1(l) = tmp3(l)
          endif
      end do

      call rott(tmp1(1),tmp1(2),tmp1(3),rtmp1(1),rtmp1(2),rtmp1(3),1)
      call rott(tmp2(1),tmp2(2),tmp2(3),rtmp2(1),rtmp2(2),rtmp2(3),1)
      call rott(tmp3(1),tmp3(2),tmp3(3),rtmp3(1),rtmp3(2),rtmp3(3),1)

      itemp = 5*cot1
      idiam1 = (rodrad/scal)*ihigh

      if (shade.eq.1) then
          if (atcol.eq.1.and..not.scnd) then
             i1 = 16 +(icol(ik)-1)*10
          else
             i1 = 16 +(iatclr(k)-1)*10
          endif
          i2 = (9-itemp)*rfac+.5d0
      else
          if (atcol.eq.1.and..not.scnd) then
             i1 = icol(ik)
          else
             i1 = iatclr(k)
          endif
      endif

      do l=1,3
          tmp1(l) = rtmp2(l) - rtmp1(l)
      end do

      clll = dsqrt(tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2))
      if (clll.lt.tol) clll = tol
      irodl = (clll/scal)*ihigh + 1

      call impsc(tmp1,zvect,cosb)

      idiam2 = (dabs(cosb)*rodrad/scal)*ihigh
      cosa = -tmp1(1)/clll
      sina = tmp1(2)/clll

      icy = (0.5d0 -(rtmp2(1) - xv)/scal)*ihigh
      icx = (0.5d0 -(rtmp2(2) - yv)/scal)*ihigh

      if (shade.eq.1) then
          call plrod3(icx,icy,idiam1,idiam2,idiam1,irodl,
     &                i1,i2,cosa,sina)
      else
          call plrodx(icx,icy,idiam1,idiam2,idiam1,irodl,
     &                i1,cosa,sina)
      endif

      if (twostk) then

          xs = (rzp(m)/scali-1.0d0)/2.0d0
          rfac = (1.0d0-1.d0*xs*xs)
          if (rfac.lt.0.0d0) rfac = 0.d0

          if (shade.eq.1) then
             if (atcol.eq.1.and..not.scnd) then
                i1 = 16 +(icol(im)-1)*10
             else
                i1 = 16 +(iatclr(m)-1)*10
             endif
          else
             if (atcol.eq.1.and..not.scnd) then
                i1 = icol(im)
             else
                i1 = iatclr(m)
             endif
          endif

          i2 = (9-itemp)*rfac+.5d0

          do l=1,3
             tmp1(l) = rtmp1(l) - rtmp3(l)
          end do

          clll = dsqrt(tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2))

          if (clll.lt.tol) clll = tol
          irodl = (clll/scal)*ihigh + 1

          call impsc(tmp1,zvect,cosb)

          idiam2 = (dabs(cosb)*rodrad/scal)*ihigh

          cosa = -tmp1(1)/clll
          sina = tmp1(2)/clll

          icy = (0.5d0 -(rtmp1(1) - xv)/scal)*ihigh
          icx = (0.5d0 -(rtmp1(2) - yv)/scal)*ihigh

          if (shade.eq.1) then
              call plrod3(icx,icy,idiam1,idiam2,idiam1,irodl,
     &                    i1,i2,cosa,sina)
          else
              call plrodx(icx,icy,idiam1,idiam2,idiam1,irodl,
     &                    i1,cosa,sina)
          endif
      endif

      return
      end

      subroutine fndmzz(ixyz,istat,imap)

c this is really fndmap

      implicit double precision (a-h,o-z)
      common /zmfrst/ ihaszm, nz, mxzat
      dimension imap(*)

      istat = 0
      do i=1,nz
         if (imap(i).eq.ixyz) then
             istat = i
             return
         endif
      end do

      return
      end

      subroutine pldst(k,m,ixt,iyt,ipltyp)
      implicit double precision (a-h,o-z)
      parameter (maxdm=20)
      common /distmn/ rdm(maxdm),idmon(2,maxdm),ndm
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character vstr*8

      if (ipltyp.eq.0) call setcol(15)
      do i=1,ndm
         if ((idmon(1,i).eq.k.and.idmon(2,i).eq.m).or.
     &       (idmon(2,i).eq.k.and.idmon(1,i).eq.m)) then
             write(vstr,'(f8.3)') rdm(i)
             if (ipltyp.eq.0) then
                call drwstr(ixt,iyt,vstr,8,-1)
             else
                write(iun4,*) ixt,' ',iyt,'  m (',vstr,') show'
             endif
         endif
      end do

      return
      end

      subroutine pldstd(ixp,iyp)
      implicit double precision (a-h,o-z)
      parameter (maxdm=20)
      common /distmn/ rdm(maxdm),idmon(2,maxdm),ndm
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character vstr*8
      dimension ixp(*),iyp(*)

      call setcol(15)

      do i=1,ndm
         ixt = (ixp(idmon(1,i)) + ixp(idmon(2,i))) / 2
         iyt = (iyp(idmon(1,i)) + iyp(idmon(2,i))) / 2
         write(vstr,'(f8.3)') rdm(i)
         call drwstr(ixt,iyt,vstr,8,-1)
      end do

      return
      end

      subroutine pllad(ixp,iyp,ianz,k,qat,inr,iqon,ires,ipost,
     &                 ityp,ipdbt,iamino,reson)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (mxppmf=16)
      parameter (mxlpmf=26)
      parameter (mxmol2=41)
      parameter (mxmm3=164)
      parameter (mxchm=82)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxchtp=136)
      parameter (mxmsf=235)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxamo=201)
      parameter (mxgmx=53)
      parameter (mxgmx2=57)
      parameter (mxg43=49)
      parameter (mxres=42)
      character str*19
      character*5 hstr 
      character*8 zst
      integer fndmap
      character*2 elemnt
      common /elem/elemnt(mxel)
      integer*2 ityp,ipdbt
      common /types/ iff
      character*2 ppmf, lpmf
      character*5 mol2
      character*19 mm3
      character*20 chmtnk
      character*20 ambstr
      character*20 amostr
      character*4 chmsf
      common /ftypes/ihasl(11),mol2(mxmol2),mm3(mxmm3),chmtnk(mxchtp),
     &               chmsf(mxmsf),ambstr(mxamb),amostr(mxamo),
     &               ppmf(mxppmf),lpmf(mxlpmf)
      character*3 pdbsym,hsym,chtnk,ambtnk
      character*2 amotnk,gffstr
      common /symbol/ pdbsym(mxsym),hsym(mxhsym),chtnk(mxchtp),
     &                ambtnk(mxamb),amotnk(mxamo),gffstr(mxgff)
      character*5 gro43,grogmx,grog2x,gro43l,grogmxl,grog2xl
      character*35 gro43s, grogms, grog2s
      common /symgro/ gro43(mxg43),grogmx(mxgmx),grog2x(mxgmx2),
     &                gro43l(mxg43),grogmxl(mxgmx),grog2xl(mxgmx2),
     &                gro43s(mxg43),grogms(mxgmx),grog2s(mxgmx2)
      integer reson
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /nmr/    shlnuc(numatm),ihsnmr
      dimension ityp(*),ipdbt(*),iamino(*),reson(*)

      if (ianz.eq.100) return

      nstr = 0
      str = "                  "

      if (iqon.eq.0) then

         str = elemnt(ianz)
         nstr = 2

      elseif (iqon.eq.1) then

         if (inr.eq.1) then
            call zzrstr(fndmap(k),zst,6,nstr)
            str =  elemnt(ianz)//zst(1:nstr)
            nstr = nstr + 2
         else
            call zzrstr(k,zst,6,nstr)
            str =  elemnt(ianz)//zst(1:nstr)
            nstr = nstr + 2
         endif

      elseif (iqon.eq.2) then

         il = ityp(k)
         ift = iff + 1
         if (ihasl(ift).gt.0.and.(il.gt.0.or.iff.eq.7)) then
            nstr = ihasl(ift)
            if (iff.eq.1) then
               str = mm3(il)
            elseif (iff.eq.2) then
               str = chtnk(il)
            elseif (iff.eq.3) then
               str = ambtnk(il)
            elseif (iff.eq.4) then
               str = amotnk(il)
            elseif (iff.eq.5) then
               str = mol2(il)
            elseif (iff.eq.6) then
               str = chmsf(il)
            elseif (iff.eq.7) then
               if (il.gt.0) then
                  str = ambtnk(il)
               elseif (il.lt.0) then
                  str = gffstr(iabs(il))
               endif
            elseif (iff.eq.8) then
               if (ires.gt.0) then
                  str = ppmf(il)
               else
                  str = lpmf(il)
               endif
            elseif (iff.eq.9) then
               str = grogmxl(il)
            elseif (iff.eq.10) then
               str = grog2xl(il)
            elseif (iff.eq.11) then
               str = gro43l(il)
            endif
         elseif (ihasl(ift).eq.0.and.il.gt.0) then
            str = hstr(il)
            nstr = 5
         else
            str = ' '
            nstr = 1
         endif

      elseif (iqon.eq.3) then

         if (ipost.eq.1) then
            write(str,'(a2,a1,f7.4)') elemnt(ianz),' ',qat
         endif
         nstr = 10

      elseif (iqon.eq.4) then

         if (ires.gt.0) then
            iam = iamino(ires)
            if (iam.gt.0.and.iam.le.mxres) then
               if (reson(ires).eq.1) then
                  if (ianz.eq.1) then
c                     ih = (ipdbt(k)-1)/3
c                     str = hsym(ih+1)
c                     nstr = 3
            ik = ipdbt(k)
            ih = (ipdbt(k)-1)/3
            call zzrstr(ik,zst,6,nstr)
            str =  hsym(ih+1)//zst(1:nstr)
            nstr = nstr + 3
                  else
c                     str = pdbsym(ipdbt(k))
c                     nstr = 3
            ik = ipdbt(k)
            call zzrstr(ik,zst,6,nstr)
            str =  pdbsym(ipdbt(k))//zst(1:nstr)
            nstr = nstr + 3
                  endif
               endif
            endif
         endif

      elseif (iqon.eq.6) then

         if (ipost.eq.1) then
            write(str,'(a2,a1,f7.4)') elemnt(ianz),' ',qat
         endif
         nstr = 10

      endif

      if (nstr.gt.0) then
         if (ipost.eq.1) then
            write(iun4,'(''labelcol'')')
            write(iun4,*) iyp,' ',ixp,' m (',str(1:nstr),') show'
         else
            if (iqon.eq.3) then
c this exception is to avoid recursive IO in fortran
               call drwqstr(ixp,iyp,ianz,qat,k)
            elseif (iqon.eq.6) then
               shl = shlnuc(k)
               call drwqstr(ixp,iyp,ianz,shl,k)
            else
               call drwstr(ixp,iyp,str,nstr,k)
            endif
         endif
      endif

      return
      end

      subroutine plalad(iqon,rzp,ixp,iyp,
     &                  icalf,ianf,islu,nchain,iamino,reson,
     &                  irsnr,achain,scali)
      implicit double precision (a-h,o-z)
      logical dozme
      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme
      integer reson
      integer dolabs,fancy,shade,atcol,persp,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      parameter (mxel=100)
      parameter (mxres=42)
      character*3 aminos
      common /amino/aminos(mxres)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      character str*13
      character*1 achain
      dimension rzp(*),ixp(*),iyp(*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),reson(*),
     &          irsnr(*),achain(*)

      if (iqon.eq.5.and.ipdbon.eq.1.and.dolabs.eq.1) then

         nstr = 8
         if (nchain.gt.1) nstr = 10

         do k=1,nchain

            do j=ianf(k),islu(k)

               jj = j - ianf(k) + 1

               if (reson(j).eq.1) then

                  l = icalf(1,j)

                  if (shade.eq.1) then
                     xs = (rzp(l)/scali-1.0d0)/2.0d0
                     rfac = (1.0d0-1.d0*xs*xs)
                     if (rfac.lt.0.0d0) rfac = 0.d0
                     kcol = 156 + 9*rfac
                  else
                     kcol = 15
                  endif

                  call setcol(kcol)

                  str(1:3) = aminos(iamino(j))
                  str(4:4) = ' '
                  write(str(5:8),'(i4)') irsnr(j)

                  if (nchain.gt.1) then
                     str(9:9) = '.'
                     str(10:10) = achain(j)
                  endif

                  call drwstr(ixp(l),iyp(l),str,nstr,l)

               endif

            end do

         end do

      end if

      return
      end

      subroutine reslad(rzp,ixp,iyp,
     &                  icalf,nchain,iamino,reson,
     &                  irsnr,lab,ihets,
     &                  achain,ncalf,scali)
      implicit double precision (a-h,o-z)
      integer reson
      parameter (mxres=42)
      parameter (mxheta=150)
      character*3 aminos
      common /amino/aminos(mxres)
      character*3 hetz
      common /clfstr/ ihashz,ihetq(mxheta),ihqset(mxheta),ihhadd(mxheta)
     &                ,labhet(mxheta),ilcset,ligcat(mxheta),hetz(mxheta)
      character str*13
      character*1 achain
      dimension rzp(*),ixp(*),iyp(*)
      dimension icalf(6,*),iamino(*),reson(*),
     &          irsnr(*),lab(*),ihets(*),achain(*)


      nstr = 8
      if (nchain.gt.1) nstr = 10

      do j=1,ncalf

         if (reson(j).eq.1.and.lab(j).eq.7) then

            l = icalf(1,j)

            if (shade.eq.1) then
               xs = (rzp(l)/scali-1.0d0)/2.0d0
               rfac = (1.0d0-1.d0*xs*xs)
               if (rfac.lt.0.0d0) rfac = 0.d0
               kcol = 156 + 9*rfac
            else
               kcol = 15
            endif

            call setcol(kcol)

            str(1:3) = aminos(iamino(j))
            str(4:4) = ' '
            write(str(5:8),'(i4)') irsnr(j)

            if (nchain.gt.1) then
               str(9:9) = '.'
               str(10:10) = achain(j)
            endif

            call drwstr(ixp(l),iyp(l),str,nstr,-1)

         endif

      end do

      nstr = 4
      do i=1,mxheta
         if (ligcat(i).ne.-1.and.ihets(i).ge.1.and.
     &       labhet(i).eq.6) then
            l = ligcat(i) + 1

            if (shade.eq.1) then
               xs = (rzp(l)/scali-1.0d0)/2.0d0
               rfac = (1.0d0-1.d0*xs*xs)
               if (rfac.lt.0.0d0) rfac = 0.d0
               kcol = 156 + 9*rfac
            else
               kcol = 15
            endif

            call setcol(kcol)

            str(1:3) = hetz(i)
            str(4:4) = ' '

            call drwstr(ixp(l),iyp(l),str,nstr,-1)
         endif
      end do

      return
      end


      subroutine plald(iami,ixp,iyp,
     &                 icalf,iamino,irsnr)
      implicit double precision (a-h,o-z)
      parameter (mxres=42)
      character*3 aminos
      common /amino/aminos(mxres)
      character str*13
      dimension ixp(*),iyp(*),icalf(6,*),iamino(*),irsnr(*)

      nstr = 8
      l = icalf(1,iami)
      kcol = 15

      call setcol(kcol)

      str(1:3) = aminos(iamino(iami))
      str(4:4) = ' '
      write(str(5:8),'(i4)') irsnr(iami)

      call drwstr(ixp(l),iyp(l),str,nstr,l)

      return
      end

      subroutine stlab(irs,ilab,idres,iqon,lab,labhet)
      implicit double precision (a-h,o-z)
      dimension lab(*),labhet(*)

      idres = 0
      ilab = iqon
      if (irs.gt.0) then
         if (lab(irs).gt.0) then
            idres = 1
            ilab = lab(irs)-2
         endif
      else
         if (labhet(abs(irs)+1).gt.0) then
            idres = 1
            ilab = labhet(abs(irs)+1)-2
         endif
      endif

      return
      end

