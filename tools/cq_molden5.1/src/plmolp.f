      subroutine plmold(coo,qat,rzp,ixp,iyp,
     &                 iconn,ianz,iaton,iatclr,iresid,inat,
     &                 xv,yv,zv,c0,
     &                 icalf,ncalf,ianf,islu,nchain,iamino,reson,
     &                 icxp,icyp,irsnr,achain,scali)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
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
      integer reson
      character*1 achain
      character*3 aminos
      common /amino/aminos(mxres)
      logical dsurf,scnd,dcell,doq
      integer dolabs,fancy,shade,atcol,persp,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo

      real xx
      integer*2 ixx(24)
      character str*13
      dimension zvect(3)
      dimension coo(3,*),rzp(*),ixp(*),iyp(*),iconn(mxcon+1,*),
     &          ianz(*),iaton(*),iatclr(*),iresid(*),inat(*),qat(*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),reson(*),
     &          icxp(*),icyp(*),irsnr(*)

      toang  = 0.52917706d0
      doq = (ihasq.eq.1.and.iqon.eq.3)
      roddef = 0.7d0*vrad(1) 
      zvect(1) = 0.0d0
      zvect(2) = 0.0d0
      zvect(3) = 1.0d0

      scalp = 20.0d0
c      scalp = scali / 2.0d0
      icltan = 8
      colscp = colscpd*uscl

      call gethei(ihigh)

c     coordinates viewpoint

      do j=1,iatoms
          call rott(coo(1,j),coo(2,j),coo(3,j),xc,yc,zc,1)
          rzp(j) = zc
          rp = c0/(zv - rzp(j) + c0)
          if (iaton(j).ge.1.and.(zv-rzp(j)).gt.0.0d0) then
             ixp(j) = (0.5d0 -rp*(xc - xv)/scalp)*ihigh
             iyp(j) = (0.5d0 -rp*(yc - yv)/scalp)*ihigh
          else
             ixp(j) = 0
             iyp(j) = 0
          endif
      end do

      if (ipdbon.eq.1) then
         do k=1,ncalf
             j = icalf(1,k)
             call rott(coo(1,j),coo(2,j),coo(3,j),xc,yc,zc,1)
             rp = c0/(zv - zc + c0)
             if ((zv-zc).gt.0.0d0) then
                icxp(k) = (0.5d0 -rp*(xc - xv)/scalp)*ihigh
                icyp(k) = (0.5d0 -rp*(yc - yv)/scalp)*ihigh
             else
                icxp(k) = 0
                icyp(k) = 0
             endif
         end do
      endif

      inr = 0
      call shsort(iatoms,rzp,inat)

      inr = 0
      call labnr(inr)

      if (fancy.eq.1) then

         do i=1,iatoms
c Max number of shades per color is 10:(0-9)
            k = inat(i)
            dsurf = .false.
            if (k.gt.natorg.and.natorg.ne.0) dsurf = .true.
            scnd = .false.
            if (backb.eq.1.and.iresid(k).le.0.and.iresid(k).ge.-3) 
     &          scnd = .true.
            if (iaton(k).ge.1.and.(zv-rzp(k)).gt.0.0d0) then
            tmp = c0 / (zv - rzp(k) + c0)
            xs = (zv-rzp(k))/scali
            rfac = (1.0d0-colscp*xs*xs)
            if (rfac.lt.0.0d0) rfac = 0.d0
            if (.not.dsurf) then
               itemp = 9*rfac
               ia = 0
               if (iaton(k).ge.2) ia = iaton(k)
               if (ivdwpl.eq.1.and..not.scnd) then
                  isph = (vdwr(ianz(k))*tmp/(toang*scalp))*ihigh*2.4d0
     &                   + 0.5d0
               else
                  isph = (vrad(ianz(k))*tmp/scalp)*ihigh*2 + 0.5d0
               endif
               if (shade.eq.1) then
                  if (atcol.eq.1) then
                     i1 = 16 +(icol(ianz(k))-1)*10
                  else
                     i1 = 16 +(iatclr(k)-1)*10
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
               
               if (dolabs.eq.1) then
                   if (shade.eq.1) then
                      kcol = 56 + 9*rfac
                   else
                      kcol = 15
                   endif
                   call setcol(kcol)
                   call pllab(ixp(k),iyp(k),ianz(k),k,
     &                        qat(k),inr,iqon,iresid(k),0)

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
                     kcol = 56 + 9*rfac
                  else
                     kcol = 15
                  endif
                  call setcol(kcol)
                  call drwstr(ixp(k),iyp(k),ccell,2)
               endif
            endif


            xx = 3.0
            call cwidth(xx)
            call sollin

            if(iconn(1,k).gt.mxcon)
     &          write(iun3,*) 'plmol bigger than mxcon'
            do j=1,iconn(1,k)
              im = iconn(j+1,k)
              m = abs(im)
              if (iaton(m).ge.1.and.rzp(k).le.rzp(m)
     &            .and.(zv-rzp(m)).gt.0.0d0) then
                  if (dsurf) then
                     if (atcol.eq.1) then
                         if (iatclr(k).ne.iatclr(m)) then
                            call astck
     &                  (zvect,colscp,scalp,icltan,m,k,ihigh,shade,0,1)
                         else
                            call mstck(zvect,colscp,icltan,m,k,shade,0)
                         endif
                     else
                         call astck(zvect,colscp,scalp,icltan,
     &                              m,k,ihigh,shade,0,0)
                     endif
                  else
                     if (im.lt.0) then
                       call astck(zvect,colscp,scalp,icltan,
     &                              m,k,ihigh,shade,1,0)
                     else
                      if (ivdwpl.eq.0.or.scnd) 
     &                 call fstck(roddef,rfac,scalp,zvect,colscp,icltan,
     &                          m,k,ihigh,shade,atcol,scnd)
          
                     endif
                  endif
              endif
            end do
            if (fyesno.eq.1) 
     &          call plfcp(shade,ixx,k,ihigh,colscp,icltan,zvect,scalp)
            if (idipon.eq.1) call pldipp(shade,ihigh,colscp,icltan,zvect
     &                           ,dipo,scalp,coo,ianz,xv,yv,zv,c0,scali)
            endif
         end do
      else

        if (atcol.eq.1) then
         do i=1,iatoms
            k = inat(i)
            if (iaton(k).ge.1.and.(zv-rzp(k)).gt.0.0d0) then
            scnd = .false.
            if (backb.eq.1.and.iresid(k).le.0.and.iresid(k).ge.-3) 
     &          scnd = .true.
            dsurf = .false.
            if (k.gt.natorg.and.natorg.ne.0) dsurf = .true.
            ia = 0
            if (iatoms.eq.1) ia = 2
            if (iaton(k).ge.2) ia = iaton(k)
            if (ia.ne.0) call plsel(iyp(k),ixp(k),ia)
            if ((dolabs.eq.1.and..not.dsurf).or.
     &          (dsurf.and.iclon.eq.1)) then
                if (shade.eq.1) then
                   xs = (zv-rzp(k))/scali
                   rfac = (1.0d0-colscp*xs*xs)
                   if (rfac.lt.0.0d0) rfac = 0.d0
                   kcol = 56 + 9*rfac
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
     &                call drwstr(ixp(k),iyp(k),ccell,2)
                else
                   call pllab(ixp(k),iyp(k),ianz(k),k,
     &                        qat(k),inr,iqon,iresid(k),0)
                endif
            endif
            nc = 0
            do j=1,iconn(1,k)
              im = iconn(j+1,k)
              m = abs(im)
              if (iaton(m).ge.1.and.rzp(k).le.rzp(m)
     &            .and.(zv-rzp(m)).gt.0.0d0) then
                if (scnd.and..not.dsurf) then
                  call mstck(zvect,colscp,icltan,m,k,shade,0)
                else
                  id = 0
                  if (im.lt.0) id = 1
                  call astck
     &                 (zvect,colscp,scalp,icltan,m,k,ihigh,shade,id,0)
                endif
              endif
              if (iaton(m).ge.1.and.im.gt.0) nc = nc + 1
            end do
            if (nc.eq.0) call snglat(zvect,scalp,colscp,icltan,k,ihigh,
     &                               shade,atcol,1,0,0)
            if (fyesno.eq.1) 
     &          call plfcp(shade,ixx,k,ihigh,colscp,icltan,zvect,scalp)
            if (idipon.eq.1) call pldipp(shade,ihigh,colscp,icltan,zvect
     &                           ,dipo,scalp,coo,ianz,xv,yv,zv,c0,scali)
            endif
         end do
        else
         do i=1,iatoms
            k = inat(i)
            if (iaton(k).ge.1.and.(zv-rzp(k)).gt.0.0d0) then
            ia = 0
            if (iatoms.eq.1) ia = 2
            if (iaton(k).ge.2) ia = iaton(k)
            if (ia.ne.0) call plsel(iyp(k),ixp(k),ia)
            dsurf = .false.
            if (k.gt.natorg.and.natorg.ne.0) dsurf = .true.
            if ((dolabs.eq.1.and..not.dsurf).or.
     &          (dsurf.and.iclon.eq.1)) then
                if (shade.eq.1) then
                   xs = (zv-rzp(k))/scali
                   rfac = (1.0d0-colscp*xs*xs)
                   if (rfac.lt.0.0d0) rfac = 0.d0
                   kcol = 56 + 9*rfac
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
     &                call drwstr(ixp(k),iyp(k),ccell,2)
                else
                   call pllab(ixp(k),iyp(k),ianz(k),k,
     &                        qat(k),inr,iqon,iresid(k),0)
                endif
            endif
            nc = 0
            do j=1,iconn(1,k)
              im = iconn(j+1,k)
              m = abs(im)
              if (iaton(m).ge.1.and.rzp(k).le.rzp(m)
     &            .and.(zv-rzp(m)).gt.0.0d0) then
                 if (im.lt.0) then
                     call astck
     &                 (zvect,colscp,scalp,icltan,m,k,ihigh,shade,1,1)
                 else
                    if (iatclr(k).ne.iatclr(m)) then
                       call astck
     &                 (zvect,colscp,scalp,icltan,m,k,ihigh,shade,0,1)
                    else
                       call mstck(zvect,colscp,icltan,m,k,shade,0)
                    endif
                 endif
              endif
              if (iaton(m).ge.1.and.im.gt.0) nc = nc + 1
            end do
            if (nc.eq.0) call snglat(zvect,scalp,colscp,icltan,k,ihigh,
     &                               shade,atcol,1,0,0)
            if (fyesno.eq.1) 
     &          call plfcp(shade,ixx,k,ihigh,colscp,icltan,zvect,scalp)
            if (idipon.eq.1) call pldipp(shade,ihigh,colscp,icltan,zvect
     &                           ,dipo,scalp,coo,ianz,xv,yv,zv,c0,scali)
            endif
         end do
        endif

      endif

      xx = 1.0
      call cwidth(xx)
      call sollin

      if (iqon.eq.5.and.ipdbon.eq.1.and.dolabs.eq.1) then
         nstr = 8
         if (nchain.gt.1) nstr = 10
         do k=1,nchain
            do j=ianf(k),islu(k)
               jj = j - ianf(k) + 1
               if (reson(j).eq.1) then
                  l = icalf(1,j)
                  if (zv-rzp(l).gt.0.0d0) then
                     if (shade.eq.1) then
                        xs = (zv-rzp(l))/scali
                        rfac = (1.0d0-colscp*xs*xs)
                        if (rfac.lt.0.0d0) rfac = 0.d0
                        kcol = 156 + 9*rfac
                     else
                        kcol = 15
                     endif
                     call setcol(kcol)
                     str(1:3) = aminos(iamino(j))
                     str(4:4) = ' '
c                     write(str(5:8),'(i4)') jj
                     write(str(5:8),'(i4)') irsnr(j)
                     if (nchain.gt.1) then
                        str(9:9) = '.'
c                        if (k.lt.10) str(10:10) = char(k+48)
                        if (k.lt.10) str(10:10) = achain(j)
                     endif
                     call drwstr(ixp(l),iyp(l),str,nstr)
                  endif
               endif
            end do
         end do
      end if
      
      return
      end

      subroutine plfcd(shade,ixx,k,ihigh,colsc,icltan,zvect,scalp,
     &                 coo,rzp,ixp,iyp,
     &                 xv,yv,zv,c0,scali)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat)
      dimension tmp1(3),tmp2(3),zvect(3)
      real xx
      integer*2 ixx(24)
      integer shade
      dimension coo(3,*),rzp(*),ixp(*),iyp(*)

      do l=1,3
          tmp1(l) = fc(l,k) - coo(l,k)
      end do

      call rott(tmp1(1),tmp1(2),tmp1(3),tmp2(1),tmp2(2),tmp2(3),0)
      call impsc(tmp2,zvect,cot1)

      itemp = abs(icltan*cot1)
      xs = (zv-rzp(k))/scali
      rfac = (1.0d0-colsc*xs*xs)

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

      rp = c0/(zv - zc + c0)
      ixx(4) = (0.5d0 - rp*(xc-xv)/scalp)*ihigh
      ixx(3) = (0.5d0 - rp*(yc-yv)/scalp)*ihigh

      call drawseg(ixx,1,0)

      return
      end

      subroutine pldipp(shade,ihigh,colsc,icltan,zvect,dipo,scalp,
     &                 coo,ianz,xv,yv,zv,c0,scali)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      common /athlp/ iatoms, mxnat
      real xx
      integer*2 ixx(4)
      integer shade
      dimension tmp1(3),tmp2(3),tmp3(3),zvect(3)
      dimension coo(3,*),ianz(*),dipo(*)

      call cntvec(tmp1,coo,ianz,iatoms)

      do l=1,3
          tmp2(l) = dipo(l) - tmp1(l)
      end do

      call rott(tmp2(1),tmp2(2),tmp2(3),tmp3(1),tmp3(2),tmp3(3),0)
      call impsc(tmp3,zvect,cot1)

      call rott(tmp1(1),tmp1(2),tmp1(3),xc,yc,zc,1)
      rp = c0/(zv - zc + c0)
      ixx(2) = (0.5d0 - rp*(xc-xv)/scalp)*ihigh
      ixx(1) = (0.5d0 - rp*(yc-yv)/scalp)*ihigh

      itemp = abs(icltan*cot1)
      xs = (zv-zc)/scali

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

      rp = c0/(zv - zc + c0)
      ixx(4) = (0.5d0 - rp*(xc-xv)/scalp)*ihigh
      ixx(3) = (0.5d0 - rp*(yc-yv)/scalp)*ihigh

      call drawseg(ixx,1,0)

      return
      end

      subroutine mstcd(zvect,colsc,icltan,m,k,shade,idash,
     &                 coo,rzp,ixp,iyp,iatclr,zv,scali)
      implicit double precision (a-h,o-z)
      integer shade
      real xx
      integer*2 ixx(24)
      dimension tmp1(3),tmp2(3),zvect(3)
      dimension coo(3,*),rzp(*),ixp(*),iyp(*),iatclr(*)

      do l=1,3
         tmp1(l) = coo(l,m) - coo(l,k)
      end do

      call rott(tmp1(1),tmp1(2),tmp1(3),tmp2(1),tmp2(2),tmp2(3),0)
      call impsc(tmp2,zvect,cot1)

      itemp = icltan*cot1
      xs = (zv -(rzp(k) + 0.5d0*(rzp(m) - rzp(k))))/scali
      rfac = (1.0d0-colsc*xs*xs)

      if (rfac.lt.0.0d0) rfac = 0.d0
      if (shade.eq.1) then
         kcol = 16 +(iatclr(k)-1)*10+(9-itemp)*rfac
      else
         kcol = iatclr(k)
      endif

      call setcol(kcol)

      itmp = 6*rfac
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

      subroutine astcd(zvect,colsc,scalp,icltan,m,k,ihigh,shade,
     &                 idash,imon,coo,rzp,ixp,iyp,ianz,iatclr,
     &                 xv,yv,zv,c0,scali)
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
         tmp1(l) = coo(l,m) - coo(l,k)
      end do

      call rott(tmp1(1),tmp1(2),tmp1(3),tmp2(1),tmp2(2),tmp2(3),0)
      call impsc(tmp2,zvect,cot1)

      itemp = icltan*cot1

      do l=1,3
         temp(l) = (coo(l,m) - coo(l,k))/2.0d0 + coo(l,k)
      end do

      call rott(temp(1),temp(2),temp(3),xc,yc,zc,1)

      rzt = zc
      rp = c0/(zv - rzt + c0)

      ixt = (0.5d0 -rp*(xc - xv)/scalp)*ihigh
      iyt = (0.5d0 -rp*(yc - yv)/scalp)*ihigh

      xs = (zv -(rzp(k) + 0.5d0*(rzt - rzp(k))))/scali
      rfac = (1.0d0-colsc*xs*xs)

      if (rfac.lt.0.0d0) rfac = 0.d0
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

      itmp = 6*rfac
      xx = (itmp+5)/2+0.5
      if (xx.lt.0.5) xx = 1.0

      call cwidth(xx)
      call dash(idash)

      ixx(1) = iyp(k)
      ixx(2) = ixp(k)
      ixx(3) = iyt
      ixx(4) = ixt

      call drawseg(ixx,1,0)

      xs = (zv -(rzp(m) + 0.5d0*(rzt - rzp(m))))/scali
      rfac = (1.0d0-colsc*xs*xs)
      if (rfac.lt.0.0d0) rfac = 0.d0

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

      itmp = 6*rfac
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

      subroutine fstcd(roddef,rfac,scalp,zvect,colsc,icltan,m,k,
     &                 ihigh,shade,atcol,scnd,coo,rzp,ianz,iatclr,
     &                 xv,yv,zv,c0,scali)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /helpar/ helrod, ihtype
      logical scnd,twostk
      integer shade,atcol
      real xx, cosa, sina
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
            tmp1(l) = coo(l,m)
         endif
      end do

      call rott(tmp1(1),tmp1(2),tmp1(3),rtmp1(1),rtmp1(2),rtmp1(3),1)
      call rott(tmp2(1),tmp2(2),tmp2(3),rtmp2(1),rtmp2(2),rtmp2(3),1)
      call rott(tmp3(1),tmp3(2),tmp3(3),rtmp3(1),rtmp3(2),rtmp3(3),1)

      rp1 = c0/(zv - rtmp2(3) + c0)
      rp2 = c0/(zv - rtmp1(3) + c0)

      idiam1 = (rp1*rodrad/scalp)*ihigh
      idiam3 = (rp2*rodrad/scalp)*ihigh

      itemp = 5*cot1

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

      icy1 = (0.5d0 -rp1*(rtmp2(1) - xv)/scalp)*ihigh
      icx1 = (0.5d0 -rp1*(rtmp2(2) - yv)/scalp)*ihigh
      icy2 = (0.5d0 -rp2*(rtmp1(1) - xv)/scalp)*ihigh
      icx2 = (0.5d0 -rp2*(rtmp1(2) - yv)/scalp)*ihigh

      do l=1,3
         tmp1(l) = rtmp2(l) - rtmp1(l)
      end do

      call impsc(tmp1,zvect,cosb)

      idiam2 = (rp1*dabs(cosb)*rodrad/scalp)*ihigh
      tmp1(1) =  icy1 - icy2
      tmp1(2) =  icx1 - icx2

      clll = dsqrt(tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2))

      if (clll.lt.tol) clll = tol
      irodl = clll +1

      cosa = tmp1(1)/clll
      sina = -tmp1(2)/clll

      if (shade.eq.1) then
          call plrod3(icx1,icy1,idiam1,idiam2,idiam3,irodl,
     &                i1,i2,cosa,sina)
      else
          call plrodx(icx1,icy1,idiam1,idiam2,idiam3,irodl,
     &                i1,cosa,sina)
      endif

      if (twostk) then
         idiam1 = (rp1*rodrad/scalp)*ihigh
         xs = (zv-rzp(m))/scali

         rfac = (1.0d0-colsc*xs*xs)
         if (rfac.lt.0.0d0) rfac = 0.d0

         if (shade.eq.1) then
            if (atcol.eq.1.and..not.scnd) then
               i1 = 16 +(icol(im)-1)*10
            else
               i1 = 16 +(iatclr(m)-1)*10
            endif
            i2 = (9-itemp)*rfac+.5d0
         else
            if (atcol.eq.1.and..not.scnd) then
               i1 = icol(im)
            else
               i1 = iatclr(m)
            endif
         endif

         rp3 = c0/(zv - rtmp3(3) + c0)
         idiam3 = (rp3*rodrad/scalp)*ihigh

         icy3 = (0.5d0 -rp3*(rtmp3(1) - xv)/scalp)*ihigh
         icx3 = (0.5d0 -rp3*(rtmp3(2) - yv)/scalp)*ihigh

         do l=1,3
            tmp1(l) = rtmp1(l) - rtmp3(l)
         end do

         call impsc(tmp1,zvect,cosb)

         idiam2 = (rp1*dabs(cosb)*rodrad/scalp)*ihigh

         tmp1(1) =  icy3 - icy2
         tmp1(2) =  icx3 - icx2

         clll = dsqrt(tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2))

         if (clll.lt.tol) clll = tol
         irodl = clll 

         cosa = -tmp1(1)/clll
         sina = tmp1(2)/clll

         if (shade.eq.1) then
             call plrod3(icx2,icy2,idiam1,idiam2,idiam3,irodl,
     &                   i1,i2,cosa,sina)
         else
             call plrodx(icx2,icy2,idiam1,idiam2,idiam3,irodl,
     &                   i1,cosa,sina)
         endif

      endif

      return
      end
