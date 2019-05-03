      subroutine plposd(backb,dolabs,icolps,fancy,atcol,persp,shade,
     &                  idelx,coo,ianz,iaton,iresid,iatclr,iconn,
     &                  ixpp,iypp,rzpp,inat,qat,
     &                  xv,yv,zv,
     &                  icalf,ncalf,ianf,islu,nchain,iamino,reson,
     &                  scal)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (mxcon=10)
      parameter (small=1.0d-4)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /athlp/ iatoms, mxnat
      common /surf/  natorg,noscnd
      common /cllab/ iclon,iclpnt(4)
      character*2 clstr,ccell
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /clchr/  clstr(4)
      common /gracom/ uscl,colscd,colscpd,ivdwpl
      logical dozme
      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme

      logical dosol,isbck,isbck1,dsurf,dcell,dclpnt
      logical doq
      integer dolabs,fancy,persp,shade,atcol,coltp,backb
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      parameter (mxres=42)
      integer reson
      character*3 aminos
      common /amino/aminos(mxres)
      character*2 elemnt
      common /elem/ elemnt(mxel)
c
c     Add normal coordinate arrows
c
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      real frmul
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
c
      character*13 str
      dimension ixpp(*),iypp(*),rzpp(*)
      dimension ctemp1(3),ctemp2(3),ctemp3(3),ctemp4(3),zvect(3)
      dimension coo(3,*),ianz(*),iaton(*),iresid(*),iatclr(*),
     &          iconn(mxcon+1,*),inat(*),qat(*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),reson(*)

      if (iatoms.le.0) return

      toang  = 0.52917706d0
      doq = (ihasq.eq.1.and.iqon.eq.3)
      tol = small
      ipost = 2000
      roddef = 0.7d0*vrad(1)
      scle = scal/1.2d0
      zvect(1) = 0.0d0
      zvect(2) = 0.0d0
      zvect(3) = 1.0d0


      rzmax = 0.0d0

      do j=1,iatoms
         call rott(coo(1,j),coo(2,j),coo(3,j),xc,yc,zc,1)
         ixpp(j) = (0.5d0 +(xc - xv)/scle)*ipost
         iypp(j) = (0.5d0 -(yc - yv)/scle)*ipost
         rzpp(j) = zc
         rtmp = dabs(rzpp(j))
         if (rtmp.gt.rzmax) rzmax = rtmp
      end do

      call shsort(iatoms,rzpp,inat)

      inr = 0
      call labnr(inr)

      do i=1,iatoms
         k = inat(i)
         if (persp.eq.1.and.zv-rzpp(k).lt.0.0d0) goto 100
         dsurf = .false.
         if (k.gt.natorg.and.natorg.ne.0) dsurf = .true.
         if (iaton(k).ge.1.and..not.(ianz(k).eq.99.and.idelx.eq.1)) then
         isbck1 = .false.
         if (backb.eq.1.and.iresid(k).gt.0) then
            nn = 6
            if (iamino(iresid(k)).le.23) nn=3
            do ii=1,nn
               call fndcal(k,ii,idum,istat,icalf,ncalf)
               if (istat.eq.1) isbck1 = .true.
            end do
         endif
         if ((backb.eq.0.and.fancy.eq.1).or.
     &      (backb.eq.1.and.fancy.eq.1.and.(iresid(k).lt.-3.or.
     &      (iresid(k).gt.0.and.reson(iresid(k)).eq.1)))) then
            if (ivdwpl.eq.1) then
               isph = vdwr(ianz(k))/(toang*scle)*ipost
               if (dsurf) isph = 0
            else
               isph = (vrad(ianz(k))/scle)*ipost
            endif
            write(iun4,*)'/rx ',iypp(k),' def'
            write(iun4,*)'/ry ',ixpp(k),' def'
            if (icolps.eq.1) then
               if (atcol.eq.1) then
                  write(iun4,*)'/col ',icol(ianz(k)),' def'
               else
                  write(iun4,*)'/col ',iatclr(k),' def'
               endif
            else
               write(iun4,*)'/col 0 def'
            endif
            write(iun4,*)'/rad ',isph,' def'
            write(iun4,*)'doatom'
         endif
         dcell = .false.
         dclpnt = .false.
         if (dsurf.and.iclon.eq.1) then
            do ik=1,4
              if (k.eq.iclpnt(ik)) then
                  dcell = .true.
                  ccell = clstr(ik)
              endif
            end do
            if (k.le.natorg+8) dclpnt = .true.
         endif
        
         nc = 0
         do j=1,iconn(1,k)
            im = iconn(j+1,k)
            m = abs(im)
            if (iaton(m).ge.1.and.im.gt.0) nc = nc + 1
            if (iaton(m).ge.1.and.rzpp(k).le.rzpp(m).and.
     &          .not.(ianz(m).eq.99.and.idelx.eq.1)) then
             isbck = .false.
             if (backb.eq.1.and.iresid(m).gt.0.and.isbck1) then
                nn = 6
                if (iamino(iresid(m)).le.23) nn=3
                do ii=1,nn
                   call fndcal(m,ii,idum,istat,icalf,ncalf)
                   if (istat.eq.1) isbck = .true.
                end do
             endif
             if (((backb.eq.0.and.fancy.eq.1).or.
     &            (backb.eq.1.and.(iresid(k).lt.-3.or.
     &           (iresid(k).gt.0.and..not.isbck))
     &           .and.fancy.eq.1)).and..not.dsurf) then
c
c stick between displayed atoms
c
              
              v1 = vrad(ianz(k))
              v2 = vrad(ianz(m))
              if (ianz(k).ne.1.and.ianz(m).ne.1) then
c                  rodrad = roddef*2
                  rodrad = 0.46d0*v1
                  if (v2.lt.v1) rodrad = 0.46d0*v2
              else
                  rodrad = roddef
              endif
              do l=1,3
                 ctemp1(l) = coo(l,k) - coo(l,m)
              end do
              sc1 = dsqrt(dabs(vrad(ianz(k))**2-rodrad**2))
              call vsc1(ctemp1,sc1,tol)
              do l=1,3
                 ctemp2(l) = coo(l,k) - ctemp1(l)
              end do
              sc2 = dsqrt(dabs(vrad(ianz(m))**2-rodrad**2))
              call vsc1(ctemp1,sc2,tol)
              do l=1,3
                 ctemp3(l) = coo(l,m) + ctemp1(l)
              end do
              sct = 1.0d0
             else
c              rodrad = roddef*2
              rodrad = 0.304d0
              do l=1,3
                 ctemp2(l) = coo(l,k)
                 ctemp3(l) = coo(l,m)
              end do

c this 1.2 is only there for esthetics get it out

              sct = 1.2d0*(rzpp(k) + rzpp(m))/(2.0d0*rzmax) + 2.0d0
             endif

             call rott(ctemp2(1),ctemp2(2),ctemp2(3),xc,yc,zc,1)
             iy1 = (0.5d0 +(xc - xv)/scle)*ipost
             ix1 = (0.5d0 -(yc - yv)/scle)*ipost
             call rott(ctemp3(1),ctemp3(2),ctemp3(3),xc,yc,zc,1)
             iy2 = (0.5d0 +(xc - xv)/scle)*ipost
             ix2 = (0.5d0 -(yc - yv)/scle)*ipost
             do l=1,3
                ctemp4(l) = (coo(l,m) - coo(l,k))/2.0d0 + coo(l,k)
             end do
             call rott(ctemp4(1),ctemp4(2),ctemp4(3),xc,yc,zc,1)
             iy3 = (0.5d0 +(xc - xv)/scle)*ipost
             ix3 = (0.5d0 -(yc - yv)/scle)*ipost


             ihd = (sct*rodrad/scle)*ipost

             if (backb.eq.1) then
                if (isbck.or.(iresid(k).le.0.and.iresid(k).ge.-3))
     &          then
                   dosol = .true.
                else
                   if (fancy.eq.1.and..not.dsurf) then
                      dosol = .true.
                   else
                      dosol = .false.
                   endif
                endif
             else 
                dosol = .false.
                if (fancy.eq.1.and..not.dsurf) dosol = .true.
             endif
             if (im.lt.0) dosol = .false.

             ctemp1(1) = 1.0d0*(ixpp(k)-ixpp(m))
             ctemp1(2) = 1.0d0*(iypp(k)-iypp(m))
             ctemp1(3) = ((rzpp(k)-rzpp(m))/scal)*ipost
             call impsc(ctemp1,zvect,cosb)
             cosb = dabs(cosb)
             if (cosb.lt.1.0d-07) cosb = 1.0d-07

             if (dosol) then
               call plbnd(iun4,ix1,iy1,ix2,iy2,ihd,im,icolps,iresid(k),
     &          iatclr(k),cosb,backb,atcol,shade,isbck,dosol,dclpnt)
             else
               coltp = atcol
               if (fancy.eq.1.and.dsurf) then
                   if (coltp.eq.1) then
                      coltp = 0
                   else
                      coltp = 1
                   endif
               endif

               if (coltp.eq.1.and..not.dclpnt) then
                  ic1 = icol(ianz(k))
                  ic2 = icol(ianz(m))
                  if (ianz(k).eq.6) ic1 = 10
                  if (ianz(m).eq.6) ic2 = 10
               else
                  ic1 = iatclr(k)
                  ic2 = iatclr(m)
               endif
               if (icolps.eq.1.and.ianz(k).ne.ianz(m)) then
                call plbnd(iun4,ix1,iy1,ix3,iy3,ihd,im,icolps,iresid(k),
     &          ic1,cosb,backb,atcol,shade,isbck,dosol,dclpnt)
                call plbnd(iun4,ix3,iy3,ix2,iy2,ihd,im,icolps,iresid(m),
     &          ic2,cosb,backb,atcol,shade,isbck,dosol,dclpnt)
               else
                call plbnd(iun4,ix1,iy1,ix2,iy2,ihd,im,icolps,iresid(k),
     &          ic1,cosb,backb,atcol,shade,isbck,dosol,dclpnt)
               endif
               if (im.lt.0) call pldst(k,m,ix3,iy3,1)
             endif

            endif
         end do
         if (nc.eq.0.and.fancy.eq.0) call snglat(zvect,scle,0.0d0,0,
     &                               k,ipost, shade,atcol,0,1,icolps)
         if (dcell) then
             write(iun4,'(''labelcol'')')
             write(iun4,*) iypp(k),' ',ixpp(k),
     &         ' m (',ccell,') show'
         endif


         if (dolabs.ge.1.and..not.dsurf) then
            call pllab(ixpp(k),iypp(k),ianz(k),k,qat(k),inr,iqon,
     &                 iresid(k),1)
         endif
         endif
100      continue
      end do

      if (ifrq.gt.0) then

         amplit = 15.0d0*frmul
         arrow = 0.2d0

         do j=1,iatoms

c     normal vector

            call rott(coo(1,j)+amplit*a(1,j),
     &           coo(2,j)+amplit*a(2,j),
     &           coo(3,j)+amplit*a(3,j),xc,yc,zc,1)
            ixfp = (0.5d0 +(xc - xv)/scle)*ipost
            iyfp = (0.5d0 -(yc - yv)/scle)*ipost
            dsq = dsqrt(dble((ixpp(j)-ixfp)**2+(iypp(j)-iyfp)**2))
            if (dsq.ge.(50*frthr)) then
               write(iun4,*)'/x1',iypp(j),' def'
               write(iun4,*)'/y1',ixpp(j),' def'
               write(iun4,*)'/x2',iyfp,' def'
               write(iun4,*)'/y2',ixfp,' def'
               write(iun4,*)'0 setcol'
               write(iun4,*)'dostick'

c     add arrows

               length = 50
               dx = ixpp(j)-ixfp
               if (dx.lt.0.0d0) length = -length
               if (dabs(dx).gt.1.0d-10) then
                  ang = datan((iypp(j)-iyfp)/dx)
               else
                  ang = 1.57079632679489661922
                  if ((iypp(j)-iyfp).le.0.0d0) ang=-ang
               endif
               write(iun4,*)'/x1',iyfp,' def'
               write(iun4,*)'/y1',ixfp,' def'
               write(iun4,*)'/x2',iyfp+length*sin(ang+arrow),' def'
               write(iun4,*)'/y2',ixfp+length*cos(ang+arrow),' def'
               write(iun4,*)'0 setcol'
               write(iun4,*)'dostick'
               write(iun4,*)'/x1',iyfp,' def'
               write(iun4,*)'/y1',ixfp,' def'
               write(iun4,*)'/x2',iyfp+length*sin(ang-arrow),' def'
               write(iun4,*)'/y2',ixfp+length*cos(ang-arrow),' def'
               write(iun4,*)'0 setcol'
               write(iun4,*)'dostick'

            endif

         end do

      endif

      if (iqon.eq.5.and.ipdbon.eq.1.and.dolabs.eq.1) then
         write(iun4,'(''labelcol'')')
         do k=1,nchain
            do j=ianf(k),islu(k)
               jj = j - ianf(k) + 1
               if (reson(j).eq.1) then
                  l = icalf(1,j)
                  str(1:3) = aminos(iamino(j))
                  str(4:4) = ' '
                  write(str(5:8),'(i4)') jj
                  if (nchain.gt.1) then
                     str(9:9) = '.'
                     if (k.lt.10) str(10:10) = char(k+48)
                  endif
                  write(iun4,*) iypp(l),' ',ixpp(l),' m (',str,') show'
               endif
            end do
         end do
      end if

      return
      end

      subroutine plbnd(iun4,ix1,iy1,ix2,iy2,ihd,im,icolps,ires,iatclr,
     &                 cosb,backb,atcol,shade,isbck,dosol,dclpnt)
      implicit double precision (a-h,o-z)
      parameter (cosm=0.1d0)
      logical isbck,dosol,dclpnt
      integer shade,atcol,backb

      write(iun4,*)'/x1 ',ix1,' def'
      write(iun4,*)'/y1 ',iy1,' def'
      write(iun4,*)'/x2 ',ix2,' def'
      write(iun4,*)'/y2 ',iy2,' def'
      if (.not.dosol.and.(icolps.eq.0.or.iatclr.eq.15)
     &    .and.dabs(cosb).lt.cosm) then
         rsgn = 1.0d0
         if (cosb.lt.0.0d0) rsgn = -1.0d0
         cosb = cosm*rsgn 
      endif
      write(iun4,*)'/cosb ',real(cosb),' def'

      if (dosol) then
          if (icolps.eq.1) then
             if ((atcol.eq.1.and..not.(ires.le.0.and.
     &           ires.ge.-3.and.backb.eq.1)).or.isbck) then
                write(iun4,*)'/col 0 def'
             else
                write(iun4,*)'/col ',iatclr,' def'
             endif
          else
                write(iun4,*)'/col 0 def'
          endif
          write(iun4,*)'/hd ',ihd,' def'
          write(iun4,*)'dorod'
      else
          if (icolps.eq.1) then
              if (dclpnt) then
                 write(iun4,*) 'cellcol'
              else
                 if (iatclr.eq.15.and.shade.eq.0) then
                    write(iun4,*) 0,' setcol'
                 else
                    write(iun4,*) iatclr,' setcol'
                 endif
              endif
          endif
          if (im.lt.0) then
              if (icolps.eq.1) write(iun4,*)'1 setcol'
              write(iun4,*) '0 setlinecap [8 14] 0 setdash'
              write(iun4,*) '5 setlinewidth'
          endif
          if (shade.eq.1) then
             write(iun4,*)'doshadedstick'
          else
             write(iun4,*)'dostick'
          endif
          if (im.lt.0) then
              if (icolps.eq.1) write(iun4,*)'0 setgray'
              write(iun4,*) '2 setlinecap [] 0 setdash'
              write(iun4,*) '2 setlinewidth'
          endif
      endif

      return
      end
