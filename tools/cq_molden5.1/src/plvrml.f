      subroutine plvrmd(iun,fancy,atcol,dolabs,ihnd,backb,dohead,
     &                  rz,t,coo,ianz,iaton,iatclr,iresid,iconn,
     &                  icalf,ianf,islu,nchain,iamino,reson,
     &                  scal,scali)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      parameter (mxcon=10)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /athlp/ iatoms, mxnat

      common /cllab/ iclon,iclpnt(4)
      character*2 clstr,ccell
      common /clchr/ clstr(4)
      common /gracom/uscl,colscd,colscpd,ivdwpl

      logical doscnd, dohead
      integer reson
      parameter (mxres=42)
      character*3    aminos
      common /amino/ aminos(mxres)
      character*2    elemnt
      common /elem/  elemnt(mxel)
      common /surf/  natorg,noscnd
      common /vrcol/ jcol(3,16)
      common /vropt/ ivtwo,ihand,ivadd
      common /helpar/helrod, ihtype

      character str*13
      logical dsurf
      integer dolabs,fancy,atcol,backb
      dimension rz(3),t(3),tmp(3), tlpos(3), nat(mxel)
      dimension coo(3,*),ianz(*),iaton(*),iatclr(*),iresid(*),
     &          iconn(mxcon+1,*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),reson(*)

      toang  = 0.52917706d0
      roddef = 0.7d0*vrad(1)
      doscnd = .true.

      if (iatoms.eq.0) return

      ihand = ihnd

      call plcini

      if (ivtwo.eq.2) then

c povray

          xc = 0.0d0
          yc = 0.0d0
          zc = 1.0d0
          tmp(1) = zc*rz(1)
          tmp(2) = zc*rz(2)
          tmp(3) = zc*rz(3)

          do i=1,3
             tmp(i) = tmp(i)*scali*2.5d0 + t(i)
          end do

          tlpos(1) = 9.0d0
          tlpos(2) = 9.0d0
          tlpos(3) = 9.0d0

          call plphd(iun,ihand,tmp,t,tlpos)

      elseif (ivtwo.eq.1) then
         if (dohead) then
            write(iun,'(''#VRML V2.0 utf8'')')
            write(iun,*) 'NavigationInfo { type ',char(34),
     &                   'EXAMINE',char(34),' }'
            write(iun,*) 'Viewpoint { position 0 0 ',scal,
     &        ' description ',char(34),'MoldenVRML',char(34),' }'
         else
            write(iun,*) 'Group {'
            write(iun,*) ' children ['
         endif
      else
         write(iun,'(''#VRML V1.0 ascii'')')
      endif

      if (fancy.eq.1.or.ivtwo.eq.2) then
         do i=1,iatoms
           ia = ianz(i)
           if (iaton(i).ge.1.and.ia.ne.100) then
            if (atcol.eq.1) then
               ic = icol(ia)
            else
               ic = iatclr(i)
            endif
            if (ivdwpl.eq.1) then
               call plvsph(iun,jcol,ic,coo(1,i),vdwr(ia)/toang)
            else
               call plvsph(iun,jcol,ic,coo(1,i),vrad(ia))
            endif
            n = iconn(1,i)
            if (n.ne.0) then
               do j=1,n
                 im = iconn(j+1,i)
                 m = abs(im)
                 ma = ianz(m)
                 if (iaton(m).ge.1.and.ma.ne.100.and.im.gt.0) then
c                    if (ia.ne.1.and.im.ne.1) then
                        rodrad = roddef*2
c                    else
c                        rodrad = roddef
c                    endif
                    ido = 1
                    if (ma.eq.ia) then
                       if (m.gt.i) then
                          do l=1,3
                             tmp(l) = coo(l,m)
                          end do
                       else
                          ido = 0
                       endif
                    else
                       do l=1,3
                          tmp(l) = 
     &                   (coo(l,m) - coo(l,i))/2.0d0 + coo(l,i)
                       end do
                    endif
                    if (ido.eq.1) then
                           call plvrod(iun,1,jcol,ic,
     &                          coo(1,i),tmp,rodrad)
                    endif
                 endif
               end do
            endif
           endif
         end do

         if (ivtwo.eq.2) goto 100
         if (ivtwo.eq.1) then
            write(iun,*) 'Transform {'
            write(iun,*) '  children ['
         else
            write(iun,*) 'Separator {'
         endif
         do k=1,15
            nat(k) = 0
         end do
         do i=1,iatoms
            if (iaton(i).ge.1.and.(i.gt.natorg.and.natorg.ne.0)) then
               k = iatclr(i)
               nat(k) = nat(k) + 1
            endif
         end do
         do k=1,15
            if (nat(k).gt.0) then
               call plvst(iun,jcol,k,0,.true.)
            endif
         end do
         if (ivtwo.eq.1) write(iun,*) '  ]'
         write(iun,*) '}'
100      continue
      else
c
c        multi colored sticks
c
         if (ivtwo.eq.1) then
            write(iun,*) 'Transform {'
            write(iun,*) '  children ['
         else
            write(iun,*) 'Separator {'
         endif

         if (atcol.eq.1) then
            kmax = mxel
         else
            kmax = 15
         endif
         do k=1,kmax
            nat(k) = 0
         end do
         do i=1,iatoms
            if (iaton(i).ge.1) then
               if (atcol.eq.1) then
                  k = ianz(i)
               else
                  k = iatclr(i)
               endif
               nat(k) = nat(k) + 1
            endif
         end do
         do k=1,kmax
            if (nat(k).gt.0) call plvst(iun,jcol,k,atcol,.false.)
         end do
         if (ivtwo.eq.1) write(iun,*) '  ]'
         write(iun,*) '}'
      endif

      if (doscnd) then
         do i=1,iatoms
           dsurf = .false.
           if (i.gt.natorg.and.natorg.ne.0) dsurf = .true.
           if (iaton(i).ge.1.and.ianz(i).eq.100.and..not.dsurf) then
              n = iconn(1,i)
              if (n.ne.0) then
                 do j=1,n
                    im = iconn(j+1,i)
                    m = abs(im)
                    if (iaton(m).ge.1.and.ianz(m).eq.100) then
                       ic = iatclr(i)
                       if (m.gt.i) then
                           call plvrod(iun,1,jcol,ic,
     &                          coo(1,i),coo(1,m),helrod)
                       else
                           if (n.eq.1) 
     &                     call plvsph(iun,jcol,ic,coo(1,i),helrod)
                       endif
                    endif
                 end do
              endif
           endif
         end do
      endif

      if (iclon.eq.1) then
          ic = 15
          do i=1,4
             ccell = clstr(i)
             do j=1,3
                 tmp(j) = coo(j,iclpnt(i))
             end do
             call pltxt(iun,jcol,ic,tmp,ccell,2)
          end do
      endif

      if (backb.eq.1.and.dolabs.eq.1) then
          ic = 15
          do k=1,nchain
             do j=ianf(k),islu(k)
                jj = j - ianf(k) + 1
                if (reson(j).eq.1) then
                   l = icalf(1,j)
                   do i=1,3
                      tmp(i) = coo(i,l)
                   end do
                   str = '             '
                   str(1:3) = aminos(iamino(j))
                   str(4:4) = ' '
                   write(str(5:8),'(i4)') jj
                   if (nchain.gt.1) then
                       str(9:9) = '.'
                       if (k.lt.10) str(10:10) = char(k+48)
                   endif
                   call pltxt(iun,jcol,ic,tmp,str,1)
                endif
             end do
          end do
      endif

      if (ivtwo.eq.1.and..not.dohead) write(iun,*) ']},'

c
c --- Close out the object, and display it --- FPA 2/27/2000
c
      if (ivtwo .eq. 2) then
        write(iun,*) '}'
        write(iun,*) 'object {molecule}'
      end if
c
      ihand = 0

      return
      end

      subroutine pltxt(iun,jcol,ic,c,str,isize)
      implicit double precision (a-h,o-z)
      common /vropt/ ivtwo,ihand,ivadd
      dimension jcol(3,16), c(3)
      character*(*) str

      if (ivtwo.eq.2) then
      elseif (ivtwo.eq.1) then
         write(iun,*) 'Transform {'
         write(iun,'(''  translation '',3f12.5)') (c(j),j=1,3)
         write(iun,*) '  children ['
         write(iun,*) '    Shape {'
         write(iun,*) '      appearance Appearance {'
         call plvcol(iun,jcol,ic,.false.)
         write(iun,*) '      }'
         write(iun,*) '      geometry Text {'
         write(iun,*) '        string ',
     &                char(34),str,char(34)
         write(iun,*) '        fontStyle FontStyle { size ',isize,' }'
         write(iun,*) '      }'
         write(iun,*) '    }'
         write(iun,*) '  ]'
         write(iun,*) '}'
      else
         write(iun,*) 'Separator {'
         call plvcol(iun,jcol,ic,.false.)
         write(iun,*) '   Transform {'
         write(iun,'(''      translation '',3f12.5)') (c(j),j=1,3)
         write(iun,*) '   }'
         write(iun,*) '   AsciiText { string ',
     &                char(34),str,char(34),' }'
         write(iun,*) '}'
      endif

      return
      end

      subroutine plvcol(iun,jcol,ic,emit)
      implicit double precision (a-h,o-z)
      logical emit
      common /vropt/ ivtwo,ihand,ivadd
      dimension jcol(3,16)

      if (ivtwo.eq.2) then

c        write(iun,*) '  texture {'
c        write(iun,
c     &  '(''pigment { color <'',f12.5,'','',f12.5,'','',f12.5,''> }'')')
c     &              (float(jcol(l,ic))/255,l=1,3)
c        write(iun,*) 
c     &  ' finish { ambient BSAMBI diffuse BSDIFF specular BSSPEC}'
c        write(iun,*) '  }'
         write(iun,*) 'texture { color',char(ic+64),' }'

      elseif (ivtwo.eq.1) then
         write(iun,*) '        material Material {'
      else
         write(iun,*) 'Material {'
      endif

      if (ivtwo.ne.2) then

           write(iun,'(''          diffuseColor  '',3f12.5)') 
     &              (float(jcol(l,ic))/255,l=1,3)
           if ((ivtwo.eq.0.or.ivtwo.eq.1).and.emit) then
              write(iun,'(''          emissiveColor '',3f12.5)') 
     &              (float(jcol(l,ic))/255,l=1,3)
           endif
           write(iun,*) '        }'

      endif


      return
      end

      subroutine plvcoo(iun,t)
      implicit double precision (a-h,o-z)
      dimension t(3)

      write(iun,'(f12.5,'' '',f12.5,'' '',f12.5,'','')') (t(l),l=1,3)

      return
      end

      subroutine plvrod(iun,iaddsp,jcol,ic,p1,p2,rad)
      implicit double precision (a-h,o-z)
      common /vropt/ ivtwo,ihand,ivadd
      dimension p1(3),p2(3),p3(3),v1(3),v2(3),v3(3)
      dimension jcol(3,16)

      do i=1,3
         v1(i) = p2(i) - p1(i)
         v2(i) = 0.0d0
         p3(i) = v1(i)/2.0d0 + p1(i)
      end do
      v2(2) = 1.0d0

      call impsc(v1,v2,cosa)
      if (dabs(cosa).eq.1.0d0) then
         angle = 0.0d0
         do i=1,3
            v3(i) = 0.0d0
         end do
         v3(2) = 1.0d0
      else
         angle = dacos(cosa)
         call crprod(v2,v1,v3)
         vl = vlen(v3)
         do i=1,3
            v3(i) = v3(i) / vl
         end do
      endif

      if (iaddsp.ne.0) call plvsph(iun,jcol,ic,p1,rad*.999d0)

      if (ivtwo.eq.2) then
         write(iun,*) 'cylinder {'
         if (ihand.eq.1) then
         write(iun,
     &   '(''<'',f12.5,'','',f12.5,'','',f12.5,''>, '',
     &     ''<'',f12.5,'','',f12.5,'','',f12.5,''>, '',f12.5)')
     &   -p1(2),p1(1),-p1(3),-p2(2),p2(1),-p2(3),rad
         else
         write(iun,
     &   '(''<'',f12.5,'','',f12.5,'','',f12.5,''>, '',
     &     ''<'',f12.5,'','',f12.5,'','',f12.5,''>, '',f12.5)')
     &   (p1(i),i=1,3),(p2(j),j=1,3),rad
         endif
         call plvcol(iun,jcol,ic,.false.)
         write(iun,*) '}'
      elseif (ivtwo.eq.1) then
         write(iun,*) 'Transform {'
         write(iun,'(''  translation '',3f12.5)') (p3(i),i=1,3)
         write(iun,'(''  rotation '',4f12.5)') (v3(i),i=1,3),angle
         write(iun,*) '  children ['
         write(iun,*) '    Shape {'
         write(iun,*) '      appearance Appearance {'
         call plvcol(iun,jcol,ic,.false.)
         write(iun,*) '      }'
         write(iun,*) '      geometry Cylinder {'
         write(iun,*) '         radius ',real(rad)
         write(iun,*) '         height ',real(vlen(v1))
         write(iun,*) '         top    FALSE'
         write(iun,*) '         bottom FALSE'
         write(iun,*) '      }'
         write(iun,*) '    }'
         write(iun,*) '  ]'
         write(iun,*) '}'
      elseif (ivtwo.eq.0) then
         write(iun,*) 'Separator {'
         call plvcol(iun,jcol,ic,.false.)
         write(iun,*) '   Transform {'
         write(iun,'(''      translation '',3f12.5)') (p3(i),i=1,3)
         write(iun,'(''      rotation '',4f12.5)') (v3(i),i=1,3),angle
         write(iun,*) '   }'
         write(iun,*) '   Cylinder {'
         write(iun,*) '      parts SIDES'
         write(iun,*) '      radius ',real(rad)
         write(iun,*) '      height ',real(vlen(v1))
         write(iun,*) '   }'
         write(iun,*) '}'
      endif

      return
      end

      subroutine plvsph(iun,jcol,ic,p1,rad)
      implicit double precision (a-h,o-z)
      common /vropt/ ivtwo,ihand,ivadd
      dimension p1(3)
      dimension jcol(3,16)

      if (ivtwo.eq.2) then
         write(iun,'(''sphere { '')')
         if (ihand.eq.1) then
         write(iun,
     &   '(''<'',f12.5,'','',f12.5,'','',f12.5,''>, '',f12.5)')
     &   -p1(2),p1(1),-p1(3),rad
         else
         write(iun,
     &   '(''<'',f12.5,'','',f12.5,'','',f12.5,''>, '',f12.5)')
     &   (p1(i),i=1,3),rad
         endif
         call plvcol(iun,jcol,ic,.false.)
         write(iun,'('' }'')')
      elseif (ivtwo.eq.1) then
         write(iun,*) 'Transform {'
         write(iun,'(''  translation '',3f12.5)') (p1(i),i=1,3)
         write(iun,*) '  children ['
         write(iun,*) '    Shape {'
         write(iun,*) '      appearance Appearance {'
         call plvcol(iun,jcol,ic,.false.)
         write(iun,*) '      }'
         write(iun,'(''      geometry Sphere { radius '',f10.5,'' }'')')
     &                 rad
         write(iun,*) '    }'
         write(iun,*) '  ]'
         write(iun,*) '}'
      elseif (ivtwo.eq.0) then
         write(iun,*) 'Separator {'
         call plvcol(iun,jcol,ic,.false.)
         write(iun,*) '   Transform {'
         write(iun,'(''      translation '',3f12.5)') (p1(i),i=1,3)
         write(iun,*) '   }'
         write(iun,'(''   Sphere { radius '',f10.5,'' }'')') rad
         write(iun,*) '}'
      endif

      return
      end

      subroutine plvsd(iun,jcol,k,atcol,hndexl,
     &                 coo,ianz,iaton,iatclr,iconn)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      parameter (mxcon=10)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /athlp/ iatoms, mxnat
      common /vropt/ ivtwo,ihand,ivadd
      logical hndexl
      integer atcol
      dimension jcol(3,16), tmp(3)
      dimension coo(3,*),ianz(*),iaton(*),iatclr(*),iconn(mxcon+1,*)

      if (atcol.eq.1) then
         ic = icol(k)
         if (k.eq.6) ic = 10
      else
         ic = k
      endif

      if (ivtwo.eq.1) then
         write(iun,*) '    Shape {'
         write(iun,*) '      appearance Appearance {'
         call plvcol(iun,jcol,ic,.true.)
         write(iun,*) '      }'
         write(iun,*) '      geometry IndexedLineSet {'
         write(iun,*) '        coord Coordinate { point ['
      else
         call plvcol(iun,jcol,ic,.true.)
         write(iun,*) 'Coordinate3 { point ['
      endif

c
      do i=1,iatoms
        iaa = ianz(i)
        if (atcol.eq.1) then
           ia = iaa
        else
           ia = iatclr(i)
           if (hndexl.and.iaa.ne.100) ia = 0
        endif
        if (iaton(i).ge.1.and.ia.eq.k) then
         n = iconn(1,i)
         if (n.ne.0) then
            call plvcoo(iun,coo(1,i))
            do j=1,n
              im = iconn(j+1,i)
              m = abs(im)
              if (atcol.eq.1) then
                 ma = ianz(m)
              else
                 ma = iatclr(m)
              endif
              if (iaton(m).ge.1) then
                 ido = 1
                 if (ma.eq.k) then
                    if (m.gt.i) then
                       do l=1,3
                          tmp(l) = coo(l,m)
                       end do
                    else
                       ido = 0
                    endif
                 else
                    do l=1,3
                       tmp(l) = 
     &                (coo(l,m) - coo(l,i))/2.0d0 + coo(l,i)
                    end do
                 endif
                 if (ido.eq.1) call plvcoo(iun,tmp)
              endif
            end do
         endif
        endif
      end do
      write(iun,*) ']}'
c
      if (ivtwo.eq.1) then
         write(iun,*) '        coordIndex ['
      else
         write(iun,*) 'IndexedLineSet { coordIndex ['
      endif

      icnt = 0
      do i=1,iatoms
        iaa = ianz(i)
        if (atcol.eq.1) then
           ia = iaa
        else
           ia = iatclr(i)
           if (hndexl.and.iaa.ne.100) ia = 0
        endif
        if (iaton(i).ge.1.and.ia.eq.k) then
         n = iconn(1,i)
         if (n.ne.0) then
            itmp = icnt
            icnt = icnt + 1
            do j=1,n
              im = iconn(j+1,i)
              m = abs(im)
              if (atcol.eq.1) then
                 ma = ianz(m)
              else
                 ma = iatclr(m)
              endif
              if (iaton(m).ge.1) then
                 ido = 1
                 if (ma.eq.k.and.m.le.i) ido = 0
                 if (ido.eq.1) then
                  write(iun,'(i5,'','',i5,'',-1,'')') 
     &                  itmp,icnt
                  icnt = icnt + 1
                 endif
              endif
            end do
         endif
        endif
      end do
      if (ivtwo.eq.1) then
         write(iun,*) '        ]'
      else
         write(iun,*) ']}'
      endif

      if (ivtwo.eq.1) then
         write(iun,*) '      }'
         write(iun,*) '    }'
      endif

      return
      end

      subroutine plvmol(iun,npts,adjus)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)
      common /atopla/ xsym(numatm),ysym(numatm),zsym(numatm),
     &                isym(numatm)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /vrcol/ jcol(3,16)
      logical opfil
      dimension tmp1(3),tmp2(3),tmp3(3)

      roddef = (0.115d0/0.52917706d0)/r(1)
      if (adjus.eq.1.0d0) roddef = roddef*0.52917706d0

      if (opfil(43,'molden_connect',14,1,1,1)) then

         call messg(7)
         do i=1,natoms
            ia = nat(i)
            ic = icol(ia)
            tmp1(1) = -ysym(i) / r(1)
            tmp1(2) = -xsym(i) / r(1)
            tmp1(3) = -zsym(i) / r(1)
            call plvsph(iun,jcol,ic,tmp1,roddef)
         end do

         do while (.true.)

            read(43,*,end=100) i1,i2
            ia = nat(i1)
            ic = icol(ia)
            tmp1(1) = -ysym(i1) / r(1)
            tmp1(2) = -xsym(i1) / r(1)
            tmp1(3) = -zsym(i1) / r(1)
            if (i2.lt.0) then
               read(43,*,end=100) j1,j2,j3,j4,j5
               ja = nat(j1)
               nc = 0
               tmp2(1) = 0.0d0
               tmp2(2) = 0.0d0
               tmp2(3) = 0.0d0
               if (j1.ne.0) then
                  tmp2(1) = tmp2(1) - ysym(j1) / r(1)
                  tmp2(2) = tmp2(2) - xsym(j1) / r(1)
                  tmp2(3) = tmp2(3) - zsym(j1) / r(1)
                  nc = nc + 1
               endif
               if (j2.ne.0) then
                  tmp2(1) = tmp2(1) - ysym(j2) / r(1)
                  tmp2(2) = tmp2(2) - xsym(j2) / r(1)
                  tmp2(3) = tmp2(3) - zsym(j2) / r(1)
                  nc = nc + 1
               endif
               if (j3.ne.0) then
                  tmp2(1) = tmp2(1) - ysym(j3) / r(1)
                  tmp2(2) = tmp2(2) - xsym(j3) / r(1)
                  tmp2(3) = tmp2(3) - zsym(j3) / r(1)
                  nc = nc + 1
               endif
               if (j4.ne.0) then
                  tmp2(1) = tmp2(1) - ysym(j4) / r(1)
                  tmp2(2) = tmp2(2) - xsym(j4) / r(1)
                  tmp2(3) = tmp2(3) - zsym(j4) / r(1)
                  nc = nc + 1
               endif
               if (j5.ne.0) then
                  tmp2(1) = tmp2(1) - ysym(j5) / r(1)
                  tmp2(2) = tmp2(2) - xsym(j5) / r(1)
                  tmp2(3) = tmp2(3) - zsym(j5) / r(1)
                  nc = nc + 1
               endif
               do i=1,3
                  tmp2(i) = tmp2(i) / dble(nc)
               end do
               call plvsph(iun,jcol,icol(ja),tmp2,roddef)
            else
               ja = nat(i2)
               tmp2(1) = -ysym(i2) / r(1)
               tmp2(2) = -xsym(i2) / r(1)
               tmp2(3) = -zsym(i2) / r(1)
            endif
            if (ja.eq.ia) then
               do l=1,3
                  tmp3(l) = tmp2(l)
               end do
               call plvrod(iun,0,jcol,ic,tmp1,tmp3,roddef)
            else
               do l=1,3
                  tmp3(l) = 
     &           (tmp2(l) - tmp1(l))/2.0d0 + tmp1(l)
               end do
               call plvrod(iun,0,jcol,ic,tmp1,tmp3,roddef)
               ic = icol(ja)
               do l=1,3
                  tmp1(l) = tmp2(l)
               end do
               call plvrod(iun,0,jcol,ic,tmp1,tmp3,roddef)
            endif

         end do
100      close(43)
         return
      endif

      do i=1,natoms
         ia = nat(i)
         ic = icol(ia)
         tmp1(1) = -ysym(i) / r(1)
         tmp1(2) = -xsym(i) / r(1)
         tmp1(3) = -zsym(i) / r(1)
         call plvsph(iun,jcol,ic,tmp1,roddef)
         do j=1,natoms
            ja = nat(j)
            tmp2(1) = -ysym(j) / r(1)
            tmp2(2) = -xsym(j) / r(1)
            tmp2(3) = -zsym(j) / r(1)
            dmaxsq = (vdwr(ia) + vdwr(ja))**2
            dijsq = ((xyz(1,i)-xyz(1,j))*adjus)**2
     &             +((xyz(2,i)-xyz(2,j))*adjus)**2
     &             +((xyz(3,i)-xyz(3,j))*adjus)**2
            if (dijsq.lt.dmaxsq) then
                 ido = 1
                 if (ja.eq.ia) then
                    if (j.gt.i) then
                       do l=1,3
                          tmp3(l) = tmp2(l)
                       end do
                    else
                       ido = 0
                    endif
                 else
                    do l=1,3
                       tmp3(l) = 
     &                (tmp2(l) - tmp1(l))/2.0d0 + tmp1(l)
                    end do
                 endif
                 if (ido.eq.1) then
                        call plvrod(iun,0,jcol,ic,
     &                       tmp1,tmp3,roddef)
                 endif
              endif
         end do
      end do

      return
      end

      subroutine calcgr(ind,jnd,knd,fvals,griddim,xydim,scale,gr)
      implicit double precision (a-h,o-z), integer ( i-n)
      integer griddim(3), xydim
      integer dir(3), dir1(3), dir2(3)
      dimension fvals(*),scale(3),gr(3)

      
      dir(1) = ind
      dir(2) = jnd
      dir(3) = knd

      do i=1,3
        fact = 1.0d0
        do j=1,3
            dir1(j) = dir(j)
            dir2(j) = dir(j)
        end do
        if (dir(i).eq.0) then
            dir2(i) = dir2(i) + 1
        elseif (dir(i).eq.(griddim(i)-1)) then
            dir1(i) = dir1(i) - 1
        else
            dir1(i) = dir1(i) - 1
            dir2(i) = dir2(i) + 1
            fact = 0.5d0
        endif
        f2 = fvals(dir2(1) + dir2(2)*griddim(1) + dir2(3)*xydim)
        f1 = fvals(dir1(1) + dir1(2)*griddim(1) + dir1(3)*xydim)
        gr(i) = fact * (f2 - f1) / scale(i)

c        f3 = fvals(dir(1) + dir(2)*griddim(1) + dir(3)*xydim)
c        if (dir(i).eq.0.or.dir(i).eq.(griddim(i)-1)) then
c           gr(i) = fact * (f2 - f1) / scale(i)
c        else
c           gr(i) = fact * ((f3 - f1) + (f2-f3)) / scale(i)
c        endif
      end do

      return
      end



      subroutine mcubes(nx,ny,nz,fvals,fmap,mapit,origin,cntval1,
     &                  cntval2,rad,orb,dolap,lapdbl,trans,iun)
c>>>
c>>> Added Common Block /vrcol/ to get the colors  FPA 3/07/2000
c>>>
      implicit double precision (a-h,o-z), integer ( i-n)
      integer trans, griddim(3), xydim, casemask(8)
      integer cases
      integer orb,iptr,pptr,wpass,npass,tptr,caseptr,findex,posneg
      logical dolap,lapdbl,cneg,utrian
      character*80 lstr
      common /cubes/ cases(16,256),lines(2,12)
      common /vrcol/ jcol(3,16)
      common /vropt/ ivtwo,ihand,ivadd
      common /grdhlp/ mx3d,mx3d2
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      logical valenc,bonds,ovrlap,atomic,doori,dum
      common /option/ idum,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dum
      common /psa/   psa, tsa, exs, pol, pol2, epmin, epmax, ipsa,
     &               icpsa, idtpsa
c>>>
c --- variables for POV macros --- FPA 3/7/2000
      integer colors
      double precision colrgb
      dimension colrgb(3)
c ---                          ---
c>>>
      dimension index(3), linepts(2)
      dimension fvals(*),fmap(*),origin(3),scale(3),cntvls(2)
      dimension grd(3),xyz(3),txyz(3,3),tgrd(3,3),cval(3),g(3),c(3)
      dimension cubval(8), cbmval(8), cubxyz(8,3), cubgr(3,8),rad(3)
      dimension tc(3,3),col(3)
      data casemask /1,2,4,8,16,32,64,128/
      data lines /1,2,2,3,4,3,1,4,5,6,6,7,8,7,5,8,1,5,2,6,4,8,3,7/
      data ((cases(i,j),i=1,16),j=1,50) /
     & -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,8,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,1,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 1,8,3,9,8,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 1,2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,8,3,1,2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 9,2,11,0,2,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 2,8,3,2,11,8,11,9,8,-1,-1,-1,-1,-1,-1,-1,
     & 3,10,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,10,2,8,10,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 1,9,0,2,3,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 1,10,2,1,9,10,9,8,10,-1,-1,-1,-1,-1,-1,-1,
     & 3,11,1,10,11,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,11,1,0,8,11,8,10,11,-1,-1,-1,-1,-1,-1,-1,
     & 3,9,0,3,10,9,10,11,9,-1,-1,-1,-1,-1,-1,-1,
     & 9,8,11,11,8,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 4,7,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 4,3,0,7,3,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,1,9,8,4,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 4,1,9,4,7,1,7,3,1,-1,-1,-1,-1,-1,-1,-1,
     & 1,2,11,8,4,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 3,4,7,3,0,4,1,2,11,-1,-1,-1,-1,-1,-1,-1,
     & 9,2,11,9,0,2,8,4,7,-1,-1,-1,-1,-1,-1,-1,
     & 2,11,9,2,9,7,2,7,3,7,9,4,-1,-1,-1,-1,
     & 8,4,7,3,10,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 10,4,7,10,2,4,2,0,4,-1,-1,-1,-1,-1,-1,-1,
     & 9,0,1,8,4,7,2,3,10,-1,-1,-1,-1,-1,-1,-1,
     & 4,7,10,9,4,10,9,10,2,9,2,1,-1,-1,-1,-1,
     & 3,11,1,3,10,11,7,8,4,-1,-1,-1,-1,-1,-1,-1,
     & 1,10,11,1,4,10,1,0,4,7,10,4,-1,-1,-1,-1,
     & 4,7,8,9,0,10,9,10,11,10,0,3,-1,-1,-1,-1,
     & 4,7,10,4,10,9,9,10,11,-1,-1,-1,-1,-1,-1,-1,
     & 9,5,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 9,5,4,0,8,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,5,4,1,5,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 8,5,4,8,3,5,3,1,5,-1,-1,-1,-1,-1,-1,-1,
     & 1,2,11,9,5,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 3,0,8,1,2,11,4,9,5,-1,-1,-1,-1,-1,-1,-1,
     & 5,2,11,5,4,2,4,0,2,-1,-1,-1,-1,-1,-1,-1,
     & 2,11,5,3,2,5,3,5,4,3,4,8,-1,-1,-1,-1,
     & 9,5,4,2,3,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,10,2,0,8,10,4,9,5,-1,-1,-1,-1,-1,-1,-1,
     & 0,5,4,0,1,5,2,3,10,-1,-1,-1,-1,-1,-1,-1,
     & 2,1,5,2,5,8,2,8,10,4,8,5,-1,-1,-1,-1,
     & 11,3,10,11,1,3,9,5,4,-1,-1,-1,-1,-1,-1,-1,
     & 4,9,5,0,8,1,8,11,1,8,10,11,-1,-1,-1,-1,
     & 5,4,0,5,0,10,5,10,11,10,0,3,-1,-1,-1,-1,
     & 5,4,8,5,8,11,11,8,10,-1,-1,-1,-1,-1,-1,-1,
     & 9,7,8,5,7,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 9,3,0,9,5,3,5,7,3,-1,-1,-1,-1,-1,-1,-1/
      data ((cases(i,j),i=1,16),j=51,100) /
     & 0,7,8,0,1,7,1,5,7,-1,-1,-1,-1,-1,-1,-1,
     & 1,5,3,3,5,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 9,7,8,9,5,7,11,1,2,-1,-1,-1,-1,-1,-1,-1,
     & 11,1,2,9,5,0,5,3,0,5,7,3,-1,-1,-1,-1,
     & 8,0,2,8,2,5,8,5,7,11,5,2,-1,-1,-1,-1,
     & 2,11,5,2,5,3,3,5,7,-1,-1,-1,-1,-1,-1,-1,
     & 7,9,5,7,8,9,3,10,2,-1,-1,-1,-1,-1,-1,-1,
     & 9,5,7,9,7,2,9,2,0,2,7,10,-1,-1,-1,-1,
     & 2,3,10,0,1,8,1,7,8,1,5,7,-1,-1,-1,-1,
     & 10,2,1,10,1,7,7,1,5,-1,-1,-1,-1,-1,-1,-1,
     & 9,5,8,8,5,7,11,1,3,11,3,10,-1,-1,-1,-1,
     & 5,7,0,5,0,9,7,10,0,1,0,11,10,11,0,-1,
     & 10,11,0,10,0,3,11,5,0,8,0,7,5,7,0,-1,
     & 10,11,5,7,10,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 11,6,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,8,3,5,11,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 9,0,1,5,11,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 1,8,3,1,9,8,5,11,6,-1,-1,-1,-1,-1,-1,-1,
     & 1,6,5,2,6,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 1,6,5,1,2,6,3,0,8,-1,-1,-1,-1,-1,-1,-1,
     & 9,6,5,9,0,6,0,2,6,-1,-1,-1,-1,-1,-1,-1,
     & 5,9,8,5,8,2,5,2,6,3,2,8,-1,-1,-1,-1,
     & 2,3,10,11,6,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 10,0,8,10,2,0,11,6,5,-1,-1,-1,-1,-1,-1,-1,
     & 0,1,9,2,3,10,5,11,6,-1,-1,-1,-1,-1,-1,-1,
     & 5,11,6,1,9,2,9,10,2,9,8,10,-1,-1,-1,-1,
     & 6,3,10,6,5,3,5,1,3,-1,-1,-1,-1,-1,-1,-1,
     & 0,8,10,0,10,5,0,5,1,5,10,6,-1,-1,-1,-1,
     & 3,10,6,0,3,6,0,6,5,0,5,9,-1,-1,-1,-1,
     & 6,5,9,6,9,10,10,9,8,-1,-1,-1,-1,-1,-1,-1,
     & 5,11,6,4,7,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 4,3,0,4,7,3,6,5,11,-1,-1,-1,-1,-1,-1,-1,
     & 1,9,0,5,11,6,8,4,7,-1,-1,-1,-1,-1,-1,-1,
     & 11,6,5,1,9,7,1,7,3,7,9,4,-1,-1,-1,-1,
     & 6,1,2,6,5,1,4,7,8,-1,-1,-1,-1,-1,-1,-1,
     & 1,2,5,5,2,6,3,0,4,3,4,7,-1,-1,-1,-1,
     & 8,4,7,9,0,5,0,6,5,0,2,6,-1,-1,-1,-1,
     & 7,3,9,7,9,4,3,2,9,5,9,6,2,6,9,-1,
     & 3,10,2,7,8,4,11,6,5,-1,-1,-1,-1,-1,-1,-1,
     & 5,11,6,4,7,2,4,2,0,2,7,10,-1,-1,-1,-1,
     & 0,1,9,4,7,8,2,3,10,5,11,6,-1,-1,-1,-1,
     & 9,2,1,9,10,2,9,4,10,7,10,4,5,11,6,-1,
     & 8,4,7,3,10,5,3,5,1,5,10,6,-1,-1,-1,-1,
     & 5,1,10,5,10,6,1,0,10,7,10,4,0,4,10,-1,
     & 0,5,9,0,6,5,0,3,6,10,6,3,8,4,7,-1,
     & 6,5,9,6,9,10,4,7,9,7,10,9,-1,-1,-1,-1,
     & 11,4,9,6,4,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 4,11,6,4,9,11,0,8,3,-1,-1,-1,-1,-1,-1,-1,
     & 11,0,1,11,6,0,6,4,0,-1,-1,-1,-1,-1,-1,-1,
     & 8,3,1,8,1,6,8,6,4,6,1,11,-1,-1,-1,-1/
      data ((cases(i,j),i=1,16),j=101,150) /
     & 1,4,9,1,2,4,2,6,4,-1,-1,-1,-1,-1,-1,-1,
     & 3,0,8,1,2,9,2,4,9,2,6,4,-1,-1,-1,-1,
     & 0,2,4,4,2,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 8,3,2,8,2,4,4,2,6,-1,-1,-1,-1,-1,-1,-1,
     & 11,4,9,11,6,4,10,2,3,-1,-1,-1,-1,-1,-1,-1,
     & 0,8,2,2,8,10,4,9,11,4,11,6,-1,-1,-1,-1,
     & 3,10,2,0,1,6,0,6,4,6,1,11,-1,-1,-1,-1,
     & 6,4,1,6,1,11,4,8,1,2,1,10,8,10,1,-1,
     & 9,6,4,9,3,6,9,1,3,10,6,3,-1,-1,-1,-1,
     & 8,10,1,8,1,0,10,6,1,9,1,4,6,4,1,-1,
     & 3,10,6,3,6,0,0,6,4,-1,-1,-1,-1,-1,-1,-1,
     & 6,4,8,10,6,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 7,11,6,7,8,11,8,9,11,-1,-1,-1,-1,-1,-1,-1,
     & 0,7,3,0,11,7,0,9,11,6,7,11,-1,-1,-1,-1,
     & 11,6,7,1,11,7,1,7,8,1,8,0,-1,-1,-1,-1,
     & 11,6,7,11,7,1,1,7,3,-1,-1,-1,-1,-1,-1,-1,
     & 1,2,6,1,6,8,1,8,9,8,6,7,-1,-1,-1,-1,
     & 2,6,9,2,9,1,6,7,9,0,9,3,7,3,9,-1,
     & 7,8,0,7,0,6,6,0,2,-1,-1,-1,-1,-1,-1,-1,
     & 7,3,2,6,7,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 2,3,10,11,6,8,11,8,9,8,6,7,-1,-1,-1,-1,
     & 2,0,7,2,7,10,0,9,7,6,7,11,9,11,7,-1,
     & 1,8,0,1,7,8,1,11,7,6,7,11,2,3,10,-1,
     & 10,2,1,10,1,7,11,6,1,6,7,1,-1,-1,-1,-1,
     & 8,9,6,8,6,7,9,1,6,10,6,3,1,3,6,-1,
     & 0,9,1,10,6,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 7,8,0,7,0,6,3,10,0,10,6,0,-1,-1,-1,-1,
     & 7,10,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 7,6,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 3,0,8,10,7,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,1,9,10,7,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 8,1,9,8,3,1,10,7,6,-1,-1,-1,-1,-1,-1,-1,
     & 11,1,2,6,10,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 1,2,11,3,0,8,6,10,7,-1,-1,-1,-1,-1,-1,-1,
     & 2,9,0,2,11,9,6,10,7,-1,-1,-1,-1,-1,-1,-1,
     & 6,10,7,2,11,3,11,8,3,11,9,8,-1,-1,-1,-1,
     & 7,2,3,6,2,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 7,0,8,7,6,0,6,2,0,-1,-1,-1,-1,-1,-1,-1,
     & 2,7,6,2,3,7,0,1,9,-1,-1,-1,-1,-1,-1,-1,
     & 1,6,2,1,8,6,1,9,8,8,7,6,-1,-1,-1,-1,
     & 11,7,6,11,1,7,1,3,7,-1,-1,-1,-1,-1,-1,-1,
     & 11,7,6,1,7,11,1,8,7,1,0,8,-1,-1,-1,-1,
     & 0,3,7,0,7,11,0,11,9,6,11,7,-1,-1,-1,-1,
     & 7,6,11,7,11,8,8,11,9,-1,-1,-1,-1,-1,-1,-1,
     & 6,8,4,10,8,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 3,6,10,3,0,6,0,4,6,-1,-1,-1,-1,-1,-1,-1,
     & 8,6,10,8,4,6,9,0,1,-1,-1,-1,-1,-1,-1,-1,
     & 9,4,6,9,6,3,9,3,1,10,3,6,-1,-1,-1,-1,
     & 6,8,4,6,10,8,2,11,1,-1,-1,-1,-1,-1,-1,-1,
     & 1,2,11,3,0,10,0,6,10,0,4,6,-1,-1,-1,-1/
      data ((cases(i,j),i=1,16),j=151,200) /
     & 4,10,8,4,6,10,0,2,9,2,11,9,-1,-1,-1,-1,
     & 11,9,3,11,3,2,9,4,3,10,3,6,4,6,3,-1,
     & 8,2,3,8,4,2,4,6,2,-1,-1,-1,-1,-1,-1,-1,
     & 0,4,2,4,6,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 1,9,0,2,3,4,2,4,6,4,3,8,-1,-1,-1,-1,
     & 1,9,4,1,4,2,2,4,6,-1,-1,-1,-1,-1,-1,-1,
     & 8,1,3,8,6,1,8,4,6,6,11,1,-1,-1,-1,-1,
     & 11,1,0,11,0,6,6,0,4,-1,-1,-1,-1,-1,-1,-1,
     & 4,6,3,4,3,8,6,11,3,0,3,9,11,9,3,-1,
     & 11,9,4,6,11,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 4,9,5,7,6,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,8,3,4,9,5,10,7,6,-1,-1,-1,-1,-1,-1,-1,
     & 5,0,1,5,4,0,7,6,10,-1,-1,-1,-1,-1,-1,-1,
     & 10,7,6,8,3,4,3,5,4,3,1,5,-1,-1,-1,-1,
     & 9,5,4,11,1,2,7,6,10,-1,-1,-1,-1,-1,-1,-1,
     & 6,10,7,1,2,11,0,8,3,4,9,5,-1,-1,-1,-1,
     & 7,6,10,5,4,11,4,2,11,4,0,2,-1,-1,-1,-1,
     & 3,4,8,3,5,4,3,2,5,11,5,2,10,7,6,-1,
     & 7,2,3,7,6,2,5,4,9,-1,-1,-1,-1,-1,-1,-1,
     & 9,5,4,0,8,6,0,6,2,6,8,7,-1,-1,-1,-1,
     & 3,6,2,3,7,6,1,5,0,5,4,0,-1,-1,-1,-1,
     & 6,2,8,6,8,7,2,1,8,4,8,5,1,5,8,-1,
     & 9,5,4,11,1,6,1,7,6,1,3,7,-1,-1,-1,-1,
     & 1,6,11,1,7,6,1,0,7,8,7,0,9,5,4,-1,
     & 4,0,11,4,11,5,0,3,11,6,11,7,3,7,11,-1,
     & 7,6,11,7,11,8,5,4,11,4,8,11,-1,-1,-1,-1,
     & 6,9,5,6,10,9,10,8,9,-1,-1,-1,-1,-1,-1,-1,
     & 3,6,10,0,6,3,0,5,6,0,9,5,-1,-1,-1,-1,
     & 0,10,8,0,5,10,0,1,5,5,6,10,-1,-1,-1,-1,
     & 6,10,3,6,3,5,5,3,1,-1,-1,-1,-1,-1,-1,-1,
     & 1,2,11,9,5,10,9,10,8,10,5,6,-1,-1,-1,-1,
     & 0,10,3,0,6,10,0,9,6,5,6,9,1,2,11,-1,
     & 10,8,5,10,5,6,8,0,5,11,5,2,0,2,5,-1,
     & 6,10,3,6,3,5,2,11,3,11,5,3,-1,-1,-1,-1,
     & 5,8,9,5,2,8,5,6,2,3,8,2,-1,-1,-1,-1,
     & 9,5,6,9,6,0,0,6,2,-1,-1,-1,-1,-1,-1,-1,
     & 1,5,8,1,8,0,5,6,8,3,8,2,6,2,8,-1,
     & 1,5,6,2,1,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 1,3,6,1,6,11,3,8,6,5,6,9,8,9,6,-1,
     & 11,1,0,11,0,6,9,5,0,5,6,0,-1,-1,-1,-1,
     & 0,3,8,5,6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 11,5,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 10,5,11,7,5,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 10,5,11,10,7,5,8,3,0,-1,-1,-1,-1,-1,-1,-1,
     & 5,10,7,5,11,10,1,9,0,-1,-1,-1,-1,-1,-1,-1,
     & 11,7,5,11,10,7,9,8,1,8,3,1,-1,-1,-1,-1,
     & 10,1,2,10,7,1,7,5,1,-1,-1,-1,-1,-1,-1,-1,
     & 0,8,3,1,2,7,1,7,5,7,2,10,-1,-1,-1,-1,
     & 9,7,5,9,2,7,9,0,2,2,10,7,-1,-1,-1,-1,
     & 7,5,2,7,2,10,5,9,2,3,2,8,9,8,2,-1/
      data ((cases(i,j),i=1,16),j=201,250) /
     & 2,5,11,2,3,5,3,7,5,-1,-1,-1,-1,-1,-1,-1,
     & 8,2,0,8,5,2,8,7,5,11,2,5,-1,-1,-1,-1,
     & 9,0,1,5,11,3,5,3,7,3,11,2,-1,-1,-1,-1,
     & 9,8,2,9,2,1,8,7,2,11,2,5,7,5,2,-1,
     & 1,3,5,3,7,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,8,7,0,7,1,1,7,5,-1,-1,-1,-1,-1,-1,-1,
     & 9,0,3,9,3,5,5,3,7,-1,-1,-1,-1,-1,-1,-1,
     & 9,8,7,5,9,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 5,8,4,5,11,8,11,10,8,-1,-1,-1,-1,-1,-1,-1,
     & 5,0,4,5,10,0,5,11,10,10,3,0,-1,-1,-1,-1,
     & 0,1,9,8,4,11,8,11,10,11,4,5,-1,-1,-1,-1,
     & 11,10,4,11,4,5,10,3,4,9,4,1,3,1,4,-1,
     & 2,5,1,2,8,5,2,10,8,4,5,8,-1,-1,-1,-1,
     & 0,4,10,0,10,3,4,5,10,2,10,1,5,1,10,-1,
     & 0,2,5,0,5,9,2,10,5,4,5,8,10,8,5,-1,
     & 9,4,5,2,10,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 2,5,11,3,5,2,3,4,5,3,8,4,-1,-1,-1,-1,
     & 5,11,2,5,2,4,4,2,0,-1,-1,-1,-1,-1,-1,-1,
     & 3,11,2,3,5,11,3,8,5,4,5,8,0,1,9,-1,
     & 5,11,2,5,2,4,1,9,2,9,4,2,-1,-1,-1,-1,
     & 8,4,5,8,5,3,3,5,1,-1,-1,-1,-1,-1,-1,-1,
     & 0,4,5,1,0,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 8,4,5,8,5,3,9,0,5,0,3,5,-1,-1,-1,-1,
     & 9,4,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 4,10,7,4,9,10,9,11,10,-1,-1,-1,-1,-1,-1,-1,
     & 0,8,3,4,9,7,9,10,7,9,11,10,-1,-1,-1,-1,
     & 1,11,10,1,10,4,1,4,0,7,4,10,-1,-1,-1,-1,
     & 3,1,4,3,4,8,1,11,4,7,4,10,11,10,4,-1,
     & 4,10,7,9,10,4,9,2,10,9,1,2,-1,-1,-1,-1,
     & 9,7,4,9,10,7,9,1,10,2,10,1,0,8,3,-1,
     & 10,7,4,10,4,2,2,4,0,-1,-1,-1,-1,-1,-1,-1,
     & 10,7,4,10,4,2,8,3,4,3,2,4,-1,-1,-1,-1,
     & 2,9,11,2,7,9,2,3,7,7,4,9,-1,-1,-1,-1,
     & 9,11,7,9,7,4,11,2,7,8,7,0,2,0,7,-1,
     & 3,7,11,3,11,2,7,4,11,1,11,0,4,0,11,-1,
     & 1,11,2,8,7,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 4,9,1,4,1,7,7,1,3,-1,-1,-1,-1,-1,-1,-1,
     & 4,9,1,4,1,7,0,8,1,8,7,1,-1,-1,-1,-1,
     & 4,0,3,7,4,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 4,8,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 9,11,8,11,10,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 3,0,9,3,9,10,10,9,11,-1,-1,-1,-1,-1,-1,-1,
     & 0,1,11,0,11,8,8,11,10,-1,-1,-1,-1,-1,-1,-1,
     & 3,1,11,10,3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 1,2,10,1,10,9,9,10,8,-1,-1,-1,-1,-1,-1,-1,
     & 3,0,9,3,9,10,1,2,9,2,10,9,-1,-1,-1,-1,
     & 0,2,10,8,0,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 3,2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 2,3,8,2,8,11,11,8,9,-1,-1,-1,-1,-1,-1,-1,
     & 9,11,2,0,9,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/
      data ((cases(i,j),i=1,16),j=251,256) /
     & 2,3,8,2,8,11,0,1,8,1,11,8,-1,-1,-1,-1,
     & 1,11,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 1,3,8,9,1,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,9,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & 0,3,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     & -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/
  
      emax  = -1000000.0d0
      emin  =  1000000.0d0
      epmax = -1000000.0d0
      epmin =  1000000.0d0
      psa   = 0.0d0
      tsa   = 0.0d0
      exs   = 0.0d0

c ityp
c  0   orbital
c  1   normal density
c  2   laplacian
c  3   mapped electrostatic potential
c  4   electrostatic potential
c  5   difference density
c  6   unknown

      ityp = 0
      if (orb.eq.0) then
         if (trans.eq.1) then
            ityp = 1
         elseif (dolap) then
            ityp = 2
         elseif (mapit.eq.1) then
            ityp = 3
         elseif (molpot.or.elpot.or.chpot) then
            ityp = 4
         elseif (bonds) then
            ityp = 5
         else
            ityp = 6
         endif
      else
         ityp = 0
      endif

      icolver = 0
      if (mapit.eq.1.and.ivtwo.eq.1) then
         open(unit=47,form='formatted',status='scratch',err=10)
         icolver = 1
      endif
10    continue
      
      griddim(1) = nx
      griddim(2) = ny
      griddim(3) = nz

c rad(2) slaat op scale(1) en griddim(1)

      do i=1,3
          scale(i) = 1.0d0/(griddim(i)-1.0d0)
      end do
      scale(1) = scale(1) * rad(2) / rad(1)
      scale(3) = scale(3) * rad(3) / rad(1)

      cntvls(1) = cntval1
      cntvls(2) = cntval2

      npass = 3
      if (trans.eq.1.or.(dolap.and..not.lapdbl).or.
     &    mapit.eq.1) npass = 2
      if (ivtwo.eq.4) then
          npass = 3
          if (trans.eq.1.or.(dolap.and..not.lapdbl).or.
     &        mapit.eq.1) npass = 1
      endif

      do wpass=0,npass
         if (wpass.eq.0) then
             if (ivtwo.eq.3.or.ivtwo.eq.4) then
             elseif (ivtwo.eq.2) then
                write(iun,'(''camera {'')')
                write(iun,*) 'location <0.0, 0.0, 3.0>'
                write(iun,*) 'look_at <0.0, 0.0, 0.0>'
                write(iun,*) 'angle 15'
                write(iun,'(''}'')')
                write(iun,*) 
     &            'light_source { <9, 9, 9> color rgb<1, 1, 1> }'
                write(iun,*) '#declare MeshGold = texture {'
                write(iun,*) 
     &              'pigment { color rgbf<0.96, 0.82, 0.65, 0.7> }'
                write(iun,*) 
     &             'finish { ambient 0.4 diffuse 0.4 specular 0.9}'
                write(iun,*) '}'
                write(iun,*) '#declare MeshRed = texture {'
                write(iun,*) 
     &             'pigment { color rgbf<0.8, 0.2, 0.2, 0.7> }'
                write(iun,*) 
     &             'finish { ambient 0.2 diffuse 0.6 specular 0.9}'
                write(iun,*) '}'
                write(iun,*) '#declare MeshBlue = texture {'
                write(iun,*) 
     &             'pigment { color rgbf<0.2, 0.2, 0.8, 0.7> }'
                write(iun,*) 
     &             'finish { ambient 0.2 diffuse 0.6 specular 0.9}'
                write(iun,*) '}'
                write(iun,*) '#declare BSAMBI = 0.2;'
                write(iun,*) '#declare BSDIFF = 0.8;'
                write(iun,*) '#declare BSSPEC = 0.8;'
c>>>
c>>>  Redefine the Colors for the Atoms, etc, here.
c>>>  FPA 3/07/2000
c>>>
                do colors = 1,15
c
                    write(iun,'(a14,a1,a3)') '#declare color',
     &                                       char(colors+64),' = '
                    write(iun,*) 'texture { '
                    do n = 1,3
                       colrgb(n) = (dble(jcol(n,colors))/255.0d0)
                    end do
                    write(iun,'(a15,3f8.4,a3)') ' pigment { rgb<',
     &                         (colrgb(n),n=1,3),' >}'
                    write(iun,*)
     &        'finish {ambient BSAMBI diffuse BSDIFF specular BSSPEC}'
                    write(iun,*)'}'
                 end do
c
c --- Declare the molecule object ---- FPA 2/27/2000
c
                 write(iun,*) '# declare molecule = union { '
c>>
c>>
             elseif (ivtwo.eq.1) then
                write(iun,'(''#VRML V2.0 utf8'')')
                write(iun,*) 'NavigationInfo { type ',char(34),
     &                       'EXAMINE',char(34),' }'
                write(iun,*) 'Viewpoint { position 0.0 0.0 1.0 ',
     &                       'description ',char(34),'MoldenVRML',
     &                        char(34),' }'
                if (ivadd.eq.1) then
                   write(iun,*) 'Background {'
                   write(iun,*) 
     &             'skyColor [1.0 1.0 0.9,0.5 0.5 0.4,0.2 0.2 0.1]'
                   write(iun,*)
     &             'skyAngle [1.571,2.771]}'
                   write(iun,*) 'PointLight { on TRUE'
                   write(iun,*) '             location 0.0 0.5 0.5'
                   write(iun,*) '             radius 100.0'
                   write(iun,*) '             intensity 1.0'
                   write(iun,*) '             ambientIntensity 0.3'
                   write(iun,*) '             color 1.0 1.0 1.0'
                   write(iun,*) '             attenuation 1.0 0.0 0.0'
                   write(iun,*) '           }'
                endif

                write(iun,*) 'WorldInfo {'
                write(iun,*) '  title ',char(34),'MoldenVRML',char(34)
                write(iun,*) '  info ',char(34)
             else
                write(iun,'(''#VRML V1.0 ascii'')')
                write(iun,*) 'Separator {'
                write(iun,*) 'Info { string ',char(34),
     &                     'Created by Molden'
             endif
             if (ivtwo.lt.2) then
             if (orb.eq.0) then
                if (trans.eq.1) then
                   write(iun,*) '               Density'
                else
                   if (dolap) then
                      write(iun,*) 
     &                '               Laplacian'

                   else
                      write(iun,*) 
     &                '               Molecular Difference Density'
                   endif
                endif
             else
                write(iun,'(''               orbital '',i4)') orb
             endif
             write(iun,'(''               Contour Value: '',f12.5)') 
     &             cntval1
             write(iun,'(''               Grid Size    : '',f12.5)') 
     &             rad(1)
             write(iun,*) char(34),' }'
             endif

             if (ivtwo.eq.3.or.ivtwo.eq.4) then
                i0 = 1
                if (mapit.eq.1) i0 = 3
                if (trans.eq.1) i0 = -1*i0
                if (ivtwo.eq.3) then
                   call ogbeg(i0,orb,cntval1)
                else
                   call ogbegg(i0,i0,ityp,orb,cntval1,mapit,lstr)
                endif
             elseif (ivtwo.eq.2) then
                write(iun,*) 'mesh {'
             elseif (ivtwo.eq.1) then
                write(iun,*) '   PROTO Surface ['
                write(iun,*) '   field SFFloat aInt 0.0'
                write(iun,*) '   field SFColor dCol 1.0 0.1 0.1'
                write(iun,*) '   field SFColor sCol 0.8 0.8 0.8'
                write(iun,*) '   field SFFloat shi 0.75'
                write(iun,*) '   field SFFloat trans 0.0'
                write(iun,*) '   field MFInt32 CIndex [ ]'
                write(iun,*) '   ]'
                write(iun,*) '   {'
                write(iun,*) '    Group {'
                write(iun,*) '    children ['
                write(iun,*) '    DEF Trigger TouchSensor { }'
                write(iun,*) '    Shape {'
                write(iun,*) '      appearance Appearance {'
                write(iun,*) '       material DEF CT Material {'
                write(iun,*) '        ambientIntensity IS aInt'
                write(iun,*) '        diffuseColor IS dCol'
                write(iun,*) '        specularColor IS sCol'
                write(iun,*) '        shininess IS shi'
                write(iun,*) '        transparency IS trans'
                write(iun,*) '       }'
                write(iun,*) '      }'
                write(iun,*) '      geometry IndexedFaceSet {'
                write(iun,*) '        coord Coordinate { point [ '
             else
                write(iun,*) 
     &              'NormalBinding { value PER_VERTEX_INDEXED }'
                write(iun,*) 'Coordinate3 { point [ '
             endif
         endif

         if (wpass.eq.1) then
             if (ivtwo.eq.2.or.ivtwo.eq.3.or.ivtwo.eq.4) then
             elseif (ivtwo.eq.1) then
                write(iun,*) '        normalPerVertex TRUE'
                write(iun,*) '        normal Normal { vector [ '
                if (icolver.eq.1.and.ivtwo.eq.1) then
                  write(47,*) '        colorPerVertex TRUE'
                  write(47,*) '        color Color { color [ '
                endif
             else
                write(iun,*) 'Normal { vector [ '
             endif
         endif

         if (wpass.eq.2) then
             if (ivtwo.eq.3.or.ivtwo.eq.4) then
             elseif (ivtwo.eq.2) then
             elseif (ivtwo.eq.1) then
                write(iun,*) '   Surface {'
                if (orb.eq.0.and..not.dolap) then
                   write(iun,*) '     aInt 0.0'
                   write(iun,*) '     dCol 1.0 1.0 0.0'
                   write(iun,*) '     sCol 0.9 0.8 0.8'
                   write(iun,*) '     shi 0.55'
                   if (trans.eq.1) write(iun,*) '     trans 0.5'
                else
                   if (ivadd.eq.1) then
                      write(iun,*) '     aInt 0.7'
                   else
                      write(iun,*) '     aInt 0.0'
                   endif
                   write(iun,*) '     dCol 0.15 0.15 1.0'
                   write(iun,*) '     sCol 0.8 0.8 0.8'
                   write(iun,*) '     shi 0.75'
                   if (trans.eq.1) write(iun,*) '     trans 0.6'
                endif
                write(iun,*) '     CIndex ['
             else
                call plvmat(iun,1,orb,trans,dolap)
                write(iun,*) 'IndexedFaceSet { coordIndex [ '
             endif
         endif

         if (wpass.eq.3) then
             if (ivtwo.eq.3.or.ivtwo.eq.4) then
                i1 = 2
                if (ivtwo.eq.3) then
                   call ogbeg(i1,orb,cntval2)
                else
                   call ogbegg(i1,i1,ityp,orb,cntval2,mapit,lstr)
                endif
             elseif (ivtwo.eq.2) then
                write(iun,*) 'mesh {'
             elseif (ivtwo.eq.1) then
                write(iun,*) '   Surface {'
                if (ivadd.eq.1) then
                   write(iun,*) '     aInt 0.7'
                else
                   write(iun,*) '     aInt 0.0'
                endif
                write(iun,*) '     dCol 1.0 0.1 0.1'
                write(iun,*) '     sCol 0.8 0.8 0.8'
                write(iun,*) '     shi 0.75'
                if (trans.eq.1) write(iun,*) '     trans 0.6'
                write(iun,*) '     CIndex ['
             else
                call plvmat(iun,2,orb,trans,dolap)
                write(iun,*) 'IndexedFaceSet { coordIndex [ '
             endif
         endif

         iptr = 0
         pptr = 0
c         xydim = griddim(1) * griddim(2)
         xydim = mx3d2
         do posneg=1,2
            contval = cntvls(posneg)
            cneg = .false.
            if (contval.lt.0.0d0) cneg = .true.
      

            do k=0,griddim(3)-2
               cubxyz(1,3) = origin(3) + k*scale(3)
               ztmp = origin(3) + (k+1)*scale(3)

               do j=0,griddim(2)-2
                  cubxyz(1,2) = origin(2) + j*scale(2)
                  ytmp = origin(2) + (j+1)*scale(2)

                  do i=0,griddim(1)-2
                     findex = i + j*griddim(1) + k*xydim + 1

c                    create Cube Function Values

                     cubval(1) = fvals(findex)
                     cubval(2) = fvals(findex+1)
                     cubval(3) = fvals(findex+1 + griddim(1))
                     cubval(4) = fvals(findex + griddim(1))
                     cubval(5) = fvals(findex + xydim)
                     cubval(6) = fvals(findex+1 + xydim)
                     cubval(7) = fvals(findex+1 + griddim(1) + xydim)
                     cubval(8) = fvals(findex + griddim(1) + xydim)

                     if (mapit.eq.1) then
                        cbmval(1) = fmap(findex)
                        cbmval(2) = fmap(findex+1)
                        cbmval(3) = fmap(findex+1 + griddim(1))
                        cbmval(4) = fmap(findex + griddim(1))
                        cbmval(5) = fmap(findex + xydim)
                        cbmval(6) = fmap(findex+1 + xydim)
                        cbmval(7) = fmap(findex+1 + griddim(1) + xydim)
                        cbmval(8) = fmap(findex + griddim(1) + xydim)
                     endif

c                    which of the 256 cases is it ?

                     caseptr = 1
                     do l=1,8
                         if (cneg) then
                            if (cubval(l).ge.contval) 
     &                      caseptr = caseptr + casemask(l)
                         else
                            if (cubval(l).le.contval) 
     &                      caseptr = caseptr + casemask(l)
                         endif
                     end do

                     if (caseptr.eq.1.or.caseptr.eq.256 ) goto 100


c                    create Cube points

                     cubxyz(1,1) = origin(1) + i*scale(1)
                     xtmp = origin(1) + (i+1)*scale(1)

                     cubxyz(2,1) = xtmp
                     cubxyz(2,2) = cubxyz(1,2)
                     cubxyz(2,3) = cubxyz(1,3)

                     cubxyz(3,1) = xtmp
                     cubxyz(3,2) = ytmp
                     cubxyz(3,3) = cubxyz(1,3)

                     cubxyz(4,1) = cubxyz(1,1)
                     cubxyz(4,2) = ytmp
                     cubxyz(4,3) = cubxyz(1,3)

                     cubxyz(5,1) = cubxyz(1,1)
                     cubxyz(5,2) = cubxyz(1,2)
                     cubxyz(5,3) = ztmp

                     cubxyz(6,1) = xtmp
                     cubxyz(6,2) = cubxyz(1,2)
                     cubxyz(6,3) = ztmp

                     cubxyz(7,1) = xtmp
                     cubxyz(7,2) = ytmp
                     cubxyz(7,3) = ztmp

                     cubxyz(8,1) = cubxyz(1,1)
                     cubxyz(8,2) = ytmp
                     cubxyz(8,3) = ztmp

c                    create Cube gradients

                     call calcgr(i,j,k,fvals,griddim,xydim,
     &                           scale,cubgr(1,1))
                     call calcgr(i+1,j,k,fvals,griddim,xydim,
     &                           scale,cubgr(1,2))
                     call calcgr(i+1,j+1,k,fvals,griddim,xydim,
     &                           scale,cubgr(1,3))
                     call calcgr(i,j+1,k,fvals,griddim,xydim,
     &                           scale,cubgr(1,4))
                     call calcgr(i,j,k+1,fvals,griddim,xydim,
     &                           scale,cubgr(1,5))
                     call calcgr(i+1,j,k+1,fvals,griddim,xydim,
     &                           scale,cubgr(1,6))
                     call calcgr(i+1,j+1,k+1,fvals,griddim,xydim,
     &                           scale,cubgr(1,7))
                     call calcgr(i,j+1,k+1,fvals,griddim,xydim,
     &                           scale,cubgr(1,8))

                     tptr = 1
                     do while (cases(tptr,caseptr).gt.-1)
                        do l=0,2
                           linepts(1) = 
     &                       lines(1,cases(tptr+l,caseptr)+1)
                           linepts(2) = 
     &                       lines(2,cases(tptr+l,caseptr)+1)
                           factor = (contval - cubval(linepts(1))) / 
     &                       (cubval(linepts(2)) - cubval(linepts(1)))
                             
                           do m=1,3
                               xyz(m) = cubxyz(linepts(1),m) + factor*
     &                   (cubxyz(linepts(2),m) - cubxyz(linepts(1),m))
                               txyz(m,l+1) = xyz(m)
                               grd(m) = cubgr(m,linepts(1)) + factor*
     &                   (cubgr(m,linepts(2)) - cubgr(m,linepts(1)))
                           end do
                           if (mapit.eq.1) then
                               cval(l+1) = cbmval(linepts(1)) + factor*
     &                        (cbmval(linepts(2)) - cbmval(linepts(1)))
                               if (cval(l+1).gt.emax) emax = cval(l+1)
                               if (cval(l+1).lt.emin) emin = cval(l+1)
                           endif
                           if ((wpass.eq.0).and.ivtwo.lt.2)
     &                        write(iun,'(3f12.5,a)') (xyz(n),n=1,3),','
               
                           index(l+1) = iptr
                           call normal(grd)
                           if (contval.gt.0.0) call arrsgn(grd,3)
                           if (ivtwo.ge.2) then
                              do m=1,3
                                 tgrd(m,l+1) = grd(m)
                              end do
                           endif
                           if (wpass.eq.1.and.ivtwo.lt.2)
     &                        write(iun,'(3f10.6,a)') (grd(n),n=1,3),','
                           if (wpass.eq.1.and.icolver.eq.1) then
                              call parcol(cval(l+1),col)
                              write(47,'(3f7.3,a)') (col(n),n=1,3),','
                           endif
                           iptr = iptr + 1
                        end do

                        if (((wpass.eq.0.and.posneg.eq.1).or.(
     &                        wpass.eq.3.and.posneg.eq.2))
     &                      .and.ivtwo.ge.2) then
                         if (utrian(txyz)) then
                          if (ivtwo.eq.2)  then
                              write(iun,*) 'smooth_triangle {'
                           do l=1,3
                             if (l.eq.3) then
                              write(iun,'(2(a,3(f10.6,a)))') 
     &              '<',txyz(1,l),',',txyz(2,l),',',txyz(3,l),'>, ',
     &              '<',tgrd(1,l),',',tgrd(2,l),',',tgrd(3,l),'>'
                             else
                              write(iun,'(2(a,3(f10.6,a)))') 
     &              '<',txyz(1,l),',',txyz(2,l),',',txyz(3,l),'>, ',
     &              '<',tgrd(1,l),',',tgrd(2,l),',',tgrd(3,l),'>,'
                             endif
                           end do
                           write(iun,*) '}'
                          else
                          endif

                          if (ivtwo.eq.3) then

                           do l=1,3
                            if (mapit.eq.1) then
                               call parcol(cval(l),col)
                               call ogcol(col(1),col(2),col(3))
                            endif
                            call ognorm(tgrd(1,l),tgrd(2,l),tgrd(3,l))
                            call ogvert(txyz(1,l),txyz(2,l),txyz(3,l))
                           end do

                          elseif (ivtwo.eq.4) then

                           do l=1,3
                            if (mapit.eq.1) then
                               call parcol(cval(l),col)
                               call tstpsa(cval,pol,pol2,ir)
                               if (ir.eq.3.and.icpsa.eq.1) then
                                 cav = cval(1)+cval(2)+cval(3)
                                 if (cav.gt.0) then
                                    call ogcoll(0.0d0,0.0d0,1.0d0)
                                 else
                                    call ogcoll(1.0d0,0.0d0,0.0d0)
                                 endif
                               else
                                 call ogcoll(col(1),col(2),col(3))
                               endif
                            endif
                            call rtgbck(tgrd(1,l),tgrd(2,l),tgrd(3,l),g)
                            call ognrm(g(1),g(2),g(3))
                            call rotbck(txyz(1,l),txyz(2,l),txyz(3,l),c)
                            call ogvrt(c(1),c(2),c(3))
                            do m=1,3
                               tc(m,l) = c(m)
                               if (cval(m).gt.epmax) epmax = cval(m)
                               if (cval(m).lt.epmin) epmin = cval(m)
                            end do
                           end do

                           if (mapit.eq.1.and.ipsa.eq.1) then 
                               call calpsa(tc,cval,pol,pol2,
     &                                     psa,tsa,exs,3)
                           endif

                          endif
                         endif
                        endif

                        if (ivtwo.lt.2) then
                          if (wpass.eq.2.and.posneg.eq.1) then
                               write(iun,'(3(i7,'',''),a)') 
     &                         (index(n),n=1,3),'-1,'
                          endif
                          if (wpass.eq.3.and.posneg.eq.2)then
                               write(iun,'(3(i7,'',''),a)') 
     &                         (index(n),n=1,3),'-1,'
                          endif
                        endif
                        pptr = pptr + 1
                        tptr = tptr + 3

                     end do
100               continue
                  end do
               end do
            end do
         end do

         if (ivtwo.eq.3.or.ivtwo.eq.4) then
            if (wpass.eq.0.or.wpass.eq.3) then
                if (ivtwo.eq.3) then
                   call ogend
                else
                   call ogendd(-1)
                endif
            endif
         elseif (ivtwo.eq.2) then
            if (wpass.eq.0) then
               if (npass.eq.2) then
                  write(iun,*) ' texture { MeshGold }'
                  write(iun,*) '}'
               else
                  write(iun,*) ' texture { MeshBlue }'
                  write(iun,*) '}'
               endif
            elseif (wpass.eq.3) then
               write(iun,*) ' texture { MeshRed }'
               write(iun,*) '}'
            endif
         else
            write(iun,*) ']}'
            if (icolver.eq.1.and.wpass.eq.1) then
               rewind(47)
               do while (.true.)
                  read(47,'(a)',end=200) lstr
                  write(iun,'(a)') lstr(1:linlen(lstr))
               end do
200            close(47)
               write(iun,*) ']}'
            endif
         endif

         if (wpass.eq.1.and.ivtwo.eq.1) then
            write(iun,*) '        coordIndex IS CIndex'
            write(iun,*) '  }'
            write(iun,*) ' }'
            write(iun,*) ' ]}'
            write(iun,*) 'DEF TTScript Script {'
            write(iun,*) '        eventIn SFTime  touchTime'
            write(iun,*) '        eventOut        SFFloat Trans'
            write(iun,*) '        field           SFBool  opaque TRUE'
            write(iun,*) '        url'
            write(iun,*) '          ',char(34),'vrmlscript:'
            write(iun,*) '               function touchTime (timeVal) {'
            write(iun,*) '                 if (opaque == true) {'
            write(iun,*) '                        Trans = 0.4;'
            write(iun,*) '                        opaque = false;'
            write(iun,*) '                 } else {'
            write(iun,*) '                        Trans = 0.0;'
            write(iun,*) '                        opaque = true;'
            write(iun,*) '                 }'
            write(iun,*) '                }'
            write(iun,*) '          ',char(34)
            write(iun,*) '}'
            write(iun,*) 'ROUTE Trigger.touchTime TO TTScript.touchTime'
            write(iun,*) 'ROUTE TTScript.Trans TO CT.transparency'
            write(iun,*) '}'
         endif

      end do

      if (ivtwo.eq.0) write(iun,*) '}'

      if (mapit.eq.1) then
         print*,' '
         print*,'Emax = ',emax
         print*,'Emin = ',emin
         print*,' '
         if (ipsa.eq.1) call prtpsa
      endif

      return
      end


      subroutine normal(x)
      implicit double precision (a-h,o-z)
      dimension x(3)

      den = dsqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
      if (den.ne.0.0d0) then
         do i=1,3
            x(i) = x(i) / den
         end do
      endif

      return
      end

      subroutine arrsgn(x,n)
      implicit double precision (a-h,o-z)
      dimension x(3)

      do i=1,3
         x(i) = -x(i)
      end do

      return
      end

      logical function utrian(txyz)
      implicit double precision (a-h,o-z)
      dimension txyz(3,3)

      utrian = .true.

      if (txyz(1,1).eq.txyz(1,2).and.txyz(2,1).eq.txyz(2,2).and.
     & txyz(3,1).eq.txyz(3,2)) utrian = .false.

      if (txyz(1,1).eq.txyz(1,3).and.txyz(2,1).eq.txyz(2,3).and.
     & txyz(3,1).eq.txyz(3,3)) utrian = .false.

      if (txyz(1,2).eq.txyz(1,3).and.txyz(2,2).eq.txyz(2,3).and.
     & txyz(3,2).eq.txyz(3,3)) utrian = .false.

      return
      end

      subroutine epvrmd(vdwr,moddma,natoms,norbs,idops,iaton)
      implicit double precision (a-h,o-z)
      logical valenc,bonds,ovrlap,atomic,doori,denok,dolap,opfil
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap 
      character vrmlfil*256
      common /vrunit/ vrmlfil

      dimension vdwr(*),iaton(*)

      if ((iftyp.ge.2.and.iftyp.le.4).and.norbs.ne.0) then
          iopt = 2
          if (elpot) iopt = 1
          call xyzcoo(1,0,0)

          do i=1,natoms
              iaton(i) = 1
          end do

          call doconn
          ipsi = 0
          dolap = .false.
          call denmak(denok)
          if (iopt.eq.2)
     &          call muldma(vdwr,moddma,.true.,.true.)
          call allsrf(0,iopt,1,4)
          call doscal
          if (idops.eq.1) then
             call plpost(0,0,1,1,1,0,1,0)
             call qupd
          else
             iun = 46
             if (opfil(iun,vrmlfil,linlen(vrmlfil),1,0,0)) then
                call plvrml(iun,1,1,0,1,0,.true.)
                close(iun)
             endif
          endif
      else
          call inferr('No Electrostatic Potential !',1)
      endif
      
      stop

      return
      end

      subroutine plvhd(iun)
      implicit double precision (a-h,o-z)

      write(iun,'(''#VRML V2.0 utf8'')')
      write(iun,*) 'NavigationInfo { type ',char(34),
     &             'EXAMINE',char(34),' }'
      write(iun,*) 'DEF V0 Viewpoint { position 0 0 0 }'
      write(iun,*) ' Background {'
      write(iun,*) 'skyColor [0 0.2 0.7,0 0.5 1,1 1 1 ]'
      write(iun,*) 'skyAngle [ 1.309, 1.571 ]'
      write(iun,*) 'groundColor [0.1 0.1 0,0.4 0.25 0.2,0.6 0.6 0.6]'
      write(iun,*) 'groundAngle [1.309,1.571]}'

      write(iun,*) 'DEF ANIME Transform {'
      write(iun,*) '  children ['
      write(iun,*) '    DEF Sensor SphereSensor { },'
      write(iun,*) '    DEF MOLDENS Switch {'
      write(iun,*) '      whichChoice 0'
      write(iun,*) '      choice ['

      return
      end

      subroutine plvedd(iun,loop,epoints,scal)
      implicit double precision (a-h,o-z)
      common /pnthlp/ ipoints,ipnt
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      integer loop
      dimension epoints(*)


      write(iun,*) '      ]'
      write(iun,*) '    }'
      write(iun,*) '  ]'
      write(iun,*) '}'
      write(iun,*) ' '
      if (loop.eq.1) then
         write(iun,*) 'DEF MOLDEN_TIMER TimeSensor {'
         write(iun,*) 'cycleInterval 1 loop TRUE stopTime 1}'
      else
         write(iun,*) 
     &   'DEF MOLDEN_TIMER TimeSensor { cycleInterval 10 }'
      endif
      write(iun,*) ' '

      write(iun,*) 'EXTERNPROTO SwitchInterpolator ['
      write(iun,*) 'eventIn SFFloat set_fraction'
      write(iun,*) 'eventIn SFVec3f hud_translation'
      write(iun,*) 'eventIn SFRotation hud_rotation'
      write(iun,*) 'eventOut SFInt32 value_changed'
      write(iun,*) 'eventOut SFTime touchTime'
      write(iun,*) 'eventOut SFVec3f v0_position'
      write(iun,*) 'field MFFloat key'
      write(iun,*) 'field MFFloat keyValue'
      write(iun,*) 'field MFFloat Energies'
      write(iun,*) 'field MFInt32 Window'
      write(iun,*) 'field SFVec3f vpnt'
      write(iun,*) '] ',char(34),
     &    'http://www.cmbi.kun.nl/~schaft/molden/wrl/animator.wrl',
     &             char(34)
      write(iun,*) ' '


      write(iun,*) 'DEF MOLDEN_ANIMATOR SwitchInterpolator {'
      write(iun,*) '  key [0,1]'
      write(iun,*) '  keyValue [0,',ipoints-1,']'
      write(iun,*) '  Energies ['
      do i=1,ipoints-1
         write(iun,'(f12.6,a)') epoints(i),','
      end do
      write(iun,'(f12.6)') epoints(ipoints)
      write(iun,*) '  ]'
      if (loop.eq.1) then
         write(iun,*) '  Window [0,0]'
      else
         write(iun,*) '  Window [2,1]'
      endif
      write(iun,*) '  vpnt 0 0 ',scal
      write(iun,*) '}'

      write(iun,*) 
     &'DEF P0 ProximitySensor{ center 0 0 0 size 100 100 100 }'
      write(iun,*) 
     &'ROUTE Sensor.rotation_changed TO ANIME.set_rotation'
      write(iun,*) 
     &'ROUTE P0.position_changed TO MOLDEN_ANIMATOR.hud_translation'
      write(iun,*) 
     &'ROUTE P0.orientation_changed TO MOLDEN_ANIMATOR.hud_rotation'
      write(iun,*) 
     &'ROUTE MOLDEN_ANIMATOR.v0_position TO V0.set_position'
      write(iun,*) 
     &'ROUTE MOLDEN_ANIMATOR.touchTime TO MOLDEN_TIMER.startTime'
      write(iun,*) 'ROUTE MOLDEN_TIMER.fraction_changed TO ',
     &             'MOLDEN_ANIMATOR.set_fraction'
      write(iun,*) 
     &'ROUTE MOLDEN_ANIMATOR.value_changed TO MOLDENS.whichChoice'

      return
      end

      subroutine plvmat(iun,iopt,orb,trans,dolap)
      implicit double precision (a-h,o-z)
      integer orb,trans
      logical dolap

      write(iun,*) 'Material {'

      if (iopt.eq.1) then
         if (orb.eq.0.and..not.dolap) then
            write(iun,*) '          ambientColor 1.000 1.000 1.000'
            write(iun,*) '          diffuseColor 1.000 1.000 0.000'
            write(iun,*) '          specularColor 0.900 0.800 0.800'
            write(iun,*) '          shininess 0.550'
         else
            write(iun,*) '          ambientColor 0.150 0.150 1.000'
            write(iun,*) '          diffuseColor 0.150 0.150 1.000'
            write(iun,*) '          specularColor 0.800 0.800 0.800'
            write(iun,*) '          shininess 0.750'
         endif
      else
         write(iun,*) '          ambientColor 1.000 0.000 0.000'
         write(iun,*) '          diffuseColor 1.000 0.100 0.100'
         write(iun,*) '          specularColor 0.800 0.800 0.800'
         write(iun,*) '          shininess 0.750'
      endif

      if (trans.eq.1)
     &     write(iun,*) '          transparency 0.55'
      write(iun,*) '}'

      return
      end

      subroutine parcol(cval,col)
      implicit double precision (a-h,o-z)
      common /mapcls/ valcol(5),colmap(3,5)
      dimension col(3)

      if (cval.lt.valcol(1)) then
         i1 = 1
         i2 = 1
      elseif (cval.lt.valcol(2)) then
         i1 = 1
         i2 = 2
      elseif (cval.lt.valcol(3)) then
         i1 = 2
         i2 = 3
      elseif (cval.lt.valcol(4)) then
         i1 = 3
         i2 = 4
      elseif (cval.lt.valcol(5)) then
         i1 = 4
         i2 = 5
      else
         i1 = 5
         i2 = 5
      endif

      if (i1.ne.i2) then
         factor = (cval - valcol(i1)) / (valcol(i2) - valcol(i1))
         do j=1,3
            col(j) = colmap(j,i1) + 
     &               factor*(colmap(j,i2) - colmap(j,i1))
         end do
      else
         do j=1,3
            col(j) = colmap(j,i1)
         end do
      endif

      return
      end

      subroutine plcini
      implicit double precision (a-h,o-z)
      common /vrcol/ jcol(3,16)

      data jcol(1,1),jcol(2,1),jcol(3,1) /255,0,0/
      data jcol(1,2),jcol(2,2),jcol(3,2) /255,159,9/
      data jcol(1,3),jcol(2,3),jcol(3,3) /0,255,0/
      data jcol(1,4),jcol(2,4),jcol(3,4) /78,255,187/
      data jcol(1,5),jcol(2,5),jcol(3,5) /0,255,255/
      data jcol(1,6),jcol(2,6),jcol(3,6) /255,191,0/
      data jcol(1,7),jcol(2,7),jcol(3,7) /132,193,214/
      data jcol(1,8),jcol(2,8),jcol(3,8) /115,115,115/
      data jcol(1,9),jcol(2,9),jcol(3,9) /255,0,255/
      data jcol(1,10),jcol(2,10),jcol(3,10) /16,176,16/
      data jcol(1,11),jcol(2,11),jcol(3,11) /239,202,140/
      data jcol(1,12),jcol(2,12),jcol(3,12) /255,122,0/
      data jcol(1,13),jcol(2,13),jcol(3,13) /230,214,92/
      data jcol(1,14),jcol(2,14),jcol(3,14) /184,56,6/
      data jcol(1,15),jcol(2,15),jcol(3,15) /255,255,255/

      idum = 1

      return
      end

      subroutine plphd(iun,ihand,tmp,t,tlpos)
      implicit double precision (a-h,o-z)
      common /vrcol/ jcol(3,16)
c
c --- variables for POV macros ---
c
      integer colors
      double precision colrgb
      dimension colrgb(3),tmp(3),t(3),tlpos(3)

      write(iun,'(''camera {'')')
      if (ihand.eq.1) then

         write(iun,
     &     '(''location <'',f12.5,'','',f12.5,'','',f12.5,''>'')')
     &     -tmp(2),tmp(1),-tmp(3)
         write(iun,
     &     '(''look_at <'',f12.5,'','',f12.5,'','',f12.5,''>'')')
     &      -t(2),t(1),-t(3)

      else

         write(iun,
     &     '(''location <'',f12.5,'','',f12.5,'','',f12.5,''>'')')
     &      (tmp(j),j=1,3)
         write(iun,'(''right < -1.33,  0,  0>'')')
         write(iun,
     &     '(''look_at <'',f12.5,'','',f12.5,'','',f12.5,''>'')')
     &      (t(j),j=1,3)

      endif

      write(iun,*) 'angle 45'
      write(iun,'(''}'')')

      if (ihand.eq.1) then
         write(iun,'(''sky_sphere {'')')
         write(iun,'(''  pigment {'')')
         write(iun,'(''    gradient y'')')
         write(iun,'(''    color_map {'')')
         write(iun,'(''      [0 color rgb <1, 0, 0>]'')')
         write(iun,'(''      [1 color rgb <0, 0, 1>]'')')
         write(iun,'(''    }'')')
         write(iun,'(''    scale 2'')')
         write(iun,'(''    translate -1'')')
         write(iun,'(''  }'')')
         write(iun,'(''}'')')
      else 
         write(iun,'(''background { color <0.6, 0.6, 0.6> }'')')
      endif

      write(iun,
     & '(''light_source { <'',f12.5,'','',f12.5,'','',f12.5,''>'')')
     &  (tlpos(j),j=1,3)
      write(iun,'('' color rgb<1, 1, 1> }'')')

      write(iun,*) '#declare BSAMBI = 0.2;'
      write(iun,*) '#declare BSDIFF = 0.8;'
      write(iun,*) '#declare BSSPEC = 0.8;'
c
c----- Add macros for colors ----- FPA 2/27/2000 
c
c general form is:
c  #declare color1 = 
c           texture {
c           pigment {color rgb <r1, r2, r3> }
c           finish {BSAMBI, BSDIFF, BSSPEC }
c           }
c
c
      do colors = 1,15
c              
           write(iun,'(a14,a1,a3)') '#declare color',
     1                      char(colors+64),' = '
           write(iun,*) 'texture { '
           do n = 1,3
                colrgb(n) = (real(jcol(n,colors))/255.0)
           end do
           write(iun,'(a15,3f8.4,a3)') ' pigment { rgb<',
     1                     (colrgb(n),n=1,3),
     2                     ' >}'
           write(iun,*)
     1        'finish {ambient BSAMBI diffuse BSDIFF specular BSSPEC}'
           write(iun,*)'}'

       end do
c
c --- Declare the molecule object ---- FPA 2/27/2000
c
       write(iun,*) '# declare molecule = union { '
c ---

       return
       end

      subroutine prtpsa
      implicit double precision (a-h,o-z)
      common /psa/ psa, tsa, exs, pol, pol2, epmin, epmax, ipsa,
     &             icpsa, idtpsa

      write(*,'(a,f8.5)') "Apolar Potential: < ",pol
      write(*,'(a,f8.5)') "Apolar Potential: > ",pol2
      write(*,'(a,f8.3,a)') "Polar Surface area: ",psa," Angstrom**2"
      write(*,'(a,f8.3,a)') "Total Surface area: ",tsa," Angstrom**2"
      write(*,'(a,f5.2,a)') "PSA percentage: ",100.0d0*psa/tsa," %"
      write(*,'(a,f8.3)') "Minimum Potential on surface: ",epmin
      write(*,'(a,f8.3)') "Maximum Potential on surface: ",epmax
      write(*,'(a,f10.6,a,f10.6)') "EXS: ",exs," ",exs/tsa

      if (idtpsa.eq.1) call tpsa

      return
      end

      subroutine wcubev(nx,ny,nz,fvals,origin,rad)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (numatm=2000)
      parameter (mxel=100)
      integer griddim(3), xydim, findex
      character*80 lstr
      common /grdhlp/ mx3d,mx3d2
      common /coord / xyz(3,numatm)
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      dimension fvals(*),origin(3),scale(3)
      dimension cubxyz(3), rad(3), c(3)
  
      open(unit=80,form='formatted',file='points',err=10)
      
      griddim(1) = nx
      griddim(2) = ny
      griddim(3) = nz

c rad(2) slaat op scale(1) en griddim(1)

      do i=1,3
          scale(i) = 1.0d0/(griddim(i)-1.0d0)
      end do
      scale(1) = scale(1) * rad(2) / rad(1)
      scale(3) = scale(3) * rad(3) / rad(1)

      xydim = mx3d2

      write(80,'(a)') '$CONTRL UNITS=BOHR $END'
      write(80,'(a)') '$DATA'
      write(80,'(a)') 'title'
      write(80,'(a)') 'c1'
      write(80,'(a)') ''

      do i=1,natoms
          write(80,'(a2,1x,i3,1x,3(f9.4,1x))') 
     &       elemnt(nat(i)),nat(i),(xyz(l,i),l=1,3)
      end do

      write(80,'(a)') '$END'
      write(80,'(a)') '$ELPOT IEPOT=1 WHERE=POINTS OUTPUT=BOTH $END'
      write(80,'(a)') '$POINTS'
      
      write(80,'(a,i8)') 'BOHR ',nx*ny*nz
      do k=0,griddim(3)-1
         cubxyz(3) = origin(3) + k*scale(3)

         do j=0,griddim(2)-1
            cubxyz(2) = origin(2) + j*scale(2)

            do i=0,griddim(1)-1
               cubxyz(1) = origin(1) + i*scale(1)

               findex = i + j*griddim(1) + k*xydim + 1
               cubval = fvals(findex)

               call rotbck(cubxyz(1),cubxyz(2),cubxyz(3),c)
               write(80,'(3(f9.4,1x))') (c(l),l=1,3)

            end do
         end do
      end do

      write(80,'(a)') '$END'
      close(80)

      return

10    print*,'error writing file'
      return
      end
