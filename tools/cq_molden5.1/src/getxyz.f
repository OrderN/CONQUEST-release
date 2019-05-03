      subroutine getxyd(igetxy,heat,iaddprv,
     &                  coo,ianz,iaton,iatclr,iconn,qat,
     &                  nat,norg,icent,nspg,ichx,nopr,ir,it,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,o-z), integer ( i-n)
      integer getlin
      character*137 str
      character*2 catom, catomt, tolowf,iel
      parameter (numatm=2000)
      parameter (mxcon=10)
      parameter (maxsym=108)
      common /athlp/  iatoms, mxnat
      common /zmfrst/ ihaszm, nz, mxzat
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /surf/   natorg,noscnd
      character*137 line
      common /curlin/ line
      logical dozme
      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb
      common /align/  vecs(3,3),nscnd,iscst,ialtyp,iocnt
      common /cllab/  iclon,iclpnt(4)
      common /cllhlp/ ico(3,8),icn(4,8),mcol(32)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      integer*2 ir,it
      logical addcel,gnreal,dall
      dimension iel(maxsym),icel(8),rr(3,3)
      dimension coo(3,*),ianz(*),iaton(*),iatclr(*),iconn(mxcon+1,*),
     &          qat(*),ir(3,3,192),it(3,192)
      data iel/'bq',
     &         'h ', 'he',
     &         'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne',
     &         'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar',
     &         'k ', 'ca',
     &                     'sc', 'ti', 'v ', 'cr', 'mn',
     &                     'fe', 'co', 'ni', 'cu', 'zn',
     &                     'ga', 'ge', 'as', 'se', 'br', 'kr',
     & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     & 'in','sn','sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     & 'pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     & 'ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po',
     & 'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm',
     & 'bk','cf','x ','oo','oa','ob','oc','ab','bc','ac','zz'/

      igetxy = 1
      addcel = .false.
      dall = .false.
      ncell = 0
      toang = 0.52917706d0
      todeg = 45.0d0 / datan(1.0d0)
      igusxyz = 0

c get number of atoms

      if (ircus.eq.1) then
         call searchq(line,'***** BEGIN IRC',
     &      'TIME                Q          P',
     &      'KINETIC       POTENTIAL       TOTAL',
     &      '$data',istat)
         if (istat.ne.0) then
            igust = 0
            if (index(line,'KINETIC').ne.0.or.
     &          index(line,'TIME').ne.0) then
               igust = 1
               call readel(line,1)
            endif
            if (icdex(line,'$data').ne.0) then
               igusxyz = 1
            endif
         else
            goto 100
         endif
      else
         if (getlin(0).eq.1) then
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.2) then
               if (iaddprv.eq.1) then
                   iscst = iatoms
                   nscnd = itype
                   natoms = nscnd
                   if (iatoms+nscnd.gt.mxnat) then
                       dall = .true.
                       goto 100
                   endif
               else
                   iscst = 0
                   iatoms = itype
                   natoms = iatoms
                   if (iatoms.gt.mxnat) then
                       dall = .true.
                       goto 100
                   endif
               endif
            else
               goto 100
            endif
         else
            goto 100
         endif
      endif

c get Heat of formation

      heat = 0.0d0
      heato = 0.0d0
      iefnd = 0
      if (igusxyz.eq.0) then
         if (getlin(1).eq.1) then
             do while(.true.)
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ircus.eq.1) then
                   if (ktype.eq.1) then
                       if (nstr.eq.2) then
                          if (str.eq.'E=') iefnd = 1
                       endif
                   endif
                endif
                if (ktype.eq.3) then
                   heato = heat
                   heat = rtype
                endif
                if (ircus.eq.1) then
                   if (ktype.eq.0.or.(ktype.eq.3.and.iefnd.eq.1)) then
                      if (igust.eq.1) heat = heato
                      goto 110
                   endif
                else
                   if (ktype.eq.0.or.ktype.eq.3) goto 110
                endif
             end do
         else
            goto 100
         endif
      else
         call readel(line,1)
      endif
      
110   continue

      if (ircus.eq.1) then
         if (igust.eq.1) then
            call search(line,'CARTESIAN COORDINATES',istat)
         else
            call readel(line,1)
         endif
         iatoms = 0
         do while (getlin(0).eq.1)

            if (line(1:4).eq.' ---') goto 120
            if (igusxyz.ne.0) then
               if (icdex(line,'$end').ne.0) goto 120
            endif
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.1) then
               if (nstr.ge.4) then
                  if (str(1:4).eq.'****') goto 120
                  if (str(1:4).eq.'MASS') goto 120
               endif
            endif
            iatoms = iatoms + 1

            if (igust.eq.0.or.igusxyz.eq.1) 
     &           ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.3) then
               ianz(iatoms) = int(rtype)
            else
               goto 100
            endif

            if (.not.gnreal(coo(1,iatoms),3,.false.)) goto 100
c            print*,(coo(i,iatoms),i=1,3)
         end do

120      continue
         backspace iun2
         call doconn
         ihaszm = 0
         return
      endif

      do iat=iscst+1,iscst+natoms
         if (getlin(0).eq.1) then
           ktype = nxtwrd(str,nstr,itype,rtype)
           if (ktype.eq.1) then
                 iatmp = 0
                 if (nstr.eq.1) then
                    catomt(1:1) = str(1:1)
                    catomt(2:2) = ' '
                 else
                    catomt = str(1:2)
                 endif
                 catom = tolowf(catomt)
                 do j=1,maxsym
                    if (catom .eq. iel(j)) iatmp = j - 1
                 end do
                 if (catom.eq.'xx'.or.catom.eq.'bq') iatmp = 99
                 if (iatmp.le.0.or.iatmp.gt.maxsym-1) goto 100
                 if (iatmp.eq.106.and.ncell.eq.0) iatmp = 89
                 if (iatmp.gt.99) then
                    icel(iatmp-99) = iat
                    addcel = .true.
                    ncell = ncell + 1
                 else
                    if (addcel) then
                        call inferr('Cell points before atoms',1)
                        goto 100
                    endif
                 endif
                 ianz(iat) = iatmp
           else
              goto 100
           endif
           do i=1,4
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.3) then
                 if (i.eq.4) then
                    qat(iat) = rtype
                    ihasq = 1
                 else
                    coo(i,iat) = rtype / toang
                 endif
              elseif (ktype.eq.2.and.i.ne.4) then
                 coo(i,iat) = dble(itype) / toang
              else
                 if (i.ne.4) goto 100
              endif
           end do
        else
          goto 100
        endif
      end do

      iatoms = natoms + iscst

c look for cell data

      if (getlin(0).ne.1) goto 130
      if (getlin(0).ne.1) then
         backspace iun2
         goto 130
      endif

      if (icdex(line,'Unit cell').eq.0) then

         backspace iun2
         backspace iun2

      else

         do i=1,3
           if (getlin(0).eq.1) then
              if (.not.gnreal(rr(1,i),3,.false.)) goto 100
           endif
         end do
         nspg = 1
         a = vlen(rr(1,1))
         b = vlen(rr(1,2))
         c = vlen(rr(1,3))
         call impsc(rr(1,2),rr(1,3),csa)
         call impsc(rr(1,1),rr(1,3),csb)
         call impsc(rr(1,1),rr(1,2),csc)
         alpha = dacos(csa)*todeg
         beta = dacos(csb)*todeg
         gamma = dacos(csc)*todeg
         call prcell(nspg,a,b,c,alpha,beta,gamma)
         call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)
         call cprot(nspg,nopr,icent,ir,it,.false.)
         call prop(nopr,ir,it)
         call cpmol(nat,norg,a,b,c,alpha,beta,gamma,
     &              coo,ianz,iatclr,iconn)

         iftyp = 6
         ichx = 1
         call fdat(1,0,0,0,0,0)

      endif

130   continue

      if (iaddprv.eq.1) then
         do i=iscst+1,iscst+nscnd
             do j=1,3
                iaton(i) = 1
             end do
         end do
      endif

      if (addcel) iatoms = iatoms - 8
      if (iaddprv.eq.1) then
         call convar(coo,iconn,ianz,nscnd,iscst)
      else
         call cooxyz(ianz,iatoms)
         call doconn
      endif
      if (addcel) then
         natorg = iatoms
         iclon = 1
         iclpnt(1) = icel(1)
         iclpnt(2) = icel(2)
         iclpnt(3) = icel(3)
         iclpnt(4) = icel(4)
         do i=1,8
            l = icel(i)
            iconn(1,l) = icn(1,i)
            do j=2,4
               iconn(j,l) = iatoms + icn(j,i)
            end do
            ianz(l) = 100
            iatclr(l) = 11
         end do
         iatoms = iatoms + 8
      endif

      if (iaddprv.eq.0) ihaszm = 0

      return

100   if (dall) then
         rewind iun2
         iatoms = iscst
         igetxy = -1
         return
      endif
      igetxy = 0
      return
      end

      subroutine aln2md(iopt,istat,
     &                  coo,ianz,iaton,iatclr,iresid)
      implicit double precision (a-h,o-z)
      logical fndsnd,ocen
      parameter (mxpsav=100)
      common /athlp/ iatoms, mxnat
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt
      common /cllmat/rrx(3),rry(3),rrz(3),tr(3),tz(3),tzorg(3),itz
      common /cllsav/rdum(12,mxpsav),npsav
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      dimension vec1(3),vec2(3),vec3(3)
      dimension coo(3,*),ianz(*),iaton(*),iatclr(*),iresid(*)

      istat = 1
      fndsnd = .false.
      ocen = .true.

      ialtyp = iopt

      if (ialtyp.eq.0) then
         call getxyz(igetxy,heat,1)
         if (igetxy.eq.1) then
             fndsnd = .true.
         else
             fndsnd = .false.
         endif
         if (.not.fndsnd) then
             call rdsrf(iun2,istats,iesp,1,0)
             if (istats.eq.1) fndsnd = .true.
         endif
         if (.not.fndsnd) then
            print*,'This is not an xyz or molden surface file !'
            istat = 0
            return
         endif
      else if (ialtyp.eq.1) then
         ioatms = iatoms
         do i=1,iatoms
            if (iresid(i).le.0.and.iresid(i).ge.-3) then
               ioatms = i - 1
               goto 3000
            endif
         end do
3000     namols = namols + 1
         iatoms = ioatms
         iscst = ioatms
         
         call pdbstd(istat,.true.,1)
         if (istat.eq.-1) then
            call allcoo(50000,0)
            call pdbstd(istat,.true.,1)
         endif
         if (istat.eq.1) then
            fndsnd = .true.
            call xyzcoo(0,1,1)
c
c commented out calfa because it added hydrogens to structure1 in 
c the space of structure 2: it SHOULD NOT do that
c            call calfa(istat,1,1,ioatms+1,1,1)
c
            call dohcon(2)
            ocen = .false.
            nscnd = iatoms - ioatms
         else
            istat = 0
            print*,'This is not a valid PDB file !'
            return
         endif
      else if (ialtyp.eq.2) then
         namols = namols + 1
         fndsnd = .true.
         ialtyp = 1
         ocen = .false.
      else if (ialtyp.eq.3) then
         namols = namols + 1
         fndsnd = .true.
         ialtyp = 0
      endif

      if (fndsnd) then

         call cntvec(vec1,coo,ianz,iscst)
         call cntvec(vec2,coo(1,iscst+1),ianz(iscst+1),nscnd)

         do i=1,3
             vec3(i) = vec1(i) - vec2(i) 
         end do

         if (ialtyp.ne.1) then

            do i=1,nscnd
               call trcoo(vec3,coo(1,iscst+i))
            end do

            do i=1,3
               tr(i) = vec1(i)
            end do

         else

c set center of mol2  (tr)
c set pos mol2 to center of mol1 (tz)

            do i=1,3
               tr(i) = vec2(i)
               if (iocnt.eq.1) then
                  tz(i) = vec3(i)
               else
                  tz(i) = 0.0d0
               endif
            end do
c
c           copy coordinates of second molecule to the top of the
c           coordinate array
c
            nstrt = mxnat-nscnd
            do i=1,nscnd
               do j=1,3
                  coo(j,nstrt+i) = coo(j,iscst+i)
               end do
               ianz(nstrt+i) = ianz(iscst+i)
            end do
         endif

         if (ocen) call docent


         do i=1,iscst
             if (ianz(i).ne.100) then
                iatclr(i) = 1
             endif
         end do
         irs = iresid(iscst+1)
         if (irs.le.0.and.irs.ge.-3) irs = -4
         do i=iscst+1,iscst+nscnd
             iatclr(i) = 5
             iaton(i) = 1
c             iresid(i) = -4
             if (irs.eq.-4) iresid(i) = irs
         end do

         if (ialtyp.eq.1) then
c vec1 is just a dummy in this case, doesnt do anything
c            call inirot
            call clini
            call alnrot(vec1,1)
            call ligzmt
            call pmfass(0,0)
            npsav = 0
         endif

         call mtinv3
         call rarbxi
         call rarbyi
         call rarbzi

      endif

      return
      end

      subroutine alnrod(vec,irot,coo)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      common /cllmat/rx(3),ry(3),rz(3),t(3),tz(3),tzorg(3),itz
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt
      dimension vec(3),tzt(3)
      dimension coo(3,*)

      nstrt = mxnat-nscnd

      if (irot.eq.1) then

c do rotation

         if (ialtyp.eq.1) then

            do i=1,nscnd
               x = coo(1,nstrt+i)
               y = coo(2,nstrt+i)
               z = coo(3,nstrt+i)
               coo(1,iscst+i) = 
     &           (x-t(1))*rx(1)+(y-t(2))*rx(2)+(z-t(3))*rx(3)+t(1)
               coo(2,iscst+i) = 
     &           (x-t(1))*ry(1)+(y-t(2))*ry(2)+(z-t(3))*ry(3)+t(2)
               coo(3,iscst+i) = 
     &           (x-t(1))*rz(1)+(y-t(2))*rz(2)+(z-t(3))*rz(3)+t(3)
               do j=1,3
                  coo(j,iscst+i) = coo(j,iscst+i) + tz(j)
               end do
            end do

         else

            do i=iscst+1,iscst+nscnd
               x = coo(1,i)
               y = coo(2,i)
               z = coo(3,i)
               coo(1,i) = 
     &            (x-t(1))*rx(1)+(y-t(2))*rx(2)+(z-t(3))*rx(3)+t(1)
               coo(2,i) = 
     &            (x-t(1))*ry(1)+(y-t(2))*ry(2)+(z-t(3))*ry(3)+t(2)
               coo(3,i) = 
     &            (x-t(1))*rz(1)+(y-t(2))*rz(2)+(z-t(3))*rz(3)+t(3)
            end do

         endif

      else

c do translation

         if (ialtyp.eq.1) then
            do j=1,3
               tz(j) = tz(j) - vec(j)
            end do
            do i=1,nscnd
               x = coo(1,nstrt+i)
               y = coo(2,nstrt+i)
               z = coo(3,nstrt+i)
               coo(1,iscst+i) = 
     &           (x-t(1))*rx(1)+(y-t(2))*rx(2)+(z-t(3))*rx(3)+t(1)
               coo(2,iscst+i) = 
     &           (x-t(1))*ry(1)+(y-t(2))*ry(2)+(z-t(3))*ry(3)+t(2)
               coo(3,iscst+i) = 
     &           (x-t(1))*rz(1)+(y-t(2))*rz(2)+(z-t(3))*rz(3)+t(3)
               do j=1,3
                  coo(j,iscst+i) = coo(j,iscst+i) + tz(j)
               end do
            end do
            do j=1,3
               tzt(j) = tzorg(j) - tz(j)
            end do
            vl = vlen(tzt)
            if (vl.gt.1.0) call updres
         else
            do i=iscst+1,iscst+nscnd
               do j=1,3
                  coo(j,i) = coo(j,i) - vec(j)
               end do
            end do
            do j=1,3
               t(j) = t(j) + vec(j)
            end do
         endif

      endif

      return
      end

      subroutine alnsorg(iatom,coo)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      common /cllmat/rx(3),ry(3),rz(3),t(3),tz(3),tzorg(3),itz
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt
      dimension coo(3,*)

      nstrt = mxnat-nscnd

      if (ialtyp.eq.1) then

         if (iatom.le.iscst) then
             call inferr('This is not an atom of docked molecule !',0)
             return
         endif
         im = iatom - iscst

         x = coo(1,nstrt+im)
         y = coo(2,nstrt+im)
         z = coo(3,nstrt+im)

         xn = 
     &      (x-t(1))*rx(1)+(y-t(2))*rx(2)+(z-t(3))*rx(3)-(x-t(1))
         yn = 
     &      (x-t(1))*ry(1)+(y-t(2))*ry(2)+(z-t(3))*ry(3)-(y-t(2))
         zn = 
     &      (x-t(1))*rz(1)+(y-t(2))*rz(2)+(z-t(3))*rz(3)-(z-t(3))

         t(1) = x
         t(2) = y
         t(3) = z

         tz(1) = tz(1) + xn
         tz(2) = tz(2) + yn
         tz(3) = tz(3) + zn

      else

         do i=1,3
            t(i) = coo(i,iatom)
         end do

      endif


      return
      end

      subroutine alnsav()
      implicit double precision (a-h,o-z)
      parameter (mxpsav=100)
      common /cllmat/rx(3),ry(3),rz(3),t(3),tz(3),tzorg(3),itz
      common /cllsav/rxs(3,mxpsav),rys(3,mxpsav),rzs(3,mxpsav),
     &               tzs(3,mxpsav),npsav

      if (npsav.lt.mxpsav) then
         npsav = npsav + 1

         do j=1,3
            rxs(j,npsav) = rx(j)
            rys(j,npsav) = ry(j)
            rzs(j,npsav) = rz(j)
            tzs(j,npsav) = tz(j)
         end do

      endif

      return
      end

      subroutine alnwrd(coo,ianz)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      common /athlp/ iatoms, mxnat
      parameter (mxpsav=100)
      common /cllmat/rx(3),ry(3),rz(3),t(3),tz(3),tzorg(3),itz
      common /cllsav/rxs(3,mxpsav),rys(3,mxpsav),rzs(3,mxpsav),
     &               tzs(3,mxpsav),npsav
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt
      character*2 ggstr,tocapf, elemnt
      common /elem/elemnt(mxel)
      logical opfil
      dimension xyzt(3)
      dimension coo(3,*),ianz(*)

      toang = 0.52917706d0
      iun = 46
      ires = 1
      nstrt = mxnat-nscnd

      do l=1,npsav
         if (opfil(iun,'pose.'//ggstr(l-1),7,1,0,0)) then
            write(iun,'(a)') 'HEADER'
            do i=1,nscnd
               x = coo(1,nstrt+i)
               y = coo(2,nstrt+i)
               z = coo(3,nstrt+i)
               xyzt(1) = (x-t(1))*rxs(1,l)+(y-t(2))*rxs(2,l)+
     &                   (z-t(3))*rxs(3,l)+t(1)
               xyzt(2) = (x-t(1))*rys(1,l)+(y-t(2))*rys(2,l)+
     &                   (z-t(3))*rys(3,l)+t(2)
               xyzt(3) = (x-t(1))*rzs(1,l)+(y-t(2))*rzs(2,l)+
     &                   (z-t(3))*rzs(3,l)+t(3)
               do j=1,3
                  xyzt(j) = xyzt(j) + tzs(j,l)
               end do
               ian = ianz(nstrt+i)
               write(iun,'(a6,i5,1x,a4,1x,i3,2x,i4,4x,3f8.3)')
     &            'HETATM',i,tocapf(elemnt(ian))//'  ',
     &             ires,ires,(xyzt(j)*toang,j=1,3)
            end do
            close(iun)
         endif
      end do

      return
      end

      subroutine alnsed(isel,coo)
      implicit double precision (a-h,o-z)
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt
      integer dolabs,fancy,persp,shade,atcol,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      dimension isel(6),isel1(3),isel2(3),vec(3),coo(3,*)

      do i=1,3
         j = isel((i-1)*2+1)
         k = isel((i-1)*2+2)
         if (j.le.iscst.and.k.gt.iscst) then
            isel1(i) = j
            isel2(i) = k - iscst
         elseif (k.le.iscst.and.j.gt.iscst) then
            isel1(i) = k
            isel2(i) = j - iscst
            if (i.eq.1) then
                idum = isel(1)
                isel(1) = isel(2)
                isel(2) = idum
            endif
         else
            call inferr('atom pair of same molecule ',0)
            return
         endif
      end do

      do i=1,3
          vec(i) = coo(i,isel(1)) - coo(i,isel(2))
      end do

      do i=1,nscnd
         call trcoo(vec,coo(1,iscst+i))
      end do
   
      call alntwo(coo,isel1,coo(1,iscst+1),nscnd,isel2)

      call setorg(isel(2))

      irtcel = 0

      return
      end

      subroutine getchd(igetca,qat)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (numatm=2000)
      integer getlin
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      character*137 str
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      character*137 line
      common /curlin/ line
      logical gnreal
      dimension r(3),qat(*)

      igetca = 1
      if (getlin(0).eq.1) then
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            iatoms = itype
            if (iatoms.ne.natoms) goto 100
         else
            goto 100
         endif
      else
         goto 100
      endif
      call readel(line,1)

      do i=1,iatoms
         if (getlin(0).eq.1) then
           ktype = nxtwrd(str,nstr,itype,rtype)
           if (ktype.ne.1) goto 100
           if (.not.gnreal(r,3,.false.)) goto 100
           ktype = nxtwrd(str,nstr,itype,rtype)
           if (ktype.eq.3) then
               qat(i) = rtype
           else
               goto 100
           endif
         endif
      end do

      ihasq = 1

      return
100   igetca = 0
      return
      end

      double precision function deter(r)
      implicit double precision (a-h,o-z)
      dimension r(3,3)
c  calculate the determinant of a 3x3 matrix

      deter = r(1,1)*r(2,2)*r(3,3) - r(1,3)*r(2,2)*r(3,1) +
     &        r(2,1)*r(3,2)*r(1,3) - r(2,3)*r(3,2)*r(1,1) +
     &        r(3,1)*r(1,2)*r(2,3) - r(3,3)*r(1,2)*r(2,1)
      return
      end

      subroutine invmat3(r,ri)
      implicit double precision (a-h,o-z)
      dimension r(3,3),ri(3,3)

c  invert matrix M

      det = deter(r)
      ri(1,1) = (r(2,2)*r(3,3)-r(2,3)*r(3,2)) / det
      ri(1,2) = (r(1,3)*r(3,2)-r(1,2)*r(3,3)) / det
      ri(1,3) = (r(1,2)*r(2,3)-r(1,3)*r(2,2)) / det
      ri(2,1) = (r(2,3)*r(3,1)-r(2,1)*r(3,3)) / det
      ri(2,2) = (r(1,1)*r(3,3)-r(1,3)*r(3,1)) / det
      ri(2,3) = (r(1,3)*r(2,1)-r(1,1)*r(2,3)) / det
      ri(3,1) = (r(2,1)*r(3,2)-r(2,2)*r(3,1)) / det
      ri(3,2) = (r(1,2)*r(3,1)-r(1,1)*r(3,2)) / det
      ri(3,3) = (r(1,1)*r(2,2)-r(1,2)*r(2,1)) / det

      return
      end

      subroutine getmod(igetmo,coo,qat,ianz,iconn)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (numatm=2000)
      parameter (mxcon=10)
      parameter (maxsym=108)
      character*80 title
      character*2 catom, catomt, tolowf,iel
      character*5 vers
      common /athlp/ iatoms, mxnat
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      character*137 line
      common /curlin/ line
      character*137 cas,nsc
      common /regis/  cas,nsc
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      logical ocorin
      dimension iel(maxsym),ichg(8),iqtmp(8)
      dimension coo(3,*),ianz(*),qat(*),iconn(mxcon+1,*)
      data iel/'lp',
     &         'h ', 'he',
     &         'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne',
     &         'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar',
     &         'k ', 'ca',
     &                     'sc', 'ti', 'v ', 'cr', 'mn',
     &                     'fe', 'co', 'ni', 'cu', 'zn',
     &                     'ga', 'ge', 'as', 'se', 'br', 'kr',
     & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     & 'in','sn','sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     & 'pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     & 'ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po',
     & 'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm',
     & 'bk','cf','x ','oo','oa','ob','oc','ab','bc','ac','zz'/

      igetmo = 1
      ocorin = .false.

      iatoms = 0
      nbnds = 0

      toang = 0.52917706d0

      call rewmf

c get title line

      read (iun2,'(a)',end=100,err=100) title
      if (title(1:4).eq.'$$$$') then
         call readel(line,1)
         read (iun2,'(a)',end=100,err=100) title
      endif

      read (iun2,'(a)',end=100,err=100) line

      read (iun2,'(a)',end=100,err=100) line
      if (line(1:6).eq.'Corina') ocorin = .true.

      read (iun2,'(a)',end=100,err=100) line
      if (linlen(line).ge.39) then
         read (line,'(i3,i3,6x,i1,21x,a5)',end=100,err=100) 
     &   iatoms,nbnds,ister,vers
         if (vers.ne.'V2000'.and.vers.ne.'V3000') goto 100
      else
         if (linlen(line).ge.6) then
            read (line,'(i3,i3)',end=100,err=100) 
     &            iatoms,nbnds
         else
            goto 100
         endif
      endif

      do i=1,mxnat
         qat(i) = 0.0d0
      end do
      ichgrec = 0

      do i=1,iatoms
         read (iun2,'(3f10.4,1x,a2,6x,i3)',end=100,err=100) 
     &       (coo(j,i),j=1,3),catomt,itchg
         do j=1,3
            coo(j,i) = coo(j,i) / toang
         end do
         iatmp = 0
         catom = tolowf(catomt)
         do j=1,maxsym
            if (catom .eq. iel(j)) iatmp = j - 1
         end do
         if (catom.eq.'lp') iatmp = 99
         if (iatmp.le.0.or.iatmp.gt.maxsym-1) goto 100
         ianz(i) = iatmp
         qat(i) = dble(itchg)
         iconn(1,i) = 0
      end do

      do i=1,nbnds
         read (iun2,'(3i3)',end=100,err=100) 
     &      iat1,iat2,ibndtyp
         if (iconn(1,iat1).lt.mxcon) then
            iconn(1,iat1) = iconn(1,iat1) + 1
            iconn(1+iconn(1,iat1),iat1) = iat2
         endif
         if (iconn(1,iat2).lt.mxcon) then
            iconn(1,iat2) = iconn(1,iat2) + 1
            iconn(1+iconn(1,iat2),iat2) = iat1
         endif
      end do


      do while (.true.) 
         read (iun2,'(a)',end=200,err=100) line
         if (linlen(line).ge.1) then
            if (line(1:6).eq.'M  CHG') then
               ihasq = 1
               if (ichgrec.eq.0) then
                  do i=1,iatoms
                     qat(i) = 0.0d0
                  end do
                  ichgrec = 1
               endif
               nchg = 0
               read(line,'(6x,i3)') nchg
               read(line,'(9x,8(1x,i3,1x,i3))') 
     &               (ichg(i),iqtmp(i),i=1,nchg)
               do i=1,nchg
                 qat(ichg(i)) = dble(iqtmp(i))
               end do
            else
               if (index(line,'<NSC>').ne.0) then
                  read (iun2,'(a)',end=200,err=100) nsc
               endif
               if (index(line,'<CAS_RN>').ne.0) then
                  read (iun2,'(a)',end=200,err=100) cas
               endif
            endif
         endif
      end do

      call cooxyz(ianz,iatoms)
200   continue
      ihaszm = 0

      return
100   igetmo = 0
      return
      end

      subroutine outmod(iun,coo,ianz,iconn)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (mxcon=10)
      parameter (maxsym=108)
      parameter (mxel=100)
      character*2 elemnt,tocapf
      character*3 atom,atomt
      common /athlp/  iatoms, mxnat
      common /elem/   elemnt(mxel)
      character*137 line
      common /curlin/ line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension coo(3,*),ianz(*),iconn(mxcon+1,*)

      toang = 0.52917706d0
      natoms = iatoms
      nbnds = 0
      idochg = 0

      do i=1,natoms
         do j=1,iconn(1,i)
            if (iconn(1+j,i).gt.0) then
               if (iconn(1+j,i).gt.i) nbnds = nbnds + 1
            endif
         end do
      end do

      write(iun,'(a)') ' '
      write(iun,'(a)') '  -MOLDEN-'
      write(iun,'(a)') 'Molden generated mol file'
      write(iun,'(i3,i3,a)') natoms,nbnds,
     &          '  0  0  0  0  0  0  0  0999 V2000'

      do i=1,natoms
         atomt = '   '
         atomt(1:2) = tocapf(elemnt(ianz(i)))
         call leftj(atomt,atom)
         write(iun,'(3f10.4,a1,a3,a)') (coo(j,i)*toang,j=1,3),' ',
     &      atom, ' 0  0  0  0  0  0  0  0  0  0  0  0'
      end do

      do i=1,natoms
         do j=1,iconn(1,i)
            k = iconn(1+j,i)
            if (k.gt.0) then
               if (k.gt.i) then
                  ibt = ibtyp(i,k,idochg,0,ianz)
                  write(iun,'(3i3,a)') i,k,ibt,
     &               '  0  0  0  0'
               endif
            endif
         end do
      end do

      write(iun,'(a)') 'M  END'

      return
      end
