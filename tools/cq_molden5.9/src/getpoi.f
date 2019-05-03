      subroutine getpod(icom,ifd,iff,doscl,ioatms,ioadd,
     &                  coo,ianz,iaton,iconn,
     &                  scal,scali,smag,natc,ichx,icrtp)

      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (mxel=100)
      parameter (mxonh=100)
      parameter (numatm=2000)
      parameter (small=1.0d-4)
      common /pnthlp/ ipoints,ipnt
      common /athlp/ iatoms, mxnat
      common /gamori/imap(numat1),cstand(3,numat1),istand(numat1),
     &               msucc
      common /oniomh/ ioni,nion,ionih(mxonh),natonh(mxonh),
     &                iomap(numatm),icntat(mxonh),fct(mxonh),
     &                xyzi(3,mxonh)
      common /gauori/ nzm,nso,nio,nzo,ioropt,ifor,
     &                ixyz98,iopr,isymm,irc,imp2,icntp,itd
      common /align/  vecs(3,3),nscnd,iscst,ialtyp,iocnt
      common /zmpart/ ipart,imn,imx,idcur
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /geocnv/ fdum3(10),idum(7),ieav,ifrav,mxpnt
      common /mopver/ mopopt
      common /animo/  movie

      character lstr*137
      character*6 astr
      character*8 refcod
      integer first,prev,next,last,getmf,doscl
      logical dummy,dozme
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /cllab/  iclon,iclpnt(4)
      common /ecce/iecce
      character*23 stmp
      dimension dum(3),coo(3,*),ianz(*),iaton(*),iconn(mxcon+1,*)

c iftyp      file type
c
c 1          mopac
c 2          gamess
c 3          gamus
c 4          gauss
c 5          adfin
c 6          chemx
c 7          cpmd
c 8          qchem
c 9          orca
c 10         xyz
c 11         tinker + small
c 12         tinker protien
c 13         mol file
c 15         nwchem

cd     write(iun3,'(a)')'enter subroutine getpoi'

      first = -1
      next  = -2
      last  = -3
      prev  = -4
      tol = small
      idoprt = 0

      icomm = icom

c     parse variables to C used in getpoi

      call parpoi(nzm,nso,nio,nzo,ioropt,ifor,ixyz98,iopr,isymm,irc,
     &            imp2,icntp,msucc,ioni,mopopt,isbin,irtype,
     &            ipdbgro,ifav,ioxyz,iconv,ircus,nscnd,iscst,ialtyp)

c for gaussian
 
      stmp = 'Input orientation:'
      ns = 18
      if (ixyz98.gt.0.and.(ioropt.eq.2)) then
          stmp = 'Standard orientation:'
          ns = 22
      endif

      if (dozme) then
c since dumzm resets the ipart variable, we need the help variable dopart
c this is because otherwise conpdb downstream goes wrong
          if (ipart.eq.1) idoprt = ipart
          ifav = 0
          call dumzm(coo,ianz,iatoms)
          goto 100
      endif

c     mopac binary .gpt or GRAPH file (msf excluded by
c                                      iftyp.ne.6: .not.chemx)

      if (isbin.eq.1.and.iftyp.ne.6) iatoms = 0

c     convert mopac coordinates to a.u.
c     vdwr are in angstrom, vrad in a.u.

      if (iftyp.eq.1.and.iconv.eq.1) call xyzcoo(1,1,0)

c     CPMD

      if (iftyp.eq.7.and.ipoints.eq.0.and.
     &    (ioxyz.eq.1.or.irtype.eq.4)) then
          call xyzcoo(1,0,0)
          call haszm(.false.)
      endif

c     GAMESS xyz optimise or FREQ job

      if (iftyp.eq.2.and.ipoints.eq.0.and.
     &    (ioxyz.eq.1.or.irtype.eq.4)) then
          call xyzcoo(1,0,0)
          call haszm(.false.)
      endif

c     GAUSSIAN single point or FREQ job, the first point

      if (iftyp.eq.4.and.
     &    (ipoints.le.1.or.irtype.eq.4).and.icomm.eq.-1) then
          call rewfil
          call gaupoi(istat)
          if (istat.eq.-1) then
             ifav = 0
          else
             ifav = 1
          endif
          if ((istat.eq.0.or.irtype.eq.4).and.ioni.ne.1) 
     &          call xyzcoo(1,0,0)
          if (istat.eq.0) call haszm(.false.)
      endif
      
      if (icomm.eq.prev) then
          if ( ipnt.le.1 ) return
          ipnt = ipnt -1
          icomm = ipnt
      endif

      ifav = 0

c     NOT {GAMESS-UK,GAMESS-US,GAUSSIAN or CPMD}
c     AND iatoms=0: EXIT

      if (iatoms.eq.0.and..not.
     &    ((iftyp.ge.2.and.iftyp.le.4).or.iftyp.eq.7)) then
          call inferr('No coordinates found !',0)
cd         write(iun3,'(a)')'leave subroutine getpoi'
          return
      endif

c Mopac2007 AUX file

      if (iftyp.eq.1.and.mopopt.eq.4.and.ipoints.ge.1) then

          if (icomm .eq. first) then
             call rewfil
             call search(lstr,'HEAT_OF_FORM_UPDATED',istat)
             call redel(lstr,2)
             if (istat.eq.1) call gtapnt(iatoms,coo,istat) 
             if (istat.eq.1) ipnt = 1
          elseif (icomm .eq. next) then
             call search(lstr,'HEAT_OF_FORM_UPDATED',istat)
             call redel(lstr,2)
             if (istat.eq.1) call gtapnt(iatoms,coo,istat) 
             if (istat.eq.1) then
                ipnt = ipnt +1
             else
                call rewfil
                call search(lstr,'HEAT_OF_FORM_UPDATED',istat)
                call redel(lstr,2)
                if (istat.eq.1) call gtapnt(iatoms,coo,istat) 
                if (istat.eq.1) ipnt = 1
             endif
          elseif (icomm.gt.0) then
             call rewfil
             do i=1,icomm
                call search(lstr,'HEAT_OF_FORM_UPDATED',istat)
                call redel(lstr,2)
                if (istat.eq.1) call gtapnt(iatoms,coo,istat) 
                if (istat.eq.1) then
                    dummy = .true.
                else
                    dummy = .false.
                endif
             end do
             ipnt = icomm
          endif
c         convert coo from angs to au
          call xyzcoo(0,1,0)
      endif

      if (iftyp.eq.1.and.iconv.eq.0.and.ipoints.ge.1) then
          if (icomm .eq. first) then
             call rewfil
             call getmop(iatoms,heat,0,1,istat) 
             if (istat.eq.1) ipnt = 1
          elseif (icomm .eq. next) then
             call getmop(iatoms,heat,0,0,istat) 
             if (istat.eq.1) then
                ipnt = ipnt +1
             else
                call rewfil
                call getmop(iatoms,heat,0,0,istat) 
                if (istat.eq.1) ipnt = 1
             endif
          elseif (icomm.gt.0) then
             call rewfil
             idomopf = 1
             do i=1,icomm
                call getmop(iatoms,heat,0,idomopf,istat) 
                if (istat.eq.1) then
                    dummy = .true.
                else
                    dummy = .false.
                endif
                idomopf = 0
             end do
             ipnt = icomm
          endif
      endif

      if ((iftyp.ge.10.or.(iftyp.eq.5.and.ihasg.ge.1)).and.
     &    ipoints.gt.1) then

          if (icomm .eq. first) then

             if (iftyp.eq.5) then
                call rewmf
                call srchmf(lstr,'[GEOMETRIES]',istat)
             else
                call rewfil
             endif

             if (iftyp.eq.5.and.ihasg.eq.2) then
                call getzm(iatoms,0,0,istat)
                if (istat.eq.1) ipnt = 1
             else
                if (iftyp.eq.11.or.iftyp.eq.12) then
                   call gettnk(igttnk,0,ipdbon,iff,iheat,heat)
                   if (igttnk.eq.1) ipnt = 1
                elseif (iftyp.eq.13) then
                   call getmol(igetmo,0)
                   if (igetmo.eq.1) ipnt = 1
                else
                   call getxyz(igetxy,heat,0)
                   if (igetxy.eq.1) ipnt = 1
                endif
             endif

             if (iftyp.eq.5) then
                iolin = getmf()
                call rdfc(ipnt,istat)
                if (istat.eq.1) then
                   ifav = 1
                else
                   ifav = 0
                endif
                call putmf(iolin)
             endif

          elseif (icomm .eq. next.and..not.
     &            (iftyp.eq.5.and.ifrav.eq.1)) then
             if (iftyp.eq.5.and.ihasg.eq.2) then
                call getzm(iatoms,0,0,istat)
                if (istat.eq.1) then
                   dummy = .true.
                else
                   dummy = .false.
                endif
             else
                if (iftyp.eq.11.or.iftyp.eq.12) then
                   call gettnk(igttnk,0,ipdbon,iff,iheat,heat)
                   if (igttnk.eq.1) then
                      dummy = .true.
                   else
                      dummy = .false.
                   endif
                elseif (iftyp.eq.13) then
                   call getmol(igetmo,0)
                   if (igetmo.eq.1) then
                      dummy = .true.
                   else
                      dummy = .false.
                   endif
                else
                   call getxyz(igetxy,heat,0)
                   if (igetxy.eq.1) then
                       dummy = .true.
                   else
                       dummy = .false.
                   endif
                endif
             endif
             if (dummy) then
                ipnt = ipnt +1
             else
                if (iftyp.eq.5) then
                   call rewmf
                   call srchmf(lstr,'[GEOMETRIES]',istat)
                else
                   call rewfil
                endif
                if (iftyp.eq.5.and.ihasg.eq.2) then
                   call getzm(iatoms,0,0,istat)
                   if (istat.eq.1) then
                      dummy = .true.
                   else
                      dummy = .false.
                   endif
                else
                   if (iftyp.eq.11.or.iftyp.eq.12) then
                      call gettnk(igttnk,0,ipdbon,iff,iheat,heat)
                      if (igttnk.eq.1) then
                         dummy = .true.
                      else
                         dummy = .false.
                      endif
                   elseif (iftyp.eq.13) then
                      call getmol(igetmo,0)
                      if (igetmo.eq.1) then
                         dummy = .true.
                      else
                         dummy = .false.
                      endif
                   else
                      call getxyz(igetxy,heat,0)
                      if (igetxy.eq.1) then
                          dummy = .true.
                      else
                          dummy = .false.
                      endif
                   endif
                endif
                if (dummy) ipnt = 1
             endif
          elseif (icomm.gt.0.or.
     &           (icomm.eq.next.and.iftyp.eq.5.and.ifrav.eq.1)) then
             if (icomm.eq.next) then
                ipntt = ipnt + 1
                if (ipntt.gt.ipoints) ipntt = ipoints
             endif
             if (icomm.gt.0) ipntt = icomm
             if (iftyp.eq.5) then
                call rewmf
                call srchmf(lstr,'[GEOMETRIES]',istat)
             else
                call rewfil
             endif
             do i=1,ipntt
                if (iftyp.eq.5.and.ihasg.eq.2) then
                   call getzm(iatoms,0,0,istat)
                   if (istat.eq.1) then
                      dummy = .true.
                   else
                      dummy = .false.
                   endif
                else
                   if (iftyp.eq.11.or.iftyp.eq.12) then
                      if (iff.eq.7) then

c                        speed things up ambfor MD trajectories

                         if (i.ne.ipntt) then
                            call gtheat(igttnk,iheat,heat)
                         else
                            call search(lstr,'[AMBFOR]',istat)
                            if (istat.eq.1) then
                                call bckfil
                                call gettnk(
     &                               igttnk,0,ipdbon,iff,iheat,heat)
                            else
                                igttnk = 0
                            endif
                         endif
                      else
                         call gettnk(igttnk,0,ipdbon,iff,iheat,heat)
                      endif
                      if (igttnk.eq.1) then
                         dummy = .true.
                      else
                         dummy = .false.
                      endif
                   elseif (iftyp.eq.13) then
                      call getmol(igetmo,0)
                      if (igetmo.eq.1) then
                         dummy = .true.
                      else
                         dummy = .false.
                      endif
                   else
                      call getxyz(igetxy,heat,0)
                      if (igetxy.eq.1) then
                          dummy = .true.
                      else
                          dummy = .false.
                      endif
                   endif
                endif
             end do
             if (iftyp.eq.5) then
                call rdfc(ipntt,istat)
                if (istat.eq.1) then
                   ifav = 1
                else
                   ifav = 0
                endif
             endif
             ipnt = ipntt
          endif
      endif

      if (iftyp.eq.8) then
          if ( icomm .eq. first ) then
             call rewmf
             ipnt = 1
          elseif ( icomm .eq. next ) then
             if ( ipnt.ge.ipoints ) return
             ipnt = ipnt +1
          elseif ( icomm .eq. last ) then
             call rewmf
             ipnt = ipoints
          elseif (icomm.gt.0) then
             call rewmf
             ipnt = icomm
          endif
          call qcxyz(0,.false.,ipnt,istat)
          if (istat.eq.-1) then
               ifav = 0
          else
               ifav = 1
          endif
      endif

      if (iftyp.eq.9) then
          if ( icomm .eq. first ) then
             call rewmf
             ipnt = 1
          elseif ( icomm .eq. next ) then
             if ( ipnt.ge.ipoints ) return
             ipnt = ipnt +1
          elseif ( icomm .eq. last ) then
             call rewmf
             ipnt = ipoints
          elseif (icomm.gt.0) then
             call rewmf
             ipnt = icomm
          endif
          call orcxyz(0,ipnt,istat)
          if (istat.eq.-1) then
               ifav = 0
          else
               ifav = 1
          endif
      endif

      if (iftyp.eq.15) then
          if ( icomm .eq. first ) then
             if (iecce.eq.1) then
                call rewfil
             else
                call rewmf
             endif
             ipnt = 1
          elseif ( icomm .eq. next ) then
             if ( ipnt.ge.ipoints ) return
             ipnt = ipnt +1
          elseif ( icomm .eq. last ) then
             if (iecce.eq.1) then
                call rewfil
             else
                call rewmf
             endif
             ipnt = ipoints
          elseif (icomm.gt.0) then
             if (iecce.eq.1) then
                call rewfil
             else
                call rewmf
             endif
             ipnt = icomm
          endif
          if (iecce.eq.1) then
             call enwxyz(0,ipnt,istat)
          else
             call nwxyz(0,ipnt,istat)
          endif
      endif

      if (iftyp.eq.2.and.ipoints.gt.0) then
          if ( icomm .eq. first ) then
             call rewfil
             ipnt = 1
          elseif ( icomm .eq. next ) then
             if ( ipnt.ge.ipoints ) return
             ipnt = ipnt +1
          elseif ( icomm .eq. last ) then
             call rewfil
             ipnt = ipoints
          elseif (icomm.gt.0) then
             call rewfil
             ipnt = icomm
          endif
          if (ioxyz.eq.1) then
             itmp = ipnt - 1
             call gampoi(itmp,istat,ioxyz)
          else
             call gampoi(ipnt,istat,ioxyz)
          endif
          if ( istat .eq. -1 ) then
               ifav = 0
          else
               ifav = 1
          endif
c          if (msucc.eq.0.and.ioxyz.eq.0) ifav = 0
      endif

      if (iftyp.eq.4.and.ipoints.gt.1) then
          if ( icomm .eq. first ) then
             call rewfil
             ipnt = 1
          elseif ( icomm .eq. next ) then
             if ( ipnt.ge.ipoints ) return
             ipnt = ipnt +1
          elseif ( icomm .eq. last ) then
             call rewfil
             do i=1,ipoints
                 if ((nzm.eq.nzo.or.2*nzm.eq.nzo).and.nzm.gt.0) then
                    call search(lstr,'Z-MATRIX (ANGSTROMS',istat)
                 elseif (nzo.gt.0) then
                    call search(lstr,'Z-Matrix orientation:',istat)
                 else
                    call search(lstr,stmp(1:ns),istat)
                 endif
             end do
             call bckfil
             ipnt = ipoints
          elseif (icomm.gt.0) then
             call rewfil
             do i=1,icomm
                 if ((nzm.eq.nzo.or.2*nzm.eq.nzo).and.nzm.gt.0) then
                    call search(lstr,'Z-MATRIX (ANGSTROMS',istat)
                 elseif (nzo.gt.0) then
                    call search(lstr,'Z-Matrix orientation:',istat)
                 else
                    call search(lstr,stmp(1:ns),istat)
                 endif
             end do
             call bckfil
             ipnt = icomm
          endif
          call gaupoi(istat)
          if ( istat .eq. -1 ) then
               ifav = 0
          else
               ifav = 1
          endif
      endif

      if (iftyp.eq.3.and..not.(ipoints.le.1.and.iatoms.gt.0) ) then
          if ( icomm .eq. first ) then
             call rewfil
             ipnt = 1
             call gamupt(ipnt,istat)
          elseif ( icomm .eq. next ) then
             if ( ipnt.ge.ipoints ) return
             ipnt = ipnt +1
             call gamupt(ipnt,istat)
          elseif ( icomm .eq. last ) then
             call rewfil
             do i=1,ipoints
                call gamupt(ipoints,istat)
             end do
             call bckfil
             ipnt = ipoints
          elseif (icomm.gt.0) then
             call rewfil
             do i=1,icomm
                call gamupt(icomm,istat)
             end do
             call bckfil
             ipnt = icomm
          endif
          ifav = 1
      endif

      if (iftyp.eq.7.and..not.(ipoints.le.1.and.iatoms.gt.0)
     &     .and.irtype.ne.3 ) then
c          iunt = iun2
c          iun2 = iun5
          if ( icomm .eq. first ) then
             call rewfil
             ipnt = 1
             call cpmdpt(ipnt,istat)
          elseif ( icomm .eq. next ) then
             if ( ipnt.ge.ipoints ) return
             ipnt = ipnt +1
             call cpmdpt(ipnt,istat)
          elseif ( icomm .eq. last ) then
             call rewfil
             do i=1,ipoints
                call cpmdpt(ipoints,istat)
             end do
             call bckfil
             ipnt = ipoints
          elseif (icomm.gt.0) then
             call rewfil
             do i=1,icomm
                call cpmdpt(icomm,istat)
             end do
             call bckfil
             ipnt = icomm
          endif
          ifav = 1
c          iun2 = iunt
      endif

      if (iftyp.eq.7.and..not.(ipoints.le.1.and.iatoms.gt.0)
     &     .and.irtype.eq.3 ) then
          if (icomm.ne.first) then
             iunt = iun2
             iun2 = iun5
          endif
          if ( icomm .eq. first ) then
             call rewfil
             ipnt = 1
             call cpmdpt(ipnt,istat)
          elseif ( icomm .eq. next ) then
             if ( ipnt.ge.ipoints ) return
             ipnt = ipnt +1
             call cpmdptdyn(ipnt,istat)
          elseif ( icomm .eq. last ) then
             call rewfil
             do i=1,ipoints
                call cpmdptdyn(ipoints,istat)
             end do
             call bckfil
             ipnt = ipoints
          elseif (icomm.gt.0) then
             call rewfil
             do i=1,icomm
                call cpmdptdyn(icomm,istat)
             end do
             call bckfil
             ipnt = icomm
          endif
          ifav = 1
          if (icomm.ne.first) iun2 = iunt
      endif

      if (((icrtp.le.3.and.icrtp.ge.1).or.icrtp.eq.5)
     &     .and.ipoints.gt.1) then
          if (icomm.eq.first) then
             call rewfil
             if (icrtp.eq.3) then
                call rfbio(0,1,istat)
             elseif (icrtp.eq.1) then
                call rdchx(0,3,0,0,0,istat,1)
             elseif (icrtp.eq.5) then
                call getxdt(1,natc,istat)
             else
                call rfdat(0,istat,refcod)
             endif
             ipnt = 1
          elseif ( icomm .eq. next ) then
             if ( ipnt.ge.ipoints ) return
             ipnt = ipnt +1
             if (icrtp.eq.3) then
                call rfbio(0,0,istat)
             elseif (icrtp.eq.1) then
                call rdchx(0,3,0,0,0,istat,1)
             elseif (icrtp.eq.5) then
                call getxdt(ipnt,natc,istat)
             else
                call rfdat(0,istat,refcod)
             endif
          elseif (icomm.gt.0) then
             call rewfil
             if (icrtp.eq.3) then
                call rfbio(0,1,istat)
             elseif (icrtp.eq.1) then
                call rdchx(0,3,0,0,0,istat,1)
             elseif (icrtp.eq.5) then
                call getxdt(icomm,natc,istat)
             else
                call rfdat(0,istat,refcod)
             endif
             if (icrtp.ne.5) then
                do i=1,icomm-1
                   if (icrtp.eq.3) then
                      call rfbio(0,0,istat)
                   elseif (icrtp.eq.1) then
                      call rdchx(0,3,0,0,0,istat,1)
                   else
                      call rfdat(0,istat,refcod)
                   endif
                end do
             endif
             ipnt = icomm
          endif
          if (.not.(icrtp.eq.3.and.istat.eq.1)) 
     &       call fdat(ifd,0,0,0,0,0)
      endif

      if (ipdbon.eq.1.and.ipdbgro.eq.1) then
          if (icomm.eq.first) then
             ipnt = 1
             call gtfrm(ipnt)
          elseif ( icomm .eq. next ) then
             if ( ipnt.ge.ipoints ) return
             ipnt = ipnt + 1
             call gtfrm(ipnt)
          elseif (icomm.gt.0) then
             call rewfil
             call gtfrm(icomm)
             ipnt = icomm
          endif
      endif

100   continue

      if (ioadd.eq.0) then
        if (.not.(
     &    ( (iftyp.eq.2.or.iftyp.eq.3).and.ipoints.eq.0
     &      .and.(ioxyz.eq.1.or.irtype.eq.4) ).or.
     &    (iftyp.eq.4.and.ipoints.le.1).or.
     &    (icrtp.eq.3.and.ipoints.gt.1).or.
     &    (ipdbon.eq.1.and.dozme))) call docent
      endif

      if (iftyp.eq.2.and.ipoints.gt.0.and.ioxyz.eq.0.and.msucc.eq.1
     &    .and.(.not.dozme)) 
     &    call rotmom(ipnt,ifav)


      ibeg = 1
      if (ioadd.eq.1) ibeg = ioatms + 1

      do i=ibeg,iatoms
         if (ipdbgro.eq.1) then
            if (ianz(i).eq.100.and.ipnt.gt.1) then
               if (i.gt.iclpnt(1)+7) iaton(i) = 0
            else
               iaton(i) = 1
            endif
         else
            iaton(i) = 1
         endif
      end do

      if (ipdbon.eq.1) then
         if (ialtyp.eq.1) then
            if (dozme) then
               call convar(coo,iconn,ianz,nscnd,iscst)
               call pmfass(1,0)
               call totpmf(dum(1))
               call upsco()
            endif
         else
            if (dozme) then
               if (idoprt.eq.0) call conpdb
               call chkbck(0)
            else
               if (((iftyp.eq.11.or.iftyp.eq.12).and.iff.eq.7)
     &             .or.ipdbgro.eq.1) call chkbck(1)
            endif
         endif
      else
         if ((.not.(iftyp.eq.6.or.ichx.ne.0).and.
     &        .not.(iftyp.eq.5.and.ihasg.eq.1).and.
     &        .not.(iftyp.ge.10.and.iftyp.lt.15)
     &       ).or.(dozme.and..not.ichx.eq.1)) call doconn
      endif

      if (ipdbon.eq.0) then
         if (.not.ichx.eq.1.and.iff.ne.7) call dohcon(0)
      else
         call domcon(1,1)
      endif
c
c  get scale
c

      if (ioadd.eq.0) then
         if (doscl.eq.1.and..not.(icrtp.eq.3.and.ipoints.gt.1).and.
     &       .not.(ipdbon.eq.1.and.dozme)) then
            call doscal
            if (dozme.and.iatoms.lt.7) scali = 3.5d0
            scal = scali * 2.4d0 * smag
         endif
      endif

      if (ifav.eq.1) call parfc

      call drwgeo

      if (movie.eq.0.or.(movie.ne.0.and.ipnt.eq.ipoints)) then
         if (ctoz) call allzmt(ipdbon)
         if (.not.dozme) call upzme
      endif

      dozme = .false.


      if (icrtp.eq.2.and.ipoints.gt.1) then
          lstr = refcod
          call inferr(lstr,0)
cd         write(iun3,'(a)')'leave subroutine getpoi'
          return
      elseif (icrtp.eq.5.and.ipoints.gt.1) then
          call zerstr(ipnt,astr,6,1)
          lstr = 'point '//astr
          call inferr(lstr,0)
          return
      endif

      if (.not.( iftyp.eq.7.or.
     &          (iftyp.ge.2.and.iftyp.le.4).or.
     &          (iftyp.eq.1.and.iconv.eq.0).or.
     &          (iftyp.ge.10.and.ipoints.gt.1)
     &         ).or.dozme) then
cd        write(iun3,'(a)')'leave subroutine getpoi'
         return
      endif

      if (ipoints.eq.0) then
         lstr = ' '
         if (iftyp.ne.1) then
             lstr = 'Single point'
         endif
      elseif (ipnt.eq.1) then
         lstr = 'First point'
      elseif (ipnt.eq.ipoints) then
         lstr = 'Last point'
      else
         if (ioxyz.eq.1) then
            call zerstr(ipnt-1,astr,6,1)
         else
            call zerstr(ipnt,astr,6,1)
         endif
         lstr = 'point '//astr
      endif

      if (iftyp.eq.7) lstr=lstr(1:22)//'CPMD output file'

      call inferr(lstr,0)

cd     write(iun3,'(a)')'leave subroutine getpoi'
      return
      end

      subroutine docond(coo,iconn,ianz)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      dimension coo(3,*),ianz(*),iconn(mxcon+1,*)

      call convar(coo,iconn,ianz,iatoms,0)

      return
      end

      subroutine convar(coo,iconn,ianz,iatoms,ioff)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (mxel=100)
      parameter (mxcon=10)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      dimension coo(3,*),iconn(mxcon+1,*),ianz(*),ctemp(3)

      adjus=0.52917706d0

      do i=ioff+1,ioff+iatoms

         k = 0
         iconn(1,i) = 0
         iani = ianz(i)

         do j=ioff+1,ioff+iatoms

            ianj = ianz(j)

            if (iani.gt.0.and.ianj.gt.0.and.
     &          iani.le.100.and.ianj.le.100) then

               dmaxsq = (vdwr(iani) + vdwr(ianj))/adjus
               dmaxsq = dmaxsq*dmaxsq

               dijsq = 0
               do l=1,3
                  ctemp(l) = coo(l,i) - coo(l,j)
                  dijsq = dijsq + ctemp(l)*ctemp(l)
               end do

               if (i.ne.j.and.iani.ne.99.and.ianj.ne.99) then
                 if (dijsq.lt.dmaxsq) then
                     k = k + 1
                     if (k.le.mxcon) then
                        iconn(1,i) = k
                        iconn(k+1,i) = j
                     else
                        write(iun3,*) 
     &                     'more than mxconn connections found'
                     endif
                 endif
               endif

            endif

         end do
      end do

      return
      end

      subroutine dohcod(ihflag,coo,ianz,iconn,isurf)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxel=100)
      parameter (nacc=6)
      parameter (ndon=4)
      common /athlp/ iatoms, mxnat
      logical hon
      common /hcon/ iacc(nacc),idon(ndon),hdmin,hdmax,hamin,hamax,hon
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      logical isacc,isdon
      real angle,hamin,hamax
      dimension ctemp(3),isel(4)
      dimension coo(3,*),ianz(*),iconn(mxcon+1,*),isurf(*)

      call domcon(1,1)

      if (.not.hon) return

      hdmin2 = hdmin*hdmin
      hdmax2 = hdmax*hdmax

c
c Hydrogen bonds
c
      do i=1,iatoms
         do j=i+1,iatoms
           if (ihflag.eq.0.or.
     &     (ihflag.gt.0.and.(isurf(i).eq.1.or.isurf(j).eq.1))) then
            iai = ianz(i)
            iaj = ianz(j)
            dijsq = 0
            do l=1,3
               ctemp(l) = coo(l,i) - coo(l,j)
               dijsq = dijsq + ctemp(l)*ctemp(l)
            end do
            if ((iai.eq.1.and.isacc(iaj)).or.
     &          (iaj.eq.1.and.isacc(iai))) then
               if (dijsq.gt.hdmin2.and.dijsq.lt.hdmax2) then
                  isel(1) = 0
                  if (iai.eq.1) then
                      if (iconn(1,i).gt.0) then
                         iai = ianz(iconn(2,i))
                         if (isdon(iai)) then
                            isel(1) = iconn(2,i)
                            isel(2) = i
                            isel(3) = j
                         endif
                      endif
                  else
                      if (iconn(1,j).gt.0) then
                         iaj = ianz(iconn(2,j))
                         if (isdon(iaj)) then
                            isel(1) = iconn(2,j)
                            isel(2) = j
                            isel(3) = i
                         endif
                      endif
                  endif
                  if (isel(1).ne.0) then
                    call intcor(intc,angle,isel,3)
                    if (intc.eq.1) then
                     if (abs(angle).gt.hamin.and.
     &                    abs(angle).lt.hamax) then
                         k = iconn(1,i)
                         if (k.lt.mxcon) then
                            iconn(1,i) = iconn(1,i) + 1
                            iconn(iconn(1,i)+1,i) = -j
                         else
                            write(iun3,*) 
     &                        'more than mxconn connections found'
                         endif
                         k = iconn(1,j)
                         if (k.lt.mxcon) then
                            iconn(1,j) = iconn(1,j) + 1
                            iconn(iconn(1,j)+1,j) = -i
                         else
                            write(iun3,*) 
     &                        'more than mxconn connections found'
                         endif
                     endif
                    endif
                  endif
               endif
            endif
           endif
         end do
      end do

      if (ihflag.eq.1) call clrsrf

      return
      end

      subroutine nohcod(iconn)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      dimension icn(mxcon+1),iconn(mxcon+1,*)

c
c Deactivate Hydrogen bonds
c
      do i=1,iatoms
          inn = 0
          do j=1,iconn(1,i)
             if (iconn(1+j,i).gt.0) then
                inn = inn + 1
                icn(inn) = iconn(1+j,i)
             endif
          end do  
          iconn(1,i) = inn
          do j=1,inn
             iconn(1+j,i) = icn(j)
          end do
      end do

      return
      end

      logical function isacc(iatnr)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (nacc=6)
      parameter (ndon=4)
      logical hon
      common /hcon/ iacc(nacc),idon(ndon),hdmin,hdmax,hamin,hamax,hon
      real hamin,hamax

      isacc = .false.

      do i=1,nacc
          if (iatnr.eq.iacc(i)) isacc = .true.
      end do

      return
      end

      logical function isdon(iatnr)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (nacc=6)
      parameter (ndon=4)
      logical hon
      common /hcon/ iacc(nacc),idon(ndon),hdmin,hdmax,hamin,hamax,hon
      real hamin,hamax

      isdon = .false.

      do i=1,ndon
          if (iatnr.eq.idon(i)) isdon = .true.
      end do

      return
      end

      subroutine doscad(t,coo,zv,pincr,scal,scali,smag,iscupd)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /athlp/ iatoms, mxnat
      common /gracom/  uscl,colscd,colscpd,ivdwpl
      dimension t(3),coo(3,*)

      scali = 0.0d0
      do i=1,iatoms
         dist = 0.0d0
         do j=1,3
             dt = coo(j,i)-t(j)
             dist = dist + dt*dt
         end do
         if (dist.gt.scali) scali = dist
      end do
      if (scali.eq.0.0d0) then
         scali = 1.0d0
      else
         scali = dsqrt(scali)
      endif
      scal = scali * 2.4d0 * smag

      if (scali.gt.15.0d0) then
         colscpd = 1.3d0
      else
         colscpd = .4d0
      endif

      colscd = 2.0d0

      if (iscupd.eq.1) then
         zv = scali
         iscupd = 0
      endif

      pincr = 0.02d0*scali

      return
      end

      subroutine redcod(iconn)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      dimension icnn(mxcon),iconn(mxcon+1,*)

      do j=1,iatoms
         ibnds = 0
         do i=1,iconn(1,j)
            ioke = 1
            do k=1,i-1
               if (iconn(k+1,j).eq.iconn(i+1,j)) ioke = 0
            end do
            if (ioke.eq.1.and.ibnds.lt.mxcon) then
               ibnds = ibnds + 1
               icnn(ibnds) = iconn(i+1,j)
            endif
         end do
         iconn(1,j) = ibnds
         do i=1,ibnds
            iconn(i+1,j) = icnn(i)
         end do
      end do

      return
      end

      subroutine domcod(ipt,idoall,monmod,iconn)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (maxdm=20)
      common /distmn/ rdm(maxdm),idmon(2,maxdm),ndm
      integer srfmap, srfloc
      common /hlpsrf/ npts1,npts2,npts3,srfmap,srfloc,ifogl,itsrf
      real rout
      dimension iconn(mxcon+1,*)


      if (idoall.eq.1) then
         ibeg = 1
         iend = ndm
      else
         ibeg = ipt
         iend = ipt
      endif

      do i=ibeg,iend
         if (monmod.eq.2) then
            call intcor(intc,rout,idmon(1,i),2)
         else
            intc = 1
         endif
         if (intc.eq.1) then
            if (monmod.eq.2) rdm(i) = dble(rout)*0.52917706d0
            if (ifogl.ne.1) then
               i1 = idmon(1,i)
               i2 = idmon(2,i)
               k = iconn(1,i1)
               if (k.lt.mxcon) then
                  iconn(1,i1) = k + 1
                  iconn(k+2,i1) = -i2
               endif
               k = iconn(1,i2)
               if (k.lt.mxcon) then
                  iconn(1,i2) = k + 1
                  iconn(k+2,i2) = -i1
               endif
            endif
         endif
      end do

      call ogmon

      return
      end

      subroutine chkbcd(ncalf,ihet,reson,iaton,iatclr,ision)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      integer fancy,shade,atcol,dolabs,persp,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      common /cllab/  iclon,iclpnt(4)
      common /boxion/ ionon,ionclr
      parameter (mxheta=150)
      integer reson
      dimension ihet(*),reson(*),iaton(*),iatclr(*)

      if (backb.eq.1) then
         call actcal(1)
         call ribbs
         if (ihet(1).eq.1) call acthel(1,0,0,0)
         if (ihet(2).eq.1) call acthel(1,1,0,0)
         if (ihet(3).eq.1) call acthel(1,2,0,0)
         if (ihet(4).eq.1) call acthel(1,3,0,0)
         do i=5,mxheta
            j = -i+1
            if (ihet(i).eq.1) then
               call actami(j,0,1,0)
            endif
         end do
         if (ionon.eq.0) then
            call actami(ision,ionclr,1,0)
         endif
         do i=1,ncalf
            if (reson(i).eq.1) then
               call actami(i,0,1,0)
            endif
         end do
         if (iclon.eq.1) then
            icllow = iclpnt(1) 
            do i=0,7
               iaton(icllow+i) = 1
               iatclr(icllow+i) = 11
            end do
         endif
      endif

      return
      end
