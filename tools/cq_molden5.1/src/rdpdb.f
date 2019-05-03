      subroutine rdpdd(istat,coo,ianz,ihashb)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      common /athlp/ iatoms, mxnat
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character*2 elemnt
      common /elem/elemnt(mxel)
      character lstr*137
      character*2 tstr,tocapf
      dimension coo(3,*),ianz(*)

      ihashb = 0

      rewind iun2
      istat = 1
      
10    call searchd(lstr,'ATOM','HETATM',istat)
      if (istat.eq.0) then
           call inferr('no ATOM/HETATM found !',0)
           goto 100
      endif
      if (lstr(1:4).ne.'ATOM'.and.lstr(1:6).ne.'HETATM') goto 10
      backspace iun2

      iatoms = 0
      do while (.true.)
15        read(iun2,'(a)',end=20) lstr
          if (lstr(1:4).ne.'ATOM'.and.lstr(1:6).ne.'HETATM') goto 15
          if (iatoms.ge.mxnat) then
               call inferr('exceeded maximum atoms !',1)
               goto 100
          endif
          iatoms = iatoms + 1
          ianz(iatoms) = 0
          read(lstr,'(12x,a2,16x,3f8.3)') tstr,
     &        coo(1,iatoms),coo(2,iatoms),coo(3,iatoms)
          it = ichar(tstr(1:1))
          if (.not.(it.ge.65.and.it.le.90).and..not.
     &             (it.ge.97.and.it.le.122)) tstr(1:1) = ' '
          do j=1,99
             if (tstr.eq.tocapf(elemnt(j))) ianz(iatoms) = j
          end do
          if (ianz(iatoms).eq.0.and.tstr(1:1).ne.' ') then
             tstr(1:1) = ' '
             do j=1,99
                if (tstr.eq.tocapf(elemnt(j))) ianz(iatoms) = j
             end do
          endif
          if (ianz(iatoms).le.0) then
             if (tstr(2:2).eq.'D') then
                ianz(iatoms) = 1
             else
                write(iun3,*) 'Unclassified atom =',tstr
                write(iun3,*) lstr
                ianz(iatoms) = 99
             endif
          endif
      end do

20    continue

      return
100   istat = 0
      return
      end

      subroutine pdbsiz(iszhnt)
      implicit double precision (a-h,o-z)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character lstr*137
      rewind iun2

      iszhnt = 0
      do while (.true.)
         read(iun2,'(a)',end=20) lstr
         if (lstr(1:4).eq."ATOM".or.lstr(1:6).eq."HETATM") 
     &       iszhnt = iszhnt + 1
      end do

20    iszhnt = 3*iszhnt

      return
      end

      subroutine pdbsdd(istat,doscnd,ioadd,
     &                  coo,ianz,iaton,iresid,iconn,ityp,
     &                  ncalf,ianf,islu,nchain,iamino,ihet,
     &                  isal,irsnr,achain,ihashb,ishoh)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxchtp=136)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxamo=201)
      integer*2 ityp
      character*3 chtnk,pdbsym,hsym,ambtnk
      character*2 amotnk,gffstr
      common /symbol/ pdbsym(mxsym),hsym(mxhsym),chtnk(mxchtp),
     &                ambtnk(mxamb),amotnk(mxamo),gffstr(mxgff)
      common /athlp/ iatoms, mxnat
      common /surf / natorg,nosncd

      parameter (numcal=50000)
      parameter (mxres=42)
      parameter (mxheta=150)

      parameter (mxmmul=100)
      character*1 achain
      logical isatm,ishet,hasalt,doscnd,skip,dall,goods
      character*3 aminos
      common /amino/ aminos(mxres)
      common /multim/imulm, nmulm,ihasqm(mxmmul)
      character*3 hetz
      common /clfstr/ ihashz,ihetq(mxheta),ihqset(mxheta),ihhadd(mxheta)
     &                ,hetz(mxheta)
      character*4 pdbc
      common /pdbcod/ pdbc

      parameter (mxel=100)
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character lstr*137
      character*2 tocapf,hetred
      character*4 tstr,tmpstr
      character*3 resi,resit
      character*1 restmp,resa,altloc,aloc,ach
      dimension ipdb(mxsym),hetred(4),
     &          ihpdb(mxhsym*3),tmp(3),isc(numcal)
      dimension coo(3,*),ianz(*),iaton(*),iresid(*),iconn(mxcon+1,*),
     &          ityp(*)
      dimension ianf(*),islu(*),iamino(*),ihet(*),isal(*),irsnr(*),
     &          achain(*)

      data hetred/
     &   'AC','AN','AO','AP'/

      imulm = 1
      if (ioadd.eq.0) then
         call parsfn('Helix',5,1)
         call parsfn('Beta',4,1)
         call parsfn('RNA/DNA',7,1)
         call parsfn('Coil',4,1)
         hetz(1) = 'Hel'
         hetz(2) = 'Bet'
         hetz(3) = 'DNA'
         hetz(4) = 'Coi'
      endif

      ihashb = 0


      ioatms = 0
      if (ioadd.eq.1) then
         natorg = iatoms
         do i=ncalf+1,numcal
             isal(i) = 3
             isc(i) = 3
         end do
         noff = ncalf
         ioatms = iatoms
      else
         do i=1,numcal
             isal(i) = 3
             isc(i) = 3
         end do
         noff = 0
      endif

      if (.not.doscnd) goto 5

      rewind iun2
      call search(lstr,'HEADER',istat)
      if (istat.eq.1) then
         pdbc = lstr(63:66)
         call tolow(pdbc,4)
         call parstr(pdbc,11)
      else
         rewind iun2
      endif
c
c     Proces Helix/Sheet information
c
      rewind iun2

      do while (.true.)
          call searchd(lstr,'HELIX','SHEET',istat)
          if (istat.eq.0) goto 5
          if (index(lstr,'HELIX').eq.1) then
             read(lstr,'(21x,i4,8x,i4)',err=4,end=4) i1,i2
             ihashb = 1
             if (noff+i1.gt.numcal) i1 = numcal
             if (noff+i2.gt.numcal) i2 = numcal
             i3 = i2 - i1 + 1
             if (i1.ne.0.and.i3.gt.1) then
                do i=i1,i2
                   isc(noff+i) = 0
                end do
             endif
          elseif (index(lstr,'SHEET').eq.1) then
             read(lstr,'(22x,i4,7x,i4)',err=4,end=4) i1,i2
             ihashb = 1
             if (noff+i1.gt.numcal) i1 = numcal
             if (noff+i2.gt.numcal) i2 = numcal
             i3 = i2 - i1 + 1
             if (i1.ne.0.and.i3.gt.1) then
                do i=i1,i2
                   isc(noff+i) = 1
                end do
             endif
          endif
      end do
4     call inferr('Error reading HELIX/SHEET record !',0)
5     continue

      rewind iun2
      istat = 1
      hasalt = .false.
      dall = .false.
      
10    call searchd(lstr,'HETATM','ATOM',istat)
      if (istat.eq.0) then
           call inferr('no ATOM/HETATM record found !',0)
           rewind iun2
           goto 210
      endif

      if (lstr(1:4).ne.'ATOM'.and.lstr(1:6).ne.'HETATM') goto 10
      backspace iun2

      if (ioadd.eq.1) then
         call numhet(nhmol)
         nmulm = nhmol - 3
         nmulmt = nmulm
         icres = ncalf
         ianf(nchain+1) = 1
         islu(nchain+1) = 0
      else
         nhmol = 3
         iatoms = 0
         icres = 0
         nchain = 1
         ianf(1) = 1
         islu(1) = 0
         nmulm = 0
      endif

      jres = -1
      irtmp = -10000
      ipdb(38) = 0
      restmp = ' '
      ihashy = 0
      skip = .false.
      aloc = ' '

      do while (.true.)
15        read(iun2,'(a)',end=20) lstr
          isatm = (lstr(1:4).eq.'ATOM'.or.lstr(1:6).eq.'HETATM')
          if (lstr(1:6).eq.'ENDMDL') goto 20
          if (.not.isatm) goto 15
c
c found atom
c
          read(lstr,'(12x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3)') tstr,altloc,
     &       resi,ach,iresnr,resa,(tmp(l),l=1,3)
           
          call tocap(tstr,4)
          call tocap(resi,3)

          jrest = 0
          if (resi(2:3).eq.'  ') resi = '  '//resi(1:1)

          if (resi(2:2).eq.'D'.and.resi(1:1).eq.' ') then
             if (resi(3:3).eq.'A'.or.resi(3:3).eq.'C'.or.
     &           resi(3:3).eq.'G'.or.resi(3:3).eq.'T') then
                 resi(2:2) = ' '
             endif
          endif

          do j=1,mxres
             if (resi.eq.aminos(j)) jrest = j
          end do
          if (resi.eq.'HIP'.or.resi.eq.'HID'.or.
     &            resi.eq.'HIE') then
             jrest = 17
             resi = 'HIS'
          endif
          if (jrest.eq.0) goto 15

          do i=1,4
             if (ichar(tstr(i:i)).eq.39) tstr(i:i) = '*'
          end do
          if (tstr(3:4).eq.'**'.and.tstr(2:2).ne.'O') then
             tstr(2:4) = tstr(1:3)
             tstr(1:1) = ' '
          endif

          if (tstr(2:4).eq.'OT ') tstr(2:4) = 'OXT'
          if (tstr(2:4).eq.'OCT') tstr(2:4) = 'OXT'
          if (tstr(2:4).eq.'OP2') tstr(2:4) = 'O2P'
          if (tstr(2:4).eq.'OP1') tstr(2:4) = 'O1P'

          if (iresnr.ne.irtmp.or.restmp.ne.resa) then
c
c next amino acid
c
             if (iatoms.gt.mxnat-30) goto 20
             if (skip) then
                skip = .false.
                if (iresnr.eq.irtmp.and.restmp.ne.resa) skip = .true.
             else
c                if (icres.gt.0.and.jres.ne.0) then
                if (icres.gt.0.and.jres.gt.0) then
                   call mkcon(ipdb,jres,icres,ihpdb,ihashy,0)
                   call pdbtyp(ipdb,ihpdb,jres,1)
                endif
                if (jres.ne.0) icres  = icres + 1
                if (restmp.ne.resa) skip = .true. 
             endif
             if (icres.gt.0.and.icres.le.numcal) then
                iamino(icres) = jrest
                irsnr(icres) = iresnr
                achain(icres) = ach
                if (ihashb.eq.1.and.iresnr.gt.0.and.iresnr.le.numcal) 
     &              isal(icres) = isc(iresnr)
             endif
             irtmp = iresnr
             restmp = resa
             jres = jrest
             do j=1,mxsym
                ipdb(j) = 0
             end do
             do j=1,mxhsym*3
                ihpdb(j) = 0
             end do
             ihashy = 0
             if (skip) goto 15
          else
             if (skip) goto 15
          endif

          if (iatoms.ge.mxnat-1000) then
             dall = .true.
             goto 100
          endif

          iatoms = iatoms + 1
          ianz(iatoms) = 0
          iconn(1,iatoms) = 0
          iresid(iatoms) = icres

          do l=1,3
             coo(l,iatoms) = tmp(l)
          end do

c determine ianz, see if hydrogen atom

          call detanz(resi,tstr,lstr,ifnd,ish,ianz(iatoms),ihashy)

          ihstyp = 0
          if (ish.eq.1) then
            do j=1,mxhsym
               if (tstr(2:4).eq.hsym(j)) then
                  ihstyp = 1
                  if (ihpdb((j-1)*3+1).eq.0) then
                     aloc = altloc
                  elseif (altloc.ne.aloc) then
                     hasalt = .true.
                     iatoms = iatoms - 1
                     goto 15
                  endif
                  do l=1,3
                     if (ihpdb((j-1)*3+l).eq.0) then
                        ihpdb((j-1)*3+l) = iatoms
                        goto 1000
                     endif
                  end do
1000              continue
               endif
            end do
            if (tstr(2:4).eq.'HG3'.and.ihpdb(16).ne.0.and.
     &          ihpdb(13).eq.0) then
                ihstyp = 1
                ihpdb(13) = ihpdb(16)
                ihpdb(16) = iatoms
            endif
            if (tstr(2:4).eq.'HD3'.and.ihpdb(25).ne.0.and.
     &          ihpdb(22).eq.0) then
                ihstyp = 1
                ihpdb(22) = ihpdb(25)
                ihpdb(25) = iatoms
            endif
            if (tstr(2:4).eq.'HE3'.and.ihpdb(37).ne.0.and.
     &          ihpdb(31).eq.0) then
                ihstyp = 1
                ihpdb(31) = ihpdb(37)
                ihpdb(37) = iatoms
            endif
          else
            do j=1,mxsym
               if (tstr(2:4).eq.pdbsym(j)) then
                  ihstyp = 1
c
c Allow only the first of alternate locations
c
                  if (ipdb(j).eq.0) then
                     ipdb(j) = iatoms
                  else
                     hasalt = .true.
                     iatoms = iatoms - 1
                     goto 15
                  endif
               endif
            end do
          endif
          if (ihstyp.eq.0) then
              nhmol = nhmol + 1
              iresid(iatoms) = -nhmol
              irtmp = -10000
              resi = 'UNK'
              call parsfn(resi,3,1)
              hetz(nhmol+1) = 'UNK'
          endif
          iaton(iatoms)  = 1
      end do


20    continue

c round up last amino acid

      
      if (jres.gt.0) then
         call mkcon(ipdb,jres,icres,ihpdb,ihashy,0)
         call pdbtyp(ipdb,ihpdb,jres,1)
      endif
      islu(nchain) = icres
      ncalf = icres
      namatm = iatoms
c
c do hetatm
c
      ishoh = 0
      irtmp = -10000
      resit = '   '

      rewind iun2
210   call searchd(lstr,'HETATM','ATOM',istat)
      if (istat.eq.0) goto 220
      if (lstr(1:6).ne.'HETATM'.and.lstr(1:4).ne.'ATOM') goto 210
      backspace iun2

      do while (.true.)
215        read(iun2,'(a)',end=220) lstr
           ishet = (lstr(1:6).eq.'HETATM')
           isatm = (lstr(1:4).eq.'ATOM')
           if (lstr(1:3).eq.'TER') irtmp = -10000
           if (lstr(1:6).eq.'ENDMDL') goto 220
           resi = lstr(18:20)
           if (ishet.or.isatm) then
              jrest = 0
              ishet = .false.
              if (resi(2:3).eq.'  ') resi = '  '//resi(1:1)
              if (resi(2:2).eq.'D'.and.resi(1:1).eq.' ') then
                 if (resi(3:3).eq.'A'.or.resi(3:3).eq.'C'.or.
     &               resi(3:3).eq.'G'.or.resi(3:3).eq.'T') then
                     resi(2:2) = ' '
                 endif
              endif
              do j=1,mxres
                 if (resi.eq.aminos(j)) jrest = j
              end do
              if (resi.eq.'HIP'.or.resi.eq.'HID'.or.
     &            resi.eq.'HIE') then
                 jrest = 17
                 resi = 'HIS'
              endif
              if (jrest.eq.0) ishet = .true.
           endif
           if (.not.ishet) goto 215

           if (iatoms.ge.mxnat-1000) then
               dall = .true.
               goto 100
           endif

           iatoms = iatoms + 1
           ianz(iatoms) = 0
           read(lstr,'(12x,a4,a1,a3,2x,i4,4x,3f8.3)') tstr,altloc,
     &        resi,iresnr, coo(1,iatoms),coo(2,iatoms),coo(3,iatoms)
           
           call tocap(tstr,4)

           if (altloc.ne.' '.and.altloc.ne.'A'.and.altloc.ne.'a'.and.
     &         altloc.ne.'1') then
              iatoms = iatoms - 1
              hasalt = .true.
              goto 215
           endif

           if (resi.ne.'HOH') then
              tmpstr = tstr(1:4)
              isl = len(tmpstr)
              call spatrm(tmpstr,isl)
              call prslab(tmpstr,isl,iatoms)
           endif

           ic1 = 0
           ic2 = 0
           ls = len(lstr)
           goods = .false.
           if (ls.ge.78) then
              ic0 = ichar(lstr(76:76))
              ic1 = ichar(lstr(77:77))
              ic2 = ichar(lstr(78:78))
              goods = ( ((ic1.ge.65.and.ic1.le.90).or.ic1.eq.32) .and.
     &                  ((ic2.ge.65.and.ic2.le.90).or.ic2.eq.32) .and.
     &                  ic0.eq.32)
              if (goods) then
                  goods = (.not.(ic1.eq.32.and.ic2.eq.32))
              endif
           endif

           if (goods) then

              read(lstr,'(76x,a2)') tstr(1:2)

              do j=1,99
                if (tstr(1:2).eq.tocapf(elemnt(j))) ianz(iatoms) = j
              end do
            
           else

              if (resi.eq.'HEM') then
                 if (tstr(1:2).eq.'HA'.or.tstr(1:2).eq.'HB'.or.
     &               tstr(1:2).eq.'HM') tstr(1:4) = ' H  '
                 if (tstr(2:3).eq.'FE') tstr(1:3) = 'FE '
              endif

              if (resi.eq.'NDP') then
                 if (tstr(1:2).eq.'NP') tstr(1:4) = ' P  '
              endif

              do j=1,4
                 if (tocapf(tstr(1:2)).eq.hetred(j)) then
                     tstr(1:1) = ' '
                 endif
              end do

              it = ichar(tstr(1:1))
              if (.not.(it.ge.65.and.it.le.90).and..not.
     &           (it.ge.97.and.it.le.122)) tstr(1:1) = ' '
              it2 = ichar(tstr(2:2))
              if (.not.(it2.ge.65.and.it2.le.90).and..not.
     &           (it2.ge.97.and.it2.le.122)) tstr(2:2) = ' '
              if (tstr(2:2).eq.' '.and.tstr(1:1).ne.' ') then
                 tstr(2:2) = tstr(1:1)
                 tstr(1:1) = ' '
              endif 

              do j=1,99
                if (tstr(1:2).eq.tocapf(elemnt(j))) ianz(iatoms) = j
              end do

              if (ianz(iatoms).eq.0.and.tstr(1:1).ne.' ') then
                 tstr(1:1) = ' '
                 do j=1,99
                     if (tstr(1:2).eq.tocapf(elemnt(j))) 
     &                  ianz(iatoms) = j
                 end do
              endif

           endif

           if (ianz(iatoms).le.0) then
              write(iun3,*) 'Unclassified atom =',tstr
              write(iun3,*) lstr
              ianz(iatoms) = 99
           endif
           iconn(1,iatoms) = 0
           iaton(iatoms)  = 0
           if (iresnr.ne.irtmp.or.resi(1:3).ne.resit(1:3)) then
              irtmp = iresnr
              resit = resi
              if (resi(1:3).eq.'HOH'.or.resi(1:3).eq.'hoh') then
                  if (ishoh.eq.0) then
                     nhmol = nhmol + 1
                     if (nmulm.lt.mxmmul) nmulm = nmulm + 1
                     call parsfn(resi,3,1)
                     hetz(nhmol+1) = resi(1:3)
                     ishoh = nhmol
                  endif  
              else
                  nhmol = nhmol + 1
                  if (nmulm.lt.mxmmul) nmulm = nmulm + 1
                  call parsfn(resi,3,1)
                  hetz(nhmol+1) = resi(1:3)
              endif
           endif
           if (resi(1:3).eq.'HOH'.or.resi(1:3).eq.'hoh') then
              iresid(iatoms) = -ishoh
           else
              iresid(iatoms) = -nhmol
           endif
      end do
220   continue
c
c do connectivity between all except amino-amino
c
      nhet = iatoms - namatm
      do ii=1,nhet
         i = namatm + ii
         do j=ioatms+1,namatm + ii - 1
            call connij(idum,i,j,0)
         end do
      end do

c
c make s-s connections
c
      do i=1,ncalf
         if (iamino(i).eq.4) then
            do j=1,i-1
               if (iamino(j).eq.4) then
                  ii = 0
                  jj = 0
                  do k=1,iatoms
                     if (ianz(k).eq.16) then
                        if (iresid(k).eq.i) ii = k
                        if (iresid(k).eq.j) jj = k
                     endif
                  end do
                  if (ii.ne.0.and.jj.ne.0) then
                     call connij(idum,ii,jj,0)
                     if (idum.eq.1) then
                        ityp(ii) = 82
                        ityp(jj) = 82
                        do l=1,iconn(1,ii)
                           ll = iabs(iconn(l+1,ii))
                           if (ianz(ll).eq.6) ityp(ll) = 38
                        end do
                        do l=1,iconn(1,jj)
                           ll = iabs(iconn(l+1,jj))
                           if (ianz(ll).eq.6) ityp(ll) = 38
                        end do
                     endif
                  endif
               endif
            end do
         endif
      end do

      do i=1,ncalf
         if (iamino(i).eq.21) iamino(i) = 10
         if (iamino(i).eq.22) iamino(i) = 14
      end do
      if (islu(1).lt.ianf(1)) nchain = 0

      if (ioadd.eq.1) then
         do i=nmulmt,nmulm
            ihasqm(i) = 0
         end do
         if (nmulmt.eq.0) nmulmt = 1
         do i=nmulmt+4,nhmol+1
            ihet(i) = 1
         end do
      endif

      ihashz = 1

      if (hasalt) call inferr('Alternate location(s) detected',0)
      return

100   istat = 0
      if (dall) then
         iatoms = natorg
         ncalf = noff
         istat = -1
      else
         if (islu(1).lt.ianf(1)) nchain = 0
      endif
      if (hasalt) call inferr('Alternate location(s) detected',0)
      return
      end

      subroutine conpdd(ianz,iconn,iresid,ncalf,iamino)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      common /athlp/ iatoms, mxnat
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension iconn(mxcon+1,*),ianz(*),iresid(*),iamino(*)

      namatm = 0

      do i=1,iatoms
         iconn(1,i) = 0
         if (iresid(i).le.0.and.namatm.eq.0) namatm = i-1
      end do

      if (namatm.eq.0) namatm = iatoms

      do i=1,ncalf
         call getpdb(i,ipdb,ihpdb)
         call mkcon(ipdb,iamino(i),i,ihpdb,1,1)
      end do

c
c do connectivity between all except amino-amino
c
      nhet = iatoms - namatm

      do ii=1,nhet

         i = namatm + ii
         ir1 = iresid(i)

         do j=1,namatm + ii - 1
            ir2 = iresid(j)
            if (.not.(ir1.gt.0.and.ir2.gt.0)) then
               call connij(idum,i,j,1)
            endif
         end do

      end do

c
c make s-s connections
c
      do i=1,ncalf
         if (iamino(i).eq.4) then
            do j=1,i-1
               if (iamino(j).eq.4) then
                  ii = 0
                  jj = 0
                  do k=1,iatoms
                     if (ianz(k).eq.16) then
                        if (iresid(k).eq.i) ii = k
                        if (iresid(k).eq.j) jj = k
                     endif
                  end do
                  if (ii.ne.0.and.jj.ne.0) then
                     call connij(idcon,ii,jj,1)
                     if (idcon.eq.1) then
                        do l=1,iconn(1,ii)
                           ll = iabs(iconn(l+1,ii))
                        end do
                        do l=1,iconn(1,jj)
                           ll = iabs(iconn(l+1,jj))
                        end do
                     endif
                  endif
               endif
            end do
         endif
      end do

      return
      end

      subroutine convpdd(ianz,iconn,iresid,ncalf,iamino)
c same as conpdb, only before conversion to angstrom
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      common /athlp/ iatoms, mxnat
      common /cllab/ iclon,iclpnt(4)
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension iconn(mxcon+1,*),ianz(*),iresid(*),iamino(*)

      namatm = 0
      do i=1,iatoms
         iconn(1,i) = 0
         if (iresid(i).le.0.and.namatm.eq.0) namatm = i-1
      end do
      if (namatm.eq.0) namatm = iatoms

      do i=1,ncalf
         call getpdb(i,ipdb,ihpdb)
         call mkcon(ipdb,iamino(i),i,ihpdb,1,1)
      end do

c
c do connectivity between all except amino-amino
c
      nhet = iatoms - namatm
      do ii=1,nhet
         i = namatm + ii
         ir1 = iresid(i)
         do j=1,namatm + ii - 1
            ir2 = iresid(j)
            if (ir1.eq.ir2.and..not.(ir1.gt.0.and.ir2.gt.0)) then
               call connij(idum,i,j,0)
            endif
         end do
      end do

c
c make s-s connections
c
      do i=1,ncalf
         if (iamino(i).eq.4) then
            do j=1,i-1
               if (iamino(j).eq.4) then
                  ii = 0
                  jj = 0
                  do k=1,iatoms
                     if (ianz(k).eq.16) then
                        if (iresid(k).eq.i) ii = k
                        if (iresid(k).eq.j) jj = k
                     endif
                  end do
                  if (ii.ne.0.and.jj.ne.0) then
                     call connij(idcon,ii,jj,0)
                     if (idcon.eq.1) then
                        do l=1,iconn(1,ii)
                           ll = iabs(iconn(l+1,ii))
                        end do
                        do l=1,iconn(1,jj)
                           ll = iabs(iconn(l+1,jj))
                        end do
                     endif
                  endif
               endif
            end do
         endif
      end do

      return
      end

      subroutine consld(iconn,iresid)
c same as conpdb, only before conversion to angstrom, non residue
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      dimension iconn(mxcon+1,*),iresid(*)

      namatm = 0
      do i=1,iatoms
         if (iresid(i).lt.-3) then
            iconn(1,i) = 0
         endif
         if (iresid(i).lt.-3.and.namatm.eq.0) namatm = i-1
      end do
      if (namatm.eq.0) namatm = iatoms

c
c do connectivity between all except amino-amino
c
      nhet = iatoms - namatm
      do ii=1,nhet
         i = namatm + ii
         ir1 = iresid(i)
         do j=1,namatm + ii - 1
            ir2 = iresid(j)
            if (ir1.eq.ir2) then
               call connij(idcon,i,j,0)
            endif
         end do
      end do

      return
      end
c
      subroutine connid(idcon,i,j,idoconv,iconn,ianz,coo)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      parameter (mxcon=10)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      logical coni,conj
      dimension ctemp(3),iconn(mxcon+1,*),ianz(*),coo(3,*)

      toang = 0.52917706d0
      toang2 = toang*toang

      idcon = 0
      coni = .true.
      conj = .true.

      k = iconn(1,i)
      m = iconn(1,j)

      do l=1,k
         if (iconn(1+l,i).eq.j) coni = .false.
      end do

      do l=1,m
         if (iconn(1+l,j).eq.i) conj = .false.
      end do

      if (.not.coni.and..not.conj) return

      dmaxsq = vdwr(ianz(i)) + vdwr(ianz(j))
      dmaxsq = dmaxsq * dmaxsq

      do l=1,3
         ctemp(l) = coo(l,i) - coo(l,j)
      end do

      dijsq = ctemp(1)*ctemp(1) + ctemp(2)*ctemp(2) 
     &             + ctemp(3)*ctemp(3)

      if (idoconv.eq.1) dmaxsq = dmaxsq / toang2
      if (dijsq.lt.dmaxsq) then
          idocon = 1
          k = k + 1
          m = m + 1
          if (k.le.mxcon.and.coni) then
             iconn(1,i) = k
             iconn(k+1,i) = j
          else
             write(iun3,*)'more than mxconn connections found'
          endif
          if (m.le.mxcon.and.conj) then
             iconn(1,j) = m
             iconn(m+1,j) = i
          else
             write(iun3,*)'more than mxconn connections found'
          endif
      endif

      return
      end

      subroutine mkcon(ipdb,jres,icres,ihpdb,ihashy,idoconv)
      implicit double precision (a-h,o-z)
      dimension ipdb(*), ihpdb(*)

      if (jres.gt.23) then
         call mknbck(ipdb,ihpdb,jres,icres,ihashy,idoconv)
      else
         call mkback(ipdb,ihpdb,jres,icres,ihashy,idoconv)
      endif
          
      goto (10,10,30,40,50,60,70,80,90,100,110,120,130,140,150,
     &      160,170,180,190,200,100,140,150,210,220,230,240,250,
     &      210,220,220,230,230,230,230,230,230,230,250,250,250,
     &      250) jres
      return
c
c glycine, alanine
c
10    goto 1000
c
c serine
c
30    call conat(ipdb,5,31,0)
      if (ihashy.eq.1) call conath(ipdb,ihpdb,31,10)
      goto 1000
c
c cysteine
c
40    call conat(ipdb,5,37,0)
c this H is dubious
      if (ihashy.eq.1) call conath(ipdb,ihpdb,37,10)
      goto 1000
c
c threonine
c
50    call conat(ipdb,5,32,0)
      call conat(ipdb,5,8,0)
      if (ihashy.eq.1) then
          call conath(ipdb,ihpdb,32,10)
          call conath(ipdb,ihpdb,32,13)
          call conath(ipdb,ihpdb,8,16)
          call conath(ipdb,ihpdb,8,17)
          call conath(ipdb,ihpdb,8,18)
      endif
      goto 1000
c
c isoleucine
c
60    call conat(ipdb,5,7,0)
      call conat(ipdb,5,8,0)
      call conat(ipdb,7,10,0)
      if (ihashy.eq.1) then
          call conath(ipdb,ihpdb,7,13)
          call conath(ipdb,ihpdb,7,14)
          call conath(ipdb,ihpdb,8,16)
          call conath(ipdb,ihpdb,8,17)
          call conath(ipdb,ihpdb,8,18)
          call conath(ipdb,ihpdb,10,22)
          call conath(ipdb,ihpdb,10,23)
          call conath(ipdb,ihpdb,10,24)
      endif
      goto 1000
c
c valine
c
70    call conat(ipdb,5,7,0)
      call conat(ipdb,5,8,0)
      if (ihashy.eq.1) then
          call conath(ipdb,ihpdb,7,13)
          call conath(ipdb,ihpdb,7,14)
          call conath(ipdb,ihpdb,7,15)
          call conath(ipdb,ihpdb,8,16)
          call conath(ipdb,ihpdb,8,17)
          call conath(ipdb,ihpdb,8,18)
      endif
      goto 1000
c
c methionine
c
80    call conat(ipdb,5,6,0)
      call conat(ipdb,6,36,0)
      call conat(ipdb,36,12,0)
      if (ihashy.eq.1) then
          call conath(ipdb,ihpdb,6,10)
          call conath(ipdb,ihpdb,6,11)
c        try HG1 and HG2 too
          call conath(ipdb,ihpdb,6,13)
          call conath(ipdb,ihpdb,6,16)
          call conath(ipdb,ihpdb,12,28)
          call conath(ipdb,ihpdb,12,29)
          call conath(ipdb,ihpdb,12,30)
          call conath(ipdb,ihpdb,12,31)
          call conath(ipdb,ihpdb,12,34)
          call conath(ipdb,ihpdb,12,37)
      endif
      goto 1000
c
c aspartic acid
c
90    call conat(ipdb,5,6,0)
      call conat(ipdb,6,29,0)
      call conat(ipdb,6,30,0)
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,29,22)
         call conath(ipdb,ihpdb,30,25)
      endif
      goto 1000
c
c asparagine
c
100   call conat(ipdb,5,6,0)
      call conat(ipdb,6,29,0)
      call conat(ipdb,6,21,0)
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,21,25)
         call conath(ipdb,ihpdb,21,26)
      endif
      goto 1000
c
c leucine
c
110   call conat(ipdb,5,6,0)
      call conat(ipdb,6,10,0)
      call conat(ipdb,6,11,0)
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,6,10)
         call conath(ipdb,ihpdb,10,22)
         call conath(ipdb,ihpdb,10,23)
         call conath(ipdb,ihpdb,10,24)
         call conath(ipdb,ihpdb,11,25)
         call conath(ipdb,ihpdb,11,26)
         call conath(ipdb,ihpdb,11,27)
      endif
      goto 1000
c
c lysine
c
120   call conat(ipdb,5,6,0)
      call conat(ipdb,6,9,0)
      call conat(ipdb,9,12,0)
      call conat(ipdb,12,27,0)
      if (ihashy.eq.1) then

         call conath(ipdb,ihpdb,6,10)
         call conath(ipdb,ihpdb,6,11)
c        try HG1 and HG2 too
         call conath(ipdb,ihpdb,6,13)
         call conath(ipdb,ihpdb,6,16)

         call conath(ipdb,ihpdb,9,19)
         call conath(ipdb,ihpdb,9,20)
c        try HD1 and HD2 too
         call conath(ipdb,ihpdb,9,22)
         call conath(ipdb,ihpdb,9,25)

         call conath(ipdb,ihpdb,12,28)
         call conath(ipdb,ihpdb,12,29)
c        try HE1 and HE2 , HE3 too
         call conath(ipdb,ihpdb,12,31)
         call conath(ipdb,ihpdb,12,34)
         call conath(ipdb,ihpdb,12,37)

         call conath(ipdb,ihpdb,27,40)
         call conath(ipdb,ihpdb,27,41)
         call conath(ipdb,ihpdb,27,42)
         call conath(ipdb,ihpdb,27,43)
         call conath(ipdb,ihpdb,27,46)
         call conath(ipdb,ihpdb,27,49)
      endif
      goto 1000
c
c glutamic acid
c
130   call conat(ipdb,5,6,0)
      call conat(ipdb,6,9,0)
      call conat(ipdb,9,34,0)
      call conat(ipdb,9,35,0)
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,6,10)
         call conath(ipdb,ihpdb,6,11)
c        try HG1 and HG2 too
         call conath(ipdb,ihpdb,6,13)
         call conath(ipdb,ihpdb,6,16)
         call conath(ipdb,ihpdb,34,31)
         call conath(ipdb,ihpdb,35,34)
      endif
      goto 1000
c
c glutamine
c
140   call conat(ipdb,5,6,0)
      call conat(ipdb,6,9,0)
      call conat(ipdb,9,34,0)
      call conat(ipdb,9,24,0)
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,6,10)
         call conath(ipdb,ihpdb,6,11)
c        try HG1 and HG2 too
         call conath(ipdb,ihpdb,6,13)
         call conath(ipdb,ihpdb,6,16)
         call conath(ipdb,ihpdb,24,34)
         call conath(ipdb,ihpdb,24,35)
      endif
      goto 1000
c
c proline
c 
150   call conat(ipdb,5,6,0)
      call conat(ipdb,6,9,0)
c     bit tricky here
      call conat(ipdb,9,1,1)
      if (ipdb(28).ne.0) call conat(ipdb,6,28,0)
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,6,10)
         call conath(ipdb,ihpdb,6,11)
c        try HG1 and HG2 too
         call conath(ipdb,ihpdb,6,13)
         call conath(ipdb,ihpdb,6,16)
         call conath(ipdb,ihpdb,9,19)
         call conath(ipdb,ihpdb,9,20)
c        try HD1 and HD2 too
         call conath(ipdb,ihpdb,9,22)
         call conath(ipdb,ihpdb,9,25)
      endif
      goto 1000
c
c arginine
c
160   call conat(ipdb,5,6,0)
      call conat(ipdb,6,9,0)
      call conat(ipdb,9,22,0)
      call conat(ipdb,22,17,0)
      call conat(ipdb,17,25,0)
      call conat(ipdb,17,26,0)
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,6,10)
         call conath(ipdb,ihpdb,6,11)
         call conath(ipdb,ihpdb,6,13)
         call conath(ipdb,ihpdb,6,16)
         call conath(ipdb,ihpdb,9,19)
         call conath(ipdb,ihpdb,9,20)
         call conath(ipdb,ihpdb,9,22)
         call conath(ipdb,ihpdb,9,25)
         call conath(ipdb,ihpdb,22,28)
         call conath(ipdb,ihpdb,25,55)
         call conath(ipdb,ihpdb,25,56)
         call conath(ipdb,ihpdb,26,58)
         call conath(ipdb,ihpdb,26,59)
      endif
      goto 1000
c
c histidine
c
170   call conat(ipdb,5,6,0)
      call conat(ipdb,6,11,0)
      call conat(ipdb,6,20,0)
      call conat(ipdb,11,24,0)
      call conat(ipdb,20,13,0)
c tricky
      call conat(ipdb,13,24,1)
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,20,22)
         call conath(ipdb,ihpdb,11,25)
         call conath(ipdb,ihpdb,13,31)
         call conath(ipdb,ihpdb,24,34)
      endif
      goto 1000
c
c phenylalanine
c
180   call conat(ipdb,5,6,0)
      call conat(ipdb,6,10,0)
      call conat(ipdb,6,11,0)
      call conat(ipdb,10,13,0)
      call conat(ipdb,11,14,0)
      call conat(ipdb,13,17,0)
      call conat(ipdb,14,17,1)
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,10,22)
         call conath(ipdb,ihpdb,11,25)
         call conath(ipdb,ihpdb,13,31)
         call conath(ipdb,ihpdb,14,34)
         call conath(ipdb,ihpdb,17,40)
      endif
      goto 1000
c
c tyrosine
c
190   call conat(ipdb,5,6,0)
      call conat(ipdb,6,10,0)
      call conat(ipdb,6,11,0)
      call conat(ipdb,10,13,0)
      call conat(ipdb,11,14,0)
      call conat(ipdb,13,17,0)
      call conat(ipdb,14,17,1)
      call conat(ipdb,17,33,0)
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,10,22)
         call conath(ipdb,ihpdb,11,25)
         call conath(ipdb,ihpdb,13,31)
         call conath(ipdb,ihpdb,14,34)
         call conath(ipdb,ihpdb,33,52)
      endif
      goto 1000
c
c tryptophan
c
200   call conat(ipdb,5,6,0)
      call conat(ipdb,6,10,0)
      call conat(ipdb,6,11,0)
      call conat(ipdb,10,23,0)
      call conat(ipdb,23,14,0)
      call conat(ipdb,14,11,1)
      call conat(ipdb,11,15,0)
      call conat(ipdb,14,18,0)
      call conat(ipdb,15,19,0)
      call conat(ipdb,18,16,0)
      call conat(ipdb,16,19,1)
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,10,22)
         call conath(ipdb,ihpdb,23,31)
         call conath(ipdb,ihpdb,15,37)
         call conath(ipdb,ihpdb,18,46)
         call conath(ipdb,ihpdb,19,49)
         call conath(ipdb,ihpdb,16,58)
      endif
      goto 1000
c
c adenosine
c
210   call conat(ipdb,54,66,0)
      call conat(ipdb,66,59,0)
      call conat(ipdb,66,56,0)
      call conat(ipdb,59,65,0)
      call conat(ipdb,65,57,0)
      call conat(ipdb,57,58,0)
      call conat(ipdb,57,56,0)
      call conat(ipdb,58,64,0)
      call conat(ipdb,58,60,0)
      call conat(ipdb,60,55,0)
      call conat(ipdb,55,62,0)
      call conat(ipdb,56,62,1)
      call conat(ipdb,56,66,1)
      if (jres.eq.29) then
          call conat(ipdb,60,71,0)
          call conat(ipdb,60,79,0)
      endif
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,55,82)
         call conath(ipdb,ihpdb,59,112)
         call conath(ipdb,ihpdb,64,106)
         call flagh(ihpdb,106)
         call conath(ipdb,ihpdb,64,109)
         call flagh(ihpdb,109)
         call conath(ipdb,ihpdb,64,136)
         call flagh(ihpdb,136)
         call conath(ipdb,ihpdb,64,137)
         call flagh(ihpdb,137)
         call conath(ipdb,ihpdb,64,103)
         call flagh(ihpdb,103)
         call conath(ipdb,ihpdb,64,104)
         call flagh(ihpdb,104)
         if (jres.eq.29) then
             call conath(ipdb,ihpdb,79,148)
             call conath(ipdb,ihpdb,79,149)
             call conath(ipdb,ihpdb,79,150)
             call conath(ipdb,ihpdb,58,103)
         endif
      endif
      goto 1000
c
c cytidine
c
220   call conat(ipdb,54,60,0)
      call conat(ipdb,60,55,0)
      call conat(ipdb,55,67,0)
      call conat(ipdb,55,62,0)
      call conat(ipdb,62,56,0)
      call conat(ipdb,56,63,0)
      call conat(ipdb,56,57,0)
      call conat(ipdb,57,58,0)
      call conat(ipdb,60,58,1)
      if (jres.eq.30) then
         call conat(ipdb,57,72,0)
         call conat(ipdb,57,83,0)
      endif
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,63,94)
         call flagh(ihpdb,94)
         call conath(ipdb,ihpdb,63,97)
         call flagh(ihpdb,97)
         call conath(ipdb,ihpdb,63,130)
         call flagh(ihpdb,130)
         call conath(ipdb,ihpdb,63,131)
         call flagh(ihpdb,131)
         call conath(ipdb,ihpdb,58,103)
         if (jres.eq.30) then
             call conath(ipdb,ihpdb,83,160)
             call conath(ipdb,ihpdb,83,161)
             call conath(ipdb,ihpdb,83,162)
         else
             call conath(ipdb,ihpdb,57,100)
         endif
      endif
      goto 1000
c
c guanosine
c
230   call conat(ipdb,54,66,0)
      call conat(ipdb,66,59,0)
      call conat(ipdb,66,56,0)
      call conat(ipdb,59,65,0)
      call conat(ipdb,65,57,0)
      call conat(ipdb,57,58,0)
      call conat(ipdb,57,56,0)
      call conat(ipdb,58,69,0)
      call conat(ipdb,58,60,0)
      call conat(ipdb,60,55,0)
      call conat(ipdb,55,62,0)
      if (jres.ne.38) call conat(ipdb,55,61,0)
      call conat(ipdb,56,62,1)
      call conat(ipdb,56,66,1)
      if (jres.eq.32) then
          call conat(ipdb,60,71,0)
          call conat(ipdb,60,79,0)
      endif
      if (jres.eq.33.or.jres.eq.34) then
          call conat(ipdb,61,73,0)
          call conat(ipdb,61,80,0)
      endif
      if (jres.eq.34) then
          call conat(ipdb,61,74,0)
          call conat(ipdb,61,79,0)
      endif
      if (jres.eq.35) then
          call conat(ipdb,65,75,0)
          call conat(ipdb,65,85,0)
      endif
      if (jres.eq.37) then
          call conat(ipdb,61,90,0)
          call conat(ipdb,62,88,0)
          call conat(ipdb,90,89,0)
          call conat(ipdb,90,91,0)
          call conat(ipdb,91,60,1)
          call conat(ipdb,91,92,0)
          call conat(ipdb,92,93,0)
          call conat(ipdb,93,94,0)
          call conat(ipdb,94,95,0)
          call conat(ipdb,95,96,0)
          call conat(ipdb,95,97,0)
          call conat(ipdb,97,98,0)
          call conat(ipdb,94,99,0)
          call conat(ipdb,99,100,0)
          call conat(ipdb,100,101,0)
          call conat(ipdb,100,102,0)
          call conat(ipdb,102,103,0)
      endif

      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,59,112)
         if (jres.ne.32.and.jres.ne.37) then
            call conath(ipdb,ihpdb,60,79)
            call flagh(ihpdb,79)
            call conath(ipdb,ihpdb,60,121)
            call flagh(ihpdb,121)
         endif
         if (jres.ne.33.and.jres.ne.34.and.jres.ne.37
     &       .and.jres.ne.38) then
            call conath(ipdb,ihpdb,61,85)
            call flagh(ihpdb,85)
            call conath(ipdb,ihpdb,61,125)
            call flagh(ihpdb,125)
            call conath(ipdb,ihpdb,61,82)
            call flagh(ihpdb,82)
            call conath(ipdb,ihpdb,61,83)
            call flagh(ihpdb,83)
         endif
         if (jres.ne.34.and.jres.ne.37.and.jres.ne.38) then
            call conath(ipdb,ihpdb,61,88)
            call flagh(ihpdb,88)
            call conath(ipdb,ihpdb,61,124)
            call flagh(ihpdb,124)
         endif
         if (jres.eq.33.or.jres.eq.34) then
            call conath(ipdb,ihpdb,80,151)
            call conath(ipdb,ihpdb,80,152)
            call conath(ipdb,ihpdb,80,153)
         endif
         if (jres.eq.31.or.jres.eq.34) then
            call conath(ipdb,ihpdb,79,148)
            call conath(ipdb,ihpdb,79,149)
            call conath(ipdb,ihpdb,79,150)
         endif
         if (jres.eq.35) then
            call conath(ipdb,ihpdb,85,166)
            call conath(ipdb,ihpdb,85,167)
            call conath(ipdb,ihpdb,85,168)
            call conath(ipdb,ihpdb,59,113)
         endif
         if (jres.eq.37) then
            call conath(ipdb,ihpdb,88,91)
            call conath(ipdb,ihpdb,88,92)
            call conath(ipdb,ihpdb,88,93)
            call conath(ipdb,ihpdb,89,175)
            call conath(ipdb,ihpdb,89,176)
            call conath(ipdb,ihpdb,89,177)
            call conath(ipdb,ihpdb,92,178)
            call conath(ipdb,ihpdb,92,179)
            call conath(ipdb,ihpdb,93,181)
            call conath(ipdb,ihpdb,93,182)
            call conath(ipdb,ihpdb,94,184)
            call conath(ipdb,ihpdb,98,187)
            call conath(ipdb,ihpdb,98,188)
            call conath(ipdb,ihpdb,98,189)
            call conath(ipdb,ihpdb,99,124)
            call conath(ipdb,ihpdb,103,190)
            call conath(ipdb,ihpdb,103,191)
            call conath(ipdb,ihpdb,103,192)
         endif
      endif
      goto 1000
c
c thymidine
c
240   call conat(ipdb,54,60,0)
      call conat(ipdb,60,55,0)
      call conat(ipdb,55,67,0)
      call conat(ipdb,55,62,0)
      call conat(ipdb,62,56,0)
      call conat(ipdb,56,68,0)
      call conat(ipdb,56,57,0)
      call conat(ipdb,57,58,0)
      call conat(ipdb,60,58,1)
      call conat(ipdb,57,70,1)
      call conat(ipdb,57,83,1)
      call conat(ipdb,57,75,1)
      call conat(ipdb,57,85,1)
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,58,103)
         call conath(ipdb,ihpdb,62,91)
         call flagh(ihpdb,91)
         call conath(ipdb,ihpdb,62,127)
         call flagh(ihpdb,127)
         call conath(ipdb,ihpdb,83,160)
         call conath(ipdb,ihpdb,83,161)
         call conath(ipdb,ihpdb,83,162)
         call conath(ipdb,ihpdb,70,160)
         call conath(ipdb,ihpdb,70,161)
         call conath(ipdb,ihpdb,70,162)
         call conath(ipdb,ihpdb,75,166)
         call conath(ipdb,ihpdb,75,167)
         call conath(ipdb,ihpdb,75,168)
         call conath(ipdb,ihpdb,85,166)
         call conath(ipdb,ihpdb,85,167)
         call conath(ipdb,ihpdb,85,168)
      endif
      goto 1000
c
c uridine
c
250   if (jres.eq.42) then
         call conat(ipdb,54,57,0)
      else
         call conat(ipdb,54,60,0)
      endif
      call conat(ipdb,60,55,0)
      call conat(ipdb,55,67,0)
      call conat(ipdb,55,62,0)
      call conat(ipdb,62,56,0)
      call conat(ipdb,56,68,0)
      call conat(ipdb,56,57,1)
      call conat(ipdb,57,58,0)
      call conat(ipdb,60,58,1)
      if (jres.eq.41) then
         call conat(ipdb,57,72,0)
         call conat(ipdb,57,83,0)
      endif
      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,58,103)
         if (jres.ne.41.and.jres.ne.42) call conath(ipdb,ihpdb,57,100)
         if (jres.eq.40) then
            call conath(ipdb,ihpdb,58,104)
            call conath(ipdb,ihpdb,57,101)
         endif
         call conath(ipdb,ihpdb,62,91)
         call flagh(ihpdb,91)
         call conath(ipdb,ihpdb,62,127)
         call flagh(ihpdb,127)
         if (jres.eq.41) then
            call conath(ipdb,ihpdb,83,160)
            call conath(ipdb,ihpdb,83,161)
            call conath(ipdb,ihpdb,83,162)
         endif
         if (jres.eq.42) call conath(ipdb,ihpdb,60,121)
      endif

1000  continue

      return
      end

      subroutine pdbtyd(ipdb,ihpdb,jres,ihashy,ipdbt)
      implicit double precision (a-h,o-z)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      integer*2 ipdbt
      dimension ipdb(mxsym), ihpdb(mxhsym*3), ipdbt(*)

      do i=1,mxsym
         it = i
         if (i.eq.70) it = 83
         if (i.eq.71) it = 79
         if (i.eq.72) it = 83
         if (i.eq.73) it = 80
         if (i.eq.74) it = 79
         if (i.eq.75) it = 85
         if (jres.gt.23.and.i.eq.38) it = 76
         if (ipdb(i).gt.0) ipdbt(ipdb(i)) = it
      end do

      if (ihashy.eq.1) then
         do i=1,mxhsym*3
            it = i
            if (i.eq.85) it = 124
            if (i.eq.88) it = 125
            if (i.eq.94) it = 130
            if (i.eq.97) it = 131
            if (i.eq.106) it = 136
            if (i.eq.109) it = 137
            if (ihpdb(i).gt.0) ipdbt(ihpdb(i)) = it
         end do
      endif

      return
      end

      subroutine typeid(ipdb,jres,ihpdb,ihashy,
     &                  ianz,iconn,ityp)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxrss=20)
      parameter (mxatt=10)
      parameter (mxath=8)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      logical nterm
      integer*2 ityp
      common /types/ iff
      common /chmtyp/ ncc(mxrss),icc(2,mxatt,mxrss),
     &                nhh(mxrss),ihh(2,mxath,mxrss)
      dimension ipdb(*), ihpdb(*),ianz(*),iconn(mxcon+1,*),ityp(*)

      iff = 2

c     currently only works for aminoacids, not nucleic acids 
c     (mxrss=20)

c     Backbone, no N-Term, C-Term yet
c               (Pro,Gly N-Term) CYS, sulfur bridged

      do i=1,mxsym
         call settyp(ityp,ipdb(i),0)
      end do
      do i=1,mxhsym*3
         call settyp(ityp,ihpdb(i),0)
      end do

      if (jres.le.0.or.jres.gt.mxrss) return

      if (jres.eq.17) then
c ihis
c 
c  1  HIS+ (HIP)
c  2  HISD (HID)
c  3  HISE (HIE)
c
         ihis = 1
         if (ihashy.eq.1) then
            if (ihpdb(34).eq.0) ihis = 2
            if (ihpdb(22).eq.0) ihis = 3
         endif
      endif

      call settyp(ityp,ipdb(1),63)
      call settyp(ityp,ipdb(2),23)
      call settyp(ityp,ipdb(3),20)
      call settyp(ityp,ipdb(4),74)

      if (ipdb(38).ne.0) then
          ityp(ipdb(38)) = 79
          call settyp(ityp,ipdb(4),79)
          call settyp(ityp,ipdb(3),22)
      endif

      do i=1,ncc(jres)
         call settyp(ityp,ipdb(icc(1,i,jres)),icc(2,i,jres))
      end do

      if (jres.eq.4) then
         do i=1,iconn(1,ipdb(37))
            if (ianz(abs(iconn(1+i,ipdb(37)))).eq.16) 
     &          call settyp(ityp,ipdb(37),82)
         end do
      endif

      if (jres.eq.17) then
         if (ihis.eq.2) then
             call settyp(ityp,ipdb(6),43)
             call settyp(ityp,ipdb(11),44)
             call settyp(ityp,ipdb(13),46)
             call settyp(ityp,ipdb(20),69)
             call settyp(ityp,ipdb(24),68)
         endif
         if (ihis.eq.3) then
             call settyp(ityp,ipdb(6),44)
             call settyp(ityp,ipdb(11),43)
             call settyp(ityp,ipdb(13),46)
             call settyp(ityp,ipdb(20),68)
             call settyp(ityp,ipdb(24),69)
         endif
      endif
      
      if (ihashy.eq.1) then

         n = 0
         do i=1,3
            if (ihpdb(i).ne.0) n = n + 1
         end do
         nterm = .false.
         if (n.eq.3.or.(jres.eq.15.and.n.eq.2)) nterm = .true.

         if (nterm) then
c H
            nt = 6
            if (jres.eq.15) nt = 7
            do i=1,n
               call settyp(ityp,ihpdb(i),nt)
            end do

c N
            nt = 65
            if (jres.eq.15.and.n.eq.2) nt = 67
            call settyp(ityp,ipdb(1),nt)

c CA
            nt = 23
            if (jres.eq.1) nt = 29
            if (jres.eq.15) nt = 33
            call settyp(ityp,ipdb(2),nt)

c CD
            if (jres.eq.15) call settyp(ityp,ipdb(9),34)

         else
            call settyp(ityp,ihpdb(1),3)
         endif
c HA
         nt = 4
         if (nterm.and.jres.ne.1) nt = 5
         call settyp(ityp,ihpdb(4),nt)
         call settyp(ityp,ihpdb(5),nt)

c HB, HG
         do i=1,3
            call settyp(ityp,ihpdb(6+i),1)
            call settyp(ityp,ihpdb(9+i),1)
         end do

         do i=1,nhh(jres)
            call settyp(ityp,ihpdb(ihh(1,i,jres)),ihh(2,i,jres))
         end do

c check for lysine alternative HZ1,HZ2,HZ3 instead of 1HZ,2HZ,3HZ

         if (jres.eq.12) then
            call settyp(ityp,ihpdb(43),6)
            call settyp(ityp,ihpdb(46),6)
            call settyp(ityp,ihpdb(49),6)
         endif

         if (jres.eq.17) then
            if (ihis.eq.2) then
                call settyp(ityp,ihpdb(22),11)
                call settyp(ityp,ihpdb(25),14)
                call settyp(ityp,ihpdb(31),12)
            endif
            if (ihis.eq.3) then
                call settyp(ityp,ihpdb(25),15)
                call settyp(ityp,ihpdb(31),12)
                call settyp(ityp,ihpdb(34),11)
            endif
         endif
      endif

      return
      end

      subroutine chkrna(ncalf,iamino)
      implicit double precision (a-h,o-z)
      dimension iamino(*)

      modres = 0
      do i=1,ncalf
         if (iamino(i).eq.32) modres = 1
         if (iamino(i).eq.38) modres = 1
         if (iamino(i).eq.39) modres = 1
      end do

      if (modres.eq.1) then
         print*,' '
         print*,'The following Modified nucleosides '
         print*,'are not supported by Ambfor !:'
         print*,' '
         print*,'(1MG,I,+U)'
         print*,' '
      endif

      return
      end

      subroutine typamd(ipdb,jres,ihpdb,ihashy,
     &                  ianz,iconn,ityp)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxres=42)
      parameter (mxrss=20)
      parameter (mxrsa=23)
      parameter (mxata=9)
      parameter (mxatha=9)
      parameter (mxrso=24)
      parameter (mxato=10)
      parameter (mxatho=11)
      parameter (mxrsn=19)
      parameter (mxnucl=39)
      parameter (mxhnuc=27)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxambc=49)
      logical nterm,cterm
      integer*2 ityp
      common /types/ iff
      common /athlp/ iatoms, mxnat
      integer amboff,ntca,ctca
      common /ambtyp/ amboff(mxrss*3),ntca(mxrss),ctca(mxrss),
     &                ncca(mxrsa),icca(2,mxata,mxrsa),
     &                nhha(mxrsa),ihha(2,mxatha,mxrsa),
     &                ncco(mxrso),icco(2,mxato,mxrso),
     &                nhho(mxrso),ihho(2,mxatho,mxrso),
     &                nnuc(mxrsn),nhnuc(mxrsn),irna(mxres-23),
     &                inuc(2,mxnucl,mxrsn),ihnuc(2,mxhnuc,mxrsn)
      integer ambvdt
      common /fcharg/ ambchg(mxamb),ambvw1(mxambc),ambvw2(mxambc),
     &                gfvdw(2,mxgff),ambvdt(mxamb),cysneg(9)
      logical ismet

      dimension ipdb(*), ihpdb(*),ianz(*),iconn(mxcon+1,*),ityp(*),
     &          icn(4)

      iff = 3

c     currently only works for aminoacids, not nucleic acids 
c     (mxrss=20) (HIP,HID,HIE and sulfer bridge covered)
c                (ORN, mALA and pyroglut not covered)


      if (jres.le.0.or.jres.gt.mxres) return
      if (jres.gt.mxrss) goto 100

      do i=1,mxsym
         call settyp(ityp,ipdb(i),0)
         call setchg(ipdb(i),0)
      end do
      do i=1,mxhsym*3
         call settyp(ityp,ihpdb(i),0)
         call setchg(ihpdb(i),0)
      end do

      nterm = .false.
      n = 0
      do i=1,3
         if (ihpdb(i).ne.0) n = n + 1
      end do
      if (n.eq.0) then
         if (ihpdb(79).ne.0) n = n + 1
         if (ihpdb(82).ne.0) n = n + 1
         if (ihpdb(91).ne.0) n = n + 1
      endif
      if (n.eq.3.or.(jres.eq.15.and.n.eq.2)) nterm = .true.
      cterm = .false.
      if (ipdb(38).ne.0) cterm = .true.

      if (jres.eq.17) then
c ihis
c 
c  1  HIS+ (HIP)
c  2  HISD (HID)
c  3  HISE (HIE)
c
         ihis = 1
         if (ihashy.eq.1) then
            if (ihpdb(34).eq.0) ihis = 2
            if (ihpdb(22).eq.0) ihis = 3
         endif
      endif

      if (jres.eq.4) then
c iss = 0 (cysteine -SH), iss = 1 (cystine -SS-)
         iss = 0
         if (ipdb(37).ne.0) then
            do i=1,iconn(1,ipdb(37))
               if (ianz(abs(iconn(1+i,ipdb(37)))).eq.16) 
     &             iss = 1
            end do
         endif
      endif

      if (nterm) then
          ioff = amboff(jres+20)
      else if(cterm) then
          ioff = amboff(jres+40)
      else
        ioff = amboff(jres)
      endif

      if (jres.eq.17) then
         if (ihis.eq.2) ioff = 178
         if (ihis.eq.3) ioff = 194
      endif

      if (jres.eq.4.and.iss.eq.1) ioff = 87
      if (jres.eq.9.and.(ihpdb(22).ne.0.or.ihpdb(25).ne.0)) ioff = 660
      if (jres.eq.13.and.(ihpdb(31).ne.0.or.ihpdb(34).ne.0)) ioff = 672

      ihtot = 0
      do i=0,2
         if (ihpdb(40+i).ne.0) ihtot = ihtot + 1
      end do

      if (jres.eq.12.and.ihtot.eq.2) ioff = 686

      if (nterm) then
         if (jres.eq.17 .and. ihis.eq.2) ioff = 436
         if (jres.eq.17 .and. ihis.eq.3) ioff = 442
         if (jres.eq.4 .and. iss.eq.1) ioff = 398
      endif
      if (cterm) then
         if (jres.eq.17 .and. ihis.eq.2) ioff = 584
         if (jres.eq.17 .and. ihis.eq.3) ioff = 590
         if (jres.eq.4 .and. iss.eq.1) ioff = 549
      endif

      call settyp(ityp,ipdb(1),ioff)
      call settyp(ityp,ipdb(2),ioff+1)
      call settyp(ityp,ipdb(3),ioff+2)

      if (jres.eq.15.and..not.nterm) then
         call settyp(ityp,ipdb(4),ioff+3)
         if (cterm) call settyp(ityp,ipdb(38),ioff+3)
      else
         call settyp(ityp,ipdb(4),ioff+4)
         if (cterm) call settyp(ityp,ipdb(38),ioff+4)
      endif

      call settyp(ityp,ihpdb(1),ioff+3)
      call settyp(ityp,ihpdb(2),ioff+3)
      call settyp(ityp,ihpdb(3),ioff+3)
      call settyp(ityp,ihpdb(79),ioff+3)
      call settyp(ityp,ihpdb(82),ioff+3)
      call settyp(ityp,ihpdb(91),ioff+3)

      if (jres.eq.15) then
         if (nterm) then
            call settyp(ityp,ihpdb(4),409)
         else
            call settyp(ityp,ihpdb(4),ioff+4)
            call settyp(ityp,ihpdb(5),ioff+4)
         endif
      else
         call settyp(ityp,ihpdb(4),ioff+5)
         call settyp(ityp,ihpdb(5),ioff+5)
      endif

      if (nterm.or.cterm) then
         ioff = amboff(jres)
         if (jres.eq.17) then
            if (ihis.eq.2) ioff = 178
            if (ihis.eq.3) ioff = 194
         endif
         if (jres.eq.4.and.iss.eq.1) ioff = 87
      endif

      if (jres.eq.15) then
         call settyp(ityp,ipdb(5),ioff+5)
         call settyp(ityp,ihpdb(7),ioff+6)
         call settyp(ityp,ihpdb(8),ioff+6)
         call settyp(ityp,ihpdb(9),ioff+6)
      else
         call settyp(ityp,ipdb(5),ioff+6)
         call settyp(ityp,ihpdb(7),ioff+7)
         call settyp(ityp,ihpdb(8),ioff+7)
         call settyp(ityp,ihpdb(9),ioff+7)
      endif

      jtmp = jres
      if (jres.eq.17) then
         if (ihis.eq.2) jtmp = 22
         if (ihis.eq.3) jtmp = 23
      endif
      if (jres.eq.4.and.iss.eq.1) jtmp = 21
    
      do i=1,ncca(jtmp)
         call settyp(ityp,ipdb(icca(1,i,jtmp)),icca(2,i,jtmp))
      end do

      do i=1,nhha(jtmp)
         call settyp(ityp,ihpdb(ihha(1,i,jtmp)),ihha(2,i,jtmp))
      end do

      if (jres.eq.8) then
         call settyp(ityp,ihpdb(13),267)
         call settyp(ityp,ihpdb(16),267)
         call settyp(ityp,ihpdb(31),270)
         call settyp(ityp,ihpdb(34),270)
         call settyp(ityp,ihpdb(37),270)
      endif
      if (jres.eq.13) then
         call settyp(ityp,ihpdb(13),241)
         call settyp(ityp,ihpdb(16),241)
      endif
      if (jres.eq.14) then
         call settyp(ityp,ihpdb(13),253)
         call settyp(ityp,ihpdb(16),253)
      endif
      if (jres.eq.16) then
         call settyp(ityp,ihpdb(13),296)
         call settyp(ityp,ihpdb(16),296)
         call settyp(ityp,ihpdb(22),298)
         call settyp(ityp,ihpdb(25),298)
      endif

      if (jres.eq.15) then
         if (nterm) then
            call settyp(ityp,ipdb(9),410)
            call settyp(ityp,ihpdb(19),411)
            call settyp(ityp,ihpdb(20),411)
         endif
CNF  Nothing different between PRO and other residues for the C-Term case
CNF         if (cterm) then
CNF            call settyp(ityp,ipdb(9),412)
CNF            call settyp(ityp,ihpdb(19),413)
CNF            call settyp(ityp,ihpdb(20),413)
CNF         endif
      endif

c check for lysine alternative HZ1,HZ2,HZ3 instead of 1HZ,2HZ,3HZ

      if (jres.eq.12) then
         call settyp(ityp,ihpdb(13),280)
         call settyp(ityp,ihpdb(16),280)
         call settyp(ityp,ihpdb(22),282)
         call settyp(ityp,ihpdb(25),282)
         call settyp(ityp,ihpdb(34),284)
         call settyp(ityp,ihpdb(37),284)
         call settyp(ityp,ihpdb(43),286)
         call settyp(ityp,ihpdb(46),286)
         call settyp(ityp,ihpdb(49),286)
      endif

c check asph, RECTIFICATION: asph, gluh neutral lys, neg cys
c are not in the tinker amber set.
c we could add our own, but it has not been done yet
c Original AMBER itself does not work with so many atom types.
c Tinker creates the extra atom types to incrporate the charges
c in there as well. Charges deviate slightly between
c cterm,nterm and regular residues
c see ~schaft/compile/linux/molden4.6/
c     forf/amber9.ffparms/dat/leap/prep/all_amino94.in
c for the charges
c we should add our own extra types for these residues

      if (jres.eq.9) then
         if (ihpdb(22).ne.0) then
            call settyp(ityp,ipdb(6),668)
            call settyp(ityp,ipdb(29),670)
            call settyp(ityp,ipdb(30),669)
            call settyp(ityp,ihpdb(22),671)
         endif
         if (ihpdb(25).ne.0) then
            call settyp(ityp,ipdb(6),668)
            call settyp(ityp,ipdb(29),669)
            call settyp(ityp,ipdb(30),670)
            call settyp(ityp,ihpdb(25),671)
         endif
      endif

      if (jres.eq.13) then
         if (ihpdb(31).ne.0) then
            call settyp(ityp,ipdb(6),680)
            call settyp(ityp,ipdb(9),682)
            call settyp(ityp,ipdb(34),684)
            call settyp(ityp,ipdb(35),683)
            call settyp(ityp,ihpdb(10),681)
            call settyp(ityp,ihpdb(11),681)
            call settyp(ityp,ihpdb(31),685)
         endif
         if (ihpdb(34).ne.0) then
            call settyp(ityp,ipdb(6),680)
            call settyp(ityp,ipdb(9),682)
            call settyp(ityp,ipdb(34),683)
            call settyp(ityp,ipdb(35),684)
            call settyp(ityp,ihpdb(10),681)
            call settyp(ityp,ihpdb(11),681)
            call settyp(ityp,ihpdb(34),685)
         endif
      endif

      if (jres.eq.12.and.ihtot.eq.2) then
            call settyp(ityp,ipdb(6),694)
            call settyp(ityp,ipdb(9),696)
            call settyp(ityp,ipdb(12),698)
            call settyp(ityp,ipdb(27),700)
            call settyp(ityp,ihpdb(10),695)
            call settyp(ityp,ihpdb(11),695)
            call settyp(ityp,ihpdb(19),697)
            call settyp(ityp,ihpdb(20),697)
            call settyp(ityp,ihpdb(28),699)
            call settyp(ityp,ihpdb(29),699)
            call settyp(ityp,ihpdb(40),701)
            call settyp(ityp,ihpdb(41),701)
            call settyp(ityp,ihpdb(42),701)
            call settyp(ityp,ihpdb(43),701)
            call settyp(ityp,ihpdb(46),701)
            call settyp(ityp,ihpdb(49),701)
      endif

      do i=1,mxsym
         call setchg(ipdb(i),1)
      end do
      do i=1,mxhsym*3
         call setchg(ihpdb(i),1)
      end do

      if (jres.eq.4) then

c ic = 0 (negative cysteine), apply different charges

         if (ipdb(37).ne.0) then
            ic = 0
            do i=1,iconn(1,ipdb(37))
               ii = iconn(1+i,ipdb(37))
               if (ii.gt.0.and.ii.ne.ipdb(5)) then
                   if (.not.ismet(ianz(ii))) ic = ic + 1
               endif
            end do

            if (ic.eq.0) then

               do i=1,mxsym
                  if (ipdb(i).ne.0) then
                     it = ityp(ipdb(i)) 
                     if (it.ge.77.and.it.le.85) then
                        it = 76 - it
                        call setchg(ipdb(i),it)
                     endif
                  endif
               end do

               do i=1,mxhsym*3
                  if (ihpdb(i).ne.0) then
                     it = ityp(ihpdb(i)) 
                     if (it.ge.77.and.it.le.85) then
                        it = 76 - it
                        call setchg(ihpdb(i),it)
                     endif
                  endif
               end do

            endif
         endif
      endif

      return

100   continue

c     NUCLEOTIDES

      isrna = 1
      jtmp = irna(jres - 23)

      if (ipdb(53).eq.0.and.jres.ne.27) then
         isrna = 0
         jtmp = jtmp + 4
      endif

      do i=1,nnuc(jtmp)
         call settyp(ityp,ipdb(inuc(1,i,jtmp)),inuc(2,i,jtmp))
      end do

      do i=1,nhnuc(jtmp)
         call settyp(ityp,ihpdb(ihnuc(1,i,jtmp)),ihnuc(2,i,jtmp))
      end do

      if (ihpdb(115).ne.0) then
c        3'-Hydroxyl
         if (isrna.eq.1) then
            call settyp(ityp,ipdb(51),1237)
            call settyp(ityp,ihpdb(115),1238)
         else
            call settyp(ityp,ipdb(51),1249)
            call settyp(ityp,ihpdb(115),1250)
         endif
      endif

      if (ihpdb(118).ne.0) then
c        5'-Hydroxyl
         if (isrna.eq.1) then
            call settyp(ityp,ipdb(46),1232)
            call settyp(ityp,ihpdb(118),1233)
         else
            call settyp(ityp,ipdb(46),1244)
            call settyp(ityp,ihpdb(118),1245)
         endif
      endif

c     check P connected to O3* (ipdb(51)), if nconn O = 4 => 3'-Phosphate

      call chkpo4(ipdb(51),isrna,1239,1240,1241,1251,1252,1253,
     &                  ianz,iconn,ityp)

c     check P connected to O5* (ipdb(46)), if nconn O = 4 => 5'-Phosphate

      call chkpo4(ipdb(46),isrna,1234,1235,1236,1246,1247,1248,
     &                  ianz,iconn,ityp)

      do i=1,mxsym
         call setchg(ipdb(i),1)
      end do

      do i=1,mxhsym*3
         call setchg(ihpdb(i),1)
      end do

      return
      end

      subroutine chkpo4(io35,isrna,ir1,ir2,ir3,id1,id2,id3,
     &                  ianz,iconn,ityp)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      integer*2 ityp
      dimension icn(3),ianz(*),iconn(mxcon+1,*),ityp(*)

      ip = 0
      if (io35.eq.0) return

      do i=1,iconn(1,io35)
         j = iconn(1+i,io35)
         if (j.gt.0) then
            if (ianz(j).eq.15) ip = j
         endif
      end do

      if (ip.gt.0) then
         no = 0
         do i=1,iconn(1,ip)
            j = iconn(1+i,ip)

            if (j.gt.0) then
               if (ianz(j).eq.8.and.j.ne.io35) then
                  ido = 1

                  do k=1,iconn(1,j)
                     l = iconn(1+k,j)
                     if (l.gt.0) then
                        if (l.gt.1.and.l.ne.ip) ido = 0
                     endif
                  end do

                  if (ido.eq.1) then
                     no = no + 1
                     icn(no) = j
                  endif

               endif
            endif
         end do

         if (no.eq.3) then

            if (isrna.eq.1) then

c     R-3'-Phosphate: O3* 1239, P 1240 OP 1241

               call settyp(ityp,io35,ir1)
               call settyp(ityp,ip,ir2)
               call setchg(ip,1)

               do i=1,3
                  call settyp(ityp,icn(i),ir3)
                  call setchg(icn(i),1)
               end do

            else

c     D-3'-Phosphate: O3* 1251, P 1252 OP 1253

               call settyp(ityp,io35,id1)
               call settyp(ityp,ip,id2)
               call setchg(ip,1)

               do i=1,3
                  call settyp(ityp,icn(i),id3)
                  call setchg(icn(i),1)
               end do

            endif
         endif
      endif

      return
      end

      subroutine typado(ipdb,jres,ihpdb,ihashy,
     &                  ianz,iconn,ityp)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxres=42)
      parameter (mxrss=20)
      parameter (mxrsa=23)
      parameter (mxata=9)
      parameter (mxatha=9)
      parameter (mxrso=24)
      parameter (mxato=10)
      parameter (mxatho=11)
      parameter (mxrsn=19)
      parameter (mxnucl=39)
      parameter (mxhnuc=27)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxambc=49)
      logical nterm,cterm
      integer*2 ityp
      common /types/ iff
      integer amboff,ntca,ctca
      common /ambtyp/ amboff(mxrss*3),ntca(mxrss),ctca(mxrss),
     &                ncca(mxrsa),icca(2,mxata,mxrsa),
     &                nhha(mxrsa),ihha(2,mxatha,mxrsa),
     &                ncco(mxrso),icco(2,mxato,mxrso),
     &                nhho(mxrso),ihho(2,mxatho,mxrso),
     &                nnuc(mxrsn),nhnuc(mxrsn),irna(mxres-23),
     &                inuc(2,mxnucl,mxrsn),ihnuc(2,mxhnuc,mxrsn)

      integer ambvdt
      common /fcharg/ ambchg(mxamb),ambvw1(mxambc),ambvw2(mxambc),
     &                gfvdw(2,mxgff),ambvdt(mxamb),cysneg(9)

      dimension ipdb(*), ihpdb(*),ianz(*),iconn(mxcon+1,*),ityp(*)

      iff = 4

c     currently only works for aminoacids, not nucleic acids 
c     (mxrss=20) (HIP,HID,HIE and sulfer bridge covered)
c                (ORN, mALA and pyroglut not covered)


      if (jres.le.0.or.jres.gt.mxrss) return

      do i=1,mxsym
         call settyp(ityp,ipdb(i),0)
         call setchg(ipdb(i),0)
      end do
      do i=1,mxhsym*3
         call settyp(ityp,ihpdb(i),0)
         call setchg(ihpdb(i),0)
      end do

      nterm = .false.
      n = 0
      do i=1,3
         if (ihpdb(i).ne.0) n = n + 1
      end do
      if (n.eq.3.or.(jres.eq.15.and.n.eq.2)) nterm = .true.
      cterm = .false.
      if (ipdb(38).ne.0) cterm = .true.

      ihis = 0
      if (jres.eq.17) then
c ihis
c 
c  1  HIS+ (HIP)
c  2  HISD (HID)
c  3  HISE (HIE)
c
         ihis = 1
         if (ihashy.eq.1) then
            if (ihpdb(34).eq.0) ihis = 2
            if (ihpdb(22).eq.0) ihis = 3
         endif
      endif

      iss = 0
      if (jres.eq.4) then
c iss = 0 (cysteine -SH), iss = 1 (cystine -SS-)
         do i=1,iconn(1,ipdb(37))
            if (ianz(abs(iconn(1+i,ipdb(37)))).eq.16) 
     &          iss = 1
         end do
      endif

      iasph = 0
      iasphh = 0
      if (jres.eq.9) then
c iasph is it aspartate or aspartic acid
         if (ipdb(28).ne.0) then
            do i=1,iconn(1,ipdb(28))
               if (ianz(abs(iconn(1+i,ipdb(28)))).eq.1) then
                   iasph = 28
                   iasphh = iabs(iconn(1+i,ipdb(28)))
               endif
            end do
         endif
         if (ipdb(29).ne.0) then
            do i=1,iconn(1,ipdb(29))
               if (ianz(abs(iconn(1+i,ipdb(29)))).eq.1) then
                   iasph = 29
                   iasphh = iabs(iconn(1+i,ipdb(29)))
               endif
            end do
         endif
         if (ipdb(30).ne.0) then
            do i=1,iconn(1,ipdb(30))
               if (ianz(abs(iconn(1+i,ipdb(30)))).eq.1) then
                   iasph = 30
                   iasphh = iabs(iconn(1+i,ipdb(30)))
               endif
            end do
         endif
      endif

      ioff = 7

      if (jres.eq.1) ioff = 1
      if (jres.eq.15) ioff = 50

c N
      if (nterm) then
         if (jres.eq.15) then
            call settyp(ityp,ipdb(1),194)
         else
            call settyp(ityp,ipdb(1),190)
         endif
      else
         call settyp(ityp,ipdb(1),ioff)
      endif
c CA
      if (jres.eq.3) then
         call settyp(ityp,ipdb(2),33)
      elseif (jres.eq.4) then
         call settyp(ityp,ipdb(2),44)
      elseif (jres.eq.15.and.nterm) then
         call settyp(ityp,ipdb(2),196)
      else
         call settyp(ityp,ipdb(2),ioff+1)
      endif
         
c C
      call settyp(ityp,ipdb(3),ioff+2)
      if (cterm) call settyp(ityp,ipdb(3),192)
      if (jres.eq.15.and.nterm) call settyp(ityp,ipdb(3),197)

c HN
      if (nterm) then
         if (jres.eq.15) then
            call settyp(ityp,ihpdb(1),195)
            call settyp(ityp,ihpdb(2),195)
            call settyp(ityp,ihpdb(3),195)
         else
            call settyp(ityp,ihpdb(1),191)
            call settyp(ityp,ihpdb(2),191)
            call settyp(ityp,ihpdb(3),191)
         endif
         
      else
         call settyp(ityp,ihpdb(1),ioff+3)
         call settyp(ityp,ihpdb(2),ioff+3)
         call settyp(ityp,ihpdb(3),ioff+3)
      endif

c O
      if (jres.eq.15.and..not.cterm) then
         if (nterm) then
            call settyp(ityp,ipdb(4),198)
         else
            call settyp(ityp,ipdb(4),ioff+3)
         endif
      else
         if (cterm) then
            call settyp(ityp,ipdb(4),193)
            call settyp(ityp,ipdb(38),193)
         else
            call settyp(ityp,ipdb(4),ioff+4)
         endif
      endif


c HA
      if (jres.eq.15) then
         if (nterm) then
            call settyp(ityp,ihpdb(4),199)
            call settyp(ityp,ihpdb(5),199)
         else
            call settyp(ityp,ihpdb(4),ioff+4)
            call settyp(ityp,ihpdb(5),ioff+4)
         endif
      else
         call settyp(ityp,ihpdb(4),ioff+5)
         call settyp(ityp,ihpdb(5),ioff+5)
      endif

      jtmp = jres
      if (jres.eq.17) then
         if (ihis.eq.2) jtmp = 22
         if (ihis.eq.3) jtmp = 23
      endif
      if (jres.eq.4.and.iss.eq.1) jtmp = 21
      if (jres.eq.9.and.iasph.ne.0) jtmp = 24
    
      do i=1,ncco(jtmp)
         call settyp(ityp,ipdb(icco(1,i,jtmp)),icco(2,i,jtmp))
      end do

      do i=1,nhho(jtmp)
         call settyp(ityp,ihpdb(ihho(1,i,jtmp)),ihho(2,i,jtmp))
      end do

      if (jres.eq.15) then
         if (nterm) then
            call settyp(ityp,ipdb(9),200)
            call settyp(ityp,ihpdb(19),201)
            call settyp(ityp,ihpdb(20),201)
         endif
      endif

      if (jres.eq.9.and.iasph.ne.0) then
         call settyp(ityp,ipdb(iasph),136)
         if (iasphh.ne.0) call settyp(ityp,iasphh,137)
c         call settyp(ityp,ihpdb(19),137)
c         call settyp(ityp,ihpdb(20),137)
c         call settyp(ityp,ihpdb(21),137)
      endif

c check for lysine alternative HZ1,HZ2,HZ3 instead of 1HZ,2HZ,3HZ

      if (jres.eq.12) then
         call settyp(ityp,ihpdb(43),168)
         call settyp(ityp,ihpdb(46),168)
         call settyp(ityp,ihpdb(49),168)
      endif

c      do i=1,mxsym
c         call setchg(ipdb(i),1)
c      end do
c      do i=1,mxhsym*3
c         call setchg(ihpdb(i),1)
c      end do

      return
      end

      subroutine settyp(itypa,iat,ityp)
      implicit double precision (a-h,o-z)
      integer*2 itypa
      dimension itypa(*)
      common /typoni/ ioniad

      if (iat.gt.0) then
         if (ioniad.gt.0) then
            itypa(iat) =  itypa(iat) + ityp
         else
            itypa(iat) =  ityp
         endif
      endif

      return
      end

      subroutine setchd(iat,iopt,qat,ityp)
      implicit double precision (a-h,o-z)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxambc=49)
      integer*2 ityp,it,it10000
      common /typoni/ ioniad
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      integer ambvdt
      common /fcharg/ ambchg(mxamb),ambvw1(mxambc),ambvw2(mxambc),
     &                gfvdw(2,mxgff),ambvdt(mxamb),cysneg(9)
      dimension qat(*),ityp(*)

      it10000=10000


      if (iat.gt.0) then
         if (iopt.eq.1) then
            if (ioniad.gt.0) then
               it = mod(ityp(iat),it10000)
            else
               it = ityp(iat)
            endif
            if (it.gt.0) then
               qat(iat) = ambchg(int(it))
            endif
            ihasq = 1
         elseif (iopt.eq.0) then
            qat(iat) = 0.0d0
         elseif (iopt.lt.0) then
            qat(iat) = cysneg(abs(iopt))
         endif
      endif

      return
      end

      subroutine flagd(ihpdb,iat,isurf)
      implicit double precision (a-h,o-z)
      dimension ihpdb(*),isurf(*)

      if (ihpdb(iat).ne.0) isurf(ihpdb(iat)) = 1

      return
      end

      subroutine conad(ipdb,iat1,iat2,iop,iconn)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      dimension ipdb(*),iconn(mxcon+1,*)

      ia1 = ipdb(iat1)
      ia2 = ipdb(iat2)

      if (ia1.ne.0.and.ia2.ne.0) then
         if (iconn(1,ia1).lt.mxcon) then
            iconn(1,ia1) = iconn(1,ia1) + 1
            iconn(iconn(1,ia1)+1,ia1) = ia2
         endif
         if (iop.eq.0) then
            iconn(1,ia2) = 1
            iconn(2,ia2) = ia1
         else
            if (iconn(1,ia2).lt.mxcon) then
               iconn(1,ia2) = iconn(1,ia2) + 1
               iconn(iconn(1,ia2)+1,ia2) = ia1
            endif
         endif
      endif

      return
      end

      subroutine conatd(ipdb,ihpdb,iat1,iat2,iconn)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      dimension ipdb(*), ihpdb(*), iconn(mxcon+1,*)

      ia1 = ipdb(iat1)
      ia2 = ihpdb(iat2)

      if (ia1.ne.0.and.ia2.ne.0) then
         if (iconn(1,ia1).lt.mxcon) then
            iconn(1,ia1) = iconn(1,ia1) + 1
            iconn(iconn(1,ia1)+1,ia1) = ia2
         endif
         iconn(1,ia2) = 1
         iconn(2,ia2) = ia1
      endif

      return
      end

      subroutine mkbacd(ipdb,ihpdb,jres,icres,ihashy,idoconv,
     &                  iconn,coo,
     &                  icalf,ianf,islu,nchain,iamino)

      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxcon=10)
      logical newch
      dimension ipdb(*),ihpdb(*),tmp(3),coo(3,*),iconn(mxcon+1,*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*)


      fct = 1.0d0
      if (idoconv.eq.1) fct = 1.0d0 / (0.52917706d0*0.52917706d0)

c N-CA
      if (ipdb(1).ne.0.and.ipdb(2).ne.0) then
         iconn(1,ipdb(1)) = 1
         iconn(2,ipdb(1)) = ipdb(2)
      endif

c CA-N, CA-C
      if (ipdb(2).ne.0) then
         iconn(1,ipdb(2)) = 0
         if (ipdb(1).ne.0) then
            iconn(2,ipdb(2)) = ipdb(1)
            iconn(1,ipdb(2)) = 1
            if (ipdb(3).ne.0) then
               iconn(3,ipdb(2)) = ipdb(3)
               iconn(1,ipdb(2)) = 2
            endif
         else
            if (ipdb(3).ne.0) then
               iconn(2,ipdb(2)) = ipdb(3)
               iconn(1,ipdb(2)) = 1
            endif
         endif
      endif

c C-CA,C-O
      if (ipdb(4).ne.0.and.ipdb(3).ne.0) then
         iconn(1,ipdb(3)) = 0
         if (ipdb(2).ne.0) then
            iconn(2,ipdb(3)) = ipdb(2)
            iconn(3,ipdb(3)) = ipdb(4)
            iconn(1,ipdb(3)) = 2
         else
            iconn(2,ipdb(3)) = ipdb(4)
            iconn(1,ipdb(3)) = 1
         endif
c O-C
         iconn(1,ipdb(4)) = 1
         iconn(2,ipdb(4)) = ipdb(3)

      elseif (ipdb(3).ne.0.and.ipdb(2).ne.0) then
         iconn(1,ipdb(3)) = 1
         iconn(2,ipdb(3)) = ipdb(2)
      elseif (ipdb(3).ne.0) then
         iconn(1,ipdb(3)) = 0
      elseif (ipdb(4).ne.0) then
         iconn(1,ipdb(4)) = 0
      endif


c Possible OXT
      if (ipdb(38).ne.0.or.ipdb(76).ne.0) then
         if (ipdb(38).ne.0) then
            iox = ipdb(38)
         else
            iox = ipdb(76)
         endif
         iconn(1,ipdb(3)) = iconn(1,ipdb(3)) + 1
         iconn(iconn(1,ipdb(3))+1,ipdb(3)) = iox
         iconn(1,iox) = 1
         iconn(2,iox) = ipdb(3)
      endif
c C beta
      if (jres.gt.1.and.ipdb(5).ne.0) then
         iconn(1,ipdb(2)) = 3
         iconn(4,ipdb(2)) = ipdb(5)
         iconn(1,ipdb(5)) = 1
         iconn(2,ipdb(5)) = ipdb(2)
      endif
 
      if (icres.le.numcal) then
         icalf(1,icres) = ipdb(2)
         icalf(2,icres) = ipdb(1)
         icalf(3,icres) = ipdb(3)
         icalf(4,icres) = 0
      endif

      if (icres.gt.1.and.icres.le.numcal) then
         newch = .false.
         ic = icalf(3,icres-1)
         in = ipdb(1)
         if (iamino(icres-1).le.23.and.(ic.gt.0.and.in.gt.0)) then
            do i=1,3
               tmp(i) = coo(i,in) - coo(i,ic)
            end do
            distsq = tmp(1)*tmp(1) + tmp(2)*tmp(2) +
     &               tmp(3)*tmp(3)
            if (distsq.lt.3.1684d0*fct) then
c
c connect N current residue to C previous
c
               iconn(1,in) = iconn(1,in) + 1 
               iconn(iconn(1,in)+1,in) = ic
c
c connect C previous to N current residue 
c
               iconn(1,ic) = iconn(1,ic) + 1 
               iconn(iconn(1,ic)+1,ic) = in
            else
               newch = .true.
            endif
         else 
            newch = .true.
         endif

         if (newch.and.idoconv.eq.0) then
            if (nchain.lt.mxchai) then
               islu(nchain) = icres-1
               nchain = nchain + 1
               ianf(nchain) = icres
            endif
         endif

      endif

      if (ihashy.eq.1) then
          do i=1,3
             call conath(ipdb,ihpdb,1,i)
             call conath(ipdb,ihpdb,2,3+i)
             call conath(ipdb,ihpdb,5,6+i)
          end do
          call conath(ipdb,ihpdb,1,79)
          call conath(ipdb,ihpdb,1,82)
          call conath(ipdb,ihpdb,1,91)
      endif

      return
      end

      subroutine mknbcd(ipdb,ihpdb,jres,icres,ihashy,idoconv,
     &                  iconn,coo,
     &                  icalf,ianf,islu,nchain,iamino)
      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxcon=10)
      logical newch
      dimension ipdb(*),ihpdb(*),tmp(3),coo(3,*),iconn(mxcon+1,*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*)

      fct = 1.0d0
      if (idoconv.eq.1) fct = 1.0d0 / (0.52917706d0*0.52917706d0)

      if (ipdb(43).gt.0) then
         nc = 0
         if (ipdb(44).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(43)) = ipdb(44)
         endif
         if (ipdb(45).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(43)) = ipdb(45)
         endif
         if (ipdb(46).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(43)) = ipdb(46)
         endif
         iconn(1,ipdb(43)) = nc
      endif

      if (ipdb(44).gt.0) iconn(1,ipdb(44)) = 0
      if (ipdb(45).gt.0) iconn(1,ipdb(45)) = 0
      if (ipdb(46).gt.0) iconn(1,ipdb(46)) = 0

      if (ipdb(43).gt.0) then

         if (ipdb(44).gt.0) then
            iconn(2,ipdb(44)) = ipdb(43)
            iconn(1,ipdb(44)) = 1
         endif

         if (ipdb(45).gt.0) then
            iconn(2,ipdb(45)) = ipdb(43)
            iconn(1,ipdb(45)) = 1
         endif

         if (ipdb(46).gt.0) then
            iconn(2,ipdb(46)) = ipdb(43)
            iconn(1,ipdb(46)) = 1
         endif
      endif


      if (ipdb(47).gt.0) then
         if (ipdb(46).gt.0) then
            iconn(1,ipdb(46)) = iconn(1,ipdb(46)) + 1
            iconn(1+iconn(1,ipdb(46)),ipdb(46)) = ipdb(47)
            iconn(2,ipdb(47)) = ipdb(46)
            iconn(1,ipdb(47)) = 1
         endif
         if (ipdb(48).gt.0) then
            if (ipdb(46).gt.0) then
               iconn(1,ipdb(47)) = 2
               iconn(3,ipdb(47)) = ipdb(48)
            else
               iconn(1,ipdb(47)) = 1
               iconn(2,ipdb(47)) = ipdb(48)
            endif
         endif
      endif

      if (ipdb(48).gt.0) then
         nc = 0
         if (ipdb(47).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(48)) = ipdb(47)
         endif
         if (ipdb(49).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(48)) = ipdb(49)
         endif
         if (ipdb(50).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(48)) = ipdb(50)
         endif
         iconn(1,ipdb(48)) = nc
      endif

      if (ipdb(49).gt.0) then
         nc = 0
         if (ipdb(48).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(49)) = ipdb(48)
         endif
         if (ipdb(54).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(49)) = ipdb(54)
         endif
         iconn(1,ipdb(49)) = nc
      endif

      if (ipdb(49).gt.0) then
         nc = 0
         if (ipdb(48).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(50)) = ipdb(48)
         endif
         if (ipdb(51).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(50)) = ipdb(51)
         endif
         if (ipdb(52).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(50)) = ipdb(52)
         endif
         iconn(1,ipdb(50)) = nc
      endif

      if (ipdb(50).gt.0.and.ipdb(51).gt.0) then
         iconn(1,ipdb(51)) = 1
         iconn(2,ipdb(51)) = ipdb(50)
      endif

      if (ipdb(52).gt.0) then
         nc = 0
         if (ipdb(50).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(52)) = ipdb(50)
         endif
         if (ipdb(54).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(52)) = ipdb(54)
         endif
         iconn(1,ipdb(52)) = nc
      endif

      if (ipdb(53).gt.0.and.ipdb(52).gt.0) then
         nc = iconn(1,ipdb(52))
         nc = nc + 1
         iconn(1+nc,ipdb(52)) = ipdb(53)
         iconn(1,ipdb(52)) = nc
         iconn(1,ipdb(53)) = 1
         iconn(2,ipdb(53)) = ipdb(52)
      endif

      if (ipdb(54).ne.0) then
         nc = 0
         if (ipdb(54).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(54)) = ipdb(49)
         endif
         if (ipdb(54).gt.0) then
            nc = nc + 1
            iconn(1+nc,ipdb(54)) = ipdb(52)
         endif
         iconn(1,ipdb(54)) = nc
      endif

c Possible OXT

      if (ipdb(38).ne.0.and.ipdb(43).ne.0) then
         nc = iconn(1,ipdb(43))
         nc = nc + 1
         iconn(1+nc,ipdb(43)) = ipdb(38)
         iconn(1,ipdb(43)) = nc
         iconn(1,ipdb(38)) = 1
         iconn(2,ipdb(38)) = ipdb(43)
      endif

c Possible O3P (OXT is old O3P new)

      if (ipdb(76).ne.0.and.ipdb(43).ne.0) then
         nc = iconn(1,ipdb(43))
         nc = nc + 1
         iconn(1+nc,ipdb(43)) = ipdb(76)
         iconn(1,ipdb(43)) = nc
         iconn(1,ipdb(76)) = 1
         iconn(2,ipdb(76)) = ipdb(43)
      endif

      if (jres.eq.36.or.jres.eq.31) then
          call conat(ipdb,53,73,0)
          call conat(ipdb,53,80,0)
      endif

      if (ihashy.eq.1) then
         call conath(ipdb,ihpdb,47,61)
         call conath(ipdb,ihpdb,47,62)
         call conath(ipdb,ihpdb,48,64)
         call conath(ipdb,ihpdb,50,67)
         call conath(ipdb,ihpdb,52,70)
         call conath(ipdb,ihpdb,52,71)
         call conath(ipdb,ihpdb,54,73)
         if (ipdb(53).ne.0) call conath(ipdb,ihpdb,53,76)
         if ((jres.eq.36.or.jres.eq.31).and.ipdb(80).ne.0) then
             call conath(ipdb,ihpdb,80,151)
             call conath(ipdb,ihpdb,80,152)
             call conath(ipdb,ihpdb,80,153)
         endif
         call conath(ipdb,ihpdb,51,115)
         call conath(ipdb,ihpdb,46,118)
      endif

      if (icres.le.numcal) then
         icalf(1,icres) = ipdb(43)
         icalf(2,icres) = ipdb(46)
         icalf(3,icres) = ipdb(47)
         icalf(4,icres) = ipdb(48)
         icalf(5,icres) = ipdb(50)
         icalf(6,icres) = ipdb(51)
      endif

      if (ipdb(43).eq.0.and.icres.le.numcal) then
         if (ihpdb(118).ne.0) then
             icalf(1,icres) = ihpdb(118)
         elseif (ipdb(46).ne.0) then
             icalf(1,icres) = ipdb(46)
         endif
      endif

      if (icres.gt.1.and.icres.le.numcal) then
         newch = .false.
         if (iamino(icres-1).gt.23) then
            ip = ipdb(43)
            io = icalf(6,icres-1)
            if (io.gt.0.and.ip.gt.0) then
               do i=1,3
                  tmp(i) = coo(i,ip) - coo(i,io)
               end do
               distsq = tmp(1)*tmp(1) + tmp(2)*tmp(2) +
     &                  tmp(3)*tmp(3)
            else
               distsq = 100000.0d0
            endif
            if (distsq.lt.3.1684d0*fct) then
c
c connect P current residue to O previous
c
               iconn(1,ip) = iconn(1,ip) + 1 
               iconn(iconn(1,ip)+1,ip) = io
c
c connect O previous to P current residue 
c
               iconn(1,io) = iconn(1,io) + 1 
               iconn(iconn(1,io)+1,io) = ip
            else
               newch = .true.
            endif
         else
            newch = .true.
         endif

         if (newch.and.idoconv.eq.0) then
            if (nchain.lt.mxchai) then
               islu(nchain) = icres-1
               nchain = nchain + 1
               ianf(nchain) = icres
            endif
         endif

      endif

      return
      end

      subroutine addhd(ires,jres,ipdb,ihpdb,nterm,
     &                 ianz,iaton,iatclr,iresid,iconn,isurf,ipdbt,ityp,
     &                 ncalf,icalf,coo)
c
c Add hydrogens
c Only the Hydrogens attached to the backbone N are not covered
c This is done by routine hcoord
c
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxres=42)
      parameter (mxhs=15)
      parameter (mxyg=26)
      common /hcom/  iat(4,mxhs,mxres),nrs(mxres),rint(3,mxhs,mxres),
     &               iatyg(4,mxyg),rintyg(3,mxyg)
      common /athlp/ iatoms, mxnat
      integer*2 ipdbt,ityp
      logical addat, o5pcon
      dimension ipdb(*), ihpdb(*), it(4), rt(4), icalf(6,*)
      dimension ianz(*),iaton(*),iatclr(*),iconn(mxcon+1,*),
     &          iresid(*),coo(3,*),ipdbt(*),ityp(*),isurf(*)
      data nrs /2,4,4,4,6,10,8,8,3,5,10,12,5,7,7,12,6,8,8,9,3*0,
     &          11,11,11,12,10,13,13,13,13,13,15,15,13,0,10,0,12,
     &          12,10/
c gly
      data ((iat(i,j,1),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 3,1,2,5, 52*0/  
      data ((rint(i,j,1),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 0.9779d0,107.66d0,-118.0d0, 
     &     39*0.0d0/  
c ala
      data ((iat(i,j,2),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 3,2,5,7, 3,2,5,8, 3,2,5,9,
     &     44*0/  
      data ((rint(i,j,2),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 0.990d0,109.8d0,180.0d0, 
     &     0.990d0,109.8d0,60.0d0, 0.990d0,109.8d0,-60.0d0,
     &     33*0.0d0/  
c ser
      data ((iat(i,j,3),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 2,5,31,10, 31,2,5,7, 31,2,5,8,
     &     44*0/  
      data ((rint(i,j,3),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 0.867d0,107.06d0,180.0d0, 
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     33*0.0d0/  
c cys
      data ((iat(i,j,4),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 2,5,37,10, 37,2,5,7, 37,2,5,8,
     &     44*0/  
      data ((rint(i,j,4),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 0.867d0,107.06d0,180.0d0, 
     &     0.990d0,112.5d0,120.0d0, 0.990d0,112.5d0,-120.0d0,
     &     33*0.0d0/  
c thr
      data ((iat(i,j,5),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 32,8,5,7, 2,5,32,13, 
     &     2,5,8,16, 2,5,8,17, 2,5,8,18,
     &     36*0/  
      data ((rint(i,j,5),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 0.990d0,109.08d0,-120.0d0, 
     &     0.990d0,112.5d0,-120.0d0, 0.978d0,110.0d0,180.0d0,
     &     0.978d0,110.0d0,60.0d0, 0.978d0,110.0d0,-60.0d0,
     &     27*0.0d0/  
c ile
      data ((iat(i,j,6),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 7,8,5,7,  7,5,8,16, 7,5,8,17, 7,5,8,18,
     &     10,5,7,13, 10,5,7,14, 5,7,10,22, 5,7,10,23, 5,7,10,24,
     &     20*0/  
      data ((rint(i,j,6),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 0.990d0,109.08d0,-120.0d0, 
     &     0.978d0,110.0d0,180.0d0,
     &     0.978d0,110.0d0,60.0d0, 0.978d0,110.0d0,-60.0d0,
     &     0.978d0,110.0d0,120.0d0, 0.978d0,110.0d0,-120.0d0,
     &     0.990d0,109.8d0,180.0d0, 
     &     0.990d0,109.8d0,60.0d0, 0.990d0,109.8d0,-60.0d0,
     &     15*0.0d0/  
c val
      data ((iat(i,j,7),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 7,8,5,7,  7,5,8,16, 7,5,8,17, 7,5,8,18,
     &     8,5,7,13, 8,5,7,14, 8,5,7,15,
     &     28*0/  
      data ((rint(i,j,7),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 0.990d0,109.08d0,120.0d0, 
     &     0.978d0,110.0d0,180.0d0,
     &     0.978d0,110.0d0,60.0d0, 0.978d0,110.0d0,-60.0d0,
     &     0.990d0,109.8d0,180.0d0, 
     &     0.990d0,109.8d0,60.0d0, 0.990d0,109.8d0,-60.0d0,
     &     21*0.0d0/  
c met
      data ((iat(i,j,8),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7,  6,2,5,8,
     &     36,5,6,10, 36,5,6,11, 
     &     6,36,12,28, 6,36,12,29, 6,36,12,30,
     &     28*0/  
      data ((rint(i,j,8),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 0.990d0,109.08d0,120.0d0, 
     &     0.990d0,109.08d0,-120.0d0,
     &     0.978d0,110.0d0,120.0d0, 0.978d0,110.0d0,-120.0d0,
     &     0.990d0,109.8d0,180.0d0, 
     &     0.990d0,109.8d0,60.0d0, 0.990d0,109.8d0,-60.0d0,
     &     21*0.0d0/  
c asp
      data ((iat(i,j,9),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7, 6,2,5,8,
     &     48*0/  
      data ((rint(i,j,9),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     36*0.0d0/  
c asn
      data ((iat(i,j,10),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7, 6,2,5,8,
     &     29,6,21,25, 29,6,21,26,
     &     40*0/  
      data ((rint(i,j,10),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     0.950d0,120.0d0,180.0d0, 0.950d0,120.0d0,0.0d0,
     &     30*0.0d0/  
c leu
      data ((iat(i,j,11),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7, 6,2,5,8,
     &     5,10,6,10, 5,6,10,22, 5,6,10,23, 5,6,10,24,
     &     5,6,11,25, 5,6,11,26, 5,6,11,27,
     &     20*0/  
      data ((rint(i,j,11),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     0.9779d0,107.63d0,120.0d0,
     &     0.990d0,109.8d0,180.0d0, 
     &     0.990d0,109.8d0,60.0d0, 0.990d0,109.8d0,-60.0d0,
     &     0.990d0,109.8d0,180.0d0, 
     &     0.990d0,109.8d0,60.0d0, 0.990d0,109.8d0,-60.0d0,
     &     15*0.0d0/  
c lys
      data ((iat(i,j,12),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7, 6,2,5,8, 9,5,6,10, 9,5,6,11,
     &     12,6,9,19, 12,6,9,20, 27,9,12,28, 27,9,12,29,
     &     9,12,27,40, 9,12,27,41, 9,12,27,42,
     &     12*0/
      data ((rint(i,j,12),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     0.990d0,109.8d0,180.0d0, 
     &     0.990d0,109.8d0,60.0d0, 0.990d0,109.8d0,-60.0d0,
     &     9*0/
c glu
      data ((iat(i,j,13),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7, 6,2,5,8, 9,5,6,10, 9,5,6,11,
     &     40*0/
      data ((rint(i,j,13),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     30*0/
c gln
      data ((iat(i,j,14),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7, 6,2,5,8, 9,5,6,10, 9,5,6,11,
     &     34,9,24,34, 34,9,24,35,
     &     32*0/
      data ((rint(i,j,14),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     0.950d0,120.0d0,180.0d0, 0.950d0,120.0d0,0.0d0,
     &     24*0/
c pro
      data ((iat(i,j,15),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7, 6,2,5,8, 9,5,6,10, 9,5,6,11,
     &     1,6,9,19, 1,6,9,20,
     &     32*0/
      data ((rint(i,j,15),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 
     &     0.990d0,109.4d0,-120.0d0, 0.990d0,109.4d0,120.0d0,
     &     0.990d0,109.4d0,-120.0d0, 0.990d0,109.4d0,120.0d0,
     &     0.990d0,109.4d0,-120.0d0, 0.990d0,109.4d0,120.0d0,
     &     24*0/
c arg
      data ((iat(i,j,16),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7, 6,2,5,8, 9,5,6,10, 9,5,6,11,
     &     22,6,9,19, 22,6,9,20, 9,17,22,28,
     &     22,17,25,55, 22,17,25,56, 22,17,26,58, 22,17,26,59,
     &     12*0/
      data ((rint(i,j,16),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     0.990d0,112.5d0,-120.0d0, 0.990d0,112.5d0,120.0d0,
     &     1.020d0,120.0d0,180.0d0,
     &     1.020d0,120.0d0,0.0d0, 1.020d0,120.0d0,180.0d0,
     &     1.020d0,120.0d0,0.0d0, 1.020d0,120.0d0,180.0d0,
     &     9*0/
c his
      data ((iat(i,j,17),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7, 6,2,5,8, 
     &     6,13,20,22, 20,24,13,31, 24,6,11,25,
     &     36*0/
      data ((rint(i,j,17),i=1,3),j=1,mxhs) /
     &     0.9779d0,107.63d0,118.0d0, 
     &     1.110d0,109.4d0,-120.0d0, 1.110d0,109.4d0,120.0d0,
     &     1.020d0,126.0d0,180.0d0, 1.100d0,126.0d0,180.0d0,
     &     1.020d0,126.0d0,180.0d0, 1.100d0,126.0d0,180.0d0,
     &     24*0/
c phe
      data ((iat(i,j,18),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7, 6,2,5,8, 
     &     6,13,10,22, 10,17,13,31, 13,14,17,40, 17,11,14,34,
     &     14,6,11,25,
     &     28*0/
      data ((rint(i,j,18),i=1,3),j=1,mxhs) /
     &     1.110d0,107.9d0,118.0d0, 
     &     1.110d0,109.4d0,-120.0d0, 1.110d0,109.4d0,120.0d0,
     &     1.100d0,120.0d0,180.0d0, 1.100d0,120.0d0,180.0d0,
     &     1.100d0,120.0d0,180.0d0, 1.100d0,120.0d0,180.0d0,
     &     1.100d0,120.0d0,180.0d0,
     &     21*0/
c tyr
      data ((iat(i,j,19),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7, 6,2,5,8, 
     &     6,13,10,22, 10,17,13,31, 13,17,33,52, 17,11,14,34,
     &     14,6,11,25,
     &     28*0/
      data ((rint(i,j,19),i=1,3),j=1,mxhs) /
     &     1.110d0,107.9d0,118.0d0, 
     &     1.110d0,109.4d0,-120.0d0, 1.110d0,109.4d0,120.0d0,
     &     1.100d0,120.0d0,180.0d0, 1.100d0,120.0d0,180.0d0,
     &     0.970d0,108.0d0,0.0d0, 1.100d0,120.0d0,180.0d0,
     &     1.100d0,120.0d0,180.0d0,
     &     21*0/
c trp
      data ((iat(i,j,20),i=1,4),j=1,mxhs) /
     &     3,1,2,4, 6,2,5,7, 6,2,5,8, 
     &     6,23,10,22, 10,14,23,31, 14,16,18,46, 18,19,16,58,
     &     16,15,19,49, 11,19,15,37,
     &     24*0/
      data ((rint(i,j,20),i=1,3),j=1,mxhs) /
     &     1.110d0,107.9d0,118.0d0, 
     &     1.110d0,109.4d0,-120.0d0, 1.110d0,109.4d0,120.0d0,
     &     1.100d0,124.0d0,180.0d0, 1.050d0,124.0d0,180.0d0,
     &     1.100d0,120.0d0,180.0d0, 1.100d0,120.0d0,180.0d0,
     &     1.100d0,120.0d0,180.0d0, 1.100d0,120.0d0,180.0d0,
     &     18*0/

c asx,glx,hyp not yet

      data (((iat(i,j,k),i=1,4),j=1,mxhs),k=21,23) /180*0/
      data (((rint(i,j,k),i=1,3),j=1,mxhs),k=21,23) /135*0.0d0/

c adenosine

      data ((iat(i,j,24),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     66,65,59,112,62,60,55,82, 60,58,64,-136, 60,58,64,-137,
     &     16*0/
      data ((rint(i,j,24),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,122.0d0,180.0d0, 1.080d0,116.0d0,180.0d0,
     &     1.010d0,120.0d0,180.0d0, 1.010d0,120.0d0,0.0d0,
     &     12*0/

c cytidine

      data ((iat(i,j,25),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     60,57,58,103,58,56,57,100, 57,56,63,-130, 57,56,63,-131,
     &     16*0/
      data ((rint(i,j,25),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,118.0d0,180.0d0, 1.080d0,123.4d0,180.0d0,
     &     1.010d0,120.0d0,180.0d0, 1.010d0,120.0d0,0.0d0,
     &     12*0/

c guanosine

      data ((iat(i,j,26),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     66,65,59,112,58,55,60,-121, 60,55,61,-124, 60,55,61,-125,
     &     16*0/
      data ((rint(i,j,26),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,122.0d0,180.0d0, 1.008d0,118.5d0,180.0d0,
     &     1.009d0,120.0d0,180.0d0, 1.009d0,120.0d0,0.0d0,
     &     12*0/

c thimidine

      data ((iat(i,j,27),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     60,57,58,103, 56,55,62,-127, 
     &     56,57,83,160, 56,57,83,161, 56,57,83,162,
     &     12*0/
      data ((rint(i,j,27),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,118.0d0,180.0d0, 1.010d0,117.0d0,180.0d0, 
     &     1.090d0,109.8d0,180.0d0,1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0,
     &     9*0/

c uridine

      data ((iat(i,j,28),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     60,57,58,103,58,56,57,100, 56,55,62,-127, 
     &     20*0/
      data ((rint(i,j,28),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,118.0d0,180.0d0, 1.080d0,119.4d0,180.0d0,
     &     1.010d0,117.0d0,180.0d0, 
     &     15*0/

c 1MA

      data ((iat(i,j,29),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     66,65,59,112,62,60,55,82, 
     &     60,58,64,-137, 58,60,79,148, 58,60,79,149, 58,60,79,150, 
     &     0,0,0,0, 0,0,0,0/
      data ((rint(i,j,29),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0, 
     &     1.080d0,120.0d0,180.0d0,1.080d0,120.0d0,180.0d0,
     &     1.010d0,120.0d0,180.0d0,
     &     1.090d0,109.8d0,180.0d0, 1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0,
     &     1.080d0,90.0d0,90.0d0,0.0d0,.0d0,0.0d0/

c 5MC

      data ((iat(i,j,30),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     60,57,58,103, 57,56,63,-130, 57,56,63,-131, 
     &     56,57,83,160, 56,57,83,161, 56,57,83,162,
     &     8*0/
      data ((rint(i,j,30),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0, 1.080d0,118.0d0,180.0d0,
     &     1.010d0,120.0d0,180.0d0, 1.010d0,120.0d0,0.0d0,
     &     1.090d0,109.8d0,180.0d0,1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0,
     &     6*0/

c OMC

      data ((iat(i,j,31),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 
     &     52,53,80,151, 52,53,80,152, 52,53,80,153,
     &     60,57,58,103,58,56,57,100, 57,56,63,-130, 57,56,63,-131,
     &     8*0/
      data ((rint(i,j,31),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,180.0d0, 1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0,
     &     1.080d0,118.0d0,180.0d0, 1.080d0,123.4d0,180.0d0,
     &     1.010d0,120.0d0,180.0d0, 1.010d0,120.0d0,0.0d0,
     &     6*0/

c 1MG

      data ((iat(i,j,32),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     66,65,59,112,
     &     58,60,79,148, 58,60,79,149, 58,60,79,150,
     &     60,55,61,-124, 60,55,61,-125,
     &     8*0/
      data ((rint(i,j,32),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0, 1.080d0,122.0d0,180.0d0, 
     &     1.090d0,109.8d0,180.0d0, 1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0,
     &     1.009d0,120.0d0,180.0d0, 1.009d0,120.0d0,0.0d0,
     &     6*0/

c 2MG

      data ((iat(i,j,33),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     66,65,59,112,58,55,60,-121, 60,55,61,-124, 
     &     55,61,80,151, 55,61,80,152, 55,61,80,153,
     &     8*0/
      data ((rint(i,j,33),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,122.0d0,180.0d0, 1.008d0,118.5d0,180.0d0,
     &     1.009d0,120.0d0,0.0d0, 
     &     1.090d0,109.8d0,180.0d0, 1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0,
     &     6*0/

c M2G

      data ((iat(i,j,34),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     66,65,59,112,58,55,60,-121, 
     &     55,61,79,148, 55,61,79,149, 55,61,79,150,
     &     55,61,80,151, 55,61,80,152, 55,61,80,153/
      data ((rint(i,j,34),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,122.0d0,180.0d0, 1.008d0,118.5d0,180.0d0,
     &     1.090d0,109.8d0,180.0d0, 1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0,
     &     1.090d0,109.8d0,180.0d0, 1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0/

c 7MG

      data ((iat(i,j,35),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     66,65,59,112,58,55,60,-121, 60,55,61,-124, 
     &     60,55,61,-125,57,65,85,166, 57,65,85,167, 57,65,85,168,
     &     0,0,0,0/
      data ((rint(i,j,35),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,122.0d0,180.0d0, 
     &     1.008d0,118.5d0,180.0d0,
     &     1.009d0,120.0d0,180.0d0, 1.009d0,120.0d0,0.0d0,
     &     1.090d0,109.8d0,180.0d0, 1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0,0.0d0,0.0d0,0.0d0/

c OMG

      data ((iat(i,j,36),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 
     &     52,53,80,151, 52,53,80,152, 52,53,80,153,
     &     66,65,59,112,58,55,60,-121, 60,55,61,-124, 60,55,61,-125,
     &     8*0/
      data ((rint(i,j,36),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,180.0d0, 1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0,
     &     1.080d0,122.0d0,180.0d0, 1.008d0,118.5d0,180.0d0,
     &     1.009d0,120.0d0,180.0d0, 1.009d0,120.0d0,0.0d0,
     &     6*0/

c NO YG

      data ((iat(i,j,37),i=1,4),j=1,mxhs) /60*0/
      data ((rint(i,j,37),i=1,3),j=1,mxhs) /45*0.0d0/

c inosine

      data ((iat(i,j,38),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     66,65,59,112,58,55,60,-121,62,60,55,82, 
     &     20*0/
      data ((rint(i,j,38),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,122.0d0,180.0d0, 1.008d0,118.5d0,180.0d0,
     &     1.080d0,122.0d0,180.0d0,
     &     15*0/

c NO +U

      data ((iat(i,j,39),i=1,4),j=1,mxhs) /60*0/
      data ((rint(i,j,39),i=1,3),j=1,mxhs) /45*0.0d0/

c H2U

      data ((iat(i,j,40),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     60,57,58,103,60,57,58,104,58,56,57,100,58,56,57,101,
     &     56,55,62,-127, 
     &     12*0/
      data ((rint(i,j,40),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,118.0d0,120.0d0, 1.080d0,119.4d0,-120.0d0,
     &     1.080d0,118.0d0,120.0d0, 1.080d0,119.4d0,-120.0d0,
     &     1.010d0,117.0d0,180.0d0, 
     &     9*0/

c 5MU

      data ((iat(i,j,41),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     60,57,58,103,56,55,62,-127, 
     &     56,57,83,160, 56,57,83,161, 56,57,83,162,
     &     12*0/
      data ((rint(i,j,41),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,118.0d0,180.0d0, 1.010d0,117.0d0,180.0d0, 
     &     1.090d0,109.8d0,180.0d0,1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0,
     &     9*0/


c PSU

      data ((iat(i,j,42),i=1,4),j=1,mxhs) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     60,57,58,103,55,58,60,-121, 56,55,62,-127, 
     &     20*0/
      data ((rint(i,j,42),i=1,3),j=1,mxhs) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,118.0d0,180.0d0, 1.010d0,117.0d0,180.0d0,
     &     1.010d0,117.0d0,180.0d0, 
     &     15*0/

c YG

      data ((iatyg(i,j),i=1,4),j=1,mxyg) /
     &     46,48,47,61, 46,48,47,62, 47,49,48,64,
     &     48,52,50,67, 50,54,52,70, 52,49,54,73, 50,52,53,76,
     &     66,65,59,112,56,62,88,91, 56,62,88,92, 56,62,88,93,
     &     61,90,89,175, 61,90,89,176, 61,90,89,177,
     &     91,93,92,178, 91,93,92,179,
     &     92,94,93,181, 92,94,93,182, 93,95,94,184,
     &     95,97,98,187, 95,97,98,188, 95,97,98,189,
     &     94,100,99,-124,
     &     100,102,103,190, 100,102,103,191, 100,102,103,192/
      data ((rintyg(i,j),i=1,3),j=1,mxyg) /
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,-120.0d0, 1.090d0,109.8d0,120.0d0,
     &     0.957d0,108.4d0,49.8d0,
     &     1.080d0,122.0d0,180.0d0, 1.090d0,109.8d0,180.0d0,
     &     1.090d0,109.8d0,60.0d0, 1.090d0,109.8d0,-60.0d0,
     &     1.090d0,109.8d0,180.0d0,1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,120.0d0, 1.090d0,109.8d0,-120.0d0,
     &     1.090d0,109.8d0,120.0d0,
     &     1.090d0,109.8d0,180.0d0,1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0, 1.010d0,120.0d0,180.0d0,
     &     1.090d0,109.8d0,180.0d0,1.090d0,109.8d0,60.0d0,
     &     1.090d0,109.8d0,-60.0d0/

c
c for nucleotides the H's that make hydrogen bonds are added
c by subroutine mkcon, addhy1 and addhy2, these atoms are flagged
c via array isurf, which is used by dohcon to create hydrogen bonds
c this routine does not do that
c
c OXT = O3P , do we have to add H's to O2P,O3P, at present NO
c

      if (jres.le.0.or.jres.gt.mxres) return

      if (jres.eq.37) then
         nh = mxyg
      else
         nh = nrs(jres)
      endif

      do i=1,nh

          if (jres.eq.37) then
              do j=1,4
                 it(j) = iatyg(j,i)
              end do
              do j=1,3
                 rt(j) = rintyg(j,i)
              end do
          else
              do j=1,4
                 it(j) = iat(j,i,jres)
              end do
              do j=1,3
                 rt(j) = rint(j,i,jres)
              end do
          endif

          ih = abs(it(4))
          ihp = ihpdb(ih)

          if (jres.eq.4.and.i.eq.2) then
              if (ipdb(37).gt.0.and.ipdb(37).lt.mxnat) then
                 if (iconn(1,ipdb(37)).gt.1) ihp = 1
              endif
          endif
          if (jres.eq.3.and.i.eq.2) then
              if (ipdb(31).gt.0.and.ipdb(31).lt.mxnat) then
                 if (iconn(1,ipdb(31)).gt.1) ihp = 1
              endif
          endif
          if (jres.eq.12.and.i.ge.10) then
              n = 0
              if (ipdb(27).gt.0.and.ipdb(27).lt.mxnat) then
                 n = iconn(1,ipdb(27)) - 1
              endif
              if (n.gt.2) ihp = 1
          endif
          if (jres.eq.17.and.i.eq.4) then
              if (ipdb(22).gt.0) then
                 if (iconn(1,ipdb(22)).gt.2) ihp = 1
              endif
          endif
          if (ihp.eq.0) then
              if (addat(ipdb(it(1)),ipdb(it(2)),
     &            ipdb(it(3)),1,rt(1),rt(2),rt(3),ihpdb(ih),
     &            ianz,iaton,iatclr,iconn,iresid,ityp,ipdbt,
     &            ncalf,icalf,coo)) then
                  iresid(ihpdb(ih)) = ires
                  ipdbt(ihpdb(ih)) = ih
                  if (it(4).lt.0) isurf(ihpdb(ih)) = 1
              endif
          endif
      end do

      if (jres.gt.23) then
c
c deal with DNA as well
c
         if (ipdb(53).eq.0) then
             if (addat(ipdb(50),ipdb(54),ipdb(52),1,
     &        1.090d0,109.8d0,120.0d0,ihpdb(71),
     &            ianz,iaton,iatclr,iconn,iresid,ityp,ipdbt,
     &            ncalf,icalf,coo)) then
                 iresid(ihpdb(71)) = ires
                 ipdbt(ihpdb(71)) = 71
             endif
         endif

c HO3

         if (ipdb(51).ne.0) then
             n = icred(ipdb(51),inoh,ih,ianz,iconn)
             if (inoh.eq.1.and.ih.eq.0) then
                if (addat(ipdb(52),ipdb(50),ipdb(51),1,
     &           0.960d0,109.8d0,120.0d0,ihpdb(115),
     &            ianz,iaton,iatclr,iconn,iresid,ityp,ipdbt,
     &            ncalf,icalf,coo)) then
                    iresid(ihpdb(115)) = ires
                    ipdbt(ihpdb(115)) = 115
                endif
             endif
         endif

         if (ipdb(43).eq.0.and.ipdb(46).ne.0) then
             if (.not.o5pcon(ipdb(46),ianz,iconn)) then
                if (addat(ipdb(48),ipdb(47),ipdb(46),1,
     &           0.960d0,109.8d0,120.0d0,ihpdb(118),
     &            ianz,iaton,iatclr,iconn,iresid,ityp,ipdbt,
     &            ncalf,icalf,coo)) then
                    iresid(ihpdb(118)) = ires
                    ipdbt(ihpdb(118)) = 118
                    icalf(1,ires) = ihpdb(118)
                endif
             endif
         endif

      else

         if (nterm.eq.1) then
            n = icred(ipdb(1),inoh,ih,ianz,iconn)
            if (inoh.eq.1.and.ih.eq.1.and.ihpdb(1).ne.0) then
               if (addat(ihpdb(1),ipdb(2),ipdb(1),1,
     &          1.010d0,109.8d0,120.0d0,ihpdb(2),
     &            ianz,iaton,iatclr,iconn,iresid,ityp,ipdbt,
     &            ncalf,icalf,coo)) then
                   iresid(ihpdb(2)) = ires
                   ipdbt(ihpdb(2)) = 2
               endif
               if (addat(ihpdb(1),ipdb(2),ipdb(1),1,
     &          1.010d0,109.8d0,-120.0d0,ihpdb(3),
     &            ianz,iaton,iatclr,iconn,iresid,ityp,ipdbt,
     &            ncalf,icalf,coo)) then
                   iresid(ihpdb(3)) = ires
                   ipdbt(ihpdb(3)) = 3
               endif
            endif
         endif

      endif

      return
      end

      logical function o5pcon(iat,ianz,iconn)
      parameter (mxcon=10)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      dimension ianz(*),iconn(mxcon+1,*)
   
      o5pcon = .false.

      do i=1,iconn(1,iat)
         ii = iconn(1+i,iat)
         if (ii.gt.0) then
            if (ianz(ii).eq.15) o5pcon = .true.
         endif
      end do
      
      return
      end

      subroutine rdprot
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxres=42)
      parameter (maxam=1000)
      character*137 line,string
      common /curlin/ line
      character*3 aminos
      common /amino/  aminos(mxres)
      integer getlin
      dimension iamin(maxam),angs(7,maxam)

      mxam = maxam
      namin = 0
      do while (getlin(0).eq.1)
         iam = 0
         ktype = nxtwrd(string,nstr,itype,rtype)
         if (ktype.eq.1.and.nstr.eq.3) then
             call tocap(string,3)
             do i=1,20
                if (string(1:3).eq.aminos(i)) iam = i
             end do
             if (iam.ne.0) then
                if (namin.lt.maxam) then
                   namin = namin + 1
                   iamin(namin) = iam - 1
                   do i=1,7
                      angs(i,namin) = -1
                      ktype = nxtwrd(string,nstr,itype,rtype)
                      if (ktype.eq.3) angs(i,namin) = rtype
                      if (ktype.eq.2) angs(i,namin) = dble(itype)
                   end do
                else
                   call inferr('To many RESIDUES !',0)
                endif
             endif
         endif
      end do

      call readsq(iamin,angs,namin)

      return
      end

      subroutine numhed(numhet,iresid)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      dimension iresid(*)

      nhmtmp = -3
      do i=1,iatoms
         if (iresid(i).lt.nhmtmp) nhmtmp = iresid(i)
      end do
      numhet = iabs(nhmtmp)

      return
      end

      logical function chkpdb(ipdb,jres,icres,irsnr)
      implicit double precision (a-h,o-z)
      parameter (mxres=42)
      parameter (mxares=20)
      parameter (mxaat=13)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxchtp=136)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxamo=201)
      character*3 aminos
      common /amino/aminos(mxres)
      integer resn,resat
      common /resa/ resat(mxares,mxaat),resn(mxares)
      character*3 pdbsym,hsym,chtnk,ambtnk
      character*2 amotnk,gffstr
      common /symbol/ pdbsym(mxsym),hsym(mxhsym),chtnk(mxchtp),
     &                ambtnk(mxamb),amotnk(mxamo),gffstr(mxgff)
      dimension ipdb(*),irsnr(*)

      data resn/4,4,6,6,7,8,7,8,8,8,8,9,9,9,8,11,10,11,12,13/
c glycine
      data (resat(1,i),i=1,mxaat)/
     &          1,2,3,4,0,0,0,0,0,0,0,0,0/
c alanine
      data (resat(2,i),i=1,mxaat)/
     &          1,2,3,4,0,0,0,0,0,0,0,0,0/
c serine
      data (resat(3,i),i=1,mxaat)/
     &          1,2,3,4,5,31,0,0,0,0,0,0,0/
c cysteine
      data (resat(4,i),i=1,mxaat)/
     &          1,2,3,4,5,37,0,0,0,0,0,0,0/
c threonine
      data (resat(5,i),i=1,mxaat)/
     &          1,2,3,4,5,32,8,0,0,0,0,0,0/
c isoleucine
      data (resat(6,i),i=1,mxaat)/
     &          1,2,3,4,5,7,8,10,0,0,0,0,0/
c valine
      data (resat(7,i),i=1,mxaat)/
     &          1,2,3,4,5,7,8,0,0,0,0,0,0/
c methionine
      data (resat(8,i),i=1,mxaat)/
     &          1,2,3,4,5,6,36,12,0,0,0,0,0/
c aspartic acid
      data (resat(9,i),i=1,mxaat)/
     &          1,2,3,4,5,6,29,30,0,0,0,0,0/
c asparagine
      data (resat(10,i),i=1,mxaat)/
     &          1,2,3,4,5,6,29,21,0,0,0,0,0/
c leucine
      data (resat(11,i),i=1,mxaat)/
     &          1,2,3,4,5,6,10,11,0,0,0,0,0/
c lysine
      data (resat(12,i),i=1,mxaat)/
     &          1,2,3,4,5,6,9,12,27,0,0,0,0/
c glutamic acid
      data (resat(13,i),i=1,mxaat)/
     &          1,2,3,4,5,6,9,34,35,0,0,0,0/
c glutamine
      data (resat(14,i),i=1,mxaat)/
     &          1,2,3,4,5,6,9,24,34,0,0,0,0/
c proline 28 ?
      data (resat(15,i),i=1,mxaat)/
     &          1,2,3,4,5,6,9,1,0,0,0,0,0/
c arginine
      data (resat(16,i),i=1,mxaat)/
     &          1,2,3,4,5,6,9,17,22,25,26,0,0/
c histidine 
      data (resat(17,i),i=1,mxaat)/
     &          1,2,3,4,5,6,11,13,20,24,0,0,0/
c phenylalanine
      data (resat(18,i),i=1,mxaat)/
     &          1,2,3,4,5,6,10,11,13,14,17,0,0/
c tyrosine
      data (resat(19,i),i=1,mxaat)/
     &          1,2,3,4,5,6,10,11,13,14,17,33,0/
c tryptophan
      data (resat(20,i),i=1,mxaat)/
     &          1,2,3,4,5,6,10,11,14,15,16,18,23/

      iresa = -1
      if (jres.le.20) then
         do i=1,resn(jres)
            if (ipdb(resat(jres,i)).le.0) then
               chkpdb = .false.
               print*,'incomplete residue ',irsnr(icres),' ',
     &                 aminos(jres),' ',pdbsym(resat(jres,i))
               return
            endif
         end do
      endif

      chkpdb = .true.
      return
      end
          
      subroutine detanz(resi,tstr,lstr,ifnd,ish,ianz,ihashy)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character lstr*137
      character*3 resi, atmp
      character*4 tstr
      character*2 tocapf

          ifnd = 1
          ish = 0
          atmp = tstr(2:4)
c
c Take care of ambiguity; asx->asn ,glx->gln ,his
c
          if (resi.eq.'ASX') then
              if (atmp.eq.'AD1') then
                  atmp = 'OD1'
              elseif (atmp.eq.'AD2') then
                  atmp = 'ND2'
              endif
          elseif (resi.eq.'GLX') then
              if (atmp.eq.'AE1') then
                  atmp = 'OE1'
              elseif (atmp.eq.'AE2') then
                  atmp = 'NE2'
              endif
          elseif (resi.eq.'HIS') then
              if (atmp.eq.'AD1') then
                  atmp = 'ND1'
              elseif (atmp.eq.'AD2') then
                  atmp = 'CD2'
              elseif (atmp.eq.'AE1') then
                  atmp = 'CE1'
              elseif (atmp.eq.'AE2') then
                  atmp = 'NE2'
              endif
          elseif (resi.eq.'ILE') then
              if (atmp.eq.'CD ') then
                  atmp = 'CD1'
              elseif (atmp.eq.'HD1') then
                  atmp = 'HD1'
                  tstr(1:4) = '1HD1'
              elseif (atmp.eq.'HD2') then
                  atmp = 'HD1'
                  tstr(1:4) = '2HD1'
              elseif (atmp.eq.'HD3') then
                  atmp = 'HD1'
                  tstr(1:4) = '3HD1'
              endif
          elseif (resi.eq.'SER'.or.resi.eq.'CYS') then
              if (atmp.eq.'HG1') atmp = 'HG '
          elseif (resi.eq.'  T'.or.resi.eq.' DT') then
              if (atmp(1:2).eq.'H7') atmp = 'HM7'
          endif
          if (atmp.ne.tstr(2:4)) tstr(2:4) = atmp


          if (tstr(1:2).eq.'HH'.or.tstr(1:2).eq.'HD'.or.
     &        tstr(1:2).eq.'HE'.or.tstr(1:2).eq.'HG') then
              atmp = tstr(1:3)
              tstr = ' '//atmp
          endif

c to unscrew whatif specific types

          if (tstr(2:4).eq.'O**') tstr(2:4) = 'OXT'
          if (tstr(2:4).eq.'O* ') tstr(2:4) = 'O  '
          if (tstr(2:4).eq.'H5M') tstr(2:4) = 'HM5'
          if (tstr(2:4).eq.'H3T') tstr(2:4) = 'HO3'
          if (tstr(1:4).eq.'1H4 ') tstr(1:4) = ' H41'
          if (tstr(1:4).eq.'2H4 ') tstr(1:4) = ' H42'

          if (tstr(2:4).eq.'OT1') tstr(2:4) = 'O  '
          if (tstr(2:4).eq.'OT2') tstr(2:4) = 'OXT'

          it = ichar(tstr(2:2))
          it4 = ichar(tstr(4:4))
          if (it.eq.67) then
             ianz = 6
          elseif (it.eq.78) then
             ianz = 7
          elseif (it.eq.79) then
             ianz = 8
          elseif (it.eq.80) then
             ianz = 15
          elseif (it.eq.83) then
             ianz = 16
          elseif (it.eq.72.or.it.eq.68) then
             ianz = 1
             if (it.eq.68) tstr(2:2) = 'H'
             if (tstr(2:4).eq.'HN ') tstr(2:4) = 'H  '
             if ((tstr(2:3).eq.'HA'.or.tstr(2:3).eq.'HB').and.
     &           (tstr(4:4).ne.' ')) tstr(4:4) = ' '
             if (tstr(2:3).eq.'HT'.and.(it4.ge.49.and.it4.le.51))
     &       then
                tstr(3:4) = '  '
                tstr(1:1) = tstr(4:4)
             endif
             ihashy = 1
             ish = 1
c Below convert 1H8 to H81, which not recognized by molden
c it probably served to read in force field H label variant
c             if (tstr(4:4).eq.' ') then
c                it1 = ichar(tstr(1:1))
c                it3 = ichar(tstr(3:3))
c                if ((it1.ge.49.and.it1.le.57).and.
c     &              (it3.ge.49.and.it3.le.57)) tstr(4:4) = tstr(1:1)
c             endif
          else
             ifnd = 0 
          endif
          
          if (ifnd.eq.0) then
            it = ichar(tstr(1:1))
            if (.not.(it.ge.65.and.it.le.90).and..not.
     &             (it.ge.97.and.it.le.122)) tstr(1:1) = ' '
            do j=1,99
             if (tstr(1:2).eq.tocapf(elemnt(j))) ianz = j
            end do
            if (ianz.eq.0.and.tstr(1:1).ne.' ') then
               tstr(1:1) = ' '
               do j=1,99
                  if (tstr(1:2).eq.tocapf(elemnt(j))) 
     &                ianz = j
               end do
            endif
            if (ianz.le.0) then
               if (tstr(2:2).eq.'D'.and.tstr(2:3).ne.'DY') then
                  ianz = 1
               else
                  write(iun3,*) 'Unclassified atom =',tstr
                  write(iun3,*) lstr
                  ianz = 99
               endif
            endif
          endif

      return
      end

c
c H-Bond geometry
c Thornton, MacDonald 1994
c
c D...A < 3.9 angstrom
c H...A < 2.5 angstrom
c D-H..A > 90 degrees
c
      subroutine opthdd(iupres,nupres,irsnr,
     &                  ncalf,icalf,coo,q,iconn,ianz,iresid,iamino,
     &                  ityp,ipdbt)
      implicit double precision (a-h,o-z)
      parameter (mxres=42)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      common /types/ iff
      integer*2 ityp,ipdbt
      character*3 aminos
      common /amino/ aminos(mxres)
      dimension copt(3),coo(3,*),q(*),iresid(*),ityp(*),iupres(*)
      dimension iamino(*),icalf(6,*),iconn(mxcon+1,*),ipdbt(*)
      dimension irsnr(*),ianz(*)

      print*," "
      print*,"Optimising O-H orientation of residues:"
      print*," "

      nupres = 0
      
      do i=1,iatoms
         it = ityp(i)
         ires = iresid(i)

         if (it.eq.64.or.it.eq.74.or.it.eq.137) then

            if (ires.gt.0) print*,aminos(iamino(ires)),
     &                     " ",irsnr(ires)
            if (it.eq.64) ang = 72.94d0
            if (it.eq.74) ang = 67.5d0
            if (it.eq.137) ang = 72.0d0

            call fndoh(i,ang,ires,copt,istat)
            if (istat.eq.1) then
               nupres = nupres + 1
               iupres(nupres) = iresid(i)
               do j=1,3
                  coo(j,i) = copt(j)
               end do
            endif
         endif
      end do


      ifftmp = iff

      call flpasn(nupres,iupres,ncalf,icalf,iamino,coo,q,iconn,
     &            ianz,iresid,irsnr,ityp)
      
      call hiseva(nupres,iupres,ncalf,icalf,iamino,coo,q,iconn,
     &            iresid,irsnr,ityp,ipdbt)

      iff = ifftmp

      call evwat

      return
      end

      logical function chk2flp(icl,ctmp,coo,iamino)
      implicit double precision (a-h,o-z)
      parameter (mxdbl=20)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (numcal=50000)
      parameter (mxliga=200)
      common /pmfrk/ pmlig(mxliga),ipmlig(mxliga),npmlig,
     &               ipmres(numcal),npmres
      dimension iamino(*),coo(3,*),idbl(mxdbl),cgd(3),ctmp(*)
      dimension ipdb(mxsym),ihpdb(mxhsym*3)

      toang = 0.52917706d0
      toang2 = toang*toang

      ndbl = 0
      chk2flp = .false.

      do i=1,npmres
         im = ipmres(i)
         if (im.gt.0) then
            ia = iamino(im)
            if (ia.eq.10.or.ia.eq.14) then
               if (ndbl.lt.mxdbl) then
                  ndbl = ndbl + 1
                  idbl(ndbl) = im
               endif
            endif
         endif
      end do

      icl = 0
      dmn = 10000.0d0
      dthres =  25.0d0/toang2

      do i=1,ndbl
         im = idbl(i)
         call getpdb(im,ipdb,ihpdb)
         ia = iamino(im)
         if (ia.eq.10) then
            icgd = ipdb(6)
         else
            icgd = ipdb(9)
         endif
         do j=1,3
            cgd(j) = coo(j,icgd)
         end do
         d2 = dist2(ctmp,cgd)
         if (d2.lt.dmn.and.d2.lt.dthres) then
            dmn = d2
            icl = im
         endif
      end do

      if (icl.ne.0) chk2flp = .true.

      return
      end

      subroutine sngflp(ires,iamino,coo)
      implicit double precision (a-h,o-z)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension coo(3,*),iamino(*)
      dimension chn1(3),chn2(3),cot(3),cnt(3)

      call getpdb(ires,ipdb,ihpdb)

      if (iamino(ires).eq.10) then
         iat1 = ipdb(21)
         iat2 = ipdb(6)
         iat3 = ipdb(29)
         ih1 =  ihpdb(25)
         ih2 =  ihpdb(26)
      else
         iat1 = ipdb(24)
         iat2 = ipdb(9)
         iat3 = ipdb(34)
         ih1 =  ihpdb(34)
         ih2 =  ihpdb(35)
      endif

      do j=1,3
         cnt(j)  = coo(j,iat1)
         cot(j)  = coo(j,iat3)
      end do

c swap O and N

      do j=1,3
         coo(j,iat1) = cot(j)
         coo(j,iat3) = cnt(j)
      end do

      call fliph(iat3,iat2,iat1,istat,
     &           0.950d0,120.0d0,180.0d0,chn1,coo)
      call fliph(iat3,iat2,iat1,istat,
     &           0.950d0,120.0d0,0.0d0,chn2,coo)

c New H coordinates

      do j=1,3
         coo(j,ih1) = chn1(j)
         coo(j,ih2) = chn2(j)
      end do

      return
      end

      subroutine flpasn(nupres,iupres,ncalf,icalf,iamino,coo,q,
     &                  iconn,ianz,iresid,irsnr,ityp)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      character*3 resn,resd
      integer*2 ityp
      logical chk2flp,dblflp,flp1,flp2,ck2het
      dimension ipdb(mxsym),ihpdb(mxhsym*3),iupres(*)
      dimension jpdb(mxsym),jhpdb(mxhsym*3),iflip(4)
      dimension irsnr(*),iamino(*),icalf(6,*),iresid(*),ityp(*)
      dimension coo(3,*),q(*),iconn(mxcon+1,*),ianz(*),
     &          chn1(3),chn2(3),cho1(3),cho2(3),cot(3),cnt(3),vec(3)
      dimension dhn1(3),dhn2(3),dho1(3),dho2(3),
     &          dot(3),dnt(3)

      do i=1,ncalf

         if (iamino(i).eq.10.or.iamino(i).eq.14) then

            call getpdb(i,ipdb,ihpdb)

            if (.not.ck2het(ipdb,coo,iconn,ianz,iresid)) then

               nupres = nupres + 1
               iupres(nupres) = i

c iat1 = ND2,NE2 
c iat2 = CG,CD
c iat3 = OD1,OE1

               flp1 = .false.
               flp2 = .false.

               if (iamino(i).eq.10) then
                  iat1 = ipdb(21)
                  iat2 = ipdb(6)
                  iat3 = ipdb(29)
                  ih1 =  ihpdb(25)
                  ih2 =  ihpdb(26)
                  resn = 'ASN'
               else
                  iat1 = ipdb(24)
                  iat2 = ipdb(9)
                  iat3 = ipdb(34)
                  ih1 =  ihpdb(34)
                  ih2 =  ihpdb(35)
                  resn = 'GLN'
               endif

               iflip(1) = iat1
               iflip(2) = iat3
               iflip(3) = ih1
               iflip(4) = ih2

               do j=1,3
                  cho1(j) = coo(j,ih1)
                  cho2(j) = coo(j,ih2)
                  cnt(j)  = coo(j,iat1)
                  cot(j)  = coo(j,iat3)
               end do

c update interacting residues

               do j=1,3
                  vec(j) = coo(j,iat2)
               end do

               call intres(0,i,vec,coo,iresid,irsnr,ncalf,icalf)

c check if other asn/gln is near: DOUBLE flip

               dblflp = chk2flp(idbl,vec,coo,iamino)

               if (dblflp.and.idbl.gt.i) then
                  call getpdb(idbl,jpdb,jhpdb)
                  if (iamino(idbl).eq.10) then
                     jat1 = jpdb(21)
                     jat2 = jpdb(6)
                     jat3 = jpdb(29)
                     jh1 =  jhpdb(25)
                     jh2 =  jhpdb(26)
                     resd = 'ASN'
                  else
                     jat1 = jpdb(24)
                     jat2 = jpdb(9)
                     jat3 = jpdb(34)
                     jh1 =  jhpdb(34)
                     jh2 =  jhpdb(35)
                     resd = 'GLN'
                  endif

                  do j=1,3
                     dho1(j) = coo(j,jh1)
                     dho2(j) = coo(j,jh2)
                     dnt(j)  = coo(j,jat1)
                     dot(j)  = coo(j,jat3)
                  end do

               endif

               if ((dblflp.and.idbl.gt.i).or..not.dblflp) then

c SINGLE flip res1

                  call evflp(iflip,score1,coo,q,iresid,ityp)
                  call evqvdw(i,ih1,ih2,scor,coo,q,ityp)
                  score1 = score1 + scor

                  if (dblflp) then
                     call evqvdw(idbl,jh1,jh2,scor,coo,q,ityp)
                     score1 = score1 + scor
                  endif

c swap O and N

                  do j=1,3
                     coo(j,iat1) = cot(j)
                     coo(j,iat3) = cnt(j)
                  end do

                  call fliph(iat3,iat2,iat1,istat,
     &                       0.950d0,120.0d0,180.0d0,chn1,coo)
                  call fliph(iat3,iat2,iat1,istat,
     &                       0.950d0,120.0d0,0.0d0,chn2,coo)

c New H coordinates

                  do j=1,3
                     coo(j,ih1) = chn1(j)
                     coo(j,ih2) = chn2(j)
                  end do

                  call evflp(iflip,score2,coo,q,iresid,ityp)
                  call evqvdw(i,ih1,ih2,scor,coo,q,ityp)
                  score2 = score2 + scor

                  if (dblflp) then
                     call evqvdw(idbl,jh1,jh2,scor,coo,q,ityp)
                     score2 = score2 + scor
                  endif

                  if (.not.dblflp) then
                     if (score2.lt.score1) flp1 = .true.
                  endif

               endif

               if (dblflp.and.idbl.gt.i) then

c double flip, flip res2

c swap O and N res 2

                  do j=1,3
                     coo(j,jat1) = dot(j)
                     coo(j,jat3) = dnt(j)
                  end do

                  call fliph(jat3,jat2,jat1,istat,
     &                       0.950d0,120.0d0,180.0d0,dhn1,coo)
                  call fliph(jat3,jat2,jat1,istat,
     &                       0.950d0,120.0d0,0.0d0,dhn2,coo)

c New H coordinates res2

                  do j=1,3
                     coo(j,jh1) = dhn1(j)
                     coo(j,jh2) = dhn2(j)
                  end do

                  call evflp(iflip,score3,coo,q,iresid,ityp)
                  call evqvdw(i,ih1,ih2,scor,coo,q,ityp)
                  score3 = score3 + scor
                  call evqvdw(idbl,jh1,jh2,scor,coo,q,ityp)
                  score3 = score3 + scor

                  do j=1,3
                     coo(j,iat1) = cnt(j)
                     coo(j,iat3) = cot(j)
                     coo(j,ih1)  = cho1(j)
                     coo(j,ih2)  = cho2(j)
                  end do

                  call evflp(iflip,score4,coo,q,iresid,ityp)
                  call evqvdw(i,ih1,ih2,scor,coo,q,ityp)
                  score4 = score4 + scor
                  call evqvdw(idbl,jh1,jh2,scor,coo,q,ityp)
                  score4 = score4 + scor

                  if (score2.lt.score1.and.score2.lt.score3.and.
     &                score2.lt.score4) flp1 = .true.

                  if (score4.lt.score1.and.score4.lt.score2.and.
     &                score4.lt.score3) flp2 = .true.


                  if (score3.lt.score1.and.score3.lt.score2.and.
     &                score3.lt.score4) then
                      flp1 = .true.
                      flp2 = .true.
                  endif

               endif

               if (flp1) then
c flip res 1
                  do j=1,3
                     coo(j,iat1) = cot(j)
                     coo(j,iat3) = cnt(j)
                     coo(j,ih1)  = chn1(j)
                     coo(j,ih2)  = chn2(j)
                  end do

                  print*,' ' 
                  print*,' Residue ',irsnr(i),' ',resn,' Flipped'
                  print*,' ' 
                  if (.not.dblflp) then
                     print*,' Score unflipped ',score1
                     print*,' Score   flipped ',score2
                     print*,' ' 
                  endif

               else

                  do j=1,3
                     coo(j,iat1) = cnt(j)
                     coo(j,iat3) = cot(j)
                     coo(j,ih1)  = cho1(j)
                     coo(j,ih2)  = cho2(j)
                  end do

               endif

               if (dblflp.and.idbl.gt.i) then

                  if (flp2) then
c flip res 2
                     do j=1,3
                        coo(j,jat1) = dot(j)
                        coo(j,jat3) = dnt(j)
                        coo(j,jh1)  = dhn1(j)
                        coo(j,jh2)  = dhn2(j)
                     end do

                     print*,' ' 
                     print*,' Residue ',irsnr(idbl),' ',resd,' Flipped'
                     print*,' ' 

                  else

                     do j=1,3
                        coo(j,jat1) = dnt(j)
                        coo(j,jat3) = dot(j)
                        coo(j,jh1)  = dho1(j)
                        coo(j,jh2)  = dho2(j)
                     end do

                  endif

                  print*,' Two flippable res in contact !'
                  print*,' ' 
                  print*,'   ',irsnr(i),'    ',irsnr(idbl),
     &                      '   score'
                  print*,' ------------------------'
                  print*,' unflip  unflip ',score1
                  print*,'   flip  unflip ',score2
                  print*,'   flip    flip ',score3
                  print*,' unflip    flip ',score4
                  print*,' ------------------------'
                  print*,' ' 

               endif

            endif

         endif

      end do

      return
      end

      subroutine fndod(itar,ang,ires,copt,istat,
     &                 coo,q,iconn,ityp,iresid,irsnr,ncalf,icalf)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxambc=49)
      integer*2 ityp
      integer ambvdt
      common /fcharg/ ambchg(mxamb),ambvw1(mxambc),ambvw2(mxambc),
     &                gfvdw(2,mxgff),ambvdt(mxamb),cysneg(9)
      dimension coo(3,*),iconn(mxcon+1,*),q(*),ityp(*),icalf(6,*)
      dimension iresid(*),irsnr(*)
      dimension a(3),b(3),c(3),copt(3),u(3),v(3),w(3),or(3),ar(3),vec(3)

      istat = 0
      todeg = 45.0d0 / datan(1.0d0)
      toang = 0.52917706d0
      sang = dsin(ang/todeg)
      cang = dcos(ang/todeg)
      oh = 1.002d0/toang
      soh = sang*oh
      coh = cang*oh

      do i=1,3
         a(i) = 0.0d0
         b(i) = 0.0d0
         vec(i) = coo(i,itar)
      end do

      i2 = 0
      i3 = 0
      in = iconn(1,itar)
      do i=1,in
         inn = iconn(1+i,itar)
         if (inn.gt.0) then
            i2 = inn
            do j=1,3
               a(j) = coo(j,i2)
            end do
         endif
      end do

      if (i2.ne.0) then
         in = iconn(1,i2)
         do i=1,in
            inn = iconn(1+i,i2)
            if (inn.gt.0.and.inn.ne.itar) then
               i3 = inn
               do j=1,3
                  b(j) = coo(j,i3)
               end do
            endif
         end do
      endif

      if (i2.eq.0.or.i3.eq.0) then
         print*,'problem with connections of Hydrogen to be optimised'
         return
      endif

      call intres(itar,ires,vec,coo,iresid,irsnr,ncalf,icalf)

      do i=1,3
         v(i) = a(i) - b(i)
      end do
      
      vl = vlen(v)

      call arbper(v,ar)

      do i=1,3
         v(i) = v(i)/vl
      end do

      call crprod(v,ar,u)
      call crprod(v,u,w)

      do i=1,3
         or(i) = a(i) + v(i)*coh
      end do
      
      qh = q(itar)
      iptp = ityp(itar)
      il = ambvdt(iptp)
      vdr = ambvw1(il)
      vde = ambvw2(il)

      scoro = 10000.0d0

      do i=0,360,2
         theta = dble(i) / todeg
         st = dsin(theta)
         ct = dcos(theta)

         do j=1,3 
            c(j) = soh*ct*u(j) + soh*st*w(j) + or(j)
         end do

         call evpos(c,qh,vdr,vde,itar,i2,i3,scor,coo,q,iresid,ityp)

         if (scor.lt.scoro) then
            scoro = scor
            do j=1,3 
               copt(j) = c(j)
            end do
         endif

      end do

      istat = 1

      return
      end

      subroutine intres(itar,ires,vec,coo,iresid,irsnr,ncalf,icalf)
      implicit double precision (a-h,o-z)
      parameter (mxliga=200)
      parameter (numcal=50000)
c (mis)use pmf work array to store close residues
      common /pmfrk/ pmlig(mxliga),ipmlig(mxliga),npmlig,
     &               ipmres(numcal),npmres
      common /athlp/ iatoms, mxnat
      dimension coo(3,*),icalf(6,*),vec(3),iresid(*),irsnr(*)


      erad = 19.0d0*19.0d0

c max radius amino acid estim. 10 Angs

      npmres = 0
      do i=1,ncalf

         if ((itar.eq.0.and.ires.ne.i).or.itar.ne.0) then

            if (itar.eq.0) then
               d2 = dist2(coo(1,icalf(1,i)),vec)
            else
               d2 = dist2(coo(1,icalf(1,i)),coo(1,itar))
            endif

            if (d2.lt.erad) then
               if (npmres.lt.numcal) then
                  npmres = npmres + 1
                  ipmres(npmres) = i
               endif
            endif

         endif
      end do

      do i=1,iatoms

         if (iresid(i).lt.-3) then

            if (itar.eq.0) then
               d2 = dist2(coo(1,i),vec)
            else
               d2 = dist2(coo(1,i),coo(1,itar))
            endif

            if (d2.lt.100.0d0) then

               ido = 1
               do j=1,npmres
                  if (ipmres(j).eq.iresid(i)) ido = 0
               end do

               if (ido.eq.1) then
                  if (npmres.lt.numcal) then
                     npmres = npmres + 1
                     ipmres(npmres) = iresid(i)
                  endif
               endif
            endif

        endif
      end do

      return
      end

      subroutine evqvdw(ires,i2,i3,scor,coo,q,ityp)

c score hydrogens of NH2 group asn/gln against other H same residue
c vdw part only

      implicit double precision (a-h,o-z)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxambc=49)
      integer*2 ityp,iptp,ihtp
      integer ambvdt
      common /fcharg/ ambchg(mxamb),ambvw1(mxambc),ambvw2(mxambc),
     &                gfvdw(2,mxgff),ambvdt(mxamb),cysneg(9)
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension coo(3,*),q(*),ityp(*)

      econv = 332.05382d0
      toang = 0.52917706d0
      scor = 0.0d0
      qh = 0.4223d0

      call getpdb(ires,ipdb,ihpdb)

      do j=1,mxhsym*3

         ipat = ihpdb(j)
         if (ipat.ne.0.and.ipat.ne.i2.and.ipat.ne.i3) then
            d2 = dist2(coo(1,ipat),coo(1,i2))
            d2 = dsqrt(d2)*toang
            eqt = (econv * q(ipat) * qh / d2)
            iptp = ityp(ipat)
            ihtp = ityp(i2)

            if (iptp.gt.0.and.ihtp.gt.0) then
               il = ambvdt(ihtp)
               vdwr1 = ambvw1(il)
               vdwe1 = ambvw2(il)
               il = ambvdt(iptp)
               vdwr2 = ambvw1(il)
               vdwe2 = ambvw2(il)
               rsum = vdwr1 + vdwr2
               p6 = (rsum/d2)**6.0d0
               p12 = p6 * p6
               epsm = dsqrt(vdwe1 * vdwe2)
               scor    = scor + eqt + epsm * (p12 - 2.0d0*p6)
            endif

            d2 = dist2(coo(1,ipat),coo(1,i3))
            d2 = dsqrt(d2)*toang
            eqt = (econv * q(ipat) * qh / d2)
            ihtp = ityp(i3)

            if (iptp.gt.0.and.ihtp.gt.0) then
               il = ambvdt(ihtp)
               vdwr1 = ambvw1(il)
               vdwe1 = ambvw2(il)
               il = ambvdt(iptp)
               vdwr2 = ambvw1(il)
               vdwe2 = ambvw2(il)
               rsum = vdwr1 + vdwr2
               p6 = (rsum/d2)**6.0d0
               p12 = p6 * p6
               epsm = dsqrt(vdwe1 * vdwe2)
               scor    = scor + eqt + epsm * (p12 - 2.0d0*p6)
            endif
         endif

      end do

      return
      end 

      subroutine evpos(vec,qh,vdwr1,vdwe1,itar,i2,i3,scor,coo,q,
     &                 iresid,ityp)
      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxambc=49)
      parameter (mxliga=200)
      common /pmfrk/ pmlig(mxliga),ipmlig(mxliga),npmlig,
     &               ipmres(numcal),npmres
      integer*2 ityp,iptp
      integer ambvdt
      common /fcharg/ ambchg(mxamb),ambvw1(mxambc),ambvw2(mxambc),
     &                gfvdw(2,mxgff),ambvdt(mxamb),cysneg(9)
      common /athlp/ iatoms, mxnat
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension coo(3,*),q(*),ityp(*),iresid(*),vec(*)

      toang = 0.52917706d0
      econv = 332.05382d0

      eq = 0.0d0
      ev = 0.0d0

      eqmn = 10000.0d0
      eqmx = -10000.0d0
      iqmn = 0
      iqmx = 0
      irn = 0
      irx = 0

      do i=1,npmres

         if (ipmres(i).gt.0) then
            call getpdb(ipmres(i),ipdb,ihpdb)

            do j=1,mxsym

               ip = ipdb(j)
               if (ip.ne.0.and.ip.ne.i2.and.ip.ne.i3) then

                  ipat = ip
                  iptp = ityp(ipat)

                  if (iptp.gt.0) then
                     d2 = dist2(coo(1,ipat),vec)
                     d2 = dsqrt(d2)*toang
                     eqt = (econv * q(ipat) * qh / d2)
                     if (eqt.lt.eqmn) then
                         eqmn = eqt
                         iqmn = ipat
                         irn = ipmres(i)
                     endif
                     if (eqt.gt.eqmx) then
                         eqmx = eqt
                         iqmx = ipat
                         irx = ipmres(i)
                     endif
                     eq = eq + eqt

                     il = ambvdt(iptp)
                     vdwr2 = ambvw1(il)
                     vdwe2 = ambvw2(il)

                     rsum = vdwr1 + vdwr2
                     p6 = (rsum/d2)**6.0d0
                     p12 = p6 * p6
                     epsm = dsqrt(vdwe1 * vdwe2)
                     ev    = ev + epsm * (p12 - 2.0d0*p6)

                  endif

               endif
            end do

            do j=1,mxhsym*3

               if (ihpdb(j).ne.0.and.ihpdb(j).ne.itar) then
                  ipat = ihpdb(j)
                  d2 = dist2(coo(1,ipat),vec)
                  d2 = dsqrt(d2)*toang
                  eqt = (econv * q(ipat) * qh / d2)
                  if (eqt.lt.eqmn) then
                      eqmn = eqt
                      iqmn = ipat
                      irn = ipmres(i)
                  endif
                  if (eqt.gt.eqmx) then
                      eqmx = eqt
                      iqmx = ipat
                      irx = ipmres(i)
                  endif
                  eq = eq + eqt

                  iptp = ityp(ipat)
                  if (iptp.gt.0) then
                     il = ambvdt(iptp)
                     vdwr2 = ambvw1(il)
                     vdwe2 = ambvw2(il)

                     rsum = vdwr1 + vdwr2
                     p6 = (rsum/d2)**6.0d0
                     p12 = p6 * p6
                     epsm = dsqrt(vdwe1 * vdwe2)
                     ev    = ev + epsm * (p12 - 2.0d0*p6)
                  endif

               endif
            end do

         else

            do j=1,iatoms
               if (iresid(j).eq.ipmres(i)) then
                  iptp = ityp(j)
                  il = abs(iptp)

                  d2 = dist2(coo(1,j),vec)
                  d2 = dsqrt(d2)*toang
                  eqt = (econv * q(j) * qh / d2)
                  eq = eq + eqt

                  vdwr2 = gfvdw(1,il)
                  vdwe2 = gfvdw(2,il)

                  rsum = vdwr1 + vdwr2
                  p6 = (rsum/d2)**6.0d0
                  p12 = p6 * p6
                  epsm = dsqrt(vdwe1 * vdwe2)
                  ev    = ev + epsm * (p12 - 2.0d0*p6)

               endif
            end do

         endif

      end do

      scor = eq + ev

      return
      end

      logical function ck2het(ipdb,coo,iconn,ianz,iresid)
      implicit double precision (a-h,o-z)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxcon=10)
      logical ismet
      dimension ipdb(mxsym)
      dimension coo(3,*),iresid(*),ianz(*),iconn(mxcon+1,*)

      ck2het = .false.

      do i=1,mxsym
         ii = ipdb(i)
         if (ii.ne.0) then
            do j=1,iconn(1,ii)
               jj = iconn(1+j,ii)
               if (jj.gt.0) then
                  if (iresid(jj).lt.-3) then
                     if (.not.ismet(ianz(jj))) ck2het = .true.
                  endif
               endif
            end do
         endif
      end do

      return
      end

      subroutine evres(ires,ipdb,ihpdb,iretyp,score,
     &                 coo,q,iresid,ityp,ncalf,icalf)
      implicit double precision (a-h,o-z)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxambc=49)
      integer*2 ityp
      integer ambvdt
      common /fcharg/ ambchg(mxamb),ambvw1(mxambc),ambvw2(mxambc),
     &                gfvdw(2,mxgff),ambvdt(mxamb),cysneg(9)
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension coo(3,*),q(*),ityp(*),iresid(*),icalf(6,*)

      i0 = 0
      score = 0.0d0
      scorp = 0.0d0
      scorh = 0.0d0

      if (iretyp.eq.1) then
         call typamb(ipdb,17,ihpdb,1)
      endif

      do j=1,mxsym
         i2 = 0
         i3 = 0
         if (j.eq.1) then
            if (ires.gt.1) then
               i2 = icalf(3,ires-1)
               i3 = icalf(1,ires-1)
            endif
         elseif (j.eq.2) then
            if (ires.gt.1) then
               i2 = icalf(3,ires-1)
            endif
            if (ires.lt.ncalf) then
               i3 = icalf(2,ires+1)
            endif
         elseif (j.eq.3) then
            if (ires.lt.ncalf) then
               i2 = icalf(1,ires+1)
               i3 = icalf(2,ires+1)
            endif
         endif
         if (ipdb(j).ne.0) then
            ipat  = ipdb(j)
            iptp  = ityp(ipat)
            if (iptp.gt.0) then
               il    = ambvdt(iptp)
               vdr = ambvw1(il)
               vde = ambvw2(il)
               qh    = q(ipat)
               call evpos(coo(1,ipat),qh,
     &                 vdr,vde,i0,i2,i3,scor1,coo,q,iresid,ityp)
               scorp = scorp + scor1
            endif
         endif
      end do

      do j=1,mxhsym*3
         ipat  = ihpdb(j)
         if (ipat.ne.0) then
            iptp  = ityp(ipat)
            if (iptp.gt.0) then
               il    = ambvdt(iptp)
               vdr = ambvw1(il)
               vde = ambvw2(il)
               qh    = q(ipat)
               call evpos(coo(1,ipat),qh,
     &                    vdr,vde,i0,i0,i0,scor2,coo,q,iresid,ityp)
               scorh = scorh + scor2
            endif
         endif
      end do
    
      score = scorp + scorh

      return
      end

      subroutine evflp(iflip,score,coo,q,iresid,ityp)
      implicit double precision (a-h,o-z)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxambc=49)
      integer*2 ityp
      integer ambvdt
      common /fcharg/ ambchg(mxamb),ambvw1(mxambc),ambvw2(mxambc),
     &                gfvdw(2,mxgff),ambvdt(mxamb),cysneg(9)
      dimension iflip(4)
      dimension coo(3,*),q(*),ityp(*),iresid(*)

      i0 = 0
      score = 0.0d0

      do j=1,4
         jj = iflip(j)
         if (jj.ne.0) then
            iptp  = ityp(jj)
            if (iptp.gt.0) then
               il    = ambvdt(iptp)
               vdr = ambvw1(il)
               vde = ambvw2(il)
               qh    = q(jj)
               call evpos(coo(1,jj),qh,
     &                 vdr,vde,i0,i0,i0,scor1,coo,q,iresid,ityp)
               score = score + scor1
            endif
         endif
      end do

      return
      end

      subroutine hiseva(nupres,iupres,
     &                  ncalf,icalf,iamino,coo,q,
     &                  iconn,iresid,irsnr,ityp,ipdbt)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxambc=49)
      logical chkat
      integer*2 ityp,ipdbt
      common /athlp/ iatoms, mxnat
      integer ambvdt
      common /fcharg/ ambchg(mxamb),ambvw1(mxambc),ambvw2(mxambc),
     &                gfvdw(2,mxgff),ambvdt(mxamb),cysneg(9)
      dimension ipdb(mxsym),ihpdb(mxhsym*3),iconn(mxcon+1,*),iupres(*)
      dimension coo(3,*),q(*),iamino(*),icalf(6,*),ityp(*),ipdbt(*)
      dimension iresid(*),vhd(3),vhe(3),vec3(3),icnn(mxcon),irsnr(*)

c  1  HIS+ (HIP)
c  2  HISD (HID)
c  3  HISE (HIE)

      print*,' '
      print*,'Evaluation of Histidines:'
      print*,' '

      do i=1,ncalf

         if (iamino(i).eq.17) then

            nupres = nupres + 1
            iupres(nupres) = i

c when chkat is false, you should not take average
c vhd and vhe !!

            call getpdb(i,ipdb,ihpdb)
            if (.not.chkat(ipdb(6),ipdb(13),ipdb(20),
     &                1.020d0,126.0d0,180.0d0,vhd,coo)) then
                call vclr(vhd,1,3)
            endif
            if (.not.chkat(ipdb(13),ipdb(11),ipdb(24),
     &                1.020d0,126.0d0,180.0d0,vhe,coo)) then
                call vclr(vhe,1,3)
            endif

            do j=1,3
               vec3(j) = (vhd(j) + vhe(j)) / 2.0d0
            end do

            call intres(0,i,vec3,coo,iresid,irsnr,ncalf,icalf)

c save original HD1 pointer

            ihpdb22 = ihpdb(22)
            ihpdb34 = ihpdb(34)

c HIP
            ihpdb(22) = iatoms + 1
            ihpdb(34) = iatoms + 2
            do j=1,3
               coo(j,iatoms+1) = vhd(j)
               coo(j,iatoms+2) = vhe(j)
            end do
            call evres(i,ipdb,ihpdb,1,scorep,coo,q,iresid,ityp,
     &                 ncalf,icalf)

c HID
            ihpdb(22) = iatoms + 1
            ihpdb(34) = 0
            call evres(i,ipdb,ihpdb,1,scored,coo,q,iresid,ityp,
     &                 ncalf,icalf)

c HIE
            ihpdb(22) = 0
            ihpdb(34) = iatoms + 2
            call evres(i,ipdb,ihpdb,1,scoree,coo,q,iresid,ityp,
     &                 ncalf,icalf)

            write(*,'(a,i4,a,f9.3,a,f9.3,a,f9.3)') ' his residue ',
     &        irsnr(i),' score + ',scorep,' D ',scored,' E ',scoree

c decide HID or HIE: 
c assumes HID as default

            ihis = 2
            if (scoree.lt.scored) ihis = 3

            if (ihis.eq.1) then

                ihpdb(22) = ihpdb22
                ihpdb(34) = iatoms + 1

                do j=1,3
                   coo(j,ihpdb22)  = vhd(j)
                   coo(j,iatoms+1) = vhe(j)
                end do
   
                iconn(1,iatoms+1) = 1
                iconn(2,iatoms+1) = ipdb(24)

                iconn(1,ipdb(24)) = iconn(1,ipdb(24)) + 1
                ncon = iconn(1,ipdb(24))
                iconn(1+ncon,ipdb(24)) = iatoms + 1

c dont forget iatclr etc. (ipdbt also)

                iatoms = iatoms + 1


            else if (ihis.eq.2) then

                if (ihpdb22.gt.0.and.ihpdb22.le.mxnat) then
                   ihpdb(22) = ihpdb22
                   ihpdb(34) = 0
                   ipdbt(ihpdb(22)) = 22 
                   do j=1,3
                      coo(j,ihpdb22) = vhd(j)
                   end do
                endif

            else if(ihis.eq.3) then

c connect new HE2 to NE2


                if (ihpdb22.gt.0.and.ihpdb22.le.mxnat.and.
     &             ipdb(24).ne.0) then

                   ihpdb(22) = 0
                   if (ihpdb34.ne.0) then
                      ihpdb(34) = ihpdb34
                      ipdbt(ihpdb(34)) = 34 
                   else
                      ihpdb(34) = ihpdb22
                      ipdbt(ihpdb(34)) = 34 
   
                      do j=1,3
                         coo(j,ihpdb22) = vhe(j)
                      end do
   
                      iconn(2,ihpdb22) = ipdb(24)
   
                      iconn(1,ipdb(24)) = iconn(1,ipdb(24)) + 1
                      ncon = iconn(1,ipdb(24))
                      iconn(1+ncon,ipdb(24)) = ihpdb22
   
c disconnect HD1 from ND1
   
                      ic = 0
                      ncon = iconn(1,ipdb(20))

                      do j=1,ncon
                         icon = iconn(1+j,ipdb(20))
                         if (icon.ne.ihpdb22) then
                            ic = ic + 1
                            icnn(ic) = icon
                         endif
                      end do

                      iconn(1,ipdb(20)) = ic

                      do j=1,ic
                         iconn(1+j,ipdb(20)) = icnn(j)
                      end do

                   endif

                endif

            endif

            call typamb(ipdb,17,ihpdb,1)

         endif
      end do

      print*,' '
      print*,'Molden only chooses between HID and HIE.'
      print*,'Only the score of HIP is calculated'
      print*,' '

      return
      end

      subroutine sethis(ires,ihis,
     &                  coo,q,iresid,iatclr,iaton,iconn,ianz,
     &                  ncalf,icalf,ityp,ipdbt)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      logical chkat
      integer*2 ipdbt,ityp
      common /athlp/ iatoms, mxnat
      dimension ipdb(mxsym),ihpdb(mxhsym*3),iconn(mxcon+1,*)
      dimension coo(3,*),q(*),iresid(*),ianz(*)
      dimension iatclr(*),iaton(*),ityp(*),ipdbt(*),icalf(6,*)
      dimension vhd(3),vhe(3),vec3(3),icnn(mxcon)

c  1  HIS+ (HIP)
c  2  HISD (HID)
c  3  HISE (HIE)

      call getpdb(ires,ipdb,ihpdb)

      if (.not.chkat(ipdb(6),ipdb(13),ipdb(20),
     &                1.020d0,126.0d0,180.0d0,vhd,coo)) then
          call vclr(vhd,1,3)
      endif

      if (.not.chkat(ipdb(13),ipdb(11),ipdb(24),
     &                1.020d0,126.0d0,180.0d0,vhe,coo)) then
          call vclr(vhe,1,3)
      endif

      do j=1,3
         vec3(j) = (vhd(j) + vhe(j)) / 2.0d0
      end do


c save original HD1 pointer

      ihpdb22 = ihpdb(22)
      ihpdb34 = ihpdb(34)

      if (ihis.eq.1) then

c HIP
          idel = 0
          ipnt1 = ihpdb22
          ipnt2 = ihpdb34
          if (ipnt1.eq.0) ipnt1 = iatoms + 1
          if (ipnt2.eq.0) ipnt2 = iatoms + 1

          ihpdb(22) = ipnt1
          ihpdb(34) = ipnt2

          do j=1,3
             coo(j,ipnt1) = vhd(j)
             coo(j,ipnt2) = vhe(j)
          end do
   
          iresid(ipnt1) = ires
          iresid(ipnt2) = ires
          iaton(ipnt1) = 1
          iaton(ipnt2) = 1
          ianz(ipnt1) = 1
          ianz(ipnt2) = 1
          ipdbt(ipnt1) = 22 
          ipdbt(ipnt2) = 34 
          iatclr(ipnt1) = iatclr(ipdb(20))
          iatclr(ipnt2) = iatclr(ipdb(24))

          if (ihpdb34.eq.0) then

             iconn(1,ipnt2) = 1
             iconn(2,ipnt2) = ipdb(24)

             iconn(1,ipdb(24)) = iconn(1,ipdb(24)) + 1
             ncon = iconn(1,ipdb(24))
             iconn(1+ncon,ipdb(24)) = ipnt2

          endif

          if (ihpdb22.eq.0) then

             iconn(1,ipnt1) = 1
             iconn(2,ipnt1) = ipdb(20)

             iconn(1,ipdb(20)) = iconn(1,ipdb(20)) + 1
             ncon = iconn(1,ipdb(20))
             iconn(1+ncon,ipdb(20)) = ipnt1

          endif

          if (ihpdb22.eq.0.or.ihpdb34.eq.0) then
             iatoms = iatoms + 1
          endif


      else if (ihis.eq.2) then

c HID

          idel = ihpdb34
          ihpdb(34) = 0

          if (ihpdb22.ne.0) then
             ihpdb(22) = ihpdb22
             ipdbt(ihpdb(22)) = 22 
          endif

          if (ihpdb22.eq.0) then

             ihpdb(22) = ihpdb34
             ipdbt(ihpdb(22)) = 22 

             do j=1,3
                coo(j,ihpdb34) = vhd(j)
             end do

c connect new HD1 to ND1

             iconn(2,ihpdb34) = ipdb(20)

             iconn(1,ipdb(20)) = iconn(1,ipdb(20)) + 1
             ncon = iconn(1,ipdb(20))
             iconn(1+ncon,ipdb(20)) = ihpdb34

          endif

          if (ihpdb34.ne.0) then
c disconnect HE2 from NE2

             ic = 0
             ncon = iconn(1,ipdb(24))

             do j=1,ncon
                icon = iconn(1+j,ipdb(24))
                if (icon.ne.ihpdb34) then
                   ic = ic + 1
                   icnn(ic) = icon
                endif
             end do

             iconn(1,ipdb(24)) = ic

             do j=1,ic
                iconn(1+j,ipdb(24)) = icnn(j)
             end do

          endif


      else if(ihis.eq.3) then

c HIE

c connect new HE2 to NE2

          idel = ihpdb22
          ihpdb(22) = 0
          if (ihpdb34.ne.0) then
             ihpdb(34) = ihpdb34
             ipdbt(ihpdb(34)) = 34 
          endif

          if (ihpdb34.eq.0) then
             ihpdb(34) = ihpdb22
             ipdbt(ihpdb(34)) = 34 

             do j=1,3
                coo(j,ihpdb22) = vhe(j)
             end do

             iconn(2,ihpdb22) = ipdb(24)

             iconn(1,ipdb(24)) = iconn(1,ipdb(24)) + 1
             ncon = iconn(1,ipdb(24))
             iconn(1+ncon,ipdb(24)) = ihpdb22
          endif

c disconnect HD1 from ND1

          if (ihpdb22.ne.0) then
             ic = 0
             ncon = iconn(1,ipdb(20))

             do j=1,ncon
                icon = iconn(1+j,ipdb(20))
                if (icon.ne.ihpdb22) then
                   ic = ic + 1
                   icnn(ic) = icon
                endif
             end do

             iconn(1,ipdb(20)) = ic

             do j=1,ic
                iconn(1+j,ipdb(20)) = icnn(j)
             end do

          endif

      endif

      call typamb(ipdb,17,ihpdb,1)

      if (ihpdb22.ne.0.and.ihpdb34.ne.0.and.idel.ne.0) then
          call atdel(idel,coo,q,iresid,iatclr,iaton,iconn,ianz,
     &               ncalf,icalf,ityp,ipdbt)
      endif

c two problems:
c going from hid,hie to hip, a new atom is added. This is on top of the list
c above the fake atoms for the secundary structure calculation.
c best solution would be to insert it in the residue range.
c
c going the other way, we have a atom position we nolonger should use.
c However it is still there. Decreasing iatoms, is not enough. We have to move
c everything one up. There are such routines in xwin.c i think.

      return
      end

      subroutine evwad(coo,q,iresid,irsnr,iatclr,iaton,iconn,ianz,
     &                 ncalf,icalf,ityp,ipdbt)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      integer*2 ityp,ipdbt
      common /athlp/ iatoms, mxnat
      parameter (mxliga=200)
      parameter (numcal=50000)
c (mis)use pmf work array to store close residues
      common /pmfrk/ pmlig(mxliga),ipmlig(mxliga),npmlig,
     &               ipmres(numcal),npmres
      logical ismet
      dimension coo(3,*),q(*),iresid(*),irsnr(*),ianz(*),
     &          iconn(mxcon+1,*)
      dimension iatclr(*),iaton(*),ityp(*),ipdbt(*),icalf(6,*)
      dimension idona(2),rdon(2),vec(3),copt(3)
      

      toang = 0.52917706d0
      toang2 = toang*toang
      d2hb = 3.9d0*3.9d0/toang2
      d2ha = 2.5d0*2.5d0/toang2

      print*," "
      print*,"Optimising waters:"
      print*,"Deleting waters not in contact with protein,"
      print*,"Adding hydrogens to the ones that are."
      print*," "

      natoms = iatoms

      i = 1

      do while(i.le.iatoms)

         if (iresid(i).lt.-3.and.ianz(i).eq.8) then

            ic = 0
            do j=1,iconn(1,i)
               icon = iconn(1+j,i)
               if (icon.gt.0) then
                  if (.not.ismet(ianz(icon))) ic = ic + 1
               endif
            end do

c           only if O is not connected, 
c           we know this is water without H's

            if (ic.eq.0) then

               idon = 0
               iacc = 0

               do j=1,natoms
                  ia = ianz(j)
                  if (iresid(j).gt.0.and.(ia.eq.7.or.ia.eq.8)) then
                     d2 = dist2(coo(1,i),coo(1,j))

                     if (d2.le.d2hb) then

c found O,N closer than 3.9 angs from water O

                        iacct = 0
                        do k=1,iconn(1,j)
                           icon = iconn(1+k,j)
                           if (icon.gt.0.and.ianz(iabs(icon)).eq.1) then
                              d2h = dist2(coo(1,i),coo(1,icon))

c the water O is an acceptor: iacc = 1

                              if (d2h.le.d2ha) iacct = 1

                           endif
                        end do


                        if (iacct.eq.1) then
                           iacc = 1
                        else
                           if (idon.eq.0) then
                              idon = 1
                              idona(1) = j
                              rdon(1) = d2
                           else if (idon.eq.1) then
                              idon = 2
                              idont = idona(1)
                              rdont = rdon(1)
                              if (d2.lt.rdon(1)) then
                                 idona(1) = j
                                 rdon(1) = d2
                                 idona(2) = idont
                                 rdon(2) = rdont
                              else
                                 idona(2) = j
                                 rdon(2) = d2
                              endif
                           else if (idon.eq.2) then
                              idont = idona(1)
                              rdont = rdon(1)
                              if (d2.lt.rdon(1)) then
                                 idona(1) = j
                                 rdon(1) = d2
                                 idona(2) = idont
                                 rdon(2) = rdont
                              else if (d2.lt.rdon(2)) then
                                 idona(2) = j
                                 rdon(2) = d2
                              endif
                           endif
                        endif

                     endif

                  endif
               end do

c we require Owat as: 1 acceptor and 1 donor / 2 donors or more

               if (iatoms+2.gt.mxnat) then
                  print*,'No room to add new hydrogens !'
                  return
               endif

               if (iacc.ge.2.or.idon.eq.2.or.
     &             (iacc.eq.1.and.idon.eq.1)) then

c              calculate vector from donor to Owat

                  do k=1,3
                     vec(k) = coo(k,i) - coo(k,idona(1))
                  end do

                  vl = vlen(vec)

                  do k=1,3
                     coo(k,iatoms+1) = coo(k,i) - 
     &                                 0.957*vec(k)/(vl*toang)
                  end do

                  iconn(1,i) = 1
                  iconn(2,i) = iatoms + 1

                  iconn(1,iatoms+1) = 1
                  iconn(2,iatoms+1) = i

                  ianz(iatoms+1) = 1
                  iatclr(iatoms+1) = 1
                  iaton(iatoms+1) = iaton(i)
                  q(iatoms+1) =  0.4710d0
                  ityp(iatoms+1) = 650
                  iresid(iatoms+1) = iresid(i)

                  iatoms = iatoms + 1

                  q(i) = -0.8340d0
                  ityp(i) = 649

                  ang = 75.48
                  ires = iresid(i)

                  iconn(1,i) = 2
                  iconn(3,i) = iatoms + 1
 
                  iconn(1,iatoms+1) = 1
                  iconn(2,iatoms+1) = i

                  ianz(iatoms+1) = 1
                  iatclr(iatoms+1) = 1
                  iaton(iatoms+1) = iaton(i)
                  q(iatoms+1) =  0.4710d0
                  ityp(iatoms+1) = 650
                  iresid(iatoms+1) = iresid(i)

c set temporary coordinates to water  O 

                  do j=1,3
                        coo(j,iatoms+1) = coo(j,i)
                  end do

                  call fndod(iatoms+1,ang,ires,copt,istat,
     &                       coo,q,iconn,ityp,iresid,irsnr,ncalf,icalf)

                  if (istat.eq.1) then

                     do j=1,3
                        coo(j,iatoms+1) = copt(j)
                     end do

                     iatoms = iatoms + 1

                  endif

                  i = i + 1

               else

c throw away this water Oxygen

                  call atdel(i,coo,q,iresid,iatclr,iaton,iconn,ianz,
     &                       ncalf,icalf,ityp,ipdbt)

               endif

            else
               i = i + 1
            endif

         else
            i = i + 1
         endif

      end do

      return
      end

      subroutine atdel(iat,coo,q,iresid,iatclr,iaton,iconn,ianz,
     &                 ncalf,icalf,ityp,ipdbt)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      integer*2 ityp,ipdbt
      common /athlp/ iatoms, mxnat
      parameter (numcal=50000)
      dimension coo(3,*),q(*),iresid(*),ianz(*),
     &          iconn(mxcon+1,*)
      dimension iatclr(*),iaton(*),ityp(*),ipdbt(*),icalf(6,*)

      do j=1,iatoms-1

         if (j.lt.iat) then

            do k=1,iconn(1,j)
               icon = iconn(1+k,j)
               if (icon.lt.0) then
                  if (iabs(icon).lt.iat) then
                     iconn(1+k,j) = icon
                  else
                     iconn(1+k,j) = icon + 1
                  endif
               else
                  if (icon.lt.iat) then
                     iconn(1+k,j) = icon
                  else
                     iconn(1+k,j) = icon - 1
                  endif
               endif
            end do

         else

            j1 = j + 1

            do k=1,3
               coo(k,j) = coo(k,j1)
            end do

            iconn(1,j) = iconn(1,j1)

            do k=1,iconn(1,j1)
               icon = iconn(1+k,j1)
               if (icon.lt.0) then
                  if (iabs(icon).lt.iat) then
                     iconn(1+k,j) = icon
                  else
                     iconn(1+k,j) = icon + 1
                  endif
               else
                  if (icon.lt.iat) then
                     iconn(1+k,j) = icon
                  else
                     iconn(1+k,j) = icon - 1
                  endif
               endif
            end do
   
            ianz(j)   = ianz(j1)
            iatclr(j) = iatclr(j1)
            iaton(j)  = iaton(j1)
            q(j)      = q(j1)
            ityp(j)   = ityp(j1)
            ipdbt(j)  = ipdbt(j1)
            iresid(j) = iresid(j1)
                     
         endif
      end do

      do k=1,ncalf
         if (icalf(1,k).gt.iat) 
     &       icalf(1,k) = icalf(1,k) - 1
         if (icalf(4,k).gt.iat) 
     &       icalf(4,k) = icalf(4,k) - 1
      end do

      iatoms = iatoms - 1

      return
      end

      subroutine quwad(nwat,iresid,iconn,ianz)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      dimension iresid(*),ianz(*),iconn(mxcon+1,*)
      
      nwat = 0

      do i=1,iatoms

         if (iresid(i).lt.-3.and.ianz(i).eq.8) then

            ic = 0
            do j=1,iconn(1,i)
               icon = iconn(1+j,i)
               if (icon.gt.0) ic = ic + 1
            end do

            if (ic.eq.0) nwat = nwat + 1

         endif

      end do

      return
      end

