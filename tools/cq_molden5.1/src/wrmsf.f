      subroutine wrmsfd(iun,
     &                  coo,qat,ianz,iaton,iatclr,iresid,iconn,
     &                  lring,ityp,nat,nspg,icel,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxel=100)
      parameter (numatm=2000)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (mxsg=238)
      common /athlp/  iatoms, mxnat
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      integer*2 ityp
      common /types/  iff
      character*7 spnam
      common /spgrnm/ spnam(mxsg)
      character*2 elemnt
      common /elem/   elemnt(mxel)
      character*200 header
      character*100 fname
      character*80  col
      character*10  version
      character*6   mklab
      integer*2 nbc(numat1),ibc(numat1,mxcon),ictyp
      integer*2 i2dum,ibc2,iorder(numat1*mxcon)
      real rdum
      logical ochg
      dimension ibc2(2,numat1*mxcon),rr(3,3)
      dimension coo(3,*),ianz(*),iaton(*),iatclr(*),iresid(*),
     &          iconn(mxcon+1,*),qat(*),ityp(*),lring(*),
     &          ir(3,3,192),it(3,192)

      istat = 1
      toang = 0.52917706d0
      torad = datan(1.0d0) / 45.0d0
      idochg = 0

      fname = 'kull'
      natoms = iatoms

      if (icel.eq.1) then

         if (ochg(idum,ianz)) then
             idochg = 1
         else
             idochg = 0
         endif

         natoms = nat
         nstor = mxnat-natoms

         do i=1,natoms
            do j=1,3
               coo(j,i) = coo(j,nstor+i)
            end do
            ianz(i) = ianz(nstor+i)
            do j=1,iconn(1,nstor+i)+1
               iconn(j,i) = iconn(j,nstor+i)
            end do
            iatclr(i) = iatclr(nstor+i)
         end do

      else
         if (ihasq.eq.1) idochg = 1
      endif

      do i=1,natoms
         lring(i) = 0
         iaton(i) = 2
      end do

      ns = 1
      nr = 1
      version = 'QUANTAR3.3'
      header = 'Molden generated MSF'
      write(iun) ns,nr,natoms,version,header
      col = 'END'
      write(iun) col

c
c write segment information
c
      i2dum = 1
      write(iun) 'NEWS'
      write(iun) i2dum
      write(iun) i2dum

c
c write residue information
c
      i2dum = 1
      write(iun) 'RES1  '
      write(iun) 'RES1'
      write(iun) i2dum
      i2dum = natoms
      write(iun) i2dum
      i2dum = 1
      write(iun) i2dum

      rdum = 0.0e0
      write(iun) rdum
      write(iun) rdum
      write(iun) rdum
      rdum = 10.0e0
      write(iun) rdum
c
c write atom info
c
      if (icel.eq.1) then

         call setrr(alpha,beta,gamma,a,b,c,rr)

         do j=1,3
            write(iun) (real(trc(coo(1,i),rr,j)),i=1,natoms)
         end do
      else
         do j=1,3
            write(iun) (real(coo(j,i)*toang),i=1,natoms)
         end do
      endif
 
      write(iun) (mklab(ianz(i),i),i=1,natoms)

      i2dum = 1
      write(iun) (i2dum,i=1,natoms)

      if (iff.eq.6) then
         write(iun) (ityp(i),i=1,natoms)
      else
         write(iun) (ictyp(i,ianz(i),idochg,ianz,iconn),i=1,natoms)
      endif

      write(iun) (real(qat(i)),i=1,natoms)

c
c write the extra data
c

      rdum = 0.0e0
      write(iun) 'BVALUE    ','REAL',natoms,fname
      write(iun) (rdum,i=1,natoms)

      write(iun) 'CONNECT   ','BOND',natoms,fname

      do i=1,natoms
         nbc(i) = 0
         do j=1,iconn(1,i)
            if (iconn(1+j,i).gt.0) then
               nbc(i) = nbc(i) + 1
               ibc(i,nbc(i)) = iconn(1+j,i)
            endif
         end do
      end do
      write(iun) (nbc(i),(ibc(i,j),j=1,nbc(i)),i=1,natoms)

      nbnds = 0
      do i=1,natoms
         do j=1,iconn(1,i)
            jj = iconn(1+j,i)
            if (jj.gt.0) then
               if (jj.gt.i) then
                  nbnds = nbnds + 1
                  iorder(nbnds) = ibtyp(i,jj,idochg,1,ianz)
                  if (iorder(nbnds).eq.4) iorder(nbnds) = 7
                  ibc2(1,nbnds) = i
                  ibc2(2,nbnds) = jj
               endif
            endif
         end do
      end do

      write(iun) 'ORDER     ','BOND',nbnds,fname
      write(iun) (iorder(i),(ibc2(j,i),j=1,2),i=1,nbnds)

      if (icel.eq.1) then
          write(iun) spnam(nspg),'   ','SYMM',2,fname
          write(col,'(a5,6(f9.3),i2,i6,3x,a7)') 'CELL ',a,b,c,
     &        alpha/torad,beta/torad,gamma/torad,1,nspg,spnam(nspg)
          write(iun) col
          col = 'END'
          write(iun) col
          call fdat(1,0,0,0,0,0)
      endif

      do i=1,natoms
         iaton(i) = 1
      end do

      return
      end

      subroutine getrcn(iat,iconn,ianz)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /conrl/ ibnds,icnn(mxcon),io,in,ic,ih,ian1,ian2,ian3,ian4
      dimension iconn(mxcon+1,*),ianz(*)

      ibnds = 0
      io = 0
      in = 0
      ic = 0
      ih = 0
      ian1 = 0
      ian2 = 0
      ian3 = 0
      ian4 = 0

      if (iat.eq.0) return

      do i=1,iconn(1,iat)
         if (iconn(i+1,iat).gt.0) then
            ibnds = ibnds + 1
            icnn(ibnds) = iconn(i+1,iat)
            ia = ianz(icnn(ibnds))
            if (ia.eq.1) ih = ih + 1
            if (ia.eq.6) ic = ic + 1
            if (ia.eq.7) in = in + 1
            if (ia.eq.8) io = io + 1
         endif
      end do

      if (ibnds.gt.0) ian1 = ianz(icnn(1))
      if (ibnds.gt.1) ian2 = ianz(icnn(2))
      if (ibnds.gt.2) ian3 = ianz(icnn(3))
      if (ibnds.gt.4) ian4 = ianz(icnn(4))

      return
      end


      subroutine addcoo(ic,icc,rcoo)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (numatm=2000)
      common /extchg/ exchg(3,3),iexchg(3),nexchg
      common /coord / xyz(3,numatm)
      dimension tmp(3),iflg(numatm)

      nexchg = nexchg + 1
      iexchg(nexchg) = 1

      do i=1,3
         tmp(i) = xyz(i,ic) - xyz(i,icc)
      end do

      rc = vlen(tmp)

      do i=1,3
         tmp(i) = (tmp(i)*rcoo)/rc
      end do

      do i=1,3
         exchg(i,nexchg) = xyz(i,ic) + tmp(i)
      end do

      print*,"Added extra positive charge at coordinates (bohr):"
      print*,(exchg(jj,nexchg),jj=1,3)

      return
      end

      subroutine chkcod(kcoo,kcooh,ianz,iconn)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (numatm=2000)
      common /athlp/ iatoms, mxnat
      common /totchg/ itot
      dimension iconn(mxcon+1,*),ianz(*),icnn(mxcon),ifl(numatm)

      call qmtot

      do i=1,iatoms
         ifl(i) = 0
      end do

      kcoo = 0
      kcooh = 0

      do i=1,iatoms
         call ispn(irs,i,irng,idochg,1)
         if (ifl(i).eq.0.and.irs.eq.10) then
            kcoo = 1
            ic  = 0
            io  = i
            io2 = 0
            icc = 0
            ifl(i) = 1

            icooh = 0
            do j=1,iconn(1,i)
               jj = iconn(1+j,i)
               ian = ianz(iconn(1+j,i))
               if (ian.eq.6) ic = jj
               if (ian.eq.1) icooh = 1
            end do
            if (icooh.eq.0) then
               do j=1,iconn(1,ic)
                  jj = iconn(1+j,ic)
                  ian = ianz(jj)
                  if (ian.eq.8.and.jj.ne.i) io2 = jj
                  if (ian.eq.6) icc = jj
               end do
               if (io2.ne.0) then
                  ifl(io2) = 1
                  do j=1,iconn(1,io2)
                     ian = ianz(iconn(1+j,io2))
                     if (ian.eq.1) icooh = 1
                  end do
               endif
            endif
            if (icooh.eq.0.and.itot.lt.0) then
               kcooh = 0
               call addcoo(ic,icc,2.281d0)
            endif
            
         endif
      end do

      return
      end

      integer*2 function ictyp(iat,ian,idochg,ianz,iconn)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxt=14)
      parameter (mxmsf=235)
      character*4 atype
      common /atypes/ ihbt(mxt),atype(mxt)
      common /msft/   ltyp(mxmsf)
      common /conrl/ ibnds,icnn(mxcon),io,in,ic,ih,ian1,ian2,ian3,ian4
      dimension iconn(mxcon+1,*),ianz(*)

      call getrcn(iat,iconn,ianz)

      call ispn(irs,iat,irng,idochg,1)
      ihb = ihbt(irs)
      
      ictyp = 1
      if (ian.eq.1) then
         ictyp = 1
         if (ibnds.eq.1) then
            if (ian1.eq.6.or.ian1.eq.14) ictyp = 3
            if (ian1.eq.7) then
               ib = icred(icnn(1),idum1,idum2,ianz,iconn)
               if (ib.eq.4) ictyp = 2
            endif
            if (ian1.eq.8.or.ian1.eq.14) then
               if (ian1.eq.14) ictyp = 8
               if (ian1.eq.8) then
                  ib = 0
                  ic = 0
                  kat = icnn(1)
                  do j=1,iconn(1,kat)
                     jat = iconn(j+1,kat)
                     if (jat.gt.0) then
                        ib = ib + 1
                        if (ianz(jat).eq.6) ic = jat
                     endif
                  end do
                  if (ib.eq.2.and.ic.ne.0) then
                     ictyp = 8
                     ib = 0
                     io = 0
                     do j=1,iconn(1,ic)
                        jat = iconn(j+1,ic)
                        if (jat.gt.0.and.jat.ne.icnn(1)) then
                           ib = ib + 1
                           if (ianz(jat).eq.8) io = jat
                        endif
                     end do
                     if (ib.eq.2.and.io.ne.0) ictyp = 1
                  endif
               endif
            endif
         endif
      elseif (ian.eq.6) then
         if (ihb.eq.1) then
            ictyp = 18
         elseif (ihb.eq.2) then
            ictyp = 14
            if (ibnds.eq.3) then
               do i=1,3
                  call ispn(irs,icnn(i),irng,idochg,1)
                  if (ianz(icnn(i)).eq.6.and.
     &                ihbt(irs).eq.2)
     &                ictyp = 16
               end do
               if (io.gt.0) ictyp = 14
            endif
         elseif (ihb.eq.3) then
            ictyp = 10
         elseif (ihb.eq.4) then
            if (irng.eq.1) then
               ictyp = 22
            elseif (irng.eq.2) then
               ictyp = 21
            elseif (irng.eq.3) then
               ictyp = 27
            elseif (irng.eq.4) then
               ictyp = 25
            elseif (irng.eq.5) then
               ictyp = 26
            endif
         endif
      elseif (ian.eq.7) then
         if (ihb.eq.1) then
            ictyp = 31
         elseif (ihb.eq.2) then
            ictyp = 32
         elseif (ihb.eq.3) then
            ictyp = 36
            do i=1,ibnds
               call ispn(irs,icnn(i),jrng,idochg,1)
               iht = ihbt(irs)
               if (ianz(icnn(i)).eq.6.and.
     &             (iht.eq.2.or.iht.eq.4))
     &             ictyp = 32
            end do
            if (ibnds.eq.3.and.io.eq.2) ictyp = 38
         elseif (ihb.eq.4) then
            ictyp = 35
            if (irng.eq.2) ictyp = 34
         endif
      elseif (ian.eq.8) then
         if (ihb.eq.1) then
            ictyp = 48
         elseif (ihb.eq.2) then
            ictyp = 40
            if (ibnds.eq.1.and.ian1.ne.1) then
               jat = icnn(1)
               if (iconn(1,jat).eq.3.or.
     &           (ian1.eq.15.and.iconn(1,jat).eq.4)) then
                  ih = 0
                  ic = 0
                  ipo = 0
                  do i=1,iconn(1,jat)
                     kat = iconn(1+i,jat)
                     if (kat.gt.0.and.kat.ne.iat) then
                        if (ianz(kat).eq.8) then
                           if (ian1.eq.15) ipo = ipo + 1
                           ib = icred(kat,idum1,idum2,ianz,iconn)
                           if (ib.eq.2.and.ian1.eq.6) ictyp = 51
                           if (ib.eq.1) ictyp = 43
                        elseif (ianz(kat).eq.1) then
                           ih = ih + 1
                        elseif (ianz(kat).eq.6) then
                           ic = ic + 1
                        endif
                        if (ib.eq.0) then
                           if (ic.eq.2.and.ian1.eq.6) ictyp = 42
                           if (ic.eq.1.and.ih.eq.1.and.ian1.eq.6) 
     &                        ictyp = 41
                        endif
                     endif
                  end do
                  if (ipo.eq.4) ictyp = 43
               endif
            endif
         elseif (ihb.eq.3) then
            ictyp = 45
            if (ibnds.eq.2) then
               if (ian1.ne.1.and.ian2.ne.1) then
                  ictyp = 50
                  do i=1,2
                     jat = icnn(i)
                     if ((ianz(jat).ne.1.and.iconn(1,jat).eq.3).or.
     &               (ianz(jat).eq.15.and.iconn(1,jat).eq.4)) then
                        do j=1,iconn(1,jat)
                           kat = iconn(1+j,jat)
                           if (kat.gt.0.and.kat.ne.iat) then
                              if (ianz(kat).eq.8) then
                                 ib = icred(kat,idum1,idum2,ianz,iconn)
                                 if (ib.eq.1) ictyp = 49
                              endif
                           endif
                        end do
                     endif
                  end do
                  if (ian1.eq.14.and.ian1.eq.14) ictyp = 55
                  if (ian1.eq.13.and.ian1.eq.14) ictyp = 56
                  if (ian1.eq.14.and.ian1.eq.13) ictyp = 56
                  if (ian1.eq.13.and.ian1.eq.13) ictyp = 56
               endif
            endif
         elseif (ihb.eq.4) then
            ictyp = 53
            if (irng.eq.2) ictyp = 52
         endif
      elseif (ian.eq.15) then
         if (ihb.eq.1) then
            ictyp = 233
         elseif (ihb.eq.2) then
            ictyp = 62
         elseif (ihb.eq.3) then
            ictyp = 60
            if (ibnds.eq.4.or.ibnds.eq.3) then
               if (io.eq.3) ictyp = 61
               if (io.eq.4) ictyp = 62
            endif
         elseif (ihb.eq.4) then
            ictyp = 64
         endif
      elseif (ian.eq.16) then
         if (ihb.eq.1) then
            ictyp = 70
         elseif (ihb.eq.2) then
            ictyp = 75
         elseif (ihb.eq.3) then
            ictyp = 70
            if (ibnds.eq.2) then
               if (ian1.ne.1.and.ian2.ne.1.and.
     &             ian1.ne.16.and.ian2.ne.16) ictyp = 74
            endif
         elseif (ihb.eq.4) then
            ictyp = 73
            if (irng.eq.2) ictyp = 72
         endif
         if (io.eq.4) then
            ictyp = 79
         elseif (io.eq.3) then
            ictyp = 78
         elseif (io.eq.2) then
            ictyp = 77
         elseif (io.eq.1) then
            ictyp = 76
         endif
      elseif (ian.eq.99) then
         ictyp = 489
      else
         do i=1,mxmsf
            if (ltyp(i).eq.ian) ictyp = i
         end do
      endif
      
      return
      end


      integer*2 function igtyp(iat,ian,idochg,ianz,iconn)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxt=14)
      character*4 atype
      common /atypes/ ihbt(mxt),atype(mxt)
      dimension icnn(mxcon),iccn(mxcon),icnn2(mxcon),iconn(mxcon+1,*),
     &          ianz(*)

      ibnds = 0
      ih = 0
      ihc = 0
      io = 0
      in = 0
      ic = 0
      is = 0
      ian1 = 0
      ian2 = 0
      ian3 = 0

      do i=1,iconn(1,iat)
         if (iconn(i+1,iat).gt.0) then
            ibnds = ibnds + 1
            icnn(ibnds) = iconn(i+1,iat)
            ia = ianz(icnn(ibnds))
            if (ia.eq.1) ih = ih + 1
            if (ia.eq.6) ic = ic + 1
            if (ia.eq.6.or.ia.eq.1) ihc = ihc + 1
            if (ia.eq.7) in = in + 1
            if (ia.eq.8) io = io + 1
            if (ia.eq.16) is = is + 1
         endif
      end do
      if (ibnds.gt.0) ian1 = ianz(icnn(1))
      if (ibnds.gt.1) ian2 = ianz(icnn(2))
      if (ibnds.gt.2) ian3 = ianz(icnn(3))

      call ispn(irs,iat,irng,idochg,1)
      
c irs    1      2      3      4      5      6       7
c      '    ','.1  ','.2  ','.3  ','.4  ','.ar ','.cat'
c irs    8      9     10     11     12     13      14
c      '.am ','.pl3','.co2','.spc','.t3p','.O  ','.O2 '

      igtyp = 1
      if (ian.eq.1) then

         igtyp = 1
         if (ibnds.eq.1) then
            if (ian1.eq.6) then
               igtyp = 25
               call ispn(itmp,icnn(1),irng,idochg,1)
               if (itmp.eq.6.or.itmp.eq.3) igtyp = 24

               call flth(icnn(1),icnn2,ibnds2,iconn)

               nhcp = 0
               do i=1,ibnds2
                   ia = ianz(icnn2(i))
                   if (.not.(ia.eq.6.or.ia.eq.1.or.ia.eq.15)) then
                      nhcp = nhcp + 1
                      iccn(nhcp) = icnn2(i)
                   endif
               end do

               if (nhcp.ne.0) then
                  if (igtyp.eq.25) then
                     igtyp = 18 + nhcp
                  elseif (igtyp.eq.24) then
                     igtyp = 21 + nhcp
                  endif
               endif
               if (nhcp.eq.1.and.ianz(iccn(1)).eq.8.and.
     &             (itmp.eq.6.or.itmp.eq.3)) then
                  igtyp = 24
               endif
            endif
            if (ian1.eq.7)  igtyp = 26
            if (ian1.eq.8)  then
               igtyp = 27
               call flth(icnn(1),icnn2,ibnds2,iconn)
               if (ibnds2.eq.2) then
                  if (ianz(icnn2(1)).eq.1.and.ianz(icnn2(2)).eq.1)
     &               igtyp = 30
               endif
            endif
            if (ian1.eq.15) igtyp = 28
            if (ian1.eq.16) igtyp = 29
         endif

      elseif (ian.eq.6) then

         if (irs.eq.2) igtyp = 3
         if (irs.eq.3) igtyp = 4
         if (irs.eq.4) igtyp = 5
         if (irs.eq.6) igtyp = 6
         if (irs.eq.7) then
            igtyp = 6
            if (in.eq.3) igtyp = 72
         endif

         if ((irs.eq.3.or.irs.eq.6).and.ibnds.eq.3.and.
     &      (io.ge.1.or.is.ge.1)) then
            icon = 0
            if (ian1.eq.8.or.ian1.eq.16) icon = icnn(1)
            if (ian2.eq.8.or.ian2.eq.16) icon = icnn(2)
            if (ian3.eq.8.or.ian3.eq.16) icon = icnn(3)
            if (icon.gt.0) then
               call ispn(itmp,icon,irng,idochg,1)
               if (itmp.eq.10.or.itmp.eq.3) igtyp = 2
            endif
         endif

      elseif (ian.eq.7) then

c n1     (.1)
         if (irs.eq.2) igtyp = 37
c n2, nc (.2,.ar) (nc contributes 1 el to pi system)
         if ((irs.eq.3.or.irs.eq.6).and.ibnds.eq.2) then
            igtyp = 38
            if (irs.eq.6) igtyp = 43
         endif
c na     (.2,.ar,.pl) (.pl contributes 2 el to pi system)
         if ((irs.eq.3.or.irs.eq.6.or.irs.eq.9).and.ibnds.eq.3) 
     &      igtyp = 41
c n3     (.3)
         if (irs.eq.4.and.ibnds.eq.3) igtyp = 39
c no     (.2,.3,.pl)
         if ((irs.eq.3.or.irs.eq.4.or.irs.eq.9).and.
     &      ibnds.eq.3.and.io.eq.2) igtyp = 48
c n4     (.4)
         if (irs.eq.5) igtyp = 40
c n      (.am)
         if (irs.eq.8) igtyp = 36
c nh     (.3)
         if ((irs.eq.4.or.irs.eq.9).and.ibnds.eq.3.and.ihc.eq.3) then
            iar = 0
            if (ian1.ne.1) then
               call ispn(itmp,icnn(1),irng,idochg,1)
               if (itmp.eq.6) iar = iar + 1
            endif
            if (ian2.ne.1) then
               call ispn(itmp,icnn(2),irng,idochg,1)
               if (itmp.eq.6) iar = iar + 1
            endif
            if (ian3.ne.1) then
               call ispn(itmp,icnn(3),irng,idochg,1)
               if (itmp.eq.6) iar = iar + 1
            endif
            if (iar.gt.0) igtyp = 47
         endif
         if (ibnds.eq.3.and.ic.gt.0) then
            do i=1,ibnds
               nat = ianz(icnn(i))
               if (nat.eq.6) then
                  call flth(icnn(i),icnn2,ibnds2,iconn)
                  if (ibnds2.eq.3) then
                      iscat = 1
                      do j=1,3
                         if (ianz(icnn2(j)).ne.7) iscat = 0
                      end do
                      if (iscat.eq.1) igtyp = 47
                  endif
               endif
            end do
         endif

      elseif (ian.eq.8) then

         if (irs.eq.3.or.irs.eq.10) igtyp = 49
         if ((irs.eq.4.or.irs.eq.6).and.ibnds.eq.2) then
            if (ih.eq.2) then 
               igtyp = 52
            elseif (ih.eq.1) then
               igtyp = 50
            else
               igtyp = 51
            endif
         endif

      elseif (ian.eq.15) then

         if (ibnds.le.2) igtyp = 53
         if (ibnds.eq.3) then

c difference between p3 and p4 may need some adjustment

            iar = 0
            if (ian1.ne.1) then
               call ispn(itmp,icnn(1),irng,idochg,1)
               if (itmp.eq.3) iar = iar + 1
            endif
            if (ian2.ne.1) then
               call ispn(itmp,icnn(2),irng,idochg,1)
               if (itmp.eq.3) iar = iar + 1
            endif
            if (ian3.ne.1) then
               call ispn(itmp,icnn(3),irng,idochg,1)
               if (itmp.eq.3) iar = iar + 1
            endif
            if (iar.gt.0) then
               igtyp = 55
            else
               igtyp = 54
            endif
         endif
         if (ibnds.ge.4) igtyp = 56

      elseif (ian.eq.16) then

         if (ibnds.eq.1) igtyp = 64
c ibnds = 2 may need some refinement, sh clear, ss clear, but s2 ?
         if (ibnds.eq.2) then
            if (ih.ge.1) then
               igtyp = 68
            else
               if ((irs.eq.4.or.irs.eq.6).or.is.eq.1) then
                  igtyp = 69
               else
                  igtyp = 65
               endif
            endif
         endif
         if (ibnds.eq.3) igtyp = 66
         if (ibnds.ge.4) igtyp = 67

      elseif (ian.eq.9) then
         if (ibnds.eq.1) igtyp = 32
      elseif (ian.eq.17) then
         if (ibnds.eq.1) igtyp = 33
      elseif (ian.eq.35) then
         if (ibnds.eq.1) igtyp = 34
      elseif (ian.eq.53) then
         if (ibnds.eq.1) igtyp = 35
      endif
      
      return
      end

      integer function icred(iat,inoh,ih,ianz,iconn)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      dimension iconn(mxcon+1,*),ianz(*)

      icred = 0
      inoh = 0
      ih = 0
      do j=1,iconn(1,iat)
         jat = iconn(j+1,iat)
         if (jat.gt.0) then
            icred = icred + 1
            if (ianz(jat).eq.1) then
               ih = ih + 1
            else
               inoh = inoh + 1
            endif
         endif
      end do

      return
      end

      character*6 function mklab(ian,lab)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxel=100)
      character*2 elemnt,ggstr
      common /elem/elemnt(mxel)

      mklab = '      '
      mklab(1:2) = elemnt(ian)
      if (lab.gt.99) return
      if (mklab(2:2).eq.' ') then
          mklab(2:3) = ggstr(lab)
      else
          mklab(3:4) = ggstr(lab)
      endif

      return
      end

      integer function ifmt(in)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      integer teni

      teni = 1
      ifmt = 1

      do while(.true.)
         teni = teni * 10
         ifmt = ifmt + 1
         if (mod(in,teni) .eq. in)  return
      end do

      return
      end

      subroutine wrtnd(iun,ianz,iaton,iconn,lring,ityp,coo)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxchtp=136)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxamo=201)
      parameter (mxel=100)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      character*2 elemnt
      common /elem/  elemnt(mxel)
      common /cllab/ iclon,iclpnt(4)
      integer*2 ityp
      common /types/ iff
      character*3 pdbsym,hsym,chtnk,chtmp,ambtnk
      character*2 amotnk,gffstr
      common /symbol/ pdbsym(mxsym),hsym(mxhsym),chtnk(mxchtp),
     &                ambtnk(mxamb),amotnk(mxamo),gffstr(mxgff)
      logical isqmmm
      integer*2 itypi,it10000
      character*40 fmtstr
      dimension icnn(mxcon),lring(*)
      dimension coo(3,*),ianz(*),iaton(*),iconn(mxcon+1,*),ityp(*)


      it10000=10000

      if (iun.gt.100) then
         isqmmm = .true.
         iun = iun - 100
      else
         isqmmm = .false.
      endif

c write tinker .xyz file, charmm atom types are written

      toang = 0.52917706d0

      natoms = 0
      ioffmx = mxnat
      ioffmn = 0
      do i=1,iatoms
         if (ianz(i).lt.100.and.ianz(i).gt.0) natoms = natoms + 1
         if (ianz(i).eq.100) then
             ioffmx = i
             if (ioffmn.eq.0) ioffmn = i
         endif
      end do
      ioffmn = ioffmx - ioffmn + 1

      do i=1,iatoms
         lring(i) = 0
         iaton(i) = 2
      end do

      fmtstr = '(i5,1x,a2,1x,3(f12.6),1x,i6,1x,8(i5,1x))'
      ifm = ifmt(natoms)
      fmtstr(3:3) = char(48+ifm)
      fmtstr(35:35) = char(48+ifm)
      if (iff.ge.2.and.iff.le.4) then
         fmtstr(9:9) = char(51)
      endif

c tinker

      if (iff.eq.2) then
         write(iun,*) natoms,
     &   ' molden generated tinker .xyz (charmm param.)'
      elseif (iff.eq.3) then
         write(iun,*) natoms,
     &   ' molden generated tinker .xyz (amber param.)'
      elseif (iff.eq.4) then
         write(iun,*) natoms,
     &   ' molden generated tinker .xyz (amoeba param.)'
      else
         write(iun,*) natoms,
     &   ' molden generated tinker .xyz (mm3 param.)'
      endif

      do i=1,iatoms

         if (ianz(i).ne.100) then
            ibnds = 0
            do j=1,iconn(1,i)
               if (iconn(j+1,i).gt.0) then
                  ibnds = ibnds + 1
                  if (iconn(j+1,i).gt.ioffmx) then
                     icnn(ibnds) = iconn(j+1,i) - ioffmn
                  else
                     icnn(ibnds) = iconn(j+1,i)
                  endif
               endif
            end do

CNF 	QM atoms are at the H level, while MM atoms are at M and L levels               
            if (isqmmm. and. ityp(i).lt.10000) then
                itypi = ityp(i) + 20000
            else if(isqmmm) then
                itypi = ityp(i) - (ityp(i)/10000)*10000
            else
                itypi = mod(ityp(i),it10000)
            endif

            if (iff.ge.2.and.iff.le.4) then

               chtmp = elemnt(ianz(i))//' '
               if (iff.eq.2) then
                  if (ityp(i).gt.0.and.(mod(itypi,it10000)).le.mxchtp)
     &               chtmp = chtnk(mod(itypi,it10000))
               endif
               if (iff.eq.3) then
                  if (ityp(i).gt.0.and.(mod(itypi,it10000)).le.mxamb)
     &            then
                     chtmp = ambtnk(mod(itypi,it10000))
                     if (itypi.gt.mxamb) itypi = itypi + (2000-mxamb)
                  endif
               endif
               if (iff.eq.4) then
                  if (ityp(i).gt.0.and.(mod(itypi,it10000)).le.mxamo)
     &               chtmp = amotnk(mod(itypi,it10000))
               endif
c               write(iun,'(i6,2x,a3,1x,3(f12.6),1x,i6,1x,8i5)')
               write(iun,fmtstr)
     &            i,chtmp,(coo(j,i)*toang,j=1,3),
     &            itypi,(icnn(j),j=1,ibnds)

            elseif (iff.eq.1) then

c               write(iun,'(i5,1x,a2,1x,3(f12.6),1x,i6,1x,8i5)')
               write(iun,fmtstr)
     &            i,elemnt(ianz(i)),(coo(j,i)*toang,j=1,3),
     &            itypi,(icnn(j),j=1,ibnds)
            else

               itypi = mmtyp(i,ianz(i),0,ianz,iaton,iconn)
               if (isqmmm. and. ityp(i).lt.10000) then
                  itypi = itypi + 20000
               elseif(isqmmm) then
                  itypi = itypi - (itypi/10000)*10000
               else
                  itypi = mod(itypi,it10000)
               endif
c               write(iun,'(i5,1x,a2,1x,3(f12.6),1x,i6,1x,8i5)')
               write(iun,fmtstr)
     &            i,elemnt(ianz(i)),(coo(j,i)*toang,j=1,3),
     &            itypi,(icnn(j),j=1,ibnds)

            endif


         endif

      end do

      do i=1,iatoms
         iaton(i) = 1
      end do

      return
      end

      subroutine reordc(iatom,ibnds,icnn,iconn,nmet,imet,lring,ityp)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      logical doit
      integer*2 ityp
      dimension icnn(mxcon), iconn(mxcon+1,*),imet(*),lring(*),ityp(*)

      ibnds = 0
      do i=1,nmet
          if (iatom.eq.imet(i)) return
      end do

      do j=1,iconn(1,iatom)
         jj = iconn(j+1,iatom)

         if (jj.gt.0) then
            doit = .true.
            do i=1,nmet
               if (jj.eq.imet(i)) doit = .false.
            end do
            if (ityp(iatom).eq.649.and.ityp(jj).ne.650) doit = .false.
            if (ityp(iatom).ne.650.and.ityp(jj).eq.649) doit = .false.
            if (doit) then
               ibnds = ibnds + 1
               icnn(ibnds) = lring(iconn(j+1,iatom))
            endif
         endif
      end do

      return
      end

      subroutine wraone(iun,iatom,ibnds,icnn,ianz,isurf,ityp,coo,q)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxchtp=136)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxamo=201)
      character*3 pdbsym,hsym,chtnk,chtmp,ambtnk
      character*2 amotnk,gffstr
      common /symbol/ pdbsym(mxsym),hsym(mxhsym),chtnk(mxchtp),
     &                ambtnk(mxamb),amotnk(mxamo),gffstr(mxgff)
      character*80 tnknm
      integer tnkbg,tnkit,tnkarc,tnkarf,tnkprg
      common /tnkopt/ rmsgrd,tnkbg,tnkit,tnkarc,tnkarf,tnkprg,tnknm,icst
      parameter (mxel=100)
      character*2 elemnt
      common /elem/  elemnt(mxel)
      integer*2 ityp
      parameter (mxcon=10)
      dimension icnn(mxcon)
      dimension coo(3,*),q(*),ianz(*),ityp(*),isurf(*)

      toang = 0.52917706d0

      if (icst.eq.1) then
         iopt = isurf(iatom)
      else
         iopt = 1
      endif

      chtmp = elemnt(ianz(iatom))//' '
      itypi = ityp(iatom)

      if (itypi.gt.0) then


c Amber
          if (itypi.le.mxamb) then
            chtmp = ambtnk(itypi)
            if (itypi.gt.mxamb) itypi = itypi + (2000-mxamb)
          endif

          write(iun,'(i6,2x,a3,1x,3(f9.3),1x,i4,1x,8(i6,1x))')
     &      iopt,chtmp,(coo(j,iatom)*toang,j=1,3),
     &      itypi,(icnn(j),j=1,ibnds)

      elseif (itypi.lt.0) then
c GAFF
          if (iabs(itypi).le.mxgff) then
            chtmp = gffstr(iabs(itypi))//' '
          endif

          write(iun,
     &      '(i6,2x,a3,1x,3(f9.3),1x,i3,1x,f6.3,1x,8(i6,1x))')
     &      iopt,chtmp,(coo(j,iatom)*toang,j=1,3),
     &      itypi,q(iatom),(icnn(j),j=1,ibnds)

      elseif (itypi.eq.0) then

          write(iun,
     &      '(i6,2x,a3,1x,3(f9.3),1x,i3,1x,f6.3,1x,8(i6,1x))')
     &      iopt,chtmp,(coo(j,iatom)*toang,j=1,3),
     &      itypi,q(iatom),(icnn(j),j=1,ibnds)

      endif


      return
      end

      logical function ismet(ian)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mexcl=10)
      common /metexc/ qexcl(mexcl),ianexc(mexcl)

      ismet = .false.

      do i=1,mexcl
         if (ian.eq.ianexc(i)) ismet = .true.
      end do

      return
      end

      subroutine wrgfd(iun,ianz,iaton,iconn,isurf,lring,lwrit,ncalf,
     &                 ishoh,iresid,ityp,coo,q,
     &                 xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxchtp=136)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxamo=201)
      parameter (mxcon=10)
      parameter (numcal=50000)
      parameter (mxmet=500)
      common /athlp/ iatoms, mxnat
      integer*2 ityp
      logical ismet
      common /types/ iff
      common /pbc/ abc(3),ibox,icell,igfmap
      character*3 pdbsym,hsym,chtnk,chtmp,ambtnk
      character*2 amotnk,gffstr
      common /symbol/ pdbsym(mxsym),hsym(mxhsym),chtnk(mxchtp),
     &                ambtnk(mxamb),amotnk(mxamo),gffstr(mxgff)
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb
      logical odupl
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension icnn(mxcon),ihcnn(mxcon)
      dimension coo(3,*),q(*),ianz(*),iaton(*),iconn(mxcon+1,*),
     &          ityp(*),isurf(*),lring(*),lwrit(*),iresid(*)
      dimension iresb(numcal),iresn(numcal),imet(mxmet)

c write AMBFOR .xyz file, amber/gaff atom types are written


      if (iff.ne.7) then
         iff = 7
         call dotyp(0)
      endif

      torad = datan(1.0d0) / 45.0d0
      nuntyp = 0
      nmet = 0

      do i=1,iatoms

         if (ianz(i).lt.100.and.ianz(i).gt.0) then
             if (ityp(i).eq.0.or.ityp(i).eq.-1) then
                print*,'atom ',i,' unknown ',ianz(i)
                nuntyp = nuntyp + 1
             endif
         endif
         if (nmet.lt.mxmet) then
             if (ismet(ianz(i))) then
                 nmet = nmet + 1
                 imet(nmet) = i
             endif
         endif

      end do

      if (nuntyp.ne.0) then
         print*,nuntyp,' untyped atoms '
         print*,' '
         print*,'Remember, the only metals supported by amber/gaff'
         print*,'are:'
         print*,' '
         do i=651,658
            print*,ambtnk(i)
         end do
      endif


c     determine new order of atoms to be written (kept in lring)

      do i=1,iatoms
         lwrit(i) = 0
      end do

      natoms = 0
      do i=1,ncalf
         
         call getpdb(i,ipdb,ihpdb)
         do j=1,mxsym
             if (ipdb(j).ne.0) then
                natoms = natoms + 1
                lring(ipdb(j)) = natoms
                lwrit(ipdb(j)) = 1
             endif
         end do
         do j=1,mxhsym*3
             if (ihpdb(j).ne.0) then
                natoms = natoms + 1
                lring(ihpdb(j)) = natoms
                lwrit(ihpdb(j)) = 1
             endif
         end do

      end do

      do i=1,iatoms
         if (ianz(i).ne.100.and.lwrit(i).eq.0) then
            if (ianz(i).eq.8.and.ityp(i).eq.649) then
               natoms = natoms + 1 
               lring(i) = natoms
               lwrit(i) = 1
               do j=1,iconn(1,i)
                  jj = iconn(j+1,i)
                  if (jj.gt.0) then
                     if (ianz(jj).eq.1.and.ityp(jj).eq.650) then
                        natoms = natoms + 1 
                        lring(jj) = natoms
                        lwrit(jj) = 1
                     endif
                  endif
               end do
            else
               if (lwrit(i).eq.0.and.
     &            odupl(i,coo,xa,ya,yb,za,zb,zc)) then
                  natoms = natoms + 1 
                  lring(i) = natoms
                  lwrit(i) = 1
               endif
            endif
         endif
      end do

c ambfor amber + gaff

      nwramb = natoms

      if (ibox.eq.1) then
         write(iun,'(a,3(f9.3,1x))')
     &      '[AMBFOR] box ',abc(1),abc(2),abc(3)
      else
         if (icell.eq.1) then
             write(iun,'(a,6(f9.3,1x))')
     &          '[AMBFOR] cell ',a,b,c,
     &          alpha/torad,beta/torad,gamma/torad
         else
             write(iun,*) '[AMBFOR]'
         endif
      endif

      write(iun,*) natoms,
     &   ' molden generated ambfor .xyz '

      do i=1,iatoms
         lwrit(i) = 0
      end do

      ires = 0
      do i=1,ncalf

         iresb(i) = 0
         call getpdb(i,ipdb,ihpdb)

         do j=1,mxsym
             if (ipdb(j).ne.0) then

                 call reordc(ipdb(j),ibnds,icnn,iconn,nmet,imet,lring,
     &                       ityp)
                 call wraone(iun,ipdb(j),ibnds,icnn,ianz,isurf,
     &                       ityp,coo,q)

                 lwrit(ipdb(j)) = 1
                 if (iresb(i).eq.0) then
                     iresb(i) = lring(ipdb(j))
                     iresn(i) = i
                     ires = iresid(ipdb(j))
                 endif
             endif
         end do

         do j=1,mxhsym*3
             if (ihpdb(j).ne.0) then
                 call reordc(ihpdb(j),ibnds,icnn,iconn,nmet,imet,lring,
     &                       ityp)
                 call wraone(iun,ihpdb(j),ibnds,icnn,ianz,isurf,
     &                       ityp,coo,q)

                 lwrit(ihpdb(j)) = 1
             endif
         end do

      end do

      nres = ncalf
      nhoh = 0
      do i=1,iatoms

         if (ianz(i).ne.100.and.lwrit(i).eq.0) then

            call reordc(i,ibnds,icnn,iconn,nmet,imet,lring,ityp)
            call wraone(iun,i,ibnds,icnn,ianz,isurf,ityp,coo,q)
            lwrit(i) = 1

            if (iresid(i).ne.ires) then
               nres = nres + 1
               iresb(nres) = lring(i)
               iresn(nres) = iresid(i)
               ires = iresid(i)
            elseif (iresid(i).eq.-ishoh.and.ityp(i).eq.649) then
               nres = nres + 1
               nhoh = nhoh + 1
               iresb(nres) = lring(i)
               iresn(nres) = -(ishoh+nhoh)
            endif

            if (ianz(i).eq.8.and.ityp(i).eq.649) then
               do j=1,iconn(1,i)
                  jj = iconn(j+1,i)
                  if (jj.gt.0) then
                     if (ianz(jj).eq.1.and.ityp(jj).eq.650) then
                        call reordc(jj,ibnds,icnn,iconn,nmet,imet,lring,
     &                              ityp)
                        call wraone(iun,jj,ibnds,icnn,ianz,isurf,
     &                              ityp,coo,q)
                        lwrit(jj) = 1
                     endif
                  endif
               end do

            endif

         endif

      end do

      write(iun,'(a)') '[RESIDUES]'

      do i=1,nres
         write(iun,'(i5,1x,i5)') iresn(i),iresb(i)
      end do

      do i=1,iatoms
         lwrit(i) = lring(i)
      end do

      do i=1,iatoms
         il = lwrit(i)
         if (il.gt.0.and.il.le.mxnat) then
            lring(lwrit(i)) = i
         endif
      end do

      do i=natoms+1,iatoms
         lring(i) = 0
      end do

      return
      end

      subroutine wrogd(iun,rx,ry,rz,
     &                 ianz,iaton,iatclr,iresid,iconn,lring,coo,rzp,
     &                 reson,ianf,nchain,ncalf)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxel=100)
      parameter (mxcon=10)
      parameter (ncmx=32)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /athlp/ iatoms, mxnat
      common /surf/  natorg,noscnd
      character*2 elemnt
      common /elem/  elemnt(mxel)
      logical dozme
      common /getpnt/irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &               iconv,ircus,dozme
      integer dolabs,fancy,shade,atcol,persp,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      common /vropt/ ivtwo,ihand,ivadd
      integer reson
      common /cllab/ iclon,iclpnt(4)
      logical doit
      common /gracom/ uscl,colscd,colscpd,ivdwpl
      common /strips/ qnormo(3),crpnto(3,ncmx),crnrmo(3,ncmx),
     &                numcir,nquad
      common /vrcol/ jcol(3,16)
      common /setogl/ idirogl
      dimension rx(3),ry(3),rz(3)
      dimension icnn(mxcon),t(3),tmp1(3),tmp2(3),tmp3(3),tlpos(3)
      dimension rxyzt(3,3),ianz(*),iresid(*),iaton(*),iatclr(*),
     &          coo(3,*),rzp(*),iconn(mxcon+1,*),lring(*),reson(*),
     &          ianf(*)

      ihand = 0

c write moldenogl .ogl file or ribbon povray file

      toang = 0.52917706d0

      do i=1,iatoms
         lring(i) = 0
      end do

      natoms = 0
      do i=1,iatoms

         doit = .false.
         if (ipdbon.eq.1) then
             if (iresid(i).gt.0) then
                if (reson(iresid(i)).eq.1) 
     &            doit = .true. 
             elseif (iresid(i).lt.-3) then
                  doit = .true. 
             endif
         else
             doit = .true. 
         endif

         if (ianz(i).lt.100.and.ianz(i).gt.0.and.iaton(i).gt.0
     &       .and.doit) then
             natoms = natoms + 1
             lring(i) = natoms
         endif

      end do


      if (ivtwo.eq.3) then

c ogl file

        if (atcol.eq.1) then
         if (ivdwpl.eq.1.and.fancy.eq.1) then
            if (numcir.gt.8) then
               write(iun,'(a)') 
     &           '[MOLECULE] UNSCALED SPACEFILL CONN HIGH'
            else
               write(iun,'(a)') '[MOLECULE] UNSCALED SPACEFILL CONN'
            endif
         else
            if (numcir.gt.8) then
               write(iun,'(a)') '[MOLECULE] UNSCALED CONN HIGH'
            else
               write(iun,'(a)') '[MOLECULE] UNSCALED CONN'
            endif
         endif
        else
         if (ivdwpl.eq.1.and.fancy.eq.1) then
            if (numcir.gt.8) then
               write(iun,'(a)') 
     &          '[MOLECULE] UNSCALED GRPCOL SPACEFILL CONN HIGH'
            else
               write(iun,'(a)') 
     &          '[MOLECULE] UNSCALED GRPCOL SPACEFILL CONN'
            endif
         else
            if (numcir.gt.8) then
               write(iun,'(a)') 
     &          '[MOLECULE] UNSCALED GRPCOL CONN HIGH'
            else
               write(iun,'(a)') '[MOLECULE] UNSCALED GRPCOL CONN'
            endif
         endif
        endif
        write(iun,*) natoms

      elseif (ivtwo.eq.2) then

        call cntvec(t,coo,ianz,iatoms)

c        sc1 = vlen(t)*1.4d0

        rzpmax = -10000.d0
        do i=1,iatoms
           if (rzp(i).gt.rzpmax) rzpmax = rzp(i)
        end do

        sc1 = rzpmax*1.8d0

c transpose rotation matrix

        do i=1,3
           rxyzt(1,i) = rx(i)
           rxyzt(2,i) = ry(i)
           rxyzt(3,i) = rz(i)
        end do

        g = rxyzt(1,2) 
        rxyzt(1,2) = rxyzt(2,1)
        rxyzt(2,1) = g

        g = rxyzt(1,3) 
        rxyzt(1,3) = rxyzt(3,1)
        rxyzt(3,1) = g

        g = rxyzt(2,3) 
        rxyzt(2,3) = rxyzt(3,2)
        rxyzt(3,2) = g


c position of the eye in tmp1

        tmp2(1) = 0.0d0
        tmp2(2) = 0.0d0
        tmp2(3) = sc1

        tmp1(1) = 
     &    tmp2(1)*rxyzt(1,1)+tmp2(2)*rxyzt(1,2)+tmp2(3)*rxyzt(1,3)
        tmp1(2) = 
     &    tmp2(1)*rxyzt(2,1)+tmp2(2)*rxyzt(2,2)+tmp2(3)*rxyzt(2,3)
        tmp1(3) = 
     &    tmp2(1)*rxyzt(3,1)+tmp2(2)*rxyzt(3,2)+tmp2(3)*rxyzt(3,3)


c eye looks at center tmp2

        tmp2(1) = 0.0d0
        tmp2(2) = 0.0d0
        tmp2(3) = 0.0d0
        
c position of the light in tlpos

        tmp3(1) = dsqrt(3.0d0)/3.0d0*sc1
        tmp3(2) = dsqrt(3.0d0)/3.0d0*sc1
        tmp3(3) = dsqrt(3.0d0)/3.0d0*sc1

        tlpos(1) = 
     &    tmp3(1)*rxyzt(1,1)+tmp3(2)*rxyzt(1,2)+tmp3(3)*rxyzt(1,3)
        tlpos(2) = 
     &    tmp3(1)*rxyzt(2,1)+tmp3(2)*rxyzt(2,2)+tmp3(3)*rxyzt(2,3)
        tlpos(3) = 
     &    tmp3(1)*rxyzt(3,1)+tmp3(2)*rxyzt(3,2)+tmp3(3)*rxyzt(3,3)


        call plcini
        call plphd(iun,ihand,tmp1,tmp2,tlpos)

      endif
      call cntvec(t,coo,ianz,iatoms)

c  the molecule

      ilin = 0
      do i=1,iatoms

         if (i.gt.natorg.and.natorg.ne.0) ilin = 1

         if (ianz(i).lt.100) then

           ibnds = 0
           do j=1,iconn(1,i)
              if (iconn(j+1,i).gt.0) then
                 if (iaton(iconn(j+1,i)).gt.0) then
                    if (lring(iconn(j+1,i)).gt.0) then
                       ibnds = ibnds + 1
                       icnn(ibnds) = lring(iconn(j+1,i))
                    endif
                 endif
              endif
           end do

           doit = .false.
           if (ipdbon.eq.1) then
               if (iresid(i).gt.0) then
                  if (reson(iresid(i)).eq.1) 
     &              doit = .true. 
               elseif (iresid(i).lt.-3) then
                    doit = .true. 
               endif
           else
               doit = .true. 
           endif

           if (iaton(i).gt.0.and.doit) then

             do j=1,3
                 tmp1(j) = (coo(j,i)-t(j))*toang
             end do

             ia = ianz(i)
             if (atcol.eq.1) then

              if (ivtwo.eq.3) then

               write(iun,'(i3,1x,3(f12.6),1x,i2,1x,8i5)')
     &         ianz(i),(tmp1(j),j=1,3),ibnds,
     &         (icnn(k),k=1,ibnds)

              elseif (ivtwo.eq.2) then

               ic = icol(ia)
               call plvsph(iun,jcol,ic,tmp1,vdwr(ia)/toang)

              endif

             else

              if (ivtwo.eq.3) then

               write(iun,'(i3,1x,i2,1x,3(f12.6),1x,i2,1x,8i5)')
     &         ianz(i),iatclr(i),(tmp1(j),j=1,3),
     &         ibnds,(icnn(k),k=1,ibnds)

              elseif (ivtwo.eq.2) then

               ic = iatclr(i)
               call plvsph(iun,jcol,ic,tmp1,vdwr(ia)/toang)

              endif
             endif

           endif


         endif

      end do

      if (ivtwo.eq.2) then
         write(iun,*) '}'
         write(iun,*) 'molecule'
      endif

c lines

      if (ilin.eq.1.and.ivtwo.eq.3) then
         write(iun,'(a)') '[LINES]'
         do i=1,iatoms
            if (i.gt.natorg.and.natorg.ne.0) then
               do j=1,iconn(1,i)
                  jj = abs(iconn(j+1,i))
                  if (jj.lt.i) then
                    write(iun,'(i2,1x,3(f12.6),1x,3(f12.6))')
     &              iatclr(i),((coo(k,i)-t(k))*toang,k=1,3),
     &              ((coo(k,jj)-t(k))*toang,k=1,3)
                  endif
               end do
            endif
         end do
      endif

c ribbons

100   if (ipdbon.eq.1) then

         if (ivtwo.eq.3) then

          write(iun,'(a)') '[COL STRANDTOP] 1.0 0.0 1.0'
          write(iun,'(a)') '[COL STRANDBOTTOM] 1.0 0.0 1.0'
          write(iun,'(a)') '[COL HELIXOUT] 0.0 1.0 0.0'
          write(iun,'(a)') '[COL HELIXIN] 0.6 0.6 0.6'
          write(iun,'(a)') '[COL RNA] 0.5 1.0 0.5'
          write(iun,'(a)') '[COL COIL] 1.0 1.0 1.0'

          idirogl = 0
          call ribgll(iun,ianf,nchain,ncalf,iatoms)


         elseif (ivtwo.eq.2) then

          write(iun,*) '#declare STRANDTOP = texture {'
          write(iun,*) 'pigment { color rgb<0.0, 0.0, 1.0> }'
          write(iun,*)
     &       'finish { ambient 0.4 diffuse 0.4 specular 0.9}'
          write(iun,*) '}'

          write(iun,*) '#declare STRANDBOTTOM = texture {'
          write(iun,*) 'pigment { color rgb<1.0, 0.0, 1.0> }'
          write(iun,*)
     &       'finish { ambient 0.4 diffuse 0.4 specular 0.9}'
          write(iun,*) '}'

          write(iun,*) '#declare HELIXOUT = texture {'
          write(iun,*) 'pigment { color rgb<0.0, 1.0, 0.0> }'
          write(iun,*)
     &       'finish { ambient 0.4 diffuse 0.4 specular 0.9}'
          write(iun,*) '}'

          write(iun,*) '#declare HELIXIN = texture {'
          write(iun,*) 'pigment { color rgb<0.6, 0.6, 0.6> }'
          write(iun,*)
     &       'finish { ambient 0.4 diffuse 0.4 specular 0.9}'
          write(iun,*) '}'

          write(iun,*) '#declare RNA = texture {'
          write(iun,*) 'pigment { color rgb<0.5, 1.0, 0.5> }'
          write(iun,*)
     &       'finish { ambient 0.4 diffuse 0.4 specular 0.9}'
          write(iun,*) '}'

          write(iun,*) '#declare COIL = texture {'
          write(iun,*) 'pigment { color rgb<1.0, 1.0, 1.0> }'
          write(iun,*)
     &       'finish { ambient 0.4 diffuse 0.4 specular 0.9}'
          write(iun,*) '}'

          call ribbon(0,2,1,iun,0,0,0)
          call ribbon(0,2,2,iun,0,0,0)
          call ribbon(0,2,3,iun,0,0,0)
          call ribbon(0,2,4,iun,0,0,0)
          call ribbon(1,2,1,iun,0,0,0)
          call ribbon(1,2,2,iun,0,0,0)
          call ribbon(1,2,3,iun,0,0,0)
          call ribbon(1,2,4,iun,0,0,0)
          call ribbon(3,2,1,iun,0,0,0)
          call ribbon(2,2,1,iun,0,0,0)

         endif

      endif

      ihand = 1

      do i=1,iatoms
         lring(i) = 0
      end do

      return
      end

      integer function mmtyp(iat,ian,idochg,ianz,iaton,iconn)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxt=14)
      character*4 atype
      common /atypes/ ihbt(mxt),atype(mxt)
      logical smrng
      integer*2 ictyp
      dimension icnn(mxcon),iring(4),ianz(*),iaton(*),iconn(mxcon+1,*)

      call ispn(isp,iat,irng,idochg,0)
      ihb = ihbt(isp)
      ict = ictyp(iat,ian,idochg,ianz,iconn)

      ibnds = 0
      ih0 = 0
      ic0 = 0
      in0 = 0
      io0 = 0
      ip0 = 0
      is0 = 0
      ili = 0
      do i=1,iconn(1,iat)
         if (iconn(i+1,iat).gt.0) then
            ibnds = ibnds + 1
            icnn(ibnds) = iconn(i+1,iat)
            ian0 = ianz(iconn(i+1,iat))
            if (ian0.eq.1) then
               ih0 = ih0 + 1
            elseif (ian0.eq.3) then
               ili = ili + 1
            elseif (ian0.eq.6) then
               ic0 = ic0 + 1
            elseif (ian0.eq.7) then
               in0 = in0 + 1
            elseif (ian0.eq.8) then
               io0 = io0 + 1
            elseif (ian0.eq.15) then
               ip0 = ip0 + 1
            elseif (ian0.eq.16) then
               is0 = is0 + 1
            endif
         endif
      end do

      mmtyp = 0

      if (ian.eq.1) then

         mmtyp = 5

         if (ibnds.eq.1) then

            ian1 = ianz(icnn(1))

            if (ian1.eq.6) then

               mmtyp = 5
CNF Add the case of acetylene -> mmtyp = 124
               if (icred(icnn(1),idum1,idum2,ianz,iconn)
     &             .eq.2) mmtyp = 124

            elseif (ian1.eq.7) then

               mmtyp = 23
               if (icred(icnn(1),idum1,idum2,ianz,iconn)
     &             .eq.4) mmtyp = 48
               call ispn(irs,icnn(1),irng,idochg,0)
               if (irs.eq.8) mmtyp = 28

            elseif (ian1.eq.8) then

CNF Add the case of phenol/enol -> mmtyp = 73
               mmtyp = 21

               if (icred(icnn(1),idum1,idum2,ianz,iconn)
     &             .eq.2) then
                  do i=1,iconn(1,icnn(1))
                     l = abs(iconn(i+1,icnn(1)))
                     if (ianz(l).eq.6) then
                        if (icred(l,idum1,idum2,ianz,iconn)
     &                      .eq.3) then
                           io = 0
                           ic = 0
                           do j=1,iconn(1,l)
                              k = abs(iconn(j+1,l))
                              if (ianz(k).eq.8) io = io + 1
                              if (ianz(k).eq.6) ic = ic + 1
                           end do
                           if (io.eq.2) mmtyp = 24
                           if (ictyp(l,ianz(l),idochg,ianz,iconn)
     &                          .eq.16) mmtyp = 24
                           if (ic.ge.1) mmtyp = 73
                        endif
                     endif
                  end do
               endif
            elseif (ian1.eq.16) then
               mmtyp = 44
            endif
         endif
      elseif (ian.eq.6) then
         mmtyp = 1
         if (ihb.eq.1) then
            mmtyp = 4
         elseif(ihb.eq.2.or.ihb.eq.4) then
            mmtyp = 2
            do i=1,ibnds
               if (ianz(icnn(i)).eq.8) then
                  if (icred(icnn(i),idum1,idum2,ianz,iconn)
     &                .eq.1) mmtyp = 3
               endif
            end do
            if (smrng(iat,iring,nring,ianz,iaton,iconn)) then
               if (nring.eq.3) mmtyp = 38
               if (nring.eq.4) mmtyp = 57
            endif
         elseif(ihb.eq.3) then
            mmtyp = 1
            if (smrng(iat,iring,nring,ianz,iaton,iconn)) then
               if (nring.eq.3) mmtyp = 22
               if (nring.eq.4) mmtyp = 56
            endif
         endif
      elseif (ian.eq.7) then
         mmtyp = 8
         if (ihb.eq.1) then
            mmtyp = 10
         elseif(ihb.eq.2.or.ihb.eq.4) then
            mmtyp = 37
            if (isp.eq.8) mmtyp = 9
         elseif(ihb.eq.3) then
            mmtyp = 8
            if (isp.eq.5) mmtyp = 39
            if (isp.eq.8) mmtyp = 9
            if (ict.eq.34) mmtyp = 40
            if (ict.eq.38) mmtyp = 46
            if (ili.gt.0) mmtyp = 164
         endif
      elseif (ian.eq.8) then
         mmtyp = 6
         if (ihb.eq.2) then
            mmtyp = 7
            if (ibnds.eq.1.and.ic0.eq.1) then
               do i=1,iconn(1,icnn(1))
                  iat2 = iconn(i+1,icnn(1))
                  if (ianz(iat2).eq.8.and.iat2.ne.iat) then
                     icc = 0
                     do j=1,iconn(1,iat2)
                       iat3 = iconn(j+1,iat2)
                       if (iat3.gt.0) then
                          if (ianz(iat3).eq.6) icc = icc + 1
                       endif
                     end do
                     if (icc.eq.1) mmtyp = 77
                     if (icc.eq.2) mmtyp = 78
                  endif
               end do
            endif
            if (ict.eq.52) mmtyp = 41
            if (isp.eq.10.and.ianz(icnn(1)).eq.6) mmtyp = 47
         elseif(ihb.eq.3) then
            mmtyp = 6
            if (ic0.ge.1) then
               do i=1,ibnds
                  if (ianz(icnn(i)).eq.6.and.icnn(i).ne.iat) then
                     do j=1,iconn(1,icnn(i))
                        jj = iconn(1+j,icnn(i))
                        if (ianz(iabs(jj)).eq.8.and.iabs(jj).ne.iat)
     &                  then
                          if (icred(icnn(i),idum1,idum2,ianz,iconn)
     &                        .eq.3) mmtyp = 75 
                        endif
                     end do
                  endif
               end do
            endif
            if (ict.eq.52) mmtyp = 41
            if (smrng(iat,iring,nring,ianz,iaton,iconn)) then
               if (nring.eq.3) mmtyp = 49
            endif
         endif
      elseif (ian.eq.15) then
         mmtyp = 60
         if (ibnds.eq.3) mmtyp = 25
      elseif (ian.eq.16) then
         mmtyp = 15
         if (isp.eq.13) mmtyp = 17
         if (isp.eq.14) mmtyp = 18
         if (ict.eq.72.or.ibnds.eq.1) mmtyp = 42
         if (ibnds.eq.2.and.is0.eq.1) mmtyp = 104
      else
         if (ian.eq.2) mmtyp = 51
         if (ian.eq.3) mmtyp = 163
         if (ian.eq.5) then
            mmtyp = 27
            if (ibnds.eq.3) mmtyp = 26
         endif
         if (ian.eq.9)  mmtyp = 11
         if (ian.eq.10) mmtyp = 52
         if (ian.eq.12) mmtyp = 59
         if (ian.eq.14) mmtyp = 19
         if (ian.eq.17) mmtyp = 12
         if (ian.eq.18) mmtyp = 53
         if (ian.eq.20) mmtyp = 125
         if (ian.eq.26) mmtyp = 61
         if (ian.eq.27) mmtyp = 65
         if (ian.eq.28) mmtyp = 63
         if (ian.eq.32) mmtyp = 31
         if (ian.eq.34) mmtyp = 34
         if (ian.eq.35) mmtyp = 13
         if (ian.eq.36) mmtyp = 54
         if (ian.eq.38) mmtyp = 126
         if (ian.eq.50) mmtyp = 32
         if (ian.eq.52) mmtyp = 35
         if (ian.eq.53) mmtyp = 14
         if (ian.eq.54) mmtyp = 55
         if (ian.ge.56.and.ian.le.71) mmtyp = ian + 71
         if (ian.eq.82) mmtyp = 33
      endif

      if (mmtyp.eq.0) print*,'no mm3 type for this atom'

      return
      end

      logical function smrng(i,iring,nring,ianz,iaton,iconn)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      logical ocnos
      dimension iring(4),ianz(*),iaton(*),iconn(mxcon+1,*)

c finds 3,4 membered rings

      smrng = .false.

      if (ocnos(i,ianz,iaton)) then
         do j=1,iconn(1,i)
            jj = abs(iconn(j+1,i))
            if (ocnos(jj,ianz,iaton)) then
               do k=1,iconn(1,jj)
                  kk = abs(iconn(k+1,jj))
                  if (ocnos(kk,ianz,iaton).and.kk.ne.i) then
                     do l=1,iconn(1,kk)
                        ll = abs(iconn(l+1,kk))
                        if (ocnos(ll,ianz,iaton).and.ll.ne.jj) then
                           if (ll.eq.i) then
                              nring = 3
                              iring(1) = i
                              iring(2) = jj
                              iring(3) = kk
                              smrng = .true.
                              return
                           else
                              do m=1,iconn(1,ll)
                                 mm = abs(iconn(m+1,ll))
                                 if (ocnos(mm,ianz,iaton).and.mm.ne.kk) 
     &                           then
                                    if (mm.eq.i) then
                                       nring = 4
                                       iring(1) = i
                                       iring(2) = jj
                                       iring(3) = kk
                                       iring(4) = ll
                                       smrng = .true.
                                       return
                                    endif
                                 endif
                              end do
                           endif
                        endif
                     end do
                  endif
               end do
            endif
         end do
      endif

      return
      end

      subroutine gtheat(igttnk,iheat,heat)
      implicit double precision (a-h,o-z), integer ( i-n)
      integer getlin
      character*137 str
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character*137 line
      common /curlin/ line

      igttnk = 1
      idebug = 0

      if (idebug.eq.1) print*,'subroutine gtheat'

      call search(line,'[AMBFOR]',istat)

      if (istat.eq.1) then

         if (getlin(0).eq.1) idum = 1

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.2) goto 100
         
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.1) then
            if (nstr.eq.4) then
               if (icdex(str,"emin").ne.0) then
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.eq.2) then
                       iheat = 1
                       heat = dble(itype)
                   endif
                   if (ktype.eq.3) then
                       iheat = 1
                       heat = rtype
                   endif
               endif
            endif
         endif
      else
         goto 100
      endif
 
      return

100   igttnk = 0
      return
      end

      subroutine gettyd(ires,iat1,iat2,ityp,ipdbt,
     &                  ianz,iresid,iamino,icalf,ncalf)
      implicit double precision (a-h,o-z)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxgff=72)
      parameter (mxamb=1590)
      parameter (mxrsn=19)
      parameter (mxcat=79)
      parameter (numcal=50000)
      common /ambtmp/ ianamb(mxamb),iangff(mxgff),itpdbt(mxamb),
     &                iamca(2,mxcat),iamop(2,mxrsn)
      integer*2 ipdbt,ityp
      dimension icalf(6,*)
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension ianz(*),iresid(*),ipdbt(*),ityp(*),iamino(*)

      do i=iat1,iat2
         iresid(i) = ires
      end do

      if (ires.le.0) return

      do i=1,mxsym
          ipdb(i) = 0
      end do

      do i=1,mxhsym*3
          ihpdb(i) = 0
      end do

      iam = 0

      do i=iat1,iat2

         it = ityp(i)
         if (it.ne.0.and.it.le.mxamb) then
             if (ianz(i).eq.1) then

                ip = itpdbt(it)

                if (ip.ne.0) then
                   ipinc = 1
                   if (it.eq.117.or.it.eq.119.or.it.eq.132.or.it.eq.134)
     &                ipinc = 3
                   if (it.eq.303) then
                      if (ihpdb(ip).eq.0) then
                         ipdbt(i) =  ip
                         ihpdb(ip) = i
                      elseif (ihpdb(ip+1).eq.0) then
                         ipdbt(i) =  ip+1
                         ihpdb(ip+1) = i
                      elseif (ihpdb(ip+3).eq.0) then
                         ipdbt(i) =  ip+3
                         ihpdb(ip+3) = i
                      else
                         ipdbt(i) =  ip+4
                         ihpdb(ip+4) = i
                      endif
                   else
                      if (ihpdb(ip).eq.0) then
                         ipdbt(i) =  ip
                         ihpdb(ip) = i
                      elseif (ihpdb(ip+ipinc).eq.0) then
                         ipdbt(i) =  ip+ipinc
                         ihpdb(ip+ipinc) = i
                      else
                         ipdbt(i) =  ip+2
                         ihpdb(ip+2) = i
                      endif
                   endif
                endif

             else

                ip = itpdbt(it)
                if (ip.ne.0) then
                   if (ip.eq.4) then
                      if (ipdb(4).eq.0) then
                         ipdbt(i) =  4
                         ipdb(4) = i
                      elseif (ipdb(38).eq.0) then
                         ipdbt(i) = 38
                         ipdb(38) = i
                      endif
                   else
                      if (ipdb(ip).eq.0) then
                         ipdbt(i) =  ip
                         ipdb(ip) = i
                      elseif (ipdb(ip+1).eq.0) then
                         ipdbt(i) =  ip+1
                         ipdb(ip+1) = i
                      endif
                   endif
                endif

                if (ipdbt(i).eq.2) then
                   do j=1,mxcat
                      if (iamca(1,j).eq.it) iam = iamca(2,j)
                   end do
                elseif (ipdbt(i).eq.46) then
                   do j=1,mxrsn
                      if (iamop(1,j).eq.it) iam = iamop(2,j)
                   end do
                   if (iam.eq.0) then
                      inc = 0
                      if (it.eq.1244.or.it.eq.1232) inc = 1
                      if (it.eq.1246.or.it.eq.1234) inc = 1
                      it1 = ityp(i+inc)
                      do j=1,mxrsn
                         if (iamop(1,j)+inc.eq.it1) iam = iamop(2,j)
                      end do
                   endif
                   
                endif

             endif
         endif

      end do

      iamino(ires) = iam
  
      if (iam.gt.0.and.iam.le.20) then

         icalf(1,ires) = ipdb(2)
         icalf(2,ires) = ipdb(1)
         icalf(3,ires) = ipdb(3)
         icalf(4,ires) = ihpdb(1)

      elseif (iam.gt.23) then

         icalf(1,ires) = ipdb(43)
         icalf(2,ires) = ipdb(46)
         icalf(3,ires) = ipdb(47)
         icalf(4,ires) = ipdb(48)
         icalf(5,ires) = ipdb(50)
         icalf(6,ires) = ipdb(51)

c         if (ipdb(43).eq.0.and.ipdb(46).ne.0.and.ires.le.numcal)
c     &       icalf(1,ires) = ipdb(46)

         if (ipdb(43).eq.0.and.ires.le.numcal) then
             if (ihpdb(118).ne.0) then
                 icalf(1,ires) = ihpdb(118)
             elseif (ipdb(46).ne.0) then
                 icalf(1,ires) = ipdb(46)
             endif
         endif

      endif

      ncalf = ires

      return
      end

      subroutine gettnd(igttnk,idebug,ipdbon,iffset,iheat,heat,
     &                  ianz,iconn,iatclr,ityp,coo,q,
     &                  isurf,issdon,iclon,ichx,ishoh,nspg,nat,norg,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,o-z), integer ( i-n)
      integer getlin
      character*137 str, tstr
      character*2 catom, catomt, tolowf,iel, ggstr
      parameter (numatm=2000)
      parameter (mxcon=10)
      parameter (numcal=50000)
      parameter (maxsym=108)
      parameter (mxrsn=19)
      parameter (mxgff=72)
      parameter (mxamb=1590)
      parameter (mxcat=79)
      common /athlp/ iatoms, mxnat
      integer*2 ityp
      common /types/ iff
      common /zmfrst/ ihaszm, nz, mxzat
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character*137 line,tline
      common /curlin/ line
      character*80 tnknm
      integer tnkbg,tnkit,tnkarc,tnkarf,tnkprg
      common /tnkopt/ rmsgrd,tnkbg,tnkit,tnkarc,tnkarf,tnkprg,tnknm,icst
      common /tnkpro/ iresrd
      common /pbc/ abc(3),ibox,icell,igfmap
      common /ambtmp/ ianamb(mxamb),iangff(mxgff),itpdbt(mxamb),
     &                iamca(2,mxcat),iamop(2,mxrsn)
      logical first,dall,box,gnreal
      dimension iel(maxsym)
      dimension coo(3,*),q(*),ianz(*),iconn(mxcon+1,*),ityp(*),
     &          isurf(*),iatclr(*),v1(3),v2(3),v3(3)
      dimension ires(numcal),ibeg(numcal),iend(numcal)

      data iangff /
     & 99, 6, 6, 6, 6, 6, 
     &  6, 6, 6, 6, 6, 6, 
     &  6, 6, 6, 6, 6, 6, 
     &  1, 1, 1, 1, 1, 1, 
     &  1, 1, 1, 1, 1, 1,
     &  1, 9,17,35,53, 7, 
     &  7, 7, 7, 7, 7, 7, 
     &  7, 7, 7, 7, 7, 7, 
     &  8, 8, 8, 8,15,15, 
     & 15,15,15,15,15,15,
     & 15,15,15,16,16,16, 
     & 16,16,16,16,16,6/

      data (ianamb(i),i=1,60) /
     &   7, 6, 6, 1, 8, 1,
     &   7, 6, 6, 1, 8, 1,
     &   6, 1, 7, 6, 6, 1,
     &   8, 1, 6, 1, 6, 1,
     &   6, 1, 7, 6, 6, 1,
     &   8, 1, 6, 1, 6, 1,
     &   6, 1, 6, 1, 7, 6,
     &   6, 1, 8, 1, 6, 1,
     &   6, 1, 6, 1, 6, 1,
     &   7, 6, 6, 1, 8, 1/
      data (ianamb(i),i=61,120) /
     &   6, 1, 8, 1, 7, 6,
     &   6, 1, 8, 1, 6, 1,
     &   8, 1, 6, 1, 7, 6,
     &   6, 1, 8, 1, 6, 1,
     &  16, 1, 7, 6, 6, 1,
     &   8, 1, 6, 1,16, 7,
     &   6, 6, 8, 1, 6, 1,
     &   6, 1, 6, 1, 7, 6,
     &   6, 1, 8, 1, 6, 1,
     &   6, 6, 1, 6, 1, 6/
      data (ianamb(i),i=121,180) /
     &   1, 7, 6, 6, 1, 8,
     &   1, 6, 1, 6, 6, 1,
     &   6, 1, 6, 8, 1, 7,
     &   6, 6, 1, 8, 1, 6,
     &   1, 6, 6, 1, 6, 7,
     &   1, 6, 6, 1, 6, 1,
     &   6, 1, 6, 1, 7, 6,
     &   6, 1, 8, 1, 6, 1,
     &   6, 7, 1, 6, 1, 6,
     &   1, 7, 1, 7, 6, 6/
      data (ianamb(i),i=181,240) /
     &   1, 8, 1, 6, 1, 6,
     &   7, 1, 6, 1, 6, 1,
     &   7, 7, 6, 6, 1, 8,
     &   1, 6, 1, 6, 7, 6,
     &   1, 6, 1, 7, 1, 7,
     &   6, 6, 1, 8, 1, 6,
     &   1, 6, 8, 7, 6, 6,
     &   1, 8, 1, 6, 1, 6,
     &   8, 7, 1, 7, 6, 6,
     &   1, 8, 1, 6, 1, 6/
      data (ianamb(i),i=241,300) /
     &   1, 6, 8, 7, 6, 6,
     &   1, 8, 1, 6, 1, 6,
     &   1, 6, 8, 7, 1, 7,
     &   6, 6, 1, 8, 1, 6,
     &   1, 6, 1,16, 6, 1,
     &   7, 6, 6, 1, 8, 1,
     &   6, 1, 6, 1, 6, 1,
     &   6, 1, 7, 1, 7, 6,
     &   6, 1, 8, 1, 6, 1,
     &   6, 1, 6, 1, 7, 1/
      data (ianamb(i),i=301,360) /
     &   6, 7, 1, 7, 6, 6,
     &   1, 8, 1, 6, 1, 6,
     &   1, 6, 1, 7, 1, 7,
     &   6, 6, 1, 8, 6, 1,
     &   7, 6, 6, 1, 8, 1,
     &   6, 1, 6, 1, 6, 8,
     &   6, 1, 8, 6, 1, 6,
     &   8, 7, 1, 7, 1, 6,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1,8/
      data (ianamb(i),i=361,420) /
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 6, 1, 7, 6, 6,
     &   1, 8, 1, 7, 6, 6/
      data (ianamb(i),i=421,480) /
     &   1, 8, 1, 7, 6, 6,
     &   1, 8, 1, 7, 6, 6,
     &   1, 8, 1, 7, 6, 6,
     &   1, 8, 1, 7, 6, 6,
     &   1, 8, 1, 7, 6, 6,
     &   1, 8, 1, 7, 6, 6,
     &   1, 8, 1, 7, 6, 6,
     &   1, 8, 1, 7, 6, 6,
     &   1, 8, 1, 7, 6, 6,
     &   1, 8, 1, 7, 6, 6/
      data (ianamb(i),i=481,540) /
     &   1, 8, 1, 7, 6, 6,
     &   1, 8, 1, 7, 6, 6,
     &   1, 8, 1, 7, 6, 6,
     &   1, 8, 7, 6, 6, 1,
     &   8, 1, 7, 6, 6, 1,
     &   8, 1, 7, 6, 6, 1,
     &   8, 1, 7, 6, 6, 1,
     &   8, 1, 7, 6, 6, 1,
     &   8, 1, 7, 6, 6, 1,
     &   8, 1, 7, 6, 6, 1/
      data (ianamb(i),i=541,600) /
     &   8, 1, 7, 6, 6, 1,
     &   8, 1, 7, 6, 6, 1,
     &   8, 1, 7, 6, 6, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8/
      data (ianamb(i),i=601,659) /
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   1, 7, 6, 6, 1, 8,
     &   8, 1, 3,11,19,37,
     &  55,12,20,30,17/
      data (ianamb(i),i=660,671) /
     &   7, 6, 6, 1, 8, 1,
     &   6, 1, 6, 8, 8, 1/
      data (ianamb(i),i=672,685) /
     &   7, 6, 6, 1, 8, 1,
     &   6, 1, 6, 1, 6, 8,
     &   8, 1/
      data (ianamb(i),i=686,701) /
     &   7, 6, 6, 1, 8, 1,
     &   6, 1, 6, 1, 6, 1,
     &   6, 1, 7, 1/

      data (ianamb(i),i=702,1000) /299*0/

      data (ianamb(i),i=1001,1253) /
     & 8,6,1,1,6,1,8,6,1,6,1,6,1,8,1,8,7,6,6,7,
     & 6,7,6,7,6,1,7,1,1,1,8,6,1,1,6,1,8,6,1,6,
     & 1,6,1,8,1,8,7,6,6,7,6,7,6,7,6,1,7,1,1,8,
     & 1,8,6,1,1,6,1,8,6,1,6,1,6,1,8,1,8,7,6,7,
     & 6,6,6,8,7,1,1,1,1,8,6,1,1,6,1,8,6,1,6,1,
     & 6,1,8,1,8,7,6,7,6,6,6,8,1,8,1,1,8,6,1,1,
     & 6,1,8,6,1,6,1,6,1,1,8,7,6,6,7,6,7,6,7,6,
     & 1,7,1,1,1,8,6,1,1,6,1,8,6,1,6,1,6,1,1,8,
     & 7,6,6,7,6,7,6,7,6,1,7,1,1,8,1,8,6,1,1,6,
     & 1,8,6,1,6,1,6,1,1,8,7,6,7,6,6,6,8,7,1,1,
     & 1,1,8,6,1,1,6,1,8,6,1,6,1,6,1,1,8,7,6,7,
     & 6,6,6,8,1,8,6,1,1,15,8,8,1,8,15,8,8,1,8,15,
     & 8,15,8,8,1,8,15,8,8,1,8,15,8/

      data (ianamb(i),i=1254,1590) /
     & 8,6,1,1,6,1,8,6,1,6,1,6,1,8,1,8,7,6,6,7,
     & 6,7,6,7,6,1,7,1,1,6,1,8,6,1,6,1,8,6,1,6,
     & 1,6,1,8,1,8,7,6,7,6,6,6,1,8,7,1,6,1,8,6,
     & 1,6,1,8,6,1,6,1,6,1,8,8,7,6,8,7,6,7,1,6,
     & 1,6,1,6,1,8,6,1,6,1,8,6,1,6,1,6,1,8,1,8,
     & 7,6,6,7,6,1,7,6,7,1,6,8,6,1,6,1,8,6,1,6,
     & 1,8,6,1,6,1,6,1,8,1,8,7,6,6,7,6,7,6,7,6,
     & 8,1,1,7,6,1,8,6,1,6,1,8,6,1,6,1,6,1,8,1,
     & 8,7,6,6,7,6,7,6,7,6,8,1,1,7,1,6,1,8,6,1,
     & 6,1,8,6,1,6,1,6,1,8,8,7,6,6,7,6,7,6,7,6,
     & 8,1,1,7,1,6,1,8,6,1,6,1,8,6,1,6,1,6,1,8,
     & 1,8,7,6,6,7,6,7,6,7,6,8,7,6,6,1,6,1,6,1,
     & 6,1,6,1,6,1,6,8,8,6,1,7,1,6,8,8,6,1,8,6,
     & 1,6,1,8,6,1,6,1,6,1,8,1,8,7,6,7,6,6,6,8,
     & 8,1,1,1,8,6,1,6,1,8,6,1,6,1,6,1,8,1,8,7,
     & 6,7,6,6,6,8,8,1,1,6,1,8,6,1,6,1,8,6,1,6,
     & 1,6,1,8,1,8,6,6,7,6,7,6,8,8,1,1,1/

      data (itpdbt(i),i=1,150) /
     & 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 5, 7, 1, 2,
     & 3, 1, 4, 4, 5, 7, 7, 13, 8, 16, 1, 2, 3, 1, 4,
     & 4, 5, 7, 6, 10, 10, 22, 11, 25, 1, 2, 3, 1, 4, 4,
     & 5, 7, 7, 13, 8, 16, 10, 22, 1, 2, 3, 1, 4, 4, 5,
     & 7, 31, 10, 1, 2, 3, 1, 4, 4, 5, 7, 32, 13, 8, 16,
     & 1, 2, 3, 1, 4, 4, 5, 7, 37, 10, 1, 2, 3, 1, 4,
     & 4, 5, 7, 37, 1, 2, 3, 4, 4, 5, 7, 6, 10, 9, 19,
     & 1, 2, 3, 1, 4, 4, 5, 7, 6, 10, 22, 13, 31, 17, 40,
     & 1, 2, 3, 1, 4, 4, 5, 7, 6, 10, 22, 13, 31, 17, 33,
     & 52, 1, 2, 3, 1, 4, 4, 5, 7, 6, 10, 22, 11, 23/

      data (itpdbt(i),i=151,300) /
     & 31, 14, 15, 37, 18, 46, 19, 49, 16, 58, 1, 2, 3, 1, 4,
     & 4, 5, 7, 6, 20, 22, 11, 25, 13, 31, 24, 34, 1, 2, 3,
     & 1, 4, 4, 5, 7, 6, 20, 22, 11, 25, 13, 31, 24, 1, 2,
     & 3, 1, 4, 4, 5, 7, 6, 20, 11, 25, 13, 31, 24, 34, 1,
     & 2, 3, 1, 4, 4, 5, 7, 6, 29, 1, 2, 3, 1, 4, 4,
     & 5, 7, 6, 29, 21, 25, 1, 2, 3, 1, 4, 4, 5, 7, 6,
     & 10, 9, 34, 1, 2, 3, 1, 4, 4, 5, 7, 6, 10, 9, 34,
     & 24, 34, 1, 2, 3, 1, 4, 4, 5, 7, 6, 10, 36, 12, 28,
     & 1, 2, 3, 1, 4, 4, 5, 7, 6, 10, 9, 19, 12, 28, 27,
     & 40, 1, 2, 3, 1, 4, 4, 5, 7, 6, 10, 9, 19, 22, 28/

      data (itpdbt(i),i=301,450) /
     & 17, 25, 55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 1, 2, 3, 1, 4, 5, 9, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4,
     & 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2,
     & 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4,
     & 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2,
     & 3, 1, 4, 4, 9, 19, 1, 2, 3, 1, 4, 4, 1, 2, 3,
     & 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4,
     & 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3/

      data (itpdbt(i),i=451,600) /
     & 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4,
     & 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3,
     & 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4,
     & 1, 2, 3, 1, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1,
     & 4, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1,
     & 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1,
     & 4, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1,
     & 2, 3, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4,
     & 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2,
     & 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4/

      data (itpdbt(i),i=601,750) /
     & 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2,
     & 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4,
     & 4, 1, 2, 3, 1, 4, 4, 1, 2, 3, 1, 4, 4, 1, 2,
     & 3, 1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     & 2, 3, 1, 4, 4, 5, 8, 6, 29, 30, 25, 1, 2, 3, 1,
     & 4, 4, 5, 8, 6, 11, 9, 34, 35, 34, 1, 2, 3, 1, 4,
     & 4, 5, 8, 6, 11, 9, 20, 12, 29, 27, 42, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/

      data (itpdbt(i),i=751,900) /
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/

      data (itpdbt(i),i=901,1050) /
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 46, 47, 61, 62, 48,
     & 64, 49, 54, 73, 50, 67, 52, 70, 53, 76, 51, 66, 56, 57, 65,
     & 59, 62, 55, 60, 58, 82, 64, 136, 137, 112, 46, 47, 61, 62, 48,
     & 64, 49, 54, 73, 50, 67, 52, 70, 53, 76, 51, 66, 56, 57, 65/

      data (itpdbt(i),i=1051,1200) /
     & 59, 62, 55, 60, 58, 121, 61, 124, 125, 69, 112, 46, 47, 61, 62,
     & 48, 64, 49, 54, 73, 50, 67, 52, 70, 53, 76, 51, 60, 55, 62,
     & 56, 57, 58, 67, 63, 130, 131, 100, 103, 46, 47, 61, 62, 48, 64,
     & 49, 54, 73, 50, 67, 52, 70, 53, 76, 51, 60, 55, 62, 56, 57,
     & 58, 67, 127, 68, 100, 103, 46, 47, 61, 62, 48, 64, 49, 54, 73,
     & 50, 67, 52, 70, 71, 51, 66, 56, 57, 65, 59, 62, 55, 60, 58,
     & 82, 64, 136, 137, 112, 46, 47, 61, 62, 48, 64, 49, 54, 73, 50,
     & 67, 52, 70, 71, 51, 66, 56, 57, 65, 59, 62, 55, 60, 58, 121,
     & 61, 124, 125, 69, 112, 46, 47, 61, 62, 48, 64, 49, 54, 73, 50,
     & 67, 52, 70, 71, 51, 60, 55, 62, 56, 57, 58, 67, 63, 130, 131/

      data (itpdbt(i),i=1201,1253) /
     & 100, 103, 46, 47, 61, 62, 48, 64, 49, 54, 73, 50, 67, 52, 70,
     & 71, 51, 60, 55, 62, 56, 57, 58, 67, 127, 68, 83, 160, 103, 43,
     & 44, 46, 118, 46, 43, 44, 51, 115, 51, 43, 44, 43, 44, 46, 118,
     & 46, 43, 44, 51, 115, 51, 43, 44/

c 1MA

      data (itpdbt(i),i=1254,1284) /
     & 46, 47, 61, 62, 48, 64, 49, 54, 73, 50, 67, 52, 70, 53,
     & 76, 51, 66, 56, 57, 65, 59, 62, 55, 60, 58, 82, 64, 136,
     & 112, 79, 148/

c 5MC

      data (itpdbt(i),i=1285,1311) /
     & 46, 47, 61, 48, 64, 49, 54, 73, 50, 67, 52, 70, 53, 76,
     & 51, 60, 55, 62, 56, 57, 58, 103, 67, 63, 130, 83, 160/

c OMC

      data (itpdbt(i),i=1312,1338) /
     & 46, 47, 61, 48, 64, 49, 54, 73, 50, 67, 52, 70, 53, 51,
     & 60, 55, 67, 62, 56, 63, 130, 57, 100, 58, 103, 80, 151/
     
c 2MG

      data (itpdbt(i),i=1339,1369) /
     & 46, 47, 61, 48, 64, 49, 54, 73, 50, 67, 52, 70, 53, 76,
     & 51, 66, 56, 57, 65, 59, 112, 62, 55, 60, 79, 58, 69, 61,
     & 124, 80, 151/

c M2G

      data (itpdbt(i),i=1370,1399) /
     & 46, 47, 61, 48, 64, 49, 54, 73, 50, 67, 52, 70, 53, 76,
     & 51, 66, 56, 57, 65, 59, 62, 55, 60, 58, 69, 121, 112, 61,
     & 79, 148/
c 7MG

      data (itpdbt(i),i=1400,1430) /
     & 46, 47, 61, 48, 64, 49, 54, 73, 50, 67, 52, 70, 53, 76,
     & 51, 66, 56, 57, 65, 59, 62, 55, 60, 58, 69, 121, 112, 61,
     & 125, 85, 166/

c OMG

      data (itpdbt(i),i=1431,1460) /
     & 46, 47, 61, 48, 64, 49, 54, 73, 50, 67, 52, 70, 53, 51,
     & 66, 56, 57, 65, 59, 62, 55, 60, 58, 69, 121, 112, 61, 124,
     & 80, 151/

c YG

      data (itpdbt(i),i=1461,1511) /
     & 46, 47, 61, 48, 64, 49, 54, 73, 50, 67, 52, 70, 53, 76,
     & 51, 66, 56, 57, 65, 59, 62, 55, 60, 58, 69, 61, 90, 91,
     & 112, 88, 91, 89, 175, 92, 178, 93, 181, 94, 184, 95, 96,
     & 97, 98, 187, 99, 124, 100, 101, 102, 103, 190/
    
c H2U

      data (itpdbt(i),i=1512,1537) /
     & 46, 47, 61, 48, 64, 49, 54, 73, 50, 67, 52, 70, 53, 76,
     & 51, 60, 55, 62, 56, 57, 58, 67, 68, 127, 100, 103/ 

c 5MU
 
      data (itpdbt(i),i=1538,1564) /
     & 46, 47, 61, 48, 64, 49, 54, 73, 50, 67, 52, 70, 53, 76,
     & 51, 60, 55, 62, 56, 57, 58, 67, 68, 127, 103, 83, 160/

c PSU
 
      data (itpdbt(i),i=1565,1590) /
     & 46, 47, 61, 48, 64, 49, 54, 73, 50, 67, 52, 70, 53, 76,
     & 51, 57, 58, 60, 55, 62, 56, 67, 68, 103, 121, 127/

c CA amber types

       data ((iamca(i,j),i=1,2),j=1,mxcat) /
     & 2,1, 8,2, 16,7, 28,11, 42,6, 56,3, 66,5, 78,4, 
     & 88,4, 97,15, 108,18, 123,19, 139,20, 162,17, 179,17, 195,17, 
     & 211,9, 221,10, 233,13, 245,14, 259,8, 272,12, 288,16, 305,0, 
     & 319,2, 326,13, 351,1, 357,2, 363,7, 369,11, 375,6, 381,3,
     & 387,5, 393,4, 399,4, 405,15, 413,18, 419,19, 425,20, 431,17,
     & 437,17, 443,17, 449,9, 455,10, 461,13, 467,14, 473,8, 479,12,
     & 485,16, 491,0, 497,0, 502,1, 508,2, 514,7, 520,11, 526,6,
     & 532,3, 538,5, 544,4, 550,4, 556,15, 561,18, 567,19, 573,20,
     & 579,17, 585,17, 591,17, 597,9, 603,10, 609,13, 615,14, 621,8,
     & 627,12, 633,16, 639,0, 645,0, 661,9, 673,13, 687,12/

c o5* amber nucleotides

       data ((iamop(i,j),i=1,2),j=1,mxrsn) /
     & 1001,24, 1031,26, 1062,25, 1090,28, 
     & 1117,24, 1146,26, 1176,25, 1203,27,
     & 1254,29, 1285,30, 1312,31, 1339,33,
     & 1370,34, 1400,35, 1431,36, 1461,37,
     & 1512,40, 1538,41, 1565,42/

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


      igttnk = 1
      first = .true.
      dall = .false.
      ifndhn = 0
      igaff = 0
      iprot = 0
      toang = 0.52917706d0
      box = .false.

      if (idebug.eq.1) print*,'subroutine gettnk'

c get number of atoms

10    continue

      if (getlin(0).eq.1) then

         tline = line

         if (icdex(line,'[AMBMD]').ne.0) then
            call readel(line,1)
         endif

         if (icdex(line,'[AMBFOR]').ne.0) then
            igaff = 1
            first = .false.
            ic = icdex(line,'box')
            if (ic.ne.0) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)
               box = .true.
               do i=1,3
                  v1(i) = 0.0d0
                  v2(i) = 0.0d0
                  v3(i) = 0.0d0
               end do
               if (gnreal(abc,3,.false.)) then
                  ibox = 1
                  v1(1) = abc(1)/toang
                  v2(2) = abc(2)/toang
                  v3(3) = abc(3)/toang
               else
                  box = .false.
               endif
            endif

            ic = icdex(line,'cell')
            if (ic.ne.0) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)
               icell = 1

               if (gnreal(abc,3,.false.)) then
                  a = abc(1)
                  b = abc(2)
                  c = abc(3)
                  if (gnreal(abc,3,.false.)) then
                     alpha = abc(1)
                     beta  = abc(2)
                     gamma = abc(3)
                  else
                     icell = 0
                  endif
               else
                  icell = 0
               endif
               if (icell.eq.1) then
                  call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)
                  iclon = 1
                  ichx = 1
                  nspg = 1
               endif
            endif

            if (getlin(0).eq.1) idum = 1
         endif

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            iatoms = itype
            if (iatoms.gt.mxnat) then
               dall = .true.
               goto 100
            endif
         else
            goto 100
         endif

         
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.1) then
            if (nstr.eq.4) then
               if (icdex(str,"emin").ne.0) then
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.eq.2) then
                       iheat = 1
                       heat = dble(itype)
                   endif
                   if (ktype.eq.3) then
                       iheat = 1
                       heat = rtype
                   endif
               endif
            endif
         endif
      else
         goto 100
      endif

110   continue

      do iat=1,iatoms
         if (getlin(0).eq.1) then

           tline = line

           ktype = nxtwrd(str,nstr,itype,rtype)
           if (ktype.ne.2) goto 100
           if (igaff.eq.1) isurf(iat) = itype

           ktype = nxtwrd(str,nstr,itype,rtype)
           if (ktype.eq.1) then
              tstr = str
              nstrt = nstr
           else
              goto 100
           endif

           do i=1,3
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.3) then
                 coo(i,iat) = rtype
              else
                 goto 100
              endif
           end do

           ktype = nxtwrd(str,nstr,itype,rtype)
           if (ktype.eq.2) then
              ityp(iat) = itype
              if (itype.lt.0) then
                 igaff = 1
c this is a gaff atom type, so read in an extra charge field

                 ktype = nxtwrd(str,nstr,itype,rtype)
                 if (ktype.eq.3) then
                    q(iat) = rtype
                    ihasq = 1
                 else
                    goto 100
                 endif

              endif
           else
              goto 100
           endif
           nc = 0
           do i=1,mxcon
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.2) then
                 iconn(i+1,iat) = itype
                 nc = nc + 1
              endif
           end do
           iconn(1,iat) = nc

c determine elemnt

           if (igaff.eq.1) then

c AMBFOR atomnumber

              i1 = abs(ityp(iat))
              if (ityp(iat).lt.0) then
                 if (i1.gt.0.and.i1.le.mxgff) iatmp = iangff(i1)
              else
                 iprot = 1
                 if (i1.gt.0.and.i1.le.mxamb) iatmp = ianamb(i1)
              endif

           else

c other force fields atomnumber

              if (nstrt.eq.1) then
                 catomt(1:1) = tstr(1:1)
                 catomt(2:2) = ' '
              else
                 catomt = tstr(1:2)
                 if (catomt(2:2).eq.'+'.or.catomt(2:2).eq.'-'
     &              .or.catomt(2:2).eq.'*') catomt(2:2) = ' '
              endif
              catom = tolowf(catomt)
              iatmp = 0

              if (first.and.iffset.le.1) then

                 do j=1,maxsym
                    if (catom .eq. iel(j)) iatmp = j - 1
                 end do
                 if (catom.eq.'xx') iatmp = 99
                 if (catom.eq.'lp') iatmp = 99

              else

                 if (nstrt.eq.3) then
                     if (tstr(1:3).eq.'SOD'.or.tstr(1:3).eq.'sod') 
     &                   iatmp = 11
                     if (tstr(1:3).eq.'CAL'.or.tstr(1:3).eq.'cal') 
     &                   iatmp = 20
                     if (tstr(1:3).eq.'Li+') iatmp = 3
                     if (tstr(1:3).eq.'Na+') iatmp = 11
                     if (tstr(1:3).eq.'Rb+') iatmp = 37
                     if (tstr(1:3).eq.'Cs+') iatmp = 55
                     if (tstr(1:3).eq.'Mg+') iatmp = 12
                     if (tstr(1:3).eq.'Ca+') iatmp = 20
                     if (tstr(1:3).eq.'Zn+') iatmp = 30
                     if (tstr(1:3).eq.'Cl+') iatmp = 17
                 elseif (nstrt.eq.2) then
                     if (tolowf(tstr(1:2)).eq.'mg') iatmp = 12
                     if (tolowf(tstr(1:2)).eq.'fe') iatmp = 26
                     if (tolowf(tstr(1:2)).eq.'zn') iatmp = 30
                     if (tolowf(tstr(1:2)).eq.'h1') ifndhn = 1
                     if (tolowf(tstr(1:2)).eq.'hn') ifndhn = 2
                     if (tolowf(tstr(1:2)).eq.'k+') iatmp = 19
                 endif

                 if (iatmp.eq.0) then
                     if (catom(1:1).eq.'h') iatmp = 1
                     if (catom(1:1).eq.'c') iatmp = 6
                     if (catom(1:1).eq.'n') iatmp = 7
                     if (catom(1:1).eq.'o') iatmp = 8
                     if (catom(1:1).eq.'f') iatmp = 9
                     if (catom(1:1).eq.'p') iatmp = 15
                     if (catom(1:1).eq.'s') iatmp = 16
                     if (catom(1:2).eq.'cl') iatmp = 17
                     if (catom(1:2).eq.'br') iatmp = 35
                     if (catom(1:1).eq.'i') iatmp = 53
                     if (catom(1:1).eq.'x') iatmp = 99
                 endif
              endif
              if (iatmp.le.0.or.iatmp.gt.maxsym-1) goto 100

           endif

           ianz(iat) = iatmp

        else
          print*,'Number of atoms read ',iat,
     &           ' less than specified in the header of file ',
     &           iatoms
          goto 100
        endif

      end do

      if (getlin(0).eq.1) then

         if (icdex(line,'[RESIDUES]').ne.0) then

             ihsres = 0
             do while (getlin(0).eq.1)

                if (icdex(line,'[AMBFOR]').ne.0) then
                   backspace(iun2)
                   goto 200
                endif

                if (iresrd.eq.0) then

                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.eq.2) then
                      ihsres = ihsres + 1
                      irs = itype
                      ires(ihsres) = irs
                   endif

                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.eq.2) then
                      istrt = itype
                      ibeg(ihsres) = istrt
                      if (ihsres.gt.1) iend(ihsres-1) = istrt - 1
                   endif

                endif

             end do

200          if (iresrd.eq.0) then

                iend(ihsres) = iatoms

                do i=1,ihsres
                   call gettyp(ires(i),ibeg(i),iend(i))
                end do

                call chkbrk

                call parsfn('Helix',5,1)
                call parsfn('Beta',4,1)
                call parsfn('RNA/DNA',7,1)
                call parsfn('Coil',4,1)

                ihet = 0
                iwat = 0
                iset = 0

                do i=1,ihsres
                   if (ires(i).lt.0) then
                      ihet = ihet + 1
                      if (iend(i)-ibeg(i)+1.eq.3) then
                         ib = ibeg(i)
                         it1 = ityp(ib)
                         it2 = ityp(ib+1)
                         it3 = ityp(ib+2)
                         if (it1.eq.649.and.it2.eq.650.and.it3.eq.650)
     &                      iwat = 1
                      endif
                      if (iwat.eq.1) then
                         if (iset.eq.0) then
                            ishoh = iabs(ires(i))
                            iset = 1
                         endif
                         call parsfn('HOH',3,1)
                      else
                         call parsfn('HET'//ggstr(ihet),5,1)
                      endif
                   endif
                end do

                iresrd = 1

             endif

         else

             backspace(iun2)

         endif


      endif


      call xyzcoo(0,1,0)
      call cooxyz(ianz,iatoms)

      ihaszm = 0
      issdon = 0

      if (iffset.eq.0) then
         if (first) then
            iff = 1
            if (igaff.eq.1) iff = 7
         else
            ipdbon = 1
            if (igaff.eq.1) then
               iff = 7
               if (iprot.ne.1) ipdbon = 0
            else
               if (ifndhn.eq.1) then
                  iff = 3
               elseif (ifndhn.eq.2) then
                  iff = 4
               else
                  iff = 2
               endif
            endif
         endif
         if (iff.eq.3.or.iff.eq.7) then
            do iat=1,iatoms
               call setchg(iat,1)
            end do
         endif
      endif

      if (icell.eq.1) then
         call cpmol2(nat,norg,xa,ya,yb,za,zb,zc,
     &              coo,ianz,iatclr,iconn)
         call addc(coo,ianz,iconn,iatclr,iatoms,xa,ya,yb,za,zb,zc)
      endif

      if (box) then
         iclon = 1
         call addtbx(v1,v2,v3)
      endif

      return

100   if (dall) then
         rewind iun2
         igttnk = -1
         return
      endif

      if (first.and.iffset.eq.0) then
         rewind iun2
         first = .false.
         ifndhn = 0
         goto 10
      else
         if (idebug.eq.1) print*,'ERROR:',tline
      endif
      igttnk = 0

      return
      end

      subroutine dotyd(icel,ianz,iaton,iatclr,iconn,iresid,
     &                 lwrit,lring,ityp,coo,qat,icont,
     &                 icalf,ncalf,ianf,islu,nchain,iamino,ishoh,
     &                 nat,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)

      parameter (mxel=100)
      parameter (numatm=2000)
      parameter (mxcon=10)
      parameter (mxt=14)
      parameter (mxppmf=16)
      parameter (mxlpmf=26)
      parameter (mxmol2=41)
      parameter (mxmm3=164)
      parameter (mxchtp=136)
      parameter (mxmsf=235)
      parameter (mxamb=1590)
      parameter (mxamo=201)
      parameter (mxsym=103)
      parameter (mxhsym=64)

      common /athlp/ iatoms, mxnat
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      integer*2 ityp
      common /types/ iff
      character*2  ppmf, lpmf
      character*5  mol2
      character*19 mm3
      character*20 chmtnk
      character*20 ambstr
      character*20 amostr
      character*4  chmsf
      common /ftypes/ihasl(11),mol2(mxmol2),mm3(mxmm3),chmtnk(mxchtp),
     &               chmsf(mxmsf),ambstr(mxamb),amostr(mxamo),
     &               ppmf(mxppmf),lpmf(mxlpmf)
      character*6 atmp
      character*4 atype
      common /atypes/ ihbt(mxt),atype(mxt)
      character*2 elemnt,atom
      common /elem/elemnt(mxel)
      integer dolabs,fancy,persp,shade,atcol,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      common /typoni/ ioniad
      integer*2 ictyp,igtyp
      logical ochg
      dimension rr(3,3),tr1(3)
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension coo(3,*),ianz(*),iaton(*),iatclr(*),iconn(mxcon+1,*),
     &          iresid(*),ityp(*),lwrit(*),lring(*),qat(*),icont(*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*)

      natoms = iatoms

      if (icel.eq.1) then

         if (ochg(idum,ianz)) then
            idochg = 1
         else
            idochg = 0
         endif

         natoms = nat
         nstor = mxnat-natoms

         do i=1,natoms

            do j=1,3
               coo(j,i) = coo(j,nstor+i)
            end do

            ianz(i) = ianz(nstor+i)

            do j=1,iconn(1,nstor+i)+1
               iconn(j,i) = iconn(j,nstor+i)
            end do

            iatclr(i) = iatclr(nstor+i)

         end do

      else
         if (ihasq.eq.1) idochg = 1
      endif

      do i=1,natoms
         lwrit(i) = 0
         lring(i) = 0
         icont(i) = iaton(i)
         iaton(i) = 2
         if (ioniad.eq.0) ityp(i) = 0
      end do

      if (icel.eq.1) then

         call setrr(alpha,beta,gamma,a,b,c,rr)
  
         do i=1,nat
            do k=1,3
               tr1(k) = trc(coo(1,i),rr,k)
            end do
            do k=1,3
               coo(k,i) = tr1(k)
            end do
         end do

      endif

      if (iff.eq.2) then

c Tinker Charmm
         do i=1,ncalf
             call getpdb(i,ipdb,ihpdb)
             call typeit(ipdb,iamino(i),ihpdb,1)
         end do

      elseif (iff.eq.3.or.iff.eq.7) then

         ifftmp = iff

c Tinker Amber

         call chkrna(ncalf,iamino)

         do i=1,ncalf
             call getpdb(i,ipdb,ihpdb)
             call typamb(ipdb,iamino(i),ihpdb,1)
         end do

         do i=1,nchain
c n-term cap
            call chkcap(ianf(i),1,iconn,ianz,iresid,qat,icalf,iamino)
c c-term cap
            call chkcap(islu(i),2,iconn,ianz,iresid,qat,icalf,iamino)
         end do

         iff = ifftmp
         do i=1,natoms
            if (iresid(i).le.-4) then
               if (iresid(i).eq.-ishoh) then
                  if (ianz(i).eq.8) ityp(i) = 649
                  if (ianz(i).eq.1) ityp(i) = 650
               else
                  ncnn = 0
                  do j=1,iconn(1,i)
                     if (iconn(1+j,i).gt.0) ncnn = ncnn + 1
                  end do

                  iset = 0
                  if (ianz(i).eq.3) then
                     ityp(i) = 651
                     qat(i) = 1.0d0
                     iset = 1
                  else if (ianz(i).eq.11) then
                     ityp(i) = 652
                     qat(i) = 1.0d0
                     iset = 1
                  else if (ianz(i).eq.19) then
                     ityp(i) = 653
                     qat(i) = 1.0d0
                     iset = 1
                  else if (ianz(i).eq.37) then
                     ityp(i) = 654
                     qat(i) = 1.0d0
                     iset = 1
                  else if (ianz(i).eq.55) then
                     ityp(i) = 655
                     qat(i) = 1.0d0
                     iset = 1
                  else if (ianz(i).eq.12) then
                     ityp(i) = 656
                     qat(i) = 2.0d0
                     iset = 1
                  else if (ianz(i).eq.20) then
                     ityp(i) = 657
                     qat(i) = 2.0d0
                     iset = 1
                  else if (ianz(i).eq.30) then
                     ityp(i) = 658
                     qat(i) = 2.0d0
                     iset = 1
                  else if (ianz(i).eq.17.and.ncnn.eq.0) then
                     ityp(i) = 659
                     qat(i) = -1.0d0
                     iset = 1
                  endif

                  if (iset.ne.1) then
                     if (iff.eq.7) then
                        ityp(i) = igtyp(i,ianz(i),idochg,ianz,iconn)
                        ityp(i) = -ityp(i)
                     endif
                  endif

               endif
            endif
         end do

      elseif (iff.eq.4) then

c Tinker Amoeba
         do i=1,ncalf
             call getpdb(i,ipdb,ihpdb)
             call typamo(ipdb,iamino(i),ihpdb,1)
         end do

      else

         do i=1,natoms
            ityp(i) = 0

            if (iff.eq.1) then
c Tinker MM3
                ityp(i) = mmtyp(i,ianz(i),0,ianz,iaton,iconn)

            elseif (iff.eq.5) then

c Sybyl Mol2
                ityp(i) = 1
                atom = elemnt(ianz(i))
                call ispn(irs,i,irng,idochg,0)
                atmp = atom//atype(irs)
                if (atmp(1:1).eq.' ') atmp(1:5) = atmp(2:6)

                do j=1,mxmol2
                   if (atmp(1:5).eq.mol2(j)) ityp(i) = j
                end do

            elseif (iff.eq.6) then

c Quanta Charmm
                ityp(i) = ictyp(i,ianz(i),idochg,ianz,iconn)

            elseif (iff.eq.8) then

c PMF scoring
c                call ipmtyp(ipmt,i,ianz(i),idochg)
c                ityp(i) = ipmt

            endif
         end do

      endif

      do i=1,iatoms
         iaton(i) = icont(i)
      end do

      if (icel.eq.1) call fdat(ifd,0,0,0,0,0)

c      do i=1,iatoms
c         iaton(i) = icont(i)
c      end do

      return
      end

      subroutine chkcap(irs,iterm,iconn,ianz,iresid,qat,
     &                  icalf,iamino)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      common /conrl/ ibnds,icnn(mxcon),io,in,ic,ih,ian1,ian2,ian3,ian4
      dimension iconn(mxcon+1,*),ianz(*),iresid(*),qat(*)
      dimension nace(6),qace(6),nme(6),qnme(6),qnh2(3)
      dimension icalf(6,*),iamino(*)
      data qace /0.59720,-0.56790,-0.36620,0.11230,0.11230,0.11230/
      data qnme /-0.41570,0.27190,-0.14900,0.09760,0.09760,0.09760/
      data qnh2 /-0.46300,0.23150,0.23150/

      if (iterm.eq.1) then
         ncterm = icalf(2,irs)
      else
         ncterm = icalf(3,irs)
      endif

      call getrcn(ncterm,iconn,ianz)

      if (iterm.eq.1) then
c N-term
c        establish carbonyl C of ACE

         do i=1,6
             nace(i) = 0
         end do

         iace = 0
         if ((ic.eq.2.and.iamino(irs).ne.15).or.
     &       (ic.eq.3.and.iamino(irs).eq.15)) then
             
             do i=1,ibnds
                if (ianz(icnn(i)).eq.6.and.iresid(icnn(i)).ne.irs)
     &            iace = icnn(i)
             end do
             if (iace.eq.0) return

         else
            return
         endif

         nace(1) = iace
         call getrcn(iace,iconn,ianz)
         if (io.eq.1) then
             do i=1,ibnds
                if (ianz(icnn(i)).eq.8) nace(2) = icnn(i)
             end do
         else
             return
         endif

         if (ic.eq.1) then
             do i=1,ibnds
                if (ianz(icnn(i)).eq.6) nace(3) = icnn(i)
             end do
         else
             return
         endif

         call getrcn(nace(3),iconn,ianz)
         if (ih.eq.3) then
             itel = 4
             do i=1,ibnds
                if (ianz(icnn(i)).eq.1) then
                    nace(itel) = icnn(i)
                    itel = itel + 1
                endif
             end do
         else
             return
         endif

c        check if all atoms ACE found

         do i=1,6
            if (nace(i).le.0) return
         end do

         do i=1,6
            qat(nace(i)) = qace(i)
         end do

      else

c C-term
c        establish n-methyl N of NME/NH2

         do i=1,6
             nme(i) = 0
         end do

         inme = 0

         if (in.eq.1) then
             
             do i=1,ibnds
                if (ianz(icnn(i)).eq.7.and.iresid(icnn(i)).ne.irs)
     &            inme = icnn(i)
             end do
             if (inme.eq.0) return

         else
            return
         endif

         nme(1) = inme
         inh2 = 0
         call getrcn(inme,iconn,ianz)
         if (ih.eq.1) then
             do i=1,ibnds
                if (ianz(icnn(i)).eq.1) nme(2) = icnn(i)
             end do
         else if (ih.eq.2) then
             inh2 = 1
             itel = 0
             do i=1,ibnds
                if (ianz(icnn(i)).eq.1) then
                   nme(2+itel) = icnn(i)
                   itel = itel + 1
                endif
             end do
         else
             return
         endif

         if (inh2.eq.1) then

c           this is NH2 cap

            do i=1,3
               if (nme(i).le.0) return
            end do

            do i=1,3
               qat(nme(i)) = qnh2(i)
            end do

            return
         endif

c        this is NME cap

         if (ic.eq.2) then
             do i=1,ibnds
                if (ianz(icnn(i)).eq.6.and.iresid(icnn(i)).ne.irs) 
     &             nme(3) = icnn(i)
             end do
         else
             return
         endif

         call getrcn(nme(3),iconn,ianz)
         if (ih.eq.3) then
             itel = 4
             do i=1,ibnds
                if (ianz(icnn(i)).eq.1) then
                    nme(itel) = icnn(i)
                    itel = itel + 1
                endif
             end do
         else
             return
         endif

c        check if all atoms NME found

         do i=1,6
            if (nme(i).le.0) return
         end do

         do i=1,6
            qat(nme(i)) = qnme(i)
         end do

      endif

      return
      end

      subroutine fixchg(ich,qat)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxel=100)
      common /athlp/ iatoms, mxnat
      dimension qat(*)

      addup = dble(ich)

      totch = 0.0d0
      do i=1,iatoms
          totch = totch + qat(i)
      end do 

      totch = totch - addup
      if (dabs(totch).gt.1.0d-7) then
           totch = totch/dble(iatoms)
           do i=1,iatoms
               qat(i) = qat(i) - totch
           end do
      endif

      return
      end
