      subroutine rdchd(idebug,iop,istdbd,iuseab,moddma,istat,
     &                 icssr,coo,iconn,ianz,iatclr,ityp,qat,
     &                 nat,icent,inorm,ncon,nspg,kz,ichx,icrtp,
     &                 nopr,ir,it,
     &                 xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
c      implicit double precision (a-h,o-z)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      integer*2 ir,it
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      integer*2 ityp
      common /types/ iff
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      parameter (mxel=100)
      character*2 elemnt
      common /elem/elemnt(mxel)
                                                                                
      character lstr*137
      character*4 atname,tspg
      character*2 tstr,tocapf
      character*3 spg
      character*8 refcod
      logical datlin,ochg,hascon
      dimension icon(8)
      dimension coo(3,*),ityp(*),qat(*),iconn(mxcon+1,*),
     &          ianz(*),iatclr(*),ir(3,3,192),it(3,192)
c
c     Read Chemx file
c
      istat = 1
      nptdum = 61
      if (icssr.ne.1) call rewfil

      toang = 0.52917706d0

c      a = 1.0d0
c      b = 1.0d0
c      c = 1.0d0
c      alpha = 90.0d0
c      beta = 90.0d0
c      gamma = 90.0d0

      hascon = .false.

      
      call nxtlin(lstr,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 100

      if (icdex(lstr,'[Molden Format]').ne.0) then
c check if not Molden format, not nice 
         istat = 0
         return
      endif

      if (icdex(lstr,'@<TRIPOS>').ne.0) then
         istat = 0
         return
      endif

      if (linlen(lstr).ge.62) then
         lstr = lstr(39:62)
         if (datlin(lstr)) then
            read(lstr,'(3f8.3)',end=100,err=100) at,bt,ct
         else
            goto 100
         endif
      else
         goto 100
      endif

      a = at
      b = bt
      c = ct

      call nxtlin(lstr,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 100
      if (linlen(lstr).ge.45) then
         if (linlen(lstr).ge.57) then
            tspg = lstr(55:58)
            if (tspg(1:1).eq.'=') then
               spg = tspg(2:4)
            else
               spg = tspg(1:3)
            endif
            read(spg,'(i3)',end=100,err=100) nspg
         else
            nspg = 1
         endif
         if (linlen(lstr).eq.46) then
            lstr = lstr(23:46)
         else
            lstr = lstr(22:45)
         endif
         if (datlin(lstr)) then
            read(lstr,'(3f8.3)',end=100,err=100) alphat,betat,gammat
         else
            goto 100
         endif
      else
         goto 100
      endif

      alpha = alphat
      beta = betat
      gamma = gammat

      call nxtlin(lstr,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 100
      if (linlen(lstr).ge.8) then
         lstr = lstr(1:8)
         if (datlin(lstr)) then
            read(lstr,'(2i4)',end=100,err=100) iatoms,inorm
         else
            goto 100
         endif
      else
         goto 100
      endif

      icent = 1
      nopr = 0
      kz = 0

      if (iatoms.lt.1) goto 100

      call redel(lstr,1)
      if (inorm.eq.0) then
C
C Set up matrix to convert fractional to Angstrom coordinates.
C
         call prcell(nspg,a,b,c,alpha,beta,gamma)
         call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)
         call cprot(nspg,nopr,icent,ir,it,.false.)
         if (idebug.eq.1) call prop(nopr,ir,it)
      else
         call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)
      endif
c
c     Atom loop
c
      nstor = mxnat-iatoms

1000  format(5x,a4,2x,3(f9.5,1x),8i4,f8.3,1x,i3)

      do i=1,iatoms
         iconn(1,i) = 0
         call nxtlin(lstr,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 100
         if (linlen(lstr).gt.81) then
            read(lstr,1000,
     &      end=100,err=100) atname,coo(1,i),coo(2,i),coo(3,i),
     &      (icon(j),j=1,8),qat(nstor+i),ityp(i)
            iff = 1
         elseif (linlen(lstr).gt.72) then
            read(lstr,1000,
     &      end=100,err=100) atname,coo(1,i),coo(2,i),coo(3,i),
     &      (icon(j),j=1,8),qat(nstor+i)
         elseif (linlen(lstr).gt.39) then
            read(lstr,1000,
     &      end=100,err=100) atname,coo(1,i),coo(2,i),coo(3,i),
     &      (icon(j),j=1,8)
            qat(nstor+i) = 0.0d0
         else
            read(lstr,1000,
     &      end=100,err=100) atname,coo(1,i),coo(2,i),coo(3,i)
            do j=1,8
               icon(j) = 0
            end do
            qat(nstor+i) = 0.0d0
         endif
         tstr = atname(1:2)
         j = ichar(tstr(2:2))
         if (j.lt.65.or.j.gt.122.or.(j.lt.97.and.j.gt.90)) then
            tstr(2:2) = tstr(1:1)
            tstr(1:1) = ' '
         endif
         do j=1,99
             if (tocapf(tstr).eq.tocapf(elemnt(j))) ianz(i) = j
         end do
         if (ianz(i).eq.0) goto 100

         do j=1,min(mxcon,8)
             if (icon(j).ne.0) then
                hascon = .true.
                iconn(1,i) = iconn(1,i) + 1
                iconn(iconn(1,i)+1,i) = icon(j)
             endif
         end do

      end do

c if read in charges are zero and old atoms array matches up with
c with new retain old charges

      ihasq = 1
      iq = 1

      do i=1,iatoms
         if (qat(nstor+i).ne.0.0d0) iq = 0
      end do

      if (.not.(iq.eq.1.and.ochg(idum,ianz))) then

         do i=1,iatoms
            qat(i) = qat(nstor+i)
         end do

      endif

      if (iq.eq.1.and..not.ochg(idum,ianz)) ihasq = 0

      do i=1,iatoms

         do j=1,3
            coo(j,nstor+i) = coo(j,i)
         end do

         call fr2crt(coo(1,i),xa,ya,yb,za,zb,zc)

         do j=1,3
             coo(j,i) = coo(j,i)/toang
         end do

         ianz(nstor+i) = ianz(i)
         iatclr(nstor+i) = 1

      end do

      if (.not.hascon) call doconn
      call dohcon(0)

      do i=1,iatoms
         do j=1,iconn(1,i)+1
            iconn(j,nstor+i) = iconn(j,i)
         end do
      end do

      nat = iatoms
      ncon = 1

      istat = 2
      ichx = 1
      call fdat(iop,1,istdbd,iuseab,moddma,idebug)
      icrtp = 1
      return

100   istat = 0
      if (icssr.eq.1) return

      iatoms = 0
      call rewfil
      call rfbio(idebug,1,istat)
      if (istat.eq.0) then
         call rdshlx(idebug,istat)
         if (istat.eq.0) then
            call rewfil
            call rdvasp(nptdum,nptdum,nptdum,0,istat,
     &           nstr,0,idebug)
            if (istat.eq.0) then
               call rewfil
               call rfdat(idebug,istat,refcod)
               if (istat.eq.0) then
                  call rewfil
                  call rdmsi(idebug,istat)
                  if (istat.eq.0) then
                     call rewfil
                     call rdcif(idebug,istat)
                     if (istat.eq.0) then
                        call rewfil
                        call rdconquest(nptdum,nptdum,nptdum,0,
     &                       istat,nstr,0,idebug)
                        
                        if (istat.ne.0) icrtp = 8
                        
                     else
                        icrtp = 7
                     end if
                  else
                     icrtp = 6
                  endif
               else
                  icrtp = 2
               endif
            else
               icrtp = 5
            endif
         else
            icrtp = 4
         endif
      else
         icrtp = 3
      endif
      if (istat.gt.0) iftyp = 6
      if (istat.eq.2) then
          ichx = 1
          call fdat(iop,1,istdbd,iuseab,moddma,idebug)
      endif

      return
      end

      subroutine rfbid(idebug,ifrst,istat,coo,qat,ianz,iatclr,iconn,
     &                 nat,norg,icent,inorm,ncon,nspg,kz,nopr,ir,it,
     &                 xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxel=100)
      parameter (mxsg=238)
      common /athlp/ iatoms, mxnat
      character*137 line
      character*2 elemnt,tolowf,atom,atomt
      common /elem/ elemnt(mxel)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      integer*2 ir,it
      character*7 spnam
      common /spgrnm/ spnam(mxsg)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /curlin/ line
      logical pbcon
      integer getlin
      character*137 string
      character*7 spgrt,spgru
      dimension coo(3,*),qat(*),ianz(*),iatclr(*),iconn(mxcon+1,*)
      dimension ir(3,3,192),it(3,192)

      pbcon = .false.
      toang = 0.52917706d0
      istat = 1

      if (getlin(0).eq.1) then
         if (index(line,'BIOSYM').ne.0) then
            call redel(line,3)
         else
            if (ifrst.eq.1) goto 100
            call redel(line,1)
         endif
      else
         goto 100
      endif

      iatoms = 0
      do while (getlin(0).eq.1)
          ktype = nxtwrd(string,nstr,itype,rtype)
          if (ktype.eq.1) then
             if (string(1:nstr).eq.'PBC') then
                 pbcon = .true.
                 ktype = nxtwrd(string,nstr,itype,rtype)
                 if (ktype.eq.3) a = rtype
                 ktype = nxtwrd(string,nstr,itype,rtype)
                 if (ktype.eq.3) b = rtype
                 ktype = nxtwrd(string,nstr,itype,rtype)
                 if (ktype.eq.3) c = rtype
                 ktype = nxtwrd(string,nstr,itype,rtype)
                 if (ktype.eq.3) alpha = rtype
                 ktype = nxtwrd(string,nstr,itype,rtype)
                 if (ktype.eq.3) beta = rtype
                 ktype = nxtwrd(string,nstr,itype,rtype)
                 if (ktype.eq.3) gamma = rtype
                 ktype = nxtwrd(string,nstr,itype,rtype)
                 nspg = 1
                 if (ktype.eq.1) then
                    if (nstr.ge.3) then
                       string = string(2:nstr-1)
                       call tocap(string,nstr-2)
                       spgrt = '       '
                       spgrt(1:nstr-2) = string(1:nstr-2) 
                       ifound = 0
                       do i=1,mxsg
                          spgru = spnam(i)
                          call tocap(spgru,7)
                          if (spgrt.eq.spgru) then
                             nspg = i
                             ifound = 1
                          endif
                       end do
                       if (ifound.eq.0)
     &                    call inferr('No spacegroup match',0)
                    endif
                 endif
             elseif (string(1:nstr).eq.'end'.or.
     &               string(1:nstr).eq.'END') then
                 goto 10
             else
c take care of silly space in label when count gt 99
                 ic = ichar(string(nstr:nstr))
                 if (.not.(ic.ge.48.and.ic.le.57)) then
                    ktype = nxtwrd(string,nstr,itype,rtype)
                 endif

                 iatoms = iatoms + 1
                 do i=1,3
                    ktype = nxtwrd(string,nstr,itype,rtype)
                    if (ktype.eq.3) coo(i,iatoms) = rtype/toang
                 end do
                 do i=1,3
                    ktype = nxtwrd(string,nstr,itype,rtype)
                 end do
                 ktype = nxtwrd(string,nstr,itype,rtype)
                 if (ktype.eq.1.and.(nstr.eq.1.or.nstr.eq.2)) then
                    if (nstr.eq.1) then
                       atomt(1:1) = ' '
                       atomt(2:2) = string(1:1)
                    else
                       atomt = string(1:2)
                    endif
                    atom = tolowf(atomt)
                    do j=1,mxel
                       if (atom .eq. tolowf(elemnt(j))) ianz(iatoms) = j
                    end do
                 endif
                 ktype = nxtwrd(string,nstr,itype,rtype)
                 if (ktype.eq.3) qat(iatoms) = rtype
                 iatclr(iatoms) = 1
             endif
          endif
      end do

10    continue
      call redel(line,1)

      ihasq = 1
      call doconn
      ncon = 1

      if (pbcon) then
         istat = 2
         nopr = 0
         call prcell(nspg,a,b,c,alpha,beta,gamma)
         call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)
         call cpmol(nat,norg,a,b,c,alpha,beta,gamma,
     &              coo,ianz,iatclr,iconn)
         call cprot(nspg,nopr,icent,ir,it,.false.)
         if (idebug.eq.1) call prop(nopr,ir,it)

         nat = iatoms

      else

         call doscal

      endif

      return

100   istat = 0
c      print*,'rfbio: error at line:',line
      return
      end

      subroutine geobio(formax,forrms,dismax,disrms,epoints,isav)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      common /curlin/ line
      character line*137
      character string*137
      integer getlin
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)

      igcvav = 1
      ifmxav = 0
      ifrmav = 0
      idmxav = 0
      idrmav = 0
      ieav = 1
      nepnts = 0
      ngeoms = 0

      call rewfil

      ngeoms = 0
      ept = 0.0d0

      do while (getlin(0).eq.1)

         if (index(line,'BIOSYM').ne.0) then
            call redel(line,1)
            if (getlin(0).eq.1) idum = 1
         endif

         do i=1,6
            ktype = nxtwrd(string,nstr,itype,rtype)
            if (ktype.eq.3) then
               ept = rtype
               goto 10
            endif
         end do
10       if (getlin(0).eq.0) goto 30

         do while (getlin(0).eq.1)
            ktype = nxtwrd(string,nstr,itype,rtype)
            if (ktype.eq.1) then
               if (string(1:nstr).eq.'end'.or.
     &             string(1:nstr).eq.'END') then
                      ngeoms = ngeoms + 1
                      epoints(ngeoms) = ept
                      ept = 0.0d0
                   goto 20
               endif
            endif
         end do

20       if (getlin(0).eq.0) goto 30

      end do

30    if (ngeoms.eq.0) then
          igcvav = 0
          ieav = 0
      else
          nepnts = ngeoms
      endif

      return
      end

      subroutine geofdt(formax,forrms,dismax,disrms,epoints,isav)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      character line*137
      character*8 refcode
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)

      igcvav = 1
      ifmxav = 0
      ifrmav = 0
      idmxav = 0
      idrmav = 0
      ieav = 0
      nepnts = 0
      ngeoms = 0

      call rewfil

      ngeoms = 0
      do while (.true.)

         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 100
         read(line,'(1x,a8,14x,i3)',end=100,err=100) refcode,ncards
         if (ncards.gt.1) then
             call redel(line,ncards-1)
         else
             goto 100
         endif

         ngeoms = ngeoms + 1

      end do

100   if (ngeoms.eq.0) then
          igcvav = 0
          ieav = 0
      else
          nepnts = ngeoms
      endif

      return
      end

      subroutine getxdd(ipnt,nat,istat,coo)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      integer getlin
      logical gnreal
      common /vasp/nheadx,natx

      character*137 line
      common /athlp/ iatoms, mxnat
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension tmp(3),coo(3,*)

      call rewfil

      if (nheadx.eq.0) nheadx=6

      call redel(line,nheadx+(ipnt-1)*(nat+1))

      nstor = mxnat - nat

      do i=1,nat
         if (getlin(0).eq.1) then
            if (gnreal(tmp,3,.false.)) then
               do j=1,3
                  coo(j,nstor+i) = tmp(j)
               end do
            else
               print*,'error reading XDATCAR'
            endif
         else
            print*,'error reading XDATCAR'
         endif
      end do

      return
      end

      subroutine geoxdt(formax,forrms,dismax,disrms,epoints,isav)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      logical opfil
      integer getlin
      character*137 str
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)
      common /vasp/nheadx,natx

      igcvav = 1
      ifmxav = 0
      ifrmav = 0
      idmxav = 0
      idrmav = 0
      ieav   = 0
      nepnts = 0
      ngeoms = 0

      nheadx = 6
      
      iuntmp = iun2
      iun2 = 48

      if (opfil(iun2,'XDATCAR',7,1,1,1)) then

         call confrm(1,istat)
         if (istat.eq.1) then
            if (getlin(1).eq.1) then
               ktype = nxtwrd(str,nstr,itype,rtype)

               if (ktype.eq.1) then

c found a string => new format XDATCAR
c find the number of lines ... 

                  nheadx = 8
                  nlines = 8
                  natx   = 0

                  do iline=1,6
                     itmp = getlin(1)
                     call gtplin(str)
                  end do
                  
                  do while (nxtwrd(str,nstr,itype,rtype).eq.2)
                     natx = natx + itype
                  end do

                  do while (getlin(1).eq.1)  
                     nlines = nlines + 1
                  end do

                  ngeoms = (nlines-nheadx)/(natx+1)

               elseif (ktype.eq.2) then

c we do not believe VASP for the number of geometries :p

                  natx = itype
                  nlines = 2

                  do while (getlin(1).eq.1)  
                     nlines = nlines + 1
                  end do

                  ngeoms = (nlines-nheadx)/(natx+1)

               end if
            endif
         else
            close(iun2) 
            iun2 = iuntmp
         endif

      endif

100   if (ngeoms.eq.0) then
          igcvav = 0
          ieav = 0
      else
          nepnts = ngeoms
      endif

      return
      end

      subroutine gtplin(str)
c     Put parsed line into str
      character*(*) str
      character*137 line
      common /curlin/ line

      str = line

      return
      end

      subroutine geochx(formax,forrms,dismax,disrms,epoints,isav)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      character line*137
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)

      igcvav = 1
      ifmxav = 0
      ifrmav = 0
      idmxav = 0
      idrmav = 0
      ieav = 0
      nepnts = 0
      ngeoms = 0

      call rewfil

      ngeoms = 0
      do while (.true.)

         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 100
         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 100
         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 100
         read(line,'(i4)',end=100,err=100) ncards
         ngeoms = ngeoms + 1
         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 100
         id = index(line,'Energy:')
         if (id.ne.0) then
            epoints(ngeoms) = reada(line,id+7,len(line))
            ieav = 1
         endif
         if (ncards.gt.1) then
             call redel(line,ncards)
         else
             goto 100
         endif


      end do

100   if (ngeoms.eq.0) then
          igcvav = 0
          ieav = 0
      else
          nepnts = ngeoms
      endif

      return
      end

      subroutine rfdad(idebug,istat,refcod,
     &                 coo,ianz,iatclr,iconn,
     &                 nat,norg,icent,inorm,ncon,nspg,kz,nopr,ir,it,
     &                 xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxel=100)
      parameter (mxtyp=8)
      parameter (mxop=2822)
      parameter (mxrop=60)
      parameter (mxtop=56)
      parameter (mxspg=232)
      parameter (mxspg6=158*6)
      parameter (mxsg=238)
      common /athlp/ iatoms, mxnat
      character*137 line
      character*1 tstr
      character*2 atom,elemnt,tocapf
      character*5 atmtmp(3)
      character*8 spg,refcod
      common /elem/   elemnt(mxel)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /symgrp/ irr(3,3,mxrop),itt(3,mxtop),nop(mxspg),iri(mxop),
     &                iti(mxop),ipt(mxspg),icen(mxspg),ispm(mxsg)
      common /settng/ jset(3,9),iset(6,mxspg),nset
      character*7 spnam
      common /spgrnm/ spnam(mxsg)
      integer*2 ir,it
      common /cllhlp/ ico(3,8),icn(4,8),mcol(32)
      common /exmass/ emass(mxel),cutt(4,4)
      character*2 labtyp
      common /buck/ alk(mxtyp,mxtyp),blk(mxtyp,mxtyp),clk(mxtyp,mxtyp),
     &              alkp(mxtyp),blkp(mxtyp),clkp(mxtyp),
     &              iptyp(mxtyp),labtyp(mxtyp)
      common /buick/ aij(mxel),bij(mxel),cij(mxel),aijhp,bijhp,cijhp
      character*9 fmt
      integer*2 irr,itt
      dimension ii(3,3),tr1(3)
      dimension icelld(6),iprec(6),icesd(6),icon(40)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*)
      dimension ir(3,3,192),it(3,192)
      
c mcol
c 1       white
c 2       red
c 3       green
c 4       blueish grey
c 5       blue
c 6       purple
c 7       ecru (cell axes)
c 8       light orange

      data mcol /15,1,3,7,4,9,13,6,2,5,8,10,11,12,14,14,16*15/
      data ico(1,1),ico(2,1),ico(3,1) /0,0,0/
      data ico(1,2),ico(2,2),ico(3,2) /1,0,0/
      data ico(1,3),ico(2,3),ico(3,3) /0,1,0/
      data ico(1,4),ico(2,4),ico(3,4) /0,0,1/
      data ico(1,5),ico(2,5),ico(3,5) /1,1,0/
      data ico(1,6),ico(2,6),ico(3,6) /0,1,1/
      data ico(1,7),ico(2,7),ico(3,7) /1,0,1/
      data ico(1,8),ico(2,8),ico(3,8) /1,1,1/
      data icn(1,1),icn(2,1),icn(3,1),icn(4,1) /3,2,3,4/
      data icn(1,2),icn(2,2),icn(3,2),icn(4,2) /3,1,5,7/
      data icn(1,3),icn(2,3),icn(3,3),icn(4,3) /3,1,5,6/
      data icn(1,4),icn(2,4),icn(3,4),icn(4,4) /3,1,6,7/
      data icn(1,5),icn(2,5),icn(3,5),icn(4,5) /3,2,3,8/
      data icn(1,6),icn(2,6),icn(3,6),icn(4,6) /3,3,4,8/
      data icn(1,7),icn(2,7),icn(3,7),icn(4,7) /3,2,4,8/
      data icn(1,8),icn(2,8),icn(3,8),icn(4,8) /3,5,6,7/

      data ispm /1,2,19,3,12,20,21,22,11,23,14,15,17,4,5,225,226,
     &           16,6,24,25,26,27,28,227,29,30,31,13,32,33,34,9,35,
     &           36,37,38,39,40,41,42,43,44,45,46,47,228,48,49,50,51,
     &           52,53,54,55,56,57,58,59,10,7,8,60,61,62,63,64,65,66,
     &           67,68,69,70,71,72,18,73,74,75,76,77,78,79,80,81,82,
     &           83,84,229,85,86,87,88,89,90,91,92,93,230,94,95,96,
     &           97,98,231,99,100,101,102,
     &           103,104,105,106,107,108,109,110,111,112,113,114,115,
     &           116,117,118,119,120,121,122,123,124,125,126,127,128,
     &           129,130,131,132,133,134,135,136,137,138,218,140,219,
     &           142,143,232,144,145,146,147,148,149,150,151,221,222,
     &           154,155,156,157,223,224,160,161,162,163,164,165,166,
     &           167,168,0,169,170,171,172,173,0,174,175,176,0,177,
     &           178,179,180,181,182,183,184,185,186,187,188,0,189,
     &           190,191,192,193,194,195,0,196,197,198,199,200,0,
     &           201,202,203,204,205,206,207,208,209,210,211,212,213,
     &           214,215,216,217,139,141,220,152,153,158,159/


      data (itt(i,1),i=1,3)           /0,0,0/
      data (itt(i,2),i=1,3)           /6,0,0/
      data (itt(i,3),i=1,3)           /0,6,0/
      data (itt(i,4),i=1,3)           /0,0,6/
      data (itt(i,5),i=1,3)           /6,6,0/
      data (itt(i,6),i=1,3)           /6,0,6/
      data (itt(i,7),i=1,3)           /0,6,6/
      data (itt(i,8),i=1,3)           /6,6,6/
      data (itt(i,9),i=1,3)           /0,0,3/
      data (itt(i,10),i=1,3)          /0,0,9/
      data (itt(i,11),i=1,3)          /3,3,3/
      data (itt(i,12),i=1,3)          /3,9,9/
      data (itt(i,13),i=1,3)          /9,3,9/
      data (itt(i,14),i=1,3)          /9,9,3/
      data (itt(i,15),i=1,3)          /0,3,3/
      data (itt(i,16),i=1,3)          /0,9,9/
      data (itt(i,17),i=1,3)          /6,3,9/
      data (itt(i,18),i=1,3)          /6,9,3/
      data (itt(i,19),i=1,3)          /3,0,3/
      data (itt(i,20),i=1,3)          /3,6,9/
      data (itt(i,21),i=1,3)          /9,0,9/
      data (itt(i,22),i=1,3)          /9,6,3/
      data (itt(i,23),i=1,3)          /3,9,6/
      data (itt(i,24),i=1,3)          /3,3,0/
      data (itt(i,25),i=1,3)          /9,9,0/
      data (itt(i,26),i=1,3)          /9,3,6/
      data (itt(i,27),i=1,3)          /0,6,3/
      data (itt(i,28),i=1,3)          /6,0,9/
      data (itt(i,29),i=1,3)          /9,3,3/
      data (itt(i,30),i=1,3)          /9,9,9/
      data (itt(i,31),i=1,3)          /6,6,3/
      data (itt(i,32),i=1,3)          /6,6,9/
      data (itt(i,33),i=1,3)          /0,6,9/
      data (itt(i,34),i=1,3)          /6,0,3/
      data (itt(i,35),i=1,3)          /3,9,3/
      data (itt(i,36),i=1,3)          /3,3,9/
      data (itt(i,37),i=1,3)          /0,0,4/
      data (itt(i,38),i=1,3)          /0,0,8/
      data (itt(i,39),i=1,3)          /4,8,8/
      data (itt(i,40),i=1,3)          /8,4,4/
      data (itt(i,41),i=1,3)          /4,8,2/
      data (itt(i,42),i=1,3)          /8,4,10/
      data (itt(i,43),i=1,3)          /0,0,10/
      data (itt(i,44),i=1,3)          /0,0,2/
      data (itt(i,45),i=1,3)          /0,3,9/
      data (itt(i,46),i=1,3)          /9,0,3/
      data (itt(i,47),i=1,3)          /3,9,0/
      data (itt(i,48),i=1,3)          /0,9,3/
      data (itt(i,49),i=1,3)          /6,3,3/
      data (itt(i,50),i=1,3)          /6,9,9/
      data (itt(i,51),i=1,3)          /9,6,9/
      data (itt(i,52),i=1,3)          /3,0,9/
      data (itt(i,53),i=1,3)          /3,6,3/
      data (itt(i,54),i=1,3)          /3,3,6/
      data (itt(i,55),i=1,3)          /9,9,6/
      data (itt(i,56),i=1,3)          /9,3,0/

      data ((irr(i,j,1),j=1,3),i=1,3) /2,1,1,1,2,1,1,1,2/
      data ((irr(i,j,2),j=1,3),i=1,3) /0,1,1,1,2,1,1,1,2/
      data ((irr(i,j,3),j=1,3),i=1,3) /2,1,1,1,0,1,1,1,2/
      data ((irr(i,j,4),j=1,3),i=1,3) /2,1,1,1,0,1,1,1,0/
      data ((irr(i,j,5),j=1,3),i=1,3) /0,1,1,1,2,1,1,1,0/
      data ((irr(i,j,6),j=1,3),i=1,3) /0,1,1,1,0,1,1,1,2/
      data ((irr(i,j,7),j=1,3),i=1,3) /1,0,1,2,1,1,1,1,2/
      data ((irr(i,j,8),j=1,3),i=1,3) /1,2,1,0,1,1,1,1,2/
      data ((irr(i,j,9),j=1,3),i=1,3) /0,1,1,1,0,1,1,1,0/
      data ((irr(i,j,10),j=1,3),i=1,3) /2,1,1,1,2,1,1,1,0/
      data ((irr(i,j,11),j=1,3),i=1,3) /1,0,1,0,1,1,1,1,2/
      data ((irr(i,j,12),j=1,3),i=1,3) /1,0,1,2,1,1,1,1,0/
      data ((irr(i,j,13),j=1,3),i=1,3) /1,2,1,0,1,1,1,1,0/
      data ((irr(i,j,14),j=1,3),i=1,3) /1,2,1,2,1,1,1,1,0/
      data ((irr(i,j,15),j=1,3),i=1,3) /1,0,1,0,1,1,1,1,0/
      data ((irr(i,j,16),j=1,3),i=1,3) /1,2,1,2,1,1,1,1,2/
      data ((irr(i,j,17),j=1,3),i=1,3) /1,0,1,2,0,1,1,1,2/
      data ((irr(i,j,18),j=1,3),i=1,3) /0,2,1,0,1,1,1,1,2/
      data ((irr(i,j,19),j=1,3),i=1,3) /0,2,1,1,2,1,1,1,0/
      data ((irr(i,j,20),j=1,3),i=1,3) /2,1,1,2,0,1,1,1,0/
      data ((irr(i,j,21),j=1,3),i=1,3) /2,0,1,1,0,1,1,1,0/
      data ((irr(i,j,22),j=1,3),i=1,3) /0,1,1,0,2,1,1,1,0/
      data ((irr(i,j,23),j=1,3),i=1,3) /1,1,2,2,1,1,1,2,1/
      data ((irr(i,j,24),j=1,3),i=1,3) /1,2,1,1,1,2,2,1,1/
      data ((irr(i,j,25),j=1,3),i=1,3) /0,1,1,1,1,0,1,0,1/
      data ((irr(i,j,26),j=1,3),i=1,3) /1,1,0,1,0,1,0,1,1/
      data ((irr(i,j,27),j=1,3),i=1,3) /0,2,1,1,2,1,1,1,2/
      data ((irr(i,j,28),j=1,3),i=1,3) /2,1,1,2,0,1,1,1,2/
      data ((irr(i,j,29),j=1,3),i=1,3) /0,1,1,0,2,1,1,1,2/
      data ((irr(i,j,30),j=1,3),i=1,3) /2,0,1,1,0,1,1,1,2/
      data ((irr(i,j,31),j=1,3),i=1,3) /1,2,1,0,2,1,1,1,2/
      data ((irr(i,j,32),j=1,3),i=1,3) /2,0,1,2,1,1,1,1,2/
      data ((irr(i,j,33),j=1,3),i=1,3) /1,0,1,2,0,1,1,1,0/
      data ((irr(i,j,34),j=1,3),i=1,3) /0,2,1,0,1,1,1,1,0/
      data ((irr(i,j,35),j=1,3),i=1,3) /1,1,2,0,1,1,1,0,1/
      data ((irr(i,j,36),j=1,3),i=1,3) /1,2,1,1,1,0,0,1,1/
      data ((irr(i,j,37),j=1,3),i=1,3) /1,1,0,2,1,1,1,0,1/
      data ((irr(i,j,38),j=1,3),i=1,3) /1,0,1,1,1,2,0,1,1/
      data ((irr(i,j,39),j=1,3),i=1,3) /1,1,0,0,1,1,1,2,1/
      data ((irr(i,j,40),j=1,3),i=1,3) /1,0,1,1,1,0,2,1,1/
      data ((irr(i,j,41),j=1,3),i=1,3) /1,1,2,2,1,1,1,0,1/
      data ((irr(i,j,42),j=1,3),i=1,3) /1,2,1,1,1,2,0,1,1/
      data ((irr(i,j,43),j=1,3),i=1,3) /1,1,0,2,1,1,1,2,1/
      data ((irr(i,j,44),j=1,3),i=1,3) /1,0,1,1,1,2,2,1,1/
      data ((irr(i,j,45),j=1,3),i=1,3) /1,1,2,0,1,1,1,2,1/
      data ((irr(i,j,46),j=1,3),i=1,3) /1,2,1,1,1,0,2,1,1/
      data ((irr(i,j,47),j=1,3),i=1,3) /0,1,1,1,1,2,1,2,1/
      data ((irr(i,j,48),j=1,3),i=1,3) /1,1,0,1,2,1,2,1,1/
      data ((irr(i,j,49),j=1,3),i=1,3) /2,1,1,1,1,2,1,0,1/
      data ((irr(i,j,50),j=1,3),i=1,3) /1,1,2,1,2,1,0,1,1/
      data ((irr(i,j,51),j=1,3),i=1,3) /2,1,1,1,1,0,1,2,1/
      data ((irr(i,j,52),j=1,3),i=1,3) /1,1,2,1,0,1,2,1,1/
      data ((irr(i,j,53),j=1,3),i=1,3) /0,1,1,1,1,2,1,0,1/
      data ((irr(i,j,54),j=1,3),i=1,3) /1,1,0,1,2,1,0,1,1/
      data ((irr(i,j,55),j=1,3),i=1,3) /0,1,1,1,1,0,1,2,1/
      data ((irr(i,j,56),j=1,3),i=1,3) /1,1,0,1,0,1,2,1,1/
      data ((irr(i,j,57),j=1,3),i=1,3) /2,1,1,1,1,2,1,2,1/
      data ((irr(i,j,58),j=1,3),i=1,3) /1,1,2,1,2,1,2,1,1/
      data ((irr(i,j,59),j=1,3),i=1,3) /2,1,1,1,1,0,1,0,1/
      data ((irr(i,j,60),j=1,3),i=1,3) /1,1,2,1,0,1,0,1,1/

c               P21 P21/c C2/c P212121 pbca pnma  pna21  pbcn   Cc

      data (iti(i),i=1,397) /
     &          3, 7, 5,4,8, 5,7,6, 5,7,6, 8,3,6, 8,5,4, 5,4,8, 5,4,8, 
     &          5,1,5, 6,2,4, 3, 5,1,5, 5,5,1, 4, 9,4,10, 1, 3, 4,
     &          5,1,5,  1,  5,1,5,4,8,4,8, 5,1,5,1,5,1,5, 
     &          7,6,5,1,7,6,5,1,7,6,5,1,7,6,5, 8,1,8,1,8,1,8,
     &          6,5,7,8,3,4,2, 1,4,4, 4,4,1, 2,2,1, 7,7,1,
     &          1,6,6, 5,5,1, 8,8,1, 5,1,5,1,5,1,5, 5,1,5,4,8,4,8, 
     &          5,4,8,4,8,1,5, 7,1,7,1,7,1,7, 7,3,4,3,4,1,7,
     &          7,2,8,2,8,1,7, 7,5,6,5,6,1,7, 
     &          7,6,5,1,7,6,5,1,7,6,5,1,7,6,5, 
     &          7,6,5,11,12,13,14,11,12,13,14,6,5,1,7, 8,1,8,1,8,1,8,
     &          8,5,4,5,4,1,8, 8,2,7,2,7,1,8, 1,8,8,1,1,8,8,
     &          4,4,1, 3,2,5, 2,1,2, 7,8,2, 1,6,6, 6,4,2, 5,5,1,
     &          6,7,5, 3,7,4, 8,8,1, 2,3,5, 5,1,5,4,8,4,8,
     &          5,1,5,7,6,7,6, 5,1,5,1,5,1,5, 5,4,8,4,8,1,5,
     &          5,1,5,3,2,3,2, 5,6,7,4,8,2,3,
     &          7,6,5,1,7,6,5,1,7,6,5,1,7,6,5,
     &          7,6,5,15,16,17,18,19,20,21,22,23,24,25,26,8,1,8,1,8,1,8,
     &          8,5,4,5,4,1,8, 8,5,4,7,2,6,3, 8,1,8,3,6,3,6, 1,1,1,
     &          4,1,4, 10,4,9, 8,1,8,1,8,1,8, 8,27,28,1,8,27,28, 1,1,1,
     &          8,1,8,1,8,1,8, 1,1,1, 4,1,4, 5,5,1,5,1,5,1/

      data (iti(i),i=398,753) /
     &          7,5,6, 8,1,8,1,8,1,8, 8,29,12,6,3,30,11, 5,1,5,1,5,1,5,
     &          9,4,10,1,10,4,9, 31,4,32,1,32,4,31, 4,1,4,1,4,1,4,
     &          8,1,8,1,8,1,8, 10,4,9,1,9,4,10, 32,4,31,1,31,4,32,
     &          8,1,1,8,8,1,1,8,8,1,1,8,8,1,8,
     &          27,27,8,8,28,28,1,8,28,28,1,1,27,27,8, 1,1,1,5,5,5,5,
     &          4,1,4,4,1,4,1, 8,1,8,8,1,8,1, 1,1,1,4,4,4,4,
     &          1,1,1,8,8,8,8, 4,1,4,5,8,5,8, 
     &          8,1,8,1,8,1,8,1,8,1,8,1,8,1,8,
     &          8,1,8,1,8,1,8,4,5,4,5,4,5,4,5,
     &          8,27,28,1,8,27,28,1,8,27,28,1,8,27,28,
     &          8,27,28,1,8,27,28,4,5,33,34,4,5,33,34,
     &          1,1,1,1,1,1,1, 1,1,1,4,4,4,4, 1,1,1,5,5,5,5,
     &          1,1,1,8,8,8,8, 1,1,1,1,1,1,1, 1,1,1,4,4,4,4,
     &          1,1,1,5,5,5,5, 1,1,1,8,8,8,8,
     &          8,1,8,1,8,1,8,1,8,1,8,1,8,1,8,
     &          8,1,8,1,8,1,8,4,5,4,5,4,5,4,5,
     &          8,1,8,1,8,1,8,1,8,1,8,1,8,1,8,
     &          8,1,8,1,8,1,8,27,28,27,28,27,28,27,28,
     &          1,1,1,1,1,1,1, 1,1,1,4,4,4,4, 2,5,3,3,5,2,1,
     &          2,5,3,7,8,6,4, 1,1,1,5,5,5,5, 1,1,1,8,8,8,8/

      data (iti(i),i=754,1083) /
     &          1,5,5,1,5,5,1,5,5,5,5,1,1,1,1,
     &          2,5,3,6,4,7,8, 4,1,4,1,4,1,4, 4,1,4,4,1,4,1,
     &          6,5,7,3,8,2,4, 6,5,7,7,5,6,1, 4,1,4,5,8,5,8,
     &          8,1,8,8,1,8,1, 6,5,7,2,4,3,8, 6,5,7,6,1,7,5,
     &          8,1,8,1,8,1,8,1,8,1,8,1,8,1,8,
     &          1,4,1,4,1,4,4,8,8,5,8,5,8,5,5,
     &          8,13,35,6,3,14,36,1,8,36,14,6,3,35,13,
     &          8,13,35,6,3,14,36,4,5,11,30,2,7,12,29,
     &          1,1, 37,38, 38,37, 39,40,1,1,39,39,40,40, 1,1,
     &          39,40,1,1,39,39,40,40, 1,1,1,1,1, 1,1,1,1,1,
     &          37,38,1,38,37, 38,37,37,38,1 ,38,37,1,37,38, 
     &          1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 1,4,1,4,4,
     &          1,1,4,4,4, 
     &          39,40,1,1,39,39,40,40,1,39,40,1,1,39,39,40,40,
     &          39,40,1,1,39,39,40,40,4,41,42,4,4,41,41,42,42,
     &          1,1,1,1,1, 1,1,4,4,4, 1,1,1,1,1, 1,1,4,4,4,
     &          39,40,1,1,39,39,40,40,1,39,40,1,1,39,39,40,40,
     &          39,40,1,1,39,39,40,40,4,41,42,4,4,41,41,42,42,
     &          1,1,1,1,1, 37,38,4,43,44, 38,37,4,44,43,
     &          38,37,1,38,37, 37,38,1,37,38, 1,1,4,4,4/

      data (iti(i),i=1084,1483) /
     &          1,1,1,1,1, 1,1,1,1,1, 1,1,4,4,4,
     &          37,38,37,1,38,4,43,44,43,4,44,
     &          43,44,38,1,4,43,37,38,44,4,37,
     &          38,37,38,1,37,1,38,37,38,1,37,
     &          37,38,37,1,38,1,37,38,37,1,38,
     &          1,1,4,4,4,1,1,1,4,4,4, 1,4,1,4,1,4,1,4,1,4,4,
     &          1,1,4,4,4,4,4,4,1,1,1, 1,1,4,4,4,1,1,1,4,4,4,
     &          1,1,4,4,4,4,4,4,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,
     &          1,1,4,4,4,1,1,1,4,4,4, 1,1,1,1,1,1,1,1,1,1,1,
     &          1,4,1,4,1,4,1,4,1,4,4, 4,1,1,4,4,1,1,4,4,1,4,
     &          1,1,4,4,4,1,1,1,4,4,4, 1,1,1,1,1,1,1,1,1,1,1,
     &          1,1,1,1,1,1,1,1,1,1,1,7,7,7,7,7,7,7,7,7,7,7,7,
     &          6,6,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,
     &          8,1,1,8,8,1,8,1,1,8,8,1,8,1,1,8,8,1,8,1,1,8,8,
     &          1,1,5,5,5,7,7,7,6,6,6,
     &          8,1,1,8,8,5,4,5,5,4,4,7,2,7,7,2,2,6,3,6,6,3,3,
     &          1,1,5,5,5,7,7,7,6,6,6,
     &          1,1,1,1,1,1,1,1,1,1,1,7,7,7,7,7,7,7,7,7,7,7,7,
     &          6,6,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,
     &          7,6,5,1,1,6,5,5,7,7,6,15,16,17,18,15,15,17,18,18,
     &          16,16,17,19,20,21,22,19,19,21,22,22,20,20,21,23,24,
     &          25,26,23,23,25,26,26,24,24,25/

      data (iti(i),i=1484,1869) /
     &          8,1,1,8,8,1,8,1,1,8,8,1,8,1,1,8,8,1,8,1,1,8,8,
     &          1,1,6,6,6,5,5,5,7,7,7,
     &          8,1,1,8,8,6,3,6,6,3,3,5,4,5,5,4,4,7,2,7,7,2,2,
     &          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     &          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     &          7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
     &          6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
     &          5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
     &          7,6,5,1,1,6,5,5,7,7,6,11,13,12,14,11,11,12,14,14,
     &          13,13,12,5,6,7,1,5,5,7,1,1,6,6,7,13,11,14,12,13,
     &          13,14,12,12,11,11,14,7,1,5,6,7,7,5,6,6,1,1,5,14,
     &          12,13,11,14,14,13,11,11,12,12,13,6,5,1,7,6,6,1,7,
     &          7,5,5,1,12,14,11,13,12,12,11,13,13,14,14,11,
     &          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     &          8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
     &          1,1,13,13,13,5,5,5,11,11,11,7,7,7,12,12,12,6,6,6,
     &          14,14,14,
     &          1,1,35,35,35,5,5,5,30,30,30,7,7,7,29,29,29,6,6,6,
     &          36,36,36,
     &          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/

      data (iti(i),i=1870,2363) /
     &            7,6,5,1,1,6,5,5,7,7,6,1,6,7,5,1,1,7,5,5,6,6,7,
     &          1,7,6,5,1,1,6,5,5,7,7,6,1,6,7,5,1,1,7,5,5,6,6,7,
     &          1,7,6,5,1,1,6,5,5,7,7,6,1,6,7,5,1,1,7,5,5,6,6,7,
     &          1,7,6,5,1,1,6,5,5,7,7,6,1,6,7,5,1,1,7,5,5,6,6,7,
     &            8,1,1,8,8,1,8,1,1,8,8,1,8,1,1,8,8,1,8,1,1,8,8,
     &          1,8,1,1,8,8,1,8,1,1,8,8,1,8,1,1,8,8,1,8,1,1,8,8,
     &          1,1,1,1,1,1,1,1,1,1,1,8,8,8,8,8,8,8,8,8,8,8,8,
     &            7,6,5,1,1,6,5,5,7,7,6,1,7,6,5,1,1,6,5,5,7,7,6,
     &          1,7,6,5,1,1,6,5,5,7,7,6,1,7,6,5,1,1,6,5,5,7,7,6,
     &          4,2,3,8,4,4,3,8,8,2,2,3,4,2,3,8,4,4,3,8,8,2,2,3,
     &          4,2,3,8,4,4,3,8,8,2,2,3,4,2,3,8,4,4,3,8,8,2,2,3,
     &            8,1,1,8,8,5,4,5,5,4,4,7,2,7,7,2,2,6,3,6,6,3,3,
     &          11,30,11,11,30,30,14,36,14,14,36,36,13,35,13,13,
     &          35,35,12,29,12,12,29,29,
     &          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     &          1,1,5,5,5,7,7,7,6,6,6,8,8,8,4,4,4,3,3,3,2,2,2,
     &          1,1,1,1,1,1,1,1,1,1,1,8,8,8,8,8,8,8,8,8,8,8,8,
     &          1,1,5,5,5,7,7,7,6,6,6,1,1,1,5,5,5,6,6,6,7,7,7,
     &            7,6,5,1,1,6,5,5,7,7,6,1,7,6,5,1,1,6,5,5,7,7,6,
     &          1,7,6,5,1,1,6,5,5,7,7,6,1,7,6,5,1,1,6,5,5,7,7,6,
     &          1,6,7,5,1,1,7,5,5,6,6,7,1,6,7,5,1,1,7,5,5,6,6,7,
     &          1,6,7,5,1,1,7,5,5,6,6,7,1,6,7,5,1,1,7,5,5,6,6,7/

      data (iti(i),i=2364,2742) /
     &            7,6,5,1,1,6,5,5,7,7,6,1,7,6,5,1,1,6,5,5,7,7,6,
     &          1,7,6,5,1,1,6,5,5,7,7,6,1,7,6,5,1,1,6,5,5,7,7,6,
     &          4,2,3,8,4,4,3,8,8,2,2,3,4,2,3,8,4,4,3,8,8,2,2,3,
     &          4,2,3,8,4,4,3,8,8,2,2,3,4,2,3,8,4,4,3,8,8,2,2,3,
     &            1,1,1,1,1,19,19,19,19,19,19,16,16,16,16,16,16,
     &          25,25,25,25,25,25,7,7,7,7,7,7,20,20,20,20,20,20,
     &          15,15,15,15,15,15,26,26,26,26,26,26,6,6,6,6,6,6,
     &          21,21,21,21,21,21,18,18,18,18,18,18,23,23,23,23,
     &          23,23,5,5,5,5,5,5,22,22,22,22,22,22,17,17,17,17,
     &          17,17,24,24,24,24,24,24,
     &            1,1,4,4,4,15,15,15,45,45,45,19,19,19,46,46,46,
     &          24,24,24,47,47,47,7,6,5,16,17,18,20,21,22,23,26,
     &          25,7,6,5,16,17,18,20,21,22,23,26,25,7,6,5,16,17,
     &          18,20,21,22,23,26,25,3,2,8,48,49,50,51,52,53,54,
     &          55,56,3,2,8,48,49,50,51,52,53,54,55,56,3,2,8,48,
     &          49,50,51,52,53,54,55,56,
     &            8,1,1,8,8,1,8,1,1,8,8,1,8,1,1,8,8,1,8,1,1,8,8,
     &          1,8,1,1,8,8,1,8,1,1,8,8,1,8,1,1,8,8,1,8,1,1,8,8,
     &            8,1,1,8,8,6,3,6,6,3,3,5,4,5,5,4,4,7,2,7,7,2,2,
     &          11,30,11,11,30,30,12,29,12,12,29,29,14,36,14,14,
     &          36,36,13,35,13,13,35,35/

      data (iti(i),i=2743,2822) /
     &          8, 1,1, 1,1, 
     &          39,40,1,1,39,39,40,40,1,40,39,1,1,40,40,39,39,
     &          1,1,1,1,1, 1,8,1,8,8, 1,1,1,1,1, 1,8,1,8,8,
     &          1,1,1, 1,4,4, 1,1,1, 1,1,1, 1,1,1,1,1,1,1,
     &          1,1,1,1,1,1,1, 1,4,4,1,1,4,4, 37,38,1,37,38/

c               C2     Pca21 P21/m C2/m P21212 P2/c P41  P2  Pm Pc
c               Cm     P2/m     C2221        C222
c               F222                              I222
c               I212121        Pmc21  Pcc2   Pma2   Pnc2
c               Pmn21  Pba2   Pnn2        Cmm2           Cmc21
c                   Ccc2          Amm2           Abm2
c                   Ama2            Aba2
c                            Fmm2
c                            Fdd2                          Imm2
c                    Iba2            Ima2         Pnnn
c               Pccm   Pban   Pmma   Pnna   Pmna   Pcca   Pbam
c               Pccn   Pbcm   Pnnm   Pmmn      Cmcm
c               Cmca             Cmmm          Cccm
c               Cmma             Ccca
c                       Fmmm
c                       Fddd                   Immm
c               Ibam             Ibca          Imma         P4
c               P42   P43        I4            I41          P-4
c               I-4              P4/m          P42/m        P4/n
c               P42/n       I4/m          I41/a            P4212
c                      P4122          P41212              P4222
c                      P42212         P4322               P43212
c               I422
c               I4122                    P4bm

      data (iri(i),i=1,397) /
     &          5, 5, 1,5,5, 4,5,6, 2,3,6, 2,3,6, 2,3,6, 2,3,6, 1,3,3,
     &          1,5,5, 2,3,6, 5, 1,5,5, 4,5,6, 5,  7,6,8, 5, 3, 3,
     &          1,3,3,  5,  1,4,4,5,5,6,6, 1,4,4,5,5,6,6,
     &          1,1,1,4,4,4,4,5,5,5,5,6,6,6,6, 1,4,4,5,5,6,6,
     &          6,4,5,1,6,4,5, 2,3,6, 2,3,6, 2,3,6, 2,3,6, 2,3,6,
     &          2,3,6, 2,3,6, 1,2,2,3,3,6,6, 1,2,2,3,3,6,6,
     &          1,2,2,3,3,6,6, 1,2,2,3,3,6,6, 1,2,2,3,3,6,6, 
     &          1,2,2,3,3,6,6, 1,2,2,3,3,6,6, 
     &          1,1,1,2,2,2,2,3,3,3,3,6,6,6,6,
     &          1,1,1,2,2,2,2,3,3,3,3,6,6,6,6, 1,2,2,3,3,6,6,
     &          1,2,2,3,3,6,6, 1,2,2,3,3,6,6, 6,9,10,4,5,2,3,
     &          2,3,6, 2,3,6, 2,3,6, 2,3,6, 2,3,6, 2,3,6, 2,3,6,
     &          2,3,6, 2,3,6, 2,3,6, 2,3,6, 1,2,2,3,3,6,6,
     &          1,2,2,3,3,6,6, 1,2,2,3,3,6,6, 1,2,2,3,3,6,6,
     &          1,2,2,3,3,6,6, 1,2,2,3,3,6,6, 
     &          1,1,1,2,2,2,2,3,3,3,3,6,6,6,6,
     &          1,1,1,2,2,2,2,3,3,3,3,6,6,6,6, 1,2,2,3,3,6,6,
     &          1,2,2,3,3,6,6, 1,2,2,3,3,6,6, 1,2,2,3,3,6,6, 7,6,8,
     &          7,6,8, 7,6,8, 1,7,7,6,6,8,8, 1,7,7,6,6,8,8, 12,6,13,
     &          1,12,12,6,6,13,13, 7,6,8, 7,6,8, 7,8,6,9,12,10,13/

      data (iri(i),i=398,725) /
     &          7,6,8, 1,7,7,6,6,8,8, 1,7,7,6,6,8,8, 7,6,8,14,4,15,5,
     &          7,6,8,5,14,4,15, 7,6,8,14,4,15,5, 7,6,8,5,14,4,15,
     &          7,6,8,14,4,15,5, 7,6,8,5,14,4,15, 7,6,8,14,4,15,5,
     &          1,7,4,7,4,6,14,6,14,8,5,8,5,15,15,
     &          7,4,6,14,8,5,15,1,7,4,6,14,8,5,15, 7,6,8,2,16,3,11,
     &          7,6,8,2,16,3,11, 7,6,8,2,16,3,11, 7,6,8,2,16,3,11,
     &          7,6,8,2,16,3,11, 7,6,8,2,16,3,11, 
     &          1,7,7,6,6,8,8,2,2,16,16,3,3,11,11,
     &          1,7,7,6,6,8,8,2,2,16,16,3,3,11,11,
     &          1,7,7,6,6,8,8,2,2,16,16,3,3,11,11,
     &          1,7,7,6,6,8,8,2,2,16,16,3,3,11,11,
     &          12,6,13,5,16,4,11, 12,6,13,5,16,4,11,
     &          12,6,13,5,16,4,11, 12,6,13,5,16,4,11,
     &          12,6,13,2,14,3,15, 12,6,13,2,14,3,15,
     &          12,6,13,2,14,3,15, 12,6,13,2,14,3,15,
     &          1,12,12,6,6,13,13,2,2,14,14,3,3,15,15,
     &          1,12,12,6,6,13,13,2,2,14,14,3,3,15,15,
     &          1,12,12,6,6,13,13,5,5,16,16,4,4,11,11,
     &          1,12,12,6,6,13,13,16,16,4,4,11,11,5,5,
     &          7,6,8,2,16,3,11, 7,6,8,2,16,3,11/

      data (iri(i),i=726,1068) /
     &          7,6,8,2,16,3,11, 7,6,8,2,16,3,11,
     &          7,6,8,2,16,3,11, 7,6,8,2,16,3,11,
     &          13,4,9,6,11,7,12,5,10,16,8,2,15,3,14,
     &          7,6,8,2,16,3,11, 7,6,8,2,16,3,11, 7,6,8,2,16,3,11,
     &          7,6,8,2,16,3,11, 7,6,8,2,16,3,11, 7,6,8,2,16,3,11,
     &          7,6,8,2,16,3,11, 7,6,8,2,16,3,11, 7,6,8,2,16,3,11,
     &          1,7,7,6,6,8,8,2,2,16,16,3,3,11,11,
     &          7,4,6,14,8,5,15,1,7,4,6,14,8,5,15,
     &          1,7,7,6,6,8,8,2,2,16,16,3,3,11,11,
     &          1,7,7,6,6,8,8,2,2,16,16,3,3,11,11,
     &          17,18, 17,18, 17,18, 1,1,17,18,17,18,17,18,
     &          17,18, 1,1,17,18,17,18,17,18, 17,18,15,19,20,
     &          17,18,14,21,22, 17,18,14,21,22, 17,18,15,19,20,
     &          17,18,14,21,22, 23,24,15,25,26, 17,18,11,27,28,
     &          17,16,18,29,30, 17,11,18,28,27, 17,18,16,30,29,
     &          1,1,17,18,17,18,17,18,11,11,11,27,28,27,28,27,28,
     &          1,1,17,18,17,18,17,18,11,11,11,27,28,27,28,27,28,
     &          17,18,16,30,29, 17,18,16,30,29, 17,18,11,27,28,
     &          17,18,11,27,28, 
     &          1,1,17,18,17,18,17,18,11,11,11,27,28,27,28,27,28,
     &          1,1,17,18,17,18,17,18,11,11,11,27,28,27,28,27,28,
     &          17,18,6,31,32, 17,18,6,31,32, 17,18,6,31,32/

      data (iri(i),i=1069,1312) /
     &          17,18,6,31,32, 17,18,6,31,32, 17,18,6,31,32,
     &          17,18,10,33,34, 17,18,6,31,32, 17,18,6,31,32,
     &          17,18,14,21,22,6,31,32,15,19,20,
     &          32,15,17,21,6,20,18,14,31,19,22,
     &          17,18,14,21,22,6,31,32,15,19,20,
     &          17,18,14,21,22,6,31,32,15,19,20,
     &          17,18,6,31,32,14,21,22,15,19,20,
     &          32,16,17,27,6,29,18,11,31,30,28,
     &          17,18,6,31,32,11,27,28,16,30,29,
     &          17,18,6,31,32,11,27,28,16,30,29,
     &          17,18,10,33,34,11,27,28,15,19,20,
     &          17,18,10,33,34,14,21,22,16,30,29,
     &          17,18,10,33,34,14,21,22,16,30,29,
     &          17,18,6,31,32,11,27,28,16,30,29,
     &          32,15,17,21,6,20,18,14,31,19,22,
     &          32,15,17,21,6,20,18,14,31,19,22,
     &          17,18,6,31,32,11,27,28,16,30,29,
     &          23,24,4,35,36,5,37,38,6,39,40,
     &          6,4,23,5,39,35,37,24,36,38,40,1,6,4,23,5,39,35,37,
     &          24,36,38,40,1,6,4,23,5,39,35,37,24,36,38,40,1,6,4/

      data (iri(i),i=1313,1563) /
     &          23,5,39,35,37,24,36,38,40,
     &          1,23,24,23,24,4,4,35,36,35,36,5,5,37,38,37,38,6,6,
     &          39,40,39,40,
     &          23,24,4,35,36,5,37,38,6,39,40,
     &          1,23,24,23,24,4,4,35,36,35,36,5,5,37,38,37,38,6,6,
     &          39,40,39,40,
     &          23,24,10,41,42,2,43,44,5,37,38,
     &          6,4,23,5,39,35,37,24,36,38,40,1,6,4,23,5,39,35,37,
     &          24,36,38,40,1,6,4,23,5,39,35,37,24,36,38,40,1,6,4,
     &          23,5,39,35,37,24,36,38,40,
     &          1,1,1,23,24,23,24,23,24,23,24,2,2,2,2,43,44,43,44,
     &          43,44,43,44,3,3,3,3,45,46,45,46,45,46,45,46,6,6,6,
     &          6,39,40,39,40,39,40,39,40,
     &          1,23,24,23,24,10,10,41,42,41,42,2,2,43,44,43,44,5,
     &          5,37,38,37,38,
     &          23,24,10,41,42,2,43,44,5,37,38,
     &          1,23,24,23,24,10,10,41,42,41,42,2,2,43,44,43,44,5,
     &          5,37,38,37,38,
     &          23,24,7,47,48,4,35,36,15,25,26,5,37,38,14,49,50,6,
     &          39,40,8,51,52/

      data (iri(i),i=1564,1869) /
     &          6,5,4,23,35,39,37,24,38,36,40,14,15,8,7,49,47,25,51,
     &          50,52,48,26,1,6,5,4,23,35,39,37,24,38,36,40,14,15,8,
     &          7,49,47,25,51,50,52,48,26,1,6,5,4,23,35,39,37,24,38,
     &          36,40,14,15,8,7,49,47,25,51,50,52,48,26,1,6,5,4,23,
     &          35,39,37,24,38,36,40,14,15,8,7,49,47,25,51,50,52,48,26,
     &          1,1,1,23,24,23,24,23,24,23,24,7,7,7,7,47,48,47,48,
     &          47,48,47,48,4,4,4,4,35,36,35,36,35,36,35,36,15,15,
     &          15,15,25,26,25,26,25,26,25,26,5,5,5,5,37,38,37,38,
     &          37,38,37,38,14,14,14,14,49,50,49,50,49,50,49,50,6,
     &          6,6,6,39,40,39,40,39,40,39,40,8,8,8,8,51,52,51,52,
     &          51,52,51,52,
     &          7,4,23,6,14,47,8,5,39,15,51,35,25,37,49,24,48,50,36,
     &          26,38,40,52,1,7,4,23,6,14,47,8,5,39,15,51,35,25,37,
     &          49,24,48,50,36,26,38,40,52,
     &          23,24,7,47,48,4,35,36,15,25,26,5,37,38,14,49,50,6,
     &          39,40,8,51,52,
     &          23,24,7,47,48,4,35,36,15,25,26,5,37,38,14,49,50,6,
     &          39,40,8,51,52,
     &          23,24,12,53,54,4,35,36,11,55,56,5,37,38,16,57,58,6,
     &          39,40,13,59,60/

      data (iri(i),i=1870,2176) /
     &          1,1,1,23,24,23,24,23,24,23,24,12,12,12,12,53,54,53,
     &          54,53,54,53,54,4,4,4,4,35,36,35,36,35,36,35,36,11,11,
     &          11,11,55,56,55,56,55,56,55,56,5,5,5,5,37,38,37,38,37,
     &          38,37,38,16,16,16,16,57,58,57,58,57,58,57,58,6,6,6,6,
     &          39,40,39,40,39,40,39,40,13,13,13,13,59,60,59,60,59,
     &          60,59,60,
     &          1,23,24,23,24,12,12,53,54,53,54,4,4,35,36,35,36,11,11,
     &          55,56,55,56,5,5,37,38,37,38,16,16,57,58,57,58,6,6,39,
     &          40,39,40,13,13,59,60,59,60,
     &          23,24,4,35,36,5,37,38,6,39,40,16,57,58,12,53,54,13,59,
     &          60,11,55,56,
     &          1,1,1,23,24,23,24,23,24,23,24,4,4,4,4,35,36,35,36,35,
     &          36,35,36,5,5,5,5,37,38,37,38,37,38,37,38,6,6,6,6,39,
     &          40,39,40,39,40,39,40,16,16,16,16,57,58,57,58,57,58,57,
     &          58,12,12,12,12,53,54,53,54,53,54,53,54,13,13,13,13,59,
     &          60,59,60,59,60,59,60,11,11,11,11,55,56,55,56,55,56,55,
     &          56,
     &          1,23,24,23,24,4,4,35,36,35,36,5,5,37,38,37,38,6,6,39,
     &          40,39,40,16,16,57,58,57,58,12,12,53,54,53,54,13,13,59,
     &          60,59,60,11,11,55,56,55,56/

      data (iri(i),i=2177,2458) /
     &          23,24,10,41,42,2,43,44,5,37,38,16,57,58,14,49,50,8,51,
     &          52,13,59,60,
     &          23,24,10,41,42,2,43,44,5,37,38,16,57,58,14,49,50,8,51,
     &          52,13,59,60,
     &          23,24,10,41,42,2,43,44,5,37,38,16,57,58,14,49,50,8,51,
     &          52,13,59,60,
     &          23,24,10,41,42,2,43,44,5,37,38,16,57,58,14,49,50,8,51,
     &          52,13,59,60,
     &          1,1,1,23,24,23,24,23,24,23,24,10,10,10,10,41,42,41,42,
     &          41,42,41,42,2,2,2,2,43,44,43,44,43,44,43,44,5,5,5,5,
     &          37,38,37,38,37,38,37,38,16,16,16,16,57,58,57,58,57,58,
     &          57,58,14,14,14,14,49,50,49,50,49,50,49,50,8,8,8,8,51,
     &          52,51,52,51,52,51,52,13,13,13,13,59,60,59,60,59,60,59,
     &          60,
     &          1,1,1,23,24,23,24,23,24,23,24,10,10,10,10,41,42,41,42,
     &          41,42,41,42,2,2,2,2,43,44,43,44,43,44,43,44,5,5,5,5,
     &          37,38,37,38,37,38,37,38,16,16,16,16,57,58,57,58,57,58,
     &          57,58,14,14,14,14,49,50,49,50,49,50,49,50,8,8,8,8,51,
     &          52,51,52,51,52,51,52,13,13,13,13,59,60,59,60,59,60,59,
     &          60/

      data (iri(i),i=2459,2742) /
     &          23,24,57,16,58,5,37,38,53,12,54,2,43,44,47,7,48,10,41,
     &          42,49,14,50,1,23,24,57,16,58,5,37,38,53,12,54,2,43,44,
     &          47,7,48,10,41,42,49,14,50,1,23,24,57,16,58,5,37,38,53,
     &          12,54,2,43,44,47,7,48,10,41,42,49,14,50,1,23,24,57,16,
     &          58,5,37,38,53,12,54,2,43,44,47,7,48,10,41,42,49,14,50,
     &          23,24,57,16,58,4,35,36,59,13,60,5,37,38,53,12,54,6,39,
     &          40,55,11,56,1,1,1,4,4,4,5,5,5,6,6,6,23,23,23,35,35,35,
     &          37,37,37,39,39,39,24,24,24,36,36,36,38,38,38,40,40,40,
     &          57,57,57,59,59,59,53,53,53,55,55,55,16,16,16,13,13,13,
     &          12,12,12,11,11,11,58,58,58,60,60,60,54,54,54,56,56,56,
     &          1,23,24,23,24,10,10,41,42,41,42,2,2,43,44,43,44,5,5,
     &          37,38,37,38,16,16,57,58,57,58,14,14,49,50,49,50,8,8,
     &          51,52,51,52,13,13,59,60,59,60,
     &          1,23,24,23,24,10,10,41,42,41,42,2,2,43,44,43,44,5,5,
     &          37,38,37,38,16,16,57,58,57,58,14,14,49,50,49,50,8,8,
     &          51,52,51,52,13,13,59,60,59,60/

      data (iri(i),i=2743,2822) /
     &          5, 23,24, 23,24, 
     &          1,1,17,18,17,18,17,18,14,14,14,21,22,21,22,21,22,
     &          23,24,16,57,58, 24,16,23,57,58, 23,24,16,57,58,
     &          23,24,16,57,58, 6,4,5, 4,6,5, 6,3,2, 6,4,5,
     &          6,5,4,15,14,8,7, 6,2,3,16,11,7,8, 6,7,8,2,3,16,11,
     &          17,18,20,19,15/

      data nop  /0,0,1,1,3,3,3,3,3,3,3,3,3,1,3,3,1,3,1,1,1,3,1,7,
     &           7,15,7,7,3,3,3,3,3,3,3,7,7,7,7,7,7,
     &           7,15,15,7,7,7,7,3,3,3,3,3,3,3,3,3,3,3,
     &           7,7,7,7,7,7,15,15,7,7,7,7,3,3,3,7,7,
     &           3,7,3,3,7,3,7,7,7,7,7,7,7,7,7,15,15,7,
     &           7,7,7,7,7,15,15,15,15,7,7,7,7,7,7,7,7,
     &           15,15,15,15,7,7,7,7,7,7,15,7,7,7,7,7,7,7,7,7,
     &           15,15,15,15,2,2,2,8,2,8,5,5,5,5,5,5,5,5,5,5,17,17,
     &           5,5,5,5,17,17,5,5,5,5,5,5,5,5,5,11,11,11,11,11,11,
     &           11,11,11,11,11,11,11,11,11,11,47,23,11,23,11,47,47,
     &           23,11,23,23,95,95,47,23,23,23,95,47,23,95,47,23,
     &           23,23,23,95,95,95,95,47,47,1,2,2,17,5,5,5,5,
     &           3,3,3,3,7,7,7,5/

      data icen /0,1,2,1,1,2,1,1,2,1,2,2,2,1,1,2,1,2,2,2,2,2,1,2,
     &           2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
     &           2, 2, 2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,
     &           1,1,1,1,1,1, 1, 1,1,1,1,1,2,2,2,2,2,
     &           2,2,1,1,2,1,1,1,2,2,2,2,2,2,2,2,2,2,
     &           2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
     &           2,2,2,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,
     &           1,1,1,1,1,2,2,2,2,1,1,2,2,2,2,2,2,2,
     &           2,2,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2,2,
     &           1,1,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,2,
     &           2,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2,2,2,
     &           2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,2,1,2,2,2,1,1,
     &           2,2,2,1,2,2,2,2/

      data ipt  /0,0,1,2,3,6,9,12,15,18,21,24,27,30,31,34,37,38,41,
     &           42,43,44,47,48,55,62,77,84,91,94,97,100,103,106,109,
     &           112,119,126,133,140,147,154,161,176,191,198,205,212,
     &           219,222,225,228,231,234,237,240,243,246,249,252,259,
     &           266,273,280,287,294,309,324,331,338,345,352,355,358,
     &           361,368,375,378,385,388,391,398,401,408,415,422,429,
     &           436,443,450,457,464,479,494,501,508,515,522,529,536,
     &           551,566,581,596,603,610,617,624,631,638,645,652,667,
     &           682,697,712,719,726,733,740,747,754,769,776,783,790,
     &           797,804,811,818,825,832,847,862,877,892,894,896,898,
     &           906,908,916,921,926,931,936,941,946,951,956,961,966,
     &           983,1000,1005,1010,1015,1020,1037,1054,1059,1064,
     &           1069,1074,1079,1084,1089,1094,1099,1110,1121,1132,
     &           1143,1154,1165,1176,1187,1198,1209,1220,1231,1242,
     &           1253,1264,1275,1322,1345,1356,1379,1390,1437,1484,
     &           1507,1518,1541,1564,1659,1754,1801,1824,1847,1870,
     &           1965,2012,2035,2130,2177,2200,2223,2246,2269,2364,
     &           2459,2554,2649,2696,2743,2744,2746,2748,2765,2770,
     &           2775,2780,2785,2788,2791,2794,2797,2804,2811,2818/

c a,b,c
      data (jset(i,1),i=1,3)          /1,2,3/
c c,b,a
      data (jset(i,2),i=1,3)          /3,2,1/
c a,c,-b
      data (jset(i,3),i=1,3)          /1,3,-2/
c -b,c,a
      data (jset(i,4),i=1,3)          /-2,3,1/
c c,a,b
      data (jset(i,5),i=1,3)          /3,1,2/
c b,c,a
      data (jset(i,6),i=1,3)          /2,3,1/
c a,-c,b
      data (jset(i,7),i=1,3)          /1,-3,2/
c b,a,-c
      data (jset(i,8),i=1,3)          /2,1,-3/
c -c,b,a
      data (jset(i,9),i=1,3)          /-3,2,1/

      data (iset(i,1),i=1,6) / 1,1,1,1,1,1 /
      data (iset(i,2),i=1,6) / 1,1,1,1,1,1 /
      data (iset(i,3),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,4),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,5),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,6),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,7),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,8),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,9),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,10),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,11),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,12),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,13),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,14),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,15),i=1,6) / 1,2,3,4,1,1 /
      data (iset(i,16),i=1,6) / 1,1,1,1,1,1 /
      data (iset(i,17),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,18),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,19),i=1,6) / 1,1,1,1,1,1 /
      data (iset(i,20),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,21),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,22),i=1,6) / 1,1,1,1,1,1 /
      data (iset(i,23),i=1,6) / 1,1,1,1,1,1 /
      data (iset(i,24),i=1,6) / 1,1,1,1,1,1 /
      data (iset(i,25),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,26),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,27),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,28),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,29),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,30),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,31),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,32),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,33),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,34),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,35),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,36),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,37),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,38),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,39),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,40),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,41),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,42),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,43),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,44),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,45),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,46),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,47),i=1,6) / 1,1,1,1,1,1 /
      data (iset(i,48),i=1,6) / 1,1,1,1,1,1 /
      data (iset(i,49),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,50),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,51),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,52),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,53),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,54),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,55),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,56),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,57),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,58),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,59),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,60),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,61),i=1,6) / 1,1,6,1,1,1 /
      data (iset(i,62),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,63),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,64),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,65),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,66),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,67),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,68),i=1,6) / 1,5,6,7,8,9 /
      data (iset(i,69),i=1,6) / 1,1,1,1,1,1 /
      data (iset(i,70),i=1,6) / 1,1,1,1,1,1 /
      data (iset(i,71),i=1,6) / 1,1,1,1,1,1 /
      data (iset(i,72),i=1,6) / 1,5,6,1,1,1 /
      data (iset(i,73),i=1,6) / 1,1,6,1,1,1 /
      data (iset(i,74),i=1,6) / 1,5,6,7,8,9 /
c still have to do tetragonal 75-142
      data ((iset(i,j),i=1,6),j=75,mxspg) / mxspg6*1 /

      data emass /1.008d0,4.003d0,6.939d0,9.012d0,10.811d0,12.011d0,
     & 14.007d0,16.000d0,19.000d0,20.183d0,22.990d0,24.312d0,26.982d0,
     & 28.086d0,30.974d0,32.064d0,35.453d0,39.948d0,82*0.0d0/

      data cutt /1.15d0,1.3d0,1.2d0,1.15d0,
     &           1.3d0,1.6d0,1.6d0,1.5d0,
     &           1.2d0,1.6d0,1.4d0,1.3d0,
     &           1.15d0,1.5d0,1.3d0,1.3d0/

      data iptyp /1,1,6,7,8,15,16,17/
      data labtyp /'HY','HP','CA','NI','OX','PH','SU','CL'/

c aij unit = eV
c bij unit = 1/angstrom
c cij unit = eV/(angstrom*6)

c
c Chlorine: Ley-Yeh Hsu and Donald E. Williams, Acta Cryst. 
c          (1980). A36,277-281
c
c Phosfor no validated data
c
      data aij /124.07167d0,4*0.0d0,3832.147d0,2638.0285d0,
     &          2384.4658d0,6*0.0d0,9583.327d0,0.0d0,9583.327d0,
     &          83*0.0d0/
      data bij /3.740d0,4*0.0d0,3.600d0,3.780d0,
     &          3.960d0,6*0.0d0,3.510d0,0.0d0,3.510d0,
     &          83*0.0d0/
      data cij /1.413698d0,4*0.0d0,25.286949d0,14.286224d0,
     &          11.645288d0,6*0.0d0,80.222339d0,0.0d0,80.222339d0,
     &          83*0.0d0/

c     hp met hy hp ca ni ox ph su cl
      data (alkp(i),i=1,mxtyp) /
     &         94.99328416d0,72.72993021d0,527.9316116d0,
     &         438.0224024d0,416.4397106d0,
     &         834.86209d0,0.0d0,834.86209d0/

c     hp met hy hp ca ni ox ph su cl
      data (blkp(i),i=1,mxtyp) /
     &         0.238095238d0,0.214592274d0,0.242130751d0,
     &         0.236966825d0,0.232018561d0,
     &         0.25940337d0,0.0d0,0.25940337d0/

c     hp met hy hp ca ni ox ph su cl
      data (clkp(i),i=1,mxtyp) /
     &         0.489668882d0,0.169608822d0,2.070963451d0,
     &         1.556621244d0,1.405397999d0,
     &         3.6830222d0,0.0d0,3.6830222d0/


      toang = 0.52917706d0

      istat = 2

c record type 1, Directory Information

      call nxtlin(line,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 1000
      read(line,1010,end=1000,err=1000) tstr,refcod,isys,ncards,
     &                                nrfac,nrem,ndis,nerr,nopr,nrad,
     &                                nat,nsat,ncon,icell,iatfor,icent

      if (tstr.ne.'#') goto 1000
      if (nat.eq.0) goto 1000
      ntot = nrfac+nrem+ndis+nerr
      ncomm = ntot/80
      if (ntot-ncomm*80.gt.0) ncomm = ncomm + 1

      print*,'Refcode            = ',refcod
      print*,'Number of Atoms    = ',nat
      print*,'Number of SAtoms   = ',nsat
      print*,'Centre of Symmetry = ',icent
      if (idebug.eq.1) then
         print*,'Number of connectivity integers = ',ncon
      endif

c record type 2, Unit Cell Parameters

      inorm = 0
      if (icell.gt.0) then
         torad = datan(1.0d0) / 45.0d0
         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 1000
         read(line,1020,end=1000,err=1000) (icelld(i),i=1,6),
     &       (iprec(i),i=1,6),(icesd(i),i=1,6),idm,idx,nspg,spg,kz,itol
         a = dfloat(icelld(1))/10**iprec(1)
         b = dfloat(icelld(2))/10**iprec(2)
         c = dfloat(icelld(3))/10**iprec(3)
         alpha = dfloat(icelld(4))/10**iprec(4)
         beta  = dfloat(icelld(5))/10**iprec(5)
         gamma = dfloat(icelld(6))/10**iprec(6)
         dm = dfloat(idm)/100
         dx = dfloat(idx)/100
         tol = dfloat(itol)/100

c check for alternative settings
         
         if (nspg.eq.14.and.spg(5:5).eq.'n'.or.spg(5:5).eq.'N') 
     &      nspg = 231

         if (nspg.eq.146.and.nopr.eq.9) nspg = 232
         if (nspg.eq.148.and.nopr.eq.9) nspg = 233
         if (nspg.eq.155.and.nopr.eq.18) nspg = 234
         if (nspg.eq.160.and.nopr.eq.18) nspg = 235
         if (nspg.eq.161.and.nopr.eq.18) nspg = 236
         if (nspg.eq.166.and.nopr.eq.18) nspg = 237
         if (nspg.eq.167.and.nopr.eq.18) nspg = 238

         call prcell(nspg,a,b,c,alpha,beta,gamma)

         print*,' '
         print*,'No. of formula units per cell = ',kz

         call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)

      endif


c record type 3, Text Information

      do i=1,ncomm
         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 1000
      end do

c record type 4, Symmetry Positions

      ncomm = nopr/5
      if (nopr-ncomm*5.gt.0) ncomm = ncomm + 1

      ntot = 0
      do l=1,ncomm
         jleft = nopr - ntot
         if (jleft.ge.5) jleft = 5
         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 1000
         read(line,1040,end=1000,err=1000) 
     &     (((ir(i,j,ntot+k),j=1,3),it(i,ntot+k),i=1,3),k=1,jleft)
         ntot = ntot + 5
      end do

      if (idebug.eq.1) call prop(nopr,ir,it)
      if (idebug.eq.1) call symrec(nopr,icent,ir,it,iun3,.true.)

c      call cprop(nopr,ir,it)

c record type 5, Radius Values

      if (nrad.gt.0) then
         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 1000
      endif
c      read(line,1050,end=1000,err=1000) (atom(i),irad(i),i=1,nrad)

c record type 6, Atomic Coordinates

      nat = nat+nsat
      norg = nat
      nsum = nat
      iatoms = nsum
      if (iatfor.gt.0) then
         k = 1
         do while (.true.)
             jleft = iatoms - k + 1
             if (jleft.ge.3) jleft = 3
             call nxtlin(line,jstat)
             if (jstat.eq.1.or.jstat.eq.2) goto 1000
             read(line,1060,end=1000,err=1000) 
     &           (atmtmp(i),(ii(i,j),j=1,3),i=1,jleft)
             do i=1,jleft
                l = k-1+i
                atom = atmtmp(i)(1:2)
                j = ichar(atom(2:2))
                if (j.lt.65.or.j.gt.122.or.(j.lt.97.and.j.gt.90)) 
     &          then
                   atom(2:2) = atom(1:1)
                   atom(1:1) = ' '
                endif
                do j=1,99
                    if (tocapf(atom).eq.tocapf(elemnt(j))) ianz(l) = j
                end do
                iatclr(l) = 1

                do j=1,3
                   coo(j,l) = dfloat(ii(i,j))/100000
                end do
             end do
             k = k + 3
             if (k.gt.iatoms) goto 50
         end do
50       continue
      endif

      if (iatoms.eq.0) goto 1000
c record type 7, Connection Table

      if (nat.le.100) then
          fmt = '(40i2)'
          num = 40
      else
          fmt = '(26i3,2x)'
          num = 26
      endif

      do i=1,iatoms
         iconn(1,i) = 0
      end do

      if (ncon.le.0) goto 200

      jcon = 1
      isave = 0
      do while (.true.)
         jleft = ncon - jcon + 1
         if (jleft.ge.num) jleft = num
         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 1000
         read(line,fmt,end=1000,err=1000) (icon(i),i=1,jleft)
         i = 1
         do while (.true.)
            if (i.gt.jleft) goto 100
            if (jcon.le.nsum) then
               iat = jcon
               jat = icon(i)
               i = i + 1
               jcon = jcon + 1
            else
               if (i.eq.jleft) then
                  isave = icon(i)
                  i = i + 1
               else
                  if (isave.ne.0) then
                     iat = isave
                     jat = icon(i)
                     i = i + 1
                     isave = 0
                  else
                     iat = icon(i)
                     jat = icon(i+1)
                     i = i + 2
                  endif
                  jcon = jcon + 2
               endif
            endif
            if (jat.ne.0) then
               if (iconn(1,iat).lt.mxcon) then
                  iconn(1,iat) = iconn(1,iat) + 1
                  iconn(iconn(1,iat)+1,iat) = jat
               else
                  write(iun3,*) 'more than mxconn connections found'
               endif
               if (iconn(1,jat).lt.mxcon) then
                  iconn(1,jat) = iconn(1,jat) + 1
                  iconn(iconn(1,jat)+1,jat) = iat
               else
                  write(iun3,*) 'more than mxconn connections found'
               endif
            endif
         end do
100      continue
         if (jcon.gt.ncon) goto 200
      end do
200   continue
      call redcon

      do i=1,iatoms
         if (inorm.eq.0) call fr2crt(coo(1,i),xa,ya,yb,za,zb,zc)
         do j=1,3
             coo(j,i) = coo(j,i)/toang
         end do
      end do

      call dohcon(0)

      do i=1,iatoms
         do j=1,3
            tr1(j) = coo(j,i) * toang
         end do
         call crt2fr(tr1,coo(1,i),xa,ya,yb,za,zb,zc)
      end do

c
c     Copy original atoms to top of coordinates array
c
      nstor = mxnat-iatoms
      do i=1,iatoms
         do j=1,3
            coo(j,nstor+i) = coo(j,i)
         end do
         ianz(nstor+i) = ianz(i)
         do j=1,iconn(1,i)+1
            iconn(j,nstor+i) = iconn(j,i)
         end do
         iatclr(nstor+i) = iatclr(i)
      end do

      if (idebug.eq.1) then
         do i=1,iatoms
            print*,i,' ',elemnt(ianz(i)),iconn(1,i),
     &               (iconn(1+j,i),j=1,iconn(1,i))
         end do
      endif

      return
1010  format(a1,a8,i1,13x,i3,8i3,3x,i3,i1,1x,2i1)
1020  format(6i6,6i1,6i2,2i3,i3,a8,i3,i2,4x)
1040  format(5(3i1,i2,3i1,i2,3i1,i2),5x)
1050  format(16(a2,i3))
1060  format(a5,3i7,1x,a5,3i7,1x,a5,3i7)
1070  format(40i2)
1075  format(26i3,2x)

1000  istat = 0
      iatoms = 0
      return
      end

      subroutine fdad(iop,ifrst,istdbd,iuseab,moddma,idebug,
     &                coo,qat,ianz,iaton,iatclr,iresid,iconn,
     &                lwrit,lring,ityp,scal,scali,smag,
     &                nat,norg,icent,inorm,ncon,nspg,nopr,ir,it,
     &                xa,ya,yb,za,zb,zc,a,b,c)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxdma=5000)
      parameter (mxel=100)
      parameter (mtel=8)
      parameter (mxtyp=8)
      parameter (mxsite=300)
      parameter (lmax=4)
      parameter (lmaxf=(lmax+1)*(lmax+1))
      parameter (mxplev=5)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /athlp/  iatoms, mxnat
      character*8   ctag
      common /multip/ qmom(lmaxf,mxsite),car(3,mxsite),ctag(mxsite),
     &                nsites
      common /mulmap/ mulivt(lmaxf)
      character*2 elemnt
      common /elem/   elemnt(mxel)
      common /exmass/ emass(mxel),cutt(4,4)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /surf/   natorg,noscnd
      integer*2 ir,it
      common /cllhlp/ ico(3,8),icn(4,8),mcol(32)
      common /cllab/  iclon,iclpnt(4)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      integer*2 ityp
      character*2 labtyp
      common /buck/ alk(mxtyp,mxtyp),blk(mxtyp,mxtyp),clk(mxtyp,mxtyp),
     &              alkp(mxtyp),blkp(mxtyp),clkp(mxtyp),
     &              iptyp(mxtyp),labtyp(mxtyp)
      common /potlev/grdw, plevel(mxplev),ipcol(mxplev+1),nplev
      character*4 label,lab1,lab2,lab3,ends,dmalab
      character*80 tnkrt
      character*137 string
      logical applop,symmex,ldupl,rescel,doinv,addcel,expcel,
     &        doivrt,multok,ringg,obenz,doint,opfil
      integer getlin
      dimension tr1(3),tr2(3),vec(3),vect(3),veci(3),iorg(mxdma),
     &          iacn(mxdma),ido(mxdma),iatel(mtel),iytyp(mxtyp),
     &          iring(6),i3(3),j3(3),qrot(13),iaxdef(mxsite,3)
      dimension coo(3,*),qat(*),ianz(*),iaton(*),iatclr(*),
     &          iconn(mxcon+1,*),ityp(*),iresid(*),lwrit(*),lring(*)
      dimension ir(3,3,192),it(3,192),tmp(3)
      data mulivt / 1,
     &              1, 1,-1,
     &              1, 1,-1, 1,-1,
     &              1, 1,-1, 1,-1, 1,-1,
     &              1, 1,-1, 1,-1, 1,-1, 1,-1/
      
      if (nat.gt.mxdma) return

      doint = .false.

      do i=1,mtel
         iatel(i) = 0
      end do

      do i=1,mxtyp
         iytyp(i) = 0
      end do

      idun = 25
      ends = 'ENDS'
      doivrt = .false.
      multok = .true.
      tol = 1.0d-6
      ctol = 1.0d-6

      iclon = 0
      toang = 0.52917706d0

      addcel = .true.
      rescel = .false.
      symmex = .true.
      expcel = .false.

      if (iop.eq.1.or.iop.eq.15.or.iop.eq.16.or.iop.eq.17) then
         addcel = .false.
         symmex = .false.
      elseif (iop.eq.2) then
         symmex = .false.
      elseif (iop.eq.4.or.iop.eq.20) then
         rescel = .true.
      elseif (iop.eq.6) then
         rescel = .true.
      elseif (iop.ge.7) then
         expcel = .true.
      endif

      if (expcel) then
         iamin = 0
         iamax = 0
         ibmin = 0
         ibmax = 0
         icmin = 0
         icmax = 0
         if (iop.eq.7.or.iop.eq.8.or.iop.eq.11.or.iop.eq.12) then
          iamin = -1
          iamax =  1
         endif
         if (iop.eq.7.or.iop.eq.9.or.iop.eq.11.or.iop.eq.13) then
          ibmin = -1
          ibmax =  1
         endif
         if (iop.eq.7.or.iop.eq.10.or.iop.eq.12.or.iop.eq.13) then
          icmin = -1
          icmax =  1
         endif
      endif
      tl = 1.0d-4

      nsum = nat
      iatoms = nsum

c
c     Copy original atoms from top of coordinates array
c
c     iorg = nr. original atom
c     iacn = counter for each type of atom
c     iatclr = molecule number, or nr. of symmetry copy
c              if neg. inversion has taken place
c
      nstor = mxnat-iatoms
      do i=1,nsum
         do j=1,3
            coo(j,i) = coo(j,nstor+i)
            if (abs(coo(j,i)).lt.ctol) coo(j,i) = 0.0d0
         end do
         ianz(i) = ianz(nstor+i)
         iorg(i) = i
         ido(i) = 1
         iacn(i) = itell(iatel,ianz(i))
         k = 0
         do j=1,iconn(1,nstor+i)
c            if (iconn(j+1,nstor+i).gt.0) then
               iconn(j+1,i) = iconn(j+1,nstor+i)
               k = k + 1
c            endif
         end do
         iconn(1,i) = k
         iatclr(i) = iatclr(nstor+i)
      end do


      ncol = 2
      call cntvec(vec,coo,ianz,iatoms)

      if (symmex) then

         do j=2,nopr

            if (applop(ir(1,1,j),it(1,j),vec,vect,0,invrt)) then
               icol1 = ncol
               ncol = ncol + 1
               call tvec(vect,tr1,idebug)

               do i=1,nsum
                   if (applop(ir(1,1,j),it(1,j),coo(1,i),coo(1,iatoms+1)
     &                 ,0,invrt)) idum = 1
                   iatoms = iatoms + 1
                   ianz(iatoms) = ianz(i)
                   iatclr(iatoms) = icol1*invrt

                   call trcoo(tr1,coo(1,iatoms))

                   if (iatoms.lt.mxdma) then
                      iorg(iatoms) = i
                      iacn(iatoms) = iacn(i)
                   endif

                   qat(iatoms) = qat(i)
                   ityp(iatoms) = ityp(i)

               end do
            endif

         end do

         do j=1,nopr
            doinv = .false.
            if (applop(ir(1,1,j),it(1,j),vec,vect,0,invrt)) idum = 1

            if (icent.eq.1) then
               sum = 0.0d0
               do i=1,3
                  veci(i) = 0.0d0
               end do
               do i=1,nsum
                  if (applop(ir(1,1,j),it(1,j),coo(1,i),tmp,1,invrt)) 
     &                idum = 1
                  wat = dble(ianz(i))
                  do k=1,3
                     veci(k) = veci(k) + wat*tmp(k)
                  end do
               end do
		
               if (sum.gt.0.0d0) then
                  do k=1,3
                     veci(k) = veci(k) / sum
                  end do
               endif

               do k=1,3
                  tr2(k) = gsmod(veci(k) - vect(k))
               end do

               if (vlen(tr2).gt.tl) then
                  doinv = .true.
                  icol2 = ncol
                  ncol = ncol + 1
                  if (applop(ir(1,1,j),it(1,j),vec,veci,1,invrt)) then
                     do k=1,3
                        tr2(k) = gsmod(veci(k) - vect(k))
                     end do
                  endif
                  call tvec(veci,tr2,idebug)
               endif
            endif

            if (doinv) then

               do i=1,nsum
                   if (applop(ir(1,1,j),it(1,j),coo(1,i),
     &                 coo(1,iatoms+1),1,invrt)) then
                      iatoms = iatoms + 1
                      ianz(iatoms) = ianz(i)
                      iatclr(iatoms) = icol2*invrt
                      call trcoo(tr2,coo(1,iatoms))
                      if (iatoms.lt.mxdma) then
                         iorg(iatoms) = i
                         iacn(iatoms) = iacn(i)
                      endif
                      qat(iatoms) = qat(i)
                      ityp(iatoms) = ityp(i)
                   endif
               end do
            endif

         end do
c
c Put Original in the Cell
c
         call tvec(vec,tr1,idebug)
         do i=1,nsum
            call trcoo(tr1,coo(1,i))
         end do
      endif

      if (idebug.eq.1) then
         do i=1,iatoms
            print*,i,' ',elemnt(ianz(i)),iconn(1,i),
     &               (iconn(1+j,i),j=1,iconn(1,i))
         end do
      endif

      if (symmex) then
         if (rescel) then
            do i=1,iatoms
               do j=1,3
                  if (coo(j,i).gt.1.0d0.or.coo(j,i).lt.0.0d0) then
                      coo(j,i) = dmod(coo(j,i),1.0d0)
                  endif
                  if (coo(j,i).lt.0.0d0) coo(j,i) = 1.0d0 + coo(j,i)
               end do
            end do
         endif
c
c Clean up duplicates
c
60       continue
         do i=1,iatoms
            do j=i+1,iatoms
                ldupl = .true.
                do k=1,3
                   ldupl = (ldupl.and.dabs(coo(k,i)-coo(k,j)).lt.tl)
                end do
                if (ldupl) then
                   call delat(j,iacn,iorg,mxdma)
                   goto 60
                endif
            end do
         end do
      endif

      if (iop.eq.18) then
         nat = iatoms
         nsum = iatoms
         nstor = mxnat-iatoms
         do i=1,iatoms
            do j=1,3
               coo(j,nstor+i) = coo(j,i)
            end do
            ianz(nstor+i) = ianz(i)
            iatclr(nstor+i) = 1
            iatclr(i) = 1
         end do
      endif

      do i=1,iatoms
         if (inorm.eq.0) call fr2crt(coo(1,i),xa,ya,yb,za,zb,zc)
         do j=1,3
             coo(j,i) = coo(j,i)/toang
         end do
      end do

      if (symmex) then
         do i=1,iatoms
            iconn(1,i) = 0
         end do
         call doconn
         if (iop.eq.18) then
            do i=1,nsum
               do j=1,iconn(1,i)+1
                  iconn(j,nstor+i) = iconn(j,i)
               end do
            end do
         endif
      elseif (iop.eq.1.and.ncon.le.0) then
         call doconn
         do i=1,nsum
            do j=1,iconn(1,i)+1
               iconn(j,nstor+i) = iconn(j,i)
            end do
         end do
      endif

      if (iop.eq.1.and.istdbd.gt.0) then
         radian = (120.0d0/45.0d0)*datan(1.0d0)
         do i=1,nsum
             lring(i) = 0
             iaton(i) = 2
         end do
         do i=1,nsum
            if (ianz(i).eq.6.or.ianz(i).eq.7) then
               obenz = .false.
               if (ianz(i).eq.6) then
                   chsc = 1.08d0 / toang
                   if (istdbd.gt.1) then
                     if (.not.ringg(i,iring,nring,.true.,
     &                       ianz,iaton,iconn,lwrit,lring)) idum = 0
                     if (nring.eq.6) then
                        icons = 0
                        do k=1,6
                           icons = icons + iconn(1,iring(k))
                        end do
                        if (icons.eq.18) then
                           obenz = .true.
                           chsc = 1.01d0 / toang
                        endif
                     endif
                   endif
               elseif (ianz(i).eq.7) then
                   chsc = 1.01d0 / toang
               endif
               do j=1,iconn(1,i)
                  n = iconn(1+j,i)
                  if (ianz(n).eq.1) then
c
c correct x-h bondlength
c
                     do k=1,3
                        tr1(k) = coo(k,n) - coo(k,i)
                     end do 
                     call vsc1(tr1,chsc,tol)
                     do k=1,3
                        coo(k,n) = coo(k,i) + tr1(k)
                        tr1(k) = coo(k,n) * toang
                     end do
                     call crt2fr(tr1,coo(1,nstor+n),xa,ya,yb,za,zb,zc)
c
c correct benzene hcc angle and dihedral in one go
c
                     if (obenz) then
                        if (iconn(1,i).eq.3) then
                           iat1 = i
                           iat2 = iconn(2,i)
                           if (iat2.eq.n) then
                               iat2 = iconn(3,i)
                               iat3 = iconn(4,i)
                           else
                               iat3 = iconn(3,i)
                               if (iat3.eq.n) iat3 = iconn(4,i)
                           endif
                           do k=1,3
                              tr1(k) = coo(k,iat2) - coo(k,iat1)
                              tr2(k) = coo(k,iat3) - coo(k,iat1)
                           end do 
                           do k=1,3
                              vect(k) = -(tr1(k) + tr2(k))
                           end do
                           call vsc1(vect,chsc,tol)
                           do k=1,3
                              coo(k,n) = coo(k,iat1) + vect(k)
                              vect(k) = coo(k,n) * toang
                           end do
                           call crt2fr(vect,coo(1,nstor+n),
     &                                 xa,ya,yb,za,zb,zc)
                        endif
                     endif
                  endif
               end do
            endif
         end do
      endif

      if (iop.eq.15) then
c
c For DMAREL
c
         if (iuseab.eq.1) then
c
c Use Ab Initio coordinates and orient
c
            do j=1,3
               i3(j) = j
               j3(j) = j
            end do
c
c align mol2 (car) with 3 atoms (j3) of mol1 (coo)
c
            call alntwo(coo,i3,car,nsum,j3)
c
c translate center of mass mol2 to center of mass mol1
c
            call trtwo(coo,ianz,nsum,car,ianz,nsites,.false.)
         endif
         if (moddma.eq.1) then
c
c do bond centers
c
            nsumn = nsum
            do i=1,nsum
               do j=i+1,nsum
                   ttot = 0.0d0
                   do k=1,3
                      vect(k) = coo(k,j) - coo(k,i)
                      ttot = ttot + vect(k)*vect(k)
                   end do
                   dvdw = (vdwr(ianz(i))+vdwr(ianz(j))) / toang
                   dvdw2 = dvdw*dvdw
                   if (ttot.lt.dvdw2) then
                       if (nsumn.lt.mxnat) then
                          nsumn = nsumn + 1
                          if (iuseab.eq.0) then
                             do k=1,3
                               coo(k,nsumn) = coo(k,i) + vect(k)/2.0d0
                             end do
                             ianz(nsumn) = 99
                             iconn(1,nsumn) = 0
                             iatclr(nsumn) = 1
                          endif
                          iconn(2,nsumn) = i
                          iconn(3,nsumn) = j
                        else
                          write(iun3,*) 'exceded max. # of sites'
                        endif
                   endif
               end do
            end do
            if (nsumn.ne.nsites) write(iun3,*) 'problem'
         endif
c
c Move molecule in the backup in the top of the coo,iconn arrays
c to allow for bondcenters
c
         if (nsites.gt.nsum) then
             nstoro = nstor
             nstor = mxnat - nsites
             do i=1,nsum
                ianz(nstor+i) = ianz(nstoro+i)
                do j=1,iconn(1,nstoro+i)+1
                   iconn(j,nstor+i) = iconn(j,nstoro+i)
                end do
                iatclr(nstor+i) = iatclr(nstoro+i)
             end do
             do i=nstor+nsum+1,mxnat
                ianz(i) = 99
                iconn(1,i) = 0
                iatclr(i) = 1
             end do
             norg = nsum
             nsum = nsites
             nat = nsum
         else 
             norg = nsum
         endif

         if (iuseab.eq.1) then
            do i=1,nsum
               do j=1,3
                  coo(j,i) = car(j,i)
               end do
            end do
         endif

c
c cartesians back to fractional coordinates in the backup
c
         do i=1,nsum
            do j=1,3
               tr1(j) = coo(j,i) * toang
            end do
            call crt2fr(tr1,coo(1,nstor+i),xa,ya,yb,za,zb,zc)
            if (i.gt.norg) then
               iconn(2,nstor+i) = iconn(2,i)
               iconn(3,nstor+i) = iconn(3,i)
            endif
         end do
      endif

      itemp = iatoms
c
c Expand Cell
c
      if (expcel) then
         do i=iamin,iamax
            do j=ibmin,ibmax
               do k=icmin,icmax
                  if (.not.(i.eq.0.and.j.eq.0.and.k.eq.0)) then
                      istrt = iatoms
                      tr1(1) = 1.0d0*i
                      tr1(2) = 1.0d0*j
                      tr1(3) = 1.0d0*k
                      call fr2crt(tr1,xa,ya,yb,za,zb,zc)
                      do m=1,3
                          tr1(m) = tr1(m)/toang
                      end do
                      do l=1,itemp
                        if (ianz(l).ne.100) then
                         iatoms = iatoms + 1
                         do m=1,3
                            coo(m,iatoms) = coo(m,l) + tr1(m)
                         end do
                         ianz(iatoms) = ianz(l)
                         iatclr(iatoms) = iatclr(l)
                         iconn(1,iatoms) = iconn(1,l)
                         qat(iatoms) = qat(l)
                         ityp(iatoms) = ityp(l)
                         do m=1,iconn(1,l)
                            if (iconn(m+1,l).lt.0) then
                               iconn(m+1,iatoms) = iconn(m+1,l) - istrt
                            else
                               iconn(m+1,iatoms) = iconn(m+1,l) + istrt
                            endif
                         end do
                        endif
                      end do
                  endif
               end do
            end do
         end do
      endif

      if (symmex) call dohcon(0)
c
c Add Cell
c
      if (addcel) then
         call addc(coo,ianz,iconn,iatclr,iatoms,xa,ya,yb,za,zb,zc)
      else
         natorg = 0
      endif

      if (doint) then
         call calcij(0,coo,ianz,iresid,iconn,qat)
      endif

      if (iop.eq.20.or.iop.eq.16.or.iop.eq.17) then
c
c Write DMAREL / Tinker  input
c
         if (iop.eq.20) then
            open(unit=idun,file='dmain',form='formatted',
     &        status='unknown',err=1000)
         else
            open(unit=idun,file='tnk.key',form='formatted',
     &        status='unknown',err=1000)
         endif
         goto 1010
1000     idun = 6
1010     continue

         if (iop.eq.16) goto 110
c
c Check if multipoles are alright
c
         if (nsum.eq.nsites) then
            do i=1,norg
               if (ctag(i)(1:2).ne.elemnt(ianz(i))) multok = .false.        
            end do
            if (.not.multok) 
     &         call inferr('No multipole match',0)
            if (iop.eq.20) then
                if (dabs(car(1,1)-car(1,2)).gt.tol.or.
     &           dabs(car(2,1)-car(2,2)).gt.tol) multok = .false.
                if (dabs(car(2,2)-car(2,3)).gt.tol.or.
     &           dabs(car(2,1)-car(2,2)).gt.tol) multok = .false.
                if (.not.multok) 
     &           call inferr('Problem with multipole axes',0)
            endif
         else
            multok = .false.
         endif

         if (multok) then
c
c Fix any minute charge leaks
c
            totch = 0.0d0
            do i=1,nsites
                totch = totch + qmom(1,i)
            end do
            if (idebug.eq.1) print*,'total charge =',totch
            if (dabs(totch).gt.1.0d-7) then
               totch = totch / dble(nsites)
               totfx = 0.0d0
               do i=1,nsites
                   qmom(1,i) = qmom(1,i) - totch
                   totfx = totfx + qmom(1,i)
               end do
               if (idebug.eq.1) print*,'total fixed charge =',totfx
            endif

         endif

110      continue
111      format(1x,7(f10.6,1x))

      endif
         
      if (iop.eq.20) then
c
c Write DMAREL input
c
         cvect = c
         cuto = 20.0d0
         write(idun,500) 'TITLE'
         write(idun,500) 'Molden generated DMAREL input'
         write(idun,500) ends
         write(idun,500) 'DIMENSION 500000'
         write(idun,500) 'JOBT 500'
         write(idun,500) 'SCALE 0.9'
         write(idun,500) 'CHGC'
         write(idun,500) 'THERMAL'
         write(idun,'(a,2f10.6)') 'CUTO ',cvect,cuto/cvect
         write(idun,'(a,f10.6)') 'RDMA ',cuto/cvect
         write(idun,'(a,a)') 'PRIN BASI 1 PLUT 10012 MOLE 1 ',
     &             'BOND 2 TORS 2 MINI 2 GEOM 1'

         write(idun,500) 'LATT'
         do i=2,4
            write(idun,'(10x,3(f17.13,2x))') (coo(j,iatoms-8+i)
     &                                    *toang/cvect,j=1,3)
         end do
         write(idun,'(a)') 'BASI'
         do i=1,iatoms-8
            label = dmalab(nat,ianz(i),iorg(i),iacn(i),iatclr(i),iytyp)
            write(idun,'(2(a4,1x),3(f17.13,4x),i3)') label,'CODA',
     &            (coo(j,i)*toang/cvect,j=1,3),iabs(iatclr(i))
            if (iatclr(i).lt.0) then
               doivrt = .true.
               if (ido(iorg(i)).ne.-1) then
                  isite = iorg(i)
                  if (multok) then
                     write(idun,*) 'LEVEL 4'
                     write(idun,111) qmom(1,isite)
                     write(idun,111) 
     &                    (qmom(j,isite)*dble(mulivt(j)),j=2,4)
                     write(idun,111) 
     &                    (qmom(j,isite)*dble(mulivt(j)),j=5,9)
                     write(idun,111) 
     &                    (qmom(j,isite)*dble(mulivt(j)),j=10,16)
                     write(idun,111) 
     &                    (qmom(j,isite)*dble(mulivt(j)),j=17,23)
                     write(idun,111) 
     &                    (qmom(j,isite)*dble(mulivt(j)),j=24,25)
                  else
                     if (ihasq.eq.1) then
                        write(idun,*) 'LEVEL 0'
                        write(idun,111) qat(i)
                     else
                        write(idun,*) 'LEVEL 1'
                        write(idun,*) '0.0'
                        write(idun,*) '0.0 0.0 0.0'
                     endif
                  endif
                  ido(iorg(i)) = -1
               else
                  write(idun,*) 'DUPL'
               endif
            else
               if (ido(iorg(i)).gt.0) then
                  if (multok) then
                     write(idun,*) 'LEVEL 4'
                     isite = iorg(i)
                     write(idun,111) qmom(1,isite)
                     write(idun,111) (qmom(j,isite),j=2,4)
                     write(idun,111) (qmom(j,isite),j=5,9)
                     write(idun,111) (qmom(j,isite),j=10,16)
                     write(idun,111) (qmom(j,isite),j=17,23)
                     write(idun,111) (qmom(j,isite),j=24,25)
                  else
                     if (ihasq.eq.1) then
                        write(idun,*) 'LEVEL 0'
                        write(idun,111) qat(i)
                     else
                        write(idun,*) 'LEVEL 1'
                        write(idun,*) '0.0'
                        write(idun,*) '0.0 0.0 0.0'
                     endif
                  endif
                  ido(iorg(i)) = 0
               else
                  write(idun,*) 'DUPL'
               endif
            endif
         end do
         write(idun,500) ends
         write(idun,500) 'POTE'
         write(idun,500) 'SPEC'
         do i=1,nsum
            label = dmalab(nat,ianz(i),iorg(i),iacn(i),1,iytyp)
            write(idun,'(1x,2(a4,2x),f4.1,2x,f7.3)') 
     &            label,'CODA',0.0,emass(ianz(i))
         end do
         if (doivrt) then
            do i=1,nsum
               label = dmalab(nat,ianz(i),iorg(i),iacn(i),-1,iytyp)
               write(idun,'(1x,2(a4,2x),f4.1,2x,f7.3)') 
     &               label,'CODA',0.0,emass(ianz(i))
            end do
         endif
         write(idun,500) ends
         do i=1,mxtyp
            if (iytyp(i).eq.1) then
               do j=i,mxtyp
                  if (iytyp(j).eq.1) then
                     call wrbuck(idun,labtyp(i),labtyp(j),
     &                           alk(i,j),blk(i,j),clk(i,j))
                  endif
               end do
               if (moddma.eq.1) then
                     call wrbuck(idun,labtyp(i),'XX',
     &                           0.0d0,1.0d0,0.0d0)
               endif
            endif
         end do
         write(idun,500) ends
         write(idun,500) 'MOLE'
         write(idun,500) ' NFXP'
         write(idun,500) ' CUTM 10.0'
         write(idun,500) ' NBUR 20'
         write(idun,500) 'NNCU'
         invval = 1
600      continue
         if (moddma.eq.1) then
            do i=1,nsum
               if (i.gt.norg) then
                   j = abs(iconn(2,nstor+i))
                   k = abs(iconn(3,nstor+i))
                   lab1 = dmalab(nat,ianz(i),iorg(i),iacn(i),
     &                           invval,iytyp)
                   lab2 = dmalab(nat,ianz(j),iorg(j),iacn(j),
     &                           invval,iytyp)
                   lab3 = dmalab(nat,ianz(k),iorg(k),iacn(k),
     &                           invval,iytyp)
                   call getdis(vdwr(ianz(j)),vdwr(ianz(k)),
     &                        ianz(j),ianz(k),dmaxs)
                   write(idun,'(1x,4(a4,2x),f6.3)') 
     &                     lab2,'CODA',lab1,'CODA',dmaxs/2.0d0
                   write(idun,'(1x,4(a4,2x),f6.3)') 
     &                     lab3,'CODA',lab1,'CODA',dmaxs/2.0d0
               endif
            end do
         else
            do i=1,nsum
               n = iconn(1,nstor+i)
               do j=1,n
                  nb = iconn(1+j,nstor+i)
                  na = abs(nb)
                  if (na.gt.i) then
                     lab1 = dmalab(nat,ianz(i),iorg(i),iacn(i),
     &                             invval,iytyp)
                     lab2 = dmalab(nat,ianz(na),iorg(na),iacn(na),
     &                             invval,iytyp)
                     call getdis(vdwr(ianz(i)),vdwr(ianz(na)),
     &                        ianz(i),ianz(na),dmaxs)
                     if (nb.lt.0) dmaxs = 2.0d0
                     write(idun,'(1x,4(a4,2x),f6.3)') 
     &                     lab1,'CODA',lab2,'CODA',dmaxs
                  endif
               end do
            end do
         endif
         if (invval.eq.1.and.doivrt) then
            invval = -1
            goto 600
         endif
         write(idun,500) ends
         iaxes = 1
         if (doivrt) iaxes = 2
         write(idun,'(a,1x,i1)') 'MOLX',iaxes
         lab1 = dmalab(nat,ianz(1),iorg(1),iacn(1),1,iytyp)
         lab2 = dmalab(nat,ianz(2),iorg(2),iacn(2),1,iytyp)
         lab3 = dmalab(nat,ianz(3),iorg(3),iacn(3),1,iytyp)
         call ifnn(ic1,nstor,1,2)
         call ifnn(ic2,nstor,1,3)
         if (moddma.eq.1) then
            ic1 = ic1 * 2
            ic2 = ic2 * 2
         endif
         write(idun,'(5(a,1x),i2)') 
     &           'Z LINE',lab1,'CODA',lab2,'CODA',ic1
         write(idun,'(5(a,1x),i2,1x,a,1x,a,1x,i2)') 
     &           'X PLANE',lab1,'CODA',lab3,'CODA',ic2,lab2,'CODA',ic1
         if (doivrt) then
            lab1 = dmalab(nat,ianz(1),iorg(1),iacn(1),-1,iytyp)
            lab2 = dmalab(nat,ianz(2),iorg(2),iacn(2),-1,iytyp)
            lab3 = dmalab(nat,ianz(3),iorg(3),iacn(3),-1,iytyp)
            write(idun,'(5(a,1x),i2)') 
     &           'Z LINE',lab1,'CODA',lab2,'CODA',ic1
            write(idun,'(5(a,1x),i2,1x,a,1x,a,1x,i2)') 
     &           'X PLANE',lab1,'CODA',lab3,'CODA',ic2,lab2,'CODA',ic1
         endif
         write(idun,500) ends
         write(idun,500) 'STAR PLUT'
         write(idun,500) 'CONP'
         write(idun,500) 'VIDEO 100'
         write(idun,500) 'MAXI 300'
         write(idun,500) 'MAXD 0.5'
         write(idun,500) 'MAXT 0.5'
         write(idun,500) 'START'
         write(idun,500) 'STOP'

         if (idun.ne.6) then
             close(idun)
             call inferr('Wrote file: '//'dmain',0)
         endif
      endif

      if (iop.eq.17) then

         do i=1,mxsite
            iaxdef(i,1) = 0
         end do

         if (opfil(46,'axes.def',8,1,1,1)) then
            iuntmp = iun2
            iun2 = 46
            do while (getlin(0).eq.1)
               ktype = nxtwrd(string,nstr,itype,rtype)
               if (ktype.eq.2) then
                  iax = itype
                  ktype = nxtwrd(string,nstr,itype,rtype)
                  if (ktype.eq.2) iaxdef(iax,1) = itype
                  ktype = nxtwrd(string,nstr,itype,rtype)
                  if (ktype.eq.2) then
                     iaxdef(iax,2) = iabs(itype)
                     if (itype.lt.0) then
                        iaxdef(iax,3) = 2
                     else
                        iaxdef(iax,3) = 1
                     endif
                  endif
               endif
            end do
            close(46)
            iun2 = iuntmp
         endif
      endif

      if (iop.eq.16.or.iop.eq.17) then

         call getenv('XTINKER',tnkrt)
         if (tnkrt(1:1).eq.' ') then
            call getenv('TNK_ROOT',tnkrt)
         endif
         itklen = linlen(tnkrt)
         write(idun,'(a)') ' '
         write(idun,'(a)') 
     &     'parameters '//tnkrt(1:itklen)
     &    //'/params/mm3dipoleintra.prm'
         write(idun,'(a)') ' '
         ierr = 0

      endif

      if (iop.eq.17) then

         if (multok) then
             do i=1,nsites
               if (iaxdef(i,1).ne.0) then
                iz1 = iaxdef(i,1)
                ix = iaxdef(i,2)
                latyp = iaxdef(i,3)
               else
                iz1 = 0
                iz2 = 0
                ih1 = 0
                ih2 = 0
                latyp = 1
                do j=1,iconn(1,i)
                   jj = iconn(j+1,i)
                   if (jj.gt.0) then
                      if (ianz(jj).ne.1) then
                         if (iz1.eq.0) then
                            iz1 = jj
                         else
                            iz2 = jj
                         endif
                      else
                         if (ih1.eq.0) then
                            ih1 = jj
                         else
                            ih2 = jj
                         endif
                      endif
                   endif
                end do
                if (iz1.eq.0.and.ih1.ne.0) iz1 = ih1
                if (iz2.eq.0) then
                   if (ianz(i).eq.8.and.ih1.ne.0) then
                      ix = ih1
                      if (iz1.eq.ih1.and.ih2.ne.0) ix = ih2
                      latyp = 2
                   else
                      ix = 0
                      ih3 = 0
                      do k=1,iconn(1,iz1)
                         kk = iconn(k+1,iz1)
                         if (kk.gt.0.and.kk.ne.i) then
                            if (ianz(kk).ne.1) then
                               ix = kk
                            else
                               ih3 = kk
                            endif
                         endif
                      end do
                      if (ix.eq.0) then
                         if (iz1.ne.ih1.and.ih1.ne.0) then
                             ix = ih1
                         elseif (ih2.ne.0) then
                             ix = ih2
                         elseif (ih3.ne.0) then
                             ix = ih3
                         endif
                      endif
                   endif
                else
                   ix = iz2
                endif
               endif
               if (iz1.ne.0.and.ix.ne.0) then
                  call rotloc(i,iz1,ix,latyp,qrot)
                  call wrpol(idun,i,iz1,ix,latyp,qrot)
               else
                  ierr = 1
                  print*,'Atom ',i,' Error defining axes'
                  call inferr('Error defining local axes!',0)
               endif
             end do


         else
             ierr = 1
             call inferr('Multipoles # cites != # atoms!',0)
         endif
      endif

      if (iop.eq.16.or.iop.eq.17) then
         if (idun.ne.6) close(idun)
      endif
c
c put colors right again, we abused this array for dmarel writing
c
      iats = iatoms
      if (addcel) iats = iats - 8
      do i=1,iats
         mc = mcol(iabs(iatclr(i)))
         if (mc.gt.15) mc=15
         iatclr(i) = mc
      end do

      if (ifrst.ne.1) then
         do i=1,iatoms
            iaton(i) = 1
         end do
         if (iop.eq.5.or.iop.eq.6) then
            if (addcel) then
               itemp = natorg
               natorg = iatoms
            endif
            grdwt = grdw
            grdw = 0.9d0
            call allsrf(0,3,1,4)
c            call allsrf(0,0,1,1)
            grdw = grdwt
            if (addcel) natorg = itemp
         endif
         call docent
         call doscal
         scal = scali * 2.4d0 * smag
      endif
 
      if (iop.eq.18) then
         nspg = 1
         nopr = 1
      endif

      if (idebug.eq.1) then
         do i=1,iatoms
            print*,i,' ',elemnt(ianz(i)),(coo(j,i),j=1,3)
         end do
      endif

500   format(a)
      return
      end

      logical function applop(ir,it,a,b,icent,invrt)
      implicit double precision (a-h,o-z)
      integer*2 ir,it
      dimension ir(3,3),it(3),a(3),b(3)

      applop = .true.
      tol = 1.0d-5
      invrt = 1
      itrace = 0

      do i=1,3
         b(i) = dfloat(it(i))/12.0d0
         do j=1,3
            if (i.eq.j) itrace = itrace + (ir(i,j)-1)
            b(i) = b(i) + dfloat((ir(i,j)-1))*a(j)
         end do
         if (abs(b(i)).lt.tol) b(i) = 0.0d0
      end do

      if (itrace.eq.1.or.itrace.eq.-3) invrt = -1
  
      if (icent.eq.1) then
         do i=1,3
            b(i) = -b(i)
         end do
         invrt = -invrt
      endif

      dtot = dabs(a(1)-b(1))+dabs(a(2)-b(2))+dabs(a(3)-b(3))
      if (dtot.lt.tol) applop = .false.

      return
      end

      subroutine delad(iat,iacn,iorg,mxdma,coo,ianz,iatclr,qat,ityp)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      integer*2 ityp
      dimension iacn(*),iorg(*)
      dimension coo(3,*),ianz(*),iatclr(*),qat(*),ityp(*)

      do m=iat,iatoms-1

         do k=1,3
            coo(k,m) = coo(k,m+1)
         end do

         ianz(m) = ianz(m+1)
         iatclr(m) = iatclr(m+1)
         qat(m) = qat(m+1)
         ityp(m) = ityp(m+1)

         if (m.le.mxdma) then
            iacn(m) = iacn(m+1)
            iorg(m) = iorg(m+1)
         endif

      end do

      iatoms = iatoms - 1

      return
      end

      subroutine tvec(vec,vect,idebug)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      dimension vec(3),vv(3),vect(3),ivec(3)

c     find unit vector translation that end up in cell

      do j=1,3
         vv(j) = 0.0d0
      end do

      do i=0,10
         do j=1,3
            vt = vec(j)+dble(i)*1.0d0
            if (vt.ge.0.0d0.and.vt.lt.1.0d0) then
               vv(j) = vt
               ivec(j) = i
            endif
            vt = vec(j)-dble(i)*1.0d0
            if (vt.ge.0.0d0.and.vt.lt.1.0d0) then
               vv(j) = vt
               ivec(j) = -i
            endif
         end do
      end do

      do i=1,3
         vect(i) = vv(i) - vec(i)
      end do

      if (idebug.eq.1) then
          print*,'tvec ',(ivec(i),i=1,3)
      endif

      return
      end

      subroutine trcoo(tr,veco)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      dimension tr(3),veco(3)

c     translation

      do i=1,3
         veco(i) = veco(i) + tr(i)
      end do

      return
      end

      double precision function gsmod(f)
      implicit double precision (a-h,o-z)

      r1 = dint(f)
      if (f.lt.0.0d0) then
         r2 = dint(f)-1.0d0
      else
         r2 = dint(f)+1.0d0
      endif
 
      r = r1
      gsmod = dabs(f-r1)
      if (dabs(f-r2).lt.gsmod) r = r2

      gsmod = f-r

      return
      end

      subroutine fr2crt(vec,xa,ya,yb,za,zb,zc)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /frc/ ifr2crt
      dimension vec(3)

      if (ifr2crt.eq.0) then
         vec(3) = za*vec(1) + zb*vec(2) + zc*vec(3)
         vec(2) = ya*vec(1) + yb*vec(2)
         vec(1) = xa*vec(1)
      else
         vec(1) = xa*vec(1) + ya*vec(2) + za*vec(3)
         vec(2) =           + yb*vec(2) + zb*vec(3)
         vec(3) =                         zc*vec(3)
      endif

      return
      end

      subroutine crt2fr(veci,veco,xa,ya,yb,za,zb,zc)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /frc/ ifr2crt
      dimension veci(3),veco(3)

      if (ifr2crt.eq.0) then
         veco(1) =  veci(1) / xa
         veco(2) = (veci(2) - ya*veco(1)) / yb
         veco(3) = (veci(3) - za*veco(1) - zb*veco(2)) / zc 
      else
         veco(1) = (veci(1) - veci(2)*(ya/yb) + veci(3)*
     &                        ((zb/zc)*(ya/yb) - (za/zc))) / xa
         veco(2) =           (veci(2)         - veci(3)*(zb/zc)) / yb
         veco(3) =                              veci(3) / zc
      endif

      return
      end

      integer function itell(iatel,ianz)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      dimension iatel(*)

      if (ianz.eq.1) then
         indx = 1
      elseif (ianz.eq.6) then
         indx = 2
      elseif (ianz.eq.7) then
         indx = 3
      elseif (ianz.eq.8) then
         indx = 4
      elseif (ianz.eq.15) then
         indx = 5
      elseif (ianz.eq.16) then
         indx = 6
      elseif (ianz.eq.17) then
         indx = 7
      else
         indx = 8
      endif

      iatel(indx) = iatel(indx) + 1
      itell = iatel(indx)

      return
      end

      subroutine pold(nat,ipolh,iatom,ianz,iconn)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      dimension iconn(mxcon+1,*),ianz(*)

      ipolh = 0

      nstor = mxnat - nat
      do i=1,iconn(1,nstor+iatom)
         jatom = iconn(1+i,nstor+iatom)
         if (ianz(nstor+jatom).eq.7.or.ianz(nstor+jatom).eq.8) 
     &       ipolh = 1
      end do

      return
      end

      character*4 function dmalab(nat,ianz,indx,iacn,invrt,iytyp)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxel=100)
      character*2 elemnt
      common /elem/elemnt(mxel)
      character*2 ggstr,agstr
      dimension iytyp(*)

      if (ianz.eq.1) then
         dmalab(1:2) = 'HY'
         iytyp(1) = 1
         call polh(nat,ipolh,indx)
         if (ipolh.eq.1) then
            dmalab(1:2) = 'HP'
            iytyp(2) = 1
         endif
      elseif (ianz.eq.6) then
         dmalab(1:2) = 'CA'
         iytyp(3) = 1
      elseif (ianz.eq.7) then
         dmalab(1:2) = 'NI'
         iytyp(4) = 1
      elseif (ianz.eq.8) then
         dmalab(1:2) = 'OX'
         iytyp(5) = 1
      elseif (ianz.eq.15) then
         dmalab(1:2) = 'PH'
         iytyp(6) = 1
         print*,'WARNING Phosphor parameters are bogous'
      elseif (ianz.eq.16) then
         dmalab(1:2) = 'SU'
         iytyp(7) = 1
      elseif (ianz.eq.17) then
         dmalab(1:2) = 'CL'
         iytyp(8) = 1
      elseif (ianz.eq.11) then
         dmalab(1:2) = 'NA'
      else
         dmalab(1:2) = elemnt(ianz)
      endif
      lind = iacn
      if (ianz.eq.99) then
         if (invrt.lt.0) then
            dmalab = dmalab(1:2)//agstr(lind,.true.)
         else
            dmalab = dmalab(1:2)//agstr(lind,.false.)
         endif
      else
         if (invrt.lt.0) lind = lind + 50
         dmalab = dmalab(1:2)//ggstr(lind)
      endif

      return
      end

      
      subroutine getdis(vdwr1,vdwr2,ianz1,ianz2,dist)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxel=100)
      common /exmass/ emass(mxel),cutt(4,4)

      ind1 = 0
      ind2 = 0
      if (ianz1.eq.1) then
         ind1 = 1
      elseif (ianz1.ge.6.and.ianz1.le.8) then
         ind1 = ianz1 - 4
      endif
      if (ianz2.eq.1) then
         ind2 = 1
      elseif (ianz2.ge.6.and.ianz2.le.8) then
         ind2 = ianz2 - 4
      endif

      if (ind1.ne.0.and.ind2.ne.0) then
         dist = cutt(ind1,ind2)
      else
         dist = (vdwr1 + vdwr2) * 0.9d0
      endif

      return
      end

      subroutine ifnd(ifnn,noff,ia1,ia2,iconn)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      dimension iconn(mxcon+1,*)


      ifnn = 0

      do i=1,iconn(1,noff+ia1)
         in = iconn(i+1,noff+ia1)
         if (in.eq.ia2) then
            ifnn = 1
            return
         endif
         do j=1,iconn(1,noff+in)
            jn = iconn(j+1,noff+in)
            if (jn.ne.ia1) then
               if (jn.eq.ia2) then
                  ifnn = 2
                  return
               endif
               do k=1,iconn(1,noff+jn)
                  kn = iconn(k+1,noff+jn)
                  if (kn.ne.in) then
                     if (kn.eq.ia2) then
                        ifnn = 3
                        return
                     endif
                     do l=1,iconn(1,noff+kn)
                        ln = iconn(l+1,noff+kn)
                        if (ln.ne.jn) then
                           if (ln.eq.ia2) then
                              ifnn = 4
                              return
                           endif
                           do m=1,iconn(1,noff+ln)
                              mn = iconn(l+1,noff+ln)
                              if (mn.ne.kn) then
                                 if (mn.eq.ia2) then
                                    ifnn = 5
                                    return
                                 endif
                                 do n=1,iconn(1,noff+mn)
                                    nn = iconn(l+1,noff+mn)
                                    if (nn.ne.ln) then
                                       if (nn.eq.ia2) then
                                          ifnn = 6
                                          return
                                       endif
                                    endif
                                 end do
                              endif
                           end do
                        endif
                     end do
                  endif
               end do
            endif
         end do
      end do

      return
      end

      subroutine alntwo(coo,i3,car,nsites,j3)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      dimension coo(3,*),car(3,*),v1(3),v2(3),v3(3),vct1(3),vct2(3)
      dimension i3(3),j3(3)

      tol = 1.0d-6

      do j=1,3
         v1(j) = coo(j,i3(2)) - coo(j,i3(1))
         vct1(j) = car(j,j3(2)) - car(j,j3(1))
      end do

      call vsc1(v1,1.0d0,tol)
      call vsc1(vct1,1.0d0,tol)
      call crprod(v1,vct1,v2)

      if (vlen(v2).gt.tol) then
         call vsc1(v2,1.0d0,tol)
         call crprod(v1,v2,v3)
         call vsc1(v3,1.0d0,tol)
         call impsc(v1,vct1,cosa)
         call impsc(v3,vct1,sina)

         do i=1,nsites
            if (i.ne.j3(1)) then

               do j=1,3
                  vct1(j) = car(j,i) - car(j,j3(1))
               end do

               call timpsc(vct1,v1,vct2(1))
               call timpsc(vct1,v2,vct2(2))
               call timpsc(vct1,v3,vct2(3))

               do j=1,3
                  car(j,i) = ( cosa*vct2(1) + sina*vct2(3))*v1(j) +
     &                       (-sina*vct2(1) + cosa*vct2(3))*v3(j) +
     &                        v2(j)*vct2(2) + car(j,j3(1))
               end do

            endif

         end do

      endif
c
      do j=1,3
         v1(j) = coo(j,i3(2)) - coo(j,i3(1))
         v2(j) = coo(j,i3(3)) - coo(j,i3(1))
      end do

      call crprod(v1,v2,vct1)

      do j=1,3
         v1(j) = car(j,j3(2)) - car(j,j3(1))
         v3(j) = car(j,j3(3)) - car(j,j3(1))
      end do

      call crprod(v1,v3,v2)

      if (vlen(vct1).gt.tol.and.vlen(v2).gt.tol) then

         call vsc1(vct1,1.0d0,tol)
         call vsc1(v1,1.0d0,tol)
         call vsc1(v2,1.0d0,tol)
         call crprod(v1,v2,v3)

         if (vlen(v3).gt.tol) then

            call vsc1(v3,1.0d0,tol)
            call impsc(vct1,v2,cosa)
            call impsc(vct1,v3,sina)

            sina = -sina

            do i=1,nsites
               if (i.ne.j3(1).and.i.ne.j3(2)) then

                  do j=1,3
                     vct1(j) = car(j,i) - car(j,j3(1))
                  end do

                  call timpsc(vct1,v1,vct2(1))
                  call timpsc(vct1,v2,vct2(2))
                  call timpsc(vct1,v3,vct2(3))

                  do j=1,3
                     car(j,i) = ( cosa*vct2(2) + sina*vct2(3))*v2(j) +
     &                          (-sina*vct2(2) + cosa*vct2(3))*v3(j) +
     &                           v1(j)*vct2(1) + car(j,j3(1))
                  end do

               endif

            end do

         endif

      endif

      return
      end

      subroutine trtwo(coo,ianz1,nsum,car,ianz2,nsites,allsit)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      logical allsit
      dimension coo(3,*),car(3,*),ianz1(*),ianz2(*),vct1(3),vct2(3)

c
c center of mass
c
      call cntvec(vct1,coo,ianz1,nsum)

      if (allsit) then
         call cntvec(vct2,car,ianz2,nsites)
      else
         call cntvec(vct2,car,ianz2,nsum)
      endif

      do i=1,nsites
         do j=1,3
            car(j,i) = car(j,i) - vct2(j) + vct1(j)
         end do
      end do

      return
      end

      subroutine wrbuck(idun,lab1,lab2,alk,blk,clk)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      character*2 lab1,lab2

      write(idun,'(a4,2x,2(a2,4x,a4,4x))') 
     &      'BUCK',lab1,'CODA',lab2,'CODA'
      write(idun,'(1x,5(f14.8,2x))') 
     &      alk,blk,clk,0.0d0,70.0d0
      write(idun,'(a)') 'ENDS'

      return
      end

      subroutine chkbrd(iconn,icalf,ianf,islu,iamino,isal,reson,
     &                  ncalf,nchain)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxchai=50)
      integer reson
      dimension iconn(mxcon+1,*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),reson(*),isal(*)

      nchain = 1
      ianf(1) = 1

      do i=1,ncalf-1
         inewch = 0
         if (iamino(i).le.23) then
            isami = 1
         else
            isami = 0
         endif

         if (iamino(i+1).le.23) then
            isnew = 1
         else
            isnew = 0
         endif

         if (isami.ne.isnew) then
            inewch = 1
         else
            ifnd = 0
            if (isami.eq.1) then
               ii = icalf(3,i)
               jj = icalf(2,i+1)
               if (ii.gt.0) then
                  do j=1,iconn(1,ii)
                     if (iconn(1+j,ii).eq.jj) ifnd = 1
                  end do
               endif
            else
               ii = icalf(6,i)                   
               jj = icalf(1,i+1)
               if (ii.gt.0) then
                  do j=1,iconn(1,ii)
                     if (iconn(1+j,ii).eq.jj) ifnd = 1
                  end do
               endif
            endif
            if (ifnd.eq.0) inewch = 1
         endif

         if (inewch.eq.1.and.nchain.lt.mxchai) then
            if (nchain.gt.0) islu(nchain) = i
            nchain = nchain + 1
            ianf(nchain) = i + 1
         endif

         reson(i) = 1
         isal(i) = 3
      end do

      reson(ncalf) = 1
      islu(nchain) = ncalf

      return
      end 

      subroutine rdgrod(idebug,istat,
     &                  coo,qat,ianz,iatclr,iresid,iconn,ityp,ipdbt,
     &                  icalf,ncalf,ianf,islu,nchain,iamino,reson,
     &                  isal,irsnr,ishoh)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxel=100)
      parameter (mxop=2822)
      parameter (mxrop=60)
      parameter (mxtop=56)
      parameter (mxsg=238)
      parameter (mxspg=232)
      parameter (mxppmf=16)
      parameter (mxlpmf=26)
      parameter (mxmol2=41)
      parameter (mxmm3=164)
      parameter (mxchtp=136)
      parameter (mxmsf=235)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxamo=201)
      parameter (mxchai=50)
      parameter (mxres=42)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxmmul=100)

      common /athlp/ iatoms, mxnat
      character*137 line
      character*2 elemnt
      common /elem/elemnt(mxel)
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
      integer*2 ityp,ipdbt
      common /types/  iff
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /settng/ jset(3,9),iset(6,mxspg),nset
      integer reson
      character*3 pdbsym,hsym,chtnk,ambtnk
      character*2 amotnk,gffstr
      common /symbol/ pdbsym(mxsym),hsym(mxhsym),chtnk(mxchtp),
     &                ambtnk(mxamb),amotnk(mxamo),gffstr(mxgff)
      character*3 aminos
      common /amino/  aminos(mxres)
      common /cllab/  iclon,iclpnt(4)
      common /surf/   natorg,noscnd

      common /multim/ imulm, nmulm,ihasqm(mxmmul)
      character*137 ftopf
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character fniun*256
      common /fnunit/ fniun
      common /curlin/ line
      integer getlin
      character*137 string
      character*4 tstr
      character*5 resnm,atnm
      logical opfil, boxadd, dall
      dimension ihpdb(mxhsym*3)
      dimension v1(3),v2(3),v3(3)
      dimension coo(3,*),ianz(*),iatclr(*),iresid(*),iconn(mxcon+1,*),
     &          ityp(*),ipdbt(*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),reson(*),isal(*),
     &          irsnr(*),qat(*)

      dall = .false.
      toang = 0.52917706d0
      istat = 1
      nhmol = 3
      imulm = 1
      nmulm = 0

      call rewfil

      if (idebug.eq.1) print*,'rdgro'

      call redel(line,1)

      call inferr('.GRO file',0)

      if (getlin(0).eq.1) then
          ktype = nxtwrd(string,nstr,itype,rtype)
          if (ktype.eq.2) then
              iatoms = itype
              if (iatoms.gt.mxnat) then
                 dall = .true.
                 goto 100
              endif
          else
              goto 100
          endif
      else
          goto 100
      endif

      call parsfn('Helix',5,1)
      call parsfn('Beta',4,1)
      call parsfn('RNA/DNA',7,1)
      call parsfn('Coil',4,1)

      ncalf = 0
      ishoh = 0
      ihohl = 0
      ihohu = 0
      irstmp = 0
      irsold = 0

      do i=1,iatoms
         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 100

         read(line,'(i5,2a5,i5,3f8.3,3f8.4)',err=100) iresid(i),
     &      resnm,atnm,idum,(coo(j,i),j=1,3)

c nanometer -> angstrom -> bohr

         do j=1,3
            coo(j,i) = coo(j,i)*10.0d0
         end do

         ianz(i) = 0
         iconn(1,i) = 0
         iatclr(i) = 1

c new residue
         if (iresid(i).ne.irsold) then
            irsold = iresid(i)
            jatyp = 0

            do j=1,mxres
               if (resnm(1:3).eq.aminos(j)) jatyp = j
            end do

            if (jatyp.ne.0) then
c amino acid
               ncalf = ncalf + 1
               iamino(ncalf) = jatyp
               irsnr(ncalf) = iresid(i)
            else
c other residue
               nhmol = nhmol + 1
               if (nmulm.lt.mxmmul) then
                  nmulm = nmulm + 1
                  ihasqm(nmulm) = 0
               endif

               if (resnm.eq.'WATER'.or.resnm.eq.'HOH'.or.
     &             resnm.eq.'SOL') then
                  if (ishoh.eq.0) then
                     ishoh = nhmol
                  endif
                  if (nhmol.gt.ihohu) ihohu = nhmol
                  iresid(i) = -nhmol
               else
                  call parsfn(resnm,5,1)
                  iresid(i) = -nhmol
               endif
            endif
            irstmp = iresid(i)

         else

c next atom of the same residue
            iresid(i) = irstmp
         endif

c process atom label
c gromacs string is right justified, we need left justified

         tstr(1:4) = atnm(2:5) 
         if (tstr(2:3).eq.'  ') then
            tstr(2:4) = tstr(4:4)//'  '
         elseif (tstr(2:2).eq.' ') then
            tstr(2:4) = tstr(3:4)//' '
         endif
         if (resnm(1:3).eq.'HEM') then
            if (tstr(2:3).eq.'FE') tstr(1:3) = 'FE '
         endif

         call detanz(resnm(1:3),tstr,line,ifnd,ish,ianz(i),ihashy)

         tstr(1:4) = tstr(2:4)//' ' 
         if (tstr(1:3).eq.'O1 ') tstr(1:3) = 'O  '
         if (tstr(1:3).eq.'O2 ') tstr(1:3) = 'OXT'
         if (iamino(ncalf).le.20) then
             if (tstr(1:3).eq.'H1 ') tstr(1:3) = 'H  '
             if (tstr(1:3).eq.'H2 ') tstr(1:3) = 'H  '
             if (tstr(1:3).eq.'H3 ') tstr(1:3) = 'H  '
         endif

         ipdbt(i) = 0
         do j=1,mxsym
            if (tstr(1:3).eq.pdbsym(j)) ipdbt(i) = j
         end do

         if (ipdbt(i).eq.0) then
c  check H atomname
            if (tstr(1:1).eq.'H') then
               do j=1,mxhsym
                  if (tstr(1:3).eq.hsym(j)) then
                     ipdbt(i) = (j-1)*3 + 1
                  endif
               end do
               if (ipdbt(i).eq.1)  icalf(4,ncalf) = i
            endif
         else
c known non H atomname
            if (ipdbt(i).eq.1)  icalf(2,ncalf) = i
            if (ipdbt(i).eq.2)  icalf(1,ncalf) = i
            if (ipdbt(i).eq.3)  icalf(3,ncalf) = i

            if (ipdbt(i).eq.43) icalf(1,ncalf) = i
            if (ipdbt(i).eq.46) icalf(2,ncalf) = i
            if (ipdbt(i).eq.47) icalf(3,ncalf) = i
            if (ipdbt(i).eq.48) icalf(4,ncalf) = i
            if (ipdbt(i).eq.50) icalf(5,ncalf) = i
            if (ipdbt(i).eq.51) icalf(6,ncalf) = i

         endif

      end do

c read box parameters
      do i=1,3
         v1(i) = 0.0d0
         v2(i) = 0.0d0
         v3(i) = 0.0d0
      end do
      v1(1) = 1.0d0
      v2(2) = 1.0d0
      v3(3) = 1.0d0

      boxadd = .false.

      call nxtlin(line,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 60
      read(line,*,end=55,err=55) v1(1),v2(2),v3(3),
     &                    v1(2),v1(3),v2(1),v2(3),v3(1),v3(2)
55    read(line,*,end=60,err=60) v1(1),v2(2),v3(3)
c
c succesfully read in box parameters
c
      boxadd = .true.

60    continue

      if (ishoh.ne.0) then
         call parsfn('Solvent',7,1)
      endif


c adjust H types

      do j=1,ncalf
         do i=1,mxhsym*3
             ihpdb(i) = 0
         end do

         do i=1,iatoms
            if (iresid(i).eq.j) then
               if (ianz(i).eq.1) then
                  if (ihpdb(ipdbt(i)).eq.0) then
                     ihpdb(ipdbt(i)) = i
                  elseif (ihpdb(ipdbt(i)+1).eq.0) then
                     ipdbt(i) = ipdbt(i) + 1
                     ihpdb(ipdbt(i)) = i
                  elseif (ihpdb(ipdbt(i)+2).eq.0) then
                     ipdbt(i) = ipdbt(i) + 2
                     ihpdb(ipdbt(i)) = i
                  endif
               endif
            endif
         end do

         ia = iamino(j)

         do i=1,mxhsym*3
            if (ihpdb(i).ne.0) then

               itmp = i
               iatmp = ihpdb(i)

               if (ia.ge.17.and.ia.le.20) then
c HIS PHE TYR TRP
c 1HD,1HE
                  if (i.eq.19.or.i.eq.28) itmp = i + 3
c 2HD,2HE
                  if (i.eq.20.or.i.eq.29) itmp = i + 5
c TRP: 2HE -> HE3
                  if (ia.eq.20.and.i.eq.29) itmp = 37
                           
               endif

               if (ia.eq.20) then
c TRP: 1HZ -> HZ2, 2HZ -> HZ3
                  if (i.eq.40) itmp = 46
                  if (i.eq.41) itmp = 49
               endif

c THR: 1HG -> HG1
               if (ia.eq.5.and.i.eq.10) itmp = 13

               if (itmp.ne.i) then
                   ihpdb(i) = 0
                   ihpdb(itmp) = iatmp
                   ipdbt(iatmp) = itmp
               endif

            endif
         end do
               
      end do

c Before we can go on we need to read connectivity or create it 
 
10    continue
      
      iuncon = 0
      i1 = icdex(fniun,'.gro')
      call parsfn(fniun,linlen(fniun)-4,18)
      if (i1.ne.0) then
         ftopf = fniun(1:i1-1)//'.top'
         nlen = i1-1+4
         if (opfil(46,ftopf,nlen,1,1,1)) then
            iuncon = 46
         else 
            print*,' '
            print*,'Didnt find topology file: '//ftopf
            print*,'Molden will generate connectivity'
            print*,' '
         endif
      endif

      if (iuncon.ne.0) then
         do i=1,nmulm
            ihasqm(1) = 1
         end do
c if there is a .top file get connectivity from there and charges
         call gettop(qat,iconn,iuncon,ityp)
         call conslv
      else
         call convpdb
      endif

      do i=1,iatoms
         irs = iresid(i)
         if (irs.le.-ishoh.and.irs.ge.-ihohu) then
            iresid(i) = -ishoh
         endif
      end do
c
c process box parameters
c
      if (boxadd) call addbox(v1,v2,v3)
      natorg = iatoms

c check for chain breaks

      call chkbrd(iconn,icalf,ianf,islu,iamino,isal,reson,
     &                  ncalf,nchain)

c     get trajectory file

      nstep = 0
      call gettrr(nstep,istats)
      if (istats.eq.0) then
         call getene(nstep,istats)
         close(iun2)
      endif
         
      return

100   istat = 0
      iatoms = 0
      iff = 0
      if (dall) then
         istat = -1
         call rewfil
      endif

      return
      end

      subroutine addbod(v1,v2,v3,coo,ianz,iatclr,iconn)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /cllab/ iclon,iclpnt(4)
      common /cllhlp/ ico(3,8),icn(4,8),mcol(32)
      common /athlp/ iatoms, mxnat
      dimension v1(*),v2(*),v3(*)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*)

      iclon = 1
      iclpnt(1) = iatoms + 1
      do j=1,3
         coo(j,iatoms+1) = 0.0d0
      end do

      iclpnt(2) = iatoms + 2
      do j=1,3
         coo(j,iatoms+2) = v1(j)*10.0d0
      end do

      iclpnt(3) = iatoms + 3

      do j=1,3
         coo(j,iatoms+3) = v2(j)*10.0d0
      end do

      iclpnt(4) = iatoms + 4
      do j=1,3
         coo(j,iatoms+4) = v3(j)*10.0d0
      end do

      do j=1,3
         coo(j,iatoms+5) = (v1(j)+v2(j))*10.0d0
      end do

      do j=1,3
         coo(j,iatoms+6) = (v2(j)+v3(j))*10.0d0
      end do

      do j=1,3
         coo(j,iatoms+7) = (v1(j)+v3(j))*10.0d0
      end do

      do j=1,3
         coo(j,iatoms+8) = (v1(j)+v2(j)+v3(j))*10.0d0
      end do

      do i=1,8
         l = iatoms + i
         iconn(1,l) = icn(1,i)
         do j=2,4
            iconn(j,l) = iatoms + icn(j,i)
         end do
         ianz(l) = 100
         iatclr(l) = 11
      end do
      iatoms = iatoms + 8

      return
      end

      subroutine addtbd(v1,v2,v3,coo,ianz,iatclr,iconn,iaton,iresid)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /cllab/ iclon,iclpnt(4)
      common /cllhlp/ ico(3,8),icn(4,8),mcol(32)
      common /athlp/ iatoms, mxnat
      dimension v1(*),v2(*),v3(*)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*),iaton(*),
     &          iresid(*)

      iclpnt(1) = iatoms + 1
      do j=1,3
         coo(j,iatoms+1) = -0.5d0*(v1(j)+v2(j)+v3(j))
      end do

      iclpnt(2) = iatoms + 2
      do j=1,3
         coo(j,iatoms+2) = -0.5d0*(-v1(j)+v2(j)+v3(j))
      end do

      iclpnt(3) = iatoms + 3

      do j=1,3
         coo(j,iatoms+3) = -0.5d0*(v1(j)-v2(j)+v3(j))
      end do

      iclpnt(4) = iatoms + 4
      do j=1,3
         coo(j,iatoms+4) = -0.5d0*(v1(j)+v2(j)-v3(j))
      end do

      do j=1,3
         coo(j,iatoms+5) = -0.5d0*(-v1(j)-v2(j)+v3(j))
      end do

      do j=1,3
         coo(j,iatoms+6) = -0.5d0*(v1(j)-v2(j)-v3(j))
      end do

      do j=1,3
         coo(j,iatoms+7) = -0.5d0*(-v1(j)+v2(j)-v3(j))
      end do

      do j=1,3
         coo(j,iatoms+8) = 0.5d0*(v1(j)+v2(j)+v3(j))
      end do

      do i=1,8
         l = iatoms + i
         iconn(1,l) = icn(1,i)
         do j=2,4
            iconn(j,l) = iatoms + icn(j,i)
         end do
         ianz(l) = 100
         iatclr(l) = 11
         iaton(l) = 1
         iresid(l) = -4
      end do
      iatoms = iatoms + 8

      return
      end

      subroutine gettop(qat,iconn,iuncon,ityp)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxgmx=53)
      parameter (mxgmx2=57)
      parameter (mxg43=49)
      parameter (mxcon=10)
      character*137 line
      character*5 styp
      integer*2 ityp
      common /types/ iff
      common /athlp/ iatoms, mxnat
      character*5 gro43,grogmx,grog2x,gro43l,grogmxl,grog2xl
      character*35 gro43s, grogms, grog2s
      common /symgro/ gro43(mxg43),grogmx(mxgmx),grog2x(mxgmx2),
     &                gro43l(mxg43),grogmxl(mxgmx),grog2xl(mxgmx2),
     &                gro43s(mxg43),grogms(mxgmx),grog2s(mxgmx2)
      dimension qat(*), iconn(mxcon+1,*),ityp(*)

      rewind iuncon
      do while (.true.)
         read(iuncon,'(a)',err=80,end=80) line
         if (icdex(line,';').eq.0) then
            if (icdex(line,'[').ne.0) then
               if (icdex(line,'atoms').ne.0) then
                  ici = 0
                  iat = 0
                  do while (.true.)
                     read(iuncon,'(a)',err=80,end=80) line
                     if (icdex(line,'[').ne.0) then
                         backspace iuncon
                         goto 70
                     endif
                     if (ici.eq.0) then
                        ici = icdex(line,'charge')
                        if (ici.ne.0) then
                           read(iuncon,'(a)',err=80,end=80) line
                        endif
                     endif
                     
                     ic1 = icdex(line,';')
                     if (ic1.ne.0) line = line(1:ic1-1)
                     
                     if (ici.ne.0) then
                         styp = line(13:17)
                         ll = len(line)
                         line = line(ici:ll)
                         read(line,'(f6.3)',err=80) q
                         iat = iat + 1
                         qat(iat) = q
                         ityp(iat) = 0
                         if (iff.eq.9) then
                            do l=1,mxgmx
                               if (styp.eq.grogmx(l)) ityp(iat) = l
                            end do
                         elseif (iff.eq.10) then
                            do l=1,mxgmx2
                               if (styp.eq.grog2x(l)) ityp(iat) = l
                            end do
                         elseif (iff.eq.11) then
                            do l=1,mxg43
                               if (styp.eq.gro43(l)) ityp(iat) = l
                            end do
                         endif
                     endif
                  end do
70                continue
               elseif (icdex(line,'bonds').ne.0) then
                  do while (.true.)
                     read(iuncon,'(a)',err=80,end=80) line
                     if (icdex(line,'[').ne.0) goto 80
                     if (icdex(line,';').eq.0) then
                         read(line,'(i5,x,i5)',err=80) i1,i2
                         if ((i1.gt.0.and.i1.lt.mxnat).and.
     &                       (i2.gt.0.and.i2.lt.mxnat)) then
                            if (iconn(1,i1).lt.mxcon) then
                               iconn(2+iconn(1,i1),i1) = i2
                               iconn(1,i1) = iconn(1,i1) + 1
                            endif
                            if (iconn(1,i2).lt.mxcon) then
                               iconn(2+iconn(1,i2),i2) = i1
                               iconn(1,i2) = iconn(1,i2) + 1
                            endif
                         endif
                     endif
                  end do
               endif
            elseif (index(line,'#include '//char(34)//'ff').ne.0)
     &         then
               if (icdex(line,'ffG43').ne.0) then
                   nfftyp = mxg43
                   iff = 11
               endif
               if (icdex(line,'ffgmx2').ne.0) then
                   iff = 10
                   nfftyp = mxgmx2
               elseif (icdex(line,'ffgmx').ne.0) then
                   iff = 9
                   nfftyp = mxgmx
               endif
            endif
         endif
      end do

80    continue

      close(iuncon)

      return
      end

      subroutine gettrra(fniun,iuncon,nstep,istat)
      implicit double precision (a-h,o-z)
      parameter (MAXPNT=2000)
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      common /geocnv2/ formax(MAXPNT),forrms(MAXPNT),dismax(MAXPNT),
     &                 disrms(MAXPNT),epoints(MAXPNT),isav(MAXPNT)
      character*(*) fniun
      character*137 ftopf,line
      logical opfil,extest

      istat = 0
      igcvav = 1
      ifmxav = 0
      ifrmav = 0
      idmxav = 0
      idrmav = 0
      ieav = 0
      ngeoms = 0
      nepnts = 0
      dmaxt = 0.0d0
      fmaxt = 0.0d0
      drmst = 0.0d0
      frmst = 0.0d0
      istp1 = 0
      istp2 = 0

      i1 = icdex(fniun,'.gro')
      if (i1.ne.0) then
         ftopf = fniun(1:i1-1)//'.trra'
         inquire(file=ftopf,exist=extest)
         if (.not.extest) then
            ftopf = fniun(1:i1-1)//'.trr'
            inquire(file=ftopf,exist=extest)
            if (extest) then
               ftopf = 
     &         'gmxdump -f '//fniun(1:i1-1)//'.trr > '//
     &                        fniun(1:i1-1)//'.trra'
               nlen = 2*(i1-1)+23
               iextmp = 0
               print*,
     &         'Running gmxdump to create ascii trajectory file.'
               print*,' '
               call exstr(ftopf,nlen,iextmp)
               ftopf = fniun(1:i1-1)//'.trra'
            else
               print*,' '
               print*,'Didnt find binary trajectory file: '//ftopf
               print*,'Didnt find ascii trajectory file: '
     &                 //ftopf(1:i1+3)//'a'
               print*,' '
               goto 100
            endif
         endif
 
         if (opfil(iuncon,ftopf,i1+4,1,1,0)) then
            rewind iuncon
            istat = 1
            ngeoms = 0 
            do while (.true.)
               read(iuncon,'(a)',err=100,end=100) line
               if (icdex(line,'frame').ne.0) then
                  ngeoms = ngeoms + 1
                  if (ngeoms.eq.1) then
                     read(iuncon,'(a)',err=100,end=100) line
                     it1 = icdex(line,'step=')
                     it2 = icdex(line,'time=')
                     read(line(it1+5:it2-1),*) istp1
                  elseif (ngeoms.eq.2) then
                     read(iuncon,'(a)',err=100,end=100) line
                     it1 = icdex(line,'step=')
                     it2 = icdex(line,'time=')
                     read(line(it1+5:it2-1),*) istp2
                  endif
               endif
            end do
         else
            print*,' '
            print*,
     &        'Didnt find ascii trajectory file: '//ftopf
            print*,' '
         endif

      endif
   
100   continue
      nstep = istp2 - istp1 
      return
      end

      subroutine getenea(fniun,iuncon,nstep,nestp,istat)
      implicit double precision (a-h,o-z)
      parameter (MAXPNT=2000)
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      common /geocnv2/ formax(MAXPNT),forrms(MAXPNT),dismax(MAXPNT),
     &                 disrms(MAXPNT),epoints(MAXPNT),isav(MAXPNT)
      character*(*) fniun
      character*137 ftopf,line
      logical opfil,extest

      istat = 0
      ieav = 0
      nepnts = 0
      istp1 = 0
      istp2 = 0
      nestp = 0
      eprepos = 1.54321d-20

      i1 = icdex(fniun,'.gro')
      if (i1.ne.0) then
         ftopf = fniun(1:i1-1)//'.edra'
         inquire(file=ftopf,exist=extest)
         if (.not.extest) then
            ftopf = fniun(1:i1-1)//'.edr'
            inquire(file=ftopf,exist=extest)
            if (extest) then
               ftopf = 
     &         'gmxdump -e '//fniun(1:i1-1)//'.edr > '//
     &                        fniun(1:i1-1)//'.edra'
               nlen = 2*(i1-1)+23
               iextmp = 0
               print*,
     &         'Running gmxdump to create ascii energy file.'
               print*,' '
               call exstr(ftopf,nlen,iextmp)
               ftopf = fniun(1:i1-1)//'.edra'
            else
               print*,' '
               print*,'Didnt find binary energy file: '//ftopf
               print*,'Didnt find ascii energy file: '//
     &                 ftopf(1:i1+3)//'a'
               print*,' '
               goto 200
            endif
         endif
 
         if (opfil(iuncon,ftopf,i1+4,1,1,0)) then

c finds step in energy points
            rewind iuncon
            nepnts = 0 
            do while (.true.)
               read(iuncon,'(a)',err=100,end=100) line
               ind1 = index(line,'step:')
               if (ind1.gt.8) then
                  nepnts = nepnts + 1
                  if (nepnts.eq.1) then
                     read(line(ind1+5:),*) istp1
                  elseif (nepnts.eq.2) then
                     read(line(ind1+5:),*) istp2
                  else
                     goto 50
                  endif
               endif
            end do

50          rewind iuncon
            nestp = istp2 - istp1
            rfrac = dble(nestp) / dble(nstep)
            do i=1,MAXPNT
               epoints(i) = eprepos
            end do
            istat = 1
            nepnts = 0 
            l = 0
            do while (.true.)
               read(iuncon,'(a)',err=100,end=100) line
               ind1 = index(line,'Total Energy')
               if (ind1.gt.8) then
                  line = line(ind1+12:ind1+26)
                  l = l + 1
                  nepnts = int(dble(l-1)*rfrac) + 1
                  if (nepnts.le.MAXPNT) then
                     read(line,*) epoints(nepnts)
c               print*,'epoints=',epoints(nepnts),' nepnts=',nepnts
                  endif
                  ieav = 1
               endif
            end do
         else
            print*,' '
            print*,
     &        'Didnt find ascii energy file: '//ftopf
            print*,' '
            goto 200
         endif

      endif
   
100   continue
      if (ieav.eq.1) then
         iold = 1
         do i=1,nepnts
            if (epoints(i).eq.eprepos) then
               epoints(i) = epoints(iold)
            else
               iold = i
            endif
         end do
c         call gmmcnv
         close(iuncon)
      endif

200   return
      end

      subroutine gropd(ipnt,coo)
      implicit double precision (a-h,o-z)
      character*137 line
      logical getd3
      common /athlp/ iatoms, mxnat
      common /cllab/ iclon,iclpnt(4)
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      dimension d(3),v(3,3),coo(3,*)

c      print*,'gropt ',ipnt
      istat = 0
      toang = 0.52917706d0

      do while (.true.)
         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 100
         read(line,'(a)',err=100,end=100) line
         ind1 = icdex(line,'frame')
         if (ind1.ne.0) then
            ind2 = index(line,':') 
            if (ind2.ne.0) then
               read(line(ind1+6:ind2-1),*) iframe
c               print*,'iframe=',iframe,' ipnt=',ipnt
               if (iframe.eq.ipnt) then
                  iat = 0
c read updated box
                  if (iclon.eq.1) then
                     ioff = iclpnt(1)
                     call redel(line,2)
                     do i=1,3
                        call nxtlin(line,jstat)
                        if (jstat.eq.1.or.jstat.eq.2) goto 100
                        if (getd3(d,line)) then
                           do j=1,3
                              v(j,i) = d(j)*10.0d0/toang
                           end do
                        else
                           goto 100
                        endif
                     end do

                     do i=1,3
                        do j=1,3
                           coo(j,ioff+i) = v(j,i)
                        end do
                     end do

                     do j=1,3
                        coo(j,ioff+4) = (v(j,1)+v(j,2))
                     end do

                     do j=1,3
                        coo(j,ioff+5) = (v(j,2)+v(j,3))
                     end do

                     do j=1,3
                        coo(j,ioff+6) = (v(j,1)+v(j,3))
                     end do

                     do j=1,3
                        coo(j,ioff+7) = (v(j,1)+v(j,2)+v(j,3))
                     end do

                     call redel(line,1)
                  else
                     call redel(line,6)
                  endif

                  do while (.true.)
                     call nxtlin(line,jstat)
                     if (jstat.eq.1.or.jstat.eq.2) goto 100
                     if (getd3(d,line)) then
                        iat = iat + 1
                        do j=1,3
                           coo(j,iat) = d(j)*10.0d0/toang
                        end do
                     else
                        goto 100
                     endif
                  end do
               endif
            endif
         endif
      end do

   
100   continue
      call bckfil

      return
      end

      logical function getd3(d,line)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      character*(*) line
      dimension d(3)

      getd3 = .true.
      i1 = index(line,'{')
      i2 = index(line,'}')
      if (i1.gt.0.and.i2.gt.0) then
         line = line(i1+1:i2-1)
         read(line,*,err=100) (d(j),j=1,3)
      else
         goto 100
      endif

      return

100   getd3 = .false.
      return
      end

      subroutine rdmod(idebug,ipdbon,ioadd,istat,
     &                 coo,qat,ianz,iatclr,iresid,iconn,ityp,ipdbt,
     &                 icalf,ncalf,ianf,islu,nchain,iamino,reson,
     &                 isal,ishoh,ihashb,
     &                 nat,norg,icent,ncon,nspg,kz,ichx,nopr,ir,it,
     &                 xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxel=100)
      parameter (mxspg=232)
      parameter (mxppmf=16)
      parameter (mxlpmf=26)
      parameter (mxmol2=41)
      parameter (mxmm3=164)
      parameter (mxchtp=136)
      parameter (mxmsf=235)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxamo=201)
      parameter (mxchai=50)
      parameter (mxres=42)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxmmul=100)

      common /athlp/ iatoms, mxnat
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      character*137 line
      character*2 elemnt,tolowf,atom,atomt
      common /elem/  elemnt(mxel)
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
      integer*2 ityp,ipdbt
      common /types/  iff
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      integer*2 ir,it
      common /settng/ jset(3,9),iset(6,mxspg),nset
      integer reson
      character*3 pdbsym,hsym,chtnk,ambtnk
      character*2 amotnk,gffstr
      common /symbol/ pdbsym(mxsym),hsym(mxhsym),chtnk(mxchtp),
     &                ambtnk(mxamb),amotnk(mxamo),gffstr(mxgff)
      character*3 aminos
      common /amino/aminos(mxres)
      common /multim/ imulm, nmulm,ihasqm(mxmmul)

      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /curlin/ line
      integer getlin
      character*137 string
      integer*2 irr,itt
      character*3 tstr
      character*2 gstr
      character*40 rstr
      logical trypdb,dall
      dimension icont(mxcon+1)
      dimension ihpdb(mxhsym*3)
      dimension coo(3,*),ianz(*),iatclr(*),iresid(*),iconn(mxcon+1,*),
     &          qat(*),ityp(*),ipdbt(*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),reson(*),isal(*)
      dimension ir(3,3,192),it(3,192)


c Hydrogen symbols used by MOL2 in relation with Molden:
c
c Besides the H, a hydrogen can have three additional symbols;
c
c - a CHARACTER specifiying its position; example G stands for Gamma
c
c - a numeral, counting the number of equivalent positions (n); example
c   if there are two Gamma carbons (CG1,CG2 with HG1,HG2) 
c   if there is only ONE it is left out
c
c - a numeral, countings the number of hydrogens connected to a
c   position (m), if there is only ONE it CAN (but not nescessarilly)
c   be left out
c
c MOLDEN writes always mHn, while reading it first throws away the 'm'
c        and then starts counting internally
c
c Sybyl MOl2 may write anything of the following: Hn, Hm, Hnm
c
c example Hnm: VAL HG11,HG12,HG13
c example Hn:  PHE HE1
c example Hm:  PHE HB1,HB2
c
c Now when is it Hn and when is it Hm ?
c

c ipdbon is set independent of ioadd

      ipdbont = ipdbon
      dall = .false.

      trypdb = .true.
      if (istat.eq.-1) then
         trypdb = .false.
      endif

      toang = 0.52917706d0
      istat = 1

      call rewmf

      if (idebug.eq.1) print*,'rdmol'

      if (idebug.eq.1) print*,'@<TRIPOS>MOLECULE'
      call srchmf(line,'@<TRIPOS>MOLECULE',istat1)
      if (istat1.ne.1) goto 100

      imulm = 1

      if (ioadd.eq.1) then
         noff = iatoms
         ncoff = ncalf
         call numhet(nhmol)
         nmulm = nhmol - 4
         nmulmt = nmulm
         nhmol = nhmol + 1
         if (nmulmt.lt.1) nmulmt = 1
         ihasqt = ihasq
      else
         noff = 0
         ncoff = 0
         iatoms = 0
         ibnds = 0
         nhmol = 4
         nmulm = 0
         ihashb = 0
      endif

      ihasq = 0
      iff = 5
      ipdbon = 0
      inucl = 0
      ishoh = 0
      ncalf = 0
      iback = 0

      call inferr('Mol2 file',0)

      call redel(line,1)
      
      if (getlin(0).eq.1) then
          ktype = nxtwrd(string,nstr,itype,rtype)
          if (ktype.eq.2) then
              iatoms = itype
              if (iatoms.ge.mxnat) then
                 dall = .true.
                 goto 100
              endif
          else
              goto 100
          endif
          ktype = nxtwrd(string,nstr,itype,rtype)
          if (ktype.eq.2) then
              ibnds = itype
          else
              goto 100
          endif
          ktype = nxtwrd(string,nstr,itype,rtype)
          if (ktype.eq.2.and.trypdb) then
              ncalf = itype
          endif
      else
          goto 100
      endif

      if (getlin(0).eq.1) then
          ktype = nxtwrd(string,nstr,itype,rtype)
          if (ktype.eq.1) then
             if (index(string,'PROTEIN').ne.0) ipdbon = 1
             if (index(string,'BIOPOLYMER').ne.0) ipdbon = 1
             if (index(string,'NUCLEIC_ACID').ne.0) then
                 ipdbon = 1
                 inucl = 1
             endif
          endif
          if (ipdbon.eq.0) ncalf = 0
      endif

      if (.not.trypdb) then
          ipdbon = 0
          inucl = 0
      endif

      if (getlin(0).eq.1) then
          ktype = nxtwrd(string,nstr,itype,rtype)
          if (ktype.eq.1) then
             if (index(string,'NO_CHARGES').eq.0) ihasq = 1
          endif
      else
          goto 100
      endif

      if (ipdbon.eq.1) then
          call srchmf(line,'BACKBONE',istat1)
          if (istat1.eq.1) iback = 1
          call rewmf
      endif

      if (idebug.eq.1) print*,'@<TRIPOS>ATOM'
      call srchmf(line,'@<TRIPOS>ATOM',istat1)
      if (istat1.ne.1) goto 100
      if (idebug.eq.1.and.ipdbon.eq.1) print*,'@<TRIPOS> ipdbon'

      nfnd = 0
      irest = 0
      irst = -4
      do while (getlin(0).eq.1)
          if (idebug.eq.1) print*,'@',line
          ktype = nxtwrd(string,nstr,itype,rtype)
          if (ktype.eq.2) then

             nfnd = nfnd + 1
             iat = noff+itype
             ktype = nxtwrd(string,nstr,itype,rtype)

             if (ipdbon.eq.1) then
                if (ktype.eq.1) then
                   idonxt = 0
                   if (index(line,'HEM').ne.0) then
                      if (nstr.eq.1.and.string(1:1).eq.'N') then
                         idonxt = 1
                      endif
                   endif
                   if (nstr.eq.1) then
                      string(2:3) = '  '
                   elseif (nstr.eq.2) then
                      string(3:3) = ' '
                   elseif (nstr.eq.4) then
                      if (string.eq.'HCAP') string = 'HO3'
c HTER connected to O3P (or OXT) there is no H for it
                   endif
                   if (inucl.eq.1) then
                      if (ichar(string(3:3)).eq.39) string(3:3) = '*'
                   endif
                   ipdbt(iat) = 0
                   do j=1,mxsym
                      if (string(1:3).eq.pdbsym(j)) ipdbt(iat) = j
                   end do
                   if (ipdbt(iat).eq.0) then
                          
                      ic1 = ichar(string(1:1))

                      if (string(1:1).eq.'H'.or.(string(2:2).eq.'H'.and.
     &                    (ic1.ge.49.and.ic1.le.57))) then
                         if (ic1.ge.49.and.ic1.le.57) then
                            string = string(2:)
                            nstr = nstr - 1
                            ic2 = ichar(string(2:2))
                            ic3 = ichar(string(3:3))
                         else
c                            ihashb = 1
                            if (nstr.eq.3) then
                               if ((ic3.ge.49.and.ic3.le.57).and..not.
     &                             (ic2.ge.49.and.ic2.le.57)) then
                                  if (ic2.ne.72.and.ic2.ne.104) 
     &                               string(3:3) = ' '
                               endif
                            endif
                         endif
                         do j=1,mxhsym
                            if (string(1:3).eq.hsym(j)) then
                               ipdbt(iat) = (j-1)*3 + 1
                            endif
                         end do
                      endif
                   endif
                else
                   ipdbt(iat) = 0
                endif
                if (idonxt.eq.1) ktype = nxtwrd(string,nstr,itype,rtype)
             endif

             do i=1,3
                ktype = nxtwrd(string,nstr,itype,rtype)
                if (ktype.eq.3) then
                   coo(i,iat) = rtype/toang
                   iatclr(iat) = 1
                else
                   goto 100
                endif
             end do
             if (idebug.eq.1) print*,(coo(i,iat),i=1,3)
             iconn(1,iat) = 0

             string(1:5) = '     '
             ktype = nxtwrd(string,nstr,itype,rtype)
             if (ktype.eq.1) then
                if (nstr.eq.1) then
                   atomt(1:1) = ' '
                   atomt(2:2) = string(1:1)
                else
                   ic = ichar(string(2:2))
                   if (ic.ge.48.and.ic.le.57.or.string(2:2).eq.'.') 
     &             then
                      atomt(1:1) = ' '
                      atomt(2:2) = string(1:1)
                   else
                      atomt = string(1:2)
                   endif
                endif
                atom = tolowf(atomt)
                do j=1,mxel
                   if (atom .eq. tolowf(elemnt(j))) ianz(iat) = j
                end do
                if (atom.eq.'du') ianz(iat) = 99
                if (atom.eq.'lp') ianz(iat) = 99
                ityp(iat) = 0
                do j=1,mxmol2
                   if (string(1:5).eq.mol2(j)) ityp(iat) = j
                end do
                
             else
                goto 100
             endif

             ktype = nxtwrd(string,nstr,itype,rtype)
             if (ktype.eq.2) then

                if (ipdbon.eq.1) then
                   itype = itype + ncoff
                   iresid(iat) = itype
                   isres = 0

                   if ((iback.eq.1.and.icdex(line,'BACKBONE').ne.0).or.
     &                  iback.eq.0) then

                     if (ianz(iat).eq.1) then

                        if (ipdbt(iat).eq.1) icalf(4,itype) = iat

                     else

                        if (ipdbt(iat).eq.1) isres = 1
                        if (ipdbt(iat).eq.1) icalf(2,itype) = iat
                        if (ipdbt(iat).eq.2) icalf(1,itype) = iat
                        if (ipdbt(iat).eq.3) icalf(3,itype) = iat
                        if (ipdbt(iat).eq.43) icalf(1,itype) = iat
                        if (ipdbt(iat).eq.46) icalf(2,itype) = iat
                        if (ipdbt(iat).eq.47) icalf(3,itype) = iat
                        if (ipdbt(iat).eq.48) icalf(4,itype) = iat
                        if (ipdbt(iat).eq.50) icalf(5,itype) = iat
                        if (ipdbt(iat).eq.51) icalf(6,itype) = iat

                     endif

                   endif

                   if (icdex(line,'WATER').ne.0) then
                      if (ishoh.eq.0) then
                         if (nmulm.lt.mxmmul) nmulm = nmulm + 1
                         ishoh = nhmol
                         nhmol = nhmol + 1
                      endif
                      iresid(iat) = -ishoh
                   elseif (icdex(line,'DICT').ne.0) then
                      idum = 0
                   else
                      if (icdex(line,'<0>').ne.0) then
                         iresid(iat) = irst
                      else
                         if (iresid(iat).ne.irest) then
                            if (isres.eq.1) then
                               irst = itype
                            else
                               if (nmulm.lt.mxmmul) nmulm = nmulm + 1
                               iresid(iat) = -nhmol
                               irst = -nhmol
                               nhmol = nhmol + 1
                            endif
                         else
                               iresid(iat) = irst
                         endif
                      endif
                   endif

                   irest = itype
                else
                   iresid(iat) = itype
                   if (ioadd.eq.1.and.ipdbont.eq.1) then
                      if (iresid(iat).ne.irest) then
                         if (nmulm.lt.mxmmul) nmulm = nmulm + 1
                         call parsfn(gstr(nhmol),2,1)
                         iresid(iat) = -nhmol
                         irst = -nhmol
                         nhmol = nhmol + 1
                      else
                         iresid(iat) = irst
                      endif
                   else
                      iresid(iat) = -4
                   endif

                   irest = itype
                endif
             else
                goto 5
             endif

             ktype = nxtwrd(string,nstr,itype,rtype)
             if (ktype.eq.0) goto 5

             if (ihasq.eq.1) then

                ktype = nxtwrd(string,nstr,itype,rtype)
                if (ktype.eq.3) then
                   qat(iat) = rtype
                else
                   goto 100
                endif

             else
                qat(iat) = 0.0d0
             endif

5            continue

          elseif (ktype.eq.1.or.ktype.eq.0) then

             if (nfnd.ne.iatoms) goto 100
             if (ktype.eq.1) then
                if (index(string,'@<TRIPOS>').ne.0) then
                    call bckfil
                    goto 10
                endif
             endif

          endif
      end do

10    continue

      if (idebug.eq.1) print*,'@<TRIPOS>BOND'
      call srchmf(line,'@<TRIPOS>BOND',istat1)
      if (istat1.ne.1) goto 100

      if (idebug.eq.1) print*,'found bonds'

      ncon = 0
      ibnd = 0
      do i=1,ibnds
         if (getlin(0).eq.1) then
            ktype = nxtwrd(string,nstr,itype,rtype)
            ktype = nxtwrd(string,nstr,itype,rtype)
            if (ktype.eq.2) then
               iat1 = noff+itype
            else
               goto 100
            endif
            ktype = nxtwrd(string,nstr,itype,rtype)
            if (ktype.eq.2) then
               iat2 = noff+itype
            else
               goto 100
            endif

            ibnd = ibnd + 1
            if (iconn(1,iat1).lt.mxcon) then
               iconn(1,iat1) = iconn(1,iat1) + 1
               iconn(iconn(1,iat1)+1,iat1) = iat2
            endif
            if (iconn(1,iat2).lt.mxcon) then
               iconn(1,iat2) = iconn(1,iat2) + 1
               iconn(iconn(1,iat2)+1,iat2) = iat1
            endif

         else
            goto 100
         endif
      end do

      ncon = 1

      if (idebug.eq.1) print*,'redcon'
      if (ioadd.eq.0) call redcon

      call bckfil

      if (idebug.eq.1) print*,'@<TRIPOS>CRYSIN'
      call srchmf(line,'@<TRIPOS>CRYSIN',istat1)

      if (istat1.eq.1) then
         if (getlin(0).eq.1) then
             ktype = nxtwrd(string,nstr,itype,rtype)
             if (ktype.eq.3) a = rtype
             ktype = nxtwrd(string,nstr,itype,rtype)
             if (ktype.eq.3) b = rtype
             ktype = nxtwrd(string,nstr,itype,rtype)
             if (ktype.eq.3) c = rtype
             ktype = nxtwrd(string,nstr,itype,rtype)
             if (ktype.eq.3) alpha = rtype
             ktype = nxtwrd(string,nstr,itype,rtype)
             if (ktype.eq.3) beta = rtype
             ktype = nxtwrd(string,nstr,itype,rtype)
             if (ktype.eq.3) gamma = rtype
             ktype = nxtwrd(string,nstr,itype,rtype)
             if (ktype.eq.2) then
                if (itype.ge.0) then
                   nspg = itype
                else
                   goto 100
                endif
                call cprot(nspg,nopr,icent,ir,it,.false.)
                call prcell(nspg,a,b,c,alpha,beta,gamma)
                kz = 0
             endif
             ktype = nxtwrd(string,nstr,itype,rtype)
             if (ktype.eq.2) then
                if (itype.ge.1.and.itype.le.6) then
                   nset = itype
                   call aplset(nspg,nopr,ir,it)
                else
                   goto 100
                endif
             endif
         else
             return
         endif
          
         istat = 2
         ichx = 1
         if (idebug.eq.1) call prop(nopr,ir,it)
         call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)
         call cpmol(nat,norg,a,b,c,alpha,beta,gamma,
     &              coo,ianz,iatclr,iconn)

      endif

      if (ipdbon.eq.1) then

         call rewmf

         if (idebug.eq.1) print*,'@<TRIPOS>SUBSTRUCTURE'
         call srchmf(line,'@<TRIPOS>SUBSTRUCTURE',istat1)
         if (istat1.ne.0) then

            if (ioadd.eq.0) then
               call parsfn('Helix',5,1)
               call parsfn('Beta',4,1)
               call parsfn('RNA/DNA',7,1)
               call parsfn('Coil',4,1)
            endif

            ido = 1
c            do i=1,ncalf
            ncalf = 0
            do while (ido.eq.1)
               if (getlin(0).eq.1) then
                  ktype = nxtwrd(string,nstr,itype,rtype)
                  if (ktype.eq.1.and.nstr.ge.1) then
                     if (string(1:1).eq.'@') then
                        ido = 0
                        goto 200
                     endif
                  endif
                  if (ktype.eq.2) then
                     irstmp = ncoff + itype
                  else
                     ido = 0
                     goto 200
c                     goto 100
                  endif
                  iamino(irstmp) = 0

                  ktype = nxtwrd(string,nstr,itype,rtype)
                  if (ktype.eq.1) then
                     if (nstr.le.40) then
                        rstr = string(1:nstr)
                        nrstr = nstr
                     else
                        rstr = '****'
                        nrstr = 4
                     endif
                  endif

                  do j=1,2
                     ktype = nxtwrd(string,nstr,itype,rtype)
                  end do
                  if (ktype.eq.1) then
                     if (string(1:nstr).eq.'RESIDUE') then
                        do j=1,3
                           ktype = nxtwrd(string,nstr,itype,rtype)
                        end do
                        if (ktype.eq.1) then
                           if (nstr.ge.3) then
                              tstr = string(1:3)
                           elseif (nstr.eq.2) then
                              tstr(2:3) = string(1:2)
                              tstr(1:1) = ' '
                           elseif (nstr.eq.1) then
                              tstr(3:3) = string(1:1)
                              tstr(1:2) = '  '
                           endif
                           do j=1,mxres
                              if (tstr.eq.aminos(j)) 
     &                            iamino(irstmp) = j
                           end do
                           ncalf = ncalf + 1
                        else
                           print*,'Unknown Residue'
c                           goto 100
                        endif
                     else
                        call parsfn(rstr,nrstr,1)
                     endif
                  else
                     goto 100
                  endif

               else
                  ido = 0
               endif

200            continue
            end do

            if (ishoh.ne.0) then
               call parsfn('WATER',5,1)
            endif

c clean up, anything not a recognised residue

            ncalft = ncalf

            do i=ncoff+1,ncoff+ncalft
               if (iamino(i).le.0) then
                  ncalf = ncalf - 1
                  nhmol = nhmol + 1
                  if (nmulm.lt.mxmmul) nmulm = nmulm + 1
                  do j=ncoff+i,ncoff+ncalf
                     iamino(j) = iamino(j+1)
                     do k=1,6
                        icalf(k,j) = icalf(k,j+1)
                     end do
                  end do
                  iamino(ncoff+ncalf+1) = -1
                  do j=noff+1,noff+iatoms
                     if (iresid(j).eq.i) then
                        iresid(j) = -nhmol
                     elseif (iresid(j).gt.i) then
                        iresid(j) = iresid(j) - 1
                     endif
                  end do
               endif
            end do


c clean up, lp an du

            i = 0
            id = 0
            do while (i.le.iatoms)
               if (ianz(noff+i).eq.99) then
                  id = id + 1
                  iatoms = iatoms - 1

c loop over lp connecties haal by iconn de referenties naar lp weg

                  do j=1,iconn(1,noff+i)
                     jj = iconn(1+j,noff+i)
                     n = 0
                     do k=1,iconn(1,jj)
                        if (iconn(1+k,jj).ne.noff+i) then
                           n = n + 1
                           icont(n) = iconn(1+k,jj)
                        endif
                     end do
                     iconn(1,jj) = n
                     do k=1,n
                        iconn(1+k,jj) = icont(k)
                     end do
                  end do

                  do j=noff+i,noff+iatoms
                     ianz(j) = ianz(j+1)
                     ipdbt(j) = ipdbt(j+1)
                     iresid(j) = iresid(j+1)
                     ityp(j) = ityp(j+1)
                     qat(j) = qat(j+1)
                     do k=1,3
                        coo(k,j) = coo(k,j+1)
                     end do
                     do k=1,mxcon+1
                        iconn(k,j) = iconn(k,j+1)
                     end do
                  end do

c voor iconn alles groter dan i wordt i-1 

                  do j=noff+1,noff+iatoms
                     do k=1,iconn(1,j)
                        if (iconn(1+k,j).gt.noff+i) then
                           iconn(1+k,j) = iconn(1+k,j) - 1
                        endif
                     end do
                  end do

                  ianz(noff+iatoms+1) = 0

                  do j=1,ncalf
                     do k=1,6
                        if (icalf(k,j).gt.noff+i) 
     &                     icalf(k,j) = icalf(k,j) - 1
                     end do
                  end do

               else
                  i = i + 1
               endif
            end do

            if (id.gt.0) print*,'Cleaned up ',id,' dummy atoms'

c check for chain breaks

            if (ioadd.eq.1) then
               nchain = nchain + 1
               ianf(nchain) = ncoff+1
            else
               if (trypdb) then
                  nchain = 1
                  ianf(1) = 1
               else
                  nchain = 0
               endif
            endif

            do i=ncoff+1,ncoff+ncalf-1
               inewch = 0
               if (iamino(i).le.23) then
                  isami = 1
               else
                  isami = 0
               endif
               if (iamino(i+1).le.23) then
                  isnew = 1
               else
                  isnew = 0
               endif
               if (isami.ne.isnew) then
                  inewch = 1
               else
                  ifnd = 0
                  if (isami.eq.1) then
                     ii = icalf(3,i)
                     jj = icalf(2,i+1)
                     do j=1,iconn(1,ii)
                        if (iconn(1+j,ii).eq.jj) ifnd = 1
                     end do
                  else
                     ii = icalf(6,i)                   
                     jj = icalf(1,i+1)
                     do j=1,iconn(1,ii)
                        if (iconn(1+j,ii).eq.jj) ifnd = 1
                     end do
                  endif
                  if (ifnd.eq.0) inewch = 1
               endif
               if (inewch.eq.1.and.nchain.lt.mxchai) then
                  if (nchain.gt.0) islu(nchain) = i
                  nchain = nchain + 1
                  ianf(nchain) = i + 1
               endif
               reson(i) = 1
               isal(i) = 3
            end do
            islu(nchain) = ncoff+ncalf

c adjust H types

            do j=ncoff+1,ncoff+ncalf
               do i=1,mxhsym*3
                   ihpdb(i) = 0
               end do

               do i=noff+1,noff+iatoms
                  if (iresid(i).eq.j) then
                     if (ianz(i).eq.1) then
                        if (ihpdb(ipdbt(i)).eq.0) then
                           ihpdb(ipdbt(i)) = i
                        elseif (ihpdb(ipdbt(i)+1).eq.0) then
                           ipdbt(i) = ipdbt(i) + 1
                           ihpdb(ipdbt(i)) = i
                        elseif (ihpdb(ipdbt(i)+2).eq.0) then
                           ipdbt(i) = ipdbt(i) + 2
                           ihpdb(ipdbt(i)) = i
                        endif
                     endif
                  endif
               end do

               ia = iamino(j)

               do i=1,mxhsym*3
                  if (ihpdb(i).ne.0) then

                     itmp = i
                     iatmp = ihpdb(i)

                     if (ia.ge.17.and.ia.le.20) then
c HIS PHE TYR TRP
c 1HD,1HE
                        if (i.eq.19.or.i.eq.28) itmp = i + 3
c 2HD,2HE
                        if (i.eq.20.or.i.eq.29) itmp = i + 5
c TRP: 2HE -> HE3
                        if (ia.eq.20.and.i.eq.29) itmp = 37
                           
                     endif

                     if (ia.eq.20) then
c TRP: 1HZ -> HZ2, 2HZ -> HZ3
                        if (i.eq.40) itmp = 46
                        if (i.eq.41) itmp = 49
                     endif

c THR: 1HG -> HG1
                     if (ia.eq.5.and.i.eq.10) itmp = 13

                     if (itmp.ne.i) then
                         ihpdb(i) = 0
                         ihpdb(itmp) = iatmp
                         ipdbt(iatmp) = itmp
                     endif

                  endif
               end do
               
            end do

         else
            ipdbon = 0
         endif

      endif

      if (idebug.eq.1) print*,'almost done'

      if (ioadd.eq.1) then
         iatoms = iatoms + noff
         ncalf = ncalf + ncoff
         do i=nmulmt,nmulm
            ihasqm(i) = ihasq
         end do
         ihasq = ihasqt
      else
         do i=1,nmulm
            ihasqm(i) = ihasq
         end do
      endif

c      call hetnum(nhmolc)
c Make sure the total number of hetatm molecules is in sequence
c with the number parsed to C with parsfn, otherwise you get into trouble
c when adding new molecules
c
      call numhet(nhmolt)
      nhmold = nhmolt - nhmol
      if (nhmold.gt.0) then
         do i=1,nhmold
            call parsfn('?',1,1)
         end do
      endif

c      if (iatoms.gt.100.and.ipdbon.eq.0) call haszm(.true.)
      
      if (istat.gt.0.and.ipdbon.eq.1.and.ncalf.le.1) istat = 3

      iftyp = 6

      if (ihasq.eq.1) then
         qtot = 0.0d0
         do i=noff+1,iatoms
            qtot = qtot + qat(i)
         end do
         write(6,'(a,f6.3)') 'sum of charges =',qtot
      endif

      return

100   istat = 0
      if (ioadd.eq.0) then
         iatoms = 0
         iff = 0
      endif
      if (dall) then
         iatoms = noff
         ncalf = ncoff
         istat = -1
      endif

      return
      end


      subroutine chkmld(iok,ianz,ityp)
      parameter (mxmol2=41)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      save mol2ck
      integer*2 ityp
      common /athlp/ iatoms, mxnat
      dimension mol2ck(mxmol2),ianz(*),ityp(*)
      data mol2ck / 0,0,0,0,
     &    6,6,6,6,6,
     &    7,7,7,7,7,7,7,
     &    8,8,8,8,8,
     &    16,16,16,16,15,
     &    1,1,1,
     &    9,17,35,53,14,0,0,
     &    11,19,20,3,13/

      iok = 1

      imiss = 0
      do i=1,iatoms
         it = ityp(i)
         if (it.le.0) then
            iat1 = 0
         else
            iat1 = mol2ck(ityp(i))
         endif
         iat2 = ianz(i)
         if (iat1.ne.0) then
            if (iat1.ne.iat2) imiss = imiss + 1
         endif
      end do

      if (imiss.gt.0) iok = 0

      return
      end

      subroutine cprot(nspg,nopr,icent,ir,it,doprt)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxop=2822)
      parameter (mxrop=60)
      parameter (mxtop=56)
      parameter (mxsg=238)
      parameter (mxspg=232)
      integer*2 irr,itt
      common /symgrp/ irr(3,3,mxrop),itt(3,mxtop),nop(mxspg),iri(mxop),
     &                iti(mxop),ipt(mxspg),icen(mxspg),ispm(mxsg)
      character*7 spnam
      common /spgrnm/ spnam(mxsg)
      integer*2 ir,it
      dimension ir(3,3,192),it(3,192)
      logical doprt


      if (nspg.eq.0) then
         ispg = 0
      else
         ispg = ispm(nspg)
      endif

      if (ispg.ne.0) then
         nopr = nop(ispg)
         ip = ipt(ispg)

         do i=1,3
            it(i,1) = itt(i,1)
            do j=1,3
               ir(i,j,1) = irr(i,j,1)
            end do
         end do

         do k=1,nopr
            do i=1,3
               it(i,1+k) = itt(i,iti(ip+k-1))
               do j=1,3
                  ir(i,j,1+k) = irr(i,j,iri(ip+k-1))
               end do
            end do
         end do
         nopr = nopr + 1
         icent = icen(ispg)
         if (doprt) call inferr(spnam(nspg),0)
      else 
         if (nopr.eq.0) then
            nopr = 0
            call inferr('NO SpaceGroup Info',0)
         endif
      endif

      return
      end

      subroutine aplset(nspg,nopr,ir,it)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxsg=238)
      parameter (mxspg=232)
      integer*2 ir,it
      common /settng/ jset(3,9),iset(6,mxspg),nset
      character*7 spnam
      common /spgrnm/ spnam(mxsg)
      dimension itmp(3),irtmp(3,3),ir(3,3,192),it(3,192)


      do k=2,nopr
            
c copy operators

         do i=1,3
            itmp(i) = it(i,k)
            do j=1,3
               irtmp(i,j) = ir(i,j,k)
            end do
         end do

         do i=1,3
            inew = jset(i,iset(nset,nspg))
            it(i,k) = itmp(iabs(inew))
            
            do j=1,3
               jnew = jset(j,iset(nset,nspg))
               ir(i,j,k) = irtmp(iabs(inew),iabs(jnew))
               if (inew.lt.0) then
                  if (ir(i,j,k).eq.2) then
                      ir(i,j,k) = 0
                  elseif (ir(i,j,k).eq.0) then
                      ir(i,j,k) = 2
                  endif
               endif
               if (jnew.lt.0) then
                  if (ir(i,j,k).eq.2) then
                      ir(i,j,k) = 0
                  elseif (ir(i,j,k).eq.0) then
                      ir(i,j,k) = 2
                  endif
               endif
            end do
         end do
      end do

      return
      end

      subroutine cpmol(nat,norg,a,b,c,alpha,beta,gamma,
     &                 coo,ianz,iatclr,iconn)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      dimension tr1(3),rr(3,3)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*)

      toang = 0.52917706d0

      ca = dcos(alpha)
      cb = dcos(beta)
      cc = dcos(gamma)
      sc = dsin(gamma)
      sa = dsin(alpha)
      cfac = dsqrt(1.0d0-ca*ca-cb*cb-cc*cc+2.0d0*ca*cb*cc)
      rr(1,1) = 1.0d0 / a
      rr(1,2) = -cc/(a*sc)
      rr(1,3) = (ca*cc-cb) / (a*sc*cfac)
      rr(2,1) = 0.0d0
      rr(2,2) = 1.0d0 / (b*sc)
      rr(2,3) = (cb*cc-ca) / (b*sc*cfac)
      rr(3,1) = 0.0d0
      rr(3,2) = 0.0d0
      rr(3,3) = sc / (c*cfac)
   
      nat = iatoms
      norg = nat
      nstor = mxnat-iatoms

c      call dohcon(0)

      do i=1,nat

         do j=1,3
            tr1(j) = coo(j,i) * toang
         end do

         do k=1,3
           coo(k,nstor+i) = 0.0d0
           do j=1,3
              coo(k,nstor+i) = coo(k,nstor+i) + rr(k,j)*tr1(j)
           end do
         end do

         ianz(nstor+i) = ianz(i)

         do j=1,iconn(1,i)+1
            iconn(j,nstor+i) = iconn(j,i)
         end do

         iatclr(nstor+i) = iatclr(i)

      end do

      return
      end

      subroutine cpmol2(nat,norg,xa,ya,yb,za,zb,zc,
     &                 coo,ianz,iatclr,iconn)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      dimension tr1(3)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*)

      toang = 0.52917706d0

      nat = iatoms
      norg = nat
      nstor = mxnat-iatoms

c
c     Copy original atoms to top of coordinates array
c
      nstor = mxnat-iatoms
      do i=1,iatoms
         do j=1,3
            tr1(j) = coo(j,i) * toang
         end do
         call crt2fr(tr1,coo(1,nstor+i),xa,ya,yb,za,zb,zc)

         ianz(nstor+i) = ianz(i)
         do j=1,iconn(1,i)+1
            iconn(j,nstor+i) = iconn(j,i)
         end do
         iatclr(nstor+i) = iatclr(i)
      end do

      return
      end

      subroutine wrmod(iun,
     &               coo,qat,ianz,iaton,iatclr,iconn,iresid,
     &               lring,iamap,ityp,ipdbt,
     &               icalf,ncalf,iamino,ishoh,
     &               nat,nspg,icel,
     &               a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxel=100)
      parameter (mxt=14)
      parameter (mxppmf=16)
      parameter (mxlpmf=26)
      parameter (mxmol2=41)
      parameter (mxmm3=164)
      parameter (mxchtp=136)
      parameter (mxmsf=235)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxamo=201)
      parameter (mxres=42)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      common /athlp/ iatoms, mxnat
      character*2 elemnt,tocapf,atom
      common /elem/elemnt(mxel)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      character*6 atmp
      character*4 atype,stmp
      common /atypes/ ihbt(mxt),atype(mxt)
      integer*2 ityp,ipdbt
      common /types/ iff
      character*2 ppmf, lpmf
      character*5 mol2
      character*19 mm3
      character*20 chmtnk
      character*20 ambstr,dtmp
      character*20 amostr
      character*4 chmsf
      common /ftypes/ihasl(11),mol2(mxmol2),mm3(mxmm3),chmtnk(mxchtp),
     &               chmsf(mxmsf),ambstr(mxamb),amostr(mxamo),
     &               ppmf(mxppmf),lpmf(mxlpmf)
      logical dozme
      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme
      character*3 chtnk,pdbsym,hsym,ambtnk
      character*2 amotnk,gffstr
      common /symbol/ pdbsym(mxsym),hsym(mxhsym),chtnk(mxchtp),
     &                ambtnk(mxamb),amotnk(mxamo),gffstr(mxgff)
      character*3 aminos
      common /amino/aminos(mxres)
      real energy
      common /ener/iener,energy
      logical ochg
      dimension tr1(3),rr(3,3)
      dimension iamap(*)
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension coo(3,*),qat(*),ianz(*),iaton(*),iatclr(*),
     &          iconn(mxcon+1,*),iresid(*),ityp(*),ipdbt(*),lring(*)
      dimension icalf(6,*),iamino(*)
      data atype /'    ','.1  ','.2  ','.3  ','.4  ','.ar ','.cat',
     &            '.am ','.pl3','.co2','.spc','.t3p','.O  ','.O2 '/
      data ihbt /   3   ,  1   ,  2   ,  3   ,  3   ,  4   ,  4   ,
     &              3   ,  3   ,  2   ,  3   ,  3   ,  2   ,  2   /

      toang = 0.52917706d0
      torad = datan(1.0d0) / 45.0d0

c     eerst geen h-bonds

      jflg = 0
      natoms = iatoms
      nbnds = 0
      idochg = 0

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

      iatred = 0
      do i=1,natoms
         if (.not.(ianz(i).eq.100.and.
     &        (iresid(i).le.0.and.iresid(i).ge.-3))) then
            iatred = iatred + 1
            do j=1,iconn(1,i)
               if (iconn(1+j,i).gt.0) then
                  if (iconn(1+j,i).gt.i) nbnds = nbnds + 1
               endif
            end do
         endif
      end do
      natoms = iatred

      do i=1,natoms
         lring(i) = 0
         iaton(i) = 2
      end do


      write(iun,'(a)') '@<TRIPOS>MOLECULE'
      if (iener.eq.1) then
         write(iun,'(a,f9.3)') 'E=',energy
      else
         write(iun,*) 'Molden generated mol2'
      endif
      if (ipdbon.eq.1) then
         irs = ncalf
         ireso = 0
         do i=1,natoms
            ires = iresid(i)
            if (ires.lt.-3.and.ires.ne.-ishoh) then
               if (ires.ne.ireso) then
                  ireso = ires
                  irs = irs + 1
               endif
            endif
         end do
         write(iun,'(3(i6,1x))') natoms,nbnds,irs
         write(iun,*) 'PROTEIN'
      else
         write(iun,'(3(i6,1x))') natoms,nbnds,1
         write(iun,*) 'SMALL'
      endif
      if (idochg.eq.1) then
         write(iun,*) 'USER_CHARGES'
      else
         write(iun,*) 'NO_CHARGES'
      endif
      write(iun,*) '****'
      write(iun,*) '****'

      if (ipdbon.eq.1) then
         write(iun,'(a)') '@<TRIPOS>DICT'
         write(iun,*) 'PROTEIN PROTEIN'
      endif

      write(iun,'(a)') '@<TRIPOS>ATOM'

      if (icel.eq.1) then

         call setrr(alpha,beta,gamma,a,b,c,rr)
   
         do i=1,nat
            do k=1,3
               tr1(k) = trc(coo(1,i),rr,k)
            end do
            atom = tocapf(elemnt(ianz(i)))
            if (idochg.eq.1) then
               q = qat(i)
            else
               q = 0.0d0
            endif
            call ispn(irs,i,irng,idochg,0)
            atmp = atom//atype(irs)
            if (iff.eq.5) then
               if (ityp(i).le.0.or.ityp(i).gt.mxmol2) then
                  jflg = 1
                  atmp = ' '//mol2(1)
               else
                  atmp = ' '//mol2(ityp(i))
               endif
            endif
            write(iun,1000) i,atom,(tr1(j),j=1,3),atmp,q
         end do

      else

         if (ipdbon.eq.1) then

            do i=1,mxnat
               iamap(i) = 0
            end do 

            iat = 0
            irs = 0
            do i=1,ncalf
               irs = irs + 1

               call getpdb(i,ipdb,ihpdb)

c all non hydrogen residue atoms

               do j=1,mxsym

                ip = ipdb(j)

                if (ip.ne.0) then
                   iat = iat + 1
                   iamap(ip) = iat
                   if (idochg.eq.1) then
                      q = qat(ip)
                   else
                      q = 0.0d0
                   endif

                   atom = tocapf(elemnt(ianz(ip)))
                   call ispn(irss,ip,irng,idochg,0)
                   atmp = atom//atype(irss)

                   if (iff.eq.5) then
                      if (ityp(ip).le.0.or.ityp(ip).gt.mxmol2) then
                         jflg = 1
                         atmp = ' '//mol2(1)
                      else
                         atmp = ' '//mol2(ityp(ip))
                      endif
                   endif

                   iam = iamino(i)
                   if (iam.gt.0.and.iam.le.mxres) then
                      stmp = pdbsym(j)//' '
                   endif

                   if (j.ge.1.and.j.le.4) then
                       dtmp = 'BACKBONE|DICT|DIRECT'
                   else
                       dtmp = 'DICT'
                   endif

                   write(iun,1001) 
     &               iat,stmp,(coo(k,ip)*toang,k=1,3),atmp,i,
     &               aminos(iamino(i)),q,dtmp

                endif

               end do

c all hydrogen residue atoms

               do j=1,mxhsym*3

                  ihp = ihpdb(j)

                  if (ihp.ne.0) then
                     iat = iat + 1
                     iamap(ihp) = iat
                     if (idochg.eq.1) then
                        q = qat(ihp)
                     else
                        q = 0.0d0
                     endif
                     atom = tocapf(elemnt(ianz(ihp)))
                     call ispn(irss,ihp,irng,idochg,0)
                     atmp = atom//atype(irss)
                     if (iff.eq.5) then
                        if (ityp(ihp).le.0.or.ityp(ihp).gt.mxmol2) then
                           jflg = 1
                           atmp = ' '//mol2(1)
                        else
                           atmp = ' '//mol2(ityp(ihp))
                        endif
                     endif

                     ih = (j-1)/3 
                     il = j - ih*3
                     if (ihpdb(ih*3+2).eq.0) then
                        stmp = hsym(ih+1)//' '
                     else
                        if (hsym(ih+1)(3:3).ne.' ') then
                           stmp = hsym(ih+1)//char(48+il)
                        else
                           stmp = hsym(ih+1)//' '
                        endif
                     endif

                     dtmp = 'DICT'

                     write(iun,1001) 
     &               iat,stmp,(coo(k,ihp)*toang,k=1,3),atmp,i,
     &               aminos(iamino(i)),q,dtmp

                  endif

               end do

            end do

            igrp = 0
            ireso = 0
            do i=1,iatoms
               ires = iresid(i)
               if (ires.lt.-3.and.ianz(i).ne.100) then
                  iat = iat + 1
                  iamap(i) = iat
                  atom = tocapf(elemnt(ianz(i)))
                  if (idochg.eq.1) then
                     q = qat(i)
                  else
                     q = 0.0d0
                  endif
                  call ispn(irss,i,irng,idochg,0)
                  atmp = atom//atype(irss)
                  if (iff.eq.5) then
                     if (ityp(i).le.0.or.ityp(i).gt.mxmol2) then
                        jflg = 1
                        atmp = ' '//mol2(1)
                     else
                        atmp = ' '//mol2(ityp(i))
                     endif
                  endif
                  dtmp = '    '
                  if (ires.eq.-ishoh) dtmp = 'WATER'
                  if (ires.ne.ireso) then
                     ireso = ires
                     irs = irs + 1
                     igrp = igrp + 1
                     stmp = 'GRP'//char(48+igrp)
                  endif
                  write(iun,1002) 
     &              iat,atom,(coo(k,i)*toang,k=1,3),atmp,irs,
     &               stmp,q,dtmp
               endif
            end do

         else
            do i=1,iatoms
               atom = tocapf(elemnt(ianz(i)))
               if (idochg.eq.1) then
                  q = qat(i)
               else
                  q = 0.0d0
               endif
               call ispn(irs,i,irng,idochg,0)
               atmp = atom//atype(irs)
               if (iff.eq.5) then
                  if (ityp(i).le.0.or.ityp(i).gt.mxmol2) then
                     jflg = 1
                     atmp = ' '//mol2(1)
                  else
                     atmp = ' '//mol2(ityp(i))
                  endif
               endif
               if (.not.(ianz(i).eq.100.and.
     &              (iresid(i).le.0.and.iresid(i).ge.-3))) then
                  if (ipdbon.eq.1) then
                     ires = iresid(i)
                     if (ires.gt.0) then
                        iam = iamino(ires)
                        if (iam.gt.0.and.iam.le.mxres) then
                           if (ianz(i).eq.1) then
                              ih = (ipdbt(i)-1)/3
                              stmp = ' '//hsym(ih+1)
                           else
                              stmp = ' '//pdbsym(ipdbt(i))
                           endif
                        endif
                        if (ianz(i).ne.1.and.
     &                  (ipdbt(i).ge.1.and.ipdbt(i).le.4)) then
                            dtmp = 'BACKBONE|DICT|DIRECT'
                        else
                            dtmp = 'DICT'
                        endif
                        write(iun,1001) 
     &                  i,stmp,(coo(j,i)*toang,j=1,3),atmp,ires,
     &                  aminos(iamino(ires)),0.0d0,dtmp
                     else
                        write(iun,1000) 
     &                  i,atom,(coo(j,i)*toang,j=1,3),atmp,q
                     endif
                  else
                     write(iun,1000) 
     &                  i,atom,(coo(j,i)*toang,j=1,3),atmp,q
                  endif
               endif
            end do

         endif

      endif


      write(iun,'(a)') '@<TRIPOS>BOND'

      ibnds = 1

      if (ipdbon.eq.1) then

         do i=1,natoms
            do j=1,iconn(1,i)
               k = iconn(1+j,i)
               if (k.gt.0) then
                  if (k.gt.i) then
                     ibt = ibtyp(i,k,idochg,0,ianz)
                     kt = iamap(k)
                     it = iamap(i)
                     if (ibt.eq.2) then
                        ifl1 = 0
                        ifl2 = 0
                        do jj=1,iconn(1,i)
                           kk = iconn(1+jj,i)
                           if (kk.gt.0.and.kk.ne.k) then
                               ibth = ibtyp(i,kk,idochg,0,ianz)
                               if (ibth.eq.2.or.ibth.eq.4) ifl1 = 1
                           endif
                        end do
                        do jj=1,iconn(1,k)
                           kk = iconn(1+jj,k)
                           if (kk.gt.0.and.kk.ne.i) then
                               ibth = ibtyp(k,kk,idochg,0,ianz)
                               if (ibth.eq.2.or.ibth.eq.4) ifl2 = 1
                           endif
                        end do
                        if (ifl1.eq.1.and.ifl2.eq.1) ibt = 1
                     endif
                     if (ibt.eq.4) then
                        write(iun,'(3(i6,1x),a)') ibnds,it,kt,'    ar'
                     else
                        write(iun,'(4(i6,1x))') ibnds,it,kt,ibt
                     endif
                     ibnds = ibnds + 1
                  endif
               endif
            end do
         end do

      else

         do i=1,natoms
            do j=1,iconn(1,i)
               k = iconn(1+j,i)
               if (k.gt.0) then
                  if (k.gt.i) then
                     ibt = ibtyp(i,k,idochg,0,ianz)
                     if (ibt.eq.2) then
                        ifl1 = 0
                        ifl2 = 0
                        do jj=1,iconn(1,i)
                           kk = iconn(1+jj,i)
                           if (kk.gt.0.and.kk.ne.k) then
                               ibth = ibtyp(i,kk,idochg,0,ianz)
                               if (ibth.eq.2.or.ibth.eq.4) ifl1 = 1
                           endif
                        end do
                        do jj=1,iconn(1,k)
                           kk = iconn(1+jj,k)
                           if (kk.gt.0.and.kk.ne.i) then
                               ibth = ibtyp(k,kk,idochg,0,ianz)
                               if (ibth.eq.2.or.ibth.eq.4) ifl2 = 1
                           endif
                        end do
                        if (ifl1.eq.1.and.ifl2.eq.1) ibt = 1
                     endif
                     if (ibt.eq.4) then
                        write(iun,'(3(i6,1x),a)') ibnds,i,k,'    ar'
                     else
                        write(iun,'(4(i6,1x))') ibnds,i,k,ibt
                     endif
                     ibnds = ibnds + 1
                  endif
               endif
            end do
         end do

      endif

      write(iun,'(a)') '@<TRIPOS>SUBSTRUCTURE'
      if (ipdbon.eq.1) then
         irs = 0
         do i=1,ncalf
            irs = irs + 1
            write(iun,'(i5,1x,a4,i5,a,a3)') 
     &       i,aminos(iamino(i))//' ',iamap(icalf(1,i)),
     &       ' RESIDUE    1 A   ',aminos(iamino(i))
         end do
         igrp = 0
         ireso = 0
         do i=1,natoms
            ires = iresid(i)
            if (ires.lt.-3.and.ires.ne.-ishoh) then
               if (ires.ne.ireso) then
                  ireso = ires
                  irs = irs + 1
                  igrp = igrp + 1
                  stmp = 'GRP'//char(48+igrp)
                  write(iun,'(i5,1x,a4,i5,a)') 
     &             irs,stmp,iamap(i),
     &             ' GROUP      0 A     ****    0 ROOT'
               endif
            endif
         end do
      else
         write(iun,*) '     1 RES1       1'
      endif

      if (icel.eq.1) then
         write(iun,'(a)') '@<TRIPOS>CRYSIN'
         write(iun,'(6(f10.4,1x),i3,1x,i1)') 
     &         a,b,c,alpha/torad,beta/torad,gamma/torad,nspg,1
      endif

      call chkmol2(iok)
      if (iok.eq.0) call messg(17)
      if (jflg.eq.1) call messg(15)

      if (icel.eq.1) call fdat(1,0,0,0,0,0)

      do i=1,natoms
         iaton(i) = 1
      end do

1000  format(i6,1x,a2,1x,3(f10.4,1x),a6,' 1 RES1',f10.4)
1001  format(i6,1x,a4,1x,3(f10.4,1x),a6,1x,i4,1x,a3,f10.4,1x,a)
1002  format(i6,1x,a4,1x,3(f10.4,1x),a6,1x,i4,1x,a4,f10.4,1x,a)
      return
      end

      subroutine flth(iat,icnn,ibnds,iconn)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      dimension iconn(mxcon+1,*),icnn(mxcon)

      ibnds = 0
      do i=1,iconn(1,iat)
         if (iconn(i+1,iat).gt.0) then
            ibnds = ibnds + 1
            icnn(ibnds) = iconn(i+1,iat)
         endif
      end do

      return
      end

      subroutine ispnd(ispn,iat,irng,idochg,ifive,
     &                      qat,ianz,iaton,iconn,lwrit,lring)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/  iatoms, mxnat
      logical ringg
      real angle,dihed
      dimension iring(6),isel(4),icnn(mxcon),icnn2(mxcon),ioke(2)
      dimension qat(*),iconn(mxcon+1,*),ianz(*),iaton(*),
     &          lwrit(*),lring(*)

c ifive not only signals look for five membered aromatic rings
c but also c==0 counts for no Pi electrons
c
c determine number of bonds and filter out hydrogen bonds

      ibnds = 0
      inn = 0
      do i=1,iconn(1,iat)
         if (iconn(i+1,iat).gt.0) then
            ibnds = ibnds + 1
            icnn(ibnds) = iconn(i+1,iat)
            if (ianz(iconn(i+1,iat)).eq.7) inn = inn + 1
         endif
      end do

      ispn = 1
      irng = 0
      ian = ianz(iat)
      if (ian.eq.6) then
         ispn = 4
         if (ibnds.eq.3) then
            ispn = 3
            if (inn.eq.3) ispn = 7
         elseif (ibnds.eq.2) then
            ispn = 2
         endif
      elseif (ian.eq.7) then
         ispn = 4
         if (ibnds.eq.4) then
            ispn = 5
         elseif (ibnds.eq.1) then
            ispn = 2
         elseif (ibnds.eq.2) then
            isel(1) = icnn(1)
            isel(2) = iat
            isel(3) = icnn(2)
            ispn = 2
            call intcor(intc,angle,isel,3)
            if (intc.eq.1) then
                if (abs(angle).lt.170.0e0) ispn = 3
            endif
         elseif (ibnds.eq.3) then
            ispn = 4
            isel(1) = iat
            isel(2) = icnn(1)
            isel(3) = icnn(2)
            isel(4) = icnn(3)
            call intcor(intc,dihed,isel,4)
            if (intc.eq.1) then
                if (abs(dihed).lt.3.0e0) ispn = 9
            endif
            inco = 0
            do i=1,3
               jat = icnn(i)
               call flth(jat,icnn2,ibnds2,iconn)
               if (ianz(jat).eq.6.and.ibnds2.eq.3) then
                   do j=1,3
                      kat = icnn2(j)
                      if (ianz(kat).eq.8) then
                         call flth(kat,icnn2,ibnds2,iconn)
                         if (ibnds2.eq.1) inco = jat
                      endif
                   end do
               endif
            end do
            if (inco.ne.0) then
               iamid = 1
               do i=1,3
                  jat = icnn(i)
                  if (jat.ne.inco) then
                     if (.not.(ianz(jat).eq.1.or.ianz(jat).eq.6.or.
     &                         ianz(jat).eq.7))
     &                 iamid = 0
                  endif
               end do
               if (iamid.eq.1) ispn = 8
            endif
         endif
      elseif (ian.eq.8.or.ian.eq.16) then
         ispn = 4
         if (ibnds.eq.1) then
            ispn = 3
            jat = icnn(1)
            if (ian.eq.8.and.
     &         (ianz(jat).eq.15.or.ianz(jat).eq.6)) then
               ipo2 = 0
               call flth(jat,icnn2,ibnds2,iconn)
               do i=1,ibnds2
                  kat = icnn2(i)
                  if (kat.gt.0) then
                     if (ianz(kat).eq.8.and.
     &                icred(kat,idum1,idum2,ianz,iconn).eq.1) 
     &                      ipo2 = ipo2 + 1
                  endif
               end do
               if (ipo2.eq.2) ispn = 10
            endif
         endif
         if (ian.eq.16) then
            iso2 = 0
            do i=1,ibnds
               call flth(icnn(i),icnn2,ibnds2,iconn)
               if (ianz(icnn(i)).eq.8.and.ibnds2.eq.1)
     &            iso2 = iso2 + 1
            end do
            if (iso2.eq.1) ispn = 13
            if (iso2.eq.2) ispn = 14
         endif
      elseif (ian.eq.15) then
         ispn = 4
      endif

      if (ian.ne.6.and.ian.ne.7.and.ian.ne.8.and.ian.ne.15
     &    .and.ian.ne.16) return

      do i=1,iatoms
         lring(i) = 0
      end do

      do ii=1,2
         ioke(ii) = 0
         narel = 0
         if (ringg(iat,iring,nring,.true.,ianz,iaton,iconn,
     &             lwrit,lring)) then
            irng = -1
            if (nring.eq.6.or.nring.eq.5) then
               do i=1,nring
                  jat = iring(i)
                  iann = ianz(jat)
                  if (iann.eq.6) then
                     narel = narel + 1
                     if (ifive.eq.1) then
                        call flth(jat,icnn2,ibnds2,iconn)
                        do j=1,ibnds2
                           kat = icnn2(j)
                           if (ianz(kat).eq.8.or.
     &                         ianz(kat).eq.16) then
                              if (icred(kat,idum1,idum2,ianz,iconn)
     &                            .eq.1) narel = narel - 1
                           endif
                        end do
                     endif
                  elseif (iann.eq.8.or.iann.eq.16) then
                     narel = narel + 2
                  elseif (iann.eq.15) then
                     narel = narel + 1
                  elseif (iann.eq.7) then
                     call flth(jat,icnn2,inb,iconn)
                     if (inb.eq.2.or.inb.eq.3) then
                        if (inb.eq.3) then
                           if (idochg.eq.1) then
                              if (abs(qat(jat)-1.0d0).gt.0.1d0) then
                                 narel = narel + 2
                              else
                                 narel = narel + 1
                              endif
                           else
                              narel = narel + 2
                           endif
                        elseif (inb.eq.2) then
                           narel = narel + 1
                        endif
                     else
                        ioke(ii) = 0
                     endif
                  endif
               end do
            endif
            if ((nring.eq.6.or.(nring.eq.5.and.ifive.eq.1))
     &         .and.narel.eq.6) ioke(ii) = nring
            do i=1,nring
               lring(iring(i)) = 1
            end do
         endif
      end do
      if (ioke(1).gt.0.or.ioke(2).gt.0) then
          if (ian.ne.8) ispn = 6
          if (ioke(1).eq.6.or.ioke(2).eq.6) irng = 1
          if (ioke(1).eq.5.or.ioke(2).eq.5) then
              irng = 2
              if (ian.eq.16) ispn = 3
          endif
          if (ioke(1).eq.6.and.ioke(2).eq.6) then
              irng = 3
          elseif (ioke(1).eq.5.and.ioke(2).eq.5) then
              irng = 4
          elseif (ioke(1).gt.0.and.ioke(2).gt.0) then
              irng = 5
          endif
      endif

      return
      end 

      subroutine setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,ipr)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /frc/ ifr2crt

      torad = datan(1.0d0) / 45.0d0
      alpha = alpha*torad
      beta = beta*torad
      gamma = gamma*torad

      ca = dcos(alpha)
      cb = dcos(beta)
      cc = dcos(gamma)
      sc = dsin(gamma)
      sa = dsin(alpha)
      cvol = a*b*c*
     &   dsqrt(1-ca**2-cb**2-cc**2+2.0d0*ca*cb*cc)

      if (ifr2crt.eq.0) then
         xa = cvol / (b*c*sa)
         ya = (a*(cc-(ca*cb)))/sa
         yb = (b*sa)
         za = (a*cb)
         zb = (b*ca)
         zc = c
      else
         xa = a
         ya = b*cc
         za = c*cb
         yb = b*sc
         zb = (c*(ca-(cb*cc)))/sc
         zc = cvol / (a*b*sc)
      endif

      if (ipr.eq.1) then
         print*,' '
         write(*,'(a,f9.3)') ' Cell Volume = ',cvol
      endif

      return
      end

      subroutine prop(nopr,ir,it)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      integer*2 ir,it
      dimension ir(3,3,192),it(3,192)

      do i=1,nopr
         print*,'ir',i
         do j=1,3
            print*,(ir(j,k,i),k=1,3)
         end do
         print*,'it',i
         print*,(it(j,i),j=1,3)
      end do

      return
      end

      subroutine cprop(nopr,ir,it)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxop=2822)
      parameter (mxrop=60)
      parameter (mxtop=56)
      parameter (mxsg=238)
      parameter (mxspg=232)
      integer*2 ir,it
      integer*2 irr,itt
      common /symgrp/ irr(3,3,mxrop),itt(3,mxtop),nop(mxspg),iri(mxop),
     &                iti(mxop),ipt(mxspg),icen(mxspg),ispm(mxsg)
      dimension irs(100),its(100),ir(3,3,192),it(3,192)

      do i=1,nopr
       irs(i) = 0
       do l=1,mxrop
         ifnd = l
         do j=1,3
            do k=1,3
               if (ir(j,k,i).ne.irr(j,k,l)) ifnd = 0
            end do
         end do
         if (ifnd.ne.0) then
            irs(i) = l
            goto 10
         endif
       end do
10     if (irs(i).eq.0) then
            print*,'ir',i
            do j=1,3
               print*,(ir(j,k,i),k=1,3)
            end do
       endif

       its(i) = 0
       do l=1,mxtop
         ifnd = l
         do k=1,3
             if (it(k,i).ne.itt(k,l)) ifnd = 0
         end do
         if (ifnd.ne.0) then
            its(i) = l
            goto 20
         endif
       end do
20     if (its(i).eq.0) then
         print*,'it',i
         print*,(it(j,i),j=1,3)
       endif
      end do

      print*,'nopr=',nopr-1
      print*,'irs'
      write(6,'(100i3)')(irs(i),i=2,nopr)
      print*,'its'
      write(6,'(100i3)')(its(i),i=2,nopr)

      return
      end

      logical function chkrec(nopr,ir,it)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      integer*2 ir,it
      dimension ir(3,3,192),it(3,192)


      chkrec = .true.
      do j=1,3
         do k=1,3
               if (j.eq.k.and.ir(j,k,nopr).ne.2) chkrec =.false. 
               if (j.ne.k.and.ir(j,k,nopr).ne.1) chkrec =.false. 
         end do
         if (it(j,nopr).ne.0) chkrec =.false.
      end do
      if (chkrec) return

      return
      end

      subroutine recsym(nopr,ir,it,keywin,ilatt)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      integer*2 ir,it
      character*(*) keywin
      character*80 ky
      character*80 keywrd
      character*1 ksym
      integer*2 itr
      logical chkrec
      dimension ky(3),n(3),ksym(3),itr(3),ir(3,3,192),it(3,192)

      keywrd = keywin(5:min0(len(keywin),80))
  
      ky(1) = ' '
      ky(2) = ' '
      ky(3) = ' '

      do i=1,2
         i1 = index(keywrd,',')
         if (i1.gt.1.and.i1.lt.len(keywrd)) then
           ky(i) = keywrd(1:min0(i1-1,80))
           keywrd = keywrd(i1+1:)
         else
           goto 100
         endif
      end do
      ky(3) = keywrd(1:min0(80,len(keywrd)))

      nopr = nopr + 1

      ksym(1) = 'X'
      ksym(2) = 'Y'
      ksym(3) = 'Z'

      do i=1,3
         call spatrm(ky(i),n(i))
         do j=1,3
            ir(i,j,nopr) = 1
            k = index(ky(i)(1:n(i)),ksym(j))
            if (k.ne.0) then
               ky(i)(k:k) = ' '
               ir(i,j,nopr) = 2
               if (k.gt.1) then
                  if (ky(i)(k-1:k-1).eq.'-') then
                     ir(i,j,nopr) = 0
                     ky(i)(k-1:k-1) = ' '
                  endif
                  if (ky(i)(k-1:k-1).eq.'+') then
                     ky(i)(k-1:k-1) = ' '
                  endif
               endif
            endif
         end do
         ks = index(ky(i)(1:n(i)),'/')
         if (ks.ne.0) then
            i1 = ichar(ky(i)(ks-1:ks-1))
            i2 = ichar(ky(i)(ks+1:ks+1))
            if ((i1.ge.48.and.i1.le.57).and.
     &          (i2.ge.48.and.i2.le.57)) then
                i1 = i1 - 48
                i2 = i2 - 48
                tt = dble(i1)/dble(i2)*12.0d0
                it(i,nopr) = krnd(tt)
            endif
         else
            tt = reada(ky(i),1,n(i))*12.0d0
            it(i,nopr) = krnd(tt)
         endif
      end do

      if (chkrec(nopr,ir,it)) then
         nopr = nopr - 1
         return
      endif

      latt = iabs(ilatt)

      if (latt.eq.2) then
         itr(1) = 6
         itr(2) = 6
         itr(3) = 6
         call symcpr(ir,nopr,nopr+1)
         call symcpt(it,itr,nopr,nopr+1)
         nopr = nopr + 1
      elseif (latt.eq.4) then
         noprt = nopr
         do i=1,3
            itr(1) = 6
            itr(2) = 6
            itr(3) = 6
            itr(i) = 0
            call symcpr(ir,noprt,nopr+1)
            call symcpt(it,itr,noprt,nopr+1)
            nopr = nopr + 1
         end do
      elseif (latt.ge.5) then
         itel = latt - 4
         itr(1) = 6
         itr(2) = 6
         itr(3) = 6
         itr(itel) = 0
         call symcpr(ir,nopr,nopr+1)
         call symcpt(it,itr,nopr,nopr+1)
         nopr = nopr + 1
      endif

100   return
      end

      subroutine symrec(nopr,icent,ir,it,iun,doall)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      character*2 sopx,sopy,sopz
      common /opcom/ sopx(3),sopy(3),sopz(3)
      character*12 st
      logical setcen,doall,ldum
      integer*2 itr,ir,it
      dimension st(3),rst(3),lst(3),itr(3),iceneq(192)
      dimension ir(3,3,192),it(3,192)
      data sopx/'-X','  ','+X'/
      data sopy/'-Y','  ','+Y'/
      data sopz/'-Z','  ','+Z'/

      ilatt = 1
      ldum = setcen(iceneq,nopr,0,ndum)

      itr(1) = 6
      itr(2) = 6
      itr(3) = 6

      call symiop(iceneq,itr,nopr,ir,it)

      if (setcen(iceneq,nopr,2,nred)) then
c
c I ?
c
         ilatt = 2

      else


c
c F ?
c
         ldum = setcen(iceneq,nopr,0,ndum)
         do i=1,3
            itr(1) = 6
            itr(2) = 6
            itr(3) = 6
            itr(i) = 0
            call symiop(iceneq,itr,nopr,ir,it)
         end do
         if (setcen(iceneq,nopr,2,nred)) then
            if (nred.eq.4) then
               ilatt = 4
               goto 10
            endif
         endif
c
c A,B,C ?
c
         do i=1,3
            ldum = setcen(iceneq,nopr,0,ndum)
            itr(1) = 6
            itr(2) = 6
            itr(3) = 6
            itr(i) = 0
            call symiop(iceneq,itr,nopr,ir,it)
            if (setcen(iceneq,nopr,2,nred)) then
                ilatt = i + 4
                goto 10
            endif
         end do

c
c P (Dont know how to do 3 yet)
c
         ilatt = 1
         ldum = setcen(iceneq,nopr,1,ndum)

      endif
      
10    continue

      if (icent.eq.2) then
         ilatt = -1*ilatt
      endif
      write(iun,'(a,i2)') 'LATT ',ilatt

      ist = 2
      if (doall) then
         ldum = setcen(iceneq,nopr,1,ndum)
         ist = 1
      endif

      do i=ist,nopr
         
        if (iceneq(i).eq.1) then

         do k=1,3
            st(k) = '            '
            st(k) = 
     &        sopx(ir(k,1,i)+1)//sopy(ir(k,2,i)+1)//sopz(ir(k,3,i)+1)

            do j=1,6
               if (st(k)(j:j).ne.' ') then
                  st(k) = st(k)(j:12)
                  goto 100
               endif
            end do
100         if (st(k)(1:1).eq.'+') st(k) = st(k)(2:12)
            lst(k) = 6
            do j=12,1,-1
               if (st(k)(j:j).ne.' ') then
                  lst(k) = j
                  goto 200
               endif
            end do

200         continue
            if (it(k,i).ne.0) then
               rst(k) = dfloat(it(k,i))/12.0d0
               if (rst(k).lt.0.0d0) then
                  st(k)(lst(k)+1:lst(k)+1) = '-'
               else
                  st(k)(lst(k)+1:lst(k)+1) = '+'
               endif
               write(st(k)(lst(k)+2:lst(k)+6),'(f5.3)') abs(rst(k))
               lst(k) = lst(k)+6
            endif
         end do

         write(iun,'(a)') 'SYMM '//st(1)(1:lst(1))//','//
     &                 st(2)(1:lst(2))//','//st(3)(1:lst(3))

        endif

      end do

      return
      end

      logical function setcen(iceneq,nopr,iop,nred)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      dimension iceneq(*)

      setcen = .true.

      if (iop.eq.0) then
          do i=1,nopr
             iceneq(i) = 0
          end do
      elseif (iop.eq.1) then
          do i=1,nopr
             iceneq(i) = 1
          end do
      elseif (iop.eq.2) then
          nunq = 0
          do i=1,nopr
             if (iceneq(i).eq.0) setcen = .false.
             if (iceneq(i).eq.1) nunq = nunq + 1
          end do
          nred = 1
          if (nunq.ne.0) nred = nopr/nunq
      endif

      return
      end 

      subroutine symiop(iceneq,itr,nopr,ir,it)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      integer*2 ir,it,itr
      logical symeqr,symeqt
      dimension itr(3),iceneq(*),ir(3,3,192),it(3,192)

c     iceneq
c
c      1  has equal
c      0  has no equal
c     -1  centering equal
      
      do i=1,nopr
         if (iceneq(i).ge.0) then
            do j=i+1,nopr
               if (symeqr(ir,i,j)) then
                  if (symeqt(it,itr,i,j)) then
                     iceneq(i) = 1
                     iceneq(j) = -1
                  endif
               endif
            end do
         endif
      end do

      return
      end

      logical function symeqr(ir,iopr,jopr)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      integer*2 ir
      dimension ir(3,3,192)

      symeqr = .false.

      do k=1,3
         do j=1,3
            if (ir(k,j,iopr).ne.ir(k,j,jopr)) return
         end do
      end do

      symeqr = .true.

      return
      end

      logical function symeqt(it,itr,iopr,jopr)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      integer*2 it,jt,itr
      dimension itr(3),it(3,192)

      symeqt = .false.

      do i=1,3
         jt = (it(i,iopr) + itr(i)) - it(i,jopr)
         if (jt.ne.0.and.jt.ne.12) return
      end do

      symeqt = .true.

      return
      end

      subroutine symcpr(ir,iopr,jopr)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      integer*2 ir
      dimension ir(3,3,192)

      do k=1,3
         do j=1,3
            ir(k,j,jopr) = ir(k,j,iopr)
         end do
      end do

      return
      end

      subroutine symcpt(it,itr,iopr,jopr)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      integer*2 it,jt,itr
      dimension itr(3),it(3,192)

      do i=1,3
         jt = it(i,iopr) + itr(i)
         if (jt.eq.12) jt = 0
         it(i,jopr) = jt
      end do

      return
      end

      subroutine prcell(nspg,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxsg=238)
      character*7 spnam
      common /spgrnm/ spnam(mxsg)

      print*,' '
      print*,'Unit Cell Parameters:'
      print*,' '
      write(*,1080) ' a     = ',a
      write(*,1080) ' b     = ',b
      write(*,1080) ' c     = ',c
      write(*,1080) ' alpha = ',alpha
      write(*,1080) ' beta  = ',beta
      write(*,1080) ' gamma = ',gamma
      print*,' '
      print*,'Space Group Number = ',nspg
      if (nspg.gt.0) print*,'Space Group Symbol  = ',spnam(nspg)
      print*,' '

1080  format(a,f7.3)
      return
      end

      logical function ochg(idum,ianz)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      dimension ianz(*)

      ochg = .true.
      if (natoms.eq.0) ochg = .false.
      do i=1,natoms
         if (nat(i).ne.ianz(i)) ochg = .false.
      end do

      return
      end

      integer function ibtyp(iat,jat,idochg,ifive,ianz)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxt=14)
      character*4 atype
      common /atypes/ ihbt(mxt),atype(mxt)
      dimension ianz(*)

      call ispn(irs,iat,irng,idochg,ifive)
      ihb1 = ihbt(irs)
      call ispn(irs,jat,irng,idochg,ifive)
      ihb2 = ihbt(irs)
      ibtyp = 1
      if (ihb1.eq.2.and.ihb2.eq.2) ibtyp = 2
      if (ihb1.eq.1.and.ihb2.eq.1) ibtyp = 3
      if (ihb1.eq.4.and.ihb2.eq.4) ibtyp = 4
      if ((ihb1.eq.4.and.ihb2.eq.2).or.(ihb1.eq.2.and.ihb2.eq.4)) 
     &   ibtyp = 2
      if (ihb1.eq.3.and.ianz(iat).eq.15
     &    .and.ihb2.eq.2.and.ianz(jat).eq.8) ibtyp = 2
      if (ihb2.eq.3.and.ianz(jat).eq.15
     &    .and.ihb1.eq.2.and.ianz(iat).eq.8) ibtyp = 2

      return
      end

      double precision function trc(coo,rr,k)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      dimension coo(3),rr(3,3)

      trc = 0.0d0
      do j=1,3
         trc = trc + rr(k,j)*coo(j)
      end do

      return
      end

      subroutine setrr(alpha,beta,gamma,a,b,c,rr)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      dimension rr(3,3)

      ca = dcos(alpha)
      cb = dcos(beta)
      cc = dcos(gamma)
      sc = dsin(gamma)
      sa = dsin(alpha)
      cfac = dsqrt(1.0d0-ca*ca-cb*cb-cc*cc+2.0d0*ca*cb*cc)
      rr(1,1) = a
      rr(1,2) = b*cc
      rr(1,3) = c*cb
      rr(2,1) = 0.0d0
      rr(2,2) = b*sc
      rr(2,3) = c*(ca-cb*cc)/sc
      rr(3,1) = 0.0d0
      rr(3,2) = 0.0d0
      rr(3,3) = c*cfac/sc

      return
      end

      subroutine calcij(nuse,coo,ianz,iresid,iconn,qat)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxel=100)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /buick/ aij(mxel),bij(mxel),cij(mxel),aijhp,bijhp,cijhp
      logical dopair,debug
      dimension coo(3,*),iresid(*),ianz(*),iconn(mxcon+1,*),qat(*)

      if (ihasq.eq.0) then
         print*,'No charges available'
         return
      endif

      debug = .true.
c      ev2kcl = 0.043363d0 old value (1/23.061d0)
      ev2kcl = 23.061d0
      hr2kcl = 627.5095d0
      toang = 0.52917706d0


c nuse = 4, exclude 1-2,1-3,1-4 
c nuse = 0, only pairs from different residues

      elec = 0.0d0
      elec1 = 0.0d0
      evdw = 0.0d0

      print*,' '
      print*,'======================================='

      do i=1,iatoms

         do j=i+1,iatoms

             dopair = .true.
             if (nuse.eq.0) then
                if (iresid(i).eq.iresid(j).or.
     &          (iresid(i).eq.0.or.iresid(j).eq.0)) dopair = .false.
             else
c 1-2
                do k=1,iconn(1,j)
                   kk = iconn(1+k,j)
                   if (kk.eq.i.and.nuse.ge.2) dopair = .false.
                   if (kk.gt.0) then
c 1-3
                      do l=1,iconn(1,kk)
                         ll = iconn(1+l,kk)
                         if (ll.eq.i.and.nuse.ge.3) dopair = .false.
                         if (ll.gt.0.and.ll.ne.j) then
c 1-4
                            do m=1,iconn(1,ll)
                               mm = iconn(1+m,ll)
                               if (mm.eq.i.and.nuse.ge.4) 
     &                            dopair = .false.
                            end do
                         endif
                      end do
                   endif
                end do
             endif

             if (dopair) then

                r2 = 0.0d0
                do k=1,3
                   tmp = coo(k,i) - coo(k,j)
                   r2 = r2 +tmp*tmp
                end do
                if (r2.gt.0.0d0) then
                   r1 = dsqrt(r2)
c
c electrostatic
c
                   if (ihasq.eq.1) elec = elec + qat(i)*qat(j)/r1
 
c vdw
                   r1 = r1 * toang
                   r2 = r2 * toang * toang

                   alk = dsqrt(aij(ianz(i))*aij(ianz(j)))
                   blk = (bij(ianz(i)) + bij(ianz(j))) / 2.0d0
                   clk = dsqrt(cij(ianz(i))*cij(ianz(j)))
                   evdw = evdw + alk*dexp(-blk*r1) - clk/(r2*r2*r2)
                endif

             endif

         end do

      end do

      elec = elec * hr2kcl
      evdw = evdw * ev2kcl

      if (nuse.eq.0) then
         print*,'Only pair-energy of atoms from different residues'
      elseif (nuse.eq.2) then
         print*,'Exclude pair-energy of 1-2 neighbour atoms'
      elseif (nuse.eq.3) then
         print*,'Exclude pair-energy of 1-2,1-3 neighbour atoms'
      elseif (nuse.eq.4) then
         print*,'Exclude pair-energy of 1-2,1-3,1-4 neighbour atoms'
      endif
      print*,' '

      print*,'Electrostatic energy = ',elec
      print*,'Van der Waals energy = ',evdw
      print*,'Total         energy = ',evdw+elec

      return
      end

      subroutine setpp(idoest)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxel=100)
      parameter (mxtyp=8)
      character*2 labtyp
      common /buck/ alk(mxtyp,mxtyp),blk(mxtyp,mxtyp),clk(mxtyp,mxtyp),
     &              alkp(mxtyp),blkp(mxtyp),clkp(mxtyp),
     &              iptyp(mxtyp),labtyp(mxtyp)
      common /buick/ aij(mxel),bij(mxel),cij(mxel),aijhp,bijhp,cijhp

      aijhp = 52.128991d0
c      bijhp = 0.214592d0
      bijhp = 4.660d0
      cijhp = 0.222819d0

      do i=1,mxtyp
         do j=1,mxtyp

            ai = aij(iptyp(i))
            if (i.eq.2) ai = aijhp
            aj = aij(iptyp(j))
            if (j.eq.2) aj = aijhp
            if (ai.gt.0.0d0.and.aj.gt.0.0d0) then
               alk(i,j) = dsqrt(ai*aj)
            else
               alk(i,j) = 0.0d0
            endif

            ci = cij(iptyp(i))
            if (i.eq.2) ci = cijhp
            cj = cij(iptyp(j))
            if (j.eq.2) cj = cijhp
            if (ci.gt.0.0d0.and.cj.gt.0.0d0) then
               clk(i,j) = dsqrt(ci*cj)
            else
               clk(i,j) = 0.0d0
            endif

            bi = bij(iptyp(i))
            if (i.eq.2) bi = bijhp
            bj = bij(iptyp(j))
            if (j.eq.2) bj = bijhp
            if (bi.gt.0.0d0.and.bj.gt.0.0d0) then
               blk(i,j) = 2.0d0 / (bi+bj)
            else
               blk(i,j) = 0.0d0
            endif

         end do
      end do

      if (idoest.eq.1) then
         do i=1,mxtyp
            alk(2,i) = alkp(i)
            blk(2,i) = blkp(i)
            clk(2,i) = clkp(i)
         end do
         do i=1,mxtyp
            alk(i,2) = alk(2,i)
            blk(i,2) = blk(2,i)
            clk(i,2) = clk(2,i)
         end do
      endif

      return
      end

      subroutine wrcryd(iun,ianz,coo,
     &                  nat,nspg,icel,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (tol1=1.0d-3)
      parameter (tol2=1.0d-2)
      common /athlp/ iatoms, mxnat
      character*80 cllstr
      character*7 cella,cellb,cellc
      dimension ianz(*),coo(3,*)

      toang = 0.52917706d0
      torad = datan(1.0d0) / 45.0d0


      if (icel.ne.1) then
         call inferr('No Cell Data !',0)
         return
      endif

      natoms = iatoms
      natoms = nat
      nstor = mxnat-natoms

      write(iun,'(a)') 'Molden generated Crystal95 input'
      write(iun,'(a)') 'CRYSTAL'
      write(iun,'(a)') '0 0 0'
      write(iun,'(i3)') nspg
     
      write(cella,'(f7.4)') a
      write(cellb,'(f7.4)') b
      write(cellc,'(f7.4)') c

      cllstr = cella
      ncll = 7
      if (dabs(a-b).gt.tol1) then
         cllstr = cllstr(1:ncll)//' '//cellb 
         ncll = ncll + 8
      endif

      if (dabs(a-c).gt.tol1.and.dabs(b-c).gt.tol1) then
         cllstr = cllstr(1:ncll)//' '//cellc 
         ncll = ncll + 8
      endif

      dalp = alpha/torad
      dbet = beta/torad
      dgam = gamma/torad
      write(cella,'(f7.3)') dalp
      write(cellb,'(f7.3)') dbet
      write(cellc,'(f7.3)') dgam

      if (dabs(dalp-90.0d0).gt.tol2) then
         cllstr = cllstr(1:ncll)//' '//cella 
         ncll = ncll + 8
      endif

      if (dabs(dbet-90.0d0).gt.tol2) then
         cllstr = cllstr(1:ncll)//' '//cellb 
         ncll = ncll + 8
      endif

      if (dabs(dgam-90.0d0).gt.tol2) then
         cllstr = cllstr(1:ncll)//' '//cellc 
         ncll = ncll + 8
      endif

      write(iun,'(a)') cllstr(1:ncll)

      write(iun,'(i4)') natoms

      do i=1,natoms
         write(iun,'(i2,1x,3(f12.6,1x))') 
     &         ianz(nstor+i),(coo(j,nstor+i),j=1,3)
      end do

      write(iun,'(a)') 'END'

      call basprt(iun,.false.,.true.)
      write(iun,'(a)') 'END'

c     gen. inf. 
      write(iun,'(a)') 'RHF'
      write(iun,'(a)') 'MULTISCF'
      write(iun,'(a)') 'END'

c     scf input
      write(iun,'(a)') 'GUESSPAT'
      write(iun,'(a)') 'END'

      call inferr('Wrote file: '//'crystal95.in',0)

      return
      end

      subroutine rdshld(idebug,istat,
     &                  coo,ianz,iconn,iatclr,
     &                  nat,icent,ncon,nspg,kz,nopr,ir,it,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxel=100)
      parameter (mxsg=238)
      common /athlp/ iatoms, mxnat
      integer*2 ir,it
      character*2 elemnt
      common /elem/elemnt(mxel)
      integer getlin
      character*137 line,str
      common /curlin/ line
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*7 spnam,spntmp,spgru
      common /spgrnm/ spnam(mxsg)
      logical gnreal,sfac,chkrec
      dimension isfac(mxel),ctmp(7),tmp(3)
      dimension ianz(*),coo(3,*),iconn(mxcon+1,*),iatclr(*)
      dimension ir(3,3,192),it(3,192)

      call rewfil
      istat = 2

      toang = 0.52917706d0
      sfac = .false.

      if (getlin(0).eq.1) then
         ic = icdex(line,'TITL')
         ic1 = icdex(line,'REM')
         if (ic.eq.0.and.ic1.eq.0) goto 100
         if (line(ic+4:ic+5).eq.'e]'.or.line(ic+4:ic+5).eq.'E]')
     &       goto 100
      else
         goto 100
      endif

      nspg = 0

      kz = 0
      ilatt = 0
      isymm = 0
      ispgr = 0
      icent = 1
      icnt = 0
      nopr = 0
      nat = 0

      call recsym(nopr,ir,it,'SYMM X,Y,Z',1)

      do while (getlin(0).eq.1)

         if (icdex(line,'LATT').ne.0) then

             ktype = nxtwrd(str,nstr,itype,rtype)
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.2) then
                if (itype.lt.0) icent = 2
                ilatt = iabs(itype)
                if (ilatt.eq.3) call inferr(
     &          'WARNING: LATT 3 not supported, missing operators!',1) 
             else
                goto 100
             endif
        
         elseif (icdex(line,'SYMM').ne.0) then

             if (ilatt.eq.0) ilatt = 1
c             if (isymm.eq.0) then
c                call recsym(nopr,ir,it,'SYMM X,Y,Z',ilatt)
c                isymm = 1
c             endif
             call tocap(line,len(line))
             call recsym(nopr,ir,it,line,ilatt)

         elseif (icdex(line,'CELL').ne.0) then

             ktype = nxtwrd(str,nstr,itype,rtype)

             i7 = 1
             do i=1,7
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.3) then
                   ctmp(i) = rtype
                elseif (ktype.eq.2) then
                   ctmp(i) = dble(itype)
                elseif (ktype.eq.0) then
                   i7 = 0
                endif
             end do

             a = ctmp(1+i7)
             b = ctmp(2+i7)
             c = ctmp(3+i7)
             alpha = ctmp(4+i7)
             beta  = ctmp(5+i7)
             gamma = ctmp(6+i7)


         elseif (icdex(line,'SPGR').ne.0) then

             ispgr = 1
             ktype = nxtwrd(str,nstr,itype,rtype)
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.1) then
                 spntmp = '       '
                 spntmp(1:nstr) = str(1:nstr)
                 call tocap(spntmp,7)
                 ifound = 0
                 do i=1,mxsg
                    spgru = spnam(i)
                    call tocap(spgru,7)
                    if (spntmp.eq.spgru) then
                       nspg = i
                       ifound = 1
                    endif
                 end do
                 if (ifound.eq.0)
     &              call inferr('No spacegroup match',0)

             endif

         elseif (icdex(line,'SFAC').ne.0) then

c             sfac = .true.
c problem with ICSD generated spf files
c
             sfac = .false.
             ktype = nxtwrd(str,nstr,itype,rtype)
             do while (nxtwrd(str,nstr,itype,rtype).eq.1)
                if (nstr.eq.1.or.nstr.eq.2) then
                   icnt = icnt + 1
                   isfac(icnt) = iatnum(str,nstr)
                endif
             end do

         elseif (icdex(line,'REM').ne.0) then
             idum = 1
         else

            iattmp = 0
            jatom = 0
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.1) then
               if (nstr.eq.4) then
                  if (str(1:4).eq.'ATOM'.or.str(1:4).eq.'atom') 
     &            then
                     ktype = nxtwrd(str,nstr,itype,rtype)
                     jatom = 1
                  endif
               endif
               if (.not.sfac) then
                  if (nstr.gt.1) then
                     ic = ichar(str(2:2))
                     if ((ic.ge.97.and.ic.le.122).or.
     &                   (ic.ge.65.and.ic.le.90)) then
                        iattmp = iatnum(str,2)
                     else
                        iattmp = iatnum(str,1)
                     endif
                  else
                     iattmp = iatnum(str,1)
                  endif
                  if (nstr.gt.2) then
                     ic = ichar(str(3:3))
                     if ((ic.ge.97.and.ic.le.122).or.
     &                   (ic.ge.65.and.ic.le.90)) iattmp = 0
                  endif
               endif
               if (jatom.eq.1) then
                  ktype = 2
               else
                  ktype = nxtwrd(str,nstr,itype,rtype)
               endif
               if (ktype.eq.2.and.(
     &            (sfac.and.itype.le.icnt.and.itype.gt.0) .or.
     &            (.not.sfac.and.iattmp.gt.0.and.iattmp.lt.mxel)
     &            )) then
                  if (gnreal(tmp,3,.false.)) then
                     nat = nat + 1
                     if (sfac) then
                        ianz(nat) = isfac(itype)
                     else
                        ianz(nat) = iattmp
                     endif
                     
                     do j=1,3
                        coo(j,nat) = tmp(j)
                     end do
                  endif
               elseif (ktype.eq.3.and.(
     &            (sfac.and.itype.le.icnt.and.itype.gt.0) .or.
     &            (.not.sfac.and.iattmp.gt.0.and.iattmp.lt.mxel)
     &            )) then
                  nat = nat + 1
                  ianz(nat) = iattmp
                  tmp(1) = rtype
                  if (gnreal(tmp(2),2,.false.)) then
                     do j=1,3
                        coo(j,nat) = tmp(j)
                     end do
                  else
                     nat = nat - 1
                  endif
               endif
            endif

         endif

      end do

      if (nat.eq.0) goto 100

      call prcell(nspg,a,b,c,alpha,beta,gamma)
      call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)

      if ((ilatt.eq.0.or.ilatt.eq.1).and.ispgr.eq.0) then
         ispgr = 1
         nspg = 0
      endif

      if (ispgr.eq.1) then
           icent = 2
           call cprot(nspg,nopr,icent,ir,it,.false.)
           if (idebug.eq.1) call prop(nopr,ir,it)
      endif

      nstor = mxnat-nat

      do i=1,nat
         do j=1,3
            coo(j,nstor+i) = coo(j,i)
         end do
         iconn(1,i) = 0
         call fr2crt(coo(1,i),xa,ya,yb,za,zb,zc)
         do j=1,3
             coo(j,i) = coo(j,i)/toang
         end do
      end do

      iatoms = nat

      call doconn
      call dohcon(0)

      do i=1,nat
         ianz(nstor+i) = ianz(i)
         do j=1,iconn(1,i)+1
            iconn(j,nstor+i) = iconn(j,i)
         end do
         iatclr(nstor+i) = 1
      end do

      ncon = 1

      if (idebug.eq.1) call prop(nopr,ir,it)

      return

100   istat = 0
      return
      end

      subroutine wrshld(iun,idospf,coo,ianz,
     &                  nat,icent,nspg,kz,icel,nopr,ir,it,
     &                  a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxel=100)
      parameter (mxsg=238)
      common /athlp/ iatoms, mxnat
      integer*2 ir,it
      character*7 spnam
      common /spgrnm/ spnam(mxsg)
      character*2 elemnt,tocapf
      character*3 tstr
      common /elem/elemnt(mxel)
      character*80 sfac,ustr
      dimension isfac(2,mxel),coo(3,*),ianz(*)
      dimension ir(3,3,192),it(3,192)

      toang = 0.52917706d0
      torad = datan(1.0d0) / 45.0d0


      if (icel.ne.1) then
         call inferr('No Cell Data !',0)
         return
      endif

      call fdat(4,0,0,0,0,0)

      nstor = mxnat-nat

      dalp = alpha/torad
      dbet = beta/torad
      dgam = gamma/torad

      if (idospf.eq.1) then
         write(iun,'(a)') 'TITL Molden generated SPF file'
         write(iun,'(a,3(f7.4,1x),3(f7.2,1x))') 
     &                 'CELL ',a,b,c,dalp,dbet,dgam
      else
         write(iun,'(a)') 'TITL Molden generated SHELX input'
         write(iun,'(a,3(f7.4,1x),3(f7.2,1x))') 
     &                 'CELL 0.0 ',a,b,c,dalp,dbet,dgam
      endif

      if (idospf.eq.1) then
         write(iun,'(a)') 'SPGR '//spnam(nspg)
      else
         kzz = kz
         if (kz.eq.0) then
             kzz = krnd(dfloat((iatoms-4))/dfloat(nat))
         endif
         write(iun,'(a,i2,a)') 'ZERR ',kzz,' 0.0 0.0 0.0 0.0 0.0 0.0'

         call symrec(nopr,icent,ir,it,iun,.false.)

         do i=1,mxel
            isfac(1,i) = 0
         end do
         do i=1,iatoms
            if (ianz(i).ne.100) then
               isfac(1,ianz(i)) = isfac(1,ianz(i)) + 1
            endif
         end do
         sfac = 'SFAC'
         ustr = 'UNIT'
         icnt = 0
         do i=1,mxel
            if (isfac(1,i).ne.0) then
                icnt = icnt + 1
                isfac(2,i) = icnt
                sfac =  sfac(1:5+(icnt-1)*3)//' '//tocapf(elemnt(i))
                write(tstr,'(i3)') isfac(1,i)
                ustr =  ustr(1:5+(icnt-1)*4)//' '//tstr
            endif
         end do
         write(iun,'(a)') sfac
         write(iun,'(a)') ustr

      endif

      do i=1,nat
         if (idospf.eq.1) then
            write(iun,'(a,a2,1x,3(f12.6,1x))') 'ATOM ',
     &       tocapf(elemnt(ianz(nstor+i))),
     &       (coo(j,nstor+i),j=1,3)
         else
            write(iun,'(a2,1x,i2,3(f12.6,1x))') 
     &       tocapf(elemnt(ianz(nstor+i))),isfac(2,ianz(nstor+i)),
     &       (coo(j,nstor+i),j=1,3)
         endif
      end do

      if (idospf.eq.0) write(iun,'(a)') 'HKLF 1'
      write(iun,'(a)') 'END'

      if (idospf.eq.1) then
        call inferr('Wrote file: pluton.spf',0)
      else 
        call inferr('Wrote file: '//'shelx.ins',0)
      endif

      return
      end

      subroutine wrcifd(iun,coo,ianz,
     &                  nat,nspg,icel,
     &                  a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxel=100)
      parameter (mxsg=238)
      common /athlp/ iatoms, mxnat
      character*7 spnam
      common /spgrnm/ spnam(mxsg)
      character*2 elemnt,tocapf,el
      character*3 tstr
      common /elem/elemnt(mxel)
      dimension coo(3,*),ianz(*)

      toang = 0.52917706d0
      torad = datan(1.0d0) / 45.0d0


      if (icel.ne.1) then
         call inferr('No Cell Data !',0)
         return
      endif

      call fdat(4,0,0,0,0,0)

      nstor = mxnat-nat

      dalp = alpha/torad
      dbet = beta/torad
      dgam = gamma/torad

      write(iun,'(a)') 'data_Molden'
      write(iun,'(a)') '_symmetry_space_group_name_H-M '//spnam(nspg)
      write(iun,'(a,i4)') '_symmetry_Int_Tables_number ',nspg 
      write(iun,'(a,f7.4)') '_cell_length_a ',a 
      write(iun,'(a,f7.4)') '_cell_length_b ',b 
      write(iun,'(a,f7.4)') '_cell_length_c ',c 
      write(iun,'(a,f7.4)') '_cell_angle_alpha ',dalp 
      write(iun,'(a,f7.4)') '_cell_angle_beta ',dbet 
      write(iun,'(a,f7.4)') '_cell_angle_gamma ',dgam 

      write(iun,'(a)') 'loop_'
      write(iun,'(a)') '_atom_site_label'
      write(iun,'(a)') '_atom_site_type_symbol'
      write(iun,'(a)') '_atom_site_fract_x'
      write(iun,'(a)') '_atom_site_fract_y'
      write(iun,'(a)') '_atom_site_fract_z'

      do i=1,nat
           el = elemnt(ianz(nstor+i))
           write(iun,'(a2,1x,a2,1x,3(f10.4,1x))') 
     &       el,el,(coo(j,nstor+i),j=1,3)
      end do

      write(iun,'(a)') '#END'

      call inferr('Wrote file: mol.cif',0)

      return
      end

      subroutine clini
      implicit double precision (a-h,o-z)
      common /cllmat/rx(3),ry(3),rz(3),t(3),tz(3),tzorg(3),itz

      do i=1,3
         rx(i) = 0.0d0
         ry(i) = 0.0d0
         rz(i) = 0.0d0
      end do

      rx(1) = 1.0d0
      ry(2) = 1.0d0
      rz(3) = 1.0d0

      return
      end

      subroutine clmat(inct,theang,iopt)
c
c actually inirot +xyzrot in rott.f
c
      implicit double precision (a-h,o-z)
      common /cllmat/rx(3),ry(3),rz(3),t(3),tz(3),tzorg(3),itz

      todeg = 45.0d0 / datan(1.0d0)

      sa = dsin(theang/todeg)
      ca = dcos(theang/todeg)

      if (iopt.eq.1) call clini

      if (inct.eq.-3) then
         call rarbx(theang/todeg)
      elseif (inct.eq.-2) then
         call rarby(theang/todeg)
      elseif (inct.eq.-1) then
         call rarbz(theang/todeg)
      endif

      return

c This is the old way of rotating the ligand/mol2, but the rotation
c axes of the ligand/mol2, rotated along with the overall rotation

      if (inct.eq.-3) then
c rotation around x-axis (A-axis)
         do i=1,3
            y = ry(i)
            z = rz(i)
            ry(i) = ca*y + sa*z
            rz(i) = ca*z - sa*y
         end do
      elseif (inct.eq.-2) then
c rotation around y-axis (B-axis)
         do i=1,3
            x = rx(i)
            z = rz(i)
            rx(i) = ca*x + sa*z
            rz(i) = ca*z - sa*x
         end do
      elseif (inct.eq.-1) then
c rotation around z-axis (C-axis in cell)
         do i=1,3
            x = rx(i)
            y = ry(i)
            rx(i) = ca*x - sa*y
            ry(i) = ca*y + sa*x
         end do
      endif

      return
      end

      subroutine cllrod(vec,irot,ifd,coo,ianz,
     &                  nat,xa,ya,yb,za,zb,zc)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      common /cllmat/rx(3),ry(3),rz(3),t(3),tz(3),tzorg(3),itz
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt
      dimension vec(3),coo(3,*),ianz(*)

      ctol = 1.0d-6
      iatoms = nat

c copy original from top

      nstor = mxnat-iatoms

      do i=1,iatoms
         do j=1,3
            coo(j,i) = coo(j,nstor+i)
            if (abs(coo(j,i)).lt.ctol) coo(j,i) = 0.0d0
         end do
         ianz(i) = ianz(nstor+i)
      end do

c convert fractional to cartesian

      do i=1,iatoms
         call fr2crt(coo(1,i),xa,ya,yb,za,zb,zc)
      end do

      if (irot.eq.1) then

c calculate center of mass

         call cntvec(t,coo,ianz,iatoms)

c do rotation

         do i=1,iatoms
            x = coo(1,i)
            y = coo(2,i)
            z = coo(3,i)
            coo(1,i) = (x-t(1))*rx(1)+(y-t(2))*rx(2)+(z-t(3))*rx(3)+t(1)
            coo(2,i) = (x-t(1))*ry(1)+(y-t(2))*ry(2)+(z-t(3))*ry(3)+t(2)
            coo(3,i) = (x-t(1))*rz(1)+(y-t(2))*rz(2)+(z-t(3))*rz(3)+t(3)
         end do

      else

c do translation

         do i=1,iatoms
            do j=1,3
               coo(j,i) = coo(j,i) - vec(j)
            end do
         end do

      endif

      do i=1,iatoms
         call crt2fr(coo(1,i),coo(1,nstor+i),xa,ya,yb,za,zb,zc)
      end do

      call fdat(ifd,1,0,0,0,0)
      call docent
      if (iatoms.eq.iscst) iatoms = iatoms + nscnd

      return
      end

      subroutine zm2fd(cm,ctmp,imkeep,coo,ianz,iatclr,iconn,
     &                 nat,ichx,xa,ya,yb,za,zb,zc)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxsite=300)
      parameter (lmax=4)
      parameter (lmaxf=(lmax+1)*(lmax+1))
      common /athlp/ iatoms, mxnat
      integer dolabs,fancy,persp,shade,atcol,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb,nions,ntota,nresi
      dimension imkeep(3),i3(3),tr1(3)
      dimension cm(3,3),ctmp(3,*)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*)

c      if (.not.ichx.eq.1) return

c  misuse common multip as scratch

      toang = 0.52917706d0
      ctol = 1.0d-6
      nstor = mxnat-nat

      if (ichx.eq.1) then
         do i=1,nat
            do j=1,3
               ctmp(j,i) = coo(j,nstor+i)
               if (abs(ctmp(j,i)).lt.ctol) ctmp(j,i) = 0.0d0
            end do
            call fr2crt(ctmp(1,i),xa,ya,yb,za,zb,zc)

            do j=1,3
                ctmp(j,i) = ctmp(j,i)/toang
            end do

         end do
      endif

      do j=1,3
         i3(j) = j
      end do

      if (ichx.eq.1) then
         do i=1,3
             tr1(i) = ctmp(i,imkeep(1)) - coo(i,1)
         end do
      else
         do i=1,3
             tr1(i) = cm(i,1) - coo(i,1)
         end do
      endif

      do i=1,iatoms
         call trcoo(tr1,coo(1,i))
      end do

      if (ichx.eq.1) then
         call alntwo(ctmp,imkeep,coo,iatoms,i3)
      else
         call alntwo(cm,i3,coo,iatoms,i3)
         return
      endif

      call doconn
      call dohcon(0)

      nat = iatoms
      nstor = mxnat-nat

      do i=1,iatoms

         do j=1,3
            tr1(j) = coo(j,i) * toang
         end do
         call crt2fr(tr1,coo(1,nstor+i),xa,ya,yb,za,zb,zc)

         ianz(nstor+i) = ianz(i)

         do j=1,iconn(1,i)+1
            iconn(j,nstor+i) = iconn(j,i)
         end do
         iatclr(nstor+i) = 1

      end do

      call fdat(ifd,1,0,0,0,0)
      call docent
      call upzme

      return
      end


      subroutine chgpar(ianz,coo,
     &                  nat,icent,nspg,ichx,nopr,ir,it,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,o-z)
      integer*2 ir,it
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb,nions,ntota,nresi
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /athlp/ iatoms, mxnat
      dimension vec(3),vecfr(3),coo(3,*),ianz(*)
      dimension ir(3,3,192),it(3,192)

      if (.not.ichx.eq.1) return

c  misuse common multip as scratch

c      toang = 0.52917706d0
      ctol = 1.0d-6
      nstor = mxnat-nat

      do i=1,nat
         do j=1,3
            coo(j,i) = coo(j,nstor+i)
            if (abs(coo(j,i)).lt.ctol) coo(j,i) = 0.0d0
         end do
         call fr2crt(coo(1,i),xa,ya,yb,za,zb,zc)

      end do

      call cntvec(vec,coo,ianz,nat)
      call crt2fr(vec,vecfr,xa,ya,yb,za,zb,zc)

      if (nspg.gt.0) call cprot(nspg,nopr,icent,ir,it,.true.)

      if (nspg.ge.1.and.nspg.le.2) then

c triclinic
          idum = 1

      elseif ((nspg.ge.3.and.nspg.le.15).or.nspg.eq.231) then

c monoclinic (B axis unique)
          alpha = 90.0d0
          gamma = 90.0d0

      elseif (nspg.ge.16.and.nspg.le.74) then

c orthorhombic
          alpha = 90.0d0
          beta = 90.0d0
          gamma = 90.0d0

      elseif (nspg.ge.75.and.nspg.le.142) then

c tetragonal
          b = a
          alpha = 90.0d0
          beta = 90.0d0
          gamma = 90.0d0

      elseif ((nspg.ge.143.and.nspg.le.167).or.
     &        (nspg.ge.232.and.nspg.le.238)) then

c trigonal (hexagonal cell)

          b = a
          alpha = 90.0d0
          beta = 90.0d0
          gamma = 120.0d0

c      elseif 146,148,155,160,161,166,167 then (7 in all)
c all R groups in the range 143-167 can have an alternative
c (rhombohedral cell), spacegroup name is appended with an small r
c so you have R3 (hexagonal cell) and R3r (rhombohedral cell)
c
cc trigonal (rhombohedral cell)
c          b = a
c          c = a
c          beta = alpha
c          gamma = alpha

      elseif (nspg.ge.168.and.nspg.le.194) then

c hexagonal
          b = a
          alpha = 90.0d0
          beta = 90.0d0
          gamma = 120.0d0

      elseif (nspg.ge.195.and.nspg.le.230) then

c cubic
          b = a
          c = a
          alpha = 90.0d0
          beta = 90.0d0
          gamma = 90.0d0

      endif

      call prcell(nspg,a,b,c,alpha,beta,gamma)
      call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)
      print*,' '
c     call symrec(nopr,icent,ir,it,iun3,.true.)
c what you need is printing of centre of symmtry operations + unity

      call fr2crt(vecfr,xa,ya,yb,za,zb,zc)

      do j=1,3
         vec(j) = vecfr(j) - vec(j)
      end do

      do i=1,nat
         call trcoo(vec,coo(1,i))
         call crt2fr(coo(1,i),coo(1,nstor+i),xa,ya,yb,za,zb,zc)
      end do


      return
      end

      subroutine cllvec(isel,inum,coo)
      implicit double precision (a-h,o-z)
      common /cllmat/rx(3),ry(3),rz(3),t(3),tz(3),tzorg(3),itz
      dimension isel(4),coo(3,*)

      if (inum.ge.2) then
         do i=1,3
            tz(i) = coo(i,isel(2)) - coo(i,isel(1))
            if (inum.eq.4) 
     &         tz(i) = tz(i) + (coo(i,isel(4)) - coo(i,isel(3)))
         end do
         call vsc1(tz,0.001d0*0.52917706d0,1.0d-4)
         itz = 1
      endif
      
      return
      end

      subroutine mkcell(coo,ianz,iatclr,iconn,
     &                  nat,icent,inorm,nspg,nopr,ir,it,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      integer*2 ir,it
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      dimension vec(3),abc(3)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*)
      dimension ir(3,3,192),it(3,192)

      toang = 0.52917706d0

      nstor = mxnat-iatoms
      nat = iatoms

      call cntvec(vec,coo,ianz,iatoms)
      do i=1,3
         abc(i) = 0.0d0
      end do
      do i=1,iatoms
         do j=1,3
            t = dabs(coo(j,i) - vec(j))
            if (t.gt.abc(j)) abc(j) = t
         end do
      end do

      do i=1,3
         abc(i) = 2.0d0*abc(i) + 3.0d0
      end do
      a = abc(1)*toang
      b = abc(2)*toang
      c = abc(3)*toang
      alpha = 90.0d0
      beta = 90.0d0
      gamma = 90.0d0
      call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)

      inorm = 0
      nspg = 1
      call cprot(nspg,nopr,icent,ir,it,.false.)
      
      do i=1,iatoms

         do j=1,3
             coo(j,nstor+i) = 
     &         (coo(j,i) - vec(j) + abc(j)/2.0d0)/abc(j)
         end do
         ianz(nstor+i) = ianz(i)
         do j=1,iconn(1,i)+1
            iconn(j,nstor+i) = iconn(j,i)
         end do
         iatclr(nstor+i) = iatclr(i)

      end do

      iftyp = 6

      return
      end

      subroutine wrchd(iun,coo,ianz,iatclr,iconn,qat,ityp,
     &                 nat,inorm,nspg,icel,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      integer*2 ityp
      common /types/ iff
      parameter (mxel=100)
      character*2 elemnt
      common /elem/elemnt(mxel)
      logical ochg,dochg
      character*2 tocapf
      character*4 atname
      dimension ich(8),vec(3),abc(3)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*),
     &          qat(*),ityp(*)

      torad = datan(1.0d0) / 45.0d0
      toang = 0.52917706d0

      dochg = .false.

      if (icel.eq.1) then

         dochg = ochg(idum,ianz)
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
         if (ihasq.eq.1) dochg = .true.
         natoms = iatoms
         inorm = 0
         call cntvec(vec,coo,ianz,iatoms)
         do i=1,3
            abc(i) = 0.0d0
         end do
         do i=1,iatoms
            do j=1,3
               t = dabs(coo(j,i) - vec(j))
               if (t.gt.abc(j)) abc(j) = t
            end do
         end do
         do i=1,3
            abc(i) = 2.0d0*abc(i) + 3.0d0
         end do
      endif

      if (dochg) call fxleak(natoms,1,qat)

      if (icel.eq.1) then
         write(iun,'(38x,3f8.3)') a,b,c
         write(iun,'(21x,3f8.3,9x,i3)') 
     &     alpha/torad,beta/torad,gamma/torad,nspg
      else
         write(iun,'(38x,3f8.3)') (abc(i)*toang,i=1,3)
         write(iun,'(21x,a)') 
     &     '  90.000  90.000  90.000           1'
      endif
      write(iun,'(2i4)') natoms,inorm
      write(iun,'(8x,a)') 'molden generated cssr'

      do i=1,natoms
          do j=1,8
             ich(j) = 0
          end do
          k = 0
          do j=1,iconn(1,i)
             if (iconn(1+j,i).gt.0) then
                k = k + 1
                ich(k) = iconn(1+j,i)
             endif
          end do

          if (dochg) then
             iq = int(qat(i)*1000)
             q = dfloat(iq)/1000.0d0
          else
             q = 0.0d0
          endif

          if (iff.ne.0) then
             it = ityp(i)
          else
             it = 0
          endif

          atname = tocapf(elemnt(ianz(i)))//'  '
          if (atname(1:1).eq.' ') then
              atname(1:1) = atname(2:2)
              atname(2:2) = ' '
          endif

          if (icel.eq.1) then
             if (iff.eq.0) then
                write(iun,100)
     &          i,atname,(coo(j,i),j=1,3),(ich(j),j=1,8),q
             else
                write(iun,100)
     &          i,atname,(coo(j,i),j=1,3),(ich(j),j=1,8),q,it
             endif
          else
             if (iff.eq.0) then
                write(iun,100) i,atname,
     &          (((coo(j,i)-vec(j)+abc(j)/2.0d0)/abc(j)),j=1,3),
     &          (ich(j),j=1,8),q
             else
                write(iun,100) i,atname,
     &          (((coo(j,i)-vec(j)+abc(j)/2.0d0)/abc(j)),j=1,3),
     &          (ich(j),j=1,8),q,it
             endif
          endif
      end do

      if (icel.eq.1) call fdat(1,0,0,0,0,0)

100   format(i4,1x,a4,2x,3(f9.5,1x),8i4,f8.3,i4)

      return
      end

      subroutine wrpol(iun,i,j,k,itype,qrot)      
      implicit double precision (a-h,o-z)
      dimension qrot(*)

      kk = k
      if (itype.eq.2) kk = -kk

      write(iun,'(a,3i5,3x,3f12.5)') 'multipole     ',
     &      -i,j,kk,qrot(1)
      write(iun,100) (qrot(l+1),l=1,3)
      write(iun,100) qrot(5)
      write(iun,100) (qrot(l+7),l=1,2)
      write(iun,100) (qrot(l+10),l=1,3)

100   format(32x,3f12.5)

      return
      end

      subroutine rotloc(i,j,k,itype,qrot)      
      implicit double precision (a-h,o-z)

c     Input:   i=multipole site  I
c              j=z-axis defining atom  
c              k=x-axis definng atom
c              itype=1 ('Z-then-X') or 2 ('Bisector')
c     Output : qrot(1-13): cartesian multipoles in local axes

      parameter (mxsite=300)
      parameter (lmax=4)
      parameter (lmaxf=(lmax+1)*(lmax+1))
      character*8   ctag
      common /multip/ qmom(lmaxf,mxsite),car(3,mxsite),ctag(mxsite),
     &                nsites
      dimension qrot(*),a(3,3),b(3,3)

C Convert to cartesian multipoles

      w34 = dsqrt(3.0d0/4.0d0)

c             q=Q_00
      qrot(1) = qmom(1,i)

c             m_x=Q_1c
      qrot(2) = qmom(3,i)

c             m_y=Q_11s
      qrot(3) = qmom(4,i)

c             m_z=Q_110
      qrot(4) = qmom(2,i)

c             Q_xx=V(3/4)Q_22c - 1/2 Q_20
      qrot(5) = w34*qmom(8,i) - qmom(5,i)/2.0d0

c             Q_xy=V(3/4)Q_22s
      qrot(8) = w34*qmom(9,i)

c             Q_yy=-V(3/4)Q_22c - 1/2 Q_20
      qrot(9) = -w34*qmom(8,i) - qmom(5,i)/2.0d0

c             Q_xz=V(3/4)Q_21c
      qrot(11) = w34*qmom(6,i)

c             Q_yz=V(3/4)Q_21s
      qrot(12) = w34*qmom(7,i)

c             Q_zz=Q_20
      qrot(13) = qmom(5,i)

c             make quadrupole matrix symmetric
      qrot(6) = qrot(8)
      qrot(7) = qrot(11)
      qrot(10) = qrot(12)

c
c Get rotation matrix local to global
c
      call rtmat(i,j,k,itype,a)
c      
c Global to local is inverse:
c
      do l=1,3
         do m=1,3
            b(l,m) = a(m,l)
         end do
      end do
      
c Rotate moments:

      call rotpole(b,qrot)

      end

      subroutine fxleak(nchrg,idebug,q)      
      implicit double precision (a-h,o-z)
      dimension q(*)
c
c Fix any minute charge leaks, fix to zero
c
      totch = 0.0d0

      do i=1,nchrg
          totch = totch + q(i)
      end do

      if (idebug.eq.1) print*,'total charge =',totch

      if (dabs(totch).gt.1.0d-7) then

         totch = totch / dble(nchrg)
         totfx = 0.0d0

         do i=1,nchrg
             q(i) = q(i) - totch
             totfx = totfx + q(i)
         end do

         if (idebug.eq.1) print*,'total fixed charge =',totfx

      endif

      return
      end

      subroutine rdmsd(idebug,istat,coo,ianz,iatclr,iconn,qat,
     &                 nat,icent,ncon,nspg,ichx,nopr,ir,it,
     &                 xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxmsi=500)
      parameter (mxcon=10)
      common /athlp/  iatoms, mxnat
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      character*137 line,string
      common /curlin/ line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      integer*2 ir,it
      integer getlin
      logical gnreal,ochg
      dimension tmp(3),iatptr(mxmsi),imsi(2,4*mxmsi),qtmp(mxmsi)
      dimension a3(3),b3(3),c3(3)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*),qat(*)
      dimension ir(3,3,192),it(3,192)

      toang = 0.52917706d0
      iatoms = 0
      istat = 1
      ihasq = 0
      nspg = 0

      do i=1,mxmsi
         qtmp(i) = 0.0d0
      end do

      call nxtlin(line,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 100
      if (line(1:5).ne.'# MSI') goto 100

      todeg = 45.0d0 / datan(1.0d0)
      inewat = 0
      idm = 0
      nbnd = 0
      imodel = 0
      ihcl = 0

      do while (getlin(2).eq.1)

         ktype = nxtwrd(string,nstr,itype,rtype)
         if (ktype.eq.2) iptt = itype

         ktype = nxtwrd(string,nstr,itype,rtype)
         if (ktype.eq.1) then
            if (string(1:nstr).eq.'Atom') then
               ipt = iptt
               inewat = 0
               idm = 0
            elseif (string(1:nstr).eq.'DistanceMonitor') then
               idm = 1
            elseif (string(1:nstr).eq.'HBond') then
               idm = 1
            elseif (string(1:nstr).eq.'Bond') then
               idm = 1
            elseif (string(1:nstr).eq.'Anchor') then
               idm = 1
            elseif (string(1:nstr).eq.'Model') then
               imodel = imodel + 1
               if (imodel.gt.1) goto 10
            endif
         endif

         ktype = nxtwrd(string,nstr,itype,rtype)
         if (ktype.eq.1) then

            if (string(1:nstr).eq.'ACL') then

               ktype = nxtwrd(string,nstr,itype,rtype)
               if (ktype.eq.2) then
                  iatnmr = itype
                  if (iatnmr.eq.0) iatnmr = 1
               endif

            elseif (string(1:nstr).eq.'Charge') then

               ktype = nxtwrd(string,nstr,itype,rtype)
               if (ktype.eq.3) then
                  ihasq = 1
                  qt = rtype
               endif

            elseif (string(1:nstr).eq.'XYZ'.and.idm.eq.0) then

               if (gnreal(tmp,3,.false.)) then
                  if (iatoms.lt.mxmsi) then
                     iatoms = iatoms + 1
                     iconn(1,iatoms) = 0
                     iatptr(iatoms) = ipt
                     ianz(iatoms) = iatnmr
                     iatclr(iatoms) = 1
                     if (ihasq.eq.1) qtmp(iatoms) = qt
                     do j=1,3
                        coo(j,iatoms) = tmp(j) / toang
                     end do
                     inewat = 1
                  endif
               endif

            elseif (string(1:nstr).eq.'Atom1') then

               ktype = nxtwrd(string,nstr,itype,rtype)
               if (ktype.eq.2) then
                  iat1 = itype
                  if (getlin(2).eq.1) then

                     do i=1,3
                        ktype = nxtwrd(string,nstr,itype,rtype)
                     end do

                     if (ktype.eq.1) then
                       if (string(1:nstr).eq.'Atom2') then
                        ktype = nxtwrd(string,nstr,itype,rtype)
                        if (ktype.eq.2.and.inewat.eq.1) then
                           if (nbnd.lt.mxmsi*4) then
                              nbnd = nbnd + 1
                              imsi(1,nbnd) = iat1
                              imsi(2,nbnd) = itype
                           endif
                        endif
                       endif
                     endif

                  endif
               endif

            elseif (string(1:nstr).eq.'A3') then

               if (gnreal(tmp,3,.false.)) then
                  do j=1,3
                     a3(j) = tmp(j)
                  end do
                  a = vlen(a3)
                  ihcl = ihcl + 1
               endif
               
            elseif (string(1:nstr).eq.'B3') then

               if (gnreal(tmp,3,.false.)) then
                  do j=1,3
                     b3(j) = tmp(j)
                  end do
                  b = vlen(b3)
                  ihcl = ihcl + 1
               endif
               
            elseif (string(1:nstr).eq.'C3') then

               if (gnreal(tmp,3,.false.)) then
                  do j=1,3
                     c3(j) = tmp(j)
                  end do
                  c = vlen(c3)
                  ihcl = ihcl + 1
               endif
               
            elseif (string(1:nstr).eq.'SpaceGroup') then

               ktype = nxtwrd(string,nstr,itype,rtype)
               if (ktype.eq.2) then
                   nspg = itype
               endif

            endif
         endif
      end do

10    continue

      do i=1,nbnd
         iat1 = 0
         do j=1,iatoms
            if (iatptr(j).eq.imsi(1,i)) iat1 = j
         end do
         iat2 = 0
         do j=1,iatoms
            if (iatptr(j).eq.imsi(2,i)) iat2 = j
         end do
         if (iat1.ne.0.and.iat2.ne.0) then
            if (iconn(1,iat1).lt.mxcon) then
               iconn(1,iat1) = iconn(1,iat1) + 1
               iconn(1+iconn(1,iat1),iat1) = iat2
            endif
            if (iconn(1,iat2).lt.mxcon) then
               iconn(1,iat2) = iconn(1,iat2) + 1
               iconn(1+iconn(1,iat2),iat2) = iat1
            endif
         endif
      end do
      
c if read in charges are zero and old atoms array matches up with
c with new retain old charges

      iq = 0
      do i=1,iatoms
         if (qtmp(i).ne.0.0d0) iq = 1
      end do

      if (iq.eq.1) then
         ihasq = 1
         do i=1,iatoms
            qat(i) = qtmp(i)
         end do
      else
         if (ochg(idum,ianz)) then
            ihasq = 1
         else
            ihasq = 0
            do i=1,iatoms
               qat(i) = 0.0d0
            end do
         endif
      endif

      if (ihcl.eq.3) then

         if (nspg.eq.0) nspg = 1

         call impsc(b3,c3,csa)
         call impsc(a3,c3,csb)
         call impsc(a3,b3,csc)
         alpha = dacos(csa)*todeg
         beta = dacos(csb)*todeg
         gamma = dacos(csc)*todeg

         call prcell(nspg,a,b,c,alpha,beta,gamma)
         call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)
         call cprot(nspg,nopr,icent,ir,it,.false.)
         nat = iatoms
         ncon = 0
         natoms = nat
         nstor = mxnat-natoms
         do i=1,natoms
            do j=1,3
               tmp(j) = coo(j,i) * toang
            end do
            call crt2fr(tmp,coo(1,nstor+i),xa,ya,yb,za,zb,zc)
            ianz(nstor+i) = ianz(i)
            do j=1,iconn(1,i)+1
               iconn(j,nstor+i) = iconn(j,i)
            end do
            iatclr(nstor+i) = iatclr(i)
         end do

         istat = 2
         ichx = 1

      endif

      return

100   istat = 0
      return
      end

      subroutine wrmsd(iun,coo,ianz,iconn,qat,
     &                 nat,nspg,icel,xa,ya,yb,za,zb,zc)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxel=100)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      character*2 elemnt
      common /elem/elemnt(mxel)
      dimension tmp(3)
      dimension a3(3),b3(3),c3(3),icont(mxcon+1)
      dimension coo(3,*),ianz(*),iconn(mxcon+1,*),qat(*)

      toang = 0.52917706d0
      write(iun,'(a)') '# MSI CERIUS2 DataModel File Version 3 9'
      write(iun,'(a)') '(1 Model'

      if (icel.eq.1) then

         do i=1,3
            a3(i) = 0.0d0
            b3(i) = 0.0d0
            c3(i) = 0.0d0
         end do
         a3(1) = 1.0d0
         b3(2) = 1.0d0
         c3(3) = 1.0d0

         natoms = nat
         nstor = mxnat-natoms
         do i=1,natoms
            do j=1,iconn(1,nstor+i)+1
               iconn(j,i) = iconn(j,nstor+i)
            end do
         end do

         write(iun,'(a)') ' (A I PeriodicType 100)'
         call fr2crt(a3,xa,ya,yb,za,zb,zc)
         write(iun,'(a,3f12.6,a)') ' (A D A3 (',(a3(i),i=1,3),'))'
         call fr2crt(b3,xa,ya,yb,za,zb,zc)
         write(iun,'(a,3f12.6,a)') ' (A D B3 (',(b3(i),i=1,3),'))'
         call fr2crt(c3,xa,ya,yb,za,zb,zc)
         write(iun,'(a,3f12.6,a)') ' (A D C3 (',(c3(i),i=1,3),'))'
         write(iun,'(a,a1,i3,a,a1,a)') ' (A C SpaceGroup ',char(34),
     &        nspg,' 1',char(34),')'

      else
         natoms = iatoms
      endif

      write(iun,'(a)') ' (2 Sequence'
      write(iun,'(a)') '  (A O SequenceList (3))'
      write(iun,'(a,a1,a,a1,a)') '  (A C Label ',char(34),
     &        'NEWS',char(34),')'
      write(iun,'(a)') '  (3 Subunit'
      write(iun,'(a,a1,a,a1,a)') '   (A C SubId ',char(34),
     &        'RES1',char(34),')'
      write(iun,'(a,a1,a,a1,a)') '   (A C Type ',char(34),
     &        '/RES1',char(34),')'
      write(iun,'(a,a1,a,a1,a)') '   (A C Label ',char(34),
     &        'RES1',char(34),')'

      ntel = 3
      do i=1,natoms
         if (icel.eq.1) then
            do j=1,3
               tmp(j) = coo(j,nstor+i)
            end do
            call fr2crt(tmp,xa,ya,yb,za,zb,zc)
            ian = ianz(nstor+i)
         else
            do j=1,3
               tmp(j) = coo(j,i) * toang
            end do
            ian = ianz(i)
         endif
         write(iun,'(a,i4,a)') '   (',ntel+i,' Atom'
         if (ian.eq.1) then
            write(iun,'(a,a1,a,a2,a1,a)') '    (A C ACL ',char(34),
     &         '  0 ',elemnt(ian),char(34),')'
         else
            write(iun,'(a,a1,i3,a,a2,a1,a)') '    (A C ACL ',char(34),
     &         ian,' ',elemnt(ian),char(34),')'
         endif
         if (ihasq.eq.1) 
     &       write(iun,'(a,f12.6,a)') '    (A F Charge ',qat(i),')'
         write(iun,'(a,3f12.6,a)') '    (A D XYZ (',
     &         (tmp(j),j=1,3),'))'
         write(iun,'(a,i3,a)') '    (A I Id ',i,')'
         write(iun,'(a)') '   )'
      end do

      ntel = ntel + natoms
      do i=1,natoms
         if (icel.eq.1) then
            do j=1,iconn(1,nstor+i)+1
               icont(j) = iconn(j,nstor+i)
            end do
         else
            do j=1,iconn(1,i)+1
               icont(j) = iconn(j,i)
            end do
         endif
         do j=1,icont(1)
            ntel = ntel + 1
            write(iun,'(a,i4,a)') '   (',ntel,' Bond'
            jj = iabs(icont(1+j))
            if (jj.gt.i.and.icont(1+j).gt.0) then
               write(iun,'(a,i4,a)') '    (A O Atom1 ',3+i,')'
               write(iun,'(a,i4,a)') '    (A O Atom2 ',3+jj,')'
            endif
            write(iun,'(a)') '   )'
         end do
      end do

      write(iun,'(a)') '  )'
      write(iun,'(a)') ' )'
      write(iun,'(a)') ')'


      return
      end

      subroutine zm3rot(cm,imkeep,coo)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      dimension imkeep(3),i3(3),tr1(3)
      dimension cm(3,3),coo(3,*)

      do j=1,3
         i3(j) = j
      end do

      do i=1,3
         tr1(i) = cm(i,1) - coo(i,1)
      end do

      do i=1,iatoms
         call trcoo(tr1,coo(1,i))
      end do

      call alntwo(cm,i3,coo,iatoms,i3)

      return
      end

      subroutine rdcifd(idebug,istat,
     &                  coo,ianz,iconn,iatclr,
     &                  nat,icent,ncon,nspg,nopr,ir,it,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      integer*2 ir,it
      integer getlin
      character*137 line,str
      common /curlin/ line
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      logical chkrec
      dimension ianz(*),coo(3,*),iconn(mxcon+1,*),iatclr(*)
      dimension ir(3,3,192),it(3,192)

      call rewfil
      istat = 2

      toang = 0.52917706d0

      nspg = 0
      ilatt = 0

      icent = 0
      nat = 0
      iloop = -1
      ilab = 0
      ifrx = 0
      ifry = 0
      ifrz = 0

      call recsym(nopr,ir,it,'SYMM X,Y,Z',1)

      call search(line,'data_',istat1)
      if (istat1.eq.0) goto 100

      do while (getlin(3).eq.1)

         if (icdex(line,'_symmetry_Int_Tables_number').ne.0) then

             ktype = nxtwrd(str,nstr,itype,rtype)
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.2) then
                nspg = iabs(itype)
             else
                goto 100
             endif
        
         elseif (icdex(line,'_cell_length_a').ne.0) then

             ktype = nxtwrd(str,nstr,itype,rtype)
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.3) then
                a = rtype
             else
                goto 100
             endif

         elseif (icdex(line,'_cell_length_b').ne.0) then

             ktype = nxtwrd(str,nstr,itype,rtype)
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.3) then
                b = rtype
             else
                goto 100
             endif

         elseif (icdex(line,'_cell_length_c').ne.0) then

             ktype = nxtwrd(str,nstr,itype,rtype)
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.3) then
                c = rtype
             else
                goto 100
             endif

         elseif (icdex(line,'_cell_angle_alpha').ne.0) then

             ktype = nxtwrd(str,nstr,itype,rtype)
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.3) then
                alpha = rtype
             elseif (ktype.eq.2) then
                alpha = dble(itype)
             else
                goto 100
             endif

         elseif (icdex(line,'_cell_angle_beta').ne.0) then

             ktype = nxtwrd(str,nstr,itype,rtype)
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.3) then
                beta = rtype
             elseif (ktype.eq.2) then
                beta = dble(itype)
             else
                goto 100
             endif

         elseif (icdex(line,'_cell_angle_gamma').ne.0) then

             ktype = nxtwrd(str,nstr,itype,rtype)
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.3) then
                gamma = rtype
             elseif (ktype.eq.2) then
                gamma = dble(itype)
             else
                goto 100
             endif

         elseif (icdex(line,' _symmetry_equiv_pos_as_xyz').ne.0) then

             ilatt = 1

             do while (.true.)
                if (getlin(2).eq.0) goto 90
                if (linlen(line).eq.0) goto 90
                call tocap(line,len(line))
                call recsym(nopr,ir,it,'SYMM '//line,ilatt)
             end do

90          continue

         elseif (icdex(line,'loop_').ne.0) then
             iloop = 0
             ilab = 0
             ifrx = 0
             ifry = 0
             ifrz = 0
         elseif (iloop.ge.0) then
             iloop = iloop + 1
c             if (icdex(line,'_atom_site_label').ne.0)   ilab = iloop
             if (icdex(line,'_atom_site_type_symbol').ne.0) ilab = iloop
             if (icdex(line,'_atom_site_fract_x').ne.0) ifrx = iloop
             if (icdex(line,'_atom_site_fract_y').ne.0) ifry = iloop
             if (icdex(line,'_atom_site_fract_z').ne.0) ifrz = iloop
             if (index(line,'#').ne.0) goto 40
             if (linlen(line).eq.0.and.nat.ne.0) goto 50
             if (index(line,'_').eq.0.and.line(1:1).ne.'#') then
                if (ilab.ne.0.and.
     &              ifrx.ne.0.and.ifry.ne.0.and.ifrz.ne.0) then
c proces data line
                  nat = nat + 1

                  do i=1,ifrz

                     ktype = nxtwrx(str,nstr,itype,rtype)

                     if (i.eq.ilab) then
                        if (ktype.eq.1) then
                           call parlab(str,nstr,iat)
                           if (iat.gt.0) then
                              ianz(nat) = iat
                           else
                              goto 100
                           endif
                        else
                           goto 100
                        endif
                     elseif (i.eq.ifrx) then
                        if (ktype.eq.2) then
                           coo(1,nat) = dble(itype)
                        elseif (ktype.eq.3) then
                           coo(1,nat) = rtype
                        else
                           goto 100
                        endif
                     elseif (i.eq.ifry) then
                        if (ktype.eq.2) then
                           coo(2,nat) = dble(itype)
                        elseif (ktype.eq.3) then
                           coo(2,nat) = rtype
                        else
                           goto 100
                        endif
                     elseif (i.eq.ifrz) then
                        if (ktype.eq.2) then
                           coo(3,nat) = dble(itype)
                        elseif (ktype.eq.3) then
                           coo(3,nat) = rtype
                        else
                           goto 100
                        endif
                     endif
                  end do
                endif
             endif

         endif

40       continue

      end do

50    continue

      if (nat.eq.0) goto 100

      call prcell(nspg,a,b,c,alpha,beta,gamma)
      call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)

      call cprot(nspg,nopr,icent,ir,it,.false.)
      if (idebug.eq.1) call prop(nopr,ir,it)

      nstor = mxnat-nat

      do i=1,nat
         do j=1,3
            coo(j,nstor+i) = coo(j,i)
         end do
         iconn(1,i) = 0
         call fr2crt(coo(1,i),xa,ya,yb,za,zb,zc)
         do j=1,3
             coo(j,i) = coo(j,i)/toang
         end do
      end do

      iatoms = nat

      call doconn
      call dohcon(0)

      do i=1,nat
         ianz(nstor+i) = ianz(i)
         do j=1,iconn(1,i)+1
            iconn(j,nstor+i) = iconn(j,i)
         end do
         iatclr(nstor+i) = 1
      end do

      ncon = 1

      if (idebug.eq.1) call prop(nopr,ir,it)

      return

100   istat = 0
      if (idebug.eq.1) print*,'rdcif: Error line: ',line
      return
      end

      
      subroutine parlab(str,nstr,iat)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      character*137 str
      character*2 catom,ctmp,tolowf,iel
      dimension iel(100)
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
     & 'bk','cf','x '/


      iat = 0
      if (nstr.eq.1) then
         ctmp(1:1) = str(1:1)
         ctmp(2:2) = ' '
      else
         ctmp = str(1:2)
      endif

      catom = tolowf(ctmp)
      ic = ichar(catom(2:2))
      if (ic.lt.97.or.ic.gt.122) catom(2:2) = ' '

      do j=1,100
         if (catom .eq. iel(j)) iat = j - 1
      end do

      return
      end

      logical function odupl(it,coo,xa,ya,yb,za,zb,zc)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      logical odup
      common /pbc/ abc(3),ibox,icell,igfmap
      dimension coo(3,*),crt1(3),crt2(3),fr1(3),fr2(3)

      odupl = .true.

      if (icell.ne.1) return

      toang = 0.52917706d0
      tl = 1.0d-4

      if (it.gt.1) then

         do j=1,3
            crt1(j) = coo(j,it) * toang
         end do

         call crt2fr(crt1,fr1,xa,ya,yb,za,zb,zc)

         odupl = .false.

         do i=1,it-1

            do j=1,3
               crt2(j) = coo(j,i) * toang
            end do

            call crt2fr(crt2,fr2,xa,ya,yb,za,zb,zc)

            odup = .true.
            do j=1,3
               dd = dmod(dabs(fr2(j)-fr1(j)),1.0d0)
               odup = (odup.and.dd.lt.tl)
            end do

            if (odup) return

         end do
      endif

      odupl = .true.

      return
      end

      subroutine addc(coo,ianz,iconn,iatclr,iatoms,xa,ya,yb,za,zb,zc)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /cllhlp/ ico(3,8),icn(4,8),mcol(32)
      common /cllab/  iclon,iclpnt(4)
      common /surf/   natorg,noscnd
      dimension ianz(*),coo(3,*),iconn(mxcon+1,*),iatclr(*)

      toang = 0.52917706d0
      natorg = iatoms
      iclon = 1
      iclpnt(1) = iatoms + 1
      iclpnt(2) = iatoms + 2
      iclpnt(3) = iatoms + 3
      iclpnt(4) = iatoms + 4

      do i=1,8
         l = iatoms + i
         do j=1,3
            coo(j,l) = dfloat(ico(j,i))
         end do
         call fr2crt(coo(1,l),xa,ya,yb,za,zb,zc)
         do j=1,3
            coo(j,l) = coo(j,l)/toang
         end do
         iconn(1,l) = icn(1,i)
         do j=2,4
            iconn(j,l) = iatoms + icn(j,i)
         end do
         ianz(l) = 100
         iatclr(l) = 11
      end do

      iatoms = iatoms + 8

      return
      end

      subroutine updc(coo,xa,ya,yb,za,zb,zc)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /cllhlp/ ico(3,8),icn(4,8),mcol(32)
      common /surf/   natorg,noscnd
      dimension coo(3,*)

      toang = 0.52917706d0

      do i=1,8
         l = natorg + i
         do j=1,3
            coo(j,l) = dfloat(ico(j,i))
         end do
         call fr2crt(coo(1,l),xa,ya,yb,za,zb,zc)
         do j=1,3
            coo(j,l) = coo(j,l)/toang
         end do
      end do

      return
      end

