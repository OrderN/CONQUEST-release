      subroutine rdgdud(idebug,ibefo,istatio,irtype,istats,
     &                 focc,focb,nocc,nocb,ncols,ncolb,coo,ianz)

c THIS IS REALLY rdgaus

      implicit double precision ( a-h,o-z)
      parameter (numat1=20000)
      parameter (mxonh=100)
      common /athlp/ iatoms, mxnat
      common /qmchar/ qch(numat1),ihasesp
      parameter (numatm=2000)
      parameter (mxel=100)
      character*137 line,lint
      character*2 els
      common /coord / xyz(3,numatm)
      common /pseudo /ipseud,ivale(numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /orbhlp/ mxorb,iuhf,ispd
      common /gauori/ nzm,nso,nio,nzo,ioropt,ifor,
     &                ixyz98,iopr,isymm,irc,imp2,icntp,itd
      common /gauver/ ivers
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /oniomh/ ioni,nion,ionih(mxonh),natonh(mxonh),
     &                iomap(numatm),icntat(mxonh),fct(mxonh),
     &                xyzi(3,mxonh)
      common /nmr/    shlnuc(numatm),ihsnmr
      common /uvspec/ ihasex
      logical rohf
      dimension focc(*),focb(*),v1(3)
      dimension coo(3,*),ianz(*)

      istats = 1
      irtype = 0
      istatio = 0
      isymm = 1
      rohf = .false.
      ig94 = 0
      ixyz98 = 0
      iopr = 0
      irc = 0
      imp2 = 0
      itd = 0
      icntp = 0
      ihasex = 0

      call search(line,' ***********************',istat)
      read(iun2,'(a)') line
      if (line(11:12).eq.'DV') then
          ivers = 2003
      else
          read(line,'(10x,i2)') ivers
          if (ivers.eq.3) ivers = 2003
          if (ivers.eq.9) ivers = 2009
          if (ivers.ge.94) ig94 = 1
      endif
      if (ivers.ge.98) then
          ixyz98 = 1
          call search(line,'will use Cartesian coordinates',istat)
          if (istat.ne.0) ixyz98 = 2
          rewind iun2
      endif

      if (ig94.eq.0) then
         call search(line,'Restricted open shell SCF:',istat)
         if (istat.ne.0) rohf = .true.
         rewind iun2
      endif
      
      call searchd(line,'Frequencies --',
     &                  'Excitation energies and oscillator strengths:',
     &                  istat)
      if (istat.ne.0) then
          if (index(line,'Frequencies --').ne.0) irtype = 4
          if (index(line,'Excitation').ne.0) then
              ihasex = 1
              call gttrns(ist)
          endif
      endif
      rewind iun2

      call search(line,'Symmetry turned off',istat)
      if (istat.ne.0) isymm = 0
      rewind iun2

      call search(line,'IRC-IRC-IRC',istat)
      if (istat.ne.0) then
         irc = 1
         call search(line,'Integration scheme',istat)
         if (istat.ne.0) then
             if (icdex(line,'HPC').ne.0) irc = 2
         endif
      endif
      rewind iun2

      call searchd(line,'Computing MP2 derivatives',
     &                  'Computing MP2/KS-MP2 derivatives',istat)
      if (istat.ne.0) imp2 = 1
      rewind iun2

      call search(line,'Computing CIS/TD-HF/TD-KS derivatives',istat)
      if (istat.ne.0) itd = 1
      rewind iun2


      call search(line,'Counterpoise: corrected',istat)
      if (istat.ne.0) icntp = 1
      rewind iun2

      call searchd(line,'ONIOM: saving','ONIOM: Cut',istat)
      if (istat.ne.0) then
         ioni = 1
         nion = 0
         do while (index(line,'Cut').ne.0)
            if (nion.lt.mxonh) then
               nion = nion + 1
               id = index(line,'/') + 1
               els = line(id:id+1)
               if (els(2:2).eq.' ') then
                  els = ' '//els(1:1)
               endif
               do i=1,99
                  if (els.eq.elemnt(i)) natonh(nion) = i
               end do
               if (line(id+8:id+8).eq.' ') then
                  read(line(id+2:id+38),'(i6,7x,i6,8x,f10.6)') 
     &              ionih(nion),icntat(nion),fct(nion)
               else if (line(id+7:id+7).eq.' ') then
                  read(line(id+2:id+35),'(i5,7x,i5,8x,f9.6)') 
     &              ionih(nion),icntat(nion),fct(nion)
               else
                  print*,'error reading ONIOM: Cut atoms'
               endif
            endif
            call readel(line,1)
         end do
         if (idebug.eq.1) then
            print*,'Oniom boundary atoms:'
            do i=1,nion
               print*,natonh(i),' ',ionih(i)
            end do
         endif
      else 
         ioni = 0
      endif
      
      if (ioni.eq.1) then
         call searchd(line,'high level on model system',
     &                     'generating new system at layer 1',istat)
         if (istat.eq.0) then
            call inferr('no high level on model system found!',1)
            goto 1000
         endif
      endif
c     
c     look for number of orbitals and electrons
c
      call search(line,'primitive gaussians',istat)
      if (istat.eq.0) then
         call inferr('no primitive gaussian found!',1)
         goto 1000
      endif
      if (index(line,'primitive gaussians').eq.30) then
         read(line,'(1x,i3)',err=1000,end=1000) norbs
      else
         read(line,'(1x,i5)',err=1000,end=1000) norbs
      endif
      if (norbs.gt.mxorb) 
     &   call inferr('Exceeding MaxNum of Orbitals!',1)
      call search(line,'alpha electrons',istat)
      if (istat.eq.0) then
         call inferr('no alpha electrons found!',1)
          goto 1000
      endif
      if (index(line,'alpha electrons').eq.6) then
         read(line,'(1x,i3,21x,i4)',err=1000,end=1000) neleca,nelecb
      else
         read(line,'(1x,i5,20x,i5)',err=1000,end=1000) neleca,nelecb
      endif
      nelecs = neleca + nelecb
      iorbs = 0
      call search(line,'One-electron integrals computed using',istat)
      if (istat.eq.0) then
         rewind iun2
         call search(line,'alpha electrons',istat)
      endif
      do i=1,30
         read(iun2,'(a)',end=999) line
         if (index(line,'RelInt:').eq.0) then
            if (index(line,'NBsUse=').ne.0) then
               read(line,'(9x,i5)',err=1000,end=1000) iorbs
            endif
         endif
      end do
999   continue
      if (ig94.eq.1) then
         ncols  = norbs
         ncolb  = norbs
         if (iorbs.ne.0) then
            ncols = iorbs
            ncolb = iorbs
         endif
      else
         ncols  = max0(neleca+5,nelecb+5)
         ncols  = min0(ncols,norbs)
         if (rohf) ncols = norbs
         ncolb  = ncols
      endif
c
      do i=1,min0(mxorb,norbs)
         focc(i) = 0.0d0
      end do
      do i=1,min0(mxorb,neleca)
         focc(i) = 1.0d0
      end do

      rewind iun2
      if (ioni.eq.1) then
         call searchd(line,'high level on model system',
     &                     'generating new system at layer 1',istat)
         if (istat.eq.0) then
            call inferr('no high level on model system found!',1)
            goto 1000
         endif
      endif
c      call searchd(line,'UHF open shell SCF','E(UHF)',istat)
      call search(line,'Beta Molecular Orbital Coefficients',istat)
      if (istat.ne.0) then
         iuhf=1
         nocc=neleca
         nocb=nelecb
         do i=1,min0(mxorb,norbs)
            focb(i) = 0.0d0
         end do
         do i=1,min0(mxorb,nelecb)
            focb(i) = 1.0d0
         end do
      else
         nocc=max0(neleca,nelecb)
         do i=1,min0(mxorb,nelecb)
            focc(i) = focc(i) + 1.0d0
         end do
      endif
c
      call search(line,'-- Stationary point found',istat)
      rewind iun2
      ihssta = istat

c      if (ivers.ge.2009) then
c          ibefo = 1
c          write(iun3,*)
c     &          'Warning G09: using orbitals input geometry!'
c          write(iun3,*) ' '
c      endif

      if (istat.eq.0.or.ibefo.eq.1.or.(istat.eq.1.and.irtype.eq.4)) 
     &then
          write(iun3,*)'First standard orientation encountered used'
          write(iun3,*)'for density calculations'
          call search(line,'Z-MATRIX (ANGSTROMS AND DEGREES)',istat)
          if (istat.eq.0) then
              call haszm(.false.)
          else
              call readel(line,2)
              if (irtype.eq.4) then
                 call convzmat(coo,ianz,iatoms,1,0,1)
              else
                 call convzmat(coo,ianz,iatoms,1,1,1)
              endif
          endif
          call rdcor(idebug,istat)
          if (istat.eq.0) goto 1000
          if (idebug.eq.1) 
     &    call inferr('Succesfully read coordinates',0)
      else
          istatio = 1

          write(iun3,*)'Coordinates of stationary point used'
          write(iun3,*)'for density calculations'

10        call rdcor(idebug,istat)
          if (istat.eq.0) goto 20
          if (idebug.eq.1) 
     &    call inferr('found coordinates',0)
          goto 10
20        continue
      endif
c
      rewind iun2 
      call search(line,'Charges from ESP fit',istat)
      if (ihssta.eq.1) call search(line,'Charges from ESP fit',istat)
      if (istat.eq.1) then
          call readel(line,2)
          do i=1,natoms
              read(iun2,'(a)',err=1000,end=1000) line
              if (linlen(line).eq.18) then
                 read(line,'(7x,f11.6)',err=1000,end=1000) qch(i)
              else if (linlen(line).eq.21) then
                 read(line,'(11x,f11.6)',err=1000,end=1000) qch(i)
              endif
          end do
          ihasesp = 1
      endif

      rewind iun2
      call fndor(idebug)

      if (ixyz98.gt.0.and.natoms.gt.50) then

c check for IOP(2/11=1)

         rewind iun2
         call search(line,' 2/',istat)
         if (istat.eq.1) then
            if (ichar(line(1:1)).eq.32.and.
     &          ichar(line(2:2)).eq.50.and.
     &          ichar(line(3:3)).eq.47) then
               if (index(line,'11=1').ne.0) then
                  iopr = 1
               endif
            endif
         endif
      endif

      if (ioni.eq.1) then

c coordinates are read, calculate the coordinates of the oniom link 
c atoms

         do j=1,nion
            do k=1,3
                v1(k) = xyz(k,ionih(j))-xyz(k,icntat(j))
            end do
            do k=1,3
               xyzi(k,j) = xyz(k,icntat(j)) + 
     &          fct(j)*v1(k)
            end do
         end do
         call xyzcoo(1,0,0)

      endif

      rewind iun2
      if (ioni.eq.1) then
         call searchd(line,'high level on model system',
     &                     'generating new system at layer 1',istat)
         if (istat.eq.0) then
            call inferr('no high level on model system found!',1)
            goto 1000
         endif
      endif

c if ioni=1, rdbas transforms the xyz array to those of H layer only

      call rdbas(idebug,0,istat)
      call norml

      if (istat.eq.0) goto 1000
      if (idebug.eq.1) 
     &    call inferr('Succesfully read basis-set',0)
      if (idebug.eq.1) call basprt(iun3,.false.,.false.)
      if (norbs.gt.mxorb) goto 1000

      call search(line,'Pseudopotential Parameters',istat)
      if (istat.eq.1) then
         ipseud = 1
         call readel(line,4)
         do i=1,natoms
            read(iun2,'(a)') line
            read(iun2,'(a)') lint
            if (icdex(lint,'No pseudo').ne.0) then
               read(line,'(14x,i3)') ivale(i) 
            else
               read(line,'(26x,i3)') ivale(i) 
               do while(.true.)
                  call readel(line,1)
                  if (line(5:5).ne.' ') then
                     backspace iun2
                     goto 100
                  endif
               end do
            endif
100         continue
         end do
      endif

      natbck = natoms
      if (ioni.eq.1) then

c     rdbas sets the array iomap, when ioni=1

         natmp = 0
         do i=1,numatm
            if (iomap(i).ne.0) natmp = natmp + 1
         end do
         natoms = natmp
      endif

      rewind iun2
      if (ioni.eq.1.and.istatio.eq.1) then
         call search(line,'-- Stationary point found',istat)
         if (istat.eq.0) then
            call inferr('no Stationary point found!',1)
            goto 1000
         endif
      endif

      call nmrshl
      if (ihsnmr.eq.2) call nmrcpl(idebug)
      rewind iun2

      if (ioni.eq.1) then
         call searchd(line,'high level on model system',
     &                     'generating new system at layer 1',istat)
         if (istat.eq.0) then
            call inferr('error reading vectors: '//
     &      'no high level on model system found for Stationary point!'
     &      ,1)
            call inferr(
     &   'no vectors for high level on model system '//
     &   'of stationary point: try molden -b',1)
            goto 1000
         endif
      endif
      if (istatio.eq.1) then
         call rdvect(idebug,ig94,istat)
         istato = istat

         if (ioni.ne.1) then

            if (ivers.lt.2009) then

               call rdvect(idebug,ig94,istat)
               if (istat.eq.0.and.istato.eq.0) goto 1000
               
            endif

         else if (ioni.eq.1.and.istato.eq.0) then

            call inferr(
     &   'no vectors for high level on model system '//
     &   'of stationary point: use molden -b',1)
            goto 1000

         endif

      else

         call rdvect(idebug,ig94,istat)
         if (istat.eq.0) goto 1000

      endif

      if (idebug.eq.1) call inferr('Succesfully read vectors',0)


      if (idebug.eq.1) 
     &         call inferr('Succesfully read Gaussian output',0)

      return

1000  if (idebug.eq.1)  then
          call inferr('Error reading Gaussian output!',1)
          print*,line
      endif

      istats = 0
      return

      end

c     ====================================================================      
      subroutine cubtst(iun,ijag)
      implicit double precision ( a-h,o-z)
      character*137 line,str
      common /curlin/ line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /grdhlp/ mx3d,mx3d2
      integer currec,bigend,buflen
      common /rdrec/ iund,currec
      common /isend/ bigend,nlen
      character keywrd*320, keyori*320
      common /keywrd/ keywrd,keyori
      integer getlin
      character*4 buff(40)
      integer*2 buf2(256)
      integer*2 idum

      if (ijag.eq.1) then
c Jaguar cube
         iuntmp = iun2
         iun2 = iun

         do while(getlin(1).eq.1)
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.1) then
               if (nstr.ge.4) then
                  if (icdex(str,'npts').ne.0) then
                      ktype = nxtwrd(str,nstr,itype,rtype)
                      if (ktype.eq.2) npts1 = itype
                      ktype = nxtwrd(str,nstr,itype,rtype)
                      if (ktype.eq.2) npts2 = itype
                      ktype = nxtwrd(str,nstr,itype,rtype)
                      if (ktype.eq.2) npts3 = itype
                  endif
               endif
            endif
         end do

      elseif (ijag.eq.2) then

         iund = 10
         bigend = 1
         nlen = 1
         currec = 0
         buflen = 40

10       open(unit=iund,file=keywrd,access="direct",
     &        recl=nlen)

         
         call getrec(buff,buflen,1,ierr)
         if (ierr.eq.1) then
            bigend = 0
            currec = 0
            rewind(iund)
            call getrec(buff,buflen,1,ierr)
            if (ierr.eq.1.and.nlen.eq.1) then
               currec = 0
               bigend = 1
               nlen = 4
               close(iund)
               goto 10
            endif
         endif

         if (ierr.eq.1) then
            ijag = -1
            return
         endif

         call byter(buff(26),npts1)
         call byter(buff(27),npts2)
         call byter(buff(28),npts3)

         rewind(iund)
         currec = 0

      elseif (ijag.eq.3) then

c O map (DSN6)

         iund = 10
         bigend = 1
         nlen = 1
         currec = -1

12       open(unit=iund,file=keywrd(1:linlen(keywrd)),access="direct",
     &        recl=nlen)

         call getrc2(buf2,ierr)
         if (ierr.eq.1) then

            bigend = 0
            currec = -1
            call getrc2(buf2,ierr)
            if (ierr.eq.1.and.nlen.eq.1) then

               currec = -1
               bigend = 1
               nlen = 2
               close(iund)
               goto 12

            endif

         endif

         if (ierr.eq.1) then
            ijag = -1
            return
         endif

         npts1 = int(buf2(7))
         npts2 = int(buf2(8))
         npts3 = int(buf2(9))

c         print*,'npts ',npts1,npts2,npts3
         currec = -1

      else

         read(iun,'(a)') line
         read(iun,'(a)') line
         read(iun,'(a)') line

         read(iun,'(i5)') npts1
         read(iun,'(i5)') npts2
         read(iun,'(i5)') npts3

      endif

      nptsmx = npts1
      if (npts2.gt.nptsmx) nptsmx = npts2
      if (npts3.gt.nptsmx) nptsmx = npts3

      if (nptsmx.gt.mx3d) then
         call allgrd(nptsmx)
      endif

      if (ijag.eq.1) then
         rewind iun2
         iun2 = iuntmp
      else
         rewind iun
      endif

      return

100   print*,'error reading grid file'
      return
      end
      
c     ====================================================================
      subroutine rdcubd(npts1,npts2,npts3,iposng,ipsi,istat,iun,
     &                  idebug,denn,pmnn)
      implicit double precision ( a-h,o-z)

c THIS IS REALLY rdcube

      parameter (numatm=2000)
      parameter (mxel=100)
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /grdhlp/ mx3d,mx3d2
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /eulx/   ca,cb,sa,sb,cc,sc
      logical oclose,osshad,ofill,ofrst
      common /spacom/ grey,iscol,oclose,osshad,ofill,ofrst
      common /zproj/  zval(numatm),nindx(numatm),nconn(numatm)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      character*137 line,lstr,str
      common /curlin/ line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      logical gnreal,oampac,ojag
      integer getlin
      dimension orig(3),tmp(3),dtmp(6)
      dimension denn(*),pmnn(*)

      ojag = .false.
      if (istat.eq.1) ojag = .true.

      istat = 1
      itype = 0
      toang  = 0.52917706d0
      rneg = 0.0d0

      if (ojag) then

c Jaguar cube file

         iuntmp = iun2
         iun2 = iun
         do while(.true.)
            if (getlin(1).eq.1) then
               if (icdex(line,'&end').ne.0) goto 10
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.1) then
                  if (nstr.eq.4) then
                     if (icdex(str,'npts').ne.0) then
                        ktype = nxtwrd(str,nstr,itype,rtype)
                        if (ktype.eq.2) npts1 = itype
                        ktype = nxtwrd(str,nstr,itype,rtype)
                        if (ktype.eq.2) npts2 = itype
                        ktype = nxtwrd(str,nstr,itype,rtype)
                        if (ktype.eq.2) npts3 = itype
                     endif
                  elseif (nstr.eq.6) then
                     if (icdex(str,'orig').ne.0) then
                        if (.not.gnreal(orig,3,.false.)) goto 100
                     endif
                  elseif (nstr.eq.7) then
                     if (icdex(str,'extent').ne.0) then
                        if (str(7:7).eq.'x'.or.str(7:7).eq.'X') then
                           if (.not.gnreal(v1,3,.false.)) goto 100
                        endif
                        if (str(7:7).eq.'y'.or.str(7:7).eq.'Y') then
                           if (.not.gnreal(v2,3,.false.)) goto 100
                        endif
                        if (str(7:7).eq.'z'.or.str(7:7).eq.'Z') then
                           if (.not.gnreal(tmp,3,.false.)) goto 100
                        endif
                     endif
                  endif
               endif
            endif 
         end do
10       natoms = 0
         nat(1) = 99

      else

c Gaussian Cube

         read(iun,'(a)') line
         read(iun,'(a)') line
         if (index(line,'Density').ne.0) itype = 0
         if (index(line,'MO coefficients').ne.0) itype = 1
         ipsi = itype

         read(iun,'(a)') line
         oampac = (linlen(line).eq.44)
         if (oampac) then
            read(line,'(i5,3f13.6)') natoms,(orig(i),i=1,3)
         else
            read(line,'(i5,3f12.6)') natoms,(orig(i),i=1,3)
         endif
         if (natoms.lt.0) itype = 1
         natoms = iabs(natoms)

         if (oampac) then
            read(iun,'(i5,3f13.6)') npts1,(v1(i),i=1,3)
            read(iun,'(i5,3f13.6)') npts2,(v2(i),i=1,3)
            read(iun,'(i5,3f13.6)') npts3,(tmp(i),i=1,3)
         else
            read(iun,'(i5,3f12.6)') npts1,(v1(i),i=1,3)
            read(iun,'(i5,3f12.6)') npts2,(v2(i),i=1,3)
            read(iun,'(i5,3f12.6)') npts3,(tmp(i),i=1,3)
         endif

      endif

      if (npts1.gt.mx3d.or.npts2.gt.mx3d.or.npts3.gt.mx3d) then
          lstr = 'dimension greater than maximum !'
          print*,'npts1 ',npts1,' npts2 ',npts2,' npts3 ',npts3,
     &           ' maxdim ',mx3d 
          goto 100
      endif

      if (ojag) then

         r(1) = vlen(v1)
         r(2) = vlen(v2)
         r(3) = vlen(tmp)

         px = orig(1) + (v1(1) + v2(1) + tmp(1))/2.0d0
         py = orig(2) + (v1(2) + v2(2) + tmp(2))/2.0d0
         pz = orig(3) + (v1(3) + v2(3) + tmp(3))/2.0d0

      else

         r(1) = vlen(v1)*(npts1-1)
         r(2) = vlen(v2)*(npts2-1)
         r(3) = vlen(tmp)*(npts3-1)

         px = orig(1) + ((npts1-1)*v1(1) + (npts2-1)*v2(1) + 
     &                   (npts3-1)*tmp(1))/2.0d0
         py = orig(2) + ((npts1-1)*v1(2) + (npts2-1)*v2(2) + 
     &                   (npts3-1)*tmp(2))/2.0d0
         pz = orig(3) + ((npts1-1)*v1(3) + (npts2-1)*v2(3) + 
     &                   (npts3-1)*tmp(3))/2.0d0

      endif

      call vsc1(v1,1.0d0,1.0d-4)
      call vsc1(v2,1.0d0,1.0d-4)
      call vsc1(tmp,1.0d0,1.0d-4)

      cx = tmp(1)
      cy = tmp(2)
      cz = tmp(3)

      if (.not.ojag) then
         do i=1,natoms
            if (oampac) then
               read(iun,'(i5,4f13.6)') nat(i),tdum,(xyz(j,i),j=1,3)
            else
               read(iun,'(i5,4f12.6)') nat(i),tdum,(xyz(j,i),j=1,3)
            endif
         end do
      endif

      if (idebug.eq.1) then
         print*,'r=',(r(i),i=1,3)
         print*,'px py pz ',px,py,pz
         print*,'cx cy cz ',cx,cy,cz
         print*,'v1=',(v1(i),i=1,3)
         print*,'v2=',(v2(i),i=1,3)
         if (ojag) then
            do i=1,natoms
               print*, nat(i),(xyz(j,i),j=1,3)
            end do
         endif
      endif

      call crprod(v1,v2,dtmp)
      call impsc(dtmp,tmp,cosb)
      if (1.0d0-dabs(cosb).gt.0.001d0) then
          lstr = 'Only rectangular grids supported !'
          print*,'ERROR: Only rectangular grids supported ',cosb
          goto 100
      endif

      if (.not.ojag.and.itype.eq.1) then
         read(iun,'(i5)') no
         if (no.gt.1) then
            lstr = 'ONLY single orbital cube files are supported !'
            goto 100
         endif
         nlines = (no + 1) / 10
         if (no+1 - nlines*10.gt.0) nlines = nlines + 1
         do i=1,nlines-1
            read(iun,'(a)') line
         end do
      endif

      ij = 0
      iposng = 0
      
      if (ojag) then

c Jaguar cube file

         do i=1,npts1
            do j=1,npts2
               ij = ij + 1
               do k=1,npts3

                 if (getlin(0).eq.1) then
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.eq.3) then
                      if (rtype.lt.0.0d0) iposng = 1
                      denn((k-1)*mx3d2 + ij) = rtype 
                   endif
                 else
                   lstr = 'Not enough lines'
                   goto 100
                 endif

               end do
            end do
         end do

      else

c Gaussian Cube

         do i=1,npts1
            do j=1,npts2
               ij = ij + 1
               kk = 0
               do while (kk.lt.npts3)
                  n = 6
                  if (npts3-kk.lt.n) n = npts3-kk
                  read(iun,'(6E13.5)') (dtmp(k),k=1,n)
                  do k=1,n
                     denn((npts3-(kk+k))*mx3d2 + ij) = dtmp(k)
                     if (dtmp(k).lt.rneg) rneg = dtmp(k)
                  end do
                  kk = kk + 6
               end do
            end do
         end do
      endif

c allow tolerance for density, since gaussian density cubes often
c have small negative values

      if (dabs(rneg).gt.1.0d08) iposng = 1

      iplat = 0

      ca = v2(2)
      sa = -v2(1)
      sb = -v1(3)
      if (dabs(ca).gt.0.001d0) then
         cb = v1(1) / ca
      else
         cb = v1(2) / sa
      endif

      cc = 1.0d0
      sc = 0.0d0


      call parrat
      call proato

      do i=1,npts3
          pmnn(i) = 100000.d0
      end do

      do i=1,natoms
         nconn(i) = 0
         do j=1,natoms
           if (i.ne.j) then
            dmaxsq = (vdwr(nat(i)) + vdwr(nat(j)))**2
            dijsq = ((xyz(1,i)-xyz(1,j))*toang)**2
     &             +((xyz(2,i)-xyz(2,j))*toang)**2
     &             +((xyz(3,i)-xyz(3,j))*toang)**2
            if (dijsq.lt.dmaxsq) then
                nconn(i) = nconn(i) + 1
            endif
           endif
         end do
      end do

      call xyzcoo(1,0,0)

      ofrst = .false.

      if (ojag) iun2 = iuntmp

      lstr = 'found cube file'
      call inferr(lstr,0)
      return

100   continue
110   istat = 0
      if (ojag) iun2 = iuntmp
      call inferr(lstr,0)

      return
      end

c     ====================================================================      
      subroutine rdvasd(npts1,npts2,npts3,iposng,istat,lenf,
     &                  idocub,idebug,denn,pmnn,bucket,
     &                  coo,ianz,iatclr,iconn,
     &                  nnat,norg,icent,inorm,ncon,nspg,ichx,
     &                  nopr,ir,it,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision ( a-h,o-z)

c THIS IS REALLY rdvasp

      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (mxcon=10)
      parameter (lnbuck=10)
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /grdhlp/ mx3d,mx3d2
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /eulx/   ca,cb,sa,sb,cc,sc
      integer*2 ir,it
      logical oclose,osshad,ofill,ofrst
      common /spacom/ grey,iscol,oclose,osshad,ofill,ofrst
      common /zproj/  zval(numatm),nindx(numatm),nconn(numatm)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /athlp/ iatoms, mxnat
      character*2 elemnt
      common /elem/elemnt(mxel)

      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line,lstr
      logical opfil,gnreal
      integer getlin
      character*137 str
      character*2 tocapf
      character keywrd*320, keyori*320
      common /keywrd/ keywrd,keyori
      dimension tmp(3),dtmp(10),xyzt(3)
      dimension ieltyp(20),ielnum(20)
      dimension denn(*),pmnn(*),bucket(*)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*),
     &          ir(3,3,192),it(3,192)

      logical felem
      character*137 strtmp
      common /vasp/nheadx,natx
      
      istat = 1
      itype = 0
      toang  = 0.52917706d0
      todeg = 45.0d0 / datan(1.0d0)
      lstr = ' '

      if (idocub.eq.1) then
         iuntmp = iun2
         iun2 = 21

         if (.not.opfil(21,keywrd(1:lenf),lenf,1,1,0)) then
            lstr = keywrd(1:lenf)//' non-existent'
            goto 110
         endif
      endif

c     title

      natt = 0

      if (idocub.eq.1) then
         read(iun2,'(a)') line
      else

         felem = .false.
         if (getlin(0).eq.1) then
            do while(.true.)
               ktype = nxtwrd(str,nstr,itype,rtype)

               if (ktype.eq.1.and.(nstr.eq.1.or.nstr.eq.2)) then
                  if (nstr.eq.1) then
                     str(2:2) = str(1:1)
                     str(1:1) = ' '
                  endif
                  ity = 0
                  do i=1,99
                     if (tocapf(str).eq.tocapf(elemnt(i))) ity = i
                  end do
                  if (ity.ne.0) then
                     felem = .true.
                     natt = natt + 1
                     ieltyp(natt) = ity
                  endif 
               else
                  goto 10
               endif
            end do
         endif
      endif

10    continue

c     scale factor

      if (getlin(0).eq.1) then
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.3) then
            scl = rtype
         endif
      endif

c     unit vectors

      if (getlin(0).eq.1) then
         if (.not.gnreal(v1,3,.false.)) goto 100
      else
         goto 100
      endif

      if (getlin(0).eq.1) then
         if (.not.gnreal(v2,3,.false.)) goto 100
      else
         goto 100
      endif

      if (getlin(0).eq.1) then
         if (.not.gnreal(tmp,3,.false.)) goto 100
      else
         goto 100
      endif

      r(1) = vlen(v1)*scl
      r(2) = vlen(v2)*scl
      r(3) = vlen(tmp)*scl

      if (r(1).le.0.0d0.or.r(2).le.0.0d0.or.r(3).le.0.0d0) goto 100

      a = r(1)
      b = r(2)
      c = r(3)
      call impsc(v2,tmp,csa)
      call impsc(v1,tmp,csb)
      call impsc(v1,v2,csc)
      alpha = dacos(csa)*todeg
      beta = dacos(csb)*todeg
      gamma = dacos(csc)*todeg
      nspg = 1

      call prcell(nspg,a,b,c,alpha,beta,gamma)
      call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)
      call cprot(nspg,nopr,icent,ir,it,.false.)
      if (idebug.eq.1) call prop(nopr,ir,it)

c In new POSCARs there might be a line here to give the elements names

      if (getlin(0).eq.1) then
         
         call gtplin(strtmp)
         ktype = nxtwrd(str,nstr,itype,rtype)
         
         if (ktype.eq.1) then
            
            if (.not.felem) then
               call setlin(strtmp,0)
               natt = 0
               do while(.true.)
                  ktype = nxtwrd(str,nstr,itype,rtype)
                  if (ktype.eq.1.and.(nstr.eq.1.or.nstr.eq.2)) then
                     if (nstr.eq.1) then
                        str(2:2) = str(1:1)
                        str(1:1) = ' '
                     endif
                     ity = 0
                     do i=1,99
                        if (tocapf(str).eq.tocapf(elemnt(i))) ity = i
                     end do
                     if (ity.ne.0) then
                        felem = .true.
                        natt = natt + 1
                        ieltyp(natt) = ity
                     endif 
                  else
                     goto 11
                  endif
               end do
 11            continue
            endif

c     we now read the number of atoms of each sort

            natoms = 0
            ity = 0
            if (getlin(0).eq.1) then
               do i=1,50
                  ktype = nxtwrd(str,nstr,itype,rtype)
                  if (ktype.eq.2) then
                     natoms = natoms + itype
                     ity = ity + 1
                     ielnum(ity) = itype
                  else
                     goto 20
                  endif
               end do
            else
               goto 100
            endif
            
 20         continue

         elseif (ktype.eq.2) THEN

c in fact, it is already the number of atoms of each sort   
c we read them from the stored line

            call setlin(strtmp,0)
            natoms = 0
            ity = 0
            do i=1,50
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.2) then
                  natoms = natoms + itype
                  ity = ity + 1
                  ielnum(ity) = itype
               else
                  goto 21
               endif
            end do
          
 21         continue
         endif
      endif

      natx = natoms

c      if (natt.eq.ity) then
      if (natt.ge.ity) then
         natt=ity
         k = 0
         do i=1,natt
            do j=1,ielnum(i)
               k = k + 1
               ianz(k) = ieltyp(i)
               nat(k)  = ieltyp(i)
            end do
         end do
      endif

c Check for Direct (=fractional coordinates)

      if (getlin(0).eq.1) then
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.1) then
            call tocap(str,nstr)
            if (str(1:1).eq.'S') then
               if (getlin(0).eq.1) then
                  ktype = nxtwrd(str,nstr,itype,rtype)
                  if (ktype.ne.1) goto 100
               endif
            endif

            if (str(1:1).ne.'C'.and.str(1:1).ne.'c'.and.
     &          str(1:1).ne.'K'.and.str(1:1).ne.'k') then
               str  = 'DIRECT'
               nstr = 6
            endif
         else
            goto 100
         endif
      end if

c read in fractional coordinates

      do i=1,natoms
         if (natt.ne.ity) then
            nat(i) = 99
            ianz(i) = 99
         endif
         if (getlin(0).eq.1) then
            if (.not.gnreal(xyzt,3,.false.)) goto 100
            do j=1,3
               xyz(j,i) = v1(j)*scl*xyzt(1) + v2(j)*scl*xyzt(2) +
     &                    tmp(j)*scl*xyzt(3)
               coo(j,i) = xyzt(j)
            end do
            !print*, nat(i),(xyz(j,i),j=1,3)
         else
            goto 100            
         endif
      end do

c      if (idocub.eq.0) then

         ncon = 1
         inorm = 0 
         iatoms = natoms
         nnat = natoms
         norg = natoms
         nstor = mxnat-iatoms
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

         call doconn
         call dohcon(0)
   
         do i=1,iatoms
            do j=1,iconn(1,i)+1
               iconn(j,nstor+i) = iconn(j,i)
            end do
         end do
c      endif

      ichx = 1

      if (idocub.eq.0) then
         istat = 2
         return
      endif
      
      read(iun2,'(a)') line

c x,y,z dimension of grid

      if (getlin(0).eq.1) then
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            npts1 = itype
         else
            goto 100
         endif
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            npts2 = itype
         else
            goto 100
         endif
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            npts3 = itype
         else
            goto 100
         endif
      else
         goto 100
      endif

      if (npts1.gt.mx3d.or.npts2.gt.mx3d.or.npts3.gt.mx3d) then
          lstr = 'dimension greater than maximum !'
          print*,'npts1 ',npts1,' npts2 ',npts2,' npts3 ',npts3,
     &           ' maxdim ',mx3d 
          goto 100
      endif

      px = (scl*v1(1) + scl*v2(1) + scl*tmp(1))/2.0d0
      py = (scl*v1(2) + scl*v2(2) + scl*tmp(2))/2.0d0
      pz = (scl*v1(3) + scl*v2(3) + scl*tmp(3))/2.0d0

c      do i=1,natoms
c         xyz(1,i) = xyz(1,i) + px
c         xyz(2,i) = xyz(2,i) + py
c         xyz(3,i) = xyz(3,i) + pz
c      end do

      call vsc1(v1,1.0d0,1.0d-4)
      call vsc1(v2,1.0d0,1.0d-4)
      call vsc1(tmp,1.0d0,1.0d-4)


      cx = tmp(1)
      cy = tmp(2)
      cz = tmp(3)

      if (idebug.eq.1) then
         print*,'r=',(r(i),i=1,3)
         print*,'px py pz ',px,py,pz
         print*,'cx cy cz ',cx,cy,cz
         print*,'v1=',(v1(i),i=1,3)
         print*,'v2=',(v2(i),i=1,3)
         do i=1,natoms
            print*, nat(i),(xyz(j,i),j=1,3)
         end do
      endif

      call crprod(v1,v2,dtmp)
      call impsc(dtmp,tmp,cosb)
      if (1.0d0-dabs(cosb).gt.0.001d0) then
          lstr = 'Only rectangular grids supported !'
          print*,'ERROR: Only rectangular grids supported ',cosb
          goto 100
      endif

      j = 1
      k = 1
      ijkt = 0
      iposng = 0
      nall = npts1*npts2*npts3
      ib = 0
      nb = (npts1/lnbuck)*lnbuck
      if (npts1.gt.nb) nb = nb + lnbuck

      do while (ijkt.lt.nall)
c fill bucket
c          ijkt = ijk
          do while (ib.lt.nb.and.ijkt.lt.nall)
             n = lnbuck
             if (nall-ijkt.lt.n) n = nall-ijkt
c             print*,'n=',n,' bucket',bucket(ib),' ',nall,' ',ijkt
             read(21,'(10f8.3)') (bucket(ib+l),l=1,n)
             ib = ib + n
             ijkt = ijkt + n
          end do
c full bucket
          ib = ib - npts1
c fill density array
          do i=1,npts1
             denn((k-1)*mx3d2 + (i-1)*npts2+j) = bucket(i)
          end do
c          print*,j,' ',k
c          print*,(bucket(l),l=1,npts1)
c move access bucket to beginning of bucket
          do l=1,ib
             bucket(l) = bucket(npts1+l)
          end do
          j = j + 1
          if (j.gt.npts2) then
             j = 1
             k = k + 1
          endif
c          ijk = ijk + npts1
      end do

      iplat = 0

      ca = v2(2)
      sa = -v2(1)
      sb = -v1(3)
      if (dabs(ca).gt.0.001d0) then
         cb = v1(1) / ca
      else
         cb = v1(2) / sa
      endif

      cc = 1.0d0
      sc = 0.0d0


      call parrat
      call proato

      do i=1,npts3
          pmnn(i) = 100000.d0
      end do

      do i=1,natoms
         nconn(i) = 0
         do j=1,natoms
           if (i.ne.j) then
            dmaxsq = (vdwr(nat(i)) + vdwr(nat(j)))**2
            dijsq = ((xyz(1,i)-xyz(1,j))*toang)**2
     &             +((xyz(2,i)-xyz(2,j))*toang)**2
     &             +((xyz(3,i)-xyz(3,j))*toang)**2
            if (dijsq.lt.dmaxsq) then
                nconn(i) = nconn(i) + 1
            endif
           endif
         end do
      end do

      if (idocub.eq.1) then
          close(21)
          call xyzcoo(1,0,0)
      endif

      ofrst = .false.

      if (idocub.eq.1) iun2 = iuntmp

      lstr = 'found cube file'
      call inferr(lstr,0)
      return

100   if (idocub.eq.1) close(21)
110   if (idocub.eq.1) iun2 = iuntmp
      istat = 0
      call inferr(lstr,0)

      return
      end
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$      
c     ====================================================================      
      subroutine rdcqd(npts1,npts2,npts3,iposng,istat,lenf,
     &                  idocub,idebug,denn,pmnn,bucket,
     &                  coo,ianz,iatclr,iconn,
     &                  nnat,norg,icent,inorm,ncon,nspg,ichx,
     &                  nopr,ir,it,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision ( a-h,o-z)

c THIS IS REALLY rdconquest

      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (mxcon=10)
      parameter (lnbuck=10)
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /grdhlp/ mx3d,mx3d2
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /eulx/   ca,cb,sa,sb,cc,sc
      integer*2 ir,it
      logical oclose,osshad,ofill,ofrst
      common /spacom/ grey,iscol,oclose,osshad,ofill,ofrst
      common /zproj/  zval(numatm),nindx(numatm),nconn(numatm)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /athlp/ iatoms, mxnat
      character*2 elemnt
      common /elem/elemnt(mxel)

      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line,lstr
      logical opfil,gnreal
      integer getlin
      character*137 str
      character*2 tocapf
      character keywrd*320, keyori*320
      common /keywrd/ keywrd,keyori
      dimension tmp(3),dtmp(10),xyzt(3)
      dimension ieltyp(20),ielnum(20)
      dimension denn(*),pmnn(*),bucket(*)
      dimension coo(3,*),ianz(*),iatclr(*),iconn(mxcon+1,*),
     &          ir(3,3,192),it(3,192)

      logical felem
      character*137 strtmp
      !common /vasp/nheadx,natx

      istat = 1
      itype = 0
      toang  = 0.52917706d0
      todeg = 45.0d0 / datan(1.0d0)
      lstr = ' '

      print*, '*** begin of Conquest READ ***'   
      print*, '    WARNINGS:'
      print*, '     > coordinates must be fractional'
      print*, '     > cell vectors in Bohr unit' 

      if (idocub.eq.1) then
         iuntmp = iun2
         iun2 = 21
         
         if (.not.opfil(21,keywrd(1:lenf),lenf,1,1,0)) then
            lstr = keywrd(1:lenf)//' non-existent'
            goto 110
         endif
      endif


c     unit vectors
      
      if (getlin(0).eq.1) then
         if (.not.gnreal(v1,3,.false.)) goto 100
         r(1) = vlen(v1)
      else
         goto 100
      endif
      
      if (getlin(0).eq.1) then
         if (.not.gnreal(v2,3,.false.)) goto 100
         r(2) = vlen(v2)
      else
         goto 100
      endif
      
      if (getlin(0).eq.1) then
         if (.not.gnreal(tmp,3,.false.)) goto 100
         r(3) = vlen(tmp)
      else
         goto 100
      endif

      r(1) = vlen(v1)
      r(2) = vlen(v2)
      r(3) = vlen(tmp)
      
      if (r(1).le.0.0d0.or.r(2).le.0.0d0.or.r(3).le.0.0d0) goto 100

      a = r(1)
      b = r(2)
      c = r(3)
      call impsc(v2,tmp,csa)
      call impsc(v1,tmp,csb)
      call impsc(v1,v2,csc)
      alpha = dacos(csa)*todeg
      beta  = dacos(csb)*todeg
      gamma = dacos(csc)*todeg
      nspg  = 1

      call prcell(nspg,a,b,c,alpha,beta,gamma)
      call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)
      call cprot(nspg,nopr,icent,ir,it,.false.)
      if (idebug.eq.1) call prop(nopr,ir,it)
      
c     get atom number and types

      natt  = 0      
      felem = .false.      
      !
      if (getlin(0).eq.1) then
         !
         ktype = nxtwrd(str,nstr,itype,rtype)
         !
         if (ktype.eq.2) then
            natoms = itype
            !
         end if
         !
         do while(.true.)
            !
            ktype = nxtwrd(str,nstr,itype,rtype)
            !
            if (ktype.eq.1.and.(nstr.eq.1.or.nstr.eq.2)) then
               if (nstr.eq.1) then
                  str(2:2) = str(1:1)
                  str(1:1) = ' '
                  !
               endif
               !
               ity = 0
               do i=1,99
                  if (tocapf(str).eq.tocapf(elemnt(i))) ity = i
               end do
               !
               if (ity.ne.0) then
                  felem = .true.
                  natt = natt + 1
                  ieltyp(natt) = ity
                  !
               endif
               !
            else
               goto 10
               !
            endif
            !
         end do
         !
      endif
      
 10   continue

c     get atom positions: fractional coordinate assumed!

      do i = 1, natoms
         
         if (natt.ne.ity) then
            nat(i)  = 99
            ianz(i) = 99
         endif
         
         if (getlin(0).eq.1) then
            
            if (.not.gnreal(xyzt,3,.false.)) goto 100
            
            do j = 1, 3
               xyz(j,i) = v1(j)*xyzt(1) + v2(j)*xyzt(2) +
     &                    tmp(j)*xyzt(3)
               coo(j,i) = xyzt(j)
            end do

            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype .ne. 2) goto 100

            do j = 1, natt
               if (itype .eq. j) then
                  !print*, itype, ie  ltyp(j)
                  ianz(i) = ieltyp(j)
                  nat(i)  = ieltyp(j)
               end if
            end do
            
            print*, nat(i), (xyzt(j), j=1,3 ), itype, ktype
            
         else
            goto 100
         endif
         
      end do

      ncon   = 1
      inorm  = 0 
      iatoms = natoms
      nnat   = natoms
      norg   = natoms
      nstor  = mxnat-iatoms
      
      do i = 1, iatoms
         
         do j=1,3
            coo(j,nstor+i) = coo(j,i)
         end do
         
         call fr2crt(coo(1,i),xa,ya,yb,za,zb,zc)
         do j=1,3
            coo(j,i) = coo(j,i)/toang
         end do
         ianz(nstor+i)   = ianz(i)
         iatclr(nstor+i) = 1
         
      end do
      
      call doconn
      call dohcon(0)
      
      do i=1,iatoms
         do j=1,iconn(1,i)+1
            iconn(j,nstor+i) = iconn(j,i)
         end do
      end do
         
      ichx = 1

      if (idocub.eq.0) then
         print*, '*** end of Conquest READ ***' 
         istat = 2
         return
      endif

      
100   if (idocub.eq.1) close(21)
110   if (idocub.eq.1) iun2 = iuntmp

      istat = 0
      call inferr(lstr,0)

      return
      end

      
c     ====================================================================      
      subroutine wrcqd(iun,coo,ianz,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision ( a-h,o-z)
      parameter (mxel=100)
      parameter (mxetyp=20)
      common /athlp/ iatoms, mxnat
      character*2 elemnt
      common /elem/elemnt(mxel)
      logical fndel
      dimension tmp(3),tr(3),rr(3,3),ieltyp(mxetyp),ielnum(mxetyp)
      dimension coo(3,*),ianz(*)

      toang  = 0.52917721092d0

      print*, '*** begin of Conquest WRITE ***'
      print*, '    WARNINGS:'
      print*, '     > coordinates will be fractional'
      print*, '     > cell vectors in Bohr unit'

      natt = 0
      do i=1,mxetyp
         ielnum(i) = 0
      end do

      natoms = 0
      do i=1,iatoms

         if (ianz(i).gt.0.and.ianz(i).lt.99) then

            fndel = .false.
            do j=1,natt
               if (ianz(i).eq.ieltyp(j)) then
                  ielnum(j) = ielnum(j) + 1
                  fndel = .true.
               endif
            end do

            if (.not.fndel) then
               natt = natt + 1
               ieltyp(natt) = ianz(i)
               ielnum(natt) = 1
            endif

            natoms = natoms + 1
         endif

      end do

c      write(iun,'(20(a2,1x))') (elemnt(ieltyp(i)),i=1,natt)

      scl = 1.0d0
c      write(iun,'(f12.6)') scl

      call setrr(alpha,beta,gamma,a,b,c,rr)

      tmp(1) = 1.0d0
      tmp(2) = 0.0d0
      tmp(3) = 0.0d0
c      call fr2crt(tmp,xa,ya,yb,za,zb,zc)
      do k=1,3
         tr(k) = trc(tmp,rr,k)
      end do

c      write(iun,'(3(f12.6,1x))') (tmp(i),i=1,3)
      write(iun,'(3(f20.12,1x))') (tr(i),i=1,3)
      
      tmp(1) = 0.0d0
      tmp(2) = 1.0d0
      tmp(3) = 0.0d0
c      call fr2crt(tmp,xa,ya,yb,za,zb,zc)
      do k=1,3
         tr(k) = trc(tmp,rr,k)
      end do
c      write(iun,'(3(f12.6,1x))') (tmp(i),i=1,3)
      write(iun,'(3(f20.12,1x))') (tr(i),i=1,3)

      tmp(1) = 0.0d0
      tmp(2) = 0.0d0
      tmp(3) = 1.0d0
c      call fr2crt(tmp,xa,ya,yb,za,zb,zc)
      do k=1,3
         tr(k) = trc(tmp,rr,k)
      end do
c      write(iun,'(3(f12.6,1x))') (tmp(i),i=1,3)
      write(iun,'(3(f20.12,1x))') (tr(i),i=1,3)
 
      write(iun,'(i12,4x,20(a2,1x))') natoms,
     &     (elemnt(ieltyp(i)),i=1,natt)


c      write(iun,'(20(i3,1x))') (ielnum(i),i=1,natt)
c      write(iun,'(20(a2,1x))') (elemnt(ieltyp(i)),i=1,natt)


c      write(iun,'(a)') 'Direct'

      do k=1,natt
         !i_elnum = 0
         !n_elnum = ielnum(k)
         do i=1,iatoms

            if (ianz(i).eq.ieltyp(k)) then
               !i_elnum = i_elnum + 1
               !
               do j=1,3
                  tmp(j) = coo(j,i) * toang
               end do

               call crt2fr(tmp,tr,xa,ya,yb,za,zb,zc)

               write(iun,'(3(f16.10,1x),4x,i4,2x,a8)') (tr(j),j=1,3),
     &              k,' T  T  T'

            endif

         end do

      end do

      call inferr('Wrote file: Conquest_coord',0)

      print*, '*** end of Conquest WRITE ***'

      return
      end
      
c     ====================================================================
      subroutine wrvasd(iun,coo,ianz,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision ( a-h,o-z)
      parameter (mxel=100)
      parameter (mxetyp=20)
      common /athlp/ iatoms, mxnat
      character*2 elemnt
      common /elem/elemnt(mxel)
      logical fndel
      dimension tmp(3),tr(3),rr(3,3),ieltyp(mxetyp),ielnum(mxetyp)
      dimension coo(3,*),ianz(*)

      toang  = 0.52917706d0

      natt = 0
      do i=1,mxetyp
         ielnum(i) = 0
      end do

      do i=1,iatoms

         if (ianz(i).gt.0.and.ianz(i).lt.99) then

            fndel = .false.
            do j=1,natt
               if (ianz(i).eq.ieltyp(j)) then
                  ielnum(j) = ielnum(j) + 1
                  fndel = .true.
               endif
            end do

            if (.not.fndel) then
               natt = natt + 1
               ieltyp(natt) = ianz(i)
               ielnum(natt) = 1
            endif

         endif

      end do

      write(iun,'(20(a2,1x))') (elemnt(ieltyp(i)),i=1,natt)

      scl = 1.0d0
      write(iun,'(f12.6)') scl

      call setrr(alpha,beta,gamma,a,b,c,rr)

      tmp(1) = 1.0d0
      tmp(2) = 0.0d0
      tmp(3) = 0.0d0
c      call fr2crt(tmp,xa,ya,yb,za,zb,zc)
      do k=1,3
         tr(k) = trc(tmp,rr,k)
      end do

c      write(iun,'(3(f12.6,1x))') (tmp(i),i=1,3)
      write(iun,'(3(f12.6,1x))') (tr(i),i=1,3)
      
      tmp(1) = 0.0d0
      tmp(2) = 1.0d0
      tmp(3) = 0.0d0
c      call fr2crt(tmp,xa,ya,yb,za,zb,zc)
      do k=1,3
         tr(k) = trc(tmp,rr,k)
      end do
c      write(iun,'(3(f12.6,1x))') (tmp(i),i=1,3)
      write(iun,'(3(f12.6,1x))') (tr(i),i=1,3)

      tmp(1) = 0.0d0
      tmp(2) = 0.0d0
      tmp(3) = 1.0d0
c      call fr2crt(tmp,xa,ya,yb,za,zb,zc)
      do k=1,3
         tr(k) = trc(tmp,rr,k)
      end do
c      write(iun,'(3(f12.6,1x))') (tmp(i),i=1,3)
      write(iun,'(3(f12.6,1x))') (tr(i),i=1,3)

      write(iun,'(20(i3,1x))') (ielnum(i),i=1,natt)

      write(iun,'(a)') 'Direct'

      do k=1,natt

         do i=1,iatoms
            if (ianz(i).eq.ieltyp(k)) then
               do j=1,3
                  tmp(j) = coo(j,i) * toang
               end do
               call crt2fr(tmp,tr,xa,ya,yb,za,zb,zc)
               write(iun,'(3(f12.6,1x))') (tr(j),j=1,3)
            endif
         end do

      end do

      call inferr('Wrote file: POSCAR',0)

      return
      end

c     ====================================================================      
      subroutine wrmopd(iun,coo,ianz,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision ( a-h,o-z)
      parameter (mxel=100)
      parameter (mxetyp=20)
      common /athlp/ iatoms, mxnat
      character*2 elemnt
      common /elem/elemnt(mxel)
      dimension tmp(3)
c      dimension tr(3),rr(3,3)
      dimension coo(3,*),ianz(*)

      toang  = 0.52917706d0

      write(iun,'(a)') 'AM1 NOINTER XYZ'
      write(iun,'(a)') ' '
      write(iun,'(a)') ' '

c      call setrr(alpha,beta,gamma,a,b,c,rr)



      do i=1,iatoms
          do j=1,3
             tmp(j) = coo(j,i) * toang
          end do
          if (ianz(i).gt.0.and.ianz(i).lt.99) then
             if (i.eq.1) then
                write(iun,'(a2,1x,3(f12.6,a3))') elemnt(ianz(i)),
     &            tmp(1),' 0 ',tmp(2),' 0 ',tmp(3),' 0 '
             else
                write(iun,'(a2,1x,3(f12.6,a3))') elemnt(ianz(i)),
     &            tmp(1),' 1 ',tmp(2),' 1 ',tmp(3),' 1 '
             endif
          endif
      end do



      tmp(1) = 1.0d0
      tmp(2) = 0.0d0
      tmp(3) = 0.0d0
      call fr2crt(tmp,xa,ya,yb,za,zb,zc)
c      do k=1,3
c         tr(k) = trc(tmp,rr,k)
c      end do

      write(iun,'(a3,3(f12.6,a3))') 
     &    'Tv ',tmp(1),' 1 ',tmp(2),' 0 ',tmp(3),' 0 '
c      write(iun,'(3(f12.6,1x))') (tr(i),i=1,3)
      
      tmp(1) = 0.0d0
      tmp(2) = 1.0d0
      tmp(3) = 0.0d0
      call fr2crt(tmp,xa,ya,yb,za,zb,zc)
c      do k=1,3
c         tr(k) = trc(tmp,rr,k)
c      end do
c      write(iun,'(3(f12.6,1x))') (tmp(i),i=1,3)
      write(iun,'(a3,3(f12.6,a3))') 
     &    'Tv ',tmp(1),' 0 ',tmp(2),' 1 ',tmp(3),' 0 '
c      write(iun,'(3(f12.6,1x))') (tr(i),i=1,3)

      tmp(1) = 0.0d0
      tmp(2) = 0.0d0
      tmp(3) = 1.0d0

      call fr2crt(tmp,xa,ya,yb,za,zb,zc)

c      do k=1,3
c         tr(k) = trc(tmp,rr,k)
c      end do
c      write(iun,'(3(f12.6,1x))') (tmp(i),i=1,3)

      write(iun,'(a3,3(f12.6,a3))') 
     &    'Tv ',tmp(1),' 0 ',tmp(2),' 0 ',tmp(3),' 1 '
c      write(iun,'(3(f12.6,1x))') (tr(i),i=1,3)

      call inferr('Wrote file: mopac.dat',0)

      return
      end

c     ====================================================================
      subroutine rdgrdd(npts1,npts2,npts3,iun,istat,
     &                  denn,dens,pmnn,buff)
      implicit double precision (a-h,o-z)
      parameter (mxdir=500)
      logical oclose,osshad,ofill,ofrst
      common /spacom/ grey,iscol,oclose,osshad,ofill,ofrst
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /grdhlp/ mx3d,mx3d2
      common /comsrf/  vo(3), rr(3),vv1(3),vv2(3),vv3(3), wo(3),
     &                 sl(3),isl
      common /cnthlp/ rng1,rng2,vlcnt,vlcnt2
      real ra,orig,den
      character*137 lstr
      character*4 buff(*)
      character*72 header
      integer currec,bigend,buflen
      common /rdrec/ iund,currec
      common /isend/ bigend,nlen
      dimension orig(3)
      dimension denn(*),pmnn(*),dens(*)

c nnx,nny,nnz    the number of points in x, y and z direction
c ra             real step along x,y and z direction per point
c rx,ry,rz       origin of the grid

      toang  = 0.52917706d0
      istat = 1
      lstr = ' '

      iplat = 0

      buflen = mx3d2
      rewind(iund)
      call getrec(buff,buflen,0,ierr)

      do i =1,18
         header(4*(i-1)+1:4*i) = buff(i)(1:4)
      end do

      print*,'header = ',header

      call byter(buff(26),npts2)
      call byter(buff(27),npts1)
      call byter(buff(28),npts3)

      print*,'npts1=',npts1,' npts2=',npts2,' npts3=',npts3

      call byter(buff(29),ra)

      ra = ra / toang

      print*,'ra=',ra

      call byter(buff(30),orig(2))
      call byter(buff(31),orig(1))
      call byter(buff(32),orig(3))

      print*,'orgx,orgy,orgz ',orig(1),' ',orig(2),' ',orig(3)

      v1(1) = 0.0d0
      v1(2) = 1.0d0
      v1(3) = 0.0d0

      v2(1) = 1.0d0
      v2(2) = 0.0d0
      v2(3) = 0.0d0

      r(1) = dble(ra)*dble(npts1-1)
      r(2) = dble(ra)*dble(npts2-1)
      r(3) = dble(ra)*dble(npts3-1)

      do i=1,3
         orig(i) = orig(i) / toang
         rr(i)  = r(i)
         vv1(i) = v1(i)
         vv2(i) = v2(i)
         vv3(i) = 0.0d0
         sl(i) = 1.0d0
      end do

      vv3(3) = 1.0d0

      isl = 1
      sl(1) = -1.0d0
      sl(2) = -r(2)/r(1)
      sl(3) = -r(3)/r(1)

      call vnrm(sl)

      px =  dble(orig(2))
      py =  dble(orig(1))
      pz =  dble(orig(3))


      vo(1) = px
      vo(2) = py
      vo(3) = pz

      cx = 0.0d0
      cy = 0.0d0
      cz = 1.0d0

      pmax = -1000000.d0
      pmin = 1000000.d0

      do k=1,npts3

c get plate

         call getrec(buff,buflen,0,ierr)
         if (ierr.eq.1) goto 100

         call byter(buff(1),nz)
         call byter(buff(2),nx)
         call byter(buff(3),ny)

c         print*,'nz,nx,ny ',nz,' ',nx,' ',ny

         call getrec(buff,buflen,0,ierr)

         ij = 0
         do i=1,ny
            do j=1,nx
               ij = ij + 1
               call byter(buff(ij),den)
               dens((j-1)*ny+i) = dble(den)
            end do
         end do
         do i=1,nx*ny
            denn((nz-1)*mx3d2 + i) = dens(i)
            if (dens(i).gt.pmax) pmax = dens(i)
            if (dens(i).lt.pmin) pmin = dens(i)
         end do

      end do

      rng1 = pmin
      rng2 = pmax

      call messg(22)

      ofrst = .false.

      lstr = 'found grid file'
      call inferr(lstr,0)
      return

100   istat = 0
      lstr = 'error reading GRID file'
      call inferr(lstr,0)
      return
      end

c     ====================================================================      
      subroutine rdomad(npts1,npts2,npts3,iun,istat,
     &                  denn,dens,pmnn,ichx)
      implicit double precision (a-h,o-z)
      parameter (mxdir=500)
      logical dolap, lapdbl
      common /vropt/ ivtwo,ihand,ivadd
      common /comsrf/  vo(3), r(3),v1(3),v2(3),v3(3), wo(3),
     &                 sl(3),isl
      common /grdhlp/ mx3d,mx3d2
      common /cnthlp/ rng1,rng2,vlcnt,vlcnt2
      character*137 lstr
      integer*2 buff(256),i1to2
      integer*1 bff(512)
      integer currec,bigend
      common /rdrec/ iund,currec
      common /isend/ bigend,nlen
      equivalence (buff,bff)
      dimension denn(*),pmnn(*),dens(*), fdum(3)

c nnx,nny,nnz    the number of points in x, y and z direction
c ra             real step along x,y and z direction per point
c rx,ry,rz       origin of the grid

      istat = 1
      lstr = ' '

      toang  = 0.52917706d0
      iplat = 0

      currec = -1
      call getrc2(buff,ierr)

c 7-9 number of grid points in x,y,z

      npts1 = int(buff(7))
      npts2 = int(buff(8))
      npts3 = int(buff(9))

      print*,' '
      print*,'O/DSN6 map file:'
      print*,' '
      print*,'npts1 npts2 npts3 ',npts1,npts2,npts3

c 1-3 lower limits for x,y,z

      ior1 = int(buff(1)) 
      ior2 = int(buff(2)) 
      ior3 = int(buff(3)) 

      print*,'orgx,orgy,orgz ',ior1,ior2,ior3

c 4-6 limit range for x,y,z

      iex1 = int(buff(4)) 
      iex2 = int(buff(5)) 
      iex3 = int(buff(6)) 

      print*,'iex1,iex2,iex3 ',iex1,iex2,iex3

c     the grid does not cover the whole unit cell, but just the molecule
c     the grid therefor does not run from 1 to npts in each direction but
c     from ior1 to ior1 + iex1 -1 (in the x direction)
c     we need npts1,npts2,npts3 however to calculate the spacing between
c     two consequetive grid points in all three directions 
c     (together with a,b,c)
c
c     the problem with existing molden maps is they have equally many
c     points to the right and left of the origin (vo)
c     in these new maps the origin is not the center of the map, but
c     the start of the map in all three directions (npts1=npts2=npts3=1)


c 18 F1

      sc  =  1.0d0 / dble(buff(18))
      sca  =  1.0d0 / (dble(buff(18)) * toang)

c     the spacing between two consequetive grid points in all three directions 

c 10-12 cell dimensions (times F1) in x,y,z
c sca gets F1 (buff(18) out again 

      a = sca * dble(buff(10)) / dble(npts1)
      b = sca * dble(buff(11)) / dble(npts2)
      c = sca * dble(buff(12)) / dble(npts3)

      print*,'a b c',sc*dble(buff(10)),sc*dble(buff(11)),
     *               sc*dble(buff(12))

c 13-15 cell angles (times F1) in x,y,z
c sc gets F1 (buff(18) out again 

      alpha = sc * dble(buff(13)) * 3.14159265d0 / 180.0d0
      beta = sc * dble(buff(14)) * 3.14159265d0 / 180.0d0
      gamma = sc * dble(buff(15)) * 3.14159265d0 / 180.0d0

c 16 multiplicative term (times F2) for density value to map to [0,255]
c 17 additive term (times F2) for density value to map to [0,255]
c 19 F2

      prod = dble(buff(16)) / dble(buff(19))
      plus = dble(buff(17)) / dble(buff(19))

      print*,'alpha beta gamma ',sc*dble(buff(13)),sc*dble(buff(14)),
     &                           sc*dble(buff(15))

      v1(1) = a
      v1(2) = 0.0d0
      v1(3) = 0.0d0
    
      v2(1) = dcos(gamma) * b
      v2(2) = dsin(gamma) * b
      v2(3) = 0.0d0
    
      z1 = dcos(beta)
      z2 = (dcos(alpha) - dcos(beta)*dcos(gamma)) / dsin(gamma)
      z3 = dsqrt(1.0 - z1*z1 - z2*z2);

      v3(1) = z1 * c
      v3(2) = z2 * c
      v3(3) = z3 * c
    
c     Convert the origin from grid space to cartesian coordinates

      vo(1) = v1(1) * dble(ior1) + 
     &        v2(1) * dble(ior2) + v3(1) * dble(ior3)
      vo(2) = v2(2) * dble(ior2) + v3(2) * dble(ior3)
      vo(3) =                      v3(3) * dble(ior3)
     
c     normalise basis vectors

      call vsc1(v1,1.0d0,1.0d-4)
      call vsc1(v2,1.0d0,1.0d-4)
      call vsc1(v3,1.0d0,1.0d-4)

      isl = 1
      sl(1) = 1.0d0
      sl(2) = b/a
      sl(3) = c/a
      call vnrm(sl)

      r(1) = a*dble(iex1-1)
      r(2) = b*dble(iex2-1)
      r(3) = c*dble(iex3-1)

      print*,'r ',(r(i),i=1,3)
      npts1 = iex1
      npts2 = iex2
      npts3 = iex3

c     spacing between consequetive points in the three directions
c     rf1 = r(1) / dble(iex1-1)
c     rf2 = r(2) / dble(iex2-1)
c     rf3 = r(3) / dble(iex3-1)
c
c     do kc=1,npts3 do ic=1,npts2 do jc=1,npts3
c     coordinates of each grid point (jc,ic,kc) rt(1,2,3)
c     rt(1) = v1(1)*(jc-1)*rf1 + v2(1)*(ic-1)*rf2 + v3(1)*(kc-1)*rf3 + vo(1)
c     rt(2) = v1(2)*(jc-1)*rf1 + v2(2)*(ic-1)*rf2 + v3(2)*(kc-1)*rf3 + vo(2)
c     rt(3) = v1(3)*(jc-1)*rf1 + v2(3)*(ic-1)*rf2 + v3(3)*(kc-1)*rf3 + vo(3)
c
c     we have to change rotbck to:
c      do i=1,3
c         wn(i) = v1(i)*(w1)*r(1) +
c     &           v2(i)*(w2)*r(2) +
c     &           v3(i)*(w3)*r(3) + vo(i)
c      end do
c
c alternatively we could make it a general purpose routine by specifying
c in common comsrf what is the location of the origin in grid coordinates
c

      pmax = -1000000.d0
      pmin = 1000000.d0

      ibrik = npts1/8
      if (npts1.gt.ibrik*8) ibrik = ibrik + 1
      jbrik = npts2/8
      if (npts2.gt.jbrik*8) jbrik = jbrik + 1
      kbrik = npts3/8
      if (npts3.gt.kbrik*8) kbrik = kbrik + 1

      do k=1,kbrik
         do j=1,jbrik
            do i=1,ibrik

               call getrc2(buff,ierr)

               do iz=1,8
                  if ((iz + (k-1)*8).le.npts3) then
                     do iy=1,8
                        if ((iy + (j-1)*8).le.npts2) then
                           do ix=1,8

                              if ((ix + (i-1)*8).le.npts1) then

                                ipnt = ix        +
     &                                 (iy-1)*8  +
     &                                 (iz-1)*64 
                                den = dble(i1to2(bff(ipnt))) - plus
                                den = den / prod

                                denn( (ix+(i-1)*8)          +
     &                               ((iy+(j-1)*8)-1)*npts1 +
     &                               ((iz+(k-1)*8)-1)*mx3d2) 
     &                          = den

                                if (den.gt.pmax) pmax = den
                                if (den.lt.pmin) pmin = den

                              endif
                        
                           end do
                        endif
                     end do
                  endif
               end do

            end do
         end do
      end do

      do i=1,npts3
         pmnn(i) = 100000.d0
      end do

      print*,'Dens Max=',pmax,' Dens Min=',pmin
      rng1 = pmin
      rng2 = pmax

      vlcnt = 1.00d0
      lstr = 'found omap file'
      call inferr(lstr,0)

      ichx = 1

      close(iun)
      call messg(16)

      return

100   istat = 0
      lstr = 'error reading OMAP/DSN6 file'
      call inferr(lstr,0)
      return
      end

c     ====================================================================      
      subroutine dpomad(iopt,denn)
      implicit double precision (a-h,o-z)
      logical dolap, lapdbl
      common /comsrf/  vo(3), r(3),v1(3),v2(3),v3(3), wo(3),
     &                 sl(3),isl
      common /cnthlp/ rng1,rng2,vlcnt,vlcnt2
      integer srfmap, srfloc
      common /hlpsrf/ npts1,npts2,npts3,srfmap,srfloc,ifogl,itsrf
      common /vropt/ ivtwo,ihand,ivadd
      dimension origin(3),rt(3)
      dimension denn(*),fdum(3)

      call curs(1)
      wo(1) = 0.0d0
      wo(2) = 0.0d0
      wo(3) = 0.0d0

      do i=1,3
         origin(i) = 0.0d0
         rt(i) = 1.0d0
      end do

      ipsi   = 0
      dolap  = .false.
      lapdbl = .false.
      if (iopt.eq.1) then
         itrans = 0
      else
         itrans = 1
      endif
      mapit  = 0

      ivtwot = ivtwo
      ivtwo = 4

      call mcubes(npts1,npts2,npts3,denn,fdum,mapit,origin,vlcnt,
     &                      vlcnt2,rt,ipsi,dolap,lapdbl,itrans,iun)

      ivtwo = ivtwot

      call curs(0)

      return
      end

      subroutine byter(intin,intout)
      implicit none
      integer intin, intout,bigend,nlen
      common /isend/ bigend,nlen

      if (bigend.eq.1) then
         intout = intin
      else
         call mvbits( intin, 24, 8, intout, 0  )
         call mvbits( intin, 16, 8, intout, 8  )
         call mvbits( intin,  8, 8, intout, 16 )
         call mvbits( intin,  0, 8, intout, 24 )
      endif

      return
      end

c     ====================================================================      
      subroutine getrec(buff,buflen,isil,ierr)
      implicit real (a-h,o-z)
      character*4 buff(*)
      integer recln,currec,ierr,buflen
      common /rdrec/ iund,currec

      ierr = 0
      currec = currec + 1

      read(iund,rec=currec,err=100) intin

      call byter(intin,intout)

      recln = intout/4

      if (recln.gt.buflen) goto 100

      do i =1,recln
         read(iund,rec=currec+i,err=100) buff(i)
      end do
   
      currec = currec + recln
      currec = currec + 1

      read(iund,rec=currec,err=100) intin


      call byter(intin,intchk)

      if (intout.ne.intchk) then
         ierr = 1
         if (isil.eq.0) print*,'getrec: error reading file'
      endif

      return

100   ierr = 1
      return
      end

c     ====================================================================      
      subroutine bytr2(intin)
      implicit none
      integer bigend,nlen,i
      integer*2 intin(256),intout
      common /isend/ bigend,nlen

      do i=1,256
         if (bigend.ne.1) then
            call mvbits( intin(i),  8, 8, intout, 0 )
            call mvbits( intin(i),  0, 8, intout, 8 )
            intin(i) = intout
         endif
      end do

      return
      end

      integer*2 function i1to2(intin)
      integer*1 intin

      i1to2 = intin
      if (i1to2.lt.0) i1to2 = abs(intin) + 128
      return
      end

      subroutine getrc2(buff,ierr)
      implicit real (a-h,o-z)
      integer*2 buff(256)
      integer currec,ierr
      integer*2 into
      common /rdrec/ iund,currec

      ierr = 0
      currec = currec + 1

      do i =1,256
         read(iund,rec=currec+i,err=100) buff(i)
      end do
   
      call bytr2(buff)
      if (currec.eq.0) then
         if (buff(19).ne.100) goto 100
      endif

      currec = currec + 255

      return

100   ierr = 1
      return
      end

c     ====================================================================      
      subroutine fndor(idebug)
      implicit double precision ( a-h,o-z)
      character*137 line
      common /gauori/ nzm,nso,nio,nzo,ioropt,ifor,
     &                ixyz98,iopr,isymm,irc,imp2,icntp,itd

      if (idebug.eq.1) print*,'start find orientations'
      ifor = 0
      nzm = 0
      nso = 0
      nio = 0
      nzo = 0 
      istat = 1
      do while(istat.eq.1)
         call searchq(line,'Z-MATRIX (ANGSTROMS',
     &                     'Standard orientation:',
     &                     'Input orientation:',
     &                     'Z-Matrix orientation:',istat)
         if (icdex(line,'Z-MATRIX (ANGSTROMS').ne.0) nzm = nzm + 1
         if (icdex(line,'Standard orientation:').ne.0) nso = nso + 1
         if (icdex(line,'Input orientation:').ne.0) nio = nio + 1
         if (icdex(line,'Z-Matrix orientation:').ne.0) nzo = nzo + 1
      end do

      if (nso.eq.nio) ioropt = 1
      if (nso.eq.0.and.nio.eq.0) ioropt = 0
      if (nso.gt.nio) ioropt = 2
      if (nso.lt.nio) ioropt = 3

      if (idebug.eq.1) then
         print*,'nzm=',nzm,' nso=',nso,' nio=',nio,' nzo=',nzo
      endif

      return
      end

c     ====================================================================
      subroutine wrcubd(npts1,npts2,npts3,ipsi,denn)

c THIS IS REALLY wrcube

      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)
      common /grdhlp/ mx3d,mx3d2
      dimension v3(3),o(3),p(3),dtmp(6),denn(*)

      iun = 21

c     normalize the perpendicular vector

      rn = dsqrt(cx*cx+cy*cy+cz*cz)

      v3(1) = cx / rn
      v3(2) = cy / rn
      v3(3) = cz / rn

      p(1) = px
      p(2) = py
      p(3) = pz

c     calculate the origin for the cube

      do i=1,3
         o(i) = p(i) - 0.5d0*(v1(i)*r(1) + v2(i)*r(2) + v3(i)*r(3))
      end do

      write(iun,'(a)') 'Molden generated cube file'

      if (ipsi.eq.0) then
         write(iun,*) 'Density'
      else
         if (ipsi.lt.0) then
            write(iun,'(''Beta Orbital '',i4)') iabs(ipsi)
         else
            write(iun,'(''Orbital '',i4)') ipsi
         endif
      endif

c     number of atoms and origin of the cube

      write(iun,'(i5,4(f12.6))') natoms,(o(i),i=1,3)

c     number of points in each dimension, and voxel size

      write(iun,'(i5,4(f12.6))') npts1,
     &            (v1(i)*r(1)/dble(npts1-1),i=1,3)
      write(iun,'(i5,4(f12.6))') npts2,
     &            (v2(i)*r(2)/dble(npts2-1),i=1,3)
      write(iun,'(i5,4(f12.6))') npts3,
     &            (v3(i)*r(3)/dble(npts3-1),i=1,3)

c     atomic numbers and coordinates of the atoms (in bohrs)

      do i=1,natoms
        write(iun,'(i5,4(f12.6))') nat(i),0.0,(xyz(j,i),j=1,3)
      end do

c     density data

      ij = 0 

      do i=1,npts1
         do j=1,npts2
            ij = ij + 1
            kk = 0
            do while (kk.lt.npts3)
               n = 6
               if (npts3-kk.lt.n) n = npts3-kk
               do k=1,n
                  dtmp(k) = denn((npts3-(kk+k))*mx3d2 + ij)
                  if (dtmp(k).lt.0.0d0) then
                     if (dtmp(k).gt.-1.0d-99) dtmp(k) = -1.0d-99
                  else
                     if (dtmp(k).lt.1.0d-99) dtmp(k) = 1.0d-99
                  endif
               end do
               write(iun,'(6E13.5)') (dtmp(k),k=1,n)
               kk = kk + 6
            end do
         end do
      end do

      return
      end

c     ====================================================================      
      subroutine nmrshl
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      character*137 line
      common /nmr/    shlnuc(numatm),ihsnmr
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      call searchd(line,'GIAO Magnetic shielding',
     &                 'Magnetic shielding (ppm)',istat)
      if (istat.ne.0) then
          do i=1,natoms
             call readel(line,1)
             i1 = icdex(line,'Isotropic =')
             if (i1.ne.0) then
                 read(line(i1+11:i1+22),'(f11.4)') shlnuc(i)
             endif
             call readel(line,4)
          end do
          ihsnmr = 1

          call search(line,'Total nuclear spin-spin coupling J',istat)
          if (istat.ne.0) ihsnmr = 2
      else
          rewind iun2
          call searchd(line,'Total nuclear spin-spin coupling J',
     &                      'Fermi Contact (FC) contribution to J',
     &                      istat)
          if (istat.ne.0) ihsnmr = 2
      endif

      return
      end

c     ====================================================================      
      subroutine nmrcpd(idebug,couplj)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      character*137 line
      common /nmr/    shlnuc(numatm),ihsnmr
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension couplj(*)

      ntimes = natoms/5
      nt = mod(natoms,5)
      if (nt.ne.0) ntimes = ntimes + 1
      
      ndone = 0
      do l=1,ntimes
         call readel(line,1)
         nitems = 5
         ndo = ndone + nitems
         if (ndo.gt.natoms) nitems = nitems - (ndo - natoms)
         do i=(l-1)*5+1,natoms
            call readel(line,1)
            read(line,'(7x,5f14.10)') 
     &          (couplj((i-1)*natoms+(l-1)*5+j),j=1,nitems)
         end do
         ndone = ndone + nitems
      end do

      if (idebug.eq.1) then
         do i=1,natoms
            print*,'atom ',i,' ',(couplj((i-1)*natoms+j),j=1,i)
         end do
      endif

      return
      end

c     ====================================================================      
      subroutine qmtot
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /totchg/ itot

      cor = 0.0d0
      do i=1,natoms
         cor = cor + nat(i)
      end do
      itot = cor - nelecs

      return
      end

