      subroutine rdnwcd(idebug,irtype,istats,
     &                 ianz,iatoms,focc,focb,nocc,nocb,ncols,ncolb)

c THIS IS REALLY rdnwch

      implicit double precision ( a-h,o-z)
      parameter (numat1=20000)
      parameter (mxonh=100)
      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (MAXPNT=2000)
      common /qmchar/ qch(numat1),ihasesp
      character*137 line,cline
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /orbhlp/ mxorb,iuhf,ispd
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /oniomh/ ioni,nion,ionih(mxonh),natonh(mxonh),
     &                iomap(numatm),icntat(mxonh),fct(mxonh),
     &                xyzi(3,mxonh)
      common /pntlin/ lpnt(MAXPNT)
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      dimension focc(*),focb(*),ianz(*)

      istats = 1
      irtype = 0
      istatio = 0
      toang=0.52917706d0
      irtype = 1
      ihbas = 0
      ihorb = 0
      ioni = 0
      do i=1,MAXPNT
         lpnt(i) = -1
      end do

      call rewmf

      call srcdmf(line,'XYZ format geometry','number of quantum atoms',
     &            istat)
      if (istat.ne.0) then
         if (icdex(line,'quantum').ne.0) then
            natoms = dint(reada(line,25,len(line)))
         else
            call readel(line,2)
            natoms = dint(reada(line,1,len(line)))
         endif
      endif

      call rewmf
      irtype = 1
      call srcdmf(line,'NWChem Geometry Optimization',
     &                 'NWChem Nuclear Hessian and Frequency',istat)
      if (istat.ne.0) then
         irtype = 2
         if (icdex(line,'Frequency').ne.0) irtype = 4
      else
         call rewmf
      endif

      if (irtype.eq.1) then
         ipnt = 1
         call nwxyz(idebug,ipnt,istat)
      elseif (irtype.eq.4) then
         ipnt = 1
         call fnwxyz(idebug,istat)
      elseif (irtype.eq.2) then
         ipnt = 1
         do while(.true.)
            call nwxyz(idebug,ipnt,istat)
            if (istat.eq.0) goto 50
            ipnt = ipnt + 1
         end do
50       continue
      endif

      call rewmf
      call srchmf(line,'Chemical Shielding Tensors',istat)
      if (istat.ne.0) then
         backspace iun2
         call nmcshl
      endif

      if (idebug.eq.1) 
     &         call inferr('Succesfully read NWchem output',0)

      return

1000  if (idebug.eq.1) 
     &         call inferr('Error reading NWchem output!',1)
      istats = 0
      return

      end

      subroutine nwxyd(idebug,ipnt,istat,
     &                 ianz,coo)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (MAXPNT=2000)
      integer getlin
      logical gnreal
      character*137 line,lline,str
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /athlp/  iatoms, mxnat
      common /coord / xyz(3,numatm)
      common /curlin/ line
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character*2 elemnt,tstr,tocapf
      common /elem/elemnt(mxel)
      common /pntlin/ lpnt(MAXPNT)
      common /pntcnv/ icv(MAXPNT)
      integer getmf
      dimension r(3),coo(3,*),ianz(*)

      toang=0.52917706d0
      istat = 1
      call haszm(.false.)


      if (lpnt(ipnt).ne.-1) then
         call putmf(lpnt(ipnt))
      else
         icv(ipnt) = 0
         call rewmf
         do i=1,ipnt
             call srchmf(line,'Output coordinates',istats)
             if (icdex(line,'angstroms').ne.0) icv(ipnt) = 1
             if (istats.eq.0) goto 100
             lpnt(ipnt) = getmf() + 1
         end do
      endif

      if (icv(ipnt).eq.0) toang = 1.0d0
      if (idebug.eq.1)  print*,'coordinates'

      do i=1,3
         call rdmf(line,lline,istat)
      end do

      do i=1,natoms
          if (getlin(1).eq.1) then
                ktype = nxtwrd(str,nstr,itype,rtype)
                ktype = nxtwrd(str,nstr,itype,rtype)
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.3) then
                   nat(i) = dint(rtype)
                   ianz(i) = nat(i)
                else
                   goto 100
                endif

                if (gnreal(r,3,.false.)) then
                   if (idebug.eq.1)  then
                       print*,elemnt(nat(i)),(r(j),j=1,3)
                   endif
                   do j=1,3
                       xyz(j,i) = r(j) / toang
                       coo(j,i) = xyz(j,i)
                   end do
                else
                   goto 100
                endif
          else 
             goto 100
          endif
      end do

      iatoms = natoms

      return

100   istat = 0
      call rewmf
      return

      end

      subroutine fnwxyd(idebug,istat,
     &                  ianz,coo)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (MAXPNT=2000)
      parameter (maxfat=1000)
      integer getlin
      logical gnreal
      character*137 line,lline,str
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /athlp/  iatoms, mxnat
      common /coord / xyz(3,numatm)
      common /curlin/ line
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character*2 elemnt,tstr,tocapf
      common /elem/elemnt(mxel)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat)
      integer getmf
      dimension r(3),coo(3,*),ianz(*)

      toang=0.52917706d0
      istat = 1
      call haszm(.false.)

      call rewmf
      call srchmf(line,'Atom information',istats)
      if (istats.eq.0) goto 100

      call readel(line,2)

      natoms = 0 
      do while (.true.)
          if (getlin(1).eq.1) then
                if (line(2:6).eq.'-----') goto 200

                natoms = natoms + 1
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.1) then
                   if (nstr.eq.1) then
                      tstr(2:2) = str(1:1)
                      tstr(1:1) = ' '
                   elseif (nstr.eq.2) then
                      tstr(1:2) = str(1:2)
                   endif
                   do k=1,mxel
                      if (tocapf(tstr).eq.tocapf(elemnt(k))) then
                         nat(natoms) = k
                      endif
                   end do
                   ianz(natoms) = nat(natoms)
                else
                   goto 100
                endif

                ktype = nxtwrd(str,nstr,itype,rtype)

                if (gnreal(r,3,.false.)) then
                   if (idebug.eq.1)  then
                       print*,elemnt(nat(natoms)),(r(j),j=1,3)
                   endif
                   do j=1,3
                       xyz(j,natoms) = r(j) / toang
                       coo(j,natoms) = xyz(j,natoms) * toang
                   end do
                else
                   goto 100
                endif
          else 
             goto 100
          endif
      end do

200   continue

      iatoms = natoms

      return

100   istat = 0
      call rewmf
      return

      end

      subroutine geonwc(formax,forrms,dismax,disrms,epoints,isav)
      implicit double precision (a-h,o-z)
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      character*137 line,str
      common /curlin/ line
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      integer getlin
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)

      call rewmf

      igcvav = 1
      ifmxav = 1
      ifrmav = 1
      idmxav = 1
      idrmav = 1
      ieav = 1
      nepnts = 0
      ngeoms = 0
      dmaxt = 0.0d0
      fmaxt = 0.0d0
      drmst = 0.0d0
      frmst = 0.0d0

      do while(.true.)

         call srcdmf(line,'@ Step       Energy',
     &                    '  Step       Energy',istats)
         if (istats.ne.1) goto 100
         call readel(line,1)
         if (getlin(1).ne.1) goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            if (itype.lt.mxpnt) then
               ngeoms = itype + 1
            else
            endif
         else
            goto 100
         endif

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.3) then
            epoints(ngeoms) = rtype
         else
            goto 100
         endif
         
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.3) then
            formax(ngeoms) = rtype
         else
            goto 100
         endif
         
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.3) then
            forrms(ngeoms) = rtype
         else
            goto 100
         endif
         
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.3) then
            disrms(ngeoms) = rtype
         else
            goto 100
         endif
         
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.3) then
            dismax(ngeoms) = rtype
         else
            goto 100
         endif
         
         isav(ngeoms) = 1
      end do

100   nepnts = ngeoms

      return
      end

      subroutine cnvnwc
      parameter (MAXITER=1000)
      implicit double precision (a-h,o-z)
      common /convrg/ convg1(MAXITER),convg2(MAXITER),cnvmax,cnvmin,
     &jstrt1,jend1,jstrt2,jend2,icvav1,icvav2
      character*137 line,str
      common /curlin/ line
      integer getlin

      call rewmf

      icvav1 = 0
      icvav2 = 0
      jstrt1 = 1
      jstrt2 = 1
      jend1 = 0
      jend2 = 0
      ifrst = 1


c     SCF convergence 

      do while(.true.)
         call srchmf(line,'SCF ITERATIONS',istat)
         if (istat.eq.0) goto 100
         if (getlin(1).ne.1) goto 100
         if (getlin(1).ne.1) goto 100
         do while(.true.)
            if (ifrst.eq.1.and.jend1.ge.MAXITER) goto 10
            if (ifrst.eq.0.and.jend2.ge.MAXITER) goto 100
            if (getlin(1).ne.1) goto 100
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.0) then
                goto 10
            else
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.3) then
                  if (ifrst.eq.1) then
                      jend1 = jend1 + 1
                      convg1(jend1) = rtype
                  else 
                      jend2 = jend2 + 1
                      convg2(jend2) = rtype
                  endif
               endif
            endif
         end do
10       if (ifrst.eq.1) then
            ifrst = 0
            icvav1 = 1
         else
            icvav2 = 1
         endif
      end do

100   return
      end

      subroutine getnfd(istat,coo)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /athlp/ iatoms, mxnat
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      character*137 line,str
      common /curlin/ line
      integer getlin
      dimension coo(3,*)

      istat = 1

      ivibs = 0
      ihasi = 0
      idirct = 1
      nframe = 5
      iframe = 0
      call rewmf

      call iatnox(natoms)

      do i=1,natoms
        do j=1,3
           fcoo(j,i) = coo(j,i)
        end do
      end do

c     Read in NWchem Frequencies

      nvibs = 3*natoms - 6
      if (natoms.eq.1) nvibs = 0
      if (natoms.eq.2) nvibs = 1

      call srchmf(line,'Eckart algorithm',istat)
      if (istat.eq.0) goto 10
      call srchmf(line,'NORMAL MODE EIGENVECTORS',istat)
      if (istat.eq.0) goto 10

      call srchmf(line,'Infra Red Intensities',istat)
      if (istat.eq.0) goto 10
      if (getlin(1).ne.1) goto 100
      if (getlin(1).ne.1) goto 100

      do i=1,natoms*3
         if (getlin(1).ne.1) goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.3) then
            ivibs = ivibs + 1
            freq(ivibs) = rtype
         else
            goto 100
         endif
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.3) then
            ihasi = 1
            frint(ivibs) = rtype
         else
            goto 100
         endif
      end do

10    nfreq = ivibs
      call parptr(1,freq,freq,nfreq)
      call parptr(112,frint,ramint,ihasi)

      do i=1,6
         frint(i) = 0.0d0
      end do

      return

100   istat = 0
      return
      end 

      subroutine nwcord(idebug,ifreq,istat)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      character*137 lstr
      character*2 gstr
      character*4 tstr
      real freqt
      dimension freqt(9)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      character*137 line,str
      common /curlin/ line
      integer getlin
      logical gnreal
      dimension r(3)

c     Get ifreq 'th Norm. Cordinates from NWchem Output

      istat = 1

      call rewmf
      call iatnox(iatoms)
      nvibs = nfreq

      ioff = ifreq - ((ifreq-1)/6)*6
      call srchmf(line,'Eckart algorithm',istat)
      if (istat.eq.0) goto 100
      call srchmf(line,'NORMAL MODE EIGENVECTORS',istat)
      if (istat.eq.0) goto 100

      call readel(line,3)

      ivibs = 0
      do while(.true.)
         call readel(line,1)
         if (ivibs.gt.nvibs) goto 20
         if (linlen(line).le.1) goto 20
         
         tstr = ' '//gstr(ifreq)//' '

         if (index(line,tstr).ne.0) then
             call readel(line,3)
             do i=1,iatoms
                do j=1,3
                    read(iun2,'(a)') line
                    read(line,'(12x,6f12.5)',err=100,end=100) 
     &                   (freqt(k),k=1,6)
                    a(j,i) = freqt(ioff)
                end do
             end do
         else
             call readel(line,3)
             call readel(line,iatoms*3)
         endif
         call readel(line,1)
         ivibs = ivibs + 6
      end do

20    if (idebug.eq.1) call prtfr(ifreq)
      return

100   istat = 0
      call inferr('Error reading Norm. Coords. !',0)
      return
      end

      subroutine nmcshd(couplj)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      character*137 line
      common /nmr/    shlnuc(numatm),ihsnmr
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension couplj(*)

      call rewmf
      call srchmf(line,'Chemical Shielding Tensors',istat)
      if (istat.ne.0) then
          do i=1,natoms
             call srchmf(line,' isotropic =',istat)
             if (istat.ne.0) then
                 i1 = icdex(line,'=')
                 if (i1.ne.0) then
                    shlnuc(i) = reada(line,i1+1,len(line))
                 endif
             endif
          end do
          ihsnmr = 1

          call search(line,'Indirect Spin-Spin Tensors',istat)
          if (istat.ne.0) then
             ihsnmr = 2
             do i=1,natoms
                do j=1,natoms
                   if (i.ne.j) then
                      call search(line,
     &                   'Isotropic Spin-Spin Coupling',istat)
                      if (istat.ne.0) then
                         i1 = icdex(line,'=')
                         if (i1.ne.0) then
                            couplj((i-1)*natoms + j) = 
     &                            reada(line,i1+1,len(line))
                         endif
                      endif
                   endif
                end do
             end do
          endif
      endif

      return
      end

