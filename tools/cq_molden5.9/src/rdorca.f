      subroutine rdorcd(idebug,irtype,istats,
     &                 ianz,iatoms,focc,focb,nocc,nocb,ncols,ncolb)

c THIS IS REALLY rdorca

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
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /oniomh/ ioni,nion,ionih(mxonh),natonh(mxonh),
     &                iomap(numatm),icntat(mxonh),fct(mxonh),
     &                xyzi(3,mxonh)
      common /pntlin/ lpnt(MAXPNT)
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
      call srctmf(line,'* Geometry Optimization Run *',
     &                 '* Single Point Calculation *',
     &                 '* Energy+Gradient Calculation *',istat)
      if (istat.eq.0) goto 1000
      if (index(line,'Single').ne.0) irtype = 1
      if (index(line,'Geometry').ne.0) irtype = 2
      if (index(line,'Gradient').ne.0) irtype = 1
      call srcdmf(line,'ORCA NUMERICAL FREQUENCIES',
     &                 'VIBRATIONAL FREQUENCIES',
     &                 istat)
      if (istat.ne.0) irtype = 4

c      call srchmf(line,'$rem',istat)
c      if (istat.eq.0) goto 1000
c
c      do while (.true.)
c          call rdmf(line,cline,istat)
c          if (istat.eq.0) goto 10
c          if (index(cline,'$END').ne.0) goto 10
c          if (index(cline,' SP').ne.0) irtype = 1
c          if (index(cline,' OPT').ne.0) irtype = 2
c          if (index(cline,'FREQ').ne.0) irtype = 4
c          if (index(cline,'PRINT_').ne.0) then
c             if (index(cline,'ORBITALS').ne.0) then
c                 if (index(cline,'TRUE').ne.0) ihbas = 1
c             endif
c             if (index(cline,'GENERAL_BASIS').ne.0) then
c                 if (index(cline,'TRUE').ne.0) ihorb = 1
c             endif
c          endif
c      end do
c
c10    continue
c
c      call rewmf
c
c     
c     look for number of orbitals and electrons
c
c
c      call srchmf(line,'beta electrons',istat)
c      if (istat.eq.0) then
c         call inferr('no electrons line found!',1)
c          goto 1000
c      endif
c
c      read(line,'(11x,i8,11x,i8)',err=1000,end=1000) neleca,nelecb
c
c      nelecs = neleca + nelecb
c
c      call srchmf(line,'basis functions',istat)
c      if (istat.eq.0) then
c         call inferr('no basis functions found!',1)
c         goto 1000
c      endif
c      
c      i1 = index(line,'basis functions')
c      i2 = index(line,'and')
c      line = line(i2+3:i1-1)
c      read(line,*,err=1000,end=1000) norbs
c
c      if (norbs.gt.mxorb) 
c     &   call inferr('Exceeding MaxNum of Orbitals!',1)
c
c      ncols  = norbs
c      ncolb  = norbs
cc
c      do i=1,min0(mxorb,norbs)
c         focc(i) = 0.0d0
c      end do
c      do i=1,min0(mxorb,neleca)
c         focc(i) = 1.0d0
c      end do
c
c      call rewmf
c
c      call srchmf(line,'Final Beta MO Coefficients',istat)
c      if (istat.ne.0) then
c         iuhf=1
c         nocc=neleca
c         nocb=nelecb
c         do i=1,min0(mxorb,norbs)
c            focb(i) = 0.0d0
c         end do
c         do i=1,min0(mxorb,nelecb)
c            focb(i) = 1.0d0
c         end do
c      else
c         nocc=max0(neleca,nelecb)
c         do i=1,min0(mxorb,nelecb)
c            focc(i) = focc(i) + 1.0d0
c         end do
c      endif
c
      call rewmf

      if (irtype.eq.1.or.irtype.eq.4) then
         ipnt = 1
         call orcxyz(idebug,ipnt,istat)
      elseif (irtype.eq.2) then
         ipnt = 1
         do while(.true.)
            call orcxyz(idebug,ipnt,istat)
            if (istat.eq.0) goto 50
            ipnt = ipnt + 1
         end do
50       continue
      endif

      call cooxyz(ianz,iatoms)

c      call rewmf
c      if (ihbas.eq.1) then
c         call rdbas(idebug,2,istat)
c         call norml
c         if (istat.eq.0) goto 1000
c      endif
c
c      if (idebug.eq.1) 
c     &    call inferr('Succesfully read basis-set',0)
c      if (idebug.eq.1) call basprt(iun3,.false.,.false.)
c      if (norbs.gt.mxorb) goto 1000
c
c      if (ihorb.eq.1) then
c
c          call rdqvec(idebug,istat)
c          if (istat.eq.0) goto 1000
c          if (idebug.eq.1) call inferr('Succesfully read vectors',0)
c       
c      endif

      if (idebug.eq.1) 
     &         call inferr('Succesfully read Orca output',0)

      return

1000  if (idebug.eq.1) 
     &         call inferr('Error reading Orca output!',1)
      istats = 0
      return

      end

      subroutine prsomd(molin)
      implicit double precision (a-h,o-z)
      common /mfdata/ nmols,imol,iendmf,ielin,maxmol,mollin
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line,caps
      character*40 title
      character*2 gstr
      dimension molin(*)
      
      call rewfil

      imol = 1
      molin(1) = 1
      nmols = 0
      ilin = 0


1        call nxtlin(line,jstat)
         if (jstat.eq.1) goto 100
         if (jstat.eq.2) goto 200
         caps = line
         len1  = len(line)
         call tocap(caps,len1)
         ilin = ilin + 1
         if (index(caps,'* O   R   C   A *').eq.0) goto 1

         if (nmols.lt.maxmol) then
            nmols = nmols + 1
            molin(nmols) = ilin
            title = 'molecule '//gstr(nmols)
            if (nmols.gt.1) call parsfn(title,linlen(title),13)
         endif

         goto 1


100   ielin = ilin
      if (nmols.eq.0) then
         nmols = 1
         title = 'molecule '//gstr(nmols)
      endif
      call parsfn(title,linlen(title),13)

      return

200   return
      end

      subroutine orcxyd(idebug,ipnt,istat,
     &                 ianz,coo)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (MAXPNT=2000)
      parameter (maxfat=1000)
      integer getlin
      logical gnreal
      character*137 line,lline,cline,str
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /athlp/  iatoms, mxnat
      common /coord / xyz(3,numatm)
      common /curlin/ line
      common /pntlin/ lpnt(MAXPNT)
      character*2 elemnt,tstr,tocapf
      common /elem/elemnt(mxel)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat)
      integer getmf
      dimension r(3),coo(3,*),ianz(*)

      toang=0.52917706d0
      istat = 1
      call haszm(.false.)

      if (idebug.eq.1) print*,'ipnt=',ipnt,' lpnt=',lpnt(ipnt)

      if (lpnt(ipnt).ne.-1) then
         call putmf(lpnt(ipnt))
      else
         call rewmf
         do i=1,ipnt
            if (ipnt.eq.1) then
               call srctmf(line,'OPTIMIZATION CYCLE',
     &                 '* Single Point Calculation *',
     &                 '* Energy+Gradient Calculation *',istats)
               if (istats.eq.0) goto 100
            else
               call srchmf(line,'OPTIMIZATION CYCLE',istats)
               if (istats.eq.0) goto 100
            endif
            call srchmf(lline,'CARTESIAN COORDINATES (ANGSTROEM)',
     &                            istats)
            if (istats.eq.0) goto 100
            lpnt(ipnt) = getmf()
         end do
      endif

      call rdmf(lline,cline,istats)

      i = 0
      if (idebug.eq.1)  print*,'coordinates'

      do while(.true.)
          if (getlin(1).eq.1) then
             if (line(1:4).ne.'----') then
                ktype = nxtwrd(str,nstr,itype,rtype)
                tstr = '  '
                if (ktype.eq.0) goto 50
                if (ktype.eq.1.and.nstr.le.2) then
                   i = i + 1
                   ianz(i) = 0
                   if (nstr.eq.1) tstr(2:2) = str(1:1)
                   if (nstr.eq.2) tstr(1:2) = str(1:2)
                   do k=1,mxel
                     if (tocapf(tstr).eq.tocapf(elemnt(k))) then
                           ianz(i) = k
                     endif
                   end do
                endif
                if (gnreal(r,3,.false.)) then
                   if (idebug.eq.1)  then
                       print*,elemnt(ianz(i)),(r(j),j=1,3)
                   endif
                   do j=1,3
                       coo(j,i) = r(j) / toang
                   end do
                else
                   goto 100
                endif
             endif
          else 
             goto 100
          endif
      end do

50    iatoms = i

      call srcdmf(lline,'CARTESIAN GRADIENT',
     &                  'INTERNAL COORDINATES (ANGSTROEM)',istats)
      if (istats.eq.0) goto 200

      if (icdex(lline,'INTERNAL').ne.0) then
         if (getlin(1).ne.1) goto 200
         call convzmat(coo,ianz,iatoms,2,0,0)
         call srchmf(lline,'CARTESIAN GRADIENT',istats)
         if (istats.eq.0) goto 200
      endif

      if (getlin(1).ne.1) goto 200
      if (idebug.eq.1)  print*,'forces'
      call rdmf(lline,cline,istats)

      do i=1,iatoms
         if (getlin(1).ne.1) goto 200
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (gnreal(r,3,.false.)) then
            if (idebug.eq.1)  print*,elemnt(ianz(i)),(r(j),j=1,3)
            do j=1,3
               fxyz(j,i) = r(j)
            end do
         endif
      end do

      return

100   istat = 0
      call rewmf
      return

200   istat = -1
      return
      end

c      subroutine rdqvcd(idebug,istats,
c     &                  vectrs,vectrb,eiga,eigb,ncols,ncolb)
c      implicit double precision (a-h,o-z)
c      parameter (numatm=2000)
c      character*137 line
c      common /orbhlp/ mxorb,iuhf,ispd
c      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
c      common /moldat/ natoms, norbs, nelecs,nat(numatm)
c      common /rdwr/   iun1,iun2,iun3,iun4,iun5
c      common /curlin/ line
c      integer getlin
c      character*137 lline,str
c      real eiga,eigb
c      dimension vectrs(*),vectrb(*),eiga(*),eigb(*)
c
c      istats = 1
c      ido5d = 1
c      ido7f = 1
c      ido9g = 1
c
c      call rewmf
c
c      call srchmf(lline,'Final Alpha MO Eigenvalues',istat)
c      if (istat.eq.0) goto 100
c
c      ncols = 0
c      nlines = 0
c      do i=1,norbs
c
c          if (getlin(1).ne.1) goto 100
c          nlines = nlines + 1
c          ktype = nxtwrd(str,nstr,itype,rtype)
c          if (ktype.eq.1) goto 10
c
c          if (getlin(1).eq.1) then
c             ktype = nxtwrd(str,nstr,itype,rtype)
c             do j=1,6
c                ktype = nxtwrd(str,nstr,itype,rtype)
c                if (ktype.eq.3) then
c                    ncols = ncols + 1
c                    eiga(ncols) = real(rtype)
c                elseif (ktype.eq.0) then
c                    goto 10
c                else
c                    goto 100
c                endif
c             end do
c          else
c             goto 100
c          endif
c      end do
c      
c10    ncol = 0
c
c      if (getlin(1).ne.1) goto 100
c      do i=1,nlines
c         if (getlin(1).ne.1) goto 100
c         do j=1,norbs
c             if (getlin(1).ne.1) goto 100
c             ktype = nxtwrd(str,nstr,itype,rtype)
c             ncolt = ncol  
c             do k=1,6
c                ktype = nxtwrd(str,nstr,itype,rtype)
c                if (ktype.eq.3) then
c                    ncolt = ncolt + 1
c                    vectrs((ncolt-1)*mxorb+j) = rtype
c                elseif (ktype.eq.0) then
c                    goto 20
c                else
c                    goto 100
c                endif
c             end do
c20           continue
c         end do
c         ncol = ncol + 6
c      end do
c
c      if (iuhf.eq.1) then
c         call srchmf(lline,'Final Beta MO Eigenvalues',istat)
c         if (istat.eq.0) goto 100
c      endif
c
c      return
c
c100   istats = 0
c      return
c      end

      subroutine geoorc(formax,forrms,dismax,disrms,epoints,isav)
      implicit double precision (a-h,o-z)
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      character*137 line,str
      common /curlin/ line
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
         call srchmf(line,'OPTIMIZATION CYCLE',istat)
         if (istat.eq.0) goto 100

         if (ngeoms.lt.mxpnt) then
            ngeoms = ngeoms + 1
         else
            goto 100
         endif

         call srchmf(line,'Total Energy       :',istat)
         if (istat.eq.0) goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3) goto 100
         epoints(ngeoms) = rtype

         call srchmf(line,'Geometry convergence',istats)
         if (istats.eq.0) goto 100
         if (getlin(1).ne.1) goto 100
         if (getlin(1).ne.1) goto 100

         if (getlin(1).ne.1) goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3) goto 100
         forrms(ngeoms) = rtype

         if (getlin(1).ne.1) goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3) goto 100
         formax(ngeoms) = rtype

         if (getlin(1).ne.1) goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3) goto 100
         disrms(ngeoms) = rtype

         if (getlin(1).ne.1) goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3) goto 100
         dismax(ngeoms) = rtype

         isav(ngeoms) = 1
      end do

100   nepnts = ngeoms

      return
      end

      subroutine cnvorc
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

      logical function zreado(nz,ianz,iz,bl,alpha,beta)
c
c      z-matrix reading routine for orca output.
c
      implicit double precision (a-h,o-z)
      character *2 iel,tstr
      character *2 tocapf
      character*137 line,str
      common /curlin/ line
      integer getlin
      logical gnreal, gnint

      dimension ianz(*), bl(*), alpha(*), beta(*)
      dimension iz(4,*), izz(3), r(3)
      dimension iel(100)

      data iel/'bq',
     $         'h ', 'he',
     $         'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne',
     $         'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar',
     $         'k ', 'ca',
     $                     'sc', 'ti', 'v ', 'cr', 'mn',
     $                     'fe', 'co', 'ni', 'cu', 'zn',
     $                     'ga', 'ge', 'as', 'se', 'br', 'kr',
     $ 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     $ 'in','sn','sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     $ 'pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     $ 'ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po',
     $ 'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm',
     $ 'bk','cf','x '  /

      zreado = .true.

      nz = 0

      do while (.true.)

         if (getlin(1).ne.1) goto 200

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.0) then
             return
         elseif (ktype.eq.1) then

             nz = nz + 1

             ianz(nz) = 0
             if (nstr.eq.1) tstr(2:2) = str(1:1)
             if (nstr.eq.2) tstr(1:2) = str(1:2)

             do i=1,100
                if (tocapf(tstr).eq.tocapf(iel(i))) ianz(nz) = i - 1
             end do

             if (gnint(izz,3,.false.)) then
                do i=1,3
                   iz(i,nz) = izz(i)
                end do
             else
                goto 200
             endif

             if (gnreal(r,3,.false.)) then
                bl(nz)    = r(1)
                alpha(nz) = r(2)
                beta(nz)  = r(3)
                do i=1,3
                   iz(i,nz) = izz(i)
                end do
             else
                goto 200
             endif

         endif

      end do

200   call inferr('error in z-matrix !',0)
      zreado = .false.

      return
      end

      subroutine getofd(istat,coo)
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
      call rewmf

      call iatnox(natoms)

      do i=1,natoms
        do j=1,3
           fcoo(j,i) = coo(j,i)
        end do
      end do

c     Read in Orca Frequencies

      nvibs = 3*natoms - 6
      if (natoms.eq.1) nvibs = 0
      if (natoms.eq.2) nvibs = 1

      call srchmf(line,'VIBRATIONAL FREQUENCIES',istat)
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
            frint(ivibs) = 0.0d0
         endif
      end do

      call srchmf(line,'IR SPECTRUM',istat)
      if (istat.eq.0) goto 10
      call redel(line,4)

      do i=1,natoms*3
         if (getlin(1).ne.1) goto 100
         if (linlen(line).le.1) goto 10
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.1) then
            l = index(str,':')
            if (l.gt.0) then
               ifrq  = dint(reada(str,1,l-1))
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.3) then
                  frint(ifrq+1) = rtype
                  ihasi = 1
               endif
            endif
         else
            goto 100
         endif
      end do

10    nfreq = ivibs
      call parptr(1,freq,freq,nfreq)
      call parptr(112,frint,ramint,ihasi)

      return

100   istat = 0
      return
      end 

      subroutine ocoord(idebug,ifreq,istat)
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

c     Get ifreq 'th Norm. Cordinates from Orca Output

      istat = 1

      call rewmf
      call iatnox(iatoms)
      nvibs = nfreq

      ifreqt = ifreq - 1
      ioff = ifreqt - ((ifreqt-1)/5)*5
      call srchmf(line,'NORMAL MODES',istat)
      if (istat.eq.0) goto 100

      call redel(line,6)

      ivibs = 0
      do while(.true.)
         call redel(line,1)
         if (ivibs.gt.nvibs) goto 20
         if (linlen(line).le.1) goto 20
         
         tstr = ' '//gstr(ifreqt)//' '

         if (index(line,tstr).ne.0) then
             do i=1,iatoms
                do j=1,3
                    call nxtlin(line,jstat)
                    if (jstat.eq.1.or.jstat.eq.2) goto 100
                    read(line,'(11x,6f11.6)',err=100,end=100) 
     &                   (freqt(k),k=1,6)
                    a(j,i) = freqt(ioff)
                end do
             end do
         else
             call redel(line,iatoms*3)
         endif
         ivibs = ivibs + 6
      end do

20    if (idebug.eq.1) call prtfr(ifreq)
      return

100   istat = 0
      call inferr('Error reading Norm. Coords. !',0)
      return
      end

