      subroutine rdqchd(idebug,irtype,istats,
     &                 focc,focb,nocc,nocb,ncols,ncolb)

c THIS IS REALLY rdqchm

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
      dimension focc(*),focb(*)

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
      call srchmf(line,'$rem',istat)
      if (istat.eq.0) goto 1000

      do while (.true.)
          call rdmf(line,cline,istat)
          if (istat.eq.0) goto 10
          if (index(cline,'$END').ne.0) goto 10
          if (index(cline,' SP').ne.0) irtype = 1
          if (index(cline,' OPT').ne.0) irtype = 2
          if (index(cline,'FREQ').ne.0) irtype = 4
          if (index(cline,'PRINT_').ne.0) then
             if (index(cline,'ORBITALS').ne.0) then
                 if (index(cline,'TRUE').ne.0) ihbas = 1
             endif
             if (index(cline,'GENERAL_BASIS').ne.0) then
                 if (index(cline,'TRUE').ne.0) ihorb = 1
             endif
          endif
      end do

10    continue

      call rewmf

c     
c     look for number of orbitals and electrons
c

      call srchmf(line,'beta electrons',istat)
      if (istat.eq.0) then
         call inferr('no electrons line found!',1)
          goto 1000
      endif

      read(line,'(11x,i8,11x,i8)',err=1000,end=1000) neleca,nelecb

      nelecs = neleca + nelecb

      call srchmf(line,'basis functions',istat)
      if (istat.eq.0) then
         call inferr('no basis functions found!',1)
         goto 1000
      endif
      
      i1 = index(line,'basis functions')
      i2 = index(line,'and')
      line = line(i2+3:i1-1)
      read(line,*,err=1000,end=1000) norbs

      if (norbs.gt.mxorb) 
     &   call inferr('Exceeding MaxNum of Orbitals!',1)

      ncols  = norbs
      ncolb  = norbs
c
      do i=1,min0(mxorb,norbs)
         focc(i) = 0.0d0
      end do
      do i=1,min0(mxorb,neleca)
         focc(i) = 1.0d0
      end do

      call rewmf

      call srchmf(line,'Final Beta MO Coefficients',istat)
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
      call rewmf
      do while(.true.)
         call qcxyz(idebug,.true.,0,istat)
         if (istat.eq.0) goto 50
      end do
50    continue
      call xyzcoo(1,0,0)

      call rewmf
      if (ihbas.eq.1) then
         call rdbas(idebug,2,istat)
         call norml
         if (istat.eq.0) goto 1000
      endif

      if (idebug.eq.1) 
     &    call inferr('Succesfully read basis-set',0)
      if (idebug.eq.1) call basprt(iun3,.false.,.false.)
      if (norbs.gt.mxorb) goto 1000

      if (ihorb.eq.1) then

          call rdqvec(idebug,istat)
          if (istat.eq.0) goto 1000
          if (idebug.eq.1) call inferr('Succesfully read vectors',0)
       
      endif

c      call rewfil
c      call search(line,'Charges from ESP fit',istat)
c      if (istatio.eq.1) call search(line,'Charges from ESP fit',istat)
c      if (istat.eq.1) then
c          call redel(line,2)
c          do i=1,natbck
c              call nxtlin(line,jstat)
c              if (jstat.eq.1.or.jstat.eq.2) goto 1000
c              if (linlen(line).eq.18) then
c                 read(line,'(7x,f11.6)',err=1000,end=1000) qch(i)
c              else if (linlen(line).eq.21) then
c                 read(line,'(11x,f11.6)',err=1000,end=1000) qch(i)
c              endif
c          end do
c          ihasesp = 1
c      endif
c
c      call rewfil



      if (idebug.eq.1) 
     &         call inferr('Succesfully read QChem output',0)

      return

1000  if (idebug.eq.1) 
     &         call inferr('Error reading QChem output!',1)
      istats = 0
      return

      end

      subroutine prsqmd(molin)
      implicit double precision (a-h,o-z)
      common /mfdata/ nmols,imol,iendmf,ielin,maxmol,mollin
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line,caps
      character*40 title
      character*2 gstr
      dimension molin(*)
      
      call rewfil

      imol = 1
      molin(1) = 0
      nmols = 0
      ilin = 0


1        call nxtlin(line,jstat)
         if (jstat.eq.1) goto 100
         if (jstat.eq.2) goto 200

         caps = line
         len1  = len(line)
         call tocap(caps,len1)
         ilin = ilin + 1
         if (index(caps,'$COMMENT').ne.0) then
            call nxtlin(line,jstat)
            if (jstat.eq.1) goto 100
            if (jstat.eq.2) goto 200
            title = line(1:40)
            ilin = ilin + 1
         endif
         if (index(caps,'WELCOME TO Q-CHEM').eq.0) goto 1

         if (nmols.lt.maxmol) then
            nmols = nmols + 1
            molin(nmols) = ilin
            if (nmols.gt.1) call parsfn(title,linlen(title),13)
            title = 'molecule '//gstr(nmols)
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

      subroutine qcxyd(idebug,nuclear,ipnt,istat,
     &                 ianz,coo)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (MAXPNT=2000)
      parameter (maxfat=1000)
      integer getlin
      logical gnreal,nuclear
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

      if (nuclear) then
         call srchmf(lline,'Standard Nuclear Orientation',istats)
         if (istats.eq.0) goto 100
         call rdmf(lline,cline,istats)
         call rdmf(lline,cline,istats)
      else
         if (idebug.eq.1) print*,'ipnt=',ipnt,'lpnt=',lpnt(ipnt)
         if (lpnt(ipnt).ne.-1) then
            call putmf(lpnt(ipnt))
         else
            call rewmf
            do i=1,ipnt
               call srchmf(lline,'Optimization Cycle:',istats)
               if (istats.eq.0) goto 100
               lpnt(ipnt) = getmf()
            end do
         endif
         call rdmf(lline,cline,istats)
         call rdmf(lline,cline,istats)
         call rdmf(lline,cline,istats)
      endif

      i = 0
      if (idebug.eq.1)  print*,'coordinates'
      do while(.true.)
          if (getlin(1).eq.1) then
             if (icdex(line,'ATOM').eq.0) then
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.2) goto 50
                i = i + 1
                ianz(i) = 0
                ktype = nxtwrd(str,nstr,itype,rtype)
                tstr = '  '
                if (ktype.eq.1) then
                   if (nstr.eq.1) tstr(2:2) = str(1:1)
                   if (nstr.eq.2) tstr(1:2) = str(1:2)
                   do k=1,mxel
                     if (tocapf(tstr).eq.tocapf(elemnt(k))) then
                        if (nuclear) then
                           nat(i) = k
                        else 
                           ianz(i) = k
                        endif
                     endif
                   end do
                endif
                if (gnreal(r,3,.false.)) then
                   if (idebug.eq.1)  then
                      if (nuclear) then
                         print*,elemnt(nat(i)),(r(j),j=1,3)
                      else
                         print*,elemnt(ianz(i)),(r(j),j=1,3)
                      endif
                   endif
                   do j=1,3
                      if (nuclear) then
                         xyz(j,i) = r(j) / toang
                      else
                         coo(j,i) = r(j) / toang
                      endif
                   end do
                else
                   goto 100
                endif
             endif
          else 
             goto 100
          endif
      end do

50    if (nuclear) then
         natoms = i
      else
         iatoms = i
      endif

      if (.not.nuclear) then
          call srchmf(lline,'Cartesian Gradient',istats)
          if (istats.eq.0) goto 200
          if (getlin(1).ne.1) goto 200
          if (idebug.eq.1)  print*,'forces'
          do i=1,iatoms
             if (getlin(1).ne.1) goto 200
             ktype = nxtwrd(str,nstr,itype,rtype)
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (gnreal(r,3,.false.)) then
                if (idebug.eq.1)  print*,elemnt(ianz(i)),(r(j),j=1,3)
                do j=1,3
                   fxyz(j,i) = r(j)
                end do
             endif
          end do
      endif

      return

100   istat = 0
      call rewmf
      return
200   istat = -1
      return
      end

      subroutine rdqvcd(idebug,istats,
     &                  vectrs,vectrb,eiga,eigb,ncols,ncolb)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      character*137 line
      common /orbhlp/ mxorb,iuhf,ispd
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /curlin/ line
      integer getlin
      character*137 lline,str
      real eiga,eigb
      dimension vectrs(*),vectrb(*),eiga(*),eigb(*)

      istats = 1
      ido5d = 1
      ido7f = 1
      ido9g = 1

      call rewmf

      call srchmf(lline,'Final Alpha MO Eigenvalues',istat)
      if (istat.eq.0) goto 100

      ncols = 0
      nlines = 0
      do i=1,norbs

          if (getlin(1).ne.1) goto 100
          nlines = nlines + 1
          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.1) goto 10

          if (getlin(1).eq.1) then
             ktype = nxtwrd(str,nstr,itype,rtype)
             do j=1,6
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.3) then
                    ncols = ncols + 1
                    eiga(ncols) = real(rtype)
                elseif (ktype.eq.0) then
                    goto 10
                else
                    goto 100
                endif
             end do
          else
             goto 100
          endif
      end do
      
10    ncol = 0

      if (getlin(1).ne.1) goto 100
      do i=1,nlines
         if (getlin(1).ne.1) goto 100
         do j=1,norbs
             if (getlin(1).ne.1) goto 100
             ktype = nxtwrd(str,nstr,itype,rtype)
             ncolt = ncol  
             do k=1,6
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.3) then
                    ncolt = ncolt + 1
                    vectrs((ncolt-1)*mxorb+j) = rtype
                elseif (ktype.eq.0) then
                    goto 20
                else
                    goto 100
                endif
             end do
20           continue
         end do
         ncol = ncol + 6
      end do

      if (iuhf.eq.1) then
         call srchmf(lline,'Final Beta MO Eigenvalues',istat)
         if (istat.eq.0) goto 100
      endif

      return

100   istats = 0
      return
      end

      subroutine geoqcm(formax,forrms,dismax,disrms,epoints,isav)
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
      ifrmav = 0
      idmxav = 1
      idrmav = 0
      ieav = 1
      nepnts = 0
      ngeoms = 0
      dmaxt = 0.0d0
      fmaxt = 0.0d0
      drmst = 0.0d0
      frmst = 0.0d0

      do while(.true.)
         call srchmf(line,'Optimization Cycle:',istat)
         if (istat.eq.0) goto 100

         if (ngeoms.lt.mxpnt) then
            ngeoms = ngeoms + 1
         else
            goto 100
         endif

         call srchmf(line,'Energy is',istat)
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3) goto 100
         epoints(ngeoms) = rtype

         call srchmf(line,'Cnvgd?',istats)
         if (getlin(1).ne.1) goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3) goto 100
         formax(ngeoms) = rtype
         if (getlin(1).ne.1) goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3) goto 100
         dismax(ngeoms) = rtype
         isav(ngeoms) = 1
      end do

100   nepnts = ngeoms

      return
      end

      subroutine cnvqcm
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
         call srchmf(line,'Cycle       Energy',istat)
         if (istat.eq.0) goto 100
         if (getlin(1).ne.1) goto 100
         do while(.true.)
            if (ifrst.eq.1.and.jend1.ge.MAXITER) goto 10
            if (ifrst.eq.0.and.jend2.ge.MAXITER) goto 100
            if (getlin(1).ne.1) goto 100
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.1) goto 10
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.3) goto 100
            if (ifrst.eq.1) then
                jend1 = jend1 + 1
                convg1(jend1) = rtype
            else 
                jend2 = jend2 + 1
                convg2(jend2) = rtype
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

      subroutine getqfd(istat,coo)
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
      ihasi = 1
      call rewmf

      call iatnox(natoms)

      do i=1,natoms
        do j=1,3
           fcoo(j,i) = coo(j,i)
        end do
      end do

c     Read in Qchem Frequencies

      nvibs = 3*natoms - 6
      if (natoms.eq.1) nvibs = 0
      if (natoms.eq.2) nvibs = 1

      do while(.true.)
         call srchmf(line,'Frequency:',istat)
         if (istat.eq.0) goto 10

         ivibst = ivibs
         ktype = nxtwrd(str,nstr,itype,rtype)
         do j=1,3
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.3) goto 20
            ivibs = ivibs + 1
            freq(ivibs) = rtype
         end do

20       if (getlin(1).ne.1) goto 100
         if (getlin(1).ne.1) goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)
         ktype = nxtwrd(str,nstr,itype,rtype)
         do j=1,3
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.3) goto 30
            ivibst = ivibst + 1
            frint(ivibst) = rtype
         end do
30       continue

      end do

10    nfreq = ivibs
      call parptr(1,freq,freq,nfreq)
      call parptr(112,frint,ramint,ihasi)

      return

100   istat = 0
      return
      end 

      subroutine qcoord(idebug,ifreq,istat)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      character*2 gstr
      character*5 tstr
      real freqt
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      character*137 line,str
      common /curlin/ line
      integer getlin
      logical gnreal
      dimension r(3)

c     Get ifreq 'th Norm. Coordinates from Qchem Output

      istat = 1

      toang=0.52917706d0
      call rewmf

      call iatnox(iatoms)

      nvibs = nfreq
      ifrh = (ifreq-1) / 3
      ifrl = ifreq - ifrh*3
      do i=1,ifrh+1
         call srchmf(line,'Frequency:',istat)
         if (istat.eq.0) goto 100
      end do
      call srchmf(line,'X      Y      Z',istat)
      if (istat.eq.0) goto 100

      do i=1,iatoms
         if (getlin(1).ne.1) goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)
         do j=1,ifrl
            if (gnreal(r,3,.false.)) then
               do k=1,3
                   a(k,i) = r(k) / toang
               end do
            else
               goto 100
            endif
         end do
      end do


      if (idebug.eq.1) call prtfr(ifreq)
      return

100   istat = 0
      call inferr('Error reading Norm. Coords. !',0)
      return
      end
