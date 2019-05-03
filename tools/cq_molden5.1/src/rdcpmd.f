C This soubroutine reads in information from {}.out files 
C of CPMD (car-parrinello molecular dynamics code) of J.Hutter et al.
C This module has been written by TEODORO LAINO 
C NEST Laboratories - INFM - Scuola Normale Superiore - Sept 2002
C Distributed with the permission of TEODORO LAINO

      subroutine rdcpmddd(idebug,ibefo,istatio,ioxyz,irtype,ihsend,
     &                  istats,vectrs,vectrb,focc,focb,eiga,eigb,
     &                  nocc,nocb,ncols,ncolb)

c THIS IS REALLY cpmd subroutine interface to molden

      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs, nat(numatm)
      common /coord / xyz(3,numatm)
      common /orbhlp/ mxorb,iuhf,ispd
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line
      real eiga,eigb
      dimension vectrs(*),vectrb(*),focc(*),focb(*),eiga(*),eigb(*)

      istats = 1
      ioxyz = 1
      istatio = 0
      iuhf = 0

      if (idebug.eq.1) write(iun3,'(a)')'subroutine rdcpmd'
c 
c..  CPMD gives information on number of populated states
c..
      call search(line,'NUMBER OF STATES:',istat)
      if (istat.eq.1) then 
         read(line(50:),*)norbs
         call search(line,'NUMBER OF ELECTRONS',istat)
         read(line(50:),*) xelec
         nelecs=int(xelec)
      endif   
      rewind(iun2)

      if (norbs .gt. mxorb) then

c we'll never reach this condition with CPMD
c at least at the moment...
c Next thing to do is to read orbital information from restart file...

         call inferr('Exceeding MaxNum of Orbitals!',1)
         goto 1000
      endif

c==== Establish RUNTYPE  ====
C..
c.. ! irtype=3 ! MOLECULAR DYNAMICS
c.. ! irtype=2 ! OPTIMIZATION GEOMETRY
c.. ! irtype=1 ! OPTIMIZATION WAVEFUNCTION
c.. ! irtype=4 ! VIBRATIONAL ANALYSIS
c.. ! irtype=5 ! Reading just geometry ! Unknown flag
c..            ! Here you can add whatever you want to implement
c..
      call search(line,'MOLECULAR DYNAMICS',istat)
      if (istat.ne.0) then 
         irtype = 3

c In molecular dynamics you have to read coordinates from
c TRAJECTORY file and information from CPMD input file

      else 
         rewind(iun2)
         call search(line,'GEOMETRY OPTIMIZATION',istat)
         if (istat.ne.0) then 
            irtype = 2

c In geometry optimization the coordinates and info must be
c read from input CPMD file

         else 
            rewind(iun2)
            call search(line,'VIBRATIONAL ANALYSIS',istat)
            if (istat.ne.0) then 
               irtype = 4

c In vibrational analysis the normal coordinates must be read
c from MOLVIB file.

            else 
               rewind(iun2)
               call searchd(line,'WAVEFUNCTION OPTIMIZATION',
     *              'PROPERTIES',istat)
               if (istat.ne.0) then 
                  irtype = 1

c In wavefunction optimization (single point) or properties 
c calculation we read just molecular geometry and display it!!

               else
                  call inferr('Reading just Geometry!',1)                  
                  irtype = 5
               endif
            endif
         endif
      endif
      
      call search(line,'END OF GEOMETRY OPTIMIZATION',istat)
      if (istat.eq.1) istatio = 1

      rewind(iun2)
c
c==== First read coordinates 
c
      if (irtype.eq.1.or.irtype.eq.5) then
c
c ====== Stationary Point ================
c
         call rdcpmolu(istat)

         if (istat.eq.0) goto 1000 

         call xyzcoo(1,0,0)

      endif

      return

1000  istats = 0
      call inferr('ERROR reading CPMD output file!',1)
      return
      end


      subroutine rdcpmold(istats,coo,ianz)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (numatm=2000,mxel=100)
      common /moldat/ natoms, norbs, nelecs, nat(numatm)
      common /coord / xyz(3,numatm)
      common /athlp/  iatoms, mxnat
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /curlin/ line
      character*137 line,str
      integer getlin
      logical dobohr,docoo,excyc
      character elemnt*2
      common /elem/elemnt(mxel)
      dimension coo(3,*),ianz(*)
c
c====== read molecular geometry ==================
c
c     sline = search this line first
c     iemlin = then skip this many lines
c
      istats = 1
      toang = 0.52917706d0
      docoo = .false.
      dobohr = .false.

c      write(*,*)'debug: rdcpmolu!'
      call search(line,'********* ATOMS *********',istat)
      if (istat.eq.0) goto 1000
      call readel(line,1)

      if (docoo) then
         iatoms = 0
      else
         natoms = 0
      endif

      do while (getlin(0).eq.1)
         if (linlen(line).le.1) return
         if (index(line,'********').ne.0) return
         if (docoo) then
            iatoms = iatoms + 1
         else
            natoms = natoms + 1
         endif
         ktype = nxtwrd(str,nstr,itype,rtype) 
         if (ktype.ne.2) goto 1000
         ktype = nxtwrd(str,nstr,itype,rtype) 
         if (ktype.ne.1) goto 1000

         indexb = index(str(1:2),' ')
c         write(*,'(A,4I5)')str(1:2),ichar(str(1:1)),ichar(str(2:2)),
c     &        ichar(' '),indexb
         if (indexb.eq.2) str = ' '//str
         itype = 0
         excyc = .false.
         do while (.not.excyc)
            itype = itype + 1
c            write(*,'(2A)')str(1:2),elemnt(itype)
            if (index(str(1:2),elemnt(itype)).ne.0) excyc = .true.
            if (itype.gt.mxel) then
               write(*,'(A)')'Error determining atomic number!'
               goto 1000
            endif
         end do

c itype atomic number

         if (docoo) then
            ianz(iatoms) = itype
         else
            nat(natoms) = itype
         endif

         do i=1,3
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.3) goto 1000
            if (dobohr) rtype = rtype / toang
            if (docoo) then
               coo(i,iatoms) = rtype
            else
               xyz(i,natoms) = rtype
            endif
         end do
      end do

1000  istats = 0
      call inferr('CPMD: ERROR reading molecular geometry!',1)
      stop
      return
      end

      subroutine cpmdpd(ipnt,istat,coo,ianz)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      parameter (maxfat=1000)
      parameter (mxel=100)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat) 
      common /athlp/ iatoms, mxnat
      integer getlin
      logical excyc
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character*137 line,str
      common /curlin/ line
      character elemnt*2
      common /elem/elemnt(mxel)
      dimension coo(3,*),ianz(*)
c
c     Read Next Point in geometry optimization
c
      toang=0.52917706d0
      iatoms = 0
      jatoms = 0

      call search(line,
     &     'ATOM          COORDINATES                GRADIENTS',
     &     istat)

      if (istat.eq.0) return

      call haszm(.false.)

      do while (getlin(0).eq.1)
         if (linlen(line).le.1) return
         if (index(line,'********').ne.0) return
         iatoms=iatoms+1
         jatoms=jatoms+1
         ktype = nxtwrd(str,nstr,itype,rtype) 
         if (ktype.ne.2) goto 1000
         ktype = nxtwrd(str,nstr,itype,rtype) 
         if (ktype.ne.1) goto 1000

         indexb = index(str(1:2),' ')
         if (indexb.eq.2) str = ' '//str
         itype = 0
         excyc = .false.
         do while (.not.excyc)
            itype = itype + 1
            if (index(str(1:2),elemnt(itype)).ne.0) excyc = .true.
            if (itype.gt.mxel) then
               write(*,'(A)')'Error determining atomic number!'
               goto 1000
            endif
         end do

c itype atomic number

         ianz(iatoms) = itype

c read geometries...

         do i=1,3
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.3) goto 1000
            coo(i,iatoms) = rtype 
         end do

c read gradients...

         do i=1,3
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.3) goto 1000
            fxyz(i,jatoms) = rtype
         end do
      end do
      
      return

1000  istat = 0
      return

      end

      subroutine geocpmd(formax,forrms,dismax,disrms,epoints,isav,
     &                   istatc)
      implicit double precision (a-h,o-z)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      character lstr*137
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)

      istatc = -1
c
c     Read Geometry Convergence data
c
      rewind iun2
      igcvav = 1
      ifmxav = 1
      ifrmav = 1
      idmxav = 0
      idrmav = 0
      ieav = 1
      nepnts = 0
      ngeoms = 0
      dmaxt = 0.0d0
      fmaxt = 0.0d0
      drmst = 0.0d0
      frmst = 0.0d0

      i = 0
      do while (.true.)
         call search(lstr,' GEOMETRY STEP NR. ',istat)
         if (istat.eq.0) goto 100
         i = i + 1
         if (i.le.mxpnt) then
            ngeoms =  ngeoms + 1
            isav(i) = 1
            read(iun2,'(a)')lstr
            i1 = index(lstr,'ETOT=')
            epoints(i) = reada(lstr,i1+6,len(lstr))
            i1 = index(lstr,'GNMAX=')
            formax(i) = reada(lstr,i1+7,len(lstr))
            read(iun2,'(a)')lstr
            i1 = index(lstr,'GNORM=')
            forrms(i) = reada(lstr,i1+7,len(lstr))
c            write(*,*)epoints(i),formax(i),forrms(i)
         else
            istatc = i
         endif
      end do

100   continue
      if (ngeoms.eq.0) then
          igcvav = 0
          ifmxav = 0
          ifrmav = 0
          idmxav = 0
          idrmav = 0
          ieav = 0
      else
          nepnts = ngeoms
      endif

      return
      end

      subroutine cnvcpmd
      parameter (MAXITER=1000)
      implicit double precision (a-h,o-z)
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      common /convrg/ convg1(MAXITER),convg2(MAXITER),cnvmax,cnvmin,
     &jstrt1,jend1,jstrt2,jend2,icvav1,icvav2
      common /curlin/ line
      character*137 line,str
      integer getlin
      logical last,datlin

c
c     Read Scf Convergence data type... In cpmd we don't have
c     a real scf calculation but an optimization of the wavefunction
c     Here we read how energy and gradient of WF converges...
c
      rewind iun2
      icvav1 = 1
      icvav2 = 1

      call searchd(line,' NFI      GEMAX',' ETOT        DETOT',istat)
      if (istat.eq.0) goto 1000

      jstrt1 = 0

      do while (getlin(0).eq.1)
         if (linlen(line).le.1) goto 100
         indiis = index('ODIIS|',line)
         if (indiis.eq.0) then
            if (datlin(line)) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.ne.2) goto 1000
               jend1 = itype
               if (jend1.gt.MAXITER) then
                  write(*,'(a)')'Reached maximum number of iterations!'
                  goto 100
               endif
               do i=1,3
                  ktype = nxtwrd(str,nstr,itype,rtype)
                  if ((ktype.eq.3).and.(i.eq.3)) then
                     convg1(jend1) = rtype
                  endif
               end do
               if (jstrt1.eq.0) jstrt1 = jend1
            endif
         endif 
      end do

100   continue

      last = .false.

c we try to understand how many WFoptimization have been performed...

      NUMOPT = 0
      do while (istat.ne.0)
         NUMOPT = NUMOPT+1
         call search(line,' GEOMETRY STEP NR.',istat)
      end do
      
      NUMOPT = NUMOPT - 2
      IF (NUMOPT.lt.0) goto 2000 

      rewind iun2
      do i=1,NUMOPT
         call search(line,' GEOMETRY STEP NR.',istat)
      end do

      do i=1,4 
         read(iun2,'(a)')line
      end do

c     there could be somethingelse appended to file 
c     before wavefunction optimization cycle 
      call checkdummylines_scf(iun2,istat)

      if (istat.eq.0) goto 2000

      last = .true.
      jstrt2 = 0

      do while (getlin(0).eq.1)
         if (linlen(line).le.1) goto 2000
         indiis = index('ODIIS|',line)
         if (indiis.eq.0) then
            if (datlin(line)) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.ne.2) goto 2000
               jend2 = itype
               if (jend2.gt.MAXITER) then
                  write(*,'(A)')'Reached maximum number of iterations!'
                  goto 999
               endif
               do i=1,3
                  ktype = nxtwrd(str,nstr,itype,rtype)
                  if ((ktype.eq.3).and.(i.eq.3)) then
                     convg2(jend2) = rtype
                  endif
               end do
               if (jstrt2.eq.0) jstrt2 = jend2
            endif
         endif 
      end do

 999  continue

      return

1000  icvav1 = 0
      icvav2 = 0
      return

2000  if (.not.last) then
         icvav2 = 0
      endif

      return

      end

      subroutine cpmdgetfd(istat,coo)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      common /curlin/ line
      character*137 line,str
      integer getlin
      dimension coo(3,*)

      istat = 1

      ivibs = 0
      ihasi = 0
      idirct = 1
      nframe = 5
      iframe = 0
      rewind iun2

      call iatnox(natoms)

      do i=1,natoms
        do j=1,3
	  fcoo(j,i) = coo(j,i)
        end do
      end do

c     Read in CPMD  Frequencies

      call search(line,'PURIFICATION OF DYNAMICAL MATRIX',istat)
      if (istat.eq.0) goto 100
      call readel(line,4)
 10   if (getlin(0).eq.0) goto 100
      do i=1,4
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3.and.ktype.ne.0) then
            write(*,'(A,I6)')'Number of vibrational frequencies read: ',
     &           ivibs
            goto 200
         elseif (ktype.eq.0) then
            goto 100
         else
            ivibs = ivibs + 1
            freq(ivibs) = rtype 
         endif
      end do
20    continue
      goto 10

100   if (ivibs.eq.0) istat = 0

      nfreq = ivibs
      call parptr(1,freq,freq,nfreq)

      return

 200  call inferr('Error reading CPMD frequencies!',1)
      istat=0
      return
      end 

      subroutine cpmdcoorg(idebug,ifreq,istat)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /curlin/ line
      character*137 line,str
      integer getlin

c     Get ifreq 'th Norm. Cordinates from CPMD MOLVIB FILE...

      iunt = iun2
      iun2 = iun5
      ncol = 8
      istat = 1
      call iatnox(iatoms)
      ndim = iatoms*3
      rewind iun2

      call search(line,'>>>>>>> NEW SET',istat)
      if (istat.eq.0) goto 1000

c     columns have been stored in order of 8 eigenvectors...

      n1 = ifreq/ncol
      n2 = mod(ifreq,ncol)
      if (n2.eq.0) then 
         n1 = n1 - 1
         n2 = ncol
      endif

      numjump = ((iatoms*3)+3)*n1
      do i=1,numjump
         read(iun2,'(a)')line
      end do
c
      do i=1,3
         read(iun2,'(a)')line
      end do

c     put the pointer at the beginning of the data file...

      do i=1,iatoms
         do k=1,3
            if (getlin(0).eq.0) goto 1000
            do j=1,n2
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.ne.3) goto 1000
            enddo
            a(k,i) = rtype
         end do
      end do
     
      if (idebug.eq.1) call prtfr(ifreq)
      iun2 = iunt
      return

1000  istat = 0
      call inferr('Error reading Norm. Coords. in MOLVIB!',0)
      iun2 = iunt
      return
      end


      subroutine dyncpdd(ipoints,isav,epoints)
      implicit double precision (a-h,o-z)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      character lstr*137,line*137
      common /curlin/ line
      dimension epoints(*),isav(*)

c
c     Read Molecular Dynamics data
c
      rewind iun2
      ipoints = 0
      nread = 6
      igcvav = 1
      ieav = 1
      nepnts = 0
      ngeoms = 0
      ifmxav = 0
      ifrmav = 0
      idmxav = 0
      idrmav = 0

      call search(lstr,
     &     'NFI    EKINC   TEMPP           EKS      ECLASSIC    ',istat)
      
      do i=1,mxpnt
 10      read(iun2,'(a)',end=100)line
         leng=linlen(line)
         ktype = nxtwrd(lstr,nstr,itype,rtype)
         if (leng.eq.0) goto 100
         if (ktype.ne.2) goto 10
         if (i.ne.itype) goto 1000
         isav(i) = 1
         ngeoms =  ngeoms + 1
         do j=2,nread
            ktype = nxtwrd(lstr,nstr,itype,rtype)
         end do
         if (ktype.ne.3) goto 100
         epoints(i) = rtype
      end do

100   continue
      if (ngeoms.eq.0) then
          igcvav = 0
          ifmxav = 0
          ifrmav = 0
          idmxav = 0
          idrmav = 0
          ieav   = 0
      else
          nepnts = ngeoms
      endif

      ipoints = ngeoms
      return

 1000 call inferr('Error reading molecular dynamics data!',0)
      return

      end

      subroutine cpmdptdyd(ipnt,istat,coo)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      parameter (maxfat=1000)
      parameter (mxel=100)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat) 
      common /athlp/ iatoms, mxnat
      integer getlin
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character*137 line,str
      common /curlin/ line
      character elemnt*2
      common /elem/elemnt(mxel)
      dimension coo(3,*)
c
c     Read Next Point in geometry optimization
c
      toang=0.52917706d0
      rewind iun2

      call search(line,'NEW',istat)

      if (istat.eq.0) rewind iun2
      numlinetojump = iatoms*(ipnt-1)
      do i=1,numlinetojump
         read(iun2,'(A)')line
      end do

      call haszm(.false.)

      do i=1,iatoms
         if (getlin(0).eq.0) return
         if (linlen(line).le.1) return
         ktype = nxtwrd(str,nstr,itype,rtype) 
         if (ktype.ne.2) goto 1000

c read geometries...

         do j=1,3
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.3) goto 1000
            coo(j,i) = rtype 
         end do
      end do
      
      return

1000  istat = 0
      return
      end

      subroutine wrcpdd(iun,ianz,coo,
     &                  nat,icel,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (numat1=20000)
      parameter (mxel=100)
      parameter (mxsg=238)
      common /athlp/ iatoms, mxnat
      character*7 spnam
      common /spgrnm/ spnam(mxsg)
      character*2 elemnt,tocapf
      common /elem/elemnt(mxel)
      character*80 pseudopot,string
      dimension iatomdone(numat1),ianz(*),coo(3,*)

      toang = 0.52917706d0
      torad = datan(1.0d0) / 45.0d0

      do i=1,iatoms
         iatomdone(i)=0
      end do

      if (icel.ne.1) then
         call inferr('No Cell Data !',0)
         return
      endif

      call fdat(2,0,0,0,0,0)

      nstor = mxnat-nat
      pseudopot='_Pseudopotential KLEINMAN-BYLANDER'
      dalp = alpha
      dbet = beta
      dgam = gamma
c     dump input to file cpmd.inp
      write(iun,'(a/)') 
     &     '# Car-Parrinello Molecular Dynamics Input File. '
      write(iun,'(A)')'&CPMD'
      write(iun,'(A)')' MOLECULAR DYNAMICS'
      write(iun,'(A)')' QUENCH BO '
      write(iun,'(A)')' MAXSTEP'
      write(iun,'(A)')' 10'
      write(iun,'(A)')' EMASS'
      write(iun,'(A)')' 280.'
      write(iun,'(A)')' SCALED MASSES'
      write(iun,'(A)')' TIMESTEP '
      write(iun,'(A)')' 6.'
      write(iun,'(A)')' ISOLATED MOLECULE'
      write(iun,'(A/)')'&END'
      write(iun,'(A)')'&SYSTEM'
      write(iun,'(A)')' SYMMETRY'
      write(iun,'(A)')'  0'
      write(iun,'(A)')' CELL'
      write(iun,'(1x,3(f7.4,4x),3(f7.2,4x))') a/toang,b/a,c/a,
     &     dcos(dalp),dcos(dbet),dcos(dgam)
      write(iun,'(A)')' CHARGE'
      write(iun,'(A)')' 0 '
      write(iun,'(A)')' CUTOFF'
      write(iun,'(A)')' 60'
      write(iun,'(A/)')'&END'
      write(iun,'(A)')'&DFT'
      write(iun,'(A)')' FUNCTIONAL BLYP'
      write(iun,'(A)')' NEWCODE'
      write(iun,'(A/)')'&END'
      write(iun,'(A)')'&ATOMS'
c     scrittura della geometria molecolare....
      natoms=iatoms-8 

      do i=1,natoms
         if (iatomdone(i).eq.0) then
c     count how many elements of this kind we have...
            tocapf = elemnt(ianz(i))    
            numatkind = 1
            do j=i+1,natoms
               jj = ianz(j)
               if (tocapf(1:2).eq.elemnt(jj)(1:2)) numatkind=numatkind+1
            end do

c            print *,numatkind,natoms
c     now we write what we have to....

            indb = 1
            if (tocapf(1:1).eq.' ') indb = 2
            string = '*'//tocapf(indb:2)//pseudopot
            write(iun,'(a)')string
            if (tocapf.ne.' H') then
               write(iun,'(a)')' LMAX=P'
            else
               write(iun,'(a)')' LMAX=S'
            endif
            write(iun,'(i4)') numatkind
            do j=i,natoms
               jj=ianz(j)
               if (tocapf(1:2).eq.elemnt(jj)(1:2)) then
                  iatomdone(j)=1
                  write(iun,'(2x,3f12.6,i10)')(coo(k,j),k=1,3),j
               end if
            end do
            write(iun,'(a)')' ' 
         endif
      end do
      write(iun,'(a)')'&END'

      call inferr('Wrote file: cpmd.inp',0)      
      return
      end

      subroutine cpmdceld(nspg, 
     &                    xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (mxel=100)
      parameter (mxsg=238)
      character*137 line
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*7 spnam
      common /spgrnm/ spnam(mxsg)
      common /curlin/ line
      logical pbcon
      character*137 string

      pbcon=.true.
      toang   = 0.52917706d0
      toangle = 45.d0/datan(1.d0)
      istat = 1
      rewind iun2

c     In the future fix this approximation...

      nspg = 1
      write(*,'(2A)')'CPMD interface can, at the moment,  work',
     &     ' only with P1 point group!'

      call search(line,'CELL DIMENSION:',istat)
      if (istat.eq.0) goto 100

      ind1=index(line,'CELL DIMENSION:')
      line=line(ind1+16:)
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.eq.3) a = rtype
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.eq.3) b = rtype * a
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.eq.3) c = rtype * a
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.eq.3) alpha = acos(rtype)*toangle
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.eq.3) beta = acos(rtype)*toangle
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.eq.3) gamma = acos(rtype)*toangle

      call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)

      return

100   istat = 0
      call inferr('cpmdcell: error reading cell parameter.',0)
      return
      end

      subroutine checkdummylines_scf(iun2,istat)
      implicit double precision(a-h,o-z)
      common /curlin/ line
      character*137 line,string


 5    read(iun2,'(a)',err=100,end=100)line
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.ne.2) goto  5
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.ne.3) goto  5
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.ne.3) goto  5
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.ne.3) goto  5
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.ne.3) goto  5
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.ne.3) goto  5
      ktype = nxtwrd(string,nstr,itype,rtype)
      if (ktype.ne.0) goto  5

      backspace iun2
      istat = 1
      return

 100  continue
      istat = 0
      return
      end
