      subroutine rdgadd(idebug,ibefo,istatio,irtype,ihsend,
     &                  istats,vectrs,vectrb,focc,focb,eiga,eigb,
     &                  nocc,nocb,ncols,ncolb)

c THIS IS REALLY rdgamu

      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs, nat(numatm)
      common /coord / xyz(3,numatm)
      common /orbhlp/ mxorb,iuhf,ispd
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line
      real eiga,eigb
      dimension vectrs(*),vectrb(*),focc(*),focb(*),eiga(*),eigb(*)
      parameter (one=1.0d+00, two=2.0d+00)

      i0 = 0
      i1 = 1
      i2 = 2
      istats = 1
      istatio = 0
      iuhf = 0

      if (idebug.eq.1) write(iun3,'(a)')'subroutine rdgamu'

c
c==== Read # of basis functions, electrons
c
C**      call search(line,'NUMBER OF BASIS FUNCTIONS',istat)
c.. New GAMESS can handle both cartesian and polar coords
c.. However, it prints eigenvectors in terms of 6d
c..
      call searchd(line,'NUMBER OF CARTESIAN GAUSSIAN BASIS',
     &                  'NUMBER OF BASIS FUNCTIONS',istat)
      if (istat.eq.1) then 
         if (icdex(line,'CARTESIAN GAUSSIAN').ne.0) then
c
c This is a new version of GAMESS.  Lines are wider
c
            read(line,'(47x,i5)') norbs
            call search(line,'NUMBER OF ELECTRONS',istat)
            read(line,'(47x,i5)') nelecs
            call search(line,'NUMBER OF OCCUPIED ORBITALS (ALPHA',istat)
            read(line,'(47x,i5)') nalpha
c           write(*,*) ' ... nalpha = ',nalpha
            call search(line,'NUMBER OF OCCUPIED ORBITALS (BETA ',istat)
            read(line,'(47x,i5)') nbeta
c           write(*,*) ' ... nbeta  = ',nbeta 
            read(line,'(47x,i5)') nbeta

c       ... check if PP used
            call search(line,'PP    =',istat)
c           write(*,*) line(9:12)
            if (line(9:12).eq.'NONE') then
              is_pp=0
            else
              is_pp=1
            endif

c       ... read PP valence electrons
            if (is_pp.eq.1) then
              call search(line,'NUMBER OF ELECTRONS KEPT IN THE ',istat)
              read(line,'(50x,i4)') nelecs
              call readel(line,1)
              read(line,'(50x,i4)') nalpha
              call readel(line,1)
              read(line,'(50x,i4)') nbeta
c             write(*,*) nelecs,nalpha,nbeta
            endif
c
c======= Occupy Orbitals - just in case, done here
c
            do i=1,norbs
              focc(i) = 0.0d0
            end do

            if (nalpha.eq.nbeta) then

              do mel=1,nalpha
                focc(mel) = two
              end do

            else

              do mel=1,nalpha
                focc(mel) = one
              end do

              do mel=1,nbeta
                focb(mel) = one
              end do

            endif
         else 
c
c.. Old GAMESS (Pre-Jan10 2000)
c
            read(line,'(38x,i5)') norbs
            read(iun2,'(38x,i5)') nelecs
         endif   
      endif   
      rewind(iun2)

      if (norbs .gt. mxorb) then
         call inferr('Exceeding MaxNum of Orbitals!',1)
         goto 1000
      endif

c==== Establish Scftype and Runtype ====

      irtype = 1

      call search(line,'CONTRL OPTIONS',istat)
      call search(line,'SCFTYP=',istat)
      if (index(line,'UHF').ne.0) iuhf=1
c      call search(line,'RUNTYP=',istat)
      if (index(line,'HESSIAN').ne.0) irtype=4

      call search(line,'STATIONARY POINT LOCATION RUN',istat)
      if (istat.eq.1) then
          call search(line,'HSSEND =',istat)
          if (istat.eq.1) then
              if (index(line,' T').ne.0) then
                 if (ihsend.eq.1) then
                    irtype=4
                 else
                    call inferr(
     &               'Use commandlineflag -H for normal modes',0)
                 endif
              endif
          endif
          call search(line,'EQUILIBRIUM GEOMETRY LOCATED',istat)
          if (istat.eq.1) istatio = 1
          if (irtype.ne.4) irtype = 2
      endif

      rewind iun2
c
c==== First read coordinates, and vectors
c
      if (istatio.eq.1.and.ibefo.eq.0) then
c
c ====== Stationary Point ================
c
         call search(line,'EQUILIBRIUM GEOMETRY LOCATED',istat)

         call rdmolu(i1,i2,i0,i1,istat)
         call searcht(line,'MOLECULAR ORBITALS',
     &                     '-MCHF- NATURAL ORBITALS',
     &                     'MCSCF NATURAL ORBITALS',istat)
c     &                     'NATURAL ORBITALS IN ATOMIC', istat)
         if (iuhf.eq.1) then
            call readel(line,3)
         else
            call readel(line,1)
         endif
         call rdvecu(vectrs,eiga,focc,norbs,ncols,istats,.false.,
     &               .false.)

         if (iuhf.eq.1) then
            backspace iun2
            backspace iun2
            backspace iun2
            call search(line,'BETA SET',istat)
            call rdvecu(vectrb,eigb,focc,norbs,ncolb,istats,.false.,
     &                  .false.)
         endif

         call inferr('Using Density of stationary point',0)

      else
c
c ====== First Point ===========
c
         call rdmolu(i2,i1,i0,i0,istat)
c STRIPPED a space character because PC-UNIX output misses first column
c         call searcht(line,'          EIGENVECTORS',
         call searcht(line,'         EIGENVECTORS',
     &                     '-MCHF- NATURAL ORBITALS',
     &                     'MCSCF NATURAL ORBITALS',istat)
c     &                     'NATURAL ORBITALS IN ATOMIC',istat)
         call readel(line,1)
         call rdvecu(vectrs,eiga,focc,norbs,ncols,istats,.false.,
     &               .false.)
         if (iuhf.eq.1) then
            backspace iun2
            backspace iun2
c            call search(line,'          EIGENVECTORS',istat)
            call search(line,'         EIGENVECTORS',istat)
            call readel(line,1)
            call rdvecu(vectrb,eigb,focc,norbs,ncolb,istats,.false.,
     &                  .false.)
         endif
         call inferr('Using Density of first point',0)

      endif

      rewind iun2
c
c======= Read Basis Set ===========
c
c     mdebug=idebug
c     idebug=1
      call rdbasg(idebug,.true.,istat)
c     idebug=mdebug

      if (istat.eq.0) goto 1000

      call search(line,'MULLIKEN AND LOWDIN POPULATION ANALYSES',
     &            istat)

c     Alpha

      call search(line,'MULLIKEN ATOMIC POPULATION',istat)
      if (iuhf.eq.1) call readel(line,1)
      call rdpopu(focc,natoms,istats)

c     Beta

      if (iuhf.eq.1) then
         rewind iun2
         call search(line,'MULLIKEN ATOMIC POPULATION',istat)
         call search(line,'MULLIKEN ATOMIC POPULATION',istat)
         call readel(line,1)
         call rdpopu(focb,natoms,istats)
      endif

      call xyzcoo(1,0,0)

c>>> Added Pipek-Mezey and Edmiston-Ruedenberg  FPA 3/15/2000

      rewind iun2
      call searcht(line,'BOYS ORBITAL LOCALIZATION',
     &                  'LOCALIZED BY THE POPULATION',
     &                  'EDMISTON-RUEDENBERG',istat)
      if (istat.ne.0) then

          if (index(line,'BOYS ORBITAL LOCALIZATION').ne.0) then
              print*,'***** READING BOYS  LOCALIZED ORBITALS *****'
              call search(line,'THE BOYS LOCALIZED ORBITALS ARE',
     &                    istat)
              call rdvecu(vectrs,eiga,focc,norbs,ncols,istats,.true.,
     &                    .false.)

              if (iuhf.eq.1) then
                call search(line,'THE BOYS LOCALIZED ORBITALS ARE',
     &                      istat)
                call rdvecu(vectrb,eigb,focc,norbs,ncolb,istats,.true.,
     &                      .false.)
              endif

          elseif (index(line,'LOCALIZED BY THE POPULATION').ne.0) then

              print*,
     &            '***** READING Pipek-Mezey LOCALIZED ORBITALS *****'
              call search(line,'POPULATION LOCALIZED ORBITALS ARE',
     &                    istat)
              call rdvecu(vectrs,eiga,focc,norbs,ncols,istats,.true.,
     &                    .false.)

              if (iuhf.eq.1) then
                call search(line,'POPULATION LOCALIZED ORBITALS ARE',
     &                     istat)
                call rdvecu(vectrb,eigb,focc,norbs,ncolb,istats,.true.,
     &                      .false.)
              endif

          elseif (index(line,'EDMISTON-RUEDENBERG').ne.0) then

              print*,
     &          '** READING Edmiston-Ruedenberg LOCALIZED ORBITALS**'
              call search(line,'ENERGY LOCALIZED ORBITALS',istat)
              call rdvecu(vectrs,eiga,focc,norbs,ncols,istats,.true.,
     &                    .false.)

              if (iuhf.eq.1) then
                call search(line,'ENERGY LOCALIZED ORBITALS',istat)
                call rdvecu(vectrb,eigb,focc,norbs,ncolb,istats,.true.,
     &                      .false.)
              endif
          endif

c no orbital energies associated with localised orbitals

          do i=1,norbs
              eiga(i) = 0.0e0
              eigb(i) = 0.0e0
          end do

      endif

c-------- Try for GUGA Natural Orbitals

      rewind iun2
      call search(line,'NATURAL ORBITALS IN ATOMIC',istat)
      if (istat.ne.0) then
          call rdvecu(vectrs,eiga,focc,norbs,ncols,istats,.false.,
     &                .true.)
c no orbital energies associated with localised orbitals

          do i=1,norbs
              eiga(i) = 0.0e0
              eigb(i) = 0.0e0
          end do

      endif

      rewind iun2
c      call gamunmr

      return

1000  istats = 0
      call inferr('ERROR reading GAMESS output file!',1)
      return
      end

      subroutine rdvecu(v,e,focc,norbs,ncol,istats,boys,nos)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (mxcol=10)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /curlin/ line
      common /orbhlp/ mxorb,iuhf,ispd
      character*137 line,str
      integer getlin
      real e
      logical boys,reo,nos
      dimension v(*), e(*), focc(*)
      dimension vt(6,mxcol),irord(6)
c
c====== read vectors ==================
c
      istats = 1

      reo = .false.
      irord(1) = 3
      irord(2) = 1
      irord(3) = 2
      irord(4) = 5
      irord(5) = 6
      irord(6) = 4

      lmin = 0
      do while (.true.)

         if (boys) then
             call readel(line,1)
         else
             if (nos) then
                call readel(line,3)
             else
                call readel(line,2)
             endif
         endif

c
c    Get Eigenvalues
c
         if (getlin(0).ne.1) goto 1000
         if (linlen(line).le.1) then
             if (getlin(0).ne.1) goto 1000
         endif
         nc = 0
         do i=1,mxcol
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.0) goto 5
             if (ktype.eq.3) then
                nc = nc + 1
                if (nos) then
                   focc(lmin+i) = rtype
                else
                   e(lmin+i) = rtype
                endif
             elseif (ktype.eq.2.and.boys) then
                nc = nc + 1
             else
                goto 100
             endif
         end do

5        continue
         call readel(line,1)
         
c
c    Get Eigenvectors
c
         istrt = 0
         itel = 0

         do j=1,norbs

            if (getlin(0).ne.1) goto 1000
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.2) goto 100
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.1) goto 100
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.1.and.ktype.ne.2) goto 100
            if (ktype.eq.2) then
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.1) goto 100
            endif
            if (istrt.eq.0.and.nstr.eq.3) then
                if (str(1:3).eq.'XXY') then
                   istrt = j
                   itel = 0
                   reo = .true.
                endif
            endif

            do i=1,mxcol
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.0) goto 10
                if (ktype.ne.3) goto 100
                v((lmin+i-1)*mxorb+j) = rtype
            end do
          
10          continue

           if (istrt.ne.0) then
              itel = itel + 1
              do k=1,nc
                 vt(itel,k) = v((lmin+k-1)*mxorb+j)
              end do
           endif

           if (itel.eq.6) then
              do l=1,6
                 do k=1,nc
                    v((lmin+k-1)*mxorb + (istrt+l-1)) = vt(irord(l),k)
                 end do
              end do
              istrt = 0
              itel = 0
           endif

         end do

         lmin = lmin + nc
      end do

100   continue
      ncol = lmin

      if (reo) then
          print*,'========================================'
          print*,'Changed order of F functions:'
          print*,' '
          print*,'xxy,xxz,yyx,yyz,zzx,zzy ->'
          print*,'xyy,xxy,xxz,xzz,yzz,yyz'
          print*,'========================================'
      endif

      return

1000  istats = 0
      call inferr('ERROR reading vectors!',1)
      return
      end

      subroutine rdpopu(focc,natoms,istats)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (mxcol=10)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /curlin/ line
      character*137 line,str
      integer getlin
      dimension focc(*)
c
c====== read populations ==================
c
      istats = 1

      
      lmin = 0
      do while (.true.)

         call readel(line,3)

c
c    Get Populations
c
         if (getlin(0).ne.1) goto 1000
         nc = 0
         do i=1,mxcol
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.0) goto 5
             if (ktype.ne.3) goto 100
             nc = nc + 1
             focc(lmin+i) = rtype
         end do

5        continue
         call readel(line,1)
         
c
c    Dummy reads
c
         do j=1,natoms

            if (getlin(0).ne.1) goto 1000
          
         end do

         lmin = lmin + nc
      end do

100   continue
      ncol = lmin

      return

1000  istats = 0
      call inferr('ERROR reading Populations!',1)
      return
      end

      subroutine rdmodu(ilino,iemlin,idocoo,idobohr,istats,
     &                  coo,ianz)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs, nat(numatm)
      common /coord / xyz(3,numatm)
      common /athlp/  iatoms, mxnat
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /curlin/ line
      character*137 line,str
      integer getlin
      dimension coo(3,*),ianz(*)
c
c====== read molecular geometry ==================
c
c     sline = search this line first
c     iemlin = then skip this many lines
c
      istats = 1
      toang = 0.52917706d0

      if (ilino.eq.1) then
         call search(line,
     &        'COORDINATES OF ALL ATOMS ARE (ANGS)',istat)
      elseif (ilino.eq.2) then
         call search(line,
     &        'ATOM      ATOMIC                      COORDINATES',
     &             istat)
      endif

      if (istat.eq.0) goto 1000
      call readel(line,iemlin)

      if (idocoo.eq.1) then
         iatoms = 0
      else
         natoms = 0
      endif

      do while (getlin(0).eq.1)
         if (linlen(line).le.1) return
         if (index(line,'......').ne.0) return
         if (idocoo.eq.1) then
            iatoms = iatoms + 1
         else
            natoms = natoms + 1
         endif
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.gt.2) goto 1000
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3) goto 1000
         if (idocoo.eq.1) then
            ianz(iatoms) = rtype
         else
            nat(natoms) = rtype
         endif
         do i=1,3
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.3) goto 1000
            if (idobohr.eq.1) rtype = rtype / toang
            if (idocoo.eq.1) then
               coo(i,iatoms) = rtype
            else
               xyz(i,natoms) = rtype
            endif
         end do
      end do

1000  istats = 0
      call inferr('ERROR reading molecular geometry!',1)
      return
      end

      subroutine gamupd(ipnt,istat,coo,ianz)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      parameter (maxfat=1000)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat) 
      common /athlp/ iatoms, mxnat
      common /gammus/iold

      integer getlin
      logical getzmu,usnew
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character*137 line,str
      common /curlin/ line
      dimension coo(3,*),ianz(*)

      i0 = 0
      i1 = 1
      i2 = 2
c
c     Read Next Point in geometry optimization
c
      toang=0.52917706d0
      iatoms = 0
      jatoms = 0
      usnew = .false.

      if (iold.eq.1) then
         if (ipnt.eq.1) then
            call searchd(line,
     &       'ATOM      ATOMIC                      COORDINATES',
     &                'GRADIENT (HARTREE/BOHR)',istat)
         else
            call searchd(line,'COORDINATES OF ALL ATOMS ARE (ANGS)',
     &                'GRADIENT (HARTREE/BOHR)',istat)
         endif
         
         if (index(line,'ALL ATOMS').ne.0.or.
     &       index(line,'ATOM      ATOMIC').ne.0) then
            backspace iun2
            if (ipnt.eq.1) then
                call rdmolu(i2,i1,i1,i0,istat)
            else
                call rdmolu(i1,i2,i1,i1,istat)
            endif
            call search(line,'GRADIENT (HARTREE/BOHR)',istat)
         endif
         if (istat.eq.0) return
         call haszm(.false.)
      else
         if (ipnt.eq.1) then
           call searcht(line,
     &       'ATOM      ATOMIC                      COORDINATES',
     &             'GRADIENT (HARTREE/BOHR)',
     &             'CURRENT FULLY SUBSTITUTED Z-MATRIX',istat)
         else
           call searcht(line,'COORDINATES OF ALL ATOMS ARE (ANGS)',
     &             'GRADIENT (HARTREE/BOHR)',
     &             'CURRENT FULLY SUBSTITUTED Z-MATRIX',istat)
         endif
         if (istat.eq.0) return

         if (index(line,'ALL ATOMS').ne.0.or.
     &       index(line,'ATOMIC').ne.0) then
             backspace iun2
             if (ipnt.eq.1) then
                call rdmolu(i2,i1,i1,i0,istat)
             else
                call rdmolu(i1,i2,i1,i1,istat)
             endif
             call searchd(line,'GRADIENT (HARTREE/BOHR)',
     &             'CURRENT FULLY SUBSTITUTED Z-MATRIX',istat)
             if (istat.eq.0) return
         endif

         if (index(line,'CURRENT FULLY').ne.0) then
             if (getzmu(0)) then
                 call haszm(.true.)
             endif
             call search(line,'GRADIENT (HARTREE/BOHR)',istat)
             if (istat.eq.0) return
         else
             call haszm(.false.)
         endif
      endif

      call readel(line,3)

      usnew = .false.
      if (line(2:4).eq.'---') usnew = .true.

      if (.not.usnew) then
         call readel(line,2)
         iatoms = 0
      endif

      do while (getlin(0).eq.1)
         if (linlen(line).le.1) return
         jatoms = jatoms + 1
         if (.not.usnew) iatoms = iatoms + 1
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.2) goto 1000
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.1) goto 1000
         ktype = nxtwrd(str,nstr,itype,rtype)

         if (ktype.ne.3) goto 1000
         if (.not.usnew) ianz(iatoms) = rtype

         if (.not.usnew) then
            do i=1,3
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.ne.3) goto 1000
               coo(i,iatoms) = rtype
            end do
         endif

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

      subroutine geogus(formax,forrms,dismax,disrms,epoints,isav)
      implicit double precision (a-h,o-z)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      character lstr*137
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)

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

      do i=1,mxpnt
         call search(lstr,'  NSERCH=',istat)
         if (istat.eq.0) goto 100
         if (index(lstr,'FAILURE TO LOCATE').ne.0) goto 100
         ngeoms =  ngeoms + 1
         isav(i) = 1
         i1 = index(lstr,'ENERGY=')
         epoints(i) = reada(lstr,i1+7,len(lstr))
         call search(lstr,'MAXIMUM GRADIENT',istat)
         if (istat.eq.0) goto 100
         i1 = index(lstr,'MAXIMUM GRADIENT =')
         formax(i) = reada(lstr,i1+18,len(lstr))
         i1 = index(lstr,'RMS GRADIENT =')
         forrms(i) = reada(lstr,i1+14,len(lstr))

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

      subroutine cnvgus
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
c     Read Scf Convergence data
c
      rewind iun2
      icvav1 = 1
      icvav2 = 1

c     SCF convergence first point

      call searchd(line,' ITER EX',' ITER    TOTAL ENERGY',istat)
      if (istat.eq.0) goto 1000

      jstrt1 = 0

      do while (getlin(0).eq.1)
         if (linlen(line).le.1) goto 100
         if (datlin(line)) then
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.2) goto 1000
            jend1 = itype
            do i=1,3
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.3) then
                  convg1(jend1) = rtype
                  goto 50
               endif
            end do
50          continue
            if (jstrt1.eq.0) jstrt1 = jend1
         endif
      end do

100   continue


      last = .false.

200   call searchd(line,' ITER EX',' ITER    TOTAL ENERGY',istat)
      if (istat.eq.0) goto 2000

      last = .true.

      jstrt2 = 0

      do while (getlin(0).eq.1)
         if (linlen(line).le.1) goto 200
         if (datlin(line)) then
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.ne.2) goto 2000
            if (itype.gt.MAXITER) goto 2000
            jend2 = itype
            do i=1,3
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.3) then
                  convg2(jend2) = rtype
                  goto 500
               endif
            end do
500         continue
            if (jstrt2.eq.0) jstrt2 = jend2
         endif
      end do

      return

1000  icvav1 = 0
      icvav2 = 0
      return

2000  if (.not.last) then
         icvav2 = 0
      endif
      return

      end

      subroutine ugetfd(istat,coo)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /athlp/ iatoms, mxnat
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      common /curlin/ line
      character*137 line,str
      integer getlin
      dimension coo(3,*)

      istat = 1

      ivibs = 0
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

c     Read in Gamess US Frequencies

10    itel = ivibs
      call search(line,'FREQUENCY:',istat)
      if (istat.eq.0) goto 100
      backspace iun2
      if (getlin(0).eq.0) goto 100
      ktype = nxtwrd(str,nstr,itype,rtype)
      if (ktype.ne.1) goto 100
      do i=1,9
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3) then
            if (ktype.eq.1.and.nstr.eq.1) then
                if (str(1:1).eq.'I') then
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.ne.3) goto 20
                endif
            endif
            if (ktype.eq.0) goto 20
         endif
         ivibs = ivibs + 1
         freq(ivibs) = rtype
      end do
20    continue

c skip symmetry line if present
      if (getlin(0).eq.0) goto 100
      if (index(line,'SYMMETRY').ne.0) then
c        write(*,*) line(1:20)
      else
         backspace iun2
      endif

c get intensities

      if (getlin(0).eq.0) goto 100
      if (index(line,'REDUCED').ne.0) then
         if (getlin(0).eq.0) goto 100
      endif
      ktype = nxtwrd(str,nstr,itype,rtype)
      if (ktype.ne.1) goto 100
      if (nstr.eq.2.and.str.eq.'IR') 
     &   ktype = nxtwrd(str,nstr,itype,rtype)
      do i=1,9
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3) then
            if (ktype.eq.1.and.nstr.eq.1) then
                if (str(1:1).eq.'I') then
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.ne.3) goto 30
                endif
            endif
            if (ktype.eq.0) goto 30
         endif
         ihasi = 1
         itel = itel + 1
         frint(itel) = rtype
         ramint(itel) = 0.0d0
      end do
30    continue
      goto 10

100   if (ivibs.eq.0) istat = 0

      nfreq = ivibs
      call parptr(1,freq,freq,nfreq)
      call parptr(112,frint,ramint,ihasi)

      return
      end 

      subroutine ucoorg(idebug,ifreq,istat)
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

c     Get ifreq 'th Norm. Cordinates from Gamess US Output

      istat = 1
      rewind iun2
      call iatnox(iatoms)

10    call search(line,'FREQUENCY:',istat)
      if (istat.eq.0) goto 100
      ivibs = 0
      ihfor = 0
      backspace iun2
      backspace iun2
      if (getlin(0).eq.0) goto 1000
      if (index(line,'FORCE CONST.:').ne.0) then
         ihfor = 1
         backspace iun2
         backspace iun2
         if (getlin(0).eq.0) goto 1000
      endif
      do i=1,9
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.2) goto 20
         ivibs = ivibs + 1
         if (itype.eq.ifreq) goto 100
      end do
20    continue
      if (ihfor.eq.1) then
         call readel(line,3)
      else
         call readel(line,2)
      endif
      goto 10

100   if (ihfor.eq.1) then
         call readel(line,3)
      else
         call readel(line,2)
      endif
c     Skip next line in Gamess US 20 JUN 2002 Output
      if (index(line,'SYMMETRY:').ne.0) call readel(line,1)
      if (index(line,'REDUCED MASS:').ne.0) call readel(line,1)
      if (index(line,'INTENSITY:').ne.0) call readel(line,1)

      do i=1,iatoms
         do j=1,3
            if (getlin(0).eq.0) goto 1000
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.2) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)
            endif
            if (ktype.ne.1) goto 1000
            do k=1,ivibs
               ktype = nxtwrd(str,nstr,itype,rtype)
            end do
            if (ktype.ne.3) goto 1000
            a(j,i) = rtype
         end do
      end do

      if (idebug.eq.1) call prtfr(ifreq)
      return

1000  istat = 0
      call inferr('Error reading Norm. Coords. !',0)
      return
      end

      subroutine getzzz(istatz,
     &              bl,alph,bet,ibl,ialph,ibet,imap,ianz,iz)

c this is really getzmu

      implicit double precision (a-h,o-z), integer ( i-n)
      integer getlin
      character*137 str
      character*2 catom,catomt,tolowf,iel
      common /zmfrst/ ihaszm, nz, mxzat
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      dimension iel(100)
      character*137 line
      common /curlin/ line
      dimension bl(*),alph(*),bet(*),ibl(*),ialph(*),ibet(*),
     &          imap(*),ianz(*),iz(4,*)
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

      istatz = 0
      maxnz = mxzat
      do i=1,3
         do j=1,4
           iz(j,i) = 0
          end do
      end do
      nz = 0

      do while (getlin(0).eq.1)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.1) then
           if (nstr.le.8) then
              nz = nz + 1
              imap(nz) = 0
              numv = 3
              if (nz.le.3) numv = nz - 1
c
c Atom String
c
              if (nstr.eq.1) then
                 catomt(1:1) = str(1:1)
                 catomt(2:2) = ' '
              else
                 catomt = str(1:2)
              endif
              catom = tolowf(catomt)
              do j=1,100
                 if (catom .eq. iel(j)) ianz(nz) = j - 1
              end do
              do j=1,numv
c
c Connectivity
c
                 if (nxtwrd(str,nstr,itype,rtype).ne.2) then
                    goto 100
                 else
                    iz(j,nz) = itype
                 endif
c
c Variable
c
                 ktype = nxtwrd(str,nstr,itype,rtype)
                 if (ktype.eq.0.or.ktype.eq.1) then
                    goto 100
                 elseif (ktype.eq.2) then
                     tmpvar = 1.0d0*itype 
                 elseif (ktype.eq.3) then
                     tmpvar = rtype 
                 endif
                 if (j.eq.1) then
                    bl(nz) = tmpvar
                 elseif (j.eq.2) then
                    alph(nz) = tmpvar 
                 elseif (j.eq.3) then
                    bet(nz) = tmpvar
                 endif
              end do
c
c Check for Gamess ITYPE
c
              ktype = nxtwrd(str,nstr,itype,rtype)
              iz(4,nz) = 0
              if (ktype.ne.0) then
                 if (ktype.eq.2.and.itype.eq.1.or.itype.eq.0) then
                     iz(4,nz) = itype
                 elseif (ktype.eq.3.and.rtype.eq.-1.0d0) then
                     iz(4,nz) = -1
                 else
                     goto 100
                 endif
              endif
           endif
c
c Empty Line
c
         elseif (ktype.eq.0.and.nz.gt.1) then
           goto 200
         endif
      end do
c
c Out of Lines, Didnt Find Zmat
c
100   continue
      
      return

200   continue
c
c Found Zmat
c

      istatz = 1
      do i=1,nz
         ibl(i) = 1
         ialph(i) = 1
         ibet(i) = 1
      end do
c      do i=1,nz
c         print*,ianz(i),bl(i),alph(i),bet(i),(iz(j,i),j=1,4)
c      end do

      return
      end

      subroutine gamunmr
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      character*137 line,str
      common /curlin/ line
      integer getlin
      common /nmr/    shlnuc(numatm),ihsnmr
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      call search(line,'GIAO CHEMICAL SHIELDING TENSOR (PPM):',istat)

      if (istat.ne.0) then
          call readel(line,1)
          do i=1,natoms
             call readel(line,6)
             if (getlin(0).ne.1) goto 1000
             print*,'i=',i,line
             ktype = nxtwrd(str,nstr,itype,rtype)
             print*,'ktype=',ktype
             if (ktype.eq.3) then
                shlnuc(i) = rtype
             else
                goto 1000
             endif
          end do
          ihsnmr = 2

      endif

      return

1000  call inferr('Error reading Isotropical Shielding !!',1)
      return
      end

