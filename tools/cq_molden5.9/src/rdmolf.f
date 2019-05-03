      subroutine rdmodd(idebug,istatio,irtype,iesp,istats,
     &                 vectrs,vectrb,focc,focb,eiga,eigb,
     &                 nocc,nocb,ncols,ncolb,
     &                 stoalfa,stobnorm,istos,naorbs,
     &                 formax,forrms,dismax,disrms,epoints,isav)

c THIS IS REALLY rdmolf
 
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (MAXITER=1000)
c
c     *******************************************************************
c     F.Mariotti: ADF common for store AO data
c     *******************************************************************
c
c
      common /athlp/ iatoms, mxnat

c     molecula data
      common /moldat/ natoms,norbs,nelecs,nat(numatm) 

c     atoms coordinates
      common /coord / xyz(3,numatm)

c     orbitals data, via arguments
 
      common /orbhlp/ mxorb,iuhf,ispd

c     scf convergence
      common /convrg/ convg1(MAXITER),convg2(MAXITER),cnvmax,cnvmin,
     &jstrt1,jend1,jstrt2,jend2,icvav1,icvav2

c     geo convergence
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt

      common /pseudo / ipseud,ivale(numatm)

c     reading and writing routines
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      character*137 line,str,tstr
      common /curlin/ line
      integer getlin,genaos
      logical doconv
      real eiga,eigb,stoalfa,stobnorm
      dimension vectrs(*),vectrb(*),focc(*),focb(*),eiga(*),eigb(*)
      dimension istos(5,*),stoalfa(*),stobnorm(*)
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)

c
c     ***************
c     local variables
c     ***************
c

c     Title is read but not saved in molden

      character*137 MOLTitle

c     Cell data are read but not saved in molden

      double precision dCellAng(3), dCellVec(3)

      istatio = 1
      doconv = .true.
      toang = 0.52917706d0
      isgau = 0
      istats = 0
      iatfnd = 1

c
      nelecs = 0
      relecs = 0.0d0
      norbs = 0

      if (idebug.eq.1) call inferr('subroutine rdmolf',0)
c
      irtype = 0
c
c     *********************
c     search for title line
c     *********************
c
      call rewmf
      call srchmf(line,'[Title]',istat)
c
c
      if (istat.eq.1) then
          MOLTitle = line
      else
          if (idebug.eq.1) call inferr('NO Title card found',0)
      endif

      call rewmf

c
c     ********************
c     search for Cell data
c     ********************
c

c
c     F.Mariotti COMMENT
c     Now cell data are not used. Only stored in un-common variables
c

      call srchmf(line,'[Cell]',istat)
c
      if (istat.eq.0) then
          if (idebug.eq.1) call inferr('No Cell card found',0)
          goto 5
      endif

      if (getlin(0).eq.1) then
         do i=1,3
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.2) then
               dCellAng(i) = dble(itype)
            elseif (ktype.eq.3) then
               dCellAng(i) = rtype
            else
               goto 5
            endif
         end do
         do i=1,3
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.2) then
               dCellVec(i) = dble(itype)
            elseif (ktype.eq.3) then
               dCellVec(i) = rtype
            else
               goto 5
            endif
         end do
      else
         goto 5
      endif

      if (idebug.eq.1) call inferr('MOLF: Cell data read',0)
c
5     call rewmf

c
c     *********************
c     search for Atoms data
c     *********************
c
c     Default in Angstroms
c     ADF/STO GTO works with xyz in au's
c

      call srchmf(line,'[Atoms]',istat)
      if (istat.eq.1) then
         if (idebug.eq.1) call inferr('ATOMS card found',0)
         if (icdex(line,'Angs').ne.0) doconv = .true.
         if (icdex(line,'AU').ne.0) doconv = .false.
      else
         call inferr('No ATOMS card found',1)
         iatfnd = 0
         goto 100
      end if
c
      natoms = 0
c
c     get coordinates line
c

      do while (getlin(0).eq.1)

c     Atom data

          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.1) then
             if (index(str,'[').ne.0) goto 10
             tstr = line(1:8)
          else
             print*,line
             call inferr('string expected',1)
             goto 100
          endif

          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.2) then
             na = itype
             if (na.gt.numatm) then
                call inferr('Exceeded max. number of atoms',1)
                goto 100
             endif
          else
             print*,line
             call inferr('integer expected',1)
             goto 100
          endif

c     Atoms type

          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.2) then
             nat(na) = itype
          else
             print*,line
             call inferr('integer expected',1)
             goto 100
          endif

          do i=1,3
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.3) then
                xyz(i,na) = rtype
             else
                print*,line
                call inferr('real expected',1)
                goto 100
             endif
          end do
      end do
c
10    if (idebug.eq.1) call inferr('MOLF: Atoms data read',0)

c
      natoms = na
      iatoms = na

      if (doconv) then
         do i=1,natoms
             do j=1,3
                xyz(j,i) = xyz(j,i) / toang  
             end do
         end do
      endif

      call xyzcoo(1,0,0)

c search for Pseudopotential records

      call rewmf
      call srchmf(line,'[PSEUDO]',istat)
      if (istat.eq.1) then
         if (idebug.eq.1) call inferr('PSEUDO card found',0)
         ipseud = 1
         do while (getlin(0).eq.1)

            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.1) then
               if (index(str,'[').ne.0) goto 12
            else
               print*,line
               call inferr('string expected',1)
               goto 100
           endif

           ktype = nxtwrd(str,nstr,itype,rtype)
           if (ktype.eq.2) then
              na = itype
              if (na.gt.numatm) then
                 call inferr('Exceeded max. number of atoms',1)
                 goto 100
              endif
           else
              print*,line
              call inferr('integer expected',1)
              goto 100
           endif

          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.2) then
             ivale(na) = itype
          else
             print*,line
             call inferr('integer expected',1)
             goto 100
          endif

        end do

      end if

12    continue
c
c     *******************************
c     search for Atomic Orbitals data
c     *******************************
c
      call rewmf
      ibasst = 1

      call srchmf(line,'[STO]',istat)
      if (istat.eq.0) then
         if (idebug.eq.1) call inferr('No STO card found',0)

         call rdbas(idebug,1,istat)
         if (istat.eq.0) then
            print*,line
            call inferr('Error reading Basis set',1)
            goto 100
         else
            call norml
            isgau = 1
            if (idebug.eq.1) call basprt(iun3,.false.,.false.)
            call rewmf
            call srchmf(line,'[5D',istat)
            if (istat.eq.1) then
               if (icdex(line,'[5D]').ne.0.or.
     &             icdex(line,'[5D7F]').ne.0) then
                   ido5d = 1
                   ido7f = 1
               elseif (icdex(line,'[5D10F]').ne.0) then
                   ido5d = 1
               endif
            endif
            call rewmf
            call srchmf(line,'[7F]',istat)
            if (istat.eq.1) ido7f = 1
            call rewmf
            call srchmf(line,'[9G]',istat)
            if (istat.eq.1) ido9g = 1
            nbas = genaos(.true.)
            goto 20
         endif
      else
         isgau = 0
         if (idebug.eq.1) call inferr('STO card found',0)
      end if
c
      nbas = 0
      do while (getlin(0).eq.1)

c     Atom data

          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.1) then
             if (index(str,'[').ne.0) goto 20
          else
             nbas = nbas + 1
             if (nbas.gt.mxorb) then
                ibasst = 0
                goto 15
             endif

             if (ktype.eq.2) then
                istos(1,nbas) = itype
             else
                print*,line
                call inferr('real expected',1)
                goto 100
             endif
             do i=2,5
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.2) then
                   istos(i,nbas) = itype
                else
                   print*,line
                   call inferr('real expected',1)
                   goto 100
                endif
             end do

             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.3) then
                stoalfa(nbas) = rtype
             else
                goto 100
             endif

             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.3) then
                stobnorm(nbas) = rtype
             else
                print*,line
                print*,'expected a real number'
                goto 100
             endif

          endif
15        continue
      end do

20    continue
c
      if (idebug.eq.1) call inferr('MOLF: AO data read',0)
c
      naorbs = nbas

      if (ibasst.eq.0) then
          call inferr('AO basis larger than parameter',1)
          goto 100
      endif


c
c     **********************************
c     search for Molecular Orbitals data
c     **********************************
c
      call rewmf

      call srchmf(line,'[MO]',istat)
      if (istat.eq.0) then
         call inferr('No MO card found',0)
         goto 100
      else
         if (idebug.eq.1) call inferr('MO card found',0)
      endif
c
      nmos = 0
      iuhf = 0
      ncols = 0
      ncolb = 0

      do i=1,mxorb
         focc(i) = 0.0d0
         focb(i) = 0.0d0
      end do
      do i=1,mxorb
        ii = (i-1)*mxorb
        do j=1,mxorb
           vectrs(ii+j) = 0.0d0
           vectrb(ii+j) = 0.0d0
        end do
      end do

      do while (getlin(1).eq.1)

         ktype = nxtwrd(str,nstr,itype,rtype)

         if (ktype.eq.1) then

            if (icdex(str,'Sym').ne.0) then
               iab = 0
               ktype = nxtwrd(str,nstr,itype,rtype)
            endif

            if (icdex(str,'Ene').ne.0) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.3) ener = rtype
            endif

            if (icdex(str,'Spin').ne.0) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.1) then
                  if (nstr.eq.5) then
                      iab = 0
                      ncols = ncols + 1
                      if (ncols.gt.mxorb) then
                         call inferr('N of MO larger than parameter',1)
                         goto 100
                      else
                          eiga(ncols) = real(ener)
                      endif
                  elseif (nstr.eq.4) then
                      iab = 1
                      iuhf = 1
                      ncolb = ncolb + 1
                      if (ncolb.gt.mxorb) then
                         call inferr('N of MO larger than parameter',1)
                         goto 100
                      else
                          eigb(ncolb) = real(ener)
                      endif
                  endif
               endif
            endif

            if (icdex(str,'Occup').ne.0) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.3) then
                  if (iab.eq.0) then
                     if (ncols.le.mxorb) focc(ncols) = rtype
                     relecs = relecs + rtype
                  else
                     if (ncolb.le.mxorb) focb(ncolb) = rtype
                     relecs = relecs + rtype
                  endif
               elseif (ktype.eq.2) then
                  if (iab.eq.0) then
                     if (ncols.le.mxorb) focc(ncols) = dble(itype)
                     relecs = relecs + dble(rtype)
                  else
                     if (ncolb.le.mxorb) focb(ncolb) = dble(itype)
                     relecs = relecs + dble(rtype)
                  endif
               endif
            endif

            if (index(str,'[').ne.0) goto 30

         elseif (ktype.eq.2) then
            ii = itype
            if (ii.gt.mxorb) then
               call inferr('N of MO larger than parameter',0)
               goto 100
            endif

            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.3) then
               if (iab.eq.0) then
                  vectrs((ncols-1)*mxorb+ii) = rtype
               else
                  vectrb((ncolb-1)*mxorb+ii) = rtype
               endif
            else
               print*,line
               call inferr('real expected',1)
            endif
         else
            if (ktype.eq.0) goto 30
            print*,line
            call inferr('integer expected',1)
            goto 100
         endif

      end do

      nelecs = dnint(relecs)

30    norbs = nbas
      ihasd = 1
      istats = 1

100   continue
 
      if (idebug.eq.1) call inferr('MOLF: MO data read',0)
      
      call rewmf

      call srchmf(line,'[SCFCONV]',istat)

      if (istat.eq.0) then

         call inferr('No SCFCONV card found',0)

      else
         if (idebug.eq.1) call inferr('SCFCONV card found',0)

         i1or2 = 0
         do while (getlin(0).eq.1)

            ktype = nxtwrd(str,nstr,itype,rtype)

            if (ktype.eq.1) then
                
               if (icdex(str,'scf-first').ne.0) then

                  i1or2 = 1
                  icvav1 = 1
                  ktype = nxtwrd(str,nstr,itype,rtype)
                  if (ktype.eq.2) jstrt1 = itype
                  ktype = nxtwrd(str,nstr,itype,rtype)
                  ktype = nxtwrd(str,nstr,itype,rtype)
                  if (ktype.eq.2) jend1 = itype
                  icnt = jstrt1

               elseif (icdex(str,'scf-last').ne.0) then

                  i1or2 = 2
                  icvav2 = 1
                  ktype = nxtwrd(str,nstr,itype,rtype)
                  if (ktype.eq.2) jstrt2 = itype
                  ktype = nxtwrd(str,nstr,itype,rtype)
                  ktype = nxtwrd(str,nstr,itype,rtype)
                  if (ktype.eq.2) jend2 = itype
                  icnt = jstrt2

               elseif (index(str,'[').ne.0) then
                  goto 40
               endif

            elseif (ktype.eq.3) then

               if (i1or2.eq.1) then
                   convg1(icnt) = rtype
               elseif (i1or2.eq.2) then
                   convg2(icnt) = rtype
               endif
               icnt = icnt + 1

            endif

         end do
      endif

40    continue
      call mmcnv

      call parptr(4,convg1,fdum,idum)
 
      call rewmf

      call srchmf(line,'[GEOCONV]',istat)

      if (istat.eq.0) then

         call inferr('No GEOCONV card found',0)

      else
         if (idebug.eq.1) call inferr('GEOCONV card found',0)

         nepnts = 0
         do while (getlin(0).eq.1)
            ktype = nxtwrd(str,nstr,itype,rtype)

            if (ktype.eq.1) then

               if (index(str,'[').ne.0) goto 50

               if (icdex(str,'energy').ne.0) then
                  iopt = 1
                  ieav = 1
               endif
               if (icdex(str,'max-f').ne.0) then
                  iopt = 2
                  ifmxav = 1
                  ngeoms = 0
               endif
               if (icdex(str,'rms-f').ne.0) then
                  iopt = 3
                  ifrmav = 1
                  ngeoms = 0
               endif
               if (icdex(str,'max-s').ne.0) then
                  iopt = 4
                  idmxav = 1
                  ngeoms = 0
               endif
               if (icdex(str,'rms-s').ne.0) then
                  iopt = 5
                  idrmav = 1
                  ngeoms = 0
               endif
 
               if (str(1:1).eq.'-') then
                  ngeoms = ngeoms + 1
                  isav(ngeoms) = 0
               endif

            elseif (ktype.eq.3) then

               if (iopt.eq.1) nepnts = nepnts + 1
               if (iopt.gt.1) then
                  ngeoms = ngeoms + 1
                  isav(ngeoms) = 1
               endif

               if (iopt.eq.1) then
                  epoints(nepnts) = rtype
               elseif (iopt.eq.2) then
                  formax(ngeoms) = rtype
               elseif (iopt.eq.3) then
                  forrms(ngeoms) = rtype
               elseif (iopt.eq.4) then
                  dismax(ngeoms) = rtype
               elseif (iopt.eq.5) then
                  disrms(ngeoms) = rtype
               endif

            endif

         end do
      endif

50    continue


      call rewmf

      call srchmf(line,'[GEOMETRIES]',istat)

      if (istat.eq.0) then

         call inferr('No GEOMETRIES card found',0)

      else
         if (idebug.eq.1) call inferr('GEOMETRIES card found',0)
         ihasg = 1
         if (icdex(line,'ZMAT').ne.0) ihasg = 2
         igcvav = 1
      endif

      call rewmf

      call srchmf(line,'[FORCES]',istat)
      if (istat.eq.0) then
         ifrav = 0 
         call inferr('No FORCES card found',0)
      else
         ifrav = 1 
         if (idebug.eq.1) call inferr('FORCES card found',0)
      endif

      call rewmf

      call srchmf(line,'[FREQ]',istat)

      if (istat.eq.0) then

         call inferr('No FREQ card found',0)

      else
         irtype = 4
         if (idebug.eq.1) call inferr('FREQ card found',0)
         if (iatfnd.ne.1) then
             call getfra(istat)
             call getint(istat)
         endif
      endif

      call rdsrf(iun2,istat,iesp,0,idebug)

      return
      end

      subroutine prtmold(iun,ihaszm,ipoints,
     &                 vectrs,vectrb,focc,focb,eiga,eigb,
     &                 nocc,nocb,ncols,ncolb,
     &                 stoalfa,stobnorm,istos,naorbs,
     &                 formax,forrms,dismax,disrms,epoints,isav,ianz)
c
c THIS IS REALLY prtmolf
c
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      parameter (MAXITER=1000)
      parameter (MAXPNT=2000)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)

c [GTO]
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /moldat/ natoms,norbs,nelecs,nat(numatm) 
c [ATOMS]
      common /coord / xyz(3,numatm)
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
c [SCFCONV]
      common /convrg/ convg1(MAXITER),convg2(MAXITER),cnvmax,cnvmin,
     &jstrt1,jend1,jstrt2,jend2,icvav1,icvav2
c [GEOCONV]
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      logical dozme
      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme
c [FREQ]
      common /athlp/ iatoms, mxnat
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /pseudo / ipseud,ivale(numatm)
c
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
c [STO],[MO]
      common /orbhlp/ mxorb,iuhf,ispd
      real eiga,eigb,stoalfa,stobnorm
      dimension vectrs(*),vectrb(*),focc(*),focb(*),eiga(*),eigb(*)
      dimension istos(5,*),stoalfa(*),stobnorm(*)
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)
      dimension ifrcs(MAXPNT),ianz(*)


      iunt = iun3
      iun3 = iun
      toang = 0.52917706d0

      write(iun3,'(a)') '[Molden Format]'

      if (ihasd.eq.0) goto 50

c [ATOMS]
      write(iun3,'(a)') '[Atoms] Angs'

      fact = toang
      do i=1,natoms
          write(iun3,'(a2,1x,i5,1x,i2,3(1x,f12.6))')
     &       elemnt(nat(i)),i,nat(i),(xyz(j,i)*fact,j=1,3)
      end do
      
      
      if (isgau.eq.1) then
         if (ido5d.eq.1.and.ido7f.eq.1) then
             write(iun3,'(a)') '[5D7F]'
         elseif (ido5d.eq.1.and.ido7f.eq.0) then
             write(iun3,'(a)') '[5D10F]'
         elseif (ido5d.eq.0.and.ido7f.eq.1) then
             write(iun3,'(a)') '[7F]'
         elseif (ido9g.eq.1) then
             write(iun3,'(a)') '[9G]'
         endif
c [PSEUDO]
         if (ipseud.eq.1) then
             write(iun3,'(a)') '[PSEUDO]'
             do i=1,natoms
                write(iun3,'(a2,1x,i5,1x,i3)')
     &                elemnt(nat(i)),i,ivale(i)
             end do
         endif
c [GTO]
         write(iun3,'(a)') '[GTO]'
         call basprt(iun3,.true.,.false.)
      else
c [STO]
         write(iun3,'(a)') '[STO]'
         do i=1,norbs
             write(iun3,'(5(i5,1x),2(f10.4,1x))')
     &         (istos(j,i),j=1,5),stoalfa(i),stobnorm(i)
         end do
      endif

100   format(i4,1x,f10.6)

c [MO]
      write(iun3,'(a)') '[MO]'
      do i=1,ncols
         write(iun3,'(a,f10.4)') ' Ene= ',eiga(i)
         write(iun3,'(a)') ' Spin= Alpha'
         write(iun3,'(a,f10.6)') ' Occup= ',focc(i)
         do j=1,norbs
            write(iun3,100) j,vectrs((i-1)*mxorb+j)
         end do
      end do

      if (iuhf.eq.1) then
         do i=1,ncolb
            write(iun3,'(a,f10.4)') ' Ene= ',eigb(i)
            write(iun3,'(a)') ' Spin= Beta'
            write(iun3,'(a,f10.6)') ' Occup= ',focb(i)
            do j=1,norbs
               write(iun3,100) j,vectrb((i-1)*mxorb+j)
            end do
         end do
      endif

50    continue

c [SCFCONV]
      if (icvav1.eq.1.or.icvav2.eq.1) then
         write(iun3,'(a)') '[SCFCONV]'
         if (icvav1.eq.1) then
            write(iun3,'(a,i4,a,i4)') 
     &          'scf-first ',jstrt1,' THROUGH ',jend1
            do i=jstrt1,jend1
                write(iun3,'(f14.6)') convg1(i)
            end do
         endif
         if (icvav2.eq.1) then
            write(iun3,'(a,i4,a,i4)') 
     &          'scf-last ',jstrt2,' THROUGH ',jend2
            do i=jstrt2,jend2
                write(iun3,'(f14.6)') convg2(i)
            end do
         endif
      endif

c [GEOCONV]
      i = ieav + ifmxav + ifrmav + idmxav + idrmav
      if (i.ne.0) then
         write(iun3,'(a)') '[GEOCONV]'

         if (ieav.eq.1) then
            write(iun3,'(a)') 'energy'
            do i=1,nepnts
               write(iun,'(f14.6)') epoints(i)
            end do
         endif

         if (ifmxav.eq.1) 
     &      call prtarr(iun3,'max-force',formax,ngeoms,isav)

         if (ifrmav.eq.1) 
     &      call prtarr(iun3,'rms-force',forrms,ngeoms,isav)

         if (idmxav.eq.1) 
     &      call prtarr(iun3,'max-step',dismax,ngeoms,isav)

         if (idrmav.eq.1) 
     &      call prtarr(iun3,'rms-step',disrms,ngeoms,isav)

      endif

      call rewmf

      if (ipoints.gt.1.or.ihasd.eq.0) then

         if (ihaszm.eq.1) then
             write(iun3,'(a)') '[GEOMETRIES] ZMAT'
         else
             write(iun3,'(a)') '[GEOMETRIES] XYZ' 
         endif

         if (ipoints.le.1.and.ihasd.eq.0) then
            if (ihaszm.eq.1) then
               call wrzmat(iun3,1)
            else
               call wrcart(iun3,0,0,0)
            endif
         else
            do i=1,ipoints
                call getpoi(i,0,0,1,idum,0)
                ifrcs(i) = 0
                if (ifav.eq.1) ifrcs(i) = 1
                if (ihaszm.eq.1) then
                   call wrzmat(iun3,1)
                else
                   call wrcart(iun3,0,0,0)
                endif
            end do
         endif
      endif

c [FORCES]
      if (ipoints.gt.1) then
         ifpts = 0
         do i=1,ipoints
            if (ifrcs(i).eq.1) ifpts = ifpts + 1
         end do
         if (ifpts.gt.0) then
            call rewmf
            if (iftyp.eq.3) then
               write(iun3,'(a)') '[FORCES] COORD'
            else
               write(iun3,'(a)') '[FORCES]'
            endif
            do i=1,ipoints
               call getpoi(i,0,0,1,idum,0)
               if (ifav.eq.1) then
                  write(iun3,'(a,i4)') 'point ',i
                  call wrfc(iun3)
               endif
            end do
         endif
      endif

c [FREQ]
      if (irtype.eq.4) then

         write(iun3,'(a)') '[FREQ]'
         do i=1,nfreq
            write(iun3,'(f10.4)') freq(i)
         end do

         if (ihasi.ne.0) then
            write(iun3,'(a)') '[INT]'
            do i=1,nfreq
               if (iabs(ihasi).eq.1) then
                  write(iun3,'(f10.4)') frint(i)
               elseif (iabs(ihasi).eq.2) then
                  write(iun3,'(f10.4,1x,f10.4)') frint(i),ramint(i)
               endif
            end do
         endif

         write(iun3,'(a)') '[FR-COORD]'
         do i=1,iatoms
             write(iun3,'(a2,1x,3(f12.6,1x))') 
     &        elemnt(ianz(i)),(fcoo(j,i),j=1,3)
         end do

         write(iun3,'(a)') '[FR-NORM-COORD]'
         do i=1,nfreq

            write(iun3,'(a,i5)') 'vibration ',i
            if (iftyp.eq.1) then
                call mcoord(0,i,istat)
            elseif (iftyp.eq.2.or.iftyp.eq.3) then
                if (iftyp.eq.3) then
                   call ucoorg(0,i,istat)
                else
                   call ncoorg(0,i,istat)
                endif
            elseif (iftyp.eq.4) then
                call ncoord(0,i,istat)
            endif

            do k=1,iatoms
                write(iun3,'(3(f12.6,1x))') (a(j,k),j=1,3)
            end do

         end do

      endif

      iun3 = iunt

      return
      end

      subroutine prtarr(iun,str,arr,narr,isav)
      implicit double precision (a-h,o-z)
      character*(*) str
      dimension arr(*),isav(*)

      write(iun,'(a)') str
      do i=1,narr
         if (isav(i).eq.1) then
            write(iun,'(f14.6)') arr(i)
         else
            write(iun,'(a)') '-'
         endif
      end do

      return
      end

      subroutine parsmd(molin)
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

c      call srchmf(line,'[MULTIMOL]',istat)
c      if (istat.ge.1) then


1        call nxtlin(line,jstat)
         if (jstat.eq.1) goto 100
         if (jstat.eq.2) goto 200
         caps = line
         len1  = len(line)
         call tocap(caps,len1)
         ilin = ilin + 1
         if (index(caps,'[TITL').ne.0) then
            call nxtlin(line,jstat)
            if (jstat.eq.1) goto 100
            if (jstat.eq.2) goto 200
            title = line(1:40)
            ilin = ilin + 1
         endif
         if (index(caps,'[MOLE').eq.0) goto 1

         if (nmols.lt.maxmol) then
            nmols = nmols + 1
            molin(nmols) = ilin
            if (nmols.gt.1) call parsfn(title,linlen(title),13)
            title = 'molecule '//gstr(nmols)
         endif

         goto 1

c      endif

100   ielin = ilin
      if (nmols.eq.0) then
         nmols = 1
         title = 'molecule '//gstr(nmols)
      endif
      call parsfn(title,linlen(title),13)

      return

200   return
      end

      subroutine srchmf(line,name,istat)
      implicit double precision (a-h,o-z)
      common /mfdata/ nmols,imol,iendmf,ielin,maxmol,mollin
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      common /mflin/ linmf
      character*(*) line,name
      character*80 caps
      character*137 caps2

      istat = 0
      caps  = name  
      len1  = len(name)
      call tocap(caps,len1)
1     call nxtlin(line,jstat)
      if (jstat.eq.1) goto 100
      if (jstat.eq.2) goto 200

      linmf = linmf + 1
      if (linmf.gt.iendmf) goto 100
      caps2 = line
      len2 = len(line)
      call tocap(caps2,len2)
      if (index(caps2(1:len2),caps(1:len1)).eq.0) goto 1
 
      istat = 1
      return

100   call rewmf
      return

200   continue
      print*,'error search mf'

      return
      end 

      subroutine srcdmf(line,name1,name2,istat)
      implicit double precision (a-h,o-z)
      common /mfdata/ nmols,imol,iendmf,ielin,maxmol,mollin
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      common /mflin/ linmf
      character*(*) line,name1,name2
      character*80 caps
      character*137 caps1,caps2

      istat = 0
      caps1  = name1
      caps2  = name2
      len1  = len(name1)
      len2  = len(name2)
      call tocap(caps1,len1)
      call tocap(caps2,len2)
1     call nxtlin(line,jstat)
      if (jstat.eq.1) goto 100
      if (jstat.eq.2) goto 200

      linmf = linmf + 1
      if (linmf.gt.iendmf) goto 100
      if (index(line,name1).eq.0.and.index(line,name2).eq.0
     &.and.index(line,caps1(1:len1)).eq.0
     &.and.index(line,caps2(1:len2)).eq.0) go to 1

 
      istat = 1
      return

100   call rewmf
      return

200   continue
      print*,'error searchd mf'

      return
      end 

      subroutine srctmf(line,name1,name2,name3,istat)
      implicit double precision (a-h,o-z)
      common /mfdata/ nmols,imol,iendmf,ielin,maxmol,mollin
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      common /mflin/ linmf
      character*(*) line,name1,name2,name3
      character*80 caps
      character*137 caps1,caps2,caps3

      istat = 0
      caps1  = name1
      caps2  = name2
      caps3  = name3
      len1  = len(name1)
      len2  = len(name2)
      len3  = len(name3)
      call tocap(caps1,len1)
      call tocap(caps2,len2)
      call tocap(caps3,len3)
1     call nxtlin(line,jstat)
      if (jstat.eq.1) goto 100
      if (jstat.eq.2) goto 200

      linmf = linmf + 1
      if (linmf.gt.iendmf) goto 100
      if (index(line,name1).eq.0.and.index(line,name2).eq.0
     &    .and.index(line,name3).eq.0
     &.and.index(line,caps1(1:len1)).eq.0
     &.and.index(line,caps2(1:len2)).eq.0
     &.and.index(line,caps3(1:len3)).eq.0) go to 1

 
      istat = 1
      return

100   call rewmf
      return

200   continue
      print*,'error searchd mf'

      return
      end 

      subroutine rdmf(line,cline,istat)
      implicit double precision (a-h,o-z)
      common /mfdata/ nmols,imol,iendmf,ielin,maxmol,mollin
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      common /mflin/ linmf
      character*(*) line,cline

      istat = 0
      call nxtlin(line,jstat)
      if (jstat.eq.1) goto 100
      if (jstat.eq.2) goto 200

      linmf = linmf + 1
      if (linmf.gt.iendmf) goto 100
      cline = line
      len2 = len(line)
      call tocap(cline,len2)

      istat = 1
      return

100   call rewmf
      return

200   continue
      print*,'error search mf'
      return
c
c
      end 

      subroutine rewmd(molin)
      implicit double precision (a-h,o-z)
      common /mfdata/ nmols,imol,iendmf,ielin,maxmol,mollin
      common /mflin/ linmf
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*137 line
      dimension molin(*)

      call rewfil
      linmf = 0

      imol2 = imol
      if (imol.eq.0) imol2=1
      
      do i=1,molin(imol2)
          call nxtlin(line,jstat)
          if (jstat.eq.1) goto 100
          if (jstat.eq.2) goto 200

          linmf = linmf + 1
      end do
      call bckfil
      linmf = linmf - 1

      if (imol.ge.nmols) then
         iendmf = ielin
      else
         iendmf = molin(imol+1)
      endif

      return

100   continue
200   continue
      print*,'error rewind mf'
      return

      end 

      subroutine setmf(im)
      implicit double precision (a-h,o-z)
      common /mfdata/ nmols,imol,iendmf,ielin,maxmol,mollin

      imol = im
      if (im.gt.nmols) imol = nmols

      return
      end 

      integer function getmf()
      integer linmf
      common /mflin/ linmf

      getmf = linmf

      return
      end 

      subroutine putmf(iline)
      implicit double precision (a-h,o-z)
      common /mflin/ linmf
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*137 line

      call rewfil

      do i=1,iline
          call nxtlin(line,jstat)
          if (jstat.eq.1) goto 100
          if (jstat.eq.2) goto 100
      end do
      call bckfil

      linmf = iline

100   return
      end 

      subroutine wrfd(iun3,coo)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      parameter (maxfat=1000)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat) 
      common /athlp/ iatoms, mxnat
      logical ctoz,molpot,elpot,chpot
      common /choic/ iftyp,isbin,ctoz,molpot,elpot,chpot
      dimension coo(3,*)


      fsgn = 1.0d0
      if (iftyp.eq.2.or.iftyp.eq.3) fsgn = -1.0d0

      write(iun3,'(i4)') iatoms
      do i=1,iatoms
         if (iftyp.eq.3) then
            write(iun3,'(6(1x,f12.6))') (fsgn*fxyz(j,i),j=1,3),
     &            (coo(j,i),j=1,3)
         else
            write(iun3,'(3(1x,f12.6))') (fsgn*fxyz(j,i),j=1,3)
         endif
      end do

      return
      end

      subroutine rdfd(ipnt,istats,coo)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      parameter (mxel=100)
      parameter (maxfat=1000)
      logical gnreal,addcoo
      integer getlin
      character*137 str
      character*137 line
      common /curlin/ line
      common /athlp/ iatoms, mxnat
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat) 
      dimension ftmp(3),coo(3,*)

      istats = 1
      addcoo = .false.

c      print*,'rdfc point ',ipnt
      call rewmf
      call srchmf(line,'[FORCES]',istat)
      if (istat.eq.1) then
         if (icdex(line,'COORD').ne.0) addcoo = .true.
      else
         goto 100
      endif

      do while(.true.)
         call srchmf(line,'point',istat)
         if (istat.eq.1) then
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.1) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.2) then
                  if (itype.eq.ipnt) goto 10
               else
                  goto 100
               endif
            else
               goto 100
            endif
         else
            goto 100
         endif
      end do

10    if (getlin(0).eq.1) then
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            iatoms = itype
         else
            goto 100
         endif
      else
         goto 100
      endif

      do i=1,iatoms
         if (getlin(0).eq.1) then
            if (.not.gnreal(ftmp,3,.false.)) goto 100
            do j=1,3
               fxyz(j,i) = ftmp(j)
            end do
c            print*,'fxyz ',i,' ',(fxyz(j,i),j=1,3)
            if (addcoo) then
               if (.not.gnreal(coo(1,i),3,.false.)) goto 100
            endif
         else
            goto 100
         endif
      end do

      return

100   istats = 0
      return
      end

      subroutine prsgmd(iopt,molin)
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

         if ((iopt.eq.3.or.iopt.eq.1).and.nmols.eq.0) then
            ilin = ilin + 1
            nmols = 1
            molin(nmols) = ilin
            call parsfn(line,linlen(line),13)
            call nxtlin(line,jstat)
            if (jstat.eq.1) goto 100
            if (jstat.eq.2) goto 200

         endif
         caps = line
         len1  = len(line)
         call tocap(caps,len1)
         ilin = ilin + 1
         if (iopt.eq.1) then
            if (index(caps,'@<TRIPOS>MOLECULE').eq.0) goto 1
         elseif (iopt.eq.2) then
            if (index(caps,'* O   R   C   A *').eq.0) goto 1
         elseif (iopt.eq.3) then
            if (index(caps,'$$$$').eq.0) goto 1
         elseif (iopt.eq.4) then
            if (index(caps,'PREPARE MODULE').ne.0) then
                if (nmols.ge.0) nmols = nmols - 1
            endif
            if (index(caps,'QM/MM INTERFACE MODULE').ne.0) then
                if (nmols.ge.0) nmols = nmols - 1
            endif
            if (index(caps,'NWCHEM INPUT MODULE').eq.0) goto 1
         endif

         if (nmols.lt.maxmol) then
            nmols = nmols + 1
            molin(nmols) = ilin
            if (iopt.eq.1) then
               call nxtlin(line,jstat)
               if (jstat.eq.1) goto 100
               if (jstat.eq.2) goto 200

               if (line(1:16).eq.'Molden generated') then
                  title = line(18:40)
               else
                  title = line(1:40)
               endif
               ilin = ilin + 1
               call parsfn(title,linlen(title),13)
            elseif (iopt.eq.2) then
               title = 'molecule '//gstr(nmols)
               if (nmols.gt.1) call parsfn(title,linlen(title),13)
            elseif (iopt.eq.3) then

               call nxtlin(line,jstat)
               if (jstat.eq.1) goto 100
               if (jstat.eq.2) goto 200

               ilin = ilin + 1
               molin(nmols) = ilin
               call parsfn(line,linlen(line),13)
            elseif (iopt.eq.4) then
               call redel(line,4)
               ilin = ilin + 4
               call lsparm(line,ll)
               if (index(line,'WARN').ne.0) ll = 0
               do while(ll.eq.0.or.ll.eq.1)
                  call redel(line,1)
                  ilin = ilin + 1
                  call lsparm(line,ll)
               end do 
               title = line(1:40)
               molin(nmols) = ilin
               call parsfn(title,linlen(title),13)
            endif
         endif

         goto 1


100   ielin = ilin
      if (nmols.eq.0) then
         nmols = 1
         title = 'molecule '//gstr(nmols)
      endif
      if (iopt.eq.3) then
         nmols = nmols + 1
         molin(nmols) = ilin
      endif
      if (iopt.eq.3.or.iopt.eq.4) then
         nmols = nmols - 1
      else
         if (iopt.ne.1) then
            call parsfn(title,linlen(title),13)
         endif
      endif

      return

200   return
      end


      subroutine rotfil
      implicit double precision (a-h,o-z)
      integer getlin
      logical gnreal,opfil
      common /filrot/ctmp(3,3),ir3(3)
      common /curlin/ line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line,str,tstr

      toang = 0.52917706d0
      iun2t = iun2
      iun2 = 46

      if (opfil(46,'rotfil',6,1,1,0)) then
         i = 0
         do while (.true.)
         if (getlin(0).eq.1) then
               i = i + 1
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.2) then
                  ir3(i) = itype
               else
                  goto 100
               endif
               if (.not.gnreal(ctmp(1,i),3,.false.)) goto 100
               do j=1,3
                  ctmp(j,i) = ctmp(j,i) / toang
               end do
         else
            goto 10
         endif
         end do
10       print*,'ir3=',(ir3(i),i=1,3)
         print*,'ctmp(1,1)=',(ctmp(i,1),i=1,3)
         print*,'ctmp(1,2)=',(ctmp(i,2),i=1,3)
         print*,'ctmp(1,3)=',(ctmp(i,3),i=1,3)
      endif

      iun2 = iun2t
      return

100   print*,'error parsing rotfil'
      iun2 = iun2t
      return
      end
