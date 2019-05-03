      subroutine rdgdd(idebug,ibefo,istatio,ioxyz,irtype,istats,
     &                 vectrs,vectrb,focc,focb,eiga,eigb,
     &                 nocc,nocb,ncols,ncolb,coo,ianz)

c THIS IS REALLY rdgam

      implicit double precision (a-h,o-z), integer (i-n)
      common /athlp/ iatoms, mxnat
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs, nat(numatm)
      common /coord / xyz(3,numatm)
      common /orbhlp/ mxorb,iuhf,ispd
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      parameter (mxel=100)
      character*2 elemnt
      common /elem/elemnt(mxel)
      character*137 line
      logical local
      real eiga,eigb
      dimension vectrs(*),vectrb(*),focc(*),focb(*),eiga(*),eigb(*)
      dimension coo(3,*),ianz(*)

      istats = 1
      ioxyz = 0
      istatio = 0
cd     write(iun3,'(a)')'enter subroutine rdgam'

      irtype = 0
      call search(line,'run type',istat)
      if (istat.eq.0) goto 1000
      call tocap(line,136)
      if (index(line,'SCF').ne.0)      irtype = 1
      if (index(line,'OPTIMIZE').ne.0) irtype = 2
      if (index(line,'OPTXYZ').ne.0)   irtype = 3
      if (index(line,'FORCE').ne.0)    irtype = 4
      if (index(line,'HESSIAN').ne.0)  irtype = 4
      rewind iun2

      if (ibefo.eq.1.or.irtype.eq.1.or.irtype.eq.4) goto 17
c
c======= read nuclear coordinates ================
c
      ncnt = 0
      rewind iun2
111   continue
         call search(line,'nuclear coordinates',istat)
         if (istat.eq.0) then
            if (ioxyz.eq.1) then
               call search(line,'molecular orbitals',istat)
               if (istat.eq.0) ncnt = 0
            endif
            if ((ioxyz.eq.1.and.ncnt.gt.1).or.
     &          (ioxyz.eq.0.and.ncnt.ge.1)) then
               istatio = 1
               goto 12
            else
               goto 17
            endif
         endif
         if (ncnt.eq.0.and.line(2:6).eq.'point') ioxyz = 1
         ncnt = ncnt + 1
            
         call readel(line,2)
         call rdgeom(numatm,xyz,nat,natoms,istat)
      goto 111

17    istatio = 0
      call search(line,'input z-matrix',istat)
      if (istat.eq.0) ioxyz = 1
      if (ioxyz.eq.0) then
         call readel(line,1)
         call search(line,'z-matrix',istat)
         if (istat.eq.0) goto 1000
         call readel(line,2)
         call convzmat(coo,ianz,iatoms,0,1,0)
      endif
      call rdmolg(istat)
      if (istat.eq.0) goto 1000

12    continue
      if (istatio.eq.1) then
         call inferr('Using Density of stationary point',0)
      else
         call inferr('Using Density of first point',0)
      endif

      call rdbasg(idebug,.false.,istat)
      if (istat.eq.0) goto 1000
c
c============ read norbs nelecs ==================
c
      do i=1,norbs
         focc(i) = 0.0d0
         eiga(i) = 0.0d0
         focb(i) = 0.0d0
         eigb(i) = 0.0d0
      end do
c      call search(line,'total number of basis func',istat)
      call search(line,'basis functions',istat)
      if (index(line,'cartesian').ne.0) then
         read(line,'(45x,i5)',err=1000) norbs
         read(iun2,'(a)',err=1000) line
         read(iun2,'(a)',err=1000) line
         read(line,'(45x,i5)',err=1000) nelecs
      else
         read(line,'(38x,i5)',err=1000) norbs
         read(iun2,'(a)',err=1000) line
         read(line,'(38x,i5)',err=1000) nelecs
      endif
      if (norbs .gt. mxorb) then
         call inferr('Exceeding MaxNum of Orbitals!',1)
         goto 1000
      endif
      call searchd(line,'scf type','scf mode',istat)
      if (index(line,'uhf').ne.0) goto 7000
c
c     read occupancies
c
      if (index(line,'casscf').ne.0.and.
     &    (irtype.eq.1.or.irtype.eq.4)) then
         call search(line,'active space',istat)
         if (istat.ne.0) then
            read(line,'(24x,i3,3x,i3)',err=1000)ncore,nact
            nact=nact-ncore
         endif
         call search(line,
     &        'active orbitals to diagonalise density matrix',istat)
         if (istat.ne.0) then
            call readel(line,6)
            call rdcasocc(ncore,nact,focc,nocc,istat)
            if (istat.eq.0) then
               ncols = norbs
               rewind iun2
            else
               ncols = min0(nocc+12,norbs)
            endif
         else
            ncols = norbs
            call inferr('Missing occupation numbers use CANON 1',1)
         endif
      else
         call rdocc(0,norbs,focc,eiga,nocc,istat)
         if (istat.eq.0) then
            ncols = norbs
            rewind iun2
         else
            ncols = min0(nocc+5,norbs)
         endif
      endif
c
c     read vectors
c
      local = .false.
      call searchq(line,'gvb natural orbital','eigenvectors',
     &             'alpha spin lmo',
     &             'casscf mos (canonicalised)',istat)

      if (istat.eq.0) then
         if (istatio.eq.1) then
            rewind iun2
         else
            goto 1000
         endif
      else
         if (index(line,'gvb natural orbital').ne.0) then
           call readel(line,11)
         else if (index(line,'alpha spin lmo').ne.0) then
           call readel(line,5)
           local = .true.
         else if (index(line,'casscf mos (canonicalised)').ne.0) then
           call readel(line,4)
           local = .true.
         else
           call readel(line,9)
         endif
         call readvv(vectrs,norbs,nocc,local)
         if (ibefo.eq.1) then
cd           write(iun3,'(a)')'leave subroutine rdgam'
            return
         endif
      endif
c
c======= try for the vectors of an optimise or saddle job====
c
6000  call searcht(line,'m.o. irrep        orbital',
     &                  'm.o.  irrep  orbital',
     &                  'm.o. irrep  orbital',istat)
      if (istat.eq.1) then
         backspace iun2
         call rdocc(0,norbs,focc,eiga,nocc,istat)
         if (istat.eq.0) goto 1000
         ncols = min0(nocc+5,norbs)
         call search(line,'molecular orbitals',istat)
         if (istat.eq.0) goto 1000
         call readel(line,9)
         call readvv(vectrs,norbs,nocc,.false.)
      endif
c
c try for spinfree mp2
c
c      call search(line,'spinfree mp2 natural orbitals',istat)
c      if (istat.eq.1) then
c      endif

cd     write(iun3,'(a)')'leave subroutine rdgam'
      return
c
c====== uhf part ==============================
c
7000  iuhf = 1
      call rdocc(1,norbs,focc,eiga,nocc,istat)
      if (istat.eq.0) goto 1000
      ncols = min0(nocc+5,norbs)
      call search(line,'eigenvectors',istat)
      if (istat.eq.0) then
          call inferr('Missing Alpha eigenvectors',1)
          goto 1000
      endif
      call readel(line,9)
      call readvv(vectrs,norbs,nocc,.false.)
c
      call rdocc(2,norbs,focb,eigb,nocb,istat)
      if (istat.eq.0) goto 1000
      ncolb = min0(nocb+5,norbs)
      call search(line,'eigenvectors',istat)
      if (istat.eq.0) then
          call inferr('Missing Beta eigenvectors',1)
          goto 1000
      endif
      call readel(line,9)
      call readvv(vectrb,norbs,nocb,.false.)
      if (ibefo.eq.1) then
cd        write(iun3,'(a)')'leave subroutine rdgam'
         return
      endif
c
c======= try for the vectors of an optimise or saddle job====
c
      if (istatio.eq.0) then
cd        write(iun3,'(a)')'leave subroutine rdgam'
         return
      endif

      call rdocc(1,norbs,focc,eiga,nocc,istat)
      if (istat.eq.0) goto 1000
      ncols = min0(nocc+5,norbs)
      call search(line,'molecular orbitals',istat)
      if (istat.eq.1) then
         call readel(line,9)
         call readvv(vectrs,norbs,nocc,.false.)
         call rdocc(2,norbs,focb,eigb,nocb,istat)
         if (istat.eq.0) goto 1000
         ncolb = min0(nocb+5,norbs)
         call search(line,'molecular orbitals',istat)
         if (istat.eq.0) then
            call inferr('missing Beta orbitals or occupancies',1)
            goto 1000
         endif
         call readel(line,9)
         call readvv(vectrb,norbs,nocb,.false.)
      endif
cd     write(iun3,'(a)')'leave subroutine rdgam'
      return

1000  istats = 0
      call inferr('ERROR reading GAMESS output file!',1)
cd     write(iun3,'(a)')'leave subroutine rdgam'
      return
      end

      subroutine rdbasg(idebug,us,istats)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      character*137 line
      character*8  dchar
      character*4  dchar1, ylabel(26)
      common /moldat/ natoms, norbs, nelecs, nat(numatm)
      common /coord / xyz(3,numatm)
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      logical onewat, us
      dimension nangm(26),ncons(26)
      data ylabel /'1s','2s','2p','2sp','3s','3p','3sp',
     +     '3d','4s','4p','4sp','4d','5s','5p','5d',
     +     's','p','d','f','g','k','l','m','5sp','****','****'/
      data nangm /0,0,1,1,0,1,1,2,0,1,1,2,0,1,2,0,1,2,3,4,0,1,
     +            2,1,0,0/
      data ncons /0,0,1,0,0,1,0,2,0,1,0,2,0,1,2,0,1,2,3,4,0,0,
     +            0,0,0,0/
c ?? what about k,l,m ??
c answ k = s, l = sp, and m = spd; we dont handle spd yet

c
c====== read basis set info ======================
c
      istats = 1
      if (idebug.eq.1) write(iun3,'(a)')'enter subroutine rdbasg'
      rewind iun2

      iguspc = 0
      if (us) then
         call search(line,'ATOMIC BASIS SET',istat)
         if (istat.eq.0) then
            call inferr('no ATOMIC BASIS SET found !',1)
            goto 1000
         endif
         if (index(line,'ATOMIC BASIS SET').eq.5) iguspc = 1
         call readel(line,5)
      else
         call search(line,'molecular basis set',istat)
         if (istat.eq.0) then
            call inferr('no molecular basis set found !',1)
            goto 1000
         endif
         call readel(line,9)
      endif

c start loop

      nshold = 0
      iprold = 1
      idf = 0
      iatcnt = 0
      igusn = 0

100   read(iun2,'(a)',err=1000) line
      if (iguspc.eq.1) line = ' '//line
      read(line,'(1x,a)',err=1000) dchar
      if (dchar.eq.'========'.or.dchar.eq.'TOTAL NU') 
     &   goto 200

      if (index(line,'.').eq.0.and.linlen(line).gt.1) then

c new atom

         onewat = .true.
         iatcnt = iatcnt + 1
      else
         if (linlen(line).gt.1) then
            if (us) then
               dchar1 = '    '
               if (line(8:10).eq.'   ') then
                  read(line,'(10x,a)') dchar1(1:2)
               else
                  read(line,'(7x,a)') dchar1(1:2)
               endif
            else
               read(line,'(21x,a)') dchar1
            endif
            call tolow(dchar1,4)

c data line

            if (index(dchar1,'sp').ne.0.or.index(dchar1,'l').ne.0) then

              if (us) then

c gamess-us
               if (index(line,'(').eq.0 ) then
                  igusn = 1
                  if (index(line,'.').eq.34) then
                     read(line,9005,err=1000) 
     &                  nshell,dchar1,iprim,ex,cd1,cd2
                  elseif (index(line,'.').eq.28) then
                     read(line,9004,err=1000) 
     &                  nshell,dchar1,iprim,ex,cd1,cd2
                  else
                     read(line,9003,err=1000) 
     &                  nshell,dchar1,iprim,ex,cd1,cd2
                  endif
               else
                  if (index(line,'(').eq.48) then
                     read(line,9002,err=1000) 
     &                  nshell,dchar1,iprim,ex,cd1,cd2
                  else
                     read(line,9001,err=1000) 
     &                  nshell,dchar1,iprim,ex,cd1,cd2
                  endif
               endif

              else

c gamess-uk

               read(line,9000,err=1000) nshell,dchar1,iprim,ex,cd1,cd2

              endif

            else

              if (us) then

c gamess-us
               if (index(line,'(').eq.0) then
                  igusn = 1
                  if (index(line,'.').eq.34) then
                     read(line,9005,err=1000) nshell,dchar1,iprim,ex,cd1
                  elseif (index(line,'.').eq.28) then
                     read(line,9004,err=1000) nshell,dchar1,iprim,ex,cd1
                  else
                     read(line,9003,err=1000) nshell,dchar1,iprim,ex,cd1
                  endif
               else
                  if (index(line,'(').eq.48) then
                     read(line,8002,err=1000) nshell,dchar1,iprim,ex,cd1
                  else
                     read(line,8001,err=1000) nshell,dchar1,iprim,ex,cd1
                  endif
               endif

              else

c gamess-uk
               read(line,8000,err=1000) nshell,dchar1,iprim,ex,cd1

              endif

            endif

            call tolow(dchar1,4)
            exx(iprim) = ex
            ityp = locatc(ylabel,26,dchar1)
      go to (1220,1220,1240,1280,1220,1240,1280,1260,1220,1240,1280,
     +      1260,1220,1240,1260,1220,1240,1260,1250,1270,
     +      1220,1280,1260,1280),ityp
c S
 1220       if (igusn.eq.1) call renorm(ex,cd1,0)
            c1(iprim) = cd1
            go to 1300
c P
 1240       if (igusn.eq.1) call renorm(ex,cd1,1)
            c2(iprim) = cd1 
            go to 1300
c D
 1260       if (igusn.eq.1) call renorm(ex,cd1,2)
            idf=idf+1
            c3(idf) = cd1 
            go to 1300
c F
 1250       if (igusn.eq.1) call renorm(ex,cd1,3)
            idf=idf+1
            c4(idf) = cd1 
            go to 1300
 1270       if (igusn.eq.1) call renorm(ex,cd1,4)
            idf=idf+1
            c5(idf) = cd1 
            go to 1300
c SP
 1280       if (igusn.eq.1) then
                call renorm(ex,cd1,0)
                call renorm(ex,cd2,1)
            endif
            c1(iprim) = cd1 
            c2(iprim) = cd2 
 1300       continue

            if (nshell.ne.nshold) then

c   starting new shell

              if (nshold.ne.0) then

                 shelln(nshold) = iprim-iprold
                 iprold = iprim

                 if (ngap.ne.0.and.onewat) then

c   copy read shells to skipped shells

                    nskip = nshold - nfirst + 1
                    natsk = ngap / nskip
                    nanci = nfirst - ngap
                    do 2000 j=0,natsk
                    do 2000 i=0,nskip-1
                       shella(nanci+j*nskip+i) = shella(nfirst+i)
                       shelln(nanci+j*nskip+i) = shelln(nfirst+i)
                       shellc(nanci+j*nskip+i) = shellc(nfirst+i)
                       shellt(nanci+j*nskip+i) = shellt(nfirst+i)
                       shladf(nanci+j*nskip+i) = shladf(nfirst+i)
                       gx(nanci+j*nskip+i) = xyz(1,iatcnt-1+j)
                       gy(nanci+j*nskip+i) = xyz(2,iatcnt-1+j)
                       gz(nanci+j*nskip+i) = xyz(3,iatcnt-1+j)
2000                continue
                    iatcnt = iatcnt + natsk
                 endif
              endif
              if (onewat) then
                  if (nshell-1.ne.nshold) then

c  have equivalent centres been skipped

                     ngap = nshell - nshold - 1
                     nfirst= nshell
                  else
                     ngap = 0
                  endif
              endif
              shella(nshell) = iprim
              shladf(nshell) = idf
              shellt(nshell) = nangm(ityp)
              shellc(nshell) = ncons(ityp)
              gx(nshell) = xyz(1,iatcnt)
              gy(nshell) = xyz(2,iatcnt)
              gz(nshell) = xyz(3,iatcnt)
            endif
            nshold = nshell
            onewat = .false.
         endif
      endif
      goto 100
200   continue

c     round up last shell

      shelln(nshell) = iprim - iprold + 1
      if (ngap.ne.0) then
c     copy read shells to skipped shells
         nskip = nshold - nfirst + 1
         natsk = ngap / nskip
         nanci = nfirst - ngap
         do 3000 j=0,natsk
         do 3000 i=0,nskip-1
            shella(nanci+j*nskip+i) = shella(nfirst+i)
            shelln(nanci+j*nskip+i) = shelln(nfirst+i)
            shellc(nanci+j*nskip+i) = shellc(nfirst+i)
            shellt(nanci+j*nskip+i) = shellt(nfirst+i)
            shladf(nanci+j*nskip+i) = shladf(nfirst+i)
            gx(nanci+j*nskip+i) = xyz(1,iatcnt+j)
            gy(nanci+j*nskip+i) = xyz(2,iatcnt+j)
            gz(nanci+j*nskip+i) = xyz(3,iatcnt+j)
3000     continue
      endif
      do 520 i=1,nshell
         do 530 j=1,natoms
            d1 = (gx(i) - xyz(1,j))**2
            d2 = (gy(i) - xyz(2,j))**2
            d3 = (gz(i) - xyz(3,j))**2
            dtot=d1+d2+d3
            if (dtot.lt.1.0d-6) jan(i)=j
530      continue
520   continue
      if (idebug.eq.1) call basprt(iun3,.false.,.false.)
      if (idebug.eq.1) write(iun3,'(a)')'leave subroutine rdbasg'
      return

1000  istats = 0
      call inferr('ERROR reading GAMESS Basis set!',1)
      if (idebug.eq.1) write(iun3,'(a)')'leave subroutine rdbasg'
      return
9000  format(15x,i3,3x,a4,3x,i3,1x,2f15.6,18x,f15.6)
9001  format(1x,i3,3x,a2,i4,1x,f15.6,f12.6,14x,f12.6)
9002  format(1x,i3,3x,a2,i4,1x,f15.6,f18.12,20x,f18.12)
9003  format(1x,i3,3x,a2,i4,f20.7,2f20.12)
9004  format(1x,i3,3x,a2,i4,f22.7,2f20.12)
9005  format(1x,i6,3x,a2,i7,f22.7,2f18.12)
8000  format(15x,i3,3x,a4,3x,i3,1x,2f15.6)
8001  format(1x,i3,3x,a2,i4,1x,f15.6,f12.6)
8002  format(1x,i3,3x,a2,i4,1x,f15.6,f18.12)
      end

      subroutine rdmolg(istats)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs, nat(numatm)
      common /coord / xyz(3,numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line
c
c====== read molecular geometry ==================
c
cd     write(iun3,*)'enter subroutine rdmolg'
      istats = 1
      rewind iun2
      call search(line,'molecular geometry',istat)
      if (istat.eq.0) goto 1000
      call search(line,'atom   atomic',istat)
      if (istat.eq.0) goto 1000
      call readel(line,3)

         natoms = 0
5        read(iun2,'(a)',err=1000) line
         if (line(10:14).eq.'*****') then
cd           write(iun3,*)'leave subroutine rdmolg'
            return
         endif
         if (line(11:70).eq.' ') goto 5
         natoms = natoms + 1
         read(line,'(22x,f5.1,3(f12.7,3x),i5)',err=1000) charge,
     &   xyz(1,natoms),xyz(2,natoms),xyz(3,natoms),idum
         nat(natoms) = charge
         if (index(line,'bq').ne.0.or.index(line,'BQ').ne.0)
     &   natoms = natoms - 1
         goto 5

1000  istats = 0
cd     write(iun3,*)'leave subroutine rdmolg'
      call inferr('ERROR reading molecular geometry!',1)
      return
      end

      subroutine rdocc(imode,norbs,focc,eig,nocc,istat)
      implicit double precision (a-h,o-z), integer ( i-n)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line
      real eig
      dimension eig(*),focc(*)

      istat = 1

      if (imode.eq.0) then
         call searcht(line,'m.o. irrep        orbital',
     &                     'm.o.  irrep  orbital',
     &                     'm.o. irrep  orbital',istat)
         call readel(line,2)
      elseif (imode.eq.1) then
         call search(line,'alpha set',istat)
         if (istat.eq.0) call inferr('missing - alpha set - orbitals',1)
         iels = 7
         if (index(line,'--').ne.0) iels = 5
         call readel(line,iels)
      elseif (imode.eq.2) then
         call search(line,'beta set',istat)
         if (istat.eq.0) call inferr('missing - beta set - orbitals',1)
         iels = 7
         if (index(line,'--').ne.0) iels = 5
         call readel(line,iels)
      endif
      if (istat.eq.0) return
      istat = 1

      nocc = 0
      do i=1,norbs
         read(iun2,'(a)',err=1000) line
         if (line(39:39).eq.'.') then
            read(line,'(11x,f16.8,16x,f14.4)',err=1000) 
     &        eig(i),focc(i)
         else
            read(line,'(11x,f16.8,f20.7)',err=1000) eig(i),focc(i)
         endif
         if (focc(i).gt.0.0d0) nocc = nocc + 1
      end do

      return

1000  if (nocc.eq.0) istat = 0
      return
      end

      subroutine rdcasocc(ncore,nact,focc,nocc,istat)
      implicit double precision (a-h,o-z), integer ( i-n)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line
      dimension focc(*)

      istat = 1
      max=8

      nocc = 0
 10   if (nocc.lt.ncore+nact) then
         ibot=nocc+1
         itop=min0(nocc+max,ncore+nact)
         read(iun2,'(a)',err=1000) line
         read(line,'(10x,8f14.7)',err=1000) (focc(i),i=ibot,itop)
         read(iun2,'(a)',err=1000) line
         nocc=itop
         goto 10
      endif

      do i=1,ncore
         focc(i)=2.0d0
      enddo

      return

1000  if (nocc.eq.0) istat = 0
      return
      end

      subroutine rdgeom(numatm,coo,ianz,iatoms,istat)
      implicit double precision (a-h,p-w), integer (i-n), logical (o)
      dimension coo(3,numatm)
      dimension ianz(numatm)
      character*80 line
      character*2   tstr
      character*2   tocapf
      parameter (mxel=100)
      character*2   elemnt
      common /elem/elemnt(mxel)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
c
c     Read a GAMESS-UK geometry. If no geometry with in the current
c     format is found we'll try older formats.
c
c     We assume the first line we read contains something like
c     "x              y              z            chg  tag",
c     "atom             x              y              z",      or
c     "atom             x              y              z            chg"
c
c
cd     write(iun3,'(a)')'enter subroutine rdgeom'
      istat = 0
      iatoms = 0
      call readel(line,1)
      if (index(line,'tag').ne.0) then
c
c...     found a geometry in the current format
c
         call readel(line,1)
100      call readel(line,1)
         if (index(line,'====').eq.0) then
            iatoms = iatoms + 1
            read(line,'(2x,3f15.7,10x,a2)',err=900)
     &           coo(1,iatoms),coo(2,iatoms),coo(3,iatoms),tstr
            if (tstr(2:2).eq.' ') then
               tstr(2:2) = tstr(1:1)
               tstr(1:1) = ' '
            endif
            do jj=1,99
               if(tocapf(tstr).eq.tocapf(elemnt(jj))) ianz(iatoms) = jj
            enddo
            if (tocapf(tstr).eq.'BQ') iatoms=iatoms-1
            goto 100
         endif
      else 
c
c        Found one of old formats
c
c...     skip a potentially present empty line
c
         call readel(line,1)
         if (index(line,'.').ne.0) then
            backspace iun2
         elseif (index(line,'======').ne.0) then
            call readel(line,1)
         endif
c
c...     old formats the same from here on
c
200      call readel(line,1)
         if (index(line,'.').ne.0) then
            iatoms=iatoms+1
            read(line,'(20x,a2,8x,3f15.6)',err=900)
     &           tstr,coo(1,iatoms),coo(2,iatoms),coo(3,iatoms)
            if (tstr(2:2).eq.' ') then
               tstr(2:2) = tstr(1:1)
               tstr(1:1) = ' '
            endif
            do jj=1,99
               if(tocapf(tstr).eq.tocapf(elemnt(jj))) ianz(iatoms) = jj
            enddo
            if (tocapf(tstr).eq.'BQ') iatoms=iatoms-1
            goto 200
         endif
      endif
cd     write(iun3,'(a)')'leave subroutine rdgeom'
      return

900   istat=-1
cd     write(iun3,'(a)')'leave subroutine rdgeom'
      return
      end
