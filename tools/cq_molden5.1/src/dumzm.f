      subroutine dumzz(cc,ianc,natoms,
     &              bl,alph,bet,ibl,ialph,ibet,imap,iann,iz,
     &              c,cz,alpha,beta,ian,coo,iresid,issdon)

c this is really dumzm

      implicit double precision (a-h,o-z)
      logical ottest, oerror
      common /athlp/ iatoms, mxnat
      common /zmfrst/ ihaszm, nz, mxzat
      common /zmpart/ ipart,imn,imx,idcur
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt
      common /cllab/ iclon,iclpnt(4)
      dimension ian(*),c(3,*),cz(3,*),alpha(*),beta(*)
      dimension cc(3,*),ianc(*),imkeep(3)
      dimension bl(*),alph(*),bet(*),ibl(*),ialph(*),ibet(*),
     &          imap(*),iann(*),iz(4,*)
      dimension cm(3,3),dum(3)
      dimension coo(3,*),iresid(*)

c
c     Convert Z-matrix to coordinates
c
      if (ialtyp.eq.1) then
         nstrt = mxnat-nscnd

         do i=1,3
            imkeep(i) = imap(i) - iscst
            do j=1,3
               cm(j,i) = coo(j,nstrt + imkeep(i))
               cz(j,i) = coo(j,nstrt + imkeep(i))
            end do
         end do
      else
         do i=1,3
            imkeep(i) = imap(i)
            do j=1,3
               cm(j,i) = coo(j,imap(i))
               cz(j,i) = coo(j,imap(i))
            end do
         end do
      endif

      if (idcur.eq.1) call curs(1)
      maxnz = mxzat
      ottest=.true.

      call stoc(maxnz,nz,ipart,imn,imx,iann,iz,bl,alph,bet,
     &        ottest,jatoms,ian,c,cz,imap,alpha,beta,
     &        oerror,.true.,.true.)

      if (.not.oerror) then

         issdon = 0
         iclon = 0
         ihaszm = 1

         if (ialtyp.eq.1) then

            nscnd = jatoms
            nstrt = mxnat-nscnd
            natoms = iscst + nscnd

c        Determine whether imap refers to protein+lig or lig

            idorem = 1
            do i=1,3
               if (imap(i).gt.iscst) idorem = 0
            end do

            if (idorem.eq.1) then
               do i=1,jatoms
                  imap(i) = imap(i) + iscst
               end do
            endif

         else

c            if (ipart.eq.1) then
c               jatoms = natoms
c            elseif (ipart.eq.2) then
c
c distinction between ipart=1 and 2 lost, gives problems
c since because of fake sec.struc atoms, natoms also contain fake atoms
c
            if (ipart.ne.0) then
               jatoms = nz
               natoms = nz
            else
               natoms = jatoms
            endif
            nstrt = 0
         endif

         do i=1,jatoms

            do j=1,3
               cc(j,nstrt+i) = c(j,i)
            end do

            ianc(nstrt+i) = ian(i)

            if (ialtyp.eq.1) then
               ianc(iscst+i) = ian(i)
               iresid(iscst+i) = -4
            endif

         end do
         
         call rdmapf(cm,imkeep,ierr)

         if (ialtyp.eq.1) then
            call alnrot(dum,1)
         else
            call zm2fr(cm,c,imkeep)
         endif

      else

         call inferr('ERROR Zmat NOT parsed !!',1)

      endif

c set to whole zmatrix parsing next time around

      ipart = 0
      if (idcur.eq.1) call curs(0)

      return
      end

      subroutine prtzz(bl,alph,bet,ibl,ialph,ibet,imap,ianz,iz)

c this is really prtzm

      parameter (mxel=100)
      implicit double precision (a-h,o-z)
c
      common /zmfrst/ ihaszm, nz, mxzat

      character*2 elemnt
      common /elem/elemnt(mxel)
      dimension bl(*),alph(*),bet(*),ibl(*),ialph(*),ibet(*),
     &          imap(*),ianz(*),iz(4,*)

      print*,'nz=',nz
      do i=1,nz
         print*,elemnt(ianz(i)),iz(1,i),bl(i),iz(2,i),alph(i),iz(3,i),
     &          bet(i)
      end do

      return
      end

      subroutine chkmzz(istat,qat,rzp,ianzz,imap,ianz,lring,ityp,ipdbt)

c this is really chkmap

      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      common /zmfrst/ ihaszm, nz, mxzat
      integer*2 ityp,ipdbt
      common /types/ iff
      dimension ianzz(*),imap(*), qat(*), rzp(*)
      dimension ianz(*),ityp(*),ipdbt(*),lring(*)
c
c check if zmatrix ordering corresponds to coo ordering
c if not convert z-matrix to coo
c
c (mis)use lring array
c

      istat = 0

      if (ihaszm.eq.0) goto 100

      istat = 1

      if (nz.ne.iatoms) return

      do i=1,iatoms
        if (ianz(i).ne.ianzz(i)) goto 200
      end do

100   istat = 0
      return

200   if (iff.ne.0) then
         do i=1,iatoms
            lring(i) = ityp(i)
         end do
         do i=1,nz
            ityp(i) = lring(imap(i))
         end do
         do i=1,iatoms
            lring(i) = ipdbt(i)
         end do
         do i=1,nz
            ipdbt(i) = lring(imap(i))
         end do
         do i=1,iatoms
            rzp(i) = qat(i)
         end do
         do i=1,nz
            qat(i) = rzp(imap(i))
         end do
c have misused the rzp array, so we must restore it
         call qupd
      endif

      return
      end

      subroutine rdmapf(cm,imap,ierr)
      implicit double precision (a-h,o-z)
      logical opfil
      integer getlin
      character*137 line
      character*137 string
      common /curlin/ line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /gtmfil/ igtfil
      dimension cm(3,3),imap(3)
      dimension cmt(3,3),imapt(3)

c     igtfil.eq.1 triggered by a keypress 'c'

      if (igtfil.eq.0) return

      toang = 0.52917706d0
      ierr = 0
      iform = 1
      if (opfil(46,'mapfile',7,iform,1,1)) then
         print*,'opened file mapfile'
         iuntmp = iun2
         iun2 = 46
         do i=1,3
             ii = getlin(0)
             ktype = nxtwrd(string,nstr,itype,rtype)
             if (ktype.eq.2) then
                imapt(i) = itype
             else
                ierr = 1
                goto 100
             endif
             do j=1,3
                ktype = nxtwrd(string,nstr,itype,rtype)
                if (ktype.eq.3) then
                   cmt(j,i) = rtype/toang
                else
                   ierr = 1
                   goto 100
                endif
             end do
         end do
100      close(46)
         iun2 = iuntmp
         do i=1,3
             imap(i) = imapt(i)
             do j=1,3
                cm(j,i) = cmt(j,i)
             end do
         end do
         if (ierr.eq.1) print*,'error reading mapfile'
      else
         print*,'could read mapfile'
      endif

      return
      end

