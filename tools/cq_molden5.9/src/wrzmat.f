      subroutine allzmd(ipdbon,ianz,iaton,coo)
c
c This routine sets up the conversion of cartesian coordinates 
c to zmatrix, the actual converting is done by intzmt
c
      implicit double precision (a-h,o-z)
      common /athlp/  iatoms, mxnat
      common /zmfrst/ ihaszm, nz, mxzat
      common /zmatst/ nwrit,nvar,isimpl
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      logical opfil
      dimension coo(3,*),ianz(*),iaton(*)

      if (ihaszm.eq.0.and.ipdbon.eq.0) then
         do i=1,iatoms
             if (ianz(i).lt.100) then
                iaton(i) = 2
             else
                iaton(i) = 1
             endif
         end do
         if (isimpl.eq.3) then
            if (opfil(46,'mapfile',7,1,1,1)) then
               iuntmp = iun2
               iun2 = 46
               call getzm(iatoms,0,1,istat)
               if (istat.eq.1) then
                   ihaszm = 1
                   call dumzm(coo,ianz,iatoms)
               endif
               iun2 = iuntmp
               close(46)
            endif
            isimpl = 0
         else
            if (iatoms.gt.mxzat-100) then
c ask if you want a zmatrix of this molecule > 100 atoms
               call dozmt(istat)
               if (istat.eq.0) then
                  do i=1,iatoms
                      iaton(i) = 1
                  end do
                  return
               endif
            endif
            call intzmt(0)
         endif
         if (ihaszm.eq.1) ihaszm = 2
         call setarr(9,idum,idum)
      endif

      return
      end

      subroutine ligzmd(ianz,iaton)
c
c this a dedicated version of allzmt, it only converts the ligand
c part of a receptor - ligand complex to a zmatrix, used by manual
c docking
c
      implicit double precision (a-h,o-z)
      common /zmfrst/ ihaszm, nz, mxzat
      common /align/  vecs(3,3),nscnd,iscst,ialtyp,iocnt
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension ianz(*),iaton(*)

      do i=1,iscst
          iaton(i) = 1
      end do

      do i=iscst+1,iscst+nscnd
          if (ianz(i).lt.100) then
             iaton(i) = 2
          else
             iaton(i) = 1
          endif
      end do

      call intzmt(0)
      if (ihaszm.eq.1) ihaszm = 2

      do i=iscst+1,iscst+nscnd
          iaton(i) = 1
      end do

      return
      end

      subroutine pdbzmd(ianz,iaton,iresid,iconn,ishoh)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/  iatoms, mxnat
      common /zmfrst/ ihaszm, nz, mxzat
      dimension ianz(*),iaton(*),iresid(*),iconn(mxcon+1,*)

      if (ihaszm.eq.0) then

         ino = 0
         do i=1,iatoms
             if (ianz(i).lt.100) then
                if (iresid(i).eq.-ishoh) then
                   if (ianz(i).eq.8) then
                      ido = 0
                      do j=1,iconn(1,i)
                         jj = iconn(1+j,i)
                         if (jj.gt.0.and.ianz(iabs(jj)).gt.1) ido = 1
                      end do
                      if (ido.eq.1) then
                         iaton(i) = 2
                      else
                         ino = ino + 1
                         iaton(i) = 1
                      endif
                   elseif (ianz(i).eq.1) then
                      ido = 0
                      do j=1,iconn(1,i)
                         jj = iconn(1+j,i)
                         if (jj.gt.0) then
                            if (ianz(iabs(jj)).ne.8) then
                               ido = 1
                            else
                               io = iabs(jj)
                               do ll=1,iconn(1,io)
                                  iio = iconn(1+ll,io)
                                  if (iio.gt.0.and.iio.ne.i) then
                                     if (ianz(iio).ne.1) ido = 1
                                  endif
                               end do
                            endif
                         endif
                      end do
                      if (ido.eq.1) then
                         iaton(i) = 2
                      else
                         ino = ino + 1
                         iaton(i) = 1
                      endif
                   else
                      iaton(i) = 2
                   endif
                else
                   iaton(i) = 2
                endif
             else
                iaton(i) = 1
             endif
         end do

         call stowat(ino)

         call intzmt(1)
         ihaszm = 1

         do i=1,iatoms
             iaton(i) = 1
         end do

      endif

      return
      end

      subroutine haswad(ino,ianz,iaton,iresid,iconn)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/  iatoms, mxnat
      common /zmfrst/ ihaszm, nz, mxzat
      dimension ianz(*),iaton(*),iresid(*),iconn(mxcon+1,*)

      ino = 0

      if (ihaszm.eq.0) then

         do i=1,iatoms
             if (ianz(i).lt.100) then
                if (iresid(i).le.0) then

                   if (ianz(i).eq.8) then

                      ido = 0

                      do j=1,iconn(1,i)
                         jj = iconn(1+j,i)
                         if (jj.gt.0.and.ianz(iabs(jj)).gt.1) ido = 1
                      end do

                      if (ido.ne.1) ino = ino + 1

                   elseif (ianz(i).eq.1) then

                      ido = 0

                      do j=1,iconn(1,i)
                         jj = iconn(1+j,i)
                         if (jj.gt.0) then
                            if (ianz(iabs(jj)).ne.8) then
                               ido = 1
                            else
                               io = iabs(jj)
                               do ll=1,iconn(1,io)
                                  iio = iconn(1+ll,io)
                                  if (iio.gt.0.and.iio.ne.i) ido = 1
                               end do
                            endif
                         endif
                      end do

                      if (ido.ne.1) ino = ino + 1

                   endif
                endif
             endif
         end do
      endif

      return
      end

      subroutine haszm(zmok)
      implicit double precision (a-h,o-z)
      common /zmfrst/ ihaszm, nz, mxzat

      logical zmok

      if (zmok) then
         ihaszm = 1
      else
         ihaszm = 0
         nz = 0
      endif
      
      return
      end

      logical function zmqok(idum)
      implicit double precision (a-h,o-z)
      common /zmfrst/ ihaszm, nz, mxzat
      character*80 glin1,glin2,gtitl,rungam
      character*15 jname,qname
      common /gauopt/ ito,imo,ibo,itotc,imult,ibatch,ihess,itime,
     &                iwxyz,ichh,ichm,ichl,imh,imm,iml,iexk,
     &                glin1,glin2,gtitl,jname,qname,rungam

      zmqok = .false.
      if (ihaszm.ne.0.or.iwxyz.ne.0) zmqok = .true.
      
      return
      end

      subroutine setzm(newnz)
      implicit double precision (a-h,o-z)
      common /zmfrst/ ihaszm, nz, mxzat

      if (ihaszm.ne.1) then
         ihaszm = 1
         nz = newnz
      endif
      
      return
      end

      logical function chkres(ires,icalf,iamino)
      implicit double precision (a-h,o-z)
      dimension icalf(6,*),iamino(*)

      chkres = .true.

      nat = 3
      if (iamino(ires).gt.23) nat = 6

      do i=1,nat
         do j=i+1,nat
            if (icalf(i,ires).eq.icalf(j,ires).or.icalf(i,ires).le.0)
     &      chkres = .false.
         end do
      end do

      return
      end

      logical function addiz(ipdb,k1,k2,l1,l2,lwrit,lring)
      implicit double precision (a-h,o-z)
      parameter (mxzz=115)
      common /pdbzm/ izz(4,mxzz)
      logical pcklin,dorng
      dimension isel(4)
      dimension ipdb(*),lwrit(*),lring(*)


      addiz = .true.
      if (k1.gt.mxzz.or.k2.gt.mxzz) return

      do k=k1,k2
         dorng = .false.
         if (l1.ne.0.and.k.ge.l1.and.k.le.l2) dorng = .true.
         do l=1,4
            it = ipdb(izz(l,k))
            if (it.eq.0) return
            if (l.eq.1) then
               it1 = it
            else
               if (lwrit(it).eq.0) return
            endif
            isel(l) = it
         end do
         
c         if (.not.pcklin(isel)) goto 100
         if (.not.pcklin(isel)) return
         if (dorng) lring(it1) = 1
      end do

      return

100   addiz = .false.
      return

      end

      logical function addone(i1,i2,i3,i4,iadd,lwrit)
      implicit double precision (a-h,o-z)
      logical pcklin
      dimension isel(4),lwrit(*)

      addone = .true.

      isel(1) = i1
      isel(2) = i2
      isel(3) = i3
      isel(4) = i4

      do l=1,4
         if (isel(l).eq.0) return
         if (l.gt.1) then
            if (lwrit(isel(l)).eq.0) return
         endif
      end do
         
      if (.not.pcklin(isel)) goto 100

      iadd = iadd + 1

      return

100   addone = .false.
      return

      end

      subroutine intzmd(ispdb,ianz,iaton,iresid,iconn,lwrit,lring,
     &                  icalf,ianf,islu,nchain,iamino)
      implicit double precision (a-h,o-z)
      parameter (ncmax=1000)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxzz=115)

      common /athlp/  iatoms, mxnat
      common /selatm/ jring(numat1)
      common /zmatst/ nwrit,nvar,isimpl
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /surf/   natorg,noscnd
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /pdbzm/ izz(4,mxzz)

      logical ocnos,parlea,pcklin,ringg,wring,debug,opfil,chkres,
     &        addiz,addone,parleh
      integer getlin
      character*137 str
      dimension isel(4),icalft(6),iring(6),idisc(ncmax)
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension ianz(*),iaton(*),iresid(*),iconn(mxcon+1,*),
     &          lwrit(*),lring(*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*)

      data ((izz(i,j),i=1,4),j=1,92) /
     &          49,48,47,46, 54,49,48,47, 52,54,49,48, 53,52,54,49,
     &          73,53,52,54, 44,43,46,47, 45,43,46,47,
     &          66,54,49,48, 59,66,54,49, 65,59,66,54, 57,65,59,66,
     &          56,57,65,59, 62,56,57,65, 55,62,56,57, 60,55,62,56,
     &          58,60,55,62, 60,54,49,48, 58,60,54,49, 57,58,60,54,
     &          56,57,58,60, 62,56,57,58, 55,62,56,57, 67,55,62,56,
     &          63,56,57,58, 68,56,57,58, 72,55,58,60, 83,57,58,60,
     &          64,58,57,65, 69,58,57,65, 61,55,60,58, 31, 5, 2, 1,
     &          37, 5, 2, 1, 32, 5, 2, 1,  8, 5, 2,32,  7, 5, 2, 1,
     &           8, 5, 2, 7, 10, 7, 5, 2,  7, 5, 2, 1,  8, 5, 2, 7,
     &           6, 5, 2, 1, 36, 6, 5, 2, 12,36, 6, 5,  6, 5, 2, 1,
     &          29, 6, 5, 2, 30, 6, 5,29,  6, 5, 2, 1, 29, 6, 5, 2,
     &          21, 6, 5,29,  6, 5, 2, 1, 10, 6, 5, 2, 11, 6, 5,10,
     &           6, 5, 2, 1,  9, 6, 5, 2, 12, 9, 6, 5, 27,12, 9, 6,
     &           6, 5, 2, 1,  9, 6, 5, 2, 34, 9, 6, 5, 35, 9, 6,34,
     &           6, 5, 2, 1,  9, 6, 5, 2, 34, 9, 6, 5, 24, 9, 6,34,
     &           6, 5, 2, 1,  9, 6, 5, 2,  6, 5, 2, 1,  9, 6, 5, 2,
     &          22, 9, 6, 5, 17,22, 9, 6, 25,17,22, 9, 26,17,22,25,
     &           6, 5, 2, 1, 20, 6, 5, 2, 13,20, 6, 5, 24,13,20, 6,
     &          11,24,13,20,  6, 5, 2, 1, 10, 6, 5, 2, 13,10, 6, 5,
     &          17,13,10, 6, 14,17,13,10, 11,14,17,13, 33,17,13,10,
     &           6, 5, 2, 1, 10, 6, 5, 2, 23,10, 6, 5, 14,23,10, 6,
     &          11,14,23,10, 15,11,14,23, 19,15,11,14, 16,19,15,11,
     &          18,16,19,15/

      data ((izz(i,j),i=1,4),j=93,98) /
     &          57,54,49,48, 58,57,54,49, 60,58,57,54,
     &          55,60,58,57, 62,55,60,58, 56,62,55,60/

      data ((izz(i,j),i=1,4),j=99,115) /
     &          90,61,55,60, 91,90,61,55, 92,91,90,61,
     &          93,92,91,90, 94,93,92,91, 95,94,93,92, 96,95,93,92,
     &          97,95,94,96, 98,97,95,96, 99,94,93,95, 100,99,94,93,
     &          101,100,99,94, 102,100,99,101, 103,102,100,101,
     &          88,62,55,60, 89,90,61,55, 76,43,46,47/

      debug = .false.
      if (debug) print*,'routine intzmt'
c
c     Convert Cartesian coordinates to Z-matrix
c
      call haszm(.true.)

      natorg = 0

      nvar = 0
      ndisc = 0
      nwrit = 0
      nwrith = 0

      noth = 0
      nsel = 0
      nconn = 0
      nzpdb = 0

      ilead = 0
      ilold = 0

      do i=1,iatoms
         lring(i) = 0
         lwrit(i) = 0
         if (ianz(i).ne.1.and.iaton(i).eq.2) noth = noth + 1
         if (iaton(i).eq.2) nsel = nsel + 1
         if (iaton(i).eq.2.and.iconn(1,i).ne.0) nconn = nconn + 1
      end do

      if (nconn.lt.3) then
         do i=1,iatoms
            if (iaton(i).eq.2.and.iconn(1,i).eq.0) then
                lwrit(i) = 1
                call icrcon(icrcn,isel,idisc,ndisc,1,0)
                if (icrcn.ne.0) then
                   jj = isel(1)
                   lwrit(jj) = 1
                   call icrcon(icrcn,isel,idisc,ndisc,1,0)
                   if (icrcn.ne.0) then
                       lwrit(i) = 0
                       lwrit(jj) = 0
                       goto 5
                   else 
                       lwrit(i) = 0
                       lwrit(jj) = 0
                   endif
                else
                   lwrit(i) = 0
                endif
            endif
         end do
      endif

      if (ispdb.eq.1) then

         ist = 1
         do i=1,6
           icalft(i) = 0
         end do
         irst = 0

         do i=1,nchain
           iaw = 0
           do j=ianf(i),islu(i)
           
            if (chkres(j,icalf,iamino)) then
             call getpdb(j,ipdb,ihpdb)
c
c Amino Acids
c
             if (iamino(j).le.23) then
c n (first residue)
               if (ist.eq.1) then
                  isel(1) = icalf(2,1)
                  isel(2) = 0
                  isel(3) = 0
                  isel(4) = 0
                  if (.not.pcklin(isel)) goto 100
               elseif (iaw.eq.0) then
                  isel(1) = icalf(2,j)
                  if (irst.gt.23) then
                    isel(2) = icalft(6)
                    isel(3) = icalft(5)
                    isel(4) = icalft(4)
                  else
                    isel(2) = icalft(3)
                    isel(3) = icalft(1)
                    isel(4) = icalft(2)
                  endif
                  if (.not.pcklin(isel)) goto 100
               endif
c ca
               isel(1) = icalf(1,j)
               isel(2) = icalf(2,j)
               if (ist.eq.1) then
                  isel(3) = 0
                  isel(4) = 0
               else
                  if (iaw.eq.0) then
                     if (irst.gt.23) then
                       isel(3) = icalft(6)
                       isel(4) = icalft(5)
                     else
                       isel(3) = icalft(3)
                       isel(4) = icalft(1)
                     endif
                     ifd = 0 
                     do k=1,iconn(1,isel(2))
                       kk = abs(iconn(k+1,isel(2)))
                       if (kk.eq.isel(3)) ifd = 1
                     end do
                     if (ifd.eq.0) 
     &               idum = jcrcon(isel(2),isel(3),idisc,ndisc,iconn)
                  else
                     if (iamino(j-1).gt.23) then
                       isel(3) = icalf(6,j-1)
                       isel(4) = icalf(5,j-1)
                     else
                       isel(3) = icalf(3,j-1)
                       isel(4) = icalf(1,j-1)
                     endif
                  endif
               endif
               if (.not.pcklin(isel)) goto 100
               if (iamino(j).eq.15) then
                  if (lwrit(icalf(1,j)).ne.0) lring(icalf(1,j)) = 1
               endif
c co
               isel(1) = icalf(3,j)
               isel(2) = icalf(1,j)
               isel(3) = icalf(2,j)
               if (ist.eq.1) then
                  isel(4) = 0
               else
                  if (iaw.eq.0) then
                     if (irst.gt.23) then
                       isel(4) = icalft(6)
                     else
                       isel(4) = icalft(3)
                     endif
                  else
                     if (iamino(j-1).gt.23) then
                       isel(4) = icalf(6,j-1)
                     else
                       isel(4) = icalf(3,j-1)
                     endif
                  endif
               endif
               if (.not.pcklin(isel)) goto 100
c n or oxt
               if (ipdb(38).ne.0) then
                  inxt = ipdb(38)
               elseif (ipdb(76).ne.0) then
                  inxt = ipdb(76)
               else
                  if (j.ne.islu(i)) then
                     if (iamino(j+1).gt.23) then
                       inxt = icalf(1,j+1)
                     else
                       inxt = icalf(2,j+1)
                     endif
                  else
                     inxt = 0
                  endif
               endif

               if (inxt.ne.0) then
                  isel(1) = inxt
                  isel(2) = icalf(3,j)
                  isel(3) = icalf(1,j)
                  isel(4) = icalf(2,j)
                  if (.not.pcklin(isel)) goto 100
                  if (iamino(j+1).eq.15) then
                     if (inxt.eq.icalf(2,j+1)) then
                        if (lwrit(icalf(2,j+1)).ne.0)
     &                     lring(icalf(2,j+1)) = 1
                     endif
                  endif
               endif
c o
               isel(1) = ipdb(4)
               isel(2) = icalf(3,j)
               isel(3) = icalf(1,j)
               if (inxt.eq.0) then
                  isel(4) = icalf(2,j)
               else
                  isel(4) = inxt
               endif
               if (.not.pcklin(isel)) goto 100

               do k=1,6
                 icalft(k) = icalf(k,j)
               end do
               irst = iamino(j)
               iaw = 1
c nh
               
               if (ist.eq.1) then
                  is = icalf(3,j)
               else
                  is = icalf(3,j-1)
               endif
               if (.not.addone(icalf(4,j),icalf(2,j),icalf(1,j),
     &             is,nwrith,lwrit)) goto 100

c c beta
               if (.not.addone(ipdb(5),icalf(1,j),icalf(3,j),
     &             icalf(2,j),idum,lwrit)) goto 100
               if (iamino(j).eq.15) then
                   if (ipdb(5).ne.0) then
                      if (lwrit(ipdb(5)).ne.0) lring(ipdb(5)) = 1
                   endif
               endif

c ser
               if (.not.addiz(ipdb,31,31,0,0,lwrit,lring)) goto 100
c cys
               if (.not.addiz(ipdb,32,32,0,0,lwrit,lring)) goto 100
c thr
               if (iamino(j).eq.5) then
                  if (.not.addiz(ipdb,33,34,0,0,lwrit,lring)) goto 100
               endif
c ile
               if (iamino(j).eq.6) then
                  if (.not.addiz(ipdb,35,37,0,0,lwrit,lring)) goto 100
               endif
c val
               if (iamino(j).eq.7) then
                  if (.not.addiz(ipdb,38,39,0,0,lwrit,lring)) goto 100
               endif
c met
               if (iamino(j).eq.8) then
                  if (.not.addiz(ipdb,40,42,0,0,lwrit,lring)) goto 100
               endif
c asp
               if (iamino(j).eq.9) then
                  if (.not.addiz(ipdb,43,45,0,0,lwrit,lring)) goto 100
               endif
c asn
               if (iamino(j).eq.10) then
                  if (.not.addiz(ipdb,46,48,0,0,lwrit,lring)) goto 100
               endif
c leu
               if (iamino(j).eq.11) then
                  if (.not.addiz(ipdb,49,51,0,0,lwrit,lring)) goto 100
               endif
c lys
               if (iamino(j).eq.12) then
                  if (.not.addiz(ipdb,52,55,0,0,lwrit,lring)) goto 100
               endif
c glu
               if (iamino(j).eq.13) then
                  if (.not.addiz(ipdb,56,59,0,0,lwrit,lring)) goto 100
               endif
c gln
               if (iamino(j).eq.14) then
                  if (.not.addiz(ipdb,60,63,0,0,lwrit,lring)) goto 100
               endif
c pro
               if (iamino(j).eq.15) then
                  if (.not.addiz(ipdb,64,65,64,65,lwrit,lring)) goto 100
               endif
c lys
               if (iamino(j).eq.16) then
                  if (.not.addiz(ipdb,66,71,0,0,lwrit,lring)) goto 100
               endif
c his
               if (iamino(j).eq.17) then
                  if (.not.addiz(ipdb,72,76,72,76,lwrit,lring)) goto 100
               endif
c phe
               if (iamino(j).eq.18) then
                  if (.not.addiz(ipdb,77,82,77,82,lwrit,lring)) goto 100
               endif
c tyr
               if (iamino(j).eq.19) then
                  if (.not.addiz(ipdb,77,83,77,82,lwrit,lring)) goto 100
               endif
c trp
               if (iamino(j).eq.20) then
                  if (.not.addiz(ipdb,84,92,84,92,lwrit,lring)) goto 100
               endif

               if (.not.addone(ihpdb(4),icalf(1,j),icalf(3,j),
     &             icalf(2,j),nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(5),icalf(1,j),icalf(3,j),
     &             icalf(2,j),nwrith,lwrit)) goto 100

               do m=6,mxhsym*3
                  l = ihpdb(m)
                  if (l.ne.0) then
                     if (ianz(l).eq.1.and.lwrit(l).eq.0
     &                    .and.iaton(l).eq.2) then
                        if (.not.parleh(l,isel,0)) then
                           if (ifcon(l,idisc,ndisc,0,ianz,iaton,iconn,
     &                               lwrit) .ne.0) then
                              if (.not.parleh(l,isel,0)) idum = 0
                           endif
                        endif
                        if (lwrit(l).eq.1) nwrith = nwrith + 1
                     endif
                  endif
               end do

c first residue more n-h ?

               if (ist.eq.1) then
                  do l=1,3
                     if (ihpdb(l).ne.icalf(4,j).and.ihpdb(l).ne.0) 
     &               then
                        if (.not.addone(ihpdb(l),icalf(2,j),
     &                     icalf(1,j),icalf(4,j),nwrith,lwrit)) goto 100
                     endif
                  end do
               endif

             else
c
c Nucleic Acids
c

c P
               isel(1) = icalf(1,j)
               if (ist.eq.1) then
                  isel(2) = 0
                  isel(3) = 0
                  isel(4) = 0
               else
                  if (iaw.eq.0) then
                     if (irst.le.23) then
                       isel(2) = icalft(3)
                       isel(3) = icalft(1)
                       isel(4) = icalft(2)
                     else
                       isel(2) = icalft(6)
                       isel(3) = icalft(5)
                       isel(4) = icalft(4)
                     endif
                     ifd = 0 
                     do k=1,iconn(1,isel(1))
                       kk = abs(iconn(k+1,isel(1)))
                       if (kk.eq.isel(2)) ifd = 1
                     end do
                     if (ifd.eq.0) 
     &               idum = jcrcon(isel(1),isel(2),idisc,ndisc,iconn)
                  else
                     if (iamino(j-1).le.23) then
                        isel(2) = icalf(3,j-1)
                        isel(3) = icalf(1,j-1)
                        isel(4) = icalf(2,j-1)
                     else
                        isel(2) = icalf(6,j-1)
                        isel(3) = icalf(5,j-1)
                        isel(4) = icalf(4,j-1)
                     endif
                  endif
               endif
               if (.not.pcklin(isel)) goto 100
c O5*
               isel(1) = icalf(2,j)
               isel(2) = icalf(1,j)
               if (ist.eq.1) then
                  isel(3) = 0
                  isel(4) = 0
               else
                  if (iaw.eq.0) then
                     if (irst.le.23) then
                        isel(3) = icalft(3)
                        isel(4) = icalft(1)
                     else
                        isel(3) = icalft(6)
                        isel(4) = icalft(5)
                     endif
                  else
                     if (iamino(j-1).le.23) then
                        isel(3) = icalf(3,j-1)
                        isel(4) = icalf(1,j-1)
                     else
                        isel(3) = icalf(6,j-1)
                        isel(4) = icalf(5,j-1)
                     endif
                  endif
               endif
               if (.not.pcklin(isel)) goto 100
c C5*
               isel(1) = icalf(3,j)
               isel(2) = icalf(2,j)
               isel(3) = icalf(1,j)
               if (ist.eq.1) then
                  isel(4) = 0
               else
                  if (iaw.eq.0) then
                     if (irst.le.23) then
                        isel(4) = icalft(3)
                     else
                        isel(4) = icalft(6)
                     endif
                  else
                     if (iamino(j-1).le.23) then
                        isel(4) = icalf(3,j-1)
                     else
                        isel(4) = icalf(6,j-1)
                     endif
                  endif
               endif
               if (.not.pcklin(isel)) goto 100
c C4*
               isel(1) = icalf(4,j)
               isel(2) = icalf(3,j)
               isel(3) = icalf(2,j)
               isel(4) = icalf(1,j)
               if (.not.pcklin(isel)) goto 100
               if (lwrit(isel(1)).ne.0) lring(isel(1)) = 1
c C3*
               isel(1) = icalf(5,j)
               isel(2) = icalf(4,j)
               isel(3) = icalf(3,j)
               isel(4) = icalf(2,j)
               if (.not.pcklin(isel)) goto 100
               if (lwrit(isel(1)).ne.0) lring(isel(1)) = 1
c O3*
               isel(1) = icalf(6,j)
               isel(2) = icalf(5,j)
               isel(3) = icalf(4,j)
               isel(4) = icalf(3,j)
               if (.not.pcklin(isel)) goto 100

               do k=1,6
                 icalft(k) = icalf(k,j)
               end do
               irst = iamino(j)
               iaw = 1
c O4*
               if (.not.addiz(ipdb,1,1,1,1,lwrit,lring)) goto 100
c C1*
               if (.not.addiz(ipdb,2,2,2,2,lwrit,lring)) goto 100
c C2*
               if (.not.addiz(ipdb,3,3,3,3,lwrit,lring)) goto 100
c O2*
               if (.not.addiz(ipdb,4,4,0,0,lwrit,lring)) goto 100
               if (.not.addiz(ipdb,5,5,0,0,lwrit,lring)) goto 100

c O1P, O2P
               if (.not.addiz(ipdb,6,6,0,0,lwrit,lring)) goto 100
               if (.not.addiz(ipdb,7,7,0,0,lwrit,lring)) goto 100
c OXT
               if (addiz(ipdb,115,115,0,0,lwrit,lring)) then
                   ldum = 1
               endif
               

               k1 = 0
c Adenosine Guanosine
               if (iamino(j).eq.24.or.iamino(j).eq.26.or.
     &             iamino(j).eq.29.or.(iamino(j).ge.32.and.
     &             iamino(j).le.38)) then
                  k1 = 8
                  k2 = 16
                  l1 = 8
                  l2 = 16
               endif
c cytidine Uridine Thymine
               if (iamino(j).eq.25.or.iamino(j).eq.27.or.
     &             iamino(j).eq.28.or.iamino(j).eq.30.or.
     &             iamino(j).eq.31.or.(iamino(j).ge.39.and.
     &             iamino(j).le.41)) then
                  k1 = 17
                  k2 = 23
                  l1 = 17
                  l2 = 22
               endif
c PSU 
               if (iamino(j).eq.42) then
                  k1 = 93
                  k2 = 98
                  l1 = 93
                  l2 = 98
               endif

               if (k1.ne.0) then
                  if (.not.addiz(ipdb,k1,k2,l1,l2,lwrit,lring)) goto 100
               endif

               if (iamino(j).eq.42) then
                   if (.not.addiz(ipdb,23,23,0,0,lwrit,lring)) goto 100
               endif

               if (.not.addiz(ipdb,24,24,0,0,lwrit,lring)) goto 100
               if (.not.addiz(ipdb,25,25,0,0,lwrit,lring)) goto 100
               if (.not.addiz(ipdb,26,26,0,0,lwrit,lring)) goto 100
               if (.not.addiz(ipdb,27,27,0,0,lwrit,lring)) goto 100
               if (.not.addiz(ipdb,28,28,0,0,lwrit,lring)) goto 100
               if (.not.addiz(ipdb,29,29,0,0,lwrit,lring)) goto 100
               if (.not.addiz(ipdb,30,30,0,0,lwrit,lring)) goto 100

c - YG
               if (iamino(j).eq.37) then
                  k1 = 99
                  k2 = 114
                  l1 = 99
                  l2 = 114
                  if (.not.addiz(ipdb,k1,k2,l1,l2,lwrit,lring)) goto 100
               endif
c - 1MA
               if (iamino(j).eq.29.or.iamino(j).eq.32) then
                  if (.not.addone(ipdb(79),ipdb(60),ipdb(58),
     &                ipdb(57),nwrith,lwrit)) goto 100
               endif

c - OMC, OMG
               if (iamino(j).eq.31.or.iamino(j).eq.36) then
                  if (.not.addone(ipdb(80),ipdb(53),ipdb(52),
     &                ipdb(54),nwrith,lwrit)) goto 100
               endif

c - 2MG, M2G
               if (iamino(j).eq.33.or.iamino(j).eq.34) then
                  if (.not.addone(ipdb(80),ipdb(61),ipdb(55),
     &                ipdb(60),nwrith,lwrit)) goto 100
               endif

c - M2G
               if (iamino(j).eq.34) then
                  if (.not.addone(ipdb(79),ipdb(61),ipdb(55),
     &                ipdb(60),nwrith,lwrit)) goto 100
               endif

c - 7MG
               if (iamino(j).eq.35) then
                  if (.not.addone(ipdb(85),ipdb(65),ipdb(59),
     &                ipdb(66),nwrith,lwrit)) goto 100
               endif

c backbone H
               if (.not.addone(ihpdb(61),ipdb(47),ipdb(46),ipdb(43),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(62),ipdb(47),ipdb(46),ipdb(43),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(64),ipdb(48),ipdb(47),ipdb(46),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(67),ipdb(50),ipdb(48),ipdb(47),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(115),ipdb(51),ipdb(50),ipdb(52),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(70),ipdb(52),ipdb(50),ipdb(48),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(71),ipdb(52),ipdb(50),ipdb(48),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(73),ipdb(54),ipdb(52),ipdb(50),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(76),ipdb(53),ipdb(52),ipdb(50),
     &            nwrith,lwrit)) goto 100

c cytidine
c 101 and 104 are for H2U

               if (.not.addone(ihpdb(130),ipdb(63),ipdb(56),ipdb(57),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(131),ipdb(63),ipdb(56),ipdb(57),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(100),ipdb(57),ipdb(58),ipdb(56),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(101),ipdb(57),ipdb(58),ipdb(56),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(103),ipdb(58),ipdb(60),ipdb(57),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(104),ipdb(58),ipdb(60),ipdb(57),
     &            nwrith,lwrit)) goto 100

c Guanosine
               if (.not.addone(ihpdb(112),ipdb(59),ipdb(66),ipdb(65),
     &            nwrith,lwrit)) goto 100
               if (iamino(j).eq.35) then
                  if (.not.addone(ihpdb(113),ipdb(59),ipdb(66),
     &                ipdb(65),nwrith,lwrit)) goto 100
               endif
               if (.not.addone(ihpdb(121),ipdb(60),ipdb(58),ipdb(57),
     &            nwrith,lwrit)) goto 100
               if (.not.iamino(j).eq.37) then
                  if (.not.addone(ihpdb(124),ipdb(61),ipdb(55),
     &                ipdb(60),nwrith,lwrit)) goto 100
               endif
               if (.not.addone(ihpdb(125),ipdb(61),ipdb(55),ipdb(60),
     &            nwrith,lwrit)) goto 100
c Adenosine
               if (.not.addone(ihpdb(82),ipdb(55),ipdb(60),ipdb(58),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(136),ipdb(64),ipdb(58),ipdb(57),
     &            nwrith,lwrit)) goto 100
               if (.not.addone(ihpdb(137),ipdb(64),ipdb(58),ipdb(57),
     &            nwrith,lwrit)) goto 100
c Uridine
               if (.not.addone(ihpdb(127),ipdb(62),ipdb(56),ipdb(57),
     &            nwrith,lwrit)) goto 100

c Thymidine CM5 (= C5A,C5M) hydrogens, 5MC, 5MU H

               if (iamino(j).eq.27.or.iamino(j).eq.30
     &             .or.iamino(j).eq.41) then
                  if (ipdb(83).ne.0) then
                     if (.not.addone(ihpdb(160),ipdb(83),ipdb(57),
     &                   ipdb(56),nwrith,lwrit)) goto 100
                     if (.not.addone(ihpdb(161),ipdb(83),ipdb(57),
     &                   ihpdb(160),nwrith,lwrit)) goto 100
                     if (.not.addone(ihpdb(162),ipdb(83),ipdb(57),
     &                   ihpdb(160),nwrith,lwrit)) goto 100
                  elseif (ipdb(85).ne.0) then
                     if (.not.addone(ipdb(85),ipdb(57),ipdb(58),
     &                   ipdb(60),nwrith,lwrit)) goto 100
                     if (.not.addone(ihpdb(166),ipdb(85),ipdb(57),
     &                   ipdb(56),nwrith,lwrit)) goto 100
                     if (.not.addone(ihpdb(167),ipdb(85),ipdb(57),
     &                   ihpdb(160),nwrith,lwrit)) goto 100
                     if (.not.addone(ihpdb(168),ipdb(85),ipdb(57),
     &                   ihpdb(160),nwrith,lwrit)) goto 100
                  endif
               endif

c - 1MA H
               if (iamino(j).eq.29.or.iamino(j).eq.32) then
                  if (.not.addone(ihpdb(148),ipdb(79),ipdb(60),
     &                ipdb(58),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(149),ipdb(79),ipdb(60),
     &                ihpdb(148),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(150),ipdb(79),ipdb(60),
     &                ihpdb(148),nwrith,lwrit)) goto 100
c    H6 already covered
               endif

c - OMC, OMG H
               if (iamino(j).eq.31.or.iamino(j).eq.36) then
                  if (.not.addone(ihpdb(151),ipdb(80),ipdb(53),
     &                ipdb(52),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(152),ipdb(80),ipdb(53),
     &                ihpdb(151),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(153),ipdb(80),ipdb(53),
     &                ihpdb(151),nwrith,lwrit)) goto 100
               endif

c - 2MG H
               if (iamino(j).eq.33.or.iamino(j).eq.34) then
                  if (.not.addone(ihpdb(151),ipdb(80),ipdb(61),
     &                ipdb(55),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(152),ipdb(80),ipdb(61),
     &                ihpdb(151),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(153),ipdb(80),ipdb(61),
     &                ihpdb(151),nwrith,lwrit)) goto 100
               endif

c - M2G H
               if (iamino(j).eq.34) then
                  if (.not.addone(ihpdb(148),ipdb(79),ipdb(61),
     &                ipdb(55),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(149),ipdb(79),ipdb(61),
     &                ihpdb(148),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(150),ipdb(79),ipdb(61),
     &                ihpdb(148),nwrith,lwrit)) goto 100
               endif

c - 7MG H
               if (iamino(j).eq.35) then
                  if (.not.addone(ihpdb(166),ipdb(85),ipdb(65),
     &                ipdb(59),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(167),ipdb(85),ipdb(65),
     &                ihpdb(166),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(168),ipdb(85),ipdb(65),
     &                ihpdb(166),nwrith,lwrit)) goto 100
c the others ones are covered under guanosine
               endif

c - YG H
               if (iamino(j).eq.37) then
                  if (.not.addone(ihpdb(91),ipdb(88),ipdb(62),
     &                ipdb(56),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(92),ipdb(88),ipdb(62),
     &                ihpdb(91),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(93),ipdb(88),ipdb(62),
     &                ihpdb(91),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(175),ipdb(89),ipdb(90),
     &                ipdb(61),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(176),ipdb(89),ipdb(90),
     &                ihpdb(175),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(177),ipdb(89),ipdb(90),
     &                ihpdb(175),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(178),ipdb(92),ipdb(91),
     &                ipdb(93),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(179),ipdb(92),ipdb(91),
     &                ihpdb(178),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(181),ipdb(93),ipdb(92),
     &                ipdb(94),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(182),ipdb(93),ipdb(92),
     &                ihpdb(181),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(184),ipdb(94),ipdb(93),
     &                ipdb(95),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(187),ipdb(98),ipdb(97),
     &                ipdb(95),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(188),ipdb(98),ipdb(97),
     &                ihpdb(187),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(189),ipdb(98),ipdb(97),
     &                ihpdb(187),nwrith,lwrit)) goto 100
c 124 = HN2, this should actually be HN20, (or else H20)
c also wrong at routine addhs
                  if (.not.addone(ihpdb(124),ipdb(99),ipdb(100),
     &                ipdb(94),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(190),ipdb(103),ipdb(102),
     &                ipdb(100),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(191),ipdb(103),ipdb(102),
     &                ihpdb(190),nwrith,lwrit)) goto 100
                  if (.not.addone(ihpdb(192),ipdb(103),ipdb(102),
     &                ihpdb(190),nwrith,lwrit)) goto 100
               endif

             endif

             ist = 0
            endif

100         continue
           end do
         end do
         if (ist.eq.0) then
            nzpdb = nwrit
            call haszm(.true.)
c            goto 5
            goto 1000
         endif
      endif
      if (debug) call prtzm

c if the z-matrix is to be created using a predefined ordering of
c atoms

      if (isimpl.eq.2) then
         if (opfil(46,'mapfile',7,1,1,1)) then
            iuntmp = iun2
            iun2 = 46
            do i=1,iatoms
               if (getlin(0).eq.1) then
                  ktype = nxtwrd(str,nstr,itype,rtype)
                  if (ktype.eq.2) then
                     if (i.le.numat1) jring(i) = itype
                  else
                     isimpl = 1
                  endif
               else
                  isimpl = 1
               endif
            end do
            iun2 = iuntmp
            close(46)
         else
            isimpl = 1
         endif
      endif

5     if (isimpl.ge.1) goto 910
      
      do i=1,4
         isel(i) = 0
      end do

      if (ilead.eq.0) then
         do i=1,iatoms
            if (.not.ringg(i,iring,nring,.false.,ianz,iaton,iconn,
     &                     lwrit,lring)) goto 2000
            if (nring.ge.5) ilead = i
         end do
      endif
      if (ilead.eq.0) then
c        get other lead, eindstandig iets
         do i=1,4
            isel(i) = 0
         end do
         icmin = 100 
         do i=1,iatoms
            if (ianz(i).eq.6.and.iaton(i).eq.2) then
               icnos = 0
               do j=1,iconn(1,i)
                  jj = abs(iconn(j+1,i))
                  if (ocnos(jj,ianz,iaton)) icnos = icnos + 1
               end do
               if (icnos.lt.icmin) then
                   icmin = icnos
                   icat  = i
               endif
            endif
         end do
         if (icmin.ne.100) ilead = icat
         if (ilead.eq.0) then
            icmin = 100 
            do i=1,iatoms
               if (ocnos(i,ianz,iaton)) then
                  icnos = 0
                  do j=1,iconn(1,i)
                     jj = abs(iconn(j+1,i))
                     if (ocnos(jj,ianz,iaton)) icnos = icnos + 1
                  end do
                  if (icnos.lt.icmin) then
                      icmin = icnos
                      icat  = i
                  endif
               endif
            end do
            if (icmin.ne.100) ilead = icat
         endif
         if (ilead.eq.0) then
            mconn = -1 
            do i=1,iatoms
               if (ianz(i).ne.1.and.iaton(i).eq.2) then
                   if (iconn(1,i).gt.mconn) then
                      mconn = iconn(1,i)
                      ilead = i
                   endif
               endif
            end do
         endif
         if (ilead.eq.0) ilead = 1
      endif

      if (debug) print*,'start lead ',ilead
      
20    continue
      if (ilead.eq.ilold) goto 2000

      ilold = ilead
      if (.not.ringg(ilead,iring,nring,.false.,ianz,iaton,iconn,
     &               lwrit,lring)) goto 2000
      if (nring.ge.5) then
         if (debug) print*,'found ring n=',nring
c
c     ilead is a ring member
c
         if (.not.wring(iring,nring,isel,ianz,iaton,iconn,lwrit,lring)) 
     &      goto 2000
         if (debug) call prtzm
c
c        check for connected rings
c
10       do i=1,iatoms
            if (lring(i).eq.1) then
               if (.not.ringg(i,iring,nring,.false.,ianz,iaton,iconn,
     &                        lwrit,lring)) goto 2000
               if (nring.ge.5) then
                  if (.not.wring(iring,nring,isel,ianz,iaton,iconn,
     &                           lwrit,lring)) goto 2000
                  goto 10
               endif
            endif
         end do
      else
c
c     ilead is not a ring member
c
         if (isel(1).eq.0) then
c
c           this a new lead, so get 3 atoms to which it is connected
c
            idum = 0
            if (parlea(ilead,isel,ispdb)) idum = 1
            if (debug) then
               if (idum.eq.1) then
                  print*,'parse new lead'
               else
                  print*,'failed parsing new lead'
               endif
            endif
         else
            if (debug) print*,'parse extended lead'
            if (.not.pcklin(isel)) goto 2000
         endif
         if (debug) call prtzm
c
c        get atom connected to prev lead and ne isel(2)
c
         ilead = 0
         do j=1,iconn(1,isel(1))
            jj = abs(iconn(j+1,isel(1)))
            if (ianz(jj).ne.1.and.jj.ne.isel(2).and.lwrit(jj).eq.0
     &          .and.iaton(jj).eq.2)   ilead = jj
         end do
         if (ilead.eq.0) goto 1000
         if (debug) print*,'extend old lead ',ilead
         isel(4) = isel(3)
         isel(3) = isel(2)
         isel(2) = isel(1)
         isel(1) = ilead
         goto 20
      endif

1000  continue
      if (debug) print*,'make brand new lead'
      
c
c     get brand new lead, first ring connected then other ?
c
      ilead = 0
c      if (ispdb.eq.1) goto 1010
c
c     Check for atoms connected to rings
c
      do i=1,iatoms
         if (lring(i).ne.0) then
             do j=1,iconn(1,i)
                jj = abs(iconn(j+1,i))
                if (ianz(jj).ne.1.and.lwrit(jj).eq.0
     &          .and.iaton(jj).eq.2) then
                    if (.not.(ispdb.eq.1.and.iresid(i).gt.0)) then
                       ilead = jj
                       if (debug) print*,'ring atom'
                       goto 900
                    endif
                endif
             end do
         endif
      end do
      if (debug) print*,'not a ring atom'
c
c     other connected atoms
c
      if (ilead.eq.0) then
         do i=1,iatoms
            if (lwrit(i).ne.0.and.iaton(i).eq.2) then
                do j=1,iconn(1,i)
                   jj = abs(iconn(j+1,i))
                   if (ianz(jj).ne.1.and.ianz(jj).lt.99
     &             .and.lwrit(jj).eq.0.and.iaton(jj).eq.2) then
                      if (.not.(ispdb.eq.1.and.iresid(i).gt.0)) then
                       ilead = jj
                       if (debug) 
     &                     print*,'connected atom ',ilead,' to ',i,
     &                            ' written nr. ',lwrit(i)
                       goto 900
                      endif
                   endif
                end do
            endif
         end do
      endif
1010  continue
      if (debug) print*,'not a connected atom'
c
c     Check for non bonded atoms
c
      if (ilead.eq.0) then
         nanz = 1
         if (nwrit-nwrith.lt.noth.or.ispdb.eq.1) then
            ispdbt = ispdb
            if (ispdb.eq.1.and.nchain.le.0) ispdbt = 0
            call icrcon(ilead,isel,idisc,ndisc,nanz,ispdbt)
            if (ilead.ne.0.and.ispdbt.eq.1) then
c               ilead = 0
               goto 20
            endif
c         else
c            print*,'nwrit=',nwrit,' nwrith=',nwrith,' noth=',noth
         endif
      endif
      if (debug) print*,'non bonded atom'
900   if (ilead.ne.0) then
         if (debug) 
     &       print*,'found new lead , ilead=',ilead,' ',ianz(ilead)
         do i=1,4
            isel(i) = 0
         end do
         goto 20
      endif

910   continue
c      print*,'910 nwrit=',nwrit,' nwrith=',nwrith,' noth=',noth
c
c     Take care of hydrogens
c
      ilo = -1
      do l=1,iatoms
         il = l
         if (isimpl.ge.2) then
            if (l.le.numat1) il = jring(l)
         endif

         if (isimpl.ge.1.and.il.eq.0) then
            if (ilo.eq.-1) goto 2000

            do ii=1,iatoms
               if (ii.le.numat1) jring(ii) = 0
            end do

            isimpl = 0
c           No more atoms preferred ordering, continue with none
            ilead = 0

            do j=1,iconn(1,isel(1))
               jj = abs(iconn(j+1,isel(1)))
               if (ianz(jj).ne.1.and.jj.ne.isel(2).and.lwrit(jj).eq.0
     &             .and.iaton(jj).eq.2)   ilead = jj
            end do

            if (ilead.eq.0) goto 1000
            if (debug) print*,'extend old lead ',ilead

            isel(4) = isel(3)
            isel(3) = isel(2)
            isel(2) = isel(1)
            isel(1) = ilead

            if (debug) print*,'extend isel',(isel(ii),ii=1,4)
            goto 20
         endif

         if ((ianz(il).eq.1.or.isimpl.ge.1).and.lwrit(il).eq.0
     &        .and.iaton(il).eq.2) then
            ilo = il
            ilead = il
            if (.not.parlea(il,isel,0)) then
               if (ifcon(il,idisc,ndisc,0,ianz,iaton,iconn,lwrit)
     &                   .ne.0) then
                  if (.not.parlea(il,isel,0)) then
                      idum = 0
                  else
                    if (isimpl.ge.1) then
                      if (.not.ringg(il,iring,nring,.false.,
     &                               ianz,iaton,iconn,lwrit,lring)) then
                          goto 2000
                      else
                         lring(il) = 1
                      endif
                    endif
                  endif
               endif
            else
               if (isimpl.ge.1) then
                  if (.not.ringg(il,iring,nring,.false.,
     &                           ianz,iaton,iconn,lwrit,lring)) then
                      goto 2000
                  else
                     lring(il) = 1
                  endif
               endif
            endif
            if (debug) call prtzm
         endif
      end do

2000  continue

c
c     Disconnect Fake connections
c
      do i=1,ndisc,2

         j = idisc(i)
         k = idisc(i+1)

         do m=1,iconn(1,j)
c            if (abs(iconn(1+m,j)).eq.k) then
            if (iconn(1+m,j).eq.k) then
               do n=m+1,iconn(1,j)
                  iconn(n,j) = iconn(n+1,j)
               end do
               iconn(1,j) = iconn(1,j)-1
            endif
         end do

         do m=1,iconn(1,k)
c            if (abs(iconn(1+m,k)).eq.j) then
            if (iconn(1+m,k).eq.j) then
               do n=m+1,iconn(1,k)
                  iconn(n,k) = iconn(n+1,k)
               end do
               iconn(1,k) = iconn(1,k)-1
            endif
         end do

      end do

      if (ispdb.eq.1.and.nzpdb.ne.0) then
         call setzm(nzpdb)
      endif

      if (nwrit.lt.nsel) then
          print*,'Error in automatic cartesian -> z-matrix conversion'
          print*,'No. of atoms written ',nwrit,' not equal to total ',
     &           'No. of selected atoms ',nsel
          if (.true.) then
              do i=1,iatoms
                 if (iaton(i).eq.2.and.lwrit(i).eq.0) then
                    print*,i,' ',ianz(i),(iconn(j+1,i),j=1,iconn(1,i))
                 endif
              end do
          endif
      endif

      if (isimpl.ge.2) isimpl = 0

      return
      end

      subroutine icrcod(icrcon,isel,idisc,ndisc,nanz,ispdb,
     &                        ianz,iaton,iconn,lwrit,coo,
     &                        icalf,ncalf)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      parameter (ncmax=1000)
      logical pline,chkclf
      real bl,rmin
      dimension isel(4),idisc(ncmax)
      dimension ianz(*),iaton(*),iconn(mxcon+1,*),coo(3,*),lwrit(*)
      dimension icalf(6,*)

      icrcon = 0

      do i=1,4
         isel(i) = 0
      end do

      rmin = 1.0e5
      im1 = 0
      im2 = 0
      jtarg = 0
      do i=1,iatoms
         if (ianz(i).ne.nanz.and.lwrit(i).eq.0.and.
     &       iaton(i).eq.2) then
            if (ispdb.eq.1) then
               jtarg = 1
               dmn = 10000000.d0
               do j=1,ncalf
                  if (chkclf(j,icalf).and.
     &               iconn(1,icalf(1,j)).lt.mxcon) then
                     d2 = dist2(coo(1,icalf(1,j)),coo(1,i))
                     if (d2.lt.dmn
     &                  .and.(i.ne.icalf(1,j))
     &                  .and.(i.ne.icalf(2,j))
     &                  .and.(i.ne.icalf(3,j))
     &                  .and.(i.ne.icalf(4,j))) then
                        dmn = d2
                        jtarg = j
                     endif
                  endif
               end do
               isel(1) = i
               isel(2) = icalf(1,jtarg)
               call intcor(intc,bl,isel,2)
               if (intc.eq.0) then
                  call haszm(.false.)
                  return
               endif
               if (bl.lt.rmin) then
                  rmin = bl
                  im1 = i
                  im2 = icalf(1,jtarg)
               endif
               goto 10
            else
               do j=1,iatoms
                  if (ianz(j).ne.nanz.and.ianz(j).lt.99
     &            .and.iconn(1,j).lt.mxcon
     &            .and.lwrit(j).ne.0.and.iaton(j).eq.2.and.i.ne.j) then
                     isel(1) = i
                     isel(2) = j
                     call intcor(intc,bl,isel,2)
                     if (intc.eq.0) then
                        call haszm(.false.)
                        return
                     endif
                     if (bl.lt.rmin) then
                        rmin = bl
                        im1 = i
                        im2 = j
                     endif
                  endif
               end do
            endif
         endif
      end do

10    icrcon = jcrcon(im1,im2,idisc,ndisc,iconn)
      if (ispdb.eq.1.and.jtarg.ne.0) then
         isel(3) = icalf(2,jtarg)
         isel(4) = icalf(3,jtarg)
c         if (.not.pline(isel)) goto 100
      endif

      return
100   print*,'error'
      return
      end

      logical function chkclf(ires,icalf)
      implicit double precision (a-h,o-z)
      dimension icalf(6,*)

      chkclf = .true.
      do i=1,3
         if (icalf(1,ires).eq.icalf(1+i,ires)) chkclf = .false.
      end do
      do i=1,2
         if (icalf(2,ires).eq.icalf(2+i,ires)) chkclf = .false.
      end do
      if (icalf(3,ires).eq.icalf(4,ires)) chkclf = .false.

c      print*,'chkclf ires=',ires,' ',chkclf
c      if (.not.chkclf) print*,(icalf(i,ires),i=1,4)

      return
      end

      integer function ifcon(ilead,idisc,ndisc,nanz,
     &                       ianz,iaton,iconn,lwrit)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (ncmax=1000)
      common /athlp/ iatoms, mxnat
      real bl,rmin
      dimension isel(4),idisc(ncmax)
      dimension ianz(*),iaton(*),iconn(mxcon+1,*),lwrit(*)

      ifcon = 0

      rmin = 1.0e5
      im1 = 0
      im2 = 0
      do j=1,iatoms
         if (ianz(j).ne.nanz.and.ianz(j).lt.99
     &   .and.iconn(1,j).lt.mxcon
     &   .and.lwrit(j).ne.0.and.iaton(j).eq.2.and.ilead.ne.j) then
            isel(1) = ilead
            isel(2) = j
            call intcor(intc,bl,isel,2)
            if (intc.eq.0) then
               call haszm(.false.)
               return
            endif
            if (bl.lt.rmin) then
               rmin = bl
               im1 = ilead
               im2 = j
            endif
         endif
      end do

      ifcon = jcrcon(im1,im2,idisc,ndisc,iconn)

      return
      end

      integer function jcrcon(im1,im2,idisc,ndisc,iconn)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (ncmax=1000)
      dimension idisc(ncmax),iconn(mxcon+1,*)

      jcrcon = 0

      icim1 = 0
      icim2 = 0
      icimt = 0

      if (im1.ne.0.and.im2.ne.0) then
c
c Create fake connection
c

c check if they are not connected

          
          do i=1,iconn(1,im1)
             if (iconn(1+i,im1).eq.im2) icim1 = 1
          end do

          do i=1,iconn(1,im2)
             if (iconn(1+i,im2).eq.im1) icim2 = 1
          end do

          icimt = icim1 + icim2

          if (icimt.eq.2) then
              jcrcon = im1
              return
          endif

          if (iconn(1,im1).lt.mxcon.and.
     &        iconn(1,im2).lt.mxcon) then

             jcrcon =  im1
             iconn(2+iconn(1,im1),im1) = im2
             iconn(1,im1) = iconn(1,im1) + 1
             iconn(2+iconn(1,im2),im2) = im1
             iconn(1,im2) = iconn(1,im2) + 1

             if (ndisc.lt.ncmax-1) then
                 ndisc = ndisc + 1
                 idisc(ndisc) = im1
                 ndisc = ndisc + 1
                 idisc(ndisc) = im2
             else
                 print*,'intzmt: array to hold fake conn. full'
             endif

          else
              print*,'intzmt: couldnt create fake connection'
          endif

      endif

      return
      end

      subroutine getrng(iat,iout,ianz,iconn)
      implicit double precision (a-h,o-z)
      common /athlp/  iatoms, mxnat
      parameter (mxcon=10)
      dimension ianz(*),iconn(mxcon+1,*)
   
      iout = 0
      call ring5(nring,iat,ianz,iconn)

      if (nring.eq.5) iout = 1
   
      return
      end

      subroutine ring5(nring,iat,ianz,iconn)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      real dihed
      dimension isel(4),iring(6),ianz(*),iconn(mxcon+1,*)
      
c oflat finds flat 6,5 membered rings
  
      do i=1,6
         iring(i) = 0
      end do

      nring = 0
      ia = ianz(iat)
      if (ia.eq.6.or.ia.eq.7.or.ia.eq.8.or.ia.eq.16) then
         do j=1,iconn(1,iat)
            jj = abs(iconn(j+1,iat))
            ja = ianz(jj)
            if (ja.eq.6.or.ja.eq.7.or.ja.eq.8.or.ja.eq.16) then
               do k=1,iconn(1,jj)
                  kk = abs(iconn(k+1,jj))
                  ka = ianz(kk)
                  if ((ka.eq.6.or.ka.eq.7.or.ka.eq.8.or.ka.eq.16)
     &                 .and.kk.ne.iat) then
                     do l=1,iconn(1,kk)
                        ll = abs(iconn(l+1,kk))
                        la = ianz(ll)
                        if ((la.eq.6.or.la.eq.7.or.la.eq.8.or.la.eq.16)
     &                      .and.ll.ne.iat.and.ll.ne.jj) then
                           isel(1) = iat
                           isel(2) = jj
                           isel(3) = kk
                           isel(4) = ll
                           call intcor(intc,dihed,isel,4)
                           if (intc.eq.0) then
                              nring = 0
                              return
                           endif

                           if (abs(dihed).lt.5.0e0) then
                              do m=1,iconn(1,iat)
                                 mm = abs(iconn(m+1,iat))
                                 ma = ianz(mm)
                                 if ((ma.eq.6.or.ma.eq.7.or.ma.eq.8.or.
     &                                ma.eq.16).and.mm.ne.jj) then
                                    do n=1,iconn(1,ll)
                                       nn = abs(iconn(n+1,ll))
                                       na = ianz(nn)
                                       if ((na.eq.6.or.na.eq.7.or.
     &                                      na.eq.8.or.na.eq.16)
     &                                      .and.nn.ne.kk) then
                                          if (mm.eq.nn) then
                                             nring = 5
                                             iring(1) = iat
                                             iring(2) = jj
                                             iring(3) = kk
                                             iring(4) = ll
                                             iring(5) = mm
                                             call intcor(intc,dihed,
     &                                          iring(2),4)
                                             if (intc.eq.0) then
                                                nring = 0
                                                return
                                             endif

                                             if (abs(dihed).lt.5.0e0)
     &                                       then
                                                return
                                             endif

                                          endif
                                       endif
                                    end do
                                 endif
                              end do
                           endif
                        endif
                     end do
                  endif
               end do
            endif
         end do
      endif

      return
      end

      subroutine ringd(i,iring,nring,iflat,ianz,iaton,iconn,
     &                 lwrit,lring,iret)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      logical oflat,ringg
      dimension iring(6),ianz(*),iaton(*),iconn(mxcon+1,*)
      dimension lwrit(*),lring(*)

      if (iflat.gt.0) then
          oflat = .true.
      else 
          oflat = .false.
      endif

      iret = 0

      if (ringg(i,iring,nring,oflat,ianz,iaton,iconn,lwrit,lring)) then
         iret = 1
      endif

      return
      end

      logical function ringg(i,iring,nring,oflat,ianz,iaton,iconn,
     &                       lwrit,lring)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      logical ocnos,oring,oflat
      real dihed
      dimension isel(4),iring(6),ianz(*),iaton(*),iconn(mxcon+1,*)
      dimension lwrit(*),lring(*)
      
c oflat finds flat 6,5 membered rings
  
      ringg = .true.

      if (ocnos(i,ianz,iaton)) then
         do j=1,iconn(1,i)
            jj = abs(iconn(j+1,i))
            if (ocnos(jj,ianz,iaton)) then
               do k=1,iconn(1,jj)
                  kk = abs(iconn(k+1,jj))
                  if (ocnos(kk,ianz,iaton).and.kk.ne.i) then
                     do l=1,iconn(1,kk)
                        ll = abs(iconn(l+1,kk))
                        if (ocnos(ll,ianz,iaton).and.ll.ne.jj
     &                      .and.ll.ne.i) then
                           isel(1) = i
                           isel(2) = jj
                           isel(3) = kk
                           isel(4) = ll
                           call intcor(intc,dihed,isel,4)
                           if (intc.eq.0) then
                              call haszm(.false.)
                              ringg = .false.
                              return
                           endif
                           if ((.not.oflat.and.abs(dihed).lt.90.0e0)
     &                        .or.(oflat.and.abs(dihed).lt.5.0e0)) then
                              do m=1,iconn(1,i)
                                 mm = abs(iconn(m+1,i))
                                 if (ocnos(mm,ianz,iaton).and.mm.ne.jj) 
     &                           then
                                    do n=1,iconn(1,ll)
                                       nn = abs(iconn(n+1,ll))
                                       if (ocnos(nn,ianz,iaton).and.
     &                                           nn.ne.kk) then
                                          if (mm.eq.nn) then
                                             nring = 5
                                             iring(1) = i
                                             iring(2) = jj
                                             iring(3) = kk
                                             iring(4) = ll
                                             iring(5) = mm
                                             if (oflat) then
                                                call intcor(intc,dihed,
     &                                          iring(2),4)
                                                if (intc.eq.0) then
                                                call haszm(.false.)
                                                ringg = .false.
                                                return
                                                endif
                                             endif
                                             if (.not.oflat.or.(oflat
     &                                       .and.abs(dihed).lt.5.0e0))
     &                                       then
                                                if (oring(nring,iring,
     &                                              lwrit,lring)) then
                                                 return
                                                endif
                                             endif
                                          else
                                          do ii=1,iconn(1,nn)
                                             ij=abs(iconn(ii+1,nn))
                                             if (ocnos(ij,ianz,iaton)
     &                                           .and.ij.ne.ll)
     &                                       then
                                               if (ij.eq.mm.and.i.ne.nn) 
     &                                         then
                                                   nring = 6
                                                   iring(1) = i
                                                   iring(2) = jj
                                                   iring(3) = kk
                                                   iring(4) = ll
                                                   iring(5) = nn
                                                   iring(6) = ij

                                             if (oflat) then
                                                call intcor(intc,dihed,
     &                                                      iring(3),4)
                                                if (intc.eq.0) then
                                                    call haszm(.false.)
                                                    ringg = .false.
                                                    return
                                                endif
                                             endif

                                             if (.not.oflat.or.(oflat
     &                                       .and.abs(dihed).lt.5.0e0))
     &                                       then
                                                if (oring(nring,iring,
     &                                              lwrit,lring)) then
                                                  return
                                                endif
                                             endif
                                               endif
                                             endif
                                          end do
                                          endif
                                       endif
                                    end do
                                 endif
                              end do
                           endif
                        endif
                     end do
                  endif
               end do
            endif
         end do
      endif

      nring = 0

      return
      end

      logical function wring(iring,nring,isel,ianz,iaton,iconn,
     &                       lwrit,lring)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /zmatst/ nwrit,nvar,isimpl
      logical onorng,pline
      real angle
      dimension isel(4),isl(4),iring(6),kring(18)
      dimension ianz(*),iaton(*),iconn(mxcon+1,*),lwrit(*),lring(*)


c     first reorder

      wring = .true.
      tol = 1.0d-10

      do i=1,nring
         kring(6+i) = iring(i)
         kring(6-nring+i) = iring(i)
         kring(6+nring+i) = iring(i)
      end do

      iseed = 1
      isens = 1

      do i=1,nring

         if (lring(kring(6+i)).eq.1.and.
     &   (lring(kring(5+i)).eq.0.or.lring(kring(7+i)).eq.0)) then

            if (lring(kring(7+i)).eq.0) isens = -1
            iseed = i
            goto 10
         endif

      end do

10    do i=1,nring
         iring(i) = kring(6+iseed+isens*(i-1))
      end do

      do i=1,nring
         kring(3+i) = iring(i)
      end do

      itot = 0

      do i=1,nring
        itot = itot + lring(iring(i)) 
      end do

      kring(3) = 0
      kring(2) = 0
      kring(1) = 0
c
c Search for atoms this ring is connected to
c
c      if (itot.le.2.and.isel(1).eq.0) then
      if (itot.le.2.and.(isel(1).ne.iring(1))) then
        iat = iring(1)
50      do i=1,iconn(1,iat)
           ii = abs(iconn(i+1,iat))
           if (ianz(ii).ne.1.and.onorng(iring,nring,ii,lwrit).and.
     &         iaton(ii).eq.2) then
              kring(3) = ii
              do j=1,iconn(1,kring(3))
                 jj = abs(iconn(j+1,kring(3)))
                 if (ianz(jj).ne.1.and.onorng(iring,nring,jj,lwrit)
     &              .and.iaton(jj).eq.2) then
                    if (kring(2).ne.0) kring(1) = kring(2)
                    kring(2) = jj
                    isl(1) = iat
                    isl(2) = kring(3)
                    isl(3) = kring(2)
                    call tomold(angle,isl,3)
                    if (.not.(abs(angle).lt.tol.or.
     &                   abs(angle).gt.180.e0-tol)) then
                       do k=1,iconn(1,kring(2))
                          kk = abs(iconn(k+1,kring(2)))
                          if (ianz(kk).ne.1.and.
     &                        onorng(iring,nring,kk,lwrit).and.
     &                        kk.ne.kring(3).and.iaton(kk).eq.2) then
                             kring(1) = kk
                             goto 100
                          endif
                       end do
                    endif
                 endif
              end do
           endif
        end do
100     continue
        if ((kring(3).eq.0.and.nwrit.gt.2).and.iat.ne.iring(2)) then
           iat = iring(2)
           goto 50
        endif
      else
         kring(3) = isel(2)
         kring(2) = isel(3)
         kring(1) = isel(4)
      endif

      do i=1,nring
         j = iring(i)
         if (lwrit(j).eq.0) then
            lring(j) = 1
            isel(1) = kring(3+i)
            isel(2) = kring(2+i)
            isel(3) = kring(1+i)
            isel(4) = kring(i)
            if (.not.pline(isel)) then
                isel(1) = 0
                wring = .false.
                return
            endif
         endif
      end do

      isel(1) = 0
      return
      end

      logical function ocnos(inum,ianz,iaton)
      implicit double precision (a-h,o-z)
      dimension ianz(*),iaton(*)
      
      ocnos = .false.
      iatnr = ianz(inum)
      if ((iatnr.eq.6.or.iatnr.eq.7.or.iatnr.eq.8.or.iatnr.eq.16)
     &    .and.iaton(inum).eq.2) ocnos = .true.

      return
      end

      logical function oring(nring,iring,lwrit,lring)
      implicit double precision (a-h,o-z)
      dimension iring(6),lwrit(*),lring(*)

      oring = .false.
      itot = 0
      jtot = 0
      do i=1,nring
        itot = itot + lring(iring(i)) 
        if (lwrit(iring(i)).ne.0) jtot = jtot + 1
      end do
      if (itot.lt.nring.and.jtot.lt.nring) oring = .true.

      return
      end

      logical function onorng(iring,nring,isub,lwrit)
      implicit double precision (a-h,o-z)
      dimension iring(6),lwrit(*)

      onorng = .true.
      do i=1,nring
        if (isub.eq.iring(i)) onorng = .false.
      end do
      if (lwrit(isub).eq.0) onorng = .false.

      return
      end

      subroutine plinz(isel,istat,
     &              bl,alph,bet,ibl,ialph,ibet,imap,ianzz,iz,
     &              lwrit,ianz)

c this is really logical function pline

      implicit double precision (a-h,o-z)
      common /zmfrst/ ihaszm, nz, mxzat
      logical lineok
      real bond,angle,dihed
      common /zmatst/ nwrit,nvar,isimpl
      dimension isel(4),ian(4)
      dimension bl(*),alph(*),bet(*),ibl(*),ialph(*),ibet(*),
     &          imap(*),ianzz(*),iz(4,*),lwrit(*),ianz(*)

      toang = 0.52917706d0

      istat = 1
      lineok = .true.

      nchk = nwrit + 1
      if (nchk.gt.4) nchk = 4
      do i=1,nchk
         lineok = (lineok.and.isel(i).ne.0)
      end do
      if (.not.lineok) goto 100

      nwrit = nwrit + 1
      lwrit(isel(1)) = nwrit
      
      do i=1,nchk
         ian(i) = ianz(isel(i))
      end do

      if (nwrit.ge.1) then
          ianzz(nwrit) = ian(1)
          imap(nwrit) = isel(1)
          iz(1,nwrit) = 0
          iz(2,nwrit) = 0
          iz(3,nwrit) = 0
          iz(4,nwrit) = 0
      endif
      if (nwrit.ge.2) then
          call intcor(intc,bond,isel,2)
          if (intc.eq.0) goto 90
          bond = bond * toang
          bl(nwrit) = bond
          ibl(nwrit) = 1
          iz(1,nwrit) = lwrit(isel(2))
      endif
      if (nwrit.ge.3) then
          call intcor(intc,angle,isel,3)
          if (intc.eq.0) goto 90
          alph(nwrit) = angle
          ialph(nwrit) = 1
          iz(2,nwrit) = lwrit(isel(3))
      endif
      if (nwrit.gt.3) then
          call intcor(intc,dihed,isel,4)
          if (intc.eq.0) goto 90
          bet(nwrit) = dihed
          ibet(nwrit) = 1
          iz(3,nwrit) = lwrit(isel(4))
          iz(4,nwrit) = 0
      endif

      nz = nwrit

      return

90    nwrit = nwrit - 1
      lwrit(isel(1)) = 0

100   istat = 0
      call haszm(.false.)

      return
      end

      subroutine dumliz(isel,blv,alphv,betv,
     &                  bl,alph,bet,ibl,ialph,ibet,imap,ianzz,iz,
     &                  lwrit,ianz)

c this is really dumlin

      implicit double precision (a-h,o-z)
      common /zmfrst/ ihaszm, nz, mxzat
      common /zmatst/ nwrit,nvar,isimpl
      dimension isel(4)
      dimension bl(*),alph(*),bet(*),ibl(*),ialph(*),ibet(*),
     &          imap(*),ianzz(*),iz(4,*),lwrit(*),ianz(*)

c
c     Bond length expected in Angstroms
c
      nwrit = nwrit + 1
      lwrit(isel(1)) = nwrit

      ianzz(nwrit) = ianz(isel(1))

      do i=1,3
         iz(i,nwrit) = lwrit(isel(i+1))
      end do

      iz(4,nwrit) = 0

      bl(nwrit) = blv
      ibl(nwrit) = 1
      alph(nwrit) = alphv
      ialph(nwrit) = 1
      bet(nwrit) = betv
      ibet(nwrit) = 1

      imap(nwrit) = isel(1)

      nz = nwrit

      return
      end

      subroutine calcd(ical,isel,nx,ianz,iaton,iconn,coo)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      common /zmatst/ nwrit,nvar,isimpl
      dimension v32(3),v43(3),c1(3),c2(3),isel(4)
      dimension coo(3,*),ianz(*),iaton(*),iconn(mxcon+1,*)

c
c     Calculate coordinates dummy atom X, and add
c
      ical = 1
      
      tol = 1.0d-10
      toang = 0.52917706d0

      do i=1,3
          v32(i) = coo(i,isel(3)) - coo(i,isel(2))
          if (isel(4).ne.0) then
             v43(i) = coo(i,isel(4)) - coo(i,isel(3))
          else
             v43(i) = 0.0d0
          endif
      end do
      if (isel(4).eq.0) then
          v43(1) = 1.0d0
          call impsc(v32,v43,cosb)
          if (dabs(dabs(cosb)-1.0d0).lt.tol) then
              v43(1) = 0.0d0
              v43(2) = 1.0d0
          endif
          
      endif

      do i=1,4
         if ((isel(i).lt.1.and.nwrit.gt.2).or.isel(i).gt.iatoms) 
     &   goto 100
         do j=i+1,4
             if (isel(i).eq.isel(j)) goto 100
         end do
      end do

      call crprod(v32,v43,c1)
      call crprod(c1,v32,c2)
      call vsc1(c2,1.0d0/toang,tol)
      
      iatoms = iatoms + 1
      nx = iatoms

      do i=1,3
         coo(i,nx) = coo(i,isel(2)) + c2(i)
      end do

      iconn(1,nx) = 1
      iconn(2,nx) = isel(2)
      ianz(nx) = 99
      iaton(nx) = 2

      if (iconn(1,isel(2)).lt.mxcon) then
         iconn(1,isel(2)) = iconn(1,isel(2)) + 1
         iconn(iconn(1,isel(2))+1,isel(2)) = nx
      endif


      return

100   ical = 0
      return
      end

      subroutine atadd(iat1,iat2,iat3,ian,bl,alpha,dih,iret,
     &                 ianz,iaton,iatclr,iconn,iresid,ityp,ipdbt,
     &                 ncalf,icalf,coo)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      logical addat
      integer*2 ityp,ipdbt
      dimension coo(3,*),ianz(*),iaton(*),iatclr(*),iresid(*),
     &          iconn(mxcon+1,*),ityp(*),ipdbt(*),icalf(6,*)

      if (.not.addat(iat1,iat2,iat3,ian,bl,alpha,dih,iret,0,
     &               ianz,iaton,iatclr,iconn,iresid,ityp,ipdbt,
     &               ncalf,icalf,coo)) then
          print*,'error adding atom'
      endif
     
      return
      end

      logical function addat(iat1,iat2,iat3,ian,bl,alpha,dih,iret,ichk,
     &                       ianz,iaton,iatclr,iconn,iresid,ityp,ipdbt,
     &                       ncalf,icalf,coo)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      integer*2 ityp,ipdbt
      common /athlp/ iatoms, mxnat
      dimension v21(3),v32(3),c1(3),c2(3),c3(3),c4(3),ctmp(3)
      dimension coo(3,*),ianz(*),iaton(*),iatclr(*),iresid(*),
     &          iconn(mxcon+1,*),ityp(*),ipdbt(*),icalf(6,*)

c
c     Add atom given bondlength, bondangle and dihedral with
c     respect to three atoms
c
      iret = 0
      addat = .true.
      
      if (iat1.le.0.or.iat1.gt.mxnat) goto 100
      if (iat2.le.0.or.iat2.gt.mxnat) goto 100
      if (iat3.le.0.or.iat3.gt.mxnat) goto 100

      tol = 1.0d-10
      toang = 0.52917706d0
      todeg = 45.0d0 / datan(1.0d0)
      sa = dsin(alpha/todeg)
      ca = dcos(alpha/todeg)
      sd = dsin(dih/todeg)
      cd = dcos(dih/todeg)

      do i=1,3
          v21(i) = coo(i,iat2) - coo(i,iat1)
          v32(i) = coo(i,iat3) - coo(i,iat2)
      end do
      call vsc1(v32,1.0d0,tol)

c     check for three atoms on a line

      call impsc(v21,v32,cosb)
      if (dabs(cosb).lt.tol) goto 100

      call crprod(v21,v32,c1)
      call vsc1(c1,1.0d0,tol)
      call crprod(c1,v32,c2)
      call vsc1(c2,1.0d0,tol)

      do i=1,3
          c3(i) = cd*c2(i) + sd*c1(i)
      end do
      
      do i=1,3
          c4(i) = -ca*v32(i) + sa*c3(i)
      end do

      call vsc1(c4,bl/toang,tol)

      if (iatoms.ge.mxnat) goto 100

      do i=1,3
         ctmp(i) = coo(i,iat3) + c4(i)
      end do

      iflg = 0
      do i=1,iatoms
          d2 = dist2(ctmp,coo(1,i))
          if (d2*toang*toang.lt.0.7d0) then
             iflg = 1
             d2t = d2
             i2t = i
          endif
      end do

c      if (iflg.eq.1) then
c          print*,'i2t=',i2t,'d2=',d2t*toang*toang,
c     &  ' iatoms=',iatoms,' bl=',bl
c          print*,'ctmp ',(ctmp(i)*toang,i=1,3)
c          print*,'ci2t ',(coo(i,i2t)*toang,i=1,3)
c      endif
      if (iflg.eq.1.and.ichk.eq.1) goto 100

      natoms = 0
      if (iresid(iat3).lt.-3) then
         do i=1,iatoms
            if (iresid(i).eq.iresid(iat3)) then
               natoms = i + 1
            endif
         end do
      endif

      if (natoms.eq.0) then
          natoms = iatoms + 1
      else


         do i=1,natoms-1
            do j=1,iconn(1,i)
               it = iconn(1+j,i)
               iat = iabs(it)
               if (it.ge.0) then
                  if (it.ge.natoms) then
                     iconn(1+j,i) = it + 1
                  endif
               else
                  if (iat.ge.natoms) then
                     iconn(1+j,i) = it - 1
                  else
                     iconn(1+j,i) = it
                  endif
               endif
            end do
         end do

         do i=iatoms,natoms,-1
            iconn(1,i+1) = iconn(1,i)
            do j=1,iconn(1,i)
               it = iconn(1+j,i)
               iat = iabs(it)
               if (it.ge.0) then
                  if (it.ge.natoms) then
                     iconn(1+j,i+1) = it + 1
                  else
                     iconn(1+j,i+1) = it
                  endif
               else
                  if (iat.ge.natoms) then
                     iconn(1+j,i+1) = it - 1
                  else
                     iconn(1+j,i+1) = it
                  endif
               endif
            end do

            ianz(i+1)    = ianz(i)
            iaton(i+1)   = iaton(i)
            iatclr(i+1)  = iatclr(i)
            iresid(i+1)  = iresid(i)
            ityp(i+1)    = ityp(i)
            ipdbt(i+1)   = ipdbt(i)

            do j=1,3
               coo(j,i+1) = coo(j,i)
            end do

            call sftlab(i)

         end do

         do i=1,ncalf
            if (icalf(1,i).ge.natoms) icalf(1,i) = icalf(1,i) + 1
            if (icalf(4,i).ge.natoms) icalf(4,i) = icalf(4,i) + 1
         end do

      endif

      iatoms = iatoms + 1
      iret = natoms

      do i=1,3
         coo(i,natoms) = coo(i,iat3) + c4(i)
      end do

      iconn(1,natoms) = 1
      iconn(2,natoms) = iat3
      ianz(natoms)    = ian
      iaton(natoms)   = 1
      iatclr(natoms)  = iatclr(iat3)
      iresid(natoms)  = iresid(iat3)
      ityp(natoms)    = 0
      ipdbt(natoms)   = 0

      if (iconn(1,iat3).lt.mxcon) then
         iconn(1,iat3) = iconn(1,iat3) + 1
         iconn(iconn(1,iat3)+1,iat3) = natoms
      endif

      return

100   addat = .false.
c      print*,'ierr=',ierr
      return
      end

      subroutine fliph(iat1,iat2,iat3,istat,bl,alpha,dih,cret,coo)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      dimension v21(3),v32(3),c1(3),c2(3),c3(3),c4(3)
      dimension coo(3,*),cret(*)

      istat = 0
      
      if (iat1.le.0.or.iat1.gt.mxnat) goto 100
      if (iat2.le.0.or.iat2.gt.mxnat) goto 100
      if (iat3.le.0.or.iat3.gt.mxnat) goto 100

      tol = 1.0d-10
      toang = 0.52917706d0
      todeg = 45.0d0 / datan(1.0d0)
      sa = dsin(alpha/todeg)
      ca = dcos(alpha/todeg)
      sd = dsin(dih/todeg)
      cd = dcos(dih/todeg)

      do i=1,3
          v21(i) = coo(i,iat2) - coo(i,iat1)
          v32(i) = coo(i,iat3) - coo(i,iat2)
      end do
      call vsc1(v32,1.0d0,tol)

c     check for three atoms on a line

      call impsc(v21,v32,cosb)
      if (dabs(cosb).lt.tol) goto 100

      call crprod(v21,v32,c1)
      call vsc1(c1,1.0d0,tol)
      call crprod(c1,v32,c2)
      call vsc1(c2,1.0d0,tol)

      do i=1,3
          c3(i) = cd*c2(i) + sd*c1(i)
      end do
      
      do i=1,3
          c4(i) = -ca*v32(i) + sa*c3(i)
      end do

      call vsc1(c4,bl/toang,tol)

      do i=1,3
         cret(i) = coo(i,iat3) + c4(i)
      end do

      return

100   istat = 1
      return
      end

      logical function chkat(iat1,iat2,iat3,bl,alpha,dih,vec,coo)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      dimension v21(3),v32(3),c1(3),c2(3),c3(3),c4(3)
      dimension vec(3),coo(3,*)

c
c     create coordinates of hypothetical added atom given 
c     bondlength, bondangle and dihedral with respect to three atoms
c
      chkat = .true.
      
      if (iat1.le.0.or.iat1.gt.mxnat) goto 100
      if (iat2.le.0.or.iat2.gt.mxnat) goto 100
      if (iat3.le.0.or.iat3.gt.mxnat) goto 100

      tol = 1.0d-10
      toang = 0.52917706d0
      todeg = 45.0d0 / datan(1.0d0)
      sa = dsin(alpha/todeg)
      ca = dcos(alpha/todeg)
      sd = dsin(dih/todeg)
      cd = dcos(dih/todeg)

      do i=1,3
          v21(i) = coo(i,iat2) - coo(i,iat1)
          v32(i) = coo(i,iat3) - coo(i,iat2)
      end do
      call vsc1(v32,1.0d0,tol)

c     check for three atoms on a line

      call impsc(v21,v32,cosb)
      if (dabs(cosb).lt.tol) goto 100

      call crprod(v21,v32,c1)
      call vsc1(c1,1.0d0,tol)
      call crprod(c1,v32,c2)
      call vsc1(c2,1.0d0,tol)

      do i=1,3
          c3(i) = cd*c2(i) + sd*c1(i)
      end do
      
      do i=1,3
          c4(i) = -ca*v32(i) + sa*c3(i)
      end do
      call vsc1(c4,bl/toang,tol)

      if (iatoms.ge.mxnat) goto 100

      do i=1,3
         vec(i) = coo(i,iat3) + c4(i)
      end do

      return

100   chkat = .false.
      return
      end

      subroutine wlinz(iun,iopt,igamb,
     &                 bl,alph,bet,ibl,ialph,ibet,imap,ianzz,iz,
     &                 iconn,ianz,ityp,qat)

c this is really wline

      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      parameter (mxel=100)
      character*2 elemnt
      common /elem/  elemnt(mxel)
      common /zmfrst/ ihaszm, nz, mxzat

      character*2 tolowf
      character*8  bstr,tbstr
      character*10 astr,tastr
      character*7  dstr,tdstr
      common /zmatst/ nwrit,nvar,isimpl
      integer*2 ityp
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      logical nolink

CNF   Common to get the type of the atoms
CNF   0 < H < 10000 ; 10000 <= M < 20000 ; L >= 20000
      integer iqmmm
      character*3 oniom(6)
      character*3 link
      character*14 amblab

      dimension isel(4),ian(4)
      dimension bl(*),alph(*),bet(*),ibl(*),ialph(*),ibet(*),
     &          imap(*),ianzz(*),iz(4,*),
     &          ianz(*),ityp(*),iconn(mxcon+1,*), qat(*)

      data oniom /'H  ','M  ','L  ','H H','M H','L H'/

      toang = 0.52917706d0
      zero = 0.0d0
      i0 = 0
      i1 = 1

      nwrit = nwrit + 1

      isel(1) = nwrit

      do i=1,3
         isel(i+1) = iz(i,nwrit)
      end do

      do i=1,4
         ian(i) = 0
         if (isel(i).gt.0.and.isel(i).le.mxnat) then
            ia = ianzz(isel(i))
            if (ia.gt.0.and.ia.le.mxel) ian(i) = ia
         endif
      end do

CNF   iqmmm (H: 0; M: 1; L: 2)

      iqmmm = ityp(isel(1))/10000
      link = oniom(iqmmm+1)
      do i = 2,iconn(1,isel(1))+1
          j = iconn(i,isel(1))
          ij = ityp(j)/10000
          if (j.ne.0 .and. ij.lt.iqmmm) link = oniom(iqmmm+4)
      end do
      call wrlab(isel(1),amblab,ianz,ityp,qat)

      iop1 = ibl(nwrit)
      iop2 = ialph(nwrit)
      iop3 = ibet(nwrit)
      j1 = iop1
      j2 = iop2
      j3 = iop3
      if (abs(j1).gt.1) j1 = 0
      if (abs(j2).gt.1) j2 = 0
      if (abs(j3).gt.1) j3 = 0

      call mkvar(nwrit,ianzz,iz,bstr,astr,dstr)
      if (nwrit.ge.2) then
          if (iop1.eq.0.and.nolink(nwrit,ibl,nz)) then
             bstr = '        '
             write(bstr,'(f8.5)') bl(nwrit)
          elseif (abs(iop1).gt.1) then
             tbstr = '        '
             call mkvar(iabs(ibl(nwrit)),ianzz,iz,tbstr,tastr,tdstr)
             bstr = tbstr
          else
             nvar = nvar + 1
          endif
      endif

      if (nwrit.ge.3) then
          if (iop2.eq.0.and.nolink(nwrit,ialph,nz)) then
             astr = '          '
             write(astr,'(f8.3)') alph(nwrit)
          elseif (abs(iop2).gt.1) then
             tastr = '        '
             call mkvar(iabs(ialph(nwrit)),ianzz,iz,tbstr,tastr,tdstr)
             astr = tastr
          else
             nvar = nvar + 1
          endif
      endif

      if (nwrit.ge.4) then
          if (iop3.eq.0.and.nolink(nwrit,ibet,nz)) then
             dstr = 'dih    '
             if (bet(nwrit).lt.0.0d0) then
                write(dstr,'(f7.2)') bet(nwrit)
             else
                write(dstr,'(f7.3)') bet(nwrit)
             endif
          elseif (abs(iop3).gt.1) then
             tdstr = '        '
             call mkvar(iabs(ibet(nwrit)),ianzz,iz,tbstr,tastr,tdstr)
             if (iop3.lt.0) then
                dstr(1:1) = '-'
                dstr(2:7) = tdstr(1:6)
             else
                dstr = tdstr
             endif
          else
             nvar = nvar + 1
          endif
      endif


      if (nwrit.eq.1) then
          if (isel(1).eq.0) print*,'isel(1)=0'
          if (iopt.le.2) then
             if (iopt.eq.2 .and. iqmmm.ne.0)
     &        print*,'The first atom must be at the H level only !'
             if (igamb.eq.1) then
                write(iun,'(a14)')amblab
             else
                write(iun,'(a2)')tolowf(elemnt(ian(1)))
             endif
          else
             write(iun,3000) elemnt(ian(1)),zero,
     &              i0,zero,i0,zero,i0,i0,i0,i0
          endif
      elseif (nwrit.eq.2) then
          if (isel(1).eq.0) print*,'isel(1)=0'
          if (isel(2).eq.0) print*,'isel(2)=0'
          if (iopt.le.2) then

             if (iopt.eq.2 .and. iqmmm.ne.0) then
2011            format(a2,1x,i3,1x,a8,1x,a3)
20110           format(a14,1x,i3,1x,a8,1x,a3)
                if (igamb.eq.1) then
                   write(iun,20110) amblab,isel(2), bstr(1:8), link
                else
                   write(iun,2011) tolowf(elemnt(ian(1))),
     &                             isel(2), bstr(1:8), link
                endif
             else
2001            format(a2,1x,i3,1x,a8)
20010           format(a14,1x,i3,1x,a8)
                if (igamb.eq.1) then
                   write(iun,20010) amblab,isel(2), bstr(1:8)
                else
                   write(iun,2001)tolowf(elemnt(ian(1))),
     &                            isel(2), bstr(1:8)
                endif
             endif

             if (iop1.eq.1) 
     &           write(47,'(a8,3x,f9.6)') bstr(1:8),bl(nwrit)
             if (iop1.eq.0.and..not.nolink(nwrit,ibl,nz)) 
     &           write(50,'(a8,3x,f9.6)') bstr(1:8),bl(nwrit)
          else
             write(iun,3000) elemnt(ian(1)),bl(nwrit),
     &              j1,zero,i0,zero,i0,isel(2),i0,i0
             if (abs(iop1).gt.1) 
     &           write(47,'(i4,a3,i4)') abs(iop1),',1,',nwrit
          endif
      elseif (nwrit.eq.3) then
          if (isel(1).eq.0) print*,'isel(1)=0'
          if (isel(2).eq.0) print*,'isel(2)=0'
          if (isel(3).eq.0) print*,'isel(3)=0'
         
          if (iopt.le.2) then

             if(iopt.eq.2 .and. iqmmm.ne.0) then
2012            format(a2,1x,i4,1x,a8,1x,i4,1x,a10,1x,a3)
20120           format(a14,1x,i4,1x,a8,1x,i4,1x,a10,1x,a3)
                if (igamb.eq.1) then
                   write(iun,20120) amblab,
     &               isel(2),bstr(1:8),isel(3),astr(1:10),link
                else
                   write(iun,2012) tolowf(elemnt(ian(1))),
     &               isel(2),bstr(1:8),isel(3),astr(1:10),link
                endif
             else
2002            format(a2,1x,i4,1x,a8,1x,i4,1x,a10)
20020           format(a14,1x,i4,1x,a8,1x,i4,1x,a10)
                if (igamb.eq.1) then
                   write(iun,20020) amblab,
     &               isel(2),bstr(1:8),isel(3),astr(1:10)
                else
                   write(iun,2002) tolowf(elemnt(ian(1))),
     &               isel(2),bstr(1:8),isel(3),astr(1:10)
                endif
             endif

             if (iop1.eq.1) 
     &           write(47,'(a8,3x,f9.6)') bstr(1:8),bl(nwrit)
             if (iop1.eq.0.and..not.nolink(nwrit,ibl,nz)) 
     &           write(50,'(a8,3x,f9.6)') bstr(1:8),bl(nwrit)
             if (iop2.eq.1) 
     &           write(47,'(a10,1x,f8.3)') astr(1:10),alph(nwrit)
             if (iop2.eq.0.and..not.nolink(nwrit,ialph,nz)) 
     &           write(50,'(a10,1x,f8.3)') astr(1:10),alph(nwrit)
          else
             write(iun,3000) elemnt(ian(1)),bl(nwrit),
     &               j1,alph(nwrit),j2,zero,i0,isel(2),isel(3),i0
             if (abs(iop1).gt.1) 
     &           write(47,'(i4,a3,i4)') abs(iop1),',1,',nwrit
             if (abs(iop2).gt.1) 
     &           write(47,'(i4,a3,i4)') abs(iop2),',2,',nwrit
          endif
      else
          if (isel(1).eq.0) print*,'isel(1)=0'
          if (isel(2).eq.0) print*,'isel(2)=0'
          if (isel(3).eq.0) print*,'isel(3)=0'
          if (isel(4).eq.0) print*,'isel(4)=0'
          if (iopt.le.2) then

             if (iopt.eq.2.and.iqmmm.ne.0) then
2013            format(a2,1x,i4,1x,a8,1x,i4,1x,a10,1x,i4,1x,a7,i2,1x,a3)
20130          format(a14,1x,i4,1x,a8,1x,i4,1x,a10,1x,i4,1x,a7,i2,1x,a3)
                if (igamb.eq.1) then
                   write(iun,20130) amblab,
     &              isel(2),bstr(1:8),isel(3),
     &              astr(1:10),isel(4),dstr(1:7),
     &              i0,link
                else
                   write(iun,2013) tolowf(elemnt(ian(1))),
     &              isel(2),bstr(1:8),isel(3),
     &              astr(1:10),isel(4),dstr(1:7),
     &              i0,link
                endif
             else
2003            format(a2,1x,i4,1x,a8,1x,i4,1x,a10,1x,i4,1x,a7)
20030           format(a14,1x,i4,1x,a8,1x,i4,1x,a10,1x,i4,1x,a7)
                if (igamb.eq.1) then
                   write(iun,20030) amblab,
     &              isel(2),bstr(1:8),isel(3),
     &              astr(1:10),isel(4),dstr(1:7)
                else
                   write(iun,2003) tolowf(elemnt(ian(1))),
     &              isel(2),bstr(1:8),isel(3),
     &              astr(1:10),isel(4),dstr(1:7)
                endif
             endif

             if (iop1.eq.1) 
     &           write(47,'(a8,3x,f9.6)') bstr(1:8),bl(nwrit)
             if (iop2.eq.1) 
     &           write(47,'(a10,1x,f8.3)') astr(1:10),alph(nwrit)
             if (iop3.eq.1) 
     &           write(47,'(a7,4x,f8.3)') dstr(1:7),bet(nwrit)
             if (iop1.eq.0.and..not.nolink(nwrit,ibl,nz)) 
     &           write(50,'(a8,3x,f9.6)') bstr(1:8),bl(nwrit)
             if (iop2.eq.0.and..not.nolink(nwrit,ialph,nz)) 
     &           write(50,'(a10,1x,f8.3)') astr(1:10),alph(nwrit)
             if (iop3.eq.0.and..not.nolink(nwrit,ibet,nz)) 
     &           write(50,'(a7,4x,f8.3)') dstr(1:7),bet(nwrit)
          else
             write(iun,3000) elemnt(ian(1)),bl(nwrit),
     &          j1,alph(nwrit),j2,bet(nwrit),j3,
     &          isel(2),isel(3),isel(4)
             if (abs(iop1).gt.1) 
     &           write(47,'(i4,a3,i4)') abs(iop1),',1,',nwrit
             if (abs(iop2).gt.1) 
     &           write(47,'(i4,a3,i4)') abs(iop2),',2,',nwrit
             if (iop3.gt.1) 
     &           write(47,'(i4,a3,i4)') abs(iop3),',3,',nwrit
             if (iop3.lt.0) 
     &           write(47,'(i4,a4,i4)') abs(iop3),',14,',nwrit
          endif
      endif

3000  format(a2,1x,3(f11.6,1x,i1,1x),3x,3(i4,1x))
      return
      end

      logical function nolink(nwrit,ibl,nz)
      implicit double precision (a-h,o-z)
      dimension ibl(*)

      nolink = .true.

      do i=nwrit+1,nz
         if (iabs(ibl(i)).eq.nwrit) nolink = .false.
      end do

      return
      end

      subroutine mkvar(nwrit,ianzz,iz,bstr,astr,dstr)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      character*2 elemnt
      common /elem/elemnt(mxel)
      dimension ianzz(*),iz(4,*)

      character*2 tolowf
      character*(*) bstr
      character*(*) astr
      character*(*) dstr

      dimension isel(4),ian(4)


      isel(1) = nwrit
      do i=1,3
         isel(i+1) = iz(i,nwrit)
      end do
      do i=1,4
         ian(i) = ianzz(isel(i))
      end do

      if (nwrit.ge.2) then
             bstr = '        '
             n = 2
             bstr(1:n) = tolowf(elemnt(ian(1)))
             if (bstr(1:1).eq.' ') then
                bstr(1:1) = bstr(2:2)
                n = 1
             endif
c no axp  bstr(1:n+2) = bstr(1:n)//tolowf(elemnt(ian(2)))
             bstr = bstr(1:n)//tolowf(elemnt(ian(2)))
             n = n + 2
             if (bstr(n-1:n-1).eq.' ') then
                bstr(n-1:n-1) = bstr(n:n)
                n = n - 1
             endif
             if (nwrit.ge.3) then
                astr = '          '
                astr(1:n) = bstr(1:n)
             endif
             n = n + 1
             if (nwrit.lt.10) then
                write(bstr(n:n),'(i1)') nwrit
             elseif (nwrit.lt.100) then
                write(bstr(n:n+1),'(i2)') nwrit
             elseif (nwrit.lt.1000) then
                write(bstr(n:n+2),'(i3)') nwrit
             else
                write(bstr(n:n+3),'(i4)') nwrit
             endif
      endif

      if (nwrit.ge.3) then
             n = n - 1
c no axp  astr(1:n+2) = astr(1:n)//tolowf(elemnt(ian(3)))
             astr = astr(1:n)//tolowf(elemnt(ian(3)))
             n = n + 2
             if (astr(n-1:n-1).eq.' ') then
                astr(n-1:n-1) = astr(n:n)
                n = n - 1
             endif
             n = n + 1
             if (nwrit.lt.10) then
                write(astr(n:n),'(i1)') nwrit
             elseif (nwrit.lt.100) then
                write(astr(n:n+1),'(i2)') nwrit
             elseif (nwrit.lt.1000) then
                write(astr(n:n+2),'(i3)') nwrit
             else
                write(astr(n:n+3),'(i4)') nwrit
             endif
      endif

      if (nwrit.ge.4) then
             dstr = 'dih    '
             n = 4
             if (nwrit.lt.10) then
                write(dstr(n:n),'(i1)') nwrit
             elseif (nwrit.lt.100) then
                write(dstr(n:n+1),'(i2)') nwrit
             elseif (nwrit.lt.1000) then
                write(dstr(n:n+2),'(i3)') nwrit
             else
                write(dstr(n:n+3),'(i4)') nwrit
             endif
      endif

      return
      end

      logical function parlea(ilead,isel,ispdb)
      implicit double precision (a-h,o-z)
      dimension isel(4)
      logical pcklin,pline
      logical debug

      debug = .false.
c
c     Find three atoms connected to ilead (Function prelea)
c     And parse zmatrix line (pline)
c
      parlea = .true.

c 
c to do: - current handeling of 3 atoms on a line is only correct
c          if they are exactly on a line angle123 = 180.0
c          otherwise you need to calculate another angle(or dihdral ?)
c        - after one parse in zmat-editor dummy is gone again, and
c          center of molecule not correct
c

      call prelea(iprel,ilead,isel,ispdb,0)
      if (iprel.eq.1) then
         
         if (debug) print*,'try lead without dummy'
         if (.not.pline(isel)) goto 100
         if (debug) print*,'lead without dummy'
         return
      else
         if (debug) print*,'try lead with dummy'
         call prelea(iprel,ilead,isel,ispdb,1)
         if (iprel.eq.0) goto 100
         if (debug) print*,'lead with dummy'
      endif
c
c     parse zmat line, check if dummy needs to be inserted
c
      if (.not.pcklin(isel)) goto 100

      return
100   parlea = .false.
      return
      end

      logical function parleh(ilead,isel,ispdb)
      implicit double precision (a-h,o-z)
      dimension isel(4)
      logical pcklin,pline
      logical debug

      debug = .false.
c
c     Find three atoms connected to ilead (Function preleh)
c     And parse zmatrix line (pline)
c     this version of the routine is for Hydrogen only
c
      parleh = .true.

c 
c to do: - current handeling of 3 atoms on a line is only correct
c          if they are exactly on a line angle123 = 180.0
c          otherwise you need to calculate another angle(or dihdral ?)
c        - after one parse in zmat-editor dummy is gone again, and
c          center of molecule not correct
c

      call preleh(iprel,ilead,isel,ispdb,0)
      if (iprel.eq.1) then
         
         if (debug) print*,'try lead without dummy'
         if (.not.pline(isel)) goto 100
         if (debug) print*,'lead without dummy'
         return
      else
         if (debug) print*,'try lead with dummy'
         call preleh(iprel,ilead,isel,ispdb,1)
         if (iprel.eq.0) goto 100
         if (debug) print*,'lead with dummy'
      endif
c
c     parse zmat line, check if dummy needs to be inserted
c
      if (.not.pcklin(isel)) goto 100

      return
100   parleh = .false.
      return
      end

      logical function pcklin(isel)
      implicit double precision (a-h,o-z)
      dimension isel(4)
      logical pline
      real angle,bl
      common /zmatst/ nwrit,nvar,isimpl

      pcklin = .true.

      tol = 1.0d-10
      toang = 0.52917706d0

      if (nwrit.gt.1) then
         call tomold(angle,isel,3)
         if ((abs(angle).lt.tol.or.abs(angle).gt.180.e0-tol)) then
            call calcx(ical,isel,nx)
            if (ical.eq.0) goto 100
            itmp = isel(1)
            isel(1) = nx
            call dumlin(isel,1.0d0,90.0d0,0.0d0)
            isel(1) = itmp
            isel(4) = isel(3)
            isel(3) = nx
            call tomold(bl,isel,2)
            call dumlin(isel,bl*toang,90.0d0,180.0d0)
         else
            if (.not.pline(isel)) goto 100
         endif
      else
         if (.not.pline(isel)) goto 100
      endif

      return

100   pcklin = .false.
      return
      end

      subroutine prelead(iprel,ilead,isel,ispdb,ithree,
     &                        ianz,iaton,iresid,iconn,lwrit)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      common /zmatst/ nwrit,nvar,isimpl
      real angle,tol
      logical doit
      dimension isel(4),iscr(4)
      dimension ianz(*),iaton(*),iresid(*),iconn(mxcon+1,*),lwrit(*)


      tol = 0.1e0

      isel(1) = ilead
      isel(2) = 0
      isel(3) = 0
      isel(4) = 0

      itmp = 0
      jtmp = 0

      nanz = 1
      nfnd = 1

10    if (nwrit.gt.0) then
         do i=1,iconn(1,isel(1))
            ii = abs(iconn(i+1,isel(1)))
            if (ianz(ii).ne.nanz.and.lwrit(ii).ne.0
     &      .and.ianz(ii).lt.99.and.iaton(ii).eq.2.and.
     &      (ispdb.eq.0.or..not.(ispdb.eq.1.and.iresid(ii).gt.0)))
     &      then
               if (nwrit.eq.1) then
                  isel(2) = ii
                  nfnd = 2
                  goto 100
               endif
               do j=1,iconn(1,ii)
                  jj = abs(iconn(j+1,ii))
                  if (ianz(jj).ne.nanz.and.lwrit(jj).ne.0
     &                .and.jj.ne.isel(1).and.iaton(jj).eq.2.and.
     &      (ispdb.eq.0.or..not.(ispdb.eq.1.and.iresid(jj).gt.0)))
     &               then

                     doit = .false.
                     if (ithree.eq.0) then
                        iscr(1) = isel(1)
                        iscr(2) = ii
                        iscr(3) = jj
                        iscr(4) = 0
                        call tomold(angle,iscr,3)
                        if (abs(angle).ge.tol.and.
     &                      abs(angle).lt.180.e0-tol) then
                            doit = .true.
                        endif
                     else
                        doit = .true.
                     endif
                        
                     if (doit) then

                       if (nwrit.eq.2) then
                          isel(2) = ii
                          isel(3) = jj
                          nfnd = 3
                          goto 100
                       else
                          itmp = ii
                          jtmp = jj
                       endif
                     
                       do k=1,iconn(1,jj)
                          kk = abs(iconn(k+1,jj))
                          if (ianz(kk).ne.nanz.and.lwrit(kk).ne.0
     &                        .and.kk.ne.ii.and.iaton(kk).eq.2.and.
     &      (ispdb.eq.0.or..not.(ispdb.eq.1.and.iresid(kk).gt.0)))
     &                    then
                             iscr(1) = ii
                             iscr(2) = jj
                             iscr(3) = kk
                             iscr(4) = 0
                             call tomold(angle,iscr,3)
                             if (abs(angle).ge.tol.and.
     &                           abs(angle).lt.180.e0-tol) then
                                isel(2) = ii
                                isel(3) = jj
                                isel(4) = kk
                                nfnd = 4
                                goto 100
                             endif
                          endif
                       end do

                     endif
                  endif
               end do
            endif
         end do
      endif
100   if (nanz.eq.1.and.nfnd.lt.4.and.nfnd.lt.nwrit+1) then
c
c     go back once more and allow hydrogens
c
         nanz = 0
         goto 10
      endif
      if (isel(4).eq.0.and.itmp.ne.0) then
         do i=1,iconn(1,itmp)
            ii = abs(iconn(i+1,itmp))
            if (lwrit(ii).ne.0.and.iaton(ii).eq.2.and.ii.ne.jtmp) 
     &      then
                if (ithree.eq.1) then
                   isel(4) = ii
                else
                   iscr(1) = itmp
                   iscr(2) = jtmp
                   iscr(3) = ii
                   iscr(4) = 0
                   call tomold(angle,iscr,3)
                   if (abs(angle).ge.tol.and.
     &                 abs(angle).lt.180.e0-tol) then
                      isel(2) = itmp
                      isel(3) = jtmp
                      isel(4) = ii
                   endif
                endif
            endif
         end do
      endif
      if ((nwrit.eq.3.or.nwrit.eq.4).and.isel(4).eq.0.and.ithree.eq.1) 
     & then
         isel(2) = itmp
         isel(3) = jtmp
         do i=1,iatoms
            if (lwrit(i).ne.0) then
               k = 0
               do j=1,4
                  if (isel(j).eq.i) k = 1
               end do
               if (k.eq.0) isel(4) = i
            endif
         end do 
      endif

      iprel = 1

      nchk = nwrit + 1
      if (nchk.gt.4) nchk = 4
      do i=1,nchk
         if (iprel.eq.1.and.isel(i).ne.0) then
             iprel = 1
         else
             iprel = 0
         endif
      end do

      return
      end

      subroutine prelehd(iprel,ilead,isel,ispdb,ithree,
     &                        ianz,iaton,iresid,iconn,lwrit)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      common /zmatst/ nwrit,nvar,isimpl
      real angle,tol
      logical doit
      dimension isel(4),iscr(4),nn(mxcon+1)
      dimension inat(mxcon+1)
      dimension ianz(*),iaton(*),iresid(*),iconn(mxcon+1,*),lwrit(*)


      tol = 0.1e0

      isel(1) = ilead
      isel(2) = 0
      isel(3) = 0
      isel(4) = 0

      itmp = 0
      jtmp = 0

      nanz = 1
      nfnd = 1

10    if (nwrit.gt.0) then
         do i=1,iconn(1,isel(1))
            ii = abs(iconn(i+1,isel(1)))
            mm = iconn(1,ii)
            if (mm.gt.2) then
               call cnvcon(iconn,inat,lwrit,ii,nn)
            else
               nn(1) = mm
               do jj=1,mm
                  nn(1+jj) = iconn(1+jj,ii)
               end do
            endif

            if (ianz(ii).ne.nanz.and.lwrit(ii).ne.0
     &      .and.ianz(ii).lt.99.and.iaton(ii).eq.2.and.
     &      (ispdb.eq.0.or..not.(ispdb.eq.1.and.iresid(ii).gt.0)))
     &      then
               if (nwrit.eq.1) then
                  isel(2) = ii
                  nfnd = 2
                  goto 100
               endif
               do j=1,nn(1)
                  jj = abs(nn(j+1))
                  if (ianz(jj).ne.nanz.and.lwrit(jj).ne.0
     &                .and.jj.ne.isel(1).and.iaton(jj).eq.2.and.
     &      (ispdb.eq.0.or..not.(ispdb.eq.1.and.iresid(jj).gt.0)))
     &               then

                     doit = .false.
                     if (ithree.eq.0) then
                        iscr(1) = isel(1)
                        iscr(2) = ii
                        iscr(3) = jj
                        iscr(4) = 0
                        call tomold(angle,iscr,3)
                        if (abs(angle).ge.tol.and.
     &                      abs(angle).lt.180.e0-tol) then
                            doit = .true.
                        endif
                     else
                        doit = .true.
                     endif
                        
                     if (doit) then

                       if (nwrit.eq.2) then
                          isel(2) = ii
                          isel(3) = jj
                          nfnd = 3
                          goto 100
                       else
                          itmp = ii
                          jtmp = jj
                       endif
                     
                       do k=1,iconn(1,ii)
                          kk = abs(iconn(k+1,ii))
                          if (ianz(kk).ne.nanz.and.lwrit(kk).ne.0
     &                        .and.kk.ne.jj.and.kk.ne.isel(1)
     &                        .and.iaton(kk).eq.2.and.
     &      (ispdb.eq.0.or..not.(ispdb.eq.1.and.iresid(kk).gt.0)))
     &                    then
                             iscr(1) = ii
                             iscr(2) = jj
                             iscr(3) = kk
                             iscr(4) = 0
                             call tomold(angle,iscr,3)
                             if (abs(angle).ge.tol.and.
     &                           abs(angle).lt.180.e0-tol) then
                                isel(2) = ii
                                isel(3) = jj
                                isel(4) = kk
                                nfnd = 4
                                goto 100
                             endif
                          endif
                       end do

                       do k=1,iconn(1,jj)
                          kk = abs(iconn(k+1,jj))
                          if (ianz(kk).ne.nanz.and.lwrit(kk).ne.0
     &                        .and.kk.ne.ii.and.
     &                        jj.ne.kk.and.iaton(kk).eq.2.and.
     &      (ispdb.eq.0.or..not.(ispdb.eq.1.and.iresid(kk).gt.0)))
     &                    then
                             iscr(1) = ii
                             iscr(2) = jj
                             iscr(3) = kk
                             iscr(4) = 0
                             call tomold(angle,iscr,3)
                             if (abs(angle).ge.tol.and.
     &                           abs(angle).lt.180.e0-tol) then
                                isel(2) = ii
                                isel(3) = jj
                                isel(4) = kk
                                nfnd = 4
                                goto 100
                             endif
                          endif
                       end do

                     endif
                  endif
               end do
            endif
         end do
      endif
100   if (nanz.eq.1.and.nfnd.lt.4.and.nfnd.lt.nwrit+1) then
c
c     go back once more and allow hydrogens
c
         nanz = 0
         goto 10
      endif
      if (isel(4).eq.0.and.itmp.ne.0) then
         do i=1,iconn(1,itmp)
            ii = abs(iconn(i+1,itmp))
            if (lwrit(ii).ne.0.and.iaton(ii).eq.2.and.ii.ne.jtmp) 
     &      then
                if (ithree.eq.1) then
                   isel(4) = ii
                else
                   iscr(1) = itmp
                   iscr(2) = jtmp
                   iscr(3) = ii
                   iscr(4) = 0
                   call tomold(angle,iscr,3)
                   if (abs(angle).ge.tol.and.
     &                 abs(angle).lt.180.e0-tol) then
                      isel(2) = itmp
                      isel(3) = jtmp
                      isel(4) = ii
                   endif
                endif
            endif
         end do
      endif
      if ((nwrit.eq.3.or.nwrit.eq.4).and.isel(4).eq.0.and.ithree.eq.1) 
     & then
         isel(2) = itmp
         isel(3) = jtmp
         do i=1,iatoms
            if (lwrit(i).ne.0) then
               k = 0
               do j=1,4
                  if (isel(j).eq.i) k = 1
               end do
               if (k.eq.0) isel(4) = i
            endif
         end do 
      endif

      iprel = 1

      nchk = nwrit + 1
      if (nchk.gt.4) nchk = 4
      do i=1,nchk
         if (iprel.eq.1.and.isel(i).ne.0) then
             iprel = 1
         else
             iprel = 0
         endif
      end do

      return
      end

      subroutine wrzmaz(iun,iopt,
     &          bl,alph,bet,ibl,ialph,ibet,imap,janz,iz,epoints,ianz)

c this is really wrzmat

      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      parameter (mxbas=17)
      parameter (mxtsk=13)
      parameter (mxmth=8)
      parameter (mxvec=7)
      common /athlp/  iatoms, mxnat
      common /zmfrst/ ihaszm, nz, mxzat
      common /zmatst/ nwrit,nvar,isimpl
      common /pnthlp/ ipoints,ipnt
      character*2 elemnt,tolowf,tstr
      common /elem/elemnt(mxel)
      common /geocnv/ fdum3(10),idum(7),ieav,ifrav,mxpnt
      logical hassym
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      character*80 glin1,glin2,gtitl,rungam,freez
      character*15 jname,qname
      common /gauopt/ ito,imo,ibo,itotc,imult,ibatch,ihess,itime,
     &                iwxyz,ichh,ichm,ichl,imh,imm,iml,iexk,
     &                glin1,glin2,gtitl,jname,qname,rungam
      character*80 lstr
      character*5 hstr,htstr
      character*3 sr
      character*7 bases
      character*8 basnw
      character*3 met
      common /sets/ nbas(6,mxbas),bases(mxbas)

      character*14 tasks
      common /tsk/ tasks(mxtsk)
      character*10 meth
      common /methods/ meth(mxmth)
      character*8 vecstr
      common /vstrs/ vecstr(mxvec)
      logical dblvec,domap
      common /zmtyp/ igztyp
      common /typoni/ ioniad
      common /types/ iff

      dimension bl(*),alph(*),bet(*),ibl(*),ialph(*),ibet(*),
     &          imap(*),janz(*),iz(4,*)
      dimension epoints(*),ianz(*),basnw(8)

      data bases/'sto3g  ','3-21g  ','4-31g  ','6-31g  ','6-31g* ',
     &           '6-31g**','mini1  ','mini4  ','midi1  ','midi4  ',
     &           'dz     ','tzv    ','tzvp   ','ecpmin ','ecpdz  ',
     &           'ecptzv ','ecptzvp'/
      data basnw/'sto3g   ','3-21g   ','4-31g   ','6-31g   ','6-31g*  ',
     &           '6-31g** ','6-31+G* ','6-311G**'/
      data (nbas(i,1),i=1,6) /1,5,9,13,19,19/
      data (nbas(i,2),i=1,6) /2,9,13,17,29,23/
      data (nbas(i,3),i=1,6) /2,9,13,17,29,23/
      data (nbas(i,4),i=1,6) /2,9,13,17,29,23/
      data (nbas(i,5),i=1,6) /2,15,19,0,0,0/
      data (nbas(i,6),i=1,6) /5,15,19,0,0,0/
      data (nbas(i,7),i=1,6) /1,5,9,10,0,19/
      data (nbas(i,8),i=1,6) /1,5,9,10,0,19/
      data (nbas(i,9),i=1,6) /2,9,13,11,23,23/
      data (nbas(i,10),i=1,6) /2,9,13,14,26,29/
      data (nbas(i,11),i=1,6) /2,10,18,0,0,0/
      data (nbas(i,12),i=1,6) /3,14,21,0,46,0/
      data (nbas(i,13),i=1,6) /6,20,27,0,52,0/
      data (nbas(i,14),i=1,6) /1,4,4,4,10,4/
      data (nbas(i,15),i=1,6) /2,8,8,8,20,8/
      data (nbas(i,16),i=1,6) /3,12,12,0,0,0/
      data (nbas(i,17),i=1,6) /6,18,18,0,0,0/

      data tasks/'scf           ','optimise      ','optxyz        ',
     &           'saddle        ','force         ','hessian       ',
     &           'ci            ',
     &           'polarisability','hyper         ','magnet        ',
     &           'raman         ','infrared      ','analyse       '/
      data meth/ 'rhf       ','direct rhf','uhf       ',
     &           'gvb 1     ','mp2       ','direct mp2',
     &           'mp3       ','casscf    '/
      data vecstr/ 'atoms   ','hcore   ','minguess',
     &             'extguess','alphas  ','1       ','1 2     '/

      open(unit=47,form='formatted',status='scratch')
      open(unit=50,form='formatted',status='scratch')

      nwrit = 0
      nvar = 0
      hassym = .false.
      domap = .false.

      nglin1 = linlen(glin1)
      nglin2 = linlen(glin2)
      ngtitl = linlen(gtitl)

c MOPAC

      if (iopt.eq.3.or.iopt.eq.21) then
        do i=1,nz
           if (ibl(i).gt.1.or.ialph(i).gt.1.or.abs(ibet(i)).gt.1)
     &         hassym = .true.
        end do
        if (iopt.eq.21) then
           if (glin1(nglin1:nglin1).ne.'+') then
              write(iun,'(a)') glin1(1:nglin1)//' +'
           else
              write(iun,'(a)') glin1(1:nglin1)
           endif
           if (hassym) then
              write(iun,'(a)') glin2(1:nglin2)//' SYMMETRY'
           else
              write(iun,'(a)') glin2(1:nglin2)
           endif
           write(iun,'(a)') gtitl(1:ngtitl)
        else
           if (hassym) then
               write(iun,'(a)') 'SYMMETRY'
           else
               write(iun,'(a)') ' '
           endif
           if (ieav.eq.1) then
            if (iftyp.ge.2.and.iftyp.le.4) then
               write(iun,'(a,f12.6)') 'scf done: ',epoints(ipnt)
            else
               write(iun,'(a,f12.6)') 'FINAL HEAT OF FORMATION = ',
     &              epoints(ipnt)
            endif
           else
               write(iun,'(a)') ' '
           endif
        endif
        write(iun,'(a)') ' '
      endif

C NWCHEM

      if (iopt.eq.22.or.iopt.eq.7) then
        write(iun,'(a)') 'title '//char(34)//gtitl(1:ngtitl)//char(34)
        write(iun,'(a)') 'start '//gtitl(1:ngtitl)
        write(iun,'(a)') 'echo'
        write(iun,'(a,i3)') 'charge ',itotc
        write(iun,'(a)') 'geometry'
        if (iwxyz.eq.1) then
           if (icdex(glin1,'amber').ne.0) then
                 call wrcart(iun,0,2,1)
           else
                 call wrcart(iun,0,1,1)
           endif
           goto 889
        else
           write(iun,'(a)') 'zmatrix'
        endif
      endif

C GAMESS-UK

      if (iopt.eq.19) then
         if (ibatch.eq.0) write(iun,'(a,i5)') 'time ',itime
         if (imo.eq.6) write(iun,'(a)') 'memory 2000000'
         if (ito.ge.2.and.ito.le.4.and.ihess.eq.0) then
            write(iun,'(a)') 'dumpfile ed3 1500'
         endif
         if (ito.eq.13) write(iun,'(a)') 'restart'
         if (ito.eq.7) write(iun,'(a)') 'super off nosym'

         write(iun,'(a)') 'title'
         write(iun,'(a)') gtitl(1:ngtitl)
         write(iun,'(a,i3)') 'charge ',itotc
         write(iun,'(a,i2)') 'multiplicity ',imult
      endif

      if (iopt.eq.1.and.igztyp.eq.1) then
         nzt = nz*3 - 6
         if (nz.eq.1) nzt = 0
         if (nz.eq.2) nzt = 1
         if (nz.eq.3) nzt = 3

         write(iun,'(a,i4,a)') 
     &     ' $contrl runtyp=optimize coord=zmt nzvar=',nzt,' $end'
         write(iun,'(a)') 
     &     ' $basis  gbasis=N31 ngauss=6 ndfunc=1 npfunc=1 $end'

         nfrez = 0
         do i=1,nz
            if (ibl(i).eq.0.and.i.gt.1) nfrez = nfrez + 1
            if (ialph(i).eq.0.and.i.gt.2) nfrez = nfrez + 1
            if (ibet(i).eq.0.and.i.gt.3) nfrez = nfrez + 1
         end do

         if (nfrez.ne.0) then
            write(iun,'(a)') ' $statpt ifreez(1)= '
            freez = ' '
            nfrez = 0
            lfrez = 0
            do i=1,nz
               nv = (i-1)*3
               if (i.le.2) nv = 0
               if (i.eq.3) nv = 1
               if (i.ge.4) nv = nv - 6 
               if (ibl(i).eq.0.and.i.gt.1) then
                  nvt = nv + 1
                  htstr = hstr(nvt)
                  l = index(htstr,')')
                  sr = htstr(2:l-1)
                  freez = freez(1:lfrez) // sr(1:l-2) // ','
                  nfrez = nfrez + 1
                  lfrez = lfrez + l -1
               endif
               if (ialph(i).eq.0.and.i.gt.2) then
                  nvt = nv + 2
                  htstr = hstr(nvt)
                  l = index(htstr,')')
                  sr = htstr(2:l-1)
                  freez = freez(1:lfrez) // sr(1:l-2) // ','
                  nfrez = nfrez + 1
                  lfrez = lfrez + l -1
               endif
               if (ibet(i).eq.0.and.i.gt.3) then
                  nvt = nv + 3
                  htstr = hstr(nvt)
                  l = index(htstr,')')
                  sr = htstr(2:l-1)
                  freez = freez(1:lfrez) // sr(1:l-2) // ','
                  nfrez = nfrez + 1
                  lfrez = lfrez + l -1
               endif
               if (nfrez.gt.17) then
                  call spatrm(freez,l)
                  if (i.eq.nz) then
                      if (freez(l:l).eq.',') freez(l:l) = ' '
                  endif
                  write(iun,'(a)') freez(1:l)
                  freez = ' '
                  nfrez = 0
                  lfrez = 0
               endif
            end do

            call spatrm(freez,l)
            if (freez(l:l).eq.',') freez(l:l) = ' '
            write(iun,'(a)') freez(1:l)// '$end'
         endif
         write(iun,'(a)') ' $data'
         write(iun,'(a)') 'title'
         write(iun,'(a)') 'C1'
      endif

      if ((iopt.eq.1.or.iopt.eq.19).and.igztyp.ne.1) 
     &   write(iun,'(a)') 'zmat angstroms'

C GAUSSIAN

      if (iopt.eq.20) then
         write(iun,'(a)') '$ RunGauss'
         write(iun,'(a)') glin1(1:nglin1)
         write(iun,'(a)') glin2(1:nglin2)
         write(iun,'(a)') ' '
         write(iun,'(a)') gtitl(1:ngtitl)
         write(iun,'(a)') ' '

         ioniom = 1
         do i = 1,nglin1
           if (glin1(i:i).eq.':') ioniom = ioniom + 1
         end do
         if (ioniom .eq. 1) then
            write(iun,'(i3,1x,i2)') itotc,imult
         elseif (ioniom.eq.2) then
            write(iun,'(i3,1x,i2,1x,i3,1x,i2)') ichh,imh,ichl,iml
         else
            write(iun,'(i3,1x,i2,1x,i3,1x,i2,1x,i3,1x,i2)') ichh,imh,
     &          ichm,imm,ichl,iml
         endif
         if (iwxyz.eq.1) then
              if (icdex(glin1,'amber').ne.0) then
                 call wrcart(iun,0,2,1)
              else
                 call wrcart(iun,0,1,1)
              endif
              write(iun,'(a)') ' '
              goto 889
         endif
      endif

      iiopt = iopt
      if (iopt.eq.19) iiopt = 1
      if (iopt.eq.20) iiopt = 2
      if (iopt.eq.21) iiopt = 3
      if (iopt.eq.22) iiopt = 1

      igamb = 0
      if (icdex(glin1,'amber').ne.0) then
          igamb = 1
          if (iff.ne.3) then
             iff = 3
             ioniad = 1
             call dotyp(0)
             ioniad = 0
          endif
      endif

      do i=1,nz
          call wline(iun,iiopt,igamb)
      end do

      rewind(47)
      if (iiopt.eq.1) then
         if (igztyp.eq.1) then
            write(iun,'(a)') ' '
         else
            if (.not.(iopt.eq.22.and.nz.eq.1)) 
     &           write(iun,'(a)') 'variables'
         endif
      endif
      if (iiopt.eq.2) write(iun,'(a)') ' '
      if (iiopt.eq.3.and.hassym) write(iun,'(a)') ' '
      do while (.true.)
         read(47,'(a)',end=888) lstr
         if (igztyp.eq.1) then
            ix = index(lstr,' ')
            if (ix.gt.0) lstr(ix:ix) = '='
         endif
         ilen = ifblen(lstr)
         if (ilen.ge.1) write(iun,'(a)') lstr(1:ilen)
      end do

888   rewind(50)
      if (iiopt.eq.1) then
         if (igztyp.eq.1) then
            write(iun,'(a)') ' '
            write(iun,'(a)') ' '
            write(iun,'(a)') ' $end'
         else
            write(iun,'(a)') 'constants'
         endif
      endif
      if (iiopt.eq.2) write(iun,'(a)') ' '
      do while (.true.)
         read(50,'(a)',end=889) lstr
         if (igztyp.eq.1) then
            ix = index(lstr,' ')
            if (ix.gt.0) lstr(ix:ix) = '='
         endif
         ilen = ifblen(lstr)
         if (ilen.ge.1) write(iun,'(a)') lstr(1:ilen)
      end do
889   continue

      if (iiopt.eq.1.and.igztyp.ne.1) write(iun,'(a)') 'end'
      if (iiopt.ge.2.or.iopt.eq.3) write(iun,'(a)') ' '
      close(47)
      close(50)

      if (domap) then
         write(iun,'(a)') 'map'
         do i=1,nz
            write(iun,'(i5)') imap(i)
         end do
      endif

      if (iopt.eq.22) then
        write(iun,'(a)') 'end'
        if (iexk.eq.0) write(iun,'(a)') 'ecce_print ecce.out'
        write(iun,'(a,a,a,a)') 'basis ',char(34)//"ao basis"//
     &                           char(34)//' cartesian print'
        if (ibo.eq.1) then
           write(iun,'(a)') '* library sto-3g'
        else
           write(iun,'(a,a,a,a)') '* library '//char(34)//
     &                          basnw(ibo)//char(34)
        endif
        write(iun,'(a)') 'end'

        if (imo.gt.2) then
           write(iun,'(a)') 'dft'
           if (imult.eq.2) write(iun,'(a)') 'mult 2'
           if (imult.eq.3) write(iun,'(a)') 'mult 3'
           if (imo.eq.3) write(iun,'(a)') 'XC slater perdew81'
           if (imo.eq.4) write(iun,'(a)') 'XC xpbe96 cpbe96'
           if (imo.eq.5) write(iun,'(a)') 'XC b3lyp'
           write(iun,'(a)') 'end'
        endif

        if (imo.eq.1) then
            met = 'scf'
        else if (imo.eq.2) then
            met = 'uhf'
        else if (imo.gt.2) then
            met = 'dft'
        endif

        if (imult.ne.1.and.imo.eq.1) then
          write(iun,'(a)') 'scf'
          if (imult.eq.2) write(iun,'(a)') 'doublet'
          if (imult.eq.3) write(iun,'(a)') 'triplet'
          write(iun,'(a)') 'end'
        endif

        if (ito.eq.1) then
          write(iun,'(a,a)') 'task '//met//' energy'
        else if (ito.eq.2) then
          write(iun,'(a,a)') 'task '//met//' optimize'
        else if (ito.eq.3) then
          write(iun,'(a,a)') 'task '//met//' saddle'
        else if (ito.eq.4) then
          write(iun,'(a,a)') 'task '//met//' freq'
        endif

      endif

      if (iopt.eq.19) then
         write(iun,'(a,a)') 'basis ',bases(ibo)
         if (ibo.ge.14.and.ibo.le.17) then
            write(iun,'(a)') 'pseudo ecp'
            do i=1,iatoms
               tstr = tolowf(elemnt(ianz(i)))
               if (ianz(i).ne.1.and.ianz(i).ne.99)
     &            write(iun,'(a)') tstr//' '//tstr
            end do
         endif

         if (ito.ge.2.and.ito.le.4.and.ihess.eq.0) then
            write(iun,'(a,a,a)') 'runtype ',tasks(ito),' ed3 1'
         else
            write(iun,'(a,a)') 'runtype ',tasks(ito)
         endif

         if (ito.ne.7) write(iun,'(a,a)') 'scftype ',meth(imo)

         nelec = 0
         nnoth = 0
         do i=1,iatoms
            if (ianz(i).ne.99) nelec = nelec + ianz(i)
            if (ianz(i).ne.99.and.ianz(i).ne.1) nnoth = nnoth + 1
         end do
         nelec = nelec - itotc
         nelec2 = nelec / 2
         norb = nelec2
         if (nelec - nelec2*2.ne.0) norb = norb + 1

         if (imo.eq.8.and..not.ito.eq.7) then
            write(iun,'(a)') 'config'
            if (nnoth.gt.0) write(iun,'(a,i3,a,i3)') 
     &          'fzc ',1,' to ',nnoth
            write(iun,'(a,i3,a,i3)') 
     &          'doc ',nnoth+1,' to ',nelec2
            norb = nelec2
            if (nelec - nelec2*2.ne.0) then
                write(iun,'(a,i3)') 'alp ',norb
            endif
            write(iun,'(a,i3,a,i3)') 'uoc ',norb+1,' ',norb+2
            write(iun,'(a)') 'end'
            write(iun,'(a)') 'superci 1 to 4'
            write(iun,'(a)') 'newton 5 to 20'
            write(iun,'(a)') 'hessian 5 9 13 17'
         endif

         if (imo.eq.1.or.imo.eq.2.or.imo.eq.4.or.ito.eq.7) then
            if (imult.ne.1) then
                if (imult.eq.3) then
                   write(iun,'(a)') 'open 2 2'
                else
                   write(iun,'(a)') 'open 1 1'
                endif
            endif
         endif

         if (ito.eq.7) then
            nbasis = 0
            do i=1,iatoms
               ia = ianz(i)
               if (ia.ne.99) then
                   if (ia.ge.1.and.ia.le.2) then
                       igrp = 1
                   elseif (ia.ge.3.and.ia.le.10) then
                       igrp = 2
                   elseif (ia.ge.11.and.ia.le.18) then
                       igrp = 3
                   elseif (ia.ge.19.and.ia.le.20) then
                       igrp = 4
                   elseif (ia.ge.21.and.ia.le.30) then
                       igrp = 5
                   elseif (ia.ge.31) then
                       igrp = 6
                   endif
                   nbasis = nbasis + nbas(igrp,ibo)
               endif
            end do

            write(iun,'(a,i3,i3,i3)') 
     &         'direct ',nelec,norb,nbasis-norb
            if (imult.eq.2) write(iun,'(a)') 'spin doublet'
            if (imult.eq.3) write(iun,'(a)') 'spin triplet'
            write(iun,'(a)') 'conf'
            idm1 = 1
            idm2 = 2
            if (imult.eq.1) then
               write(iun,'(100(i3,1x))') (idm2,i=1,nelec2)
            endif
            if (imult.eq.2) then
               write(iun,'(100(i3,1x))') (idm2,i=1,nelec2),idm1
            endif
            if (imult.eq.3) then
               write(iun,'(100(i3,1x))') (idm2,i=1,nelec2-1),idm1,idm1
            endif
         endif

         if (ito.eq.13) then
            write(iun,'(a)') 'local'
            write(iun,'(a,i3,a)') '1 to ',norb,' end'
         endif

         dblvec = .false.
         if (imo.eq.3.or.imo.eq.4.or.imo.eq.8)
     &       dblvec = .true.
         if ((imo.eq.1.or.imo.eq.2).and.imult.ne.1) dblvec = .true.
         if (ito.eq.7.and.imult.ne.1) dblvec = .true.

         ivo = 2
         if (ibo.eq.1) ivo = 3
         if (ibo.ge.2.and.ibo.le.5) ivo = 4
         if (ito.eq.13) then
             ivo = 6
             if (dblvec) ivo = 7
         endif
         write(iun,'(a,a)') 'vectors ',vecstr(ivo)

         if (dblvec) then
            if (ito.eq.13) then
               write(iun,'(a,i3)') 'enter 3 4'
            else
               write(iun,'(a,i3)') 'enter 1 2'
            endif
         else
            if (ito.eq.13) then
               write(iun,'(a,i3)') 'enter 2'
            else
               write(iun,'(a,i3)') 'enter 1'
            endif
         endif

      endif
      
      return
      end

      subroutine wrcard(iun,dopdb,idogau,ipdbwh,epoints,
     &                  coo,qat,ianz,iaton,iresid,iconn,ityp,
     &                  ncalf,ianf,islu,nchain,iamino,reson,
     &                  irsnr,achain,ishoh)
c This is really wrcart
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxcon=10)
      parameter (mxel=100)
      parameter (mxres=42)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxchtp=136)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxamo=201)
      parameter (mxheta=150)
      common /athlp/  iatoms, mxnat
      common /pnthlp/ ipoints,ipnt
      character*2 elemnt,tocapf
      common /elem/   elemnt(mxel)
      common /geocnv/ fdum3(10),idum(7),ieav,ifrav,mxpnt
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /typoni/ ioniad
      common /cllab/  iclon,iclpnt(4)
      integer reson
      integer*2 ityp
      character*1 achain
      common /types/ iff
      character*3 chtnk,pdbsym,hsym,ambtnk
      character*2 amotnk,gffstr
      common /symbol/ pdbsym(mxsym),hsym(mxhsym),chtnk(mxchtp),
     &                ambtnk(mxamb),amotnk(mxamo),gffstr(mxgff)
      character*3 hetz
      common /clfstr/ ihashz,ihetq(mxheta),ihqset(mxheta),ihhadd(mxheta)
     &                ,labhet(mxheta),ilcset,ligcat(mxheta),hetz(mxheta)
      character*3 aminos
      common /amino/aminos(mxres)
      character*3 hez
      character*4 stmp
      character*14 ttmp
      integer ionilv,ism,isl,dopdb
      character*3 oniom(6), link

      dimension ipdb(mxsym),ihpdb(mxhsym*3),ich(3)
      dimension epoints(*),qat(*),ityp(*)
      dimension coo(3,*),ianz(*),iaton(*),iresid(*),iconn(mxcon+1,*)
      dimension ianf(*),islu(*),iamino(*),reson(*),irsnr(*),achain(*)

      data oniom /'H  ','M  ','L  ','H H','M H','L H'/

      ism = 0
      isl = 0

c idogau = 2, try to write amber type

      igaut = idogau
      if (igaut.eq.2.and.iff.ne.3) then
         iff = 3
         ioniad = 1
         call dotyp(0)
         ioniad = 0
      endif


      toang  = 0.52917706d0
      
      jatoms = iatoms
      if (iclon.eq.1) jatoms = iatoms - 8
      natoms = 0
      do i=1,jatoms
          if (ianz(i).lt.100) natoms = natoms + 1
          if ((ityp(i)/10000) .eq. 1) ism = 1
          if ((ityp(i)/10000) .eq. 2) isl = 1
      end do

      ionilv = 1 + ism + isl

      if (dopdb.eq.1) then
          irsnrt = 1
          do i=1,ncalf
             if (irsnr(i).le.0) irsnrt = 0
          end do
          if (irsnrt.eq.0) then
             do i=1,ncalf
                irsnr(i) = i
             end do
          endif
          irsnrt = 1
          do i=1,ncalf
             if (len(achain(i)).le.0) then
                irsnrt = 0
             else
                ic = ichar(achain(i))
                if (.not.((ic.ge.97.and.ic.le.122).or.
     &                    (ic.ge.65.and.ic.le.90))) irsnrt = 0
             endif
          end do
          if (irsnrt.eq.0) then
             do k=1,nchain
                do j=ianf(k),islu(k)
                   achain(j) = char(64+k)
                end do
             end do
          endif
      endif

      if (dopdb.eq.1) then
         write(iun,'(a)') 'HEADER'
      else
         if (igaut.eq.0) then
            write(iun,'(i5)') natoms
            if (ieav.eq.1) then
               iipnt = ipnt
               if (ipnt.eq.0) iipnt = 1
               if (iftyp.ge.2.and.iftyp.le.4) then
                  write(iun,'(a,f12.6)') 'scf done: ',epoints(iipnt)
               else
                  write(iun,'(a,f12.6)') 'FINAL HEAT OF FORMATION = ',
     &                                   epoints(iipnt)
               endif

            elseif (ionilv.ne.1) then
               write(iun,'(a5,1x,i1,1x,a6)') 'Oniom',ionilv,'levels'
            else
               write(iun,'(a)') ' '
            endif
         endif
      endif

      if (dopdb.eq.1) then

         iat = 0
         do i=1,ncalf
            if (ipdbwh.ne.2.or.(ipdbwh.eq.2.and.reson(i).eq.1)) then
               call getpdb(i,ipdb,ihpdb)
               do j=1,mxsym
                if (ipdb(j).ne.0) then
                   iat = iat + 1
                   stmp = ' '//pdbsym(j)
                   write(iun,'(a4,2x,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)') 
     &             'ATOM',iat,stmp,aminos(iamino(i)),achain(i),irsnr(i),
     &             (coo(k,ipdb(j))*toang,k=1,3)
                endif
               end do
               if (ipdbwh.eq.0.or.(ipdbwh.eq.2.and.reson(i).eq.1)) then
                  do j=1,mxhsym*3
                     if (ihpdb(j).ne.0) then
                        iat = iat + 1
                        ih = (j-1)/3 
                        il = j - ih*3
                        if (ihpdb(ih*3+2).eq.0) then
                           stmp = ' '//hsym(ih+1)
                        else
                           stmp = char(48+il)//hsym(ih+1)
                        endif
                        write(iun,
     &                  '(a4,2x,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)') 
     &                  'ATOM',iat,stmp,aminos(iamino(i)),achain(i),
     &                  irsnr(i),(coo(k,ihpdb(j))*toang,k=1,3)
                     endif
                  end do
               endif
            endif
         end do

         call hetnum(numhet)
         if (numhet.eq.0) numhet = 4

         do l=4,numhet
          do i=1,iatoms
            ir = iresid(i)
            if (ir.eq.-l.and.ir.ne.-ishoh.and.ianz(i).ne.100) then
              if (.not.((ipdbwh.eq.1.and.ianz(i).eq.1).or.
     &           (ipdbwh.eq.2.and.iaton(i).eq.0))) then
                 iat = iat + 1
                 ires = abs(iresid(i)) - 2
                 if (ihashz.eq.1) then
                    hez = hetz(ires-1+4)
                    if (hez.eq."   ") then
                       call gtht(ich,ires-1+4)
                       hez(1:1) = char(ich(1))
                       hez(2:2) = char(ich(2))
                       hez(3:3) = char(ich(3))
                    endif
                 else
                    write(hez,'(i3)') ires
                 endif
		 stmp = "    "
		 call gethet(i,nstat,stmp)
		 if (nstat.eq.0) stmp = tocapf(elemnt(ianz(i)))//'  '
                 write(iun,'(a6,i5,1x,a4,1x,a3,2x,i4,4x,3f8.3)') 
     &            'HETATM',iat,stmp,
     &             hez,ncalf+ires,(coo(j,i)*toang,j=1,3)
              endif
            endif
          end do
         end do

         do i=1,iatoms
            ir = iresid(i)
            if (ir.eq.-ishoh.and.ianz(i).ne.100) then
              if (.not.((ipdbwh.eq.1.and.ianz(i).eq.1).or.
     &           (ipdbwh.eq.2.and.iaton(i).eq.0))) then
                 iat = iat + 1
                 ires = abs(iresid(i)) - 2
                 if (ihashz.eq.1) then
                    hez = hetz(ires-1+4)
                    if (hez.eq."   ") then
                       call gtht(ich,ires-1+4)
                       hez(1:1) = char(ich(1))
                       hez(2:2) = char(ich(2))
                       hez(3:3) = char(ich(3))
                    endif
                 else
                    write(hez,'(i3)') ires
                 endif
		 stmp = "    "
		 call gethet(i,nstat,stmp)
		 if (nstat.eq.0) stmp = tocapf(elemnt(ianz(i)))//'  '
                 write(iun,'(a6,i5,1x,a4,1x,a3,2x,i4,4x,3f8.3)') 
     &            'HETATM',iat,stmp,
     &             hez,ncalf+ires,(coo(j,i)*toang,j=1,3)
              endif
            endif
         end do

      else

         do i=1,iatoms
           if (ianz(i).ne.100) then
              if (ihasq.eq.1.and.igaut.eq.0) then

                 write(iun,'(a2,1x,3(f12.6,1x),f9.6)') 
     &           elemnt(ianz(i)),(coo(j,i)*toang,j=1,3),qat(i)

              else

                 if (ionilv.eq.1) then

                    if (igaut.le.1) then
                       write(iun,'(a2,1x,3(f12.6,1x),1x,a1)')
     &                 elemnt(ianz(i)),(coo(j,i)*toang,j=1,3)
                    else
                       call wrlab(i,ttmp,ianz,ityp,qat)
                       write(iun,'(a14,1x,3(f12.6,1x),1x,a1)')
     &                 ttmp,(coo(j,i)*toang,j=1,3)
                    endif

                 else

                    iqmmm = ityp(i)/10000
                    link = oniom(iqmmm+1)
                    do l = 2,iconn(1,i)+1
                        j = iconn(l,i)
                        ij = ityp(j)/10000
                        if (j.ne.0 .and. ij.lt.iqmmm) 
     &                     link = oniom(iqmmm+4)
                    end do

                    if (igaut.le.1) then

                       write(iun,'(a2,1x,3(f12.6,1x),1x,a3)')
     &                 elemnt(ianz(i)),(coo(j,i)*toang,j=1,3),
     &                 link

                    elseif (igaut.eq.2) then

                       call wrlab(i,ttmp,ianz,ityp,qat)
                       write(iun,'(a14,1x,3(f12.6,1x),1x,a3)')
     &                 ttmp,(coo(j,i)*toang,j=1,3),
     &                 link

                    endif
                 endif
              endif
           endif
         end do

      endif

      return
      end

      subroutine wrpnt(filenm,lenfn,iopt,iapp,ipoints,
     &                 fancy,atcol,dolabs,backb)
      implicit double precision (a-h,o-z)
      common /zmfrst/ ihaszm, nz, mxzat
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb,nions,ntota,nresi
      common /vropt/  ivtwo,ihand,ivadd
      character*(*) filenm
      character*80 ftmp
      logical zmqok, zok
      integer dolabs,fancy,atcol,backb

c
c     iapp
c     -------
c     0    open, write and close
c     1    open, write
c     2          write and close
c     3          write
c
c
c     iopt
c     -------
c     1    zmat gamess
c     2    zmat gaussian
c     3    zmat mopac
c     4    cartesian coordinates
c     5    VRML
c     6    moldenogl
c     7    zmat nwchem
c
c     ixyz
c     -------
c     0    xyz file
c     1    mol file
c     2    msf file
c     3    tinker xyz file
c     4    tinker xyz file ?
c     5    pdb file
c     6    chemx format
c     7    msi format
c     8    molden format file 
c     9    mol2 file 
c     11   ambfor xyz file 
c     12   mopac xyz input file 

      iun = 51

      zok = zmqok(idum)

      if (lenfn.eq.0.or.filenm.eq.' ') then
         call zmterr('Invalid Filename !',1,0,1)
      else
         if (iapp.le.1) then
             if (iopt.eq.4.and.ixyz.eq.2) then
                open(unit=iun,form='unformatted',
     &          file=filenm,status='unknown',err=2000)
             else
                open(unit=iun,form='formatted',
     &          file=filenm,status='unknown',err=2000)
             endif
         endif

         if (iopt.eq.4) then

             if (ixyz.eq.1) then
                call wrmol(iun)
             elseif (ixyz.eq.2) then
                call wrmsf(iun)
             elseif (ixyz.eq.3) then
                call wrtnk(iun)
             elseif (ixyz.eq.4) then
                call wrtnk(iun+100)
             elseif (ixyz.eq.5) then
                call wrcart(iun,1,0,ipdbwh)
             elseif (ixyz.eq.6) then
                call wrchx(iun)
             elseif (ixyz.eq.7) then
                call wrmsi(iun)
             elseif (ixyz.eq.8) then
                call prtmolf(iun,ihaszm,ipoints)
             elseif (ixyz.eq.9) then
                call outmol(iun)
             elseif (ixyz.eq.11) then
                call appchg
                call wrgff(iun)
             elseif (ixyz.eq.12) then
                call mopxyz(iun)
             else
                call wrcart(iun,0,0,0)
             endif

         elseif ((iopt.ge.1.and.iopt.le.3).or.iopt.eq.7.or.
     &           (iopt.ge.19.and.iopt.le.22)) then

             if (zok) call wrzmat(iun,iopt)

         elseif (iopt.eq.5) then

             if (iapp.le.1) call plvhd(iun)
             call plvrml(iun,fancy,atcol,dolabs,1,backb,.false.)
             if (iapp.eq.0.or.iapp.eq.2) call plvend(iun,0)
             if (iapp.eq.4) call plvend(iun,1)

         elseif (iopt.eq.6) then

             if (iapp.le.1.and.ivtwo.eq.3) then
                 write(iun,'(a)') '[MOLDENOGL]'
             endif
             call wrogl(iun)

         endif
      endif

      if (iapp.eq.0.or.iapp.eq.2.or.iapp.eq.4) then
         close(iun)
         ftmp = 'Succesfully wrote file: '//filenm(1:lenfn)
         if (iopt.lt.19) call zmterr(ftmp,1,0,0)
      endif
      return

2000  call zmterr('Error Opening File !',1,0,1)
      return
      end

      subroutine mapzzz(iun,imod,iff,izmtmp,istatz,
     &              bl,alph,bet,ibl,ialph,ibet,imap,ianzz,iz,ianz)

c this is really logical function mapxyz

      implicit double precision (a-h,o-z)
      common /athlp/  iatoms, mxnat
      common /zmfrst/ ihaszm, nz, mxzat
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real bond,angle,dihed
      dimension isel(4)
      dimension bl(*),alph(*),bet(*),ibl(*),ialph(*),ibet(*),
     &          imap(*),ianzz(*),iz(4,*),ianz(*)

      istatz = 1
      toang = 0.52917706d0
      iuntmp = iun2
      iun2 = iun
      ierr = 0
      idummy = 0
      
      if (imod.eq.1) then
         call rdmol(0,itmp,0,istat)
         if (istat.eq.0) then
            ierr = 1
            goto 100
         endif
      elseif (imod.eq.2) then
         call rdmsf(0,istat)
         if (istat.eq.0) then
            ierr = 2
            goto 100
         endif
      elseif (imod.eq.3) then
         call gettnk(igttnk,0,idummy,iff,iheat,heat)
         if (igttnk.eq.0) then
            ierr = 3
            goto 100
         else
            call dohcon(0)
         endif
      elseif (imod.eq.8) then
         call tnkfst(igttnk,0,0)
         if (igttnk.eq.0) then
            ierr = 3
            goto 100
         else
            call dohcon(0)
         endif
      elseif (imod.eq.7) then
         call rdchx(0,1,0,0,0,istat,1)
         if (istat.eq.0) then
            ierr = 4
            goto 100
         endif
      elseif (imod.eq.0) then
         call getxyz(igetxy,heat,0)
         if (igetxy.eq.0) then
            ierr = 4
            goto 100
         endif
      endif

      if (izmtmp.eq.0) then
         iun2 = iuntmp
         return
      endif

      ii = 0
      do i=1,nz
         if (imap(i).gt.0) then
            if (ianz(imap(i)).ne.ianzz(i)) then
               ierr = 5
               goto 100
            endif
            ii = ii + 1
         endif
      end do

      if (ii.ne.iatoms) then
         ierr = 4
         goto 100
      endif

      do i=1,nz
      
        isel(1) = imap(i)
        do j=1,3
           isel(j+1) = imap(iz(j,i))
        end do

        if (i.ge.2) then
            call intcor(intc,bond,isel,2)
            if (intc.eq.0) goto 100
            bond = bond * toang
            bl(i) = bond
        endif
        if (i.ge.3) then
            call intcor(intc,angle,isel,3)
            if (intc.eq.0) goto 100
            alph(i) = angle
        endif
        if (i.gt.3) then
            call intcor(intc,dihed,isel,4)
            if (intc.eq.0) goto 100
            bet(i) = dihed
        endif
      
      end do

      iun2 = iuntmp
      ihaszm = 1
      call upzme
      return

100   ihaszm = 0
      iun2 = iuntmp
      call inferr('ERROR mapping/reading XYZ file !',0)
      call zmterr('ERROR mapping/reading XYZ file !',1,0,1)
      print*,'mapping ERROR ',ierr
      istatz = 0
      return
      end

c========================= For the convenience of C =======
c========================= Map functions to subroutines ===

      logical function chkmap(idum)
      implicit double precision (a-h,o-z)

      call chkmaz(istat)

      if (istat.eq.0) chkmap = .false.
      if (istat.eq.1) chkmap = .true.

      return
      end

      integer function fndmap(ixyz)
      implicit double precision (a-h,o-z)

      call fndmaz(ixyz,istat)

      fndmap = istat

      return
      end

      logical function pline(isel)
      implicit double precision (a-h,o-z)
      dimension isel(4)

      call plinzz(isel,istat)

      if (istat.eq.0) pline = .false.
      if (istat.eq.1) pline = .true.

      return
      end

      logical function mapxyz(iun,imod,iff,izmtmp)
      implicit double precision (a-h,o-z)

      call mapxzz(iun,imod,iff,izmtmp,istat)

      if (istat.eq.0) mapxyz = .false.
      if (istat.eq.1) mapxyz = .true.

      return
      end

      logical function getzmu(idum)
      implicit double precision (a-h,o-z), integer (i-n)

      call getzmz(istat)

      if (istat.eq.0) getzmu = .false.
      if (istat.eq.1) getzmu = .true.

      return
      end

      subroutine atmd(mopac,ipsi,adfin)
      implicit double precision (a-h,o-z)
      logical mopac,adfin

      if (mopac) then
         imo = 1
      else
         imo = 0
      endif
      if (adfin) then
         iao = 1
      else
         iao = 0
      endif

      call atdd(imo,ipsi,iao)

      return
      end

      subroutine denmak(denok)
      implicit double precision (a-h,o-z)
      logical denok

      call denmad(ido)

      if (ido.eq.1) then
         denok = .true.
      else
         denok = .false.
      endif

      return
      end

      subroutine muldma(vdwr,moddma,domul,dodma)
      implicit double precision(a-h,o-z)
      logical domul,dodma
      dimension vdwr(*)

      if (domul) then
         idm = 1
      else
         idm = 0
      endif
      if (dodma) then
         idd = 1
      else
         idd = 0
      endif

      call muldmd(vdwr,moddma,idm,idd)

      return
      end

      subroutine rdgam(idebug,befo,statio,ioxyz,irtype,istats)
      implicit double precision (a-h,o-z)
      logical befo,statio

      if (befo) then
         ibefo = 1
      else
         ibefo = 0
      endif

      call rdgad(idebug,ibefo,istatio,ioxyz,irtype,istats)

      if (istatio.eq.1) then
         statio = .true.
      else
         statio = .false.
      endif

      return
      end

      subroutine rdgamu(idebug,befo,statio,irtype,hesend,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      logical befo,statio,hesend

      if (befo) then
         ibefo = 1
      else
         ibefo = 0
      endif

      if (hesend) then
         ihsend = 1
      else
         ihsend = 0
      endif

      call rdgamd(idebug,ibefo,istatio,irtype,ihsend,istats)

      if (istatio.eq.1) then
         statio = .true.
      else
         statio = .false.
      endif

      return
      end

      subroutine rdcpmd(idebug,befo,statio,ioxyz,irtype,hesend,istats)
      implicit double precision (a-h,o-z), integer (i-n)
      logical befo,statio,hesend

      if (befo) then
         ibefo = 1
      else
         ibefo = 0
      endif

      if (hesend) then
         ihsend = 1
      else
         ihsend = 0
      endif

      call rdcpmdd(idebug,ibefo,istatio,ioxyz,irtype,ihsend,istats)

      if (istatio.eq.1) then
         statio = .true.
      else
         statio = .false.
      endif

      return
      end

      subroutine rdgaus(idebug,befo,statio,irtype,istats)
      implicit double precision (a-h,o-z)
      logical befo,statio

      if (befo) then
         ibefo = 1
      else
         ibefo = 0
      endif

      call rdgaud(idebug,ibefo,istatio,irtype,istats)

      if (istatio.eq.1) then
         statio = .true.
      else
         statio = .false.
      endif

      return
      end

      subroutine rdmaux(idebug,statio,istats)
      implicit double precision (a-h,o-z)
      logical statio

      call rdmaud(idebug,istatio,istats)

      if (istatio.eq.1) then
         statio = .true.
      else
         statio = .false.
      endif

      return
      end

      subroutine rdmolf(idebug,statio,irtype,iesp,istats)
      implicit double precision (a-h,o-z)
      logical statio

      call rdmold(idebug,istatio,irtype,iesp,istats)

      if (istatio.eq.1) then
         statio = .true.
      else
         statio = .false.
      endif

      return
      end

      subroutine rdnorb(naorbs,impas)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs, nat(numatm)
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /orbhlp/ mxorb,iuhf,ispd
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line
      logical datlin


      if (iftyp.eq.3) then

         call search(line,'NUMBER OF CARTESIAN GAUSSIAN BASIS',istat)
         if (istat.eq.1) then
             read(line,'(47x,i5)',err=1000) norbs
         else
             call rewfil
             call search(line,'NUMBER OF BASIS FUNCTIONS',istat)
             read(line,'(38x,i5)',err=1000) norbs
         endif

      elseif (iftyp.eq.2) then

         call search(line,'basis functions',istat)
         if (index(line,'cartesian').ne.0) then
            read(line,'(45x,i5)',err=1000) norbs
         else
            read(line,'(38x,i5)',err=1000) norbs
         endif

      elseif (iftyp.eq.4) then

         call search(line,'primitive gaussians',istat)
         if (istat.ne.0) then
            if (index(line,'primitive gaussians').eq.30) then
               read(line,'(1x,i3)',err=1000,end=1000) norbs
            else
               read(line,'(1x,i5)',err=1000,end=1000) norbs
            endif
         else
            goto 1000
         endif

      elseif (iftyp.eq.5) then

         norbs = naorbs

      elseif (iftyp.eq.1) then

         if (isbin.eq.1) then

             read(iun2,end=1000,err=1000)
     &       idum1,norbs,idum2,((rdum,j=1,idum1),i=1,3)
             if (idum1.gt.numatm) goto 1000 

         else

             call nxtlin(line,jstat)
             if (jstat.eq.1.or.jstat.eq.2) goto 1000

             if (.not.datlin(line)) goto 1000
             call bckfil
             if (index(line,"START OF MOPAC").ne.0) then
                call search(line,'AO_ZETA[',istat)
                if (istat.eq.1) then
                   norbs = intlin(line,istat)
                else
                   goto 1000
                endif
             else
                if (impas.eq.1) then
                   
                   call nxtlin(line,jstat)
                   if (jstat.eq.1.or.jstat.eq.2) goto 1000
                   read(line,*,end=1000,err=1000)
     &             idum1,norbs,idum2,((rdum,j=1,idum1),i=1,3)
                   if (idum1.gt.numatm) goto 1000 
                else
                   call nxtlin(line,jstat)
                   if (jstat.eq.1.or.jstat.eq.2) goto 1000
                   read(line,'(3i8)',err=1000)
     &                  idum1,norbs,idum2
                   do i=1,idum1
                      call nxtlin(line,jstat)
                      if (jstat.eq.1.or.jstat.eq.2) goto 1000
                      read(line,'(1x,3f16.8)',err=1000)
     &                   (rdum,j=1,3)
                   end do
                   if (idum1.gt.numatm) goto 1000 
                endif
             endif
         endif

      endif

      if (norbs.gt.mxorb) then
         nsiz = norbs/256 
         lsiz = norbs - nsiz*256
         if (lsiz.gt.0) nsiz = nsiz + 1
         nsiz = nsiz*256
         call allorb(nsiz,0)
         print*,'Allocating memory for orbitals !'
      endif

1000  call rewfil
      return
      end

      subroutine wrlab(i,ttmp,ianz,ityp,qat)
      implicit double precision (a-h,o-z)

      parameter (mxel=100)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxchtp=136)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      parameter (mxamo=201)

      character*2 elemnt
      character*(*) ttmp
      common /elem/   elemnt(mxel)
      integer*2 ityp,it10000
      character*3 chtnk,pdbsym,hsym,ambtnk
      character*2 amotnk,gffstr
      character*7 chgstr
      common /symbol/ pdbsym(mxsym),hsym(mxhsym),chtnk(mxchtp),
     &                ambtnk(mxamb),amotnk(mxamo),gffstr(mxgff)
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb,nions,ntota,nresi
      dimension ityp(*),ianz(*),qat(*)

      it10000 = 10000

      ttmp = '              '
      ttmp = elemnt(ianz(i))
c      if (ianz(i).eq.16) ttmp(2:2) = 'p'
      im = mod(ityp(i),it10000)

      if (im.gt.0) then
         ttmp(3:3) = '-'
         ttmp = ttmp(1:3)//ambtnk(im)
         if (iambch.eq.1) then
            ttmp(7:7) = '-'
            write(chgstr,'(f7.4)') qat(i)
            ttmp = ttmp(1:7)//chgstr
         endif
      endif

      call tocap(ttmp,14)
      if (ttmp(4:5).eq.'H1'.and.
     &   (im.ge.350.and.im.le.500)) ttmp(5:5) = 'P'

      l = 14
      call spatrm(ttmp,l)

      return
      end

      subroutine cnvcon(iconn,inat,lwrit,ii,nn)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      dimension nn(mxcon+1),icnn(mxcon+1),lcnn(mxcon+1),lb(mxcon+1)
      dimension inat(mxcon+1),iconn(mxcon+1,*),lwrit(*)

      ncnn = iconn(1,ii)
      mcnn = 0

      do j=1,ncnn

         icnn(j) = abs(iconn(j+1,ii))
         jj = lwrit(icnn(j))

         if (jj.gt.0) then
             mcnn = mcnn + 1
             lcnn(mcnn) = jj
             lb(mcnn) = j
         endif

      end do

      call srti(mcnn,lcnn,inat)

      nn(1) = mcnn
      do j=1,mcnn
         nn(j+1) = icnn(lb(inat(j)+1))
      end do

      return
      end
