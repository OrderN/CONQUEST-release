      subroutine hcoodd(istat,ioatms,nstrt,
     &                  ipdbt,coo,ianz,iaton,iresid,iconn,
     &                  icalf,ncalf,ianf,islu,nchain,iamino)
c hcoord, hbond en vadar hebben geen ioadd provisie
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat

      parameter (numcal=50000)
      common /hbonds/ hbd(2,numcal),ihb(2,numcal)
      integer*2 ipdbt
      logical hashy
      dimension tmp(3),ipdbt(*)
      dimension coo(3,*),ianz(*),iaton(*),iresid(*),iconn(mxcon+1,*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*)

      istat = 1

c     nhbond in A.U. Atomic Units

      if (iatoms+ncalf.le.mxnat) then
         nhatm = iatoms
         do k=1,nchain
c
c           first residue
c
          if (iamino(ianf(k)).le.23.and.ianf(k).ge.nstrt) then
            hashy = .false.
            do n=1,iconn(1,icalf(2,ianf(k)))
               nn = abs(iconn(n+1,icalf(2,ianf(k))))
               if (ianz(nn).eq.1.and.nn.gt.ioatms) then
                  hashy = .true.
                  icalf(4,ianf(k)) = nn
                  iaton(nn) = 0
                  iresid(nn) = ianf(k)
               endif
            end do
            if (.not.hashy.and.iamino(ianf(k)).ne.15.and.
     &         iamino(ianf(k)).ne.0) then
               call bckok(ibckok,ianf(k),2)
               if (ibckok.eq.1) then
                  nhatm = nhatm + 1
                  icalf(4,ianf(k)) = nhatm
                  do j=1,3
                     tmp(j) = coo(j,icalf(1,ianf(k))) 
     &                      - coo(j,icalf(3,ianf(k)))
                  end do
                  tmpl = vlen(tmp)
                  do j=1,3
                     coo(j,nhatm) = coo(j,icalf(2,ianf(k))) 
     &                              +(tmp(j)/tmpl)*1.89d0
                  end do
                  iaton(nhatm) = 0
                  ianz(nhatm) = 1
                  iconn(1,nhatm) = 1
                  m = icalf(2,ianf(k))
                  iconn(2,nhatm) = m
                  iconn(1,m) = iconn(1,m) + 1
                  iconn(1+iconn(1,m),m) = nhatm
                  iresid(nhatm) = ianf(k)
                  ipdbt(nhatm) = 1
               endif
            endif
            ihb(1,ianf(k)) = 0
            ihb(2,ianf(k)) = 0
           endif
c
c           other residues
c
           if (ianf(k).ge.nstrt) then
            do i=ianf(k)+1,islu(k)
             if (iamino(i).le.23) then
               icalf(4,i) = 0
               hashy = .false.
               do n=1,iconn(1,icalf(2,i))
                  nn = abs(iconn(n+1,icalf(2,i)))
                  if (ianz(nn).eq.1.and.nn.gt.ioatms) then
                     hashy = .true.
                     icalf(4,i) = nn
                     iaton(nn) = 0
                     iresid(nn) = i
                  endif
               end do

               if (.not.hashy.and.iamino(i).ne.15) then

                  m = icalf(2,i)
                  idxc = icalf(3,i-1)
                  idxo = 0

                  do j=1,iconn(1,idxc)
                     kk = abs(iconn(j+1,idxc))
                     if (ianz(kk).eq.8.and.kk.gt.ioatms) idxo = kk
                  end do

                  call bckok(ibckok,i,2)

                  if (idxo.ne.0.and.ibckok.eq.1) then

                     nhatm = nhatm + 1
                     icalf(4,i) = nhatm

                     do j=1,3
                        tmp(j) = coo(j,idxc) - coo(j,idxo)
                     end do

                     tmpl = vlen(tmp)

                     do j=1,3
                        coo(j,nhatm) = coo(j,m) +
     &                             (tmp(j) / tmpl)*1.89d0
                     end do

                     iaton(nhatm) = 0
                     ianz(nhatm) = 1
                     iconn(1,nhatm) = 1
                     iconn(2,nhatm) = m
                     iconn(1,m) = iconn(1,m) + 1
                     iconn(1+iconn(1,m),m) = nhatm
                     iresid(nhatm) = i
                     ipdbt(nhatm) = 1
                  else
                     icalf(4,i) = 0
                  endif

               endif

               ihb(1,i) = 0
               ihb(2,i) = 0
              endif
            end do
           endif
         end do
         iatoms = nhatm
      else
         istat = 0
         do i=1,ncalf
            if (iamino(i).le.23) icalf(4,i) = 0
         end do
      endif

      return
      end

      double precision function phi(idx,icalf,ncalf,ianf,nchain)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      real dihed
      dimension isel(4)
      dimension icalf(6,*),ianf(*)

      phi = 0.0d0

      do i=1,nchain
        if (idx.eq.ianf(i)) return
      end do

      if (idx.le.ncalf.and.idx.gt.1) then
         call bckok(ibckok,idx,2)
         if (ibckok.eq.1) then
            isel(1) = icalf(3,idx-1)
            isel(2) = icalf(2,idx)
            isel(3) = icalf(1,idx)
            isel(4) = icalf(3,idx)
            call tomold(dihed,isel,4)
            phi = dihed
         endif
      endif

      return
      end

      double precision function psi(idx,icalf,ncalf,islu,nchain)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      real dihed
      dimension isel(4)
      dimension icalf(6,*),islu(*)

      psi = 0.0d0

      do i=1,nchain
        if (idx.eq.islu(i)) return
      end do

      if (idx.lt.ncalf.and.idx.ge.1) then
         isel(1) = icalf(2,idx)
         isel(2) = icalf(1,idx)
         isel(3) = icalf(3,idx)
         isel(4) = icalf(2,idx+1)
         call tomold(dihed,isel,4)
         psi = dihed
      endif

      return
      end

      subroutine hand(idx1,idx2,hang,coo,ianz,iconn,icalf)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (numcal=50000)
      common /hbonds/ hbd(2,numcal),ihb(2,numcal)
      dimension tmp(3),tem(3),coo(3,*),ianz(*),iconn(mxcon+1,*)
      dimension icalf(6,*)

      idxc = icalf(3,idx1)
      idxo = 0
      do j=1,iconn(1,idxc)
         k = abs(iconn(j+1,idxc))
         if (ianz(k).eq.8) idxo = k
      end do
      if (idxo.eq.0) then
          print*,'subroutine hang: c without o'
          hang = 0.0d0
      else
         do j=1,3
            tmp(j) = coo(j,idxo) - coo(j,idxc)
         end do
   
         idxn = icalf(2,idx2)
         do j=1,3
            tem(j) = coo(j,idxn) - coo(j,icalf(4,idx2))
         end do

         call impsc(tmp,tem,cosb)
         hang = cosb
      endif
      
      return
      end

      integer function indhb(idx,distnw)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (numcal=50000)
      common /hbonds/ hbd(2,numcal),ihb(2,numcal)

      dismax = 0.0d0
      do j=1,2
         if (ihb(j,idx).eq.0) then
            indhb = j
            return
         elseif (hbd(j,idx).gt.dismax) then
            dismax = hbd(j,idx)
            maxid = j
         endif
      end do
      indhb = -1
      if (distnw.lt.hbd(maxid,idx)) indhb = maxid

      return
      end

      subroutine hbodd(ioatms,nstrt,coo,ianz,iconn,
     &                 icalf,ncalf,iamino)
c hcoord, hbond en vadar hebben geen ioadd provisie
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (mxcon=10)
      parameter (numcal=50000)
      common /hbonds/ hbd(2,numcal),ihb(2,numcal)
      dimension tmp(3),tem(3)
      dimension coo(3,*),ianz(*),iconn(mxcon+1,*)
      dimension icalf(6,*),iamino(*)
  
      todeg=45.0d0/datan(1.0d0)

      do k=nstrt,ncalf
         dismin = 4.0d0/0.52917706d0
         do i=nstrt,ncalf
           if ((icalf(4,i).ne.0.and.icalf(3,k).ne.0)
     &      .and.iamino(i).le.23) then
            idxc = icalf(3,k)
            idxo = 0
            do j=1,iconn(1,idxc)
               m = abs(iconn(j+1,idxc))
               if (ianz(m).eq.8.and.m.gt.ioatms) idxo = m
            end do
            do j=1,3
               tmp(j) = coo(j,idxo) - coo(j,icalf(4,i))
            end do
            dho = vlen(tmp)
            do j=1,3
               tem(j) = coo(j,idxo) - coo(j,icalf(2,i))
            end do
            dno = vlen(tem)
c check if angstrom
            if (dho.lt.6.6d0.and.dho.lt.dno.and.i.ne.k) then
               call hang(k,i,hand)
               cosb = hand
               if (dabs(dacos(cosb)*todeg).lt.95.d0.and.dho.lt.dismin
     &         .and.dho.lt.(5.1d0+1.89d0*dabs(cosb))) then
                  dismin = dho
                  i1 = indhb(k,dho)
                  if (i1.ne.-1) then
                     if (i1.eq.1.or.(i1.eq.2.and.ihb(1,k).ne.i)) then
                        ihb(i1,k) = i
                        hbd(i1,k) = dho
                     endif
                  endif
                  i2 = indhb(i,dho)
                  if (i2.ne.-1) then
                     if (i2.eq.1.or.(i2.eq.2.and.ihb(1,i).ne.k)) then
                        ihb(i2,i) = k
                        hbd(i2,i) = dho
                     endif
                  endif
               endif
            endif
          endif
         end do
      end do

      return
      end

      subroutine vadar(ioatms,nstrt,
     &                 icalf,ncalf,ianf,islu,nchain,iamino,isal)
c hcoord, hbond en vadar hebben geen ioadd provisie
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (numcal=50000)
      common /hbonds/ hbd(2,numcal),ihb(2,numcal)
      logical legitb
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),isal(*)

      phlast = 0.0d0

c hoe is isal geinitialiseerd -1 ?

c 0 Helix
c 1 Beta
c 2 RNA/DNA
c 3 Coil

      do i=nstrt,ncalf
        if (iamino(i).le.23) then
         ph = phi(i,icalf,ncalf,ianf,nchain)
         ps = psi(i,icalf,ncalf,islu,nchain)
         if (ihb(1,i).ne.0) then
            if (legitb(i,0).and.ph.lt.-34.d0.and.ph.gt.-118.d0
     &          .and.ps.gt.-95.d0.and.ps.lt.6.d0) then
                isal(i) = 0
c hoe zit het met i=1, hieronder
                if (i.gt.1.and.isal(i-1).ne.0.and.isal(i-2).ne.
     &              0.and.phlast.lt.-55.d0.and.phlast.gt.-90.d0) 
     &              isal(i-1) = 0
                
            elseif (legitb(i,4).and.isal(i-1).eq.0) then
                isal(i) = 0
            elseif (legitb(i,1).and.ph.lt.-45.d0.and.ph.gt.-180.d0
     &          .and.((ps.gt.95.d0.and.ps.lt.180.d0).or.
     &          (ps.lt.-170.d0.and.ps.gt.-180.d0))) then
                isal(i) = 1
                if (i.gt.1.and.isal(i-1).ne.1.and.isal(i-2).ne.
     &          1.and.phlast.lt.-100.d0) isal(i-1) = 1
            elseif (legitb(i,5).and.ph.lt.-95.d0) then
                isal(i) = 1
            else
                if (i.gt.1.and.isal(i-1).eq.0.and.ph.lt.-55.d0
     &              .and.ph.gt.-100.d0) then
                    isal(i) = 0
                else
                    isal(i) = 3
                endif
            endif
         else
            isal(i) = 3
         endif
         phlast = ph
        else
           isal(i) = 2
        endif
      end do

c     smoothing

      do i=nstrt+2,ncalf-1
         if (isal(i-1).eq.1.and.isal(i).eq.3.and.isal(i+1).eq.1)
     &       isal(i) = 1
      end do
      do i=nstrt+2,ncalf-1
         if (isal(i-1).ne.1.and.isal(i).eq.1.and.isal(i+1).ne.1)
     &       isal(i) = 3
      end do
      do i=nstrt+2,ncalf-1
         if (isal(i-1).ne.0.and.isal(i).eq.0.and.isal(i+1).ne.0)
     &       isal(i) = 3
      end do

c     filter

      icount = 1
      ilast = 0
      do i=nstrt+2,ncalf
         if (isal(i).eq.isal(i-1)) then
            icount = icount + 1
         else
            if (ilast.eq.0.and.icount.lt.4) then
               do j=i-1,i-icount,-1
                  isal(j) = 3
               end do
            endif
            if (ilast.eq.1.and.icount.lt.4) then
               do j=i-1,i-icount,-1
                  isal(j) = 3
               end do
            endif
            icount = 1
         endif
         ilast = isal(i)
      end do

c      do i=1,nchain
c         if (iamino(ianf(i)).gt.23) isal(ianf(i)) = 3
c      end do

      return
      end

      logical function legitb(idx,itype)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (numcal=50000)
      common /hbonds/ hbd(2,numcal),ihb(2,numcal)
      logical b3, b4

      b3 = .false.
      b4 = .false.
      
      legitb = .true.
      do i=1,2
         idist = abs(ihb(i,idx) - idx)
         if (ihb(i,idx).eq.0) idist = 0

         if (itype.eq.0) then
            if (idist.gt.2.and.idist.lt.6) return
         elseif (itype.eq.4) then
            if (idist.eq.3) b3 = .true.
            if (idist.eq.4) b4 = .true.
            if (b3.and.b4) return
         elseif (itype.eq.1) then
            if (idist.gt.2) return
         elseif (itype.eq.5) then
            if (idist.ge.6) return
         endif

      end do

      legitb = .false.

      return
      end

      subroutine acthd(iopt,hbfilt,coo,iconn,ianz,icalf,ncalf)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (numcal=50000)
      common /hbonds/ hbd(2,numcal),ihb(2,numcal)
      dimension tmp(3)
      dimension coo(3,*),iconn(mxcon+1,*),ianz(*)
      dimension icalf(6,*)

      do i=1,ncalf
        if (icalf(4,i).ne.0) then
         if (ihb(1,i).ne.0.and.hbd(1,i).lt.hbfilt) then

            m = icalf(4,i)
            idxc = icalf(3,ihb(1,i))
            ioxy = 0

            do j=1,iconn(1,idxc)
               n = abs(iconn(j+1,idxc))
               if (ianz(n).eq.8) ioxy = n
            end do

            do j=1,3
               tmp(j) = coo(j,m) - coo(j,ioxy)
            end do

            tmpl = vlen(tmp)

            if (dabs(tmpl-hbd(1,i)).gt.1.0d-3) then

               m = icalf(4,ihb(1,i))
               idxc = icalf(3,i)
               ioxy = 0

               do j=1,iconn(1,idxc)
                  n = abs(iconn(j+1,idxc))
                  if (ianz(n).eq.8) ioxy = n
               end do

            endif

            call hbconn(iopt,ioxy,m)

            if (ihb(2,i).ne.0.and.hbd(2,i).lt.hbfilt) then
               m = icalf(4,i)
               idxc = icalf(3,ihb(2,i))
               ioxy = 0

               do j=1,iconn(1,idxc)
                  n = abs(iconn(j+1,idxc))
                  if (ianz(n).eq.8) ioxy = n
               end do

               do j=1,3
                  tmp(j) = coo(j,m) - coo(j,ioxy)
               end do

               tmpl = vlen(tmp)

               if (dabs(tmpl-hbd(2,i)).gt.1.0d-3) then
                  m = icalf(4,ihb(2,i))
                  idxc = icalf(3,i)
                  ioxy = 0

                  do j=1,iconn(1,idxc)
                     n = abs(iconn(j+1,idxc))
                     if (ianz(n).eq.8) ioxy = n
                  end do

               endif

               call hbconn(iopt,ioxy,m)

            endif
         endif
       endif
      end do

      return
      end

      subroutine disabd(iopt,ianz,iaton,iconn)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      dimension ianz(*),iaton(*),iconn(mxcon+1,*)

      do i=1,iatoms
          if (ianz(i).eq.1) then
             if (iopt.eq.0) then
                iaton(i) = 1
             elseif (iopt.eq.1) then
                iaton(i) = 0
             elseif (iopt.eq.2) then
                if (iconn(1,i).eq.1) then
                    if (iconn(2,i).gt.0) iaton(i) = 0
                elseif (iconn(1,i).eq.2) then
                    if (iconn(2,i).gt.0.and.iconn(3,i).gt.0) 
     &                  iaton(i) = 0
               
                endif
             endif
          endif
      end do

      return
      end

      subroutine hbcond(iopt,iat1,iat2,ianz,iatclr,iaton,iconn)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxel=100)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      dimension ianz(*),iatclr(*),iaton(*),iconn(mxcon+1,*)


      if (iat1.ne.0) then
         if (iopt.eq.1) then
            iaton(iat2) = 1
            iaton(iat1) = 1
            iatclr(iat2) = 15
            iatclr(iat1) = icol(ianz(iat1))
            iconn(1,iat2) = iconn(1,iat2) + 1
            iconn(1+iconn(1,iat2),iat2) = -iat1
            iconn(1,iat1) = iconn(1,iat1) + 1
            iconn(1+iconn(1,iat1),iat1) = -iat2
         else
            iaton(iat2) = 0
            iaton(iat1) = 0
            iconn(1,iat2) = iconn(1,iat2) - 1
            iconn(1,iat1) = iconn(1,iat1) - 1
         endif
      endif
 
      return
      end

      subroutine ribbd(issdon)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      common /surf/  natorg,noscnd

c      call getnat(natoms)
c      iatoms = natoms
      noscnd = iatoms
c      natorg = iatoms

      call ribbon(0,0,1,0,0,0,0)
      call ribbon(1,0,1,0,0,0,0)
      call ribbon(2,0,1,0,0,0,0)
      call ribbon(3,0,1,0,0,0,0)

      call setarr(8,noscnd,idum)

      issdon = 1

      return
      end

      subroutine getcdh(ires,ichain,ianf,islu,nchain)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      dimension ianf(*),islu(*)

      ichain = 1

      do i=1,nchain
         if (ires.ge.ianf(i).and.ires.le.islu(i)) ichain = i
      end do

      return
      end

      subroutine ribbod(istat,iscnd,dogl,nr,iungl,
     &                  ipart,ist,ichain,tori,
     &                  coo,ianz,iaton,iatclr,iresid,iconn,
     &                  icalf,ncalf,ianf,islu,nchain,iamino,
     &                  isal,ihashb)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (alpha = 32.0d0, beta = -11.0d0)
      parameter (torad = 3.1415926536/180.0)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      parameter (numcal=50000)
      parameter (ncmx=32)
      common /strips/ qnormo(3),crpnto(3,ncmx),crnrmo(3,ncmx),
     &                numcir,nquad
      common /nfhlp/ nfree
      integer dogl
      dimension vec1(3),vec2(3),quad(3,4),cooh(3,numcal),
     &          coot(3,numcal),cooca(3,numcal),hint(3),cvec(3),rvec(3)
      dimension tori(3),qnrm(3),hintb(3,3)
      dimension coo(3,*),ianz(*),iaton(*),iatclr(*),iresid(*),
     &          iconn(mxcon+1,*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),isal(*)

      istat = 1

      if (iscnd.eq.0) then
         iheldf = 3
      elseif (iscnd.eq.1) then
         iheldf = 4
      elseif (iscnd.eq.2) then
         iheldf = 10
      endif

      if (ipart.eq.1) then
         istrt  = ist
      else
         istrt  = 1
         nfree = iatoms + 1
      endif
      iatorg = iatoms
      cilrad = 0.2d0 / 0.52917706d0
      hlxwid = 2.4d0 / 0.52917706d0
      if (iscnd.eq.0) then
          hermhx = 4.7d0 / 0.52917706d0
      else
          hermhx = 3.0d0 / 0.52917706d0
      endif
      cos1 = dcos(torad*alpha)
      sin1 = dsin(torad*alpha)
      cos2 = dcos(torad*beta)
      sin2 = dsin(torad*beta)
 
100   continue

      if (istrt.ge.ncalf) then
          if (ipart.eq.1) ist = istrt
          return
      endif

      lhelx = 0
      ifrst = 0
      ilast = 0
      ia = 0
      is = 0

      do i=istrt,ncalf
         if (isal(i).eq.iscnd.and.ifrst.eq.0) then
             ifrst = i
             do nn=1,nchain
                if (i.ge.ianf(nn).and.i.le.islu(nn)) then
                    ia = ianf(nn)
                    is = islu(nn)
                endif
             end do
         endif

         if (iscnd.eq.3.and.
     &      (isal(i).eq.3.and.isal(i+1).ne.3)) then
             do nn=1,nchain
                if (i.ge.ianf(nn).and.i.le.islu(nn)) then
                    ia = ianf(nn)
                    is = islu(nn)
                endif
             end do
             if (i.eq.ia) then
                ifrst = i
                lhelx = 2
                ilast = ifrst + 1
                istrt = ilast + 1
                goto 20
             endif
         endif

         if (iscnd.eq.3.and.
     &      (isal(i).ne.3.and.isal(i+1).eq.3)) then
c             print*,i,islu(nchain),isal(i),isal(islu(nchain))
             do nn=1,nchain
                if (i.ge.ianf(nn).and.i.le.islu(nn)) then
                    ia = ianf(nn)
                    is = islu(nn)
                endif
             end do
             if (i+1.eq.is) then
                ifrst = i
                lhelx = 2
                ilast = ifrst + 1
                istrt = ilast + 1
                goto 20
             endif
         endif

         if (iscnd.eq.3.and.
     &      ((isal(i).eq.0.and.isal(i+1).eq.1).or.
     &      (isal(i).eq.1.and.isal(i+1).eq.0))) then
             iia = 0
             iis = 0
             do nn=1,nchain
                if (i.ge.ianf(nn).and.i.le.islu(nn)) then
                    iia = ianf(nn)
                    iis = islu(nn)
                endif
             end do
             if (i.ne.iia.and.i.ne.iis.and.i+1.ne.iis) then
                ifrst = i
                lhelx = 2
                ilast = ifrst + 1
                istrt = ilast + 1
                goto 20
             endif
         endif

         if (i.eq.ncalf.and.ifrst.eq.0) then
             if (ipart.eq.1) ist = i
             return
         endif
             
         if (isal(i).ne.iscnd.and.ifrst.ne.0) goto 10

         if (isal(i).eq.iscnd) lhelx = lhelx + 1
         if (i.eq.is.and.ifrst.ne.0) goto 10

      end do
10    ilast = ifrst + lhelx - 1
      istrt = ilast + 1
 

      if (iscnd.eq.3) then
         if (ifrst.gt.ia) then
            ifrst = ifrst - 1
            lhelx = lhelx + 1
         endif
         if (ifrst+lhelx-1.lt.is) lhelx = lhelx + 1
         if (lhelx.eq.1) goto 100
      endif

20    iqopt = 1

      if (ipart.eq.1.and.iscnd.ge.2.and.ifrst.gt.islu(ichain)) then
          ist = ilast
          return
      endif

c minimim Helix/Beta length of Four
 
      if (lhelx.le.0) then
          if (ipart.eq.1) ist = ilast
          return
      endif

      if (lhelx.gt.numcal) lhelx = numcal
      if (iscnd.ne.3.and.lhelx.lt.4.and.iamino(ifrst).le.23
     &    .and.ihashb.eq.0) goto 100
      if (iscnd.ne.3.and.lhelx.lt.3.and.iamino(ifrst).le.23) goto 100
      if (iscnd.eq.2.and.lhelx.eq.1) goto 100

      ihlx = 0
      do i=1,lhelx
         call bckok(ibckok,ifrst+i-1,1)
         if (ibckok.eq.1) then
            ihlx = ihlx + 1
            do j=1,3
               cooca(j,i) = coo(j,icalf(1,ifrst+i-1))
            end do
         endif
      end do
      lhelx = ihlx
 
c mag alleen in combinatie met chain break voor of na
c anders trek een recht lijntje
c      if (lhelx.le.2) goto 100

c compute tangent vectors for ca's in 
 
      do ica = 2,lhelx-1

        if (iscnd.eq.1) then
c           do j=1,3
c first try
c              cvec(j) = cooca(j,ica) - 
c     &                 (cooca(j,ica+1) + cooca(j,ica-1))*0.5d0
c           end do
c           call vsc1(cvec,1.0d0,1.0d-4)
           do j=1,3
              vec1(j) = cooca(j,ica) - cooca(j,ica-1)
              vec2(j) = cooca(j,ica+1) - cooca(j,ica)
           end do
           call crprod(vec1,vec2,cvec)
           call vsc1(cvec,1.0d0,1.0d-4)
           
           do j=1,3
              cooh(j,ica) = cvec(j)
           end do

           do j=1,3
              cvec(j) = cooca(j,ica+1) - cooca(j,ica-1)
           end do
           call vsc1(cvec,1.0d0,1.0d-4)

           do j=1,3
              vec1(j) = cooca(j,ica) - cooca(j,ica-1)
              vec2(j) = cooca(j,ica+1) - cooca(j,ica)
           end do
           call crprod(vec1,vec2,rvec)
           call vsc1(rvec,1.0d0,1.0d-4)
 
           do j=1,3
              coot(j,ica) = (cos2*cvec(j) + sin2*rvec(j))*hermhx
           end do

        else
           do j=1,3
              cvec(j) = cooca(j,ica+1) - cooca(j,ica-1)
           end do
           call vsc1(cvec,1.0d0,1.0d-4)

           do j=1,3
              vec1(j) = cooca(j,ica) - cooca(j,ica-1)
              vec2(j) = cooca(j,ica+1) - cooca(j,ica)
           end do
           call crprod(vec1,vec2,rvec)
           call vsc1(rvec,1.0d0,1.0d-4)
 
           do j=1,3
              cooh(j,ica) = cos1*rvec(j) + sin1*cvec(j)
              coot(j,ica) = (cos2*cvec(j) + sin2*rvec(j))*hermhx
           end do
        endif
 
      end do
 
c find ca before and after chain
c Later check for chain break
 
      if (ifrst.gt.ia) ifrst = ifrst - 1
      if (ilast.lt.is) ilast = ilast + 1
 
c start initialization

      do j=1,3
           cooh(j,1) = cooh(j,2)
           coot(j,1) = cooca(j,2) - coo(j,icalf(1,ifrst))
           if (lhelx.eq.2) cooh(j,1) = -coot(j,1)
      end do
      call vsc1(coot(1,1),hermhx,1.0d-4)
 
c end initialization
 
      do j=1,3
         cooh(j,lhelx) = cooh(j,lhelx-1)
         if (iscnd.eq.1) then
            coot(j,lhelx) = coot(j,lhelx-1)
         else
            if (icalf(1,ilast).ne.0) then
               coot(j,lhelx) = coo(j,icalf(1,ilast)) - cooca(j,lhelx-1)
            endif
         endif
      end do
      call vsc1(coot(1,lhelx),hermhx,1.0d-4)

      if (iscnd.eq.1) then

c strand align and smooth
c alternating vector flip for subsequent ca's

        do ica=2,lhelx
           call timpsc(cooh(1,ica-1),cooh(1,ica),cimp)
           if (cimp.lt.0.0d0) then
              do j=1,3
                  cooh(j,ica) = -1.0d0 * cooh(j,ica)
              end do
           endif
        end do
        do ica=2,lhelx-1
           do j=1,3
               cooh(j,ica) = 
     &           cooh(j,ica-1) + cooh(j,ica) + cooh(j,ica+1)
           end do
           call vsc1(cooh(1,ica),1.0d0,1.0d-4)
        end do
        do ica=2,lhelx-1
           do j=1,3
               cooca(j,ica) = 0.5d0*(cooca(j,ica) +
     &           (cooca(j,ica-1) + cooca(j,ica+1))*0.5d0)
           end do
        end do

        if (dogl.eq.1.and.nr.eq.1) 
     &     write(iungl,'(a)') '[RIBBON] STRANDTOP'
c        if (dogl.eq.1.and.nr.ge.2) 
c     &     write(iungl,'(a)') '[RIBBON] STRANDBOTTOM'
      else
        if (dogl.eq.1.and.nr.eq.1.and.iscnd.lt.2) 
     &     write(iungl,'(a)') '[RIBBON] HELIXOUT'
c        if (dogl.eq.1.and.nr.ge.2.and.iscnd.lt.2) 
c     &     write(iungl,'(a)') '[RIBBON] HELIXIN'
        if (dogl.eq.1.and.nr.ge.2.and.iscnd.lt.2) 
     &     write(iungl,'(a)') '[HELIXIN]'
        if (dogl.eq.1.and.nr.eq.1.and.iscnd.eq.2) 
     &   write(iungl,'(a)') '[RIBBON] RNA'
        if (dogl.eq.1.and.nr.eq.1.and.iscnd.eq.3) 
     &   write(iungl,'(a)') '[RIBBON] COIL'
      endif

      if (dogl.eq.2) write(iungl,'(a)') 'mesh {'

 
c helix first segment edge

      do j=1,3
         if (dogl.ge.1) then
            if (iscnd.eq.1) then
               quad(j,1) = cooca(j,1) + hlxwid*cooh(j,1)/2.0d0
               quad(j,4) = cooca(j,1) - hlxwid*cooh(j,1)/2.0d0
            else
               quad(j,1) = cooca(j,1) + cilrad*cooh(j,1)
               quad(j,4) = cooca(j,1) - cilrad*cooh(j,1)
            endif
            hintb(j,3) = cooca(j,1)
         else
            if (nfree.le.mxnat) then
               coo(j,nfree)   = cooca(j,1)
            endif
         endif
      end do

      if (dogl.eq.0) then
         if (nfree.le.mxnat) then
            iconn(1,nfree) = 1
            iconn(2,nfree) = nfree+1
            iaton(nfree) = 1
            iatclr(nfree) = iheldf
            ianz(nfree) = 100
            iresid(nfree) = -iscnd
            nfree = nfree + 1
         else
            goto 1000
         endif
      endif
 
c helix/strand start 
 
      do iquad = 1, nquad
        t = dfloat(iquad)/dfloat(nquad)
 
        call intpol(hint,t,cooca(1,1),cooca(1,2),coot(1,1),coot(1,2))

        if (iscnd.eq.1) then
           rscal = hlxwid/2.0d0
        else
           rscal = cilrad + (hlxwid/2.0d0 - cilrad)*
     &          0.5d0*(-dcos(180.0d0*torad*t) + 1.0d0)
        endif

        do j=1,3
           if (dogl.ge.1) then
              vec1(j) = rscal*cooh(j,1)
              quad(j,2) = hint(j) + vec1(j)
              quad(j,3) = hint(j) - vec1(j)
              hintb(j,1) = hintb(j,2)
              hintb(j,2) = hintb(j,3)
              hintb(j,3) = hint(j)
           else
              if (nfree.le.mxnat) then
                 coo(j,nfree)   = hint(j)
              else
                 goto 1000
              endif
           endif
        end do

        if (dogl.ge.1) then

c compute normal

           call vsc1(vec1,1.0d0,1.0d-8)
           do j=1,3
              vec2(j) = (quad(j,1) + quad(j,4))*0.5d0 -
     &                  (quad(j,2) + quad(j,3))*0.5d0
           end do
           call vsc1(vec2,1.0d0,1.0d-8)
           call crprod(vec1,vec2,qnrm)
           call vsc1(qnrm,1.0d0,1.0d-8)

c put out quadrilateral

           if (iscnd.ge.2) then
              if (iquad.ge.2) then
                  call wrcoil(hintb(1,1),hintb(1,2),
     &                        hintb(1,3),iungl,tori,iqopt,dogl)
                  iqopt = 2
              endif
           else
              call wrquad(quad,qnrm,iungl,tori,iqopt,nr,dogl)
              iqopt = 2
           endif
           do j=1,3
              quad(j,1) = quad(j,2)
              quad(j,4) = quad(j,3)
           end do

        else

           if (nfree.le.mxnat) then
              iconn(1,nfree) = 2
              iconn(2,nfree) = nfree+1
              iconn(3,nfree) = nfree-1
              iaton(nfree) = 1
              iatclr(nfree) = iheldf
              ianz(nfree) = 100
              iresid(nfree) = -iscnd
              nfree = nfree + 1
           else
              goto 1000
           endif

        endif

      end do
 
c helix/strand main 
 
      do ica = 2, lhelx - 2
        do iquad = 1, nquad

          t = dfloat(iquad)/dfloat(nquad)
 
          call intpol(hint,t,cooca(1,ica),cooca(1,ica+1),
     &                coot(1,ica),coot(1,ica+1))

          rscal = hlxwid/2.0d0
          do j=1,3
             vec1(j) = (1.0d0-t)*cooh(j,ica) + t*cooh(j,ica+1)
          end do

          call vsc1(vec1,rscal,1.0d-4)

          do j=1,3
              if (dogl.ge.1) then
                 quad(j,2) = hint(j) + vec1(j)
                 quad(j,3) = hint(j) - vec1(j)
                 hintb(j,1) = hintb(j,2)
                 hintb(j,2) = hintb(j,3)
                 hintb(j,3) = hint(j)
              else
                 if (nfree.le.mxnat) then
                    coo(j,nfree)   = hint(j)
                 else
                    goto 1000
                 endif
              endif
          end do

          if (dogl.ge.1) then

c compute normal

             call vsc1(vec1,1.0d0,1.0d-8)
             do j=1,3
                vec2(j) = (quad(j,1) + quad(j,4))*0.5d0 -
     &                    (quad(j,2) + quad(j,3))*0.5d0
             end do
             call vsc1(vec2,1.0d0,1.0d-8)
             call crprod(vec1,vec2,qnrm)
             call vsc1(qnrm,1.0d0,1.0d-8)

c put out quadrilateral

             if (iscnd.ge.2) then
                call wrcoil(hintb(1,1),hintb(1,2),
     &                      hintb(1,3),iungl,tori,iqopt,dogl)
             else
                call wrquad(quad,qnrm,iungl,tori,iqopt,nr,dogl)
             endif
             do j=1,3
                quad(j,1) = quad(j,2)
                quad(j,4) = quad(j,3)
             end do

          else

             if (nfree.le.mxnat) then
                iconn(1,nfree) = 2
                iconn(2,nfree) = nfree+1
                iconn(3,nfree) = nfree-1
                iaton(nfree) = 1
                iatclr(nfree) = iheldf
                ianz(nfree) = 100
                iresid(nfree) = -iscnd
                nfree = nfree + 1
             else
                goto 1000
             endif

          endif


        end do
      end do
 
c helix/strand end
 
      if (lhelx.eq.2 ) then
        if (dogl.ge.1) then
           goto 200
        else
           if (nfree-1.le.mxnat) then
              iconn(1,nfree-1) = 1
              iconn(2,nfree-1) = iconn(3,nfree-1)
           else
              goto 1000
           endif
           goto 100
        endif
      endif

      do iquad = 1, nquad
        t = dfloat(iquad)/dfloat(nquad)
 
        call intpol(hint,t,cooca(1,lhelx-1),cooca(1,lhelx),
     &              coot(1,lhelx-1),coot(1,lhelx))

        if (iscnd.eq.1) then
           rscal = 1.8d0*hlxwid*(1.01d0-t)/2.0d0
        else
           rscal = cilrad + (hlxwid/2.0d0 - cilrad)*
     &          0.5d0*(dcos(180.0d0*torad*t) + 1.0d0)
        endif
        do j=1,3
           if (dogl.ge.1) then
              vec1(j) = rscal*cooh(j,lhelx)
              quad(j,2) = hint(j) + vec1(j)
              quad(j,3) = hint(j) - vec1(j)
              hintb(j,1) = hintb(j,2)
              hintb(j,2) = hintb(j,3)
              hintb(j,3) = hint(j)
           else
              if (nfree.le.mxnat) then
                 coo(j,nfree)   = hint(j)
              else
                 goto 1000
              endif
           endif
        end do

        if (dogl.ge.1) then

c compute normal

           call vsc1(vec1,1.0d0,1.0d-8)
           do j=1,3
              vec2(j) = (quad(j,1) + quad(j,4))*0.5d0 -
     &                  (quad(j,2) + quad(j,3))*0.5d0
           end do
           call vsc1(vec2,1.0d0,1.0d-8)
           call crprod(vec1,vec2,qnrm)
           call vsc1(qnrm,1.0d0,1.0d-8)


c put out quadrilateral

           if (iscnd.ge.2) then
              call wrcoil(hintb(1,1),hintb(1,2),
     &                      hintb(1,3),iungl,tori,iqopt,dogl)
           else
              call wrquad(quad,qnrm,iungl,tori,iqopt,nr,dogl)
           endif
           do j=1,3
              quad(j,1) = quad(j,2)
              quad(j,4) = quad(j,3)
           end do

        else

           if (nfree.le.mxnat) then
              iconn(1,nfree) = 2
              iconn(2,nfree) = nfree-1
              iconn(3,nfree) = nfree+1
              iaton(nfree) = 1
              iatclr(nfree) = iheldf
              ianz(nfree) = 100
              iresid(nfree) = -iscnd
              nfree = nfree + 1
           else
              goto 1000
           endif

        endif
 

      end do
 
200   continue

      if (dogl.ge.1) then

         if (iscnd.ge.2) then

            do j=1,3
              hintb(j,1) = hintb(j,2)
              hintb(j,2) = hintb(j,3)
c              hintb(j,3) = cooca(j,lhelx+1)
            end do

            do j=1,3
              hintb(j,3) = hintb(j,2) + (hintb(j,2) - hintb(j,1))
            end do

            call wrcoil(hintb(1,1),hintb(1,2),
     &                  hintb(1,3),iungl,tori,iqopt,dogl)
         else
            call wrquad(quad,qnrm,iungl,tori,3,nr,dogl)
         endif

      endif

      if (dogl.eq.2) then

         if (iscnd.eq.1) then
           if (nr.eq.1)
     &        write(iungl,'(a)') 'texture { STRANDTOP }'
           if (nr.ge.2)
     &        write(iungl,'(a)') 'texture { STRANDBOTTOM }'
         else
           if (nr.eq.1.and.iscnd.lt.2)
     &        write(iungl,'(a)') 'texture { HELIXOUT }'
           if (nr.ge.2.and.iscnd.lt.2)
     &        write(iungl,'(a)') 'texture { HELIXIN }'
           if (nr.eq.1.and.iscnd.eq.2)
     &      write(iungl,'(a)') 'texture { RNA }'
           if (nr.eq.1.and.iscnd.eq.3)
     &      write(iungl,'(a)') 'texture { COIL }'
         endif

         write(iungl,'(a)') '}'

      endif

      if (ipart.eq.1) then

         if (iscnd.ge.2) then
            if (ilast.gt.islu(ichain)) then
               ist = ilast
               return
            endif
         else
            ist = ilast
            return
         endif
      endif

      if (dogl.ge.1) goto 100


      if (nfree-1.le.mxnat) then
          iatoms = nfree-1
          iconn(1,nfree-1) = 1
      else
          goto 1000
      endif

      goto 100

1000  istat = 0
      iatoms = iatorg
      istart = istrt
      return
 
      end

      subroutine intpol(cnew,g,c1,c2,v1,v2)
      implicit double precision (a-h,o-z)
      dimension cnew(3),c1(3),c2(3),v1(3),v2(3)
 
      do i = 1, 3
        cnew(i) = c1(i)*(2.0*g*g*g - 3.0*g*g + 1.0)
     &          + c2(i)*(-2.0*g*g*g + 3.0*g*g)
     &          + v1(i)*(g*g*g - 2.0*g*g + g) + v2(i)*(g*g*g - g*g)
      end do
 
      return
      end

      subroutine wrquad(quad,qnrm,iungl,tori,iopt,nr,dogl)
      implicit double precision (a-h,o-z)
      parameter (ncmx=32)
      common /strips/ qnormo(3),crpnto(3,ncmx),crnrmo(3,ncmx),
     &                numcir,nquad
      integer dogl
      dimension quad(3,4),v1(3),v2(3),tori(3)
      dimension quadt(3,4),qnew(3),qnrm(3)

 
      ribbw = 0.4d0
      qs = 1.0

      call vsc1(qnrm,1.0d0,1.0d-4)

      if (iopt.eq.1.and.nr.eq.4) then

c in between side quad of begin

         do j=1,3
            quadt(j,1) = quad(j,4) + qnrm(j)*ribbw
            quadt(j,2) = quad(j,4) - qnrm(j)*ribbw
            quadt(j,3) = quad(j,1) - qnrm(j)*ribbw
            quadt(j,4) = quad(j,1) + qnrm(j)*ribbw
         end do
         do j=1,3
            qnew(j) = (quad(j,2) + quad(j,3))*0.5d0 -
     &                (quad(j,1) + quad(j,4))*0.5d0
         end do
         call vsc1(qnew,1.0d0,1.0d-4)

         if (dogl.eq.2) write(iungl,'(a)') 'smooth_triangle {'
         call wrvert(qnew,qs,quadt(1,1),tori,iungl,dogl)
         call wrvert(qnew,qs,quadt(1,4),tori,iungl,dogl)
         if (dogl.eq.2) then
            call wrvert(qnew,qs,quadt(1,3),tori,iungl,dogl)
            write(iungl,'(a)') '}'
            write(iungl,'(a)') 'smooth_triangle {'
         endif

         call wrvert(qnew,qs,quadt(1,3),tori,iungl,dogl)
         call wrvert(qnew,qs,quadt(1,2),tori,iungl,dogl)
         if (dogl.eq.2) then
            call wrvert(qnew,qs,quadt(1,1),tori,iungl,dogl)
            write(iungl,'(a)') '}'
         endif

      endif

      if (iopt.eq.1) then
         do j=1,3
            qnew(j) = qnrm(j)
         end do
      elseif (iopt.eq.2) then
         do j=1,3
            qnew(j) = qnrm(j) + qnormo(j)
         end do
      elseif (iopt.eq.3) then
         do j=1,3
            qnew(j) = qnormo(j)
         end do
      endif

      call vsc1(qnew,1.0d0,1.0d-4)

c top and bottom strip

      if (nr.eq.1.or.nr.eq.2) then
         do i=1,4
            do j=1,3
               if (nr.eq.1) then
                  qs = 1.0
                  quadt(j,i) = quad(j,i) + qnew(j)*ribbw
               elseif (nr.eq.2) then
                  qs = -1.0
                  quadt(j,5-i) = quad(j,i) - qnew(j)*ribbw
               endif
            end do
         end do
      endif

c left side strip

      if (nr.eq.3) then
         qs = -1.0
         do j=1,3
            quadt(j,1) = quad(j,1) - qnew(j)*ribbw
            quadt(j,2) = quad(j,2) - qnew(j)*ribbw
         end do
         do j=1,3
            quadt(j,4) = quad(j,1) + qnew(j)*ribbw
            quadt(j,3) = quad(j,2) + qnew(j)*ribbw
         end do
         do j=1,3
            qnew(j) = quad(j,4) - quad(j,1)
         end do
         call vsc1(qnew,1.0d0,1.0d-4)
      endif

c right side strip

      if (nr.eq.4) then
         qs = 1.0
         do j=1,3
            quadt(j,1) = quad(j,4) + qnew(j)*ribbw
            quadt(j,2) = quad(j,3) + qnew(j)*ribbw
         end do
         do j=1,3
            quadt(j,4) = quad(j,4) - qnew(j)*ribbw
            quadt(j,3) = quad(j,3) - qnew(j)*ribbw
         end do
         do j=1,3
            qnew(j) = quad(j,4) - quad(j,1)
         end do
         call vsc1(qnew,1.0d0,1.0d-4)
      endif


      if (iopt.eq.2.or.iopt.eq.3) then

         if (dogl.eq.2) then
            call wrvert(qnew,qs,quadt(1,4),tori,iungl,dogl)
            write(iungl,'(a)') '}'
            write(iungl,'(a)') 'smooth_triangle {'
         endif

         call wrvert(qnew,qs,quadt(1,4),tori,iungl,dogl)
         call wrvert(qnew,qs,quadt(1,1),tori,iungl,dogl)
         if (dogl.eq.2) then
            call wrvert(crnrmo(1,1),qs,crpnto(1,1),tori,iungl,dogl)
            write(iungl,'(a)') '}'
         endif
      endif

      if (iopt.eq.1.or.iopt.eq.2) then
         if (dogl.eq.2) write(iungl,'(a)') 'smooth_triangle {'

         call wrvert(qnew,qs,quadt(1,1),tori,iungl,dogl)
         call wrvert(qnew,qs,quadt(1,4),tori,iungl,dogl)

         if (dogl.eq.2) then
            do i=1,4
               do j=1,3
                  crpnto(j,1) = quadt(j,1)
                  crnrmo(j,1) = qnew(j)
               end do
            end do
         endif

      endif

      if (iopt.eq.3.and.nr.eq.4) then

c in between side quad of begin

         qs = -1.0d0;

         do j=1,3
            quadt(j,1) = quad(j,1) + qnrm(j)*ribbw
            quadt(j,4) = quad(j,4) + qnrm(j)*ribbw
            quadt(j,2) = quad(j,1) - qnrm(j)*ribbw
            quadt(j,3) = quad(j,4) - qnrm(j)*ribbw
         end do
         do j=1,3
            v1(j) = quadt(j,4) - quadt(j,1)
            v2(j) = quadt(j,2) - quadt(j,1)
         end do
         call crprod(v1,v2,qnew)
         call vsc1(qnew,1.0d0,1.0d-8)

         if (dogl.eq.2) write(iungl,'(a)') 'smooth_triangle {'
         call wrvert(qnew,qs,quadt(1,1),tori,iungl,dogl)
         call wrvert(qnew,qs,quadt(1,4),tori,iungl,dogl)
         if (dogl.eq.2) then
            call wrvert(qnew,qs,quadt(1,3),tori,iungl,dogl)
            write(iungl,'(a)') '}'
            write(iungl,'(a)') 'smooth_triangle {'
         endif

         call wrvert(qnew,qs,quadt(1,3),tori,iungl,dogl)
         call wrvert(qnew,qs,quadt(1,2),tori,iungl,dogl)
         if (dogl.eq.2) then
            call wrvert(qnew,qs,quadt(1,1),tori,iungl,dogl)
            write(iungl,'(a)') '}'
         endif

      endif

      do j=1,3
         qnormo(j) = qnrm(j)
      end do


      return
      end

      subroutine wrcoil(c1,c2,c3,iungl,tori,iopt,dogl)
      implicit double precision (a-h,o-z)
      parameter (ncmx=32)
      common /strips/ qnormo(3),crpnto(3,ncmx),crnrmo(3,ncmx),
     &                numcir,nquad
      integer dogl
      dimension c1(3),c2(3),c3(3),v1(3),v2(3),vt1(3),vt2(3),vt3(3)
      dimension tori(3),cirpnt(3,ncmx),cirnrm(3,ncmx)
      dimension ids(2)
 
      pi = 3.141592654d0
      twopi = 2*pi
      coilw = 0.5d0
      ids(1) = 1
      ids(2) = -1


c iopt = 1, first segment, only c2 and c3 defined

      do j=1,3
         vt1(j) = c1(j) - c2(j)
         vt2(j) = c3(j) - c2(j)
      end do
      call crprod(vt1,vt2,v2)
      call vsc1(v2,1.0d0,1.0d-8)

      call vsc1(vt1,1.0d0,1.0d-8)
      call vsc1(vt2,1.0d0,1.0d-8)

      do j=1,3
         v1(j) = vt1(j) + vt2(j)
      end do
      call vsc1(v1,1.0d0,1.0d-8)

      if (vlen(v1).lt.1.0d-8.and.iopt.ne.1) then
         do j=1,3
            v2(j) = crnrmo(j,numcir)
         end do
         call crprod(vt1,v2,v1)
      endif

      call crprod(v1,v2,vt3)

      
      do i=1,numcir
         s = twopi * dble(i)/dble(numcir)
         do j=1,3
            cirnrm(j,i) = dsin(s)*v1(j) + dcos(s)*v2(j)
            cirpnt(j,i) = c2(j) + cirnrm(j,i)*coilw
         end do
      end do

      if (iopt.eq.1) then

         call crprod(v2,vt1,v1)
         call vsc1(v1,1.0d0,1.0d-8)

         do i=1,numcir
            s = twopi * dble(i)/dble(numcir)
            do j=1,3
               crnrmo(j,i) = dsin(s)*v1(j) + dcos(s)*v2(j)
               crpnto(j,i) = c1(j) + crnrmo(j,i)*coilw
            end do
         end do

      endif

      idir = 1
      dismin = 10000.0d0
      do l=1,2
         do iof=0,numcir-1
            dsq = 0.0d0
            do i=1,numcir
               i1 = i
               i2 = iof + ids(l)*i
               if (i2.gt.numcir) then
                  i2 = i2 - numcir
               elseif (i2.lt.1) then
                  i2 = i2 + numcir
               endif
               dsq = dsq + dist2(crpnto(1,i1),cirpnt(1,i2))
            end do
            if (dsq.lt.dismin) then
               dismin = dsq
               ioff = iof
               idir = ids(l)
            endif
         end do
      end do

c write quads

      do i=1,numcir
         i1 = i
         i2 = i + 1
         if (i2.gt.numcir) i2 = i2 - numcir
         l1 = ioff + idir*i
         if (l1.gt.numcir) then
            l1 = l1 - numcir
         elseif (l1.lt.1) then
            l1 = l1 + numcir
         endif
         l2 = ioff + idir*(i+1)
         if (l2.gt.numcir) then
            l2 = l2 - numcir
         elseif (l2.lt.1) then
            l2 = l2 + numcir
         endif

c dogl.eq.2 (povray), split quad up in to triangles

         if (dogl.eq.2) write(iungl,'(a)') 'smooth_triangle {'
         call wrvert(crnrmo(1,i1),1.0d0,crpnto(1,i1),tori,iungl,dogl)
         call wrvert(crnrmo(1,i2),1.0d0,crpnto(1,i2),tori,iungl,dogl)
         if (dogl.eq.2) then
            call wrvert(cirnrm(1,l2),1.0d0,cirpnt(1,l2),tori,iungl,dogl)
            write(iungl,'(a)') '}'
            write(iungl,'(a)') 'smooth_triangle {'
         endif

         call wrvert(cirnrm(1,l2),1.0d0,cirpnt(1,l2),tori,iungl,dogl)
         call wrvert(cirnrm(1,l1),1.0d0,cirpnt(1,l1),tori,iungl,dogl)
         if (dogl.eq.2) then
            call wrvert(crnrmo(1,i1),1.0d0,crpnto(1,i1),tori,iungl,dogl)
            write(iungl,'(a)') '}'
         endif

      end do

      do i=1,numcir
         do j=1,3
            crnrmo(j,i) = cirnrm(j,i)
            crpnto(j,i) = cirpnt(j,i)
         end do
      end do

      return
      end

      subroutine wrvert(qnrm,qs,qvert,tori,iungl,dogl)
      implicit double precision (a-h,o-z)
      integer dogl
      common /setogl/ idirogl
      dimension qnrm(3),qvert(3),tori(3)

      toang = 0.52917706d0

      if (idirogl.eq.1) then
c        call ognrm(qnrm(1),qnrm(2),qnrm(3))
        call ognrm(qnrm(1)*qs,qnrm(2)*qs,qnrm(3)*qs)
        call ogvrt(qvert(1),qvert(2),qvert(3))
      else
        if (dogl.eq.1) then
         write(iungl,'(3(f9.6,1x))') qnrm(1)*qs,qnrm(2)*qs,qnrm(3)*qs
         write(iungl,'(3(f12.6,1x))') (qvert(1)-tori(1))*toang,
     &     (qvert(2)-tori(2))*toang,(qvert(3)-tori(3))*toang
        elseif (dogl.eq.2) then
         write(iungl,'(2(a,3(f12.6,a)))')
     &   '<',(qvert(1)-tori(1))*toang,',',(qvert(2)-tori(2))*toang,
     &   ',',(qvert(3)-tori(3))*toang,'>, ',
     &   '<',qnrm(1),',',qnrm(2),',',qnrm(3),'>'
        endif
      endif

      return
      end

      subroutine ribgl(ianf,nchain,ncalf,iatoms)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /setogl/ idirogl
      dimension ianf(*)

      idirogl = 1
      iun = 47
      open(unit=iun,form='formatted',status='scratch')

      call ribgll(iun,ianf,nchain,ncalf,iatoms)

      close(iun)
      idirogl = 0

      return
      end

      subroutine ribgll(iun,ianf,nchain,ncalf,iatoms)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /nfhlp/ nfree
      common /setogl/ idirogl
      dimension ianf(*)

      nfree = iatoms + 1

c helix

      istart = 1

      if (idirogl.eq.1) call ribpnt(1,0)
      
      do while (istart.lt.ncalf)

         if (idirogl.eq.1) then
            call ogribb(0)
            call sribcol(2)
         endif

         ist = istart
         call ribbon(0,1,1,iun,1,ist,0)
         if (idirogl.eq.1) then
            call sribcol(3)
            call setcll
         endif
         ist = istart
         call ribbon(0,1,2,iun,1,ist,0)
         ist = istart
         call ribbon(0,1,3,iun,1,ist,0)
         ist = istart
         call ribbon(0,1,4,iun,1,ist,0)

         istart = ist

         call getchn(ist,ichain)
         if (idirogl.eq.1) call ogendd(ichain)

      end do

      if (idirogl.eq.1) call ribpnt(0,0)

c strand

      istart = 1

      if (idirogl.eq.1) call ribpnt(1,1)

      do while (istart.lt.ncalf)

         if (idirogl.eq.1) then
            call ogribb(1)
            call sribcol(1)
         endif

         ist = istart
         call ribbon(1,1,1,iun,1,ist,0)
         ist = istart
         call ribbon(1,1,2,iun,1,ist,0)
         ist = istart
         call ribbon(1,1,3,iun,1,ist,0)
         ist = istart
         call ribbon(1,1,4,iun,1,ist,0)

         istart = ist

         call getchn(ist,ichain)
         if (idirogl.eq.1) call ogendd(ichain)

      end do

      if (idirogl.eq.1) then
          call ribpnt(0,1)
          call ribpnt(1,2)
      endif

      do i=1,nchain
          ist = ianf(i)
c rna
          if (idirogl.eq.1) then
             call ogribb(2)
             call sribcol(4)
          endif
          call ribbon(2,1,1,iun,1,ist,i)
          if (idirogl.eq.1) then
             call ogendd(i)
          endif
      end do

      if (idirogl.eq.1) then
          call ribpnt(0,2)
          call ribpnt(1,3)
      endif

      do i=1,nchain
          ist = ianf(i)
c coil
          if (idirogl.eq.1) then
             call ogribb(3)
             call sribcol(5)
          endif
          call ribbon(3,1,1,iun,1,ist,i)
          if (idirogl.eq.1) then
             call ogendd(i)
          endif

      end do

      if (idirogl.eq.1) call ribpnt(0,3)
  
      close(iun)

      return
      end

