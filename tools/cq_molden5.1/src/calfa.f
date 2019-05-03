      subroutine calfd(istat,istpdb,iaddh,ioatms,nstrt,ioadd,
     &                 ianz,iconn,ityp,
     &                 icalf,ncalf,ianf,islu,nchain,iamino,
     &                 isal,irsnr,ihashb,rphi,rpsi)

      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      parameter (mxsym=103)
      parameter (mxhsym=64)
      logical chkpdb
      common /surf/ natorg,noscnd
      logical cn,cco,ccb,cco1,chkhs
      common /types/ iff
      integer*2 ityp
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension ianz(*),iconn(mxcon+1,*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),isal(*),irsnr(*),
     &          ityp(*),rphi(*),rpsi(*)

c icalf(1,*) is Calfa
c icalf(2,*) is N connected to Calfa
c icalf(3,*) is C=O connected to Calfa

      istat = 1

      if (istpdb.eq.1) goto 1000

      if (ioadd.eq.0) ncalf = 0

      do i=ioatms+1,iatoms

         cn  = .false.
         cco = .false.
         ccb = .false.
         cco1 = .false.

         if (ianz(i).eq.6.and..not.(iff.eq.7.and.ityp(i).lt.0)) then

            icnso = 0

            do j=1,iconn(1,i)

               jj = abs(iconn(j+1,i))

               if (jj.gt.ioatms) then
                  if (ianz(jj).eq.6 ) icnso = icnso + 1
                  if (ianz(jj).eq.7 ) icnso = icnso + 1
                  if (ianz(jj).eq.8 ) icnso = icnso + 1
                  if (ianz(jj).eq.16) icnso = icnso + 1
               endif

            end do

           if (icnso.eq.2.or.icnso.eq.3) then
c is c connected to nitrogen ?
c
            do j=1,iconn(1,i)
               jj = abs(iconn(j+1,i))
               if (jj.gt.ioatms) then
                  if (ianz(jj).eq.7) then
                      cn = .true.
                      icn = jj
                  endif
               endif
            end do

            if (.not.cn) goto 100
c
c is c connected to carbonyl ?
c
            do j=1,iconn(1,i)
               l = abs(iconn(j+1,i))
               cco1 = .false.
               if (ianz(l).eq.6.and.l.gt.ioatms) then
                  if (iconn(1,l).eq.3.or.iconn(1,l).eq.2) then

                     do k=1,iconn(1,l)
                        kk = abs(iconn(k+1,l))
                        if (ianz(kk).eq.8.and.kk.gt.ioatms) then
                           cco1 = .true.
                        end if
                     end do

                     if (cco1) then
                       if (iconn(1,l).eq.3) then
                          do k=1,iconn(1,l)
                             kk = abs(iconn(k+1,l))
                             if (ianz(kk).eq.7.and.kk.gt.ioatms) then
                                cco = .true.
                                icco = l
                             end if
                          end do
                       else
                          cco = .true.
                          icco = l
                       endif
                     endif

                  endif 
               endif
            end do

            if (.not.cco) goto 100
c
c is c connected to carbon other than those found?
c
            if (icnso.eq.3) then

               do j=1,iconn(1,i)
                  l = abs(iconn(j+1,i))
                  if (l.gt.ioatms) then
                    if (ianz(l).eq.6.and.l.ne.icco) ccb = .true.
                  endif
               end do
            endif

            if (icnso.eq.2) ccb = .true.

            if (ccb) then
               ncalf = ncalf + 1
               icalf(1,ncalf) = i
               icalf(2,ncalf) = icn
               icalf(3,ncalf) = icco
            endif

          endif
         endif
100      continue
      end do

c     round up terminal C alfa

      if (ioadd.eq.1) then
         nstr = nstrt
      else 
         nchain = 0
         nstr = 1
      endif

      ncalft = ncalf

5     continue 

      ns = nstr

      do i=ns,ncalft
         do j=1,iconn(1,icalf(3,i))
            l = abs(iconn(j+1,icalf(3,i)))

c l = N connected to CO , but not yet in list

            if (ianz(l).eq.7.and.l.gt.ioatms) then

               call fndcal(l,2,idum,istat,icalf,ncalft)

               if (istat.eq.0) then

                  do k=1,iconn(1,l)
                    m = abs(iconn(k+1,l))

                    if (ianz(m).eq.6.and.m.ne.icalf(3,i).and.
     &                  m.gt.ioatms) then

c m = C(alfa?) connected to N and is not previous CO

                      do n=1,iconn(1,m)
                        mm = abs(iconn(n+1,m))
                        if (ianz(mm).eq.6.and.mm.gt.ioatms) then

                        if (iconn(1,mm).eq.3) then

                          ioo = 0
                          do in=1,3
                             if (ianz(abs(iconn(in+1,mm))).eq.8) 
     &                           ioo = ioo + 1
                          end do

c Calfa is connected to a CO2

                          if (ioo.eq.2) then
                             ncalf = ncalf + 1
                             icalf(1,ncalf) =  m
                             icalf(2,ncalf) =  l
                             icalf(3,ncalf) = mm
                             nstr = nstr + 1
                             goto 5
                          endif

                        endif

c Calfa is connected to a CO

                        if (iconn(1,mm).eq.2.and.(ianz(abs(iconn(2,mm)))
     &                  .eq.8.or.ianz(abs(iconn(3,mm))).eq.8)) then
                             ncalf = ncalf + 1
                             icalf(1,ncalf) =  m
                             icalf(2,ncalf) =  l
                             icalf(3,ncalf) = mm
                             nstr = nstr + 1
                             goto 5
                        endif

                        endif
                      end do
                    endif
                  end do
               endif
            end if
         end do
         nstr = nstr + 1
      end do

      if (ioadd.eq.1) then
         nstr = nstrt
      else 
         nchain = 0
         nstr = 1
      endif

10    continue

      do i=nstr,ncalf

         ico = 0

         do j=1,iconn(1,icalf(2,i))

c jj = N alfa not connected to a CO alfa of previous residue
c this an N-term

            jj = abs(iconn(j+1,icalf(2,i)))
            if (ianz(jj).eq.6.and.jj.gt.ioatms) then
               call fndcal(jj,3,indx,istat,icalf,ncalf)
               if (istat.eq.1) ico = 1
            endif

         end do

         if (ico.eq.0) then

             do k=1,3
                itemp = icalf(k,nstr)
                icalf(k,nstr) = icalf(k,i)
                icalf(k,i) = itemp
             end do

             nchain = nchain + 1
             ianf(nchain) = nstr
             goto 20

         endif

      end do

      goto 30

20    continue

          do j=1,iconn(1,icalf(3,nstr))
             jj = abs(iconn(j+1,icalf(3,nstr)))

c jj = N alfa connected to an CO alfa

             if (ianz(jj).eq.7.and.jj.gt.ioatms) then

                call fndcal(jj,2,indx,istat,icalf,ncalf)
                if (istat.eq.1) then

c                  go along chain untill C-Term found

                   nstr = nstr + 1

                   do k=1,3
                      itemp = icalf(k,nstr)
                      icalf(k,nstr) = icalf(k,indx)
                      icalf(k,indx) = itemp
                   end do

                   goto 20

                endif

             endif

          end do

c found C-Term; C=O, not connected to N 

          islu(nchain) = nstr
          nstr = nstr + 1

          goto 10
30    continue

      if (nchain.eq.0.and.ncalf.gt.0) then
          nchain = 1
          ianf(nchain) = 1
          islu(nchain) = ncalf
      endif

c determine what residue, doesnt work for nucls

      call samino(istat)

      call parsfn('Helix',5,1)
      call parsfn('Beta',4,1)
      call parsfn('RNA/DNA',7,1)
      call parsfn('Coil',4,1)
      call parsfn('HET',3,1)

1000  continue

c hcoord, hbond en vadar hebben geen ioadd provisie

      call hcoord(ioatms,nstrt)
      call hbond(ioatms,nstrt)
      if (ihashb.eq.0) then
          call vadar(ioatms,nstrt,icalf,ncalf,ianf,islu,nchain,iamino,
     &               isal)
          call inferr('Secondary Structure calculated',0)
      else
          isold = 10000
          do i=1,ncalf
             if (iamino(i).gt.23) isal(i) = 2
             if (isal(i).ne.isold) then
                 if (i.gt.1) then
                    if (lng.le.2) isal(i-1) = 3
                 endif
                 if (i.gt.2) then
                    if (lng.eq.2) isal(i-2) = 3
                 endif
                 lng = 1
             else
                 lng = lng + 1
             endif
             isold = isal(i)
          end do
      endif
      
      if (iff.eq.0) then
         do i=1,nchain
            do j=ianf(i),islu(i)
               call getpdb(j,ipdb,ihpdb)
               if (.not.chkpdb(ipdb,iamino(j),j,irsnr)) then
                   print*,'incomplete residue: internal no. ',j
               endif
               if (iaddh.eq.1.and..not.chkhs(ihpdb)) then
                  if (j.eq.ianf(i)) then
                     call addhs(j,iamino(j),ipdb,ihpdb,1)
                  else
                     call addhs(j,iamino(j),ipdb,ihpdb,0)
                  endif
               endif
               call typeit(ipdb,iamino(j),ihpdb,1)
            end do
         end do
      else
         do i=1,nchain
            do j=ianf(i),islu(i)
               call getpdb(j,ipdb,ihpdb)
               if (.not.chkpdb(ipdb,iamino(j),j,irsnr)) then
                   print*,'incomplete residue: internal no.',j
               endif
            end do
         end do
      endif

      call prtcal(rphi,rpsi,icalf,ncalf,ianf,islu,nchain,iamino,isal)

      noscnd = iatoms

      return
      end

      subroutine setcdd(iopt,nchain,ianf,islu,irsnr,achain)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      character*1 achain
      dimension ianf(*),islu(*),irsnr(*),achain(*)

      if (iopt.eq.1) then

           do k=1,nchain
              do j=ianf(k),islu(k)
                 achain(j) = char(64+k)
                 irsnr(j) = j
              end do
           end do

      else 

           do k=1,nchain
              do j=ianf(k),islu(k)
                 achain(j) = char(64+k)
                 jj = j - ianf(k) + 1
                 irsnr(j) = jj
              end do
           end do

      endif

      return
      end

