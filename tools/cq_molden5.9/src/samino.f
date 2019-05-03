      subroutine samind(istat,ianz,iresid,iconn,ipdbt,
     &                  icalf,ncalf,iamino)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      integer*2 ipdbt
      dimension icd(3),icdh(3),icdhn(3)
      dimension ianz(*),iresid(*),iconn(mxcon+1,*),ipdbt(*)
      dimension icalf(6,*),iamino(*)
  
      istat = 0

      do i=1,iatoms
         iresid(i) = -4
      end do

      do k=1,ncalf

      icbet = 0
      icalp = icalf(1,k)
      call setvar(iresid,mxnat,icalf(1,k),k)
      call setvr2(ipdbt,mxnat,icalf(1,k),2)
      call setvr2(ipdbt,mxnat,icalf(2,k),1)
      call setvr2(ipdbt,mxnat,icalf(3,k),3)

      ich = 0
      do i=1,iconn(1,icalf(1,k))
         ii = abs(iconn(i+1,icalf(1,k)))
         call setvar(iresid,mxnat,ii,k)
         if (ianz(ii).eq.6.and.ii.ne.icalf(3,k)) icbet = ii
         if (ianz(ii).eq.1) then
            ich = ich + 1
            call setvr2(ipdbt,mxnat,ii,3 + ich)
         endif
      end do
      inh = 0
      do i=1,iconn(1,icalf(2,k))
         ii = abs(iconn(i+1,icalf(2,k)))
         if (ianz(ii).eq.1) then
            inh = inh + 1
            call setvar(iresid,mxnat,ii,k)
            call setvr2(ipdbt,mxnat,ii,inh)
         endif
      end do
      io = 0
      do i=1,iconn(1,icalf(3,k))
         ii = abs(iconn(i+1,icalf(3,k)))
         if (ianz(ii).eq.8) then
            io = io + 1
            call setvar(iresid,mxnat,ii,k)
            if (io.eq.1) then
               call setvr2(ipdbt,mxnat,ii,4)
            else
               call setvr2(ipdbt,mxnat,ii,38)
            endif
         endif
      end do
      iamino(k) = 0
      if (icbet.eq.0) then
c Glycine
         iamino(k) = 1
         goto 100
      endif

      call setvr2(ipdbt,mxnat,icbet,5)
      icc = 0
      icn = 0
      ico = 0
      ics = 0
      ich = 0
      do i=1,iconn(1,icbet)
         ii = abs(iconn(i+1,icbet))
         call setvar(iresid,mxnat,ii,k)
         if (ianz(ii).eq.6) icc = icc + 1
         if (ianz(ii).eq.7) icn = icn + 1
         if (ianz(ii).eq.8) ico = ico + 1
         if (ianz(ii).eq.16) ics = ics + 1
         if (ianz(ii).eq.1) then
            ich = ich + 1
            call setvr2(ipdbt,mxnat,ii,6 + ich)
         endif
      end do
      itot = icc + icn + ico + ics

      if (itot.eq.1) then
c Alanine
         iamino(k) = 2
         goto 100
      endif

      if (itot.eq.2.and.icc.eq.1.and.ico.eq.1) then
c Serine
         iamino(k) = 3
         do i=1,iconn(1,icbet)
            ii = abs(iconn(i+1,icbet))
            if (ianz(ii).eq.8) then
               call setvr2(ipdbt,mxnat,ii,31)
               do n=1,iconn(1,ii)
                  nn = abs(iconn(n+1,ii))
                  if (ianz(nn).eq.1) then 
                      call setvar(iresid,mxnat,nn,k)
                      call setvr2(ipdbt,mxnat,nn,10)
                  endif
               end do
            endif
         end do
         goto 100
      endif

      if (itot.eq.2.and.icc.eq.1.and.ics.eq.1) then
c Cysteine
         iamino(k) = 4
         do i=1,iconn(1,icbet)
            ii = abs(iconn(i+1,icbet))
            if (ianz(ii).eq.16) then
               call setvr2(ipdbt,mxnat,ii,37)
               do n=1,iconn(1,ii)
                  nn = abs(iconn(n+1,ii))
                  if (ianz(nn).eq.1) then
                     call setvr2(ipdbt,mxnat,nn,10)
                     call setvar(iresid,mxnat,nn,k)
                  endif
               end do
            endif
         end do
         goto 100
      endif

      if (itot.eq.3.and.icc.eq.2.and.ico.eq.1) then
c Threonine
         iamino(k) = 5
         do i=1,iconn(1,icbet)
            ii = abs(iconn(i+1,icbet))
            ia = ianz(ii)
            ich = 0
            if ((ia.eq.6.and.ii.ne.icalp)
     &          .or.ia.eq.8) then
               if (ia.eq.8) call setvr2(ipdbt,mxnat,ii,32)
               if (ia.eq.6) call setvr2(ipdbt,mxnat,ii,8)
               do n=1,iconn(1,ii)
                  nn = abs(iconn(n+1,ii))
                  if (ianz(nn).eq.1) then
                     if (ia.eq.8) then
                        call setvr2(ipdbt,mxnat,nn,13)
                     else
                        ich = ich + 1
                        call setvr2(ipdbt,mxnat,nn,15 + ich)
                     endif
                  endif
                  call setvar(iresid,mxnat,nn,k)
               end do
            endif
         end do
         goto 100
      endif

      if (itot.eq.3.and.icc.eq.3) then

         icg1 = 0
         do i=1,iconn(1,icbet)
            ii = abs(iconn(i+1,icbet))
            if (ianz(ii).eq.6.and.ii.ne.icalp) then
               do n=1,iconn(1,ii)
                  nn = abs(iconn(n+1,ii))
                  ia = ianz(nn)
                  if (ia.eq.6.and.nn.ne.icbet) icg1 = ii
               end do
               if (icg1.eq.0) icg1 = ii
            endif
         end do

         do i=1,iconn(1,icbet)
            ii = abs(iconn(i+1,icbet))
            if (ianz(ii).eq.6.and.ii.ne.icalp) then
               call setvr2(ipdbt,mxnat,ii,8)
               if (ii.eq.icg1) call setvr2(ipdbt,mxnat,ii,7)
               ih = 0
               do n=1,iconn(1,ii)
                  nn = abs(iconn(n+1,ii))
                  ia = ianz(nn)
                  call setvar(iresid,mxnat,nn,k)
                  if (ia.eq.6.and.nn.ne.icbet) then
c Isoleucine
                     iamino(k) = 6
                     call setvr2(ipdbt,mxnat,nn,10)
                     ich2 = 0
                     do l=1,iconn(1,nn)
                       ll = abs(iconn(l+1,nn))
                       call setvar(iresid,mxnat,ll,k)
                       if (ianz(ll).eq.1) then
                          ich2 = ich2 + 1
                          call setvr2(ipdbt,mxnat,ll,21+ich2)
                       endif
                     end do
 
                  elseif (ia.eq.1) then

                     ih = ih + 1
                     call setvr2(ipdbt,mxnat,nn,15+ih)
                     if (ii.eq.icg1) call setvr2(ipdbt,mxnat,nn,12+ih)

                  endif
               end do

            endif
         end do
c Valine
         if (iamino(k).ne.6) iamino(k) = 7
         goto 100
      endif

      if (itot.ne.2) print*, 'error in classifying aminoacids'

      icgam = 0
      do i=1,iconn(1,icbet)
         ii = abs(iconn(i+1,icbet))
         if (ianz(ii).eq.6.and.ii.ne.icbet.and.ii.ne.icalp) then
            icgam = ii
            call setvr2(ipdbt,mxnat,icgam,6)
         endif
      end do
      icc = 0
      icn = 0
      ico = 0
      ics = 0
      ich = 0
      do i=1,iconn(1,icgam)
         ii = abs(iconn(i+1,icgam))
         ia = ianz(ii)
         call setvar(iresid,mxnat,ii,k)
         if (ia.eq.1) then
             ich = ich + 1
             call setvr2(ipdbt,mxnat,ii,9+ich)
         endif
         if (ia.eq.6) icc = icc + 1
         if (ia.eq.7) icn = icn + 1
         if (ia.eq.8) ico = ico + 1
         if (ia.eq.16) ics = ics + 1
      end do
      itot = icc + icn + ico + ics

      if (itot.eq.2.and.icc.eq.1.and.ics.eq.1) then
c Methionine
          iamino(k) = 8
          do i=1,iconn(1,icgam)
             ii = abs(iconn(i+1,icgam))
             if (ianz(ii).eq.16) then
                call setvr2(ipdbt,mxnat,ii,36)
                do n=1,iconn(1,ii)
                   nn = abs(iconn(n+1,ii))
                   if (ianz(nn).eq.6.and.nn.ne.icgam) then
                     call setvr2(ipdbt,mxnat,nn,12)
                     call setvar(iresid,mxnat,nn,k)
                     ich = 0
                     do l=1,iconn(1,nn)
                       ll = abs(iconn(l+1,nn))
                       if (ianz(ll).eq.1) then
                          call setvar(iresid,mxnat,ll,k)
                          ich = ich + 1
                          call setvr2(ipdbt,mxnat,ll,27+ich)
                       endif
                     end do
                   endif
                end do
             endif
          end do
          goto 100
      endif

      if (itot.eq.2.and.icc.eq.2) then
         icdel = 0
         do i=1,iconn(1,icgam)
            ii = abs(iconn(i+1,icgam))
            if(ianz(ii).eq.6.and.ii.ne.icbet) icdel = ii
         end do
         icc = 0
         icn = 0
         ico = 0
         ics = 0
         do i=1,iconn(1,icdel)
            ii = abs(iconn(i+1,icdel))
            call setvar(iresid,mxnat,ii,k)
            if (ianz(ii).eq.6) icc = icc + 1
            if (ianz(ii).eq.7) icn = icn + 1
            if (ianz(ii).eq.8) ico = ico + 1
            if (ianz(ii).eq.16) ics = ics + 1
         end do
         itot = icc + icn + ico + ics

         if (itot.eq.2.and.icc.eq.2) then
c Lysine
          iamino(k) = 12
          call setvr2(ipdbt,mxnat,icdel,9)
          ich = 0
          do i=1,iconn(1,icdel)
             ii = abs(iconn(i+1,icdel))
             ia = ianz(ii)
             if (ia.eq.6.and.ii.ne.icgam) then
                call setvr2(ipdbt,mxnat,ii,12)
                ich2 = 0
                do n=1,iconn(1,ii)
                   nn = abs(iconn(n+1,ii))
                   call setvar(iresid,mxnat,nn,k)
                   if (ianz(nn).eq.7) then
                     call setvr2(ipdbt,mxnat,nn,27)
                     inh = 0
                     do l=1,iconn(1,nn)
                       ll = abs(iconn(l+1,nn))
                       call setvar(iresid,mxnat,ll,k)
                       if (ianz(ll).eq.1) then
                         inh = inh + 1
                         call setvr2(ipdbt,mxnat,ll,39+inh)
                       endif
                     end do
                   endif
                   if (ianz(nn).eq.1) then
                     ich2 = ich2 + 1
                     call setvr2(ipdbt,mxnat,nn,27+ich2)
                   endif
                end do
             endif
             if (ia.eq.1) then
                ich = ich + 1
                call setvr2(ipdbt,mxnat,ii,18+ich)
             endif
          end do
          goto 100
         endif

         if (itot.eq.2.and.icc.eq.1.and.icn.eq.1) then

            call setvr2(ipdbt,mxnat,icdel,9)
            ich = 0
            do i=1,iconn(1,icdel)
               ii = abs(iconn(i+1,icdel))
               ia = ianz(ii)
               if (ia.eq.1) then
                  ich = ich + 1
                  call setvr2(ipdbt,mxnat,ii,18+ich)
               endif
            end do

            do i=1,iconn(1,icdel)
               ii = abs(iconn(i+1,icdel))
               ia = ianz(ii)
               if (ia.eq.7) then
                  if (icalf(2,k).eq.ii) then
c Proline
                      iamino(k) = 15
                      goto 100
                  else
c Arginine
                     call setvr2(ipdbt,mxnat,ii,22)
                     iamino(k) = 16
                  endif
                  do n=1,iconn(1,ii)
                     nn = abs(iconn(n+1,ii))
                     call setvar(iresid,mxnat,nn,k)
                     if (ianz(nn).eq.6.and.nn.ne.icdel) then
                        call setvr2(ipdbt,mxnat,nn,17)
                        inn = 0
                        do l=1,iconn(1,nn)
                           ll = abs(iconn(l+1,nn))
                           call setvar(iresid,mxnat,ll,k)
                           if (ianz(ll).eq.7.and.ll.ne.ii) then
                              inn = inn + 1
                              call setvr2(ipdbt,mxnat,ll,24+inn)
                              inh = 0
                              do m=1,iconn(1,ll)
                                 mm = abs(iconn(m+1,ll))
                                 call setvar(iresid,mxnat,mm,k)
                                 if (ianz(mm).eq.1) then
                                    inh = inh + 1
                                    if (inn.eq.1) then
                                       call 
     &                                 setvr2(ipdbt,mxnat,mm,54+inh)
                                    else
                                       call 
     &                                 setvr2(ipdbt,mxnat,mm,57+inh)
                                    endif
                                 endif
                              end do
                           endif
                        end do
                     endif
                     if (ianz(nn).eq.1) call setvr2(ipdbt,mxnat,nn,28)
                  end do
                  goto 100
               endif
            end do
         endif

         if (itot.eq.3.and.icc.eq.1.and.ico.eq.2) then
c Glutamate
            call setvr2(ipdbt,mxnat,icdel,9)
            iamino(k) = 13
            io = 0
            do i=1,iconn(1,icdel)
               ii = abs(iconn(i+1,icdel))
               if (ianz(ii).eq.8) then
                  io = io + 1
                  call setvr2(ipdbt,mxnat,ii,33+io)
                  do n=1,iconn(1,ii)
                     nn = abs(iconn(n+1,ii))
                     call setvar(iresid,mxnat,nn,k)
                  end do
               endif
            end do
            goto 100
         endif

         if (itot.eq.3.and.icc.eq.1.and.ico.eq.1.and.icn.eq.1) then
c Glutamine
            call setvr2(ipdbt,mxnat,icdel,9)
            iamino(k) = 14
            do i=1,iconn(1,icdel)
               ii = abs(iconn(i+1,icdel))
               ia = ianz(ii)
               if (ia.eq.8) then
                  call setvr2(ipdbt,mxnat,ii,34)
               elseif (ia.eq.7) then
                  call setvr2(ipdbt,mxnat,ii,24)
                  inh = 0
                  do n=1,iconn(1,ii)
                     nn = abs(iconn(n+1,ii))
                     call setvar(iresid,mxnat,nn,k)
                     if (ianz(nn).eq.1) then
                        inh = inh + 1
                        call setvr2(ipdbt,mxnat,nn,33+inh)
                     endif
                  end do
               endif
            end do
            goto 100
         endif
      endif

      if (itot.eq.3.and.icc.eq.2.and.ico.eq.1) then
c HydroxyProline
         iamino(k) = 15
         do i=1,iconn(1,icgam)
            ii = abs(iconn(i+1,icgam))
            if (ianz(ii).eq.8.or.ianz(ii).eq.6) then
               do n=1,iconn(1,ii)
                  nn = abs(iconn(n+1,ii))
                  call setvar(iresid,mxnat,nn,k)
               end do
            endif
         end do
         goto 100
      endif

      if (itot.eq.3.and.icc.eq.1.and.ico.eq.2) then
c Aspartate
         iamino(k) = 9
         io = 0
         do i=1,iconn(1,icgam)
            ii = abs(iconn(i+1,icgam))
            if (ianz(ii).eq.8) then
               io = io + 1
               call setvr2(ipdbt,mxnat,ii,28+io)
               do n=1,iconn(1,ii)
                  nn = abs(iconn(n+1,ii))
                  call setvar(iresid,mxnat,nn,k)
               end do
            endif
         end do
         goto 100
      endif

      if (itot.eq.3.and.icc.eq.2.and.icn.eq.1) then
c Histidine
         iamino(k) = 17
         do i=1,iconn(1,icgam)
            ii = abs(iconn(i+1,icgam))
            ia = ianz(ii)
            if ((ia.eq.6.and.ii.ne.icbet).or.ia.eq.7) 
     &      then
               if (ia.eq.6) call setvr2(ipdbt,mxnat,ii,11)
               if (ia.eq.7) call setvr2(ipdbt,mxnat,ii,20)
               do n=1,iconn(1,ii)
                  nn = abs(iconn(n+1,ii))
                  na = ianz(nn)
                  if (nn.ne.icgam) then
                     if (na.eq.1.or.na.eq.6.or.na.eq.7) 
     &                  call setvar(iresid,mxnat,nn,k)
                     if (na.eq.1) then
                        if (ia.eq.7) then
                           call setvr2(ipdbt,mxnat,nn,22)
                        else
                           call setvr2(ipdbt,mxnat,nn,25)
                        endif
                     endif
                     if (na.eq.6) call setvr2(ipdbt,mxnat,nn,13)
                     if (na.eq.7) call setvr2(ipdbt,mxnat,nn,24)
                     do l=1,iconn(1,nn)
                        ll = abs(iconn(l+1,nn))
                        if (ianz(ll).eq.1) then
                           call setvar(iresid,mxnat,ll,k)
                           if (na.eq.6) then
                              call setvr2(ipdbt,mxnat,ll,31)
                           else
                              call setvr2(ipdbt,mxnat,ll,34)
                           endif
                        endif
                     end do
                  endif
               end do
            endif
         end do
         goto 100
      endif

      if (itot.eq.3.and.icc.eq.1.and.ico.eq.1.and.icn.eq.1) then
c Asparagine
         iamino(k) = 10
         do i=1,iconn(1,icgam)
            ii = abs(iconn(i+1,icgam))
            ia = ianz(ii)
            if (ia.eq.8.or.ia.eq.7) then
               if (ia.eq.8) call setvr2(ipdbt,mxnat,ii,29)
               if (ia.eq.7) then
                  call setvr2(ipdbt,mxnat,ii,21)
                  inh = 0
                  do n=1,iconn(1,ii)
                     nn = abs(iconn(n+1,ii))
                     call setvar(iresid,mxnat,nn,k)
                     if (ianz(nn).eq.1) then
                        inh = inh + 1
                        call setvr2(ipdbt,mxnat,nn,24+inh)
                     endif
                  end do
               endif
            endif
         end do
         goto 100
      endif

      if (itot.eq.3.and.icc.eq.3) then

         itrp = 0
         ice2 = 0
         inn = 0
         do i=1,iconn(1,icgam)
            ii = abs(iconn(i+1,icgam))
            if (ianz(ii).eq.6.and.ii.ne.icbet) then
               do n=1,iconn(1,ii)
                  nn = abs(iconn(n+1,ii))
                  if (nn.ne.icgam.and.ianz(nn).eq.7) then
c Tryptophan
                      iamino(k) = 20
                      itrp = 1
                      inn = nn
                  endif
               end do
               if (itrp.eq.1) then
                  do n=1,iconn(1,ii)
                     nn = abs(iconn(n+1,ii))
                     if (nn.ne.icgam.and.ianz(nn).eq.6) then
                        do m=1,iconn(1,nn)
                          mm = abs(iconn(m+1,nn))
                          if (mm.ne.ii.and.mm.eq.inn) ice2 = nn
                        end do
                     endif
                  end do
               endif
            endif
         end do
    
         l = 0
         icc = 0
         do i=1,2
            icd(i) = 0
            icdh(i) = 0
            icdhn(i) = 0
         end do

         do i=1,iconn(1,icgam)
            ii = abs(iconn(i+1,icgam))
            if (ianz(ii).eq.6.and.ii.ne.icbet) then
               icc = icc + 1
               icd(icc) = ii
               if (icc.eq.1) then
                  call setvr2(ipdbt,mxnat,ii,10)
               else
                  call setvr2(ipdbt,mxnat,ii,11)
               endif
               do n=1,iconn(1,ii)
                  nn = abs(iconn(n+1,ii))
                  na = ianz(nn)
                  if (nn.ne.icgam) then

                     call setvar(iresid,mxnat,nn,k)
                  
                     if (na.eq.6) then
                         l = 1
                         if (itrp.eq.1) then
                            if (nn.eq.ice2) then
                               call setvr2(ipdbt,mxnat,nn,14)
                            else
                               call setvr2(ipdbt,mxnat,nn,15)
                            endif
                         else
                            if (icc.eq.1) then
                               call setvr2(ipdbt,mxnat,nn,13)
                            else
                               call setvr2(ipdbt,mxnat,nn,14)
                            endif
                         endif
                     elseif (na.eq.1) then
                         icdh(icc) = nn
                         icdhn(icc) = icdhn(icc) + 1
                         if (icc.eq.1) then
                            call setvr2(ipdbt,mxnat,nn,22+icdhn(1)-1)
                         else
                            call setvr2(ipdbt,mxnat,nn,25+icdhn(2)-1)
                         endif
                     elseif (na.eq.7) then
c Tryptophan
                         iamino(k) = 20
                         call setvr2(ipdbt,mxnat,nn,23)
                         if (ii.ne.icd(1)) then
                            call setvr2(ipdbt,mxnat,icd(1),11)
                            call setvr2(ipdbt,mxnat,icd(2),10)
                            call setvr2(ipdbt,mxnat,icdh(1),25)
                            call setvr2(ipdbt,mxnat,icdh(2),22)
                         endif
                     endif

                     do m=1,iconn(1,nn)
                       mm = abs(iconn(m+1,nn))
                       ma = ianz(mm)
                       call setvar(iresid,mxnat,mm,k)
                       if (mm.ne.ii) then
                        if (ma.eq.6.and.na.ne.7) then
                        
                           if (itrp.eq.1) then
                              if (nn.eq.ice2) then
                                 call setvr2(ipdbt,mxnat,mm,18)
                              else
                                 call setvr2(ipdbt,mxnat,mm,19)
                              endif
                           else
                              call setvr2(ipdbt,mxnat,mm,17)
                           endif

                           do j=1,iconn(1,mm)
                              jj = abs(iconn(j+1,mm))
                              ja = ianz(jj)
                              if (jj.ne.nn) then
                              call setvar(iresid,mxnat,jj,k)
c Tyrosine
                              if (ja.eq.8) then
                                 iamino(k) = 19
                                 call setvr2(ipdbt,mxnat,jj,33)
                              elseif (ja.eq.6.and.itrp.eq.1) then
                                 call setvr2(ipdbt,mxnat,jj,16)
                              elseif (ja.eq.1) then
                                 if (itrp.eq.1) then
                                    if (nn.eq.ice2) then
                                       call setvr2(ipdbt,mxnat,jj,46)
                                    else
                                       call setvr2(ipdbt,mxnat,jj,49)
                                    endif
                                 else
                                    call setvr2(ipdbt,mxnat,jj,40)
                                 endif
                              endif
                              do l=1,iconn(1,jj)
                                 ll = abs(iconn(l+1,jj))
                                 if (ianz(ll).eq.1.and.ll.ne.mm) then
                                    call setvar(iresid,mxnat,ll,k)
                                    if (ja.eq.8) then
                                       call setvr2(ipdbt,mxnat,ll,52)
                                    else
                                       if (itrp.eq.1)
     &                                 call setvr2(ipdbt,mxnat,ll,58)
                                    endif
                                 endif
                              end do
                              endif
                           end do

                        elseif (ma.eq.1) then

                           if (na.eq.7) then
                              call setvr2(ipdbt,mxnat,mm,31)
                           else
                              if (itrp.eq.1) then
                                 call setvr2(ipdbt,mxnat,mm,37)
                              else
                                 if (icc.eq.1) then
                                    call setvr2(ipdbt,mxnat,mm,31)
                                 else
                                    call setvr2(ipdbt,mxnat,mm,34)
                                 endif
                              endif
                           endif

                        endif

                       endif
                     end do
                  endif
               end do
            endif
         end do
         if (l.eq.0) then
c Leucine
            iamino(k) = 11
            goto 100
         endif
      endif
c Phenylalanine
      if (iamino(k).ne.20.and.iamino(k).ne.19) iamino(k) = 18

100   continue
      end do

      istat = 1

      return
      end

      subroutine setvar(iarr,nelm,ielm,ival)
      implicit double precision (a-h,o-z)
      dimension iarr(*)

      if (ielm.gt.0.and.ielm.le.nelm) iarr(ielm) = ival
         
      return
      end

      subroutine setvr2(iarr,nelm,ielm,ival)
      implicit double precision (a-h,o-z)
      integer*2 iarr
      dimension iarr(*)

      if (ielm.gt.0.and.ielm.le.nelm) iarr(ielm) = ival
         
      return
      end

      logical function chkhs(ihpdb)
      implicit double precision (a-h,o-z)
      parameter (mxhsym=64)
      dimension ihpdb(*)

      chkhs = .false.
      do i=2,mxhsym*3
         if (ihpdb(i).ne.0) chkhs = .true.
      end do

      return
      end

      subroutine getpdd(ires,ipdb,ihpdb,
     &                  ianz,iresid,ipdbt)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat

      parameter (mxsym=103)
      parameter (mxhsym=64)
      integer*2 ipdbt
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension ianz(*),iresid(*),ipdbt(*)

      do i=1,mxsym
          ipdb(i) = 0
      end do
      do i=1,mxhsym*3
          ihpdb(i) = 0
      end do

      do i=1,iatoms
          if (iresid(i).eq.ires) then
              if (ipdbt(i).ne.0) then
                  if (ianz(i).eq.1) then
                     ihpdb(ipdbt(i)) = i
                  else
                     ipdb(ipdbt(i)) = i
                  endif
              endif
          endif
      end do

      return
      end

      subroutine bckod(ibckok,ires,iop,
     &                  ianz,iresid,ipdbt,iamino)
      implicit double precision (a-h,o-z)

      parameter (mxsym=103)
      common /athlp/ iatoms, mxnat
      integer*2 ipdbt
      dimension ipdb(mxsym)
      dimension ianz(*),iresid(*),ipdbt(*),iamino(*)

c checks if Calpha is present (iop=1)

      ibckok = 0

      do i=1,mxsym
          ipdb(i) = 0
      end do

      do i=1,iatoms
          if (iresid(i).eq.ires) then
              if (ipdbt(i).ne.0) then
                  if (ianz(i).ne.1) then
                     ipdb(ipdbt(i)) = i
                  endif
              endif
          endif
      end do

c has Calpha
      if (iamino(ires).gt.23) then
         if (iop.eq.1.and.(ipdb(43).ne.0.or.ipdb(46).ne.0)) 
     &       ibckok = 1
      else
         if (iop.eq.1.and.ipdb(2).ne.0) ibckok = 1
         if (iop.eq.2.and.
     &      (ipdb(1).ne.0.and.ipdb(2).ne.0.and.ipdb(3).ne.0)
     &      ) ibckok = 1
      endif

      return
      end
