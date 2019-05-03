      subroutine proxid(itarg,thresh,ikleur,idosrf,coo,iresid,ishoh)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      logical uniq,debug
      parameter (mxnb=1000)
      dimension ineigh(mxnb)
      dimension coo(3,*),iresid(*)

      debug = .false.

      do i=1,mxnb
         ineigh(i) = 0
      end do
      numnb = 0

      do i=1,iatoms
         if (iresid(i).eq.itarg) then
            do j=1,iatoms
               if (iresid(j).ne.itarg.and.iresid(j).ne.-ishoh) then
                  dist = dsqrt(dist2(coo(1,i),coo(1,j)))
                  if (dist.lt.thresh) then
                     uniq = .true.
                     do k=1,numnb
                        if (ineigh(k).eq.iresid(j)) uniq = .false.
                     end do
                     if (uniq) then
                        numnb = numnb + 1
                        ineigh(numnb) = iresid(j)
                     endif
                  endif
               endif
            end do
         endif
      end do
      
      if (debug) then
         print*,'neighbouring residues'
         do i=1,numnb
            print*,i,' ',ineigh(i)
         end do
      endif

      do i=1,numnb
         call actami(ineigh(i),ikleur,1,idosrf)
      end do

      return
      end

      subroutine proxd(itarg,backb,adds,idocom,thresh,
     &                 coo,iresid,iaton,iconn,iams,ishoh)
      implicit double precision (a-h,o-z)
      parameter (maxdm=20)
      parameter (maxrs=100)
      parameter (mxcon=10)
      parameter (mxra=100)
      logical dozme,notcon,ispdb
      common /getpnt/irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &               iconv,ircus,dozme
      common /athlp/ iatoms, mxnat
      common /distmn/ rdm(maxdm),idmon(2,maxdm),ndm
      integer backb,adds
      dimension irs(mxra),irss(maxrs)
      dimension coo(3,*),iresid(*),iaton(*),iconn(mxcon+1,*),iams(*)

      call clrmon

      do i=1,maxdm
         rdm(i) = 1.0d10
      end do
      ndm = maxdm

      nrss = 0
      toang  = 0.52917706d0

      ispdb = (ipdbon.eq.1)

      ir1 = iresid(itarg)
      nrs = 0
      if (ispdb) then
         do i=1,iatoms
            ir2 = iresid(i)
            if (ir1.eq.ir2) then
               if (nrs.lt.mxra) then
                  nrs = nrs + 1
                  irs(nrs) = i
               endif
            endif
         end do
      else
         nrs = 1
         irs(1) = itarg
      endif

      do i=1,iatoms
        ir2 = iresid(i)
        if (iresid(i).gt.0.or.iresid(i).lt.-3) then
          if ((.not.ispdb.and.i.ne.itarg).or.
     &             (ispdb.and.ir1.ne.ir2)) then
c     &             (backb.eq.1.and.ir1.ne.ir2.and.ir2.gt.0)) then
           do n=1,nrs
             in = irs(n)
             dist = dsqrt(dist2(coo(1,i),coo(1,in)))
             if (dist.lt.thresh/toang) then
                ihrss = 0
                do j=1,nrss
                   if (irss(j).eq.ir2) ihrss = 1
                end do
                if (ihrss.eq.0) then
                   if (nrss.lt.maxrs) then
                      nrss = nrss + 1
                      irss(nrss) = ir2
                   endif
                endif
             endif
             do j=1,ndm
                if (dist.lt.rdm(j)) then
                    notcon = .true.
                    do k=1,iconn(1,in)
                       kk = iconn(1+k,in)
                       if (kk.eq.i) then
                           notcon = .false.
                       else
                          do l=1,iconn(1,iabs(kk))
                             ll = iconn(1+l,iabs(kk))
                             if (ll.eq.i) then
                                notcon = .false.
                             elseif (ll.ne.in) then
                                do m=1,iconn(1,iabs(ll))
                                   mm = iconn(1+m,iabs(ll))
                                   if (mm.eq.i) then
                                      notcon = .false.
                                   elseif (mm.ne.kk.and.backb.eq.1) then
                                      do nm=1,iconn(1,iabs(mm))
                                        km = iconn(1+nm,iabs(mm))
                                        if (km.eq.i) notcon = .false.
                                      end do
                                   endif
                                end do
                             endif
                          end do
                       endif
                    end do
                    if (notcon) then
                       do m=ndm-1,j,-1
                          idmon(1,m+1) = idmon(1,m)
                          idmon(2,m+1) = idmon(2,m)
                          rdm(m+1) = rdm(m)
                       end do
                       idmon(1,j) = in
                       idmon(2,j) = i
                       rdm(j) = dist
                       goto 100
                    endif
                endif
             end do
100          continue

           end do
          endif
        endif
      end do
      

      if (ispdb) then
         do i=1,nrss
            if (irss(i).ne.-ishoh) call actami(irss(i),0,1,0)
            if (adds.eq.1.and.irss(i).gt.0) iams(irss(i)) = 1
         end do
         do i=1,maxdm
            ir = iresid(idmon(2,i))
            if (ir.eq.-ishoh) then
               iaton(idmon(2,i)) = 1
            else
c               call actami(iresid(idmon(2,i)),0,1,0)
            endif
         end do
      else
         do i=1,maxdm
            iaton(idmon(2,i)) = 1
         end do
      endif

      if (idocom.eq.1) then
         call domcon(ndm,1)
      else
         call clrmon
      endif

      return
      end
