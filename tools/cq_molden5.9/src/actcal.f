      subroutine actcad(iopt,
     &                  iaton,iatclr,iresid,
     &                  icalf,ianf,islu,nchain,iamino,ibck)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      parameter (mxres=42)
      parameter (mxchai=50)
      common /clfhlp/ isndcl(4),iamicl(mxres),ichcol(mxchai)
      common /cllab/  iclon,iclpnt(4)
      dimension iaton(*),iatclr(*),iresid(*),iamino(*)
      dimension icalf(6,*),ianf(*),islu(*),ibck(*)

      icllow = 0
      if (iclon.eq.1) icllow = iclpnt(1) 

      if (iopt.eq.1) then
         do i=1,iatoms
             iaton(i) = 0
         end do

         do i=1,nchain
             do j=ianf(i),islu(i)
                if (iamino(j).gt.23) then
                   n = 6
                else
                   n = 3
                endif
                do k=1,n
                   if (icalf(k,j).ne.0) then
                       iaton(icalf(k,j)) = 1
                       iatclr(icalf(k,j)) = ichcol(i)
                   endif
                end do
             end do
         end do
      end if

      if (iopt.le.0) then
         do i=1,iatoms
             iaton(i) = 1
             if (iclon.eq.0.or.(iclon.eq.1.and.i.lt.icllow)) then
                if (iresid(i).eq.0) then
                   iatclr(i) =  isndcl(1)
                elseif (iresid(i).eq.-1) then
                   iatclr(i) =  isndcl(2)
                elseif (iresid(i).eq.-2) then
                   iatclr(i) =  isndcl(3)
                elseif (iresid(i).eq.-3) then
                   iatclr(i) =  isndcl(4)
                elseif (iopt.eq.0.and.iresid(i).lt.-3) then
                   iatclr(i) =  1
                else
                   if (iopt.eq.0) iatclr(i) =  12
                endif
             endif
         end do
      end if

      do i=1,4
         ibck(i) = 1
      end do

      return
      end

      subroutine acthed(iopt,iscnd,jcolsp,inclbb,
     &                  iaton,iatclr,iresid,iconn,
     &                  icalf,ianf,islu,nchain,iamino,ihet,
     &                  reson,isal,ibck)
c     switches off the backbone and activates helix, beta sheet, coil
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (mxel=100)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      common /athlp/ iatoms, mxnat
      parameter (mxres=42)
      parameter (mxchai=50)
      common /clfhlp/ isndcl(4),iamicl(mxres),ichcol(mxchai)
      integer reson
      logical edge
      dimension iaton(*),iatclr(*),iresid(*),iconn(mxcon+1,*)
      dimension icalf(6,*),ianf(*),islu(*),iamino(*),ihet(*),reson(*),
     &          isal(*),ibck(*)

      if (iscnd.lt.0.or.iscnd.gt.3) return


      jcol = jcolsp
      if (jcol.eq.0) then
         if (iscnd.eq.0) then
            jcol = isndcl(1)
         elseif (iscnd.eq.1) then
            jcol = isndcl(2)
         elseif (iscnd.eq.2) then
            jcol = isndcl(3)
         elseif (iscnd.eq.3) then
            jcol = isndcl(4)
         else
            jcol = 12
         endif
      endif

      icol(100) = 10
      ifnd = 0
      do i=1,iatoms
         if (iresid(i).eq.-iscnd) then
            iaton(i) = iopt
            iatclr(i) = jcol
            ifnd = 1
         endif
      end do

      if (ifnd.eq.0) then
         if (iscnd.eq.0) then
             call inferr('No Helix found !',0)
         elseif (iscnd.eq.1) then
             call inferr('No Beta found !',0)
         elseif (iscnd.eq.2) then
             call inferr('No RNA/DNA backbone found !',0)
         elseif (iscnd.eq.3) then
             call inferr('No Coil found !',0)
         endif
      endif

      iset = inclbb
      if (iopt.eq.0) iset = 1
      ibck(iscnd+1) = iset

      do i=1,nchain
          istr = 0
          ichlen = islu(i) - ianf(i) + 1
          do k=1,ichlen
             j = ianf(i) + k - 1
             edge = .false.
             ido = 0
             if (iscnd.eq.3) then

                if (isal(j).eq.iscnd.and.reson(j).eq.0) ido = 1

                if (isal(j).ne.iscnd.and.reson(j).eq.0) then
                 if (iamino(j).le.23) then
                    if (.not.(k.eq.1.or.k.eq.ichlen)) then
                       if (isal(j-1).eq.iscnd.or.
     &                     isal(j+1).eq.iscnd) then
                         ido = 1
                         if ((ihet(isal(j)+1).eq.0.and.iset.eq.0) 
     &                   .or.(ihet(isal(j)+1).eq.1.and.iset.eq.1))
     &                   edge = .true.
                       endif
                    endif
                 endif
                endif

                if ((isal(j).eq.0.and.isal(j+1).eq.1).or.
     &              (isal(j).eq.1.and.isal(j+1).eq.0).or.
     &              (isal(j).eq.0.and.isal(j-1).eq.1).or.
     &              (isal(j).eq.1.and.isal(j-1).eq.0)) then
                   if (.not.(k.eq.1.or.k.eq.ichlen)) then
                      if ((isal(j).eq.0.and.isal(j+1).eq.1).or.
     &                    (isal(j).eq.1.and.isal(j+1).eq.0)) ido = 2
                      if ((isal(j).eq.0.and.isal(j-1).eq.1).or.
     &                    (isal(j).eq.1.and.isal(j-1).eq.0)) ido = 3
                      if (ihet(isal(j)+1).eq.0.and.iset.eq.0) 
     &                    edge = .true.
                      if (ihet(isal(j)+1).eq.1.and.iset.eq.1) 
     &                    edge = .true.
                   endif
                endif

             else

                if (isal(j).eq.iscnd.and.reson(j).eq.0) then
                 ido = 1
                 if (iamino(j).le.23) then
                    if (.not.(k.eq.1.or.k.eq.ichlen)) then

                       if (.not.(isal(j-1).eq.iscnd.and.
     &                     isal(j+1).eq.iscnd)) then
                         if (isal(j-1).ne.iscnd) then
                            if (isal(j-1).eq.3) then
                              if (ihet(isal(j-1)+1).eq.0.and.iset.eq.0)
     &                            edge = .true.
                              if (ihet(isal(j-1)+1).eq.1.and.iset.eq.1)
     &                            edge = .true.
                            else
                              if (ihet(4).eq.0.and.iset.eq.0) 
     &                            edge = .true.
                              if (ihet(4).eq.1.and.iset.eq.1) 
     &                            edge = .true.
                            endif
                         endif
                         if (isal(j+1).ne.iscnd) then
                            if (isal(j+1).eq.3) then
                              if (ihet(isal(j+1)+1).eq.0.and.iset.eq.0)
     &                            edge = .true.
                              if (ihet(isal(j+1)+1).eq.1.and.iset.eq.1)
     &                            edge = .true.
                            else
                              if (ihet(4).eq.0.and.iset.eq.0) 
     &                            edge = .true.
                              if (ihet(4).eq.1.and.iset.eq.1) 
     &                            edge = .true.
                            endif
                         endif
                       endif

                    endif
                 endif
                endif
             endif

             if (ido.gt.0) then
                if (edge) then

                   ioff = 1
                   if (ido.eq.1.and.isal(j-1).eq.iscnd) ioff = -1
                   if (ido.eq.3) ioff = -1
                   if (iamino(j).le.23) then
                      if (iset.eq.1) iaton(icalf(1,j)) = 1
                      do l=2,3
                         i1 = icalf(l,j)
                         do n=1,iconn(1,i1)
                            do ll=2,3
                               i2 = icalf(ll,j+ioff)
                               if (abs(iconn(n+1,i1)).eq.i2) 
     &                             iaton(i1) = iset
                            end do
                         end do
                      end do
                   endif

                else

                   n=6
                   if (iamino(j).le.23) n=3
                   do l=1,n
                      if (icalf(l,j).ne.0) then
                         iaton(icalf(l,j)) = iset
                      endif
                   end do

                endif
             endif

          end do
      end do

      ihet(iscnd+1) = iopt

      return
      end
