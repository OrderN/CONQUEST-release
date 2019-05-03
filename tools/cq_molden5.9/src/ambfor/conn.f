      subroutine cond34(istat,nac,iac,nad,iad,iconn)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxac=3*mxcon)
      parameter (mxad=9*mxcon)
      common /athlp/  iatoms, mxnat
      dimension iconn(mxcon+1,*), 
     &          nac(*),nad(*),iac(mxac,*),iad(mxad,*)

      istat = 1

      do i=1,iatoms
         nac(i) = 0

         do j=1,iconn(1,i)
            jj = iconn(1+j,i)

            if (jj.gt.0) then

               do k=1,iconn(1,jj)
                  kk = iconn(1+k,jj)

                  if (kk.ne.i.and.kk.gt.0) then
                     ido = 1

c check if 1-3 atom is not already 1-2

                     do l=1,iconn(1,i)
                        if (kk .eq. iconn(1+l,i)) ido = 0
                     end do

                     if (nac(i).lt.mxac) then
                        if (ido.eq.1) then
                           nac(i) = nac(i) + 1
                           iac(nac(i),i) = kk
                        endif
                     else
                        istat = 0
                        return
                     endif

                  endif

               end do

            endif

         end do

      end do

      do i=1,iatoms
         nad(i) = 0

         do j=1,nac(i)
            jj = iac(j,i)

            do k =1,iconn(1,jj)
               kk = iconn(1+k,jj)

               if (kk.ne.i.and.kk.gt.0) then
                  ido = 1

c check if 1-4 atom is not already 1-3, or 1-2

                  do l=1,iconn(1,i)
                     if (kk.eq.iconn(1+l,i)) ido = 0
                  end do
                  do l=1,nac(i)
                     if (kk.eq.iac(l,i)) ido = 0
                  end do
                  do l=1,nad(i)
                     if (kk.eq.iad(l,i)) ido = 0
                  end do

                  if (nad(i).lt.mxad) then
                     if (ido.eq.1) then
                        nad(i) = nad(i) + 1
                        iad(nad(i),i) = kk
                     endif
                  else
                     istat = 0
                     return
                  endif

               endif

            end do

         end do

      end do

      return
      end

