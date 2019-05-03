      subroutine fndcal(itarg,isel,indx,istat,icalf,ncalf)
      implicit double precision (a-h,p-x),integer (i-n),logical (o)
      dimension icalf(6,*)

      istat = 1
      do i=1,ncalf
         if (icalf(isel,i).eq.itarg) then
            indx = i
            return
         endif
      end do

      istat = 0
      return
      end
