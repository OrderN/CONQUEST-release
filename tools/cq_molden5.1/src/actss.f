      subroutine actsd(iopt,
     &                 ianz,iaton,iatclr,iresid,iconn)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      character*2 gstr
      character lstr*137
      dimension ianz(*),iaton(*),iatclr(*),iresid(*),iconn(mxcon+1,*)

c     shouldnt issdon be set here ?

      iss = 0
      do i=1,iatoms
         if (ianz(i).eq.16) then
             do j=1,iconn(1,i)
                if (ianz(abs(iconn(j+1,i))).eq.16) then
                   ires1 = iresid(i)
                   ires2 = iresid(abs(iconn(j+1,i)))
                   if (iopt.eq.1) then
                      call actami(ires1,0,1,0)
                      call actami(ires2,0,1,0)
                      iss = iss + 1
                      iaton(i) = 1
                      iaton(abs(iconn(j+1,i))) = 1
                      iatclr(i) = 2
                      iatclr(abs(iconn(j+1,i))) = 2
                   else
                      call actami(ires1,0,0,0)
                      call actami(ires2,0,0,0)
                      iaton(i) = 0
                      iaton(abs(iconn(j+1,i))) = 0
                   endif
                endif
             end do
         endif
      end do

      if (iopt.eq.1) then
         if (iss.eq.0) then
            call inferr('No Sulfur Bridges Found !',0)
         else
            lstr = gstr(iss/2)//' Sulfur Bridges Found'
            call inferr(lstr,0)
         endif
      endif

      return
      end
