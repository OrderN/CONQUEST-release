      subroutine parfd(coo,fscal)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      parameter (maxfat=1000)
      parameter (mxel=100)
      parameter (mxcon=10)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat) 
      common /athlp/ iatoms, mxnat

      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      dimension coo(3,*)

      fdefsc = 10.0d0
      if (iftyp.eq.2.or.iftyp.eq.3) fdefsc = -fdefsc

      do i=1,iatoms
         do j=1,3
             fc(j,i) = coo(j,i) + fdefsc*fscal*fxyz(j,i)
         end do
      end do

      call parptr(135,fc,fdum,idum)

      return
      end
