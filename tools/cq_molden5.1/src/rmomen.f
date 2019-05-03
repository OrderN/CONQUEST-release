      double precision function rmomen(b)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (maxfat=1000)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat) 
      common /athlp/ iatoms, mxnat

      dimension temp(3),rmom(3),b(3,numat1)

      call rotcor(b)

      do j=1,3
          rmom(j) = 0.0d0
      end do

      do i=1,iatoms
        call crpros(b(1,i),fxyz(1,i),temp)
        do j=1,3
           rmom(j) = rmom(j) + temp(j)
        end do
      end do

      rmomen = vlen(rmom)

      return
      end
