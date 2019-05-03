      subroutine mulprt
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (mxsite=300)
      parameter (lmax=4)
      parameter (lmaxf=(lmax+1)*(lmax+1))
      character*8   ctag
      common /multip/ q(lmaxf,mxsite),car(3,mxsite),ctag(mxsite),
     &                nsites

      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      write(iun3,*)' '
      write(iun3,*)'========== multipoles =========='
      write(iun3,*)' '

      do ij=1,nsites
         write(iun3,*)'***** site ',ij,' *****'
         write(iun3,*)' '
         write(iun3,*)'tag ',ctag(ij)
         write(iun3,*)' '
         write(iun3,*)'coordinates ',(car(j,ij),j=1,3)
         write(iun3,*)' '
         write(iun3,*)'q(',ij,')=',q(1,ij)
         do i=1,3
           write(iun3,*)'dip(',ij,',',i,')=', q(i+1,ij)
         end do
         do i=1,5
           write(iun3,*)'qua(',ij,',',i,')=', q(i+4,ij)
         end do
         do i=1,7
           write(iun3,*)'octa(',ij,',',i,')=',q(i+9,ij)
         end do
      end do

      write(iun3,*)' '
      write(iun3,*)'================================'
      write(iun3,*)' '

      return
      end
