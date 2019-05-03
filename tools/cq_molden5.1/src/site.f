      subroutine site(line)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (mxsite=300)
      parameter (lmax=4)
      parameter (lmaxf=(lmax+1)*(lmax+1))
      common /multip/ q(lmaxf,mxsite),car(3,mxsite),ctag(mxsite),
     &                nsites
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*8   ctag,attag
      character*137 line

10    read(line,11) nsites,attag 
11    format(5x,i3,7x,a8)
      ctag(nsites) = attag
      
      do i=1,lmaxf
         q(i,nsites) = 0.0d0
      end do

      read(iun2,'(a)') line
      read(iun2,'(a)') line
      read(iun2,'(a)') line

      read(line,12) car(1,nsites),car(2,nsites),car(3,nsites)
12    format(6x,f10.6,5x,f10.6,5x,f10.6)

15    read(iun2,'(a)') line
      if( index(line,'site').ne.0.) goto 10
      if(index(line,'total multipoles referred to origin').ne.0) return
      call getreal(line,'q00 ',q(1,nsites))
      call getreal(line,'q10 ',q(2,nsites))
      call getreal(line,'q11c',q(3,nsites))
      call getreal(line,'q11s',q(4,nsites))
      call getreal(line,'q20 ',q(5,nsites))
      call getreal(line,'q21c',q(6,nsites))
      call getreal(line,'q21s',q(7,nsites))
      call getreal(line,'q22c',q(8,nsites))
      call getreal(line,'q22s',q(9,nsites))
      call getreal(line,'q30 ',q(10,nsites))
      call getreal(line,'q31c',q(11,nsites))
      call getreal(line,'q31s',q(12,nsites))
      call getreal(line,'q32c',q(13,nsites))
      call getreal(line,'q32s',q(14,nsites))
      call getreal(line,'q33c',q(15,nsites))
      call getreal(line,'q33s',q(16,nsites))
      goto 15
      end
