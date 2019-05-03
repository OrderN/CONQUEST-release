      subroutine readel(line,n)
      character*(*) line
      common /rdwr/ iun1,iun2,iun3,iun4,iun5

      do i=1,n
        read(iun2,'(a)',end=100) line
      end do
c
      return
c
100   call inferr('Premature end of file !',1)
c
      end 
