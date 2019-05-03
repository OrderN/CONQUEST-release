      subroutine readel(line,n)
      character*(*) line
      character*137 lstr
      common /rdwr/ iun1,iun2,iun3,iun4,iun5

      do i=1,n
        read(iun2,'(a)',err=100,end=100) lstr
      end do
      line(1:137) = lstr(1:137)
c
      return
c
100   call inferr('Premature end of file !',1)
c
      end 
