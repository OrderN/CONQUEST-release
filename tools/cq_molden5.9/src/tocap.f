      subroutine tocap(line,ncars)
      integer ncars
      character*(*) line

      do i=1,ncars
         j = ichar(line(i:i))
         if (j.ge.97.and.j.le.122) line(i:i) = char(j-32)
      end do

      return
      end

      subroutine tolow(line,ncars)
      integer ncars
      character*(*) line

      do i=1,ncars
         j = ichar(line(i:i))
         if (j.ge.65.and.j.le.90) line(i:i) = char(j+32)
      end do

      return
      end
