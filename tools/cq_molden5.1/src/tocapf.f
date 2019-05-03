      character*2 function tocapf(line)
      character*2 line

      tocapf = line

      do i=1,2
         j = ichar(tocapf(i:i))
         if (j.ge.97.and.j.le.122) tocapf(i:i) = char(j-32)
      end do

      return
      end

      character*2 function tolowf(line)
      character*2 line

      tolowf = line

      do i=1,2
         j = ichar(tolowf(i:i))
         if (j.ge.65.and.j.le.90) tolowf(i:i) = char(j+32)
      end do

      return
      end
