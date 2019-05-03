      function ifblen(string)
      implicit double precision (a-h,o-z)
      character*(*) string

      do i=len(string),1,-1
         ifblen = i
         if (string(i:i).ne.' ') return
      end do

      ifblen = 0

      return
      end
