      subroutine getreal(line,stri,vari)
      implicit double precision (a-h,o-z), integer ( i-n)
      character*137 line
      character*4   stri
      indx = index(line,stri)
      if (indx.ne.0) then
         vari = reada(line,indx+6,len(line))
      endif

      return
      end
