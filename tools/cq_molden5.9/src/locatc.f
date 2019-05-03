      function locatc(label,nf,itext)
      implicit double precision (a-h,p-w),integer (i-n)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y), logical (o)
      character*(*) label,itext
      dimension label(*)

      do i=1,nf
         if (label(i).eq.itext) then
            locatc = i
            return
         endif
      end do

      locatc = 0

      return
      end
