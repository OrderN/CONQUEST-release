      subroutine settc(cnam,i,j)
      character*3 cnam,cstr(7)
      character*80 tmps
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension iang(7),icon(7)
      data cstr/'S  ','P  ','D  ','F  ','SP ','SPD','G  '/
      data iang/0,1,2,3,1,2,4/
      data icon/0,1,2,3,0,0,4/

      do k=1,7
         if (cnam.eq.cstr(k)) goto 20
      end do
      tmps = 'Unrecognised type of shell >>'//cnam//'<<'
      call inferr(tmps,1)
20    i = iang(k)
      j = icon(k)

      return
      end
