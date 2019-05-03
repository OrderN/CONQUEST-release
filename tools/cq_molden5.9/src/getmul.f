      subroutine getmul
      implicit double precision (a-h,o-z), integer ( i-n)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*137 line
c     Dirk Huckriede got Bug out 25 Nov 93

      call search(line,'distributed multipole analysis module',
     &            istat)
      do i=1,13
         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 100

      end do
3        call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 100

         if (index(line,'site').ne.0) call site(line)
         if (index(line,'total multipoles referred to origin').ne.0) 
     &   goto 2
         goto 3
2     continue
      return
100   call inferr('DMA: Premature End Of File !',1)
      end
