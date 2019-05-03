      subroutine scback(line,name,istat)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*(*) line,name
      character*137 caps

      istat = 0
      caps  = name
      len1  = len(name)
      call tocap(caps,len1)
1     call bckfil
      call nxtlin(line,jstat)
      if (jstat.eq.1) goto 100
      if (jstat.eq.2) goto 200
      call bckfil
      if (index(line,name).eq.0.and.index(line,caps(1:len1)).eq.0)
     &  go to 1
c
      istat=1
      return
100   call rewfil
      return
200   continue
      return
c
c
      end
