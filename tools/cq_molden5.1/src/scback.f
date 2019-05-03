      subroutine scback(line,name,istat)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*(*) line,name
      character*137 caps

      istat = 0
      caps  = name
      len1  = len(name)
      call tocap(caps,len1)
1     backspace iun2
      read(iun2,'(a)',end=100,err=200) line
      backspace iun2
      if (index(line,name).eq.0.and.index(line,caps(1:len1)).eq.0)
     &  go to 1
c
      istat=1
      return
100   rewind iun2
      return
200   continue
      return
c
c
      end
