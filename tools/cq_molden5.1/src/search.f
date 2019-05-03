      subroutine search(line,name,istat)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*(*) line,name
      character*80 caps

      istat = 0
      caps  = name  
      len1  = len(name)
      call tocap(caps,len1)
1     read(iun2,'(a)',end=100,err=200) line
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

      subroutine srclit(line,name,istat)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*(*) line,name
      character*80 caps

      istat = 0
      caps  = name  
      len1  = len(name)
      call tocap(caps,len1)
1     read(iun2,'(a)',end=100,err=200) line
      if (index(line,name).ne.1.and.index(line,caps(1:len1)).ne.1)
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

      integer function icdex(line,name)
      character*(*) line,name
      character*80 caps
      character*137 cline

      il1 = len(line)
      if (il1.gt.137) il1 = 137
      il2 = len(name)
      if (il2.gt.80) il2 = 80

      cline(1:il1) = line
      caps(1:il2) = name
      call tocap(cline,il1)
      call tocap(caps,il2)

      icdex = index(cline(1:il1),caps(1:il2))

      return
      end

      subroutine leftj(line,name)
      character*(*) line,name

      name = line
      len1  = len(line)

      nsp = 0
      do i=1,len1
         if (ichar(name(i:i)).ne.32) then
            nsp = i
            goto 10
         endif
      end do

10    continue

      if (nsp.ne.0) then
         name = name(nsp:len1)
      endif

      return
      end

      subroutine searchu(line,name,istat)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*(*) line,name
      character*80 caps
      character*137 caps2

      istat = 0
      caps  = name  
      len1  = len(name)
      call tocap(caps,len1)
1     read(iun2,'(a)',end=100,err=200) line
      caps2 = line
      len2 = len(line)
      call tocap(caps2,len2)
      if (index(caps2(1:len2),caps(1:len1)).eq.0) goto 1
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
