      subroutine searchd(line,name1,name2,istat)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*(*) line,name1,name2
      character*137 caps1, caps2

      istat = 0
      caps1 = name1
      caps2 = name2
      len1  = len(name1)
      len2  = len(name2)
      call tocap(caps1,len1)
      call tocap(caps2,len2)
1     call nxtlin(line,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 100
      if (index(line,name1).eq.0.and.index(line,name2).eq.0
     &.and.index(line,caps1(1:len1)).eq.0
     &.and.index(line,caps2(1:len2)).eq.0) go to 1
c
      istat=1
      return
c
100   call rewfil
c
      return
      end 

      subroutine seardu(line,name1,name2,istat)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*(*) line,name1,name2
      character*137 caps1, caps2, caps3

      istat = 0
      len1  = min(137,len(name1))
      len2  = min(137,len(name2))
      caps1 = name1(1:len1)
      caps2 = name2(1:len2)
      call tocap(caps1,len1)
      call tocap(caps2,len2)

1     call nxtlin(line,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 100
      len3  = min(137,len(line))
      caps3 = line(1:len3)
      call tocap(caps3,len3)
      if (index(caps3,caps1(1:len1)).eq.0.and.
     &    index(caps3,caps2(1:len2)).eq.0) goto 1
 
      istat = 1
      return
 
100   call rewfil
      return
 
      end 

      subroutine seartu(line,name1,name2,name3,istat)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*(*) line,name1,name2,name3
      character*137 caps1, caps2, caps3, capst

      istat = 0
      len1  = min(137,len(name1))
      len2  = min(137,len(name2))
      len3  = min(137,len(name3))
      caps1 = name1(1:len1)
      caps2 = name2(1:len2)
      caps3 = name3(1:len3)
      call tocap(caps1,len1)
      call tocap(caps2,len2)
      call tocap(caps3,len3)

1     call nxtlin(line,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 100
      lent  = min(137,len(line))
      capst = line(1:lent)
      call tocap(capst,lent)
      if (index(capst,caps1(1:len1)).eq.0.and.
     &    index(capst,caps2(1:len2)).eq.0.and.
     &    index(capst,caps3(1:len3)).eq.0) goto 1
 
      istat = 1
      return
 
100   call rewfil
      return
 
      end 

      subroutine searcht(line,name1,name2,name3,istat)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*(*) line,name1,name2,name3
      character*137 caps1, caps2, caps3

      istat = 0
      caps1 = name1
      caps2 = name2
      caps3 = name3
      len1  = len(name1)
      len2  = len(name2)
      len3  = len(name3)
      call tocap(caps1,len1)
      call tocap(caps2,len2)
      call tocap(caps3,len3)
1     call nxtlin(line,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 100
      if (index(line,name1).eq.0.and.index(line,name2).eq.0
     &.and.index(line,name3).eq.0
     &.and.index(line,caps1(1:len1)).eq.0
     &.and.index(line,caps2(1:len2)).eq.0
     &.and.index(line,caps3(1:len3)).eq.0) go to 1
c
      istat=1
      return
c
100   call rewfil
      return
c
      end 

      subroutine searchq(line,name1,name2,name3,name4,istat)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*(*) line,name1,name2,name3,name4
      character*137 caps1, caps2, caps3, caps4

      istat = 0
      caps1 = name1
      caps2 = name2
      caps3 = name3
      caps4 = name4
      len1  = len(name1)
      len2  = len(name2)
      len3  = len(name3)
      len4  = len(name4)
      call tocap(caps1,len1)
      call tocap(caps2,len2)
      call tocap(caps3,len3)
      call tocap(caps4,len4)
1     call nxtlin(line,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 100
      if (index(line,name1).eq.0.and.index(line,name2).eq.0
     &.and.index(line,name3).eq.0.and.index(line,name4).eq.0
     &.and.index(line,caps1(1:len1)).eq.0
     &.and.index(line,caps2(1:len2)).eq.0
     &.and.index(line,caps3(1:len3)).eq.0
     &.and.index(line,caps4(1:len4)).eq.0) go to 1
c
      istat=1
      return
c
100   call rewfil
      return
c
      end 

      subroutine searchv(line,name1,name2,name3,name4,name5,istat)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*(*) line,name1,name2,name3,name4,name5
      character*137 caps1, caps2, caps3, caps4, caps5

      istat = 0
      caps1 = name1
      caps2 = name2
      caps3 = name3
      caps4 = name4
      caps5 = name5
      len1  = len(name1)
      len2  = len(name2)
      len3  = len(name3)
      len4  = len(name4)
      len5  = len(name5)
      call tocap(caps1,len1)
      call tocap(caps2,len2)
      call tocap(caps3,len3)
      call tocap(caps4,len4)
      call tocap(caps5,len5)
1     call nxtlin(line,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 100
      if (index(line,name1).eq.0.and.index(line,name2).eq.0
     &.and.index(line,name3).eq.0.and.index(line,name4).eq.0
     &.and.index(line,name5).eq.0.and.index(line,caps1(1:len1)).eq.0
     &.and.index(line,caps2(1:len2)).eq.0
     &.and.index(line,caps3(1:len3)).eq.0
     &.and.index(line,caps4(1:len4)).eq.0
     &.and.index(line,caps5(1:len5)).eq.0) go to 1
c
      istat=1
      return
c
100   call rewfil
      return
c
      end 
