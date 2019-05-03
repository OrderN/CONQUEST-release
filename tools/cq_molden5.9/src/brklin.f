      integer function getlin(ieq)
      character*137 line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /curlin/ line
      common /mflin/  linmf

      getlin = 1

      call nxtlin(line,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 100

      linmf = linmf + 1
      do i=1,137
         if (ichar(line(i:i)).eq.9) line(i:i) = ' '
      end do

      if (ieq.eq.1.or.ieq.eq.2) then
         do i=1,137
            if (ichar(line(i:i)).eq.61) line(i:i) = ' '
         end do
      endif

      if (ieq.eq.2) then
         do i=1,137
            ii = ichar(line(i:i))
            if (ii.eq.40.or.ii.eq.41.or.ii.eq.34.or.ii.eq.39) 
     &         line(i:i) = ' '
         end do
      endif

      if (ieq.eq.3) then
         do while (.true.)
            i1 = index(line,'(')
            i2 = index(line,')')
            if (i1.gt.0.and.i2.gt.0) then
               line = line(1:i1-1)//line(i2+1:)
            else
               return
            endif
         end do
      endif

      return
100   getlin = 0
      return
      end

      subroutine setlin(str,ic)
      character*(*) str
      character*137 line
      common /curlin/ line

      line = str
      if (ic.ne.0) then
         do i=1,linlen(str)
            if (ichar(line(i:i)).eq.ic) line(i:i) = ' '
         end do
      endif

      return
      end

      integer function nxtwrd(string,strlen,itype,rtype)
c
c      string           nxtwrd = 1 
c      integer          nxtwrd = 2 
c      real             nxtwrd = 3 
c      no word          nxtwrd = 0
c
      character*(*) string
      integer itype,strlen
      double precision rtype
      double precision reada
      logical chkstr
      character*137 line
      common /curlin/ line

      nxtwrd = 0
      
      llen = linlen(line)
      if (llen.eq.0) return

      do while (line(1:1).eq.' ')
          line = line(2:)
      end do

      iend = index(line,' ')
      if (iend.eq.0) then
         iend = llen
      else
         iend = iend - 1
      endif
      if (chkstr(line,iend)) then
           nxtwrd = 1
           string = line(1:iend)
           strlen = iend
      elseif (index(line(1:iend),'.').ne.0) then
           nxtwrd = 3
           rtype = reada(line,1,iend)
      else
           nxtwrd = 2
           itype = reada(line,1,iend)
      endif

      line = line(iend+1:)

      return
      end

      integer function nxtwrz(string,strlen,itype,rtype)
c
c      string           nxtwrd = 1 
c      integer          nxtwrd = 2 
c      real             nxtwrd = 3 
c      n*int            nxtwrd = 4 
c      no word          nxtwrd = 0
c
      character*(*) string
      integer itype,strlen
      double precision rtype
      double precision reada
      logical chkstd
      character*137 line
      common /curlin/ line

      nxtwrz = 0
      
      llen = linlen(line)
      if (llen.eq.0) return

      do while (line(1:1).eq.' ')
          line = line(2:)
      end do

      iend = index(line,' ')
      if (iend.eq.0) then
         iend = llen
      else
         iend = iend - 1
      endif
      if (chkstd(line,iend)) then
           nxtwrz = 1
           string = line(1:iend)
           strlen = iend
      elseif (index(line(1:iend),'.').ne.0) then
           if (index(line(1:iend),'*').ne.0) then
              ied = index(line,'*')
              if (ied.eq.0) then
                 ied = llen
              else
                 ied = ied - 1
              endif
              itype = reada(line,1,ied)
              nxtwrz = 4
           else
              nxtwrz = 3
              rtype = reada(line,1,iend)
           endif
      else
           nxtwrz = 2
           itype = reada(line,1,iend)
      endif

      line = line(iend+1:)

      return
      end

      integer function nxtwrx(string,strlen,itype,rtype)
c
c      string           nxtwrd = 1 
c      integer          nxtwrd = 2 
c      real             nxtwrd = 3 
c      no word          nxtwrd = 0
c
      character*(*) string
      integer itype,strlen
      double precision rtype
      double precision reada
      logical chkstr
      character*137 line
      common /curlin/ line

      nxtwrx = 0
      
      nine  = ichar('9')
      izero = ichar('0')

      llen = linlen(line)
      if (llen.eq.0) return

      do while (line(1:1).eq.' ')
          line = line(2:)
      end do

      if (llen.gt.3) then
          do i=1,llen-2
              ii = ichar(line(i+1:i+1))
              if (line(i:i).eq.'('.and.line(i+2:i+2).eq.')'.and.
     &           (ii.ge.izero.and.ii.le.nine)) then
                 line = line(i-1:)//line(i+3:)
              endif
          end do
      endif

      llen = linlen(line)
      if (llen.eq.0) return

      iend = index(line,' ')
      if (iend.eq.0) then
         iend = llen
      else
         iend = iend - 1
      endif

      if (chkstr(line,iend)) then
           nxtwrx = 1
           string = line(1:iend)
           strlen = iend
      elseif (index(line(1:iend),'.').ne.0) then
           nxtwrx = 3
           rtype = reada(line,1,iend)
      else
           nxtwrx = 2
           itype = reada(line,1,iend)
      endif

      line = line(iend+1:)

      return
      end

      logical function chkstr(line,iend)
      character*(*) line
      chkstr = .false.

      ie    = ichar('e')
      iee   = ichar('E')
      id    = ichar('d')
      idd   = ichar('D')
      nine  = ichar('9')
      izero = ichar('0')
      minus = ichar('-')
      iplus = ichar('+')
      idot  = ichar('.')
      icomma = ichar(',')
      islash  = ichar('/')

      ihase = 0
      idig = 0
      do i=1,iend
         n = ichar(line(i:i))
         if ((n.eq.ie.or.n.eq.iee.or.n.eq.id.or.n.eq.idd)
     &      .and.ihase.eq.0.and.idig.eq.1) then
             n = izero
             ihase = 1
         endif
         if (n.lt.iplus.or.n.gt.nine.or.n.eq.islash
     &       .or.n.eq.icomma) goto 100
         idig = 1
      end do

      n = ichar(line(1:1))
      n2 = ichar(line(2:2))
      if (iend.eq.1) then
         if (n.eq.minus) goto 100
         if (n.eq.iplus) goto 100
         if (n.eq.ie.or.n.eq.iee) goto 100
         if (n.eq.id.or.n.eq.idd) goto 100
      elseif (iend.gt.1) then
         if (n.eq.minus.and.n2.eq.minus) goto 100
      endif

      return
100   chkstr = .true.
      return
      end

      logical function chkstd(line,iend)
      character*(*) line
      chkstd = .false.

      ie    = ichar('e')
      iee   = ichar('E')
      id    = ichar('d')
      idd   = ichar('D')
      nine  = ichar('9')
      izero = ichar('0')
      minus = ichar('-')
      iplus = ichar('+')
      idot  = ichar('.')
      icomma = ichar(',')
      islash = ichar('/')
      istar  = ichar('*')

      ihase = 0
      idig = 0
      do i=1,iend
         n = ichar(line(i:i))
         if ((n.eq.ie.or.n.eq.iee.or.n.eq.id.or.n.eq.idd)
     &      .and.ihase.eq.0.and.idig.eq.1) then
             n = izero
             ihase = 1
         endif
         if ((n.lt.iplus.or.n.gt.nine.or.n.eq.islash
     &       .or.n.eq.icomma).and.n.ne.istar) goto 100
         idig = 1
      end do

      n = ichar(line(1:1))
      n2 = ichar(line(2:2))
      if (iend.eq.1) then
         if (n.eq.minus) goto 100
         if (n.eq.iplus) goto 100
         if (n.eq.ie.or.n.eq.iee) goto 100
         if (n.eq.id.or.n.eq.idd) goto 100
      elseif (iend.gt.1) then
         if (n.eq.minus.and.n2.eq.minus) goto 100
      endif

      return
100   chkstd = .true.
      return
      end

      integer function linlen(line)
      character*(*) line
      integer i,n

      linlen = 0

      do i=len(line),1,-1
         n = ichar(line(i:i))
         if (n.gt.32.and.n.le.126) goto 100
      end do

      return
100   linlen = i
      return
      end

      logical function dat3ln(lin)
      integer i,itype,ktype,nstr
      double precision rtype
      character*(*) lin
      character*137 str
      character*137 line
      common /curlin/ line

      dat3ln = .true.

      line = lin

      do i=1,3
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.ne.3.and.ktype.ne.2) goto 100
      end do

      return
100   dat3ln = .false.
      return
      end

      logical function datlin(line)
      character*(*) line

      datlin = .true.

      do i=1,linlen(line)
         n = ichar(line(i:i))
         if ((n.lt.43.or.n.gt.57).and.n.ne.32.and.n.ne.68.
     &       and.n.ne.100.and.n.ne.69.and.n.ne.101) goto 100
      end do

      return
100   datlin = .false.
      return
      end

      logical function gnreal(r,n,doget)
      implicit double precision (a-h,o-z)
      character*137 line,str
      common /curlin/ line
      integer getlin
      logical doget
      dimension r(*)

      gnreal = .true.

      if (doget) then
         if (getlin(0).ne.1) gnreal = .false.
      endif

      if (gnreal) then
          do i=1,n
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.3) then
                r(i) = rtype
             elseif (ktype.eq.2) then
                r(i) = dble(itype)
             else
                gnreal = .false.
             endif
          end do
      endif

      return
      end
    
      logical function gnint(iarr,n,doget)
      implicit double precision (a-h,o-z)
      character*137 line,str
      common /curlin/ line
      integer getlin
      logical doget
      dimension iarr(*)

      gnint = .true.

      if (doget) then
         if (getlin(0).ne.1) gnint = .false.
      endif


      if (gnint) then
          do i=1,n
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.2) then
                iarr(i) = itype
             else
                gnint = .false.
             endif
          end do
      endif

      return
      end
    
      subroutine lsparm(str,l)
      character*(*) str

      l = len(str)
      do i=1,l
         if (str(i:i).ne.' ') goto 10
      end do
10    str = str(i:)
      l = linlen(str)
      return
      end

      subroutine spatrm(str,l)
      character*(*) str

      j = 1
      l = len(str)
      do while (j.le.l) 
         if (str(j:j).eq.' ') then
           if (l.eq.1) then
              return
           else
              if (j.ne.l) then
                 if (j.eq.1) then
                    str(1:l-1) = str(2:l)
                    str(l:l) = ' '
                 else
c                    str(j:l-1) = str(j:j-1)//str(j+1:l)
                    str(j:l-1) = str(j+1:l)
                    str(l:l) = ' '
                 endif
              endif
              l = l - 1
           endif
         else
           j = j + 1
         endif
      end do

c      if (l.lt.len(str)) str(l+1:l+1) = char(0)

      return
      end

      integer function krnd(r)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)

      krnd = int(r)
      if (r-dfloat(krnd).ge.0.5d0) krnd = krnd + 1

      return
      end

      subroutine rmnull(line)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      character*137 line

      ii = 0
      do while (ii.le.137)
         ii = ii + 1
         jj = ichar(line(ii:ii))
         if (jj.eq.0) then
            do kk=ii,136
               line(kk:kk) = line(kk+1:kk+1)
            end do
            if (ii.eq.138) return
            ii = ii - 1
         endif
         if (jj.eq.10.or.jj.eq.13) return
      end do

      return
      end

      subroutine rwfile()
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      rewind iun2

      return
      end

      subroutine bcfile()
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      backspace iun2

      return
      end

      subroutine nxline(line,istat)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      character*137 line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      istat = 0

      read(iun2,'(a)',end=100,err=200) line

      return
100   istat = 1
      return
200   istat = 2
      return
      end

