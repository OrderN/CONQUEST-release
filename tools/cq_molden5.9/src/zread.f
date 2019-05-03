      subroutine zread(nz,ianz,iz,bl,alpha,beta)
c
c
c
c      z-matrix reading routine.
c      It expects to be positioned after the string point #
c      adapted by g.schaftenaar from sprint gamess
c
c
      implicit double precision (a-h,o-z)
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character * 2 iel,jel
      character *137 string
      dimension ianz(*), bl(*), alpha(*), beta(*)
      dimension iz(4,*)
      dimension iel(100)
c
c      note that lower case letters are used in the atomic symbols
c      and in some of the format statements in hollereith strings.
c
      data iel/'bq',
     $         'h ', 'he',
     $         'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne',
     $         'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar',
     $         'k ', 'ca',
     $                     'sc', 'ti', 'v ', 'cr', 'mn',
     $                     'fe', 'co', 'ni', 'cu', 'zn',
     $                     'ga', 'ge', 'as', 'se', 'br', 'kr',
     $ 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     $ 'in','sn','sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     $ 'pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     $ 'ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po',
     $ 'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm',
     $ 'bk','cf','x '  /
 2110 format(1x,i2,2x,i2,3x,a2)
 2111 format(1x,i3,2x,i3,3x,a2)
 2210 format(1x,i2,2x,i2,3x,a2,3x,i2,2x,f9.6,2x,i3)
 2211 format(1x,i3,2x,i3,3x,a2,2x,i3,2x,f9.6,2x,i4)
 2310 format(1x,i2,2x,i2,3x,a2,3x,i2,2x,f9.6,2x,i3,2x,
     $          i2,1x,f8.3,2x,i3)
 2311 format(1x,i3,2x,i3,3x,a2,2x,i3,2x,f9.6,2x,i4,2x,
     $          i3,1x,f8.3,2x,i4)
 2410 format(1x,i2,2x,i2,3x,a2,3x,i2,2x,f9.6,2x,i3,2x,
     $          i2,1x,f8.3,2x,i3,2x,i2,1x,f8.3,2x,i3,2x,i2)
 2411 format(1x,i3,2x,i3,3x,a2,2x,i3,2x,f9.6,2x,i4,2x,
     $          i3,1x,f8.3,2x,i4,2x,i3,1x,f8.3,2x,i4,2x,i2)
c
c
c
c     print the heading.
c
c
c     first card of z-matrix.
c
      jmode = 0
      call nxtlin(string,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 200

      if (string(4:4).ne.' ') jmode = 1
      if ( index(string(2:5),'====') .ne. 0 ) return
      if (jmode.eq.0) then
         read(string,2110,err=200) icard,icent,jel
      else
         read(string,2111,err=200) icard,icent,jel
      endif

      do i=1,100
         if ( jel .eq. iel(i) ) ianz(1) = i - 1
      end do

      nz = 1    
  
c
c     second card of the z-matrix.
c
      call nxtlin(string,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 200
      if ( index(string(2:5),'====') .ne. 0 ) return
      if (jmode.eq.0) then
         read(string,2210,err=200) icard,icent,jel,iz(1,2),pbl,np1
      else
         read(string,2211,err=200) icard,icent,jel,iz(1,2),pbl,np1
      endif

      do i=1,100
         if ( jel .eq. iel(i) ) ianz(2) = i - 1
      end do

      bl(2) = pbl
      nz = 2

c
c     third card.
c
      call nxtlin(string,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 200

      if ( index(string(2:5),'====') .ne. 0 ) return
      if (jmode.eq.0) then
         read(string,2310,err=200) 
     &   icard,icent,jel,iz(1,3),pbl,np1,iz(2,3),pa,np2
      else
         read(string,2311,err=200) 
     &   icard,icent,jel,iz(1,3),pbl,np1,iz(2,3),pa,np2
      endif

      do i=1,100
         if ( jel .eq. iel(i) ) ianz(3) = i - 1
      end do

      bl(3) = pbl
      alpha(3) = pa
      nz=3
c
c     cards 4 through nz.
c
      do while (.true.)
         call nxtlin(string,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 200

         if ( index(string(2:5),'====') .ne. 0 ) return
         nz = nz + 1
         if (jmode.eq.0) then
            read(string,2410,err=200) icard,icent,jel,iz(1,nz),pbl,np1,
     &                        iz(2,nz),pa,np2,iz(3,nz),pb,np3,iz(4,nz)
         else
            read(string,2411,err=200) icard,icent,jel,iz(1,nz),pbl,np1,
     &                        iz(2,nz),pa,np2,iz(3,nz),pb,np3,iz(4,nz)
         endif

         do i=1,100
            if ( jel .eq. iel(i) ) ianz(nz) = i - 1
         end do

         bl(icard) = pbl
         alpha(icard) = pa
         beta(icard) = pb

      end do

200   call inferr('error in z-matrix !',0)
      call inferr(string,0)

      return
      end

      logical function zreadg(nz,ianz,iz,bl,alpha,beta)
c
c
c
c      z-matrix reading routine for gaussian output.
c      It expects to be positioned after the string point #
c      adapted by g.schaftenaar from sprint gamess
c
c
      implicit double precision (a-h,o-z)
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character *2 iel,jel
      character *100 string
      character *2 tocapf
      dimension ianz(*), bl(*), alpha(*), beta(*)
      dimension iz(4,*)
      dimension iel(100)
c
c      note that lower case letters are used in the atomic symbols
c      and in some of the format statements in hollereith strings.
c
      data iel/'bq',
     $         'h ', 'he',
     $         'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne',
     $         'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar',
     $         'k ', 'ca',
     $                     'sc', 'ti', 'v ', 'cr', 'mn',
     $                     'fe', 'co', 'ni', 'cu', 'zn',
     $                     'ga', 'ge', 'as', 'se', 'br', 'kr',
     $ 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     $ 'in','sn','sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     $ 'pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     $ 'ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po',
     $ 'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm',
     $ 'bk','cf','x '  /
 2110 format(1x,i3,1x,i3,2x,a2)
 2210 format(1x,i3,1x,i3,2x,a2,2x,i3,f11.6)
 2310 format(1x,i3,1x,i3,2x,a2,2x,i3,f11.6,6x,i3,1x,f8.3)
 2410 format(1x,i3,1x,i3,2x,a2,2x,i3,f11.6,6x,i3,1x,f8.3,
     $          6x,i3,1x,f8.3,6x,i3)
 2111 format(1x,i6,1x,i6,2x,a2)
 2211 format(1x,i6,1x,i6,2x,a2,2x,i6,f11.6)
 2311 format(1x,i6,1x,i6,2x,a2,2x,i6,f11.6,9x,i6,1x,f8.3)
 2411 format(1x,i6,1x,i6,2x,a2,2x,i6,f11.6,9x,i6,1x,f8.3,
     $          9x,i6,1x,f8.3,9x,i6)
      zreadg = .true.
      iform = 1

c     first card of z-matrix.

      call nxtlin(string,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 200

      il = linlen(string)
      if (il.eq.12.or.il.eq.11) then
         iform = 1
      else
         iform = 2
      endif
      if ( index(string(2:5),'----') .ne. 0 ) return

      if (iform.eq.1) then
         read(string,2110,err=200) icard,icent,jel
      else
         read(string,2111,err=200) icard,icent,jel
      endif

      do i=1,100
         if (tocapf(jel).eq.tocapf(iel(i))) ianz(1) = i - 1
      end do
      nz = 1    
  

c     second card of the z-matrix.

      call nxtlin(string,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 200

      if ( index(string(2:5),'----') .ne. 0 ) return
      if (iform.eq.1) then
         read(string,2210,err=200) icard,icent,jel,iz(1,2),pbl
      else
         read(string,2211,err=200) icard,icent,jel,iz(1,2),pbl
      endif

      do i=1,100
         if (tocapf(jel).eq.tocapf(iel(i))) ianz(2) = i - 1
      end do

      if (iz(1,2).le.0) then
          zreadg = .false.
          return
      endif

      bl(2) = pbl
      nz = 2


c     third card.

      call nxtlin(string,jstat)
      if (jstat.eq.1.or.jstat.eq.2) goto 200

      if ( index(string(2:5),'----') .ne. 0 ) return

      if (iform.eq.1) then
         read(string,2310,err=200) icard,icent,jel,iz(1,3),pbl,
     &                             iz(2,3),pa
      else
         read(string,2311,err=200) icard,icent,jel,iz(1,3),pbl,
     &                             iz(2,3),pa
      endif

      do i=1,100
         if (tocapf(jel).eq.tocapf(iel(i))) ianz(3) = i - 1
      end do

      bl(3) = pbl
      alpha(3) = pa
      nz=3

c     cards 4 through nz.

      do while (.true.)
         call nxtlin(string,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 200

         if ( index(string(2:5),'----') .ne. 0 ) return
         nz = nz + 1
         if (iform.eq.1) then
            read(string,2410,err=200) icard,icent,jel,iz(1,nz),pbl,
     &                    iz(2,nz),pa,iz(3,nz),pb,iz(4,nz)
         else
            read(string,2411,err=200) icard,icent,jel,iz(1,nz),pbl,
     &                    iz(2,nz),pa,iz(3,nz),pb,iz(4,nz)
         endif

         do i=1,100
            if (tocapf(jel).eq.tocapf(iel(i))) ianz(nz) = i - 1
         end do

         bl(icard) = pbl
         alpha(icard) = pa
         beta(icard) = pb
      end do

200   call inferr('error in z-matrix !',0)
      call inferr(string,0)

      return
      end
