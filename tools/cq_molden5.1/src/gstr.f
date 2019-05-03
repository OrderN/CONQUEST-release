      character*2 function gstr(i)
      character*1 sl,sh
      integer i,il,ih

      gstr = '**'
      if ( i .gt. 99 .or. i.lt.0) return

      ih = i / 10
      sh = char(ih+48)
      if (ih.eq.0) sh = char(32)
      il = i - ih*10
      sl = char(il+48)
      
      gstr = sh//sl

      return
      end      

      character*2 function ggstr(i)
      character*1 sl,sh
      integer i,il,ih

      ggstr = '**'
      if ( i .gt. 99 .or. i.lt.0) return

      ih = i / 10
      sh = char(ih+48)
      if (ih.eq.0) sh = char(48)
      il = i - ih*10
      sl = char(il+48)
      
      ggstr = sh//sl

      return
      end      

      character*2 function agstr(i,upp)
      character*1 sl,sh
      integer i,il,ih
      logical upp

      agstr = '**'
      if (i.gt.26*26.or.i.lt.0) return

      ioff = 97
      if (upp) ioff = 65
      ih = i / 26
      sh = char(ih+ioff)
      il = i - ih*26
      sl = char(il+ioff)
      
      agstr = sh//sl

      return
      end      

      character*5 function hstr(i)
      character*1 sl,sh,shh
      integer i,il,ih,ihh

      hstr = '(***)'
      if ( i .gt. 999 .or. i.lt.0) return


      ihh = i / 100
      shh = char(ihh+48)

      ih = (i - ihh*100) / 10 
      sh = char(ih+48)

      il = i - ihh*100 - ih*10
      sl = char(il+48)

      if (ihh.eq.0.and.ih.eq.0) then
          hstr = '('//sl//')  '
      elseif (ihh.eq.0) then
          hstr = '('//sh//sl//') '
      else
          hstr = '('//shh//sh//sl//')'
      endif

      return
      end      

      subroutine zerstr(inum,str,isiz,iopt)
      integer inum,i,j,isiz,iwrk,ibig,icnt
      character*(*) str

      do i=1,isiz
         str(i:i) = char(48)
      end do

      iwrk = inum
      do i=isiz,1,-1
         ibig = 10**(i-1)
         if (ibig.le.iwrk) then
            icnt = iwrk/ibig
            j = isiz-i+1
            str(j:j) = char(48+icnt)
            iwrk = iwrk - icnt*ibig
         endif
      end do
      
      if (iopt.eq.1) then
         i = 1
         do while (str(i:i).eq.char(48) )
            str(i:i) = char(32)
            i = i + 1
         end do
      endif

      return
      end      

      subroutine zzrstr(inum,str,isiz,isout)
      integer inum,i,j,isiz,iwrk,ibig,icnt,isout
      character*(*) str

      do i=1,isiz
         str(i:i) = char(48)
      end do
      str(1:1) = '('

      istrt = 0
      iwrk = inum
      j = 1
      do i=isiz,1,-1
         ibig = 10**(i-1)
         if (ibig.le.iwrk.or.istrt.eq.1) then
            istrt = 1
            j = j + 1
            icnt = iwrk/ibig
            str(j:j) = char(48+icnt)
            iwrk = iwrk - icnt*ibig
         endif
      end do
      j = j + 1
      str(j:j) = ')'
      isout = j
      
      return
      end      

