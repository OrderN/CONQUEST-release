      subroutine crprod(a,b,c)
c
c       calculates cross product:   a x b = c
c                                   -   -   -
      implicit real (a-h,o-z)
      dimension a(3),b(3),c(3)
c
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)

      return
      end
      subroutine crpros(a,b,c)
c
c       calculates cross product:   a x b = c
c                                   -   -   -
      implicit real (a-h,o-z)
      dimension a(3),b(3),c(3)
c
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)

      return
      end
      subroutine vsc1(a,scale,tol)
      implicit real (a-h,o-z)
      dimension a(3)

      rlen = vlen(a)

      if (rlen.gt.tol) then
         do i=1,3
            a(i) = a(i)*scale/rlen
         end do
      endif

      return
      end

      subroutine vscal(a,n,scale)
      implicit real (a-h,o-z), integer (i-n)
      dimension a(n)


      do i=1,n
         a(i) = a(i)*scale
      end do

      return
      end

      subroutine vcpy(a,b,n)
      implicit real (a-h,o-z), integer (i-n)
      dimension a(n),b(n)

      do i=1,n
         a(i) = b(i)
      end do

      return
      end

      subroutine vadd(a,b,c,n)
      implicit real (a-h,o-z), integer (i-n)
      dimension a(n),b(n),c(n)

      do i=1,n
         c(i) = a(i) + b(i)
      end do

      return
      end

      subroutine vsub(a,b,c,n)
      implicit real (a-h,o-z), integer (i-n)
      dimension a(n),b(n),c(n)

      do i=1,n
         c(i) = a(i) - b(i)
      end do

      return
      end

      subroutine vseti(ia,n,ival)
      implicit real (a-h,o-z), integer (i-n)
      dimension ia(n)

      do i=1,n
         ia(i) = ival
      end do

      return
      end

      subroutine vsetr(a,n,rval)
      implicit real (a-h,o-z), integer (i-n)
      dimension a(n)

      do i=1,n
         a(i) = rval
      end do

      return
      end

      real function vlen(a)
      implicit real (a-h,o-z)
      dimension a(3)

      tot = a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
      vlen = 0.0e0
      if (tot.gt.0.0e0) vlen = sqrt(tot)

      return
      end
      subroutine impsc(a,b,c)
      implicit real (a-h,o-z)
      dimension a(3),b(3)

      rimp = 0.0e0
    
      do i=1,3
         rimp = rimp + a(i)*b(i)
      end do

      al = vlen(a)
      bl = vlen(b)

      if (al.gt.0.0e0.and.bl.gt.0.0e0) then
         c = rimp/(al*bl)
      else
         c = 0.0e0
      endif

      return
      end

      subroutine timpsc(a,b,c)
      implicit real (a-h,o-z)
      dimension a(3),b(3)

      c = 0.0e0
    
      do i=1,3
         c = c + a(i)*b(i)
      end do

      return
      end

      character*2 function tocapf(line)
      character*2 line

      tocapf = line

      do i=1,2
         j = ichar(tocapf(i:i))
         if (j.ge.97.and.j.le.122) tocapf(i:i) = char(j-32)
      end do

      return
      end

      character*2 function tolowf(line)
      character*2 line

      tolowf = line

      do i=1,2
         j = ichar(tolowf(i:i))
         if (j.ge.65.and.j.le.90) tolowf(i:i) = char(j+32)
      end do

      return
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

      subroutine tocap(line,ncars)
      integer ncars
      character*(*) line

      do i=1,ncars
         j = ichar(line(i:i))
         if (j.ge.97.and.j.le.122) line(i:i) = char(j-32)
      end do

      return
      end

      subroutine zerstr(inum,str,isiz)
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
      
      return
      end      

      subroutine fndchr(str,ilen,chr,ic)
      implicit real (a-h,o-z), integer (i-n)
      character*(*) str
      character*1 chr
      
      ic = 0
      do i=ilen,1,-1
         if (str(i:i).eq.chr) then
            ic = i
            return
         endif
      end do

      return
      end      

       logical function gargpl(fl,n,strng1,strng2)
       integer n,fln
       character*(*) fl
       character*(*) strng1,strng2

       gargpl = .false.
       
       fln = len(fl)
       if (strng1(1:fln).ne.fl) return 
       if (strng1(fln+1:fln+1).eq.' ') then
           n = n + 1
           call getarg(n,strng2)
           if (strng2(1:1).eq.'-'.or.strng2(1:1).eq.' ') return
           strng2 = strng2(1:index(strng2,' ')-1)
       else
           strng2 = strng1(fln+1:index(strng1,' ')-1)
       endif
       gargpl = .true.

       return
       end
