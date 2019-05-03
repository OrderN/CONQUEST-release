      logical function obin(line)
      implicit double precision (a-h,o-z)
      character line*(*)

      obin = .false.
      jmin = 256
      jmax = 0
      nbin = 0
      ncr = 0

      do j=1,len(line)
         k = ichar(line(j:j))
         if (k.gt.jmax.and.k.ne.9.and.k.ne.27) jmax = k
         if (k.lt.jmin.and.k.ne.9.and.k.ne.27) jmin = k
         if (k.gt.126.or.k.lt.32) nbin = nbin + 1
         if (k.eq.13) ncr = ncr + 1
      end do
           
      if (jmax.gt.126.or.jmin.lt.32) then
         if (ncr.ne.nbin) obin = .true.
      endif

      return
      end

      logical function odos(fniun)
      implicit double precision (a-h,o-z)
      character*(*) fniun
      integer*1 k ,kr
      dimension kr(800)

      odos = .false.

      close(48)
      open(unit=48,form='unformatted',file=fniun
     &     ,status='old',err=1000)

      rewind(48)

      jmin = 256
      jmax = 0
      nbin = 0
      ncr = 0
      ncrlf = 0

      do i=1,800
         kr(i) = 0
      end do

      read(48,end=10,err=100) kr

10    continue

      do j=1,800
         k = kr(j)
         if (k.gt.jmax.and.k.ne.9.and.k.ne.27) jmax = k
         if (k.lt.jmin.and.k.ne.9.and.k.ne.27) jmin = k
         if (k.gt.126.or.k.lt.32) nbin = nbin + 1
         if (k.eq.13) ncr = ncr + 1
         if (k.eq.10.and.ko.eq.1) ncrlf = ncrlf + 1
         ko = 0
         if (k.eq.13) ko = 1
      end do
           
100   close(48)
      if (jmax.gt.126.or.jmin.lt.32) then
         if (ncrlf.gt.0) odos = .true.
      endif

1000  continue
      return
      end
