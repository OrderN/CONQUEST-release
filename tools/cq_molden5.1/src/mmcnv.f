      subroutine mmcnv
      parameter (MAXITER=1000)
      implicit double precision (a-h,o-z)
      common /convrg/ convg1(MAXITER),convg2(MAXITER),cnvmax,cnvmin,
     &jstrt1,jend1,jstrt2,jend2,icvav1,icvav2
      
      cnvmax = -1000000.d0
      cnvmin =  1000000.d0

      if (icvav1.eq.1) then
         do i=jstrt1,jend1
            if (convg1(i).lt.cnvmin) cnvmin = convg1(i)
            if (convg1(i).gt.cnvmax) cnvmax = convg1(i)
         end do
      endif

      if (icvav2.eq.1) then
         do i=jstrt2,jend2
            if (convg2(i).lt.cnvmin) cnvmin = convg2(i)
            if (convg2(i).gt.cnvmax) cnvmax = convg2(i)
         end do
      endif

      return
      end
