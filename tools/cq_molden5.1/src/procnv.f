      subroutine procnv
      parameter (MAXITER=1000)
      implicit double precision (a-h,o-z)
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /convrg/ convg1(MAXITER),convg2(MAXITER),cnvmax,cnvmin,
     &jstrt1,jend1,jstrt2,jend2,icvav1,icvav2
      logical dozme
      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme

      if (iftyp.eq.2.or.iftyp.eq.3) then
          if (iftyp.eq.3) then
             call cnvgus
          else
             call cnvgam
          endif
      elseif (iftyp.eq.4) then
          call cnvgau
      elseif (iftyp.eq.8) then
          call cnvqcm
      elseif (iftyp.eq.9) then
          call cnvorc
      elseif (iftyp.eq.7) then
          if (irtype.le.3) call cnvcpmd
      else
          if (iftyp.eq.5) then
             icvav1 = 0
             icvav2 = 0
          endif
      endif

      if (.not.(iftyp.eq.1.or.ipdbon.eq.1.or.iftyp.eq.6.or.
     &          isbin.eq.1.or.iftyp.eq.5)) 
     &     call mmcnv

      call parptr(4,convg1,fdum,idum)

      return
      end
