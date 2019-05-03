      subroutine fcij(itype,istrt,isp,idf,cij)
      implicit double precision (a-h,o-z)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      integer shella,shelln,shellt,shellc,shladf,aos
      dimension cij(35)
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp

      if (itype-3.lt.0) then
         cij(1) = c1(isp)
         if (itype.gt.0) then
            cij(2) = c2(isp)
            cij(3) = c2(isp)
            cij(4) = c2(isp)
            if (itype-1.gt.0) then
               iidf = idf + isp - istrt
               cij( 5) = c3(iidf)
               cij( 6) = c3(iidf)
               cij( 7) = c3(iidf)
               cij( 8) = c3(iidf)*dsqrt(3.0d0)
               cij( 9) = cij(8)
               cij(10) = cij(8)
            endif
         endif
      elseif (itype.eq.3) then
         iidf = idf + isp - istrt
         cij(11) = c4(iidf)
         cij(12) = c4(iidf)
         cij(13) = c4(iidf)
         cij(14) = c4(iidf)*dsqrt(5.0d0)
         cij(15) = cij(14)
         cij(16) = cij(14)
         cij(17) = cij(14)
         cij(18) = cij(14)
         cij(19) = cij(14)
         cij(20) = c4(iidf)*dsqrt(15.0d0)
      else
         iidf = idf + isp - istrt
         cij(21) = c5(iidf)
         cij(22) = c5(iidf)
         cij(23) = c5(iidf)
         cij(24) = c5(iidf)*dsqrt(7.0d0)
         cij(25) = cij(24)
         cij(26) = cij(24)
         cij(27) = cij(24)
         cij(28) = cij(24)
         cij(29) = cij(24)
         cij(30) = c5(iidf)*dsqrt(35.0d0/3.0d0)
         cij(31) = cij(30)
         cij(32) = cij(30)
         cij(33) = c5(iidf)*dsqrt(35.0d0)
         cij(34) = cij(33)
         cij(35) = cij(33)
      endif
 
      return
      end
