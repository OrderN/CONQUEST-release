      subroutine gmmcnv(formax,forrms,dismax,disrms,epoints,isav)
      implicit double precision (a-h,o-z)
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)
      
      fgmax = -1000000.d0
      fgmin =  1000000.d0
      dgmax = -1000000.d0
      dgmin =  1000000.d0
      enmax = -1000000.d0
      enmin =  1000000.d0

      do i=1,ngeoms
        if(isav(i).eq.1) then
         if (formax(i).lt.fgmin) fgmin = formax(i)
         if (formax(i).gt.fgmax) fgmax = formax(i)
         if (forrms(i).lt.fgmin) fgmin = forrms(i)
         if (forrms(i).gt.fgmax) fgmax = forrms(i)
         if (dismax(i).lt.dgmin) dgmin = dismax(i)
         if (dismax(i).gt.dgmax) dgmax = dismax(i)
         if (disrms(i).lt.dgmin) dgmin = disrms(i)
         if (disrms(i).gt.dgmax) dgmax = disrms(i)
        endif
      end do

      do i=1,nepnts
         if (epoints(i).lt.enmin) enmin = epoints(i)
         if (epoints(i).gt.enmax) enmax = epoints(i)
      end do

      return
      end
