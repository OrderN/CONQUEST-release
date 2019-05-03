      subroutine geomop(formax,forrms,dismax,disrms,epoints,isav,
     &                  coo,ianz)
      implicit double precision (a-h,o-z)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      common /athlp/ iatoms, mxnat
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)
      dimension coo(3,*),ianz(*)

      rewind iun2

      igcvav = 1
      ifmxav = 0
      ifrmav = 0
      idmxav = 0
      idrmav = 0
      ieav = 1
      nepnts = 0
      ngeoms = 0
      dmaxt = 0.0d0
      fmaxt = 0.0d0
      drmst = 0.0d0
      frmst = 0.0d0
      enmax = -1000000.0d0
      enmin = 1000000.0d0
      do i=1,mxpnt
         epoints(i) = 0.0d0
      end do
      ihasep = 0

      idomopf = 1
      do while (.true.)
         call getmop(iatoms,heat,0,idomopf,istat)
         if (istat.eq.0) goto 100
         idomopf = 0
         if (ngeoms.eq.mxpnt) then
             call inferr('exceeded MAXPNT !',1)
             goto 100
         endif
         ngeoms = ngeoms + 1
         isav(ngeoms) = 1
         if (heat.ne.0.0d0) then
            ihasep = 1
            epoints(ngeoms) = heat
            if (heat.gt.enmax) enmax = heat
            if (heat.lt.enmin) enmin = heat
         endif
      end do

100   continue
      nepnts =  ngeoms
      if (ihasep.eq.0) then
         igcvav = 0
         ieav = 0
      endif

      return
      end

      subroutine geoxyz(ihasg,iff,
     &      formax,forrms,dismax,disrms,epoints,isav,coo,ianz,istat)
      implicit double precision (a-h,o-z)
      logical goon
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      logical dozme
      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /athlp/ iatoms, mxnat
      character*137 line
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)
      dimension coo(3,*),ianz(*)

      istat = -1
      ngeoms = 0
      if (ihasg.gt.0) then
         call rewmf
         call srchmf(line,'[GEOMETRIES]',istat1)
      else
         rewind iun2
         igcvav = 1
         ifmxav = 0
         ifrmav = 0
         idmxav = 0
         idrmav = 0
         ieav = 1
         nepnts = 0
         dmaxt = 0.0d0
         fmaxt = 0.0d0
         drmst = 0.0d0
         frmst = 0.0d0
         enmax = -1000000.0d0
         enmin = 1000000.0d0
         ihasep = 0
         do i=1,mxpnt
            epoints(i) = 0.0d0
         end do
      endif

      do while (.true.)

         iheat = 0
         if (ihasg.le.1) then
             if (iftyp.eq.11.or.iftyp.eq.12) then
                if (iff.eq.7) then
                   call gtheat(igttnk,iheat,heat)
                else
                   call gettnk(igttnk,0,ipdbon,iff,iheat,heat)
                endif
                if (igttnk.eq.1) then
                   goon = .true.
                else
                   goon = .false.
                endif
             elseif (iftyp.eq.13) then
                call getmol(igetmo)
                if (igetmo.eq.1) then
                   goon = .true.
                else
                   goon = .false.
                endif
             else
                call getxyz(igetxy,heat,0)
                if (igetxy.eq.1) then
                   goon = .true.
                else
                   goon = .false.
                endif
                iheat = 1
             endif
         else
             call getzm(iatoms,0,0,istat1)
             if (istat1.eq.1) then
                 goon = .true.
             else
                 goon = .false.
             endif
         endif
         if (.not.goon) goto 80

         ngeoms = ngeoms + 1
         if (ihasg.eq.0.and.ngeoms.le.mxpnt) then
            isav(ngeoms) = 1
            if (iheat.eq.1) then
               ihasep = 1
               epoints(ngeoms) = heat
               if (heat.gt.enmax) enmax = heat
               if (heat.lt.enmin) enmin = heat
            endif
         endif
      end do

80    continue

      if (ngeoms.gt.mxpnt) then
          istat = ngeoms
          return
      endif

100   continue
      nepnts =  ngeoms
      if (enmax.eq.enmin) enmax = enmin + 0.00001d0
      if (ihasg.eq.0) then
         if (ihasep.eq.0) then
            igcvav = 0
            ieav = 0
         endif
      endif

      return
      end
