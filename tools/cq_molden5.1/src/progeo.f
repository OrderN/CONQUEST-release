      subroutine pargeo
      parameter (MAXPNT=2000)
      implicit double precision ( a-h,o-z )
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      common /geocnv2/ formax(MAXPNT),forrms(MAXPNT),dismax(MAXPNT),
     &                 disrms(MAXPNT),epoints(MAXPNT),isav(MAXPNT)

      mxpnt = MAXPNT
      call parptr(3,fmaxt,fdum,idum)
      call parptr(151,formax,fdum,idum)

      return
      end

      subroutine proged(ipoints,iff,istat,
     &                  formax,forrms,dismax,disrms,epoints,isav,
     &                  coo,ianz,icrtp)
C this is really progeo
      implicit double precision ( a-h,o-z )
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      logical dozme
      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb
      common /mopver/ mopopt
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*),coo(3,*),ianz(*)

      istat = -1

c can NOT set ngeoms = 0 here because some reading routines have
c already set ngeoms (for example:  mopac2007 aux files)
c                    (iftyp = 1, mopopt 4)

      if (iftyp.eq.1.and.iconv.eq.0) then
          call geomop(
     &            formax,forrms,dismax,disrms,epoints,isav,coo,ianz)
      elseif (iftyp.eq.2) then
          call geogam(formax,forrms,dismax,disrms,epoints,isav)
      elseif (iftyp.eq.3) then
          call geogus(formax,forrms,dismax,disrms,epoints,isav)
      elseif (iftyp.eq.4) then
          call geogau(formax,forrms,dismax,disrms,epoints,isav)
      elseif (iftyp.eq.5) then
          if (ihasg.gt.0) then
             call geoxyz(ihasg,iff,
     &       formax,forrms,dismax,disrms,epoints,isav,coo,ianz,needm)
             if (needm.ne.-1) then
                istat = needm
                return
             else
                call gmmcnv(formax,forrms,dismax,disrms,epoints,isav)
             endif
          else
             ngeoms = 0
          endif
      elseif (iftyp.eq.7) then
          call geocpmd(formax,forrms,dismax,disrms,epoints,isav,
     &                 needm)
          if (needm.ne.-1) istat = needm
      elseif (iftyp.eq.8) then
          call geoqcm(formax,forrms,dismax,disrms,epoints,isav)
      elseif (iftyp.eq.9) then
          call geoorc(formax,forrms,dismax,disrms,epoints,isav)
      elseif (iftyp.eq.10) then
          call geoxyz(0,iff,
     &      formax,forrms,dismax,disrms,epoints,isav,coo,ianz,needm)
          if (needm.ne.-1) then
             istat = needm
             return
          endif
      elseif (iftyp.eq.11) then
          if (iff.eq.7) then
             call geoxyz(0,iff,
     &        formax,forrms,dismax,disrms,epoints,isav,coo,ianz,needm)
          else
             call geoxyz(-1,iff,
     &        formax,forrms,dismax,disrms,epoints,isav,coo,ianz,needm)
          endif
          if (needm.ne.-1) then
             istat = needm
             return
          endif
      elseif (iftyp.eq.12) then
          call geoxyz(0,iff,
     &        formax,forrms,dismax,disrms,epoints,isav,coo,ianz,needm)
          if (needm.ne.-1) then
             istat = needm
             return
          endif
      elseif (iftyp.eq.14) then
          ngeoms = 1
      elseif (iftyp.eq.15) then
          call geonwc(formax,forrms,dismax,disrms,epoints,isav)
      elseif (icrtp.eq.5) then
          call geoxdt(formax,forrms,dismax,disrms,epoints,isav)
      elseif (icrtp.eq.3) then
          call geobio(formax,forrms,dismax,disrms,epoints,isav)
          call gmmcnv(formax,forrms,dismax,disrms,epoints,isav)
      elseif (icrtp.eq.2) then
          call geofdt(formax,forrms,dismax,disrms,epoints,isav)
          call gmmcnv(formax,forrms,dismax,disrms,epoints,isav)
      elseif (icrtp.eq.1) then
          call geochx(formax,forrms,dismax,disrms,epoints,isav)
          call gmmcnv(formax,forrms,dismax,disrms,epoints,isav)
      elseif (ipdbgro.eq.1) then
c done this already
          idum = 1
      elseif (iftyp.eq.1.and.mopopt.eq.4) then
          idum = 1
      else
          igcvav = 0
      endif

      if ((iftyp.ge.2.and.iftyp.le.4).or.(iftyp.ge.7.and.iftyp.le.9)
     &    .or.(iftyp.eq.15)
     &    .or.(iftyp.eq.1.and.(iconv.eq.0.or.mopopt.eq.4)) ) 
     &     call gmmcnv(formax,forrms,dismax,disrms,epoints,isav)

      ipoints = ngeoms

      return
      end
