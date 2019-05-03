      subroutine geogam(formax,forrms,dismax,disrms,epoints,isav)
      implicit double precision (a-h,o-z)
      logical optxyz
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /geocnv/ fmaxt,frmst,dmaxt,drmst,fgmin,fgmax,dgmin,dgmax,
     &                enmax,enmin,ngeoms,nepnts,igcvav,ifmxav,ifrmav,
     &                idmxav,idrmav,ieav,ifrav,mxpnt
      character lstr*137
      dimension formax(*),forrms(*),dismax(*),disrms(*),epoints(*)
      dimension isav(*)

cd     write(iun3,'(a)')'enter subroutine geogam'
      rewind iun2
      igcvav = 1
      ifmxav = 1
      ifrmav = 1
      idmxav = 1
      idrmav = 1
      ieav = 1
      nepnts = 0
      ngeoms = 0
      dmaxt = 0.0d0
      fmaxt = 0.0d0
      drmst = 0.0d0
      frmst = 0.0d0

      optxyz = .false.
      call search(lstr,'optxyz',istat)
      if (istat.ne.0) optxyz = .true.
      if (optxyz) then
         ifrmav = 0
         idmxav = 0
         idrmav = 0
      endif

      do i=1,mxpnt
         if (optxyz) then
            call search(lstr,'nuclear coordinates',istat)
         else
            call search(lstr,'-  point',istat)
         endif
         if (istat.eq.0) goto 100

10       ibeg = index(lstr,'point') + 5
         lstr = lstr(ibeg:)
c
c...     Try using i4 format first if that fails there might
c...     be a space less so try i3 format, if that fails as well
c...     give up
c
         read(lstr,'(i4)',end=100,err=20) ntemp
         goto 30
20       read(lstr,'(i3)',end=100,err=100) ntemp
30       continue
         if (optxyz) ntemp = ntemp + 1

         if (optxyz) then
             call searchd(lstr,'total energy',
     &              'largest component of gradient',istat)
         else
             call searcht(lstr,'total energy','information on',
     &                     '-  point',istat)
         endif
         if (istat.eq.0) goto 100
         if (index(lstr,'total energy').eq.21) then
             if (linlen(lstr).eq.61) then
                read(lstr,'(43x,f18.10)',end=100,err=100) 
     &              epoints(ntemp)
             else if (linlen(lstr).eq.55) then
                read(lstr,'(37x,f18.10)',end=100,err=100) 
     &              epoints(ntemp)
             else if (linlen(lstr).eq.54) then
                read(lstr,'(38x,f16.10)',end=100,err=100) 
     &              epoints(ntemp)
             else if (linlen(lstr).eq.52) then
                read(lstr,'(37x,f16.10)',end=100,err=100) 
     &              epoints(ntemp)
             endif
             nepnts = ntemp
             if (optxyz) then
                 call search(lstr,'largest component of gradient',istat)
             else
                 call searchd(lstr,'information on',
     &                     '-  point',istat)
                 if (istat.eq.0) goto 100
             endif
         elseif (index(lstr,'total energy').eq.11) then
             if (linlen(lstr).eq.69) then
                read(lstr,'(33x,f18.10)',end=100,err=100) 
     &              epoints(ntemp)
             else if (linlen(lstr).eq.57) then
                read(lstr,'(39x,f18.10)',end=100,err=100) 
     &              epoints(ntemp)
             else if (linlen(lstr).eq.53) then
                read(lstr,'(37x,f16.10)',end=100,err=100) 
     &              epoints(ntemp)
             else if (linlen(lstr).eq.52) then
                read(lstr,'(37x,f15.8)',end=100,err=100) 
     &              epoints(ntemp)
             else
                read(lstr,'(28x,f18.10)',end=100,err=100)
     &              epoints(ntemp)
             endif
             nepnts = ntemp
             if (optxyz) then
                 call search(lstr,'largest component of gradient',istat)
             else
                 call searcht(lstr,'information on','information for',
     &                     '-  point',istat)
                 if (istat.eq.0) goto 100
             endif
         else
             epoints(ntemp) = 0.0d0
             nepnts = ntemp
         endif

         isav(ntemp) = 1
         if ( index(lstr,'-  point').ne.0 .or.
     &        index(lstr,'-  POINT').ne.0 ) then
            isav(ntemp) = 0
            ngeoms = ntemp
            goto 10
         endif

         if (optxyz) then
            ngeoms = ntemp
            read(lstr,'(61x,f15.7)',end=100,err=100)
     &            formax(ngeoms)
            forrms(ngeoms) = 0.0d0
            dismax(ngeoms) = 0.0d0
            disrms(ngeoms) = 0.0d0
         else
            call readel(lstr,2)

            read(iun2,'(21x,f14.8,15x,f14.8)',end=100,err=100) 
     &           tmp3,dmaxt
            read(iun2,'(21x,f14.8,15x,f14.8)',end=100,err=100) 
     &           tmp4,drmst
            read(iun2,'(21x,f14.8,15x,f14.8)',end=100,err=100) 
     &           tmp1,fmaxt
            read(iun2,'(21x,f14.8,15x,f14.8)',end=100,err=100) 
     &           tmp2,frmst
            ngeoms = ntemp
            formax(ngeoms) = tmp1
            forrms(ngeoms) = tmp2
            dismax(ngeoms) = tmp3
            disrms(ngeoms) = tmp4
         endif

      end do
100   continue
      if (ngeoms.eq.0) then
          igcvav = 0
          ifmxav = 0
          ifrmav = 0
          idmxav = 0
          idrmav = 0
          ieav = 0
      endif

cd     write(iun3,*)'*** ngeom=',ngeom
cd     write(iun3,'(a)')'leave subroutine geogam'
      return
      end
