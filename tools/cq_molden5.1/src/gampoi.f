      subroutine gampod(ipoint,istat,ioxyz,coo,ianz)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      parameter (maxfat=1000)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat) 
      common /athlp/ iatoms, mxnat

      character*2 gstr
      character*3 tstr
      character*5 hstr,hstrt
      character lstr*137
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      dimension coo(3,*),ianz(*)

cd     write(iun3,'(a)')'enter subroutine gampoi'

      if (ioxyz.eq.1) then
c
c        cartesian optim.
c
         if (ipoint.gt.99) then
            hstrt = hstr(ipoint)
            tstr = hstrt(2:4)
            call search(lstr,' point '//tstr,istat)
         else
            call search(lstr,' point  '//gstr(ipoint),istat)
         endif
         if ( istat .eq. 0 ) then
cd           write(iun3,'(a)')'leave subroutine gampoi'
            return
         endif
         call readel(lstr,2)
         call rdgeom(mxnat,coo,ianz,iatoms,istat)
         call haszm(.false.)
      else
c
c standard optim
c
         call searchd(lstr,'-  point '//gstr(ipoint),
     &                     '-  point  '//gstr(ipoint), istat)
         if ( istat .eq. 0 ) then
cd           write(iun3,'(a)')'leave subroutine gampoi'
            return
         endif
         call readel(lstr,7)
         call convzmat(coo,ianz,iatoms,0,1,0)
      endif

      call searchd(lstr,'gradient of the energy',
     &             '-  point',istat)
      if ( istat .eq. 0 ) return
      if ( index(lstr,'-  point').ne.0 .or. 
     &     index(lstr,'-  POINT').ne.0 ) then
           istat = -1
           backspace iun2
cd          write(iun3,'(a)')'leave subroutine gampoi'
           return
      endif

      call search(lstr,' atom',istmp)
      backspace iun2
      backspace iun2
      backspace iun2
      max = 0
140   min = max + 1
      max = max + 8
      if ( max .gt. iatoms ) max = iatoms
      call readel(lstr,5)
      do n=1,3
         read(iun2,'(6x,8(f15.7))',err=200)(fxyz(n,i),i=min,max)
      end do
      if ( max .lt. iatoms ) goto 140

cd     write(iun3,'(a)')'leave subroutine gampoi'
      return

200   istat =-1 
cd     write(iun3,'(a)')'leave subroutine gampoi'
      return
      end
