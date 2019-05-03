      subroutine readvv(v,norbs,nocc,local)
      implicit double precision (a-h,o-z), integer ( i-n)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      character*137 line
      logical conti,local,reo,gamnew
      common /orbhlp/ mxorb,iuhf,ispd
      dimension v(*)
      dimension vt(6,12),irord(6)

      irord(1) = 3
      irord(2) = 1
      irord(3) = 2
      irord(4) = 5
      irord(5) = 6
      irord(6) = 4

      reo = .false.
      conti = .true.
      gamnew = .true.
      nlin = 8
      if (local) nlin = 3
c==== is format high or low ? =========
      max   = 7
      lprnt = norbs

      read(iun2,'(a)') line

      if (line(20:20).eq.'.') gamnew = .false.
      if (gamnew) then
         if (line(27:27).eq.' ') max = 10
      else
         if (line(25:25).eq.' ') max = 12
      endif
      if (max.ge.10.and..not.local) lprnt = min0(nocc+5,norbs)

      backspace iun2
      imax = 0
c======================================

100   imin = imax + 1
      imax = imax + max
      if ( imax .ge. lprnt ) then
         imax  = lprnt
         conti =.false.
      endif

      istrt = 0
      itel = 0
      nvec = imax - imin + 1
      do j=1,norbs
           read(iun2,'(a)',end=200) line
           if (max.ge.10) then
c=========== format low  =============
              if (gamnew) then
                 read(line,'(17x,10f9.4)',err=150)
     &            (v((i-1)*mxorb+j),i=imin,imax)
              else
                 read(line,'(15x,12f9.4)',err=150)
     &            (v((i-1)*mxorb+j),i=imin,imax)
              endif
           else
c=========== format high =============
              if (gamnew) then
                 read(line,'(17x,7f15.10)',err=150)
     &            (v((i-1)*mxorb+j),i=imin,imax)
              else
                 read(line,'(15x,7f15.10)',err=150)
     &            (v((i-1)*mxorb+j),i=imin,imax)
              endif
           endif
c
c Reorder the F functions, to fit the gaussian order
c
           if (istrt.eq.0.and.line(12:15).eq.' xxy') then
              reo =.true.
              istrt = j
              itel = 0
           endif
           if (istrt.ne.0) then
              itel = itel + 1
              do k=1,nvec
                 vt(itel,k) = v((imin+k-1-1)*mxorb+j)
              end do
           endif
           if (itel.eq.6) then
              do l=1,6
                 do k=1,nvec
                    v((imin+k-1-1)*mxorb+(istrt+l-1)) = vt(irord(l),k)
                 end do
              end do
              istrt = 0
              itel = 0
           endif
      end do

50    if (.not. conti) then
         if (reo) then
             print*,'========================================'
             print*,'     Changed order of F functions:'
             print*,' '
             print*,'     xxy,xxz,xyy,yyz,xzz,yzz ->'
             print*,' '
             print*,'     xyy,xxy,xxz,xzz,yzz,yyz'
             print*,'========================================'
         endif
         return
      endif

      do j=1,nlin
         read(iun2,'(a)',end=200) line
      end do
      goto 100
150   conti = .false.
      goto 50

200   call inferr('Error while reading vectors',1)
      end
