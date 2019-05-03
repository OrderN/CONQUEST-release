      subroutine tomold(rout,isel,inum)
      implicit double precision (a-h,o-z)
      real rout
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      dimension isel(4)

      call intcor(intc,rout,isel,inum)
      if (intc.eq.0) then
          write(iun3,*) 'Error calculating internal coordinate'
c          write(iun3,*) 'isel ',(isel(i),i=1,inum)
      endif

      return
      end

      subroutine setmon(isel,inum)
      implicit double precision (a-h,o-z)
      parameter (maxdm=20)
      real rout
      common /distmn/ rdm(maxdm),idmon(2,maxdm),ndm
      dimension isel(4)

      if (inum.eq.2) then
          if (ndm.lt.maxdm) then
             call intcor(intc,rout,isel,inum)
             if (intc.eq.1) then
                ndm = ndm + 1
                idmon(1,ndm) = isel(1)
                idmon(2,ndm) = isel(2)
                rdm(ndm) = dble(rout)*0.52917706d0
             endif
             call domcon(ndm,0)
          else
             call inferr('To many Distance Monitors',0)
          endif
      endif

      return
      end

      subroutine stjmon(coupl,isel,inum)
      implicit double precision (a-h,o-z)
      parameter (maxdm=20)
      common /distmn/ rdm(maxdm),idmon(2,maxdm),ndm
      dimension isel(4)

      if (inum.eq.2) then
          if (ndm.lt.maxdm) then
             ndm = ndm + 1
             idmon(1,ndm) = isel(1)
             idmon(2,ndm) = isel(2)
             rdm(ndm) = coupl
             call domcon(ndm,0)
          else
             call inferr('To many Monitors',0)
          endif
      endif

      return
      end

      subroutine clrmod(iconn)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (maxdm=20)
      common /distmn/ rdm(maxdm),idmon(2,maxdm),ndm
      logical found
      dimension icnn(mxcon),iconn(mxcon+1,*)

      do i=1,ndm

         k = idmon(1,i)
         n = iconn(1,k)
         icnt = 0
         found = .false.
         do j=n,1,-1
            if (iconn(1+j,k).ne.-idmon(2,i).or.found) then
                icnt = icnt + 1
                icnn(icnt) = iconn(1+j,k)
            else
                found = .true.
            endif
         end do
         do j=1,icnt
            iconn(1+(icnt-j+1),k) = icnn(j)
         end do
         iconn(1,k) = icnt

         k = idmon(2,i)
         n = iconn(1,k)
         icnt = 0
         found = .false.
         do j=n,1,-1
            if (iconn(1+j,k).ne.-idmon(1,i).or.found) then
                icnt = icnt + 1
                icnn(icnt) = iconn(1+j,k)
            else
                found = .true.
            endif
         end do
         do j=1,icnt
            iconn(1+(icnt-j+1),k) = icnn(j)
         end do
         iconn(1,k) = icnt

      end do
      
      ndm = 0

      return
      end

      subroutine intcod(intc,rout,isel,inum,coo)
      implicit double precision (a-h,o-z)
      parameter (small=1.0d-5)
      common /athlp/ iatoms, mxnat
      real rout
      dimension v1(3),v2(3),v3(3),c12(3),c23(3),isel(4),coo(3,*)

      intc = 1

      todeg = 45.0d0 / datan(1.0d0)

      do i=1,3
          v1(i) = coo(i,isel(1)) - coo(i,isel(2))
          if (inum.ge.3) v2(i) = coo(i,isel(3)) - coo(i,isel(2))
          if (inum.ge.4) v3(i) = coo(i,isel(4)) - coo(i,isel(3))
      end do

      if (inum.eq.2) then
          rout = vlen(v1)
      endif

      if (inum.eq.3) then
          call impsc(v1,v2,cosb)
          if (cosb.gt.1.0d0) cosb = 1.0d0
          if (cosb.lt.-1.0d0) cosb = -1.0d0
          rout = dacos(cosb)*todeg
      endif

      if (inum.eq.4) then

          rout = 0.0

          do i=1,4
             if (isel(i).lt.1.or.isel(i).gt.iatoms) goto 100
             do j=i+1,4
                if (isel(i).eq.isel(j)) goto 100
             end do
          end do

          call crprod(v1,v2,c12)

          do i=1,3
             v2(i) = -v2(i)
          end do

          call crprod(v2,v3,c23)
          if (vlen(c12).eq.0.0d0.or.vlen(c23).eq.0.0d0) goto 100

          call impsc(c12,c23,cosb)
          call crprod(c12,c23,v1)

          if (vlen(v1) .lt. small ) then
              dihsgn = 1.0d0
          else
              call impsc(v1,v2,dihsgn)
              if (dihsgn.ge.0.0d0) then
                 dihsgn = -1.0d0
              else
                 dihsgn = 1.0d0
              endif
          endif

          if (cosb.gt.1.0d0) cosb = 1.0d0
          if (cosb.lt.-1.0d0) cosb = -1.0d0
          rout = dacos(cosb)*todeg*dihsgn
      endif

      return

100   intc = 0
      return
      end
