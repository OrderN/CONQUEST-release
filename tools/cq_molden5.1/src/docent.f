      subroutine docend(t,coo,ianz)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      common /athlp/ iatoms, mxnat
      dimension coo(3,*),ianz(*),t(*)

      call cntvec(t,coo,ianz,iatoms)

      return
      end

      subroutine doccd(coo,ianz)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      common /athlp/ iatoms, mxnat
      dimension vec(3)
      dimension coo(3,*),ianz(*)

      call cntvec(vec,coo,ianz,iatoms)

      do i=1,iatoms
         coo(1,i) = coo(1,i) - vec(1)
         coo(2,i) = coo(2,i) - vec(2)
         coo(3,i) = coo(3,i) - vec(3)
      end do

      return
      end

      subroutine cntvec(vec,coo,ianz,iatoms)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      dimension vec(3), coo(3,*), ianz(*)

      do i=1,3
         vec(i) = 0.0d0
      end do

      if (iatoms.le.0) return
      sum = 0.0d0
      do i=1,iatoms
         wat = dble(ianz(i))
         sum = sum + wat
         vec(1) = vec(1) + wat*coo(1,i)
         vec(2) = vec(2) + wat*coo(2,i)
         vec(3) = vec(3) + wat*coo(3,i)
      end do
 
      if (sum.gt.0.0d0) then
         vec(1) = vec(1) / sum
         vec(2) = vec(2) / sum
         vec(3) = vec(3) / sum
      endif

      return
      end

      subroutine cntvc2(vec,coo,iatoms)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      dimension vec(3), coo(3,*)

      do i=1,3
         vec(i) = 0.0d0
      end do

      if (iatoms.le.0) return
      do i=1,iatoms
         do j=1,3
            vec(j) = vec(j) + coo(j,i)
         end do
      end do
 
      do j=1,3
         vec(j) = vec(j) / dble(iatoms)
      end do

      return
      end

      subroutine dodcnt(td)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (numatm=2000)
      common /coord /xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      dimension td(3)

      do i=1,3
         td(i) = 0.0d0
      end do

      if (natoms.le.0) return

      sum = 0.0d0
      do i=1,natoms
         wat = dble(nat(i))
         sum = sum + wat
         do j=1,3
            td(j) = td(j) + wat*xyz(j,i)
         end do
      end do
 
      if (sum.gt.0.0d0) then
         do j=1,3
            td(j) = td(j) / sum
         end do
      endif

      return
      end
