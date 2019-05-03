      subroutine xyzcod(idocopy,idoconv,ioadd,
     &                  ianz,coo)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /athlp/  iatoms, mxnat
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /surf/   natorg,noscnd
      dimension coo(3,*),ianz(*)

      toang = 0.52917706d0

      if (idocopy.eq.1) iatoms = natoms

      noff = 1
      if (ioadd.eq.1) noff = natorg + 1

      do i=noff,iatoms
         if (idocopy.eq.1) ianz(i) = nat(i)
         do j=1,3
           if (idoconv.eq.1) then
              if (idocopy.eq.1) then
                 coo(j,i) = xyz(j,i) / toang
              else
                 coo(j,i) = coo(j,i) / toang
              endif
           else
              coo(j,i) = xyz(j,i)
           endif
         end do
      end do

      if (ioadd.eq.1) natorg = 0

      return
      end

      subroutine cooxyz(ianz,iatoms)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      dimension ianz(*)
      
      natoms = iatoms
      if (natoms.gt.numatm) natoms = numatm

      do i=1,natoms
         nat(i) = ianz(i)
      end do

      return
      end
