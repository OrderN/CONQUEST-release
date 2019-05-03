      subroutine parori
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common/orihlp/ori(numatm),iuser(numatm),oalpha(numatm),
     $              obeta(numatm),norien,ibal

      norien = 0
      do k=1,natoms
          if (nat(k).eq.8.or.nat(k).eq.9.or.nat(k).eq.16
     &       .or.nat(k).eq.17) then
             norien        = norien + 1
             ori(norien)   = k
             iuser(norien) = 0
          endif
      end do

      return
      end
