module dynpro_math

  use kinds, only      : DP
  use constants,  only : pi, eps12
  use constants,  only : ps   => AU_PS
  use constants,  only : c    => C_SI
  use constants,  only : kb   => K_BOLTZMANN_SI
  use constants,  only : hp   => H_PLANCK_SI
  use constants,  only : bohr => BOHR_RADIUS_ANGS
  use constants,  only : ry   => AUTOEV

  !USE dynpro_psd, ONLY : 

  implicit none

contains

  subroutine check_power2(arg)

    implicit none

    integer, intent(IN)    :: arg
    integer :: pow2(31)
    integer :: i, iostat

    data pow2 / 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,    &
         2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, &
         524288, 1048576, 2097152, 4194304, 8388608, 16777216,  &
         33554432, 67108864, 134217728, 268435456, 536870912,   &
         1073741824 / !, 2147483648, 4294967296, 8589934592, 17179869184 


    i = 1 ; iostat = -1
    do while ( i < 32 )  
       if ( arg == pow2(i) ) iostat = 0
       i = i + 1
    end do

    if ( iostat == -1  ) then
       print*, '      ERROR: frame number must be an int^2 for Fourier/PSD'
       print*, ' Dynpro aborted'
       stop
    else
       return
    end if


    return
  end subroutine check_power2

  subroutine inverse(at, atinv)

    ! compute inverse of 3*3 matrix

    implicit none

    real(DP), intent(in)  :: at(3, 3)
    real(DP), intent(out) :: atinv(3, 3)

    real(DP) :: det

    atinv(1, 1) = at(2, 2) * at(3, 3) - at(2, 3) * at(3, 2)
    atinv(2, 1) = at(2, 3) * at(3, 1) - at(2, 1) * at(3, 3)
    atinv(3, 1) = at(2, 1) * at(3, 2) - at(2, 2) * at(3, 1)
    atinv(1, 2) = at(1, 3) * at(3, 2) - at(1, 2) * at(3, 3)
    atinv(2, 2) = at(1, 1) * at(3, 3) - at(1, 3) * at(3, 1)
    atinv(3, 2) = at(1, 2) * at(3, 1) - at(1, 1) * at(3, 2)
    atinv(1, 3) = at(1, 2) * at(2, 3) - at(1, 3) * at(2, 2)
    atinv(2, 3) = at(1, 3) * at(2, 1) - at(1, 1) * at(2, 3)
    atinv(3, 3) = at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1)

    det = at(1, 1) * atinv(1, 1) + at(1, 2) * atinv(2, 1) + &
         at(1, 3) * atinv(3, 1)
    atinv(:,:) = atinv(:,:) / det

    return
  end subroutine inverse


  subroutine vec_product(a1, a2, a3)

    implicit none

    real(DP), intent(in)         :: a1(3), a2(3)
    real(DP), intent(out)        :: a3(3)

    a3(1) = a1(2)*a2(3) - a1(3)*a2(2)
    a3(2) = a1(3)*a2(1) - a1(1)*a2(3)
    a3(3) = a1(1)*a2(2) - a1(2)*a2(1)

    return

  end subroutine VEC_PRODUCT

end module dynpro_math
