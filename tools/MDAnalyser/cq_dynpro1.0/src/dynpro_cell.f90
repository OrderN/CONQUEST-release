module dynpro_cell
  !
  use kinds, only      : DP
  use constants,   only : pi, eps12
  use constants,   only : ps   => AU_PS
  use constants,   only : c    => C_SI
  use constants,   only : kb   => K_BOLTZMANN_SI
  use constants,   only : hp   => H_PLANCK_SI
  use constants,   only : bohr => BOHR_RADIUS_ANGS
  use constants,   only : ry   => AUTOEV
  !
  use dynpro_math, only : inverse
  !
  implicit none
  !
contains
  !
  subroutine xyz_to_abc(n, m, xyz, abc, cel)
    
    implicit none
    
    integer,  intent(in)         :: n, m
    real(DP), intent(in)         :: xyz(3,n*m)
    real(DP), intent(in)         :: cel(3,3*m)
    real(DP), intent(inout)      :: abc(3,n*m)
    
    real(DP)                     :: at(3,3), at_inv(3,3), sigma(3,n)
    integer                      :: i, j, k, l, q
    !
    k = 0
    l = 0
    do i = 1, m
       at(:,:) = cel(1:3,(i+k):(i*3))
       call inverse(at, at_inv)
       abc(:,(i+l):(i*n)) = matmul(at_inv(:,:), xyz(:,(i+l):(i*n)))
       k = k + 2
       l = l + (n - 1)
    end do
    !
  end subroutine XYZ_TO_ABC
  !
  !
  !
  subroutine center_cell(n, m, ct_atom, cel, abc, ct_abc, ct_xyz)
    
    implicit none
    
    integer,  intent(in)  :: n, m, ct_atom
    real(DP), intent(in)  :: abc(3,n*m), cel(3,3*m)
    real(DP), intent(out) :: ct_abc(3,n*m) , ct_xyz(3,n*m)
    
    real(DP)      :: sigma(3, n), sigma_catom(3, n)
    integer       :: i, j, k, l, q
    !
    k = 0
    l = 0
    !
    do q = 1, m
       !
       sigma(:,:)  = abc(:,(q+l):(q*n))
       !
       sigma_catom = 0.0d0
       !
       do i = 1, 3
          sigma_catom(i,:) =  sigma(i,:) - sigma(i,ct_atom)
          do j = 1, n
             if (sigma_catom(i,j) < -0.5d0) then
                sigma_catom(i,j) = sigma_catom(i,j) + 1.0d0
             else
                if (sigma_catom(i,j) > 0.5d0) then
                   sigma_catom(i,j) = sigma_catom(i,j) - 1.0d0
                endif
             end if
          end do
       end do
       !
       ct_abc(:,(q+l):(q*n)) = sigma_catom(:,1:n)
       ct_xyz(:,(q+l):(q*n)) = matmul(cel(:,(k+q):(q*3)), &
            ct_abc(:,(q+l):(q*n)))
       !
       k = k + 2
       l = l + (n - 1)
       !
    end do
    !
    return
  end subroutine center_cell




end module dynpro_cell
