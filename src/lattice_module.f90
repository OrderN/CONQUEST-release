! $Id$
! -----------------------------------------------------------
! Module lattice_module
! -----------------------------------------------------------

!!****h* Conquest/lattice_module
!!  NAME
!!   lattice_module
!!  PURPOSE
!!   This module is preared for non-orthorhombic calculations.
!!   We now assume that cell_vector can be non-orthogonal to each other.
!!   (we may merge this module to some other modules in the future.)
!!
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   10/Nov/2025  
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
module lattice_module

  use datatypes
  use global_module, only: cell_vector, inv_cell_vector
  implicit none
  real(double) :: cell_length(1:3)
  real(double) :: volume

contains
 subroutine set_cell_parameters
  use mat33, only: inver3, deter3
  implicit none
  real(double):: bc1,bc2,bc3
  integer :: ii
 
  !cell length
  do ii=1,3
   cell_length(ii) = cell_vector(1,ii)**2 + cell_vector(2,ii)**2 + cell_vector(3,ii)**2
   cell_length(ii) = sqrt(cell_length(ii))
  enddo 

  !volume 
   volume = deter3(cell_vector)

  !inv_cell_vector
   call inver3(cell_vector, inv_cell_vector)

  return
 end subroutine set_cell_parameters

 subroutine get_pos_cart(pos_frac, pos_cart)
  use datatypes
  implicit none
  real(double), intent(in) :: pos_frac(3)
  real(double), intent(out):: pos_cart(3)

  !
   pos_cart(1) = cell_vector(1,1)*pos_frac(1) + cell_vector(1,2)*pos_frac(2)+cell_vector(1,3)*pos_frac(3)
   pos_cart(2) = cell_vector(2,1)*pos_frac(1) + cell_vector(2,2)*pos_frac(2)+cell_vector(2,3)*pos_frac(3)
   pos_cart(3) = cell_vector(3,1)*pos_frac(1) + cell_vector(3,2)*pos_frac(2)+cell_vector(3,3)*pos_frac(3)

  return
 end subroutine get_pos_cart

 subroutine get_pos_frac(pos_cart, pos_frac)
  use datatypes
  implicit none
  real(double), intent(in) :: pos_cart(3)
  real(double), intent(out):: pos_frac(3)

  !
   pos_frac(1) = inv_cell_vector(1,1)*pos_cart(1) + inv_cell_vector(1,2)*pos_cart(2)+inv_cell_vector(1,3)*pos_cart(3)
   pos_frac(2) = inv_cell_vector(2,1)*pos_cart(1) + inv_cell_vector(2,2)*pos_cart(2)+inv_cell_vector(2,3)*pos_cart(3)
   pos_frac(3) = inv_cell_vector(3,1)*pos_cart(1) + inv_cell_vector(3,2)*pos_cart(2)+inv_cell_vector(3,3)*pos_cart(3)

  return
 end subroutine get_pos_frac
 
end module lattice_module
