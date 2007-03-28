! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_alloc
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_alloc *
!!
!!NAME
!! cq_ut_psp2pao_alloc
!!PURPOSE
!! Allocation of arrays
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 03/10/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_alloc

  implicit none

  interface cq_allocate
     module procedure cq_allocate_int
     module procedure cq_allocate_real
     module procedure cq_allocate_real_2
     module procedure cq_allocate_logical
     module procedure cq_allocate_orbital_info
     module procedure cq_allocate_basis_data
     module procedure cq_allocate_psp_in_table_fhi
  end interface

contains

! -----------------------------------------------------------------------------
! Subroutine allocate_writable_globals
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_alloc/allocate_writable_globals *
!!
!!NAME
!! allocate_writable_globals
!!USAGE
!!
!!PURPOSE
!! Allocates several global arrays that will be printed to the wavefunction 
!! file and will be used during the basis calculation
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine allocate_writable_globals(no_orbitals, points_mesh, psp_comps, partial_core)

  use cq_ut_psp2pao_global

  implicit none 

  ! Passed variables
  integer :: no_orbitals
  integer :: points_mesh
  integer :: psp_comps
  logical :: partial_core

  call cq_allocate(gl_orbitals, no_orbitals, 'allocate_writable_globals', 'gl_orbitals')
  call cq_allocate(gl_occ, 2 * no_orbitals, 'allocate_writable_globals', 'gl_occ')         ! Occupancies
  call cq_allocate(gl_r, points_mesh, 'allocate_writable_globals', 'gl_r')
  call cq_allocate(gl_rho, 2 * points_mesh, 'allocate_writable_globals', 'gl_rho')
  call cq_allocate(gl_ul, 2 * no_orbitals*points_mesh, 'allocate_writable_globals', 'gl_ul')
  call cq_allocate(gl_eigenvalues, 2 * no_orbitals, 'allocate_writable_globals', 'gl_eigenvalues')
  call cq_allocate(gl_v_nuclear, points_mesh, 'allocate_writable_globals', 'gl_v_nuclear')

  if(psp_comps > 0) then
     call cq_allocate(gl_v_nonlocal, psp_comps * points_mesh, 'allocate_writable_globals', 'gl_v_nonlocal')
  end if
  if(partial_core) then
     call cq_allocate(gl_rhopc, points_mesh, 'allocate_writable_globals', 'gl_rhopc')
  end if

  end subroutine allocate_writable_globals

! -----------------------------------------------------------------------------
! Subroutine allocate_psp_globals
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/allocate_psp_globals *
!!
!!NAME
!! allocate_psp_globals
!!USAGE
!!
!!PURPOSE
!! Allocates several global arrays that will be used to store the pseudopotential
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine allocate_psp_globals(points_psp, psp_comp, partial_core)

  use cq_ut_psp2pao_global
  use cq_ut_psp2pao_types

  implicit none 

  ! Passed variables
  integer :: points_psp
  integer :: psp_comp
  real(double) :: partial_core

  call cq_allocate(gl_psp_in%r, points_psp, 'allocate_psp_globals', 'gl_psp_in%r')
  call cq_allocate(gl_psp_in%v_nuclear, points_psp, 'allocate_psp_globals', 'gl_psp_in%v_nuclear')
  call cq_allocate(gl_psp_in%v_nuclear2, points_psp, 'allocate_psp_globals', 'gl_psp_in%v_nuclear2')
  call cq_allocate(gl_psp_in%l_nonlocal, psp_comp, 'allocate_psp_globals', 'gl_psp_in%l_nonlocal')
  call cq_allocate(gl_psp_in%sign_nonlocal, psp_comp, 'allocate_psp_globals', 'gl_psp_in%sign_nonlocal')
  call cq_allocate(gl_psp_in%v_nonlocal, psp_comp * points_psp, 'allocate_psp_globals', 'gl_psp_in%v_nonlocal')
  call cq_allocate(gl_psp_in%v_nonlocal2, psp_comp * points_psp, 'allocate_psp_globals', 'gl_psp_in%v_nonlocal2')
  call cq_allocate(gl_psp_in%v_nonlocal_cutoff, psp_comp, 'allocate_psp_globals', 'gl_psp_in%v_nonlocal_cutoff')

  if(partial_core > 0.0_double) then
     call cq_allocate(gl_psp_in%rhopc, points_psp, 'allocate_psp_globals', 'gl_psp_in%rhopc')
     call cq_allocate(gl_psp_in%rhopc2, points_psp, 'allocate_psp_globals', 'gl_psp_in%rhopc2')
  end if

  end subroutine allocate_psp_globals


! -----------------------------------------------------------------------------
! Subroutine allocate_wavefunction_set
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/allocate_wavefunction_set *
!!
!!NAME
!! allocate_wavefunction_set
!!USAGE
!!
!!PURPOSE
!! Allocates one wavefunction set
!! NOTE that the global array of set MUST have been allocated
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine allocate_wavefunction_set(no_set, no_orbitals, points_mesh, psp_comp, partial_core)

  use cq_ut_psp2pao_global

  implicit none 

  ! Passed variables
  integer :: no_set
  integer :: no_orbitals
  integer :: points_mesh
  integer :: psp_comp
  logical :: partial_core

  if(.not.associated(gl_wf_set)) then
    write(*,*) 'The array gl_wf_set has not been allocated when calling allocate_wavefunction_set'
    stop
  end if

  !*ast* TO DO: Fix the labels so they indicate the actual no_set

  call cq_allocate(gl_wf_set(no_set)%r, points_mesh, 'allocate_wavefunction_set', 'gl_wf_set(no_set)%r')
  call cq_allocate(gl_wf_set(no_set)%rho, 2*points_mesh, 'allocate_wavefunction_set', 'gl_wf_set(no_set)%rho')
  call cq_allocate(gl_wf_set(no_set)%l_nonlocal, psp_comp, &
                   'allocate_wavefunction_set', 'gl_wf_set(no_set)%l_nonlocal')
  call cq_allocate(gl_wf_set(no_set)%sign_nonlocal, psp_comp, &
                   'allocate_wavefunction_set', 'gl_wf_set(no_set)%sign_nonlocal')
  call cq_allocate(gl_wf_set(no_set)%v_local, points_mesh, 'allocate_wavefunction_set', 'gl_wf_set(no_set)%v_local')
  call cq_allocate(gl_wf_set(no_set)%v_nonlocal, points_mesh*psp_comp, &
                   'allocate_wavefunction_set', 'gl_wf_set(no_set)%v_nonlocal')
  call cq_allocate(gl_wf_set(no_set)%v_nonlocal_cutoff, psp_comp, &
                   'allocate_wavefunction_set', 'gl_wf_set(no_set)%v_nonlocal_cutoff')

  ! Second derivatives for splines
  call cq_allocate(gl_wf_set(no_set)%v_local2, points_mesh, 'allocate_wavefunction_set', 'gl_wf_set(no_set)%v_local2')
  call cq_allocate(gl_wf_set(no_set)%v_nonlocal2, points_mesh*psp_comp, &
                   'allocate_wavefunction_set', 'gl_wf_set(no_set)%v_nonlocal2')
  call cq_allocate(gl_wf_set(no_set)%orb_n, no_orbitals, 'allocate_wavefunction_set', 'gl_wf_set(no_set)%orb_n')
  call cq_allocate(gl_wf_set(no_set)%orb_l, no_orbitals, 'allocate_wavefunction_set', 'gl_wf_set(no_set)%orb_l')
  call cq_allocate(gl_wf_set(no_set)%orb_occ, no_orbitals, 'allocate_wavefunction_set', 'gl_wf_set(no_set)%orb_occ')
  call cq_allocate(gl_wf_set(no_set)%orb_keep, no_orbitals, 'allocate_wavefunction_set', 'gl_wf_set(no_set)%orb_keep')
  call cq_allocate(gl_wf_set(no_set)%orb_cutoff_radius, no_orbitals, &
                   'allocate_wavefunction_set', 'gl_wf_set(no_set)%orb_cutoff_radius')
  call cq_allocate(gl_wf_set(no_set)%orb_eigenvalues, no_orbitals, &
                   'allocate_wavefunction_set', 'gl_wf_set(no_set)%orb_eigenvalues')
  call cq_allocate(gl_wf_set(no_set)%orb_ul, no_orbitals, points_mesh, &
                   'allocate_wavefunction_set', 'gl_wf_set(no_set)%orb_ul')
  call cq_allocate(gl_wf_set(no_set)%orb_ul2, no_orbitals, points_mesh, &
                   'allocate_wavefunction_set', 'gl_wf_set(no_set)%orb_ul2')
  if(gl_wf_set(no_set)%partial_core) then
     call cq_allocate(gl_wf_set(no_set)%rhopc, points_mesh, 'allocate_wavefunction_set', 'gl_wf_set(no_set)%rhopc')
     call cq_allocate(gl_wf_set(no_set)%rhopc2, points_mesh, 'allocate_wavefunction_set', 'gl_wf_set(no_set)%rhopc2')
  end if

  end subroutine allocate_wavefunction_set


! -----------------------------------------------------------------------------
! Subroutine allocate_basis
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/allocate_basis *
!!
!!NAME
!! allocate_basis
!!USAGE
!!
!!PURPOSE
!! Allocates the basis
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine allocate_basis(no_orbitals, points_mesh, psp_comp, partial_core)

  use cq_ut_psp2pao_global

  implicit none 

  ! Passed variables
  integer :: no_orbitals
  integer :: points_mesh
  integer :: psp_comp
  logical :: partial_core

  call cq_allocate(gl_basis%r, points_mesh, 'allocate_basis', 'gl_basis%r')
  call cq_allocate(gl_basis%rho, 2 * points_mesh, 'allocate_basis', 'gl_basis%rho')
  call cq_allocate(gl_basis%l_nonlocal, psp_comp, 'allocate_basis', 'gl_basis%l_nonlocal')
  call cq_allocate(gl_basis%sign_nonlocal, psp_comp, 'allocate_basis', 'gl_basis%sign_nonlocal')
  call cq_allocate(gl_basis%v_local, points_mesh, 'allocate_basis', 'gl_basis%v_local')
  call cq_allocate(gl_basis%v_nonlocal, points_mesh * psp_comp, 'allocate_basis', 'gl_basis%v_nonlocal')
  call cq_allocate(gl_basis%v_nonlocal_cutoff, psp_comp, 'allocate_basis', 'gl_basis%v_nonlocal_cutoff')
  call cq_allocate(gl_basis%orb_n, no_orbitals, 'allocate_basis', 'gl_basis%orb_n')
  call cq_allocate(gl_basis%orb_l, no_orbitals, 'allocate_basis', 'gl_basis%orb_l')
  call cq_allocate(gl_basis%orb_cutoff_radius, no_orbitals, 'allocate_basis', 'gl_basis%orb_cutoff_radius')
  call cq_allocate(gl_basis%orb_keep, no_orbitals, 'allocate_basis', 'gl_basis%orb_keep')
  call cq_allocate(gl_basis%orb_occ, no_orbitals, 'allocate_basis', 'gl_basis%orb_occ')
  call cq_allocate(gl_basis%orb_eigenvalues, no_orbitals, 'allocate_basis', 'gl_basis%orb_eigenvalues')
  call cq_allocate(gl_basis%orb_ul, no_orbitals, points_mesh, 'allocate_basis', 'gl_basis%orb_ul')

  if(partial_core) then
     call cq_allocate(gl_basis%rhopc, points_mesh, 'allocate_basis', 'gl_basis%rhopc')
  end if

  end subroutine allocate_basis


! -----------------------------------------------------------------------------
! Subroutine deallocate_writable_globals
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_alloc/deallocate_writable_globals *
!!
!!NAME
!! deallocate_writable_globals
!!USAGE
!!
!!PURPOSE
!! Deallocates the globals that were allocated in allocate_writable_globals
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine deallocate_writable_globals(psp_comps, partial_core)

  use cq_ut_psp2pao_global

  implicit none 

  ! Passed variables
  integer :: psp_comps
  logical :: partial_core

  if(associated(gl_orbitals)) deallocate(gl_orbitals)
  if(associated(gl_occ))      deallocate(gl_occ)
  if(associated(gl_r))        deallocate(gl_r)
  if(associated(gl_rho))      deallocate(gl_rho)

  if(associated(gl_ul))          deallocate(gl_ul)
  if(associated(gl_eigenvalues)) deallocate(gl_eigenvalues)

  if(associated(gl_v_nuclear))     deallocate(gl_v_nuclear)
  if(psp_comps > 0) then
     if(associated(gl_v_nonlocal)) deallocate(gl_v_nonlocal)
  end if
  if(partial_core) then
     if(associated(gl_rhopc))      deallocate(gl_rhopc)
  end if

  end subroutine deallocate_writable_globals


! -----------------------------------------------------------------------------
! Subroutine deallocate_psp_globals
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/deallocate_psp_globals *
!!
!!NAME
!! deallocate_psp_globals
!!USAGE
!!
!!PURPOSE
!! Deallocates the globals that were allocated in allocate_psp_globals
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine deallocate_psp_globals(partial_core)

  use cq_ut_psp2pao_global

  implicit none 

  ! Passed variables
  logical :: partial_core

  if(associated(gl_psp_in%r))                    deallocate(gl_psp_in%r)
  if(associated(gl_psp_in%v_nuclear))            deallocate(gl_psp_in%v_nuclear)
  if(associated(gl_psp_in%v_nuclear2))           deallocate(gl_psp_in%v_nuclear2)
  if(associated(gl_psp_in%l_nonlocal))           deallocate(gl_psp_in%l_nonlocal)
  if(associated(gl_psp_in%sign_nonlocal))        deallocate(gl_psp_in%sign_nonlocal)
  if(associated(gl_psp_in%v_nonlocal))           deallocate(gl_psp_in%v_nonlocal)
  if(associated(gl_psp_in%v_nonlocal2))          deallocate(gl_psp_in%v_nonlocal2)
  if(associated(gl_psp_in%v_nonlocal_cutoff))    deallocate(gl_psp_in%v_nonlocal_cutoff)

  if(partial_core) then
     if(associated(gl_psp_in%rhopc))      deallocate(gl_psp_in%rhopc)
     if(associated(gl_psp_in%rhopc2))     deallocate(gl_psp_in%rhopc2)
  end if

  end subroutine deallocate_psp_globals


! -----------------------------------------------------------------------------
! Subroutine deallocate_wavefunction_sets
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/deallocate_wavefunction_sets *
!!
!!NAME
!! deallocate_wavefunction_sets
!!USAGE
!!
!!PURPOSE
!! Deallocates all wavefunction set and the global array
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine deallocate_wavefunction_sets(no_sets)

  use cq_ut_psp2pao_global

  implicit none 

  ! Passed variables
  integer :: no_sets

  ! Local variables
  integer :: i

  do i=1,no_sets
     if(associated(gl_wf_set(i)%r))                  deallocate(gl_wf_set(i)%r)
     if(associated(gl_wf_set(i)%rho))                deallocate(gl_wf_set(i)%rho)
     if(associated(gl_wf_set(i)%l_nonlocal))         deallocate(gl_wf_set(i)%l_nonlocal)
     if(associated(gl_wf_set(i)%sign_nonlocal))      deallocate(gl_wf_set(i)%sign_nonlocal)
     if(associated(gl_wf_set(i)%v_local))            deallocate(gl_wf_set(i)%v_local)
     if(associated(gl_wf_set(i)%v_nonlocal))         deallocate(gl_wf_set(i)%v_nonlocal)
     if(associated(gl_wf_set(i)%v_nonlocal_cutoff))  deallocate(gl_wf_set(i)%v_nonlocal_cutoff)
     if(associated(gl_wf_set(i)%v_local2))           deallocate(gl_wf_set(i)%v_local2)
     if(associated(gl_wf_set(i)%v_nonlocal2))        deallocate(gl_wf_set(i)%v_nonlocal2)
     if(associated(gl_wf_set(i)%orb_n))              deallocate(gl_wf_set(i)%orb_n)
     if(associated(gl_wf_set(i)%orb_l))              deallocate(gl_wf_set(i)%orb_l)
     if(associated(gl_wf_set(i)%orb_occ))            deallocate(gl_wf_set(i)%orb_occ)
     if(associated(gl_wf_set(i)%orb_keep))           deallocate(gl_wf_set(i)%orb_keep)
     if(associated(gl_wf_set(i)%orb_cutoff_radius))  deallocate(gl_wf_set(i)%orb_cutoff_radius)
     if(associated(gl_wf_set(i)%orb_eigenvalues))    deallocate(gl_wf_set(i)%orb_eigenvalues)
     if(associated(gl_wf_set(i)%orb_ul))             deallocate(gl_wf_set(i)%orb_ul)
     if(associated(gl_wf_set(i)%orb_ul2))            deallocate(gl_wf_set(i)%orb_ul2)

     if(gl_wf_set(i)%partial_core) then
        if(associated(gl_wf_set(i)%rhopc))           deallocate(gl_wf_set(i)%orb_ul2)
        if(associated(gl_wf_set(i)%rhopc2))          deallocate(gl_wf_set(i)%orb_ul2)
     end if
   end do

   if(associated(gl_wf_set))  deallocate(gl_wf_set)

  end subroutine deallocate_wavefunction_sets


! -----------------------------------------------------------------------------
! Subroutine deallocate_basis
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/deallocate_basis *
!!
!!NAME
!! allocate_basis
!!USAGE
!!
!!PURPOSE
!! Deallocates the basis
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine deallocate_basis(partial_core)

  use cq_ut_psp2pao_global

  implicit none 

  ! Passed variables

  logical :: partial_core

  if(associated(gl_basis%r))                 deallocate(gl_basis%r)
  if(associated(gl_basis%rho))               deallocate(gl_basis%rho)
  if(associated(gl_basis%l_nonlocal))        deallocate(gl_basis%l_nonlocal)
  if(associated(gl_basis%sign_nonlocal))     deallocate(gl_basis%sign_nonlocal)
  if(associated(gl_basis%v_local))           deallocate(gl_basis%v_local)
  if(associated(gl_basis%v_nonlocal))        deallocate(gl_basis%v_nonlocal)
  if(associated(gl_basis%v_nonlocal_cutoff)) deallocate(gl_basis%v_nonlocal_cutoff)
  if(associated(gl_basis%orb_n))             deallocate(gl_basis%orb_n)
  if(associated(gl_basis%orb_l))             deallocate(gl_basis%orb_l)
  if(associated(gl_basis%orb_cutoff_radius)) deallocate(gl_basis%orb_cutoff_radius)
  if(associated(gl_basis%orb_keep))          deallocate(gl_basis%orb_keep)
  if(associated(gl_basis%orb_occ))           deallocate(gl_basis%orb_occ)
  if(associated(gl_basis%orb_eigenvalues))   deallocate(gl_basis%orb_eigenvalues)
  if(associated(gl_basis%orb_ul))            deallocate(gl_basis%orb_ul)
  if(partial_core) then
     if(associated(gl_basis%rhopc))          deallocate(gl_basis%rhopc)
  end if

  end subroutine deallocate_basis


! -----------------------------------------------------------------------------
! Subroutine deallocate_nonwritable_globals
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_alloc/deallocate_nonwritable_globals *
!!
!!NAME
!! deallocate_nonwritable_globals
!!USAGE
!!
!!PURPOSE
!! Deallocates the globals that are not used for writing to file
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine deallocate_nonwritable_globals

  use cq_ut_psp2pao_global

  implicit none 

  if(associated(gl_norm))               deallocate(gl_norm)
  if(associated(gl_v_nonlocal_int_in))  deallocate(gl_v_nonlocal_int_in)
  if(associated(gl_v_nonlocal_int_out)) deallocate(gl_v_nonlocal_int_out)
  if(associated(gl_v_ext))              deallocate(gl_v_ext)
  if(associated(gl_v_hartree))          deallocate(gl_v_hartree)
  if(associated(gl_v_xc))               deallocate(gl_v_xc)
  if(associated(gl_v))                  deallocate(gl_v)
  if(associated(gl_tmp1))               deallocate(gl_tmp1)
  if(associated(gl_tmp2))               deallocate(gl_tmp2)
  if(associated(gl_tmp3))               deallocate(gl_tmp3)
  if(associated(gl_dtmp1))              deallocate(gl_dtmp1)
  if(associated(gl_dtmp2))              deallocate(gl_dtmp2)
  if(associated(gl_dtmp3))              deallocate(gl_dtmp3)

  end subroutine deallocate_nonwritable_globals


! -----------------------------------------------------------------------------
! Subroutine cq_allocate_int
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/cq_allocate_int *
!!
!!NAME
!! cq_allocate_int
!!USAGE
!!
!!PURPOSE
!! Allocate an array of integers
!!INPUTS
!!
!!USES
!!
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine cq_allocate_int(array, size, routine, message)

  implicit none

  ! Passed variables
  integer, dimension(:), pointer :: array
  integer :: size
  character(len=*) :: routine
  character(len=*) :: message

  ! Local variables
  integer :: stat

  nullify(array)
  allocate(array(size), STAT=stat)
  if(stat /= 0) then
     write(*,*) routine,': Error allocating ',message
     stop
  end if

  end subroutine cq_allocate_int

! -----------------------------------------------------------------------------
! Subroutine cq_allocate_real
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/cq_allocate_real *
!!
!!NAME
!! cq_allocate_real
!!USAGE
!!
!!PURPOSE
!! Allocate a real(double) array
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine cq_allocate_real(array, size, routine, message)

  use cq_ut_psp2pao_types

  implicit none

  ! Passed variables
  real(double), dimension(:), pointer :: array
  integer :: size
  character(len=*) :: routine
  character(len=*) :: message

  ! Local variables
  integer :: stat

  nullify(array)
  allocate(array(size), STAT=stat)
  if(stat /= 0) then
     write(*,*) routine,': Error allocating ',message
     stop
  end if

  end subroutine cq_allocate_real

! -----------------------------------------------------------------------------
! Subroutine cq_allocate_real_2
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/cq_allocate_real_2 *
!!
!!NAME
!! cq_allocate_real_2
!!USAGE
!!
!!PURPOSE
!! Allocate a real(double) 2D-array
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine cq_allocate_real_2(array, size1, size2, routine, message)

  use cq_ut_psp2pao_types

  implicit none

  ! Passed variables
  real(double), dimension(:,:), pointer :: array
  integer :: size1
  integer :: size2
  character(len=*) :: routine
  character(len=*) :: message

  ! Local variables
  integer :: stat

  nullify(array)
  allocate(array(size1,size2), STAT=stat)
  if(stat /= 0) then
     write(*,*) routine,': Error allocating ',message
     stop
  end if

  end subroutine cq_allocate_real_2

! -----------------------------------------------------------------------------
! Subroutine cq_allocate_logical
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/cq_allocate_logical *
!!
!!NAME
!! cq_allocate_logical
!!USAGE
!!
!!PURPOSE
!! Allocate an array of logicals
!!INPUTS
!!
!!USES
!!
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine cq_allocate_logical(array, size, routine, message)

  implicit none

  ! Passed variables
  logical, dimension(:), pointer :: array
  integer :: size
  character(len=*) :: routine
  character(len=*) :: message

  ! Local variables
  integer :: stat

  nullify(array)
  allocate(array(size), STAT=stat)
  if(stat /= 0) then
     write(*,*) routine,': Error allocating ',message
     stop
  end if

  end subroutine cq_allocate_logical

! -----------------------------------------------------------------------------
! Subroutine cq_allocate_orbital_info
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/cq_allocate_orbital_info *
!!
!!NAME
!! cq_allocate_orbital_info
!!USAGE
!!
!!PURPOSE
!! Allocate an orbital_info array
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine cq_allocate_orbital_info(array, size, routine, message)

  use cq_ut_psp2pao_types

  implicit none

  ! Passed variables
  type(orbital_info), dimension(:), pointer :: array
  integer :: size
  character(len=*) :: routine
  character(len=*) :: message

  ! Local variables
  integer :: stat

  nullify(array)
  allocate(array(size), STAT=stat)
  if(stat /= 0) then
     write(*,*) routine,': Error allocating ',message
     stop
  end if

  end subroutine cq_allocate_orbital_info

! -----------------------------------------------------------------------------
! Subroutine cq_allocate_basis_data
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/cq_allocate_basis_data *
!!
!!NAME
!! cq_allocate_basis_data
!!USAGE
!!
!!PURPOSE
!! Allocate a basis_data array
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine cq_allocate_basis_data(array, size, routine, message)

  use cq_ut_psp2pao_types

  implicit none

  ! Passed variables
  type(basis_data), dimension(:), pointer :: array
  integer :: size
  character(len=*) :: routine
  character(len=*) :: message

  ! Local variables
  integer :: stat

  nullify(array)
  allocate(array(size), STAT=stat)
  if(stat /= 0) then
     write(*,*) routine,': Error allocating ',message
     stop
  end if

  end subroutine cq_allocate_basis_data

! -----------------------------------------------------------------------------
! Subroutine cq_allocate_psp_in_table_fhi
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_alloc/cq_allocate_psp_in_table_fhi *
!!
!!NAME
!! cq_allocate_psp_in_table_fhi
!!USAGE
!!
!!PURPOSE
!! Allocate a psp_in_table_fhi array
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 09/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine cq_allocate_psp_in_table_fhi(array, size, routine, message)

  use cq_ut_psp2pao_types

  implicit none

  ! Passed variables
  type(psp_in_table_fhi), dimension(:), pointer :: array
  integer :: size
  character(len=*) :: routine
  character(len=*) :: message

  ! Local variables
  integer :: stat

  nullify(array)
  allocate(array(size), STAT=stat)
  if(stat /= 0) then
     write(*,*) routine,': Error allocating ',message
     stop
  end if

  end subroutine cq_allocate_psp_in_table_fhi

end module cq_ut_psp2pao_alloc

