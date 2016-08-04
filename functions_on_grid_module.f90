! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module 
! ------------------------------------------------------------------------------
! Code area 8: indexing
! ------------------------------------------------------------------------------

!!****h* Conquest/functions_on_grid *
!!  NAME
!!   functions_on_grid
!!  PURPOSE
!!   Holds various variables relating to functions (e.g. support
!!   functions, NL pseudopotential projectors, PAOs) on the grid and
!!   all associated routines and variables (like sizes etc)
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2006/03/04
!!  MODIFICATION HISTORY
!!   Added overflow checks and headers
!!   2008/05/16 ast
!!    Added timer
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2016/08/01 17:30 nakata
!!    Introduced atomf instead of sf and paof
!!  SOURCE
!!
module functions_on_grid

  use datatypes
  use global_module,          only: atomf, nlpf
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_allocation

  implicit none
  save

!!***

  !!****s* functions_on_grid/fn_on_grid *
  !!  NAME
  !!   fn_on_grid
  !!  PURPOSE
  !!   Contains data relating to functions on grid
  !!  AUTHOR
  !!   D.R.Bowler
  !!  MODIFICATION HISTORY
  !!   2011/09/19 L.Tong
  !!     Added new function on grid H_dn_on_supportfns As the name
  !!     suggests it should be used to store the spin down component
  !!     of H acting on supportfns.  Because pseudofns is optional,
  !!     H_dn_on_supportfns will not be set as a parameter, but
  !!     instead will be set as 3 or 4 depending on whether pseudofns
  !!     is present or not
  !!   2012/03/13 L.Tong
  !!   - Changed spin implementation, now H_on_supportfns is a
  !!     constant array of dimension 2, corresponding to two spin
  !!     components.
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomf
  !!  SOURCE
  !!
  type fn_on_grid

     real(double), dimension(:), pointer :: griddata
     integer :: size
     integer :: type

  end type fn_on_grid
  !!***

  integer, dimension(3) :: gridsize
  type(fn_on_grid), dimension(:), allocatable :: gridfunctions

  integer, parameter :: supportfns = 1
  integer, parameter, dimension(2) :: H_on_atomf = (/2, 3/)
  integer, parameter :: pseudofns = 4
  integer :: current_fn_on_grid
  integer, parameter :: mx_fns_on_grid = 20

  !real(double), target, dimension(SUPPORT_SIZE*NSF) :: support, workspace_support, workspace2_support
  
  ! RCS tag for object file identification
  character(len=80), private :: &
       RCSid = "$Id$"

contains

!!****f* functions_on_grid/associate_fn_on_grid *
!!
!!  NAME 
!!   associate_fn_on_grid
!!  USAGE
!!   
!!  PURPOSE
!!   Allocates space for functions on grid
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2006 sometime
!!  MODIFICATION HISTORY
!!   2007/02/01 14:02 dave
!!    Added if checks for zero size to stop underflow errors on Hitachi
!!   2008/05/16 ast
!!    Added timer
!!   2011/09/19 L.Tong
!!    Added routines for H_dn_on_supportfns, this is used only if
!!    doing spin polarised calculations. Note that H_dn_on_supportfns
!!    = 3 if pseudofns not present, otherwise it is set to be 4
!!   2012/03/13 L.Tong
!!    Rewrote spin implementation. Now H_on_supportfns is defined as
!!    an array (/2, 3/). H_on_supportfns(spin) will then denote the
!!    array for the given spin channel.
!!   2016/07/13 18:30 nakata
!!    Renamed H_on_supportfns -> H_on_atomf
!!   2016/08/01 17:30 nakata
!!    Introduced atomf instead of sf
!!  SOURCE
!!  
  subroutine associate_fn_on_grid

    use GenComms, only: cq_abort
    use global_module, only: flag_basis_set, blips, area_index, &
                             nspin, flag_analytic_blip_int
    use memory_module, only: reg_alloc_mem, type_dbl
    use numbers, only: zero

    implicit none

    integer :: stat, i, spin

    call start_timer(tmr_std_allocation)
    if (.not. allocated(gridfunctions)) then
       allocate(gridfunctions(mx_fns_on_grid), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating gridfunctions ", mx_fns_on_grid)
    end if

    ! Support functions
    gridfunctions(supportfns)%size = gridsize(atomf)
    allocate(gridfunctions(supportfns)%griddata(gridsize(atomf)))
    if (gridsize(atomf) > 0) then
       do i = 1, gridsize(atomf)
          gridfunctions(supportfns)%griddata(i) = zero
       end do
    end if
    call reg_alloc_mem(area_index, gridsize(atomf), type_dbl)
    gridfunctions(supportfns)%type = atomf

    ! H acting on support functions
    do spin = 1, nspin
       gridfunctions(H_on_atomf(spin))%size = gridsize(atomf)
       allocate(gridfunctions(H_on_atomf(spin))%griddata(gridsize(atomf)))
       if(gridsize(atomf)>0) then
          do i=1,gridsize(atomf)
             gridfunctions(H_on_atomf(spin))%griddata(i) = zero
          end do
       end if
       call reg_alloc_mem(area_index, gridsize(atomf), type_dbl)
       gridfunctions(H_on_atomf(spin))%type = atomf
    end do

    if (flag_basis_set == blips .and. (.not. flag_analytic_blip_int)) then
       ! Non-local projector functions
       gridfunctions(pseudofns)%size = gridsize(nlpf)
       allocate(gridfunctions(pseudofns)%griddata(gridsize(nlpf)))
       if(gridsize(nlpf)>0) then
          do i=1,gridsize(nlpf)
             gridfunctions(pseudofns)%griddata(i) = zero
          end do
       end if
       call reg_alloc_mem(area_index,gridsize(nlpf),type_dbl)
       gridfunctions(pseudofns)%type = nlpf
       ! Matrix acting on support functions
       current_fn_on_grid = 4
    else
       current_fn_on_grid = 3
    end if

    call stop_timer(tmr_std_allocation)
  end subroutine associate_fn_on_grid
!!***

  !!****f* functions_on_grid/dissociate_fn_on_grid *
  !!
  !!  NAME 
  !!   dissociate_fn_on_grid
  !!  USAGE
  !!   
  !!  PURPOSE
  !!   Frees space for functions on grid
  !!  INPUTS
  !!   
  !!   
  !!  USES
  !!   
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2006 sometime
  !!  MODIFICATION HISTORY
  !!   2007/02/01 14:02 dave
  !!    Added if checks for zero size to stop underflow errors on Hitachi
  !!   2008/05/16 ast
  !!    Added timer
  !!   2011/09/19 L.Tong
  !!    Added routines for H_dn_on_supportfns
  !!   2012/03/13 L.Tong
  !!    Changed implementation for spin polarisation. 
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomf
  !!  SOURCE
  !!  
  subroutine dissociate_fn_on_grid

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: flag_basis_set, blips, area_index, &
                             nspin, flag_analytic_blip_int
    use memory_module, ONLY: reg_dealloc_mem, type_dbl

    implicit none

    integer :: stat, spin

    call start_timer(tmr_std_allocation)
    ! Support functions
    deallocate(gridfunctions(supportfns)%griddata)
    call reg_dealloc_mem (area_index,                               &
                          size(gridfunctions(supportfns)%griddata), &
                          type_dbl)

    ! H acting on support functions 
    do spin = 1, nspin
       deallocate(gridfunctions(H_on_atomf(spin))%griddata)
       call reg_dealloc_mem(area_index,                                    &
                            size(gridfunctions(H_on_atomf(spin))%griddata), &
                            type_dbl)
    end do

    if (flag_basis_set == blips .and. &
        (.not. flag_analytic_blip_int)) then
       ! Non-local projector functions
       deallocate(gridfunctions(pseudofns)%griddata)
       call reg_dealloc_mem(area_index,                              &
                            size(gridfunctions(pseudofns)%griddata), &
                            type_dbl)
    end if
    
    current_fn_on_grid = 0
    call stop_timer(tmr_std_allocation)
  end subroutine dissociate_fn_on_grid
  !!***

!!****f* functions_on_grid/allocate_temp_fn_on_grid *
!!
!!  NAME 
!!   allocate_temp_fn_on_grid
!!  USAGE
!!   
!!  PURPOSE
!!   Allocates space for functions on grid (temporary)
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2006 sometime
!!  MODIFICATION HISTORY
!!   2007/02/01 14:02 dave
!!    Added if checks for zero size to stop underflow errors on Hitachi
!!   2008/05/16 ast
!!    Added timer
!!  SOURCE
!!  
  function allocate_temp_fn_on_grid(type)

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_index
    use memory_module, ONLY: reg_alloc_mem, type_dbl
    use numbers, ONLY: zero

    implicit none

    ! Passed variables
    integer :: type
    ! Result
    integer :: allocate_temp_fn_on_grid
    ! Local variables
    integer :: size, stat, i

    call start_timer(tmr_std_allocation)
    ! Increment variable
    current_fn_on_grid = current_fn_on_grid + 1
    if(current_fn_on_grid>mx_fns_on_grid) call cq_abort("allocate_temp_fn: increase fns on grid ",&
         current_fn_on_grid,mx_fns_on_grid)
    ! Set size
    size = gridsize(type)
    call reg_alloc_mem(area_index,size,type_dbl)
    ! Allocate space and assign variables
    stat = 0 
    allocate(gridfunctions(current_fn_on_grid)%griddata(size),STAT=stat)
    if(size>0) then
       do i=1,size
          gridfunctions(current_fn_on_grid)%griddata(i) = zero
       end do
    end if
    if(stat/=0) call cq_abort("Error allocating temporary grid function ",current_fn_on_grid,size)
    gridfunctions(current_fn_on_grid)%size = size
    gridfunctions(current_fn_on_grid)%type = type
    allocate_temp_fn_on_grid = current_fn_on_grid
    call stop_timer(tmr_std_allocation)
  end function allocate_temp_fn_on_grid
!!***

!!****f* functions_on_grid/free_temp_fn_on_grid *
!!
!!  NAME 
!!   free_temp_fn_on_grid
!!  USAGE
!!   
!!  PURPOSE
!!   Frees space for functions on grid (temporary)
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2006 sometime
!!  MODIFICATION HISTORY
!!   2007/02/01 14:02 dave
!!    Added if checks for zero size to stop underflow errors on Hitachi
!!   2008/05/16 ast
!!    Added timer
!!  SOURCE
!!  
  subroutine free_temp_fn_on_grid(num)

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_index
    use memory_module, ONLY: reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    integer :: num

    call start_timer(tmr_std_allocation)
    if(num/=current_fn_on_grid) call cq_abort("Out-of-order deallocation of temp grid fn",num)
    call reg_dealloc_mem(area_index,gridfunctions(current_fn_on_grid)%size,type_dbl)
    deallocate(gridfunctions(current_fn_on_grid)%griddata)
    gridfunctions(current_fn_on_grid)%size = 0
    gridfunctions(current_fn_on_grid)%type = 0
    current_fn_on_grid = current_fn_on_grid - 1
    call stop_timer(tmr_std_allocation)
  end subroutine free_temp_fn_on_grid
!!***

end module functions_on_grid
