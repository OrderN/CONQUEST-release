! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: functions_on_grid_module.f90,v 1.1.2.2 2006/03/31 12:17:35 drb Exp $
! ------------------------------------------------------------------------------
! Module 
! ------------------------------------------------------------------------------
! Code area 8: indexing
! ------------------------------------------------------------------------------

!!****h* Conquest/functions_on_grid *
!!  NAME
!!   functions_on_grid
!!  PURPOSE
!!   Holds various variables relating to functions (e.g. support functions, NL pseudopotential projectors, PAOs) 
!!   on the grid and all associated routines and variables (like sizes etc)
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2006/03/04
!!  MODIFICATION HISTORY
!!   Added overflow checks and headers
!!  SOURCE
!!
module functions_on_grid

  use datatypes
  use global_module, ONLY: sf, nlpf, paof

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
  integer, parameter :: H_on_supportfns = 2
  integer, parameter :: pseudofns = 3
  integer :: current_fn_on_grid
  integer, parameter :: mx_fns_on_grid = 20

  !real(double), target, dimension(SUPPORT_SIZE*NSF) :: support, workspace_support, workspace2_support
  
  ! RCS tag for object file identification
  character(len=80), private :: RCSid = "$Id: functions_on_grid_module.f90,v 1.1.2.2 2006/03/31 12:17:35 drb Exp $"

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
!!  SOURCE
!!  
  subroutine associate_fn_on_grid

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: flag_basis_set, blips, area_index
    use memory_module, ONLY: reg_alloc_mem, type_dbl
    use numbers, ONLY: zero

    implicit none

    integer :: stat, i

    if(.not.allocated(gridfunctions)) then
       allocate(gridfunctions(mx_fns_on_grid),STAT=stat)
       if(stat/=0) call cq_abort("Error allocating gridfunctions ",mx_fns_on_grid)
    end if
    ! Support functions
    gridfunctions(supportfns)%size = gridsize(sf)
    allocate(gridfunctions(supportfns)%griddata(gridsize(sf)))
    if(gridsize(sf)>0) then
       do i=1,gridsize(sf)
          gridfunctions(supportfns)%griddata(i) = zero
       end do
    end if
    call reg_alloc_mem(area_index,gridsize(sf),type_dbl)
    gridfunctions(supportfns)%type = sf
    ! H acting on support functions 
    gridfunctions(H_on_supportfns)%size = gridsize(sf)
    allocate(gridfunctions(H_on_supportfns)%griddata(gridsize(sf)))
    if(gridsize(sf)>0) then
       do i=1,gridsize(sf)
          gridfunctions(H_on_supportfns)%griddata(i) = zero
       end do
    end if
    call reg_alloc_mem(area_index,gridsize(sf),type_dbl)
    gridfunctions(H_on_supportfns)%type = sf
    if(flag_basis_set==blips) then
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
       current_fn_on_grid = 3
    else
       current_fn_on_grid = 2
    end if
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
!!  SOURCE
!!  
  subroutine dissociate_fn_on_grid

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: flag_basis_set, blips, area_index
    use memory_module, ONLY: reg_dealloc_mem, type_dbl

    implicit none

    integer :: stat

    ! Support functions
    call reg_dealloc_mem(area_index,size(gridfunctions(supportfns)%griddata),type_dbl)
    deallocate(gridfunctions(supportfns)%griddata)
    ! H acting on support functions 
    call reg_dealloc_mem(area_index,size(gridfunctions(H_on_supportfns)%griddata),type_dbl)
    deallocate(gridfunctions(H_on_supportfns)%griddata)
    if(flag_basis_set==blips) then
       ! Non-local projector functions
       call reg_dealloc_mem(area_index,size(gridfunctions(pseudofns)%griddata),type_dbl)
       deallocate(gridfunctions(pseudofns)%griddata)
    end if
    current_fn_on_grid = 0
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
!!  SOURCE
!!  
  integer function allocate_temp_fn_on_grid(type)

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_index
    use memory_module, ONLY: reg_alloc_mem, type_dbl
    use numbers, ONLY: zero

    implicit none

    ! Passed variables
    integer :: type
    ! Local variables
    integer :: size, stat, i

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
!!  SOURCE
!!  
  subroutine free_temp_fn_on_grid(num)

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_index
    use memory_module, ONLY: reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    integer :: num

    if(num/=current_fn_on_grid) call cq_abort("Out-of-order deallocation of temp grid fn",num)
    call reg_dealloc_mem(area_index,gridfunctions(current_fn_on_grid)%size,type_dbl)
    deallocate(gridfunctions(current_fn_on_grid)%griddata)
    gridfunctions(current_fn_on_grid)%size = 0
    gridfunctions(current_fn_on_grid)%type = 0
    current_fn_on_grid = current_fn_on_grid - 1
  end subroutine free_temp_fn_on_grid
!!***

end module functions_on_grid
