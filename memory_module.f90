! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module memory_module
! ------------------------------------------------------------------------------

!!****h* Conquest/memory_module *
!!  NAME
!!   memory_module
!!  PURPOSE
!!   Store details of memory allocated by areas of code
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2006/09/27
!!  MODIFICATION HISTORY
!!   2008/02/06 08:24 dave
!!    Changed for output to file not stdout
!!   2011/12/01 L.Tong
!!    Added type_cplx, and no_bytes associated to it
!!  SOURCE
!!
module memory_module

  use datatypes
  use global_module, ONLY: io_lun

  implicit none
  save

  integer, allocatable, dimension(:) :: max_alloc_area, tot_alloc_area
  integer :: tot_alloc, max_alloc

  integer, parameter, dimension(3) :: no_bytes = (/4,8,16/)
  integer, parameter :: type_int = 1
  integer, parameter :: type_dbl = 2
  integer, parameter :: type_cplx = 3
  
  interface reg_dealloc_mem
     module procedure reg_dealloc_mem_32
     module procedure reg_dealloc_mem_64
  end interface

  ! RCS tag for object file identification
  character(len=80), private :: RCSid = "$Id$"
!!***

contains

! -----------------------------------------------------------
! Subroutine reg_alloc_mem
! -----------------------------------------------------------

!!****f* memory_module/reg_alloc_mem *
!!
!!  NAME 
!!   reg_alloc_mem
!!  USAGE
!! 
!!  PURPOSE
!!   Registers the allocation of memory in a particular code area
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2006/09/27
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine reg_alloc_mem(area, amount, type)

    use units, ONLY: m_units, mem_units, mem_conv
    use global_module, ONLY: iprint_gen
    use GenComms, ONLY: myid

    implicit none

    ! Passed variables
    integer, intent(in) :: area
    integer, intent(in) :: amount
    integer, intent(in) :: type

    if(iprint_gen>4.AND.myid==0) write(io_lun,fmt='(10x,"Allocating in area ",i3,f10.3," ",a2)') &
         area, real(amount*no_bytes(type),double)*mem_conv,mem_units(m_units)
    tot_alloc_area(area) = tot_alloc_area(area) + amount*no_bytes(type)
    max_alloc_area(area) = max(max_alloc_area(area),tot_alloc_area(area))
    tot_alloc = tot_alloc + amount*no_bytes(type)
    max_alloc = max(max_alloc, tot_alloc)

  end subroutine reg_alloc_mem
!!***

! -----------------------------------------------------------
! Subroutine reg_dealloc_mem
! -----------------------------------------------------------

!!****f* memory_module/reg_dealloc_mem *
!!
!!  NAME 
!!   reg_dealloc_mem
!!  USAGE
!! 
!!  PURPOSE
!!   Registers the deallocation of memory in a particular code area
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2006/09/27
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine reg_dealloc_mem_32(area, amount,type)

    implicit none

    ! Passed variables
    integer :: area
    integer :: amount
    integer :: type

    tot_alloc_area(area) = tot_alloc_area(area) - amount*no_bytes(type)
    tot_alloc = tot_alloc - amount*no_bytes(type)

  end subroutine reg_dealloc_mem_32

  subroutine reg_dealloc_mem_64(area, amount,type)

    implicit none

    ! Passed variables
    integer :: area
    integer(8):: amount
    integer :: type

    tot_alloc_area(area) = tot_alloc_area(area) - amount*no_bytes(type)
    tot_alloc = tot_alloc - amount*no_bytes(type)

  end subroutine reg_dealloc_mem_64
!!***

! -----------------------------------------------------------
! Subroutine init_reg_mem
! -----------------------------------------------------------

!!****f* memory_module/init_reg_mem *
!!
!!  NAME 
!!   init_reg_mem
!!  USAGE
!! 
!!  PURPOSE
!!   Initialises variables used in tracking memory allocation
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2006/09/27
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine init_reg_mem

    use global_module, ONLY: n_areas
    use GenComms, ONLY: cq_abort

    implicit none

    ! Local variables
    integer :: stat

    allocate(max_alloc_area(n_areas), tot_alloc_area(n_areas), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating max_alloc_area: ",n_areas)
    max_alloc_area = 0
    tot_alloc_area = 0
    max_alloc = 0
    tot_alloc = 0
    return
  end subroutine init_reg_mem
!!***

  subroutine write_mem_use

    use datatypes
    use global_module, ONLY: n_areas
    use units, ONLY: m_units, mem_units, mem_conv
    use GenComms, ONLY: inode, ionode

    implicit none

    integer :: i

    if(inode==ionode) then
       do i=1,n_areas
          write(io_lun,'(2x,"Max mem use for area ",i4," is ",f10.3," ",a2)')&
               i,real(max_alloc_area(i))*mem_conv,mem_units(m_units)
       end do
       write(io_lun,'(2x,"Max total mem use is ",f10.3," ",a2)')&
            real(max_alloc)*mem_conv,mem_units(m_units)
    end if
  end subroutine write_mem_use

end module memory_module
