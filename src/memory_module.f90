! -*- mode: F90; mode: font-lock -*-
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
!!   2014/01/17 lat 
!!    Added optional what and lun arguments in:
!!      reg_alloc_mem
!!      reg_dealloc_mem
!!      write_mem_use 
!!  SOURCE
!!
module memory_module

  use datatypes
  use global_module, only: io_lun

  implicit none
  save

  integer, allocatable, dimension(:) :: max_alloc_area, tot_alloc_area
  integer :: tot_alloc, max_alloc

  integer, parameter, dimension(3) :: no_bytes = (/4,8,16/)
  integer, parameter :: type_int  = 1
  integer, parameter :: type_dbl  = 2
  integer, parameter :: type_cplx = 3

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
  !!  USES
  !!
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2006/09/27
  !!  MODIFICATION HISTORY
  !!   2014/01/17 lat 
  !!    Added optional printing options
  !!  SOURCE
  !!
  subroutine reg_alloc_mem(area, amount, type, what, lun)

    use units,         only: m_units, mem_units, mem_conv
    use global_module, only: iprint_gen
    use GenComms,      only: myid

    implicit none

    ! Passed variables
    integer, intent(in) :: area
    integer, intent(in) :: amount
    integer, intent(in) :: type
    character(*), optional, intent(in) :: what
    integer,      optional, intent(in) :: lun

    if (iprint_gen>4 .and. present(what) .and. present(lun)) then
       write(lun,fmt='(10x,"Allocating ",a12," in area ",i3,f10.3," ",a2)') &
            adjustl(what), area, real(amount*no_bytes(type),double)*mem_conv, &
            mem_units(m_units)
    end if

    if (iprint_gen>4.AND.myid==0) then
       if (present(what)) then
          write(io_lun,fmt='(10x,"Allocating ",a12," in area ",i3,f10.3," ",a2)') &
               adjustl(what), area, real(amount*no_bytes(type),double)*mem_conv, &
               mem_units(m_units)
       else
          write(io_lun,fmt='(10x,"Allocating in area ",i3,f10.3," ",a2)') &
               area, real(amount*no_bytes(type),double)*mem_conv,mem_units(m_units)
       end if
    end if
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
  !!  USES
  !!
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2006/09/27
  !!  MODIFICATION HISTORY
  !!   2014/01/17 lat 
  !!    Added optional printing options
  !!  SOURCE
  !!
  subroutine reg_dealloc_mem(area, amount, type, what, lun)

    use units, only: m_units, mem_units, mem_conv
    use global_module, only: iprint_gen

    implicit none

    ! Passed variables
    integer :: area
    integer :: amount
    integer :: type

    character(*), optional, intent(in) :: what
    integer,      optional, intent(in) :: lun

    if (iprint_gen>4 .and. present(what).and.present(lun)) then
       write(lun,fmt='(10x,"Delocating ",a12," in area ",i3,f10.3," ",a2)') &
            adjustl(what), area, real(amount*no_bytes(type),double)*mem_conv, &
            mem_units(m_units)
    end if

    tot_alloc_area(area) = tot_alloc_area(area) - amount*no_bytes(type)
    tot_alloc = tot_alloc - amount*no_bytes(type)

  end subroutine reg_dealloc_mem

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
  !!  USES
  !!
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2006/09/27
  !!  MODIFICATION HISTORY
  !!   2014/01/17 lat 
  !!    Added optional printing options
  !!  SOURCE
  !!
  subroutine init_reg_mem

    use global_module, only: n_areas
    use GenComms,      only: cq_abort

    implicit none

    ! Local variables
    integer :: stat

    allocate(max_alloc_area(n_areas), tot_alloc_area(n_areas), STAT=stat)
    if (stat /= 0) call cq_abort("Error allocating max_alloc_area: ", n_areas)
    max_alloc_area = 0
    tot_alloc_area = 0
    max_alloc = 0
    tot_alloc = 0

    return
  end subroutine init_reg_mem
  !!***
  subroutine write_mem_use(lun, area)

    use datatypes
    use global_module, only: n_areas, iprint
    use units,         only: m_units, mem_units, mem_conv
    use GenComms,      only: inode, ionode

    implicit none

    integer :: i
    integer, optional, intent(in) :: lun
    integer, optional, intent(in) :: area

    if (present(lun).and.present(area)) then
       write(lun,'(4x,"Max mem use for area ",i4," is ",f10.3," ",a2)') &
            area,real(max_alloc_area(area))*mem_conv,mem_units(m_units)
       !write(io_lun,'(4x,"Max mem use for area ",i4," is ",f10.3," ",a2)') &
       !     area,real(max_alloc_area(area))*mem_conv,mem_units(m_units)
    else
       if(inode==ionode) then
          if(iprint>2) then
             do i=1,n_areas
                write(io_lun,'(4x,"Max mem use for area ",i4," is ",f10.3," ",a2)') &
                     i,real(max_alloc_area(i))*mem_conv,mem_units(m_units)
             end do
          end if
          write(io_lun,'(/4x,"Max total mem use is         ",f10.3," ",a2)') &
               real(max_alloc)*mem_conv,mem_units(m_units)
       end if
    end if

    return
  end subroutine write_mem_use

end module memory_module
