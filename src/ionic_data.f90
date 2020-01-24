! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module ionic_data
! ------------------------------------------------------------------------------
! Code area 1: initialisation
! ------------------------------------------------------------------------------

!!****h* Conquest/ionic_data *
!!  NAME
!!   ionic_data
!!  PURPOSE
!!   Collects the reading and processing of ionic data into one module
!!  USES
!!
!!  AUTHOR
!!   M.J.Gillan
!!  CREATION DATE
!!   June 2002
!!  MODIFICATION HISTORY
!!   14:29, 17/09/2002 drb 
!!    Added module header and a check for gauss/pao before calling read_pao
!!   13:56, 2003/09/22 dave
!!    Added read for basis set: defaults to blips
!!   2008/02/06 08:08 dave
!!    Changed for output to file not stdout
!!   2015/06/08 lat
!!    Added experimental backtrace
!!  SOURCE
!!
module ionic_data

  use global_module, only: io_lun
  use timer_module,  only: start_timer, stop_timer, cq_timer 
  use timer_module,  only: start_backtrace, stop_backtrace

  implicit none

  ! Area identification
  integer, parameter, private :: area = 1

!!***

contains

  !!****f* ionic_data/get_ionic_data *
  !!
  !!  NAME 
  !!   get_ionic_data
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Gets ionic data - will read PAOs, charge density and 
  !!   pseudopotentials
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   M.J.Gillan
  !!  CREATION DATE
  !!   June 2002
  !!  MODIFICATION HISTORY
  !!   14:31, 17/09/2002 drb 
  !!    Added headers and a check for PAOs before reading them
  !!   11:16, 24/09/2002 mjg & drb 
  !!    Added lines to check on whether PAOS are needed
  !!   12:01, 24/09/2002 mjg & drb 
  !!    Added calls to make or read atomic densities (and used
  !!    atomic_density module)
  !!   15:33, 25/09/2002 mjg & drb 
  !!    Added flag_no_atomic_densities to specify that atomic
  !!    densities have not been made
  !!   07:56, 2003/09/22 dave
  !!    Added flag to choose basis set
  !!   08:54, 26/05/2005 dave 
  !!    Added various bits to stop PAO reading if SIESTA pseudos are
  !!    used - the PAOs are now read during .ion file reading.
  !!    **NOTE** we assume that the FIRST call to init_pseudo_siesta
  !!    is during read_and_write and that the tables are read in there
  !!   10:55, 13/02/2006 drb 
  !!    Line length change
  !!   16:00, 2017/02/21 nakata
  !!    Commented out get_support_pao_rep which is no longer used
  !!   15:00, 2017/03/08 nakata
  !!    Removed get_support_pao_rep which is no longer used
  !!  SOURCE
  !!
  subroutine get_ionic_data(inode,ionode,level)

    use datatypes
    use atomic_density,         only: make_atomic_density_from_paos, &
                                      spline_atomic_density
    use species_module,         only: n_species
    use GenComms,               only: cq_abort

    implicit none

    ! Passed variables
    integer, optional    :: level
    integer, intent(in)  :: inode, ionode

    ! Local variables
    character(len=10) :: init_blip_method
    logical           :: flag_blips_from_pao
    type(cq_timer)    :: backtrace_timer
    integer           :: backtrace_level

!****lat<$    
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='get_ionic_data',&
         where=area,level=backtrace_level)
!****lat>$

    call make_atomic_density_from_paos(inode, ionode, n_species)
    call spline_atomic_density(n_species)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_ionic_data')
!****lat>$

    return
  end subroutine get_ionic_data
  !!***

end module ionic_data





