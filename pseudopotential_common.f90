!!****h* Conquest/pseudopotential_common
!!  NAME
!!   pseudopotential_common
!!  PURPOSE
!!   Contains data which are common to old pseudo and tm pseudo
!!   This module should be used in pseudo_tm_module and 
!!   pseudopotential_data.
!!  USES
!!   datatypes, maxima_module
!!  AUTHOR
!!   Tsuyoshi MIYAZAKI
!!  CREATION DATE
!!   15/11/2002
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
module pseudopotential_common
  use datatypes, ONLY: double

 implicit none
 save

  public :: pseudo_type, non_local, &
            pseudopotential, core_radius, &
            OLDPS, SIESTA, STATE, ABINIT, core_correction

  ! RCS tag for object file identification
  character(len=80), private :: RCSid = "$Id$"

  logical :: non_local
  logical :: flag_angular_new

  integer :: pseudo_type
  integer, parameter:: OLDPS=0, SIESTA=1, STATE=2, ABINIT=3
  
  real(double), allocatable, dimension(:) :: pseudopotential
  real(double), allocatable, dimension(:) :: core_radius
  real(double) :: core_correction

end module pseudopotential_common
!!***
