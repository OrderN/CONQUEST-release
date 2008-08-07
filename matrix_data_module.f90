! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module matrix_data
! ------------------------------------------------------------------------------
! Code area 2: matrices
! ------------------------------------------------------------------------------

!!****h* Conquest/matrix_data *
!!  NAME
!!   matrix_data
!!  PURPOSE
!!   Gives access to ALL the indices (i.e. the matrix
!!   types) and data (i.e. the elements) for the matrices in the system.
!!  USES
!!   datatypes, matrix_module, maxima_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   12/04/00
!!  MODIFICATION HISTORY
!!   18/03/2002 dave
!!    Added ROBODoc header, RCS Id and Log tags and static for object file id
!!   29/07/2002 drb
!!    Changed temp_T so that is has dimension (NSF,NSF,....) instead of (NCF,NCF,....)
!!    to fix a problem when NCF/=NSF
!!   31/07/2002 dave
!!    Added data_M12 (for Pulay force in exact diagonalisation); removed misleading 
!!    statement in purpose tag above
!!   14:36, 26/02/2003 drb 
!!    Reordered and removed integ references during bug search
!!   15:01, 12/03/2003 drb 
!!    Added data_KE and data_NL
!!   2006/01/25 16:57 dave
!!    Starting to add new variables and routines for encapsulated matrix data
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines: specifically declared matrix_pointer type, made dS and dH not allocatable and
!!    all data_ variables targets
!!   2006/08/30 08:12 dave
!!    Finalised conversion of matrix indexing to allocatable
!!  SOURCE
!!
module matrix_data

  use datatypes
  use matrix_module, ONLY: matrix, matrix_halo ! Think about this - how do we want to arrange ?
  
  implicit none
  save

  ! This will need to change if the above parameters are changed
  integer, parameter :: mx_matrices = 20

  ! Store ALL indices in a large array
  type(matrix), allocatable, dimension(:,:), target :: mat
  type(matrix_halo), allocatable, dimension(:), target :: halo
  real(double) :: rcut(mx_matrices)
  character(len=3), dimension(mx_matrices) :: mat_name

  integer, dimension(:), pointer :: Smatind, dSmatind, Lmatind, LTrmatind, Hmatind, dHmatind, &
       SPmatind, PAOPmatind, PSmatind, LSmatind, SLmatind, LHmatind, HLmatind, LSLmatind, &
       SLSmatind, Tmatind, TTrmatind, TSmatind, THmatind, TLmatind  

  ! Parameters for the different matrix ranges
  integer, parameter :: Srange   = 1   ! STS,TST,TS.TS
  integer, parameter :: Lrange   = 2   ! SLS,HLSLS,SLSLS,THT,LSLSL
  integer, parameter :: Hrange   = 3   ! K,LSLSL
  integer, parameter :: SPrange  = 4   ! PS,U   ! Maybe rename ?
  integer, parameter :: LSrange  = 5   ! SL
  integer, parameter :: LHrange  = 6   ! HL
  integer, parameter :: LSLrange = 7   ! HLS
  integer, parameter :: SLSrange = 8   ! LHL
  integer, parameter :: Trange   = 9   ! T is S^-1
  integer, parameter :: TSrange  = 10  
  integer, parameter :: THrange  = 11
  integer, parameter :: TLrange  = 12
  integer, parameter :: PSrange  = 13   ! PS,U   ! Maybe rename ?
  integer, parameter :: LTrrange = 14
  integer, parameter :: SLrange  = 15
  integer, parameter :: TTrrange = 16
  integer, parameter :: dSrange  = 17   ! STS,TST,TS.TS
  integer, parameter :: dHrange  = 18   ! K,LSLSL
  integer, parameter :: PAOPrange= 19  ! PS,U   ! Maybe rename ?
  integer, parameter :: HLrange = 20
!!***

!!****s* multiply_module/matrix_pointer *
!!  NAME
!!   matrix_pointer
!!  PURPOSE
!!   Contains pointer to matrix data and length of matrix
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type matrix_pointer
     integer :: length
     integer :: sf1_type, sf2_type
     real(double), pointer, dimension(:) :: matrix
  end type matrix_pointer
!!***

  ! RCS tag for object file identification 
  character(len=80), private :: RCSid = "$Id$"

end module matrix_data
