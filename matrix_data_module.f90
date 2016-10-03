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
!!   2008/03/12 06:12 dave
!!    Added max_range parameter to find longest range matrix
!!   2012/01/18 16:55 dave
!!    Added blip transfer type for analytic blip integration
!!   2014/01/17 13:20 lat
!!    Added new matrix X and SX for exchange
!!   2016/08/09 18:00 nakata
!!    Added parameters Satomf_range, Hatomf_range, AP_range, PA_range, HTr_range,
!!                     SFcoeff_range, SFcoeffTr_range, Satomfsf_range, Hatomfsf_range, LD_range,
!!                     SFcoeff_matind, SFcoeffTr_matind, Satomf_matind, Hatomf_matind, 
!!                     SPatomf_matind, PSatomf_matind, Satomfsf_matind, Hatomfsf_matind, HTr_matind, LD_matind
!!    for PAO-based calculations
!!  SOURCE
!!
module matrix_data

  use datatypes
  use matrix_module, only: matrix, matrix_halo, blip_transfer ! Think about this - how do we want to arrange ?
  
  implicit none
  save

  ! This will need to change if the above parameters are changed
  integer, parameter :: mx_matrices = 32   ! nakata3

  ! Store ALL indices in a large array
  type(matrix),      allocatable, dimension(:,:), target :: mat
  type(matrix_halo), allocatable, dimension(:),   target :: halo
  real(double),      dimension(mx_matrices) :: rcut
  character(len=3),  dimension(mx_matrices) :: mat_name

  integer, dimension(:), pointer :: Smatind, dSmatind, Lmatind, LTrmatind, Hmatind, dHmatind, &
       SPmatind, PAOPmatind, PSmatind, LSmatind, SLmatind, LHmatind, HLmatind, LSLmatind,     &
       SLSmatind, Tmatind, TTrmatind, TSmatind, THmatind, TLmatind, Xmatind, SXmatind   
  integer, dimension(:), pointer :: SFcoeff_matind, SFcoeffTr_matind, Satomf_matind, Hatomf_matind, &
                                    SPatomf_matind, PSatomf_matind, Satomfsf_matind, Hatomfsf_matind,   &
                                    HTr_matind, LD_matind

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
  integer, parameter :: HLrange  = 20

  integer, parameter :: Xrange   = 21
  integer, parameter :: SXrange  = 22

!!! nakata3
  ! The indices for ATOMF-based-matrix ranges will be set later.
  integer :: Satomf_range    ! 23
  integer :: Hatomf_range    ! 24
  integer :: AP_range        ! 25
  integer :: PA_range        ! 26
  integer :: HTr_range       ! 27
  integer :: SFcoeff_range   ! 28
  integer :: SFcoeffTr_range ! 29
  integer :: Satomfsf_range  ! 30
  integer :: Hatomfsf_range  ! 31
  integer :: LD_range        ! 32
!!! end nakata3

  integer :: max_range ! Indexes matrix with largest range

  type(blip_transfer) :: blip_trans
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
