! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: gto_format.f90 lat $
! ------------------------------------------------------------------------------
! Module gto_format
! ------------------------------------------------------------------------------
! Code area 11: basis functions
! ------------------------------------------------------------------------------

!!****h* Conquest/gto_format *
!!  NAME
!!   gto_format
!!  PURPOSE
!!   Creates a defined type to hold data for GTO-atomic orbitals,
!!  USES
!!   datatypes
!!  AUTHOR
!!   Lionel Truflandier 
!!  CREATION DATE
!!   20/04/2012
!!  MODIFICATION HISTORY
!!  SOURCE
module gto_format

  use datatypes 

  implicit none
  save


  character(len=80), private :: RCSid = "$Id: gto_format.f90 lat $"

  type shell
     character(len=2) :: t    ! "type" (= s,p,sp,d,f...)
     integer          :: n    ! number of gto (= number of primitive functions for one "type" basis function)
     real(double)     :: f    ! scaling factor     

     ! Cartesian Gaussian type orbital info
     integer :: nn ! total angular momentum (nn=nx+ny+nz)
     integer :: nf ! number of "type" functions

     integer, allocatable :: nx(:) ! x angular momentum index (also called l)
     integer, allocatable :: ny(:) ! y ...  (...m)
     integer, allocatable :: nz(:) ! z ...  (...n)
     real(double), allocatable :: norm(:) 

     character(len=16), allocatable :: nt(:) ! "type" (= s,px,py...)

     ! Gaussian primitive functions info
     real(double), allocatable :: ak(:)   ! orbital exponent of index k (also called zeta_k)
     real(double), allocatable :: dk(:)   ! contraction coefficient of index k (also called c_k)
     real(double), allocatable :: Nk(:,:) ! norm of the primitive of index k 
  end type shell

  type atom_gto
     character(len=2) :: label  ! atom label
     integer          :: nshell ! number of shell
     integer          :: nang   ! total number of angular functions
     integer          :: nprim  ! total number of primitives
     type(shell), allocatable :: shell(:)  ! shell derived type
  end type atom_gto

  type(atom_gto), allocatable, dimension(:) :: gto

end module gto_format
!
!!****h* Conquest/gto_format_new *
!!  NAME
!!   gto_format
!!  PURPOSE
!!   Creates a defined type to hold data for GTO-atomic orbitals
!!   using the pao_format as reference
!!  USES
!!   datatypes
!!  AUTHOR
!!   Lionel Truflandier 
!!  CREATION DATE
!!   20/04/2012
!!  MODIFICATION HISTORY
!!  SOURCE
module gto_format_new

  use datatypes 

  implicit none

  save

  type gto_zeta
     integer      :: n    ! principal quantum number
     integer      :: kind ! regular zeta or polarisation
     real(double) :: occ  ! occ. of the state when generated PAOs
     !
     integer      :: ngto ! number of GTO primitives
     ! Gaussian primitive functions info
     real(double), pointer :: a(:)   ! orbital exponent 
     real(double), pointer :: d(:)   ! contraction coefficient
     real(double), pointer :: c(:)   ! origin 
  end type gto_zeta

  type sph_hrmnc
     !
     integer :: size
     ! exponent
     integer, pointer :: nx(:)
     integer, pointer :: ny(:)
     integer, pointer :: nz(:)
     ! coefficient     
     real(double), pointer ::  c(:)
     ! name
     character(len=12) :: nt
  end type sph_hrmnc
  
  type angmom_gto
     !
     integer          :: l_value ! (0,1,2...)
     character(len=1) :: l_name  ! (S,P,D...)
     integer          :: n_zeta_in_angmom
     type(gto_zeta),  pointer, dimension(:) :: zeta
     !
     ! cartesian Gaussian
     ! 
     integer :: nf_gto          ! number of l-type function
     integer, pointer :: nx(:) ! x angular momentum index (also called l)
     integer, pointer :: ny(:) ! y ...  (...m)
     integer, pointer :: nz(:) ! z ...  (...n)
     real(double), pointer :: norm(:)
     character(len=12), pointer :: nt(:) !  (s,px,py...)
     !
     ! real spherical harmonics 
     !
     integer :: nf_sph
     type(sph_hrmnc),  pointer :: transform_sph(:)
     real(double),     pointer :: norm_sph(:)
     !
  end type angmom_gto
  
  type gto_to_sf
     !
     ! cartesian Gaussian
     !      
     integer  :: ngto
     integer  :: nx ! x angular momentum index (also called l)
     integer  :: ny ! y ...  (...m)
     integer  :: nz ! z ...  (...n)
     !
     real(double)      :: norm
     character(len=12) :: nt      !  (s,px,py...)
     !
     ! real spherical harmonics
     !
     integer           :: sph_size 
     integer, pointer  :: sph_nx(:) ! x angular momentum index (also called l)
     integer, pointer  :: sph_ny(:) ! y ...  (...m)
     integer, pointer  :: sph_nz(:) ! z ...  (...n)
     real(double),      pointer :: sph_c(:)
     
     ! GTO data
     real(double), pointer, dimension(:) :: a   ! orbital exponent 
     real(double), pointer, dimension(:) :: d   ! contraction coefficient
     real(double), pointer, dimension(:) :: c   ! other if needed
  end type gto_to_sf


  type species_gto
     character(len=2) :: label
     integer :: greatest_angmom
     integer :: n_zeta_tot
     integer :: nsf_gto_tot
     integer :: nsf_sph_tot     
     type(angmom_gto), pointer, dimension(:) :: angmom
     type(gto_to_sf),  pointer, dimension(:) :: sf
  end type species_gto

  type(species_gto), allocatable, dimension(:) :: gto
  !
end module gto_format_new

!!***
