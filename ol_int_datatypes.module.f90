! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module ol_int_datatypes
! ------------------------------------------------------------------------------
! Code area 11: basis functions
! ------------------------------------------------------------------------------

!!****h* Conquest/ol_int_datatypes
!!  NAME
!!   ol_int_datatypes
!!  PURPOSE
!!   store data structures for overlap integral evaluation
!!   tables : soon to be redundant generic table type, with
!!   1) npnts - no of points, 2)del_x - spacing of points,
!!   3) arr_vals : array containing table values
!!
!!   ol_integral : data type for storing radial tables which
!!   are obtained for calculation of overlap integral
!!   1)rad_tbls - set of tables corresponding to an l1,l2,l3 triplet
!!   2)l_values : array containing l3 values for each of rad_tbls
!!   3)no_of_lvals : multiplicity of radial tables
!! 
!!   ang_coeff : data type for storing integrals of spherical
!!   harmonic triple products
!!   1) n_m_combs : no of m combinations corresponding to given l1,l2,l3
!!   2) m_combs : pointer containing numerical values of integrals
!!
!!   rad_tables : array to hold values of radial tables
!!
!!   ol_index : array to index the position of a given set of radial tables within
!!   rad_tables
!!
!!   coefficients: array to hold values of triple spherical harmonic integrals.
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!! 
!!  SOURCE
!!

module ol_int_datatypes
  use datatypes
  implicit none
  save
  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id$"


  !datatypes used in calculating overlap integrals

    type tables
       integer :: npnts
       real(double) :: del_x
       real(double) , pointer, dimension(:) :: arr_vals
       !RC 18/09/03 adding tables to hold 2nd derivatives
       real(double), pointer, dimension(:) :: arr_vals2
    end type tables
    
    type ol_integral
       type(tables), pointer, dimension(:) :: rad_tbls
       integer, pointer, dimension(:) :: l_values
       integer :: no_of_lvals
    end type ol_integral
    
    type ang_coeff
       integer :: n_m_combs
       real(double), pointer, dimension(:) :: m_combs
    end type ang_coeff
    
    integer, allocatable, dimension(:,:,:,:,:,:) :: ol_index
    integer, allocatable, dimension(:,:,:,:,:,:) :: ol_index_nlpf_pao
    integer, allocatable, dimension(:,:,:,:,:,:) :: ol_index_napf_pao ! NA projectors
    integer, allocatable, dimension(:,:,:,:,:,:) :: ol_index_ke
 
    type(ang_coeff), allocatable, dimension(:) :: coefficients

    type(ol_integral), allocatable, dimension(:) :: rad_tables
    type(ol_integral), allocatable, dimension(:) :: rad_tables_nlpf_pao
    type(ol_integral), allocatable, dimension(:) :: rad_tables_ke
    ! Neutral atom Projector functions
    ! One- and two-centre integrals i=k and/or j=k <phi_i|V_k|phi_j>
    type(ol_integral), allocatable, dimension(:) :: rad_tables_paoNApao
    ! PAO-Projector for k/=i and k/=j
    type(ol_integral), allocatable, dimension(:) :: rad_tables_napf_pao

end module ol_int_datatypes
!!***
