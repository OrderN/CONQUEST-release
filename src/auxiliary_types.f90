! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module auxiliary_types
! ------------------------------------------------------------------------------
! Code area : 
! ------------------------------------------------------------------------------

!****h* Conquest/auxiliary_types
! NAME
!  auxiliary_types
! PURPOSE
!  
! AUTHOR
!  M.Arita
! CREATION DATE
!  2013/12/05
! MODIFICATION HISTORY
! SOURCE
!
module auxiliary_types
  ! Module usage
  use datatypes
  implicit none

  !****s* auxiliary_types/
  ! NAME
  !  group_aux
  ! PURPOSE
  !  Define all variables related to auxiliary groups
  ! AUTHOR
  ! SOURCE
  !
  type group_aux
   ! Schalor
   integer :: n_grp                        ! No. of groups
   character(20) :: filename               ! auxiliary file name
   ! Arrays
   character(20),pointer :: grp_name(:)    ! group name
   integer,pointer :: n_atom_in_grp(:)     ! No. of atoms composing group
   integer,pointer :: n_subgrp(:)          ! no. of subgroups for each group
   integer,pointer :: iatom_beg(:)         ! initial address for glob_atom w.r.t subgroups
   integer,pointer :: ibeg_grp(:)          ! initial address for glob_atom w.r.t groups
   integer,pointer :: glob_atom(:)         ! global labels of the atoms in the group
  end type group_aux
  
end module auxiliary_types
