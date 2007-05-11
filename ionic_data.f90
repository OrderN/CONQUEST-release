! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!  SOURCE
!!
module ionic_data

  implicit none

  ! RCS tag for object file identification
  character(len=80), private, save :: RCSid = "$Id$"
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
!!    Added calls to make or read atomic densities (and used atomic_density module)
!!   15:33, 25/09/2002 mjg & drb 
!!    Added flag_no_atomic_densities to specify that atomic densities have not
!!    been made
!!   07:56, 2003/09/22 dave
!!    Added flag to choose basis set
!!   08:54, 26/05/2005 dave 
!!    Added various bits to stop PAO reading if SIESTA pseudos are used - the PAOs are now
!!    read during .ion file reading.  **NOTE** we assume that the FIRST call to init_pseudo_siesta
!!    is during read_and_write and that the tables are read in there
!!   10:55, 13/02/2006 drb 
!!    Line length change
!!  SOURCE
!!
  subroutine get_ionic_data(inode,ionode,flag_no_atomic_densities)

    use datatypes
    use read_pao_info, ONLY : read_pao
    use atomic_density, ONLY: read_atomic_density, make_atomic_density_from_paos, spline_atomic_density, &
         flag_atomic_density_from_pao, atomic_density_method
    use species_module, ONLY : n_species
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: flag_basis_set, blips, PAOs, iprint_init
    use make_rad_tables, ONLY : get_support_pao_rep
    use pseudopotential_common, ONLY: pseudo_type, SIESTA, ABINIT
    use blip, ONLY: init_blip_flag

    implicit none

    ! Passed variables
    integer, intent(IN) :: inode, ionode
    logical, intent(OUT) :: flag_no_atomic_densities

    ! Local variables
    character(len=10) :: init_blip_method
    logical :: flag_blips_from_pao!, flag_atomic_density_from_pao
    logical, external :: leqi ! External subroutine from Siesta FDF
    
    ! Decide whether blips are to be initialised from PAOs
    flag_blips_from_pao = .false.
    if(leqi(init_blip_flag,'pao')) flag_blips_from_pao = .true.
    ! Decide whether atomic densities are to be initialised from PAOs
    flag_atomic_density_from_pao = .false.
    if(leqi(atomic_density_method,'pao')) flag_atomic_density_from_pao = .true.
    ! If PAOs are needed, then read them
    ! Note that we now read PAOs from the .ion file if we're using that pseudopotential
    if((pseudo_type/=SIESTA.AND.pseudo_type/=ABINIT).AND. &
         (flag_blips_from_pao.OR.flag_atomic_density_from_pao.OR.flag_basis_set==PAOs)) then
       if(inode == ionode.AND.iprint_init>1) &
            write(unit=*,fmt='(//" **** get_ionic_data: about to call read_pao_info")')
       call read_pao(inode,ionode,n_species)
    else if((pseudo_type==SIESTA.OR.pseudo_type==ABINIT).AND. &
         (flag_blips_from_pao.OR.flag_atomic_density_from_pao.OR.flag_basis_set==PAOs)) then
       if(inode == ionode.AND.iprint_init>1) &
            write(unit=*,fmt='(//" **** get_ionic_data: PAO input from init_pseudo_tm already done")')
    end if
    ! If we're using PAOs as a basis, need coefficients
    if(flag_basis_set==PAOs) call get_support_pao_rep(inode,ionode)
    ! Get the atomic densities
    if(flag_atomic_density_from_pao) then ! Use PAOs for atomic densities
       call make_atomic_density_from_paos(inode,ionode,n_species)
    else if(.NOT.flag_atomic_density_from_pao) then
       if(leqi(atomic_density_method,'read')) then ! User has specified a file
          call read_atomic_density(inode,ionode,n_species)
       else ! Signal to the calling routine that no atomic densities have been created
          flag_no_atomic_densities = .true.
       end if
    end if
    if(.NOT.flag_no_atomic_densities) then ! Spline tables
       call spline_atomic_density(n_species)
    end if
    return
  end subroutine get_ionic_data
!!***

end module ionic_data





