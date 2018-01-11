! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module pseudo_tm_module
! ------------------------------------------------------------------------------
! Code area 10: pseudopotentials
! ------------------------------------------------------------------------------

!!****h* Conquest/pseudo_tm_module *
!!  NAME
!!   pseudo_tm_module
!!  PURPOSE
!!   sets up 
!!
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   Jul. 2002
!!  MODIFICATION HISTORY
!!   11:04, 2003/02/28 dave
!!    Added deallocation for structures
!!   14:03, 2003/12/19 dave & rc
!!    Added new angular routines for PAO compatibility
!!   2008/02/06 08:33 dave
!!    Changed for output to file not stdout
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2015/05/01 13:41 dave and sym
!!    Adding stress
!!  SOURCE
!!
module pseudo_tm_module

  use datatypes,              only: double
  use pseudo_tm_info,         only: pseudo_info, pseudo, setup_pseudo_info, loc_pot, loc_chg
  use pseudopotential_common, only: pseudopotential, core_radius, flag_angular_new
  use global_module,          only: iprint_pseudo, io_lun
  use timer_module,           only: start_timer, stop_timer, cq_timer
  use timer_module,           only: start_backtrace, stop_backtrace
  use timer_stdclocks_module, only: tmr_std_pseudopot, tmr_std_basis, tmr_std_allocation

  implicit none

  !arrays
  private
  public :: init_pseudo_tm, loc_pp_derivative_tm, &
       nonloc_pp_derivative_tm, set_tm_pseudo, deallocate_pseudo_tm, &
       loc_HF_stress, loc_G_stress, make_neutral_atom
  !public
  ! pseudofunctions is arranged like 
  !  (n_pts in block, ncf, naba atoms, iblock)
  !The definition of the following arrays are moved to
  ! pseudopotential_common because of consistency with
  ! old pseudopotentials.
  !real(double),allocatable,save :: pseudofunctions(:)
  !real(double),allocatable,save :: pseudopotential(:)
  !real(double),allocatable,save :: core_radius(:)
  !private
  !%%! real(double), allocatable,save :: spherical_harmonic_norm(:)
  !offset_mcomp(nl, ispecies)
  ! we now assume multi references for each l. 
  ! (in the following, we write ... as nl)
  ! offset_mcomp is prepared to calculate the 
  ! index of m-components for each nl. 
  integer, allocatable, save :: offset_mcomp(:,:)
  !the followings are just workarrays
  integer, allocatable :: ip_store(:)
  real(double), allocatable :: x_store(:)
  real(double), allocatable :: y_store(:)
  real(double), allocatable :: z_store(:)
  real(double), allocatable :: r_store(:)
  real(double), dimension(3) :: loc_HF_stress, loc_G_stress
  !interface
  !   interface nonloc_pp_derivative_tm
  !     module procedure nonloc_pp_derivative_tm
  !   end interface nonloc_pp_derivative_tm

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: &
       RCSid = "$Id$"
  !!***

contains

  ! -----------------------------------------------------------
  ! Subroutine init_pseudo_tm
  ! -----------------------------------------------------------

!!****f* pseudo_tm_module/init_pseudo_tm *
!!
!!  NAME 
!!   init_pseudo_tm
!!  USAGE
!! 
!!  PURPOSE
!!   Initialises Troullier-Martin pseudopotentials
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   Jul. 2002
!!  MODIFICATION HISTORY
!!   18:25, 2003/02/27 dave
!!    Added ROBODoc header
!!   2011/09/16 11:01 dave
!!    Tidied to use species_module variables only
!!   2015/06/08 lat
!!    Added experimental backtrace
!!   2017/04/24 11:52 dave
!!    Updated rcutmax for neutral atom
!!  SOURCE
!!
  subroutine init_pseudo_tm(ecore, ncf_out) 
    
    use datatypes,      only: double
    use global_module,  only: iprint_pseudo, flag_neutral_atom
    use block_module,   only: n_pts_in_block
    use GenComms,       only: cq_abort, myid, inode, ionode
    use species_module, only: species, n_species, species_label, &
         nlpf_species
    
    implicit none

    real(double), intent(out) :: ecore

    integer :: ncf_max   ! calculated ncf
    integer :: nl_max    ! calculated n_L
    integer :: lcomp_max !
    integer :: n_mcomp   ! size of spherical_harmonic_norm
    integer :: stat
    integer :: ispecies, nlm, nl, the_l
    ! For PreConquest or calcluates maximum  16/05/2003 TM
    integer, intent(out), OPTIONAL :: ncf_out
    logical :: flag_calc_max = .false.
    type(cq_timer)    :: backtrace_timer

!****lat<$
    call start_backtrace(t=backtrace_timer,who='init_pseudo_tm',where=10,level=2)
!****lat>$

    call start_timer(tmr_std_pseudopot)
    if(PRESENT(ncf_out)) then
       flag_calc_max = .true.
    else
       flag_calc_max = .false.
    endif

    if(myid==0.AND.iprint_pseudo>2) &
         write(io_lun,fmt='(10x,"Entering init_pseudo_tm")')
    call setup_pseudo_info
    ! in the moduel pseudo_tm_info
    if(flag_calc_max) then
       call calc_max                 !contained in this subroutine
       ncf_out= ncf_max 
    else
       call check_max                !contained in this subroutine
    endif
    call allocate_pseudo_tm    !contained in this subroutine
    !%%! call setup_spherical_harmonic  !contained in this subroutine
    call setup_offset_mcomp        !contained in this subroutine
    call setup_core_radius         !contained in this subroutine
    if(flag_calc_max) then
       call deallocate_pseudo_tm  !contained in this module
       call stop_timer(tmr_std_pseudopot)    ! This is timed inside set_tm_pseudo
!****lat<$
       call stop_backtrace(t=backtrace_timer,who='init_pseudo_tm')
!****lat>$
       return
    endif

    call stop_timer(tmr_std_pseudopot)    ! This is timed inside set_tm_pseudo
    call set_tm_pseudo         !contained in this module
    call start_timer(tmr_std_pseudopot)
    call get_energy_shift(ecore)   !contained in this module
    call stop_timer(tmr_std_pseudopot)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='init_pseudo_tm')
!****lat>$

    return

  contains

    !sbrt check_max  -------------------------------
    subroutine check_max
      ncf_max = 1; nl_max=1 ; lcomp_max = 0
      do ispecies = 1, n_species
         nl_max = max(nl_max, pseudo(ispecies)%n_pjnl)
         nlm =  0  ! # of projector functions
        if(pseudo(ispecies)%n_pjnl>0) then
         do nl= 1, pseudo(ispecies)%n_pjnl
            the_l = pseudo(ispecies)%pjnl_l(nl)
            nlm = nlm + 2* the_l +1
            lcomp_max = max(lcomp_max, the_l)
         enddo
        endif
         !  (we may have multi refs. for each l)
         ncf_max = max(ncf_max, nlm)     ! # of nlm 
      enddo
      n_mcomp = (lcomp_max + 1) ** 2  ! # of lm for l = 0-lcomp_max
      return
    end subroutine check_max
    !sbrt calc_max  -------------------------------
    subroutine calc_max
      !calculates ncf_max, nl_max, lcomp_max
      if(n_species <= 0 ) &
           call cq_abort ('calc_max in init_pseudo: n_species error',&
           n_species)
      !ncf_max, nl_max, lcomp_max
      ncf_max = 1; nl_max=1 ; lcomp_max = 0
      do ispecies = 1, n_species
         nl_max = max(nl_max, pseudo(ispecies)%n_pjnl)
         nlm =  0  ! # of projector functions
        if(pseudo(ispecies)%n_pjnl > 0) then
         do nl= 1, pseudo(ispecies)%n_pjnl
            the_l = pseudo(ispecies)%pjnl_l(nl)
            nlm = nlm + 2* the_l +1
            lcomp_max = max(lcomp_max, the_l)
         enddo
        endif
         !  (we may have multi refs. for each l)
         ncf_max = max(ncf_max, nlm)     ! # of nlm
         nlpf_species(ispecies)= nlm
      enddo
      n_mcomp = (lcomp_max + 1) ** 2  ! # of lm for l = 0-lcomp_max
      return
    end subroutine calc_max

    !sbrt allocate_pseudo  -------------------------------
    subroutine allocate_pseudo_tm
      !allocate(pseudofunctions(ncf*core_size), STAT=stat)
      ! if(stat /= 0) call cq_abort &
      !   ('ERROR! init_pseudo: alloc pseudofunctions', stat, ncf*core_size)
      !allocate(pseudopotential(N_GRID_MAX), STAT=stat)
      ! if(stat /= 0) call cq_abort &
      !   ('ERROR! init_pseudo: alloc pseudopotential', stat, N_GRID_MAX)
      call start_timer(tmr_std_allocation)
      if(.NOT.allocated(core_radius)) then
         allocate(core_radius(n_species), STAT=stat)
         if(stat /= 0) &
              call cq_abort ('ERROR! init_pseudo: alloc core_radius', &
              stat, n_species)
      end if

      !%%! allocate(spherical_harmonic_norm(n_mcomp), STAT=stat)
      !%%! if(stat /= 0) call cq_abort &
      !%%!      ('ERROR! init_pseudo: alloc spherical_harmonic_norm', stat, n_mcomp)
      allocate(offset_mcomp(nl_max,n_species), STAT=stat)
      if(stat /= 0) call cq_abort &
           ('ERROR! init_pseudo: alloc offset_mcomp', stat, nl_max*n_species)
      if(.not.flag_calc_max) then
         allocate(ip_store(n_pts_in_block), STAT=stat)
         if(stat /= 0) call cq_abort &
              ('ERROR! init_pseudo: alloc ip_store', stat, n_pts_in_block)
         allocate(x_store(n_pts_in_block), STAT=stat)
         if(stat /= 0) call cq_abort &
              ('ERROR! init_pseudo: alloc x_store', stat, n_pts_in_block)
         allocate(y_store(n_pts_in_block), STAT=stat)
         if(stat /= 0) call cq_abort &
              ('ERROR! init_pseudo: alloc y_store', stat, n_pts_in_block)
         allocate(z_store(n_pts_in_block), STAT=stat)
         if(stat /= 0) call cq_abort &
              ('ERROR! init_pseudo: alloc z_store', stat, n_pts_in_block)
         allocate(r_store(n_pts_in_block), STAT=stat)
         if(stat /= 0) call cq_abort &
              ('ERROR! init_pseudo: alloc r_store', stat, n_pts_in_block)
      else
         allocate(ip_store(1), x_store(1), y_store(1), z_store(1), r_store(1),&
              STAT = stat)
         if(stat /= 0) call cq_abort &
              ('ERROR! init_pseudo: alloc for calc_max', stat)
      endif
      call stop_timer(tmr_std_allocation)
      return
    end subroutine allocate_pseudo_tm

    !%%! !sbrt setup_spherical_harmonic  -------------------------------
    !%%! subroutine setup_spherical_harmonic
    !%%!   use numbers
    !%%!   implicit none
    !%%!   ! norm for the l=0, m=0
    !%%!   if(n_mcomp >= 1) then
    !%%!      spherical_harmonic_norm(1) = one / sqrt( four * pi )
    !%%!   endif
    !%%!   ! norm for the x, y and z type functions
    !%%!   if(n_mcomp >= 4) then
    !%%!      spherical_harmonic_norm(2) = sqrt( three / ( four * pi ) )
    !%%!      spherical_harmonic_norm(3) = sqrt( three / ( four * pi ) )
    !%%!      spherical_harmonic_norm(4) = sqrt( three / ( four * pi ) )
    !%%!   endif
    !%%!   ! norm for the x^2-y^2 and 3z^2-1 type functions
    !%%!   if(n_mcomp == 9) then
    !%%!      spherical_harmonic_norm(5) = sqrt( ( three*five ) / ( four*four*pi ) )
    !%%!      spherical_harmonic_norm(6) = sqrt( five / ( four * four * pi ) )
    !%%!      ! norm for the x * y, x * z and y * z type functions
    !%%!      spherical_harmonic_norm(7) = sqrt( ( three*five ) / ( four * pi ) )
    !%%!      spherical_harmonic_norm(8) = sqrt( ( three*five ) / ( four * pi ) )
    !%%!      spherical_harmonic_norm(9) = sqrt( ( three*five ) / ( four * pi ) )
    !%%!   else if(n_mcomp/=4.AND.n_mcomp/=9) then
    !%%!      if(inode==ionode) &
    !%%!           write(io_lun,fmt='(10x,"Warning! Possible problem with ang. mom. in  pseudo_tm_module: ",i4)') n_mcomp
    !%%!   endif
    !%%!   return
    !%%! end subroutine setup_spherical_harmonic

    !sbrt setup_offset_mcomp  -------------------------------
    subroutine setup_offset_mcomp
      do ispecies = 1, n_species
         nlm =  0  ! # of projector functions
         offset_mcomp(1, ispecies) = 0
         if(pseudo(ispecies)%n_pjnl > 1) then
            do nl= 1, pseudo(ispecies)%n_pjnl-1
               the_l = pseudo(ispecies)%pjnl_l(nl)
               offset_mcomp(nl+1,ispecies) =  &
                    offset_mcomp(nl,ispecies) + 2* the_l + 1
            enddo
         endif
      enddo
      return
    end subroutine setup_offset_mcomp

    !sbrt setup_core_radius  -------------------------------
    subroutine setup_core_radius
      use numbers, only: zero
      real(double) :: rcutmax
      do ispecies = 1, n_species
         rcutmax = zero
         if(pseudo(ispecies)%tm_loc_pot==loc_pot) then
            rcutmax = max(rcutmax, pseudo(ispecies)%vlocal%cutoff)
         else
            rcutmax = max(rcutmax, pseudo(ispecies)%chlocal%cutoff)
         end if
         if(flag_neutral_atom) rcutmax = max(rcutmax, pseudo(ispecies)%vna%cutoff)
         if(pseudo(ispecies)%flag_pcc) &
              rcutmax = max(rcutmax, pseudo(ispecies)%chpcc%cutoff)
        if(pseudo(ispecies)%n_pjnl > 0) then
         do nl= 1, pseudo(ispecies)%n_pjnl
            rcutmax = max(rcutmax, pseudo(ispecies)%pjnl(nl)%cutoff)
         enddo
        endif
         core_radius(ispecies) = rcutmax
      enddo
      return
    end subroutine setup_core_radius

  end subroutine init_pseudo_tm
!!***

! -----------------------------------------------------------
! Subroutine deallocate_pseudo_tm
! -----------------------------------------------------------

!!****f* pseudo_tm_module/deallocate_pseudo_tm *
!!
!!  NAME 
!!   deallocate_pseudo_tm
!!  USAGE
!! 
!!  PURPOSE
!!   Deallocates memory for Troullier-Martin pseudopotentials
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   18:19, 2003/02/27 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine deallocate_pseudo_tm

    use GenComms, only: cq_abort

    implicit none

    integer :: stat

    !deallocate(spherical_harmonic_norm, STAT=stat)
    !if(stat /= 0) call cq_abort('ERROR! init_pseudo: dealloc spherical_harmonic_norm', stat)
    call start_timer(tmr_std_allocation)
    deallocate(offset_mcomp, STAT=stat)
    if(stat /= 0) call cq_abort('ERROR! init_pseudo: dealloc offset_mcomp', stat)
    deallocate(ip_store, STAT=stat)
    if(stat /= 0) call cq_abort('ERROR! init_pseudo: dealloc ip_store', stat)
    deallocate(x_store, STAT=stat)
    if(stat /= 0) call cq_abort('ERROR! init_pseudo: dealloc x_store', stat)
    deallocate(y_store, STAT=stat)
    if(stat /= 0) call cq_abort('ERROR! init_pseudo: dealloc y_store', stat)
    deallocate(z_store, STAT=stat)
    if(stat /= 0) call cq_abort('ERROR! init_pseudo: dealloc z_store', stat)
    deallocate(r_store, STAT=stat)
    if(stat /= 0) call cq_abort('ERROR! init_pseudo: dealloc r_store', stat)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_pseudo_tm
  !!***

! -----------------------------------------------------------
! Subroutine set_tm_pseudo
! -----------------------------------------------------------

!!****f* pseudo_tm_module/set_tm_pseudo *
!!
!!  NAME 
!!   set_tm_pseudo
!!  USAGE
!! 
!!  PURPOSE
!!   Puts the local pseudopotential on the grid points belonging
!!   to this processor and puts the non-local projector functions
!!   similarly (though in an array of size NCF*CORE_SIZE rather 
!!   than the NSF*SUPPORT_SIZE)
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   20/06/2002
!!  MODIFICATION HISTORY
!!   2007/07/30 08:06 dave
!!    Tidied up cutoff finding: re-using core_radius instead of finding max every time
!!   2008/03/03 18:50 dave
!!    Changed float to real
!!   2011/12/22 13:55 dave
!!    Removed calculation of pseudo_functions for analytic blips
!!   2016/01/07 11:58 dave
!!    Added switch for loop over neighbours for neutral atom, tidied
!!   2016/06/23 15:53 dave
!!    Bug fix: added a new loop for NLPF with blips, neutral atom but NOT analytic blip integrals
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks
!!   2016/08/08 15:30 nakata
!!    Removed unused sf in global_module
!!   2017/11/27 15:27 dave   
!!    Added if loop for NA projectors to avoid making pseudopotential on grid
!!  SOURCE
!!
  subroutine set_tm_pseudo

    use datatypes
    use numbers
    use global_module, only: rcellx,rcelly,rcellz,id_glob, ni_in_cell, &
                             iprint_pseudo, species_glob, nlpf,        &
                             flag_basis_set, blips,                    &
                             IPRINT_TIME_THRES3, flag_analytic_blip_int, &
                             flag_neutral_atom, dens
    use species_module, only: species, nlpf_species, n_species
    !  At present, these arrays are dummy arguments.
    use block_module, only : nx_in_block,ny_in_block,nz_in_block, &
         n_pts_in_block
    use group_module, only : blocks, parts
    use primary_module, only: domain
    use cover_module, only: DCS_parts
    use set_blipgrid_module, only : naba_atoms_of_blocks
    use GenBlas, only: axpy, copy
    use GenComms, only: my_barrier, cq_abort, inode, ionode, myid
    use angular_coeff_routines, only : pp_elem
    use hartree_module, only: hartree
    use functions_on_grid, only: gridfunctions, pseudofns
    use dimens, only: n_my_grid_points
    use maxima_module, only: maxngrid
    use timer_module, only: cq_timer, start_timer, stop_print_timer, WITH_LEVEL
    use pseudopotential_common, only: flag_neutral_atom_projector
    
    implicit none 
    !local
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom
    integer :: the_species
    integer :: j,iblock,the_l,ipoint, igrid
    integer :: no_of_ib_ia, offset_position
    integer :: position,iatom
    integer :: stat, nl, npoint, ip
    integer :: i,m, pseudo_neighbour

    real(double):: dcellx_block,dcelly_block,dcellz_block
    real(double):: xatom,yatom,zatom,alpha,step
    real(double):: xblock,yblock,zblock
    real(double) :: r_from_i
    real(double) :: rr,a,b,c,d,x,y,z,nl_potential
    real(double) :: coulomb_energy
    real(double) :: rcut
    real(double) :: r1, r2, r3, r4, core_charge, gauss_charge
    real(double) :: val
    real(double),allocatable :: chlocal_density(:), coulomb_potential(:)

    logical :: local_charge = .false.

    type(cq_timer) :: tmr_l_tmp1

    ! --  Start of subroutine  ---
    call start_timer(tmr_std_pseudopot)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    ! Set up neighbour list depending on whether we use VNA
    if(flag_neutral_atom) then
       pseudo_neighbour = dens
    else
       pseudo_neighbour = nlpf
    end if
    ! Local charge
    do i=1,n_species
       if(pseudo(i)%tm_loc_pot==loc_chg) local_charge = .true.
    end do
    !allocates allocatable arrays
    call start_timer(tmr_std_allocation)
    allocate(chlocal_density(maxngrid), STAT=stat)
    if(stat/=0) call cq_abort('set_ps: alloc c_d ', stat)
    if(.NOT.local_charge) then
       allocate(coulomb_potential(maxngrid), STAT=stat)
       if(stat/=0) call cq_abort('set_ps: alloc c_p ', stat)
    end if
    if(stat/=0) call cq_abort('set_ps: alloc c_d ', stat)
    call stop_timer(tmr_std_allocation)
    chlocal_density   = zero

    pseudopotential = zero
    if(flag_basis_set==blips.AND.(.NOT.flag_analytic_blip_int)) then
       gridfunctions(pseudofns)%griddata = zero
    end if


    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz

    !  This subroutine assumes Troullier-Martin's pseudopotential.
    ! For Siesta-styly, the local part of the pseudopotentials is made by the
    ! distribution of artificial core charge (local charge) density. 
    ! As this core charge distribution is smooth enough, 
    ! we don't need to consider the short range term of the local
    ! pseudopotential. 
    !      29/06/2002  Tsuyoshi MIYAZAKI
    ! For ABINIT style, the local potential is tabulated as are the projectors
    !      30/07/2007 David Bowler (belated note)

    call my_barrier()

    no_of_ib_ia = 0

    ! loop arround grid points in the domain, and for each
    ! point, get the core charge (local charge) density from the atoms 
    ! whose distances from the grid point are within the cutoff. 
    the_species = 1
    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
       if(.NOT.flag_neutral_atom_projector) then
       if(naba_atoms_of_blocks(pseudo_neighbour)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom=0
          do ipart=1,naba_atoms_of_blocks(pseudo_neighbour)%no_of_part(iblock)
             jpart=naba_atoms_of_blocks(pseudo_neighbour)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('set_ps: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atoms_of_blocks(pseudo_neighbour)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atoms_of_blocks(pseudo_neighbour)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)

                if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                   call cq_abort('set_ps: globID ERROR ', &
                        ii,parts%icell_beg(ind_part))
                endif
                if(icover > DCS_parts%mx_mcover) then
                   call cq_abort('set_ps: icover ERROR ', &
                        icover,DCS_parts%mx_mcover)
                endif

                !-- I (TM) am not sure where the next line should be.
                !   (no_of_ib_ia = no_of_ib_ia +1)
                !  This determines how to store 
                !   pseudofunctions(ipoints,ncf,naba_atm,iblock).
                !  Now, I assume we should consider all naba_atm whose 
                ! distance from the block is within the maximum of core_radius.
                ! This is needed to keep the consistency with <set_bucket>.
                ! However, we can change this strategy by changing naba_atoms_of_blocks(pseudo_neighbour).
                !no_of_ib_ia = no_of_ib_ia +1

                xatom=DCS_parts%xcover(icover)
                yatom=DCS_parts%ycover(icover)
                zatom=DCS_parts%zcover(icover)

                !the_species=species(ig_atom)
                the_species = species_glob(ig_atom)

                !calculates distances between the atom and integration grid points
                !in the block and stores which integration grids are neighbours.
                if( flag_neutral_atom ) then ! for Neutral atom potential
                   rcut = pseudo(the_species)%vna%cutoff + very_small  
                else
                   rcut = core_radius(the_species) + very_small   !!   2007/07/30 drb
                end if
                call check_block &
                     (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
                     npoint,ip_store,r_store,x_store,y_store,z_store) !out
                !r_from_i = sqrt((xatom-xblock)**2+(yatom-yblock)**2+ &
                !     (zatom-zblock)**2 )

                !Local part  ----
                ! construct the local part of the pseudopotential.
                ! local part of Siesta's pseudopotential is generated
                ! from the pseudo core charge (local charge) distribution
                if(npoint > 0) then
                   if(flag_neutral_atom) then
                      step = pseudo(the_species)%vna%delta
                      do ip=1, npoint
                         ipoint=ip_store(ip)
                         igrid=n_pts_in_block*(iblock-1)+ipoint
                         if(igrid > n_my_grid_points) call cq_abort &
                              ('set_ps: igrid error ', igrid, n_my_grid_points) 
                         r_from_i = r_store(ip)
                         x = x_store(ip)
                         y = y_store(ip)
                         z = z_store(ip)
                         j = aint( r_from_i / step ) + 1
                         if(j+1 <= pseudo(the_species)%vna%n) then
                            rr = real(j,double) * step
                            a = ( rr - r_from_i ) / step
                            b = one - a
                            c = a * ( a * a - one ) * step * step / six
                            d = b * ( b * b - one ) * step * step / six

                            r1=pseudo(the_species)%vna%f(j)
                            r2=pseudo(the_species)%vna%f(j+1)
                            r3=pseudo(the_species)%vna%d2(j)
                            r4=pseudo(the_species)%vna%d2(j+1)
                            !if(abs(a * r1 + b * r2 + c * r3 + d * r4)>1e5_double) then
                            !!   write(*,*) 'Error!'
                            !   write(*,*) j,a,b,r1,r2,r3,r4
                            !end if
                            pseudopotential(igrid) = &
                                 pseudopotential(igrid) &
                                 + a * r1 + b * r2 + c * r3 + d * r4
                         end if
                      end do
                   else if(pseudo(the_species)%tm_loc_pot==loc_pot) then
                      step = pseudo(the_species)%vlocal%delta

                      do ip=1, npoint

                         ipoint=ip_store(ip)
                         igrid=n_pts_in_block*(iblock-1)+ipoint
                         if(igrid > n_my_grid_points) call cq_abort &
                              ('set_ps: igrid error ', igrid, n_my_grid_points) 
                         r_from_i = r_store(ip)
                         x = x_store(ip)
                         y = y_store(ip)
                         z = z_store(ip)
                         j = aint( r_from_i / step ) + 1
                         !As we use the maximum of cutoff for check_block
                         ! overrun should occur in some cases.
                         !if(j > N_TAB-1) call cq_abort &
                         !     ('set_ps: overrun problem1',j)
                         ! check j (j+1 =< N_TAB)
                         ! Use the spline interpolation tables
                         !  cutoff for this function might be smaller
                         !  than cut off used in check_block
                         if(j+1 <= pseudo(the_species)%vlocal%n) then
                            ! rr is BEYOND the point we're interested in
                            rr = real(j,double) * step
                            a = ( rr - r_from_i ) / step
                            b = one - a
                            c = a * ( a * a - one ) * step * step / six
                            d = b * ( b * b - one ) * step * step / six

                            r1=pseudo(the_species)%vlocal%f(j)
                            r2=pseudo(the_species)%vlocal%f(j+1)
                            r3=pseudo(the_species)%vlocal%d2(j)
                            r4=pseudo(the_species)%vlocal%d2(j+1)

                            pseudopotential(igrid) =  pseudopotential(igrid)  + a * r1 + b * r2 + c * r3 + d * r4
                            !if(r_from_i>RD_ERR) then
                            !   pseudopotential(igrid) = pseudopotential(igrid) + &
                            !        pseudo(the_species)%zval * derf(sqrt(pseudo(the_species)%alpha)*r_from_i)/r_from_i
                            !else
                            !   pseudopotential(igrid) = pseudopotential(igrid) + &
                            !        pseudo(the_species)%zval * two*sqrt(pseudo(the_species)%alpha/pi)
                            !end if
                            r2 = r_from_i*r_from_i
                            gauss_charge = pseudo(the_species)%prefac* exp(-pseudo(the_species)%alpha * r2)
                            chlocal_density(igrid) = chlocal_density(igrid) + gauss_charge * pseudo(the_species)%zval
                         endif ! (j+1 < pseudo(the_species)%vlocal(nl)%n) then
                      enddo ! ip=1, npoint
                   else
                      step = pseudo(the_species)%chlocal%delta

                      do ip=1, npoint

                         ipoint=ip_store(ip)
                         igrid=n_pts_in_block*(iblock-1)+ipoint
                         if(igrid > n_my_grid_points) call cq_abort &
                              ('set_ps: igrid error ', igrid, n_my_grid_points) 
                         r_from_i = r_store(ip)
                         x = x_store(ip)
                         y = y_store(ip)
                         z = z_store(ip)
                         j = aint( r_from_i / step ) + 1
                         !As we use the maximum of cutoff for check_block
                         ! overrun should occur in some cases.
                         !if(j > N_TAB-1) call cq_abort &
                         !     ('set_ps: overrun problem1',j)
                         ! check j (j+1 =< N_TAB)
                         ! Use the spline interpolation tables
                         !  cutoff for this function might be smaller
                         !  than cut off used in check_block
                         if(j+1 <= pseudo(the_species)%chlocal%n) then
                            rr = real(j,double) * step
                            a = ( rr - r_from_i ) / step
                            b = one - a
                            c = a * ( a * a - one ) * step * step / six
                            d = b * ( b * b - one ) * step * step / six

                            r1=pseudo(the_species)%chlocal%f(j)
                            r2=pseudo(the_species)%chlocal%f(j+1)
                            r3=pseudo(the_species)%chlocal%d2(j)
                            r4=pseudo(the_species)%chlocal%d2(j+1)

                            core_charge =  a * r1 + b * r2 + c * r3 + d * r4
                            ! Note that (core_charge < 0)

                            chlocal_density(igrid)=chlocal_density(igrid) + core_charge

                         endif ! (j+1 < pseudo(the_species)%chlocal(nl)%n) then
                      enddo ! ip=1, npoint
                   end if ! tm_loc_pot==loc_pot
                endif! (npoint > 0) then
             enddo ! naba_atoms
          enddo ! naba_part
       endif !(naba_atoms_of_blocks(pseudo_neighbour)%no_of_part(iblock) > 0) !naba atoms?
       end if ! NAP
       ! -----------------------------------------------------------------------------------
       ! NB This loop is VERY rarely used: only when non-analytic blip operations are chosen
       ! DRB 2016/01/07
       ! -----------------------------------------------------------------------------------
       ! Bug fix starts here: loop over naba_atoms_of_blocks(dens) for NA gave array over-run on the
       ! grid for projector functions DRB 2016/06/23
       ! -----------------------------------------------------------------------------------
       !Projector Functions  ------------------------------------------
       ! I removed some fourteen year old comments (no longer relevant) DRB 2016/06/23
       if(flag_basis_set==blips.AND.(.NOT.flag_analytic_blip_int)) then
          if(naba_atoms_of_blocks(nlpf)%no_of_part(iblock) > 0) then ! if there are naba atoms
             iatom=0
             do ipart=1,naba_atoms_of_blocks(nlpf)%no_of_part(iblock)
                jpart=naba_atoms_of_blocks(nlpf)%list_part(ipart,iblock)
                if(jpart > DCS_parts%mx_gcover) then 
                   call cq_abort('set_ps: JPART ERROR ',ipart,jpart)
                endif
                ind_part=DCS_parts%lab_cell(jpart)
                do ia=1,naba_atoms_of_blocks(nlpf)%no_atom_on_part(ipart,iblock)
                   iatom=iatom+1
                   ii = naba_atoms_of_blocks(nlpf)%list_atom(iatom,iblock)
                   icover= DCS_parts%icover_ibeg(jpart)+ii-1
                   ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)

                   if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                      call cq_abort('set_ps: globID ERROR ', &
                           ii,parts%icell_beg(ind_part))
                   endif
                   if(icover > DCS_parts%mx_mcover) then
                      call cq_abort('set_ps: icover ERROR ', &
                           icover,DCS_parts%mx_mcover)
                   endif

                   !-- I (TM) am not sure where the next line should be.
                   !   (no_of_ib_ia = no_of_ib_ia +1)
                   !  This determines how to store 
                   !   pseudofunctions(ipoints,ncf,naba_atm,iblock).
                   !  Now, I assume we should consider all naba_atm whose 
                   ! distance from the block is within the maximum of core_radius.
                   ! This is needed to keep the consistency with <set_bucket>.
                   ! However, we can change this strategy by changing naba_atoms_of_blocks(pseudo_neighbour).
                   !no_of_ib_ia = no_of_ib_ia +1

                   xatom=DCS_parts%xcover(icover)
                   yatom=DCS_parts%ycover(icover)
                   zatom=DCS_parts%zcover(icover)

                   !the_species=species(ig_atom)
                   the_species = species_glob(ig_atom)

                   rcut = core_radius(the_species) + very_small   !!   2007/07/30 drb

                   call check_block &
                        (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
                        npoint,ip_store,r_store,x_store,y_store,z_store) !out
                   if(npoint>0) then
                      do nl= 1, pseudo(the_species)%n_pjnl

                         the_l=pseudo(the_species)%pjnl_l(nl)
                         offset_position = no_of_ib_ia + &
                              offset_mcomp(nl, the_species) * n_pts_in_block
                         step = pseudo(the_species)%pjnl(nl)%delta

                         do ip=1,npoint
                            ipoint=ip_store(ip)
                            position= offset_position + ipoint
                            if(position > gridfunctions(pseudofns)%size) call cq_abort &
                                 ('set_ps: position error ', position, gridfunctions(pseudofns)%size)

                            r_from_i = r_store(ip)
                            x = x_store(ip)
                            y = y_store(ip)
                            z = z_store(ip)
                            j = aint( r_from_i / step ) + 1
                            !As we use the maximum of cutoff for check_block
                            ! overrun should occur in some cases.
                            !if(j > N_TAB-1) then
                            !   call cq_abort('set_ps: overrun problem2',j)
                            !endif
                            ! check j (j+1 =< N_TAB)
                            ! Use the spline interpolation tables 
                            !  cutoff for this function might be smaller
                            !  than cut off used in check_block
                            if(j+1 <= pseudo(the_species)%pjnl(nl)%n) then
                               rr = real(j,double) * step
                               a = ( rr - r_from_i ) / step
                               b = one - a
                               c = a * ( a * a - one ) * step * step / six
                               d = b * ( b * b - one ) * step * step / six

                               r1=pseudo(the_species)%pjnl(nl)%f(j)
                               r2=pseudo(the_species)%pjnl(nl)%f(j+1)
                               r3=pseudo(the_species)%pjnl(nl)%d2(j)
                               r4=pseudo(the_species)%pjnl(nl)%d2(j+1)

                               nl_potential =  a * r1 + b * r2 + c * r3 + d * r4

                               !pjnl = chi_nl(r)/r**l ----
                               if(the_l > 0) then
                                  nl_potential= nl_potential *  r_from_i**(real(the_l,double))
                               endif
                               ! for r_from_i < RD_ERR
                               !   x, y, z are set to be 0 in check_block

                               ! Removed flag_angular_new switch DRB 2016/01/07
                               if (the_l>0) then
                                  i = 0
                                  do m = -the_l, the_l
                                     call pp_elem(nl_potential,the_l,m,x,y,z,r_from_i,val)
                                     gridfunctions(pseudofns)%griddata(position+i*n_pts_in_block) = val
                                     i = i+1
                                  enddo
                               else
                                  call pp_elem(nl_potential,the_l,0,x,y,z,r_from_i,val)
                                  gridfunctions(pseudofns)%griddata(position) = val
                               endif
                            endif !(j+1 < pseudo(the_species)%pjnl(nl)%n) then
                         enddo ! ip=1,npoint
                      enddo ! nl= 1, pseudo(the_species)%n_pjnl
                   endif! (npoint > 0) then
                   no_of_ib_ia = no_of_ib_ia + nlpf_species(the_species)*n_pts_in_block
                enddo ! naba_atoms
             enddo ! naba_part
          endif !(naba_atoms_of_blocks(pseudo_neighbour)%no_of_part(iblock) > 0) !naba atoms?
       end if ! blips and NOT analytic blip integrals
    enddo ! iblock : primary set of blocks
    ! now we must use FFT to transform the core charge density into
    ! reciprocal space to construct the pseudopotential in 
    ! reciprocal space. This is then transformed back into real space 
    ! and returned in pseudopotential from subroutine hartree()
    call my_barrier()

    ! local potential
    if(.NOT.flag_neutral_atom) then
       if(local_charge) then
          call hartree( chlocal_density, pseudopotential, maxngrid, coulomb_energy )
       else
          call hartree( chlocal_density, coulomb_potential, maxngrid, coulomb_energy )
          call axpy( n_my_grid_points, -one, coulomb_potential, 1, &
               pseudopotential, 1 )
          call start_timer(tmr_std_allocation)
          deallocate(coulomb_potential, STAT=stat)
          if(stat/=0) call cq_abort('set_ps: dealloc c_p ', stat)
          call stop_timer(tmr_std_allocation)
       end if
    end if
    
    !deallocates allocatable arrays
    call start_timer(tmr_std_allocation)
    deallocate(chlocal_density, STAT=stat)
    if(stat/=0) call cq_abort('set_ps: dealloc c_d ', stat)
    call stop_timer(tmr_std_allocation)

    call my_barrier()
    call stop_print_timer(tmr_l_tmp1,"initialising pseudopotential (TM)",IPRINT_TIME_THRES3)
    call stop_timer(tmr_std_pseudopot)
    return
  end subroutine set_tm_pseudo
!!***

! -----------------------------------------------------------
! Subroutine loc_pp_derivative_tm
! -----------------------------------------------------------

!!****f* pseudo_tm_module/loc_pp_derivative_tm *
!!
!!  NAME
!!   loc_pp_derivative_tm
!!  USAGE
!!
!!  PURPOSE
!!   Gets the Hellman-Feynman contribution to the
!!   atomic forces. For this both the electron density
!!   and Hartree potential are needed.
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   02/07/2002
!!  MODIFICATION HISTORY
!!   10:21, 2004/01/09 dave
!!    Added check on j (to stop array overrun)
!!   12:14, 17/10/2005 drb 
!!    Added deallocate for h_potential (from TM)
!!   2007/02/19 08:32 dave
!!    Added differential of gaussian charges
!!   2015/05/01 13:44 dave and sym
!!    Adding stress
!!   2015/09/04 07:53 dave
!!    Small changes to avoid messing up Hartree stress
!!   2016/01/14 16:47 dave
!!    Implementing neutral atom forces
!!   2016/01/21 08:23 dave
!!    Adding neutral atom stress and tidying
!!   2016/01/29 14:25 dave
!!    Bug fix for local G stress: accumulate loc_charge (was just storing !)
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks
!!   2016/08/08 15:30 nakata
!!    Removed unused sf in global_module
!!   2017/11/27 15:28 dave
!!    Remove calculation of HF forces for NA projectors
!!  SOURCE
!!
  subroutine loc_pp_derivative_tm ( hf_force, density, size )

    use datatypes
    use numbers
    use dimens, only: grid_point_volume, n_my_grid_points
    use global_module, only: rcellx,rcelly,rcellz,id_glob, iprint_pseudo, &
         species_glob, nlpf,ni_in_cell, flag_neutral_atom, dens
    use block_module, only : n_pts_in_block
    use group_module, only : blocks, parts
    use primary_module, only: domain
    use cover_module, only: DCS_parts
    use set_blipgrid_module, only : naba_atoms_of_blocks

    use species_module, only: species
    use GenComms, only: gsum, cq_abort, inode, ionode
    use hartree_module, only: hartree, Hartree_stress
    use maxima_module, only: maxngrid
    use density_module, only: density_atom
    use atomic_density, only: atomic_density_table ! for Neutral atom potential
    use pseudopotential_common, only: flag_neutral_atom_projector

    implicit none   

    ! Passed variables
    integer :: size
    real(double),intent(in)  :: density( size )
    real(double),intent(out) :: hf_force( 3, ni_in_cell )

    ! Local variables
    integer :: i, j, pseudo_neighbour

    real(double) :: a, b, c, d, alpha, beta, gamma, delta, derivative 
    real(double) :: fx_1, fy_1, fz_1
    real(double) :: fx_2, fy_2, fz_2
    real(double) :: h_energy, r_from_i, x, y, z, step, elec_here, front, gauss
    real(double):: r1, r2, r3, r4, da, db, dc, dd, gauss_charge

    real(double):: dcellx_block,dcelly_block,dcellz_block
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom
    real(double):: xatom,yatom,zatom
    real(double):: xblock,yblock,zblock
    integer :: the_species
    integer :: iblock,ipoint,igrid
    real(double) :: rr
    integer ::  iatom, stat, ip, npoint, nl
    real(double) :: rcut, temp(3)

    ! allocatable
    real(double),allocatable, dimension(:) :: h_potential, drho_tot
    ! Store potential scaled by 2GaGb/G^4 for stress
    real(double),allocatable :: loc_charge(:)

    call start_timer(tmr_std_pseudopot)
    if(iprint_pseudo>2.AND.inode==ionode) write(io_lun,fmt='(4x,"Doing TM force with pseudotype: ",i3)') pseudo(1)%tm_loc_pot
    ! the structure of this subroutine is similar to set_tm_pseudo et.
    HF_force = zero
    loc_HF_stress = zero
    loc_G_stress = zero

    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz

    ! get Hartree potential
    call start_timer(tmr_std_allocation)
    if(flag_neutral_atom) then
       ! NB we will use the variable h_potential even though it is the drho potential
       allocate(h_potential(maxngrid),STAT=stat)
       if(stat /= 0) call cq_abort &
            ('loc_pp_derivative_tm: alloc h_potential',stat)
       allocate(drho_tot(maxngrid),STAT=stat)
       if(stat /= 0) call cq_abort &
            ('loc_pp_derivative_tm: alloc drho_tot',stat)
       drho_tot = zero
       drho_tot(1:size) = density(1:size) - density_atom(1:size)
    else
       allocate(h_potential(maxngrid),STAT=stat)
       if(stat /= 0) call cq_abort &
            ('loc_pp_derivative_tm: alloc h_potential',stat)
       allocate(loc_charge(maxngrid),STAT=stat)
       if(stat /= 0) call cq_abort &
            ('loc_pp_derivative_tm: alloc loc_charge',stat)
    end if
    call stop_timer(tmr_std_allocation)
    ! Save existing Hartree stress
    temp = Hartree_stress
    if(flag_neutral_atom) then
       call hartree( drho_tot, h_potential, maxngrid, h_energy )
       pseudo_neighbour = dens
    else
       call hartree( density, h_potential, maxngrid, h_energy )
       pseudo_neighbour = nlpf
       loc_charge = zero
    end if
    Hartree_stress = temp

    ! now loop over grid points and accumulate HF part of the force
    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
       if(naba_atoms_of_blocks(pseudo_neighbour)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom=0
          do ipart=1,naba_atoms_of_blocks(pseudo_neighbour)%no_of_part(iblock)
             jpart=naba_atoms_of_blocks(pseudo_neighbour)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then
                call cq_abort('set_ps: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atoms_of_blocks(pseudo_neighbour)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atoms_of_blocks(pseudo_neighbour)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)

                if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                   call cq_abort('set_ps: globID ERROR ', &
                        ii,parts%icell_beg(ind_part))
                endif
                if(ig_atom > ni_in_cell) then
                   call cq_abort('HF_for: igatom ERROR ', &
                        ii,ig_atom)
                endif
                if(icover > DCS_parts%mx_mcover) then
                   call cq_abort('set_ps: icover ERROR ', &
                        icover,DCS_parts%mx_mcover)
                endif

                xatom=DCS_parts%xcover(icover)
                yatom=DCS_parts%ycover(icover)
                zatom=DCS_parts%zcover(icover)

                !the_species=species(ig_atom)
                the_species = species_glob(ig_atom)

                !calculates distances between the atom and integration grid points
                !in the block and stores which integration grids are neighbours.
                if( flag_neutral_atom ) then ! for Neutral atom potential
                   rcut = pseudo(the_species)%vna%cutoff + very_small   !!   30072007 drb
                else ! ABINIT or SIESTA PP
                   rcut = core_radius(the_species) + very_small   !!   03032003TM
                end if
                call check_block &
                     (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
                     npoint,ip_store,r_store,x_store,y_store,z_store) !out

                !  When we use SIESTA's pseudopotential, the local part of 
                ! pseudopotentials is generated by the pseudo core charge
                ! (local charge) distribution rho_chlocal(r).
                !  Thus,
                !   HF_force = - d/dR_{mu} 
                !       ( \int d3r \int d3r' rho_chlocal(r) rho(r')/ |r-r'| )
                !  Here, 
                !   rho_chlocal(r) = \sum_{mu} rho^{mu}_chlocal(|r-R_mu|)
                !  and
                !   \int d3r' rho(r')/|r-r'| is calculated by using FFT 
                !  (by calling hartree) and is stored as h_potential(:)
                ! 
                !  - d/dR_mu (rho_chlocal (r) ) = 
                !        d(rho^{mu}_chlocal(s))/ds |(s=|r-R_mu|)
                !         * (r-R_mu)/|r-R_mu|
                !  the derivative in l.h.s will be calculated by using
                !  spline interpolations, as in the following.
                !        02/07/2002  T. Miyazaki
                !
                if(npoint > 0) then
                   do ip=1, npoint

                      ipoint=ip_store(ip)
                      igrid=n_pts_in_block*(iblock-1)+ipoint
                      if(igrid > n_my_grid_points) call cq_abort &
                           ('set_ps: igrid error ', igrid, n_my_grid_points)
                      r_from_i = r_store(ip)
                      x = x_store(ip)
                      y = y_store(ip)
                      z = z_store(ip)
                      if(flag_neutral_atom) then
                         step = pseudo(the_species)%vna%delta
                      else if(pseudo(the_species)%tm_loc_pot==loc_pot) then
                         step = pseudo(the_species)%vlocal%delta
                      else
                         step = pseudo(the_species)%chlocal%delta
                      end if
                      j = aint( r_from_i / step ) + 1
                      !As we use the maximum of cutoff for check_block
                      ! overrun should occur in some cases.
                      !if(j > N_TAB-1) call cq_abort &
                      !     ('set_ps: overrun problem3',j)
                      ! check j (j+1 =< N_TAB)
                      ! Use the spline interpolation tables
                      if(flag_neutral_atom) then
                         if(j+1 <= pseudo(the_species)%vna%n.AND.(.NOT.flag_neutral_atom_projector)) then
                            rr = real(j,double) * step
                            a = ( rr - r_from_i ) / step
                            b = one - a

                            alpha = -one / step
                            beta =  one / step
                            gamma = -step * ( three * a * a - one ) / six
                            delta =  step * ( three * b * b - one ) / six

                            r1=pseudo(the_species)%vna%f(j)
                            r2=pseudo(the_species)%vna%f(j+1)
                            r3=pseudo(the_species)%vna%d2(j)
                            r4=pseudo(the_species)%vna%d2(j+1)

                            derivative = &
                                 alpha * r1 + beta * r2 + gamma * r3 + delta * r4

                            fx_2 = x * derivative * density( igrid )
                            fy_2 = y * derivative * density( igrid )
                            fz_2 = z * derivative * density( igrid )
                            i=ig_atom
                            HF_force(1,i) = HF_force(1,i) + fx_2 * grid_point_volume
                            HF_force(2,i) = HF_force(2,i) + fy_2 * grid_point_volume
                            HF_force(3,i) = HF_force(3,i) + fz_2 * grid_point_volume
                            loc_HF_stress(1) = loc_HF_stress(1) + fx_2 * grid_point_volume * x*r_from_i
                            loc_HF_stress(2) = loc_HF_stress(2) + fy_2 * grid_point_volume * y*r_from_i
                            loc_HF_stress(3) = loc_HF_stress(3) + fz_2 * grid_point_volume * z*r_from_i
                         end if ! j+1<pseudo(the_species)%vna%n
                         ! Now the atomic density interacting with the potential from drho
                         ! Note that this is the force from the SCREENED Hartree potential
                         step = atomic_density_table(the_species)%delta
                         j = aint( r_from_i / step ) + 1

                         if(j+1 <= atomic_density_table(the_species)%length) then
                            rr = real(j,double) * step
                            a = ( rr - r_from_i ) / step
                            b = one - a

                            alpha = -one / step
                            beta =  one / step
                            gamma = -step * ( three * a * a - one ) / six
                            delta =  step * ( three * b * b - one ) / six

                            r1=atomic_density_table(the_species)%table(j)
                            r2=atomic_density_table(the_species)%table(j+1)
                            r3=atomic_density_table(the_species)%d2_table(j)
                            r4=atomic_density_table(the_species)%d2_table(j+1)

                            derivative =  &
                                 alpha * r1 + beta * r2 + gamma * r3 + delta * r4

                            fx_2 = minus_one * x * derivative * h_potential( igrid )
                            fy_2 = minus_one * y * derivative * h_potential( igrid )
                            fz_2 = minus_one * z * derivative * h_potential( igrid )
                            i=ig_atom
                            HF_force(1,i) = HF_force(1,i) + fx_2 * grid_point_volume
                            HF_force(2,i) = HF_force(2,i) + fy_2 * grid_point_volume
                            HF_force(3,i) = HF_force(3,i) + fz_2 * grid_point_volume
                            loc_HF_stress(1) = loc_HF_stress(1) + fx_2 * grid_point_volume * x*r_from_i
                            loc_HF_stress(2) = loc_HF_stress(2) + fy_2 * grid_point_volume * y*r_from_i
                            loc_HF_stress(3) = loc_HF_stress(3) + fz_2 * grid_point_volume * z*r_from_i

                         end if ! j+1<atomic_density_table(the_species)%length

                      else if(pseudo(the_species)%tm_loc_pot==loc_pot) then
                         if(j+1 <= pseudo(the_species)%vlocal%n) then
                            elec_here = density(igrid) * grid_point_volume
                            gauss = exp( -pseudo(the_species)%alpha * r_from_i*r_from_i )
                            rr = real(j,double) * step
                            a = ( rr - r_from_i ) / step
                            b = one - a
                            da = -one / step
                            db =  one / step
                            dc = -step * ( three * a * a - one ) / six
                            dd =  step * ( three * b * b - one ) / six
                            derivative =  &
                                 da * pseudo(the_species)%vlocal%f(j) +  &
                                 db * pseudo(the_species)%vlocal%f(j+1) + &
                                 dc * pseudo(the_species)%vlocal%d2(j) +  &
                                 dd * pseudo(the_species)%vlocal%d2(j+1)

!                            write(20+inode,fmt='(a15,i4,4f12.4)') 'Atom, x, y, z: ',ig_atom,x,y,z,derivative
                            ! This is a little obscure, but the multiples below have had 
                            ! the minus sign which should be there removed to correct an 
                            ! earlier on dropped in the da, db, dc, dd expressions.DRB 27.11.97
                            fx_1 = x * derivative
                            fy_1 = y * derivative 
                            fz_1 = z * derivative 
                            front = two * pseudo(the_species)%alpha * pseudo(the_species)%prefac * &
                                 gauss * h_potential( igrid ) * pseudo(the_species)%zval* r_from_i * grid_point_volume
                            fx_2 = front* x 
                            fy_2 = front* y
                            fz_2 = front* z
                            i=ig_atom
                            HF_force(1,ig_atom) = HF_force(1,ig_atom) + fx_1 * elec_here + fx_2
                            HF_force(2,ig_atom) = HF_force(2,ig_atom) + fy_1 * elec_here + fy_2
                            HF_force(3,ig_atom) = HF_force(3,ig_atom) + fz_1 * elec_here + fz_2
                            loc_HF_stress(1) = loc_HF_stress(1) + (fx_1 * elec_here + fx_2) * x*r_from_i
                            loc_HF_stress(2) = loc_HF_stress(2) + (fy_1 * elec_here + fy_2) * y*r_from_i
                            loc_HF_stress(3) = loc_HF_stress(3) + (fz_1 * elec_here + fz_2) * z*r_from_i
                            r2 = r_from_i*r_from_i
                            gauss_charge = pseudo(the_species)%prefac* exp(-pseudo(the_species)%alpha * r2)
                            loc_charge(igrid) = loc_charge(igrid) + gauss_charge * pseudo(the_species)%zval
                         end if ! j+1<pseudo(the_species)%chlocal%n
                      else
                         if(j+1 <= pseudo(the_species)%chlocal%n) then
                            rr = real(j,double) * step
                            a = ( rr - r_from_i ) / step
                            b = one - a
                            c = a * ( a * a - one ) * step * step / six
                            d = b * ( b * b - one ) * step * step / six

                            alpha = -one / step
                            beta =  one / step
                            gamma = -step * ( three * a * a - one ) / six
                            delta =  step * ( three * b * b - one ) / six

                            r1=pseudo(the_species)%chlocal%f(j)
                            r2=pseudo(the_species)%chlocal%f(j+1)
                            r3=pseudo(the_species)%chlocal%d2(j)
                            r4=pseudo(the_species)%chlocal%d2(j+1)

                            derivative =  &
                                 alpha * r1 + beta * r2 + gamma * r3 + delta * r4

                            !elec_here = density(igrid) * grid_point_volume

                            fx_2 = x * derivative * h_potential( igrid )
                            fy_2 = y * derivative * h_potential( igrid )
                            fz_2 = z * derivative * h_potential( igrid )


                            i=ig_atom
                            ! *grid_point_volume is added  03032003 T. MIYAZAKI
                            HF_force(1,i) = HF_force(1,i) + fx_2 * grid_point_volume
                            HF_force(2,i) = HF_force(2,i) + fy_2 * grid_point_volume
                            HF_force(3,i) = HF_force(3,i) + fz_2 * grid_point_volume

                            ! SYM, DRB 2014/09/02 and 2015/03: Accumulate HF stress and local charge (for reciprocal space stress)
                            loc_HF_stress(1) = loc_HF_stress(1) + fx_2 * grid_point_volume * x*r_from_i
                            loc_HF_stress(2) = loc_HF_stress(2) + fy_2 * grid_point_volume * y*r_from_i
                            loc_HF_stress(3) = loc_HF_stress(3) + fz_2 * grid_point_volume * z*r_from_i
                            loc_charge(igrid) = loc_charge(igrid)  + a*r1+b*r2+c*r3+d*r4
                         end if ! j+1<pseudo(the_species)%chlocal%n
                      end if
                   enddo ! ip=1, npoint
                endif  !(npoint > 0) then

             enddo ! ia=1,naba_atoms_of_blocks(pseudo_neighbour)%no_atom_on_part(ipart,iblock)
          enddo ! ipart=1,naba_atoms_of_blocks(pseudo_neighbour)%no_of_part(iblock)
       endif    ! (naba_atoms_of_blocks(pseudo_neighbour)%no_of_part(iblock) > 0) then ! if there are naba atoms
    enddo ! iblock = 1, domain%groups_on_node ! primary set of blocks
    ! Deallocate added by TM, 2005/08/11
    if(flag_neutral_atom) then
       call start_timer(tmr_std_allocation)
       deallocate(drho_tot,STAT=stat)
       if(stat /= 0) call cq_abort &
            ('loc_pp_derivative_tm: dealloc drho_tot',stat)
    else
       ! Save existing Hartree stress
       temp = Hartree_stress
       call hartree( density, h_potential, maxngrid, h_energy, loc_charge,loc_G_stress )
       Hartree_stress = temp
       if(pseudo(the_species)%tm_loc_pot==loc_pot) loc_G_stress = -loc_G_stress
       call start_timer(tmr_std_allocation)
       deallocate(loc_charge,STAT=stat)
       if(stat /= 0) call cq_abort &
            ('loc_pp_derivative_tm: dealloc loc_charge',stat)
    end if
    deallocate(h_potential,STAT=stat)
    if(stat /= 0) call cq_abort &
         ('loc_pp_derivative_tm: dealloc h_potential',stat)
    call stop_timer(tmr_std_allocation)
    ! and add contributions from all nodes
    !  In the future this should be replaced by summation with local communication
    !    Tsuyoshi Miyazaki
    call gsum(HF_force,3,ni_in_cell)
    call gsum(loc_HF_stress,3)
    ! Don't gsum loc_G_stress - that's done in hartree
    call stop_timer(tmr_std_pseudopot)

    return
  end subroutine loc_pp_derivative_tm
!!***

! -----------------------------------------------------------
! Subroutine nonloc_pp_derivative_tm
! -----------------------------------------------------------

!!****f* pseudo_tm_module/nonloc_pp_derivative_tm *
!!
!!  NAME 
!!   nonloc_pp_derivative_tm
!!  USAGE
!! 
!!  PURPOSE
!!   Calculates the derivatives of projector functions 
!!   w.r.t. positions of the atoms.
!!     d(chi^mu_ln(r))/d(R_mu,direction)
!!  INPUTS
!!   direction: 1=x, 2=y, 3=z
!!  OUTPUTS
!!   d_pseudofunctions =  d(chi^mu_ln(r))/d(R_mu,direction)
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   01/07/2002
!!  MODIFICATION HISTORY
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks
!!  SOURCE
!!
  subroutine nonloc_pp_derivative_tm(direction, dpseudofns)

    use datatypes
    use GenComms, only: inode,ionode
    use numbers
    use global_module, only: rcellx,rcelly,rcellz,id_glob,ni_in_cell,iprint_pseudo, species_glob, nlpf
    use species_module, only: species, nlpf_species
    !  At present, these arrays are dummy arguments.
    use block_module, only : nx_in_block,ny_in_block,nz_in_block, &
         n_pts_in_block
    use group_module, only : blocks, parts
    use primary_module, only: domain
    use cover_module, only: DCS_parts
    use set_blipgrid_module, only : naba_atoms_of_blocks
    use GenComms, only: my_barrier, cq_abort
    use angular_coeff_routines, only : pp_gradient
    use functions_on_grid, only: gridfunctions, fn_on_grid

    implicit none 

    !Passed
    integer, intent(in) :: direction
    integer :: dpseudofns

    ! local
    real(double):: dcellx_block,dcelly_block,dcellz_block
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom
    real(double):: xatom,yatom,zatom,step
    real(double):: xblock,yblock,zblock
    integer :: the_species
    integer :: j,iblock,the_l,ipoint, ll
    real(double) :: r_from_i, r_the_l
    real(double) :: rr,a,b,c,d,x,y,z,nl_potential,nl_potential_derivative
    real(double) :: alpha,beta,gamma,delta
    integer :: no_of_ib_ia, offset_position
    integer :: position,iatom
    integer :: nl, npoint, ip
    real(double) :: r1, r2, r3, r4
    real(double) :: rcut, rl, rl1
    real(double) :: nl_potential_derivative_new, val, nl_potential_new
    integer :: m, i

    ! --  Start of subroutine  ---
    !  This subroutine calculates the derivatives of 
    ! projector functions w.r.t. atomic positions.
    ! The structure of this subroutine is almost the same as in
    ! set_pseudopotential_tm
    
    call start_timer(tmr_std_pseudopot)
    gridfunctions(dpseudofns)%griddata = zero

    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz

    no_of_ib_ia = 0

    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
       if(naba_atoms_of_blocks(nlpf)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom=0
          do ipart=1,naba_atoms_of_blocks(nlpf)%no_of_part(iblock)
             jpart=naba_atoms_of_blocks(nlpf)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('ps_derivative: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atoms_of_blocks(nlpf)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atoms_of_blocks(nlpf)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)

                if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                   call cq_abort('ps_derivative: globID ERROR ', &
                        ii,parts%icell_beg(ind_part))
                endif
                if(icover > DCS_parts%mx_mcover) then
                   call cq_abort('ps_derivative: icover ERROR ', &
                        icover,DCS_parts%mx_mcover)
                endif

                !-- I (TM) am not sure where the next line should be.
                !   See the comment in set_pseudopotential_tm
                !no_of_ib_ia = no_of_ib_ia +1

                xatom=DCS_parts%xcover(icover)
                yatom=DCS_parts%ycover(icover)
                zatom=DCS_parts%zcover(icover)

                !the_species=species(ig_atom)
                the_species = species_glob(ig_atom)

                !calculates distances between the atom and integration grid points
                !in the block and stores which integration grids are neighbours.
                rcut = core_radius(the_species) + very_small   !!   03032003TM
                call check_block &
                     (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &   ! in
                     npoint,ip_store,r_store,x_store,y_store,z_store) !out

                !Projector Functions  ------------------------------------------
                ! 1D array pseudofunctions is like
                ! pseudofunction(n_pts_in_block,ncf,naba_atm,iblock) 

                if(npoint > 0) then
                   do nl= 1, pseudo(the_species)%n_pjnl

                      the_l=pseudo(the_species)%pjnl_l(nl)
                      offset_position = no_of_ib_ia + &
                           offset_mcomp(nl, the_species) * n_pts_in_block
                      step = pseudo(the_species)%pjnl(nl)%delta

                      do ip=1,npoint
                         ipoint=ip_store(ip)
                         position= offset_position + ipoint
                         if(position > gridfunctions(dpseudofns)%size) call cq_abort &
                              ('ps_derivative: position error ', position, gridfunctions(dpseudofns)%size)
                         
                         r_from_i = r_store(ip)
                         x = x_store(ip)  
                         y = y_store(ip)
                         z = z_store(ip)
                         j = aint( r_from_i / step ) + 1
                         !As we use the maximum of cutoff for check_block
                         ! overrun should occur in some cases.
                         !if(j > N_TAB) then
                         !   call cq_abort('ps_derivative: overrun problem4',j)
                         ! check j (j+1 =< N_TAB)
                         !endif
                         ! Use the spline interpolation tables 
                         !  cutoff for this function might be smaller
                         !  than cut off used in check_block
                         ! FIX 2007/03/19 DRB
                         if(j+1 <= pseudo(the_species)%pjnl(nl)%n) then
                            !if(j+1 < pseudo(the_species)%chlocal%n) then
                            rr = real(j,double) * step
                            a = ( rr - r_from_i ) / step
                            b = one - a
                            c = a * ( a * a - one ) * step * step / six
                            d = b * ( b * b - one ) * step * step / six

                            alpha = -one / step
                            beta =  one / step
                            gamma = -step * ( three * a * a - one ) / six
                            delta =  step * ( three * b * b - one ) / six

                            r1=pseudo(the_species)%pjnl(nl)%f(j)
                            r2=pseudo(the_species)%pjnl(nl)%f(j+1)
                            r3=pseudo(the_species)%pjnl(nl)%d2(j)
                            r4=pseudo(the_species)%pjnl(nl)%d2(j+1)

                            nl_potential =  a * r1 + b * r2 + c * r3 + d * r4
                            nl_potential_derivative =  &
                                 alpha * r1 + beta * r2 + gamma * r3 + delta * r4

                            ! Note that  pjnl = chi_ln(r)/ r**l
                            !   d(chi_ln(r)/r**l )/dr = 
                            !      -l r**(-l-1) * chi_ln + r**(-l) * d(chi_ln)/dr
                            !   d(chi_ln(r))/dr = r**l d(chi_ln(r)/r**l)/dr + 
                            !                      l*r**(l-1) (chi_ln(r)/r**l)

                            if(r_from_i > RD_ERR) then
                               !OLD if(the_l == 1) then
                               !OLD    nl_potential_derivative = r_from_i * nl_potential_derivative &
                               !OLD         + nl_potential
                               !OLD    nl_potential =  nl_potential/r_from_i
                               !OLD elseif(the_l > 1) then
                               !OLDif(the_l > 0) then
                               !OLD   rl1= r_from_i**(the_l-1)
                               !OLD   rl = rl1 * r_from_i
                               !OLD   nl_potential_derivative = rl * nl_potential_derivative &
                               !OLD        + the_l * rl1 * nl_potential
                               !OLD   nl_potential =  nl_potential/rl
                               !OLDendif
                            else
                               nl_potential_derivative = zero
                               nl_potential= zero   !!   03032003TM
                               r_from_i = RD_ERR
                            endif
                            
                            !RC putting in call to new spherical harmonic scheme here;
                            if(flag_angular_new) then !using new spherical harmonics
                               !forming nl_potential derivative
                               if(the_l>0) then
                                  r_the_l = real(the_l,double)
                                  ! More efficient way of doing r**l and r**(l-1)
                                  rl = r_from_i
                                  rl1 = one
                                  do ll=1,the_l-1
                                     rl = rl*r_from_i
                                     rl1 = rl1*r_from_i
                                  end do
                                  nl_potential_derivative_new = rl*nl_potential_derivative&
                                       &+r_the_l*rl1*nl_potential
                                  nl_potential_new = nl_potential*rl
                               else
                                  nl_potential_derivative_new = nl_potential_derivative
                                  nl_potential_new = nl_potential
                               end if
                               if(the_l.gt.0) then
                                  i = 0
                                  do m = -the_l, the_l
                                     !i = 0
                                     !call pp_elem_derivative(direction,nl_potential_new,& 
                                     ! Timing in area 11 (basis) is necessary here, 
                                     !   to avoid double start of timer
                                     call start_timer(tmr_std_basis)
                                     call pp_gradient(direction,nl_potential_new,&
                                          &nl_potential_derivative_new,x,y,z,r_from_i,the_l,m,val)
                                     call stop_timer(tmr_std_basis)
                                     gridfunctions(dpseudofns)%griddata(position+i*n_pts_in_block) = val
                                     i = i+1
                                     !RC for s/p value comparisons
                                     !if(the_l.lt.2) then
                                     !write(45+inode,*) direction,m,the_l, val, 'val'
                                     !else
                                     !   continue
                                     !endif
                                  enddo
                               else
                                  !call pp_elem_derivative(direction,nl_potential_new,&
                                  call pp_gradient(direction,nl_potential_new,&
                                       &nl_potential_derivative_new,x,y,z,r_from_i,the_l,0,val)
                                  gridfunctions(dpseudofns)%griddata(position) = val
                                  !!  write(45+inode,*)  0,the_l, val, 'val'
                               endif
                            else
                               call cq_abort("Must have FlagNewAngular T for siesta/abinit pseudos")
                               !%%! ! S, P or D components
                               !%%! select case(the_l)
                               !%%! case (0)     ! S component
                               !%%!    if(direction == 1) then
                               !%%!       gridfunctions(dpseudofns)%griddata(position) =   &
                               !%%!            nl_potential_derivative * spherical_harmonic_norm(1) * x
                               !%%!    elseif(direction == 2) then
                               !%%!       gridfunctions(dpseudofns)%griddata(position) =   &
                               !%%!            nl_potential_derivative * spherical_harmonic_norm(1) * y
                               !%%!    else
                               !%%!       gridfunctions(dpseudofns)%griddata(position) =   &
                               !%%!            nl_potential_derivative * spherical_harmonic_norm(1) * z
                               !%%!    endif
                               !%%! case (1)     !  P components
                               !%%!    if(direction == 1) then
                               !%%!       !Px
                               !%%!       gridfunctions(dpseudofns)%griddata(position) =   &
                               !%%!            spherical_harmonic_norm(2) *   &
                               !%%!            ( nl_potential_derivative * x * x * r_from_i &
                               !%%!            + nl_potential)
                               !%%!       !Py
                               !%%!       gridfunctions(dpseudofns)%griddata(position+1*n_pts_in_block) =   &
                               !%%!            spherical_harmonic_norm(3) *   &
                               !%%!            nl_potential_derivative * y * x * r_from_i
                               !%%!       !Pz
                               !%%!       gridfunctions(dpseudofns)%griddata(position+2*n_pts_in_block) =   &
                               !%%!            spherical_harmonic_norm(4) *   &
                               !%%!            nl_potential_derivative * z * x * r_from_i
                               !%%!    elseif(direction == 2) then
                               !%%!       !Px
                               !%%!       gridfunctions(dpseudofns)%griddata(position) =   &
                               !%%!            spherical_harmonic_norm(2) *   &
                               !%%!            nl_potential_derivative * x * y * r_from_i
                               !%%!       !Py
                               !%%!       gridfunctions(dpseudofns)%griddata(position+1*n_pts_in_block) =   &
                               !%%!            spherical_harmonic_norm(3) *   &
                               !%%!            ( nl_potential_derivative * y * y * r_from_i &
                               !%%!            + nl_potential)
                               !%%!       !Pz
                               !%%!       gridfunctions(dpseudofns)%griddata(position+2*n_pts_in_block) =   &
                               !%%!            spherical_harmonic_norm(4) *   &
                               !%%!            nl_potential_derivative * z * y * r_from_i
                               !%%!    else
                               !%%!       !Px
                               !%%!       gridfunctions(dpseudofns)%griddata(position) =   &
                               !%%!         spherical_harmonic_norm(2) *   &
                               !%%!         nl_potential_derivative * x * z * r_from_i
                               !%%!       !Py
                               !%%!       gridfunctions(dpseudofns)%griddata(position+1*n_pts_in_block) =   &
                               !%%!            spherical_harmonic_norm(3) *   &
                               !%%!            nl_potential_derivative * y * z * r_from_i
                               !%%!       !Pz
                               !%%!       gridfunctions(dpseudofns)%griddata(position+2*n_pts_in_block) =   &
                               !%%!            spherical_harmonic_norm(4) *   &
                               !%%!            ( nl_potential_derivative * z * z * r_from_i &
                               !%%!            + nl_potential)
                               !%%!    endif
                               !%%! case (2)     !  D components
                               !%%!    if(direction == 1) then
                               !%%!       !D(x^2-y^2)
                               !%%!       gridfunctions(dpseudofns)%griddata(position) =                               &
                               !%%!            spherical_harmonic_norm(5) * r_from_i *             &
                               !%%!            ( nl_potential_derivative * r_from_i * (x*x-y*y) * x &
                               !%%!            + nl_potential * two * x )
                               !%%!       !D(3z^2-1)
                               !%%!       gridfunctions(dpseudofns)%griddata(position+1*n_pts_in_block) =              &
                               !%%!            spherical_harmonic_norm(6) * r_from_i *             &
                               !%%!            ( nl_potential_derivative * r_from_i * (three*z*z-one)*x &
                               !%%!            + nl_potential * (-two) * x )
                               !%%!       !D(xy)
                               !%%!       gridfunctions(dpseudofns)%griddata(position+2*n_pts_in_block) =              &
                               !%%!            spherical_harmonic_norm(7) * r_from_i *             &
                               !%%!            ( nl_potential_derivative * r_from_i * x * y * x  &
                               !%%!            + nl_potential * y )
                               !%%!       !D(xz)
                               !%%!       gridfunctions(dpseudofns)%griddata(position+3*n_pts_in_block) =              &
                               !%%!            spherical_harmonic_norm(8) * r_from_i *             &
                               !%%!            ( nl_potential_derivative * r_from_i * x * z * x  &
                               !%%!            + nl_potential * z )
                               !%%!       !D(yz)
                               !%%!       gridfunctions(dpseudofns)%griddata(position+4*n_pts_in_block) =              &
                               !%%!            spherical_harmonic_norm(9) * r_from_i *             &
                               !%%!            ( nl_potential_derivative * r_from_i * y * z * x )
                               !%%!    elseif(direction == 2) then
                               !%%!       !D(x^2-y^2)
                               !%%!       gridfunctions(dpseudofns)%griddata(position) =                             &
                               !%%!            spherical_harmonic_norm(5) * r_from_i *              &
                               !%%!            ( nl_potential_derivative * r_from_i * (x*x-y*y) * y &
                               !%%!            + nl_potential * (-two) * y )
                               !%%!       !D(3z^2-1)
                               !%%!       gridfunctions(dpseudofns)%griddata(position+1*n_pts_in_block) =            &
                               !%%!            spherical_harmonic_norm(6) * r_from_i *              &
                               !%%!            ( nl_potential_derivative * r_from_i * (three*z*z-one)*y &
                               !%%!            + nl_potential * (-two) * y )
                               !%%!       !D(xy)
                               !%%!       gridfunctions(dpseudofns)%griddata(position+2*n_pts_in_block) =              &
                               !%%!            spherical_harmonic_norm(7) * r_from_i *             &
                               !%%!            ( nl_potential_derivative * r_from_i * x * y * y  &
                               !%%!            + nl_potential * x )
                               !%%!       !D(xz)
                               !%%!       gridfunctions(dpseudofns)%griddata(position+3*n_pts_in_block) =              &
                               !%%!            spherical_harmonic_norm(8) * r_from_i *             &
                               !%%!         ( nl_potential_derivative * r_from_i * x * z * y )
                               !%%!       !D(yz)
                               !%%!       gridfunctions(dpseudofns)%griddata(position+4*n_pts_in_block) =              &
                               !%%!            spherical_harmonic_norm(9) * r_from_i *             &
                               !%%!            ( nl_potential_derivative * r_from_i * y * z * y  &
                               !%%!            + nl_potential * z )
                               !%%!    else
                               !%%!       !D(x^2-y^2)
                               !%%!       gridfunctions(dpseudofns)%griddata(position) =                              &
                               !%%!            spherical_harmonic_norm(5) *                            &
                               !%%!            nl_potential_derivative * (x*x-y*y)*z  &
                               !%%!            * r_from_i **2
                               !%%!       !D(3z^2-1)
                               !%%!       gridfunctions(dpseudofns)%griddata(position+1*n_pts_in_block) =             &
                               !%%!            spherical_harmonic_norm(6) * r_from_i *                 &
                               !%%!            ( nl_potential_derivative * r_from_i * (three*z*z-one)*z &
                               !%%!            + nl_potential * four * z ) 
                               !%%!       !D(xy)
                               !%%!       gridfunctions(dpseudofns)%griddata(position+2*n_pts_in_block) =             &
                               !%%!            spherical_harmonic_norm(7) * r_from_i *               &
                               !%%!            ( nl_potential_derivative *r_from_i * x * y * z )
                               !%%!       !D(xz)
                               !%%!       gridfunctions(dpseudofns)%griddata(position+3*n_pts_in_block) =             &
                               !%%!            spherical_harmonic_norm(8) * r_from_i * &
                               !%%!            ( nl_potential_derivative * r_from_i* x * z * z     &
                               !%%!            + nl_potential * x ) 
                               !%%!       !D(yz)
                               !%%!       gridfunctions(dpseudofns)%griddata(position+4*n_pts_in_block) =             &
                               !%%!            spherical_harmonic_norm(9) * r_from_i * &
                               !%%!            ( nl_potential_derivative * r_from_i* y * z * z     &
                               !%%!            + nl_potential * y ) 
                               !%%!    endif
                               !%%! end select
                            endif !flag_angular_new
                         endif ! (j+1 < pseudo(the_species)%chlocal%n) then
                      enddo ! ip=1,npoint
                   enddo ! nl= 1, pseudo(the_species)%n_pjnl
                endif !(npoint > 0) then
                no_of_ib_ia = no_of_ib_ia + nlpf_species(the_species)*n_pts_in_block
             enddo ! naba_atoms
          enddo ! naba_part
       endif !(naba_atoms_of_blocks(nlpf)%no_of_part(iblock) > 0) !naba atoms?
    enddo ! iblock : primary set of blocks
    call stop_timer(tmr_std_pseudopot)
    
    return
  end subroutine nonloc_pp_derivative_tm
!!***

! -----------------------------------------------------------
! Subroutine get_energy_shift
! -----------------------------------------------------------

!!****f* pseudo_tm_module/get_energy_shift *
!!
!!  NAME
!!   get_energy_shift
!!  USAGE
!!
!!  PURPOSE
!!   calculates the contribution from G=0 term, as
!!  INPUTS
!!
!!  OUTPUTS
!!   real(double) :: e_core - energy calculated
!!  USES
!!   datatypes, numbers, common, dimens, atoms, species_module
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   08/07/2002
!!  MODIFICATION HISTORY
!!   2008/03/03 18:51 dave
!!    Changed float to real
!!   2016/02/09 08:17 dave
!!    Changed to use erfc from functions module
!!  SOURCE
!!
  subroutine get_energy_shift(e_core)

    use datatypes
    use numbers
    use GenComms, only: inode,ionode
    use dimens, only: volume
    use species_module, only: species, n_species, type_species
    use global_module, only: ni_in_cell, iprint_pseudo
    use GenComms, only: cq_abort

    implicit none

    ! Passed variables
    real(double), intent(OUT) :: e_core

    ! Local variables
    integer :: iatom, ispecies, stat
    integer :: number_of_electrons
    integer, allocatable :: n_atoms_of_species(:)
    real(double), allocatable :: eshift(:)
    real(double), parameter :: eps = 1.0e-5_double
    real(double) :: eshift1, eshift2, beta, z, eshift3, vlocG0_1, vlocG0_2, vlocG0_3
    integer :: nfac1, nfac2, nfac3
    integer :: n_atoms
     n_atoms = ni_in_cell

    call start_timer(tmr_std_allocation)
    allocate(n_atoms_of_species(n_species), STAT=stat)
    if(stat /= 0) call cq_abort &
         ('Error in allocating n_atoms_of_species',stat, n_species)
    allocate(eshift(n_species), STAT=stat)
    if(stat /= 0) call cq_abort &
         ('Error in allocating eshift',stat, n_species)
    call stop_timer(tmr_std_allocation)

    ! first calculate the number of atoms of each species
    n_atoms_of_species = 0
    number_of_electrons = 0
    do iatom=1, n_atoms
       n_atoms_of_species(species(iatom)) =  &
            n_atoms_of_species(species(iatom)) + 1
       number_of_electrons = number_of_electrons + &
            pseudo(species(iatom))%zval
    end do

    ! loop over species and calculate the correction factor for each one
    e_core = zero

    !calculates the shift for each species
    !the accuracy of v_local(r) = \int_{r
    nfac1= 6
    nfac2= 8
    nfac3= 10

    eshift(:) =zero
    do ispecies = 1, n_species
       beta = pseudo(ispecies)%alpha
       z = pseudo(ispecies)%zval
     if(type_species(ispecies) > 0) then
       if(pseudo(ispecies)%tm_loc_pot == loc_pot) then
          call calc_energy_shift(pseudo(ispecies),nfac1, eshift1, vlocG0_1)
          call calc_energy_shift(pseudo(ispecies),nfac2, eshift2, vlocG0_2)
          call calc_energy_shift(pseudo(ispecies),nfac3, eshift3, vlocG0_3)
          !if(abs(eshift2-eshift3) > eps) write(io_lun,fmt='(10x,"WARNING!!!   eshift")')
          if(inode == ionode.AND.iprint_pseudo>0) then
             write (io_lun, &
                    '(10x,a11,i2,a12,3i3,a19,d15.7,/,63x,d15.7,/,63x,d15.7)') &
                   'ispecies = ', ispecies, '  eshift for', &
                   nfac1, nfac2, nfac3, ' times fine mesh = ', &
                   eshift1, eshift2, eshift3
             write (io_lun,'(55x,a8,d15.7)') ' diff = ', eshift2 - eshift3
          end if
          if(inode == ionode.AND.iprint_pseudo>0) then
             write (io_lun, &
                    '(10x,a11,i2,a12,3i3,a19,d15.7,/,63x,d15.7,/,63x,d15.7)') &
                   'ispecies = ', ispecies, '  eshift for', &
                   nfac1, nfac2, nfac3, ' times fine mesh = ', &
                   vlocG0_1, vlocG0_2, vlocG0_3
             write (io_lun,'(55x,a8,d15.7)') ' diff = ', vlocG0_2 - vlocG0_3
          end if
          eshift(ispecies) = eshift3 - vlocG0_3
       else
          call calc_energy_shift(pseudo(ispecies),nfac1, eshift1)
          call calc_energy_shift(pseudo(ispecies),nfac2, eshift2)
          call calc_energy_shift(pseudo(ispecies),nfac3, eshift3)
          !if(abs(eshift2-eshift3) > eps) write(io_lun,fmt='(10x,"WARNING!!!   eshift")')
          if(inode == ionode.AND.iprint_pseudo>0) then
             write (io_lun, &
                    '(10x,a11,i2,a12,3i3,a19,d15.7,/,63x,d15.7,/,63x,d15.7)') &
                   'ispecies = ', ispecies, '  eshift for', &
                   nfac1, nfac2, nfac3, ' times fine mesh = ', &
                   eshift1, eshift2, eshift3
             write (io_lun,'(55x,a8,d15.7)') ' diff = ', eshift2 - eshift3
          end if
          eshift(ispecies) = eshift3
       end if
     endif ! (type_species(ispecies) > 0) 
    enddo

    do ispecies=1, n_species
       e_core = e_core + real( n_atoms_of_species(ispecies), double ) *  &
            eshift(ispecies)
    end do
    !  I (TM) am not sure whether number_of_electrons should be
    ! used in the following because in Conquest, we may 
    ! perform constant chemical-potential calculations.
    ! (Even in constatnt chemical-potential calculations, we
    !  might ought to assume the system is neutral, otherwise
    !  the total energy would diverge. However, assuming the
    !  extra uniform field, we can calculate the total energy.
    !  I think I should consider it carefully later.)

    e_core = e_core *  number_of_electrons / volume

    call start_timer(tmr_std_allocation)
    deallocate(n_atoms_of_species,eshift, STAT=stat)
    if(stat /= 0) call cq_abort &
         ('deallocation in get_energy_shift',stat)
    call stop_timer(tmr_std_allocation)
    return
  contains

    !sbrt calc_energy_shift -----------------------------------------
    subroutine calc_energy_shift(pseudo, nfac, eshift, vlocG0)

      use datatypes
      use numbers
      use GenComms, only: inode, ionode
      use pseudo_tm_info, only: rad_func, rad_alloc, rad_dealloc
      use global_module, only: iprint_pseudo
      use functions, only: erfc

      implicit none

      !passed
      type(pseudo_info), intent(in) :: pseudo
      integer, intent(in) :: nfac
      real(double), intent(out) :: eshift
      real(double), intent(out), optional :: vlocG0
      !local
      real(double) :: z, beta
      integer :: nmesh, nmesh_vloc
      real(double) :: step, step_fine, rr, rr_imesh
      real(double) :: a,b,c,d, r1,r2,r3,r4
      type(rad_func) :: func1, func2 ! , func, vloc
      integer :: imesh, imesh_rr, imesh_int
      real(double) :: rho, vlocal, intvloc
      real(double), allocatable :: wos(:)
      integer :: stat
      real(double) :: zcore, erfarg
      integer, save:: icall = 0

      z = pseudo%zval
      beta = pseudo%alpha
      icall = icall+1

      if(pseudo%tm_loc_pot==loc_pot) then
         ! Commented out: I think that here we want z not zcore
         !%%! ! First, check integral of local charge
         !%%! nmesh = pseudo%chlocal%n
         !%%! nmesh_vloc = (nmesh-1) * nfac + 1
         !%%! 
         !%%! call rad_alloc(func1, nmesh_vloc)
         !%%! allocate(wos(nmesh_vloc), STAT = stat)
         !%%! func1%cutoff = pseudo%chlocal%cutoff ; func1%delta = pseudo%chlocal%cutoff/(nmesh_vloc-1)
         !%%! 
         !%%! !sets weight of simpson integration
         !%%! call set_wos(nmesh_vloc, func1%delta, wos)
         !%%! !Interpolates core_charge multiplied by 4pi*r**2 and 4pi*r
         !%%! step = pseudo%chlocal%delta
         !%%! zcore = zero
         !%%! do imesh = 1, func1%n
         !%%!    rr = (imesh-1)*func1%delta
         !%%!    imesh_rr = aint(rr/step+very_small) +1
         !%%!    rr_imesh = step * imesh_rr
         !%%! 
         !%%!    if(imesh < func1%n) then
         !%%!       a = ( rr_imesh - rr ) / step
         !%%!       b = one - a
         !%%!       c = a * ( a * a - one ) * step * step / six
         !%%!       d = b * ( b * b - one ) * step * step / six
         !%%! 
         !%%!       r1= pseudo%chlocal%f(imesh_rr)
         !%%!       r2= pseudo%chlocal%f(imesh_rr+1)
         !%%!       r3= pseudo%chlocal%d2(imesh_rr)
         !%%!       r4= pseudo%chlocal%d2(imesh_rr+1)
         !%%! 
         !%%!       rho = a* r1 + b*r2 + c*r3 + d* r4
         !%%!    else
         !%%!       rho = pseudo%chlocal%f(imesh_rr)
         !%%!    endif
         !%%! 
         !%%!    func1%f(imesh) = rho* four * pi * rr **2
         !%%!    zcore=zcore + func1%f(imesh)*wos(imesh)
         !%%! enddo
         !%%! if(inode==ionode.AND.iprint_pseudo>1) write(io_lun,fmt='(10x,"Integral of interpolated chlocal ",f8.3)') zcore
         !%%! zcore = - zcore   ! makes zcore positive
         eshift = zero
         intvloc = zero
         nmesh = pseudo%vlocal%n
         nmesh_vloc = (nmesh-1) * nfac + 1
         call start_timer(tmr_std_allocation)
         allocate(wos(nmesh_vloc), STAT = stat)
         if(stat /= 0) call cq_abort ('ERROR in allocating wos in calc_energy_shift', stat)
         call stop_timer(tmr_std_allocation)
         step = pseudo%vlocal%delta
         step_fine = pseudo%vlocal%cutoff/(nmesh_vloc-1)
         call set_wos(nmesh_vloc, step_fine, wos)
         do imesh = 1, nmesh_vloc
            rr = (imesh-1) * step_fine
            imesh_rr = aint(rr/step+very_small) +1
            rr_imesh = step * imesh_rr

            if(imesh_rr < pseudo%vlocal%n) then ! This is OK - will stop at n-1
               a = ( rr_imesh - rr ) / step
               b = one - a
               c = a * ( a * a - one ) * step * step / six
               d = b * ( b * b - one ) * step * step / six

               r1=pseudo%vlocal%f(imesh_rr)
               r2=pseudo%vlocal%f(imesh_rr+1)
               r3=pseudo%vlocal%d2(imesh_rr)
               r4=pseudo%vlocal%d2(imesh_rr+1)

               vlocal = a* r1 + b*r2 + c*r3 + d* r4
            else
               vlocal = pseudo%vlocal%f(pseudo%vlocal%n) ! Safer this way DRB 2017/04/24
            endif
            ! We need to scale by r
            if(rr > RD_ERR) then
               vlocal = vlocal*rr
            else
               vlocal=zero
            end if


            !Choose whether we use gaussian distribution
            ! for the calculation of energy shift
            eshift = eshift + (vlocal + z ) * four *pi * rr * wos(imesh)
            !intvloc = intvloc + vlocal * four *pi * rr * wos(imesh)
            ! z or zcore?
            ! eshift = eshift + (vlocal + z ) * four *pi * rr * wos(imesh)
            !erf eshift = eshift +  four *pi * rr * wos(imesh) * &
            !erf            (vlocal + zcore * erf(sqrt(beta) *rr) )
            erfarg = one - erfc(sqrt(beta) *rr)
            !erfarg = derf(sqrt(beta) *rr)
            !write(52,*) rr, erf(sqrt(beta) *rr), erfarg
            intvloc = intvloc +  four *pi * rr * wos(imesh) * &
                 (vlocal + z * erfarg )
            !(vlocal + z * erf(sqrt(beta) *rr) )
            if(imesh==nmesh_vloc.AND.inode==ionode.AND.iprint_pseudo>1) &
                 write(io_lun,fmt='(2x,"Sum of local potential and core charge: ",3f12.8)') &
                 vlocal + z, vlocal, z
         enddo
         call start_timer(tmr_std_allocation)
         deallocate(wos, STAT = stat)
         if(stat /= 0) call cq_abort ('ERROR in deallocating wos in calc_energy_shift', stat)
         call stop_timer(tmr_std_allocation)
         !write(io_lun,*) 'Int. vloc: ',intvloc
         !intvloc = intvloc + pi*z/beta
         if(PRESENT(vlocG0)) vlocG0 = intvloc
         !write(io_lun,*) 'Int. vloc, eshift: ',intvloc,eshift, beta, pi*z/beta, z/beta
      else
         nmesh = pseudo%chlocal%n
         nmesh_vloc = (nmesh-1) * nfac + 1

         call rad_alloc(func1, nmesh_vloc)
         call rad_alloc(func2, nmesh_vloc)

         call start_timer(tmr_std_allocation)
         allocate(wos(nmesh_vloc), STAT = stat)
         call stop_timer(tmr_std_allocation)
         !check the integral of chlocal ---------------
         call set_wos(nmesh, pseudo%chlocal%delta, wos)
         rho = zero
         do imesh = 1, pseudo%chlocal%n
            rr = (imesh-1)*pseudo%chlocal%delta
            rho = rho + four*pi*rr**2* pseudo%chlocal%f(imesh) * wos(imesh)
            !         if(inode == ionode.and.icall==1) write(40,*) rr, pseudo%chlocal%f(imesh)
         enddo
         if(inode == ionode.AND.iprint_pseudo>1) write(io_lun,fmt='(2x,"Integral of chlocal = ",2f8.3)') rho, z
         !check the integral of chlocal ---------------

         if(stat /= 0) call cq_abort &
              ('ERROR in allocating wos in calc_energy_shift', stat, nmesh_vloc)

         func1%cutoff = pseudo%chlocal%cutoff ; func1%delta = pseudo%chlocal%cutoff/(nmesh_vloc-1)
         func2%cutoff = pseudo%chlocal%cutoff ; func2%delta = pseudo%chlocal%cutoff/(nmesh_vloc-1)

         !sets weight of simpson integration
         call set_wos(nmesh_vloc, func1%delta, wos)
         !Interpolates core_charge multiplied by 4pi*r**2 and 4pi*r
         step = pseudo%chlocal%delta
         zcore = zero
         do imesh = 1, func1%n
            rr = (imesh-1)*func1%delta
            imesh_rr = aint(rr/step+very_small) +1
            rr_imesh = step * imesh_rr

            if(imesh < func1%n) then
               a = ( rr_imesh - rr ) / step
               b = one - a
               c = a * ( a * a - one ) * step * step / six
               d = b * ( b * b - one ) * step * step / six

               r1= pseudo%chlocal%f(imesh_rr)
               r2= pseudo%chlocal%f(imesh_rr+1)
               r3= pseudo%chlocal%d2(imesh_rr)
               r4= pseudo%chlocal%d2(imesh_rr+1)

               rho = a* r1 + b*r2 + c*r3 + d* r4
            else
               rho = pseudo%chlocal%f(imesh_rr)
            endif

            !         if(inode == ionode.and.icall==1) write(41,1001) rr, rho, r1,r2,r3,r4
1001        format(2d15.6, 5x, 4d12.4)

            func1%f(imesh) = rho* four * pi * rr **2
            func2%f(imesh) = rho* four * pi * rr
            zcore=zcore + func1%f(imesh)*wos(imesh)
         enddo
         if(inode == ionode.AND.iprint_pseudo>1) write(io_lun,fmt='(2x,"Integral of interpolated chlocal ",f8.3)') zcore
         zcore = - zcore   ! makes zcore positive

         !Makes vlocal & integrates vlocal
         ! Vlocal(r) = {\int_{0}^{r} ds 4pi s**2 rho(s)}/ r + &
         !                  \int_{r}^{\inf} ds 4pi s rho(s)
         !Here, I define the 
         eshift = zero
         intvloc = zero
         do imesh = 1, nmesh_vloc
            rr = (imesh-1) * func1%delta

            vlocal = zero    ! vlocal = v_local(r) * r
            do imesh_int = 1, nmesh_vloc
               if(imesh_int < imesh) then
                  vlocal = vlocal + func1%f(imesh_int) * wos(imesh_int)
               else
                  vlocal = vlocal + func2%f(imesh_int) *rr * wos(imesh_int)
               endif
            enddo
            !if(rr>RD_ERR) then
            !   write(51,*) rr,vlocal/rr
            !else
            !   write(51,*) rr,func2%f(1)*wos(1)
            !end if
            !Choose whether we use gaussian distribution
            ! for the calculation of energy shift
            eshift = eshift + (vlocal + zcore ) * four *pi * rr * wos(imesh)
            intvloc = intvloc + vlocal * four *pi * rr * wos(imesh)
            ! z or zcore?
            ! eshift = eshift + (vlocal + z ) * four *pi * rr * wos(imesh)
            !erf eshift = eshift +  four *pi * rr * wos(imesh) * &
            !erf            (vlocal + zcore * erf(sqrt(beta) *rr) )
            if(imesh==nmesh_vloc.AND.inode==ionode.AND.iprint_pseudo>1) &
                 write(io_lun,fmt='(2x,"Sum of local potential and core charge: ",3f8.3)') &
                 vlocal + zcore, vlocal, zcore
         enddo
         !Choose whether we use gaussian distribution
         ! for the calculation of energy shift
         !erf eshift = eshift + pi*zcore/beta

         call rad_dealloc(func1)
         call rad_dealloc(func2)
         call start_timer(tmr_std_allocation)
         deallocate(wos, STAT = stat)
         if(stat /= 0) call cq_abort &
              ('ERROR in deallocating wos in calc_energy_shift', stat)
         call stop_timer(tmr_std_allocation)
         !write(io_lun,*) 'Int. vloc, eshift: ',intvloc,eshift
      end if ! tm_loc_pot == loc_pot
      return
    end subroutine calc_energy_shift

    !sbrt set_wos ----------------------
    subroutine set_wos(n, delta, wos)
      implicit none
      integer, intent(in):: n
      real(double), intent(in) :: delta
      real(double), intent(out) :: wos(1:n)
      real(double) :: one_third, two_third, four_third
      integer :: imesh

      if(n < 4) call cq_abort('ERROR n in set_wos = ',n)
      one_third = one/three
      two_third = two/three
      four_third = four/three

      do imesh = 2, n-1, 2
         wos(imesh  ) = four_third * delta
         wos(imesh+1) =  two_third * delta
      enddo
      wos( 1     ) =  one_third * delta
      wos( n     ) =  one_third * delta
      return
    end subroutine set_wos

  end subroutine get_energy_shift
!!***

! -----------------------------------------------------------
! Subroutine make_neutral_atom
! -----------------------------------------------------------

!!****f* pseudo_tm_module/make_neutral_atom *
!!
!!  NAME 
!!   make_neutral_atom
!!  USAGE
!!   
!!  PURPOSE
!!   Prepares arrays and matrix elements for neutral atom potential
!!  INPUTS
!!   None
!!   
!!  USES
!!   
!!  AUTHOR
!!   N. Watanabe (Mizuho) with TM and DRB
!!  CREATION DATE
!!   2014
!!  MODIFICATION HISTORY
!!   2015/11/09 08:29 dave
!!    Moved to pseudo_tm_module and included in main trunk
!!   2016/01/07 08:37 dave
!!    Changed splint calls into explicit calculations (faster, removed some bug)  
!!  SOURCE
!!  
  subroutine make_neutral_atom

    use datatypes
    use numbers
    use dimens, only: volume
    use species_module, ONLY: n_species
    use pseudo_tm_info, ONLY: rad_func, pseudo, loc_pot, rad_alloc, rad_dealloc
    use atomic_density, ONLY: atomic_density_table
    use spline_module, ONLY: spline
    use functions, ONLY: j0
    use pseudopotential_common, ONLY: flag_neutral_atom_projector, maxL_neutral_atom_projector, &
         numN_neutral_atom_projector
    use species_module, ONLY: napf_species, npao_species
    use pao_format, ONLY: pao, paoVNA, paopao
    use bessel_integrals, only: general_bessel

    implicit none

    integer      :: npoint
    integer :: ik
    integer :: isp
    integer :: ir, ir1, ir2, j, i_pjna, l, n, jj
    
    real(double) :: cutoff
    real(double) :: delta
    real(double) :: r, r_1, r_2, yp1, ypn, r1, r2, r3, r4
    real(double) :: rr, a, b, c, d, step, alpha, beta, gamma
    real(double) :: chatom_at, chlocal_at, chna_at
    real(double) :: vlocal_at, vchatom_at, vna_here, pao_here, dvna_here, dpao_here
    real(double) :: integ_inside, integ_outside
    real(double) :: k, overlap, norm, sum_chna
    real(double) :: sqrt_two_pi ! sqrt(two/pi)

    logical :: range_flag
    
    ! parameters of integration, it can be changed.
    integer, parameter :: k_length = 2048 ! size of radial function in reciprocal space
    real(double), parameter :: k_cutoff = 16.0_double

    real(double), allocatable :: Rad(:)
    real(double), allocatable :: Radbar(:,:)
    
    integer :: j_pjna
    integer :: m

    ! Bessel function variables
    integer :: zeroind, bessel_order
    real(double), allocatable, dimension(:,:) :: zeroloc
    real(double), allocatable, dimension(:,:) :: qvec
    real(double) :: rm, rn, jl, jm, jn, x, fac
    logical :: not_done

    sqrt_two_pi = sqrt(two)/sqrt_pi
    if(flag_neutral_atom_projector) then
       allocate(napf_species(n_species),paoVNA(n_species),paopao(n_species))
       allocate(zeroloc(15,0:6),qvec(15,0:8))
       zeroloc = zero
       qvec = zero
       ! Find Bessel function zeroes
       npoint = 1001
       do bessel_order=0,6!lmax_ps+lmax_pao ! Bessel order
          cutoff = 100.0_double
          delta = cutoff/real(npoint-1,double)
          zeroind = 1
          do ir = 2,npoint
             r = ir*delta
             rm = (ir-1)*delta
             jl = general_bessel(r,bessel_order)
             jm = general_bessel(rm,bessel_order)
             !jj = spherical_bessels(r,bessel_order)
             !jm = spherical_bessels(rm,bessel_order)
             if(jl*jm<zero) then ! Bracketed a zero
                not_done = .true.
                do while(not_done)
                   rn = half*(r+rm)
                   jn = general_bessel(rn,bessel_order)
                   !jn = spherical_bessels(rn,bessel_order)
                   if(jn*jm<zero) then !
                      r = rn
                      jl = jn
                   else if(jn*jl<zero) then
                      rm = rn
                      jm = jn
                   else if(abs(jn)<1e-8_double) then
                      not_done = .false.
                   end if
                   if(abs(r-rm)<1e-8) not_done = .false.
                end do
                zeroloc(zeroind,bessel_order) = rn
                !write(io_lun,*) '# zero is at ',zeroloc(zeroind)
                zeroind = zeroind+1
                if(zeroind>15) exit
             end if
          end do
       end do
    end if
    ! Loop over species
    do isp=1, n_species
       if(flag_neutral_atom_projector) then
          cutoff = pseudo(isp)%vna%cutoff
          do bessel_order=0,6
             do ir = 1, zeroind-1
                qvec(ir, bessel_order) = zeroloc(ir,bessel_order)/cutoff
             end do
          end do
       end if
       if( pseudo(isp)%tm_loc_pot==loc_pot ) then ! ABINIT PP
          cutoff = atomic_density_table(isp)%cutoff
          delta  = atomic_density_table(isp)%delta
       else ! SIESTA PP
          cutoff = max( atomic_density_table(isp)%cutoff, pseudo(isp)%chlocal%cutoff )
          delta  = min( atomic_density_table(isp)%delta, pseudo(isp)%chlocal%delta )
       end if
       npoint = int(cutoff/delta) + 1
       !-------------------------------------------------------
       ! UNNECESSARY ROUTINES 
       !-------------------------------------------------------
       if(.false.) then ! remove this for now
          !-------------------------------------------------------
          ! calculate chna(r) - atomic charge
          pseudo(isp)%chna%cutoff = cutoff
          pseudo(isp)%chna%delta = delta
          call rad_alloc( pseudo(isp)%chna, npoint )

          norm = zero
          sum_chna = zero
          do ir=1, pseudo(isp)%chna%n
             r = (ir-1)*pseudo(isp)%chna%delta
             chatom_at = zero
             if( r<atomic_density_table(isp)%cutoff ) then

                step = atomic_density_table(isp)%delta
                j = aint( r / step ) + 1
                if(j+1<=atomic_density_table(isp)%length) then
                   rr = real(j,double) * step
                   a = ( rr - r ) / step
                   b = one - a
                   c = a * ( a * a - one ) * step * step / six
                   d = b * ( b * b - one ) * step * step / six

                   r1=atomic_density_table(isp)%table(j)
                   r2=atomic_density_table(isp)%table(j+1)
                   r3=atomic_density_table(isp)%d2_table(j)
                   r4=atomic_density_table(isp)%d2_table(j+1)

                   chatom_at =  a * r1 + b * r2 + c * r3 + d * r4
                end if
                !call splint(atomic_density_table(isp)%delta, &
                !     atomic_density_table(isp)%table, &
                !     atomic_density_table(isp)%d2_table, &
                !     atomic_density_table(isp)%length, &
                !     r, chatom_at, range_flag )
             end if

             if( pseudo(isp)%tm_loc_pot==loc_pot ) then ! ABINIT PP
                chlocal_at = minus_one * pseudo(isp)%zval &
                     * pseudo(isp)%prefac * exp(-pseudo(isp)%alpha*r*r)
             else ! SIESTA PP
                chlocal_at = zero
                if( r<pseudo(isp)%chlocal%cutoff ) then
                   step = pseudo(isp)%chlocal%delta
                   j = aint( r / step ) + 1
                   if(j+1<=pseudo(isp)%chlocal%n) then
                      rr = real(j,double) * step
                      a = ( rr - r ) / step
                      b = one - a
                      c = a * ( a * a - one ) * step * step / six
                      d = b * ( b * b - one ) * step * step / six

                      r1=pseudo(isp)%chlocal%f(j)
                      r2=pseudo(isp)%chlocal%f(j+1)
                      r3=pseudo(isp)%chlocal%d2(j)
                      r4=pseudo(isp)%chlocal%d2(j+1)

                      chlocal_at =  a * r1 + b * r2 + c * r3 + d * r4
                   end if
                   !call splint( pseudo(isp)%chlocal%delta, &
                   !     pseudo(isp)%chlocal%f, &
                   !     pseudo(isp)%chlocal%d2, &
                   !     pseudo(isp)%chlocal%n, &
                   !     r, chlocal_at, range_flag )
                end if
             end if

             pseudo(isp)%chna%f(ir) = chatom_at + chlocal_at
             norm = norm + r*r* pseudo(isp)%chna%f(ir)
          end do
          sum_chna = norm*four*pi* delta
          yp1 = zero!(pseudo(isp)%chna%f(2)-pseudo(isp)%chna%f(1)) &
          !/pseudo(isp)%chna%delta
          ypn = zero!(pseudo(isp)%chna%f(pseudo(isp)%chna%n)-pseudo(isp)%chna%f(pseudo(isp)%chna%n-1)) &
          !/pseudo(isp)%chna%delta

          call spline( pseudo(isp)%chna%n, pseudo(isp)%chna%delta, &
               pseudo(isp)%chna%f, yp1, ypn, pseudo(isp)%chna%d2 )


          !-------------------------------------------------------
          ! calculate Vna(r) - neutral atom potential
          pseudo(isp)%vna%cutoff = cutoff
          pseudo(isp)%vna%delta = delta
          call rad_alloc( pseudo(isp)%vna, npoint )
          do ir1=1, pseudo(isp)%chna%n
             r_1 = real(ir1-1,double)*pseudo(isp)%chna%delta

             integ_inside = zero
             do ir2=1, ir1
                r_2 = real(ir2-1,double)*pseudo(isp)%chna%delta

                step = pseudo(isp)%chna%delta
                j = aint( r_2 / step ) + 1
                chna_at = zero
                if(j+1<=pseudo(isp)%chna%n) then
                   rr = real(j,double) * step
                   a = ( rr - r_2 ) / step
                   b = one - a
                   c = a * ( a * a - one ) * step * step / six
                   d = b * ( b * b - one ) * step * step / six

                   r1=pseudo(isp)%chna%f(j)
                   r2=pseudo(isp)%chna%f(j+1)
                   r3=pseudo(isp)%chna%d2(j)
                   r4=pseudo(isp)%chna%d2(j+1)

                   chna_at =  a * r1 + b * r2 + c * r3 + d * r4
                end if
                !call splint( pseudo(isp)%chna%delta, pseudo(isp)%chna%f, &
                !     pseudo(isp)%chna%d2, pseudo(isp)%chna%n, r2, chna_at, range_flag )
                if( ir2 == 1 .or. ir2 == ir1 ) then
                   chna_at = chna_at*half
                end if
                integ_inside = integ_inside + r_2*r_2 * chna_at
             end do
             if( r_1 == zero ) then
                integ_inside = zero
             else
                integ_inside = integ_inside * (four*pi/r_1)*pseudo(isp)%chna%delta
             end if

             integ_outside = zero
             do ir2=ir1, pseudo(isp)%chna%n
                r_2 = real(ir2-1,double)*pseudo(isp)%chna%delta

                step = pseudo(isp)%chna%delta
                j = aint( r_2 / step ) + 1
                chna_at = zero
                if(j+1<=pseudo(isp)%chna%n) then
                   rr = real(j,double) * step
                   a = ( rr - r_2 ) / step
                   b = one - a
                   c = a * ( a * a - one ) * step * step / six
                   d = b * ( b * b - one ) * step * step / six

                   r1=pseudo(isp)%chna%f(j)
                   r2=pseudo(isp)%chna%f(j+1)
                   r3=pseudo(isp)%chna%d2(j)
                   r4=pseudo(isp)%chna%d2(j+1)

                   chna_at =  a * r1 + b * r2 + c * r3 + d * r4
                end if
                !call splint( pseudo(isp)%chna%delta, pseudo(isp)%chna%f, &
                !     pseudo(isp)%chna%d2, pseudo(isp)%chna%n, r2, chna_at, range_flag )
                if( ir2 == ir1 .or. ir2 == pseudo(isp)%chna%n ) then
                   chna_at = chna_at*half
                end if
                integ_outside = integ_outside + r_2 * chna_at
             end do
             integ_outside = integ_outside * four*pi*pseudo(isp)%chna%delta

             if( pseudo(isp)%tm_loc_pot==loc_pot ) then ! ABINIT PP
                vlocal_at = zero
                if( r_1<pseudo(isp)%vlocal%cutoff ) then

                   step = pseudo(isp)%vlocal%delta
                   j = aint( r_1 / step ) + 1
                   if(j+1<=pseudo(isp)%vlocal%n) then
                      rr = real(j,double) * step
                      a = ( rr - r_1 ) / step
                      b = one - a
                      c = a * ( a * a - one ) * step * step / six
                      d = b * ( b * b - one ) * step * step / six

                      r1=pseudo(isp)%vlocal%f(j)
                      r2=pseudo(isp)%vlocal%f(j+1)
                      r3=pseudo(isp)%vlocal%d2(j)
                      r4=pseudo(isp)%vlocal%d2(j+1)

                      vlocal_at =  a * r1 + b * r2 + c * r3 + d * r4
                   end if
                   !call splint( pseudo(isp)%vlocal%delta, pseudo(isp)%vlocal%f, &
                   !     pseudo(isp)%vlocal%d2, &
                   !     pseudo(isp)%vlocal%n, &
                   !     r1, vlocal_at, range_flag )
                end if
                pseudo(isp)%vna%f(ir1) = integ_inside+integ_outside + vlocal_at
             else ! SIESTA PP
                pseudo(isp)%vna%f(ir1) = integ_inside+integ_outside
             end if
          end do

          yp1 = (pseudo(isp)%vna%f(2)-pseudo(isp)%vna%f(1)) &
               /pseudo(isp)%vna%delta
          ypn = (pseudo(isp)%vna%f(pseudo(isp)%vna%n)-pseudo(isp)%vna%f(pseudo(isp)%vna%n-1)) &
               /pseudo(isp)%vna%delta

          call spline( pseudo(isp)%vna%n, pseudo(isp)%vna%delta, &
               pseudo(isp)%vna%f, yp1, ypn, pseudo(isp)%vna%d2 )

          !-------------------------------------------------------
          ! calculate Vatom(r) - atomic potential
          norm = zero
          do ir1=1, atomic_density_table(isp)%length
             r_1 = (ir1-1)*atomic_density_table(isp)%delta

             integ_inside = zero
             do ir2=1, ir1
                r_2 = (ir2-1)*atomic_density_table(isp)%delta

                chatom_at = zero
                step = atomic_density_table(isp)%delta
                j = aint( r_2 / step ) + 1
                if(j+1<=atomic_density_table(isp)%length) then
                   rr = real(j,double) * step
                   a = ( rr - r_2 ) / step
                   b = one - a
                   c = a * ( a * a - one ) * step * step / six
                   d = b * ( b * b - one ) * step * step / six

                   r1=atomic_density_table(isp)%table(j)
                   r2=atomic_density_table(isp)%table(j+1)
                   r3=atomic_density_table(isp)%d2_table(j)
                   r4=atomic_density_table(isp)%d2_table(j+1)

                   chatom_at =  a * r1 + b * r2 + c * r3 + d * r4
                end if
                !write(io_lun,*) 'splint6'
                !call splint(atomic_density_table(isp)%delta, &
                !     atomic_density_table(isp)%table, &
                !     atomic_density_table(isp)%d2_table, &
                !     atomic_density_table(isp)%length, &
                !     r2, chatom_at, range_flag )

                if( ir2 == 1 .or. ir2 == ir1 ) then
                   chatom_at = chatom_at*half
                end if
                integ_inside = integ_inside + r_2*r_2 * chatom_at
             end do
             if( r_1 == zero ) then
                integ_inside = zero
             else
                integ_inside = integ_inside * (four*pi/r_1)*atomic_density_table(isp)%delta
             end if

             integ_outside = zero
             do ir2=ir1, atomic_density_table(isp)%length
                r_2 = (ir2-1)*atomic_density_table(isp)%delta

                chatom_at = zero
                step = atomic_density_table(isp)%delta
                j = aint( r_2 / step ) + 1
                if(j+1<=atomic_density_table(isp)%length) then
                   rr = real(j,double) * step
                   a = ( rr - r_2 ) / step
                   b = one - a
                   c = a * ( a * a - one ) * step * step / six
                   d = b * ( b * b - one ) * step * step / six

                   r1=atomic_density_table(isp)%table(j)
                   r2=atomic_density_table(isp)%table(j+1)
                   r3=atomic_density_table(isp)%d2_table(j)
                   r4=atomic_density_table(isp)%d2_table(j+1)

                   chatom_at =  a * r1 + b * r2 + c * r3 + d * r4
                end if
                !write(io_lun,*) 'splint7'
                !call splint(atomic_density_table(isp)%delta, &
                !     atomic_density_table(isp)%table, &
                !     atomic_density_table(isp)%d2_table, &
                !     atomic_density_table(isp)%length, &
                !     r2, chatom_at, range_flag )

                if( ir2 == ir1 .or. ir2 == atomic_density_table(isp)%length ) then
                   chatom_at = chatom_at*half
                end if
                integ_outside = integ_outside + r_2 * chatom_at
             end do
             integ_outside = integ_outside * four*pi*atomic_density_table(isp)%delta

             vchatom_at = integ_inside+integ_outside
             norm = norm + (vchatom_at*r_1 - pseudo(isp)%zval)*four*pi*r_1*atomic_density_table(isp)%delta
          end do
          pseudo(isp)%eshift = norm
       end if ! .false. - removing unnecessary routines
       !-------------------------------------------------------
       ! END OF UNNECESSARY ROUTINES 
       !-------------------------------------------------------

       if( flag_neutral_atom_projector ) then
          !-------------------------------------------------------
          ! calculate pjna(r)
          !
          ! We need (at least) the same number of projectors as PAOs
          ! Following simple 1- and 2-centre tests we may add projectors
          ! for the 3-centre integrals; it would be easiest to add these
          ! on-top of the other functions
          !
          paoVNA(isp)%greatest_angmom = pao(isp)%greatest_angmom
          paoVNA(isp)%count = pao(isp)%count
          paopao(isp)%greatest_angmom = pao(isp)%greatest_angmom
          paopao(isp)%count = pao(isp)%count
          i_pjna = 0
          napf_species(isp) = 0
          do l=0, maxL_neutral_atom_projector
             do n=1, numN_neutral_atom_projector(l)
                i_pjna = i_pjna + 1
             end do
          end do
          napf_species(isp) = i_pjna
          pseudo(isp)%n_pjna = i_pjna
          npoint = pseudo(isp)%vna%n          
          allocate(paoVNA(isp)%angmom(0:pao(isp)%greatest_angmom))
          allocate(paopao(isp)%angmom(0:pao(isp)%greatest_angmom))
          allocate( pseudo(isp)%pjna(pseudo(isp)%n_pjna) )
          allocate( pseudo(isp)%pjna_l(pseudo(isp)%n_pjna) )
          allocate( pseudo(isp)%pjna_n(pseudo(isp)%n_pjna) )
          allocate( pseudo(isp)%pjna_ekb(pseudo(isp)%n_pjna) )
          allocate( Rad(npoint) )
          allocate( Radbar(npoint,maxval(numN_neutral_atom_projector(:))) )

          ! 1- and 2-centre integrals
          i_pjna = 0
          do l=0, pao(isp)%greatest_angmom
             paoVNA(isp)%angmom(l)%n_zeta_in_angmom = pao(isp)%angmom(l)%n_zeta_in_angmom
             paopao(isp)%angmom(l)%n_zeta_in_angmom = pao(isp)%angmom(l)%n_zeta_in_angmom
             allocate(paoVNA(isp)%angmom(l)%zeta(pao(isp)%angmom(l)%n_zeta_in_angmom))
             allocate(paopao(isp)%angmom(l)%zeta(pao(isp)%angmom(l)%n_zeta_in_angmom))
             do n=1, pao(isp)%angmom(l)%n_zeta_in_angmom
                i_pjna = i_pjna + 1
                paoVNA(isp)%angmom(l)%zeta(n)%length = pao(isp)%angmom(l)%zeta(n)%length
                paoVNA(isp)%angmom(l)%zeta(n)%cutoff = pao(isp)%angmom(l)%zeta(n)%cutoff
                paoVNA(isp)%angmom(l)%zeta(n)%delta = pao(isp)%angmom(l)%zeta(n)%delta
                allocate(paoVNA(isp)%angmom(l)%zeta(n)%table(paoVNA(isp)%angmom(l)%zeta(n)%length))
                allocate(paoVNA(isp)%angmom(l)%zeta(n)%table2(paoVNA(isp)%angmom(l)%zeta(n)%length))
                paoVNA(isp)%angmom(l)%zeta(n)%table = zero
                paoVNA(isp)%angmom(l)%zeta(n)%table2 = zero
                ! pao-pao
                paopao(isp)%angmom(l)%zeta(n)%length = pseudo(isp)%vna%n !pao(isp)%angmom(l)%zeta(n)%length
                paopao(isp)%angmom(l)%zeta(n)%cutoff = pseudo(isp)%vna%cutoff!pao(isp)%angmom(l)%zeta(n)%cutoff
                paopao(isp)%angmom(l)%zeta(n)%delta  = pseudo(isp)%vna%delta!pao(isp)%angmom(l)%zeta(n)%delta
                allocate(paopao(isp)%angmom(l)%zeta(n)%table(paopao(isp)%angmom(l)%zeta(n)%length))
                allocate(paopao(isp)%angmom(l)%zeta(n)%table2(paopao(isp)%angmom(l)%zeta(n)%length))
                paopao(isp)%angmom(l)%zeta(n)%table = zero
                paopao(isp)%angmom(l)%zeta(n)%table2 = zero
             end do
          end do
          ! 3-centre integrals via projectors
          i_pjna = 0
          napf_species(isp) = 0
          do l=0, maxL_neutral_atom_projector
             do n=1, numN_neutral_atom_projector(l)
                i_pjna = i_pjna + 1
                pseudo(isp)%pjna_l(i_pjna) = l
                pseudo(isp)%pjna_n(i_pjna) = n
                pseudo(isp)%pjna(i_pjna)%cutoff = pseudo(isp)%vna%cutoff
                pseudo(isp)%pjna(i_pjna)%delta  = pseudo(isp)%vna%delta
                call rad_alloc( pseudo(isp)%pjna(i_pjna), pseudo(isp)%vna%n )
                napf_species(isp) = napf_species(isp) + (2*l+1)
             end do
          end do
          ! Create Vna phi for each PAO for each species
          do l=0, pao(isp)%greatest_angmom
             do n=1, pao(isp)%angmom(l)%n_zeta_in_angmom
                !   ! |Vna_ln> = Vna(r) |PAO'_ln>
                step = pseudo(isp)%vna%delta
                do j=1,paoVNA(isp)%angmom(l)%zeta(n)%length
                   r = paoVNA(isp)%angmom(l)%zeta(n)%delta * real(j-1,double)
                   jj = aint(r/step) + 1
                   if(jj+1<=pseudo(isp)%vna%n) then
                      rr = real(jj,double) * step
                      a = ( rr - r ) / step
                      b = one - a
                      c = a * ( a * a - one ) * step * step / six
                      d = b * ( b * b - one ) * step * step / six                   
                      r1=pseudo(isp)%vna%f(jj)
                      r2=pseudo(isp)%vna%f(jj+1)
                      r3=pseudo(isp)%vna%d2(jj)
                      r4=pseudo(isp)%vna%d2(jj+1)
                      vna_here = a * r1 + b * r2 + c * r3 + d * r4
                      paoVNA(isp)%angmom(l)%zeta(n)%table(j) = pao(isp)%angmom(l)%zeta(n)%table(j) * vna_here
                   else
                      paoVNA(isp)%angmom(l)%zeta(n)%table(j) = zero
                   end if
                end do
                ! Spline paoVNA
                delta = pao(isp)%angmom(l)%zeta(n)%delta
                npoint = pao(isp)%angmom(l)%zeta(n)%length
                yp1 = (paoVNA(isp)%angmom(l)%zeta(n)%table(2) &
                     - paoVNA(isp)%angmom(l)%zeta(n)%table(1)       )/delta
                ypn = (paoVNA(isp)%angmom(l)%zeta(n)%table(npoint) &
                     - paoVNA(isp)%angmom(l)%zeta(n)%table(npoint-1))/delta
                call spline( npoint, delta, &
                     paoVNA(isp)%angmom(l)%zeta(n)%table, yp1, ypn, &
                     paoVNA(isp)%angmom(l)%zeta(n)%table2 )
                ! Create paopao: all PAOs on the same grid (for ease of multiplication)
                step = pao(isp)%angmom(l)%zeta(n)%delta
                do j=1,pseudo(isp)%vna%n
                   r = pseudo(isp)%vna%delta * real(j-1,double)
                   jj = aint(r/step) + 1
                   if(jj+1<=pao(isp)%angmom(l)%zeta(n)%length) then
                      rr = real(jj,double) * step
                      a = ( rr - r ) / step
                      b = one - a
                      c = a * ( a * a - one ) * step * step / six
                      d = b * ( b * b - one ) * step * step / six                   
                      r1=pao(isp)%angmom(l)%zeta(n)%table(jj)
                      r2=pao(isp)%angmom(l)%zeta(n)%table(jj+1)
                      r3=pao(isp)%angmom(l)%zeta(n)%table2(jj)
                      r4=pao(isp)%angmom(l)%zeta(n)%table2(jj+1)
                      pao_here = a * r1 + b * r2 + c * r3 + d * r4
                      paopao(isp)%angmom(l)%zeta(n)%table(j) = pao_here
                   else
                      paopao(isp)%angmom(l)%zeta(n)%table(j) = zero
                   end if
                end do
                delta = paopao(isp)%angmom(l)%zeta(n)%delta
                npoint = paopao(isp)%angmom(l)%zeta(n)%length
                yp1 = (paopao(isp)%angmom(l)%zeta(n)%table(2) &
                     - paopao(isp)%angmom(l)%zeta(n)%table(1)       )/delta
                ypn = (paopao(isp)%angmom(l)%zeta(n)%table(npoint) &
                     - paopao(isp)%angmom(l)%zeta(n)%table(npoint-1))/delta
                call spline( npoint, delta, &
                     paopao(isp)%angmom(l)%zeta(n)%table, yp1, ypn, &
                     paopao(isp)%angmom(l)%zeta(n)%table2 )
             end do
          end do
          ! 3-centre integrals
          delta = pseudo(isp)%vna%delta
          cutoff = pseudo(isp)%vna%cutoff
          npoint = pseudo(isp)%vna%n
          Radbar = zero
          do i_pjna=1, pseudo(isp)%n_pjna
             n = pseudo(isp)%pjna_n(i_pjna)
             l = pseudo(isp)%pjna_l(i_pjna)
             step = pseudo(isp)%pjna(i_pjna)%delta
             a = (four/cutoff)**2
             Rad = zero
             do ir=1, npoint
                r = (ir-1)*delta
                !if(l>0.AND.ir>1) then
                !   fac = one/(r**l)
                !else
                   fac = one
                !end if
                Rad(ir) = fac*general_bessel(qvec(n,l)*r,l)!r**real(n-1,double) * exp(-a*r**2) ! n-1
                !Rad(ir) =  r**real(n-1,double)*exp(-a*r**2)
             end do
             ! Normalise
             norm = overlapRadials(Rad,Rad,delta,l)
             Rad(:) = Rad(:) * sqrt(one/norm)

             ! |R'_ln>=|R_ln>
             Radbar(:,n) = Rad(:)

             do j_pjna=1, i_pjna-1
                if( pseudo(isp)%pjna_l(j_pjna) /= l ) then
                   cycle
                end if
                m = pseudo(isp)%pjna_n(j_pjna)

                ! calc <R_ln|VNA_lm>
                overlap = overlapRadials(Rad,pseudo(isp)%pjna(j_pjna)%f,delta,l)
                ! subtract |R_ln> by |R'_lm>*(<R_ln|VNA_lm>*(1/C_lm)
                Radbar(:,n) = Radbar(:,n) - Radbar(:,m) * &
                     (pseudo(isp)%pjna_ekb(j_pjna)*overlap)  !DRB 2017/06/13
             end do

             ! |Vna_ln> = Vna(r) |R'_ln>
             pseudo(isp)%pjna(i_pjna)%f(:) = pseudo(isp)%vna%f(:) * Radbar(:,n)
             ! calc <R'_ln|Vna_i>
             overlap = overlapRadials(Radbar(:,n),pseudo(isp)%pjna(i_pjna)%f(:),delta,l)

             ! 1/C_ln = 1/<R'_i|Vna_ln>
             pseudo(isp)%pjna_ekb(i_pjna) = one/overlap
          end do
          ! make splines
          do i_pjna=1, pseudo(isp)%n_pjna
             yp1 = (pseudo(isp)%pjna(i_pjna)%f(2) &
                  - pseudo(isp)%pjna(i_pjna)%f(1)       )/delta
             ypn = (pseudo(isp)%pjna(i_pjna)%f(npoint) &
                  - pseudo(isp)%pjna(i_pjna)%f(npoint-1))/delta
             call spline( npoint, delta, &
                  pseudo(isp)%pjna(i_pjna)%f, yp1, ypn, &
                  pseudo(isp)%pjna(i_pjna)%d2 )
          end do

          deallocate( Rad, Radbar )

       end if
       !-------------------------------------------------------
       ! calculate chatom(k) - atomic charge
       atomic_density_table(isp)%k_length = k_length
       atomic_density_table(isp)%k_delta  = k_cutoff/atomic_density_table(isp)%k_length
       allocate( atomic_density_table(isp)%k_table(0:atomic_density_table(isp)%k_length) )
       atomic_density_table(isp)%k_table(:) = zero

       do ir=1, npoint
          r = (ir-1)*delta

          chatom_at = zero
          step = atomic_density_table(isp)%delta
          j = aint( r / step ) + 1
          if(j+1<=atomic_density_table(isp)%length) then
             rr = real(j,double) * step
             a = ( rr - r ) / step
             b = one - a
             c = a * ( a * a - one ) * step * step / six
             d = b * ( b * b - one ) * step * step / six

             r1=atomic_density_table(isp)%table(j)
             r2=atomic_density_table(isp)%table(j+1)
             r3=atomic_density_table(isp)%d2_table(j)
             r4=atomic_density_table(isp)%d2_table(j+1)

             chatom_at =  a * r1 + b * r2 + c * r3 + d * r4
          end if
          !write(io_lun,*) 'splint8'
          !call splint( atomic_density_table(isp)%delta, &
          !     atomic_density_table(isp)%table, &
          !     atomic_density_table(isp)%d2_table, &
          !     atomic_density_table(isp)%length, &
          !     r, chatom_at,range_flag)

          do ik=0, atomic_density_table(isp)%k_length
             k = ik*atomic_density_table(isp)%k_delta

             atomic_density_table(isp)%k_table(ik) = &
                  atomic_density_table(isp)%k_table(ik) &
                  + sqrt_two_pi*delta*r*r*j0(k*r) * chatom_at
          end do
       end do
       !-------------------------------------------------------
       ! UNNECESSARY ROUTINES 
       !-------------------------------------------------------
       if(.false.) then
          ! Local
          !-------------------------------------------------------
          ! calculate FT
          pseudo(isp)%chlocal%k_length = k_length
          pseudo(isp)%chlocal%k_delta  = k_cutoff/pseudo(isp)%chlocal%k_length
          allocate( pseudo(isp)%chlocal%k_table(0:pseudo(isp)%chlocal%k_length) )
          pseudo(isp)%chlocal%k_table(:) = zero

          do ir=1, npoint
             r = (ir-1)*delta

             chatom_at = zero
             step = pseudo(isp)%chlocal%delta
             j = aint( r / step ) + 1
             if(j+1<=pseudo(isp)%chlocal%n) then
                rr = real(j,double) * step
                a = ( rr - r ) / step
                b = one - a
                c = a * ( a * a - one ) * step * step / six
                d = b * ( b * b - one ) * step * step / six

                r1=pseudo(isp)%chlocal%f(j)
                r2=pseudo(isp)%chlocal%f(j+1)
                r3=pseudo(isp)%chlocal%d2(j)
                r4=pseudo(isp)%chlocal%d2(j+1)

                chatom_at =  a * r1 + b * r2 + c * r3 + d * r4
             end if
             do ik=0, pseudo(isp)%chlocal%k_length
                k = ik*pseudo(isp)%chlocal%k_delta

                pseudo(isp)%chlocal%k_table(ik) = &
                     pseudo(isp)%chlocal%k_table(ik) &
                     + sqrt_two_pi*delta*r*r*j0(k*r) * chatom_at
             end do
          end do
          ! NA
          !-------------------------------------------------------
          ! calculate FT
          pseudo(isp)%chna%k_length = k_length
          pseudo(isp)%chna%k_delta  = k_cutoff/pseudo(isp)%chna%k_length
          allocate( pseudo(isp)%chna%k_table(0:pseudo(isp)%chna%k_length) )
          pseudo(isp)%chna%k_table(:) = zero

          do ir=1, npoint
             r = (ir-1)*delta

             chatom_at = zero
             step = pseudo(isp)%chna%delta
             j = aint( r / step ) + 1
             if(j+1<=pseudo(isp)%chna%n) then
                rr = real(j,double) * step
                a = ( rr - r ) / step
                b = one - a
                c = a * ( a * a - one ) * step * step / six
                d = b * ( b * b - one ) * step * step / six

                r1=pseudo(isp)%chna%f(j)
                r2=pseudo(isp)%chna%f(j+1)
                r3=pseudo(isp)%chna%d2(j)
                r4=pseudo(isp)%chna%d2(j+1)

                chatom_at =  a * r1 + b * r2 + c * r3 + d * r4
             end if
             do ik=0, pseudo(isp)%chna%k_length
                k = ik*pseudo(isp)%chna%k_delta

                pseudo(isp)%chna%k_table(ik) = &
                     pseudo(isp)%chna%k_table(ik) &
                     + sqrt_two_pi*delta*r*r*j0(k*r) * chatom_at
             end do
          end do
       end if ! End of if(.false.) commenting out unnecessary routines
       !-------------------------------------------------------
       ! END OF UNNECESSARY ROUTINES 
       !-------------------------------------------------------
    end do ! n_species

  end subroutine make_neutral_atom
  !!***
  
  ! -----------------------------------------------------------
  ! Function overlapRadials
  ! -----------------------------------------------------------

  !!****f* pseudo_tm_module/overlapRadials *
  !!
  !!  NAME
  !!   overlapRadials
  !!  USAGE
  !!   
  !!  PURPOSE
  !!   Calculates the overlap between two radial functions after
  !!   scaling by the appropriate power of r
  !!  INPUTS
  !!   
  !!  AUTHOR
  !!   N. Watanabe (Mizuho) with TM and DRB
  !!  CREATION DATE
  !!   2016
  !!  MODIFICATION HISTORY
  !!   2017/06/01 16:53 dave
  !!    Changed to remove result argument
  !!  SOURCE
  !!
  real(double) function overlapRadials( radial1, radial2, delta, l )
    use datatypes
    use numbers
    implicit none
    real(double), intent(in) :: radial1(:)
    real(double), intent(in) :: radial2(:)
    real(double), intent(in) :: delta
    integer, intent(in) :: l

    integer :: ir, nr
    real(double) :: r
    nr = size(radial1)
    
    overlapRadials = zero
    do ir=1, nr
       r = (ir-1)*delta
       overlapRadials = overlapRadials + delta*r*r*radial1(ir)*radial2(ir)
       !overlapRadials = overlapRadials + delta*r**(2*l+2)*radial1(ir)*radial2(ir)
    end do

    return
  end function overlapRadials
  !!***

! -----------------------------------------------------------
! Subroutine check_block
! -----------------------------------------------------------

!!****f* pseudo_tm_module/check_block *
!!
!!  NAME
!!   check_block
!!  USAGE
!!
!!  PURPOSE
!!   finds integration grids in a given block
!!   whose distances from the atom (xatom, yatom, zatom)
!!   are within a cutoff 'rcut'.
!! 
!!  INPUTS
!! (xblock, yblock, zblock): the position of the l.h.s.
!!                             of the block
!! (xatom, yatom, zatom)   : the position of the atom
!!             rcut          : cutoff 
!!
!!  OUTPUTS
!!    npoint  : # of integration grids whose distances 
!!              from the atom are within the cutoff rcut.
!!   for the following,  0 < ii < n_pts_in_block.
!!    ipoint(ii)  : index of found points ii,
!!    r_store(ii) : distance between the found point ii and the atom
!!    x_store(ii) : direction cosine x/r, for r=0, x/r = 0
!!    y_store(ii) : direction cosine y/r, for r=0, y/r = 0
!!    z_store(ii) : direction cosine z/r, for r=0, z/r = 0
!!  USES
!!
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   Jul. 2002
!!  MODIFICATION HISTORY
!!   2008/03/03 18:37 dave
!!    Removed dsqrt
!!  SOURCE
!!
  subroutine check_block &
       (xblock,yblock,zblock,xatom,yatom,zatom,rcut, &
       npoint, ip_store, r_store, x_store, y_store, z_store) 

    use numbers
    use global_module, only: rcellx,rcelly,rcellz
    use group_module,  only: blocks
    use block_module,  only: nx_in_block,ny_in_block,nz_in_block, &
         n_pts_in_block


    implicit none
    !Passed 
    real(double), intent(in):: xblock, yblock, zblock
    real(double), intent(in):: xatom, yatom, zatom, rcut
    integer, intent(out) :: npoint, ip_store(n_pts_in_block)
    real(double),intent(out) :: r_store(n_pts_in_block)
    real(double),intent(out) :: x_store(n_pts_in_block)
    real(double),intent(out) :: y_store(n_pts_in_block)
    real(double),intent(out) :: z_store(n_pts_in_block)
    !Local
    real(double):: dcellx_block,dcelly_block,dcellz_block
    real(double):: dcellx_grid, dcelly_grid, dcellz_grid
    real(double):: dx, dy, dz
    integer :: ipoint, iz, iy, ix
    real(double) ::  r2, r_from_i, rx, ry, rz, x, y, z, rcut2


    rcut2 = rcut* rcut
    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nx_in_block
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/ny_in_block
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nz_in_block

    ipoint=0
    npoint=0
    do iz=1,nz_in_block
       do iy=1,ny_in_block
          do ix=1,nx_in_block
             ipoint=ipoint+1

             dx=dcellx_grid*(ix-1)
             dy=dcelly_grid*(iy-1)
             dz=dcellz_grid*(iz-1)

             rx=xblock+dx-xatom
             ry=yblock+dy-yatom
             rz=zblock+dz-zatom
             r2 = rx * rx + ry * ry + rz * rz

             if(r2 < rcut2) then
                npoint=npoint+1
                r_from_i = sqrt( r2 )

                ! direction cosines needed for spher harms
                if ( r_from_i > RD_ERR ) then
                   x = rx / r_from_i
                   y = ry / r_from_i
                   z = rz / r_from_i
                else
                   x = zero ; y = zero ; z = zero
                   !!   r_from_i = RD_ERR !!   04/07/2002 TM
                end if
                ip_store(npoint)=ipoint
                r_store(npoint)=r_from_i
                x_store(npoint)=x
                y_store(npoint)=y
                z_store(npoint)=z
             endif  ! (r2 < rcut2) then

          enddo ! ix=1,nx_in_block
       enddo ! iy=1,ny_in_block
    enddo ! iz=1,nz_in_block
    return
  end subroutine check_block
!!***
end module pseudo_tm_module
