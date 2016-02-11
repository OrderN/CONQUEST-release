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
       nonloc_pp_derivative_tm, set_tm_pseudo, deallocate_pseudo_tm, loc_HF_stress, loc_G_stress
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
!!  SOURCE
!!
  subroutine init_pseudo_tm(ecore, ncf_out) 
    
    use datatypes,      only: double
    use global_module,  only: iprint_pseudo
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
! Subroutine set_tm_pseudo
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
!!  SOURCE
!!
  subroutine set_tm_pseudo

    use datatypes
    use numbers
    use global_module, only: rcellx,rcelly,rcellz,id_glob, ni_in_cell, &
                             iprint_pseudo, species_glob, nlpf, sf,    &
                             flag_basis_set, blips,                    &
                             IPRINT_TIME_THRES3, flag_analytic_blip_int
    use species_module, only: species, nlpf_species, n_species
    !  At present, these arrays are dummy arguments.
    use block_module, only : nx_in_block,ny_in_block,nz_in_block, &
         n_pts_in_block
    use group_module, only : blocks, parts
    use primary_module, only: domain
    use cover_module, only: DCS_parts
    use set_blipgrid_module, only : naba_atm
    use GenBlas, only: axpy, copy
    use GenComms, only: my_barrier, cq_abort, inode, ionode, myid
    use angular_coeff_routines, only : pp_elem
    use hartree_module, only: hartree
    use functions_on_grid, only: gridfunctions, pseudofns
    use dimens, only: n_my_grid_points
    use maxima_module, only: maxngrid
    use timer_module, only: cq_timer, start_timer, stop_print_timer, WITH_LEVEL

    implicit none 
    !local
    real(double):: dcellx_block,dcelly_block,dcellz_block
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom
    real(double):: xatom,yatom,zatom,alpha,step
    real(double):: xblock,yblock,zblock
    integer :: the_species
    integer :: j,iblock,the_l,ipoint, igrid
    real(double) :: r_from_i
    real(double) :: rr,a,b,c,d,x,y,z,nl_potential
    integer :: no_of_ib_ia, offset_position
    integer :: position,iatom
    integer :: stat, nl, npoint, ip
    integer :: i,m
    real(double) :: coulomb_energy
    real(double) :: rcut
    real(double) :: r1, r2, r3, r4, core_charge, gauss_charge
    real(double) :: val
    !allocatable
    real(double),allocatable :: chlocal_density(:), coulomb_potential(:)
    logical :: local_charge = .false.
    type(cq_timer) :: tmr_l_tmp1

    ! --  Start of subroutine  ---
    call start_timer(tmr_std_pseudopot)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
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
    !write(51,*) 'Pseudo: ',pseudo(the_species)%tm_loc_pot,pseudo(the_species)%vlocal%delta,pseudo(the_species)%vlocal%n
    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
       if(naba_atm(nlpf)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom=0
          do ipart=1,naba_atm(nlpf)%no_of_part(iblock)
             jpart=naba_atm(nlpf)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('set_ps: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atm(nlpf)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atm(nlpf)%list_atom(iatom,iblock)
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
                ! However, we can change this strategy by changing naba_atm(nlpf).
                !no_of_ib_ia = no_of_ib_ia +1

                xatom=DCS_parts%xcover(icover)
                yatom=DCS_parts%ycover(icover)
                zatom=DCS_parts%zcover(icover)

                !the_species=species(ig_atom)
                the_species = species_glob(ig_atom)

                !calculates distances between the atom and integration grid points
                !in the block and stores which integration grids are neighbours.
                rcut = core_radius(the_species) + very_small   !!   30072007 drb
                ! write(io_lun,*) ' rcut for check_block = ', rcut
                call check_block &
                     (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
                     npoint,ip_store,r_store,x_store,y_store,z_store) !out
                r_from_i = sqrt((xatom-xblock)**2+(yatom-yblock)**2+ &
                     (zatom-zblock)**2 )

                !Local part  ----
                ! construct the local part of the pseudopotential.
                ! local part of Siesta's pseudopotential is generated
                ! from the pseudo core charge (local charge) distribution
                ! 
                !write(51,*) 'Pseudo points: ',npoint
                if(npoint > 0) then
                   if(pseudo(the_species)%tm_loc_pot==loc_pot) then
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
                         !write(51,*) 'Pseudo dist: ',r_from_i, j
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

                            !write(51,*) 'Pseudo: ',igrid,r_from_i,j,a * r1 + b * r2 + c * r3 + d * r4
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

                !Projector Functions  ------------------------------------------
                ! pseudofunction(n_pts_in_block,ncf,naba_atm,iblock) 
                !  ncf is fixed for all atoms in this version, 
                !  though it is not necessary.  
                !  Now, I keep the strategy of the old version  for ncf, 
                !  but we can pack pseudofunction with respect to naba_atm 
                !  by using species-dependent core radius.
                !  (needs to be checked!  01/07/2002 TM)

                if(flag_basis_set==blips.AND.npoint > 0.AND.(.NOT.flag_analytic_blip_int)) then
                   do nl= 1, pseudo(the_species)%n_pjnl

                      the_l=pseudo(the_species)%pjnl_l(nl)
                      !offset_position = (no_of_ib_ia-1) * ncf * n_pts_in_block + &
                      !     offset_mcomp(nl, the_species) * n_pts_in_block
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
                             
                            !RC putting in call for alternative Y_lm's here
                            if(flag_angular_new) then
                               if (the_l>0) then
                                  i = 0
                                  do m = -the_l, the_l
                                     call pp_elem(nl_potential,the_l,m,x,y,z,r_from_i,val)
                                     gridfunctions(pseudofns)%griddata(position+i*n_pts_in_block) = val
                                     !write(io_lun,*) val, 'val', the_l, m,myid
                                     i = i+1
                                  enddo
                               else
                                  call pp_elem(nl_potential,the_l,0,x,y,z,r_from_i,val)
                                  gridfunctions(pseudofns)%griddata(position) = val
                                  !write(io_lun,*) val, 'val', the_l, m,myid
                               endif
                            else
                               call cq_abort("Must have FlagNewAngular T for siesta/abinit pseudos")
                               
                               !%%! ! S, P or D components
                               !%%! if ( the_l .eq. 0 ) then      ! NonLocal : S component
                               !%%!    gridfunctions(pseudofns)%griddata(position) =   &
                               !%%!         nl_potential*spherical_harmonic_norm(1)
                               !%%!    !write(io_lun,*) pseudofunctions(position), 'val_orig',myid
                               !%%!    
                               !%%! else if ( the_l .eq. 1 ) then ! Nonlocal : P components
                               !%%!    gridfunctions(pseudofns)%griddata(position) =  & 
                               !%%!         nl_potential*spherical_harmonic_norm(2) * x 
                               !%%!    gridfunctions(pseudofns)%griddata(position+1*n_pts_in_block) =  &
                               !%%!         nl_potential*spherical_harmonic_norm(3) * y 
                               !%%!    gridfunctions(pseudofns)%griddata(position+2*n_pts_in_block) =   &
                               !%%!         nl_potential*spherical_harmonic_norm(4) * z
                               !%%!    !write(io_lun,*) pseudofunctions(position),the_l, -1,myid 
                               !%%!    !write(io_lun,*) pseudofunctions(position+1*n_pts_in_block),the_l,0,myid
                               !%%!    !write(io_lun,*) pseudofunctions(position+2*n_pts_in_block),the_l,+1,myid
                               !%%!    
                               !%%! else if ( the_l .eq. 2 ) then ! Nonlocal : D components
                               !%%!    gridfunctions(pseudofns)%griddata(position) =  &
                               !%%!         nl_potential*spherical_harmonic_norm(5) * (x*x - y*y) 
                               !%%!    gridfunctions(pseudofns)%griddata(position+1*n_pts_in_block) = &
                               !%%!         nl_potential*spherical_harmonic_norm(6) * (three*z*z-one)
                               !%%!    gridfunctions(pseudofns)%griddata(position+2*n_pts_in_block) = &
                               !%%!         nl_potential*spherical_harmonic_norm(7) * x * y 
                               !%%!    gridfunctions(pseudofns)%griddata(position+3*n_pts_in_block) = &
                               !%%!      nl_potential*spherical_harmonic_norm(8) * x * z 
                               !%%!    gridfunctions(pseudofns)%griddata(position+4*n_pts_in_block) = &
                               !%%!         nl_potential*spherical_harmonic_norm(9) * y * z 
                               !%%!    
                               !%%! else
                               !%%!    write(io_lun,fmt='(10x," ERROR in set_tm_pseudo! l_component = ",i3)') the_l
                               !%%! end if
                            endif ! flag_angular_new
                         endif !(j+1 < pseudo(the_species)%pjnl(nl)%n) then
                      enddo ! ip=1,npoint
                   enddo ! nl= 1, pseudo(the_species)%n_pjnl
                endif! (npoint > 0) then
                no_of_ib_ia = no_of_ib_ia + nlpf_species(the_species)*n_pts_in_block
             enddo ! naba_atoms
          enddo ! naba_part
       endif !(naba_atm(nlpf)%no_of_part(iblock) > 0) !naba atoms?
    enddo ! iblock : primary set of blocks
    !RC putting in stop for testing
    !stop
    

    ! now we must use FFT to transform the core charge density into
    ! reciprocal space to construct the pseudopotential in 
    ! reciprocal space. This is then transformed back into real space 
    ! and returned in pseudopotential from subroutine hartree()
    call my_barrier()

    ! local potential 
    if(local_charge) then
       call hartree( chlocal_density, pseudopotential, maxngrid, coulomb_energy )
    else
       !write(io_lun,*) 'Calling hartree to do FFT of loc pot'
       call hartree( chlocal_density, coulomb_potential, maxngrid, coulomb_energy )
       call axpy( n_my_grid_points, -one, coulomb_potential, 1, &
            pseudopotential, 1 )
       call start_timer(tmr_std_allocation)
       deallocate(coulomb_potential, STAT=stat)
       if(stat/=0) call cq_abort('set_ps: dealloc c_p ', stat)
       call stop_timer(tmr_std_allocation)
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
!!   2016/01/29 14:25 dave
!!    Bug fix for local G stress: accumulate loc_charge (was just storing !)
!!  SOURCE
!!
  subroutine loc_pp_derivative_tm ( hf_force, density, size )

    use datatypes
    use numbers
    use dimens, only: grid_point_volume, n_my_grid_points
    use global_module, only: rcellx,rcelly,rcellz,id_glob, iprint_pseudo, species_glob, nlpf,ni_in_cell, sf
    use block_module, only : n_pts_in_block
    use group_module, only : blocks, parts
    use primary_module, only: domain
    use cover_module, only: DCS_parts
    use set_blipgrid_module, only : naba_atm

    use species_module, only: species
    use GenComms, only: gsum, cq_abort, inode, ionode
    use hartree_module, only: hartree, Hartree_stress
    use maxima_module, only: maxngrid

    implicit none   

    ! Passed variables
    integer :: size
    real(double),intent(in)  :: density( size )
    real(double),intent(out) :: hf_force( 3, ni_in_cell )

    ! Local variables
    integer :: i, j 

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
    real(double),allocatable :: h_potential(:)
    ! Store potential scaled by 2GaGb/G^4 for stress
    real(double),allocatable :: loc_charge(:)

    call start_timer(tmr_std_pseudopot)
    if(iprint_pseudo>2.AND.inode==ionode) write(io_lun,fmt='(4x,"Doing TM force with pseudotype: ",i3)') pseudo(1)%tm_loc_pot
    ! the structure of this subroutine is similar to set_tm_pseudo et.
    HF_force = 0
    loc_HF_stress = zero
    loc_G_stress = zero

    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz

    ! get Hartree potential
    call start_timer(tmr_std_allocation)
    allocate(h_potential(maxngrid),STAT=stat)
    if(stat /= 0) call cq_abort &
         ('loc_pp_derivative_tm: alloc h_potential',stat)
    allocate(loc_charge(maxngrid),STAT=stat)
    if(stat /= 0) call cq_abort &
         ('loc_pp_derivative_tm: alloc loc_charge',stat)
    call stop_timer(tmr_std_allocation)
    ! Save existing Hartree stress
    temp = Hartree_stress
    call hartree( density, h_potential, maxngrid, h_energy )
    Hartree_stress = temp
    loc_charge = zero

    ! now loop over grid points and accumulate HF part of the force
    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
       if(naba_atm(nlpf)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom=0
          do ipart=1,naba_atm(nlpf)%no_of_part(iblock)
             jpart=naba_atm(nlpf)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then
                call cq_abort('set_ps: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atm(nlpf)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atm(nlpf)%list_atom(iatom,iblock)
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
                rcut = core_radius(the_species) + very_small   !!   03032003TM
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
                if(pseudo(the_species)%tm_loc_pot==loc_pot) then
                   step = pseudo(the_species)%vlocal%delta
                else
                   step = pseudo(the_species)%chlocal%delta
                end if
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
                      j = aint( r_from_i / step ) + 1
                      !As we use the maximum of cutoff for check_block
                      ! overrun should occur in some cases.
                      !if(j > N_TAB-1) call cq_abort &
                      !     ('set_ps: overrun problem3',j)
                      ! check j (j+1 =< N_TAB)
                      ! Use the spline interpolation tables
                      if(pseudo(the_species)%tm_loc_pot==loc_pot) then
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
                            loc_charge(igrid) = loc_charge(igrid) + a*r1+b*r2+c*r3+d*r4
                         end if ! j+1<pseudo(the_species)%chlocal%n
                      end if
                   enddo ! ip=1, npoint
                endif  !(npoint > 0) then

             enddo ! ia=1,naba_atm(nlpf)%no_atom_on_part(ipart,iblock)
          enddo ! ipart=1,naba_atm(nlpf)%no_of_part(iblock)
       endif    ! (naba_atm(nlpf)%no_of_part(iblock) > 0) then ! if there are naba atoms
    enddo ! iblock = 1, domain%groups_on_node ! primary set of blocks
    ! Save existing Hartree stress
    temp = Hartree_stress
    call hartree( density, h_potential, maxngrid, h_energy, loc_charge,loc_G_stress )
    Hartree_stress = temp
    if(pseudo(the_species)%tm_loc_pot==loc_pot) loc_G_stress = -loc_G_stress

    ! and add contributions from all nodes
    !  In the future this should be replaced by summation with local communication
    !    Tsuyoshi Miyazaki
    call gsum(HF_force,3,ni_in_cell)
    call gsum(loc_HF_stress,3)
    ! Don't gsum loc_G_stress - that's done in hartree
    ! Deallocate added by TM, 2005/08/11
    call start_timer(tmr_std_allocation)
    deallocate(loc_charge,STAT=stat)
    if(stat /= 0) call cq_abort &
         ('loc_pp_derivative_tm: dealloc loc_charge',stat)
    deallocate(h_potential,STAT=stat)
    if(stat /= 0) call cq_abort &
         ('loc_pp_derivative_tm: dealloc h_potential',stat)
    call stop_timer(tmr_std_allocation)
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
!!  
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
    use set_blipgrid_module, only : naba_atm
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
       if(naba_atm(nlpf)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom=0
          do ipart=1,naba_atm(nlpf)%no_of_part(iblock)
             jpart=naba_atm(nlpf)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('ps_derivative: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atm(nlpf)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atm(nlpf)%list_atom(iatom,iblock)
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
       endif !(naba_atm(nlpf)%no_of_part(iblock) > 0) !naba atoms?
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
      use ewald_module, only: erfc

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
               vlocal = pseudo%vlocal%f(imesh_rr)
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
