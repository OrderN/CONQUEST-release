! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module density_module
! ------------------------------------------------------------------------------
! Code area 5: charge density and self-consistency
! ------------------------------------------------------------------------------

!!****h* Conquest/density_module *
!!  NAME
!!   density_module
!!  PURPOSE
!!   Density initialization routines
!!  AUTHOR
!!   L.K.Dash
!!  CREATION DATE
!!   02/05/02
!!  MODIFICATION HISTORY
!!   17/07/02 lkd
!!    moved declaration of density table variables to derived types in
!!    density_read_module.f90
!!   14/08/2002 drb and mjg
!!    changed naba table to supp from proj and added header
!!   15:21, 25/09/2002 mjg & drb
!!    Tidied and updated to use dens tables from set_blipgrid, and
!!    changed to linear interpolation (for now)
!!   15:25, 27/02/2003 drb & tm
!!    Added flag_no_atomic_densities
!!   2006/03/04 09:10 dave
!!    Added get_electronic_density
!!   2008/02/04 17:13 dave
!!    Changdes for output to file not stdout
!!   2008/05/23 ast
!!    Added timers
!!   2009/01/23 14:16 dave
!!    Added Becke weights and Becke atomic charges
!!   2011/03/30 19:15 M.Arita
!!    Added the charge density for core electrons
!!   2011/04/13 L.Tong
!!    Added spin polarisation
!!    Added density_dn module variable for electron density for spin_down component
!!   2011/09/16 11:16 dave
!!    Update and bug fix for atomic radii in Becke weights
!!   2011/07/20 12:00 dave
!!    Changes for cDFT: Becke weight matrix
!!   2011/10/06 13:57 dave
!!    Added cDFT routines
!!   2011/11/25 L.Tong
!!    Major change to spin implementation, now density refers to total
!!    density, which is NOT allocated for spin polarised calculations
!!    as default, and density_up and density_dn are arrays for spin up
!!    and down components, which are allocated at the start for spin
!!    polarised calculations. For spin non-polarised calculations,
!!    only density is allocated at the start
!!   2012/04/21 L.Tong
!!   - Moved get_cdft_constraint to cdft_module, since it actually
!!     does not use variables in density module. And this removes
!!     energy module dependence here.
!!   2012/05/29 L.Tong
!!   - Added subroutine for calculating electron numbers
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2015/11/09 08:34 dave with TM and NW from Mizuho
!!    Changes to implement neutral atom potential - mainly the set_density_atom routine
!!   2016/01/12 08:26 dave
!!    Merging set_density_atom and set_density and renaming to set_atomic_density
!!  SOURCE
module density_module

  use datatypes
  use global_module,          only: io_lun, iprint_SC
  use timer_module,           only: start_timer, stop_timer, cq_timer, stop_print_timer
  use timer_module,           only: start_backtrace, stop_backtrace
  use timer_stdclocks_module, only: tmr_std_chargescf

  implicit none
  save

  logical :: flag_no_atomic_densities
  ! total charge density
  ! values at gridpoints, to be calculated. second dim is spin
  real(double), allocatable, dimension(:,:) :: density
  ! for Neutral atom potential
  real(double), allocatable, dimension(:)   :: density_atom
  ! values at gridpoints, to be calculated for P.C.C.
  real(double), allocatable, dimension(:)   :: density_pcc
  ! Scaling factor for atomic densities, array of spin
  real(double), allocatable, dimension(:)   :: density_scale

  ! Becke weights
  real(double), allocatable, dimension(:)   :: bw
  real(double), allocatable, dimension(:,:) :: bwgrid
  ! total atomcharge, second dim is spin
  real(double), allocatable, dimension(:,:) :: atomcharge

  logical :: weights_set = .false.
  real(double), dimension(95) :: atrad

  ! for spin polarised calculations
  ! ! density for spin up component
  ! real(double) :: density_scale_up
  ! real(double), allocatable, dimension(:) :: density_up
  ! real(double), allocatable, dimension(:) :: atomcharge_up
  ! ! density for spin down component
  ! real(double) :: density_scale_dn
  ! real(double), allocatable, dimension(:) :: density_dn
  ! real(double), allocatable, dimension(:) :: atomcharge_dn

  ! Area identification
  integer, parameter, private :: area = 5

  ! RCS tag for object file identification
  character(len=80), private :: &
       RCSid = "$Id$"

!!***

contains

  ! -----------------------------------------------------------
  ! Subroutine set_atomic_density
  ! -----------------------------------------------------------

  !!****f* density_module/set_atomic_density *
  !!
  !!  NAME
  !!   set_density
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Puts the atomic charge density on the grid points belonging
  !!   to this processor.  Based largely on the set_pseudopotential
  !!   routine from pseudopotential.module.f90
  !!  INPUTS
  !! (none)
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   L.K.Dash
  !!  CREATION DATE
  !!   02.05.02
  !!  MODIFICATION HISTORY
  !!   14/08/2002 drb and mjg
  !!    Changed naba table to be use table based on support region radius
  !!   08:14, 2003/03/19 dave
  !!    Tidied and re-instated splint call
  !!   11:50, 30/09/2003 drb
  !!    Added entry print statement
  !!   2004/10/29 drb
  !!    Bug fix for species
  !!   2007/05/08 17:00 dave
  !!    Added scaling to fix electron number (specifically to ensure
  !!    correct charge in cell)
  !!   2008/05/23 ast
  !!    Added timers
  !!   2011/05/16 L. Tong
  !!    - Added spin polarisation
  !!    - For fixed spin case the densities are renormalised to give the
  !!      spin populations, which sums to give total electrons in cell
  !!    - For variable spin case the densities are renormalised to give the
  !!      total electrons in cell
  !!   2012/03/13 L.Tong
  !!    - Rewriting spin implementation. Now the spin non-polarised
  !!      case are treated the same way as fixed spin case.
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2016/01/12 08:26 dave
  !!    - Changed name and included VNA possibility; also added temporary array to store density
  !!   2016/07/20 16:30 nakata
  !!    Renamed naba_atm -> naba_atoms_of_blocks
  !!  SOURCE
  !!
  subroutine set_atomic_density(flag_set_density,level)

    use datatypes
    use numbers
    use global_module,       only: rcellx, rcelly, rcellz, id_glob, &
                                   ni_in_cell, iprint_SC,           &
                                   species_glob, dens, ne_in_cell,  &
                                   ne_spin_in_cell,                 &
                                   IPRINT_TIME_THRES3, nspin,       &
                                   spin_factor,                     &
                                   flag_fix_spin_population, &
                                   flag_neutral_atom
    use block_module,        only: nx_in_block, ny_in_block,        &
                                   nz_in_block, n_pts_in_block
    use group_module,        only: blocks, parts
    use primary_module,      only: domain
    use cover_module,        only: DCS_parts
    use set_blipgrid_module, only: naba_atoms_of_blocks
    use GenComms,            only: my_barrier, cq_abort, inode, ionode, gsum
    use atomic_density,      only: atomic_density_table
    use spline_module,       only: splint
    use dimens,              only: n_my_grid_points, grid_point_volume
    use GenBlas,             only: rsum, scal
    use timer_module,        only: WITH_LEVEL
    use maxima_module,  only: maxngrid

    implicit none

    ! Passed Variables
    integer, optional :: level
    logical :: flag_set_density

    ! Local Variables
    integer :: ipart, jpart, ind_part, ia, ii, icover, ig_atom
    integer :: the_species
    integer :: ix, iy, iz, j, iblock, ipoint, igrid
    integer :: no_of_ib_ia, offset_position
    integer :: position, iatom, icheck, spin

    real(double) :: dcellx_block, dcelly_block, dcellz_block
    real(double) :: dcellx_grid, dcelly_grid, dcellz_grid
    real(double) :: xatom, yatom, zatom, step, loc_cutoff
    real(double) :: xblock, yblock, zblock, alpha
    real(double) :: dx, dy, dz, rx, ry, rz, r2, r_from_i
    ! local charge density returned from splint routine
    real(double), dimension(:), allocatable :: store_density
    real(double)   :: local_density, scale
    type(cq_timer) :: tmr_l_tmp1
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level, stat

    ! logical flag to warn if splint routine called out of the
    ! tabulated range.
    logical :: range_flag

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='set_density', &
         where=area,level=backtrace_level)
!****lat>$
    allocate(store_density(maxngrid), STAT=stat)
    if (stat/=0) call cq_abort('Error allocating store_density: ', maxngrid)

    if (inode == ionode .and. iprint_SC >= 2) &
         write (io_lun, fmt='(2x,"Entering set_density")')

    call start_timer(tmr_std_chargescf)
    call start_timer(tmr_l_tmp1, WITH_LEVEL)

    !call sub_enter_output('set_density',2,'5')

    ! initialize density
    if(flag_set_density) then
       density = zero
    end if
    store_density = zero

    ! determine the block and grid spacing
    dcellx_block = rcellx / blocks%ngcellx
    dcelly_block = rcelly / blocks%ngcelly
    dcellz_block = rcellz / blocks%ngcellz
    dcellx_grid = dcellx_block / nx_in_block
    dcelly_grid = dcelly_block / ny_in_block
    dcellz_grid = dcellz_block / nz_in_block

    ! loop around grid points in this domain, and for each
    ! point, get contributions to the charge density from atoms which are
    ! within the cutoff distance to that grid point

    call my_barrier()

    do iblock = 1, domain%groups_on_node ! loop over blocks of grid points
       !write(io_lun,*) 'Block ',iblock
       ! determine the position of this block
       xblock = (domain%idisp_primx(iblock) + domain%nx_origin - 1) * dcellx_block
       yblock = (domain%idisp_primy(iblock) + domain%ny_origin - 1) * dcelly_block
       zblock = (domain%idisp_primz(iblock) + domain%nz_origin - 1) * dcellz_block
       ! if there are neighbour partitions
       if (naba_atoms_of_blocks(dens)%no_of_part(iblock) > 0) then
          iatom = 0
          ! loop over neighbour partitions of this block
          do ipart = 1, naba_atoms_of_blocks(dens)%no_of_part(iblock)
             jpart = naba_atoms_of_blocks(dens)%list_part(ipart,iblock)
             if (jpart > DCS_parts%mx_gcover) then
                call cq_abort('set_ps: JPART ERROR ', ipart, jpart)
             endif
             ind_part = DCS_parts%lab_cell(jpart)
             ! .... then atoms in this partition.
             do ia = 1, naba_atoms_of_blocks(dens)%no_atom_on_part(ipart, iblock)
                iatom = iatom + 1
                ii = naba_atoms_of_blocks(dens)%list_atom(iatom, iblock)
                icover = DCS_parts%icover_ibeg(jpart) + ii - 1
                ig_atom = id_glob(parts%icell_beg(ind_part) + ii - 1)

                if (parts%icell_beg(ind_part) + ii - 1 > ni_in_cell) then
                   call cq_abort('set_ps: globID ERROR ', &
                                 ii, parts%icell_beg(ind_part))
                endif
                if (icover > DCS_parts%mx_mcover) then
                   call cq_abort('set_ps: icover ERROR ', &
                                 icover, DCS_parts%mx_mcover)
                endif

                ! determine the position of the current atom
                xatom = DCS_parts%xcover(icover)
                yatom = DCS_parts%ycover(icover)
                zatom = DCS_parts%zcover(icover)

                ! determine the type of the current atom
                the_species = species_glob(ig_atom)
                loc_cutoff = atomic_density_table(the_species)%cutoff
                ! step in the density table
                step = loc_cutoff / &
                       real(atomic_density_table(the_species)%length - 1, &
                            double)
                icheck = 0
                ipoint = 0
                ! loop over the grid points in the current block
                do iz = 1, nz_in_block
                   do iy = 1, ny_in_block
                      do ix = 1, nx_in_block
                         ipoint = ipoint + 1
                         igrid = n_pts_in_block * (iblock - 1) + ipoint
                         ! position= offset_position + ipoint
                         if (igrid > n_my_grid_points) &
                              call cq_abort('set_density: igrid error ', &
                                            igrid, n_my_grid_points)
                         dx = dcellx_grid * (ix - 1)
                         dy = dcelly_grid * (iy - 1)
                         dz = dcellz_grid * (iz - 1)
                         ! determine separation between the current
                         ! grid point and atom
                         rx = xblock + dx - xatom
                         ry = yblock + dy - yatom
                         rz = zblock + dz - zatom
                         r2 = rx * rx + ry * ry + rz * rz
                         if (r2 < loc_cutoff * loc_cutoff) then
                            r_from_i = sqrt(r2)
                            call splint(step, &
                                 atomic_density_table(the_species)%table(:),&
                                 atomic_density_table(the_species)%d2_table(:),&
                                 atomic_density_table(the_species)%length,&
                                 r_from_i, local_density, range_flag)
                            if (range_flag) &
                                 call cq_abort('set_density: overrun problem')
                            ! Store the density for this grid point
                            store_density(igrid) = store_density(igrid) + local_density
                         end if ! if this point is within cutoff
                      end do !ix gridpoints
                   end do  !iy gridpoints
                end do   !iz gridpoints
             end do ! end loop over atoms in this partition
          end do ! end loop over partitions
       end if ! end if there are neighbour atoms
    end do ! end loop over blocks

    ! renormalise density
    local_density = zero
    do iz = 1, n_my_grid_points
       local_density = local_density + store_density(iz)
    end do
    local_density = local_density * grid_point_volume
    call gsum(local_density)
    scale = ne_in_cell/local_density
    store_density(:) = store_density(:) * scale
    
    if(flag_neutral_atom) density_atom(:) = store_density(:)
    if(flag_set_density) then
       ! first assume atomic density is evenly
       ! divided among the spin channels, this
       ! applies also to spin non-polarised
       ! calculations
       local_density = half*local_density
       do spin = 1, nspin
          density_scale(spin) = ne_spin_in_cell(spin) / local_density
          density(:,spin) = density_scale(spin) * half* store_density(:) ! Check later
          if (inode == ionode .and. iprint_SC > 0) &
               write (io_lun, &
               fmt='(10x,"In set_density, electrons (spin=",i1,"): ",f20.12)') &
               spin, density_scale(spin) * local_density
       end do
    end if

    call my_barrier()
    deallocate(store_density, STAT=stat)
    if (stat/=0) call cq_abort('Error deallocating store_density: ', maxngrid)

    call stop_print_timer(tmr_l_tmp1, "set_density", IPRINT_TIME_THRES3)
    call stop_timer(tmr_std_chargescf)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='set_density')
!****lat>$

    return
  end subroutine set_atomic_density
  !!***
  
  ! -----------------------------------------------------------
  ! Subroutine set_density_pcc
  ! -----------------------------------------------------------

  !!****f* density_module/set_densityi_pcc *
  !!
  !!  NAME
  !!   set_density_pcc
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Puts the P.C.C. charge density on the grid points belonging
  !!   to this processor.  Based largely on the set_pseudopotential
  !!   routine from pseudopotential.module.f90
  !!   This subroutine is based upon sbrt: set_density.
  !!  INPUTS
  !! (none)
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   M.Arita
  !!  CREATION DATE
  !!   2011/03/30
  !!  MODIFICATION HISTORY
  !!   2016/07/20 16:30 nakata
  !!    Renamed naba_atm -> naba_atoms_of_blocks
  !!  SOURCE
  !!
  subroutine set_density_pcc()

    use datatypes
    use numbers,             only: zero, one
    use global_module,       only: rcellx, rcelly, rcellz, id_glob, &
                                   ni_in_cell, iprint_SC,           &
                                   species_glob, dens, ne_in_cell,  &
                                   IPRINT_TIME_THRES3
    use block_module,        only: nx_in_block, ny_in_block,        &
                                   nz_in_block, n_pts_in_block
    use group_module,        only: blocks, parts
    use primary_module,      only: domain
    use cover_module,        only: DCS_parts
    use set_blipgrid_module, only: naba_atoms_of_blocks
    use GenComms,            only: my_barrier, cq_abort, inode, ionode, gsum
    use pseudo_tm_info,      only: pseudo
    use spline_module,       only: splint
    use dimens,              only: n_my_grid_points, grid_point_volume
    use GenBlas,             only: rsum, scal
    use timer_module

    implicit none

    ! Local Variables
    integer :: ipart, jpart, ind_part, ia, ii, icover, ig_atom
    integer :: the_species
    integer :: ix, iy, iz, j, iblock, ipoint, igrid
    integer :: no_of_ib_ia, offset_position
    integer :: position, iatom, icheck

    real(double) :: dcellx_block, dcelly_block, dcellz_block
    real(double) :: dcellx_grid, dcelly_grid, dcellz_grid
    real(double) :: xatom, yatom, zatom
    real(double) :: pcc_cutoff, pcc_step
    real(double) :: xblock, yblock, zblock,  alpha
    real(double) :: dx, dy, dz, rx, ry, rz, r2, r_from_i
    real(double) :: pcc_density !P.C.C. charge density returned from splint routine
    type(cq_timer) :: tmr_l_tmp1
    ! logical flag to warn if splint routine called out of the tabulated range.
    logical :: range_flag

    if (inode == ionode .and. iprint_SC >= 2) &
         write (io_lun, fmt='(2x,"Entering set_density_pcc")')

    call start_timer(tmr_std_chargescf)
    call start_timer(tmr_l_tmp1, WITH_LEVEL)
    density_pcc = zero  ! initialize density
    ! write(io_lun,*) 'Size of density: ',size(density)
    ! call scal(n_my_grid_points,zero,density,1)

    ! determine the block and grid spacing
    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nx_in_block
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/ny_in_block
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nz_in_block

    ! loop around grid points in this domain, and for each
    ! point, get contributions to the charge density from atoms which are
    ! within the cutoff distance to that grid point

    call my_barrier()

    do iblock = 1, domain%groups_on_node ! loop over blocks of grid points
       !write(io_lun,*) 'Block ',iblock
       ! determine the position of this block
       xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
       if(naba_atoms_of_blocks(dens)%no_of_part(iblock) > 0) then ! if there are neighbour partitions
          iatom=0
          ! loop over neighbour partitions of this block
          do ipart=1,naba_atoms_of_blocks(dens)%no_of_part(iblock)   !for now, dens is used even for P.C.C.
             jpart=naba_atoms_of_blocks(dens)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then
                call cq_abort('set_ps: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             ! .... then atoms in this partition.
             do ia=1,naba_atoms_of_blocks(dens)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atoms_of_blocks(dens)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)
                ! determine the type of the current atom
                the_species=species_glob(ig_atom)
                ! for P.C.C.
                if (.NOT. pseudo(the_species)%flag_pcc) then
                  cycle
                endif

                if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                   call cq_abort('set_ps: globID ERROR ', &
                        ii,parts%icell_beg(ind_part))
                endif
                if(icover > DCS_parts%mx_mcover) then
                   call cq_abort('set_ps: icover ERROR ', &
                        icover,DCS_parts%mx_mcover)
                endif

                ! determine the position of the current atom
                xatom=DCS_parts%xcover(icover)
                yatom=DCS_parts%ycover(icover)
                zatom=DCS_parts%zcover(icover)

                pcc_cutoff = pseudo(the_species)%chpcc%cutoff
                ! step in the density table
                pcc_step = pseudo(the_species)%chpcc%delta
                icheck=0
                ipoint=0
                ! loop over the grid points in the current block
                do iz=1,nz_in_block
                   do iy=1,ny_in_block
                      do ix=1,nx_in_block
                         ipoint=ipoint+1
                         igrid=n_pts_in_block*(iblock-1)+ipoint
                         ! position= offset_position + ipoint
                         if(igrid > n_my_grid_points) &
                              call cq_abort('set_density: igrid error ', &
                              igrid, n_my_grid_points)
                         dx=dcellx_grid*(ix-1)
                         dy=dcelly_grid*(iy-1)
                         dz=dcellz_grid*(iz-1)
                         ! determine separation between the current
                         ! grid point and atom
                         rx=xblock+dx-xatom
                         ry=yblock+dy-yatom
                         rz=zblock+dz-zatom
                         r2 = rx * rx + ry * ry + rz * rz
                         if(r2 < pcc_cutoff * pcc_cutoff) then
                            r_from_i = sqrt(r2)
                            ! j = aint(r_from_i/step) + 1
                            ! ! check j (j+1 =< N_TAB_MAX)
                            ! if(j+1 > atomic_density_table(the_species)%&
                            !      length) then
                            !    call cq_abort('set_density: overrun &
                            !         &problem',j)
                            ! endif
                            ! ! Linear interpolation for now
                            ! alpha = (r_from_i/step)-real(j-1,double)
                            ! local_density = (one-alpha)*&
                            !      atomic_density_table(the_species)%&
                            !      table(j) + alpha*&
                            !      atomic_density_table(the_species)%&
                            !      table(j+1)
                            call splint( pcc_step,pseudo(the_species)%&
                                 chpcc%f(:), pseudo(the_species)%&
                                 chpcc%d2(:), pseudo(the_species)%&
                                 chpcc%n, r_from_i,pcc_density,&
                                 range_flag)
                            if(range_flag) &
                                 call cq_abort('set_density_pcc: overrun problem')
                            ! recalculate the density for this grid point
                            density_pcc(igrid) = density_pcc(igrid) + pcc_density
                         endif ! if this point is within cutoff
                         ! test output of densities
                         ! print '(4(f10.6,a))', xblock+dx, " *", &
                         !      yblock+dy, " *", zblock+dz, " *", &
                         !      density(igrid), " *"
                      enddo !ix gridpoints
                   enddo  !iy gridpoints
                enddo   !iz gridpoints
             enddo ! end loop over atoms in this partition
          enddo ! end loop over partitions
       endif ! end if there are neighbour atoms
    enddo ! end loop over blocks
    ! local_density = zero
    ! do iz=1,n_my_grid_points
    !   local_density = local_density + density(iz)
    ! end do
    ! local_density = local_density*grid_point_volume
    ! local_density = grid_point_volume * rsum( n_my_grid_points, density, 1 )

    ! call gsum(local_density)
    ! Correct electron density
    ! density_scale = ne_in_cell/local_density
    ! density = density_scale*density
    ! if(inode.eq.ionode.AND.iprint_SC>0) &
    !      write(io_lun,fmt='(10x,"In set_density, electrons: ",f20.12)') &
    !      density_scale*local_density
    call my_barrier()
    call stop_print_timer(tmr_l_tmp1,"set_density",IPRINT_TIME_THRES3)
    call stop_timer(tmr_std_chargescf)
    return
  end subroutine set_density_pcc
  !!***

  ! -----------------------------------------------------------
  ! Subroutine get_electronic_density
  ! -----------------------------------------------------------

  !!****f* density_module/get_electronic_density *
  !!
  !!  NAME
  !!   get_electronic_density
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Gets electronic density on grid
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   E.H.Hernandez
  !!  CREATION DATE
  !!   16/04/96
  !!  MODIFICATION HISTORY
  !!   14/08/2000 TM
  !!    Added calls for new act_on_vector
  !!   15/05/2001 dave
  !!    Converted to F90 and ROBODoc and reduced passed variables
  !!   11/06/2001 dave
  !!    Added RCS Id and Log tags and GenComms dependence
  !!   15:00, 12/03/2003 drb
  !!    Removed reference to workspace module
  !!   14:40, 29/07/2003 drb
  !!    Changed density = 0 to density = 0.0_double
  !!   13:53, 2003/09/22 dave
  !!    Removed references to old variables MAX_?_ELEMENTS
  !!   11:57, 30/09/2003 drb
  !!    Added print statement on entry
  !!   10:09, 13/02/2006 drb
  !!    Removed all explicit references to data_ variables and rewrote
  !!    in terms of new matrix routines
  !!   2006/03/04 09:08 dave
  !!    Changed to use new grid ideas, put into density_module
  !!   2008/05/23 ast
  !!    Added timers
  !!   2011/05/20 L.Tong
  !!    Added Spin polarisation
  !!   Tuesday, 2011/08/02 L.Tong
  !!    Corrected the bug that spin == 2 is evaluated at the same time
  !!    as test for present (spin), this will give a segmentation fault
  !!    if spin is not present with some compilers. Moved the spin == 2
  !!    conditional inside the present (spin) conditional
  !!   2011/08/25 L.Tong
  !!    Removed the dependence on global module variable
  !!    flag_spin_polarisation, and uses present (spin) to determine if
  !!    we are in spin polarisation mode.
  !!    - if spin is present, then denout gives either density_up or
  !!      density_dn, and electrons gives either electrons_up and
  !!      electrons_dn
  !!   2012/03/18 L.Tong
  !!    - removed optional parameter ~spin~
  !!    - both ~spin~ components are done in the subroutine, using loop
  !!      over spin index.
  !!    - There is no longer a factor of two difference for spin
  !!      non-polarised case. Spin non-polarised case will just store the
  !!      spin up component in density(1). This is for consistency.
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2016/07/15 18:30 nakata
  !!    Renamed sf_H_sf_rem -> atomf_H_atomf_rem
  !!   2016/07/20 16:30 nakata
  !!    Renamed naba_atm -> naba_atoms_of_blocks
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf instead of sf
  !!   2016/08/09 14:00 nakata
  !!    Renamed support -> atom_fns
  !!   2016/08/09 18:00 nakata
  !!    Introduced matKatomf
  !!  SOURCE
  !!
  subroutine get_electronic_density(denout, electrons, atom_fns, &
                                    atom_fns_K, inode, ionode, size, level)

    use datatypes
    use numbers
    use GenBlas,                     only: scal, rsum
    use mult_module,                 only: matKatomf   ! nakata3
    use dimens,                      only: n_my_grid_points, grid_point_volume
    use block_module,                only: n_pts_in_block
    use set_bucket_module,           only: rem_bucket, atomf_H_atomf_rem
    use calc_matrix_elements_module, only: act_on_vectors_new
    use primary_module,              only: domain
    use set_blipgrid_module,         only: naba_atoms_of_blocks
    use GenComms,                    only: gsum
    use global_module,               only: iprint_SC, ni_in_cell, &
                                           flag_Becke_weights, nspin, &
                                           spin_factor, &
                                           sf, paof, atomf   ! nakata3
    use functions_on_grid,           only: gridfunctions, fn_on_grid

    implicit none

    ! Passed variables
    real(double), dimension(:)   :: electrons
    real(double), dimension(:,:) :: denout
    integer :: inode, ionode, size
    integer :: atom_fns, atom_fns_K
    integer, optional :: level

    ! Local variables
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level
    integer        :: blk, i_count_alpha, n, n_i, n_point, spin

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='get_electronic_density',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$

    call start_timer(tmr_std_chargescf)

    if (inode == ionode .and. iprint_SC >= 2) &
         write (io_lun,fmt='(2x,"Entering get_electronic_density")')

    do spin = 1, nspin
       gridfunctions(atom_fns_K)%griddata = zero
       call act_on_vectors_new(inode-1, rem_bucket(atomf_H_atomf_rem), &
                               matKatomf(spin), atom_fns_K, atom_fns)
       denout(:,spin) = zero
       i_count_alpha = 0
       do blk = 1, domain%groups_on_node
          n_point = (blk - 1) * n_pts_in_block
          i_count_alpha = (naba_atoms_of_blocks(atomf)%ibegin_blk_orb(blk)-1) * n_pts_in_block
          !if(blk == 1) then
          !   i_count_alpha = 0
          !else
          !   i_count_alpha= i_count_alpha+ &
          !        naba_atoms_of_blocks(atomf)%no_of_atom(blk-1)*nsf*n_pts_in_block
          !endif
          if (naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0) then   !TM 30/Jun/2003
             !do n_i=1, naba_atoms_of_blocks(atomf)%no_of_atom(blk)*NSF*n_pts_in_block, &
             do n_i = 1, naba_atoms_of_blocks(atomf)%no_of_orb(blk)*n_pts_in_block, n_pts_in_block
                do n=1, n_pts_in_block
                   denout(n_point+n, spin) =                                       &
                        denout(n_point+n, spin) +                                  &
                        gridfunctions(atom_fns_K)%griddata(i_count_alpha+n_i+n-1) * &
                        gridfunctions(atom_fns)%griddata(i_count_alpha+n_i+n-1)
                end do
             end do
          end if !(naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0)   !TM 30/Jun/2003
       end do ! blk

       ! FOR DEBUGGING   T. MIYAZAKI 30/Jun/2003
       !do n=1, size
       !   if (denout(n,spin) < -RD_ERR) then
       !      write(io_lun,*) ' WARNING!!!   density < 0 ',n, size, denout(n,spin)
       !   endif
       !enddo
       ! FOR DEBUGGING   T. MIYAZAKI 30/Jun/2003

       ! and calculate the electron number, note in the spin_polarised case
       ! the electrons corresponds to the number of electrons corresponding
       ! to different spin components. In the spin non-polarised form electrons
       ! corresponds to the total number of electrons
       electrons(spin) = grid_point_volume * &
                         rsum(n_my_grid_points, denout(:,spin), 1)
       call gsum(electrons(spin))

       if (inode == ionode .and. iprint_SC > 1) &
            write (io_lun, '(2x,"Electrons (spin=",i1,"): ",f25.15)') &
                  spin, electrons(spin)
    end do ! spin

    if (flag_Becke_weights) &
         call build_Becke_charges(atomcharge, denout, size)

    ! atom_fns_K is using the same memory as h_on_atomfns, (in other
    ! words we are using h_on_atomfns as temorary storage), so lets be
    ! safe and set it back to zero
    gridfunctions(atom_fns_K)%griddata = zero

    call stop_timer(tmr_std_chargescf)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_electronic_density',echo=.true.)
!****lat>$

    return
  end subroutine get_electronic_density
  !!***

  ! -----------------------------------------------------------
  ! Subroutine get_band_density
  ! -----------------------------------------------------------

  !!****f* density_module/get_band_density *
  !!
  !!  NAME
  !!   get_band_density
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Gets electronic density for one band on grid (closely follows
  !!   get_electronic_density, above)
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2015/06/05
  !!  MODIFICATION HISTORY
  !!   2016/07/15 18:30 nakata
  !!    Renamed sf_H_sf_rem -> atomf_H_atomf_rem
  !!   2016/07/20 16:30 nakata
  !!    Renamed naba_atm -> naba_atoms_of_blocks
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf instead of sf
  !!   2016/08/09 14:00 nakata
  !!    Renamed support -> atom_fns
  !!   2016/09/15 18:00 nakata
  !!    Renamed matBand -> matBand_atomfns
  !!  SOURCE
  !!
  subroutine get_band_density(denout, spin, atom_fns, atom_fns_K, matBand_atomfns, size)

    use datatypes
    use numbers
    use GenBlas,                     only: scal, rsum
    use dimens,                      only: n_my_grid_points, grid_point_volume
    use block_module,                only: n_pts_in_block
    use set_bucket_module,           only: rem_bucket, atomf_H_atomf_rem
    use calc_matrix_elements_module, only: act_on_vectors_new
    use primary_module,              only: domain
    use set_blipgrid_module,         only: naba_atoms_of_blocks
    use GenComms,                    only: gsum, inode, ionode
    use global_module,               only: iprint_SC, atomf, ni_in_cell
    use functions_on_grid,           only: gridfunctions, fn_on_grid

    implicit none

    ! Passed variables

    integer :: size, spin, matBand_atomfns
    integer :: atom_fns, atom_fns_K
    real(double), dimension(size) :: denout

    ! Local variables
    integer :: blk, i_count_alpha, n, n_i, n_point
    real(double) :: electrons

    if (inode == ionode .and. iprint_SC >= 2) &
         write (io_lun,fmt='(2x,"Entering get_band_density")')

    gridfunctions(atom_fns_K)%griddata = zero
    call act_on_vectors_new(inode-1, rem_bucket(atomf_H_atomf_rem), &
         matBand_atomfns, atom_fns_K, atom_fns)
    denout(:) = zero
    i_count_alpha = 0
    do blk = 1, domain%groups_on_node
       n_point = (blk - 1) * n_pts_in_block
       i_count_alpha = (naba_atoms_of_blocks(atomf)%ibegin_blk_orb(blk)-1) * n_pts_in_block
       if (naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0) then   !TM 30/Jun/2003
          do n_i = 1, naba_atoms_of_blocks(atomf)%no_of_orb(blk)*n_pts_in_block, n_pts_in_block
             do n=1, n_pts_in_block
                denout(n_point+n) =                                       &
                     denout(n_point+n) +                                  &
                     gridfunctions(atom_fns_K)%griddata(i_count_alpha+n_i+n-1) * &
                     gridfunctions(atom_fns)%griddata(i_count_alpha+n_i+n-1)
             end do
          end do
       end if !(naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0)   !TM 30/Jun/2003
    end do ! blk
    electrons = grid_point_volume * rsum(n_my_grid_points, denout(:), 1)
    call gsum(electrons)

    if (inode == ionode .and. iprint_SC > 1) &
         write (io_lun, '(2x,"Electrons: ",f25.15)') &
         electrons
    call stop_timer(tmr_std_chargescf)

    return
  end subroutine get_band_density
  !!***

  ! -----------------------------------------------------------
  ! Subroutine build_Becke_weights
  ! -----------------------------------------------------------

  !!****f* density_module/build_Becke_weights *
  !!
  !!  NAME
  !!   build_Becke_weights
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Builds Becke weights on grid for atoms
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D. R. Bowler & A. M. P. Sena
  !!  CREATION DATE
  !!   2009
  !!  MODIFICATION HISTORY
  !!   2011/08/03 10:44 dave
  !!    Moving all detailed cDFT work into cdft_module
  !!   2016/07/20 16:30 nakata
  !!    Renamed naba_atm -> naba_atoms_of_blocks
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf
  !!   2016/09/16 17:00 nakata
  !!    Introduced RadiusPAO
  !!  SOURCE
  !!
  subroutine build_Becke_weights

    use datatypes
    use GenBlas,             only: scal, rsum
    use numbers,             only: zero, RD_ERR, one, half
    use block_module,        only: n_pts_in_block
    use primary_module,      only: domain
    use set_blipgrid_module, only: naba_atoms_of_blocks
    use GenComms,            only: gsum, cq_abort, inode, ionode
    use global_module,       only: rcellx, rcelly, rcellz, iprint_SC, &
                                   atomf, paof, sf,                   &   ! nakata3
                                   ni_in_cell, species_glob,          &
                                   id_glob, io_lun,                   &
                                   flag_Becke_atomic_radii,           &
                                   flag_perform_cdft, flag_cdft_atom
    use cover_module,        only: DCS_parts
    use group_module,        only: blocks, parts
    use dimens,              only: RadiusSupport, RadiusPAO, atomicnum

    implicit none

    ! Local variables
    integer      :: blk, no_of_ib_ia, n, n_i, ipos, jpos, stat, at,   &
                    atomi, atomj, point
    integer      :: ipart, jpart, ind_part, ia, ii, icover, ig_atom,  &
                    iatom, iblock, the_species
    real(double) :: sum, mu, sij, sji, Rij
    real(double) :: dcellx_block, dcelly_block, dcellz_block
    real(double) :: xblock, yblock, zblock, aij, chi, uij, nu
    integer,      dimension(:),   allocatable :: npoint, globatom
    integer,      dimension(:,:), allocatable :: ip_store
    real(double), dimension(:),   allocatable :: xatom, yatom, zatom, &
                                                 rcut, rad
    real(double), dimension(:,:), allocatable :: r_store, x_store,    &
                                                 y_store, z_store

    call start_timer(tmr_std_chargescf)
    if(inode==ionode.AND.iprint_SC>=2) &
         write(io_lun,fmt='(2x,"Entering build_Becke_weights")')
    if(.NOT.weights_set) then
       call assign_atomic_radii
       weights_set = .true.
    end if
    no_of_ib_ia=0
    if(allocated(bw)) deallocate(bw)
    do blk=1, domain%groups_on_node
       no_of_ib_ia = no_of_ib_ia + naba_atoms_of_blocks(atomf)%no_of_atom(blk)
    end do
    no_of_ib_ia = no_of_ib_ia * n_pts_in_block
    allocate(bw(no_of_ib_ia),STAT=stat)
    bw = -one
    no_of_ib_ia=0
    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz
    ! atomcharge = zero   ! LT: seems to be redundant 2011/06/18
    if(flag_perform_cdft) bwgrid = zero
    do blk=1, domain%groups_on_node
       xblock=(domain%idisp_primx(blk)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(blk)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(blk)+domain%nz_origin-1)*dcellz_block
       if(naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0) then
          allocate(xatom(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                   &
                   yatom(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                   &
                   zatom(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                   &
                   rcut(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                    &
                   npoint(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                  &
                   globatom(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                &
                   ip_store(n_pts_in_block,naba_atoms_of_blocks(atomf)%no_of_atom(blk)), &
                   rad(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                     &
                   x_store(n_pts_in_block,naba_atoms_of_blocks(atomf)%no_of_atom(blk)),  &
                   y_store(n_pts_in_block,naba_atoms_of_blocks(atomf)%no_of_atom(blk)),  &
                   z_store(n_pts_in_block,naba_atoms_of_blocks(atomf)%no_of_atom(blk)),  &
                   r_store(n_pts_in_block,naba_atoms_of_blocks(atomf)%no_of_atom(blk)),  &
                   STAT=stat)
          ip_store = 0
          x_store = zero
          y_store = zero
          z_store = zero
          r_store = zero
          at = 0
          do ipart=1,naba_atoms_of_blocks(atomf)%no_of_part(blk)
             jpart=naba_atoms_of_blocks(atomf)%list_part(ipart,blk)
             if(jpart > DCS_parts%mx_gcover) then
                call cq_abort('build_Becke_weights: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart,blk)
                at=at+1
                ii = naba_atoms_of_blocks(atomf)%list_atom(at,blk)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)
                xatom(at)=DCS_parts%xcover(icover)
                yatom(at)=DCS_parts%ycover(icover)
                zatom(at)=DCS_parts%zcover(icover)
                the_species=species_glob(ig_atom)
                globatom(at) = ig_atom
!!! 2016.9.16 nakata3
                if (atomf.eq.sf) then
                   rcut(at) = RadiusSupport(the_species)
                else if (atomf.eq.paof) then
                   rcut(at) = RadiusPAO(the_species)
                endif
!!! nakata3 end
                rad(at) = atrad(atomicnum(the_species))
             end do
          end do
          call check_block(xblock,yblock,zblock,xatom,yatom,zatom,&
               rcut, npoint,ip_store, x_store, y_store, z_store,&
               r_store,n_pts_in_block,naba_atoms_of_blocks(atomf)%no_of_atom(blk))
          do point = 1, n_pts_in_block
             ! Sum is normalisation for grid point, Z = \sum_i bw_i
             sum = zero
             ! Actually we want a flag set in check_block so that this whole double loop
             ! is skipped if r_in criterion is met and the normalisation
             ! if(flag) then set weight for atom = 1, others = 0
             ! else
             do atomi = 1, naba_atoms_of_blocks(atomf)%no_of_atom(blk)
                ipos = no_of_ib_ia + n_pts_in_block*(atomi-1)
                if(ip_store(point,atomi)>0) then
                   ! Here we accumulate bw which is \prod_{j\ne i} s(mu_{ij})
                   if(naba_atoms_of_blocks(atomf)%no_of_atom(blk)==1) then
                      bw(ipos+point) = one
                   else if(r_store(point,atomi)<1.0e-6) then
                      bw(ipos+point) = one
                      do atomj = atomi+1,naba_atoms_of_blocks(atomf)%no_of_atom(blk)
                         jpos = no_of_ib_ia + n_pts_in_block*(atomj-1)
                         bw(jpos+point) = zero
                      end do
                   else
                      do atomj = atomi+1,naba_atoms_of_blocks(atomf)%no_of_atom(blk)
                         jpos = no_of_ib_ia + n_pts_in_block*(atomj-1)
                         if(ip_store(point,atomj)>0) then
                            if(r_store(point,atomj)<1.0e-6) then
                               bw(ipos+point) = zero
                               bw(jpos+point) = one
                            else
                               ! Do we want a pre-calculated set of these ?
                               Rij = sqrt ((xatom(atomi) - xatom(atomj)) * &
                                    (xatom(atomi) - xatom(atomj)) + &
                                    (yatom(atomi) - yatom(atomj)) * &
                                    (yatom(atomi) - yatom(atomj)) + &
                                    (zatom(atomi) - zatom(atomj)) * &
                                    (zatom(atomi) - zatom(atomj)))
                               ! make mu
                               mu = (r_store(point,atomi) - &
                                    r_store(point,atomj)) / Rij
                               if(flag_Becke_atomic_radii) then
                                  chi = rad(atomi) / rad(atomj)
                                  uij = (chi - one) / (chi + one)
                                  aij = uij / (uij * uij - one)
                                  if(aij<-half) aij = - half
                                  if(aij>half) aij = half
                                  !sij = s(mu)
                                  nu = mu+aij*(one-mu*mu)
                                  sij = s(nu)
                                  sji = s(-nu)
                               else
                                  sij = s(mu)
                                  sji = s(-mu)
                               end if
                               if(bw(ipos+point)<zero) bw(ipos+point) = one
                               if(bw(jpos+point)<zero) bw(jpos+point) = one
                               bw(ipos + point) = bw(ipos + point) * sij
                               bw(jpos + point) = bw(jpos + point) * sji
                            end if
                         else
                            bw(jpos+point) = zero
                         end if
                      end do ! atomj
                      if(bw(ipos+point)<zero) bw(ipos+point) = zero
                   end if
                   ! Accumulate w_i into p
                   sum = sum + bw(ipos+point)
                else
                   bw(ipos+point) = zero
                end if ! ip_store(point,atomi)>0
             end do ! atomi
             do atomi = 1, naba_atoms_of_blocks(atomf)%no_of_atom(blk)
                ipos = no_of_ib_ia + n_pts_in_block*(atomi-1)
                ! Normalise w
                if(abs(sum)>1.0e-8_double) &
                     bw(ipos+point) = bw(ipos+point) / sum
                if(flag_perform_cdft) then
                   if(flag_cdft_atom(globatom(atomi)) > 0) then
                      bwgrid(point + (blk - 1) * n_pts_in_block, &
                           flag_cdft_atom(globatom(atomi))) = &
                           bwgrid(point + (blk - 1) * n_pts_in_block,&
                           flag_cdft_atom(globatom(atomi))) + &
                           bw(ipos + point)
                   end if
                end if
             end do ! atomi
             ! end if
          end do ! point
          no_of_ib_ia = no_of_ib_ia + n_pts_in_block*naba_atoms_of_blocks(atomf)%no_of_atom(blk)
          deallocate(xatom, yatom, zatom, rcut, npoint, globatom, ip_store, &
               x_store, y_store, z_store, r_store, rad, STAT=stat)
       endif !(naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0)   !TM 30/Jun/2003
    end do ! blk
  end subroutine build_Becke_weights
  !!***

  ! -----------------------------------------------------------
  ! Subroutine build_Becke_weight_forces
  ! -----------------------------------------------------------

  !!****f* density_module/build_Becke_weight_forces *
  !!
  !!  NAME
  !!   build_Becke_weight_forces
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Builds contributions to forces from cDFT using Becke weights
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D. R. Bowler and A. M. P. Sena
  !!  CREATION DATE
  !!   2009
  !!  MODIFICATION HISTORY
  !!    2011/05/20 L.Tong
  !!      * Removed input chden and size, redundant, as the density data is
  !!        read from the module variable density.
  !!      * Added correction for spin polarisation: for spin polarisation
  !!      density is contribution from spin up and density_dn is contribution
  !!      from spin down and hence elec_here should be the sum of both
  !!      spin components.
  !!   2011/10/06 13:58 dave
  !!    Completed correct force form
  !!   2012/03/13 L.Tong
  !!    Rewrote spin implementation
  !!   2016/07/20 16:30 nakata
  !!    Renamed naba_atm -> naba_atoms_of_blocks
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf
  !!   2016/09/16 17:00 nakata
  !!    Introduced RadiusPAO
  !!  SOURCE
  !!
  subroutine build_Becke_weight_forces(weight_force)

    use datatypes
    use GenBlas,             only: scal, rsum
    use numbers,             only: zero, RD_ERR, two, one, half
    use block_module,        only: n_pts_in_block
    use primary_module,      only: domain
    use set_blipgrid_module, only: naba_atoms_of_blocks
    use GenComms,            only: gsum, cq_abort, inode, ionode
    use global_module,       only: rcellx, rcelly, rcellz, iprint_SC, &
                                   atomf, paof, sf,                   &
                                   ni_in_cell, species_glob,          &
                                   id_glob, io_lun, flag_cdft_atom,   &
                                   flag_Becke_atomic_radii,           &
                                   nspin, spin_factor
    use cover_module,        only: DCS_parts
    use group_module,        only: blocks, parts
    use dimens,              only: RadiusSupport, RadiusPAO,          &
                                   grid_point_volume,                 &
                                   atomicnum
    use cdft_data,           only: cDFT_Type, cDFT_Fix_Charge,        &
                                   cDFT_Fix_ChargeDifference, cDFT_Vc

    implicit none

    ! Passed variables
    real(double), dimension(3,ni_in_cell) :: weight_force
    ! Local variables
    integer      :: blk, no_of_ib_ia, n, n_i, ipos, jpos, stat, at,    &
                    atomi, atomj, point
    integer      :: ipart, jpart, ind_part, ia, ii, icover, ig_atom,   &
                    iatom, iblock, the_species, gi, spin
    real(double) :: sum, mu, sij, sji, Rij, dwjbydix, dwjbydiy,        &
                    dwjbydiz, elec_here
    real(double) :: dmux, dmuy, dmuz, tij, tji, sum1x, sum2x, sum1y,   &
                    sum2y, sum1z, sum2z
    real(double) :: dcellx_block, dcelly_block, dcellz_block
    real(double) :: xblock, yblock, zblock, checkE, aij, chi, uij, nu, &
                    dnu
    integer,      dimension(:),   allocatable :: npoint, globatom
    integer,      dimension(:,:), allocatable :: ip_store
    real(double), dimension(:),   allocatable :: xatom, yatom, zatom,  &
                                                 rcut, rad
    real(double), dimension(:,:), allocatable :: r_store, x_store,     &
                                                 y_store, z_store
    real(double), dimension(:),   allocatable :: sum0x, sum0y, sum0z



    call start_timer(tmr_std_chargescf)
    weight_force = zero
    if (inode == ionode .AND. iprint_SC >= 2)&
         write(io_lun, fmt='(2x,"Entering build_Becke_weight_forces")')
    no_of_ib_ia = 0
    dcellx_block = rcellx / blocks%ngcellx
    dcelly_block = rcelly / blocks%ngcelly
    dcellz_block = rcellz / blocks%ngcellz
    ! atomcharge = zero       ! LT seems redundant 2011/06/18
    do blk=1, domain%groups_on_node
       xblock = (domain%idisp_primx(blk) + domain%nx_origin-1) * dcellx_block
       yblock = (domain%idisp_primy(blk) + domain%ny_origin-1) * dcelly_block
       zblock = (domain%idisp_primz(blk) + domain%nz_origin-1) * dcellz_block
       if (naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0) then
          allocate(xatom(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                   &
                   yatom(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                   &
                   zatom(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                   &
                   rcut(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                    &
                   npoint(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                  &
                   globatom(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                &
                   rad(naba_atoms_of_blocks(atomf)%no_of_atom(blk)),                     &
                   ip_store(n_pts_in_block,naba_atoms_of_blocks(atomf)%no_of_atom(blk)), &
                   x_store(n_pts_in_block,naba_atoms_of_blocks(atomf)%no_of_atom(blk)),  &
                   y_store(n_pts_in_block,naba_atoms_of_blocks(atomf) %no_of_atom(blk)), &
                   z_store(n_pts_in_block,naba_atoms_of_blocks(atomf)%no_of_atom(blk)),  &
                   r_store(n_pts_in_block,naba_atoms_of_blocks(atomf)%no_of_atom(blk)),  &
                   STAT=stat)
          ip_store = 0
          x_store = zero
          y_store = zero
          z_store = zero
          r_store = zero
          ! Replace this with part/atom loop !
          ! If there's only one atom there's no force
          if (naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 1) then
             at = 0
             do ipart = 1, naba_atoms_of_blocks(atomf)%no_of_part(blk)
                jpart = naba_atoms_of_blocks(atomf)%list_part(ipart,blk)
                if (jpart > DCS_parts%mx_gcover) then
                   call cq_abort('build_Becke_weight_forces: JPART ERROR ', &
                                 ipart, jpart)
                endif
                ind_part=DCS_parts%lab_cell(jpart)
                do ia = 1, naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart,blk)
                   at = at + 1
                   ii = naba_atoms_of_blocks(atomf)%list_atom(at,blk)
                   icover = DCS_parts%icover_ibeg(jpart)+ii-1
                   ig_atom = id_glob(parts%icell_beg(ind_part)+ii-1)
                   xatom(at) = DCS_parts%xcover(icover)
                   yatom(at) = DCS_parts%ycover(icover)
                   zatom(at) = DCS_parts%zcover(icover)
                   the_species = species_glob(ig_atom)
                   globatom(at) = ig_atom
!!! 2016.9.16 nakata3
                   if (atomf.eq.sf) then
                      rcut(at) = RadiusSupport(the_species)
                   else if (atomf.eq.paof) then
                      rcut(at) = RadiusPAO(the_species)
                   endif
!!! nakata3 end
                   rad(at) = atrad(atomicnum(the_species))
                end do
             end do
             call check_block(xblock, yblock, zblock, xatom, yatom,   &
                              zatom, rcut, npoint, ip_store, x_store, &
                              y_store, z_store,r_store,               &
                              n_pts_in_block,naba_atoms_of_blocks(atomf)% &
                              no_of_atom(blk))
             allocate(sum0x(naba_atoms_of_blocks(atomf)%no_of_atom(blk)), &
                      sum0y(naba_atoms_of_blocks(atomf)%no_of_atom(blk)), &
                      sum0z(naba_atoms_of_blocks(atomf)%no_of_atom(blk)))
             do point = 1, n_pts_in_block
                elec_here = zero
                do spin = 1, nspin
                   elec_here = elec_here + spin_factor * &
                               density(n_pts_in_block*(blk-1)+point,spin) * &
                               grid_point_volume
                end do
                ! Actually we want a flag set in check_block so that
                ! this whole double loop is skipped if r_in criterion
                ! is met and the normalisation
                ! if(flag) then set weight for atom = 1, others = 0
                ! else
                ! First time around we need to store the sums over t
                ! and \nabla \mu as well as storing t and nabla
                do atomi = 1, naba_atoms_of_blocks(atomf)%no_of_atom(blk)
                   ipos = no_of_ib_ia + n_pts_in_block*(atomi-1)
                   if (ip_store(point,atomi) > 0) then
                      ! Here we accumulate bw which is \prod_{j\ne i} s(mu_{ij})
                      if (naba_atoms_of_blocks(atomf)%no_of_atom(blk) == 1 .or. &
                          r_store(point, atomi) < 1.0e-6) then
                         sum0x = zero
                         sum0y = zero
                         sum0z = zero
                         sum1x = zero
                         sum1y = zero
                         sum1z = zero
                         sum2x = zero
                         sum2y = zero
                         sum2z = zero
                      else
                         sum0x = zero
                         sum0y = zero
                         sum0z = zero
                         sum1x = zero
                         sum1y = zero
                         sum1z = zero
                         sum2x = zero
                         sum2y = zero
                         sum2z = zero
                         do atomj = 1, naba_atoms_of_blocks(atomf)%no_of_atom(blk)
                            if (atomj /= atomi) then
                               jpos = no_of_ib_ia + n_pts_in_block*(atomj-1)
                               if (ip_store(point,atomj) > 0) then
                                  if (abs(bw(jpos+point)-one) < 1.0e-8_double .or. &
                                      abs(bw(jpos+point))<1.0e-8_double) then
                                     sum0x(atomj) = zero
                                     sum0y(atomj) = zero
                                     sum0z(atomj) = zero
                                  else
                                     if (r_store(point,atomi) < 1.0e-6)       &
                                          write (*, *) 'Zero distance: ',     &
                                                        point, atomi,         &
                                                        r_store(point,atomi), &
                                                        bw(ipos+point),       &
                                                        bw(jpos+point)
                                     ! Do we want a pre-calculated set of these ?
                                     Rij = sqrt((xatom(atomi)-xatom(atomj)) * &
                                                (xatom(atomi)-xatom(atomj)) + &
                                                (yatom(atomi)-yatom(atomj)) * &
                                                (yatom(atomi)-yatom(atomj)) + &
                                                (zatom(atomi)-zatom(atomj)) * &
                                                (zatom(atomi)-zatom(atomj)))
                                     ! make mu
                                     mu = (r_store(point,atomi) -&
                                          & r_store(point,atomj)) / Rij
                                     if (flag_Becke_atomic_radii) then
                                        chi = rad(atomi) / rad(atomj)
                                        uij = (chi-one) / (chi+one)
                                        aij = uij / (uij*uij-one)
                                        if (aij < -half) aij = -half
                                        if (aij > half) aij = half
                                        !sij = s(mu)
                                        nu = mu+aij*(one-mu*mu)
                                        dnu = one-two*aij*mu
                                        dmux = (x_store(point,atomi) /             &
                                                (Rij*r_store(point,atomi)) -       &
                                                (xatom(atomi)-xatom(atomj)) * mu / &
                                                (Rij * Rij)) * (one - two*mu*aij)
                                        dmuy = (y_store(point,atomi) /             &
                                                (Rij*r_store(point,atomi)) -       &
                                                (yatom(atomi)-yatom(atomj)) * mu / &
                                                (Rij*Rij)) * (one - two*mu*aij)
                                        dmuz = (z_store(point,atomi) /             &
                                                (Rij*r_store(point,atomi)) -       &
                                                (zatom(atomi)-zatom(atomj)) * mu / &
                                                (Rij*Rij)) * (one - two*mu*aij)
                                        tij = t( nu)!*dnu
                                        tji = t(-nu)!*dnu
                                     else
                                        dmux = x_store(point,atomi) /             &
                                               (Rij*r_store(point,atomi)) -       &
                                               (xatom(atomi)-xatom(atomj)) * mu / &
                                               (Rij*Rij)
                                        dmuy = y_store(point,atomi) /             &
                                               (Rij*r_store(point,atomi)) -       &
                                               (yatom(atomi)-yatom(atomj)) * mu / &
                                               (Rij*Rij)
                                        dmuz = z_store(point,atomi) /             &
                                               (Rij*r_store(point,atomi)) -       &
                                               (zatom(atomi)-zatom(atomj)) * mu / &
                                               (Rij*Rij)
                                        tij = t( mu)
                                        tji = t(-mu)
                                     end if
                                     sum0x(atomj) = tji*dmux
                                     sum0y(atomj) = tji*dmuy
                                     sum0z(atomj) = tji*dmuz
                                     sum1x = sum1x + tij*dmux
                                     sum1y = sum1y + tij*dmuy
                                     sum1z = sum1z + tij*dmuz
                                     sum2x = sum2x + bw(jpos+point)*tji*dmux
                                     sum2y = sum2y + bw(jpos+point)*tji*dmuy
                                     sum2z = sum2z + bw(jpos+point)*tji*dmuz
                                  end if
                               end if
                            end if
                         end do ! atomj
                         do atomj = 1,naba_atoms_of_blocks(atomf)%no_of_atom(blk)
                            jpos = no_of_ib_ia + n_pts_in_block*(atomj-1)
                            if (atomi == atomj) then
                               dwjbydix = bw(jpos+point) * sum1x - &
                                          bw(jpos+point) * &
                                          (bw(ipos+point) * sum1x - sum2x)
                               dwjbydiy = bw(jpos+point) * sum1y - &
                                          bw(jpos+point) * &
                                          (bw(ipos+point) * sum1y - sum2y)
                               dwjbydiz = bw(jpos+point) * sum1z - &
                                          bw(jpos+point) * &
                                          (bw(ipos+point) * sum1z - sum2z)
                            else
                               dwjbydix = - bw(jpos+point) * &
                                          sum0x(atomj) - bw(jpos+point) * &
                                          (bw(ipos+point) * sum1x - sum2x)
                               dwjbydiy = - bw(jpos+point) * &
                                          sum0y(atomj) - bw(jpos+point) * &
                                          (bw(ipos+point) * sum1y - sum2y)
                               dwjbydiz = - bw(jpos+point) * &
                                          sum0z(atomj) - bw(jpos+point) * &
                                          (bw(ipos+point) * sum1z - sum2z)
                            end if
                            gi = globatom(atomi)
                            if (cDFT_Type == cDFT_Fix_ChargeDifference) then
                               if (flag_cdft_atom(globatom(atomj)) == 1) then
                                  weight_force(1,gi) =      &
                                       weight_force(1,gi) - &
                                       cDFT_Vc(1) * elec_here * dwjbydix
                                  weight_force(2,gi) =      &
                                       weight_force(2,gi) - &
                                       cDFT_Vc(1) * elec_here * dwjbydiy
                                  weight_force(3,gi) =      &
                                       weight_force(3,gi) - &
                                       cDFT_Vc(1) * elec_here * dwjbydiz
                               else if (flag_cdft_atom(globatom(atomj)) == 2) then
                                  weight_force(1,gi) =      &
                                       weight_force(1,gi) + &
                                       cDFT_Vc(1) * elec_here * dwjbydix
                                  weight_force(2,gi) =      &
                                       weight_force(2,gi) + &
                                       cDFT_Vc(1) * elec_here * dwjbydiy
                                  weight_force(3,gi) =      &
                                       weight_force(3,gi) + &
                                       cDFT_Vc(1) * elec_here * dwjbydiz
                               end if
                            else
                               if (flag_cdft_atom(globatom(atomi)) > 0) then
                                  weight_force(1,gi) =               &
                                       weight_force(1,gi) -          &
                                       cDFT_Vc(flag_cdft_atom(gi)) * &
                                       elec_here * dwjbydix
                                  weight_force(2,gi) =               &
                                       weight_force(2,gi) -          &
                                       cDFT_Vc(flag_cdft_atom(gi)) * &
                                       elec_here * dwjbydiy
                                  weight_force(3,gi) =               &
                                       weight_force(3,gi) -          &
                                       cDFT_Vc(flag_cdft_atom(gi)) * &
                                       elec_here * dwjbydiz
                               end if
                            end if
                         end do ! atom j
                      end if
                   end if
                end do ! atomi
             end do ! point
             deallocate (sum0x, sum0y, sum0z)
          end if ! naba_atoms_of_blocks%no_of_atom(blk)>1
          no_of_ib_ia = no_of_ib_ia + n_pts_in_block * &
                        naba_atoms_of_blocks(atomf)%no_of_atom(blk)
          deallocate(xatom, yatom, zatom, rcut, npoint, globatom, &
                     ip_store, rad, x_store, y_store, z_store,    &
                     r_store, STAT=stat)
       endif !(naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0)   !TM 30/Jun/2003
    end do ! blk
    call gsum(weight_force, 3, ni_in_cell)
    return
  end subroutine build_Becke_weight_forces
  !!***

  ! -----------------------------------------------------------
  ! Subroutine build_Becke_weight_matrix
  ! -----------------------------------------------------------

  !!****f* density_module/build_Becke_weight_matrix *
  !!
  !!  NAME
  !!   build_Becke_weight_matrix
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Builds cDFT weight matrix for constrained atoms using Becke weights
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D. R. Bowler and A. M. P. Sena
  !!  CREATION DATE
  !!   2009
  !!  MODIFICATION HISTORY
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomfns
  !!   2016/07/15 18:30 nakata
  !!    Renamed sf_H_sf_rem -> atomf_H_atomf_rem
  !!   2016/07/20 16:30 nakata
  !!    Renamed naba_atm -> naba_atoms_of_blocks
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf instead of sf
  !!    Removed unused pao_H_sf_rem
  !!   2016/08/08 15:30 nakata
  !!    Renamed supportfns -> atomfns
  !!  SOURCE
  !!
  subroutine build_Becke_weight_matrix(matWc, ngroups)

    use datatypes
    use numbers
    use global_module,               only: atomf
    use block_module,                only: n_blocks, n_pts_in_block
    use primary_module,              only: domain
    use set_blipgrid_module,         only: naba_atoms_of_blocks
    use functions_on_grid,           only: gridfunctions, fn_on_grid, &
                                           H_on_atomfns, atomfns
    use calc_matrix_elements_module, only: norb
    use set_bucket_module,           only: rem_bucket, atomf_H_atomf_rem
    use calc_matrix_elements_module, only: get_matrix_elements_new
    use GenComms,                    only: inode

    implicit none

    ! input parameters
    integer :: ngroups
    integer, dimension(ngroups) :: matWc
    ! local variables
    integer :: n, m, nb, atom, nsf1, point, i

    ! Create w_c(r)|chi_j>
    do i=1,ngroups
       gridfunctions(H_on_atomfns(1))%griddata = zero
       n = 0
       m = 0
       do nb=1,domain%groups_on_node
          if(naba_atoms_of_blocks(atomf)%no_of_atom(nb)>0) then
             do atom = 1, naba_atoms_of_blocks(atomf)%no_of_atom(nb)
                do nsf1 = 1, norb(naba_atoms_of_blocks(atomf),atom,nb)
                   do point = 1, n_pts_in_block
                      n = n + 1
                      gridfunctions(H_on_atomfns(1))%griddata(n) = &
                           gridfunctions(atomfns)%griddata(n) * &
                           bwgrid(m+point,i)
                   end do
                end do
             end do
          end if ! (naba_atoms_of_blocks(atomf)%no_of_atom(nb)>0)
          m = m + n_pts_in_block
       end do
       call get_matrix_elements_new(inode-1, rem_bucket(atomf_H_atomf_rem), &
                                    matWc(i), atomfns, H_on_atomfns(1))
    end do
  end subroutine build_Becke_weight_matrix
  !!***

  ! -----------------------------------------------------------
  ! Subroutine build_Becke_charges
  ! -----------------------------------------------------------

  !!****f* density_module/build_Becke_charges *
  !!
  !!  NAME
  !!   build_Becke_charges
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Builds Becke charges for atoms
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2009
  !!  MODIFICATION HISTORY
  !!    2011/05/20 L.Tong
  !!      Added correction for spin polarisation, added optional spin
  !!      for controlling which atomcharge or atomcharge_dn is calculated.
  !!   Wednesday, 2011/08/03 L.Tong
  !!    Fixed the bug caused by the conditional if (present (spin)
  !!    .and. spin == 2), on some compilers which does not short circuit
  !!    and evaluates both sides of .and. and when the subroutine is
  !!    called without the optional spin parameter, then spin == 2 gives
  !!    a segmentation fault. Fixed this by introducing logical variable
  !!    do_spin_down, we work out this variable correctly first, and
  !!    then apply conditional on it in place of the previous present
  !!    (spin) .and. spin == 2 conditional. Also removed the redoundant
  !!    flag_spin_polarisation. The presense of spin parameter is
  !!    sufficient in telling the subroutine it is calculating in spin
  !!    polarised case
  !!   2011/11/27 L.Tong
  !!    - Made atomcharge densities subroutine input parameter rather than
  !!      dependent on module global. This makes implementation for spin
  !!      polarised calculations easier.
  !!    - Changed chden and defined atomch to assumed shape array
  !!   2012/03/13 L.Tong
  !!    - Rewrote the spin implementation
  !!    - now the subroutine calculates both spin channels. The spin
  !!      non-polarised case is treated just as the single spin up
  !!      channel in fixed spin population case.
  !!   2013/03/06 17:03 dave
  !!   - spin_factor added in output
  !!   2016/07/20 16:30 nakata
  !!    Renamed naba_atm -> naba_atoms_of_blocks
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf instead of sf
  !!  SOURCE
  !!
  subroutine build_Becke_charges(atomch, chden, size)

    use datatypes
    use GenComms,            only: gsum, cq_abort, inode, ionode
    use numbers,             only: zero, RD_ERR, two, one
    use global_module,       only: atomf, ni_in_cell, id_glob, io_lun, &
                                   iprint_SC, nspin, spin_factor
    use primary_module,      only: domain
    use set_blipgrid_module, only: naba_atoms_of_blocks
    use cover_module,        only: DCS_parts
    use dimens,              only: grid_point_volume
    use block_module,        only: n_pts_in_block
    use group_module,        only: parts

    implicit none

    ! Passed variables
    integer                      :: size
    real(double), dimension(:,:) :: atomch
    real(double), dimension(:,:) :: chden

    ! Local variables
    integer :: blk, no_of_ib_ia, at, point, ipos, spin
    integer :: ipart, jpart, ind_part, ia, ii, icover, ig_atom, iatom, &
               iblock, the_species
    logical :: do_spin_down

    atomch = zero
    no_of_ib_ia = 0
    do blk = 1, domain%groups_on_node
       if (naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0) then
          at = 0
          do ipart = 1, naba_atoms_of_blocks(atomf)%no_of_part(blk)
             jpart = naba_atoms_of_blocks(atomf)%list_part(ipart, blk)
             if (jpart > DCS_parts%mx_gcover) then
                call cq_abort('PAO_to_grid_global: JPART ERROR ', &
                              ipart, jpart)
             endif
             ind_part = DCS_parts%lab_cell(jpart)
             do ia = 1, naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart, blk)
                at = at + 1
                ii = naba_atoms_of_blocks(atomf)%list_atom(at, blk)
                icover= DCS_parts%icover_ibeg(jpart) + ii - 1
                ig_atom= id_glob(parts%icell_beg(ind_part) + ii - 1)
                do point = 1, n_pts_in_block
                   ipos = no_of_ib_ia + n_pts_in_block * (at - 1)
                   do spin = 1, nspin
                      atomch(ig_atom,spin) = atomch(ig_atom,spin) + &
                           bw(ipos + point) * &
                           chden(point+(blk-1)*n_pts_in_block,spin)
                   end do ! spin
                end do ! point
             end do ! ia
          end do ! ipart
          no_of_ib_ia = no_of_ib_ia + n_pts_in_block * &
                        naba_atoms_of_blocks(atomf)%no_of_atom(blk)
       end if !(naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0)   !TM 30/Jun/2003
    end do ! blk
    !write(*,*) 'Done blocks'
    do spin = 1, nspin
       atomch(:,spin) = atomch(:,spin) * grid_point_volume
       call gsum(atomch(:,spin), ni_in_cell)
       ! print out the atom charge density information
       if (inode == ionode .and. iprint_SC > 2) then
          if (nspin == 1) then
             do blk = 1, ni_in_cell
                write (io_lun, &
                       fmt='(2x,"Atom ",i4," Becke charge ",f20.12)') &
                       blk, spin_factor*atomch(blk,spin)
             end do
          else
             do blk = 1, ni_in_cell
                write (io_lun, &
                       fmt='(2x,"Atom ",i4," Becke charge (spin=",i1,") ",f20.12)') &
                       blk, spin, atomch(blk,spin)
             end do
          end if
       end if
    end do ! spin

    return
  end subroutine build_Becke_charges
  !!***

  ! -----------------------------------------------------------
  ! Function s
  ! -----------------------------------------------------------

  !!****f* density_module/s *
  !!
  !!  NAME
  !!   s
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates s function for Becke weights
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2009
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  real(double) function s(mu)

    use numbers, only: zero, half, one

    implicit none

    real(double) :: mu, mua, mua2
    real(double), parameter :: a = 0.64_double
    real(double), parameter :: c1 = 2.1875_double
    real(double), parameter :: c3 = 2.1875_double
    real(double), parameter :: c5 = 1.3125_double
    real(double), parameter :: c7 = 0.3125_double

    if(mu<-a) then
       s = one
    else if(mu>a) then
       s = zero
    else
       mua = mu/a
       mua2 = mua*mua
       s = half*(one -  mua*(c1 - mua2*(c3 - mua2*(c5 - c7*mua2))))
    end if
  end function s
  !!***

  ! -----------------------------------------------------------
  ! Function t
  ! -----------------------------------------------------------

  !!****f* density_module/t *
  !!
  !!  NAME
  !!   t
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates t function for Becke weights
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2009
  !!  MODIFICATION HISTORY
  !!   2013/03/06 17:04 dave
  !!   - Added check for size of s
  !!  SOURCE
  !!
  real(double) function t(mu)

    use numbers, only: zero, half, one, three, five, seven, RD_ERR

    implicit none

    real(double) :: mu, mua, mua2, s
    real(double), parameter :: a = 0.64_double
    real(double), parameter :: reca = 1.5625_double
    real(double), parameter :: c1 = 2.1875_double
    real(double), parameter :: c3 = 2.1875_double
    real(double), parameter :: c5 = 1.3125_double
    real(double), parameter :: c7 = 0.3125_double

    if(mu<=-a) then
       t = zero
    else if(mu>=a) then
       t = zero
    else
       mua = mu/a
       mua2 = mua*mua
       t = -half*(reca*c1 - mua2*(three*reca*c3 -&
           mua2*(five*reca*c5 - seven*reca*c7*mua2)))
       s = half*(one -  mua*(c1 - mua2*(c3 - mua2*(c5 - c7*mua2))))
       if(abs(s)<RD_ERR) then
          !write(io_lun,fmt='(2x,"Warning: s too small in t(mu)",2f12.5)') s,mu
          t = zero
       else
          t = t/s
       end if
    end if
  end function t
  !!***


  !!****f* density_module/assign_atomic_radii *
  !!
  !!  NAME
  !!   assign_atomic_radii
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Assigns atomic radii; these are from J. C. Slater,
  !!   J. Chem. Phys. 41, 3199 (1964) H altered to 0.35 A following
  !!   Becke, J. Chem. Phys. 88, 2547 (1988) Noble gas and At, Rn, Fr
  !!   radii from Dalton Trans., 2008, 2832 (we may want to use all
  !!   from here ?)
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2009
  !!  MODIFICATION HISTORY
  !!   2011/07/20 17:47 dave
  !!    Extended table list significantly
  !!  SOURCE
  !!
  subroutine assign_atomic_radii

    use datatypes
    use numbers, only: zero
    use units, only: AngToBohr

    implicit none

    atrad = zero
    atrad(1)=0.35_double*AngToBohr
    atrad(2)=0.28_double*AngToBohr

    atrad(3)=1.45_double*AngToBohr ! Li
    atrad(4)=1.05_double*AngToBohr
    atrad(5)=0.85_double*AngToBohr
    atrad(6)=0.70_double*AngToBohr
    atrad(7)=0.65_double*AngToBohr
    atrad(8)=0.60_double*AngToBohr
    atrad(9)=0.50_double*AngToBohr ! F
    atrad(10)=0.58_double*AngToBohr

    atrad(11)=1.80_double*AngToBohr ! Na
    atrad(12)=1.50_double*AngToBohr
    atrad(13)=1.25_double*AngToBohr
    atrad(14)=1.10_double*AngToBohr
    atrad(15)=1.00_double*AngToBohr
    atrad(16)=1.00_double*AngToBohr
    atrad(17)=1.00_double*AngToBohr ! Cl
    atrad(18)=1.06_double*AngToBohr

    atrad(19)=2.20_double*AngToBohr ! K
    atrad(20)=1.80_double*AngToBohr
    atrad(21)=1.60_double*AngToBohr
    atrad(22)=1.40_double*AngToBohr
    atrad(23)=1.35_double*AngToBohr
    atrad(24)=1.40_double*AngToBohr
    atrad(25)=1.40_double*AngToBohr
    atrad(26)=1.40_double*AngToBohr ! Fe
    atrad(27)=1.35_double*AngToBohr
    atrad(28)=1.35_double*AngToBohr
    atrad(29)=1.35_double*AngToBohr
    atrad(30)=1.35_double*AngToBohr
    atrad(31)=1.30_double*AngToBohr ! Ga
    atrad(32)=1.25_double*AngToBohr
    atrad(33)=1.15_double*AngToBohr
    atrad(34)=1.15_double*AngToBohr
    atrad(35)=1.15_double*AngToBohr ! Br
    atrad(36)=1.16_double*AngToBohr

    atrad(37)=2.35_double*AngToBohr ! Rb
    atrad(38)=2.00_double*AngToBohr
    atrad(39)=1.80_double*AngToBohr
    atrad(40)=1.55_double*AngToBohr ! Zr
    atrad(41)=1.45_double*AngToBohr
    atrad(42)=1.45_double*AngToBohr
    atrad(43)=1.35_double*AngToBohr
    atrad(44)=1.30_double*AngToBohr ! Ru
    atrad(45)=1.35_double*AngToBohr
    atrad(46)=1.40_double*AngToBohr
    atrad(47)=1.60_double*AngToBohr ! Ag
    atrad(48)=1.55_double*AngToBohr
    atrad(49)=1.55_double*AngToBohr ! In
    atrad(50)=1.45_double*AngToBohr
    atrad(51)=1.45_double*AngToBohr
    atrad(52)=1.40_double*AngToBohr
    atrad(53)=1.40_double*AngToBohr ! I
    atrad(54)=1.40_double*AngToBohr

    atrad(55)=2.60_double*AngToBohr ! Cs
    atrad(56)=2.15_double*AngToBohr
    atrad(57)=1.95_double*AngToBohr ! La
    atrad(58)=1.85_double*AngToBohr
    atrad(59)=1.85_double*AngToBohr
    atrad(60)=1.85_double*AngToBohr
    atrad(61)=1.85_double*AngToBohr
    atrad(62)=1.85_double*AngToBohr
    atrad(63)=1.85_double*AngToBohr ! Eu
    atrad(64)=1.80_double*AngToBohr
    atrad(65)=1.75_double*AngToBohr
    atrad(66)=1.75_double*AngToBohr
    atrad(67)=1.75_double*AngToBohr
    atrad(68)=1.75_double*AngToBohr
    atrad(69)=1.75_double*AngToBohr
    atrad(70)=1.75_double*AngToBohr
    atrad(71)=1.75_double*AngToBohr ! Lu
    atrad(72)=1.55_double*AngToBohr
    atrad(73)=1.45_double*AngToBohr
    atrad(74)=1.35_double*AngToBohr ! W
    atrad(75)=1.35_double*AngToBohr
    atrad(76)=1.30_double*AngToBohr
    atrad(77)=1.35_double*AngToBohr ! Ir
    atrad(78)=1.35_double*AngToBohr
    atrad(79)=1.35_double*AngToBohr
    atrad(80)=1.50_double*AngToBohr ! Hg
    atrad(81)=1.90_double*AngToBohr
    atrad(82)=1.80_double*AngToBohr
    atrad(83)=1.60_double*AngToBohr
    atrad(84)=1.90_double*AngToBohr ! Po
    atrad(85)=1.50_double*AngToBohr
    atrad(86)=1.50_double*AngToBohr

    atrad(87)=2.60_double*AngToBohr
    atrad(88)=2.15_double*AngToBohr ! Ra
    atrad(89)=1.95_double*AngToBohr
    atrad(90)=1.80_double*AngToBohr
    atrad(91)=1.80_double*AngToBohr
    atrad(92)=1.75_double*AngToBohr
    atrad(93)=1.75_double*AngToBohr
    atrad(94)=1.75_double*AngToBohr
    atrad(95)=1.75_double*AngToBohr
    return
  end subroutine assign_atomic_radii
  !!***

  ! -----------------------------------------------------------
  ! Subroutine check_block
  ! -----------------------------------------------------------

  !!****f* density_module/check_block *
  !!
  !!  NAME
  !!   check_block
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Checks distances for blocks
  !!   N.B. This has been constructed for THIS MODULE and specifically
  !!   for the Becke weight scheme; the minus sign in the x, y, z below
  !!   are needed for this.
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   TM
  !!  CREATION DATE
  !!   2005 ?
  !!  MODIFICATION HISTORY
  !!   2011/10/06 14:00 dave
  !!    Finalising for cDFT
  !!  SOURCE
  !!
  subroutine check_block(xblock, yblock, zblock, xatom, yatom, zatom, &
                         rcut, npoint, ip_store, x_store, y_store,    &
                         z_store, r_store, blocksize, natoms)

    use numbers
    use global_module, only: rcellx,rcelly,rcellz
    use group_module,  only: blocks
    use block_module,  only: nx_in_block, ny_in_block, nz_in_block!, &
    ! n_pts_in_block

    implicit none

    ! Passed
    integer :: blocksize, natoms
    real(double), intent(in) :: xblock, yblock, zblock
    integer, intent(out) :: npoint(natoms), ip_store(blocksize,natoms)
    real(double), dimension(natoms), intent(in)  :: &
         xatom, yatom, zatom, rcut
    real(double), dimension(blocksize,natoms), intent(out) :: &
         r_store, x_store, y_store, z_store
    ! Local
    real(double) :: dcellx_block, dcelly_block, dcellz_block
    real(double) :: dcellx_grid, dcelly_grid, dcellz_grid
    real(double) :: dx, dy, dz
    integer      :: ipoint, iz, iy, ix, at
    real(double) :: r2, r_from_i, rx, ry, rz, x, y, z, rcut2

    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz

    dcellx_grid=dcellx_block/nx_in_block
    dcelly_grid=dcelly_block/ny_in_block
    dcellz_grid=dcellz_block/nz_in_block

    ! Add loop to find R_in (shortest atom-atom distance)
    ipoint=0
    npoint=0
    do iz=1,nz_in_block
       do iy=1,ny_in_block
          do ix=1,nx_in_block
             ipoint=ipoint+1

             dx=dcellx_grid*(ix-1)
             dy=dcelly_grid*(iy-1)
             dz=dcellz_grid*(iz-1)
             do at = 1,natoms
                rcut2 = rcut(at)* rcut(at)
                rx=xblock+dx-xatom(at)
                ry=yblock+dy-yatom(at)
                rz=zblock+dz-zatom(at)
                r2 = rx * rx + ry * ry + rz * rz
                ! Add test and flag to find points within 0.18*R_in of atom ?
                if(r2 < rcut2) then
                   npoint(at)=npoint(at)+1
                   !ip_store(npoint(at),at)=ipoint
                   ip_store(ipoint,at) = npoint(at)
                   x_store(ipoint,at) = -rx
                   y_store(ipoint,at) = -ry
                   z_store(ipoint,at) = -rz
                   r_store(ipoint,at)=sqrt(r2)
                   !write(55,*) ipoint,at,r_store(ipoint,at)
                else
                   ip_store(ipoint,at) = 0
                   r_store(ipoint,at)=zero
                   x_store(ipoint,at)=zero
                   y_store(ipoint,at)=zero
                   z_store(ipoint,at)=zero
                end if  ! (r2 < rcut2) then
             end do! at = 1,natoms
          end do ! ix=1,nx_in_block
       end do ! iy=1,ny_in_block
    end do ! iz=1,nz_in_block
    return
  end subroutine check_block
  !!***

  !!****f* density_module/electron_number
  !! PURPOSE
  !!   calculates electron number. If optional spin component exsists
  !!   then calculate the electron number of the given spin component
  !! USAGE
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/05/29
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine electron_number(electrons)
    use datatypes
    use numbers
    use global_module, only: nspin
    use dimens,        only: grid_point_volume, n_my_grid_points
    use GenComms,      only: gsum
    use GenBlas,       only: rsum
    implicit none
    ! passed parameters
    real(double), dimension(nspin), intent(out) :: electrons
    ! local variables
    integer :: spin
    do spin = 1, nspin
       electrons(spin) = grid_point_volume * &
                         rsum(n_my_grid_points, density(:,spin), 1)
       call gsum(electrons(spin))
    end do
  end subroutine electron_number
  !!*****

end module density_module
