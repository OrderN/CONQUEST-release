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
!!    moved declaration of density table variables to derived types in density_read_module.f90
!!   14/08/2002 drb and mjg
!!    changed naba table to supp from proj and added header
!!   15:21, 25/09/2002 mjg & drb 
!!    Tidied and updated to use dens tables from set_blipgrid, and changed to linear interpolation (for now)
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
!!  SOURCE
module density_module
  
  use datatypes
  use global_module, ONLY: io_lun
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_chargescf

  implicit none
  save

  logical :: flag_no_atomic_densities
  ! charge density
  real(double), allocatable, dimension(:) :: density     ! values at gridpoints, to be calculated.
  real(double), allocatable, dimension(:) :: density_pcc ! values at gridpoints, to be calculated for P.C.C.
  real(double) :: density_scale ! Scaling factor for atomic densities

  ! Becke weights
  real(double), allocatable, dimension(:) :: bw, bwgrid
  real(double), allocatable, dimension(:) :: atomcharge
  logical :: weights_set = .false.
  real(double), dimension(35) :: atrad
  ! RCS tag for object file identification
  character(len=80), private :: RCSid = "$Id$"

!!***

contains

! -----------------------------------------------------------
! Subroutine set_density
! -----------------------------------------------------------

!!****f* density_module/set_density *
!!
!!  NAME 
!!   set_density
!!  USAGE
!! 
!!  PURPOSE
!!   Puts the local charge density on the grid points belonging
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
!!    Added scaling to fix electron number (specifically to ensure correct charge in cell)
!!   2008/05/23 ast
!!    Added timers
!!  SOURCE
!!
  subroutine set_density()

    use datatypes
    use numbers, ONLY: zero, one
    use global_module, ONLY: rcellx,rcelly,rcellz,id_glob,ni_in_cell, iprint_SC, &
                             species_glob, dens, ne_in_cell, IPRINT_TIME_THRES3
    use block_module, ONLY : nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
    use group_module, ONLY : blocks, parts
    use primary_module, ONLY: domain
    use cover_module, ONLY: DCS_parts
    use set_blipgrid_module, ONLY : naba_atm
    use GenComms, ONLY: my_barrier, cq_abort, inode, ionode, gsum
    use atomic_density, ONLY: atomic_density_table
    use spline_module, ONLY: splint
    use dimens, ONLY: n_my_grid_points, grid_point_volume
    use GenBlas, ONLY: rsum, scal
    use timer_module

    implicit none

    ! Local Variables
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom
    integer :: the_species
    integer :: ix,iy,iz,j,iblock,ipoint,igrid
    integer :: no_of_ib_ia,offset_position
    integer :: position,iatom,icheck

    real(double):: dcellx_block,dcelly_block,dcellz_block
    real(double):: dcellx_grid, dcelly_grid, dcellz_grid
    real(double):: xatom,yatom,zatom,step, loc_cutoff
    real(double):: xblock,yblock,zblock, alpha
    real(double) :: dx,dy,dz,rx,ry,rz,r2,r_from_i
    real(double) :: local_density !local charge density returned from splint routine
    type(cq_timer) :: tmr_l_tmp1

    logical :: range_flag ! logical flag to warn if splint routine called out of the tabulated range.

    if(inode==ionode.AND.iprint_SC>=2) write(io_lun,fmt='(2x,"Entering set_density")')

    call start_timer(tmr_std_chargescf)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    density = zero  ! initialize density
    !write(io_lun,*) 'Size of density: ',size(density)
    !call scal(n_my_grid_points,zero,density,1)
    
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
       if(naba_atm(dens)%no_of_part(iblock) > 0) then ! if there are neighbour partitions
          iatom=0
          ! loop over neighbour partitions of this block
          do ipart=1,naba_atm(dens)%no_of_part(iblock)
             jpart=naba_atm(dens)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('set_ps: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             ! .... then atoms in this partition.
             do ia=1,naba_atm(dens)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atm(dens)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)

                if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                   call cq_abort('set_ps: globID ERROR ', ii,parts%icell_beg(ind_part))
                endif
                if(icover > DCS_parts%mx_mcover) then
                   call cq_abort('set_ps: icover ERROR ', icover,DCS_parts%mx_mcover)
                endif

                ! determine the position of the current atom
                xatom=DCS_parts%xcover(icover)
                yatom=DCS_parts%ycover(icover)
                zatom=DCS_parts%zcover(icover)

                ! determine the type of the current atom
                the_species=species_glob(ig_atom)
                loc_cutoff = atomic_density_table(the_species)%cutoff
                ! step in the density table
                step = loc_cutoff/real(atomic_density_table(the_species)%length-1,double)
                icheck=0
                ipoint=0
                ! loop over the grid points in the current block
                do iz=1,nz_in_block
                   do iy=1,ny_in_block
                      do ix=1,nx_in_block
                         ipoint=ipoint+1
                         igrid=n_pts_in_block*(iblock-1)+ipoint
                         ! position= offset_position + ipoint
                         if(igrid > n_my_grid_points) call cq_abort('set_density: igrid error ', igrid, n_my_grid_points)
                         dx=dcellx_grid*(ix-1)
                         dy=dcelly_grid*(iy-1)
                         dz=dcellz_grid*(iz-1)
                         ! determine separation between the current grid point and atom
                         rx=xblock+dx-xatom
                         ry=yblock+dy-yatom
                         rz=zblock+dz-zatom
                         r2 = rx * rx + ry * ry + rz * rz
                         if(r2 < loc_cutoff*loc_cutoff) then
                            r_from_i = sqrt(r2)
!                            j = aint(r_from_i/step) + 1
!                            ! check j (j+1 =< N_TAB_MAX)
!                            if(j+1 > atomic_density_table(the_species)%length) then
!                               call cq_abort('set_density: overrun problem',j)
!                            endif
!                            ! Linear interpolation for now
!                            alpha = (r_from_i/step)-real(j-1,double)
!                            local_density = (one-alpha)*atomic_density_table(the_species)%table(j) + &
!                                 alpha*atomic_density_table(the_species)%table(j+1)
                            call splint(step,atomic_density_table(the_species)%table(:), & 
                                 atomic_density_table(the_species)%d2_table(:), &
                                 atomic_density_table(the_species)%length, & 
                                 r_from_i,local_density,range_flag)
                            if(range_flag) call cq_abort('set_density: overrun problem')
                            ! recalculate the density for this grid point
                            density(igrid) = density(igrid) + local_density
                         endif ! if this point is within cutoff
                         ! test output of densities
!                         print '(4(f10.6,a))', xblock+dx, " *", yblock+dy, " *", zblock+dz, " *", density(igrid), " *"
                      enddo !ix gridpoints
                   enddo  !iy gridpoints
                enddo   !iz gridpoints
             enddo ! end loop over atoms in this partition
          enddo ! end loop over partitions
       endif ! end if there are neighbour atoms
    enddo ! end loop over blocks
    
    local_density = zero
    do iz=1,n_my_grid_points
       local_density = local_density + density(iz)
    end do
    local_density = local_density*grid_point_volume
    !local_density = grid_point_volume * rsum( n_my_grid_points, density, 1 )

    call gsum(local_density)
    ! Correct electron density
    density_scale = ne_in_cell/local_density
    density = density_scale*density
    if(inode.eq.ionode.AND.iprint_SC>0) write(io_lun,fmt='(10x,"In set_density, electrons: ",f20.12)') density_scale*local_density
    call my_barrier()
    call stop_print_timer(tmr_l_tmp1,"set_density",IPRINT_TIME_THRES3)
    call stop_timer(tmr_std_chargescf)
    return
  end subroutine set_density
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
!!  SOURCE
!!
  subroutine set_density_pcc()

    use datatypes
    use numbers, ONLY: zero, one
    use global_module, ONLY: rcellx,rcelly,rcellz,id_glob,ni_in_cell, iprint_SC, &
                             species_glob, dens, ne_in_cell, IPRINT_TIME_THRES3
    use block_module, ONLY : nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
    use group_module, ONLY : blocks, parts
    use primary_module, ONLY: domain
    use cover_module, ONLY: DCS_parts
    use set_blipgrid_module, ONLY : naba_atm
    use GenComms, ONLY: my_barrier, cq_abort, inode, ionode, gsum
    use pseudo_tm_info, ONLY : pseudo
    use spline_module, ONLY: splint
    use dimens, ONLY: n_my_grid_points, grid_point_volume
    use GenBlas, ONLY: rsum, scal
    use timer_module

    implicit none

    ! Local Variables
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom
    integer :: the_species
    integer :: ix,iy,iz,j,iblock,ipoint,igrid
    integer :: no_of_ib_ia,offset_position
    integer :: position,iatom,icheck

    real(double) :: dcellx_block,dcelly_block,dcellz_block
    real(double) :: dcellx_grid, dcelly_grid, dcellz_grid
    real(double) :: xatom,yatom,zatom
    real(double) :: pcc_cutoff, pcc_step
    real(double) :: xblock,yblock,zblock, alpha
    real(double) :: dx,dy,dz,rx,ry,rz,r2,r_from_i
    real(double) :: pcc_density !P.C.C. charge density returned from splint routine
    type(cq_timer) :: tmr_l_tmp1

    logical :: range_flag ! logical flag to warn if splint routine called out of the tabulated range.

    if(inode==ionode.AND.iprint_SC>=2) write(io_lun,fmt='(2x,"Entering set_density_pcc")')

    call start_timer(tmr_std_chargescf)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    density_pcc = zero  ! initialize density
    !write(io_lun,*) 'Size of density: ',size(density)
    !call scal(n_my_grid_points,zero,density,1)
    
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
       if(naba_atm(dens)%no_of_part(iblock) > 0) then ! if there are neighbour partitions
          iatom=0
          ! loop over neighbour partitions of this block
          do ipart=1,naba_atm(dens)%no_of_part(iblock)   !for now, dens is used even for P.C.C.
             jpart=naba_atm(dens)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('set_ps: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             ! .... then atoms in this partition.
             do ia=1,naba_atm(dens)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atm(dens)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)
                ! determine the type of the current atom
                the_species=species_glob(ig_atom)
                ! for P.C.C.
                if (.NOT. pseudo(the_species)%flag_pcc) then
                  cycle
                endif

                if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                   call cq_abort('set_ps: globID ERROR ', ii,parts%icell_beg(ind_part))
                endif
                if(icover > DCS_parts%mx_mcover) then
                   call cq_abort('set_ps: icover ERROR ', icover,DCS_parts%mx_mcover)
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
                         if(igrid > n_my_grid_points) call cq_abort('set_density: igrid error ', igrid, n_my_grid_points)
                         dx=dcellx_grid*(ix-1)
                         dy=dcelly_grid*(iy-1)
                         dz=dcellz_grid*(iz-1)
                         ! determine separation between the current grid point and atom
                         rx=xblock+dx-xatom
                         ry=yblock+dy-yatom
                         rz=zblock+dz-zatom
                         r2 = rx * rx + ry * ry + rz * rz
                         if(r2 < pcc_cutoff * pcc_cutoff) then
                            r_from_i = sqrt(r2)
!                            j = aint(r_from_i/step) + 1
!                            ! check j (j+1 =< N_TAB_MAX)
!                            if(j+1 > atomic_density_table(the_species)%length) then
!                               call cq_abort('set_density: overrun problem',j)
!                            endif
!                            ! Linear interpolation for now
!                            alpha = (r_from_i/step)-real(j-1,double)
!                            local_density = (one-alpha)*atomic_density_table(the_species)%table(j) + &
!                                 alpha*atomic_density_table(the_species)%table(j+1)
                            call splint( pcc_step,pseudo(the_species)%chpcc%f(:), & 
                                         pseudo(the_species)%chpcc%d2(:), &
                                         pseudo(the_species)%chpcc%n, & 
                                         r_from_i,pcc_density,range_flag)
                            if(range_flag) call cq_abort('set_density_pcc: overrun problem')
                            ! recalculate the density for this grid point
                            density_pcc(igrid) = density_pcc(igrid) + pcc_density
                         endif ! if this point is within cutoff
                         ! test output of densities
!                         print '(4(f10.6,a))', xblock+dx, " *", yblock+dy, " *", zblock+dz, " *", density(igrid), " *"
                      enddo !ix gridpoints
                   enddo  !iy gridpoints
                enddo   !iz gridpoints
             enddo ! end loop over atoms in this partition
          enddo ! end loop over partitions
       endif ! end if there are neighbour atoms
    enddo ! end loop over blocks
    
    !local_density = zero
    !do iz=1,n_my_grid_points
    !   local_density = local_density + density(iz)
    !end do
    !local_density = local_density*grid_point_volume
    !local_density = grid_point_volume * rsum( n_my_grid_points, density, 1 )

    !call gsum(local_density)
    ! Correct electron density
    !density_scale = ne_in_cell/local_density
    !density = density_scale*density
    !if(inode.eq.ionode.AND.iprint_SC>0) write(io_lun,fmt='(10x,"In set_density, electrons: ",f20.12)') density_scale*local_density
    call my_barrier()
    call stop_print_timer(tmr_l_tmp1,"set_density",IPRINT_TIME_THRES3)
    call stop_timer(tmr_std_chargescf)
    return
  end subroutine set_density_pcc
!!***

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
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2006/03/04 09:08 dave
!!    Changed to use new grid ideas, put into density_module
!!   2008/05/23 ast
!!    Added timers
!!  SOURCE
!!
  subroutine get_electronic_density(denout, electrons, support, support_K, inode, ionode, size)

    use datatypes
    use GenBlas, ONLY: scal, rsum
    use numbers, ONLY: zero, very_small, two
    use mult_module, ONLY: matK
    use dimens, ONLY: n_my_grid_points, grid_point_volume
    use block_module, ONLY: n_pts_in_block 
    use set_bucket_module,     ONLY: rem_bucket, sf_H_sf_rem
    use calc_matrix_elements_module, ONLY: act_on_vectors_new
    use primary_module, ONLY: domain
    use set_blipgrid_module, ONLY: naba_atm
    use GenComms, ONLY: gsum
    use global_module, ONLY: iprint_SC, sf, ni_in_cell, flag_Becke_weights
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid

    implicit none

    ! Passed variables

    integer :: inode, ionode, size
    integer :: support, support_K
    real(double) :: electrons
    real(double), dimension(size) :: denout

    ! Local variables

    integer :: blk, i_count_alpha, n, n_i, n_point 

    call start_timer(tmr_std_chargescf)
    if(inode==ionode.AND.iprint_SC>=2) write(io_lun,fmt='(2x,"Entering get_electronic_density")')
    gridfunctions(support_K)%griddata = zero
    call act_on_vectors_new(inode-1,rem_bucket(sf_H_sf_rem),matK,support_K,support)
    ! we can now calculate the density at each point

    denout = zero

    i_count_alpha=0
    do blk=1, domain%groups_on_node
       n_point = (blk - 1) * n_pts_in_block
       i_count_alpha = (naba_atm(sf)%ibegin_blk_orb(blk)-1)*n_pts_in_block
       !if(blk == 1) then
       !   i_count_alpha = 0
       !else
       !   i_count_alpha= i_count_alpha+ &
       !        naba_atm(sf)%no_of_atom(blk-1)*nsf*n_pts_in_block
       !endif
       if(naba_atm(sf)%no_of_atom(blk) > 0) then   !TM 30/Jun/2003
          !do n_i=1, naba_atm(sf)%no_of_atom(blk)*NSF*n_pts_in_block, &
          do n_i=1, naba_atm(sf)%no_of_orb(blk)*n_pts_in_block, n_pts_in_block
             do n=1, n_pts_in_block
                denout(n_point+n) = denout(n_point+n) + &
                     gridfunctions(support_K)%griddata(i_count_alpha+n_i+n-1) * &
                     gridfunctions(support)%griddata(i_count_alpha+n_i+n-1)
             end do
          end do
       endif !(naba_atm(sf)%no_of_atom(blk) > 0)   !TM 30/Jun/2003
    end do

    ! FOR DEBUGGING   T. MIYAZAKI 30/Jun/2003
    !do n=1, size
    !   if (denout(n) < -very_small) then
    !      write(io_lun,*) ' WARNING!!!   density < 0 ',n, size, denout(n)
    !   endif
    !enddo
    ! FOR DEBUGGING   T. MIYAZAKI 30/Jun/2003
    ! scale the density by two for spin

    call scal( size, two, denout, 1 )

    ! and calculate the electron number

    electrons = grid_point_volume * rsum( n_my_grid_points, denout, 1 )

    call gsum(electrons)
    if(inode.eq.ionode.AND.iprint_SC>1) write(io_lun,*) 'Electrons: ',electrons
    if(flag_Becke_weights) call build_Becke_charges(denout,size)
    ! support_K is using the same memory as h_on_support, so lets be safe
    ! and set it back to zero
    gridfunctions(support_K)%griddata = zero
    call stop_timer(tmr_std_chargescf)

    return
  end subroutine get_electronic_density
!!***

!!****f* / *
!!
!!  NAME 
!!   
!!  USAGE
!!   
!!  PURPOSE
!!   
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   
!!  CREATION DATE
!!   
!!  MODIFICATION HISTORY
!!  
!!  SOURCE
!!  
  subroutine build_Becke_weights!(chden,size)

    use datatypes
    use GenBlas, ONLY: scal, rsum
    use numbers, ONLY: zero, very_small, two, one, half
    use block_module, ONLY: n_pts_in_block 
    use primary_module, ONLY: domain
    use set_blipgrid_module, ONLY: naba_atm
    use GenComms, ONLY: gsum, cq_abort, inode, ionode
    use global_module, ONLY: rcellx,rcelly,rcellz,iprint_SC, sf, ni_in_cell, species_glob, id_glob, io_lun, &
         flag_Becke_atomic_radii
    use cover_module, ONLY: DCS_parts
    use group_module, ONLY : blocks, parts
    use dimens, ONLY: RadiusSupport, atomicrad
    use maxima_module, ONLY: maxngrid

    implicit none

    !! Passed variables
    !integer :: size
    !real(double), dimension(size) :: chden

    ! Local variables

    integer :: blk, no_of_ib_ia, n, n_i, ipos, jpos, stat, at, atomi, atomj, point
    integer, dimension(:), allocatable :: npoint, globatom
    integer, dimension(:,:), allocatable :: ip_store
    real(double) :: sum, mu, sij, sji, Rij
    real(double), dimension(:), allocatable :: xatom, yatom, zatom, rcut, rad
    real(double), dimension(:,:), allocatable :: r_store, x_store, y_store, z_store

    real(double):: dcellx_block,dcelly_block,dcellz_block
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom, iatom, iblock, the_species
    real(double):: xblock,yblock,zblock, aij, chi, uij, nu

    call start_timer(tmr_std_chargescf)
    if(inode==ionode.AND.iprint_SC>=2) write(io_lun,fmt='(2x,"Entering build_Becke_weights")')
    if(.NOT.weights_set) then
       call assign_atomic_radii
       weights_set = .true.
    end if
    no_of_ib_ia=0
    if(allocated(bw)) deallocate(bw)
    if(allocated(bwgrid)) deallocate(bwgrid)
    do blk=1, domain%groups_on_node
       no_of_ib_ia = no_of_ib_ia + naba_atm(sf)%no_of_atom(blk)
    end do
    no_of_ib_ia = no_of_ib_ia * n_pts_in_block
    allocate(bw(no_of_ib_ia),STAT=stat)
    allocate(bwgrid(maxngrid),STAT=stat)
    bw = -one
    no_of_ib_ia=0
    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz
    atomcharge = zero
    do blk=1, domain%groups_on_node
       !write(io_lun,*) 'Block: ', blk,naba_atm(sf)%no_of_atom(blk)
       xblock=(domain%idisp_primx(blk)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(blk)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(blk)+domain%nz_origin-1)*dcellz_block
       if(naba_atm(sf)%no_of_atom(blk) > 0) then  
          allocate(xatom(naba_atm(sf)%no_of_atom(blk)), yatom(naba_atm(sf)%no_of_atom(blk)), &
               zatom(naba_atm(sf)%no_of_atom(blk)), rcut(naba_atm(sf)%no_of_atom(blk)), &
               npoint(naba_atm(sf)%no_of_atom(blk)), globatom(naba_atm(sf)%no_of_atom(blk)), &
               ip_store(n_pts_in_block,naba_atm(sf)%no_of_atom(blk)), rad(naba_atm(sf)%no_of_atom(blk)),&
               x_store(n_pts_in_block,naba_atm(sf)%no_of_atom(blk)), &
               y_store(n_pts_in_block,naba_atm(sf)%no_of_atom(blk)), &
               z_store(n_pts_in_block,naba_atm(sf)%no_of_atom(blk)), &
               r_store(n_pts_in_block,naba_atm(sf)%no_of_atom(blk)),STAT=stat)
          ip_store = 0
          x_store = zero
          y_store = zero
          z_store = zero
          r_store = zero
          ! Replace this with part/atom loop !
          at = 0
          do ipart=1,naba_atm(sf)%no_of_part(blk)
             jpart=naba_atm(sf)%list_part(ipart,blk)
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('PAO_to_grid_global: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atm(sf)%no_atom_on_part(ipart,blk)
                at=at+1
                ii = naba_atm(sf)%list_atom(at,blk)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)
                xatom(at)=DCS_parts%xcover(icover)
                yatom(at)=DCS_parts%ycover(icover)
                zatom(at)=DCS_parts%zcover(icover)
                the_species=species_glob(ig_atom)
                globatom(at) = ig_atom
                rcut(at) = RadiusSupport(the_species)
                rad(at) = atomicrad(the_species)
             end do
          end do
          ! ** NEED rcut ** ! 
          !write(55,*) blk
          call check_block(xblock,yblock,zblock,xatom,yatom,zatom,rcut,&
               npoint,ip_store, x_store, y_store, z_store,r_store,n_pts_in_block,naba_atm(sf)%no_of_atom(blk))
          do point = 1, n_pts_in_block
             ! Sum is normalisation for grid point, Z = \sum_i bw_i
             sum = zero
             ! Actually we want a flag set in check_block so that this whole double loop
             ! is skipped if r_in criterion is met and the normalisation
             ! if(flag) then set weight for atom = 1, others = 0
             ! else
             do atomi = 1, naba_atm(sf)%no_of_atom(blk)
                ipos = no_of_ib_ia + n_pts_in_block*(atomi-1)
                if(ip_store(point,atomi)>0) then
                   ! Here we accumulate bw which is \prod_{j\ne i} s(mu_{ij})
                   if(naba_atm(sf)%no_of_atom(blk)==1) then
                      bw(ipos+point) = one
                   else
                      do atomj = atomi+1,naba_atm(sf)%no_of_atom(blk)
                         jpos = no_of_ib_ia + n_pts_in_block*(atomj-1)
                         if(ip_store(point,atomj)>0) then
                            ! Do we want a pre-calculated set of these ?
                            Rij = sqrt((xatom(atomi)-xatom(atomj))*(xatom(atomi)-xatom(atomj)) + &
                                 (yatom(atomi)-yatom(atomj))*(yatom(atomi)-yatom(atomj)) + &
                                 (zatom(atomi)-zatom(atomj))*(zatom(atomi)-zatom(atomj)))
                            ! make mu
                            mu = (r_store(point,atomi) - r_store(point,atomj))/Rij
                            if(flag_Becke_atomic_radii) then
                               chi = rad(atomi)/rad(atomj)
                               uij = (chi-one)/(chi+one)
                               aij = uij/(uij*uij-one)
                               if(aij<-half) aij = -half
                               if(aij>half) aij = half
                               !sij = s(mu)
                               nu = mu+aij*(1-mu*mu)
                               sij = s(nu)
                               sji = s(-nu)
                            else
                               sij = s(mu)
                               sji = s(-mu)
                            end if
                            if(bw(ipos+point)<zero) bw(ipos+point) = one
                            if(bw(jpos+point)<zero) bw(jpos+point) = one
                            bw(ipos+point) = bw(ipos+point)*sij
                            !write(54,fmt='(2i8,3f20.12)') ipos,point,mu,sij,bw(ipos+point)
                            bw(jpos+point) = bw(jpos+point)*sji
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
             !write(52,fmt='(2x,f20.12)') sum
             do atomi = 1, naba_atm(sf)%no_of_atom(blk)
                ipos = no_of_ib_ia + n_pts_in_block*(atomi-1)
                ! Normalise w
                if(abs(sum)>1.0e-8_double) bw(ipos+point) = bw(ipos+point)/sum
                ! ALEX: change constraint to appropriate array and uncomment lines below
                ! if(constraint(globatom(atomi))) bwgrid(point+(blk-1)*n_pts_in_block) = &
                !    bwgrid(point+(blk-1)*n_pts_in_block) + bw(ipos+point)
             end do ! atomi
             ! end if
          end do ! point
          no_of_ib_ia = no_of_ib_ia + n_pts_in_block*naba_atm(sf)%no_of_atom(blk)
          deallocate(xatom, yatom, zatom, rcut, npoint, globatom, ip_store, &
               x_store, y_store, z_store, r_store, rad, STAT=stat)
       endif !(naba_atm(sf)%no_of_atom(blk) > 0)   !TM 30/Jun/2003
    end do ! blk
    !!write(*,*) 'Done blocks'
    !atomcharge = atomcharge*grid_point_volume
    !call gsum(atomcharge,ni_in_cell)
    !!write(*,*) 'Done gsum'
    !if(inode==ionode) then
    !   do blk = 1, ni_in_cell
    !      write(io_lun,fmt='(2x,"Atom ",i4," Becke charge ",f20.12)') blk,atomcharge(blk)
    !   end do
    !end if
  end subroutine build_Becke_weights
!!***

!!****f* / *
!!
!!  NAME 
!!   
!!  USAGE
!!   
!!  PURPOSE
!!   
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   
!!  CREATION DATE
!!   
!!  MODIFICATION HISTORY
!!  
!!  SOURCE
!!  
  subroutine build_Becke_weight_forces(weight_force,chden,size)

    use datatypes
    use GenBlas, ONLY: scal, rsum
    use numbers, ONLY: zero, very_small, two, one
    use block_module, ONLY: n_pts_in_block 
    use primary_module, ONLY: domain
    use set_blipgrid_module, ONLY: naba_atm
    use GenComms, ONLY: gsum, cq_abort, inode, ionode
    use global_module, ONLY: rcellx,rcelly,rcellz,iprint_SC, sf, ni_in_cell, species_glob, id_glob, io_lun
    use cover_module, ONLY: DCS_parts
    use group_module, ONLY : blocks, parts
    use dimens, ONLY: RadiusSupport, grid_point_volume
    use maxima_module, ONLY: maxngrid

    implicit none

    ! Passed variables
    real(double), dimension(3,ni_in_cell) :: weight_force
    integer :: size
    real(double), dimension(size) :: chden

    ! Local variables

    integer :: blk, no_of_ib_ia, n, n_i, ipos, jpos, stat, at, atomi, atomj, point
    integer, dimension(:), allocatable :: npoint, globatom
    integer, dimension(:,:), allocatable :: ip_store
    real(double) :: sum, mu, sij, sji, Rij, dwjbydix, dwjbydiy, dwjbydiz, elec_here
    real(double), dimension(:), allocatable :: xatom, yatom, zatom, rcut
    real(double), dimension(:,:), allocatable :: r_store, x_store, y_store, z_store, &
         sum1x,sum2x,sum1y,sum2y,sum1z,sum2z,dmux,dmuy,dmuz,tij

    real(double):: dcellx_block,dcelly_block,dcellz_block
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom, iatom, iblock, the_species
    real(double):: xblock,yblock,zblock

    call start_timer(tmr_std_chargescf)
    if(inode==ionode.AND.iprint_SC>=2) write(io_lun,fmt='(2x,"Entering build_Becke_weights")')
    no_of_ib_ia=0
    if(allocated(bw)) deallocate(bw)
    if(allocated(bwgrid)) deallocate(bwgrid)
    do blk=1, domain%groups_on_node
       no_of_ib_ia = no_of_ib_ia + naba_atm(sf)%no_of_atom(blk)
    end do
    no_of_ib_ia = no_of_ib_ia * n_pts_in_block
    allocate(bw(no_of_ib_ia),STAT=stat)
    allocate(bwgrid(maxngrid),STAT=stat)
    bw = -one
    no_of_ib_ia=0
    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz
    atomcharge = zero
    do blk=1, domain%groups_on_node
       !write(io_lun,*) 'Block: ', blk,naba_atm(sf)%no_of_atom(blk)
       xblock=(domain%idisp_primx(blk)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(blk)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(blk)+domain%nz_origin-1)*dcellz_block
       if(naba_atm(sf)%no_of_atom(blk) > 0) then  
          allocate(xatom(naba_atm(sf)%no_of_atom(blk)), yatom(naba_atm(sf)%no_of_atom(blk)), &
               zatom(naba_atm(sf)%no_of_atom(blk)), rcut(naba_atm(sf)%no_of_atom(blk)), &
               npoint(naba_atm(sf)%no_of_atom(blk)), globatom(naba_atm(sf)%no_of_atom(blk)), &
               ip_store(n_pts_in_block,naba_atm(sf)%no_of_atom(blk)), &
               x_store(n_pts_in_block,naba_atm(sf)%no_of_atom(blk)), &
               y_store(n_pts_in_block,naba_atm(sf)%no_of_atom(blk)), &
               z_store(n_pts_in_block,naba_atm(sf)%no_of_atom(blk)), &
               r_store(n_pts_in_block,naba_atm(sf)%no_of_atom(blk)),STAT=stat)
          ip_store = 0
          r_store = zero
          ! Replace this with part/atom loop !
          if(naba_atm(sf)%no_of_atom(blk)>1) then ! If there's only one atom there's no force
             at = 0
             do ipart=1,naba_atm(sf)%no_of_part(blk)
                jpart=naba_atm(sf)%list_part(ipart,blk)
                if(jpart > DCS_parts%mx_gcover) then 
                   call cq_abort('PAO_to_grid_global: JPART ERROR ',ipart,jpart)
                endif
                ind_part=DCS_parts%lab_cell(jpart)
                do ia=1,naba_atm(sf)%no_atom_on_part(ipart,blk)
                   at=at+1
                   ii = naba_atm(sf)%list_atom(at,blk)
                   icover= DCS_parts%icover_ibeg(jpart)+ii-1
                   ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)
                   xatom(at)=DCS_parts%xcover(icover)
                   yatom(at)=DCS_parts%ycover(icover)
                   zatom(at)=DCS_parts%zcover(icover)
                   the_species=species_glob(ig_atom)
                   globatom(at) = ig_atom
                   rcut(at) = RadiusSupport(the_species)
                end do
             end do
             ! ** NEED rcut ** ! 
             !write(55,*) blk
             call check_block(xblock,yblock,zblock,xatom,yatom,zatom,rcut,&
                  npoint,ip_store, x_store, y_store, z_store,r_store,n_pts_in_block,naba_atm(sf)%no_of_atom(blk))
             allocate(sum1x(naba_atm(sf)%no_of_atom(blk),naba_atm(sf)%no_of_atom(blk)), &
                  sum1y(naba_atm(sf)%no_of_atom(blk),naba_atm(sf)%no_of_atom(blk)), &
                  sum1z(naba_atm(sf)%no_of_atom(blk),naba_atm(sf)%no_of_atom(blk)), &
                  sum2x(naba_atm(sf)%no_of_atom(blk),naba_atm(sf)%no_of_atom(blk)), &
                  sum2y(naba_atm(sf)%no_of_atom(blk),naba_atm(sf)%no_of_atom(blk)), &
                  sum2z(naba_atm(sf)%no_of_atom(blk),naba_atm(sf)%no_of_atom(blk)), &
                  dmux(naba_atm(sf)%no_of_atom(blk),naba_atm(sf)%no_of_atom(blk)), &
                  dmuy(naba_atm(sf)%no_of_atom(blk),naba_atm(sf)%no_of_atom(blk)), &
                  dmuz(naba_atm(sf)%no_of_atom(blk),naba_atm(sf)%no_of_atom(blk)), &
                  tij(naba_atm(sf)%no_of_atom(blk),naba_atm(sf)%no_of_atom(blk)))
             do point = 1, n_pts_in_block
                elec_here = density(n_pts_in_block*(blk-1)+point) * grid_point_volume
                ! Actually we want a flag set in check_block so that this whole double loop
                ! is skipped if r_in criterion is met and the normalisation
                ! if(flag) then set weight for atom = 1, others = 0
                ! else
                sum1x = zero
                sum1y = zero
                sum1z = zero
                sum2x = zero
                sum2y = zero
                sum2z = zero
                dmux = zero
                dmuy = zero
                dmuz = zero
                tij = zero
                ! First time around we need to store the sums over t and \nabla \mu
                ! as well as storing t and nabla
                do atomi = 1, naba_atm(sf)%no_of_atom(blk)
                   ipos = no_of_ib_ia + n_pts_in_block*(atomi-1)
                   if(ip_store(point,atomi)>0) then
                      ! Here we accumulate bw which is \prod_{j\ne i} s(mu_{ij})
                      do atomj = 1,naba_atm(sf)%no_of_atom(blk)
                         if(atomj/=atomi) then
                            jpos = no_of_ib_ia + n_pts_in_block*(atomj-1)
                            if(ip_store(point,atomj)>0) then
                               ! Do we want a pre-calculated set of these ?
                               Rij = sqrt((xatom(atomi)-xatom(atomj))*(xatom(atomi)-xatom(atomj)) + &
                                    (yatom(atomi)-yatom(atomj))*(yatom(atomi)-yatom(atomj)) + &
                                    (zatom(atomi)-zatom(atomj))*(zatom(atomi)-zatom(atomj)))
                               ! make mu
                               mu = (r_store(point,atomi) - r_store(point,atomj))/Rij
                               dmux(atomi,atomj) = x_store(point,atomi)/Rij - (xatom(atomi)-xatom(atomj))*mu/Rij
                               dmuy(atomi,atomj) = y_store(point,atomi)/Rij - (yatom(atomi)-yatom(atomj))*mu/Rij
                               dmuz(atomi,atomj) = z_store(point,atomi)/Rij - (zatom(atomi)-zatom(atomj))*mu/Rij
                               tij(atomi,atomj) = t(mu)
                               tij(atomj,atomi) = t(-mu)
                               sum1x(atomi,atomj) = sum1x(atomi,atomj) + tij(atomi,atomj)*dmux(atomi,atomj)
                               sum1y(atomi,atomj) = sum1y(atomi,atomj) + tij(atomi,atomj)*dmuy(atomi,atomj)
                               sum1z(atomi,atomj) = sum1z(atomi,atomj) + tij(atomi,atomj)*dmuz(atomi,atomj)
                               sum2x(atomi,atomj) = sum2x(atomi,atomj) + bw(jpos+point)*tij(atomj,atomi)*dmux(atomi,atomj)
                               sum2y(atomi,atomj) = sum2y(atomi,atomj) + bw(jpos+point)*tij(atomj,atomi)*dmuy(atomi,atomj)
                               sum2z(atomi,atomj) = sum2z(atomi,atomj) + bw(jpos+point)*tij(atomj,atomi)*dmuz(atomi,atomj)
                            end if
                         end if
                      end do ! atomj
                      if(bw(ipos+point)<zero) bw(ipos+point) = zero
                   end if
                end do ! atomi
                ! Now accumulate \nabla_i w_j
                do atomi = 1, naba_atm(sf)%no_of_atom(blk)
                   ipos = no_of_ib_ia + n_pts_in_block*(atomi-1)
                   if(ip_store(point,atomi)>0) then
                      ! Here we accumulate bw which is \prod_{j\ne i} s(mu_{ij})
                      do atomj = 1,naba_atm(sf)%no_of_atom(blk)
                         !if(constraint(globatom(j))) then
                         jpos = no_of_ib_ia + n_pts_in_block*(atomj-1)
                         if(ip_store(point,atomj)>0) then
                            dwjbydix = -bw(jpos+point)*tij(atomj,atomi)*dmux(atomi,atomj) - &
                                 bw(jpos+point)*(bw(ipos+point)*sum1x(atomi,atomj) - sum2x(atomi,atomj))
                            dwjbydiy = -bw(jpos+point)*tij(atomj,atomi)*dmuy(atomi,atomj) - &
                                 bw(jpos+point)*(bw(ipos+point)*sum1y(atomi,atomj) - sum2y(atomi,atomj))
                            dwjbydiz = -bw(jpos+point)*tij(atomj,atomi)*dmuz(atomi,atomj) - &
                                 bw(jpos+point)*(bw(ipos+point)*sum1z(atomi,atomj) - sum2z(atomi,atomj))
                            weight_force(1,globatom(atomi)) = weight_force(1,globatom(atomi)) + elec_here * dwjbydix
                            weight_force(2,globatom(atomi)) = weight_force(2,globatom(atomi)) + elec_here * dwjbydiy
                            weight_force(3,globatom(atomi)) = weight_force(3,globatom(atomi)) + elec_here * dwjbydiz
                         end if
                         !end if
                      end do
                   end if
                end do
             end do ! point
          end if
          no_of_ib_ia = no_of_ib_ia + n_pts_in_block*naba_atm(sf)%no_of_atom(blk)
          deallocate(xatom, yatom, zatom, rcut, npoint, globatom, ip_store, r_store,STAT=stat)
       endif !(naba_atm(sf)%no_of_atom(blk) > 0)   !TM 30/Jun/2003
    end do ! blk
    call gsum(weight_force,3,ni_in_cell)
  end subroutine build_Becke_weight_forces
!!***

  subroutine build_Becke_weight_matrix

    use datatypes
    use numbers
    use global_module, ONLY: sf
    use block_module, ONLY: n_blocks, n_pts_in_block
    use primary_module, ONLY: domain
    use set_blipgrid_module, ONLY: naba_atm
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid, H_on_supportfns, supportfns
    use calc_matrix_elements_module, ONLY: norb
    use set_bucket_module, ONLY: rem_bucket, sf_H_sf_rem, pao_H_sf_rem
    use calc_matrix_elements_module, ONLY: get_matrix_elements_new
    use GenComms, ONLY: inode, ionode

    implicit none

    integer :: n, m, nb, atom, nsf1, point, matW

    ! Create w_c(r)|chi_j>
    gridfunctions(H_on_supportfns)%griddata = zero
    n = 0
    m = 0
    do nb=1,domain%groups_on_node
       if(naba_atm(sf)%no_of_atom(nb)>0) then
          do atom = 1, naba_atm(sf)%no_of_atom(nb)
             do nsf1 = 1, norb(naba_atm(sf),atom,nb)
                do point = 1, n_pts_in_block
                   n = n + 1
                   gridfunctions(H_on_supportfns)%griddata(n) = &
                        gridfunctions(supportfns)%griddata(n)*bwgrid(m+point)
                end do
             end do
          end do
       end if ! (naba_atm(sf)%no_of_atom(nb)>0)
       m = m + n_pts_in_block
    end do
    call get_matrix_elements_new(inode-1,rem_bucket(sf_H_sf_rem),matW,supportfns,H_on_supportfns)

  end subroutine build_Becke_weight_matrix

  subroutine build_Becke_charges(chden,size)

    use datatypes
    use GenComms, ONLY: gsum, cq_abort, inode, ionode
    use numbers, ONLY: zero, very_small, two, one
    use global_module, ONLY: sf, ni_in_cell, id_glob, io_lun
    use primary_module, ONLY: domain
    use set_blipgrid_module, ONLY: naba_atm
    use cover_module, ONLY: DCS_parts
    use dimens, ONLY: grid_point_volume
    use block_module, ONLY: n_pts_in_block 
    use group_module, ONLY : parts

    implicit none

    ! Passed variables
    integer :: size
    real(double), dimension(size) :: chden

    ! Local variables

    integer :: blk, no_of_ib_ia, at, point, ipos
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom, iatom, iblock, the_species

    atomcharge = zero
    no_of_ib_ia=0
    do blk=1, domain%groups_on_node
       if(naba_atm(sf)%no_of_atom(blk) > 0) then  
          at = 0
          do ipart=1,naba_atm(sf)%no_of_part(blk)
             jpart=naba_atm(sf)%list_part(ipart,blk)
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('PAO_to_grid_global: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atm(sf)%no_atom_on_part(ipart,blk)
                at=at+1
                ii = naba_atm(sf)%list_atom(at,blk)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)
                do point = 1, n_pts_in_block
                   ipos = no_of_ib_ia + n_pts_in_block*(at-1)
                   atomcharge(ig_atom) = atomcharge(ig_atom) + &
                        bw(ipos+point)*chden(point+(blk-1)*n_pts_in_block)
                end do ! point
             end do ! ia
          end do ! ipart
          no_of_ib_ia = no_of_ib_ia + n_pts_in_block*naba_atm(sf)%no_of_atom(blk)
       endif !(naba_atm(sf)%no_of_atom(blk) > 0)   !TM 30/Jun/2003
    end do ! blk
    !write(*,*) 'Done blocks'
    atomcharge = atomcharge*grid_point_volume
    call gsum(atomcharge,ni_in_cell)
    if(inode==ionode) then
       do blk = 1, ni_in_cell
          write(io_lun,fmt='(2x,"Atom ",i4," Becke charge ",f20.12)') blk,atomcharge(blk)
       end do
    end if
  end subroutine build_Becke_charges

!!****f* / *
!!
!!  NAME 
!!   
!!  USAGE
!!   
!!  PURPOSE
!!   
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   
!!  CREATION DATE
!!   
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
   
  real(double) function t(mu)

    use numbers, only: zero, half, one, three, five, seven

    implicit none

    real(double) :: mu, mua, mua2, s
    real(double), parameter :: a = 0.64_double
    real(double), parameter :: reca = 1.5625_double
    real(double), parameter :: c1 = 2.1875_double
    real(double), parameter :: c3 = 2.1875_double
    real(double), parameter :: c5 = 1.3125_double
    real(double), parameter :: c7 = 0.3125_double

    if(mu<-a) then
       t = zero
    else if(mu>a) then
       t = zero
    else
       mua = mu/a
       mua2 = mua*mua
       t = -half*(reca*c1 - mua2*(three*reca*c3 - mua2*(five*reca*c5 - seven*reca*c7*mua2)))
       s = half*(one -  mua*(c1 - mua2*(c3 - mua2*(c5 - c7*mua2))))
       t = t/s
    end if
  end function t

  ! Assigns atomic radii; these are from J. C. Slater, J. Chem. Phys. 41, 3199 (1964)
  ! H altered to 0.35 A, F to 0.9 A
  subroutine assign_atomic_radii

    use datatypes
    use units, ONLY: AngToBohr

    implicit none

    atrad(1)=0.35_double*AngToBohr

    atrad(3)=1.45_double*AngToBohr
    atrad(4)=1.05_double*AngToBohr
    atrad(5)=0.85_double*AngToBohr

    atrad(6)=0.70_double*AngToBohr
    atrad(7)=0.65_double*AngToBohr
    atrad(8)=0.60_double*AngToBohr
    atrad(9)=0.90_double*AngToBohr


    atrad(11)=1.80_double*AngToBohr
    atrad(12)=1.50_double*AngToBohr
    atrad(13)=1.25_double*AngToBohr
    atrad(14)=1.10_double*AngToBohr
    atrad(15)=1.00_double*AngToBohr
    atrad(16)=1.00_double*AngToBohr
    atrad(17)=1.00_double*AngToBohr


    atrad(32)=1.25_double*AngToBohr
    atrad(33)=1.15_double*AngToBohr
    atrad(34)=1.15_double*AngToBohr
    atrad(35)=1.15_double*AngToBohr

  end subroutine assign_atomic_radii

  subroutine check_block &
       (xblock,yblock,zblock,xatom,yatom,zatom,rcut, &
       npoint, ip_store, x_store, y_store, z_store, r_store,blocksize,natoms) 

    use numbers
    use global_module, ONLY: rcellx,rcelly,rcellz
    use group_module,  ONLY: blocks
    use block_module,  ONLY: nx_in_block,ny_in_block,nz_in_block!, &
!         n_pts_in_block


    implicit none
    !Passed 
    integer :: blocksize, natoms
    real(double), intent(in):: xblock, yblock, zblock
    real(double), intent(in), dimension(natoms):: xatom, yatom, zatom, rcut
    integer, intent(out) :: npoint(natoms), ip_store(blocksize,natoms)
    real(double), dimension(blocksize,natoms), intent(out) :: r_store, x_store, y_store, z_store
    !Local
    real(double):: dcellx_block,dcelly_block,dcellz_block
    real(double):: dcellx_grid, dcelly_grid, dcellz_grid
    real(double):: dx, dy, dz
    integer :: ipoint, iz, iy, ix, at
    real(double) ::  r2, r_from_i, rx, ry, rz, x, y, z, rcut2


    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nx_in_block
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/ny_in_block
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nz_in_block

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
                   x_store(ipoint,at) = rx
                   y_store(ipoint,at) = ry
                   z_store(ipoint,at) = rz
                   r_store(ipoint,at)=sqrt(r2)
                   !write(55,*) ipoint,at,r_store(ipoint,at)
                else
                   ip_store(ipoint,at) = 0
                   r_store(ipoint,at)=zero
                endif  ! (r2 < rcut2) then
             end do! at = 1,natoms
          enddo ! ix=1,nx_in_block
       enddo ! iy=1,ny_in_block
    enddo ! iz=1,nz_in_block
    return
  end subroutine check_block
!!***

end module density_module


