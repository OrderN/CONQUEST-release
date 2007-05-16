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
!!  SOURCE
module density_module
  
  use datatypes
  
  implicit none
  save

  logical :: flag_no_atomic_densities
  ! charge density
  real(double), allocatable, dimension(:) :: density ! values at gridpoints, to be calculated.

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
!!  SOURCE
!!
  subroutine set_density()

    use datatypes
    use numbers, ONLY: zero, one
    use global_module, ONLY: rcellx,rcelly,rcellz,id_glob,ni_in_cell, iprint_SC, species_glob, dens, ne_in_cell
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

    logical :: range_flag ! logical flag to warn if splint routine called out of the tabulated range.

    if(inode==ionode.AND.iprint_SC>=2) write(*,fmt='(2x,"Entering set_density")')

    density = zero  ! initialize density
    !write(*,*) 'Size of density: ',size(density)
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
       !write(*,*) 'Block ',iblock
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
    density = ne_in_cell*density/local_density
    if(inode.eq.ionode.AND.iprint_SC>0) write(*,*) 'In set_density, electrons: ',local_density
    call my_barrier()
    return
  end subroutine set_density
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
    use global_module, ONLY: iprint_SC, sf
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid

    implicit none

    ! Passed variables

    integer :: inode, ionode, size
    integer :: support, support_K
    real(double) :: electrons
    real(double), dimension(size) :: denout

    ! Local variables

    integer :: blk, i_count_alpha, n, n_i, n_point 

    if(inode==ionode.AND.iprint_SC>=2) write(*,fmt='(2x,"Entering get_electronic_density")')
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
    !      write(*,*) ' WARNING!!!   density < 0 ',n, size, denout(n)
    !   endif
    !enddo
    ! FOR DEBUGGING   T. MIYAZAKI 30/Jun/2003
    ! scale the density by two for spin

    call scal( size, two, denout, 1 )

    ! and calculate the electron number

    electrons = grid_point_volume * rsum( n_my_grid_points, denout, 1 )

    call gsum(electrons)
    if(inode.eq.ionode.AND.iprint_SC>1) write(*,*) 'Electrons: ',electrons

    ! support_K is using the same memory as h_on_support, so lets be safe
    ! and set it back to zero
    gridfunctions(support_K)%griddata = zero

    return
  end subroutine get_electronic_density
!!***

end module density_module


