! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module pseudopotential_data
! ------------------------------------------------------------------------------
! Code area 10: pseudopotentials
! ------------------------------------------------------------------------------

!!****h* Conquest/pseudopotential_data
!!  NAME
!!   pseudopotential_data
!!  PURPOSE
!!   Contains all pseudopotential data and routines
!!  USES
!!   GenBlas, GenComms, atoms, block_module, common, cover_module, datatypes, dimens, 
!!   global_module, group_module, maxima_module, numbers, primary_module, set_blipgrid_module,
!!   species_module, splines
!!  AUTHOR
!!   IJB ?
!!  CREATION DATE
!! ?
!!  MODIFICATION HISTORY
!!   11/06/2001 dave
!!    Added ROBODoc header, RCS Id and Log tags and included
!!    all pseudopotential routines
!!   12/06/2001 dave
!!    Included pseudopotential_derivatives
!!   18/03/2002 dave
!!    Added RCS Id and Log tags and static tag for object file id
!!   24/06/2002 dave
!!    Changed raw open to io_assign and io_close from fdf
!!   12:21, 04/02/2003 drb 
!!    Various small changes
!!   09:39, 11/11/2004 dave 
!!    Removed N_ATOMS_MAX reference and changed to mx_icell
!!   2008/02/06 08:34 dave
!!    Changed for output to file not stdout
!!   2008/05/25 ast
!!    Added timers
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!  SOURCE
!!
module pseudopotential_data

  use datatypes
  use global_module,          only: io_lun
  use pseudopotential_common, only: pseudopotential, core_radius, &
                                    non_local
  use GenComms,               only: cq_abort
  use species_module,         only: non_local_species
  use timer_module,           only: start_timer, stop_timer, cq_timer
  use timer_module,           only: start_backtrace, stop_backtrace
  use timer_stdclocks_module, only: tmr_std_pseudopot, tmr_std_allocation

  implicit none

  save

  real(double), allocatable, dimension(:) :: core_radius_2, radius_max, radius_max_2, ps_exponent, &
       spherical_harmonic_norm, loc_radius
  real(double), allocatable, dimension(:,:) :: scale_functions, recip_scale, local_pseudopotential, &
       d2_local_pseudopotential
  real(double), allocatable, dimension(:,:,:) :: nl_pseudopotential, d2_nl_pseudopotential

  integer, allocatable, dimension(:) :: n_projectors, n_points_max, n_l_components
  integer, allocatable, dimension(:,:) :: l_core, l_component

!!***

contains

! -----------------------------------------------------------
! Subroutine init_pseudo
! -----------------------------------------------------------

!!****f* pseudopotential_data/init_pseudo *
!!
!!  NAME 
!!   init_pseudo
!!  USAGE
!! 
!!  PURPOSE
!!   Initialises pseudopotentials - reads, splines, gets core correction
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   11/06/2001 dave
!!  MODIFICATION HISTORY
!!   2008/05/25 ast
!!    Added timers
!!   2011/12/09 L.Tong
!!    Removed obsolete parameter number_of_bands
!!  SOURCE
!!
  subroutine init_pseudo(e_core)

    use datatypes
    use global_module,  only: IPRINT_TIME_THRES3
    use species_module, only: ps_file
    use timer_module
    
    ! Passed variables
    real(double)   :: e_core
    type(cq_timer) :: tmr_l_tmp1

    call start_timer(tmr_std_pseudopot)
    call start_timer(tmr_l_tmp1, WITH_LEVEL)
    call spline_pseudopotential()
    call set_pseudopotential()
    call get_core_correction(e_core)
    call stop_print_timer(tmr_l_tmp1, "initialising pseudopotentials", &
         IPRINT_TIME_THRES3)
    call stop_timer(tmr_std_pseudopot)
    return
  end subroutine init_pseudo
!!***

! -----------------------------------------------------------
! Subroutine set_pseudopotential
! -----------------------------------------------------------

!!****f* pseudopotential_data/set_pseudopotential *
!!
!!  NAME 
!!   set_pseudopotential
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
!!   E.H.Hernandez/D.R.Bowler/T.Miyazaki
!!  CREATION DATE
!!   15/03/95
!!  MODIFICATION HISTORY
!!   12/9/95 by Chris Goringe to use new data distribution
!!   22/2/96 by E. Hernandez to implement non-local
!!         Kleinman-Bylander pseudopotentials.
!!   17/4/96 by E. Hernandez to use block-data distribution
!!   21/5/96 by E. Hernandez for multiple species.
!!   08/2/01 by T. Miyazaki for new set_grid
!!               Loop structure has been changed.
!!   This subroutine has been rewritten to be consistent with new set_grid.
!!   I also change the loop structure.
!!  (Old)  do all_integration_grids in my domain
!!            do all_atoms in the sim. cell
!!             do l-component
!!               pseudofunctions
!!             enddo
!!               pseudo_density & pseudopotential
!!          .....
!!  (New)  do primary_blocks (my domain)
!!            do naba_atm
!!             do integration_grids in the block
!!               do l-component
!!                pseudofunctions
!!              (in the last component of l) pseudo_density & pot
!!               enddo
!!          .....
!!   25/05/2001 dave
!!    Stripped subroutine call, added ROBODoc header, indented
!!   11/06/2001 dave
!!    Added RCS Id and Log tags and GenComms
!!   12:21, 04/02/2003 drb 
!!    Small changes: improved definition of gauss and zeroed offset_position
!!   2008/03/03 18:37 dave
!!    Removed dsqrt, changed float to real
!!   2008/05/25 ast
!!    Added timers
!!   2011/12/12 17:26 dave
!!    Removed calculation of pseudofunctions for analytic blips
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks
!!  SOURCE
!!
  subroutine set_pseudopotential()

    use datatypes
    use GenBlas, only: axpy
    use numbers
    use global_module, only: rcellx,rcelly,rcellz,id_glob,ni_in_cell, &
                             species_glob, nlpf, flag_basis_set,      &
                             blips, flag_analytic_blip_int
    use species_module, only: charge, species, nlpf_species
    !  At present, these arrays are dummy arguments.
    use block_module, only : nx_in_block, ny_in_block, nz_in_block, &
         n_pts_in_block
    use group_module, only : blocks, parts
    use primary_module, only: domain
    use cover_module, only: DCS_parts
    use set_blipgrid_module, only : naba_atoms_of_blocks
    use GenComms, only: my_barrier, cq_abort, inode, ionode
    use hartree_module, only: hartree
    use functions_on_grid, only: gridfunctions, pseudofns
    use dimens, only: n_my_grid_points
    use maxima_module, only: maxngrid

    implicit none

    ! local
    real(double), allocatable, dimension(:) :: pseudo_density
    real(double):: dcellx_block,dcelly_block,dcellz_block
    real(double):: dcellx_grid, dcelly_grid, dcellz_grid
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom, stat
    real(double):: xatom,yatom,zatom,q,alpha,step,gauss
    real(double):: xblock,yblock,zblock
    integer :: the_species,l_comp_max,lcomp,the_ncf
    logical :: non_local_flag
    integer :: ix,iy,iz,j,iblock,the_l,ipoint,igrid
    real(double) :: dx,dy,dz,rx,ry,rz,r2,r_from_i
    real(double) :: rr,a,b,c,d,x,y,z,nl_potential,local_potential
    integer :: no_of_ib_ia, no_of_mcomp,offset_position
    integer :: position,iatom,icheck

    real(double) :: coulomb_potential(blocks%mx_ngonn*n_pts_in_block) !automatic array
    real(double) :: coulomb_energy
    ! --  Start of subroutine  ---

    call start_timer(tmr_std_allocation)
    allocate(pseudo_density(maxngrid),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating pseudo_density: ",maxngrid,stat)
    call stop_timer(tmr_std_allocation)

    pseudopotential = zero
    pseudo_density = zero
    if(flag_basis_set==blips.AND.(.NOT.flag_analytic_blip_int)) then
       gridfunctions(pseudofns)%griddata = zero
    end if


    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nx_in_block
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/ny_in_block
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nz_in_block

    ! loop arround grid points in this domain, and for each
    ! point, get contributions to the short-range
    ! potential and the core charge density from atoms which are 
    ! within the cutoff distance to that grid point

    call my_barrier()

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

                !--- I (TM) am not sure where the next line should be.
                !     (no_of_ib_ia = no_of_ib_ia +1)
                !    This determines how to store 
                !      pseudofunctions(ipoints,ncf,naba_atm,iblock).
                !    Now, I assume we should consider all naba_atm whose 
                !    distance from the block is within the maximum of core_radius.
                !    This is needed to keep the consistency with <set_bucket>.
                !    However, we can change this strategy by changing naba_atoms_of_blocks(nlpf).
                !no_of_ib_ia = no_of_ib_ia +1

                xatom=DCS_parts%xcover(icover)
                yatom=DCS_parts%ycover(icover)
                zatom=DCS_parts%zcover(icover)

                !the_species=species(ig_atom)
                the_species = species_glob(ig_atom)
                q =charge(the_species)
                alpha = ps_exponent(the_species)
                step = radius_max(the_species)/real(n_points_max(the_species)-1,double)
                non_local_flag = non_local_species(the_species)
                ! pseudofunction(n_pts_in_block,ncf,naba_atm,iblock)
                !  ncf is fixed for all atoms in this version,
                !  although it is not needed.
                !  Now, I keep the strategy of the old version  for ncf, 
                !  but I pack pseudofunction with respect to naba_atm
                !  by using species-dependent core radius.
                !    T. Miyazaki 8/Feb/2001
                if(non_local_flag.AND.flag_basis_set==blips) then
                   l_comp_max=n_l_components(the_species) ! s,p,d except local
                   the_ncf=nlpf_species(the_species)
                   !do lcomp=1,l_comp_max
                   !   the_l=l_component(lcomp,the_species)
                   !   the_ncf=the_ncf+2*the_l+1
                   !enddo
                else
                   l_comp_max=1
                   the_l = -1
                   the_ncf=0
                endif
                no_of_mcomp=0
                do lcomp=1,l_comp_max
                   if(non_local_flag) then
                      the_l=l_component(lcomp,the_species)
                      !offset_position = (no_of_ib_ia-1) * ncf * n_pts_in_block + &
                      !     no_of_mcomp * n_pts_in_block
                      offset_position = no_of_ib_ia * n_pts_in_block + &
                           no_of_mcomp * n_pts_in_block                     
                   else
                      offset_position = 0
                   endif
!!!!  CDIR NODEP            
                   icheck=0
                   !VEC  do ipoint=1,n_pts_in_block
                   !VEC    iz=(ipoint-1)/(nx_in_block*ny_in_block)+1
                   !VEC    iy=(ipoint-1-(iz-1)*nx_in_block*ny_in_block)/nx_in_block+1
                   !VEC    ix=ipoint-(iz-1)*nx_in_block*ny_in_block-(iy-1)*nx_in_block
                   ipoint=0
                   do iz=1,nz_in_block
                      do iy=1,ny_in_block
                         do ix=1,nx_in_block
                            ipoint=ipoint+1
                            igrid=n_pts_in_block*(iblock-1)+ipoint
                            if(.NOT.flag_analytic_blip_int) then
                               position= offset_position + ipoint
                               if(position > gridfunctions(pseudofns)%size) then
                                  call cq_abort('set_ps: position error ', &
                                       position, gridfunctions(pseudofns)%size)
                               endif
                            end if
                            if(igrid > n_my_grid_points) then
                               call cq_abort('set_ps: igrid error ', &
                                    igrid, n_my_grid_points)
                            endif
                            dx=dcellx_grid*(ix-1)
                            dy=dcelly_grid*(iy-1)
                            dz=dcellz_grid*(iz-1)
                            rx=xblock+dx-xatom
                            ry=yblock+dy-yatom
                            rz=zblock+dz-zatom
                            r2 = rx * rx + ry * ry + rz * rz
                            ! if within non-local cutoff, make the non-local
                            ! projections and store them in 'pseudofunctions'
                            if(r2 < core_radius_2(the_species)) then
                               icheck=1
                               r_from_i = sqrt( r2 )
                               j = floor( r_from_i / step ) + 1
                               ! check j (j+1 =< n_points_max)
                               if(j > n_points_max(the_species)) then
                                  call cq_abort('set_ps: overrun problem',j,n_points_max(the_species))
                               endif
                               ! Use the spline interpolation tables to build
                               ! the local part of the pseudopotl and 
                               ! pseudo-functions involved
                               ! in the non-local part on the grid
                               rr = real(j,double) * step
                               a = ( rr - r_from_i ) / step
                               b = one - a
                               c = a * ( a * a - one ) * step * step / six
                               d = b * ( b * b - one ) * step * step / six
                               !--- NONLOCAL PART START ---------
                               if((the_l >= 0).AND.(.NOT.flag_analytic_blip_int)) then  ! NONLOCAL case
                                  ! direction cosines needed for spher harms
                                  if ( r_from_i .gt. RD_ERR ) then
                                     x = rx / r_from_i
                                     y = ry / r_from_i
                                     z = rz / r_from_i
                                  else
                                     x = zero ; y = zero ; z = zero
                                  end if
                                  nl_potential =  &
                                       a * nl_pseudopotential(j,lcomp,the_species) + &
                                       b * nl_pseudopotential(j+1,lcomp,the_species) + &
                                       c * d2_nl_pseudopotential(j,lcomp,the_species) + &
                                       d * d2_nl_pseudopotential(j+1,lcomp,the_species)
                                  ! S, P or D components
                                  if ( the_l .eq. 0 ) then      ! NonLocal : S component
                                     gridfunctions(pseudofns)%griddata(position) =   &
                                          nl_potential * spherical_harmonic_norm(1) 

                                  else if ( the_l .eq. 1 ) then ! Nonlocal : P components
                                     gridfunctions(pseudofns)%griddata(position) =  & 
                                          nl_potential * spherical_harmonic_norm(2) * x 
                                     gridfunctions(pseudofns)%griddata(position+1*n_pts_in_block) =  &
                                          nl_potential * spherical_harmonic_norm(3) * y 
                                     gridfunctions(pseudofns)%griddata(position+2*n_pts_in_block) =   &
                                          nl_potential * spherical_harmonic_norm(4) * z 

                                  else if ( the_l .eq. 2 ) then ! Nonlocal : D components
                                     gridfunctions(pseudofns)%griddata(position) =  &
                                          nl_potential * spherical_harmonic_norm(5) * &
                                          (x*x - y*y ) 
                                     gridfunctions(pseudofns)%griddata(position+1*n_pts_in_block) = &
                                          nl_potential * spherical_harmonic_norm(6) * &
                                          ( three*z*z-one ) 
                                     gridfunctions(pseudofns)%griddata(position+2*n_pts_in_block) = &
                                          nl_potential * spherical_harmonic_norm(7) * &
                                          x * y 
                                     gridfunctions(pseudofns)%griddata(position+3*n_pts_in_block) = &
                                          nl_potential * spherical_harmonic_norm(8) * &
                                          x * z 
                                     gridfunctions(pseudofns)%griddata(position+4*n_pts_in_block) = &
                                          nl_potential * spherical_harmonic_norm(9) * &
                                          y * z 
                                  end if
                               endif
                               !--- NONLOCAL PART END ---------

                               !--- LOCAL PART START ---------
                               !  In the old version, the following tasks were performed
                               ! after the loop of l-components. The present version 
                               ! changes the order of loops, thus we have put this part below.
                               !  For the atoms which have non-local parts, the following
                               ! should be done once.
                               if(lcomp == l_comp_max) then 
                                  !construct the local part of the pseudopotential
                                  gauss = (alpha / pi)**1.5_double * exp(-alpha * r2)
                                  pseudo_density(igrid) = pseudo_density(igrid) + gauss * q

                                  local_potential =   &
                                       a * local_pseudopotential(j,the_species)+  &
                                       b * local_pseudopotential(j+1,the_species)+ &
                                       c * d2_local_pseudopotential(j,the_species)+  &
                                       d * d2_local_pseudopotential(j+1,the_species)

                                  pseudopotential(igrid) = pseudopotential(igrid) + &
                                       local_potential
                               endif
                               !--- LOCAL PART END ---------
                            endif ! if this point is naba points of the atom
                            ! To be consistent with the old version -----   

                         enddo !ix
                      enddo  !iy
                   enddo   !iz
                   !ORI enddo ! ipoints : points in the block

                   if(the_l >= 0) then
                      no_of_mcomp=no_of_mcomp+2*the_l+1
                   endif

                enddo ! l-components :

                no_of_ib_ia = no_of_ib_ia + the_ncf*n_pts_in_block
             enddo ! naba_atoms
          enddo ! naba_part
       endif !(naba_atoms_of_blocks(nlpf)%no_of_part(iblock) > 0) !naba atoms?
    enddo ! iblock : primary set of blocks

    ! now we must use FFT to transform the core charge density into
    ! reciprocal space to construct the long range part of the 
    ! pseudopotential in reciprocal space. This is then transformed
    ! back into real space and returned in coulomb_potential from
    ! subroutine hartree()
    call my_barrier()

    coulomb_potential = 0
    call hartree( pseudo_density, coulomb_potential, maxngrid, coulomb_energy )

    ! now add the two contributions together

    call axpy( n_my_grid_points, -one, coulomb_potential, 1, &
         pseudopotential, 1 )

    call start_timer(tmr_std_allocation)
    deallocate(pseudo_density, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating pseudo_density: ",n_my_grid_points,stat)
    call stop_timer(tmr_std_allocation)
    call my_barrier()
    return
  end subroutine set_pseudopotential
!!***

! -----------------------------------------------------------
! Subroutine pseudopotential_derivatives
! -----------------------------------------------------------

!!****f* pseudopotential_data/pseudopotential_derivatives *
!!
!!  NAME 
!!   pseudopotential_derivatives
!!  USAGE
!! 
!!  PURPOSE
!!   Calculates the derivatives of the non-local part
!!   of the pseudopotentials on the grid with respect to the 
!!   position of the atoms.
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   16/02/2001
!!  MODIFICATION HISTORY
!!   12/06/2001 dave
!!    Included in pseudopotential_data and added ROBODoc header
!!   2008/03/03 18:38 dave
!!    Removed dsqrt, changed float to real
!!   2008/05/25 ast
!!    Added timers
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks
!!  SOURCE
!!
  subroutine pseudopotential_derivatives(direction, dpseudofns )

    use datatypes
    use numbers
    use species_module, only: species, nlpf_species
    use global_module,  only: rcellx,rcelly,rcellz,id_glob,ni_in_cell, &
                              species_glob, nlpf
    !  At present, these arrays are dummy arguments.
    use block_module,   only : nx_in_block,ny_in_block,nz_in_block,    &
                               n_pts_in_block
    use group_module,   only : blocks, parts
    use primary_module, only: domain
    use cover_module,   only: DCS_parts
    use set_blipgrid_module, only : naba_atoms_of_blocks
    use functions_on_grid,   only: gridfunctions, fn_on_grid
    use dimens,              only: n_my_grid_points
    use GenComms,            only: my_barrier, cq_abort, inode,ionode

    implicit none
    ! dummy arguments
    integer,intent(in) :: direction
    integer :: dpseudofns
    ! local
    real(double):: dcellx_block,dcelly_block,dcellz_block
    real(double):: dcellx_grid, dcelly_grid, dcellz_grid
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom
    real(double):: xatom,yatom,zatom,step
    real(double):: xblock,yblock,zblock
    integer :: the_species,lcomp
    logical :: non_local_flag
    integer :: ix,iy,iz,j,iblock,the_l,ipoint,igrid
    real(double) :: dx,dy,dz,rx,ry,rz,r2,r_from_i
    real(double) :: rr,a,b,c,d,x,y,z,nl_potential
    real(double) :: nl_potential_derivative
    integer :: no_of_ib_ia, no_of_mcomp,offset_position, the_ncf
    integer :: position,iatom,icheck
    real(double) :: alpha,beta,gamma,delta

    call start_timer(tmr_std_pseudopot)
    gridfunctions(dpseudofns)%griddata = zero

    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nx_in_block
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/ny_in_block
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nz_in_block

    ! loop arround grid points in this domain, and for each
    ! point, get contributions to the short-range
    ! potential and the core charge density from atoms which are 
    ! within the cutoff distance to that grid point
    call my_barrier()
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
                call cq_abort("Error in pseudopotential derivatives: jpart overflow ",jpart,DCS_parts%mx_gcover)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atoms_of_blocks(nlpf)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atoms_of_blocks(nlpf)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)
                if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                   call cq_abort("Error in pseudopotential derivatives: atom overflow ", &
                        parts%icell_beg(ind_part),ni_in_cell)
                endif
                if(icover > DCS_parts%mx_mcover) then
                   call cq_abort("Error in pseudopotential derivatives: cover overflow ",&
                        icover,DCS_parts%mx_mcover)
                endif
                !--- I (TM) am not sure whether the next line should be here.
                !     (no_of_ib_ia = no_of_ib_ia +1)
                !    This determines how to store
                !      pseudofunctions(ipoints,ncf,naba_atm,iblock).
                !    Now, I assume we should consider all naba_atm whose
                !    distance from the block is within the maximum of core_radius.
                !    This is needed to keep the consistency with <set_bucket>.
                !    However, we can change this strategy by changing naba_atoms_of_blocks(nlpf).
                !no_of_ib_ia = no_of_ib_ia +1

                xatom=DCS_parts%xcover(icover)
                yatom=DCS_parts%ycover(icover)
                zatom=DCS_parts%zcover(icover)

                the_species = species_glob(ig_atom)
                !the_species=species(ig_atom)
                step = radius_max(the_species)/real(n_points_max(the_species)-1,double)
                non_local_flag = non_local_species(the_species)
                !the_ncf= ncf !!   This part should be changed when we consider
                !             !!   multiple ncf cases.
                ! pseudofunction(n_pts_in_block,ncf,naba_atm,iblock)
                !  ncf is fixed for all atoms in this version,
                !  although it is not needed.
                !  Now, I keep the strategy of the old version  for ncf, 
                !  but I pack pseudofunction with respect to naba_atm
                !  by using species-dependent core radius.
                !    T. Miyazaki 8/Feb/2001
                if(non_local_flag) then
                   the_ncf=nlpf_species(the_species)
                else
                   the_ncf=0
                end if
                if(non_local_flag) then
                   no_of_mcomp=0

                   do lcomp=1,n_l_components(the_species)
                      the_l=l_component(lcomp,the_species)
                      offset_position = no_of_ib_ia * n_pts_in_block  + no_of_mcomp * n_pts_in_block
!!!!  CDIR NODEP            
                      icheck=0
                      !ORI  do ipoint=1,n_pts_in_block
                      !ORI    iz=(ipoint-1)/(nx_in_block*ny_in_block)+1
                      !ORI    iy=(ipoint-1-(iz-1)*nx_in_block*ny_in_block)/nx_in_block+1
                      !ORI    ix=ipoint-(iz-1)*nx_in_block*ny_in_block-(iy-1)*nx_in_block
                      ipoint=0
                      do iz=1,nz_in_block
                         do iy=1,ny_in_block
                            do ix=1,nx_in_block
                               ipoint=ipoint+1

                               igrid=n_pts_in_block*(iblock-1)+ipoint
                               position= offset_position + ipoint

                               if(position > gridfunctions(dpseudofns)%size) then
                                  call cq_abort('pseudo_derivs: posn error ', &
                                       position, gridfunctions(dpseudofns)%size)
                               endif
                               if(igrid > n_my_grid_points) then
                                  call cq_abort('pseudo_derivs: grid error ', &
                                       igrid, n_my_grid_points)
                               endif
                               dx=dcellx_grid*(ix-1)
                               dy=dcelly_grid*(iy-1)
                               dz=dcellz_grid*(iz-1)
                               rx=xblock+dx-xatom
                               ry=yblock+dy-yatom
                               rz=zblock+dz-zatom
                               r2 = rx * rx + ry * ry + rz * rz
                               ! if we are within the non-local cutoff, build
                               ! the derivatives of the non-local projectors
                               if(r2 < core_radius_2(the_species)) then
                                  icheck=1

                                  r_from_i = sqrt( r2 )
                                  j = floor( r_from_i / step ) + 1

                                  ! OLD version : the following lines are
                                  if(j > n_points_max(the_species)-1) then
                                     call cq_abort('pseudo_derivs: j_error ',&
                                          j, n_points_max(the_species)-1)
                                  endif
                                  ! use splines to find derivatives
                                  rr = real(j,double) * step

                                  a = ( rr - r_from_i ) / step
                                  b = one - a
                                  c = a * ( a * a - one ) * step * step / six
                                  d = b * ( b * b - one ) * step * step / six

                                  alpha = -one / step
                                  beta =  one / step
                                  gamma = -step * ( three * a * a - one ) / six
                                  delta = step * ( three * b * b - one ) / six

                                  ! direction cosines for spherical harmonics
                                  if ( r_from_i .gt. RD_ERR ) then
                                     x = rx / r_from_i
                                     y = ry / r_from_i
                                     z = rz / r_from_i
                                  else
                                     x = zero ; y = zero ; z = zero
                                     r_from_i = RD_ERR
                                  end if

                                  nl_potential =  &
                                       a * nl_pseudopotential(j,lcomp,the_species) + &
                                       b * nl_pseudopotential(j+1,lcomp,the_species) + &
                                       c * d2_nl_pseudopotential(j,lcomp,the_species) + &
                                       d * d2_nl_pseudopotential(j+1,lcomp,the_species)

                                  nl_potential_derivative = &
                                       alpha * nl_pseudopotential(j,lcomp,the_species) + &
                                       beta * nl_pseudopotential(j+1,lcomp,the_species) + &
                                       gamma * d2_nl_pseudopotential(j,lcomp,the_species) + &
                                       delta * d2_nl_pseudopotential(j+1,lcomp,the_species)

                                  ! S, P or D components
                                  select case(the_l)
                                  case (0)     ! S component
                                     if(direction == 1) then
                                        gridfunctions(dpseudofns)%griddata(position) =   &
                                             nl_potential_derivative * spherical_harmonic_norm(1) * x
                                     elseif(direction == 2) then
                                        gridfunctions(dpseudofns)%griddata(position) =   &
                                             nl_potential_derivative * spherical_harmonic_norm(1) * y
                                     else
                                        gridfunctions(dpseudofns)%griddata(position) =   &
                                             nl_potential_derivative * spherical_harmonic_norm(1) * z
                                     endif

                                     !check
                                  case (1)     !  P components
                                     if(direction == 1) then
                                        !Px
                                        gridfunctions(dpseudofns)%griddata(position) =   &
                                             spherical_harmonic_norm(2) *   &
                                             ( nl_potential_derivative * x * x +   &
                                             nl_potential * ( one - x * x ) / r_from_i )
                                        !Py
                                        gridfunctions(dpseudofns)%griddata(position+1*n_pts_in_block) =   &
                                             spherical_harmonic_norm(3) *   &
                                             ( nl_potential_derivative * y * x -   &
                                             nl_potential * y * x / r_from_i )
                                        !Pz
                                        gridfunctions(dpseudofns)%griddata(position+2*n_pts_in_block) =   &
                                             spherical_harmonic_norm(4) *   &
                                             ( nl_potential_derivative * z * x -   &
                                             nl_potential * z * x / r_from_i )
                                     elseif(direction == 2) then
                                        !Px
                                        gridfunctions(dpseudofns)%griddata(position) =   &
                                             spherical_harmonic_norm(2) *   &
                                             ( nl_potential_derivative * x * y -   &
                                             nl_potential * x * y / r_from_i )
                                        !Py
                                        gridfunctions(dpseudofns)%griddata(position+1*n_pts_in_block) =   &
                                             spherical_harmonic_norm(3) *   &
                                             ( nl_potential_derivative * y * y +   &
                                             nl_potential * ( one - y * y ) / r_from_i )
                                        !Pz
                                        gridfunctions(dpseudofns)%griddata(position+2*n_pts_in_block) =   &
                                             spherical_harmonic_norm(4) *   &
                                             ( nl_potential_derivative * z * y -   &
                                             nl_potential * z * y / r_from_i )
                                     else
                                        !Px
                                        gridfunctions(dpseudofns)%griddata(position) =   &
                                             spherical_harmonic_norm(2) *   &
                                             ( nl_potential_derivative * x * z -   &
                                             nl_potential * x * z / r_from_i )
                                        !Py
                                        gridfunctions(dpseudofns)%griddata(position+1*n_pts_in_block) =   &
                                             spherical_harmonic_norm(3) *   &
                                             ( nl_potential_derivative * y * z -   &
                                             nl_potential * y * z / r_from_i )
                                        !Pz
                                        gridfunctions(dpseudofns)%griddata(position+2*n_pts_in_block) =   &
                                             spherical_harmonic_norm(4) *   &
                                             ( nl_potential_derivative * z * z +   &
                                             nl_potential * ( one - z * z ) / r_from_i )
                                     endif
                                  case (2)     !  D components
                                     if(direction == 1) then
                                        !D(x^2-y^2)
                                        gridfunctions(dpseudofns)%griddata(position) =                               &
                                             spherical_harmonic_norm(5) *                             &
                                             ( nl_potential_derivative * ( x * x - y * y ) * x +      &
                                             nl_potential * two * x *                               &
                                             ( one - x * x + y * y ) / r_from_i )
                                        !D(3z^2-1)
                                        gridfunctions(dpseudofns)%griddata(position+1*n_pts_in_block) =              &
                                             spherical_harmonic_norm(6) *                             &
                                             ( nl_potential_derivative * ( three * z * z - one ) * x -&
                                             nl_potential * six * z * z * x / r_from_i )
                                        !D(xy)
                                        gridfunctions(dpseudofns)%griddata(position+2*n_pts_in_block) =              &
                                             spherical_harmonic_norm(7) *                             &
                                             ( nl_potential_derivative * x * y * x +                  &
                                             nl_potential * y * ( one - two * x * x ) / r_from_i )

                                        !D(xz)
                                        gridfunctions(dpseudofns)%griddata(position+3*n_pts_in_block) =              &
                                             spherical_harmonic_norm(8) *                             &
                                             ( nl_potential_derivative * x * z * x +                  &
                                             nl_potential * z * ( one - two * x * x ) / r_from_i )
                                        !D(yz)
                                        gridfunctions(dpseudofns)%griddata(position+4*n_pts_in_block) =              &
                                             spherical_harmonic_norm(9) *                             &
                                             ( nl_potential_derivative * y * z * x -                  &
                                             nl_potential * two * x * z * y / r_from_i )
                                     elseif(direction == 2) then
                                        !D(x^2-y^2)
                                        gridfunctions(dpseudofns)%griddata(position) =                               &
                                             spherical_harmonic_norm(5) *                             &
                                             ( nl_potential_derivative * ( x * x - y * y ) * x +      &
                                             nl_potential * two * x *                               &
                                             ( one - x * x + y * y ) / r_from_i )
                                        !D(3z^2-1)
                                        gridfunctions(dpseudofns)%griddata(position+1*n_pts_in_block) =              &
                                             spherical_harmonic_norm(6) *                             &
                                             ( nl_potential_derivative * ( three * z * z - one ) * y -&
                                             nl_potential * six * z * z * y / r_from_i )
                                        !D(xy)
                                        gridfunctions(dpseudofns)%griddata(position+2*n_pts_in_block) =              &
                                             spherical_harmonic_norm(7) *                             &
                                             ( nl_potential_derivative * x * y * y +                  &
                                             nl_potential * x * ( one - two * y * y ) / r_from_i )
                                        !D(xz)
                                        gridfunctions(dpseudofns)%griddata(position+3*n_pts_in_block) =              &
                                             spherical_harmonic_norm(8) *                             &
                                             ( nl_potential_derivative * x * z * y -                  &
                                             nl_potential * two * x * z * y / r_from_i )
                                        !D(yz)
                                        gridfunctions(dpseudofns)%griddata(position+4*n_pts_in_block) =              &
                                             spherical_harmonic_norm(9) *                             &
                                             ( nl_potential_derivative * y * z * y +                  &
                                             nl_potential * z * ( one - two * y * y ) / r_from_i )
                                     else
                                        !D(x^2-y^2)
                                        gridfunctions(dpseudofns)%griddata(position) =                              &
                                             spherical_harmonic_norm(5) *                            &
                                             ( nl_potential_derivative * ( x * x - y * y ) * z +     &
                                             nl_potential * two * z *                              &
                                             ( y * y - x * x ) / r_from_i )
                                        !D(3z^2-1)
                                        gridfunctions(dpseudofns)%griddata(position+1*n_pts_in_block) =             &
                                             spherical_harmonic_norm(6) *                            &
                                             ( nl_potential_derivative * ( three * z * z - one ) * z - &
                                             nl_potential * six * z * ( one - z * z ) / r_from_i )
                                        !D(xy)
                                        gridfunctions(dpseudofns)%griddata(position+2*n_pts_in_block) =             &
                                             spherical_harmonic_norm(7) *                            &
                                             ( nl_potential_derivative * x * y * y +                 &
                                             nl_potential * x * ( one - two * y * y ) / r_from_i )
                                        !D(xz)
                                        gridfunctions(dpseudofns)%griddata(position+3*n_pts_in_block) =             &
                                             spherical_harmonic_norm(8) *                            &
                                             ( nl_potential_derivative * x * z * z +                 &
                                             nl_potential * x * ( one - two * z * z ) / r_from_i )
                                        !D(yz)
                                        gridfunctions(dpseudofns)%griddata(position+4*n_pts_in_block) =             &
                                             spherical_harmonic_norm(9) *                            &
                                             ( nl_potential_derivative * y * z * z +                 &
                                             nl_potential * y * ( one - two * z * z ) / r_from_i )
                                     endif
                                  end select

                               endif ! if this point is naba points of the atom

                            enddo !ix
                         enddo  !iy
                      enddo   !iz
                      !ORI enddo ! ipoints : points in the block

                      no_of_mcomp=no_of_mcomp+2*the_l+1

                   enddo ! l-components :

                endif ! if this atom has non-local parts
                no_of_ib_ia = no_of_ib_ia + the_ncf*n_pts_in_block
             enddo ! naba_atoms
          enddo ! naba_part
       endif !(naba_atoms_of_blocks(nlpf)%no_of_part(iblock) > 0) !naba atoms?
    enddo ! iblock : primary set of blocks

    call my_barrier()
    call stop_timer(tmr_std_pseudopot)
    return
  end subroutine pseudopotential_derivatives
!!***

! -----------------------------------------------------------
! Subroutine read_pseudopotential
! -----------------------------------------------------------

!!****f* pseudopotential_data/read_pseudopotential *
!!
!!  NAME 
!!   read_pseudopotential
!!  USAGE
!! 
!!  PURPOSE
!!   Reads in and distributes pseudopotentials
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   E.H.Hernandez and D.R.Bowler
!!  CREATION DATE
!!   23/1/96 EHH
!!  MODIFICATION HISTORY
!!   Frequent by DRB to correct non-local pseudopotentials
!!   28/05/2001 dave
!!    Stripped subroutine call
!!   08/06/2001 dave
!!    Changed exitpm to cq_abort, added abort if NCF too
!!    small and added RCS Id and Log tags
!!   11/06/2001 dave
!!    Included in pseudopotential_data
!!   24/06/2002 dave
!!    Changed raw open/close to io_assign and io_close from FDF
!!   12:20, 04/02/2003 drb 
!!    Small changes: reversed arguments to gcopy and initialised lcore to zero
!!   2008/05/25 ast
!!    Added timers
!!   2008/09/01 08:21 dave
!!    Added io_ routines from input_module
!!   2015/06/08 lat
!!    Added experimental backtrace
!!  SOURCE
!!
  subroutine read_pseudopotential(inode, ionode)

    use datatypes
    use numbers
    use global_module,  only: iprint_pseudo
    use species_module, only: n_species, nlpf_species, ps_file
    use GenComms,       only: gcopy, cq_abort
    use input_module,   only: io_assign, io_close

    implicit none

    ! Passed variables
    integer        :: inode, ionode

    ! Local variables
    type(cq_timer) :: backtrace_timer
    integer        :: i, ios, j, jnode
    integer        :: l, m, n_l, n, lun, nl_max, npts_max

    real(double), parameter :: ln10 = 2.302585092994_double
    real(double)            :: tmpr

!****lat<$
    call start_backtrace(t=backtrace_timer,who='read_pseudopotential',where=10,level=3)
!****lat>$

    call start_timer(tmr_std_pseudopot)
    ! loop over the different species of atoms in the simulation box,
    ! and for each species read the relevant pseudopotential data
    ! we only do this in the ionode; the remaining nodes receive the data
    ! from this one
    call start_timer(tmr_std_allocation)
    allocate(n_projectors(n_species) )
    allocate(n_points_max(n_species) )
    allocate(n_l_components(n_species) )
    allocate(core_radius(n_species) )
    allocate(core_radius_2(n_species) )
    allocate(loc_radius(n_species) )      
    allocate(radius_max(n_species) )      
    allocate(radius_max_2(n_species) )      
    allocate(ps_exponent(n_species) )      
    call stop_timer(tmr_std_allocation)
    n_projectors = 0
    if(inode==ionode) then
       nl_max = 0
       npts_max = 0
       do n=1, n_species
          call io_assign(lun)
          open( unit=lun, file=ps_file(n), iostat=ios )
          if ( ios .ne. 0 ) then
             write(io_lun,1) ps_file(n)
             call cq_abort('read_pseudopotential: file error')
          end if
          read(lun,*) n_points_max(n)
          if(n_points_max(n)>npts_max) npts_max = n_points_max(n)
          read(lun,*) loc_radius(n)
          read(lun,*) core_radius(n)
          read(lun,*) radius_max(n)
          ! then read the local part of the pseudopotential
          do i=1,n_points_max(n)
             read(lun,*) tmpr
          end do
          ! for consistency, the non-local pseudopotential data is read in
          ! a la CASTEP, the only difference being that we allow for the 
          ! data corresponding to different species to be read from 
          ! different files. The first value read is the l value. 
          ! n_componets keeps track of the total number of
          ! projections (a projection being specified by a value of l and m)
          ! for each species.
          if ( non_local_species(n) ) then
             read(lun,*) n_l_components(n)
             n_l = 0
             do l=1, n_l_components(n)
                read(lun,*) i
                n_l = n_l + 2*i+1
                read(lun,*) tmpr
                do i=1,n_points_max(n)
                   read(lun,*) tmpr
                end do
             end do
             if(n_l>nl_max) nl_max = n_l
          end if ! non_local
          call io_close(lun)
          ! now end loop over species
       end do ! species
    end if ! inode==ionode
    ! and copy the data just read in to the remaining nodes
    call gcopy(nl_max)
    call gcopy(npts_max)
    call gcopy(n_points_max, n_species )
    call gcopy(core_radius, n_species )
    call gcopy(loc_radius, n_species )      
    call gcopy(radius_max, n_species )      
    call gcopy(n_l_components, n_species )
    call start_timer(tmr_std_allocation)
    allocate(l_core(nl_max,n_species))
    allocate(l_component(nl_max, n_species ))
    allocate(recip_scale(nl_max, n_species ))
    allocate(scale_functions(nl_max, n_species ))
    allocate(local_pseudopotential(npts_max, n_species ))
    allocate(d2_local_pseudopotential(npts_max, n_species ))
    allocate(nl_pseudopotential(npts_max,nl_max,n_species))
    allocate(d2_nl_pseudopotential(npts_max,nl_max,n_species))
    call stop_timer(tmr_std_allocation)
    if(inode==ionode) then
       if(inode==ionode.AND.iprint_pseudo>=2) write(io_lun,2)
       do n=1, n_species
          call io_assign(lun)
          open( unit=lun, file=ps_file(n), iostat=ios )
          if ( ios .ne. 0 ) then
             write(io_lun,1) ps_file(n)
             call cq_abort('read_pseudopotential: file error')
          end if
          read(lun,*) i
          read(lun,*) i
          read(lun,*) i
          read(lun,*) i
          ! then read the local part of the pseudopotential
          read(lun,*) ( local_pseudopotential(i,n), i=1, n_points_max(n) )
          ! for consistency, the non-local pseudopotential data is read in
          ! a la CASTEP, the only difference being that we allow for the 
          ! data corresponding to different species to be read from 
          ! different files. The first value read is the l value. 
          ! n_componets keeps track of the total number of
          ! projections (a projection being specified by a value of l and m)
          ! for each species.
          if ( non_local_species(n) ) then
             read(lun,*) i
             do n_l=1, n_l_components(n)
                read(lun,*) l 
                l_component(n_l,n) = l
                ! this is the scaling value for componets of this angular momentum
                read(lun,*) scale_functions(n_l,n)
                ! and this is the pseudopotential data itself
                read(lun,*) (nl_pseudopotential(i,n_l,n),i=1,n_points_max(n))
             end do
          end if ! non_local
          call io_close(lun)
          ! now end loop over species
       end do ! species
    end if ! inode==ionode
    call gcopy( local_pseudopotential, npts_max, n_species )
    call gcopy( l_component, nl_max, n_species)
    call gcopy( nl_pseudopotential,npts_max,nl_max,n_species)
    call gcopy( scale_functions, nl_max, n_species )
    ! now we need to set the ps_exponent for the Gaussian pseudo-charge
    ! distribution for each species. This is done by requiring that the
    ! Gaussian falls off to a value of less than 10^-6 at the core radius.
    ! Also, set the auxiliary array recip_scale, which will be useful for
    ! constructing the non-local pseudopotentials on the grid
    l_core = 0
    do n=1, n_species
       core_radius_2(n) = core_radius(n) * core_radius(n)
       radius_max_2(n) = radius_max(n) * radius_max(n)
       ps_exponent(n) = six * ln10 / (loc_radius(n)*loc_radius(n))
       if ( non_local_species(n) ) then
          do j=1, n_l_components(n)
             recip_scale(j,n) = one / scale_functions(j,n) 
          end do
       end if
       n_l = 0
       do l=1, n_l_components(n)
          do m=1,2*l_component(l,n)+1
             n_l = n_l + 1
             l_core(n_l,n) = l
          end do
       end do
       n_projectors(n) = n_l
       nlpf_species(n) = n_l
       if (inode==ionode.AND.iprint_pseudo>=2) write(io_lun,3) n, ps_exponent(n)
    end do
    call stop_timer(tmr_std_pseudopot)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='read_pseudopotential')
!****lat>$

    return

1   format(10x,'An error occurred while trying to open file ', &
         a40,/,10x,'File does not exist or contains less data than needed')
2   format(2x,'Reading pseudopotentials'/)
3   format(2x,'Species ',i3,' Atomic exponent: ',f15.8)

  end subroutine read_pseudopotential
!!***

! -----------------------------------------------------------
! Subroutine spline_pseudopotential
! -----------------------------------------------------------

!!****f* pseudopotential_data/spline_pseudopotential *
!!
!!  NAME 
!!   spline_pseudopotential
!!  USAGE
!! 
!!  PURPOSE
!!   The task of this subroutine is to construct
!!   spline interpolation tables for the local and non-local
!!   parts of the pseudopotentials for the different species
!!   of atoms present in the system. This subroutine also
!!   calculates the values of the array spherical_harmonic_norm()
!!   containing the normalisation factors of the real functions
!!   resulting from the linear combination of spherical
!!   harmonics of equal l and |m|, which are needed for the
!!   construction of the non-local part of the pseudopotential
!!   on the integration grid.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   E.H.Hernandez/D.R.Bowler
!!  CREATION DATE
!!   01/04/96
!!  MODIFICATION HISTORY
!!   23/05/2001 dave
!!    F90, ROBODoc, indented
!!   11/06/2001 dave
!!    Included in pseudopotential_data
!!   2008/03/03 18:54 dave
!!    Changed float to real
!!   2008/05/25 ast
!!    Added timers
!!  SOURCE
!!
  subroutine spline_pseudopotential()

    use datatypes
    use numbers
    use species_module, only:  n_species, species
    use splines, only: spline

    implicit none

    ! Local variables

    integer :: i, j, n, n_l, stat

    real(double) :: d_end, d_origin, delta_r

    ! first of all construct the spherical_harmonic_norm() array

    call start_timer(tmr_std_allocation)
    allocate(spherical_harmonic_norm(9),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating spherical harmonic nrom: ",stat)
    call stop_timer(tmr_std_allocation)
    ! norm for the l=0, m=0
    spherical_harmonic_norm(1) = one / sqrt( four * pi )
    ! norm for the x, y and z type functions
    spherical_harmonic_norm(2) = sqrt( three / ( four * pi ) )
    spherical_harmonic_norm(3) = sqrt( three / ( four * pi ) )
    spherical_harmonic_norm(4) = sqrt( three / ( four * pi ) )
    ! norm for the x^2-y^2 and 3z^2-1 type functions
    spherical_harmonic_norm(5) = sqrt( ( three*five ) / ( four*four*pi ) )
    spherical_harmonic_norm(6) = sqrt( five / ( four * four * pi ) )
    ! norm for the x * y, x * z and y * z type functions
    spherical_harmonic_norm(7) = sqrt( ( three*five ) / ( four * pi ) )
    spherical_harmonic_norm(8) = sqrt( ( three*five ) / ( four * pi ) )
    spherical_harmonic_norm(9) = sqrt( ( three*five ) / ( four * pi ) )

    ! loop over species and do the interpolation
    do n=1, n_species
       ! do the splining for the local part of the pseudopotential
       delta_r = radius_max(n) / real( n_points_max(n) - 1, double)
       d_origin = ( local_pseudopotential(2,n) -  &
            local_pseudopotential(1,n) ) / delta_r
       d_end = ( local_pseudopotential(n_points_max(n),n) - &
            local_pseudopotential(n_points_max(n)-1,n) ) / &
            delta_r
       call spline( n_points_max(n), delta_r, local_pseudopotential(1,n),  &
            d_origin, d_end, d2_local_pseudopotential(1,n) )
       if ( non_local_species(n) ) then
          ! now calculate the second derivatives of these pseudopotentials for
          ! spline purposes; d_origin and d_end are the 1st derivatives at
          ! the origin and the end of the table respectively. These are 
          ! approximated in order to get continuous potentials
          delta_r = radius_max(n) / real( n_points_max(n) - 1, double)
          do n_l=1, n_l_components(n)
             d_origin = ( nl_pseudopotential(2,n_l,n) -  &
                  nl_pseudopotential(1,n_l,n) ) / delta_r
             d_end = ( nl_pseudopotential(n_points_max(n),n_l,n) - &
                  nl_pseudopotential(n_points_max(n)-1,n_l,n) ) / &
                  delta_r
             call spline( n_points_max(n), delta_r, &
                  nl_pseudopotential(1,n_l,n),  &
                  d_origin, d_end, d2_nl_pseudopotential(1,n_l,n) )
          end do
       end if
    end do
    return
  end subroutine spline_pseudopotential
!!***

! -----------------------------------------------------------
! Subroutine get_core_correction
! -----------------------------------------------------------

!!****f* pseudopotential_data/get_core_correction *
!!
!!  NAME 
!!   get_core_correction
!!  USAGE
!!   get_core_correction(number_of_bands, e_core)
!!  PURPOSE
!!   Calculates the correction to the pseudopotential
!!   energy due to the non-Coulomb nature of the
!!   pseudopotential in the core region. This 
!!   contribution is constant and thus needs to be
!!   calculated only once at the begining of the
!!   run. The correction arises only from the 
!!   Gaussian pseudo-charge that is used to construct
!!   the long-range Coulomb tail of the 
!!   pseudopotential (see subroutine set_pseudopotential())
!!   The short range part of the pseudopotential is
!!   dealt with in real space, and thus does not give
!!   rise to a correction.
!!  INPUTS
!!   real(double) :: e_core - energy calculated
!!   real(double) :: number_of_bands - obvious, really
!!  USES
!!   datatypes, numbers, dimens, atoms, species_module
!!  AUTHOR
!!   E.H.Hernandez
!!  CREATION DATE
!!   19/06/96
!!  MODIFICATION HISTORY
!!   21/05/2001 dave
!!    F90, ROBODoc header
!!   25/05/2001 dave
!!    Stripped subroutine call
!!   11/06/2001 dave
!!    Included in pseudopotential_data
!!   2011/12/09 L.Tong
!!    Removed obsolete parameter number_of_bands, replaced two *
!!    number_of_bands by ne_in_cell
!!  SOURCE
!!
  subroutine get_core_correction(e_core)

    use datatypes
    use numbers
    use dimens, only: volume
    use species_module, only: charge, species, n_species
    use global_module, only: ni_in_cell, ne_in_cell

    implicit none

    ! Passed variables
    real(double), intent(OUT) :: e_core

    ! Local variables
    integer :: i, n
    integer, allocatable, dimension(:) :: n_atoms_of_species

    ! first calculate the number of atoms of each species
    call start_timer (tmr_std_allocation)
    allocate (n_atoms_of_species(n_species))
    call stop_timer (tmr_std_allocation)
    n_atoms_of_species = 0
    do i=1, ni_in_cell
       n_atoms_of_species(species(i)) =  &
            n_atoms_of_species(species(i)) + 1
    end do

    ! loop over species and calculate the correction factor for each one
    e_core = zero
    do n=1, n_species
       e_core = e_core + real (n_atoms_of_species(n), double) * &
            charge(n) / ps_exponent(n)
    end do
    e_core = e_core * pi * ne_in_cell / volume
    call start_timer (tmr_std_allocation)
    deallocate (n_atoms_of_species)
    call stop_timer (tmr_std_allocation)
    return
  end subroutine get_core_correction
!!***

end module pseudopotential_data
