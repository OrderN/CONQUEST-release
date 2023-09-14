! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module PAO_grid_transform_module
! ------------------------------------------------------------------------------
! Code area 11: basis operations
! ------------------------------------------------------------------------------

!!****h* Conquest/PAO_grid_transform_module *
!!  NAME
!!   PAO_grid_transform_module
!!  PURPOSE
!!   Transforms support functions represented in a PAO basis onto integration grid
!!
!!   This follows rather closely the blip_grid_transform_module
!!  USES
!!
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   16:37, 2003/09/22
!!  MODIFICATION HISTORY
!!   2005/07/11 10:18 dave
!!    Tidied up use statements (only one occurence) for ifort
!!   10:31, 13/02/2006 drb
!!    Changed lines to conform to F90 standard
!!   2006/06/20 08:15 dave
!!    Various changes for variable NSF
!!   2007/05/11 07:56 dave
!!    Tidying up timing calls throughout
!!   2008/02/01 17:45 dave
!!    Changes to write output to file not stdout
!!   2008/05/28 ast
!!    Added timers
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2016/11/10 16:15 nakata
!!    Added subroutine "single_PAO_to_grad"
!!   2017/03/08 15:00 nakata
!!    Removed PAO_to_grid, do_PAO_transform, PAO_to_grad, do_PAO_grad_transform,
!!    PAO_to_grid_global and PAO_to_grad_global, which are no longer used.
!!  SOURCE
!!
module PAO_grid_transform_module

  use global_module,          only: io_lun
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_basis, tmr_std_allocation

  implicit none

!!***

contains

!!****f* PAO_grid_transform_module/PAO_or_gradPAO_to_grid *
!!
!!  NAME
!!   PAO_or_gradPAO_to_grid
!!  USAGE
!!
!!  PURPOSE
!!   Projects the PAO functions onto the grid (rather than support functions)
!!   Used for gradients of energy wrt PAO coefficients
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2004 sometime, probably
!!  MODIFICATION HISTORY
!!   2006/06/20 08:28 dave
!!    Various changes for variable NSF
!!   2008/05/28 ast
!!    Added timers
!!   2009/07/08 16:48 dave
!!    Added code for one-to-one PAO to SF assignment
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks
!!   2016/07/29 18:30 nakata
!!    Renamed supports_on_atom -> blips_on_atom
!!   2016/08/01 18:30 nakata
!!    Removed unused sf
!!   2016/08/09 14:00 nakata
!!    Renamed support -> pao_fns
!!   2016/09/30 18:30 nakata
!!    Changed paof to atomf (atomf is sf for one_to_one PAOs but quivalent to paof)
!!   2016/12/29 19:30 nakata
!!    Removed support_spec_format (blips_on_atom and flag_one_to_one) and this_atom,
!!    which are no longer needed
!!    Removed unused npao1 and atom_species
!!   2023/08/23 11:20 tkoskela
!!    Added OMP threading, merged single_PAO_to_grad and single_PAO_to_grid into
!!    single subroutine PAO_or_gradPAO_to_grid
!!  SOURCE
!!
  subroutine PAO_or_gradPAO_to_grid(pao_fns, evaluate, direction)

    use datatypes, only: double
    use primary_module, ONLY: bundle
    use GenComms, ONLY: my_barrier, cq_abort, mtime
    use dimens, ONLY: r_h
    use GenComms, ONLY: inode, ionode
    use numbers
    use global_module, ONLY: rcellx,rcelly,rcellz,id_glob,ni_in_cell, iprint_basis, species_glob, atomf
    use species_module, ONLY: species, npao_species
    !  At present, these arrays are dummy arguments.
    use block_module, ONLY : nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
    use group_module, ONLY : blocks
    use primary_module, ONLY: domain
    use cover_module, ONLY: DCS_parts
    use set_blipgrid_module, ONLY : naba_atoms_of_blocks
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid
    use pao_format

    implicit none
    integer, intent(in) :: pao_fns
    integer, intent(in) :: direction

    !local
    integer :: iblock,ia,ipart,ip,l1,acz,m1 ! Loop index variables
    integer :: nblock, npart, natom ! Array dimensions
    integer :: naba_part_label,ind_part,icover,iatom ! indices for calculating species. TODO:Need better names
    integer :: position,next_offset_position ! Temporary indices to gridfunctions%griddata
    integer :: my_species ! Temporary variables to reduce indirect accesses
    integer :: npoint ! outputs of check_block
    integer :: count1 ! incremented counter, maps from (l1, acz, m1) to linear index of gridfunctions%griddata
    real(double):: dcellx_block,dcelly_block,dcellz_block ! grid dimensions, should be moved
    real(double) :: x,y,z ! Temporary variables to reduce indirect accesses
    real(double) :: rcut ! Input to check_block
    real(double) :: val ! output, written into gridfunctions%griddata
    real(double) :: xblock,yblock,zblock ! inputs to check_block
    real(double) :: xatom,yatom,zatom ! inputs to check_block
    integer, allocatable, dimension(:) :: ip_store ! outputs of check_block
    integer, allocatable, dimension(:,:,:) :: offset_position ! precomputed offsets
    real(double), allocatable, dimension(:) :: x_store, y_store, z_store, r_store ! outputs of check_block

    interface
       ! Interface to return a value val given arguments
       ! direction,species,l,acz,m,x,y,z. Implemented by
       ! evaluate_pao() and pao_elem_derivative_2().
       subroutine evaluate(direction,species,l,acz,m,x,y,z,val)
         use datatypes, only: double
         integer, intent(in) :: species,l,acz,m
         integer, intent(in) :: direction
         real(kind=double), intent(in) :: x,y,z
         real(kind=double), intent(out) :: val
       end subroutine evaluate
    end interface

    nblock = domain%groups_on_node
    npart = maxval(naba_atoms_of_blocks(atomf)%no_of_part)
    natom = maxval(naba_atoms_of_blocks(atomf)%no_atom_on_part)

    call start_timer(tmr_std_basis)
    call start_timer(tmr_std_allocation)
    allocate(ip_store(n_pts_in_block ))
    allocate(x_store( n_pts_in_block ))
    allocate(y_store( n_pts_in_block ))
    allocate(z_store( n_pts_in_block ))
    allocate(r_store( n_pts_in_block ))
    allocate(offset_position(natom, npart, nblock))
    call stop_timer(tmr_std_allocation)
    ! --  Start of subroutine  ---

    ! No need to compute these here, since they are derived from globals
    ! Store with rcellx, rcelly, rcellz?
    ! Reduces code duplication
    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz

    call my_barrier()

    gridfunctions(pao_fns)%griddata = zero
    next_offset_position = 0
    offset_position = huge(0)
    rcut = r_h + RD_ERR

    ! loop arround grid points in the domain, and for each
    ! point, get the PAO values

    ! Note: Using OpenMP in this loop requires some redesign because there is a loop
    !       carrier dependency in next_offset_position.
    blocks_loop: do iblock = 1, domain%groups_on_node ! primary set of blocks
       iatom = 0
       parts_loop: do ipart=1,naba_atoms_of_blocks(atomf)%no_of_part(iblock)
          naba_part_label = naba_atoms_of_blocks(atomf)%list_part(ipart,iblock)
          ind_part = DCS_parts%lab_cell(naba_part_label)
          atoms_loop: do ia=1,naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart,iblock)

             call get_species(iblock, naba_part_label, ind_part, iatom, icover, my_species)

             offset_position(ia, ipart, iblock) = next_offset_position
             next_offset_position = offset_position(ia, ipart, iblock) + &
                  npao_species(my_species) * n_pts_in_block
          end do atoms_loop
       end do parts_loop
    end do blocks_loop

    !$omp parallel do default(none) &
    !$omp             schedule(dynamic) &
    !$omp             shared(domain, naba_atoms_of_blocks, offset_position, pao_fns, atomf, &
    !$omp                    dcellx_block, dcelly_block, dcellz_block, dcs_parts, &
    !$omp                    rcut, n_pts_in_block, pao, gridfunctions, direction) &
    !$omp             private(ia, ipart, iblock, l1, acz, m1, count1, x, y, z, val, position, &
    !$omp                     npoint, r_store, ip_store, x_store, y_store, z_store, my_species, &
    !$omp                     xblock, yblock, zblock, iatom, xatom, yatom, zatom, naba_part_label, ind_part, icover)
    blocks_loop_omp: do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock = ( domain%idisp_primx(iblock) + domain%nx_origin - 1 ) * dcellx_block
       yblock = ( domain%idisp_primy(iblock) + domain%ny_origin - 1 ) * dcelly_block
       zblock = ( domain%idisp_primz(iblock) + domain%nz_origin - 1 ) * dcellz_block
       iatom = 0
       parts_loop_omp: do ipart=1,naba_atoms_of_blocks(atomf)%no_of_part(iblock)
          naba_part_label = naba_atoms_of_blocks(atomf)%list_part(ipart,iblock)
          ind_part = DCS_parts%lab_cell(naba_part_label)
          atoms_loop_omp: do ia=1,naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart,iblock)

             call get_species(iblock, naba_part_label, ind_part, iatom, icover, my_species)

             xatom = DCS_parts%xcover(icover)
             yatom = DCS_parts%ycover(icover)
             zatom = DCS_parts%zcover(icover)

             !calculates distances between the atom and integration grid points
             !in the block and stores which integration grids are neighbours.
             call check_block (xblock, yblock, zblock, xatom, yatom, zatom, rcut, &  ! in
                  npoint, ip_store, r_store, x_store, y_store, z_store, & !out
                  n_pts_in_block) ! in

             points_loop_omp: do ip=1,npoint
                position = offset_position(ia, ipart, iblock) + ip_store(ip)
                x = x_store(ip)
                y = y_store(ip)
                z = z_store(ip)
                ! For this point-atom offset, we accumulate the PAO on the grid
                l_loop: do l1 = 0,pao(my_species)%greatest_angmom
                   z_loop: do acz = 1,pao(my_species)%angmom(l1)%n_zeta_in_angmom
                      m_loop: do m1=-l1,l1
                         call evaluate(direction,my_species,l1,acz,m1,x,y,z,val)
                         gridfunctions(pao_fns)%griddata(position) = val
                         position = position + n_pts_in_block
                      end do m_loop
                   end do z_loop
                end do l_loop
             enddo points_loop_omp
          enddo atoms_loop_omp
       enddo parts_loop_omp
    enddo blocks_loop_omp
    !$omp end parallel do
    call my_barrier()
    call start_timer(tmr_std_allocation)
    deallocate(ip_store,x_store,y_store,z_store,r_store,offset_position)
    call stop_timer(tmr_std_allocation)
    call stop_timer(tmr_std_basis)
    return
  end subroutine PAO_or_gradPAO_to_grid
  !!***

! -----------------------------------------------------------
! Subroutine check_block
! -----------------------------------------------------------

!!****f* PAO_grid_transform_module/check_block *
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
!!
!!  SOURCE
!!
  subroutine check_block &
       (xblock,yblock,zblock,xatom,yatom,zatom,rcut, &
       npoint, ip_store, r_store, x_store, y_store, z_store,blocksize)

    use numbers
    use global_module, ONLY: rcellx,rcelly,rcellz
    use group_module,  ONLY: blocks
    use block_module,  ONLY: nx_in_block,ny_in_block,nz_in_block!, &
!         n_pts_in_block


    implicit none
    !Passed
    integer :: blocksize
    real(double), intent(in):: xblock, yblock, zblock
    real(double), intent(in):: xatom, yatom, zatom, rcut
    integer, intent(out) :: npoint, ip_store(blocksize)
    real(double),intent(out) :: r_store(blocksize)
    real(double),intent(out) :: x_store(blocksize)
    real(double),intent(out) :: y_store(blocksize)
    real(double),intent(out) :: z_store(blocksize)
    !Local
    real(double):: dcellx_block,dcelly_block,dcellz_block
    real(double):: dcellx_grid, dcelly_grid, dcellz_grid
    real(double):: dx, dy, dz
    integer :: ipoint, iz, iy, ix
    real(double) ::  r2, r_from_i, rx, ry, rz, x, y, z, rcut2


    rcut2 = rcut* rcut
    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz

    dcellx_grid=dcellx_block/nx_in_block
    dcelly_grid=dcelly_block/ny_in_block
    dcellz_grid=dcellz_block/nz_in_block

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
                if(npoint>blocksize) write(io_lun,*) 'ALERT: Too many points ! ',npoint,blocksize
                r_from_i = sqrt( r2 )

!                ! direction cosines needed for spher harms
!                if ( r_from_i > RD_ERR ) then
!                   x = rx / r_from_i
!                   y = ry / r_from_i
!                   z = rz / r_from_i
!                else
!                   x = zero ; y = zero ; z = zero
!                   !!   r_from_i = RD_ERR !!   04/07/2002 TM
!                end if
                x = rx
                y = ry
                z = rz
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

  !< Helper routine to get my_species and icover from block, particle and atom indices
  !>
  subroutine get_species(iblock, naba_part_label, ind_part, iatom, icover, my_species)

    use global_module, ONLY: id_glob, species_glob, atomf
    use set_blipgrid_module, ONLY : naba_atoms_of_blocks
    use cover_module, ONLY: DCS_parts
    use group_module, ONLY : blocks, parts

    integer, intent(in) :: iblock, naba_part_label, ind_part
    integer, intent(inout) :: iatom
    integer, intent(out) :: icover, my_species

    integer :: tmp_index

    iatom = iatom + 1

    tmp_index = naba_atoms_of_blocks(atomf)%list_atom(iatom,iblock) - 1

    icover= DCS_parts%icover_ibeg(naba_part_label) + tmp_index
    my_species = species_glob(id_glob(parts%icell_beg(ind_part) + tmp_index))

  end subroutine get_species

end module PAO_grid_transform_module
