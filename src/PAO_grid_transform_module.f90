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

!!****f* PAO_grid_transform_module/single_PAO_to_grid *
!!
!!  NAME
!!   single_PAO_to_grid
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
!!  SOURCE
!!
  subroutine single_PAO_to_grid(pao_fns)

    use datatypes
    use primary_module, ONLY: bundle
    use GenComms, ONLY: my_barrier, cq_abort, mtime
    use dimens, ONLY: r_h
    use GenComms, ONLY: inode, ionode
    use numbers
    use global_module, ONLY: rcellx,rcelly,rcellz,id_glob,ni_in_cell, iprint_basis, species_glob, atomf
    use species_module, ONLY: species, npao_species
    !  At present, these arrays are dummy arguments.
    use block_module, ONLY : nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
    use group_module, ONLY : blocks, parts
    use primary_module, ONLY: domain
    use cover_module, ONLY: DCS_parts
    use set_blipgrid_module, ONLY : naba_atoms_of_blocks
    use angular_coeff_routines, ONLY : evaluate_pao
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid
    use pao_format

    implicit none
    integer,intent(in) :: pao_fns

    !local
    integer :: iblock,ia,ipart,ip,l1,acz,m1 ! Loop index variables
    integer :: nblock, npart, natom ! Array dimensions
    integer :: jpart,ind_part,ii,icover,ig_atom,iatom ! Temporary indices for calculating species
    integer :: position,next_offset_position ! Temporary indices to gridfunctions%griddata
    integer :: my_species ! Temporary variables to reduce indirect accesses
    integer :: npoint ! outputs of check_block
    integer :: count1 ! incremented counter, maps from (l1, acz, m1) to linear index of gridfunctions%griddata
    real(double):: dcellx_block,dcelly_block,dcellz_block ! grid dimensions, should be moved
    real(double) :: x,y,z ! Temporary variables to reduce indirect accesses
    real(double) :: rcut ! Input to check_block
    real(double) :: val ! output, written into gridfunctions%griddata
    real(double) xblock,yblock,zblock ! inputs to check_block
    integer, allocatable, dimension(:) :: ip_store ! outputs of check_block
    integer, allocatable, dimension(:,:,:) :: offset_position, species_store ! precomputed variables
    real(double), allocatable, dimension(:) :: x_store, y_store, z_store, r_store ! outputs of check_block
    real(double), allocatable, dimension(:,:,:) :: xatom,yatom,zatom ! inputs to check_block, precomputed

    nblock = domain%groups_on_node
    npart = maxval(naba_atoms_of_blocks(atomf)%no_of_part)
    natom = maxval(naba_atoms_of_blocks(atomf)%no_atom_on_part)

    call start_timer(tmr_std_basis)
    call start_timer(tmr_std_allocation)
    allocate(xatom(natom, npart, nblock))
    allocate(yatom(natom, npart, nblock))
    allocate(zatom(natom, npart, nblock))
    allocate(ip_store(n_pts_in_block ))
    allocate(x_store( n_pts_in_block ))
    allocate(y_store( n_pts_in_block ))
    allocate(z_store( n_pts_in_block ))
    allocate(r_store( n_pts_in_block ))
    allocate(species_store(    natom, npart, nblock))
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
       part_in_block: if(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom=0
          parts_loop: do ipart=1,naba_atoms_of_blocks(atomf)%no_of_part(iblock)
             jpart=naba_atoms_of_blocks(atomf)%list_part(ipart,iblock)
             ind_part=DCS_parts%lab_cell(jpart)
             atoms_loop: do ia=1,naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atoms_of_blocks(atomf)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)

                xatom(ia, ipart, iblock) = DCS_parts%xcover(icover)
                yatom(ia, ipart, iblock) = DCS_parts%ycover(icover)
                zatom(ia, ipart, iblock) = DCS_parts%zcover(icover)
                species_store(ia, ipart, iblock) = species_glob(ig_atom)

                offset_position(ia, ipart, iblock) = next_offset_position
                next_offset_position = offset_position(ia, ipart, iblock) + &
                     npao_species(species_store(ia, ipart, iblock)) * n_pts_in_block
             end do atoms_loop
          end do parts_loop
       end if part_in_block
    end do blocks_loop

    !$omp parallel do default(none) &
    !$omp             schedule(dynamic) &
    !$omp             shared(domain, naba_atoms_of_blocks, offset_position, pao_fns, atomf, &
    !$omp                    dcellx_block, dcelly_block, dcellz_block, species_store, &
    !$omp                    xatom, yatom, zatom, rcut, n_pts_in_block, pao, gridfunctions) &
    !$omp             private(ia, ipart, iblock, l1, acz, m1, count1, x, y, z, val, position, &
    !$omp                     npoint, r_store, ip_store, x_store, y_store, z_store, my_species, &
    !$omp                     xblock, yblock, zblock)
    blocks_loop_omp: do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock = ( domain%idisp_primx(iblock) + domain%nx_origin - 1 ) * dcellx_block
       yblock = ( domain%idisp_primy(iblock) + domain%ny_origin - 1 ) * dcelly_block
       zblock = ( domain%idisp_primz(iblock) + domain%nz_origin - 1 ) * dcellz_block
       part_if_omp: if(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) then ! if there are naba atoms
          parts_loop_omp: do ipart=1,naba_atoms_of_blocks(atomf)%no_of_part(iblock)
             atoms_loop_omp: do ia=1,naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart,iblock)

                !calculates distances between the atom and integration grid points
                !in the block and stores which integration grids are neighbours.
                call check_block (xblock, yblock, zblock, &
                     xatom(ia, ipart, iblock), &
                     yatom(ia, ipart, iblock), &
                     zatom(ia, ipart, iblock), &
                     rcut, &  ! in
                     npoint, ip_store, r_store, x_store, y_store, z_store, & !out
                     n_pts_in_block) ! in

                my_species = species_store(ia, ipart, iblock)

                npoint_if_omp : if(npoint > 0) then
                   points_loop_omp: do ip=1,npoint
                      position = offset_position(ia, ipart, iblock) + ip_store(ip)
                      x = x_store(ip)
                      y = y_store(ip)
                      z = z_store(ip)
                      ! For this point-atom offset, we accumulate the PAO on the grid
                      count1 = 0
                      l_loop: do l1 = 0,pao(my_species)%greatest_angmom
                         z_loop: do acz = 1,pao(my_species)%angmom(l1)%n_zeta_in_angmom
                            m_loop: do m1=-l1,l1
                               call evaluate_pao(my_species,l1,acz,m1,x,y,z,val)
                               gridfunctions(pao_fns)%griddata(position + count1 * n_pts_in_block) = val
                               count1 = count1+1
                            end do m_loop
                         end do z_loop
                      end do l_loop
                   enddo points_loop_omp
                endif npoint_if_omp
             enddo atoms_loop_omp
          enddo parts_loop_omp
       endif part_if_omp
    enddo blocks_loop_omp
    !$omp end parallel do
    call my_barrier()
    call start_timer(tmr_std_allocation)
    ! Could just let these go out of scope at the end?
    deallocate(ip_store,x_store,y_store,z_store,r_store,species_store,offset_position)
    deallocate(xatom,yatom,zatom)
    call stop_timer(tmr_std_allocation)
    call stop_timer(tmr_std_basis)
    return
  end subroutine Single_PAO_to_grid
!!***

!!****f* PAO_grid_transform_module/single_PAO_to_grad *
!!
!!  NAME
!!   single_PAO_to_grad
!!  USAGE
!!
!!  PURPOSE
!!   Projects the gradient of PAO functions onto the grid (rather than support functions)
!!   Used for gradients of energy wrt atomic coordinates
!!
!!   This subroutine is based on sub:single_PAO_to_grid in PAO_grid_transform_module.f90.
!!   TODO: There is a lot of code duplication between single_PAO_to_grid and single_PAO_to_grad
!!
!!  INPUTS
!!
!!  USES
!!
!!  AUTHOR
!!   A. Nakata
!!  CREATION DATE
!!   2016/11/10
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine single_PAO_to_grad(direction, pao_fns)

    use datatypes
    use primary_module, ONLY: bundle
    use GenComms, ONLY: my_barrier, cq_abort, mtime
    use dimens, ONLY: r_h
    use GenComms, ONLY: inode, ionode
    use numbers
    use global_module, ONLY: rcellx,rcelly,rcellz,id_glob,ni_in_cell, iprint_basis, species_glob, atomf
    use species_module, ONLY: species, npao_species
    !  At present, these arrays are dummy arguments.
    use block_module, ONLY : nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
    use group_module, ONLY : blocks, parts
    use primary_module, ONLY: domain
    use cover_module, ONLY: DCS_parts
    use set_blipgrid_module, ONLY : naba_atoms_of_blocks
    use angular_coeff_routines, ONLY : pao_elem_derivative_2
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid
    use pao_format

    implicit none
    integer,intent(in) :: pao_fns
    integer,intent(in) :: direction

    !local
    real(double):: dcellx_block,dcelly_block,dcellz_block
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom
    real(double):: xatom,yatom,zatom
    real(double):: xblock,yblock,zblock
    integer :: the_species
    integer :: j,iblock,the_l,ipoint, igrid
    real(double) :: r_from_i
    real(double) :: rr,a,b,c,d,x,y,z,nl_potential
    integer :: no_of_ib_ia, offset_position
    integer :: position,iatom
    integer :: stat, nl, npoint, ip
    integer :: i,m, m1min, m1max,acz,m1,l1,count1
    integer     , allocatable :: ip_store(:)
    real(double), allocatable :: x_store(:)
    real(double), allocatable :: y_store(:)
    real(double), allocatable :: z_store(:)
    real(double), allocatable :: r_store(:)
    real(double) :: rcut
    real(double) :: val
    integer :: max_num_blocks, current_num_blocks
    real(double), dimension(:), allocatable :: temp_block_storage

    call start_timer(tmr_std_basis)
    call start_timer(tmr_std_allocation)
    allocate(ip_store(n_pts_in_block),x_store(n_pts_in_block),y_store(n_pts_in_block),z_store(n_pts_in_block), &
         r_store(n_pts_in_block))
    call stop_timer(tmr_std_allocation)
    ! --  Start of subroutine  ---

    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz

    call my_barrier()

    no_of_ib_ia = 0
    gridfunctions(pao_fns)%griddata = zero

    max_num_blocks = maxval(npao_species)*n_pts_in_block
    allocate(temp_block_storage(max_num_blocks))

    ! loop arround grid points in the domain, and for each
    ! point, get the d_PAO/d_R values

    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
       if(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom=0
          do ipart=1,naba_atoms_of_blocks(atomf)%no_of_part(iblock)
             jpart=naba_atoms_of_blocks(atomf)%list_part(ipart,iblock)
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atoms_of_blocks(atomf)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)

                xatom=DCS_parts%xcover(icover)
                yatom=DCS_parts%ycover(icover)
                zatom=DCS_parts%zcover(icover)
                the_species=species_glob(ig_atom)

                !calculates distances between the atom and integration grid points
                !in the block and stores which integration grids are neighbours.
                rcut = r_h + RD_ERR
                call check_block (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
                     npoint,ip_store,r_store,x_store,y_store,z_store,n_pts_in_block) !out

                if(npoint > 0) then
                   ! Temporary storage
                   current_num_blocks = npao_species(the_species)*n_pts_in_block
                   temp_block_storage = zero
                   offset_position = no_of_ib_ia
                   !$omp parallel do default(none) &
                   !$omp             schedule(dynamic) &
                   !$omp             reduction(+: temp_block_storage) &
                   !$omp             shared(n_pts_in_block, direction, pao, the_species, offset_position, &
                   !$omp                    npoint, no_of_ib_ia, ip_store, r_store, x_store, y_store, z_store) &
                   !$omp             private(r_from_i, ip, position, ipoint, &
                   !$omp                     x, y, z, l1, acz, m1, count1, val)
                   do ip=1,npoint
                      ipoint=ip_store(ip)
                      position= offset_position + ipoint

                      r_from_i = r_store(ip)
                      x = x_store(ip)
                      y = y_store(ip)
                      z = z_store(ip)
                      ! For this point-atom offset, we accumulate the PAO on the grid
                      count1 = 1

                      do l1 = 0,pao(the_species)%greatest_angmom
                         do acz = 1,pao(the_species)%angmom(l1)%n_zeta_in_angmom
                            do m1=-l1,l1
                               call pao_elem_derivative_2(direction,the_species,l1,acz,m1,x,y,z,val)
                               ! Each loop iteration accesses consequtive elements of gridfunctions%griddata
                               ! so this should be safe as a shared variable. Incrementing the index
                               ! inside the loop might not be great for vectorization, but with the function
                               ! call here, it probably won't vectorize anyways.
                               ! This loop nest really should be inside pao_elem_derivative_2, see issue #198
                               !gridfunctions(pao_fns)%griddata(position+(count1-1)*n_pts_in_block) = val
                               temp_block_storage(ipoint + (count1-1)*n_pts_in_block) = val
                               count1 = count1+1
                            end do ! m1
                         end do ! acz
                      end do ! l1
                   enddo ! ip=1,npoint
                   !$omp end parallel do
                   gridfunctions(pao_fns)%griddata(no_of_ib_ia+1:no_of_ib_ia+current_num_blocks) = &
                        temp_block_storage(1:current_num_blocks)
                endif! (npoint > 0) then
                no_of_ib_ia = no_of_ib_ia + npao_species(the_species)*n_pts_in_block
             enddo ! naba_atoms
          enddo ! naba_part
       endif !(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) !naba atoms?
    enddo ! iblock : primary set of blocks
    deallocate(temp_block_storage)
    call my_barrier()
    call start_timer(tmr_std_allocation)
    deallocate(ip_store,x_store,y_store,z_store,r_store)
    call stop_timer(tmr_std_allocation)
    call stop_timer(tmr_std_basis)
    return
  end subroutine Single_PAO_to_grad
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

end module PAO_grid_transform_module
