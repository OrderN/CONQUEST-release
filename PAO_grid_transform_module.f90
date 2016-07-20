! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!  SOURCE
!!
module PAO_grid_transform_module

  use global_module,          only: io_lun
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_basis, tmr_std_allocation

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
!!***

contains

!!****f* PAO_grid_transform_module/PAO_to_grid *
!!
!!  NAME 
!!   PAO_to_grid
!!  USAGE
!! 
!!  PURPOSE
!!   Takes support functions represented as PAOs and returns their values on the integration
!!   grid.
!!  INPUTS
!!   integer,intent(in) :: myid
!!   integer,intent(in) :: NSF, SUPPORT_SIZE
!!
!!   real(double),intent(out):: support(SUPPORT_SIZE*NSF)
!!  USES
!!   datatypes, numbers, primary_module, GenComms, blip_grid_transform_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   16:40, 2003/09/22
!!  MODIFICATION HISTORY
!!   11:54, 2003/10/03 dave
!!    Changed writes, removed stop
!!   2006/06/19 07:57 dave
!!    Made variable-NSF compatible
!!   2008/05/28 ast
!!    Added timers
!!  SOURCE
!!
  subroutine PAO_to_grid(myid, support)

    use datatypes
    use numbers, ONLY: zero
    use primary_module, ONLY: bundle
    use GenComms, ONLY: my_barrier, cq_abort, mtime
    use blip_grid_transform_module, ONLY: distribute_result
    use global_module, ONLY: iprint_basis, IPRINT_TIME_THRES2
    use functions_on_grid, ONLY: gridfunctions
    use species_module, ONLY: nsf_species
    use support_spec_format, ONLY: flag_paos_atoms_in_cell
    use timer_module

    implicit none

    integer,intent(in) :: myid

    integer :: support

    integer :: iprim, nsf_send
    integer :: ii
    real(double) :: now, nthen
    type(cq_timer) :: tmr_l_tmp1

    call start_timer(tmr_std_basis)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    gridfunctions(support)%griddata = zero

    nthen = mtime()
    if(bundle%mx_iprim < 1) call cq_abort("PAO_to_grid: no primary atoms ! ",myid,bundle%mx_iprim)
    if(flag_paos_atoms_in_cell) then
       call PAO_to_grid_global(myid,support)
    else
       do iprim = 1,bundle%mx_iprim
          call my_barrier()
          if( iprim <= bundle%n_prim ) then
             nsf_send = nsf_species(bundle%species(iprim))
             call do_PAO_transform(iprim,nsf_send)
             if(iprint_basis>1) then
                now = mtime()
                if(myid==0) write(io_lun,fmt='(2x,"Proc ",i5," Time for PAO_grid transform (in ms): ",f15.8)') myid,now-nthen
                nthen = now
             end if
          else
             nsf_send = 0
          endif
          call distribute_result(myid,iprim,nsf_send,support)
          if(iprint_basis>1) then
             now = mtime()
             if(myid==0) write(io_lun,fmt='(2x,"Proc ",i5," Time for PAO_grid distribute (in ms): ",f15.8)') myid,now - nthen
             nthen = now
          end if
       end do
    end if
    call stop_print_timer(tmr_l_tmp1,"PAO to grid transform",IPRINT_TIME_THRES2)
    call stop_timer(tmr_std_basis)
    return
  end subroutine PAO_to_grid
!!***

!!****f* PAO_grid_transform_module/do_PAO_transform *
!!
!!  NAME 
!!   do_PAO_transform
!!  USAGE
!! 
!!  PURPOSE
!!   Does the actual transform onto a local grid
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki/D.R.Bowler
!!  CREATION DATE
!!   21/6/2000/16:41, 2003/09/22
!!  MODIFICATION HISTORY
!!   16:42, 2003/09/22 dave
!!    Copied and modified do_blip_transform
!!   11:54, 2003/10/03 dave
!!    Changed writes, removed stop
!!   13:18, 03/10/2003 drb 
!!    Bug fix - added cq_abort to use GenComms
!!   16:01, 30/10/2003 drb 
!!    Changed calls, now works
!!   07:57, 2004/07/27 dave
!!    Changed to use new PAO/SF structure
!!   2008/05/28 ast
!!    Added timers
!!   2009/07/08 16:48 dave
!!    Added code for one-to-one PAO to SF assignment
!!   2016/07/14 16:30 nakata
!!    Renamed naba_blk_supp -> naba_blocks_of_atoms
!!  SOURCE
!!
  subroutine do_PAO_transform(iprim,this_nsf)

    use datatypes
    use numbers
    use global_module,   ONLY: rcellx,rcelly,rcellz, iprint_basis
    use group_module,    ONLY: blocks
    use primary_module,  ONLY: bundle
    use cover_module,    ONLY: BCS_blocks
    use set_blipgrid_module, ONLY: naba_blocks_of_atoms
    use comm_array_module,   ONLY: send_array
    use block_module,    ONLY: nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
    use support_spec_format, ONLY: supports_on_atom, flag_paos_atoms_in_cell, flag_one_to_one
    use GenComms, ONLY: myid, cq_abort
    use dimens, ONLY: r_h
    use angular_coeff_routines, ONLY: evaluate_pao

    implicit none

    integer,intent(in)     :: iprim
    integer,intent(in)     :: this_nsf

    ! Local variables
    real(double):: rec_sgs, dx, dy, dz
    integer     :: x, y, z
    integer     :: nsf1

    integer     :: nblkx,nblky,nblkz  ! nos. of integ. pts in a block 
    ! along x,y and z direction
    integer     :: naba_blk,nx_blk,ny_blk,nz_blk,igrid
    real(double) :: dcellx_block,dcelly_block,dcellz_block
    real(double) :: dcellx_grid, dcelly_grid, dcellz_grid
    real(double):: dsum(this_nsf), val

    integer     :: ix,iy,iz,acz, l1, m1, count1
    integer     :: ind_blk,ind_blip
    integer     :: nxmin_grid,nxmax_grid,nymin_grid,nymax_grid,nzmin_grid,nzmax_grid
    integer     :: ncover_yz, atom_species, this_atom, stat

    ! Species of atom
    atom_species = bundle%species(iprim)
    ! How are we storing PAO coefficients ?
    if(flag_one_to_one) then
       this_atom = atom_species
    else
       if(flag_paos_atoms_in_cell) then
          this_atom = bundle%ig_prim(iprim)
       else
          this_atom = iprim
       end if
    end if
    nblkx=nx_in_block
    nblky=ny_in_block
    nblkz=nz_in_block

    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nblkx
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/nblky
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nblkz

    ncover_yz=BCS_blocks%ncovery*BCS_blocks%ncoverz
    igrid=0

    call start_timer(tmr_std_allocation)
    allocate(send_array(naba_blocks_of_atoms%no_naba_blk(iprim)*n_pts_in_block*this_nsf),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating send_array in do_PAO_grid: ",&
         naba_blocks_of_atoms%no_naba_blk(iprim)*n_pts_in_block*this_nsf)
    call stop_timer(tmr_std_allocation)
    send_array(:) = zero

    ! Loop over neighbour blocks
    DO naba_blk=1,naba_blocks_of_atoms%no_naba_blk(iprim)! naba blks in NOPG order
       ! iprim : primary seq. no. of the atom
       ind_blk=naba_blocks_of_atoms%list_naba_blk(naba_blk,iprim) !CC in BCS_blocks
       ind_blk=ind_blk-1

       nx_blk=ind_blk/ncover_yz
       ny_blk=(ind_blk-nx_blk*ncover_yz)/BCS_blocks%ncoverz
       nz_blk=ind_blk-nx_blk*ncover_yz-ny_blk*BCS_blocks%ncoverz

       nx_blk=nx_blk-BCS_blocks%nspanlx+BCS_blocks%nx_origin
       ny_blk=ny_blk-BCS_blocks%nspanly+BCS_blocks%ny_origin
       nz_blk=nz_blk-BCS_blocks%nspanlz+BCS_blocks%nz_origin

       do iz=1,nblkz    ! z-direction in a block
          z=nblkz*(nz_blk-1)+iz-1
          do iy=1,nblky    ! y-direction in a block
             y=nblky*(ny_blk-1)+iy-1
             do ix=1,nblkx    ! x-direction in a block
                x=nblkx*(nx_blk-1)+ix-1
                ! I want to find the position of the grid point relative to the atom in its primary set.
                ! The grid-point position _would_ be x + xadd, with xadd = origin + cover block x - FSC block x
                ! cover block x is nx_blk BEFORE the manipulation above
                ! 
                ! Actually, isn't nx_blk AFTER the manipulation just what I want ? We subtract the left-span 
                ! and add the origin
                ! 
                ! OK - let's assume that x, y and z are EXACTLY what we want for the grid point.
                igrid=igrid+1  ! seq. no. of integration grids
                !if(this_nsf*igrid > mx_send_real) call cq_abort("do_PAO_transform: overflow ",igrid,mx_send_real)
                ! generate the offset vector                
                dx = x*dcellx_grid - bundle%xprim(iprim)
                dy = y*dcelly_grid - bundle%yprim(iprim)
                dz = z*dcellz_grid - bundle%zprim(iprim)
                if(dx*dx+dy*dy+dz*dz<r_h*r_h) then
                   ! loop over PAOs
                   dsum = zero
                   count1 = 1
                   do l1 = 0,supports_on_atom(this_atom)%lmax
                      do acz = 1,supports_on_atom(this_atom)%naczs(l1)
                         do m1=-l1,l1
                            call evaluate_pao(atom_species,l1,1,m1,dx,dy,dz,val)
                            if(flag_one_to_one) then
                                  dsum(count1) = dsum(count1) + val
                            else
                               do nsf1 = 1,supports_on_atom(this_atom)%nsuppfuncs
                                  dsum(nsf1) = dsum(nsf1) + &
                                       supports_on_atom(this_atom)%supp_func(nsf1)%coefficients(count1)* val
                               end do ! nsf
                            end if
                            count1 = count1+1
                         end do ! m1
                      end do ! acz
                   end do ! l1
                   ! store in send_array
                   do nsf1 = 1,this_nsf
                      send_array(this_nsf*igrid-this_nsf+nsf1) = dsum(nsf1)
                   enddo ! nsf1=this_nsf
                end if ! if(dx*dx...<r_h*r_h)
             enddo    ! x-direction in a block
          enddo    ! y-direction in a block
       enddo    ! z-direction in a block
    end do    ! Loop over naba blocks for this primary atom
    return
  end subroutine do_PAO_transform
!!***

!!****f* PAO_grid_transform_module/PAO_to_grad *
!!
!!  NAME 
!!   PAO_to_grad
!!  USAGE
!! 
!!  PURPOSE
!!   Takes support functions represented as PAOs and returns their differentials in a 
!!   specified direction on the integration grid.
!!  INPUTS
!!   integer,intent(in) :: myid
!!   integer,intent(in) :: direction
!!   integer,intent(in) :: NSF, SUPPORT_SIZE
!!
!!   real(double),intent(out):: support(SUPPORT_SIZE*NSF)
!!  USES
!!   datatypes, numbers, primary_module, GenComms, blip_grid_transform_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   08:01, 2003/10/01
!!  MODIFICATION HISTORY
!!   11:53, 2003/10/03 dave
!!    Changed write statements, removed stop
!!   02:30, 2003/11/13 dave
!!    Added direction to call to do_PAO_grad
!!   2008/05/28 ast
!!    Added timers
!!  SOURCE
!!
  subroutine PAO_to_grad(myid, direction, support)

    use datatypes
    use numbers, ONLY: zero
    use global_module, ONLY: iprint_basis
    use primary_module, ONLY: bundle
    use GenComms, ONLY: my_barrier, cq_abort, mtime
    use blip_grid_transform_module, ONLY: distribute_result
    use functions_on_grid, ONLY: gridfunctions
    use species_module, ONLY: nsf_species
    use support_spec_format, ONLY: flag_paos_atoms_in_cell

    implicit none

    integer,intent(in) :: myid, direction
    integer :: support

    integer :: iprim, nsf_send
    integer :: ii
    real(double) :: now, nthen

    call start_timer(tmr_std_basis)
    gridfunctions(support)%griddata = zero

    nthen = mtime()
    if(bundle%mx_iprim < 1) call cq_abort("PAO_to_grad: no primary atoms ! ",myid,bundle%mx_iprim)
    if(flag_paos_atoms_in_cell) then
       call PAO_to_grad_global(myid,direction,support)
    else
       do iprim = 1,bundle%mx_iprim
          call my_barrier()
          if( iprim <= bundle%n_prim ) then
             nsf_send = nsf_species(bundle%species(iprim))
             call do_PAO_grad_transform(direction,iprim,nsf_send)
		     if(iprint_basis>1) then
                now = mtime()
                if(myid==0) write(io_lun,fmt='(2x,"Proc ",i5," Time for PAO_grad transform (in ms): ",f15.8)') myid,now - nthen
                nthen = now
             end if
          else
             nsf_send = 0
          endif
          call distribute_result(myid,iprim,nsf_send,support)
          if(iprint_basis>1) then
             now = mtime()
             if(myid==0) write(io_lun,fmt='(2x,"Proc ",i5," Time for PAO_grad distribute (in ms): ",f15.8)') myid,now - nthen
             nthen = now
          end if
       end do
    end if
    call stop_timer(tmr_std_basis)
    return
  end subroutine PAO_to_grad
!!***

!!****f* PAO_grid_transform_module/do_PAO_grad_transform *
!!
!!  NAME 
!!   do_PAO_grad_transform
!!  USAGE
!! 
!!  PURPOSE
!!   Does the actual transform of the derivative of the PAO onto a local grid
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   11:53, 2003/10/03 
!!  MODIFICATION HISTORY
!!   16:42, 2003/09/22 dave
!!    Copied and modified do_blip_transform
!!   11:53, 2003/10/03 dave
!!    Changed write statements, removed stop
!!   16:08, 30/10/2003 drb 
!!    Brought up-to-date with do_PAO_grid_transform: still needs DIRECTION !
!!   02:30, 2003/11/13 dave
!!    Added direction
!!   07:57, 2004/07/27 dave
!!    Changed to use new PAO/SF structure
!!   2008/05/28 ast
!!    Added timers
!!   2009/07/08 16:48 dave
!!    Added code for one-to-one PAO to SF assignment
!!   2016/07/14 16:30 nakata
!!    Renamed naba_blk_supp -> naba_blocks_of_atoms
!!  SOURCE
!!
  subroutine do_PAO_grad_transform(direction,iprim,this_nsf)

    use datatypes
    use numbers
    use global_module,   ONLY: rcellx,rcelly,rcellz,iprint_basis
    use group_module,    ONLY: blocks
    use primary_module,  ONLY: bundle
    use cover_module,    ONLY: BCS_blocks
    use set_blipgrid_module, ONLY: naba_blocks_of_atoms
    use comm_array_module,   ONLY: send_array
    use block_module,    ONLY: nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
    use support_spec_format, ONLY: supports_on_atom, flag_paos_atoms_in_cell, flag_one_to_one
    use GenComms, ONLY: myid, cq_abort
    use dimens, ONLY: r_h
    use angular_coeff_routines, ONLY: pao_elem_derivative_2, evaluate_pao

    implicit none

    integer,intent(in)     :: iprim, direction
    integer,intent(in)     :: this_nsf

    ! Local variables
    real(double):: rec_sgs, dx, dy, dz
    integer     :: x, y, z
    integer     :: nsf1

    integer     :: nblkx,nblky,nblkz  ! nos. of integ. pts in a block 
    ! along x,y and z direction
    integer     :: naba_blk,nx_blk,ny_blk,nz_blk,igrid
    real(double) :: dcellx_block,dcelly_block,dcellz_block
    real(double) :: dcellx_grid, dcelly_grid, dcellz_grid
    real(double):: dsum(this_nsf),val,val2,pao_val,pao_val2

    integer     :: ix,iy,iz,acz,l1,m1,count1
    integer     :: ind_blk,ind_blip
    integer     :: nxmin_grid,nxmax_grid,nymin_grid,nymax_grid,nzmin_grid,nzmax_grid
    integer     :: ncover_yz, atom_species, this_atom, stat

    ! Species of atom
    atom_species = bundle%species(iprim)
    ! How are we storing PAO coefficients ?
    if(flag_one_to_one) then
       this_atom = atom_species
    else
       if(flag_paos_atoms_in_cell) then
          this_atom = bundle%ig_prim(iprim)
       else
          this_atom = iprim
       end if
    end if
    nblkx=nx_in_block
    nblky=ny_in_block
    nblkz=nz_in_block

    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nblkx
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/nblky
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nblkz

    ncover_yz=BCS_blocks%ncovery*BCS_blocks%ncoverz
    igrid=0

    call start_timer(tmr_std_allocation)
    allocate(send_array(naba_blocks_of_atoms%no_naba_blk(iprim)*n_pts_in_block*this_nsf),&
             STAT=stat)
    if(stat/=0) call cq_abort("Error allocating send_array in do_PAO_grad: ",&
         naba_blocks_of_atoms%no_naba_blk(iprim)*n_pts_in_block*this_nsf)
    call stop_timer(tmr_std_allocation)
    send_array(:) = zero

    ! Loop over neighbour blocks
    DO naba_blk=1,naba_blocks_of_atoms%no_naba_blk(iprim)! naba blks in NOPG order
       ! iprim : primary seq. no. of the atom
       ind_blk=naba_blocks_of_atoms%list_naba_blk(naba_blk,iprim) !CC in BCS_blocks
       ind_blk=ind_blk-1

       nx_blk=ind_blk/ncover_yz
       ny_blk=(ind_blk-nx_blk*ncover_yz)/BCS_blocks%ncoverz
       nz_blk=ind_blk-nx_blk*ncover_yz-ny_blk*BCS_blocks%ncoverz

       nx_blk=nx_blk-BCS_blocks%nspanlx+BCS_blocks%nx_origin
       ny_blk=ny_blk-BCS_blocks%nspanly+BCS_blocks%ny_origin
       nz_blk=nz_blk-BCS_blocks%nspanlz+BCS_blocks%nz_origin

       do iz=1,nblkz    ! z-direction in a block
          z=nblkz*(nz_blk-1)+iz-1
          do iy=1,nblky    ! y-direction in a block
             y=nblky*(ny_blk-1)+iy-1
             do ix=1,nblkx    ! x-direction in a block
                x=nblkx*(nx_blk-1)+ix-1
                ! I want to find the position of the grid point relative to the atom in its primary set.
                ! The grid-point position _would_ be x + xadd, with xadd = origin + cover block x - FSC block x
                ! cover block x is nx_blk BEFORE the manipulation above
                ! 
                ! Actually, isn't nx_blk AFTER the manipulation just what I want ? We subtract the left-span and add the origin
                ! 
                ! OK - let's assume that x, y and z are EXACTLY what we want for the grid point.
                igrid=igrid+1  ! seq. no. of integration grids
                if(this_nsf*igrid > naba_blocks_of_atoms%no_naba_blk(iprim)*n_pts_in_block*this_nsf) &
                     call cq_abort("do_PAO_grad_transform: overflow ",this_nsf*igrid, &
                     naba_blocks_of_atoms%no_naba_blk(iprim)*n_pts_in_block*this_nsf)
                ! generate the offset vector                
                dx = x*dcellx_grid - bundle%xprim(iprim)
                dy = y*dcelly_grid - bundle%yprim(iprim)
                dz = z*dcellz_grid - bundle%zprim(iprim)
                if(dx*dx+dy*dy+dz*dz<r_h*r_h) then
                   ! loop over PAOs
                   dsum = zero
                   count1 = 1
                   do l1 = 0,supports_on_atom(this_atom)%lmax
                      do acz = 1,supports_on_atom(this_atom)%naczs(l1)
                         do m1=-l1,l1
                            call pao_elem_derivative_2(direction,atom_species,l1,acz,m1,dx,dy,dz,val)
                            if(flag_one_to_one) then
                               dsum(count1) = dsum(count1) + val
                            else
                               do nsf1 = 1,supports_on_atom(this_atom)%nsuppfuncs
                                  dsum(nsf1) = dsum(nsf1) + &
                                       supports_on_atom(this_atom)%supp_func(nsf1)%coefficients(count1)* val
                               end do ! nsf
                            end if
                            count1 = count1+1
                         end do ! m1
                      end do ! acz
                   end do ! l1
                   ! store in send_array
                   do nsf1 = 1,this_nsf
                      send_array(this_nsf*igrid-this_nsf+nsf1) = dsum(nsf1)
                   enddo
                end if
             enddo    ! x-direction in a block
          enddo    ! y-direction in a block
       enddo    ! z-direction in a block
    end do    ! Loop over naba blocks for this primary atom
    return
  end subroutine do_PAO_grad_transform
!!***

!!****f* PAO_grid_transform_module/PAO_to_grid_global *
!!
!!  NAME 
!!   PAO_to_grid_global
!!  USAGE
!! 
!!  PURPOSE
!!   Transforms PAOs to grid assuming global PAO coefficient storage: no communication involved
!!  INPUTS
!!   myid: processor number
!!   support: index for grid function
!!  USES
!! 
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2004 sometime probably
!!  MODIFICATION HISTORY
!!   2006/06/20 08:24 dave
!!    Various changes for variable NSF
!!   2008/05/28 ast
!!    Added timers
!!   2009/07/08 16:48 dave
!!    Added code for one-to-one PAO to SF assignment
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks
!!  SOURCE
!!
  subroutine PAO_to_grid_global(myid, support)

    use datatypes
    use primary_module, ONLY: bundle
    use GenComms, ONLY: my_barrier, cq_abort, mtime
    use dimens, ONLY: r_h
    use GenComms, ONLY: inode, ionode, cq_abort
    use numbers
    use global_module, ONLY: rcellx,rcelly,rcellz,id_glob,ni_in_cell, iprint_basis, species_glob, sf
    use species_module, ONLY: species, nsf_species
    !  At present, these arrays are dummy arguments.
    use block_module, ONLY : nx_in_block,ny_in_block,nz_in_block, &
         n_pts_in_block
    use group_module, ONLY : blocks, parts
    use primary_module, ONLY: domain
    use cover_module, ONLY: DCS_parts
    use set_blipgrid_module, ONLY : naba_atoms_of_blocks
    use angular_coeff_routines, ONLY : evaluate_pao
    use support_spec_format, ONLY: supports_on_atom, flag_paos_atoms_in_cell, flag_one_to_one
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid

    implicit none 
    integer,intent(in) :: myid
    integer :: support


    !local
    real(double):: dcellx_block,dcelly_block,dcellz_block
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom, this_atom
    real(double):: xatom,yatom,zatom,alpha,step
    real(double):: xblock,yblock,zblock
    integer :: the_species
    integer :: j,iblock,the_l,ipoint, igrid
    real(double) :: r_from_i
    real(double) :: rr,a,b,c,d,x,y,z,nl_potential
    integer :: no_of_ib_ia, offset_position
    integer :: position,iatom
    integer :: stat, nl, npoint, ip, this_nsf
    integer :: i,m, m1min, m1max,acz,m1,l1,nsf1,atom_species, count1
    integer     , allocatable :: ip_store(:)
    real(double), allocatable :: x_store(:)
    real(double), allocatable :: y_store(:)
    real(double), allocatable :: z_store(:)
    real(double), allocatable :: r_store(:)
    real(double) :: coulomb_energy
    real(double) :: rcut
    real(double) :: r1, r2, r3, r4, core_charge, gauss_charge
    real(double) :: val, theta, phi, r_tmp
    real(double), dimension(:), allocatable :: sfsum

    if(.NOT.flag_paos_atoms_in_cell) call cq_abort("Calling global PAO-to-grid with local PAO storage")
    call my_barrier()

    call start_timer(tmr_std_allocation)
    allocate(ip_store(n_pts_in_block),x_store(n_pts_in_block),y_store(n_pts_in_block),z_store(n_pts_in_block), &
         r_store(n_pts_in_block))
    call stop_timer(tmr_std_allocation)
    ! --  Start of subroutine  ---

    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz

    no_of_ib_ia = 0
    gridfunctions(support)%griddata = zero

    ! loop arround grid points in the domain, and for each
    ! point, get the core charge (local charge) density from the atoms 
    ! whose distances from the grid point are within the cutoff. 

    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
       if(naba_atoms_of_blocks(sf)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom=0
          do ipart=1,naba_atoms_of_blocks(sf)%no_of_part(iblock)
             jpart=naba_atoms_of_blocks(sf)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('PAO_to_grid_global: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atoms_of_blocks(sf)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atoms_of_blocks(sf)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)

                if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                   call cq_abort('PAO_to_grid_global: globID ERROR ', &
                        ii,parts%icell_beg(ind_part))
                endif
                if(icover > DCS_parts%mx_mcover) then
                   call cq_abort('PAO_to_grid_global: icover ERROR ', &
                        icover,DCS_parts%mx_mcover)
                endif

                !-- I (TM) am not sure where the next line should be.
                !   (no_of_ib_ia = no_of_ib_ia +1)
                !  This determines how to store 
                !   pseudofunctions(ipoints,ncf,naba_atm,iblock).
                !  Now, I assume we should consider all naba_atm whose 
                ! distance from the block is within the maximum of core_radius.
                ! This is needed to keep the consistency with <set_bucket>.
                ! However, we can change this strategy by changing naba_atoms_of_blocks(sf).
                !no_of_ib_ia = no_of_ib_ia +1

                xatom=DCS_parts%xcover(icover)
                yatom=DCS_parts%ycover(icover)
                zatom=DCS_parts%zcover(icover)

                the_species=species_glob(ig_atom)
                this_nsf = nsf_species(the_species)
                call start_timer(tmr_std_allocation)
                allocate(sfsum(this_nsf))
                call stop_timer(tmr_std_allocation)
                
                !calculates distances between the atom and integration grid points
                !in the block and stores which integration grids are neighbours.
                rcut = r_h + RD_ERR
                call check_block (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
                     npoint,ip_store,r_store,x_store,y_store,z_store,n_pts_in_block) !out
                !r_from_i = sqrt((xatom-xblock)**2+(yatom-yblock)**2+ &
                !     (zatom-zblock)**2 )
                if(flag_one_to_one) then
                   this_atom = the_species
                else
                   this_atom = ig_atom
                end if
                if(npoint > 0) then
                   !offset_position = (no_of_ib_ia-1) * this_nsf * n_pts_in_block
                   offset_position = no_of_ib_ia
                   do ip=1,npoint
                      ipoint=ip_store(ip)
                      position= offset_position + ipoint
                      if(position > gridfunctions(support)%size) &
                           call cq_abort('PAO_to_grid_global: position error 1 ', &
                           position, gridfunctions(support)%size)
                      r_from_i = r_store(ip)
                      x = x_store(ip)
                      y = y_store(ip)
                      z = z_store(ip)
                      sfsum = zero
                      ! For this point-atom offset, we loop over ALL PAOs, and accumulate each PAO multiplied by
                      ! the appropriate coefficient in a local store before putting it in place later
                      count1 = 1
                      do l1 = 0,supports_on_atom(this_atom)%lmax
                         do acz = 1,supports_on_atom(this_atom)%naczs(l1)
                            do m1=-l1,l1
                               call evaluate_pao(the_species,l1,acz,m1,x,y,z,val)
                               if(flag_one_to_one) then
                                  sfsum(count1) = sfsum(count1) + val
                               else
                                  do nsf1 = 1,supports_on_atom(this_atom)%nsuppfuncs
                                     sfsum(nsf1) = sfsum(nsf1) + &
                                          supports_on_atom(this_atom)%supp_func(nsf1)%coefficients(count1)* val
                                  end do ! nsf1
                               end if
                               count1 = count1+1
                            end do ! m1
                         end do ! acz
                      end do ! l1
                      ! Now store data
                      if(position+(this_nsf-1)*n_pts_in_block > gridfunctions(support)%size) &
                           call cq_abort('PAO_to_grid_global: position error 2 ', &
                           position+(this_nsf-1)*n_pts_in_block, gridfunctions(support)%size)
                      do nsf1 = 1,this_nsf
                         gridfunctions(support)%griddata(position+(nsf1-1)*n_pts_in_block) = sfsum(nsf1)
                      end do ! nsf
                   enddo ! ip=1,npoint
                endif! (npoint > 0) then
                call start_timer(tmr_std_allocation)
                deallocate(sfsum)
                call stop_timer(tmr_std_allocation)
                no_of_ib_ia = no_of_ib_ia + this_nsf*n_pts_in_block
             enddo ! naba_atoms
          enddo ! naba_part
       endif !(naba_atoms_of_blocks(sf)%no_of_part(iblock) > 0) !naba atoms?
    enddo ! iblock : primary set of blocks
    call my_barrier()
    call start_timer(tmr_std_allocation)
    deallocate(ip_store,x_store,y_store,z_store,r_store)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine PAO_to_grid_global
!!***

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
!!  SOURCE
!!
  subroutine single_PAO_to_grid(support)

    use datatypes
    use primary_module, ONLY: bundle
    use GenComms, ONLY: my_barrier, cq_abort, mtime
    use dimens, ONLY: r_h
    use GenComms, ONLY: inode, ionode
    use numbers
    use global_module, ONLY: rcellx,rcelly,rcellz,id_glob,ni_in_cell, iprint_basis, species_glob, sf, paof
    use species_module, ONLY: species, npao_species
    !  At present, these arrays are dummy arguments.
    use block_module, ONLY : nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
    use group_module, ONLY : blocks, parts
    use primary_module, ONLY: domain
    use cover_module, ONLY: DCS_parts
    use set_blipgrid_module, ONLY : naba_atoms_of_blocks
    use angular_coeff_routines, ONLY : evaluate_pao
    use support_spec_format, ONLY: supports_on_atom, flag_one_to_one
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid

    implicit none 
    integer,intent(in) :: support

    !local
    real(double):: dcellx_block,dcelly_block,dcellz_block
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom, this_atom
    real(double):: xatom,yatom,zatom,alpha,step
    real(double):: xblock,yblock,zblock
    integer :: the_species
    integer :: j,iblock,the_l,ipoint, igrid
    real(double) :: r_from_i
    real(double) :: rr,a,b,c,d,x,y,z,nl_potential
    integer :: no_of_ib_ia, offset_position
    integer :: position,iatom
    integer :: stat, nl, npoint, ip, this_nsf
    integer :: i,m, m1min, m1max,acz,m1,l1,npao1,atom_species, count1
    integer     , allocatable :: ip_store(:)
    real(double), allocatable :: x_store(:)
    real(double), allocatable :: y_store(:)
    real(double), allocatable :: z_store(:)
    real(double), allocatable :: r_store(:)
    real(double) :: coulomb_energy
    real(double) :: rcut
    real(double) :: r1, r2, r3, r4, core_charge, gauss_charge
    real(double) :: val, theta, phi, r_tmp

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
    gridfunctions(support)%griddata = zero

    ! loop arround grid points in the domain, and for each
    ! point, get the PAO values

    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
       if(naba_atoms_of_blocks(paof)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom=0
          do ipart=1,naba_atoms_of_blocks(paof)%no_of_part(iblock)
             jpart=naba_atoms_of_blocks(paof)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('single_PAO_to_grid: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atoms_of_blocks(paof)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atoms_of_blocks(paof)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)

                if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                   call cq_abort('single_PAO_to_grid: globID ERROR ', &
                        ii,parts%icell_beg(ind_part))
                endif
                if(icover > DCS_parts%mx_mcover) then
                   call cq_abort('single_PAO_to_grid: icover ERROR ', &
                        icover,DCS_parts%mx_mcover)
                endif

                xatom=DCS_parts%xcover(icover)
                yatom=DCS_parts%ycover(icover)
                zatom=DCS_parts%zcover(icover)
                the_species=species_glob(ig_atom)

                !calculates distances between the atom and integration grid points
                !in the block and stores which integration grids are neighbours.
                rcut = r_h + RD_ERR
                call check_block (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
                     npoint,ip_store,r_store,x_store,y_store,z_store,n_pts_in_block) !out
                r_from_i = sqrt((xatom-xblock)**2+(yatom-yblock)**2+ &
                     (zatom-zblock)**2 )

                if(flag_one_to_one) then
                   this_atom = the_species
                else
                   this_atom = ig_atom
                end if
                if(npoint > 0) then
                   !offset_position = (no_of_ib_ia-1) * npao * n_pts_in_block
                   offset_position = no_of_ib_ia
                   do ip=1,npoint
                      ipoint=ip_store(ip)
                      position= offset_position + ipoint
                      if(position > gridfunctions(support)%size) call cq_abort &
                           ('single_pao_to_grid: position error ', position, gridfunctions(support)%size)

                      r_from_i = r_store(ip)
                      x = x_store(ip)
                      y = y_store(ip)
                      z = z_store(ip)
                      ! For this point-atom offset, we accumulate the PAO on the grid
                      count1 = 1
                      do l1 = 0,supports_on_atom(this_atom)%lmax
                         do acz = 1,supports_on_atom(this_atom)%naczs(l1)
                            do m1=-l1,l1
                               call evaluate_pao(the_species,l1,acz,m1,x,y,z,val)
                               if(position+(count1-1)*n_pts_in_block > gridfunctions(support)%size) &
                                    call cq_abort('single_pao_to_grid: position error ', &
                                    position, gridfunctions(support)%size)
                               gridfunctions(support)%griddata(position+(count1-1)*n_pts_in_block) = val
                               count1 = count1+1
                            end do ! m1
                         end do ! acz
                      end do ! l1
                   enddo ! ip=1,npoint
                endif! (npoint > 0) then
                no_of_ib_ia = no_of_ib_ia + npao_species(the_species)*n_pts_in_block
             enddo ! naba_atoms
          enddo ! naba_part
       endif !(naba_atoms_of_blocks(paof)%no_of_part(iblock) > 0) !naba atoms?
    enddo ! iblock : primary set of blocks
    call my_barrier()
    call start_timer(tmr_std_allocation)
    deallocate(ip_store,x_store,y_store,z_store,r_store)
    call stop_timer(tmr_std_allocation)
    call stop_timer(tmr_std_basis)
    return
  end subroutine Single_PAO_to_grid
!!***

!!****f* PAO_grid_transform_module/PAO_to_grad_global *
!!
!!  NAME 
!!   PAO_to_grad_Global
!!  USAGE
!! 
!!  PURPOSE
!!   Transforms gradient of support functions to grid using PAO basis functions
!!   This version assumes all PAO coeffients on all atoms are stored on all processors, 
!!   so no communication is required
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2004 sometime
!!  MODIFICATION HISTORY
!!   2006/06/20 08:30 dave
!!    Variable NSF preparation, tidying, name change
!!   2008/05/28 ast
!!    Added timers
!!   2009/07/08 16:48 dave
!!    Added code for one-to-one PAO to SF assignment
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks
!!  SOURCE
!!
  subroutine PAO_to_grad_global(myid, direction, support)

    use datatypes
    use primary_module, ONLY: bundle
    use GenComms, ONLY: my_barrier, cq_abort, mtime
    use dimens, ONLY: r_h
    use GenComms, ONLY: inode, ionode
    use numbers
    use global_module, ONLY: rcellx,rcelly,rcellz,id_glob,ni_in_cell, iprint_basis, species_glob, sf
    use species_module, ONLY: species, nsf_species
    !  At present, these arrays are dummy arguments.
    use block_module, ONLY : nx_in_block,ny_in_block,nz_in_block, &
         n_pts_in_block
    use group_module, ONLY : blocks, parts
    use primary_module, ONLY: domain
    use cover_module, ONLY: DCS_parts
    use set_blipgrid_module, ONLY : naba_atoms_of_blocks
    use angular_coeff_routines, ONLY : pao_elem_derivative_2
    use support_spec_format, ONLY: supports_on_atom, flag_one_to_one
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid

    implicit none 
    integer,intent(in) :: myid, direction

    integer :: support


    !local
    real(double):: dcellx_block,dcelly_block,dcellz_block
    integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom, this_atom
    real(double):: xatom,yatom,zatom,alpha,step
    real(double):: xblock,yblock,zblock
    integer :: the_species
    integer :: j,iblock,the_l,ipoint, igrid
    real(double) :: r_from_i
    real(double) :: rr,a,b,c,d,x,y,z,nl_potential
    integer :: no_of_ib_ia, offset_position
    integer :: position,iatom
    integer :: stat, nl, npoint, ip, this_nsf
    integer :: i,m, m1min, m1max,acz,m1,l1,nsf1,atom_species, count1
    integer     , allocatable :: ip_store(:)
    real(double), allocatable :: x_store(:)
    real(double), allocatable :: y_store(:)
    real(double), allocatable :: z_store(:)
    real(double), allocatable :: r_store(:)
    real(double) :: coulomb_energy
    real(double) :: rcut
    real(double) :: r1, r2, r3, r4, core_charge, gauss_charge
    real(double) :: val, theta, phi, r_tmp
    real(double), dimension(:), allocatable :: sfsum

    call start_timer(tmr_std_allocation)
    allocate(ip_store(n_pts_in_block),x_store(n_pts_in_block),y_store(n_pts_in_block),z_store(n_pts_in_block), &
         r_store(n_pts_in_block), STAT=stat)
    if(stat /= 0) call cq_abort(' Error allocating store in PAO_to_grad_global: ',n_pts_in_block)
    call stop_timer(tmr_std_allocation)
    ! --  Start of subroutine  ---

    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz

    call my_barrier()

    no_of_ib_ia = 0
    gridfunctions(support)%griddata = zero

    ! loop arround grid points in the domain, and for each
    ! point, get the core charge (local charge) density from the atoms 
    ! whose distances from the grid point are within the cutoff. 

    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
       if(naba_atoms_of_blocks(sf)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom=0
          do ipart=1,naba_atoms_of_blocks(sf)%no_of_part(iblock)
             jpart=naba_atoms_of_blocks(sf)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('PAO_to_grad_global: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atoms_of_blocks(sf)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atoms_of_blocks(sf)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)

                if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                   call cq_abort('PAO_to_grad_global: globID ERROR ', &
                        ii,parts%icell_beg(ind_part))
                endif
                if(icover > DCS_parts%mx_mcover) then
                   call cq_abort('PAO_to_grad_global: icover ERROR ', &
                        icover,DCS_parts%mx_mcover)
                endif

                xatom=DCS_parts%xcover(icover)
                yatom=DCS_parts%ycover(icover)
                zatom=DCS_parts%zcover(icover)

                the_species=species_glob(ig_atom)
                this_nsf = nsf_species(the_species)
                call start_timer(tmr_std_allocation)
                allocate(sfsum(this_nsf))
                call stop_timer(tmr_std_allocation)

                !calculates distances between the atom and integration grid points
                !in the block and stores which integration grids are neighbours.
                rcut = r_h + RD_ERR
                call check_block (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
                     npoint,ip_store,r_store,x_store,y_store,z_store,n_pts_in_block) !out
                r_from_i = sqrt((xatom-xblock)**2+(yatom-yblock)**2+ &
                     (zatom-zblock)**2 )
                if(flag_one_to_one) then
                   this_atom = the_species
                else
                   this_atom = ig_atom
                end if

                if(npoint > 0) then
                   !offset_position = (no_of_ib_ia-1) * nsf * n_pts_in_block
                   offset_position = no_of_ib_ia
                   do ip=1,npoint
                      ipoint=ip_store(ip)
                      position= offset_position + ipoint
                      if(position > gridfunctions(support)%size) call cq_abort &
                           ('PAO_to_grad_global: position error ', position, gridfunctions(support)%size)

                      r_from_i = r_store(ip)
                      x = x_store(ip)
                      y = y_store(ip)
                      z = z_store(ip)
                      sfsum = zero
                      ! For this point-atom offset, we loop over ALL PAOs, and accumulate each PAO multiplied by
                      ! the appropriate coefficient in a local store before putting it in place later
                      ! Actually, the zeta should be *inside* the l loop
                      count1 = 1
                      do l1 = 0,supports_on_atom(this_atom)%lmax
                         do acz = 1,supports_on_atom(this_atom)%naczs(l1)
                            do m1=-l1,l1
                               call pao_elem_derivative_2(direction,the_species,l1,acz,m1,x,y,z,val)
                               if(flag_one_to_one) then
                                  sfsum(count1) = sfsum(count1) + val
                               else
                                  do nsf1 = 1,supports_on_atom(this_atom)%nsuppfuncs
                                     sfsum(nsf1) = sfsum(nsf1) + &
                                          supports_on_atom(this_atom)%supp_func(nsf1)%coefficients(count1)* val
                                  end do ! nsf1
                               end if
                               count1 = count1+1
                            end do ! m1
                         end do ! acz
                      end do ! l1
                      ! Now store data
                      do nsf1 = 1,this_nsf
                         gridfunctions(support)%griddata(position+(nsf1-1)*n_pts_in_block) = sfsum(nsf1)
                      end do ! nsf
                   enddo ! ip=1,npoint
                endif! (npoint > 0) then
                call start_timer(tmr_std_allocation)
                deallocate(sfsum)
                call stop_timer(tmr_std_allocation)
                no_of_ib_ia = no_of_ib_ia + this_nsf*n_pts_in_block
             enddo ! naba_atoms
          enddo ! naba_part
       endif !(naba_atoms_of_blocks(sf)%no_of_part(iblock) > 0) !naba atoms?
    enddo ! iblock : primary set of blocks
    call start_timer(tmr_std_allocation)
    deallocate(ip_store,x_store,y_store,z_store, r_store, STAT=stat)
    if(stat /= 0) call cq_abort(' Error in deallocation in PAO_to_grad_global')
    call stop_timer(tmr_std_allocation)
    call my_barrier()
    return
  end subroutine PAO_to_grad_global
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
