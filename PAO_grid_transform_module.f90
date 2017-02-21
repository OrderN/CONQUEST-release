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
!!   2016/11/10 16:15 nakata
!!    Added subroutine "single_PAO_to_grad"
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

!MS2 !!****f* PAO_grid_transform_module/PAO_to_grid *
!MS2 !!
!MS2 !!  NAME 
!MS2 !!   PAO_to_grid
!MS2 !!  USAGE
!MS2 !! 
!MS2 !!  PURPOSE
!MS2 !!   Takes support functions represented as PAOs and returns their values on the integration
!MS2 !!   grid.
!MS2 !!  INPUTS
!MS2 !!   integer,intent(in) :: myid
!MS2 !!   integer,intent(in) :: NSF, SUPPORT_SIZE
!MS2 !!
!MS2 !!   real(double),intent(out):: support(SUPPORT_SIZE*NSF)
!MS2 !!  USES
!MS2 !!   datatypes, numbers, primary_module, GenComms, blip_grid_transform_module
!MS2 !!  AUTHOR
!MS2 !!   D.R.Bowler
!MS2 !!  CREATION DATE
!MS2 !!   16:40, 2003/09/22
!MS2 !!  MODIFICATION HISTORY
!MS2 !!   11:54, 2003/10/03 dave
!MS2 !!    Changed writes, removed stop
!MS2 !!   2006/06/19 07:57 dave
!MS2 !!    Made variable-NSF compatible
!MS2 !!   2008/05/28 ast
!MS2 !!    Added timers
!MS2 !!   2016/08/09 14:00 nakata
!MS2 !!    Renamed support -> pao_fns
!MS2 !!  SOURCE
!MS2 !!
!MS2   subroutine PAO_to_grid(myid, pao_fns)
!MS2 
!MS2     use datatypes
!MS2     use numbers, ONLY: zero
!MS2     use primary_module, ONLY: bundle
!MS2     use GenComms, ONLY: my_barrier, cq_abort, mtime
!MS2     use blip_grid_transform_module, ONLY: distribute_result
!MS2     use global_module, ONLY: iprint_basis, IPRINT_TIME_THRES2
!MS2     use functions_on_grid, ONLY: gridfunctions
!MS2     use species_module, ONLY: nsf_species
!MS2     use support_spec_format, ONLY: flag_paos_atoms_in_cell
!MS2     use timer_module
!MS2 
!MS2     implicit none
!MS2 
!MS2     integer,intent(in) :: myid
!MS2 
!MS2     integer :: pao_fns
!MS2 
!MS2     integer :: iprim, nsf_send
!MS2     integer :: ii
!MS2     real(double) :: now, nthen
!MS2     type(cq_timer) :: tmr_l_tmp1
!MS2 
!MS2     call start_timer(tmr_std_basis)
!MS2     call start_timer(tmr_l_tmp1,WITH_LEVEL)
!MS2     gridfunctions(pao_fns)%griddata = zero
!MS2 
!MS2     nthen = mtime()
!MS2     if(bundle%mx_iprim < 1) call cq_abort("PAO_to_grid: no primary atoms ! ",myid,bundle%mx_iprim)
!MS2     if(flag_paos_atoms_in_cell) then
!MS2        call PAO_to_grid_global(myid,pao_fns)
!MS2     else
!MS2        do iprim = 1,bundle%mx_iprim
!MS2           call my_barrier()
!MS2           if( iprim <= bundle%n_prim ) then
!MS2              nsf_send = nsf_species(bundle%species(iprim))
!MS2              call do_PAO_transform(iprim,nsf_send)
!MS2              if(iprint_basis>1) then
!MS2                 now = mtime()
!MS2                 if(myid==0) write(io_lun,fmt='(2x,"Proc ",i5," Time for PAO_grid transform (in ms): ",f15.8)') myid,now-nthen
!MS2                 nthen = now
!MS2              end if
!MS2           else
!MS2              nsf_send = 0
!MS2           endif
!MS2           call distribute_result(myid,iprim,nsf_send,pao_fns)
!MS2           if(iprint_basis>1) then
!MS2              now = mtime()
!MS2              if(myid==0) write(io_lun,fmt='(2x,"Proc ",i5," Time for PAO_grid distribute (in ms): ",f15.8)') myid,now - nthen
!MS2              nthen = now
!MS2           end if
!MS2        end do
!MS2     end if
!MS2     call stop_print_timer(tmr_l_tmp1,"PAO to grid transform",IPRINT_TIME_THRES2)
!MS2     call stop_timer(tmr_std_basis)
!MS2     return
!MS2   end subroutine PAO_to_grid
!MS2 !!***
!MS2 
!MS2 !!****f* PAO_grid_transform_module/do_PAO_transform *
!MS2 !!
!MS2 !!  NAME 
!MS2 !!   do_PAO_transform
!MS2 !!  USAGE
!MS2 !! 
!MS2 !!  PURPOSE
!MS2 !!   Does the actual transform onto a local grid
!MS2 !!  INPUTS
!MS2 !! 
!MS2 !! 
!MS2 !!  USES
!MS2 !! 
!MS2 !!  AUTHOR
!MS2 !!   T.Miyazaki/D.R.Bowler
!MS2 !!  CREATION DATE
!MS2 !!   21/6/2000/16:41, 2003/09/22
!MS2 !!  MODIFICATION HISTORY
!MS2 !!   16:42, 2003/09/22 dave
!MS2 !!    Copied and modified do_blip_transform
!MS2 !!   11:54, 2003/10/03 dave
!MS2 !!    Changed writes, removed stop
!MS2 !!   13:18, 03/10/2003 drb 
!MS2 !!    Bug fix - added cq_abort to use GenComms
!MS2 !!   16:01, 30/10/2003 drb 
!MS2 !!    Changed calls, now works
!MS2 !!   07:57, 2004/07/27 dave
!MS2 !!    Changed to use new PAO/SF structure
!MS2 !!   2008/05/28 ast
!MS2 !!    Added timers
!MS2 !!   2009/07/08 16:48 dave
!MS2 !!    Added code for one-to-one PAO to SF assignment
!MS2 !!   2016/07/14 16:30 nakata
!MS2 !!    Renamed naba_blk_supp -> naba_blocks_of_atoms
!MS2 !!   2016/07/29 18:30 nakata
!MS2 !!    Renamed supports_on_atom -> blips_on_atom
!MS2 !!   2016/08/05 10:27 dave
!MS2 !!    Bug fix: evaluate_pao was passed 1 instead of acz
!MS2 !!  SOURCE
!MS2 !!
!MS2   subroutine do_PAO_transform(iprim,this_nsf)
!MS2 
!MS2     use datatypes
!MS2     use numbers
!MS2     use global_module,   ONLY: rcellx,rcelly,rcellz, iprint_basis
!MS2     use group_module,    ONLY: blocks
!MS2     use primary_module,  ONLY: bundle
!MS2     use cover_module,    ONLY: BCS_blocks
!MS2     use set_blipgrid_module, ONLY: naba_blocks_of_atoms
!MS2     use comm_array_module,   ONLY: send_array
!MS2     use block_module,    ONLY: nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
!MS2     use support_spec_format, ONLY: blips_on_atom, flag_paos_atoms_in_cell, flag_one_to_one
!MS2     use GenComms, ONLY: myid, cq_abort
!MS2     use dimens, ONLY: r_h
!MS2     use angular_coeff_routines, ONLY: evaluate_pao
!MS2 
!MS2     implicit none
!MS2 
!MS2     integer,intent(in)     :: iprim
!MS2     integer,intent(in)     :: this_nsf
!MS2 
!MS2     ! Local variables
!MS2     real(double):: rec_sgs, dx, dy, dz
!MS2     integer     :: x, y, z
!MS2     integer     :: nsf1
!MS2 
!MS2     integer     :: nblkx,nblky,nblkz  ! nos. of integ. pts in a block 
!MS2     ! along x,y and z direction
!MS2     integer     :: naba_blk,nx_blk,ny_blk,nz_blk,igrid
!MS2     real(double) :: dcellx_block,dcelly_block,dcellz_block
!MS2     real(double) :: dcellx_grid, dcelly_grid, dcellz_grid
!MS2     real(double):: dsum(this_nsf), val
!MS2 
!MS2     integer     :: ix,iy,iz,acz, l1, m1, count1
!MS2     integer     :: ind_blk,ind_blip
!MS2     integer     :: nxmin_grid,nxmax_grid,nymin_grid,nymax_grid,nzmin_grid,nzmax_grid
!MS2     integer     :: ncover_yz, atom_species, this_atom, stat
!MS2 
!MS2     ! Species of atom
!MS2     atom_species = bundle%species(iprim)
!MS2     ! How are we storing PAO coefficients ?
!MS2     if(flag_one_to_one) then
!MS2        this_atom = atom_species
!MS2     else
!MS2        if(flag_paos_atoms_in_cell) then
!MS2           this_atom = bundle%ig_prim(iprim)
!MS2        else
!MS2           this_atom = iprim
!MS2        end if
!MS2     end if
!MS2     nblkx=nx_in_block
!MS2     nblky=ny_in_block
!MS2     nblkz=nz_in_block
!MS2 
!MS2     dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nblkx
!MS2     dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/nblky
!MS2     dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nblkz
!MS2 
!MS2     ncover_yz=BCS_blocks%ncovery*BCS_blocks%ncoverz
!MS2     igrid=0
!MS2 
!MS2     call start_timer(tmr_std_allocation)
!MS2     allocate(send_array(naba_blocks_of_atoms%no_naba_blk(iprim)*n_pts_in_block*this_nsf),STAT=stat)
!MS2     if(stat/=0) call cq_abort("Error allocating send_array in do_PAO_grid: ",&
!MS2          naba_blocks_of_atoms%no_naba_blk(iprim)*n_pts_in_block*this_nsf)
!MS2     call stop_timer(tmr_std_allocation)
!MS2     send_array(:) = zero
!MS2 
!MS2     ! Loop over neighbour blocks
!MS2     DO naba_blk=1,naba_blocks_of_atoms%no_naba_blk(iprim)! naba blks in NOPG order
!MS2        ! iprim : primary seq. no. of the atom
!MS2        ind_blk=naba_blocks_of_atoms%list_naba_blk(naba_blk,iprim) !CC in BCS_blocks
!MS2        ind_blk=ind_blk-1
!MS2 
!MS2        nx_blk=ind_blk/ncover_yz
!MS2        ny_blk=(ind_blk-nx_blk*ncover_yz)/BCS_blocks%ncoverz
!MS2        nz_blk=ind_blk-nx_blk*ncover_yz-ny_blk*BCS_blocks%ncoverz
!MS2 
!MS2        nx_blk=nx_blk-BCS_blocks%nspanlx+BCS_blocks%nx_origin
!MS2        ny_blk=ny_blk-BCS_blocks%nspanly+BCS_blocks%ny_origin
!MS2        nz_blk=nz_blk-BCS_blocks%nspanlz+BCS_blocks%nz_origin
!MS2 
!MS2        do iz=1,nblkz    ! z-direction in a block
!MS2           z=nblkz*(nz_blk-1)+iz-1
!MS2           do iy=1,nblky    ! y-direction in a block
!MS2              y=nblky*(ny_blk-1)+iy-1
!MS2              do ix=1,nblkx    ! x-direction in a block
!MS2                 x=nblkx*(nx_blk-1)+ix-1
!MS2                 ! I want to find the position of the grid point relative to the atom in its primary set.
!MS2                 ! The grid-point position _would_ be x + xadd, with xadd = origin + cover block x - FSC block x
!MS2                 ! cover block x is nx_blk BEFORE the manipulation above
!MS2                 ! 
!MS2                 ! Actually, isn't nx_blk AFTER the manipulation just what I want ? We subtract the left-span 
!MS2                 ! and add the origin
!MS2                 ! 
!MS2                 ! OK - let's assume that x, y and z are EXACTLY what we want for the grid point.
!MS2                 igrid=igrid+1  ! seq. no. of integration grids
!MS2                 !if(this_nsf*igrid > mx_send_real) call cq_abort("do_PAO_transform: overflow ",igrid,mx_send_real)
!MS2                 ! generate the offset vector                
!MS2                 dx = x*dcellx_grid - bundle%xprim(iprim)
!MS2                 dy = y*dcelly_grid - bundle%yprim(iprim)
!MS2                 dz = z*dcellz_grid - bundle%zprim(iprim)
!MS2                 if(dx*dx+dy*dy+dz*dz<r_h*r_h) then
!MS2                    ! loop over PAOs
!MS2                    dsum = zero
!MS2                    count1 = 1
!MS2                    do l1 = 0,blips_on_atom(this_atom)%lmax
!MS2                       do acz = 1,blips_on_atom(this_atom)%naczs(l1)
!MS2                          do m1=-l1,l1
!MS2                             call evaluate_pao(atom_species,l1,acz,m1,dx,dy,dz,val)
!MS2                             if(flag_one_to_one) then
!MS2                                   dsum(count1) = dsum(count1) + val
!MS2                             else
!MS2                                do nsf1 = 1,blips_on_atom(this_atom)%nsuppfuncs
!MS2                                   dsum(nsf1) = dsum(nsf1) + &
!MS2                                        blips_on_atom(this_atom)%supp_func(nsf1)%coefficients(count1)* val
!MS2                                end do ! nsf
!MS2                             end if
!MS2                             count1 = count1+1
!MS2                          end do ! m1
!MS2                       end do ! acz
!MS2                    end do ! l1
!MS2                    ! store in send_array
!MS2                    do nsf1 = 1,this_nsf
!MS2                       send_array(this_nsf*igrid-this_nsf+nsf1) = dsum(nsf1)
!MS2                    enddo ! nsf1=this_nsf
!MS2                 end if ! if(dx*dx...<r_h*r_h)
!MS2              enddo    ! x-direction in a block
!MS2           enddo    ! y-direction in a block
!MS2        enddo    ! z-direction in a block
!MS2     end do    ! Loop over naba blocks for this primary atom
!MS2     return
!MS2   end subroutine do_PAO_transform
!MS2 !!***
!MS2 
!MS2 !!****f* PAO_grid_transform_module/PAO_to_grad *
!MS2 !!
!MS2 !!  NAME 
!MS2 !!   PAO_to_grad
!MS2 !!  USAGE
!MS2 !! 
!MS2 !!  PURPOSE
!MS2 !!   Takes support functions represented as PAOs and returns their differentials in a 
!MS2 !!   specified direction on the integration grid.
!MS2 !!  INPUTS
!MS2 !!   integer,intent(in) :: myid
!MS2 !!   integer,intent(in) :: direction
!MS2 !!   integer,intent(in) :: NSF, SUPPORT_SIZE
!MS2 !!
!MS2 !!   real(double),intent(out):: support(SUPPORT_SIZE*NSF)
!MS2 !!  USES
!MS2 !!   datatypes, numbers, primary_module, GenComms, blip_grid_transform_module
!MS2 !!  AUTHOR
!MS2 !!   D.R.Bowler
!MS2 !!  CREATION DATE
!MS2 !!   08:01, 2003/10/01
!MS2 !!  MODIFICATION HISTORY
!MS2 !!   11:53, 2003/10/03 dave
!MS2 !!    Changed write statements, removed stop
!MS2 !!   02:30, 2003/11/13 dave
!MS2 !!    Added direction to call to do_PAO_grad
!MS2 !!   2008/05/28 ast
!MS2 !!    Added timers
!MS2 !!   2016/08/09 14:00 nakata
!MS2 !!    Renamed support -> pao_fns
!MS2 !!  SOURCE
!MS2 !!
!MS2   subroutine PAO_to_grad(myid, direction, pao_fns)
!MS2 
!MS2     use datatypes
!MS2     use numbers, ONLY: zero
!MS2     use global_module, ONLY: iprint_basis
!MS2     use primary_module, ONLY: bundle
!MS2     use GenComms, ONLY: my_barrier, cq_abort, mtime
!MS2     use blip_grid_transform_module, ONLY: distribute_result
!MS2     use functions_on_grid, ONLY: gridfunctions
!MS2     use species_module, ONLY: nsf_species
!MS2     use support_spec_format, ONLY: flag_paos_atoms_in_cell
!MS2 
!MS2     implicit none
!MS2 
!MS2     integer,intent(in) :: myid, direction
!MS2     integer :: pao_fns
!MS2 
!MS2     integer :: iprim, nsf_send
!MS2     integer :: ii
!MS2     real(double) :: now, nthen
!MS2 
!MS2     call start_timer(tmr_std_basis)
!MS2     gridfunctions(pao_fns)%griddata = zero
!MS2 
!MS2     nthen = mtime()
!MS2     if(bundle%mx_iprim < 1) call cq_abort("PAO_to_grad: no primary atoms ! ",myid,bundle%mx_iprim)
!MS2     if(flag_paos_atoms_in_cell) then
!MS2        call PAO_to_grad_global(myid,direction,pao_fns)
!MS2     else
!MS2        do iprim = 1,bundle%mx_iprim
!MS2           call my_barrier()
!MS2           if( iprim <= bundle%n_prim ) then
!MS2              nsf_send = nsf_species(bundle%species(iprim))
!MS2              call do_PAO_grad_transform(direction,iprim,nsf_send)
!MS2              if(iprint_basis>1) then
!MS2                 now = mtime()
!MS2                 if(myid==0) write(io_lun,fmt='(2x,"Proc ",i5," Time for PAO_grad transform (in ms): ",f15.8)') myid,now - nthen
!MS2                 nthen = now
!MS2              end if
!MS2           else
!MS2              nsf_send = 0
!MS2           endif
!MS2           call distribute_result(myid,iprim,nsf_send,pao_fns)
!MS2           if(iprint_basis>1) then
!MS2              now = mtime()
!MS2              if(myid==0) write(io_lun,fmt='(2x,"Proc ",i5," Time for PAO_grad distribute (in ms): ",f15.8)') myid,now - nthen
!MS2              nthen = now
!MS2           end if
!MS2        end do
!MS2     end if
!MS2     call stop_timer(tmr_std_basis)
!MS2     return
!MS2   end subroutine PAO_to_grad
!MS2 !!***
!MS2 
!MS2 !!****f* PAO_grid_transform_module/do_PAO_grad_transform *
!MS2 !!
!MS2 !!  NAME 
!MS2 !!   do_PAO_grad_transform
!MS2 !!  USAGE
!MS2 !! 
!MS2 !!  PURPOSE
!MS2 !!   Does the actual transform of the derivative of the PAO onto a local grid
!MS2 !!  INPUTS
!MS2 !! 
!MS2 !! 
!MS2 !!  USES
!MS2 !! 
!MS2 !!  AUTHOR
!MS2 !!   D.R.Bowler
!MS2 !!  CREATION DATE
!MS2 !!   11:53, 2003/10/03 
!MS2 !!  MODIFICATION HISTORY
!MS2 !!   16:42, 2003/09/22 dave
!MS2 !!    Copied and modified do_blip_transform
!MS2 !!   11:53, 2003/10/03 dave
!MS2 !!    Changed write statements, removed stop
!MS2 !!   16:08, 30/10/2003 drb 
!MS2 !!    Brought up-to-date with do_PAO_grid_transform: still needs DIRECTION !
!MS2 !!   02:30, 2003/11/13 dave
!MS2 !!    Added direction
!MS2 !!   07:57, 2004/07/27 dave
!MS2 !!    Changed to use new PAO/SF structure
!MS2 !!   2008/05/28 ast
!MS2 !!    Added timers
!MS2 !!   2009/07/08 16:48 dave
!MS2 !!    Added code for one-to-one PAO to SF assignment
!MS2 !!   2016/07/14 16:30 nakata
!MS2 !!    Renamed naba_blk_supp -> naba_blocks_of_atoms
!MS2 !!   2016/07/29 18:30 nakata
!MS2 !!    Renamed supports_on_atom -> blips_on_atom
!MS2 !!  SOURCE
!MS2 !!
!MS2   subroutine do_PAO_grad_transform(direction,iprim,this_nsf)
!MS2 
!MS2     use datatypes
!MS2     use numbers
!MS2     use global_module,   ONLY: rcellx,rcelly,rcellz,iprint_basis
!MS2     use group_module,    ONLY: blocks
!MS2     use primary_module,  ONLY: bundle
!MS2     use cover_module,    ONLY: BCS_blocks
!MS2     use set_blipgrid_module, ONLY: naba_blocks_of_atoms
!MS2     use comm_array_module,   ONLY: send_array
!MS2     use block_module,    ONLY: nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
!MS2     use support_spec_format, ONLY: blips_on_atom, flag_paos_atoms_in_cell, flag_one_to_one
!MS2     use GenComms, ONLY: myid, cq_abort
!MS2     use dimens, ONLY: r_h
!MS2     use angular_coeff_routines, ONLY: pao_elem_derivative_2, evaluate_pao
!MS2 
!MS2     implicit none
!MS2 
!MS2     integer,intent(in)     :: iprim, direction
!MS2     integer,intent(in)     :: this_nsf
!MS2 
!MS2     ! Local variables
!MS2     real(double):: rec_sgs, dx, dy, dz
!MS2     integer     :: x, y, z
!MS2     integer     :: nsf1
!MS2 
!MS2     integer     :: nblkx,nblky,nblkz  ! nos. of integ. pts in a block 
!MS2     ! along x,y and z direction
!MS2     integer     :: naba_blk,nx_blk,ny_blk,nz_blk,igrid
!MS2     real(double) :: dcellx_block,dcelly_block,dcellz_block
!MS2     real(double) :: dcellx_grid, dcelly_grid, dcellz_grid
!MS2     real(double):: dsum(this_nsf),val,val2,pao_val,pao_val2
!MS2 
!MS2     integer     :: ix,iy,iz,acz,l1,m1,count1
!MS2     integer     :: ind_blk,ind_blip
!MS2     integer     :: nxmin_grid,nxmax_grid,nymin_grid,nymax_grid,nzmin_grid,nzmax_grid
!MS2     integer     :: ncover_yz, atom_species, this_atom, stat
!MS2 
!MS2     ! Species of atom
!MS2     atom_species = bundle%species(iprim)
!MS2     ! How are we storing PAO coefficients ?
!MS2     if(flag_one_to_one) then
!MS2        this_atom = atom_species
!MS2     else
!MS2        if(flag_paos_atoms_in_cell) then
!MS2           this_atom = bundle%ig_prim(iprim)
!MS2        else
!MS2           this_atom = iprim
!MS2        end if
!MS2     end if
!MS2     nblkx=nx_in_block
!MS2     nblky=ny_in_block
!MS2     nblkz=nz_in_block
!MS2 
!MS2     dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nblkx
!MS2     dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/nblky
!MS2     dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nblkz
!MS2 
!MS2     ncover_yz=BCS_blocks%ncovery*BCS_blocks%ncoverz
!MS2     igrid=0
!MS2 
!MS2     call start_timer(tmr_std_allocation)
!MS2     allocate(send_array(naba_blocks_of_atoms%no_naba_blk(iprim)*n_pts_in_block*this_nsf),&
!MS2              STAT=stat)
!MS2     if(stat/=0) call cq_abort("Error allocating send_array in do_PAO_grad: ",&
!MS2          naba_blocks_of_atoms%no_naba_blk(iprim)*n_pts_in_block*this_nsf)
!MS2     call stop_timer(tmr_std_allocation)
!MS2     send_array(:) = zero
!MS2 
!MS2     ! Loop over neighbour blocks
!MS2     DO naba_blk=1,naba_blocks_of_atoms%no_naba_blk(iprim)! naba blks in NOPG order
!MS2        ! iprim : primary seq. no. of the atom
!MS2        ind_blk=naba_blocks_of_atoms%list_naba_blk(naba_blk,iprim) !CC in BCS_blocks
!MS2        ind_blk=ind_blk-1
!MS2 
!MS2        nx_blk=ind_blk/ncover_yz
!MS2        ny_blk=(ind_blk-nx_blk*ncover_yz)/BCS_blocks%ncoverz
!MS2        nz_blk=ind_blk-nx_blk*ncover_yz-ny_blk*BCS_blocks%ncoverz
!MS2 
!MS2        nx_blk=nx_blk-BCS_blocks%nspanlx+BCS_blocks%nx_origin
!MS2        ny_blk=ny_blk-BCS_blocks%nspanly+BCS_blocks%ny_origin
!MS2        nz_blk=nz_blk-BCS_blocks%nspanlz+BCS_blocks%nz_origin
!MS2 
!MS2        do iz=1,nblkz    ! z-direction in a block
!MS2           z=nblkz*(nz_blk-1)+iz-1
!MS2           do iy=1,nblky    ! y-direction in a block
!MS2              y=nblky*(ny_blk-1)+iy-1
!MS2              do ix=1,nblkx    ! x-direction in a block
!MS2                 x=nblkx*(nx_blk-1)+ix-1
!MS2                 ! I want to find the position of the grid point relative to the atom in its primary set.
!MS2                 ! The grid-point position _would_ be x + xadd, with xadd = origin + cover block x - FSC block x
!MS2                 ! cover block x is nx_blk BEFORE the manipulation above
!MS2                 ! 
!MS2                 ! Actually, isn't nx_blk AFTER the manipulation just what I want ? We subtract the left-span and add the origin
!MS2                 ! 
!MS2                 ! OK - let's assume that x, y and z are EXACTLY what we want for the grid point.
!MS2                 igrid=igrid+1  ! seq. no. of integration grids
!MS2                 if(this_nsf*igrid > naba_blocks_of_atoms%no_naba_blk(iprim)*n_pts_in_block*this_nsf) &
!MS2                      call cq_abort("do_PAO_grad_transform: overflow ",this_nsf*igrid, &
!MS2                      naba_blocks_of_atoms%no_naba_blk(iprim)*n_pts_in_block*this_nsf)
!MS2                 ! generate the offset vector                
!MS2                 dx = x*dcellx_grid - bundle%xprim(iprim)
!MS2                 dy = y*dcelly_grid - bundle%yprim(iprim)
!MS2                 dz = z*dcellz_grid - bundle%zprim(iprim)
!MS2                 if(dx*dx+dy*dy+dz*dz<r_h*r_h) then
!MS2                    ! loop over PAOs
!MS2                    dsum = zero
!MS2                    count1 = 1
!MS2                    do l1 = 0,blips_on_atom(this_atom)%lmax
!MS2                       do acz = 1,blips_on_atom(this_atom)%naczs(l1)
!MS2                          do m1=-l1,l1
!MS2                             call pao_elem_derivative_2(direction,atom_species,l1,acz,m1,dx,dy,dz,val)
!MS2                             if(flag_one_to_one) then
!MS2                                dsum(count1) = dsum(count1) + val
!MS2                             else
!MS2                                do nsf1 = 1,blips_on_atom(this_atom)%nsuppfuncs
!MS2                                   dsum(nsf1) = dsum(nsf1) + &
!MS2                                        blips_on_atom(this_atom)%supp_func(nsf1)%coefficients(count1)* val
!MS2                                end do ! nsf
!MS2                             end if
!MS2                             count1 = count1+1
!MS2                          end do ! m1
!MS2                       end do ! acz
!MS2                    end do ! l1
!MS2                    ! store in send_array
!MS2                    do nsf1 = 1,this_nsf
!MS2                       send_array(this_nsf*igrid-this_nsf+nsf1) = dsum(nsf1)
!MS2                    enddo
!MS2                 end if
!MS2              enddo    ! x-direction in a block
!MS2           enddo    ! y-direction in a block
!MS2        enddo    ! z-direction in a block
!MS2     end do    ! Loop over naba blocks for this primary atom
!MS2     return
!MS2   end subroutine do_PAO_grad_transform
!MS2 !!***
!MS2 
!MS2 !!****f* PAO_grid_transform_module/PAO_to_grid_global *
!MS2 !!
!MS2 !!  NAME 
!MS2 !!   PAO_to_grid_global
!MS2 !!  USAGE
!MS2 !! 
!MS2 !!  PURPOSE
!MS2 !!   Transforms PAOs to grid assuming global PAO coefficient storage: no communication involved
!MS2 !!  INPUTS
!MS2 !!   myid: processor number
!MS2 !!   pao_fns: index for grid function
!MS2 !!  USES
!MS2 !! 
!MS2 !!  AUTHOR
!MS2 !!   D. R. Bowler
!MS2 !!  CREATION DATE
!MS2 !!   2004 sometime probably
!MS2 !!  MODIFICATION HISTORY
!MS2 !!   2006/06/20 08:24 dave
!MS2 !!    Various changes for variable NSF
!MS2 !!   2008/05/28 ast
!MS2 !!    Added timers
!MS2 !!   2009/07/08 16:48 dave
!MS2 !!    Added code for one-to-one PAO to SF assignment
!MS2 !!   2016/07/20 16:30 nakata
!MS2 !!    Renamed naba_atm -> naba_atoms_of_blocks
!MS2 !!   2016/07/29 18:30 nakata
!MS2 !!    Renamed supports_on_atom -> blips_on_atom
!MS2 !!   2016/08/01 17:30 nakata
!MS2 !!    Introduced atomf instead of sf
!MS2 !!   2016/08/09 14:00 nakata
!MS2 !!    Renamed support -> pao_fns
!MS2 !!  SOURCE
!MS2 !!
!MS2   subroutine PAO_to_grid_global(myid, pao_fns)
!MS2 
!MS2     use datatypes
!MS2     use primary_module, ONLY: bundle
!MS2     use GenComms, ONLY: my_barrier, cq_abort, mtime
!MS2     use dimens, ONLY: r_h
!MS2     use GenComms, ONLY: inode, ionode, cq_abort
!MS2     use numbers
!MS2     use global_module, ONLY: rcellx,rcelly,rcellz,id_glob,ni_in_cell, iprint_basis, species_glob, atomf
!MS2     use species_module, ONLY: species, nsf_species
!MS2     !  At present, these arrays are dummy arguments.
!MS2     use block_module, ONLY : nx_in_block,ny_in_block,nz_in_block, &
!MS2          n_pts_in_block
!MS2     use group_module, ONLY : blocks, parts
!MS2     use primary_module, ONLY: domain
!MS2     use cover_module, ONLY: DCS_parts
!MS2     use set_blipgrid_module, ONLY : naba_atoms_of_blocks
!MS2     use angular_coeff_routines, ONLY : evaluate_pao
!MS2     use support_spec_format, ONLY: blips_on_atom, flag_paos_atoms_in_cell, flag_one_to_one
!MS2     use functions_on_grid, ONLY: gridfunctions, fn_on_grid
!MS2 
!MS2     implicit none 
!MS2     integer,intent(in) :: myid
!MS2     integer :: pao_fns
!MS2 
!MS2 
!MS2     !local
!MS2     real(double):: dcellx_block,dcelly_block,dcellz_block
!MS2     integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom, this_atom
!MS2     real(double):: xatom,yatom,zatom,alpha,step
!MS2     real(double):: xblock,yblock,zblock
!MS2     integer :: the_species
!MS2     integer :: j,iblock,the_l,ipoint, igrid
!MS2     real(double) :: r_from_i
!MS2     real(double) :: rr,a,b,c,d,x,y,z,nl_potential
!MS2     integer :: no_of_ib_ia, offset_position
!MS2     integer :: position,iatom
!MS2     integer :: stat, nl, npoint, ip, this_nsf
!MS2     integer :: i,m, m1min, m1max,acz,m1,l1,nsf1,atom_species, count1
!MS2     integer     , allocatable :: ip_store(:)
!MS2     real(double), allocatable :: x_store(:)
!MS2     real(double), allocatable :: y_store(:)
!MS2     real(double), allocatable :: z_store(:)
!MS2     real(double), allocatable :: r_store(:)
!MS2     real(double) :: coulomb_energy
!MS2     real(double) :: rcut
!MS2     real(double) :: r1, r2, r3, r4, core_charge, gauss_charge
!MS2     real(double) :: val, theta, phi, r_tmp
!MS2     real(double), dimension(:), allocatable :: sfsum
!MS2 
!MS2     if(.NOT.flag_paos_atoms_in_cell) call cq_abort("Calling global PAO-to-grid with local PAO storage")
!MS2     call my_barrier()
!MS2 
!MS2     call start_timer(tmr_std_allocation)
!MS2     allocate(ip_store(n_pts_in_block),x_store(n_pts_in_block),y_store(n_pts_in_block),z_store(n_pts_in_block), &
!MS2          r_store(n_pts_in_block))
!MS2     call stop_timer(tmr_std_allocation)
!MS2     ! --  Start of subroutine  ---
!MS2 
!MS2     dcellx_block=rcellx/blocks%ngcellx
!MS2     dcelly_block=rcelly/blocks%ngcelly
!MS2     dcellz_block=rcellz/blocks%ngcellz
!MS2 
!MS2     no_of_ib_ia = 0
!MS2     gridfunctions(pao_fns)%griddata = zero
!MS2 
!MS2     ! loop arround grid points in the domain, and for each
!MS2     ! point, get the core charge (local charge) density from the atoms 
!MS2     ! whose distances from the grid point are within the cutoff. 
!MS2 
!MS2     do iblock = 1, domain%groups_on_node ! primary set of blocks
!MS2        xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
!MS2        yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
!MS2        zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
!MS2        if(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) then ! if there are naba atoms
!MS2           iatom=0
!MS2           do ipart=1,naba_atoms_of_blocks(atomf)%no_of_part(iblock)
!MS2              jpart=naba_atoms_of_blocks(atomf)%list_part(ipart,iblock)
!MS2              if(jpart > DCS_parts%mx_gcover) then 
!MS2                 call cq_abort('PAO_to_grid_global: JPART ERROR ',ipart,jpart)
!MS2              endif
!MS2              ind_part=DCS_parts%lab_cell(jpart)
!MS2              do ia=1,naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart,iblock)
!MS2                 iatom=iatom+1
!MS2                 ii = naba_atoms_of_blocks(atomf)%list_atom(iatom,iblock)
!MS2                 icover= DCS_parts%icover_ibeg(jpart)+ii-1
!MS2                 ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)
!MS2 
!MS2                 if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
!MS2                    call cq_abort('PAO_to_grid_global: globID ERROR ', &
!MS2                         ii,parts%icell_beg(ind_part))
!MS2                 endif
!MS2                 if(icover > DCS_parts%mx_mcover) then
!MS2                    call cq_abort('PAO_to_grid_global: icover ERROR ', &
!MS2                         icover,DCS_parts%mx_mcover)
!MS2                 endif
!MS2 
!MS2                 !-- I (TM) am not sure where the next line should be.
!MS2                 !   (no_of_ib_ia = no_of_ib_ia +1)
!MS2                 !  This determines how to store 
!MS2                 !   pseudofunctions(ipoints,ncf,naba_atm,iblock).
!MS2                 !  Now, I assume we should consider all naba_atm whose 
!MS2                 ! distance from the block is within the maximum of core_radius.
!MS2                 ! This is needed to keep the consistency with <set_bucket>.
!MS2                 ! However, we can change this strategy by changing naba_atoms_of_blocks(atomf).
!MS2                 !no_of_ib_ia = no_of_ib_ia +1
!MS2 
!MS2                 xatom=DCS_parts%xcover(icover)
!MS2                 yatom=DCS_parts%ycover(icover)
!MS2                 zatom=DCS_parts%zcover(icover)
!MS2 
!MS2                 the_species=species_glob(ig_atom)
!MS2                 this_nsf = nsf_species(the_species)
!MS2                 call start_timer(tmr_std_allocation)
!MS2                 allocate(sfsum(this_nsf))
!MS2                 call stop_timer(tmr_std_allocation)
!MS2                 
!MS2                 !calculates distances between the atom and integration grid points
!MS2                 !in the block and stores which integration grids are neighbours.
!MS2                 rcut = r_h + RD_ERR
!MS2                 call check_block (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
!MS2                      npoint,ip_store,r_store,x_store,y_store,z_store,n_pts_in_block) !out
!MS2                 !r_from_i = sqrt((xatom-xblock)**2+(yatom-yblock)**2+ &
!MS2                 !     (zatom-zblock)**2 )
!MS2                 if(flag_one_to_one) then
!MS2                    this_atom = the_species
!MS2                 else
!MS2                    this_atom = ig_atom
!MS2                 end if
!MS2                 if(npoint > 0) then
!MS2                    !offset_position = (no_of_ib_ia-1) * this_nsf * n_pts_in_block
!MS2                    offset_position = no_of_ib_ia
!MS2                    do ip=1,npoint
!MS2                       ipoint=ip_store(ip)
!MS2                       position= offset_position + ipoint
!MS2                       if(position > gridfunctions(pao_fns)%size) &
!MS2                            call cq_abort('PAO_to_grid_global: position error 1 ', &
!MS2                            position, gridfunctions(pao_fns)%size)
!MS2                       r_from_i = r_store(ip)
!MS2                       x = x_store(ip)
!MS2                       y = y_store(ip)
!MS2                       z = z_store(ip)
!MS2                       sfsum = zero
!MS2                       ! For this point-atom offset, we loop over ALL PAOs, and accumulate each PAO multiplied by
!MS2                       ! the appropriate coefficient in a local store before putting it in place later
!MS2                       count1 = 1
!MS2                       do l1 = 0,blips_on_atom(this_atom)%lmax
!MS2                          do acz = 1,blips_on_atom(this_atom)%naczs(l1)
!MS2                             do m1=-l1,l1
!MS2                                call evaluate_pao(the_species,l1,acz,m1,x,y,z,val)
!MS2                                if(flag_one_to_one) then
!MS2                                   sfsum(count1) = sfsum(count1) + val
!MS2                                else
!MS2                                   do nsf1 = 1,blips_on_atom(this_atom)%nsuppfuncs
!MS2                                      sfsum(nsf1) = sfsum(nsf1) + &
!MS2                                           blips_on_atom(this_atom)%supp_func(nsf1)%coefficients(count1)* val
!MS2                                   end do ! nsf1
!MS2                                end if
!MS2                                count1 = count1+1
!MS2                             end do ! m1
!MS2                          end do ! acz
!MS2                       end do ! l1
!MS2                       ! Now store data
!MS2                       if(position+(this_nsf-1)*n_pts_in_block > gridfunctions(pao_fns)%size) &
!MS2                            call cq_abort('PAO_to_grid_global: position error 2 ', &
!MS2                            position+(this_nsf-1)*n_pts_in_block, gridfunctions(pao_fns)%size)
!MS2                       do nsf1 = 1,this_nsf
!MS2                          gridfunctions(pao_fns)%griddata(position+(nsf1-1)*n_pts_in_block) = sfsum(nsf1)
!MS2                       end do ! nsf
!MS2                    enddo ! ip=1,npoint
!MS2                 endif! (npoint > 0) then
!MS2                 call start_timer(tmr_std_allocation)
!MS2                 deallocate(sfsum)
!MS2                 call stop_timer(tmr_std_allocation)
!MS2                 no_of_ib_ia = no_of_ib_ia + this_nsf*n_pts_in_block
!MS2              enddo ! naba_atoms
!MS2           enddo ! naba_part
!MS2        endif !(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) !naba atoms?
!MS2     enddo ! iblock : primary set of blocks
!MS2     call my_barrier()
!MS2     call start_timer(tmr_std_allocation)
!MS2     deallocate(ip_store,x_store,y_store,z_store,r_store)
!MS2     call stop_timer(tmr_std_allocation)
!MS2     return
!MS2   end subroutine PAO_to_grid_global
!MS2 !!***

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
    integer :: stat, nl, npoint, ip, this_nsf
    integer :: i,m, m1min, m1max,acz,m1,l1,count1
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
    gridfunctions(pao_fns)%griddata = zero

    ! loop arround grid points in the domain, and for each
    ! point, get the PAO values

    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
       yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
       zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
       if(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom=0
          do ipart=1,naba_atoms_of_blocks(atomf)%no_of_part(iblock)
             jpart=naba_atoms_of_blocks(atomf)%list_part(ipart,iblock)
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('single_PAO_to_grid: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atoms_of_blocks(atomf)%list_atom(iatom,iblock)
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

                if(npoint > 0) then
                   !offset_position = (no_of_ib_ia-1) * npao * n_pts_in_block
                   offset_position = no_of_ib_ia
                   do ip=1,npoint
                      ipoint=ip_store(ip)
                      position= offset_position + ipoint
                      if(position > gridfunctions(pao_fns)%size) call cq_abort &
                           ('single_pao_to_grid: position error ', position, gridfunctions(pao_fns)%size)

                      r_from_i = r_store(ip)
                      x = x_store(ip)
                      y = y_store(ip)
                      z = z_store(ip)
                      ! For this point-atom offset, we accumulate the PAO on the grid
                      count1 = 1
                      do l1 = 0,pao(the_species)%greatest_angmom
                         do acz = 1,pao(the_species)%angmom(l1)%n_zeta_in_angmom
                            do m1=-l1,l1
                               call evaluate_pao(the_species,l1,acz,m1,x,y,z,val)
                               if(position+(count1-1)*n_pts_in_block > gridfunctions(pao_fns)%size) &
                                    call cq_abort('single_pao_to_grid: position error ', &
                                    position, gridfunctions(pao_fns)%size)
                               gridfunctions(pao_fns)%griddata(position+(count1-1)*n_pts_in_block) = val
                               count1 = count1+1
                            end do ! m1
                         end do ! acz
                      end do ! l1
                   enddo ! ip=1,npoint
                endif! (npoint > 0) then
                no_of_ib_ia = no_of_ib_ia + npao_species(the_species)*n_pts_in_block
             enddo ! naba_atoms
          enddo ! naba_part
       endif !(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) !naba atoms?
    enddo ! iblock : primary set of blocks
    call my_barrier()
    call start_timer(tmr_std_allocation)
    deallocate(ip_store,x_store,y_store,z_store,r_store)
    call stop_timer(tmr_std_allocation)
    call stop_timer(tmr_std_basis)
    return
  end subroutine Single_PAO_to_grid
!!***

!MS2 !!****f* PAO_grid_transform_module/PAO_to_grad_global *
!MS2 !!
!MS2 !!  NAME 
!MS2 !!   PAO_to_grad_Global
!MS2 !!  USAGE
!MS2 !! 
!MS2 !!  PURPOSE
!MS2 !!   Transforms gradient of support functions to grid using PAO basis functions
!MS2 !!   This version assumes all PAO coeffients on all atoms are stored on all processors, 
!MS2 !!   so no communication is required
!MS2 !!  INPUTS
!MS2 !! 
!MS2 !! 
!MS2 !!  USES
!MS2 !! 
!MS2 !!  AUTHOR
!MS2 !!   D. R. Bowler
!MS2 !!  CREATION DATE
!MS2 !!   2004 sometime
!MS2 !!  MODIFICATION HISTORY
!MS2 !!   2006/06/20 08:30 dave
!MS2 !!    Variable NSF preparation, tidying, name change
!MS2 !!   2008/05/28 ast
!MS2 !!    Added timers
!MS2 !!   2009/07/08 16:48 dave
!MS2 !!    Added code for one-to-one PAO to SF assignment
!MS2 !!   2016/07/20 16:30 nakata
!MS2 !!    Renamed naba_atm -> naba_atoms_of_blocks
!MS2 !!   2016/07/29 18:30 nakata
!MS2 !!    Renamed supports_on_atom -> blips_on_atom
!MS2 !!   2016/08/01 17:30 nakata
!MS2 !!    Introduced atomf instead of sf
!MS2 !!   2016/08/09 14:00 nakata
!MS2 !!    Renamed support -> pao_fns
!MS2 !!  SOURCE
!MS2 !!
!MS2   subroutine PAO_to_grad_global(myid, direction, pao_fns)
!MS2 
!MS2     use datatypes
!MS2     use primary_module, ONLY: bundle
!MS2     use GenComms, ONLY: my_barrier, cq_abort, mtime
!MS2     use dimens, ONLY: r_h
!MS2     use GenComms, ONLY: inode, ionode
!MS2     use numbers
!MS2     use global_module, ONLY: rcellx,rcelly,rcellz,id_glob,ni_in_cell, iprint_basis, species_glob, atomf
!MS2     use species_module, ONLY: species, nsf_species
!MS2     !  At present, these arrays are dummy arguments.
!MS2     use block_module, ONLY : nx_in_block,ny_in_block,nz_in_block, &
!MS2          n_pts_in_block
!MS2     use group_module, ONLY : blocks, parts
!MS2     use primary_module, ONLY: domain
!MS2     use cover_module, ONLY: DCS_parts
!MS2     use set_blipgrid_module, ONLY : naba_atoms_of_blocks
!MS2     use angular_coeff_routines, ONLY : pao_elem_derivative_2
!MS2     use support_spec_format, ONLY: blips_on_atom, flag_one_to_one
!MS2     use functions_on_grid, ONLY: gridfunctions, fn_on_grid
!MS2 
!MS2     implicit none 
!MS2     integer,intent(in) :: myid, direction
!MS2 
!MS2     integer :: pao_fns
!MS2 
!MS2 
!MS2     !local
!MS2     real(double):: dcellx_block,dcelly_block,dcellz_block
!MS2     integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom, this_atom
!MS2     real(double):: xatom,yatom,zatom,alpha,step
!MS2     real(double):: xblock,yblock,zblock
!MS2     integer :: the_species
!MS2     integer :: j,iblock,the_l,ipoint, igrid
!MS2     real(double) :: r_from_i
!MS2     real(double) :: rr,a,b,c,d,x,y,z,nl_potential
!MS2     integer :: no_of_ib_ia, offset_position
!MS2     integer :: position,iatom
!MS2     integer :: stat, nl, npoint, ip, this_nsf
!MS2     integer :: i,m, m1min, m1max,acz,m1,l1,nsf1,atom_species, count1
!MS2     integer     , allocatable :: ip_store(:)
!MS2     real(double), allocatable :: x_store(:)
!MS2     real(double), allocatable :: y_store(:)
!MS2     real(double), allocatable :: z_store(:)
!MS2     real(double), allocatable :: r_store(:)
!MS2     real(double) :: coulomb_energy
!MS2     real(double) :: rcut
!MS2     real(double) :: r1, r2, r3, r4, core_charge, gauss_charge
!MS2     real(double) :: val, theta, phi, r_tmp
!MS2     real(double), dimension(:), allocatable :: sfsum
!MS2 
!MS2     call start_timer(tmr_std_allocation)
!MS2     allocate(ip_store(n_pts_in_block),x_store(n_pts_in_block),y_store(n_pts_in_block),z_store(n_pts_in_block), &
!MS2          r_store(n_pts_in_block), STAT=stat)
!MS2     if(stat /= 0) call cq_abort(' Error allocating store in PAO_to_grad_global: ',n_pts_in_block)
!MS2     call stop_timer(tmr_std_allocation)
!MS2     ! --  Start of subroutine  ---
!MS2 
!MS2     dcellx_block=rcellx/blocks%ngcellx
!MS2     dcelly_block=rcelly/blocks%ngcelly
!MS2     dcellz_block=rcellz/blocks%ngcellz
!MS2 
!MS2     call my_barrier()
!MS2 
!MS2     no_of_ib_ia = 0
!MS2     gridfunctions(pao_fns)%griddata = zero
!MS2 
!MS2     ! loop arround grid points in the domain, and for each
!MS2     ! point, get the core charge (local charge) density from the atoms 
!MS2     ! whose distances from the grid point are within the cutoff. 
!MS2 
!MS2     do iblock = 1, domain%groups_on_node ! primary set of blocks
!MS2        xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
!MS2        yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
!MS2        zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
!MS2        if(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) then ! if there are naba atoms
!MS2           iatom=0
!MS2           do ipart=1,naba_atoms_of_blocks(atomf)%no_of_part(iblock)
!MS2              jpart=naba_atoms_of_blocks(atomf)%list_part(ipart,iblock)
!MS2              if(jpart > DCS_parts%mx_gcover) then 
!MS2                 call cq_abort('PAO_to_grad_global: JPART ERROR ',ipart,jpart)
!MS2              endif
!MS2              ind_part=DCS_parts%lab_cell(jpart)
!MS2              do ia=1,naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart,iblock)
!MS2                 iatom=iatom+1
!MS2                 ii = naba_atoms_of_blocks(atomf)%list_atom(iatom,iblock)
!MS2                 icover= DCS_parts%icover_ibeg(jpart)+ii-1
!MS2                 ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)
!MS2 
!MS2                 if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
!MS2                    call cq_abort('PAO_to_grad_global: globID ERROR ', &
!MS2                         ii,parts%icell_beg(ind_part))
!MS2                 endif
!MS2                 if(icover > DCS_parts%mx_mcover) then
!MS2                    call cq_abort('PAO_to_grad_global: icover ERROR ', &
!MS2                         icover,DCS_parts%mx_mcover)
!MS2                 endif
!MS2 
!MS2                 xatom=DCS_parts%xcover(icover)
!MS2                 yatom=DCS_parts%ycover(icover)
!MS2                 zatom=DCS_parts%zcover(icover)
!MS2 
!MS2                 the_species=species_glob(ig_atom)
!MS2                 this_nsf = nsf_species(the_species)
!MS2                 call start_timer(tmr_std_allocation)
!MS2                 allocate(sfsum(this_nsf))
!MS2                 call stop_timer(tmr_std_allocation)
!MS2 
!MS2                 !calculates distances between the atom and integration grid points
!MS2                 !in the block and stores which integration grids are neighbours.
!MS2                 rcut = r_h + RD_ERR
!MS2                 call check_block (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
!MS2                      npoint,ip_store,r_store,x_store,y_store,z_store,n_pts_in_block) !out
!MS2                 r_from_i = sqrt((xatom-xblock)**2+(yatom-yblock)**2+ &
!MS2                      (zatom-zblock)**2 )
!MS2                 if(flag_one_to_one) then
!MS2                    this_atom = the_species
!MS2                 else
!MS2                    this_atom = ig_atom
!MS2                 end if
!MS2 
!MS2                 if(npoint > 0) then
!MS2                    !offset_position = (no_of_ib_ia-1) * nsf * n_pts_in_block
!MS2                    offset_position = no_of_ib_ia
!MS2                    do ip=1,npoint
!MS2                       ipoint=ip_store(ip)
!MS2                       position= offset_position + ipoint
!MS2                       if(position > gridfunctions(pao_fns)%size) call cq_abort &
!MS2                            ('PAO_to_grad_global: position error ', position, gridfunctions(pao_fns)%size)
!MS2 
!MS2                       r_from_i = r_store(ip)
!MS2                       x = x_store(ip)
!MS2                       y = y_store(ip)
!MS2                       z = z_store(ip)
!MS2                       sfsum = zero
!MS2                       ! For this point-atom offset, we loop over ALL PAOs, and accumulate each PAO multiplied by
!MS2                       ! the appropriate coefficient in a local store before putting it in place later
!MS2                       ! Actually, the zeta should be *inside* the l loop
!MS2                       count1 = 1
!MS2                       do l1 = 0,blips_on_atom(this_atom)%lmax
!MS2                          do acz = 1,blips_on_atom(this_atom)%naczs(l1)
!MS2                             do m1=-l1,l1
!MS2                                call pao_elem_derivative_2(direction,the_species,l1,acz,m1,x,y,z,val)
!MS2                                if(flag_one_to_one) then
!MS2                                   sfsum(count1) = sfsum(count1) + val
!MS2                                else
!MS2                                   do nsf1 = 1,blips_on_atom(this_atom)%nsuppfuncs
!MS2                                      sfsum(nsf1) = sfsum(nsf1) + &
!MS2                                           blips_on_atom(this_atom)%supp_func(nsf1)%coefficients(count1)* val
!MS2                                   end do ! nsf1
!MS2                                end if
!MS2                                count1 = count1+1
!MS2                             end do ! m1
!MS2                          end do ! acz
!MS2                       end do ! l1
!MS2                       ! Now store data
!MS2                       do nsf1 = 1,this_nsf
!MS2                          gridfunctions(pao_fns)%griddata(position+(nsf1-1)*n_pts_in_block) = sfsum(nsf1)
!MS2                       end do ! nsf
!MS2                    enddo ! ip=1,npoint
!MS2                 endif! (npoint > 0) then
!MS2                 call start_timer(tmr_std_allocation)
!MS2                 deallocate(sfsum)
!MS2                 call stop_timer(tmr_std_allocation)
!MS2                 no_of_ib_ia = no_of_ib_ia + this_nsf*n_pts_in_block
!MS2              enddo ! naba_atoms
!MS2           enddo ! naba_part
!MS2        endif !(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) !naba atoms?
!MS2     enddo ! iblock : primary set of blocks
!MS2     call start_timer(tmr_std_allocation)
!MS2     deallocate(ip_store,x_store,y_store,z_store, r_store, STAT=stat)
!MS2     if(stat /= 0) call cq_abort(' Error in deallocation in PAO_to_grad_global')
!MS2     call stop_timer(tmr_std_allocation)
!MS2     call my_barrier()
!MS2     return
!MS2   end subroutine PAO_to_grad_global
!MS2 !!***

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
    real(double):: xatom,yatom,zatom,alpha,step
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
    gridfunctions(pao_fns)%griddata = zero

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
             if(jpart > DCS_parts%mx_gcover) then 
                call cq_abort('single_PAO_to_grad: JPART ERROR ',ipart,jpart)
             endif
             ind_part=DCS_parts%lab_cell(jpart)
             do ia=1,naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart,iblock)
                iatom=iatom+1
                ii = naba_atoms_of_blocks(atomf)%list_atom(iatom,iblock)
                icover= DCS_parts%icover_ibeg(jpart)+ii-1
                ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)

                if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                   call cq_abort('single_PAO_to_grad: globID ERROR ', &
                        ii,parts%icell_beg(ind_part))
                endif
                if(icover > DCS_parts%mx_mcover) then
                   call cq_abort('single_PAO_to_grad: icover ERROR ', &
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

                if(npoint > 0) then
                   !offset_position = (no_of_ib_ia-1) * npao * n_pts_in_block
                   offset_position = no_of_ib_ia
                   do ip=1,npoint
                      ipoint=ip_store(ip)
                      position= offset_position + ipoint
                      if(position > gridfunctions(pao_fns)%size) call cq_abort &
                           ('single_pao_to_grad: position error ', position, gridfunctions(pao_fns)%size)

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
                               if(position+(count1-1)*n_pts_in_block > gridfunctions(pao_fns)%size) &
                                    call cq_abort('single_pao_to_grad: position error ', &
                                    position, gridfunctions(pao_fns)%size)
                               gridfunctions(pao_fns)%griddata(position+(count1-1)*n_pts_in_block) = val
                               count1 = count1+1
                            end do ! m1
                         end do ! acz
                      end do ! l1
                   enddo ! ip=1,npoint
                endif! (npoint > 0) then
                no_of_ib_ia = no_of_ib_ia + npao_species(the_species)*n_pts_in_block
             enddo ! naba_atoms
          enddo ! naba_part
       endif !(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) !naba atoms?
    enddo ! iblock : primary set of blocks
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
