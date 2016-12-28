! $Id$
! -----------------------------------------------------------
! Module S_matrix_module
! -----------------------------------------------------------
! Code area 3: Operators
! -----------------------------------------------------------

!!****h* Conquest/S_matrix_module
!!  NAME
!!   S_matrix_module
!!  PURPOSE
!!   Collects together the routines needed to get the S matrix
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/05/2001
!!  MODIFICATION HISTORY
!!   16:03, 04/02/2003 drb 
!!    Created get_onsite_S (mainly taken from kinetic energy)
!!   14:26, 26/02/2003 drb 
!!    Corrected get_S_matrix
!!   07:52, 2003/09/22 dave
!!    Added flags to choose between blips and PAOs as basis set
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2008/02/01 17:53 dave
!!    Changes for output to file not stdout
!!   2008/05/22 ast
!!    Added timers
!!   2013/03/06 L.Tong
!!   - Added module globals InvSMaxSteps and InvSDeltaOmegaTolerance
!!     which can be set by the user, instead of original hard-wired
!!     numbers of 100 and 0.0001
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2015/04/29 17:23 dave and kane
!!    Added get_r_on_support for calculations of stress, polarisation and TDDFT in periodic boundaries
!!   2016/08/05 15:00 nakata
!!    Renamed get_r_on_support -> get_r_on_atomfns
!!  SOURCE
!!
module S_matrix_module

  use datatypes
  use global_module,          only: io_lun
  use timer_module,           only: start_timer, stop_timer
  use timer_module,           only: start_backtrace, stop_backtrace
  use timer_stdclocks_module, only: tmr_std_smatrix

  implicit none

  real(double) :: InvSTolerance
  integer      :: InvSMaxSteps
  real(double) :: InvSDeltaOmegaTolerance

  ! RCS tag for object file identification
  character(len=80), save, private :: &
       RCSid = "$Id$"
!!***

contains

! -----------------------------------------------------------
! Subroutine get_S_matrix
! -----------------------------------------------------------

!!****f* S_matrix_module/get_S_matrix *
!!
!!  NAME 
!!   get_S_matrix
!!  USAGE
!! 
!!  PURPOSE
!!   Gets a new S matrix by doing blip_to_support, integrating
!!   and updating inverse S matrix
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/05/2001
!!  MODIFICATION HISTORY
!!   08:35, 2003/02/05 dave
!!    Added lines and call to get_onsite_S to replace onsite elements with analytical values
!!   14:25, 26/02/2003 drb
!!    Removed - breaking code for some reason
!!   07:53, 2003/09/22 dave
!!    Added choice between blips and PAOs
!!   08:18, 2003/10/01 dave
!!    A little more for PAOs
!!   12:52, 31/10/2003 drb 
!!    Changed call to assemble
!!   08:31, 2004/07/23 dave
!!    Added call to allow calculation of derivative of S wrt PAO coefficient
!!   2008/05/22 ast
!!    Added timer
!!   2009/10/27 07:10 dave
!!    Analytic on-site integrals added
!!   2011/07/21 11:48 dave
!!    Changes for cDFT (added call to Becke weight matrix)
!!   2012/04/03 11:57 dave
!!    Analytic integrals for S and KE added
!!   2013/03/06 L.Tong
!!   - Changed the hard-wired max number of Hottelling steps (100) to
!!     use module variable InvSMaxSteps
!!   - Changed the hard-wired deltaOmega tolerance for Hottelling
!!     iterations (0.0001_double) to use module variable
!!     InvSDeltaOmegaTolerance
!!   2015/06/08 lat
!!    - Added experimental backtrace
!!   2016/07/29 18:30 nakata
!!    Renamed supports_on_atom -> blips_on_atom
!!   2016/08/08 15:30 nakata
!!    Renamed supportfns -> atomfns
!!   2016/09/28 15:30 nakata
!!    Introduce get_S_matrix_blips and get_S_matrix_PAO
!!    Introduce PAO->SF transformation for contracted SFs
!!  TODO
!!    
!!  SOURCE
!!
  subroutine get_S_matrix(inode, ionode)

    use datatypes
    use numbers
    use global_module,               only: iprint_ops, flag_basis_set, &
                                           blips, PAOs,                &
                                           ni_in_cell,                 &
                                           flag_perform_cdft,          &
                                           IPRINT_TIME_THRES1,         &
                                           flag_do_SFtransform
    use matrix_data,                 only: Srange
    use mult_module,                 only: matS, matSatomf, transform_ATOMF_to_SF
    use io_module,                   only: dump_matrix
    use timer_module,                only: cq_timer, start_timer,      &
                                           stop_print_timer,           &
                                           WITH_LEVEL
    use cdft_module,                 only: make_weights

    implicit none

    ! Passed variables
    integer :: inode, ionode
    
    ! Local variables
    type(cq_timer) :: tmr_l_tmp1
    type(cq_timer) :: backtrace_timer


!****lat<$
    call start_backtrace(t=backtrace_timer,who='get_S_matrix',where=3,level=2)
!****lat>$

    call start_timer(tmr_std_smatrix)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)

    ! get S in atomic-function basis
    if (flag_basis_set == blips) call get_S_matrix_blips(inode, ionode)
    if (flag_basis_set == PAOs)  call get_S_matrix_PAO(inode, ionode)

    ! For blips and one_to_one PAOs, atomic functions are SFs so that no transformation is needed.
    ! For contracted SFs, transform S from atomic-function basis to SF basis
    if (flag_do_SFtransform) call transform_ATOMF_to_SF(matS, matSatomf, 1, Srange)

    call dump_matrix("NS1",matS,inode)

    ! get the new InvS matrix
    call  Iter_Hott_InvS(iprint_ops, InvSMaxSteps, &
                         InvSDeltaOmegaTolerance, ni_in_cell, &
                         inode, ionode)

    if (flag_perform_cdft) call make_weights

    call stop_print_timer(tmr_l_tmp1, "get_S_matrix", &
                          IPRINT_TIME_THRES1)
    call stop_timer(tmr_std_smatrix)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_S_matrix')
!****lat>$

    return
  end subroutine get_S_matrix
!!***

! -----------------------------------------------------------
! Subroutine get_S_matrix_PAO
! -----------------------------------------------------------

!!****f* S_matrix_module/get_S_matrix_PAO *
!!
!!  NAME 
!!   get_S_matrix_PAO
!!  USAGE
!! 
!!  PURPOSE
!!   Gets a new S matrix in PAO-basis by assemble_2
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler, A. Nakata
!!  CREATION DATE
!!   28/09/2016
!!  MODIFICATION HISTORY
!!   2016/09/28 nakata
!!   This subroutine is the part for PAOs in previous get_S_matrix
!!   
!!  TODO
!!    
!!  SOURCE
!!
  subroutine get_S_matrix_PAO(inode, ionode)

    use global_module,               only: iprint_ops
    use matrix_data,                 only: aSa_range
    use mult_module,                 only: matSatomf
    use build_PAO_matrices,          only: assemble_2
    use io_module,                   only: dump_matrix

    implicit none

    ! Passed variables
    integer :: inode, ionode


    ! Get S matrix with assemble
    if (inode == ionode .and. iprint_ops > 2) &
         write (io_lun, *) 'Calling assemble_2 for Satomf: ', matSatomf
    call assemble_2(aSa_range, matSatomf, 1)
    !call dump_matrix("NS_atomf",matSatomf,inode)

    return
  end subroutine get_S_matrix_PAO
!!***

! -----------------------------------------------------------
! Subroutine get_S_matrix_blips
! -----------------------------------------------------------

!!****f* S_matrix_module/get_S_matrix_blips *
!!
!!  NAME 
!!   get_S_matrix_blips
!!  USAGE
!! 
!!  PURPOSE
!!   Gets a new S matrix by doing blip_to_support and integrating
!!
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler, A. Nakata
!!  CREATION DATE
!!   28/09/2016
!!  MODIFICATION HISTORY
!!   2016/09/28 nakata
!!   This subroutine is the part for blips in previous get_S_matrix
!!
!!  TODO
!!    
!!  SOURCE
!!
  subroutine get_S_matrix_blips(inode, ionode)

    use datatypes
    use numbers
    use global_module,               only: iprint_ops,            &
                                           flag_onsite_blip_ana,  &
                                           id_glob, species_glob, &
                                           flag_analytic_blip_int, nspin
    use matrix_data,                 only: Srange, halo, blip_trans
    use mult_module,                 only: matS, ltrans, matrix_scale, &
                                           matK, matM12, return_matrix_block_pos,&
                                           matKE, matrix_pos
    use set_bucket_module,           only: rem_bucket
    use calc_matrix_elements_module, only: get_matrix_elements_new
    use blip_grid_transform_module,  only: blip_to_support_new
    use primary_module ,             only: bundle
    use functions_on_grid,           only: atomfns
    use support_spec_format,         only: blips_on_atom, support_function, &
                                           coefficient_array, support_gradient, &
                                           grad_coeff_array
    use species_module,              only: nsf_species
    use comms_module,                ONLY: start_blip_transfer, fetch_blips
    use group_module,                ONLY: parts
    use blip,                        ONLY: blip_info
    use cover_module,                ONLY: BCS_parts
    use GenComms,                    ONLY: myid, cq_abort, my_barrier, mtime
    use mpi

    implicit none

    ! Passed variables
    integer :: inode, ionode
    
    ! Local variables
    integer :: np, ni, iprim, spec, this_nsf, icall, jpart, ind_part, j_in_halo, pb_len, ist, nab
    integer :: i, j,jseq,specj,i_in_prim,speci, pb_st, nblipsj, gcspart, nod, nsfj,sends,ierr, spin
    integer :: neigh_global_num, neigh_global_part, neigh_species, neigh_prim, wheremat, this_nsfi, this_nsfj
    integer, allocatable, dimension(:) :: nreqs
    real(double) :: dx, dy, dz, time0, time1
    real(double), allocatable, dimension(:), target :: part_blips
    type(support_function) :: supp_on_j
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    real(double), allocatable, dimension(:,:,:) :: this_data_K
    real(double), allocatable, dimension(:,:,:) :: this_data_M12


    ! Project support functions onto grid
    if (inode == ionode .and. iprint_ops > 2) &
         write (io_lun, *) 'Doing blip-to-support ', atomfns
    call blip_to_support_new(inode-1, atomfns)

    if (inode == ionode .and. iprint_ops > 2) &
         write (io_lun, *) 'Doing integration ', atomfns
    ! Integrate
    if(flag_analytic_blip_int) then
       call matrix_scale(zero,matS)
       call matrix_scale(zero,matKE)
       grad_coeff_array = zero
       allocate(nreqs(blip_trans%npart_send))
       ! For speed, we should have a blip_trans%max_len and max_nsf and allocate once
       time0 = mtime()
       call start_blip_transfer(nreqs,sends,parts%mx_ngonn)
       time1 = mtime()
       do jpart=1,halo(Srange)%np_in_halo
          pb_len = blip_trans%len_recv(jpart)
          ind_part = halo(Srange)%lab_hcell(jpart)
          nod = parts%i_cc2node(ind_part)
          ! Fetch remote blip coefficients for partition
          ! or copy local blip coefficients
          time0 = mtime()
          if(jpart>1) then
             if(ind_part/=halo(Srange)%lab_hcell(jpart-1)) then
                allocate(part_blips(pb_len))
                part_blips = zero
                if(nod==myid+1) then
                   pb_st = blip_trans%partst(parts%i_cc2seq(ind_part))
                   part_blips(1:pb_len) = coefficient_array(pb_st:pb_st+pb_len-1)
                else
                   call fetch_blips(part_blips,pb_len,nod-1,(myid)*parts%mx_ngonn + parts%i_cc2seq(ind_part))
                end if
             end if
          else
             allocate(part_blips(pb_len))
             part_blips = zero
             if(nod==myid+1) then
                pb_st = blip_trans%partst(parts%i_cc2seq(ind_part))
                part_blips(1:pb_len) = coefficient_array(pb_st:pb_st+pb_len-1)
             else
                call fetch_blips(part_blips,pb_len,nod-1,(myid)*parts%mx_ngonn + parts%i_cc2seq(ind_part))
             end if
          endif
          time1 = mtime()
          gcspart = halo(Srange)%i_hbeg(halo(Srange)%lab_hcover(jpart))
          pb_st = 1
          time0 = mtime()
          do j=1,halo(Srange)%nh_part(jpart) ! Loop over atoms j in partition
             j_in_halo = halo(Srange)%j_beg(jpart)+j-1
             jseq = halo(Srange)%j_seq(j_in_halo)

             specj = species_glob( id_glob( parts%icell_beg(halo(Srange)%lab_hcell(jpart))+jseq-1) )
             nblipsj = blip_info(specj)%NBlipsRegion
             this_nsfj = nsf_species(specj)
             allocate(supp_on_j%supp_func(this_nsfj))
             do nsfj=1,this_nsfj
                supp_on_j%supp_func(nsfj)%ncoeffs = nblipsj
                supp_on_j%supp_func(nsfj)%coefficients => part_blips(pb_st:pb_st+nblipsj-1)
                pb_st = pb_st+nblipsj
             end do
             do i=1,ltrans(Srange)%n_hnab(j_in_halo) ! Loop over atoms i: primary set neighbours of j
                i_in_prim=ltrans(Srange)%i_prim(ltrans(Srange)%i_beg(j_in_halo)+i-1)
                speci = bundle%species(i_in_prim)
                this_nsfi = nsf_species(speci)
                dx = BCS_parts%xcover(gcspart+jseq-1)-bundle%xprim(i_in_prim)
                dy = BCS_parts%ycover(gcspart+jseq-1)-bundle%yprim(i_in_prim)
                dz = BCS_parts%zcover(gcspart+jseq-1)-bundle%zprim(i_in_prim)
                allocate(this_data_K(this_nsfi,this_nsfj,nspin),this_data_M12(this_nsfi,this_nsfj,nspin))
                this_data_K = zero
                this_data_M12 = zero
                do spin=1,nspin
                   wheremat = matrix_pos(matK(spin),i_in_prim,j_in_halo,1,1)
                   call return_matrix_block_pos(matK(spin),wheremat,this_data_K(:,:,spin),this_nsfi*this_nsfj)
                   wheremat = matrix_pos(matM12(spin),i_in_prim,j_in_halo,1,1)
                   call return_matrix_block_pos(matM12(spin),wheremat,this_data_M12(:,:,spin),this_nsfi*this_nsfj)
                end do
                call get_S_analytic(blips_on_atom(i_in_prim),supp_on_j,support_gradient(i_in_prim), &
                     matS,matKE,this_data_M12,this_data_K, i_in_prim, &
                     j_in_halo,dx,dy,dz,speci,specj,this_nsfi,this_nsfj)
                deallocate(this_data_M12,this_data_K)
             end do
             do nsfj=1,this_nsfj
                nullify(supp_on_j%supp_func(nsfj)%coefficients)
             end do
             deallocate(supp_on_j%supp_func)
          end do
          if(jpart<halo(Srange)%np_in_halo) then
             if(ind_part/=halo(Srange)%lab_hcell(jpart+1)) then
                deallocate(part_blips)
             end if
          else
             deallocate(part_blips)
          end if
          time1 = mtime()
       end do
       if(sends>0) then
          do i=1,sends
             call MPI_Wait(nreqs(i),mpi_stat,ierr)
             if(ierr/=0) call cq_abort("Error waiting for blip send to finish",i)
          end do
       end if
       call my_barrier
       deallocate(nreqs)
    else
       call get_matrix_elements_new(inode-1, rem_bucket(1), matS, &
            atomfns, atomfns)
       ! Do the onsite elements analytically
       if (flag_onsite_blip_ana) then
          iprim = 0
          do np = 1, bundle%groups_on_node
             if (bundle%nm_nodgroup(np) > 0) then
                do ni = 1, bundle%nm_nodgroup(np)
                   iprim = iprim + 1
                   spec = bundle%species(iprim)
                   this_nsf = nsf_species(spec)
                   call get_onsite_S(blips_on_atom(iprim), matS, &
                        np, ni, iprim, this_nsf, spec)
                end do
             end if
          end do
       end if
    end if  ! flag_analytic_blip_int

    return
  end subroutine get_S_matrix_blips
!!***
    
! -----------------------------------------------------------
! Subroutine get_r_on_atomfns
! -----------------------------------------------------------

!!****f* S_matrix_module/get_r_on_atomfns *
!!
!!  NAME 
!!   get_r_on_atomfns
!!  USAGE
!!   
!!  PURPOSE
!!   Scales the support functions on the grid by x/y/z
!!  INPUTS
!!   direction chooses x/y/z (1-3) or all directions (0)
!!   inputgf is the grid function to be scaled by x/y/z  
!!   gridfunc1 should contain support functions on entry;
!!   on exit will be scaled by x/y/z
!!   If direction=0 we do all three directions, and gridfunc1-3 should all contain support functions
!!   We can choose memory vs storage
!!  USES
!!   
!!  AUTHOR
!!   D. R. Bowler with J. Kane Shenton
!!  CREATION DATE
!!   2015/04/29
!!  MODIFICATION HISTORY
!!   2015/05/01 09:01 dave
!!    Changed position-finding routines to use check-block (removes some problems with z-component)  
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks
!!   2016/08/01 17:30 nakata
!!    Introduced atomf instead of sf
!!   2016/08/05 15:00 nakata
!!    Renamed get_r_on_support -> get_r_on_atomfns
!!   2016/11/09 21:00 nakata
!!    Used natomf_species instead of nsf_species
!!  SOURCE
!!  
  subroutine get_r_on_atomfns(direction,inputgf,gridfunc1,gridfunc2,gridfunc3)

    use datatypes
    use numbers
    use global_module,               only: rcellx,rcelly,rcellz, atomf, species_glob, ni_in_cell, id_glob
    use cover_module,                only: DCS_parts
    use block_module,                only: nx_in_block,ny_in_block,nz_in_block, &
                                           n_blocks, n_pts_in_block
    use group_module,                only: blocks, parts
    use primary_module,              only: domain
    use set_blipgrid_module,         only: naba_atoms_of_blocks
    use functions_on_grid,           only: gridfunctions, fn_on_grid
    use dimens,                      only: n_my_grid_points, r_h
    use GenComms,                    only: cq_abort
    use species_module,              only: natomf_species
    use PAO_grid_transform_module,   only: check_block

    implicit none

    ! Passed variables
    integer :: direction
    integer :: inputgf,gridfunc1
    integer, OPTIONAL :: gridfunc2, gridfunc3

    ! Local variables
    integer        :: ipoint, iblock, nsf1, ix, iy, iz, igrid, position
    integer        :: iatom, this_nsf
    integer        :: ipart, jpart, ind_part, ia, ii, icover, ig_atom, no_of_ib_ia
    real(double)   :: rx, ry, rz, dx, dy, dz
    real(double)   :: dcellx_block, dcelly_block, dcellz_block
    real(double)   :: dcellx_grid, dcelly_grid, dcellz_grid
    real(double)   :: xatom, yatom, zatom
    real(double)   :: xblock,yblock,zblock
    real(double)   :: sfni

    integer     , allocatable :: ip_store(:)
    real(double), allocatable :: x_store(:)
    real(double), allocatable :: y_store(:)
    real(double), allocatable :: z_store(:)
    real(double), allocatable :: r_store(:)
    integer :: offset_position, ip, npoint, stat
    real(double) :: x, y, z, rcut
    
    ! Test for direction=0 and all three grid functions
    if(direction==0.AND.(.NOT.PRESENT(gridfunc2).OR..NOT.PRESENT(gridfunc3))) then
       call cq_abort("Error in get_r_on_atomfns: three directions requested but not passed")
    end if
    allocate(ip_store(n_pts_in_block),x_store(n_pts_in_block),y_store(n_pts_in_block),z_store(n_pts_in_block), &
         r_store(n_pts_in_block), STAT=stat)
    if(stat /= 0) call cq_abort(' Error allocating store in get_r_on_atomfns: ',n_pts_in_block)

    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nx_in_block
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/ny_in_block
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nz_in_block
    no_of_ib_ia = 0
    do iblock = 1, domain%groups_on_node
       if (iblock*n_pts_in_block > n_my_grid_points) call cq_abort('get_nonSC_force: igrid error ', igrid,n_my_grid_points)
       xblock = (domain%idisp_primx(iblock) + domain%nx_origin - 1) * dcellx_block
       yblock = (domain%idisp_primy(iblock) + domain%ny_origin - 1) * dcelly_block
       zblock = (domain%idisp_primz(iblock) + domain%nz_origin - 1) * dcellz_block
       if (naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom = 0
          do ipart = 1, naba_atoms_of_blocks(atomf)%no_of_part(iblock)
             jpart = naba_atoms_of_blocks(atomf)%list_part(ipart,iblock)
             if (jpart > DCS_parts%mx_gcover) call cq_abort('set_ps: JPART ERROR ', ipart, jpart)
             ind_part = DCS_parts%lab_cell(jpart)
             do ia = 1, naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart, iblock)
                iatom = iatom + 1
                ii = naba_atoms_of_blocks(atomf)%list_atom(iatom, iblock)
                icover = DCS_parts%icover_ibeg(jpart) + ii - 1
                if (parts%icell_beg(ind_part) + ii - 1 > ni_in_cell) &
                     call cq_abort('set_ps: globID ERROR ', ii, parts%icell_beg(ind_part))
                ig_atom = id_glob(parts%icell_beg(ind_part) + ii - 1)
                if (icover > DCS_parts%mx_mcover) &
                     call cq_abort ('set_ps: icover ERROR ', icover, DCS_parts%mx_mcover)
                xatom = DCS_parts%xcover(icover)
                yatom = DCS_parts%ycover(icover)
                zatom = DCS_parts%zcover(icover)
                this_nsf = natomf_species(species_glob(ig_atom))
                rcut = r_h + RD_ERR
                call check_block (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
                     npoint,ip_store,r_store,x_store,y_store,z_store,n_pts_in_block) !out
                if(npoint > 0) then
                   !offset_position = (no_of_ib_ia-1) * nsf * n_pts_in_block
                   offset_position = no_of_ib_ia
                   do ip=1,npoint
                      ipoint=ip_store(ip)
                      position= offset_position + ipoint
                      if(position > gridfunctions(inputgf)%size) call cq_abort &
                           ('PAO_to_grad_global: position error ', position, gridfunctions(inputgf)%size)

                      !r_from_i = r_store(ip)
                      x = x_store(ip)
                      y = y_store(ip)
                      z = z_store(ip)
                      if(direction==0) then
                         rx = x
                         ry = y
                         rz = z
                         do nsf1=1,this_nsf
                            sfni = gridfunctions(inputgf)%griddata(position+(nsf1-1)*n_pts_in_block)
                            gridfunctions(gridfunc1)%griddata(position+(nsf1-1)*n_pts_in_block) = sfni * rx
                            gridfunctions(gridfunc2)%griddata(position+(nsf1-1)*n_pts_in_block) = sfni * ry
                            gridfunctions(gridfunc3)%griddata(position+(nsf1-1)*n_pts_in_block) = sfni * rz
                         enddo ! nsf1
                      else ! Store position in rx and scale appropriate direction
                         if(direction==1) then
                            rx = x
                         else if(direction==2) then
                            rx = y
                         else if(direction==3) then
                            rx = z
                         end if
                         do nsf1=1,this_nsf
                            sfni = gridfunctions(inputgf)%griddata(position+(nsf1-1)*n_pts_in_block)
                            gridfunctions(gridfunc1)%griddata(position+(nsf1-1)*n_pts_in_block) = sfni * rx
                         enddo ! nsf1
                      end if
                   end do
                end if
                no_of_ib_ia=no_of_ib_ia+this_nsf*n_pts_in_block
             end do ! naba_atoms
          end do ! naba_part
       end if !(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) !naba atoms?
    end do ! iblock : primary set of blocks
    deallocate(ip_store,x_store,y_store,z_store, r_store, STAT=stat)
    if(stat /= 0) call cq_abort(' Error in deallocation in PAO_to_grad_global')
    return
  end subroutine get_r_on_atomfns
!!***

! -----------------------------------------------------------
! Subroutine Iter_Hott_InvS
! -----------------------------------------------------------

!!****f* S_matrix_module/Iter_Hott_InvS *
!!
!!  NAME 
!!   Iter_Hott_InvS
!!  USAGE
!! 
!!  PURPOSE
!!   This subroutine finds the inverse of S by Hotelling's
!!   method - see Numerical Recipes (ed 2), pp49-50.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   16/10/98 DRB
!!  MODIFICATION HISTORY
!!   04/05/00 by D.R.Bowler to use new matrix mults
!!   22/05/2001 dave
!!    Added ROBODoc header, indented
!!   25/05/2001 dave
!!    Included in S_matrix_module, relocated HotInvSmm
!!   11/06/2001 dave
!!    Changed to use GenComms
!!   14:54, 28/08/2003 drb 
!!    Added a logical test for diagonalisation (which doesn't need S^-1)
!!   2004/10/28 drb
!!    Added Told so that we can revert to the best T at exit, and changed omega definition 
!!    so that it's divided by total number of orbitals; added user-set criterion on omega 
!!    tolerance
!!   2013/08/20 M.Arita
!!    Added variables and calls for reusing T-matrix
!!  SOURCE
!!
  subroutine Iter_Hott_InvS(output_level, n_L_iterations, tolerance,n_atoms,&
       inode, ionode)

    use datatypes
    use numbers
    use global_module, ONLY: IPRINT_TIME_THRES1,flag_MDold,flag_TmatrixReuse,restart_T, &
                             runtype,n_proc_old
    use matrix_data, ONLY: Trange, TSrange, mat, Srange
    use mult_module, ONLY: allocate_temp_matrix, free_temp_matrix, store_matrix_value, matrix_scale, matrix_sum, &
         matT, matS, return_matrix_value, T_trans
    use primary_module, only: bundle
    use GenComms, ONLY: gsum, my_barrier
    use DiagModule, ONLY: diagon
    use species_module, ONLY: nsf_species, species
    use timer_module, ONLY: cq_timer,start_timer,stop_print_timer,WITH_LEVEL
    use input_module, ONLY: leqi
    use io_module2, ONLY: grab_matrix2,dump_matrix2,InfoT
    use UpdateInfo_module, ONLY: Matrix_CommRebuild

    implicit none

    ! Passed variables
    integer :: output_level, inode, ionode,n_atoms, n_L_iterations
    real(double) ::  tolerance

    ! Local variables
    integer :: n_iterations, nn, np, nb, ist, ip, i,j,nsf1,nsf2
    integer :: matI, matT1, matTM, matTold, nfile, symm
    real(double) :: step, tot, eps, x, omega, oldomega, deltaomega, n_orbs
    type(cq_timer) :: tmr_l_tmp1
    logical,save :: flag_readT = .false.

    matI = allocate_temp_matrix(TSrange,0)
    matT1 = allocate_temp_matrix(Trange,0)
    matTold = allocate_temp_matrix(Trange,0)
    matTM = allocate_temp_matrix(TSrange,0)

    if(diagon) then ! Don't need S^-1 for diagonalisation
       call matrix_scale(zero,matT)
       ip = 1
       nb = 1
       do np = 1,bundle%groups_on_node
          if(bundle%nm_nodgroup(np)>0) then
             do i=1,bundle%nm_nodgroup(np)
                do j = 1, nsf_species(bundle%species(ip))
                   call store_matrix_value(matT,np,i,ip,nb,j,j,one,1)
                enddo
                ip = ip+1
             enddo
          end if ! bundle%nm_nodgroup(np)>0
       enddo
    else
       call start_timer(tmr_l_tmp1,WITH_LEVEL)
       n_orbs = zero
       do i=1,n_atoms
          n_orbs = n_orbs + real(nsf_species(species(i)),double)
       end do
       ! First construct the identity
       if (inode.eq.ionode.and.output_level>=2) write(io_lun,*) 'Zeroing data'
       call matrix_scale(zero,matT1)
       call matrix_scale(zero,matI)
       call matrix_scale(zero,matT)
       call matrix_scale(zero,matTold)
       call matrix_scale(zero,matTM)
       if (inode.eq.ionode.and.output_level>=2) write(io_lun,*) 'Creating I'
       ip = 1
       nb = 1
       do np = 1,bundle%groups_on_node
          if(bundle%nm_nodgroup(np)>0) then
             do i=1,bundle%nm_nodgroup(np)
                do j = 1, nsf_species(bundle%species(ip))
                   call store_matrix_value(matT,np,i,ip,nb,j,j,one,1)
                   call store_matrix_value(matI,np,i,ip,nb,j,j,one,1)
                enddo
                ip = ip+1
             enddo
          end if ! bundle%nm_nodgroup(np)>0
       enddo
       call my_barrier


       if (.NOT. flag_readT .AND. .NOT. restart_T) then ! Initialising invS
         ! Construct the initial guess for T as epsilon.S^T, where epsilon is
         ! given as 1/(sum_jk S^2_jk)
         tot = 0.0_double
         ip = 1
         do np = 1,bundle%groups_on_node
            if(bundle%nm_nodgroup(np)>0) then
               do i=1,bundle%nm_nodgroup(np)
                  do nb=1,mat(np,Srange)%n_nab(i)
                     ist = mat(np,Srange)%i_acc(i)+nb-1
                     do nsf1 = 1,mat(np,Srange)%ndimi(i)
                        do nsf2 = 1,mat(np,Srange)%ndimj(ist)
                           eps = return_matrix_value(matS,np,i,ip,nb,nsf1,nsf2)
                           tot = tot + eps*eps
                        enddo
                     enddo
                  enddo
                  ip = ip+1
               enddo
            end if ! bundle%nm_nodgroup(np)>0
         enddo
         call gsum(tot)
         eps = 1.0_double/(tot)
         if(output_level>1.and.inode==ionode) write(io_lun,*) 'Eps, tot: ',eps,tot
         call matrix_scale(zero,matT)
         call matrix_sum(zero,matT,eps,matS)
         call stop_print_timer(tmr_l_tmp1,"inverse S preliminaries",IPRINT_TIME_THRES1)
       elseif (.NOT. flag_MDold .AND. (flag_readT .OR. restart_T)) then
         !Reusing previously computed inverse S-matrix
         call grab_matrix2('T',inode,nfile,InfoT)
         call my_barrier()
         call Matrix_CommRebuild(InfoT,Trange,T_trans,matT,nfile,symm)
       endif
       ! and evaluate the current value of the functional and its gradient
       deltaomega = zero
       oldomega = zero
       if (inode==ionode.and.output_level>=2) write(io_lun,*) 'Starting loop'
       do n_iterations=1,n_L_iterations
          call start_timer(tmr_l_tmp1,WITH_LEVEL)
          if (inode==ionode.and.output_level>=2) &
               write(io_lun,2) n_iterations
          deltaomega = deltaomega * half
          ! check for convergence
          if(n_iterations<3.or.abs(deltaomega)>tolerance) then
             call HotInvS_mm( matI, matS, matT, matT1, matTM, omega,n_iterations)
             deltaomega = omega - oldomega
             if(inode==ionode.and.output_level>=1) then
                write(io_lun,*) 'Omega is ',omega/n_orbs
                if(omega>zero) write(io_lun,*) 'R is ',sqrt(omega)/n_orbs
                write(io_lun,*) 'deltaomega is ',n_iterations,deltaomega
             endif
             if ( omega>oldomega.and.oldomega/=0.0_double) then
                if(inode==ionode) write(io_lun,*) 'Truncation error reached !'
                call matrix_sum(zero,matT,one,matTold)
                call stop_print_timer(tmr_l_tmp1,"an inverse S iteration",IPRINT_TIME_THRES1)
                exit
             endif
             oldomega = omega 
             call matrix_sum(zero,matTold,one,matT)
             call matrix_sum(zero,matT,one,matT1)
          else
             call stop_print_timer(tmr_l_tmp1,"an inverse S iteration",IPRINT_TIME_THRES1)
             exit  ! Leave the do loop
          endif
          call stop_print_timer(tmr_l_tmp1,"an inverse S iteration",IPRINT_TIME_THRES1)
       end do
       ! If this isn't a good guess, then reset to I
       if((omega/n_orbs)>InvSTolerance) then
          if(inode==ionode) write(io_lun,*) 'Setting InvS to I'
          call matrix_scale(zero,matT)
          ip = 1
          nb = 1
          do np = 1,bundle%groups_on_node
             if(bundle%nm_nodgroup(np)>0) then
                do i=1,bundle%nm_nodgroup(np)
                   do j = 1, nsf_species(bundle%species(ip))
                      call store_matrix_value(matT,np,i,ip,nb,j,j,one,1)
                   enddo
                   ip = ip+1
                enddo
             end if ! bundle%nm_nodgroup(np)>0
          enddo
       endif
    end if! End if NOT diagon
    call free_temp_matrix(matTM)
    call free_temp_matrix(matTold)
    call free_temp_matrix(matT1)
    call free_temp_matrix(matI)
    ! Dump T-matrix
    if (.NOT. flag_MDold .AND. .NOT. leqi(runtype,'static') .AND. flag_TmatrixReuse) then
      call dump_matrix2('T',matT,inode,Trange)
      flag_readT = .true.
    endif

1   format(20x,'Starting functional value: ',f15.7,' a.u.')
2   format(/,20x,'Conjugate Gradients InvS iteration:',i5)
3   format(/,20x,'Functional value reached after ',i5,' InvS iterations: ',&
         /,20x,' Omega: ', f15.7, ' DeltaOmega: ', f15.7)
4   format('InvS is ',4i5,f15.7)
5   format('T0S,A is ',4i5,2f15.7)
    return

  end subroutine Iter_Hott_InvS
!!***

!!****f* S_matrix_module/HotInvS_mm *
!!
!!  NAME 
!! 
!!  USAGE
!! 
!!  PURPOSE
!!   To perform all the required matrix multiplications to produce
!!   an inverse S matrix, which will give "tensorially correct" 
!!   gradients for the density matrix
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   15/10/98
!!  MODIFICATION HISTORY
!!   11/06/2001 dave
!!    Changed to use GenComms
!!   12:59, 16/11/2004 dave & tsuyoshi
!!    Removed TS_TS_T multiplication and replaced omega with Frobenius norm (basically same !)
!!   12:12, 03/02/2006 drb 
!   Reworked for new scheme
!!  SOURCE
!!
  subroutine HotInvS_mm( matI, matS, matT0, matT1, matA, omega,n)

    use datatypes
    use numbers
    use mult_module, ONLY: matrix_scale, matrix_product, matrix_sum, free_temp_matrix, allocate_temp_matrix, &
         mult, T_S_TS, TS_T_T, matrix_product_trace
    use matrix_data, ONLY: TSrange, Trange
    use GenComms, ONLY: gsum, inode, my_barrier

    implicit none

    ! Passed variables
    integer :: matI, matS, matT0, matT1, matA, n
    real( double ) :: omega

    !     Local Variables
    integer :: matT0S, matGrad

    matT0S = allocate_temp_matrix(TSrange,0)
    matGrad = allocate_temp_matrix(Trange,0)

    call matrix_scale(zero,matT0S)
    call matrix_scale(zero,matGrad)
    ! Create T0.S
    call matrix_product(matT0, matS, matT0S, mult( T_S_TS ) )
    ! Now A=I-TS, as a diagnostic
    call my_barrier()
    call matrix_sum(zero,matA,one,matI)
    call matrix_sum(one,matA,-one,matT0S)
    ! Create T0S.T0
    call matrix_product(matT0S, matT0, matGrad, mult( TS_T_T ) )
    ! T1 = 2T0 - Grad
    call matrix_sum(zero,matT1,two,matT0)
    call matrix_sum(one,matT1,-one,matGrad)
    omega = matrix_product_trace(matA,matA)
    call free_temp_matrix(matGrad)
    call free_temp_matrix(matT0S)

  end subroutine HotInvS_mm
!!***

! -----------------------------------------------------------
! Subroutine get_onsite_S
! -----------------------------------------------------------

!!****f* S_matrix_module/get_onsite_S *
!!
!!  NAME 
!!   get_onsite_S
!!  USAGE
!! 
!!  PURPOSE
!!   This routine evalutes the overlay matrix elements for the block
!!   diagonal. This can be done analytically relatively quickly, because the
!!   support functions are represented on the _same_ blip grid, as they are
!!   associated with a single atom. 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler/C.M.Goringe
!!  CREATION DATE
!!   16:03, 04/02/2003 drb 
!!  MODIFICATION HISTORY
!!   2009/10/27 07:10 dave
!!    Correctly coded
!!   2011/11/15 08:03 dave
!!    Changes to blip data
!!  SOURCE
!!
  subroutine get_onsite_S(blip_co, matS, np, nn, ip, this_nsf, spec)

    use datatypes
    use numbers
    use GenBlas, ONLY: axpy, copy, scal, gemm
    use blip, ONLY: blip_info
    use support_spec_format, ONLY: support_function
    use GenComms, ONLY: cq_abort
    use mult_module, ONLY: store_matrix_value, scale_matrix_value, return_matrix_value

    implicit none

    ! Shared Variables
    integer :: this_nsf, spec, np, nn, ip, matS

    type(support_function) :: blip_co

    ! Local Variables
    integer, parameter :: MAX_D = 3

    real(double) ::  FAC(0:MAX_D)

    real(double), allocatable, dimension(:) :: work1, work2, work4, work6
    real(double) :: temp(this_nsf,this_nsf)

    integer :: dx, dy, dz, offset, l, at, nsf1, stat, i1, i2

    allocate(work1(blip_info(spec)%FullArraySize*this_nsf), &
         work2(blip_info(spec)%FullArraySize*this_nsf), &
         work4(blip_info(spec)%FullArraySize*this_nsf), &
         work6(blip_info(spec)%FullArraySize*this_nsf), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating arrays for onsite S elements: ", &
         blip_info(spec)%FullArraySize,this_nsf)
    ! first, we copy the blip functions for this atom onto a cubic grid;
    ! we make this grid 'too big' in order to have a fast routine below.

    FAC(0) = 151.0_double/140.0_double
    FAC(1) = 1191.0_double/2240.0_double
    FAC(2) = 3.0_double/56.0_double
    FAC(3) = 1.0_double/2240.0_double

    work1 = zero
    offset = blip_info(spec)%BlipArraySize+1

    do dx = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
       do dy = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
          do dz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             l = blip_info(spec)%blip_number(dx,dy,dz)
             if (l.ne.0) then
                at = (((dz+offset)*blip_info(spec)%OneArraySize + (dy+offset))*blip_info(spec)%OneArraySize + &
                     (dx+offset)) * this_nsf
                do nsf1 = 1,this_nsf
                   work1(nsf1+at) = blip_co%supp_func(nsf1)%coefficients(l)
                enddo
             end if
          end do
       end do
    end do

    ! now, for each direction in turn, we need to apply a number of
    ! 'spreading operations'. Do z first... put blip(z) in 2, and
    ! del2blip(z) in 3

    call copy(blip_info(spec)%FullArraySize*this_nsf,work1,1,work2,1)
    call scal(blip_info(spec)%FullArraySize*this_nsf,FAC(0),work2,1)
    do dz = 1, MAX_D
       offset = dz * blip_info(spec)%OneArraySize * blip_info(spec)%OneArraySize * this_nsf
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dz), &
            work1(1:), 1, work2(1+offset:), 1 )
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dz), &
            work1(1+offset:), 1, work2(1:), 1 )
    end do

    ! now do y : put blip(y).blip(z) in 4,
    ! blip(y).del2blip(z) + del2blip(y).blip(z)  in 5

    call copy(blip_info(spec)%FullArraySize*this_nsf,work2,1,work4,1)
    call scal(blip_info(spec)%FullArraySize*this_nsf,FAC(0),work4,1)
    do dy = 1, MAX_D
       offset = dy * blip_info(spec)%OneArraySize * this_nsf
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dy), &
            work2(1:), 1, work4(1+offset:), 1 )
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dy), &
            work2(1+offset:), 1, work4(1:), 1 )
    end do

    work6 = zero
    call axpy(blip_info(spec)%FullArraySize*this_nsf,FAC(0),work4,1,work6,1)
    do dx = 1, MAX_D
       offset = dx * this_nsf
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dx), &
            work4(1:), 1, work6(1+offset:), 1 )
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dx), &
            work4(1+offset:), 1, work6(1:), 1 )
    end do

    ! and now get the matrix elements by multiplication...

    temp = zero
    call gemm('n','t',this_nsf,this_nsf,&
         blip_info(spec)%OneArraySize*blip_info(spec)%OneArraySize*blip_info(spec)%OneArraySize, &
         one,work1,this_nsf,work6,this_nsf,zero,temp,this_nsf )
    do i1 = 1,this_nsf
       do i2=1,this_nsf
          call scale_matrix_value(matS,np,nn,ip,0,i2,i1,zero,1)
          call store_matrix_value(matS,np,nn,ip,0,i2,i1, &
               blip_info(spec)%SupportGridSpacing*blip_info(spec)%SupportGridSpacing* &
               blip_info(spec)%SupportGridSpacing*temp(i2,i1),1)
       end do
    end do
    deallocate(work1,work2, work4,work6, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating arrays for onsite S blip elements: ",&
         blip_info(spec)%FullArraySize,this_nsf)
    return
  end subroutine get_onsite_S
!!***

!!****f* S_matrix_module/get_S_analytic *
!!
!!  NAME 
!!   get_S_analytic
!!  USAGE
!!   
!!  PURPOSE
!!   Calculates various matrix elements and gradients analytically from blips.  The overlap
!!   and kinetic energy matrix elements and the associated contributions to the gradient of 
!!   energy wrt blip coefficients are found (if integrals are performed numerically, these 
!!   are found in different places: S is found in this module, KE in H_matrix_module, and the
!!   gradients in blip_gradient.module and specifically the S pulay term (in get_support_gradient)
!!   and the KE term.
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2012/02 (roughly)
!!  MODIFICATION HISTORY
!!  
!!  SOURCE
!!  
  subroutine get_S_analytic(blipL_co, blipR_co, blip_grad, matS, matT, dataM12, dataK, ip, j_in_halo, &
       dx, dy, dz, speci, specj,this_nsfL,this_nsfR)

    use datatypes
    use numbers
    use global_module, ONLY: nspin, spin_factor
    use GenBlas, ONLY: axpy, copy, scal, gemm
    use blip, ONLY: blip_info
    use support_spec_format, ONLY: support_function
    use GenComms, ONLY: cq_abort, mtime, inode, ionode
    use mult_module, ONLY: store_matrix_value, scale_matrix_value, return_matrix_value_pos, store_matrix_value_pos, matrix_pos
    use species_module, ONLY: nsf_species
    use nlpf2blip, ONLY: do_blip_integrals, do_d2blip_integrals

    implicit none

    ! Shared Variables
    integer :: speci, specj, np, nn, ip, matS, j_in_halo,num, matT

    type(support_function) :: blipL_co
    type(support_function) :: blipR_co
    type(support_function) :: blip_grad
    real(double) :: dx, dy, dz
    integer ::  this_nsfL, this_nsfR
    real(double), dimension(this_nsfL,this_nsfR,nspin) :: dataK, dataM12

    ! Local Variables
    integer, parameter :: MAX_D = 4

    real(double) ::  FAC(-MAX_D:MAX_D,3), D2FAC(-MAX_D:MAX_D,3)

    real(double), allocatable, dimension(:) :: work1, work2, work3, work4, work5, work6
    real(double), allocatable, dimension(:,:) :: temp, temp2
    real(double) :: tmp, facx, facy, facz, factot, t0, t1, mat_val, bgv

    integer :: ix, iy, iz, offset, l, at, nsf1, nsf2, stat, i1, i2, m, x,y,z, xmin, xmax, ymin, ymax, zmin, zmax, &
         idx, idy, idz, wheremat, nx, ny, nz, ir, spin

    bgv = blip_info(speci)%SupportGridSpacing*blip_info(speci)%SupportGridSpacing*blip_info(speci)%SupportGridSpacing
    dx=-dx
    dy=-dy
    dz=-dz
    allocate(temp(this_nsfR,this_nsfL),temp2(this_nsfR,this_nsfL))
    temp = zero
    temp2 = zero
    idx = floor(abs(dx/blip_info(speci)%SupportGridSpacing))
    if(dx<zero) then
       idx=-idx
       if(abs(dx-real(idx,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idx=idx-1
    else
       if(abs(dx-real(idx,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idx=idx+1
    end if
    idy = floor(abs(dy/blip_info(speci)%SupportGridSpacing))
    if(dy<zero) then
       idy=-idy
       if(abs(dy-real(idy,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idy=idy-1
    else
       if(abs(dy-real(idy,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idy=idy+1
    end if
    idz = floor(abs(dz/blip_info(speci)%SupportGridSpacing))
    if(dz<zero) then
       idz=-idz
       if(abs(dz-real(idz,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idz=idz-1
    else
       if(abs(dz-real(idz,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idz=idz+1
    end if
    call do_blip_integrals(FAC,dx/blip_info(speci)%SupportGridSpacing-real(idx,double), &
         dy/blip_info(speci)%SupportGridSpacing-real(idy,double),dz/blip_info(speci)%SupportGridSpacing-real(idz,double))
    call do_d2blip_integrals(D2FAC,dx/blip_info(speci)%SupportGridSpacing-real(idx,double), &
         dy/blip_info(speci)%SupportGridSpacing-real(idy,double),dz/blip_info(speci)%SupportGridSpacing-real(idz,double))
    allocate(work1(blip_info(speci)%FullArraySize*this_nsfR),work2(blip_info(speci)%FullArraySize*this_nsfR), &
         work3(blip_info(speci)%FullArraySize*this_nsfR),work4(blip_info(speci)%FullArraySize*this_nsfR), &
         work5(blip_info(speci)%FullArraySize*this_nsfR),work6(blip_info(speci)%FullArraySize*this_nsfR), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating arrays for onsite S elements: ",blip_info(speci)%FullArraySize,this_nsfR)

    ! first, we copy the blip functions for this atom onto a cubic grid;
    ! we make this grid 'too big' in order to have a fast routine below.
    work1 = zero
    offset = blip_info(speci)%BlipArraySize+1+3
    xmin = -blip_info(specj)%BlipArraySize
    xmax =  blip_info(specj)%BlipArraySize
    ymin = -blip_info(specj)%BlipArraySize
    ymax =  blip_info(specj)%BlipArraySize
    zmin = -blip_info(specj)%BlipArraySize
    zmax =  blip_info(specj)%BlipArraySize
    do ix = xmin,xmax
       do iy = ymin,ymax
          do iz = zmin,zmax
             l = blip_info(specj)%blip_number(ix,iy,iz)
             if (l.ne.0) then
                at = (((iz+offset)*blip_info(speci)%OneArraySize + (iy+offset))*blip_info(speci)%OneArraySize + &
                     (ix+offset)) * this_nsfR
                do nsf1 = 1,this_nsfR
                   work1(nsf1+at) = blipR_co%supp_func(nsf1)%coefficients(l)
                enddo
             end if
          end do
       end do
    end do
    
    ! now, for each direction in turn, we need to apply a number of
    ! 'spreading operations'. Do z first... put blip(z) in 2, blip.d2blip into 3
    call copy(blip_info(speci)%FullArraySize*this_nsfR,work1,1,work2,1)
    call scal(blip_info(speci)%FullArraySize*this_nsfR,FAC(0,3),work2,1)
    call copy(blip_info(speci)%FullArraySize*this_nsfR,work1,1,work3,1)
    call scal(blip_info(speci)%FullArraySize*this_nsfR,D2FAC(0,3),work3,1)
    do iz = 1, MAX_D
       offset = iz * blip_info(speci)%OneArraySize * blip_info(speci)%OneArraySize * this_nsfR
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), FAC(iz,3), &
            work1(1:), 1, work2(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), FAC(-iz,3), &
            work1(1+offset:), 1, work2(1:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), D2FAC(iz,3), &
            work1(1:), 1, work3(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), D2FAC(-iz,3), &
            work1(1+offset:), 1, work3(1:), 1 )
    end do
    ! now do y : put blip(y).blip(z) in 4blip.d2blip into 5
    call copy(blip_info(speci)%FullArraySize*this_nsfR,work2,1,work4,1)
    call scal(blip_info(speci)%FullArraySize*this_nsfR,FAC(0,2),work4,1)
    call copy(blip_info(speci)%FullArraySize*this_nsfR,work2,1,work5,1)
    call scal(blip_info(speci)%FullArraySize*this_nsfR,D2FAC(0,2),work5,1)
    call axpy(blip_info(speci)%FullArraySize*this_nsfR,FAC(0,2),work3,1,work5,1)    
    do iy = 1, MAX_D
       offset = iy * blip_info(speci)%OneArraySize * this_nsfR
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), FAC(iy,2), &
            work2(1:), 1, work4(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), FAC(-iy,2), &
            work2(1+offset:), 1, work4(1:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), D2FAC(iy,2), &
            work2(1:), 1, work5(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), D2FAC(-iy,2), &
            work2(1+offset:), 1, work5(1:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), FAC(iy,2), &
            work3(1:), 1, work5(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), FAC(-iy,2), &
            work3(1+offset:), 1, work5(1:), 1 )
    end do
    ! Now x
    work6 = zero
    call axpy(blip_info(speci)%FullArraySize*this_nsfR,FAC(0,1),work4,1,work6,1)    
    ! Reuse array space for kinetic energy
    work3 = zero
    call axpy(blip_info(speci)%FullArraySize*this_nsfR,FAC(0,1),work5,1,work3,1)    
    call axpy(blip_info(speci)%FullArraySize*this_nsfR,D2FAC(0,1),work4,1,work3,1)    
    do ix = 1, MAX_D
       offset = ix * this_nsfR
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), FAC(ix,1), &
            work4(1:), 1, work6(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), FAC(-ix,1), &
            work4(1+offset:), 1, work6(1:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), FAC(ix,1), &
            work5(1:), 1, work3(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), FAC(-ix,1), &
            work5(1+offset:), 1, work3(1:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), D2FAC(ix,1), &
            work4(1:), 1, work3(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsfR-offset), D2FAC(-ix,1), &
            work4(1+offset:), 1, work3(1:), 1 )
    end do
    ! work6 now holds the results of S integral, work3 holds the results of the T integral
    ! and now get the matrix elements by multiplication...
    ! Now set work1 to left SF coeffs
    deallocate(work1)
    allocate(work1(blip_info(speci)%FullArraySize*this_nsfL))
    work1 = zero
    offset = blip_info(speci)%BlipArraySize+1+3
    xmin = -blip_info(specj)%BlipArraySize
    xmax =  blip_info(specj)%BlipArraySize
    ymin = -blip_info(specj)%BlipArraySize
    ymax =  blip_info(specj)%BlipArraySize
    zmin = -blip_info(specj)%BlipArraySize
    zmax =  blip_info(specj)%BlipArraySize
    if(idx<0) then 
       xmin = xmin -3
    else if(idx>0) then 
       xmax = xmax +3
    end if
    if(idy<0) then 
       ymin = ymin -3
    else if(idy>0) then 
       ymax = ymax +3
    end if
    if(idz<0) then 
       zmin = zmin -3
    else if(idz>0) then 
       zmax = zmax +3
    end if
    ! Set coeffs for left SF and calculate blip grad contribution at the same time
    do ix = xmin,xmax
       do iy = ymin,ymax
          do iz = zmin,zmax
             if(abs(ix-idx)<=blip_info(specj)%BlipArraySize.AND.abs(iy-idy)<=blip_info(specj)%BlipArraySize.AND. &
                  abs(iz-idz)<=blip_info(specj)%BlipArraySize) then
                l = blip_info(specj)%blip_number(ix-idx,iy-idy,iz-idz)
             else
                l = 0
             end if
             if (l/=0) then
                at = (((iz+offset)*blip_info(specj)%OneArraySize + (iy+offset))*blip_info(specj)%OneArraySize + &
                     (ix+offset)) * this_nsfL
                do nsf1 = 1,this_nsfL
                   work1(nsf1+at) = blipL_co%supp_func(nsf1)%coefficients(l)
                enddo
                !if((abs(dx)>RD_ERR).AND.(abs(dy)>RD_ERR).AND.(abs(dz)>RD_ERR)) then
                at = (((iz+offset)*blip_info(specj)%OneArraySize + (iy+offset))*blip_info(specj)%OneArraySize + &
                     (ix+offset)) * this_nsfR
                do nsf2 =1,this_nsfR
                   do nsf1 = 1,this_nsfL
                      do spin=1,nspin
                         blip_grad%supp_func(nsf1)%coefficients(l) = &
                              blip_grad%supp_func(nsf1)%coefficients(l) &
                              -two*bgv*spin_factor*dataM12(nsf1,nsf2,spin)*work6(nsf2+at) &
                              +blip_info(speci)%SupportGridSpacing*spin_factor*dataK(nsf1,nsf2,spin)*work3(nsf2+at)
                      end do
                   end do
                enddo
                !end if
             end if
          end do
       end do
    end do
    call gemm('n','t',this_nsfR,this_nsfL,blip_info(speci)%OneArraySize*blip_info(speci)%OneArraySize* &
         blip_info(speci)%OneArraySize, one,work6,this_nsfR,work1,this_nsfL,zero,temp,this_nsfR )
    call gemm('n','t',this_nsfR,this_nsfL,blip_info(specj)%OneArraySize*blip_info(specj)%OneArraySize* &
         blip_info(specj)%OneArraySize, one,work3,this_nsfR,work1,this_nsfL,zero,temp2,this_nsfR )
    do i1 = 1,this_nsfL
       do i2=1,this_nsfR
          wheremat = matrix_pos(matS,ip,j_in_halo,i1,i2)
          mat_val = bgv*temp(i2,i1)
          call store_matrix_value_pos(matS,wheremat,mat_val)
          wheremat = matrix_pos(matT,ip,j_in_halo,i1,i2)
          ! Minus sign for KE: half is put in H_matrix_module
          mat_val = -blip_info(speci)%SupportGridSpacing*temp2(i2,i1)
          call store_matrix_value_pos(matT,wheremat,mat_val)
       end do
    end do
    deallocate(work1,work2,work3,work4,work5,work6, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating arrays for onsite S blip elements: ",blip_info(specj)%FullArraySize,this_nsfL)
    return
  end subroutine get_S_analytic
!!***

!!****f* S_matrix_module/get_dS_analytic_oneL *
!!
!!  NAME 
!!   get_S_analytic
!!  USAGE
!!   
!!  PURPOSE
!!   Calculates force contribution for S and KE matrices one direction at a time (hence one in name)
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
  subroutine get_dS_analytic_oneL(blipL_co, blipR_co, forS, forT, dataM12, dataK, ip, j_in_halo, &
       dx, dy, dz, speci, specj,this_nsfL,this_nsfR,direction)

    use datatypes
    use numbers
    use global_module, ONLY: nspin, spin_factor
    use GenBlas, ONLY: axpy, copy, scal, gemm
    use blip, ONLY: blip_info
    use support_spec_format, ONLY: support_function
    use GenComms, ONLY: cq_abort, mtime, inode, ionode
    use mult_module, ONLY: store_matrix_value, scale_matrix_value, return_matrix_value_pos, store_matrix_value_pos, matrix_pos
    use species_module, ONLY: nsf_species
    use nlpf2blip, ONLY: do_blip_integrals, do_dblip_integrals, do_d2blip_integrals, do_d3blip_integrals

    implicit none

    ! Shared Variables
    integer :: speci, specj, np, nn, ip, j_in_halo,num, direction

    type(support_function) :: blipL_co
    type(support_function) :: blipR_co
    real(double) :: dx, dy, dz
    real(double) :: forS, forT
    integer ::  this_nsfL, this_nsfR
    real(double), dimension(this_nsfL,this_nsfR,nspin) :: dataK, dataM12

    ! Local Variables
    integer, parameter :: MAX_D = 4

    real(double) ::  FAC(-MAX_D:MAX_D,3), DFAC(-MAX_D:MAX_D,3), D2FAC(-MAX_D:MAX_D,3), D3FAC(-MAX_D:MAX_D,3)

    real(double), allocatable, dimension(:) :: work1, work2, work3, work4, work5, work6
    real(double), allocatable, dimension(:,:) :: temp, temp2
    real(double) :: tmp, facx, facy, facz, factot, t0, t1, mat_val

    integer :: ix, iy, iz, offset, l, at, nsf1, nsf2, stat, i1, i2, m, x,y,z, xmin, xmax, ymin, ymax, zmin, zmax, &
         idx, idy, idz, wheremat, nx, ny, nz, ir, dir, spin

    allocate(temp(this_nsfR,this_nsfL),temp2(this_nsfR,this_nsfL))
    temp = zero
    temp2 = zero
    idx = floor(abs(dx/blip_info(speci)%SupportGridSpacing))
    if(dx<zero) then
       idx=-idx
       if(abs(dx-real(idx,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idx=idx-1
    else
       if(abs(dx-real(idx,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idx=idx+1
    end if
    idy = floor(abs(dy/blip_info(speci)%SupportGridSpacing))
    if(dy<zero) then
       idy=-idy
       if(abs(dy-real(idy,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idy=idy-1
    else
       if(abs(dy-real(idy,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idy=idy+1
    end if
    idz = floor(abs(dz/blip_info(speci)%SupportGridSpacing))
    if(dz<zero) then
       idz=-idz
       if(abs(dz-real(idz,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idz=idz-1
    else
       if(abs(dz-real(idz,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idz=idz+1
    end if
    call do_blip_integrals(FAC,dx/blip_info(speci)%SupportGridSpacing-real(idx,double), &
         dy/blip_info(speci)%SupportGridSpacing-real(idy,double),dz/blip_info(speci)%SupportGridSpacing-real(idz,double))
    call do_dblip_integrals(DFAC,dx/blip_info(speci)%SupportGridSpacing-real(idx,double), &
         dy/blip_info(speci)%SupportGridSpacing-real(idy,double),dz/blip_info(speci)%SupportGridSpacing-real(idz,double))
    call do_d2blip_integrals(D2FAC,dx/blip_info(speci)%SupportGridSpacing-real(idx,double), &
         dy/blip_info(speci)%SupportGridSpacing-real(idy,double),dz/blip_info(speci)%SupportGridSpacing-real(idz,double))
    call do_d3blip_integrals(D3FAC,dx/blip_info(speci)%SupportGridSpacing-real(idx,double), &
         dy/blip_info(speci)%SupportGridSpacing-real(idy,double),dz/blip_info(speci)%SupportGridSpacing-real(idz,double))
    allocate(work1(blip_info(speci)%FullArraySize*this_nsfL),work2(blip_info(speci)%FullArraySize*this_nsfL), &
         work3(blip_info(speci)%FullArraySize*this_nsfL),work4(blip_info(speci)%FullArraySize*this_nsfL), &
         work5(blip_info(speci)%FullArraySize*this_nsfL),work6(blip_info(speci)%FullArraySize*this_nsfL), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating arrays for onsite S elements: ",blip_info(speci)%FullArraySize,this_nsfL)
    ! first, we copy the blip functions for this atom onto a cubic grid;
    ! we make this grid 'too big' in order to have a fast routine below.
    work1 = zero
    offset = blip_info(speci)%BlipArraySize+1+3
    xmin = -blip_info(specj)%BlipArraySize
    xmax =  blip_info(specj)%BlipArraySize
    ymin = -blip_info(specj)%BlipArraySize
    ymax =  blip_info(specj)%BlipArraySize
    zmin = -blip_info(specj)%BlipArraySize
    zmax =  blip_info(specj)%BlipArraySize
    do ix = xmin,xmax
       do iy = ymin,ymax
          do iz = zmin,zmax
             l = blip_info(specj)%blip_number(ix,iy,iz)
             if (l.ne.0) then
                at = (((iz+offset)*blip_info(speci)%OneArraySize + (iy+offset))*blip_info(speci)%OneArraySize + &
                     (ix+offset)) * this_nsfL
                do nsf1 = 1,this_nsfL
                   work1(nsf1+at) = blipL_co%supp_func(nsf1)%coefficients(l)
                enddo
             end if
          end do
       end do
    end do
    ! now, for each direction in turn, we need to apply a number of
    ! 'spreading operations'. Do z first... put blip(z) in 2, and
    ! del2blip(z) in 3
    if(direction==3) then
       call copy(blip_info(speci)%FullArraySize*this_nsfL,work1,1,work2,1)
       call scal(blip_info(speci)%FullArraySize*this_nsfL,DFAC(0,3),work2,1)
       call copy(blip_info(speci)%FullArraySize*this_nsfL,work1,1,work3,1)
       call scal(blip_info(speci)%FullArraySize*this_nsfL,D3FAC(0,3),work3,1)
       do iz = 1, MAX_D
          offset = iz * blip_info(speci)%OneArraySize * blip_info(speci)%OneArraySize * this_nsfL
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), DFAC(iz,3), &
               work1(1:), 1, work2(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), DFAC(-iz,3), &
               work1(1+offset:), 1, work2(1:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), D3FAC(iz,3), &
               work1(1:), 1, work3(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), D3FAC(-iz,3), &
               work1(1+offset:), 1, work3(1:), 1 )
       end do
    else
       call copy(blip_info(speci)%FullArraySize*this_nsfL,work1,1,work2,1)
       call scal(blip_info(speci)%FullArraySize*this_nsfL,FAC(0,3),work2,1)
       call copy(blip_info(speci)%FullArraySize*this_nsfL,work1,1,work3,1)
       call scal(blip_info(speci)%FullArraySize*this_nsfL,D2FAC(0,3),work3,1)
       do iz = 1, MAX_D
          offset = iz * blip_info(speci)%OneArraySize * blip_info(speci)%OneArraySize * this_nsfL
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), FAC(iz,3), &
               work1(1:), 1, work2(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), FAC(-iz,3), &
               work1(1+offset:), 1, work2(1:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), D2FAC(iz,3), &
               work1(1:), 1, work3(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), D2FAC(-iz,3), &
               work1(1+offset:), 1, work3(1:), 1 )
       end do
    end if
    ! now do y : put blip(y).blip(z) in 4,
    ! blip(y).del2blip(z) + del2blip(y).blip(z)  in 5
    if(direction==2) then
       call copy(blip_info(speci)%FullArraySize*this_nsfL,work2,1,work4,1)
       call scal(blip_info(speci)%FullArraySize*this_nsfL,DFAC(0,2),work4,1)
       call copy(blip_info(speci)%FullArraySize*this_nsfL,work2,1,work5,1)
       call scal(blip_info(speci)%FullArraySize*this_nsfL,D3FAC(0,2),work5,1)
       call axpy(blip_info(speci)%FullArraySize*this_nsfL,DFAC(0,2),work3,1,work5,1)    
       do iy = 1, MAX_D
          offset = iy * blip_info(speci)%OneArraySize * this_nsfL
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), DFAC(iy,2), &
               work2(1:), 1, work4(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), DFAC(-iy,2), &
               work2(1+offset:), 1, work4(1:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), D3FAC(iy,2), &
               work2(1:), 1, work5(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), D3FAC(-iy,2), &
               work2(1+offset:), 1, work5(1:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), DFAC(iy,2), &
               work3(1:), 1, work5(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), DFAC(-iy,2), &
               work3(1+offset:), 1, work5(1:), 1 )
       end do
    else
       call copy(blip_info(speci)%FullArraySize*this_nsfL,work2,1,work4,1)
       call scal(blip_info(speci)%FullArraySize*this_nsfL,FAC(0,2),work4,1)
       call copy(blip_info(speci)%FullArraySize*this_nsfL,work2,1,work5,1)
       call scal(blip_info(speci)%FullArraySize*this_nsfL,D2FAC(0,2),work5,1)
       call axpy(blip_info(speci)%FullArraySize*this_nsfL,FAC(0,2),work3,1,work5,1)    
       do iy = 1, MAX_D
          offset = iy * blip_info(speci)%OneArraySize * this_nsfL
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), FAC(iy,2), &
               work2(1:), 1, work4(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), FAC(-iy,2), &
               work2(1+offset:), 1, work4(1:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), D2FAC(iy,2), &
               work2(1:), 1, work5(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), D2FAC(-iy,2), &
               work2(1+offset:), 1, work5(1:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), FAC(iy,2), &
               work3(1:), 1, work5(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), FAC(-iy,2), &
               work3(1+offset:), 1, work5(1:), 1 )
       end do
    end if
    ! Now x
    work6 = zero  ! x
    ! Reuse array space for kinetic energy
    work3 = zero  ! w6,  x
    if(direction==1) then
       call axpy(blip_info(speci)%FullArraySize*this_nsfL,DFAC(0,1),work4,1,work6,1)    
       call axpy(blip_info(speci)%FullArraySize*this_nsfL,DFAC(0,1),work5,1,work3,1)    
       call axpy(blip_info(speci)%FullArraySize*this_nsfL,D3FAC(0,1),work4,1,work3,1)    
       do ix = 1, MAX_D
          offset = ix * this_nsfL
          ! dS
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), DFAC(ix,1), &
               work4(1:), 1, work6(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), DFAC(-ix,1), &
               work4(1+offset:), 1, work6(1:), 1 )
          ! dT
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), DFAC(ix,1), &
               work5(1:), 1, work3(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), DFAC(-ix,1), &
               work5(1+offset:), 1, work3(1:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), D3FAC(ix,1), &
               work4(1:), 1, work3(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), D3FAC(-ix,1), &
               work4(1+offset:), 1, work3(1:), 1 )
       end do
    else
       call axpy(blip_info(speci)%FullArraySize*this_nsfL,FAC(0,1),work4,1,work6,1)    
       call axpy(blip_info(speci)%FullArraySize*this_nsfL,FAC(0,1),work5,1,work3,1)    
       call axpy(blip_info(speci)%FullArraySize*this_nsfL,D2FAC(0,1),work4,1,work3,1)    
       do ix = 1, MAX_D
          offset = ix * this_nsfL
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), FAC(ix,1), &
               work4(1:), 1, work6(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), FAC(-ix,1), &
               work4(1+offset:), 1, work6(1:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), FAC(ix,1), &
               work5(1:), 1, work3(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), FAC(-ix,1), &
               work5(1+offset:), 1, work3(1:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), D2FAC(ix,1), &
               work4(1:), 1, work3(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsfL-offset), D2FAC(-ix,1), &
               work4(1+offset:), 1, work3(1:), 1 )
       end do
    end if
    ! work6 now holds the results of S integral, work3 holds the results of the T integral
    ! and now get the matrix elements by multiplication...
    ! Now set work1 to Projector coeffs
    deallocate(work1)
    allocate(work1(blip_info(speci)%FullArraySize*this_nsfR))
    work1 = zero
    offset = blip_info(speci)%BlipArraySize+1+3
    xmin = -blip_info(specj)%BlipArraySize
    xmax =  blip_info(specj)%BlipArraySize
    ymin = -blip_info(specj)%BlipArraySize
    ymax =  blip_info(specj)%BlipArraySize
    zmin = -blip_info(specj)%BlipArraySize
    zmax =  blip_info(specj)%BlipArraySize
    if(idx<0) then 
       xmin = xmin -3
    else if(idx>0) then 
       xmax = xmax +3
    end if
    if(idy<0) then 
       ymin = ymin -3
    else if(idy>0) then 
       ymax = ymax +3
    end if
    if(idz<0) then 
       zmin = zmin -3
    else if(idz>0) then 
       zmax = zmax +3
    end if
    do ix = xmin,xmax
       do iy = ymin,ymax
          do iz = zmin,zmax
             if(abs(ix-idx)<=blip_info(specj)%BlipArraySize.AND.abs(iy-idy)<=blip_info(specj)%BlipArraySize.AND. &
                  abs(iz-idz)<=blip_info(specj)%BlipArraySize) then
                l = blip_info(specj)%blip_number(ix-idx,iy-idy,iz-idz)
             else
                l = 0
             end if
             if (l/=0) then
                at = (((iz+offset)*blip_info(specj)%OneArraySize + (iy+offset))*blip_info(specj)%OneArraySize + &
                     (ix+offset)) * this_nsfR
                do nsf1 = 1,this_nsfR
                   work1(nsf1+at) = blipR_co%supp_func(nsf1)%coefficients(l)
                enddo
             end if
          end do
       end do
    end do
    call gemm('n','t',this_nsfR,this_nsfL,blip_info(speci)%OneArraySize*blip_info(speci)%OneArraySize* &
         blip_info(speci)%OneArraySize, one,work1,this_nsfR,work6,this_nsfL,zero,temp(:,:),this_nsfR )
    call gemm('n','t',this_nsfR,this_nsfL,blip_info(specj)%OneArraySize*blip_info(specj)%OneArraySize* &
         blip_info(specj)%OneArraySize, one,work1,this_nsfR,work3,this_nsfL,zero,temp2(:,:),this_nsfR )
    forS = zero
    forT = zero
    t1 = mtime()
    t0 = t1
    tmp = zero
    do i1 = 1,this_nsfL
       do i2=1,this_nsfR
          do spin =1,nspin
             forS = forS -two*blip_info(speci)%SupportGridSpacing*blip_info(speci)%SupportGridSpacing* &
                  temp(i2,i1)*spin_factor*dataM12(i1,i2,spin)
             forT = forT +temp2(i2,i1)*spin_factor*dataK(i1,i2,spin)
          end do
       end do
    end do
    deallocate(work1,work2,work3,work4,work5,work6, STAT=stat)
    deallocate(temp,temp2, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating arrays for onsite S blip elements: ",blip_info(specj)%FullArraySize,this_nsfL)
    return
  end subroutine get_dS_analytic_oneL
!!***
end module S_matrix_module
