! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: exx_kernel_default.f90 XX year-mm-dd XX:XX:XXZ Lionel $
! -----------------------------------------------------------
! Module exx_kernel_default.f90
! -----------------------------------------------------------
! Code area 13: EXX
! -----------------------------------------------------------

!!****h* Conquest/exx_kernel_default *
!!  NAME
!!   exx_kernel_default
!!
!!  PURPOSE
!!   Holds the routines to compute exact exchange matrix (matX)
!!   Parallelized version deriving from multiply_kernel_default
!!
!!  AUTHOR
!!   L.A.Truflandier/D.R.Bowler
!!
!!  USES
!!   mpi, GenComms, comms_module
!!   timer_module, timer_stdclocks_module
!!   exx_types, exx_io, exx_memory
!!   exx_module, exx_evalpao
!!   global_module, basic_types, memory_module
!!   species_module, primary_module, cover_module
!!   group_module, energy, datatypes, numbers
!!   multiply_module, matrix_module, matrix_data
!!   mult_module, matrix_comms_module
!!
!!  CREATION DATE
!!   2014/05/29
!!
!!  MODIFICATION HISTORY
!!   2015/06/08 lat
!!    - Added spin 
!!
!!  SOURCE
!!
module exx_kernel_default

  use datatypes

  use GenComms,                  only: my_barrier, cq_abort, mtime
  use GenComms,                  only: inode, ionode, myid, root
  use timer_module,              only: start_timer, stop_timer, print_timer, lun_tmr
  use timer_module,              only: start_backtrace, stop_backtrace, cq_timer
  use timer_stdclocks_module,    only: tmr_std_exx

  !**<lat>** ISF Poisson solver Will be available in the forthcoming version 
  !use Poisson_Solver,            only: PSolver, createKernel, gequad
 
  use exx_types,                 only: reckernel_3d, fftwrho3d, exx_debug
  use exx_io

  implicit none 

  ! Area identification
  integer, parameter, private :: area = 13

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = &
       "$Id: exx_kernel_default.f90 XX year-mm-dd XX:XX:XXZ lionel $"

!!***
  
contains
 
  !!****f* exx_kernel_default/get_X_matrix *
  !!
  !!  NAME
  !!   get_X_matrix
  !!
  !!  PURPOSE
  !!   Compute X matrix 
  !! 
  !!  INPUTS
  !!
  !!  AUTHOR
  !!   L.A.Truflandier
  !!
  !!  CREATION DATE
  !!   2013/10/30
  !!
  !!  MODIFICATION HISTORY
  !!   2018/10/08 07:59 dave
  !!    Adding recv_part variable to count partitions received from processes
  !!    to follow update of multiplies (to conform to MPI standard for
  !!    tags)
  !!   2018/11/13 18:00 nakata
  !!    Removed matS which was passed but not used
  !!  SOURCE
  !!
  subroutine get_X_matrix( exxspin, level )

    use numbers
    use global_module
    use matrix_module
    use basic_types
    use matrix_comms_module
    use comms_module
    use mpi
    use GenComms,        only: my_barrier, cq_abort, mtime
    use multiply_module, only: prefetch
    !
    use energy,         only: exx_energy
    !
    use primary_module, only: bundle
    use cover_module,   only: BCS_parts
    use matrix_data,    only: mat, Hrange, Srange, Xrange, SXrange, halo, rcut
    use mult_module,    only: S_X_SX, mat_p, mult 
    use mult_module,    only: matX, matK, matrix_scale, matrix_trace
    use mult_module,    only: matrix_product_trace, matrix_product_trace_length 
    use mult_module,    only: return_matrix_value,  matrix_pos
    use mult_module,    only: store_matrix_value,   store_matrix_value_pos
    use memory_module,  only: write_mem_use
    use species_module, only: nsf_species
    !
    use exx_memory,     only: exx_mem_alloc
    !
    use exx_types, only: prim_atomic_data, neigh_atomic_data,      &
         tmr_std_exx_evalpao, &
         tmr_std_exx_setup,   tmr_std_exx_fetch_K, &
         tmr_std_exx_fetch,   tmr_std_exx_accumul, &
         tmr_std_exx_matmult, tmr_std_exx_poisson, &
         tmr_std_exx_allocat, tmr_std_exx_dealloc, &
         tmr_std_exx_fetch,   tmr_std_exx_kernel,  &
         unit_matrix_write, unit_output_write,     & 
         unit_timers_write, unit_screen_write,     &
         unit_memory_write

    use exx_types, only: grid_spacing, r_int, extent, &
         ewald_alpha, ewald_charge, ewald_rho, ewald_pot,   &
         pulay_radius, p_omega, p_ngauss, p_gauss, w_gauss, &
         exx_psolver,p_scheme, isf_order, ngrid, kernel,    &
         exx_alloc, exx_mem, exx_phil, exx_screen_pao,      &
         exx_total_time
    !
    use exx_module,only: exx_scal_rho_3d, exx_ewald_rho, exx_ewald_pot
    !
    !**<lat>** ISF Poisson solver Will be available in the forthcoming version 
    !use Poisson_Solver, only: createBeylkin
    !
    implicit none

    integer, optional :: level

    ! Local variables
    integer :: exxspin
    integer :: lab_const
    integer :: invdir,ierr,kpart,ind_part,ncover_yz,n_which,ipart,nnode
    integer :: icall,n_cont,kpart_next,ind_partN,k_off
    integer :: icall2,stat,ilen2,lenb_rem
    ! Remote variables to be allocated
    integer(integ),allocatable :: ibpart_rem(:)
    real(double),  allocatable :: b_rem(:)

    ! Remote variables which will point to part_array
    integer(integ),pointer :: nbnab_rem(:)
    integer(integ),pointer :: ibseq_rem(:)
    integer(integ),pointer :: ibind_rem(:)
    integer(integ),pointer :: ib_nd_acc_rem(:)
    integer(integ),pointer :: npxyz_rem(:)
    integer(integ),pointer :: ibndimj_rem(:)

    ! Arrays for remote variables to point to



    integer, target :: part_array(3*mult(S_X_SX)%parts%mx_mem_grp+ &
         5*mult(S_X_SX)%parts%mx_mem_grp*mult(S_X_SX)%bmat(  exxspin  )%mx_abs)


    integer, allocatable, dimension(:) :: recv_part

    integer, dimension(:), allocatable  :: nreqs
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    integer :: offset,sends,i,j

    logical      :: flag,call_flag
    real(double) :: t0,t1

    integer :: iprim, np, nsf1, nsf2, gcspart
    integer :: maxsuppfuncs

    type(prim_atomic_data)  :: ia !i_alpha
    type(neigh_atomic_data) :: jb !j_beta
    type(neigh_atomic_data) :: kg !k_gamma
    type(neigh_atomic_data) :: ld !l_delta
    !
    real(double)      :: xyz_ghost(3), r_ghost, tmp
    integer           :: wheremat
    !
    integer           :: unit1, unit2, unit3, unit4, unit5
    character(len=20) :: filename1, filename2, filename3, filename4, filename5, filename6, filename7

    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level
    !
    !==============================================================================================================

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='get_X_matrix',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$
 
    call start_timer(tmr_std_exx)   
    call start_timer(tmr_std_exx_setup)    
    !
    if ( iprint_exx > 3 ) then
       !call io_assign(unit_matrix_write)
       !unit1 = unit_matrix_write
       !unit1 = 100001
       !unit_matrix_write = unit1
       !call get_file_name('exx_matrix.debug',numprocs,inode,filename1)
       !open(unit1,file=filename1)

       !call io_assign(unit_output_write)
       !unit2 = unit_output_write
       !unit2 = 100002
       !unit_output_write = unit2
       !call get_file_name('exx_output.debug',numprocs,inode,filename2)
       !open(unit2,file=filename2)

       call io_assign(unit_timers_write)
       call get_file_name('exx_timers',numprocs,inode,filename3)
       open(unit_timers_write,file=filename3)
       !
       !
       call io_assign(unit_memory_write)
       call get_file_name('exx_memory',numprocs,inode,filename4)
       open(unit_memory_write,file=filename4)

       !call io_assign(unit_screen_write)
       !unit5 = unit_screen_write
       !unit5 = 100005
       !unit_screen_write = unit5
       !call get_file_name('exx_screen.debug',numprocs,inode,filename5)
       !open(unit5,file=filename5)

       !call io_assign(unit_exx_debug)
       !call get_file_name('Exx_debug.debug',numprocs,inode,filename6)
       !open(unit_exx_debug,file=filename6)       
       !call exx_write_head(unit_exx_debug,inode,bundle%groups_on_node) 
       !
       !call exx_global_write()
       !
    end if
    !
    maxsuppfuncs = maxval(nsf_species)
    !
    call exx_mem_alloc(extent,0,0,'work_3d' ,'alloc')
    !
    !DRB! This appears to set up the different Poisson solvers
    if (exx_psolver == 'fftw')     then
       call exx_mem_alloc(extent,0,0,'fftw_3d','alloc')  
       call exx_mem_alloc(extent,0,0,'reckernel_3d','alloc')

       poisson_fftw: select case(p_scheme)          

       case('default')
          call exx_scal_rho_3d(inode,extent,r_int,p_scheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss)

       case('ewald')
          call exx_mem_alloc(extent,0,0,'ewald_3d','alloc')
          call exx_ewald_rho(ewald_rho,extent,ewald_alpha,r_int)
          call exx_ewald_pot(ewald_pot,extent,ewald_alpha,r_int)          
          call exx_scal_rho_3d(inode,extent,r_int,p_scheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss)

       case('pulay')
          call exx_scal_rho_3d(inode,extent,r_int,p_scheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss)

       case('yukawa')
          call exx_scal_rho_3d(inode,extent,r_int,p_scheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss)

       case('gauss')
          call cq_abort('EXX: Gaussian representation if 1/r for solving &
              &the Poisson equation &
              &is currently under testing...')
          !call createBeylkin(p_gauss,w_gauss,r_int)     
          !call exx_scal_rho_3d(inode,extent,r_int,p_scheme,pulay_radius, &
          !     p_omega,p_ngauss,p_gauss,w_gauss)

       case default
          call exx_scal_rho_3d(inode,extent,r_int,p_scheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss)

       end select poisson_fftw

    else if (exx_psolver == 'isf') then
       call cq_abort('EXX: ISF Poisson solver is not available yet ; &
              &under optimisation...')
       !call exx_mem_alloc(extent,0,0,'isf_rho','alloc')
       !call createKernel('F',ngrid,ngrid,ngrid,grid_spacing,grid_spacing,grid_spacing,isf_order,&
       !     0,1,kernel)       

    end if
    call stop_timer(tmr_std_exx_setup,.true.)
    !    
    !DRB! Zero matrix - spin polarised possible - fix this later ? 
    !DRB! Where is matX allocated ? immi
    !call matrix_scale(zero, matX(1))
    !
    !
    call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_i','alloc')       
    call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_j','alloc')
    call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_k','alloc')
    call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_l','alloc')
    !
    call exx_mem_alloc(extent,maxsuppfuncs,0,'Phy_k','alloc')
    call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'Ome_kj','alloc')       
    call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_kj','alloc')       
    call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_kj','alloc')
    !
    !
    !==============================================================================================================
    !
    call matrix_scale(zero,matX(  exxspin  ))
    !
    call start_timer(tmr_std_exx_kernel)
    !
    call start_timer(tmr_std_exx_allocat)
    if(iprint_mat>3.AND.myid==0) t0 = mtime()
    ! Allocate memory for the elements
    allocate(ibpart_rem(mult(S_X_SX)%parts%mx_mem_grp*mult(S_X_SX)%bmat(  exxspin  )%mx_abs),STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error allocating ibpart_rem')
    call stop_timer(tmr_std_exx_allocat,.true.)
    !
    sends = 0
    do i=1,mult(S_X_SX)%comms%inode
       if(mult(S_X_SX)%comms%np_send(i)>0) then
          do j=1,mult(S_X_SX)%comms%np_send(i)
             sends = sends+1
          end do
       end if
    end do
    !
    call start_timer(tmr_std_exx_allocat)
    allocate(nreqs(sends*2),STAT=stat)
    if(stat/=0) call cq_abort("mat_mult: Error allocating nreqs",sends,stat)
    call stop_timer(tmr_std_exx_allocat,.true.)
    allocate(recv_part(0:mult(S_X_SX)%comms%inode),STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error allocating recv_part')
    recv_part = zero
    !
    sends  = 0
    invdir = 0
    call Mquest_start_send(mult(S_X_SX),mat_p(matK(  exxspin  ))%matrix,nreqs,myid,mult(S_X_SX)%prim%mx_ngonn,sends)
    !ncover_yz=mult(S_X_SX)%gcs%ncovery*mult(S_X_SX)%gcs%ncoverz
    ncover_yz=mult(S_X_SX)%gcs%ncovery*mult(S_X_SX)%gcs%ncoverz
    !
    ! #ifdef OMP_M
    ! !$omp parallel default(none) &
    ! !$omp          shared(a, b, c, a_b_c, myid, lena, lenc, tmr_std_allocation, &
    ! !$omp                 ncover_yz, ibpart_rem, atrans, usegemm) &
    ! !$omp          private(kpart, icall, ind_part, ipart, nnode, b_rem, &
    ! !$omp                  lenb_rem, n_cont, part_array, ilen2, offset, &
    ! !$omp                  nbnab_rem, ibind_rem, ib_nd_acc_rem, ibseq_rem, &
    ! !$omp                  npxyz_rem, ibndimj_rem, k_off, icall2)
    ! !$omp do
    ! #end if
    !
    if ( exx_debug ) then
       call hf_write_info(unit1,inode,10,mult(S_X_SX)%ahalo%np_in_halo,0,0,0,0, &
            0,0,'XX',xyz_ghost,r_ghost,0,0,0,zero,zero,0)
    end if
    !
    !
    !call hf_write_info(unit1,inode,10,bundle%groups_on_node,0,0,0,0, &
    !     0,0,'XX',xyz_ghost,r_ghost,0,0,0,zero,zero,0)
    !
    xyz_ghost = zero
    r_ghost   = zero
    do kpart = 1,mult(S_X_SX)%ahalo%np_in_halo  ! Main loop
       icall=1
       ind_part = mult(S_X_SX)%ahalo%lab_hcell(kpart)    
       !
       !print*, 'inode', inode,'kpart', kpart, ind_part
       if ( exx_debug ) then
          call hf_write_info(unit1,inode,20,0,kpart,mult(S_X_SX)%ahalo%np_in_halo,0,0, &
               0,0,'XX',xyz_ghost,r_ghost,0,0,0,zero,zero,0)       
       end if
       !
       if(kpart>1) then  ! Is it a periodic image of the previous partition ?
          if(ind_part.eq.mult(S_X_SX)%ahalo%lab_hcell(kpart-1)) then
             icall=0
          else ! Get the data
             ipart = mult(S_X_SX)%parts%i_cc2seq(ind_part)
             nnode = mult(S_X_SX)%comms%neigh_node_list(kpart)
             recv_part(nnode) = recv_part(nnode)+1
             !
             if(allocated(b_rem)) deallocate(b_rem)
             if(mult(S_X_SX)%parts%i_cc2node(ind_part)==myid+1) then
                lenb_rem = mult(S_X_SX)%bmat(ipart)%part_nd_nabs
             else
                lenb_rem = mult(S_X_SX)%comms%ilen3rec(ipart,nnode)
             end if
             !
             allocate(b_rem(lenb_rem))
             !
             call prefetch(kpart,mult(S_X_SX)%ahalo,mult(S_X_SX)%comms,mult(S_X_SX)%bmat,icall, &
                  n_cont,part_array,mult(S_X_SX)%bindex,b_rem,lenb_rem,mat_p(matK(  exxspin  ))%matrix,   &
                  myid,ilen2,mx_msg_per_part,mult(S_X_SX)%parts,mult(S_X_SX)%prim,mult(S_X_SX)%gcs,&
                  (recv_part(nnode)-1)*2)
             !
             offset = 0
             nbnab_rem => part_array(offset+1:offset+n_cont)
             offset = offset+n_cont
             ibind_rem => part_array(offset+1:offset+n_cont)
             offset = offset+n_cont
             ib_nd_acc_rem => part_array(offset+1:offset+n_cont)
             offset = offset+n_cont
             ibseq_rem => part_array(offset+1:offset+ilen2)
             offset = offset+ilen2
             npxyz_rem => part_array(offset+1:offset+3*ilen2)
             offset = offset+3*ilen2
             ibndimj_rem => part_array(offset+1:offset+ilen2)
             !
             if(offset+ilen2>3*mult(S_X_SX)%parts%mx_mem_grp+ &
                  5*mult(S_X_SX)%parts%mx_mem_grp*mult(S_X_SX)%bmat(  exxspin  )%mx_abs) then
                call cq_abort('mat_mult: error pointing to part_array ',kpart)
             end if
             ! Create ibpart_rem
             call end_part_comms(myid,n_cont,nbnab_rem,ibind_rem,npxyz_rem,&
                  ibpart_rem,ncover_yz,mult(S_X_SX)%gcs%ncoverz)
          end if
       else ! Get the data
          ipart = mult(S_X_SX)%parts%i_cc2seq(ind_part)
          nnode = mult(S_X_SX)%comms%neigh_node_list(kpart)
          recv_part(nnode) = recv_part(nnode)+1
          !
          !
          if(allocated(b_rem)) deallocate(b_rem)
          if(mult(S_X_SX)%parts%i_cc2node(ind_part)==myid+1) then
             lenb_rem = mult(S_X_SX)%bmat(ipart)%part_nd_nabs
          else
             write(*,*) 'Error? ',nnode,mult(S_X_SX)%parts%i_cc2node(ind_part),myid
             lenb_rem = mult(S_X_SX)%comms%ilen3rec(ipart,nnode)
          end if
          !
          !
          call start_timer(tmr_std_exx_allocat)
          allocate(b_rem(lenb_rem))
          call stop_timer(tmr_std_exx_allocat,.true.)
          !
          !
          call prefetch(kpart,mult(S_X_SX)%ahalo,mult(S_X_SX)%comms,mult(S_X_SX)%bmat,icall, &
               n_cont,part_array,mult(S_X_SX)%bindex,b_rem,lenb_rem,mat_p(matK(  exxspin  ))%matrix,   & 
               myid,ilen2,mx_msg_per_part,mult(S_X_SX)%parts,mult(S_X_SX)%prim,mult(S_X_SX)%gcs,&
               (recv_part(nnode)-1)*2)
          !
          lenb_rem = size(b_rem)
          !
          !
          offset = 0              ; nbnab_rem => part_array(offset+1:offset+n_cont)
          offset = offset+n_cont  ; ibind_rem => part_array(offset+1:offset+n_cont)
          offset = offset+n_cont  ; ib_nd_acc_rem => part_array(offset+1:offset+n_cont)
          offset = offset+n_cont  ; ibseq_rem => part_array(offset+1:offset+ilen2)
          offset = offset+ilen2   ; npxyz_rem => part_array(offset+1:offset+3*ilen2)
          offset = offset+3*ilen2 ; ibndimj_rem => part_array(offset+1:offset+ilen2)
          !
          if(offset+ilen2>3*mult(S_X_SX)%parts%mx_mem_grp+ &
               5*mult(S_X_SX)%parts%mx_mem_grp*mult(S_X_SX)%bmat(  exxspin  )%mx_abs) then
             call cq_abort('Error pointing to part_array !',kpart)
          end if
          !          
          call end_part_comms(myid,n_cont,nbnab_rem,ibind_rem,npxyz_rem,&
               ibpart_rem,ncover_yz,mult(S_X_SX)%gcs%ncoverz)
          !
       end if ! End of the "if this isn't the first partition" loop
       k_off  = mult(S_X_SX)%ahalo%lab_hcover(kpart) ! --- offset for pbcs
       icall2 = 1
       !
       call m_kern_exx( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
            ibpart_rem,ibseq_rem,ibndimj_rem, & 
                                !atrans, &
            b_rem,  &
            mat_p(matX(  exxspin  ))%matrix,      &
            mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
            mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
            mult(S_X_SX)%prim%mx_iprim, &
                                !lena,     &
            lenb_rem, &
            mat_p(matX(  exxspin  ))%length)
    end do ! End of the kpart=1,ahalo%np_in_halo loop !
    !
    ! #ifdef OMP_M
    ! !$omp end do
    ! !$omp end parallel
    ! #end if
    !
    call start_timer(tmr_std_exx_dealloc)
    if(allocated(b_rem)) deallocate(b_rem)
    call stop_timer(tmr_std_exx_dealloc,.true.)
    !
    if(sends>0) then
       do i=1,sends
          call MPI_Wait(nreqs(i),mpi_stat,ierr)
          if(ierr/=0) call cq_abort("Error waiting for send to finish",i)
       end do
    end if
    !
    call my_barrier
    !
    call start_timer(tmr_std_exx_dealloc)
    deallocate(nreqs,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating nreqs')
    deallocate(recv_part,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating ibpart_rem')
    call stop_timer(tmr_std_exx_dealloc,.true.)
    !
    call my_barrier
    !
    call start_timer(tmr_std_exx_dealloc)
    deallocate(ibpart_rem,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating recv_part')
    call stop_timer(tmr_std_exx_dealloc,.true.)
    !
    call my_barrier
    !
    if(iprint_mat>3.AND.myid==0) then
       t1 = mtime()
       write(io_lun,*) 'mult time: ',t1-t0
    end if
    !
    call stop_timer(tmr_std_exx_kernel,.true.)
    !
    !exx = matrix_product_trace(matX(1),matK(1))
    !if (inode == ionode) then
    !   print*, 'exx energy'
    !end if
    !
    if ( exx_debug ) then
       !call exx_write_tail(unit1,inode)  
    end if
    !
    !
    call exx_mem_alloc(extent,0,0,'work_3d' ,'dealloc')
    !
    call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_i','dealloc')       
    call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_j','dealloc')
    call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_k','dealloc')    
    call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_l','dealloc')
    !
    call exx_mem_alloc(extent,maxsuppfuncs,0,'Phy_k','dealloc')
    call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'Ome_kj','dealloc')   
    call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_kj','dealloc')       
    call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_kj','dealloc')
    !
    !
    if (exx_psolver == 'fftw')  then
       call exx_mem_alloc(extent,0,0,'fftw_3d',     'dealloc')  
       call exx_mem_alloc(extent,0,0,'reckernel_3d','dealloc')
       select case(p_scheme)
       case('default')         
       case('gauss')          
       case('pulay')
       case('ewald')
          call exx_mem_alloc(extent,0,0,'ewald_3d','dealloc')
       case default
       end select

    else if (exx_psolver == 'isf') then
       call exx_mem_alloc(extent,0,0,'isf_rho','dealloc')
       deallocate(kernel)
    end if
    !
    if (inode == ionode) then
       !call write_mem_use(unit_memory_write,area_exx)
    end if
    !
    call my_barrier()
    !
    ! Timers and elapse time
    exx_total_time = &
         tmr_std_exx_evalpao%t_tot + &
         tmr_std_exx_poisson%t_tot + &
         tmr_std_exx_matmult%t_tot + &
         tmr_std_exx_accumul%t_tot + &
         tmr_std_exx_allocat%t_tot + &
         tmr_std_exx_dealloc%t_tot + &
         tmr_std_exx_setup%t_tot   + &
         tmr_std_exx_write%t_tot   + &
         tmr_std_exx_fetch%t_tot  
         !tmr_std_exx_fetch_K%t_tot 

    call stop_timer(tmr_std_exx,.true.)

    if ( iprint_exx > 3 ) then
       call print_timer(tmr_std_exx_setup,  "exx_setup   time:", unit_timers_write)    
       call print_timer(tmr_std_exx_write,  "exx_write   time:", unit_timers_write)    
       call print_timer(tmr_std_exx_kernel, "exx_kernel  time:", unit_timers_write)    
       call print_timer(tmr_std_exx_fetch,  "exx_fetch   time:", unit_timers_write)    
       !call print_timer(tmr_std_exx_fetch_K,"exx_fetch_K time:", unit_timers_write)    

       call print_timer(tmr_std_exx_evalpao,"exx_evalpao time:", unit_timers_write)    
       call print_timer(tmr_std_exx_poisson,"exx_poisson time:", unit_timers_write)
       call print_timer(tmr_std_exx_matmult,"exx_matmult time:", unit_timers_write)    
       call print_timer(tmr_std_exx_accumul,"exx_accumul time:", unit_timers_write)    
       call print_timer(tmr_std_exx_allocat,"exx_allocat time:", unit_timers_write)    
       call print_timer(tmr_std_exx_dealloc,"exx_dealloc time:", unit_timers_write)    
       call print_timer(tmr_std_exx,        "exx_total   time:", unit_timers_write)    

       write(unit=unit_timers_write,fmt='("Timing: Proc ",i6,": Time spent in ", a50, " = ", &
           &f12.5," s")') inode, 'get_X_matrix',  tmr_std_exx%t_tot

       write(unit=unit_timers_write,fmt='("Timing: Proc ",i6,": Time spent in ", a50, " = ", &
           &f12.5," s")') inode, 'timer calls',  tmr_std_exx%t_tot-exx_total_time 

       !call io_close(unit_matrix_write)
       !call io_close(unit_output_write)
       !call io_close(unit_screen_write)
       call io_close(unit_memory_write)
       call io_close(unit_timers_write)
       !
       !
    end if

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_X_matrix',echo=.true.)
!****lat>$

    return
  end subroutine get_X_matrix
  !!***

  !!****f* exx_kernel_default/m_kern_exx *
  !!
  !!  NAME
  !!   m_kern_exx
  !!
  !!  PURPOSE
  !! 
  !!  INPUTS
  !!
  !!  AUTHOR
  !!   L.A.Truflandier/D.R.Bowler
  !!
  !!  CREATION DATE
  !!   2014/05/29  
  !!
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine m_kern_exx(k_off, kpart, ib_nd_acc, ibaddr, nbnab, &
       ibpart, ibseq, bndim2, b, c, ahalo, chalo, & 
       at, mx_absb, mx_part,  &
       mx_iprim, lenb, lenc,debug)

    use numbers,        only: zero, one
    use matrix_module,  only: matrix_halo, matrix_trans
    use global_module,  only: area_exx, id_glob, species_glob
    !
    use basic_types,    only: primary_set
    use primary_module, only: bundle 
    use matrix_data,    only: mat, Hrange, SXrange, Xrange, Srange, halo
    use mult_module,    only: matK, return_matrix_value, mult, S_X_SX, mat_p, matX
    use cover_module,   only: BCS_parts
    use group_module,   only: parts 
    !
    use species_module, only: nsf_species, nlpf_species 
    !
    use exx_evalpao,    only: exx_phi_on_grid
    !
    use exx_types, only: prim_atomic_data, neigh_atomic_data, &
         tmr_std_exx_evalpao, &
         tmr_std_exx_setup, tmr_std_exx_fetch_K,   &
         tmr_std_exx_fetch, tmr_std_exx_accumul,   &
         tmr_std_exx_matmult, tmr_std_exx_poisson, &
         tmr_std_exx_allocat, tmr_std_exx_dealloc, &
         tmr_std_exx_fetch,   &
                                !unit_matrix_write, unit_output_write,     & 
                                !unit_timers_write, unit_screen_write,     &
                                !unit_memory_write, &
         grid_spacing, r_int, extent, &
         ewald_alpha, ewald_charge, ewald_rho, ewald_pot,   &
         pulay_radius, p_omega, p_ngauss, p_gauss, w_gauss, &
         exx_psolver,p_scheme, isf_order, ngrid, kernel,    &
         exx_alloc, exx_mem, exx_phil, exx_screen_pao,      &
         exx_total_time
    !
    use exx_types, only: phi_i, phi_j, phi_k, phi_l, &
         Phy_k, rho_kj, Ome_kj, vhf_kj, &
         work_in_3d, work_out_3d,       &
         exx_Kkl, exx_Kij

    use exx_module,only: exx_v_on_grid, get_halodat, get_iprimdat
    !
    !
    implicit none
    !
    ! Passed variables
    type(matrix_halo)  :: ahalo, chalo
    type(matrix_trans) :: at
    integer            :: mx_absb, mx_part, mx_iprim, lena, lenb, lenc
    integer            :: kpart, k_off
    real(double)       :: b(lenb)
    real(double)       :: c(lenc)
    integer, optional  :: debug
    !
    ! Remote indices
    integer(integ) :: ib_nd_acc(mx_part)
    integer(integ) :: ibaddr(mx_part)
    integer(integ) :: nbnab(mx_part)
    integer(integ) :: ibpart(mx_part*mx_absb)
    integer(integ) :: ibseq(mx_part*mx_absb)
    integer(integ) :: bndim2(mx_part*mx_absb)
    !
    ! Local variables
    integer :: jbnab2ch(mx_absb)  ! Automatic array
    integer :: nbkbeg, k, k_in_part, k_in_halo, j, jpart, jseq
    integer :: i, nabeg, i_in_prim, icad, nbbeg, j_in_halo, ncbeg
    integer :: n1, n2, n3, nb_nd_kbeg
    integer :: nd1, nd2, nd3
    integer :: naaddr, nbaddr, ncaddr
    integer :: lbnab2ch(mx_absb)  ! Automatic array
    integer :: l, lseq, lpart, l_in_halo
    integer :: np, ni, iprim
    !
    integer :: unit_exx_debug1
    !
    real(double), dimension(3) :: xyz_ghost = zero
    real(double), dimension(3) :: xyz_zero  = zero
    real(double), dimension(3) :: xyz_delta = zero
    real(double), dimension(3) :: xyz_ij    = zero
    real(double), dimension(3) :: xyz_kl    = zero
    real(double)               ::   r_ghost = zero
    !
    real(double)               ::   dr,dv,K_val
    real(double)               ::   exx_mat_elem
    real(double)               ::   screen_ij, range_ij
    real(double)               ::   screen_kl, range_kl
    !
    type(prim_atomic_data)  :: ia !i_alpha
    type(neigh_atomic_data) :: jb !j_beta
    type(neigh_atomic_data) :: kg !k_gamma
    type(neigh_atomic_data) :: ld !l_delta
    !
    integer                 :: maxsuppfuncs, nsf1, nsf2, nsf3
    integer                 :: r, s, t, tmp    
    !
    !
    dr = grid_spacing
    dv = dr**3
    ewald_alpha  = 0.5
    maxsuppfuncs = maxval(nsf_species)    
    !range_ij        = 0.5d0
    !range_kl        = 0.5d0
    !unit_exx_debug1 = 333
    !
    !
!!$
!!$ ****[ k loop ]****
!!$
    do k = 1, ahalo%nh_part(kpart)
       !
       k_in_halo  = ahalo%j_beg(kpart) + k - 1
       k_in_part  = ahalo%j_seq(k_in_halo)
       nbkbeg     = ibaddr     (k_in_part) 
       nb_nd_kbeg = ib_nd_acc  (k_in_part)
       nd3        = ahalo%ndimj(k_in_halo)
       call get_halodat(kg,kg,k_in_part,ahalo%i_hbeg(ahalo%lab_hcover(kpart)), &
            ahalo%lab_hcell(kpart),'k',.true.,unit_exx_debug1)
       !
       !print*, 'k',k, 'global_num',kg%global_num,'spe',kg%spec
       !
       call exx_phi_on_grid(inode,kg%global_num,kg%spec,extent, &
            xyz_zero,maxsuppfuncs,phi_k,r_int,xyz_zero)             
       !
       jbnab2ch = 0
       !print*, 'nbnab: ',nbnab(k_in_part),k_in_part
       do j = 1, nbnab(k_in_part)
          jpart = ibpart(nbkbeg+j-1) + k_off
          jseq  = ibseq (nbkbeg+j-1)
          jbnab2ch(j) = chalo%i_halo(chalo%i_hbeg(jpart)+jseq-1)
          !print*, 'jbnab2ch',j, jbnab2ch(j)
       end do
       !
       nbbeg = nb_nd_kbeg
       !
       !print*, 'size jbnab2ch', size(jbnab2ch)
       !print*, 'jbnab2ch', jbnab2ch
       !print*
       !
       Phy_k = zero
       !
!!$
!!$ ****[ l do loop ]****
!!$
       do l = 1, nbnab(k_in_part)
          !l_in_halo = lbnab2ch(l)                
          lpart = ibpart(nbkbeg+l-1) + k_off
          lseq  = ibseq (nbkbeg+l-1)
          call get_halodat(ld,kg,lseq,chalo%i_hbeg(lpart),         &
               BCS_parts%lab_cell(BCS_parts%inv_lab_cover(lpart)), &
                           'l',.true.,unit_exx_debug1)
          !
          !write(*,*) 'l',l, 'global_num',ld%global_num,'spe',ld%spec
          !
!!$
!!$ ****[ <kl> screening ]****
!!$
          !xyz_kl    = kg%xyz - ld%xyz
          !screen_kl = sqrt(dot_product(xyz_kl,xyz_kl))
          !if ( screen_kl < range_kl ) then
          !
          call exx_phi_on_grid(inode,ld%global_num,ld%spec,extent,     &
                               ld%xyz,maxsuppfuncs,phi_l,r_int,xyz_zero)             
          !
          do nsf2 = 1, ld%nsup
             !
             nbaddr = nbbeg + kg%nsup * (nsf2 - 1)
             !
             do nsf1 = 1, kg%nsup                         
                !
                !call start_timer(tmr_std_exx_fetch_K)
                if (exx_Kkl) then
                   K_val = b(nbaddr+nsf1-1)
                else
                   K_val = real(1,double)
                end if
                !
                !write(*,'(7I4,F12.8)') kpart,i,j,k,l,nsf2,nsf1,K_val
                !
                !call stop_timer(tmr_std_exx_fetch_K,.true.)
                
                call start_timer(tmr_std_exx_accumul)
                Phy_k(:,:,:,nsf1) = Phy_k(:,:,:,nsf1) + K_val*phi_l(:,:,:,nsf2) 
                call stop_timer(tmr_std_exx_accumul,.true.)
                
             end do
          end do
          !
          nbbeg = nbbeg + ld%nsup*kg%nsup !nd3 * nd2                
          !
       !end if !( screen kl )
       end do ! End of l = 1, nbnab(k_in_part)
!!$
!!$ ****[ i loop ]****
!!$
       do i = 1, at%n_hnab(k_in_halo)
          i_in_prim = at%i_prim(at%i_beg(k_in_halo)+i-1)
          nd1 = ahalo%ndimi      (i_in_prim)
          ni  = bundle%iprim_seq (i_in_prim)
          np  = bundle%iprim_part(i_in_prim)
          icad  = (i_in_prim - 1) * chalo%ni_in_halo !***
          !nbbeg = nb_nd_kbeg
          !
          !print*, 'i_in_prim', i_in_prim, 'nd1', nd1, 'ni', ni, 'np', np
          !
          call get_iprimdat(ia,kg,ni,i_in_prim,np,.true.,unit_exx_debug1)          
          !
          !print*, 'i',i, 'global_num',ia%ip,'spe',ia%spec
          !
          call exx_phi_on_grid(inode,ia%ip,ia%spec,extent, &
               ia%xyz,maxsuppfuncs,phi_i,r_int,xyz_zero)             
          !
          !print*, size(chalo%i_h2d), shape(chalo%i_h2d)
          ! 
!!$
!!$ ****[ j loop ]****
!!$
          do j = 1, nbnab(k_in_part)!mat(np,Xrange)%n_nab(ni)                   
             nbbeg     = nb_nd_kbeg
             j_in_halo = jbnab2ch(j) !***
             !
             !print*, j, icad, j_in_halo
             !   
             if ( j_in_halo /= 0 ) then
                !
                ncbeg = chalo%i_h2d(icad + j_in_halo) !***
                !
                if ( ncbeg /= 0 ) then 
                   jpart = ibpart(nbkbeg+j-1) + k_off
                   jseq  = ibseq (nbkbeg+j-1)
                   call get_halodat(jb,kg,jseq,chalo%i_hbeg(jpart),                     &
                                    BCS_parts%lab_cell(BCS_parts%inv_lab_cover(jpart)), &
                           'j',.true.,unit_exx_debug1)
                   !
!!$
!!$ ****[ <ij> screening ]****
!!$
                   !xyz_ij    = ia%xyz - jb%xyz
                   !screen_ij = sqrt(dot_product(xyz_ij,xyz_ij))
                   !print*, screen_ij
                   !if ( screen_ij < range_ij ) then
                   !write(*,*) 'j',j, 'global_num',jb%global_num,'spe',jb%spec
                   !
                   call exx_phi_on_grid(inode,jb%global_num,jb%spec,extent, &
                                        jb%xyz,maxsuppfuncs,phi_j,r_int,xyz_zero)             
                   !
                   !Ome_kj = zero
                   !
                   call start_timer(tmr_std_exx_accumul)
                   rho_kj = zero
                   do nsf1 = 1, kg%nsup                         
                      do nsf2 = 1, jb%nsup                                   
                         rho_kj(:,:,:,nsf1,nsf2) = &
                              Phy_k(:,:,:,nsf1) * phi_j(:,:,:,nsf2)  
                            !do r = 1, 2*extent+1
                            !   do s = 1, 2*extent+1
                            !      do t = 1, 2*extent+1                         
                            !         if(isnan(rho_kj(r,s,t,nsf1,nsf2))) write(*,*) 'rho: ',r,s,t,nsf1,nsf2

                            !         if(isnan(Phy_k(r,s,t,nsf1))) write(*,*) 'Phy: ',r,s,t,nsf1,nsf2
                            !         if(isnan(phi_j(r,s,t,nsf2))) write(*,*) 'phi: ',r,s,t,nsf1,nsf2
                            !      end do
                            !   end do
                            !end do
                      end do
                   end do
                   call stop_timer(tmr_std_exx_accumul,.true.)                 
                   !
                   do nsf1 = 1, kg%nsup
                      do nsf2 = 1, jb%nsup
                         
                         call start_timer(tmr_std_exx_poisson)    
                         work_in_3d  = zero
                         work_out_3d = zero
                         work_in_3d  = rho_kj(:,:,:,nsf1,nsf2)
                         
                         call exx_v_on_grid(inode,extent,work_in_3d,work_out_3d,r_int,   &
                              exx_psolver,p_scheme,pulay_radius,p_omega,p_ngauss,p_gauss,&
                              w_gauss)
                         
                         vhf_kj(:,:,:,nsf1,nsf2) = work_out_3d
                         
                         !do r = 1, 2*extent+1
                         !   do s = 1, 2*extent+1
                         !      do t = 1, 2*extent+1                         
                         !         if(isnan(vhf_kj(r,s,t,nsf1,nsf2))) write(*,*) 'vhf: ',r,s,t,nsf1,nsf2
                         !      end do
                         !   end do
                         !end do
                         
                         call stop_timer(tmr_std_exx_poisson,.true.)
                      end do
                   end do
                   !
                   Ome_kj = zero
                   !
                   call start_timer(tmr_std_exx_accumul)  
                   do nsf1 = 1, kg%nsup
                      do nsf2 = 1, jb%nsup                                       
                         Ome_kj(:,:,:,nsf1,nsf2) = &
                              vhf_kj(:,:,:,nsf1,nsf2) * phi_k(:,:,:,nsf1)                   
                      end do
                   end do
                   call stop_timer(tmr_std_exx_accumul,.true.)
                   !
                   if(ia%nsup/=ahalo%ndimi(i_in_prim)) write(24,*) 'Error1: ',ia%nsup,ahalo%ndimi(i_in_prim)
                   if(jb%nsup/=bndim2(nbkbeg+j-1))     write(24,*) 'Error2: ',jb%nsup,bndim2(nbkbeg+j-1)
                   !
                   call start_timer(tmr_std_exx_matmult)
                   do nsf2 = 1, jb%nsup                                         
                      !
                      ncaddr = ncbeg + ia%nsup * (nsf2 - 1)
                      !
                      do nsf1 = 1, ia%nsup
                         !
                         do nsf3 = 1, kg%nsup
                            !
                            exx_mat_elem = zero
                            !
                            do r = 1, 2*extent+1
                               do s = 1, 2*extent+1
                                  do t = 1, 2*extent+1                         
                                     
                                     exx_mat_elem = exx_mat_elem &                                    
                                          + phi_i(r,s,t,nsf1)    &
                                          * Ome_kj(r,s,t,nsf3,nsf2) * dv
                                     
                                  end do
                               end do
                            end do
                            !
                            c(ncaddr + nsf1 - 1) = c(ncaddr + nsf1 - 1) + exx_mat_elem
                            !
                         end do ! nsf3
                         !
                         !
                      end do ! nsf1
                      !
                      !
                   end do ! nsf2
                   call stop_timer(tmr_std_exx_matmult,.true.)
                   !
!!$
!!$
!                end if !( screen ij )
!!$
!!$
                end if ! ( ncbeg /=0 )
             end if ! ( j_in_halo /=0 )
             !
             !nbbeg = nbbeg + ia%nsup * jb%nsup ! nd3 * nd2
             !
!!$
!!$ ****[ j end loop ]****
!!$
          end do ! End of j = 1, mat(np,SXrange)%n_nab(ni)           
          !          
!!$
!!$ ****[ i end loop ]****
!!$
          !
       end do ! End of i = 1, at%n_hnab(k_in_halo)
       !
!!$
!!$ ****[ k end loop ]****
!!$
       !
    end do ! End of k = 1, ahalo%nh_part(kpart)
    !
    !
    return
  end subroutine m_kern_exx
  !
  !!***
end module exx_kernel_default
