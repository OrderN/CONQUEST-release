! -*- mode: F90; mode: font-lock -*-
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

  use exx_types,                 only: reckernel_3d, exx_debug
  use exx_io

  implicit none 

  ! Area identification
  integer, parameter, private :: area = 13

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
  !!   2020/12/06 16:26 Lionel
  !!    added calculation of ERIs with storage
  !!  SOURCE
  !!
  subroutine get_X_matrix( exxspin, scheme, backup_eris, niter, level )

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
    use matrix_data,    only: mat, Hrange, Srange, Xrange, SXrange
    use mult_module,    only: S_X_SX, mat_p, mult 
    use mult_module,    only: matX, matK, matrix_scale, matrix_trace
    use mult_module,    only: matrix_product_trace, matrix_product_trace_length 
    use mult_module,    only: return_matrix_value,  matrix_pos
    use mult_module,    only: store_matrix_value,   store_matrix_value_pos
    use memory_module,  only: write_mem_use, reg_alloc_mem, reg_dealloc_mem, type_dbl, type_int
    use species_module, only: nsf_species
    !
    use exx_memory,     only: exx_mem_alloc
    !
    use fft_interface_module, only: fft3_init_wrapper
    !
    use exx_types, only: prim_atomic_data, neigh_atomic_data, store_eris, &
         tmr_std_exx_evalpao, &
         tmr_std_exx_setup,   tmr_std_exx_barrier, &
         tmr_std_exx_fetch,   tmr_std_exx_accumul, &
         tmr_std_exx_matmult, tmr_std_exx_poisson, &
         tmr_std_exx_allocat, tmr_std_exx_dealloc, &
         tmr_std_exx_fetch,   tmr_std_exx_kernel,  &
         unit_timers_write,   unit_memory_write,   &
         unit_exx_debug,      unit_eri_debug,      &
         file_exx_timers,     file_exx_memory,     &
         file_exx_debug,      file_eri_debug
         
    use exx_types, only: grid_spacing, r_int, extent, &
         ewald_alpha, ewald_rho, ewald_pot,   &
         pulay_radius, p_omega, p_ngauss, p_gauss, w_gauss, &
         exx_psolver,p_scheme, kernel,    &
         exx_total_time,  eris, exx_filter
    !
    use exx_poisson, only: exx_scal_rho_3d, exx_ewald_rho, exx_ewald_pot
    !
    !**<lat>** ISF Poisson solver Will be available in the forthcoming version 
    !use Poisson_Solver, only: createBeylkin
    !
    implicit none

    integer, intent(in), optional :: level
    integer, intent(in)           :: exxspin, scheme
    logical, intent(in)           :: backup_eris
    
    integer, intent(in)           :: niter

    ! Local variables
    integer :: invdir,ierr,kpart,ind_part,ncover_yz,ipart,nnode
    integer :: icall,n_cont,k_off
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

    logical      :: get_exx, exist
    real(double) :: t0,t1
    integer      :: maxsuppfuncs, nb_eris

    !type(prim_atomic_data)  :: ia !i_alpha
    !type(neigh_atomic_data) :: jb !j_beta
    !type(neigh_atomic_data) :: kg !k_gamma
    !type(neigh_atomic_data) :: ld !l_delta
    !
    !type(store_eris), dimension(:), allocatable :: eris
    !
    real(double)      :: xyz_ghost(3), r_ghost
    !character(len=20) :: filename3, filename4, filename5, filename6

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
       !
       if ( niter == 1 ) then
          call io_assign(unit_timers_write)
          call get_file_name('exx_timers',numprocs,inode,file_exx_timers)
          open(unit_timers_write,file=file_exx_timers)
          !
          call io_assign(unit_memory_write)
          call get_file_name('exx_memory',numprocs,inode,file_exx_memory)
          open(unit_memory_write,file=file_exx_memory)
          !
       else
          !inquire(file=filename3, exist=exist)
          !if ( exist ) then
          open(unit_timers_write,file=file_exx_timers,status='old', position='append')
          open(unit_memory_write,file=file_exx_memory,status='old', position='append')
          !else
          !   open(unit_timers_write,file=filename3,status='new')
          !end if
       end if
       !
    end if

    if ( exx_debug ) then

       if ( niter == 1 ) then
          call io_assign(unit_exx_debug)
          call get_file_name('exx_debug',numprocs,inode,file_exx_debug)
          open(unit_exx_debug,file=file_exx_debug)
          !
       else
          open(unit_exx_debug,file=file_exx_debug,status='old', position='append')
          !
       end if

       call io_assign(unit_eri_debug)
       call get_file_name('eri_debug',numprocs,inode,file_eri_debug)
       !inquire(file=filename6, exist=exist)
       !if ( exist ) then
       !   open(unit_eri_debug,file=filename6,status='old', position='append')
       !else
       open(unit_eri_debug,file=file_eri_debug)
       !end if
       !
       !call exx_write_head(unit_exx_debug,inode,bundle%groups_on_node) 
       !call exx_global_write()
       !
    end if
    !
    maxsuppfuncs = maxval(nsf_species)
    !
    if ( scheme > 0 ) then
       call exx_mem_alloc(extent,0,0,'work_3d' ,'alloc')
    end if
    !
    !DRB! This appears to set up the different Poisson solvers
    if (exx_psolver == 'fftw')     then
       call exx_mem_alloc(extent,0,0,'fftw_3d','alloc')  
       call exx_mem_alloc(extent,0,0,'reckernel_3d','alloc')

       poisson_fftw: select case(p_scheme)          

       case('default')
          call exx_scal_rho_3d(inode,extent,r_int,p_scheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)

       case('ewald')
          call exx_mem_alloc(extent,0,0,'ewald_3d','alloc')
          call exx_ewald_rho(ewald_rho,extent,ewald_alpha,r_int)
          call exx_ewald_pot(ewald_pot,extent,ewald_alpha,r_int)          
          call exx_scal_rho_3d(inode,extent,r_int,p_scheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)

       case('pulay')
          call exx_scal_rho_3d(inode,extent,r_int,p_scheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)

       case('yukawa')
          call exx_scal_rho_3d(inode,extent,r_int,p_scheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)

       case('gauss')
          call cq_abort('EXX: Gaussian representation if 1/r for solving &
               &the Poisson equation &
               &is currently under testing...')
          !call createBeylkin(p_gauss,w_gauss,r_int)     
          !call exx_scal_rho_3d(inode,extent,r_int,p_scheme,pulay_radius, &
          !     p_omega,p_ngauss,p_gauss,w_gauss)

       case default
          call exx_scal_rho_3d(inode,extent,r_int,p_scheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)

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
    if ( scheme > 0 ) then
       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_i','alloc')       
       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_j','alloc')
       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_k','alloc')
       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_l','alloc')
       !
       if ( scheme == 1 ) then
          call exx_mem_alloc(extent,maxsuppfuncs,0,'Phy_k','alloc')
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'Ome_kj','alloc')       
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_kj','alloc')       
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_kj','alloc')
          !
       else if ( scheme == 2 ) then       
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_ki','alloc')       
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_lj','alloc')
          !
       else if ( scheme == 3 ) then       
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_ki','alloc')       
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_lj','alloc')
          !
          if ( niter == 1 ) then
             allocate( eris( mult(S_X_SX)%ahalo%np_in_halo ), STAT=stat)
             if(stat/=0) call cq_abort('Error allocating memory to eris/exx !',stat)
          end if
          !for now allocated here ; should not ; better to have this before
          ! SCF calculation ; note that here it is never deallocated!
          !eris_size = int( sqrt(dble(size( mat_p(matX(  exxspin  ))%matrix ))) )
          !print*, 'inode', inode, inode, size( mat_p(matX(  exxspin  ))%matrix), eris_size, &
          !     & mat_p(matX(  exxspin  ))%length, mat_p(matX(  exxspin  ))%sf1_type, mat_p(matX(  exxspin  ))%sf2_type
          !call exx_mem_alloc(eris_size**4,0,0,'eris','alloc')
          !
          !
       end if
       !
       !
    end if
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
    allocate(recv_part(0:mult(S_X_SX)%comms%inode),STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error allocating recv_part')
    call stop_timer(tmr_std_exx_allocat,.true.)
    !
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
       call hf_write_info(unit_exx_debug,inode,10,mult(S_X_SX)%ahalo%np_in_halo,0,0,0,0, &
            0,0,'XX',xyz_ghost,r_ghost,0,0,0,zero,zero,0)
    end if
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
          call hf_write_info(unit_exx_debug,inode,20,0,kpart,mult(S_X_SX)%ahalo%np_in_halo,0,0, &
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
             call start_timer(tmr_std_exx_dealloc)
             if(allocated(b_rem)) deallocate(b_rem)
             call stop_timer(tmr_std_exx_dealloc,.true.)
             !
             if(mult(S_X_SX)%parts%i_cc2node(ind_part)==myid+1) then
                lenb_rem = mult(S_X_SX)%bmat(ipart)%part_nd_nabs
             else
                lenb_rem = mult(S_X_SX)%comms%ilen3rec(ipart,nnode)
             end if
             !
             call start_timer(tmr_std_exx_allocat)
             allocate(b_rem(lenb_rem))
             call stop_timer(tmr_std_exx_allocat,.true.)   
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
          call start_timer(tmr_std_exx_dealloc)
          if(allocated(b_rem)) deallocate(b_rem)
          call stop_timer(tmr_std_exx_dealloc,.true.)
          !
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
       if ( scheme == 1 ) then

          if( myid==0) write(io_lun,*) 'EXX: performing CRI calculation on kpart =', kpart
          
          call m_kern_exx_cri( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
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

          !print*, ' mat X ', shape( mat_p(matX(  exxspin  ))%matrix)
          !print*, ' mat K ', shape( mat_p(matK(  exxspin  ))%matrix)

       else if ( scheme == 2 ) then

          if( myid==0) write(io_lun,*) 'EXX: performing on-the-fly ERI calculation on kpart =', kpart

          call m_kern_exx_eri( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
               ibpart_rem,ibseq_rem,ibndimj_rem, & 
                                !atrans, &
               b_rem,  &
               mat_p(matX(  exxspin  ))%matrix,      &
               mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
               mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
               mult(S_X_SX)%prim%mx_iprim, &
                                !lena,     &
               lenb_rem, &
               mat_p(matX(  exxspin  ))%length, backup_eris)

       else if (scheme == 3 ) then

          if ( niter == 1 ) then
             !
             get_exx = .false.
             !
             if( myid==0 ) write(io_lun,*) 'EXX: preparing store ERI calculation on kpart =', kpart
             !
             ! First dummy call to get the number of ERIs on each proc
             !
             call m_kern_exx_dummy( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
                  ibpart_rem,ibseq_rem,ibndimj_rem, & 
                  !atrans, &
                  b_rem,  &
                  mat_p(matX(  exxspin  ))%matrix,      &
                  mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
                  mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
                  mult(S_X_SX)%prim%mx_iprim, &
                  !lena,     &
                  lenb_rem, &
                  mat_p(matX(  exxspin  ))%length, nb_eris, get_exx, .false. )
             
             if( myid==0 ) write(io_lun,*) 'EXX: allocate ERIs on kpart =', kpart
             !
             ! Allocate ERI arrays
             !
             call start_timer(tmr_std_exx_allocat)
             allocate(eris(kpart)%store_eris( nb_eris ), STAT=stat)
             call stop_timer(tmr_std_exx_allocat,.true.)
             if(stat/=0) call cq_abort('Error allocating memory toeris/exx !',stat)
             eris(kpart)%store_eris = 0.0d0
             !
             call start_timer(tmr_std_exx_allocat)
             allocate(eris(kpart)%filter_eris( nb_eris ), STAT=stat)
             call stop_timer(tmr_std_exx_allocat,.true.)
             if(stat/=0) call cq_abort('Error allocating memory toeris/exx !',stat)
             call reg_alloc_mem(area_exx, nb_eris, type_int,'eris',unit_memory_write)
             eris(kpart)%filter_eris = .true.
             !
             ! Second dummy call for poor-man filtering of ERIs
             !
             if( myid==0 ) write(io_lun,*) 'EXX: setup filtering on kpart =', kpart
             !
             if ( exx_filter ) then
                !
                if( myid==0 ) write(io_lun,*) 'EXX: setup filtering on kpart =', kpart
                !
                call m_kern_exx_dummy( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
                     ibpart_rem,ibseq_rem,ibndimj_rem, & 
                     !atrans, &
                     b_rem,  &
                     mat_p(matX(  exxspin  ))%matrix,      &
                     mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
                     mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
                     mult(S_X_SX)%prim%mx_iprim, &
                     !lena,     &
                     lenb_rem, &
                     mat_p(matX(  exxspin  ))%length, nb_eris, get_exx, exx_filter )
             end if
             !
             ! should be a single call not embeded in the kpart loop... sorry for that
             call fft3_init_wrapper( 2*extent+1  )
             !
             if( myid==0 ) write(io_lun,*) 'EXX: compute and store ERIs on kpart =', kpart
             !
             ! Third call to compute and store ERIs
             !       
             call m_kern_exx_eri( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
                  ibpart_rem,ibseq_rem,ibndimj_rem, & 
                  !atrans, &
                  b_rem,  &
                  mat_p(matX(  exxspin  ))%matrix,      &
                  mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
                  mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
                  mult(S_X_SX)%prim%mx_iprim, &
                  !lena,     &
                  lenb_rem, &
                  mat_p(matX(  exxspin  ))%length, backup_eris)
             
          else

             get_exx = .true.

             if( myid==0 ) write(io_lun,*) 'EXX: use stored ERIs to get X on kpart =', kpart

             call m_kern_exx_dummy( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
                  ibpart_rem,ibseq_rem,ibndimj_rem, & 
                  !atrans, &
                  b_rem,  &
                  mat_p(matX(  exxspin  ))%matrix,      &
                  mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
                  mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
                  mult(S_X_SX)%prim%mx_iprim, &
                  !lena,     &
                  lenb_rem, &
                  mat_p(matX(  exxspin  ))%length, nb_eris, get_exx, .false. )
          end if
          
       else if ( scheme == -1 ) then

          get_exx = .false.
          
          if( myid==0 ) write(io_lun,*) 'EXX: dummy calculation on kpart =', kpart
          
          call m_kern_exx_dummy( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
               ibpart_rem,ibseq_rem,ibndimj_rem, & 
                                !atrans, &
               b_rem,  &
               mat_p(matX(  exxspin  ))%matrix,      &
               mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
               mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
               mult(S_X_SX)%prim%mx_iprim, &
                                !lena,     &
               lenb_rem, &
               mat_p(matX(  exxspin  ))%length, nb_eris, get_exx, .false. )

       end if

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
    call start_timer(tmr_std_exx_barrier)
    call my_barrier
    call stop_timer(tmr_std_exx_barrier,.true.)
    !
    call start_timer(tmr_std_exx_dealloc)
    deallocate(nreqs,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating nreqs')
    deallocate(recv_part,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating ibpart_rem')
    call stop_timer(tmr_std_exx_dealloc,.true.)
    !
    call start_timer(tmr_std_exx_barrier)
    call my_barrier
    call stop_timer(tmr_std_exx_barrier,.true.)
    !
    call start_timer(tmr_std_exx_dealloc)
    deallocate(ibpart_rem,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating recv_part')
    call stop_timer(tmr_std_exx_dealloc,.true.)
    !
    call start_timer(tmr_std_exx_barrier)
    call my_barrier
    call stop_timer(tmr_std_exx_barrier,.true.)
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
    !   print*, 'exx energy', -exx*1.d0/2.d0*0.25
    !end if
    !
    !if ( exx_debug ) then
    !call exx_write_tail(unit1,inode)  
    !end if
    !
    !
    if ( scheme > 0 ) then
       call exx_mem_alloc(extent,0,0,'work_3d' ,'dealloc')
       !
       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_i','dealloc')       
       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_j','dealloc')
       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_k','dealloc')    
       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_l','dealloc')
       !
       if ( scheme == 1 ) then
          call exx_mem_alloc(extent,maxsuppfuncs,0,'Phy_k','dealloc')
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'Ome_kj','dealloc')   
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_kj','dealloc')       
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_kj','dealloc')
          !
       else if( scheme == 2 ) then
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_ki','dealloc')       
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_lj','dealloc')
          !
       else if( scheme == 3 ) then
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_ki','dealloc')       
          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_lj','dealloc')
          !
       end if
       !
    end if
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
    !if (inode == ionode) then
    call write_mem_use(unit_memory_write,area_exx)
    !end if
    !
    call start_timer(tmr_std_exx_barrier)
    call my_barrier()
    call stop_timer(tmr_std_exx_barrier,.true.)
    !
    call stop_timer(tmr_std_exx,.true.)
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
         tmr_std_exx_barrier%t_tot + &
         tmr_std_exx_fetch%t_tot  

    if ( iprint_exx > 3 ) then
       
       call print_timer(tmr_std_exx_setup,  "exx_setup   time:", unit_timers_write)    
       call print_timer(tmr_std_exx_write,  "exx_write   time:", unit_timers_write)    
       !call print_timer(tmr_std_exx_kernel, "exx_kernel  time:", unit_timers_write)    
       call print_timer(tmr_std_exx_fetch,  "exx_fetch   time:", unit_timers_write)    
       call print_timer(tmr_std_exx_evalpao,"exx_evalpao time:", unit_timers_write)    
       call print_timer(tmr_std_exx_poisson,"exx_poisson time:", unit_timers_write)
       call print_timer(tmr_std_exx_matmult,"exx_matmult time:", unit_timers_write)    
       call print_timer(tmr_std_exx_accumul,"exx_accumul time:", unit_timers_write)    
       call print_timer(tmr_std_exx_allocat,"exx_allocat time:", unit_timers_write)    
       call print_timer(tmr_std_exx_dealloc,"exx_dealloc time:", unit_timers_write)    
       call print_timer(tmr_std_exx_barrier,"exx_barrier time:", unit_timers_write)          
       call print_timer(tmr_std_exx_kernel, "exx_kernel  time:", unit_timers_write)    
       call print_timer(tmr_std_exx,        "exx_total   time:", unit_timers_write)    
       write(unit=unit_timers_write,fmt='("Timing: Proc ",i6,": Time spent in ", a50, " = ", &
            &f12.5," s")') inode, 'get_X_matrix',  tmr_std_exx%t_tot
       write(unit=unit_timers_write,fmt='("Timing: Proc ",i6,": Time spent in ", a50, " = ", &
            &f12.5," s")') inode, 'timer calls',  tmr_std_exx%t_tot-exx_total_time
       write(unit_timers_write,*)
       !call io_close(unit_matrix_write)

       !call io_close(unit_output_write)
       !call io_close(unit_screen_write)
       !
       !call io_close(unit_eri_debug)
       !call io_close(unit_exx_debug)
       call io_close(unit_memory_write)
       call io_close(unit_timers_write)
       !
       !
    end if

    if ( exx_debug ) then
       call io_close(unit_eri_debug)
       call io_close(unit_exx_debug)       
    end if
    
    !****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_X_matrix',echo=.true.)
    !****lat>$

    return
  end subroutine get_X_matrix
  !
  !
  !
  !!****f* exx_kernel_default/m_kern_exx_cri *
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
  subroutine m_kern_exx_cri(k_off, kpart, ib_nd_acc, ibaddr, nbnab, &
       ibpart, ibseq, bndim2, b, c, ahalo, chalo, at, mx_absb, mx_part,  &
       mx_iprim, lenb, lenc )

    use numbers,        only: zero, one
    use matrix_module,  only: matrix_halo, matrix_trans
    use global_module,  only: area_exx
    !
    use basic_types,    only: primary_set
    use primary_module, only: bundle 
    use matrix_data,    only: Hrange, SXrange, Xrange, Srange
    use mult_module,    only: return_matrix_value, S_X_SX
    use cover_module,   only: BCS_parts
    !
    use species_module, only: nsf_species
    !
    use exx_evalpao,    only: exx_phi_on_grid
    !
    use exx_types, only: prim_atomic_data, neigh_atomic_data, &
         tmr_std_exx_accumul,   &
         tmr_std_exx_matmult, tmr_std_exx_poisson, &
                                !unit_matrix_write, unit_output_write,     & 
                                !unit_timers_write, unit_screen_write,     &
                                !unit_memory_write, &
         grid_spacing, r_int, extent, &
         !ewald_alpha, ewald_charge, ewald_rho, ewald_pot,   &
         pulay_radius, p_omega, p_ngauss, p_gauss, w_gauss, &
         exx_psolver,p_scheme, &         
         unit_exx_debug
    !
    use exx_types, only: phi_i, phi_j, phi_k, phi_l, &
         Phy_k, rho_kj, Ome_kj, vhf_kj, &
         work_in_3d, work_out_3d, fftwrho3d

    use exx_poisson, only: exx_v_on_grid
    !
    use exx_module,  only: get_halodat, get_iprimdat
    !
    implicit none
    !
    ! Passed variables
    type(matrix_halo)  :: ahalo, chalo
    type(matrix_trans) :: at
    integer            :: mx_absb, mx_part, mx_iprim, lenb, lenc
    integer            :: kpart, k_off
    real(double)       :: b(lenb)
    real(double)       :: c(lenc)
    !integer, optional  :: debug
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
    integer :: i, i_in_prim, icad, nbbeg, j_in_halo, ncbeg
    integer :: nb_nd_kbeg
    integer :: nd1, nd3
    integer :: nbaddr, ncaddr
    integer :: lbnab2ch(mx_absb)  ! Automatic array
    integer :: l, lseq, lpart
    integer :: np, ni
    !
    real(double), dimension(3) :: xyz_zero  = zero
    !real(double), dimension(3) :: xyz_ij    = zero
    !real(double), dimension(3) :: xyz_kl    = zero
    real(double)               ::   dr,dv,K_val
    real(double)               ::   exx_mat_elem
    !real(double)               ::   screen_ij, range_ij
    !real(double)               ::   screen_kl, range_kl
    !
    type(prim_atomic_data)  :: ia !i_alpha
    type(neigh_atomic_data) :: jb !j_beta
    type(neigh_atomic_data) :: kg !k_gamma
    type(neigh_atomic_data) :: ld !l_delta
    !
    integer                 :: maxsuppfuncs, nsf1, nsf2, nsf3
    integer                 :: r, s, t
    !
    !
    dr = grid_spacing
    dv = dr**3
    !ewald_alpha  = 0.5
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
            ahalo%lab_hcell(kpart),'k',.true.,unit_exx_debug)
       !
       !print*, 'k',k, 'global_num',kg%global_num,'spe',kg%spec, kg%xyz
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
               'l',.true.,unit_exx_debug)
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
                K_val = b(nbaddr+nsf1-1)
                !
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
          call get_iprimdat(ia,kg,ni,i_in_prim,np,.true.,unit_exx_debug)          
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
                   call get_halodat(jb,kg,jseq,chalo%i_hbeg(jpart),         &
                        BCS_parts%lab_cell(BCS_parts%inv_lab_cover(jpart)), &
                        'j',.true.,unit_exx_debug)
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
                              w_gauss,fftwrho3d,reckernel_3d)

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
  end subroutine m_kern_exx_cri
  !

  !!****f* exx_kernel_default/m_kern_exx_eri *
  !!
  !!  NAME
  !!   m_kern_exx_eri
  !!
  !!  PURPOSE
  !! 
  !!  INPUTS
  !!
  !!  AUTHOR
  !!   L.A.Truflandier/D.R.Bowler
  !!
  !!  CREATION DATE
  !!   2020/10/27
  !!
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine m_kern_exx_eri(k_off, kpart, ib_nd_acc, ibaddr, nbnab, &
       ibpart, ibseq, bndim2, b, c, ahalo, chalo, & 
       at, mx_absb, mx_part, mx_iprim, lenb, lenc, backup_eris )

    use numbers,        only: zero, one
    use matrix_module,  only: matrix_halo, matrix_trans
    use global_module,  only: area_exx
    !
    use basic_types,    only: primary_set
    use primary_module, only: bundle 
    use matrix_data,    only: Hrange, SXrange, Xrange, Srange
    use mult_module,    only: return_matrix_value, S_X_SX
    use cover_module,   only: BCS_parts
    !
    use species_module, only: nsf_species
    !
    use exx_evalpao,    only: exx_phi_on_grid
    !
    use exx_types, only: prim_atomic_data, neigh_atomic_data, &
         tmr_std_exx_accumul, tmr_std_exx_poisson,  &
         tmr_std_exx_poisson, &
         grid_spacing, r_int, extent, eris,                 &
                                !ewald_alpha, ewald_charge, ewald_rho, ewald_pot,   &
         pulay_radius, p_omega, p_ngauss, p_gauss, w_gauss, &
         exx_psolver,p_scheme, &
         unit_exx_debug, unit_eri_debug
    !
    use exx_types,  only: phi_i, phi_j, phi_k, phi_l, eris, &
         rho_ki, vhf_lj, work_in_3d, work_out_3d, fftwrho3d

    use exx_module, only: get_halodat, get_iprimdat
    !
    use exx_poisson,only: exx_v_on_grid
    !
    implicit none
    !
    ! Passed variables
    type(matrix_halo),  intent(in)    :: ahalo, chalo
    type(matrix_trans), intent(in)    :: at
    integer,            intent(in)    :: mx_absb, mx_part, mx_iprim, lenb, lenc
    integer,            intent(in)    :: kpart, k_off
    real(double),       intent(in)    :: b(lenb)
    real(double),       intent(inout) :: c(lenc)
    logical,            intent(in)    :: backup_eris
    !integer, optional  :: debug
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
    integer :: i, i_in_prim, icad, nbbeg, j_in_halo, ncbeg
    integer :: nb_nd_kbeg
    integer :: nd1, nd3
    integer :: nbaddr, ncaddr
    integer :: lbnab2ch(mx_absb)  ! Automatic array
    integer :: l, lseq, lpart
    integer :: np, ni
    !
    real(double), dimension(3) :: xyz_zero  = zero
    !real(double), dimension(3) :: xyz_ij    = zero
    !real(double), dimension(3) :: xyz_kl    = zero
    !
    real(double)               ::   dr,dv,K_val
    real(double)               ::   exx_mat_elem
    !real(double)               ::   screen_ij, range_ij
    !real(double)               ::   screen_kl, range_kl
    !
    type(prim_atomic_data)  :: ia !i_alpha
    type(neigh_atomic_data) :: jb !j_beta
    type(neigh_atomic_data) :: kg !k_gamma
    type(neigh_atomic_data) :: ld !l_delta
    !
    integer                 :: maxsuppfuncs
    integer                 :: nsf_kg, nsf_ld, nsf_ia, nsf_jb
    integer                 :: r, s, t, count
    !
    !
    dr = grid_spacing
    dv = dr**3
    maxsuppfuncs = maxval(nsf_species)
    !
    count = 1
    !
!!$
!!$ ****[ k loop ]****
!!$
    k_loop: do k = 1, ahalo%nh_part(kpart)
       !
       k_in_halo  = ahalo%j_beg(kpart) + k - 1
       k_in_part  = ahalo%j_seq(k_in_halo)
       nbkbeg     = ibaddr     (k_in_part) 
       nb_nd_kbeg = ib_nd_acc  (k_in_part)
       nd3        = ahalo%ndimj(k_in_halo)
       call get_halodat(kg,kg,k_in_part,ahalo%i_hbeg(ahalo%lab_hcover(kpart)), &
            ahalo%lab_hcell(kpart),'k',.true.,unit_exx_debug)
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
       !Phy_k = zero
       !
!!$
!!$ ****[ l do loop ]****
!!$
       l_loop: do l = 1, nbnab(k_in_part)
          !l_in_halo = lbnab2ch(l)                
          lpart = ibpart(nbkbeg+l-1) + k_off
          lseq  = ibseq (nbkbeg+l-1)
          call get_halodat(ld,kg,lseq,chalo%i_hbeg(lpart),         &
               BCS_parts%lab_cell(BCS_parts%inv_lab_cover(lpart)), &
               'l',.true.,unit_exx_debug)
          !
          !write(*,*) 'l',l, 'global_num',ld%global_num,'spe',ld%spec
          !
          call exx_phi_on_grid(inode,ld%global_num,ld%spec,extent,     &
               ld%xyz,maxsuppfuncs,phi_l,r_int,xyz_zero)             
          !
          ld_loop: do nsf_ld = 1, ld%nsup
             !
             nbaddr = nbbeg + kg%nsup * (nsf_ld - 1)
             !
             kg_loop: do nsf_kg = 1, kg%nsup                         
                !
                if ( backup_eris ) then
                   K_val = real(1,double)
                else
                   K_val = b(nbaddr+nsf_kg-1)
                end if
!!$
!!$ ****[ i loop ]****
!!$
                i_loop: do i = 1, at%n_hnab(k_in_halo)
                   i_in_prim = at%i_prim(at%i_beg(k_in_halo)+i-1)
                   nd1 = ahalo%ndimi      (i_in_prim)
                   ni  = bundle%iprim_seq (i_in_prim)
                   np  = bundle%iprim_part(i_in_prim)
                   icad  = (i_in_prim - 1) * chalo%ni_in_halo !***
                   !nbbeg = nb_nd_kbeg
                   !
                   !print*, 'i_in_prim', i_in_prim, 'nd1', nd1, 'ni', ni, 'np', np
                   !
                   call get_iprimdat(ia,kg,ni,i_in_prim,np,.true.,unit_exx_debug)          
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
                   j_loop: do j = 1, nbnab(k_in_part)!mat(np,Xrange)%n_nab(ni)                   
                      !nbbeg     = nb_nd_kbeg
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
                            !
                            call get_halodat(jb,kg,jseq,chalo%i_hbeg(jpart),         &
                                 BCS_parts%lab_cell(BCS_parts%inv_lab_cover(jpart)), &
                                 'j',.true.,unit_exx_debug)
                            !
                            call exx_phi_on_grid(inode,jb%global_num,jb%spec,extent, &
                                 jb%xyz,maxsuppfuncs,phi_j,r_int,xyz_zero)
                            !
                            jb_loop: do nsf_jb = 1, jb%nsup                                         
                               !
                               ncaddr = ncbeg + ia%nsup * (nsf_jb - 1)
                               !
                               call start_timer(tmr_std_exx_poisson) 
                               work_out_3d = zero
                               work_in_3d  = phi_l(:,:,:,nsf_ld)*phi_j(:,:,:,nsf_jb)
                               !
                               !if ( eris(kpart)%filter_eris( count ) ) then
                               call exx_v_on_grid(inode,extent,work_in_3d,work_out_3d,r_int,   &
                                    exx_psolver,p_scheme,pulay_radius,p_omega,p_ngauss,p_gauss,&
                                    w_gauss,fftwrho3d,reckernel_3d)
                               !end if
                               !
                               vhf_lj(:,:,:,nsf_ld,nsf_jb) = work_out_3d
                               !
                               call stop_timer(tmr_std_exx_poisson,.true.)

                               ia_loop: do nsf_ia = 1, ia%nsup
                                  !
                                  exx_mat_elem = zero
                                  !
                                  !if ( eris(kpart)%filter_eris( count ) ) then
                                     !
                                     call start_timer(tmr_std_exx_accumul)
                                     !
                                     rho_ki(:,:,:,nsf_kg,nsf_ia) = phi_k(:,:,:,nsf_kg)*phi_i(:,:,:,nsf_ia)
                                     !
                                     do r = 1, 2*extent+1
                                        do s = 1, 2*extent+1
                                           do t = 1, 2*extent+1                         

                                              exx_mat_elem = exx_mat_elem &                                    
                                                   + rho_ki(r,s,t,nsf_kg,nsf_ia) * K_val   &
                                                   * vhf_lj(r,s,t,nsf_ld,nsf_jb) * dv

                                           end do
                                        end do
                                     end do
                                     !
                                     call stop_timer(tmr_std_exx_accumul,.true.)
                                     !
                                  !end if
                                  !
                                  if ( exx_debug ) then
                                     write(unit_eri_debug,'(I3,A,I2,A,I3,A,I2,A,I3,A,I2,A,I3,A,I2,A,A,F20.16,5F16.12)') &
                                          kg%global_num,'{',nsf_kg,'}', ld%global_num,'{',nsf_ld,'}',&
                                          ia%ip,'{',nsf_ia,'}', jb%global_num,'{',nsf_jb,'}', '!', exx_mat_elem, &
                                          K_val, kg%r, ld%r, ia%r, jb%r
                                  end if
                                  !
                                  if ( backup_eris ) then
                                     !
                                     eris(kpart)%store_eris( count ) = exx_mat_elem
                                     !
                                  else                                     
                                     c(ncaddr + nsf_ia - 1) = c(ncaddr + nsf_ia - 1) + exx_mat_elem
                                     !
                                  end if
                                  !
                                  count = count + 1
                                  !
                               end do ia_loop ! nsf1
                               !
                               !
                            end do jb_loop ! nsf2
                            !
                         end if
                         !
                      end if

                   end do j_loop

                end do i_loop
!!$
!!$ ****[ j end loop ]****
!!$

             end do kg_loop  ! End of l = 1, nbnab(k_in_part)
!!$
!!$ ****[ i end loop ]****
!!$
          end do ld_loop

          !
          nbbeg = nbbeg + ld%nsup*kg%nsup !nd3 * nd2                
          !
!!$
!!$ ****[ l end loop ]****
!!$         
       end do l_loop
!!$
!!$ ****[ k end loop ]****
!!$
    end do k_loop
    !
    return
  end subroutine m_kern_exx_eri
  !
  !!****f* exx_kernel_default/m_kern_exx_dummy *
  !!
  !!  NAME
  !!   m_kern_exx_dummy
  !!
  !!  PURPOSE
  !! 
  !!  INPUTS
  !!
  !!  AUTHOR
  !!   L.A.Truflandier/D.R.Bowler
  !!
  !!  CREATION DATE
  !!   2020/10/27
  !!
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine m_kern_exx_dummy(k_off, kpart, ib_nd_acc, ibaddr, nbnab, &
       ibpart, ibseq, bndim2, b, c, ahalo, chalo,  & 
       at, mx_absb, mx_part, mx_iprim, lenb, lenc, &
       count_eris, compute_exx, filter_eris)

    use numbers,        only: zero, one
    use matrix_module,  only: matrix_halo, matrix_trans
    use global_module,  only: area_exx
    !
    use fft_interface_module, only: fft3_init_wrapper
    
    use basic_types,    only: primary_set
    use primary_module, only: bundle 
    use matrix_data,    only: Hrange, SXrange, Xrange, Srange
    use mult_module,    only: matK, return_matrix_value, S_X_SX
    use cover_module,   only: BCS_parts
    !
    !use group_module,   only: parts 
    !use pao_format,     only: pao
    !
    use species_module, only: nsf_species
    !
    use angular_coeff_routines, only: calc_mat_elem_gen
    !    
    use exx_evalpao, only: exx_phi_on_grid
    use exx_types,   only: reckernel_3d_filter, fftwrho3d_filter
    use exx_types,   only: prim_atomic_data, neigh_atomic_data, &         
         grid_spacing, r_int, unit_eri_debug, unit_exx_debug, &
         pulay_radius, p_omega, p_ngauss, p_gauss, w_gauss, &
         exx_psolver, p_scheme
    !
    use exx_types,  only: eris, exx_filter_thr, exx_filter_extent
    use exx_module, only: get_halodat, get_iprimdat
    !
    use exx_poisson,only: exx_v_on_grid, exx_scal_rho_3d
    !
    implicit none
    !
    ! Passed variables
    type(matrix_halo)  :: ahalo, chalo
    type(matrix_trans) :: at
    !
    integer            :: mx_absb, mx_part, mx_iprim, lenb, lenc
    integer            :: kpart, k_off, count_eris
    real(double)       :: b(lenb)
    real(double)       :: c(lenc)
    !
    logical            :: compute_exx
    logical            :: filter_eris
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
    integer :: i, i_in_prim, icad, nbbeg, j_in_halo, ncbeg
    integer :: nb_nd_kbeg
    integer :: nd1, nd3
    integer :: nbaddr, ncaddr
    integer :: lbnab2ch(mx_absb)  ! Automatic array
    integer :: l, lseq, lpart
    integer :: np, ni
    !
    real(double), dimension(3) :: xyz_zero  = zero
    !real(double), dimension(3) :: xyz_ij    = zero
    !real(double), dimension(3) :: xyz_kl    = zero
    !
    real(double)               ::   dr, dv, K_val, S_val
    real(double)               ::   exx_mat_elem
    !real(double)               ::   screen_ij, range_ij
    !real(double)               ::   screen_kl, range_kl
    !
    type(prim_atomic_data)  :: ia !i_alpha
    type(neigh_atomic_data) :: jb !j_beta
    type(neigh_atomic_data) :: kg !k_gamma
    type(neigh_atomic_data) :: ld !l_delta
    !
    integer                 :: maxsuppfuncs
    integer                 :: nsf_kg, nsf_ld, nsf_ia, nsf_jb
    integer                 :: r, s, t, stat
    !integer                 :: count_kg, count_ld
    !
    real(double), dimension(:,:,:,:),   allocatable :: phi_i_filter,  phi_j_filter
    real(double), dimension(:,:,:,:),   allocatable :: phi_k_filter,  phi_l_filter
    real(double), dimension(:,:,:,:,:), allocatable :: rho_ki_filter, vhf_lj_filter
    real(double), dimension(:,:,:),     allocatable :: work_in, work_out
    !
    dr = grid_spacing
    dv = dr**3
    !
    count_eris = 1
    !    
    maxsuppfuncs = maxval(nsf_species)    
    !
    if ( filter_eris ) then
       allocate( phi_i_filter(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1,maxsuppfuncs), STAT=stat )
       if(stat/=0) call cq_abort('Error allocating memory to phi_i_filter/exx !',stat)
       !
       allocate( phi_j_filter(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1,maxsuppfuncs), STAT=stat )
       if(stat/=0) call cq_abort('Error allocating memory to phi_j_filter/exx !',stat)
       !
       allocate( phi_k_filter(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1,maxsuppfuncs), STAT=stat )
       if(stat/=0) call cq_abort('Error allocating memory to phi_k_filter/exx !',stat)
       !
       allocate( phi_l_filter(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1,maxsuppfuncs), STAT=stat )
       if(stat/=0) call cq_abort('Error allocating memory to phi_l_filter/exx !',stat)
       !
       allocate( rho_ki_filter(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1,&
            maxsuppfuncs,maxsuppfuncs), STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to rho_ki_filter/exx !',stat)
       !
       allocate( vhf_lj_filter(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1,&
            maxsuppfuncs,maxsuppfuncs), STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to vhf_lj_filter/exx !',stat)
       !
       allocate( work_in(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1), STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to work_in_filter/exx !',stat)
       !
       allocate(work_out(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1), STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to work_out_filter/exx !',stat)
       !
       allocate(fftwrho3d_filter%arrayin(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1), STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to fftwin3d_filter/exx !',stat)
       !
       allocate(fftwrho3d_filter%arrayout(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1), STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to fftwout3d_filter/exx !',stat)
       !
       call fft3_init_wrapper( 2*exx_filter_extent+1  )
       !
       allocate(reckernel_3d_filter(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1), STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to reckernel_3d_filter/exx !',stat)

       call exx_scal_rho_3d(inode,exx_filter_extent,r_int,p_scheme,pulay_radius, &
            p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d_filter)
       
    end if

!!$
!!$ ****[ k loop ]****
!!$
    k_loop: do k = 1, ahalo%nh_part(kpart)
       !
       k_in_halo  = ahalo%j_beg(kpart) + k - 1
       k_in_part  = ahalo%j_seq(k_in_halo)
       nbkbeg     = ibaddr     (k_in_part) 
       nb_nd_kbeg = ib_nd_acc  (k_in_part)
       nd3        = ahalo%ndimj(k_in_halo)
       call get_halodat(kg,kg,k_in_part,ahalo%i_hbeg(ahalo%lab_hcover(kpart)), &
            ahalo%lab_hcell(kpart),'k',.true.,unit_exx_debug)
       !
       !print*, 'k',k, 'global_num',kg%global_num,'spe',kg%spec
       !
       if ( filter_eris ) then
          call exx_phi_on_grid(inode,kg%global_num,kg%spec,exx_filter_extent, &
               xyz_zero,maxsuppfuncs,phi_k_filter,r_int,xyz_zero)
          !print*, maxval(abs(phi_k_filter))          
       end if
       
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
!!$
!!$ ****[ l do loop ]****
!!$
       l_loop: do l = 1, nbnab(k_in_part)
          !l_in_halo = lbnab2ch(l)                
          lpart = ibpart(nbkbeg+l-1) + k_off
          lseq  = ibseq (nbkbeg+l-1)
          call get_halodat(ld,kg,lseq,chalo%i_hbeg(lpart),         &
               BCS_parts%lab_cell(BCS_parts%inv_lab_cover(lpart)), &
               'l',.true.,unit_exx_debug)
          !
          !write(*,*) 'l',l, 'global_num',ld%global_num,'spe',ld%spec
          !
          if ( filter_eris ) then             
             call exx_phi_on_grid(inode,ld%global_num,ld%spec,exx_filter_extent, &
                  ld%xyz,maxsuppfuncs,phi_l_filter,r_int,xyz_zero)
             !print*, maxval(abs(phi_l_filter))
 
          end if
          !
          ! count_ld = 0
          ld_loop: do nsf_ld = 1, ld%nsup
             !
             nbaddr = nbbeg + kg%nsup * (nsf_ld - 1)
             !
             !count_ld = count_ld + 1
             !
             !count_kg = 0
             kg_loop: do nsf_kg = 1, kg%nsup

                !count_kg = count_kg + 1
                !
                !if (inode == ionode) write(io_lun,*) 'k',k, 'global_num',kg%global_num,'spe',kg%spec, kg%l1(count_kg),  &
                !     kg%acz1(count_kg), kg%m1(count_kg)
                !
                !if (inode == ionode) write(io_lun,*) 'l',l, 'global_num',ld%global_num,'spe',ld%spec, ld%l1(count_ld),  &
                !     ld%acz1(count_ld), ld%m1(count_ld)
                !
                !if (inode == ionode) write(io_lun,*)
                !
                !call calc_mat_elem_gen(1,kg%spec,kg%l1(count_kg),kg%acz1(count_kg),kg%m1(count_kg),ld%spec,ld%l1(count_ld),&
                !     ld%acz1(count_ld),ld%m1(count_ld),ld%xyz(1),ld%xyz(2),&
                !     ld%xyz(3),S_val)
                !
                !
                K_val = b(nbaddr+nsf_kg-1)
                !
!!$
!!$ ****[ i loop ]****
!!$
                i_loop: do i = 1, at%n_hnab(k_in_halo)
                   i_in_prim = at%i_prim(at%i_beg(k_in_halo)+i-1)
                   nd1 = ahalo%ndimi      (i_in_prim)
                   ni  = bundle%iprim_seq (i_in_prim)
                   np  = bundle%iprim_part(i_in_prim)
                   icad  = (i_in_prim - 1) * chalo%ni_in_halo !***
                   !nbbeg = nb_nd_kbeg
                   !
                   !print*, 'i_in_prim', i_in_prim, 'nd1', nd1, 'ni', ni, 'np', np
                   !
                   call get_iprimdat(ia,kg,ni,i_in_prim,np,.true.,unit_exx_debug)          
                   !
                   !print*, 'i',i, 'global_num',ia%ip,'spe',ia%spec
                   !
                   if ( filter_eris ) then                      
                      call exx_phi_on_grid(inode,ia%ip,ia%spec,exx_filter_extent, &
                           ia%xyz,maxsuppfuncs,phi_i_filter,r_int,xyz_zero)
                      !print*,'phi_i_filter', maxval(abs(phi_i_filter))
                   end if
                   ! 
!!$
!!$ ****[ j loop ]****
!!$
                   j_loop: do j = 1, nbnab(k_in_part)!mat(np,Xrange)%n_nab(ni)                   
                      !nbbeg     = nb_nd_kbeg
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
                            !
                            call get_halodat(jb,kg,jseq,chalo%i_hbeg(jpart),         &
                                 BCS_parts%lab_cell(BCS_parts%inv_lab_cover(jpart)), &
                                 'j',.true.,unit_exx_debug)
                            !
                            if ( filter_eris ) then                               
                               call exx_phi_on_grid(inode,jb%global_num,jb%spec,exx_filter_extent, &
                                    jb%xyz,maxsuppfuncs,phi_j_filter,r_int,xyz_zero)
                               !print*, 'phi_j_filter', maxval(abs(phi_j_filter))
                            end if
                            !
                            !count_ld = 0
                            !
                            jb_loop: do nsf_jb = 1, jb%nsup                                         
                               !
                               !count_ld = count_ld + 1
                               
                               ncaddr = ncbeg + ia%nsup * (nsf_jb - 1)
                               !
                               !if (inode == ionode) write(io_lun,*) 'k',k, 'global_num',kg%global_num,'spe',kg%spec, &
                               !     kg%l1(count_kg),  &
                               !     kg%acz1(count_kg), kg%m1(count_kg)
                               !
                               !if (inode == ionode) write(io_lun,*) 'l',j, 'global_num',jb%global_num,'spe',jb%spec, &
                               !     jb%l1(count_ld),  &
                               !     jb%acz1(count_ld), jb%m1(count_ld)
                               !                               
                               !call calc_mat_elem_gen(1,kg%spec,kg%l1(count_kg),kg%acz1(count_kg),kg%m1(count_kg),&
                               !     jb%spec,jb%l1(count_ld),&
                               !     jb%acz1(count_ld),jb%m1(count_ld),jb%xyz(1),jb%xyz(2),&
                               !     jb%xyz(3),S_val)
                               !
                               if ( filter_eris ) then
                                  work_out = zero
                                  work_in  = phi_l_filter(:,:,:,nsf_ld)*phi_j_filter(:,:,:,nsf_jb)
                                  !
                                  !print*, 'work_in', maxval(abs(work_in))
                                  !                                  
                                  call exx_v_on_grid(inode,exx_filter_extent,work_in,work_out,r_int, &
                                       exx_psolver,p_scheme,pulay_radius,p_omega,p_ngauss,p_gauss,   &
                                       w_gauss,fftwrho3d_filter,reckernel_3d_filter)
                                  !
                                  !print*, 'work_out', maxval(abs(work_out))
                                  !
                                  vhf_lj_filter(:,:,:,nsf_ld,nsf_jb) = work_out
                                  !
                               end if
                               !
                               ia_loop: do nsf_ia = 1, ia%nsup
                                  !
                                  if ( compute_exx ) then
                                     !print*, 'kpart =', kpart, 'count_eris =', count_eris
                                     c(ncaddr + nsf_ia - 1) = c(ncaddr + nsf_ia - 1) + &
                                          eris(kpart)%store_eris(count_eris) * K_val
                                                                          
                                     !write(unit_eri_debug,'(I3,A,I2,A,I3,A,I2,A,I3,A,I2,A,I3,A,I2,A,A,6F16.12)') &
                                     !     kg%global_num,'{',nsf_kg,'}', ld%global_num,'{',nsf_ld,'}',&
                                     !     ia%ip,'{',nsf_ia,'}', jb%global_num,'{',nsf_jb,'}', '!',   &
                                     !     eris(kpart)%store_eris(count_eris), K_val,                 &
                                     !     kg%r, ld%r, ia%r, jb%r
                                     
                                  else
                                     !
                                     if ( filter_eris ) then
                                        !
                                        rho_ki_filter(:,:,:,nsf_kg,nsf_ia) = phi_k_filter(:,:,:,nsf_kg) &
                                             *phi_i_filter(:,:,:,nsf_ia)
                                        !
                                        !print*, 'vhf_lj_filter', maxval(abs(rho_ki_filter))
                                        !print*, 'rho_ki_filter', maxval(abs(rho_ki_filter))
                                        !print*, 'shape(vhf_lj_filter)', shape(vhf_lj_filter)
                                        !print*, 'shape(rho_ki_filter)', shape(rho_ki_filter)
                                        !
                                        !
                                        exx_mat_elem = zero
                                        !
                                        do r = 1, 2*exx_filter_extent+1
                                           do s = 1, 2*exx_filter_extent+1
                                              do t = 1, 2*exx_filter_extent+1                         
                                                 
                                                 exx_mat_elem = exx_mat_elem &                                    
                                                      + rho_ki_filter(r,s,t,nsf_kg,nsf_ia)    &
                                                      * vhf_lj_filter(r,s,t,nsf_ld,nsf_jb) * dv
                                              end do
                                           end do
                                        end do
                                        !
                                        if ( abs(exx_mat_elem) < exx_filter_thr ) then
                                           !
                                           eris(kpart)%filter_eris( count_eris ) = .false.
                                           !
                                        end if
                                        !                                                                        
                                        write(unit_eri_debug,'(I3,A,I2,A,I3,A,I2,A,I3,A,I2,A,I3,A,I2,A,A,F20.16,X,L,5F16.12)') &
                                             kg%global_num,'{',nsf_kg,'}', ld%global_num,'{',nsf_ld,'}',&
                                             ia%ip,'{',nsf_ia,'}', jb%global_num,'{',nsf_jb,'}', '!',   &
                                             exx_mat_elem, eris(kpart)%filter_eris(count_eris), &
                                             K_val, kg%r, ld%r, ia%r, jb%r
                                        !
                                     end if
                                     !
                                     c(ncaddr + nsf_ia - 1) = c(ncaddr + nsf_ia - 1) + 0.0d0
                                     !
                                  end if
                                  !
                                  count_eris = count_eris + 1
                                  !
                               end do ia_loop ! nsf1
                               !
                               !
                            end do jb_loop ! nsf2
                            !
                         end if
                         !
                      end if

                   end do j_loop

                end do i_loop
!!$
!!$ ****[ j end loop ]****
!!$

             end do kg_loop  ! End of l = 1, nbnab(k_in_part)
!!$
!!$ ****[ i end loop ]****
!!$
          end do ld_loop

          !
          nbbeg = nbbeg + ld%nsup*kg%nsup !nd3 * nd2                
          !
!!$
!!$ ****[ l end loop ]****
!!$         
       end do l_loop
!!$
!!$ ****[ k end loop ]****
!!$
    end do k_loop
    !
    print*, 'proc:', inode,'kpart:', kpart, 'nb. of ERIs:', count_eris
    !

  if ( filter_eris ) then
       deallocate( phi_i_filter, STAT=stat )
       if(stat/=0) call cq_abort('Error allocating memory to phi_i_filter/exx !',stat)
       !
       deallocate( phi_j_filter, STAT=stat )
       if(stat/=0) call cq_abort('Error allocating memory to phi_j_filter/exx !',stat)
       !
       deallocate( phi_k_filter, STAT=stat )
       if(stat/=0) call cq_abort('Error allocating memory to phi_k_filter/exx !',stat)
       !
       deallocate( phi_l_filter, STAT=stat )
       if(stat/=0) call cq_abort('Error allocating memory to phi_l_filter/exx !',stat)
       !
       deallocate( rho_ki_filter, STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to rho_ki_filter/exx !',stat)
       !
       deallocate( vhf_lj_filter, STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to vhf_lj_filter/exx !',stat)
       !
       deallocate( work_in, STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to work_in_filter/exx !',stat)
       !
       deallocate(work_out, STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to work_out_filter/exx !',stat)
       !
       deallocate(fftwrho3d_filter%arrayin, STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to fftwin3d_filter/exx !',stat)
       !
       deallocate(fftwrho3d_filter%arrayout, STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to fftwout3d_filter/exx !',stat)
       !
       deallocate(reckernel_3d_filter, STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to reckernel_3d_filter/exx !',stat)

    end if
    
    return
  end subroutine m_kern_exx_dummy
  !
 
 
  !!***

end module exx_kernel_default


