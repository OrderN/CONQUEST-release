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
  !!   2020/28/12 17:40 Lionel
  !!    added numerical ERI filtering (usefull for debugg)
  !!  SOURCE
  !!
  subroutine get_X_matrix( exxspin, scheme, backup_eris, niter, siter, level )

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
    use matrix_data,    only: Hrange, Srange, Xrange, SXrange
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
         tmr_std_exx_evalpao, tmr_std_exx_nsup,    &
         tmr_std_exx_setup,   tmr_std_exx_barrier, &
         tmr_std_exx_fetch,   tmr_std_exx_accumul, &
         tmr_std_exx_matmult, tmr_std_exx_poisson, &
         tmr_std_exx_allocat, tmr_std_exx_dealloc, &
         tmr_std_exx_comms,   tmr_std_exx_kernel,  &
         unit_timers_write,   unit_memory_write,   &
         unit_exx_debug,      unit_eri_debug,      &
         unit_eri_filter_debug, &         
         file_exx_timers,     file_exx_memory,     &
         file_exx_debug,      file_eri_debug,      &
         file_eri_filter_debug

    use exx_types, only: exx_alloc, r_int, extent, &
         ewald_alpha, ewald_rho, ewald_pot,   &
         pulay_radius, p_omega, p_ngauss, p_gauss, w_gauss, &
         exx_psolver,exx_pscheme, kernel,    &
         exx_total_time, eris, exx_filter, exx_gto
    !
    use exx_poisson, only: exx_scal_rho_3d, exx_ewald_rho, exx_ewald_pot
    !
    use input_module, only: fdf_boolean
    !
    implicit none

    integer, intent(in), optional :: level
    integer, intent(in)           :: exxspin, scheme
    logical, intent(in)           :: backup_eris
    
    integer, intent(in)           :: niter, siter

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

    logical      :: get_exx
    real(double) :: t0,t1
    integer      :: maxsuppfuncs, nb_eris

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
       if (fdf_boolean('IO.WriteOutToFile',.true.)) then
          if ( (niter == 1 .and. scheme == 3) .or. (niter == siter + 1  .and. (scheme == 1 .or. scheme == 2 ) )) then
             call io_assign(unit_timers_write)
             call get_file_name('exx_timers',numprocs,inode,file_exx_timers)
             open(unit_timers_write,file=file_exx_timers)
             !
             call io_assign(unit_memory_write)
             call get_file_name('exx_memory',numprocs,inode,file_exx_memory)
             open(unit_memory_write,file=file_exx_memory)
             !
          else
   
             open(unit_timers_write,file=file_exx_timers,status='old', position='append')
             open(unit_memory_write,file=file_exx_memory,status='old', position='append')
          end if
       else 
          unit_timers_write = 6
          unit_memory_write = 6
       end if
       !
    end if
    
    if ( exx_debug ) then

       if ( (niter == 1 .and. scheme == 3) .or. (niter == siter + 1  .and. (scheme == 1 .or. scheme == 2 ) )) then
          call io_assign(unit_exx_debug)
          call get_file_name('exx_debug',numprocs,inode,file_exx_debug)
          open(unit_exx_debug,file=file_exx_debug)
          !
          call io_assign(unit_eri_debug)
          call get_file_name('eri_debug',numprocs,inode,file_eri_debug)
          open(unit_eri_debug,file=file_eri_debug)       
          !
          call io_assign(unit_eri_filter_debug)
          call get_file_name('eri_filter_debug',numprocs,inode,file_eri_filter_debug)
          open(unit_eri_filter_debug,file=file_eri_filter_debug)
          !
       else
          open(unit_exx_debug,file=file_exx_debug,status='old', position='append')
          open(unit_eri_debug,file=file_eri_debug,status='old', position='append')       
          open(unit_eri_filter_debug,file=file_eri_filter_debug,status='old', position='append')
          !
       end if
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

       poisson_fftw: select case(exx_pscheme)          

       case('default')
          call exx_scal_rho_3d(inode,extent,r_int,exx_pscheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)

       case('ewald')
          call exx_mem_alloc(extent,0,0,'ewald_3d','alloc')
          call exx_ewald_rho(ewald_rho,extent,ewald_alpha,r_int)
          call exx_ewald_pot(ewald_pot,extent,ewald_alpha,r_int)          
          call exx_scal_rho_3d(inode,extent,r_int,exx_pscheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)
       
       case('pulay')
          call exx_scal_rho_3d(inode,extent,r_int,exx_pscheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)

       case('yukawa')
          call exx_scal_rho_3d(inode,extent,r_int,exx_pscheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)

       case('gauss')
          call cq_abort('EXX Gaussian representation if 1/r for solving &
               &the Poisson equation &
               &is currently under testing...')

       case default
          call exx_scal_rho_3d(inode,extent,r_int,exx_pscheme,pulay_radius, &
               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)

       end select poisson_fftw

    else if (exx_psolver == 'isf') then
       call cq_abort('EXX with ISF Poisson solver disabled')  
       !
    end if
    call stop_timer(tmr_std_exx_setup,.true.)
    !    
    if ( .not. exx_alloc ) then
       !
       if ( scheme > 0 ) then
          call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_i','alloc')       
          call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_j','alloc')
          call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_k','alloc')
          call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_l','alloc')
          !
          if ( scheme == 1 ) then
             call exx_mem_alloc(extent,maxsuppfuncs,0,'Phy_k','alloc')
             call exx_mem_alloc(extent,0,0,'Ome_kj_1d_buffer','alloc')       
             !
          end if
          !
       end if
       !
    end if
    !
    if ( scheme == 3 ) then       
       !
       if ( niter == 1 ) then
          allocate( eris( mult(S_X_SX)%ahalo%np_in_halo ), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to eris/exx !',stat)
       end if
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
    !
    call start_timer(tmr_std_exx_comms)
    call Mquest_start_send(mult(S_X_SX),mat_p(matK(  exxspin  ))%matrix,nreqs,myid,mult(S_X_SX)%prim%mx_ngonn,sends)
    call stop_timer(tmr_std_exx_comms,.true.)
    !
    ncover_yz=mult(S_X_SX)%gcs%ncovery*mult(S_X_SX)%gcs%ncoverz
    !
    !
    do kpart = 1,mult(S_X_SX)%ahalo%np_in_halo  ! Main loop
       icall=1
       ind_part = mult(S_X_SX)%ahalo%lab_hcell(kpart)    
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
             call start_timer(tmr_std_exx_comms)
             call prefetch(kpart,mult(S_X_SX)%ahalo,mult(S_X_SX)%comms,mult(S_X_SX)%bmat,icall, &
                  n_cont,part_array,mult(S_X_SX)%bindex,b_rem,lenb_rem,mat_p(matK(  exxspin  ))%matrix,   &
                  myid,ilen2,mx_msg_per_part,mult(S_X_SX)%parts,mult(S_X_SX)%prim,mult(S_X_SX)%gcs,&
                  (recv_part(nnode)-1)*2)
             call stop_timer(tmr_std_exx_comms,.true.)
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

             call start_timer(tmr_std_exx_comms)
             call end_part_comms(myid,n_cont,nbnab_rem,ibind_rem,npxyz_rem,&
                  ibpart_rem,ncover_yz,mult(S_X_SX)%gcs%ncoverz)
             call stop_timer(tmr_std_exx_comms,.true.)
             !
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
          call start_timer(tmr_std_exx_comms)
          call prefetch(kpart,mult(S_X_SX)%ahalo,mult(S_X_SX)%comms,mult(S_X_SX)%bmat,icall, &
               n_cont,part_array,mult(S_X_SX)%bindex,b_rem,lenb_rem,mat_p(matK(  exxspin  ))%matrix,   & 
               myid,ilen2,mx_msg_per_part,mult(S_X_SX)%parts,mult(S_X_SX)%prim,mult(S_X_SX)%gcs,&
               (recv_part(nnode)-1)*2)
          call stop_timer(tmr_std_exx_comms,.true.)
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
          call start_timer(tmr_std_exx_comms)
          call end_part_comms(myid,n_cont,nbnab_rem,ibind_rem,npxyz_rem,&
               ibpart_rem,ncover_yz,mult(S_X_SX)%gcs%ncoverz)
          call stop_timer(tmr_std_exx_comms,.true.)
          !
       end if ! End of the "if this isn't the first partition" loop
       k_off  = mult(S_X_SX)%ahalo%lab_hcover(kpart) ! --- offset for pbcs
       icall2 = 1
       !
       if ( scheme == 1 ) then

          if (iprint_exx > 5) write(io_lun,*) 'Proc :', myid, &
               'EXX: performing CRI calculation on kpart =', kpart
          
          call m_kern_exx_cri( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
               ibpart_rem,ibseq_rem,ibndimj_rem, & 
               b_rem,  &
               mat_p(matX(  exxspin  ))%matrix,      &
               mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
               mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
               lenb_rem, &
               mat_p(matX(  exxspin  ))%length)

       else if ( scheme == 2 ) then
          
          if (iprint_exx > 5) write(io_lun,*) 'Proc :', myid, &
               'EXX: performing on-the-fly ERI calculation on kpart =', kpart
          call m_kern_exx_eri( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
               ibpart_rem,ibseq_rem, & 
               b_rem,  &
               mat_p(matX(  exxspin  ))%matrix,      &
               mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
               mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
               lenb_rem, &
               mat_p(matX(  exxspin  ))%length, backup_eris)

       else if (scheme == 3 ) then

          if ( niter == 1 ) then
             !
             get_exx = .false.
             !
             if (iprint_exx > 5) write(io_lun,*) 'Proc :', myid, &
                  'EXX: preparing store ERI calculation on kpart =', kpart
             !
             ! First dummy call to get the number of ERIs on each proc
             !
             call m_kern_exx_dummy( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
                  ibpart_rem,ibseq_rem, &
                  b_rem,  &
                  mat_p(matX(  exxspin  ))%matrix,      &
                  mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
                  mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
                  lenb_rem, &
                  mat_p(matX(  exxspin  ))%length, nb_eris, get_exx, .false. )
             
             
             if (iprint_exx > 5) write(io_lun,*) 'Proc :', myid, &
                  'EXX: allocate ERIs on kpart =', kpart

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
             if ( exx_filter ) then
                !
                if (iprint_exx > 5) write(io_lun,*) 'Proc :', myid, &
                     'EXX: setup filtering on kpart =', kpart
                !
                call m_kern_exx_dummy( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
                     ibpart_rem,ibseq_rem, &
                     b_rem,  &
                     mat_p(matX(  exxspin  ))%matrix,      &
                     mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
                     mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
                     lenb_rem, &
                     mat_p(matX(  exxspin  ))%length, nb_eris, get_exx, exx_filter )
             end if
             !
             if (iprint_exx > 5) write(io_lun,*) 'Proc :', myid, &
                  'EXX: compute and store ERIs on kpart =', kpart
             !
             ! Third call to compute and store ERIs
             !       
             call m_kern_exx_eri( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
                ibpart_rem,ibseq_rem, & 
                b_rem,  &
                mat_p(matX(  exxspin  ))%matrix,      &
                mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
                mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
                lenb_rem, &
                mat_p(matX(  exxspin  ))%length, backup_eris)
             
          else
             !                          
             get_exx = .true.
             
             if (iprint_exx > 5) write(io_lun,*) 'Proc :', myid, &
                  'EXX: use stored ERIs to get X on kpart =', kpart

             call m_kern_exx_dummy( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
                  ibpart_rem,ibseq_rem, &
                  b_rem,  &
                  mat_p(matX(  exxspin  ))%matrix,      &
                  mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
                  mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
                  lenb_rem, &
                  mat_p(matX(  exxspin  ))%length, nb_eris, get_exx, .false. )
          end if
          
       else if ( scheme == -1 ) then

          get_exx = .false.

          if (iprint_exx > 5) write(io_lun,*) 'Proc :', myid, &
               'EXX: dummy calculation on kpart =', kpart
          
          call m_kern_exx_dummy( k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
               ibpart_rem,ibseq_rem, &
               b_rem,  &
               mat_p(matX(  exxspin  ))%matrix,      &
               mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans, &
               mult(S_X_SX)%bmat(  exxspin  )%mx_abs,mult(S_X_SX)%parts%mx_mem_grp, &
               lenb_rem, &
               mat_p(matX(  exxspin  ))%length, nb_eris, get_exx, .false. )

       end if

    end do ! End of the kpart=1,ahalo%np_in_halo loop !
    !
    !
    call start_timer(tmr_std_exx_dealloc)
    if(allocated(b_rem)) deallocate(b_rem)
    call stop_timer(tmr_std_exx_dealloc,.true.)
    !
    call start_timer(tmr_std_exx_comms)
    if(sends>0) then
       do i=1,sends
          call MPI_Wait(nreqs(i),mpi_stat,ierr)
          if(ierr/=0) call cq_abort("Error waiting for send to finish",i)
       end do
    end if
    call stop_timer(tmr_std_exx_comms,.true.)
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
    call start_timer(tmr_std_exx_dealloc)
    deallocate(ibpart_rem,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating recv_part')
    call stop_timer(tmr_std_exx_dealloc,.true.)
    !
    if(iprint_mat>3.AND.myid==0) then
       t1 = mtime()
       write(io_lun,*) 'mult time: ',t1-t0
    end if
    !
    call stop_timer(tmr_std_exx_kernel,.true.)
    !
    if ( scheme > 0 ) then
       call exx_mem_alloc(extent,0,0,'work_3d' ,'dealloc')
    end if
    !
    if ( .not. exx_alloc ) then
       !
       if ( scheme > 0 ) then
          !
          call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_i','dealloc')       
          call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_j','dealloc')
          call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_k','dealloc')    
          call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_l','dealloc')
          !
          if ( scheme == 1 ) then
             call exx_mem_alloc(extent,maxsuppfuncs,0,'Phy_k','dealloc')
             call exx_mem_alloc(extent,0,0,'Ome_kj_1d_buffer','dealloc')   
             !
          end if
          !
       end if
       !
    end if
    !
    if (exx_psolver == 'fftw')  then
       call exx_mem_alloc(extent,0,0,'fftw_3d',     'dealloc')  
       call exx_mem_alloc(extent,0,0,'reckernel_3d','dealloc')
       select case(exx_pscheme)
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
    if ( iprint_exx > 3 ) call write_mem_use(unit_memory_write,area_exx)
    !
    call stop_timer(tmr_std_exx,.true.)
    !
    ! Timers and elapse time
    exx_total_time = &
         tmr_std_exx_evalpao%t_tot + &
         tmr_std_exx_poisson%t_tot + &
         tmr_std_exx_matmult%t_tot + &
         tmr_std_exx_accumul%t_tot + &
         tmr_std_exx_nsup%t_tot    + &
         tmr_std_exx_allocat%t_tot + &
         tmr_std_exx_dealloc%t_tot + &
         tmr_std_exx_setup%t_tot   + &
         tmr_std_exx_write%t_tot   + &
         tmr_std_exx_barrier%t_tot + &
         tmr_std_exx_fetch%t_tot   + & 
         tmr_std_exx_comms%t_tot  

    if ( iprint_exx > 3 ) then
       
       call print_timer(tmr_std_exx_setup,  "exx_setup   time:", unit_timers_write)    
       call print_timer(tmr_std_exx_write,  "exx_write   time:", unit_timers_write)        
       call print_timer(tmr_std_exx_fetch,  "exx_fetch   time:", unit_timers_write)    
       call print_timer(tmr_std_exx_evalpao,"exx_evalpao time:", unit_timers_write)    
       call print_timer(tmr_std_exx_poisson,"exx_poisson time:", unit_timers_write)
       call print_timer(tmr_std_exx_matmult,"exx_matmult time:", unit_timers_write)    
       call print_timer(tmr_std_exx_accumul,"exx_accumul time:", unit_timers_write)
       call print_timer(tmr_std_exx_nsup,   "exx_nsup    time:", unit_timers_write)
       call print_timer(tmr_std_exx_allocat,"exx_allocat time:", unit_timers_write)    
       call print_timer(tmr_std_exx_dealloc,"exx_dealloc time:", unit_timers_write)    
       call print_timer(tmr_std_exx_barrier,"exx_barrier time:", unit_timers_write)
       call print_timer(tmr_std_exx_comms,  "exx_comms   time:", unit_timers_write)       
       call print_timer(tmr_std_exx_kernel, "exx_kernel  time:", unit_timers_write)    
       call print_timer(tmr_std_exx,        "exx_total   time:", unit_timers_write)    
       write(unit=unit_timers_write,fmt='("Timing: Proc ",i6,": Time spent in ", a50, " = ", &
            &f12.5," s")') inode, 'get_X_matrix',  tmr_std_exx%t_tot
       write(unit=unit_timers_write,fmt='("Timing: Proc ",i6,": Time spent in ", a50, " = ", &
            &f12.5," s")') inode, 'timer calls',  tmr_std_exx%t_tot-exx_total_time
       write(unit_timers_write,*)
       !
       if (fdf_boolean('IO.WriteOutToFile',.true.)) then
         call io_close(unit_memory_write)
         call io_close(unit_timers_write)
       end if
       !
       !
    end if

    if ( exx_debug ) then
       call io_close(unit_eri_debug)
       call io_close(unit_exx_debug)
       call io_close(unit_eri_filter_debug)       
    end if
    
    !****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_X_matrix',echo=.true.)
    !****lat>$

    return
  end subroutine get_X_matrix
  !
  !!****f* exx_kernel_default/cri_eri_inner_calculation *
  !!
  !!  NAME
  !!   cri_eri_inner_calculation
  !!
  !!  PURPOSE
  !!   Deduplicate the inner calculations of m_kern_exx_cri and m_kern_exx_eri.
  !! 
  !!  INPUTS
  !!   To ensure thread safety, variables which are altered must be passed in 
  !!   as parameters rather than imported.
  !!
  !!  AUTHOR
  !!   Connor Aird
  !!
  !!  CREATION DATE
  !!   2024/05/21  
  !!
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine cri_eri_inner_calculation(nsf1_array, phi_i, Ome_kj, nsf1, nsf_jb, nsf_kg, dv, &
               multiplier, ncaddr, ia_nsup, ewald_charge, work_out_3d, work_in_3d, &
               c, backup_eris, store_eris_ptr) 

       use exx_poisson, only: exx_v_on_grid, exx_ewald_charge

       use exx_types, only: phi_j, phi_k, ewald_rho, p_gauss, w_gauss, reckernel_3d, ewald_pot, &
            pulay_radius, p_ngauss, r_int, p_omega, exx_psolver, exx_pscheme, extent, store_eris

       use GenBlas, only: dot
       
       use numbers, only: zero
 
       implicit none

       real(double), pointer, intent(in)            :: Ome_kj(:,:,:), phi_i(:,:,:,:)
       integer,               intent(in)            :: nsf1, nsf_jb, nsf_kg        ! The indices of the loops from which this function is called 
       integer,               intent(in)            :: ncaddr, ia_nsup
       real(double),          intent(in)            :: nsf1_array(:,:,:,:), dv, multiplier
       real(double),          intent(out)           :: ewald_charge, work_out_3d(:,:,:), work_in_3d(:,:,:)
       real(double),          intent(inout)         :: c(:)
       ! Backup eris parameters. Optional as they are only needed by eri function
       logical,               intent(in)            :: backup_eris
       real(double), pointer, intent(inout), OPTIONAL :: store_eris_ptr(:,:)
       integer      :: nsf3
       real(double) :: exx_mat_elem

       work_out_3d = zero
       !
       work_in_3d = nsf1_array(:,:,:,nsf1) * phi_j(:,:,:,nsf_jb)
       !
       if (exx_psolver=='fftw' .and. exx_pscheme=='ewald') then
          call exx_ewald_charge(work_in_3d,extent,dv,ewald_charge)
          work_in_3d = work_in_3d - ewald_rho*ewald_charge
       end if
       !
       call exx_v_on_grid(inode,extent,work_in_3d,work_out_3d,r_int,   &
          exx_psolver,exx_pscheme,pulay_radius,p_omega,p_ngauss,p_gauss,&
          w_gauss,reckernel_3d)
 
       if (exx_psolver=='fftw' .and. exx_pscheme=='ewald') then
          work_out_3d = work_out_3d + ewald_pot*ewald_charge
       end if
       !
       Ome_kj = work_out_3d * phi_k(:,:,:,nsf_kg)
       !
       do nsf3 = 1, ia_nsup
          !
          ! Can we instead always store directly into store_eris_ptr(nsf2, nsf3)?
          exx_mat_elem = dot((2*extent+1)**3, phi_i(:,:,:,nsf3), 1, Ome_kj, 1) * dv * multiplier
          !
          if ( backup_eris ) then
             !
             ! eris(kpart)%store_eris( count ) = exx_mat_elem
             store_eris_ptr(nsf3, nsf_jb) = exx_mat_elem
             !
          else
             !
             c(ncaddr + nsf3 - 1) = c(ncaddr + nsf3 - 1) + exx_mat_elem
             !
          end if
          !
       end do ! nsf3 = 1, ia%nsup
  end subroutine cri_eri_inner_calculation
  !
  !!****f* exx_kernel_default/eri_gto_inner_calculation *
  !!
  !!  NAME
  !!   eri_gto_inner_calculation
  !!
  !!  PURPOSE
  !!   simplify m_kern_exx_eri to allow a simple subroutine call inside the threaded region.
  !! 
  !!  INPUTS
  !!   To ensure thread safety, variables which are altered must be passed in 
  !!   as parameters rather than imported.
  !!
  !!  AUTHOR
  !!   Connor Aird
  !!
  !!  CREATION DATE
  !!   2024/05/24  
  !!
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine eri_gto_inner_calculation(ld, kg, jb, ia, nsf_ld, nsf_kg, nsf_jb, &
               ncaddr, c, i_nt, j_nt, k_nt, l_nt, backup_eris, store_eris_ptr, &
               filter_eris_ptr, should_compute_eri_hoh)

         use exx_poisson, only: exx_v_on_grid, exx_ewald_charge

         use exx_types, only: prim_atomic_data, neigh_atomic_data, p_ngauss, store_eris

         use GenBlas, only: dot
         
         use numbers, only: zero

         use exx_erigto, only: compute_eri_hoh

         implicit none

         type(neigh_atomic_data), intent(in)    :: ld, kg, jb 
         type(prim_atomic_data),  intent(in)    :: ia
         integer,                 intent(in)    :: nsf_ld, nsf_jb, nsf_kg 
         integer,                 intent(in)    :: ncaddr
         real(double),            intent(inout) :: c(:)
         character(len=8),        intent(inout) :: i_nt, j_nt, k_nt, l_nt
         logical,                 intent(in)    :: backup_eris, should_compute_eri_hoh
         real(double), pointer,   intent(inout) :: store_eris_ptr(:,:)
         logical,      pointer,   intent(inout) :: filter_eris_ptr(:,:)
         integer      :: nsf_ia
         real(double) :: exx_mat_elem
         !
         do nsf_ia = 1, ia%nsup
            !
            exx_mat_elem = zero
            !
            if ( should_compute_eri_hoh .and. filter_eris_ptr( nsf_ia, nsf_jb ) ) then
               !
               call compute_eri_hoh( nsf_ia, nsf_jb, nsf_kg, nsf_ld, &
                     ia%spec,   jb%spec,   kg%spec,   ld%spec,  &
                     ia%xyz_ip, jb%xyz_cv, kg%xyz_cv, ld%xyz_cv,&
                     i_nt, j_nt, k_nt, l_nt,&
                     exx_mat_elem )
               !
            end if
            !
            if ( backup_eris ) then
               !
               store_eris_ptr(nsf_ia, nsf_jb) = exx_mat_elem
               !
            else
               !
               c(ncaddr + nsf_ia - 1) = c(ncaddr + nsf_ia - 1) + exx_mat_elem
               !
            end if
            !
         end do ! nsf_ia = 1, ia%nsup
  end subroutine eri_gto_inner_calculation
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
  !!   2024/02/28 Connor
  !!    Combined three instances of nsup loops into one and added omp threading
  !!   2024/05/24 Connor
  !!    Extract inner calculation into cri_eri_inner_calculation to de-duplicate code
  !!
  !!  SOURCE
  !!
  subroutine m_kern_exx_cri(k_off, kpart, ib_nd_acc, ibaddr, nbnab, &
       ibpart, ibseq, bndim2, b, c, ahalo, chalo, at, mx_absb, mx_part,  &
       lenb, lenc )

    use numbers,        only: zero, one
    use matrix_module,  only: matrix_halo, matrix_trans
    use global_module,  only: area_exx
    use GenBlas,        only: dot
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
    use exx_types, only: prim_atomic_data, neigh_atomic_data,tmr_std_exx_accumul,  &
         grid_spacing, r_int, extent, ewald_charge,p_ngauss,unit_exx_debug,        &
         tmr_std_exx_nsup,phi_i_1d_buffer, phi_j, phi_k, phi_l,Phy_k,              &
         Ome_kj_1d_buffer,work_in_3d,work_out_3d,exx_alloc
    !
    use exx_memory,  only: exx_mem_alloc 
    use exx_poisson, only: exx_v_on_grid, exx_ewald_charge
    !
    use exx_module,  only: get_halodat, get_iprimdat
    !
    implicit none
    !
    ! Passed variables
    type(matrix_halo),  intent(in)    :: ahalo, chalo
    type(matrix_trans), intent(in)    :: at
    integer,            intent(in)    :: mx_absb, mx_part, lenb, lenc
    integer,            intent(in)    :: kpart, k_off
    real(double),       intent(in)    :: b(lenb)
    real(double),       intent(inout) :: c(lenc)
    !
    ! Remote indices
    integer(integ), intent(in) :: ib_nd_acc(mx_part)
    integer(integ), intent(in) :: ibaddr(mx_part)
    integer(integ), intent(in) :: nbnab(mx_part)
    integer(integ), intent(in) :: ibpart(mx_part*mx_absb)
    integer(integ), intent(in) :: ibseq(mx_part*mx_absb)
    integer(integ), intent(in) :: bndim2(mx_part*mx_absb)
    !
    ! Local variables
    integer :: jbnab2ch(mx_absb)  ! Automatic array
    integer :: nbkbeg, k, k_in_part, k_in_halo, j, jpart, jseq
    integer :: i, i_in_prim, icad, nbbeg, j_in_halo, ncbeg
    integer :: nb_nd_kbeg
    integer :: nd1, nd3
    integer :: nbaddr, ncaddr
    integer :: l, lseq, lpart
    integer :: np, ni
    !
    real(double), dimension(3) :: xyz_zero  = zero
    real(double)               :: dr,dv,K_val
    !
    !
    ! We allocate pointers here to point at 1D arrays later and allow contiguous access when passing to BLAS dot later
    real(double), pointer :: phi_i(:,:,:,:), Ome_kj(:,:,:)
    !
    type(prim_atomic_data)  :: ia !i_alpha
    type(neigh_atomic_data) :: jb !j_beta
    type(neigh_atomic_data) :: kg !k_gamma
    type(neigh_atomic_data) :: ld !l_delta
    !
    integer :: maxsuppfuncs, nsf_kg, nsf_ld, nsf_jb, count
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
       if ( exx_alloc ) call exx_mem_alloc(extent,kg%nsup,0,'phi_k','alloc')
    
       call exx_phi_on_grid(inode,kg%global_num,kg%spec,extent, &
            xyz_zero,kg%nsup,phi_k,r_int,xyz_zero)             
       !
       jbnab2ch = 0
       do j = 1, nbnab(k_in_part)
          jpart = ibpart(nbkbeg+j-1) + k_off
          jseq  = ibseq (nbkbeg+j-1)
          jbnab2ch(j) = chalo%i_halo(chalo%i_hbeg(jpart)+jseq-1)
       end do
       !
       nbbeg = nb_nd_kbeg
       !
       if ( exx_alloc ) call exx_mem_alloc(extent,kg%nsup,0,'Phy_k','alloc')
       !
       Phy_k = zero
       !
!!$
!!$ ****[ l do loop ]****
!!$
       do l = 1, nbnab(k_in_part)
          lpart = ibpart(nbkbeg+l-1) + k_off
          lseq  = ibseq (nbkbeg+l-1)
          call get_halodat(ld,kg,lseq,chalo%i_hbeg(lpart),         &
               BCS_parts%lab_cell(BCS_parts%inv_lab_cover(lpart)), &
               'l',.true.,unit_exx_debug)
          !
!!$
!!$ ****[ <kl> screening ]****
!!$
          !
          if ( exx_alloc ) call exx_mem_alloc(extent,ld%nsup,0,'phi_l','alloc')
          !
          call exx_phi_on_grid(inode,ld%global_num,ld%spec,extent,     &
               ld%xyz,ld%nsup,phi_l,r_int,xyz_zero)             
          !
          do nsf_ld = 1, ld%nsup
             !
             nbaddr = nbbeg + kg%nsup * (nsf_ld - 1)
             !
             do nsf_kg = 1, kg%nsup                         
                !
                K_val = b(nbaddr+nsf_kg-1)
                !
                call start_timer(tmr_std_exx_accumul)
                Phy_k(:,:,:,nsf_kg) = Phy_k(:,:,:,nsf_kg) + K_val*phi_l(:,:,:,nsf_ld) 
                call stop_timer(tmr_std_exx_accumul,.true.)

             end do
          end do
          !
          nbbeg = nbbeg + ld%nsup*kg%nsup
          !
          if ( exx_alloc ) call exx_mem_alloc(extent,ld%nsup,0,'phi_l','dealloc')
          
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
          !
          call get_iprimdat(ia,kg,ni,i_in_prim,np,.true.,unit_exx_debug)
          if (ia%nsup/=ahalo%ndimi(i_in_prim)) call cq_abort('Error1: ',ia%nsup,ahalo%ndimi(i_in_prim))
          !
          if ( exx_alloc ) call exx_mem_alloc(extent,ia%nsup,0,'phi_i_1d_buffer','alloc')
          phi_i(1:2*extent+1, 1:2*extent+1, 1:2*extent+1, 1:ia%nsup) => phi_i_1d_buffer
          !
          call exx_phi_on_grid(inode,ia%ip,ia%spec,extent, &
               ia%xyz,ia%nsup,phi_i,r_int,xyz_zero)    
          !
!!$
!!$ ****[ j loop ]****
!!$
          do j = 1, nbnab(k_in_part)
             nbbeg     = nb_nd_kbeg
             j_in_halo = jbnab2ch(j) !***
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
                   if (jb%nsup/=bndim2(nbkbeg+j-1)) call cq_abort('Error2: ',jb%nsup,bndim2(nbkbeg+j-1))
                   !
                   if ( exx_alloc ) call exx_mem_alloc(extent,jb%nsup,0,'phi_j','alloc')
                   !
                   call exx_phi_on_grid(inode,jb%global_num,jb%spec,extent, &
                        jb%xyz,jb%nsup,phi_j,r_int,xyz_zero)             
                   !
                   ! We allocate a 1D buffer here, instead of 3D, to allow contiguous access when passing to BLAS dot later
                   !
                   if ( exx_alloc ) call exx_mem_alloc(extent,0,0,'Ome_kj_1d_buffer','alloc')
                   !
                   call start_timer(tmr_std_exx_nsup)
                   !
                   ! Begin the parallel region here as earlier allocations make it difficult to do before now. 
                   ! However, this should be possible in future work. 
                   !
                   !$omp parallel default(none) reduction(+: c)            &
                   !$omp     shared(kg,jb,Phy_k,ncbeg,ia,phi_i,extent,dv) &
                   !$omp     private(nsf_kg,nsf_jb,work_out_3d,work_in_3d,ewald_charge,Ome_kj_1d_buffer, &
                   !$omp             Ome_kj,ncaddr)
                   Ome_kj(1:2*extent+1, 1:2*extent+1, 1:2*extent+1) => Ome_kj_1d_buffer
                   !$omp do schedule(dynamic) collapse(2)
                   do nsf_kg = 1, kg%nsup
                      do nsf_jb = 1, jb%nsup
                         !
                         ncaddr = ncbeg + ia%nsup * (nsf_jb - 1)
                         !
                         call cri_eri_inner_calculation(Phy_k, phi_i, Ome_kj, nsf_kg, nsf_jb, nsf_kg, dv, &
                                       1.0d0, ncaddr, ia%nsup, ewald_charge, work_out_3d, work_in_3d, c,  &
                                       .false.)
                         !
                      end do ! nsf_ld = 1, jb%nsup
                   end do ! nsf_kg = 1, kg%nsup
                   !$omp end do
                   !$omp end parallel
                   !
                   call stop_timer(tmr_std_exx_nsup,.true.)
                   !
                   if ( exx_alloc ) call exx_mem_alloc(extent,0,0,'Ome_kj_1d_buffer','dealloc')
                   if ( exx_alloc ) call exx_mem_alloc(extent,jb%nsup,0,'phi_j','dealloc')
                   !
                end if ! ( ncbeg /=0 )
             end if ! ( j_in_halo /=0 )
!!$
!!$ ****[ j end loop ]****
!!$
          end do ! End of j = 1, nbnab(k_in_part)           
          !          
!!$
!!$ ****[ i end loop ]****
!!$
          !
          if ( exx_alloc ) call exx_mem_alloc(extent,ia%nsup,0,'phi_i_1d_buffer','dealloc')
          !
       end do ! End of i = 1, at%n_hnab(k_in_halo)
       !
!!$
!!$ ****[ k end loop ]****
!!$
       !
       if ( exx_alloc ) call exx_mem_alloc(extent,kg%nsup,0,'Phy_k','dealloc')
       if ( exx_alloc ) call exx_mem_alloc(extent,kg%nsup,0,'phi_k','dealloc')
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
  !!   Compute EXX matrix using GTO-ERI engine. For now
  !!   default is the Taketa, Huzinaga, and O-ohata (HOH)
  !!   approach as given in exx_erigto module
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
  !!   2024/05/24 Connor
  !!    - Extract inner calculation into cri_eri_inner_calculation
  !!    - Combine with m_kern_exx_eri_gto
  !!    - Add omp threading 
  !!
  !!  SOURCE
  !!
  subroutine m_kern_exx_eri(k_off, kpart, ib_nd_acc, ibaddr, nbnab, &
       ibpart, ibseq, b, c, ahalo, chalo, at, mx_absb, mx_part,     &
       lenb, lenc, backup_eris )

    use numbers,        only: zero
    use matrix_module,  only: matrix_halo, matrix_trans
    !
    use primary_module, only: bundle 
    use cover_module,   only: BCS_parts
    !
    use exx_evalpao,    only: exx_phi_on_grid
    !
    use exx_evalgto,    only: exx_gto_on_grid_prim
    !
    use exx_types, only: prim_atomic_data,neigh_atomic_data,grid_spacing,r_int, &
         extent,ewald_charge,p_ngauss,unit_exx_debug,tmr_std_exx_nsup,eris,     &
         phi_i_1d_buffer,phi_j,phi_k,phi_l,exx_alloc,                           &
         Ome_kj_1d_buffer,work_in_3d,work_out_3d, exx_gto
    !
    use exx_memory, only: exx_mem_alloc
    !
    use exx_module, only: get_halodat, get_iprimdat
    !
    use exx_poisson,only: exx_v_on_grid, exx_ewald_charge
    !
    ! m_kern_exx_eri_gto imports
    use exx_erigto, only: compute_eri_hoh
    !
    implicit none
    !
    ! Passed variables
    type(matrix_halo),  intent(in)    :: ahalo, chalo
    type(matrix_trans), intent(in)    :: at
    integer,            intent(in)    :: mx_absb, mx_part, lenb, lenc
    integer,            intent(in)    :: kpart, k_off
    real(double),       intent(in)    :: b(lenb)
    real(double),       intent(inout) :: c(lenc)
    logical,            intent(in)    :: backup_eris
    !
    ! Remote indices
    integer(integ), intent(in) :: ib_nd_acc(mx_part)
    integer(integ), intent(in) :: ibaddr(mx_part)
    integer(integ), intent(in) :: nbnab(mx_part)
    integer(integ), intent(in) :: ibpart(mx_part*mx_absb)
    integer(integ), intent(in) :: ibseq(mx_part*mx_absb)
    !
    ! Local variables
    integer :: jbnab2ch(mx_absb)  ! Automatic array
    integer :: nbkbeg, k, k_in_part, k_in_halo, j, jpart, jseq
    integer :: i, i_in_prim, icad, nbbeg, j_in_halo, ncbeg
    integer :: nb_nd_kbeg
    integer :: nbaddr, ncaddr
    integer :: l, lseq, lpart
    integer :: np, ni
    !
    real(double), dimension(3) :: xyz_zero  = zero
    !
    real(double)               :: dv,K_val
    !
    ! We allocate pointers here to point at 1D arrays later and allow contiguous access when passing to BLAS dot later
    real(double), pointer :: phi_i(:,:,:,:), Ome_kj(:,:,:), store_eris_ptr(:,:)
    logical,      pointer :: filter_eris_ptr(:,:)
    !
    type(prim_atomic_data)  :: ia !i_alpha
    type(neigh_atomic_data) :: jb !j_beta
    type(neigh_atomic_data) :: kg !k_gamma
    type(neigh_atomic_data) :: ld !l_delta
    !
    integer                 :: nsf_kg, nsf_ld, nsf_jb, count
    !
    logical :: should_allocate, should_compute_eri_hoh
    !
    ! m_kern_exx_eri_gto variables
    character(len=8) :: i_nt, j_nt, k_nt, l_nt
    !
    dv = grid_spacing**3
    count = 0
    should_allocate = exx_alloc .and. (.not. exx_gto)
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
       call get_halodat(kg,kg,k_in_part,ahalo%i_hbeg(ahalo%lab_hcover(kpart)), &
            ahalo%lab_hcell(kpart),'k',.true.,unit_exx_debug)
       !
       if ( should_allocate ) call exx_mem_alloc(extent,kg%nsup,0,'phi_k','alloc')
       !
       if (.not. exx_gto) call exx_phi_on_grid(inode,kg%global_num,kg%spec,extent, &
                              kg%xyz,kg%nsup,phi_k,r_int,xyz_zero)
       !
       jbnab2ch = 0
       do j = 1, nbnab(k_in_part)
          jpart = ibpart(nbkbeg+j-1) + k_off
          jseq  = ibseq (nbkbeg+j-1)
          jbnab2ch(j) = chalo%i_halo(chalo%i_hbeg(jpart)+jseq-1)
       end do
       !
       nbbeg = nb_nd_kbeg
       !
!!$
!!$ ****[ l do loop ]****
!!$
       l_loop: do l = 1, nbnab(k_in_part)
          lpart = ibpart(nbkbeg+l-1) + k_off
          lseq  = ibseq (nbkbeg+l-1)
          call get_halodat(ld,kg,lseq,chalo%i_hbeg(lpart),         &
               BCS_parts%lab_cell(BCS_parts%inv_lab_cover(lpart)), &
               'l',.true.,unit_exx_debug)
          !
          if ( should_allocate ) call exx_mem_alloc(extent,ld%nsup,0,'phi_l','alloc')
          !
          if (.not. exx_gto) call exx_phi_on_grid(inode,ld%global_num,ld%spec,extent,     &
                                 ld%xyz,ld%nsup,phi_l,r_int,xyz_zero)
          !
          ld_loop: do nsf_ld = 1, ld%nsup
             !
             nbaddr = nbbeg + kg%nsup * (nsf_ld - 1)
             !
             kg_loop: do nsf_kg = 1, kg%nsup                         
                !
                if ( backup_eris .and. (.not. exx_gto) ) then
                   K_val = real(1,double)
                else
                   K_val = b(nbaddr+nsf_kg-1)
                end if
!!$
!!$ ****[ i loop ]****
!!$
                i_loop: do i = 1, at%n_hnab(k_in_halo)
                   i_in_prim = at%i_prim(at%i_beg(k_in_halo)+i-1)
                   ni  = bundle%iprim_seq (i_in_prim)
                   np  = bundle%iprim_part(i_in_prim)
                   icad  = (i_in_prim - 1) * chalo%ni_in_halo
                   !
                   call get_iprimdat(ia,kg,ni,i_in_prim,np,.true.,unit_exx_debug)          
                   !
                   if ( should_allocate ) call exx_mem_alloc(extent,ia%nsup,0,'phi_i_1d_buffer','alloc')
                   phi_i(1:2*extent+1, 1:2*extent+1, 1:2*extent+1, 1:ia%nsup) => phi_i_1d_buffer
                   !
                   if (.not. exx_gto) call exx_phi_on_grid(inode,ia%ip,ia%spec,extent, &
                                          ia%xyz,ia%nsup,phi_i,r_int,xyz_zero)
                   !
!!$
!!$ ****[ j loop ]****
!!$
                   j_loop: do j = 1, nbnab(k_in_part)                  
                      j_in_halo = jbnab2ch(j)
                      !
                      if ( j_in_halo /= 0 ) then
                         !
                         ncbeg = chalo%i_h2d(icad + j_in_halo)
                         !
                         if ( ncbeg /= 0 ) then 
                            jpart = ibpart(nbkbeg+j-1) + k_off
                            jseq  = ibseq (nbkbeg+j-1)
                            !
                            call get_halodat(jb,kg,jseq,chalo%i_hbeg(jpart),         &
                                 BCS_parts%lab_cell(BCS_parts%inv_lab_cover(jpart)), &
                                 'j',.true.,unit_exx_debug)
                            !
                            if ( should_allocate ) call exx_mem_alloc(extent,jb%nsup,0,'phi_j','alloc')
                            !
                            if (.not. exx_gto) call exx_phi_on_grid(inode,jb%global_num,jb%spec,extent, &
                                                   jb%xyz,jb%nsup,phi_j,r_int,xyz_zero)
                            !
                            if ( should_allocate ) call exx_mem_alloc(extent,0,0,'Ome_kj_1d_buffer','alloc')
                            !
                            call start_timer(tmr_std_exx_nsup)
                            !
                            ! Point at the next block of eris to store and update counter 
                            store_eris_ptr(1:ia%nsup, 1:jb%nsup) => eris(kpart)%store_eris(count+1:count + (jb%nsup * ia%nsup))
                            filter_eris_ptr(1:ia%nsup, 1:jb%nsup) => eris(kpart)%filter_eris(count+1:count + (jb%nsup * ia%nsup))
                            count = count + (jb%nsup * ia%nsup)
                            !
                            should_compute_eri_hoh = abs(ia%xyz_ip(3)-kg%xyz_cv(3)) < ( ia%radi + kg%radi) &
                                             .and. abs(jb%xyz_cv(3)-ld%xyz_cv(3)) < ( jb%radi + ld%radi)
                            !
                            !$omp parallel default(none) reduction(+: c)                                                   &
                            !$omp     shared(ld,kg,jb,ia,nsf_kg,nsf_ld,ncbeg,phi_k,phi_j,phi_l,phi_i,extent,dv,eris,K_val, &
                            !$omp            backup_eris,phi_i_1d_buffer,kpart,store_eris_ptr,filter_eris_ptr,exx_gto,      &
                            !$omp            should_compute_eri_hoh)                                                       &
                            !$omp     private(nsf_jb,work_out_3d,work_in_3d,ewald_charge,Ome_kj_1d_buffer,Ome_kj,ncaddr,   &
                            !$omp             i_nt,j_nt,k_nt,l_nt)
                            !
                            Ome_kj(1:2*extent+1, 1:2*extent+1, 1:2*extent+1) => Ome_kj_1d_buffer
                            !$omp do schedule(dynamic)
                            jb_loop: do nsf_jb = 1, jb%nsup
                               !
                               ncaddr = ncbeg + ia%nsup * (nsf_jb - 1)
                               !
                               if (exx_gto) then
                                 !
                                 call eri_gto_inner_calculation(ld, kg, jb, ia, nsf_ld, nsf_kg, nsf_jb, ncaddr, c, i_nt, j_nt, &
                                             k_nt, l_nt, backup_eris, store_eris_ptr, filter_eris_ptr, should_compute_eri_hoh)
                                 !
                               else 
                                 !
                                 call cri_eri_inner_calculation(phi_l, phi_i, Ome_kj, nsf_ld, nsf_jb, nsf_kg, dv, K_val, &
                                             ncaddr, ia%nsup, ewald_charge, work_out_3d, work_in_3d, c, &
                                             backup_eris, store_eris_ptr)
                                 !
                               end if 
                            end do jb_loop
                            !$omp end do
                            !$omp end parallel
                            !
                            call stop_timer(tmr_std_exx_nsup,.true.)
                            !
                            if ( should_allocate ) call exx_mem_alloc(extent,0,0,'Ome_kj_1d_buffer','dealloc')
                            if ( should_allocate ) call exx_mem_alloc(extent,jb%nsup,0,'phi_j','dealloc')                                
                            !
                         end if
                         !
                      end if
                      !
!!$
!!$ ****[ j end loop ]****
!!$                      
                      !
                   end do j_loop
                   !
                   if ( should_allocate ) call exx_mem_alloc(extent,ia%nsup,0,'phi_i_1d_buffer','dealloc')
                   !
!!$
!!$ ****[ i end loop ]****
!!$
                   !
                end do i_loop
                !
             end do kg_loop
             !
          end do ld_loop
          !
          if ( should_allocate ) call exx_mem_alloc(extent,ld%nsup,0,'phi_l','dealloc')  
          !
!!$
!!$ ****[ l end loop ]****
!!$
          !
       end do l_loop
       !
       nbbeg = nbbeg + ld%nsup*kg%nsup             
       !
       if ( should_allocate ) call exx_mem_alloc(extent,kg%nsup,0,'phi_k','dealloc')
       !
!!$
!!$ ****[ k end loop ]****
!!$
       !
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
  !!   (1) count the the number of ERIs on each proc
  !!   to allocate arrays when storage is required
  !!   (2) if filter_eris = .true. prepare EXX calculation
  !!   such as filtering and screening using coarse grid Poisson
  !!   solver for ERI evaluation. Filtering is used only if
  !!   GTO fitting and GTO-ERI engine are used.
  !!   (3) if compute_exx =.true. compute EXX matrix using
  !!   stored ERIs
  !!   
  !!  INPUTS
  !!
  !!  AUTHOR
  !!   L.A.Truflandier
  !!
  !!  CREATION DATE
  !!   2020/10/27
  !!
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine m_kern_exx_dummy(k_off, kpart, ib_nd_acc, ibaddr, nbnab, &
       ibpart, ibseq, b, c, ahalo, chalo,  & 
       at, mx_absb, mx_part, lenb, lenc, &
       count_eris, compute_exx, filter_eris)

    use global_module,  only: iprint_exx 
    use numbers,        only: zero, one
    use matrix_module,  only: matrix_halo, matrix_trans
    use global_module,  only: area_exx
    !
    use fft_interface_module, only: fft3_init_wrapper
    
    use basic_types,    only: primary_set
    use primary_module, only: bundle 
    use matrix_data,    only: Hrange, SXrange, Xrange, Srange
    use mult_module,    only: return_matrix_value, S_X_SX
    use cover_module,   only: BCS_parts
    !
    use species_module, only: nsf_species
    !
    use angular_coeff_routines, only: calc_mat_elem_gen
    !    
    use exx_evalpao, only: exx_phi_on_grid
    use exx_types,   only: reckernel_3d_filter
    use exx_types,   only: prim_atomic_data,neigh_atomic_data, &         
         grid_spacing,r_int,unit_exx_debug,pulay_radius,eris,  &
         p_omega,p_ngauss,p_gauss,w_gauss,exx_psolver,         &
         exx_pscheme,exx_filter_thr,exx_filter_extent,         &
         unit_eri_filter_debug

    use exx_module, only: get_halodat, get_iprimdat
    !
    use exx_poisson,only: exx_v_on_grid, exx_scal_rho_3d
    !
    implicit none
    !
    ! Passed variables
    type(matrix_halo),  intent(in)    :: ahalo, chalo
    type(matrix_trans), intent(in)    :: at
    integer,            intent(in)    :: mx_absb, mx_part, lenb, lenc
    integer,            intent(in)    :: kpart, k_off
    real(double),       intent(in)    :: b(lenb)
    real(double),       intent(inout) :: c(lenc)
    integer, intent(inout) :: count_eris
    logical, intent(in)    :: compute_exx
    logical, intent(in)    :: filter_eris
    !
    ! Remote indices
    integer(integ), intent(in) :: ib_nd_acc(mx_part)
    integer(integ), intent(in) :: ibaddr(mx_part)
    integer(integ), intent(in) :: nbnab(mx_part)
    integer(integ), intent(in) :: ibpart(mx_part*mx_absb)
    integer(integ), intent(in) :: ibseq(mx_part*mx_absb)
    !
    ! Local variables
    integer :: jbnab2ch(mx_absb)  ! Automatic array
    integer :: nbkbeg, k, k_in_part, k_in_halo, j, jpart, jseq
    integer :: i, i_in_prim, icad, nbbeg, j_in_halo, ncbeg
    integer :: nb_nd_kbeg
    integer :: nd1, nd3
    integer :: nbaddr, ncaddr
    integer :: l, lseq, lpart
    integer :: np, ni
    !
    real(double), dimension(3) :: xyz_zero  = zero
    !
    real(double)               ::   dr, dv, K_val
    real(double)               ::   exx_mat_elem
    !
    type(prim_atomic_data)  :: ia !i_alpha
    type(neigh_atomic_data) :: jb !j_beta
    type(neigh_atomic_data) :: kg !k_gamma
    type(neigh_atomic_data) :: ld !l_delta
    !
    integer                 :: maxsuppfuncs
    integer                 :: nsf_kg, nsf_ld, nsf_ia, nsf_jb
    integer                 :: r, s, t, stat
    !
    real(double), dimension(:,:,:,:),   allocatable :: phi_i_filter,  phi_j_filter
    real(double), dimension(:,:,:,:),   allocatable :: phi_k_filter,  phi_l_filter
    real(double), dimension(:,:,:),     allocatable :: work_in, work_out
    !
    dr = grid_spacing
    dv = dr**3
    !
    count_eris = 0
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
       allocate( work_in(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1), STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to work_in_filter/exx !',stat)
       !
       allocate(work_out(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1), STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to work_out_filter/exx !',stat)
       !
       call fft3_init_wrapper( 2*exx_filter_extent+1  )
       !
       allocate(reckernel_3d_filter(2*exx_filter_extent+1,2*exx_filter_extent+1,2*exx_filter_extent+1), STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to reckernel_3d_filter/exx !',stat)

       call exx_scal_rho_3d(inode,exx_filter_extent,r_int,exx_pscheme,pulay_radius, &
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
       if ( filter_eris ) then
          call exx_phi_on_grid(inode,kg%global_num,kg%spec,exx_filter_extent, &
               xyz_zero,maxsuppfuncs,phi_k_filter,r_int,xyz_zero)
       end if
       
       jbnab2ch = 0
       do j = 1, nbnab(k_in_part)
          jpart = ibpart(nbkbeg+j-1) + k_off
          jseq  = ibseq (nbkbeg+j-1)
          jbnab2ch(j) = chalo%i_halo(chalo%i_hbeg(jpart)+jseq-1)
       end do
       !
       nbbeg = nb_nd_kbeg
       !
!!$
!!$ ****[ l do loop ]****
!!$
       l_loop: do l = 1, nbnab(k_in_part)
          lpart = ibpart(nbkbeg+l-1) + k_off
          lseq  = ibseq (nbkbeg+l-1)
          call get_halodat(ld,kg,lseq,chalo%i_hbeg(lpart),         &
               BCS_parts%lab_cell(BCS_parts%inv_lab_cover(lpart)), &
               'l',.true.,unit_exx_debug)
          !
          if ( filter_eris ) then
             !
             call exx_phi_on_grid(inode,ld%global_num,ld%spec,exx_filter_extent, &
                  ld%xyz,maxsuppfuncs,phi_l_filter,r_int,xyz_zero)
             !
          end if
          !
          ld_loop: do nsf_ld = 1, ld%nsup
             !
             nbaddr = nbbeg + kg%nsup * (nsf_ld - 1)
             !
             kg_loop: do nsf_kg = 1, kg%nsup
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
                   !
                   call get_iprimdat(ia,kg,ni,i_in_prim,np,.true.,unit_exx_debug)          
                   !
                   if ( filter_eris ) then                      
                      call exx_phi_on_grid(inode,ia%ip,ia%spec,exx_filter_extent, &
                           ia%xyz,maxsuppfuncs,phi_i_filter,r_int,xyz_zero)
                   end if
                   ! 
!!$
!!$ ****[ j loop ]****
!!$
                   j_loop: do j = 1, nbnab(k_in_part)
                      j_in_halo = jbnab2ch(j) !***
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
                            end if
                            !
                            jb_loop: do nsf_jb = 1, jb%nsup                                         
                               !
                               ncaddr = ncbeg + ia%nsup * (nsf_jb - 1)
                               !
                               if ( filter_eris ) then
                                  work_out = zero
                                  work_in  = phi_l_filter(:,:,:,nsf_ld)*phi_j_filter(:,:,:,nsf_jb)
                                  !
                                  call exx_v_on_grid(inode,exx_filter_extent,work_in,work_out,r_int, &
                                       exx_psolver,exx_pscheme,pulay_radius,p_omega,p_ngauss,p_gauss,   &
                                       w_gauss,reckernel_3d_filter)
                                  !
                               end if
                               !
                               ia_loop: do nsf_ia = 1, ia%nsup
                                  !
                                  if ( compute_exx ) then
                                     !
                                     c(ncaddr + nsf_ia - 1) = c(ncaddr + nsf_ia - 1) + &
                                          eris(kpart)%store_eris(count_eris + 1) * K_val
                                     !
                                  else
                                     !
                                     if ( filter_eris ) then
                                        !
                                        exx_mat_elem = zero
                                        !
                                        do r = 1, 2*exx_filter_extent+1
                                           do s = 1, 2*exx_filter_extent+1
                                              do t = 1, 2*exx_filter_extent+1                         
                                                 
                                                 exx_mat_elem = exx_mat_elem &      
                                                      + phi_k_filter(t,s,r,nsf_kg)  &
                                                      * phi_i_filter(t,s,r,nsf_ia)  &
                                                      * work_out(t,s,r) * dv
                                              end do
                                           end do
                                        end do
                                        !
                                        if ( abs(exx_mat_elem) < exx_filter_thr ) then
                                           !
                                           eris(kpart)%filter_eris( count_eris + 1) = .false.
                                           !
                                        end if
                                        !
                                        if (exx_debug) then
                                           
                                           write(unit_eri_filter_debug,10) count_eris+1, exx_mat_elem, &
                                                eris(kpart)%filter_eris(count_eris + 1),               &
                                                '[',ia%ip, kg%global_num,'|',ld%global_num, jb%global_num,']', &
                                                '(',nsf_ia,nsf_kg,  '|',nsf_ld,nsf_jb,  ')' , &
                                                '[',ia%name,kg%name,'|',ld%name,jb%name,']' , &
                                                K_val, kg%r, ld%r, ia%r, jb%r
                                           
                                        end if
                                        !
                                     end if
                                     !
                                  end if
                                  !
                                  count_eris = count_eris + 1
                                  !
                               end do ia_loop
                               !
                               !
                            end do jb_loop
                            !
                         end if
                         !
                      end if
                      !
!!$
!!$ ****[ j end loop ]****
!!$
                      !
                   end do j_loop
                   !
!!$
!!$ ****[ i end loop ]****
!!$
                   !
                end do i_loop
                
             end do kg_loop
             !
          end do ld_loop
          !
          nbbeg = nbbeg + ld%nsup*kg%nsup       
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
    
    if (iprint_exx > 5) write(io_lun,*) 'Proc :', myid, &
         'kpart:', kpart, 'nb. of ERIs:', count_eris
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
       deallocate( work_in, STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to work_in_filter/exx !',stat)
       !
       deallocate(work_out, STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to work_out_filter/exx !',stat)
       !
       deallocate(reckernel_3d_filter, STAT=stat)
       if(stat/=0) call cq_abort('Error allocating memory to reckernel_3d_filter/exx !',stat)

    end if

10  format(I8,X,F16.10,X,L,X,A,2I4,A,2I4,A,4X,A,2I4,A,2I4,A,A,2A4,A,2A4,A,X,8F12.6)
    
    
    return
  end subroutine m_kern_exx_dummy
  !
   subroutine plot1d_obj(obj_on_grid,extent,r_int,unit,filename)

    use numbers, ONLY: zero, twopi
    
    implicit none
    
    ! << Passed variables >>
    integer,     intent(in)    :: extent
    integer,     intent(in)    :: unit
    real(double),intent(in)    :: r_int
    character(*),intent(in)    :: filename

    real(double),intent(inout) :: obj_on_grid(2*extent+1,2*extent+1,2*extent+1)

    ! << Local variables >>
    
    real(double) :: h
    integer      :: i
 

    h = r_int/real(extent,double) 
        
    open(unit,file=filename)
    
    do i = 1, 2*extent+1
       write(unit,'(3F24.16)') &
            real(i-extent-1,double)*h, &
            obj_on_grid(extent+1,extent+1,i)
    end do

    close(unit)    

    return
  end subroutine plot1d_obj 
  !!***

end module exx_kernel_default


