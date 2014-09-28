module exx_memory

  use datatypes

  use exx_types,         ONLY: phi_i,   phi_j,  phi_k, phi_l
  use exx_types,         ONLY: rho_kj,  rho_lj, rho_ik
  use exx_types,         ONLY: Ome_j,   Ome_i, vhf_kj, Phy_k, vhf_lj
  use exx_types,         ONLY: phi_l_m, phi_k_m, Ome_kj
  
  use exx_types,         ONLY: EXXmat,Kmat_ij,Kmat_kl,EXX_tmp,Kmat_tmp,EXX_tmp2, Imat
  use exx_types,         ONLY: work_in_3d,work_out_3d 
  use exx_types,         ONLY: reckernel_3d,ewald_rho,ewald_pot
  use exx_types,         ONLY: isf_rho,isf_pot_ion

  use exx_types,         ONLY: fftwrho3d,fftwrho2d,fftwrho1d

  use exx_types,         ONLY: tmr_std_exx_allocat, tmr_std_exx_dealloc
  use exx_types,         ONLY: unit_memory_write

  use global_module,             ONLY: area_exx, io_lun
  use timer_module,              ONLY: start_timer, stop_timer, print_timer
  use memory_module,             ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
  use GenComms,                  ONLY: cq_abort

  !use fftw_module
  use fft_interface_module,      ONLY: fft3_init_wrapper
  !use fft_interface_module,      ONLY: fft3_dest_wrapper
  !use fft_interface_module,      ONLY: FFTW_FORWARD, FFTW_BACKWARD, FFTW_ESTIMATE

  implicit none

contains

  subroutine exx_mem_alloc(extent,nsf1,nsf2,matrix,flag,unit,n_neigh,neigh)
    
    use numbers,     ONLY: zero, one
    use GenComms,    ONLY: inode, ionode

    implicit none

    integer, parameter   :: FFTW_FORWARD  = -1
    integer, parameter   :: FFTW_BACKWARD = +1
    integer, parameter   :: FFTW_ESTIMATE = 64

    ! << Passed variables >>
    character(*), intent(in) :: flag
    character(*), intent(in) :: matrix

    integer, intent(in) :: extent     ! grid size (2*extent+1)
    integer, intent(in) :: nsf1, nsf2 ! number of support functions

    integer, optional, intent(in) :: unit    ! for writing    
    integer, optional, intent(in) :: n_neigh ! number of neigh. atoms for vhf
    integer, optional, intent(in) :: neigh   ! neigh. no.

    ! << Local variables >>
    integer             :: ierr, stat, planR, planF, lun, i, j
    logical, parameter  :: accumulate = .true.

    lun  = unit_memory_write
    stat = 0

    if(flag=='alloc') then
       call start_timer(tmr_std_exx_allocat)       
       !
       select case (matrix)          
          !
          !===PAO=allocation========================================================================
       case('phi_i') ! allocate phi_i for primary atom
          allocate(phi_i(2*extent+1,2*extent+1,2*extent+1,nsf1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to phi_i/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          phi_i = zero
          !write(*,*) '\phi_{i\alpha} allocated', shape(phi_i)
          
       case('phi_j') ! allocate phi_j for neighbour atom [Srange]
          allocate(phi_j(2*extent+1,2*extent+1,2*extent+1,nsf1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to phi_j/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          phi_j = zero
          !write(*,*) '\phi_{i\beta} allocated', shape(phi_j)
          
       case('phi_k') ! allocate phi_k for neighbour atom [Xrange]
          allocate(phi_k(2*extent+1,2*extent+1,2*extent+1,nsf1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to phi_k/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          phi_k = zero
          !write(*,*) '\phi_{k\gamma} allocated', shape(phi_k)

       case('phi_l') ! allocate phi_l for neighbour atom [Xrange]
          allocate(phi_l(2*extent+1,2*extent+1,2*extent+1,nsf1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to phi_l/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          phi_l = zero
          !write(*,*) '\phi_{l\delta} allocated', shape(phi_l)
          
       case('phi_l_m') ! allocate phi_l_m and save it
          allocate(phi_l_m(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to phi_l_m/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          phi_l_m = zero
          !write(unit,*) '\phi_{l\delta}_m allocated', shape(phi_l_m)
          
       case('phi_k_m') ! allocate phi_k_m and save it
          allocate(phi_k_m(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to phi_k_m/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          phi_k_m = zero
          !write(unit,*) '\phi_{k\delta}_m allocated', shape(phi_k_m)
          !
          !===auxiliarydensity=allocation===========================================================
       case('rho_kj') ! allocate rho_kj for neighbour[k]/neighbour[j] atoms [Krange/Srange]
          allocate(rho_kj(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to rho_kj/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          rho_kj = zero
          !write(unit,*) '\rho_{k\gamma j\beta} allocated'

       case('rho_lj') ! allocate rho_lj for neighbour[l]/neighbour[j] atoms [Krange/Srange]
          allocate(rho_lj(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to rho_lj/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          rho_lj = zero
          !write(unit,*) '\rho_{l\delta j\beta} allocated'

       case('rho_ik') ! allocate rho_ik for primary[i]/neighbour[k] atoms [Krange/Srange]
          allocate(rho_ik(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to rho_ik/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          rho_ik = zero
          !write(unit,*) '\rho_{i\alpha k\gamma} allocated'
          !
          !===other=allocation======================================================================
       case('Ome_i') ! allocate Ome_i 
          allocate(Ome_i(2*extent+1,2*extent+1,2*extent+1,nsf1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Ome_i/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          Ome_i = zero
          !write(unit,*) '\Ome_{i\alpha} allocated'

       case('Ome_j') ! allocate Ome_j 
          allocate(Ome_j(2*extent+1,2*extent+1,2*extent+1,nsf1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Ome_j/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          Ome_j = zero
          !write(unit,*) '\Ome_{j\beta} allocated'

       case('Ome_kj') ! allocate Ome_kj 
          allocate(Ome_kj(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Ome_j/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          Ome_kj = zero
          !write(unit,*) '\Ome_{k\gamma}_{j\beta} allocated'

       case('vhf_kj') ! allocate vhf_kj
          allocate(vhf_kj(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to vhf_kj/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          vhf_kj = zero
          !write(unit,*) '\vhf_{k\gamma j\beta} allocated'

       case('Phy_k')! allocate Phy_k
          allocate(Phy_k(2*extent+1,2*extent+1,2*extent+1,nsf1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Phy_j/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          Phy_k = zero
          !write(unit,*) '\Phy_{k\gamma} allocated'

       case('vhf_lj') ! allocate vhf_on_grid_lj
          allocate(vhf_lj(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to vhf_lj/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1), &
               type_dbl)!,matrix,lun)
          vhf_lj = zero
          !write(unit,*) '\vhf_{l\deta j\beta} allocated'
          
       case('ewald_3d') ! allocate ewald rho/pot
          allocate(ewald_rho(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to ewald_rho/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'ewald_rho',lun)
          ewald_rho = zero          
          allocate(ewald_pot(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to ewald_pot/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'ewald_pot',lun)
          ewald_pot = zero
          !write(unit,*) 'ewald_rho/ewald_pot allocated'

       case('work_3d') ! allocate work matrix
          allocate(work_in_3d(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to work_in_3d/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'work_in ',lun)
          work_in_3d = zero          
          allocate(work_out_3d(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to work_out_3d/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'work_out',lun)
          work_out_3d = zero
          !write(unit,*) 'work_in_3d/work_out_3d allocated'

       case('reckernel_3d') ! allocate the reciprocal space kernel for scaling
          allocate(reckernel_3d(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to reckernel_3d/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),type_dbl)!,matrix,lun)
          reckernel_3d = zero
          !write(unit,*) 'reckernel_3d allocated'
          
       case('isf_rho') ! allocate isf_rho
          allocate(isf_rho(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to isf_rho/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'isf_rho',lun)
          isf_rho = zero          
          allocate(isf_pot_ion(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to isf_pot_ion/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'isf_pot',lun)
          isf_pot_ion = zero
          
       case('exx_mat')
          allocate(Imat(nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Imat/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2,type_dbl)!,'Imat',lun)          
          Imat = reshape((/(1.0,(0.0,i=1,nsf1),j=1,(nsf2-1)),1.0/),(/nsf1,nsf2/))
          allocate(EXXmat(nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to EXXmat/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2,type_dbl)!,'Exxmat',lun)
          EXXmat = zero
          allocate(Kmat_ij(nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Kmat_ij/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2,type_dbl)!,'Kmatij',lun)
          Kmat_ij = zero
          !allocate(EXX_tmp(extent,extent,nsf1,nsf2), STAT=stat)
          !if(stat/=0) call cq_abort('Error allocating memory to EXX_tmp/exx !',stat)
          !call reg_alloc_mem(area_exx,extent*extent*nsf1*nsf2,type_dbl,'EXX_tmp',lun)
          !EXX_tmp = zero
          allocate(Kmat_tmp(nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Kmat_tmp/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2,type_dbl)!,'Kmat_tmp',lun)
          Kmat_tmp = zero
          !allocate(EXX_tmp2(extent,nsf1,extent,extent,nsf1,nsf1), STAT=stat)
          !if(stat/=0) call cq_abort('Error allocating memory to EXX_tmp2/exx !',stat)
          !call reg_alloc_mem(area_exx,extent*extent*extent*nsf1*nsf1*nsf1,type_dbl,'EXX_tmp',lun)
          !EXX_tmp2 = zero
          !write(unit,*) 'EXXmat/Kmat_ij/Kmat_tmp/EXX_tmp2 allocated'          
          
       case('exx_Kkl')
          allocate(Kmat_kl(extent,extent,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Kmat_kl/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*extent*extent,type_dbl)!,'Kmatkl',lun)
          Kmat_kl = zero
          !write(unit,*) 'Kmat_kl allocated'

          !===FFT=allocation========================================================================
       case('fftw_3d') 
          ! purpose: FFT[rho_{k,j}(r) => rho_{k,j}(G)
          ! allocate arrayin/arrayout [remember: type_cplx=[2!]*type_dbl (I guess ?!)]
          allocate(fftwrho3d%arrayin(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to fftwin3d/exx !',stat)
          call reg_alloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'fftw_in ',lun)
          allocate(fftwrho3d%arrayout(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to fftwout3d/exx !',stat)
          call reg_alloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'fftw_out',lun)
          !allocate(fftwrho3d%auxin(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          !if(stat/=0) call cq_abort('Error allocating memory to auxin/exx !',stat)
          !call reg_alloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),type_dbl,matrix,lun)
          !allocate(fftwrho3d%auxout(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          !if(stat/=0) call cq_abort('Error allocating memory to auxout/exx !',stat)
          !call reg_alloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),type_dbl,matrix,lun)
          ! purpose: create FFTW plans
          !
          call  fft3_init_wrapper( 2*extent+1  )
          !
          !call dfftw_plan_dft_3d(fftwrho3d%planF,(2*extent+1),(2*extent+1),(2*extent+1), &
          !     fftwrho3d%arrayin,fftwrho3d%arrayout,FFTW_FORWARD,FFTW_ESTIMATE)          
          !call dfftw_plan_dft_3d(fftwrho3d%planR,(2*extent+1),(2*extent+1),(2*extent+1), &
          !     fftwrho3d%arrayin,fftwrho3d%arrayout,FFTW_BACKWARD,FFTW_ESTIMATE)          
          fftwrho3d%arrayin  = cmplx(zero,zero,double_cplx)
          fftwrho3d%arrayout = cmplx(zero,zero,double_cplx)
          !fftwrho3d%auxin    = cmplx(zero,zero,double_cplx)
          !fftwrho3d%auxout   = cmplx(zero,zero,double_cplx)

       case('fftw_2d') 
          ! allocate arrayin/arrayout [remember: type_cplx=[2!]*type_dbl (I guess ?!)]
          allocate(fftwrho2d%arrayin(2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to fftwin2d/exx !',stat)
          call reg_alloc_mem(area_exx,2*(2*extent+1)*(2*extent+1),type_dbl)!,matrix,lun)
          allocate(fftwrho2d%arrayout(2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to fftwout2d/exx !',stat)
          call reg_alloc_mem(area_exx,2*(2*extent+1)*(2*extent+1),type_dbl)!,matrix,lun)
          ! purpose: create FFTW plans
          call dfftw_plan_dft_2d(fftwrho2d%planF,(2*extent+1),(2*extent+1), &
               fftwrho2d%arrayin,fftwrho2d%arrayout,FFTW_FORWARD,FFTW_ESTIMATE)          
          call dfftw_plan_dft_2d(fftwrho2d%planR,(2*extent+1),(2*extent+1), &
               fftwrho2d%arrayin,fftwrho2d%arrayout,FFTW_BACKWARD,FFTW_ESTIMATE)          
          fftwrho2d%arrayin  = cmplx(zero,zero,double_cplx)
          fftwrho2d%arrayout = cmplx(zero,zero,double_cplx)

       case('fftw_1d')
          ! allocate arrayin/arrayout [remember: type_cplx=[2!]*type_dbl (I guess ?!)]
          allocate(fftwrho1d%arrayin(2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to fftwin1d/exx !',stat)
          call reg_alloc_mem(area_exx,2*(2*extent+1),type_dbl)!,matrix,lun)
          allocate(fftwrho1d%arrayout(2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to fftwout1d/exx !',stat)
          call reg_alloc_mem(area_exx,2*(2*extent+1),type_dbl)!,matrix,lun)
          ! purpose: create FFTW plans
          call dfftw_plan_dft_1d(fftwrho1d%planF,(2*extent+1), &
               fftwrho1d%arrayin,fftwrho1d%arrayout,FFTW_FORWARD,FFTW_ESTIMATE)          
          call dfftw_plan_dft_1d(fftwrho1d%planR,(2*extent+1), &
               fftwrho1d%arrayin,fftwrho1d%arrayout,FFTW_BACKWARD,FFTW_ESTIMATE)          
          fftwrho1d%arrayin  = cmplx(zero,zero,double_cplx)
          fftwrho1d%arrayout = cmplx(zero,zero,double_cplx)

       case default
          call cq_abort('Error allocation/exx_mem !',stat)          

       end select

       call stop_timer(tmr_std_exx_allocat,accumulate)
    else       
       call start_timer(tmr_std_exx_dealloc)
       select case (matrix)
          !===PAO=deallocation======================================================================
       case('phi_i')
          deallocate(phi_i,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to phi_i/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(*,*) '\phi_{i\alpha} deallocated'

       case('phi_j')
          deallocate(phi_j,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to phi_j/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(*,*) '\phi_{j\beta} deallocated'

       case('phi_k')
          deallocate(phi_k,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to phi_k/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(*,*) '\phi_{k\gamma} deallocated'

       case('phi_l')
          deallocate(phi_l,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to phi_l/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(*,*) '\phi_{l\delta} deallocated'
          
       case('phi_l_m')
          deallocate(phi_l_m,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to phi_l_m/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(unit,*) '\phi_{l\delta}_m deallocated'

       case('phi_k_m')
          deallocate(phi_k_m,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to phi_k_m/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(unit,*) '\phi_{k\delta}_m deallocated'

          !===auxiliarydensity=deallocation=========================================================
       case('rho_kj') ! allocate rho_kj for neighbour[k]/neighbour[j] atoms [Krange/Srange]
          deallocate(rho_kj, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to rho_kj/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(unit,*) '\rho_{k\gamma j\beta} deallocated'

       case('rho_lj') ! allocate rho_lj for neighbour[l]/neighbour[j] atoms [Krange/Srange]
          deallocate(rho_lj, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to rho_lj/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(unit,*) '\rho_{l\delta j\beta} deallocated'
          
       case('rho_ik') ! allocate rho_ik for primary[i]/neighbour[j] atoms [Krange/Srange]
          deallocate(rho_ik, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to rho_ik/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(unit,*) '\rho_{i\delta k\gamma} deallocated'          
          !===other=deallocation====================================================================
       case('Ome_i')
          deallocate(Ome_i,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to Ome_i/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(unit,*) '\Ome_{i\alpha} deallocated'

       case('Ome_j')
          deallocate(Ome_j,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to Ome_j/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(unit,*) '\Ome_{j\beta} deallocated'

       case('Ome_kj')
          deallocate(Ome_kj,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to Ome_j/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(unit,*) '\Ome_{k\gamma}_{j\beta} deallocated'

       case('vhf_kj')
          deallocate(vhf_kj, STAT=stat)          
          if(stat/=0) call cq_abort('Error deallocating memory to vhf_kj/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(unit,*) '\vhf_{k\gamma j\beta} deallocated'
          
       case('Phy_k')
          deallocate(Phy_k, STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Phy_k/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(unit,*) '\Phy_{k\gamma} deallocated'

       case('vhf_lj')
          deallocate(vhf_lj, STAT=stat)          
          if(stat/=0) call cq_abort('Error deallocating memory to vhf_lj/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)
          !write(unit,*) '\vhf_{l\delta j\beta} deallocated'

       case('ewald_3d')
          deallocate(ewald_rho, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to ewald_rho/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'ewald_rho',lun)
          deallocate(ewald_pot, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to ewald_pot/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'ewald_pot',lun)
          !write(unit,*) 'ewald_rho/ewald_pot deallocated'
          
       case('work_3d') 
          deallocate(work_in_3d, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to work_in_3d/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'work_in ',lun)
          deallocate(work_out_3d, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to work_out_3d/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'work_out',lun)
          !write(unit,*) 'work_in_3d/work_out_3d deallocated'          
          
       case('reckernel_3d') 
          deallocate(reckernel_3d, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to reckernel_3d/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,matrix,lun)          
          !write(unit,*) 'reckernel_3d deallocated'

       case('isf_rho') ! allocate isf_rho
          deallocate(isf_rho, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to isf_rho/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'isf_rho',lun)          
          deallocate(isf_pot_ion, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to isf_pot_ion/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'isf_pot',lun)          
          !write(unit,*) 'isf_rho/isf_pot_ion deallocated'

       case('exx_mat') ! deallocate exx matrix
          deallocate(Imat, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to Imat/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2,type_dbl)!,'Imat',lun)
          deallocate(EXXmat, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to EXXmat/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2,type_dbl)!,'Exxmat',lun)
          deallocate(Kmat_ij, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to Kmat_ij/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2,type_dbl)!,'Kmatij',lun)
          !deallocate(EXX_tmp, STAT=stat)
          !if(stat/=0) call cq_abort('Error deallocating memory to EXX_tmp/exx !',stat)
          !call reg_dealloc_mem(area_exx,1000*1000*nsf1*nsf2,type_dbl,'EXX_tmp',lun)
          deallocate(Kmat_tmp, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to Kmat_tmp/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2,type_dbl)!,'Kmat_tmp',lun)
          !deallocate(EXX_tmp2, STAT=stat)
          !if(stat/=0) call cq_abort('Error deallocating memory to EXX_tmp2/exx !',stat)
          !call reg_dealloc_mem(area_exx,extent*extent*extent*nsf1*nsf1*nsf1,type_dbl,'EXX_tmp2',lun)
          !write(unit,*) 'EXXmat/Kmat_ij/Kmat_tmp/EXX_tmp2 allocated'

       case('exx_Kkl') ! deallocate exx matrix
          deallocate(Kmat_kl, STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Kmat_kl/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*extent*extent,type_dbl)!,'Kmatkl',lun)
          !write(unit,*) 'Kmat_kl allocated'

          !===FFT=deallocation======================================================================
       case('fftw_3d') 

          !call dfftw_destroy_plan(fftwrho3d%planF)
          !call dfftw_destroy_plan(fftwrho3d%planR)
          !call fft3_dest_wrapper( )
          !call fft3_dest_wrapper( )

          deallocate(fftwrho3d%arrayin, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to fftwin3d/exx !',stat)
          call reg_dealloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'fftw_in ',lun)        
          deallocate(fftwrho3d%arrayout, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to fftwout3d/exx !',stat)
          call reg_dealloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl)!,'fftw_out',lun)
          !deallocate(fftwrho3d%auxout, STAT=stat)
          !if(stat/=0) call cq_abort('Error deallocating memory to fftw_auxout/exx !',stat)
          !call reg_dealloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),type_dbl)        
          !deallocate(fftwrho3d%auxin, STAT=stat)
          !if(stat/=0) call cq_abort('Error deallocating memory to fftw_auxin/exx !',stat)
          !call reg_dealloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),type_dbl)
          
       case('fftw_2d') 
          call dfftw_destroy_plan(fftwrho2d%planF)
          call dfftw_destroy_plan(fftwrho2d%planR)
          deallocate(fftwrho2d%arrayin, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to fftwin2d/exx !',stat)
          call reg_dealloc_mem(area_exx,2*(2*extent+1)*(2*extent+1),type_dbl)        
          deallocate(fftwrho2d%arrayout, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to fftwout2d/exx !',stat)
          call reg_dealloc_mem(area_exx,2*(2*extent+1)*(2*extent+1),type_dbl)       
          
       case('fftw_1d') 
          call dfftw_destroy_plan(fftwrho1d%planF)
          call dfftw_destroy_plan(fftwrho1d%planR)
          deallocate(fftwrho1d%arrayin, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to fftwin1d/exx !',stat)
          call reg_dealloc_mem(area_exx,2*(2*extent+1),type_dbl)        
          deallocate(fftwrho1d%arrayout, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to fftwout1d/exx !',stat)
          call reg_dealloc_mem(area_exx,2*(2*extent+1),type_dbl)       

       case default
          call cq_abort('Error deallocation/exx_mem !',stat)          

       end select
       call stop_timer(tmr_std_exx_dealloc,accumulate)
    end if
    
    return
  end subroutine exx_mem_alloc

end module exx_memory
