! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id: $
! -----------------------------------------------------------
! Module exx_memory.f90
! -----------------------------------------------------------
! Code area 13: EXX
! -----------------------------------------------------------

!!***h* Conquest/exx_memory *
!!  NAME
!!   exx_memory
!!  CREATION DATE
!!   2011/02/11
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
module exx_memory

  use datatypes

  use exx_types,         ONLY: phi_i, phi_j,  phi_k,  phi_l
  use exx_types,         ONLY: Phy_k, Ome_kj, rho_kj, vhf_kj
  use exx_types,         ONLY: rho_ki, vhf_lj

  
  use exx_types,         ONLY: work_in_3d,work_out_3d 
  use exx_types,         ONLY: reckernel_3d,ewald_rho,ewald_pot
  use exx_types,         ONLY: isf_rho,isf_pot_ion

  use exx_types,         ONLY: fftwrho3d

  use exx_types,         ONLY: eris
  
  use exx_types,         ONLY: tmr_std_exx_allocat, tmr_std_exx_dealloc
  use exx_types,         ONLY: unit_memory_write

  use global_module,             ONLY: area_exx, io_lun
  use timer_module,              ONLY: start_timer, stop_timer, print_timer
  use memory_module,             ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl, type_cplx
  use GenComms,                  ONLY: cq_abort

  use fft_interface_module,      ONLY: fft3_init_wrapper
  !**<lat>** don't touch
  !use fftw_module
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

    !print*, "EXX: exx_memory lu =", lun
    
!!$
!!$
!!$ ALLOCATION for EXX
!!$
!!$
    if(flag=='alloc') then
       call start_timer(tmr_std_exx_allocat)       
       !
       select case (matrix)          
          !
!!$
!!$
!!$ for PAO grids allocations
!!$
!!$
          !
       case('phi_i') ! allocate phi_i for primary atom
          allocate(phi_i(2*extent+1,2*extent+1,2*extent+1,nsf1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to phi_i/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          phi_i = zero
          !write(*,*) '\phi_{i\alpha} allocated', shape(phi_i)
          !
          !
       case('phi_j') ! allocate phi_j for neighbour atom [Srange]
          allocate(phi_j(2*extent+1,2*extent+1,2*extent+1,nsf1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to phi_j/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          phi_j = zero
          !write(*,*) '\phi_{i\beta} allocated', shape(phi_j)
          !
          !
       case('phi_k') ! allocate phi_k for neighbour atom [Xrange]
          allocate(phi_k(2*extent+1,2*extent+1,2*extent+1,nsf1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to phi_k/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          phi_k = zero
          !write(*,*) '\phi_{k\gamma} allocated', shape(phi_k)
          !
          !
       case('phi_l') ! allocate phi_l for neighbour atom [Xrange]
          allocate(phi_l(2*extent+1,2*extent+1,2*extent+1,nsf1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to phi_l/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          phi_l = zero
          !write(*,*) '\phi_{l\delta} allocated', shape(phi_l)
          !
          !
!!$
!!$
!!$ for auxiliary density and potential grid allocations
!!$
!!$
          !
       case('rho_kj') ! allocate rho_kj for neighbour[k]/neighbour[j] atoms [Krange/Srange]
          allocate(rho_kj(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to rho_kj/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          rho_kj = zero
          !write(unit,*) '\rho_{k\gamma j\beta} allocated'
          !
          !
       case('Ome_kj') ! allocate Ome_kj 
          allocate(Ome_kj(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Ome_j/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          Ome_kj = zero
          !write(unit,*) '\Ome_{k\gamma}_{j\beta} allocated'
          !
          !
       case('vhf_kj') ! allocate vhf_kj
          allocate(vhf_kj(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to vhf_kj/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          vhf_kj = zero
          !write(unit,*) '\vhf_{k\gamma j\beta} allocated'
          !
          !
       case('Phy_k')! allocate Phy_k
          allocate(Phy_k(2*extent+1,2*extent+1,2*extent+1,nsf1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Phy_k/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          Phy_k = zero
          !write(unit,*) '\Phy_{k\gamma} allocated'
          !
          !
       case('rho_ki') ! allocate rho_ki for neighbour[k]/neighbour[i] atoms [Krange/Srange]
          allocate(rho_ki(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to rho_ki/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          rho_ki = zero
          !write(unit,*) '\rho_{k\gamma i\alpha} allocated'
          !
          !
      case('vhf_lj') ! allocate vhf_lj
          allocate(vhf_lj(2*extent+1,2*extent+1,2*extent+1,nsf1,nsf2), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to vhf_lj/exx !',stat)
          call reg_alloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          vhf_lj = zero
          !write(unit,*) '\vhf_{l\delta j\beta} allocated'          

!!$
!!$
!!$ for working grids allocations
!!$        (related to FFT)
!!$
          !
       case('ewald_3d') ! allocate ewald rho/pot
          allocate(ewald_rho(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to ewald_rho/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'ewald_rho',lun)
          ewald_rho = zero          
          allocate(ewald_pot(2*extent+1,2*extent+1,2*extent+1), STAT=stat) 
          if(stat/=0) call cq_abort('Error allocglobating memory to ewald_pot/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'ewald_pot',lun)
          ewald_pot = zero
          !write(unit,*) 'ewald_rho/ewald_pot allocated'
          !
          !
       case('work_3d') ! allocate work matrix
          allocate(work_in_3d(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to work_in_3d/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'work_in',lun)
          work_in_3d = zero          
          allocate(work_out_3d(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to work_out_3d/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'work_out',lun)
          work_out_3d = zero
          !write(unit,*) 'work_in_3d/work_out_3d allocated'
          !
          !
       case('reckernel_3d') ! allocate the reciprocal space kernel for scaling
          allocate(reckernel_3d(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to reckernel_3d/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_cplx,matrix,lun)
          reckernel_3d = zero
          !write(unit,*) 'reckernel_3d allocated'
          !
          !
       case('isf_rho') ! allocate isf_rho
          allocate(isf_rho(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to isf_rho/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'isf_rho',lun)
          isf_rho = zero          
          allocate(isf_pot_ion(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to isf_pot_ion/exx !',stat)
          call reg_alloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'isf_pot',lun)
          isf_pot_ion = zero
          !
          !
!!$
!!$
!!$ FFT allocations
!!$
!!$
          !
       case('fftw_3d') 
          ! purpose: FFT[rho_{k,j}(r) => rho_{k,j}(G)
          ! allocate arrayin/arrayout [remember: type_cplx=[2!]*type_dbl (I guess ?!)]
          allocate(fftwrho3d%arrayin(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to fftwin3d/exx !',stat)
          call reg_alloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_cplx,'fftw_in ',lun)
          allocate(fftwrho3d%arrayout(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to fftwout3d/exx !',stat)
          call reg_alloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_cplx,'fftw_out',lun)
          !
          !**<lat>** don't touch
          !
          !allocate(fftwrho3d%auxin(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          !if(stat/=0) call cq_abort('Error allocating memory to auxin/exx !',stat)
          !call reg_alloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),type_dbl,matrix,lun)
          !allocate(fftwrho3d%auxout(2*extent+1,2*extent+1,2*extent+1), STAT=stat)
          !if(stat/=0) call cq_abort('Error allocating memory to auxout/exx !',stat)
          !call reg_alloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),type_dbl,matrix,lun)
          !
          !**<lat>**  purpose: create FFTW plans
          !
          !call dfftw_plan_dft_3d(fftwrho3d%planF,(2*extent+1),(2*extent+1),(2*extent+1), &
          !     fftwrho3d%arrayin,fftwrho3d%arrayout,FFTW_FORWARD,FFTW_ESTIMATE)          
          !call dfftw_plan_dft_3d(fftwrho3d%planR,(2*extent+1),(2*extent+1),(2*extent+1), &
          !     fftwrho3d%arrayin,fftwrho3d%arrayout,FFTW_BACKWARD,FFTW_ESTIMATE)          
          !
          call fft3_init_wrapper( 2*extent+1  )
          !
          fftwrho3d%arrayin  = cmplx(zero,zero,double_cplx)
          fftwrho3d%arrayout = cmplx(zero,zero,double_cplx)
          !
!!$
!!$
!!$ ERI allocation
!!$
!!$
       !case('eris') ! allocate electron repulsion integrals
       !   allocate(eris(extent), STAT=stat)
       !   if(stat/=0) call cq_abort('Error allocating memory to eris/exx !',stat)
       !   call reg_alloc_mem(area_exx,extent,type_dbl,matrix,lun)
       !   eris = zero
          !write(*,*) 'ERIs allocated', shape(eris)
          !
       case default
          call cq_abort('Error allocation/exx_mem !',stat)          
          !
          !
       end select

       call stop_timer(tmr_std_exx_allocat,accumulate)
    else       
!!$
!!$
!!$ DEALLOCATION for EXX
!!$
!!$
       call start_timer(tmr_std_exx_dealloc)
       select case (matrix)
       case('phi_i')
          deallocate(phi_i,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to phi_i/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          !write(*,*) '\phi_{i\alpha} deallocated'
          !
          !
       case('phi_j')
          deallocate(phi_j,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to phi_j/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          !write(*,*) '\phi_{j\beta} deallocated'
          !
          !
       case('phi_k')
          deallocate(phi_k,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to phi_k/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          !write(*,*) '\phi_{k\gamma} deallocated'
          !
          !
       case('phi_l')
          deallocate(phi_l,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to phi_l/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          !write(*,*) '\phi_{l\delta} deallocated'
          !
          !
       case('rho_kj') ! allocate rho_kj for neighbour[k]/neighbour[j] atoms [Krange/Srange]
          deallocate(rho_kj, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to rho_kj/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          !write(unit,*) '\rho_{k\gamma j\beta} deallocated'
          !
          !
       case('Ome_kj')
          deallocate(Ome_kj,STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to Ome_j/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          !write(unit,*) '\Ome_{k\gamma}_{j\beta} deallocated'

       case('vhf_kj')
          deallocate(vhf_kj, STAT=stat)          
          if(stat/=0) call cq_abort('Error deallocating memory to vhf_kj/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          !write(unit,*) '\vhf_{k\gamma j\beta} deallocated'
          !
          !
       case('Phy_k')
          deallocate(Phy_k, STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory to Phy_k/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          !write(unit,*) '\Phy_{k\gamma} deallocated'
          !
          !
       case('rho_ki') ! allocate rho_ki for neighbour[k]/neighbour[i] atoms [Krange/Srange]
          deallocate(rho_ki, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to rho_ki/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          !write(unit,*) '\rho_{k\gamma i\alpha} deallocated'
          !
          !
       case('vhf_lj')
          deallocate(vhf_lj, STAT=stat)          
          if(stat/=0) call cq_abort('Error deallocating memory to vhf_lj/exx !',stat)
          call reg_dealloc_mem(area_exx,nsf1*nsf2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)
          !write(unit,*) '\vhf_{l\delta j\beta} deallocated'
          !
          !          
       case('ewald_3d')
          deallocate(ewald_rho, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to ewald_rho/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'ewald_rho',lun)
          deallocate(ewald_pot, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to ewald_pot/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'ewald_pot',lun)
          !write(unit,*) 'ewald_rho/ewald_pot deallocated'
          !
          !
       case('work_3d') 
          deallocate(work_in_3d, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to work_in_3d/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'work_in ',lun)
          deallocate(work_out_3d, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to work_out_3d/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'work_out',lun)
          !write(unit,*) 'work_in_3d/work_out_3d deallocated'          
          !
          !
       case('reckernel_3d') 
          deallocate(reckernel_3d, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to reckernel_3d/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,matrix,lun)          
          !write(unit,*) 'reckernel_3d deallocated'
          !
          !
       case('isf_rho') ! allocate isf_rho
          deallocate(isf_rho, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to isf_rho/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'isf_rho',lun)          
          deallocate(isf_pot_ion, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to isf_pot_ion/exx !',stat)
          call reg_dealloc_mem(area_exx,(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'isf_pot',lun)          
          !write(unit,*) 'isf_rho/isf_pot_ion deallocated'
          !
          !
       case('fftw_3d')
          !call dfftw_destroy_plan(fftwrho3d%planF)
          !call dfftw_destroy_plan(fftwrho3d%planR)
          !call fft3_dest_wrapper( )
          !call fft3_dest_wrapper( )
          deallocate(fftwrho3d%arrayin, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to fftwin3d/exx !',stat)
          call reg_dealloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'fftw_in ',lun)        
          deallocate(fftwrho3d%arrayout, STAT=stat)
          if(stat/=0) call cq_abort('Error deallocating memory to fftwout3d/exx !',stat)
          call reg_dealloc_mem(area_exx,2*(2*extent+1)*(2*extent+1)*(2*extent+1),&
               type_dbl,'fftw_out',lun)
          !
       case('eris')
          deallocate(eris, STAT=stat)          
          if(stat/=0) call cq_abort('Error deallocating memory to eris/exx !',stat)
          call reg_dealloc_mem(area_exx,extent,type_dbl,matrix,lun)
          !write(unit,*) 'ERIs deallocated'  
          !
          !
       case default
          call cq_abort('Error deallocation/exx_mem !',stat)          

       end select
       call stop_timer(tmr_std_exx_dealloc,accumulate)
    end if
    !
    return
  end subroutine exx_mem_alloc
!!$
!!$
!!$
!!$
end module exx_memory
