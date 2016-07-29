! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: exx_module.f90 XX year-mm-dd XX:XX:XXZ Lionel $
! -----------------------------------------------------------
! Module exx_module
! -----------------------------------------------------------
! Code area 13: EXX
! -----------------------------------------------------------

!!***h* Conquest/exx_module *
!!  NAME
!!   exx_module
!!
!!  PURPOSE
!!   Holds routines called in exx_kernel_module
!!   to compute exact exchange
!!
!!  USES
!!   GenComms, dimens, datatypes, numbers, units 
!!   group_module, cover_module, global_module
!!   species_module, primary_module, matrix_data
!!   atomic_density, support_spec_format
!!   fft_interface_module, exx_types, exx_io 
!!   DiagModule
!!
!!  AUTHOR
!!   L.A.Truflandier (lat)
!!
!!  CREATION DATE
!!   2011/02/11
!!
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
module exx_module

  use datatypes
  use GenComms,                  only: my_barrier, cq_abort, mtime
  use GenComms,                  only: inode, ionode, myid, root
  use timer_module,              only: start_backtrace, stop_backtrace, cq_timer
  !**<lat>** ISF Poisson solver Will be available in the forthcoming version 
  !use Poisson_Solver,            only: PSolver, createKernel, gequad

  use exx_types,                 only: reckernel_3d, fftwrho3d, exx_debug
  use exx_io

  use fft_interface_module,      only: fft3_exec_wrapper, fft3_init_wrapper

  implicit none 
  
  ! Area identification
  integer, parameter, private :: area = 13

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = &
       "$Id: exx_module.f90 XX year-mm-dd XX:XX:XXZ lionel $"

  !!***

contains

  subroutine exx_scal_rho_3d(inode,extent,r_int,scheme,cutoff,omega,n_gauss, &
                             p_gauss,w_gauss)

    use datatypes
    use numbers,           ONLY: twopi, pi, sqrt_pi, very_small
    use numbers,           ONLY: zero, one, two, four, three_halves
     
    implicit none

    ! << Passed variables >>    
    character(*), intent(in) :: scheme
    integer,      intent(in) :: extent
    integer,      intent(in) :: inode
    real(double), intent(in) :: r_int

    real(double), intent(in) :: cutoff
    real(double), intent(in) :: omega

    integer,      intent(in) :: n_gauss
    real(double), optional, intent(in) :: p_gauss(n_gauss)    
    real(double), optional, intent(in) :: w_gauss(n_gauss)    
  
    ! << Local variables >>
    real(double),dimension(3) :: G
    real(double)              :: charge_G
    
    real(double), parameter   :: magic_number = 0.1864d0 ! 32/43 = nombre magic !

    real(double)              :: r_G, r_A, v_G, v_R, sphere, fftwnorm, tmp, alpha
    real(double)              :: rho_zero, Gkernel, Gkernel_zero, arg, factor, delta    
    real(double)              :: max_, min_

    integer                   :: gauss
    integer                   :: nx, ny, nz, ngrid   
    integer                   ::  x,  y,  z    
    integer                   ::  i,  j,  k

    ngrid     = 2*extent+1
    fftwnorm  = sqrt(real(ngrid**3,double)) 

    r_A    = real(2,double)*r_int
    r_G    = zero
    G      = zero

    v_G    = (twopi**3)*real(1,double)/(real(2,double)*r_int)**3
    v_R    = (real(2,double)*r_int)**3
    sphere = ((twopi*real(1,double)/(real(2,double)*r_int))**3)*(1.0d0/6.0d0*pi)

    reckernel_3d = zero

    poisson_scheme: select case(scheme)       
    case('default')
       delta  = zero
    case('ewald')
       delta  = zero
    case('pulay')
       delta  = zero
    case('yukawa')
       delta  = omega
    case('gauss')
       delta  = zero
    case default
       delta  = zero
    end select poisson_scheme

    grid_x_loop: do nx = -extent, extent
       x = nx + extent + 1
       if (nx <= 0) then
          G(1) = twopi*( real(nx+extent+delta,double))/r_A
       else 
          G(1) = twopi*(-real(extent-nx+(one-delta),double))/r_A
       end if
       
       grid_y_loop: do ny = -extent, extent
          y  = ny + extent + 1 
          if (ny <= 0) then                         
             G(2) = twopi*( real(ny+extent+delta,double))/r_A
          else
             G(2) = twopi*(-real(extent-ny+(one-delta),double))/r_A
          end if

          grid_z_loop: do nz = -extent, extent             
             z  = nz + extent + 1             
             if (nz <= 0) then
                G(3) = twopi*( real(nz+extent+delta,double))/r_A 
             else                
                G(3) = twopi*(-real(extent-nz+(one-delta),double))/r_A
             end if
             
             r_G = sqrt(dot_product(G,G))                         

             !if (scheme == 'yukawa') then
             !   if (r_G > 0.1d0) then
             !      delta = zero
             !   end if
             !end if

             if (r_G == zero) then                

                !write(*,*) x, y, z
                poisson_zero: select case(scheme)
                case('default')
                   reckernel_3d(x,y,z) = cmplx(0,0,double_cplx)          
                case('ewald')
                   reckernel_3d(x,y,z) = cmplx(0,0,double_cplx)          
                case('pulay')
                   reckernel_3d(x,y,z) = cmplx(cutoff**2/two,0,double_cplx)     
                case('yukawa')
                   reckernel_3d(x,y,z) = cmplx(0,0,double_cplx) 
                case('gauss')
                   Gkernel_zero = zero
                   do gauss = 1, n_gauss
                      Gkernel_zero = Gkernel_zero + &
                           w_gauss(gauss)*sqrt_pi/(four*p_gauss(gauss)**three_halves)
                   end do
                   reckernel_3d(x,y,z) = cmplx(Gkernel_zero,double_cplx)          
                case default 
                   reckernel_3d(x,y,z) = cmplx(zero,zero,double_cplx)          
                end select poisson_zero
             else
                poisson_xyz: select case(scheme)
                case('default')
                   reckernel_3d(x,y,z) = cmplx(one/r_G**2,zero,double_cplx)
                case('ewald')
                   reckernel_3d(x,y,z) = cmplx(one/r_G**2,zero,double_cplx)
                case('pulay')
                   reckernel_3d(x,y,z) = (one - cos(r_G*cutoff))/r_G**2
                case('yukawa')
                   reckernel_3d(x,y,z) = cmplx(one/r_G**2,zero,double_cplx)
                case('gauss')
                   Gkernel = zero
                   arg     = zero
                   factor  = zero

                   do gauss = 1, n_gauss                      
                      factor  = w_gauss(gauss)*sqrt(pi)/(four*p_gauss(gauss)**three_halves)
                      arg     = -(r_G/(two*sqrt(p_gauss(gauss))))**two
                      Gkernel = Gkernel + factor*exp(arg)
                      max_ = max(max_,arg)
                      min_ = min(min_,arg)

                   end do
                   reckernel_3d(x,y,z) = cmplx(Gkernel,zero,double_cplx)
                case default
                   reckernel_3d(x,y,z) = cmplx(one/r_G**2,zero,double_cplx)
                end select poisson_xyz
             end if

          end do grid_z_loop
       end do grid_y_loop
    end do grid_x_loop    

    return
  end subroutine exx_scal_rho_3d  
  !
  !
  !
  !
  subroutine exx_v_on_grid(inode,extent,rho,potential,r_int,poisson,scheme,&
       alpha,omega,n_gauss,p_gauss,w_gauss)
    
    use numbers,   ONLY: zero, one, fourpi
    use exx_types, ONLY: reckernel_3d                 ! FFTW
    use exx_types, ONLY: kernel, isf_rho, isf_pot_ion ! ISF  
    use exx_types, ONLY: pulay_radius, ewald_rho, ewald_pot, ewald_charge
    use exx_types, ONLY: tmr_std_exx_poisson

    implicit none

    ! << Input variables >>
    integer,      intent(in) :: inode
    integer,      intent(in) :: extent 
    real(double), intent(in) :: r_int
    character(*), intent(in) :: poisson    

    real(double), dimension(2*extent+1,2*extent+1,2*extent+1), &
         intent(in)  :: rho
    
    character(*), intent(in) :: scheme
    real(double), intent(in) :: alpha
    real(double), intent(in) :: omega        
    integer,      intent(in) :: n_gauss
    real(double), dimension(n_gauss) :: p_gauss, w_gauss

    ! << Output variables >>
    real(double), dimension(2*extent+1,2*extent+1,2*extent+1), &
         intent(out) :: potential
    
    ! << Local variables >>    
    integer       :: ng, i, j, k
    real(double)  :: fftwnorm, normr
    real(double)  :: grid_spacing
    real(double)  :: r(3)
    logical, parameter :: accumulate = .true.

    ! ... For Poisson/ISF real space FFT >>
    real(kind=8)          :: isf_eh, isf_exc, isf_vxc, isf_spacing

    ! ... Dummy variables >>
    real(double), dimension(:,:,:), allocatable :: dum_pot_ion, dum_rho

    ng           = 2*extent+1
    fftwnorm     = sqrt(real(ng**3,double))   
    grid_spacing = r_int/real(extent,double)

    !call start_timer(tmr_std_exx_poisson)    

    select case(poisson)

    case('fftw')       
       potential = zero

       ! setup[rho(r)] 
       fftwrho3d%arrayin  = cmplx(rho,zero,double_cplx)

       ! FFT_F[rho(r)] => rho(G)      
       call fft3_exec_wrapper( fftwrho3d%arrayin, ng , +1 )

       ! scale[rho_(G)] = 4pi*rho(G)/|G|^2
       fftwrho3d%arrayin = fourpi*fftwrho3d%arrayin*reckernel_3d
       
       ! FFT_B[4pi*rho(G)/|G|^2] = V(r')
       call fft3_exec_wrapper( fftwrho3d%arrayin, ng , -1 )
       
       ! Normalization
       potential = real(fftwrho3d%arrayin) / fftwnorm**2

    case('isf')       
       !
       call cq_abort('EXX: not available yet')
       !
       !isf_rho     = rho
       !potential   = zero
       !isf_pot_ion = zero
       !
       !call PSolver('F','G',0,1,ng,ng,ng,0,grid_spacing,grid_spacing,grid_spacing, &
       !     isf_rho,kernel,isf_pot_ion,isf_eh,isf_exc,isf_vxc,zero,.false.,1)       
       !potential = isf_rho
       !
    end select

    !call stop_timer(tmr_std_exx_poisson,accumulate)    

    return
  end subroutine exx_v_on_grid
  !
  !
  subroutine get_X_params(level)
    
    use exx_types, only: exx_psolver,p_scheme,p_scheme_default,p_cutoff,p_factor,p_omega
    use exx_types, only: pulay_factor,pulay_radius, magic_number,ewald_alpha,isf_order
    use exx_types, only: extent,ngrid,r_int,grid_spacing,volume,edge,screen

    use exx_types, only: exx_scheme, exx_phil, exx_mem, exx_screen, exx_alloc, exx_cutoff
    use exx_types, only: exx_cartesian, exx_overlap, exx_radius, exx_screen_pao, exx_hgrid
    use exx_types, only: unit_output_write, exx_phik, exx_gto, exx_debug
    use exx_types, only: tmr_std_exx_setup

    use atomic_density, only: atomic_density_table
    use species_module, only: n_species

    use global_module,  only: iprint_exx, exx_hgrid_medium, exx_hgrid_coarse, exx_alpha
    use numbers,        only: zero, very_small, one, two, three
    use units,          only: BohrToAng
    use dimens,         only: r_super_x, r_super_y, r_super_z, &
                              n_grid_x, n_grid_y, n_grid_z,    &
                              r_h, r_nl, r_h   
    implicit none

    !real(double), optional, intent(in) :: exx_scf_hgrid
    integer, optional :: level

    character(100) :: solver, scheme
    character(100) :: phil_scheme, eri_scheme, alloc_scheme, pao_scheme
    character(100) ::  mem_scheme, scr_scheme
    character(100) :: filename
    real(double)   :: threshold = 1.0d-4
    real(double)   :: r_shift   = 0.0_double
    real(double)   :: gs(3)     = 0.0_double
    real(double)   :: r_max, r_maxs, gs_min
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level
    integer        :: i,j,k,l,unit
   
    !
!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='get_X_params',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$
    !
    call start_timer(tmr_std_exx_setup) 
    !
    ! **<lat>** This routine need to be re-written !
    !
    exx_phil        = .false.
    exx_phik        = .false.
    exx_screen      = .false.
    exx_screen_pao  = .false.
    exx_gto         = .false.
    exx_cutoff      = zero              ! do not touch  
    unit            = unit_output_write ! do not touch  
    !if (exx_mem == 2) exx_alloc = .false.  
    !
    ! Find out the finest grid spacing from input
    ! Input grid spacing from input (Bohr unit)
    gs(1) = r_super_x/n_grid_x
    gs(2) = r_super_y/n_grid_y
    gs(3) = r_super_z/n_grid_z
    !
    ! Find out the finest grid spacing from input
    ! we use the same grid spacing for Ox,Oy,Oz   
    gs_min = minval(gs)
    !
    ! Setup grid spacing of EXX
    if ( exx_hgrid < very_small .and. exx_hgrid >= zero ) then
       grid_spacing = gs_min
    else if ( exx_hgrid >= zero ) then
       grid_spacing = exx_hgrid       
    else
       call cq_abort('EXX: unrecognised grid_spacing ',exx_hgrid)
       !
    end if
    !
    !
    ! **<lat>** Need to be adapted to GPFA ; not too difficult
    ! Find r_max for integration
    r_max = zero
    do i = 1, n_species
       r_max = max(r_max,atomic_density_table(i)%cutoff)
    end do
    !
    ! Setup EXX r_int
    if (exx_radius > very_small) then
       r_int = exx_radius
    else
       r_int = real(ceiling(r_max),double)
    end if
    !
    ! Number of grid points alongthe ]O,X] segment
    extent = ceiling(r_int/grid_spacing)    
    !
    ! Re-compute the real grid spacing with respect to extent
    grid_spacing = r_int/real(extent,double)
    !
    ! Number of grid points along [-X,X] segment
    ngrid  = 2*extent+1
    edge   = two*r_int 
    volume = edge**3    
    !
    ! Setup Poisson solver
    if (exx_psolver == 'fftw') then 
       solver = 'reciprocal space'
       scheme = trim(exx_psolver)     
    else 
       exx_psolver = 'isf'
       solver   = 'real space'
       scheme   = trim(exx_psolver)
       if (isf_order <= 0) then
          isf_order = 8
       end if
    end if
    !
    ! Below for output purpose
    if (exx_scheme==1) then
       eri_scheme = "3center reduction integrals"
    else
       eri_scheme = "4center eris"
    end if
    !
    if (exx_phil) then
       phil_scheme = 'storage'
    else
       phil_scheme = 'on-the-fly'
    end if
    !
    if (exx_alloc) then
       alloc_scheme = 'on-the-fly'
    else
       alloc_scheme = 'once'
    end if
    !
    if (exx_cartesian) then
       pao_scheme = 'cartesian'
    else
       pao_scheme = 'spherical'
    end if
    !
    if (exx_mem == 1) then
       mem_scheme = 'high'
    else if (exx_mem == 2) then
       mem_scheme = 'low'
    else if (exx_mem == 3) then
       mem_scheme = 'medium'
    else
       mem_scheme = 'low'
       exx_mem    = 2
    end if
    !
    if (exx_screen) then
       if (exx_cutoff > very_small) then
          scr_scheme     = 'user'
          exx_screen_pao = .false.
          screen         = exx_cutoff
       else
          exx_screen_pao =.true.
          scr_scheme = 'pao'
       end if
    else
       exx_screen_pao = .false.
       scr_scheme = 'n/a'
       screen     = 100.0d0
    end if
    !
    if (exx_scheme==2) then
       phil_scheme  = 'on-the-fly'
       alloc_scheme = 'on-the-fly'
       mem_scheme   = 'high'
    end if
    !
    if ( inode == ionode .and. exx_debug ) then
       write(unit,2) ('Entering in the Hartree-Fock module')
       write(unit,41) eri_scheme
       write(unit,42) phil_scheme
       write(unit,43) alloc_scheme
       write(unit,48) mem_scheme
       write(unit,44) exx_screen, scr_scheme, screen, screen*BohrToAng
       
       write(unit,46) exx_overlap
       write(unit,47) pao_scheme
       write(unit,20)
       write(unit,21) r_max, r_max*BohrToAng
       write(unit,22) r_int, r_int*BohrToAng
       write(unit,23) grid_spacing, grid_spacing*BohrToAng
       write(unit,24) ngrid
       write(unit,25) 
       write(unit,26) edge, edge*BohrToAng
       write(unit,27) volume, volume*BohrToAng**3
       write(unit,28) ngrid**3

    end if

    if (exx_psolver == 'fftw') then       
       scheme = trim(scheme)//'/'//trim(p_scheme)
       poisson_fftw: select case(p_scheme)          
          
       case('default')
          scheme = trim(scheme)//trim(p_scheme_default)
          if ( inode == ionode ) then
             !write(unit,30) solver, scheme
             !write(unit,31) ('G=0 component neglected... warning: inaccurate !')
          end if
       case('ewald')
          if (ewald_alpha < very_small) then
             ewald_alpha = real(3.0,double)
          end if
          if ( inode == ionode ) then
             !write(unit,30) solver, scheme
             !write(unit,32) ewald_alpha
          end if
          
       case('pulay')
          if (pulay_radius < very_small) then
             !pulay_radius = r_maxs*sqrt(three)
             pulay_radius = r_int*sqrt(three)
             !pulay_radius = exx_radius

          end if
          if (pulay_factor < very_small) then
             pulay_factor = one
          end if 
          pulay_radius = pulay_factor*pulay_radius
          if ( inode == ionode ) then
             !write(unit,30) solver, scheme
             !write(unit,33) pulay_factor
             !write(unit,34) pulay_radius
          end if
       case('yukawa')
          if (p_omega < very_small) then
             p_omega = magic_number
             p_omega = -log(threshold*r_int)/r_int
          end if
          !write(unit,30) solver, scheme
          !write(unit,35) p_omega
          
       case('gauss')
          if ( inode == ionode ) then
             !write(unit,30) solver, scheme
          end if
       case default
          scheme = trim(p_scheme_default)
          if ( inode == ionode ) then
             !write(unit,30) solver, scheme
             !write(unit,31) ('WARNING: G=0 component neglected !')
          end if
          
       end select poisson_fftw
       
    else if (exx_psolver == 'isf') then
       if ( inode == ionode ) then
          !write(unit,30) solver, scheme
          !write(unit,36) isf_order
       end if
    end if
   
!****lat<$
    if (inode == ionode .and. iprint_exx > 2) then
       write (io_lun,50) &
            grid_spacing,     &
            r_int,            &
            extent,           &
            trim(exx_psolver),&
            exx_alpha 
       
    end if

50 format (' EXX: gs = ',f6.4,', rc = ',f7.4,', extent = ',i4,', psolv = ',a,', alpha = ', f4.2)
!****lat>$
    !
    !
    call stop_timer(tmr_std_exx_setup,.true.)
    !
!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_X_matrix',echo=.true.)
!****lat>$
    !
    
    return

    !! EXX Formats start here
1   format(/1x,104a/)
2   format(1x,34a,8x,/)
    
    !! Grid settings
20  format(/25x,'Grid Settings: Cubic')  
21  format( 29x,'r_max_pao: ',    f12.4,' a0,',f12.4,' Ang')  
22  format( 25x,'r_integration: ',f12.4,' a0,',f12.4,' Ang')  
23  format( 26x,'grid_spacing: ', f12.8,' a0,',f12.8,' Ang')  
24  format( 25x,'n_grid_points: ',i4)  
    
25  format(/25x,'Local FFT Box: Cubic')  
26  format( 27x,'edge_length: ',f12.4,' a0,',f12.4,' Ang')
27  format( 32x,'volume: ',     f12.4,' a0,',f12.4,' Ang')
28  format( 23x,'n^3_grid_points: ',i12)  
    
    !! Poisson solver
30  format(/24x,'Poisson Solver: ',a20,/32x,'Scheme: ',a20)  
    ! Warning
31  format( 24x,a34) 
    ! Ewald scheme
32  format( 33x,'alpha: ',f8.6,' a0^{-1/2}')
    ! Pulay scheme
33  format( 27x,'cutoff_fact: ',f4.1, ' a0')
34  format( 28x,'cutoff_pot: ',f4.1, ' a0')
    ! Yukawa scheme
35  format( 33x,'omega: ',f8.6, ' a0^{-1}')
    ! ISF Psolver
36  format( 31x,'n_order: ',i3)

    !! User settings
40  format(/25x,'Exact exchange settings:')  
41  format( 32x,'method: ',a30)
42  format( 34x,'phil: ',a30)
43  format( 28x,'allocation: ',a30)
44  format( 29x,'screening: ',l,', method: ',a4, ', cutoff: ',f12.4,' a0,',f12.4,' Ang')

46  format( 27x,'overlap_box: ',l)
47  format( 27x,'pao_on_grid: ',a30)
48  format( 24x,'memory_storage: ',a30)
    
  end subroutine get_X_params
  !
  !
  subroutine exx_ewald_charge(rho,extent,dv,charge)

    use numbers,   ONLY: zero

    implicit none

    integer                   :: extent
    real(double), dimension(2*extent+1,2*extent+1,2*extent+1), &
         intent(in)           :: rho
    real(double), intent(in)  :: dv
    real(double), intent(out) :: charge
    integer                   :: i, j, k
    
    charge = zero
    
    do i = 1, 2*extent+1
       do j = 1, 2*extent+1
          do k = 1, 2*extent+1                                    
             charge = charge + rho(i,j,k)*dv
          end do
       end do
    end do
    
    return
  end subroutine exx_ewald_charge

  subroutine exx_ewald_rho(gauss,extent,alpha,r_int)
    
    use numbers,         ONLY: zero, one, two, three_halves
    use numbers,         ONLY: twopi, pi, three_halves

    implicit none
    
    integer,      intent(in)    :: extent
    real(double), intent(inout) :: alpha
    real(double), intent(in)    :: r_int
    real(double), dimension(2*extent+1,2*extent+1,2*extent+1), &
         intent(inout)          :: gauss

    real(double)             :: gnorm, dv, factor
    real(double)             :: grid_spacing
    real(double)             :: r(3)

    integer                  :: i, j, k
    
    grid_spacing = r_int/real(extent,double)     
    factor       = (alpha/pi)**three_halves
    !factor   = one
    dv           = grid_spacing**3

    gauss    = zero
    gnorm    = zero
    do i = -extent, extent      
       r(1) = real(i,double)*grid_spacing
       do j = -extent, extent
          r(2) = real(j,double)*grid_spacing
          do k = -extent, extent
             r(3) = real(k,double)*grid_spacing
             
             gauss(extent+i+1,extent+j+1,extent+k+1) = &
                  exp(-dot_product(r,r)*alpha)             
             gnorm = &
                  gnorm + gauss(extent+i+1,extent+j+1,extent+k+1)*dv

          end do
       end do
    end do

    gauss = gauss*factor
    gnorm = gnorm*factor

    return
  end subroutine exx_ewald_rho
  
  subroutine exx_ewald_pot(potential,extent,alpha,r_int)
    
    use numbers,         ONLY: zero, one, two, twopi, pi
    use functions,      ONLY: erfc

    implicit none

    integer,      intent(in)    :: extent
    real(double), intent(inout) :: alpha
    real(double), intent(in)    :: r_int
    real(double), dimension(2*extent+1,2*extent+1,2*extent+1), &
         intent(inout)          :: potential

    real(double)             :: arg, factor, dv
    real(double)             :: grid_spacing, r(3)

    integer                  :: i, j, k
    
    grid_spacing = r_int/real(extent,double)     
    dv           = grid_spacing**3   

    potential    = zero
    do i = -extent, extent      
       r(1) = real(i,double)*grid_spacing
       do j = -extent, extent      
          r(2) = real(j,double)*grid_spacing
          do k = -extent, extent      
             r(3) = real(k,double)*grid_spacing
             
             if (i == 0 .and. j == 0 .and. k==0) then
                potential(extent+i+1,extent+j+1,extent+k+1) = two*sqrt(alpha/pi)
             else
                arg    = sqrt(dot_product(r,r))*sqrt(alpha)
                factor = one/sqrt(dot_product(r,r))
                potential(extent+i+1,extent+j+1,extent+k+1) = (one - erfc(arg))*factor             
             end if
             
          end do
       end do
    end do

    return
  end subroutine exx_ewald_pot

!!$ For ISF real-space Poisson solver
!!$
!!$subroutine createBeylkin(p_gauss,w_gauss,r_int)
!!$
!!$    use numbers, ONLY: zero, one, two, twopi, pi
!!$    
!!$    implicit none
!!$
!!$    ! << Passed variables >>
!!$    real(double), intent(in) :: r_int
!!$
!!$    ! << Local variables >>
!!$    integer, parameter :: n_gauss = 89
!!$    real(double)       :: ur_gauss,dr_gauss,acc_gauss
!!$    real(double)       :: norm, vol, factor1, factor2 
!!$    real(double)       :: p_gauss(n_gauss)
!!$    real(double)       :: w_gauss(n_gauss)
!!$
!!$    p_gauss = zero
!!$    w_gauss = zero
!!$
!!$    vol     = (real(2,double)*r_int)**3
!!$    norm    = sqrt(vol)
!!$
!!$    factor1 = one/norm
!!$    factor2 = one/vol
!!$        
!!$    call gequad(n_gauss,p_gauss,w_gauss,ur_gauss,dr_gauss,acc_gauss)
!!$        
!!$    return
!!$  end subroutine createBeylkin


  !
!!$  subroutine get_neighdat(nb,ia,hl,ind,part,verbose,unit)
!!$    
!!$    use numbers,        only: three
!!$    use group_module,   only: parts
!!$    use cover_module,   only: BCS_parts
!!$    use global_module,  only: id_glob,atom_coord
!!$    use species_module, only: species_label
!!$    use atomic_density, only: atomic_density_table
!!$    use numbers,        only: zero, one, two, twopi, pi
!!$    !
!!$    use matrix_data, only: mat, Hrange, Srange
!!$    use exx_types,   only: prim_atomic_data, neigh_atomic_data
!!$    use exx_types,   only: tmr_std_exx_fetch
!!$    !
!!$    !
!!$    implicit none
!!$    !
!!$    type(neigh_atomic_data), intent(inout) :: nb
!!$    type(prim_atomic_data),  intent(inout) :: ia
!!$    type(neigh_atomic_data), intent(inout) :: hl
!!$    integer, intent(in)                    :: ind
!!$    integer, intent(in)                    :: part
!!$    logical, intent(in)                    :: verbose
!!$    !
!!$    real(double)      :: xyz_Ang(3)      
!!$    character(len=20) :: filename
!!$    integer, optional :: unit
!!$    !
!!$    call start_timer(tmr_std_exx_fetch)
!!$    ! Get u(v) data
!!$    nb%nb  = ind
!!$    nb%ist = mat(part,Hrange)%i_acc(ia%pr) + nb%nb - 1
!!$    nb%npc = mat(part,Hrange)%i_part(nb%ist)
!!$    nb%nic = mat(part,Hrange)%i_seq(nb%ist)
!!$    !
!!$    nb%global_part = BCS_parts%lab_cell(mat(part,Hrange)%i_part(nb%ist))
!!$    nb%global_num  = id_glob(parts%icell_beg(nb%global_part) +  &
!!$         mat(part,Hrange)%i_seq(nb%ist)-1)
!!$    !
!!$    nb%spec = BCS_parts%spec_cover(BCS_parts%icover_ibeg(nb%npc) + nb%nic-1)
!!$    nb%name = species_label(nb%spec)               
!!$    nb%radi = atomic_density_table(nb%spec)%cutoff
!!$    nb%nsup = mat(part,Hrange)%ndimj(nb%ist)                   
!!$    !
!!$    ! Calculate R_u
!!$    nb%xyz_nb(1) = atom_coord(1,nb%global_num)
!!$    nb%xyz_nb(2) = atom_coord(2,nb%global_num)
!!$    nb%xyz_nb(3) = atom_coord(3,nb%global_num)
!!$    !
!!$    nb%xyz_cv(1) = BCS_parts%xcover(BCS_parts%icover_ibeg(nb%npc) + &
!!$         nb%nic - 1)
!!$    nb%xyz_cv(2) = BCS_parts%ycover(BCS_parts%icover_ibeg(nb%npc) + &
!!$         nb%nic - 1)
!!$    nb%xyz_cv(3) = BCS_parts%zcover(BCS_parts%icover_ibeg(nb%npc) + &
!!$         nb%nic - 1) 
!!$    !
!!$    ! Calculate R_iu
!!$    nb%xyz(1) = - nb%xyz_nb(1) + hl%xyz_hl(1)
!!$    nb%xyz(2) = - nb%xyz_nb(2) + hl%xyz_hl(2)
!!$    nb%xyz(3) = - nb%xyz_nb(3) + hl%xyz_hl(3)
!!$    nb%r      = sqrt(dot_product(nb%xyz,nb%xyz))
!!$    !
!!$    ! Calculate D_iu
!!$    nb%d      = sqrt(three)*(nb%radi + hl%radi)
!!$    !
!!$    xyz_Ang(1) = nb%xyz_nb(1)
!!$    xyz_Ang(2) = nb%xyz_nb(2)
!!$    xyz_Ang(3) = nb%xyz_nb(3)
!!$    !
!!$    if ( exx_debug ) then
!!$       write(unit,'(I8,4X,A9,4X,A,I8,A2,I3,A,6X,A2,X,I8,X,4F12.4,3X,I3,2I8,X,F7.3,F14.8)')        &
!!$            nb%global_part,'{j\beta} ','{',nb%global_num,'\ ',nb%nsup,'}', nb%name, nb%global_num, &
!!$            xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), zero, nb%spec, ind, nb%nsup, nb%radi,zero
!!$    end if
!!$    call stop_timer(tmr_std_exx_fetch,.true.)
!!$
!!$    return
!!$  end subroutine get_neighdat

  subroutine get_halodat(hl,kl,ind,part_cover,part,which,verbose,unit)

    use numbers,        only: three
    use group_module,   only: parts
    use cover_module,   only: BCS_parts
    use global_module,  only: id_glob,atom_coord, species_glob
    use species_module, only: species_label, nsf_species
    use atomic_density, only: atomic_density_table
    use numbers,        only: zero, one, two, twopi, pi
    !
    use matrix_data, only: mat, Hrange, Srange
    use exx_types, only: neigh_atomic_data
    use exx_types, only: tmr_std_exx_fetch
    !
    implicit none
    !
    type(neigh_atomic_data), intent(inout) :: hl
    type(neigh_atomic_data), intent(inout) :: kl
    integer,          intent(in)           :: ind
    integer,          intent(in)           :: part_cover
    integer,          intent(in)           :: part
    logical,          intent(in)           :: verbose
    character(len=1), intent(in)           :: which
    !
    real(double)      :: xyz_Ang(3)      
    character(len=20) :: filename
    integer, optional :: unit
    !
    call start_timer(tmr_std_exx_fetch)
    ! Get u(v) data
    hl%npc = part_cover
    hl%nic = ind
    !
    hl%global_part = parts%icell_beg(part)
    hl%global_num  = id_glob(parts%icell_beg(part)+ind-1)
    !
    !if (which == 'k') then
       !hl%spec  = BCS_parts%spec_cover(BCS_parts%icover_ibeg(hl%npc)+hl%nic-1)
    hl%spec  = species_glob(hl%global_num)
    !else
    !   hl%spec = BCS_parts%spec_cover(part_cover+ind-1)
    !end if
    !
    !hl%spec     = BCS_parts%spec_cover(BCS_parts%icover_ibeg(hl%npc)+hl%nic-1)
    !
    hl%name = species_label(hl%spec)               
    hl%radi = atomic_density_table(hl%spec)%cutoff
    hl%nsup = nsf_species(hl%spec)
    !
    ! Calculate R_u
    hl%xyz_hl(1) = atom_coord(1,hl%global_num)
    hl%xyz_hl(2) = atom_coord(2,hl%global_num)
    hl%xyz_hl(3) = atom_coord(3,hl%global_num)
    !
    hl%xyz_cv(1) = BCS_parts%xcover(part_cover+ind-1)
    hl%xyz_cv(2) = BCS_parts%ycover(part_cover+ind-1)
    hl%xyz_cv(3) = BCS_parts%zcover(part_cover+ind-1)
    !
    ! Calculate R_iu
    hl%xyz(1) = - hl%xyz_hl(1) + kl%xyz_hl(1)
    hl%xyz(2) = - hl%xyz_hl(2) + kl%xyz_hl(2)
    hl%xyz(3) = - hl%xyz_hl(3) + kl%xyz_hl(3)
    hl%r      = sqrt(dot_product(hl%xyz,hl%xyz))
    !
    ! Calculate D_iu
    hl%d      = sqrt(three)*(hl%radi + kl%radi)
    !
    xyz_Ang(1) = hl%xyz_hl(1)
    xyz_Ang(2) = hl%xyz_hl(2)
    xyz_Ang(3) = hl%xyz_hl(3)
    !
    if ( exx_debug ) then
       if (which == 'k') then
          write(unit,'(I8,6X,A9,2X,A,I8,A2,I3,A,6X,A2,X,I8,X,4F12.4,3X,I3,2I8,X,F7.3)')     &
               part,'{k\gamma}','{',hl%global_num,'\ ',hl%nsup,'}', hl%name, hl%global_num, &
               xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), zero, hl%spec, ind, hl%nsup, hl%radi
       else
          write(unit,'(I8,6X,A9,2X,A,I8,A2,I3,A,6X,A2,X,I8,X,4F12.4,3X,I3,2I8,X,F7.3)')     &
               part,'{l\delta}','{',hl%global_num,'\ ',hl%nsup,'}', hl%name, hl%global_num, &
               xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), zero, hl%spec, ind, hl%nsup, hl%radi
       end if
    end if    
    call stop_timer(tmr_std_exx_fetch,.true.)

    return
  end subroutine get_halodat

!!   2016/07/29 18:30 nakata
!!    Renamed supports_on_atom -> blips_on_atom
  subroutine get_iprimdat(ia,hl,ind,iprim,part,verbose,unit)
    
    use units,          only: BohrToAng
    use numbers,        only: very_small, zero, three
    use primary_module, only: bundle
    use species_module, only: species_label,nsf_species
    use atomic_density, only: atomic_density_table
    use support_spec_format, only: blips_on_atom, flag_one_to_one
    !
    use matrix_data, only: mat, Hrange, Srange
    use exx_types, only: prim_atomic_data, neigh_atomic_data
    use exx_types, only: tmr_std_exx_fetch
    use group_module,only: parts 
    use GenComms,    only: inode
    !
    implicit none
    !
    type(prim_atomic_data), intent(inout) :: ia
    type(neigh_atomic_data),intent(inout) :: hl
    integer, intent(in)                   :: ind
    integer, intent(in)                   :: iprim
    integer, intent(in)                   :: part
    logical, intent(in)                   :: verbose
    !
    real(double)      :: xyz_Ang(3)      
    character(len=20) :: filename
    integer, optional :: unit
    !
    call start_timer(tmr_std_exx_fetch)
    ! Get u(v) data
    ia%pr  = ind
    ia%num = bundle%nm_nodbeg(part) + ia%pr - 1 
    ia%ip  = bundle%ig_prim(iprim)             

    ia%spec = bundle%species(ia%num)
    ia%name = species_label(ia%spec)             
    ia%radi = atomic_density_table(ia%spec)%cutoff
    ia%labcell = ia%ip 
    !if (flag_one_to_one) then
    !   ia%nsup = blips_on_atom(ia%spec)%nsuppfuncs
    !else
    !   ia%nsup = blips_on_atom(ia%ip)%nsuppfuncs
    !end if
    ia%nsup = nsf_species(ia%spec)

    ! Calculate R_u    
    ia%xyz_ip(1) = bundle%xprim(ia%num)
    ia%xyz_ip(2) = bundle%yprim(ia%num)
    ia%xyz_ip(3) = bundle%zprim(ia%num)

    ! Calculate R_iu
    ia%xyz(1) = - ia%xyz_ip(1) + hl%xyz_hl(1)
    ia%xyz(2) = - ia%xyz_ip(2) + hl%xyz_hl(2)
    ia%xyz(3) = - ia%xyz_ip(3) + hl%xyz_hl(3)
    ia%r      = sqrt(dot_product(ia%xyz,ia%xyz))
    !
    ! Calculate D_iu
    ia%d      = sqrt(three)*(ia%radi + hl%radi)
    !    
    xyz_Ang(1) = ia%xyz_ip(1)
    xyz_Ang(2) = ia%xyz_ip(2)
    xyz_Ang(3) = ia%xyz_ip(3)
    !
    if ( exx_debug ) then
       write(unit,'(I8,2X,A9,6X,A,I8,A2,I3,A,6X,A2,X,I8,X,4F12.4,3X,I3,2I8,X,F7.3,F12.8)')           &
            parts%ngnode(parts%inode_beg(inode)+part-1),'{i\alpha}','{',ia%labcell,'\ ',ia%nsup,'}', &
            ia%name, ia%labcell, xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), zero, ia%spec, ind, ia%nsup, zero
    end if
    call stop_timer(tmr_std_exx_fetch,.true.)

    return
  end subroutine get_iprimdat

!!$  subroutine get_ghostdat(hl)
!!$    
!!$    use units,          only: BohrToAng
!!$    use numbers,        only: very_small, zero, three
!!$    use primary_module, only: bundle
!!$    use species_module, only: species_label,nsf_species
!!$    use atomic_density, only: atomic_density_table
!!$    use support_spec_format, only: blips_on_atom, flag_one_to_one
!!$    !
!!$    use matrix_data, only: mat, Hrange, Srange
!!$    use exx_types, only: prim_atomic_data, neigh_atomic_data
!!$    use exx_types, only: tmr_std_exx_fetch
!!$    use group_module,only: parts 
!!$    use GenComms,    only: inode
!!$    !
!!$    implicit none
!!$    !
!!$    type(neigh_atomic_data),intent(inout) :: hl
!!$    !
!!$    !call start_timer(tmr_std_exx_fetch)
!!$    ! Get u(v) data
!!$    hl%xyz_hl(1) = 0_double
!!$    hl%xyz_hl(2) = 0_double
!!$    hl%xyz_hl(3) = 0_double
!!$    !call stop_timer(tmr_std_exx_fetch,.true.)
!!$    !
!!$    return
!!$  end subroutine get_ghostdat

  !!***
end module exx_module
