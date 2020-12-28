module exx_poisson

  use GenComms,            only: cq_abort
  use datatypes,           only: double, double_cplx
  !use exx_types,           only: fftwrho3d, exx_debug
  use fft_interface_module,only: fft3_exec_wrapper, fft3_init_wrapper
  
contains
  !
  !
  subroutine exx_scal_rho_3d(inode,extent,r_int,scheme,cutoff,omega,n_gauss, &
       p_gauss,w_gauss,reckernel)

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

    complex(double), dimension(2*extent+1,2*extent+1,2*extent+1), intent(out) :: reckernel

    
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

    reckernel = zero

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
                   reckernel(x,y,z) = cmplx(0,0,double_cplx)          
                case('ewald')
                   reckernel(x,y,z) = cmplx(0,0,double_cplx)          
                case('pulay')
                   reckernel(x,y,z) = cmplx(cutoff**2/two,0,double_cplx)     
                case('yukawa')
                   reckernel(x,y,z) = cmplx(0,0,double_cplx) 
                case('gauss')
                   Gkernel_zero = zero
                   do gauss = 1, n_gauss
                      Gkernel_zero = Gkernel_zero + &
                           w_gauss(gauss)*sqrt_pi/(four*p_gauss(gauss)**three_halves)
                   end do
                   reckernel(x,y,z) = cmplx(Gkernel_zero,double_cplx)          
                case default 
                   reckernel(x,y,z) = cmplx(zero,zero,double_cplx)          
                end select poisson_zero
             else
                poisson_xyz: select case(scheme)
                case('default')
                   reckernel(x,y,z) = cmplx(one/r_G**2,zero,double_cplx)
                case('ewald')
                   reckernel(x,y,z) = cmplx(one/r_G**2,zero,double_cplx)
                case('pulay')
                   reckernel(x,y,z) = (one - cos(r_G*cutoff))/r_G**2
                case('yukawa')
                   reckernel(x,y,z) = cmplx(one/r_G**2,zero,double_cplx)
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
                   reckernel(x,y,z) = cmplx(Gkernel,zero,double_cplx)
                case default
                   reckernel(x,y,z) = cmplx(one/r_G**2,zero,double_cplx)
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
       alpha,omega,n_gauss,p_gauss,w_gauss,fftwrho,reckernel)
    
    use numbers,   ONLY: zero, one, fourpi
    use exx_types, ONLY: fftw3d                       ! FFTW
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
         intent(in) :: rho
    complex(double), dimension(2*extent+1,2*extent+1,2*extent+1),&
         intent(in) :: reckernel
    
    character(*), intent(in) :: scheme
    real(double), intent(in) :: alpha
    real(double), intent(in) :: omega        
    integer,      intent(in) :: n_gauss
    real(double), dimension(n_gauss) :: p_gauss, w_gauss

    type(fftw3d), intent(in) :: fftwrho
    
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
       fftwrho%arrayin  = cmplx(rho,zero,double_cplx)

       ! FFT_F[rho(r)] => rho(G)      
       call fft3_exec_wrapper( fftwrho%arrayin, ng , +1 )

       ! scale[rho_(G)] = 4pi*rho(G)/|G|^2
       fftwrho%arrayin = fourpi*fftwrho%arrayin*reckernel
       
       ! FFT_B[4pi*rho(G)/|G|^2] = V(r')
       call fft3_exec_wrapper( fftwrho%arrayin, ng , -1 )
       
       ! Normalization
       potential = real(fftwrho%arrayin) / fftwnorm**2

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
    use functions,      ONLY: erfc_cq

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
                potential(extent+i+1,extent+j+1,extent+k+1) = (one - erfc_cq(arg))*factor             
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
  
end module exx_poisson
