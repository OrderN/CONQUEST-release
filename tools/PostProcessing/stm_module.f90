  ! Module for STM simulation following Paz and Soler pssb 243, 1080
module stm

  use datatypes

  implicit none

contains

  subroutine make_stm_current

    use datatypes
    use numbers
    use units
    use local, ONLY: nptsx, nptsy, nptsz, n_bands_active, nkp, nprocs, stm_bias, &
         eigenvalues, current, grid_x, grid_y, grid_z, wtk, band_active_kp, n_bands_total, efermi
    use dimens, ONLY: GridCutoff
    use fft_interface_module, ONLY: fft3_init_wrapper_nonu, fft3_exec_wrapper_nonu, fft3_dest_wrapper_nonu
    use block_module, only: n_pts_in_block
    use process, ONLY: read_domain, pao_dpao_to_grid
    use io_module, ONLY: get_file_name
    use output, ONLY: write_cube
    use read, ONLY: read_psi_coeffs, read_eigenvalues
    use angular_coeff_routines, ONLY: set_prefac_real, set_fact, set_prefac

    implicit none

    integer :: band, nk, proc, i, idum1, idum2
    real(double) :: e_tip, band_rel
    real(double), dimension(3) :: dg
    real(double) :: rbx, rby, rbz, sq
    real(double), dimension(:), allocatable :: k_vec_x, k_vec_y, k_vec_z
    real(double), dimension(:,:,:), allocatable :: delta_of_S, g_k, rho_band
    real(double), dimension(:,:,:,:), allocatable :: cvec, drho_band
    complex(double_cplx), dimension(:,:,:), allocatable :: A_of_r, B_of_r, psi, STM_current
    complex(double_cplx), dimension(:,:,:,:), allocatable :: dpsi
    character(len=50) :: ci, filename
    logical :: use_band

    call set_fact(8)
    call set_prefac(9)
    call set_prefac_real(9)
    !dg = pi/sqrt(two*GridCutoff) ! Equivalent grid spacing
    dg(1) = grid_x!/BohrToAng
    dg(2) = grid_y!/BohrToAng
    dg(3) = grid_z!/BohrToAng
    ! Read eigenvalues
    !if(stm_bias<zero) then
       call read_eigenvalues!(stm_bias, zero)
    !else
    !   call read_eigenvalues(zero, stm_bias)
    !end if
    ! Read coefficients
    call read_psi_coeffs("Process")
    ! Allocate space for variables
    allocate(delta_of_S(nptsx,nptsy,nptsz), g_k(nptsx,nptsy,nptsz))
    delta_of_S = zero
    g_k = zero
    allocate(cvec(nptsx,nptsy,nptsz,3))
    cvec = zero
    allocate(current(nptsx,nptsy,nptsz), STM_current(nptsx,nptsy,nptsz))
    current = zero
    STM_current = zero
    !call create_delta_of_S_and_c(delta_of_S, cvec)
    allocate(A_of_r(nptsx,nptsy,nptsz))
    allocate(B_of_r(nptsx,nptsy,nptsz))
    ! Space for psi
    allocate(psi(nptsx,nptsy,nptsz))
    allocate(dpsi(nptsx,nptsy,nptsz,3))
    allocate(rho_band(nptsx,nptsy,nptsz))
    allocate(drho_band(nptsx,nptsy,nptsz,3))
    ! Read charge
    do proc = 1,nprocs
       write(ci,'("chden")')
       call get_file_name(ci,nprocs,proc,filename)
       open(unit=17,file=filename)
       call read_domain(17,proc,rho_band)
       close(unit=17)
    end do
    write(ci,'("TotalCharge")')
    call write_cube(rho_band,ci)
    do i=1,3
       drho_band(:,:,:,i) = (cshift(rho_band,1,i) - cshift(rho_band,-1,i))/dg(i)
       write(ci,'("drho",I0.1)') i
       call write_cube(drho_band(:,:,:,i),ci)
    end do
    ! Create delta(S) and c(r) from band charge density
    call create_delta_of_S_and_c(rho_band, drho_band, delta_of_S, cvec)
    ! Create k vectors
    allocate(k_vec_x(nptsx),k_vec_y(nptsy),k_vec_z(nptsz))
    call create_k_vectors(k_vec_x,k_vec_y,k_vec_z)
    ! Set up FFTs
    call fft3_init_wrapper_nonu(nptsx, nptsy, nptsz)
    ! Loop over eigenstates
    do band = 1, n_bands_total
       !write(*,*) 'Band: ',band_no(band)
       do nk = 1, nkp
          !write(*,*) 'E(',band_no(band),nk,"): ", eigenvalues(nk,band,1)
          band_rel = eigenvalues(band, nk, 1) - efermi(1)
          use_band = .false.
          if(band_rel*stm_bias > zero .and. abs(band_rel) < abs(stm_bias)) use_band = .true.
          if(band_active_kp(band,nk,1)==1 .and. use_band) then
             current = zero
             psi = zero
             dpsi = zero
             call pao_dpao_to_grid(band,nk,1,psi,dpsi)
             !write(*,*) "Made band"
             ! Tip energy
             e_tip = eigenvalues(band,nk,1) - stm_bias
             ! Create Green's function in reciprocal space
             call create_g_k(g_k, e_tip, k_vec_x, k_vec_y, k_vec_z)
             !write(*,*) 'Made g_k'
             ! Create A(r) = c(r) dot grad psi(r)
             A_of_r = cmplx(zero,zero)
             do i = 1, 3
                !A_of_r = A_of_r + cvec(:,:,:,i) * (cshift(current,1,i) - cshift(current,-1,i))/dg(i)
                A_of_r = A_of_r + cvec(:,:,:,i) * dpsi(:,:,:,i)
             end do
             !write(*,*) 'Made A(r)'
             ! FFT
             call fft3_exec_wrapper_nonu(A_of_r, nptsx, nptsy, nptsz, -1) ! Forward
             ! Create B(r) by component: FFT, dot with k and accumulate
             ! x
             B_of_r(:,:,:) = cvec(:,:,:,1) * psi
             call fft3_exec_wrapper_nonu(B_of_r, nptsx, nptsy, nptsz, -1) ! Forward
             do i = 1, nptsx
                A_of_r(i,:,:) = A_of_r(i,:,:) + cmplx(zero,k_vec_x(i))*B_of_r(i,:,:)
             end do
             ! y
             B_of_r(:,:,:) = cvec(:,:,:,2) * psi
             call fft3_exec_wrapper_nonu(B_of_r, nptsx, nptsy, nptsz, -1) ! Forward
             do i = 1, nptsy
                A_of_r(:,i,:) = A_of_r(:,i,:) + cmplx(zero,k_vec_y(i))*B_of_r(:,i,:)
             end do
             ! z
             B_of_r(:,:,:) = cvec(:,:,:,3) * psi
             call fft3_exec_wrapper_nonu(B_of_r, nptsx, nptsy, nptsz, -1) ! Forward
             do i = 1, nptsz
                A_of_r(:,:,i) = A_of_r(:,:,i) + cmplx(zero,k_vec_z(i))*B_of_r(:,:,i)
             end do
             ! Scale by g_k
             B_of_r = A_of_r * g_k
             ! Transform back and accumulate into current
             call fft3_exec_wrapper_nonu(B_of_r, nptsx, nptsy, nptsz, +1) ! Backward
             current = current + wtk(nk)*B_of_r*conjg(B_of_r)
          end if
       end do ! nkp
    end do ! n_bands_active
    !current = STM_current*conjg(STM_current)
    ci = "TestImage"
    call write_cube(current,ci)
    call fft3_dest_wrapper_nonu
    return
  end subroutine make_stm_current
  
  subroutine create_delta_of_S_and_c(rho, drho, delta_of_S, cvec)

    use datatypes
    use numbers
    use units
    use dimens, ONLY: GridCutoff
    use local, ONLY: current, nptsx, nptsy, nptsz, charge_stub, nprocs, grid_x, grid_y, grid_z
    use block_module, only: n_pts_in_block
    use io_module, ONLY: get_file_name
    use output, ONLY: write_cube

    implicit none

    ! Passed variables
    real(double), dimension(nptsx, nptsy, nptsz) :: delta_of_S, rho
    real(double), dimension(nptsx, nptsy, nptsz, 3) :: cvec, drho

    ! Local variables
    real(double) :: rho0, dS, rho_min, rho_max, s_over_ds_2
    real(double), dimension(3) :: dg
    integer :: proc, i
    character(len=50) :: ci

    delta_of_S = zero
    ! These could be made user-defined
    rho0 = three/(fourpi*729.0_double) ! r_s = 9 bohr
    dS = eight*two*pi*sqrt(five/(HaToeV*GridCutoff)) ! NB 5 eV is an assumed workfunction
    write(*,*) 'deltaS is ',dS
    !dg = pi/sqrt(two*GridCutoff) ! Equivalent grid spacing
    !write(*,*) 'Effective grid spacing: ',dg
    dg(1) = grid_x!/BohrToAng
    dg(2) = grid_y!/BohrToAng
    dg(3) = grid_z!/BohrToAng
    write(*,*) 'Min and max values of rho: ',minval(rho), maxval(rho)
    write(*,*) 'Present rho0: ',rho0
    rho0 = min(rho0,0.005*maxval(rho))
    write(*,*) 'New rho0: ',rho0
    ! Set limits on values of rho
    rho_min = rho0*exp(-dS)
    rho_max = rho0*exp(dS)
    !where (rho>rho_min) ! .and. rho<rho_max)
    delta_of_S = log(rho/rho0)
    !elsewhere
    !   delta_of_S = zero
    !end where
    !delta_of_S(:,:,1:nptsz/2) = zero
    ci = "justS"
    call write_cube(delta_of_S,ci)
    ! Assemble delta(S)
    where (rho>rho_min .and. rho<rho_max) !>rho_min .and. rho < rho_max)
       delta_of_S = (one - (delta_of_S/dS)*(delta_of_S/dS))!log(rho/rho0)/dS)*(log(rho/rho0)/dS))
       delta_of_S = fifteen*delta_of_S * delta_of_S/(sixteen*dS)
    elsewhere
       delta_of_S = zero
    end where
    delta_of_S(:,:,1:nptsz/2) = zero
    ! Write out
    ci = "deltaS"
    call write_cube(delta_of_S,ci)
    ! Now make c = delta(S).grad rho/rho 
    do i = 1,3
       where (delta_of_S>RD_ERR) ! .and. abs(rho) > RD_ERR)
          cvec(:,:,:,i) = drho(:,:,:,i)*delta_of_S/rho
       elsewhere
          cvec(:,:,:,i) = zero
       end where
       if(i==1) ci = "cvec_x"
       if(i==2) ci = "cvec_y"
       if(i==3) ci = "cvec_z"
       call write_cube(cvec(:,:,:,i),ci)
    end do
    !do i = 1,3
       ! Grad rho using finite differences (centred)
       !cvec(:,:,:,i) = (cshift(rho,1,i) - cshift(rho,-1,i))/dg(i)
       ! delta(S)*grad rho/rho
       !where (delta_of_S>RD_ERR) ! .and. abs(rho) > RD_ERR)
       !   cvec(:,:,:,i) = cvec(:,:,:,i)*delta_of_S/rho
       !elsewhere
       !   cvec(:,:,:,i) = zero
       !end where
       !current = cvec(:,:,:,i)
       !if(i==1) ci = "cvec_x"
       !if(i==2) ci = "cvec_y"
       !if(i==3) ci = "cvec_z"
       !call write_cube(ci)
    !end do
  end subroutine create_delta_of_S_and_c

  subroutine create_k_vectors(kx,ky,kz)

    use datatypes
    use numbers
    use local, ONLY: nptsx, nptsy, nptsz, stm_z_min, stm_z_max, stm_x_min, stm_x_max, &
         stm_y_min, stm_y_max

    implicit none

    ! Passed variables
    real(double), dimension(nptsx) :: kx
    real(double), dimension(nptsy) :: ky
    real(double), dimension(nptsz) :: kz

    ! Local variables
    integer :: i
    real(double) :: dk_x, dk_y, dk_z

    dk_x = pi/(stm_x_max - stm_x_min)
    dk_y = pi/(stm_y_max - stm_y_min)
    dk_z = pi/(stm_z_max - stm_z_min)
    do i=1,nptsx/2+1
       kx(i) = dk_x*real(i-1,double)
    end do
    do i=nptsx/2+2,nptsx
       kx(i) = dk_x*real(i-1-nptsx,double)
    end do
    do i=1,nptsy/2+1
       ky(i) = dk_y*real(i-1,double)
    end do
    do i=nptsy/2+2,nptsy
       ky(i) = dk_y*real(i-1-nptsy,double)
    end do
    do i=1, nptsz/2+1
       kz(i) = dk_z*real(i-1,double)
    end do
    do i=nptsz/2+2,nptsz
       kz(i) = dk_z*real(i-1-nptsz,double)
    end do
  end subroutine create_k_vectors
  
  subroutine create_g_k(g_k, e_tip, kx, ky, kz)

    use datatypes
    use numbers
    use local, ONLY: nptsx, nptsy, nptsz
    use units

    implicit none

    ! Passed variables
    real(double), dimension(nptsx, nptsy, nptsz) :: g_k
    real(double) :: e_tip
    real(double), dimension(nptsx) :: kx
    real(double), dimension(nptsy) :: ky
    real(double), dimension(nptsz) :: kz

    ! Local variables
    integer :: i, j, k
    real(double) :: k_x2, k_y2, k_z2, kappa_2

    !kappa_2 = two*(five/HaToeV - e_tip)
    kappa_2 = two*(five - e_tip)/HaToeV
    do k = 1, nptsz
       k_z2 = kz(k)*kz(k)
       do j=1, nptsy
          k_y2 = ky(j)*ky(j)
          do i=1, nptsx
             k_x2 = kx(i)*kx(i)
             g_k(i,j,k) = one/(kappa_2 + k_x2 + k_y2 + k_z2)
          end do
       end do
    end do
    return
  end subroutine create_g_k

  subroutine tersoff_hamann

    use datatypes
    use numbers
    use local, ONLY: nkp, wtk, efermi, current, nptsx, nptsy, nptsz, eigenvalues, stm_bias, fermi_offset, &
         n_bands_total, root_file, grid_z, band_full_to_active
    use output, ONLY: write_dx_density, write_cube, write_dx_coords
    use global_module, only : nspin
    use read, ONLY: read_eigenvalues, read_psi_coeffs
    use process, ONLY: pao_to_grid
    use units, ONLY: HaToeV
    use angular_coeff_routines, ONLY: set_prefac_real, set_fact, set_prefac

    implicit none

    integer :: proc, band, nk, kp, ispin, i
    real(double) :: weight, rbx, rby, rbz, sq, test
    real(double), dimension(2) :: Emin, Emax
    character(len=50) :: filename, ci
    complex(double_cplx), dimension(:,:,:), allocatable :: psi

    write(*,*) 'Processing STM'
    ! Read eigenvalues
    call read_eigenvalues
    ! Read eigenvector coefficients
    call read_psi_coeffs("Process")
    call set_fact(8)
    call set_prefac(9)
    call set_prefac_real(9)
    allocate(current(nptsx,nptsy,nptsz))
    allocate(psi(nptsx,nptsy,nptsz))

    fermi_offset = fermi_offset/HaToeV   ! default of fermi_offset is zero
    if(stm_bias<0) then
       Emin = stm_bias/HaToeV + Efermi + fermi_offset
       Emax = Efermi + fermi_offset
    else
       Emin = Efermi + fermi_offset
       Emax = stm_bias/HaToeV + Efermi + fermi_offset
    end if
    if(nspin==1) then
       write(*,fmt='(2x,"Including bands between ",f7.3," and ",f7.3," eV")') Emin(1)*HaToeV, Emax(1)*HaToeV
    end if

    current = zero
    do ispin=1,nspin
       do band=1,n_bands_total
          do kp = 1,nkp
             if(eigenvalues(band,kp,ispin) >= Emin(ispin) .and. &
                  eigenvalues(band,kp,ispin) <= Emax(ispin)) then
                call pao_to_grid(band_full_to_active(band), kp, ispin, psi)
                current = current + wtk(kp)*psi*conjg(psi)
             end if
          end do ! kp
       end do ! bands = 1, n_bands_total
    end do
    ! Output STM image
    ci = trim(root_file)
    ! VESTA will draw a surface on the top of the cube if the base is non-zero...
    current(:,:,1) = zero
    call write_cube(current,ci)
    ! Output height file for colouring STM image
    current = zero
    do i=1, nptsz
       current(:,:,i) = real(i-1,double)*grid_z
    end do
    ci = trim("STMHeight")
    call write_cube(current,ci)
    deallocate(current, psi)
    return
  end subroutine tersoff_hamann

  ! Old implementation of broadening to be included
  !%%!                if(stm_bias<0) then
  !%%!                   if(eigenvalues(nk,band)>efermi) then
  !%%!                      sq = (eigenvalues(nk,band)-efermi)/stm_broad
  !%%!                      weight = wtk(nk)/(exp(sq)+one)
  !%%!                   else if(eigenvalues(nk,band)<(efermi+stm_bias)) then
  !%%!                      sq = (eigenvalues(nk,band)-efermi-stm_bias)/stm_broad
  !%%!                      weight = wtk(nk)/(exp(sq)+one)
  !%%!                   else
  !%%!                      weight = wtk(nk)
  !%%!                   end if
  !%%!                else
  !%%!                   if(eigenvalues(nk,band)<efermi) then
  !%%!                      sq = (eigenvalues(nk,band)-efermi)/stm_broad
  !%%!                      weight = wtk(nk)/(exp(sq)+one)
  !%%!                   else if(eigenvalues(nk,band)>(efermi+stm_bias)) then
  !%%!                      sq = (eigenvalues(nk,band)-efermi-stm_bias)/stm_broad
  !%%!                      weight = wtk(nk)/(exp(sq)+one)
  !%%!                   else
  !%%!                      weight = wtk(nk)
  !%%!                   end if
  !%%!                end if
  
end module stm
