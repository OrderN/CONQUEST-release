module process

  use datatypes

  implicit none

  ! Maximum possible number of spin components for simplicity
  real(double), dimension(4) :: range_offset


contains

  subroutine assign_blocks

    use datatypes
    use local, ONLY: block_store, nprocs, block_size_x, block_size_y, block_size_z, &
         stm_z_min, stm_z_max, stm_x_min, stm_x_max, stm_y_min, stm_y_max


    implicit none

    integer :: proc, iblock, ig1, ind_group, block_x, block_y, block_z, nblock
    real(double) :: rbx, rby, rbz

    nblock = 0
    do proc=1,nprocs
       block_store(proc)%active = 0
       do iblock=1, block_store(proc)%num_blocks
          ! Create starting point for block using ind_block
          nblock = nblock + 1
          ! Check min area - if the RHS of block is in LHS of area, keep it
          rbx = block_size_x*real(block_store(proc)%nx(iblock),double)
          rby = block_size_y*real(block_store(proc)%ny(iblock),double)
          rbz = block_size_z*real(block_store(proc)%nz(iblock),double)
          if(rbx>=stm_x_min.AND.rby>=stm_y_min.AND.rbz>=stm_z_min) &
               block_store(proc)%active(iblock) = block_store(proc)%active(iblock) + 1
          !write(*,*) 'RHS: ',rbx,rby,rbz
          ! Now LHS of block and RHS of area
          rbx = block_size_x*real(block_store(proc)%nx(iblock)-1,double)
          rby = block_size_y*real(block_store(proc)%ny(iblock)-1,double)
          rbz = block_size_z*real(block_store(proc)%nz(iblock)-1,double)
          if(rbx<=stm_x_max.AND.rby<=stm_y_max.AND.rbz<=stm_z_max) &
               block_store(proc)%active(iblock) = block_store(proc)%active(iblock) + 1
          !write(*,*) 'RHS: ',rbx,rby,rbz
          !write(*,*) 'Active: ',iblock!block_store(proc)%active(iblock)
          if(block_store(proc)%active(iblock)==2) then
             !write(*,*) 'Active: ',iblock!block_store(proc)%active(iblock)
             block_store(proc)%active(iblock)=1
          else
             block_store(proc)%active(iblock)=0
          end if
       end do
    end do
    return
  end subroutine assign_blocks

  subroutine process_charge

    use datatypes
    use numbers
    use local, ONLY: nprocs, n_bands_active, nkp, block_store, efermi, wtk, stm_broad, &
         nxmin, nymin, nzmin, current, nptsx, nptsy, nptsz, &
         eigenvalues, band_no, stm_bias, charge_stub, n_bands_active, band_no, nkp, flag_by_kpoint, &
         flag_output, dx, cube, gpv
    use block_module, only: n_pts_in_block, in_block_x,in_block_y,in_block_z
    use io_module, ONLY: get_file_name
    use output, ONLY: write_dx_density, write_cube, write_dx_coords
    use global_module, only : nspin

    implicit none

    integer :: proc, ispin
    real(double) :: n_elect
    character(len=50) :: filename, ci

    allocate(current(nptsx,nptsy,nptsz))
    do ispin=1,nspin
       current = zero
       if(nspin==1) then
          ci = TRIM(charge_stub)
       else
          if(ispin==1) ci = TRIM(charge_stub)//"_up"
          if(ispin==2) ci = TRIM(charge_stub)//"_dn"
       end if
       do proc = 1, nprocs
       do proc = 1, nprocs
          call get_file_name(ci,nprocs,proc,filename)
          ! Open file
          open(unit=17,file=filename)
          call read_domain(17,proc,current)
          close(unit=17)
       end do ! proc
       call write_cube(current,ci)
       n_elect = sum(current)*gpv
       write(*,*) 'Number of electrons in active area: ',n_elect
    end do ! ispin
    return
  end subroutine process_charge

  subroutine process_bands

    use datatypes
    use numbers
    use local, ONLY: nkp, efermi, current, nptsx, nptsy, nptsz, eigenvalues, flag_by_kpoint, &
         n_bands_total, band_active_kp, flag_proc_range, band_full_to_active, evec_coeff,&
         E_procwf_min, E_procwf_max, flag_procwf_range_Ef, band_proc_no, n_bands_process, &
         grid_x, grid_y, grid_z, wtk, flag_outputWF_real
    use output, ONLY: write_cube
    use global_module, only : nspin
    use read, ONLY: read_eigenvalues, read_psi_coeffs
    use pao_format, ONLY: pao
    use global_module, ONLY: ni_in_cell, atom_coord, species_glob
    use species_module, ONLY: nsf_species, n_species
    use angular_coeff_routines, ONLY: set_prefac_real, set_fact, set_prefac
    use GenComms, ONLY: cq_abort

    implicit none

    integer :: proc, band, nk, idum1, idum2, kp, ispin, i, j, k, ipao, jpao
    real(double) :: weight, rbx, rby, rbz, sq, test, gpv, integral
    real(double), dimension(2) :: Emin, Emax
    character(len=50) :: filename, ci
    complex(double_cplx), dimension(:,:,:), allocatable :: psi
    integer :: i_atom, i_spec, i_l, i_zeta,i_m,j_atom,j_spec,j_l,j_zeta,j_m
    integer :: i_band, i_pao, j_pao
    real(double), parameter :: band_integral_tol = 1e-3_double
    real(double) :: max_band_integral_deviation, integral_deviation
    ! Create arrays needed by Conquest PAO routines
    call set_fact(8)
    call set_prefac(9)
    call set_prefac_real(9)
    gpv = grid_x*grid_y*grid_z
    ! Read eigenvalues
    call read_eigenvalues
    ! Read eigenvector coefficients
    call read_psi_coeffs("Process")
    allocate(current(nptsx,nptsy,nptsz))
    allocate(psi(nptsx,nptsy,nptsz))
    max_band_integral_deviation = zero

    if (flag_outputWF_real .and. (nkp.ne.1)) &
       call cq_abort("OutputWF_real is available only for Gamma-point calculations.")

    if(flag_proc_range) then
       Emin = E_procwf_min
       Emax = E_procwf_max
       if(flag_by_kpoint) then
          write(*,fmt='(4x,"Writing bands at each k-point")')
       else
          write(*,fmt='(4x,"Summing over k-points")')
       end if
       if(flag_procwf_range_Ef) then
          Emin = efermi + Emin
          Emax = efermi + Emax
       end if
       write(*,fmt='(4x,"Writing bands between ",e12.4," and ",e12.4,"Ha as specified in input file")') &
            Emin(1),Emax(1)
       do ispin=1,nspin
          if(flag_by_kpoint) then ! Separate bands by k-point
             do band=1,n_bands_total
                current = zero
                do kp = 1,nkp
                   if(eigenvalues(band,kp,ispin) >= Emin(ispin) .and. &
                        eigenvalues(band,kp,ispin) <= Emax(ispin) .and. &
                        band_active_kp(band,kp,ispin)==1) then
                      psi = zero
                      current = zero
                      call pao_to_grid(band_full_to_active(band), kp, ispin, psi)
                      if (.not.flag_outputWF_real) then
                         current = psi*conjg(psi)   ! band density
                      else
                         current = real(psi)        ! only real-part of WF
                      endif
                      write(ci,'("Band",I0.6,"den_kp",I0.3,"S",I0.1)') band, kp, ispin
                      call write_cube(current,ci)
                      integral = gpv*sum(current)
                      integral_deviation = abs(integral - one)
                      max_band_integral_deviation = max(integral_deviation, max_band_integral_deviation)
                      ! Check for problems with band integral
                      if(integral_deviation>band_integral_tol) &
                           write(*,fmt='(4x,"Integral of band ",i5," with energy ",f17.10," is ",f17.10)') &
                           band,eigenvalues(band,kp,ispin),integral
                   end if
                end do ! kp
             end do ! bands = 1, n_bands_total
          else ! Sum over k-points
             do band=1,n_bands_total
                current = zero
                idum1=0
                do kp = 1,nkp
                   if(eigenvalues(band,kp,ispin) >= Emin(ispin) .and. &
                        eigenvalues(band,kp,ispin) <= Emax(ispin) .and. &
                        band_active_kp(band,kp,ispin)==1) then
                      call pao_to_grid(band_full_to_active(band), kp, ispin, psi)
                      integral = gpv*sum(psi*conjg(psi))
                      integral_deviation = abs(integral - one)
                      max_band_integral_deviation = max(integral_deviation, max_band_integral_deviation)
                      ! Check for problems with band integral
                      if(integral_deviation>band_integral_tol) &
                           write(*,fmt='(4x,"Integral of band ",i5," at kp ",i5," is ",f17.10)') &
                           band,kp,integral
                      if (.not.flag_outputWF_real) then
                         current = current + psi*conjg(psi)*wtk(kp)   ! band density
                      else
                         current = current + real(psi)*wtk(kp)        ! only real-part of WF
                      endif
                      idum1 = 1
                   end if
                end do ! kp
                if(idum1==1) then
                   write(ci,'("Band",I0.6,"den_totS",I0.1)') band, ispin
                   call write_cube(current,ci)
                   integral = gpv*sum(current)
                   integral_deviation = abs(integral - one)
                   max_band_integral_deviation = max(integral_deviation, max_band_integral_deviation)
                   ! Check for problems with band integral
                   if(integral_deviation>band_integral_tol) &
                        write(*,fmt='(4x,"Integral of band ",i5," is ",f17.10)') &
                        band,integral
                end if
             end do ! bands
          end if
       end do
    else ! User has provided list of bands
       write(*,fmt='(4x,"Writing ",i4," bands specified in input file")') n_bands_process
       if(flag_by_kpoint) then
          write(*,fmt='(4x,"Writing bands at each k-point")')
       else
          write(*,fmt='(4x,"Summing over k-points")')
       end if
       do ispin=1,nspin
          if(flag_by_kpoint) then ! Separate bands by k-point
             do band=1,n_bands_process
                current = zero
                do kp = 1,nkp
                   ! This clause is needed in case the user chose an energy range that only selects some k-points
                   if(band_active_kp(band_proc_no(band),kp,ispin)==1) then
                      psi = zero
                      current = zero
                      call pao_to_grid(band_full_to_active(band_proc_no(band)), kp, ispin, psi)
                      if (.not.flag_outputWF_real) then
                         current = psi*conjg(psi)   ! band density
                      else
                         current = real(psi)        ! only real-part of WF
                      endif
                      write(ci,'("Band",I0.6,"den_kp",I0.3,"S",I0.1)') band_proc_no(band), kp, ispin
                      call write_cube(current,ci)
                      integral = gpv*sum(current)
                      integral_deviation = abs(integral - one)
                      max_band_integral_deviation = max(integral_deviation, max_band_integral_deviation)
                      ! Check for problems with band integral
                      if(integral_deviation>band_integral_tol) &
                           write(*,fmt='(4x,"Integral of psi squared ",i5," with energy ",f17.10," is ",f17.10)') &
                           band,eigenvalues(band,kp,ispin),integral
                   end if
                end do ! kp
             end do ! bands = 1, n_bands_total
          else ! Sum over k-points
             do band=1,n_bands_process
                current = zero
                do kp = 1,nkp
                   ! This clause is needed in case the user chose an energy range that only selects some k-points
                   if(band_active_kp(band_proc_no(band),kp,ispin)==1) then
                      call pao_to_grid(band_full_to_active(band_proc_no(band)), kp, ispin, psi)
                      integral = gpv*sum(psi*conjg(psi))
                      integral_deviation = abs(integral - one)
                      max_band_integral_deviation = max(integral_deviation, max_band_integral_deviation)
                      ! Check for problems with band integral
                      if(integral_deviation>band_integral_tol) &
                           write(*,fmt='(4x,"Integral of band ",i5," at kp ",i5," is ",f17.10)') &
                           band,kp,integral
                      if (.not.flag_outputWF_real) then
                         current = current + psi*conjg(psi)*wtk(kp)   ! band density
                      else
                         current = current + real(psi)*wtk(kp)        ! only real-part of WF
                      endif
                   end if
                end do ! kp
                write(ci,'("Band",I0.6,"den_totS",I0.1)') band_proc_no(band), ispin
                call write_cube(current,ci)
                integral = gpv*sum(current)
                integral_deviation = abs(integral - one)
                max_band_integral_deviation = max(integral_deviation, max_band_integral_deviation)
                ! Check for problems with band integral
                if(integral_deviation>band_integral_tol) &
                     write(*,fmt='(4x,"Integral of band ",i5," is ",f17.10)') &
                     band,integral
             end do ! bands
          end if
       end do
    end if
    write(*,fmt='(4x,"Largest deviation of band integral from one is ",f8.5)') max_band_integral_deviation
    return
  end subroutine process_bands

  subroutine process_dos

    use datatypes
    use numbers, ONLY: zero, RD_ERR, twopi, half, one, two, four, six
    use local, ONLY: eigenvalues, n_bands_total, nkp, wtk, efermi, &
         flag_total_iDOS, flag_procwf_range_Ef, flag_expand_range
    use read, ONLY: read_eigenvalues, read_psi_coeffs
    use global_module, ONLY: nspin, n_DOS, E_DOS_min, E_DOS_max, sigma_DOS, &
         flag_fix_spin_population
    use units, ONLY: HaToeV

    implicit none

    ! Local variables
    integer :: i_band, i_kp, i_spin, n_DOS_wid, n_band, n_min, n_max, i
    real(double) :: Ebin, dE_DOS, a, pf_DOS, spin_fac, peak_width
    real(double), dimension(nspin) :: total_electrons, iDOS_low
    real(double), dimension(:,:), allocatable :: total_DOS, iDOS
    real(double), dimension(:,:), allocatable :: occ

    write(*,fmt='(/2x,"Calculating density of states (DOS)")')
    if(nspin==1) then
       spin_fac = two
    else if(nspin==2) then
       spin_fac = one
    end if
    ! Read eigenvalues
    call read_eigenvalues
    allocate(total_DOS(n_DOS,nspin),iDOS(n_DOS,nspin),occ(n_bands_total,nkp))
    total_DOS = zero
    iDOS = zero
    ! Define peak width
    peak_width = six*sigma_DOS
    ! Offset for the energy range used for DOS display
    range_offset = zero
    do i_spin = 1, nspin
       if(flag_procwf_range_Ef) then
          if(nspin>1.and.flag_fix_spin_population) then
             range_offset(i_spin) = efermi(1)
          else
             range_offset(i_spin) = efermi(i_spin)
          end if
       end if
    end do
    ! Set limits and broaden energy range if needed
    if(abs(E_DOS_min)<RD_ERR) then
       E_DOS_min = minval(eigenvalues(1,:,:))
       if(flag_expand_range) E_DOS_min = E_DOS_min - peak_width
    else
       if(flag_expand_range) E_DOS_min = E_DOS_min - peak_width
    end if
    if(abs(E_DOS_max)<RD_ERR) then
       E_DOS_max = maxval(eigenvalues(n_bands_total,:,:))
       if(flag_expand_range) E_DOS_max = E_DOS_max + peak_width
    else
       if(flag_expand_range) E_DOS_max = E_DOS_max + peak_width
    end if
    ! Spacing, width, prefactor
    dE_DOS = (E_DOS_max - E_DOS_min)/real(n_DOS-1,double)
    write(*,fmt='(2x,"Dividing DOS into ",i5," bins of width ",f12.6," Ha")') n_DOS, dE_DOS
    write(*,fmt='(2x,"Sigma is: ",f12.6," Ha")') sigma_DOS
    if(peak_width < dE_DOS) write(*,fmt='(4x,"Sigma is much less than bin size: this may cause errors")')
    n_DOS_wid = floor(peak_width/dE_DOS) ! How many bins either side of state we consider
    pf_DOS = one/(sigma_DOS*sqrt(twopi))
    total_electrons = zero
    iDOS_low = zero ! Integral of DOS to lowest energy bound
    ! Accumulate DOS over bands and k-points for each spin
    do i_spin = 1, nspin
       write(*,fmt='(2x,"DOS lower limit: ",f12.5," Ha")') E_DOS_min + range_offset(i_spin)
       write(*,fmt='(2x,"DOS upper limit: ",f12.5," Ha")') E_DOS_max + range_offset(i_spin)
       write(*,fmt='(2x,"Fermi level:     ",f12.5," Ha")') efermi(i_spin)
       occ = zero
       call occupy(occ,eigenvalues,efermi,i_spin)
       do i_kp = 1, nkp
          do i_band=1,n_bands_total ! All bands
             if(eigenvalues(i_band, i_kp, i_spin)>=E_DOS_min + range_offset(i_spin) .and. &
                  eigenvalues(i_band, i_kp, i_spin)<=E_DOS_max + range_offset(i_spin)) then
                n_band = floor((eigenvalues(i_band, i_kp, i_spin) - (E_DOS_min + range_offset(i_spin)))/dE_DOS) + 1
                n_min = n_band - n_DOS_wid
                if(n_min<1) n_min = 1
                n_max = n_band + n_DOS_wid
                if(n_max>n_DOS) n_max = n_DOS
                do i = n_min, n_max
                   Ebin = real(i-1,double)*dE_DOS + E_DOS_min + range_offset(i_spin)
                   a = (Ebin-eigenvalues(i_band, i_kp, i_spin))/sigma_DOS
                   total_DOS(i,i_spin) = total_DOS(i,i_spin) + wtk(i_kp)*pf_DOS*exp(-half*a*a)
                   total_electrons(i_spin) = total_electrons(i_spin) + occ(i_band,i_kp)*wtk(i_kp)*pf_DOS*exp(-half*a*a)
                   iDOS(i,i_spin) = iDOS(i,i_spin) + wtk(i_kp)*pf_DOS*exp(-half*a*a)
                end do
             else if(eigenvalues(i_band, i_kp, i_spin)<E_DOS_min + range_offset(i_spin)) then
                iDOS_low(i_spin) = iDOS_low(i_spin) + wtk(i_kp)
             end if
          end do
       end do
       ! Now integrate DOS
       !write(*,*) 'iDOS_low is ',iDOS_low
       do i = 2, n_DOS
          iDOS(i,i_spin) = iDOS(i,i_spin) + iDOS(i-1,i_spin)
       end do
    end do
    ! Include spin factor
    iDOS = iDOS*dE_DOS*spin_fac
    if(flag_total_iDOS) then
       do i_spin = 1, nspin
          iDOS(:,i_spin) = iDOS(:,i_spin) + spin_fac*iDOS_low(i_spin)
       end do
    end if
    total_electrons = total_electrons*dE_DOS*spin_fac
    total_DOS = total_DOS*spin_fac
    write(*,fmt='(2x,"Integration of DOS in terms of absolute energies")')
    if(nspin==1) then
       write(*,fmt='(2x,"DOS between ",f11.3," Ha and Ef integrates to ",f12.3," electrons")') &
            E_DOS_min + range_offset(1), total_electrons(1)
       write(*,fmt='(2x,"DOS between ",f11.3," Ha and ",f11.3," Ha integrates to ",f12.3," electrons")') &
            E_DOS_min + range_offset(1), E_DOS_max + range_offset(1),dE_DOS*sum(total_DOS(:,1))
    else
       write(*,fmt='(2x,"Spin Up DOS integrated to Ef gives ",f12.3," electrons")') total_electrons(1)
       write(*,fmt='(2x,"Spin Dn DOS integrated to Ef gives ",f12.3," electrons")') total_electrons(2)
       write(*,fmt='(2x,"Spin Up DOS between ",f11.3," Ha and ",f11.3," Ha integrates to ",f12.3," electrons")') &
            E_DOS_min + range_offset(1), E_DOS_max + range_offset(1),dE_DOS*sum(total_DOS(:,1))
       write(*,fmt='(2x,"Spin Up DOS between ",f11.3," Ha and ",f11.3," Ha integrates to ",f12.3," electrons")') &
            E_DOS_min + range_offset(2), E_DOS_max + range_offset(2),dE_DOS*sum(total_DOS(:,2))
    end if
    ! Since we write out DOS against eV we need this conversion to get the integral right
    total_DOS = total_DOS/HaToeV
    ! Write out DOS, shifted to Ef = 0
    open(unit=17, file="DOS.dat")
    do i_spin = 1, nspin
       write(17,fmt='("# Spin ",I1)') i_spin
       write(17,fmt='("# Original Fermi-level: ",f12.5," eV")') HaToeV*efermi(i_spin)
       write(17,fmt='("# DOS shifted relative to Fermi-level")')
       if(flag_total_iDOS) then
          write(17,fmt='("#  Energy (eV)     DOS (/eV)    Total iDOS")')
       else
          write(17,fmt='("#  Energy (eV)     DOS (/eV)    Local iDOS")')
       end if
       do i=1, n_DOS
          write(17,fmt='(3f14.5)') HaToeV*(E_DOS_min + range_offset(i_spin) - efermi(i_spin) + dE_DOS*real(i-1,double)), &
               total_DOS(i,i_spin), iDOS(i,i_spin)
       end do
       write(17,fmt='("&")')
    end do
    close(unit=17)
    deallocate(total_DOS,iDOS,occ)
    return
  end subroutine process_dos

  ! Initially produce DOS projected onto all atoms
  subroutine process_pdos

    use datatypes
    use numbers, ONLY: zero, RD_ERR, pi, twopi, half, one, two, four, six
    use local, ONLY: eigenvalues, n_bands_total, nkp, wtk, efermi, flag_total_iDOS, &
         evec_coeff, scaled_evec_coeff, flag_procwf_range_Ef, flag_l_resolved, flag_lm_resolved, &
         flag_rotate_pdos, flag_rotate_pdos_mode,flag_rotate_pdos_debug,  band_full_to_active, n_atoms_pDOS, pDOS_atom_index, &
         pdos_ax, pdos_ay, pdos_az, rotate_pdos_atoms_euler, euler_angles, find_neighbours, rotate_pdos_natoms, &
         nghbr_arr, U1, U2
    use read, ONLY: read_eigenvalues, read_psi_coeffs, read_nprocs_from_blocks
    use global_module, ONLY: nspin, n_DOS, E_DOS_min, E_DOS_max, sigma_DOS, ni_in_cell, species_glob
    use units, ONLY: HaToeV
    use species_module, ONLY: nsf_species, npao_species
    use pao_format, ONLY: pao

    implicit none

    ! Local variables
    integer :: i_band, i_kp, i_spin, n_DOS_wid, n_band, n_min, n_max, i, j, k, i_atom,max_nsf, i_spec, &
         i_l, nzeta, sf_offset, max_l, norbs, i_m, i_band_c, i_z
    real(double) :: Ebin, dE_DOS, a, pf_DOS, spin_fac, coeff, check_electrons, peak_width
   integer, dimension(:), allocatable :: nghbr_atoms
    real(double), dimension(:,:), allocatable :: bond ! For pDOS rotation
    real(double), dimension(:,:,:), allocatable :: pDOS
    real(double), dimension(:,:,:,:), allocatable :: pDOS_l
    real(double), dimension(:,:,:,:,:), allocatable :: pDOS_lm
    real(double), dimension(:,:), allocatable :: occ
    real(double), dimension(:,:), allocatable :: total_electrons
    real(double), dimension(:,:,:), allocatable :: total_electrons_l
    real(double) :: A1(3, 3), A2(5, 5), C1(3, 3), C2(5, 5), E1(3,3), E2(5,5), Qxyz(3,3)
    real(double) :: rod(9), angle, axis(3), temp_matrix(3,3)
    character(len=25) :: filename,fmt_dos
    complex(double_cplx),external :: zdotc

    write(*,fmt='(/2x,"Calculating projected density of states (DOS)")')
    if(flag_l_resolved .and. flag_lm_resolved) then
       write(*,fmt='(4x,"Resolving in l and m")')
    else if(flag_l_resolved) then
       write(*,fmt='(4x,"Resolving in l")')
    end if
    if(n_atoms_pDOS==ni_in_cell) then
       write(*,fmt='(4x,"Producing pDOS for all atoms in cell")')
    else
       if(n_atoms_pDOS==1) then
          write(*,fmt='(4x,"Producing pDOS for ",i7," atom in cell")') n_atoms_pDOS
       else
          write(*,fmt='(4x,"Producing pDOS for ",i7," atoms in cell")') n_atoms_pDOS
       end if
    end if
    call read_nprocs_from_blocks
    if(nspin==1) then
       spin_fac = two
    else if(nspin==2) then
       spin_fac = one
    end if
    ! DOS processing called first, so eigenvalues already read
    ! Read eigenvector coefficients scaled by Sij into variable evec_coeff
    call read_psi_coeffs("ProcessSij")
    max_nsf = maxval(nsf_species)
    max_l = maxval(pao(:)%greatest_angmom)
    ! The subroutine read_psi_coeffs allocates evec_coeff, so we make a copy and deallocate
    allocate(scaled_evec_coeff(max_nsf, ni_in_cell, n_bands_total, nkp, nspin))
    scaled_evec_coeff = evec_coeff
    deallocate(evec_coeff)
    ! Read eigenvector coefficients
    call read_psi_coeffs("Process")
    allocate(occ(n_bands_total,nkp))

    ! Call pDOS rotation subroutines if desired. Also includes debug output
    if (flag_rotate_pdos) then
      write(*,fmt='(2x,"Rotating wavefunction coefficients")')
      call initialise_A_mat(A1, A2)
      if (flag_rotate_pdos_debug .and. flag_rotate_pdos_mode /= 1) then
         write(*, fmt='(/2x,"ROTATION DEBUG OUTPUT: MODE ", (I0,1x))') &
            flag_rotate_pdos_mode

         write(*, fmt='(/2x,"ROTATION DEBUG OUTPUT: A^(l) MATRICES")')
         write(*, fmt='(/4x, "A1: ")')
         do j = 1, 3
            write(*, fmt='(/4x,3(f10.5,1X))') A1(j, :)
         end do
         write(*, fmt='(/4x, "A2: ")')
         do j = 1, 5
            write(*, fmt='(/4x,5(f10.5,1X))') A2(j, :)
         end do
         write(*, fmt='(/4x, "For l = 1, orbital basis is [|y>, |z>, |x>]")')
         write(*, fmt='(/4x, "For l = 2, orbital basis is [|xy>, |yz>, |3z^2-r^2>, |xz>, |x^2-y^2>]")')
      end if
      if (allocated(U1))  deallocate(U1)
      if (allocated(U2))  deallocate(U2)

      if (flag_rotate_pdos_mode == 0) then
         write(*,fmt='(2x,"Using user input axes")')
         allocate(U1(3,3,n_atoms_pDOS))
         allocate(U2(5,5,n_atoms_pDOS))
         call get_pdos_axes

         if(n_atoms_pDOS==ni_in_cell) then ! All atoms
            write(*, fmt='(/2x,"New local axes for all atoms in unit cell:")')
         else
            write(*, fmt='(/2x,"New local axes for specified atoms: ", *(I0,1x))') &
               pDOS_atom_index
         end if
         write(*, fmt='(/4x,"x: ",3(f10.5,1X))', advance="no") pdos_ax
         write(*, fmt='(/4x,"y: ",3(f10.5,1X))', advance="no") pdos_ay
         write(*, fmt='(/4x,"z: ",3(f10.5,1X))', advance="no") pdos_az
         call calculate_axis_angle(pdos_ax, pdos_ay, pdos_az, axis, angle)
         call construct_rodrigues(axis, angle, rod)
         ! Construct all C^l matrices from rodrigues
         call construct_C1(rod, C1)
         call construct_C2(rod, C2)
         do i = 1, n_atoms_pDOS
            call construct_Ul(1, A1, C1, U1(:,:,i))
            call construct_Ul(2, A2, C2, U2(:,:,i))
            U1(:,:,i) = transpose(U1(:,:,i))
            U2(:,:,i) = transpose(U2(:,:,i))

            if (flag_rotate_pdos_debug) then
               call euler_from_axisangle(axis, angle)
               write(*, fmt='(/4x, "C1: ")')
               do j = 1, 3
                  write(*, fmt='(/6x,3(f10.5,1X))',  advance='no') C1(j,:)
               end do
               write(*, fmt='(/4x, "C2: ")')
               do j = 1, 5
                  write(*, fmt='(/6x,5(f10.5,1X))',  advance='no') C2(j, :)
               end do
               call print_orbital_weights(i)
            end if ! end rotation debug output
          end do

      else if (flag_rotate_pdos_mode == 1) then
         write(*,fmt='(2x,"Using extrinsic Euler angles in active zyz convention")')
         allocate(U1(3,3,rotate_pdos_natoms))
         allocate(U2(5,5,rotate_pdos_natoms))
         Qxyz = 0.0
         Qxyz(1,2) = 1.0
         Qxyz(2,3) = 1.0
         Qxyz(3,1) = 1.0
         do i = 1, rotate_pdos_natoms
            write(*, fmt='(/2x,"New local axes for specified atoms: ", *(I0,1x))') &
               rotate_pdos_atoms_euler(i)
            call construct_EulerMatrices(E1, E2, i)
            U1(:,:,i) = (E1)
            U2(:,:,i) = (E2)
            ! We only explicitly construct wavefunction rotation matrices
            ! However can extrapolate new local axes using E1 since orbital basis is {y,z,x}
            write(*, fmt='(/4x,"Image of x,y,z under this active rotation, in standard basis x,y,z")', advance="no")
            write(*, fmt='(/4x,"x: ",3(f10.5,1X))', advance="no") &
            E1(3,3), E1(1,3), E1(2,3)
            write(*, fmt='(/4x,"y: ",3(f10.5,1X))', advance="no") &
            E1(3,1), E1(1,1), E1(2,1)
            write(*, fmt='(/4x,"z: ",3(f10.5,1X))', advance="no") &
            E1(3,2), E1(1,2), E1(2,2)
            temp_matrix = matmul(inv(Qxyz), matmul(E1, Qxyz))
            temp_matrix = inv(temp_matrix)
             write(*, fmt='(/4x,"Equivalent basis transformation (inverse of active coordinates)")', advance="no")
             write(*, fmt='(/4x,"x: ",3(f10.5,1X))', advance="no") &
               temp_matrix(:,1)
            write(*, fmt='(/4x,"y: ",3(f10.5,1X))', advance="no") &
               temp_matrix(:,2)
            write(*, fmt='(/4x,"z: ",3(f10.5,1X))', advance="no") &
               temp_matrix(:,3)
            if (flag_rotate_pdos_debug) then
               write(*, fmt='(/4x, "alpha(z) beta(y) gamma(z): ",3(f10.5,1X))') &
                  euler_angles(1,i), euler_angles(2,i),euler_angles(3,i)
               call print_orbital_weights(i)
            end if ! end rotation debug mode
         end do
      else if (flag_rotate_pdos_mode == 2) then
         write(*,fmt='(2x,"Using user input atom numbers and local geometry")')
         allocate(U1(3,3,rotate_pdos_natoms))
         allocate(U2(5,5,rotate_pdos_natoms))

         do i = 1, rotate_pdos_natoms
            call nearest_neighbours(find_neighbours(1, i), bond)
            write(*, fmt='(/2x,"Located neighours of atom ", I0, ": ", *(I0,1x))') &
                    find_neighbours(1,i), nghbr_arr
            call axes_from_nn(find_neighbours(1, i), bond)
            write(*, fmt='(/2x,"New local axes for atom ", I0, ": ")') &
                    find_neighbours(1,i)
            write(*, fmt='(/4x,"x: ",3(f10.5,1X))', advance="no") pdos_ax
            write(*, fmt='(/4x,"y: ",3(f10.5,1X))', advance="no") pdos_ay
            write(*, fmt='(/4x,"z: ",3(f10.5,1X))', advance="no") pdos_az
            call get_pdos_axes
            call calculate_axis_angle(pdos_ax, pdos_ay, pdos_az, axis, angle)
            call construct_rodrigues(axis, angle, rod)
            ! Construct all C^l matrices from rodrigues
            call construct_C1(rod, C1)
            call construct_C2(rod, C2)
            call construct_Ul(1, A1, C1, U1(:,:,i))
            call construct_Ul(2, A2, C2, U2(:,:,i))
            U1(:,:,i) = transpose(U1(:,:,i))
            U2(:,:,i) = transpose(U2(:,:,i))
            if (flag_rotate_pdos_debug) then
               call euler_from_axisangle(axis, angle)
               write(*, fmt='(/4x, "C1: ")')
               do j = 1, 3
                  write(*, fmt='(/6x,3(f10.5,1X))',  advance='no') C1(j,:)
               end do
               write(*, fmt='(/4x, "C2: ")')
               do j = 1, 5
                  write(*, fmt='(/6x,5(f10.5,1X))',  advance='no') C2(j, :)
               end do
               call print_orbital_weights(i)
            end if ! end rotation debug output
         end do
      end if ! pdos rotation mode
      call rotate_coefficients
      deallocate(U1)
      deallocate(U2)
    end if  ! rotate pdos
    ! Set up storage based on pDOS per atom, or l/lm resolved per atom
    if(flag_lm_resolved) then
       allocate(pDOS_lm(-max_l:max_l,0:max_l,n_atoms_pDOS,n_DOS,nspin))
       pDOS_lm = zero
       allocate(total_electrons_l(0:max_l,n_atoms_pDOS,nspin))
       total_electrons_l = zero
       ! For total pDOS
       allocate(pDOS(n_atoms_pDOS,n_DOS,nspin))
       pDOS = zero
    else if(flag_l_resolved) then
       allocate(pDOS_l(0:max_l,n_atoms_pDOS,n_DOS,nspin))
       pDOS_l = zero
       allocate(total_electrons_l(0:max_l,n_atoms_pDOS,nspin))
       total_electrons_l = zero
       ! For total pDOS
       allocate(pDOS(n_atoms_pDOS,n_DOS,nspin))
       pDOS = zero
    else
       allocate(pDOS(n_atoms_pDOS,n_DOS,nspin))
       pDOS = zero
       allocate(total_electrons(n_atoms_pDOS,nspin))
       total_electrons = zero
    end if
    ! E_DOS_min and max and sigma_DOS already set in process_dos
    ! Spacing, width, prefactor
    peak_width = six*sigma_DOS
    dE_DOS = (E_DOS_max - E_DOS_min)/real(n_DOS-1,double)
    n_DOS_wid = floor(peak_width/dE_DOS) ! How many bins either side of state we consider
    pf_DOS = one/(sigma_DOS*sqrt(twopi))
    ! Accumulate DOS over bands and k-points for each spin
    do i_spin = 1, nspin
       occ = zero
       call occupy(occ,eigenvalues,efermi,i_spin)
       if(flag_procwf_range_Ef) then
          E_DOS_min = E_DOS_min + efermi(i_spin)
          E_DOS_max = E_DOS_max + efermi(i_spin)
       end if
       do i_kp = 1, nkp
          do i_band=1,n_bands_total ! All bands
             if(eigenvalues(i_band, i_kp, i_spin)>E_DOS_min + range_offset(i_spin) .and. &
                  eigenvalues(i_band, i_kp, i_spin)<E_DOS_max + range_offset(i_spin)) then
                i_band_c = band_full_to_active(i_band)
                n_band = floor((eigenvalues(i_band, i_kp, i_spin) - (E_DOS_min + range_offset(i_spin)))/dE_DOS) + 1
                n_min = n_band - n_DOS_wid
                if(n_min<1) n_min = 1
                n_max = n_band + n_DOS_wid
                if(n_max>n_DOS) n_max = n_DOS
                do i = n_min, n_max
                   Ebin = real(i-1,double)*dE_DOS + E_DOS_min + range_offset(i_spin)
                   a = (Ebin-eigenvalues(i_band, i_kp, i_spin))/sigma_DOS
                   do i_atom = 1, n_atoms_pDOS
                      i_spec = species_glob(pDOS_atom_index(i_atom))
                      if(flag_l_resolved .and. flag_lm_resolved) then
                         sf_offset = 1
                         do i_l = 0, pao(i_spec)%greatest_angmom
                            nzeta = pao(i_spec)%angmom(i_l)%n_zeta_in_angmom
                            do i_z = 1, nzeta
                               do i_m = -i_l,i_l
                                  coeff = conjg(evec_coeff(sf_offset,pDOS_atom_index(i_atom), i_band_c,i_kp,i_spin)) * &
                                       scaled_evec_coeff(sf_offset,pDOS_atom_index(i_atom), i_band_c,i_kp,i_spin)
                                  pDOS_lm(i_m,i_l,i_atom,i,i_spin) = &
                                       pDOS_lm(i_m,i_l,i_atom,i,i_spin) + &
                                       wtk(i_kp)*pf_DOS*exp(-half*a*a)*coeff
                                  pDOS(i_atom,i,i_spin) = pDOS(i_atom,i,i_spin) + &
                                       wtk(i_kp)*pf_DOS*exp(-half*a*a)*coeff
                                  total_electrons_l(i_l,i_atom, i_spin) = &
                                       total_electrons_l(i_l,i_atom, i_spin) + &
                                       occ(i_band,i_kp)*wtk(i_kp)*pf_DOS*exp(-half*a*a)*coeff
                                  sf_offset = sf_offset + 1
                               end do
                            end do
                         end do
                      else if(flag_l_resolved) then
                         sf_offset = 0
                         do i_l = 0, pao(i_spec)%greatest_angmom
                            nzeta = pao(i_spec)%angmom(i_l)%n_zeta_in_angmom
                            norbs = nzeta*(2*i_l+1)
                            coeff = zdotc(norbs, evec_coeff(sf_offset+1:sf_offset+norbs,pDOS_atom_index(i_atom), &
                                 i_band_c,i_kp,i_spin),1, &
                                 scaled_evec_coeff(sf_offset+1:sf_offset+norbs,pDOS_atom_index(i_atom), &
                                 i_band_c,i_kp,i_spin),1)
                            pDOS_l(i_l,i_atom,i,i_spin) = pDOS_l(i_l,i_atom,i,i_spin) + &
                                 wtk(i_kp)*pf_DOS*exp(-half*a*a)*coeff
                            pDOS(i_atom,i,i_spin) = pDOS(i_atom,i,i_spin) + wtk(i_kp)*pf_DOS*exp(-half*a*a)*coeff
                            total_electrons_l(i_l,i_atom, i_spin) = total_electrons_l(i_l,i_atom, i_spin) + &
                                 occ(i_band,i_kp)*wtk(i_kp)*pf_DOS*exp(-half*a*a)*coeff
                            sf_offset = sf_offset + norbs
                         end do
                      else
                         coeff = zdotc(npao_species(i_spec),evec_coeff(1:npao_species(i_spec), &
                              pDOS_atom_index(i_atom),i_band_c,i_kp,i_spin),1, &
                              scaled_evec_coeff(1:npao_species(i_spec), &
                              pDOS_atom_index(i_atom),i_band_c,i_kp,i_spin),1)
                         pDOS(i_atom,i,i_spin) = pDOS(i_atom,i,i_spin) + wtk(i_kp)*pf_DOS*exp(-half*a*a)*coeff
                         total_electrons(i_atom, i_spin) = total_electrons(i_atom, i_spin) + &
                              occ(i_band,i_kp)*wtk(i_kp)*pf_DOS*exp(-half*a*a)*coeff
                      end if
                   end do
                end do
             end if
          end do
       end do
    end do ! do i_spin = 1, n_spin
    ! Include spin factor and convert Ha to eV
    if(flag_l_resolved .and. flag_lm_resolved) then
       pDOS_lm = pDOS_lm*spin_fac/HaToeV
       pDOS = pDOS*spin_fac/HaToeV
       total_electrons_l = total_electrons_l*dE_DOS*spin_fac
    else if(flag_l_resolved) then
       pDOS_l = pDOS_l*spin_fac/HaToeV
       pDOS = pDOS*spin_fac/HaToeV
       total_electrons_l = total_electrons_l*dE_DOS*spin_fac
    else
       pDOS = pDOS*spin_fac/HaToeV
       total_electrons = total_electrons*dE_DOS*spin_fac
    end if
    write(*,fmt='(2x,"Integration of DOS in terms of absolute energies")')
    if(nspin==1) then
       !write(*,fmt='(2x,"Results of integrating pDOS between ",f11.3," and ",f11.3," Ha (electrons per atom).")') &
       write(*,fmt='(2x,"Results of integrating pDOS between ",f11.3," and Ef (electrons per atom).")') &
            E_DOS_min + range_offset(1)
       if(flag_l_resolved) then ! l and l-m
          write(*,fmt='(4x,"   Atom       Total        l=0        l=1        l=2")')
          write(fmt_DOS,*) max_l + 2 ! Number of columns
          fmt_DOS = '(4x,i7,x,'//trim(adjustl(fmt_DOS))//'f11.3)'
          do i_atom = 1, n_atoms_pDOS
             write(*,fmt=fmt_DOS) pDOS_atom_index(i_atom), sum(total_electrons_l(:,i_atom,1)),total_electrons_l(:,i_atom,1)
          end do
          write(*,fmt='(2x,"Integrated pDOS: ",f11.3," electrons")') sum(total_electrons_l)
          if(flag_lm_resolved) then
             norbs = 1
             do i_l=0,max_l
                norbs = norbs + 2*i_l + 1
             end do
             norbs = norbs + 1 ! Total pDOS column
             write(fmt_DOS,*) norbs ! Number of columns
             fmt_DOS = '('//trim(adjustl(fmt_DOS))//'f12.5)'
          else
             write(fmt_DOS,*) max_l + 3 ! Number of columns (extra columns for energy and total pDOS)
             fmt_DOS = '('//trim(adjustl(fmt_DOS))//'f12.5)'
          end if
       else
          write(*,fmt='(4x,"   Atom   Electrons")')
          do i_atom = 1, n_atoms_pDOS
             write(*,fmt='(4x,i7,x,f11.3)') pDOS_atom_index(i_atom), total_electrons(i_atom,1)
          end do
          write(*,fmt='(2x,"Integrated pDOS: ",f11.3," electrons")') sum(total_electrons)
       end if
    else
       if(flag_l_resolved) then
          write(*,fmt='(2x,"Results of integrating pDOS between ",f11.3," and Ef (electrons per atom).")') &
               E_DOS_min + range_offset(1)
          write(*,fmt='(4x,"   Atom     Spin Up        l=0        l=1        l=2  Spin Down        l=0        l=1        l=2")')
          write(fmt_DOS,*) 2*(max_l + 2) ! Number of columns
          fmt_DOS = '(4x,i7,x,'//trim(adjustl(fmt_DOS))//'f11.3)'
          do i_atom = 1, n_atoms_pDOS
             write(*,fmt=fmt_DOS) pDOS_atom_index(i_atom), sum(total_electrons_l(:,i_atom,1)),total_electrons_l(:,i_atom,1), &
                  sum(total_electrons_l(:,i_atom,2)),total_electrons_l(:,i_atom,2)
          end do
          write(*,fmt='(2x,"Integrated spin up pDOS: ",f11.3," electrons")') sum(total_electrons_l(:,:,1))
          write(*,fmt='(2x,"Integrated spin dn pDOS: ",f11.3," electrons")') sum(total_electrons_l(:,:,2))
          if(flag_lm_resolved) then
             norbs = 1
             do i_l=0,max_l
                norbs = norbs + 2*i_l + 1
             end do
             norbs = norbs + 1 ! Total pDOS column
             write(fmt_DOS,*) norbs ! Number of columns
             fmt_DOS = '('//trim(adjustl(fmt_DOS))//'f12.5)'
          else
             write(fmt_DOS,*) max_l + 3 ! Number of columns
             fmt_DOS = '('//trim(adjustl(fmt_DOS))//'f12.5)'
          end if
       else
          write(*,fmt='(2x,"Results of integrating pDOS between ",f11.3," and Ef (electrons per atom).")') &
               E_DOS_min + range_offset(1)
          write(*,fmt='(4x,"   Atom     Spin Up  Spin Down")')
          do i_atom = 1, n_atoms_pDOS
             write(*,fmt='(4x,i7,x,2f11.3)') pDOS_atom_index(i_atom), total_electrons(i_atom,1), total_electrons(i_atom,2)
          end do
          write(*,fmt='(2x,"Integrated spin up pDOS: ",f11.3," electrons")') sum(total_electrons(:,1))
          write(*,fmt='(2x,"Integrated spin dn pDOS: ",f11.3," electrons")') sum(total_electrons(:,2))
       end if
    end if
    ! Write out DOS, shifted to Ef = 0
    do i_atom = 1, n_atoms_pDOS
       i_spec = species_glob(pDOS_atom_index(i_atom))
       if(flag_l_resolved .and. flag_lm_resolved) then
          write(filename,'("Atom",I0.7,"DOS_lm.dat")') pDOS_atom_index(i_atom)
       else if(flag_l_resolved) then
          write(filename,'("Atom",I0.7,"DOS_l.dat")') pDOS_atom_index(i_atom)
       else
          write(filename,'("Atom",I0.7,"DOS.dat")') pDOS_atom_index(i_atom)
       end if
       open(unit=17, file=filename)
       do i_spin = 1, nspin
          write(17,fmt='("# Spin ",I1)') i_spin
          write(17,fmt='("# Original Fermi-level: ",f12.5," eV")') HaToeV*efermi(i_spin)
          write(17,fmt='("# DOS shifted relative to Fermi-level")')
          if(flag_l_resolved .and. flag_lm_resolved) then
             write(17,fmt='("# Energy(eV)   pDOS(/eV)")')
             write(17,fmt='("#                  Total         l=0         l=1                                 l=2")')
             do i=1, n_DOS
                write(17,fmt=fmt_dos) HaToeV*(E_DOS_min + range_offset(i_spin) - efermi(i_spin) + &
                     dE_DOS*real(i-1,double)), pDOS(i_atom,i,i_spin), &
                     ((pDOS_lm(i_m,i_l,i_atom,i,i_spin),i_m=-i_l,i_l),i_l=0,pao(i_spec)%greatest_angmom)
             end do
          else if(flag_l_resolved) then
             write(17,fmt='("# Energy(eV)   pDOS(/eV)")')
             write(17,fmt='("#                  Total         l=0         l=1         l=2")')
             do i=1, n_DOS
                write(17,fmt=fmt_dos) HaToeV*(E_DOS_min + range_offset(i_spin) - efermi(i_spin) + &
                     dE_DOS*real(i-1,double)), pDOS(i_atom,i,i_spin), &
                     pDOS_l(0:pao(i_spec)%greatest_angmom,i_atom,i,i_spin)
             end do
          else
             write(17,fmt='("# Energy(eV)   pDOS(/eV)")')
             do i=1, n_DOS
                write(17,fmt='(2f12.5)') HaToeV*(E_DOS_min + range_offset(i_spin) - efermi(i_spin) + dE_DOS*real(i-1,double)), &
                     pDOS(i_atom,i,i_spin)
             end do
          end if
          write(17,fmt='("&")')
       end do
       close(unit=17)
    end do
    return
  end subroutine process_pdos

  subroutine print_orbital_weights(atom_index)
   use local, ONLY: U1, U2
   implicit none

   integer, intent(in) :: atom_index

   integer :: i,j,k
   real(double) :: pdos_weight, identity3(3,3), identity5(5,5)
   character(len=30), parameter :: p_orb(3) = (/ "|y>", "|z>","|x>"/)
   character(len=30), parameter :: d_orb(5) = (/ "|xy>      ", "|yz>      ", &
     "|3z^2-r^2>", "|xz>      ", "|x^2-y^2> "/)
   character(len=50) :: pdos_weight_str
   character(len=100) :: line

   write(*, fmt='(/4x, "U1: ")')
   do j = 1, 3
      write(*, fmt='(/6x,3(f10.5,1X))',  advance='no') U1(j,:,atom_index)
   end do
   write(*, fmt='(/4x, "p-orbital weight decomposition")')
   write(*, '(/6x, A10, 5A10)') "", (trim(p_orb(j)), j = 1, 3)
   do k = 1, 3
      write(*, '(6x, A10, 5F10.5)') trim(p_orb(k)), &
         (U1(k,j,atom_index)*U1(k,j,atom_index), j = 1, 3)
   end do  ! end printing l = 2 orbital weights
   write(*, fmt='(/4x, "U2: ")')
   do j = 1, 5
      write(*, fmt='(/6x,5(f10.5,1X))', advance='no') U2(j,:,atom_index)
   end do
   write(*, fmt='(/4x, "d-orbital weight decomposition")')

   write(*, '(/6x, A10, 5A10)') "", (trim(d_orb(j)), j = 1, 5)
   do k = 1, 5
      write(*, '(6x, A10, 5F10.5)') trim(d_orb(k)), &
         (U2(k,j,atom_index)*U2(k,j,atom_index), j = 1, 5)
   end do  ! end printing l = 2 orbital weights

   ! Orthogonality check
   identity3 = 0.0
   do i = 1, size(identity3(:,1))
      identity3(i,i) = 1.0
   end do
   identity5 = 0.0
   do i = 1, size(identity5(:,1))
      identity5(i,i) = 1.0
   end do
   if(all(abs(matmul(U1(:,:,atom_index), transpose(U1(:,:,atom_index))) - identity3) < 1e-7)) then
      write(*, '(/2x, A)', advance='no') 'U1 matrix is orthogonal: PASS'
   else
      write(*, '(/2x, A)', advance='no') 'U1 matrix is orthogonal: FAIL'
   endif

   if(all(abs(matmul(U2(:,:,atom_index), transpose(U2(:,:,atom_index))) - identity5) < 1e-7)) then
      write(*, '(/2x, A)', advance='no') 'U2 matrix is orthogonal: PASS'
   else
      write(*, '(/2x, A)', advance='no') 'U2 matrix is orthogonal: FAIL'
   endif
   write(*,*)
end subroutine
   ! -----------------------------------------------------------------------------
   ! Subroutine get_pdos_axes
   ! -----------------------------------------------------------------------------

   !!****f* ProcModule/get_pdos_axes *
   !!
   !!  NAME
   !!   get_pdos_axes - Create and normalise local rotation axes
   !!  USAGE
   !!   get_pdos_axes
   !!  PURPOSE
   !!   Create and normalise local rotation axes. Checks for mutual orthogonality.
   !!  INPUTS
   !!    NONE
   !!  USES
   !!   datatypes, local, GenComms
   !!  AUTHOR
   !!   C. Xu
   !!  CREATION DATE
   !!   05/03/2026
   !!  MODIFICATION HISTORY
   !!    03/2026 - C. Xu: Orthogonality check
   !!  SOURCE
   !!
  subroutine get_pdos_axes
   use datatypes
   use local, ONLY: pdos_ax, pdos_ay, pdos_az
   use GenComms, ONLY: cq_abort
   implicit none
   real(double), parameter :: tol = 1e-10

   ! Normalise new x-axis correctly
   pdos_ax = pdos_ax / norm2(pdos_ax)
   ! Normalise new y-axis correctly
   pdos_ay = pdos_ay / norm2(pdos_ay)
   ! Normalise new z-axis correctly
   pdos_az = pdos_az / norm2(pdos_az)
   if (abs(dot_product(pdos_ax, pdos_ay)) > tol) &
      call cq_abort("get_pdos_axes: pDOS_ax vector was not orthogonal to pDOS_ay.")
   if (abs(dot_product(pdos_az, pdos_ay)) > tol) &
       call cq_abort("get_pdos_axes: pDOS_az vector was not orthogonal to pDOS_ay.")
   if (abs(dot_product(pdos_az, pdos_ax)) > tol) &
       call cq_abort("get_pdos_axes: pDOS_az vector was not orthogonal to pDOS_ax.")
  end subroutine

   subroutine initialise_A_mat(A1, A2)
      use datatypes
      implicit none
      real(double), intent(out) :: A1(3, 3), A2(5, 5)
      A1 = 0.0
      A2 = 0.0
      ! A1, FOR p-orbitals, l = 1
      A1(1, 3) = 1.0
      A1(2, 1) = 1.0
      A1(3, 2) = 1.0
      ! d-orbitals, l = 2
      A2(1, 2) = 1.0
      A2(2, 4) = 1.0
      A2(3, 1) = 1.0
      A2(4, 5) = 2.0
      A2(5, 3) = 2.0*sqrt(3.0)
   end subroutine initialise_A_mat

   subroutine antisym_matrix(vector, matrix)
      use datatypes
      implicit none
      real(double), intent(in) :: vector(3)
      real(double), intent(out) :: matrix(3,3)
      matrix = 0.0
      matrix(1, 2) = -vector(3)
      matrix(1, 3) = vector(2)
      matrix(2, 1) = vector(3)
      matrix(2, 3) = -vector(1)
      matrix(3, 1) = -vector(2)
      matrix(3, 2) = vector(1)
   end subroutine
   ! -----------------------------------------------------------------------------
   ! Subroutine euler_from_axisangle
   ! -----------------------------------------------------------------------------

   !!****f* ProcModule/euler_from_axisangle *
   !!
   !!  NAME
   !!   euler_from_axisangle
   !!  USAGE
   !!   euler_from_axisangle(axis, angle)
   !!  PURPOSE
   !!   Calculate Euler angles (alpha, beta, gamma) from a rotation axis and angle.
   !!   The rotation axis is assumed to be in [0, pi]
   !!  INPUTS
   !!    real(double), intent(in) :: axis(3) - axis of rotation
   !!    real(double), intent(in) :: angle   - rotation angle in [0, pi]
   !!  USES
   !!   datatypes, numbers
   !!  AUTHOR
   !!   C. Xu
   !!  CREATION DATE
   !!   20/08/2026
   !!  MODIFICATION HISTORY
   !!
   !!  SOURCE
   !!
   subroutine euler_from_axisangle(axis, angle)
      use datatypes
      use numbers, ONLY: pi
      implicit none
      real(double), intent(in) :: axis(3), angle
      write(*, fmt= &
         '(/2x,"Equivalent Euler angles (up to a sign). Using `Process.RotatePDOSMode 1`should give the same result.")')
      write(*, fmt='(/4x,"alpha: ", (f10.5,1X))', advance="no") &
         180/pi*(datan2(axis(3)*tan(angle / 2),1.0) + datan2(axis(2), axis(1)) - (pi/2))
      write(*, fmt='(/4x,"beta: ", (f10.5,1X))',  advance="no") &
         180/pi*2*asin(sqrt(axis(1)*axis(1) + axis(2)*axis(2)) * sin(angle / 2))
      write(*, fmt='(/4x,"gamma: ", (f10.5,1X))',  advance="no") &
         180/pi*(datan2(axis(3)*tan(angle / 2),1.0) - datan2(axis(2), axis(1)) + (pi/2))
   end subroutine

   ! -----------------------------------------------------------------------------
   ! Subroutine construct_EulerMatrices
   ! -----------------------------------------------------------------------------

   !!****f* ProcModule/construct_EulerMatrices *
   !!
   !!  NAME
   !!   construct_EulerMatrices
   !!  USAGE
   !!   construct_EulerMatrices(E1, E2, atom_index)
   !!  PURPOSE
   !!   Create rotation matrices using Euler angle input: right-handed,
   !!   active Euler angles about intrinsic (fixed) cell axes.
   !!  INPUTS
   !!    integer, intent(in) :: atom_index - atom to perform rotation for
   !!  USES
   !!   datatypes, numbers, local, GenComms
   !!  AUTHOR
   !!   C. Xu
   !!  CREATION DATE
   !!   10/04/2026
   !!  MODIFICATION HISTORY
   !!    13/04/2026 - C. Xu: Add Euler matrix for d-orbitals
   !!    20/07/2026 - C. Xu: Add output matrices in debugging mode
   !!    30/07/2026 - C. Xu: Correct E2(3,2): alpha -> beta
   !!  SOURCE
   !!    Credits to R. Johnson for the mathematical expressions
   !!    Also see: Quantum Theory of Angular Momentum
   subroutine construct_EulerMatrices(E1, E2, atom_index)
      use datatypes
      use numbers, ONLY: pi
      use local, ONLY: euler_angles, flag_rotate_pdos_debug
      use GenComms, ONLY: cq_abort

      implicit none
      integer, intent(in) :: atom_index
      real(double), intent(out) :: E1(3,3), E2(5,5)
      real(double) :: ca, sa, cb, sb, cg, sg, c2a, s2a, c2b, s2b, c2g, s2g
      real(double) :: euler_alpha, euler_beta, euler_gamma
      real(double) :: Q1(3,3), Q2(5,5), identity3(3,3), identity5(5,5)
      integer :: i
      euler_alpha =  euler_angles(1,atom_index)
      euler_beta =  euler_angles(2,atom_index)
      euler_gamma =  euler_angles(3,atom_index)
      !Q matrices are change-of-basis for the orbital convention
      Q1 = 0.0
      Q2 = 0.0

      Q1(1,3) = 1.0
      Q1(2,2) = 1.0
      Q1(3,1) = 1.0

      Q2(1,5) = 1.0
      Q2(2,4) = 1.0
      Q2(3,3) = 1.0
      Q2(4,2) = 1.0
      Q2(5,1) = 1.0

      E1 = 0.0
      E2 = 0.0

      ca = cos(euler_alpha)
      sa = sin(euler_alpha)
      cb = cos(euler_beta)
      sb = sin(euler_beta)
      cg = cos(euler_gamma)
      sg = sin(euler_gamma)

      c2a = cos(2.0 * euler_alpha)
      s2a = sin(2.0 * euler_alpha)
      c2b = cos(2.0 * euler_beta)
      s2b = sin(2.0 * euler_beta)
      c2g = cos(2.0 * euler_gamma)
      s2g = sin(2.0 * euler_gamma)

      ! Rotation of l = 1 coefficients
      E1(1,1) = ca*cb*cg - sa*sg
      E1(1,2) = ca*sb
      E1(1,3) = -cg*sa - ca*cb*sg
      E1(2,1) = -cg*sb
      E1(2,2) = cb
      E1(2,3) = sb*sg
      E1(3,1) = cb*cg*sa + ca*sg
      E1(3,2) = sa*sb
      E1(3,3) = ca*cg - cb*sa*sg


      ! Rotation of l = 2 coefficients
      E2(1,1) = 0.25*c2a*c2g*(3.0+c2b) - s2a*cb*s2g
      E2(1,2) = 0.5*c2a*s2b*cg - s2a*sb*sg
      E2(1,3) = 0.5*sqrt(3.0)*c2a*sb*sb
      E2(1,4) = -sb*(c2a*cb*sg + s2a*cg)
      E2(1,5) = -cb*s2a*c2g -0.25*c2a*(3.0+c2b)*s2g

      E2(2,1) = sa*sb*s2g - 0.5*ca*c2g*s2b
      E2(2,2) = ca*c2b*cg - sa*cb*sg
      E2(2,3) = sqrt(3.0)*ca*sb*cb
      E2(2,4) = -cb*sa*cg -ca*c2b*sg
      E2(2,5) = sb*(ca*cb*s2g + c2g*sa)

      E2(3,1) = 0.5*sqrt(3.0)*sb*sb*c2g
      E2(3,2) = -sqrt(3.0)*sb*cb*cg
      E2(3,3) = 0.25*(3.0*c2b+1.0)
      E2(3,4) = sqrt(3.0)*sb*cb*sg
      E2(3,5) = -sqrt(3.0)*sb*sb*sg*cg

      E2(4,1) = -sb*(sa*cb*c2g + ca*s2g)
      E2(4,2) = sa*c2b*cg + ca*cb*sg
      E2(4,3) = sqrt(3.0)*sa*sb*cb
      E2(4,4) = ca*cb*cg - sa*c2b*sg
      E2(4,5) = sb*(sa*cb*s2g - ca*c2g)

      E2(5,1) = 0.25*s2a*(c2b + 3.0)*c2g + c2a*cb*s2g
      E2(5,2) = sb*(s2a*cb*cg + c2a*sg)
      E2(5,3) = sqrt(3.0)*sa*sb*sb*ca
      E2(5,4) = sb*(c2a*cg - 2*sa*ca*cb*sg)
      E2(5,5) = c2a*cb*c2g - sa*ca*(3.0 + c2b)*sg*cg

      ! Print out matrices before basis conversion
      ! This basis is in descending order of m
      if (flag_rotate_pdos_debug) then
         write(*, fmt='(/4x, "E1 (in old basis): ")')
         do i = 1, 3
            write(*, fmt='(/6x,3(f10.5,1X))', advance='no') E1(i,:)
         end do
         write(*, fmt='(/4x, "E2 (in old basis): ")')
         do i = 1, 5
            write(*, fmt='(/6x,5(f10.5,1X))', advance='no') E2(i,:)
         end do
      end if
      ! Require change of basis to correct orbital convention
      E1 = matmul(inv(Q1), matmul(E1, Q1))
      E2 = matmul(inv(Q2), matmul(E2, Q2))

      if (flag_rotate_pdos_debug) then
         identity3 = 0.0
         identity5 = 0.0
         do i = 1,3
            identity3(i,i) = 1.0
         end do
         do i = 1,5
            identity5(i,i) = 1.0
         end do

         if(all(abs(matmul(E1, transpose(E1)) - identity3) < 1e-7)) then
            write(*, '(/2x, "E1 matrix is orthogonal: PASS")', advance='no')
         else
            write(*, '(/2x, "E1 matrix is orthogonal: FAIL")', advance='no')
         endif
         if(all(abs(matmul(E2, transpose(E2)) - identity5) < 1e-7)) then
            write(*, '(/2x, "E2 matrix is orthogonal: PASS")', advance='no')
         else
            write(*, '(/2x, "E2 matrix is orthogonal: FAIL")', advance='no')
         endif
         write(*,*)
      end if
   end subroutine construct_EulerMatrices
   ! -----------------------------------------------------------------------------
   ! Subroutine calculate_axis_angle
   ! -----------------------------------------------------------------------------

   !!****f* ProcModule/calculate_axis_angle *
   !!
   !!  NAME
   !!   calculate_axis_angle
   !!  USAGE
   !!   call calculate_axis_angle(w1, w2, w3, axis, angle)
   !!  PURPOSE
   !!   Given a set of 3 orthonormal vectors, calculate the axis of rotation and angle needed to go from the simulation
   !!   cell axes to this basis (order matters).
   !!   This subroutine forms the basis change matrix and finds its normalised eigenvector with real eigenvalue one
   !!   This eigenvector is the rotation axis. The rotation angle is computed as the angle satisfying
   !!   cos \theta = (Tr(R) - 1) /2; sin\theta = -Tr(K_n R) / 2, \Theta = datan2(sin, cos) and shifted to the interval
   !!    [0, 2pi]
   !!  INPUTS
   !!    None
   !!  USES
   !!   datatypes, local, global, pao_format, GenComms
   !!  AUTHOR
   !!   C. Xu
   !!  CREATION DATE
   !!   02/03/2026
   !!  MODIFICATION HISTORY
   !!    04/03/2026 - C. Xu: use LAPACK to find rotation axis
   !!    05/03/2026 - C. Xu: Correct use of lapack and rotation angle
   !!    09/03/2026 - C. Xu: use atan2 for rotation angle
   !!    19/03/2026 - C. Xu: initial attempt to fix 90 degree rotation
   !!    24/03/2026 - C. Xu: add error checking
   !!    10/04/2026 - C. Xu: fixes to determinant check, negative angles and euler matrix
   !!    26/05/2026, 08/06 - C. Xu: comments and clearer error messages
   !!    20/08/2026 - C. Xu: Clean debug and negate axis/angle when necessary
   !!  SOURCE
   !!
   subroutine calculate_axis_angle(w1, w2, w3, axis, angle)
      use datatypes
      use numbers, ONLY: pi, one
      use GenComms, ONLY: cq_abort
      use local, ONLY: flag_rotate_pdos_debug
      implicit none
      external DGEEV

      ! Allow the user to define new coordinate system
      ! Rotate from simulation axes (standard Cartesian) to local axes
      real(double), intent(in) ::  w1(3), w2(3), w3(3)
      real(double), intent(out) :: axis(3), angle
      real(double) :: pos_diff(3), det, sin_angle, cos_angle, angle2
      real(double) :: basis_matrix(3,3), temp(3,3), antisym_axis(3,3)
      ! Local variables
      integer :: i, j, N, LDA, LDVL, LDVR, INFO, LDWORK
      real(double), allocatable::  A(:,:), WR(:),   WI(:),   VL(:,:), VR(:,:), WORK(:), tra
      real(double), parameter :: tol = 1e-10
      real(double) :: real_eigenvalue,imag_eigenvalue
      N = 3
      LDA = N
      LDVL = LDA
      LDVR = LDA
      allocate(A(N,N))
      allocate(WR(N))
      allocate(WI(N))
      allocate(VL(LDVL, N))
      allocate(VR(LDVR, N))
      allocate(WORK(5*N))
      LDWORK = 5*N
      ! Define change of basis matrix;
      ! Columns of new basis in terms of old basis
      basis_matrix(:,1) = w1 / norm2(w1)
      basis_matrix(:,2) = w2 / norm2(w2)
      basis_matrix(:,3) = w3 / norm2(w3)
      ! Make copy because DGEEV destroys matrix
      basis_matrix = transpose(basis_matrix)
      A = basis_matrix
      ! If determinant flips sign, orientation of the coordinate system has changed
      det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
      - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
      + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

      if (det < 1.0 - tol .or. det > 1.0 + tol) &
         call cq_abort("calculate_axis_angle: determinant is not +1. Make sure input axes form a right-handed basis.")
      call DGEEV("N", "V", N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, &
           WORK, LDWORK, INFO)
      if (INFO /= 0) then
         write(*, fmt='(A,I0)') 'DGEEV INFO = ', INFO
         call cq_abort('calculate_axis_angle: DGEEV, Eigenvalues of basis change matrix could not be found!')
      end if
      ! Require eigenvector with eigenvalue 1, no complex, to find axis
      do i = 1, N
         if (abs(WR(i) - one) < tol .and. abs(WI(i)) < tol) then
            axis  = VR(:, i)   ! Column i of VR is the i-th right eigenvector
            exit
         end if
      end do
      if (flag_rotate_pdos_debug) then
         real_eigenvalue = 0.0
         imag_eigenvalue = 0.0
         write(*, '(/2x, "Eigenvalues/Eigenvectors of change of base matrix: ")')
         do j = 1, N
            ! Print eigenvalue
            if (abs(WI(j)) < 1e-7) then
               write(*,'(/4x, A,I3,A,F12.6)') 'Eigenvalue ', j, ': ', WR(j)
            else
               write(*,'(/4x, A,I3,A,F12.6,SP,F12.6,A)') 'Eigenvalue ', j, ': ', WR(j), WI(j), 'i'
            end if
            write(*,'(/6x, 2X,A)', advance="no") 'Eigenvector:'
            do i = 1, N
               if (WI(j) == 0.0d0) then
                     write(*,'(/7x, 4X,A,I3,A,F10.5)', advance="no") 'component ', i, ': ', VR(i,j)
               else if (WI(j) > 0.0d0) then
                     ! First complex conjugate pair: v = VR(:,j) + i*VR(:,j+1)
                     write(*,'(/7x, 4X,A,I3,A,F10.5,SP,F10.5,A)', advance="no") &
                        'component ', i, ': ', VR(i,j), VR(i,j+1), 'i'
               else
                     ! Second complex conjugate pair: v = VR(:,j-1) - i*VR(:,j)
                     write(*,'(/7x, 4XA,I3,A,F10.5,SP,F10.5,A)', advance="no") &
                        'component ', i, ': ', VR(i,j-1), -VR(i,j), 'i'
               end if
            end do
            write(*,*)
         end do
      end if ! End debug: print all eigenvalues/vectors
      if (any(axis /= axis)) &
         call cq_abort("calculate_axis_angle: NaN in rotation axis.")
      if (all(abs(axis) < tol)) &
         call cq_abort("calculate_axis_angle: Rotation axis was (0,0,0). Cannot perform rotation.")
      tra = basis_matrix(1,1) + basis_matrix(2,2) + basis_matrix(3,3)
      antisym_axis = 0.0
      axis = axis / norm2(axis)
      cos_angle = (tra - 1.0)
      call antisym_matrix(axis, antisym_axis)
      temp = matmul(antisym_axis, basis_matrix)
      sin_angle = -(temp(1,1) + temp(2,2) + temp(3,3))
      if (abs(sin_angle) < tol .and. abs(cos_angle) < tol) &
            call cq_abort("calculate_axis_angle: Both arguments to datan2 are zero.")
      angle = datan2(sin_angle, cos_angle)
      if (angle /= angle) &
            call cq_abort("calculate_axis_angle: NaN rotation angle.")
      if (flag_rotate_pdos_debug) &
            write(*, fmt='(/4x,"Rotation angle (rad/deg) [-pi, pi]: ",2(f10.5,1X))') &
                  angle, angle * 180/pi
            write(*, fmt='(/4x,"Rotation axis from LAPACK : ",3(f10.5,1X))') &
               axis
      if (angle < 0.0) then
         ! Keep angle between [0, pi] and negate both axis and angle
         angle = -angle
         axis = -axis !R(n, theta) = R(-n, -theta)
      end if
      write(*, fmt='(/4x,"Rotation angle (rad/deg) [0, pi]: ",2(f10.5,1X))') &
               angle, angle * 180/pi
      write(*, fmt='(/4x,"Rotation axis: ", 3(f10.5,1X))') axis
   end subroutine calculate_axis_angle
   subroutine construct_rodrigues(axis, angle, matrix)
      use datatypes
      implicit none
! We will compute the rodrigues matrix, specifically R^T(axis, -angle)
      real(double), intent(in) :: axis(3), angle
      real(double), intent(out) :: matrix(9)
      real(double) :: K(3, 3), KT(3, 3), identity(3, 3), temp(3, 3)
      real(double) :: norm_axis(3)
      K = 0.0
      identity = 0.0
      identity(1, 1) = 1.0
      identity(2, 2) = 1.0
      identity(3, 3) = 1.0
      norm_axis = axis / norm2(axis)
      ! Define K
      call antisym_matrix(norm_axis, K)

      KT = transpose(K)
      ! Rodrigues rotation matrix but with -angle
      temp = identity + (sin(angle)*KT) + (1 - cos(angle))*matmul(KT, KT)
      ! Check the transpose with a simple utility
      matrix = reshape(transpose(temp), shape=(/9/))
   end subroutine construct_rodrigues

   subroutine construct_C1(rod, C1)
      use datatypes
      use local, ONLY: flag_rotate_pdos_debug
      implicit none
      real(double), intent(in) :: rod(9)
      real(double), intent(out) :: C1(3, 3)
      real(double) :: identity(3,3)
      integer :: i
      C1(1, 1) = rod(1)
      C1(1, 2) = rod(2)
      C1(1, 3) = rod(3)
      C1(2, 1) = rod(4)
      C1(2, 2) = rod(5)
      C1(2, 3) = rod(6)
      C1(3, 1) = rod(7)
      C1(3, 2) = rod(8)
      C1(3, 3) = rod(9)
      if (flag_rotate_pdos_debug) then
         identity = 0.0
         identity(1, 1) = 1.0
         identity(2, 2) = 1.0
         identity(3, 3) = 1.0
         if(all(abs(matmul(C1, transpose(C1)) - identity) < 1e-7)) then
            write(*, '(/2x, "C1 matrix is orthogonal: PASS")', advance='no')
         else
            write(*, '(/2x, "C1 matrix is orthogonal: FAIL")', advance='no')
         endif
      end if
   end subroutine construct_C1
   subroutine construct_C2(r, C2)
      use datatypes
      use local, ONLY: flag_rotate_pdos_debug
      implicit none
      real(double), intent(in) :: r(9)
      real(double), intent(out) :: C2(5, 5)
      real(double) :: identity(5,5)

      C2(1, 1) = r(6)*r(8) + r(5)*r(9)
      C2(1, 2) = r(6)*r(7) + r(4)*r(9)
      C2(1, 3) = r(5)*r(7) + r(4)*r(8)
      C2(1, 5) = r(6)*r(9)*0.5
      C2(1, 4) = r(4)*r(7) + r(6)*r(9)*0.5
! Row 2
      C2(2, 1) = r(3)*r(8) + r(2)*r(9)
      C2(2, 2) = r(3)*r(7) + r(1)*r(9)
      C2(2, 3) = r(2)*r(7) + r(1)*r(8)
      C2(2, 5) = r(3)*r(9)*0.5
      C2(2, 4) = r(1)*r(7) + r(3)*r(9)*0.5
! Row 3
      C2(3, 1) = r(3)*r(5) + r(2)*r(6)
      C2(3, 2) = r(3)*r(4) + r(1)*r(6)
      C2(3, 3) = r(2)*r(4) + r(1)*r(5)
      C2(3, 5) = r(3)*r(6)*0.5
      C2(3, 4) = r(1)*r(4) + r(3)*r(6)*0.5

! Row 4
      C2(4, 1) = 2*(r(2)*r(3) - r(5)*r(6))
      C2(4, 2) = 2*(r(1)*r(3) - r(4)*r(6))
      C2(4, 3) = 2*(r(1)*r(2) - r(4)*r(5))
      C2(4, 5) = 0.5*(r(3)*r(3) - r(6)*r(6))
      C2(4, 4) = r(1)*r(1) - r(4)*r(4) + 0.5*(r(3)*r(3) - r(6)*r(6))

! Row 5
      C2(5, 1) = (4*r(8)*r(9)) - 2*(r(2)*r(3) + r(5)*r(6))
      C2(5, 2) = (4*r(7)*r(9)) - 2*(r(1)*r(3) + r(4)*r(6))
      C2(5, 3) = (4*r(7)*r(8)) - 2*(r(1)*r(2) + r(4)*r(5))
      C2(5, 5) = r(9)*r(9) - 0.5*(r(3)*r(3) + r(6)*r(6))
      C2(5, 4) = C2(5, 5) - r(1)*r(1) - r(4)*r(4) + 2*r(7)*r(7)
      if (flag_rotate_pdos_debug) then
         identity = 0.0
         identity(1, 1) = 1.0
         identity(2, 2) = 1.0
         identity(3, 3) = 1.0
         identity(4, 4) = 1.0
         identity(5, 5) = 1.0
         if(all(abs(matmul(C2, transpose(C2)) - identity) < 1e-7))  then
            write(*, '(/2x, "C2 matrix is orthogonal: PASS")', advance='yes')
         else
            write(*, '(/2x, "C2 matrix is orthogonal: FAIL")', advance='yes')
         endif
      end if
   end subroutine construct_C2
   function inv(A) result(Ainv)
      use datatypes
      use GenComms, ONLY: cq_abort
      implicit none

      real(double), dimension(:, :), intent(in) :: A
      real(double), dimension(size(A, 1), size(A, 2)) :: Ainv

      real(double), dimension(size(A, 1)) :: work  ! work array for LAPACK
      integer, dimension(size(A, 1)) :: ipiv   ! pivot indices
      integer :: n, info

      external DGETRF
      external DGETRI

      Ainv = A
      n = size(A, 1)

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      call DGETRF(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
         call cq_abort("inv: Matrix is numerically singular!")
      end if

      ! DGETRI computes the inverse using the LU factorization by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
         call cq_abort("inv: Matrix inversion failed!")
      end if
   end function inv

   subroutine construct_Ul(angmom, Al, Cl, Ul)
      use datatypes
      implicit none

      integer, intent(in) :: angmom
      real(double), intent(in) ::   Al(2*angmom + 1, 2*angmom + 1), Cl(2*angmom + 1, 2*angmom + 1)
      real(double), intent(out) :: Ul(:,:)
      real(double) :: Al_inv(2*angmom + 1, 2*angmom + 1)
      Al_inv = inv(Al)
      Ul = matmul(Al_inv, matmul(Cl, Al))

   end subroutine construct_Ul
   ! -----------------------------------------------------------------------------
   ! Subroutine rotate_coefficients
   ! -----------------------------------------------------------------------------

   !!****f* ProcModule/rotate_coefficients *
   !!
   !!  NAME
   !!   rotate_coefficients
   !!  USAGE
   !!   call rotate_coefficients
   !!  PURPOSE
   !!   Apply rotation matrices to wavefunction coefficients. Formally, for each l, each 2l + 1 coefficients
   !!   corresponding to a (m, zeta) are multiplied by the corresponding U^{(l)} rotation matrix.
   !!   As different rotation modes use different arrays and lookups, these are taken care of at the start of the
   !!   subroutine
   !!   The loop over all values is pasted and modified from process_pdos.
   !!  INPUTS
   !!    None
   !!  USES
   !!   datatypes, local, global, pao_format, GenComms
   !!  AUTHOR
   !!   C. Xu
   !!  CREATION DATE
   !!   02/03/2026
   !!  MODIFICATION HISTORY
   !!    14/04/2026 - C. Xu: loop over active bands only
   !!    03/05/2026 - C. Xu: correct g_atom_lookup
   !!    12/05/2026 - C. Xu: correct de/allocation behaviour. Edit desc
   !!  SOURCE
   !!
   subroutine rotate_coefficients
      use datatypes
      use local, ONLY: n_bands_total, nkp, n_atoms_pDOS, evec_coeff, scaled_evec_coeff, &
      pDOS_atom_index, band_full_to_active, rotate_pdos_natoms, find_neighbours, &
      flag_rotate_pdos_mode, rotate_pdos_atoms_euler, U1, U2
      use global_module, ONLY: nspin, species_glob
      use pao_format,    ONLY: pao
      use GenComms, ONLY: cq_abort

      implicit none
      ! Local variables
      integer :: i_atom, i_spec, i_band, i_band_c, i_kp, i_spin, g_atom
      integer :: i_l, i_z, nzeta, norbs, sf_offset, rotate_counter
      integer, dimension(:), allocatable :: g_atom_lookup
      ! Currently input axes depends on n_atoms_pDOS
      if (flag_rotate_pdos_mode == 0) then
         rotate_counter = n_atoms_pDOS
      else
         rotate_counter = rotate_pdos_natoms
      end if
      if (.not. allocated(g_atom_lookup)) allocate(g_atom_lookup(rotate_counter))
      !  Define the lookups correctly for each mode
      if (flag_rotate_pdos_mode == 2) then
         ! Atom counter should be the order the user input pDOSNeighbours block
            g_atom_lookup = find_neighbours(1, :)
      else if (flag_rotate_pdos_mode == 1) then
         ! Atom counter should be the order the user input pDOSEuler block
            g_atom_lookup = rotate_pdos_atoms_euler
      else
            g_atom_lookup = pDOS_atom_index
      end if
      do i_spin = 1, nspin
         do i_kp = 1, nkp
            do i_band = 1, n_bands_total
               i_band_c = band_full_to_active(i_band)
               if(i_band_c > 0) then
               do i_atom = 1, rotate_counter
                  ! Get global atom number from input
                  g_atom = g_atom_lookup(i_atom)
                  i_spec = species_glob(g_atom)
                  sf_offset = 0
                  ! Include l = 0 to correctly calculate offset
                  do i_l = 0, pao(i_spec)%greatest_angmom
                     nzeta = pao(i_spec)%angmom(i_l)%n_zeta_in_angmom
                     norbs = 2*i_l + 1
                     do i_z = 1, nzeta
                        select case(i_l)
                        case(1)
                           evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin) = &
                              matmul(U1(:,:,i_atom), evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin))
                           scaled_evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin) = &
                              matmul(U1(:,:,i_atom), scaled_evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin))
                        case(2)
                           evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin) = &
                              matmul(U2(:,:,i_atom), evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin))
                           scaled_evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin) = &
                              matmul(U2(:,:,i_atom), scaled_evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin))
                        end select
                        sf_offset = sf_offset + norbs
                     end do ! i_z
                  end do ! i_l
               end do ! i_atom
            end if ! if band is active
            end do ! i_band
         end do ! i_kp
      end do ! i_spin
      deallocate(g_atom_lookup)
   end subroutine rotate_coefficients

   ! -----------------------------------------------------------------------------
   ! Subroutine nearest_neighbours
   ! -----------------------------------------------------------------------------

   !!****f* ProcModule/nearest_neighbours *
   !!
   !!  NAME
   !!   nearest_neighbours - Find nearest neighbours and their bond vectors
   !!  USAGE
   !!   nearest_neighbours(atomno, bond)
   !!  PURPOSE
   !!   Evaluates the nearest neighbours using periodic boundary conditions
   !!    and returns their bond vectors
   !!
   !!  INPUTS
   !!   integer, intent(in) :: atomno !atom to find neighbours of
   !!   real(double), intent(out), allocatable :: bond(:,:) ! array to hold bond
   !!       information regarding an atom
   !!  USES
   !!   datatypes, dimens, local, global, GenComms
   !!  AUTHOR
   !!   C. Xu
   !!  CREATION DATE
   !!   15/04/2026
   !!  MODIFICATION HISTORY
   !!    12/05/2026 - C. Xu: correct de/allocation behaviour. Edit desc
   !!  SOURCE
   !!
   subroutine nearest_neighbours(atomno, bond)
      use datatypes
      use dimens, ONLY: r_super_x, r_super_y, r_super_z
      use local, ONLY: find_neighbours, nghbr_arr, flag_rotate_pdos_debug
      use global_module, ONLY: ni_in_cell, atom_coord, species_glob
      use GenComms, ONLY: cq_abort
      use species_module, ONLY: species_label

      implicit none

      integer, intent(in) :: atomno
      real(double), intent(out), allocatable :: bond(:,:)




      ! Local variables
      real(double), parameter :: err = 1e-5
      real(double) :: distances(1:ni_in_cell), temp(1:ni_in_cell)
      real(double) ::cx, cy, cz, dx, dy, dz, atomno_pos(3)
      integer :: i, j, min_idx, ibond_max_len, idx_direction
      ! Reset allocatables as this subroutine is called in a loop
      ! and populated with new data every time
      ! and populated with new data every time
      if (allocated(nghbr_arr)) deallocate(nghbr_arr)
      if (allocated(bond)) deallocate(bond)

      cx = atom_coord(1, atomno)
      cy = atom_coord(2, atomno)
      cz = atom_coord(3, atomno)
      ! ACCOUNT FOR PERIODICITY
      do i = 1, ni_in_cell

         if (i == atomno) then
            distances(i) = huge(1.0d0)
            cycle
         end if

         dx = atom_coord(1,i) - cx
         dy = atom_coord(2,i) - cy
         dz = atom_coord(3,i) - cz

         dx = dx - r_super_x * nint(dx / r_super_x)
         dy = dy - r_super_y * nint(dy / r_super_y)
         dz = dz - r_super_z * nint(dz / r_super_z)

         distances(i) =  dx*dx + dy*dy + dz*dz
      end do

      temp = distances
      ! Eliminate zeros - self-distance
      where (abs(temp) < err) temp = huge(1.0)
      select case (find_neighbours(2, findloc(find_neighbours(1,:), atomno, dim=1)))
      case (0)
         allocate(nghbr_arr(4))
      case (1)
         allocate(nghbr_arr(6))
      case default
         call cq_abort("nearest_neighbours: Did not correctly allocate nghbr_arr.")
      end select
      do j = 1, size(nghbr_arr)
         min_idx = minloc(temp, 1)
         nghbr_arr(j) = min_idx
         temp(min_idx) = huge(1.0)
      end do

      allocate(bond(3, size(nghbr_arr)))

      atomno_pos = (/cx, cy, cz/)
      do i = 1, size(nghbr_arr)
         bond(1, i) = atom_coord(1, nghbr_arr(i)) - atomno_pos(1)
         bond(2, i) = atom_coord(2, nghbr_arr(i)) - atomno_pos(2)
         bond(3, i) = atom_coord(3, nghbr_arr(i)) - atomno_pos(3)

         bond(1, i) = bond(1, i) - r_super_x * nint(bond(1, i) / r_super_x)
         bond(2, i) = bond(2, i) - r_super_y * nint(bond(2, i) / r_super_y)
         bond(3, i) = bond(3, i) - r_super_z * nint(bond(3, i) / r_super_z)

      end do

       if (flag_rotate_pdos_debug) then
         write(*, fmt='(/2x,"Unit cell dimensions (x,y,z)", 3(f10.5))', advance="no") &
            r_super_x, r_super_y, r_super_z
         write(*, fmt='(/2x,"Outputting cell-periodic distances ", &
            "of all atoms relative to atom", I0)', advance="no") atomno
         do i = 1, ni_in_cell
            if (i /= atomno) then
               write(*, fmt='(/4x, I0, A, f10.5)', advance="no") &
                  i, " " // species_label(species_glob(i)), distances(i)
            end if
         end do
         write(*, fmt='(/2x,"Outputting cell-periodic unnormalised bond lengths")')

         do i = 1, size(nghbr_arr)
            write(*, fmt='(/4x, I0, A, 3f10.5)', advance="no") &
               nghbr_arr(i), " " // species_label(species_glob(nghbr_arr(i))), bond(:,i)
         end do

      end if

   end subroutine
   ! -----------------------------------------------------------------------------
   ! Subroutine axes_from_nn
   ! -----------------------------------------------------------------------------

   !!****f* ProcModule/axes_from_nn *
   !!
   !!  NAME
   !!   axes_from_nn - Gets rotation axes from bond vectors
   !!  USAGE
   !!   axes_from_nn(atomno, bond)
   !!  PURPOSE
   !! This subroutine will use computed nearest_neighbours
   !! This subroutine will use computed nearest_neighbours
   !!  and specified local geometry to construct  a new set of local
   !!  axes to rotate the pDOS into
   !!
   !!  For square-planar geometry:
   !!
   !!  For square-planar geometry:
   !!   Longest bond is chosen as x_hat.
   !!   Bonds which are closest to orthogonality are
   !!      projected onto a plane defined by x_hat
   !!   y_hat is chosen from the projected direction
   !!      which has the minimal difference to its unprojected bond
   !!   z_hat is computed as the cross-product from these 2 directions, and thus
   !!   z_hat is computed as the cross-product from these 2 directions, and thus
   !!      point perpendicular to the planar geometry
   !!  For octahedra: choose z_hat as longest or shortest bond
   !!   Use this direction to define a plane - project the 4 atoms onto it
   !!   Call y_hat the bond which is the closest to the plane
   !!    Find x_hat = z cross y
   !!  INPUTS
   !!   integer, intent(in) :: atomno - atom to find neighbours of
   !!   real(double), intent(in) :: bond(:,:) - array to hold bond vectors
   !!  USES
   !!   datatypes, local, global, GenComms
   !!  AUTHOR
   !!   C. Xu
   !!  CREATION DATE
   !!   15/04/2026
   !!  MODIFICATION HISTORY
   !!   20/04/2026 C. Xu - allow user to specify neighbour as secondary direction
   !!  SOURCE
   !!
   subroutine axes_from_nn(atomno, bond)
      use datatypes
      use local, ONLY: find_neighbours, nghbr_arr, pdos_ax, pdos_ay, pdos_az
      use global_module, ONLY: ni_in_cell, atom_coord
      use GenComms, ONLY: cq_abort
      implicit none
      ! This subroutine will use computed nearest_neighbours
      ! and specified local geometry to construct  a new set of local
      ! axes to rotate the pDOS into
      !
      ! For square-planar geometry:
      !  Longest bond is chosen as x_hat.
      !  Bonds which are closest to orthogonality are
      !     projected onto a plane defined by x_hat
      !  y_hat is chosen from the projected direction
      !     which has the minimal difference to its unprojected bond
      !  z_hat is computed as the cross-product from these 2 directions, and thus
      !     point perpendicular to the planar geometry
      ! For octahedra: choose z_hat as longest or shortest bond
      !  Use this direction to define a plane - project the 4 atoms onto it
      !  Call y_hat the bond which is the closest to the plane
      ! Find x_hat by x_hat = z cross y
      integer, intent(in) :: atomno
      real(double), intent(in) :: bond(:,:)
      ! Local variables
      real(double) ::atomno_pos(3), plane_normal(3), proj_vector(3), cos_angle
      real(double), allocatable :: dots(:), bond_lengths(:), norm_bond(:,:)
      integer :: i, ibond_principal, ibond_sec, idx_direction, minormax, find_atomno

      allocate(dots(size(bond(1,:))))
      allocate(norm_bond(3, size(bond(1,:))))
      allocate(bond_lengths(size(bond(1,:))))
      do i = 1, size(bond_lengths)
         bond_lengths(i) = sqrt(dot_product(bond(:, i), bond(:, i)))
         if (bond_lengths(i) > 1e-10) then
              norm_bond(:, i) = bond(:, i) / bond_lengths(i)
          else
              call cq_abort("axes_from_nn: Zero-length bond")
          end if
      end do
      ! find_atomno: index of user supplied central atom, atomno
      find_atomno = findloc(find_neighbours(1,:), atomno, dim=1)
      minormax = find_neighbours(3, find_atomno)
      if (minormax < 0) then
        ibond_principal = minloc(bond_lengths,1) ! gets neighbour with min bond length
      else if (minormax .eq. 0) then
         ibond_principal = maxloc(bond_lengths,1) ! gets neighbour with maximum bond length
      else
         ibond_principal = findloc(nghbr_arr, minormax, dim=1)
      end if
      if (ibond_principal == 0) &
      call cq_abort('axes_from_nn: principal neighbour not found for atom ', find_atomno)

      write(*,fmt='(/2x,"Principal neighbour found: atom ", (I0,1x))') nghbr_arr(ibond_principal)
      dots = matmul(transpose(norm_bond), norm_bond(:, ibond_principal))  ! dot prod with chosen normal
      ! select index of closest to perpendicular bond - there may be two
      idx_direction = minloc(abs(dots), 1)
      if (find_neighbours(4, find_atomno) == 0) then
         ! If 0, we choose second direction by closest projection
         call project_onto_plane(bond(:, ibond_principal),bond(:, idx_direction), proj_vector)
      else
         ! This entry > 0 -> corresponds to neighbour
         ibond_sec = findloc(nghbr_arr, find_neighbours(4, find_atomno), dim=1)
         ! Handle when input neighbour is not found
         if (ibond_sec == 0) &
            call cq_abort('axes_from_nn: second neighbour not found for atom ', find_atomno)

         if (abs(dots(ibond_sec)) > 0.5) &
            print *, "WARNING: Chosen secondary neighbour appears to not be very perpendicular to principal."
         write(*,fmt='(/2x,"Secondary neighbour found: atom ", (I0,1x))') nghbr_arr(ibond_sec)

         call project_onto_plane(bond(:, ibond_principal),bond(:, ibond_sec), proj_vector)
   end if ! choice of 2nd direction
   if ( find_neighbours(1, find_atomno) == find_neighbours(4, find_atomno)) &
            call cq_abort("axes_from_nn: cannot have the principal neighbour also as the second direction")
      select case (find_neighbours(2, find_atomno))
      case (0)
         pdos_ax = norm_bond(:, ibond_principal)
         pdos_ay = proj_vector
         pdos_az = cross_product(pdos_ax, pdos_ay)
      case (1)
         pdos_ay = proj_vector
         pdos_az = norm_bond(:, ibond_principal)
         pdos_ax = cross_product(pdos_ay, pdos_az)

      end if
      pdos_ax = pdos_ax / norm2(pdos_ax)
      pdos_ay = pdos_ay / norm2(pdos_ay)
      pdos_az = pdos_az / norm2(pdos_az)
   end subroutine axes_from_nn
   function cross_product(a, b) result(cross)
      use datatypes
      implicit none
      real(double), dimension(3), intent(in) :: a, b
      real(double), dimension(3) :: cross

      cross(1) = a(2) * b(3) - a(3) * b(2)
      cross(2) = a(3) * b(1) - a(1) * b(3)
      cross(3) = a(1) * b(2) - a(2) * b(1)
   end function cross_product

   subroutine project_onto_plane(plane_normal, vector, proj_vector)
      use datatypes
      implicit none
      real(double), intent(in) :: plane_normal(3), vector(3)
      real(double), intent(out) :: proj_vector(3)
      real(double) :: norm_plane_normal(3)
      ! Subroutine to project vector onto plane defined by plane_normal vector
      norm_plane_normal = plane_normal / norm2(plane_normal)

      proj_vector = vector - ((dot_product(vector, norm_plane_normal))*norm_plane_normal)
   end subroutine project_onto_plane
  subroutine process_band_structure

    use datatypes
    use numbers, ONLY: zero, RD_ERR, twopi, half, one, two, four, six
    use local, ONLY: eigenvalues, n_bands_total, nkp, wtk, efermi, flag_total_iDOS, &
         flag_procwf_range_Ef, kx,  ky, kz, flag_proc_band_str
    use read, ONLY: read_eigenvalues, read_psi_coeffs
    use global_module, ONLY: nspin, n_DOS, E_DOS_min, E_DOS_max, sigma_DOS
    use units, ONLY: HaToeV

    implicit none

    ! Local variables
    integer :: i_band, i_kp, i_spin, n_DOS_wid, n_band, n_min, n_max, i
    real(double) :: Ebin, dE_DOS, a, pf_DOS, spin_fac, dE
    real(double), dimension(nspin) :: total_electrons, iDOS_low
    real(double), dimension(:,:), allocatable :: total_DOS, iDOS
    real(double), dimension(:,:), allocatable :: occ

    write(*,fmt='(/2x,"Processing eigenvalues to write band structure")')
    if(nspin==1) then
       spin_fac = two
    else if(nspin==2) then
       spin_fac = one
    end if
    ! Read eigenvalues
    call read_eigenvalues
    if(abs(E_DOS_min)<RD_ERR) then
       E_DOS_min = minval(eigenvalues(1,:,:))
       write(*,fmt='(2x,"Band structure lower limit set automatically: ",f12.5," Ha")') E_DOS_min
    else
       write(*,fmt='(2x,"Band structure lower limit set by user: ",f12.5," Ha")') E_DOS_min
    end if
    if(abs(E_DOS_max)<RD_ERR) then
       E_DOS_max = maxval(eigenvalues(n_bands_total,:,:))
       write(*,fmt='(2x,"Band structure upper limit set automatically: ",f12.5," Ha")') E_DOS_max
    else
       write(*,fmt='(2x,"Band structure upper limit set by user: ",f12.5," Ha")') E_DOS_max
    end if
    write(*,fmt='(2x,"Writing band structure files")')
    if(flag_proc_band_str==4) then
       write(*,fmt='(4x,"X-axis will be k-point index")')
    else if(flag_proc_band_str==0) then
       write(*,fmt='(4x,"All k-point coordinates provided")')
    else if(flag_proc_band_str==1) then
       write(*,fmt='(4x,"X-axis will be kx")')
    else if(flag_proc_band_str==2) then
       write(*,fmt='(4x,"X-axis will be ky")')
    else if(flag_proc_band_str==3) then
       write(*,fmt='(4x,"X-axis will be kz")')
    end if
    open(unit=17, file="BandStructure.dat")
    do i_spin = 1, nspin
       dE = zero
       if(flag_procwf_range_Ef) dE = -efermi(i_spin)
       write(17,fmt='("# Spin ",I1)') i_spin
       if(flag_procwf_range_Ef) then
          write(17,fmt='("# Original Fermi-level: ",f12.5," eV")') HaToeV*efermi(i_spin)
          write(17,fmt='("# Bands shifted relative to Fermi-level")')
       else
          write(17,fmt='("# Fermi-level: ",f12.5," eV")') HaToeV*efermi(i_spin)
       end if
       do i_band=1,n_bands_total ! All bands
          if(minval(eigenvalues(i_band, :, i_spin))+dE>=E_DOS_min .and. &
               maxval(eigenvalues(i_band, :, i_spin))+dE<=E_DOS_max) then
             write(17,fmt='("# Band ",i6)') i_band
             do i_kp = 1, nkp
                if(flag_procwf_range_Ef) then
                   if(flag_proc_band_str==4) then
                      write(17,fmt='(i4,f20.12)') i_kp, &
                           HaToeV*(eigenvalues(i_band, i_kp, i_spin) - efermi(i_spin))
                   else if(flag_proc_band_str==0) then
                      write(17,fmt='(4f20.12)') kx(i_kp), ky(i_kp), kz(i_kp), &
                           HaToeV*(eigenvalues(i_band, i_kp, i_spin) - efermi(i_spin))
                   else if(flag_proc_band_str==1) then
                      write(17,fmt='(2f20.12)') kx(i_kp), &
                           HaToeV*(eigenvalues(i_band, i_kp, i_spin) - efermi(i_spin))
                   else if(flag_proc_band_str==2) then
                      write(17,fmt='(2f20.12)') ky(i_kp), &
                           HaToeV*(eigenvalues(i_band, i_kp, i_spin) - efermi(i_spin))
                   else if(flag_proc_band_str==3) then
                      write(17,fmt='(2f20.12)') kz(i_kp), &
                           HaToeV*(eigenvalues(i_band, i_kp, i_spin) - efermi(i_spin))
                   end if
                else
                   if(flag_proc_band_str==4) then
                      write(17,fmt='(i4,f20.12)') i_kp, HaToeV*eigenvalues(i_band, i_kp, i_spin)
                   else if(flag_proc_band_str==0) then
                      write(17,fmt='(4f20.12)') kx(i_kp), ky(i_kp), kz(i_kp), &
                           HaToeV*eigenvalues(i_band, i_kp, i_spin)
                   else if(flag_proc_band_str==1) then
                      write(17,fmt='(2f20.12)') kx(i_kp), HaToeV*eigenvalues(i_band, i_kp, i_spin)
                   else if(flag_proc_band_str==2) then
                      write(17,fmt='(2f20.12)') ky(i_kp), HaToeV*eigenvalues(i_band, i_kp, i_spin)
                   else if(flag_proc_band_str==3) then
                      write(17,fmt='(2f20.12)') kz(i_kp), HaToeV*eigenvalues(i_band, i_kp, i_spin)
                   end if
                end if
             end do ! i_kp = nkp
             write(17,fmt='("&")')
          end if
       end do
    end do
    close(unit=17)
    return
  end subroutine process_band_structure

  ! Important note
  !
  ! Formally we have: PAO(\mathbf{r}) = f(r) r^l Y_{lm}(\hat{\mathbf{r}})
  !
  ! However it is much easier when dealing with the explicit Cartesian form
  ! of spherical harmonics to scale the Y_{lm} by r^l (because the Cartesian
  ! form has 1/r^l as part of it) so this is what we do.  The routine
  ! spherical_harmonic_rl returns this (and its differential)
  !
  ! So we have dPAO/dalpha = df/dr.dr/dalpha.(r^lY) + f.d(r^lY)/dalpha where alpha
  ! is x/y/z and dr/dalpha = alpha/r.
  subroutine pao_dpao_to_grid(i_band, i_kp, i_spin, psi, dpsi)

    use datatypes
    use numbers
    use units
    use global_module, ONLY: ni_in_cell, atom_coord, species_glob
    use pao_format, ONLY: pao
    use local, ONLY: nptsx, nptsy, nptsz, grid_x, grid_y, grid_z, stm_x_min, stm_x_max, &
         stm_y_min, stm_y_max, stm_z_min, stm_z_max, evec_coeff, kx, ky, kz, i_job
    use dimens, ONLY: RadiusAtomf, r_super_x, r_super_y, r_super_z
    use angular_coeff_routines, ONLY: evaluate_pao, pao_elem_derivative_2

    implicit none

    ! Passed variables
    integer :: i_band, i_kp, i_spin
    complex(double_cplx), dimension(nptsx, nptsy, nptsz) :: psi
    complex(double_cplx), dimension(nptsx, nptsy, nptsz, 3) :: dpsi

    ! Local variables
    integer :: i_atom, i_grid_x, i_grid_y, i_grid_z, i_l, i_zeta, i_m, ix, iy, iz
    integer :: i_spec, j, npao, i_mult
    integer :: minx, maxx, miny, maxy, minz, maxz
    real(double) :: dr, dx, dy, dz, sph_rl, f_r, df_r, dx_dr, dy_dr, dz_dr, del_r
    real(double) :: a, b, c, d, r1, r2, r3, r4, rr, kr, krx, kry, krz
    real(double), dimension(3) :: dsph_rl, dg
    complex(double_cplx) :: phase, phase_shift

    psi = zero
    dpsi = zero
    ! Grid spacing
    dg(1) = grid_x!/BohrToAng
    dg(2) = grid_y!/BohrToAng
    dg(3) = grid_z!/BohrToAng
    ! Loop over atoms
    do i_atom = 1, ni_in_cell
       i_spec = species_glob(i_atom)
       if(atom_coord(3, i_atom) + RadiusAtomf(i_spec) >= stm_z_min) then ! Is the atom in STM region?
          kr = kx(i_kp)*atom_coord(1, i_atom) + ky(i_kp)*atom_coord(2, i_atom) + kz(i_kp)*atom_coord(3, i_atom)
          ! Find grid limits
          minx = floor( (atom_coord(1, i_atom) - RadiusAtomf(i_spec))/dg(1) )
          maxx = floor( (atom_coord(1, i_atom) + RadiusAtomf(i_spec))/dg(1) ) + 1
          miny = floor( (atom_coord(2, i_atom) - RadiusAtomf(i_spec))/dg(2) )
          maxy = floor( (atom_coord(2, i_atom) + RadiusAtomf(i_spec))/dg(2) ) + 1
          minz = floor( (atom_coord(3, i_atom) - RadiusAtomf(i_spec))/dg(3) )
          maxz = floor( (atom_coord(3, i_atom) + RadiusAtomf(i_spec))/dg(3) ) + 1
          if(i_job==4.or.i_job==5) then ! STM not band density, so no z periodicity
             if(stm_z_min>zero) then
                minz = minz - floor(stm_z_min/dg(3))
             end if
             if(minz<1) minz = 1
             maxz = min(maxz, nptsz)
          end if
          ! Loop over grid points
          do i_grid_z = minz, maxz
             ! Find z grid position and dz
             dz = dg(3)*real(i_grid_z,double)+stm_z_min - atom_coord(3,i_atom)
             ! Wrap if we're making band densities, but not for STM simulation
             krz = zero
             if(i_job==3) then
                iz = modulo(i_grid_z,nptsz) + 1
                ! Extra phase if we extend outside simulation cell
                if(i_grid_z<1 .or. i_grid_z>nptsz) then
                   i_mult = -(i_grid_z - modulo(i_grid_z,nptsz))/nptsz
                   krz = kz(i_kp)*r_super_z*i_mult
                end if
             else
                iz = i_grid_z + 1
             end if
             do i_grid_y = miny, maxy
                ! Find y grid position and dy and wrap grid point
                dy = dg(2)*real(i_grid_y,double)+stm_y_min - atom_coord(2,i_atom)
                kry = zero
                iy = modulo(i_grid_y,nptsy)+1
                ! Extra phase if we extend outside simulation cell
                if(i_grid_y<1.or.i_grid_y>nptsy) then
                   i_mult = -(i_grid_y - modulo(i_grid_y,nptsy))/nptsy
                   kry = ky(i_kp)*r_super_y*i_mult
                end if
                do i_grid_x = minx, maxx
                   ! Find x grid position and dx and wrap grid point
                   dx = dg(1)*real(i_grid_x,double)+stm_x_min - atom_coord(1,i_atom)
                   krx = zero
                   ix = modulo(i_grid_x,nptsx)+1
                   ! Extra phase if we extend outside simulation cell
                   if(i_grid_x<1 .or. i_grid_x>nptsx) then
                      i_mult = -(i_grid_x - modulo(i_grid_x,nptsx))/nptsx
                      krx = kx(i_kp)*r_super_x*i_mult
                   end if
                   ! Calculate dr
                   dr = sqrt(dx*dx + dy*dy + dz*dz)
                   if(dr<=RadiusAtomf(i_spec)) then
                      phase = cmplx(cos(kr+krx+kry+krz),sin(kr+krx+kry+krz))
                      ! dr/dx = x/r etc.  Are the variable names confusing?
                      dx_dr = dx/dr
                      dy_dr = dy/dr
                      dz_dr = dz/dr
                      npao = 1
                      ! Loop over l
                      do i_l = 0, pao(i_spec)%greatest_angmom
                         ! Loop over zeta
                         do i_zeta = 1, pao(i_spec)%angmom(i_l)%n_zeta_in_angmom
                            ! Find f(r), df/dr
                            ! Loop over m
                            do i_m = -i_l, i_l
                               call evaluate_pao(i_spec,i_l,i_zeta,i_m,dx,dy,dz,f_r)
                               ! Accumulate into psi
                               psi(ix, iy, iz) = psi(ix, iy, iz) + &
                                    phase*evec_coeff(npao,i_atom,i_band, i_kp, i_spin) * f_r
                               call pao_elem_derivative_2(1,i_spec,i_l,i_zeta,i_m,dx,dy,dz,df_r)
                               dpsi(ix, iy, iz, 1) = dpsi(ix, iy, iz, 1) + &
                                    phase*evec_coeff(npao,i_atom,i_band, i_kp, i_spin) * df_r
                               call pao_elem_derivative_2(2,i_spec,i_l,i_zeta,i_m,dx,dy,dz,df_r)
                               dpsi(ix, iy, iz, 2) = dpsi(ix, iy, iz, 2) + &
                                    phase*evec_coeff(npao,i_atom,i_band, i_kp, i_spin) * df_r
                               call pao_elem_derivative_2(3,i_spec,i_l,i_zeta,i_m,dx,dy,dz,df_r)
                               dpsi(ix, iy, iz, 3) = dpsi(ix, iy, iz, 3) + &
                                    phase*evec_coeff(npao,i_atom,i_band, i_kp, i_spin) * df_r
                               npao = npao + 1
                            end do ! i_m
                         end do ! i_zeta
                      end do ! i_l
                   end if ! dr <= RadiusAtomf
                end do ! i_grid_x
             end do ! i_grid_y
          end do ! i_grid_z
       end if ! Atom is in STM region
    end do ! i_atom
    return
  end subroutine pao_dpao_to_grid

  subroutine pao_to_grid(i_band, i_kp, i_spin, psi)

    use datatypes
    use numbers
    use units
    use global_module, ONLY: ni_in_cell, atom_coord, species_glob
    use pao_format, ONLY: pao
    use local, ONLY: nptsx, nptsy, nptsz, grid_x, grid_y, grid_z, stm_x_min, stm_x_max, &
         stm_y_min, stm_y_max, stm_z_min, stm_z_max, evec_coeff, kx, ky, kz, i_job
    use dimens, ONLY: RadiusAtomf, r_super_x, r_super_y, r_super_z
    use angular_coeff_routines, ONLY: evaluate_pao

    implicit none

    ! Passed variables
    integer :: i_band, i_kp, i_spin
    complex(double_cplx), dimension(nptsx, nptsy, nptsz) :: psi

    ! Local variables
    integer :: i_atom, i_grid_x, i_grid_y, i_grid_z, i_l, i_zeta, i_m, ix, iy, iz
    integer :: i_spec, j, npao, i_mult
    integer :: minx, maxx, miny, maxy, minz, maxz
    real(double) :: dr, dx, dy, dz, sph_rl, f_r, df_r, dx_dr, dy_dr, dz_dr, del_r
    real(double) :: a, b, c, d, r1, r2, r3, r4, rr, kr, krx, kry, krz
    real(double), dimension(3) :: dsph_rl, dg
    complex(double_cplx) :: phase, phase_shift

    psi = zero
    ! Grid spacing
    dg(1) = grid_x!/BohrToAng
    dg(2) = grid_y!/BohrToAng
    dg(3) = grid_z!/BohrToAng
    ! Loop over atoms
    do i_atom = 1, ni_in_cell
       i_spec = species_glob(i_atom)
       if(atom_coord(3, i_atom) + RadiusAtomf(i_spec) >= stm_z_min) then ! Is the atom in STM region?
          kr = kx(i_kp)*atom_coord(1, i_atom) + ky(i_kp)*atom_coord(2, i_atom) + kz(i_kp)*atom_coord(3, i_atom)
          !phase = cmplx(cos(kr),sin(kr))
          ! Find grid limits
          minx = floor( (atom_coord(1, i_atom) - RadiusAtomf(i_spec))/dg(1) )
          maxx = floor( (atom_coord(1, i_atom) + RadiusAtomf(i_spec))/dg(1) ) + 1
          miny = floor( (atom_coord(2, i_atom) - RadiusAtomf(i_spec))/dg(2) )
          maxy = floor( (atom_coord(2, i_atom) + RadiusAtomf(i_spec))/dg(2) ) + 1
          minz = floor( (atom_coord(3, i_atom) - RadiusAtomf(i_spec))/dg(3) )
          maxz = floor( (atom_coord(3, i_atom) + RadiusAtomf(i_spec))/dg(3) ) + 1
          ! Account for STM limits
          if(i_job==4.or.i_job==5) then ! STM not band density
             if(stm_z_min>zero) then
                minz = minz - floor(stm_z_min/dg(3))
             end if
             if(minz<1) minz = 1
             maxz = min(maxz, nptsz)
          end if
          ! Loop over grid points
          do i_grid_z = minz, maxz
             ! Find z grid position and dz
             dz = dg(3)*real(i_grid_z,double)+stm_z_min - atom_coord(3,i_atom)
             ! Wrap if we're making band densities, but not for STM simulation
             krz = zero
             if(i_job==3) then
                iz = modulo(i_grid_z,nptsz) + 1
                ! Extra phase if we extend outside simulation cell
                if(i_grid_z<1 .or. i_grid_z>nptsz) then
                   i_mult = -(i_grid_z - modulo(i_grid_z,nptsz))/nptsz
                   krz = kz(i_kp)*r_super_z*i_mult
                end if
             else
                iz = i_grid_z + 1
             end if
             do i_grid_y = miny, maxy
                ! Find y grid position and dy and wrap grid point
                dy = dg(2)*real(i_grid_y,double)+stm_y_min - atom_coord(2,i_atom)
                kry = zero
                iy = modulo(i_grid_y,nptsy)+1
                ! Extra phase if we extend outside simulation cell
                if(i_grid_y<1.or.i_grid_y>nptsy) then
                   i_mult = -(i_grid_y - modulo(i_grid_y,nptsy))/nptsy
                   kry = ky(i_kp)*r_super_y*i_mult
                end if
                do i_grid_x = minx, maxx
                   ! Find x grid position and dx and wrap grid point
                   dx = dg(1)*real(i_grid_x,double)+stm_x_min - atom_coord(1,i_atom)
                   krx = zero
                   ix = modulo(i_grid_x,nptsx)+1
                   ! Extra phase if we extend outside simulation cell
                   if(i_grid_x<1 .or. i_grid_x>nptsx) then
                      i_mult = -(i_grid_x - modulo(i_grid_x,nptsx))/nptsx
                      krx = kx(i_kp)*r_super_x*i_mult
                   end if
                   ! Calculate dr
                   dr = sqrt(dx*dx + dy*dy + dz*dz)
                   if(dr<=RadiusAtomf(i_spec)) then
                      phase = cmplx(cos(kr+krx+kry+krz),sin(kr+krx+kry+krz))
                      npao = 1
                      ! Loop over l
                      do i_l = 0, pao(i_spec)%greatest_angmom
                         ! Loop over zeta
                         do i_zeta = 1, pao(i_spec)%angmom(i_l)%n_zeta_in_angmom
                            ! Loop over m
                            do i_m = -i_l, i_l
                               call evaluate_pao(i_spec,i_l,i_zeta,i_m,dx,dy,dz,f_r)
                               ! Accumulate into psi
                               psi(ix, iy, iz) = psi(ix, iy, iz) + &
                                    phase*evec_coeff(npao,i_atom,i_band, i_kp, i_spin) * f_r
                               npao = npao + 1
                            end do ! i_m
                         end do ! i_zeta
                      end do ! i_l
                   end if ! dr <= RadiusAtomf
                end do ! i_grid_x
             end do ! i_grid_y
          end do ! i_grid_z
       end if ! Atom is in STM region
    end do ! i_atom
    return
  end subroutine pao_to_grid

  subroutine read_domain(lun,proc,data)

    use datatypes
    use numbers
    use local, ONLY: block_store, nxmin, nymin, nzmin, current, nptsx, nptsy, nptsz
    use block_module, only: n_pts_in_block, in_block_x,in_block_y,in_block_z

    implicit none

    ! Passed
    integer :: lun, proc
    real(double), dimension(nptsx,nptsy,nptsz) :: data

    ! Local
    real(double), dimension(n_pts_in_block) :: local_grid
    integer :: iblock, point, ptx, pty, ptz, npx, npy, npz

    do iblock=1,block_store(proc)%num_blocks
       local_grid = zero
       ! Read block
       do point = 1,n_pts_in_block
          read(lun,*) local_grid(point)
       end do
       point = 0
       if(block_store(proc)%active(iblock)==1) then
          do ptz = 1, in_block_z
             npz = in_block_z*(block_store(proc)%nz(iblock)-1) + ptz - nzmin
             if(npz>0.AND.npz<=nptsz) then
                do pty = 1, in_block_y
                   npy = in_block_y*(block_store(proc)%ny(iblock)-1) + pty - nymin
                   if(npy>0.AND.npy<=nptsy) then
                      do ptx = 1, in_block_x
                         npx = in_block_x*(block_store(proc)%nx(iblock)-1) + ptx - nxmin
                         point = point + 1
                         if(npx>0.AND.npx<=nptsx) then
                            data(npx,npy,npz) = data(npx,npy,npz) + local_grid(point)
                         end if
                      end do
                   else
                      do ptx = 1, in_block_x
                         point = point + 1
                      end do
                   end if
                end do
             else
                do pty = 1, in_block_y
                   do ptx = 1, in_block_x
                      point = point + 1
                   end do
                end do
             end if
          end do
       end if ! if active block
    end do ! iblock
  end subroutine read_domain

  subroutine occupy(occu, ebands, Ef, spin)

    use datatypes
    use numbers
    use local, ONLY: nkp, n_bands_total, kT, flag_smear_type, iMethfessel_Paxton, wtk

    implicit none

    ! Passed variables
    real(double), dimension(:) :: Ef
    real(double), dimension(:,:) :: occu
    real(double), dimension(:,:,:) :: ebands
    integer :: spin

    ! Local variables
    integer :: i_kp, i_band

    do i_kp = 1, nkp
       do i_band = 1, n_bands_total
          select case(flag_smear_type)
          case (0) ! Fermi smearing
             occu(i_band,i_kp) = &
                  fermi(ebands(i_band,i_kp,spin) - Ef(spin), kT)
          case (1) ! Methfessel Paxton smearing
             occu(i_band,i_kp) = &
                  MP_step(ebands(i_band,i_kp,spin) - Ef(spin), &
                  iMethfessel_Paxton, kT)
          end select
       end do
    end do
    return
  end subroutine occupy
  ! -----------------------------------------------------------------------------
  ! Function fermi
  ! -----------------------------------------------------------------------------

  !!****f* DiagModule/fermi *
  !!
  !!  NAME
  !!   fermi - evaluate fermi function
  !!  USAGE
  !!   fermi(E,kT)
  !!  PURPOSE
  !!   Evaluates the fermi occupation of an energy
  !!
  !!   I'm assuming (for the sake of argument) that if both the energy and
  !!   the smearing (kT) are zero then we get an occupation of 0.5 - this is
  !!   certainly the limit if E and kT are equal and heading to zero, or if
  !!   E is smaller than kT and both head for zero.
  !!  INPUTS
  !!   real(double), intent(in) :: E - energy
  !!   real(double), intent(in) :: kT - smearing energy
  !!  USES
  !!   datatypes, numbers
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   23/04/2002
  !!  MODIFICATION HISTORY
  !!   2006/10/02 17:54 dave
  !!    Small fix to prevent maths overflows by only calculating
  !!    exponential if x well bounded
  !!   2012/01/22 L.Tong
  !!    Small change to use FORTRAN 90 function declaration notation
  !!    This works better with etags
  !!  SOURCE
  !!
  function fermi(E,kT)

    use datatypes
    use numbers, only: zero, one, half

    implicit none

    ! result
    real(double) :: fermi

    ! Passed variables
    real(double), intent(in) :: E
    real(double), intent(in) :: kT

    ! Local variables
    real(double) :: x
    real(double), parameter :: cutoff = 10.0_double

    if(kT==zero) then
       if(E>zero) then
          fermi = zero
       else if(E<zero) then
          fermi = one
       else if(E==zero) then
          fermi = half
       end if
    else
       x = E/kT
       if(x > cutoff) then
          fermi = zero
       elseif(x < -cutoff) then
          fermi = one
       else
          fermi = one/(one + exp(x))
       endif
    end if
  end function fermi
  !!***


  ! -----------------------------------------------------------------------------
  ! Function MP_step
  ! -----------------------------------------------------------------------------

  !!****f* DiagModule/MP_step *
  !!
  !!  NAME
  !!   MP_step - evaluate Methfessel-Paxton step function
  !!  USAGE
  !!   MP_step(E,order,smear)
  !!  PURPOSE
  !!   Evaluates the order (order) Methfessel-Paxton approximation to
  !!   step function
  !!  INPUTS
  !!   real(double), intent(in) :: E - energy
  !!   integer, intent(in) :: order - order of Methfessel expansion
  !!   real(double), intent(in) :: smear - smearing energy, nothing to
  !!                               do with physical temperature
  !!  USES
  !!   datatypes, numbers
  !!  AUTHOR
  !!   L.Tong (lt)
  !!  CREATION DATE
  !!   2010/06/15 00:17
  !!  MODIFICATION HISTORY
  !!   2012/01/22 L.Tong
  !!     - Small change to use FORTRAN 90 function declaration notation
  !!       This works better with etags
  !!  SOURCE
  !!
  function MP_step(E,order,smear)

    use datatypes
    use numbers, only: zero, one, half, two, four, pi
    use functions, ONLY: erfc_cq

    implicit none

    ! Result
    real(double) :: MP_step

    ! Passed variables
    real(double), intent(in) :: E
    integer,      intent(in) :: order
    real(double), intent(in) :: smear

    ! Internal variables
    real(double) :: x, A, H0, H1, H2, nd, x2
    integer      :: n

    ! in case of smear==0, we have the exact step function
    if(smear==zero) then
       if(E>zero) then
          MP_step = zero
       else if(E<zero) then
          MP_step = one
       else if(E==zero) then
          MP_step = half
       end if
    else if(smear>zero) then
       x = E/smear
       if(order==0) then
          MP_step = half*erfc_cq(x)
       else
          x2 = x*x
          A = one/sqrt(pi)
          H0 = one
          H1 = two*x
          MP_step = half*erfc_cq(x)
          nd = one
          do n=1,order
             A = A/((-four)*real(n,double))
             MP_step = MP_Step + A*H1*exp(-x2)
             H2 = two*x*H1 - two*nd*H0
             H0 = H1
             H1 = H2
             nd = nd + one
             H2 = two*x*H1 - two*nd*H0
             H0 = H1
             H1 = H2
             nd = nd + one
          end do
       end if
    end if
  end function MP_step
  !!***

end module process
