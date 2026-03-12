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
         flag_rotate_pdos, band_full_to_active, n_atoms_pDOS, pDOS_atom_index, &
         pdos_ax, pdos_ay, pdos_az
    use read, ONLY: read_eigenvalues, read_psi_coeffs, read_nprocs_from_blocks
    use global_module, ONLY: nspin, n_DOS, E_DOS_min, E_DOS_max, sigma_DOS, ni_in_cell, species_glob
    use units, ONLY: HaToeV
    use species_module, ONLY: nsf_species, npao_species
    use pao_format, ONLY: pao

    implicit none
    
    ! Local variables
    integer :: i_band, i_kp, i_spin, n_DOS_wid, n_band, n_min, n_max, i, i_atom,max_nsf, i_spec, &
         i_l, nzeta, sf_offset, max_l, norbs, i_m, i_band_c, i_z
    real(double) :: Ebin, dE_DOS, a, pf_DOS, spin_fac, coeff, check_electrons, peak_width
    real(double), dimension(:,:,:), allocatable :: pDOS
    real(double), dimension(:,:,:,:), allocatable :: pDOS_l
    real(double), dimension(:,:,:,:,:), allocatable :: pDOS_lm
    real(double), dimension(:,:), allocatable :: occ
    real(double), dimension(:,:), allocatable :: total_electrons
    real(double), dimension(:,:,:), allocatable :: total_electrons_l
    real(double) :: A1(3, 3), A2(5, 5), A1inv(3, 3), A2inv(5, 5)
    real(double) :: U1(3, 3), U2(5, 5), U1d(3, 3), U2d(5, 5)
    real(double) :: C1(3, 3), C2(5, 5)
    real(double) :: rod(3, 3), angle, axis(3)


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
    ! Call rotation subroutine if desired
    allocate(occ(n_bands_total,nkp))
    if (flag_rotate_pdos) then
      write(*,fmt='(2x, "Initialising A matrices")')
      call initialise_A_mat(A1, A2)
      ! ! Receive user input axes, calculate axis from eigenvector
      ! Normalise vectors, then check they form an orthonormal basis
      call get_pdos_axes
      call calculate_axis_angle(pdos_ax, pdos_ay, pdos_az, axis, angle)
      ! Make Rodrigues
      !axis = (/real(double) :: 1.0, 1.0, 1.0/)
      !angle = pi * 0.437547648
      call construct_rodrigues(axis, angle, rod)

      ! Construct all C^l matrices from rodrigues
      call construct_C1(rod, C1)
      call construct_C2(rod, C2)
      ! compute U^l
      call construct_Ul(1, A1, C1, U1)
      call construct_Ul(2, A2, C2, U2)
      U1d = transpose(U1) ! UU^T = I
      U2d = transpose(U2) ! UU^T = I
      write(*,fmt='(2x,"Rotating wavefunction coefficients")')
      call rotate_coefficients(U1d, U2d)
    end if 
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
                write(17,fmt=fmt_dos) HaToeV*(E_DOS_min + range_offset(i_spin) - efermi(i_spin) + dE_DOS*real(i-1,double)), &
                     pDOS(i_atom,i,i_spin), ((pDOS_lm(i_m,i_l,i_atom,i,i_spin),i_m=-i_l,i_l),i_l=0,pao(i_spec)%greatest_angmom)
             end do
          else if(flag_l_resolved) then
             write(17,fmt='("# Energy(eV)   pDOS(/eV)")')
             write(17,fmt='("#                  Total         l=0         l=1         l=2")')
             do i=1, n_DOS
                write(17,fmt=fmt_dos) HaToeV*(E_DOS_min + range_offset(i_spin) - efermi(i_spin) + dE_DOS*real(i-1,double)), &
                     pDOS(i_atom,i,i_spin), pDOS_l(0:pao(i_spec)%greatest_angmom,i_atom,i,i_spin)
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

  subroutine get_pdos_axes
   use datatypes
   use local, ONLY: pdos_ax, pdos_ay, pdos_az
   use dimens, ONLY: r_super_x, r_super_y, r_super_z
   implicit none

   real(double), parameter :: tol = 1e-10

   ! Normalise new x-axis correctly
   pdos_ax = pdos_ax / norm2(pdos_ax)
   ! Normalise new y-axis correctly
   pdos_ay = pdos_ay / norm2(pdos_ay)
   ! Normalise new y-axis correctly
   pdos_az = pdos_az / norm2(pdos_az)
   if (abs(dot_product(pdos_ax, pdos_ay)) > tol) &
      stop "pDOS_ax vector was not orthogonal to pDOS_ay"
   if (abs(dot_product(pdos_az, pdos_ay)) > tol) &
       stop "pDOS_az vector was not orthogonal to pDOS_ay"
   if (abs(dot_product(pdos_az, pdos_ax)) > tol) &
       stop "pDOS_az vector was not orthogonal to pDOS_ax"
  end subroutine

  subroutine initialise_A_mat(A1, A2)
   use datatypes
      implicit none
      real(double), intent(out) :: A1(3, 3), A2(5, 5)
      A1 = 0.0
      A2 = 0.0
      !A3 = 0.0
      ! A1, FOR p-orbitals, l = 1
      A1(1, 3) = 1.0
      A1(2, 1) = 1.0
      A1(3, 2) = 1.0
      ! d-orbitals, l = 2
      A2(1, 2) = 1.0
      A2(2, 4) = 1.0
      A2(3, 1) = 1.0
      A2(4, 5) = 2.0
      A2(5, 3) = 2.0*sqrt(3.0_double)
      ! f-orbitals, l = 3
      ! A3(1, 5) = 1
      ! A3(2, 3) = 1
      ! A3(3, 4) = 0.5*sqrt(6.0)
      ! A3(4, 2) = 0.05*sqrt(10.0)
      ! A3(5, 1) = 0.2*sqrt(15.0)
      ! A3(6, 7) = 0.2*sqrt(15.0)
      ! A3(7, 6) = 0.05*sqrt(10.0)
   end subroutine initialise_A_mat
   subroutine calculate_axis_angle(w1, w2, w3, axis, angle)
      use datatypes
      implicit none
      external DGEEV

      ! Allow the user to define new coordinate system
      ! Rotate from simulation axes (standard Cartesian) to local axes
      real(double), intent(in) ::  w1(3), w2(3), w3(3)
      real(double), intent(out) :: axis(3), angle
      real(double) :: pos_diff(3), cross_prod(3), det, sin_angle, cos_angle
      real(double) :: basis_matrix(3,3), temp(3,3), antisym_axis(3,3)
      ! Local variables
      integer :: i, j, N, LDA, LDVL, LDVR, INFO, LDWORK
      real(double), allocatable::  A(:,:), WR(:),   WI(:),   VL(:,:), VR(:,:), WORK(:), tra
      real(double), parameter :: tol = 1e-10
      N = 3
      LDA = N
      LDVL = N
      LDVR = N
      allocate(A(N,N))
      allocate(WR(N))
      allocate(WI(N))
      allocate(VL(LDVL, N))
      allocate(VR(LDVR, N))
      allocate(WORK(5*N))
      LDWORK = 5*N
      ! Define change of basis matrix;
      ! Columns of new basis in terms of old basis
      basis_matrix(:,1) = w1
      basis_matrix(:,2) = w2
      basis_matrix(:,3) = w3
      ! Make copy because DGEEV destroys matrix
      A = basis_matrix
      ! If determinant flips sign, orientation of the coordinate system has changed
      det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
      - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
      + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

      if (det  - 1.000 > tol) &
         stop "ERROR: determinant of change-of-basis matrix is negative"
      call DGEEV("N", "V", N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, &
           WORK, LDWORK, INFO)
      if (INFO /= 0) then
         write(*, fmt='(A,I0)') 'DGEEV INFO = ', INFO
         stop 'DGEEV, Eigenvalues of basis change matrix could not be found!'
      end if
      ! Require eigenvector with eigenvalue 1, no complex, to find axis
      do i = 1, N
      if (abs(WR(i) - 1.00000) < tol .and. abs(WI(i)) < tol) then
         axis  = VR(:, i)   ! Column i of VR is the i-th right eigenvector
         ! write(*,'(A,I0,A,F12.8,A,F12.8,A)') &
         ! '  Eigenvalue ', i, ' = (', WR(i), ', ', WI(i), 'i)'
         ! write(*,'(A,3F12.8)') '  Eigenvector = ', axis
         exit
      end if
   end do
   ! Calculate rotation angle using trace in radians
   ! Since only using acos gives -theta, theta, use atan2 to get correct signed angle
   ! by constructing sin theta as well
   ! Must check sign is correct later on
   tra = basis_matrix(1,1) + basis_matrix(2,2) + basis_matrix(3,3)
   
   antisym_axis = 0.0
   cos_angle = 0.5 * (tra - 1)
   antisym_axis(1, 2) = -axis(3)
   antisym_axis(1, 3) = axis(2)
   antisym_axis(2, 1) = axis(3)
   antisym_axis(2, 3) = -axis(1)
   antisym_axis(3, 1) = -axis(2)
   antisym_axis(3, 2) = axis(1)
   temp = matmul(antisym_axis, basis_matrix)
   sin_angle = -0.5 *(temp(1,1) + temp(2,2) + temp(3,3))
   ! Now angle in -[pi, pi] -> need to check how it affects Rodrigues
   angle = datan2(sin_angle, cos_angle)
   !print *, angle * 57.295779513
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
      !Normalise axis to unit vector 
      norm_axis = axis / norm2(axis) 
      ! Define K
      K(1, 2) = -norm_axis(3)
      K(1, 3) = norm_axis(2)
      K(2, 1) = norm_axis(3)
      K(2, 3) = -norm_axis(1)
      K(3, 1) = -norm_axis(2)
      K(3, 2) = norm_axis(1)

      KT = transpose(K)
      ! Rodrigues rotation matrix but with -angle
      temp = identity - (sin(angle)*KT) + (1 - cos(angle))*matmul(KT, KT)
      matrix = reshape(transpose(temp), shape=(/9/))
   end subroutine construct_rodrigues

   subroutine construct_C1(rod, C1)
      use datatypes
      implicit none
      real(double), intent(in) :: rod(9)
      real(double), intent(out) :: C1(3, 3)
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
   end subroutine construct_C1
   subroutine construct_C2(rod, C2)
      use datatypes
      implicit none
      real(double), intent(in) :: rod(9)
      real(double), intent(out) :: C2(5, 5)
      real(double) :: r(9)
      integer :: i
      do i = 1, size(r)
         r(i) = rod(i)
      end do
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
end subroutine construct_C2
   function inv(A) result(Ainv)
      use datatypes
      real(double), dimension(:, :), intent(in) :: A
      real(double), dimension(size(A, 1), size(A, 2)) :: Ainv

      real(double), dimension(size(A, 1)) :: work  ! work array for LAPACK
      integer, dimension(size(A, 1)) :: ipiv   ! pivot indices
      integer :: n, info

      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A, 1)

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
         stop 'Matrix inversion failed!'
      end if
   end function inv

   subroutine construct_Ul(angmom, Al, Cl, Ul)
      use datatypes
      implicit none
! U has dimensions(2l +1, 2l+1)
      integer, intent(in) :: angmom
      real(double), intent(in) ::   Al(2*angmom + 1, 2*angmom + 1), Cl(2*angmom + 1, 2*angmom + 1)
      real(double), intent(out) :: Ul(2*angmom + 1, 2*angmom + 1)
      real(double) :: Al_inv(2*angmom + 1, 2*angmom + 1)
      Al_inv = inv(Al)
      Ul = matmul(Al_inv, matmul(Cl, Al))
      
   end subroutine construct_Ul
   
   subroutine rotate_coefficients(U1, U2)
   use datatypes
   use local, ONLY: n_bands_total, nkp, n_atoms_pDOS, evec_coeff, scaled_evec_coeff, &
   n_atoms_pDOS, pDOS_atom_index, band_full_to_active
   use global_module, ONLY: nspin, species_glob
   use pao_format,    ONLY: pao

   implicit none

   real(double), intent(in) :: U1(3, 3)   ! rotation matrix for l=1
   real(double), intent(in) :: U2(5, 5)   ! rotation matrix for l=2

   ! Local variables
   integer :: i_atom, i_spec, i_band, i_kp, i_spin, i_band_c, g_atom
   integer :: i_l, i_z, nzeta, norbs, sf_offset

   do i_spin = 1, nspin
      do i_kp = 1, nkp
         do i_band = 1, n_bands_total
            i_band_c = band_full_to_active(i_band)
            do i_atom = 1, n_atoms_pDOS
               g_atom = pDOS_atom_index(i_atom)
               i_spec = species_glob(g_atom)
               sf_offset = 0
               ! Include l = 0 to correctly calculate offset
               do i_l = 0, pao(i_spec)%greatest_angmom
                  nzeta = pao(i_spec)%angmom(i_l)%n_zeta_in_angmom
                  norbs = 2*i_l + 1
                  do i_z = 1, nzeta
                  !evec_coeff(sf_offset,pDOS_atom_index(i_atom), i_band_c,i_kp,i_spin)
                     select case(i_l)
                     case(1)
                        evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin) = &
                           matmul(U1, evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin))
                        scaled_evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin) = &
                           matmul(U1, scaled_evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin))
                     case(2)
                        evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin) = &
                           matmul(U2, evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin))
                        scaled_evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin) = &
                           matmul(U2, scaled_evec_coeff(sf_offset+1:sf_offset+norbs, g_atom, i_band_c, i_kp, i_spin))
                     end select
                     sf_offset = sf_offset + norbs
                  end do ! i_z
               end do ! i_l
            end do ! i_atom
         end do ! i_band
      end do ! i_kp
   end do ! i_spin
   end subroutine rotate_coefficients

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
