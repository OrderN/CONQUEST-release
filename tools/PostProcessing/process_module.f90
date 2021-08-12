module process

  implicit none
  
contains

  subroutine assign_blocks

    use datatypes
    use local, ONLY: block_store, nprocs, block_size_x, block_size_y, block_size_z, &
         stm_z_min, stm_z_max, stm_x_min, stm_x_max, stm_y_min, stm_y_max
    
    implicit none

    integer :: proc, iblock, ig1, ind_group, block_x, block_y, block_z, nblock
    real(double) :: rbx, rby, rbz

    write(*,*) 'Range for x: ',stm_x_min, stm_x_max
    write(*,*) 'Range for y: ',stm_y_min, stm_y_max
    write(*,*) 'Range for z: ',stm_z_min, stm_z_max
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
    write(*,*) "Blocks assigned"
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

    write(*,*) 'Processing charge'
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
         n_bands_total, band_active_kp, flag_proc_range, band_full_to_active, &
         E_procwf_min, E_procwf_max, flag_procwf_range_Ef, band_proc_no, n_bands_process
    use output, ONLY: write_cube
    use global_module, only : nspin, flag_Multisite
    use read, ONLY: read_eigenvalues, read_psi_coeffs, read_MSSF_coeffs

    implicit none

    integer :: proc, band, nk, idum1, idum2, kp, ispin
    real(double) :: weight, rbx, rby, rbz, sq, test
    real(double), dimension(2) :: Emin, Emax
    character(len=50) :: filename, ci
    complex(double_cplx), dimension(:,:,:), allocatable :: psi, psi_sum

    write(*,*) 'Processing bands'
    ! Read eigenvalues
    call read_eigenvalues
    ! Read eigenvector coefficients
    call read_psi_coeffs
    if(flag_Multisite) call read_MSSF_coeffs
    !if(flag_Multisite) call write_one_MSSF
    
    allocate(current(nptsx,nptsy,nptsz))
    allocate(psi(nptsx,nptsy,nptsz))
    !allocate(psi_sum(nptsx,nptsy,nptsz))
    !psi_sum = zero
    if(flag_proc_range) then
       Emin = E_procwf_min
       Emax = E_procwf_max
       if(flag_procwf_range_Ef) then
          Emin = efermi + Emin
          Emax = efermi + Emax
       end if
       do ispin=1,nspin
          if(flag_by_kpoint) then ! Separate bands by k-point
             do band=1,n_bands_total
                current = zero
                do kp = 1,nkp
                   if(eigenvalues(band,kp,ispin) >= Emin(ispin) .and. &
                        eigenvalues(band,kp,ispin) <= Emax(ispin) .and. &
                        band_active_kp(band,kp,ispin)==1) then
                      if(flag_Multisite) then
                         call mssf_to_grid_alt(band_full_to_active(band), kp, ispin, psi)
                      else
                         call pao_to_grid(band_full_to_active(band), kp, ispin, psi)
                      end if
                      current = psi*conjg(psi)
                      write(ci,'("Band",I0.6,"den_kp",I0.3,"S",I0.1)') band, kp, ispin
                      call write_cube(current,ci)
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
                      !write(*,*) 'Band: ',band,eigenvalues(band,kp,ispin)
                      if(flag_Multisite) then
                         call mssf_to_grid_alt(band_full_to_active(band), kp, ispin, psi)
                      else
                         call pao_to_grid(band_full_to_active(band), kp, ispin, psi)
                      end if
                      !psi_sum = psi_sum + psi*conjg(psi)
                      current = current + psi*conjg(psi)
                      idum1 = 1
                   end if
                end do ! kp
                if(idum1==1) then
                   write(ci,'("Band",I0.6,"den_totS",I0.1)') band, ispin
                   call write_cube(current,ci)
                end if
             end do ! bands
          end if
       end do
    else ! User has provided list of bands
       do ispin=1,nspin
          if(flag_by_kpoint) then ! Separate bands by k-point
             do band=1,n_bands_process
                current = zero
                do kp = 1,nkp
                   ! This clause is needed in case the user chose an energy range that only selects some k-points
                   if(band_active_kp(band_proc_no(band),kp,ispin)==1) then
                      if(flag_Multisite) then
                         call mssf_to_grid_alt(band_full_to_active(band_proc_no(band)), kp, ispin, psi)
                      else
                         call pao_to_grid(band_full_to_active(band_proc_no(band)), kp, ispin, psi)
                      end if
                      current = psi*conjg(psi)
                      write(ci,'("Band",I0.6,"den_kp",I0.3,"S",I0.1)') band_proc_no(band), kp, ispin
                      call write_cube(current,ci)
                   end if
                end do ! kp
             end do ! bands = 1, n_bands_total
          else ! Sum over k-points
             do band=1,n_bands_process
                current = zero
                do kp = 1,nkp
                   ! This clause is needed in case the user chose an energy range that only selects some k-points
                   if(band_active_kp(band_proc_no(band),kp,ispin)==1) then
                      if(flag_Multisite) then
                         call mssf_to_grid_alt(band_full_to_active(band_proc_no(band)), kp, ispin, psi)
                      else
                         call pao_to_grid(band_full_to_active(band_proc_no(band)), kp, ispin, psi)
                      end if
                      current = current + psi*conjg(psi)
                   end if
                end do ! kp
                write(ci,'("Band",I0.6,"den_totS",I0.1)') band_proc_no(band), ispin
                call write_cube(current,ci)
             end do ! bands
          end if
       end do
    end if
    !current = zero
    !current = real(psi_sum)!*conjg(psi_sum)
    !write(ci,'("BandSum")') 
    !call write_cube(current,ci)
    return
  end subroutine process_bands

  subroutine process_mssf

    use datatypes
    use numbers
    use local, ONLY: mssf_proc_no, n_mssf_process
    use output, ONLY: write_cube
    use global_module, only : nspin, flag_Multisite
    use read, ONLY: read_eigenvalues, read_psi_coeffs, read_MSSF_coeffs

    implicit none

    integer :: i_mssf
    character(len=50) :: filename, ci
    complex(double_cplx), dimension(:,:,:), allocatable :: psi, psi_sum

    write(*,*) 'Processing MSSF'
    call read_MSSF_coeffs
    do i_mssf = 1, n_mssf_process
       write(*,*) 'Processing atom: ',mssf_proc_no(i_mssf)
       call write_one_MSSF(mssf_proc_no(i_mssf))
    end do
    return
  end subroutine process_mssf

  subroutine process_dos

    use datatypes
    use numbers, ONLY: zero, RD_ERR, twopi, half, one, two, four, six
    use local, ONLY: eigenvalues, n_bands_total, nkp, wtk, efermi, flag_total_iDOS
    use read, ONLY: read_eigenvalues, read_psi_coeffs
    use global_module, ONLY: nspin, n_DOS, E_DOS_min, E_DOS_max, sigma_DOS
    use units, ONLY: HaToeV

    implicit none
    
    ! Local variables
    integer :: i_band, i_kp, i_spin, n_DOS_wid, n_band, n_min, n_max, i
    real(double) :: Ebin, dE_DOS, a, pf_DOS, spin_fac
    real(double), dimension(nspin) :: total_electrons, iDOS_low
    real(double), dimension(:,:), allocatable :: total_DOS, iDOS
    real(double), dimension(:,:), allocatable :: occ

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
    ! Set limits on DOS output
    if(abs(E_DOS_min)<RD_ERR) then
       E_DOS_min = minval(eigenvalues(1,:,:))
       write(*,fmt='(2x,"DOS lower limit set automatically: ",f12.5," Ha")') E_DOS_min
    else
       write(*,fmt='(2x,"DOS lower limit set by user: ",f12.5," Ha")') E_DOS_min
    end if
    if(abs(E_DOS_max)<RD_ERR) then
       E_DOS_max = maxval(eigenvalues(n_bands_total,:,:))
       write(*,fmt='(2x,"DOS upper limit set automatically: ",f12.5," Ha")') E_DOS_max
    else
       write(*,fmt='(2x,"DOS upper limit set by user: ",f12.5," Ha")') E_DOS_max
    end if
    ! Spacing, width, prefactor
    dE_DOS = (E_DOS_max - E_DOS_min)/real(n_DOS-1,double)
    ! Set sigma automatically
    if(abs(sigma_DOS)<RD_ERR) then
       sigma_DOS = four*dE_DOS
       write(*,fmt='(2x,"Sigma set automatically: ",f12.6," Ha")') sigma_DOS
    else
       write(*,fmt='(2x,"Sigma set by user: ",f12.6," Ha")') sigma_DOS
       if(six*sigma_DOS < dE_DOS) write(*,fmt='(4x,"Sigma is much less than bin size: this may cause errors")')
    end if
    ! Adjust limits to allow full peak to be seen
    E_DOS_min = E_DOS_min - two*sigma_DOS
    E_DOS_max = E_DOS_max + two*sigma_DOS
    ! Recalculate dE_DOS now that we've broadened it
    dE_DOS = (E_DOS_max - E_DOS_min)/real(n_DOS-1,double)
    write(*,fmt='(2x,"Dividing DOS into ",i5," bins of width ",f12.6," Ha")') n_DOS, dE_DOS
    n_DOS_wid = floor(six*sigma_DOS/dE_DOS) ! How many bins either side of state we consider
    pf_DOS = one/(sigma_DOS*sqrt(twopi))
    total_electrons = zero
    iDOS_low = zero ! Integral of DOS to lowest energy bound
    ! Accumulate DOS over bands and k-points for each spin
    do i_spin = 1, nspin
       occ = zero
       call occupy(occ,eigenvalues,efermi,i_spin)
       do i_kp = 1, nkp
          do i_band=1,n_bands_total ! All bands
             if(eigenvalues(i_band, i_kp, i_spin)>E_DOS_min .and. &
                  eigenvalues(i_band, i_kp, i_spin)<E_DOS_max) then
                n_band = floor((eigenvalues(i_band, i_kp, i_spin) - E_DOS_min)/dE_DOS) + 1
                n_min = n_band - n_DOS_wid
                if(n_min<1) n_min = 1
                n_max = n_band + n_DOS_wid
                if(n_max>n_DOS) n_max = n_DOS
                do i = n_min, n_max
                   Ebin = real(i-1,double)*dE_DOS + E_DOS_min
                   a = (Ebin-eigenvalues(i_band, i_kp, i_spin))/sigma_DOS
                   total_DOS(i,i_spin) = total_DOS(i,i_spin) + wtk(i_kp)*pf_DOS*exp(-half*a*a)
                   total_electrons(i_spin) = total_electrons(i_spin) + occ(i_band,i_spin)*wtk(i_kp)*pf_DOS*exp(-half*a*a)
                   iDOS(i,i_spin) = iDOS(i,i_spin) + wtk(i_kp)*pf_DOS*exp(-half*a*a)
                end do
             else if(eigenvalues(i_band, i_kp, i_spin)<E_DOS_min) then
                iDOS_low(i_spin) = iDOS_low(i_spin) + wtk(i_kp)
             end if
          end do
       end do
       ! Now integrate DOS
       write(*,*) 'iDOS_low is ',iDOS_low
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
    if(nspin==1) then
       write(*,fmt='(2x,"DOS between ",f11.3," and ",f11.3," Ha integrates to ",f12.3," electrons")') &
            E_DOS_min, E_DOS_max, total_electrons(1)
    else
       write(*,fmt='(2x,"Spin Up DOS integrates to ",f12.3," electrons")') total_electrons(1)
       write(*,fmt='(2x,"Spin Dn DOS integrates to ",f12.3," electrons")') total_electrons(2)
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
          write(17,fmt='(3f14.5)') HaToeV*(E_DOS_min + dE_DOS*real(i-1,double)-efermi(i_spin)), &
               total_DOS(i,i_spin), iDOS(i,i_spin)
       end do
       write(17,fmt='("&")')
    end do
    close(unit=17)
    return
  end subroutine process_dos

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

    implicit none

    ! Passed variables
    integer :: i_band, i_kp, i_spin
    complex(double_cplx), dimension(nptsx, nptsy, nptsz) :: psi
    complex(double_cplx), dimension(nptsx, nptsy, nptsz, 3) :: dpsi

    ! Local variables
    integer :: i_atom, i_grid_x, i_grid_y, i_grid_z, i_l, i_zeta, i_m, ix, iy, iz
    integer :: i_spec, j, npao
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
       !write(*,*) 'Atom, spec, pos: ',i_atom, i_spec, atom_coord(3, i_atom) + RadiusAtomf(i_spec), stm_z_min
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
          !if(stm_x_min>zero) minx = max(minx, floor(stm_x_min/dg(1)) )
          !if(stm_x_max<r_super_x) maxx = min(maxx, floor(stm_x_max/dg(1)+1) )
          !if(stm_y_min>zero) miny = max(miny, floor(stm_y_min/dg(2)) )
          !if(stm_y_max<r_super_y) maxy = min(maxy, floor(stm_y_max/dg(2)+1) )
          if(i_job==4.or.i_job==5) then ! STM not band density, so no z periodicity
             if(stm_z_min>zero) then
                minz = minz - floor(stm_z_min/dg(3))
             end if
             if(minz<1) minz = 1
             maxz = min(maxz, nptsz)
          end if
          !end if
          ! Loop over grid points
          do i_grid_z = minz, maxz
             ! Find z grid position and dz
             dz = dg(3)*real(i_grid_z-1,double)+stm_z_min - atom_coord(3,i_atom)
             iz = i_grid_z
             ! Wrap if we're making band densities, but not for STM simulation
             krz = zero
             if(i_job==3) then
                if(i_grid_z<1) then
                   iz = i_grid_z + nptsz
                   krz = kz(i_kp)*r_super_z
                end if
                if(i_grid_z>nptsz) then
                   iz = i_grid_z - nptsz
                   krz = -kz(i_kp)*r_super_z
                end if
             end if
             do i_grid_y = miny, maxy
                ! Find y grid position and dy and wrap grid point
                dy = dg(2)*real(i_grid_y-1,double)+stm_y_min - atom_coord(2,i_atom)
                iy = i_grid_y
                kry = zero
                if(i_grid_y<1) then
                   iy = i_grid_y + nptsy
                   kry = ky(i_kp)*r_super_y
                end if
                if(i_grid_y>nptsy) then
                   iy = i_grid_y - nptsy
                   kry = -ky(i_kp)*r_super_y
                end if
                do i_grid_x = minx, maxx
                   ! Find x grid position and dx and wrap grid point
                   dx = dg(1)*real(i_grid_x-1,double)+stm_x_min - atom_coord(1,i_atom)
                   ix = i_grid_x
                   krx = zero
                   if(i_grid_x<1) then
                      ix = i_grid_x + nptsx
                      krx = kx(i_kp)*r_super_x
                   end if
                   if(i_grid_x>nptsx) then
                      ix = i_grid_x - nptsx
                      krx = -kx(i_kp)*r_super_x
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
                            del_r = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%delta
                            j = floor(dr/del_r) + 1
                            if(j+1<=pao(i_spec)%angmom(i_l)%zeta(i_zeta)%length) then
                               rr = real(j,double)*del_r
                               a = (rr - dr)/del_r
                               b = one - a
                               c = a * ( a * a - one ) * del_r * del_r / six
                               d = b * ( b * b - one ) * del_r * del_r / six
                               r1 = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%table(j)
                               r2 = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%table(j+1)
                               r3 = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%table2(j)
                               r4 = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%table2(j+1)
                               f_r = a*r1 + b*r2 + c*r3 + d*r4
                               ! Re-use a and b before redefining
                               c = -del_r*(three*a*a-one)/six
                               d =  del_r*(three*b*b-one)/six
                               a = -one/del_r
                               b = -a
                               df_r = a*r1 + b*r2 + c*r3 + d*r4
                            else
                               f_r = zero
                               df_r = zero
                            end if
                            ! Loop over m
                            do i_m = -i_l, i_l
                               ! Get spherical harmonic and gradient
                               ! NB returns r^l times spherical harmonic
                               call spherical_harmonic_rl(dx, dy, dz, i_l, i_m, sph_rl, dsph_rl)
                               ! Accumulate into psi, grad psi
                               psi(ix, iy, iz) = psi(ix, iy, iz) + &
                                    phase*evec_coeff(npao,i_atom,i_band, i_kp, i_spin) * sph_rl * f_r
                               dpsi(ix, iy, iz, 1) = dpsi(ix, iy, iz, 1) + &
                                    phase*evec_coeff(npao,i_atom,i_band, i_kp, i_spin) * &
                                    (sph_rl * df_r * dx_dr + dsph_rl(1) * f_r)
                               dpsi(ix, iy, iz, 2) = dpsi(ix, iy, iz, 2) + &
                                    phase*evec_coeff(npao,i_atom,i_band, i_kp, i_spin) * &
                                    (sph_rl * df_r * dy_dr + dsph_rl(2) * f_r)
                               dpsi(ix, iy, iz, 3) = dpsi(ix, iy, iz, 3) + &
                                    phase*evec_coeff(npao,i_atom,i_band, i_kp, i_spin) * &
                                    (sph_rl * df_r * dz_dr + dsph_rl(3) * f_r)
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

  ! This is really band_to_grid as it loops over all atoms and PAOs
  subroutine pao_to_grid(i_band, i_kp, i_spin, psi)

    use datatypes
    use numbers
    use units
    use global_module, ONLY: ni_in_cell, atom_coord, species_glob
    use pao_format, ONLY: pao
    use local, ONLY: nptsx, nptsy, nptsz, grid_x, grid_y, grid_z, stm_x_min, stm_x_max, &
         stm_y_min, stm_y_max, stm_z_min, stm_z_max, evec_coeff, kx, ky, kz, i_job
    use dimens, ONLY: RadiusAtomf, r_super_x, r_super_y, r_super_z

    implicit none

    ! Passed variables
    integer :: i_band, i_kp, i_spin
    complex(double_cplx), dimension(nptsx, nptsy, nptsz) :: psi

    ! Local variables
    integer :: i_atom, i_grid_x, i_grid_y, i_grid_z, i_l, i_zeta, i_m, ix, iy, iz
    integer :: i_spec, j, npao
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
       !write(*,*) 'Atom, spec, pos: ',i_atom, i_spec, atom_coord(3, i_atom) + RadiusAtomf(i_spec), stm_z_min
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
          !if(stm_x_min>zero.and.minx<0) then
          !   if(minx+nptsx>floor(stm_x_min
          !if(stm_x_min>zero) minx = max(minx, floor(stm_x_min/dg(1)) )
          !if(stm_x_max<r_super_x) maxx = min(maxx, floor(stm_x_max/dg(1)+1) )
          !if(stm_y_min>zero) miny = max(miny, floor(stm_y_min/dg(2)) )
          !if(stm_y_max<r_super_y) maxy = min(maxy, floor(stm_y_max/dg(2)+1) )
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
             dz = dg(3)*real(i_grid_z-1,double)+stm_z_min - atom_coord(3,i_atom)
             iz = i_grid_z
             ! Wrap if we're making band densities, but not for STM simulation
             krz = zero
             if(i_job==3) then
                if(i_grid_z<1) then
                   iz = i_grid_z + nptsz
                   krz = kz(i_kp)*r_super_z
                end if
                if(i_grid_z>nptsz) then
                   iz = i_grid_z - nptsz
                   krz = -kz(i_kp)*r_super_z
                end if
             end if
             do i_grid_y = miny, maxy
                ! Find y grid position and dy and wrap grid point
                dy = dg(2)*real(i_grid_y-1,double)+stm_y_min - atom_coord(2,i_atom)
                iy = i_grid_y
                kry = zero
                if(i_grid_y<1) then
                   iy = i_grid_y + nptsy
                   kry = ky(i_kp)*r_super_y
                end if
                if(i_grid_y>nptsy) then
                   iy = i_grid_y - nptsy
                   kry = -ky(i_kp)*r_super_y
                end if
                do i_grid_x = minx, maxx
                   ! Find x grid position and dx and wrap grid point
                   dx = dg(1)*real(i_grid_x-1,double)+stm_x_min - atom_coord(1,i_atom)
                   ix = i_grid_x
                   krx = zero
                   if(i_grid_x<1) then
                      ix = i_grid_x + nptsx
                      krx = kx(i_kp)*r_super_x
                   end if
                   if(i_grid_x>nptsx) then
                      ix = i_grid_x - nptsx
                      krx = -kx(i_kp)*r_super_x
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
                            ! Find f(r), df/dr
                            del_r = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%delta
                            j = floor(dr/del_r) + 1
                            if(j+1<=pao(i_spec)%angmom(i_l)%zeta(i_zeta)%length) then
                               rr = real(j,double)*del_r
                               a = (rr - dr)/del_r
                               b = one - a
                               c = a * ( a * a - one ) * del_r * del_r / six
                               d = b * ( b * b - one ) * del_r * del_r / six
                               r1 = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%table(j)
                               r2 = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%table(j+1)
                               r3 = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%table2(j)
                               r4 = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%table2(j+1)
                               f_r = a*r1 + b*r2 + c*r3 + d*r4
                            else
                               f_r = zero
                            end if
                            ! Loop over m
                            do i_m = -i_l, i_l
                               ! Get spherical harmonic and gradient
                               ! NB returns r^l times spherical harmonic
                               call spherical_harmonic_rl(dx, dy, dz, i_l, i_m, sph_rl, dsph_rl)
                               ! Accumulate into psi, grad psi
                               psi(ix, iy, iz) = psi(ix, iy, iz) + &
                                    phase*evec_coeff(npao,i_atom,i_band, i_kp, i_spin) * sph_rl * f_r
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

  ! This is really band_to_grid as it loops over all atoms and MSSFs
  subroutine mssf_to_grid(i_band, i_kp, i_spin, psi)

    use datatypes
    use numbers
    use units
    use global_module, ONLY: ni_in_cell, atom_coord, species_glob
    use pao_format, ONLY: pao
    use local, ONLY: nptsx, nptsy, nptsz, grid_x, grid_y, grid_z, stm_x_min, stm_x_max, &
         stm_y_min, stm_y_max, stm_z_min, stm_z_max, evec_coeff, kx, ky, kz, i_job, &
         nab_glob, MSSF_coeffs
    use dimens, ONLY: RadiusAtomf, r_super_x, r_super_y, r_super_z, RadiusSupport

    implicit none

    ! Passed variables
    integer :: i_band, i_kp, i_spin
    complex(double_cplx), dimension(nptsx, nptsy, nptsz) :: psi

    ! Local variables
    integer :: i_atom, i_grid_x, i_grid_y, i_grid_z, i_l, i_zeta, i_m, ix, iy, iz
    integer :: i_spec, j, npao, i_mssf, i_nsf, i_glob_nabi, i_nabi
    integer :: minx, maxx, miny, maxy, minz, maxz
    real(double) :: dr, dx, dy, dz, sph_rl, f_r, df_r, dx_dr, dy_dr, dz_dr, del_r
    real(double) :: a, b, c, d, r1, r2, r3, r4, rr, kr, krx, kry, krz, b_ia_jb
    real(double) :: dx_mssf, dy_mssf, dz_mssf, dr_mssf, rem
    real(double), dimension(3) :: dsph_rl, dg
    complex(double_cplx) :: phase, phase_shift, c_ia_nk
    
    psi = zero
    ! Grid spacing
    dg(1) = grid_x!/BohrToAng
    dg(2) = grid_y!/BohrToAng
    dg(3) = grid_z!/BohrToAng
    ! Loop over atoms
    do i_atom = 1, ni_in_cell
       write(*,*) 'Atom: ',i_atom
       i_spec = species_glob(i_atom)
       !phase = cmplx(cos(kr),sin(kr))
       ! Find grid limits
       minx = floor( (atom_coord(1, i_atom) - RadiusAtomf(i_spec))/dg(1) )    
       maxx = floor( (atom_coord(1, i_atom) + RadiusAtomf(i_spec))/dg(1) ) + 1
       miny = floor( (atom_coord(2, i_atom) - RadiusAtomf(i_spec))/dg(2) )    
       maxy = floor( (atom_coord(2, i_atom) + RadiusAtomf(i_spec))/dg(2) ) + 1
       minz = floor( (atom_coord(3, i_atom) - RadiusAtomf(i_spec))/dg(3) )    
       maxz = floor( (atom_coord(3, i_atom) + RadiusAtomf(i_spec))/dg(3) ) + 1
       ! Loop over grid points
       do i_grid_z = minz, maxz
          ! Find z grid position and dz
          dz = dg(3)*real(i_grid_z-1,double) - atom_coord(3,i_atom)
          iz = i_grid_z
          krz = zero
          if(i_grid_z<1 .or. i_grid_z>nptsz) then
             rem = mod(i_grid_z,nptsz)
             if(i_grid_z<1) rem = rem + nptsz
             iz = rem
             krz = kz(i_kp)*r_super_z*(rem - i_grid_z)/nptsz
          end if
          do i_grid_y = miny, maxy
             ! Find y grid position and dy and wrap grid point
             dy = dg(2)*real(i_grid_y-1,double) - atom_coord(2,i_atom)
             iy = i_grid_y
             kry = zero
             if(i_grid_y<1 .or. i_grid_y>nptsy) then
                rem = mod(i_grid_y,nptsy)
                if(i_grid_y<1) rem = rem + nptsy
                iy = rem
                kry = ky(i_kp)*r_super_y*(rem - i_grid_y)/nptsy
             end if
             do i_grid_x = minx, maxx
                ! Find x grid position and dx and wrap grid point
                dx = dg(1)*real(i_grid_x-1,double) - atom_coord(1,i_atom)
                ix = i_grid_x
                krx = zero
                if(i_grid_x<1 .or. i_grid_x>nptsx) then
                   rem = mod(i_grid_x,nptsx)
                   if(i_grid_x<1) rem = rem + nptsx
                   ix = rem
                   krx = kx(i_kp)*r_super_x*(rem - i_grid_x)/nptsx
                end if
                ! Calculate dr
                dr = sqrt(dx*dx + dy*dy + dz*dz)
                if(dr<=RadiusAtomf(i_spec)) then
                   npao = 1
                   ! Loop over l
                   do i_l = 0, pao(i_spec)%greatest_angmom
                      ! Loop over zeta
                      do i_zeta = 1, pao(i_spec)%angmom(i_l)%n_zeta_in_angmom
                         ! Find f(r), df/dr
                         del_r = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%delta
                         j = floor(dr/del_r) + 1
                         if(j+1<=pao(i_spec)%angmom(i_l)%zeta(i_zeta)%length) then
                            rr = real(j,double)*del_r
                            a = (rr - dr)/del_r
                            b = one - a
                            c = a * ( a * a - one ) * del_r * del_r / six
                            d = b * ( b * b - one ) * del_r * del_r / six
                            r1 = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%table(j)
                            r2 = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%table(j+1)
                            r3 = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%table2(j)
                            r4 = pao(i_spec)%angmom(i_l)%zeta(i_zeta)%table2(j+1)
                            f_r = a*r1 + b*r2 + c*r3 + d*r4
                         else
                            f_r = zero
                         end if
                         ! Loop over m
                         do i_m = -i_l, i_l
                            c_ia_nk = evec_coeff(npao,i_atom,i_band, i_kp, i_spin)
                            ! Get spherical harmonic and gradient
                            ! NB returns r^l times spherical harmonic
                            call spherical_harmonic_rl(dx, dy, dz, i_l, i_m, sph_rl, dsph_rl)
                            ! Loop over MSSF we contribute to
                            do i_mssf = 1, nab_glob(i_atom)%i_count
                               i_glob_nabi = nab_glob(i_atom)%i_glob(i_mssf)
                               kr = kx(i_kp)*atom_coord(1, nab_glob(i_atom)%i_glob(i_mssf)) + &
                                    ky(i_kp)*atom_coord(2, nab_glob(i_atom)%i_glob(i_mssf)) + &
                                    kz(i_kp)*atom_coord(3, nab_glob(i_atom)%i_glob(i_mssf))
                               phase = cmplx(cos(kr+krx+kry+krz),sin(kr+krx+kry+krz))
                               i_nabi = nab_glob(i_atom)%neigh(i_mssf)
                               ! Calculate displacement using dx, dy, dz and check range
                               dx_mssf = dx - nab_glob(i_atom)%disp(i_mssf,1)
                               dy_mssf = dy - nab_glob(i_atom)%disp(i_mssf,2)
                               dz_mssf = dz - nab_glob(i_atom)%disp(i_mssf,3)
                               dr_mssf = sqrt(dx_mssf*dx_mssf + dy_mssf*dy_mssf + dz_mssf*dz_mssf)
                               !write(*,*) 'Loc: ',dx_mssf, dy_mssf, dz_mssf
                               if(dr_mssf < RadiusSupport(species_glob(i_glob_nabi))) then
                                  ! Loop over MSSFs on atom
                                  do i_nsf = 1, MSSF_coeffs(i_glob_nabi)%n_sf
                                     ! Coeff
                                     b_ia_jb = MSSF_coeffs(i_glob_nabi)% &
                                          neigh_coeff(i_nabi)%coeff(npao,i_nsf)
                                     psi(ix, iy, iz) = psi(ix, iy, iz) + &
                                          phase * c_ia_nk * b_ia_jb * sph_rl * f_r
                                  end do ! i_nsf
                               end if
                            end do ! i_mssf
                            npao = npao + 1
                         end do ! i_m
                      end do ! i_zeta
                   end do ! i_l
                end if ! dr <= RadiusAtomf
             end do ! i_grid_x
          end do ! i_grid_y
       end do ! i_grid_z
    end do ! i_atom
    return
  end subroutine mssf_to_grid

  ! This is really band_to_grid as it loops over all atoms and MSSFs
  ! Loop over c_ia and then over neighbours jb - less efficient but
  ! should be simpler indexing
  subroutine mssf_to_grid_alt(i_band, i_kp, i_spin, psi)

    use datatypes
    use numbers
    use units
    use global_module, ONLY: ni_in_cell, atom_coord, species_glob
    use pao_format, ONLY: pao
    use local, ONLY: nptsx, nptsy, nptsz, grid_x, grid_y, grid_z, stm_x_min, stm_x_max, &
         stm_y_min, stm_y_max, stm_z_min, stm_z_max, evec_coeff, kx, ky, kz, i_job, &
         nab_glob, MSSF_coeffs
    use dimens, ONLY: RadiusAtomf, r_super_x, r_super_y, r_super_z, RadiusSupport

    implicit none

    ! Passed variables
    integer :: i_band, i_kp, i_spin
    complex(double_cplx), dimension(nptsx, nptsy, nptsz) :: psi

    ! Local variables
    integer :: i_atom, i_grid_x, i_grid_y, i_grid_z, i_l, i_zeta, i_m, ix, iy, iz
    integer :: i_spec, j, npao, i_mssf, i_nsf, i_glob_nabi, i_nabi
    integer :: minx, maxx, miny, maxy, minz, maxz, j_spec
    real(double) :: dr, dx, dy, dz, sph_rl, f_r, df_r, dx_dr, dy_dr, dz_dr, del_r
    real(double) :: a, b, c, d, r1, r2, r3, r4, rr, kr, krx, kry, krz, b_ia_jb
    real(double) :: dx_mssf, dy_mssf, dz_mssf, dr_mssf
    real(double), dimension(3) :: dsph_rl, dg
    complex(double_cplx) :: phase, phase_shift, c_ia_nk
    
    psi = zero
    ! Grid spacing
    dg(1) = grid_x!/BohrToAng
    dg(2) = grid_y!/BohrToAng
    dg(3) = grid_z!/BohrToAng
    ! Loop over atoms
    do i_atom = 1, ni_in_cell
       i_spec = species_glob(i_atom)
       kr = kx(i_kp)*atom_coord(1, i_atom) + ky(i_kp)*atom_coord(2, i_atom) + kz(i_kp)*atom_coord(3, i_atom)
       ! Find grid limits
       minx = floor( (atom_coord(1, i_atom) - RadiusSupport(i_spec))/dg(1) )
       maxx = floor( (atom_coord(1, i_atom) + RadiusSupport(i_spec))/dg(1) ) + 1
       miny = floor( (atom_coord(2, i_atom) - RadiusSupport(i_spec))/dg(2) )    
       maxy = floor( (atom_coord(2, i_atom) + RadiusSupport(i_spec))/dg(2) ) + 1
       minz = floor( (atom_coord(3, i_atom) - RadiusSupport(i_spec))/dg(3) )    
       maxz = floor( (atom_coord(3, i_atom) + RadiusSupport(i_spec))/dg(3) ) + 1
       ! Loop over grid points
       do i_grid_z = minz, maxz
          ! Find z grid position and dz
          dz = dg(3)*real(i_grid_z-1,double) - atom_coord(3,i_atom)
          iz = i_grid_z
          krz = zero
          if(i_grid_z<1) then
             iz = i_grid_z + nptsz
             krz = kz(i_kp)*r_super_z
          end if
          if(i_grid_z>nptsz) then
             iz = i_grid_z - nptsz
             krz = -kz(i_kp)*r_super_z
          end if
          do i_grid_y = miny, maxy
             ! Find y grid position and dy and wrap grid point
             dy = dg(2)*real(i_grid_y-1,double) - atom_coord(2,i_atom)
             iy = i_grid_y
             kry = zero
             if(i_grid_y<1) then
                iy = i_grid_y + nptsy
                kry = ky(i_kp)*r_super_y
             end if
             if(i_grid_y>nptsy) then
                iy = i_grid_y - nptsy
                kry = -ky(i_kp)*r_super_y
             end if
             do i_grid_x = minx, maxx
                ! Find x grid position and dx and wrap grid point
                dx = dg(1)*real(i_grid_x-1,double) - atom_coord(1,i_atom)
                ix = i_grid_x
                krx = zero
                if(i_grid_x<1) then
                   ix = i_grid_x + nptsx
                   krx = kx(i_kp)*r_super_x
                end if
                if(i_grid_x>nptsx) then
                   ix = i_grid_x - nptsx
                   krx = -kx(i_kp)*r_super_x
                end if
                ! Calculate dr
                dr = sqrt(dx*dx + dy*dy + dz*dz)
                if(dr < RadiusSupport(i_spec)) then
                   phase = cmplx(cos(kr+krx+kry+krz),sin(kr+krx+kry+krz))
                   do i_mssf = 1, MSSF_coeffs(i_atom)%n_neighbours
                      j_spec = species_glob(MSSF_coeffs(i_atom)%neigh_coeff(i_mssf)%i_glob)
                      dx_mssf = dx - MSSF_coeffs(i_atom)%neigh_coeff(i_mssf)%dx
                      dy_mssf = dy - MSSF_coeffs(i_atom)%neigh_coeff(i_mssf)%dy
                      dz_mssf = dz - MSSF_coeffs(i_atom)%neigh_coeff(i_mssf)%dz
                      dr_mssf = sqrt(dx_mssf*dx_mssf + dy_mssf*dy_mssf + dz_mssf*dz_mssf)
                      if(dr_mssf<=RadiusAtomf(j_spec)) then
                         npao = 1
                         ! Loop over l
                         do i_l = 0, pao(j_spec)%greatest_angmom
                            ! Loop over zeta
                            do i_zeta = 1, pao(j_spec)%angmom(i_l)%n_zeta_in_angmom
                               ! Find f(r), df/dr
                               del_r = pao(j_spec)%angmom(i_l)%zeta(i_zeta)%delta
                               j = floor(dr_mssf/del_r) + 1
                               if(j+1<=pao(j_spec)%angmom(i_l)%zeta(i_zeta)%length) then
                                  rr = real(j,double)*del_r
                                  a = (rr - dr_mssf)/del_r
                                  b = one - a
                                  c = a * ( a * a - one ) * del_r * del_r / six
                                  d = b * ( b * b - one ) * del_r * del_r / six
                                  r1 = pao(j_spec)%angmom(i_l)%zeta(i_zeta)%table(j)
                                  r2 = pao(j_spec)%angmom(i_l)%zeta(i_zeta)%table(j+1)
                                  r3 = pao(j_spec)%angmom(i_l)%zeta(i_zeta)%table2(j)
                                  r4 = pao(j_spec)%angmom(i_l)%zeta(i_zeta)%table2(j+1)
                                  f_r = a*r1 + b*r2 + c*r3 + d*r4
                               else
                                  f_r = zero
                               end if
                               ! Loop over m
                               do i_m = -i_l, i_l
                                  ! Get spherical harmonic and gradient
                                  ! NB returns r^l times spherical harmonic
                                  call spherical_harmonic_rl(dx_mssf, dy_mssf, dz_mssf, i_l, i_m, sph_rl, dsph_rl)
                                  ! Loop over MSSF we contribute to
                                  ! Loop over MSSFs on atom
                                  do i_nsf = 1, MSSF_coeffs(i_atom)%n_sf
                                     ! Coeff
                                     c_ia_nk = evec_coeff(i_nsf,i_atom,i_band, i_kp, i_spin)
                                     b_ia_jb = MSSF_coeffs(i_atom)%neigh_coeff(i_mssf)%coeff(npao,i_nsf)
                                     psi(ix, iy, iz) = psi(ix, iy, iz) + &
                                          phase * c_ia_nk * b_ia_jb * sph_rl * f_r
                                  end do ! i_nsf
                                  npao = npao + 1
                               end do ! i_m
                            end do ! i_zeta
                         end do ! i_l
                      end if ! dr <= RadiusAtomf
                   end do
                end if ! dr_mssf < RadiusSupport
             end do ! i_grid_x
          end do ! i_grid_y
       end do ! i_grid_z
    end do ! i_atom
    return
  end subroutine mssf_to_grid_alt

  ! This is really band_to_grid as it loops over all atoms and MSSFs
  subroutine write_one_MSSF(i_atom)

    use datatypes
    use numbers
    use units
    use global_module, ONLY: ni_in_cell, atom_coord, species_glob
    use pao_format, ONLY: pao
    use local, ONLY: nptsx, nptsy, nptsz, grid_x, grid_y, grid_z, stm_x_min, stm_x_max, &
         stm_y_min, stm_y_max, stm_z_min, stm_z_max, evec_coeff, kx, ky, kz, i_job, &
         nab_glob, MSSF_coeffs
    use dimens, ONLY: RadiusAtomf, r_super_x, r_super_y, r_super_z, RadiusMS, RadiusSupport
    use output, ONLY: write_mssf_cube
    
    implicit none

    ! Passed variables
    integer :: i_atom

    ! Local variables
    integer :: i_grid_x, i_grid_y, i_grid_z, i_l, i_zeta, i_m, ix, iy, iz
    integer :: i_spec, j, npao, i_mssf, i_nsf, i_nab, j_spec, npx, npy, npz
    integer :: minx, maxx, miny, maxy, minz, maxz
    real(double), dimension(:,:,:), allocatable :: psi
    real(double) :: dr, dx, dy, dz, sph_rl, f_r, df_r, dx_dr, dy_dr, dz_dr, del_r
    real(double) :: a, b, c, d, r1, r2, r3, r4, rr, kr, krx, kry, krz, c_ia_nk, b_ia_jb
    real(double) :: dx_mssf, dy_mssf, dz_mssf
    real(double), dimension(3) :: dsph_rl, dg
    complex(double_cplx) :: phase, phase_shift
    character(len=50) :: ci

    ! Grid spacing
    dg(1) = grid_x!/BohrToAng
    dg(2) = grid_y!/BohrToAng
    dg(3) = grid_z!/BohrToAng
    ! Find grid limits
    i_spec = species_glob(i_atom)
    minx = floor( - RadiusSupport(i_spec))/dg(1)    
    maxx = floor(   RadiusSupport(i_spec))/dg(1) + 1
    miny = floor( - RadiusSupport(i_spec))/dg(2)    
    maxy = floor(   RadiusSupport(i_spec))/dg(2) + 1
    minz = floor( - RadiusSupport(i_spec))/dg(3)    
    maxz = floor(   RadiusSupport(i_spec))/dg(3) + 1
    npx = maxx - minx + 1
    npy = maxy - miny + 1
    npz = maxz - minz + 1
    allocate(psi(npx,npy,npz))
    write(*,*) 'Atom, SF: ',i_atom,MSSF_coeffs(i_atom)%n_sf
    do i_nsf = 1, MSSF_coeffs(i_atom)%n_sf
       psi = zero
       ! Loop over grid points
       do i_grid_z = minz, maxz
          ! Find z grid position and dz
          dz = dg(3)*real(i_grid_z-1,double)
          iz = i_grid_z - minz + 1
          do i_grid_y = miny, maxy
             ! Find y grid position and dy and wrap grid point
             dy = dg(2)*real(i_grid_y-1,double)
             iy = i_grid_y - miny + 1
             do i_grid_x = minx, maxx
                ! Find x grid position and dx and wrap grid point
                dx = dg(1)*real(i_grid_x-1,double)
                ix = i_grid_x - minx + 1
                ! Calculate dr
                dr = sqrt(dx*dx + dy*dy + dz*dz)
                if(dr<=RadiusSupport(i_spec)) then
                   do i_nab = 1, MSSF_coeffs(i_atom)%n_neighbours
                      dx_dr = dx - MSSF_coeffs(i_atom)%neigh_coeff(i_nab)%dx
                      dy_dr = dy - MSSF_coeffs(i_atom)%neigh_coeff(i_nab)%dy
                      dz_dr = dz - MSSF_coeffs(i_atom)%neigh_coeff(i_nab)%dz
                      dr = sqrt(dx_dr*dx_dr + dy_dr*dy_dr + dz_dr*dz_dr)
                      j_spec = species_glob(MSSF_coeffs(i_atom)%neigh_coeff(i_nab)%i_glob)
                      npao = 1
                      ! Loop over l
                      do i_l = 0, pao(j_spec)%greatest_angmom
                         ! Loop over zeta
                         do i_zeta = 1, pao(j_spec)%angmom(i_l)%n_zeta_in_angmom
                            ! Find f(r), df/dr
                            del_r = pao(j_spec)%angmom(i_l)%zeta(i_zeta)%delta
                            j = floor(dr/del_r) + 1
                            if(j+1<=pao(j_spec)%angmom(i_l)%zeta(i_zeta)%length) then
                               rr = real(j,double)*del_r
                               a = (rr - dr)/del_r
                               b = one - a
                               c = a * ( a * a - one ) * del_r * del_r / six
                               d = b * ( b * b - one ) * del_r * del_r / six
                               r1 = pao(j_spec)%angmom(i_l)%zeta(i_zeta)%table(j)
                               r2 = pao(j_spec)%angmom(i_l)%zeta(i_zeta)%table(j+1)
                               r3 = pao(j_spec)%angmom(i_l)%zeta(i_zeta)%table2(j)
                               r4 = pao(j_spec)%angmom(i_l)%zeta(i_zeta)%table2(j+1)
                               f_r = a*r1 + b*r2 + c*r3 + d*r4
                            else
                               f_r = zero
                            end if
                            ! Loop over m
                            do i_m = -i_l, i_l
                               ! Get spherical harmonic and gradient
                               ! NB returns r^l times spherical harmonic
                               call spherical_harmonic_rl(dx_dr, dy_dr, dz_dr, i_l, i_m, sph_rl, dsph_rl)
                               ! Accumulate into psi, grad psi
                               ! Coeff
                               b_ia_jb = MSSF_coeffs(i_atom)%neigh_coeff(i_nab)%coeff(npao,i_nsf)
                               psi(ix, iy, iz) = psi(ix, iy, iz) + b_ia_jb * sph_rl * f_r
                               npao = npao + 1
                            end do ! i_m
                         end do ! i_zeta
                      end do ! i_l
                   end do ! i_nab
                end if ! dr <= RadiusAtomf
             end do ! i_grid_x
          end do ! i_grid_y
       end do ! i_grid_z
       write(ci,'("Atom",I0.6,"MSSF",I0.3)') i_atom, i_nsf
       call write_mssf_cube(psi,npx, npy, npz, ci)
    end do
    return
  end subroutine write_one_MSSF

  ! Return spherical harmonic times r^l (and derivative of this)
  subroutine spherical_harmonic_rl(x, y, z, l, m, sph, dsph)
    
    use datatypes
    use numbers

    implicit none

    ! Passed variables
    integer :: l, m
    real(double) :: x, y, z, sph
    real(double), dimension(3) :: dsph

    ! Local variables
    real(double) :: r
    real(double) :: prefac

    dsph = zero
    if(l==0) then
       prefac = one/sqrt(fourpi)
       sph = prefac
       dsph = zero
    else if(l==1) then
       prefac = sqrt(3/fourpi)
       select case(m)
       case(-1) ! py
          sph = prefac*y
          dsph(2) = prefac
       case(0)  ! pz
          sph = prefac*z
          dsph(3) = prefac
       case(1)  ! px
          sph = prefac*x
          dsph(1) = prefac
       end select
    else if(l==2) then
       select case(m)
       case(-2) ! xy
          prefac = sqrt(fifteen/fourpi)
          sph = prefac*x*y
          dsph(1) = prefac*y
          dsph(2) = prefac*x
       case(-1) ! yz
          prefac = sqrt(fifteen/fourpi)
          sph = prefac*y*z
          dsph(2) = prefac*z
          dsph(3) = prefac*y
       case(0) ! 3z^2 - r^2
          prefac = sqrt(five/(sixteen*pi))
          sph = prefac*(two*z*z - x*x - y*y)
          dsph(1) = -prefac*two*x
          dsph(2) = -prefac*two*y
          dsph(3) =  prefac*four*z
       case(1) ! xz
          prefac = sqrt(fifteen/fourpi)
          sph = prefac*z*x
          dsph(1) = prefac*z
          dsph(3) = prefac*x
       case(2) ! x^2 - y^2
          prefac = sqrt(fifteen/(sixteen*pi))
          sph = prefac*(x*x - y*y)
          dsph(1) =  prefac*two*x
          dsph(2) = -prefac*two*y
       end select
    end if
  end subroutine spherical_harmonic_rl

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
