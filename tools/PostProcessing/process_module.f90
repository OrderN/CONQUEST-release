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
         n_bands_total, band_active_kp, flag_proc_range, &
         E_procwf_min, E_procwf_max, flag_procwf_range_Ef, band_proc_no, n_bands_process
    use output, ONLY: write_cube
    use global_module, only : nspin
    use read, ONLY: read_eigenvalues, read_psi_coeffs

    implicit none

    integer :: proc, band, nk, idum1, idum2, kp, ispin
    real(double) :: weight, rbx, rby, rbz, sq, test
    real(double), dimension(2) :: Emin, Emax
    character(len=50) :: filename, ci
    complex(double_cplx), dimension(:,:,:), allocatable :: psi

    write(*,*) 'Processing bands'
    ! Read eigenvalues
    call read_eigenvalues
    ! Read eigenvector coefficients
    call read_psi_coeffs
    allocate(current(nptsx,nptsy,nptsz))
    allocate(psi(nptsx,nptsy,nptsz))
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
                        eigenvalues(band,kp,ispin) <= Emax(ispin)) then
                      call pao_to_grid(band, kp, ispin, psi)
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
                        eigenvalues(band,kp,ispin) <= Emax(ispin)) then
                      call pao_to_grid(band, kp, ispin, psi)
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
                      call pao_to_grid(band_proc_no(band), kp, ispin, psi)
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
                      call pao_to_grid(band_proc_no(band), kp, ispin, psi)
                      current = current + psi*conjg(psi)
                   end if
                end do ! kp
                write(ci,'("Band",I0.6,"den_totS",I0.1)') band_proc_no(band), ispin
                call write_cube(current,ci)
             end do ! bands
          end if
       end do
    end if
    return
  end subroutine process_bands

  subroutine process_dos

    use datatypes
    use numbers, ONLY: zero, RD_ERR, twopi, half, one
    use local, ONLY: eigenvalues, n_bands_total, nkp, wtk, efermi
    use read, ONLY: read_eigenvalues, read_psi_coeffs
    use global_module, ONLY: nspin, n_DOS, E_DOS_min, E_DOS_max, sigma_DOS
    use units, ONLY: HaToeV

    implicit none
    
    ! Local variables
    integer :: i_band, i_kp, i_spin, n_DOS_wid, n_band, n_min, n_max, i
    real(double) :: Ebin, dE_DOS, a, pf_DOS
    real(double), dimension(:,:), allocatable :: total_DOS
    
    ! Read eigenvalues
    call read_eigenvalues
    allocate(total_DOS(n_DOS,nspin))
    total_DOS = zero
    ! Set limits on DOS output
    if(abs(E_DOS_min)<RD_ERR .and. abs(E_DOS_max)<RD_ERR) then
       E_DOS_min = 1e30_double
       E_DOS_max = -1e30_double
       do i=1,nkp
          if(E_DOS_min>eigenvalues(1,i,1)) E_DOS_min = eigenvalues(1,i,1)
          if(E_DOS_max<eigenvalues(n_bands_total,i,1)) E_DOS_max = eigenvalues(n_bands_total,i,1)
       end do
       write(*,fmt='(2x,"DOS limits set automatically: ",2f12.5)') E_DOS_min, E_DOS_max
    end if
    ! Spacing, width, prefactor
    dE_DOS = (E_DOS_max - E_DOS_min)/real(n_DOS-1,double)
    n_DOS_wid = floor(6.0_double*sigma_DOS/dE_DOS) ! How many bins either side of state we consider
    pf_DOS = one/(sigma_DOS*sqrt(twopi))
    ! Accumulate DOS over bands and k-points for each spin
    do i_spin = 1, nspin
       do i_kp = 1, nkp
          do i_band=1,n_bands_total ! All bands
             n_band = floor((eigenvalues(i_band, i_kp, i_spin) - E_DOS_min)/dE_DOS) + 1
             n_min = n_band - n_DOS_wid
             if(n_min<1) n_min = 1
             n_max = n_band + n_DOS_wid
             if(n_max>n_DOS) n_max = n_DOS
             do i = n_min, n_max
                Ebin = real(i-1,double)*dE_DOS + E_DOS_min
                a = (Ebin-eigenvalues(i_band, i_kp, i_spin))/sigma_DOS
                total_DOS(i,i_spin) = total_DOS(i,i_spin) + wtk(i_kp)*pf_DOS*exp(-half*a*a)
             end do
          end do
       end do
    end do
    ! Write out DOS, shifted to Ef = 0
    open(unit=17, file="DOS.dat")
    do i_spin = 1, nspin
       write(17,fmt='("# Spin ",I1)') i_spin
       write(17,fmt='("# Original Fermi-level: ",f12.5," eV")') HaToeV*efermi(i_spin)
       write(17,fmt='("# DOS shifted relative to Fermi-level")')
       do i=1, n_DOS
          write(17,fmt='(2f12.5)') HaToeV*(E_DOS_min + dE_DOS*real(i-1,double)-efermi(i_spin)), total_DOS(i,i_spin)
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
  
end module process
