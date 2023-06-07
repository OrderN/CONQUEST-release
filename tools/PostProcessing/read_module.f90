! Utilities associated with reading input for Conquest charge/STM manipulation
module read

  implicit none

  character(len=80) :: block_file
  
contains

  ! Read Conquest_input file for parameters from simulation, output parameters and coordinates
  subroutine read_input

    use global_module, ONLY: flag_assign_blocks, flag_fractional_atomic_coords, nspin, &
         flag_wf_range_Ef, E_DOS_min, E_DOS_max, sigma_DOS, n_DOS
    use local
    use input_module
    use numbers
    use io_module, ONLY: pdb_format, pdb_template, read_atomic_positions, flag_MatrixFile_BinaryFormat
    use dimens, ONLY: r_super_x, r_super_y, r_super_z, GridCutoff
    use species_module, ONLY: n_species, species_label, species_file, mass, type_species, charge, nsf_species
    use units, ONLY: HaToeV
    use block_module, only: n_pts_in_block, in_block_x,in_block_y,in_block_z, blocks_raster, blocks_hilbert
    use pseudo_tm_info, only: setup_pseudo_info
    use GenComms,       only: cq_abort
    use pseudopotential_common, only: pseudo_type, ABINIT, OLDPS, SIESTA
    
    implicit none

    character(len=80) :: input_string, proc_coords
    integer :: i, j, n_grid_x, n_grid_y, n_grid_z
    integer :: n_kp_lines
    logical :: flag_kp_lines, flag_spin_polarisation, flag_Multisite
    real(double) :: dk
    character(len=5) :: ps_type !To find which pseudo we use
    character(len=3) :: job

    ! Load the Conquest_input files
    call load_input
    ! Now scan for parameters
    flag_spin_polarisation   = fdf_boolean('Spin.SpinPolarised', .false.)
    nspin = 1
    if(flag_spin_polarisation) nspin = 2
    ! Grid spacing
    n_grid_x   = fdf_integer('Grid.PointsAlongX',0)
    n_grid_y   = fdf_integer('Grid.PointsAlongY',0)
    n_grid_z   = fdf_integer('Grid.PointsAlongZ',0)    
    if(n_grid_x>0.AND.n_grid_y>0.AND.n_grid_z>0) then
       dk = pi/min(n_grid_x, n_grid_y, n_grid_z)
       GridCutoff = half*dk*dk
    else
       GridCutoff = fdf_double('Grid.GridCutoff',50.0_double)
    end if
    ! Block size in grid points
    in_block_x = fdf_integer('Grid.InBlockX',4)
    in_block_y = fdf_integer('Grid.InBlockY',4)
    in_block_z = fdf_integer('Grid.InBlockZ',4)
    n_pts_in_block = in_block_x*in_block_y*in_block_z
    ! Distribution of grid points
    input_string = fdf_string(10,'General.BlockAssign','Hilbert')
    if(leqi(input_string,'Raster')) then
       flag_assign_blocks = blocks_raster
       block_file = 'raster_make_blk.dat'
    else if(leqi(input_string,'Hilbert')) then
       flag_assign_blocks = blocks_hilbert
       block_file = 'hilbert_make_blk.dat'
    end if
    ! If we don't read blocks, we'll need to have a rastering routine
    !flag_read_blocks = fdf_boolean('Grid.ReadBlocksCharge',.false.)
    ! Atomic coordinates - file, format etc
    flag_fractional_atomic_coords = &
         fdf_boolean('IO.FractionalAtomicCoords',.true.)
    input_string = fdf_string(80,'IO.Coordinates',' ')
    proc_coords = fdf_string(80,'Process.Coordinates',input_string)
    pdb_format = fdf_boolean ('IO.PdbIn',.false.)
    if ( pdb_format ) then
       pdb_template = fdf_string(80,'IO.PdbTemplate',proc_coords)
    else
       pdb_template = fdf_string(80,'IO.PdbTemplate',' ')
    end if
    ! Format of wavefunction coefficient files
    flag_MatrixFile_BinaryFormat = fdf_boolean('IO.MatrixFile.BinaryFormat', .true.)
    ! Number of species
    n_species = fdf_integer('General.NumberOfSpecies',1)
    ! And read the positions
    call read_atomic_positions(trim(proc_coords))
    flag_Multisite = fdf_boolean('Basis.MultisiteSF', .false.)
    ! Find job to perform
    job = fdf_string(3,'Process.Job','pos') ! Default to converting output
    if(leqi(job,'pos').or.leqi(job,'coo')) then
       i_job = 1
       ! Allow user to specify output filename
       root_file = fdf_string(50,'Process.RootFile',trim(proc_coords))
       input_string = fdf_string(3,'Process.CoordFormat','xyz')
       if(leqi(input_string,'xyz')) then
          coord_format=1
       else if(leqi(input_string,'cel')) then
          coord_format=2
       else if(leqi(input_string,'xsf')) then
          coord_format=3
          flag_write_forces = fdf_boolean('Process.WriteForces',.false.)
          flag_write_spin_moments = fdf_boolean('Process.WriteSpinMoments',.false.)
          if(flag_write_forces .and. flag_write_spin_moments) call cq_abort("Cannot have both forces and spin moments output")
       else
          call cq_abort("Unrecognised output format: "//trim(input_string))
       end if
    else if(leqi(job,'chg').or.leqi(job,'cha').or.leqi(job,'den')) then
       i_job = 2
    else if(leqi(job,'ban')) then
       i_job = 3
       if(flag_Multisite) call cq_abort("Not yet compatible with multi-site support functions")
    else if(leqi(job,'ter').or.leqi(job,'th')) then
       i_job = 4
       if(flag_Multisite) call cq_abort("Not yet compatible with multi-site support functions")
       ! Allow user to specify output filename
       root_file = fdf_string(50,'Process.RootFile','STM')
    else if(leqi(job,'stm')) then
       i_job = 5
       if(flag_Multisite) call cq_abort("Not yet compatible with multi-site support functions")
       ! Allow user to specify output filename
       root_file = fdf_string(50,'Process.RootFile','STM')
    else if(leqi(job,'dos')) then
       i_job = 6
    else if(leqi(job,'pdo').or.leqi(job,'pro')) then
       i_job = 7
    end if
    ! 
    charge_stub = fdf_string(80,'Process.ChargeStub','chden')
    ! STM parameters
    ! NB Bias will be in volts
    stm_bias = fdf_double('STM.BiasVoltage',zero)
    ! Allow for Fermi level offset, also in volts
    fermi_offset = fdf_double('STM.FermiOffset',zero)
    ! Restrict the area that we consider
    stm_x_min = zero
    stm_y_min = zero
    stm_x_max = r_super_x
    stm_y_max = r_super_y
    if(i_job==4.or.i_job==5) then
       stm_z_min = fdf_double('Process.MinZ',zero)
       stm_z_max = fdf_double('Process.MaxZ',r_super_z)
    else
       stm_z_min = zero
       stm_z_max = r_super_z
    end if
    !stm_x_min = fdf_double('Process.MinX',zero)
    !stm_x_max = fdf_double('Process.MaxX',r_super_x)
    !stm_y_min = fdf_double('Process.MinY',zero)
    !stm_y_max = fdf_double('Process.MaxY',r_super_y)
    ! Broadening of the energy levels - in Ha
    stm_broad = fdf_double('STM.Sigma',0.001_double)
    stm_broad = stm_broad*HaToeV
    ! Read number of k-points for consistency
    nkp = fdf_integer('Diag.NumKpts',1)
    flag_kp_lines = fdf_boolean('Diag.KspaceLines',.false.)
    if(flag_kp_lines) then
       n_kp_lines = fdf_integer('Diag.NumKptLines',1)
       nkp = nkp*n_kp_lines
    end if
    ! Smearing parameters from Conquest
    kT = fdf_double('Diag.kT',0.001_double)
    ! Method to approximate step function for occupation number
    flag_smear_type = fdf_integer('Diag.SmearingType',0)
    iMethfessel_Paxton = fdf_integer('Diag.MPOrder',0)
    ! Repeat the cell
    nrptx = fdf_integer('Process.RepeatX',1)
    nrpty = fdf_integer('Process.RepeatY',1)
    nrptz = fdf_integer('Process.RepeatZ',1)
    nsampx = fdf_integer('Process.SampleX',1)
    nsampy = fdf_integer('Process.SampleY',1)
    nsampz = fdf_integer('Process.SampleZ',1)
    ! Energy limits (relative to Ef or absolute)
    E_wf_min = fdf_double('IO.min_wf_E',zero)
    E_wf_max = fdf_double('IO.max_wf_E',zero)
    if(i_job==3 .or. i_job==4 .or. i_job==5) then ! Band-resolved charge or STM
       ! Read in details of bands output from Conquest
       n_bands_active=fdf_integer('IO.maxnoWF',0)
       if(n_bands_active>0) then
          allocate(band_no(n_bands_active))
          if (fdf_block('WaveFunctionsOut')) then
             if(1+block_end-block_start<n_bands_active) then
                write(*,*) "Too few wf no in WaveFunctionsOut: ",1+block_end-block_start,n_bands_active
                stop
             end if
             do i=1,n_bands_active
                read(unit=input_array(block_start+i-1),fmt=*) band_no(i)
             end do
             call fdf_endblock
          else
             call cq_abort("You specified bands with IO.maxnoWF but didn't give the WaveFunctionsOut block")
          end if
          flag_wf_range = .false.
       else
          if(abs(E_wf_max-E_wf_min)>1e-8_double) then
             flag_wf_range = .true.
             flag_wf_range_Ef = fdf_boolean('IO.WFRangeRelative',.true.)
          else
             call cq_abort("No bands specified!")
          end if
       end if
       ! Now read details of bands to output from processing
       ! Energy limits
       E_procwf_min = fdf_double('Process.min_wf_E',E_wf_min)
       E_procwf_max = fdf_double('Process.max_wf_E',E_wf_max)
       ! Is the range relative to Ef (T) or absolute (F)
       flag_procwf_range_Ef = fdf_boolean('Process.WFRangeRelative',.true.)
       if(abs(E_procwf_max-E_procwf_min)>1e-8_double) then
          flag_proc_range = .true.
       else
          n_bands_process = fdf_integer('Process.noWF',0)
          if(n_bands_process>0) then
             allocate(band_proc_no(n_bands_process))
             if (fdf_block('WaveFunctionsProcess')) then
                if(1+block_end-block_start<n_bands_process) then
                   write(*,*) "Too few wf no in WaveFunctionsOut: ",1+block_end-block_start,n_bands_process
                   stop
                end if
                do i=1,n_bands_process
                   read(unit=input_array(block_start+i-1),fmt=*) band_proc_no(i)
                end do
                call fdf_endblock
             end if
             flag_proc_range = .false.
          end if
       end if
    end if ! i_job == 3 or 4 or 5
    ! Output format
    ps_type = fdf_string(4,'Process.OutputFormat','cube')
    if(leqi(ps_type,'cube')) then
       flag_output = cube
    !else if(leqi(ps_type,'dx')) then
    !   flag_output = dx
    else
       write(*,*) "Output format ",ps_type," not recognised.  Only cube available at present."
       stop
    end if
    flag_by_kpoint = fdf_boolean('Process.outputWF_by_kpoint',.false.)
    ! DOS
    ! Add flag for window relative to Fermi level
    E_DOS_min = fdf_double('Process.min_DOS_E',E_wf_min)
    E_DOS_max = fdf_double('Process.max_DOS_E',E_wf_max)
    sigma_DOS = fdf_double('Process.sigma_DOS',zero) ! Adjust to minimum of 4*energy spacing
    n_DOS = fdf_integer('Process.n_DOS',1001)
    flag_total_iDOS = fdf_boolean('Process.TotalIntegratedDOS',.false.)
    if(i_job==7) then
       ! If no limits specified, cover whole range
       if(abs(E_wf_max-E_wf_min)<1e-8_double) then
          E_wf_min = -BIG
          E_wf_max =  BIG
       end if
       flag_wf_range = .true.
       flag_procwf_range_Ef = fdf_boolean('Process.WFRangeRelative',.false.)
       flag_l_resolved = fdf_boolean('Process.pDOS_l_resolved',.false.)
       flag_lm_resolved = fdf_boolean('Process.pDOS_lm_resolved',.false.)
       if(flag_lm_resolved .and. (.not.flag_l_resolved)) flag_l_resolved = .true.
    end if
    ! Now read PS files for atomic information
    call allocate_species_vars
    ps_type = fdf_string(5,'General.PseudopotentialType','haman') 
    if(leqi(ps_type,'siest')) then
       pseudo_type = SIESTA
    else if(leqi(ps_type,'plato').OR.leqi(ps_type,'haman')) then
       pseudo_type = ABINIT
    else
       pseudo_type = OLDPS
    endif
    if(fdf_block('ChemicalSpeciesLabel')) then
       if(1+block_end-block_start<n_species) & 
            call cq_abort("Too few species in ChemicalSpeciesLabel: ",&
            1+block_end-block_start,n_species)
       do i=1,n_species
          read (unit=input_array(block_start+i-1),fmt=*) &
               j, mass(j),       &
               species_label(j)
          type_species(j)=j
       end do
       call fdf_endblock
    end if
    do i=1,n_species
       if(fdf_block(species_label(i))) then
          charge(i)        = fdf_double ('Atom.ValenceCharge',zero)
          nsf_species(i)   = fdf_integer('Atom.NumberOfSupports',0)
          call fdf_endblock
       end if
    end do
    call setup_pseudo_info ! Get atomic number of species
    call io_close(fdf_out)
    return
  end subroutine read_input

  ! We read in the block positions
  subroutine read_block_input

    use datatypes
    use local, ONLY: nblockx, nblocky, nblockz, block_store, nprocs, block_size_x, block_size_y, block_size_z, &
         grid_x, grid_y, grid_z, nptsx, nptsy, nptsz, stm_z_min, stm_z_max, stm_x_min, stm_x_max, &
         stm_y_min, stm_y_max, nxmin, nymin, nzmin, gpv
    use dimens, ONLY: r_super_x, r_super_y, r_super_z, volume
    use block_module, only: n_pts_in_block, in_block_x,in_block_y,in_block_z
    !use units, ONLY: BohrToAng
    
    implicit none

    integer :: proc, idum, idum2, iblock, blocks_on_proc, ind_group
    real(double) :: BohrToAng

    BohrToAng = 1.0_double

    write(*,fmt='(2x,"Opening block file: ",a)') block_file
    open(unit=17,file=block_file)
    ! Numbers of blocks in x/y/z
    read(17,*) nblockx,nblocky,nblockz
    write(*,fmt='(4x,"Blocks in x, y, z: ",3i5)') nblockx,nblocky,nblockz
    ! Physical size of blocks in bohr
    block_size_x = BohrToAng*r_super_x/real(nblockx,double)
    block_size_y = BohrToAng*r_super_y/real(nblocky,double)
    block_size_z = BohrToAng*r_super_z/real(nblockz,double)
    write(*,'(4x,"Block sizes in bohr: ",3f12.5)') block_size_x, block_size_y, block_size_z
    ! Grid spacing in bohr
    grid_x = BohrToAng*r_super_x/real(nblockx*in_block_x,double)
    grid_y = BohrToAng*r_super_y/real(nblocky*in_block_y,double)
    grid_z = BohrToAng*r_super_z/real(nblockz*in_block_z,double)
    ! Cell volume
    volume = BohrToAng*BohrToAng*BohrToAng*r_super_x*r_super_y*r_super_z
    ! Number of points in the area that the user has specified
    nptsx = floor(BohrToAng*(stm_x_max - stm_x_min)/grid_x)
    nptsy = floor(BohrToAng*(stm_y_max - stm_y_min)/grid_y)
    nptsz = floor(BohrToAng*(stm_z_max - stm_z_min)/grid_z)
    write(*,fmt='(4x,"Points in user-specified area: ",3i5)') nptsx,nptsy,nptsz
    gpv = volume/real(nptsx*nptsy*nptsz,double)
    ! Work out starting grid points in user area
    nxmin = floor(BohrToAng*stm_x_min/grid_x)
    if(stm_x_min - grid_x*real(nxmin,double)>1e-3) nxmin = nxmin + 1
    nymin = floor(BohrToAng*stm_y_min/grid_y)
    if(stm_y_min - grid_y*real(nymin,double)>1e-3) nymin = nymin + 1
    nzmin = floor(BohrToAng*stm_z_min/grid_z)
    if(stm_z_min - grid_z*real(nzmin,double)>1e-3) nzmin = nzmin + 1
    read(17,*) nprocs
    write(*,fmt='(2x,"Original run on ",i8," processors")') nprocs
    ! Allocate array of types
    allocate(block_store(nprocs))
    ! Loop over processors
    do proc = 1,nprocs
       ! Read from make_blk
       read(17,*) idum,blocks_on_proc,idum2
       if(idum/=proc) write(*,*) 'Proc mismatch in block file: ',idum,proc
       allocate(block_store(proc)%nx(blocks_on_proc),block_store(proc)%ny(blocks_on_proc), &
            block_store(proc)%nz(blocks_on_proc),block_store(proc)%num(blocks_on_proc), &
            block_store(proc)%active(blocks_on_proc))
       block_store(proc)%process = idum
       block_store(proc)%num_blocks = blocks_on_proc
       ! Allocate num, x, y, z for proc
       do iblock=1, blocks_on_proc
          read(17,*) idum,ind_group
          block_store(proc)%num(iblock) = ind_group
          block_store(proc)%nx(iblock) = 1+(ind_group-1)/(nblocky*nblockz)
          block_store(proc)%ny(iblock) = 1+(ind_group-1-(block_store(proc)%nx(iblock)-1)*nblocky*nblockz)/nblockz
          block_store(proc)%nz(iblock) = ind_group - (block_store(proc)%nx(iblock)-1)*nblocky*nblockz - &
               (block_store(proc)%ny(iblock)-1)*nblockz
       end do
    end do
    close(unit=17)
    return
  end subroutine read_block_input

  ! We read in the block positions
  subroutine read_nprocs_from_blocks

    use datatypes
    use local, ONLY: nprocs
    
    implicit none

    integer :: proc, idum, idum2, iblock, blocks_on_proc, ind_group, nblockx,nblocky,nblockz

    write(*,fmt='(2x,"Opening block file: ",a)') block_file
    open(unit=17,file=block_file)
    ! Numbers of blocks in x/y/z
    read(17,*) nblockx,nblocky,nblockz
    read(17,*) nprocs
    write(*,fmt='(2x,"Original run on ",i8," processors")') nprocs
    close(unit=17)
    return
  end subroutine read_nprocs_from_blocks
  
  subroutine read_eigenvalues

    use datatypes
    use numbers
    use local, ONLY: nkp, n_eval_window, efermi, kx, ky, kz, wtk, efermi, stm_bias, stm_broad, &
         n_bands_active, eigenvalues, band_no, flag_only_charge, flag_by_kpoint, E_wf_min, E_wf_max, &
         flag_wf_range, n_bands_total, band_active_kp, band_active_all, n_bands_process, band_proc_no, &
         band_full_to_active, i_job
    use units, ONLY: HaToeV
    use global_module, only: flag_wf_range_Ef, nspin
    
    implicit none

    ! Local variables
    integer :: idum, n_evals,i_kp,i_eval, i, i_spin, i_low, i_high
    real(double), dimension(2) :: Emin, Emax
    real(double) :: eval, Eband
    integer, dimension(:), allocatable :: active_bands
    real(double), dimension(:,:),allocatable :: tmp_evals ! Read and store
    character(len=80) :: str, str2
    
    open(unit=17,file='eigenvalues.dat')
    read(17,*) str,n_evals,str2,idum
    n_bands_total = n_evals
    if(i_job==7 .and. n_bands_active==0) then
       n_bands_active = n_bands_total
    end if
    if(idum/=nkp) then
       write(*,fmt='(4x,"Reading k-points from eigenvalues file, not setting from input file")')
       nkp = idum
    end if
    ! Fermi level
    efermi = zero
    if(nspin==1) then
       read(17,fmt='(a6,f18.10)') str,efermi(1)
       write(*,fmt='(4x,"Fermi level: ",f12.5," Ha")') efermi(1)
    else
       read(17,fmt='(a6,2f18.10)') str,efermi(1), efermi(2)
       write(*,fmt='(4x,"Fermi levels: ",2f12.5," Ha")') efermi
    end if
    read(17,*) str
    ! Allocate memory
    allocate(kx(nkp),ky(nkp),kz(nkp),wtk(nkp))
    allocate(eigenvalues(n_evals,nkp,nspin), band_active_kp(n_evals,nkp,nspin), band_active_all(n_evals,nspin))
    allocate(band_full_to_active(n_evals))
    band_full_to_active = 0
    eigenvalues = zero
    band_active_kp = 0
    band_active_all = 0
    ! Active bands lie within the energy window specified for STM images
    Emin = E_wf_min
    Emax = E_wf_max
    if(flag_wf_range_Ef) then
       Emin = efermi + Emin
       Emax = efermi + Emax
    end if
    if(flag_wf_range) then
       if(nspin==1) then
          write(*,fmt='(2x,"Reading bands between ",e12.3,"Ha and ",e12.3,"Ha")') Emin(1), Emax(1)
       else
          write(*,fmt='(2x,"SpinU reading bands between ",e12.3,"Ha and ",e12.3,"Ha")') Emin(1), Emax(1)
          write(*,fmt='(2x,"SpinD reading bands between ",e12.3,"Ha and ",e12.3,"Ha")') Emin(2), Emax(2)
       end if
    end if
    ! Now loop over k-points and read eigenvalues
    i_low = n_evals+1
    i_high = 0
    do i_spin = 1, nspin
       do i_kp = 1,nkp
          read(17,*) idum,kx(i_kp),ky(i_kp),kz(i_kp),wtk(i_kp)
          do i_eval = 1,n_evals
             read(17,*) idum, eval
             eigenvalues(i_eval,i_kp,i_spin) = eval
             if(flag_wf_range) then
                if(eval>=Emin(i_spin).AND.eval<=Emax(i_spin)) then
                   band_active_kp(i_eval, i_kp, i_spin) = 1
                   band_active_all(i_eval, i_spin) = 1
                   if(i_eval<i_low) i_low = i_eval
                   if(i_eval>i_high) i_high = i_eval
                end if
             end if
          end do
       end do
       if(flag_wf_range) then
          write(*,fmt='(2x,"Bands ",I5," to ",I5," were written out by CONQUEST")') i_low, i_high
          n_bands_active = i_high - i_low + 1
       else
          do i_eval=1,n_bands_active
             band_active_kp(band_no(i_eval),:,i_spin) = 1
             band_active_all(band_no(i_eval),i_spin) = 1
             band_full_to_active(band_no(i_eval)) = i_eval
          end do
          write(*,fmt='(2x,"Active bands: ",(I5))') band_no
       end if
    end do
    if(flag_wf_range) then
       allocate(band_no(n_bands_active))
       do i_eval = 1, n_bands_active
          band_no(i_eval) = i_low + i_eval - 1
          band_full_to_active(i_low+i_eval-1) = i_eval
       end do
    end if
    close(unit=17)
    !n_bands_active = sum(band_active_all)
  end subroutine read_eigenvalues

  ! Here we will read ALL the bands written out by Conquest (n_bands_active)
  subroutine read_psi_coeffs(stub)

    use datatypes
    use numbers
    use species_module, ONLY: nsf_species
    use global_module, ONLY: ni_in_cell, species_glob, nspin
    use local, ONLY: nkp, nprocs, n_bands_active, band_active_kp, eigenvalues, evec_coeff, n_bands_active, band_no
    use io_module, ONLY: flag_MatrixFile_BinaryFormat

    implicit none

    ! Passed variables
    character(len=*) :: stub

    ! Local variables
    integer :: i_sf, i_proc, i_kp, i, i_band, i_prim, n_prim, i_glob, i_atom, i_spin
    real(double) :: eval
    complex(double_cplx) :: coeff
    character(len=80) :: filename, str

    ! Allocate space
    allocate(evec_coeff(maxval(nsf_species), ni_in_cell, n_bands_active, nkp, nspin))
    evec_coeff = zero
    if(flag_MatrixFile_BinaryFormat) then
       ! Read coefficients
       do i_spin = 1, nspin
          do i_proc = 1, nprocs
             if(nspin==1) then
                write(filename,'(a,I0.7,"WF.dat")') trim(stub),i_proc
             else
                write(filename,'(a,I0.7,"WF",I0.1,".dat")') trim(stub),i_proc, i_spin
             end if
             open(unit=17,file=filename,form='unformatted')
             do i_kp = 1, nkp
                read(17) n_prim ! Number of atoms on process
                do i_band = 1, n_bands_active
                   if(band_active_kp(band_no(i_band), i_kp, i_spin) == 1) then
                      read(17) i, eval
                      ! Loop over primary atoms
                      do i_prim = 1, n_prim
                         read(17) i_atom ! Global number of atom
                         ! Loop over SFs
                         do i_sf = 1, nsf_species(species_glob(i_atom))
                            read(17) evec_coeff(i_sf, i_atom, i_band, i_kp, i_spin)
                         end do ! nsf
                      end do ! n_prim primary atoms
                   end if
                end do ! bands
             end do ! nkp kpoints
             close(unit=17)
          end do! nprocs processes
       end do ! i_spin = nspin
    else
       ! Read coefficients
       do i_spin = 1, nspin
          do i_proc = 1, nprocs
             if(nspin==1) then
                write(filename,'(a,I0.7,"WF.dat")') trim(stub),i_proc
             else
                write(filename,'(a,I0.7,"WF",I0.1,".dat")') trim(stub),i_proc, i_spin
             end if
             open(unit=17,file=filename)
             do i_kp = 1, nkp
                read(17,*) n_prim ! Number of atoms on process
                do i_band = 1, n_bands_active
                   if(band_active_kp(band_no(i_band), i_kp, i_spin) == 1) then
                      read(17,*) i, eval
                      ! Loop over primary atoms
                      do i_prim = 1, n_prim
                         read(17,*) i_atom ! Global number of atom
                         ! Loop over SFs
                         do i_sf = 1, nsf_species(species_glob(i_atom))
                            read(17,*) evec_coeff(i_sf, i_atom, i_band, i_kp, i_spin)
                         end do ! nsf
                      end do ! n_prim primary atoms
                   end if
                end do ! bands
             end do ! nkp kpoints
             close(unit=17)
          end do! nprocs processes
       end do ! i_spin = nspin
    end if ! Binary format
  end subroutine read_psi_coeffs
 
  subroutine allocate_species_vars

    use numbers
    use dimens,         only: atomicnum, RadiusSupport, RadiusAtomf, InvSRange
    use memory_module,  only: reg_alloc_mem, type_dbl
    use species_module, only: charge, charge_up, charge_dn
    use species_module, only: mass, npao_species, type_species, nsf_species
    use species_module, only: species_label, species_file, n_species
    use GenComms,       only: cq_abort

    implicit none

    ! Local variables
    integer        :: stat

    allocate(npao_species(n_species),STAT=stat)
    allocate(nsf_species(n_species),STAT=stat)
    nsf_species = zero
    allocate(type_species(n_species),STAT=stat)
    allocate(atomicnum(n_species),STAT=stat)
    allocate(RadiusSupport(n_species),STAT=stat)
    RadiusSupport = zero
    allocate(InvSRange(n_species),STAT=stat)
    allocate(RadiusAtomf(n_species),STAT=stat)
    RadiusAtomf = zero
    if(stat/=0) call cq_abort("Error allocating atomicnum in allocate_species_vars: ",n_species,stat)
    allocate(charge(n_species),STAT=stat)
    allocate(charge_up(n_species),STAT=stat)
    allocate(charge_dn(n_species),STAT=stat)
    charge_up = zero
    charge_dn = zero
    if(stat/=0) call cq_abort("Error allocating charge in allocate_species_vars: ",                  n_species,stat)
    allocate(mass(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating mass in allocate_species_vars: ",                    n_species,stat)
    allocate(species_label(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating species_label in allocate_species_vars: ",           n_species,stat)
    allocate(species_file(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating species_file in allocate_species_vars: ",            n_species,stat)
    return
  end subroutine allocate_species_vars
  
end module read
