! Utilities associated with reading input for Conquest charge/STM manipulation
module read

  implicit none

  character(len=80) :: block_file
  
contains

  ! Read Conquest_input file for parameters from simulation, output parameters and coordinates
  subroutine read_input

    use global_module, ONLY: flag_read_blocks, flag_assign_blocks, flag_fractional_atomic_coords, nspin, &
         flag_wf_range_Ef
    use local
    use input_module
    use numbers
    use io_module, ONLY: pdb_format, pdb_template, read_atomic_positions
    use dimens, ONLY: r_super_x, r_super_y, r_super_z
    use species_module, ONLY: n_species, species_label, species_file, mass, type_species, charge, nsf_species
    use units, ONLY: HaToeV
    use block_module, only: n_pts_in_block, in_block_x,in_block_y,in_block_z, blocks_raster, blocks_hilbert
    use pseudo_tm_info, only: setup_pseudo_info
    use GenComms,       only: cq_abort
    use pseudopotential_common, only: pseudo_type, ABINIT, OLDPS, SIESTA
    
    implicit none

    character(len=80) :: input_string
    integer :: i, j
    integer :: n_kp_lines
    logical :: flag_kp_lines, flag_spin_polarisation
    character(len=5)  :: ps_type !To find which pseudo we use

    ! Load the Conquest_input files
    call load_input
    ! Now scan for parameters
    flag_spin_polarisation   = fdf_boolean('Spin.SpinPolarised', .false.)
    nspin = 1
    if(flag_spin_polarisation) nspin = 2
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
    flag_read_blocks = fdf_boolean('Grid.ReadBlocksCharge',.false.)
    ! Atomic coordinates - file, format etc
    flag_fractional_atomic_coords = &
         fdf_boolean('IO.FractionalAtomicCoords',.true.)
    input_string = fdf_string(80,'IO.Coordinates',' ')
    pdb_format = fdf_boolean ('IO.PdbIn',.false.)
    if ( pdb_format ) then
       pdb_template = fdf_string(80,'IO.PdbTemplate',input_string)
    else
       pdb_template = fdf_string(80,'IO.PdbTemplate',' ')
    end if
    n_species = fdf_integer('General.NumberOfSpecies',1)
    ! And read the positions
    call read_atomic_positions(trim(input_string))
    flag_only_charge = fdf_boolean('Process.OnlyCharge',.true.)
    if(flag_only_charge) charge_stub = fdf_string(80,'Process.ChargeStub','chden')
    ! STM parameters
    ! NB Bias will be in volts
    stm_bias = fdf_double('STM.BiasVoltage',zero)
    ! Allow for Fermi level offset, also in volts
    fermi_offset = fdf_double('STM.FermiOffset',zero)
    ! Restrict the area that we consider
    stm_z_min = fdf_double('Process.MinZ',zero)
    stm_z_max = fdf_double('Process.MaxZ',r_super_z)
    stm_x_min = fdf_double('Process.MinX',zero)
    stm_x_max = fdf_double('Process.MaxX',r_super_x)
    stm_y_min = fdf_double('Process.MinY',zero)
    stm_y_max = fdf_double('Process.MaxY',r_super_y)
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
    ! Allow user to specify output filename
    root_file = fdf_string(50,'Process.RootFile','STM')
    ! Repeat the cell
    nrptx = fdf_integer('Process.RepeatX',1)
    nrpty = fdf_integer('Process.RepeatY',1)
    nrptz = fdf_integer('Process.RepeatZ',1)
    nsampx = fdf_integer('Process.SampleX',1)
    nsampy = fdf_integer('Process.SampleY',1)
    nsampz = fdf_integer('Process.SampleZ',1)
    if(flag_only_charge) then
       ! Read in details of band output
       flag_by_kpoint = fdf_boolean('IO.outputWF_by_kpoint',.false.)
       ! Have a specific number of bands been output ?
       n_bands_active=fdf_integer('IO.maxnoWF',0)
       !if(n_bands_active==0) then
       !   write(*,*) 'Please specify band numbers you want converted using'
       !   write(*,*) 'IO.maxnoWF for number of bands and %block WaveFunctionsOut'
       !   write(*,*) 'to specify bands'
       !   stop
       !end if
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
          end if
       end if
       flag_range = .false.
       E_wf_min = fdf_double('IO.min_wf_E',zero)
       E_wf_max = fdf_double('IO.max_wf_E',zero)
       if(abs(E_wf_max-E_wf_min)>1e-8_double) flag_range = .true.
       ! Is the range relative to Ef (T) or absolute (F)
       flag_wf_range_Ef = fdf_boolean('IO.WFRangeRelative',.true.)
    end if
    ps_type = fdf_string(4,'Process.OutputFormat','cube')
    if(leqi(ps_type,'cube')) then
       flag_output = cube
    else if(leqi(ps_type,'dx')) then
       flag_output = dx
    else
       write(*,*) "Output format ",root_file," not recognised.  Choose cube or dx."
       stop
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
          write(*,*) 'Species ',j,species_label(j)
       end do
       call fdf_endblock
    end if
    do i=1,n_species
       if(fdf_block(species_label(i))) then
          charge(i)        = fdf_double ('Atom.ValenceCharge',zero)
          nsf_species(i)   = fdf_integer('Atom.NumberOfSupports',0)
          call fdf_endblock
       else
          write(*,*) "Can't find species block: ",species_label(i)
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

  subroutine read_eigenvalues_bands

    use datatypes
    use numbers
    use local, ONLY: nkp, n_eval_window, efermi, kx, ky, kz, wtk, efermi, stm_bias, stm_broad, &
         n_bands_active, eigenvalues, band_no, flag_only_charge, flag_by_kpoint, E_wf_min, E_wf_max, &
         flag_range
    use units, ONLY: HaToeV
    use global_module, only: flag_wf_range_Ef
    
    implicit none

    integer :: idum, n_evals,nk,nev, i
    real(double) :: Emin, Emax, eval, Eband
    integer, dimension(:), allocatable :: active_bands
    real(double), dimension(:,:),allocatable :: tmp_evals ! Read and store
    character(len=80) :: str, str2
    
    open(unit=17,file='eigenvalues.dat')
    read(17,*) str,n_evals,str2,idum
    if(idum/=nkp) then
       write(*,fmt='("Reading k-points from eigenvalues file, not setting from input file")')
       nkp = idum
    end if
    ! Fermi level
    read(17,fmt='(a6,f12.5)') str,efermi
    write(*,fmt='(4x,"Fermi level: ",f12.5," Ha")') efermi
    read(17,*) str
    ! Allocate memory
    allocate(kx(nkp),ky(nkp),kz(nkp),wtk(nkp))
    allocate(tmp_evals(n_evals,nkp))
    ! Active bands lie within the energy window specified for STM images
    n_bands_active = 0
    Emin = E_wf_min
    Emax = E_wf_max
    if(flag_wf_range_Ef) then
       Emin = efermi + Emin
       Emax = efermi + Emax
    end if
    write(*,fmt='(2x,"Searching for bands between ",f9.3,"Ha and ",f9.3,"Ha")') Emin, Emax
    ! Now loop over k-points and read eigenvalues
    do nk = 1,nkp
       read(17,*) idum,kx(nk),ky(nk),kz(nk),wtk(nk)
       do nev = 1,n_evals
          read(17,*) idum,eval
          tmp_evals(nev,nk) = eval
          if(eval>Emin.AND.eval<Emax) n_bands_active = n_bands_active+1
       end do
    end do
    close(unit=17)
    allocate(band_no(n_bands_active))
    i = 1
    do nk = 1,nkp
       do nev = 1,n_evals
          if(tmp_evals(nev,nk)>Emin.AND.tmp_evals(nev,nk)<Emax) then
             band_no(i) = nev
             i = i+1
          end if
       end do
    end do
    deallocate(tmp_evals)
  end subroutine read_eigenvalues_bands

  ! For STM output, we need to read in eigenvalues
  ! Added band-by-band allocation for charge only 2016/04/19 08:10 dave
  ! Changing to STM version (see band version above)
  subroutine read_eigenvalues_stm

    use datatypes
    use numbers
    use local, ONLY: nkp, n_eval_window, efermi, kx, ky, kz, wtk, efermi, stm_bias, stm_broad, &
         n_bands_active, eigenvalues, band_no, flag_only_charge, flag_by_kpoint, E_wf_min, E_wf_max, &
         flag_range, fermi_offset
    use units, ONLY: HaToeV
    use global_module, only: flag_wf_range_Ef
    
    implicit none

    integer :: idum, n_evals,nk,nev
    real(double) :: Emin, Emax, eval, Eband
    integer, dimension(:), allocatable :: active_bands
    real(double), dimension(:,:),allocatable :: tmp_evals ! Read and store
    character(len=80) :: str, str2
    
    open(unit=17,file='eigenvalues.dat')
    read(17,*) str,n_evals,str2,idum
    if(idum/=nkp) then
       write(*,fmt='("Warning ! Kpoints from input incompatible with eigenvalues: ",2i5)') idum,nkp
       write(*,fmt='("This is not a problem if you set k-points automatically.")')
       nkp = idum
    end if
    ! Fermi level
    read(17,fmt='(a6,f12.5)') str,efermi
    efermi = efermi*HaToeV
    write(*,fmt='(4x,"Fermi level: ",f12.5," eV")') efermi
    if(abs(fermi_offset)>1e-8_double) then
       efermi = efermi + fermi_offset
       write(*,fmt='(4x,"After applying offset, Fermi level: ",f12.5," eV")') efermi
    end if
    read(17,*) str
    ! Allocate memory
    allocate(kx(nkp),ky(nkp),kz(nkp),wtk(nkp))
    allocate(tmp_evals(n_evals,nkp))
    ! Active bands lie within the energy window specified for STM images
    allocate(active_bands(n_evals))
    active_bands = 0
    ! Allow for thermal broadening
    if(stm_bias<0) then
       Emax = efermi + two*stm_broad
       Emin = efermi + stm_bias - two*stm_broad
    else
       Emin = efermi - two*stm_broad
       Emax = efermi + stm_bias + two*stm_broad
    end if
    write(*,fmt='(4x,"Limits on integration: ",2f12.5," eV")') Emin+two*stm_broad, Emax-two*stm_broad
    write(*,fmt='(4x,"Thermal broadening at edges: ",f12.5," eV")') stm_broad
    ! Now loop over k-points and read eigenvalues
    do nk = 1,nkp
       read(17,*) idum,kx(nk),ky(nk),kz(nk),wtk(nk)
       do nev = 1,n_evals
          read(17,*) idum,eval
          ! Convert to eV (for easy comparison with STM)
          eval = eval*HaToeV
          tmp_evals(nev,nk) = eval
          if(eval>Emin.AND.eval<Emax) active_bands(nev) = 1
          if(flag_only_charge.AND.flag_range.AND.(eval>E_wf_min.AND.eval<E_wf_max)) active_bands(nev) = 1
       end do
    end do
    close(unit=17)
    ! Allocate space for active bands
    idum = 0
    do nev = 1, n_evals
       if(active_bands(nev)==1) idum = idum+1
    end do
    n_bands_active = idum
    allocate(eigenvalues(nkp,n_bands_active))
    allocate(band_no(n_bands_active))
    ! Store active bands
    eigenvalues = -1e30_double ! Clearly not active
    idum = 0
    do nev=1,n_evals
       if(active_bands(nev)==1) then
          idum = idum + 1
          band_no(idum) = nev
          if(.NOT.flag_only_charge) then
             if(flag_wf_range_Ef) then
                do nk=1,nkp
                   Eband = tmp_evals(nev,nk)-efermi
                   if(Eband>Emin.AND.Eband<Emax) &
                        eigenvalues(nk,idum) = tmp_evals(nev,nk)
                end do
             else
                do nk=1,nkp
                   if(tmp_evals(nev,nk)>Emin.AND.tmp_evals(nev,nk)<Emax) &
                        eigenvalues(nk,idum) = tmp_evals(nev,nk)
                end do
             end if
          end if
       end if
    end do
    deallocate(active_bands,tmp_evals)
  end subroutine read_eigenvalues_stm

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
    allocate(type_species(n_species),STAT=stat)
    allocate(atomicnum(n_species),STAT=stat)
    allocate(RadiusSupport(n_species),STAT=stat)
    allocate(InvSRange(n_species),STAT=stat)
    allocate(RadiusAtomf(n_species),STAT=stat)
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
