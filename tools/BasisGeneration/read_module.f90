module read

  use datatypes
  use GenComms, ONLY: cq_abort
  use global_module, ONLY: iprint
  
  implicit none

  character(len=80) :: block_file
  character(len=80), save :: pseudo_file_name, vps_file_name, vkb_file_name

  ! Local parameters 
  integer :: basis_size
  integer, parameter :: minimal = 1
  integer, parameter :: small = 2
  integer, parameter :: medium = 3
  integer, parameter :: full = 4

  integer :: energy_units ! Local for reading; 1 is Ha, 2 is eV
  real(double) :: energy_conv ! Define a factor for energy conversion (often 1.0)
  real(double) :: energy_semicore ! Threshold for semi-core states
  real(double) :: shallow_state_energy ! Energy to define shallow states (maybe too large)
  real(double), save :: gen_energy_semicore ! System-wide threshold for semi-core states
  real(double) :: width, prefac ! Defaults
  logical, save :: flag_gen_use_Vl
  logical :: flag_adjust_deltaE = .false.

  character(len=30), dimension(:), allocatable, save :: species_label
      
contains

  ! Read Conquest_ion_input into memory, and obtain general parameters
  subroutine read_general_input

    use numbers
    use species_module, only: n_species
    use input_module, only: load_input, fdf_integer, fdf_double, fdf_boolean, &
         fdf_block, block_start, block_end, input_array, fdf_endblock
    use pseudo_tm_info, ONLY: pseudo
    use pseudo_atom_info, ONLY: flag_default_cutoffs, flag_plot_output

    implicit none

    integer :: i, j, ios
    
    !
    ! Load the Conquest_ion_input file into memory
    !
    call load_input
    !
    ! Scan for general parameters
    !
    ! IO
    ! 
    iprint = fdf_integer('IO.Iprint',0)
    flag_plot_output = fdf_boolean('IO.PlotOutput',.false.)
    !
    ! Species
    !
    n_species = fdf_integer('General.NumberOfSpecies',1)
    allocate(pseudo(n_species)) ! Preserve compatibility with Conquest
    !
    ! Read species labels
    !
    allocate(species_label(n_species))
    if(fdf_block('SpeciesLabels')) then
       if(1+block_end-block_start<n_species) & 
            call cq_abort("Too few species in SpeciesLabels: ",&
            1+block_end-block_start,n_species)
       do i=1,n_species
          read (unit=input_array(block_start+i-1),fmt=*, iostat=ios) j, species_label(j)
          if(ios/=0) call line_read_error('SpeciesLabels',ios,i)
       end do
       call fdf_endblock
    else
       call cq_abort("You need to specify SpeciesLabels block")
    end if
    !
    ! Retain previous general flags? Read and use as defaults for atom
    !
    flag_gen_use_Vl = fdf_boolean('General.UseVl',.false.)
    gen_energy_semicore = fdf_double('General.SemicoreEnergy',-one)
    if(gen_energy_semicore>zero) &
         write(*,fmt='(4x,"Error: your semi-core threshold is positive ! ",f6.3)') gen_energy_semicore
    return
  end subroutine read_general_input

  ! Read input from Conquest_ion_input for a species
  ! This routine also reads input and output from Hamann code
  subroutine read_species_input(species)

    use numbers
    use input_module, ONLY: fdf_block, fdf_string, leqi, fdf_integer, fdf_double, fdf_boolean, &
         fdf_endblock
    use pseudo_atom_info, ONLY: paos, flag_use_Vl, val, deltaE_large_radius, deltaE_small_radius, &
         flag_default_cutoffs, pao_cutoff_energies, pao_cutoff_radii, pao_cutoff_default
    use units, ONLY: HaToeV
    use mesh, ONLY: mesh_type, hamann, siesta, alpha, beta, delta_r_reg

    ! Passed variables
    integer :: species

    ! Local variables
    character(len=80) :: input_string
    integer :: i
    logical :: flag_default
    
    if(fdf_block(species_label(species))) then
       !
       ! Set semi-core energy
       ! 
       energy_semicore = fdf_double('Atom.SemicoreEnergy',gen_energy_semicore)
       if(energy_semicore>zero) &
            write(*,fmt='(4x,"Error: your semi-core threshold is positive ! ",f6.3)') energy_semicore
       ! Energy which detects a state near zero which might have too large a radius
       shallow_state_energy = fdf_double('Atom.ShallowEnergy',-0.152_double)
       !
       ! Get Hamann input and output file names and read files
       !
       pseudo_file_name = fdf_string(80,'Atom.PseudopotentialFile',' ')
       call read_hamann_input(species)
       vkb_file_name    = fdf_string(80,'Atom.VKBFile',' ')
       call read_vkb(species)
       !
       ! Mesh
       !
       mesh_type = hamann
       alpha = fdf_double("Mesh.Alpha",alpha)
       beta = fdf_double("Mesh.Beta",beta)
       delta_r_reg = fdf_double("Atom.RegularSpacing",0.01_double)       
       !
       ! Form for zetas: compress or split norm
       !
       input_string = fdf_string(80,'Atom.ZetaForm','split')
       if(leqi(input_string(1:3),'spl')) then
          paos%flag_zetas = 1 ! Split-norm approach
       else if(leqi(input_string(1:3),'com')) then
          paos%flag_zetas = 2 ! Compression approach
       else
          call cq_abort("Unrecognised zeta form flag "//input_string(1:8))
       end if
       !
       ! Confinement
       !
       width = fdf_double("Atom.WidthConfine",one)
       prefac = fdf_double("Atom.PrefacConfine",zero)
       !
       ! Cutoffs
       !
       input_string = fdf_string(8,"Atom.Cutoffs","default")
       if(leqi(input_string(1:2),'en')) then ! Take reasonable keyword
          paos%flag_cutoff = pao_cutoff_energies
       else if(leqi(input_string(1:2),'ra')) then ! Again
          paos%flag_cutoff = pao_cutoff_radii
       else if(leqi(input_string(1:7),'default')) then
          paos%flag_cutoff = pao_cutoff_radii !default
       else
          call cq_abort("Unrecognised atomic cutoff flag "//input_string(1:8))
       end if
       !
       ! Energy units
       !
       energy_units = 1
       energy_conv = one
       input_string = fdf_string(2,"Atom.EnergyUnits","Ha")
       if(leqi(input_string(1:2),"eV")) then
          energy_units = 2
          energy_conv = one / HaToeV
          deltaE_large_radius = fdf_double("Atom.dE_large_radius",0.02_double)
          deltaE_small_radius = fdf_double("Atom.dE_small_radius",two)
          deltaE_large_radius = deltaE_large_radius / HaToeV
          deltaE_small_radius = deltaE_small_radius / HaToeV
       else if(leqi(input_string(1:2),"Ha")) then
          deltaE_large_radius = fdf_double("Atom.dE_large_radius",0.00073498_double)
          deltaE_small_radius = fdf_double("Atom.dE_small_radius",0.073498_double)
       end if
       if(flag_adjust_deltaE) then
          deltaE_large_radius = deltaE_large_radius*two
          write(*,fmt='(2x,"Shallow energy state present ", &
               "so adjusting large radius energy shift to ",f7.5," Ha")') deltaE_large_radius
          write(*,fmt='(2x,"Adjust Atom.ShallowEnergy if needed")')
          flag_adjust_deltaE = .false. ! Allow for other elements not to have this
       end if
       !
       ! Basis size
       !
       input_string = fdf_string(7,'Atom.BasisSize','none')
       if(leqi(input_string(1:7),'minimal')) then
          basis_size = minimal
       else if(leqi(input_string(1:5),'small')) then
          basis_size = small
       else if(leqi(input_string,'medium').OR.leqi(input_string(1:4),'none')) then ! Default
          basis_size = medium
       else if(leqi(input_string(1:4),'full').OR.leqi(input_string(1:5),'large')) then
          basis_size = full
       else
          call cq_abort("Error ! Unknown Atom.BasisSize specification: "//input_string(1:7))
       end if
       !
       ! Basis block
       !
       input_string = fdf_string(80,'Atom.BasisBlock','none')
       if(.NOT.(leqi(input_string(1:4),'none'))) then !.AND.paos%flag_cutoff==3)) then
          paos%n_shells = fdf_integer("Atom.PAO_N_Shells",0)
          if(paos%n_shells==0) call cq_abort("Number of PAO shells not specified in atom block")
          if(val%n_occ>paos%n_shells) &
               call cq_abort("More valence shells in PP than in Atom.BasisBlock ",val%n_occ,paos%n_shells)
          ! Close species block
          call fdf_endblock
          call read_basis_block(input_string)
          flag_default_cutoffs = .false.
       else
          ! Close species block
          call fdf_endblock
          flag_default_cutoffs = .true.
       end if
       !
       ! Polarisation
       !
       flag_default = .true.
       ! If there are no polarisation orbitals, change default
       if(basis_size == minimal .or. paos%n_shells==val%n_occ) flag_default = .false.
       paos%flag_perturb_polarise = fdf_boolean("Atom.Perturbative_Polarised",flag_default)
       if(paos%flag_perturb_polarise .and. flag_default) then
          ! This is for compatibility but should be false
          flag_use_Vl = fdf_boolean('Atom.UseVl',flag_gen_use_Vl)
          ! Do we have state to polarise? 
          paos%polarised_n = fdf_integer("Atom.PolarisedN",0)
          paos%polarised_l = fdf_integer("Atom.PolarisedL",-1)
          paos%polarised_shell = 0
          !
          ! Default to outer shell (or previous if outer shell has l=2)
          !
          if(paos%polarised_n==0.AND.paos%polarised_l==-1) then
             ! Outer-most occupied shell
             i = val%n_occ
             ! If it is l=2 (or more!) then default to one lower
             if(val%l(i)>1) i = i-1
             if(i==0) &
                  call cq_abort("We have one shell with l=2; I can't polarise this automatically.")
             paos%polarised_n = val%n(i)
             paos%polarised_l = val%l(i)
             paos%polarised_shell = i
             write(*,fmt='(4x,"Polarising shell ",i2," with n=",i2," and l=",i2)') &
                  i, paos%polarised_n, paos%polarised_l
          end if
       else if(paos%flag_perturb_polarise .and. (.not.flag_default)) then
          call cq_abort("Can't set Atom.Perturbative_Polarised T with no polarisation orbitals")
       end if
    else
       call cq_abort("Can't find species block for label "//species_label(species))
    end if ! if fdf_block(species(species)) - is this species defined ?!
    return
  end subroutine read_species_input

  ! Read the basis specification for an atom
  subroutine read_basis_block(basis_block_string)

    use numbers
    use pseudo_atom_info, ONLY: paos, pao_cutoff_energies, pao_cutoff_radii, pao_cutoff_default, &
         allocate_pao, allocate_pao_z
    use input_module

    implicit none

    character(len=80) :: basis_block_string

    integer :: j, k, ios, n_shells_read
    integer :: maxl, n_paos
    
    if(fdf_block(basis_block_string)) then
       !
       ! Initialise
       ! 
       call allocate_pao(paos%n_shells)
       paos%width = width
       paos%prefac = prefac
       maxl = 0
       n_paos = 0
       !
       ! Set number of shells to read
       !
       if(paos%flag_perturb_polarise) then
          n_shells_read = paos%n_shells-1
       else
          n_shells_read = paos%n_shells
       end if
       !
       ! Read zetas for each shell
       !
       do j=1,n_shells_read
          read (unit=input_array(block_start-1+j),fmt=*, iostat=ios) paos%n(j), paos%l(j), paos%nzeta(j)
          if(ios/=0) call line_read_error(basis_block_string,ios,j)
          if(paos%nzeta(j)>maxl) maxl = paos%nzeta(j)
          n_paos = n_paos + paos%nzeta(j)
       end do
       if(paos%flag_perturb_polarise) then
          call set_polarisation
          paos%nzeta(paos%n_shells) = paos%nzeta(paos%polarised_shell)
       end if
       !
       ! Allocate
       !
       call allocate_pao_z(maxl)
       !paos%total_paos = n_paos
       !
       ! Check for format error
       !
       if(block_start-1+2*n_shells_read>block_end) then
          write(*,fmt='(2x,"Not enough lines in atom block ",(a))') basis_block_string
          call cq_abort("Basis block error")
       end if
       !
       ! Read radii/energies
       !
       do j=1,n_shells_read
          ios = 0
          if(paos%flag_cutoff==pao_cutoff_energies) then ! Energies
             read (unit=input_array(block_start-1+n_shells_read+j),fmt=*, iostat=ios) &
                  (paos%energy(k,j),k=1,paos%nzeta(j))
             if(paos%energy(1,j)>two) write(*,fmt='("Warning: energy shift of ",f6.3, &
                  & "Ha is rather large.  Are these radii ?")') paos%energy(1,j)
          else
             read (unit=input_array(block_start-1+n_shells_read+j),fmt=*, iostat=ios) &
                  (paos%cutoff(k,j),k=1,paos%nzeta(j))
             if(paos%cutoff(1,j)<two) write(*,fmt='("Warning: radius of ",f6.3, &
                  & "a0 is rather small.  Are these energies ?")') paos%cutoff(1,j)
          end if
          if(ios/=0) call line_read_error(basis_block_string,ios,j)
       end do
       if(paos%flag_perturb_polarise) then
          paos%cutoff(:,paos%n_shells) = paos%cutoff(:,paos%polarised_shell)
          paos%energy(:,paos%n_shells) = paos%energy(:,paos%polarised_shell)
       end if
       !
       ! Check for confinement potentials
       !
       if(1+block_end-block_start>2*n_shells_read) then
          do j=1,n_shells_read
             read (unit=input_array(block_start-1+2*n_shells_read+j),fmt=*,iostat=ios) &
                  paos%width(j),paos%prefac(j)
          if(ios/=0) call line_read_error(basis_block_string,ios,j)
          end do
          if(paos%flag_perturb_polarise) then
             paos%width(paos%n_shells) = paos%width(paos%polarised_shell)
             paos%prefac(paos%n_shells) = paos%prefac(paos%polarised_shell)
          end if
       end if
       !
       ! Sort the energies/cutoffs so that largest is first
       !
       if(paos%flag_cutoff==pao_cutoff_radii) then
          ! We want to sort cutoffs large -> small
          do j=1,paos%n_shells
             call lsortr(paos%cutoff(:,j),paos%nzeta(j))
          end do
       else
          ! We want to sort energies small -> large
          do j=1,paos%n_shells
             call lsort(paos%energy(:,j),paos%nzeta(j))
          end do
       end if
       !
       ! Close block
       !
       call fdf_endblock
    else
       call cq_abort("Can't find species basis block "//basis_block_string)
    end if
    return
  end subroutine read_basis_block

  ! Find inner shell and set npao
  subroutine set_polarisation

    use numbers
    use pseudo_atom_info, ONLY: paos, val, flag_use_Vl
    
    implicit none

    ! Local variables
    integer :: i

    !
    ! Locate perturbed shell index, if necessary
    !
    if(paos%polarised_shell==0) then
       do i=val%n_occ,1,-1
          if(val%n(i)==paos%polarised_n.AND.val%l(i)==paos%polarised_l) then
             paos%polarised_shell=i
             write(*,fmt='(4x,"Polarising shell ",i3)') i
             exit
          end if
       end do
       if(paos%polarised_shell==0) call cq_abort("Can't find shell to polarise ",&
            paos%n(paos%polarised_shell), paos%l(paos%polarised_shell))
    end if
    ! Set n, l
    paos%l(paos%n_shells) = paos%l(paos%polarised_shell)+1
    paos%n(paos%n_shells) = max(paos%n(paos%polarised_shell),paos%l(paos%n_shells)+1)
    if(iprint>2) write(*,fmt='("For polarisation, we will perturb shell with n=",i2," and l=",i2)') &
         paos%n(paos%polarised_shell), paos%l(paos%polarised_shell)
    ! Check for inner shell
    paos%inner(paos%n_shells) = 0
    do i=val%n_occ,1,-1
       if(paos%l(i)==paos%l(paos%n_shells)) paos%inner(paos%n_shells) = i
    end do
    !
    ! Check whether there is an inner shell
    !
    !if(val%inner(paos%polarised_shell)>0) &
    !     call cq_abort("Can't use perturbative polarisation on shell with nodes: ", &
    !     paos%polarised_n, paos%polarised_l)
    !
    ! Set npao for polarisation - use same as the shell being perturbed
    !
    paos%npao(paos%n_shells) = val%npao(paos%polarised_shell)
    !
    ! Checks for perturbing s orbital
    !
    ! Li, Be (npao set to one because l=0)
    if(paos%l(paos%n_shells)==1 .AND. paos%inner(paos%n_shells)==0) paos%npao(paos%n_shells) = 1
    ! We need sign change if s orbital is outermost and there is a semi-core p orbital
    if(paos%l(paos%n_shells)==1 .AND. paos%polarised_shell==val%n_occ &
         .AND. paos%inner(paos%n_shells)>0) paos%pol_pf = -one
    ! Use Vl if perturbing a p orbital with a valence or semi-core d orbital
    if(paos%l(paos%n_shells)==2) then
       do i=1,val%n_occ
          if(paos%l(i)==2) flag_use_Vl = .true.
       end do
    end if
    return
  end subroutine set_polarisation
  
  ! Error reading line
  subroutine line_read_error(basis_block_string,ios,j)

    use GenComms, only: cq_abort

    implicit none
    
    character(len=*) :: basis_block_string
    integer :: ios, j
    
    write(*,fmt='(2x,"Problem with reading atom block ",(a))') basis_block_string
    write(*,fmt='(2x,"Problem with shell ",i2)') j
    if(ios<0) then
       write(*,fmt='(2x,"End of record or end of file; please check block format")')
       call cq_abort("Basis block error")
    else if(ios>0) then
       call cq_abort("System-dependent read error: ",ios)
    endif
  end subroutine line_read_error
  
  ! Read semi-local potentials output by a DRB patch to Hamann's code
  ! Now read KB potentials
  subroutine read_vkb(i_species)

    use numbers
    use pseudo_tm_info, ONLY: pseudo
    use input_module, ONLY: io_assign, io_close
    use mesh, ONLY: siesta
    use pseudo_atom_info, ONLY: val, allocate_val, local_and_vkb, allocate_vkb, hamann_version, deltaE_large_radius
    
    implicit none

    ! Passed variables
    integer :: i_species

    ! Local variables
    integer :: ios, lun, ngrid, ell, en, i, j, n_occ, n_read
    integer :: n_shells, n_nl_proj, this_l, number_of_this_l
    integer, dimension(0:4) :: count_func
    integer, dimension(3,0:4) :: index_count_func
    character(len=2) :: char_in
    character(len=80) :: line
    logical :: flag_core_done = .false.
    real(double) :: dummy, dummy2, highest_energy

    !
    ! Zero arrays
    !
    count_func = 0
    index_count_func = 0
    !
    ! Open file
    !
    call io_assign(lun)
    open(unit=lun, file=vkb_file_name, status='old', iostat=ios)
    if ( ios > 0 ) call cq_abort('Error opening pseudopotential file: '//vkb_file_name)
    !
    ! Find number of shells and allocate space
    !
    read(lun,'(a)') line
    if(len(trim(line))<9) then ! No version string
       read(line,*) char_in,n_shells
       hamann_version = 0
    else
       read(line,*) char_in,n_shells, hamann_version
       write(*,*) 'Hamann code version ',hamann_version
    end if
    if(iprint>4) write(*,fmt='("Reading Hamann output, with ",i3," valence shells")') n_shells
    call allocate_val(n_shells)
    !
    ! Read in table of shells, energies and occupancies and test for ordering
    !
    do i = 1,n_shells
       read(lun,*) char_in,val%n(i),val%l(i),val%occ(i),val%en_ps(i)
       !
       ! Ensure that shells are ordered; relatively simplistic but effective
       !
       if(i>1) then
          do j=i,2,-1
             if(val%en_ps(j)<val%en_ps(j-1).AND.val%occ(j)>RD_ERR) then ! Only re-order occupied
                en = val%n(j-1)
                ell = val%l(j-1)
                dummy = val%occ(j-1)
                dummy2 = val%en_ps(j-1)
                val%n(j-1) = val%n(j)
                val%l(j-1) = val%l(j)
                val%occ(j-1) = val%occ(j)
                val%en_ps(j-1) = val%en_ps(j)
                val%n(j) = en
                val%l(j) = ell
                val%occ(j) = dummy
                val%en_ps(j) = dummy2
             end if
          end do
       end if
    end do
    ! 
    n_occ = 0
    if(iprint>3) then
       write(*,fmt='("Valence shells from pseudopotential calculation")')
       write(*,fmt='("  n  l    occ    energy")')
    end if
    !
    ! Check for inner shells and assign pseudo-n value (for nodes)
    !
    highest_energy = -ten ! For comparison
    do i = 1,n_shells
       if(iprint>3) write(*,fmt='(2i3,f7.2,f10.4)') val%n(i),val%l(i),val%occ(i),val%en_ps(i)
       ! Check for inner shells: count shells with this l and store
       this_l = val%l(i)
       number_of_this_l = count_func(this_l) + 1
       count_func(this_l) = number_of_this_l
       index_count_func(number_of_this_l, this_l) = i
       if(number_of_this_l>1) val%inner(i) = index_count_func(number_of_this_l-1,this_l)
       ! Set n for PAO (sets number of nodes)
       val%npao(i) = this_l + number_of_this_l 
       ! Check for semi-core state
       if(val%en_ps(i)<energy_semicore) val%semicore(i) = 1
       ! Count occupied states
       if(val%occ(i)>RD_ERR.OR.val%en_ps(i)<-1e-6_double) then
          n_occ = n_occ + 1
          highest_energy = max(val%en_ps(i),highest_energy)
       end if
    end do
    if(iprint>3) write(*,fmt='(i2," valence shells, with ",i2," occupied")') n_shells, n_occ
    if(highest_energy>shallow_state_energy) then
       flag_adjust_deltaE = .true.
    end if
    val%n_occ = n_occ
    !
    ! Read grid, charge, partial core, local potential, semilocal potentials
    !
    read(lun,*) char_in, ngrid
    call allocate_vkb(ngrid,i_species)
    !
    ! Read and store first block: local data
    !
    if(pseudo(i_species)%flag_pcc) then
       ! Read logarithmic mesh, atomic charge, partial core charge, local potential
       do i=1,ngrid
          read(lun,*) local_and_vkb%rr(i),local_and_vkb%charge(i), &
               local_and_vkb%pcc(i), local_and_vkb%local(i), &
               (local_and_vkb%semilocal_potential(i,j),j=0,pseudo(i_species)%lmax)
       end do
    else
       ! Read logarithmic mesh, atomic charge, dummy for core, local potential 
       do i=1,ngrid
          read(lun,*) local_and_vkb%rr(i),local_and_vkb%charge(i), &
               dummy, local_and_vkb%local(i), &
               (local_and_vkb%semilocal_potential(i,j),j=0,pseudo(i_species)%lmax)
       end do
    end if
    !
    ! Read number of VKB projectors
    !
    read(lun,*) line!char_in,i,j,ios ! We could this as a check on the number of projectors
    if(iprint>4) write(*,fmt='("VKB projector energies"/,"  l  energies")')
    n_read = 0
    do ell=0,pseudo(i_species)%lmax
       read(lun,*) char_in,j,pseudo(i_species)%pjnl_ekb(n_read+1:n_read+local_and_vkb%n_proj(ell))
       if(iprint>4) write(*,fmt='(i3,4f10.5)') ell, &
            pseudo(i_species)%pjnl_ekb(n_read+1:n_read+local_and_vkb%n_proj(ell))
       n_read = n_read+local_and_vkb%n_proj(ell)
    end do
    read(lun,*) char_in, ngrid
    local_and_vkb%r_vkb = local_and_vkb%rr(ngrid)
    local_and_vkb%ngrid_vkb = ngrid
    !
    ! Read VKB projectors
    !
    do i=1,ngrid
       read(lun,*) ((local_and_vkb%projector(i,j,ell), &
            j=1,local_and_vkb%n_proj(ell)), ell=0,pseudo(i_species)%lmax)
    end do
    call io_close(lun)
    return
  end subroutine read_vkb

  subroutine read_hamann_input(i_species)

    use numbers
    use input_module, ONLY: io_assign, io_close
    use pseudo_tm_info, ONLY: alloc_pseudo_info, pseudo
    use spline_module, ONLY: spline
    use radial_xc, ONLY: flag_functional_type, init_xc, functional_lda_pz81, functional_gga_pbe96, &
         functional_description
    use pseudo_atom_info, ONLY: local_and_vkb, val, input_file_length
    
    implicit none

    ! Passed variables
    integer :: i_species

    ! Local variables
    integer :: lun, nc, nv, iexc, i_shell, en, ell, icmod, i, j, zeta, ios, n_nl_proj
    integer :: input_lines
    real(double) :: fill, zval, z, zcore
    character(len=80) :: a
    character(len=2) :: sym
    character(len=4) :: file_format
    character(len=132) :: line
    logical :: done

    !
    ! Open appropriate file
    !
    call io_assign(lun)
    open(unit=lun, file=pseudo_file_name, status='old', iostat=ios)
    if ( ios > 0 ) call cq_abort('Error opening Hamann input file: '//pseudo_file_name)
    input_file_length = 0
    pseudo(i_species)%filename = pseudo_file_name
    a = get_hamann_line(lun)
    read(a,*) sym,z,nc,nv,iexc,file_format
    write(*,fmt='(/"Information about pseudopotential for species: ",a2/)') sym
    pseudo(i_species)%z = z
    !
    ! Assign and initialise XC functional for species
    !
    if(iexc==3) then
       flag_functional_type = functional_lda_pz81
    else if(iexc==4) then
       flag_functional_type = functional_gga_pbe96
    else if(iexc<0) then
       flag_functional_type = iexc
    else
       call cq_abort("Error: unrecognised iexc value: ",iexc)
    end if
    call init_xc
    write(*,fmt='("There are ",i2," core and ",i2," valence shells")') nc,nv
    a = get_hamann_line(lun)
    !
    ! Read n, l, filling for core
    !
    zcore = zero
    do i_shell = 1, nc
       read(a,*) en,ell,fill
       zcore = zcore + fill
       a = get_hamann_line(lun)
    end do
    pseudo(i_species)%zcore = zcore
    ! Read n, l, filling for valence
    zval = zero
    do i_shell = 1, nv
       read(a,*) en,ell,fill
       zval = zval + fill
       a = get_hamann_line(lun)
    end do
    pseudo(i_species)%zval = zval
    write(*,fmt='("The atomic number is",f6.2,", with valence charge ",f5.2," (core electrons: ",f5.2,")")') z,zval,zcore
    !
    ! lmax
    !
    read(a,*) pseudo(i_species)%lmax
    write(*,fmt='("Maximum angular momentum for pseudopotential is l=",i1)') pseudo(i_species)%lmax
    !
    ! Projector radii etc
    !
    local_and_vkb%core_radius = zero
    a = get_hamann_line(lun)
    do ell = 0, pseudo(i_species)%lmax
       read(a,*) en, local_and_vkb%core_radius(ell)
       write(*,'("l=",i1," core radius ",f6.3," bohr")') ell, local_and_vkb%core_radius(ell)
       a = get_hamann_line(lun)
    end do
    !
    ! Local potential read above in final step of loop, so leave
    !
    !a = get_hamann_line(lun)
    !
    ! Numbers of projectors
    !
    local_and_vkb%n_proj = 0
    local_and_vkb%n_nl_proj = 0
    a = get_hamann_line(lun)
    do ell = 0, pseudo(i_species)%lmax
       read(a,*) en,local_and_vkb%n_proj(en),fill
       local_and_vkb%n_nl_proj = local_and_vkb%n_nl_proj + local_and_vkb%n_proj(en)
       a = get_hamann_line(lun)
    end do
    write(*,fmt='("Total number of VKB projectors: ",i2)') local_and_vkb%n_nl_proj
    call alloc_pseudo_info(pseudo(i_species),local_and_vkb%n_nl_proj)
    i=0
    do ell=0,pseudo(i_species)%lmax
       zeta=0
       do j=1,local_and_vkb%n_proj(ell)
          i=i+1
          zeta = zeta+1
          pseudo(i_species)%pjnl_l(i) = ell
          pseudo(i_species)%pjnl_n(i) = zeta
       end do
    end do
    !
    ! Partial core correction ?
    !
    pseudo(i_species)%flag_pcc = .false.
    read(a,*) icmod
    if(icmod>0) then
       pseudo(i_species)%flag_pcc = .true.
       write(*,fmt='("This pseudopotential includes partial core corrections"/)') 
    else
       write(*,fmt='("This pseudopotential does not include partial core corrections"/)') 
    end if
    ! Last two categories: tests and grid size
    a = get_hamann_line(lun)
    a = get_hamann_line(lun)
    call io_close(lun)
  end subroutine read_hamann_input
  
  ! Check for comment markers
  function get_hamann_line(lun)

    use pseudo_atom_info, ONLY: input_file_length, input_file
    implicit none
   
    integer :: lun
    character(len=80) :: get_hamann_line

    character(len=80) :: a

    read(lun,'(a)') a
    input_file_length = input_file_length+1
    input_file(input_file_length) = a
    get_hamann_line = adjustl(a)
    do while(get_hamann_line(1:1).eq.'#')
       read(lun,'(a)') a
       input_file_length = input_file_length+1
       input_file(input_file_length) = a
       get_hamann_line = adjustl(a)
    end do
    return
  end function get_hamann_line

  ! Set up PAO basis in case of default, or check user-specified basis
  subroutine set_pao_initial(i_species)

    use numbers
    use species_module, ONLY: n_species
    use pseudo_atom_info, ONLY: paos, allocate_pao, allocate_pao_z, flag_default_cutoffs, val
    use pseudo_tm_info, ONLY: pseudo
    use periodic_table, ONLY: pte

    implicit none

    ! Passed variables
    integer :: i_species

    ! Local variabls
    integer :: j, i_highest_occ, this_l, number_of_this_l, i_end
    integer :: i_shell, ell, en, n_paos, n_shells, ell_hocc, max_zeta, n_pol_zeta
    integer, dimension(0:4) :: count_func
    logical :: flag_confine

    flag_confine = .false.
    write(*,fmt='(/"Species ",i2," is ",a2)') i_species, pte(nint(pseudo(i_species)%z))
    !
    ! Defaults
    !
    if(flag_default_cutoffs) then 
       write(*,fmt='(/"Default basis chosen"/)')
       !
       ! Set numbers of functions
       !
       n_shells = val%n_occ+1 ! Default
       if(basis_size==minimal) then ! One zeta for valence
          max_zeta = 1
          n_shells = val%n_occ
       else if(basis_size==small) then ! As above plus polarisation
          max_zeta = 1
          n_pol_zeta = 1
       else if(basis_size==medium) then ! Two zeta for valence plus polarisation
          max_zeta = 2
          n_pol_zeta = 1
       else if(basis_size==full) then ! Three zeta for all channels - default
          max_zeta = 3
          n_pol_zeta = 3
       end if
       !
       ! Allocate space
       !
       count_func = 0
       write(*,fmt='("  n  l  nodes  zetas")')
       call allocate_pao(n_shells)
       call allocate_pao_z(max_zeta)
       paos%width = width
       paos%prefac = prefac
       !
       ! Set information about PAOs: valence
       !
       do i_shell = 1,val%n_occ          
          this_l = val%l(i_shell)
          number_of_this_l = count_func(this_l) + 1
          count_func(this_l) = number_of_this_l
          paos%n(i_shell) = val%n(i_shell)
          paos%l(i_shell) = this_l
          paos%npao(i_shell) = val%npao(i_shell)
          if(val%semicore(i_shell)>0) then
             paos%nzeta(i_shell) = 1 
             write(*,fmt='(2i3,2i7, "  semi-core state")') paos%n(i_shell), paos%l(i_shell),&
                  paos%npao(i_shell) - paos%l(i_shell) - 1, &
                  paos%nzeta(i_shell)
          else
             paos%nzeta(i_shell) = max_zeta 
             write(*,fmt='(2i3,2i7)') paos%n(i_shell), paos%l(i_shell),&
                  paos%npao(i_shell) - paos%l(i_shell) - 1, &
                  paos%nzeta(i_shell)
          end if
          if(val%inner(i_shell)>0) then
             paos%inner(i_shell)=val%inner(i_shell)
          end if
       end do
       !
       ! Set information about PAOs: polarisation
       !
       if(basis_size>minimal) then
          i_shell = paos%n_shells
          paos%nzeta(i_shell) = n_pol_zeta
          !
          ! Determine polarisation shell: assume l(highest occupied)+1 unless that is l=2
          !
          i_highest_occ = val%n_occ
          if(paos%l(i_highest_occ)>1) then
             i_highest_occ = i_highest_occ - 1
             if(i_highest_occ==0) &
                  call cq_abort("Only one shell, with l=2; I can't polarise this automatically.")
          end if
          !
          ! Check for equivalent l in pseudopotentials - switch to perturbation if necessary
          !
          if(paos%l(i_highest_occ)+1>pseudo(i_species)%lmax.and.(.not.paos%flag_perturb_polarise)) then
             write(*,fmt='("Warning ! l for polarisation is greater than lmax for pseudopotential.&
                  &  We have to use perturbation.")')
             paos%flag_perturb_polarise = .true.
             paos%polarised_l = paos%l(i_highest_occ)
             paos%polarised_n = paos%n(i_highest_occ)
          end if
          !
          ! Set l, n, npao (for nodes)
          !
          if(paos%flag_perturb_polarise) then
             call set_polarisation
             write(*,fmt='(2i3,2i7, "  perturbed polarisation state")') paos%n(i_shell), &
                  paos%l(i_shell), paos%npao(i_shell) - paos%l(i_shell), paos%nzeta(i_shell)
          else
             paos%l(i_shell) = paos%l(i_highest_occ)+1
             paos%n(i_shell) = max(paos%n(i_highest_occ),paos%l(i_highest_occ)+1)
             paos%npao(i_shell) = paos%l(i_shell)+1
             ! Check for inner shell with same l value
             if(count_func(paos%l(i_shell))>0) then
                paos%npao(i_shell) = paos%npao(i_shell)+1
                do j = i_shell-1,1,-1
                   if(paos%l(j)==paos%l(i_shell)) paos%inner(i_shell) = j
                end do
             end if
             write(*,fmt='(2i3,2i7, "  polarisation state")') paos%n(i_shell), paos%l(i_shell),&
                  paos%npao(i_shell) - paos%l(i_shell) - 1, &
                  paos%nzeta(i_shell)
          end if
       end if
       !
       ! User-specified basis
       !
    else 
       write(*,fmt='("User-specified basis")')
       count_func = 0
       write(*,fmt='("  n  l  nodes  zetas")')
       !
       ! Create npao - used for number of nodes - and check consistency for n/l
       !
       paos%inner = 0
       do i_shell = 1,val%n_occ
          ! Checks
          if( paos%n(i_shell) /= val%n(i_shell) ) &
               call cq_abort("Mismatch in valence n value: ", paos%n(i_shell),val%n(i_shell))
          if( paos%l(i_shell) /= val%l(i_shell) ) &
               call cq_abort("Mismatch in valence l value: ", paos%l(i_shell),val%l(i_shell))
          if(paos%l(i_shell)>pseudo(i_species)%lmax) then
             write(*,fmt='(4x,"Max l value from pseudopotential: ",i2)') pseudo(i_species)%lmax
             write(*,fmt='(4x,"Max l value specified for PAO   : ",i2)') paos%l(i_shell)
             call cq_abort("Angular momentum mismatch between PS and PAO; please correct.")
          end if
          ! Confine?
          if(paos%prefac(i_shell)>RD_ERR) flag_confine = .true.
          this_l = val%l(i_shell)
          count_func(this_l) = count_func(this_l) + 1
          ! Set npao and write out
          paos%npao(i_shell) = val%npao(i_shell)
          write(*,fmt='(2i3,2i7)') paos%n(i_shell), paos%l(i_shell),&
                  paos%npao(i_shell) - paos%l(i_shell) - 1, &
                  paos%nzeta(i_shell)
          ! Check for inner shell and reasonable l
          if(val%inner(i_shell)>0) then
             paos%inner(i_shell)=val%inner(i_shell)
          end if
       end do ! i_shell = 1,val%n_occ
       !
       ! And again for unoccupied (polarisation) shells
       !
       if(paos%n_shells>val%n_occ) then
          if(paos%flag_perturb_polarise) then
             i_end = paos%n_shells - 1
          else
             i_end = paos%n_shells
          end if
          do i_shell = val%n_occ+1,i_end
             if(paos%l(i_shell)>pseudo(i_species)%lmax) then
                write(*,fmt='(4x,"Max l value from pseudopotential: ",i2)') pseudo(i_species)%lmax
                write(*,fmt='(4x,"Max l value specified for PAO   : ",i2)') paos%l(i_shell)
                call cq_abort("Angular momentum mismatch between PS and PAO; please correct.")
             end if
             ! Set npao and write out
             paos%npao(i_shell) = paos%l(i_shell)+1
             ! Check for inner shell with same l value
             if(count_func(paos%l(i_shell))>0) then
                paos%npao(i_shell) = paos%npao(i_shell)+1
                do j = i_shell-1,1,-1
                   if(paos%l(j)==paos%l(i_shell)) paos%inner(i_shell) = j
                end do
             end if
             write(*,fmt='(2i3,2i7)') paos%n(i_shell), paos%l(i_shell),&
                  paos%npao(i_shell) - paos%l(i_shell) - 1, &
                  paos%nzeta(i_shell)
          end do
          ! Add perturbative polarisation shell if selected
          if(paos%flag_perturb_polarise) then
             i_shell = paos%n_shells
             call set_polarisation
             paos%nzeta(i_shell) = paos%nzeta(paos%polarised_shell)
             write(*,fmt='(2i3,2i7, "  perturbed polarisation state")') paos%n(i_shell), &
                  paos%l(i_shell), paos%npao(i_shell) - paos%l(i_shell), paos%nzeta(i_shell)
          end if
       end if
       if(flag_confine) then
          write(*,fmt='("Using exponential confinement for PAOs")')
          write(*,fmt='("   Prefactor       Width")')
          do i_shell = 1,paos%n_shells
             write(*,fmt='(2f12.3)') paos%prefac(i_shell), paos%width(i_shell)
          end do
       end if
    end if ! Default cutoffs
  end subroutine set_pao_initial

  ! These routines were written to read the output from the Hamann code before realising that
  ! (a) we needed the semi-local potentials and (b) that the logarithmic radial mesh gave
  ! a much more accurate result.  They're kept below for reference but aren't used.
  subroutine read_abinit

    use datatypes
    use numbers
    use species_module, ONLY: n_species
    use input_module, ONLY: io_assign, io_close
    use pseudo_tm_info, ONLY: alloc_pseudo_info, pseudo, rad_alloc
    use spline_module, ONLY: spline
    
    implicit none

    character(len=80) :: a
    integer :: i_species, ios, lun, loc_lmax
    integer :: pseudo_format, pseudo_xc, lmax, l_local, grid_size
    integer :: i,j,l, n_nl_proj, i_proj, i_last, zeta
    real(double) :: zatom, zion, pcc_radius, pcc_flag, r, rlast, delta, yp1, ypn
    real(double), dimension(:,:), allocatable :: temp_array
    integer, dimension(0:6) :: n_proj ! Temporary - this is now in local_and_vkb
    
    allocate(pseudo(n_species))
    do i_species = 1,n_species
       ! open file
       call io_assign(lun)
       open(unit=lun, file=pseudo_file_name, status='old', iostat=ios)
       if ( ios > 0 ) call cq_abort('Error opening pseudopotential file: '//pseudo_file_name)
       pseudo(i_species)%filename = pseudo_file_name
       ! Read title line (and print ?)
       read (lun, '(a)') a
       write(*,fmt='(2x,"Title line for file: ",(a),/,4x,a80)') pseudo_file_name,a
       ! Read z atom, zion
       read(lun,*) zatom, zion
       write(*,fmt='(4x,"Atomic number: ",f6.1,/,4x,"Ionic charge: ",f6.1)') zatom, zion
       pseudo(i_species)%z = zatom
       pseudo(i_species)%zval = zion
       ! Line 3: format, XC, lmax, local, grid, extra
       read(lun,*) pseudo_format, pseudo_xc, lmax, l_local, grid_size
       pseudo(i_species)%lmax = lmax
       n_proj = 0
       ! PCC line
       read(lun,*) pcc_radius, pcc_flag
       if(pcc_flag>zero) pseudo(i_species)%flag_pcc = .true.
       ! Numbers of projectors for l=0 to llmax; Now we can work out n_pjnl (number of NL projectors)
       read(lun,*) (n_proj(l),l=0,lmax)
       n_nl_proj = 0
       do l=0,lmax
          n_nl_proj = n_nl_proj + n_proj(l)
          if(n_proj(l)==0) loc_lmax = lmax-1
       end do
       if(loc_lmax<lmax) then
          write(*,*) '# There is one value of l with zero KB projectors; changing lmax to ',loc_lmax
          lmax = loc_lmax
       end if
       pseudo(i_species)%lmax = lmax
       call alloc_pseudo_info(pseudo(i_species),n_nl_proj)
       i=0
       write(*,*) 'Grid size: ',grid_size
       do l=0,lmax
          zeta=0
          do j=1,n_proj(l)
             i=i+1
             zeta = zeta+1
             !call rad_alloc(pseudo(i_species)%pjnl(i),grid_size)
             pseudo(i_species)%pjnl_l(i) = l
             pseudo(i_species)%pjnl_n(i) = zeta
          end do
       end do
       write(*,*) 'Number of NL projectors: ',n_nl_proj
       !call alloc_pseudo_info(pseudo(i_species),n_nl_proj)
       ! Extension (ignore ?)
       read(lun,*) i
       ! Loop over l
       n_nl_proj = 0
       do l=0,lmax
          ! read l, energies
          write(*,*) 'l, nproj: ',l,n_proj(l)
          if(n_proj(l)>0) then
             allocate(temp_array(n_proj(l),grid_size))
             read(lun,*) j, (pseudo(i_species)%pjnl_ekb(n_nl_proj+i),i=1,n_proj(l))
             write(*,*) 'KB energies for l=',l,pseudo(i_species)%pjnl_ekb(n_nl_proj+1:n_nl_proj+n_proj(l))
             ! Read index, grid, projectors
             rlast = zero
             delta = zero
             i_last = 0
             do i=1,grid_size
                !read(lun,*) j,r,(pseudo(i_species)%pjnl(n_nl_proj+i_proj)%f(i),i_proj=1,nproj(l))
                read(lun,*) j,r,(temp_array(i_proj,i),i_proj=1,n_proj(l))
                do i_proj = 1,n_proj(l)
                   if(abs(temp_array(i_proj,i))>RD_ERR) then
                      i_last = i
                   end if
                end do
                if(i==2) then
                   delta = r-rlast
                   write(*,*) 'Delta: ',delta
                else if(i>2) then
                   if((abs(r-rlast)-delta)>RD_ERR) write(*,*) 'Possible delta error: ',r-rlast
                end if
                rlast = r
             end do
             write(*,*) 'Last non-zero entry for l=',l,' is ',i_last
             do i_proj=1,n_proj(l)
                call rad_alloc(pseudo(i_species)%pjnl(n_nl_proj+i_proj),i_last)
             end do
             do i=1,i_last
                do i_proj=1,n_proj(l)
                   pseudo(i_species)%pjnl(n_nl_proj+i_proj)%f(i) = temp_array(i_proj,i)
                end do
             end do
             deallocate(temp_array)
             pseudo(i_species)%pjnl(n_nl_proj+1:n_nl_proj+n_proj(l))%delta = delta
             ! NB we use i_last-1 because the first entry is r=0
             pseudo(i_species)%pjnl(n_nl_proj+1:n_nl_proj+n_proj(l))%cutoff = real(i_last-1,double)*delta
             write(*,*) 'Cutoff: ',pseudo(i_species)%pjnl(n_nl_proj+1)%cutoff
             yp1= zero
             ypn= zero
             do i=1,n_proj(l)
                call spline(pseudo(i_species)%pjnl(n_nl_proj+i)%n, pseudo(i_species)%pjnl(n_nl_proj+i)%delta, &
                     pseudo(i_species)%pjnl(n_nl_proj+i)%f, yp1, ypn, pseudo(i_species)%pjnl(n_nl_proj+i)%d2)
             end do
             n_nl_proj = n_nl_proj+n_proj(l)
          end if
       end do ! l=0,lmax
       ! Read lloc
       read(lun,*) i
       ! Read index, grid, potential
       call rad_alloc(pseudo(i_species)%vlocal,grid_size)
       pseudo(i_species)%vlocal%delta = delta
       pseudo(i_species)%vlocal%cutoff = grid_size*delta
       do i=1,grid_size
          read(lun,*) j,r,pseudo(i_species)%vlocal%f(i)
       end do
       ! PCC
       if(pseudo(i_species)%flag_pcc) then
          ! Find an appropriate cutoff
          allocate(temp_array(grid_size,1))
          i_last = 0
          do i=1,grid_size
             read(lun,*) j,r,temp_array(i,1)
             if(abs(temp_array(i,1))>RD_ERR) then
                i_last = i
             end if
          end do
          call rad_alloc(pseudo(i_species)%chpcc,i_last)
          pseudo(i_species)%chpcc%delta = delta
          pseudo(i_species)%chpcc%cutoff = i_last*delta
          do i=1,i_last
             pseudo(i_species)%chpcc%f(i) = temp_array(i,1)!/fourpi
          end do
          deallocate(temp_array)
          write (*,fmt='(10x,a)') "P.C.C. is taken into account."
       end if
       !call rad_alloc( pseudo(i_species)%vna, grid_size )
       !pseudo(i_species)%vna%delta = delta
       !pseudo(i_species)%vna%cutoff = grid_size*delta
       call io_close(lun)
    end do
    
  end subroutine read_abinit

  subroutine read_upf

    use datatypes
    use numbers
    use species_module, ONLY: n_species
    use input_module, ONLY: io_assign, io_close
    use pseudo_tm_info, ONLY: alloc_pseudo_info, pseudo, rad_alloc
    use spline_module, ONLY: spline
    use periodic_table, ONLY: pte
    
    implicit none

    character(len=80) :: a, str
    integer :: i_species, ios, lun, lines, st, rem, i_last
    integer :: pseudo_format, pseudo_xc, lmax, l_local, grid_size, verlen
    integer :: i,j,l, n_nl_proj, i_proj, ele, nproj_per_l
    real(double) :: zatom, zion, pcc_radius, pcc_flag, r, rlast, delta, yp1, ypn, temp(100)
    integer, dimension(0:6) :: n_proj

    allocate(pseudo(n_species))
    do i_species = 1,n_species
       ! open file
       call io_assign(lun)
       open(unit=lun, file=pseudo_file_name, status='old', iostat=ios)
       if ( ios > 0 ) call cq_abort('Error opening pseudopotential file: '//pseudo_file_name)
       pseudo(i_species)%filename = pseudo_file_name
       ! Read opening line
       str = get_upf_line(lun)
       if(str(1:4)/='<UPF') then
          call cq_abort('File error: not UPF format: '//str)
       else
          verlen = len(trim(str))-16
          write(*,*) 'UPF pseudopotential, version: ',str(15:14+verlen)
       end if
       str = get_upf_line(lun)
       if (str(1:9) .eq. '<PP_INFO>') then
          str = get_upf_line(lun)
          write(*,*) 'PP_INFO: ',str
          do while(str .ne. '</PP_INFO>')
             str = get_upf_line(lun)
             write(*,*) 'PP_INFO: ',str
          end do
       else
          write(*,*) 'No PP_INFO: ',trim(str)
       end if
       ! PP_HEADER
       str = get_upf_line(lun)
       if(str(1:10).eq.'<PP_HEADER') then
          write(*,*) 'Header'
          ! Read in a bunch of useful data
          ! Element
          str = get_upf_line(lun)
          do while(str(1:7).ne.'element')
             str = get_upf_line(lun)
          end do
          write(*,*) 'Element: ',str(10:11)
          ! We need to work out Zion from the element !
          do ele = 1, 103
             if(str(10:11).eq.pte(ele)) then
                pseudo(i_species)%z = ele
                write(*,*) 'Element number: ',ele
             end if
          end do
          ! Pseudo type
          str = get_upf_line(lun)
          do while(str(1:11).ne.'pseudo_type')
             str = get_upf_line(lun)
          end do
          if(str(14:15).ne.'NC') then
             call cq_abort("Wrong PS type in UPF: "//str(14:15))
          end if
          ! We may need to come back for relativistic
          ! core correction
          str = get_upf_line(lun)
          do while(str(1:15).ne.'core_correction')
             str = get_upf_line(lun)
          end do
          if((str(18:18).eq.'T').or.(str(18:18).eq.'t')) then
             pseudo(i_species)%flag_pcc = .true.
             write (*,fmt='(10x,a)') "P.C.C. is taken into account."
          end if
          ! functional
          str = get_upf_line(lun)
          do while(str(1:10).ne.'functional')
             str = get_upf_line(lun)
          end do
          write(*,*) 'Functional: ',str(13:15)
          ! We really need (a) a way to record and output the functional and (b) account for functional in CQ
          ! Valence
          str = get_upf_line(lun)
          do while(str(1:9).ne.'z_valence')
             str = get_upf_line(lun)
          end do
          read(str(12:19),*) pseudo(i_species)%zval
          write(*,*) 'Valence charge: ',pseudo(i_species)%zval
          ! rho_cutoff is a suggested energy for charge density expansion (?)
          ! l max
          str = get_upf_line(lun)
          do while(str(1:5).ne.'l_max')
             str = get_upf_line(lun)
          end do
          read(str(8:8),*) pseudo(i_species)%lmax
          write(*,*) 'Max l value: ',pseudo(i_species)%lmax
          ! mesh
          str = get_upf_line(lun)
          do while(str(1:9).ne.'mesh_size')
             str = get_upf_line(lun)
          end do
          read(str(12:17),*) grid_size
          write(*,*) 'Grid size: ',grid_size
          ! Projectors
          str = get_upf_line(lun)
          do while(str(1:14).ne.'number_of_proj')
             str = get_upf_line(lun)
          end do
          read(str(17:17),*) n_nl_proj
          write(*,*) 'Number of projectors: ',n_nl_proj
          ! We assume that there are the same number of projectors for each l
          ! If this is NOT the case we're in trouble...  We would then need to
          ! read the entire file, recording the projector and l for each beta
          ! For now I assume that we don't need to do this DRB 2016/06/10
          nproj_per_l = n_nl_proj/(pseudo(i_species)%lmax+1)
          if(nproj_per_l*(pseudo(i_species)%lmax+1)/=n_nl_proj) write(*,*) 'Problem with projectors ?'
          n_proj(0:pseudo(i_species)%lmax) = nproj_per_l
          call alloc_pseudo_info(pseudo(i_species),n_nl_proj)
          i=0
          write(*,*) 'Grid size: ',grid_size
          do l=0,pseudo(i_species)%lmax
             do j=1,nproj_per_l
                i=i+1
                call rad_alloc(pseudo(i_species)%pjnl(i),grid_size)
                pseudo(i_species)%pjnl_l(i) = l
                pseudo(i_species)%pjnl_n(i) = i
                write(*,*) 'n,l: ',pseudo(i_species)%pjnl_n(i),pseudo(i_species)%pjnl_l(i)
             end do
          end do
          ! We need to read the energies later - they're after all the tables
          if(str(19:20).eq.'/>') write(*,*) 'Found header end'
       end if
       write(*,*) 'Done header'
       ! Mesh: PP_MESH
       str = get_upf_line(lun)
       if (str(1:9) .eq. '<PP_MESH>') then
          str = get_upf_line(lun) ! We already have mesh_size
          str = get_upf_line(lun)
          read(str,*) rlast,r
          delta = r-rlast
          pseudo(i_species)%pjnl(1:n_nl_proj)%delta = delta
          write(*,*) 'Delta: ',delta
          str = get_upf_line(lun)
          do while(str(1:10).ne.'</PP_MESH>')
             str = get_upf_line(lun)
          end do
          write(*,*) 'Done mesh'
       end if
       pseudo(i_species)%pjnl(:)%delta = delta
       ! Local
       str = get_upf_line(lun)
       if (str(1:9) .eq. '<PP_LOCAL') then
          call rad_alloc(pseudo(i_species)%vlocal,grid_size)
          write(*,*) str(40:50)
          read(str(45:45),*) j
          write(*,*) 'Columns: ',j
          lines = (grid_size/j)
          rem = grid_size - lines*j
          write(*,*) 'Lines, rem: ',lines, rem
          st = 1
          do while(st<lines*j)
             !write(*,*) 'Position: ',st, st+j-1
             str = get_upf_line(lun)
             !write(*,*) 'Line: ',str
             read(str,*) (pseudo(i_species)%vlocal%f(i),i=st,st+j-1)
             st = st + j
          end do
          str = get_upf_line(lun)
          read(str,*) (pseudo(i_species)%vlocal%f(i),i=st,st+rem-1)
          write(*,*) 'Final number: ',st+rem-1,pseudo(i_species)%vlocal%f(st+rem-1)
          str = get_upf_line(lun)
          pseudo(i_species)%vlocal%delta = delta
          pseudo(i_species)%vlocal%cutoff = delta*grid_size
          !call rad_alloc( pseudo(i_species)%vna, grid_size )
          !pseudo(i_species)%vna%delta = delta
          !pseudo(i_species)%vna%cutoff = delta*grid_size
          if(str(1:11).ne.'</PP_LOCAL>') call cq_abort('Error in local')
       end if
       write(*,*) 'Starting non-local'
       ! Non-local
       str = get_upf_line(lun)
       if (str(1:12) .eq. '<PP_NONLOCAL') then
          do i_proj = 1,n_nl_proj
             ! Loop over projectors ? 
             str = get_upf_line(lun)
             if(str(1:8).eq.'<PP_BETA') then
                if(str(9:9).eq.'.') then
                   write(*,*) 'Projector index and read: ',i_proj,str(10:10)
                   ! I need total size and columns here
                   str = get_upf_line(lun)
                   do while(str(1:4).ne.'size')
                      str = get_upf_line(lun)
                   end do
                   read(str(7:10),*) grid_size
                   write(*,*) 'Grid size ',grid_size
                   ! Columns
                   str = get_upf_line(lun)
                   do while(str(1:7).ne.'columns')
                      str = get_upf_line(lun)
                   end do
                   read(str(10:10),*) j
                   write(*,*) 'Columns: ',j
                   str = get_upf_line(lun)
                   do while(str(1:19).ne.'cutoff_radius_index')
                      str = get_upf_line(lun)
                   end do
                   read(str(22:25),*) i_last!pseudo(i_species)%pjnl(n_nl_proj+1:n_nl_proj+nproj(l))%cutoff
                   call rad_alloc(pseudo(i_species)%pjnl(i_proj),i_last)
                   write(*,*) 'Max index: ',i_last
                   str = get_upf_line(lun)
                   do while(str(1:13).ne.'cutoff_radius')
                      str = get_upf_line(lun)
                   end do
                   write(*,*) str(16:34)
                   read(str(16:34),*) pseudo(i_species)%pjnl(i_proj)%cutoff
                else
                end if
                ! Read the projector
                lines = (grid_size/j)
                rem = grid_size - lines*j
                write(*,*) 'Lines, rem: ',lines, rem
                st = 1
                do while(st<lines*j)
                   str = get_upf_line(lun)
                   if(st<i_last) read(str,*) (pseudo(i_species)%pjnl(i_proj)%f(i),i=st,st+j-1)
                   st = st + j
                end do
                ! Read final line
                str = get_upf_line(lun)
                ! Get closing tag
                str = get_upf_line(lun)
                if(str(1:9).ne.'</PP_BETA') call cq_abort('Error reading beta: '//str(1:9))
             end if
          end do ! Loop over projectors i_proj=1,n_nl_proj
          ! Read projectors
          str = get_upf_line(lun)
          if (str(1:7) .eq. '<PP_DIJ') then
             do i_proj = 1,n_nl_proj*n_nl_proj/4
                str = get_upf_line(lun)
                read(str,*) temp((i_proj-1)*4 + 1:(i_proj-1)*4 + 4)
             end do
             do i_proj = 1,n_nl_proj
                pseudo(i_species)%pjnl_ekb(i_proj) = temp(1+(i_proj-1)*(n_nl_proj+1))
                write(*,*) 'Proj, Energy: ',i_proj,temp(1+(i_proj-1)*(n_nl_proj+1)) 
             end do
             str = get_upf_line(lun)
             if (str(1:8) .ne. '</PP_DIJ') call cq_abort('Error reading DIJ close tag: '//str(1:8))
          end if
          str = get_upf_line(lun)
          if (str(1:13) .ne. '</PP_NONLOCAL') call cq_abort('Error reading NONLOCAL close tag: '//str(1:13))
       else
          call cq_abort("Error reading non-local pseudopotentials: "//str(1:12))
       end if ! if(PP_NON_LOCAL
       ! PSWFC (ignore - should be empty)
       str = get_upf_line(lun)
       if (str(1:9) .eq. '<PP_PSWFC') then
          str = get_upf_line(lun)
          do while(str(1:10).ne.'</PP_PSWFC')
             str = get_upf_line(lun)
          end do
          write(*,*) 'Skipped PP_PSWFC'
       end if
       ! NLCC
       if(pseudo(i_species)%flag_pcc) then
          str = get_upf_line(lun)
          if(str(1:8) .eq. '<PP_NLCC') then
             call rad_alloc(pseudo(i_species)%chpcc,grid_size)
             lines = (grid_size/j)
             rem = grid_size - lines*j
             write(*,*) 'Lines, rem: ',lines, rem
             st = 1
             do while(st<lines*j)
                !write(*,*) 'Position: ',st, st+j-1
                str = get_upf_line(lun)
                !write(*,*) 'Line: ',str
                read(str,*) (pseudo(i_species)%chpcc%f(i),i=st,st+j-1)
                st = st + j
             end do
             str = get_upf_line(lun)
             read(str,*) (pseudo(i_species)%chpcc%f(i),i=st,st+rem-1)
             write(*,*) 'Final number: ',st+rem-1,pseudo(i_species)%chpcc%f(st+rem-1)
             str = get_upf_line(lun)
             if(str(1:9).ne.'</PP_NLCC>') call cq_abort('Error in local')
          else
             call cq_abort('Failed to find NLCC tag in str: '//str(1:47))
          end if
       end if
       ! PP_RHOATOM
       str = get_upf_line(lun)
       if (str(1:11) .eq. '<PP_RHOATOM') then
          if(str(38:44).eq.'columns') then
             read(str(47:47),*) j
             write(*,*) 'Columns: ',j
             call rad_alloc(pseudo(i_species)%chlocal,grid_size) ! OK, this is ATOMIC density but read as local
             lines = (grid_size/j)
             rem = grid_size - lines*j
             write(*,*) 'Lines, rem: ',lines, rem
             st = 1
             do while(st<lines*j)
                !write(*,*) 'Position: ',st, st+j-1
                str = get_upf_line(lun)
                !write(*,*) 'Line: ',str
                read(str,*) (pseudo(i_species)%chlocal%f(i),i=st,st+j-1)
                st = st + j
             end do
             str = get_upf_line(lun)
             read(str,*) (pseudo(i_species)%chlocal%f(i),i=st,st+rem-1)
             write(*,*) 'Final number: ',st+rem-1,pseudo(i_species)%chlocal%f(st+rem-1)
             str = get_upf_line(lun)
             if(str(1:13).ne.'</PP_RHOATOM>') call cq_abort('Error in local')
          else
             call cq_abort('Failed to find columns tag in str: '//str(1:47))
          end if
       else
          write(*,*) 'Expecting PP_RHOATOM tag - not found'
          call cq_abort('Wrong tag: '//str(1:11))
       end if
       ! End UPF tag
       str = get_upf_line(lun)
       if (str(1:6) .eq. '</UPF>') then
          write(*,*) 'UPF file correctly closed'
       else
          write(*,*) 'Missing close UPF tag; found instead: ',str(1:6)
       end if
       call io_close(lun)
    end do ! i=1,n_species
    return
  end subroutine read_upf

  ! Check for comment markers
  function get_upf_line(lun)

    implicit none
   
    integer :: lun
    character(len=80) :: get_upf_line

    character(len=80) :: a

    read(lun,'(a)') a
    get_upf_line = adjustl(a)
    do while(get_upf_line(1:4).eq.'<!--')
       read(lun,'(a)') a
       get_upf_line = adjustl(a)
    end do
    return
  end function get_upf_line

  ! Simple sorting routine - uses heap sort for n>2
  subroutine lsort(a,n)

    use datatypes

    implicit none

    ! Passed variables
    integer :: n
    real(double), dimension(n) :: a

    ! Local variables
    real(double), dimension(n) :: b
    real(double) :: temp
    integer, dimension(n) :: indx
    integer :: i, j, m, ir, indxt
    
    if(n==1) then
       return
    else if(n==2) then
       if(a(1)>a(2)) then
          temp = a(2)
          a(2) = a(1)
          a(1) = temp
       end if
       return
    else
       b = a ! Copy array
       do i=1,n
          indx(i) = i
       end do
       m=1+n/2
       ir = n
       do while(.true.)
          if(m>1) then
             m=m-1
             indxt = indx(m)
             temp = a(indxt)
          else
             indxt = indx(ir)
             temp = a(indxt)
             indx(ir)=indx(1)
             ir=ir-1
             if(ir==1) then
                indx(1) = indxt
                ! Copy back
                do i=1,n
                   a(i) = b(indx(i))
                end do
                return
             end if
          end if
          i = m
          j = m+m
          do while(j<=ir)
             if(j<ir) then
                if(a(indx(j))<a(indx(j+1))) j=j+1
             end if
             if(temp<a(indx(j))) then
                indx(i) = indx(j)
                i=j
                j=j+j
             else
                j=ir+1
             end if
          end do
          indx(i) = indxt
       end do ! Infinite ?!
    end if
    return
  end subroutine lsort

  ! Simple sorting routine - uses heap sort for n>2
  ! This sorts large -> small
  subroutine lsortr(a,n)

    use datatypes

    implicit none

    ! Passed variables
    integer :: n
    real(double), dimension(n) :: a

    ! Local variables
    real(double), dimension(n) :: b
    real(double) :: temp
    integer, dimension(n) :: indx
    integer :: i, j, m, ir, indxt
    
    if(n==1) then
       return
    else if(n==2) then
       if(a(1)<a(2)) then
          temp = a(2)
          a(2) = a(1)
          a(1) = temp
       end if
       return
    else
       b = a ! Copy array
       do i=1,n
          indx(i) = i
       end do
       m=1+n/2
       ir = n
       do while(.true.)
          if(m>1) then
             m=m-1
             indxt = indx(m)
             temp = a(indxt)
          else
             indxt = indx(ir)
             temp = a(indxt)
             indx(ir)=indx(1)
             ir=ir-1
             if(ir==1) then
                indx(1) = indxt
                ! Copy back
                do i=1,n
                   a(i) = b(indx(n+1-i))
                end do
                return
             end if
          end if
          i = m
          j = m+m
          do while(j<=ir)
             if(j<ir) then
                if(a(indx(j))<a(indx(j+1))) j=j+1
             end if
             if(temp<a(indx(j))) then
                indx(i) = indx(j)
                i=j
                j=j+j
             else
                j=ir+1
             end if
          end do
          indx(i) = indxt
       end do ! Infinite ?!
    end if
    return
  end subroutine lsortr
end module read
