module read

  use datatypes
  use GenComms, ONLY: cq_abort
  use global_module, ONLY: iprint
  
  implicit none

  character(len=80) :: block_file
  character(len=80), allocatable, dimension(:) :: pseudo_file_name, vps_file_name, vkb_file_name

  ! Local parameters 
  integer :: basis_size
  integer, parameter :: minimal = 1
  integer, parameter :: small = 2
  integer, parameter :: medium = 3
  integer, parameter :: full = 4

  integer :: energy_units ! Local for reading; 1 is Ha, 2 is eV
  real(double) :: energy_semicore ! Threshold for semi-core states
contains

  ! Read Conquest_input file for parameters from simulation, output parameters and coordinates
  subroutine read_input

    use input_module
    use numbers
    use species_module, ONLY: n_species
    use pseudopotential_common, ONLY: ABINIT, UPF, pseudo_file_format_oncv
    use mesh, ONLY: mesh_type, hamann, siesta, alpha, beta, nmesh_reg, delta_r_reg
    use schro, ONLY: pseudo_type, flag_user_specified, flag_default_cutoffs, flag_run_debug, &
         pao_cutoff_energies, pao_cutoff_radii, pao_cutoff_default, deltaE_large_radius, deltaE_small_radius, &
         flag_plot_output, local_and_vkb, val, n_debug_run, l_debug_run, E_debug_run
    use pao_info, ONLY: paos
    use units, ONLY: HaToeV
    
    implicit none

    character(len=80) :: input_string
    character(len=4) :: ps_format
    integer :: i, j, k, maxl, n_paos
    character(len=30), dimension(:), allocatable :: species_label
    real(double) :: temp
    
    ! Load the Conquest_input files
    call load_input
    ! Now scan for parameters
    iprint = fdf_integer('IO.Iprint',0)
    n_species = fdf_integer('General.NumberOfSpecies',1)
    ps_format = fdf_string(4,'General.PseudoFileFormatONCV','psp8')
    if((ps_format(1:3).eq.'upf').OR.(ps_format(1:3).eq.'UPF')) then
       pseudo_file_format_oncv = UPF
    else if((ps_format(1:4).eq.'psp8').OR.(ps_format(1:4).eq.'PSP8')) then
       pseudo_file_format_oncv = ABINIT
    end if
    input_string = fdf_string(6,'PS.PseudopotentialType','hamann')
    if(input_string(1:6)=='hamann') then
       pseudo_type = hamann
    else if(input_string(1:6)=='siesta') then
       pseudo_type = siesta
    end if
    flag_plot_output = fdf_boolean('IO.PlotOutput',.false.)
    input_string = fdf_string(6,'Mesh.MeshType',input_string(1:6)) ! Take mesh-type default as semilocal type
    if(input_string(1:6)=='hamann') then
       mesh_type = hamann
       !write(*,*) 'Mesh type: ',hamann
    else if(input_string(1:6)=='siesta') then
       mesh_type = siesta
       !write(*,*) 'Mesh type: ',hamann
    else
       call cq_abort("Unrecognised mesh type: "//input_string(1:6))
    end if
    alpha = fdf_double("Mesh.Alpha",alpha)
    beta = fdf_double("Mesh.Beta",beta)
    nmesh_reg = fdf_integer("Mesh.NMeshRegular",500)
    delta_r_reg = fdf_double("Mesh.RegularSpacing",0.01_double)
    allocate(pseudo_file_name(n_species),vps_file_name(n_species),species_label(n_species),&
         paos(n_species),vkb_file_name(n_species))
    allocate(local_and_vkb(n_species),val(n_species))
    ! Loop over species, finding their labels in the file
    ! All other species-dependent data will then be in the block with that label
    if(fdf_block('SpeciesLabels')) then
       if(1+block_end-block_start<n_species) & 
            call cq_abort("Too few species in SpeciesLabels: ",&
            1+block_end-block_start,n_species)
       do i=1,n_species
          read (unit=input_array(block_start+i-1),fmt=*) j, species_label(j)
       end do
       call fdf_endblock
    else
       call cq_abort("You need to specify SpeciesLabels block")
    end if
    ! Is this a debug run ? 
    flag_run_debug = fdf_boolean('General.RunDebug',.false.)
    energy_semicore = fdf_double('General.SemicoreEnergy',-one)
    if(energy_semicore>zero) write(*,fmt='(4x,"Possible error: your semi-core threshold is positive ! ",f6.3)') energy_semicore
    ! Energy cutoffs
    allocate(deltaE_large_radius(n_species), deltaE_small_radius(n_species))
    ! Now read the species-dependent data
    do i=1,n_species
       if(fdf_block(species_label(i))) then
          pseudo_file_name(i) = fdf_string(80,'Atom.PseudopotentialFile',' ')
          vps_file_name(i) = fdf_string(80,'Atom.SemilocalFile',' ')
          vkb_file_name(i) = fdf_string(80,'Atom.VKBFile',' ')
          paos(i)%flag_perturb_polarise = fdf_boolean("Atom.Perturbative_Polarised",.false.)
          if(paos(i)%flag_perturb_polarise) then
             paos(i)%polarised_n = fdf_integer("Atom.PolarisedN",0)
             if(paos(i)%polarised_n==0) call cq_abort("Must specify polarised shell n when using perturbation")
             paos(i)%polarised_l = fdf_integer("Atom.PolarisedL",-1)
             if(paos(i)%polarised_l==-1) call cq_abort("Must specify polarised shell l when using perturbation")
             paos(i)%polarised_shell = 0
          end if
          input_string = fdf_string(80,'Atom.ZetaForm','split')
          if(leqi(input_string(1:3),'spl')) then
             paos(i)%flag_zetas = 1 ! Split-norm approach
          else if(leqi(input_string(1:3),'com')) then
             paos(i)%flag_zetas = 2 ! Compression approach
          else
             call cq_abort("Unrecognised zeta form flag "//input_string(1:8))
          end if
          paos(i)%width = fdf_double("Atom.WidthConfine",one)
          paos(i)%prefac = fdf_double("Atom.PrefacConfine",zero)
          input_string = fdf_string(8,"Atom.Cutoffs","default")
          if(leqi(input_string(1:2),'en')) then ! Take reasonable keyword
             paos(i)%flag_cutoff = pao_cutoff_energies
          else if(leqi(input_string(1:2),'ra')) then ! Again
             paos(i)%flag_cutoff = pao_cutoff_radii
          else if(leqi(input_string(1:7),'default')) then
             paos(i)%flag_cutoff = pao_cutoff_default
          else
             call cq_abort("Unrecognised atomic cutoff flag "//input_string(1:8))
          end if
          energy_units = 1
          input_string = fdf_string(2,"Atom.EnergyUnits","Ha")
          if(leqi(input_string(1:2),"eV")) then
             energy_units = 2
             deltaE_large_radius(i) = fdf_double("Atom.dE_large_radius",0.02_double)
             deltaE_small_radius(i) = fdf_double("Atom.dE_small_radius",two)
             deltaE_large_radius(i) = deltaE_large_radius(i) / HaToeV
             deltaE_small_radius(i) = deltaE_small_radius(i) / HaToeV
          else if(leqi(input_string(1:2),"Ha")) then
             deltaE_large_radius(i) = fdf_double("Atom.dE_large_radius",0.00073498_double)
             deltaE_small_radius(i) = fdf_double("Atom.dE_small_radius",0.073498_double)
          end if
          input_string = fdf_string(7,'Atom.BasisSize','none')
          if(leqi(input_string(1:7),'minimal')) then
             basis_size = minimal
          else if(leqi(input_string(1:5),'small')) then
             basis_size = small
          else if(leqi(input_string,'medium')) then
             basis_size = medium
          else if(leqi(input_string(1:4),'full').OR.leqi(input_string(1:4),'none').OR.leqi(input_string(1:5),'large')) then ! Default
             basis_size = full
          else
             call cq_abort("Error ! Unknown Atom.BasisSize specification: "//input_string(1:7))
          end if
          ! If debug, check for specific n/l/energy
          if(flag_run_debug) then
             n_debug_run = fdf_integer("Atom.Debug_N",0)
             l_debug_run = fdf_integer("Atom.Debug_L",-1)
             E_debug_run = fdf_double("Atom.Debug_E",100.0_double)
          end if
          input_string = fdf_string(80,'Atom.BasisBlock','none')
          if(.NOT.(leqi(input_string(1:4),'none'))) then !.AND.paos(i)%flag_cutoff==3)) then
             paos(i)%n_shells = fdf_integer("Atom.PAO_N_Shells",0)
             if(paos(i)%n_shells==0) call cq_abort("Number of PAO shells not specified in atom block")
          end if
          call fdf_endblock
       else
          call cq_abort("Can't find species block for label "//species_label(i))
       end if
       if(input_string(1:4)=='none') then ! Default radii and zetas
          write(*,fmt='(2x,"Using default PAO specification")')
          flag_default_cutoffs = .true.
          ! We allocate the memory for PAOs elsewhere - after we've read the PS and worked out the number of valence shells
       else
          flag_default_cutoffs = .false.
          if(fdf_block(input_string)) then
             flag_user_specified = .true.
             !if(1+block_end-block_start<n_species) & 
             !     call cq_abort("Too few species in SpeciesLabels: ",&
             !     1+block_end-block_start,n_species)
             allocate(paos(i)%nzeta(paos(i)%n_shells),paos(i)%l(paos(i)%n_shells), &
                  paos(i)%n(paos(i)%n_shells),paos(i)%npao(paos(i)%n_shells), &
                  paos(i)%has_semicore(paos(i)%n_shells),paos(i)%inner(paos(i)%n_shells))
             maxl = 0
             n_paos = 0
             if(paos(i)%flag_perturb_polarise) then
                do j=1,paos(i)%n_shells-1
                   read (unit=input_array(block_start+j-1),fmt=*) paos(i)%n(j), paos(i)%l(j), paos(i)%nzeta(j)
                   if(paos(i)%nzeta(j)>maxl) maxl = paos(i)%nzeta(j)
                   !write(*,*) '# '//input_array(block_start+j-1)
                   n_paos = n_paos + paos(i)%nzeta(j)
                   if(paos(i)%n(j)==paos(i)%polarised_n.AND.paos(i)%l(j)==paos(i)%polarised_l) then
                      paos(i)%polarised_shell = j
                      if(iprint>3) write(*,fmt='("Found shell to polarise: ",i2)') j
                   end if
                end do
                if(paos(i)%polarised_shell==0) call cq_abort("Can't find shell to polarise corresponding to n, l: ", &
                     paos(i)%polarised_n,paos(i)%polarised_l)
                if(iprint>2) write(*,fmt='("For polarisation, we will perturb shell with n=",i2," and l=",i2)') &
                     paos(i)%n(paos(i)%polarised_shell), paos(i)%l(paos(i)%polarised_shell)
                paos(i)%n(paos(i)%n_shells) = paos(i)%n(paos(i)%polarised_shell)
                paos(i)%l(paos(i)%n_shells) = paos(i)%l(paos(i)%polarised_shell)+1
                paos(i)%npao(paos(i)%n_shells) = paos(i)%l(paos(i)%polarised_shell)+1
                do j=paos(i)%n_shells-1,1,-1
                   if(paos(i)%l(j)==paos(i)%l(paos(i)%n_shells)) then ! Semi-core with this l
                      paos(i)%npao(paos(i)%n_shells) = paos(i)%npao(paos(i)%n_shells) + 1
                      write(*,*) 'Found semi-core state ! ',paos(i)%l(j),paos(i)%npao(paos(i)%n_shells)
                      exit
                   end if
                end do
                paos(i)%nzeta(paos(i)%n_shells) = paos(i)%nzeta(paos(i)%polarised_shell)
                n_paos = n_paos + paos(i)%nzeta(paos(i)%n_shells)
             else
                do j=1,paos(i)%n_shells
                   read (unit=input_array(block_start+j-1),fmt=*) paos(i)%n(j), paos(i)%l(j), paos(i)%nzeta(j)
                   if(paos(i)%nzeta(j)>maxl) maxl = paos(i)%nzeta(j)
                   !write(*,*) '# '//input_array(block_start+j-1)
                   n_paos = n_paos + paos(i)%nzeta(j)
                   if(paos(i)%flag_perturb_polarise) then
                      if(paos(i)%polarised_n == paos(i)%n(j).AND.paos(i)%polarised_l == paos(i)%l(j)) paos(i)%polarised_shell = j
                   end if
                end do
             end if
             allocate(paos(i)%cutoff(maxl,paos(i)%n_shells))
             paos(i)%cutoff = zero
             allocate(paos(i)%energy(maxl,paos(i)%n_shells))
             paos(i)%energy = zero
             allocate(paos(i)%psi(maxl,paos(i)%n_shells))
             allocate(paos(i)%psi_reg(maxl,paos(i)%n_shells))
             !if(paos(i)%flag_perturb_polarise) then
             !   allocate(paos(i)%pol(maxl))
             !   allocate(paos(i)%pol_reg(maxl))
             !   n_paos = n_paos + paos(i)%nzeta(paos(i)%n_shells)
             !end if
             paos(i)%total_paos = n_paos
             if(paos(i)%flag_cutoff==pao_cutoff_energies.OR.paos(i)%flag_cutoff==pao_cutoff_default) then ! Energies
                if(paos(i)%flag_perturb_polarise) then
                   do j=1,paos(i)%n_shells-1
                      read (unit=input_array(block_start+paos(i)%n_shells-1+j-1),fmt=*) (paos(i)%energy(k,j),k=1,paos(i)%nzeta(j))
                      if(energy_units==2) paos(i)%energy(:,j) = paos(i)%energy(:,j)/ HaToeV
                      !write(*,*) '# '//input_array(block_start+paos(i)%n_shells-1+j-1)
                   end do
                   paos(i)%energy(:,paos(i)%n_shells) = paos(i)%energy(:,paos(i)%polarised_shell)
                else
                   do j=1,paos(i)%n_shells
                      read (unit=input_array(block_start+paos(i)%n_shells+j-1),fmt=*) (paos(i)%energy(k,j),k=1,paos(i)%nzeta(j))
                      if(energy_units==2) paos(i)%energy(:,j) = paos(i)%energy(:,j)/ HaToeV
                      !write(*,*) '# '//input_array(block_start+paos(i)%n_shells+j-1)
                   end do
                end if
                do j=1,paos(i)%n_shells
                   if(paos(i)%energy(1,j)>two) &
                        write(*,fmt='("Warning: energy shift of ",f6.3,"Ha is rather large.  Are these radii ?")') &
                        paos(i)%energy(1,j)
                end do
             else if(paos(i)%flag_cutoff==pao_cutoff_radii) then
                if(paos(i)%flag_perturb_polarise) then
                   do j=1,paos(i)%n_shells-1
                      read (unit=input_array(block_start+paos(i)%n_shells-1+j-1),fmt=*) (paos(i)%cutoff(k,j),k=1,paos(i)%nzeta(j))
                      !write(*,*) '# '//input_array(block_start+paos(i)%n_shells-1+j-1)
                   end do
                   paos(i)%cutoff(:,paos(i)%n_shells) = paos(i)%cutoff(:,paos(i)%polarised_shell)
                else
                   do j=1,paos(i)%n_shells
                      read (unit=input_array(block_start+paos(i)%n_shells+j-1),fmt=*) (paos(i)%cutoff(k,j),k=1,paos(i)%nzeta(j))
                      !write(*,*) '# '//input_array(block_start+paos(i)%n_shells+j-1)
                   end do
                end if
                do j=1,paos(i)%n_shells
                   if(paos(i)%cutoff(1,j)<two) &
                        write(*,fmt='("Warning: radius of ",f6.3,"a0 is rather small.  Are these energies ?")') paos(i)%cutoff(1,j)
                end do
             end if
             ! Sort the energies/cutoffs so that largest is first
             do j=1,paos(i)%n_shells
                if(paos(i)%flag_cutoff==pao_cutoff_radii) then
                   if(paos(i)%nzeta(j)==2) then
                      if(paos(i)%cutoff(1,j)<paos(i)%cutoff(2,j)) then
                         temp = paos(i)%cutoff(2,j)
                         paos(i)%cutoff(2,j) = paos(i)%cutoff(1,j)
                         paos(i)%cutoff(1,j) = temp
                      end if
                   else if(paos(i)%nzeta(j)==3) then
                      ! Not the most elegant, but simple; 1&2; 2&3; 1&2 again
                      if(paos(i)%cutoff(1,j)<paos(i)%cutoff(2,j)) then
                         temp = paos(i)%cutoff(2,j)
                         paos(i)%cutoff(2,j) = paos(i)%cutoff(1,j)
                         paos(i)%cutoff(1,j) = temp
                      end if
                      if(paos(i)%cutoff(2,j)<paos(i)%cutoff(3,j)) then
                         temp = paos(i)%cutoff(3,j)
                         paos(i)%cutoff(3,j) = paos(i)%cutoff(2,j)
                         paos(i)%cutoff(2,j) = temp
                      end if
                      if(paos(i)%cutoff(1,j)<paos(i)%cutoff(2,j)) then
                         temp = paos(i)%cutoff(2,j)
                         paos(i)%cutoff(2,j) = paos(i)%cutoff(1,j)
                         paos(i)%cutoff(1,j) = temp
                      end if
                   end if
                else ! Energies
                   if(paos(i)%nzeta(j)==2) then
                      if(paos(i)%energy(1,j)>paos(i)%energy(2,j)) then
                         temp = paos(i)%energy(2,j)
                         paos(i)%energy(2,j) = paos(i)%energy(1,j)
                         paos(i)%energy(1,j) = temp
                      end if
                   else if(paos(i)%nzeta(j)==3) then
                      ! Not the most elegant, but simple; 1&2; 2&3; 1&2 again
                      if(paos(i)%energy(1,j)>paos(i)%energy(2,j)) then
                         temp = paos(i)%energy(2,j)
                         paos(i)%energy(2,j) = paos(i)%energy(1,j)
                         paos(i)%energy(1,j) = temp
                      end if
                      if(paos(i)%energy(2,j)>paos(i)%energy(3,j)) then
                         temp = paos(i)%energy(3,j)
                         paos(i)%energy(3,j) = paos(i)%energy(2,j)
                         paos(i)%energy(2,j) = temp
                      end if
                      if(paos(i)%energy(1,j)>paos(i)%energy(2,j)) then
                         temp = paos(i)%energy(2,j)
                         paos(i)%energy(2,j) = paos(i)%energy(1,j)
                         paos(i)%energy(1,j) = temp
                      end if
                   end if
                end if
             end do
             call fdf_endblock
          else
             call cq_abort("Can't find species basis block "//input_string)
          end if
       end if
    end do
    call io_close(fdf_out)
    return
  end subroutine read_input

    ! Read semi-local potentials output by a DRB hack of Hamann's code
  ! Now read KB potentials
  subroutine read_vkb

    use numbers
    use species_module, ONLY: n_species
    use schro, ONLY: local_and_vkb, pseudo_type, val, flag_default_cutoffs, n_proj, pao_cutoff_energies, pao_cutoff_radii
    use pseudo_tm_info, ONLY: pseudo
    use input_module, ONLY: io_assign, io_close
    use mesh, ONLY: siesta
    use global_module, ONLY: flag_pcc_global
    use pao_info, ONLY: paos
    
    implicit none

    integer :: i_species, ios, lun, ngrid, ell, en, i, j, n_occ
    integer :: n_shells, n_nl_proj
    integer, dimension(0:4) :: count_func
    character(len=2) :: char_in
    character(len=80) :: line
    logical :: flag_core_done = .false.
    real(double) :: dummy, dummy2

    !allocate(local_and_vkb(n_species),val(n_species))
    do i_species = 1,n_species
       ! open file
       call io_assign(lun)
       open(unit=lun, file=vkb_file_name(i_species), status='old', iostat=ios)
       if ( ios > 0 ) call cq_abort('Error opening pseudopotential file: '//vkb_file_name(i_species))
       ! Allocate space
       read(lun,*) char_in,n_shells
       val(i_species)%n_shells = n_shells
       if(iprint>4) write(*,fmt='("Reading Hamann output, with ",i3," valence shells")') n_shells
       allocate(val(i_species)%n(n_shells),val(i_species)%npao(n_shells),val(i_species)%l(n_shells), &
            val(i_species)%occ(n_shells),val(i_species)%semicore(n_shells),val(i_species)%en_ps(n_shells), &
            val(i_species)%en_pao(n_shells),val(i_species)%has_semicore(n_shells),val(i_species)%inner(n_shells))
       count_func = 0 ! Count how many functions per l channel
       n_occ = 0
       val(i_species)%n = 0
       val(i_species)%npao = 0
       val(i_species)%l = 0
       val(i_species)%occ = zero
       val(i_species)%semicore = 0
       val(i_species)%inner = 0
       val(i_species)%en_ps = zero
       val(i_species)%en_pao = zero
       ! Read in table of shells, energies and occupancies and test for ordering
       do i = 1,n_shells
          read(lun,*) char_in,val(i_species)%n(i),val(i_species)%l(i),val(i_species)%occ(i),val(i_species)%en_ps(i)
          ! Simplistic check for correct ordering - I think that this is all that's needed though
          if(i>1) then
             do j=i,2,-1
                if(val(i_species)%en_ps(j)<val(i_species)%en_ps(j-1).AND.val(i_species)%occ(j)>RD_ERR) then ! Only re-order occupied
                   en = val(i_species)%n(j-1)
                   ell = val(i_species)%l(j-1)
                   dummy = val(i_species)%occ(j-1)
                   dummy2 = val(i_species)%en_ps(j-1)
                   val(i_species)%n(j-1) = val(i_species)%n(j)
                   val(i_species)%l(j-1) = val(i_species)%l(j)
                   val(i_species)%occ(j-1) = val(i_species)%occ(j)
                   val(i_species)%en_ps(j-1) = val(i_species)%en_ps(j)
                   val(i_species)%n(j) = en
                   val(i_species)%l(j) = ell
                   val(i_species)%occ(j) = dummy
                   val(i_species)%en_ps(j) = dummy2
                end if
             end do
          end if
       end do
       ! Check for semi-core states, count number of occupied shells and assign pseudo-n value (for nodes)
       if(iprint>3) then
          write(*,fmt='("Valence shells from pseudopotential calculation")')
          write(*,fmt='("  n  l    occ    energy")')
       end if
       do i = 1,n_shells
          if(iprint>3) write(*,fmt='(2i3,f7.2,f10.4)') &
               val(i_species)%n(i),val(i_species)%l(i),val(i_species)%occ(i),val(i_species)%en_ps(i)
          count_func(val(i_species)%l(i)) = count_func(val(i_species)%l(i)) + 1
          val(i_species)%has_semicore(i) = .false.
          ! Do we have semi-core states ? 
          if(count_func(val(i_species)%l(i))>1) then
             val(i_species)%semicore(i) = 0 ! This shell is the valence
             ! Find the corresponding inner shell
             do j=1,i-1
                if(val(i_species)%l(j)==val(i_species)%l(i)) then
                   val(i_species)%semicore(j) = 1
                   val(i_species)%npao(j) = val(i_species)%l(j)+1
                   val(i_species)%npao(i) = val(i_species)%l(i)+2
                   val(i_species)%has_semicore(i) = .true.
                   val(i_species)%inner(i) = j
                end if
             end do
          else if(val(i_species)%en_ps(i)<energy_semicore) then ! Make the value -1 a user input
             val(i_species)%semicore(i) = 1
             val(i_species)%npao(i) = val(i_species)%l(i)+1
          else
             val(i_species)%semicore(i) = 0
             val(i_species)%npao(i) = val(i_species)%l(i)+1
          end if
          if(val(i_species)%occ(i)>RD_ERR) n_occ = n_occ + 1
       end do
       if(iprint>3) write(*,fmt='(i2," valence shells, with ",i2," occupied")') n_shells, n_occ
       val(i_species)%n_occ = n_occ
       if(paos(i_species)%n_shells -2 >= val(i_species)%n_occ.AND.paos(i_species)%flag_cutoff == pao_cutoff_energies) &
            write(*,fmt='("Warning: cutoffs for unoccupied orbitals may not be reliable. Consider setting radii manually.")')
       if(.NOT.flag_default_cutoffs) then
          if(n_occ>paos(i_species)%n_shells) &
               call cq_abort("More valence shells in PP than in Atom.BasisBlock ",n_occ,paos(i_species)%n_shells)
       end if
       ! Read grid, charge, partial core, local potential, semilocal potentials
       read(lun,*) char_in, ngrid
       local_and_vkb(i_species)%ngrid = ngrid
       allocate(local_and_vkb(i_species)%rr(ngrid))
       local_and_vkb(i_species)%rr = zero
       allocate(local_and_vkb(i_species)%charge(ngrid))
       local_and_vkb(i_species)%charge = zero
       allocate(local_and_vkb(i_species)%local(ngrid))
       local_and_vkb(i_species)%local = zero
       n_nl_proj = maxval(n_proj)
       allocate(local_and_vkb(i_species)%projector(ngrid,n_nl_proj,0:pseudo(i_species)%lmax))
       local_and_vkb(i_species)%projector = zero
       allocate(local_and_vkb(i_species)%semilocal_potential(ngrid,0:pseudo(i_species)%lmax))
       local_and_vkb(i_species)%semilocal_potential = zero
       ! Read and store first block: local data
       if(flag_pcc_global) then
          allocate(local_and_vkb(i_species)%pcc(ngrid))
          local_and_vkb(i_species)%pcc = zero
          ! Read logarithmic mesh, atomic charge, partial core charge, local potential
          do i=1,ngrid
             read(lun,*) local_and_vkb(i_species)%rr(i),local_and_vkb(i_species)%charge(i), &
                  local_and_vkb(i_species)%pcc(i), local_and_vkb(i_species)%local(i), &
                  (local_and_vkb(i_species)%semilocal_potential(i,j),j=0,pseudo(i_species)%lmax)
          end do
       else
          ! Read logarithmic mesh, atomic charge, local potential - the core charge will be zero, so read as dummy
          do i=1,ngrid
             read(lun,*) local_and_vkb(i_species)%rr(i),local_and_vkb(i_species)%charge(i), &
                  dummy, local_and_vkb(i_species)%local(i), &
                  (local_and_vkb(i_species)%semilocal_potential(i,j),j=0,pseudo(i_species)%lmax)
          end do
       end if
       ! The VKB projectors are shorter-ranged than the full mesh, so store this size
       read(lun,*) line!char_in,i,j,ios ! We could this as a check on the number of projectors
       ! Read VKB projector coefficients
       n_nl_proj = 0
       if(iprint>4) write(*,fmt='("VKB projector energies"/,"  l  energies")')
       do ell=0,pseudo(i_species)%lmax
          read(lun,*) char_in,j,pseudo(i_species)%pjnl_ekb(n_nl_proj+1:n_nl_proj+n_proj(ell))
          if(iprint>4) write(*,fmt='(i3,4f10.5)') ell,pseudo(i_species)%pjnl_ekb(n_nl_proj+1:n_nl_proj+n_proj(ell))
          n_nl_proj = n_nl_proj + n_proj(ell)
       end do
       read(lun,*) char_in, ngrid
       local_and_vkb(i_species)%r_vkb = local_and_vkb(i_species)%rr(ngrid)
       local_and_vkb(i_species)%ngrid_vkb = ngrid
       ! And read VKB projectors
       do i=1,ngrid
          read(lun,*) ((local_and_vkb(i_species)%projector(i,j,ell),j=1,n_proj(ell)),ell=0,pseudo(i_species)%lmax)
       end do
       if(pseudo_type==siesta) then
          !write(*,*) '# Local_And_Vkb ps type siesta; scaling charge by 1/r2 and potential by 0.5'
          do i=1,ngrid
             local_and_vkb(i_species)%charge(i) = local_and_vkb(i_species)%charge(i)/ &
                  (local_and_vkb(i_species)%rr(i)*local_and_vkb(i_species)%rr(i))
             local_and_vkb(i_species)%local(i) = half*local_and_vkb(i_species)%local(i)
          end do
       end if
       call io_close(lun)
    end do
  end subroutine read_vkb

  subroutine read_hamann_input

    use global_module, ONLY: flag_pcc_global, flag_functional_type, functional_lda_pz81, functional_gga_pbe96, &
         functional_description
    use numbers
    use species_module, ONLY: n_species
    use input_module, ONLY: io_assign, io_close
    use pseudo_tm_info, ONLY: alloc_pseudo_info, pseudo
    use spline_module, ONLY: spline
    use schro, ONLY: n_proj, val
    
    implicit none

    integer :: i_species, lun, z, nc, nv, iexc, i_shell, en, ell, icmod, i, j, zeta, ios, n_nl_proj
    real(double) :: fill, zval
    character(len=80) :: a
    character(len=2) :: sym
    character(len=4) :: file_format

    allocate(pseudo(n_species))
    flag_pcc_global = .false.
    do i_species = 1, n_species
       ! Open appropriate file
       call io_assign(lun)
       open(unit=lun, file=pseudo_file_name(i_species), status='old', iostat=ios)
       if ( ios > 0 ) call cq_abort('Error opening Hamann input file: '//pseudo_file_name(i_species))
       pseudo(i_species)%filename = pseudo_file_name(i_species)
       ! Basic element information
       a = get_hamann_line(lun)
       read(a,*) sym,z,nc,nv,iexc,file_format
       write(*,fmt='(/"Information about pseudopotential for species: ",a2/)') sym
       pseudo(i_species)%z = z
       if(iexc==3) then
          val(i_species)%functional = functional_lda_pz81
       else if(iexc==4) then
          val(i_species)%functional = functional_gga_pbe96
       else
          call cq_abort("Error: unrecognised iexc value: ",iexc)
       end if
       ! DRB taken from initial_read.module.f90 in Conquest
       select case(val(i_species)%functional)
       case (functional_lda_pz81)
          functional_description = 'LDA PZ81'
       !case (functional_lda_gth96)
       !   functional_description = 'LDA GTH96'
       !case (functional_lda_pw92)
       !   functional_description = 'LSDA PW92'
       case (functional_gga_pbe96)
          functional_description = 'GGA PBE96'
       !case (functional_gga_pbe96_rev98)            ! This is PBE with the parameter correction
       !   functional_description = 'GGA revPBE98'   !   in Zhang & Yang, PRL 80:4, 890 (1998)
       !case (functional_gga_pbe96_r99)              ! This is PBE with the functional form redefinition
       !   functional_description = 'GGA RPBE99'     !   in Hammer et al., PRB 59:11, 7413-7421 (1999)
       !case (functional_gga_pbe96_wc)               ! Wu-Cohen nonempirical GGA functional
       !   functional_description = 'GGA WC'         !   in Wu and Cohen, PRB 73. 235116, (2006)
       !case (functional_hyb_pbe0)                   ! This is PB0E with the functional form redefinition
       !   functional_description = 'hyb PBE0'        
       !case default
       !   functional_description = 'LSDA PW92'
       end select
       write(*,fmt='("Using XC functional: ",a12)') functional_description
       write(*,fmt='("There are ",i2," core and ",i2," valence shells")') nc,nv
       a = get_hamann_line(lun)
       ! Read n, l, filling for core
       do i_shell = 1, nc
          read(a,*) en,ell,fill
          a = get_hamann_line(lun)
       end do
       ! Read n, l, filling for valence
       zval = zero
       do i_shell = 1, nv
          read(a,*) en,ell,fill
          zval = zval + fill
          a = get_hamann_line(lun)
       end do
       pseudo(i_species)%zval = zval
       write(*,fmt='("The atomic number is",i3,", with valence charge ",f4.1)') z,zval
       ! lmax
       read(a,*) pseudo(i_species)%lmax
       write(*,fmt='("Maximum angular momentum for pseudopotential is l=",i1)') pseudo(i_species)%lmax
       ! Projector radii etc
       a = get_hamann_line(lun)
       do ell = 0, pseudo(i_species)%lmax-1
          a = get_hamann_line(lun)
       end do
       ! Local potential
       a = get_hamann_line(lun)
       ! Numbers of projectors
       n_proj = 0
       n_nl_proj = 0
       a = get_hamann_line(lun)
       do ell = 0, pseudo(i_species)%lmax
          read(a,*) en,n_proj(en),fill
          n_nl_proj = n_nl_proj + n_proj(en)
          a = get_hamann_line(lun)
       end do
       write(*,fmt='("Total number of VKB projectors: ",i2)') n_nl_proj
       call alloc_pseudo_info(pseudo(i_species),n_nl_proj)
       i=0
       do ell=0,pseudo(i_species)%lmax
          zeta=0
          do j=1,n_proj(ell)
             i=i+1
             zeta = zeta+1
             !call rad_alloc(pseudo(i_species)%pjnl(i),grid_size)
             pseudo(i_species)%pjnl_l(i) = ell
             pseudo(i_species)%pjnl_n(i) = zeta
          end do
       end do
       ! Partial core corrections ?
       pseudo(i_species)%flag_pcc = .false.
       read(a,*) icmod
       if(icmod>0) then
          pseudo(i_species)%flag_pcc = .true.
          flag_pcc_global = .true.
          write(*,fmt='("This pseudopotential includes partial core corrections"/)') 
       else
          write(*,fmt='("This pseudopotential does not include partial core corrections"/)') 
       end if
       ! From here on we don't need the data, so the lines are commented out
       !! Comment line
       !read(lun, '(a)') a
       !! Log derivative analysis
       !read(lun, '(a)') a
       !read(lun, '(a)') a
       !! Comment line
       !read(lun, '(a)') a
       !! Output grid (Hamann code)
       call io_close(lun)
    end do
  end subroutine read_hamann_input
  
  ! Check for comment markers
  function get_hamann_line(lun)

    implicit none
   
    integer :: lun
    character(len=80) :: get_hamann_line

    character(len=80) :: a

    read(lun,'(a)') a
    get_hamann_line = adjustl(a)
    do while(get_hamann_line(1:1).eq.'#')
       read(lun,'(a)') a
       get_hamann_line = adjustl(a)
    end do
    return
  end function get_hamann_line

  subroutine set_pao_initial

    use numbers
    use species_module, ONLY: n_species
    use schro, ONLY: flag_default_cutoffs, val
    use pao_info, ONLY: paos
    use pseudo_tm_info, ONLY: pseudo
    use write, ONLY: pte

    implicit none

    integer :: i_species, j, i_highest_occ
    integer :: i_shell, ell, en, n_paos, n_shells, ell_hocc, max_zeta, n_pol_zeta
    integer, dimension(0:4) :: count_func

    ! If we are using defaults, then set up the structures
    if(flag_default_cutoffs) then 
       write(*,fmt='(/"Default basis chosen"/)')
       ! Set maxima based on chosen size
       if(basis_size==minimal) then ! One zeta for valence
          max_zeta = 1
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
       do i_species = 1,n_species
          write(*,fmt='(/"Species ",i2," is ",a2)') i_species, pte(pseudo(i_species)%z)
          count_func = 0
          write(*,fmt='("  n  l  nodes  zetas")')
          if(basis_size==minimal) then
             n_shells = val(i_species)%n_occ ! Occupied only
          else
             n_shells = val(i_species)%n_occ + 1 ! Occupied and polarisation
          end if
          paos(i_species)%n_shells = n_shells
          allocate(paos(i_species)%nzeta(n_shells),paos(i_species)%l(n_shells), &
               paos(i_species)%n(n_shells),paos(i_species)%npao(n_shells),paos(i_species)%l_no(n_shells),&
               paos(i_species)%has_semicore(n_shells),paos(i_species)%inner(n_shells))
          paos(i_species)%n = 0
          paos(i_species)%npao = 0
          paos(i_species)%nzeta = 0
          paos(i_species)%l = 0
          paos(i_species)%l_no = 0
          paos(i_species)%lmax = 0
          ! For default, we have SZ (shallow core) and TZTP so maximum number of zetas is 3
          allocate(paos(i_species)%cutoff(max_zeta,n_shells))
          paos(i_species)%cutoff = zero
          allocate(paos(i_species)%energy(max_zeta,n_shells))
          paos(i_species)%energy = zero
          allocate(paos(i_species)%psi(max_zeta,n_shells))
          allocate(paos(i_species)%psi_reg(max_zeta,n_shells))
          !if(paos(i_species)%flag_perturb_polarise) then
          !   allocate(paos(i_species)%pol(max_zeta))
          !   allocate(paos(i_species)%pol_reg(max_zeta))
          !end if
          n_paos = 0
          paos(i_species)%inner_shell = 0
          paos(i_species)%inner = 0
          paos(i_species)%has_semicore(:) = .false.
          do i_shell = 1,n_shells
             if(i_shell<=val(i_species)%n_occ) then ! We have n/l
                paos(i_species)%n(i_shell) = val(i_species)%n(i_shell)
                paos(i_species)%l(i_shell) = val(i_species)%l(i_shell)
                paos(i_species)%npao(i_shell) = val(i_species)%npao(i_shell)
                if(val(i_species)%semicore(i_shell)>0) then
                   paos(i_species)%nzeta(i_shell) = 1 
                   n_paos = n_paos + 1
                else
                   paos(i_species)%nzeta(i_shell) = max_zeta 
                   n_paos = n_paos + max_zeta
                end if
                do j=i_shell-1,1,-1
                   if(paos(i_species)%l(j)==paos(i_species)%l(i_shell)) then ! Semi-core with this l
                      paos(i_species)%has_semicore(i_shell) = .true.
                      paos(i_species)%inner(i_shell) = j
                   end if
                end do
             else ! Polarisation shell, defined as not being within the occupied valence states
                i_highest_occ = val(i_species)%n_occ
                if(paos(i_species)%l(i_highest_occ)<2) then ! We polarise this
                   ell = paos(i_species)%l(i_highest_occ)+1
                   paos(i_species)%npao(i_shell) = ell+1 ! Will need to change for perturb
                   if(paos(i_species)%flag_perturb_polarise)then
                      do j=1,paos(i_species)%n_shells -1
                         if(paos(i_species)%n(j)==paos(i_species)%polarised_n.AND. &
                              paos(i_species)%l(j)==paos(i_species)%polarised_l) then
                            paos(i_species)%polarised_shell = j
                            paos(i_species)%npao(i_shell) = paos(i_species)%npao(j) !+ 1
                            paos(i_species)%n(i_shell) = paos(i_species)%n(j)
                            paos(i_species)%l(i_shell) = paos(i_species)%l(j)+1
                            exit
                         end if
                      end do
                   else
                      do j=i_shell-1,1,-1
                      if(paos(i_species)%l(j)==ell) then ! Semi-core with this l
                         paos(i_species)%npao(i_shell) = paos(i_species)%npao(i_shell) + 1
                         paos(i_species)%has_semicore(i_shell) = .true.
                         paos(i_species)%inner(i_shell) = j
                         !if(paos(i_species)%flag_perturb_polarise) &
                         !     call cq_abort("We cannot support perturbative polarisation with semi-core states")
                         exit
                      end if
                   end do
                   paos(i_species)%n(i_shell) = paos(i_species)%n(i_highest_occ)
                   paos(i_species)%l(i_shell) = ell
                end if
                   paos(i_species)%nzeta(i_shell) = n_pol_zeta 
                   n_paos = n_paos + n_pol_zeta
                else if(i_highest_occ==1) then
                   call cq_abort("I can't polarise an atom with one shell and l=2")
                else ! Polarise the i_highest_occ - 1 shell
                   ell = paos(i_species)%l(i_highest_occ-1)+1
                   paos(i_species)%npao(i_shell) = ell + 1
                   if(paos(i_species)%flag_perturb_polarise)then
                      paos(i_species)%npao(i_shell) = ell
                   end if
                   do j=i_shell-1,1,-1
                      if(paos(i_species)%l(j)==ell) then ! Semi-core with this l
                         paos(i_species)%npao(i_shell) = paos(i_species)%npao(i_shell) + 1
                         exit
                      end if
                   end do
                   paos(i_species)%n(i_shell) = paos(i_species)%n(i_highest_occ-1)
                   if((paos(i_species)%n(i_shell) - ell -1) < 0) &
                        paos(i_species)%n(i_shell) = paos(i_species)%n(i_shell) + 1
                   paos(i_species)%l(i_shell) = ell
                   paos(i_species)%nzeta(i_shell) = n_pol_zeta 
                   n_paos = n_paos + n_pol_zeta
                end if
                !! Work out l
                !if(paos(i_species)%lmax<2) then
                !   ell = paos(i_species)%lmax+1
                !   paos(i_species)%npao(i_shell) = ell+1
                !   !paos(i_species)%npao(i_shell) = paos(i_species)%npao(i_shell-1)+1
                !   paos(i_species)%l(i_shell) = ell!val(i_species)%l(i_shell-1)+1
                !   paos(i_species)%nzeta(i_shell) = n_pol_zeta 
                !   n_paos = n_paos + n_pol_zeta
                !else if(flag_lmax_semicore) then
                !   ell = paos(i_species)%l(val(i_species)%n_occ) + 1 ! Polarise the outermost orbital
                !   paos(i_species)%npao(i_shell) = ell+2 ! Orthogonal to semicore
                !   !paos(i_species)%npao(i_shell) = paos(i_species)%npao(i_shell-1)+1
                !   paos(i_species)%l(i_shell) = ell!val(i_species)%l(i_shell-1)+1
                !   paos(i_species)%nzeta(i_shell) = n_pol_zeta 
                !   n_paos = n_paos + n_pol_zeta
                !else! Polarise s
                !   ell = 1
                !   paos(i_species)%npao(i_shell) = ell+1
                !   !paos(i_species)%npao(i_shell) = paos(i_species)%npao(i_shell-1)+1
                !   paos(i_species)%l(i_shell) = ell!val(i_species)%l(i_shell-1)+1
                !   paos(i_species)%nzeta(i_shell) = n_pol_zeta 
                !   n_paos = n_paos + n_pol_zeta
                !end if
                !paos(i_species)%n(i_shell) = 0
                !! Work out n
                !do j = i_shell-1,1,-1
                !   if(paos(i_species)%l(j)==ell) paos(i_species)%npao(i_shell) = paos(i_species)%npao(i_shell)+1
                !   if(paos(i_species)%n(i_shell)==0.AND.paos(i_species)%l(j)==ell-1) &
                !        paos(i_species)%n(i_shell) = paos(i_species)%n(j)
                !end do
                if(paos(i_species)%l(i_shell)>pseudo(i_species)%lmax) then ! We will need to use perturbation
                   write(*,fmt='("Warning ! l for polarisation is greater than lmax for pseudopotential.  Using perturbation.")')
                   paos(i_species)%flag_perturb_polarise = .true.
                   do j=i_shell-1,1,-1
                      !if(paos(i_species)%l(j)==pseudo(i_species)%lmax) then
                      if(paos(i_species)%l(j)<2) then
                         paos(i_species)%polarised_shell = j
                         exit
                      end if
                   end do
                   paos(i_species)%polarised_l = paos(i_species)%l(paos(i_species)%polarised_shell)
                   paos(i_species)%polarised_n = paos(i_species)%n(paos(i_species)%polarised_shell)
                   paos(i_species)%n(i_shell) = paos(i_species)%n(paos(i_species)%polarised_shell)
                   ell = paos(i_species)%l(paos(i_species)%polarised_shell)+1
                   paos(i_species)%l(i_shell) = ell
                   if((paos(i_species)%n(i_shell) - paos(i_species)%l(i_shell) -1) < 0) &
                        paos(i_species)%n(i_shell) = paos(i_species)%n(i_shell) + 1
                   !paos(i_species)%nzeta(i_shell) = paos(i_species)%nzeta(paos(i_species)%polarised_shell)
                   paos(i_species)%npao(i_shell) = paos(i_species)%npao(paos(i_species)%polarised_shell)
                   do j=i_shell-1,1,-1
                      if(paos(i_species)%l(j)==ell) then ! Semi-core with this l
                         paos(i_species)%npao(i_shell) = paos(i_species)%npao(i_shell) + 1
                         exit
                      end if
                   end do
                end if
             end if
             if(paos(i_species)%flag_perturb_polarise.AND.i_shell==n_shells) then
                if((paos(i_species)%l(i_shell)+1)>paos(i_species)%lmax) paos(i_species)%lmax=(paos(i_species)%l(i_shell)+1)
                count_func(paos(i_species)%l(i_shell)+1) = count_func(paos(i_species)%l(i_shell)) + 1
             else
                if(paos(i_species)%l(i_shell)>paos(i_species)%lmax) paos(i_species)%lmax=paos(i_species)%l(i_shell)
                count_func(paos(i_species)%l(i_shell)) = count_func(paos(i_species)%l(i_shell)) + 1
             end if
             paos(i_species)%l_no(i_shell) = count_func(paos(i_species)%l(i_shell))
             if(i_shell<=val(i_species)%n_occ) then
                if(val(i_species)%semicore(i_shell)>0) then
                   write(*,fmt='(2i3,2i7, "  semi-core state")') paos(i_species)%n(i_shell), paos(i_species)%l(i_shell),&
                        paos(i_species)%npao(i_shell) - paos(i_species)%l(i_shell) - 1, &
                        paos(i_species)%nzeta(i_shell)
                else
                   write(*,fmt='(2i3,2i7)') paos(i_species)%n(i_shell), paos(i_species)%l(i_shell),&
                        paos(i_species)%npao(i_shell) - paos(i_species)%l(i_shell) - 1, &
                        paos(i_species)%nzeta(i_shell)
                end if
             else
                if(paos(i_species)%flag_perturb_polarise) then
                   en = paos(i_species)%n(i_shell)
                   ell = paos(i_species)%l(i_shell)
                   !if(en-ell-1<0) en = en+1
                   write(*,fmt='(2i3,2i7)') en, ell,&
                        paos(i_species)%npao(i_shell) - paos(i_species)%l(i_shell) , &
                        paos(i_species)%nzeta(i_shell)
                else
                   write(*,fmt='(2i3,2i7)') paos(i_species)%n(i_shell), paos(i_species)%l(i_shell),&
                        paos(i_species)%npao(i_shell) - paos(i_species)%l(i_shell) - 1, &
                        paos(i_species)%nzeta(i_shell)
                end if
             end if
          end do ! n_shells
          if(basis_size==minimal) then
             paos(i_species)%polarised_shell = 0
          else
             paos(i_species)%polarised_shell = n_shells
          end if
          if(count_func(paos(i_species)%l(n_shells))>1.AND.paos(i_species)%flag_perturb_polarise) &
               call cq_abort("You can't have perturbative polarisation functions with semi-core states")
          paos(i_species)%total_paos = n_paos
       end do ! n_species
    else ! Check that the user specification is consistent
       write(*,fmt='("User-specified basis")')
       do i_species = 1,n_species
          count_func = 0
          write(*,fmt='("  n  l  nodes  zetas")')
          if(paos(i_species)%n_shells<val(i_species)%n_occ) then
             call cq_abort("Not enough shells specified for valence electrons: ",paos(i_species)%n_shells,val(i_species)%n_occ)
          end if
          do i_shell = 1,val(i_species)%n_occ
             if(paos(i_species)%n(i_shell) /= val(i_species)%n(i_shell)) call cq_abort("Mismatch in valence n value: ", &
                  paos(i_species)%n(i_shell),val(i_species)%n(i_shell))
             if(paos(i_species)%l(i_shell) /= val(i_species)%l(i_shell)) call cq_abort("Mismatch in valence n value: ", &
                  paos(i_species)%l(i_shell),val(i_species)%l(i_shell))
          end do ! i_shell = 1,n_shells
          ! Create npao - used for number of nodes
          paos(i_species)%lmax = 0
          paos(i_species)%inner = 0
          paos(i_species)%has_semicore(:) = .false.
          do i_shell = 1,paos(i_species)%n_shells
             count_func(paos(i_species)%l(i_shell)) = count_func(paos(i_species)%l(i_shell))+1
             if(paos(i_species)%flag_perturb_polarise.AND.i_shell==paos(i_species)%n_shells) then
                if(paos(i_species)%l(i_shell)+1>paos(i_species)%lmax) paos(i_species)%lmax=paos(i_species)%l(i_shell)+1
                !paos(i_species)%npao(i_shell) = paos(i_species)%npao(paos(i_species)%polarised_shell)
                !do j=i_shell-1,1,-1
                !   if(paos(i_species)%l(j)==ell) then ! Semi-core with this l
                !      paos(i_species)%npao(i_shell) = paos(i_species)%npao(i_shell) + 1
                !      exit
                !   end if
                !end do
                en = paos(i_species)%n(i_shell)
                if(en<3) en = en+1
                write(*,fmt='(2i3,2i7," using perturbative polarisation")') en, paos(i_species)%l(i_shell), &
                     paos(i_species)%npao(i_shell) - paos(i_species)%l(i_shell) - 1+1, & ! Silly book-keeping
                     paos(i_species)%nzeta(i_shell)
                !write(*,fmt='(2i3,2i7,"  Shell being polarised")') paos(i_species)%n(i_shell), paos(i_species)%l(i_shell),&
                !     paos(i_species)%npao(i_shell) - paos(i_species)%l(i_shell) - 1, &
                !     paos(i_species)%nzeta(i_shell)
             else
                if(paos(i_species)%l(i_shell)>paos(i_species)%lmax) paos(i_species)%lmax=paos(i_species)%l(i_shell)
                paos(i_species)%npao(i_shell) = paos(i_species)%l(i_shell) + 1 + count_func(paos(i_species)%l(i_shell))-1
                write(*,fmt='(2i3,2i7)') paos(i_species)%n(i_shell), paos(i_species)%l(i_shell),&
                     paos(i_species)%npao(i_shell) - paos(i_species)%l(i_shell) - 1, &
                     paos(i_species)%nzeta(i_shell)
                do j=i_shell-1,1,-1
                   if(paos(i_species)%l(j)==paos(i_species)%l(i_shell)) then ! Semi-core with this l
                      paos(i_species)%has_semicore(i_shell) = .true.
                      paos(i_species)%inner(i_shell) = j
                   end if
                end do
             end if
             if(i_shell>val(i_species)%n_occ) then
                if(val(i_species)%n_shells>2) then
                   if(val(i_species)%l(i_shell-1)==2.AND.val(i_species)%l(i_shell-2)==0) then ! OK - we have d/s overlap and polarise s
                      paos(i_species)%inner_shell = i_shell-2
                   end if
                else
                   paos(i_species)%inner_shell = i_shell-1
                end if
                if(iprint>4) write(*,fmt='(4x,"Using shell ",i2," as polarised shell")') paos(i_species)%inner_shell
             end if
             if(paos(i_species)%l(i_shell)>pseudo(i_species)%lmax) then ! Potential problem
                if(paos(i_species)%flag_perturb_polarise.AND.paos(i_species)%l(i_shell)==pseudo(i_species)%lmax+1) then
                   write(*,fmt='("Maximum l value for pseudopotential is ",i2," but you have chosen PAOs with l=",i2)') &
                        pseudo(i_species)%lmax, paos(i_species)%l(i_shell)
                   write(*,fmt='("As you are using perturbative polarisation, this is OK")')
                else
                   call cq_abort("Max l from PS smaller than max l from PAOs; use perturbative polarisation if needed.", &
                        pseudo(i_species)%lmax, paos(i_species)%l(i_shell))
                end if
             end if
          end do ! i_shell = 1,n_shells
          if(paos(i_species)%n_shells>val(i_species)%n_occ.AND.paos(i_species)%flag_perturb_polarise) then
             write(*,fmt='("Using perturbative polarisation: perturbing solution for n=",i2," and l=",i2)') &
                  paos(i_species)%n(val(i_species)%n_occ+1), paos(i_species)%l(val(i_species)%n_occ+1)
             if(paos(i_species)%n_shells-val(i_species)%n_occ>1) &
                  write(*,fmt='("Using normal Schrodinger solver for further shells")')
          end if
       end do ! n_species
    end if
  end subroutine set_pao_initial

  ! These routines were written to read the output from the Hamann code before realising that (a) we needed the semi-local
  ! potentials and (b) that the logarithmic radial mesh gave a much more accurate result.  They're kept below for reference
  ! but aren't used.
  subroutine read_abinit

    use datatypes
    use global_module, ONLY: flag_pcc_global
    use numbers
    use species_module, ONLY: n_species
    use input_module, ONLY: io_assign, io_close
    use pseudo_tm_info, ONLY: alloc_pseudo_info, pseudo, rad_alloc
    use spline_module, ONLY: spline
    use schro, ONLY: n_proj
    
    implicit none

    character(len=80) :: a
    integer :: i_species, ios, lun, loc_lmax
    integer :: pseudo_format, pseudo_xc, lmax, l_local, grid_size
    integer :: i,j,l, n_nl_proj, i_proj, i_last, zeta
    real(double) :: zatom, zion, pcc_radius, pcc_flag, r, rlast, delta, yp1, ypn
    real(double), dimension(:,:), allocatable :: temp_array
    
    allocate(pseudo(n_species))
    do i_species = 1,n_species
       ! open file
       call io_assign(lun)
       open(unit=lun, file=pseudo_file_name(i_species), status='old', iostat=ios)
       if ( ios > 0 ) call cq_abort('Error opening pseudopotential file: '//pseudo_file_name(i_species))
       pseudo(i_species)%filename = pseudo_file_name(i_species)
       ! Read title line (and print ?)
       read (lun, '(a)') a
       write(*,fmt='(2x,"Title line for file: ",(a),/,4x,a80)') pseudo_file_name(i_species),a
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
       if((.NOT.flag_pcc_global).AND.pseudo(i_species)%flag_pcc) flag_pcc_global = .true.
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
       end if
       !call rad_alloc( pseudo(i_species)%vna, grid_size )
       !pseudo(i_species)%vna%delta = delta
       !pseudo(i_species)%vna%cutoff = grid_size*delta
       call io_close(lun)
    end do
    if (flag_pcc_global) &
         write (*,fmt='(10x,a)') "P.C.C. is taken into account."
    
  end subroutine read_abinit

  subroutine read_upf

    use datatypes
    use global_module, ONLY: flag_pcc_global
    use numbers
    use species_module, ONLY: n_species
    use input_module, ONLY: io_assign, io_close
    use pseudo_tm_info, ONLY: alloc_pseudo_info, pseudo, rad_alloc
    use spline_module, ONLY: spline
    use schro, ONLY: n_proj
    use write, ONLY: pte
    
    implicit none

    character(len=80) :: a, str
    integer :: i_species, ios, lun, lines, st, rem, i_last
    integer :: pseudo_format, pseudo_xc, lmax, l_local, grid_size, verlen
    integer :: i,j,l, n_nl_proj, i_proj, ele, nproj_per_l
    real(double) :: zatom, zion, pcc_radius, pcc_flag, r, rlast, delta, yp1, ypn, temp(100)

    allocate(pseudo(n_species))
    do i_species = 1,n_species
       ! open file
       call io_assign(lun)
       open(unit=lun, file=pseudo_file_name(i_species), status='old', iostat=ios)
       if ( ios > 0 ) call cq_abort('Error opening pseudopotential file: '//pseudo_file_name(i_species))
       pseudo(i_species)%filename = pseudo_file_name(i_species)
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
          flag_pcc_global = .false.
          str = get_upf_line(lun)
          do while(str(1:15).ne.'core_correction')
             str = get_upf_line(lun)
          end do
          if((str(18:18).eq.'T').or.(str(18:18).eq.'t')) then
             flag_pcc_global = .true.
             pseudo(i_species)%flag_pcc = .true.
          end if
          if (flag_pcc_global) &
               write (*,fmt='(10x,a)') "P.C.C. is taken into account."
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
       if(flag_pcc_global) then
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

!%%!   ! Read semi-local potentials output by a DRB hack of Hamann's code
!%%!   ! Now read KB potentials
!%%!   subroutine read_vps
!%%! 
!%%!     use numbers
!%%!     use species_module, ONLY: n_species
!%%!     use schro, ONLY: pseudo_type, val, flag_default_cutoffs
!%%!     use pseudo_tm_info, ONLY: pseudo
!%%!     use input_module, ONLY: io_assign, io_close
!%%!     use mesh, ONLY: siesta
!%%!     use global_module, ONLY: flag_pcc_global
!%%!     use pao_info, ONLY: paos
!%%!     
!%%!     implicit none
!%%! 
!%%!     integer :: i_species, ios, lun, ngrid, ell, i, n_occ, n_paos
!%%!     integer :: n_shells
!%%!     integer, dimension(0:4) :: count_func
!%%!     character(len=2) :: char_in
!%%!     logical :: flag_core_done = .false.
!%%!     real(double) :: dummy
!%%! 
!%%!     allocate(semilocal(n_species),val(n_species))
!%%!     do i_species = 1,n_species
!%%!        ! open file
!%%!        call io_assign(lun)
!%%!        open(unit=lun, file=vps_file_name(i_species), status='old', iostat=ios)
!%%!        if ( ios > 0 ) call cq_abort('Error opening pseudopotential file: '//vps_file_name(i_species))
!%%!        ! Allocate space
!%%!        read(lun,*) char_in,n_shells
!%%!        if(paos(i_species)%n_shells==0) then
!%%!           paos(i_species)%n_shells = n_shells
!%%!        else if(paos(i_species)%n_shells/=n_shells) then
!%%!           call cq_abort("Incompatibility of valence shells: ",paos(i_species)%n_shells,n_shells)
!%%!        end if
!%%!        if(flag_default_cutoffs) then
!%%!           allocate(paos(i_species)%nzeta(n_shells),paos(i_species)%l(n_shells),paos(i_species)%n(n_shells))
!%%!           ! For default, we have SZ (shallow core) and TZTP so maximum number of zetas is 3
!%%!           allocate(paos(i_species)%cutoff(3,n_shells))
!%%!           paos(i_species)%cutoff = zero
!%%!           allocate(paos(i_species)%energy(3,n_shells))
!%%!           paos(i_species)%energy = zero
!%%!           allocate(paos(i_species)%psi(3,n_shells))
!%%!           allocate(paos(i_species)%psi_reg(3,n_shells))
!%%!           if(paos(i_species)%flag_polarise) then
!%%!              allocate(paos(i_species)%pol(3))
!%%!              allocate(paos(i_species)%pol_reg(3))
!%%!           end if
!%%!        end if
!%%!        write(*,*) '# In VPS, number of shells: ',n_shells
!%%!        allocate(val(i_species)%n(n_shells),val(i_species)%l(n_shells), &
!%%!             val(i_species)%occ(n_shells),val(i_species)%semicore(n_shells),val(i_species)%en_ps(n_shells))
!%%!        count_func = 0 ! Count how many functions per l channel
!%%!        n_occ = 0
!%%!        do i = 1,n_shells
!%%!           read(lun,*) char_in,val(i_species)%n(i),val(i_species)%l(i),val(i_species)%occ(i),val(i_species)%en_ps(i)
!%%!           write(*,*) '# Valence : ',val(i_species)%n(i),val(i_species)%l(i),val(i_species)%occ(i),val(i_species)%en_ps(i)
!%%!           count_func(val(i_species)%l(i)) = count_func(val(i_species)%l(i)) + 1
!%%!           val(i_species)%semicore(i) = 0
!%%!           if(val(i_species)%occ(i)>RD_ERR) n_occ = n_occ + 1
!%%!           if(flag_default_cutoffs) then
!%%!              paos(i_species)%n(i) = val(i_species)%n(i)
!%%!              paos(i_species)%l(i) = val(i_species)%l(i)
!%%!              paos(i_species)%nzeta(i) = 3 ! Fix for semi-core later
!%%!           end if
!%%!        end do
!%%!        write(*,*) '# valence shells, occupied shells: ',n_shells, n_occ
!%%!        val(i_species)%n_occ = n_occ
!%%!        ! Now flag semi-core states
!%%!        if(flag_default_cutoffs) n_paos = 0
!%%!        do ell=0,4
!%%!           if(count_func(ell)>1) then ! We have semi-core
!%%!              do i = 1,n_shells ! Loop over valence shells
!%%!                 if(val(i_species)%l(i)==ell.AND.(.NOT.flag_core_done)) then ! If this is the first, flag it
!%%!                    val(i_species)%semicore(i) = 1
!%%!                    val(i_species)%n(i) = ell+1
!%%!                    flag_core_done = .true.
!%%!                    paos(i_species)%nzeta(i) = 1 ! Fix for semi-core
!%%!                    if(flag_default_cutoffs) n_paos = n_paos + 1
!%%!                 else
!%%!                    val(i_species)%n(i) = ell+2                   
!%%!                    if(flag_default_cutoffs) n_paos = n_paos + 3
!%%!                 end if
!%%!              end do
!%%!           else
!%%!              do i = 1,n_shells ! Loop over valence shells
!%%!                 if(val(i_species)%l(i)==ell) then ! If this is the first, flag it
!%%!                    val(i_species)%n(i) = ell+1
!%%!                    if(flag_default_cutoffs) n_paos = n_paos + 3
!%%!                 end if
!%%!              end do
!%%!           end if
!%%!        end do
!%%!        read(lun,*) ngrid
!%%!        semilocal(i_species)%ngrid = ngrid
!%%!        allocate(semilocal(i_species)%rr(ngrid))
!%%!        allocate(semilocal(i_species)%charge(ngrid))
!%%!        allocate(semilocal(i_species)%potential(ngrid,0:pseudo(i_species)%lmax))
!%%!        if(flag_pcc_global) then
!%%!           allocate(semilocal(i_species)%pcc(ngrid))
!%%!           do i=1,ngrid
!%%!              read(lun,*) semilocal(i_species)%rr(i),semilocal(i_species)%charge(i), semilocal(i_species)%pcc(i), &
!%%!                   (semilocal(i_species)%potential(i,ell),ell=0,pseudo(i_species)%lmax)
!%%!           end do
!%%!        else
!%%!           do i=1,ngrid
!%%!              read(lun,*) semilocal(i_species)%rr(i),semilocal(i_species)%charge(i), dummy, &
!%%!                   (semilocal(i_species)%potential(i,ell),ell=0,pseudo(i_species)%lmax)
!%%!           end do
!%%!        end if
!%%!        if(pseudo_type==siesta) then
!%%!           write(*,*) '# Semilocal ps type siesta; scaling charge by 1/r2 and potential by 0.5'
!%%!           do i=1,ngrid
!%%!              semilocal(i_species)%charge(i) = semilocal(i_species)%charge(i)/ &
!%%!                   (semilocal(i_species)%rr(i)*semilocal(i_species)%rr(i))
!%%!              semilocal(i_species)%potential(i,:) = half*semilocal(i_species)%potential(i,:)
!%%!           end do
!%%!        end if
!%%!        call io_close(lun)
!%%!     end do
!%%!   end subroutine read_vps
!%%! 
!%%!   ! Read semi-local potentials output by a DRB hack of Hamann's code
!%%!   ! Only reads the potentials (assumes that all the rest is done by read_vkb)
!%%!   subroutine read_semilocal
!%%! 
!%%!     use numbers
!%%!     use species_module, ONLY: n_species
!%%!     use schro, ONLY: semilocal, pseudo_type, val, flag_default_cutoffs
!%%!     use pseudo_tm_info, ONLY: pseudo
!%%!     use input_module, ONLY: io_assign, io_close
!%%!     use mesh, ONLY: siesta
!%%!     use global_module, ONLY: flag_pcc_global
!%%!     use pao_info, ONLY: paos
!%%!     
!%%!     implicit none
!%%! 
!%%!     integer :: i_species, ios, lun, ngrid, ell, i, n_occ, n_paos
!%%!     integer :: n_shells, en
!%%!     integer, dimension(0:4) :: count_func
!%%!     character(len=2) :: char_in
!%%!     logical :: flag_core_done = .false.
!%%!     real(double) :: dummy, dummy2
!%%! 
!%%!     allocate(semilocal(n_species))
!%%!     do i_species = 1,n_species
!%%!        ! open file
!%%!        call io_assign(lun)
!%%!        open(unit=lun, file=vps_file_name(i_species), status='old', iostat=ios)
!%%!        if ( ios > 0 ) call cq_abort('Error opening pseudopotential file: '//vps_file_name(i_species))
!%%!        ! Allocate space
!%%!        read(lun,*) char_in,n_shells
!%%!        do i = 1,n_shells
!%%!           read(lun,*) char_in,en, ell, dummy, dummy2
!%%!        end do
!%%!        read(lun,*) ngrid
!%%!        semilocal(i_species)%ngrid = ngrid
!%%!        allocate(semilocal(i_species)%rr(ngrid))
!%%!        allocate(semilocal(i_species)%charge(ngrid))
!%%!        allocate(semilocal(i_species)%potential(ngrid,0:pseudo(i_species)%lmax))
!%%!        if(flag_pcc_global) then
!%%!           allocate(semilocal(i_species)%pcc(ngrid))
!%%!           do i=1,ngrid
!%%!              read(lun,*) semilocal(i_species)%rr(i),semilocal(i_species)%charge(i), semilocal(i_species)%pcc(i), &
!%%!                   (semilocal(i_species)%potential(i,ell),ell=0,pseudo(i_species)%lmax)
!%%!           end do
!%%!        else
!%%!           do i=1,ngrid
!%%!              read(lun,*) semilocal(i_species)%rr(i),semilocal(i_species)%charge(i), dummy, &
!%%!                   (semilocal(i_species)%potential(i,ell),ell=0,pseudo(i_species)%lmax)
!%%!           end do
!%%!        end if
!%%!        if(pseudo_type==siesta) then
!%%!           write(*,*) '# Semilocal ps type siesta; scaling charge by 1/r2 and potential by 0.5'
!%%!           do i=1,ngrid
!%%!              semilocal(i_species)%charge(i) = semilocal(i_species)%charge(i)/ &
!%%!                   (semilocal(i_species)%rr(i)*semilocal(i_species)%rr(i))
!%%!              semilocal(i_species)%potential(i,:) = half*semilocal(i_species)%potential(i,:)
!%%!           end do
!%%!        end if
!%%!        call io_close(lun)
!%%!     end do
!%%!   end subroutine read_semilocal
  
end module read
