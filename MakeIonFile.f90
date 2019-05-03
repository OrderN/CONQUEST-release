! Read optimised norm-conserving Vanderbilt pseudopotentials from Don Hamann's code
! as provided by PseudoDojo and SG15 libraries (Abinit and QuantumEspresso compatible)
! and write out as Conquest .ion files.  Also generate PAOs for the atoms.
program MakeIonFiles

  ! Conquest modules
  ! Also used in other modules: global, timer, mpi, units, pao_format, spline, functions, input
  use datatypes
  use numbers
  use GenComms, ONLY: myid, numprocs, init_comms, cq_abort
  use memory_module,     only: init_reg_mem
  use species_module, ONLY: n_species
  ! Slightly changed Conquest modules
  ! Removed gcopy statements
  use pseudo_tm_info
  ! UPF and pseudo_file_format_oncv are new - are they needed ?
  use pseudopotential_common, ONLY: ABINIT, UPF, pseudo_file_format_oncv
  ! Local modules
  use pseudo_atom_info
  use read
  use write
  use schro
  use mesh
  use input_module, ONLY: input_lines

  implicit none

  integer :: i_species

  call init_comms(myid,numprocs)
  call init_reg_mem

  call write_banner
  ! Load Conquest_input, read general parameters
  call read_general_input
  ! Loop over species, read appropriate data and generate PAOs
  do i_species = 1, n_species
     !
     ! Inputs
     !
     call read_species_input(i_species)
     !
     ! Setup
     !
     call set_pao_initial(i_species)
     ! Make PAOs
     call make_paos(i_species)
     !
     ! Output
     !
     call write_header(i_species)
     call write_paos(i_species)
     call write_pseudopotential(i_species)
     call deallocate_val
     call deallocate_pao
     call deallocate_vkb(i_species)
  end do
end program MakeIonFiles
