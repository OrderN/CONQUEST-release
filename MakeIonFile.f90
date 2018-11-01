! Read optimised norm-conserving Vanderbilt pseudopotentials from QuantumEspresso or ABINIT
! and write out as Conquest .ion files.  Also generate PAOs for the atoms.
program MakeIonFiles

  use datatypes
  use numbers
  use GenComms, ONLY: myid, numprocs, init_comms, cq_abort
  use memory_module,     only: init_reg_mem
  use read
  use species_module, ONLY: n_species
  use pseudopotential_common, ONLY: ABINIT, UPF, pseudo_file_format_oncv
  use write
  use schro
  use pseudo_tm_info
  use pao_info
  use mesh

  implicit none

  integer :: i, ell, en
  real(double), allocatable, dimension(:) :: psi,r
  real(double) :: zval, energy
  
  call init_comms(myid,numprocs)
  call init_reg_mem

  call write_banner
  ! Read input (from Conquest_input file)
  call read_input
  !! Read pseudopotential files
  !if(pseudo_file_format_oncv==ABINIT) then
  !   call read_abinit
  !else if(pseudo_file_format_oncv==UPF) then
  !   call read_upf
  !else
  !   call cq_abort("Unknown pseudopotential format: ",pseudo_file_format_oncv)
  !end if
  call read_hamann_input
  call read_vkb
  call set_pao_initial
  call make_paos
  call write_header(val)
  call write_paos(val)
  call write_pseudopotential

end program MakeIonFiles
