! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module initial_read
! ------------------------------------------------------------------------------
! Code area 1: initialisation
! ------------------------------------------------------------------------------

!!****h* Conquest/initial_read *
!!  NAME
!!   initial_read
!!  PURPOSE
!!   Holds the initial read routines
!!  USES
!!   atoms, construct_module, datatypes, DiagModule, dimens, fdf, GenComms, global_module,
!!   group_module, io_module, maxima_module, numbers, parse, primary_module, 
!!   pseudopotential_data, species_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   10/05/2002
!!  MODIFICATION HISTORY
!!   31/05/2002 dave
!!    Bug fix (from TM) for read_and_write
!!   17/06/2002 dave
!!    Added read to find solution method - order N or exact diagonalisation
!!   15:55, 25/09/2002 mjg & drb 
!!    Added flag for initial charge from initial K
!!   10:41, 06/03/2003 drb 
!!    Added tags for linear mixing for self-consistency in read_input
!!   11:02, 24/03/2003 drb 
!!    Simplified read_and_write
!!   14:51, 2003/06/09 dave
!!    Added new way of reading atomic coordinates, independent of make_prt
!!   2006/11/01 17:01 dave
!!    Included MP k-mesh generation from VB
!!   11:40, 10/11/2006 Veronika
!!    Added option for reading pdb files
!!   2008/02/01 03:43 dave
!!    Changes to give output to file not stdout
!!   2008/05/23 ast
!!    Added timers
!!   2010/07/23 Lianheng
!!    Added flags for Methfessel-Paxton smearing and associated Ef search method
!!   2011/07/21 16:50 dave
!!    Changes for cDFT
!!  SOURCE
!!
module initial_read

  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_allocation

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
!!***

contains

! ------------------------------------------------------------------------------
! Subroutine read_and_write
! ------------------------------------------------------------------------------

!!****f* initial_read/read_and_write *
!!
!!  NAME 
!!   read_and_write
!!  USAGE
!! 
!!  PURPOSE
!!   Reads the input data and writes out the parameters for the run
!!  INPUTS
!!   Many and varied - mainly things to be read in
!!  USES
!!   datatypes, numbers, io_module, group_module, construct_module, maxima_module, global_module,
!!   primary_module, atoms, dimens, species_module, GenComms, pseudopotential_data
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   Late 1998, I think
!!  MODIFICATION HISTORY
!!   28/05/2001 dave
!!    Added ROBODoc header, indented, stripped subroutine calls and
!!    overall argument list
!!   08/06/2001 dave
!!    Added RCS Id and Log tags and GenComms for my_barrier and cq_abort
!!   13/06/2001 dave
!!    Changed call to set_dimensions and to use pseudopotential_data
!!   10/05/2002 dave
!!    Incorporated into initial_read module
!!   31/05/2002 dave
!!    Bug fix (TM) - changed np to ind_part in loop to build atom positions from make_prt.dat
!!   11:02, 24/03/2003 drb 
!!    Removed variables associated with old phi initialisation
!!   14:52, 2003/06/09 dave
!!    Removed atomic coordinate fiddling
!!   14:59, 02/05/2005 dave 
!!    Added calculation of number of electrons in cell
!!   2007/05/08 17:00 dave
!!    Added net charge on cell to number_of_bands for consistency
!!   2011/09/29 15:06 M. Arita
!!    Slight modification is added for DFT-D2
!!  TODO
!!   Improve calculation of number of bands 28/05/2001 dave
!!   Change input so that appropriate variables are taken from modules 10/05/2002 dave
!!  SOURCE
!!
  subroutine read_and_write(start, start_L,&
       inode, ionode, restart_file, vary_mu, mu, &
       find_chdens, read_phi, number_of_bands)

    use datatypes
    use numbers
    use io_module, ONLY: read_mult, read_atomic_positions, pdb_template, pdb_format, print_process_info
    use group_module, ONLY : parts, part_method, HILBERT, PYTHON
    use construct_module, ONLY : init_group, init_primary
    use maxima_module, ONLY: maxpartsproc, maxatomsproc
    use global_module, ONLY: id_glob,x_atom_cell,y_atom_cell, &
         z_atom_cell, numprocs, iprint_init, rcellx,rcelly,rcellz,flag_old_partitions, &
         ne_in_cell, ni_in_cell, area_moveatoms, io_lun, flag_only_dispersion
    use memory_module, ONLY: reg_alloc_mem, type_dbl
    use primary_module, ONLY: bundle, make_prim
    use dimens, ONLY: r_super_x, r_super_y, r_super_z, &
         x_grid, y_grid, z_grid, r_h, set_dimensions, find_grid
    use block_module, ONLY: in_block_x, in_block_y, in_block_z
    use species_module, ONLY: n_species, species, charge, non_local_species
    use GenComms, ONLY: my_barrier, cq_abort
    use pseudopotential_data, ONLY: non_local, read_pseudopotential
    use pseudopotential_common, ONLY: core_radius, pseudo_type, OLDPS, SIESTA, ABINIT
    use pseudo_tm_module, ONLY: init_pseudo_tm
    use input_module, ONLY: fdf_string
    use force_module, ONLY: tot_force
    use DiagModule, ONLY: diagon

    implicit none

    ! Passed variables
    logical :: vary_mu, start, start_L, read_phi
    logical :: find_chdens

    character(len=40) restart_file

    integer :: inode, ionode

    real(double) :: mu
    real(double) :: number_of_bands

    ! Local variables
    character(len=80) titles,def
    character(len=80) atom_coord_file
    character(len=80) part_coord_file

    integer :: i,ierr, nnd, np, ni, ind_part
    integer :: ncf_tmp
    real(double) :: ecore_tmp

    real(double) :: HNL_fac

    ! read input data: parameters for run
    find_chdens = .true.
    call read_input(start, start_L, titles, restart_file, vary_mu, mu, find_chdens, read_phi,HNL_fac)

    ! Initialise group data for partitions and read in partitions and atoms
    call my_barrier()
    def = ' '
    atom_coord_file = fdf_string(80,'IO.Coordinates',def)
    if ( pdb_format ) then
       pdb_template = fdf_string(80,'IO.PdbTemplate',atom_coord_file)
    else
       pdb_template = fdf_string(80,'IO.PdbTemplate',' ')
    end if    
    def = 'make_prt.dat'
    part_coord_file = fdf_string(80,'IO.Partitions',def)
    call read_atomic_positions(trim(atom_coord_file))
    if(iprint_init>4) call print_process_info()
    ! By now, we'll have unit cell sizes and grid cutoff
    call find_grid
    if(diagon) call readDiagInfo
    if(inode==ionode.AND.iprint_init>1) write(io_lun,fmt='(10x,"Partitioning method: ",i2)') part_method
    call read_mult(inode-1,parts,part_coord_file)
    !if(iprint_init>4) write(io_lun,fmt='(10x,"Proc: ",i4," done read_mult: ",2i5)') inode,maxatomsproc,maxpartsproc
    if(allocated(tot_force)) then
       write(io_lun,fmt='(10x,"WARNING! Proc: ",i4," tot_force already allocated: ",i7)')  size(tot_force)
       deallocate(tot_force)
    end if
    allocate(tot_force(3,ni_in_cell))
    call reg_alloc_mem(area_moveatoms,3*ni_in_cell,type_dbl)
    call my_barrier

    ! Create primary set for atoms: bundle of partitions
    call init_primary(bundle, maxatomsproc, maxpartsproc, .true.)
    call make_prim(parts, bundle, inode-1,id_glob,x_atom_cell, y_atom_cell, z_atom_cell, species)

    ! Skip the following statement if only calculating the dispersion
    if (flag_only_dispersion) return

    !if(iprint_init>4) write(io_lun,fmt='(10x,"Proc: ",i4," done primary")') inode
    ! Read pseudopotential data
    if(pseudo_type == OLDPS) then
       call read_pseudopotential( inode, ionode)
    elseif(pseudo_type == SIESTA.OR.pseudo_type==ABINIT) then
       !just for setting core_radius in pseudopotential_common
       ! allocation and deallocation will be performed.
       call init_pseudo_tm(ecore_tmp, ncf_tmp)
    else
       call cq_abort(' Pseudotype Error ', pseudo_type)
    endif
    !if(iprint_init>4) write(io_lun,fmt='(10x,"Proc: ",i4," done pseudo")') inode

    ! Make number of bands in an astonishingly crude way
    number_of_bands = half*ne_in_cell!zero 
    do i = 1, ni_in_cell
       number_of_bands = number_of_bands + half * charge(species(i))
       ne_in_cell = ne_in_cell + charge(species(i))
    end do

    ! Set up various lengths, volumes, reciprocals etc. for convenient use
    call set_dimensions(inode, ionode,HNL_fac, non_local, n_species, non_local_species, core_radius)

    ! write out some information on the run
    if(inode==ionode) call write_info(number_of_bands, titles, mu, vary_mu, find_chdens, read_phi,HNL_fac, numprocs )

    call my_barrier()

    return
  end subroutine read_and_write
!!***

! ------------------------------------------------------------------------------
! Subroutine read_input
! ------------------------------------------------------------------------------

!!****f* initial_read/read_input *
!!
!!  NAME 
!!   read_input
!!  USAGE
!! 
!!  PURPOSE
!!   reads the input file
!!  INPUTS
!! 
!! 
!!  USES
!!   datatypes, global_module, atoms, dimens, species_module, pseudopotential_data, GenComms,
!!   fdf, parse
!!  AUTHOR
!!   E.H.Hernandez
!!  CREATION DATE
!!   08/05/95
!!  MODIFICATION HISTORY
!!   29/10/97 by Dave Bowler to add chdens reading
!!   19/11/97 by Dave Bowler to add phi reading
!!   28/05/2001 dave
!!    Stripped argument list, added ROBODoc header, indented
!!   08/06/2001 dave
!!    Added RCS Id and Log tags and GenComms for my_barrier, cq_abort
!!   09/05/2002 drb
!!    Converted to fdf formatting
!!   10/05/2002 dave
!!    Incorporated into initial_read module
!!   29/05/2002 dave
!!    Added scan for SolutionMethod tag
!!   15:54, 25/09/2002 mjg & drb 
!!    Added scan for make_initial_charge_from_K (which sets find_chdens)
!!   16:31, 2003/02/03 dave
!!    Added scans for type of run (e.g. MD, CG etc)
!!   07:49, 2003/02/04 dave
!!    Added boolean variables to control run (vary blips, self-consistency)
!!   10:41, 06/03/2003 drb 
!!    Added tags for linear mixing in self-consistency
!!   10:09, 12/03/2003 drb 
!!    Added an end tag
!!   14:52, 2003/06/09 dave
!!    Added new way to read atomic coordinates and preconditioning flag
!!   08:31, 2003/10/01 dave
!!    Changed flag_vary_blips to flag_vary_basis
!!   13:15, 03/10/2003 drb 
!!    Bug fix - missing .
!!   2004/10/29 drb
!!    Added various new flags for self consistency
!!   14:59, 02/05/2005 dave 
!!    Added call to check for net charge - initialises ne_in_cell
!!   09:14, 11/05/2005 dave 
!!    Store max_L_iterations in global variable
!!   13:34, 2006/07/09 ast
!!    Added flag for functional
!!   2007/04/02 07:53 dave
!!    Changed defaults
!!   2007/05/18 15:18 dave
!!    Added nullify for bp to stop g95 crash
!!   2007/10/15 Veronika
!!    Added keyword MaxEfIter to stop infinite loop in findFermi
!!   2008/01/24 Veronika
!!    Changed the default of General.ManyProcessors to .true.
!!   12:18, 14/02/2008 drb 
!!    Added options for buffer around primary and covering sets
!!   2008/07/16 ast
!!    New keywords for timers
!!   2009/07/08 16:47 dave
!!    Added keyword for one-to-one PAO to SF assigment
!!   2009/07/23 ast
!!    More general routine to name timer files, so that they are not 
!!    to a maximum number of processes, as it was the case before,
!!    with only three figures (max. 999 processes)
!!   2010/06/18 lt
!!    Added control flags for Methfessel-Paxton smearing method
!!   2010/07/22 15.53 Lianheng
!!    Moved control flags associated with Diagonalisation method from 
!!    here to readDiagInfo() subroutine. Added Methfessel-Paxton 
!!    smearing and related control flags. Changed SC.MaxEfIter to 
!!    Diag.MaxEfIter (moved to readDiagInfo()) for greater consistency
!!   2011/07/21 16:50 dave
!!    Added input flags for cDFT
!!   2011/09/16 11:08 dave
!!    Bug fix to say when species block not read
!!   2011/11/17 10:36 dave
!!    Changes for new blip data
!!  TODO
!!   Fix reading of start flags (change to block ?) 10/05/2002 dave
!!   Fix rigid shift 10/05/2002 dave
!!  SOURCE
!!
  subroutine read_input( start, start_L, &
       titles, restart_file, vary_mu, mu, find_chdens, read_phi,HNL_fac)

    use datatypes
    use numbers
    use units
    use global_module, ONLY: iprint, flag_vary_basis, flag_self_consistent, flag_precondition_blips, &
         flag_fractional_atomic_coords, flag_old_partitions, ne_in_cell, max_L_iterations, flag_read_blocks, &
         runtype, restart_L, restart_rho, flag_basis_set, blips, PAOs, flag_test_forces, UseGemm, &
         flag_fractional_atomic_coords, flag_old_partitions, ne_in_cell, max_L_iterations, &
         flag_functional_type, functional_description, functional_lda_pz81, functional_lda_gth96, &
         functional_lda_pw92, functional_gga_pbe96, functional_gga_pbe96_rev98, functional_gga_pbe96_r99, &
         flag_reset_dens_on_atom_move, flag_continue_on_SC_fail, &
         iprint_init, iprint_mat, iprint_ops, iprint_DM, iprint_SC, iprint_minE, iprint_time, &
         iprint_MD, iprint_index, iprint_gen, iprint_pseudo, iprint_basis, iprint_intgn, area_general, &
         global_maxatomspart, load_balance, many_processors, flag_assign_blocks, io_lun, &
         flag_pulay_simpleStep, flag_Becke_weights, flag_Becke_atomic_radii, flag_global_tolerance, & 
         flag_mix_L_SC_min, flag_onsite_blip_ana, flag_read_velocity, flag_quench_MD, temp_ion, &
         numprocs, flag_dft_d2, flag_only_dispersion, flag_perform_cDFT
    use dimens, ONLY: r_super_x, r_super_y, r_super_z, GridCutoff, &
         n_grid_x, n_grid_y, n_grid_z, r_h, r_c, RadiusSupport, NonLocalFactor, InvSRange, min_blip_sp, &
         flag_buffer_old, AtomMove_buffer, r_dft_d2
    use block_module, ONLY: in_block_x, in_block_y, in_block_z, blocks_raster, blocks_hilbert
    use species_module, ONLY: species_label, charge, mass, n_species, &
         species, ps_file, ch_file, phi_file, nsf_species, nlpf_species, npao_species, &
         non_local_species, type_species
    use GenComms, ONLY: gcopy, my_barrier, cq_abort, inode, ionode
    use DiagModule, ONLY: diagon
    use energy, ONLY: flag_check_Diag
    use DMMin, ONLY: maxpulayDMM, LinTol_DMM, n_dumpL
!TM
    use pseudopotential_common, ONLY: pseudo_type, OLDPS, SIESTA, STATE, ABINIT, flag_angular_new
    use SelfCon, ONLY: A, flag_linear_mixing, EndLinearMixing, q0, q1, n_exact, maxitersSC, maxearlySC, &
        maxpulaySC, atomch_output, flag_Kerker, flag_wdmetric
    use atomic_density, ONLY: read_atomic_density_file, atomic_density_method
    use S_matrix_module, ONLY: InvSTolerance
    use blip, ONLY: blip_info, init_blip_flag, alpha, beta
    use maxima_module, ONLY: maxnsf
    use control, ONLY: MDn_steps, MDfreq, MDtimestep, MDcgtol, CGreset
    use ewald_module, ONLY: ewald_accuracy, flag_old_ewald
    use minimise, ONLY: UsePulay, n_L_iterations, n_support_iterations, L_tolerance, sc_tolerance, &
         energy_tolerance, expected_reduction
    use pao_format, ONLY: kcut, del_k
    use support_spec_format, ONLY : flag_paos_atoms_in_cell, read_option, symmetry_breaking, support_pao_file, &
         TestBasisGrads, TestTot, TestBoth, TestS, TestH, flag_one_to_one
    use read_pao_info, ONLY: pao_info_file, pao_norm_flag
    use read_support_spec, ONLY: support_spec_file, flag_read_support_spec
    use test_force_module, ONLY: flag_test_all_forces, flag_which_force, TF_direction, TF_atom_moved, TF_delta
    use io_module, ONLY: pdb_format, pdb_altloc, append_coords, pdb_output, banner, get_file_name
    use group_module, ONLY : part_method, HILBERT, PYTHON
    use energy, ONLY: flag_check_DFT
    use H_matrix_module, ONLY: locps_output, locps_choice
    use pao_minimisation, ONLY: InitStep_paomin
    use timer_module, ONLY: time_threshold,lun_tmr, TimingOn, TimersWriteOut
    use input_module!, ONLY: load_input, input_array, block_start, block_end, fd
    use cdft_data, ONLY: cDFT_Type, cDFT_MaxIterations, cDFT_NAtoms, cDFT_Target, cDFT_Tolerance, &
         cDFT_NumberAtomGroups, cDFT_AtomList, cDFT_BlockLabel, cDFT_Vc

    implicit none

    ! Passed variables
    logical :: vary_mu, find_chdens
    logical :: start, start_L, read_phi

    character(len=40) :: restart_file
    character(len=80) :: titles

    real(double) :: mu, HNL_fac

    ! Local variables
    !type(block), pointer :: bp      ! Pointer to a block read by fdf
    !type(parsed_line), pointer :: p ! Pointer to a line broken into tokens by fdf

    integer :: i, j, k, lun, stat
    character(len=20) :: def
    character(len=80) :: coordfile, timefile, timefileroot
    character(len=10) :: basis_string, part_mode, tmp2
    character(len=6)  :: method ! To find whether we diagonalise or use O(N)
    character(len=5)  :: ps_type !To find which pseudo we use
    character(len=8)  :: tmp
    logical :: new_format
    !logical, external :: leqi
    real(double) :: r_t

    logical :: flag_ghost, find_species

    ! Set defaults
    vary_mu = .true.
    read_phi = .false.
    find_chdens = .false.

    !*** WHO READS AND BROADCASTS ? ***!

    ! Test for existence of input file
    if(inode==ionode) then
       stat = 0
       open(unit=20,file='Conquest_input',iostat=stat,status='old')
       if(stat/=0) call cq_abort("We need Conquest_input to run !")
       close(20)
    end if
    ! Start fdf reading from the file Conquest_input
    call load_input
    ! Here we need to read output file name, and open it.
    if(fdf_boolean('IO.WriteOutToFile',.true.)) then
       if(inode==ionode) then
          coordfile = fdf_string(80,'IO.OutputFile',"Conquest_out")
          call io_assign(io_lun)
          open(unit=io_lun,file=coordfile,iostat=stat)
          if(stat/=0) call cq_abort("Failed to open Conquest output file",stat)
       end if
       call gcopy(io_lun)
    else
       io_lun = 6
    end if
    ! Where will the timers will be written?
    TimingOn = fdf_boolean('IO.TimingOn',.false.)                                                                       ! tmr_rmv001
    if(TimingOn) then                                                                                                   ! tmr_rmv001
       if(fdf_boolean('IO.TimeAllProcessors',.false.)) then                                                             ! tmr_rmv001
          TimersWriteOut = .true.                                                                                       ! tmr_rmv001
          timefileroot = fdf_string(80,'IO.TimeFileRoot',"time")                                                        ! tmr_rmv001
          call get_file_name(timefileroot,numprocs,inode,timefile)                                                      ! tmr_rmv001
          call io_assign(lun_tmr)                                                                                       ! tmr_rmv001
          open(unit=lun_tmr,file=timefile,iostat=stat)                                                                  ! tmr_rmv001
          if(stat/=0) call cq_abort("Failed to open time file",stat)                                                    ! tmr_rmv001
       else                                                                                                             ! tmr_rmv001
          TimersWriteOut = .false.                                                                                      ! tmr_rmv001
          if(inode==ionode) then                                                                                        ! tmr_rmv001
             TimersWriteOut = .true.                                                                                    ! tmr_rmv001
             if(fdf_boolean('IO.WriteTimeFile',.true.)) then                                                            ! tmr_rmv001
                timefileroot = fdf_string(80,'IO.TimeFileRoot',"time")                                                  ! tmr_rmv001
                call get_file_name(timefileroot,numprocs,inode,timefile)                                                ! tmr_rmv001
                call io_assign(lun_tmr)                                                                                 ! tmr_rmv001
                open(unit=lun_tmr,file=timefile,iostat=stat)                                                            ! tmr_rmv001
                if(stat/=0) call cq_abort("Failed to open time file",stat)                                              ! tmr_rmv001
             else                                                                                                       ! tmr_rmv001
                lun_tmr = io_lun                                                                                        ! tmr_rmv001
             end if ! Output to file                                                                                    ! tmr_rmv001
          end if ! ionode switch                                                                                        ! tmr_rmv001
       end if                                                                                                           ! tmr_rmv001
    end if                                                                                                              ! tmr_rmv001
    if(inode==ionode) call banner
    new_format = fdf_boolean('IO.NewFormat',.true.)
    if(new_format) then
       def = ' '
       iprint = fdf_integer('IO.Iprint',0)
       iprint_init   = fdf_integer('IO.Iprint_init',iprint)
       iprint_mat    = fdf_integer('IO.Iprint_mat',iprint)
       iprint_ops    = fdf_integer('IO.Iprint_ops',iprint)
       iprint_DM     = fdf_integer('IO.Iprint_DM',iprint)
       iprint_SC     = fdf_integer('IO.Iprint_SC',iprint)
       iprint_minE   = fdf_integer('IO.Iprint_minE',iprint)
       iprint_MD     = fdf_integer('IO.Iprint_MD',iprint)
       iprint_index  = fdf_integer('IO.Iprint_index',iprint)
       iprint_gen    = fdf_integer('IO.Iprint_gen',iprint)
       iprint_pseudo = fdf_integer('IO.Iprint_pseudo',iprint)
       iprint_basis  = fdf_integer('IO.Iprint_basis',iprint)
       iprint_intgn  = fdf_integer('IO.Iprint_intgn',iprint)
       iprint_time   = fdf_integer('IO.Iprint_time',iprint)
       locps_output = fdf_boolean('IO.LocalPotOutput', .false.)
       locps_choice = fdf_integer('IO.LocalPotChoice', 8)
       atomch_output = fdf_boolean('IO.AtomChargeOutput', .false.)
       tmp = fdf_string(8,'General.MemoryUnits','MB')
       if(leqi(tmp(1:2),'kB')) then
          m_units = kbytes
          mem_conv = kB
       else if(leqi(tmp(1:2),'MB')) then
          m_units = mbytes
          mem_conv = MB
       else if(leqi(tmp(1:2),'GB')) then
          m_units = gbytes
          mem_conv = GB
       end if
       time_threshold = fdf_double('General.TimeThreshold',0.001_double)
       ! Read run title
       titles = fdf_string(80,'IO.Title',def)
       ! Is this a restart run ? **NB NOT AVAILABLE RIGHT NOW**
       start = fdf_boolean('General.NewRun',.true.)
       if(start) then
          start_L = .true.
          read_option = .false.  ! N.B. Original had option for a blip model - add this later
       else ! This option loads data from a file - we need to add the option to search for this
          call cq_abort('read_input: you may not select restart for a run just now')
       endif
       init_blip_flag = fdf_string(10,'Basis.InitBlipFlag','pao')
       restart_L = fdf_boolean('General.LoadL',.false.)
       restart_rho = fdf_boolean('General.LoadRho',.false.)
       ! Is there a net charge on the cell ?
       ne_in_cell = fdf_double('General.NetCharge',zero)
       ! Read coordinates file 
       flag_fractional_atomic_coords = fdf_boolean('IO.FractionalAtomicCoords',.true.)
       call my_barrier()
       !blip_width = fdf_double('blip_width',zero)
       !support_grid_spacing = fdf_double('support_grid_spacing',zero)
       GridCutoff = fdf_double('Grid.GridCutoff',25.0_double) ! Default to 25 Ha or 50 Ry cutoff for grid
       ! Grid points
       n_grid_x = fdf_integer('Grid.PointsAlongX',0)
       n_grid_y = fdf_integer('Grid.PointsAlongY',0)
       n_grid_z = fdf_integer('Grid.PointsAlongZ',0)
       ! Number of different iterations - not well defined
       n_L_iterations = fdf_integer('DM.LVariations',50)
       max_L_iterations = n_L_iterations
       n_dumpL = fdf_integer('DM.n_dumpL',n_L_iterations+1)
       n_support_iterations = fdf_integer('minE.SupportVariations',10)
       ! Initial expected drop in energy
       expected_reduction = fdf_double('minE.ExpectedEnergyReduction',zero)
       ! Keep for fixed potential calculations
       mu = fdf_double('DM.mu',zero)
       if(fdf_defined('DM.ConstantMu')) then
          vary_mu = .false.
       else
          vary_mu = .true.
       end if
       ! Radii - again, astonishingly badly named
       r_h = fdf_double('hamiltonian_range',zero)
       r_c = fdf_double('DM.L_range',one)
       r_t = fdf_double('InvSRange',r_h)
       HNL_fac = fdf_double('non_local_factor',zero)
       ! Exponents for initial gaussian blips
       alpha = fdf_double('Basis.SGaussExponent',zero)
       beta = fdf_double('Basis.PGaussExponent',zero)
       ! Tolerance on minimisation
       energy_tolerance = fdf_double('minE.EnergyTolerance',1.0e-5_double)
       UsePulay = fdf_boolean('Basis.UsePulayForPAOs',.false.)
       ! Sizes of integration grid blocks
       in_block_x = fdf_integer('Grid.InBlockX',4)
       in_block_y = fdf_integer('Grid.InBlockY',4)
       in_block_z = fdf_integer('Grid.InBlockZ',4)
       ! Solution method - O(N) or diagonalisation ?
       method = fdf_string(6,'DM.SolutionMethod','ordern') ! Default is O(N)
       if(leqi(method,'diagon')) then
          diagon = .true.
          flag_check_Diag = .true.
       else
          diagon = .false.
          flag_check_Diag = .false.
       end if
       ! Read basis set
       basis_string = fdf_string(10,'Basis.BasisSet','PAOs')
       if(leqi(basis_string,'blips')) then
          flag_basis_set = blips
       else if(leqi(basis_string,'PAOs')) then
          flag_basis_set = PAOs
       end if
       read_option = fdf_boolean('Basis.LoadCoeffs',.false.)
       flag_onsite_blip_ana = fdf_boolean('Basis.OnsiteBlipsAnalytic',.true.)
       if(read_option.AND.flag_basis_set==blips) restart_file = fdf_string(40,'Basis.LoadBlipFile',' ')
       find_chdens = fdf_boolean('SC.MakeInitialChargeFromK',.false.)
       flag_Becke_weights = fdf_boolean('SC.BeckeWeights',.false.)
       flag_Becke_atomic_radii = fdf_boolean('SC.BeckeAtomicRadii',.false.)
       ! Number of species
       n_species = fdf_integer('General.NumberOfSpecies',1)
       call allocate_species_vars
       flag_angular_new = fdf_boolean('Basis.FlagNewAngular',.true.)
       ! Tweak DRB 2007/03/23: fix abinit as standard type
       ps_type = fdf_string(5,'General.PseudopotentialType','abinit') 
       ! Write out pseudopotential type
       if(leqi(ps_type,'siest')) then
          if(inode==ionode.AND.iprint_init>0) write(io_lun,fmt='(10x,"SIESTA pseudopotential will be used. ")')
          pseudo_type = SIESTA
       else if(leqi(ps_type,'plato').OR.leqi(ps_type,'abini')) then
          if(inode==ionode.AND.iprint_init>0) write(io_lun,fmt='(10x,"ABINIT pseudopotential will be used. ")')
          pseudo_type = ABINIT
       else
          if(inode==ionode.AND.iprint_init>0) write(io_lun,fmt='(10x,"OLD pseudopotential will be used. ")')
          pseudo_type = OLDPS
       endif
       if((.NOT.flag_angular_new).AND.(pseudo_type==SIESTA.OR.pseudo_type==ABINIT).AND.inode==ionode) then
          write(io_lun,fmt='(10x,"Setting FlagNewAngular to T for Siesta/Abinit pseudopotentials")') 
          flag_angular_new = .true.
       end if
       ! Read, using old-style fdf_block, the information about the different species
       !if(fdf_block('ChemicalSpeciesLabel',lun)) then
       ! flag_ghost=.false.
       !   do i=1,n_species
       !      read(lun,*) j,mass(i),species_label(i)
       !      if(mass(i) < -very_small) flag_ghost=.true.
       !      type_species(i)=i
       !   end do
       !end if
       if(fdf_block('ChemicalSpeciesLabel')) then
          flag_ghost=.false.
          if(1+block_end-block_start<n_species) & 
               call cq_abort("Too few species in ChemicalSpeciesLabel: ",1+block_end-block_start,n_species)
          do i=1,n_species
             read(unit=input_array(block_start+i-1),fmt=*) j,mass(i),species_label(i)
             if(mass(i) < -very_small) flag_ghost=.true.
             type_species(i)=i
          end do
          call fdf_endblock
       end if
       ! I (TM) now assume we will use type_species(j) = 0
       ! for vacancy sites.   Nov. 14, 2007 
       if(flag_ghost) then
          do i=1, n_species
             if(mass(i) < -very_small) then
              find_species=.false.
               do j=1,n_species
                if(j == i) cycle
                if(species_label(i) == species_label(j)) then
                 type_species(i)=-j; find_species=.true.
                 exit
                endif
               enddo
              if(.not.find_species) then 
                 type_species(i) = -i
              endif
             elseif(mass(i) > very_small) then
                 type_species(i) = i
             else
                if(inode==ionode) write(io_lun,*) ' There are chemical species which have zero mass ',i,mass(i),species_label(i)
                 type_species(i) = 0    !  Vacancy sites
             endif
          enddo
       endif
        
       ! Read charge, mass, pseudopotential and starting charge and blip models for the individual species
       maxnsf = 0
       min_blip_sp = 1.0e8_double
       do i=1,n_species
          charge(i) = zero
          nsf_species(i) = 0
          RadiusSupport(i) = r_h
          NonLocalFactor(i) = HNL_fac
          InvSRange(i) = r_t
          blip_info(i)%SupportGridSpacing = zero
          if(pseudo_type==SIESTA.OR.pseudo_type==ABINIT) non_local_species(i) = .true.
          ! This is new-style fdf_block
          !nullify(bp)
          !if(fdf_block(species_label(i),bp)) then
          !   do while(fdf_bline(bp,line)) ! While there are lines in the block
          if(fdf_block(species_label(i))) then
             charge(i) = fdf_double('Atom.ValenceCharge',zero)
             nsf_species(i) = fdf_integer('Atom.NumberOfSupports',0)
             RadiusSupport(i) = fdf_double('Atom.SupportFunctionRange',r_h)
             InvSRange = fdf_double('Atom.InvSRange',zero)
             if(InvSRange(i)<very_small) InvSRange(i) = RadiusSupport(i)
             NonLocalFactor(i) = fdf_double('Atom.NonLocalFactor',HNL_fac)
             if(NonLocalFactor(i)>one.OR.NonLocalFactor(i)<zero) then
                if(inode==ionode) &
                     write(io_lun,fmt='(10x,"Warning: Atom.NonLocalFactor must lie between 0.0 and 1.0: ",f9.5)') &
                     NonLocalFactor(i)
                if(NonLocalFactor(i)>one) NonLocalFactor(i) = one
                if(NonLocalFactor(i)<zero) NonLocalFactor(i) = zero
             end if
             blip_info(i)%SupportGridSpacing = fdf_double('Atom.SupportGridSpacing',zero)
             min_blip_sp = min(blip_info(i)%SupportGridSpacing,min_blip_sp)
             blip_info(i)%BlipWidth = four*blip_info(i)%SupportGridSpacing
             if(pseudo_type==OLDPS.and.fdf_defined('Atom.Pseudopotential')) then
                non_local_species(i) = fdf_defined('Atom.NonLocal')
                if(non_local_species(i)) then
                   ps_file(i) = fdf_string(50,'Atom.NonLocal',def)
                else
                   if(fdf_defined('Atom.Local')) then
                      ps_file(i) = fdf_string(50,'Atom.Local',def)
                   else
                      call cq_abort('read_input: no local/non-local specification for pseudopotential !')
                   end if
                end if
             end if
             call fdf_endblock
          else
             call cq_abort("Failure to read data for species_label: "//species_label(i))
          end if
          if(nsf_species(i)==0) call cq_abort("Number of supports not specified for species ",i)
          if(flag_basis_set==blips.AND.blip_info(i)%SupportGridSpacing<very_small) &
               call cq_abort("Error: for a blip basis set you must set SupportGridSpacing for all species")
          maxnsf = max(maxnsf,nsf_species(i))
          if(RadiusSupport(i)<very_small) &
               call cq_abort("Radius of support too small for species; increase SupportFunctionRange ",i)
       end do
       ! For ghost atoms
       if(flag_ghost) then
          do i=1, n_species
             if(type_species(i) < 0) then
                charge(i) = zero
             endif
          enddo
       endif
       !blip_width = support_grid_spacing * fdf_double('blip_width_over_support_grid_spacing',four)
       flag_global_tolerance = fdf_boolean('minE.GlobalTolerance',.true.)
       L_tolerance = fdf_double('minE.LTolerance',1.0e-7_double)
       sc_tolerance  = fdf_double('minE.SCTolerance',1.0e-6_double)
       maxpulayDMM = fdf_integer('DM.MaxPulay',5)
       LinTol_DMM = fdf_double('DM.LinTol',0.1_double)
       ! Find out what type of run we're doing
       runtype = fdf_string(20,'AtomMove.TypeOfRun','static')
       flag_buffer_old = fdf_boolean('AtomMove.OldBuffer',.false.)
       AtomMove_buffer = fdf_double('AtomMove.BufferSize',4.0_double)
       flag_pulay_simpleStep = fdf_boolean('AtomMove.PulaySimpleStep',.false.)
       CGreset = fdf_boolean('AtomMove.ResetCG',.false.)
       MDn_steps = fdf_integer('AtomMove.NumSteps',100)
       MDfreq = fdf_integer('AtomMove.OutputFreq',50)
       MDtimestep = fdf_double('AtomMove.Timestep',0.5_double)
       MDcgtol = fdf_double('AtomMove.MaxForceTol',0.0005_double)
       flag_vary_basis = fdf_boolean('minE.VaryBasis',.false.)
       if(.NOT.flag_vary_basis) then
          flag_precondition_blips = .false.
       else 
          flag_precondition_blips = fdf_boolean('minE.PreconditionBlips',.true.)
       end if
       InitStep_paomin = fdf_double('minE.InitStep_paomin',5.0_double)
       flag_self_consistent = fdf_boolean('minE.SelfConsistent',.true.)
       flag_mix_L_SC_min = fdf_boolean('minE.MixedLSelfConsistent',.false.)
       ! Tweak 2007/03/23 DRB Make Pulay mixing default
       flag_linear_mixing = fdf_boolean('SC.LinearMixingSC',.true.)
       A = fdf_double('SC.LinearMixingFactor',0.5_double)
       EndLinearMixing = fdf_double('SC.LinearMixingEnd',sc_tolerance)
       flag_Kerker = fdf_boolean('SC.KerkerPreCondition',.false.)
       q0 = fdf_double('SC.KerkerFactor',0.1_double)
       flag_wdmetric = fdf_boolean('SC.WaveDependentMetric',.false.)
       q1 = fdf_double('SC.MetricFactor',0.1_double)
       n_exact = fdf_integer('SC.LateStageReset',5)
       flag_reset_dens_on_atom_move = fdf_boolean('SC.ResetDensOnAtomMove',.false.)
       flag_continue_on_SC_fail = fdf_boolean('SC.ContinueOnSCFail',.false.)
       maxitersSC = fdf_integer('SC.MaxIters',50)
       maxearlySC = fdf_integer('SC.MaxEarly',3)
       maxpulaySC = fdf_integer('SC.MaxPulay',5)
       read_atomic_density_file = fdf_string(80,'SC.ReadAtomicDensityFile','read_atomic_density.dat')
       ! Read atomic density initialisation flag
       atomic_density_method = fdf_string(10,'SC.AtomicDensityFlag','pao')
       ! cDFT flags
       flag_perform_cDFT = fdf_boolean('cDFT.Perform_cDFT',.false.)
       if(flag_perform_cDFT) then
          if(.NOT.flag_Becke_weights) then
             flag_Becke_weights = .true.
             write(io_lun,fmt='(2x,"Warning: we require Becke  weights for cDFT ! Setting to true")')
          end if
          cDFT_Type = fdf_integer('cDFT.Type',2) ! Default to D-A charge diff
          cDFT_MaxIterations = fdf_integer('cDFT.MaxIterations',50)
          cDFT_Tolerance = fdf_double('cDFT.Tolerance',0.001_double)
          !cDFT_Target = fdf_double('cDFT.Target',1.0_double)
          cDFT_NumberAtomGroups = fdf_integer('cDFT.NumberAtomGroups',1)
          allocate(cDFT_NAtoms(cDFT_NumberAtomGroups),cDFT_Target(cDFT_NumberAtomGroups), &
               cDFT_AtomList(cDFT_NumberAtomGroups),cDFT_BlockLabel(cDFT_NumberAtomGroups))
          if(fdf_block('cDFT.AtomGroups')) then
             !write(*,*) input_array(block_start+i-1:block_start+i+10)
             do i=1,cDFT_NumberAtomGroups
                read(unit=input_array(block_start+i-1),fmt=*) j,cDFT_NAtoms(i),cDFT_Target(i),cDFT_BlockLabel(i)
             end do
             call fdf_endblock
          else
             call cq_abort('Must define block cDFT.AtomGroups')
          end if
          do i=1,cDFT_NumberAtomGroups
             if(fdf_block(cDFT_BlockLabel(i))) then
                allocate(cDFT_AtomList(i)%Numbers(cDFT_NAtoms(i)))
                do j=1,cDFT_NAtoms(i)
                   read(unit=input_array(block_start+j-1),fmt=*) cDFT_AtomList(i)%Numbers(j)
                enddo
             else
                call cq_abort('cDFT block not defined: '//cDFT_BlockLabel(i))
             end if
             call fdf_endblock
          end do
       end if
       InvSTolerance = fdf_double('DM.InvSTolerance',1e-2_double)
       basis_string = fdf_string(10,'General.BlockAssign','Hilbert')
       if(leqi(basis_string,'Raster')) then
          flag_assign_blocks = blocks_raster
       else if(leqi(basis_string,'Hilbert')) then
          flag_assign_blocks = blocks_hilbert
       end if
       flag_read_blocks = fdf_boolean('Grid.ReadBlocks',.false.)
       ewald_accuracy = fdf_double('General.EwaldAccuracy',1.0e-10_double) ! Default value of 10^-10 Ha per atom
       flag_old_ewald = fdf_boolean('General.FlagOldEwald',.false.)
       UseGemm = fdf_boolean('MM.UseGemm',.false.)
        if(flag_ghost) then
           if(inode==ionode) write(io_lun,*) ' As ghost atoms are included, UseGemm must be false.'
           UseGemm = .false.
        endif
       flag_check_DFT=fdf_boolean('General.CheckDFT',.false.)
       flag_read_velocity = fdf_boolean('AtomMove.ReadVelocity',.false.)
       flag_quench_MD = fdf_boolean('AtomMove.QuenchMD',.false.)
       temp_ion = fdf_double('AtomMove.IonTemperature',300.0_double)

       del_k = fdf_double('Basis.PaoKspaceOlGridspace',0.1_double)
       kcut = fdf_double('Basis.PaoKspaceOlCutoff', 1000.0_double)
       flag_paos_atoms_in_cell = fdf_boolean('Basis.PAOs_StoreAllAtomsInCell',.true.)
       flag_one_to_one = fdf_boolean('Basis.PAOs_OneToOne',.false.)
       symmetry_breaking = fdf_boolean('Basis.SymmetryBreaking',.false.)
       support_pao_file = fdf_string(80,'Basis.SupportPaoFile','supp_pao.dat')
       pao_info_file = fdf_string(80,'Basis.PaoInfoFile','pao.dat')
       pao_norm_flag = fdf_integer('Basis.PaoNormFlag',0)
       TestBasisGrads = fdf_boolean('Basis.TestBasisGradients',.false.)
       TestTot = fdf_boolean('Basis.TestBasisGradTot',.false.)
       TestBoth = fdf_boolean('Basis.TestBasisGradBoth',.false.)
       TestS = fdf_boolean('Basis.TestBasisGrad_S',.false.)
       TestH = fdf_boolean('Basis.TestBasisGrad_H',.false.)
       support_spec_file = fdf_string(80,'Basis.SupportSpecFile','support.dat')
       flag_read_support_spec = fdf_boolean('Basis.ReadSupportSpec',.false.)
       flag_test_forces = fdf_boolean('AtomMove.TestForces',.false.)
       flag_test_all_forces = fdf_boolean('AtomMove.TestAllForces',.true.)
       if(.NOT.flag_test_all_forces) then ! Test one force component
          flag_which_force = fdf_integer('AtomMove.TestSpecificForce',1)
          if(flag_which_force>8.OR.flag_which_force<0.AND.inode==ionode) &
               write(io_lun,fmt='(10x,"Warning ! TestSpecificForce must lie between 1 and 8: ",i3)') flag_which_force
       end if
       TF_direction = fdf_integer('AtomMove.TestForceDirection',1)
       if(TF_direction>3.OR.TF_direction<0.AND.inode==ionode) then
          write(io_lun,fmt='(10x,"Warning ! TestForceDirection must lie between 1 and 3: ",i3)') TF_direction
          TF_direction = 1
       end if
       TF_atom_moved = fdf_integer('AtomMove.TestForceAtom',1)
       TF_delta = fdf_double('AtomMove.TestForceDelta',0.00001_double)
       flag_functional_type = fdf_integer('General.FunctionalType', 1)   ! LDA PZ81
       select case(flag_functional_type)
       case (functional_lda_pz81)
          functional_description = 'LDA PZ81'
       case (functional_lda_gth96)
          functional_description = 'LDA GTH96'
       case (functional_lda_pw92)
          functional_description = 'LDA PW92'
       case (functional_gga_pbe96)
          functional_description = 'GGA PBE96'
       case (functional_gga_pbe96_rev98)            ! This is PBE with the parameter correction
          functional_description = 'GGA revPBE98'   !   in Zhang & Yang, PRL 80:4, 890 (1998)
       case (functional_gga_pbe96_r99)              ! This is PBE with the functional form redefinition
          functional_description = 'GGA RPBE99'     !   in Hammer et al., PRB 59:11, 7413-7421 (1999)
       case default
          functional_description = 'LDA PZ81'
       end select
       tmp = fdf_string(8,'General.EnergyUnits','Ha')
       if(leqi(tmp(1:2),'Ha')) then
          energy_units = har
          en_conv = one
       else if(leqi(tmp(1:2),'Ry')) then
          energy_units = ryd
          en_conv = HaToRy
       else if(leqi(tmp(1:2),'eV')) then
          energy_units = ev
          en_conv = HaToeV
       endif
       tmp = fdf_string(8,'General.DistanceUnits','bohr')
       if(leqi(tmp(1:2),'a0').OR.leqi(tmp(1:2),'bo')) then
          dist_units = bohr
          dist_conv = one
       else if(leqi(tmp(1:1),'A')) then
          dist_units = ang
          dist_conv = BohrToAng
       endif
       for_conv = en_conv/dist_conv
       pdb_format = fdf_boolean('IO.PdbIn',.false.)
       pdb_altloc = fdf_string(1,'IO.PdbAltLoc',' ')
       pdb_output = fdf_boolean('IO.PdbOut',.false.)
       part_mode = fdf_string(10,'General.PartitionMethod','Hilbert')
       if (leqi (part_mode(1:6),'File')) then
          part_method = PYTHON
       else if (leqi (part_mode(1:7), 'Hilbert')) then
          part_method = HILBERT
       else
          part_method = HILBERT
          if(inode==ionode) write(io_lun,*) 'WARNING: Unrecognised partitioning method'
       end if
       tmp2 = fdf_string(10,'General.LoadBalance','atoms')
       if (leqi (tmp2(1:10),'partitions')) then
          load_balance = 1
       else if (leqi (tmp2(1:5),'atoms')) then
          load_balance = 0
       else
          load_balance = 0
       end if
       many_processors = fdf_boolean('General.ManyProcessors',.true.)
       global_maxatomspart = fdf_integer('General.MaxAtomsPartition', 34)
       append_coords = fdf_boolean('AtomMove.AppendCoords',.true.)
       flag_dft_d2 = fdf_boolean('General.DFT_D2', .false.)                   ! for DFT-D2
       if (flag_dft_d2) r_dft_d2 = fdf_double('DFT-D2_range',23.0_double)     ! for DFT-D2
       flag_only_dispersion = fdf_boolean('General.only_Dispersion',.false.)  ! for DFT-D2
    else
       call cq_abort("Old-style CQ input no longer supported: please convert")
!%%!else
!%%!   def = ' '
!%%!   iprint = fdf_integer('iprint',0)
!%%!   iprint_init   = fdf_integer('iprint_init',iprint)
!%%!   iprint_mat    = fdf_integer('iprint_mat',iprint)
!%%!   iprint_ops    = fdf_integer('iprint_ops',iprint)
!%%!   iprint_DM     = fdf_integer('iprint_DM',iprint)
!%%!   iprint_SC     = fdf_integer('iprint_SC',iprint)
!%%!   iprint_minE   = fdf_integer('iprint_minE',iprint)
!%%!   iprint_MD     = fdf_integer('iprint_MD',iprint)
!%%!   iprint_index  = fdf_integer('iprint_index',iprint)
!%%!   iprint_gen    = fdf_integer('iprint_gen',iprint)
!%%!   iprint_pseudo = fdf_integer('iprint_pseudo',iprint)
!%%!   iprint_basis  = fdf_integer('iprint_basis',iprint)
!%%!   iprint_intgn  = fdf_integer('iprint_intgn',iprint)
!%%!   tmp = fdf_string('MemoryUnits','MB')
!%%!   if(leqi(tmp(1:2),'kB')) then
!%%!      m_units = kbytes
!%%!      mem_conv = kB
!%%!   else if(leqi(tmp(1:2),'MB')) then
!%%!      m_units = mbytes
!%%!      mem_conv = MB
!%%!   else if(leqi(tmp(1:2),'GB')) then
!%%!      m_units = gbytes
!%%!      mem_conv = GB
!%%!   end if
!%%!   ! Read run title
!%%!   titles = fdf_string('title',def)
!%%!   ! Is this a restart run ? **NB NOT AVAILABLE RIGHT NOW**
!%%!   start = .true.!fdf_defined('start')
!%%!   if(start) then
!%%!      start_L = .true.
!%%!      read_option = .false.  ! N.B. Original had option for a blip model - add this later
!%%!   else ! This option loads data from a file - we need to add the option to search for this
!%%!      call cq_abort('read_input: you may not select restart for a run just now')
!%%!   endif
!%%!   !read_option = fdf_boolean('StartBlip',.true.)
!%%!   !if(.NOT.start_blips) restart_file = fdf_string('LoadBlipFile',' ')
!%%!   init_blip_flag = fdf_string('init_blip_flag','pao')
!%%!   restart_L = fdf_boolean('RestartL',.false.)
!%%!   restart_rho = fdf_boolean('RestartRho',.false.)
!%%!   ! Is there a net charge on the cell ?
!%%!   ne_in_cell = fdf_double('NetCharge',zero)
!%%!   ! Read coordinates file 
!%%!   flag_fractional_atomic_coords = fdf_boolean('FractionalAtomicCoords',.true.)
!%%!   call my_barrier()
!%%!   !blip_width = fdf_double('blip_width',zero)
!%%!   !support_grid_spacing = fdf_double('support_grid_spacing',zero)
!%%!   GridCutoff = fdf_double('GridCutoff',20.0_double) ! Default to 20 Ha or 50 Ry cutoff for grid
!%%!   ! Grid points
!%%!   n_grid_x = fdf_integer('grid_points_along_x',0)
!%%!   n_grid_y = fdf_integer('grid_points_along_y',0)
!%%!   n_grid_z = fdf_integer('grid_points_along_z',0)
!%%!   ! Number of different iterations - not well defined
!%%!   n_L_iterations = fdf_integer('l_variations',50)
!%%!   max_L_iterations = n_L_iterations
!%%!   n_support_iterations = fdf_integer('support_variations',20)
!%%!   ! Initial expected drop in energy
!%%!   expected_reduction = fdf_double('expected_energy_reduction',zero)
!%%!   ! Keep for fixed potential calculations
!%%!   mu = fdf_double('mu',zero)
!%%!   if(fdf_defined('constant_mu')) then
!%%!      vary_mu = .false.
!%%!   else
!%%!      vary_mu = .true.
!%%!   end if
!%%!   ! Radii - again, astonishingly badly named
!%%!   r_h = fdf_double('hamiltonian_range',zero)
!%%!   r_c = fdf_double('L_range',one)
!%%!   r_t = fdf_double('InvSRange',r_h)
!%%!   HNL_fac = fdf_double('non_local_factor',zero)
!%%!   ! Exponents for initial gaussian blips
!%%!   alpha = fdf_double('s_gauss_exponent',zero)
!%%!   beta = fdf_double('p_gauss_exponent',zero)
!%%!   ! Tolerance on minimisation
!%%!   energy_tolerance = fdf_double('energy_tolerance',1.0e-5_double)
!%%!   UsePulay = fdf_boolean('UsePulayForPAOs',.false.)
!%%!   ! Sizes of integration grid blocks
!%%!   in_block_x = fdf_integer('in_block_x',4)
!%%!   in_block_y = fdf_integer('in_block_y',4)
!%%!   in_block_z = fdf_integer('in_block_z',4)
!%%!   ! Solution method - O(N) or diagonalisation ?
!%%!   method = fdf_string('SolutionMethod','ordern') ! Default is O(N)
!%%!   if(leqi(method,'diagon')) then
!%%!      diagon = .true.
!%%!   else
!%%!      diagon = .false.
!%%!   end if
!%%!   ! Read basis set
!%%!   basis_string = fdf_string('BasisSet','PAOs')
!%%!   if(leqi(basis_string,'blips')) then
!%%!      flag_basis_set = blips
!%%!   else if(leqi(basis_string,'PAOs')) then
!%%!      flag_basis_set = PAOs
!%%!   end if
!%%!   find_chdens = fdf_boolean('make_initial_charge_from_K',.false.)
!%%!   ! Number of species
!%%!   n_species = fdf_integer('NumberOfSpecies',1)
!%%!   call allocate_species_vars
!%%!   flag_angular_new = fdf_boolean('FlagNewAngular',.true.)
!%%!   ! Siesta's pseudopotential or not.  14/11/2002 T. Miyazaki
!%%!   ps_type = fdf_string('PseudopotentialType','abinit') 
!%%!   ! siesta's pseudo is not used in default
!%%!   if(leqi(ps_type,'siest')) then
!%%!      if(inode==ionode.AND.iprint_init>0) write(io_lun,fmt='(10x,"SIESTA pseudopotential will be used. ")')
!%%!      pseudo_type = SIESTA
!%%!   else if(leqi(ps_type,'plato').OR.leqi(ps_type,'abini')) then
!%%!      if(inode==ionode.AND.iprint_init>0) write(io_lun,fmt='(10x,"ABINIT pseudopotential will be used. ")')
!%%!      pseudo_type = ABINIT
!%%!   else
!%%!      if(inode==ionode.AND.iprint_init>0) write(io_lun,fmt='(10x,"OLD pseudopotential will be used. ")')
!%%!      pseudo_type = OLDPS
!%%!   endif
!%%!   if((.NOT.flag_angular_new).AND.(pseudo_type==SIESTA.OR.pseudo_type==ABINIT)) then
!%%!      write(io_lun,fmt='(10x,"Setting FlagNewAngular to T for Siesta/Abinit pseudopotentials")') 
!%%!      flag_angular_new = .true.
!%%!   end if
!%%!   ! Read, using old-style fdf_block, the information about the different species
!%%!   if(fdf_block('ChemicalSpeciesLabel',lun)) then
!%%!      do i=1,n_species
!%%!         read(lun,*) j,mass(i),species_label(i)
!%%!      end do
!%%!   end if
!%%!   ! Read charge, mass, pseudopotential and starting charge and blip models for the individual species
!%%!   maxnsf = 0
!%%!   min_blip_sp = 1.0e8_double
!%%!   do i=1,n_species
!%%!      charge(i) = zero
!%%!      nsf_species(i) = 0
!%%!      RadiusSupport(i) = r_h
!%%!      NonLocalFactor(i) = HNL_fac
!%%!      InvSRange(i) = r_t
!%%!      SupportGridSpacing(i) = zero
!%%!      if(pseudo_type==SIESTA.OR.pseudo_type==ABINIT) non_local_species(i) = .true.
!%%!      ! This is new-style fdf_block
!%%!      nullify(bp)
!%%!      if(fdf_block(species_label(i),bp)) then
!%%!         do while(fdf_bline(bp,line)) ! While there are lines in the block
!%%!            p=>digest(line)           ! Break the line up
!%%!            if(search(p,'charge',j)) then               ! Charge
!%%!               charge(i) = reals(p,1)
!%%!            else if(search(p,'NumberOfSupports',j)) then            ! NSF
!%%!               nsf_species(i) = integers(p,1)
!%%!            else if(search(p,'SupportFunctionRange',j)) then            ! Support radius
!%%!               RadiusSupport(i) = reals(p,1)
!%%!               if(InvSRange(i)<very_small) InvSRange(i) = RadiusSupport(i)
!%%!               !if(r_c<two*RadiusSupport(i)) r_c = two*RadiusSupport(i)
!%%!            else if(search(p,'InvSRange',j)) then            ! Support radius
!%%!               InvSRange(i) = reals(p,1)
!%%!            else if(search(p,'NonLocalFactor',j)) then            ! Support radius
!%%!               NonLocalFactor(i) = reals(p,1)
!%%!            else if(search(p,'SupportGridSpacing',j)) then            ! Support radius
!%%!               SupportGridSpacing(i) = reals(p,1)
!%%!               min_blip_sp = min(SupportGridSpacing(i),min_blip_sp)
!%%!               BlipWidth(i) = four*SupportGridSpacing(i)
!%%!            else if(pseudo_type==OLDPS.AND.search(p,'pseudopotential',j)) then ! Pseudopotential
!%%!               if(search(p,'non_local',j)) then         ! Non-local or local ?
!%%!                  non_local_species(i) = .true.
!%%!                  ps_file(i) = names(p,3)               ! We expect a line: pseudopotential (non_)local file 
!%%!               else
!%%!                  if(.NOT.search(p,'local',j)) then
!%%!                     call cq_abort('read_input: no local/non-local specification for pseudopotential !')
!%%!                  else
!%%!                     non_local_species(i) = .false.
!%%!                     ps_file(i) = names(p,3)            ! We expect a line: pseudopotential (non_)local file 
!%%!                  end if
!%%!               end if
!%%!            end if
!%%!         end do
!%%!         call destroy(p)  ! Remove storage
!%%!      end if
!%%!      call destroy(bp)    ! Remove storage
!%%!      if(nsf_species(i)==0) call cq_abort("Number of supports not specified for species ",i)
!%%!      if(flag_basis_set==blips.AND.SupportGridSpacing(i)<very_small) &
!%%!           call cq_abort("Error: for a blip basis set you must set SupportGridSpacing for all species")
!%%!      maxnsf = max(maxnsf,nsf_species(i))
!%%!      if(RadiusSupport(i)<very_small) &
!%%!           call cq_abort("Radius of support too small for species; increase SupportFunctionRange ",i)
!%%!   end do
!%%!   !blip_width = support_grid_spacing * fdf_double('blip_width_over_support_grid_spacing',four)
!%%!   L_tolerance = fdf_double('l_tolerance',1.0e-6_double)
!%%!   sc_tolerance  = fdf_double('sc_tolerance',0.001_double)
!%%!   maxpulayDMM = fdf_integer('DM.MaxPulay',5)
!%%!   LinTol_DMM = fdf_double('DM.LinTol',0.1_double)
!%%!   ! Find out what type of run we're doing
!%%!   runtype = fdf_string('MD.TypeOfRun','static')
!%%!   MDn_steps = fdf_integer('MD.NumSteps',100)
!%%!   MDfreq = fdf_integer('MD.OutputFreq',50)
!%%!   MDtimestep = fdf_double('MD.Timestep',0.5_double)
!%%!   MDcgtol = fdf_double('MD.MaxForceTol',0.0005_double)
!%%!   flag_vary_basis = fdf_boolean('VaryBasis',.false.)
!%%!   if(fdf_boolean('VaryBlips',.false.)) flag_vary_basis = .true.
!%%!   if(.NOT.flag_vary_basis) then
!%%!      flag_precondition_blips = .false.
!%%!   else 
!%%!      flag_precondition_blips = fdf_boolean('PreconditionBlips',.true.)
!%%!   end if
!%%!   flag_self_consistent = fdf_boolean('SelfConsistent',.true.)
!%%!   ! Tweak 2007/03/23 DRB Make Pulay mixing default
!%%!   flag_linear_mixing = fdf_boolean('LinearMixingSC',.true.)
!%%!   A = fdf_double('LinearMixingFactor',0.5_double)
!%%!   EndLinearMixing = fdf_double('LinearMixingEnd',sc_tolerance)
!%%!   q0 = fdf_double('KerkerFactor',0.1_double)
!%%!   n_exact = fdf_integer('LateStageReset',5)
!%%!   maxitersSC = fdf_integer('SC.MaxIters',50)
!%%!   maxearlySC = fdf_integer('SC.MaxEarly',3)
!%%!   maxpulaySC = fdf_integer('SC.MaxPulay',5)
!%%!   read_atomic_density_file = fdf_string('read_atomic_density_file','read_atomic_density.dat')
!%%!   ! Read atomic density initialisation flag
!%%!   atomic_density_method = fdf_string('atomic_density_flag','pao')
!%%!   InvSTolerance = fdf_double('InvSTolerance',1e-2_double)
!%%!   flag_read_blocks = fdf_boolean('ReadBlocks',.false.)
!%%!   ewald_accuracy = fdf_double('EwaldAccuracy',1.0e-10_double) ! Default value of 10^-10 Ha per atom
!%%!   flag_old_ewald = fdf_boolean('FlagOldEwald',.false.)
!%%!   UseGemm = fdf_boolean('MM.UseGemm',.false.)
!%%!   del_k = fdf_double('pao_kspace_ol_gridspace',0.1_double)
!%%!   kcut = fdf_double('pao_kspace_ol_cutoff', 1000.0_double)
!%%!   flag_paos_atoms_in_cell = fdf_boolean('BasisSet.PAOs_StoreAllAtomsInCell',.true.)
!%%!   read_option = fdf_boolean('read_pao_coeffs',.false.)
!%%!   symmetry_breaking = fdf_boolean('BasisSet.SymmetryBreaking',.false.)
!%%!   support_pao_file = fdf_string('support_pao_file','supp_pao.dat')
!%%!   pao_info_file = fdf_string('pao_info_file','pao.dat')
!%%!   pao_norm_flag = fdf_integer('pao_norm_flag',0)
!%%!   TestBasisGrads = fdf_boolean('TestPAOGrads',.false.)
!%%!   TestTot = fdf_boolean('TestPAOGradTot',.false.)
!%%!   TestBoth = fdf_boolean('TestPAOGradBoth',.false.)
!%%!   TestS = fdf_boolean('TestPAOGrad_S',.false.)
!%%!   TestH = fdf_boolean('TestPAOGrad_H',.false.)
!%%!   support_spec_file = fdf_string('support_spec_file','support.dat')
!%%!   flag_read_support_spec = fdf_boolean('Basis.ReadSupportSpec',.false.)
!%%!   flag_test_forces = fdf_boolean('TestForces',.false.)
!%%!   flag_test_all_forces = fdf_boolean('TestAllForces',.true.)
!%%!   if(.NOT.flag_test_all_forces) then ! Test one force component
!%%!      flag_which_force = fdf_integer('TestSpecificForce',1)
!%%!   end if
!%%!   TF_direction = fdf_integer('TestForceDirection',1)
!%%!   TF_atom_moved = fdf_integer('TestForceAtom',1)
!%%!   TF_delta = fdf_double('TestForceDelta',0.00001_double)
!%%!   flag_functional_type = fdf_integer('FunctionalType', 1)   ! LDA PZ81
!%%!   select case(flag_functional_type)
!%%!   case (functional_lda_pz81)
!%%!      functional_description = 'LDA PZ81'
!%%!   case (functional_lda_gth96)
!%%!      functional_description = 'LDA GTH96'
!%%!   case (functional_lda_pw92)
!%%!      functional_description = 'LDA PW92'
!%%!   case (functional_gga_pbe96)
!%%!      functional_description = 'GGA PBE96'
!%%!   case default
!%%!      functional_description = 'LDA PZ81'
!%%!   end select
!%%!   tmp = fdf_string('EnergyUnits','Ha')
!%%!   if(leqi(tmp(1:2),'Ha')) then
!%%!      energy_units = har
!%%!      en_conv = one
!%%!   else if(leqi(tmp(1:2),'Ry')) then
!%%!      energy_units = ryd
!%%!      en_conv = HaToRy
!%%!   else if(leqi(tmp(1:2),'eV')) then
!%%!      energy_units = ev
!%%!      en_conv = HaToeV
!%%!   endif
!%%!   tmp = fdf_string('DistanceUnits','bohr')
!%%!   if(leqi(tmp(1:2),'a0').OR.leqi(tmp(1:2),'bo')) then
!%%!      dist_units = bohr
!%%!      dist_conv = one
!%%!   else if(leqi(tmp(1:1),'A')) then
!%%!      dist_units = ang
!%%!      dist_conv = BohrToAng
!%%!   endif
!%%!   for_conv = en_conv/dist_conv
!%%!   pdb_format = fdf_boolean('General.PdbIn',.false.)
!%%!   pdb_altloc = fdf_string('General.PdbAltLoc',' ')
!%%!   pdb_output = fdf_boolean('General.PdbOut',.false.)
!%%!   part_mode = fdf_string('General.Partitions','Python')
!%%!   if (leqi (part_mode(1:6),'Python')) then
!%%!      part_method = PYTHON
!%%!   else if (leqi (part_mode(1:7), 'Hilbert')) then
!%%!      part_method = HILBERT
!%%!   else
!%%!      part_method = HILBERT
!%%!      write(io_lun,*) 'WARNING: Unrecognised partitioning method'
!%%!   end if
!%%!   tmp2 = fdf_string('General.LoadBalance','atoms')
!%%!   if (leqi (tmp2(1:10),'partitions')) then
!%%!      load_balance = 1
!%%!   else if (leqi (tmp2(1:5),'atoms')) then
!%%!      load_balance = 0
!%%!   else
!%%!      load_balance = 0
!%%!   end if
!%%!   many_processors = fdf_boolean('General.ManyProcessors',.false.)
!%%!   global_maxatomspart = fdf_integer('General.MaxAtomsPartition', 34)
!%%!   append_coords = fdf_boolean('MD.AppendCoords',.true.)
    end if ! new_format
    !!!   TEMPORARY ? !!!!!!!   by T. MIYAZAKI
    call my_barrier
    return
  end subroutine read_input
!!***

! ------------------------------------------------------------------------------
! Subroutine allocate_species_vars
! ------------------------------------------------------------------------------

!!****f* initial_read/allocate_species_vars *
!!
!!  NAME 
!!   allocate_species_vars
!!  USAGE
!!   
!!  PURPOSE
!!   Allocates variables which depend on number of species
!!  INPUTS
!! 
!! 
!!  USES
!!   
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   Unknown
!!  MODIFICATION HISTORY
!!   2011/09/16 11:10 dave
!!    Added header and changed atomicrad to atomicnum
!!   2011/11/16 15:51 dave
!!    Changes for new blip data
!!  SOURCE
!!
  subroutine allocate_species_vars

    use dimens, ONLY: RadiusSupport, NonLocalFactor, InvSRange, atomicnum
    use blip, ONLY: blip_info
    use memory_module, ONLY: reg_alloc_mem, type_dbl
    use species_module, ONLY: nsf_species, nlpf_species, npao_species, charge, mass, non_local_species, &
         ps_file, ch_file, phi_file, species_label, n_species, type_species
    use global_module, ONLY: area_general
    use GenComms, ONLY: cq_abort

    implicit none

    ! Local variables
    integer :: stat
    
    call start_timer(tmr_std_allocation)
    allocate(RadiusSupport(n_species),atomicnum(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating RadiusSupport, atomicnum in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(blip_info(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating blip_info in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(NonLocalFactor(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating NonLocalFactor in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(nsf_species(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating nsf_species in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(nlpf_species(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating nlpf_species in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(npao_species(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating npao_species in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(charge(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating charge in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(mass(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating mass in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(non_local_species(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating non_local_species in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(ps_file(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating ps_file in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(ch_file(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating ch_file in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(phi_file(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating phi_file in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(species_label(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating species_label in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(type_species(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating type_species in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(InvSRange(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating InvSRange in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine allocate_species_vars
!!***

! ------------------------------------------------------------------------------
! Subroutine write_info
! ------------------------------------------------------------------------------

!!****f* initial_read/write_info *
!!
!!  NAME 
!!   write_info
!!  USAGE
!! 
!!  PURPOSE
!!   Writes out information about the run
!!  INPUTS
!!   Many - all the things read in
!!  USES
!!   datatypes, dimes, species_module, pseudopotential_data, atoms
!!  AUTHOR
!!   E.H.Hernandez
!!  CREATION DATE
!!   9/11/93 (!)
!!  MODIFICATION HISTORY
!!   10/8/01 DRB 
!!    Changed to F90 format and tidied
!!   10/05/2002 dave
!!    Changed titles to be a single line
!!   10/05/2002 dave
!!    Imported into initial_read
!!   29/05/2002 dave
!!    Added write out for solution method
!!   14:47, 2006/07/17
!!    Added description of the functional
!!   2010/06/18 17:00 lt
!!    Added information on the type of smearing used if using Diagonalsation method
!!  SOURCE
!!
  subroutine write_info(number_of_bands, titles, mu, vary_mu, find_chdens, read_phi,HNL_fac, NODES)

    use datatypes
    use units
    use dimens, ONLY: r_super_x, r_super_y, r_super_z, &
         n_grid_x, n_grid_y, n_grid_z, r_h, r_c
    use block_module, ONLY: in_block_x, in_block_y, in_block_z
    use species_module, ONLY: n_species, species_label, mass, charge, &
         ps_file, ch_file, phi_file, species, nsf_species
    use pseudopotential_data, ONLY: core_radius, non_local_species
    use DiagModule, ONLY: diagon, flag_smear_type, iMethfessel_Paxton
    use blip, ONLY: blip_info
    use global_module, ONLY: flag_basis_set, PAOs,blips, functional_description, &
         flag_precondition_blips, io_lun
    use minimise, ONLY: energy_tolerance, L_tolerance, sc_tolerance, &
         n_support_iterations, n_L_iterations
    use datestamp, ONLY: datestr, commentver

    implicit none

    ! Passed variables
    logical :: vary_mu, find_chdens
    logical :: read_phi

    character(len=80) :: titles

    integer :: NODES 
    integer :: n_title

    real(double) :: mu, number_of_bands, &
         alpha, beta, expected_reduction, HNL_fac

    ! Local variables
    integer :: n
    character(len=10) :: today, the_time

    call date_and_time(today, the_time)
    write(io_lun,3) today(1:4), today(5:6), today(7:8), the_time(1:2), the_time(3:4)
    write(io_lun,'(/10x,"Code compiled on: ",a,/10x,"Version comment: ",/10x,a)') datestr,commentver
    
    write(io_lun,1)
    write(io_lun,2) titles

    if(diagon) then
       write(io_lun,30) 'diagonalisation '
       select case (flag_smear_type)
       case (0)
          write(io_lun,'(/,10x,"Using Fermi-Dirac smearing")')
       case (1)
          write(io_lun,'(/,10x,"Using order ",i2," Methfessel-Paxton smearing")') iMethfessel_Paxton
       end select
    else
       write(io_lun,30) 'order N with LNV'   
    end if
    write(io_lun,4) dist_conv*r_super_x, dist_conv*r_super_y, dist_conv*r_super_z,d_units(dist_units)

    write(io_lun,9) n_grid_x, n_grid_y, n_grid_z

    write(io_lun,17) in_block_x, in_block_y, in_block_z

    write(io_lun,15) dist_conv*(r_super_x/n_grid_x),d_units(dist_units), &
         dist_conv*(r_super_y/n_grid_y), d_units(dist_units),&
         dist_conv*(r_super_z/n_grid_z),d_units(dist_units)

    write(io_lun,18) n_species, d_units(dist_units)

    do n=1, n_species
       write(io_lun,19) n, species_label(n), mass(n), charge(n), dist_conv*core_radius(n), nsf_species(n)
       if(flag_basis_set==blips) then 
          if(flag_precondition_blips) then
             write(io_lun,'(/13x,"Blip basis with preconditioning")') 
          else
             write(io_lun,'(/13x,"Blip basis - no preconditioning")') 
          end if
          write(io_lun,14) dist_conv*blip_info(n)%SupportGridSpacing, d_units(dist_units),&
               dist_conv*blip_info(n)%BlipWidth,d_units(dist_units)
       else
          write(io_lun,'(13x,"PAO basis")') 
       end if
    end do

    write(io_lun,20)

    write(io_lun,29) energy_tolerance, L_tolerance, sc_tolerance

    !write(io_lun,26)
    !do n=1, n_species
    !   write(io_lun,21) species_label(n), ps_file(n)
    !end do
    write(io_lun,27)
    do n=1, n_species
       if ( non_local_species(n) ) then
          write(io_lun,22) species_label(n), 'Non Local'
       else 
          write(io_lun,22) species_label(n), 'Local    '
       end if
    end do

    !if(.NOT.find_chdens) then
    !   write(io_lun,261)
    !   do n=1, n_species
    !      write(io_lun,211) species_label(n), ch_file(n)
    !   end do
    !endif
    if(read_phi) then
       write(io_lun,262)
       do n=1, n_species
          write(io_lun,212) species_label(n), phi_file(n)
       end do
    endif

    write(io_lun,6) number_of_bands

    if (.not.vary_mu) then
       write(io_lun,*) '          mu is constant'
       write(io_lun,16) mu
    endif

    write(io_lun,7) NODES

    write(io_lun,8) functional_description

    write(io_lun,11) n_support_iterations, n_L_iterations

    write(io_lun,13) dist_conv*r_h, d_units(dist_units),dist_conv*r_c,d_units(dist_units)
    do n=1, n_species
       write(io_lun,131) n,dist_conv*r_h+core_radius(n)*HNL_fac,d_units(dist_units)
    enddo

1   format(/10x,'Job title: ')
2   format(10x,a80)

3   format(/10x,'This job was run on ',a4,'/',a2,'/',a2,' at ',a2,':',a2,/)

4   format(/10x,'The simulation box has the following dimensions',/, &
         10x,'a = ',f11.5,' b = ',f11.5,' c = ',f11.5,' 'a2)

5   format(/10x,'The simulation box contains ',i7,' atoms.')

6   format(/10x,'The number of bands is ',f8.2)

7   format(/10x,'The calculation will be performed on ',i5,' processors')

8   format(/10x, 'The functional used will be ', a15)

9   format(/10x,'The number of cell grid points in each direction is :',/, &
         20x,i5,' cell grid points along x',/, &
         20x,i5,' cell grid points along y',/, &
         20x,i5,' cell grid points along z')

11  format(/10x,'The maximum number of support iterations will be ',i5, &
         ' consisting of:',/, &
         20x,i5,' L matrix variations.')

13  format(/10x,'Support region radius = ',f7.4,' ',a2,' ',/, &
         10x,'Density Matrix range  = ',f7.4,' ',a2)
131 format(/10x,'Species ',i2,' Non-local Hamiltonian radius = ', &
         f7.4,' ',a2)

14  format(/10x,'Support-grid spacing =   ',f7.4,' ',a2,' ',/, &
         10x,'Width of (3D) b-spline = ',f7.4,' ',a2)

15  format(/10x,'integration grid spacing along x ',f7.5,' ',a2,/, &
         10x,'integration grid spacing along y ',f7.5,' ',a2,/, &
         10x,'integration grid spacing along z ',f7.5,' ',a2)

16  format(/10x,'The Chemical Potential mu is :',f7.4)

17  format(/10x,'The number of cell grid points in each block is :',/, &
         20x,i5,' cell grid points along x',/, &
         20x,i5,' cell grid points along y',/, &
         20x,i5,' cell grid points along z')

18  format(/10x,'The number of atomic species in the system is :', &
         i5,/,/,6x, &
         '------------------------------------------------------------------',/,6x, &
         '   #   Label     Mass (a.u.)   Charge (e)  NLPF Rad (',a2,')   NSF  ',/,6x, &
         '------------------------------------------------------------------')

19  format(6x,i4,3x,a30,3f11.5,9x,i3)

20  format(6x, &
         '------------------------------------------------------------------',/)

21  format(20x,a10,10x,a40)
211 format(20x,a10,10x,a40)
212 format(20x,a10,10x,a40)

22  format(20x,a30,10x,a10) 

23  format(/,10x, &
         '-------------------------------------------------------',/, &
         10x,'   Atom No.    Species      rx        ry       rz       ' &
         ,/,10x,'-------------------------------------------------------')

24  format(12x,i5,10x,a5,3f10.5)

25  format(/, &
         10x,'-------------------------------------------------------')

26  format(/,10x,'The pseudopotential data is read from:',/, &
         15x,'  Species ',12x,'  File  ')
261 format(/,10x,'The initial charge density data is read from:',/, &
         15x,'  Species ',12x,'  File  ')
262 format(/,10x,'The initial support function data is read from:',/, &
         15x,'  Species ',12x,'  File  ')

27  format(/,10x,'The pseudopotentials are of type:',/, &
           15x,'  Species ',12x,'  Type  ')

28  format(/,10x,'Exponent for the s-type support functions: ', &
           f10.7,/,10x,'Exponent for the p-type support functions: ', &
           f10.7)

29  format(/,10x,'Energy tolerance required:             ',f10.8, &
           /,10x,'L-matrix convergence tolerance:        ',f10.8, &
           /,10x,'Self consistent convergence tolerance: ',f10.8)
30  format(/,10x,'Solving for the K matrix using ',a16)

    return
  end subroutine write_info
!!***



! -----------------------------------------------------------------------------
! Subroutine readDiagInfo
! -----------------------------------------------------------------------------

!!****f* DiagModule/readDiagInfo *
!!
!!  NAME 
!!   readDiagInfo
!!  USAGE
!!   readDiagInfo(Number of Kpts, proc grid rows, proc grid cols, row block size, col block size)
!!  PURPOSE
!!   Reads in data to do with diagonalisation, allocates memory for data and broadcasts to all
!!   processors
!!  INPUTS
!!   integer, intent(OUT) :: proc_rows, proc_cols - rows and columns in processor grid
!!   integer, intent(OUT) :: block_size_r, block_size_c - block sizes for rows and columns of matrix
!!  USES
!!   GenComms, fdf, global_module, numbers
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   16/04/2002
!!  MODIFICATION HISTORY
!!   22/04/2002 dave
!!    Added weight of 1.0 to default k-point set (i.e. gamma point only)
!!   29/05/2002 dave
!!    Removed fdf_init (now all Conquest input is done via fdf)
!!   2006/09/22 08:20 dave
!!    Relocated into initial_read
!!   2010/07/22 15.53 Lianheng
!!    moved control flags associated with Diagonalisation method from
!!    the main read_input() subroutine to here.  added
!!    Methfessel-Paxton smearing and related control flags. Changed
!!    SC.MaxEfIter to Diag.MaxEfIter for greater consistency
!!   2011/02/13 L.Tong
!!    added Diag.KProcGroups for k-point parallelisation
!!    IMPORTANT: for now Diag.KProcGroups must be an integer factor
!!    of the total number of processes. In other words each group
!!    must contain the same number of processes.
!!   2011/09/02 L.Tong
!!    Corrected Conquest default values for ScaLAPACK processor grid
!!    definitions when k-point parallelisation is on. 
!!    Added local variable proc_per_group
!!   2011/09/05 L.Tong
!!    Added print out input information for k-point parallelisation
!!    related parameters
!!  SOURCE
!!
  subroutine readDiagInfo

    use datatypes
    use global_module, ONLY: iprint_init, rcellx, rcelly, rcellz, area_general, ni_in_cell, numprocs, &
         species_glob, io_lun
    use numbers, ONLY: zero, one, two, pi, very_small
    use GenComms, ONLY: cq_abort, gcopy, myid
    use input_module
    use ScalapackFormat, ONLY: proc_rows, proc_cols, block_size_r, block_size_c, proc_groups
    use DiagModule, ONLY: nkp, kk, wtk, kT, maxefermi, flag_smear_type, iMethfessel_Paxton, max_brkt_iterations, &
         gaussian_height, finess, NElec_less
    use energy, ONLY: SmearingType, MPOrder
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use species_module, ONLY: nsf_species

    implicit none

    ! Local variables
    integer :: stat, iunit, i, j, k, matrix_size
    real(double) :: a, sum
    integer :: proc_per_group
    ! k-point mesh type
    logical :: mp_mesh, done
    real(double), dimension(1:3) :: mp_shift
    real(double), allocatable, dimension(:,:) :: kk_tmp
    real(double), allocatable, dimension(:) :: wtk_tmp
    integer, dimension(1:3) :: mp
    integer :: nkp_tmp
    integer :: counter
    
    if(myid==0) then
       ! Read Control Flags associated to diagonalisation method
       maxefermi = fdf_integer('Diag.MaxEfIter',50)
       kT = fdf_double('Diag.kT',0.001_double)
       ! Method to approximate step function for occupation number
       flag_smear_type = fdf_integer('Diag.SmearingType',0)
       SmearingType = flag_smear_type
       iMethfessel_Paxton = fdf_integer('Diag.MPOrder',0)
       MPOrder = iMethfessel_Paxton
       max_brkt_iterations = fdf_integer('Diag.MaxBrktSteps',50)
       gaussian_height = fdf_double('Diag.GaussianHeight',0.1_double)
       finess = fdf_double('Diag.EfStepFiness',1.0_double)
       NElec_less = fdf_double('Diag.NElecLess',10.0_double)
       ! Read k-point parallelisation process-group incormation, default is 1
       proc_groups = fdf_integer ('Diag.KProcGroups', 1)
       if (iprint_init > 0) write (io_lun, 11) proc_groups 
       ! Read/choose ScaLAPACK processor grid dimensions
       if(fdf_defined('Diag.ProcRows')) then
          proc_rows = fdf_integer('Diag.ProcRows',0)
          proc_cols = fdf_integer('Diag.ProcCols',0)
          if(proc_rows*proc_cols==0) call cq_abort('Diag: proc grid product is zero: ',proc_rows,proc_cols)
          ! 2011/09/02 L.Tong: now checks correctly if we have more than one k-point groups
          proc_per_group = numprocs / proc_groups
          if(proc_rows*proc_cols > proc_per_group) &
               call cq_abort('Diag: proc grid product is too large: ', proc_rows*proc_cols, proc_per_group) 
       else
          ! default values
          if(numprocs<9) then ! Recommended by ScaLapack User Guide
             proc_rows = 1
             proc_cols = numprocs / proc_groups
          else
             done = .false.
             ! 2011/09/02 L.Tong, now works properly for more than one k-point groups
             proc_per_group = numprocs / proc_groups
             a = sqrt(real(proc_per_group, double))
             proc_rows = aint(a)+1
             do while(.NOT.done.AND.proc_rows>1) 
                proc_rows = proc_rows - 1
                proc_cols = aint(real(proc_per_group,double) / real(proc_rows,double))
                if(proc_per_group - proc_rows*proc_cols == 0) done = .true.
             end do
          end if
       end if
       ! Read/choose ScaLAPACK block sizes
       matrix_size = 0 
       do i=1,ni_in_cell
          matrix_size = matrix_size + nsf_species(species_glob(i))
       end do
       if(fdf_defined('Diag.BlockSizeR')) then
          block_size_r = fdf_integer('Diag.BlockSizeR',1)
          block_size_c = fdf_integer('Diag.BlockSizeC',1)
          a = real(matrix_size)/real(block_size_r)
          if(a - real(floor(a))>1e-8_double) &
               call cq_abort('block_size_r not a factor of matrix size ! ',matrix_size, block_size_r)
          a = real(matrix_size)/real(block_size_c)
          if(a - real(floor(a))>1e-8_double) &
               call cq_abort('block_size_c not a factor of matrix size ! ',matrix_size, block_size_c)
       else
          done = .false.
          block_size_r = matrix_size/max(proc_rows,proc_cols)+1
          do while(.NOT.done) 
             block_size_r = block_size_r - 1
             a = real(matrix_size)/real(block_size_r)
             if(a - real(floor(a))<1e-8_double) done = .true.
          end do
          if(block_size_r==1) then
             if(numprocs==1) then
                block_size_r = matrix_size
             else
                call cq_abort("Can't find a good block size: please set manually")
             end if
          end if
          block_size_c = block_size_r
       end if
       if(iprint_init>1) then
          write(io_lun,2) block_size_r, block_size_c
          write(io_lun,3) proc_rows, proc_cols
       end if
       ! Read k-point mesh type
       mp_mesh = fdf_boolean('Diag.MPMesh',.false.)
       if(.NOT.mp_mesh) then
          ! Read k-point number and allocate
          nkp = fdf_integer('Diag.NumKpts',1)
          if(iprint_init>1) write(io_lun,fmt='(8x,"Number of Kpoints: ",i4)') nkp
          if(nkp<1) call cq_abort("Need to specify how many kpoints !",nkp)
          allocate(kk(3,nkp),wtk(nkp),STAT=stat)
          if(stat/=0) call cq_abort('FindEvals: couldnt allocate kpoints',nkp)
          call reg_alloc_mem(area_general,4*nkp,type_dbl)
          sum = zero
          ! Read k-points
          if(fdf_block('Diag.Kpoints'))then
             if(1+block_end-block_start<nkp) call cq_abort("Kpoint error: ",1+block_end-block_start,nkp)
             do i=1,nkp
                read(unit=input_array(block_start+i-1),fmt=*) kk(1,i),kk(2,i),kk(3,i),wtk(i)
                ! Assume fractional kpoints and orthorhombic cell for now
                kk(1,i) = two*pi*kk(1,i)/rcellx
                kk(2,i) = two*pi*kk(2,i)/rcelly
                kk(3,i) = two*pi*kk(3,i)/rcellz
                sum = sum + wtk(i)
             end do
             call fdf_endblock
             wtk = wtk/sum
          else ! Force gamma point dependence
             write(io_lun,4)
             nkp = 1
             kk(1,1) = zero
             kk(2,1) = zero
             kk(3,1) = zero
             wtk(1) = one
          end if
       else
          ! Read Monkhorst-Pack mesh coefficients
          ! Default is Gamma point only 
          if(iprint_init>0) then
             write(io_lun,fmt='(/8x,"Reading Monkhorst-Pack Kpoint mesh"//)')
          end if
          ! giving up on "MPMesh i j k" for the moment
          ! mp_string = fdf_string('MPMesh','1 1 1')
          ! This must be posible!
          ! But maybe not using fdf
          !
          ! Check whether the mesh is specified in input
          ! if (fdf_boolean('mp_mesh_x',.false.)==.false.) write (*,8)
          ! The above causes an error message in case mp_mesh_x exists, because
          ! it's then handled as boolean below
          mp(1) = fdf_integer('Diag.MPMeshX',1)
          mp(2) = fdf_integer('Diag.MPMeshY',1)
          mp(3) = fdf_integer('Diag.MPMeshZ',1) 
          if(iprint_init>0) write(io_lun,fmt='(8x,a, 3i3)') ' Monkhorst-Pack mesh: ', (mp(i), i=1,3)
          if (mp(1) <= 0) call cq_abort('K-points: number of k-points must be > 0!')
          if (mp(2) <= 0) call cq_abort('K-points: number of k-points must be > 0!')
          if (mp(3) <= 0) call cq_abort('K-points: number of k-points must be > 0!')
          nkp_tmp = mp(1)*mp(2)*mp(3)
          if(iprint_init>0) write(io_lun,fmt='(8x,a, i4)') ' Number of k-points: ',nkp_tmp
          ! Read k-point shift, default (0.0 0.0 0.0)
          mp_shift(1) = fdf_double('Diag.MPShiftX',zero)
          mp_shift(2) = fdf_double('Diag.MPShiftY',zero)
          mp_shift(3) = fdf_double('Diag.MPShiftZ',zero)
          if (mp_shift(1) >= one) then
             write(io_lun,9)
             mp_shift(1) = mp_shift(1) - one
          end if
          if (mp_shift(2) >= one) then
             write(io_lun,9)
             mp_shift(2) = mp_shift(2) - one
          end if
          if (mp_shift(3) >= one) then
             write(io_lun,9)
             mp_shift(3) = mp_shift(3) - one
          end if
          if(iprint_init>0) write(io_lun,fmt='(8x,a, 3f11.6)') ' Monkhorst-Pack mesh shift:  ', (mp_shift(i), i=1,3)
          ! Allocate
          allocate(kk_tmp(3,nkp_tmp),wtk_tmp(nkp_tmp),STAT=stat)
          if(stat/=0) call cq_abort('FindEvals: couldnt allocate kpoints',nkp_tmp)
          call reg_alloc_mem(area_general,4*nkp_tmp,type_dbl)
          ! All k-points have weight 1 for now
          wtk_tmp(1:nkp_tmp) = one
          if(iprint_init>0) write(io_lun,*)
          ! Generate fractional k-point coordinates plus shift
          ! Assume orthorhombic cell for now
          do i = 0, mp(1) -1 ! x axis 
             do j = 0, mp(2) -1 ! y axis
                do k = 1, mp(3) ! z axis
                   counter = k + j * mp(3) + i * mp(2) * mp(3)
                   kk_tmp(1,counter) = ( (two * (i+1) - mp(1) - 1) / (2 * mp(1)) ) + mp_shift(1)
                   kk_tmp(2,counter) = ( (two * (j+1) - mp(2) - 1) / (2 * mp(2)) ) + mp_shift(2)
                   kk_tmp(3,counter) = ( (two * k - mp(3) - 1) / (2 * mp(3)) ) + mp_shift(3)
                end do
             end do
          end do
          ! Write out fractional k-points
          if(iprint_init>0) then
             write(io_lun,7) nkp_tmp
             do i=1,nkp_tmp
                write(io_lun,fmt='(8x,i5,3f15.6,f12.3)') i,kk_tmp(1,i),kk_tmp(2,i),kk_tmp(3,i),wtk_tmp(i)
             end do
          end if
          ! Check if k = -k and weed out equivalent k-points. Adjust wtk accordingly.
          do i = nkp_tmp, 1, -1 ! to get most coords > 0
             if (wtk_tmp(i) > zero) then
                do j = nkp_tmp, 1, -1
                   if (i /= j) then
                      if ((kk_tmp(1,i) == -kk_tmp(1,j)) .and. (kk_tmp(2,i) == -kk_tmp(2,j)) &
                           .and. (kk_tmp(3,i) == -kk_tmp(3,j))) then
                         wtk_tmp(j) = zero
                         wtk_tmp(i) = wtk_tmp(i) + one
                      end if
                   end if
                end do
             end if
          end do

          ! How many k-points are now left? 
          nkp = 0
          ! DRB fix 
          sum = zero
          do i = 1, nkp_tmp
             if (wtk_tmp(i) > very_small ) then
                nkp = nkp + 1 
                sum = sum + wtk_tmp(i)
             end if
          end do
          ! Fill in the real k-point array
          allocate(kk(3,nkp),wtk(nkp),STAT=stat)
          if(stat/=0) call cq_abort('FindEvals: couldnt allocate kpoints',nkp)
          call reg_alloc_mem(area_general,4*nkp,type_dbl)
          counter = 0
          do i = 1, nkp_tmp
             if (wtk_tmp(i) > very_small ) then
                counter = counter + 1
                kk(1,counter) = kk_tmp(1,i)
                kk(2,counter) = kk_tmp(2,i)
                kk(3,counter) = kk_tmp(3,i)
                !DRB fix to sum weights to one
                wtk(counter) = wtk_tmp(i)/sum
             end if
          end do
          if(counter>nkp) call cq_abort("Error in kpt symmetry ! ",counter,nkp)
          call reg_dealloc_mem(area_general,4*nkp_tmp,type_dbl)
          deallocate(kk_tmp,wtk_tmp,STAT=stat)
          if(stat/=0) call cq_abort('FindEvals: couldnt deallocate kpoints',nkp_tmp)
          if(iprint_init>0) then
             write(io_lun,*)
             write(io_lun,10) nkp
             do i=1,nkp
                write(io_lun,fmt='(8x,i5,3f15.6,f12.3)') i,kk(1,i),kk(2,i),kk(3,i),wtk(i)
             end do
          end if

          do i = 1, nkp
             kk(1,i) = two * pi * kk(1,i) / rcellx
             kk(2,i) = two * pi * kk(2,i) / rcelly
             kk(3,i) = two * pi * kk(3,i) / rcellz
          end do
       end if ! MP mesh branch

       ! Write out k-points
       if(iprint_init>0) then
          write(io_lun,*)
          write(io_lun,51) nkp
          write(io_lun,52)
          do i=1,nkp
             write(io_lun,fmt='(8x,i5,3f15.6,f12.3)') i,kk(1,i),kk(2,i),kk(3,i),wtk(i)
          end do
          write(io_lun,fmt='(/8x,"Finished reading Kpoints"/)')
       end if
       
       ! Write out smearing temperature
       if(iprint_init>0) write(io_lun,'(10x,"Temperature used for smearing: ",f10.6)') kT
    end if ! myid==0

    ! Distribute data to all processors
    call gcopy(block_size_r)
    call gcopy(block_size_c)
    call gcopy (proc_groups)
    call gcopy(proc_rows)
    call gcopy(proc_cols)
    call gcopy(nkp)
    if(myid/=0) then
       allocate(kk(3,nkp),wtk(nkp),STAT=stat)
       if(stat/=0) call cq_abort('FindEvals: couldnt allocate kpoints',nkp)
       call reg_alloc_mem(area_general,4*nkp,type_dbl)
    end if
    call gcopy(kk,3,nkp)
    call gcopy(wtk,nkp)
    call gcopy(kT)
    call gcopy (maxefermi)
    call gcopy (flag_smear_type)
    call gcopy (SmearingType)
    call gcopy (iMethfessel_Paxton)
    call gcopy (MPOrder)
    call gcopy (max_brkt_iterations)
    call gcopy (gaussian_height)
    call gcopy (finess)
    call gcopy (NElec_less)
    !if(iprint_init>=5.AND.myid/=0) then
    !   write(io_lun,*) 'Proc: ',myid
    !   write(io_lun,5)
    !   do i=1,nkp
    !      write(io_lun,6) nkp,kk(1,i),kk(2,i),kk(3,i),wtk(i)
    !   end do
    !end if
    return
2   format(8x,'Block size (row, col): ',2i5)
3   format(8x,'Proc grid (row, col): ',2i5)
4   format(/8x,'***WARNING***',/,2x,'No Kpoints block found - defaulting to Gamma point')
5   format(8x,i4,' Kpoints in Cartesian (inverse Angstrom) form: ')
51  format(8x,i4,' symmetry inequivalent Kpoints in Cartesian ')
52  format(12x,' (inverse Angstrom) form: ')
6   format(8x,'Kpt: ',i4,' : ',3f12.8,' Weight: ',f10.6)
7   format(8x,' All ',i4,' Kpoints in fractional coordinates: ')
8   format(/8x,'***WARNING***',/,8x,'No Kpoint mesh found - defaulting to Gamma point')
9   format(/8x,'***WARNING***',/,8x,'Specified Kpoint mesh shift in fractional coords >= 1.0.')
10  format(8x,i4,' symmetry inequivalent Kpoints in fractional coordinates: ')
11  format(8x, 'Number of processor groups for k-point parallelisation: ', i5)
  end subroutine readDiagInfo
!!***

!%%!   subroutine readDiagInfo
!%%! 
!%%!     use datatypes
!%%!     use global_module, ONLY: iprint_init, rcellx, rcelly, rcellz, area_general
!%%!     use numbers, ONLY: two, pi
!%%!     use GenComms, ONLY: cq_abort, gcopy, myid
!%%!     use fdf, ONLY: fdf_integer, fdf_block
!%%!     use ScalapackFormat, ONLY: proc_rows, proc_cols, block_size_r, block_size_c
!%%!     use DiagModule, ONLY: nkp, kk, wtk
!%%!     use memory_module, ONLY: reg_alloc_mem, type_dbl
!%%! 
!%%! 
!%%!     implicit none
!%%! 
!%%!     ! Local variables
!%%!     integer :: stat, iunit, i
!%%!     
!%%!     if(myid==0) then
!%%!        ! Read ScaLAPACK block sizes
!%%!        block_size_r = fdf_integer('SCBlockSizeR',1)
!%%!        block_size_c = fdf_integer('SCBlockSizeC',1)
!%%!        ! Read ScaLAPACK processor grid dimensions
!%%!        proc_rows = fdf_integer('ProcRows',0)
!%%!        proc_cols = fdf_integer('ProcCols',0)
!%%!        if(iprint_init>1) then
!%%!           write(io_lun,2) block_size_r, block_size_c
!%%!           write(io_lun,3) proc_rows, proc_cols
!%%!        end if
!%%!        if(proc_rows*proc_cols==0) call cq_abort('FindEval: error in proc grid',proc_rows,proc_cols)
!%%!        ! Read k-point number and allocate
!%%!        nkp = fdf_integer('NumKpts',1)
!%%!        if(nkp<1) call cq_abort("Need to specify how many kpoints !",nkp)
!%%!        allocate(kk(3,nkp),wtk(nkp),STAT=stat)
!%%!        if(stat/=0) call cq_abort('FindEvals: couldnt allocate kpoints',nkp)
!%%!        call reg_alloc_mem(area_general,4*nkp,type_dbl)
!%%!        ! Read k-points
!%%!        if(fdf_block('Kpoints',iunit))then
!%%!           do i=1,nkp
!%%!              read(iunit,*) kk(1,i),kk(2,i),kk(3,i),wtk(i)
!%%!              ! Assume fractional kpoints and orthorhombic cell for now
!%%!              kk(1,i) = two*pi*kk(1,i)/rcellx
!%%!              kk(2,i) = two*pi*kk(2,i)/rcelly
!%%!              kk(3,i) = two*pi*kk(3,i)/rcellz
!%%!           end do
!%%!        else ! Force gamma point dependence
!%%!           write(io_lun,4)
!%%!           nkp = 1
!%%!           kk(1,1) = 0.0_double
!%%!           kk(2,1) = 0.0_double
!%%!           kk(3,1) = 0.0_double
!%%!           wtk(1) = 1.0_double
!%%!        end if
!%%!        ! Write out k-points
!%%!        if(iprint_init>0) then
!%%!           write(io_lun,'(2x,i4," Kpoints in Cartesian (inverse Angstrom) form: ")') nkp
!%%!           do i=1,nkp
!%%!              write(io_lun,6) nkp,kk(1,i),kk(2,i),kk(3,i),wtk(i)
!%%!           end do
!%%!        end if
!%%!     end if
!%%!     ! Distribute data to all processors
!%%!     call gcopy(block_size_r)
!%%!     call gcopy(block_size_c)
!%%!     call gcopy(proc_rows)
!%%!     call gcopy(proc_cols)
!%%!     call gcopy(nkp)
!%%!     if(myid/=0) then
!%%!        allocate(kk(3,nkp),wtk(nkp),STAT=stat)
!%%!        if(stat/=0) call cq_abort('FindEvals: couldnt allocate kpoints',nkp)
!%%!        call reg_alloc_mem(area_general,4*nkp,type_dbl)
!%%!     end if
!%%!     call gcopy(kk,3,nkp)
!%%!     call gcopy(wtk,nkp)
!%%!     if(iprint_init>=5.AND.myid/=0) then
!%%!        write(io_lun,*) 'Proc: ',myid
!%%!        write(io_lun,5)
!%%!        do i=1,nkp
!%%!           write(io_lun,6) nkp,kk(1,i),kk(2,i),kk(3,i),wtk(i)
!%%!        end do
!%%!     end if
!%%!     return
!%%! 2   format(2x,'Block size (row, col): ',2i5)
!%%! 3   format(2x,'Proc grid (row, col): ',2i5)
!%%! 4   format(/2x,'***WARNING***',/,2x,'No Kpoints block found - defaulting to Gamma point')
!%%! 5   format(2x,i4,' Kpoints in Cartesian (inverse Angstrom) form: ')
!%%! 6   format(2x,'Kpt: ',i4,' : ',3f12.8,' Weight: ',f10.6)
!%%!   end subroutine readDiagInfo



end module initial_read
