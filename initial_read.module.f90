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
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2015/06/10 15:44 cor & dave
!!    Input parameters for band output
!!   2015/06/11 16:30 lat
!!    Added possibility of user to define filename for PAO/pseudo *ion
!!    Fixed EXX variable definitions
!!    Try to make blocks of variables clearer 
!!    Added experimental backtrace
!!   
!!  SOURCE
!!
module initial_read

  use timer_module,           only: start_timer, stop_timer, cq_timer
  use timer_module,           only: start_backtrace, stop_backtrace
  use timer_stdclocks_module, only: tmr_std_allocation

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: &
       RCSid = "$Id$"
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
  !!   datatypes, numbers, io_module, group_module, construct_module,
  !!   maxima_module, global_module, primary_module, atoms, dimens,
  !!   species_module, GenComms, pseudopotential_data
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
  !!    Bug fix (TM) - changed np to ind_part in loop to build atom
  !!    positions from make_prt.dat
  !!   11:02, 24/03/2003 drb 
  !!    Removed variables associated with old phi initialisation
  !!   14:52, 2003/06/09 dave
  !!    Removed atomic coordinate fiddling
  !!   14:59, 02/05/2005 dave 
  !!    Added calculation of number of electrons in cell
  !!   2007/05/08 17:00 dave
  !!    Added net charge on cell to number_of_bands for consistency
  !!   2011/09/19 L.Tong
  !!    Added calculation of ne_up_in_cell and ne_dn_in_cell for spin
  !!    polarised calculations
  !!   2011/09/29 15:06 M. Arita
  !!    Slight modification is added for DFT-D2
  !!   2011/12/11 L.Tong
  !!    Removed number_of_bands, it is obsolete, and use ne_in_cell
  !!    instead
  !!   2012/03/27 L.Tong
  !!   - Changed spin implementation
  !!   2014/02/04 M.Arita
  !!   - Added call for initial_read_aux for constraint-MD
  !!   2015/05/11 lat 
  !!   - Added optional total spin magnetization
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2015/06/10 15:44 cor & dave
  !!    - Input parameters for band output
  !!   2016/09/15 17:00 nakata
  !!    Added automatic setting of flag_one_to_one, atomf and nspin_SF
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!   2017/03/09 17:30 nakata
  !!    Changed to check consistence between flag_one_to_one and flag_Multisite
  !!   2017/07/11 16:36 dave
  !!    Bug fix (GitHub issue #36) to turn off basis optimisation if NSF=NPAO
  !!   2018/05/11 10:26 dave
  !!    Bug fix (GitHub issue #80) to set flag_SFcoeffReuse false if NSF=NPAO
  !!   2018/07/13 09:37 dave
  !!    Added test for missing coordinate file name
  !!  SOURCE
  !!
  subroutine read_and_write(start, start_L, inode, ionode,          &
                            vary_mu, mu, find_chdens, &
                            read_phi)

    use datatypes
    use numbers
    use io_module,              only: read_mult,                       &
                                      read_atomic_positions,           &
                                      pdb_template, pdb_format,        &
                                      print_process_info
    use group_module,           only: parts, part_method, HILBERT,     &
                                      PYTHON
    use construct_module,       only: init_group, init_primary
    use maxima_module,          only: maxpartsproc, maxatomsproc
    use global_module,          only: id_glob,x_atom_cell,y_atom_cell, &
                                      z_atom_cell, numprocs,           &
                                      iprint_init, rcellx, rcelly,     &
                                      rcellz, flag_old_partitions,     &
                                      nspin, spin_factor,              &
                                      flag_fix_spin_population,        &
                                      ne_in_cell, ne_spin_in_cell,     &
                                      ne_magn_in_cell,                 &
                                      ni_in_cell, area_moveatoms,      &
                                      io_lun, flag_only_dispersion,    &
                                      flag_basis_set, blips, PAOs,     &
                                      atomf, sf, paof,                 &
                                      flag_SpinDependentSF, nspin_SF,  &
                                      flag_Multisite,                  &
                                      flag_cdft_atom, flag_local_excitation, &
                                      flag_diagonalisation, flag_vary_basis, &
                                      flag_MDcontinue, flag_SFcoeffReuse
    use cdft_data, only: cDFT_NAtoms, &
                         cDFT_NumberAtomGroups, cDFT_AtomList
    use memory_module,          only: reg_alloc_mem, type_dbl
    use primary_module,         only: bundle, make_prim
    use dimens,                 only: r_super_x, r_super_y, r_super_z, &
                                      x_grid, y_grid, z_grid, r_h,     &
                                      set_dimensions, find_grid
    use block_module,           only: in_block_x, in_block_y,          &
                                      in_block_z
    use species_module,         only: n_species, species, charge,      &
                                      non_local_species,               &
                                      nsf_species, npao_species,       &
                                      natomf_species
    use GenComms,               only: my_barrier, cq_abort
    use pseudopotential_data,   only: non_local, read_pseudopotential
    use pseudopotential_common, only: core_radius, pseudo_type, OLDPS, &
                                      SIESTA, ABINIT
    use pseudo_tm_module,       only: init_pseudo_tm
    use input_module,           only: fdf_string, fdf_out, io_close, leqi
    use force_module,           only: tot_force
    !use DiagModule,             only: diagon
    use constraint_module,      only: flag_RigidBonds,constraints
    use support_spec_format,    only: flag_one_to_one, symmetry_breaking, read_option
    use pao_format

    implicit none

    ! Passed variables
    logical           :: vary_mu, start, start_L, read_phi
    logical           :: find_chdens
    integer           :: inode, ionode
    real(double)      :: mu

    ! Local variables
    type(cq_timer)    :: backtrace_timer
    character(len=80) :: titles, def
    character(len=80) :: atom_coord_file
    character(len=80) :: part_coord_file

    integer      :: i, ierr, nnd, np, ni, ind_part, j, acz, prncpl
    integer      :: ncf_tmp, stat
    integer      :: count_SZ, count_SZP
    real(double) :: ecore_tmp
    real(double) :: HNL_fac
    real(double), dimension(8) :: occ_n

    ! for checking the sum of electrons of spin channels
    real(double) :: sum_elecN_spin

!****lat<$
    call start_backtrace(t=backtrace_timer,who='read_and_write',where=1,level=1)
!****lat>$

    ! read input data: parameters for run
    call read_input(start, start_L, titles, vary_mu, mu,&
                    find_chdens, read_phi,HNL_fac)

    ! Initialise group data for partitions and read in partitions and atoms
    call my_barrier()
    def = ' '
    atom_coord_file = fdf_string(80,'IO.Coordinates',def)
    if(leqi(def,atom_coord_file)) call cq_abort("No coordinate file specified: please set with IO.Coordinates")
    if ( pdb_format ) then
       pdb_template = fdf_string(80,'IO.PdbTemplate',atom_coord_file)
    else
       pdb_template = fdf_string(80,'IO.PdbTemplate',' ')
    end if    
    def = 'make_prt.dat'
    part_coord_file = fdf_string(80,'IO.Partitions',def)
    if (.not. flag_MDcontinue) then
      call read_atomic_positions(trim(atom_coord_file))
    else
      call read_atomic_positions('cq.position')
    end if
    if(iprint_init>4) call print_process_info()
    ! By now, we'll have unit cell sizes and grid cutoff
    call find_grid
    if(flag_diagonalisation) call readDiagInfo
    if (flag_RigidBonds) call read_input_aux(constraints)
    ! Now that input from Conquest_input is done, we will close the file
    call io_close(fdf_out)
    if(inode==ionode.AND.iprint_init>1) &
         write(io_lun,fmt='(10x,"Partitioning method: ",i2)') part_method
    call read_mult(inode-1,parts,part_coord_file)
    !if(iprint_init>4) &
    !     write(io_lun,fmt='(10x,"Proc: ",i4," done read_mult: ",2i5)') &
    !          inode,maxatomsproc,maxpartsproc
    if(allocated(tot_force)) then
       write (io_lun,fmt='(10x,"WARNING! Proc: ",i4," tot_force &
                          &already allocated: ",i7)') &
             size(tot_force)
       deallocate(tot_force)
    end if
    allocate(tot_force(3,ni_in_cell))
    call reg_alloc_mem(area_moveatoms,3*ni_in_cell,type_dbl)
    if(flag_local_excitation) then
       allocate(flag_cdft_atom(ni_in_cell), STAT=stat)
       flag_cdft_atom = 0
       !if(cDFT_NumberAtomGroups==1) flag_cdft_atom = 2
       do i=1,cDFT_NumberAtomGroups
          do j = 1, cDFT_NAtoms(i)
             flag_cdft_atom(cDFT_AtomList(i)%Numbers(j)) = i
          end do
       end do
    end if
    call my_barrier

    ! Create primary set for atoms: bundle of partitions
    call init_primary(bundle, maxatomsproc, maxpartsproc, .true.)

    call make_prim(parts, bundle, inode-1, id_glob, x_atom_cell, &
                   y_atom_cell, z_atom_cell, species)

    ! Skip the following statement if only calculating the dispersion
    if (flag_only_dispersion) return

    !if(iprint_init>4) write(io_lun,fmt='(10x,"Proc: ",i4," done primary")') inode
    ! Read pseudopotential data
    if(pseudo_type == OLDPS) then
       call read_pseudopotential(inode, ionode)
    elseif(pseudo_type == SIESTA.OR.pseudo_type==ABINIT) then
       !just for setting core_radius in pseudopotential_common
       ! allocation and deallocation will be performed.
       call init_pseudo_tm(ecore_tmp, ncf_tmp)
    else
       call cq_abort(' Pseudotype Error ', pseudo_type)
    endif
    !if(iprint_init>4) write(io_lun,fmt='(10x,"Proc: ",i4," done pseudo")') inode
    !
    !
    !
    ! Set flag_one_to_one and check the number of SFs
    if (flag_basis_set==blips) then
       flag_one_to_one = .false.
       atomf = sf
       nspin_SF = 1
       flag_SpinDependentSF = .false.
    else if (flag_basis_set==PAOs) then
       flag_one_to_one = .true. ! primitive SFs
       atomf = sf
       nspin_SF = 1
       ! Check PAOs are contracted or not
       do i = 1, n_species
          if (nsf_species(i).ne.npao_species(i)) flag_one_to_one = .false.
       enddo
       if (flag_one_to_one .and. read_option) then
          if (inode==ionode) write(io_lun,'(A/A)') &
             'Warning!! nsf and npao are the same but read SF coefficients', &
             '          so that flag_one_to_one is set to be .false.'
          flag_one_to_one = .false.
       endif
       if (flag_one_to_one .and. flag_Multisite) then
          call cq_abort("flag_Multisite is .true., but the number of SFs is the same as the number of PAOs.")
       endif
       if (flag_one_to_one) then
          flag_SpinDependentSF = .false. ! spin-dependent SFs will be available only for contracted SFs
          flag_SFcoeffReuse = .false.
       else
          atomf = paof
       end if
       if (flag_SpinDependentSF) nspin_SF = nspin
       if (flag_one_to_one.AND.flag_vary_basis) then
          write(io_lun,fmt='(/2x,"************")')
          write(io_lun,fmt='(2x,"* WARNING !*")')
          write(io_lun,fmt='(2x,"************"/)')
          write(io_lun,fmt='(2x,"You have set minE.VaryBasis T")')
          write(io_lun,fmt='(2x,"This is not possible when numbers of support functions and PAOs are the same")')
          write(io_lun,fmt='(2x,"Setting the flag to false and proceeding"/)')
          flag_vary_basis = .false.
       end if
       ! Check symmetry-breaking for contracted SFs
       if (.not.flag_one_to_one) then
          do i = 1, 1, n_species
             count_SZ  = 0
             count_SZP = 0
             do j = 0, pao(i)%greatest_angmom
                if(pao(i)%angmom(j)%n_zeta_in_angmom>0) then
                   count_SZP = count_SZP + (2*j+1) ! Accumulate no. of ang. mom. components
                   occ_n(:) = zero
                   do acz = 1, pao(i)%angmom(j)%n_zeta_in_angmom
                      prncpl = pao(i)%angmom(j)%prncpl(acz)
                      occ_n(prncpl) = occ_n(prncpl) + pao(i)%angmom(j)%occ(acz)
                   enddo
                   do prncpl = 1, 8
                      if (occ_n(prncpl).gt.zero) count_SZ = count_SZ + (2*j+1)
                   enddo
                endif
             enddo
             ! If number of support functions is less than total number of ang. mom. components (ignoring
             ! for now multiple zetas) then there is a formal problem with basis set: we require the user
             ! to set an additional flag to assert that this is really desired
             if (flag_Multisite) then
                ! multisite SFs are symmetry-breaking usually
                ! multisite SFs should be in SZ-size at present
                if (count_SZ.ne.nsf_species(i)) then
                   if(inode==ionode) write(io_lun,'(A,I3,A,I3)') "Number of multi-site SFs", nsf_species(i), &
                                                                 " is not equal to single-zeta size", count_SZ
                   call cq_abort("Basis set error for species ",i)
                endif
             else 
                if(count_SZP>nsf_species(i)) then
                   if(.NOT.symmetry_breaking.OR..NOT.read_option) then
                      if(inode==ionode) then
                         write(io_lun,fmt='("You have a major problem with your basis set.")')
                         write(io_lun,fmt='("There are less support functions than the minimal angular momentum")')
                         write(io_lun,fmt='("components.  Either increase number of support functions on species ",i4)') i
                         write(io_lun,fmt='("to",i4," or set flag BasisSet.SymmetryBreaking to T")') count_SZP
                         write(io_lun,fmt='("You must also specify the basis set coefficients.")')
                         write(io_lun,fmt='("Use read_pao_coeffs T and support_pao_file <filename>.")')
                      end if
                      call cq_abort("Basis set error for species ",i)
                   else if(symmetry_breaking.AND.read_option) then
                      write(io_lun,fmt='("You have a major problem with your basis set.")')
                      write(io_lun,fmt='("There are ",i4," support functions and",i4," angular momentum")') &
                          nsf_species(i),count_SZP
                      write(io_lun,fmt='("components.  But as BasisSet.SymmetryBreaking is set T we will continue.")')
                   end if
                end if
             end if
          enddo ! n_species
       endif ! flag_one_to_one
    endif ! flag_basis_set
    if (atomf==sf)   natomf_species(:) =  nsf_species(:)
    if (atomf==paof) natomf_species(:) = npao_species(:)
    if (iprint_init.ge.3 .and. inode==ionode) write(io_lun,*) 'flag_one_to_one (T/F): ',flag_one_to_one
    if (iprint_init.ge.3 .and. inode==ionode) then
       if(atomf==sf) then
          write(io_lun,fmt='(2x,"Primitive atom functions are the support functions"/)')
       else if(atomf==paof) then
          write(io_lun,fmt='(2x,"Primitive atom functions are the pseudo-atomic orbitals"/)')
       endif
    end if
    !
    !
    !
    ! Make number of bands in an astonishingly crude way
    ! number_of_bands = half * ne_in_cell !zero
    do i = 1, ni_in_cell
       ! number_of_bands = number_of_bands + half * charge(species(i))
       ne_in_cell = ne_in_cell + charge(species(i))
    end do
    !
    !
    if (nspin == 2 .and. ne_magn_in_cell > zero) then
       ne_spin_in_cell(1) = half * (ne_magn_in_cell + ne_in_cell)
       ne_spin_in_cell(2) = ne_spin_in_cell(1) - ne_magn_in_cell  
       !
    end if
    !
    ! Calculate the number of electrons in the spin channels
    ! if (flag_fix_spin_population) then
    sum_elecN_spin = ne_spin_in_cell(1) + ne_spin_in_cell(2)
    if (sum_elecN_spin /= ne_in_cell) then
       if (nspin == 2 .and. flag_fix_spin_population) then
          call cq_abort('read_and_write: sum of number of electrons &
                         &in spin channels is different from total &
                         &number of electrons. ', &
                         sum_elecN_spin, ne_in_cell)
       else
          ne_spin_in_cell(:) = half * ne_in_cell
       end if
    end if
    !
    ! 
    !
    ! Set up various lengths, volumes, reciprocals etc. for convenient use
    call set_dimensions(inode, ionode,HNL_fac, non_local, n_species, &
                        non_local_species, core_radius)

    ! write out some information on the run
    if (inode == ionode) &
         call write_info(titles, mu, vary_mu, find_chdens, read_phi, &
                         HNL_fac, numprocs)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='read_and_write')
!****lat>$

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
  !!   datatypes, global_module, atoms, dimens, species_module,
  !!   pseudopotential_data, GenComms, fdf, parse
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
  !!   2011/09/19 L.Tong
  !!    - Added flag_spin_polarisation and flag_fix_spin_population
  !!    - Added ne_up_in_cell and ne_dn_in_cell, for fixed spin pupulation
  !!      calculations.
  !!    - Added A_dn as mixing parameter for spin down component (and A is
  !!      for up) for spin polarised calculations. The default (if A_dn
  !!      is missing from input) is A_dn = A
  !!    - Added functional_lsda_pw92 for LSDA
  !!   2011/11/17 10:36 dave
  !!    Changes for new blip data
  !!   2011/12/12 17:27 dave
  !!    Added flag for analytic blip integrals
  !!   2012/03/27 L.Tong
  !!   - Changed spin implementation
  !!   2012/05/29 L.Tong
  !!   - Added Checks for functional types. Makes sure for spin non-
  !!     polarised calculations only the one that has spin implemented
  !!     can be used.
  !!   2012/06/24 L.Tong
  !!   - Added input for flag_dump_L, which controls whether matL is
  !!     dumped or not
  !!   2013/02/28 L.Tong
  !!   - Removed the cq_abort when user choses non-SC calculations
  !!     with spin polarised PBE functionals. Conquest now will
  !!     perform the calculations as usual, but will give warnings
  !!     about non-SC forces and set non-SC forces to zero, and carry
  !!     on. This allows the user to still perform non-SC single point
  !!     calculations with spin polarisation, if the forces are not
  !!     important.
  !!   2013/02/08 15:44 dave (with UT)
  !!   - Flags for DeltaSCF
  !!   2013/04/03 L.Tong
  !!   - Added user input parameters for updated Hilbert curve
  !!     partitioning method
  !!   2013/07/01 M.Arita
  !!   - Added flags and parameters for the efficient MD scheme
  !!   2013/08/20 M.Arita
  !!   - Add a flag for reusing T-matrix
  !!   2014/01/21 lat
  !!   - Added flags for EXX and iprint_exx
  !!   2015/06/08 lat
  !!   - Added experimental backtrace
  !!   2015/06/10 15:47 cor & dave
  !!    Flags for wavefunction (band) output
  !!   2015/06/10 16:01 dave
  !!    Bug fix: added default value for r_exx when EXX not set (prevents problems)
  !!   2015/06/19 14:28 dave
  !!    Changed index of ChemicalSpeciesLabel to use number read from input file
  !!   2015/11/09 17:17 dave with TM and NW (Mizuho)
  !!    Added read for neutral atom flag
  !!   2015/11/24 08:30 dave
  !!    Removed flag_old_ewald (now redundant)
  !!   2016/01/28 16:44 dave
  !!    Updated module name to ion_electrostatic
  !!   2016/02/05 08:31 dave
  !!    Changed default pseudopotential to Siesta (necessary for now)
  !!   2016/08/05 10:08 dave
  !!    Added gap threshold parameter for initial estimation of partition numbers
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!   2017/03/14 SYM (with dave)
  !!    Fix lack of index on InvSRange when read in
  !!   2017/04/05 18:00 nakata
  !!    Added charge_up, charge_dn and flag_readAtomicSpin to initialise spin from input file
  !!   2017/04/24 dave
  !!    Changed pseudopotential output label to HAMANN (replacing ABINIT)
  !!   2017/05/09 dave
  !!    Removed restriction on L-matrix re-use and spin
  !!   2017/08/30 jack baker & dave
  !!    Adding parameters for simulation cell optimisation
  !!   2017/10/15 dave & an
  !!    Reading atomic spins: small tweak to test on net spin (uses abs and RD_ERR)
  !!   2017/10/24 zamaan
  !!    added thermostat flags for refactored NVT
  !!   2018/01/19 dave
  !!    Added test for lmax_ps and NA projector functions
  !!   2018/01/22 tsuyoshi (with dave)
  !!    Added new parameters:
  !!     - threshold for resetting charge density to atomic density
  !!     - time_max (General.MaxTime) to stop run after specified number of seconds in check_stop
  !!   2018/02/26 12:22 dave
  !!    Bug fix: turn off mixed L-SCF with diagonalisation (was breaking force tests)
  !!    Second issue (14:08) set SC.MinIters to 0 as was breaking force tests; also tweaked
  !!    linear mixing end default
  !!   2018/03/09 zamaan
  !!    Set some flags to true when restarting MD (loading matrices)
  !!   2018/04/12 22:00 nakata
  !!    Bug fix: turn on flag_LFD_MD_UseAtomicDensity only when flag_SFcoeffReuse is .false.
  !!   2018/04/24 zamaan
  !!    Changed AtomMove.FixCentreOfMass default to .true. except when running
  !!    FIRE qMD. minE.GlobalTolerance now defaults to .false. for MD.
  !!   2018/4/25 zamaan
  !!    added md_tdep flag for TDEP output dump
  !!   2018/05/17 12:54 dave with Ayako Nakata
  !!    flag_InitialAtomicSpin now comes from density_module (not global)
  !!   2018/07/16 16:02 dave
  !!    Bug fix: only read RadiusMS when flag_Multisite is true
  !!  TODO
  !!   Fix reading of start flags (change to block ?) 10/05/2002 dave
  !!   Fix rigid shift 10/05/2002 dave
  !!  SOURCE
  !!
  subroutine read_input(start, start_L, titles, vary_mu,&
                        mu, find_chdens, read_phi,HNL_fac)

    use datatypes
    use numbers
    use units
    use global_module, only: iprint, flag_vary_basis,                  &
                             flag_self_consistent,                     &
                             flag_precondition_blips,                  &
                             flag_fix_spin_population,                 &
                             flag_fractional_atomic_coords,            &
                             flag_old_partitions, ne_in_cell,          &
                             max_L_iterations, flag_read_blocks,       &
                             runtype, restart_LorK, restart_rho,          &
                             flag_basis_set, blips, PAOs,              &
                             flag_test_forces, UseGemm,                &
                             flag_fractional_atomic_coords,            &
                             flag_old_partitions, ne_in_cell,          &
                             ne_spin_in_cell, nspin, spin_factor,      &
                             ne_magn_in_cell,                          &
                             max_L_iterations,    &
                             flag_reset_dens_on_atom_move,             &
                             flag_continue_on_SC_fail, iprint_init,    &
                             iprint_mat, iprint_ops, iprint_DM,        &
                             iprint_SC, iprint_minE, iprint_time,      &
                             iprint_MD, iprint_index, iprint_gen,      &
                             iprint_pseudo, iprint_basis, iprint_exx,  &
                             iprint_intgn,iprint_MDdebug,area_general, &
                             global_maxatomspart, load_balance,        &
                             many_processors, flag_assign_blocks,      &
                             io_lun, flag_pulay_simpleStep,            &
                             flag_Becke_weights,                       &
                             flag_Becke_atomic_radii,                  &
                             flag_global_tolerance, flag_mix_L_SC_min, &
                             flag_onsite_blip_ana, flag_read_velocity, &
                             flag_quench_MD, temp_ion, numprocs,       &
                             flag_dft_d2, flag_only_dispersion,        &
                             flag_perform_cDFT, flag_analytic_blip_int,&
                             flag_vdWDFT, vdW_LDA_functional,          &
                             flag_dump_L, flag_DeltaSCF, dscf_source_level, &
                             dscf_target_level, dscf_source_spin,      &
                             dscf_target_spin, dscf_source_nfold,      &
                             dscf_target_nfold, flag_local_excitation, dscf_HOMO_thresh,   &
                             dscf_LUMO_thresh, dscf_HOMO_limit, dscf_LUMO_limit,           &
                             flag_MDcontinue,flag_MDdebug,flag_MDold,  &
                             flag_thermoDebug, flag_baroDebug, &
                             flag_LmatrixReuse,flag_TmatrixReuse,flag_SkipEarlyDM,McWFreq, &
                             restart_T,restart_X,flag_XLBOMD,flag_propagateX,              &
                             flag_propagateL,flag_dissipation,integratorXL, flag_FixCOM,   &
                             flag_exx, exx_alpha, exx_scf, exx_scf_tol, exx_siter,         &
                             flag_out_wf,flag_out_wf_by_kp,max_wf,out_wf,wf_self_con, flag_fire_qMD, &
                             flag_write_DOS, flag_write_projected_DOS, &
                             E_DOS_min, E_DOS_max, sigma_DOS, n_DOS, E_wf_min, E_wf_max, flag_wf_range_Ef, &
                             mx_temp_matrices, flag_neutral_atom, flag_diagonalisation, &
                             flag_SpinDependentSF, flag_Multisite, flag_LFD, flag_SFcoeffReuse, &
                             flag_opt_cell, cell_constraint_flag, cell_en_tol
    use dimens, only: r_super_x, r_super_y, r_super_z, GridCutoff,    &
                      n_grid_x, n_grid_y, n_grid_z, r_h, r_c,         &
                      RadiusSupport, RadiusAtomf, RadiusMS, RadiusLD, &
                      NonLocalFactor, InvSRange,                      &
                      min_blip_sp, flag_buffer_old, AtomMove_buffer,  &
                      r_dft_d2, r_exx
    use block_module, only: in_block_x, in_block_y, in_block_z, &
                            blocks_raster, blocks_hilbert
    use species_module, only: species_label, charge, mass, n_species,  &
                              charge, charge_up, charge_dn,            &
                              species, ps_file, ch_file, phi_file,     &
                              nsf_species, nlpf_species, npao_species, &
                              non_local_species, type_species,         &
                              species_file, species_from_files
    use GenComms,   only: gcopy, my_barrier, cq_abort, inode, ionode
    !use DiagModule, only: diagon
    use energy,     only: flag_check_Diag
    use DMMin,      only: maxpulayDMM, maxpulaystepDMM, minpulaystepDMM, &
                          LinTol_DMM, n_dumpL
!TM
    use pseudopotential_common, only: pseudo_type, OLDPS, SIESTA, &
         STATE, ABINIT, flag_angular_new, &
         flag_neutral_atom_projector, maxL_neutral_atom_projector, numN_neutral_atom_projector
    use SelfCon, only: A, flag_linear_mixing, EndLinearMixing, q0, q1,&
                       n_exact, maxitersSC, maxearlySC, maxpulaySC,   &
                       atomch_output, flag_Kerker, flag_wdmetric, minitersSC
    use atomic_density,  only: read_atomic_density_file, &
         atomic_density_method
    use density_module, only: flag_InitialAtomicSpin
    use S_matrix_module, only: InvSTolerance, InvSMaxSteps,&
                               InvSDeltaOmegaTolerance
    use blip,          only: blip_info, init_blip_flag, alpha, beta
    use maxima_module, only: maxnsf, lmax_ps
    use control,       only: MDn_steps, MDfreq, MDtimestep, MDcgtol, CGreset
    use ion_electrostatic,  only: ewald_accuracy
    use minimise,      only: UsePulay, n_L_iterations,          &
                             n_support_iterations, L_tolerance, &
                             sc_tolerance, energy_tolerance,    &
                             expected_reduction
    use pao_format,          only: kcut, del_k
    use support_spec_format, only: flag_paos_atoms_in_cell,          &
                                   read_option, symmetry_breaking,   &
                                   support_pao_file, TestBasisGrads, &
                                   TestTot, TestBoth, TestS, TestH
    use read_pao_info,     only: pao_info_file, pao_norm_flag
    use read_support_spec, only: support_spec_file, &
                                 flag_read_support_spec
    use test_force_module, only: flag_test_all_forces,           &
                                 flag_which_force, TF_direction, &
                                 TF_atom_moved, TF_delta
    use io_module, only: pdb_format, pdb_altloc, append_coords,  &
                         pdb_output, banner, get_file_name, time_max
    use group_module,     only: part_method, HILBERT, PYTHON
    use energy,           only: flag_check_DFT
    use H_matrix_module,  only: locps_output, locps_choice
    use pao_minimisation, only: InitStep_paomin
    use timer_module,     only: time_threshold,lun_tmr, TimingOn, &
                                TimersWriteOut, BackTraceOn
    use input_module!, only: load_input, input_array, block_start,
    ! block_end, fd
    use cdft_data, only: cDFT_Type, cDFT_MaxIterations, cDFT_NAtoms, &
                         cDFT_Target, cDFT_Tolerance,                &
                         cDFT_NumberAtomGroups, cDFT_AtomList,       &
                         cDFT_BlockLabel, cDFT_Vc
    use sfc_partitions_module, only: n_parts_user, average_atomic_diameter, gap_threshold
    use XLBOMD_module,         only: XLInitFreq,maxitersDissipation,kappa
    use constraint_module,     only: flag_RigidBonds,constraints,SHAKE_tol, &
                                     RATTLE_tol,maxiterSHAKE,maxiterRATTLE, &
                                     const_range,n_bond
    use exx_types, only: exx_scheme, exx_mem, exx_overlap, exx_alloc,       &
                         exx_cartesian, exx_radius, exx_hgrid, exx_psolver, &
                         exx_debug, exx_Kij, exx_Kkl, p_scheme
    use multisiteSF_module, only: flag_MSSF_smear, MSSF_Smear_Type,                      &
                                  MSSF_Smear_center, MSSF_Smear_shift, MSSF_Smear_width, &
                                  flag_LFD_ReadTVEC, LFD_TVEC_read,                      &
                                  LFD_kT, LFD_ChemP, flag_LFD_useChemPsub,               &
                                  flag_LFD_minimise, LFD_ThreshE, LFD_ThreshD,           &
                                  LFD_Thresh_EnergyRise, LFD_max_iteration,              &
                                  flag_LFD_MD_UseAtomicDensity
    use control,    only: md_ensemble
    use md_control, only: md_tau_T, md_n_nhc, md_n_ys, md_n_mts, md_nhc_mass, &
                          md_target_press, md_baro_type, md_tau_P, &
                          md_thermo_type, md_bulkmod_est, md_box_mass, &
                          flag_write_xsf, md_cell_nhc, md_nhc_cell_mass, &
                          md_calc_xlmass, md_berendsen_equil, &
                          md_tau_T_equil, md_tau_P_equil, md_p_drag, md_t_drag
    use md_model,   only: md_tdep
    use move_atoms,         only: threshold_resetCD, flag_stop_on_empty_bundle
    use Integrators, only: fire_alpha0, fire_f_inc, fire_f_dec, fire_f_alpha, fire_N_min, &
         fire_N_max, fire_max_step, fire_N_below_thresh
    use XC, only : flag_functional_type, functional_hartree_fock, functional_hyb_pbe0

    implicit none

    ! Passed variables
    logical           :: vary_mu, find_chdens
    logical           :: start, start_L, read_phi
    real(double)      :: mu, HNL_fac
    character(len=80) :: titles

    ! Local variables
    !type(block), pointer :: bp      ! Pointer to a block read by fdf
    !type(parsed_line), pointer :: p ! Pointer to a line broken into tokens by fdf
    type(cq_timer)    :: backtrace_timer
    type(cq_timer)    :: backtrace_timer2
    type(cq_timer)    :: backtrace_timer3


    integer           :: i, j, k, lun, stat
    integer           :: i_LFD, j_LFD, k_LFD, n_LFD, line_LFD
    character(len=20) :: def, tmp2
    character(len=80) :: coordfile, timefile, timefileroot
    character(len=10) :: basis_string, part_mode
    character(len=6)  :: method ! To find whether we diagonalise or use O(N)
    character(len=5)  :: ps_type !To find which pseudo we use
    character(len=8)  :: tmp
    logical           :: new_format
    !logical, external :: leqi
    real(double)      :: r_t, max_rc
    logical :: flag_ghost, find_species

    ! spin polarisation
    logical :: flag_spin_polarisation
    real(double) :: sum_elecN_spin

    ! Set defaults
    vary_mu  = .true.
    read_phi = .false.

    !*** WHO READS AND BROADCASTS ? ***!

!****lat<$
    call start_backtrace(t=backtrace_timer,who='read_input',where=1,level=2)
!****lat>$

    ! Test for existence of input file
    if (inode == ionode) then
       stat = 0
       open(unit=20,file='Conquest_input',iostat=stat,status='old')
       if (stat /= 0) call cq_abort("We need Conquest_input to run !")
       close(20)
    end if
    ! Start fdf reading from the file Conquest_input
    call load_input
    ! Here we need to read output file name, and open it.
    if (fdf_boolean('IO.WriteOutToFile',.true.)) then
       if (inode == ionode) then
          coordfile = fdf_string(80,'IO.OutputFile',"Conquest_out")
          call io_assign(io_lun)
          open(unit=io_lun,file=coordfile,iostat=stat)
          if (stat /= 0) &
               call cq_abort("Failed to open Conquest output file", stat)
       end if
       call gcopy(io_lun)
    else
       io_lun = 6
    end if
!!$
!!$
!!$
!!$  T I M E R S   &   O U T P U T
!!$
!!$
!!$
    BackTraceOn = fdf_boolean('IO.BackTraceOn',.false.)                         
    ! Where will the timers will be written?
    TimingOn = fdf_boolean('IO.TimingOn',.false.)                          ! tmr_rmv001
    if(TimingOn) then                                                      ! tmr_rmv001
       if(fdf_boolean('IO.TimeAllProcessors',.false.)) then                ! tmr_rmv001
          TimersWriteOut = .true.                                          ! tmr_rmv001
          timefileroot = fdf_string(80,'IO.TimeFileRoot',"time")           ! tmr_rmv001
          call get_file_name(timefileroot,numprocs,inode,timefile)         ! tmr_rmv001
          call io_assign(lun_tmr)                                          ! tmr_rmv001
          open(unit=lun_tmr,file=timefile,iostat=stat)                     ! tmr_rmv001
          if(stat/=0) call cq_abort("Failed to open time file",stat)       ! tmr_rmv001
       else                                                                ! tmr_rmv001
          TimersWriteOut = .false.                                         ! tmr_rmv001
          if(inode==ionode) then                                           ! tmr_rmv001
             TimersWriteOut = .true.                                       ! tmr_rmv001
             if(fdf_boolean('IO.WriteTimeFile',.true.)) then               ! tmr_rmv001
                timefileroot = fdf_string(80,'IO.TimeFileRoot',"time")     ! tmr_rmv001
                call get_file_name(timefileroot,numprocs,inode,timefile)   ! tmr_rmv001
                call io_assign(lun_tmr)                                    ! tmr_rmv001
                open(unit=lun_tmr,file=timefile,iostat=stat)               ! tmr_rmv001
                if(stat/=0) call cq_abort("Failed to open time file",stat) ! tmr_rmv001
             else                                                          ! tmr_rmv001
                lun_tmr = io_lun                                           ! tmr_rmv001
             end if ! Output to file                                       ! tmr_rmv001
          end if ! ionode switch                                           ! tmr_rmv001
       end if                                                              ! tmr_rmv001
    end if                                                                 ! tmr_rmv001
!!$ 
!!$ **<lat>** 
!!$ I think we may need (std) output file to be open before this point
!!$ at the very beginning in the main.f90 
!!$
!!$
!!$
!!$  I / O
!!$
!!$
    if(inode==ionode) call banner
    new_format = fdf_boolean('IO.NewFormat',.true.)
    if(new_format) then
       def = ' '
       iprint        = fdf_integer('IO.Iprint',       0     )
       iprint_init   = fdf_integer('IO.Iprint_init',  iprint)
       iprint_mat    = fdf_integer('IO.Iprint_mat',   iprint)
       iprint_ops    = fdf_integer('IO.Iprint_ops',   iprint)
       iprint_DM     = fdf_integer('IO.Iprint_DM',    iprint)
       iprint_SC     = fdf_integer('IO.Iprint_SC',    iprint)
       iprint_minE   = fdf_integer('IO.Iprint_minE',  iprint)
       iprint_MD     = fdf_integer('IO.Iprint_MD',    iprint)
       iprint_MDdebug= fdf_integer('IO.Iprint_MDdebug',1    )
       iprint_index  = fdf_integer('IO.Iprint_index', iprint)
       iprint_gen    = fdf_integer('IO.Iprint_gen'  , iprint)
       iprint_pseudo = fdf_integer('IO.Iprint_pseudo',iprint)
       iprint_basis  = fdf_integer('IO.Iprint_basis' ,iprint)
       iprint_intgn  = fdf_integer('IO.Iprint_intgn' ,iprint)
       iprint_time   = fdf_integer('IO.Iprint_time',  iprint)
       iprint_exx    = fdf_integer('IO.Iprint_exx',   iprint)
       !
       !
       flag_dump_L   = fdf_boolean('IO.DumpL',           .true. )
       locps_output  = fdf_boolean('IO.LocalPotOutput',  .false.)
       locps_choice  = fdf_integer('IO.LocalPotChoice',   8     )
       atomch_output = fdf_boolean('IO.AtomChargeOutput',.false.)
       !
       !
       tmp = fdf_string(8,'General.MemoryUnits','MB')
       if( leqi(tmp(1:2),     'kB') ) then
          m_units  = kbytes
          mem_conv = kB
       else if( leqi(tmp(1:2),'MB') ) then
          m_units  = mbytes
          mem_conv = MB
       else if( leqi(tmp(1:2),'GB') ) then
          m_units  = gbytes
          mem_conv = GB
       end if
       !
       !
       time_threshold = fdf_double('General.TimeThreshold',0.001_double)

       time_max  = fdf_double('General.MaxTime', zero)

       ! Read run title
       titles = fdf_string(80,'IO.Title',def)
       ! Is this a restart run ? **NB NOT AVAILABLE RIGHT NOW**
       start = fdf_boolean('General.NewRun',.true.)
       if(start) then
          start_L     = .true.
          read_option = .false.  ! N.B. Original had option for a blip model - add this later
       else ! This option loads data from a file - we need to add the option to search for this
          call cq_abort('read_input: you may not select restart for a run just now')
       endif
       init_blip_flag = fdf_string(10,'Basis.InitBlipFlag','pao')
       restart_rho    = fdf_boolean('General.LoadRho', .false.)
       restart_T      = fdf_boolean('General.LoadInvS',.false.)
       ! Is there a net charge on the cell ?
       ne_in_cell     = fdf_double('General.NetCharge',zero)
       ! Read coordinates file 
       flag_fractional_atomic_coords = &
            fdf_boolean('IO.FractionalAtomicCoords',.true.)
       call my_barrier()
       !
       !
       ! Spin polarised calculation
       ! - if we are doing spin polarised calculations or not
       flag_spin_polarisation   = fdf_boolean('Spin.SpinPolarised', .false.) 
       flag_fix_spin_population = fdf_boolean('Spin.FixSpin',       .false.)
       if (flag_spin_polarisation) then
          nspin       = 2
          spin_factor = one
          ne_spin_in_cell(1) = fdf_double('Spin.NeUP', zero)
          ne_spin_in_cell(2) = fdf_double('Spin.NeDN', zero)
          ne_magn_in_cell    = fdf_double('Spin.Magn', zero)
       end if
       !
       !
       !blip_width = fdf_double('blip_width',zero)
       !support_grid_spacing = fdf_double('support_grid_spacing',zero)

       ! Default to 25 Ha or 50 Ry cutoff for grid
       GridCutoff = fdf_double('Grid.GridCutoff',25.0_double)
       ! Grid points
       n_grid_x   = fdf_integer('Grid.PointsAlongX',0)
       n_grid_y   = fdf_integer('Grid.PointsAlongY',0)
       n_grid_z   = fdf_integer('Grid.PointsAlongZ',0)
       ! Number of different iterations - not well defined
       n_L_iterations       = fdf_integer('DM.LVariations',50)
       max_L_iterations     = n_L_iterations
       n_dumpL              = fdf_integer('DM.n_dumpL',n_L_iterations+1)
       n_support_iterations = fdf_integer('minE.SupportVariations',10)
       ! Initial expected drop in energy
       expected_reduction   = fdf_double('minE.ExpectedEnergyReduction',zero)
       ! Keep for fixed potential calculations
       mu = fdf_double('DM.mu',zero)
       if(fdf_defined('DM.ConstantMu')) then
          vary_mu = .false.
       else
          vary_mu = .true.
       end if
!!$
!!$
!!$
!!$
!!$
       ! Radii - again, astonishingly badly named
!****lat<$
       ! it seems that the initialisation of r_h here
       ! is obselete due to 'r_h = max(r_h,RadiusSupport(n))'
       ! in dimens.module.f90 is this what we want ?
!****lat>$
       r_h     = fdf_double('hamiltonian_range',zero)
       r_c     = fdf_double('DM.L_range',       one )
       r_t     = fdf_double('InvSRange',        r_h )
       HNL_fac = fdf_double('non_local_factor', zero)
       ! Exponents for initial gaussian blips
       alpha = fdf_double('Basis.SGaussExponent',zero)
       beta  = fdf_double('Basis.PGaussExponent',zero)
       ! Tolerance on minimisation
       energy_tolerance = fdf_double('minE.EnergyTolerance',1.0e-5_double)
       UsePulay         = fdf_boolean('Basis.UsePulayForPAOs',.false.)
       ! Sizes of integration grid blocks
       in_block_x = fdf_integer('Grid.InBlockX',4)
       in_block_y = fdf_integer('Grid.InBlockY',4)
       in_block_z = fdf_integer('Grid.InBlockZ',4)
       ! Solution method - O(N) or diagonalisation ?
       method = fdf_string(6,'DM.SolutionMethod','ordern') ! Default is O(N)
       if(leqi(method,'diagon')) then
          flag_diagonalisation = .true.
          flag_check_Diag = .true. 
       else
          flag_diagonalisation = .false.
          flag_check_Diag = .false.
       end if
       ! Read basis set
       basis_string = fdf_string(10,'Basis.BasisSet','PAOs')
       if(leqi(basis_string,'blips')) then
          flag_basis_set = blips
       else if(leqi(basis_string,'PAOs')) then
          flag_basis_set = PAOs
       end if
       read_option            = fdf_boolean('Basis.LoadCoeffs',           .false.)
       flag_onsite_blip_ana   = fdf_boolean('Basis.OnsiteBlipsAnalytic',  .true. )
       flag_analytic_blip_int = fdf_boolean('Basis.AnalyticBlipIntegrals',.false.)
       !
       !
!       find_chdens            = fdf_boolean('SC.MakeInitialChargeFromK',.false.)
       flag_Becke_weights     = fdf_boolean('SC.BeckeWeights',          .false.)
       flag_Becke_atomic_radii= fdf_boolean('SC.BeckeAtomicRadii',      .false.)
       ! Number of species
       n_species = fdf_integer('General.NumberOfSpecies',1)
       call allocate_species_vars
       flag_angular_new = fdf_boolean('Basis.FlagNewAngular',.true.)
       ! Multisite support functions
       flag_Multisite = fdf_boolean('Basis.MultisiteSF', .false.)
       if (flag_Multisite .and. flag_basis_set==blips) &
          call cq_abort("Multi-site support functions are available only with PAOs.")
!!$
!!$
!!$
!!$
!!$
       ! Tweak DRB 2007/03/23: fix abinit as standard type
       ! Tweak DRB 2016/02/05: fix siesta as standard ! 
       ps_type = fdf_string(5,'General.PseudopotentialType','siest') 
       ! Write out pseudopotential type
       if(leqi(ps_type,'siest')) then
          if(inode==ionode.AND.iprint_init>0) &
               write(io_lun,fmt='(10x,"SIESTA pseudopotential will be used. ")')
          pseudo_type = SIESTA
       else if(leqi(ps_type,'plato').OR.leqi(ps_type,'haman')) then
          if(inode==ionode.AND.iprint_init>0) &
               write(io_lun,fmt='(10x,"HAMANN pseudopotential will be used. ")')
          pseudo_type = ABINIT
       else
          if(inode==ionode.AND.iprint_init>0) &
               write(io_lun,fmt='(10x,"OLD pseudopotential will be used. ")')
          pseudo_type = OLDPS
       endif
       if((.NOT.flag_angular_new).AND.&
          (pseudo_type==SIESTA.OR.pseudo_type==ABINIT).AND.&
          inode==ionode) then
          write(io_lun,&
                fmt='(10x,"Setting FlagNewAngular to T for Siesta/&
                      &Abinit pseudopotentials")') 
          flag_angular_new = .true.
       end if
       ! Should we use neutral atom potential ?
       flag_neutral_atom  = fdf_boolean('General.NeutralAtom',.true.) ! Default for now
       if( flag_neutral_atom ) then
          flag_neutral_atom_projector = fdf_boolean('General.NeutralAtomProjector',.false.)
          if( flag_neutral_atom_projector ) then
             maxL_neutral_atom_projector = &
                  fdf_integer('General.NeutralAtomProjectorMaxL',3)
             if(maxL_neutral_atom_projector>lmax_ps) lmax_ps = maxL_neutral_atom_projector
             allocate( numN_neutral_atom_projector(0:maxL_neutral_atom_projector ) )
             if(fdf_block('NeutralAtom.ProjectorNumbers')) then
                if(1+block_end-block_start<maxL_neutral_atom_projector+1) &
                     call cq_abort("Too few NA projector numbers; found/needed: ", &
                     1+block_end-block_start,maxL_neutral_atom_projector+1)
                do i=1,maxL_neutral_atom_projector+1
                   read(unit=input_array(block_start+i-1),fmt=*) j,numN_neutral_atom_projector(i-1)
                end do
                call fdf_endblock
             else
                write(io_lun,fmt='(10x,"No NA projector numbers specified; defaulting to 5 per l")')
                numN_neutral_atom_projector(:) = 5
             end if
          end if ! NA projector
       else
          flag_neutral_atom_projector = .false.
       end if
!!$
!!$
!!$
!!$
!!$
       ! Read, using old-style fdf_block, the information about the
       ! different species
       !if(fdf_block('ChemicalSpeciesLabel',lun)) then
       ! flag_ghost=.false.
       !   do i=1,n_species
       !      read(lun,*) j,mass(i),species_label(i)
       !      if(mass(i) < -RD_ERR) flag_ghost=.true.
       !      type_species(i)=i
       !   end do
       !end if
       !
       !**<lat>** added possibility to read *.ion from file specified by the user
       species_from_files = fdf_boolean('General.PAOFromFiles',.false.)
       if(fdf_block('ChemicalSpeciesLabel')) then
          flag_ghost=.false.
          if(1+block_end-block_start<n_species) & 
               call cq_abort("Too few species in ChemicalSpeciesLabel: ",&
                             1+block_end-block_start,n_species)
          do i=1,n_species
             ! **<lat>** Added the optional filename
             if ( species_from_files ) then
                read (unit=input_array(block_start+i-1),fmt=*) &
                     j, mass(j),       &
                     species_label(j), &
                     species_file(j)
             else
                read (unit=input_array(block_start+i-1),fmt=*) &
                     j, mass(j),       &
                     species_label(j)
             end if
             if(mass(j) < -RD_ERR) flag_ghost=.true.
             type_species(j)=j
          end do
          call fdf_endblock
       end if
       ! I (TM) now assume we will use type_species(j) = 0
       ! for vacancy sites.   Nov. 14, 2007 
       if(flag_ghost) then
          do i=1, n_species
             if(mass(i) < -RD_ERR) then
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
             elseif(mass(i) > RD_ERR) then
                type_species(i) = i
             else
                if(inode==ionode) &
                     write (io_lun,*) &
                           ' There are chemical species which have zero mass ',&
                           i,mass(i),species_label(i)
                type_species(i) = 0    !  Vacancy sites
             endif
          enddo
       endif
!!$
!!$
!!$
!!$
!!$
       ! Read charge, mass, pseudopotential and starting charge and
       ! blip models for the individual species
       maxnsf      = 0
       max_rc = zero
       min_blip_sp = 1.0e8_double
       do i=1,n_species
          charge(i)         = zero
          nsf_species(i)    = 0
          RadiusSupport(i)  = r_h
          RadiusAtomf(i)    = r_h
          RadiusMS(i)       = zero
          RadiusLD(i)       = zero
          NonLocalFactor(i) = HNL_fac
          InvSRange(i)      = r_t
          blip_info(i)%SupportGridSpacing = zero
          if(pseudo_type==SIESTA.OR.pseudo_type==ABINIT) non_local_species(i) = .true.
          ! This is new-style fdf_block
          !nullify(bp)
          !if(fdf_block(species_label(i),bp)) then
          !   do while(fdf_bline(bp,line)) ! While there are lines in the block
          if(fdf_block(species_label(i))) then
             charge(i)        = fdf_double ('Atom.ValenceCharge',zero)
             charge_up(i)     = fdf_double ('Atom.SpinNeUp',zero)
             charge_dn(i)     = fdf_double ('Atom.SpinNeDn',zero)
             sum_elecN_spin   = charge_up(i)+charge_dn(i)
             if (abs(sum_elecN_spin)>RD_ERR) then
                flag_InitialAtomicSpin = .true.
                if (abs(sum_elecN_spin-charge(i))>RD_ERR) &
                   call cq_abort('read_input: sum of number of electrons &
                                 &in spin channels is different from total &
                                 &number of electrons for this species ', &
                                 sum_elecN_spin,charge(i))
             endif
             nsf_species(i)   = fdf_integer('Atom.NumberOfSupports',0)
             RadiusSupport(i) = fdf_double ('Atom.SupportFunctionRange',r_h)
             RadiusAtomf(i)   = RadiusSupport(i) ! = r_pao for (atomf=paof) or r_sf for (atomf==sf)
             ! Added DRB 2018/07/16 for safety
             if(flag_Multisite) RadiusMS(i)      = fdf_double ('Atom.MultisiteRange',zero)
             RadiusLD(i)      = fdf_double ('Atom.LFDRange',zero)
             if (flag_Multisite) RadiusSupport(i) = RadiusAtomf(i) + RadiusMS(i)
             ! DRB 2016/08/05 Keep track of maximum support radius
             if(RadiusSupport(i)>max_rc) max_rc = RadiusSupport(i)
             ! SYM 2017/03/13 Fix no index on InvSRange
             InvSRange(i)     = fdf_double ('Atom.InvSRange',zero)
             if(InvSRange(i)<RD_ERR) InvSRange(i) = RadiusSupport(i)
             NonLocalFactor(i) = fdf_double('Atom.NonLocalFactor',HNL_fac)
             if(NonLocalFactor(i)>one.OR.NonLocalFactor(i)<zero) then
                if(inode==ionode) &
                     write(io_lun,&
                           fmt='(10x,"Warning: Atom.NonLocalFactor &
                                 &must lie between 0.0 and 1.0: ",f9.5)') &
                          NonLocalFactor(i)
                if(NonLocalFactor(i)>one)  NonLocalFactor(i) = one
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
                      call cq_abort('read_input: no local/non-local &
                                     &specification for pseudopotential !')
                   end if
                end if
             end if
             call fdf_endblock
          else
             call cq_abort("Failure to read data for species_label: "&
                           //species_label(i))
          end if
          if(nsf_species(i)==0) &
               call cq_abort("Number of supports not specified for species ",i)
          if(flag_basis_set==blips.AND.blip_info(i)%SupportGridSpacing<RD_ERR) &
               call cq_abort("Error: for a blip basis set you must &
                              &set SupportGridSpacing for all species")
          maxnsf = max(maxnsf,nsf_species(i))
          if(RadiusSupport(i)<RD_ERR) &
               call cq_abort("Radius of support too small for &
                              &species; increase SupportFunctionRange ",i)
       end do
       ! For ghost atoms
       if(flag_ghost) then
          do i=1, n_species
             if(type_species(i) < 0) then
                charge(i) = zero
             endif
          enddo
       endif
!!$        
!!$        
!!$        
!!$        
!!$        
       ! Set variables for multisite support functions
       if (flag_Multisite) then
          flag_LFD = fdf_boolean('Multisite.LFD', .true.)
          flag_MSSF_smear = fdf_boolean('Multisite.Smear', .false.)
          if (flag_MSSF_smear) then
             MSSF_Smear_type   = fdf_integer('Multisite.Smear.FunctionType', 1)
             MSSF_Smear_center = fdf_double('Multisite.Smear.Center', zero) ! if zero, will be changed to r_s later
             MSSF_Smear_shift  = fdf_double('Multisite.Smear.Shift', zero)
             MSSF_Smear_width  = fdf_double('Multisite.Smear.Width', 0.1_double)
          endif
       else
          flag_LFD          = .false.
          flag_LFD_minimise = .false.
          flag_MSSF_smear   = .false.
       endif
       ! For LFD
       if (flag_LFD) then
          LFD_ChemP = fdf_double('Multisite.LFD.ChemP',zero)
          LFD_kT    = fdf_double('Multisite.LFD.kT',   0.1_double)
          flag_LFD_useChemPsub = fdf_boolean('Multisite.LFD.UseChemPsub', .true.)
          flag_LFD_ReadTVEC    = fdf_boolean('Multisite.LFD.ReadTVEC',    .false.)
          if (flag_LFD_ReadTVEC) then
             allocate(LFD_TVEC_read(n_species,25,100))
             if (fdf_block('LFDTrialVector')) then
                if(inode==ionode) write(io_lun,'(10x,A)') 'Reading TVEC for LFD'
                line_LFD = 0
                do i=1,n_species
                   if (nsf_species(i).gt.25) call cq_abort("nsf_species must be less or equal to 25 for Multi-site SFs")
                   do j = 1,nsf_species(i)
                      line_LFD = line_LFD + 1
                      ! i_LFD = species, j_LFD = sf, n_LFD = npao
                      read (unit=input_array(block_start+line_LFD-1),fmt=*) &
                            i_LFD,j_LFD,n_LFD,(LFD_TVEC_read(i,j,k_LFD),k_LFD=1,n_LFD)
                      if (i_LFD.ne.i)   call cq_abort("Error reading LFD TrialVector")
                      if (j_LFD.ne.j)   call cq_abort("Error reading LFD TrialVector")
                      if (n_LFD.gt.100) call cq_abort("n_PAO per atom must be less than 100 for Multi-site SFs")
                   enddo
                enddo
                call fdf_endblock
             else
                call cq_abort("No LFDTrialVector label in the input file.")
             endif
          endif ! flag_LFD_ReadTVEC
          flag_LFD_minimise = fdf_boolean('Multisite.LFD.Minimise',.true.)
          if (flag_LFD_minimise) then
             LFD_threshE = fdf_double('Multisite.LFD.Min.ThreshE',1.0e-6_double)
             LFD_threshD = fdf_double('Multisite.LFD.Min.ThreshD',1.0e-6_double)
             LFD_Thresh_EnergyRise = fdf_double('Multisite.LFD.Min.ThreshEnergyRise',LFD_threshE*ten)
             LFD_max_iteration = fdf_integer('Multisite.LFD.Min.MaxIteration',50)
          endif
       endif
       flag_LFD_MD_UseAtomicDensity = .false.
!!$
!!$
!!$
!!$
!!$ 
       !blip_width = support_grid_spacing *
       !             fdf_double('blip_width_over_support_grid_spacing',four)
       if (leqi(runtype,'md')) then
         flag_global_tolerance = fdf_boolean('minE.GlobalTolerance', .true.)
       else
         flag_global_tolerance = fdf_boolean('minE.GlobalTolerance', .false.)
       end if
       L_tolerance           = fdf_double ('minE.LTolerance',    1.0e-7_double)
       sc_tolerance          = fdf_double ('minE.SCTolerance',   1.0e-6_double)
       maxpulayDMM           = fdf_integer('DM.MaxPulay',        5            )
       minpulaystepDMM       = fdf_double ('DM.MinPulayStepSize',1.0e-3_double)
       maxpulaystepDMM       = fdf_double ('DM.MaxPulayStepSize',1.0e-1_double)
       LinTol_DMM = fdf_double('DM.LinTol',0.1_double)
!!$
!!$
!!$
!!$
!!$
       ! Find out what type of run we're doing
       runtype             = fdf_string(20,'AtomMove.TypeOfRun',       'static')
       flag_buffer_old       = fdf_boolean('AtomMove.OldBuffer',        .false.)
       AtomMove_buffer       = fdf_double ('AtomMove.BufferSize',    4.0_double)
       flag_pulay_simpleStep = fdf_boolean('AtomMove.PulaySimpleStep',  .false.)
       CGreset               = fdf_boolean('AtomMove.ResetCG',          .false.)
       MDn_steps             = fdf_integer('AtomMove.NumSteps',     100        )
       MDfreq                = fdf_integer('AtomMove.OutputFreq',    50        )
       MDtimestep            = fdf_double ('AtomMove.Timestep',      0.5_double)
       MDcgtol               = fdf_double ('AtomMove.MaxForceTol',0.0005_double)
       flag_opt_cell         = fdf_boolean('AtomMove.OptCell',          .false.)
       cell_constraint_flag  = fdf_string(20,'AtomMove.OptCell.Constraint','none')
       cell_en_tol           = fdf_double('AtomMove.OptCell.EnTol',0.00001_double)
       flag_stop_on_empty_bundle = fdf_boolean('AtomMove.StopOnEmptyBundle',.false.)
       !
       flag_vary_basis       = fdf_boolean('minE.VaryBasis', .false.)
       if(.NOT.flag_vary_basis) then
          flag_precondition_blips = .false.
       else 
          flag_precondition_blips = fdf_boolean('minE.PreconditionBlips',.true.)
       end if
       InitStep_paomin      = fdf_double ('minE.InitStep_paomin',  5.0_double)
       flag_self_consistent = fdf_boolean('minE.SelfConsistent',      .true. )
       flag_mix_L_SC_min    = fdf_boolean('minE.MixedLSelfConsistent',.false.)
       ! DRB 2018/02/26 turn off mixed L-SCF with diagonalisation
       if(flag_mix_L_SC_min.AND.flag_diagonalisation) flag_mix_L_SC_min = .false.
       ! Tweak 2007/03/23 DRB Make Pulay mixing default
       flag_linear_mixing   = fdf_boolean('SC.LinearMixingSC',        .true. )
       A(1)                 = fdf_double ('SC.LinearMixingFactor', 0.5_double)
       ! 2011/09/19 L.Tong, for spin polarised calculation
       A(2)                 = fdf_double ('SC.LinearMixingFactor_SpinDown', A(1))
       ! end 2011/09/19 L.Tong
       EndLinearMixing = fdf_double ('SC.LinearMixingEnd',sc_tolerance*1e-4_double) ! DRB 2018/02/26
       flag_Kerker     = fdf_boolean('SC.KerkerPreCondition',  .false.)
       q0              = fdf_double ('SC.KerkerFactor',     0.1_double)
       flag_wdmetric   = fdf_boolean('SC.WaveDependentMetric', .false.)
       q1              = fdf_double ('SC.MetricFactor',     0.1_double)
       n_exact         = fdf_integer('SC.LateStageReset',   5         )
       flag_reset_dens_on_atom_move = fdf_boolean('SC.ResetDensOnAtomMove',.false.)
       flag_continue_on_SC_fail     = fdf_boolean('SC.ContinueOnSCFail',   .false.)
       maxitersSC      = fdf_integer('SC.MaxIters',50)
       minitersSC      = fdf_integer('SC.MinIters',0) ! Changed default 2->0 DRB 2018/02/26
       maxearlySC      = fdf_integer('SC.MaxEarly',3 )
       maxpulaySC      = fdf_integer('SC.MaxPulay',5 )
       read_atomic_density_file = &
                       fdf_string(80,'SC.ReadAtomicDensityFile','read_atomic_density.dat')
       ! Read atomic density initialisation flag
       atomic_density_method = fdf_string(10,'SC.AtomicDensityFlag','pao')
       ! When constructing charge density from K matrix at last step, we check the total 
       ! number of electrons. If the error of electron number (per total electron number) 
       ! is larger than the following value, we use atomic charge density. (in update_H)
       threshold_resetCD     = fdf_double('SC.Threshold.Reset',0.1_double)
!!$
!!$
!!$
!!$
!!$
       ! read wavefunction output flags
       mx_temp_matrices = fdf_integer('General.MaxTempMatrices',100)
       flag_out_wf=fdf_boolean('IO.outputWF',.false.)
       if (flag_out_wf) then
          if (flag_diagonalisation .and. leqi(runtype,'static')) then
             flag_out_wf_by_kp=fdf_boolean('IO.outputWF_by_kpoint',.false.)
             wf_self_con=.false.
             ! The user can either specify which bands explicitly
             max_wf=fdf_integer('IO.maxnoWF',0)
             if(max_wf>0) allocate(out_wf(max_wf))
             if (fdf_block('WaveFunctionsOut')) then
                if(1+block_end-block_start<max_wf) &
                 call cq_abort("Too few wf no in WaveFunctionsOut:"&
                              ,1+block_end-block_start,max_wf)
                do i=1,max_wf
                   read(unit=input_array(block_start+i-1),fmt=*) out_wf(i)
                end do
                call fdf_endblock
             end if
             ! Or specify an energy range
             E_wf_min = fdf_double('IO.min_wf_E',zero)
             E_wf_max = fdf_double('IO.max_wf_E',zero)
             ! Is the range relative to Ef (T) or absolute (F)
             flag_wf_range_Ef = fdf_boolean('IO.WFRangeRelative',.true.)
             if(flag_wf_range_Ef.AND.E_wf_min==zero.AND.E_wf_max==zero) then
                flag_out_wf = .false.
                flag_wf_range_Ef = .false.
                if(inode==ionode) write(io_lun,'(2x,"Setting IO.outputWF F as no bands range given")')
             end if
          else
              call cq_abort("Won't output WFs for Order(N) or non-static runs")
          end if
       end if
       ! DOS output
       flag_write_DOS = fdf_boolean('IO.writeDOS',.false.)
       if(flag_write_DOS) then
          if(flag_diagonalisation) then
             flag_write_projected_DOS = fdf_boolean('IO.write_proj_DOS',.false.)
             E_DOS_min = fdf_double('IO.min_DOS_E',zero)
             E_DOS_max = fdf_double('IO.max_DOS_E',zero)
             sigma_DOS = fdf_double('IO.sigma_DOS',0.001_double)
             n_DOS = fdf_integer('IO.n_DOS',201)
          else
             flag_write_DOS = .false.
             if(inode==ionode) write(io_lun,'(2x,"Setting IO.writeDOS F as solving O(N)")')
          end if
       end if
!!$
!!$
!!$
!!$
!!$
       ! cDFT flags
       flag_perform_cDFT = fdf_boolean('cDFT.Perform_cDFT',.false.)

       if(flag_perform_cDFT) then
          if(flag_multisite) call cq_abort('cDFT is not supported with multi-site SFs.')
          if(.NOT.flag_Becke_weights) then
             flag_Becke_weights = .true.
             write(io_lun,fmt='(2x,"Warning: we require Becke  weights for cDFT ! Setting to true")')
          end if
          cDFT_Type          = fdf_integer('cDFT.Type',2) ! Default to D-A charge diff
          cDFT_MaxIterations = fdf_integer('cDFT.MaxIterations',50)
          cDFT_Tolerance     = fdf_double ('cDFT.Tolerance',0.001_double)
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
!!$
!!$
!!$  D e l t a - S C F 
!!$
!!$
       ! DeltaSCF flags
       flag_DeltaSCF = fdf_boolean('minE.DeltaSCF',.false.)
       if(flag_DeltaSCF) then
          if(flag_multisite) call cq_abort('DeltaSCF is not supported with multi-site SFs.')
          dscf_source_level     = fdf_integer('DeltaSCF.SourceLevel',  0)
          dscf_target_level     = fdf_integer('DeltaSCF.TargetLevel',  0)
          dscf_source_spin      = fdf_integer('DeltaSCF.SourceChannel',1)
          dscf_target_spin      = fdf_integer('DeltaSCF.TargetChannel',1)
          dscf_source_nfold     = fdf_integer('DeltaSCF.SourceNFold',  1)
          dscf_target_nfold     = fdf_integer('DeltaSCF.TargetNFold',  1)
          flag_local_excitation = fdf_boolean('DeltaSCF.LocalExcitation',.false.)
          if(flag_local_excitation) then
             dscf_homo_limit    = fdf_integer('DeltaSCF.HOMOLimit',0)
             if(dscf_homo_limit<0) then
                dscf_homo_limit = -dscf_homo_limit
                write(io_lun,fmt='(2x,"Changed sign on dscf_HOMO_limit: ",i5)') dscf_homo_limit
             end if
             dscf_lumo_limit = fdf_integer('DeltaSCF.LUMOLimit',0)
             if((dscf_homo_limit==0).AND.(dscf_lumo_limit==0)) then
                flag_local_excitation = .false.
                write(io_lun,fmt='(2x,"You must set DeltaSCF.HOMOLimit or DeltaSCF.LUMOLimit for local excitation")')
             end if
             dscf_HOMO_thresh      = fdf_double ('DeltaSCF.HOMOThresh',0.5_double)
             dscf_LUMO_thresh      = fdf_double ('DeltaSCF.LUMOThresh',0.5_double)
             cDFT_NumberAtomGroups = fdf_integer('cDFT.NumberAtomGroups',1)
             allocate(cDFT_NAtoms(cDFT_NumberAtomGroups),cDFT_Target(cDFT_NumberAtomGroups), &
                  cDFT_AtomList(cDFT_NumberAtomGroups),cDFT_BlockLabel(cDFT_NumberAtomGroups))
             if(fdf_block('cDFT.AtomGroups')) then
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
       end if
!!$
!!$
!!$
!!$
       InvSTolerance = fdf_double('DM.InvSTolerance',1e-2_double)
       InvSDeltaOmegaTolerance = &
                       fdf_double('DM.InvSDeltaOmegaTolerance',0.0001_double)
       InvSMaxSteps = fdf_integer('DM.InvSMaxSteps',100)
       basis_string = fdf_string(10,'General.BlockAssign','Hilbert')
       if(leqi(basis_string,'Raster')) then
          flag_assign_blocks = blocks_raster
       else if(leqi(basis_string,'Hilbert')) then
          flag_assign_blocks = blocks_hilbert
       end if
       flag_read_blocks = fdf_boolean('Grid.ReadBlocks',.false.)
       ! Default value of 10^-10 Ha per atom
       ewald_accuracy   = fdf_double ('General.EwaldAccuracy',1.0e-10_double)
       UseGemm = fdf_boolean('MM.UseGemm',.false.)
       if(flag_ghost) then
          if(inode==ionode) write(io_lun,*) ' As ghost atoms are included, UseGemm must be false.'
          UseGemm = .false.
       endif

       flag_check_DFT     = fdf_boolean('General.CheckDFT',.false.)
       flag_quench_MD     = fdf_boolean('AtomMove.QuenchMD',.false.)
       flag_fire_qMD = fdf_boolean('AtomMove.FIRE',.false.)
       ! If we're doing MD, then the centre of mass should generally be fixed,
       ! except for FIRE quenched MD
       if (leqi(runtype, 'md')) then
          if (flag_fire_qMD) then
            flag_FixCOM = fdf_boolean('AtomMove.FixCentreOfMass', .false.)
          else
            flag_FixCOM = fdf_boolean('AtomMove.FixCentreOfMass', .true.)
          end if
       end if
       if(flag_fire_qMD) then
          fire_N_min         = fdf_integer('AtomMove.FireNMin',5)
          fire_N_max         = fdf_integer('AtomMove.FireNMaxSlowQMD',10)
          fire_alpha0        = fdf_double ('AtomMove.FireAlpha',0.1_double)
          fire_f_inc         = fdf_double ('AtomMove.Fire_f_inc',1.1_double)
          fire_f_dec         = fdf_double ('AtomMove.Fire_f_dec',0.5_double)
          fire_f_alpha       = fdf_double ('AtomMove.Fire_f_alpha',0.99_double)
          fire_max_step      = fdf_double ('AtomMove.FireMaxStep',2.0_double)
          fire_N_below_thresh= fdf_integer('AtomMove.FireNBelowThresh',3)
          temp_ion           = fdf_double ('AtomMove.IonTemperature',0.0_double)
       else
          temp_ion           = fdf_double ('AtomMove.IonTemperature',300.0_double)
       end if
!!$
!!$
!!$
!!$
       del_k = fdf_double('Basis.PaoKspaceOlGridspace',0.1_double)
       kcut  = fdf_double('Basis.PaoKspaceOlCutoff', 1000.0_double)
       flag_paos_atoms_in_cell = fdf_boolean(  'Basis.PAOs_StoreAllAtomsInCell',.true. )
!       flag_one_to_one         = fdf_boolean(  'Basis.PAOs_OneToOne',           .false.)  ! will be set automatically later
       symmetry_breaking       = fdf_boolean(  'Basis.SymmetryBreaking',        .false.)
       support_pao_file        = fdf_string(80,'Basis.SupportPaoFile',   'supp_pao.dat')
       pao_info_file           = fdf_string(80,'Basis.PaoInfoFile',           'pao.dat')
       pao_norm_flag           = fdf_integer(  'Basis.PaoNormFlag',              0     )
       TestBasisGrads          = fdf_boolean(  'Basis.TestBasisGradients',      .false.)
       TestTot                 = fdf_boolean(  'Basis.TestBasisGradTot',        .false.)
       TestBoth                = fdf_boolean(  'Basis.TestBasisGradBoth',       .false.)
       TestS                   = fdf_boolean(  'Basis.TestBasisGrad_S',         .false.)
       TestH                   = fdf_boolean(  'Basis.TestBasisGrad_H',         .false.)
       support_spec_file       = fdf_string(80,'Basis.SupportSpecFile',   'support.dat')
       flag_read_support_spec  = fdf_boolean(  'Basis.ReadSupportSpec',         .false.)
       flag_SpinDependentSF    = fdf_boolean(  'Basis.SpinDependentSF',         .false.) ! Spin-dependence of SFs
       !
       !
       flag_test_forces        = fdf_boolean('AtomMove.TestForces',   .false.)
       flag_test_all_forces    = fdf_boolean('AtomMove.TestAllForces',.true. )
       if(.NOT.flag_test_all_forces) then ! Test one force component
          flag_which_force = fdf_integer('AtomMove.TestSpecificForce',1)
          if(flag_which_force>8.OR.flag_which_force<0.AND.inode==ionode) &
               write (io_lun,fmt='(10x,"Warning ! TestSpecificForce &
                                  &must lie between 1 and 8: ",i3)') &
                     flag_which_force
       end if
       TF_direction = fdf_integer('AtomMove.TestForceDirection',1)
       if(TF_direction>3.OR.TF_direction<0.AND.inode==ionode) then
          write (io_lun,fmt='(10x,"Warning ! TestForceDirection must &
                              &lie between 1 and 3: ",i3)') TF_direction
          TF_direction = 1
       end if
       TF_atom_moved = fdf_integer('AtomMove.TestForceAtom',1)
       TF_delta      = fdf_double ('AtomMove.TestForceDelta',0.00001_double)
!!$
!!$
!!$
!!$       
       flag_functional_type = fdf_integer('General.FunctionalType', 3)   ! LDA PW92
!!$
!!$
!!$  E X X 
!!$ 
!!$
       if ( flag_functional_type == functional_hyb_pbe0 ) then
          flag_exx = .true.
          exx_siter = fdf_integer('EXX.StartAfterIter', 1 )
          exx_scf   = fdf_integer('EXX.MethodSCF',      0 )
          r_exx     = fdf_double ('EXX.Krange'   ,   zero )
          !
       else if ( flag_functional_type == functional_hartree_fock ) then
          flag_exx = .true.
          exx_scf  = fdf_integer('EXX.MethodSCF', 0)
          r_exx    = fdf_double ('EXX.Krange', zero)
          !
       else
          ! don't touch we need it because matX is setup in set_dimensions 
          ! whatever is flag_exx
          flag_exx = .false.
          exx_scf  =  fdf_integer('EXX.MethodSCF', -1)
          r_exx    =  one 
          !             
       end if
       !
       if ( flag_exx ) then
          ! Setup the exact-exchange mixing factor (0 < exx_alpha <1)
          if ( flag_functional_type == functional_hyb_pbe0 ) then
             exx_alpha  = fdf_double ('EXX.Alpha',0.25_double)
             !
          else if( flag_functional_type == functional_hartree_fock ) then
             exx_alpha  = fdf_double ('EXX.Alpha',1.00_double)
             !
          end if
          if ( exx_alpha < zero .or. alpha > one ) then
             call cq_abort('EXX: exact-exchange mixing value is not &
                 &reasonable ', exx_alpha)
             !
          end if
          ! To control accuracy during scf
          exx_scf_tol   = sc_tolerance
          ! Grid spacing for PAO discretisation in EXX
          exx_hgrid  = fdf_double ('EXX.GridSpacing',zero)
          exx_radius = fdf_double ('EXX.IntegRadius',0.00_double) 
          ! debug mode
          exx_Kij       = .true.
          exx_Kkl       = .true.
          exx_cartesian = .true. 
          exx_overlap   = .true. 
          exx_alloc     = .false.
          exx_psolver   = 'fftw'
          p_scheme      = 'pulay'
          exx_scheme    = 1
          exx_mem       = 1
          exx_debug     = .false.
       end if
!!$
!!$
!!$  U N I T S
!!$
!!$
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
!!$
!!$
!!$  P D B
!!$
!!$

       for_conv = en_conv/dist_conv
       pdb_format = fdf_boolean ('IO.PdbIn',.false.)
       pdb_altloc = fdf_string(1,'IO.PdbAltLoc',' ')
       pdb_output = fdf_boolean ('IO.PdbOut',.false.)
       !
       !
!!$
!!$
!!$  P A R T I T I O N    M E T H O D
!!$
!!$
       part_mode  = fdf_string(10,'General.PartitionMethod','Hilbert')
       if (leqi(part_mode(1:6),'File')) then
          part_method = PYTHON
       else if (leqi(part_mode(1:7), 'Hilbert')) then
          part_method = HILBERT
       else
          part_method = HILBERT
          if(inode==ionode) &
               write(io_lun,*) 'WARNING: Unrecognised partitioning method'
       end if
       ! For Hilbert curve sfc partitioning
       tmp2 = fdf_string(20,'General.LoadBalance','atoms')
       if (leqi(tmp2(1:10),'partitions')) then
          load_balance = 1
       else if (leqi(tmp2(1:5),'atoms')) then
          load_balance = 0
       else if (leqi(tmp2(1:16),'supportfunctions')) then
          load_balance = 2
       else
          load_balance = 0
       end if
       ! many_processors = fdf_boolean('General.ManyProcessors',.true.)
       global_maxatomspart = fdf_integer('General.MaxAtomsPartition', 34)
       n_parts_user(1)     = fdf_integer('General.NPartitionsX', 0)
       n_parts_user(2)     = fdf_integer('General.NPartitionsY', 0)
       n_parts_user(3)     = fdf_integer('General.NPartitionsZ', 0)
       if(n_parts_user(1)*n_parts_user(2)*n_parts_user(3)>0.AND.&
            n_parts_user(1)*n_parts_user(2)*n_parts_user(3)<numprocs) &
            call cq_abort("More processors than partitions! ",numprocs,n_parts_user(1)*n_parts_user(2)*n_parts_user(3))
       average_atomic_diameter = &
            fdf_double('General.AverageAtomicDiameter', 5.0_double)
       gap_threshold = fdf_double('General.GapThreshold',two*max_rc)
       ! end sfc partitioning
       !
       !
!!$
!!$
!!$  
!!$
!!$
       append_coords = fdf_boolean('AtomMove.AppendCoords',.true.)
       !
       !
!!$
!!$
!!$  D F T - D 2
!!$
!!$
       flag_dft_d2   = fdf_boolean('General.DFT_D2', .false.)                 ! for DFT-D2
       if (flag_dft_d2) r_dft_d2 = fdf_double('DFT-D2_range',23.0_double)     ! for DFT-D2
       flag_only_dispersion = fdf_boolean('General.only_Dispersion',.false.)  ! for DFT-D2
       ! vdW XC functional flags
       flag_vdWDFT = fdf_boolean('General.vdWDFT', .false.)
       if (flag_vdWDFT) then
          vdW_LDA_functional = fdf_integer('vdWDFT.LDAFunctionalType', 3)
       end if
       !
       !
!!$
!!$
!!$  M O L E C U L A R    D Y N A M I C S
!!$
!!$
       ! Basic settings for MD
       flag_MDold        = fdf_boolean('AtomMove.OldMemberUpdates',.false.)
       flag_MDdebug      = fdf_boolean('AtomMove.Debug',.false.)
       flag_MDcontinue   = fdf_boolean('AtomMove.RestartRun',.false.)
       flag_SFcoeffReuse = fdf_boolean('AtomMove.ReuseSFcoeff',.false.)
       flag_LmatrixReuse = fdf_boolean('AtomMove.ReuseL',.false.)
       flag_write_xsf    = fdf_boolean('AtomMove.WriteXSF', .true.)
       if (flag_LFD .and. .not.flag_SFcoeffReuse) then
          ! if LFD, use atomic density in default when we don't reuse SFcoeff
          flag_LFD_MD_UseAtomicDensity = fdf_boolean('Multisite.LFD.UpdateWithAtomicDensity',.true.)
       endif
       ! DRB 2017/05/09 Removing restriction (now implemented)
       !if(flag_spin_polarisation.AND.flag_LmatrixReuse) then
       !   call cq_abort("L matrix re-use and spin polarisation not implemented !")
       !end if

       flag_TmatrixReuse = fdf_boolean('AtomMove.ReuseInvS',.false.)
       flag_SkipEarlyDM  = fdf_boolean('AtomMove.SkipEarlyDM',.false.)
       McWFreq           = fdf_integer('AtomMove.McWeenyFreq',0)
       ! XL-BOMD
       flag_XLBOMD       = fdf_boolean('AtomMove.ExtendedLagrangian',.false.)

       ! zamaan 2018/03/03 I can't imagine a case where you would want to
       ! restart a MD run without loading the various matrices, so I'm 
       ! Defaulting some flags to true. Calling fdf to ensure that input.log
       ! remains consistent
       if (flag_MDcontinue) then
         flag_read_velocity = fdf_boolean('AtomMove.ReadVelocity',.true.)
         restart_LorK   = fdf_boolean('General.LoadL', .true.)
         find_chdens    = fdf_boolean('SC.MakeInitialChargeFromK',.true.)
         if (flag_XLBOMD) restart_X=fdf_boolean('XL.LoadX', .true.)
       else
         flag_read_velocity = fdf_boolean('AtomMove.ReadVelocity',.false.)
         restart_LorK   = fdf_boolean('General.LoadL', .false.)
         find_chdens    = fdf_boolean('SC.MakeInitialChargeFromK',.false.)
         if (flag_XLBOMD) restart_X=fdf_boolean('XL.LoadX', .false.)
       end if

       if (flag_XLBOMD) then
         kappa=fdf_double('XL.Kappa',2.0_double)
         if (kappa.GT.2.0_double) then
           kappa=2.0_double
           if (inode.EQ.ionode) &
             write (io_lun,'(2x,a)') "WARNING: kappa should not be larger than 2.0 ! &
                                     &Setting to 2.0 "
         endif
         if (.NOT.flag_LmatrixReuse) then
           flag_LmatrixReuse = .true.
           if (inode.EQ.ionode) write (io_lun,'(2x,a)') &
             "WARNING: XL-BOMD requires L-matrix be reused ! Setting to true "
         endif
         flag_propagateX = fdf_boolean('XL.PropagateX',.true.)
         flag_propagateL = fdf_boolean('XL.PropagateL',.false.)
         if (flag_propagateX .AND. flag_propagateL) then
           flag_propagateX = .false.
           if (inode.EQ.ionode) write (io_lun,'(2x,a)')      &
             "WARNING: we require X-matrix not be propagated &
             &when employing the original scheme ! Setting to&
             & false "
         endif
         XLInitFreq          = fdf_integer('XL.ResetFreq',0)
         flag_dissipation    = fdf_boolean('XL.Dissipation',.false.)
         maxitersDissipation = fdf_integer('XL.MaxDissipation',5)
         if (flag_dissipation) then
           if (maxitersDissipation.LT.1 .OR. maxitersDissipation.GT.9) then
             maxitersDissipation = 5
             if (inode.EQ.ionode) write (io_lun,'(2x,a)') &
               "WARNING: K must be 1 <= K <= 9 ! Setting to 5"
           endif
         endif
         integratorXL=fdf_string(20,'XL.Integrator','velocityVerlet')
         if (integratorXL.NE.'velocityVerlet' .AND. integratorXL.NE.'Verlet') then
           call cq_abort("Integrator for XL-BOMD must be either velocity Verlet &
                         &or Verlet")
         endif
         if (integratorXL.EQ.'Verlet' .AND. .NOT.flag_dissipation) then
           integratorXL='velocityVerlet'
           if (inode.EQ.ionode) write (io_lun,'(2x,a)') &
             "WARNING: integrator must be velocity Verlet when dissipation &
             &does not apply ! Setting to velocity Verlet "
         endif
       endif ! XL-BOMD
       ! Constraints
       flag_RigidBonds=fdf_boolean('AtomMove.RigidBonds', .false.)
       if (flag_RigidBonds) then
         constraints%filename=fdf_string(20,'RigidBonds.File','constraint.aux')
         constraints%n_grp=fdf_integer('RigidBonds.NumberOfGroups', 0                )
         SHAKE_tol        =fdf_double ('RigidBonds.SHAKETol',      1.0E-8_double     )
         RATTLE_tol       =fdf_double ('RigidBonds.RATTLETol',     1.0E-10_double    )
         maxiterSHAKE     =fdf_integer('RigidBonds.MaxIterSHAKE',  500               )
         maxiterRATTLE    =fdf_integer('RigidBonds.MaxIterRATTLE', maxiterSHAKE      )
         const_range      =fdf_double ('RigidBonds.SearchRange',   9.448629943_double)
         allocate (n_bond(constraints%n_grp))
         ! Read block
         if (constraints%n_grp.LE.0) call cq_abort('Error: The number of constraints must be &
                                                   &positive integer:', constraints%n_grp)
         allocate (constraints%grp_name(constraints%n_grp)      , &
                   constraints%n_atom_in_grp(constraints%n_grp) , &
                   constraints%n_subgrp(constraints%n_grp))     
         do i = 1, constraints%n_grp
           if (fdf_block('RigidBondsGroups')) then
             read (unit=input_array(block_start+i-1),fmt=*) &
               j,constraints%grp_name(i),constraints%n_atom_in_grp(i),constraints%n_subgrp(i)
             if (constraints%n_atom_in_grp(i).LE.1) &
               call cq_abort('Error: n_atom_in_grp must be greater than one:', &
                             constraints%n_atom_in_grp(i))
           else
             call cq_abort('RigidBonds block not defined: '//'RigidBondsGroups')
           endif
           call fdf_endblock
         enddo
       endif ! Constraints

       md_ensemble        = fdf_string(3, 'MD.Ensemble', 'nve')
       md_calc_xlmass     = fdf_boolean('MD.CalculateXLMass', .true.)

       ! Thermostat
       md_thermo_type     = fdf_string(20, 'MD.Thermostat', 'none')
       if (leqi(md_thermo_type, 'berendsen')) then
         md_tau_T           = fdf_double('MD.tauT', one)
       else
         md_tau_T           = fdf_double('MD.tauT', one)
       end if
       md_tau_T_equil     = fdf_double('MD.tauTEquil', one)
       md_n_nhc           = fdf_integer('MD.nNHC', 5) 
       md_n_ys            = fdf_integer('MD.nYoshida', 1)
       md_n_mts           = fdf_integer('MD.nMTS', 1)
       flag_thermoDebug   = fdf_boolean('MD.ThermoDebug',.false.)
       md_t_drag          = fdf_double('MD.TDrag', zero)
       allocate(md_nhc_mass(md_n_nhc)) 
       allocate(md_nhc_cell_mass(md_n_nhc)) 
       md_nhc_mass = one
       md_nhc_cell_mass = one
       if (fdf_block('MD.NHCMass')) then
         read(unit=input_array(block_start), fmt=*) md_nhc_mass
       end if
       call fdf_endblock
       if (fdf_block('MD.CellNHCMass')) then
         read(unit=input_array(block_start), fmt=*) md_nhc_cell_mass
       end if
       call fdf_endblock

       ! Barostat
       md_baro_type       = fdf_string(20, 'MD.Barostat', 'none')
       md_target_press    = fdf_double('MD.TargetPressure', zero)
       md_box_mass        = fdf_double('MD.BoxMass', one)
       if (leqi(md_baro_type, 'berendsen')) then
         md_tau_P           = fdf_double('MD.tauP', 50.0_double)
       else
         md_tau_P           = fdf_double('MD.tauP', 100.0_double)
       end if
       md_tau_P_equil     = fdf_double('MD.tauPEquil', 100.0_double)
       md_bulkmod_est     = fdf_double('MD.BulkModulusEst', 100.0_double)
       md_cell_nhc        = fdf_boolean('MD.CellNHC', .true.)
       flag_baroDebug     = fdf_boolean('MD.BaroDebug',.false.)
       md_berendsen_equil = fdf_integer('MD.BerendsenEquil', 0)
       md_tdep            = fdf_boolean('MD.TDEP', .false.)
       md_p_drag          = fdf_double('MD.PDrag', zero)

    else
       call cq_abort("Old-style CQ input no longer supported: please convert")
    end if ! new_format

!**** TM 2017.Nov.3rd
    call check_compatibility
!****lat<$
    call stop_backtrace(t=backtrace_timer,who='read_input')
!****lat>$

    call my_barrier()

    return
   contains
    ! ------------------------------------------------------------------------------
    ! subroutine check_compatibility
    ! ------------------------------------------------------------------------------
    subroutine check_compatibility 
     ! this subroutine checks the compatibility between keywords defined in read_input
     ! add the
     !       2017.11(Nov).03   Tsuyoshi Miyazaki
     implicit none
      !check AtomMove.ExtendedLagrangian(flag_XLBOMD) &  AtomMove.TypeOfRun (runtype)
       if(runtype .NE. 'md') flag_XLBOMD=.false.
     
     return
    end subroutine check_compatibility 
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
!!    - Added header and changed atomicrad to atomicnum
!!   2011/11/16 15:51 dave
!!    - Changes for new blip data
!!   2014/10/12 10:58 lat
!!    - Added species_file
!!   2015/06/08 lat
!!    - Added experimental backtrace
!!   2016/09/16 17:00 nakata
!!    - Added RadiusAtomf, RadiusMS and RadiusLD
!!   2017/04/05 18:00 nakata
!!    - Added charge_up and charge_dn
!!  SOURCE
!!
  subroutine allocate_species_vars

    use dimens,         only: RadiusSupport, RadiusAtomf, RadiusMS, RadiusLD, &
                              NonLocalFactor, InvSRange, atomicnum
    use memory_module,  only: reg_alloc_mem, type_dbl
    use species_module, only: n_species, nsf_species, nlpf_species, npao_species, natomf_species, &
                              charge, charge_up, charge_dn
    use species_module, only: mass, non_local_species, ps_file, ch_file, phi_file 
    use species_module, only: species_label, species_file, type_species
    use global_module,  only: area_general
    use GenComms,       only: cq_abort
    use blip,           only: blip_info

    implicit none

    ! Local variables
    type(cq_timer) :: backtrace_timer
    integer        :: stat

!****lat<$
    call start_backtrace(t=backtrace_timer,who='allocate_species_vars',where=1,level=3)
!****lat>$

    call start_timer(tmr_std_allocation)
    !
    allocate(RadiusSupport(n_species),atomicnum(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating RadiusSupport, atomicnum in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,2*n_species,type_dbl)
    allocate(RadiusAtomf(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating RadiusAtomf in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(RadiusMS(n_species),RadiusLD(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating RadiusMS, RadiusLD in allocate_species_vars: ",n_species,stat)
    call reg_alloc_mem(area_general,2*n_species,type_dbl)
    allocate(blip_info(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating blip_info in allocate_species_vars: ",               n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(NonLocalFactor(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating NonLocalFactor in allocate_species_vars: ",          n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(nsf_species(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating nsf_species in allocate_species_vars: ",             n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(nlpf_species(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating nlpf_species in allocate_species_vars: ",            n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(npao_species(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating npao_species in allocate_species_vars: ",            n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(natomf_species(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating natomf_species in allocate_species_vars: ",          n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(charge(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating charge in allocate_species_vars: ",                  n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(charge_up(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating charge_up in allocate_species_vars: ",               n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(charge_dn(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating charge_dn in allocate_species_vars: ",               n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(mass(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating mass in allocate_species_vars: ",                    n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(non_local_species(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating non_local_species in allocate_species_vars: ",       n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(ps_file(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating ps_file in allocate_species_vars: ",                 n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(ch_file(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating ch_file in allocate_species_vars: ",                 n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(phi_file(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating phi_file in allocate_species_vars: ",                n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(species_label(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating species_label in allocate_species_vars: ",           n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(type_species(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating type_species in allocate_species_vars: ",            n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(InvSRange(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating InvSRange in allocate_species_vars: ",               n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    allocate(species_file(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating species_file in allocate_species_vars: ",            n_species,stat)
    call reg_alloc_mem(area_general,n_species,type_dbl)
    !
    call stop_timer(tmr_std_allocation)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='allocate_species_vars')
!****lat>$

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
  !!    Added information on the type of smearing used if using
  !!    Diagonalsation method
  !!   2011/12/06 17:02 dave
  !!    Small bug fix on formatting numbers
  !!   2011/12/11 L.Tong
  !!    Removed obsolete parameter number_of_bands
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!  SOURCE
  !!
  subroutine write_info(titles, mu, vary_mu, find_chdens, read_phi, &
                        HNL_fac, NODES)

    use datatypes
    use units
    use dimens,               only: r_super_x, r_super_y, r_super_z,   &
                                    n_grid_x, n_grid_y, n_grid_z, r_h, &
                                    r_c
    use block_module,         only: in_block_x, in_block_y, in_block_z
    use species_module,       only: n_species, species_label, mass,    &
                                    charge, ps_file, ch_file,          &
                                    phi_file, species, nsf_species
    use pseudopotential_data, only: core_radius, non_local_species
    use DiagModule,           only: flag_smear_type,           &
                                    iMethfessel_Paxton
    use blip,                 only: blip_info
    use global_module,        only: flag_basis_set, PAOs,blips,        &
                                    flag_precondition_blips, io_lun,   &
                                    flag_Multisite, flag_diagonalisation, flag_neutral_atom
    use minimise,             only: energy_tolerance, L_tolerance,     &
                                    sc_tolerance,                      &
                                    n_support_iterations,              &
                                    n_L_iterations
    use datestamp,            only: datestr, commentver
    use pseudopotential_common, only: flag_neutral_atom_projector, maxL_neutral_atom_projector, &
         numN_neutral_atom_projector

    implicit none

    ! Passed variables
    logical :: vary_mu, find_chdens
    logical :: read_phi
    character(len=80) :: titles
    integer :: NODES 
    integer :: n_title
    real(double) :: mu, alpha, beta, expected_reduction, HNL_fac

    ! Local variables
    integer :: n
    character(len=10) :: today, the_time

    call date_and_time(today, the_time)
    write(io_lun,3) today(1:4), today(5:6), today(7:8), the_time(1:2),&
                    the_time(3:4)
    write(io_lun,&
          '(/10x,"Code compiled on: ",a,/10x,"Version comment: ",/10x,a)') &
         datestr, commentver
    
    write(io_lun,1)
    write(io_lun,2) titles

    if(flag_diagonalisation) then
       write(io_lun,30) 'diagonalisation '
       select case (flag_smear_type)
       case (0)
          write(io_lun,'(/,10x,"Using Fermi-Dirac smearing")')
       case (1)
          write(io_lun,&
                '(/,10x,"Using order ",i2," Methfessel-Paxton smearing")') &
               iMethfessel_Paxton
       end select
    else
       write(io_lun,30) 'order N with LNV'   
    end if
    write(io_lun,4) dist_conv*r_super_x, dist_conv*r_super_y, &
                    dist_conv*r_super_z,d_units(dist_units)

    write(io_lun,9) n_grid_x, n_grid_y, n_grid_z

    write(io_lun,17) in_block_x, in_block_y, in_block_z

    write(io_lun,15) dist_conv*(r_super_x/n_grid_x), d_units(dist_units), &
                     dist_conv*(r_super_y/n_grid_y), d_units(dist_units), &
                     dist_conv*(r_super_z/n_grid_z), d_units(dist_units)

    write(io_lun,18) n_species, d_units(dist_units)

    do n=1, n_species
       write(io_lun,19) n, species_label(n), mass(n), charge(n), &
                        dist_conv*core_radius(n), nsf_species(n)
       if(flag_basis_set==blips) then 
          if(flag_precondition_blips) then
             write(io_lun,'(/13x,"Blip basis with preconditioning")') 
          else
             write(io_lun,'(/13x,"Blip basis - no preconditioning")') 
          end if
          write(io_lun,14) &
               dist_conv*blip_info(n)%SupportGridSpacing, d_units(dist_units), &
               dist_conv*blip_info(n)%BlipWidth,          d_units(dist_units)
       else
          write(io_lun,'(13x,"PAO basis")') 
       end if
    end do

    write(io_lun,20)

    if (flag_Multisite) write(io_lun,'(/10x,"PAOs are contracted to multi-site support functions")')

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

    if(flag_neutral_atom) then
       write(io_lun,fmt='(/13x,"Using neutral atom potential (NAP) formalism")')
       if(flag_neutral_atom_projector) then
          write(io_lun,fmt='(/13x,"Calculating 1- and 2-centre NAP integrals analytically")')
          write(io_lun,fmt='(13x,"Expanding 3-centre NAP integrals with ",i2," projectors up to l=",i2)') &
               maxval(numN_neutral_atom_projector),maxL_neutral_atom_projector
       end if
    end if

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

    ! write (io_lun, 6) number_of_bands

    if (.not.vary_mu) then
       write(io_lun,*) '          mu is constant'
       write(io_lun,16) mu
    endif

    write(io_lun,7) NODES

    write(io_lun,11) n_support_iterations, n_L_iterations

    write(io_lun,13) dist_conv*r_h, d_units(dist_units), &
                     dist_conv*r_c, d_units(dist_units)
    do n=1, n_species
       write(io_lun,131) n, dist_conv*r_h+core_radius(n) * HNL_fac, &
                         d_units(dist_units)
    end do

1   format(/10x,'Job title: ')
2   format(10x,a80)
3   format(/10x,'This job was run on ',a4,'/',a2,'/',a2,' at ',a2,':',a2,/)
4   format(/10x,'The simulation box has the following dimensions',/, &
           10x,'a = ',f11.5,' b = ',f11.5,' c = ',f11.5,' 'a2)
5   format(/10x,'The simulation box contains ',i7,' atoms.')
! 6   format(/10x,'The number of bands is ',f8.2)
7   format(/10x,'The calculation will be performed on ',i5,' processors')
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
15  format(/10x,'integration grid spacing along x ',f9.5,' ',a2,/, &
           10x,'integration grid spacing along y ',f9.5,' ',a2,/, &
           10x,'integration grid spacing along z ',f9.5,' ',a2)
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
29  format(/,10x,'Energy tolerance required:             ',f12.8, &
           /,10x,'L-matrix convergence tolerance:        ',f12.8, &
           /,10x,'Self consistent convergence tolerance: ',f12.8)
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
  !!   readDiagInfo(Number of Kpts, proc grid rows, proc grid cols,
  !!   row block size, col block size)
  !!  PURPOSE
  !!   Reads in data to do with diagonalisation, allocates memory for
  !!   data and broadcasts to all processors
  !!  INPUTS
  !!   integer, intent(OUT) :: proc_rows, proc_cols - rows and columns in
  !!                                                  processor grid
  !!   integer, intent(OUT) :: block_size_r, block_size_c - block sizes for
  !!                                                        rows and columns
  !!                                                        of matrix
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
  !!   2015/06/08 lat
  !!    Added experimental backtrace
  !!   2016/05/09 dave
  !!    Added code to specify lines in reciprocal space
  !!   2016/05/10 dave
  !!    Bug fix: hadn't scaled fractional k-point coordinates to reciprocal
  !!   2017/07/11 dave
  !!    Bug fix: added check for KPointGroups being smaller than number of processors
  !!  SOURCE
  !!
  subroutine readDiagInfo

    use datatypes
    use global_module,   only: iprint_init, rcellx, rcelly, rcellz,  &
                               area_general, ni_in_cell, numprocs,   &
                               species_glob, io_lun
    use numbers,         only: zero, one, two, pi, RD_ERR
    use GenComms,        only: cq_abort, gcopy, myid
    use input_module
    use ScalapackFormat, only: proc_rows, proc_cols, block_size_r,   &
                               block_size_c, proc_groups, matrix_size
    use DiagModule,      only: nkp, kk, wtk, kT, maxefermi,          &
                               flag_smear_type, iMethfessel_Paxton,  &
                               max_brkt_iterations, gaussian_height, &
                               finess, NElec_less
    use energy,          only: SmearingType, MPOrder
    use memory_module,   only: reg_alloc_mem, reg_dealloc_mem,       &
                               type_dbl
    use species_module,  only: nsf_species

    implicit none

    ! Local variables
    type(cq_timer) :: backtrace_timer
    integer        :: stat, iunit, i, j, k, nk_st, nkp_lines
    real(double)   :: a, sum, dkx, dky, dkz
    integer        :: proc_per_group

    ! k-point mesh type
    logical        :: mp_mesh, done, flag_lines_kpoints
    integer,      dimension(1:3)              :: mp
    real(double), dimension(1:3)              :: mp_shift
    real(double), allocatable, dimension(:,:) :: kk_tmp
    real(double), allocatable, dimension(:)   :: wtk_tmp
    integer :: nkp_tmp
    integer :: counter
    
!****lat<$    
    call start_backtrace(t=backtrace_timer,who='readDiagInfo',where=1,level=2)
!****lat>$

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
       if(proc_groups>numprocs) then
          write(io_lun,fmt='(/2x,"*************")')
          write(io_lun,fmt='(2x,"* WARNING ! *")')
          write(io_lun,fmt='(2x,"*************")')
          write(io_lun,fmt='(/2x,"Error setting Diag.KProcGroups.  We have ",i6, &
               &" processes and ",i6," KProcGroups.   Setting to 1 and continuing."/)') numprocs,proc_groups
          proc_groups = 1
       end if
       
       if (iprint_init > 0) write (io_lun, 11) proc_groups 
       ! Read/choose ScaLAPACK processor grid dimensions
       if(fdf_defined('Diag.ProcRows')) then
          proc_rows = fdf_integer('Diag.ProcRows',0)
          proc_cols = fdf_integer('Diag.ProcCols',0)
          if(proc_rows*proc_cols==0) &
               call cq_abort('Diag: proc grid product is zero: ',&
                             proc_rows,proc_cols)
          ! 2011/09/02 L.Tong: now checks correctly if we have more
          ! than one k-point groups
          proc_per_group = numprocs / proc_groups
          if(proc_rows*proc_cols > proc_per_group) &
               call cq_abort('Diag: proc grid product is too large: ', &
                             proc_rows*proc_cols, proc_per_group) 
       else
          ! default values
          if(numprocs<9) then ! Recommended by ScaLapack User Guide
             proc_rows = 1
             proc_cols = numprocs / proc_groups
          else
             done = .false.
             ! 2011/09/02 L.Tong, now works properly for more than one
             ! k-point groups
             proc_per_group = numprocs / proc_groups
             a = sqrt(real(proc_per_group, double))
             proc_rows = aint(a)+1
             do while(.NOT.done.AND.proc_rows>1) 
                proc_rows = proc_rows - 1
                proc_cols = aint(real(proc_per_group,double) / &
                                 real(proc_rows,double))
                if(proc_per_group - proc_rows*proc_cols == 0) &
                     done = .true.
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
               call cq_abort('block_size_r not a factor of matrix size ! ',&
                             matrix_size, block_size_r)
          a = real(matrix_size)/real(block_size_c)
          if(a - real(floor(a))>1e-8_double) &
               call cq_abort('block_size_c not a factor of matrix size ! ',&
                             matrix_size, block_size_c)
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
          ! Do we want lines in reciprocal space ?
          flag_lines_kpoints = fdf_boolean('Diag.KspaceLines',.false.)
          if(flag_lines_kpoints) then
             nkp_lines = fdf_integer('Diag.NumKptLines',1)
             if(iprint_init>1) write(io_lun,fmt='(8x,"Number of Kpoint lines: ",i4)') nkp_lines
             if(nkp_lines<1) call cq_abort("Need to specify how many kpoint lines !",nkp_lines)
             nkp = fdf_integer('Diag.NumKpts',2)
             if(iprint_init>1) write(io_lun,fmt='(8x,"Number of Kpoints in a line: ",i4)') nkp
             allocate(kk(3,nkp*nkp_lines),wtk(nkp*nkp_lines),STAT=stat)
             if(stat/=0) call cq_abort('FindEvals: couldnt allocate kpoints',nkp*nkp_lines)
             call reg_alloc_mem(area_general,4*nkp*nkp_lines,type_dbl)
             kk = zero
             wtk = one/real(nkp*nkp_lines,double)
             sum = zero
             nk_st = 1
             if(fdf_block('Diag.KpointLines'))then
                if(1+block_end-block_start<2*nkp_lines) &
                     call cq_abort("Kpoint line error: ",1+block_end-block_start,nkp_lines)
                do i=1,2*nkp_lines,2
                   read (unit=input_array(block_start+i-1),fmt=*) &
                        kk(1,nk_st),kk(2,nk_st),kk(3,nk_st)
                   read (unit=input_array(block_start+i),fmt=*) &
                        kk(1,nk_st+nkp-1),kk(2,nk_st+nkp-1),kk(3,nk_st+nkp-1)
                   dkx = (kk(1,nk_st+nkp-1) - kk(1,nk_st))/real(nkp-1,double)
                   dky = (kk(2,nk_st+nkp-1) - kk(2,nk_st))/real(nkp-1,double)
                   dkz = (kk(3,nk_st+nkp-1) - kk(3,nk_st))/real(nkp-1,double)
                   if(iprint_init>1) write(io_lun,fmt='(2x,"K-point spacing along line : ",i3,3f7.3)') i,dkx,dky,dkz
                   do j=1,nkp-2
                      kk(1,nk_st+j) = kk(1,nk_st+j-1)+dkx
                      kk(2,nk_st+j) = kk(2,nk_st+j-1)+dky
                      kk(3,nk_st+j) = kk(3,nk_st+j-1)+dkz
                   end do
                   nk_st = nk_st + nkp
                end do
             else
                call cq_abort("Must specify a block Diag.KpointLines to have lines of kpoints !")
             end if
             nkp = nkp*nkp_lines
             ! Write out fractional k-points
             if(iprint_init>0) then
                write(io_lun,7) nkp
                do i=1,nkp
                   write(io_lun,fmt='(8x,i5,3f15.6,f12.3)')&
                        i,kk(1,i),kk(2,i),kk(3,i),wtk(i)
                end do
             end if
             ! Scale from fractional to reciprocal
             do i = 1, nkp
                kk(1,i) = two * pi * kk(1,i) / rcellx
                kk(2,i) = two * pi * kk(2,i) / rcelly
                kk(3,i) = two * pi * kk(3,i) / rcellz
             end do
          else
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
                if(1+block_end-block_start<nkp) &
                     call cq_abort("Kpoint error: ",1+block_end-block_start,nkp)
                do i=1,nkp
                   read (unit=input_array(block_start+i-1),fmt=*) &
                        kk(1,i),kk(2,i),kk(3,i),wtk(i)
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
          if(iprint_init>0) &
               write (io_lun,fmt='(8x,a, 3i3)') &
                     ' Monkhorst-Pack mesh: ', (mp(i), i=1,3)
          if (mp(1) <= 0) &
               call cq_abort('K-points: number of k-points must be > 0!')
          if (mp(2) <= 0) &
               call cq_abort('K-points: number of k-points must be > 0!')
          if (mp(3) <= 0) &
               call cq_abort('K-points: number of k-points must be > 0!')
          nkp_tmp = mp(1)*mp(2)*mp(3)
          if(iprint_init>0) &
               write(io_lun,fmt='(8x,a, i4)') ' Number of k-points: ',nkp_tmp
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
          if(iprint_init>0) &
               write (io_lun,fmt='(8x,a, 3f11.6)') &
                     ' Monkhorst-Pack mesh shift:  ', &
                     (mp_shift(i), i=1,3)
          ! Allocate
          allocate(kk_tmp(3,nkp_tmp),wtk_tmp(nkp_tmp),STAT=stat)
          if(stat/=0) &
               call cq_abort('FindEvals: couldnt allocate kpoints',nkp_tmp)
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
                   kk_tmp(1,counter) = &
                        ((two * (i+1) - mp(1) - 1) / (2 * mp(1))) + mp_shift(1)
                   kk_tmp(2,counter) = &
                        ((two * (j+1) - mp(2) - 1) / (2 * mp(2))) + mp_shift(2)
                   kk_tmp(3,counter) = &
                        ((two * k - mp(3) - 1) / (2 * mp(3))) + mp_shift(3)
                end do
             end do
          end do
          ! Write out fractional k-points
          if(iprint_init>0) then
             write(io_lun,7) nkp_tmp
             do i=1,nkp_tmp
                write(io_lun,fmt='(8x,i5,3f15.6,f12.3)')&
                      i,kk_tmp(1,i),kk_tmp(2,i),kk_tmp(3,i),wtk_tmp(i)
             end do
          end if
          ! Check if k = -k and weed out equivalent k-points. Adjust
          ! wtk accordingly.
          do i = nkp_tmp, 1, -1 ! to get most coords > 0
             if (wtk_tmp(i) > zero) then
                do j = nkp_tmp, 1, -1
                   if (i /= j) then
                      if ((kk_tmp(1,i) == -kk_tmp(1,j)) .and. &
                          (kk_tmp(2,i) == -kk_tmp(2,j)) .and. &
                          (kk_tmp(3,i) == -kk_tmp(3,j))) then
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
             if (wtk_tmp(i) > RD_ERR ) then
                nkp = nkp + 1 
                sum = sum + wtk_tmp(i)
             end if
          end do
          ! Fill in the real k-point array
          allocate(kk(3,nkp),wtk(nkp),STAT=stat)
          if(stat/=0) &
               call cq_abort('FindEvals: couldnt allocate kpoints',nkp)
          call reg_alloc_mem(area_general,4*nkp,type_dbl)
          counter = 0
          do i = 1, nkp_tmp
             if (wtk_tmp(i) > RD_ERR ) then
                counter = counter + 1
                kk(1,counter) = kk_tmp(1,i)
                kk(2,counter) = kk_tmp(2,i)
                kk(3,counter) = kk_tmp(3,i)
                !DRB fix to sum weights to one
                wtk(counter) = wtk_tmp(i)/sum
             end if
          end do
          if(counter>nkp) &
               call cq_abort("Error in kpt symmetry ! ",counter,nkp)
          call reg_dealloc_mem(area_general,4*nkp_tmp,type_dbl)
          deallocate(kk_tmp,wtk_tmp,STAT=stat)
          if(stat/=0) &
               call cq_abort('FindEvals: couldnt deallocate kpoints',&
                             nkp_tmp)
          if(iprint_init>0) then
             write(io_lun,*)
             write(io_lun,10) nkp
             do i=1,nkp
                write (io_lun,fmt='(8x,i5,3f15.6,f12.3)') &
                      i,kk(1,i),kk(2,i),kk(3,i),wtk(i)
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
             write (io_lun,fmt='(8x,i5,3f15.6,f12.3)') &
                   i,kk(1,i),kk(2,i),kk(3,i),wtk(i)
          end do
          write(io_lun,fmt='(/8x,"Finished reading Kpoints"/)')
       end if
       
       ! Write out smearing temperature
       if(iprint_init>0) &
            write (io_lun,'(10x,"Temperature used for smearing: ",f10.6)') kT
    end if ! myid==0

    ! Distribute data to all processors
    call gcopy(block_size_r)
    call gcopy(block_size_c)
    call gcopy(matrix_size)
    call gcopy(proc_groups)
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

!****lat<$    
    call stop_backtrace(t=backtrace_timer,who='readDiagInfo')
!****lat>$

    return

2   format(8x,'Block size (row, col): ',2i5)
3   format(8x,'Proc grid (row, col): ',2i5)
4   format(/8x,'***WARNING***',/,2x,&
           'No Kpoints block found - defaulting to Gamma point')
5   format(8x,i4,' Kpoints in Cartesian (inverse Angstrom) form: ')
51  format(8x,i4,' symmetry inequivalent Kpoints in Cartesian ')
52  format(12x,' (inverse Angstrom) form: ')
6   format(8x,'Kpt: ',i4,' : ',3f12.8,' Weight: ',f10.6)
7   format(8x,' All ',i4,' Kpoints in fractional coordinates: ')
8   format(/8x,'***WARNING***',/,8x,&
           'No Kpoint mesh found - defaulting to Gamma point')
9   format(/8x,'***WARNING***',/,8x,&
           'Specified Kpoint mesh shift in fractional coords >= 1.0.')
10  format(8x,i4,' symmetry inequivalent Kpoints in fractional coordinates: ')
11  format(8x, 'Number of processor groups for k-point parallelisation: ', i5)

  end subroutine readDiagInfo
  !!***

!%%!   subroutine readDiagInfo
!%%! 
!%%!     use datatypes
!%%!     use global_module, only: iprint_init, rcellx, rcelly, rcellz, area_general
!%%!     use numbers, only: two, pi
!%%!     use GenComms, only: cq_abort, gcopy, myid
!%%!     use fdf, only: fdf_integer, fdf_block
!%%!     use ScalapackFormat, only: proc_rows, proc_cols, block_size_r, block_size_c
!%%!     use DiagModule, only: nkp, kk, wtk
!%%!     use memory_module, only: reg_alloc_mem, type_dbl
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

  ! ------------------------------------------------------------------------------
  ! Subroutine read_input_aux
  ! ------------------------------------------------------------------------------
  
  !!****f* initial_read/read_input_aux *
  !!  NAME 
  !!   read_input_aux
  !!  USAGE
  !!   call read_input_aux(aux)
  !!  PURPOSE
  !!   Reads auxiliary input files
  !!  INPUTS
  !!   type(group_aux) : aux
  !!  AUTHOR
  !!   M.Arita
  !!  CREATION DATE
  !!   2014/02/04
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !!
  subroutine read_input_aux(aux)
    ! Module usage
    use global_module, ONLY: io_lun
    use auxiliary_types, ONLY: group_aux
    use GenComms, ONLY: myid,cq_abort,gcopy
    use io_module, ONLY: io_assign,io_close
    use input_module, ONLY: fdf_block,fdf_endblock,input_array,block_start

    implicit none
    ! passed variables
    type(group_aux) :: aux
    ! local varables
    integer :: i,j,jj,l,natom,ibeg,stat,lun_aux
    integer :: isize,jsize
    logical :: done
    character(20) :: filename_aux
    character(100) :: line,line_blck

    ! Allocation
    isize = 0
    jsize = 0
    do i = 1, aux%n_grp
      isize = isize + aux%n_atom_in_grp(i)*aux%n_subgrp(i)
      jsize = jsize + aux%n_subgrp(i)
    enddo
    allocate (aux%ibeg_grp(aux%n_grp), STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating ibeg_grp: ',aux%n_grp)
    allocate (aux%glob_atom(isize), STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating n_grp: ',isize)
    allocate (aux%iatom_beg(jsize), STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating iatom_beg: ',jsize)

    if (myid.EQ.0) then
      ! Open file
      filename_aux=aux%filename
      call io_assign(lun_aux)
      open (lun_aux,file=filename_aux,status='old',action='read', iostat=stat)
      if (stat.NE.0) call cq_abort('Error opening auxiliary file !')
      ! Reckon addresses
      aux%ibeg_grp(1) = 1
      do i = 1, aux%n_grp-1
        aux%ibeg_grp(i+1) = aux%ibeg_grp(i) + aux%n_atom_in_grp(i)*aux%n_subgrp(i)
      enddo
      aux%iatom_beg(1) = 1
      l = 0
      natom = aux%n_atom_in_grp(1)
      do i = 1, aux%n_grp
        do j = 1, aux%n_subgrp(i)
          l = l + 1
          if (l.NE.1) then
            aux%iatom_beg(l) = aux%iatom_beg(l-1) + natom
          else
            cycle
          endif
          natom = aux%n_atom_in_grp(i)
        enddo
      enddo
      ! Read aux%filename
      do i = 1, aux%n_grp
        done = .false.
        ibeg = aux%ibeg_grp(i)
        do while (.NOT.done)
          read (lun_aux,'(a100)') line
          line_blck = trim(adjustl(line(7:100)))
          if (trim(adjustl(line_blck)).EQ.trim(aux%grp_name(i))) then
            do j = 1, aux%n_subgrp(i)
              read (lun_aux,*) jj,aux%glob_atom(ibeg:ibeg+aux%n_atom_in_grp(i)-1)
              ibeg = ibeg + aux%n_atom_in_grp(i)
            enddo
            done = .true.
          endif
        enddo !done
        rewind (lun_aux)
      enddo

      !% NOTE: I'd like to make use of fdf_block but it doesn't work...
      !% ibeg = 1
      !% do i = 1, aux%n_grp
      !%   !if (fdf_block(trim(adjustl(aux%grp_name(i))))) then
      !%   if (fdf_block(aux%grp_name(i))) then
      !%     do j = 1, aux%n_subgrp(i)
      !%       read (unit=input_array(block_start+j-1),fmt=*) &
      !%         jj,aux%glob_atom(ibeg:ibeg+aux%n_atom_in_grp(i)-1)
      !%         ibeg = ibeg + aux%n_atom_in_grp(i)
      !%     enddo
      !%   else
      !%     write (io_lun,*) trim(adjustl(aux%grp_name(i)))," not found !"
      !%     call cq_abort('Block not defined in auxiliary file:')
      !%   endif
      !%   rewind (lun_aux)
      !%   call fdf_endblock
      !% enddo

      ! Close file
      call io_close(lun_aux)
    endif ! myid

    ! Broadcast data to all processors
    call gcopy(aux%ibeg_grp,aux%n_grp)
    call gcopy(aux%glob_atom,isize)
    call gcopy(aux%iatom_beg,jsize)

    if (myid.EQ.0) write (io_lun,*) "Completed read_input_aux"

    return
  end subroutine read_input_aux

end module initial_read
