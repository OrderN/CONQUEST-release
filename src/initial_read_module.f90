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

  ! Index number for loading DM etc
  integer, save :: index_MatrixFile

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
  !!   2018/09/06 16:45 nakata
  !!    Changed to allow to use MSSFs larger than single-zeta size (i.e., not minimal)
  !!   2019/06/06 17:30 jack poulton and nakata
  !!    Changed to allow SZP-MSSFs with flag_MSSF_nonminimal when checking the number of MSSFs
  !!    Added MSSF_nonminimal_species
  !!   2019/07/04 zamaan
  !!    Added bibliography
  !!   2019/11/18 tsuyoshi
  !!    Removed flag_MDold
  !!  SOURCE
  !!
  subroutine read_and_write(start, start_L, inode, ionode,          &
       vary_mu, mu, find_chdens)

    use datatypes
    use numbers
    use io_module,              only: read_mult,                       &
         read_atomic_positions,           &
         pdb_template, pdb_format,        &
         print_process_info, titles
    use group_module,           only: parts, part_method
    use construct_module,       only: init_group, init_primary
    use maxima_module,          only: maxpartsproc, maxatomsproc
    use global_module,          only: id_glob,x_atom_cell,y_atom_cell, &
         z_atom_cell, numprocs,           &
         iprint_init, nspin,              &
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
    use dimens,                 only: set_dimensions, find_grid
    use species_module,         only: n_species, species, charge,      &
         non_local_species,               &
         nsf_species, npao_species,       &
         natomf_species, charge_up, charge_dn
    use density_module, only: flag_InitialAtomicSpin
    use GenComms,               only: my_barrier, cq_abort, cq_warn
    use pseudopotential_data,   only: non_local, read_pseudopotential
    use pseudopotential_common, only: core_radius, pseudo_type, OLDPS, &
         SIESTA, ABINIT
    use pseudo_tm_module,       only: init_pseudo_tm
    use pseudo_tm_info,         only: pseudo
    use input_module,           only: fdf_string, fdf_out, io_close, leqi
    use force_module,           only: tot_force
    use constraint_module,      only: flag_RigidBonds,constraints
    use support_spec_format,    only: flag_one_to_one, symmetry_breaking, read_option
    use multisiteSF_module,     only: flag_LFD_ReadTVEC, &
         flag_MSSF_nonminimal, &   !nonmin_mssf
         MSSF_nonminimal_species   !nonmin_mssf
    use md_control,             only: md_position_file
    use pao_format
    use XC,                     only: flag_functional_type, flag_different_functional

    implicit none

    ! Passed variables
    logical           :: vary_mu, start, start_L
    logical           :: find_chdens
    integer           :: inode, ionode
    real(double)      :: mu

    ! Local variables
    character(len=80) :: sub_name = "read_and_write"
    type(cq_timer)    :: backtrace_timer
    character(len=80) :: def
    character(len=80) :: atom_coord_file
    character(len=80) :: part_coord_file

    integer      :: i, j, acz, prncpl, prncpl0
    integer      :: ncf_tmp, stat, i_species
    integer      :: count_SZ, count_SZP
    real(double) :: ecore_tmp
    real(double) :: HNL_fac
    real(double), dimension(8) :: occ_n

    ! for checking the sum of electrons of spin channels
    real(double) :: sum_elecN_spin
    real(double) :: charge_tmp

    !****lat<$
    call start_backtrace(t=backtrace_timer,who='read_and_write',where=1,level=1)
    !****lat>$

    ! read input data: parameters for run
    call read_input(start, start_L, titles, vary_mu, mu,&
         find_chdens, HNL_fac)
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
    ! Functional consistency check
    if(flag_functional_type==0) then ! Not set, so take from pseudopotentials
       flag_functional_type = pseudo(1)%functional
       if(inode==ionode) then ! Check and warn
          do i_species = 2, n_species
             if(flag_functional_type/=pseudo(i_species)%functional) &
                  call cq_abort("Functionals differ between pseudopotential files: ",&
                  flag_functional_type, pseudo(i_species)%functional)
          end do
       end if
    else if(inode==ionode) then ! Check and warn
       do i_species = 1, n_species
          if(flag_functional_type/=pseudo(i_species)%functional) then
             if(flag_different_functional) then
                call cq_warn(sub_name, "Functional in input file differs to pseudopotential but proceeding: ",&
                     flag_functional_type, pseudo(i_species)%functional)
             else
                call cq_abort("Functional in input file differs to pseudopotential: ",&
                     flag_functional_type, pseudo(i_species)%functional)
             end if
          end if
       end do
    end if
    !if(iprint_init>4) write(io_lun,fmt='(10x,"Proc: ",i4," done pseudo")') inode

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
       call read_atomic_positions(md_position_file)
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

    !
    !
    !
    ! Set flag_one_to_one and check the number of SFs
    if (nspin.eq.1) then
       nspin_SF = 1
       flag_SpinDependentSF = .false.
    endif
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
          call cq_warn(sub_name, "Cannot set 1:1 flag, read coefficients and have NSF==NPAO; setting 1:1 false")
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
       if (flag_SpinDependentSF) then
          nspin_SF = nspin
          if (inode==ionode) write(io_lun,'(2X,A,I1)') 'Support functions are spin-dependent, nspin_SF = ',nspin_SF
       endif
       if (flag_one_to_one.AND.flag_vary_basis) then
          call cq_warn(sub_name, "minE.VaryBasis T is not compatible with NSF==NPAO.  Setting flag to F.")
          flag_vary_basis = .false.
       end if
       ! Check symmetry-breaking for contracted SFs
       if (.not.flag_one_to_one) then
          do i = 1, n_species
             count_SZ  = 0
             count_SZP = 0
             do j = 0, pao(i)%greatest_angmom
                if(pao(i)%angmom(j)%n_zeta_in_angmom>0) then
                   count_SZP = count_SZP + (2*j+1) ! Accumulate no. of ang. mom. components
                   occ_n(:) = zero
                   do acz = 1, pao(i)%angmom(j)%n_zeta_in_angmom
                      prncpl = pao(i)%angmom(j)%prncpl(acz)
                      if (acz.ne.1 .and. prncpl.ne.prncpl0) count_SZP = count_SZP + (2*j+1)
                      occ_n(prncpl) = occ_n(prncpl) + pao(i)%angmom(j)%occ(acz)
                      prncpl0 = prncpl
                   enddo
                   do prncpl = 1, 8
                      if (occ_n(prncpl).gt.zero) count_SZ = count_SZ + (2*j+1)
                   enddo
                endif
             enddo
             if (flag_Multisite) then
                ! multi-site SFs (MSSFs) are symmetry-breaking usually.
                ! The default of the number of the MSSFs is SingleZeta-size.
                ! SZP-size MSSFs with flag_MSSF_nonminimal is also available.
                ! If the user wants to use the other number of MSSFs,
                ! the user must provide the initial SFcoeffmatrix or the trial vectors for the LFD method.
                if (nsf_species(i).eq.count_SZ) then
                   MSSF_nonminimal_species(i) = 1   ! SZ-size MSSF
                else
                   if (.not.flag_MSSF_nonminimal) then
                      if(inode==ionode) write(io_lun,'(A/A,I3,A,I3/2A)') &
                           "You have a major problem with your multi-site support functions.",&
                           "Number of multi-site SFs", nsf_species(i), &
                           " is not equal to single-zeta (SZ) size", count_SZ, &
                           "You need to set Multisite.nonminimal to be .true. or ",&
                           "set Atom.NumberOfSupports to be SZ size."
                      call cq_abort("Multi-site support function error for species ",i)
                   else if (nsf_species(i).eq.count_SZP) then
                      MSSF_nonminimal_species(i) = 2   ! SZP-size MSSF
                   else if (nsf_species(i).ne.count_SZP) then
                      MSSF_nonminimal_species(i) = 3   ! other size
                      if (.not.flag_LFD_ReadTVEC .and. .not.read_option) then
                         if (inode==ionode) write(io_lun,'(A/A,I3,A,I3,A,I3/A/A/A/A)') &
                              "You have a major problem with your multi-site support functions.",&
                              "Number of multi-site SFs", nsf_species(i), &
                              " is not equal to single-zeta (SZ) size", count_SZ, &
                              " nor single-zeta plus polarisation (SZP) size", count_SZP, &
                              "Since the automatic setting of the trial vectors in the local filter diagonalisation method is ",&
                              "avaialble only for SZ or SZP size,",&
                              "you need to provide initial SFcoeffmatrix or trial vectors for the LFD method, ",&
                              "or change Atom.NumberOfSupports to be SZ or SZP size."
                         call cq_abort("Multi-site support function error for species ",i)
                      endif
                   endif
                endif
             else
                ! If number of support functions is less than total number of ang. mom. components (ignoring
                ! for now multiple zetas) then there is a formal problem with basis set: we require the user
                ! to set an additional flag to assert that this is really desired
                if(count_SZP>nsf_species(i)) then
                   if(.NOT.symmetry_breaking.OR..NOT.read_option) then
                      if(inode==ionode) then
                         write(io_lun,fmt='("You have a major problem with your basis set.")')
                         write(io_lun,fmt='("There are fewer support functions than the minimal angular momentum")')
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
    if (iprint_init > 3 .and. inode==ionode) &
         write(io_lun,fmt='(6x,a,L2)') 'PAO:SF mapping flag is ',flag_one_to_one
    if (iprint_init> 3 .and. inode==ionode) then
       if(atomf==sf) then
          write(io_lun,fmt='(6x,"Primitive atom functions are the support functions"/)')
       else if(atomf==paof) then
          write(io_lun,fmt='(6x,"Primitive atom functions are the pseudo-atomic orbitals"/)')
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
    ! Removed second condition otherwise variables not set 2022/10/31 08:36 dave
    if (nspin == 2) then ! .and. ne_magn_in_cell > zero) then
       ! Yes but should consider preceding set up via Spin.NeUP/Spin.NeDN
       if ( abs(ne_spin_in_cell(1))<RD_ERR .and. abs(ne_spin_in_cell(1))<RD_ERR) then
          ne_spin_in_cell(1) = half * (ne_magn_in_cell + ne_in_cell)
          ne_spin_in_cell(2) = ne_spin_in_cell(1) - ne_magn_in_cell
       end if
       !
    end if
    !
    ! Calculate the number of electrons in the spin channels
    ! if (flag_fix_spin_population) then
    sum_elecN_spin = ne_spin_in_cell(1) + ne_spin_in_cell(2)
    if (abs(sum_elecN_spin - ne_in_cell)>very_small) then
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
    ! IF (flag_InitialAtomicSpin) : atomic spin density will be set
    !
    if(flag_InitialAtomicSpin) then
      do i = 1, n_species
       charge_tmp = charge_up(i) + charge_dn(i)
       if(charge_tmp < RD_ERR) then   !
         charge_up(i) = half*charge(i)
         charge_dn(i) = half*charge(i)
         !if(inode .eq. ionode) write(io_lun,fmt='(6x,a,i3,a,2f15.8)') &
         ! 'ispecies = ', i,' charge_up, dn = ',charge_up(i), charge_dn(i)
       endif
      enddo !i = 1, n_species
     ne_spin_in_cell(:) = zero
     do i = 1, ni_in_cell
       ! number_of_bands = number_of_bands + half * charge(species(i))
       ne_spin_in_cell(1) = ne_spin_in_cell(1) + charge_up(species(i))
       ne_spin_in_cell(2) = ne_spin_in_cell(2) + charge_dn(species(i))
     end do
     !if(inode .eq. ionode) write(io_lun,fmt='(6x,a,2f15.8)') 'ne_spin_in_cell(1:2) = ',ne_spin_in_cell(1),ne_spin_in_cell(2)
    endif
    !
    !
    ! Set up various lengths, volumes, reciprocals etc. for convenient use
    call set_dimensions(inode, ionode,HNL_fac, non_local, n_species, &
         non_local_species, core_radius)

    ! write out some information on the run
    if (inode == ionode) &
         call write_info(titles, mu, vary_mu, HNL_fac, numprocs)

    !****lat<$
    call stop_backtrace(t=backtrace_timer,who='read_and_write')
    !****lat>$

    !**** TM 2020.Jul.30
    call check_compatibility

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
  !!   2017/11/13 18:15 nakata
  !!    Added flag_normalise to normalise each eigenstate of PDOS
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
  !!   2018/07/16 17:09 dave
  !!    Adding zero T requirement for flag_quench_MD as well as FIRE
  !!   2018/09/19 18:30 nakata
  !!    Added flag_pDOS_angmom for orbital angular momentum resolved PDOS
  !!   2018/10/22 14:26 dave & jsb
  !!    Added flag_PDOS_lm for (l,m) projected DOS
  !!   2018/10/30 11:52 dave
  !!    added flag_PDOS_include_semicore to allow inclusion/exclusion of semi-core states from PDOS
  !!   2019/02/28 zamaan
  !!    added stress and enthalpy tolerances for cell optimisation
  !!   2019/03/28 zamaan
  !!    Added flag_stress and flag_full_stress to toggle calculation of stress
  !!    and off-diagonal elements respectively
  !!   2019/05/21 zamaan
  !!    Added flag for RNG seed
  !!   2019/10/08 nakata
  !!    Fixed bug for the case when reading only Kmatrix but not SFcoeff (for MSSFs)
  !!   2019/11/14 20:00 nakata
  !!    set flag_SpinDependentSF to be .true. in default
  !!   2019/11/18 14:36 dave
  !!    Added flag_variable_cell set to true for cell optimisation or NPT dynamics
  !!   2019/12/04 08:09 dave
  !!    Made default pseudopotential type Hamann to fit with Conquest ion file generator
  !!    and made default grid 100Ha; removed check for old input file format; default method diagonalisation
  !!   2019/12/26 tsuyoshi  (2019/12/29)
  !!     General.LoadL => General.LoadLorK  => General.LoadDM
  !!     AtomMove.ReuseL => AtomMove.ReuseLorK => AtomMove.ReuseDM
  !!   2020/01/06 15:43 dave
  !!    Keywords for equilibration
  !!   2020/01/07 tsuyoshi
  !!     Default setting of MakeInitialChargeFromK has been changed
  !!   2022/10/28 15:56 lionel
  !!     Added ASE output file setup ; default is F
  !!   2022/12/14 10:01 dave and tsuyoshi
  !!     Update test for solution method (diagon vs ordern) following issue #47
  !!  TODO
  !!  SOURCE
  !!
  subroutine read_input(start, start_L, titles, vary_mu,&
       mu, find_chdens, HNL_fac)

    use datatypes
    use numbers
    use units
    use global_module, only: iprint, flag_vary_basis, &
         flag_self_consistent,                     &
         flag_precondition_blips,                  &
         flag_fix_spin_population,                 &
         flag_fractional_atomic_coords,            &
         ne_in_cell,          &
         max_L_iterations, flag_read_blocks,       &
         runtype, restart_DM, restart_rho,         &
         flag_basis_set, blips, PAOs,              &
         flag_test_forces, UseGemm,                &
         flag_fractional_atomic_coords,            &
         ne_spin_in_cell, nspin, spin_factor,      &
         ne_magn_in_cell,                          &
         max_L_iterations,    &
         flag_reset_dens_on_atom_move,             &
         flag_continue_on_SC_fail, iprint_init,    &
         iprint_mat, iprint_ops, iprint_DM,        &
         iprint_SC, iprint_minE, iprint_time,      &
         iprint_MD, iprint_index, iprint_gen,      &
         iprint_pseudo, iprint_basis, iprint_exx,  &
         iprint_intgn,iprint_MDdebug, &
         global_maxatomspart, load_balance,        &
         flag_assign_blocks, &
         io_ase, write_ase, ase_file,   &
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
         flag_MDcontinue,flag_MDdebug,             &
         flag_thermoDebug, flag_baroDebug, &
         flag_LmatrixReuse,flag_TmatrixReuse,flag_SkipEarlyDM,McWFreq, &
         restart_T,restart_X,flag_XLBOMD,flag_propagateX,              &
         flag_propagateL,flag_dissipation,integratorXL, flag_FixCOM,   &
         flag_exx, exx_alpha, exx_scf, exx_scf_tol, exx_siter,         &
         flag_out_wf,max_wf,out_wf,wf_self_con, flag_fire_qMD, &
         flag_write_DOS, flag_write_projected_DOS, &
         E_wf_min, E_wf_max, flag_wf_range_Ef, &
         mx_temp_matrices, flag_neutral_atom, flag_diagonalisation, &
         flag_SpinDependentSF, flag_Multisite, flag_LFD, flag_SFcoeffReuse, &
         flag_opt_cell, cell_constraint_flag, flag_variable_cell, &
         cell_en_tol, optcell_method, cell_stress_tol, &
         flag_stress, flag_full_stress, rng_seed, &
         flag_atomic_stress, flag_heat_flux, flag_DumpMatrices, flag_calc_pol, flag_do_pol_calc, &
         i_pol_dir, i_pol_dir_st, i_pol_dir_end
    use dimens, only: GridCutoff,    &
         n_grid_x, n_grid_y, n_grid_z, r_c,         &
         RadiusSupport, RadiusAtomf, RadiusMS, RadiusLD, &
         NonLocalFactor, InvSRange,                      &
         min_blip_sp, flag_buffer_old, AtomMove_buffer,  &
         r_dft_d2, r_exx
    use block_module, only: in_block_x, in_block_y, in_block_z, &
         blocks_raster, blocks_hilbert
    use species_module, only: species_label, charge, mass, n_species,  &
         charge, charge_up, charge_dn,            &
         ps_file,     &
         nsf_species, &
         non_local_species, type_species,         &
         species_file, species_from_files
    use GenComms,   only: gcopy, my_barrier, cq_abort, inode, ionode, cq_warn
    !use DiagModule, only: diagon
    !use DiagModule,             only: flag_pDOS_include_semicore
    use energy,     only: flag_check_Diag
    use DMMin,      only: maxpulayDMM, maxpulaystepDMM, minpulaystepDMM, &
         LinTol_DMM, n_dumpL
    !TM
    use pseudopotential_common, only: pseudo_type, OLDPS, SIESTA, &
         ABINIT, flag_angular_new, &
         flag_neutral_atom_projector, maxL_neutral_atom_projector, numN_neutral_atom_projector
    use SelfCon, only: A, flag_linear_mixing, EndLinearMixing, q0, q1,&
         n_exact, maxitersSC, maxearlySC, maxpulaySC,   &
         atomch_output, flag_Kerker, flag_wdmetric, minitersSC, &
         flag_newresidual, flag_newresid_abs, n_dumpSCF
    use density_module, only: flag_InitialAtomicSpin, flag_DumpChargeDensity
    use S_matrix_module, only: InvSTolerance, InvSMaxSteps,&
         InvSDeltaOmegaTolerance
    use blip,          only: blip_info, init_blip_flag, alpha, beta
    use maxima_module, only: lmax_ps
    use control,       only: MDn_steps, MDfreq, XSFfreq, XYZfreq,   &
         MDcgtol, CGreset, LBFGS_history, sqnm_trust_step
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
         pdb_output, banner, get_file_name, time_max, &
         flag_MatrixFile_RankFromZero, flag_MatrixFile_BinaryFormat, &
         flag_MatrixFile_BinaryFormat_Grab, flag_MatrixFile_BinaryFormat_Dump, &
         flag_MatrixFile_BinaryFormat_Dump_END, atom_output_threshold, flag_coords_xyz

    use group_module,     only: part_method, HILBERT, PYTHON
    use H_matrix_module,  only: flag_write_locps, flag_dump_locps
    use pao_minimisation, only: InitStep_paomin
    use timer_module,     only: time_threshold,lun_tmr, TimingOn, &
         TimersWriteOut, BackTraceOn
    use input_module, only: load_input, input_array, block_start, &
         block_end, fdf_boolean, fdf_integer, fdf_double, fdf_string, &
         fdf_block, fdf_defined, leqi, io_assign, fdf_endblock
    use cdft_data, only: cDFT_Type, cDFT_MaxIterations, cDFT_NAtoms, &
         cDFT_Target, cDFT_Tolerance,                &
         cDFT_NumberAtomGroups, cDFT_AtomList,       &
         cDFT_BlockLabel
    use sfc_partitions_module, only: n_parts_user, average_atomic_diameter, gap_threshold
    use XLBOMD_module,         only: XLInitFreq,kappa
    use mult_module,           only: maxiter_Dissipation, InitAtomicDistance_max, &
         InitAtomicDistance_min, flag_check_init_atomic_coords
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
         flag_LFD_nonSCF, LFD_ThreshE, LFD_ThreshD,             &
         LFD_Thresh_EnergyRise, LFD_max_iteration,              &
         flag_LFD_MD_UseAtomicDensity,  flag_MSSF_nonminimal,   &
         n_dumpSFcoeff, flag_mix_LFD_SCF,                       &
         MSSF_nonminimal_offset ! nonmin_mssf
    use md_control, only: md_tau_T, md_n_nhc, md_n_ys, md_n_mts, md_nhc_mass, &
         target_pressure, md_baro_type, md_tau_P, &
         md_thermo_type, md_bulkmod_est, md_box_mass, &
         flag_write_xsf, md_cell_nhc, md_nhc_cell_mass, &
         md_calc_xlmass, md_equil_steps, md_equil_press, &
         md_tau_T_equil, md_tau_P_equil, md_p_drag, &
         md_t_drag, md_cell_constraint, flag_write_extxyz, MDtimestep, md_ensemble
    use md_model,   only: md_tdep
    use move_atoms,         only: threshold_resetCD, &
         flag_stop_on_empty_bundle, &
         enthalpy_tolerance, cg_line_min, safe, backtrack, adapt_backtrack
    use Integrators, only: fire_alpha0, fire_f_inc, fire_f_dec, fire_f_alpha, fire_N_min, &
         fire_N_max, fire_max_step, fire_N_below_thresh
    use XC, only : flag_functional_type, functional_hartree_fock, functional_hyb_pbe0, &
         flag_different_functional
    use biblio, only: flag_dump_bib
    !2019/12/27 tsuyoshi
    use density_module,  only: method_UpdateChargeDensity,DensityMatrix,AtomicCharge
    use force_module, only: mix_input_output_XC_GGA_stress

    implicit none

    ! Passed variables
    logical           :: vary_mu, find_chdens
    logical           :: start, start_L
    real(double)      :: mu, HNL_fac
    character(len=80) :: titles

    ! Local variables
    character(len=80)  :: sub_name = "read_input"
    !type(block), pointer :: bp      ! Pointer to a block read by fdf
    !type(parsed_line), pointer :: p ! Pointer to a line broken into tokens by fdf
    type(cq_timer)    :: backtrace_timer


    integer           :: i, j, stat
    integer           :: i_LFD, j_LFD, k_LFD, n_LFD, line_LFD
    character(len=20) :: def, tmp2
    character(len=80) :: coordfile, timefile, timefileroot
    character(len=10) :: basis_string, part_mode
    character(len=6)  :: method ! To find whether we diagonalise or use O(N)
    character(len=5)  :: ps_type !To find which pseudo we use
    character(len=8)  :: tmp
    logical :: flag_ghost, find_species, test_ase

    ! spin polarisation
    logical :: flag_spin_polarisation
    real(double) :: sum_elecN_spin
    real(double) :: charge_tmp

    ! Set defaults
    vary_mu  = .true.

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
    !
    if (fdf_boolean('IO.WriteOutToASEFile',.false.)) then
       !call io_assign(io_ase) ! Reserve this unit on all processes
       !write_ase = .false.
       if (inode == ionode) then
          write_ase = .true.
          call io_assign(io_ase)
          open(unit=io_ase,file=ase_file,iostat=stat)
          if (stat /= 0) &
               call cq_abort("Failed to open ASE Conquest output file", stat)
       end if
    else
       if (inode == ionode) write_ase = .false.
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
    flag_dump_bib = fdf_boolean('IO.DumpBibTeX',      .true.)
    flag_dump_L   = fdf_boolean('IO.DumpL',           .true. )
    !flag for dumping locps
    flag_dump_locps(1) = fdf_boolean('IO.DumpHarPot',  .false.)
    flag_dump_locps(2) = fdf_boolean('IO.DumpXCPot',  .false.)
    flag_dump_locps(3) = fdf_boolean('IO.DumpPSPot',  .false.)
    flag_dump_locps(4) = fdf_boolean('IO.DumpESPot',  .false.)
    flag_dump_locps(5) = fdf_boolean('IO.DumpTotPot',  .false.)
    flag_write_locps = any(flag_dump_locps(:))
    atomch_output = fdf_boolean('IO.AtomChargeOutput',.false.)
    !
    ! Seed for RNG
    rng_seed      = fdf_integer('General.RNGSeed', -1)
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
    flag_check_init_atomic_coords = &
         fdf_boolean('IO.CheckInitialAtomicDistances',.true.)
    if (flag_check_init_atomic_coords) then
       InitAtomicDistance_max = fdf_double('IO.InitAtomicDistance_max', 50.0_double)
       InitAtomicDistance_min = fdf_double('IO.InitAtomicDistance_min',  0.5_double)
    end if
    atom_output_threshold = fdf_integer('IO.AtomOutputThreshold',200)
    flag_coords_xyz = fdf_boolean('IO.AtomCoordsXYZ',.false.)
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
    ! Default to 50 Ha cutoff for grid
    GridCutoff = fdf_double('Grid.GridCutoff',50.0_double)
    ! Grid points
    n_grid_x   = fdf_integer('Grid.PointsAlongX',0)
    n_grid_y   = fdf_integer('Grid.PointsAlongY',0)
    n_grid_z   = fdf_integer('Grid.PointsAlongZ',0)
    ! Number of different iterations - not well defined
    n_L_iterations       = fdf_integer('DM.LVariations',50)
    max_L_iterations     = n_L_iterations
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
    r_c     = fdf_double('DM.L_range',       one )
    HNL_fac = zero !fdf_double('non_local_factor', zero)
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
    ! Solution method - O(N) or diagonalisation - diagonalisation default
    method = fdf_string(6,'DM.SolutionMethod','diagon')
    if(leqi(method,'ordern')) then
       flag_diagonalisation = .false.
       flag_check_Diag = .false.
    else
       flag_diagonalisation = .true.
       flag_check_Diag = .true.
    end if
    ! Read basis set
    basis_string = fdf_string(10,'Basis.BasisSet','PAOs')
    if(leqi(basis_string,'blips')) then
       flag_basis_set = blips
    else if(leqi(basis_string,'PAOs')) then
       flag_basis_set = PAOs
    end if
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
    ! Hamann is default to fit with Conquest ion file generator
    ps_type = fdf_string(5,'General.PseudopotentialType','haman')
    ! Write out pseudopotential type
    if(leqi(ps_type,'siest')) then
       !if(inode==ionode.AND.iprint_init>0) &
       !     write(io_lun,fmt='(10x,"SIESTA pseudopotential will be used. ")')
       pseudo_type = SIESTA
    else if(leqi(ps_type,'plato').OR.leqi(ps_type,'haman')) then
       !if(inode==ionode.AND.iprint_init>0) &
       !     write(io_lun,fmt='(10x,"HAMANN pseudopotential will be used. ")')
       pseudo_type = ABINIT
    else
       !if(inode==ionode.AND.iprint_init>0) &
       !     write(io_lun,fmt='(10x,"OLD pseudopotential will be used. ")')
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
    flag_ghost=.false.
    if(fdf_block('ChemicalSpeciesLabel')) then
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
    !maxnsf      = 0
    !max_rc = zero
    min_blip_sp = 1.0e8_double
          flag_InitialAtomicSpin = .false.
    do i=1,n_species
       charge(i)         = zero
       charge_up(i)      = zero
       charge_dn(i)      = zero
       sum_elecN_spin    = zero
       nsf_species(i)    = 0
       RadiusSupport(i)  = zero
       RadiusAtomf(i)    = zero
       RadiusMS(i)       = zero
       RadiusLD(i)       = zero
       NonLocalFactor(i) = zero
       InvSRange(i)      = zero
       blip_info(i)%SupportGridSpacing = zero
       if(pseudo_type==SIESTA.OR.pseudo_type==ABINIT) non_local_species(i) = .true.
       ! This is new-style fdf_block
       !nullify(bp)
       !if(fdf_block(species_label(i),bp)) then
       !   do while(fdf_bline(bp,line)) ! While there are lines in the block
       if(fdf_block(species_label(i))) then
          !charge(i)        = fdf_double ('Atom.ValenceCharge',zero)
          charge_up(i)     = fdf_double ('Atom.SpinNeUp',zero)
          charge_dn(i)     = fdf_double ('Atom.SpinNeDn',zero)
          sum_elecN_spin   = charge_up(i)+charge_dn(i)
          if (abs(sum_elecN_spin)>RD_ERR) then
             flag_InitialAtomicSpin = .true.
             ! We will check that the sum of charge_up and charge_dn matches charge later
          endif
          nsf_species(i)   = fdf_integer('Atom.NumberOfSupports',0)
          RadiusSupport(i) = fdf_double ('Atom.SupportFunctionRange',zero)
          !RadiusAtomf(i)   = RadiusSupport(i) ! = r_pao for (atomf=paof) or r_sf for (atomf==sf)
          ! Added DRB 2018/07/16 for safety
          if(flag_Multisite) then
             RadiusMS(i)      = fdf_double ('Atom.MultisiteRange',zero)
             RadiusLD(i)      = fdf_double ('Atom.LFDRange',RadiusMS(i))
             if(RadiusLD(i)<RadiusMS(i)) call cq_warn(sub_name,"LFD range should be larger than MSSF range: ", &
                  RadiusLD(i),RadiusMS(i))
          end if
          ! Moved to ... so that RadiusAtomf is read from ion files
          !if (flag_Multisite) RadiusSupport(i) = RadiusAtomf(i) + RadiusMS(i)
          ! DRB 2016/08/05 Keep track of maximum support radius
          !if(RadiusSupport(i)>max_rc) max_rc = RadiusSupport(i)
          ! SYM 2017/03/13 Fix no index on InvSRange
          InvSRange(i)     = fdf_double ('Atom.InvSRange',zero)
          !if(InvSRange(i)<RD_ERR) InvSRange(i) = RadiusSupport(i)
          NonLocalFactor(i) = fdf_double('Atom.NonLocalFactor',zero)
          if(NonLocalFactor(i)>one.OR.NonLocalFactor(i)<zero) then
             call cq_warn(sub_name, "0.0 < Atom.NonLocalFactor < 1.0; you have ",NonLocalFactor(i))
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
       else if(flag_Multisite.OR.(flag_basis_set==blips)) then
          call cq_abort("Failure to read data for species_label: "&
               //species_label(i))
       end if
       !%%! if(nsf_species(i)==0) &
       !%%!      call cq_abort("Number of supports not specified for species ",i)
       if(flag_basis_set==blips.AND.blip_info(i)%SupportGridSpacing<RD_ERR) &
            call cq_abort("Error: for a blip basis set you must &
            &set SupportGridSpacing for all species")
       !maxnsf = max(maxnsf,nsf_species(i))
       !if(RadiusSupport(i)<RD_ERR) &
       !     call cq_abort("Radius of support too small for &
       !                    &species; increase SupportFunctionRange ",i)
    end do
    ! For ghost atoms
    if(flag_ghost) then
       !2019Dec30 tsuyoshi
       !the following do-loop should be commented out, since charge(:) will be updated in setup_pseudo_info
       do i=1, n_species
          if(type_species(i) < 0) then
             charge(i) = zero
          endif
       enddo

       ! At present, we cannot use neutral atom potential for ghost atoms.
       !    2018/Sep/26   Tsuyoshi Miyazaki
       ! Fixed for Neutral Atom Potential  2018/Nov/16  TM
       !   but not for neutral atom projectors ...
       if(flag_neutral_atom_projector) then
          if(inode==ionode) write(io_lun,'(10x,A)') &
               'Now, we cannot use neutral atom projector for ghost atoms. Please set General.NeutralAtomProjector as false.'
          call cq_abort("Error: Neutral Atom Projector cannot be used for ghost atom, at present.")
       endif
    endif
          if(flag_InitialAtomicSpin) then
            ! we may not need the following do-loop since charge(:) has not been defined yet, here.
            do i = 1, n_species
             charge_tmp = charge_up(i) + charge_dn(i)
             if(charge_tmp < RD_ERR) then   ! includes ghost atom?
               charge_up(i) = half*charge(i)
               charge_dn(i) = half*charge(i)
             endif
            enddo !i = 1, n_species
           !if(inode .eq. ionode) then
           ! write(io_lun,*) 'TMTMTM: flag_InitialAtomicSpin = ',flag_InitialAtomicSpin
           ! do i = 1, n_species
           !    write(io_lun,*) ' ispecies, charge_up and charge_dn = ',i, charge_up(i), charge_dn(i)
           ! enddo
           !endif
          endif
!!$
!!$
!!$
!!$
!!$
    ! Set variables for multisite support functions
    if (flag_Multisite) then
       flag_MSSF_nonminimal = fdf_boolean('Multisite.nonminimal', .false.) !nonmin_mssf
       MSSF_nonminimal_offset = fdf_double('Multisite.nonminimal.offset', 0.0_double) !nonmin_mssf
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
       flag_LFD_nonSCF   = .true.
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
       flag_LFD_nonSCF = fdf_boolean('Multisite.LFD.NonSCF',.false.)
       flag_mix_LFD_SCF = fdf_boolean('Multisite.LFD.MixLFDSCF',.true.)
       if (.NOT.flag_LFD_nonSCF) then ! Expected behaviour
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
    ! Find out what type of run we're doing
    runtype             = fdf_string(20,'AtomMove.TypeOfRun',       'static')
    if(leqi(runtype,'pulay')) then
       runtype = 'sqnm'
       if(inode==ionode) write(io_lun,fmt='(/4x,"Pulay relaxation superceded by SQNM; changing method.")')
    end if
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
    flag_buffer_old       = fdf_boolean('AtomMove.OldBuffer',        .false.)
    AtomMove_buffer       = fdf_double ('AtomMove.BufferSize',    4.0_double)
    flag_pulay_simpleStep = fdf_boolean('AtomMove.PulaySimpleStep',  .false.)
    CGreset               = fdf_boolean('AtomMove.ResetCG',          .false.)
    MDn_steps             = fdf_integer('AtomMove.NumSteps',     100        )
    MDfreq                = fdf_integer('AtomMove.OutputFreq',    50        )
    XSFfreq               = fdf_integer('AtomMove.XsfFreq',    MDfreq        )
    if (leqi(runtype,'md')) then
      XYZfreq             = fdf_integer('AtomMove.XyzFreq',    MDfreq        )
    else
      XYZfreq             = fdf_integer('AtomMove.XyzFreq',    1        )
    end if
    MDtimestep            = fdf_double ('AtomMove.Timestep',      0.5_double)
    MDcgtol               = fdf_double ('AtomMove.MaxForceTol',0.0005_double)
    sqnm_trust_step       = fdf_double ('AtomMove.MaxSQNMStep',0.2_double   )
    LBFGS_history         = fdf_integer('AtomMove.LBFGSHistory', 5          )
    flag_opt_cell         = fdf_boolean('AtomMove.OptCell',          .false.)
    ! At present (2023/07/26 just before v1.2 release) neutral atom is required for cell opt
    if(flag_opt_cell.and.(.not.flag_neutral_atom)) &
         call cq_abort("You must use neutral atom for cell optimisation")
    ! This can be removed when ewald update is implemented
    flag_variable_cell    = flag_opt_cell
    optcell_method        = fdf_integer('AtomMove.OptCellMethod', 1)
    cell_constraint_flag  = fdf_string(20,'AtomMove.OptCell.Constraint','none')
    cell_en_tol           = fdf_double('AtomMove.OptCell.EnTol',0.00001_double)
    ! It makes sense to use GPa here so I'm changing the default to 0.1GPa
    cell_stress_tol       = fdf_double('AtomMove.StressTolerance',0.1_double) !005_double)
    flag_stop_on_empty_bundle = fdf_boolean('AtomMove.StopOnEmptyBundle',.false.)
    enthalpy_tolerance    = fdf_double('AtomMove.EnthalpyTolerance', 1.0e-5_double)
    flag_stress           = fdf_boolean('AtomMove.CalcStress', .true.)
    flag_full_stress      = fdf_boolean('AtomMove.FullStress', .false.)
    flag_atomic_stress    = fdf_boolean('AtomMove.AtomicStress', .false.)
    mix_input_output_XC_GGA_stress = fdf_double('General.MixXCGGAInOut',half)
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
    if(flag_mix_L_SC_min .and. flag_diagonalisation) then
       flag_mix_L_SC_min = .false.
    else if(flag_mix_L_SC_min .and. flag_self_consistent) then
       flag_self_consistent = .false.
    end if
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
    flag_continue_on_SC_fail     = fdf_boolean('SC.ContinueOnSCFail',   .false.)
    maxitersSC      = fdf_integer('SC.MaxIters',50)
    minitersSC      = fdf_integer('SC.MinIters',0) ! Changed default 2->0 DRB 2018/02/26
    maxearlySC      = fdf_integer('SC.MaxEarly',3 )
    maxpulaySC      = fdf_integer('SC.MaxPulay',5 )
    ! New residual flags jtlp 08/2019
    flag_newresidual = fdf_boolean('SC.AbsResidual', .false.)
    flag_newresid_abs = fdf_boolean('SC.AbsResidual.Fractional', .true.)
    ! When constructing charge density from K matrix at last step, we check the total
    ! number of electrons. If the error of electron number (per total electron number)
    ! is larger than the following value, we use atomic charge density. (in update_H)
    threshold_resetCD     = fdf_double('SC.Threshold.Reset',0.1_double)
    tmp = fdf_string(4,'AtomMove.CGLineMin','safe')
    if(leqi(tmp,'safe')) then
       cg_line_min = safe
    else if(leqi(tmp,'back')) then
       cg_line_min = backtrack
    else if(leqi(tmp,'adap')) then
       cg_line_min = adapt_backtrack
    end if
!!$
!!$
!!$   New Parameters for Dumping Files
!!$
!!$

    n_dumpL              = fdf_integer('IO.DumpFreq.DMM',0)
    n_dumpSCF            = fdf_integer('IO.DumpFreq.SCF',0)
    n_dumpSFcoeff        = fdf_integer('IO.DumpFreq.SFcoeff',0)
!!$
    flag_DumpMatrices      = fdf_boolean('IO.DumpMatrices',.true.)
    flag_DumpChargeDensity = fdf_boolean('IO.DumpChargeDensity',.false.)

!!$
!!$
!!$
!!$
!!$
    ! Name and Format of Matrix Files  (see store_matrix)
    flag_MatrixFile_RankFromZero = fdf_boolean('IO.MatrixFile.RankFromZero', .true.)
    flag_MatrixFile_BinaryFormat = fdf_boolean('IO.MatrixFile.BinaryFormat', .true.)
    flag_MatrixFile_BinaryFormat_Grab = fdf_boolean('IO.MatrixFile.BinaryFormat.Grab', flag_MatrixFile_BinaryFormat)
    flag_MatrixFile_BinaryFormat_Dump_END = fdf_boolean('IO.MatrixFile.BinaryFormat.Dump', flag_MatrixFile_BinaryFormat)
    flag_MatrixFile_BinaryFormat_Dump = flag_MatrixFile_BinaryFormat_Grab

    ! To load Matrix Files having nonzero indices
    index_MatrixFile       = fdf_integer('IO.MatrixFile.FileIndex',0)

    ! read wavefunction output flags
    mx_temp_matrices = fdf_integer('General.MaxTempMatrices',100)
    wf_self_con=.false.
    flag_out_wf=fdf_boolean('IO.outputWF',.false.)
    if (flag_out_wf) then
       if (flag_diagonalisation .and. leqi(runtype,'static')) then
          ! The user can either specify which bands explicitly
          max_wf=fdf_integer('IO.maxnoWF',0)
          if(max_wf>0) then
             allocate(out_wf(max_wf))
             if (fdf_block('WaveFunctionsOut')) then
                if(1+block_end-block_start<max_wf) &
                     call cq_abort("Too few wf no in WaveFunctionsOut:"&
                     ,1+block_end-block_start,max_wf)
                do i=1,max_wf
                   read(unit=input_array(block_start+i-1),fmt=*) out_wf(i)
                end do
                call fdf_endblock
             else
                call cq_abort("You specified bands with IO.maxnoWF but didn't give the WaveFunctionsOut block")
             end if
             flag_wf_range_Ef = .false.
          else
             ! Or specify an energy range
             E_wf_min = fdf_double('IO.min_wf_E',zero)
             E_wf_max = fdf_double('IO.max_wf_E',zero)
             ! Is the range relative to Ef (T) or absolute (F)
             flag_wf_range_Ef = fdf_boolean('IO.WFRangeRelative',.true.)
             if(flag_wf_range_Ef.AND.abs(E_wf_min)<very_small.AND.abs(E_wf_max)<very_small) then
                flag_out_wf = .false.
                flag_wf_range_Ef = .false.
                if(inode==ionode) write(io_lun,'(2x,"Setting IO.outputWF F as no bands range given")')
             end if
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
          if(flag_write_projected_DOS) then
             E_wf_min = fdf_double('IO.min_wf_E',-BIG)
             E_wf_max = fdf_double('IO.max_wf_E',BIG)
          end if
          ! Possibly needed to decide if MSSF needs dealing with
          !flag_pDOS_angmom = fdf_boolean('IO.PDOS_Angmom',.false.)
          !if (flag_pDOS_angmom .and. flag_basis_set==blips) then
          !   flag_pDOS_angmom = .false.
          !   if(inode==ionode) write(io_lun,'(2x,"Setting IO.PDOS_Angmom F as using blips")')
          !endif
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
          call cq_warn(sub_name, "Becke weights required for cDFT ! Setting to true")
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
    ! Calculate bulk polarisation
    flag_calc_pol   = fdf_boolean('General.CalcPol', .false.)
    flag_do_pol_calc = .false.
    ! Find direction for polarisation calculation: 0 means all three
    i_pol_dir = 0
    i_pol_dir(1) = fdf_integer('General.PolDir',0)
    i_pol_dir_st = 1
    i_pol_dir_end = 1
    if(i_pol_dir(1)==0) then
       i_pol_dir_end = 3
       i_pol_dir(1) = 1
       i_pol_dir(2) = 2
       i_pol_dir(3) = 3
    else if(i_pol_dir(1)>3) then
       call cq_abort("Illegal value for General.PolDir: must lie between 0 and 3 ",i_pol_dir(1))
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
    else if(flag_quench_MD) then
       temp_ion           = fdf_double ('AtomMove.IonTemperature',0.0_double)
       fire_alpha0        = fdf_double ('AtomMove.QMDAlpha',0.5_double)
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
       if(flag_which_force>10.OR.flag_which_force<0.AND.inode==ionode) then
          call cq_warn(sub_name,&
               "AtomMove.TestSpecificForce must lie between 1 and 10 (setting to 1): ",&
               flag_which_force)
          flag_which_force = 1
       end if
    end if
    TF_direction = fdf_integer('AtomMove.TestForceDirection',1)
    if(TF_direction>3.OR.TF_direction<0.AND.inode==ionode) then
       call cq_warn(sub_name,&
            "TestForceDirection must lie between 1 and 3: ",TF_direction)
       TF_direction = 1
    end if
    TF_atom_moved = fdf_integer('AtomMove.TestForceAtom',1)
    TF_delta      = fdf_double ('AtomMove.TestForceDelta',0.00001_double)
!!$
!!$
!!$
!!$
    flag_functional_type = fdf_integer('General.FunctionalType', 0)   ! Read from pseudopotentials
    flag_different_functional = fdf_boolean('General.DifferentFunctional',.false.) ! Use different functional to ion files
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
       call cq_warn(sub_name,"Unrecognised partitioning method: "//part_mode//" Using Hilbert")
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
    gap_threshold = fdf_double('General.GapThreshold',zero) !two*max_rc)
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
    flag_MDdebug      = fdf_boolean('AtomMove.Debug',.false.)
    flag_MDcontinue   = fdf_boolean('AtomMove.RestartRun',.false.)
    flag_LmatrixReuse = fdf_boolean('AtomMove.ReuseDM',.true.)
    flag_SFcoeffReuse = fdf_boolean('AtomMove.ReuseSFcoeff',flag_LmatrixReuse)
    flag_write_xsf    = fdf_boolean('AtomMove.WriteXSF', .true.)
    flag_write_extxyz = fdf_boolean('AtomMove.WriteExtXYZ', .false.)
    ! tsuyoshi 2019/12/30
    if(flag_SFcoeffReuse .and. .not.flag_LmatrixReuse) then
       call cq_warn(sub_name,' AtomMove.ReuseDM should be true if AtomMove.ReuseSFcoeff is true.')
       flag_LmatrixReuse = .true.
    endif
    if(flag_Multisite) then
       if(.not.flag_SFcoeffReuse .and. flag_LmatrixReuse) then
          call cq_warn(sub_name,' AtomMove.ReuseSFcoeff should be true if AtomMove.ReuseDM is true.')
          flag_SFcoeffReuse = .true.
       endif
    endif

    ! tsuyoshi 2019/12/27
    !  New Keyword for the method to update the charge density after the movement of atoms
    !    DensityMatrix = 0; AtomicCharge = 1; LastStep = 2
    method_UpdateChargeDensity = fdf_integer('AtomMove.InitialChargeDensity',DensityMatrix)

    !  The keywords ( SC.ResetDensOnAtomMove and Multisite.LFD.UpdateWithAtomicDensity )
    !  should be removed in the near future.
    flag_reset_dens_on_atom_move = fdf_boolean('SC.ResetDensOnAtomMove',.false.)
    if(flag_reset_dens_on_atom_move) then
       call cq_warn(sub_name,' SC.ResetDensOnAtomMove will not be available soon. &
            &             Set AtomMove.InitialChargeDensity as 1, instead.')
       method_UpdateChargeDensity = AtomicCharge
    endif

    if (flag_LFD .and. .not.flag_SFcoeffReuse) then
       ! if LFD, use atomic density in default when we don't reuse SFcoeff
       flag_LFD_MD_UseAtomicDensity = fdf_boolean('Multisite.LFD.UpdateWithAtomicDensity',.true.)
       call cq_warn(sub_name,' Multisite.LFD.UpdateWithAtomicDensity will not be available soon. &
            &             Set AtomMove.InitialChargeDensity, instead.')
    endif
    if(flag_LFD_MD_UseAtomicDensity) method_UpdateChargeDensity = AtomicCharge

    if(method_UpdateChargeDensity == AtomicCharge) then
       flag_reset_dens_on_atom_move = .true.
       flag_LFD_MD_UseAtomicDensity = .true.
    endif


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
       restart_DM         = fdf_boolean('General.LoadDM', .true.)
       find_chdens    = fdf_boolean('SC.MakeInitialChargeFromK',.true.)
       if (flag_XLBOMD) restart_X=fdf_boolean('XL.LoadX', .true.)
       if (flag_multisite .or. flag_basis_set==blips) read_option = fdf_boolean('Basis.LoadCoeffs', .true.)
    else
       flag_read_velocity = fdf_boolean('AtomMove.ReadVelocity',.false.)
       restart_DM         = fdf_boolean('General.LoadDM', .false.)
       if(restart_DM) then
          find_chdens    = fdf_boolean('SC.MakeInitialChargeFromK',.true.)
       else
          find_chdens    = fdf_boolean('SC.MakeInitialChargeFromK',.false.)
       endif
       if (flag_XLBOMD) restart_X=fdf_boolean('XL.LoadX', .false.)
       if (flag_multisite .or. flag_basis_set==blips) read_option = fdf_boolean('Basis.LoadCoeffs', .false.)
    end if

    if (restart_DM .and. flag_Multisite .and. .not.read_option) then
       call cq_abort("When L or K matrix is read from files, SFcoeff also must be read from files for multi-site calculation.")
    endif
    if(find_chdens .and. (.not.restart_DM)) then
       call cq_warn(sub_name," Cannot make charge density from K without loading K! Starting from atomic densities.")
       find_chdens = .false.
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
       maxiter_Dissipation = fdf_integer('XL.MaxDissipation',5)
       if (flag_dissipation) then
          if (maxiter_Dissipation.LT.1 .OR. maxiter_Dissipation.GT.9) then
             maxiter_Dissipation = 5
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
    select case(md_ensemble)
    case('nve')
       md_thermo_type     = fdf_string(20, 'MD.Thermostat', 'none')
       md_baro_type       = fdf_string(20, 'MD.Barostat', 'none')
    case('nvt')
       md_thermo_type     = fdf_string(20, 'MD.Thermostat', 'nhc')
       md_baro_type       = fdf_string(20, 'MD.Barostat', 'none')
    case('npt')
       md_thermo_type     = fdf_string(20, 'MD.Thermostat', 'nhc')
       md_baro_type       = fdf_string(20, 'MD.Barostat', 'pr')
       flag_variable_cell  = .true.
    end select
    md_tau_T           = fdf_double('MD.tauT', -one)
    md_tau_T_equil     = fdf_double('MD.tauTEquil', one)
    md_n_nhc           = fdf_integer('MD.nNHC', 5)
    md_n_ys            = fdf_integer('MD.nYoshida', 1)
    md_n_mts           = fdf_integer('MD.nMTS', 1)
    flag_thermoDebug   = fdf_boolean('MD.ThermoDebug',.false.)
    md_t_drag          = fdf_double('MD.TDrag', zero)
    if (leqi(md_thermo_type, 'nhc')) then
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
    end if
    flag_heat_flux = fdf_boolean('MD.HeatFlux', .false.)

    ! Barostat
    target_pressure    = fdf_double('AtomMove.TargetPressure', zero)
    md_box_mass        = fdf_double('MD.BoxMass', one)
    md_tau_P           = fdf_double('MD.tauP', -one)
    md_tau_P_equil     = fdf_double('MD.tauPEquil', 100.0_double)
    md_bulkmod_est     = fdf_double('MD.BulkModulusEst', 100.0_double)
    md_cell_nhc        = fdf_boolean('MD.CellNHC', .true.)
    flag_baroDebug     = fdf_boolean('MD.BaroDebug',.false.)
    md_equil_steps     = fdf_integer('MD.EquilSteps', 0)
    if(.NOT.leqi(md_thermo_type,'svr').AND.md_equil_steps>0) then
       call cq_warn("read_input","MD equilibration only possible with SVR")
       md_equil_steps = 0
    end if
    md_equil_press     = fdf_double('MD.EquilPress',one) ! GPa
    md_equil_press     = md_equil_press/HaBohr3ToGPa
    md_tdep            = fdf_boolean('MD.TDEP', .false.)
    md_p_drag          = fdf_double('MD.PDrag', zero)
    md_cell_constraint = fdf_string(20, 'MD.CellConstraint', 'volume')

    !**** TM 2017.Nov.3rd
    !call check_compatibility
    !****lat<$
    call stop_backtrace(t=backtrace_timer,who='read_input')
    !****lat>$

    call my_barrier()

    return
  end subroutine read_input
  !!***

  ! ------------------------------------------------------------------------------
  ! subroutine check_compatibility
  ! ------------------------------------------------------------------------------
  !!****f* initial_read/check_compatibility *
  !!
  !!  NAME
  !!   check_compatibility
  !!  USAGE
  !!
  !!  PURPOSE
  !!   checks the compatibility between keywords mainly defined in read_input
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   T. Miyazaki
  !!  CREATION DATE
  !!   2017/11/03
  !!  MODIFICATION HISTORY
  !!   2020/07/30 tsuyoshi
  !!    - Moved from read_input
  !!  SOURCE
  !!


    ! this subroutine checks the compatibility between keywords defined in read_input
    !       2017.11(Nov).03   Tsuyoshi Miyazaki
    !
    ! we don't need to worry about which parameter is defined first.
    !
    subroutine check_compatibility
      use global_module, only: flag_move_atom, ni_in_cell, &
                               runtype, flag_XLBOMD, flag_diagonalisation, &
                               flag_LmatrixReuse, flag_SFcoeffReuse, flag_DumpMatrices, &
                               flag_FixCOM
      use input_module,  only: leqi
      use density_module,only: method_UpdateChargeDensity,DensityMatrix,AtomicCharge
      use GenComms,      only: cq_warn

      implicit none
      logical :: flag_FixedAtoms
      integer :: ig, k

      character(len=80) :: sub_name = "check_compatibility"

      !check AtomMove.ExtendedLagrangian(flag_XLBOMD) &  AtomMove.TypeOfRun (runtype)
      if((.NOT.leqi(runtype,'md')) .and. flag_XLBOMD) then
         flag_XLBOMD=.false.
         call cq_warn(sub_name,&
              '(AtomMove.ExtendedLagrangian): XLBOMD should be used only with MD.')
      endif
      !check AtomMove.ExtendedLagrangian(flag_XLBOMD) &  DM.SolutionMethod
      if(flag_diagonalisation .and. flag_XLBOMD) then
         flag_XLBOMD=.false.
         call cq_warn(sub_name,&
              'WARNING (AtomMove.ExtendedLagrangian): present XLBOMD is only for orderN')
      endif

      !flag_LmatrixReuse & method_UpdateChargeDensity
      if(.not.flag_LmatrixReuse .and. method_UpdateChargeDensity == DensityMatrix) then
         method_UpdateChargeDensity = AtomicCharge
         call cq_warn(sub_name,&
              'AtomMove.InitialChargeDensity is changed to AtomicCharge, since AtomMove.ReuseDM is false')
      endif

      !flag_DumpMatrices : at present, we need matrix files to reuse the previous matrix data
      if(.not.flag_DumpMatrices) then
         if(flag_SFcoeffReuse .or. flag_LmatrixReuse) then
            flag_DumpMatrices = .true.
            call cq_warn(sub_name,&
                 'IO.DumpMatrices must be true when AtomMove.ReuseSFcoeff or AtomMove.ReuseDM is true')
         endif
      endif

      !flag_FixCOM  &  flag_move_atom  2020/Jul/30

      flag_FixedAtoms = .false.
      do ig = 1, ni_in_cell
       do k = 1,3
        if(.NOT.flag_move_atom(k,ig)) flag_FixedAtoms = .true.
       enddo
      enddo

      if(flag_FixedAtoms .and. flag_FixCOM) then
         flag_FixCOM = .false.
         call cq_warn(sub_name,&
           'flag_FixCOM should be false when some of the atomic positions are fixed')
      endif

      return
    end subroutine check_compatibility
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
  !!   2019/06/06 18:00 nakata
  !!    - Added MSSF_nonminimal_species
  !!  SOURCE
  !!
  subroutine allocate_species_vars

    use dimens,             only: RadiusSupport, RadiusAtomf, RadiusMS, RadiusLD, &
         NonLocalFactor, InvSRange, atomicnum
    use memory_module,      only: reg_alloc_mem, type_dbl
    use species_module,     only: n_species, nsf_species, nlpf_species, npao_species, natomf_species, &
         charge, charge_up, charge_dn
    use species_module,     only: mass, non_local_species, ps_file, ch_file, phi_file
    use species_module,     only: species_label, species_file, type_species
    use global_module,      only: area_general
    use GenComms,           only: cq_abort
    use blip,               only: blip_info
    use multisiteSF_module, only: MSSF_nonminimal_species

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
    allocate(MSSF_nonminimal_species(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating MSSF_nonminimal_species in allocate_species_vars: ", n_species,stat)
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
  !!   2022/10/28 lionel
  !!    Add printing of species table in ASE output
  !!   2022/12/14 11:31 dave
  !!    Removed output of support grid spacing with blips (shouldn't be here: species dependent)
  !!  SOURCE
  !!
  subroutine write_info(titles, mu, vary_mu, HNL_fac, NODES)

    use datatypes
    use units
    use dimens,               only: r_super_x, r_super_y, r_super_z,   &
         n_grid_x, n_grid_y, n_grid_z, r_h, r_c, RadiusSupport
    use block_module,         only: in_block_x, in_block_y, in_block_z
    use species_module,       only: n_species, species_label, mass,    &
         charge,          &
         phi_file, nsf_species
    use pseudopotential_data, only: core_radius
    use DiagModule,           only: flag_smear_type,           &
         iMethfessel_Paxton
    use blip,                 only: blip_info
    use global_module,        only: flag_basis_set, blips,        &
         flag_precondition_blips, io_lun, io_ase, write_ase, flag_LFD, runtype, flag_opt_cell, optcell_method, &
         flag_Multisite, flag_diagonalisation, flag_neutral_atom, temp_ion, &
         flag_self_consistent, flag_vary_basis, iprint_init, flag_pcc_global, &
         nspin, flag_SpinDependentSF, flag_fix_spin_population, ne_spin_in_cell, flag_XLBOMD,&
         ase_file
    use SelfCon,              only: maxitersSC
    use GenComms,             only: cq_abort
    use minimise,             only: energy_tolerance, L_tolerance,     &
         sc_tolerance,                      &
         n_support_iterations,              &
         n_L_iterations
    use datestamp,            only: datestr, commentver
    use pseudopotential_common, only: flag_neutral_atom_projector, maxL_neutral_atom_projector, &
         numN_neutral_atom_projector, pseudo_type, OLDPS, SIESTA, ABINIT
    use input_module,         only: leqi, chrcap
    use control,    only: MDn_steps
    use md_control, only: md_ensemble
    use omp_module, only: init_threads
    use multiply_kernel, only: kernel_id

    implicit none

    ! Passed variables
    logical :: vary_mu
    character(len=80) :: titles
    character(len=3) :: ensemblestr
    integer :: NODES
    real(double) :: mu, HNL_fac

    ! Local variables
    integer :: n, stat
    integer :: threads
    character(len=10) :: today, the_time
    character(len=15) :: job_str
    character(len=5)  :: timezone

    call date_and_time(today, the_time, timezone)
3   format()
    write(io_lun,fmt='(/4x,"This job was run on  ",a4,"/",a2,"/",a2," at ",a2,":",a2," ",a5)') &
         today(1:4), today(5:6), today(7:8), the_time(1:2), the_time(3:4), timezone
    write(io_lun,&
         '(4x,"Code was compiled on ",a,/4x,"Version comment: ",a)') &
         datestr, commentver

    write(io_lun,fmt='(/4x,"Job title: ",a)') titles

    ! Job type
    job_str = "Job to be run: "
    if(leqi(runtype,'static') ) then
       write(io_lun, fmt='(4x,a15,"static calculation")') job_str
    else if(leqi(runtype,'cg')) then
       if(flag_opt_cell) then
          if(optcell_method==1) then
             write(io_lun, fmt='(4x,a15,"CG cell relaxation")') job_str
          else
             write(io_lun, fmt='(4x,a15,"Combined CG atomic and cell relaxation")') job_str
          end if
       else
          write(io_lun, fmt='(4x,a15,"CG atomic relaxation")') job_str
       end if
    else if(leqi(runtype,'sqnm')) then
       if(flag_opt_cell) then
          if(optcell_method==1) then
             write(io_lun, fmt='(4x,a15,"CG cell relaxation")') job_str
          else if(optcell_method==2) then
             write(io_lun, fmt='(4x,a15,"Combined SQNM atomic and CG cell relaxation")') job_str
          else
             write(io_lun, fmt='(4x,a15,"Combined CG atomic and cell relaxation")') job_str
          end if
       else
          write(io_lun, fmt='(4x,a15,"SQNM atomic relaxation")') job_str
       end if
    else if(leqi(runtype,'md')) then
       ensemblestr = md_ensemble
       call chrcap(ensemblestr,3)
       write(io_lun, fmt='(4x,a15,a3," MD run for ",i5," steps ")') job_str, ensemblestr, MDn_steps
       write(io_lun, fmt='(6x,"Initial ion temperature: ",f9.3,"K")') temp_ion
       if(flag_XLBOMD) write(io_lun, fmt='(6x,"Using extended Lagrangian formalism")')
    else if(leqi(runtype,'lbfgs')) then
       write(io_lun, fmt='(4x,a15,"L-BFGS atomic relaxation")') job_str
    end if
    ! Ground state search details
    write(io_lun,fmt='(/4x,"Ground state search:")')
    if(flag_basis_set==blips) then 
       write(io_lun,'(6x,"Support functions represented with blip basis")')
    else
       write(io_lun,'(6x,"Support functions represented with PAO basis")') 
       if(flag_Multisite) then
          if(flag_LFD) then
             write(io_lun,'(6x,"Multi-site SFs used with local filter diagonalisation")')
          else
             write(io_lun,'(6x,"Multi-site SFs used")')
          end if
          if(flag_SpinDependentSF) write(io_lun,'(6x,"SFs are spin dependent")')
       else
          write(io_lun,'(6x,"1:1 PAO to SF mapping")')
       end if
    end if
    if(nspin==2) then
       if(flag_fix_spin_population) then
          write(io_lun,'(6x,"Spin-polarised electrons; population difference fixed at ",f10.6)') &
               ne_spin_in_cell(1) - ne_spin_in_cell(2)
       else
          write(io_lun,'(6x,"Spin-polarised electrons; population difference free")')
       end if
    else
       write(io_lun,'(6x,"Non-spin-polarised electrons")')
    end if
    if(flag_diagonalisation) then
       write(io_lun,fmt='(6x,"Solving for the K matrix using ",a16)') 'diagonalisation '
       if(iprint_init>0) then
          select case (flag_smear_type)
          case (0)
             write(io_lun,'(6x,"Using Fermi-Dirac smearing")')
          case (1)
             write(io_lun,&
                  '(6x,"Using order ",i2," Methfessel-Paxton smearing")') &
                  iMethfessel_Paxton
          end select
       end if
    else
       write(io_lun,fmt='(6x,"Solving for the K matrix using ",a16)') 'order N with LNV'   
    end if
    if(iprint_init>0) write(io_lun,fmt='(/4x,"Integration grid size:   ",i4," x ",i4," x ",i4)') &
         n_grid_x, n_grid_y, n_grid_z

    if(iprint_init>1) write(io_lun,fmt='(4x, "Integration grid blocks: ",i4," x ",i4," x ",i4)') &
         in_block_x, in_block_y, in_block_z

    write(io_lun,fmt='(/4x,"Integration grid spacing: ",f6.3,a3," x",f6.3,a3," x",f6.3,a3)') &
         dist_conv*(r_super_x/n_grid_x), d_units(dist_units), & 
         dist_conv*(r_super_y/n_grid_y), d_units(dist_units), &
         dist_conv*(r_super_z/n_grid_z), d_units(dist_units)
    
    ! print in Conquest output
    write(io_lun,fmt='(/4x,"Number of species: ",i2)') n_species
    write(io_lun,fmt='(4x,a56)') '--------------------------------------------------------'
    write(io_lun,fmt='(4x,"|   #  mass (au)  Charge (e)  SF Rad (",a2,")  NSF  Label  |")') &
         d_units(dist_units)
    write(io_lun,fmt='(4x,a56)') '--------------------------------------------------------'

    do n=1, n_species
       write(io_lun,fmt='(4x,"|",i4,2x,f9.3,3x,f9.3,4x,f9.3,2x,i3,2x,a7,"|")') &
            n, mass(n), charge(n), dist_conv*RadiusSupport(n), nsf_species(n), species_label(n)
    end do
    write(io_lun,fmt='(4x,a56)') '--------------------------------------------------------'
    
    !
    ! BEGIN %%%% ASE printing %%%%
    !
    if ( write_ase ) then
       open(io_ase,file=ase_file, status='old', action='readwrite', iostat=stat, position='append')
       if (stat .ne. 0) call cq_abort('ASE/io_module error opening file !')
       !
       write(io_ase,fmt='(/4x,"Number of species: ",i2)') n_species
       write(io_ase,fmt='(4x,a56)') '--------------------------------------------------------'
       write(io_ase,fmt='(4x,"|   #  mass (au)  Charge (e)  SF Rad (",a2,")  NSF  Label  |")') &
            d_units(dist_units)
       write(io_ase,fmt='(4x,a56)') '--------------------------------------------------------'
       
       do n=1, n_species
          write(io_ase,fmt='(4x,"|",i4,2x,f9.3,3x,f9.3,4x,f9.3,2x,i3,2x,a7,"|")') &
               n, mass(n), charge(n), dist_conv*RadiusSupport(n), nsf_species(n), species_label(n)
       end do
       write(io_ase,fmt='(4x,a56)') '--------------------------------------------------------'
       write(io_ase,fmt='(4x,a)') 'end of species report'
       !
       close(io_ase)
    end if
    !
    ! END %%%% ASE printing %%%%
    !
    if (flag_Multisite) write(io_lun,'(/10x,"PAOs are contracted to multi-site support functions")')

    if(iprint_init>0) then
       if(flag_vary_basis) then
          write(io_lun,fmt='(10x,"Support function tolerance:  ",f12.8)') energy_tolerance
          write(io_lun,fmt='(10x,"Support function iterations: ",i4)') n_support_iterations
       end if
       if(flag_self_consistent) then
          write(io_lun,fmt='(10x,"SCF tolerance:               ",f12.8)') sc_tolerance
          write(io_lun,fmt='(10x,"SCF iterations:              ",i4)') maxitersSC
       end if
       if(.not.flag_diagonalisation) then
          write(io_lun,fmt='(10x,"O(N) tolerance:              ",f12.8)') L_tolerance
          write(io_lun,fmt='(10x,"O(N) iterations:             ",i4)') n_L_iterations
       end if
    end if
    if(iprint_init>1) then
       if(pseudo_type==SIESTA) write(io_lun,fmt='(4x,"SIESTA (TM) pseudopotential will be used. ")')
       if(pseudo_type==ABINIT) write(io_lun,fmt='(4x,"Hamann (ONCVPSP) pseudopotential will be used. ")')
    end if
    ! PCC
    if (iprint_init>2.AND.flag_pcc_global) &
         write (io_lun,fmt='(4x,a)') "Some species include partial core corrections (PCC)."
    if(flag_neutral_atom.and.iprint_init>2) then
       write(io_lun,fmt='(/13x,"Using neutral atom potential (NAP) formalism")')
       if(flag_neutral_atom_projector) then
          write(io_lun,fmt='(/13x,"Calculating 1- and 2-centre NAP integrals analytically")')
          write(io_lun,fmt='(13x,"Expanding 3-centre NAP integrals with ",i2," projectors up to l=",i2)') &
               maxval(numN_neutral_atom_projector),maxL_neutral_atom_projector
       end if
    end if

    if (.not.vary_mu) then
       write(io_lun,*) '          mu is constant'
       write(io_lun,fmt="(/10x,'The Chemical Potential mu is :',f7.4)") mu
    endif

    if(nodes>1) then
       write(io_lun,fmt="(/4x,'The calculation will be performed on ',i5,' processes')") NODES
    else
       write(io_lun,fmt="(/4x,'The calculation will be performed on ',i5,' process')") NODES
    end if

    call init_threads(threads)
    if(threads>1) then
       write(io_lun,fmt="(/4x,'The calculation will be performed on ',i5,' threads')") threads
    else if (threads==1) then
       write(io_lun,fmt="(/4x,'The calculation will be performed on ',i5,' thread')") threads
    end if
    write(io_lun,fmt='(/4x,"Using the ",a," matrix multiplication kernel")') kernel_id
    if(.NOT.flag_diagonalisation) &
         write(io_lun,fmt='(10x,"Density Matrix range  = ",f7.4,1x,a2)') &
         dist_conv*r_c, d_units(dist_units)

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
  !!   2019/12/04 08:25 dave
  !!    Added gamma-centred option and tweaked loops for MP mesh to be more logical
  !!    Removed iprint control on k-point output: this is always needed
  !!    Removed gcopy and myid checks
  !!   2019/12/05 08:12 dave
  !!    Bug fix: only write out on ionode
  !!   2022/14/06 lat
  !!    Added cq_warn when matrix size is a prime number
  !!    and set block_size_r = block_size_c = 1
  !!   2022/07/16 lionel
  !!    Added printing fractional k-points when read from block
  !!   2022/06/29 12:00 dave
  !!    Moved printing to capture default gamma point behaviour
  !!   2023/07/20 12:00 tsuyoshi
  !!    Implementing 1st version of Padding H and S matrices
  !!  SOURCE
  !!
  subroutine readDiagInfo

    use datatypes
    use functions,       only: is_prime
    use global_module,   only: iprint_init, rcellx, rcelly, rcellz,  &
         area_general, ni_in_cell, numprocs,   &
         species_glob, io_lun, io_ase, ase_file, write_ase, flag_calc_pol
    use numbers,         only: zero, one, two, pi, RD_ERR, half
    use GenComms,        only: cq_abort, cq_warn, gcopy
    use input_module
    use ScalapackFormat, only: proc_rows, proc_cols, block_size_r,   &
         block_size_c, proc_groups, matrix_size, flag_padH
    use DiagModule,      only: nkp, kk, wtk, kT, maxefermi,          &
         flag_smear_type, iMethfessel_Paxton,  &
         max_brkt_iterations, gaussian_height, &
         finess, NElec_less
    use energy,          only: SmearingType, MPOrder
    use memory_module,   only: reg_alloc_mem, reg_dealloc_mem,       &
         type_dbl
    use species_module,  only: nsf_species
    use units, only: en_conv, en_units, energy_units

    implicit none

    ! Local variables
    character(len=80) :: sub_name = "readDiagInfo"
    type(cq_timer) :: backtrace_timer
    integer        :: stat, i, j, k, nk_st, nkp_lines
    real(double)   :: a, sum, dkx, dky, dkz
    integer        :: proc_per_group
    logical        :: ms_is_prime
    
    ! k-point mesh type
    logical        :: mp_mesh, done, flag_lines_kpoints, flag_gamma, test_ase
    integer,      dimension(1:3)              :: mp
    real(double), dimension(1:3)              :: mp_shift
    real(double), allocatable, dimension(:,:) :: kk_tmp
    real(double), allocatable, dimension(:)   :: wtk_tmp
    integer :: nkp_tmp, nkp_in_line, inc
    integer :: counter
    character(len=2) :: suffix

    !****lat<$    
    call start_backtrace(t=backtrace_timer,who='readDiagInfo',where=1,level=2)
    !****lat>$

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
       call cq_warn(sub_name,"Too few processes for k-point groups (setting to 1): ", &
            numprocs, proc_groups)
       proc_groups = 1
    end if

    if (iprint_init > 1 .AND. inode==ionode) write (io_lun, 11) proc_groups 
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
          proc_rows = floor(a)+1
          do while(.NOT.done.AND.proc_rows>1) 
             proc_rows = proc_rows - 1
             proc_cols = floor(real(proc_per_group,double) / &
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
    
    ! Test if matrix_size is a prime number
    ms_is_prime = is_prime(matrix_size)
    if ( ms_is_prime ) call cq_warn(sub_name,'matrix size is a prime number', matrix_size)
    
    ! padH or not  :temporary?   
    flag_padH = fdf_boolean('Diag.PaddingHmatrix',.true.)

    if(fdf_defined('Diag.BlockSizeR')) then
       block_size_r = fdf_integer('Diag.BlockSizeR',1)
       block_size_c = fdf_integer('Diag.BlockSizeC',1)
       if(flag_padH) then
          if(block_size_c .ne. block_size_r) &
               call cq_abort('PaddingHmatrix: block_size_c needs to be block_size_r')
          block_size_c = block_size_r
       else
          a = real(matrix_size)/real(block_size_r)
          if(a - real(floor(a))>1e-8_double) &
               call cq_abort('block_size_r not a factor of matrix size ! ',&
               matrix_size, block_size_r)
          a = real(matrix_size)/real(block_size_c)
          if(a - real(floor(a))>1e-8_double) &
               call cq_abort('block_size_c not a factor of matrix size ! ',&
               matrix_size, block_size_c)
       endif
    else if (  ms_is_prime ) then
       block_size_r = 1
       block_size_c = block_size_r
       call cq_warn(sub_name,'prime: set block_size_c = block_size_r = 1 ')
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
    if(iprint_init>1.AND.inode==ionode) then
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
          if(iprint_init>1.AND.inode==ionode) then
             write(io_lun,fmt='(8x,"Number of Kpoint lines: ",i4)') nkp_lines
          else
             write(io_lun,fmt='(4x,"Using ",i3," lines of k-points specified by user")')
          end if
          if(nkp_lines<1) call cq_abort("Need to specify how many kpoint lines !",nkp_lines)
          nkp_in_line = fdf_integer('Diag.NumKpts',2)
          if(iprint_init>1.AND.inode==ionode) write(io_lun,fmt='(8x,"Number of Kpoints in a line: ",i4)') nkp_in_line
          ! Total number of k-points (potentially corrected later)
          nkp = nkp_in_line*nkp_lines
          ! Read start/end points for lines and correct for duplication
          allocate(kk_tmp(3,2*nkp_lines))
          kk_tmp = zero
          if(fdf_block('Diag.KpointLines'))then
             if(1+block_end-block_start<2*nkp_lines) &
                  call cq_abort("Kpoint line error: ",1+block_end-block_start,nkp_lines)
             do i=1,2*nkp_lines,2
                ! Start/end points for line
                read (unit=input_array(block_start+i-1),fmt=*) &
                     kk_tmp(1,i),kk_tmp(2,i),kk_tmp(3,i)
                read (unit=input_array(block_start+i),fmt=*) &
                     kk_tmp(1,i+1),kk_tmp(2,i+1),kk_tmp(3,i+1)
                ! If the start of this line duplicates the end of the last, reduce nkp by 1
                if(i>1) then
                   if(abs(kk_tmp(1,i)-kk_tmp(1,i-1))<RD_ERR .and. &
                        abs(kk_tmp(2,i)-kk_tmp(2,i-1))<RD_ERR .and. &
                        abs(kk_tmp(3,i)-kk_tmp(3,i-1))<RD_ERR) nkp = nkp - 1
                end if
             end do
          else
             call cq_abort("Must specify a block Diag.KpointLines to have lines of kpoints !")
          end if
          allocate(kk(3,nkp),wtk(nkp),STAT=stat)
          if(stat/=0) call cq_abort('FindEvals: couldnt allocate kpoints',nkp)
          call reg_alloc_mem(area_general,4*nkp,type_dbl)
          kk = zero
          nk_st = 1
          do i=1,2*nkp_lines,2
             ! Spacing for this line in three directions
             dkx = (kk_tmp(1,i+1) - kk_tmp(1,i))/real(nkp_in_line-1,double)
             dky = (kk_tmp(2,i+1) - kk_tmp(2,i))/real(nkp_in_line-1,double)
             dkz = (kk_tmp(3,i+1) - kk_tmp(3,i))/real(nkp_in_line-1,double)
             if(iprint_init>1.AND.inode==ionode) &
                  write(io_lun,fmt='(2x,"K-point spacing along line : ",i3,3f7.3)') i,dkx,dky,dkz
             ! Number of points in line
             inc = nkp_in_line
             if(i<2*nkp_lines-1) then ! If last point is the same as first point of next line, ignore
                if(abs(kk_tmp(1,i+2)-kk_tmp(1,i+1))<RD_ERR .and. &
                     abs(kk_tmp(2,i+2)-kk_tmp(2,i+1))<RD_ERR .and. &
                     abs(kk_tmp(3,i+2)-kk_tmp(3,i+1))<RD_ERR) then
                   inc = nkp_in_line - 1
                end if
             end if
             ! Initial point
             kk(1,nk_st) = kk_tmp(1,i)
             kk(2,nk_st) = kk_tmp(2,i)
             kk(3,nk_st) = kk_tmp(3,i)
             ! Intermediate points
             do j=1,inc-1
                kk(1,nk_st+j) = kk(1,nk_st+j-1)+dkx
                kk(2,nk_st+j) = kk(2,nk_st+j-1)+dky
                kk(3,nk_st+j) = kk(3,nk_st+j-1)+dkz
             end do
             nk_st = nk_st + inc
          end do
          deallocate(kk_tmp)
          wtk = one/real(nkp,double)
          ! Write out fractional k-points
          if(iprint_init>1.AND.inode==ionode) then
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
          if(iprint_init>1.AND.inode==ionode) then
             write(io_lun,fmt='(8x,"Number of Kpoints: ",i4)') nkp
          end if
          if(nkp<1) call cq_abort("Need to specify how many kpoints !",nkp)
          allocate(kk(3,nkp),wtk(nkp),STAT=stat)
          if(stat/=0) call cq_abort('FindEvals: couldnt allocate kpoints',nkp)
          call reg_alloc_mem(area_general,4*nkp,type_dbl)
          sum = zero
          ! Read k-points
          if(fdf_block('Diag.Kpoints'))then
             if(iprint_init==0) then
                if(nkp==1) then
                   write(io_lun,fmt='(4x,"Using ",i1," k-point specified by user")') nkp
                else
                   write(io_lun,fmt='(4x,"Using ",i3," k-points specified by user")') nkp
                endif
             end if
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
             if(inode==ionode) write(io_lun,fmt='(4x,"Default k-point sampling of Gamma point only")')
             nkp = 1
             kk(1,1) = zero
             kk(2,1) = zero
             kk(3,1) = zero
             wtk(1) = one
          end if
          ! Write out fractional k-points and weight (the easiest way not the cleverest)
          if(iprint_init>0.AND.inode==ionode) then
             write(io_lun,7) nkp
             do i=1,nkp
                write(io_lun,fmt='(8x,i5,3f15.6,f12.3)')&
                     i,kk(1,i)/(two*pi)*rcellx,kk(2,i)/(two*pi)*rcellx,kk(3,i)/(two*pi)*rcellx,wtk(i)
             end do
          end if
       end if
    else
       ! Read Monkhorst-Pack mesh coefficients
       ! Default is Gamma point only 
       if(iprint_init>1.AND.inode==ionode) &
            write(io_lun,fmt='(/8x,"Reading Monkhorst-Pack Kpoint mesh"//)')
       flag_gamma = fdf_boolean('Diag.GammaCentred',.false.)
       mp(1) = fdf_integer('Diag.MPMeshX',1)
       mp(2) = fdf_integer('Diag.MPMeshY',1)
       mp(3) = fdf_integer('Diag.MPMeshZ',1) 
       if(iprint_init>0.AND.inode==ionode) then
          if(flag_gamma) then
             write (io_lun,fmt='(/8x,a, i3," x ",i3," x ",i3," gamma-centred")') &
                  ' Monkhorst-Pack mesh: ', (mp(i), i=1,3)
          else
             write (io_lun,fmt='(/8x,a, i3," x ",i3," x ",i3)') &
                  ' Monkhorst-Pack mesh: ', (mp(i), i=1,3)
          end if
       else if(inode==ionode) then
          if(flag_gamma) then
             suffix = " G"
          else
             suffix = "  "
          end if
          write (io_lun,fmt='(4x,"Using a MP mesh for k-points: ", i3," x ",i3," x ",i3,a2)') &
               (mp(i), i=1,3), suffix
       end if
       if (mp(1) <= 0 .OR. mp(2) <= 0 .OR. mp(3) <= 0) &
            call cq_abort('K-points: number of k-points must be > 0!')
       nkp_tmp = mp(1)*mp(2)*mp(3)
       if(iprint_init>1.AND.inode==ionode) &
            write(io_lun,fmt='(8x,a, i4)') ' Number of k-points: ',nkp_tmp
       ! Read k-point shift, default (0.0 0.0 0.0)
       if(flag_gamma) then
          mp_shift = zero
          if(modulo(mp(1),2)==0) mp_shift(1) = half/mp(1)
          if(modulo(mp(2),2)==0) mp_shift(2) = half/mp(2)
          if(modulo(mp(3),2)==0) mp_shift(3) = half/mp(3)
       else
          mp_shift(1) = fdf_double('Diag.MPShiftX',zero)
          mp_shift(2) = fdf_double('Diag.MPShiftY',zero)
          mp_shift(3) = fdf_double('Diag.MPShiftZ',zero)
          if (mp_shift(1) >= one) then
             if(inode==ionode) write(io_lun,9)
             mp_shift(1) = mp_shift(1) - one
          end if
          if (mp_shift(2) >= one) then
             if(inode==ionode) write(io_lun,9)
             mp_shift(2) = mp_shift(2) - one
          end if
          if (mp_shift(3) >= one) then
             if(inode==ionode) write(io_lun,9)
             mp_shift(3) = mp_shift(3) - one
          end if
          if(iprint_init>0.AND.inode==ionode) &
               write (io_lun,fmt='(8x,a, 3f11.6)') &
               ' Monkhorst-Pack mesh shift:  ', &
               (mp_shift(i), i=1,3)
       end if
       ! Allocate
       allocate(kk_tmp(3,nkp_tmp),wtk_tmp(nkp_tmp),STAT=stat)
       if(stat/=0) &
            call cq_abort('FindEvals: couldnt allocate kpoints',nkp_tmp)
       call reg_alloc_mem(area_general,4*nkp_tmp,type_dbl)
       ! All k-points have weight 1 for now
       wtk_tmp(1:nkp_tmp) = one
       ! Generate fractional k-point coordinates plus shift
       ! Assume orthorhombic cell for now
       do i = 1, mp(1)       ! x axis 
          do j = 1, mp(2)    ! y axis
             do k = 1, mp(3) ! z axis
                counter = k + (j-1) * mp(3) + (i-1) * mp(2) * mp(3)
                kk_tmp(1,counter) = &
                     ((two * i - mp(1) - 1) / (2 * mp(1))) + mp_shift(1)
                kk_tmp(2,counter) = &
                     ((two * j - mp(2) - 1) / (2 * mp(2))) + mp_shift(2)
                kk_tmp(3,counter) = &
                     ((two * k - mp(3) - 1) / (2 * mp(3))) + mp_shift(3)
             end do
          end do
       end do
       ! Write out fractional k-points
       if(iprint_init>1.AND.inode==ionode) then
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
                   if (abs(kk_tmp(1,i) + kk_tmp(1,j))<RD_ERR .and. &
                        abs(kk_tmp(2,i) + kk_tmp(2,j))<RD_ERR .and. &
                        abs(kk_tmp(3,i) + kk_tmp(3,j))<RD_ERR) then
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
            call cq_abort('FindEvals: couldnt deallocate kpoints', nkp_tmp)
       if(iprint_init>1.AND.inode==ionode) then
          !
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
    ! Check polarisation
    if(flag_calc_pol) then
       if(nkp>1 .or. (nkp==1 .and. maxval(abs(kk))>RD_ERR)) &
            call cq_warn(sub_name, "Resta polarisation is only valid at gamma point")
    end if
    !
    ! BEGIN %%%% ASE printing %%%%
    !
    if (inode==ionode .and. write_ase ) then
       inquire(io_ase, opened=test_ase) 
       if ( .not. test_ase ) open(io_ase,file=ase_file, status='old', action='write',&
            iostat=stat, position='append')                       
       write(io_ase,*)
       write(io_ase,10) nkp
       do i=1,nkp
          write (io_ase,fmt='(8x,i5,3f15.6,f12.3)') &
               i, kk(1,i)*rcellx/(two*pi), kk(2,i)*rcelly/(two*pi), kk(3,i)*rcellz/(two*pi), wtk(i)
       end do
       close(io_ase)
       !
    end if
    !
    ! END %%%% ASE printing %%%%
    !    
    ! Write out k-points
    if(inode==ionode.AND.iprint_init>2) then
       write(io_lun,51) nkp
       do i=1,nkp
          write (io_lun,fmt='(8x,i5,3f15.6,f12.3)') &
               i,kk(1,i),kk(2,i),kk(3,i),wtk(i)
       end do
    end if
    ! Write out smearing temperature
    if(iprint_init>1.AND.inode==ionode) &
         write (io_lun,'(10x,"Value of kT used for smearing: ",f10.6, a2)') &
         en_conv * kT, en_units(energy_units)

    !****lat<$    
    call stop_backtrace(t=backtrace_timer,who='readDiagInfo')
    !****lat>$

    return

2   format(8x,'Block size (row, col): ',2i5)
3   format(8x,'Proc grid (row, col): ',2i5)
4   format(/8x,'***WARNING***',/,2x,&
         'No Kpoints block found - defaulting to Gamma point')
51  format(/8x,i4,' symmetry inequivalent Kpoints in Cartesian form (1/A): ')
7   format(8x,' All ',i4,' Kpoints in fractional coordinates: ')
9   format(/8x,'***WARNING***',/,8x,&
         'Specified Kpoint mesh shift in fractional coords >= 1.0.')
10  format(8x,i4,' symmetry inequivalent Kpoints in fractional coordinates: ')
11  format(8x, 'Number of processor groups for k-point parallelisation: ', i5)

  end subroutine readDiagInfo
  !!***

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
  !!   2020/10/30 lionel
  !!    correct for io_assign, wrong module!
  !!  SOURCE
  !!
  subroutine read_input_aux(aux)
    ! Module usage
    use global_module,   only: io_lun
    use auxiliary_types, only: group_aux
    use GenComms,        only: myid,cq_abort,gcopy
    use input_module,    only: io_assign,io_close
    use input_module,    only: fdf_block,fdf_endblock

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
  !!***
end module initial_read
