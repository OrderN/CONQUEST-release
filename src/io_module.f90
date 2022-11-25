! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module io_module
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/io_module
!!  NAME 
!!   io_module -- governs input and output of atom positions, types, 
!!     hopping and matrix writing
!!  PURPOSE
!!   Locate all input and output (bar comment writing) in one place so
!!   that MPI-type calls can be isolated from main body of code
!!  USES
!!   GenComms, basic_types, common, datatypes, global_module,
!!   group_module, matrix_module, species_module
!!  AUTHOR
!!    D.R.Bowler
!!  CREATION DATE   
!!   25/01/00 by D.R.Bowler
!!  MODIFICATION HISTORY 
!!   02/03/00 by DRB
!!    Changed the unit chosen to write out C-matrix
!!   12/04/00 by DRB
!!    write_c -> write_matrix
!!   18->19/04/00 by DRB
!!    Implementing the new generic group_set and primary_set types
!!   1/05/00 by DRB
!!    Incorporation into Conquest - removed the parameters read in 
!!    from stdin, and left it with reading make_prt.dat
!!   04/05/01 by DRB
!!    Used the ParaDens tidied output to improve readability
!!   21/06/2001 dave
!!    Added GenComms and RCS Id and Log tags and removed mpi
!!   18/03/2002 dave
!!    Added static tag for object file id
!!   24/06/2002 dave
!!    Changed read_atoms to use fdf io_assign routines (not raw
!!    open/close)
!!   14:00, 04/02/2003 drb 
!!    Added dump/grab routines for charge, matrices, blips, local and
!!    non-local pseudos for force checking
!!   16:43, 2003/04/15 dave
!!    Added read_atomic_coordinates and changed read_mult
!!   16:55, 2003/06/09 dave
!!    Further changes based on TM's work: now use atom_coord instead
!!    of various x_atom etc variables Changes also to mapping and
!!    numbering
!!   13:30, 22/09/2003 drb 
!!    Added constrained atoms
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote
!!    in terms of new matrix routines
!!   11:40, 10/11/2006 Veronika
!!    Added option for reading pdb files
!!   20/11/2006 Veronika
!!    Added option for writing out pdb files
!!   28/01/2008 Veronika
!!    Fixed writing out pdb files - the pdb-style atom names are now
!!    also written out
!!   2008/02/06 08:07 dave
!!    Changed for output to file not stdout
!!   2008/07/31 ast
!!    Added timers
!!   2008/09/01 08:20 dave
!!    Added new Conquest input routines
!!   2009/07/24 ast
!!    Filenames based on node/iteration number are now not limited to
!!    3 figures New routine get_file_name takes care of creating a
!!    sufficiently long string
!!   2015/06/25 17:18 dave
!!    Changed all file name lengths to 50 characters
!!   2016/03/15 14:02 dave
!!    Added tests for existence of files to all grab_ routines and abort if
!!    they don't exist
!!   2019/11/04 11:36 dave
!!    Removed redundant code (old SFC routines)
module io_module

  use datatypes,              only: double
  use numbers,                only: zero
  use global_module,          only: io_lun
  use GenComms,               only: cq_abort, gcopy
  use timer_module,           only: start_timer,     stop_timer,    cq_timer
  use timer_module,           only: start_backtrace, stop_backtrace
  use timer_stdclocks_module, only: tmr_std_matrices
  use input_module,           only: leqi, io_assign, io_close

  implicit none

  logical           :: pdb_format, append_coords, pdb_output 
  character(len=1)  :: pdb_altloc
  character(len=80) :: pdb_template

  !Maximum of wallclock time (in seconds): See subroutine 'check_stop' 2018.Jan.17 TM 
  real(double)      :: time_max =zero
  integer :: atom_output_threshold
   
  !Name and Format of  MatrixFile used in store_matrix
  logical          :: flag_MatrixFile_RankFromZero  ! Starting from 0 in ##### (*matrix2.i**.p#####)
  logical          :: flag_MatrixFile_BinaryFormat  ! Binary or Ascii
  logical          :: flag_MatrixFile_BinaryFormat_Grab ! Binary or Ascii for reading
  logical          :: flag_MatrixFile_BinaryFormat_Dump ! Binary or Ascii for writing
  ! flag_MatrixFile_BinaryFormat_Dump is set as flag_MatrixFile_BinaryFormat_Grab in the beginning,
  !  then, if we get a signal to finalise CONQUEST by calling "check_stop", it is set ast
  !  flag_MatrixFile_BinaryFormat_Dump_END.  
  !   (during the job, we need to keep the format of the file, prepared in the last job.)
  logical          :: flag_MatrixFile_BinaryFormat_Dump_END ! set from Conquest_input
  
  ! System signature
  ! Moved here from read_and_write so that it can be used for extended XYZ output
  ! Moved here from initial_read_module to slove the dependence problem
  character(len=80), save :: titles
  logical          :: flag_coords_xyz

!!***

contains

  ! --------------------------------------------------------------------
  ! Subroutine read_atomic_positions
  ! --------------------------------------------------------------------
  
  !!****f* io_module/read_atomic_positions *
  !!
  !!  NAME 
  !!   read_atomic_positions
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Reads positions, species and constraint type for atoms
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !!   datatypes, dimens, GenComms, global_module, species_module
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   16:26, 2003/04/15
  !!  MODIFICATION HISTORY
  !!   16:46, 2003/06/09 tm
  !!    id_glob_inv, atom_coord, species_glob is added.
  !!    and labelling in make_prt.dat is changed 
  !!   13:31, 22/09/2003 drb 
  !!    Added reading of flags to constrain atom movement in three directions
  !!   2006/10/09 08:21 dave
  !!    Tidying up output
  !!   20/11/2006 Veronika
  !!    Corrected reading of constraints from pdb file
  !!   15:05, 27/04/2007 drb 
  !!    Check to ensure right number of species added
  !!   2007/06/28, 21:37 mt + drb
  !!    Added coordinate wrapping to non-pdb coordinates.
  !!   2007/10/11 Veronika
  !!    Added coordinate wrapping to pdb coordinates
  !!    Added coordinate wrapping for coordinates outside 
  !!    the [-1*cell parameter; 2*cell parameter] range
  !!   2013/07/01 M.Arita & T.Miyazaki
  !!    Added shift_in_bohr when wrapping atoms and allocation of atom_coord_diff
  !!   2013/08/20 M.Arita
  !!    Bug fix & correct the if-statement
  !!   2015/06/08 lat
  !!    Added experimental backtrace
  !!   2018/01/22 tsuyoshi (with dave)
  !!    Allocate atom_coord_diff for all calculations
  !!   2019/04/04 14:17 dave
  !!    Correct bug in wrapping with non-fractional coordinates and Angstroms
  !!   2020/07/27 tsuyoshi
  !!    Added atom_vels  
  !!   2020/10/07 tsuyoshi
  !!    Removed allocation of atom_vels (moved to "control")
  !!  SOURCE
  !!
  subroutine read_atomic_positions(filename)

    use datatypes
    use dimens,         only: r_super_x, r_super_y, r_super_z
    use global_module,  only: x_atom_cell, y_atom_cell, z_atom_cell, &
                              ni_in_cell,                            &
                              flag_fractional_atomic_coords, rcellx, &
                              rcelly, rcellz, id_glob, iprint_init,  &
                              id_glob_inv, atom_coord, species_glob, &
                              flag_move_atom, area_init, shift_in_bohr, &
                              runtype,atom_coord_diff,id_glob_old,id_glob_inv_old
    use species_module, only: species, species_label, n_species
    use GenComms,       only: inode, ionode, cq_abort, cq_warn
    use memory_module,  only: reg_alloc_mem, type_dbl, type_int
    use units,          only: AngToBohr
    use units,          only: dist_units, ang
    use numbers,        only: RD_ERR, zero

    ! Passed variables
    character(len=*) :: filename

    ! Local variables
    type(cq_timer) :: backtrace_timer
    integer        :: lun, i, spec, stat
    logical        :: movex, movey, movez
    real(double)   :: x, y, z
    real(double), dimension(3) :: cell

    ! Local variables for reading pdb files
    integer      :: ios, j
    integer      :: num_move_atom
    real(double) :: num_move_atom_real
    real(double), dimension(3) :: angle
    character(len=80) :: pdb_line
    character(len=2)  :: atom_name

!****lat<$
    call start_backtrace(t=backtrace_timer,who='read_atomic_positions',where=9,level=2)
!****lat>$

    if(inode==ionode) then
       ! read a pdb file if General.pdb == .true.)
pdb:   if (pdb_format) then
          if(flag_fractional_atomic_coords) &
               call cq_abort('Pdb file format does not support &
                              &fractional coordinates.')
          if (iprint_init>2) &
               write(io_lun,fmt='(4x,a)') 'Entering read_atomic_positions, pdb file'
          call io_assign(lun)
          ! Go through the file and count the atoms
          open( unit=lun, file=filename, status='old', iostat=ios)
          if ( ios > 0 ) call cq_abort('Reading pdb file: file error')
          if (iprint_init > 2) &
               write (io_lun,'(5x,a)') &
                     'Counting atoms, checking for alternate locations'
          ni_in_cell = 0
first:    do
             read (lun, '(a)', iostat=ios) pdb_line
             ! End of file?
             if (ios < 0) exit first
             select case (pdb_line(1:6))
             case ('ATOM  ','HETATM')
                if (lgt( pdb_line(17:17),' ')) then 
                   ! File with alternate locations and no
                   ! General.PdbAltLoc specified?
                   !
                   ! This check fails if General.PdbAltLoc is in the
                   ! input file but has no value. In this case, ALL
                   ! entries with altloc will be ignored. User
                   ! beware. (Of course we can always sell it as a
                   ! feature.)
                   if (leqi( pdb_altloc,' ') ) then
                      call io_close(lun)
                      call cq_abort('PDB file with alternate &
                                     &locations, specify General.PdbAltLoc')
                   endif
                   ! Count the alternate location entry?
                   if (leqi( pdb_line(17:17), pdb_altloc)) &
                        ni_in_cell = ni_in_cell + 1
                else
                   ni_in_cell = ni_in_cell + 1
                end if
             end select
          end do first
          call io_close(lun)
          if (iprint_init>0) &
               write(io_lun,'(5x,a,i5)') 'Number of atoms: ', ni_in_cell

          ! Now read the file again and extracts the coordinates
          open(unit=lun, file=filename, status='old', iostat=ios)
          if (ios > 0) call cq_abort('Reading pdb file: file error')
          if (iprint_init > 2) write (*,*) 'Reading atomic positions'
          ! Allocate array for coordinates
          allocate(flag_move_atom(3,ni_in_cell),&
                   atom_coord(3,ni_in_cell),&
                   species_glob(ni_in_cell),STAT=stat)
          if(stat/=0) &
               call cq_abort("Failure to allocate coordinates: ",ni_in_cell)
          ! Are the following two lines necessary?
          call reg_alloc_mem(area_init, 3*ni_in_cell,type_dbl)
          call reg_alloc_mem(area_init, 4*ni_in_cell,type_int)

          i = 0
second:   do
             read (lun, '(a)', iostat=ios) pdb_line
             ! End of file?
             if (ios < 0) exit second
             select case (pdb_line(1:6))
             case ('ATOM  ','HETATM')
                if (lgt( pdb_line(17:17),' ')) then
                   ! Accept the alternate location?
                   if ( leqi( pdb_line(17:17), pdb_altloc)) then
                      i = i + 1
                      read (pdb_line(31:38),'(f8.3)') atom_coord(1,i)
                      read (pdb_line(39:46),'(f8.3)') atom_coord(2,i)
                      read (pdb_line(47:54),'(f8.3)') atom_coord(3,i)
                      atom_coord(1,i) = AngToBohr * atom_coord(1,i)
                      atom_coord(2,i) = AngToBohr * atom_coord(2,i)
                      atom_coord(3,i) = AngToBohr * atom_coord(3,i)
                      read (pdb_line(13:14),'(a2)') atom_name
                      atom_name = adjustl (atom_name)
                      do j = 1, n_species
                         if (leqi (atom_name,species_label(j) ) ) then
                            species_glob(i) = j
                            if(species_glob(i)>n_species) then
                               write(io_lun,&
                                     fmt='(4x,"** WARNING ! ** &
                                           &Species incompatibility between &
                                           &coordinates and input")')
                               call cq_abort("Species specified &
                                              &greater than number in input &
                                              &file: ",species_glob(i),&
                                              n_species)
                            end if
                            exit
                         end if
                      end do
                      flag_move_atom(:,i) = .false.
                      read (pdb_line(61:66),'(f6.0)') num_move_atom_real
                      num_move_atom = floor(num_move_atom_real)
                      if (num_move_atom_real - (real(num_move_atom)) > &
                          RD_ERR ) &
                           call cq_abort('The constraints on atom(s) &
                                          &have not been set. Check the pdb file.')
                      if (num_move_atom > 7) & 
                           call cq_abort('Out of range value for &
                                          &flag_move_atom. Check the pdb file.')
                      if (num_move_atom >= 4) then
                         flag_move_atom(3,i) = .true.
                         num_move_atom = num_move_atom - 4
                      end if
                      if (num_move_atom >= 2) then
                         flag_move_atom(2,i) = .true.
                         num_move_atom = num_move_atom - 2
                      end if
                      if (num_move_atom == 1) flag_move_atom(1,i) = .true.
                      !if (iprint_init > 0) &
                      !     write (io_lun,fmt='(3x, i7, 3f15.8, i3, 3L2)') &
                      !           i, atom_coord(1:3,i), species_glob(i), &
                      !           flag_move_atom(1:3,i)
                   end if
                else
                   i = i + 1
                   read (pdb_line(31:38),'(f8.3)') atom_coord(1,i)
                   read (pdb_line(39:46),'(f8.3)') atom_coord(2,i)
                   read (pdb_line(47:54),'(f8.3)') atom_coord(3,i)
                   atom_coord(1,i) = AngToBohr * atom_coord(1,i)
                   atom_coord(2,i) = AngToBohr * atom_coord(2,i)
                   atom_coord(3,i) = AngToBohr * atom_coord(3,i)
                   read (pdb_line(13:14),'(a2)') atom_name
                   atom_name = adjustl (atom_name)
                   do j = 1, n_species
                      if (leqi (atom_name,species_label(j) ) ) then
                         species_glob(i) = j
                         if(species_glob(i)>n_species) then
                            write (io_lun,&
                                   fmt='(2x,"** WARNING ! ** Species &
                                         &incompatibility between &
                                         &coordinates and input")')
                            call cq_abort("Species specified greater &
                                           &than number in input file: ",&
                                          species_glob(i),n_species)
                         end if
                         exit
                      end if
                   end do
                   flag_move_atom(:,i) = .false.
                   read (pdb_line(61:66),'(f6.0)') num_move_atom_real
                   num_move_atom = floor(num_move_atom_real)
                   if (num_move_atom_real - (real(num_move_atom)) > RD_ERR ) &
                        call cq_abort('The constraints on atom(s) &
                                       &have not been set. Check the pdb file.')
                   if (num_move_atom > 7) & 
                        call cq_abort('Out of range value for &
                                       &flag_move_atom. Check the pdb file.')
                   if (num_move_atom >= 4) then
                      flag_move_atom(3,i) = .true.
                      num_move_atom = num_move_atom - 4
                   end if
                   if (num_move_atom >= 2) then
                      flag_move_atom(2,i) = .true.
                      num_move_atom = num_move_atom - 2
                   end if
                   if (num_move_atom == 1) flag_move_atom(1,i) = .true.
                   !if (iprint_init > 0) &
                   !     write (io_lun,fmt='(3x, i7, 3f15.8, i3, 3L2)')&
                   !           i, atom_coord(1:3,i), species_glob(i), &
                   !           flag_move_atom(1:3,i)
                end if
             case ('CRYST1')
                read (pdb_line,'(6x,3f9.3,3f7.2,15x)') &
                     r_super_x, r_super_y, r_super_z, angle
                r_super_x = AngToBohr * r_super_x
                r_super_y = AngToBohr * r_super_y
                r_super_z = AngToBohr * r_super_z
                ! The following write statement is also in
                ! initial_read_module, but I think it's better to have
                ! it (also) here, because then the cell size will be
                ! printed together with the coordinates
                !write(io_lun,4) r_super_x, r_super_y, r_super_z
!4               format(/10x,'The simulation box has the following dimensions',/, &
!                       10x,'a = ',f9.5,' b = ',f9.5,' c = ',f9.5,' a.u.')
             end select
          end do second
          ! Wrap coordinates
          cell(1) = r_super_x
          cell(2) = r_super_y
          cell(3) = r_super_z
          do i = 1, ni_in_cell
             do j = 1, 3
                ! Introduce shift_in_bohr: (small shift for fractional coordinates may 
                ! cause a problem for very large systems.)
                !
                ! Originally, 'if'-statement in the next line was active, but since we need 
                ! a common and strict rule to treat the atoms on the boundary of the unit cell, 
                ! we need to do this wrapping with shift_in_bohr for all cases. 
                !
                ! Thus, I have commented out the next line.   23/Jul/2014 TM 
                !  (coordinate of the atoms on the boundary of the unit cell must be 0.000*** not 0.999***.)
                !
                !if ((atom_coord(j,i) < zero) .or. (atom_coord(j,i) > cell(j))) &
                   atom_coord(j,i) = &
                      atom_coord(j,i) - floor((atom_coord(j,i)+shift_in_bohr)/cell(j)) * cell(j)
             end do
          end do
          call io_close(lun)
       else ! Read normal coordinate file
          if(iprint_init>2) write(io_lun,'(10x,a40,a20)') 'Entering read_atomic_positions; reading ', filename
          call io_assign(lun)
          open(unit=lun,file=filename,status='old')
          ! Read supercell vector - for now it must be orthorhombic so
          ! we use x and y as dummy variables
          read(lun,*) r_super_x, x, y
          if(abs(x)>RD_ERR.OR.abs(y)>RD_ERR) call cq_warn('read_atomic_positions', &
               'Non-orthorhombic simulation cells are not supported by CONQUEST')
          read(lun,*) x,r_super_y, y
          if(abs(x)>RD_ERR.OR.abs(y)>RD_ERR) call cq_warn('read_atomic_positions', &
               'Non-orthorhombic simulation cells are not supported by CONQUEST')
          read(lun,*) x,y,r_super_z
          if(abs(x)>RD_ERR.OR.abs(y)>RD_ERR) call cq_warn('read_atomic_positions', &
               'Non-orthorhombic simulation cells are not supported by CONQUEST')
          read(lun,*) ni_in_cell
         !2010.06.25 TM (Angstrom Units in coords file, but not pdb)
          if(dist_units == ang) then
            r_super_x = r_super_x*AngToBohr
            r_super_y = r_super_y*AngToBohr
            r_super_z = r_super_z*AngToBohr
          endif
         !2010.06.25 TM (Angstrom Units in coords file, but not pdb)

          allocate(flag_move_atom(3,ni_in_cell),atom_coord(3,ni_in_cell),&
                   species_glob(ni_in_cell),STAT=stat)
          if(stat/=0) &
               call cq_abort("Failure to allocate coordinates: ",ni_in_cell)
          call reg_alloc_mem(area_init, 6*ni_in_cell,type_dbl)
          call reg_alloc_mem(area_init, 4*ni_in_cell,type_int)
          do i=1,ni_in_cell
             read(lun,*) x,y,z,species_glob(i),movex,movey,movez
             if(species_glob(i)>n_species) then
                write(io_lun,fmt='(6x,"** WARNING ! ** Species &
                                   &incompatibility between coordinates and &
                                   &input")')
                call cq_abort("Species specified greater than number &
                               &in input file: ", species_glob(i), n_species)
             end if
             if(flag_fractional_atomic_coords) then
                atom_coord(1,i) = x*r_super_x
                atom_coord(2,i) = y*r_super_y
                atom_coord(3,i) = z*r_super_z
             else
                atom_coord(1,i) = x
                atom_coord(2,i) = y
                atom_coord(3,i) = z
                !2010.06.25 TM (Angstrom Units in coords file, but not pdb)
                ! 2019/04/04 14:16 dave
                ! Moved here so that distances are corrected *before* wrapping below
                if(dist_units == ang) atom_coord(:,i)=atom_coord(:,i)*AngToBohr
             end if
             ! Wrap coordinates
             cell(1) = r_super_x
             cell(2) = r_super_y
             cell(3) = r_super_z
             do j = 1, 3
                ! Introduce shift_in_bohr: (small shift for fractional coordinates may
                ! cause a problem for very large systems.)
                !
                ! Originally, 'if'-statement in the next line was active, but since we need
                ! a common and strict rule to treat the atoms on the boundary of the unit cell, 
                ! we need to do this wrapping with shift_in_bohr for all cases. 
                !
                ! Thus, I have commented out the next line.   23/Jul/2014 TM   
                !  (coordinate of the atoms on the boundary of the unit cell must be 0.000*** not 0.999***.)
                !
                !if ((atom_coord(j,i) < zero) .or. (atom_coord(j,i) > cell(j))) &
                   atom_coord(j,i) = &
                      atom_coord(j,i) - floor((atom_coord(j,i)+shift_in_bohr)/cell(j)) * cell(j)
             end do
             flag_move_atom(1,i) = movex
             flag_move_atom(2,i) = movey
             flag_move_atom(3,i) = movez
          end do
          call io_close(lun)
       end if pdb
    end if
    if((iprint_init>0) .or. (iprint_init==0.AND.ni_in_cell<atom_output_threshold)) &
         call print_atomic_positions
    call gcopy(ni_in_cell)
    if(inode/=ionode) &
         allocate(flag_move_atom(3,ni_in_cell), atom_coord(3,ni_in_cell),&
                  species_glob(ni_in_cell), STAT=stat)
    if(stat/=0) call cq_abort("Failure to allocate coordinates: ",ni_in_cell)
    allocate(id_glob(ni_in_cell), id_glob_inv(ni_in_cell),&
             x_atom_cell(ni_in_cell), y_atom_cell(ni_in_cell),&
             z_atom_cell(ni_in_cell),species(ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Failure to allocate coordinates(2): ",ni_in_cell)
    call gcopy(r_super_x)
    call gcopy(r_super_y)
    call gcopy(r_super_z)
    call gcopy(species_glob,ni_in_cell)
    call gcopy(atom_coord, 3, ni_in_cell)
    call gcopy(flag_move_atom, 3, ni_in_cell)
    rcellx = r_super_x
    rcelly = r_super_y
    rcellz = r_super_z
    allocate(atom_coord_diff(3,ni_in_cell), STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating atom_coord_diff: ', 3, ni_in_cell)
    allocate(id_glob_old(ni_in_cell),id_glob_inv_old(ni_in_cell), STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating id_glob_old/id_glob_inv_old: ', &
         ni_in_cell)
    atom_coord_diff=zero
    !call gcopy(atom_coord_diff,3,ni_in_cell)
    id_glob_old=0
    id_glob_inv_old=0

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='read_atomic_positions')
!****lat>$

    return
  end subroutine read_atomic_positions
  !!***

!!$  ! --------------------------------------------------------------------
!!$  ! Subroutine write_atomic_positions_ase
!!$  ! --------------------------------------------------------------------
!!$  
!!$  !!****f* io_module/write_atomic_positions *
!!$  !!
!!$  !!  NAME 
!!$  !!   write_atomic_positions
!!$  !!  USAGE
!!$  !! 
!!$  !!  PURPOSE
!!$  !!   writes positions, species and constraint type for atoms
!!$  !!  INPUTS
!!$  !! 
!!$  !! 
!!$  !!  USES
!!$  !!   datatypes, dimens, GenComms, global_module, species_module
!!$  !!  AUTHOR
!!$  !!   D.R.Bowler
!!$  !!  CREATION DATE
!!$  !!   2022/11/08 L. Truflandier
!!$  !!     Copied from write_atomic_positions to print in ASE output
!!$  !      Forced to be fractional coordinates
!!$  !!  MODIFICATION HISTORY
!!$  !!  SOURCE
!!$  !!
!!$  subroutine write_atomic_positions_ase()
!!$    use datatypes
!!$    use numbers,        only: zero
!!$    use dimens,         only: r_super_x, r_super_y, r_super_z
!!$    use global_module,  only: ni_in_cell, species_glob, flag_move_atom, &
!!$                              min_layer, io_ase, ase_file, iprint_init, &
!!$                              IPRINT_TIME_THRES3, atom_coord, nspin
!!$    use species_module, only: species, species_label, n_species
!!$    use ScalapackFormat,only: matrix_size
!!$    use DiagModule,             only: nkp
!!$    use GenComms,       only: inode, ionode, cq_abort
!!$    use timer_module
!!$
!!$    ! Passed variables
!!$
!!$    ! Local variables
!!$    integer :: i, counter, stat
!!$    type(cq_timer)             :: tmr_l_tmp1
!!$
!!$    if(inode==ionode) then
!!$       
!!$       call start_timer(tmr_l_tmp1,WITH_LEVEL)       
!!$       if (iprint_init + min_layer > 2) write(io_lun,fmt='(6x,a)') 'Writing atomic positions in ASE output'
!!$
!!$       open(io_ase,file=ase_file, status='old', action='readwrite', iostat=stat, position='rewind')
!!$
!!$       !open(io_ase,file=ase_file, status='old', action='readwrite', iostat=stat, position='append')
!!$             
!!$       if (stat .ne. 0) call cq_abort('Error opening file !')
!!$       !
!!$       if ( nspin == 2 ) then
!!$          counter = nkp*3 + nspin*nkp + nspin*nkp*(matrix_size/3) + 1 + 2
!!$          if ( mod(matrix_size,3) > 0 ) counter = counter + nspin*nkp
!!$       
!!$       else
!!$          counter = nkp*3 + nkp*(matrix_size/3) + 1 + 1
!!$          if ( mod(matrix_size,3) > 0 ) counter = counter + nkp
!!$       
!!$       end if
!!$       counter = counter + 7 + n_species + 2 + nkp + 2 + 2 + ni_in_cell + 2 + 1
!!$
!!$       !%%%%%%%%%%
!!$       !counter = 0
!!$       !%%%%%%%%%%
!!$       do i = 1, counter
!!$          read (io_ase,*)
!!$       end do
!!$       
!!$       ! Read supercell vector - for now it must be orthorhombic so we use x and y as dummy variables
!!$       write(io_ase,fmt='(3f25.17)') r_super_x, zero, zero
!!$       write(io_ase,fmt='(3f25.17)') zero, r_super_y, zero
!!$       write(io_ase,fmt='(3f25.17)') zero, zero, r_super_z
!!$       write(io_ase,fmt='(i8)') ni_in_cell
!!$       do i=1,ni_in_cell
!!$          write(io_ase,fmt='(3f25.17,i6,3L3)') atom_coord(1,i)/r_super_x,atom_coord(2,i)/r_super_y,&
!!$               atom_coord(3,i)/r_super_z, species_glob(i),flag_move_atom(1,i),flag_move_atom(2,i), &
!!$               flag_move_atom(3,i)
!!$       end do
!!$       call stop_print_timer(tmr_l_tmp1,"writing atomic positions ASE",IPRINT_TIME_THRES3)
!!$    end if
!!$    
!!$    return
!!$  end subroutine write_atomic_positions_ase
  
  ! --------------------------------------------------------------------
  ! Subroutine write_atomic_positions
  ! --------------------------------------------------------------------
  
  !!****f* io_module/write_atomic_positions *
  !!
  !!  NAME 
  !!   write_atomic_positions
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   writes positions, species and constraint type for atoms
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !!   datatypes, dimens, GenComms, global_module, species_module
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   2006/11/07 08:02 dave
  !!  MODIFICATION HISTORY
  !!   20/11/2006 Veronika
  !!   Added option for writing out pdb files
  !!   2008/07/18
  !!     Added timers
  !!   2013/07/01 M.Arita & T.Miyazaki
  !!   Changed formats from f18.10 to f25.17
  !!  SOURCE
  !!
  subroutine write_atomic_positions(filename, pdb_temp)

    use datatypes
    use numbers,        only: zero
    use dimens,         only: r_super_x, r_super_y, r_super_z
    use global_module,  only: x_atom_cell, y_atom_cell, z_atom_cell,   &
                              ni_in_cell,                              &
                              flag_fractional_atomic_coords, rcellx,   &
                              rcelly, rcellz, iprint_init, atom_coord, &
                              species_glob, flag_move_atom, area_init, &
                              IPRINT_TIME_THRES3, min_layer
    use species_module, only: species, species_label
    use GenComms,       only: inode, ionode, cq_abort
    use units,          only: BohrToAng
    use timer_module

    ! Passed variables
    character(len=*) :: filename, pdb_temp

    ! Local variables
    integer                    :: lun, i, spec, stat, template, ios, j
    real(double)               :: x, y, z, beta
    logical                    :: movex, movey, movez
    character(len=80)          :: pdb_line
    character(len=2)           :: atom_name
    real(double), dimension(3) :: coords
    type(cq_timer)             :: tmr_l_tmp1

    if(inode==ionode) then
       call start_timer(tmr_l_tmp1,WITH_LEVEL)
       if (pdb_output) then
          if (iprint_init + min_layer > 2) write(io_lun,fmt='(6x,a)') 'Writing atomic positions'
          call io_assign(lun)
          ! No appending of coords for a pdb file
          open (unit = lun, file = 'output.pdb', status = 'replace', iostat = ios)
          if ( ios > 0 ) call cq_abort('Cannot create file.')
          ! Open the template file
          call io_assign(template)
          open (unit = template, file = pdb_temp, status = 'old', iostat = ios)
          if ( ios > 0 ) then
             write(io_lun,*) 'Filename was: |',pdb_temp,'|'
             call cq_abort('Reading template file: file error ',ios)
          end if
          ! Read the template file, reprint all lines that are not ATOM, HETATM, or CRYST1
          ! to the output file
          ! If the line contains coords/cell parameters, replace them by the calculated ones.
          ! NB We do not update the connectivity information (the CONECT entries).
          i = 0
          do
            read (template, '(a)', iostat=ios) pdb_line
            ! End of file?
            if (ios < 0) exit
            select case (pdb_line(1:6))
            case ('ATOM  ','HETATM')
              if ( leqi( pdb_line(17:17), pdb_altloc) .or. leqi( pdb_line(17:17), ' ')) then
                i = i + 1
                coords(1) = BohrToAng * atom_coord(1,i)
                coords(2) = BohrToAng * atom_coord(2,i)
                coords(3) = BohrToAng * atom_coord(3,i)
                atom_name = adjustr (species_label(species_glob(i))(1:2))
                beta = 0.0
                if (flag_move_atom(1,i)) beta = beta + 1
                if (flag_move_atom(2,i)) beta = beta + 2
                if (flag_move_atom(3,i)) beta = beta + 4
                write (lun,'(a12,a2,a16,3f8.3,a6,f6.2,a14)') &
                      pdb_line(1:12), atom_name, pdb_line(15:30), &
                      coords(:), pdb_line(55:60), beta, pdb_line(67:80)
              else
                ! If the entry is ATOM or HETATM but the alternate location
                ! does not match, just copy the line from the template
                write (lun,'(a)') pdb_line
              endif
            case ('CRYST1')
              coords(1) = BohrToAng * r_super_x
              coords(2) = BohrToAng * r_super_y
              coords(3) = BohrToAng * r_super_z
              ! Once we get to non-orthorhombic cells the line below has to be updated
              write (lun,'(a6,3f9.3,3f7.2,a)') pdb_line(1:6), coords(:), &
                    90.0, 90.0, 90.0, pdb_line(56:80)
            case default
              write (lun,'(a)') pdb_line
            end select
          end do
          call io_close(lun)
          call io_close(template)
       else
          if (iprint_init + min_layer > 2) write(io_lun,fmt='(6x,a)') 'Writing atomic positions'
          call io_assign(lun)
          if(append_coords) then
             open(unit=lun,file=filename,position='append')
             write(lun,*)
          else
             open(unit=lun,file=filename)
          end if
          ! Read supercell vector - for now it must be orthorhombic so we use x and y as dummy variables
          write(lun,fmt='(3f25.17)') r_super_x, zero, zero
          write(lun,fmt='(3f25.17)') zero, r_super_y, zero
          write(lun,fmt='(3f25.17)') zero, zero, r_super_z
          write(lun,fmt='(i8)') ni_in_cell
          do i=1,ni_in_cell
             if(flag_fractional_atomic_coords) then
                write(lun,fmt='(3f25.17,i6,3L3)') atom_coord(1,i)/r_super_x,atom_coord(2,i)/r_super_y,&
                     atom_coord(3,i)/r_super_z, species_glob(i),flag_move_atom(1,i),flag_move_atom(2,i), &
                     flag_move_atom(3,i)
             else
                write(lun,fmt='(3f25.17,i6,3L3)') atom_coord(1,i),atom_coord(2,i),atom_coord(3,i), &
                     species_glob(i),flag_move_atom(1,i),flag_move_atom(2,i), &
                     flag_move_atom(3,i)
             end if
          end do
          call io_close(lun)
       end if ! pdb_output
       call stop_print_timer(tmr_l_tmp1,"writing atomic positions",IPRINT_TIME_THRES3)
    end if
    return
  end subroutine write_atomic_positions
  !!***

  ! --------------------------------------------------------------------
  ! Subroutine read_mult
  ! --------------------------------------------------------------------
  
  !!****f* io_module/read_mult *
  !!
  !!  NAME 
  !!   read_mult -- reads in data for multiplications
  !!  USAGE
  !!   read_mult(Processor id, partition data structure)
  !!   read_mult(myid,parts)
  !!  PURPOSE
  !!   Reads in data for partitions and distributes to all processors
  !!  INPUTS
  !!   integer :: myid - processor id
  !!   type(group_set) :: parts - contains info on partitions
  !!  USES
  !!   datatypes, global_module, group_module, basic_types
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !! 
  !!  MODIFICATION HISTORY
  !!   21/06/2001 dave
  !!    Added cq_abort and gcopy
  !!   2006/10/19 16:41 dave
  !!    Added filename for partition reading
  !!   2013/04/03 L.Tong
  !!    Swapped the original create_sfc_partitions with the new
  !!    sfc_partitions_to_processors from sfc_partitions_module
  !!   2013/08/20 M.Arita
  !!    Correct the if-statement 
  !!   2017/08/29 jack baker & dave
  !!    Removed r_super_x references (redundant)
  !!  SOURCE
  !!
  subroutine read_mult(myid,parts,part_file)

    use datatypes
    use global_module,    only: numprocs, iprint_init, id_glob,       &
                                ni_in_cell, x_atom_cell, y_atom_cell, &
                                z_atom_cell, atom_coord, id_glob_inv, &
                                species_glob, runtype, id_glob_old, id_glob_inv_old
    use group_module,     only: make_cc2, part_method, PYTHON, HILBERT
    use basic_types
    use species_module,   only: species
    use maxima_module,    only: maxpartsproc, maxatomspart,           &
                                maxatomsproc, maxpartscell
    use construct_module,      only: init_group
    use sfc_partitions_module, only: sfc_partitions_to_processors

    implicit none

    ! Passed variables
    integer           :: myid
    type(group_set)   :: parts
    character(len=80) :: part_file

    integer           :: np_in_cell, mx_tmp_edge
    integer           :: irc,ierr,nnd,nnd1,np,np1,ind_part
    integer           :: n_cont,n_beg,ni,ni1,i
    integer           :: id_global
    integer           :: ind_file, stat

    integer           :: isendbuf(6)
    real(double)      :: rsendbuf(6)
    type(cq_timer)    :: backtrace_timer

!****lat<$
    call start_backtrace(t=backtrace_timer,who='read_mult',where=9,level=2)
!****lat>$

    !parts%ng_on_node = 0
    if(part_method == PYTHON) then
       if(myid==0) then
          if(iprint_init>1) then
             write(io_lun,'(10x,a29)')   '   Reading partition data    '
             write(io_lun,'(10x,a29/)') '-----------------------------'
             write(io_lun,'(10x,a29)')   '*** Partition Data Starts ***'
          end if
          call read_partitions(parts,part_file)
          if(iprint_init>1) write(io_lun,'(10x,a29)') '***  Partition Data Ends  ***'
       end if
       isendbuf = 0
       if(myid==0) then
          isendbuf(1) = parts%ngcellx
          isendbuf(2) = parts%ngcelly
          isendbuf(3) = parts%ngcellz
          isendbuf(4) = maxpartsproc
          isendbuf(5) = maxatomspart
          isendbuf(6) = maxatomsproc
          np_in_cell  = parts%ngcellx*parts%ngcelly*parts%ngcellz
       endif
       call gcopy(isendbuf,6)
       if(myid/=0) then
          parts%ngcellx = isendbuf(1)
          parts%ngcelly = isendbuf(2)
          parts%ngcellz = isendbuf(3)
          maxpartsproc  = isendbuf(4)
          maxatomspart  = isendbuf(5)
          maxatomsproc  = isendbuf(6)
          np_in_cell    = parts%ngcellx*parts%ngcelly*parts%ngcellz
          maxpartscell  = np_in_cell
          mx_tmp_edge   = max(parts%ngcellx,parts%ngcelly,parts%ngcellz)
          call init_group(parts, maxpartsproc, mx_tmp_edge, &
                          np_in_cell, maxatomspart, numprocs)
       endif
       call gcopy(parts%ng_on_node,numprocs)
       call gcopy(parts%inode_beg,numprocs)
       call gcopy(parts%ngnode,parts%mx_gcell)
       call gcopy(parts%nm_group,parts%mx_gcell)
       call gcopy(parts%icell_beg,parts%mx_gcell)
       call gcopy(id_glob, ni_in_cell)
       call gcopy(id_glob_inv, ni_in_cell)
       !if (leqi(runtype,'md') .OR. leqi(runtype,'cg')) then
       if (.NOT. leqi(runtype,'static')) then
         call gcopy(id_glob_old,ni_in_cell)
         call gcopy(id_glob_inv_old,ni_in_cell)
       endif
    else if (part_method == HILBERT) then
       if (iprint_init > 1.AND.myid==0) then
          write(io_lun,'(10x,a33)' ) 'Partitioning using Hilbert curves'
          write(io_lun,'(10x,a33/)') '---------------------------------'
       end if
       !call create_sfc_partitions(myid, parts)
       call sfc_partitions_to_processors(parts)
       np_in_cell = parts%ngcellx*parts%ngcelly*parts%ngcellz
       if (iprint_init > 2.AND.myid==0) write(io_lun,fmt='(10x,a)') 'Finished partitioning'
    end if
    ! inverse table to npnode
    do np=1,np_in_cell
       parts%inv_ngnode(parts%ngnode(np))=np
    end do
    call make_cc2(parts,numprocs)
    do ni = 1, ni_in_cell
       id_global= id_glob(ni)
       x_atom_cell(ni) = atom_coord(1,id_global)
       y_atom_cell(ni) = atom_coord(2,id_global)
       z_atom_cell(ni) = atom_coord(3,id_global)
       species(ni)     = species_glob(id_global)
    end do
!****lat<$
    call stop_backtrace(t=backtrace_timer,who='read_mult')
!****lat>$
    
    return
  end subroutine read_mult
  !!***

  ! --------------------------------------------------------------------
  ! Subroutine read_partitions
  ! --------------------------------------------------------------------
  
  !!****f* io_module/read_partitions *
  !!
  !!  NAME 
  !!   read_atoms -- reads in partition data
  !!  USAGE
  !!   read_partitions(Partition data structure, number of processors)
  !!  PURPOSE
  !!   Reads in all data associated with partitions  
  !!   (numbers within partition) on master processor and distributes to 
  !!     all other processors
  !!  INPUTS
  !!   type(group_set) :: parts - partition data structure
  !!   integer :: nnode - number of processors
  !!  USES
  !!   datatypes, global_module, basic_types, common, ham
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   16:50, 2003/04/15 dave
  !!  MODIFICATION HISTORY
  !!   2006/10/19 16:42 dave
  !!    Added filename for partitions (defaults to make_prt.dat)
  !!   2007/06/28 21:16 mt + drb
  !!    Added check for partition file access error.
  !!   2011/11/17 10:19 dave
  !!    Bug fix to format for 1141
  !!   2013/08/2013
  !!    Added storage of id_glob_old & id_glob_inv_old
  !!  SOURCE
  !!
  subroutine read_partitions(parts, part_file)

    ! Module usage
    use datatypes
    use global_module,    only: numprocs, iprint_init, ni_in_cell, &
                                id_glob, id_glob_inv, id_glob_old, &
                                id_glob_inv_old, runtype
    use maxima_module,    only: maxpartsproc, maxatomspart,        &
                                maxatomsproc, maxpartscell
    use basic_types
    use construct_module, only: init_group
    use GenComms,         only: cq_abort, myid

    implicit none

    ! Passed variables
    type(group_set)   :: parts
    character(len=80) :: part_file

    ! Local variables
    type(cq_timer) :: backtrace_timer
    integer :: nnode
    integer :: nnd,nnd1,np,np1,ind_part,n_cont,n_beg,ni,ni1
    integer :: irc,ierr,np_in_cell, lun, glob_count, map, ios
    integer :: ind_global, ntmpx, ntmpy, ntmpz, ntmp1, ntmp2, ntmp3, mx_tmp_edge

!****lat<$
    call start_backtrace(t=backtrace_timer,who='read_partitions',where=1,level=4)
!****lat>$  

    if(iprint_init>2.AND.myid==0) write(io_lun,fmt='(4x,a)') 'Entering read_partitions'
    call io_assign(lun)

    open(unit=lun, file=part_file, status='old', iostat=ios)
    if ( ios > 0 ) call cq_abort('Reading partition file: file error')
    !open(unit=lun,file=part_file)
    ! Read and check partitions along cell sides
    read(lun,*) ntmpx,ntmpy,ntmpz
    mx_tmp_edge = max(ntmpx,ntmpy,ntmpz)
    ! Find and check total number of partitions
    np_in_cell=ntmpx*ntmpy*ntmpz
    if(np_in_cell<=0) &
         call cq_abort('read_partitions: no partitions ! ',np_in_cell)
    ! Read and check number of processors
    read(lun,*) nnode
    if((nnode<1).or.(nnode>numprocs)) &
         call cq_abort('read_partitions: nodes wrong ',nnode,numprocs)
    ! Loop over nodes, reading bundle information
    maxpartsproc = 0
    maxatomspart = 0
    maxatomsproc = 0
    do nnd=1,nnode
       ntmp3 = 0  ! count atoms on processor
       ! First node number, partitions and starting pointer
       read(lun,*) nnd1,ntmp1, ntmp2
       if(nnd1/=nnd) call cq_abort('read_partitions: node label wrong ',nnd1,nnd)
       ! Check for correct numbers of partitions
       if(ntmp1<0) call cq_abort('read_partitions: parts on proc ', ntmp1)
       if(ntmp1>maxpartsproc) maxpartsproc = ntmp1
       ! If there are partitions, read them in
       if(ntmp1>0) then
          ! Check for correct maxima
          if(ntmp2 + ntmp1 -1>np_in_cell) &
               call cq_abort('read_partitions: too many parts in cell ',&
                             ntmp1 + ntmp2-1, np_in_cell)
          do np=1,ntmp1
             read(lun,*) np1,ind_part,n_cont,n_beg
             if(np1/=np) &
                  call cq_abort('read_partitions: part label error ',np1,np)
             if(n_cont>maxatomspart) maxatomspart = n_cont
             ntmp3 = ntmp3 + n_cont
             if(n_cont>0) then
                do ni=1,n_cont
                   read(lun,*) ni1,map
                   if(ni1/=ni) &
                        call cq_abort('read_partitions: atom label wrong ',ni1,ni)
                enddo
             endif
          enddo ! Loop over np = parts%ng_on_node
       endif ! if(parts%ng_on_node>0)
       if(ntmp3>maxatomsproc) maxatomsproc = ntmp3
    enddo ! nnd = 1,nnode
    if(iprint_init>3.AND.myid==0) &
         write(io_lun,fmt='(6x,a,2i5)') 'Atoms proc max: ',maxatomsproc, maxpartsproc
    call init_group(parts, maxpartsproc, mx_tmp_edge, np_in_cell, &
                    maxatomspart, numprocs)
    maxpartscell = np_in_cell
    close(unit=lun)

    open(unit=lun, file=part_file, status='old', iostat=ios)
    if ( ios > 0 ) call cq_abort('Reading partition file: file error')
    !open(unit=lun,file=part_file)
    ! Read and check partitions along cell sides
    read(lun,*) parts%ngcellx,parts%ngcelly,parts%ngcellz
    if(iprint_init>0.AND.myid==0) &
         write(io_lun,102) parts%ngcellx,parts%ngcelly,parts%ngcellz
    if(max(parts%ngcellx,parts%ngcelly,parts%ngcellz)>parts%mx_gedge) then
       call cq_abort('read_partitions: too many parts edge ', &
            max(parts%ngcellx,parts%ngcelly,parts%ngcellz),parts%mx_gedge)
    endif
    ! Read and check number of processors
    read(lun,*) nnode
    if(iprint_init>0.AND.myid==0) write(io_lun,110) nnode
    if((nnode<1).or.(nnode>numprocs)) then
       call cq_abort('read_partitions: nodes wrong ',nnode,numprocs)
    endif
    ! Output header for atoms
!    if(iprint_init==0) write(io_lun,*) '   Partn        x              y              z        Species   Proc'
    ! Loop over nodes, reading bundle information
    glob_count = 0
    do nnd=1,nnode
       ! First node number, partitions and starting pointer
       read(lun,*) nnd1,parts%ng_on_node(nnd),parts%inode_beg(nnd)
       if(iprint_init>0.AND.myid==0) &
            write(io_lun,111) nnd1,parts%ng_on_node(nnd),parts%inode_beg(nnd)
!       if(iprint_init==1) write(io_lun,*) '   Partn        x              y              z        Species   Proc'
       if(nnd1/=nnd) then
          call cq_abort('read_partitions: node label wrong ',nnd1,nnd)
       endif
       ! Check for correct numbers of partitions
       if((parts%ng_on_node(nnd)<0).or.&
            (parts%ng_on_node(nnd)>parts%mx_ngonn)) then
          call cq_abort('read_partitions: too many parts on proc ',&
                        parts%ng_on_node(nnd),parts%mx_ngonn)
       endif
       ! If there are partitions, read them in
       if(parts%ng_on_node(nnd).gt.0) then
          ! Check for correct maxima
          if(parts%inode_beg(nnd)+parts%ng_on_node(nnd)-1>parts%mx_gcell) then
             call cq_abort('read_partitions: too many parts in cell ',&
                           parts%inode_beg(nnd) + parts%ng_on_node(nnd) - 1, &
                           parts%mx_gcell)
          endif
          ! Loop over partitions on processor
          do np=1,parts%ng_on_node(nnd)
             read(lun,*) np1,ind_part,n_cont,n_beg
             ! Check for maxima
             if(ind_part>parts%mx_gcell) then
                call cq_abort('read_partitions: bad part label ', &
                     ind_part,parts%mx_gcell)
             endif
             ! Check for ordering of partitions
             if(np1/=np) then
                call cq_abort('read_partitions: part label error ',np1,np)
             endif
             ! Assign and write out
             parts%ngnode(parts%inode_beg(nnd)+np-1)=ind_part
             parts%nm_group(ind_part)=n_cont
             parts%icell_beg(ind_part)=n_beg
             if(iprint_init>1.AND.myid==0) &
                  write(io_lun,113) ind_part,np1,n_cont,n_beg
             ! Read atoms
             if(n_cont>0) then
                ! Check maxima
                if(parts%icell_beg(ind_part)+n_cont-1>ni_in_cell) &
                     call cq_abort('read_partitions: too many atoms ', &
                                   parts%icell_beg(ind_part)+n_cont-1,ni_in_cell)
                ! Read and write out
                do ni=1,n_cont
                   ! Define global number of this atom
                   ! (TM) Is it global number?
                   glob_count = glob_count+1
                   ! Read atom in partition, map from coordinate file to partition-based number
                   ! comment by TM,  ni1: partition labelling, map: global labelling
                   !                 glob_count : CC labelling??
                   read(lun,*) ni1,map
                   if(ni1/=ni) &
                        call cq_abort('read_partitions: atom label wrong ',ni1,ni)
                   ! Store global number
                   ! id_glob(parts%icell_beg(ind_part)+ni-1) = glob_count
                   id_glob(parts%icell_beg(ind_part)+ni-1) = map
                   ! Map from coord file to partition number
                   id_glob_inv(map) = parts%icell_beg(ind_part)+ni-1
                   !write(io_lun,*) ' Global map: ',parts%icell_beg(ind_part)+ni-1,map
                enddo
             endif
          enddo ! Loop over np = parts%ng_on_node
       endif ! if(parts%ng_on_node>0)
    enddo ! nnd = 1,nnode
    call io_close(lun)
    if (.NOT. leqi(runtype,'static')) then
      id_glob_old = id_glob
      id_glob_inv_old = id_glob_inv
    endif

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='read_partitions')
!****lat>$

    return
    ! Format statements placed here to improve readability of code
101 format(/8x,'cell lengths:',3f12.6)
102 format(/8x,'partitions along cell sides',3i5)
110 format(/8x,'no. of processors (i.e. no. of partition bundles):',i5/)
111 format(/8x,'processor no:',i5/8x,'no. of partitions:',i5,&
         ' initial partition:',i5)
112 format(/8x,'partitions on current processor:')
113 format(/9x,'Partn no: ',i5,' Local no: ',i5,2x,'contains ',i5,&
         ' atoms, starting at: ',i5/)
114 format(10x,i5,2x,i5,2x,3e15.6,i5)
1141 format(11x,i4,2x,3e15.6,4x,i2,3x,i6)

  end subroutine read_partitions
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine read_blocks
  ! --------------------------------------------------------------------
  
  !!****f* io_module/read_blocks *
  !!
  !!  NAME 
  !!   read_atoms -- reads in block data
  !!  USAGE
  !!   read_blocks(Block data structure, number of processors)
  !!  PURPOSE
  !!   Reads in all data associated with blocks  
  !!   (numbers within block) on master processor and distributes to 
  !!     all other processors
  !!  INPUTS
  !!   type(group_set) :: blocks - block data structure
  !!   integer :: nnode - number of processors
  !!  USES
  !!   datatypes, global_module, basic_types, common, ham
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   12:03, 15/11/2004 dave 
  !!  MODIFICATION HISTORY
  !!   2008/09/11 08:09 dave
  !!    Rewritten to do I/O on one processor only
  !!  SOURCE
  !!
  subroutine read_blocks(blocks)

    ! Module usage
    use datatypes
    use global_module, only: iprint_init, numprocs
    use maxima_module, only: maxngrid, maxblocks
    use basic_types, only: group_set
    use construct_module, only: init_group
    use group_module,  only: make_cc2
    use block_module, only : in_block_x, in_block_y, in_block_z, n_pts_in_block, &
         nx_in_block, ny_in_block, nz_in_block
    use dimens, only: n_grid_x, n_grid_y, n_grid_z
    use GenComms, only: inode, ionode
    use input_module, only: fdf_string

    implicit none

    ! Passed variables
    type(group_set) :: blocks

    ! Local variables
    integer :: nnd,nnd1,np,np1,ind_part,n_cont,n_beg,ni,ni1, nnode
    integer :: irc,ierr,np_in_cell, lun, glob_count, map
    integer :: ind_global, ntmpx, ntmpy, ntmpz, nmyblocks,ntmp1, ntmp2
    integer:: n_block_x, n_block_y, n_block_z, maxtmp
    character(len=80) :: blk_coord_file, def

    if(inode==ionode.AND.iprint_init>2) write(io_lun,fmt='(4x,a,3i6)') 'Entering read_blocks ', in_block_x,in_block_y,in_block_z
    !--- find numbers of blocks in each direction
    n_block_x = n_grid_x/in_block_x
    n_block_y = n_grid_y/in_block_y
    n_block_z = n_grid_z/in_block_z
    if(n_block_x*in_block_x/=n_grid_x) call cq_abort("Block and grid sizes incompatible ! ",in_block_x,n_grid_x)
    if(n_block_y*in_block_y/=n_grid_y) call cq_abort("Block and grid sizes incompatible ! ",in_block_y,n_grid_y)
    if(n_block_z*in_block_z/=n_grid_z) call cq_abort("Block and grid sizes incompatible ! ",in_block_z,n_grid_z)
    blocks%ngcellx=n_block_x
    blocks%ngcelly=n_block_y
    blocks%ngcellz=n_block_z
    nx_in_block=in_block_x
    ny_in_block=in_block_y
    nz_in_block=in_block_z
    if(inode==ionode) then
       call io_assign(lun)
       def = 'make_blk.dat'
       blk_coord_file = fdf_string(80,'IO.Blocks',def)
       open(unit=lun,file=blk_coord_file)
       ! Read and check blocks along cell sides
       read(lun,*) ntmpx, ntmpy, ntmpz
       if(ntmpx/=n_block_x.OR.ntmpy/=n_block_y.OR.ntmpz/=n_block_z) then
          write(io_lun,fmt='(4x,"In input, ReadBlocks T has been set, so the distribution of blocks is ")')
          write(io_lun,fmt='(4x,"read from a file.  There is an error specifying numbers of blocks.  ")')
          write(io_lun,fmt='(4x,"Found from the file: ",3i6)') ntmpx, ntmpy, ntmpz
          write(io_lun,fmt='(4x,"Calculated from input using grid points: ",3i6)') n_block_x, n_block_y, n_block_z
          call cq_abort("Aborting ! Please bring block input file and overall input file together.")
       end if
       ! Find and check total number of blocks
       np_in_cell=ntmpx* ntmpy* ntmpz
       if(np_in_cell<=0) call cq_abort('read_blocks: no blocks in cell ',np_in_cell)
       ! Read and check number of processors
       read(lun,*) nnode
       if((nnode<1).or.(nnode>numprocs)) then
          call cq_abort('read_blocks: nodes wrong ',nnode,numprocs)
       endif
       ! Read in block information to find maxima
       maxblocks = 0
       do nnd=1,nnode
          ! First node number, blocks and starting pointer
          read(lun,*) nnd1,ntmp1,ntmp2 
          if(nnd1/=nnd) call cq_abort('read_blocks: node label wrong ',nnd1,nnd)
          ! Check for correct numbers of blocks
          if(ntmp1<0) call cq_abort('read_blocks: too many blocks on proc ', nmyblocks)
          if(ntmp1>maxblocks) maxblocks = ntmp1
          ! If there are blocks, read them in
          if(ntmp1>0) then
             do np=1,ntmp1
                read(lun,*) np1,ind_part!,n_cont,n_beg
             enddo ! Loop over np = blocks%ng_on_node
          endif ! if(blocks%ng_on_node>0)
       enddo ! nnd = 1,nnode
       close(unit=lun)
    end if
    ! Now initialise block derived type
    call gcopy(ntmpx)
    call gcopy(ntmpy)
    call gcopy(ntmpz)
    ntmp1=max(ntmpx,ntmpy,ntmpz)
    np_in_cell=ntmpx* ntmpy* ntmpz
    call gcopy(maxblocks)
    call init_group(blocks, maxblocks, ntmp1, np_in_cell, n_pts_in_block, numprocs)
    if(inode==ionode) then
       open(unit=lun,file=blk_coord_file)
       ! Read and check blocks along cell sides
       read(lun,*) blocks%ngcellx,blocks%ngcelly,blocks%ngcellz
       !if(iprint_init>0) write(io_lun,102) blocks%ngcellx,blocks%ngcelly,blocks%ngcellz
       ! Find and check total number of blocks
       np_in_cell=blocks%ngcellx*blocks%ngcelly*blocks%ngcellz
       ! Read and check number of processors
       read(lun,*) nnode
       !if(iprint_init>=0) write(io_lun,110) nnode
       ! Loop over nodes, reading bundle information
       glob_count = 0
       do nnd=1,nnode
          ! First node number, blocks and starting pointer
          read(lun,*) nnd1,blocks%ng_on_node(nnd),blocks%inode_beg(nnd)
          !if(iprint_init>0) write(io_lun,111) nnd1,blocks%ng_on_node(nnd),blocks%inode_beg(nnd)
          if(nnd1/=nnd) then
             call cq_abort('read_blocks: node label wrong ',nnd1,nnd)
          endif
          ! Check for correct numbers of blocks
          if((blocks%ng_on_node(nnd)<0).or.&
               (blocks%ng_on_node(nnd)>maxblocks)) then
             call cq_abort('read_blocks: too many blocks on proc ',&
                  blocks%ng_on_node(nnd),maxblocks)
          endif
          ! If there are blocks, read them in
          if(blocks%ng_on_node(nnd)>0) then
             ! Check for correct maxima
             if(blocks%inode_beg(nnd)+blocks%ng_on_node(nnd)-1>np_in_cell) then
                call cq_abort('read_blocks: too many blocks in cell ',&
                     blocks%inode_beg(nnd)+blocks%ng_on_node(nnd)-1, np_in_cell)
             endif
             ! Loop over blocks on processor
             do np=1,blocks%ng_on_node(nnd)
                read(lun,*) np1,ind_part!,n_cont,n_beg
                ! Check for maxima
                if(ind_part>np_in_cell) then
                   call cq_abort('read_blocks: bad block label ', &
                        ind_part,np_in_cell)
                endif
                ! Check for ordering of blocks
                if(np1/=np) then
                   call cq_abort('read_blocks: block label error ',np1,np)
                endif
                ! Assign and write out
                blocks%ngnode(blocks%inode_beg(nnd)+np-1)=ind_part
                ! These next two lines are essentially meaningless, but added for completeness
                blocks%nm_group(ind_part)=n_pts_in_block!n_cont
                blocks%icell_beg(ind_part)=(ind_part-1)*n_pts_in_block+1!n_beg
                !if(iprint_init>1) write(io_lun,113) ind_part,np1,n_cont,n_beg
             enddo ! Loop over np = blocks%ng_on_node
          endif ! if(blocks%ng_on_node>0)
       enddo ! nnd = 1,nnode
    end if
    call gcopy(blocks%ngcellx)
    call gcopy(blocks%ngcelly)
    call gcopy(blocks%ngcellz)
    call gcopy(blocks%ng_on_node,numprocs)
    call gcopy(blocks%inode_beg,numprocs)
    call gcopy(blocks%ngnode,np_in_cell)
    if(inode/=ionode) then
       do nnd=1,numprocs
          do np=1,blocks%ng_on_node(nnd)
             ind_part = blocks%ngnode(blocks%inode_beg(nnd)+np-1)
             blocks%nm_group(ind_part)=n_pts_in_block!n_cont
             blocks%icell_beg(ind_part)=(ind_part-1)*n_pts_in_block+1!n_beg
          end do
       end do
    end if
    ! Define maximum number of grid points on any processor for FFTs
    maxngrid = n_pts_in_block * maxblocks
    ! FFT check
    maxtmp = n_grid_x * (n_grid_y*n_grid_z/numprocs + 1)
    if(maxtmp>maxngrid) maxngrid = maxtmp
    maxtmp = n_grid_y * (n_grid_z*n_grid_x/numprocs + 1)
    if(maxtmp>maxngrid) maxngrid = maxtmp
    maxtmp = n_grid_z * (n_grid_x*n_grid_y/numprocs + 1)
    if(maxtmp>maxngrid) maxngrid = maxtmp
    ! Build inverse map
    do np=1,np_in_cell
       blocks%inv_ngnode(blocks%ngnode(np))=np
    enddo
    call make_cc2(blocks,numprocs)
    call io_close(lun)
    ! Format statements placed here to improve readability of code
101 format(/2x,'cell lengths:',3f12.6)
102 format(/2x,'blocks along cell sides',3i5)
110 format(/'no. of processors (i.e. no. of block bundles):',i5/)
111 format(/'processor no:',i5/'no. of blocks:',i5,&
         ' initial block:',i5)
112 format(/'blocks on current processor:')
113 format(/1x,'Block no: ',i5,' Local no: ',i5,2x,'contains ',i5,&
         ' atoms, starting at: ',i5/)
114 format(2x,i5,2x,i5,2x,3e15.6,i5)
1141 format(3x,i4,2x,3e15.6,4x,i2,3x,i6)
    return
  end subroutine read_blocks
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine dump_charge
  ! --------------------------------------------------------------------
  
  !!****f* io_module/dump_charge *
  !!
  !!  NAME 
  !!   dump_charge
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Dumps the charge density
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   10:58, 15/10/2002 drb 
  !!  MODIFICATION HISTORY
  !!   2011/04/13 L.Tong
  !!     Added optional flag controlling whether to dump total charge
  !!     density (spin = 0) or the up/down components (spin = 1 and 2
  !!     respectively)
  !!   2011/12/19 L.Tong
  !!     Changed filename length from 10 to 20, to accommodate the added
  !!     length of spin dependent file names
  !!  SOURCE
  !!
  subroutine dump_charge(density,size,inode,spin)

    use datatypes
    use block_module, only: n_pts_in_block
    use primary_module, only: domain
    use global_module, only: numprocs
    !use set_blipgrid_module, only: naba_atoms_of_blocks, supp

    ! Passed variables
    integer :: size, inode
    real(double), dimension(size) :: density
    integer, optional :: spin

    ! Local variables
    integer :: lun, block, n_point, n_i, n
    character(len=50) :: filename

    ! Build a filename based on node number
    if (present(spin)) then
       select case (spin)
       case (0)
          call get_file_name('chden', numprocs, inode, filename)
       case (1)
          call get_file_name('chden_up', numprocs, inode, filename)
       case (2)
          call get_file_name('chden_dn', numprocs, inode, filename)
       case default
          call cq_abort("dump_charge: unrecagonised value for spin &
               &output flag, must be 0, 1, or 2, but is ", spin)
       end select
    else
       call get_file_name('chden',numprocs,inode,filename)
    end if
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename)
    ! Dump charge density
    !do block=1, domain%groups_on_node
    !   n_point = (block - 1) * n_pts_in_block
    !   do n_i=1, naba_atoms_of_blocks(supp)%no_of_atom(block)*NSF*n_pts_in_block, n_pts_in_block
    !      do n=1, n_pts_in_block
    !         write(unit=lun,fmt='(f30.15)') density(n_point+n)
    !      end do
    !   end do
    !end do
    do n=1, domain%groups_on_node*n_pts_in_block
       write(unit=lun,fmt='(g13.6)') density(n)
    end do
    call io_close(lun)
    return
  end subroutine dump_charge
  !!***


  !!****f* io_module/dump_band_charge
  !! PURPOSE
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   David Bowler
  !! CREATION DATE 
  !!   2015/07/02 08:22 dave
  !! MODIFICATION HISTORY
  !!   2015/07/02 08:23 dave
  !!    This was dump_charge2 (which served no purpose that I could see)
  !! SOURCE
  !!
  subroutine dump_band_charge(stub, density, size, inode, spin,kp,energy,num_kpts)

    use datatypes
    use block_module, only: n_pts_in_block
    use primary_module, only: domain
    use global_module, only: numprocs
    !use set_blipgrid_module, only: naba_atoms_of_blocks, supp

    ! Passed variables
    integer :: size, inode
    real(double), dimension(size) :: density
    character(len=*) :: stub
    integer :: spin 
    real(double), optional :: energy
    real(double), optional, dimension(3) :: kp
    integer, optional :: num_kpts

    ! Local variables
    integer :: lun, block, n_point, n_i, n
    character(len=50) :: filename

    ! Build a filename based on node number
    select case (spin)
    case (0)
       call get_file_name (stub//'den', numprocs, inode, filename)
    case (1)
       call get_file_name (stub//'den_up', numprocs, inode, filename)
    case (2)
       call get_file_name (stub//'den_dn', numprocs, inode, filename)
    end select
    ! Open file
    call io_assign (lun)
    open (unit = lun, file = filename, position = 'append')
    ! We need to write a header of the k-point and energy if split-down
    if(present(num_kpts)) then
       write(lun,fmt='(2i10)') num_kpts, domain%groups_on_node * n_pts_in_block
    end if
    if(present(kp).AND.present(energy)) then
       write(lun,fmt='(3f12.5,f17.11)') kp(1:3),energy
    end if
    do n=1, domain%groups_on_node * n_pts_in_block
       write (unit=lun, fmt='(g13.6)') density(n)
    end do
    call io_close (lun)
    return
  end subroutine dump_band_charge
  !!***

  !!****f* io_module/write_eigenvalues
  !! PURPOSE
  !!  Writes out eigenvalues
  !! INPUTS
  !!  eval
  !!  n_evals
  !!  nkp
  !!  nspin
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   David Bowler
  !! CREATION DATE 
  !!   2015/07/09 08:16
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine write_eigenvalues(eval,n_evals,nkp,nspin,kk,wtk,Ef)

    use datatypes
    
    implicit none

    ! Passed variables
    integer :: n_evals,nkp,nspin
    real(double), dimension(n_evals,nkp,nspin) :: eval
    real(double), dimension(3,nkp) :: kk
    real(double), dimension(nkp) :: wtk
    real(double), dimension(nspin) :: Ef

    ! Local variables
    integer :: lun,sp,kp,ev
    
    call io_assign (lun)
    open (unit = lun, file = 'eigenvalues.dat')
    write(lun,fmt='("#  ",i7," eigenvalues ",i4," kpoints")') n_evals, nkp
    if(nspin==1) then
       write(lun,fmt='("# Ef: ",f18.10)') Ef(1)
    else
       write(lun,fmt='("# Ef: ",2f18.10)') Ef(1),Ef(2)
    end if
    write(lun,fmt='("# Format: nk kx ky kz weight, followed by eigenvalues")')
    do sp = 1,nspin
       do kp = 1,nkp
          write(lun,fmt='(i4,4f12.5)') kp,kk(1,kp),kk(2,kp),kk(3,kp),wtk(kp)
          do ev = 1,n_evals
             write(lun,fmt='(i6,f18.10)') ev,eval(ev,kp,sp)
          end do
       end do
    end do
    call io_close(lun)
    return
  end subroutine write_eigenvalues
  !!***

  !!****f* io_module/write_eigenvalues_format_ase
  !! PURPOSE
  !!  Writes out eigenvalues with ASE format
  !! INPUTS
  !!  eval
  !!  n_evals
  !!  nkp
  !!  nspin
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   David Bowler/Lionel Truflandier
  !! CREATION DATE 
  !!   2022/10/29 16:00
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine write_eigenvalues_format_ase(eval,occ,n_evals,nkp,nspin,kk,Ef,io,file,skip_lines)

    use datatypes
    use units

    implicit none

    ! Passed variables
    integer,   intent(in) :: n_evals, nkp, nspin
    integer,   intent(in) :: io, skip_lines
    character(len=80), intent(in) :: file

    real(double), dimension(n_evals,nkp,nspin) :: eval
    real(double), dimension(n_evals,nkp,nspin) :: occ
    real(double), dimension(3,nkp)             :: kk
    real(double), dimension(2)                 :: Ef

    ! Local variables
    real(double),dimension(2) :: BandE
    integer :: i, j, spin, stat, counter, nlines


    !open(io,file=file, status='old', action='read', iostat=stat, position='rewind')
    !if (stat .ne. 0) call cq_abort('Error opening file !')
    !do
    !   read (io,*, END=11)
    !   nlines = nlines + 1
    !end do
    !11  call io_close(io)
    !print*, nlines
    
    open(io,file=file, status='old', action='readwrite', iostat=stat, position='rewind')
    if (stat .ne. 0) call cq_abort('ASE/eigenvalues error opening file !')
    !
    do i = 1, skip_lines
       read (io,*)
    end do
    !
    !if ( nspin == 2 ) then
    !   counter = nkp*3 + nspin*nkp + nspin*nkp*(n_evals/3) + 1
    !   if ( mod(n_evals,3) > 0 ) counter = counter + nspin*nkp
    
    !else
    !   counter = nkp*3 + nkp*(n_evals/3) + 1
    !   if ( mod(n_evals,3) > 0 ) counter = counter + nkp
       
    !end if

    !if ( nlines > skip_lines ) then
    !   do i = 1, counter - skip_lines
    !      backspace(io)
    !   end do
    !end if
    !
    ! Write bands!
    !
    write(io,*)
    !
    bandE = zero
    !
    do i = 1, nkp
       write (io, 7) i, kk(1,i), kk(2,i), kk(3,i)
       do spin = 1, nspin
          if (nspin == 2) &
               write (io, '(10x,"For spin = ",i1)') spin
          do j = 1, n_evals, 3
             
             if (j == n_evals) then
                write (io, 8) eval(j,i,spin), occ(j,i,spin)

                bandE(spin) = bandE(spin) + eval(j,i,spin) * occ(j,i,spin)
                
             else if (j == n_evals - 1) then
                write (io, 9) eval(j,i,spin), occ(j,i,spin), &
                     eval(j+1,i,spin), occ(j+1,i,spin)
                
                bandE(spin) = bandE(spin) + eval(j,i,spin) * occ(j,i,spin) + &
                     eval(j+1,i,spin) * occ(j+1,i,spin)

             else
                write (io, 10) eval(j,i,spin), occ(j,i,spin), &
                     eval(j+1,i,spin), occ(j+1,i,spin), &
                     eval(j+2,i,spin), occ(j+2,i,spin)
                
                bandE(spin) = bandE(spin) + eval(j,i,spin) * occ(j,i,spin) + &
                     eval(j+1,i,spin) * occ(j+1,i,spin) + &
                     eval(j+2,i,spin) * occ(j+2,i,spin)
             endif
          end do ! j=matrix_size
          write (io, &
               fmt='(10x,"Sum of eigenvalues for spin = ", &
               &i1, ": ", f18.11," ", a2)') &
               spin, en_conv * bandE(spin), en_units(energy_units)
       end do ! spin
       if (nspin == 2) then
          write (io, &
               fmt='(10x,"Total sum of eigenvalues: ", f18.11, " ",a2)') &
               en_conv * (bandE(1) + bandE(2)), en_units(energy_units)
       else
          write(io, 4) en_conv * two * bandE(1), en_units(energy_units)
       end if
    end do ! do i = 1, nkp
    !
    ! Write Fermi energy
    !
    do spin = 1, nspin
       write (io, 13) spin, en_conv * Ef(spin), &
            en_units(energy_units)
    end do
    !
    call io_close(io)
    !    
4   format(10x,'Sum of eigenvalues: ',f18.11,' ',a2)
7   format(10x,'Eigenvalues and occupancies for k-point ',i3,' : ',3f12.5)
8   format(10x,f12.5,f6.3,2x)
9   format(10x,f12.5,f6.3,2x,f12.5,f6.3,2x)
10  format(10x,f12.5,f6.3,2x,f12.5,f6.3,2x,f12.5,f6.3,2x)
13  format(10x,'Fermi energy for spin = ',i1,' is ',f18.11,' ',a2)
    
    return

  end subroutine write_eigenvalues_format_ase
  !!***
  
  !!****f* io_module/dump_DOS
  !! PURPOSE
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   David Bowler
  !! CREATION DATE 
  !!   2015/07/02 10:46
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine dump_DOS(DOS,Ef)

    use datatypes
    use global_module, only: n_DOS, E_DOS_max, E_DOS_min, sigma_DOS, nspin

    ! Passed variables
    real(double), dimension(n_DOS,nspin) :: DOS
    real(double), dimension(nspin) :: Ef

    ! Local variables
    integer :: lun, i
    real(double) :: dE, thisE
    character(len=50) :: filename

    filename = 'DOSout.dat'
    ! Open file
    call io_assign (lun)
    open (unit = lun, file = filename)
    ! Header
    if(nspin==1) then
       write(lun,fmt='(2x,"# Ef: ",f12.5," DOS limits: ",2f12.5," n DOS: ",i5)') &
            Ef(nspin),E_DOS_min,E_DOS_max, n_DOS
    else if(nspin==2) then
       write(lun,fmt='(2x,"# Ef(up.down) : ",2f12.5," DOS limits: ",2f12.5," n DOS: ",i5)') &
            Ef(1:nspin),E_DOS_min,E_DOS_max, n_DOS
    end if
    write(lun,fmt='(2x,"# Broadening: ",f12.5)') sigma_DOS
    dE = (E_DOS_max - E_DOS_min)/real(n_DOS - 1,double)
    thisE = E_DOS_min
    if(nspin==1) then
       do i=1,n_DOS
          write(lun,fmt='(2f17.10)') thisE,DOS(i,1)
          thisE = thisE + dE
       end do
    else if(nspin==2) then
       do i=1,n_DOS
          write(lun,fmt='(3f17.10)') thisE,DOS(i,1),-DOS(i,2)
          thisE = thisE + dE
       end do
    end if
    call io_close (lun)
    return
  end subroutine dump_DOS
  !!***

  
  !!****f* io_module/dump_projected_DOS
  !! PURPOSE
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   David Bowler
  !! CREATION DATE 
  !!   2015/07/02 17:21
  !! MODIFICATION HISTORY
  !!   2017/11/09 16:00 nakata
  !!    Added atomf to find whether the PDOS is with PAOs or MSSFs
  !!   2018/09/19 18:00 nakata
  !!    Corrected the order of dimension of pDOS
  !!   2018/09/19 18:30 nakata
  !!    Added orbital angular momentum resolved DOS (pDOS_angmom)
  !!   2018/10/22 14:22 dave & jsb
  !!    Adding (l,m)-projected DOS
  !! SOURCE
  !!
  subroutine dump_projected_DOS(pDOS,Ef,pDOS_angmom,Nangmom)

    use datatypes
    use global_module, only: n_DOS, E_DOS_max, E_DOS_min, sigma_DOS, flag_pDOS_angmom, &
                             nspin, atomf, sf, ni_in_cell, flag_pDOS_lm
    use primary_module,  only: bundle

    ! Passed variables
    real(double), dimension(n_DOS,bundle%n_prim,nspin) :: pDOS
    real(double), dimension(nspin) :: Ef
    real(double), OPTIONAL, dimension(:,:,:,:,:) :: pDOS_angmom ! bin, atom, l, m, spin
    integer, OPTIONAL :: Nangmom

    ! Local variables
    integer :: lun, i, j, k, iprim, natom, tmp_col
    real(double) :: dE, thisE
    character(len=50) :: filename, fmt_DOS
    character(len=200) :: colstr

    if (atomf==sf) natom = bundle%n_prim
    if (atomf/=sf) natom = ni_in_cell

    do iprim = 1,natom
       if (atomf==sf) write(filename,'("Atom",I0.7,"DOS.dat")') bundle%ig_prim(iprim)
       if (atomf/=sf) write(filename,'("Atom",I0.7,"DOS.dat")') iprim ! iprim is equal to global ID if atomf = paof (MSSFs)
       ! Open file
       call io_assign (lun)
       open (unit = lun, file = filename)
       ! Header
       if(nspin==1) then
          write(lun,fmt='(2x,"# Ef: ",f12.5," DOS limits: ",2f12.5," n DOS: ",i5)') &
               Ef(nspin),E_DOS_min,E_DOS_max, n_DOS
       else if(nspin==2) then
          write(lun,fmt='(2x,"# Ef(up.down) : ",2f12.5," DOS limits: ",2f12.5," n DOS: ",i5)') &
               Ef(1:nspin),E_DOS_min,E_DOS_max, n_DOS
       end if
       write(lun,fmt='(2x,"# Broadening: ",f12.5)') sigma_DOS
       ! Indicate columns
       if(flag_PDOS_angmom) then
          if(flag_pDOS_lm) then
             colstr = "  #     Energy     |   Total pDOS   |    l=0, m= 0   |    l=1, m=-1   |&
                  &    l=1, m= 0   |    l=1, m= 1   |    l=2, m=-2   |    l=2, m=-1   |&
                  &    l=2, m= 0   |    l=2, m= 1   |    l=2, m= 2   "
          else
             colstr = "  #     Energy     |   Total pDOS   |      l=0       |      l=1       |&
                  &      l=2       "
          end if
       else
          colstr =    "  #     Energy     |   Total pDOS    "
       end if
       write(lun,fmt='(a)') colstr
       dE = (E_DOS_max - E_DOS_min)/real(n_DOS - 1,double)
       thisE = E_DOS_min
       if(nspin==1) then
          if (flag_PDOS_angmom) then
             if(flag_pDOS_lm) then
                tmp_col = 0
                do i=1,Nangmom
                   tmp_col = tmp_col + 2*i-1  ! Because i is l+1 
                end do
                write(fmt_DOS,*) tmp_col+2 ! NB this is number of columns in format
                fmt_DOS = '('//trim(adjustl(fmt_DOS))//'f17.10)'
                do i=1,n_DOS
                   write(lun,fmt_DOS) thisE,pDOS(i,iprim,1),((pDOS_angmom(i,iprim,j,k,1),k=1,2*j-1),j=1,Nangmom)
                   thisE = thisE + dE
                end do
             else
                write(fmt_DOS,*) Nangmom+2 ! NB this is number of columns in format
                fmt_DOS = '('//trim(adjustl(fmt_DOS))//'f17.10)'
                do i=1,n_DOS
                   write(lun,fmt_DOS) thisE,pDOS(i,iprim,1),(pDOS_angmom(i,iprim,j,1,1),j=1,Nangmom)
                   thisE = thisE + dE
                end do
             end if
          else
             do i=1,n_DOS
                write(lun,fmt='(2f17.10)') thisE,pDOS(i,iprim,1)
                thisE = thisE + dE
             end do
          endif
       else if(nspin==2) then
          if (flag_PDOS_angmom) then
             if(flag_pDOS_lm) then
                tmp_col = 0
                do i=1,Nangmom
                   tmp_col = tmp_col + 2*i-1  ! Because i is l+1 
                end do
                write(fmt_DOS,*) tmp_col*2+3
                fmt_DOS = '('//trim(adjustl(fmt_DOS))//'f17.10)'
                do i=1,n_DOS
                   write(lun,fmt_DOS) thisE,pDOS(i,iprim,1),((pDOS_angmom(i,iprim,j,k,1),k=1,2*j-1),j=1,Nangmom), &
                        -pDOS(i,iprim,2),((-pDOS_angmom(i,iprim,j,k,2),k=1,2*j-1),j=1,Nangmom)
                   thisE = thisE + dE
                end do
             else
                write(fmt_DOS,*) Nangmom*2+3 ! Total number of columns
                fmt_DOS = '('//trim(adjustl(fmt_DOS))//'f17.10)'
                do i=1,n_DOS
                   write(lun,fmt_DOS) thisE,pDOS(i,iprim,1),(pDOS_angmom(i,iprim,j,1,1),j=1,Nangmom), &
                        -pDOS(i,iprim,2),(-pDOS_angmom(i,iprim,j,1,2),j=1,Nangmom)
                   thisE = thisE + dE
                end do
             end if
          else
             do i=1,n_DOS
                write(lun,fmt='(3f17.10)') thisE,pDOS(i,iprim,1),-pDOS(i,iprim,2)
                thisE = thisE + dE
             end do
          end if
       end if
       call io_close (lun)
    end do
    return
  end subroutine dump_projected_DOS
  !!***

  

  ! --------------------------------------------------------------------
  ! Subroutine grab_charge
  ! --------------------------------------------------------------------
  
  !!****f* io_module/grab_charge *
  !!
  !!  NAME 
  !!   grab_charge
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Grabs the charge density (dumped by dump_charge)
  !!   - Note that this subroutine simply graps the list data from a
  !!     file and stores in the array density. Any post processing of
  !!     the data should be done after the call and outside of this
  !!     subroutine.
  !!  INPUTS
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   10:58, 15/10/2002 drb 
  !!  MODIFICATION HISTORY
  !!   2011/11/27 L.Tong
  !!   - Added spin polarisation
  !!   2013/03/06 17:10 dave
  !!   - file name length increased
  !!  SOURCE
  !!
  subroutine grab_charge(density, size, inode, spin)

    use datatypes
    use block_module,   only: n_pts_in_block
    use primary_module, only: domain
    use global_module,  only: numprocs
    !use set_blipgrid_module, only: naba_atoms_of_blocks, supp

    ! Passed variables
    integer :: size, inode
    real(double), dimension(size) :: density
    integer, optional :: spin

    ! Local variables
    integer           :: lun, block, n_point, n_i, n, ios
    character(len=50) :: filename

    ! Build a filename based on node number
    if (present(spin)) then
       select case (spin)
       case (0)
          call get_file_name('chden', numprocs, inode, filename)
       case (1)
          call get_file_name('chden_up', numprocs, inode, filename)
       case (2)
          call get_file_name('chden_dn', numprocs, inode, filename)
       end select
    else
       call get_file_name('chden', numprocs, inode, filename)
    end if
    ! Open file
    call io_assign(lun)
    open(unit=lun, file=filename,status='old',iostat=ios)
    if(ios /= 0) call cq_abort('grab_charge: failed to open input file '//filename)
    ! Grab charge density
    ! do block = 1, domain%groups_on_node
    !   n_point = (block - 1) * n_pts_in_block
    !   do n_i = 1, naba_atoms_of_blocks(supp)%no_of_atom(block) * NSF * n_pts_in_block, &
    !            n_pts_in_block
    !      do n=1, n_pts_in_block
    !         read(unit=lun,fmt='(f30.15)') density(n_point+n)
    !      end do
    !   end do
    ! end do
    do n = 1, domain%groups_on_node * n_pts_in_block
       read (unit=lun, fmt='(g13.6)') density(n)
    end do
    call io_close(lun)
    return
  end subroutine grab_charge
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine dump_matrix
  ! --------------------------------------------------------------------
  
  !!****f* io_module/dump_matrix *
  !!
  !!  NAME 
  !!   dump_matrix
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Dumps a matrix
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   16:12, 15/10/2002 drb 
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine dump_matrix(stub,matA,inode)

    use datatypes
    use mult_module, only: return_matrix_value_pos, return_matrix_len,&
         matrix_index, return_matrix_value
    use primary_module, only: bundle
    use group_module, only: parts
    use matrix_data, only: mat
    use cover_module, only: BCS_parts
    use global_module, only: id_glob, numprocs

    ! Passed variables
    integer :: inode, matA
    character(len=*) :: stub

    ! Local variables
    integer :: lun, element, nf1, nf2,len
    integer :: np, ni, iprim,i, jsf, neigh, ist, gcspart,isf, Ah, globno
    character(len=50) :: filename

    ! Build a filename based on node number
    call get_file_name(stub//'matrix',numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename)
    !Ah = matrix_index(matA)
    !iprim = 0
    !do np=1,bundle%groups_on_node
    !   if(bundle%nm_nodgroup(np) > 0) then
    !      do ni=1,bundle%nm_nodgroup(np)
    !         iprim=iprim+1
    !         i=bundle%ig_prim(iprim)
    !         do jsf=1,mat(np,Ah)%ndimi(ni)
    !            do neigh = 1,mat(np,Ah)%n_nab(ni)
    !               ist = mat(np,Ah)%i_acc(ni)+neigh-1
    !               gcspart = BCS_parts%icover_ibeg(mat(np,Ah)%i_part(ist))+mat(np,Ah)%i_seq(ist)-1
    !               globno = id_glob( parts%icell_beg(BCS_parts%lab_cell(mat(np,Ah)%i_part(ist) )) +mat(np,Ah)%i_seq(ist)-1 )
    !               do isf=1,mat(np,Ah)%ndimj(ist)
    !                  write(unit=lun,fmt='(2i4,g14.6)') i,globno,return_matrix_value(matA,np,ni,iprim,neigh,jsf,isf)
    !               end do
    !            end do
    !         end do
    !      end do
    !   end if
    !end do
    len = return_matrix_len(matA)
    ! Dump matrix
    do element=1,len
       write(unit=lun,fmt='(f25.18)') return_matrix_value_pos(matA,element)
    end do
    call io_close(lun)
    return
  end subroutine dump_matrix
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine grab_matrix
  ! --------------------------------------------------------------------
  
  !!****f* io_module/grab_matrix *
  !!
  !!  NAME 
  !!   grab_matrix
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Grabs a matrix (dumped by dump_matrix)
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   16:12, 15/10/2002 drb 
  !!  MODIFICATION HISTORY
  !!   2010/11/30 15:39 dave
  !!    Fixed formatting mismatch between dump_matrix and grab_matrix
  !!   2013/03/06 17:10 dave
  !!   - file name length increased
  !!  SOURCE
  !!
  subroutine grab_matrix(stub,matA,inode)

    use datatypes
    use mult_module, only: store_matrix_value_pos, return_matrix_len
    use global_module, only: numprocs

    ! Passed variables
    integer :: inode, matA
    character(len=*) :: stub

    ! Local variables
    integer :: lun, element, nf1, nf2,len, ios
    real(double) :: val
    character(len=50) :: filename

    ! Build a filename based on node number
    call get_file_name(stub//'matrix',numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename,status='old',iostat=ios)
    if(ios /= 0) call cq_abort('grab_matrix: failed to open input file '//filename)
    len = return_matrix_len(matA)
    ! Grab matrix
    call start_timer(tmr_std_matrices)
    do element=1,len
       read(unit=lun,fmt='(f25.18)') val!matrix(nf2,nf1,element)
       call store_matrix_value_pos(matA,element,val)
    end do
    call stop_timer(tmr_std_matrices)
    call io_close(lun)
    return
  end subroutine grab_matrix
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine dump_blips
  ! --------------------------------------------------------------------
  
  !!****f* io_module/dump_blips *
  !!
  !!  NAME 
  !!   dump_blips
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Dumps a blips
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   16:12, 15/10/2002 drb 
  !!  MODIFICATION HISTORY
  !!   2013/03/06 17:10 dave
  !!   - file name length increased
  !!
  !!  SOURCE
  !!
  subroutine dump_blips(stub,support,inode)

    use datatypes
    use functions_on_grid, only: gridfunctions, fn_on_grid
    use global_module, only: numprocs

    implicit none

    ! Passed variables
    integer :: inode, support
    character(len=*) :: stub

    ! Local variables
    integer :: lun, element, nf1, nf2
    character(len=50) :: filename

    ! Build a filename based on node number
    call get_file_name(stub//'grid',numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename)
    ! Dump blips
    do element=1,gridfunctions(support)%size
       write(unit=lun,fmt='(f30.15)') gridfunctions(support)%griddata(element)
    end do
    call io_close(lun)
    return
  end subroutine dump_blips
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine grab_blips
  ! --------------------------------------------------------------------
  
  !!****f* io_module/grab_blips *
  !!
  !!  NAME 
  !!   grab_blips
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Grabs a blips (dumped by dump_blips)
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   16:12, 15/10/2002 drb 
  !!  MODIFICATION HISTORY
  !!   2013/03/06 17:10 dave
  !!   - file name length increased
  !!
  !!  SOURCE
  !!
  subroutine grab_blips(stub,support,inode)

    use datatypes
    use functions_on_grid, only: gridfunctions, fn_on_grid
    use global_module, only: numprocs

    implicit none

    ! Passed variables
    integer :: inode, support
    character(len=*) :: stub

    ! Local variables
    integer :: lun, element, nf1, nf2, ios
    character(len=50) :: filename

    ! Build a filename based on node number
    call get_file_name(stub//'grid',numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename,status='old',iostat=ios)
    if(ios /= 0) call cq_abort('grab_blips: failed to open input file '//filename)
    ! Grab blips
    do element=1,gridfunctions(support)%size
       read(unit=lun,fmt='(f30.15)') gridfunctions(support)%griddata(element)
    end do
    call io_close(lun)
    return
  end subroutine grab_blips
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine dump_locps
  ! --------------------------------------------------------------------
  
  !!****f* io_module/dump_locps *
  !!
  !!  NAME 
  !!   dump_locps
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Dumps the local pseudo
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   10:58, 15/10/2002 drb 
  !!  MODIFICATION HISTORY
  !!   2008/04/02  M. Todorovic
  !!     Added stub name into function
  !!   2013/03/06 17:10 dave
  !!   - file name length increased
  !!
  !!  SOURCE
  !!
  subroutine dump_locps(stub,locps,size,inode)

    use datatypes
    use global_module, only: numprocs

    ! Passed variables
    integer :: size, inode
    real(double), dimension(size) :: locps
    character(len=*) :: stub

    ! Local variables
    integer :: lun, block, n_point, n_i, n
    character(len=50) :: filename

    ! Build a filename based on node number
    call get_file_name('locps'//stub,numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename)
    ! Dump locps
    do block=1, size
       write(unit=lun,fmt='(f30.15)') locps(block)
    end do
    call io_close(lun)
    return
  end subroutine dump_locps
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine grab_locps
  ! --------------------------------------------------------------------
  
  !!****f* io_module/grab_locps *
  !!
  !!  NAME 
  !!   grab_locps
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Grabs the locps (dumped by dump_locps)
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   10:58, 15/10/2002 drb 
  !!  MODIFICATION HISTORY
  !!   2008/04/02  M. Todorovic
  !!     Added stub name into function
  !!   2013/03/06 17:10 dave
  !!   - file name length increased
  !!
  !!  SOURCE
  !!
  subroutine grab_locps(stub,locps,size,inode)

    use datatypes
    use global_module, only: numprocs

    ! Passed variables
    integer :: size, inode
    real(double), dimension(size) :: locps
    character(len=*) :: stub

    ! Local variables
    integer :: lun, block, n_point, n_i, n, ios
    character(len=50) :: filename

    ! Build a filename based on node number
    call get_file_name('locps'//stub,numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename,status='old',iostat=ios)
    if(ios /= 0) call cq_abort('grab_locps: failed to open input file '//filename)
    ! Grab locps
    do block=1, size
       read(unit=lun,fmt='(f30.15)') locps(block)
    end do
    call io_close(lun)
    return
  end subroutine grab_locps
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine dump_projs
  ! --------------------------------------------------------------------
  
  !!****f* io_module/dump_projs *
  !!
  !!  NAME 
  !!   dump_projs
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Dumps the projectors
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   10:58, 15/10/2002 drb 
  !!  MODIFICATION HISTORY
  !!   2013/03/06 17:10 dave
  !!   - file name length increased
  !!
  !!  SOURCE
  !!
  subroutine dump_projs(projs,size,inode)

    use datatypes
    use global_module, only: numprocs

    ! Passed variables
    integer :: size, inode
    real(double), dimension(size) :: projs

    ! Local variables
    integer :: lun, block, n_point, n_i, n
    character(len=15) :: filename

    ! Build a filename based on node number
    call get_file_name('projs',numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename)
    ! Dump projs
    do block=1, size
       write(unit=lun,fmt='(f30.15)') projs(block)
    end do
    call io_close(lun)
    return
  end subroutine dump_projs
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine grab_projs
  ! --------------------------------------------------------------------
  
  !!****f* io_module/grab_projs *
  !!
  !!  NAME 
  !!   grab_projs
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Grabs the charge projs (dumped by dump_projs)
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   10:58, 15/10/2002 drb 
  !!  MODIFICATION HISTORY
  !!   2013/03/06 17:10 dave
  !!   - file name length increased
  !!  SOURCE
  !!
  subroutine grab_projs(projs,size,inode)

    use datatypes
    use global_module, only: numprocs

    ! Passed variables
    integer :: size, inode
    real(double), dimension(size) :: projs

    ! Local variables
    integer :: lun, block, n_point, n_i, n, ios
    character(len=15) :: filename

    ! Build a filename based on node number
    call get_file_name('projs',numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename,status='old',iostat=ios)
    if(ios /= 0) call cq_abort('grab_projs: failed to open input file '//filename)
    ! Grab projs
    do block=1, size
       read(unit=lun,fmt='(f30.15)') projs(block)
    end do
    call io_close(lun)
    return
  end subroutine grab_projs
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine dump_blip_coeffs
  ! --------------------------------------------------------------------
  
  !!****f* io_module/dump_blip_coeffs *
  !!
  !!  NAME 
  !!   dump_blip_coeffs
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Dumps a blip_coeffs
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   16:12, 15/10/2002 drb 
  !!  MODIFICATION HISTORY
  !!   2006/07/18 08:19 dave
  !!    Changed to use new data storage
  !!  SOURCE
  !!
  subroutine dump_blip_coeffs(data_blip,size,inode)

    use datatypes
    use global_module, only: numprocs

    ! Passed variables
    integer :: inode, size
    real(double) :: data_blip(size)

    ! Local variables
    integer :: lun, i,j,k, nf1, nf2
    character(len=18) :: filename

    ! Build a filename based on node number
    call get_file_name('blip_coeffs',numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename)
    ! Dump blip_coeffs
    do k=1,size
       write(unit=lun,fmt='(f30.15)') data_blip(k)
    end do
    call io_close(lun)
    return
  end subroutine dump_blip_coeffs
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine grab_blip_coeffs
  ! --------------------------------------------------------------------
  
  !!****f* io_module/grab_blip_coeffs *
  !!
  !!  NAME 
  !!   grab_blip_coeffs
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Grabs a blip_coeffs (dumped by dump_blip_coeffs)
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   16:12, 15/10/2002 drb 
  !!  MODIFICATION HISTORY
  !!   2006/07/18 08:19 dave
  !!    Changed to use new data storage
  !!  SOURCE
  !!
  subroutine grab_blip_coeffs(data_blip,size,inode)

    use datatypes
    use global_module, only: numprocs

    ! Passed variables
    integer :: inode, size
    real(double) :: data_blip(size)

    ! Local variables
    integer :: lun, i,j,k, nf1, nf2, ios
    character(len=18) :: filename

    ! Build a filename based on node number
    call get_file_name('blip_coeffs',numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename,status='old',iostat=ios)
    if(ios /= 0) call cq_abort('grab_blip_coeffs: failed to open input file '//filename)
    ! Grab blip_coeffs
    do k=1,size
       read(unit=lun,fmt='(f30.15)') data_blip(k)
    end do
    call io_close(lun)
    return
  end subroutine grab_blip_coeffs
  !!***

  
  ! --------------------------------------------------------------------
  ! Subroutine write_matrix
  ! --------------------------------------------------------------------
  
  !!****f* io_module/write_matrix *
  !!
  !!  NAME 
  !!   write_matrix -- Writes out a matrix
  !!  USAGE
  !!   write_matrix(nnode,myid,cmat,c,parts,prim,set,
  !!     ndim1,ndim2,mx_dim1,mx_dim2,plus)
  !!  PURPOSE
  !!   Writes out elements of a matrix, sorted by global number of 
  !!     atoms and neighbours, and calculates periodic supercell offset
  !!  INPUTS
  !!   integer :: nnode, myid, plus
  !!   integer :: mx_dim1,mx_dim2,ndim1,ndim2
  !!   type(group_set)   :: parts
  !!   type(primary_set) :: prim
  !!   type(cover_set)   :: set
  !!   type(matrix), dimension(:) :: cmat
  !!   real(double) :: c(:,:,:)                    
  !!  USES
  !!   datatypes, global_module, basic_types, matrix_module
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !! 
  !!  MODIFICATION HISTORY
  !!   2006/02/22 11:49 dave
  !!    Completely updated to new matrix style
  !!  SOURCE
  !!
  subroutine write_matrix(myid,matA,range,parts,prim,set,plus)

    use datatypes
    use global_module, only: id_glob
    use basic_types, only: group_set, primary_set, cover_set
    use matrix_module, only: matrix, matrix_halo
    use matrix_data, only: halo, mat
    use mult_module, only: return_matrix_value_pos, matrix_pos

    implicit none

    ! Passed variables
    integer :: myid,matA,plus,range
    type(group_set)   :: parts
    type(primary_set) :: prim
    type(cover_set)   :: set

    ! Local variables
    integer :: unit1,unit2,n1,n2,nmodx,nmody,nmodz,ncx,ncy,ncz
    integer :: nnd,ind_part,i,ni,ind_cover,nb,np,icpartp,ind_partn
    integer :: inc,nx,ny,nz,npx,npy,npz,nnabx,nnaby,nnabz
    integer :: ndcellx,ndcelly,ndcellz,ist,isu, gcspart

    unit2=12+myid+plus
    open(unit=unit2)
    nnd=myid+1
    if(prim%n_prim.gt.0) then
       nmodx=1+set%ncoverx/parts%ngcellx
       nmody=1+set%ncovery/parts%ngcelly
       nmodz=1+set%ncoverz/parts%ngcellz
       call start_timer(tmr_std_matrices)
       do np=1,parts%ng_on_node(nnd)
          ind_part=parts%ngnode(parts%inode_beg(nnd)+np-1)
          ncx=1+(ind_part-1)/(parts%ngcelly*parts%ngcellz)
          ncy=1+(ind_part-1-(ncx-1)*parts%ngcelly*parts%ngcellz)/parts%ngcellz
          ncz=ind_part-(ncx-1)*parts%ngcelly*parts%ngcellz-(ncy-1)*parts%ngcellz
          if(prim%nm_nodgroup(np).gt.0) then
             do i=1,prim%nm_nodgroup(np)
                ni=prim%nm_nodbeg(np)+i-1
                if(mat(np,range)%n_nab(i).gt.0) then
                   do nb=1,mat(np,range)%n_nab(i)
                      ind_cover=mat(np,range)%i_part(mat(np,range)%i_acc(i)+nb-1)
                      icpartp=set%lab_cover(ind_cover)
                      ind_partn=set%lab_cell(ind_cover)
                      inc=parts%icell_beg(ind_partn)+&
                           mat(np,range)%i_seq(mat(np,range)%i_acc(i)+nb-1)-1
                      nx=1+(icpartp-1)/(set%ncovery*set%ncoverz)
                      ny=1+(icpartp-1-(nx-1)*set%ncovery*&
                           set%ncoverz)/set%ncoverz
                      nz=icpartp-(nx-1)*set%ncovery*set%ncoverz-&
                           (ny-1)*set%ncoverz
                      npx=nx-1-set%nspanlx-prim%idisp_primx(np)
                      npy=ny-1-set%nspanly-prim%idisp_primy(np)
                      npz=nz-1-set%nspanlz-prim%idisp_primz(np)
                      nnabx=ncx+npx
                      nnaby=ncy+npy
                      nnabz=ncz+npz
                      ndcellx=(nnabx-1+nmodx*parts%ngcellx)/parts%ngcellx-nmodx
                      ndcelly=(nnaby-1+nmody*parts%ngcelly)/parts%ngcelly-nmody
                      ndcellz=(nnabz-1+nmodz*parts%ngcellz)/parts%ngcellz-nmodz
                      write(unit2,596) prim%ig_prim(ni), id_glob(inc),ndcellx,ndcelly,ndcellz
596                   format(2i10,3i5)
                      ist=mat(np,range)%i_acc(i)+nb-1
                      gcspart = set%icover_ibeg(mat(np,range)%i_part(ist))+mat(np,range)%i_seq(ist)-1
                      do n1=1,mat(np,range)%ndimi(i)!ndim1
                         do n2=1,mat(np,range)%ndimj(ist)!ndim2
                            isu = matrix_pos(matA,ni,halo(range)%i_halo(gcspart),n1,n2)
                            write(unit2,597) n1,n2, return_matrix_value_pos(matA,isu)
597                         format(10x,2i5,5x,e20.12)
                         enddo ! n2
                      enddo ! n1
                   enddo ! ncnab
                endif ! ncnab.gt.0
             enddo ! nc_nodpart
          endif ! nc_nodpart.gt.0
       enddo ! np_on_node
       call stop_timer(tmr_std_matrices)
    endif ! n_prim.gt.0
    close(unit=unit2)
    return
  end subroutine write_matrix
  !!***

  !!****f* io_module/banner *
  !!
  !!  NAME 
  !!   banner
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Writes out an initial banner
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !! 
  !!  MODIFICATION HISTORY
  !!   21/06/2001 dave
  !!    Added ROBODoc header
  !!   2013/01/30 10:32 dave
  !!   - Added new authors to list
  !!  SOURCE
  !!
  subroutine banner

    implicit none

    write(io_lun,fmt='(4x,a72)') '________________________________________________________________________'
    write(io_lun,fmt='(4x,a72)') '                                                                        '
    write(io_lun,fmt='(4x,a72)') '                                CONQUEST                                '
    write(io_lun,fmt='(4x,a72)') '                                                                        '
    write(io_lun,fmt='(4x,a72)') '            Concurrent Order N QUantum Electronic STructure             '
    write(io_lun,fmt='(4x,a72)') '________________________________________________________________________'
    write(io_lun,fmt='(4x,a72)') '                                                                        '
    write(io_lun,fmt='(4x,a72)') ' Conquest lead developers:                                              '
    write(io_lun,fmt='(4x,a72)') '  D.R.Bowler (UCL, NIMS), T.Miyazaki (NIMS), A.Nakata (NIMS),           '
    write(io_lun,fmt='(4x,a72)') '  L. Truflandier (U. Bordeaux)                                          '
    write(io_lun,fmt='(4x,a72)') '                                                                        '
    write(io_lun,fmt='(4x,a72)') ' Developers:                                                            '
    write(io_lun,fmt='(4x,a72)') '  M.Arita (NIMS), J.S.Baker (UCL), V.Brazdova (UCL), R.Choudhury (UCL), '
    write(io_lun,fmt='(4x,a72)') '  S.Y.Mujahed (UCL), J.T.Poulton (UCL), Z.Raza (NIMS), A.Sena (UCL),    '
    write(io_lun,fmt='(4x,a72)') '  U.Terranova (UCL), L.Tong (UCL), A.Torralba (NIMS)                    '
    write(io_lun,fmt='(4x,a72)') '                                                                        '
    write(io_lun,fmt='(4x,a72)') ' Early development:                                                     '
    write(io_lun,fmt='(4x,a72)') '  I.J.Bush (STFC), C.M.Goringe (Keele), E.H.Hernandez (Keele)           '
    write(io_lun,fmt='(4x,a72)') '                                                                        '
    write(io_lun,fmt='(4x,a72)') ' Original inspiration and project oversight:                            '
    write(io_lun,fmt='(4x,a72)') '  M.J.Gillan (Keele, UCL)                                               '
    write(io_lun,fmt='(4x,a72)') '________________________________________________________________________'

  end subroutine banner
  !!***

  !!****f* io_module/write_positions *
  !!
  !!  NAME 
  !!   write_positions - output atom positions
  !!  USAGE
  !!   write_positions(partition structure)
  !!  PURPOSE
  !!   Writes out the positions of the atoms in an appropriate
  !!   file so that more work can be done
  !!  INPUTS
  !!   type(group_set) :: parts
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   04/07/2001
  !!  MODIFICATION HISTORY
  !!   19/10/2001 dave
  !!    Added io_assign call (from fdf) for safer file access
  !!    and changed all file numbers to lun. Comment added 22/02/2002 drb
  !!   11:39, 18/03/2003 drb 
  !!    Took from ParaDens into Conquest
  !!   13:20, 18/03/2003 drb 
  !!    Added integer argument to allow multiple outputs
  !!  SOURCE
  !!
  subroutine write_positions(num,parts)

    use basic_types, only: group_set
    use global_module, only: rcellx, rcelly, rcellz, ni_in_cell, &
         x_atom_cell, y_atom_cell, z_atom_cell, numprocs, id_glob
    use species_module, only: species

    implicit none

    ! Passed variables
    type(group_set) :: parts
    integer :: num

    ! Local variables
    integer :: nnd, np, ind_part, ni
    integer :: lun
    character(len=50) :: filename

    if(num>0) then
       ! Build a filename based on iteration number
       call get_file_name('make_prt.out',4,num,filename)
    else
       filename = "make_prt.out"
    end if
    call io_assign(lun)
    open(unit=lun,file=filename)
    write(lun,1) rcellx,rcelly,rcellz
    write(lun,2) ni_in_cell
    write(lun,3) parts%ngcellx,parts%ngcelly,parts%ngcellz
    write(lun,2) numprocs
    do nnd=1,numprocs
       write(lun,4) nnd,parts%ng_on_node(nnd),parts%inode_beg(nnd)
       if(parts%ng_on_node(nnd)>0) then
          do np=1,parts%ng_on_node(nnd)
             ind_part=parts%ngnode(parts%inode_beg(nnd)+np-1)
             write(lun,5) np,ind_part,parts%nm_group(ind_part),&
                  parts%icell_beg(ind_part)
             if(parts%nm_group(ind_part)>0) then
                do ni=1,parts%nm_group(ind_part)
                   write(lun,6) ni,id_glob(parts%icell_beg(ind_part)+ni-1),&
                        x_atom_cell(parts%icell_beg(ind_part)+ni-1),&
                        y_atom_cell(parts%icell_beg(ind_part)+ni-1),&
                        z_atom_cell(parts%icell_beg(ind_part)+ni-1),&
                        species(parts%icell_beg(ind_part)+ni-1)
                end do ! Loop over atoms
             end if ! Atoms in group
          end do ! Loop over groups
       end if ! Processor have groups
    end do ! Loop over processors
    call io_close(lun)
1   format(3e20.12)
2   format(i10)
3   format(3i5)
4   format(3i10)
5   format(5x,4i10)
6   format(2i8,3e20.12,i3)
    return
  end subroutine write_positions
  !!***

  !!****f* io_module/write_xsf *
  !!
  !!  NAME 
  !!   write_xsf
  !!  PURPOSE
  !!   Writes atomic positions to a .xsf file, viewable using VMD
  !!  INPUTS
  !!   type(group_set) :: parts
  !!  USES
  !! 
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/10/26
  !!  MODIFICATION HISTORY
  !!   
  !!  SOURCE
  !!
  subroutine write_xsf(filename, step)

    use datatypes
    use numbers,        only: zero
    use dimens,         only: r_super_x, r_super_y, r_super_z
    use global_module,  only: ni_in_cell, iprint_init, atom_coord, &
                              species_glob
    use species_module, only: species_label
    use GenComms,       only: inode, ionode, cq_abort
    use units,          only: BohrToAng
    use timer_module

    ! Passed variables
    character(len=*)  :: filename
    integer           :: step

    ! Local variables
    integer                    :: lun, i
    character(len=2)           :: atom_name

    if(inode==ionode) then
      if (iprint_init>3) write(io_lun,fmt='(6x,a)') 'Writing atomic positions to .xsf'
      call io_assign(lun)
      if(append_coords) then
         open(unit=lun,file=filename,position='append')
         write(lun,*)
      else
         open(unit=lun,file=filename)
      end if
      write(lun,'(a)') "CRYSTAL"
      write(lun,'("PRIMVEC   ",i8)') step
      write(lun,fmt='(3f14.8)') r_super_x*BohrToAng, zero, zero
      write(lun,fmt='(3f14.8)') zero, r_super_y*BohrToAng, zero
      write(lun,fmt='(3f14.8)') zero, zero, r_super_z*BohrToAng
      write(lun,'("PRIMCOORD ",i8)') step
      write(lun,fmt='(2i8)') ni_in_cell, 1
      do i=1,ni_in_cell
        atom_name = adjustr(species_label(species_glob(i))(1:2))
        write(lun,'(a4,3f16.8)') atom_name, atom_coord(:,i)*BohrToAng
                 ! species_glob(i),flag_move_atom(1,i),flag_move_atom(2,i), &
      end do
      call io_close(lun)
    end if
  end subroutine write_xsf
  !!***

  !!****f* io_module/write_extxyz *
  !!
  !!  NAME 
  !!   write_extxyz
  !!  PURPOSE
  !!   Writes atomic positions, including atomic forces, lattice vectors, 
  !!  energy, system signature, etc. to an extended format .xyz file   
  !!  which can be directly recognized by Python Library named ASE.
  !!  The units are Angstrom, eV and eV/Angstrom for distance, energy and
  !!  forces.
  !!  
  !!  INPUTS
  !!  
  !!  USES
  !! 
  !!  AUTHOR
  !!   Jianbo Lin
  !!  CREATION DATE
  !!   2021/10/18
  !!  MODIFICATION HISTORY
  !!   
  !!  SOURCE
  !!
  subroutine write_extxyz(filename, energy0, atom_force)

    use datatypes
    use timer_module
    use numbers,        only: zero
    use dimens,         only: r_super_x, r_super_y, r_super_z
    use global_module,  only: ni_in_cell, iprint_init, atom_coord, &
                              species_glob
    use species_module, only: species_label
    use GenComms,       only: inode, ionode, cq_abort
    use units,          only: BohrToAng, HaToeV

    ! Passed variables
    character(len=*)                      :: filename
    real(double)                          :: energy0
    real(double), dimension(3,ni_in_cell) :: atom_force

    ! Local variables
    integer                    :: lun, i, j, title_length
    character(len=2)           :: atom_name
    character(len=432)         :: comment
    character(len=45)          :: vec_a, vec_b, vec_c, energy_str
    character(len=80)          :: titles_xyz
    real(double)               :: for_conv_loc, en_conv_loc, dist_conv_loc

    if(inode==ionode) then
      if (iprint_init>2) write(io_lun, &
          '(6x,"Writing atomic positions to ",a,".xyz")') filename
      call io_assign(lun)
      if(append_coords) then
         open(unit=lun,file=filename,position='append')
      else
         open(unit=lun,file=filename)
      end if

      en_conv_loc = HaToeV
      dist_conv_loc = BohrToAng
      for_conv_loc = en_conv_loc/dist_conv_loc
      write(lun,'(i0)') ni_in_cell
	  
      ! Transfer space to underbar in the titles
      titles_xyz = TRIM(titles)
      title_length = len(TRIM(titles))
      do i=1, title_length
        if (titles_xyz(i:i) == ' ') titles_xyz(i:i) = '_' 
      end do
      ! Add information about system signature
      comment = 'config_type='//TRIM(titles_xyz)

      ! Add information about lattice and energy
      write(vec_a,fmt='(3f15.8)') r_super_x*BohrToAng, zero, zero
      write(vec_b,fmt='(3f15.8)') zero, r_super_y*BohrToAng, zero
      write(vec_c,fmt='(3f15.8)') zero, zero, r_super_z*BohrToAng
      comment=TRIM(comment)//' Lattice="'//ADJUSTL(vec_a)//ADJUSTL(vec_b)//TRIM(ADJUSTL(vec_c))//'" '
      comment=TRIM(comment)//' Properties=species:S:1:pos:R:3:forces:R:3 potential_energy='
      write(energy_str,'(f0.8)') energy0 * en_conv_loc
      comment = TRIM(comment)//TRIM(energy_str)//' pbc="T T T" '
      write(lun,'(a)') TRIM(comment)

      do i=1,ni_in_cell
        atom_name = adjustr(species_label(species_glob(i))(1:2))
        write(lun,'(a4,6f16.8)') atom_name, atom_coord(:,i)*dist_conv_loc, &
                 (for_conv_loc*atom_force(j,i), j = 1, 3)
                 ! species_glob(i),flag_move_atom(1,i),flag_move_atom(2,i), &
      end do
      call io_close(lun)
    end if
  end subroutine write_extxyz
  !!***

  !!****f* io_module/write_xyz *
  !!
  !!  NAME 
  !!   write_xyz
  !!  PURPOSE
  !!   Writes atomic positions to a .xyz file
  !!  INPUTS
  !!   
  !!  USES
  !! 
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2019/02/13
  !!  MODIFICATION HISTORY
  !!   
  !!  SOURCE
  !!
  subroutine write_xyz(filename, comment)

    use datatypes
    use numbers,        only: zero
    use dimens,         only: r_super_x, r_super_y, r_super_z
    use global_module,  only: ni_in_cell, iprint_init, atom_coord, &
                              species_glob
    use species_module, only: species_label
    use GenComms,       only: inode, ionode, cq_abort
    use units,          only: BohrToAng
    use timer_module

    ! Passed variables
    character(len=*)  :: filename, comment

    ! Local variables
    integer                    :: lun, i
    character(len=2)           :: atom_name

    if(inode==ionode) then
      if (iprint_init>2) write(io_lun, &
          '(6x,"Writing atomic positions to ",a,".xyz")') filename
      call io_assign(lun)
      open(unit=lun,file=filename)
      write(lun,fmt='(i8)') ni_in_cell
      write(lun,'(a)') comment
      do i=1,ni_in_cell
        atom_name = adjustr(species_label(species_glob(i))(1:2))
        write(lun,'(a4,3f16.8)') atom_name, atom_coord(:,i)*BohrToAng
      end do
      call io_close(lun)
    end if
  end subroutine write_xyz
  !!***

  ! Hopefully we'll never need this kludgy but portable way of flushing buffers !
  !  subroutine force_buffers(lun)
  !
  !    use GenComms, only: inode, ionode
  !
  !    implicit none
  !
  !    integer :: lun
  !
  !    integer :: i, len
  !    character(len=255) :: fname
  !
  !    if(inode==ionode) then
  !       inquire(unit=lun,name=fname)
  !       if(fname(1:1)==' ') return
  !       do i=1,255
  !          if(fname(i:i)==' ') then
  !             len=i
  !             exit
  !          end if
  !       end do
  !       close(lun)
  !       open(unit=lun,file=fname(1:len),position='APPEND')
  !    end if
  !    return
  !  end subroutine force_buffers
  
  !!****f* io_module/get_file_name *
  !!
  !!  NAME
  !!   get_file_name - creates a file name with a process number in it
  !!  USAGE
  !!   get_file_name(fileroot,numprocs,inode,filename)
  !!  PURPOSE
  !!   Returns a file name with proper formating, including 
  !!   node identification number an zero padding to reflect 
  !!   the total number of processes. It appends a dot after the root
  !!  INPUTS
  !!   character(lem=*) :: fileroot  ! The root of the filename, e.g. chden
  !!   integer :: numprocs           ! Number of processes; determines zero padding
  !!   integer :: inode              ! Who am I ?
  !!  OUTPUTS
  !!   character(len=*) :: filanem   ! The formatted filename
  !!  USES
  !!
  !!  AUTHOR
  !!   A.S. Torralba
  !!  CREATION DATE
  !!   23/07/2009
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine get_file_name(fileroot, numprocs, inode, filename)

    use datatypes

    implicit none

    ! Passed variables
    character(len=*), intent(in) :: fileroot
    integer, intent(in) :: numprocs
    integer, intent(in) :: inode
    character(len=*), intent(out) :: filename

    ! Parameters
    integer, parameter :: maxlen=80

    ! Local variables
    integer :: i
    integer :: padzeros
    character(len=maxlen) :: num

    if (LEN_TRIM (fileroot) + 1 > maxlen) &
         call cq_abort('get_file_name: error : string overflow')
    filename = TRIM (fileroot)//'.'
    padzeros = MAX (3, &
         FLOOR (1.0 + LOG10(REAL (numprocs, kind=double)))) - &
         FLOOR (1.0 + LOG10(REAL (inode, kind=double)))
    write (num,'(i80)') inode
    num = ADJUSTL (num)
    if (LEN_TRIM(fileroot) + padzeros + LEN_TRIM (num) + 1 > maxlen) &
         call cq_abort ('get_file_name: error : string overflow')
    do i = 1, padzeros
       filename = TRIM (filename)//'0'
    end do
    filename = TRIM (filename)//num
  end subroutine get_file_name
  !!***

  !!****f* io_module/get_file_name_2rank *
  !!
  !!  NAME
  !!   get_file_name_2rank - creates a file name with an index (MD step), and process number.
  !!  USAGE
  !!   get_file_name_2rank(fileroot,filename,step,inode)
  !!  PURPOSE
  !!   Returns a file name with the number of step, and inode (optional) 
  !!   At present, we assume 6 digits for each number, though it can be 
  !!   easily changed by changing the line 
  !!  INPUTS
  !!   character(lem=*) :: fileroot  ! The root of the filename, e.g. chden
  !!   integer :: step               ! MD step or some index
  !!   integer :: inode              ! Process Number
  !!  OUTPUTS
  !!   character(len=*) :: filename  ! The formatted filename
  !!  USES
  !!
  !!  AUTHOR
  !!   T. Miyazaki
  !!  CREATION DATE
  !!   02/06/2017
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine get_file_name_2rank(fileroot, filename, index, inode)
    use datatypes
    implicit none

    ! Passed variables
    character(len=*), intent(in) :: fileroot
    integer, intent(in) :: index
    integer, intent(in),optional :: inode
    character(len=*), intent(out) :: filename

    ! Parameters
    integer, parameter :: maxlen=80

    ! Local variables
    integer :: i
    character(len=maxlen) :: num_index
    character(len=maxlen) :: num_inode

    if (LEN_TRIM (fileroot) + 1 > maxlen) &
         call cq_abort('get_file_name: error : string overflow')

    write (num_index,'(i2.2)') index
     filename = TRIM (fileroot)//'.i'//num_index

    if(present(inode)) then
     write (num_inode,'(i6.6)') inode
     filename = TRIM (filename)//'.p'//num_inode
    endif
 
   return
  end subroutine get_file_name_2rank
  !!***


  !!****f* io_module/print_process_info *
  !!
  !!  NAME
  !!   print_process_info - prints MPI and unix PIDs and hostname
  !!  USAGE
  !!   get_file_name()
  !!  PURPOSE
  !!   Prints a line matching MPI ID to process ID and host,
  !!   so that output files can be identified. Useful to compare
  !!   memory use, for example.
  !!  INPUTS
  !!
  !!  OUTPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   A.S. Torralba
  !!  CREATION DATE
  !!   05/08/2009
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine print_process_info()

    use global_module, only: numprocs
    use GenComms, only: inode,ionode
    use mpi

    implicit none

    ! Parameters
    integer, parameter :: hostlen=80

    ! Local variables
    integer :: mypid
    integer :: getpid        ! Intrinsic
    character(len=hostlen) :: myhost
    character(len=hostlen), allocatable :: hostnames(:)
    integer, allocatable, dimension(:) :: pids
    integer :: stat, i
    integer :: mpistat(MPI_STATUS_SIZE), ierr

    ! Get process ID and hostname for this process
    mypid=getpid()
    call hostnm(myhost)

    ! Communicate the information to the ionode, so it can print
    if(inode==ionode) then
       allocate(hostnames(numprocs),STAT=stat)
       if(stat/=0) call cq_abort('Error allocating hostnames !')
       allocate(pids(numprocs),STAT=stat)
       if(stat/=0) call cq_abort('Error allocating pids !')

       do i=1,numprocs
          if(i /= ionode) then
             call MPI_recv(hostnames(i),hostlen,MPI_CHARACTER,i-1,1,MPI_COMM_WORLD,mpistat,ierr)
             if(ierr /= MPI_SUCCESS) call cq_abort('Error receiving hostnames ! ',i,ierr)
             call MPI_recv(pids(i),1,MPI_INTEGER,i-1,1,MPI_COMM_WORLD,mpistat,ierr)
             if(ierr /= MPI_SUCCESS) call cq_abort('Error receiving pids ! ',i,ierr)
          else
             hostnames(i)=myhost
             pids(i)=mypid
          end if
       end do

       write(io_lun,*)
       do i=1,numprocs
          write(io_lun,'(4x,a,i9,a,i9,2a)') 'ID-Info: MPI = ',i,'    Process = ',pids(i),'    Host = ',hostnames(i)
       end do
       write(io_lun,*)

       deallocate(hostnames)
       deallocate(pids)
    else
       call MPI_send(myhost,hostlen,MPI_CHARACTER,ionode-1,1,MPI_COMM_WORLD,ierr)
       if(ierr /= MPI_SUCCESS) call cq_abort('Error sending hostnames ! ',i,ierr)
       call MPI_send(mypid,1,MPI_INTEGER,ionode-1,1,MPI_COMM_WORLD,ierr)
       if(ierr /= MPI_SUCCESS) call cq_abort('Error sending pids ! ',i,ierr)
    end if

  end subroutine print_process_info
  !!***

  !!****f* io_module/print_atomic_positions *
  !!
  !!  NAME
  !!   print_atomic_positions
  !!  USAGE
  !!   print_atomic_positions
  !!  PURPOSE
  !!   Prints atomic positions to the output file
  !!  INPUTS
  !!
  !!  OUTPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2020/03/11
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine print_atomic_positions

    use global_module, only: atom_coord, iprint_MD, ni_in_cell, species_glob
    use dimens,         only: r_super_x, r_super_y, r_super_z, atomicnum
    use GenComms, only: inode, ionode
    use units, only: dist_conv, d_units, dist_units, BohrToAng, bohr
    use periodic_table, only: pte
    use pseudo_tm_info, only: pseudo

    implicit none

    integer :: i

    if(inode==ionode) then
       write(io_lun,fmt='(/6x,"Simulation cell dimensions: ",f10.4,a3," x ",f10.4,a3," x ",f10.4,a3)') &
            r_super_x*dist_conv, d_units(dist_units), r_super_y*dist_conv, d_units(dist_units), &
            r_super_z*dist_conv, d_units(dist_units)
       if(flag_coords_xyz) then
          write(io_lun,fmt='(6x,"           X         Y         Z")')
          if(dist_units==bohr) then
             write(io_lun,fmt='(/6x,"Atomic coordinates in XYZ format (",a2,")")') "A "
             do i = 1, ni_in_cell
                write (io_lun,fmt='(4x, a2, 3f10.4)') pte(atomicnum(species_glob(i))), atom_coord(1:3,i)*BohrToAng
             end do
             write(io_lun,fmt='(8x,"N.B. units above converted to Angstroms for xyz output")')
          else
             write(io_lun,fmt='(/6x,"Atomic coordinates (",a2,")")') d_units(dist_units)
             do i = 1, ni_in_cell
                write (io_lun,fmt='(4x, a2, 3f10.4)') pte(atomicnum(species_glob(i))), atom_coord(1:3,i)
             end do
          end if
       else
          write(io_lun,fmt='(/6x,"Atomic coordinates (",a2,")")') d_units(dist_units)
          write(io_lun,fmt='(6x,"   Atom         X         Y         Z  Species")')
          do i = 1, ni_in_cell
             write (io_lun,fmt='(6x, i7, 3f10.4, 6x, i3)') i,atom_coord(1:3,i), species_glob(i)
          end do
       end if
    end if
    return
    
  end subroutine print_atomic_positions
  !!***
  
  !!****f* io_module/write_velocity *
  !!
  !!  NAME
  !!   write_velocity - prints out velocities of the atoms 
  !!  USAGE
  !!   write_velocity 
  !!  PURPOSE
  !!   
  !!  INPUTS
  !!   character(len=*) :: filename  ! The formatted filename
  !!   real(double) :: velocity(1:3, ni_in_cell)
  !!  OUTPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   T. Miyazaki
  !!  CREATION DATE
  !!   15/06/2010  TM
  !!  MODIFICATION HISTORY
  !!   2011/10/20 09:59 dave
  !!    Syntax correction for rewind
  !!  SOURCE
  !!
  subroutine write_velocity(velocity, filename)
    use numbers
    use global_module, only: id_glob, ni_in_cell,id_glob_inv
    use GenComms, only: inode, ionode, cq_abort, gcopy

    implicit none

    ! Passed variables
    real(double), intent(in) :: velocity(1:3, ni_in_cell)
    character(len=50),intent(in) :: filename

    ! Local variables
    integer :: id_global, ni 
    integer :: lun

    if(inode == ionode) then
       call io_assign(lun)
       open(unit=lun,file=filename)
       rewind(unit=lun)
       do id_global=1,ni_in_cell
          ni=id_glob_inv(id_global)
          if(id_glob(ni) /= id_global) &
               call cq_abort(' ERROR in global labelling ',id_global,id_glob(ni))
          write(lun,101) id_global, ni, velocity(1:3, ni)
101       format(2i8,3e20.12)
       enddo
       call io_close(lun)
    endif  ! (inode == ionode) then

    return
  end subroutine write_velocity
  !!***


  !!****f* io_module/read_velocity *
  !!
  !!  NAME
  !!   read_velocity - reads velocities of the atoms from the file, and
  !!                  send them to all processors
  !!  USAGE
  !!   read_velocity 
  !!  PURPOSE
  !!   
  !!  INPUTS
  !!   character(len=*) :: filename  ! The formatted filename
  !!  OUTPUTS
  !!   real(double) :: velocity(1:3, ni_in_cell)
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   T. Miyazaki
  !!  CREATION DATE
  !!   15/06/2010  TM
  !!  MODIFICATION HISTORY
  !!   2011/10/20 09:59 dave
  !!    Syntax correction for rewind
  !!  SOURCE
  !!
  subroutine read_velocity(velocity, filename)
    use numbers
    use global_module, only: id_glob, ni_in_cell,id_glob_inv, iprint_init
    use GenComms, only: inode, ionode, cq_abort, gcopy

    implicit none

    ! Passed variables
    real(double), intent(out) :: velocity(1:3, ni_in_cell)
    character(len=50),intent(in) :: filename

    ! Local variables
    integer :: id_global, ni , ni2, id_tmp
    integer :: lun

    if(inode == ionode) then
       if(iprint_init>2) write(io_lun,'(10x,a40,a20)') 'Entering read_velocity; reading ', filename
       call io_assign(lun)
       open(unit=lun,file=filename)
       rewind(unit=lun)
       do id_global=1,ni_in_cell
          ni=id_glob_inv(id_global)
          if(id_glob(ni) /= id_global) &
               call cq_abort(' ERROR in global labelling ',id_global,id_glob(ni))
          read(lun,101) id_tmp, ni2, velocity(1:3, ni)
          if(ni2 /= ni) write(io_lun,fmt='(4x,a,i6,a,i6,a,i6)') &
               ' Order of atom has changed for global id (file_labelling) = ',id_global, &
               ' : corresponding labelling (NOprt labelling) used to be ',ni2,&
               ' : but now it is ',ni
101       format(2i8,3e20.12)
       enddo
       call io_close(lun)
    endif  ! (inode == ionode) then
    call gcopy(velocity,3,ni_in_cell)

    return
  end subroutine read_velocity
  !!***

  !!****f* io_module/read_fire *
  !!
  !!  NAME 
  !!   read_fire
  !!  USAGE
  !!   
  !!  PURPOSE
  !!   Reads parameters for FIRE on I/O process and distributes
  !!  INPUTS
  !!   
  !!   
  !!  USES
  !!   
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2015/06/19
  !!  MODIFICATION HISTORY
  !!    
  !!  SOURCE
  !!
  subroutine read_fire(fire_N, fire_N2, fire_P0, MDtimestep, fire_alpha)

    use datatypes
    use global_module, ONLY: iprint_MD
    use GenComms,       only: myid, gcopy

    implicit none

    integer :: fire_N, fire_N2
    real(double) :: MDtimestep, fire_P0, fire_alpha

    integer :: ios, lun

    fire_N = 0
    fire_N2 = 0
    fire_alpha = 0.0_double
    fire_P0 = 0.0_double
    if(myid==0) then
       call io_assign(lun)
       open(unit=lun,file='CQ.FIRE',status='old',iostat=ios)
       if (ios == 0) then
          read (lun,*) fire_N, fire_N2, fire_P0, MDtimestep, fire_alpha
          call io_close(lun)
       else
          fire_N = 0
          fire_P0 = 1.0_double
          fire_alpha = 0.1_double
          fire_N2 = 0
       endif
       if(iprint_MD>1) write (io_lun,fmt="(4x,'Initial FIRE parameters:', 2x,'N= ', i3, 2x, &
            &'N2= ', i3, 2x, 'P= ', f12.8, 2x, 'dt= ', f9.6, 2x, 'alpha= ', f9.6)") &
             fire_N, fire_N2, fire_P0, MDtimestep, fire_alpha
    else
       MDtimestep = 0.0_double ! So that we can distribute with gcopy
    end if
    call gcopy(fire_N)
    call gcopy(fire_N2)
    call gcopy(fire_alpha)
    call gcopy(fire_P0)
    call gcopy(MDtimestep)
    return
    
  end subroutine read_fire
  !!***
  
  !!****f* io_module/write_fire *
  !!
  !!  NAME 
  !!   write_fire
  !!  USAGE
  !!   
  !!  PURPOSE
  !!   Writes parameters for FIRE on I/O process 
  !!  INPUTS
  !!   
  !!   
  !!  USES
  !!   
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2015/06/19
  !!  MODIFICATION HISTORY
  !!    
  !!  SOURCE
  !!
  subroutine write_fire(fire_N, fire_N2, fire_P0, MDtimestep, fire_alpha)

    use datatypes
    use GenComms,       only: myid

    implicit none

    integer :: fire_N, fire_N2
    real(double) :: MDtimestep, fire_P0, fire_alpha

    integer :: lun

    if(myid==0) then
       call io_assign(lun)
       open(unit=lun,file='CQ.FIRE')
       write(lun,fmt='(2i4,3f22.15)') fire_N, fire_N2, fire_P0, MDtimestep, fire_alpha
       call io_close(lun)
    end if
    return
    
  end subroutine write_fire
  !!***
  
  !!****f* io_module/check_stop *
  !!
  !!  NAME 
  !!   check_stop
  !!  USAGE
  !!   
  !!  PURPOSE
  !!   Checks for a CQ.stop file for user-requested stop
  !!  INPUTS
  !!   logical :: flag_userstop
  !!   
  !!  USES
  !!   
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2011/10/05
  !!  MODIFICATION HISTORY
  !!   2011/10/20 09:58 dave
  !!    Bug fix to call io_close
  !!   2018/01/17 TM
  !!    Introduced the control by Elapsed time & Iteration number in CQ.stop
  !!  SOURCE
  !!  
  subroutine check_stop(flag_userstop,iter)

    use datatypes, only: double
    use numbers, only: very_small
    use GenComms, only: inode, ionode, gsum, mtime

    implicit none

    ! Passed variables
    logical :: flag_userstop
    integer, OPTIONAL :: iter
    real(double) :: present_time

    ! Local variables
    integer :: lun, ios, ios2, iter_max

    flag_userstop = .false.

    ! Option 1: Check "CQ.stop"
    if(inode==ionode) then
       call io_assign(lun)
       open(unit=lun,file='CQ.stop',status='old',iostat=ios)
       if(ios==0) then
        rewind lun
        read(lun,*,iostat=ios2) iter_max
        if(ios2 /= 0) then
          flag_userstop = .true.
        endif
        if(PRESENT(iter)) then
         if(iter_max>0 .and. iter > iter_max) flag_userstop = .true.
        endif
       endif
       if(flag_userstop) then
          if(PRESENT(iter)) then
             write(io_lun,fmt='(4x,"User requested stop at ieration ",i4)') iter
          else
             write(io_lun,fmt='(4x,"User requested stop")')
          end if
       end if
       call io_close(lun)
    endif
    call gsum(flag_userstop)

    ! Option 2: Check elapsed time < time_max
    if(time_max > very_small) then
     present_time = mtime()/1000.e0_double
     if(present_time > time_max) flag_userstop=.true.
    endif
  
    ! Set flag_MatrixFile_BinaryFormat_Dump
     if(flag_userstop) &
      flag_MatrixFile_BinaryFormat_Dump = flag_MatrixFile_BinaryFormat_Dump_END 
      !then, if we call dump_matrix2, the format of the file should be changed.
      !      2019/Nov/13 tsuyoshi

    return
  end subroutine check_stop
  !!***

  function return_prefix(name, level)

    character(len=120) :: return_prefix
    character(len=*)   :: name
    character(len=10):: prefix = "          "
    integer          :: level

    return_prefix = prefix(1:-2*level)//name
  end function return_prefix


  ! --------------------------------------------------------------------
  ! Subroutine write_output_ase
  ! --------------------------------------------------------------------
  
  !!****f* io_module/write_output_ase *
  !!
  !!  NAME 
  !!   read_atomic_positions
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Reads positions, species and constraint type for atoms
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !!   datatypes, dimens, GenComms, global_module, species_module
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   16:26, 2003/04/15
  !!  MODIFICATION HISTORY
  !!   16:46, 2003/06/09 tm
  !!    id_glob_inv, atom_coord, species_glob is added.
  !!    and labelling in make_prt.dat is changed 
  !!   13:31, 22/09/2003 drb 
  !!    Added reading of flags to constrain atom movement in three directions
  !!   2006/10/09 08:21 dave
  !!    Tidying up output
  !!   20/11/2006 Veronika
  !!    Corrected reading of constraints from pdb file
  !!   15:05, 27/04/2007 drb 
  !!    Check to ensure right number of species added
  !!   2007/06/28, 21:37 mt + drb
  !!    Added coordinate wrapping to non-pdb coordinates.
  !!   2007/10/11 Veronika
  !!    Added coordinate wrapping to pdb coordinates
  !!    Added coordinate wrapping for coordinates outside 
  !!    the [-1*cell parameter; 2*cell parameter] range
  !!   2013/07/01 M.Arita & T.Miyazaki
  !!    Added shift_in_bohr when wrapping atoms and allocation of atom_coord_diff
  !!   2013/08/20 M.Arita
  !!    Bug fix & correct the if-statement
  !!   2015/06/08 lat
  !!    Added experimental backtrace
  !!   2018/01/22 tsuyoshi (with dave)
  !!    Allocate atom_coord_diff for all calculations
  !!   2019/04/04 14:17 dave
  !!    Correct bug in wrapping with non-fractional coordinates and Angstroms
  !!   2020/07/27 tsuyoshi
  !!    Added atom_vels  
  !!   2020/10/07 tsuyoshi
  !!    Removed allocation of atom_vels (moved to "control")
  !!  SOURCE
  !!
end module io_module

