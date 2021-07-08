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
    use GenComms,       only: inode, ionode, cq_abort
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
               write(io_lun,*) 'Entering read_atomic_positions, pdb file'
          call io_assign(lun)
          ! Go through the file and count the atoms
          open( unit=lun, file=filename, status='old', iostat=ios)
          if ( ios > 0 ) call cq_abort('Reading pdb file: file error')
          if (iprint_init > 2) &
               write (io_lun,'(1x,a)') &
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
               write(io_lun,'(1x,a,i5)') 'Number of atoms: ', ni_in_cell

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
                                     fmt='(2x,"** WARNING ! ** &
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
                      if (iprint_init > 0) &
                           write (io_lun,fmt='(3x, i7, 3f15.8, i3, 3L2)') &
                                 i, atom_coord(1:3,i), species_glob(i), &
                                 flag_move_atom(1:3,i)
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
                   if (iprint_init > 0) &
                        write (io_lun,fmt='(3x, i7, 3f15.8, i3, 3L2)')&
                              i, atom_coord(1:3,i), species_glob(i), &
                              flag_move_atom(1:3,i)
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
                write(io_lun,4) r_super_x, r_super_y, r_super_z
4               format(/10x,'The simulation box has the following dimensions',/, &
                       10x,'a = ',f9.5,' b = ',f9.5,' c = ',f9.5,' a.u.')
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
       else
          if(iprint_init>2) write(io_lun,'(10x,a40,a20)') 'Entering read_atomic_positions; reading ', filename
          call io_assign(lun)
          open(unit=lun,file=filename,status='old')
          ! Read supercell vector - for now it must be orthorhombic so
          ! we use x and y as dummy variables
          read(lun,*) r_super_x, x, y
          read(lun,*) x,r_super_y, y
          read(lun,*) x,y,r_super_z
          read(lun,*) ni_in_cell
         !2010.06.25 TM (Angstrom Units in coords file, but not pdb)
          if(dist_units == ang) then
            r_super_x = r_super_x*AngToBohr
            r_super_y = r_super_y*AngToBohr
            r_super_z = r_super_z*AngToBohr
          endif
         !2010.06.25 TM (Angstrom Units in coords file, but not pdb)

          allocate(flag_move_atom(3,ni_in_cell),atom_coord(3,&
                   ni_in_cell),species_glob(ni_in_cell),STAT=stat)
          if(stat/=0) &
               call cq_abort("Failure to allocate coordinates: ",ni_in_cell)
          call reg_alloc_mem(area_init, 3*ni_in_cell,type_dbl)
          call reg_alloc_mem(area_init, 4*ni_in_cell,type_int)
          do i=1,ni_in_cell
             read(lun,*) x,y,z,species_glob(i),movex,movey,movez
             if(species_glob(i)>n_species) then
                write(io_lun,fmt='(2x,"** WARNING ! ** Species &
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
             !          id_glob(i) = i
!!$ LAT: put in write info
            if(iprint_init>0) &
                 write (io_lun,fmt='(3x, i7, 3f15.8, i3, 3L2)') &
                       i,atom_coord(1:3,i), species_glob(i), &
                       flag_move_atom(1:3,i)
          end do
          call io_close(lun)
       end if pdb
    end if
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
       write(lun,fmt='("# Ef: ",f12.5)') Ef(1)
    else
       write(lun,fmt='("# Ef: ",2f12.5)') Ef(1),Ef(2)
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

    write(io_lun,1) 

1   format(/12x, &
         '______________________________________________________',/,12x, &
         '______________________________________________________',/,12x, &
         '                                                      ',/,12x, &
         '                        CONQUEST                      ',/,12x, &
         '                                                      ',/,12x, &
         '    Concurrent Order N QUantum Electronic STructure   ',/,12x, &
         '______________________________________________________',/,12x, &
         '                                                      ',/,12x, &
         '                       Written by:                    ',/,12x, &
         '                                                      ',/,12x, &
         '       David Bowler           Tsuyoshi Miyazaki       ',/,12x, &
         '        (UCL,NIMS)                 (NIMS)             ',/,12x, &
         '                                                      ',/,12x, &
         '       Ayako Nakata              Zamaan Raza          ',/,12x, &
         '          (NIMS)                   (NIMS)             ',/,12x, &
         '                                                      ',/,12x, &
         '       Jack Poulton           Lionel Truflandier      ',/,12x, &
         '          (UCL)                  (Bordeaux)           ',/,12x, &
         '                                                      ',/,12x, &
         '      Shereif Mujahed            Jack Baker           ',/,12x, &
         '          (UCL)                    (UCL)              ',/,12x, &
         '                                                      ',/,12x, &
         '     Antonio Torralba         Veronika Brazdova       ',/,12x, &
         '        (UCL,NIMS)                 (UCL)              ',/,12x, &
         '                                                      ',/,12x, &
         '      Lianheng Tong            Michiaki Arita         ',/,12x, &
         '          (UCL)                    (NIMS)             ',/,12x, &
         '                                                      ',/,12x, &
         '        Alex Sena             Umberto Terranova       ',/,12x, &
         '          (UCL)                    (UCL)              ',/,12x, &
         '                                                      ',/,12x, &
         '     Rathin Choudhury           Mike Gillan           ',/,12x, &
         '          (UCL)                    (UCL)              ',/,12x, &
         '                                                      ',/,12x, &
         '               Early Development by                   ',/,12x, & 
         '                                                      ',/,12x, &
         '       Chris Goringe             Edward Hernandez     ',/,12x, &
         '          (Keele)                     (Keele)         ',/,12x, &
         '                                                      ',/,12x, &
         '                      Ian Bush                        ',/,12x, &
         '                     (Daresbury)                      ',/,12x, &
         '______________________________________________________',/,12x, &
         '______________________________________________________',/,/)

  end subroutine banner
  !!***
  
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

end module io_module

