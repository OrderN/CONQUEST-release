! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!
module io_module

  use global_module,          only: io_lun
  use GenComms,               only: cq_abort, gcopy
  use timer_stdclocks_module, only: start_timer, stop_timer, tmr_std_matrices
  use input_module,           only: leqi, io_assign, io_close

  implicit none

  logical           :: pdb_format, append_coords, pdb_output 
  character(len=1)  :: pdb_altloc
  character(len=80) :: pdb_template

  ! RCS tag for object file identification 
  character(len=80), save, private :: &
       RCSid = "$Id$"

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
                              flag_move_atom, area_init
    use species_module, only: species, species_label, n_species
    use GenComms,       only: inode, ionode, cq_abort
    use memory_module,  only: reg_alloc_mem, type_dbl, type_int
    use units,          only: AngToBohr
    use units,          only: dist_units, ang
    use numbers,        only: very_small, zero

    ! Passed variables
    character(len=*) :: filename

    ! Local variables
    integer :: lun, i, spec, stat
    real(double) :: x, y, z
    logical :: movex, movey, movez
    real(double), dimension(3) :: cell
    ! Local variables for reading pdb files
    integer :: ios, j
    character(len=80) :: pdb_line
    real(double), dimension(3) :: angle
    integer :: num_move_atom
    real(double) :: num_move_atom_real
    character(len=2) :: atom_name

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
               write (io_lun,'(x,a)') &
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
               write(io_lun,'(x,a,i5)') 'Number of atoms: ', ni_in_cell

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
                          very_small ) &
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
                   if (num_move_atom_real - (real(num_move_atom)) > very_small ) &
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
                if ((atom_coord(j,i) < zero) .or. (atom_coord(j,i) > cell(j))) &
                   atom_coord(j,i) = &
                      atom_coord(j,i) - floor(atom_coord(j,i)/cell(j)) * cell(j)
             end do
          end do
          call io_close(lun)
       else
          if(iprint_init>2) write(io_lun,*) 'Entering read_atomic_positions'
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
             end if
             ! Wrap coordinates
             cell(1) = r_super_x
             cell(2) = r_super_y
             cell(3) = r_super_z
             do j = 1, 3
                if ((atom_coord(j,i) < zero) .or. (atom_coord(j,i) > cell(j))) &
                   atom_coord(j,i) = &
                      atom_coord(j,i) - floor(atom_coord(j,i)/cell(j)) * cell(j)
             end do
             flag_move_atom(1,i) = movex
             flag_move_atom(2,i) = movey
             flag_move_atom(3,i) = movez
             !          id_glob(i) = i
             if(iprint_init>0) &
                  write (io_lun,fmt='(3x, i7, 3f15.8, i3, 3L2)') &
                        i,atom_coord(1:3,i), species_glob(i), &
                        flag_move_atom(1:3,i)
          end do
          !2010.06.25 TM (Angstrom Units in coords file, but not pdb)
             if(.not.flag_fractional_atomic_coords .and. dist_units == ang) &
              atom_coord(:,:)=atom_coord(:,:)*AngToBohr
          !2010.06.25 TM (Angstrom Units in coords file, but not pdb)
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
    return
  end subroutine read_atomic_positions
  !!***

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
                              IPRINT_TIME_THRES3
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
          if (iprint_init > 2) write(io_lun,*) 'Writing write_atomic_positions'
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
          if(iprint_init>2) write(io_lun,*) 'Writing read_atomic_positions'
          call io_assign(lun)
          if(append_coords) then
             open(unit=lun,file=filename,position='append')
             write(lun,*)
          else
             open(unit=lun,file=filename)
          end if
          ! Read supercell vector - for now it must be orthorhombic so we use x and y as dummy variables
          write(lun,fmt='(3f18.10)') r_super_x, zero, zero
          write(lun,fmt='(3f18.10)') zero, r_super_y, zero
          write(lun,fmt='(3f18.10)') zero, zero, r_super_z
          write(lun,fmt='(i8)') ni_in_cell
          do i=1,ni_in_cell
             if(flag_fractional_atomic_coords) then
                write(lun,fmt='(3f18.10,i6,3L3)') atom_coord(1,i)/r_super_x,atom_coord(2,i)/r_super_y,&
                     atom_coord(3,i)/r_super_z, species_glob(i),flag_move_atom(1,i),flag_move_atom(2,i), &
                     flag_move_atom(3,i)
             else
                write(lun,fmt='(3f18.10,i6,3L3)') atom_coord(1,i),atom_coord(2,i),atom_coord(3,i), &
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
  !!  SOURCE
  !!
  subroutine read_mult(myid,parts,part_file)

    use datatypes
    use global_module,    only: numprocs, iprint_init, id_glob,       &
                                ni_in_cell, x_atom_cell, y_atom_cell, &
                                z_atom_cell, atom_coord, id_glob_inv, &
                                species_glob
    use group_module,     only: make_cc2, part_method, PYTHON, HILBERT
    use basic_types
    use dimens,           only: r_super_x, r_super_y, r_super_z
    use species_module,   only: species
    use maxima_module,    only: maxpartsproc, maxatomspart,           &
                                maxatomsproc, maxpartscell
    use construct_module, only: init_group

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

    !parts%ng_on_node = 0
    if(part_method == PYTHON) then
       if(myid==0) then
          if(iprint_init>1) then
             write(io_lun,*) 'Reading partition data'
             write(io_lun,*) '--------------------------'
             write(io_lun,*) 
             write(io_lun,*) '  *** Partition Data Starts ***'
          end if
          call read_partitions(parts,part_file)
          if(iprint_init>1) write(io_lun,*) '  *** Partition Data Ends ***'
       end if
       isendbuf = 0
       if(myid==0) then
          isendbuf(1) = parts%ngcellx
          isendbuf(2) = parts%ngcelly
          isendbuf(3) = parts%ngcellz
          isendbuf(4) = maxpartsproc
          isendbuf(5) = maxatomspart
          isendbuf(6) = maxatomsproc
          np_in_cell = parts%ngcellx*parts%ngcelly*parts%ngcellz
       endif
       call gcopy(isendbuf,6)
       if(myid/=0) then
          parts%ngcellx = isendbuf(1)
          parts%ngcelly = isendbuf(2)
          parts%ngcellz = isendbuf(3)
          maxpartsproc = isendbuf(4)
          maxatomspart = isendbuf(5)
          maxatomsproc = isendbuf(6)
          np_in_cell = parts%ngcellx*parts%ngcelly*parts%ngcellz
          maxpartscell = np_in_cell
          mx_tmp_edge = max(parts%ngcellx,parts%ngcelly,parts%ngcellz)
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
    else if (part_method == HILBERT) then
       if (iprint_init > 1.AND.myid==0) then
          write(io_lun,*) 'Partitioning using Hilbert curves'
          write(io_lun,*) '---------------------------------'
          write(io_lun,*)
       end if
       call create_sfc_partitions(myid, parts)
       np_in_cell = parts%ngcellx*parts%ngcelly*parts%ngcellz
       if (iprint_init > 1.AND.myid==0) write(io_lun,*) 'Finished partitioning'
    endif
    ! inverse table to npnode
    do np=1,np_in_cell
       parts%inv_ngnode(parts%ngnode(np))=np
    enddo
    call make_cc2(parts,numprocs)
    do ni = 1, ni_in_cell
       id_global= id_glob(ni)
       x_atom_cell(ni) = atom_coord(1,id_global)
       y_atom_cell(ni) = atom_coord(2,id_global)
       z_atom_cell(ni) = atom_coord(3,id_global)
       species(ni)     = species_glob(id_global)
       !id_glob_inv(id_global) = ni  !in read_partitions
    enddo
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
  !!  SOURCE
  !!
  subroutine read_partitions(parts, part_file)

    ! Module usage
    use datatypes
    use global_module,    only: numprocs, iprint_init, ni_in_cell, &
                                id_glob, id_glob_inv
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
    integer :: nnode
    integer :: nnd,nnd1,np,np1,ind_part,n_cont,n_beg,ni,ni1
    integer :: irc,ierr,np_in_cell, lun, glob_count, map, ios
    integer :: ind_global, ntmpx, ntmpy, ntmpz, ntmp1, ntmp2, ntmp3, mx_tmp_edge

    if(iprint_init>2.AND.myid==0) write(io_lun,*) 'Entering read_partitions'
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
         write(io_lun,*) 'Atoms proc max: ',maxatomsproc, maxpartsproc
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
    return
  end subroutine read_partitions
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine create_sfc_partitions(myid, parts)
  ! --------------------------------------------------------------------
  
  !!****f* io_module/create_sfc_partitions
  !!
  !!  NAME
  !!   create_sfc_partitions -- creates partitions using SF curves
  !!  USAGE
  !!   create_sfc_partitions
  !!  PURPOSE
  !!   Determines the number of partitions to be created, maps the
  !!   real space onto a space-filling curve (Hilbert, at the moment),
  !!   partitions the atoms and distributes over processors
  !!  INPUTS
  !!   integer :: myid - processor id
  !!   type(group_set) :: parts - contains info on partitions
  !!  USES
  !!  
  !!  AUTHOR
  !!   Veronika Brazdova
  !!  CREATION DATE
  !!   2006/12/14
  !!  MODIFICATION HISTORY
  !!   15:05, 27/04/2007 drb & vb
  !!    Added further check for no atoms on a processor and corrected
  !!    small bug
  !!   28/01/2008 Veronika & Milica
  !!    Fixed memory allocation and deallocation
  !!   2009/07/08 16:56 dave
  !!    Formatting tweak
  !!   2011/12/08 15:47 dave
  !!    Bug fix: changed criterion for slab or line from very_small to
  !!    0.1 (following COR suggestion) Also changed many processor
  !!    scheme to NOT wrap atoms on partition boundary down one
  !!    partition
  !!  SOURCE
  !!  
  subroutine create_sfc_partitions(myid, parts)

    ! Module usage
    use datatypes
    use global_module,    only: numprocs, iprint_init, ni_in_cell, &
                                id_glob, id_glob_inv, atom_coord,  &
                                sorted_coord, global_maxatomspart, &
                                load_balance, many_processors
    use maxima_module,    only: maxpartsproc, maxatomspart,        &
                                maxatomsproc, maxpartscell
    use basic_types,      only: group_set
    use group_module,     only: make_cc2
    use numbers,          only: very_small, one, two, three
    use dimens,           only: r_super_x, r_super_y, r_super_z
    use construct_module, only: init_group
    use GenComms,         only: my_barrier, mtime

    implicit none

    ! Passed variables
    integer :: myid
    type (group_set) :: parts

    ! Local variables
    integer :: nnode, min_atoms_part, max_atoms_part
    integer :: min_parts_occ, max_parts_occ
    real(double) :: real_min_parts, real_max_parts, real_j_max, &
                    real_j_min, dims, time0, time1
    real(double), dimension(-3:3) :: minmax_coords
    real(double), dimension(1:3) :: occupied_cell, tmp_parts, hc_edge
    integer :: b, i0, i1, i2, i3, i4, i5, tmp, stat, start
    integer :: j_min, j_max, no_hc, parts_edge
    integer, dimension(ni_in_cell) :: atom_sfc_id
    integer(integ), allocatable, dimension(:) :: sfc_sequence, cc_part_id, cc_to_H
    integer, dimension(0:numprocs-1) :: parts_inode_tmp, no_atoms_proc
    logical :: refine, reassign, too_many_atoms
    logical, dimension(1:3) :: shift_start
    integer, dimension(1:3) :: coord_start, coord_end
    integer, allocatable, dimension(:,:) :: map

    ! Variables required for axes_to_transpose and transpose_to_axes
    integer, dimension(0:2) :: X
    integer, allocatable, dimension(:) :: H
    integer :: i, j, k
    integer :: Hilbert

    ! Variables related to partitions
    integer :: np_in_cell, maxatomsproc_tmp, parts_per_node, left_parts
    integer, dimension(0:numprocs-1) :: parts_boundaries
    integer :: total_no_atoms, total_no_parts, average, occupied_parts
    integer :: correct_parts
    integer :: mx_tmp_edge, counter, counter2


    ! Variables related to statistics
    ! integer :: proc_atoms_max, proc_atoms_min, mode, mode_value, median
    ! integer, allocatable, dimension(:) :: histogram
    ! integer :: av_atom_processor

    ! atom_sfc_id ... partition to which an atom belongs
    ! sfc_sequence ... sequence of partitions along the Hilbert curve
    ! parts_inode_tmp ... number of partitions on a processor
    ! no_atoms_proc ... number of atoms on a processor

    ! The Hilbert curve:
    ! X(0:2) ... real space coordinates
    ! axes_to_transpose takes X, representing the origin of a Hilbert cube
    ! and rewrites it with integers from which the Hilbert integer can be
    ! constructed - i.e., the position of the Hilbert cube along the Hilbert
    ! curve (this is NOT the Gray code of the position but the binary
    ! representation of the sequential number).

    !write(io_lun,*) 'Proc, starting sfc: ',myid
    !if(myid==0) then
    if (iprint_init > 2.AND.myid==0) &
         write(io_lun,*) 'Entering create_sfc_partitions'

    ! av_atom_part is not used at the moment 
    ! av_atom_part = 13
    min_atoms_part = 5
    if(min_atoms_part>ni_in_cell) min_atoms_part = ni_in_cell
    max_atoms_part = global_maxatomspart !20

    maxpartscell = 0 ! max number of partitions in the unit cell

    np_in_cell = 0 ! number of partitions

    min_parts_occ = ni_in_cell / max_atoms_part + 1
    max_parts_occ = ni_in_cell / min_atoms_part 

    if (max_parts_occ == 0) max_parts_occ = 1

    if (iprint_init > 3.AND.myid==0) &
         write(io_lun,*) 'ni_in_cell = ', ni_in_cell
    if (iprint_init > 3.AND.myid==0) &
         write(io_lun,*) 'min_parts_occ = ', min_parts_occ
    if (iprint_init > 3.AND.myid==0) &
         write(io_lun,*) 'max_parts_occ = ', max_parts_occ

    if (.not.many_processors) then
       ! Sort the atoms according to x, y, and z coordinates
       ! and store their *indices* in sorted_coord.

       time0 = mtime()

       allocate(sorted_coord(3,ni_in_cell), STAT = stat)
       if(stat/=0) &
            call cq_abort("Failure to allocate sorted_coord: ",ni_in_cell)

       do i = 1, ni_in_cell
          sorted_coord(1,i) = i
          sorted_coord(2,i) = i
          sorted_coord(3,i) = i
       end do

       ! The mergesort routine is the first thing to parallelise
       call mergesort(1,ni_in_cell,1)
       call mergesort(1,ni_in_cell,2)
       call mergesort(1,ni_in_cell,3)

       ! Check the volume which contains atoms:
       minmax_coords(0) = 0
       minmax_coords(1) = atom_coord(1,sorted_coord(1,ni_in_cell))
       minmax_coords(2) = atom_coord(2,sorted_coord(2,ni_in_cell))
       minmax_coords(3) = atom_coord(3,sorted_coord(3,ni_in_cell))
       minmax_coords(-1) = atom_coord(1,sorted_coord(1,1))
       minmax_coords(-2) = atom_coord(2,sorted_coord(2,1))
       minmax_coords(-3) = atom_coord(3,sorted_coord(3,1))

       if (iprint_init > 4.AND.myid==0)                   &
            write(io_lun,'(a,6f12.6)') "Minmax:",         &
            minmax_coords(1),  &
            minmax_coords(2),  &
            minmax_coords(3),  &
            minmax_coords(-1), &
            minmax_coords(-2), &
            minmax_coords(-3)
    else

       ! Check the volume which contains atoms:
       minmax_coords(1) = 0
       minmax_coords(2) = 0
       minmax_coords(3) = 0
       minmax_coords(-1) = r_super_x
       minmax_coords(-2) = r_super_y
       minmax_coords(-3) = r_super_z
       do i = 1, ni_in_cell
          if(minmax_coords(1) < atom_coord(1,i)) minmax_coords(1) = atom_coord(1,i) 
          if(minmax_coords(2) < atom_coord(2,i)) minmax_coords(2) = atom_coord(2,i) 
          if(minmax_coords(3) < atom_coord(3,i)) minmax_coords(3) = atom_coord(3,i) 
          if(minmax_coords(-1) > atom_coord(1,i)) minmax_coords(-1) = atom_coord(1,i) 
          if(minmax_coords(-2) > atom_coord(2,i)) minmax_coords(-2) = atom_coord(2,i) 
          if(minmax_coords(-3) > atom_coord(3,i)) minmax_coords(-3) = atom_coord(3,i) 
       end do

       if (iprint_init > 4.AND.myid==0) &
            write(io_lun,'(a,6f12.6)') "Minmax:", &
            minmax_coords(1),  &
            minmax_coords(2),  &
            minmax_coords(3),  &
            minmax_coords(-1), &
            minmax_coords(-2), &
            minmax_coords(-3)

    end if

    if (myid==0.AND.iprint_init>2) &
         write(io_lun,'(a,f23.10)') "Time for min-max", mtime()-time0

    do i = 1, 3
       occupied_cell(i) = minmax_coords(i) - minmax_coords(-i)
    end do

    ! Find the best level of recursion taking into account min/max_parts and
    ! size of the occupied_cell.
    ! Remember that # partitions in all three cart. directions must be the same.
    !
    ! Are all atoms in a plane or a line || to any of the main axes?

    dims = three

    if (occupied_cell(1) < 0.1_double) then!very_small) then
       tmp_parts(1) = one
       dims = dims - one
    else
       tmp_parts(1) = r_super_x / occupied_cell(1) 
    end if
    if (occupied_cell(2) < 0.1_double) then!very_small) then
       tmp_parts(2) = one
       dims = dims - one
    else
       tmp_parts(2) = r_super_y / occupied_cell(2) 
    end if
    if (occupied_cell(3) < 0.1_double) then!very_small) then
       tmp_parts(3) = one
       dims = dims - one 
    else
       tmp_parts(3) = r_super_z / occupied_cell(3) 
    end if
    if (dims < very_small) dims = one ! if dims is 0, we only have one atom

    if (iprint_init > 3.AND.myid==0) &
         write(io_lun,'(a,3f12.6)') "Tmp_parts:", tmp_parts(1), &
                                    tmp_parts(2), tmp_parts(3)
    real_min_parts = min_parts_occ*(tmp_parts(1) * tmp_parts(2) * &
                                    tmp_parts(3))**(one/dims)
    real_max_parts = max_parts_occ*(tmp_parts(1) * tmp_parts(2) * &
                                    tmp_parts(3))**(one/dims)

    if (iprint_init > 3.AND.myid==0) &
         write(io_lun,'("Min parts = ",f14.5)') real_min_parts
    if (iprint_init > 3.AND.myid==0) &
         write(io_lun,'("Max parts = ",f14.5)') real_max_parts

    ! j is the actual level of recursion, we will have 2^j^3 partitions
    real_j_max = log(real_max_parts)/log(two)/dims
    real_j_min = log(real_min_parts)/log(two)/dims
    if (iprint_init > 3.AND.myid==0) &
         write(io_lun,'("Real j max = ",f15.9)') real_j_max
    if (iprint_init > 3.AND.myid==0) &
         write(io_lun,'("Real j min = ",f15.9)') real_j_min
    j_max = floor(real_j_max)
    j_min = ceiling(real_j_min)
    if (iprint_init > 3.AND.myid==0) write(io_lun,'("J max = ",i5)') j_max
    if (iprint_init > 3.AND.myid==0) write(io_lun,'("J min = ",i5)') j_min

    ! This can happen when the interval real_j_min -- real_j_max is small
    if (j_min > j_max) then
       k = j_min
       j_min = j_max
       j_max = k
       if (iprint_init > 3.AND.myid==0) write(io_lun,'("J max = ",i5)') j_max
       if (iprint_init > 3.AND.myid==0) write(io_lun,'("J min = ",i5)') j_min
    end if
    ! Similar as above
    if (j_min == 0) j_min = 1

    ! Now j_min and j_max give the range of level of recursion:
    ! <2^j_min; 2^j_max> which is the "b" going into axes_to_transpose
    ! and transpose_to_axes

    ! checks: - too many partitions with too few / too many atoms?
    !         - overflowing integers?
    !         - close enough to the average number of atoms in a partition?

    b = j_min

    ! This is to select the HC scheme. Hard-coded at the moment, will later
    ! depend on the number of processors

    if (.not.many_processors) then

       time0 = mtime()

       hc:   do 
          refine = .false.
          maxatomspart = 0

          hc2:    do
             if (iprint_init > 3.AND.myid==0) &
                  write(io_lun,'(a,i2)') "Hilbert curve recursion level =", b

             ! Initialise the array of atom <--> Hilbert cube
             ! -1 means atom not assigned to any HC
             do tmp = 1, ni_in_cell
                atom_sfc_id(tmp) = -1
             end do

             no_hc = 2**(3*b)

             allocate(H(1:3*b), STAT=stat)
             if (stat /= 0) &
                  call cq_abort('Create_sfc_partitions: error allocating &
                                 &array H (Hilbert)')
             allocate(sfc_sequence(0:no_hc-1), STAT=stat)
             if (stat /= 0) &
                  call cq_abort('Create_sfc_partitions: error allocating &
                                 &array sfc_sequence')
             allocate(map(0:no_hc-1,1:global_maxatomspart), STAT=stat)
             if (stat /= 0) &
                  call cq_abort('Create_sfc_partitions: error allocating &
                                 &array map')
             parts_edge = 2**b 

             hc_edge(1) = r_super_x / parts_edge
             hc_edge(2) = r_super_y / parts_edge
             hc_edge(3) = r_super_z / parts_edge

             allocate(cc_part_id(0:no_hc-1), cc_to_H(no_hc), STAT=stat)
             if (stat /= 0) &
                  call cq_abort('Create_sfc_partitions: error allocating &
                                 &array cc_part_id')

             ! Take every Hilbert cube in turn
             ! Parallelise the loop(s) going over X(0) (and X(1), X(2))

             coord_start(1) = 1
             coord_end(1) = 1
             do i0 = 0, parts_edge - 1
                shift_start(1) = .true.
                coord_start(1) = coord_end(1)

                coord_start(2) = 1
                do i1 = 0, parts_edge - 1
                   shift_start(2) = .true.

                   coord_start(3) = 1
                   do i2 = 0, parts_edge - 1
                      shift_start(3) = .true.

                      X(0) = i0
                      X(1) = i1
                      X(2) = i2

                      if (iprint_init > 4.AND.myid==0) &
                           write(io_lun,'(a,3i7)') "Sent:", X(0), X(1), X(2)
                      call axes_to_transpose(X, b, 3)

                      ! Construct the Hilbert integer
                      do i = 1, b
                         do j = 1, 3
                            H(j+(i-1)*3) = iand(ishft(X(j-1),(i-b)),1)
                         end do
                      end do
                      if (iprint_init > 4.AND.myid==0) &
                           write(io_lun,'(a,12b1)') "Hilbert integer = ", H

                      ! Convert the binary Hilbert integer to decimal
                      Hilbert = 0
                      do i = 1, 3 * b
                         Hilbert = Hilbert + 2**(i-1)*H(b*3+1-i)
                      end do

                      if (iprint_init > 3.AND.myid==0) &
                           write(io_lun,'(a,3i3,i7)') "Result:", &
                           X(0), X(1), X(2), Hilbert

                      ! How many atoms are in this Hilbert cube?  A cunning
                      ! way to search in all 3 cartesian directions using
                      ! the arrays of (sorted) atom IDs

                      sfc_sequence(Hilbert) = 0
                      if (iprint_init > 4.AND.myid==0) &
                           write(io_lun,'(a,f9.5)') "X value:", (hc_edge(1) * (i0 + 1))
                      if (iprint_init > 4.AND.myid==0) &
                           write(io_lun,'(a,f9.5)') "Y value:", (hc_edge(2) * (i1 + 1))
                      if (iprint_init > 4.AND.myid==0) &
                           write(io_lun,'(a,f9.5)') "Z value:", (hc_edge(3) * (i2 + 1))

                      cc_part_id(Hilbert) = 1 + i0 * parts_edge**2 + i1 * parts_edge + i2
                      cc_to_H(1 + i0 * parts_edge**2 + i1 * parts_edge + i2) = Hilbert
                      if (iprint_init > 4.AND.myid==0) &
                           write(io_lun,'(a,i6)') "cc_part_id =", cc_part_id(Hilbert)

                      ! Check whether we are in the occupied region. If not,
                      ! sfc_sequence(Hilbert) will remain 0
                      if (((hc_edge(1) * (i0 + 1)) >= minmax_coords(-1)) .and. &
                           ( (hc_edge(1) * i0) <= minmax_coords(1))        .and. &
                           ( (hc_edge(2) * (i1 + 1)) >= minmax_coords(-2)) .and. &
                           ( (hc_edge(2) * i1) <= minmax_coords(2))        .and. &
                           ( (hc_edge(3) * (i2 + 1)) >= minmax_coords(-3)) .and. &
                           ( (hc_edge(3) * i2) <= minmax_coords(3)))  then

                         ! Loop over all atoms (potentially)
                         do i3 = coord_start(1), ni_in_cell

                            ! If x coord is > the HC x boundries, exit
                            if (iprint_init > 4.AND.myid==0) &
                                 write(io_lun,'(a,i5,a,f9.5)') "Atom ",&
                                 sorted_coord(1,i3), &
                                 " coord =", atom_coord(1,sorted_coord(1,i3))
                            if ((atom_coord(1,sorted_coord(1,i3)) > &
                                 (hc_edge(1) * (i0 + 1)))) then
                               if (iprint_init > 4.AND.myid==0) &
                                    write(io_lun,'(a,i5)') "Atom over x:", &
                                    sorted_coord(1,i3) 
                               coord_end(1) = i3
                               exit
                            end if

                            ! Are we in the right delta(x) range?
                            if ((i0 == 0 .and. (atom_coord(1,sorted_coord(1,i3)) >= hc_edge(1) * i0)) .or. &
                                ((i0 /= 0 .and. atom_coord(1,sorted_coord(1,i3)) > hc_edge(1) * i0))) then
                               if (iprint_init > 4.AND.myid==0) &
                                    write(io_lun,'(a,i5)') "  Entering y loop. Atom id: ", sorted_coord(1,i3)
                               if (shift_start(1)) then 
                                  coord_start(1) = i3
                                  shift_start(1) = .false.
                               end if

                               do i4 = 1, ni_in_cell

                                  ! If y coord is > the HC y boundries, exit
                                  if ((atom_coord(2,sorted_coord(2,i4)) > (hc_edge(2) * (i1 + 1)))) then
                                     if (iprint_init > 4.AND.myid==0) &
                                          write(io_lun,'(a,i5)') "  Atom over y:", sorted_coord(2,i4)
                                     exit
                                  end if

                                  ! Are we in the right delta(y) range? .and.
                                  ! Does the atom ID of the atom i3 appear in this range?
                                  if (sorted_coord(1,i3) == sorted_coord(2,i4)) then
                                     if ((i1 == 0 .and. atom_coord(2,sorted_coord(2,i4)) >= hc_edge(2) * i1) &
                                          .or. ((i1 /= 0 .and. atom_coord(2,sorted_coord(2,i4)) > hc_edge(2) * i1))) then
                                        if (iprint_init > 4.AND.myid==0) write(io_lun,'(a,i5)') &
                                             "    Entering z loop. Atom id: ", sorted_coord(1,i3)

                                        do i5 = 1, ni_in_cell

                                           ! If z coord is > the HC z boundries, exit
                                           if (atom_coord(3,sorted_coord(3,i5)) > (hc_edge(3) * (i2 + 1))) then
                                              if (iprint_init > 4.AND.myid==0) write(io_lun,'(a,i5)') "    Atom over z:", sorted_coord(3,i5)
                                              exit
                                           end if

                                           ! Are we in the right delta(z) range? .and.
                                           ! Does the atom ID of the atom i3 appear in this range?
                                           if (sorted_coord(1,i3) == sorted_coord(3,i5)) then
                                              if ((i2 == 0 .and. atom_coord(3,sorted_coord(3,i5)) >= hc_edge(3) * i2) &
                                                   .or. ((i2 /= 0 .and. atom_coord(3,sorted_coord(3,i5)) > hc_edge(3) * i2))) then

                                                 ! Number of atoms in this HC
                                                 sfc_sequence(Hilbert) = sfc_sequence(Hilbert) + 1
                                                 atom_sfc_id(sorted_coord(1,i3)) = Hilbert
                                                 if (sfc_sequence(Hilbert) > global_maxatomspart) then
                                                    refine = .true.
                                                    deallocate(H,STAT = stat)
                                                    if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array H (Hilbert)')
                                                    deallocate(sfc_sequence,STAT = stat)
                                                    if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array sfc_sequence')
                                                    deallocate(map,STAT = stat)
                                                    if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array sfc_sequence')
                                                    deallocate(cc_to_H,cc_part_id,STAT = stat)
                                                    if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array cc_part_id')
                                                    exit hc2
                                                 end if
                                                 map(Hilbert,sfc_sequence(Hilbert)) = sorted_coord(1,i3)
                                                 if (iprint_init > 4.AND.myid==0) &
                                                      write(io_lun,'(a,i5,a)') "    Atom", sorted_coord(1,i3), " matches"
                                                 if (iprint_init > 4.AND.myid==0) &
                                                      write(io_lun,'(a,i5,a,i5)') "    Sequence(", Hilbert,") atoms:",sfc_sequence(Hilbert)
                                                 if (maxatomspart < sfc_sequence(Hilbert)) maxatomspart = sfc_sequence(Hilbert)
                                                 exit
                                              else
                                                 if (iprint_init > 4.AND.myid==0) &
                                                      write(io_lun,'(a,i5,a)') "    Z range, atom",sorted_coord(3,i5)," does not match"
                                                 exit
                                              end if
                                           end if

                                        end do
                                        exit
                                     else
                                        if (iprint_init > 4.AND.myid==0) &
                                             write(io_lun,'(a,i5,a)') "  Y range, atom",sorted_coord(2,i4)," does not match"
                                        exit
                                     end if
                                  end if
                               end do
                            end if
                         end do
                      end if

                      if (iprint_init > 2.AND.myid==0) &
                           write(io_lun,'(a,i5,a,i5,a)') "Sfc_sequence(",Hilbert,") =",sfc_sequence(Hilbert)," finished"

                      !         if (sfc_sequence(Hilbert) > global_maxatomspart) then
                      !           refine = .true.
                      !           exit hc2
                      !         end if

                   end do
                end do
             end do

             ! Delete later
             if (iprint_init > 3.AND.myid==0) then
                do tmp = 1, ni_in_cell
                   write(io_lun,'(a,i5,a,i5)') "Atom number", tmp, " HC cube", atom_sfc_id(tmp)
                end do
             end if

             ! Parallelise
             do tmp = 1, ni_in_cell
                if (atom_sfc_id(tmp) == -1) call cq_abort('Atom not assigned to a partition! Aborting job.')
             end do

             deallocate(H,STAT = stat)
             if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array H (Hilbert)')

             if (.not.refine) exit hc

             deallocate(sfc_sequence,STAT = stat)
             if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array sfc_sequence')

             deallocate(map,STAT = stat)
             if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array sfc_sequence')

             deallocate(cc_to_H,cc_part_id,STAT = stat)
             if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array cc_part_id')

          end do hc2
          b = b + 1

       end do hc

       if(myid==0.AND.iprint_init>2) &
            write(io_lun,'(10x,a,f14.10)') "Time for assigning atoms to partitions, loop over partitions", mtime()-time0

    else ! Many Processors

       time0 = mtime()

       refine = .true.
       do 
          refine = .false.

          parts_edge = 2**b

          no_hc = 2**(3*b)

          allocate (H(1:3*b), STAT = stat)
          if (stat /= 0) call cq_abort('Create_sfc_partitions: error allocating array H (Hilbert)')
          allocate (sfc_sequence(0:no_hc-1), STAT = stat)
          if (stat /= 0) call cq_abort('Create_sfc_partitions: error allocating array sfc_sequence')
          allocate (map(0:no_hc-1,1:global_maxatomspart), STAT = stat)
          if (stat /= 0) call cq_abort('Create_sfc_partitions: error allocating array map')
          allocate (cc_part_id(0:no_hc-1), cc_to_H(no_hc), STAT = stat)
          if (stat /= 0) call cq_abort('Create_sfc_partitions: error allocating array cc_part_id')

          sfc_sequence = 0
          atom_sfc_id = -1

          hc_edge(1) = r_super_x / parts_edge
          hc_edge(2) = r_super_y / parts_edge
          hc_edge(3) = r_super_z / parts_edge

          do i = 1, ni_in_cell

             ! Assign atom to part of real space and
             ! check if it is on a boundary
             X(0) = floor(atom_coord(1,i) / hc_edge(1))
             !if ((X(0) /= 0) .and. ((atom_coord(1,i) / hc_edge(1)) == real(X(0)))) X(0) = X(0) - 1
             X(1) = floor(atom_coord(2,i) / hc_edge(2))
             !if ((X(1) /= 0) .and. ((atom_coord(2,i) / hc_edge(2)) == real(X(1)))) X(1) = X(1) - 1
             X(2) = floor(atom_coord(3,i) / hc_edge(3))
             !if ((X(2) /= 0) .and. ((atom_coord(3,i) / hc_edge(3)) == real(X(2)))) X(2) = X(2) - 1
             ! Get the Hilbert coordinates
             call axes_to_transpose(X, b, 3)

             ! Construct the Hilbert integer
             do k = 1, b
                do j = 1, 3
                   H(j+(k-1)*3) = iand(ishft(X(j-1),(k-b)),1)
                end do
             end do
             if (iprint_init > 4.AND.myid==0) write(io_lun,'(a,12b1)') "Hilbert integer = ", H

             ! Convert the binary Hilbert integer to decimal
             Hilbert = 0
             do k = 1, 3 * b
                Hilbert = Hilbert + 2**(k-1)*H(b*3+1-k)
             end do

             atom_sfc_id(i) = Hilbert
             sfc_sequence(Hilbert) = sfc_sequence(Hilbert) + 1
             if(sfc_sequence(Hilbert)>global_maxatomspart) then
                refine = .true. ! We use the refine switch to flag early exit
                exit
             end if
             map(Hilbert,sfc_sequence(Hilbert)) = i

             if (iprint_init > 3.AND.myid==0) write(io_lun,'(a,i5,a,i5)') "Atom number", i, " HC cube", Hilbert
             if (iprint_init > 3.AND.myid==0) write(io_lun,'(a,3i3,i7)') "Result:", X(0), X(1), X(2), Hilbert

          end do
          if(.not.refine) then
             do tmp = 1, ni_in_cell
                if (atom_sfc_id(tmp) == -1) call cq_abort('Atom not assigned to a partition! Aborting job.')
             end do

             do i = 0, no_hc - 1
                if (iprint_init > 3.AND.myid==0) write(io_lun,'(a,i5,a,i5,a)') "Sfc_sequence(",i,") =",sfc_sequence(i)," finished"
                if (sfc_sequence(i) > global_maxatomspart) then
                   refine = .true.
                   exit 
                end if
             end do
          end if
          if (.not.refine) exit

          deallocate(H,STAT = stat)
          if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array H (Hilbert)')

          deallocate(sfc_sequence,STAT = stat)
          if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array sfc_sequence')

          deallocate(map,STAT = stat)
          if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array sfc_sequence')

          deallocate(cc_part_id, cc_to_H, STAT = stat)
          if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array cc_part_id')

          b = b + 1

       end do

       maxatomspart = 0

       do i0 = 0, parts_edge - 1
          do i1 = 0, parts_edge - 1
             do i2 = 0, parts_edge - 1
                X(0) = i0
                X(1) = i1
                X(2) = i2
                call axes_to_transpose(X, b, 3)
                ! Construct the Hilbert integer
                do k = 1, b
                   do j = 1, 3
                      H(j+(k-1)*3) = iand(ishft(X(j-1),(k-b)),1)
                   end do
                end do
                if (iprint_init > 4.AND.myid==0) write(io_lun,'(a,12b1)') "Hilbert integer = ", H

                ! Convert the binary Hilbert integer to decimal
                Hilbert = 0
                do k = 1, 3 * b
                   Hilbert = Hilbert + 2**(k-1)*H(b*3+1-k)
                end do
                cc_part_id(Hilbert) = 1 + i0 * parts_edge**2 + i1 * parts_edge + i2
                cc_to_H(1 + i0 * parts_edge**2 + i1 * parts_edge + i2) = Hilbert
                if (maxatomspart < sfc_sequence(Hilbert)) maxatomspart = sfc_sequence(Hilbert)
                if (iprint_init > 4.AND.myid==0) write(io_lun,'(a,i6)') "cc_part_id =", cc_part_id(Hilbert)
             end do
          end do
       end do

       deallocate(H,STAT = stat)
       if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array H (Hilbert)')

       if(myid==0.AND.iprint_init>2) &
            write(io_lun,'(10x,a,f14.10)') "Time for assigning atoms to partitions, loop over partitions", mtime()-time0
    end if ! partitioning algorithm

    maxatomsproc_tmp = ni_in_cell / numprocs 
    if (maxatomsproc_tmp == 0) maxatomsproc_tmp = 1
    ! i.e., if we have too few atoms for the number of processors

    reassign = .true. ! Continue reassigning?
    too_many_atoms = .false. ! Too many atoms on the last processor?

    if (load_balance == 0) then ! balancing number of atoms per processor
       do while (reassign)

          no_atoms_proc = 0
          parts_inode_tmp = 0

          start = 1

          no_atoms_proc(0) = sfc_sequence(0)
          parts_inode_tmp(0) = 1

          parts_boundaries = -1
          parts_boundaries(0) = 0

          do i = 0, numprocs - 1
             do while (start < no_hc)

                if (no_atoms_proc(i) >= maxatomsproc_tmp .and. i < numprocs - 1) then
                   parts_boundaries(i+1) = start
                   exit
                end if

                no_atoms_proc(i) = no_atoms_proc(i) + sfc_sequence(start)
                parts_inode_tmp(i) = parts_inode_tmp(i) + 1

                start = start + 1

             end do
          end do

          reassign = .false.

          if (iprint_init > 4.AND.myid==0) then
             write(io_lun,*)
             write(io_lun,*) "Maxatomsproc_tmp =", maxatomsproc_tmp
             do i = 0, numprocs -1
                write(io_lun,'(a,i7,a,i4,a,i5)') "Processor",i, "  Atoms", no_atoms_proc(i), "  Partitions", parts_inode_tmp(i)
             end do
          end if

          if (too_many_atoms) exit

          ! No atoms on any processor?
          do i = 0, numprocs - 1
             if ((parts_inode_tmp(i) == 0) .and. (.not.too_many_atoms)) then
                maxatomsproc_tmp = maxatomsproc_tmp - 1
                reassign = .true.
                exit
             end if
          end do

          if (numprocs > 1) then
             ! Too few atoms on the last processor?
             average = (ni_in_cell - no_atoms_proc(numprocs-1)) / (numprocs - 1)
             if (iprint_init > 4.AND.myid==0) write(io_lun,*) average, (maxatomsproc_tmp - no_atoms_proc(numprocs-1))
             if ((average - no_atoms_proc(numprocs-1) > average - maxatomsproc_tmp + 1) &
                  .and. (.not.too_many_atoms) .and. (.not.reassign)) then
                maxatomsproc_tmp = maxatomsproc_tmp - 1
                reassign = .true.
             end if

             ! Too many atoms on the last processor?
             if (no_atoms_proc(numprocs-1) - average > 0) then
                maxatomsproc_tmp = maxatomsproc_tmp + 1
                too_many_atoms = .true.
                reassign = .true.
                occupied_parts =  0
                do i = 0, no_hc - 1
                   if (sfc_sequence(i) /= 0) occupied_parts = occupied_parts + 1
                end do
                correct_parts = ni_in_cell / occupied_parts
                if (iprint_init > 4.AND.myid==0) write(io_lun,*) "correct_parts =", correct_parts

                do i = numprocs - 1, 1, -1
                   if (no_atoms_proc(i) - no_atoms_proc(i-1) > correct_parts) then
                      tmp = no_atoms_proc(i-1)
                      do while (no_atoms_proc(i) - tmp > correct_parts .and. parts_boundaries(i) < no_hc)

                         parts_boundaries(i) = parts_boundaries(i) + 1
                         no_atoms_proc(i) = no_atoms_proc(i) - sfc_sequence(parts_boundaries(i)-1)
                         no_atoms_proc(i-1) = no_atoms_proc(i-1) + sfc_sequence(parts_boundaries(i)-1)
                         parts_inode_tmp(i) = parts_inode_tmp(i) - 1
                         parts_inode_tmp(i-1) = parts_inode_tmp(i-1) + 1

                      end do
                   else
                      exit
                   end if
                end do

                if (iprint_init > 4.AND.myid==0) write(io_lun,*) "Numprocs =", numprocs
                if (iprint_init > 4.AND.myid==0) write(io_lun,*) "Loop: Maxatomsproc_tmp =", maxatomsproc_tmp
                do i = 0, numprocs - 1
                   if (iprint_init > 4.AND.myid==0) &
                        write(io_lun,'(a,i7,a,i4,a,i5)') "Processor",i, "  Atoms", no_atoms_proc(i), "  Partitions", parts_inode_tmp(i)
                   if (no_atoms_proc(i) == 0) write(io_lun,'(a,i7,a)') "WARNING: Processor", i," has no atoms! Too many processors?"
                   ! Statistics
                   ! proc_atoms_min = ni_in_cell
                   ! proc_atoms_max = 0
                   ! if (no_atoms_proc(i) < proc_atoms_min) proc_atoms_min = no_atoms_proc(i)
                   ! if (no_atoms_proc(i) > proc_atoms_max) proc_atoms_max = no_atoms_proc(i)
                end do

                exit

             end if

             if ((maxatomsproc_tmp <= 0) .and. (.not.too_many_atoms)) then
                write(io_lun,'(a,i6,a)') "WARNING: You are using too many processors for", ni_in_cell, " atoms"
                exit
             end if

          end if

       end do
       do i = 0, numprocs - 1
          if (iprint_init > 4.AND.myid==0) &
               write(io_lun,'(a,i7,a,i4,a,i5)') "Processor",i, "  Atoms", no_atoms_proc(i), "  Partitions", parts_inode_tmp(i)
          ! The following needs to stay here because the reshuffling part relies on all processors having atoms
          if (no_atoms_proc(i) == 0) then
             write(io_lun,'(a,i7,a)') "WARNING: Processor", i," has no atoms! Too many processors?"
             call cq_abort("Automatic partitioner requires all processors to have atoms")
          end if
       end do

       ! Reshuffling empty partitions

       if (numprocs > 2) then

          ! We skip the first and last processor
          do i = 1, numprocs - 2

             if ((sfc_sequence(parts_boundaries(i)) == 0) .and. (parts_inode_tmp(i) > parts_inode_tmp(i-1))) then

                j = 0
                do while ( (sfc_sequence(parts_boundaries(i)+j) == 0) .and. (parts_boundaries(i)+j < parts_boundaries(i+1)) )
                   j = j + 1
                end do
                k = min (((parts_inode_tmp(i) - parts_inode_tmp(i-1)) / 2), j)
                parts_boundaries(i) = parts_boundaries(i) + k
                parts_inode_tmp(i-1) = parts_inode_tmp(i-1) + k
                parts_inode_tmp(i) = parts_inode_tmp(i) - k

             end if

          end do

       end if

       if (numprocs > 1) then

          ! The same for the last processor
          if ((sfc_sequence(parts_boundaries(numprocs-1)) == 0) .and. &
               (parts_inode_tmp(numprocs-1) > parts_inode_tmp(numprocs-2))) then

             j = 0
             do while ( (sfc_sequence(parts_boundaries(numprocs-1)+j) == 0) .and. (parts_boundaries(numprocs-1)+j < no_hc - 1) ) 
                j = j + 1
             end do
             k = min (((parts_inode_tmp(numprocs-1) - parts_inode_tmp(numprocs-2)) / 2), j)
             parts_boundaries(numprocs-1) = parts_boundaries(numprocs-1) + k
             parts_inode_tmp(numprocs-2) = parts_inode_tmp(numprocs-2) + k
             parts_inode_tmp(numprocs-1) = parts_inode_tmp(numprocs-1) - k

          end if


          do i = 0, numprocs - 1
             if (iprint_init > 1.AND.myid==0) &
                  write(io_lun,'(a,i7,a,i4,a,i5)') "Processor",i, "  Atoms", no_atoms_proc(i), "  Partitions", parts_inode_tmp(i)
             ! Probably unnecessary - was done earlier
             if (no_atoms_proc(i) == 0) then
                write(io_lun,'(a,i7,a)') "WARNING: Processor", i," has no atoms! Too many processors?"
                call cq_abort("Automatic partitioner requires all processors to have atoms")
             end if
          end do

       end if

       ! Statistics
       ! allocate(histogram(proc_atoms_min:proc_atoms_max), STAT = stat)
       ! if(stat/=0) call cq_abort("Failure to allocate histogram")
       ! histogram = 0
       ! do i = 0, numprocs - 1
       !   histogram(no_atoms_proc(i)) = histogram(no_atoms_proc(i)) + 1
       ! end do
       ! av_atom_processor = ni_in_cell / numprocs
       ! median = (proc_atoms_min + proc_atoms_max) / 2
       ! mode_value = histogram(proc_atoms_min)
       ! mode = proc_atoms_min
       ! do i = histogram(proc_atoms_min), histogram(proc_atoms_max)
       !   if (histogram(i) > mode_value) mode = i
       ! end do
       ! write(io_lun,'(a,i6,a)') "Average =", av_atom_processor, " atoms / processor"
       ! write(io_lun,'(a,i6,a)') "Mode =", mode_value, " atoms / processor"
       ! do i = histogram(proc_atoms_min), histogram(proc_atoms_max)
       !   if (histogram(i) == mode_value) write(io_lun,'(a,i6)') "Modal value: processor", i
       ! end do
       ! deallocate(histogram, STAT = stat)
       ! if(stat/=0) call cq_abort("Failure to deallocate histogram")

    else ! balancing number of partitions per processor
       ! can easily be parallelised

       no_atoms_proc = 0
       parts_inode_tmp = 0

       parts_boundaries = -1
       parts_boundaries(0) = 0

       parts_per_node = no_hc / numprocs
       left_parts = mod(no_hc, numprocs)

       do i = 0, numprocs - 2     

          if (i < left_parts) then
             parts_inode_tmp(i) = parts_per_node + 1
             parts_boundaries(i+1) = i * parts_per_node + 2
          else
             parts_inode_tmp(i) = parts_per_node 
             parts_boundaries(i+1) = i * parts_per_node + 1
          end if
          do j = parts_boundaries(i), parts_boundaries(i+1) - 1
             no_atoms_proc(i) = no_atoms_proc(i) + sfc_sequence(j)
          end do

       end do

       parts_boundaries(numprocs-1) = (numprocs - 1) * parts_per_node
       parts_inode_tmp(numprocs-1) = parts_per_node
       do j = parts_boundaries(numprocs-1), no_hc - 1
          no_atoms_proc(numprocs-1) = no_atoms_proc(numprocs-1) + sfc_sequence(j)
       end do

    end if

    total_no_atoms = 0
    total_no_parts = 0

    do i = 0, no_hc - 1
       total_no_atoms = total_no_atoms + sfc_sequence(i)
    end do

    ! If this ever happens, we're in big trouble
    if (total_no_atoms /= ni_in_cell) call cq_abort('Creafe_sfc_partitions: Atom number does not match')

    do i = 0, numprocs - 1
       total_no_parts = total_no_parts + parts_inode_tmp(i)
    end do

    ! Ditto
    if (total_no_parts /= no_hc) call cq_abort('Creafe_sfc_partitions: Partition number does not match')

    if (iprint_init > 1.AND.myid==0) then
       do i = 0, numprocs -1
          write(io_lun,'(a,i7,a,i4,a,i5)') "Processor",i, "  Atoms", no_atoms_proc(i), "  Partitions", parts_inode_tmp(i)
       end do
    end if

    ! Assign partitions to processors

    mx_tmp_edge = parts_edge
    np_in_cell = no_hc
    maxpartscell = np_in_cell

    maxatomsproc = 0
    maxpartsproc = 0
    do i = 0, numprocs - 1
       if (maxatomsproc < no_atoms_proc(i)) maxatomsproc = no_atoms_proc(i)
       if (maxpartsproc < parts_inode_tmp(i)) maxpartsproc = parts_inode_tmp(i)
    end do

    call init_group(parts, maxpartsproc, mx_tmp_edge, np_in_cell, maxatomspart, numprocs)
    parts%ngcellx = parts_edge
    parts%ngcelly = parts_edge
    parts%ngcellz = parts_edge
    if(max(parts%ngcellx,parts%ngcelly,parts%ngcellz)>parts%mx_gedge) then
       call cq_abort('read_partitions: too many parts edge ', &
            max(parts%ngcellx,parts%ngcelly,parts%ngcellz),parts%mx_gedge)
    endif
    if (iprint_init > 0.AND.myid==0) write(io_lun,'(2x,a,3i5)') "partitions along cell sides", &
         parts%ngcellx, parts%ngcelly, parts%ngcellz
    if (iprint_init > 2.AND.myid==0) write(io_lun,'(2x,a,i5)') "maximum number of atoms in a partition", maxatomspart
    if (iprint_init > 2.AND.myid==0) write(io_lun,'(2x,a,i5)') "maximum number of atoms on a processor", maxatomsproc
    if (iprint_init > 2.AND.myid==0) write(io_lun,'(2x,a,i5)') "maximum number of partitions on a processor", maxpartsproc

    counter = 1
    do i = 0, numprocs - 1  
       parts%ng_on_node(i+1) = parts_inode_tmp(i)
       parts%inode_beg(i+1) = counter
       counter = counter + parts_inode_tmp(i)
       if(parts%inode_beg(i+1)+parts%ng_on_node(i+1)-1>parts%mx_gcell) then
          call cq_abort('read_partitions: too many parts in cell ',&
               parts%inode_beg(i+1)+parts%ng_on_node(i+1)-1, parts%mx_gcell)
       endif
       !      parts%inode_beg(i+1) = parts_boundaries(i) + 1
    end do
    !%%!    ! Loop over partitions in CC order
    !%%!    counter = 1
    !%%!    do i=1,no_hc
    !%%!       parts%ngnode(cc_to_H(i)+1) = i
    !%%!       parts%nm_group(i) = sfc_sequence(cc_to_H(i))
    !%%!       parts%icell_beg(i) = counter
    !%%!       counter = counter + sfc_sequence(cc_to_H(i))
    !%%!       if(sfc_sequence(cc_to_H(i))>0) then
    !%%!          do j=1,sfc_sequence(cc_to_H(i))
    !%%!             write(io_lun,'(a,6i3)') 'Part, seq, map: ',i,j,map(cc_to_H(i),j),cc_to_H(i),parts%icell_beg(i),&
    !%%!                  parts%icell_beg(i)+j-1
    !%%!             id_glob(parts%icell_beg(i)+j-1) = map(cc_to_H(i),j)
    !%%!             id_glob_inv(map(cc_to_H(i),j)) = parts%icell_beg(i)+j-1
    !%%!          end do
    !%%!       end if
    !%%!    end do
    !%%!    if(myid==0) then
    !%%!     counter = 1
    !%%!     do i = 0, no_hc - 1
    !%%!        write(io_lun,*) 'ID: ',parts%ngnode(i+1),cc_part_id(i)
    !%%!        write(io_lun,*) 'In group: ',parts%nm_group(cc_part_id(i)),sfc_sequence(i)
    !%%!        write(io_lun,*) 'Start: ',parts%icell_beg(cc_part_id(i)),counter
    !%%!        counter = counter + sfc_sequence(i)
    !%%!     end do
    !%%!    end if
    !%%!    call my_barrier
    !%%!    stop

    counter = 1
    do i = 0, no_hc - 1
       parts%ngnode(i+1) = cc_part_id(i)
       parts%nm_group(cc_part_id(i)) = sfc_sequence(i)
       parts%icell_beg(cc_part_id(i)) = counter
       counter = counter + sfc_sequence(i)
       if(sfc_sequence(i) > 0) then
          do j = 1, sfc_sequence(i)
             id_glob(parts%icell_beg(cc_part_id(i))+j-1) = map(i,j)
             id_glob_inv(map(i,j)) = parts%icell_beg(cc_part_id(i)) + j - 1
          end do
       end if
    end do
    !    if(myid==0) then
    !       open(unit=76,file='make_prtSFC7.dat')
    !       write(io_lun,*) 'Partitions '
    !       write(io_lun,*) parts_edge,parts_edge,parts_edge
    !       write(io_lun,*) numprocs
    !       do i=1,numprocs
    !          write(io_lun,*) i,parts%ng_on_node(i),parts%inode_beg(i)
    !          do j=1,parts%ng_on_node(i)
    !             write(io_lun,*) j,parts%ngnode(parts%inode_beg(i)+j-1), &
    !                  parts%nm_group(parts%ngnode(parts%inode_beg(i)+j-1)), &
    !                  parts%icell_beg(parts%ngnode(parts%inode_beg(i)+j-1))
    !             do k=1,parts%nm_group(parts%ngnode(parts%inode_beg(i)+j-1))
    !                write(io_lun,*) k,id_glob(parts%icell_beg(parts%ngnode(parts%inode_beg(i)+j-1))+k-1)
    !             end do
    !          end do
    !       end do
    !       close(unit=76)
    !    endif
    !    call my_barrier
    !stop


    !%%!    counter = 1
    !%%!    k=1
    !%%!    do i=1,numprocs
    !%%!       do j=1,parts%ng_on_node(i)
    !%%!          parts%icell_beg(k) = counter
    !%%!          k=k+1
    !%%!          if(k<=np_in_cell) counter = counter + parts%nm_group(k)
    !%%!       end do
    !%%!    end do
    !%%!    counter = 1
    !%%!    k=1
    !%%!    do i = 0, no_hc - 1
    !%%!      !parts%ngnode(i+1) = cc_part_id(i)
    !%%!      !parts%nm_group(cc_part_id(i)) = sfc_sequence(i)
    !%%!      !parts%icell_beg(cc_part_id(i)) = counter
    !%%!      !counter = counter + sfc_sequence(i)
    !%%!      if(sfc_sequence(i) > 0) then
    !%%!        do j = 1, sfc_sequence(i)
    !%%!          id_glob(parts%icell_beg(cc_part_id(i))+j-1) = map(i,j)
    !%%!          id_glob_inv(map(i,j)) = parts%icell_beg(cc_part_id(i)) + j - 1
    !%%!          write(io_lun,'(a,6i3)') 'Part, seq, map: ',i,j,map(i,j),cc_part_id(i), &
    !%%!               parts%icell_beg(cc_part_id(i)),parts%icell_beg(cc_part_id(i))+j-1
    !%%!        end do
    !%%!      end if
    !%%!    end do

    !end if
    !write(io_lun,*) 'Proc, ending sfc: ',myid
    call my_barrier()

    !call gcopy(maxpartsproc)
    !call gcopy(mx_tmp_edge)
    !call gcopy(np_in_cell)
    !call gcopy(maxatomspart)
    !if(myid/=0) call init_group(parts, maxpartsproc, mx_tmp_edge, np_in_cell, maxatomspart, numprocs)
    !call gcopy(parts%ngcellx)
    !call gcopy(parts%ngcelly)
    !call gcopy(parts%ngcellz)
    !call gcopy(parts%ng_on_node,numprocs)
    !call gcopy(parts%inode_beg,numprocs)
    !call gcopy(parts%ngnode,parts%mx_gcell)
    !call gcopy(parts%nm_group,parts%mx_gcell)
    !call gcopy(parts%icell_beg,parts%mx_gcell)
    !call gcopy(id_glob, ni_in_cell)
    !call gcopy(id_glob_inv, ni_in_cell)

    ! Moved up to read_mult DRB 2007/03/27
    !%%!! inverse table to npnode
    !%%!do i=1,parts%ngcellx*parts%ngcelly*parts%ngcellz
    !%%!   parts%inv_ngnode(parts%ngnode(i))=i
    !%%!enddo
    !%%!call make_cc2(parts,numprocs)

    if (.not.many_processors) then
       deallocate(sorted_coord,STAT = stat)
       if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array sorted_coord')
    end if

    deallocate(sfc_sequence,STAT = stat)
    if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array sfc_sequence')

    deallocate(map,STAT = stat)
    if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array map')

    deallocate(cc_part_id, cc_to_H, STAT = stat)
    if (stat /= 0) call cq_abort('Create_sfc_partitions: error deallocating array cc_part_id')

  end subroutine create_sfc_partitions
  !!***

  ! --------------------------------------------------------------------
  ! Subroutine axes_to_transpose(X, b, nn)
  ! --------------------------------------------------------------------
  
  !!****f* io_module/axes_to_transpose
  !!
  !!  NAME
  !!   axes_to_transpose
  !!  USAGE
  !!   axes_to_transpose(coordinates, bits, dimension)
  !!  PURPOSE
  !!   Mapping a point in 3D space onto a (1D) Hilbert curve
  !!   Based on: J. Skilling, Programming the Hilbert Curve, Bayesian
  !!   Inference and Maximum Entropy Methods in Science and Engineering:
  !!   23rd International Workshop on Bayesian Inference and Maximum
  !!   Entropy Methods in Science and Engineering. AIP Conference
  !!   Proceedings, Volume 707, pp. 381-387 (2004)
  !!  
  !!  INPUTS
  !!  
  !!  USES
  !!  
  !!  AUTHOR
  !!   Veronika Brazdova
  !!  CREATION DATE
  !!   2006/12/14
  !!  MODIFICATION HISTORY
  !!  
  !!  SOURCE
  !!  
  subroutine axes_to_transpose(X, b, nn) ! position, #bits, dimension

    implicit none

    integer, dimension(0:2) :: X
    integer M, P, Q, t, b, nn
    integer i


    if (b <= 0 .or. nn <= 0) then
       write (*,*) "Error: b and n must be positive."
       return
    end if

    M = ishft (1, b - 1)

    ! Inverse undo
    Q = M
    do while (Q > 1)
       P = Q - 1
       do i = 0, nn - 1
          ! if(iand()) will NOT work in Fortran
          ! when the result is > 1 (it does in C)
          if ((iand (X(i), Q)) > 0) then ! invert
             X(0) = ieor (X(0), P)
          else
             t = iand (ieor (X(0), X(i)), P)
             X(0) = ieor (X(0), t)
             X(i) = ieor (X(i), t)
          end if
       end do
       Q = ishft (Q, -1)
    end do ! exchange

    ! Gray encode
    do i = 1, nn - 1, 1
       X(i) = ieor (X(i), X(i-1))
    end do
    t = 0
    Q = M
    do while (Q > 1)
       if ((iand (X(nn-1), Q)) > 0 ) t = ieor (t, Q-1)
       Q = ishft (Q, -1)
    end do
    do i = 0, nn - 1
       X(i) = ieor (X(i), t)
    end do

  end subroutine axes_to_transpose
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine transpose_to_axes(X, b, nn)
  ! --------------------------------------------------------------------
  
  !!****f* io_module/transpose_to_axes
  !!
  !!  NAME
  !!   transpose_to_axes
  !!  USAGE
  !!   transpose_to_axes(coordinates, bits, dimension)
  !!  PURPOSE
  !!   Mapping a point from the Hilbert curve onto cartesian coordinates
  !!   Based on: J. Skilling, Programming the Hilbert Curve, Bayesian
  !!   Inference and Maximum Entropy Methods in Science and Engineering:
  !!   23rd International Workshop on Bayesian Inference and Maximum
  !!   Entropy Methods in Science and Engineering. AIP Conference
  !!   Proceedings, Volume 707, pp. 381-387 (2004)
  !!  
  !!  INPUTS
  !!  
  !!  USES
  !!  
  !!  AUTHOR
  !!   Veronika Brazdova
  !!  CREATION DATE
  !!   2006/12/14
  !!  MODIFICATION HISTORY
  !!  
  !!  SOURCE
  !!  
  subroutine transpose_to_axes(X, b, nn) ! position, #bits, dimension

    implicit none

    integer, dimension(0:2) :: X
    integer N, P, Q, t, b, nn
    integer i

    if (b <= 0 .or. nn <= 0) then
       write(io_lun,*) "Error: b and n must be positive."
       return
    end if

    N = ishft(2, b - 1)

    ! Gray decode by H .ieor. H/2
    t = ishft(X(nn-1), -1)

    do i = nn-1, 1, -1
       X(i) = ieor(X(i), X(i-1))
    end do

    X(0) = ieor(X(0), t)
    ! Undo excess work
    Q = 2
    do
       if (Q == N) exit

       P = Q - 1
       do i = nn - 1, 0, -1
          if ( iand(X(i),Q) > 0) then
             X(0) = ieor(X(0), P)
          else
             t = iand( ieor( X(0), X(i)), P) ! invert
             X(0) = ieor( X(0), t)
             X(i) = ieor( X(i), t)
          end if
       end do
       Q = ishft(Q, 1)
    end do ! exchange

  end subroutine transpose_to_axes
  !!***


  ! -------------------------------------------------------------------- 
  ! Subroutine mergesort
  ! -------------------------------------------------------------------- 
  
  !!****f* io_module/mergesort *
  !!
  !!  NAME 
  !!   mergesort - sorts the array atom_coord
  !!  USAGE
  !!   mergesort(left_bound, right_bound, array_index)
  !!  PURPOSE
  !!   Sorts atom_coord using the mergesort algorithm and puts 
  !!   *the indices* of the sorted array into sorted_coord(3,ni_in_cell).
  !!   the first index is passed to the routine by array_index
  !!  INPUTS
  !!  
  !!  USES
  !!   
  !!  AUTHOR
  !!   V. Brazdova
  !!  CREATION DATE
  !!   15/01/2007 Veronika
  !!  MODIFICATION HISTORY
  !!   
  !!  SOURCE
  !!    
  recursive subroutine mergesort(low, high, a_i)

    use datatypes
    use global_module, only: ni_in_cell, atom_coord, sorted_coord

    implicit none

    integer :: low, high, middle
    integer :: i 
    integer :: a_i

    if (low < high) then
       ! To avoid integer overflow:
       middle = low + ((high - low ) / 2);
       call mergesort(low, middle, a_i)
       call mergesort(middle+1, high, a_i)
       call merge_lists(low, middle+1, high, a_i)
    else
       return
    end if

  contains

    subroutine merge_lists(low, middle, high, a_i)

      use datatypes
      use global_module, only: ni_in_cell, atom_coord, sorted_coord

      implicit none

      integer :: low, high, middle
      integer :: l_position, r_position
      integer :: i, j, a_i, index1, index2

      integer, dimension(ni_in_cell) :: temp_array

      l_position = low
      r_position = middle
      i = low

      do j = 1, ni_in_cell
         temp_array(j) = sorted_coord(a_i,j)
      end do

      do while ((l_position < middle) .and. (r_position <= high))
         index1 = sorted_coord(a_i,l_position)
         index2 = sorted_coord(a_i,r_position)
         if (atom_coord(a_i,index1) <= atom_coord(a_i,index2)) then
            temp_array(i) = index1
            l_position = l_position + 1
         else
            temp_array(i) = index2
            r_position = r_position + 1
         end if
         i = i + 1
      end do

      do while (l_position < middle)
         temp_array(i) = sorted_coord(a_i,l_position)
         l_position = l_position + 1
         i = i + 1
      end do

      do while (r_position <= high)
         temp_array(i) = sorted_coord(a_i,r_position)
         r_position = r_position + 1
         i = i + 1
      end do

      do j = 1, ni_in_cell
         sorted_coord(a_i,j) = temp_array(j)
      end do

    end subroutine merge_lists

  end subroutine mergesort
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
    use block_module, only : in_block_x, in_block_y, in_block_z, n_pts_in_block, set_block_module
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

    if(inode==ionode.AND.iprint_init>2) write(io_lun,*) 'Entering read_blocks ', in_block_x,in_block_y,in_block_z
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
    call set_block_module
    if(inode==ionode) then
       call io_assign(lun)
       def = 'make_blk.dat'
       blk_coord_file = fdf_string(80,'IO.Blocks',def)
       open(unit=lun,file=blk_coord_file)
       ! Read and check blocks along cell sides
       read(lun,*) ntmpx, ntmpy, ntmpz
       if(ntmpx/=n_block_x.OR.ntmpy/=n_block_y.OR.ntmpz/=n_block_z) then
          write(io_lun,fmt='(2x,"In input, ReadBlocks T has been set, so the distribution of blocks is ")')
          write(io_lun,fmt='(2x,"read from a file.  There is an error specifying numbers of blocks.  ")')
          write(io_lun,fmt='(2x,"Found from the file: ",3i6)') ntmpx, ntmpy, ntmpz
          write(io_lun,fmt='(2x,"Calculated from input using grid points: ",3i6)') n_block_x, n_block_y, n_block_z
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
    !use set_blipgrid_module, only: naba_atm, supp

    ! Passed variables
    integer :: size, inode
    real(double), dimension(size) :: density
    integer, optional :: spin

    ! Local variables
    integer :: lun, block, n_point, n_i, n
    character(len=20) :: filename

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
    !   do n_i=1, naba_atm(supp)%no_of_atom(block)*NSF*n_pts_in_block, n_pts_in_block
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


  !!****f* io_module/dump_charge2
  !! PURPOSE
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   David Bowler
  !! CREATION DATE 
  !!   
  !! MODIFICATION HISTORY
  !!   2011/08/30 L.Tong
  !!     Added spin polarisation
  !!     - optional integer parameter spin: 0 = total, 1 = up, 2 = dn
  !!     - if spin is not present, then assume spin non-polarised
  !!       calculation
  !! SOURCE
  !!
  subroutine dump_charge2(stub, density, size, inode, spin)

    use datatypes
    use block_module, only: n_pts_in_block
    use primary_module, only: domain
    use global_module, only: numprocs
    !use set_blipgrid_module, only: naba_atm, supp

    ! Passed variables
    integer :: size, inode
    real(double), dimension(size) :: density
    character(len=*) :: stub
    integer, optional :: spin 

    ! Local variables
    integer :: lun, block, n_point, n_i, n
    character(len=18) :: filename

    ! Build a filename based on node number
    if (present (spin)) then
       select case (spin)
       case (0)
          call get_file_name (stub//'den', numprocs, inode, filename)
       case (1)
          call get_file_name (stub//'den_up', numprocs, inode, filename)
       case (2)
          call get_file_name (stub//'den_dn', numprocs, inode, filename)
       end select
    else
       call get_file_name (stub//'den', numprocs, inode, filename)
    end if
    ! Open file
    call io_assign (lun)
    open (unit = lun, file = filename)
    ! Dump charge density
    ! do block=1, domain%groups_on_node
    !   n_point = (block - 1) * n_pts_in_block
    !   do n_i=1, naba_atm(supp)%no_of_atom(block) * NSF * n_pts_in_block, &
    !        n_pts_in_block
    !      do n=1, n_pts_in_block
    !         write(unit=lun,fmt='(f30.15)') density(n_point+n)
    !      end do
    !   end do
    ! end do
    do n=1, domain%groups_on_node * n_pts_in_block
       write (unit=lun, fmt='(g13.6)') density(n)
    end do
    call io_close (lun)
    return
  end subroutine dump_charge2


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
  !!   2011/11/27 L.Tong
  !!     Added spin polarisation
  !!  SOURCE
  !!
  subroutine grab_charge(density, size, inode, spin)

    use datatypes
    use block_module,   only: n_pts_in_block
    use primary_module, only: domain
    use global_module,  only: numprocs
    !use set_blipgrid_module, only: naba_atm, supp

    ! Passed variables
    integer :: size, inode
    real(double), dimension(size) :: density
    integer, optional :: spin

    ! Local variables
    integer           :: lun, block, n_point, n_i, n
    character(len=10) :: filename

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
    open(unit=lun, file=filename)
    ! Grab charge density
    !do block=1, domain%groups_on_node
    !   n_point = (block - 1) * n_pts_in_block
    !   do n_i=1, naba_atm(supp)%no_of_atom(block)*NSF*n_pts_in_block, n_pts_in_block
    !      do n=1, n_pts_in_block
    !         read(unit=lun,fmt='(f30.15)') density(n_point+n)
    !      end do
    !   end do
    !end do
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
    character(len=20) :: filename

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
    integer :: lun, element, nf1, nf2,len
    real(double) :: val
    character(len=15) :: filename

    ! Build a filename based on node number
    call get_file_name(stub//'matrix',numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename)
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
    character(len=15) :: filename

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
    integer :: lun, element, nf1, nf2
    character(len=15) :: filename

    ! Build a filename based on node number
    call get_file_name(stub//'grid',numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename)
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
    character(len=16) :: filename

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
    integer :: lun, block, n_point, n_i, n
    character(len=16) :: filename

    ! Build a filename based on node number
    call get_file_name('locps'//stub,numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename)
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
    character(len=10) :: filename

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
  !!
  !!  SOURCE
  !!
  subroutine grab_projs(projs,size,inode)

    use datatypes
    use global_module, only: numprocs

    ! Passed variables
    integer :: size, inode
    real(double), dimension(size) :: projs

    ! Local variables
    integer :: lun, block, n_point, n_i, n
    character(len=10) :: filename

    ! Build a filename based on node number
    call get_file_name('projs',numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename)
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
    integer :: lun, i,j,k, nf1, nf2
    character(len=18) :: filename

    ! Build a filename based on node number
    call get_file_name('blip_coeffs',numprocs,inode,filename)
    ! Open file
    call io_assign(lun)
    open(unit=lun,file=filename)
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
         '     Antonio Torralba         Veronika Brazdova       ',/,12x, &
         '          (UCL)                    (UCL)              ',/,12x, &
         '                                                      ',/,12x, &
         '      Lianheng Tong            Michiaki Arita         ',/,12x, &
         '          (UCL)                    (NIMS)             ',/,12x, &
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
    character(len=20) :: filename

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
         FLOOR (1.0 + LOG10 (REAL (numprocs, kind=double)))) - &
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
          write(io_lun,'(a,i9,a,i9,2a)') 'ID-Info: MPI = ',i,'    Process = ',pids(i),'    Host = ',hostnames(i)
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
    character(len=20),intent(in) :: filename

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
    use global_module, only: id_glob, ni_in_cell,id_glob_inv
    use GenComms, only: inode, ionode, cq_abort, gcopy

    implicit none

    ! Passed variables
    real(double), intent(out) :: velocity(1:3, ni_in_cell)
    character(len=20),intent(in) :: filename

    ! Local variables
    integer :: id_global, ni , ni2, id_tmp
    integer :: lun

    if(inode == ionode) then
       call io_assign(lun)
       open(unit=lun,file=filename)
       rewind(unit=lun)
       do id_global=1,ni_in_cell
          ni=id_glob_inv(id_global)
          if(id_glob(ni) /= id_global) &
               call cq_abort(' ERROR in global labelling ',id_global,id_glob(ni))
          read(lun,101) id_tmp, ni2, velocity(1:3, ni)
          if(ni2 /= ni) write(io_lun,*) &
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
  !!  SOURCE
  !!  
  subroutine check_stop(flag_userstop,iter)

    use GenComms, only: inode, ionode, gsum

    implicit none

    ! Passed variables
    logical :: flag_userstop
    integer, OPTIONAL :: iter

    ! Local variables
    integer :: lun, ios

    flag_userstop = .false.
    if(inode==ionode) then
       call io_assign(lun)
       open(unit=lun,file='CQ.stop',status='old',iostat=ios)
       if(ios==0) flag_userstop = .true.
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
    return
  end subroutine check_stop
  !!***

end module io_module

