! ----------------------------------------------------------------------------------
! Module atom_dispenser_module
! ----------------------------------------------------------------------------------

!!****h* Conquest/atom_dispenser_module *
!!  NAME
!!    atom_dispenser_module
!!  PURPOSE
!!    Finds out the relations between atoms and partitions
!!  AUTHOR
!!    Michiaki Arita
!!  CREATION DATE
!!    2013/07/01
!!  MODIFICATION HISTORY
!!
!!***
module atom_dispenser

  use datatypes
  use group_module, ONLY: parts
  use global_module, ONLY: flag_fractional_atomic_coords,rcellx,rcelly,rcellz, &
                           flag_MDdebug,shift_in_bohr,Iprint_MD
  use GenComms, ONLY: cq_abort

  implicit none

  ! RCS tag for object file idetification
  character(80), private :: RCSid = "$Id: atom_dispenser_module.f90 "

  logical, parameter :: flag_ortho = .true.

!!***

contains


  ! ---------------------------------------------------------------------------------
  ! Subroutine atom2part
  ! ---------------------------------------------------------------------------------

  !!****f* atom_dispenser/atom2part
  !!
  !!  NAME
  !!    atom2part
  !!  USAGE
  !!
  !!  PURPOSE
  !!    Eliminate an ambiguity in the case where atoms are on a partition
  !!    boundry. Also used in finding the partition which j-atom belongs to
  !!    at the new MD/cg step.
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !! 
  !!  CREATION DATE
  !!
  !!  MODIFICATION
  !!
  !!  SOURCE
  !!
  subroutine atom2part(x,y,z,ind_part,px,py,pz,atom_id)

    ! Module usage
    use numbers
    use input_module, ONLY: io_assign, io_close
    use global_module, ONLY: io_lun
    use GenComms, ONLY: inode, ionode
    ! DB
    use io_module, ONLY: get_file_name
    use global_module, ONLY: numprocs
    use cover_module, ONLY: BCS_parts
    use dimens, ONLY: r_super_x,r_super_y,r_super_z

    implicit none

    ! passed variables
    real(double) :: x, y, z
    integer :: px, py, pz
    integer :: ind_part, atom_id   !! NOTE: These are only for debugging. !!

    ! local variables
    integer :: i
    real(double) :: plen_x, plen_y, plen_z, px_max, px_min, py_max, py_min, &
                    pz_max, pz_min
    real(double) :: x_eps, y_eps, z_eps, eps
    logical :: flag_px, flag_py, flag_pz
 
    ! Debug
    integer :: lun, stat, ind_part2
    character(13) :: file_name
    !integer :: n_part
    !integer, allocatable :: n_atom_part(:)


    ! Set IO. print level for the following.
    !if (inode.EQ.ionode) write (io_lun, '(a)') "Enter atom2part."


    !! ----------- DEBUG ----------- !!
    ! NOTE: This file will be sizable.
    if (flag_MDdebug) then
      call get_file_name('dispenser',numprocs,inode,file_name)
      call io_assign(lun)
      open (lun,file=file_name,position='append')
    endif
    !! ----------- DEBUG ----------- !!
      

    ! Firstly, we need to shift the atoms by eps before deciding the 
    ! partition and updating parts.
    if (flag_fractional_atomic_coords) then
      !WRONG x = x * rcellx
      !WRONG y = y * rcelly
      !WRNOG z = z * rcellz
      eps = 1.0E-08     ! May be changed later.
    else
      !eps = 1.0E-03    ! May be changed later.
      eps = 1.0E-04     ! May be changed later.
    endif

    !TM  Eps should be considered with Cartesian (bohr) units
    eps = 1.0E-03_double    
    x_eps = x + eps
    y_eps = y + eps
    z_eps = z + eps

    if (flag_ortho) then ! if the system is orthorhombic.
      ! Calculate the partition lengths.
      plen_x = rcellx/real(parts%ngcellx,double)
      plen_y = rcelly/real(parts%ngcelly,double)
      plen_z = rcellz/real(parts%ngcellz,double)

      !x_eps = x_eps - floor(x_eps)
      !y_eps = y_eps - floor(y_eps)
      !z_eps = z_eps - floor(z_eps)

         px = floor(x_eps/plen_x) + 1
         py = floor(y_eps/plen_y) + 1
         pz = floor(z_eps/plen_z) + 1

    else ! if NOT orthorhombic.
      call cq_abort('ERROR in atom2part: flag_ortho ')
      ! You need to consider alpha, beta and gamma, but leave it for now.
      plen_x = rcellx/real(parts%ngcellx,double)
      plen_y = rcelly/real(parts%ngcelly,double)
      plen_z = rcellz/real(parts%ngcellz,double)
    endif

    ind_part  = (px-1) * (parts%ngcelly*parts%ngcellz) + (py-1) * parts%ngcellz + pz !! --- DEBUG -- !! 

    !! ----------- DEBUG ----------- !!
    if (flag_MDdebug) then
      write (lun, '(a,1x,i8)') "## Globel atom ID: ", atom_id
      write (lun, '(a,1x,l5,1x,f15.10)') "Fractional? / eps: ", flag_fractional_atomic_coords, eps
      write (lun, '(a,1x,3f15.10)') "cell lengths: ", rcellx, rcelly, rcellz
      !write (lun, '(a,1x,3i8)') "# of parts: ", parts%ngcellx, parts%ngcelly, parts%ngcellz
      write (lun, '(a,1x,3f15.10)') "Partition lengths: ", plen_x, plen_y, plen_z
      write (lun, '(a,1x,6f15.10)') "x,y,z; x,y,z_eps: ", x, y, z, x_eps, y_eps, z_eps
      !write (lun, *) "position (x, y, z): ", x, y, z
      write (lun, '(a,1x,5i7)') "px,py,pz, ind_part: ", px, py, pz, ind_part
      write (lun, *) ""
      call io_close(lun)
    endif
    !! ----------- DEBUG ----------- !!

    return
  end subroutine atom2part
  !!***


  ! ---------------------------------------------------------------------------------
  ! Subroutine allatom2part
  ! ---------------------------------------------------------------------------------

  !!****f* atom_dispenser/allatom2part
  !!
  !!  NAME
  !!    allatom2part
  !!  USAGE
  !!
  !!  PURPOSE
  !!    Eliminates an ambiguity in the case where atoms are on a partition
  !!    boundry. Also used in finding the partition which j-atom belongs to
  !!    at the new MD/cg step. Unlike atom2part, this subroutine is called 
  !!    when updating the member information.
  !!  INPUTS
  !!    ind_part
  !!  USES
  !!    global_module,dimens
  !!  AUTHOR
  !!    Michiaki Arita 
  !!  CREATION DATE
  !!    2013/07/01
  !!  MODIFICATION
  !!
  !!  SOURCE
  !!
  subroutine allatom2part(ind_part)

    ! Module usage
    use numbers
    use GenComms, ONLY: gcopy,inode,ionode
    use global_module, ONLY: ni_in_cell,x_atom_cell,y_atom_cell,z_atom_cell, &
                             atom_coord
    use dimens, ONLY: r_super_x,r_super_y,r_super_z
    ! DB
    use input_module, ONLY: io_assign,io_close
    use global_module, ONLY: io_lun,id_glob,numprocs

    implicit none

    ! passed variables
    integer :: ind_part(:)

    ! local variables
    integer :: n_atom,px,py,pz,i
    integer :: stat_alloc
    real(double) :: px_min,px_max,py_min,py_max,pz_min,pz_max, &
                    plen_x,plen_y,plen_z,cellx,celly,cellz
    real(double) :: x_eps,y_eps,z_eps,eps,x,y,z 
    logical :: flag_px,flag_py,flag_pz

    ! DB
    integer :: lun,stat
    character(20) :: file_name

    if (inode.EQ.ionode .AND. Iprint_MD.GT.2) &
      write (io_lun,*) "Entering allatom2part."

    !! ---------------- DEBUG ------------------ !!
    if (flag_MDdebug .AND. inode.EQ.ionode) then
      call io_assign(lun)
      open (lun,file='allatom2part.dat',action='write',iostat=stat)
      if (stat.NE.0) call cq_abort('Error: Fail in opening allatom2part.dat .')
      write (lun,*) "No. of atoms in system:", ni_in_cell
      write (lun,*) "atom_coord(1,:)", atom_coord(1,1:)
      write (lun,*) "atom_coord(2,:)", atom_coord(2,1:)
      write (lun,*) "atom_coord(3,:)", atom_coord(3,1:)
      write (lun,*) "x_atom_cell:", x_atom_cell(1:)
      write (lun,*) "y_atom_cell:", y_atom_cell(1:)
      write (lun,*) "z_atom_cell:", z_atom_cell(1:)
      write (lun,*) ""
    endif
    !! ---------------- DEBUG ------------------ !!

    cellx = r_super_x
    celly = r_super_y
    cellz = r_super_z
    ! Firstly, we need to wrap coordinates sutisfying p.b.c.
    ! Then, we need to shift the atoms by eps before deciding the
    ! partition and updating parts.
    ! TM eps should be considered for Cartesian (bohr) units
    eps = shift_in_bohr
    if (flag_ortho) then
      ! Calculate the partition lengths
      plen_x = rcellx/real(parts%ngcellx,double)
      plen_y = rcelly/real(parts%ngcelly,double)
      plen_z = rcellz/real(parts%ngcellz,double)
      if (flag_MDdebug .AND. inode.EQ.ionode) &
        write (lun,*) "plen_x,y,z:", plen_x,plen_y,plen_z !DB
      ! Get the partition (sim-cell gp (CC)) to which each atom belongs.
      flag_px=.false. ; flag_py=.false. ; flag_pz=.false.
      do n_atom = 1, ni_in_cell
        ! Wrap coordinates.
        ! NOTE: We cannot use xyz_atom_cell since they depend on parts labelling.
        atom_coord(1,n_atom) = atom_coord(1,n_atom) - floor((atom_coord(1,n_atom)+eps)/cellx)*cellx
        atom_coord(2,n_atom) = atom_coord(2,n_atom) - floor((atom_coord(2,n_atom)+eps)/celly)*celly
        atom_coord(3,n_atom) = atom_coord(3,n_atom) - floor((atom_coord(3,n_atom)+eps)/cellz)*cellz
        x_eps = atom_coord(1,n_atom)+eps
        y_eps = atom_coord(2,n_atom)+eps
        z_eps = atom_coord(3,n_atom)+eps

        !! -------------- DEBUG --------------- !!
        if (flag_MDdebug .AND. inode.EQ.ionode) then
          write (lun,*) ""
          write (lun,*) "glob ID          :", n_atom
          write (lun,*) "x_eps,y_eps,z_eps:"
          write (lun,*)  x_eps,y_eps,z_eps
          write (lun,*) "atom_coord(1:3, n_atom):"
          write (lun,*) atom_coord(1,n_atom),atom_coord(2,n_atom),atom_coord(3,n_atom)
        endif
        !! -------------- DEBUG --------------- !!

        px = floor(x_eps/plen_x) + 1
        py = floor(y_eps/plen_y) + 1
        pz = floor(z_eps/plen_z) + 1
        flag_px = .true. ; flag_py = .true. ; flag_pz = .true.
        if(px <= 0 .or. px > parts%ngcellx) then
          flag_px = .false.
          write(io_lun,*) ' ERROR : flag_px , px = ',px,' x_eps, plen_x, cellx = ', x_eps, plen_x, cellx
        endif
        if(py <= 0 .or. py > parts%ngcelly) then
          flag_py = .false.
          write(io_lun,*) ' ERROR : flag_py , py = ',py,' y_eps, plen_y, celly = ', y_eps, plen_y, celly
        endif
        if(pz <= 0 .or. pz > parts%ngcellz) then
          flag_pz = .false.
          write(io_lun,*) ' ERROR : flag_pz , pz = ',pz,' z_eps, plen_z, cellz = ', z_eps, plen_z, cellz
        endif
        if(.not.flag_px .or. .not.flag_py .or. .not.flag_pz) then
          write(io_lun,*) ' flag_px, flag_py, flag_pz = ',flag_px, flag_py, flag_pz
          call cq_abort(' ERROR : flag_pxyz ',n_atom)
        endif
        ! Get the partition in sim-cell (CC).
        ind_part(n_atom) = (px-1)*(parts%ngcelly*parts%ngcellz) + (py-1)*parts%ngcellz + pz
        if (flag_MDdebug .AND. inode.EQ.ionode) &
          write (lun,*) "px,py,pz, ind_part:", px,py,pz,ind_part(n_atom) !DB
      enddo  !(n_atom, ni_in_cell)	
    else
      call cq_abort('Error: flag_ortho ')
    endif  !(flag_ortho)
    if ((.NOT. flag_px) .OR. (.NOT. flag_py) .OR. (.NOT. flag_pz)) &
      call cq_abort('Error: Fail in finding partitions: allatom2part')

    !! ---------------- DEBUG ------------------ !!
    if (flag_MDdebug .AND. inode.EQ.ionode) then
      write (lun,*) ""
      write (lun,*) "ind_part(n_atom):", ind_part(1:)
      call io_close (lun)
    endif
    !! ---------------- DEBUG ------------------ !!

    return
  end subroutine allatom2part
  !!***

end module atom_dispenser
