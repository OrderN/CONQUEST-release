!!****h* Conquest/store_matrix
!!  NAME
!!   store_matrix
!!  PURPOSE
!!    Stores the data and information of a matrix for a set of atomic positions.
!!   to use them for a different set of atomic positions. (ex. during the MD)
!!  USES
!!
!!  AUTHOR
!!   Tsuyoshi Miyazaki
!!  CREATION DATE
!!   2016/09/29
!!  MODIFCATION
!!   - Added derived data types for XL-BOMD
!!
module store_matrix

  use datatypes
  use global_module, ONLY: numprocs
  use input_module, ONLY: io_assign, io_close
  use GenComms, ONLY: inode, ionode, cq_abort

  implicit none

! matrix_store : 
!   It includes the information of a matrix in my process.
!    (Information of Primary sets of atoms, Lists of Neighbours, ...  by global numbers.)
  type matrix_store
    character(len=80) :: name
    !integer :: inode  ! not necessary 
    integer :: n_prim
    integer :: matrix_size
   ! Size : n_prim
    integer, allocatable :: nsf_spec_i(:)
    integer, allocatable :: idglob_i(:)
    integer, allocatable :: jmax_i(:)
    integer, allocatable :: jbeta_max_i(:)
    integer, allocatable :: ibeg_data_matrix(:)  ! initial address of data_matrix for (iprim)th atom
   ! Size : (n_prim, jmax_i_max)
    integer, allocatable :: idglob_j(:,:)
    integer, allocatable :: beta_j(:,:)
   ! Size : (3, jmax_i_max, iprim)
    real(double), allocatable :: vec_Rij(:,:,:)
   ! Size : matrix_size
    real(double), allocatable :: data_matrix(:)
  end type matrix_store

! matrix_store_global : 
!   It includes the information for the whole system, including the atomic positions. 
!   Number of mpi processes may change during the multiple jobs.
!   In the future, we also assume that cell size, number of partitions may also change 
!   during a simulation.
!                  2016/10/04 Tsuyoshi Miyazaki@UCL
  type matrix_store_global
    integer :: index                     ! Will be used as "***"
    integer :: ni_in_cell, numprocs         ! number of atoms, number of MPI processes
    integer :: npcellx, npcelly, npcellz    ! partition
    real(double) :: rcellx, rcelly, rcellz  ! cell length (should be changed to 3x3 cell parameters)
    integer, allocatable :: glob_to_node(:)     ! global id of atoms -> index of MPI-process
    real(double), allocatable :: atom_coord(:,:)! atomic coordinates (3, global-id)
    real(double), allocatable :: atom_veloc(:,:)! atomic velocities  (3, global-id)
  end type matrix_store_global

  character(80),private :: RCSid = "$Id$"

contains
  ! -----------------------------------------------------------------------
  ! Subroutine dump_matrix_update
  ! -----------------------------------------------------------------------

  !!****f* store_matrix/dump_matrix_update *
  !!
  !!  NAME
  !!    dump_matrix_update
  !!  USAGE
  !!
  !!  PURPOSE
  !!    prints out the files (*matrix2.dat & InfoGlobal.dat) for updating a matrix
  !!    calls "dump_matrix2" and "InfoMatGlobal"
  !!
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Tsuyoshi Miyazaki
  !!  CREATION DATE
  !!   2016/10/05
  !!  MODIFICATION
  !!   
  !!  SOURCE
  !!

   subroutine dump_matrix_update(stub,matA,range,index,iprint_mode)
    use GenComms, ONLY: inode, ionode, cq_abort
    use global_module, ONLY: numprocs, id_glob
    use io_module, ONLY: get_file_name
    implicit none
    character(len=*),intent(in) :: stub
    integer,intent(in) :: matA
    integer,intent(in) :: range
    integer,intent(in), optional :: index
    integer,intent(inout), optional :: iprint_mode

    if(.not.present(iprint_mode)) iprint_mode = 0

     select case(iprint_mode)
      case(0)  ! Both "InfoGlobal.dat" and "*matrix2.dat" will be pinted out.
        call dump_InfoMatGlobal(index)
        call dump_matrix2(stub,matA,range)
      case(1)  ! only "*matrix2.dat" will be printed out.
        call dump_matrix2(stub,matA,range)
      case(2)  ! only "InfoGlobal.dat" will be printed out.
        call dump_InfoMatGlobal(index)
     end select ! case(iprint_mode)
 
    return
   end subroutine dump_matrix_update
  ! -----------------------------------------------------------------------
  ! Subroutine dump_matrix2
  ! -----------------------------------------------------------------------

  !!****f* store_matrix/dump_matrix2 *
  !!
  !!  NAME
  !!    dump_matrix2
  !!  USAGE
  !!
  !!  PURPOSE
  !!    Stores the information of matrix for my primary set of atoms.
  !!     Information of the atoms in my primary set (i) and their neighbour atoms (j)
  !!     
  !!
  !!    imax:           the number of primary atoms "i" in the current node "inode".
  !!    alpha_i():      the number of orbitals for each primary set of atom "i"
  !!                    in the current node "inode".
  !!    jmax():         the number of neighbour-atoms "j" for each atom "i" above.
  !!    jbeta_max():    the number of (jmax x beta) for each atom "i".
  !!    beta_j():       the number of orbitals for each neighbour-atom "j".
  !!    idglob_j():     the global labelling of neighbour-atoms "j".
  !!
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/08/21
  !!  MODIFICATION
  !!   2016/10/03 Tsuyoshi Miyazaki
  !!      1. Changed to use the derivative type (matrix_store) for tidying up MD continuation
  !!      2. Moved from io_module2
  !!  SOURCE
  !!

   subroutine dump_matrix2(stub,matA,range)
    use GenComms, ONLY: inode, ionode, cq_abort
    use global_module, ONLY: numprocs, id_glob
    use io_module, ONLY: get_file_name
    implicit none
    character(len=*),intent(in) :: stub
    integer,intent(in) :: matA
    integer,intent(in) :: range
    type(matrix_store):: tmp_matrix_store

    integer :: lun, iprim, nprim, jmax, jj, ibeg, jbeta_alpha, len
    character(20) :: file_name

    ! set_matrix_store : build tmp_matrix_store
    call set_matrix_store(stub,matA,range,tmp_matrix_store)

    ! Actual Dump (from dump_matrix2)
     ! First, make a file based upon the node ID.
     call get_file_name(stub//'matrix2',numprocs,inode,file_name)
     call io_assign(lun)
    open (lun,file=file_name)

      ! 1. node ID, no. of PS of atoms "i".
        nprim=tmp_matrix_store%n_prim
        write (lun,*) inode, nprim
      ! 2. no. of alpha for each "i".
        write (lun,*) tmp_matrix_store%nsf_spec_i(1:nprim)
      ! NOTE: The followings are written out with neighbour-labelling
      ! 3. no. of the neighbours "j" for each "i".
      ! 4. no. of "neighbour-j x beta" for each "i".
        write (lun,*) tmp_matrix_store%jmax_i(1:nprim)
        write (lun,*) tmp_matrix_store%jbeta_max_i(1:nprim)

     !I will change the order of dumping in the following, later.   2016/09/30: TM@UCL
      if(nprim .GT. 0) then
       do iprim=1,nprim
          jmax = tmp_matrix_store%jmax_i(iprim)
          write (lun,*) tmp_matrix_store%idglob_i(iprim)
          write (lun,*) tmp_matrix_store%beta_j(1:jmax,iprim)
          write (lun,*) tmp_matrix_store%idglob_j(1:jmax,iprim)
         do jj=1,jmax
          write (lun,*) tmp_matrix_store%vec_Rij(1:3,jj,iprim)
         enddo !jj=1,jmax
         if(iprim < nprim) then
           len = tmp_matrix_store%ibeg_data_matrix(iprim+1)-tmp_matrix_store%ibeg_data_matrix(iprim)
         else
           len = tmp_matrix_store%matrix_size-tmp_matrix_store%ibeg_data_matrix(iprim)+1
         endif
           ibeg = tmp_matrix_store%ibeg_data_matrix(iprim)
         do jbeta_alpha = 1, len
          write (lun,fmt='(f25.18)') tmp_matrix_store%data_matrix(ibeg+jbeta_alpha-1)
         enddo !jbeta_alpha = 1, len
       enddo !iprim=1,nprim
      endif  ! (nprim .GT. 0) 

    ! Close the file in the end.
    call io_close(lun)
    
    ! free_matrix_store : free tmp_matrix_store
    call free_matrix_store(tmp_matrix_store)

    return
   end subroutine dump_matrix2

  ! -----------------------------------------------------------------------
  ! Subroutine set_matrix_store  & free_matrix_store
  ! -----------------------------------------------------------------------
  !!****f* store_matrix/set_matrix_store *
  !!
  !!  NAME
  !!    set_matrix_store
  !!  USAGE
  !!
  !!  PURPOSE
  !!    Stores the information as for the primary atoms and its
  !!    neighbours first, then the matrix elements are followed.
  !!    Made from dump_matrix2 written by M. Arita
  !!
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Tsuyoshi Miyazaki/Michiaki Arita
  !!  CREATION DATE
  !!   2016/10/03
  !!  MODIFICATION
  !!
  !!  SOURCE
  !!

   subroutine set_matrix_store(stub,matA,range,matinfo)
    use numbers
    use global_module, ONLY: id_glob, species_glob, io_lun
    use species_module, ONLY: nsf_species
    use group_module, ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module, ONLY: BCS_parts
    use matrix_data, ONLY: mat
    use mult_module, ONLY: mat_p

    implicit none
    character(len=*),intent(in) :: stub
    integer,intent(in) :: matA
    integer,intent(in) :: range
    type(matrix_store), intent(out) :: matinfo
    !Local
    integer :: nprim, istat
    integer :: np, ni, iprim, atom_num, jmax_i_max
    integer :: ist, n_beta, neigh, gcspart, ibeg, j_global_part, len

    !name & n_prim
      matinfo%name=stub
      matinfo%n_prim=bundle%n_prim
    !allocation 1 
      nprim=matinfo%n_prim   ! local variable
      allocate(matinfo%nsf_spec_i(nprim), matinfo%idglob_i(nprim), matinfo%jmax_i(nprim), &
               matinfo%jbeta_max_i(nprim), matinfo%ibeg_data_matrix(nprim), STAT=istat)
      if(istat .NE. 0) call cq_abort('Allocation 1 in set_matrix_store',istat,nprim)
    
      matinfo%nsf_spec_i(1:nprim) =nsf_species(bundle%species(1:nprim))

   !!set up : jmax_i, jbeta_max_i

    iprim = 0
    jmax_i_max = -1
    matinfo%jbeta_max_i(:) = 0

    do np = 1, bundle%groups_on_node
      if (bundle%nm_nodgroup(np).GT.0) then
        do ni = 1, bundle%nm_nodgroup(np)
           iprim = iprim + 1
           atom_num = bundle%nm_nodbeg(np)+ni-1
           if(iprim .NE. atom_num) call cq_abort('iprim .ne. atom_num ?',iprim,atom_num)
           matinfo%idglob_i(iprim) = bundle%ig_prim(iprim)
           matinfo%jmax_i(iprim) =  mat(np,range)%n_nab(ni)
           if (matinfo%jmax_i(iprim).GT.jmax_i_max) jmax_i_max = matinfo%jmax_i(iprim)
          do neigh = 1, matinfo%jmax_i(iprim)
            ist = mat(np,range)%i_acc(ni) + neigh - 1 ! Neighbour-labeling
            n_beta = mat(np,range)%ndimj(ist)
            matinfo%jbeta_max_i(iprim) = matinfo%jbeta_max_i(iprim) + n_beta
          enddo !(neigh)
        enddo !(ni)
      endif
    enddo !(np)

    !allocation 2 
      allocate(matinfo%idglob_j(jmax_i_max,nprim), matinfo%beta_j(jmax_i_max,nprim), STAT=istat)
      if(istat .NE. 0) call cq_abort('Allocation 2 in set_matrix_store', nprim,jmax_i_max)
    !allocation 3 
      allocate(matinfo%vec_Rij(3,jmax_i_max,nprim), STAT=istat)
      if(istat .NE. 0) call cq_abort('Allocation 3 in set_matrix_store', nprim,jmax_i_max)
     
   !!set up : idglob_j, beta_j, vec_Rij
    matinfo%idglob_j(:,:) = 0
    matinfo%beta_j(:,:) = 0
    iprim = 0
    ibeg = 1
    matinfo%matrix_size = 0
    do np = 1, bundle%groups_on_node
      if (bundle%nm_nodgroup(np).GT.0) then
        do ni = 1, bundle%nm_nodgroup(np)
          atom_num=bundle%nm_nodbeg(np)+ni-1
          iprim = iprim + 1
          if(atom_num .NE. iprim) write(io_lun,*) 'Warning!! in set_matrix_store : iprim & atom_num = ',iprim,atom_num

          do neigh = 1, mat(np,range)%n_nab(ni)
            ist = mat(np,range)%i_acc(ni) + neigh - 1  ! Neighbour-labeling
            gcspart=BCS_parts%icover_ibeg(mat(np,range)%i_part(ist))+mat(np,range)%i_seq(ist)-1
            j_global_part = BCS_parts%lab_cell(mat(np,range)%i_part(ist))

            matinfo%idglob_j(neigh,iprim) = id_glob(parts%icell_beg(j_global_part) &
                                                +mat(np,range)%i_seq(ist)-1)
            matinfo%beta_j(neigh,iprim) = nsf_species(species_glob(matinfo%idglob_j(neigh,iprim)))
            matinfo%vec_Rij(1,neigh,iprim) = bundle%xprim(iprim) - BCS_parts%xcover(gcspart)
            matinfo%vec_Rij(2,neigh,iprim) = bundle%yprim(iprim) - BCS_parts%ycover(gcspart)
            matinfo%vec_Rij(3,neigh,iprim) = bundle%zprim(iprim) - BCS_parts%zcover(gcspart)

          enddo !neigh = 1, mat(np,range)%n_nab(ni)

          len = matinfo%jbeta_max_i(iprim)*nsf_species(bundle%species(iprim))
          matinfo%matrix_size = matinfo%matrix_size + len
          matinfo%ibeg_data_matrix(iprim) = ibeg
          ibeg = ibeg + len
        enddo !(ni)
      endif !(bundle%nm_nodgroup(np).GT.0) then
    enddo !(np)

    allocate(matinfo%data_matrix(matinfo%matrix_size), STAT = istat)
     if(istat .NE. 0) call cq_abort('Allocation 4 in set_matrix_store', istat, matinfo%matrix_size)
    matinfo%data_matrix(1:matinfo%matrix_size)=mat_p(matA)%matrix(1:matinfo%matrix_size)

    return
   end subroutine set_matrix_store

   subroutine free_matrix_store(matinfo)
    use GenComms, ONLY: cq_abort
    implicit none
    type(matrix_store), intent(inout) :: matinfo
    integer :: istat
   
     deallocate(matinfo%data_matrix, matinfo%vec_Rij, matinfo%beta_j, matinfo%idglob_j, &
                matinfo%ibeg_data_matrix, matinfo%jbeta_max_i, matinfo%jmax_i, matinfo%idglob_i, &
                matinfo%nsf_spec_i, STAT=istat)
     if(istat .NE. 0) call cq_abort("Error in deallocation at free_matrix_store",istat)
    return
   end subroutine free_matrix_store

  ! -----------------------------------------------------------------------
  ! Subroutine dump_InfoMatGlobal
  ! -----------------------------------------------------------------------
  !!****f* store_matrix/dump_InfoMatGlobal *
  !!
  !!  NAME
  !!   dump_InfoMatGlobal
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Writes out the global information
  !!    l.1 no. of atoms, no. of procs.
  !!    l.2 no. of partitions along a,b and c axes
  !!    l.3 cell size
  !!    l.4 relations between atoms (global label) and procs.
  !!    (optional) MD iteration
  !!    l.5 atomic position
  !!    l.6 atomic velocity
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Tsuyoshi Miyazaki/Michiaki Arita
  !!  CREATION DATE
  !!   2016/10/04 (2013/08/21)
  !!  MODIFICATION
  !!
  !!  SOURCE
  !!
  subroutine dump_InfoMatGlobal(MDstep,velocity)

    ! Module usage
    use global_module, ONLY: ni_in_cell,numprocs,rcellx,rcelly,rcellz,id_glob, io_lun
    use GenComms, ONLY: cq_abort, inode, ionode, my_barrier
    use group_module, ONLY: parts
    use io_module, ONLY: get_file_name, get_file_name_2rank

    implicit none
    type(matrix_store_global):: mat_global_tmp
    integer, intent(in), optional :: MDstep
    real(double), intent(in), optional :: velocity(1:3,ni_in_cell)
    integer :: lun, istat, iglob
    integer :: step_local
    character(len=80) :: filename

    if(present(MDstep)) then
     step_local = MDstep
    else
     step_local = 0
    endif

    if(present(velocity)) then
     call set_InfoMatGlobal(mat_global_tmp, step=step_local, velocity_in=velocity)
    else
     call set_InfoMatGlobal(mat_global_tmp, step=step_local)
    endif

    ! Open InfoGlobal.dat and write data.
    if(inode == ionode) then
     call io_assign(lun)
     call get_file_name_2rank('InfoGlobal',filename,step_local)
     write(io_lun,*) ' filename = ',filename
     open (lun,file=filename,iostat=istat)
     !open (lun,file='InfoGlobal.dat',iostat=istat)
     rewind lun
     if (istat.GT.0) call cq_abort('Fail in opening InfoGlobal.dat .')
     write (lun,*) mat_global_tmp%ni_in_cell, mat_global_tmp%numprocs
     write (lun,*) mat_global_tmp%npcellx,mat_global_tmp%npcelly,mat_global_tmp%npcellz
     write (lun,*) mat_global_tmp%rcellx,mat_global_tmp%rcelly,mat_global_tmp%rcellz
     write (lun,*) mat_global_tmp%glob_to_node(1:mat_global_tmp%ni_in_cell)
     !if (present(MDstep)) write (lun,*) mat_global_tmp%index
     write (lun,*) mat_global_tmp%index
 
     do iglob=1, mat_global_tmp%ni_in_cell
      write (lun,101) iglob, mat_global_tmp%atom_coord(1:3, iglob)
      101 format(5x,i8,3x,3e20.10)
     enddo !iglob=1, mat_global_tmp%ni_in_cell
     call io_close (lun)
    endif

    call free_InfoMatGlobal(mat_global_tmp)
    return
  end subroutine dump_InfoMatGlobal
  !!***

  ! -----------------------------------------------------------------------
  ! Subroutine set_InfoMatGlobal & free_InfoMatGlobal
  ! -----------------------------------------------------------------------
  !!****f* store_matrix/set_InfoMatGlobal *

  subroutine set_InfoMatGlobal(mat_glob, step, velocity_in)
    
    ! Module usage
    use numbers, ONLY: zero
    use global_module, ONLY: ni_in_cell,numprocs,rcellx,rcelly,rcellz,id_glob, atom_coord
    use GenComms, ONLY: cq_abort
    use group_module, ONLY: parts

    implicit none
    integer, intent(in), optional :: step
    real(double), intent(in), optional :: velocity_in(3,ni_in_cell)
    type(matrix_store_global),intent(out) :: mat_glob
    integer :: istat, ind_part, id_node, ni, id_global
   
    if(present(step))then
       mat_glob%index      = step
    else
       mat_glob%index      = 0
    endif

    mat_glob%ni_in_cell = ni_in_cell
    mat_glob%numprocs   = numprocs

    mat_glob%rcellx     = rcellx
    mat_glob%rcelly     = rcelly
    mat_glob%rcellz     = rcellz

    mat_glob%npcellx     = parts%ngcellx
    mat_glob%npcelly     = parts%ngcelly
    mat_glob%npcellz     = parts%ngcellz

    !allocation of glob_to_node, atom_coord, atom_veloc
    allocate(mat_glob%glob_to_node(ni_in_cell), mat_glob%atom_coord(3,ni_in_cell), &
             mat_glob%atom_veloc(3,ni_in_cell), STAT=istat)
    if(istat .NE. 0) call cq_abort('Error : allocation in set_InfoMatGlobal',istat,ni_in_cell)

    do ind_part = 1, parts%mx_gcell
      id_node = parts%i_cc2node(ind_part) ! CC labeling
      if (parts%nm_group(ind_part).NE.0) then
        do ni = 1, parts%nm_group(ind_part)
          id_global = id_glob(parts%icell_beg(ind_part)+ni-1)
          mat_glob%glob_to_node(id_global) = id_node
          if (id_node.EQ.0 .OR. id_node .GT. numprocs) &
             call cq_abort('Error: glob_to_node in set_InfoMatGlobal',id_node, numprocs)
        enddo
      endif
    enddo !(ind_part, mx_gcell)
      
    mat_glob%atom_coord(1:3,1:ni_in_cell) = atom_coord(1:3,1:ni_in_cell) 

   !velocity ( global_id, partition labelling)
    if(present(velocity_in)) then
     do ni=1, ni_in_cell
       id_global= id_glob(ni)
       mat_glob%atom_veloc(1:3,id_global)=velocity_in(1:3,ni)
     enddo
    else
       mat_glob%atom_veloc(1:3,:)=zero
    endif

   return
  end subroutine set_InfoMatGlobal

  !!***
  subroutine free_InfoMatGlobal(mat_glob)
    use GenComms, ONLY: cq_abort
    implicit none
    integer :: istat
    !type(matrix_store_global),intent(out) :: mat_glob
    type(matrix_store_global) :: mat_glob

    deallocate(mat_glob%atom_veloc, mat_glob%atom_coord, mat_glob%glob_to_node, STAT=istat)
    if(istat .NE. 0) call cq_abort('Error : deallocation in free_InfoMatGlobal',istat)
 
   return
  end subroutine free_InfoMatGlobal

  !!***
  ! -----------------------------------------------------------------------
  ! Subroutine grab_InfoMatGlobal
  ! -----------------------------------------------------------------------
  !!****f* store_matrix_module/grab_InfoMatGlobal *
  !!
  !!  NAME
  !!   grab_InfoMatGlobal
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Reads the global data used for communication & reconstruction 
  !!   Now, the 
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Tsuyoshi Miyazaki/Michiaki Arita
  !!  CREATION DATE
  !!   2017/10/13 made from grab_InfoGlobal
  !!  MODIFICATION
  !!
  !!  SOURCE
  !!
  subroutine grab_InfoMatGlobal(InfoGlob,index)

    ! Module usage
    use GenComms, ONLY: inode,ionode,gcopy
    use io_module, ONLY: get_file_name_2rank
    use global_module, ONLY: io_lun, ni_in_cell

    ! passed variables
    type(matrix_store_global),intent(out) :: InfoGlob
    integer,intent(in), optional :: index

    ! local variables
    integer :: lun, istat, iglob, ig
    integer :: index_local
    character(len=80) :: filename

    ! check whether members of InfoGlob has been allocated or not.
     if(.not.allocated(InfoGlob%glob_to_node)) then
       !allocation of glob_to_node, atom_coord, atom_veloc
        write(io_lun,*) "allocation in grab_InfoMatGlobal"
        allocate(InfoGlob%glob_to_node(ni_in_cell), InfoGlob%atom_coord(3,ni_in_cell), &
                 InfoGlob%atom_veloc(3,ni_in_cell), STAT=istat)
        if(istat .NE. 0) call cq_abort('Error : allocation in grab_InfoMatGlobal',istat,ni_in_cell)
     endif

    index_local=0
    if(present(index)) index_local=index

   !Reading  "InfoGlobal.ind***"
    if(inode == ionode) then
     call io_assign(lun)
     call get_file_name_2rank('InfoGlobal',filename,index_local)
     open (lun,file=filename,status='old',iostat=istat)
     if (istat.GT.0) then
      write(io_lun,*) " grab_InfoMatGlobal: Error in opening InfoGlobal",filename
      call cq_abort('Fail in opening InfoGlobal.dat')
     endif
 
      read(lun,*) InfoGlob%ni_in_cell, InfoGlob%numprocs
       if(ni_in_cell .NE. InfoGlob%ni_in_cell) &
        call cq_abort('Error in grab_InfoMatGlobal: ni_in_cell= ',ni_in_cell,InfoGlob%ni_in_cell)
      read(lun,*) InfoGlob%npcellx,InfoGlob%npcelly,InfoGlob%npcellz
      read(lun,*) InfoGlob%rcellx,InfoGlob%rcelly,InfoGlob%rcellz
      read(lun,*) InfoGlob%glob_to_node(1:InfoGlob%ni_in_cell)
      read(lun,*) InfoGlob%index
 
     do ig=1, InfoGlob%ni_in_cell
      read(lun,101) iglob, InfoGlob%atom_coord(1:3, ig)
      101 format(5x,i8,3x,3e20.10)
     enddo !ig=1, InfoGlob%ni_in_cell
  
     ! Velociy will be needed in the future ?
     !do ig=1, InfoGlob%ni_in_cell
     ! read(lun,101) iglob, InfoGlob%atom_veloc(1:3, ig)
     !enddo !ig=1, InfoGlob%ni_in_cell

     call io_close (lun)
    endif !(inode == ionode) 

   !gcopy InfoGlob to all nodes
    call gcopy(InfoGlob%ni_in_cell)
    call gcopy(InfoGlob%numprocs)

    call gcopy(InfoGlob%npcellx)
    call gcopy(InfoGlob%npcelly)
    call gcopy(InfoGlob%npcellz)

    call gcopy(InfoGlob%rcellx)
    call gcopy(InfoGlob%rcelly)
    call gcopy(InfoGlob%rcellz)

    call gcopy(InfoGlob%glob_to_node, ni_in_cell)
    call gcopy(InfoGlob%index)
    call gcopy(InfoGlob%atom_coord, 3, ni_in_cell)
    !call gcopy(InfoGlob%atom_veloc, 3, ni_in_cell)

   !gcopy InfoGlob : END

    return
  end subroutine grab_InfoMatGlobal
  !!***


  ! -----------------------------------------------------------------------
  ! Subroutine make_index_iter
  ! -----------------------------------------------------------------------
  !!****f* store_matrix/make_index_iter
  !!
  !!  NAME
  !!   make_index_iter
  !!  USAGE
  !!   make_index_iter(mx_store, iter_present, index_iter)
  !!  PURPOSE
  !!   to get the iteration numbers for the latest mx_store steps
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Tsuyoshi Miyazaki
  !!  CREATION DATE
  !!   2017 June
  !!  MODIFICATION
  !!
  !!  SOURCE
  !!

  subroutine make_index_iter(mx_store, iter_present, index_iter)
   implicit none
   integer, intent(in) :: mx_store
   integer, intent(in) :: iter_present
   integer, intent(out) :: index_iter(0:mx_store-1)
   integer :: iter

   index_iter(:) = 0
   if(iter_present < mx_store) then
    do iter = 0, iter_present-1
      index_iter(iter)=iter_present-iter
    enddo
   else
    do iter = 0, mx_store-1
      index_iter(iter)=mod(iter_present-iter-1,mx_store)+1
    enddo
   endif
  return
  end subroutine make_index_iter

end module store_matrix
