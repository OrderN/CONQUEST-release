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
  use GenComms, ONLY: inode, ionode, cq_abort, myid
  use io_module, ONLY: flag_MatrixFile_RankFromZero, flag_MatrixFile_BinaryFormat, &
                       flag_MatrixFile_BinaryFormat_Grab, flag_MatrixFile_BinaryFormat_Dump

  implicit none

! matrix_store : 
!   It includes the information of a matrix in my process.
!    (Information of Primary sets of atoms, Lists of Neighbours, ...  by global numbers.)
!  
!   To make the arrays defined with type used as those with "InfoMatrixFile"
!   ranks of beta_j(2rank) & vec_Rij (3rank) must be reduced by introducing ibeg_Pij(:)
!
  type matrix_store
    character(len=80) :: name
    !integer :: inode  ! not necessary 
    integer :: n_prim
    integer :: nspin
    integer :: matrix_size
   ! Size : n_prim
    integer, allocatable :: nsf_spec_i(:)
    integer, allocatable :: idglob_i(:)
    integer, allocatable :: jmax_i(:)
    integer, allocatable :: jbeta_max_i(:)
    integer, allocatable :: ibeg_data_matrix(:)  ! initial address of data_matrix for (iprim)th atom
    integer, allocatable :: ibeg_Rij(:)          ! initial address of ibeg_Pij
   ! Size : (n_prim)
    integer, allocatable :: idglob_j(:)
    integer, allocatable :: beta_j(:)
   ! Size : (3, # of j)
    real(double), allocatable :: vec_Rij(:,:)
   ! Size : matrix_size
    real(double), allocatable :: data_matrix(:,:)
  end type matrix_store

! matrix_store_global : 
!   It includes the information for the whole system, including the atomic positions. 
!   Number of mpi processes may change during the multiple jobs.
!   In the future, we also assume that cell size, number of partitions may also change 
!   during a simulation.
!                  2016/10/04 Tsuyoshi Miyazaki@UCL
  type matrix_store_global
    integer :: MDstep                       ! Will be used as "***"
    integer :: ni_in_cell, numprocs         ! number of atoms, number of MPI processes
    integer :: npcellx, npcelly, npcellz    ! partition
    real(double) :: rcellx, rcelly, rcellz  ! cell length (should be changed to 3x3 cell parameters)
    integer, allocatable :: glob_to_node(:)     ! global id of atoms -> index of MPI-process
    real(double), allocatable :: atom_coord(:,:)! atomic coordinates (3, global-id)
    real(double), allocatable :: atom_veloc(:,:)! atomic velocities  (3, global-id)
  end type matrix_store_global

! InfoMatrixFile : 
!   Moved from io_module2
!  The information in this type is basically the same as in "matrix_store".
!  But "InfoMatrixFile" is for storing the information from the files (IOs) at previous steps
!  or previous jobs. Since the number of (MPI)processes may change during the contiguous jobs,
!  each present (MPI) process may need to read multiple files. Thus, the arrays with this type
!  are defined as 1-rank array.  
!  
!  At present, in io_module2, they are defined by 
!  type(InfoMatrixFile), pointer :: InfoL(:),InfoT(:),InfoK(:),    &
!                                       InfoX(:),InfoXvel(:),InfoS(:), & ! for XL-BOMD
!                                       InfoX1(:),InfoX2(:),InfoX3(:), & ! for dissipation
!                                       InfoX4(:),InfoX5(:),InfoX6(:), &
!                                       InfoX7(:),InfoX8(:),InfoX9(:), &
!                                       InfoX10(:)
!
!                  2017/11/10 Tsuyoshi Miyazaki 
!!!!
!  N_MATRIX should be defined for each i_atom ?????
!!!!
  type InfoMatrixFile
    integer :: index_node
    integer :: natom_i
    integer :: nspin
    ! Size : natom_i
    integer, pointer :: alpha_i(:)
    integer, pointer :: idglob_i(:)
    integer, pointer :: jmax_i(:)
    integer, pointer :: jbeta_max_i(:)
    integer, pointer :: ibeg_dataL(:)
    integer, pointer :: ibeg_Pij(:)
    ! Size : jmax_i
    integer, pointer :: idglob_j(:)
    integer, pointer :: beta_j_i(:)
    ! Size : 3 x jmax_i
    real(double), pointer :: rvec_Pij(:,:)
    ! Size :
    real(double), pointer :: data_Lold(:,:)
  end type InfoMatrixFile

  character(80),private :: RCSid = "$Id$"

contains
  ! -----------------------------------------------------------------------
  ! Subroutine dump_pos_and_matrices
  ! -----------------------------------------------------------------------

  !!****f* store_matrix/dump_pos_and_matrices *
  subroutine dump_pos_and_matrices(index, MDstep, velocity)
    use global_module, ONLY: ni_in_cell, nspin, nspin_SF, flag_diagonalisation, flag_Multisite, &
                             flag_XLBOMD, flag_propagateX, flag_dissipation, integratorXL, flag_SFcoeffReuse
    use matrix_data, ONLY: Lrange, Hrange, SFcoeff_range, SFcoeffTr_range, HTr_range, Srange, LSrange
    use mult_module, ONLY: matL,L_trans, matK, matSFcoeff, matS
    use io_module, ONLY: append_coords, write_atomic_positions, pdb_template
!    use XLBOMD_module, only: matX, matXvel, dump_XL

    implicit none
    integer,intent(in), optional :: index
    integer,intent(in), optional :: MDstep
    integer :: both=0 , mat=1
    logical :: append_coords_bkup
    integer :: index_local, MDstep_local
    real(double), intent(in), optional :: velocity(1:3,ni_in_cell)
    

    !!! Check whether we should write out the files or not.  !!!
     !   1. check elapsed time
     !   2. check some specific file

    !!! Write out Files to restart..
     !   InfoGlob.dat
     !     PAO coefficients of supports  or Blips
     !     L matrix of K matrix
     !     (for XL-BOMD, files in previous steps should be also printed out)
     !     coodinates file ?

     index_local = 0; if(present(index)) index_local=index
     MDstep_local = 0; if(present(MDstep)) MDstep_local=MDstep

       if(flag_SFcoeffReuse .or. flag_Multisite) then
        call dump_matrix_update('SFcoeff',matSFcoeff,SFcoeff_range,&
                                 index=index_local,nspin=nspin_SF,iprint_mode=mat,MDstep=MDstep_local)
       endif

       if(flag_diagonalisation) then
        if(present(velocity)) then
         call dump_matrix_update('K',matK,Hrange,&
           index=index_local,nspin=nspin,iprint_mode=both,MDstep=MDstep_local,velocity=velocity)
        else
         call dump_matrix_update('K',matK,Hrange,index=index_local,nspin=nspin,iprint_mode=both,MDstep=MDstep_local)
        endif
       else
        if(present(velocity)) then
         call dump_matrix_update('L',matL,Lrange,&
           index=index_local,nspin=nspin, iprint_mode=both,MDstep=MDstep_local,velocity=velocity)
        else
         call dump_matrix_update('L',matL,Lrange,index=index_local,nspin=nspin, iprint_mode=both,MDstep=MDstep_local)
        endif
       endif

! Since XLBOMD_module uses the subroutines in this module, we cannot 
! call subroutines in XLBOMD_module.
!   (Now the following part is in DMMMin_module)
!
!       !XLBOMD for DMM (X=LS) 
!       if(.not.flag_diagonalisation) then
!          if (flag_XLBOMD) then
!             if (flag_propagateX) then
!                call dump_matrix2('X',matX(1),LSrange)
!                if(nspin==2) call dump_matrix2('X_2',matX(2),LSrange)
!                call dump_matrix2('S',matS   ,Srange)
!                if (integratorXL.EQ.'velocityVerlet') then
!                   call dump_matrix2('Xvel',matXvel(1),LSrange)
!                   if(nspin==2) call dump_matrix2('Xvel_2',matXvel(2),LSrange)
!                end if
!             else
!                call dump_matrix2('X',matX(1),Lrange)
!                if(nspin==2) call dump_matrix2('X_2',matX(2),LSrange)
!                if (integratorXL.EQ.'velocityVerlet') then
!                   call dump_matrix2('Xvel',matXvel(1),Lrange)
!                   if(nspin==2) call dump_matrix2('Xvel_2',matXvel(2),LSrange)
!                end if
!             endif
!             ! When dissipation applies
!             if (flag_dissipation) call dump_XL()
!          endif ! (flag_XLBOMD)
!       endif ! (.not.flag_diagonalisation) 

       append_coords_bkup = append_coords; append_coords = .false.
        call write_atomic_positions('coord_next.dat',trim(pdb_template))
       append_coords = append_coords_bkup

   return
  end subroutine dump_pos_and_matrices
  !!***
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
  !!   2019/09/02 TM : nspin is introduced
  !!   
  !!  SOURCE
  !!

   subroutine dump_matrix_update(stub,matA,range,index,nspin,iprint_mode,MDstep,velocity)
    use GenComms, ONLY: inode, ionode, cq_abort
    use global_module, ONLY: numprocs, id_glob, ni_in_cell
    use io_module, ONLY: get_file_name
    implicit none
    character(len=*),intent(in) :: stub
    integer,intent(in) :: matA(:)
    integer,intent(in) :: range
    integer,intent(in), optional :: index
    integer,intent(in), optional :: nspin
    integer,intent(inout), optional :: iprint_mode
    integer,intent(inout), optional :: MDstep
    real(double), intent(in), optional :: velocity(1:3,ni_in_cell)
    integer :: index_local, nspin_local

    index_local=0; if(present(index)) index_local=index
    nspin_local=1; if(present(nspin)) nspin_local=nspin
    if(.not.present(iprint_mode)) iprint_mode = 0

     select case(iprint_mode)
      case(0)  ! Both "InfoGlobal.dat" and "*matrix2.dat" will be pinted out.
       if(present(velocity)) then
        call dump_InfoMatGlobal(index_local,velocity=velocity,MDstep=MDstep)
       else
        call dump_InfoMatGlobal(index_local,MDstep=MDstep)
       endif
        call dump_matrix2(stub,matA,range,nspin_local,index=index_local)
      case(1)  ! only "*matrix2.dat" will be printed out.
        call dump_matrix2(stub,matA,range,nspin_local,index=index_local)
      case(2)  ! only "InfoGlobal.dat" will be printed out.
       if(present(velocity)) then
        call dump_InfoMatGlobal(index=index_local,velocity=velocity,MDstep=MDstep)
       else
        call dump_InfoMatGlobal(index=index_local,MDstep=MDstep)
       endif
     end select ! case(iprint_mode)
 
    return
   end subroutine dump_matrix_update
  !!***
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
  !!   2018/05/29 Tsuyoshi Miyazaki
  !!      introduced nmatrix for the case of dumping multiple matrices having a same cutoff range
  !!   2019/09/02 Tsuyoshi Miyazaki
  !!      nmatrix -> nspin
  !!  SOURCE
  !!

   subroutine dump_matrix2(stub,matA,range,nspin,index)
    use GenComms, ONLY: inode, ionode, cq_abort
    use global_module, ONLY: numprocs, id_glob
    use io_module, ONLY: get_file_name, get_file_name_2rank
    implicit none
    character(len=*),intent(in) :: stub
    integer,intent(in) :: matA(nspin)
    integer,intent(in) :: range
    integer,intent(in) :: nspin
    integer,optional,intent(in) :: index
    type(matrix_store):: tmp_matrix_store

    integer :: lun, iprim, nprim, jmax, jj, ibeg, jbeta_alpha, len, istat
    character(32) :: file_name
    integer :: index_local, nn

    index_local=0; if(present(index)) index_local=index

    ! set_matrix_store : build tmp_matrix_store
    call set_matrix_store(stub,matA,range,nspin,tmp_matrix_store)

    ! Actual Dump (from dump_matrix2)
     ! First, get the name of a file based upon the node ID or rank.
     if(flag_MatrixFile_RankFromZero) then
      call get_file_name_2rank(stub//'matrix2',file_name,index_local,myid)
     else
      call get_file_name_2rank(stub//'matrix2',file_name,index_local,inode)
     endif
     call io_assign(lun)

    ! Case 1 : Binary Format   2018Oct22 TM -----------------------
    if(flag_MatrixFile_BinaryFormat_Dump) then
     open (lun,file=file_name,form='unformatted')

      ! 1. node ID, no. of PS of atoms "i".
        nprim=tmp_matrix_store%n_prim
        write (lun) inode, nprim
      ! 2. no. of alpha for each "i".
        write (lun) tmp_matrix_store%nsf_spec_i(1:nprim)
      ! NOTE: The followings are written out with neighbour-labelling
      ! 3. no. of the neighbours "j" for each "i".
      ! 4. no. of "neighbour-j x beta" for each "i".
        write (lun) tmp_matrix_store%jmax_i(1:nprim)
        write (lun) tmp_matrix_store%jbeta_max_i(1:nprim)
      ! 5. no. of matrices whose elements will be printed out
        write (lun) tmp_matrix_store%nspin

     !I will change the order of dumping in the following, later.   2016/09/30: TM@UCL
      if(nprim .GT. 0) then
       do iprim=1,nprim
          jmax = tmp_matrix_store%jmax_i(iprim)
          ibeg = tmp_matrix_store%ibeg_Rij(iprim)
          write (lun) tmp_matrix_store%idglob_i(iprim)
          write (lun) tmp_matrix_store%beta_j(ibeg:ibeg+jmax-1)
          write (lun) tmp_matrix_store%idglob_j(ibeg:ibeg+jmax-1)
         do jj=1,jmax
          write (lun) tmp_matrix_store%vec_Rij(1:3,ibeg+jj-1)
         enddo !jj=1,jmax
         if(iprim < nprim) then
           len = tmp_matrix_store%ibeg_data_matrix(iprim+1)-tmp_matrix_store%ibeg_data_matrix(iprim)
         else
           len = tmp_matrix_store%matrix_size-tmp_matrix_store%ibeg_data_matrix(iprim)+1
         endif
           ibeg = tmp_matrix_store%ibeg_data_matrix(iprim)
         do nn=1,tmp_matrix_store%nspin
          do jbeta_alpha = 1, len
           write (lun) tmp_matrix_store%data_matrix(ibeg+jbeta_alpha-1,nn)
          enddo !jbeta_alpha = 1, len
         enddo
       enddo !iprim=1,nprim
      endif  ! (nprim .GT. 0)
    
    ! Case 2 : Text Format   (old version) -----------------------
    !                 I assume this part will be removed in the future...
    else
     open (lun,file=file_name,form='formatted')

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
      ! 5. no. of matrices whose elements will be printed out
        write (lun,*) tmp_matrix_store%nspin

     !I will change the order of dumping in the following, later.   2016/09/30: TM@UCL
      if(nprim .GT. 0) then
       do iprim=1,nprim
          jmax = tmp_matrix_store%jmax_i(iprim)
          ibeg = tmp_matrix_store%ibeg_Rij(iprim)
          write (lun,*) tmp_matrix_store%idglob_i(iprim)
          write (lun,*) tmp_matrix_store%beta_j(ibeg:ibeg+jmax-1)
          write (lun,*) tmp_matrix_store%idglob_j(ibeg:ibeg+jmax-1)
         do jj=1,jmax
          write (lun,*) tmp_matrix_store%vec_Rij(1:3,ibeg+jj-1)
         enddo !jj=1,jmax
         if(iprim < nprim) then
           len = tmp_matrix_store%ibeg_data_matrix(iprim+1)-tmp_matrix_store%ibeg_data_matrix(iprim)
         else
           len = tmp_matrix_store%matrix_size-tmp_matrix_store%ibeg_data_matrix(iprim)+1
         endif
           ibeg = tmp_matrix_store%ibeg_data_matrix(iprim)
         do nn=1,tmp_matrix_store%nspin
          do jbeta_alpha = 1, len
           write (lun,fmt='(e30.20)') tmp_matrix_store%data_matrix(ibeg+jbeta_alpha-1,nn)
          enddo !jbeta_alpha = 1, len
         enddo
       enddo !iprim=1,nprim
      endif  ! (nprim .GT. 0) 
    
    endif  
    ! Binary or Text  -----------------------

    ! Close the file in the end.
    call io_close(lun)
    
    ! free_matrix_store : free tmp_matrix_store
    call free_matrix_store(tmp_matrix_store)

    return
   end subroutine dump_matrix2
  !!***

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

   subroutine set_matrix_store(stub,matrices,range,nspin,matinfo)
    use numbers
    use global_module, ONLY: id_glob, species_glob, io_lun, sf, nlpf, atomf
    use species_module, ONLY: nsf_species, natomf_species, nlpf_species
    use group_module, ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module, ONLY: BCS_parts
    use matrix_data, ONLY: mat
    use mult_module, ONLY: mat_p
    use global_module, ONLY: rcellx, rcelly, rcellz  !2017.Dec.14

    implicit none
    character(len=*),intent(in) :: stub
    integer,intent(in) :: nspin
    integer,intent(in) :: matrices(1:nspin)
    integer,intent(in) :: range
    type(matrix_store), intent(out) :: matinfo
    !Local
    integer :: matA   ! used to be passed variable 
    integer :: nprim, istat
    integer :: np, ni, iprim, atom_num, jmax_i_max
    integer :: ist, n_beta, neigh, gcspart, ibeg, j_global_part, len
    integer :: sf1, sf2  ! 2017/Dec/06 for reading SF_coeff
    integer :: ibeg2, jcount, jst ! 2017/Dec/26 for introducing ibeg_Rij
    integer :: nn

    !matA
      matinfo%nspin=nspin
      matA = matrices(1)
    !name & n_prim
      matinfo%name=stub
      matinfo%n_prim=bundle%n_prim
    !allocation 1 
      nprim=matinfo%n_prim   ! local variable
      allocate(matinfo%nsf_spec_i(nprim), matinfo%idglob_i(nprim), matinfo%jmax_i(nprim), &
               matinfo%jbeta_max_i(nprim), matinfo%ibeg_data_matrix(nprim), &
               matinfo%ibeg_Rij(nprim), STAT=istat)
              !matinfo%jbeta_max_i(nprim), matinfo%ibeg_data_matrix(nprim), STAT=istat)
      if(istat .NE. 0) call cq_abort('Allocation 1 in set_matrix_store',istat,nprim)
    
    !sf1_type , sf2_type
    sf1 = mat_p(matA)%sf1_type;  sf2=mat_p(matA)%sf2_type
    if(sf1 == sf) then
      matinfo%nsf_spec_i(1:nprim) =nsf_species(bundle%species(1:nprim))
    elseif(sf1 == atomf) then
      matinfo%nsf_spec_i(1:nprim) =natomf_species(bundle%species(1:nprim))
    elseif(sf1 == nlpf) then
      matinfo%nsf_spec_i(1:nprim) =nlpf_species(bundle%species(1:nprim))
    else
      call cq_abort(" ERROR in set_matrix_store: No sf1_type : ",sf1)
    endif
     

   !!set up : jmax_i, jbeta_max_i

    iprim = 0
    jmax_i_max = -1
    matinfo%jbeta_max_i(:) = 0
    jcount = 0

    do np = 1, bundle%groups_on_node
      if (bundle%nm_nodgroup(np).GT.0) then
        do ni = 1, bundle%nm_nodgroup(np)
           iprim = iprim + 1
           atom_num = bundle%nm_nodbeg(np)+ni-1
           if(iprim .NE. atom_num) call cq_abort('iprim .ne. atom_num ?',iprim,atom_num)
           matinfo%idglob_i(iprim) = bundle%ig_prim(iprim)
           matinfo%jmax_i(iprim) =  mat(np,range)%n_nab(ni)
           jcount = jcount + matinfo%jmax_i(iprim)
           if (matinfo%jmax_i(iprim).GT.jmax_i_max) jmax_i_max = matinfo%jmax_i(iprim)
          do neigh = 1, matinfo%jmax_i(iprim)
            ist = mat(np,range)%i_acc(ni) + neigh - 1 ! Neighbour-labeling
            n_beta = mat(np,range)%ndimj(ist)
            matinfo%jbeta_max_i(iprim) = matinfo%jbeta_max_i(iprim) + n_beta
          enddo !(neigh)
        enddo !(ni)
      endif
    enddo !(np)
    if(jcount < 1 .or. jcount > nprim*jmax_i_max) call cq_abort(' ERROR: # of (j,i) = 0??',jcount)

    !allocation 2 
      allocate(matinfo%idglob_j(jcount), matinfo%beta_j(jcount), STAT=istat)
      if(istat .NE. 0) call cq_abort('Allocation 2 in set_matrix_store', jcount)
    !allocation 3 
      allocate(matinfo%vec_Rij(3,jcount), STAT=istat)
      if(istat .NE. 0) call cq_abort('Allocation 3 in set_matrix_store', jcount)
     
   !!set up : idglob_j, beta_j, vec_Rij
    matinfo%idglob_j(:) = 0
    matinfo%beta_j(:) = 0
    iprim = 0
    ibeg = 1; ibeg2 = 1
    matinfo%matrix_size = 0
    do np = 1, bundle%groups_on_node
      if (bundle%nm_nodgroup(np).GT.0) then
        do ni = 1, bundle%nm_nodgroup(np)
          atom_num=bundle%nm_nodbeg(np)+ni-1
          iprim = iprim + 1
          if(atom_num .NE. iprim) write(io_lun,*) 'Warning!! in set_matrix_store : iprim & atom_num = ',iprim,atom_num

          do neigh = 1, mat(np,range)%n_nab(ni)
            ist = mat(np,range)%i_acc(ni) + neigh - 1  ! Neighbour-labeling
            jst = ibeg2 + neigh -1      

            gcspart=BCS_parts%icover_ibeg(mat(np,range)%i_part(ist))+mat(np,range)%i_seq(ist)-1
            j_global_part = BCS_parts%lab_cell(mat(np,range)%i_part(ist))

            !OLD matinfo%idglob_j(neigh,iprim) = id_glob(parts%icell_beg(j_global_part) &
            matinfo%idglob_j(jst) = id_glob(parts%icell_beg(j_global_part) &
                                                +mat(np,range)%i_seq(ist)-1)
            if(sf2 == sf) then
             matinfo%beta_j(jst) = nsf_species(species_glob(matinfo%idglob_j(jst)))
            elseif(sf2 == atomf) then
             matinfo%beta_j(jst) = natomf_species(species_glob(matinfo%idglob_j(jst)))
            elseif(sf2 == nlpf) then
             matinfo%beta_j(jst) = nlpf_species(species_glob(matinfo%idglob_j(jst)))
            else
             call cq_abort(" ERROR in set_matrix_store: No sf2_type : ",sf2)
            endif

            matinfo%vec_Rij(1,jst) = bundle%xprim(iprim) - BCS_parts%xcover(gcspart)
            matinfo%vec_Rij(2,jst) = bundle%yprim(iprim) - BCS_parts%ycover(gcspart)
            matinfo%vec_Rij(3,jst) = bundle%zprim(iprim) - BCS_parts%zcover(gcspart)

           !vec_Rij is changed from cartesian unit to fractional coordinate  : 2017Dec14
            matinfo%vec_Rij(1,jst) = matinfo%vec_Rij(1,jst)/rcellx
            matinfo%vec_Rij(2,jst) = matinfo%vec_Rij(2,jst)/rcelly
            matinfo%vec_Rij(3,jst) = matinfo%vec_Rij(3,jst)/rcellz

          enddo !neigh = 1, mat(np,range)%n_nab(ni)

          len = matinfo%jbeta_max_i(iprim)*matinfo%nsf_spec_i(iprim)
          matinfo%matrix_size = matinfo%matrix_size + len
          matinfo%ibeg_data_matrix(iprim) = ibeg
          matinfo%ibeg_Rij(iprim) = ibeg2
          ibeg = ibeg + len
          ibeg2 = ibeg2 + matinfo%jmax_i(iprim)
        enddo !(ni)
      endif !(bundle%nm_nodgroup(np).GT.0) then
    enddo !(np)

    ! allocation of the matrix elements 
    !   2019/05/28 TM @ UCL
    !   data_matrix(1:matrix_size) ->  data_matrix(1:matrix_size, 1:nspin)
    allocate(matinfo%data_matrix(matinfo%matrix_size,1:nspin), STAT = istat)
     if(istat .NE. 0) call cq_abort('Allocation 4 in set_matrix_store', istat, matinfo%matrix_size)
    do nn=1, nspin
     matA = matrices(nn)
     matinfo%data_matrix(1:matinfo%matrix_size,nn)=mat_p(matA)%matrix(1:matinfo%matrix_size)
    enddo

    return
   end subroutine set_matrix_store

   subroutine free_matrix_store(matinfo)
    use GenComms, ONLY: cq_abort
    implicit none
    type(matrix_store), intent(inout) :: matinfo
    integer :: istat
   
     !ORI deallocate(matinfo%data_matrix, matinfo%vec_Rij, matinfo%beta_j, matinfo%idglob_j, &
     deallocate(matinfo%data_matrix, matinfo%vec_Rij, matinfo%beta_j, matinfo%idglob_j, matinfo%ibeg_Rij, &
                matinfo%ibeg_data_matrix, matinfo%jbeta_max_i, matinfo%jmax_i, matinfo%idglob_i, &
                matinfo%nsf_spec_i, STAT=istat)
     if(istat .NE. 0) call cq_abort("Error in deallocation at free_matrix_store",istat)
    return
   end subroutine free_matrix_store
  !!***

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
  subroutine dump_InfoMatGlobal(index,velocity,MDstep)

    ! Module usage
    use global_module, ONLY: ni_in_cell,numprocs,rcellx,rcelly,rcellz,id_glob, io_lun
    use GenComms, ONLY: cq_abort, inode, ionode, my_barrier, myid
    use group_module, ONLY: parts
    use io_module, ONLY: get_file_name, get_file_name_2rank

    implicit none
    type(matrix_store_global):: mat_global_tmp
    integer, intent(in), optional :: index
    integer, intent(in), optional :: MDstep
    real(double), intent(in), optional :: velocity(1:3,ni_in_cell)
    integer :: lun, istat, iglob
    integer :: index_local, step_local
    character(len=80) :: filename

    logical :: flag_index, flag_velocity, flag_MDstep

    if(present(index)) then
     flag_index=.true. 
     index_local = index
    else
     flag_index=.false. 
     index_local = 0
    endif

    if(present(MDstep)) then
     flag_MDstep=.true. 
     step_local = MDstep
    else
     flag_MDstep=.false. 
     step_local = 0
    endif

    if(present(velocity)) then
     flag_velocity = .true.
     call set_InfoMatGlobal(mat_global_tmp, step=step_local, velocity_in=velocity)
    else
     flag_velocity = .false.
     call set_InfoMatGlobal(mat_global_tmp, step=step_local)
    endif

    ! Open InfoGlobal.dat and write data.
    if(inode == ionode) then
     call io_assign(lun)
     call get_file_name_2rank('InfoGlobal',filename,index_local)
     write(io_lun,*) ' filename = ',filename
     open (lun,file=filename,iostat=istat)
     rewind lun
     if (istat.GT.0) call cq_abort('Fail in opening InfoGlobal.dat .')
     write (lun,*) flag_velocity, flag_MDstep,'  = flag_velocity, flag_MDstep'
     write (lun,*) mat_global_tmp%ni_in_cell, mat_global_tmp%numprocs, ' # of atoms, # of process '
     write (lun,*) mat_global_tmp%npcellx,mat_global_tmp%npcelly,mat_global_tmp%npcellz,' npcellx,y,z'
     write (lun,fmt='(3f25.15,a)') mat_global_tmp%rcellx,mat_global_tmp%rcelly,mat_global_tmp%rcellz,' rcellx,y,z'
     write (lun,*) mat_global_tmp%glob_to_node(1:mat_global_tmp%ni_in_cell)   
     write (lun,*) mat_global_tmp%MDstep,'  MD step'
 
     do iglob=1, mat_global_tmp%ni_in_cell
      write (lun,101) iglob, mat_global_tmp%atom_coord(1:3, iglob)
      101 format(5x,i8,3x,3e22.13)
     enddo !iglob=1, mat_global_tmp%ni_in_cell

     if(flag_velocity) then
      write (lun,*) 'velocity'
      do iglob=1, mat_global_tmp%ni_in_cell
       write (lun,101) iglob, mat_global_tmp%atom_veloc(1:3, iglob)
      enddo !iglob=1, mat_global_tmp%ni_in_cell
     endif 

     write (lun,*) index_local,' = index_local in dump_InfoMatGlobal'

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
       mat_glob%MDstep      = step
    else
       mat_glob%MDstep      = 0
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
    !if(present(velocity_in)) then
       allocate(mat_glob%glob_to_node(ni_in_cell), mat_glob%atom_coord(3,ni_in_cell), &
               mat_glob%atom_veloc(3,ni_in_cell), STAT=istat)
       if(istat .NE. 0) call cq_abort('Error : allocation in set_InfoMatGlobal1',istat,ni_in_cell)
    !else
    !   allocate(mat_glob%glob_to_node(ni_in_cell), mat_glob%atom_coord(3,ni_in_cell), STAT=istat)
    !   if(istat .NE. 0) call cq_abort('Error : allocation in set_InfoMatGlobal2',istat,ni_in_cell)
    !endif

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
    endif

   return
  end subroutine set_InfoMatGlobal

  !!***
  subroutine free_InfoMatGlobal(mat_glob)
    use GenComms, ONLY: cq_abort
    implicit none
    !logical, intent(in) :: flag_velocity_in
    integer :: istat
    !type(matrix_store_global),intent(out) :: mat_glob
    type(matrix_store_global) :: mat_glob

    !if(flag_velocity_in) then
     deallocate(mat_glob%atom_veloc, mat_glob%atom_coord, mat_glob%glob_to_node, STAT=istat)
     if(istat .NE. 0) call cq_abort('Error : deallocation in free_InfoMatGlobal1',istat)
    !else
    ! deallocate(mat_glob%atom_coord, mat_glob%glob_to_node, STAT=istat)
    ! if(istat .NE. 0) call cq_abort('Error : deallocation in free_InfoMatGlobal2',istat)
    !endif
 
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
  subroutine grab_InfoMatGlobal(InfoGlob,index,flag_velocity)

    ! Module usage
    use GenComms, ONLY: inode,ionode,gcopy, my_barrier
    use io_module, ONLY: get_file_name_2rank
    use global_module, ONLY: io_lun, ni_in_cell

    ! passed variables
    type(matrix_store_global),intent(out) :: InfoGlob
    integer,intent(in), optional :: index
    logical,intent(in), optional :: flag_velocity

    ! local variables
    integer :: lun, istat, iglob, ig
    integer :: index_local, index_in_file
    character(len=80) :: filename
    logical :: flag_velocity_local, flag_MDstep
    integer :: ni_in_cell_tmp, numprocs_tmp

    if(present(flag_velocity)) then
     flag_velocity_local=flag_velocity
    else
     flag_velocity_local=.false.
    endif

    ! check whether members of InfoGlob has been allocated or not.
     if(allocated(InfoGlob%glob_to_node)) call free_InfoMatGlobal(InfoGlob)
     !if(.not.allocated(InfoGlob%glob_to_node)) then
      !if(flag_velocity_local) then
      ! !allocation of glob_to_node, atom_coord, atom_veloc
      !  !write(io_lun,*) "allocation in grab_InfoMatGlobal1"
        allocate(InfoGlob%glob_to_node(ni_in_cell), InfoGlob%atom_coord(3,ni_in_cell), &
                 InfoGlob%atom_veloc(3,ni_in_cell), STAT=istat)
        if(istat .NE. 0) call cq_abort('Error : allocation in grab_InfoMatGlobal1',istat,ni_in_cell)
      !else
      !  !write(io_lun,*) "allocation in grab_InfoMatGlobal2"
      !  allocate(InfoGlob%glob_to_node(ni_in_cell), InfoGlob%atom_coord(3,ni_in_cell), STAT=istat)
      !  if(istat .NE. 0) call cq_abort('Error : allocation in grab_InfoMatGlobal2',istat,ni_in_cell)
      !endif
     !endif

    index_local=0
    if(present(index)) index_local=index

   !Reading  "InfoGlobal.ind***"
    if(inode == ionode) then
     call io_assign(lun)
     call get_file_name_2rank('InfoGlobal',filename,index_local)
     open (lun,file=filename,status='old',iostat=istat)
     if (istat.GT.0) then
      write(io_lun,*) " grab_InfoMatGlobal: Error in opening InfoGlobal: ",filename
      call cq_abort('Fail in opening InfoGlobal.dat')
     endif
 
     !Reading  Start !!!
      read(lun,*) flag_velocity_local, flag_MDstep   ! it is probably better to have this output for the future.
       if(present(flag_velocity)) then
        if(flag_velocity .neqv. flag_velocity_local) call cq_abort('Error in grab_InfoMatglobal: flag_velocity')
       endif
      read(lun,*) InfoGlob%ni_in_cell, InfoGlob%numprocs
       if(ni_in_cell .NE. InfoGlob%ni_in_cell) &
        call cq_abort('Error in grab_InfoMatGlobal: ni_in_cell= ',ni_in_cell,InfoGlob%ni_in_cell)
      read(lun,*) InfoGlob%npcellx,InfoGlob%npcelly,InfoGlob%npcellz
      read(lun,*) InfoGlob%rcellx,InfoGlob%rcelly,InfoGlob%rcellz
      read(lun,*) InfoGlob%glob_to_node(1:InfoGlob%ni_in_cell)
      read(lun,*) InfoGlob%MDstep
 
     do ig=1, InfoGlob%ni_in_cell
      read(lun,101) iglob, InfoGlob%atom_coord(1:3, ig)
      101 format(5x,i8,3x,3e22.13)
     enddo !ig=1, InfoGlob%ni_in_cell

     if(flag_velocity_local) then
      !if(.not.allocated(InfoGlob%atom_veloc)) then
      ! allocate(InfoGlob%atom_veloc(3,ni_in_cell),stat=istat) 
      ! if(istat .NE. 0) call cq_abort('Error : allocation3 in grab_InfoMatGlobal1',istat,ni_in_cell)
      !endif
      write(io_lun,*) 'Reading velocity from ',filename
      read(lun,*)
      do ig=1, InfoGlob%ni_in_cell
       read(lun,101) iglob, InfoGlob%atom_veloc(1:3, ig)
      enddo !ig=1, InfoGlob%ni_in_cell
     endif

     read(lun,*) index_in_file
     if(index_in_file .ne. index) then
      write(io_lun, *) 'ERROR: Index in the file InfoGlobal.*** is different. ',filename,index_in_file
      call cq_abort(' ERROR in Reading InfoGlobal ')
     endif

     call io_close (lun)
    endif !(inode == ionode) 

   !flag_velocity_local
    call gcopy(flag_velocity_local)

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
    call gcopy(InfoGlob%MDstep)
    call gcopy(InfoGlob%atom_coord, 3, ni_in_cell)
    if(flag_velocity_local) call gcopy(InfoGlob%atom_veloc, 3, ni_in_cell)

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
  !!***

  ! -----------------------------------------------------------------------
  ! Subroutine grab_matrix2
  ! -----------------------------------------------------------------------

  !!****f* store_matrix/grab_matrix2 *
  !!
  !!  NAME
  !!   grab_matrix2
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Reads all the necessary data for matrix reconstruction
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/08/21
  !!  MODIFICATION
  !!   2016/04/06 dave
  !!    Changed Info type from allocatable to pointer (fix gcc 4.4.7 compile issue)
  !!   2017/05/09 dave
  !!    Removed attempts to read both spin channels in one routine (this routine should now be called once for each channel)
  !!   2017/11/10 tsuyoshi 
  !!    Moved from io_module2 to store_matrix
  !!  SOURCE
  !!
  subroutine grab_matrix2(stub,inode,nfile,Info,index,nspin_in)

    ! Module usage
    use io_module, ONLY: get_file_name, get_file_name_2rank
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: io_lun,n_proc_old

    implicit none

    ! passed variables
    integer :: max_node  ! No. of nodes in the PREVIOUS job.
    integer :: inode, matA
    character(len=*) :: stub
    type(InfoMatrixFile), pointer :: Info(:)
    integer, optional :: index
    integer, optional :: nspin_in

    ! local variables
    integer :: lun,stat,padzeros,stat_alloc,size,size2,sizeL,i,j,jbeta_alpha,len,ifile,ibeg
    integer :: nfile, nrest_file, index_file
    integer :: proc_id, jmax_i_max,ios, index_local, nspin_local
    character(32) :: file_name
    character(80) :: num

    ! db
    integer :: lun_db
    character(20) :: file_db

    max_node = n_proc_old

    index_local=0; if(present(index)) index_local=index
    nspin_local=1; if(present(nspin_in)) nspin_local=nspin_in

    ! Open matrix files to get "natom_i", allocate the arrays of Info,
    ! and read the data.
    ! Note numprocs = No. of MPI_process_new and max_nodes = No. of MPI_process_old.

    nfile = int(max_node/numprocs)
    nfile = int(max_node/numprocs)
    nrest_file = mod(max_node,numprocs)
    if (inode.LT.nrest_file+1) nfile = nfile + 1

    ! Allocate Info when nfile is greater than 0.
    if (nfile.GT.0) then
      allocate (Info(nfile), STAT=stat_alloc)
      if (stat_alloc.NE.0) call cq_abort('Error allocating Info: ', nfile)
      do ifile = 1, nfile
        ! Open a file to start with.
       if(flag_MatrixFile_RankFromZero) then
        index_file=numprocs*(ifile-1)+myid
       else
        index_file=numprocs*(ifile-1)+inode
       endif

        call get_file_name_2rank(stub//'matrix2',file_name,index_local,index_file)

        call io_assign(lun)
    
       ! Case 1 : Binary Format   2018Oct22 TM -----------------------
       if(flag_MatrixFile_BinaryFormat_Grab) then
        open (lun,file=file_name,status='old',iostat=stat,form='unformatted')
        if (stat.NE.0) call cq_abort('Fail in opening Lmatrix file in Binary.')
        ! Get dimension, allocate arrays and store necessary data.
        read (lun,iostat=stat) proc_id, Info(ifile)%natom_i
         if(stat.NE.0) then
          write(io_lun,*) " ERROR in reading Binary File of Lmatrix : inode= ", inode
          write(io_lun,*) "  ** Set IO.MatirxFile.BinaryFormat or IO.MatrixFile.BinaryFormat.Grab as False.&
                          & if you use ASCII for Matrix Files."
          call cq_abort('Fail in reading Lmatrix File in Binary')
         endif
        size = Info(ifile)%natom_i
        allocate (Info(ifile)%alpha_i(size), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating alpha_i:', size)
        allocate (Info(ifile)%idglob_i(size), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating idglob_i:', size)
        allocate (Info(ifile)%jmax_i(size), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating jmax_i:', size)
        allocate (Info(ifile)%jbeta_max_i(size), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating jbeta_max_i:', size)
        allocate (Info(ifile)%ibeg_dataL(size), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating ibeg_dataL:', size)
        allocate (Info(ifile)%ibeg_Pij(size), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating ibeg_Pij:', size)
        read (lun) Info(ifile)%alpha_i(1:size)
        read (lun) Info(ifile)%jmax_i(1:size)
        read (lun) Info(ifile)%jbeta_max_i(1:size)
        jmax_i_max = -1
        do i = 1, size
          if (Info(ifile)%jmax_i(i).GT.jmax_i_max) jmax_i_max = Info(ifile)%jmax_i(i)
        enddo
        !db write (io_lun,*) "jmax_i_max: ", jmax_i_max
        size2 = 0
        do i = 1, size
          size2 = size2 + Info(ifile)%jmax_i(i)
        enddo
        allocate (Info(ifile)%idglob_j(size2), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating idglob_j:', size2)
        allocate (Info(ifile)%beta_j_i(size2), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating beta_j_i:', size2)
        allocate (Info(ifile)%rvec_Pij(3,size2), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating rvec_Pij:',size2)
        sizeL = 0
        do i = 1, size
          sizeL = sizeL + Info(ifile)%alpha_i(i)*Info(ifile)%jbeta_max_i(i)
        enddo
        read (lun) Info(ifile)%nspin
        if(Info(ifile)%nspin .ne. nspin_local) call cq_abort("ERROR: nspin mismatch ", Info(ifile)%nspin, nspin_local)
        allocate (Info(ifile)%data_Lold(sizeL,Info(ifile)%nspin), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating data_Lold:', sizeL, Info(ifile)%nspin)

        Info(ifile)%ibeg_dataL(1) = 1 ; Info(ifile)%ibeg_Pij(1) = 1
        ibeg = 1
        do i = 1, size
          read (lun) Info(ifile)%idglob_i(i)
          read (lun) Info(ifile)%beta_j_i(Info(ifile)%ibeg_Pij(i) : &
                                             Info(ifile)%ibeg_Pij(i)+Info(ifile)%jmax_i(i)-1)
          read (lun) Info(ifile)%idglob_j(Info(ifile)%ibeg_Pij(i) : &
                                             Info(ifile)%ibeg_Pij(i)+Info(ifile)%jmax_i(i)-1)
          do j = 1, Info(ifile)%jmax_i(i)
            read (lun) Info(ifile)%rvec_Pij(1:3,Info(ifile)%ibeg_Pij(i)+j-1)
          enddo
          len = Info(ifile)%jbeta_max_i(i)*Info(ifile)%alpha_i(i)

          ! spin
          if (Info(ifile)%nspin.EQ.1) then
            do jbeta_alpha = 1, len
              read (lun) Info(ifile)%data_Lold(Info(ifile)%ibeg_dataL(i)+jbeta_alpha-1, 1)
            enddo
          elseif (Info(ifile)%nspin.EQ.2) then
            do jbeta_alpha = 1, len
              read (lun,*) Info(ifile)%data_Lold(Info(ifile)%ibeg_dataL(i)+jbeta_alpha-1, 2)
            enddo
          endif

          if (i+1.LE.size) then
            Info(ifile)%ibeg_Pij(i+1) = Info(ifile)%ibeg_Pij(i) + Info(ifile)%jmax_i(i)
            Info(ifile)%ibeg_dataL(i+1) = Info(ifile)%ibeg_dataL(i) + len
          endif
        enddo !(i, size)

       else
       ! Case 2 : Text Format   old version  -----------------------
        open (lun,file=file_name,status='old',iostat=stat)
            write(*,*) " Trying to open the ASCII file : ",file_name
        if (stat.NE.0) call cq_abort('Fail in opening ASCII Lmatrix file.')
        ! Get dimension, allocate arrays and store necessary data.
        read (lun,*,iostat=stat) proc_id, Info(ifile)%natom_i
         if(stat.NE.0) then
          write(io_lun,*) " ERROR in reading ASCII File of Lmatrix : inode= ", inode
          write(io_lun,*) "  ** Set IO.MatirxFile.BinaryFormat or IO.MatrixFile.BinaryFormat.Grab as .True.&
                          & if you use Binary for Matrix Files."
          call cq_abort('Fail in reading Lmatrix File in ASCII')
         endif
        size = Info(ifile)%natom_i
        allocate (Info(ifile)%alpha_i(size), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating alpha_i:', size)
        allocate (Info(ifile)%idglob_i(size), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating idglob_i:', size)
        allocate (Info(ifile)%jmax_i(size), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating jmax_i:', size)
        allocate (Info(ifile)%jbeta_max_i(size), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating jbeta_max_i:', size)
        allocate (Info(ifile)%ibeg_dataL(size), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating ibeg_dataL:', size)
        allocate (Info(ifile)%ibeg_Pij(size), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating ibeg_Pij:', size)
        read (lun,*) Info(ifile)%alpha_i(1:size)
        read (lun,*) Info(ifile)%jmax_i(1:size)
        read (lun,*) Info(ifile)%jbeta_max_i(1:size)
        jmax_i_max = -1
        do i = 1, size
          if (Info(ifile)%jmax_i(i).GT.jmax_i_max) jmax_i_max = Info(ifile)%jmax_i(i)
        enddo
        !db write (io_lun,*) "jmax_i_max: ", jmax_i_max
        size2 = 0
        do i = 1, size
          size2 = size2 + Info(ifile)%jmax_i(i)
        enddo
        allocate (Info(ifile)%idglob_j(size2), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating idglob_j:', size2)
        allocate (Info(ifile)%beta_j_i(size2), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating beta_j_i:', size2)
        allocate (Info(ifile)%rvec_Pij(3,size2), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating rvec_Pij:',size2)
        sizeL = 0
        do i = 1, size
          sizeL = sizeL + Info(ifile)%alpha_i(i)*Info(ifile)%jbeta_max_i(i)
        enddo

        read (lun,*) Info(ifile)%nspin
        if(Info(ifile)%nspin .ne. nspin_local) call cq_abort("ERROR: nspin mismatch ", Info(ifile)%nspin, nspin_local)
        allocate (Info(ifile)%data_Lold(sizeL,Info(ifile)%nspin), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating data_Lold:', sizeL, Info(ifile)%nspin)

        Info(ifile)%ibeg_dataL(1) = 1 ; Info(ifile)%ibeg_Pij(1) = 1
        ibeg = 1
        do i = 1, size
          read (lun,*) Info(ifile)%idglob_i(i)
          read (lun,*) Info(ifile)%beta_j_i(Info(ifile)%ibeg_Pij(i) : &
                                             Info(ifile)%ibeg_Pij(i)+Info(ifile)%jmax_i(i)-1)
          read (lun,*) Info(ifile)%idglob_j(Info(ifile)%ibeg_Pij(i) : &
                                             Info(ifile)%ibeg_Pij(i)+Info(ifile)%jmax_i(i)-1)
          do j = 1, Info(ifile)%jmax_i(i)
            read (lun,*) Info(ifile)%rvec_Pij(1:3,Info(ifile)%ibeg_Pij(i)+j-1)
          enddo
          len = Info(ifile)%jbeta_max_i(i)*Info(ifile)%alpha_i(i)

          ! spin
          if (Info(ifile)%nspin.EQ.1) then
            do jbeta_alpha = 1, len
              read (lun,*) Info(ifile)%data_Lold(Info(ifile)%ibeg_dataL(i)+jbeta_alpha-1, 1)
            enddo
          elseif (Info(ifile)%nspin.EQ.2) then
            do jbeta_alpha = 1, len
              read (lun,*) Info(ifile)%data_Lold(Info(ifile)%ibeg_dataL(i)+jbeta_alpha-1, 2)
            enddo
          endif

          if (i+1.LE.size) then
            Info(ifile)%ibeg_Pij(i+1) = Info(ifile)%ibeg_Pij(i) + Info(ifile)%jmax_i(i)
            Info(ifile)%ibeg_dataL(i+1) = Info(ifile)%ibeg_dataL(i) + len
          endif
        enddo !(i, size)

       endif
       !Binary or Text -------

        ! Close the file.
       call io_close (lun)
      enddo !(ifile, nfile)
    endif

    return
  end subroutine grab_matrix2
  !!***

  !!****f* store_matrix/deallocate_InfoMatrixFile *
  !!
  !!  NAME
  !!   deallocate_InfoLmatrix_File
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Releases memories for the derived data type InfoMatrixFile(:)
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/08/21
  !!  MODIFICATION
  !!   2016/04/06 dave
  !!    Changed Info type from allocatable to pointer (fix gcc 4.4.7 compile issue)
  !!   2017/11/10 tsuyoshi
  !!    Moved from io_module2
  !!
  !!  SOURCE
  !!
  subroutine deallocate_InfoMatrixFile(nfile,Info)

    ! Module usage
    use numbers
    use GenComms, ONLY: cq_abort

    ! passed variables
    integer :: nfile
    type(InfoMatrixFile), pointer :: Info(:)

    ! local variables
    integer :: ifile,stat_alloc

    if (associated(Info)) then
      do ifile = 1, nfile
        deallocate (Info(ifile)%alpha_i, STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error deallocating alpha_i:')
        deallocate (Info(ifile)%idglob_i, STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error deallocating idglob_i:')
        deallocate (Info(ifile)%jmax_i, STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error deallocating jmax_i:')
        deallocate (Info(ifile)%jbeta_max_i, STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error deallocating jbeta_max_i:')
        deallocate (Info(ifile)%ibeg_Pij, STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error deallocating ibeg_Pij:')
        deallocate (Info(ifile)%ibeg_dataL, STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error deallocating ibeg_dataL:')
        deallocate (Info(ifile)%beta_j_i, STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error deallocating beta_j_i:')
        deallocate (Info(ifile)%idglob_j, STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error deallocating idglob_j:')
        deallocate (Info(ifile)%rvec_Pij, STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error deallocating rvec_Pij:')
        deallocate (Info(ifile)%data_Lold, STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error deallocating data_Lold:')
      enddo
      deallocate (Info, STAT=stat_alloc)
      if (stat_alloc.NE.0) call cq_abort('Error deallocating Info:', nfile)
    endif

    return
  end subroutine deallocate_InfoMatrixFile
  !!***

end module store_matrix
