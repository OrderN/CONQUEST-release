!!****h* Conquest/io_module2
!!  NAME
!!   io_module2
!!  PURPOSE
!!   Grabs and writes out the necessary information used in
!!   matrix reconstruction
!!  USES
!!
!!  AUTHOR
!!   Michiaki Arita
!!  CREATION DATE
!!   2013/08/21
!!  MODIFCATION
!!   2013/12/02 M.Arita
!!   - Added derived data types for XL-BOMD
!!   2016/04/06 dave
!!    Changed Info type from allocatable to pointer (fix gcc 4.4.7 compile issue)
!!
module io_module2

  use datatypes
  use global_module, ONLY: numprocs,flag_MDdebug,iprint_MDdebug
  use input_module, ONLY: io_assign, io_close
  use GenComms, ONLY: inode, ionode, cq_abort

  implicit none

  type InfoMatrixFile
    integer :: index_node
    integer :: natom_i
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
    !OLD real(double), pointer :: rvec_Pij(:,:,:)
    ! Size : 
    !ORI real(double), pointer :: data_Lold(:)
    real(double), pointer :: data_Lold(:,:)
  end type InfoMatrixFile

  type(InfoMatrixFile), pointer :: InfoL(:),InfoT(:),InfoK(:),    &
                                       InfoX(:),InfoXvel(:),InfoS(:), & ! for XL-BOMD
                                       InfoX1(:),InfoX2(:),InfoX3(:), & ! for dissipation
                                       InfoX4(:),InfoX5(:),InfoX6(:), &
                                       InfoX7(:),InfoX8(:),InfoX9(:), &
                                       InfoX10(:)
  integer :: n_matrix

  character(80),private :: RCSid = "$Id$"

contains

  ! -----------------------------------------------------------------------
  ! Subroutine dump_matrix2
  ! -----------------------------------------------------------------------

  !!****f* io_module2/dump_matrix2 *
  !!
  !!  NAME
  !!    dump_matrix2
  !!  USAGE
  !!
  !!  PURPOSE
  !!    Stores the information as for the primary atoms and its
  !!    neighbours first, then the matrix elements are followed.
  !!    Note that matrix data are the ones obtained at the previous (old)
  !!    step - e.g. L_{i.OLD,alpha | j.OLD,beta}.
  !!
  !!    imax:		the number of primary atoms "i" in the current node
  !!			"inode".
  !!
  !!    alpha_i():	the number of orbitals for each primary set of atom "i"
  !!			in the current node "inode".
  !!
  !!    jmax():		the number of neighbour-atoms "j" for each atom "i" above.
  !!
  !!    jbeta_max():	the number of (jmax x beta) for each atom "i".
  !!
  !!    beta_j():	the number of orbitals for each neighbour-atom "j".
  !!
  !!    idglob_j():	the global labelling of neighbour-atoms "j".
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/08/21
  !!  MODIFICATION
  !!
  !!  SOURCE
  !!
  !subroutine dump_matrix2(stub,inode,range,matA,matB)
!  subroutine dump_matrix2(stub,matA,inode,range)
!
!    ! Module usage
!    use numbers
!    use basic_types, ONLY: cover_set
!    use global_module, ONLY: id_glob, species_glob, atom_coord
!    use io_module, ONLY: get_file_name
!    use group_module, ONLY: parts
!    use mult_module, ONLY: return_matrix_value, matrix_index, mat_p
!    use primary_module, ONLY: bundle
!    use cover_module, ONLY: BCS_parts
!    use matrix_data, ONLY: Lrange, mat
!    use species_module, ONLY: nsf_species
!    use matrix_data, ONLY: matrix_pointer, mat, Lrange
!
!    ! passed variables
!    integer :: inode, matA
!    integer,intent(in) :: range
!    character(len=*) :: stub
!    !integer,optional :: matB
!
!    ! local variables
!    integer :: lun, stat_alloc
!    integer :: iprim, ipart, mem_ipart, nab, n_beta, ist, j, jmax_i_max, &
!               j_global_part, np, ni, jsf, neigh, isf, ibeg, jbeta_alpha, len, gcspart, &
!               ipartx_ps,iparty_ps,ipartz_ps,ipartx_cs,iparty_cs,ipartz_cs,CC, &
!               ni_cover,idlocal_i,i_in_cover,np_cover,np_cover_CC,x_CC,y_CC,z_CC,atom_num
!    integer, allocatable :: jmax_i(:), jbeta_max_i(:), j_global_num(:,:), &
!                            beta_j(:,:)
!    real(double), allocatable :: vec_Rij(:,:), Ltmp(:)
!    character(20) :: file_name
!    logical :: find_local
!
!    ! DEBUG
!    integer :: lun_db
!    real(double) :: dist
!    character(20) :: file_db
!
!
!    ! Firstly, make a file based upon the node ID.
!    ! Its name will be renamed later.
!    !ORI call get_file_name('Lmatrix2',numprocs,inode,file_name)
!    call get_file_name(stub//'matrix2',numprocs,inode,file_name)
!    call io_assign(lun)
!    open (lun,file=file_name)
!
!    !! -------- DEBUG -------- !!
!    if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
!      call get_file_name('io2',numprocs,inode,file_db)
!      call io_assign(lun_db)
!      !ORI open (lun_db,file=file_db,position='append')
!      open (lun_db,file=file_db)
!      write (lun_db,*) "mat_p(matL(1))%matrix(1:)"
!      write (lun_db,*) mat_p(matA)%matrix(1:)
!    endif
!    !! -------- DEBUG -------- !!
!
!    ! Then, write down the information necessary for reading Lmatrix files.
!    ! 1. node ID, no. of PS of atoms "i".
!    write (lun,*) inode, bundle%n_prim 
!
!    ! 2. no. of alpha for each "i".
!    write (lun,*) nsf_species(bundle%species(1:bundle%n_prim))
!
!    ! NOTE: The followings are written out with neighbour-labelling
!    ! 3. no. of the neighbours "j" for each "i".
!    ! 4. no. of "neighbour-j x beta" for each "i".
!    allocate (jmax_i(bundle%n_prim),jbeta_max_i(bundle%n_prim), STAT=stat_alloc)
!    if (stat_alloc.NE.0) call cq_abort('Error allocating jmax_i, jbeta_max_i: ', &
!                                        bundle%n_prim)
!    iprim = 0
!    jbeta_max_i = 0
!    jmax_i_max = -1
!    do np = 1, bundle%groups_on_node
!      if (bundle%nm_nodgroup(np).GT.0) then
!        do ni = 1, bundle%nm_nodgroup(np)
!          iprim = iprim + 1
!          atom_num = bundle%nm_nodbeg(np)+ni-1
!          jmax_i(atom_num) =  mat(np,range)%n_nab(ni)
!          if (jmax_i(atom_num).GT.jmax_i_max) jmax_i_max = jmax_i(atom_num)
!          do neigh = 1, jmax_i(atom_num)
!            ist = mat(np,range)%i_acc(ni) + neigh - 1 ! Neighbour-labeling
!            n_beta = mat(np,range)%ndimj(ist)
!            jbeta_max_i(iprim) = jbeta_max_i(iprim) + n_beta
!          enddo !(neigh)
!
!          !! ---- DEBUG ---- !!
!          if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
!            write (lun_db,*) "iprim,np,ni: ", iprim, np, ni
!            write (lun_db,*) ""
!          endif
!          !! ---- DEBUG ---- !!
!
!        enddo !(ni)
!      endif
!    enddo !(np)
!    write (lun,*) jmax_i(1:bundle%n_prim)
!    write (lun,*) jbeta_max_i(1:bundle%n_prim)
!
!    ! 5. beta, j in global ID, vec_Rij and L for each "j" over "i".
!    ! It does not look good to allocate arrays with jmax_i_max, 
!    ! but I hate to allocate them in a do-loop.
!    ! Be careful - Spin is NOT considered.
!    allocate (j_global_num(bundle%n_prim,jmax_i_max), &
!              beta_j(bundle%n_prim,jmax_i_max), STAT=stat_alloc)
!    if (stat_alloc.NE.0) call cq_abort('Error allocating j_global_num: ')
!    allocate (vec_Rij(3,jmax_i_max), STAT=stat_alloc)
!    if (stat_alloc.NE.0) call cq_abort('Error allocating vec_Rij: ',3,jmax_i_max)
!
!    j_global_num = 0
!    beta_j = 0
!    iprim = 0
!    ibeg = 1
!    do np = 1, bundle%groups_on_node
!      if (bundle%nm_nodgroup(np).GT.0) then
!        do ni = 1, bundle%nm_nodgroup(np)
!          atom_num=bundle%nm_nodbeg(np)+ni-1
!          iprim = iprim + 1
!  
!          !! ------------ DEBUG: ------------ !!
!          if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
!            write (lun_db,*) ""
!            write (lun_db,*) "parts. & members in bundle:", np, ni
!            write (lun_db,*) "No. of neighbour:", mat(np,range)%n_nab(ni)
!            write (lun_db,*) ""
!          endif
!          !! ------------ DEBUG: ------------ !!
!
!          !ORI do neigh = 1, jmax_i(iprim)
!          do neigh = 1, mat(np,range)%n_nab(ni)
!            ist = mat(np,range)%i_acc(ni) + neigh - 1  ! Neighbour-labeling
!            gcspart=BCS_parts%icover_ibeg(mat(np,range)%i_part(ist))+mat(np,range)%i_seq(ist)-1
!            j_global_part = BCS_parts%lab_cell(mat(np,range)%i_part(ist))
!
!            !! --------- DEBUG: --------- !!
!            if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
!              write (lun_db,*) "j_global_part:", j_global_part
!              write (lun_db,*) "parts%icell_beg:", parts%icell_beg(j_global_part)
!              write (lun_db,*) "mat(np,range)i_seq(ist)-1:", mat(np,range)%i_seq(ist)-1
!            endif
!            !! --------- DEBUG: --------- !!
!
!            j_global_num(atom_num,neigh) = id_glob(parts%icell_beg(j_global_part) &
!                                                +mat(np,range)%i_seq(ist)-1)
!            beta_j(atom_num,neigh) = nsf_species(species_glob(j_global_num(atom_num,neigh)))
!            vec_Rij(1,neigh) = bundle%xprim(atom_num) - BCS_parts%xcover(gcspart)
!            vec_Rij(2,neigh) = bundle%yprim(atom_num) - BCS_parts%ycover(gcspart)
!            vec_Rij(3,neigh) = bundle%zprim(atom_num) - BCS_parts%zcover(gcspart)
!
!            !! -------- DEBUG: -------- !!
!            if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
!              dist = vec_Rij(1,neigh)*vec_Rij(1,neigh)+vec_Rij(2,neigh)*vec_Rij(2,neigh) &
!                      +vec_Rij(3,neigh)*vec_Rij(3,neigh)
!              dist = sqrt(dist)
!              write (lun_db,*) "Neigh-label, part(in BCS, CC) &  dist.:"
!              write (lun_db,*) neigh,BCS_parts%lab_cover(mat(np,range)%i_part(ist)),dist
!            endif
!            !! -------- DEBUG: -------- !!
!
!          enddo
!          write (lun,*) bundle%ig_prim(iprim)
!          write (lun,*) beta_j(iprim,1:jmax_i(atom_num))
!          write (lun,*) j_global_num(iprim,1:jmax_i(atom_num))
!          do j = 1, jmax_i(atom_num)
!            !DEBUG write (lun,*) j, vec_Rij(1:3,j)
!            write (lun,*) vec_Rij(1:3,j)
!          enddo
!          len = jbeta_max_i(iprim)*nsf_species(bundle%species(iprim))
!
!          !! ---------- DEBUG ---------- !!
!          if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
!            write (lun_db,*) "mat_p(matL(1))%length: ", mat_p(matA)%length 
!            write (lun_db,*) "ni,len: ",ni,len
!            write (lun_db,*) "ibeg:", ibeg
!          endif
!          !! ---------- DEBUG ---------- !!
!          
!          do jbeta_alpha = 1, len
!            write (lun,fmt='(f25.18)') mat_p(matA)%matrix(ibeg+jbeta_alpha-1)
!          enddo
!          ibeg = ibeg + len
!        enddo !(ni)
!      endif
!    enddo !(np)
!
!    !! ---------- DEBUG: ---------- !!
!    if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
!      write (lun_db,'(/a)') "--- Lmatrix elements: mat_p(matA)%matrix(1:length)"
!      do ni = 1, mat_p(matA)%length
!        write (lun_db,'(f25.18)') mat_p(matA)%matrix(ni)
!      enddo
!    endif
!    !! ---------- DEBUG: ---------- !!
!
!    ! Deallocation section.
!    deallocate (jmax_i,jbeta_max_i, STAT=stat_alloc)
!    if (stat_alloc.NE.0) call cq_abort('Error deallocating jmax_i, jbeta_max_i: ', &
!                                        bundle%n_prim)
!    deallocate (j_global_num,beta_j, STAT=stat_alloc)
!    if (stat_alloc.NE.0) call cq_abort('Error deallocating j_global_num:')
!    deallocate (vec_Rij, STAT=stat_alloc)
!    if (stat_alloc.NE.0) call cq_abort('Error deallocating vec_Rij:')
!
!    ! Close the file in the end.
!    call io_close(lun)
!
!    !! ---------- DEBUG ---------- !!
!    if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
!      write (lun_db,*) ""
!      write (lun_db,*) ""
!      call io_close(lun_db)
!    endif
!    !! ---------- DEBUG ---------- !!
!
!    return
!  end subroutine dump_matrix2
!  !!***
!
  ! -----------------------------------------------------------------------
  ! Subroutine grab_matrix2
  ! -----------------------------------------------------------------------

  !!****f* io_module2/grab_matrix2 *
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
  !!  SOURCE
  !!
  subroutine grab_matrix2(stub,inode,nfile,Info)

    ! Module usage
    use io_module, ONLY: get_file_name
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: nspin,io_lun,n_proc_old

    implicit none

    ! passed variables
    integer :: max_node  ! No. of nodes in the PREVIOUS job.
    integer :: inode, matA
    character(len=*) :: stub
    type(InfoMatrixFile), pointer :: Info(:)

    ! local variables
    integer :: lun,stat,padzeros,stat_alloc,size,size2,sizeL,i,j,jbeta_alpha,len,ifile,ibeg
    integer :: nfile, nrest_file, index_file
    integer :: proc_id, jmax_i_max,ios
    character(20) :: file_name  
    character(80) :: num

    ! db
    integer :: lun_db
    character(20) :: file_db

    max_node = n_proc_old

    !! ---- DEBUG ---- !!
    if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
      call get_file_name('ReadLmatrix2',max_node,inode,file_db)
      call io_assign(lun_db)
      open (lun_db,file=file_db,status='old',iostat=ios)
      if(ios /= 0) call cq_abort('grab_matrix2: failed to open input file '//file_db)
    endif
    !! ---- DEBUG ---- !!

    ! Open matrix files to get "natom_i", allocate the arrays of Info, 
    ! and read the data.
    ! Note numprocs = No. of MPI_process_new and max_nodes = No. of MPI_process_old.
    nfile = int(max_node/numprocs)
    nfile = int(max_node/numprocs)
    nrest_file = mod(max_node,numprocs)
    if (inode.LT.nrest_file+1) nfile = nfile + 1
    !db write (io_lun,*) "inode; nfile, nrest_file: ", inode, nfile, nrest_file
    !db write (io_lun,*) "Info allocated?: ", allocated(Info)
    ! Allocate Info when nfile is greater than 0.
    if (nfile.GT.0) then
      allocate (Info(nfile), STAT=stat_alloc)
      if (stat_alloc.NE.0) call cq_abort('Error allocating Info: ', nfile)
      do ifile = 1, nfile
        ! Open a file to start with.
        index_file=numprocs*(ifile-1)+inode
        padzeros=MAX(3, &
                     FLOOR(1.0+LOG10(REAL(numprocs,kind=double)))) - &
                     FLOOR(1.0+LOG10(REAL(inode,kind=double))) 
        write (num,'(i80)') index_file
        num=adjustl(num)
        !ORI file_name='Lmatrix2.'
        file_name=stub//'matrix2.'
        do i = 1, padzeros
          file_name=trim(file_name)//'0'
        enddo
        file_name=trim(file_name)//num
        !db write (io_lun,*) "inode, file_name:", inode, file_name
        call io_assign(lun)
        open (lun,file=file_name,status='old',iostat=stat)
        if (stat.NE.0) call cq_abort('Fail in opening Lmatrix file.')
        ! Get dimension, allocate arrays and store necessary data.
        read (lun,*) proc_id, Info(ifile)%natom_i
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
        ! n_matrix depends on nspin
        n_matrix = 1
        !if (nspin.EQ.2) n_matrix = 2
        !ORI allocate (Info(ifile)%data_Lold(sizeL), STAT=stat_alloc)
        allocate (Info(ifile)%data_Lold(sizeL,n_matrix), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating data_Lold:', sizeL, n_matrix)

        !! ---- DEBUG ---- !!
        if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
          write (lun_db,*) 'n_proc_old:', max_node
          write (lun_db,*) proc_id, size
          write (lun_db,*) Info(ifile)%alpha_i(1:size)
          write (lun_db,*) Info(ifile)%jmax_i(1:size)
          write (lun_db,*) Info(ifile)%jbeta_max_i(1:size)
          !write (lun_db,*) "jmax_i_max:", jmax_i_max
        endif
        !! ---- DEBUG ---- !!

        Info(ifile)%ibeg_dataL(1) = 1 ; Info(ifile)%ibeg_Pij(1) = 1
        ibeg = 1
        do i = 1, size
          read (lun,*) Info(ifile)%idglob_i(i)
          if (flag_MDdebug .AND. iprint_MDdebug.GT.3) &
            write (lun_db,*) Info(ifile)%idglob_i(i)
          read (lun,*) Info(ifile)%beta_j_i(Info(ifile)%ibeg_Pij(i) : &
                                             Info(ifile)%ibeg_Pij(i)+Info(ifile)%jmax_i(i)-1)
          if (flag_MDdebug .AND. iprint_MDdebug.GT.3) &
            write (lun_db,*) Info(ifile)%beta_j_i(Info(ifile)%ibeg_Pij(i) : &
                             Info(ifile)%ibeg_Pij(i)+Info(ifile)%jmax_i(i)-1)
          read (lun,*) Info(ifile)%idglob_j(Info(ifile)%ibeg_Pij(i) : & 
                                             Info(ifile)%ibeg_Pij(i)+Info(ifile)%jmax_i(i)-1)
          if (flag_MDdebug .AND. iprint_MDdebug.GT.3) &
            write (lun_db,*) Info(ifile)%idglob_j(ibeg:ibeg+Info(ifile)%jmax_i(i)-1)
          do j = 1, Info(ifile)%jmax_i(i)
            read (lun,*) Info(ifile)%rvec_Pij(1:3,Info(ifile)%ibeg_Pij(i)+j-1)
            if (flag_MDdebug .AND. iprint_MDdebug.GT.3) &
              write (lun_db,*) Info(ifile)%rvec_Pij(1:3,Info(ifile)%ibeg_Pij(i)+j-1)
          enddo
          len = Info(ifile)%jbeta_max_i(i)*Info(ifile)%alpha_i(i)

          !! --------- DEBUG: --------- !!
          if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
            write (lun_db,*) "ibeg_dataL,  len: ", Info(ifile)%ibeg_dataL(i), len
            write (lun_db,*) "ibeg_Pij, jmax_i: ", Info(ifile)%ibeg_Pij(i),Info(ifile)%jmax_i(i)
          endif
          !! --------- DEBUG: --------- !!

          ! spin
          !if (nspin.EQ.1) then
            do jbeta_alpha = 1, len
              read (lun,*) Info(ifile)%data_Lold(Info(ifile)%ibeg_dataL(i)+jbeta_alpha-1, 1)
              if (flag_MDdebug .AND. iprint_MDdebug.GT.3) &
                 write (lun_db,'(f25.18)')               &
                 Info(ifile)%data_Lold(Info(ifile)%ibeg_dataL(i)+jbeta_alpha-1, 1)
            enddo
          !elseif (nspin.EQ.2) then
          !  do jbeta_alpha = 1, len
          !    read (lun,*) Info(ifile)%data_Lold(Info(ifile)%ibeg_dataL(i)+jbeta_alpha-1, 1), &
          !                 Info(ifile)%data_Lold(Info(ifile)%ibeg_dataL(i)+jbeta_alpha-1, 2)
          !    if (flag_MDdebug .AND. iprint_MDdebug.GT.3) &
          !       write (lun_db,'(3f25.18)')               &
          !         Info(ifile)%data_Lold(Info(ifile)%ibeg_dataL(i)+jbeta_alpha-1, 1), &
          !         Info(ifile)%data_Lold(Info(ifile)%ibeg_dataL(i)+jbeta_alpha-1, 2)
          !  enddo
          !endif

          !ORI if (i+1.LE.size) Info(ifile)%ibeg_dataL(i+1) = len + 1
          if (i+1.LE.size) then
            Info(ifile)%ibeg_Pij(i+1) = Info(ifile)%ibeg_Pij(i) + Info(ifile)%jmax_i(i)
            Info(ifile)%ibeg_dataL(i+1) = Info(ifile)%ibeg_dataL(i) + len
          endif
        enddo !(i, size)

        ! Close the file.
        call io_close (lun)
      enddo !(ifile, nfile)
    endif

    !! ---- DEBUG ---- !!
    if (flag_MDdebug .AND. iprint_MDdebug.GT.3) call io_close(lun_db)
    !! ---- DEBUG ---- !!

    return
  end subroutine grab_matrix2
  !!***

  ! -----------------------------------------------------------------------
  ! Subroutine dump_InfoGlobal
  ! -----------------------------------------------------------------------

  !!****f* io_module2/dump_InfoGlobal *
  !!
  !!  NAME
  !!   dump_InfoGlobal
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Writes out the global information
  !!    l.1 no. of atoms, no. of procs.
  !!    l.2 no. of partitions along a,b and c axes
  !!    l.3 cell size
  !!    l.4 relations between atoms (global label) and procs.
  !!    (optional) MD iteration
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/08/21
  !!  MODIFICATION
  !!
  !!  SOURCE
  !!
  subroutine dump_InfoGlobal(MDiter)

    ! Module usage
    use global_module, ONLY: ni_in_cell,numprocs,rcellx,rcelly,rcellz,id_glob
    use GenComms, ONLY: cq_abort
    use group_module, ONLY: parts

    ! passed variables
    integer,intent(in),optional :: MDiter

    ! local variables
    integer :: lun, stat, stat_alloc
    integer, allocatable :: glob_to_node_old_local(:)
    integer :: ind_part, id_node, ni, id_global


    ! Make a table showing atoms (glob) in nodes (CC) in the PREVIOUS job.
    allocate (glob_to_node_old_local(ni_in_cell), STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error in allocating glob_to_node_old_local: ', &
                                        ni_in_cell)
    glob_to_node_old_local = 0
    do ind_part = 1, parts%mx_gcell
      id_node = parts%i_cc2node(ind_part) ! CC labeling
      if (parts%nm_group(ind_part).NE.0) then
        do ni = 1, parts%nm_group(ind_part)
          id_global = id_glob(parts%icell_beg(ind_part)+ni-1)
          glob_to_node_old_local(id_global) = id_node
          if (glob_to_node_old_local(id_global).EQ.0) &
             call cq_abort('Error in glob_to_node_old_local.')
        enddo
      endif
    enddo !(ind_part, mx_gcell)

    ! Open InfoGlobal.dat and write data.
    call io_assign(lun)
    open (lun,file='InfoGlobal.dat',iostat=stat)
    if (stat.GT.0) call cq_abort('Fail in opening InfoGlobal.dat .')
    write (lun,*) ni_in_cell, numprocs
    write (lun,*) parts%ngcellx,parts%ngcelly,parts%ngcellz
    write (lun,*) rcellx,rcelly,rcellz
    write (lun,*) glob_to_node_old_local(1:ni_in_cell)
    if (present(MDiter)) write (lun,*) MDiter
    call io_close (lun)
    deallocate (glob_to_node_old_local, STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error deallocating glob_to_node_old_local :', &
                                       ni_in_cell)

    return
  end subroutine dump_InfoGlobal
  !!***

  ! -----------------------------------------------------------------------
  ! Subroutine grab_InfoGlobal
  ! -----------------------------------------------------------------------

  !!****f* io_module2/grab_InfoGlobal *
  !!
  !!  NAME
  !!   grab_InfoGlobal
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Reads the global data used for communication & reconstruction
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/08/21
  !!  MODIFICATION
  !!
  !!  SOURCE
  !!
  !subroutine grab_InfoGlobal(n_proc_old,glob_to_node_old)
  subroutine grab_InfoGlobal(n_proc_old,glob_to_node_old,MDiter)

    ! Module usage
    use GenComms, ONLY: inode,ionode,gcopy

    ! passed variables
    integer :: n_proc_old, glob_to_node_old(:)
    integer,optional :: MDiter

    ! local variables
    integer :: lun, stat
    integer :: n_atom,npcellx,npcelly,npcellz
    real(double) :: rx,ry,rz

    call io_assign(lun)
    open (lun,file='InfoGlobal.dat',status='old',iostat=stat)
    if (stat.GT.0) call cq_abort('Fail in opening InfoGlobal.dat .')
    read (lun,*) n_atom, n_proc_old
    read (lun,*) npcellx, npcelly, npcellz
    read (lun,*) rx, ry, rz
    read (lun,*) glob_to_node_old(1:)
    if (present(MDiter)) read (lun,*) MDiter
    call io_close(lun)

    return
  end subroutine grab_InfoGlobal
  !!***

  !!****f* io_module2/dump_idglob_old *
  !!
  !!  NAME
  !!   dump_id_glob_old
  !!  USAGE
  !!
  !!  PURPOSE
  !!
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
  ! -----------------------------------------------------------------------
  ! Subroutine dump_idglob_old
  ! -----------------------------------------------------------------------
  subroutine dump_idglob_old
    
    ! Module usage
    use global_module, ONLY: id_glob_old,id_glob_inv_old
    use input_module, ONLY: io_assign,io_close
    use GenComms, ONLY: cq_abort,inode

    ! local variables
    integer :: lun,stat

    ! Open section.
    call io_assign(lun)
    open (lun,file='InfoGlobal2.dat',iostat=stat)
    if (stat.NE.0) call cq_abort('Error in opening InfoGlobal2.dat:')

    write (lun,*) id_glob_old(1:)
    write (lun,*) id_glob_inv_old(1:)

    ! Close section.
    call io_close(lun)

    return
  end subroutine dump_idglob_old
  !!***

  !!****f* io_module2/grab_idglob_old *
  !!
  !!  NAME
  !!   grab_id_glob_old
  !!  USAGE
  !!
  !!  PURPOSE
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
  !!
  !!  SOURCE
  !!
  subroutine grab_idglob_old

    ! Module usage
    use global_module, ONLY: id_glob_old,id_glob_inv_old,ni_in_cell,io_lun
    use input_module, ONLY: io_assign,io_close
    use GenComms, ONLY: cq_abort,inode,ionode

    ! local variables
    integer :: lun,stat,stat_alloc
    
    ! Open section.
    call io_assign(lun)
    open (lun,file='InfoGlobal2.dat',status='old',iostat=stat)
    if (stat.NE.0) call cq_abort('Fail in opening InfoGlobal2.dat:')

    ! In the case where the arrays are not allocated.
    if (.NOT. allocated(id_glob_old)) then
      if (inode.EQ.ionode) write (io_lun,*) "WARNING ! Allocate call for id_glob_old."
      allocate (id_glob_old(ni_in_cell), STAT=stat_alloc)
      if (stat_alloc.NE.0) call cq_abort('Error allocating id_glob_old: ', ni_in_cell)
    endif
    if (.NOT. allocated(id_glob_inv_old)) then
      if (inode.EQ.ionode) write (io_lun,*) "WARNING ! Allocate call for id_glob_inv_old."
      allocate (id_glob_inv_old(ni_in_cell), STAT=stat_alloc)
      if (stat_alloc.NE.0) call cq_abort('Error allocating id_glob_inv_old: ', ni_in_cell)
    endif

    read (lun,*) id_glob_old(1:)
    read (lun,*) id_glob_inv_old(1:)

    ! Close section.
    call io_close(lun)

    return
  end subroutine grab_idglob_old
  !!***

  !!****f* io_module2/grab_InfoGlobal *
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

end module io_module2
