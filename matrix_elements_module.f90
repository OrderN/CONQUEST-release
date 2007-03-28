! $Id: matrix_elements_module.f90,v 1.1.1.1.2.1 2006/03/31 12:26:00 drb Exp $
! -----------------------------------------------------------
! Module matrix_elements_module
! -----------------------------------------------------------
! Code area 2: matrices
! -----------------------------------------------------------

!!****h* Conquest/matrix_elements_module
!!  NAME
!!   matrix_elements_module
!!  PURPOSE
!!   Creates all necessary indexing for matrices, haloes etc
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!   07/12/99 by D.R.Bowler
!!    Check put on mx_abs in get_naba
!!   19/02/00 by DRB
!!    Implementing transposes
!!   01/03/00 by DRB
!!    Finishing transposes
!!   12/04/00 by DRB 
!!    Moved itran_addr into trans_remote type
!!   20/04/00 by DRB 
!!    Removed the pass of pairs to send_pair_info  
!!   25/04/00 by DRB
!!    Now uses the generic group, primary and cover_set types
!!   08/06/2001 dave
!!    Added ROBODoc header and GenComms for my_barrier and cq_abort
!!***
module matrix_elements_module

  implicit none

contains

!!****f* matrix_elements_module/matrix_ini
!!
!!  NAME 
!!   matrix_ini
!!  USAGE
!! 
!!  PURPOSE
!!   Creates indexing for matrix
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!   08/06/2001 dave
!!    Added ROBODoc header and GenComms for my_barrier
!!  SOURCE
!!
!  subroutine matrix_ini(parts,prim,gcs,amat,aind,rcut,myid,&
!       ahalo,atrans,atrans_rem,apairs,apairind)
  subroutine matrix_ini(parts,prim,gcs,amat,aind,rcut,myid, ahalo,atrans)

    ! Module usage
    use datatypes
    use global_module
    use basic_types
    use matrix_module
    use trans_module
    use GenComms, ONLY: my_barrier, gmax

    implicit none

    ! Passed variables     
    ! First, the small group variables
    type(group_set) :: parts
    type(primary_set) :: prim
    type(cover_set) :: gcs
    ! Now general integers
    integer :: myid,nd1,nd2
    real(double) :: rcut
    ! These variables are generic matrix variables
    type(matrix), dimension (:) :: amat  
    type(matrix_trans),OPTIONAL :: atrans
    type(matrix_halo),OPTIONAL :: ahalo
    !type(trans_remote),OPTIONAL :: atrans_rem
    !type(pair_data),OPTIONAL    :: apairs(:)
    !integer(integ), pointer, dimension(:), OPTIONAL :: apairind
    integer(integ), pointer, dimension(:) :: aind

    ! Local variables
    integer :: nnd,nr,posn
    integer :: irc,ierr

    nnd = myid+1
    if(prim%n_prim.gt.0) then ! If we actually HAVE a primary set
       call get_naba_max(prim,gcs,amat,rcut)
       ! The matrix maxima MUST be global across processors !
       call gmax(amat(1)%mx_nab)
       amat(2:prim%groups_on_node)%mx_nab = amat(1)%mx_nab
       call gmax(amat(1)%mx_abs)
       amat(2:prim%groups_on_node)%mx_abs = amat(1)%mx_abs
       !write(*,*) 'Matrix maxima: ',amat(1)%mx_nab, amat(1)%mx_abs
       !write(*,*) 'Index length: ',parts%mx_ngonn*(3*parts%mx_mem_grp+5*parts%mx_mem_grp*amat(1)%mx_abs)
       allocate(aind(parts%mx_ngonn*(3*parts%mx_mem_grp+5*parts%mx_mem_grp*amat(1)%mx_abs)),STAT=ierr)
       if(ierr/=0) call cq_abort("Error allocating matrix index: ", &
            parts%mx_ngonn*(3*parts%mx_mem_grp+5*parts%mx_mem_grp*amat(1)%mx_abs),ierr)
       aind = 0
       call set_matrix_pointers(prim%nm_nodgroup,amat, aind, prim%groups_on_node,parts%mx_mem_grp)
       call allocate_matrix(amat, prim%groups_on_node,parts%mx_mem_grp)
       call get_naba(prim,gcs,amat,rcut)
       ! Now that we know neighbour numbers, sort out the pointers
       call set_matrix_pointers2(prim%nm_nodgroup,amat,aind, &
            prim%groups_on_node,parts%mx_mem_grp)
       ! * IMPORTANT * ! This must come next: it builds ndimj
       call make_npxyz(nnd,amat,prim,gcs,parts%mx_mem_grp)
       if(PRESENT(ahalo)) then ! make index arrays for halo 
          call make_halo_max(prim,gcs,amat,ahalo)
          call gmax(ahalo%mx_part)
          call gmax(ahalo%mx_halo)
          !write(*,*) 'Halo maxima: ',ahalo%mx_part, ahalo%mx_halo
          call allocate_halo(ahalo,prim%mx_iprim,gcs%mx_mcover,gcs%mx_gcover)
          call make_halo(prim,gcs,amat,ahalo)
       endif ! halo
       if(PRESENT(atrans)) then ! make index arrays for retro storage 
          atrans%mx_halo = ahalo%mx_halo
          atrans%mx_nab = amat(1)%mx_nab
          call allocate_trans(atrans,prim%mx_iprim)
          call get_retr(prim,ahalo,atrans)
       endif ! trans
       ! make index arrays for global transposes
       call my_barrier
    endif ! n_prim.gt.0
    return
  end subroutine matrix_ini
!!***

!!****f* matrix_elements_module/trans_ini *
!!
!!  NAME 
!!   trans_ini
!!  USAGE
!! 
!!  PURPOSE
!!   Initialises transposes
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2006/09/07 (though just relocating earlier code)
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine trans_ini(parts,prim,gcs,amat,myid,&
       ahalo,ahalo_rem,atrans,atrans_rem,apairs,apairind)

    ! Module usage
    use datatypes
    use global_module
    use basic_types
    use matrix_module
    use trans_module
    use GenComms, ONLY: my_barrier, gmax
    use maxima_module, ONLY: maxnabaprocs

    implicit none

    ! Passed variables     
    ! First, the small group variables
    type(group_set) :: parts
    type(primary_set) :: prim
    type(cover_set) :: gcs
    ! Now general integers
    integer :: myid
    ! These variables are generic matrix variables
    type(matrix), dimension (:) :: amat  
    type(matrix_trans) :: atrans
    type(matrix_halo) :: ahalo
    type(matrix_halo) :: ahalo_rem
    type(trans_remote) :: atrans_rem
    type(pair_data)    :: apairs(:)
    integer(integ), pointer, dimension(:) :: apairind

    ! Local variables
    integer :: nnd,nr,posn
    integer :: irc,ierr


    nnd = myid+1
    if(prim%n_prim.gt.0) then ! If we actually HAVE a primary set
       allocate(apairind(4*prim%mx_iprim*amat(1)%mx_nab), STAT=ierr)
       if(ierr/=0) call cq_abort("Error allocaing pair index: ",4*prim%mx_iprim*amat(1)%mx_nab,ierr)
       apairind = 0
       call allocate_trans_rem(atrans_rem,maxnabaprocs, prim%mx_iprim*amat(1)%mx_nab)
       call trans_halo_info(numprocs,nnd,atrans_rem,ahalo,atrans,parts)
       call set_trans_pointers(apairind,atrans_rem,apairs)
       call a_and_b(nnd,atrans,ahalo,atrans_rem,apairs,prim,gcs)
       call my_barrier
       call send_pair_info(nnd,apairind,atrans_rem)
       call my_barrier
       call index_transpose(prim%mx_iprim,nnd,atrans_rem, &
            ahalo,ahalo_rem,apairind,apairs,parts,prim,gcs)
    endif 
  end subroutine trans_ini
!!***

!!****f* matrix_elements_module/get_naba *
!!
!!  NAME 
!!   get_naba
!!  USAGE
!! 
!!  PURPOSE
!!   Creates neighbour lists for a matrix range rcut
!!  INPUTS
!!   type(primary_set) :: prim   
!!   type(cover_set) :: gcs
!!   type(matrix), dimension(:) :: amat
!!   real(double) :: rcut
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!   08/06/2001 dave
!!    Added ROBODoc header, GenComms and
!!    removed MPI_Abort for cq_abort
!!  SOURCE
!!
  subroutine get_naba(prim,gcs,amat,rcut)

    ! Module usage
    use datatypes
    use basic_types
    use matrix_module
    use GenComms, ONLY: cq_abort, inode
    use global_module, ONLY: id_glob, species_glob, sf, nlpf, paof
    use group_module, ONLY: parts
    use species_module, ONLY: nsf_species, nlpf_species, npao_species

    implicit none

    ! Passed variables
    type(primary_set) :: prim
    type(cover_set) :: gcs
    real(double) :: rcut
    type(matrix) :: amat(:)

    ! Local variables
    integer :: irc,ierr,inp, cumu_ndims, neigh_spec
    integer :: nn,j,np,ni,ist
    real(double) :: rcutsq,dx,dy,dz
    real(double), parameter :: tol=1.0e-8_double

    ! Check that prim and gcs are correctly set up
    if((.NOT.ASSOCIATED(gcs%xcover)).OR. &
         (.NOT.ASSOCIATED(prim%xprim))) then
       call cq_abort('get_naba: gcs or prim without members !')
    endif
    rcutsq=rcut*rcut
    ! loop over all atom pairs (atoms in primary set, max. cover set) -
    inp=1  ! Indexes primary atoms
    amat(1)%offset=0
    amat(1)%nd_offset=0
    do nn=1,prim%groups_on_node ! Partitions in primary set
       amat(nn)%part_nabs = 0    ! Cumulative neighbours of partition
       amat(nn)%part_nd_nabs = 0    ! Cumulative neighbours of partition
       if(prim%nm_nodgroup(nn).gt.0) then  ! Are there atoms ?
          amat(nn)%i_acc(1)=1
          amat(nn)%i_nd_acc(1)=1
          amat(nn)%n_atoms = prim%nm_nodgroup(nn) ! Redundant, but useful
          !write(*,*) 'Starting group with atoms: ',nn,prim%nm_nodgroup(nn)
          do j=1,prim%nm_nodgroup(nn)  ! Loop over atoms in partition
             amat(nn)%n_nab(j)=0
             select case(amat(nn)%sf1_type)
             case(sf)
                amat(nn)%ndimi(j) = nsf_species(prim%species(inp))
             case(nlpf)
                amat(nn)%ndimi(j) = nlpf_species(prim%species(inp))
             case(paof)
                amat(nn)%ndimi(j) = npao_species(prim%species(inp))
             end select
             cumu_ndims = 0
             do np=1,gcs%ng_cover  ! Loop over partitions in GCS
                if(gcs%n_ing_cover(np).gt.0) then  ! Are there atoms ?
                   if(gcs%icover_ibeg(np)+gcs%n_ing_cover(np)-1.gt.&
                        gcs%mx_mcover) then
                      call cq_abort('get_naba: overran gcs mx_mcover: ', &
                           gcs%icover_ibeg(np)+gcs%n_ing_cover(np)-1, &
                           gcs%mx_mcover)
                   endif
                   do ni=1,gcs%n_ing_cover(np)
                      dx=gcs%xcover(gcs%icover_ibeg(np)+ni-1)-prim%xprim(inp)
                      dy=gcs%ycover(gcs%icover_ibeg(np)+ni-1)-prim%yprim(inp)
                      dz=gcs%zcover(gcs%icover_ibeg(np)+ni-1)-prim%zprim(inp)
                      if(dx*dx+dy*dy+dz*dz<rcutsq-tol) then
                         amat(nn)%n_nab(j)=amat(nn)%n_nab(j)+1
                         !write(*,*) 'Neighbour: ',j,amat(nn)%n_nab(j)
                         if(amat(nn)%n_nab(j).gt.amat(nn)%mx_abs) then
                            call cq_abort('get_naba: n_nab>mx_nab: ',amat(nn)%n_nab(j), &
                                 amat(nn)%mx_abs)
                         endif
                         ist = amat(nn)%i_acc(j)+amat(nn)%n_nab(j)-1
                         amat(nn)%i_part(ist)=np
                         amat(nn)%i_seq(ist)=ni
                         amat(nn)%radius(ist)=sqrt(dx*dx+dy*dy+dz*dz)
                         if(gcs%lab_cover(np)==gcs%iprim_group(inp).AND.ni==j) then
                            amat(nn)%onsite(j)=amat(nn)%nd_offset+amat(nn)%i_nd_acc(j)+cumu_ndims
                         endif
                         neigh_spec = species_glob( id_glob( parts%icell_beg(gcs%lab_cell(np)) +ni-1 ))
                         select case(amat(nn)%sf2_type)
                         case(sf)
                            cumu_ndims = cumu_ndims + amat(nn)%ndimi(j)*nsf_species(neigh_spec)
                         case(nlpf)
                            cumu_ndims = cumu_ndims + amat(nn)%ndimi(j)*nlpf_species(neigh_spec)
                         case(paof)
                            cumu_ndims = cumu_ndims + amat(nn)%ndimi(j)*npao_species(neigh_spec)
                         end select
                      endif
                   enddo ! End n_inp_cover
                endif
             enddo ! End np_cover
             !write(*,*) 'Finishing prim atom: ',inode,inp,cumu_ndims,j,prim%nm_nodgroup(nn)
             if(j.lt.prim%nm_nodgroup(nn)) then
                amat(nn)%i_acc(j+1)=amat(nn)%i_acc(j)+amat(nn)%n_nab(j)
                amat(nn)%i_nd_acc(j+1)=amat(nn)%i_nd_acc(j)+cumu_ndims
             endif
             amat(nn)%part_nabs = amat(nn)%part_nabs+amat(nn)%n_nab(j)
             amat(nn)%part_nd_nabs = amat(nn)%part_nd_nabs+cumu_ndims
             inp=inp+1  ! Indexes primary-set atoms
          enddo ! End prim%nm_nodgroup
          if(nn.lt.prim%groups_on_node) then 
             amat(nn+1)%offset=amat(nn)%offset+amat(nn)%part_nabs
             amat(nn+1)%nd_offset=amat(nn)%nd_offset+amat(nn)%part_nd_nabs
          endif
       else
          amat(nn)%n_atoms = 0
          if(nn.lt.prim%groups_on_node) then 
             amat(nn+1)%offset=amat(nn)%offset+amat(nn)%part_nabs
             amat(nn+1)%nd_offset=amat(nn)%nd_offset+amat(nn)%part_nd_nabs
          endif
          amat(nn)%i_acc(1) = 1
          amat(nn)%i_nd_acc(1) = 1
       endif ! End if(prim%nm_nodgroup>0)
    enddo ! End part_on_node
    ! Now create length of the matrix
    amat%length=0
    do nn=1,prim%groups_on_node
       amat(1)%length=amat(1)%length+amat(nn)%part_nd_nabs
    enddo
    amat(2:prim%groups_on_node)%length = amat(1)%length
    return
  end subroutine get_naba
!!***

!!****f* matrix_elements_module/get_naba_max *
!!
!!  NAME 
!!   get_naba_max
!!  USAGE
!! 
!!  PURPOSE
!!   Creates maxima for neighbour lists
!!  INPUTS
!!   type(primary_set) :: prim   
!!   type(cover_set) :: gcs
!!   type(matrix), dimension(:) :: amat
!!   real(double) :: rcut
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2006/08/08
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine get_naba_max(prim,gcs,amat,rcut)

    ! Module usage
    use datatypes
    use basic_types
    use matrix_module
    use GenComms, ONLY: cq_abort, inode
    use global_module, ONLY: id_glob, species_glob, sf, nlpf, paof
    use group_module, ONLY: parts
    use species_module, ONLY: nsf_species, nlpf_species, npao_species

    implicit none

    ! Passed variables
    type(primary_set) :: prim
    type(cover_set) :: gcs
    real(double) :: rcut
    type(matrix) :: amat(:)

    ! Local variables
    integer :: irc,ierr,inp, cumu_ndims, neigh_spec
    integer :: nn,j,np,ni,ist
    integer :: tot_nabs, nabs_of_atom, mx_abs_nabs
    real(double) :: rcutsq,dx,dy,dz
    real(double), parameter :: tol=1.0e-8_double

    ! Check that prim and gcs are correctly set up
    if((.NOT.ASSOCIATED(gcs%xcover)).OR. &
         (.NOT.ASSOCIATED(prim%xprim))) then
       call cq_abort('get_naba: gcs or prim without members !')
    endif
    rcutsq=rcut*rcut
    ! loop over all atom pairs (atoms in primary set, max. cover set) -
    inp=1  ! Indexes primary atoms
    tot_nabs = 0
    nabs_of_atom = 0
    mx_abs_nabs = 0
    do nn=1,prim%groups_on_node ! Partitions in primary set
       if(prim%nm_nodgroup(nn).gt.0) then  ! Are there atoms ?
          do j=1,prim%nm_nodgroup(nn)  ! Loop over atoms in partition
             nabs_of_atom = 0
             do np=1,gcs%ng_cover  ! Loop over partitions in GCS
                if(gcs%n_ing_cover(np).gt.0) then  ! Are there atoms ?
                   if(gcs%icover_ibeg(np)+gcs%n_ing_cover(np)-1.gt.&
                        gcs%mx_mcover) then
                      call cq_abort('get_naba: overran gcs mx_mcover: ', &
                           gcs%icover_ibeg(np)+gcs%n_ing_cover(np)-1, &
                           gcs%mx_mcover)
                   endif
                   do ni=1,gcs%n_ing_cover(np)
                      dx=gcs%xcover(gcs%icover_ibeg(np)+ni-1)-prim%xprim(inp)
                      dy=gcs%ycover(gcs%icover_ibeg(np)+ni-1)-prim%yprim(inp)
                      dz=gcs%zcover(gcs%icover_ibeg(np)+ni-1)-prim%zprim(inp)
                      if(dx*dx+dy*dy+dz*dz<rcutsq-tol) then
                         nabs_of_atom = nabs_of_atom + 1
                         tot_nabs = tot_nabs + 1
                      endif
                   enddo ! End n_inp_cover
                endif
             enddo ! End np_cover
             if(nabs_of_atom>mx_abs_nabs) mx_abs_nabs = nabs_of_atom
             inp=inp+1  ! Indexes primary-set atoms
          enddo ! End prim%nm_nodgroup
       endif ! End if(prim%nm_nodgroup>0)
    enddo ! End part_on_node
    amat(1:prim%groups_on_node)%mx_nab = (tot_nabs+prim%mx_iprim-1)/prim%mx_iprim 
    amat(1:prim%groups_on_node)%mx_abs = mx_abs_nabs
    return
  end subroutine get_naba_max
!!***

!!****f* matrix_elements_module/get_retr *
!!
!!  NAME 
!!   get_retr
!!  USAGE
!! 
!!  PURPOSE
!!   Constructs indexing for "retro" storage - in other
!!   words for local transposes (required for matrix mults)
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   08/06/2001 dave
!!    Added ROBODoc header and changed MPI_abort to cq_abort
!!  SOURCE
!!
  subroutine get_retr(prim,ahalo,at)

    ! Module usage
    use datatypes
    use matrix_module
    use GenComms, ONLY: cq_abort

    implicit none

    type(primary_set) :: prim
    type(matrix_halo)::ahalo
    type(matrix_trans)::at
    integer :: irc,ierr,ni,i,cumu_ndims

    if((ahalo%ni_in_halo.le.0).or.(ahalo%ni_in_halo.gt.ahalo%mx_halo)) then
       call cq_abort('get_retr: ni_in_halo out of range')
    endif
    ! --- loop over all atom pairs (halo, primary set) --------------------
    at%i_beg(1)=1
    at%i_nd_beg(1)=1
    do ni=1,ahalo%ni_in_halo
       at%n_hnab(ni)=0
       cumu_ndims = 0
       do i=1,prim%n_prim
          if(ahalo%i_h2d((i-1)*ahalo%ni_in_halo+ni).ne.0) then
             at%n_hnab(ni)=at%n_hnab(ni)+1
             cumu_ndims = cumu_ndims + ahalo%ndimi(i)*ahalo%ndimj(ni)
             if(at%i_beg(ni)+at%n_hnab(ni)-1.gt.prim%mx_iprim*at%mx_nab) then
                call cq_abort('get_retr: i_prim arg. out of range')
             endif
             at%i_prim(at%i_beg(ni)+at%n_hnab(ni)-1)=i
          endif
       enddo
       if(ni.lt.ahalo%ni_in_halo) then
          at%i_beg(ni+1)=at%i_beg(ni)+at%n_hnab(ni)
          at%i_nd_beg(ni+1)=at%i_nd_beg(ni)+cumu_ndims
       endif
    enddo
    ! ---------------------------------------------------------------------
    return
  end subroutine get_retr
!!***

!!****f* matrix_elements_module/make_halo *
!!
!!  NAME 
!!   make_halo
!!  USAGE
!! 
!!  PURPOSE
!!   Constructs the "halo" - in other words, lists of
!!   all atoms and all partitions within range of ANY atom
!!   in the primary set.  So the halo atoms are a superset of
!!   the neighbours of atoms in the primary set, while the
!!   halo partitions are simply all partitions containing at
!!   least one halo atom.  Halo processors are those processors
!!   responsible for at least ONE halo partition.  These are 
!!   all found by searching over the cover set.
!!  INPUTS
!!   type(primary_set) :: prim
!!   type(cover_set) :: gcs
!!   type(matrix) :: amat(:)
!!   type(matrix_halo) :: ahalo
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   08/06/2001 dave
!!    Added ROBODoc header and cq_abort from GenComms
!!  SOURCE
!!
  subroutine make_halo(prim,gcs,amat,ahalo)

    use matrix_module
    use basic_types
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: id_glob, species_glob, sf, nlpf, paof
    use group_module, ONLY: parts
    use species_module, ONLY: nsf_species, nlpf_species, npao_species

    implicit none

    ! Passed variables
    type(primary_set) :: prim
    type(cover_set) :: gcs
    type(matrix)::amat(:)
    type(matrix_halo)::ahalo

    ! Local variables
    integer :: irc
    integer :: ierr,np,ni,i,nb,i_in_part,nh,j,ist,isu,nd1
    integer(integ) :: np_in_halo,ni_in_halo

    integer :: stat,pr, neigh_spec

    ! initialise flag for every part and atom in covering set
    do np=1,gcs%ng_cover
       if(gcs%n_ing_cover(np).gt.0) then
          if(gcs%icover_ibeg(np)+gcs%n_ing_cover(np)-1.gt.gcs%mx_mcover) then
             call cq_abort('make_halo: ahalo%i_halo arg. out of range')
          endif
          do ni=1,gcs%n_ing_cover(np)
             ahalo%i_halo(gcs%icover_ibeg(np)+ni-1)=0
          enddo
       endif
    enddo
    ! find atoms in covering set that are neighbours of primary set
    pr = 1
    do i=1,prim%groups_on_node
       if(prim%nm_nodgroup(i)>0) then
          do j=1,prim%nm_nodgroup(i)
             ahalo%ndimi(pr) = amat(i)%ndimi(j)
             do nb=1,amat(i)%n_nab(j)
                np=amat(i)%i_part(amat(i)%i_acc(j)+nb-1)
                ni=amat(i)%i_seq(amat(i)%i_acc(j)+nb-1)
                ahalo%i_halo(gcs%icover_ibeg(np)+ni-1)=1
             enddo
             pr = pr+1
          enddo
       end if ! (prim%nm_nodgroup(i)>0)
    enddo
    ahalo%np_in_halo=0
    ahalo%ni_in_halo=0
    do np=1,gcs%ng_cover
       i_in_part=0
       ahalo%i_hbeg(gcs%lab_cover(np))=gcs%icover_ibeg(np)
       if(gcs%n_ing_cover(np).gt.0) then
          do ni=1,gcs%n_ing_cover(np)
             if(ahalo%i_halo(gcs%icover_ibeg(np)+ni-1).eq.1) then
                i_in_part=i_in_part+1
                ahalo%ni_in_halo=ahalo%ni_in_halo+1
                if(ahalo%ni_in_halo.gt.ahalo%mx_halo) then
                   call cq_abort('make_halo: too many atoms in halo',&
                        ahalo%ni_in_halo,ahalo%mx_halo)
                endif
                ahalo%j_seq(ahalo%ni_in_halo)=ni
                neigh_spec = species_glob( id_glob(parts%icell_beg(gcs%lab_cell(np))+ni-1) )
                select case(amat(1)%sf2_type)
                case(sf)
                   ahalo%ndimj(ahalo%ni_in_halo) = nsf_species(neigh_spec)
                case(nlpf)
                   ahalo%ndimj(ahalo%ni_in_halo) = nlpf_species(neigh_spec)
                case(paof)
                   ahalo%ndimj(ahalo%ni_in_halo) = npao_species(neigh_spec)
                end select
                ahalo%i_halo(gcs%icover_ibeg(np)+ni-1)=ahalo%ni_in_halo
             endif
          enddo
       endif
       if(i_in_part.gt.0) then
          ahalo%np_in_halo=ahalo%np_in_halo+1
          if(ahalo%np_in_halo.gt.ahalo%mx_part) then
             call cq_abort('make_halo: too many partitions in halo',&
                  ahalo%np_in_halo)
          endif
          if(ahalo%np_in_halo.eq.1) then
             ahalo%j_beg(ahalo%np_in_halo)=1
          else
             ahalo%j_beg(ahalo%np_in_halo)=ahalo%j_beg(ahalo%np_in_halo-1)&
                  +ahalo%nh_part(ahalo%np_in_halo-1)
          endif
          ahalo%nh_part(ahalo%np_in_halo)=i_in_part
          ahalo%lab_hcell(ahalo%np_in_halo)=gcs%lab_cell(np)
          ahalo%lab_hcover(ahalo%np_in_halo)=gcs%lab_cover(np)
       endif
    enddo
    ! --- construct neighbour-halo transcription table --------------------
    if(ahalo%ni_in_halo.le.0) then
       call cq_abort('make_halo: no. of atoms in halo must be .ge. 1')
    endif
    do i=1,prim%n_prim
       do nh=1,ahalo%ni_in_halo
          ahalo%i_h2d((i-1)*ahalo%ni_in_halo+nh)=0
       enddo
    enddo
    pr = 1
    do i=1,prim%groups_on_node
       if(prim%nm_nodgroup(i)>0) then
          do j=1,prim%nm_nodgroup(i)
             ist = amat(i)%nd_offset+amat(i)%i_nd_acc(j)
             nd1 = amat(i)%ndimi(j)
             do nb=1,amat(i)%n_nab(j)
                isu = amat(i)%i_acc(j)+(nb-1)
                np=amat(i)%i_part(isu)
                ni=amat(i)%i_seq(isu)
                ahalo%i_h2d((pr-1)*ahalo%ni_in_halo +ahalo%i_halo(gcs%icover_ibeg(np)+ni-1))= ist
                ist = ist+nd1*amat(i)%ndimj(isu)
             enddo
             pr = pr+1
          enddo
       end if ! (prim%nm_nodgroup(i)>0)
    enddo
    ! ---------------------------------------------------------------------
    return
  end subroutine make_halo
!!***

!!****f* matrix_elements_module/make_halo_max *
!!
!!  NAME 
!!   make_halo_max
!!  USAGE
!! 
!!  PURPOSE
!!   Find halo maxima
!!  INPUTS
!!   type(primary_set) :: prim
!!   type(cover_set) :: gcs
!!   type(matrix) :: amat(:)
!!   type(matrix_halo) :: ahalo
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2006/08/08
!!  MODIFICATION HISTORY
!! 
!!  SOURCE
!!
  subroutine make_halo_max(prim,gcs,amat,ahalo)

    use matrix_module
    use basic_types
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: id_glob, species_glob, sf, nlpf, paof
    use group_module, ONLY: parts
    use species_module, ONLY: nsf_species, nlpf_species, npao_species

    implicit none

    ! Passed variables
    type(primary_set) :: prim
    type(cover_set) :: gcs
    type(matrix)::amat(:)
    type(matrix_halo)::ahalo

    ! Local variables
    integer :: irc, max_atoms_halo, max_part_halo
    integer :: ierr,np,ni,i,nb,i_in_part,nh,j,ist,isu,nd1
    integer :: np_in_halo,ni_in_halo

    integer :: stat,pr, neigh_spec

    allocate(ahalo%i_halo(gcs%mx_mcover),STAT=stat)
    if(stat/=0) write(*,*) 'Error allocating ihalo !'
    ! initialise flag for every part and atom in covering set
    do np=1,gcs%ng_cover
       if(gcs%n_ing_cover(np).gt.0) then
          if(gcs%icover_ibeg(np)+gcs%n_ing_cover(np)-1.gt.gcs%mx_mcover) then
             call cq_abort('make_halo: ahalo%i_halo arg. out of range')
          endif
          do ni=1,gcs%n_ing_cover(np)
             ahalo%i_halo(gcs%icover_ibeg(np)+ni-1)=0
          enddo
       endif
    enddo
    ! find atoms in covering set that are neighbours of primary set
    do i=1,prim%groups_on_node
       if(prim%nm_nodgroup(i)>0) then
          do j=1,prim%nm_nodgroup(i)
             do nb=1,amat(i)%n_nab(j)
                np=amat(i)%i_part(amat(i)%i_acc(j)+nb-1)
                ni=amat(i)%i_seq(amat(i)%i_acc(j)+nb-1)
                ahalo%i_halo(gcs%icover_ibeg(np)+ni-1)=1
             enddo
          enddo
       end if ! (prim%nm_nodgroup(i)>0)
    enddo
    max_atoms_halo = 0
    max_part_halo = 0
    do np=1,gcs%ng_cover
       i_in_part=0
       if(gcs%n_ing_cover(np).gt.0) then
          do ni=1,gcs%n_ing_cover(np)
             if(ahalo%i_halo(gcs%icover_ibeg(np)+ni-1).eq.1) then
                i_in_part=i_in_part+1
                max_atoms_halo = max_atoms_halo + 1
             endif
          enddo
       endif
       if(i_in_part>0) max_part_halo = max_part_halo + 1
    enddo
    ahalo%mx_part = max_part_halo
    ahalo%mx_halo = max_atoms_halo
    ! Dellocate memory
    deallocate(ahalo%i_halo,STAT=stat)
    if(stat/=0) write(*,*) 'Error allocating ihalo !'    
    return
  end subroutine make_halo_max
!!***

!!****f* matrix_elements_module/make_npxyz *
!!
!!  NAME 
!!   make_npxyz
!!  USAGE
!! 
!!  PURPOSE
!!   Makes the xyz parameters that tell processors
!!   how to convert from another processor's cover set
!!   partition to a local cover set partition
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   08/06/2001 dave
!!    Added ROBODoc header and GenComms for cq_abort
!!  SOURCE
!!
  subroutine make_npxyz(nnd,amat,prim,gcs,mx_part)

    ! Module usage
    use matrix_module
    use basic_types
    use GenComms, ONLY: cq_abort
    use group_module, ONLY: parts
    use species_module, ONLY: nsf_species, nlpf_species, npao_species
    use global_module, ONLY: id_glob, species_glob, sf, nlpf, paof

    implicit none

    ! Passed variables
    integer(integ) :: nnd,mx_part
    type(matrix), dimension(:) :: amat
    type(primary_set) :: prim
    type(cover_set) :: gcs

    ! Local variables
    integer(integ) :: np,i,ni,nb,ind_cover,ibpartp,ncoveryz
    integer(integ) :: nx,ny,nz,npx,npy,npz,ist,isu,neigh_spec

    ncoveryz = gcs%ncovery*gcs%ncoverz
    do np=1,prim%groups_on_node
       if(prim%nm_nodgroup(np).gt.0) then ! Are there atoms ?
          do i=1,prim%nm_nodgroup(np)
             ni=prim%nm_nodbeg(np)+i-1 ! Start of partn in list of primary atoms
             if(amat(np)%n_nab(i).gt.0) then 
                do nb=1,amat(np)%n_nab(i) ! Loop over neighbours
                   isu = amat(np)%i_acc(i)+nb-1
                   if(isu>mx_part*amat(np)%mx_abs) then
                      call cq_abort('make_npxyz: neighbour error', &
                           isu, mx_part*amat(np)%mx_abs)
                   endif
                   ind_cover=amat(np)%i_part(isu)
                   ibpartp=gcs%lab_cover(ind_cover)
                   nx=1+(ibpartp-1)/(ncoveryz)
                   ny=1+(ibpartp-1-(nx-1)*ncoveryz)/&
                        gcs%ncoverz
                   nz=ibpartp-(nx-1)*ncoveryz-&
                        (ny-1)*gcs%ncoverz
                   npx=nx-1-gcs%nspanlx-prim%idisp_primx(np)
                   npy=ny-1-gcs%nspanly-prim%idisp_primy(np)
                   npz=nz-1-gcs%nspanlz-prim%idisp_primz(np)
                   ! Store the xyz triplet locating this partition
                   if(3*(isu-1)+3>3*amat(np)%part_nabs) then
                      call cq_abort('make_mat: triplet error ',isu,amat(np)%part_nabs)
                   endif
                   neigh_spec = species_glob( id_glob( parts%icell_beg(gcs%lab_cell(ind_cover)) +amat(np)%i_seq(isu)-1 ))
                   select case(amat(np)%sf2_type)
                   case(sf)
                      amat(np)%ndimj(isu) = nsf_species(neigh_spec)
                   case(nlpf)
                      amat(np)%ndimj(isu) = nlpf_species(neigh_spec)
                   case(paof)
                      amat(np)%ndimj(isu) = npao_species(neigh_spec)
                   end select
                   amat(np)%npxyz(3*(isu-1)+1) = npx
                   amat(np)%npxyz(3*(isu-1)+2) = npy
                   amat(np)%npxyz(3*(isu-1)+3) = npz
                enddo
             endif
          enddo
       endif
    enddo
  end subroutine make_npxyz
!!***

!!****f* matrix_elements_module/trans_halo_info *
!!
!!  NAME 
!!   trans_halo_info
!!  USAGE
!! 
!!  PURPOSE
!!   Makes info about nodes in halo, numbers
!!   of halo neighbour pairs on each of these nodes, etc
!!   needed for communications in performing matrix transpose.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!! ?
!!  MODIFICATION HISTORY
!!   08/06/2001 dave
!!    Added ROBODoc header and cq_abort from GenComms
!!  SOURCE
!!
  subroutine trans_halo_info(nnode,mynode,adata,halo,trans,parts)

    use matrix_module
    use basic_types
    use trans_module
    use maxima_module, ONLY: maxnabaprocs
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    type(group_set) :: parts
    type(matrix_trans) :: trans
    type(trans_remote):: adata
    type(matrix_halo) :: halo

    integer(integ) :: nnode,mynode

    ! Local variables
    integer(integ) :: np,ind_part,nnd_rem,na,ia,nr,n

    ! Check on halo
    if(halo%np_in_halo.lt.1) then
       call cq_abort('trans_halo_info: no partitions !')
    endif
    ! Loop over partitions in halo
    adata%n_rem_node=0
    np=1
    ind_part=halo%lab_hcell(np)
    nnd_rem=parts%i_cc2node(ind_part)
    if(nnd_rem.ne.mynode) then
       call cq_abort('trans_halo_info: first halo node must be myself')
    endif
    adata%n_rem_node=adata%n_rem_node+1
    adata%list_rem_node(adata%n_rem_node)=nnd_rem
    adata%nhp_for_node(adata%n_rem_node)=1
    if(halo%np_in_halo.gt.1) then
       do np=2,halo%np_in_halo
          !write(*,*) 'Halo part: ',np,mynode
          ind_part=halo%lab_hcell(np)
          nnd_rem=parts%i_cc2node(ind_part)
          if(nnd_rem.ne.adata%list_rem_node(adata%n_rem_node)) then
             if(nnd_rem<adata%list_rem_node(adata%n_rem_node).AND. &
                  nnd_rem>mynode) then
                write(*,923) 
923             format(//'Possible error halo_com: halo nodes out of order')
             endif
             adata%n_rem_node=adata%n_rem_node+1
             if(adata%n_rem_node.gt.nnode) then
                call cq_abort('trans_halo_info: too many nodes',&
                     adata%n_rem_node)
             endif
             if(adata%n_rem_node>maxnabaprocs+1) then
                call cq_abort('trans_halo_info: too many neighbour nodes', adata%n_rem_node,maxnabaprocs+1)
             endif
             adata%list_rem_node(adata%n_rem_node)=nnd_rem
             adata%nhp_for_node(adata%n_rem_node)=1
          else
             adata%nhp_for_node(adata%n_rem_node)= &
                  adata%nhp_for_node(adata%n_rem_node)+1
          endif
       enddo
    endif
    ! determine no. of pairs on each halo-node, addresses etc. 
    na=0
    ia=0
    adata%i_pair_addr(1)=1
    do nr=1,adata%n_rem_node
       adata%n_pair(nr)=0
       do np=1,adata%nhp_for_node(nr)
          na=na+1
          do n=1,halo%nh_part(na)
             ia=ia+1
             adata%n_pair(nr)=adata%n_pair(nr)+trans%n_hnab(ia)
          enddo
       enddo ! End of np=1,adata%nhp_for_node
       if(nr.lt.adata%n_rem_node) then
          adata%i_pair_addr(nr+1)=adata%i_pair_addr(nr)+adata%n_pair(nr)
       endif
       if(adata%n_pair(nr)>adata%mx_pair) adata%mx_pair = adata%n_pair(nr)
    enddo ! End of nr-1,adata%n_rem_node
    return
  end subroutine trans_halo_info
!!***
end module matrix_elements_module
