! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module multiply_kernel
! ------------------------------------------------------------------------------
! Code area 2: matrices
! ------------------------------------------------------------------------------
!
!!****h* Conquest/multiply_kernel
!! NAME
!!   multiply_kernel
!! PURPOSE
!!   Contains the matrix multiplication kernel subroutines used for
!!   Conquest matrix multiplication. This is one of the hottest part
!!   of the code.
!!   |----------------------------------|
!!   | This is the version with OpenACC |
!!   |----------------------------------|
!! AUTHOR
!!   A.Naruse
!! CREATION DATE
!!   2020/06/04
!! MODIFICATION HISTORY
!!   2020/06/04 A.Naruse
!!   - Copied multiply_kernel_default.f90 and modified it for OpenACC
!! SOURCE
!!
! ------------------------------------------------------------------------------
! Copyright (c) 2020, NVIDIA CORPORATION. All rights reserved.
!
! Permission is hereby granted, free of charge, to any person obtaining a
! copy of this software and associated documentation files (the "Software"),
! to deal in the Software without restriction, including without limitation
! the rights to use, copy, modify, merge, publish, distribute, sublicense,
! and/or sell copies of the Software, and to permit persons to whom the
! Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.
! ------------------------------------------------------------------------------
module multiply_kernel
  use nvtx

  integer, allocatable, target :: work_IJKN_0(:,:,:,:), work_IJKN_1(:,:,:,:), work_IJKN_2(:,:,:,:)
  integer, allocatable, target :: work_N_0(:), work_N_1(:), work_N_2(:)

#define arr_num_k(N)          work_N(N)
#define arr_nd3(K,N)          work_IJKN(1,          1,K,N)
#define arr_num_i(K,N)        work_IJKN(1,          2,K,N)
#define arr_num_j(K,N)        work_IJKN(1,          3,K,N)
#define arr_nabeg(I,K,N)      work_IJKN(1+(I),      1,K,N)
#define arr_nd1(I,K,N)        work_IJKN(1+(I),      2,K,N)
#define arr_icad(I,K,N)       work_IJKN(1+(I),      3,K,N)
#define arr_nbbeg(J,K,N)      work_IJKN(1+max_i+(J),1,K,N)
#define arr_nd2(J,K,N)        work_IJKN(1+max_i+(J),2,K,N)
#define arr_j_in_halo(J,K,N)  work_IJKN(1+max_i+(J),3,K,N)

!!*****

contains

  !!****f* multiply_kernel/m_kern_max *
  !!
  !!  NAME
  !!   m_kern_max
  !!  USAGE
  !!
  !!  PURPOSE
  !!   multiplication kernel for maximal case and reductions thereof.
  !!
  !!   nah_part:       the number of atoms in the A-halo that are contained
  !!                   in the current partition K.
  !!
  !!   kseq(k):        the unpruned sequence number of atom k in the current
  !!                   partition. `Unpruned' means that the sequence
  !!                   number counts all atoms in the partition, and not just
  !!                   those in the A-halo. The atom in question is specified
  !!                   by its `pruned' sequence number k.
  !!
  !!   k_in_part:      temporary variable for current value of kseq(k).
  !!
  !!   k_halo(.):      the A-halo sequence number of atom k. The latter is
  !!                   specified by its partition number and its pruned
  !!                   sequence number k in that partition.
  !!
  !!   k_in_halo:      temporary variable for current value of k_halo(.)
  !!
  !!   kbeg(kpart):    specifies address in k_halo(.) for start of information
  !!                   about atoms in partition kpart. The label kpart
  !!                   goes over all partitions containing atoms in the
  !!                   A-halo.
  !!
  !!   nahnab(k_in_halo):
  !!                   number of atoms in primary set that are A-neighbours of
  !!                   a given atom in the A-halo, the latter being given by its
  !!                   A-halo sequence number k_in_halo.
  !!
  !!   i_prim(.):      sequence number of atom i in primary set, this
  !!                   atom being specified as ith neighbour of atom
  !!                   k in the A-halo
  !!
  !!   ibeg(k_in_halo):
  !!                   address in i_prim(.) where index data for atom k
  !!                   start, the latter being specified by k_in_halo.
  !!
  !!   iaaddr(k_in_halo):
  !!                   address in array a where data for atom k start,
  !!                   the latter being specified by k_in_halo.
  !!
  !!   nbnab(k_in_part):
  !!                   number of B-neighbours of atom k, the latter being
  !!                   specified by its (unpruned) seq. no. in partition K.
  !!
  !!   ibaddr(k_in_part):
  !!                   address in array b where data for atom k start.
  !!
  !!   jch2cad(.):     address in array c where data for (i,j) pair start.
  !!
  !!   icad:           address in jch2cad where info for atom i starts.
  !!
  !!   jbnab2ch(.):    the C-halo seq. no. of atom j, the latter being
  !!                   specified as the jth B-neighbour of atom k.
  !!
  !!   ni_in_chalo:        total number of atoms in the C-halo.
  !!
  !!   ibpart(.):      partition to which each atom j belongs, the latter
  !!                   being specified as jth B-neighbour of atom k.
  !!
  !!   ibseq(.):       unpruned sequence number in partition for each
  !!                   atom j, the latter being specified as in ibpart.
  !!
  !!   ibindaddr(k_in_part):
  !!                   address in arrays ibpart(.) and ibseq(.) where
  !!                   data for atom k start.
  !!
  !!   j_halo(.):      C-halo sequence number of atom j, the latter being
  !!                   specified in partition labelling.
  !!
  !!   jchbeg(jpart):   address in array j_halo(.) where data for
  !!                   partition jpart start.
  !!
  !!   kpart_s:        starting A-halo seq. no. of current partition K
  !!
  !!   n_kpart:        number of ??? handling the same halo part (?)
  !!
  !!   side:           ID of workspace buffer to be used
  !!
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   M.J.Gillan/D.R.Bowler
  !!  CREATION DATE
  !!   25/11/99
  !!  MODIFICATION HISTORY
  !!   20/06/2001 dave
  !!    Added ROBODoc header
  !!   2020/06/04 A.Naruse
  !!    - Copied from multiply_kernel_default.f90 and refactoring for GPU
  !!  SOURCE
  !!
  subroutine m_kern_max(kpart_s, n_kpart, ib_nd_acc, ibaddr, nbnab,   &
                        ibpart, ibseq, bndim2, a, b, c, ahalo, chalo, &
                        at, mx_absb, mx_part, mx_iprim, lena, lenb,   &
                        lenc, side, debug)
    use datatypes
    use matrix_module
    use basic_types,    only: primary_set
    use primary_module, only: bundle

    implicit none

    ! Passed variables
    type(matrix_halo)  :: ahalo, chalo
    type(matrix_trans) :: at
    integer            :: mx_absb, mx_part, mx_iprim, lena, lenb, lenc
    integer            :: kpart_s, n_kpart
    real(double)       :: a(lena)
    real(double)       :: b(lenb)
    real(double)       :: c(lenc)
    integer            :: side
    integer, optional  :: debug
    ! Remote indices
    integer(integ) :: ib_nd_acc(mx_part)
    integer(integ) :: ibaddr(mx_part)
    integer(integ) :: nbnab(mx_part)
    integer(integ) :: ibpart(mx_part*mx_absb)
    integer(integ) :: ibseq(mx_part*mx_absb)
    integer(integ) :: bndim2(mx_part*mx_absb)
    ! Local variables
    integer :: jbnab2ch(mx_absb)  ! Automatic array
    integer :: nbkbeg, k, k_in_part, k_in_halo, j, jpart, jseq
    integer :: i, nabeg, i_in_prim, icad, nbbeg, j_in_halo, ncbeg
    integer :: n1, n2, n3, nb_nd_kbeg
    integer :: nd1, nd2, nd3
    integer :: naaddr, nbaddr, ncaddr
    !
    real(double) :: tmp
    integer :: l
    integer :: max_i, max_j, max_k, max_nd1, max_nd2, max_nd3
    integer, pointer :: work_IJKN(:,:,:,:)
    integer, pointer :: work_N(:)

    call nvtxStartRange('m_kern_max', __LINE__)

    call m_kern_prep(kpart_s, n_kpart, ib_nd_acc, ibaddr, nbnab, &
                     ibpart, ibseq, bndim2, ahalo, chalo, at, &
                     mx_absb, mx_part, mx_iprim, side, &
                     max_i, max_j, max_k, max_nd1, max_nd2, max_nd3)

    if (side == 0) then
       work_IJKN => work_IJKN_0
       work_N    => work_N_0
    else if (side == 1) then
       work_IJKN => work_IJKN_1
       work_N    => work_N_1
    else if (side == 2) then
       work_IJKN => work_IJKN_2
       work_N    => work_N_2
    else
       call cq_abort('side must be 0, 1 or 2.')
    endif

    !$acc parallel vector_length(32) copy(c) copyin(a, b) copyin(chalo, chalo%i_h2d, work_ijkn, work_n) async(side)
    !$acc loop collapse(3) gang
    do l = 1, n_kpart
       ! Loop over B-neighbours of atom k
       do j = 1, max_j
          ! Loop over primary-set A-neighbours of k
          do i = 1, max_i
             !$acc loop collapse(3) vector
             ! Loop over atoms k in current A-halo partn
             do k = 1, max_k
                ! multiplication of ndim x ndim blocks
                do n2 = 1, max_nd2
                   do n1 = 1, max_nd1
                      if (k > arr_num_k(l)) cycle
                      if (j > arr_num_j(k,l)) cycle
                      if (i > arr_num_i(k,l)) cycle
                      j_in_halo = arr_j_in_halo(j,k,l)
                      if (j_in_halo == 0) cycle
                      icad = arr_icad(i,k,l)
                      ncbeg = chalo%i_h2d(icad+j_in_halo)
                      if (ncbeg == 0 ) cycle
                      nd1 = arr_nd1(i,k,l)
                      nd2 = arr_nd2(j,k,l)
                      nd3 = arr_nd3(k,l)
                      if (n1 > nd1) cycle
                      if (n2 > nd2) cycle
                      naaddr = arr_nabeg(i,k,l) + nd3 * (n1 - 1)
                      nbaddr = arr_nbbeg(j,k,l) + nd3 * (n2 - 1)
                      ncaddr = ncbeg + nd1 * (n2 - 1) + (n1 - 1)
                      tmp = 0.d0
                      !$acc loop seq
                      do n3 = 1, nd3
                         tmp = tmp + a(naaddr+n3-1) * b(nbaddr+n3-1)
                      end do
                      !$acc atomic update
                      c(ncaddr) = c(ncaddr) + tmp
                      !$acc end atomic
                   end do
                end do
             end do ! End of k = 1, max_k
          end do ! End of i = 1, max_i
       end do ! End of j = 1, max_j
    end do
    !$acc end parallel

    call nvtxEndRange
    return
  end subroutine m_kern_max
  !!*****


  !!****f* multiply_kernel/m_kern_min *
  !!
  !!  NAME
  !!   m_kern_min
  !!  USAGE
  !!
  !!  PURPOSE
  !!   multiplication kernel for minimal case and extensions thereof.
  !!
  !!   nah_part:       the number of atoms in the A-halo that are contained
  !!                   in the current partition K.
  !!
  !!   kseq(k):        the unpruned sequence number of atom k in the current
  !!                   partition. `Unpruned' means that the sequence
  !!                   number counts all atoms in the partition, and not just
  !!                   those in the A-halo. The atom in question is specified
  !!                   by its `pruned' sequence number k.
  !!
  !!   k_in_part:      temporary variable for current value of kseq(k).
  !!
  !!   k_halo(.):      the A-halo sequence number of atom k. The latter is
  !!                   specified by its partition number and its pruned
  !!                   sequence number k in that partition.
  !!
  !!   k_in_halo:      temporary variable for current value of k_halo(.)
  !!
  !!   kbeg(kpart):    specifies address in k_halo(.) for start of information
  !!                   about atoms in partition kpart. The label kpart
  !!                   goes over all partitions containing atoms in the
  !!                   A-halo.
  !!
  !!   nahnab(k_in_halo):
  !!                   number of atoms in primary set that are A-neighbours of
  !!                   a given atom in the A-halo, the latter being given by its
  !!                   A-halo sequence number k_in_halo.
  !!
  !!   i_prim(.):      sequence number of atom i in primary set, this
  !!                   atom being specified as ith neighbour of atom
  !!                   k in the A-halo
  !!
  !!   ibeg(k_in_halo):
  !!                   address in i_prim(.) where index data for atom k
  !!                   start, the latter being specified by k_in_halo.
  !!
  !!   iaaddr(k_in_halo):
  !!                   address in array a where data for atom k start,
  !!                   the latter being specified by k_in_halo.
  !!
  !!   nbnab(k_in_part):
  !!                   number of B-neighbours of atom k, the latter being
  !!                   specified by its (unpruned) seq. no. in partition K.
  !!
  !!   ibaddr(k_in_part):
  !!                   address in array b where data for atom k start.
  !!
  !!   jch2cad(.):     address in array ! where data for (i,j) pair start.
  !!
  !!   icad:           address in jch2cad where info for atom i starts.
  !!
  !!   jbnab2ch(.):    the C-halo seq. no. of atom j, the latter being
  !!                   specified as the jth B-neighbour of atom k.
  !!
  !!   ni_in_chalo:        total number of atoms in the C-halo.
  !!
  !!   ibpart(.):      partition to which each atom j belongs, the latter
  !!                   being specified as jth B-neighbour of atom k.
  !!
  !!   ibseq(.):       unpruned sequence number in partition for each
  !!                   atom j, the latter being specified as in ibpart.
  !!
  !!   ibindaddr(k_in_part):
  !!                   address in arrays ibpart(.) and ibseq(.) where
  !!                   data for atom k start.
  !!
  !!   j_halo(.):      C-halo sequence number of atom j, the latter being
  !!                   specified in partition labelling.
  !!
  !!   jchbeg(jpart):   address in array j_halo(.) where data for
  !!                   partition jpart start.
  !!
  !!   kpart_s:        starting A-halo seq. no. of current partition K
  !!
  !!   n_kpart:        number of ??? handling the same halo part (?)
  !!
  !!   side:           ID of workspace buffer to be used
  !!
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   M.J.Gillan/D.R.Bowler
  !!  CREATION DATE
  !!   25/11/99
  !!  MODIFICATION HISTORY
  !!   20/06/2001 dave
  !!    Added ROBODoc header
  !!   2020/06/04 A.Naruse
  !!    - Copied from multiply_kernel_default.f90 and refactoring for GPU
  !!  SOURCE
  !!
  subroutine m_kern_min(kpart_s, n_kpart, ib_nd_acc, ibaddr, nbnab,   &
                        ibpart, ibseq, bndim2, a, b, c, ahalo, chalo, &
                        at, mx_absb, mx_part, mx_iprim, lena, lenb,   &
                        lenc, side)
    use datatypes
    use matrix_module
    use basic_types,    only: primary_set
    use primary_module, only: bundle

    implicit none

    ! Passed variables
    type(matrix_halo)  :: ahalo, chalo
    type(matrix_trans) :: at
    integer            :: mx_absb, mx_part, mx_iprim, lena, lenb, lenc
    ! Remember that a is a local transpose
    real(double)       :: a(lena)
    real(double)       :: b(lenb)
    real(double)       :: c(lenc)
    integer            :: kpart_s, n_kpart
    integer            :: side
    ! dimension declarations
    integer :: ibaddr(mx_part)
    integer :: ib_nd_acc(mx_part)
    integer :: nbnab(mx_part)
    integer :: ibpart(mx_part*mx_absb)
    integer :: ibseq(mx_part*mx_absb)
    integer :: bndim2(mx_part*mx_absb)
    ! Local variables
    integer :: jbnab2ch(mx_absb)
    integer :: k, k_in_part, k_in_halo, nbkbeg, j, jpart, jseq
    integer :: i, nabeg, i_in_prim, icad, nbbeg, j_in_halo, ncbeg
    integer :: n1, n2, n3, nb_nd_kbeg
    integer :: nd1, nd2, nd3
    integer :: naaddr, nbaddr, ncaddr
    !
    real(double) :: tmp
    integer :: l
    integer :: max_i, max_j, max_k, max_nd1, max_nd2, max_nd3
    integer, pointer :: work_IJKN(:,:,:,:)
    integer, pointer :: work_N(:)

    call nvtxStartRange('m_kern_min', __LINE__)

    call m_kern_prep(kpart_s, n_kpart, ib_nd_acc, ibaddr, nbnab, &
                     ibpart, ibseq, bndim2, ahalo, chalo, at, &
                     mx_absb, mx_part, mx_iprim, side, &
                     max_i, max_j, max_k, max_nd1, max_nd2, max_nd3)

    if (side == 0) then
       work_IJKN => work_IJKN_0
       work_N    => work_N_0
    else if (side == 1) then
       work_IJKN => work_IJKN_1
       work_N    => work_N_1
    else if (side == 2) then
       work_IJKN => work_IJKN_2
       work_N    => work_N_2
    else
       call cq_abort('side must be 0, 1 or 2.')
    endif

    !$acc parallel vector_length(32) copy(a) copyin(b, c) copyin(chalo, chalo%i_h2d, work_ijkn, work_n) async(side)
    !$acc loop collapse(3) gang
    do l = 1, n_kpart
       ! Loop over B-neighbours of atom k
       do j = 1, max_j
          ! Loop over primary-set A-neighbours of k
          do i = 1, max_i
             !$acc loop collapse(3) vector
             ! Loop over atoms k in current A-halo partn
             do k = 1, max_k
                ! multiplication of ndim x ndim blocks
                do n3 = 1, max_nd3
                   do n1 = 1, max_nd1
                      if (k > arr_num_k(l)) cycle
                      if (j > arr_num_j(k,l)) cycle
                      if (i > arr_num_i(k,l)) cycle
                      j_in_halo = arr_j_in_halo(j,k,l)
                      if (j_in_halo == 0) cycle
                      icad = arr_icad(i,k,l)
                      ncbeg = chalo%i_h2d(icad+j_in_halo)
                      if (ncbeg == 0) cycle
                      nd1 = arr_nd1(i,k,l)
                      nd2 = arr_nd2(j,k,l)
                      nd3 = arr_nd3(k,l)
                      if (n1 > nd1) cycle
                      if (n3 > nd3) cycle
                      naaddr = arr_nabeg(i,k,l) + nd3 * (n1 - 1) + (n3 - 1)
                      nbaddr = arr_nbbeg(j,k,l) + (n3 - 1)
                      ncaddr = ncbeg + (n1 - 1)
                      tmp = 0.d0
                      !$acc loop seq
                      do n2 = 1, nd2
                         tmp = tmp + c(ncaddr+nd1*(n2-1)) * b(nbaddr+nd3*(n2-1))
                      end do
                      !$acc atomic update
                      a(naaddr) = a(naaddr) + tmp
                      !$acc end atomic
                   end do
                end do
             end do
          end do
       end do
    end do
    !$acc end parallel

    call nvtxEndRange
    return
  end subroutine m_kern_min
  !!*****

  !
  subroutine m_kern_prep(kpart_s, n_kpart, ib_nd_acc, ibaddr, nbnab, &
                         ibpart, ibseq, bndim2, ahalo, chalo, at, &
                         mx_absb, mx_part, mx_iprim, side, &
                         max_i, max_j, max_k, max_nd1, max_nd2, max_nd3)
     ! calculate array/matrix indices used in m_kern_max() and m_kern_min()
     use datatypes
     use matrix_module
     use GenComms, only: cq_abort

     implicit none

     ! Passed variables
     integer :: kpart_s, n_kpart
     integer :: mx_absb, mx_part, mx_iprim
     integer :: side
     integer(integ) :: ib_nd_acc(mx_part)
     integer(integ) :: ibaddr(mx_part)
     integer(integ) :: nbnab(mx_part)
     integer(integ) :: ibpart(mx_part*mx_absb)
     integer(integ) :: ibseq(mx_part*mx_absb)
     integer(integ) :: bndim2(mx_part*mx_absb)
     type(matrix_halo) :: ahalo, chalo
     type(matrix_trans) :: at
     integer :: max_i, max_j, max_k
     integer :: max_nd1, max_nd2, max_nd3
     ! Local variables
     integer :: kpart, k_off, k_in_halo, k_in_part
     integer :: nabeg, i_in_prim, icad, nbbeg, nbkbeg, jpart, jseq, j_in_halo
     integer :: nd1, nd2, nd3
     integer :: l, i, j, k, num_i, num_j, num_k
     integer, pointer :: work_IJKN(:,:,:,:)
     integer, pointer :: work_N(:)

     max_i = 0
     max_j = 0
     max_k = 0
     do l = 1, n_kpart
        kpart = kpart_s + l - 1

        num_k = ahalo%nh_part(kpart)
        max_k = max(max_k, num_k)
        do k = 1, num_k
           k_in_halo = ahalo%j_beg(kpart) + k - 1
           k_in_part = ahalo%j_seq(k_in_halo)

           num_i = at%n_hnab(k_in_halo)
           max_i = max(max_i, num_i)

           num_j = nbnab(k_in_part)
           max_j = max(max_j, num_j)
        enddo
     enddo

     !$acc wait(side)
     if (side == 0) then
        if (allocated(work_IJKN_0)) then
           deallocate( work_IJKN_0, work_N_0 )
        endif
        allocate( work_IJKN_0(1+max_i+max_j, 3, max_k, n_kpart) )
        allocate( work_N_0(n_kpart) )
        work_IJKN => work_IJKN_0
        work_N    => work_N_0
     else if (side == 1) then
        if (allocated(work_IJKN_1)) then
           deallocate( work_IJKN_1, work_N_1 )
        endif
        allocate( work_IJKN_1(1+max_i+max_j, 3, max_k, n_kpart) )
        allocate( work_N_1(n_kpart) )
        work_IJKN => work_IJKN_1
        work_N    => work_N_1
     else if (side == 2) then
        if (allocated(work_IJKN_2)) then
           deallocate( work_IJKN_2, work_N_2 )
        endif
        allocate( work_IJKN_2(1+max_i+max_j, 3, max_k, n_kpart) )
        allocate( work_N_2(n_kpart) )
        work_IJKN => work_IJKN_2
        work_N    => work_N_2
     else
        call cq_abort('side must be 0, 1 or 2.')
     endif
     work_IJKN(:,:,:,:) = 0
     work_N(:) = 0

     max_nd1 = 0
     max_nd2 = 0
     max_nd3 = 0
     do l = 1, n_kpart
        kpart = kpart_s + l - 1
        k_off = ahalo%lab_hcover(kpart)

        num_k = ahalo%nh_part(kpart)
        arr_num_k(l) = num_k
        !max_k = max(max_k, num_k)
        do k = 1, num_k
           k_in_halo = ahalo%j_beg(kpart) + k - 1
           k_in_part = ahalo%j_seq(k_in_halo)

           nd3 = ahalo%ndimj(k_in_halo)
           arr_nd3(k,l) = nd3
           max_nd3 = max(max_nd3, nd3)
           nabeg = at%i_nd_beg(k_in_halo)

           num_i = at%n_hnab(k_in_halo)
           arr_num_i(k,l) = num_i
           !max_i = max(max_i, num_i)
           do i = 1, num_i
              arr_nabeg(i,k,l) = nabeg
              i_in_prim = at%i_prim(at%i_beg(k_in_halo)+i-1)
              icad = (i_in_prim-1) * chalo%ni_in_halo
              arr_icad(i,k,l) = icad
              nd1 = ahalo%ndimi(i_in_prim)
              arr_nd1(i,k,l) = nd1
              max_nd1 = max(max_nd1, nd1)
              nabeg = nabeg + nd3 * nd1
           enddo

           nbkbeg = ibaddr(k_in_part)
           nbbeg = ib_nd_acc(k_in_part)

           num_j = nbnab(k_in_part)
           arr_num_j(k,l) = num_j
           !max_j = max(max_j, num_j)
           do j = 1, num_j
              arr_nbbeg(j,k,l) = nbbeg
              nd2 = bndim2(nbkbeg+j-1)
              arr_nd2(j,k,l) = nd2
              max_nd2 = max(max_nd2, nd2)
              nbbeg = nbbeg + nd3 * nd2
              jpart = ibpart(nbkbeg+j-1) + k_off
              jseq = ibseq(nbkbeg+j-1)
              j_in_halo = chalo%i_halo(chalo%i_hbeg(jpart)+jseq-1)
              arr_j_in_halo(j,k,l) = j_in_halo
           enddo
        enddo
     enddo

     return
  end subroutine m_kern_prep
  !

end module multiply_kernel
