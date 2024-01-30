! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module multiply_kernel
! ------------------------------------------------------------------------------
! Code area 2: matrices
! ------------------------------------------------------------------------------

!!****h* Conquest/multiply_kernel
!! NAME
!!   multiply_kernel
!! PURPOSE
!!   Contains the matrix multiplication kernel subroutines used for
!!   Conquest matrix multiplication. This is one of the hottest part
!!   of the code.
!!   |---------------------------------------------|
!!   | This is the version using GEMM library call |
!!   |---------------------------------------------|
!! AUTHOR
!!   L.Tong
!! CREATION DATE
!!   2012/08/28
!! MODIFICATION HISTORY
!! SOURCE
!!
module multiply_kernel

!!*****

contains

  !!****f* multiply_kernel/m_kern_max
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
  !!   kpart:          A-halo seq. no. of current partition K
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
  !!   k_off:          offset in partition labelling to account for p.b.c.
  !!
  !!   j_halo(.):      C-halo sequence number of atom j, the latter being
  !!                   specified in partition labelling.
  !!
  !!   jchbeg(jpart):   address in array j_halo(.) where data for
  !!                   partition jpart start.
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
  !!  SOURCE
  !!
  subroutine m_kern_max(k_off, kpart, ib_nd_acc, ibaddr, nbnab,       &
                        ibpart, ibseq, bndim2, a, b, c, ahalo, chalo, &
                        at, mx_absb, mx_part, mx_iprim, lena, lenb,   &
                        lenc, debug)
    use datatypes
    use matrix_module
    use basic_types,    only: primary_set
    use primary_module, only: bundle
    use numbers,        only: zero, one

    implicit none

    ! Passed variables
    type(matrix_halo)  :: ahalo, chalo
    type(matrix_trans) :: at
    integer      :: mx_absb, mx_part, mx_iprim, lena, lenb, lenc
    integer      :: kpart, k_off
    real(double), target :: a(lena)
    real(double), target :: b(lenb)
    real(double) :: c(lenc)
    integer, optional :: debug
    ! Remote indices
    integer(integ), intent(in) :: ib_nd_acc(:)
    integer(integ), intent(in) :: ibaddr(:)
    integer(integ), intent(in) :: nbnab(:)
    integer(integ), intent(in) :: ibpart(:)
    integer(integ), intent(in) :: ibseq(:)
    integer(integ), intent(in) :: bndim2(:)
    ! Local variables
    integer :: jbnab2ch(mx_absb)  ! Automatic array
    integer :: nbkbeg, k, k_in_part, k_in_halo, j, jpart, jseq
    integer :: i, nabeg, naend, i_in_prim, icad, j_in_halo, ncbeg, ncend
    integer :: tcbeg, tcend
    integer :: n1, n2, n3, nb_nd_kbeg
    integer :: nd1, nd2, nd3
    integer :: naaddr, nbaddr, ncaddr
    real(double), allocatable, dimension(:,:) :: tempb, tempc
    real(double), pointer, dimension(:,:) :: pointa, pointb
    integer :: sofar, maxnd1, maxnd2, maxnd3, maxlen
    integer :: nbbeg, nbend, tbbeg, tbend
    logical :: mask
    external :: dgemm
    ! OpenMP required indexing variables
    integer :: nd1_vector(at%mx_halo), nd2_vector(mx_absb), nd2_array(mx_absb)

    ! Allocate tempb, tempc to largest possible size outside the loop
    maxnd1 = maxval(ahalo%ndimi)
    maxnd2 = maxval(bndim2)
    maxnd3 = maxval(ahalo%ndimj)
    maxlen = maxval(nbnab) * maxnd2
    allocate(tempc(maxnd1,maxlen))
    ! We don't need to initialize these because:
    ! All elements of tempB are overwitten
    ! When the beta argument to dgemm is zero, C does not need to be set. See
    ! https://netlib.org/lapack/explore-html/d7/d2b/dgemm_8f_source.html

    ! Loop over atoms k in current A-halo partn
    do k = 1, ahalo%nh_part(kpart)
       ! These indices only depend on k and kpart
       k_in_halo = ahalo%j_beg(kpart) + k - 1
       k_in_part = ahalo%j_seq(k_in_halo)
       nbkbeg = ibaddr(k_in_part)
       nb_nd_kbeg = ib_nd_acc(k_in_part)
       nd3 = ahalo%ndimj(k_in_halo)

       ! Precompute indices that depend on i and j to avoid loop carrier
       ! dependency in the following loops
       nd1_vector(1) = 0
       do i = 2, at%n_hnab(k_in_halo)
         i_in_prim = at%i_prim(at%i_beg(k_in_halo) + i - 2)
         nd1_vector(i) = nd1_vector(i-1) + nd3 * ahalo%ndimi(i_in_prim)
       end do
       nd2_vector(1) = nb_nd_kbeg
       nd2_array(1) = 1
       do j = 2, nbnab(k_in_part)
          nd2 = bndim2(nbkbeg + j - 2)
          nd2_vector(j) = nd2_vector(j - 1) + nd3 * nd2
          nd2_array(j) = nd2_array(j - 1) + nd2
       end do

       ! transcription of j from partition to C-halo labelling
       copy_b: do j = 1, nbnab(k_in_part)
          jpart = ibpart(nbkbeg + j - 1) + k_off
          jseq = ibseq(nbkbeg + j - 1)
          jbnab2ch(j) = chalo%i_halo(chalo%i_hbeg(jpart) + jseq - 1)

          nd2 = bndim2(nbkbeg+j-1)
          nbend = nd2_vector(j) + nd3 * nd2 - 1
          tbend = nd2_array(j) + nd2 - 1
       end do copy_b

       pointb(1:nd3, 1:tbend) => b(nb_nd_kbeg:nbend)

       ! Loop over primary-set A-neighbours of k
       !$omp do schedule(runtime)
       do i = 1, at%n_hnab(k_in_halo)
          i_in_prim = at%i_prim(at%i_beg(k_in_halo)+i-1)
          icad = (i_in_prim - 1) * chalo%ni_in_halo
          nd1 = ahalo%ndimi(i_in_prim)
          nabeg = at%i_nd_beg(k_in_halo) + nd1_vector(i)
          naend = nabeg + (nd1 * nd3 - 1)

          pointa(1:nd3, 1:nd1) => a(nabeg:naend)

          call dgemm('t', 'n', nd1, tbend, nd3, one, pointa, &
               nd3, pointb, nd3, zero, tempc, maxnd1)

          ! Loop over B-neighbours of atom k
          ! Copy result back from temporary array and add to C
          copy_c: do j = 1, nbnab(k_in_part)
             nd2 = bndim2(nbkbeg+j-1)
             ncbeg = chalo%i_h2d(icad+jbnab2ch(j))
             ncend = ncbeg + (nd1 * nd2 - 1)
             mask = (jbnab2ch(j) /= 0 .and. ncbeg /= 0)
             tcbeg = nd2_array(j)
             tcend = tcbeg + nd2 - 1
             c(ncbeg:ncend) = c(ncbeg:ncend) + pack(tempc(1:nd1, tcbeg:tcend), mask)
          end do copy_c

          ! copy_c: do j = 1, nbnab(k_in_part)
          !    nd2 = bndim2(nbkbeg+j-1)
          !    ncbeg = chalo%i_h2d(icad+jbnab2ch(j))
          !    ncend = ncbeg + (nd1 * nd2 - 1)
          !    if (jbnab2ch(j) /= 0 .and. ncbeg /= 0) then
          !       tcbeg = nd2_array(j)
          !       tcend = tcbeg + nd2 - 1
          !       c(ncbeg:ncend) = c(ncbeg:ncend) + pack(tempc(1:nd1, tcbeg:tcend), mask)
          !    end if
          ! end do copy_c

          ! sofar = 0
          ! ! Loop over B-neighbours of atom k
          ! do j = 1, nbnab(k_in_part)
          !    nd2 = bndim2(nbkbeg+j-1)
          !    j_in_halo = jbnab2ch(j)
          !    if (j_in_halo /= 0) then
          !       ncbeg = chalo%i_h2d(icad+j_in_halo)
          !       if (ncbeg /= 0) then ! multiplication of ndim x ndim blocks
          !          do n2 = 1, nd2
          !             ncaddr = ncbeg + nd1 * (n2 - 1)
          !             do n1 = 1, nd1
          !                c(ncaddr+n1-1) = c(ncaddr+n1-1) + tempc(n1,sofar+n2)
          !             end do
          !          end do
          !       end if
          !    end if
          !    sofar = sofar + nd2
          ! end do ! end of j = 1, nbnab(k_in_part)

       end do ! end of i = 1, at%n_hnab
       !$omp end do
    end do ! end of k = 1, nahpart
    deallocate(tempc)
    return
  end subroutine m_kern_max
  !!*****


  !!****f* multiply_kernel/m_kern_min
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
  !!   kpart:          A-halo seq. no. of current partition K
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
  !!   k_off:          offset in partition labelling to account for p.b.c.
  !!
  !!   j_halo(.):      C-halo sequence number of atom j, the latter being
  !!                   specified in partition labelling.
  !!
  !!   jchbeg(jpart):   address in array j_halo(.) where data for
  !!                   partition jpart start.
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
  !!  SOURCE
  !!
  subroutine m_kern_min(k_off, kpart, ib_nd_acc, ibaddr, nbnab,       &
                        ibpart, ibseq, bndim2, a, b, c, ahalo, chalo, &
                        at, mx_absb, mx_part, mx_iprim, lena, lenb,   &
                        lenc)
    use datatypes
    use matrix_module
    use basic_types,    only: primary_set
    use primary_module, only: bundle
    use numbers,        only: one, zero

    implicit none

    ! Passed variables
    type(matrix_halo)  :: ahalo, chalo
    type(matrix_trans) :: at
    integer :: mx_absb, mx_part, mx_iprim, lena, lenb, lenc
    integer :: kpart, k_off
    ! Remember that a is a local transpose
    real(double) :: a(lena)
    real(double) :: b(lenb)
    real(double) :: c(lenc)
    ! dimension declarations
    integer(integ), intent(in) :: ib_nd_acc(:)
    integer(integ), intent(in) :: ibaddr(:)
    integer(integ), intent(in) :: nbnab(:)
    integer(integ), intent(in) :: ibpart(:)
    integer(integ), intent(in) :: ibseq(:)
    integer(integ), intent(in) :: bndim2(:)
    ! Local variables
    integer :: jbnab2ch(mx_absb)
    integer :: k, k_in_part, k_in_halo, nbkbeg, j, jpart, jseq
    integer :: i, nabeg, i_in_prim, icad, nbbeg, j_in_halo, ncbeg
    integer :: n1, n2, n3, nb_nd_kbeg
    integer :: nd1, nd2, nd3
    integer :: naaddr, nbaddr, ncaddr
    integer :: sofar, maxnd1, maxnd2, maxnd3, maxlen
    real(double), allocatable, dimension(:,:) :: tempb, tempc
    external :: dgemm
    ! OpenMP required indexing variables
    integer :: nd1_1st(at%mx_halo), nd2_1st(mx_absb)

    maxnd1 = maxval(ahalo%ndimi)
    maxnd2 = maxval(bndim2)
    maxnd3 = maxval(ahalo%ndimj)
    maxlen = maxval(nbnab) * maxnd2
    allocate(tempb(maxnd3,maxlen), tempc(maxlen,maxnd1))
    tempb = zero
    tempc = zero
    ! Loop over atoms k in current A-halo partn
    do k = 1, ahalo%nh_part(kpart)
       k_in_halo = ahalo%j_beg(kpart) + k - 1
       k_in_part = ahalo%j_seq(k_in_halo)
       nbkbeg = ibaddr(k_in_part)
       nb_nd_kbeg = ib_nd_acc(k_in_part)
       nd3 = ahalo%ndimj(k_in_halo)
       ! for OpenMP sub-array indexing
       nd1_1st(1) = 0
       do i = 2, at%n_hnab(k_in_halo)
         i_in_prim = at%i_prim(at%i_beg(k_in_halo)+i-2)
         nd1_1st(i) = nd1_1st(i-1) + nd3 * ahalo%ndimi(i_in_prim)
       end do
       nd2_1st(1) = 0
       do j = 2, nbnab(k_in_part)
          nd2_1st(j) = nd2_1st(j-1) + nd3 * bndim2(nbkbeg+j-2)
       end do
       ! transcription of j from partition to C-halo labelling
       do j = 1, nbnab(k_in_part)
          jpart = ibpart(nbkbeg+j-1) + k_off
          jseq = ibseq(nbkbeg+j-1)
          jbnab2ch(j) = chalo%i_halo(chalo%i_hbeg(jpart)+jseq-1)
       end do
!$omp do schedule(runtime)
       ! Loop over primary-set A-neighbours of k
       do i = 1, at%n_hnab(k_in_halo)
          ! nabeg = at%i_beg(k_in_halo) + i - 1
          i_in_prim = at%i_prim(at%i_beg(k_in_halo)+i-1)
          nd1 = ahalo%ndimi(i_in_prim)
          nabeg = at%i_nd_beg(k_in_halo) + nd1_1st(i)
          icad = (i_in_prim-1) * chalo%ni_in_halo
          sofar = 0
          ! Loop over B-neighbours of atom k
          do j = 1, nbnab(k_in_part)
             ! nbbeg = nbkbeg + j - 1
             nd2 = bndim2(nbkbeg+j-1)
             nbbeg = nb_nd_kbeg + nd2_1st(j)
             j_in_halo = jbnab2ch(j)
             if (j_in_halo /= 0) then
                ! nd2 = chalo%ndimj(j_in_halo)
                ncbeg = chalo%i_h2d(icad+j_in_halo)
                if (ncbeg /= 0) then ! multiplication of ndim x ndim blocks
!DIR$ NOPATTERN
                   do n2 = 1, nd2
                      nbaddr = nbbeg + nd3 * (n2 - 1)
                      ncaddr = ncbeg + nd1 * (n2 - 1)
                      do n3 = 1, nd3
                         tempb(n3,sofar+n2) = b(nbaddr+n3-1)
                      end do
                      do n1 = 1, nd1
                         tempc(sofar+n2,n1) = c(ncaddr+n1-1)
                      end do
                   end do
                   sofar = sofar + nd2
                end if
             end if
          end do
          if (sofar > 0) then
             ! m, n, k, alpha, a, lda, b, ldb, beta, c, ldc
             call dgemm('n', 'n', nd3, nd1, sofar, one, tempb, &
                        maxnd3, tempc, maxlen, one, a(nabeg:), nd3)
          end if
       end do
!$omp end do
    end do
    deallocate(tempb, tempc)
    return
  end subroutine m_kern_min
  !!*****

end module multiply_kernel
