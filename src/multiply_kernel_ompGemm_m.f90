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

  character(len=*), parameter :: kernel_id = "ompGemm_m"

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
    use numbers,        only: zero, one

    implicit none

    ! Passed variables
    type(matrix_halo)  :: ahalo, chalo
    type(matrix_trans) :: at
    integer      :: mx_absb, mx_part, mx_iprim, lena, lenb, lenc
    integer      :: kpart, k_off
    real(double) :: a(lena)
    real(double) :: b(lenb)
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
    integer :: i,j,k, k_in_halo, k_in_part
    integer :: nbkbeg, nabeg, naend, i_in_prim, icad, ncbeg, ncend
    integer :: tcbeg, tcend
    integer :: nb_nd_kbeg
    integer :: nd1, nd2, nd3
    real(double), allocatable, dimension(:,:) :: tempc
    integer :: maxnd1, maxnd2, maxnd3, maxlen
    integer :: nbbeg, nbend, tbend
    external :: dgemm
    ! OpenMP required indexing variables
    integer :: nd1_vector(at%mx_halo), nd2_vector(mx_absb), nd2_array(mx_absb)
    integer :: n1, n2, ncaddr

    ! Allocate tempc to largest possible size outside the loop
    maxnd1 = maxval(ahalo%ndimi)
    maxnd2 = maxval(bndim2)
    maxnd3 = maxval(ahalo%ndimj)
    maxlen = maxval(nbnab) * maxnd2
    allocate(tempc(maxnd1,maxlen))
    ! When the beta argument to dgemm is zero, C does not need to be set. See
    ! https://netlib.org/lapack/explore-html/d7/d2b/dgemm_8f_source.html

    ! Loop over atoms k in current A-halo partn
    do k = 1, ahalo%nh_part(kpart)

       call precompute_indices(jbnab2ch, nd1_vector, nd2_vector, nd2_array, &
            nd3, nbkbeg, k_in_halo, k_in_part, nb_nd_kbeg, &
            k, kpart, k_off, ahalo, chalo, at, &
            ib_nd_acc, ibaddr, nbnab, ibpart, ibseq, bndim2)

       ! Compute indices to B for dgemm call
       nd2 = bndim2(nbkbeg + nbnab(k_in_part))
       nbbeg = nd2_vector(1)
       nbend = nd2_vector(nbnab(k_in_part)) + nd3 * nd2 - 1
       tbend = nd2_array(nbnab(k_in_part)) + nd2 - 1

       ! Loop over primary-set A-neighbours of k
       !$omp do schedule(runtime)
       ANeighbours: do i = 1, at%n_hnab(k_in_halo)
          i_in_prim = at%i_prim(at%i_beg(k_in_halo) + i - 1)
          icad = (i_in_prim - 1) * chalo%ni_in_halo

          ! Compute indices to A for dgemm call
          nd1 = ahalo%ndimi(i_in_prim)
          nabeg = at%i_nd_beg(k_in_halo) + nd1_vector(i)
          naend = nabeg + (nd1 * nd3 - 1)

          ! Compute A*B using and store the result in tempc
          call dgemm('t', 'n', nd1, tbend, nd3, one, A(nabeg:naend), &
               nd3, B(nbbeg:nbend), nd3, zero, tempc, maxnd1)

          ! Copy result back from tempc and add to C
          copy_c: do j = 1, nbnab(k_in_part)
             if (jbnab2ch(j) /= 0) then
                ncbeg = chalo%i_h2d(icad+jbnab2ch(j))
                if (ncbeg /= 0) then
                   nd2 = bndim2(nbkbeg+j)
                   tcbeg = nd2_array(j)
                   tcend = tcbeg + nd2 - 1
                   do n2 = tcbeg, tcend
                      do n1 = 1, nd1
                         ncaddr = ncbeg + (n2 - tcbeg) * nd1 + n1 - 1
                         c(ncaddr) = c(ncaddr) + tempc(n1,n2)
                      end do
                   end do
                end if
             end if
          end do copy_c

       end do ANeighbours
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
    integer :: i,j,k, k_in_halo, k_in_part, nbkbeg
    integer :: nabeg, i_in_prim, icad, nbbeg, ncbeg
    integer :: n2, nbaddr, ncaddr
    integer :: n2beg, n2end
    integer :: nb_nd_kbeg
    integer :: nd1, nd2, nd3
    integer :: maxlen, maxnd1, maxnd2, maxnd3
    real(double), allocatable, dimension(:,:) :: tempb, tempc
    external :: dgemm
    ! OpenMP required indexing variables
    integer :: nd1_vector(at%mx_halo), nd2_vector(mx_absb), nd2_array(mx_absb)

    integer :: sofar

    maxnd1 = maxval(ahalo%ndimi)
    maxnd2 = maxval(bndim2)
    maxnd3 = maxval(ahalo%ndimj)
    maxlen = maxval(nbnab) * maxnd2
    allocate(tempb(maxnd3,maxlen), tempc(maxlen,maxnd1))
    ! Loop over atoms k in current A-halo partn
    do k = 1, ahalo%nh_part(kpart)

       call precompute_indices(jbnab2ch, nd1_vector, nd2_vector, nd2_array, &
            nd3, nbkbeg, k_in_halo, k_in_part, nb_nd_kbeg, &
            k, kpart, k_off, ahalo, chalo, at, &
            ib_nd_acc, ibaddr, nbnab, ibpart, ibseq, bndim2)

       ! Loop over primary-set A-neighbours of k
       !$omp do schedule(runtime)
       do i = 1, at%n_hnab(k_in_halo)
          i_in_prim = at%i_prim(at%i_beg(k_in_halo)+i-1)
          nd1 = ahalo%ndimi(i_in_prim)
          icad = (i_in_prim-1) * chalo%ni_in_halo
          nabeg = at%i_nd_beg(k_in_halo) + nd1_vector(i)
          sofar = 0
          ! Loop over B-neighbours of atom k

          do j = 1, nbnab(k_in_part)
             if (jbnab2ch(j) /= 0) then
                ncbeg = chalo%i_h2d(icad + jbnab2ch(j))
                if (ncbeg /= 0) then ! multiplication of ndim x ndim blocks
                   nd2 = bndim2(nbkbeg + j)
                   nbbeg = nd2_vector(j)

                   ! Vectorized version of temporary copies. Intel ifort compiler 2021.6
                   ! on myriad
                   ! compiles the memory copies into AVX2 instructions, but efficiency
                   ! is poor in the test cases I tried because nd1 and nd3 are small.
                   ! Tuomas Koskela (ARC) 21 Mar 2024 eCSE08 project
                   ! nbend = nbbeg + nd3 * nd2 - 1
                   ! ncend = ncbeg + nd1 * nd2 - 1
                   ! n2beg = sofar + 1
                   ! n2end = sofar + nd2
                   ! tempb(1:nd3, n2beg:n2end) = reshape(b(nbbeg:nbend), [nd3, nd2], order=[1,2])
                   ! tempc(n2beg:n2end, 1:nd1) = reshape(c(ncbeg:ncend), [nd2, nd1], order=[2,1])

                   do n2 = 1, nd2
                      nbaddr = nbbeg + nd3 * (n2 - 1)
                      ncaddr = ncbeg + nd1 * (n2 - 1)
                      tempb(1:nd3,sofar+n2) = b(nbaddr:nbaddr+nd3-1)
                      tempc(sofar+n2,1:nd1) = c(ncaddr:ncaddr+nd1-1)
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

  ! Precompute indices for parallel loop. These indices are starting indices
  ! of slices to A, B and C. The ending indices are computed on the fly.
  ! nd1_vector accumulates steps of size nd3 * ahalo%ndimi. This is used for indexing
  !    the 1d vector storing the A matrix
  ! nd2_vector accumulates steps of size nd3 * nd2. This is used for indexing
  !    1d vectors storing the B and C matrices
  ! nd2_array accumulates steps of size nd2. This is used for indexing the
  !    2d arrays storing the B and C matrices
  subroutine precompute_indices(jbnab2ch, nd1_vector, nd2_vector, nd2_array, &
       nd3, nbkbeg, k_in_halo, k_in_part, nb_nd_kbeg, &
       k, kpart, k_off, ahalo, chalo, at, &
       ib_nd_acc, ibaddr, nbnab, ibpart, ibseq, bndim2)

    use datatypes
    use matrix_module, only: matrix_halo, matrix_trans

    integer, intent(out) :: jbnab2ch(:), nd1_vector(:), nd2_vector(:), nd2_array(:)
    integer, intent(out) :: nd3, nbkbeg, k_in_halo, nb_nd_kbeg
    integer, intent(in) :: k, kpart, k_off
    integer(integ), intent(in) :: ib_nd_acc(:)
    integer(integ), intent(in) :: ibaddr(:)
    integer(integ), intent(in) :: nbnab(:)
    integer(integ), intent(in) :: ibpart(:)
    integer(integ), intent(in) :: ibseq(:)
    integer(integ), intent(in) :: bndim2(:)
    type(matrix_halo), intent(in)  :: ahalo, chalo
    type(matrix_trans), intent(in) :: at

    integer :: i, j, nd2_prev, i_in_prim_prev

    ! Indices that depend on k
    k_in_halo = ahalo%j_beg(kpart) + k - 1
    k_in_part = ahalo%j_seq(k_in_halo)
    nbkbeg = ibaddr(k_in_part) - 1
    nb_nd_kbeg = ib_nd_acc(k_in_part)
    nd3 = ahalo%ndimj(k_in_halo)

    nd1_vector(1) = 0
    nd2_vector(1) = nb_nd_kbeg
    nd2_array(1) = 1

    do i = 2, at%n_hnab(k_in_halo)
       i_in_prim_prev = at%i_prim(at%i_beg(k_in_halo) + i - 2)
       nd1_vector(i) = nd1_vector(i-1) + nd3 * ahalo%ndimi(i_in_prim_prev)
    end do

    ! transcription of j from partition to C-halo labelling
    do j = 1, nbnab(k_in_part)
       ! Also precompute jbnab2ch to be used later in the parallel loop.
       jpart = ibpart(nbkbeg + j) + k_off
       jseq = ibseq(nbkbeg + j)
       jbnab2ch(j) = chalo%i_halo(chalo%i_hbeg(jpart) + jseq - 1)

       if (j .gt. 1) then
          nd2_prev = bndim2(nbkbeg + j - 1)
          nd2_vector(j) = nd2_vector(j - 1) + nd3 * nd2_prev
          nd2_array(j) = nd2_array(j - 1) + nd2_prev
       end if

    end do

  end subroutine precompute_indices

end module multiply_kernel
