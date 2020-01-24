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
!!   |---------------------------------------|
!!   | This is the version with OpenMP Tasks |
!!   |---------------------------------------|
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
  !!  USES
  !!
  !!  AUTHOR
  !!   M.J.Gillan/D.R.Bowler
  !!  CREATION DATE
  !!   25/11/99
  !!  MODIFICATION HISTORY
  !!   20/06/2001 dave
  !!    Added ROBODoc header
  !!   2012/07/05 L.Tong
  !!   - Added OpenMP directives
  !!   2012/08/02 L.Tong
  !!   - Realsed that it is not possible to achieve maximum
  !!     efficiency with one piece of code for both OpenMP and standard
  !!     versions, so decided to use C preprocessors. The preprocessors
  !!     are supported by all major compilers and should not pose much
  !!     problem on different platiforms.
  !!   - The standard code uses less automatic arrays for indexing.
  !!     But for OpenMP calculations, extra indexing arrays has to be
  !!     allocated for different threads. While automatic arrays uses
  !!     stack and allocating and deallocating should be efficient,
  !!     memory allocations for this hot part of the calculation
  !!     should be reduced to minimum. Hence the original indexing
  !!     method in the routine is kept for non-OpenMP versions, while
  !!     a new indexing method with extra indexing arrays is used for
  !!     OpenMP version.
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
    use omp_lib

    implicit none

    ! Passed variables
    type(matrix_halo)  :: ahalo, chalo
    type(matrix_trans) :: at
    integer            :: mx_absb, mx_part, mx_iprim, lena, lenb, lenc
    integer            :: kpart, k_off
    real(double)       :: a(lena)
    real(double)       :: b(lenb)
    real(double)       :: c(lenc)
    integer, optional  :: debug
    ! Remote indices
    integer(integ) :: ib_nd_acc(mx_part)
    integer(integ) :: ibaddr(mx_part)
    integer(integ) :: nbnab(mx_part)
    integer(integ) :: ibpart(mx_part*mx_absb)
    integer(integ) :: ibseq(mx_part*mx_absb)
    integer(integ) :: bndim2(mx_part*mx_absb)
    ! Local variables
    integer :: nbkbeg, k, k_in_part, k_in_halo, j, jpart, jseq
    integer :: i, nabeg, i_in_prim, icad, nbbeg, j_in_halo, ncbeg
    integer :: n1, n2, n3, nb_nd_kbeg
    integer :: nd1, nd2, nd3
    integer :: naaddr, nbaddr, ncaddr
    ! The OpenMP firstprivate clause for some reason crashes (at least
    ! for gcc) whenever there is automatic arrays in its
    ! list. Everthing works fine with allocatable arrays. The
    ! allocation of allocatable arrays may be slightly slower than
    ! automatic arrays, but it is the only option at the moment.
    integer, dimension(:), allocatable :: nd1_1st, nd2_1st, jbnab2ch

    allocate(nd1_1st(at%mx_halo), nd2_1st(mx_absb), jbnab2ch(mx_absb))

!$omp parallel default(none)                                             &
!$omp          shared(kpart, ibaddr, ib_nd_acc, nbnab, ibpart, ibseq,    &
!$omp                 k_off, bndim2, mx_absb, mx_part, at, ahalo, chalo, &
!$omp                 a, b, c)                                           &
!$omp          private(i, j, k, j_in_halo, k_in_halo, k_in_part, nbkbeg, &
!$omp                  nb_nd_kbeg, nd1, nd2, nd3, jpart, jseq, jbnab2ch, &
!$omp                  nabeg, nbbeg, ncbeg, i_in_prim, icad, naaddr,     &
!$omp                  nbaddr, ncaddr, n1, n2, n3, nd1_1st, nd2_1st)
    ! run on single thread for generating tasks
!$omp single
    ! task for generating computation tasks
!$omp task untied                                                    &
!$omp      default(none)                                             &
!$omp      shared(kpart, ibaddr, ib_nd_acc, nbnab, ibpart, ibseq,    &
!$omp             k_off, bndim2, mx_absb, mx_part, at, ahalo, chalo, &
!$omp             a, b, c)                                           &
!$omp      private(i, j, k, j_in_halo, k_in_halo, k_in_part, nbkbeg, &
!$omp              nb_nd_kbeg, nd1, nd2, nd3, jpart, jseq, jbnab2ch, &
!$omp              nabeg, nbbeg, ncbeg, i_in_prim, icad, naaddr,     &
!$omp              nbaddr, ncaddr, n1, n2, n3, nd1_1st, nd2_1st)
    ! Loop over atoms k in current A-halo partn
    do k = 1, ahalo%nh_part(kpart)
       k_in_halo = ahalo%j_beg(kpart) + k - 1
       k_in_part = ahalo%j_seq(k_in_halo)
       nbkbeg = ibaddr(k_in_part)
       nb_nd_kbeg = ib_nd_acc(k_in_part)
       nd3 = ahalo%ndimj(k_in_halo)

       ! if(PRESENT(debug)) write(21+debug,*) 'Details1: ',k,nb_nd_kbeg

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
       ! Loop over primary-set A-neighbours of k
       do i = 1, at%n_hnab(k_in_halo)
          i_in_prim = at%i_prim(at%i_beg(k_in_halo)+i-1)
          nd1 = ahalo%ndimi(i_in_prim)
          icad = (i_in_prim-1) * chalo%ni_in_halo
          nabeg = at%i_nd_beg(k_in_halo) + nd1_1st(i)
          ! Loop over B-neighbours of atom k
          do j = 1, nbnab(k_in_part)
             nd2 = bndim2(nbkbeg+j-1)
             nbbeg = nb_nd_kbeg + nd2_1st(j)
             j_in_halo = jbnab2ch(j)
!$omp task if(j_in_halo /= 0)                                        &
!$omp      default(none)                                             &
!$omp      shared(nbnab, bndim2, chalo, a, b, c)                     &
!$omp      firstprivate(nd1, nd2, nd3, nbkbeg, nb_nd_kbeg, nd2_1st,  &
!$omp                   jbnab2ch, nabeg, icad, k_in_part, j_in_halo, &
!$omp                   nbbeg)                                       &
!$omp      private(ncbeg, n1, n2, n3, naaddr, nbaddr, ncaddr)
             if (j_in_halo /= 0) then
                ncbeg = chalo%i_h2d(icad+j_in_halo)
                if (ncbeg /= 0 ) then ! multiplication of ndim x ndim blocks
                   ! if (PRESENT(debug))                            &
                   !      write (21+debug, *) 'Details2: ', j, nd2, &
                   !      (nabeg-1)/(nd1*nd3), (ncbeg-1)/(nd1*nd2), &
                   !      (nbbeg-1)/(nd2*nd3)
!DIR$ NOPATTERN
                   do n2 = 1, nd2
                      nbaddr = nbbeg + nd3 * (n2 - 1)
                      ncaddr = ncbeg + nd1 * (n2 - 1)
                      do n1 = 1, nd1
                         naaddr = nabeg + nd3 * (n1 - 1)
                         do n3 = 1, nd3
                            c(ncaddr+n1-1) = c(ncaddr+n1-1) + &
                                             a(naaddr+n3-1) * b(nbaddr+n3-1)
                         end do
                      end do
                   end do
                end if
             end if ! End of if(j_in_halo.ne.0)
!$omp end task
          end do ! End of j = 1, nbnab
       end do ! End of i = 1, at%n_hnab
!$omp taskwait
    end do ! End of k = 1, nahpart
    ! need to wait for all the computation tasks to finish before the
    ! generating task is to end.
!$omp end task
!$omp end single
!$omp end parallel

    deallocate(nd1_1st, nd2_1st, jbnab2ch)

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
  !!  USES
  !!
  !!  AUTHOR
  !!   M.J.Gillan/D.R.Bowler
  !!  CREATION DATE
  !!   25/11/99
  !!  MODIFICATION HISTORY
  !!   20/06/2001 dave
  !!    Added ROBODoc header
  !!   2012/07/18 L.Tong
  !!   - Added OpenMP directives
  !!   2012/08/02 L.Tong
  !!   - Added C preprocessors for switching between codes for openMP
  !!     and standard, this is done to get the best efficiency. See
  !!     the notes in m_kern_max for details.
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

    implicit none

    ! Passed variables
    type(matrix_halo)  :: ahalo, chalo
    type(matrix_trans) :: at
    integer            :: mx_absb, mx_part, mx_iprim, lena, lenb, lenc
    integer            :: kpart, k_off
    ! Remember that a is a local transpose
    real(double)       :: a(lena)
    real(double)       :: b(lenb)
    real(double)       :: c(lenc)
    ! dimension declarations
    integer :: ibaddr(mx_part)
    integer :: ib_nd_acc(mx_part)
    integer :: nbnab(mx_part)
    integer :: ibpart(mx_part*mx_absb)
    integer :: ibseq(mx_part*mx_absb)
    integer :: bndim2(mx_part*mx_absb)
    ! Local variables
    integer :: k, k_in_part, k_in_halo, nbkbeg, j, jpart, jseq
    integer :: i, nabeg, i_in_prim, icad, nbbeg, j_in_halo, ncbeg
    integer :: n1, n2, n3, nb_nd_kbeg
    integer :: nd1, nd2, nd3
    integer :: naaddr, nbaddr, ncaddr
    ! The OpenMP firstprivate clause for some reason crashes (at least
    ! for gcc) whenever there is automatic arrays in its
    ! list. Everthing works fine with allocatable arrays. The
    ! allocation of allocatable arrays may be slightly slower than
    ! automatic arrays, but it is the only option at the moment.
    integer, dimension(:), allocatable :: nd1_1st, nd2_1st, jbnab2ch

    allocate(nd1_1st(at%mx_halo), nd2_1st(mx_absb), jbnab2ch(mx_absb))

!$omp parallel default(none)                                             &
!$omp          shared(kpart, ibaddr, ib_nd_acc, nbnab, ibpart, ibseq,    &
!$omp                 k_off, bndim2, mx_absb, mx_part, at, ahalo, chalo, &
!$omp                 a, b, c)                                           &
!$omp          private(i, j, k, j_in_halo, k_in_halo, k_in_part, nbkbeg, &
!$omp                  nb_nd_kbeg, nd1, nd2, nd3, jpart, jseq, jbnab2ch, &
!$omp                  nabeg, nbbeg, ncbeg, i_in_prim, icad, naaddr,     &
!$omp                  nbaddr, ncaddr, n1, n2, n3, nd1_1st, nd2_1st)
    ! run on single thread for generating tasks
!$omp single
    ! task for generating computation tasks
!$omp task untied                                                    &
!$omp      default(none)                                             &
!$omp      shared(kpart, ibaddr, ib_nd_acc, nbnab, ibpart, ibseq,    &
!$omp             k_off, bndim2, mx_absb, mx_part, at, ahalo, chalo, &
!$omp             a, b, c)                                           &
!$omp      private(i, j, k, j_in_halo, k_in_halo, k_in_part, nbkbeg, &
!$omp              nb_nd_kbeg, nd1, nd2, nd3, jpart, jseq, jbnab2ch, &
!$omp              nabeg, nbbeg, ncbeg, i_in_prim, icad, naaddr,     &
!$omp              nbaddr, ncaddr, n1, n2, n3, nd1_1st, nd2_1st)
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
       ! Loop over primary-set A-neighbours of k
       do i = 1, at%n_hnab(k_in_halo)
          i_in_prim = at%i_prim(at%i_beg(k_in_halo)+i-1)
          nd1 = ahalo%ndimi(i_in_prim)
          icad = (i_in_prim-1) * chalo%ni_in_halo
          nabeg = at%i_nd_beg(k_in_halo) + nd1_1st(i)
!$omp task default(none)                                        &
!$omp      shared(nbnab, bndim2, chalo, a, b, c)                &
!$omp      firstprivate(nd1, nd3, nbkbeg, nb_nd_kbeg, nd2_1st,  &
!$omp                   jbnab2ch, nabeg, icad, k_in_part)       &
!$omp      private(j, nd2, nbbeg, j_in_halo, ncbeg, n1, n2, n3, &
!$omp              naaddr, nbaddr, ncaddr)
          ! Loop over B-neighbours of atom k
          do j = 1, nbnab(k_in_part)
             nd2 = bndim2(nbkbeg+j-1)
             nbbeg = nb_nd_kbeg + nd2_1st(j)
             j_in_halo = jbnab2ch(j)
             if (j_in_halo /= 0) then
                ncbeg = chalo%i_h2d(icad+j_in_halo)
                if (ncbeg /= 0) then  ! multiplication of ndim x ndim blocks
!DIR$ NOPATTERN
                   do n2=1, nd2
                      nbaddr = nbbeg + nd3 * (n2 - 1)
                      ncaddr = ncbeg + nd1 * (n2 - 1)
                      do n1 = 1, nd1
                         naaddr = nabeg + nd3 * (n1 - 1)
                         do n3 = 1, nd3
                            a(naaddr+n3-1) = a(naaddr+n3-1) + &
                                             c(ncaddr+n1-1) * b(nbaddr+n3-1)
                         end do
                      end do
                   end do
                end if
             end if
          end do
!$omp end task
       end do
    end do
    ! need to wait for all the computation tasks to finish before the
    ! generating task is to end.
!$omp taskwait
!$omp end task
!$omp end single
!$omp end parallel

    deallocate(nd1_1st, nd2_1st, jbnab2ch)

    return
  end subroutine m_kern_min
  !!*****

end module multiply_kernel
