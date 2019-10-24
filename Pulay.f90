! $Id$
! -----------------------------------------------------------
! Module Pulay
! -----------------------------------------------------------
! Code area 9: general
! -----------------------------------------------------------

!!****h* Conquest/Pulay
!!  NAME
!!   Pulay
!!  PURPOSE
!!   Collects together routines used in doing Pulay-type
!!   minimisation.  At present only finds the values for
!!   mixing the different iterations
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!
!!  MODIFICATION HISTORY
!!   22/05/2001 dave
!!    ROBODoc header
!!   11/06/2001 dave
!!    Added RCS Id and Log tags and changed dsytr LAPACK calls
!!    to generic ones through GenBlas
!!   2008/02/01 17:49 dave
!!    Changes for output to file not stdout
!!   2012/03/01 L.Tong
!!    Added interface DoPulay
!!   2019/10/24 10:30 dave
!!    Added inode and ionode as use-associated from GenComms in module header
!!    for efficiency and best practice; changed function calls
!!  SOURCE
!!
module Pulay

  use global_module, ONLY: io_lun
  use GenComms, ONLY: inode, ionode

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
!!***

  !!****f* Pulay/DoPulay
  !! PURPOSE
  !!   interface for DoPulay_nospin and DoPulay_spin
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/01
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  interface DoPulay
     module procedure DoPulay2D_nospin  ! use the 2D version
     module procedure DoPulay_spin
     ! module procedure DoPulay_spinfixed
     ! module procedure DoPulay_spinrelax
  end interface DoPulay
  !!*****

contains

!!****f* Pulay/DoPulay_nospin *
!!
!!  NAME
!!   DoPulay_nospin
!!  USAGE
!!
!!  PURPOSE
!!   Given the matrix A, inverts and sums to create the
!!   coefficients alpha:
!!    alpha(i) = \sum_j A_ij^-1/\sum_jk A_jk^-1
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
!!   22/05/2001 dave
!!    ROBODoc header
!!   11/06/2001 dave
!!    Generic LAPACK calls
!!   2012/03/01 L.Tong
!!    Renamed to DoPulay_nospin
!!  SOURCE
!!
  subroutine DoPulay_nospin(Aij, al, n, max_n)

    ! Module usage
    use datatypes
    use GenBlas, ONLY: sytrf, sytri

    implicit none

    ! Passed variables
    integer :: n, max_n,len

    real(double) :: Aij(:),al(:)

    ! Local variables
    integer i,j, info, ipiv(max_n)
    real(double) :: Rtot,work(max_n), sum
    logical :: WRITE_OUT = .false.

    ! Start the Pulay procedure
    ! Invert A
    info = 0
    ipiv = 0
    work = 0.0_double
    call sytrf('U',n,Aij,n,ipiv,work,n,info)
    if(inode.eq.ionode.AND.info.ne.0) write(io_lun,*) 'info is ',info
    info = 0
    call sytri('U',n,Aij,n,ipiv,work,info)
    if(inode.eq.ionode.AND.info.ne.0) write(io_lun,*) 'info is ',info
    ! Only a triangle is returned, so complete Aij
    do i=1,n
       do j=i,n
          Aij(j+(i-1)*n) = Aij(i+(j-1)*n)
       enddo
    enddo
    ! Create the alphas
    al = 0.0_double
    sum = 0.0_double
    do i=1,n
       do j=1,n
          al(i) = al(i) + Aij(j+(i-1)*n)
          sum = sum+ Aij(j+(i-1)*n)
       enddo
    enddo
    do i=1,n
       al(i) = al(i)/sum
    enddo
    ! Write out alphas
    if(inode.eq.ionode.AND.WRITE_OUT) then
       write(io_lun,*) 'Alpha:'
       if(n.eq.2) write(io_lun,111) (al(i),i=1,n)
       if(n.eq.3) write(io_lun,112) (al(i),i=1,n)
       if(n.eq.4) write(io_lun,113) (al(i),i=1,n)
    endif
111 format(2f10.6)
112 format(3f10.6)
113 format(4f10.6)
    return
  end subroutine DoPulay_nospin
!!***

!!****f* Pulay/DoPulay2D_nospin *
!!
!!  NAME
!!   DoPulay2D_nospin
!!  USAGE
!!
!!  PURPOSE
!!   As for DoPulay_nospin, but allows the A matrix to be a proper
!!   2D array - in other words, distinguishes between a maximum
!!   size and a current size
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
!!   22/05/2001 dave
!!    Added ROBODoc header
!!   11/06/2001 dave
!!    Changed LAPACK calls dsytr to generic ones through GenBlas
!!   2011/09/18 L.Tong
!!    Removed a redundant (wrong?) reference to len in list of
!!    passed variables. len is not used in the inputs, nor in the code
!!   2012/03/01 L.Tong
!!    Renamed to DoPulay2D_nospin
!!   2012/05/26 L.Tong
!!   - added check in the case of Aij non-invertable, and we need
!!     iPulay (the current pulay history slot) to set alphas to
!!     correspond to linear mixing
!!  SOURCE
!!
  subroutine DoPulay2D_nospin(iPulay, Aij, al, n, max_n)

    ! Module usage
    use datatypes
    use numbers
    use GenBlas, ONLY: sytrf, sytri

    implicit none

    ! Passed variables
    integer :: iPulay, n, max_n

    real(double) :: Aij(max_n, max_n),al(max_n)

    ! Local variables
    integer i,j, info, ipiv(max_n)
    real(double) :: Rtot,work(max_n), sum, OA(max_n,max_n)
    logical :: WRITE_OUT = .false.
    logical :: Aij_invertable


    ! Start the Pulay procedure
    ! Invert A
    if(n==1) then
       Aij(n,n) = one/Aij(n,n)
    else
       info = 0
       ipiv = 0
       work = 0.0_double
       OA = Aij
       call sytrf('U',n,Aij,max_n,ipiv,work,max_n,info)
       if(inode==ionode.AND.info/=0) then
          write(io_lun,*) inode,' sytrf: info is ',info
          write(io_lun,*) inode,' A was: ',OA
          write(io_lun,*) inode,' A is: ',Aij
       end if
       info = 0
       call sytri('U',n,Aij,max_n,ipiv,work,info)
       if(inode==ionode.AND.info/=0) then
          write(io_lun,*) ' sytri: info is ',info
          Aij_invertable = .false.
       else
          Aij_invertable = .true.
          ! Only a triangle is returned, so complete Aij
          do i=1,n
             do j=i,n
                Aij(j,i) = Aij(i,j)
             enddo
          enddo
       end if
    end if
    ! Create the alphas
    if (Aij_invertable) then
       al = zero
       sum = zero
       do i=1,n
          do j=1,n
             al(i) = al(i) + Aij(j,i)
             sum = sum + Aij(j,i)
          enddo
       enddo
       do i=1,n
          al(i) = al(i)/sum
       enddo
    else
       al = zero
       al(iPulay) = one
    end if

    ! Write out alphas
    if (inode == ionode .and. write_out) then
       write (io_lun, '(1x,a)') 'Alpha are:'
       write (io_lun, '(1x,f12.6)', advance='no') (al(i), i=1,n)
       write (io_lun, '(/)')
    end if

    return
  end subroutine DoPulay2D_nospin
!!***


  !!****f* Pulay/DoPulay_spin
  !! PURPOSE
  !!   Calculates pulay cofficients for spin polarised calculations.
  !!   If the spin populations are fixed then
  !!
  !!                       \sum_j Aij_spin^{-1}
  !!       alpha_spin_i = -----------------------
  !!                       \sum_ij Aij_spin^{-1}
  !!
  !!   If the spin populations are allowed to vary then
  !!
  !!                         \sum_j \sum_spin Aij_spin^{-1}
  !!       alpha_spin_i = ------------------------------------
  !!                         \sum_ij \sum_spin  Aij_spin^{-1}
  !!
  !! INPUTS
  !! USES
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2011/09/18
  !! MODIFICATION HISTORY
  !!   2012/03/19 L.Tong
  !!   - Changed spin implementation
  !!   2012/05/26 L.Tong
  !!   - added check in the case of Aij non-invertable, and we need
  !!     iPulay (the current pulay history slot) to set alphas to
  !!     correspond to linear mixing
  !! SOURCE
  !!
  subroutine DoPulay_spin(iPulay, Aij, alpha, n, max_n)
    use datatypes
    use numbers
    use GenBlas,       only: sytrf, sytri
    use global_module, only: flag_fix_spin_population, ne_in_cell, nspin
    use GenComms,      only: cq_abort

    implicit none

    ! passed variables
    integer :: iPulay, n, max_n
    real(double), dimension(max_n,max_n,nspin), intent(inout) :: Aij
    real(double), dimension(max_n,nspin),       intent(out)   :: alpha
    ! local variables
    integer      :: ii, jj, info, spin, max_spin
    logical      :: write_out = .false.
    logical      :: Aij_invertable
    real(double) :: denominator
    integer,      dimension(max_n)       :: ipiv
    real(double), dimension(max_n)       :: work
    real(double), dimension(max_n,max_n) :: old_A

    ! store sum of spin parts in Aij if doing spin relaxation
    if (nspin == 2 .and. (.not. flag_fix_spin_population)) then
       do ii = 1, n
          do jj = 1, n
             Aij(ii,jj,1) = Aij(ii,jj,1) + Aij(ii,jj,2)
          end do
       end do
       max_spin = 1
    else
       max_spin = nspin
    end if

    ! Invert A
    do spin = 1, max_spin
       if (n == 1) then
          Aij(n,n,spin) = one / Aij(n,n,spin)
       else
          ! invert A_spin
          info = 0
          ipiv = 0
          work = zero
          old_A(:,:) = Aij(:,:,spin)
          call sytrf('U', n, Aij(:,:,spin), max_n, ipiv, work, max_n, info)
          if (inode == ionode .and. info /= 0) then
             write (io_lun, '(a,i5,a,i2,a,i5)') &
                   'Node:', inode, ' A(spin=', spin, '): sytrf: info is ', info
             write (io_lun, '(a,i5,a,i2,a)') &
                   'Node:', inode, ' A(spin=', spin, ') was: '
             write (io_lun, *) old_A
             write (io_lun, '(a,i5,a,i2,a)') &
                   'Node:', inode, ' A(spin=', spin, ') is: '
             write (io_lun, *) Aij(:,:,spin)
          end if
          info = 0
          call sytri('U', n, Aij(:,:,spin), max_n, ipiv, work, info)
          if (inode == ionode .and. info /= 0) then
             write (io_lun, '(a,i5,a,i2,a,i5)') &
                   'Node:', inode, 'A(spin=', spin, '): sytri: info is ', info
             Aij_invertable = .false.
          else
             Aij_invertable = .true.
             ! only a triangle is returned, so complete Aij_spin
             do ii = 1, n
                do jj = ii, n
                   Aij(jj,ii,spin) = Aij(ii,jj,spin)
                end do
             end do
          end if ! check if inversion successful
       end if ! (n == 1)
    end do ! spin

    ! create the alpha_spins
    if (Aij_invertable) then
       alpha = zero
       do spin = 1, max_spin
          denominator = zero
          do ii = 1, n
             do jj = 1, n
                ! note Aij_tot is symmetric, swap ii and jj for efficiency
                alpha(ii,spin) = alpha(ii,spin) + Aij(jj,ii,spin)
                denominator = denominator + Aij(jj,ii,spin)
             end do ! jj
          end do ! ii
          do ii = 1, n
             alpha(ii,spin) = alpha(ii,spin) / denominator
          end do
       end do ! spin
    else ! Aij is non-invertable
       ! in this case, just do linear mixing, i.e. set alpha(n,spin) = one
       ! the rest to zero
       do spin = 1, max_spin
          alpha(:,spin) = zero
          alpha(iPulay,spin) = one
       end do
    end if

    if (nspin == 2 .and. .not. flag_fix_spin_population) then
       ! copy alpha for the second spin component
       alpha(1:n,2) = alpha(1:n,1)
    end if

    ! write out alphas
    if (inode == ionode .and. write_out) then
       do spin = 1, nspin
          write (io_lun, '(1x,a,i2,a)') &
               'Alpha for spin = ', spin, ' are:'
          write (io_lun, '(1x,f12.6)', advance='no') (alpha(ii,spin), ii=1,n)
          write (io_lun, '(/)')
       end do
    end if
  end subroutine DoPulay_spin
  !!****

end module Pulay

