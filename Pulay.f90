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
!!  SOURCE
!!
module Pulay

  use global_module, ONLY: io_lun

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
  subroutine DoPulay_nospin(Aij, al, n, max_n, inode, ionode)

    ! Module usage
    use datatypes
    use GenBlas, ONLY: sytrf, sytri

    implicit none

    ! Passed variables
    integer :: n, max_n,len,inode,ionode

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
!!  SOURCE
!!
  subroutine DoPulay2D_nospin(Aij, al, n, max_n, inode, ionode)

    ! Module usage
    use datatypes
    use numbers
    use GenBlas, ONLY: sytrf, sytri

    implicit none

    ! Passed variables
    integer :: n, max_n ,inode, ionode

    real(double) :: Aij(max_n, max_n),al(max_n)

    ! Local variables
    integer i,j, info, ipiv(max_n)
    real(double) :: Rtot,work(max_n), sum, OA(max_n,max_n)
    logical :: WRITE_OUT = .false.

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
       if(inode==ionode.AND.info/=0) write(io_lun,*) 'info is ',info
       if(inode==ionode.AND.info/=0) write(io_lun,*) ' sytri: info is ',info
       ! Only a triangle is returned, so complete Aij
       do i=1,n
          do j=i,n
             Aij(j,i) = Aij(i,j)
          enddo
       enddo
    end if
    ! Create the alphas
    al = zero
    sum = zero
    do i=1,n
       do j=1,n
          al(i) = al(i) + Aij(j,i)
          sum = sum+ Aij(j,i)
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
!!                                     \sum_j Aij_spin^{-1} N_spin_j
!!       alpha_spin_i = N * ---------------------------------------------------
!!                           \sum_spin \sum_ij Aij_spin^{-1} N_spin_i N_spin_j
!!
!!   where N_spin_i is the history of spin populations from i-th step
!!   this is entered into the subroutine using the optional parameters
!!   N_up and N_dn.
!! INPUTS
!! USES
!! AUTHOR
!!   L.Tong
!! CREATION DATE
!!   2011/09/18
!! MODIFICATION HISTORY
!! SOURCE
!!
  subroutine DoPulay_spin (Aij_up, Aij_dn, alpha_up, alpha_dn, n, &
       max_n, inode, ionode, N_up, N_dn)
    use datatypes
    use numbers
    use GenBlas, only: sytrf, sytri
    use global_module, only: flag_fix_spin_population, ne_in_cell
    use GenComms, only: cq_abort

    implicit none

    ! passed variables
    integer :: n, max_n, inode, ionode
    real (double) :: Aij_up(max_n, max_n), Aij_dn(max_n, max_n),&
         & alpha_up(max_n), alpha_dn(max_n)
    real (double), intent (in), optional :: N_up(max_n), N_dn(max_n)
    
    ! local variables
    integer ii, jj, info, ipiv(max_n)
    real(double) :: Rtot, work(max_n), sum, sum_up, sum_dn, OA(max_n, max_n)
    logical :: write_out = .true.

    ! Check if the electron number histories are present
    if (.not. flag_fix_spin_population) then
       if ((.not. present (N_up)) .or. (.not. present (N_dn)))&
            & call cq_abort ('DoPulay: N_up and N_dn not&
            & present, required for evaluation of Pulay coeeficients&
            & if spin populations are allowed to vary.')
    end if

    ! Invert A_spin
    if (n == 1) then
       Aij_up(n,n) = one / Aij_up(n,n)
       Aij_dn(n,n) = one / Aij_dn(n,n)
    else
       ! invert A_up
       info = 0
       ipiv = 0
       work = zero
       OA = Aij_up
       call sytrf ('U', n, Aij_up, max_n, ipiv, work, max_n, info)
       if (inode == ionode .and. info /= 0) then
          write (io_lun, *) inode, ' A_up: sytrf: info is ', info
          write (io_lun, *) inode, ' A_up was: ', OA
          write (io_lun, *) inode, ' A_up is: ', Aij_up
       end if
       info = 0
       call sytri ('U', n, Aij_up, max_n, ipiv, work, info)
       if (inode == ionode .and. info /= 0) write (io_lun, *) 'A_up:&
            & sytri info is ', info
       ! invert A_dn
       info = 0
       ipiv = 0
       work = zero
       OA = Aij_dn
       call sytrf ('U', n, Aij_dn, max_n, ipiv, work, max_n, info)
       if (inode == ionode .and. info /= 0) then
          write (io_lun, *) inode, ' A_dn: sytrf: info is ', info
          write (io_lun, *) inode, ' A_dn was: ', OA
          write (io_lun, *) inode, ' A_dn is: ', Aij_dn
       end if
       info = 0
       call sytri ('U', n, Aij_dn, max_n, ipiv, work, info)
       if (inode == ionode .and. info /= 0) write (io_lun, *) 'A_dn:&
            & sytri info is ', info
       ! only a triangle is returned, so complete Aij_spin
       do ii = 1, n
          do jj = ii, n
             Aij_up(jj,ii) = Aij_up(ii,jj)
             Aij_dn(jj,ii) = Aij_dn(ii,jj)
          end do
       end do
    end if
    ! create the alpha_spins 
    alpha_up = zero
    alpha_dn = zero
    if (flag_fix_spin_population) then
       sum_up = zero
       sum_dn = zero
       do ii = 1, n
          do jj = 1, n
             alpha_up(ii) = alpha_up(ii) + Aij_up(ii,jj)
             sum_up = sum_up + Aij_up(ii,jj)
             alpha_dn(ii) = alpha_dn(ii) + Aij_dn(ii,jj)
             sum_dn = sum_dn + Aij_dn(ii,jj)
          end do
       end do
       do ii = 1, n
          alpha_up(ii) = alpha_up(ii) / sum_up
          alpha_dn(ii) = alpha_dn(ii) / sum_dn
       end do
    else ! if (flag_fix_spin_population)
       sum = zero
       do ii = 1, n
          do jj = 1, n
             alpha_up(ii) = alpha_up(ii) + Aij_up(ii,jj) * N_up(jj)
             alpha_dn(ii) = alpha_dn(ii) + Aij_dn(ii,jj) * N_dn(jj)
             sum = sum + Aij_up(ii,jj) * N_up(ii) * N_up(jj)
             sum = sum + Aij_dn(ii,jj) * N_dn(ii) * N_dn(jj)
          end do
       end do
       do ii = 1, n
          alpha_up(ii) = alpha_up(ii) * ne_in_cell / sum
          alpha_dn(ii) = alpha_dn(ii) * ne_in_cell / sum
       end do
    end if ! if (flag_fix_spin_population)
    ! write out alphas
    if (inode .eq. ionode .and. write_out) then 
       write (io_lun, *) 'Alpha for spin up:'
       if (n .eq. 2) write (io_lun, 211) (alpha_up(ii), ii = 1, n)
       if (n .eq. 3) write (io_lun, 212) (alpha_up(ii), ii = 1, n)
       if (n .eq. 4) write (io_lun, 213) (alpha_up(ii), ii = 1, n)
       write (io_lun, *) 'Alpha for spin down:'
       if (n .eq. 2) write (io_lun, 211) (alpha_dn(ii), ii = 1, n)
       if (n .eq. 3) write (io_lun, 212) (alpha_dn(ii), ii = 1, n)
       if (n .eq. 4) write (io_lun, 213) (alpha_dn(ii), ii = 1, n)
    endif
211 format(2f10.6)
212 format(3f10.6)
213 format(4f10.6)
    return
  end subroutine DoPulay_spin
!!****

end module Pulay
