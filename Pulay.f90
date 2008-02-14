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
!!  SOURCE
!!
module Pulay

  use global_module, ONLY: io_lun

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
!!***

contains

!!****f* Pulay/DoPulay *
!!
!!  NAME 
!!   DoPulay
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
!!  SOURCE
!!
  subroutine DoPulay(Aij,al,n,max_n,inode,ionode)

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
  end subroutine DoPulay
!!***

!!****f* Pulay/DoPulay2D *
!!
!!  NAME 
!!   DoPulay2D
!!  USAGE
!! 
!!  PURPOSE
!!   As for DoPulay, but allows the A matrix to be a proper
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
!!  SOURCE
!!
  subroutine DoPulay2D(Aij,al,n,max_n,inode,ionode)

    ! Module usage
    use datatypes
    use numbers
    use GenBlas, ONLY: sytrf, sytri

    implicit none

    ! Passed variables
    integer :: n, max_n,len,inode,ionode

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
       if(info/=0) write(io_lun,*) 'info is ',info
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
  end subroutine DoPulay2D
!!***
end module Pulay

