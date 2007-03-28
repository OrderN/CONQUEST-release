! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: comm_array_module.f90,v 1.2 2003/06/11 09:01:54 drb Exp $
! ------------------------------------------------------------------------------
! Module comm_array_module
! ------------------------------------------------------------------------------
! Code area 8: indexing
! ------------------------------------------------------------------------------

!!****h* Conquest/comm_array_module *
!!  NAME
!!   comm_array_module
!!  PURPOSE
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   2000/1
!!  MODIFICATION HISTORY
!!   10:01, 2003/06/11 dave
!!    Added headers
!!  SOURCE
!!
module comm_array_module

  use datatypes
  use numbers, ONLY: zero
  !use maxima_module, ONLY : mx_send_real, mx_recv_real, mx_send_int, mx_recv_int
!  use maxima_module, ONLY : mx_send_int, mx_recv_int
  use GenComms, ONLY: cq_abort

  implicit none

  ! 23/Mar/2001 T. Miyazaki
  ! integer,parameter :: mx_send_real=48640
  ! integer,parameter :: mx_recv_real=29440
  ! integer,parameter :: mx_send_int =6100
  ! integer,parameter :: mx_recv_int =2440

  save
  real(double), allocatable, dimension(:), target :: send_array
  real(double), allocatable, dimension(:), target :: recv_array
  integer, allocatable, dimension(:), target      :: isend_array
  integer, allocatable, dimension(:), target      :: irecv_array

  ! interface check_commarray_real
  !  module procedure check_commarray_real
  ! end interface check_commarray_real
  ! interface check_commarray_int
  !  module procedure check_commarray_int
  ! end interface check_commarray_int

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id: comm_array_module.f90,v 1.2 2003/06/11 09:01:54 drb Exp $"
!!***

!%%!contains

!%%! !!****f* comm_array_module/check_commarray_real *
!%%! !!
!%%! !!  NAME 
!%%! !!   check_commarray_real
!%%! !!  USAGE
!%%! !! 
!%%! !!  PURPOSE
!%%! !!   check size of send_array and recv_array
!%%! !! 
!%%! !!  INPUTS
!%%! !! 
!%%! !! 
!%%! !!  USES
!%%! !! 
!%%! !!  AUTHOR
!%%! !!   T.Miyazaki
!%%! !!  CREATION DATE
!%%! !!   2000/1
!%%! !!  MODIFICATION HISTORY
!%%! !!   09:45, 29/10/2002 dave
!%%! !!    Added ROBODoc headers
!%%! !!  SOURCE
!%%! !!
!%%!   subroutine check_commarray_real(name,myid,mx_send_in,mx_recv_in)
!%%! 
!%%!     use datatypes
!%%! 
!%%!     implicit none
!%%! 
!%%!     character(len=*),intent(in) :: name !Name of the calling subroutine
!%%!     integer,intent(in) :: myid,mx_send_in,mx_recv_in
!%%! 
!%%!     integer :: mynode,ierr=0
!%%!     integer :: irc
!%%! 
!%%!     mynode=myid+1
!%%!     write(*,*) ' In check_commarray_real called from ',name, &
!%%!          '  node no. = ',mynode,  &
!%%!          '  mx_send_in,mx_send_real = ',mx_send_in,mx_send_real, &
!%%!          '  mx_recv_in,mx_recv_real = ',mx_recv_in,mx_recv_real
!%%! 
!%%!     if(mx_send_in > mx_send_real) then
!%%!        write(*,*) ' ERROR in check_commarray_real called from',name
!%%!        write(*,*) '  node no. = ',mynode,'  mx_send_in,mx_send_real = ' &
!%%!             ,mx_send_in,mx_send_real
!%%!        ierr=1
!%%!     endif
!%%!     if(mx_recv_in > mx_recv_real) then
!%%!        write(*,*) ' ERROR in check_commarray_real called from',name
!%%!        write(*,*) '  node no. = ',mynode,'  mx_recv_in,mx_recv_real = ' &
!%%!             ,mx_recv_in,mx_recv_real
!%%!        ierr=ierr+2
!%%!     endif
!%%! 
!%%!     if(ierr /= 0) call cq_abort(' in check_commarray_real :ierr= ',ierr)
!%%!     send_array(1:mx_send_in)=zero
!%%!     recv_array(1:mx_recv_in)=zero
!%%! 
!%%!     return
!%%!   end subroutine check_commarray_real
!%%! !!***

!%%! !!****f* comm_array_module/check_commarray_int *
!%%! !!
!%%! !!  NAME 
!%%! !!   check_commarray_int
!%%! !!  USAGE
!%%! !! 
!%%! !!  PURPOSE
!%%! !!   check size of isend_array and irecv_array
!%%! !! 
!%%! !!  INPUTS
!%%! !! 
!%%! !! 
!%%! !!  USES
!%%! !! 
!%%! !!  AUTHOR
!%%! !!   T.Miyazaki
!%%! !!  CREATION DATE
!%%! !!   2000/1
!%%! !!  MODIFICATION HISTORY
!%%! !!
!%%! !!  SOURCE
!%%! !!
!%%!   subroutine check_commarray_int(name,myid,mx_send_in,mx_recv_in)
!%%! 
!%%!     use datatypes
!%%! 
!%%!     implicit none
!%%! 
!%%!     character(len=*),intent(in) :: name !Name of the calling subroutine
!%%!     integer,intent(in) :: myid,mx_send_in,mx_recv_in
!%%!     integer :: mynode,ierr=0
!%%!     integer :: irc
!%%! 
!%%!     mynode=myid+1
!%%!     write(*,*) ' In check_commarray_int called from ',name, &
!%%!          '  node no. = ',mynode,  &
!%%!          '  mx_send_in,mx_send_real = ',mx_send_in,mx_send_int, &
!%%!          '  mx_recv_in,mx_recv_real = ',mx_recv_in,mx_recv_int
!%%! 
!%%!     if(mx_send_in > mx_send_int) then
!%%!        write(*,*) ' ERROR in check_commarray_int called from',name
!%%!        write(*,*) '  node no. = ',mynode,'  mx_send_in,mx_send_int = ' &
!%%!             ,mx_send_in,mx_send_int
!%%!        ierr=1
!%%!     endif
!%%!     if(mx_recv_in > mx_recv_int) then
!%%!        write(*,*) ' ERROR in check_commarray_int called from',name
!%%!        write(*,*) '  node no. = ',mynode,'  mx_recv_in,mx_recv_int = ' &
!%%!             ,mx_recv_in,mx_recv_int
!%%!        ierr=ierr+2
!%%!     endif
!%%! 
!%%!     if(ierr /= 0) call cq_abort(' in check_commarray_int :ierr= ',ierr)
!%%!     isend_array(1:mx_send_in)=0
!%%!     irecv_array(1:mx_recv_in)=0
!%%! 
!%%!     return
!%%!   end subroutine check_commarray_int
!%%! !!***
end module comm_array_module
