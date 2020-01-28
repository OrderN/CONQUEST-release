! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module GenBlas
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/GenBlas
!!  NAME
!!   GenBlas
!!  PURPOSE
!!   Creates a generic interface to the BLAS routines so that different
!!   shape arrays can be easily passed to the same routine through an
!!   interface, without F90 compilers whinging about problems
!!  USES
!!   datatypes
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   30/06/00
!!  MODIFICATION HISTORY
!!   Frequent to add extra functionality
!!   30/05/2001 dave
!!    ROBODoc header
!!   11/06/2001 dave
!!    Added LAPACK calls sytrf, sytri, potrf, potri and RCS Id and Log tags
!!   17/06/2002 dave
!!    Added onezdot - dot product for 1D complex arrays, and extended copy to complex. Also improved headers a little
!!    and added static RCSid object
!!   15:47, 04/02/2003 drb 
!!    Added vdot calls and many comments
!!   2008/02/06 08:16 dave
!!    Changed for output to file not stdout
!!  SOURCE
!!
module GenBlas

  use datatypes
  use numbers
  use global_module, ONLY: io_lun

  implicit none

  interface dot
     module procedure onedot
     module procedure twodot
     module procedure threedot
     module procedure onezdot
  end interface
!!***

!!****f* GenBlas/vdot *
!!
!!  NAME 
!!   vdot
!!  USAGE
!!   vdot(len,x,stridex,y,stridey)
!!  PURPOSE
!!   Provides an interface to dot products plus sum across processors for arrays of different dimensionality
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   30/06/00
!!  MODIFICATION HISTORY
!!   25/04/2002 dave
!!    Added ROBODoc header
!!  SOURCE
!!
  interface vdot
     module procedure v1dot
     module procedure v2dot
     module procedure v3dot
  end interface
!!***
  interface scal
     module procedure onescal
     module procedure twoscal
     module procedure threescal
  end interface

  interface copy
     module procedure onecopy
     module procedure twocopy
     module procedure threecopy
     module procedure onezcopy
     module procedure twozcopy
  end interface

!!****f* GenBlas/axpy *
!!
!!  NAME 
!!   axpy
!!  USAGE
!!   axpy(len,x,stridex,y,stridey)
!!  PURPOSE
!!   Provides an interface to vector sum for arrays of different dimensionality
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   30/06/00
!!  MODIFICATION HISTORY
!!   25/04/2002 dave
!!    Added ROBODoc header
!!  SOURCE
!!
  interface axpy
     module procedure one_axpy
     module procedure two_axpy
     module procedure three_axpy
  end interface
!!***

  interface gemm
     module procedure one_gemm
     module procedure onetwo_gemm
     module procedure two_gemm
  end interface

  interface sytrf
     module procedure one_sytrf
     module procedure two_sytrf
  end interface

  interface sytri
     module procedure one_sytri
     module procedure two_sytri
  end interface

  interface potrf
     module procedure two_potrf
  end interface

  interface potri
     module procedure two_potri
  end interface

contains

!!****f* GenBlas/one_axpy *
!!
!!  NAME 
!!   one_axpy
!!  PURPOSE
!!   Vector sum and multiply for rank 1 arrays
!!  INPUTS
!!   integer :: len - length of array
!!   integer :: stride_a, stride_b - stride of arrays
!!   real(double) :: alpha - scale factor for a
!!   real(double) :: a, b - arrays
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  subroutine one_axpy(len,alpha,a,stride_a,b,stride_b)

    implicit none

    real(double) :: a(:),b(:),alpha
    integer :: len,stride_a,stride_b

    call daxpy(len,alpha,a,stride_a,b,stride_b)
  end subroutine one_axpy
!!***

!!****f* GenBlas/two_axpy *
!!
!!  NAME 
!!   two_axpy
!!  PURPOSE
!!   Vector sum and multiply for rank 2 arrays
!!  INPUTS
!!   integer :: len - length of array
!!   integer :: stride_a, stride_b - stride of arrays
!!   real(double) :: alpha - scale factor for a
!!   real(double) :: a, b - arrays
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  subroutine two_axpy(len,alpha,a,stride_a,b,stride_b)

    implicit none

    real(double) :: a(:,:),b(:,:),alpha
    integer :: len,stride_a,stride_b

    call daxpy(len,alpha,a,stride_a,b,stride_b)
  end subroutine two_axpy
!!***

!!****f* GenBlas/three_axpy *
!!
!!  NAME 
!!   three_axpy
!!  PURPOSE
!!   Vector sum and multiply for rank 3 arrays
!!  INPUTS
!!   integer :: len - length of array
!!   integer :: stride_a, stride_b - stride of arrays
!!   real(double) :: alpha - scale factor for a
!!   real(double) :: a, b - arrays
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  subroutine three_axpy(len,alpha,a,stride_a,b,stride_b)

    implicit none

    real(double) :: a(:,:,:),b(:,:,:),alpha
    integer :: len,stride_a,stride_b

    call daxpy(len,alpha,a,stride_a,b,stride_b)
  end subroutine three_axpy
!!***

  subroutine onescal(len,fac,a,stride)

    implicit none

    integer :: len, stride
    real(double) :: a(:), fac

    call dscal(len,fac,a,stride)
  end subroutine onescal

  subroutine twoscal(len,fac,a,stride)

    implicit none

    real(double) :: a(:,:), fac
    integer :: len, stride

    call dscal(len,fac,a,stride)
  end subroutine twoscal

  subroutine threescal(len,fac,a,stride)

    implicit none

    real(double) :: a(:,:,:), fac
    integer :: len, stride

    call dscal(len,fac,a,stride)
  end subroutine threescal

  subroutine onecopy(len,a,stridea,b,strideb)

    implicit none

    real(double) :: a(:), b(:)
    integer :: len, stridea,strideb

    call dcopy(len,a,stridea,b,strideb)
  end subroutine onecopy

  subroutine twocopy(len,a,stridea,b,strideb)

    implicit none

    real(double) :: a(:,:), b(:,:)
    integer :: len, stridea,strideb

    call dcopy(len,a,stridea,b,strideb)
  end subroutine twocopy

  subroutine threecopy(len,a,stridea,b,strideb)

    implicit none

    real(double) :: a(:,:,:), b(:,:,:)
    integer :: len, stridea,strideb

    call dcopy(len,a,stridea,b,strideb)
  end subroutine threecopy

  subroutine onezcopy(len,a,stridea,b,strideb)

    implicit none

    complex(double_cplx) :: a(:), b(:)
    integer :: len, stridea,strideb

    call zcopy(len,a,stridea,b,strideb)
  end subroutine onezcopy

  subroutine twozcopy(len,a,stridea,b,strideb)

    implicit none

    complex(double_cplx) :: a(:,:), b(:,:)
    integer :: len, stridea,strideb

    call zcopy(len,a,stridea,b,strideb)
  end subroutine twozcopy

  subroutine one_gemm(ta,tb,lda,ldb,ldc,al,a,na,b,nb,beta,c,nc)

    implicit none

    character :: ta,tb
    real(double) :: al,a(:),b(:),c(:),beta
    integer :: lda,ldb,ldc,na,nb,nc

    call dgemm(ta,tb,lda,ldb,ldc,al,a,na,b,nb,beta,c,nc)
  end subroutine one_gemm

  subroutine onetwo_gemm(ta,tb,lda,ldb,ldc,al,a,na,b,nb,beta,c,nc)

    implicit none

    character :: ta,tb
    real(double) :: al,a(:),b(:),c(:,:),beta
    integer :: lda,ldb,ldc,na,nb,nc

    call dgemm(ta,tb,lda,ldb,ldc,al,a,na,b,nb,beta,c,nc)
  end subroutine onetwo_gemm

  subroutine two_gemm(ta,tb,lda,ldb,ldc,al,a,na,b,nb,beta,c,nc)

    implicit none

    character :: ta,tb
    real(double) :: al,a(:,:),b(:,:),c(:,:),beta
    integer :: lda,ldb,ldc,na,nb,nc

    call dgemm(ta,tb,lda,ldb,ldc,al,a,na,b,nb,beta,c,nc)
  end subroutine two_gemm

  real(double) function asum(len,a,stride_a)

    implicit none

    real(double) :: a(:)
    real(double),external :: dasum
    integer :: len,stride_a

    asum = dasum(len,a,stride_a)
  end function asum

  ! Take sum of relative (i.e. not absolute) values of a
  real(double) function rsum(len,a,stride_a)

    implicit none

    real(double) :: a(:)
    real(double),external :: dasum
    integer :: len,stride_a

    ! Local
    integer :: i

    rsum = zero
    do i=1,len,stride_a
       rsum = rsum + a(i)
    end do
  end function rsum

!!****f* GenBlas/onedot *
!!
!!  NAME 
!!   onedot
!!  PURPOSE
!!   Dot product for rank 1 real(double) arrays
!!  INPUTS
!!   integer :: len - length of array
!!   integer :: stride_a, stride_b - stride of arrays
!!   real(double) :: a, b - arrays
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  real(double) function onedot(len,a,stride_a,b,stride_b)

    implicit none

    real(double) :: a(:),b(:)
    real(double),external :: ddot
    integer :: len,stride_a,stride_b

    onedot = ddot(len,a,stride_a,b,stride_b)
  end function onedot
!!***

!!****f* GenBlas/twodot *
!!
!!  NAME 
!!   twodot
!!  PURPOSE
!!   Dot product for rank 2 real(double) arrays
!!  INPUTS
!!   integer :: len - length of array
!!   integer :: stride_a, stride_b - stride of arrays
!!   real(double) :: a, b - arrays
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  real(double) function twodot(len,a,stride_a,b,stride_b)

    implicit none

    real(double) :: a(:,:),b(:,:)
    real(double),external :: ddot
    integer :: len,stride_a,stride_b

    twodot = ddot(len,a,stride_a,b,stride_b)
  end function twodot
!!***

!!****f* GenBlas/threedot *
!!
!!  NAME 
!!   onedot
!!  PURPOSE
!!   Dot product for rank 3 real(double) arrays
!!  INPUTS
!!   integer :: len - length of array
!!   integer :: stride_a, stride_b - stride of arrays
!!   real(double) :: a, b - arrays
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  real(double) function threedot(len,a,stride_a,b,stride_b)

    implicit none

    real(double) :: a(:,:,:),b(:,:,:)
    real(double),external :: ddot
    integer :: len,stride_a,stride_b

    threedot = ddot(len,a,stride_a,b,stride_b)
  end function threedot
!!***

!!****f* GenBlas/onezdot *
!!
!!  NAME 
!!   onezdot
!!  PURPOSE
!!   Dot product for rank 1 complex(double_cplx) arrays
!!
!!   N.B. Uses zdotc which returns x*.y not x.y !
!!  INPUTS
!!   integer :: len - length of array
!!   integer :: stride_a, stride_b - stride of arrays
!!   complex(double_cplx) :: a, b - arrays
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  complex(double_cplx) function onezdot(len,a,stride_a,b,stride_b)

    implicit none

    complex(double_cplx) :: a(:),b(:)
    complex(double_cplx),external :: zdotc
    integer :: len,stride_a,stride_b

    onezdot = zdotc(len,a,stride_a,b,stride_b)
  end function onezdot
!!***

  subroutine one_sytrf(uplo,n,a,lda,ipiv,work,lwork,info)

    implicit none

    character :: uplo
    integer :: n, lda, lwork, info

    integer, dimension(n) :: ipiv
    real(double), dimension(lda*n) :: a
    real(double), dimension(lwork) :: work
    
    call dsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
  end subroutine one_sytrf

  subroutine two_sytrf(uplo,n,a,lda,ipiv,work,lwork,info)

    implicit none

    character :: uplo
    integer :: n, lda, lwork, info

    integer, dimension(n) :: ipiv
    real(double), dimension(lda,n) :: a
    real(double), dimension(lwork) :: work
    
    call dsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
  end subroutine two_sytrf

  subroutine one_sytri(uplo,n,a,lda,ipiv,work,info)

    implicit none

    character :: uplo
    integer :: n, lda, lwork, info

    integer, dimension(n) :: ipiv
    real(double), dimension(lda*n) :: a
    real(double), dimension(n) :: work
    
    call dsytri(uplo,n,a,lda,ipiv,work,info)
  end subroutine one_sytri

  subroutine two_sytri(uplo,n,a,lda,ipiv,work,info)

    implicit none

    character :: uplo
    integer :: n, lda, lwork, info

    integer, dimension(n) :: ipiv
    real(double), dimension(lda,n) :: a
    real(double), dimension(n) :: work
    
    call dsytri(uplo,n,a,lda,ipiv,work,info)
  end subroutine two_sytri

  subroutine two_potrf(uplo,n,a,lda,info)

    character :: uplo
    integer :: info, lda, n
    real(double) :: a(lda, n)

    call dpotrf(uplo,n,a,lda,info)
  end subroutine two_potrf

  subroutine two_potri(uplo,n,a,lda,info)

    character :: uplo
    integer :: info, lda, n
    real(double) :: a(lda, n)

    call dpotri(uplo,n,a,lda,info)
  end subroutine two_potri

!!****f* GenBlas/v1dot *
!!
!!  NAME 
!!   v1dot
!!  PURPOSE
!!   Dot product for rank 1 real(double) arrays with sum across all processors
!!  INPUTS
!!   integer :: len - length of array
!!   integer :: stride_a, stride_b - stride of arrays
!!   real(double) :: a, b - arrays
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  real(double) function v1dot(len,a,stride_a,b,stride_b)

    use mpi

    implicit none

    real(double) :: a(:),b(:)
    real(double) :: temp
    integer :: len,stride_a,stride_b
    integer :: ierr

    v1dot = 0.0_double
    temp = dot(len,a,stride_a,b,stride_b)
    call MPI_allreduce(temp, v1dot, 1, MPI_double_precision, &
         MPI_sum, MPI_COMM_WORLD, ierr)
    if(ierr/=0) write(io_lun,*) 'Problem with MPI_allreduce in v1dot !'
  end function v1dot
!!***

!!****f* GenBlas/v2dot *
!!
!!  NAME 
!!   v2dot
!!  PURPOSE
!!   Dot product for rank 2 real(double) arrays with sum across all processors
!!  INPUTS
!!   integer :: len - length of array
!!   integer :: stride_a, stride_b - stride of arrays
!!   real(double) :: a, b - arrays
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  real(double) function v2dot(len,a,stride_a,b,stride_b)

    use mpi

    implicit none

    real(double) :: a(:,:),b(:,:)
    real(double) :: temp
    integer :: len,stride_a,stride_b
    integer :: ierr

    temp = dot(len,a,stride_a,b,stride_b)
    call MPI_allreduce(temp, v2dot, 1, MPI_double_precision, &
         MPI_sum, MPI_COMM_WORLD, ierr)
    if(ierr/=0) write(io_lun,*) 'Problem with MPI_allreduce in v1dot !'
  end function v2dot
!!***

!!****f* GenBlas/v3dot *
!!
!!  NAME 
!!   v3dot
!!  PURPOSE
!!   Dot product for rank 3 real(double) arrays with sum across all processors
!!  INPUTS
!!   integer :: len - length of array
!!   integer :: stride_a, stride_b - stride of arrays
!!   real(double) :: a, b - arrays
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  real(double) function v3dot(len,a,stride_a,b,stride_b)

    use mpi

    implicit none

    real(double) :: a(:,:,:),b(:,:,:)
    real(double) :: temp
    integer :: len,stride_a,stride_b
    integer :: ierr

    temp = dot(len,a,stride_a,b,stride_b)
    call MPI_allreduce(temp, v3dot, 1, MPI_double_precision, &
         MPI_sum, MPI_COMM_WORLD, ierr)
    if(ierr/=0) write(io_lun,*) 'Problem with MPI_allreduce in v1dot !'
  end function v3dot
!!***

end module GenBlas
