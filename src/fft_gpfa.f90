module fft_interface_module

  use datatypes

  implicit none

  real(double), save, dimension(:), allocatable :: trigs_x
  real(double), save, dimension(:), allocatable :: trigs_y
  real(double), save, dimension(:), allocatable :: trigs_z
  real(double), save, dimension(:,:), allocatable :: trigs_3

contains

  subroutine fftx_init_wrapper( len )
    use datatypes
    implicit none

    integer,intent(in) :: len

    allocate( trigs_x(2*len) )
    call setgpfa( trigs_x, len )

    return
  end subroutine fftx_init_wrapper

  subroutine fftx_exec_wrapper( cdata, size, cols, len, isign )
    use datatypes
    implicit none

    integer,intent(in) :: cols, len, isign, size
    complex(double_cplx), intent(inout) :: cdata(size)

    real(double) :: scale
    integer :: i,j

    real(double), allocatable, dimension(:) :: a, b

    scale = SQRT(1.0_double/REAL(len,double))

    allocate( a(len), b(len) )

    do i=0, cols-1
       do j=1,len
          a(j) = real(cdata(i*len+j),double)
          b(j) = aimag(cdata(i*len+j))
       end do
       call gpfa( a, b, trigs_x, 1, 1, len, 1, isign )
       do j=1,len
          cdata(i*len+j) = scale*cmplx(a(j),b(j),double_cplx)
       end do
    end do

    deallocate(a,b)

    return
  end subroutine fftx_exec_wrapper

  subroutine ffty_init_wrapper( len )
    use datatypes
    implicit none

    integer,intent(in) :: len

    allocate( trigs_y(2*len) )
    call setgpfa( trigs_y, len )

    return
  end subroutine ffty_init_wrapper

  subroutine ffty_exec_wrapper( cdata, size, cols, len, isign )
    use datatypes
    implicit none

    integer,intent(in) :: cols, len, isign, size
    complex(double_cplx), intent(inout) :: cdata(size)

    real(double) :: scale
    integer :: i,j

    real(double), allocatable, dimension(:) :: a, b

    scale = SQRT(1.0_double/REAL(len,double))

    allocate( a(len), b(len) )

    do i=0, cols-1
       do j=1,len
          a(j) = real(cdata(i*len+j),double)
          b(j) = aimag(cdata(i*len+j))
       end do
       call gpfa( a, b, trigs_y, 1, 1, len, 1, isign )
       do j=1,len
          cdata(i*len+j) = scale*cmplx(a(j),b(j),double_cplx)
       end do
    end do

    deallocate(a,b)

    return
  end subroutine ffty_exec_wrapper

  subroutine fftz_init_wrapper( len )
    use datatypes
    implicit none

    integer,intent(in) :: len

    allocate( trigs_z(2*len) )
    call setgpfa( trigs_z, len )

    return
  end subroutine fftz_init_wrapper

  subroutine fftz_exec_wrapper( cdata, size, cols, len, isign )
    use datatypes
    implicit none

    integer,intent(in) :: cols, len, isign, size
    complex(double_cplx), intent(inout) :: cdata(size)

    real(double) :: scale
    integer :: i,j

    real(double), allocatable, dimension(:) :: a, b

    scale = SQRT(1.0_double/REAL(len,double))

    allocate( a(len), b(len) )

    do i=0, cols-1
       do j=1,len
          a(j) = real(cdata(i*len+j),double)
          b(j) = aimag(cdata(i*len+j))
       end do
       call gpfa( a, b, trigs_z, 1, 1, len, 1, isign )
       do j=1,len
          cdata(i*len+j) = scale*cmplx(a(j),b(j),double_cplx)
       end do
    end do

    deallocate(a,b)

    return
  end subroutine fftz_exec_wrapper

  subroutine fft3_init_wrapper( nsize )
    use datatypes
    implicit none

    integer, intent(in) :: nsize

    allocate( trigs_3(2*nsize,3) )
    call setgpfa(trigs_3(1,1),nsize)
    trigs_3(:,2) = trigs_3(:,1)
    trigs_3(:,3) = trigs_3(:,1)

    return
  end subroutine fft3_init_wrapper

  subroutine fft3_exec_wrapper( cdata, nsize, isign )
    use datatypes
    implicit none

    complex(double_cplx), intent(inout) :: cdata(nsize,nsize,nsize)
    integer, intent(in) :: nsize, isign
    integer, dimension(3) :: nn

    nn(:) = nsize
    call GPF3D( cdata, trigs_3, 2*nsize, nsize, nn, isign )

    return
  end subroutine fft3_exec_wrapper

end module fft_interface_module
