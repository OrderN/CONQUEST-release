module fft_interface_module

  use datatypes

  implicit none

  complex(double), save, dimension(:), allocatable :: work_x
  complex(double), save, dimension(:), allocatable :: work_y
  complex(double), save, dimension(:), allocatable :: work_z

contains

  subroutine fftx_init_wrapper( len )
    use datatypes

    implicit none

    integer,intent(in) :: len
    complex(double_cplx) :: dummy

    allocate( work_x(2*len) )
    call zfft1d( dummy, len, 0, work_x )

    return
  end subroutine fftx_init_wrapper

  subroutine fftx_exec_wrapper( cdata, size, cols, len, isign )
    use datatypes

    implicit none

    integer,intent(in) :: cols, len, isign, size
    complex(double_cplx), intent(inout) :: cdata(size)

    integer :: i

    do i=0, cols-1
       call zfft1d( cdata(i*len+1), len, isign, work_x )
    end do

    if( isign == +1 ) then
       cdata(:) = cdata(:)*SQRT(REAL(len,double))
    else
       cdata(:) = cdata(:)*SQRT(1.0_double/REAL(len,double))
    end if

    return
  end subroutine fftx_exec_wrapper

  subroutine ffty_init_wrapper( len )
    use datatypes

    implicit none

    integer,intent(in) :: len
    complex(double_cplx) :: dummy

    allocate( work_y(2*len) )
    call zfft1d( dummy, len, 0, work_y )

    return
  end subroutine ffty_init_wrapper

  subroutine ffty_exec_wrapper( cdata, size, cols, len, isign )
    use datatypes

    implicit none

    integer,intent(in) :: cols, len, isign, size
    complex(double_cplx), intent(inout) :: cdata(size)

    integer :: i

    do i=0, cols-1
       call zfft1d( cdata(i*len+1), len, isign, work_y )
    end do

    if( isign == +1 ) then
       cdata(:) = cdata(:)*SQRT(REAL(len,double))
    else
       cdata(:) = cdata(:)*SQRT(1.0_double/REAL(len,double))
    end if

    return
  end subroutine ffty_exec_wrapper

  subroutine fftz_init_wrapper( len )
    use datatypes

    implicit none

    integer,intent(in) :: len
    complex(double_cplx) :: dummy

    allocate( work_z(2*len) )
    call zfft1d( dummy, len, 0, work_z )

    return
  end subroutine fftz_init_wrapper

  subroutine fftz_exec_wrapper( cdata, size, cols, len, isign )
    use datatypes

    implicit none

    integer,intent(in) :: cols, len, isign, size
    complex(double_cplx), intent(inout) :: cdata(size)

    integer :: i

    do i=0, cols-1
       call zfft1d( cdata(i*len+1), len, isign, work_z )
    end do

    if( isign == +1 ) then
       cdata(:) = cdata(:)*SQRT(REAL(len,double))
    else
       cdata(:) = cdata(:)*SQRT(1.0_double/REAL(len,double))
    end if

    return
  end subroutine fftz_exec_wrapper


  subroutine fft3_init_wrapper( nsize )
    use datatypes

    implicit none

    integer, intent(in) :: nsize
    complex(double_cplx) :: dummy

    call zfft3d( dummy, nsize, nsize, nsize, 0 )

    return
  end subroutine fft3_init_wrapper

  subroutine fft3_exec_wrapper( cdata, nsize, isign )
    use datatypes

    implicit none

    complex(double_cplx), intent(inout) :: cdata(nsize,nsize,nsize)
    integer, intent(in) :: nsize, isign

    call zfft3d( cdata, nsize, nsize, nsize, isign )

    return
  end subroutine fft3_exec_wrapper
end module fft_interface_module
