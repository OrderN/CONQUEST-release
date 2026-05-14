! Interface to hipfort hipFFTW
!
! Modifications
! 2019/05/29 16:54 dave
!  Added real-to-real (i.e. cosine and sine) transforms
module hipfft_interface_module
  use iso_c_binding, only : c_ptr, c_loc, c_sizeof
  implicit none

  type(c_ptr), save :: plan

contains

  subroutine hipfft3_init_wrapper( nsize )
    use hipfort_hipfft, only : hipfftplan3d, HIPFFT_Z2Z

    integer, intent(in)  :: nsize
    integer :: ierr

    ierr = hipfftplan3d(plan, nsize, nsize, nsize, HIPFFT_Z2Z)
  end subroutine hipfft3_init_wrapper

  subroutine hipfft3_exec_wrapper( cdata_h, nsize, isign )
    use hipfort, only : hipMalloc, hipMemcpy, hipMemcpyHostToDevice, hipMemcpyDeviceToHost
    use hipfort_hipfft, only : hipfftexecz2z, HIPFFT_BACKWARD, HIPFFT_FORWARD
    use datatypes, only : double_cplx

    integer :: nsize, isign
    complex(double_cplx), intent(inout) :: cdata_h(nsize,nsize,nsize)

    integer :: ierr
    complex(double_cplx), pointer, contiguous :: cdata_d(:,:,:)
    !type(c_ptr) :: cdata_d_ptr

    ierr = hipMalloc(cdata_d, [nsize,nsize,nsize])
    ierr = hipMemcpy(cdata_d, cdata_h, nsize*nsize*nsize, hipMemcpyHostToDevice)

    !cdata_d_ptr = c_loc(cdata_d(1,1,1))

    if( (-1)*isign == -1 ) then ! forward
      ierr = hipfftexecz2z(plan, cdata_d, cdata_d, HIPFFT_FORWARD)
    else if( (-1)*isign == +1 ) then ! reverse
      ierr = hipfftexecz2z(plan, cdata_d, cdata_d, HIPFFT_BACKWARD)
    end if

    ierr = hipMemcpy(cdata_h, cdata_d, nsize*nsize*nsize, hipMemcpyDeviceToHost)

  end subroutine hipfft3_exec_wrapper

  subroutine hipfft3_dest_wrapper( )
    use hipfort_hipfft, only : hipfftdestroy_
    integer :: ierr
    ierr = hipfftdestroy_(plan)
  end subroutine hipfft3_dest_wrapper

end module hipfft_interface_module
