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

  subroutine hipfft3_exec_wrapper( cdata, nsize, isign )
    use hipfort, only : hipMalloc, hipMemcpy, hipMemcpyHostToDevice, hipMemcpyDeviceToHost
    use hipfort_hipfft, only : hipfftexecz2z, HIPFFT_BACKWARD, HIPFFT_FORWARD, HIPFFT_SUCCESS
    use datatypes, only : double_cplx

    integer :: nsize, isign
    complex(double_cplx), intent(inout) :: cdata(nsize,nsize,nsize)

    integer :: ierr
    complex(double_cplx), pointer, contiguous :: cdata_d_ptr(:,:,:)

    !$omp target data map(tofrom: cdata(1:nsize, 1:nsize, 1:nsize))
      !$omp target data use_device_addr(cdata(1:nsize, 1:nsize, 1:nsize))
        ! Create pointer to device data to pass to execute kernels
        cdata_d_ptr => cdata
        
        ! Execute kernel
        if( (-1)*isign == -1 ) then ! forward
          ierr = hipfftexecz2z(plan, cdata_d_ptr, cdata_d_ptr, HIPFFT_FORWARD)
        else if( (-1)*isign == +1 ) then ! reverse
          ierr = hipfftexecz2z(plan, cdata_d_ptr, cdata_d_ptr, HIPFFT_BACKWARD)
        end if
        if (ierr /= HIPFFT_SUCCESS) stop "hipfftexecz2z failed"

        ! Execute kernels run asychronously, therefore we must synchronise here
        ierr = hipDeviceSynchronize_()
        if (ierr /= HIPFFT_SUCCESS) stop "hipDeviceSynchronize_ failed"
    !$omp end target data
  !$omp end target data

  end subroutine hipfft3_exec_wrapper

  subroutine hipfft3_dest_wrapper( )
    use hipfort_hipfft, only : hipfftdestroy_
    integer :: ierr
    ierr = hipfftdestroy_(plan)
  end subroutine hipfft3_dest_wrapper

end module hipfft_interface_module
