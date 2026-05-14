! Interface to hipfort hipFFTW
!
! Modifications
! 2019/05/29 16:54 dave
!  Added real-to-real (i.e. cosine and sine) transforms
module hipfft_interface_module
  use iso_c_binding, only : c_ptr
  implicit none

  type(c_ptr), save :: plan

contains

  subroutine hipfft3_init_wrapper( nsize )
    use hipfort_hifft, only : hipfftplan3d_, HIPFFT_C2C

    integer, intent(in)  :: nsize

    call hipfftplan3d_(plan, nsize, nsize, nsize, HIPFFT_C2C)
  end subroutine hipfft3_init_wrapper

  subroutine hipfft3_exec_wrapper( cdata_h, nsize, isign )
    use hipfort, only : hipMalloc, hipMemcpy, hipMemcpyHostToDevice
    use hipfort_hifft, only : hipfftexecc2c_, HIPFFT_BACKWARD, HIPFFT_FORWARD
    use datatypes, only : double_cplx

    integer :: nsize, isign
    complex(double_cplx), intent(inout) :: cdata_h(nsize,nsize,nsize)

    integer :: ierr
    complex(double_cplx), pointer :: cdata_d(nsize,nsize,nsize)
    type(c_ptr) :: cdata_d_ptr

    ierr = hipMalloc(cdata_d, size=[nsize,nsize,nsize])
    ierr = hipMemcpy(cdata_d, cdata_h, [nsize,nsize,nsize], hipMemcpyHostToDevice)

    cdata_d_ptr = c_loc(cdata_d(1,1,1))

    if( (-1)*isign == -1 ) then ! forward
      call hipfftexecc2c_(plan, cdata_d_ptr, cdata_d_ptr, HIPFFT_FORWARD)
    else if( (-1)*isign == +1 ) then ! reverse
      call hipfftexecc2c_(plan, cdata_d_ptr, cdata_d_ptr, HIPFFT_BACKWARD)
    end if
  end subroutine hipfft3_exec_wrapper

  subroutine hipfft3_dest_wrapper( )
    call hipfftdestroy_(plan)
  end subroutine hipfft3_dest_wrapper

end module hipfft_interface_module
