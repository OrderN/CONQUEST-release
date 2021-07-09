! Interface to FFTW
!
! Modifications
! 2019/05/29 16:54 dave
!  Added real-to-real (i.e. cosine and sine) transforms
module fft_interface_module

  use datatypes

  implicit none

  integer, parameter   :: FFTW_FORWARD  = -1
  integer, parameter   :: FFTW_BACKWARD = +1
  integer, parameter   :: FFTW_REDFT00  = 3
  integer, parameter   :: FFTW_RODFT00  = 7
  integer, parameter   :: FFTW_ESTIMATE = 64 

  integer(wide), save :: planx_for, planx_rev
  integer(wide), save :: plany_for, plany_rev
  integer(wide), save :: planz_for, planz_rev
  integer(wide), save :: plan3_for, plan3_rev, plan3_for_nonu, plan3_rev_nonu
  integer(wide), save :: cos_plan_for, cos_plan_rev
  integer(wide), save :: sin_plan

contains

  subroutine fftx_init_wrapper( len )
    use datatypes

    implicit none

    integer,intent(in)   :: len
    complex(double_cplx) :: dummy

    ! for forward
    call dfftw_plan_dft_1d( planx_for, &
         len, dummy, dummy, FFTW_FORWARD,  FFTW_ESTIMATE )
    ! for reverse
    call dfftw_plan_dft_1d( planx_rev, &
         len, dummy, dummy, FFTW_BACKWARD, FFTW_ESTIMATE )

    return
  end subroutine fftx_init_wrapper

  subroutine fftx_exec_wrapper( cdata, size, cols, len, isign )
    use datatypes

    implicit none

    integer,intent(in) :: cols, len, isign, size
    complex(double_cplx), intent(inout) :: cdata(size)

    integer :: i

    if( isign==-1 ) then ! forward
       do i=0, cols-1
          call dfftw_execute_dft( planx_for, cdata(i*len+1), cdata(i*len+1) )
       end do
    end if
    if( isign==+1 ) then ! reverse
       do i=0, cols-1
          call dfftw_execute_dft( planx_rev, cdata(i*len+1), cdata(i*len+1) )
       end do
    end if
    cdata(:) = cdata(:)*SQRT(1.0_double/REAL(len,double))

    return
  end subroutine fftx_exec_wrapper

  subroutine ffty_init_wrapper( len )
    use datatypes

    implicit none

    integer,intent(in)   :: len
    complex(double_cplx) :: dummy

    ! for forward
    call dfftw_plan_dft_1d( plany_for, &
         len, dummy, dummy, FFTW_FORWARD,  FFTW_ESTIMATE )
    ! for reverse
    call dfftw_plan_dft_1d( plany_rev, &
         len, dummy, dummy, FFTW_BACKWARD, FFTW_ESTIMATE )

    return
  end subroutine ffty_init_wrapper

  subroutine ffty_exec_wrapper( cdata, size, cols, len, isign )
    use datatypes

    implicit none

    integer :: cols, len, isign, size
    complex(double_cplx), intent(inout) :: cdata(size)

    integer :: i

    if( isign==-1 ) then ! forward
       do i=0, cols-1
          call dfftw_execute_dft( plany_for, cdata(i*len+1), cdata(i*len+1) )
       end do
    end if
    if( isign==+1 ) then ! reverse
       do i=0, cols-1
          call dfftw_execute_dft( plany_rev, cdata(i*len+1), cdata(i*len+1) )
       end do
    end if
    cdata(:) = cdata(:)*SQRT(1.0_double/REAL(len,double))

    return
  end subroutine ffty_exec_wrapper

  subroutine fftz_init_wrapper( len )
    use datatypes

    implicit none

    integer,intent(in)   :: len
    complex(double_cplx) :: dummy

    ! for forward
    call dfftw_plan_dft_1d( planz_for, &
         len, dummy, dummy, FFTW_FORWARD,  FFTW_ESTIMATE )
    ! for reverse
    call dfftw_plan_dft_1d( planz_rev, &
         len, dummy, dummy, FFTW_BACKWARD, FFTW_ESTIMATE )

    return
  end subroutine fftz_init_wrapper

  subroutine fftz_exec_wrapper( cdata, size, cols, len, isign )
    use datatypes

    implicit none

    integer,intent(in) :: cols, len, isign, size
    complex(double_cplx), intent(inout) :: cdata(size)

    integer :: i

    if( isign==-1 ) then ! forward
       do i=0, cols-1
          call dfftw_execute_dft( planz_for, cdata(i*len+1), cdata(i*len+1) )
       end do
    end if
    if( isign==+1 ) then ! reverse
       do i=0, cols-1
          call dfftw_execute_dft( planz_rev, cdata(i*len+1), cdata(i*len+1) )
       end do
    end if
    cdata(:) = cdata(:)*SQRT(1.0_double/REAL(len,double))

    return
  end subroutine fftz_exec_wrapper


  subroutine fft3_init_wrapper( nsize )
    use datatypes

    implicit none

    integer, intent(in)  :: nsize
    complex(double_cplx) :: dummy

    ! for forward
    call dfftw_plan_dft_3d( plan3_for, &
         nsize, nsize, nsize, dummy, dummy, FFTW_FORWARD, FFTW_ESTIMATE )
    ! for reverse
    call dfftw_plan_dft_3d( plan3_rev, &
         nsize, nsize, nsize, dummy, dummy, FFTW_BACKWARD, FFTW_ESTIMATE )

    return
  end subroutine fft3_init_wrapper

  subroutine fft3_exec_wrapper( cdata, nsize, isign )
    use datatypes

    implicit none

    integer :: nsize, isign
    complex(double_cplx), intent(inout) :: cdata(nsize,nsize,nsize)

    if( (-1)*isign == -1 ) then ! forward
       call dfftw_execute_dft( plan3_for, cdata, cdata )
    else if( (-1)*isign == +1 ) then ! reverse
       call dfftw_execute_dft( plan3_rev, cdata, cdata )
    end if

    return
  end subroutine fft3_exec_wrapper

  subroutine fft3_init_wrapper_nonu( nsizex, nsizey, nsizez )
    use datatypes

    implicit none

    integer, intent(in)  :: nsizex, nsizey, nsizez
    complex(double_cplx) :: dummy

    ! for forward
    call dfftw_plan_dft_3d( plan3_for_nonu, &
         nsizex, nsizey, nsizez, dummy, dummy, FFTW_FORWARD, FFTW_ESTIMATE )
    ! for reverse
    call dfftw_plan_dft_3d( plan3_rev_nonu, &
         nsizex, nsizey, nsizez, dummy, dummy, FFTW_BACKWARD, FFTW_ESTIMATE )

    return
  end subroutine fft3_init_wrapper_nonu

  subroutine fft3_exec_wrapper_nonu( cdata, nsizex, nsizey, nsizez, isign )
    use datatypes

    implicit none

    integer :: nsizex, nsizey, nsizez, isign
    complex(double_cplx), intent(inout) :: cdata(nsizex,nsizey,nsizez)

    if( (-1)*isign == -1 ) then ! forward
       call dfftw_execute_dft( plan3_for_nonu, cdata, cdata )
    else if( (-1)*isign == +1 ) then ! reverse
       call dfftw_execute_dft( plan3_rev_nonu, cdata, cdata )
    end if

    return
  end subroutine fft3_exec_wrapper_nonu

  subroutine fft3_dest_wrapper( )
    use datatypes

    implicit none

    ! for forward
    call dfftw_destroy_plan( plan3_for )
    ! for reverse
    call dfftw_destroy_plan( plan3_rev )

    return
  end subroutine fft3_dest_wrapper

  subroutine fft3_dest_wrapper_nonu( )
    use datatypes

    implicit none

    ! for forward
    call dfftw_destroy_plan( plan3_for_nonu )
    ! for reverse
    call dfftw_destroy_plan( plan3_rev_nonu )

    return
  end subroutine fft3_dest_wrapper_nonu

  subroutine cosft_init_wrapper( len )
    
    use datatypes

    implicit none

    integer,intent(in)   :: len
    real(double) :: dummy

    ! for forward
    call dfftw_plan_r2r_1d( cos_plan_for, len, dummy, dummy, FFTW_REDFT00, FFTW_ESTIMATE )
    ! for reverse
    call dfftw_plan_r2r_1d( cos_plan_rev, len, dummy, dummy, FFTW_REDFT00, FFTW_ESTIMATE )

    return
  end subroutine cosft_init_wrapper

  subroutine cosft_exec_wrapper( data, len, isign )
    
    use datatypes

    implicit none

    integer, intent(in) :: len, isign
    real(double), intent(inout) :: data(len)

    if( isign==-1 ) then ! forward
       call dfftw_execute_r2r( cos_plan_for, data, data )
    end if
    if( isign==+1 ) then ! reverse
       call dfftw_execute_r2r( cos_plan_rev, data, data )
       data = data/real(2*(len-1),double)
    end if

    return
  end subroutine cosft_exec_wrapper

  subroutine cosft_dest_wrapper( )
    
    use datatypes

    implicit none

    call dfftw_destroy_plan( cos_plan_for )
    call dfftw_destroy_plan( cos_plan_rev )

    return
  end subroutine cosft_dest_wrapper

  subroutine sinft_init_wrapper( len )
    
    use datatypes

    implicit none

    integer,intent(in)   :: len
    real(double) :: dummy

    ! for forward & reverse
    call dfftw_plan_r2r_1d( sin_plan, len, dummy, dummy, FFTW_RODFT00, FFTW_ESTIMATE )

    return
  end subroutine sinft_init_wrapper

  ! For the sine transform, the forward and reverse are the same except for a scaling factor
  subroutine sinft_exec_wrapper( data, len, isign )
    
    use datatypes

    implicit none

    integer, intent(in) :: len, isign
    real(double), intent(inout) :: data(len)

    call dfftw_execute_r2r( sin_plan, data, data )
    if( isign==+1 ) then ! reverse
       data = data/real(2*len+2,double)
    end if

    return
  end subroutine sinft_exec_wrapper

  subroutine sinft_dest_wrapper( )
    
    use datatypes

    implicit none

    call dfftw_destroy_plan( sin_plan )

    return
  end subroutine sinft_dest_wrapper

end module fft_interface_module
