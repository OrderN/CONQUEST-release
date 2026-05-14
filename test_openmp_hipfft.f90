program test_openmp_hipfft
  use iso_c_binding, only : c_int, c_ptr
  use hipfort                      
  use hipfort_hipfft               
  implicit none

  integer(c_int), parameter :: N = 1024
  integer(c_int) :: i

  type(c_ptr) :: plan
  
  integer, parameter :: double_cplx = selected_real_kind( 6, 70 )

  ! Host arrays
  complex(double_cplx), allocatable, target :: x(:)
  complex(double_cplx), allocatable :: x_ref(:)
  complex(double_cplx), pointer :: x_ptr(:)

  ! Error codes
  integer :: ierr
  integer :: fftstat
  real :: maxerr

  allocate(x(N), x_ref(N))

  ! Initialise a simple signal on host
  do i = 1, N
     x(i) = cmplx(real(i-1), -0.5, kind=double_cplx)
  end do
  x_ref = x
  
  !$omp target data map(tofrom: x(1:N))
    !$omp target data use_device_addr(x(1:N))
      x_ptr => x

      ! Create 1D FFT plan length N, complex-to-complex, batch=1 
      fftstat = hipfftPlan1d(plan, N, HIPFFT_Z2Z, 1)
      if (fftstat /= HIPFFT_SUCCESS) stop "hipfftPlan1d failed"

      ! Forward FFT in-place (idata==odata is allowed) 
      fftstat = hipfftExecZ2Z(plan, x_ptr, x_ptr, HIPFFT_FORWARD)
      if (fftstat /= HIPFFT_SUCCESS) stop "hipfftExecC2C forward failed"

      ! Backward (inverse) FFT in-place 
      fftstat = hipfftExecZ2Z(plan, x_ptr, x_ptr, HIPFFT_BACKWARD)
      if (fftstat /= HIPFFT_SUCCESS) stop "hipfftExecC2C backward failed"

      fftstat = hipDeviceSynchronize_()
      if (fftstat /= HIPFFT_SUCCESS) stop "hipDeviceSynchronize_ failed"

    !$omp end target data
  !$omp end target data

  ! Normalise 
  x = x / real(N)

  ! Compute max error
  maxerr = 0.0
  do i = 1, N
     maxerr = max(maxerr, abs(x(i) - x_ref(i)))
  end do

  write(*,'(a,i0)') "N = ", N
  write(*,'(a,es12.4)') "max |x_back - x_ref| = ", maxerr

  ! Cleanup
  fftstat = hipfftDestroy(plan)    
  !ierr    = hipFree(x_d)           

  deallocate(x, x_ref)

  write(*,*) "Finished hipfort + hipFFT forward/backward test"
end program test_openmp_hipfft