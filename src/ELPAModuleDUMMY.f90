module ELPA_module

  use datatypes
  use mpi
  use GenComms, ONLY: cq_abort, myid
  implicit none

  logical :: flag_elpa_dummy = .true. ! A marker to show no ELPA in compilation
  logical :: flag_use_elpa = .false.  ! This should be false for ELPAModuleDummy
  logical :: flag_elpa_GPU = .false.  ! Should ELPA use GPUs?
  character(len=16) :: elpa_solver = "ELPA1" ! ELPA1 or ELPA2
  character(len=16) :: elpa_kernel = "GENERIC"
  integer :: elpa_API = 20241105
  integer :: merow, mecol

  private
  public :: flag_use_elpa, elpa_solver, elpa_kernel, elpa_API, flag_elpa_dummy, flag_elpa_GPU
  public :: init_ELPA, end_ELPA, ELPA_zhegv

contains

  subroutine init_ELPA (matrix_size, row_size, col_size, desc, info)

    implicit none

    integer, intent(in) :: matrix_size, row_size, col_size
    integer, intent(in) :: desc(9)
    integer, intent(out) :: info

    if(flag_use_elpa) then
       call cq_abort("init_ELPA: init_ELPA is called though ELPA_module is not compiled")
    else
       call cq_abort("init_ELPA: init_ELPA is called even though use_elpa is false")
    endif

    return
  end subroutine init_ELPA

  subroutine end_ELPA( info )

    implicit none

    integer, intent(out) :: info

    call cq_abort("end_ELPA: end_ELPA should not be called")

    return
  end subroutine end_ELPA

  subroutine ELPA_zhegv( mode, matrix_size, row_size, col_size, &
       Hmat, Smat, Wvec, Zmat, info )

    implicit none

    character(len=1), intent(in) :: mode
    integer, intent(in) :: matrix_size, row_size, col_size
    complex(double_cplx), intent(inout) :: Hmat(row_size,col_size)
    complex(double_cplx), intent(inout) :: Smat(row_size,col_size)
    real(double),    intent(out)   :: Wvec(matrix_size)
    complex(double_cplx), intent(out)   :: Zmat(row_size,col_size)
    integer, intent(out) :: info

    call cq_abort("ELPA_zhev: CONQUEST should be compiled with ELPA")

    return
  end subroutine ELPA_zhegv

end module ELPA_module
