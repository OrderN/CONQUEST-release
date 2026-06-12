module ELPA_module

  use datatypes
  use mpi
!!$   use omp
  use elpa
  use GenComms, ONLY: cq_abort, cq_warn, myid

  implicit none

  logical :: flag_elpa_dummy = .false. ! A marker to show ELPA in compilation
  logical :: flag_use_elpa = .false.  ! whether we use ELPA or not
  logical :: flag_elpa_GPU = .false.  ! Should ELPA use GPUs?
  character(len=16) :: elpa_solver = "ELPA1" ! ELPA1 or ELPA2
  character(len=16) :: elpa_kernel = "GENERIC"
  integer :: elpa_API = 20241105
  integer :: merow, mecol

  class(elpa_t), pointer :: elp

  private
  public :: flag_use_elpa, elpa_solver, elpa_kernel, elpa_API, flag_elpa_dummy, flag_elpa_GPU
  public :: init_ELPA, end_ELPA, ELPA_zhegv

contains

  subroutine init_ELPA (matrix_size, row_size, col_size, desc, info)

    implicit none

    integer, intent(in) :: matrix_size, row_size, col_size
    integer, intent(in) :: desc(9)
    integer, intent(out) :: info

    integer :: context, block_size_r, block_size_c
    integer :: numrows, numcols, merow, mecol  ! numrows= proc_rows, numcols=proc_cols
    character(len=12) :: subname = "init_ELPA:  "

    context      = desc(2)
    block_size_r = desc(5)
    block_size_c = desc(6)

    if( block_size_r /= block_size_c ) then ! restriction for ELPA
       call cq_abort("Diag.BlockSizeR and Diag.BlockSizeC not same !", &
            block_size_r, block_size_c )
    end if

    !To get the information of blacs grid
    call blacs_gridinfo( context, numrows, numcols, merow, mecol )

    if( mod(numcols,numrows) /= 0 ) then ! restriction for ELPA
       call cq_warn(subname,"Diag.ProcRows is not a factor of Diag.ProcCols !",numrows,numcols)
    end if

    if( matrix_size <= block_size_r*numrows ) then ! restriction for ELPA
       call cq_abort("Diag.BlockSizeR should be less than or equal to", (matrix_size-1)/numrows )
    end if
    if( matrix_size <= block_size_c*numcols ) then ! restriction for ELPA
       call cq_abort("Diag.BlockSizeC should be less than or rqual to", (matrix_size-1)/numcols )
    end if

    info = elpa_init(elpa_API)
    if( info /= ELPA_OK ) call cq_abort("ELPA_Init: ELPA API version not supported")

    elp => elpa_allocate(info)

    call elp%set( "na",  matrix_size, info )
    if( info /= ELPA_OK ) call cq_abort("ELPA_Init: Could not set parameter na")
    call elp%set( "nev", matrix_size, info )
    if( info /= ELPA_OK ) call cq_abort("ELPA_Init: Could not set parameter nev")
    call elp%set( "local_nrows", row_size, info )
    if( info /= ELPA_OK ) call cq_abort("ELPA_Init: Could not set parameter local_nrows")
    call elp%set( "local_ncols", col_size, info )
    if( info /= ELPA_OK ) call cq_abort("ELPA_Init: Could not set parameter local_ncols")
    call elp%set( "nblk", block_size_r, info )
    if( info /= ELPA_OK ) call cq_abort("ELPA_Init: Could not set parameter nblk")
    call elp%set( "mpi_comm_parent", MPI_COMM_WORLD, info )
    if( info /= ELPA_OK ) call cq_abort("ELPA_Init: Could not set parameter mpi_comm_parent")
    call elp%set( "process_row", merow, info )
    if( info /= ELPA_OK ) call cq_abort("ELPA_Init: Could not set parameter process_row")
    call elp%set( "process_col", mecol, info )
    if( info /= ELPA_OK ) call cq_abort("ELPA_Init: Could not set parameter process_col")
    call elp%set( "blacs_context", context, info )
    if( info /= ELPA_OK ) call cq_abort("ELPA_Init: Could not set parameter blacs_cotext")

    if(flag_elpa_GPU) info = elp%setup()
    select case( elpa_solver )
    case("ELPA1")
       call elp%set( "solver", ELPA_SOLVER_1STAGE, info ) ! ELPA1
    case("ELPA2")
       call elp%set( "solver", ELPA_SOLVER_2STAGE, info ) ! ELPA2

       select case( elpa_kernel )
       case("GENERIC")
          call elp%set( "complex_kernel", ELPA_2STAGE_COMPLEX_GENERIC, info )
       case("GENERIC_SIMPLE")
          call elp%set( "complex_kernel", ELPA_2STAGE_COMPLEX_GENERIC_SIMPLE, info )
       case("SSE_ASSEMBLY")
          call elp%set( "complex_kernel", ELPA_2STAGE_COMPLEX_SSE_ASSEMBLY, info )
       case("SSE_BLOCK1")
          call elp%set( "complex_kernel", ELPA_2STAGE_COMPLEX_SSE_BLOCK1, info )
       case("SSE_BLOCK2")
          call elp%set( "complex_kernel", ELPA_2STAGE_COMPLEX_SSE_BLOCK2, info )
       case("AVX_BLOCK1")
          call elp%set( "complex_kernel", ELPA_2STAGE_COMPLEX_AVX_BLOCK1, info )
       case("AVX_BLOCK2")
          call elp%set( "complex_kernel", ELPA_2STAGE_COMPLEX_AVX_BLOCK2, info )
       case("AVX2_BLOCK1")
          call elp%set( "complex_kernel", ELPA_2STAGE_COMPLEX_AVX2_BLOCK1, info )
       case("AVX2_BLOCK2")
          call elp%set( "complex_kernel", ELPA_2STAGE_COMPLEX_AVX2_BLOCK2, info )
       case default
          call cq_abort("Invalid Diag.ELPA2Kernal " // trim(elpa_kernel) )
       end select
    case default
       call cq_abort("Invalid Diag.ELPASolver " // trim(elpa_solver) )
    end select

!!$      call elp%set( "omp_threads", omp_get_max_threads(), info )

    if(flag_elpa_GPU) then
       call elp%set("nvidia-gpu",1,info)
       call elp%set("solver",ELPA_SOLVER_1STAGE,info)
       info = elp%setup_gpu()
    else
       info = elp%setup()
    end if
    if( info /= ELPA_OK )  call cq_abort("something wrong in ELPA !")
  end subroutine init_ELPA

  subroutine end_ELPA( info )

    implicit none

    integer, intent(out) :: info

    call elpa_deallocate( elp, info )
    if( info /= ELPA_OK )  call cq_abort("end_ELPA: deallocation error")
    call elpa_uninit( info )
    if( info /= ELPA_OK )  call cq_abort("end_ELPA: uninit error")
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
    integer :: i, j

    if( mode =='N' ) then
       call elp%generalized_eigenvalues( &
            Hmat, Smat, Wvec, .false., info )
    end if
    if( mode =='V' ) then
       call elp%generalized_eigenvectors( &
            Hmat, Smat, Wvec, Zmat, .false., info )
    end if

  end subroutine ELPA_zhegv

end module ELPA_module
