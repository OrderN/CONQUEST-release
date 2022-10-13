! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module GenComms
! -----------------------------------------------------------
! Code area 9: general
! -----------------------------------------------------------

!!****h* Conquest/GenComms *
!!  NAME
!!   GenComms
!!  PURPOSE
!!   Generalises and isolates machine or implementation
!!   specific communications
!!  USES
!!   datatypes, global_module, mpi
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   04/06/2001
!!  MODIFICATION HISTORY
!!   06/06/2001 dave
!!    Corrected bugs - commented out old routines and removed
!!    dependence on pm_mpi.inc in mtmini and mtime. Also fixed
!!    bugs in cq_init, int_gsumv and indented mtmini, mtime
!!   06/06/2001 dave
!!    Implemented gsumv properly
!!   07/06/2001 dave
!!    Added gcopy and gcopy_diff as well as finishing exch
!!    Debugged
!!   08/06/2001 dave
!!    Changed the interface for the gsum - included the vector
!!    and the single variable in the same interface as that's 
!!    what it's there for !
!!    Fixed double_two_gsumv
!!   21/06/2001 dave
!!    Altered cq_abort a little to improve formatting
!!   18/03/2002 dave
!!    Added a little to the header and a static tag for object file id
!!   25/06/2002 drb
!!    Added int, double, char_one, logical and logical_one to gcopy
!!   15:45, 04/02/2003 drb 
!!    Added headers, max and min finders, changed cq_init to init_comms and mtime to real from int
!!   14:30, 26/02/2003 drb 
!!    Added logical_gsum with MPI_LOR (so if ANY of the variables are true, all are)
!!   13:33, 22/09/2003 drb 
!!    Added logical_two_copy
!!   2008/02/06 08:17 dave
!!    Changed for output to file not stdout
!!   2019/12/05 08:17 dave
!!    Added warning output
!!   2020/01/21 17:11 dave
!!    Tidying and removing non-standard FORTRAN
!!   2020/04/21 15:05 dave
!!    Update cq_warn so that warnings only go to Conquest_warnings at iprint=0
!!    Add flag_warnings to alert user to check file
!!  SOURCE
!!
module GenComms

  use datatypes
  use global_module, ONLY: numprocs, io_lun, iprint
  use mpi

  implicit none

  integer, save :: myid, root
  integer, save :: inode, ionode
  integer, save :: warning_lun = 9 ! Check this with input_module to avoid clashes
  logical, save :: flag_warnings = .false.

  integer, parameter :: xn_tag = 106
  integer, parameter :: int_gsum_buf = 10000
  integer, parameter :: double_gsum_buf = 5000
  integer, parameter :: double_cplx_gsum_buf = 2500

  real(double), save :: secnd0

  real(double) :: t1, t2, timings(13)
!!***

!!****m* GenComms/cq_abort *
!!  NAME
!!   cq_abort
!!  PURPOSE
!!   Interface to different abort routines
!!  AUTHOR
!!   D.R. Bowler
!!  MODIFCATION
!!   2014/09/10 LAT
!!    Added cq_abort for boleans
!!  SOURCE
!!  
  interface cq_abort
     module procedure cq_abort_none
     module procedure cq_abort_int
     module procedure cq_abort_real
     module procedure cq_abort_logi
  end interface
!!***

!!****m* GenComms/cq_warn *
!!  NAME
!!   cq_warn
!!  PURPOSE
!!   Interface to different warning routines
!!  AUTHOR
!!   D.R. Bowler
!!  MODIFCATION
!!  SOURCE
!!  
  interface cq_warn
     module procedure cq_warn_none
     module procedure cq_warn_int
     module procedure cq_warn_real
  end interface
!!***

!!****m* GenComms/exchv *
!!  NAME
!!   exchv
!!  PURPOSE
!!   Interface to different exchv routines
!!  AUTHOR
!!   D. R. Bowler
!!  SOURCE
!!  
  interface exchv ! Exchange data with all processors
     module procedure dcplx_exchv
  end interface
!!***

!!****m* GenComms/exch *
!!  NAME
!!   exch
!!  PURPOSE
!!   Interface to different exch routines
!!  AUTHOR
!!   D. R. Bowler
!!  SOURCE
!!  
  interface exch
     module procedure int_exch
  end interface
!!***

  !!****m* GenComms/gcopy *
  !!  NAME
  !!   gcopy
  !!  PURPOSE
  !!   Interface to different gcopy routines
  !!  AUTHOR
  !!   D. R. Bowler
  !!  MODIFICATION HISTORY
  !!   2012/04/03 L.Tong
  !!   - Added interface for double_four_copy
  !!  SOURCE
  !!  
  interface gcopy ! Copy data from root processor to all
     module procedure int_copy
     module procedure int_one_copy
     module procedure int_two_copy
     module procedure double_copy
     module procedure double_one_copy
     module procedure double_two_copy
     module procedure double_three_copy
     module procedure double_four_copy
     module procedure char_one_copy
     module procedure logical_copy
     module procedure logical_one_copy
     module procedure logical_two_copy
  end interface
  !!***

!!****m* GenComms/gcopy_diff *
!!  NAME
!!   
!!  PURPOSE
!!   Interface to different gcopy_diff routines
!!  AUTHOR
!!   D. R. Bowler
!!  SOURCE
!!  
  interface gcopy_diff ! Copies data from specified processor to all
     module procedure integer_copy_diff
     module procedure int_one_copy_diff
  end interface
!!***

!!****m* GenComms/gsum *
!!  NAME
!!   gsum
!!  PURPOSE
!!   Interface to different gsum routines
!!  AUTHOR
!!   D. R. Bowler
!!  SOURCE
!!  
  interface gsum ! Global sum on variables
     module procedure int_gsum
     module procedure double_gsum
     module procedure dcplx_gsum
     module procedure int_one_gsumv
     module procedure double_one_gsumv
     module procedure double_two_gsumv
     module procedure double_three_gsumv
     module procedure double_four_gsumv
     module procedure dcplx_one_gsumv
     module procedure logical_gsum
  end interface
!!***

!!****m* GenComms/gmin *
!!  NAME
!!   gmin
!!  PURPOSE
!!   Interface to different gmin routines
!!  AUTHOR
!!   D. R. Bowler
!!  SOURCE
!!  
  interface gmin
     module procedure double_gmin
  end interface
!!***

!!****m* GenComms/gmax *
!!  NAME
!!   gmax
!!  PURPOSE
!!   Interface to different gmax routines
!!  AUTHOR
!!   D. R. Bowler
!!  SOURCE
!!  
  interface gmax
     module procedure integer_gmax
     module procedure double_gmax
  end interface
!!***

contains

! -----------------------------------------------------------
! Subroutine init_comms
! -----------------------------------------------------------

!!****f* GenComms/init_comms *
!!
!!  NAME 
!!   init_comms
!!  USAGE
!! 
!!  PURPOSE
!!   Initialises communications
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!   26/06/2001 dave
!!    Added ROBODoc header
!!   30/08/2001 dave
!!    Moved to GenComms module
!!   2007/05/11 08:00 dave
!!    Added initialisation call for mtmini
!!  SOURCE
!!
  subroutine init_comms(myid,number_of_procs)

    ! Module usage
    use mpi

    implicit none

    ! Passed variables
    integer :: myid,number_of_procs

    ! Local variables
    integer :: ierr

    call MPI_INIT(ierr)
    if(ierr.ne.0) write(io_lun,*) 'ierr is ',ierr
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, number_of_procs, ierr )
    inode  = myid+1
    ionode = 1
    if(inode==ionode) open(unit=warning_lun, file='Conquest_warnings')
    call mtmini()
    return
  end subroutine init_comms
!!***

! -----------------------------------------------------------
! Subroutine end_comms
! -----------------------------------------------------------

!!****f* GenComms/end_comms *
!!
!!  NAME 
!!   end_comms
!!  USAGE
!! 
!!  PURPOSE
!!   Shuts down communications
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!   26/06/2001 dave
!!    ROBODoc header
!!   30/08/2001 dave
!!    Included in generic_comms and changed to end_comms
!!   2007/05/11 07:58 dave
!!    Added final call to mtime to give overall run time
!!  SOURCE
!!
  subroutine end_comms

    ! Module usage
    use mpi

    implicit none

    ! Local variables
    integer :: ierr

    if(inode==ionode) close(warning_lun)
    if(myid==0) write(io_lun,fmt='(/4x,"Total run time was: ",f19.3," seconds")') mtime()*0.001_double
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Finalize(ierr)
    return
  end subroutine end_comms
!!***

! -----------------------------------------------------------
! Subroutine cq_abort_none
! -----------------------------------------------------------

!!****f* GenComms/cq_abort_none *
!!
!!  NAME 
!!   cq_abort_none - aborts Conquest and prints a message
!!  USAGE
!!   cq_abort_none(message indicating routine and problem)
!!  PURPOSE
!!   Aborts Conquest in a reasonably graceful manner, outputting
!!   a message indicating routine and problem
!!  INPUTS
!!   character(len=80) :: message - the abort message
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   05/06/2001 dave
!!  MODIFICATION HISTORY
!!   14/06/2001 dave
!!    Created cq_abort_int and cq_abort_dble
!!   15/06/2001 dave
!!    Created cq_abort_none - no numbers passed
!!   15/06/2001 dave
!!    Merged back into a single routine with two optional arguments
!!   19/06/2001 dave
!!    Fixed formatting for second write
!!   19/07/2001 dave
!!    Reordered ifs so that BOTH integers are printed if present
!!   2006/07/06 08:00 dave
!!    Changed to take no extra arguments for new interface cq_abort
!!  SOURCE
!!
  subroutine cq_abort_none(message)

    implicit none

    ! Passed variables
    character(len=*) :: message

    ! Local variables
    integer :: ierror
    
    write(*,1)       myid+1
    write(*,4)       message
    write(io_lun,1)  myid+1
    write(io_lun,4)  message
1   format(2x,'Error in process ',i4)
4   format(2x,a)
    call flush(io_lun)
    call flush(warning_lun)
    call MPI_abort( MPI_comm_world, 1, ierror )
    stop

    return
  end subroutine cq_abort_none
!!***

! -----------------------------------------------------------
! Subroutine cq_abort_int
! -----------------------------------------------------------

!!****f* GenComms/cq_abort_int *
!!
!!  NAME 
!!   cq_abort_int - aborts Conquest and prints a message
!!  USAGE
!!   cq_abort_int(message indicating routine and problem,int1,int2)
!!  PURPOSE
!!   Aborts Conquest in a reasonably graceful manner, outputting
!!   a message indicating routine and problem and one or two integers
!!  INPUTS
!!   character(len=80) :: message - the abort message
!!   integer :: int1
!!   integer, OPTIONAL :: int2
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   05/06/2001 dave
!!  MODIFICATION HISTORY
!!   14/06/2001 dave
!!    Created cq_abort_int and cq_abort_dble
!!   15/06/2001 dave
!!    Created cq_abort_none - no numbers passed
!!   15/06/2001 dave
!!    Merged back into a single routine with two optional arguments
!!   19/06/2001 dave
!!    Fixed formatting for second write
!!   19/07/2001 dave
!!    Reordered ifs so that BOTH integers are printed if present
!!   2006/07/06 08:01 dave
!!    Changed so that first integer is required (for new interface)
!!  SOURCE
!!
  subroutine cq_abort_int(message,int1,int2)

    implicit none

    ! Passed variables
    character(len=*)  :: message
    integer           :: int1
    integer, OPTIONAL :: int2

    ! Local variables
    integer :: ierror
    
    write (*,1)       myid+1
    write (io_lun,1)  myid+1
    if(PRESENT(int2)) then
       write(*,3)      message, int1, int2
       write(io_lun,3) message, int1, int2
    else 
       write(*,2)      message, int1
       write(io_lun,2) message, int1
    endif

1   format(2x,'Error in process ',i4)
2   format(2x,a,i10)
3   format(2x,a,2i10)

    call flush(io_lun)
    call flush(warning_lun)
    call MPI_abort( MPI_comm_world, 1, ierror )
    stop

    return
  end subroutine cq_abort_int
!!***

! -----------------------------------------------------------
! Subroutine cq_abort_logi
! -----------------------------------------------------------

!!****f* GenComms/cq_abort_logi *
!!
!!  NAME 
!!   cq_abort_logi - aborts Conquest and prints a message
!!  USAGE
!!   cq_abort_int(message indicating routine and problem,logi1,logi2)
!!  PURPOSE
!!   Aborts Conquest in a reasonably graceful manner, outputting
!!   a message indicating routine and problem and one or two boleans
!!  INPUTS
!!   character(len=80) :: message - the abort message
!!   integer           :: logi1
!!   integer, OPTIONAL :: logi2
!!  USES
!! 
!!  AUTHOR
!!   L.A. Truflandier
!!  CREATION DATE
!!   2014/09/10 lat
!!  MODIFICATION
!!  SOURCE
!!
  subroutine cq_abort_logi(message,logi1,logi2)

    implicit none

    ! Passed variables
    character(len=*)  :: message
    logical           :: logi1
    logical, optional :: logi2

    ! Local variables
    integer :: ierror
    
    write (*,1)       myid+1
    write (io_lun,1)  myid+1
    if(present(logi2)) then
       write(*,3)      message, logi1, logi2
       write(io_lun,3) message, logi1, logi2
    else 
       write(*,2)      message, logi1
       write(io_lun,2) message, logi1
    endif

1   format(2x,'Error in process ',i4)
2   format(2x,a, l2)
3   format(2x,a,2l2)

    call flush(io_lun)
    call flush(warning_lun)
    call MPI_abort( MPI_comm_world, 1, ierror )
    stop

    return
  end subroutine cq_abort_logi
!!***

! -----------------------------------------------------------
! Subroutine cq_abort_real
! -----------------------------------------------------------

!!****f* GenComms/cq_abort_real *
!!
!!  NAME 
!!   cq_abort_real - aborts Conquest and prints a message
!!  USAGE
!!   cq_abort(message indicating routine and problem,real1,real2)
!!  PURPOSE
!!   Aborts Conquest in a reasonably graceful manner, outputting
!!   a message indicating routine and problem and optionally two reals
!!  INPUTS
!!   character(len=80) :: message - the abort message
!!   integer :: real1
!!   integer, OPTIONAL :: real2
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   05/06/2001 dave
!!  MODIFICATION HISTORY
!!   14/06/2001 dave
!!    Created cq_abort_int and cq_abort_dble
!!   15/06/2001 dave
!!    Created cq_abort_none - no numbers passed
!!   15/06/2001 dave
!!    Merged back into a single routine with two optional arguments
!!   19/06/2001 dave
!!    Fixed formatting for second write
!!   19/07/2001 dave
!!    Reordered ifs so that BOTH integers are printed if present
!!   2006/07/06 08:04 dave
!!    Changed to take reals for new interface
!!  SOURCE
!!
  subroutine cq_abort_real(message,real1,real2)

    implicit none

    ! Passed variables
    character(len=*)       :: message
    real(double)           :: real1
    real(double), OPTIONAL :: real2

    ! Local variables
    integer :: ierror
    
    write (*,1)       myid+1
    write (io_lun,1)  myid+1
    if(PRESENT(real2)) then 
       write(*,3)      message, real1, real2
       write(io_lun,3) message, real1, real2
    else
       write(*,2)      message, real1
       write(io_lun,2) message, real1
    endif

1   format(2x,'Error in process ',i4)
2   format(2x,a,f20.12)
3   format(2x,a,2f20.12)

    call flush(io_lun)
    call flush(warning_lun)
    call MPI_abort( MPI_comm_world, 1, ierror )
    stop

    return
  end subroutine cq_abort_real
!!***

! -----------------------------------------------------------
! Subroutine cq_warn_none
! -----------------------------------------------------------

!!****f* GenComms/cq_warn_none *
!!
!!  NAME 
!!   cq_warn_none - writes out a warning message to file
!!  USAGE
!!   cq_warn_none(message indicating routine and problem)
!!  PURPOSE
!!   Allows writing of warning messages to separate file
!!   output file  
!!  INPUTS
!!   character(len=80) :: message - the abort message
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2019/12/05 08:21 dave
!!  MODIFICATION HISTORY
!!   2020/04/21 15:05 dave
!!    Update so that warnings only go to Conquest_warnings at iprint=0
!!  SOURCE
!!
  subroutine cq_warn_none(sub_name,message)

    implicit none

    ! Passed variables
    character(len=*) :: sub_name,message

    ! Local variables
    
    if(inode==ionode) then
       if(iprint>0) write(io_lun,fmt='(/"WARNING: ",a/)')  message
       write(warning_lun,fmt='(a,": ",a)')  trim(sub_name), message
       flag_warnings = .true.
    end if
    return
  end subroutine cq_warn_none
!!***

! -----------------------------------------------------------
! Subroutine cq_warn_int
! -----------------------------------------------------------

!!****f* GenComms/cq_warn_int *
!!
!!  NAME 
!!   cq_warn_int - writes out a warning message to file
!!  USAGE
!!   cq_warn_int(message indicating routine and problem, int1, int2)
!!  PURPOSE
!!   Allows writing of warning messages to separate file
!!   output file  
!!  INPUTS
!!   character(len=80) :: message - the abort message
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2019/12/05 08:21 dave
!!  MODIFICATION HISTORY
!!   2020/04/21 15:05 dave
!!    Update so that warnings only go to Conquest_warnings at iprint=0
!!  SOURCE
!!
  subroutine cq_warn_int(sub_name,message, int1, int2)

    implicit none

    ! Passed variables
    character(len=*) :: sub_name, message
    integer           :: int1
    integer, optional :: int2

    ! Local variables
    
    if(inode==ionode) then
       if(present(int2)) then
          if(iprint>0) write(io_lun,fmt='(/"WARNING: ",a,2i8/)')  message, int1, int2
          write(warning_lun,fmt='(a,": ",a,2i8)')  trim(sub_name), message, int1, int2
       else
          if(iprint>0) write(io_lun,fmt='(/"WARNING: ",a,i8/)')  message, int1
          write(warning_lun,fmt='(a,": ",a,i8)')  trim(sub_name), message, int1
       end if
       flag_warnings = .true.
    end if
    return
  end subroutine cq_warn_int
!!***

! -----------------------------------------------------------
! Subroutine cq_warn_real
! -----------------------------------------------------------

!!****f* GenComms/cq_warn_real *
!!
!!  NAME 
!!   cq_warn_real - writes out a warning message to file
!!  USAGE
!!   cq_warn_real(message indicating routine and problem, real1, real2)
!!  PURPOSE
!!   Allows writing of warning messages to separate file
!!   output file  
!!  INPUTS
!!   character(len=80) :: message - the abort message
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2019/12/05 08:21 dave
!!  MODIFICATION HISTORY
!!   2020/04/21 15:05 dave
!!    Update so that warnings only go to Conquest_warnings at iprint=0
!!  SOURCE
!!
  subroutine cq_warn_real(sub_name, message, real1, real2)

    implicit none

    ! Passed variables
    character(len=*) :: sub_name, message
    real(double)           :: real1
    real(double), optional :: real2

    ! Local variables
    
    if(inode==ionode) then
       if(present(real2)) then
          if(iprint>0) write(io_lun,fmt='(/"WARNING: ",a,2f20.12/)')  message, real1, real2
          write(warning_lun,fmt='(a,": ",a,2f20.12)')  trim(sub_name), message, real1, real2
       else
          if(iprint>0) write(io_lun,fmt='(/"WARNING: ",a,f20.12/)')  message, real1
          write(warning_lun,fmt='(a,": ",a,f20.12)')  trim(sub_name), message, real1
       end if
       flag_warnings = .true.
    end if
    return
  end subroutine cq_warn_real
!!***

! -----------------------------------------------------------
! Subroutine int_gsum
! -----------------------------------------------------------

!!****f* GenComms/int_gsum *
!!
!!  NAME 
!!   int_gsum - global sum on an integer
!!  USAGE
!!   int_gsum(variable)
!!  PURPOSE
!!   Performs a global sum on a specified integer - in
!!   other words, sum variable over all nodes and ensure
!!   that all processors get the results
!!  INPUTS
!!   integer :: variable
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   05/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine int_gsum(variable)

    implicit none

    ! Passed variables
    integer :: variable

    ! Local variables
    integer :: temp, ierr

    t1 = MPI_wtime()
    temp = variable
    call MPI_allreduce(temp, variable, 1, MPI_Integer, &
         MPI_sum, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_Success) then
       call cq_abort('int_gsum: Problem with allreduce')
    endif
    t2 = MPI_wtime()
    timings(1) = timings(1) + t2 - t1
    return
  end subroutine int_gsum
!!***

! -----------------------------------------------------------
! Subroutine logical_gsum
! -----------------------------------------------------------

!!****f* GenComms/logical_gsum *
!!
!!  NAME 
!!   logical_gsum - global sum on a logical
!!  USAGE
!!   logical_gsum(variable)
!!  PURPOSE
!!   Performs a global sum on a specified logical - in
!!   other words, OR variable over all processors and ensure
!!   that all processors get the results
!!  INPUTS
!!   logical :: variable
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   12:07, 26/02/2003 drb 
!!  MODIFICATION HISTORY
!!   12:16, 17/10/2005 drb 
!!    Bug fix from TM: should have had MPI_logical below
!!  SOURCE
!!
  subroutine logical_gsum(variable)

    implicit none

    ! Passed variables
    logical :: variable

    ! Local variables
    logical :: temp
    integer :: ierr

    t1 = MPI_wtime()
    temp = variable
    call MPI_allreduce(temp, variable, 1, MPI_logical, &
         MPI_LOR, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_Success) then
       call cq_abort('int_gsum: Problem with allreduce')
    endif
    t2 = MPI_wtime()
    timings(1) = timings(1) + t2 - t1
    return
  end subroutine logical_gsum
!!***

!!****f* GenComms/double_gsum *
!!
!!  NAME 
!!   double_gsum - global sum on a real(double)
!!  USAGE
!!   double_gsum(variable)
!!  PURPOSE
!!   Performs a global sum on a specified real(double) - in
!!   other words, sum variable over all nodes and ensure
!!   that all processors get the results
!!  INPUTS
!!   real(double) :: variable
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   05/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine double_gsum(variable)

    use datatypes

    implicit none

    ! Passed variables
    real(double) :: variable

    ! Local variables
    integer :: ierr

    real(double) :: temp

    t1 = MPI_wtime()
    temp = variable
    call MPI_allreduce(temp, variable, 1, MPI_Double_precision, &
         MPI_sum, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_Success) then
       call cq_abort('int_gsum: Problem with allreduce')
    endif
    t2 = MPI_wtime()
    timings(2) = timings(2) + t2 - t1
    return
  end subroutine double_gsum
!!***

!!****f* GenComms/dcplx_gsum *
!!
!!  NAME 
!!   dcplx_gsum - global sum on a complex(double_cplx)
!!  USAGE
!!   dcplx_gsum(variable)
!!  PURPOSE
!!   Performs a global sum on a specified complex(double_cplx) - in
!!   other words, sum variable over all nodes and ensure
!!   that all processors get the results
!!  INPUTS
!!   complex(double_cplx) :: variable
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   05/06/2001 dave
!!  MODIFICATION HISTORY
!!   06/06/2001 dave
!!    Minor alterations to headers above
!!  SOURCE
!!
  subroutine dcplx_gsum(variable)

    use datatypes

    implicit none

    ! Passed variables
    complex(double_cplx) :: variable

    ! Local variables
    integer :: ierr

    complex(double_cplx) :: temp

    t1 = MPI_wtime()
    temp = variable
    call MPI_allreduce(temp, variable, 1, MPI_Double_complex, &
         MPI_sum, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_Success) then
       call cq_abort('int_gsum: Problem with allreduce')
    endif
    t2 = MPI_wtime()
    timings(3) = timings(3) + t2 - t1
    return
  end subroutine dcplx_gsum
!!***

!!****f* GenComms/int_one_gsumv *
!!
!!  NAME 
!!   int_one_gsumv - global sum on an integer vector
!!  USAGE
!!   int_one_gsumv(variable,len)
!!  PURPOSE
!!   Performs a global sum on a specified integer 1D vector - in
!!   other words, sum variable over all nodes and ensure
!!   that all processors get the results
!!
!!   Since we don't know in advance how big this will be, we'll
!!   have a parameter in the module header int_gsum_buf which 
!!   will allow the sum to be done piecewise (for memory conservation)
!!
!!   The shenanigans with temp are necessary as there's no in-place
!!   global sum in MPI
!!  INPUTS
!!   integer :: len
!!   integer, dimension(len) :: variable
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   05/06/2001 dave
!!  MODIFICATION HISTORY
!!   06/06/2001 dave
!!    Declared j
!!    Corrected minor errors in headers above
!!  SOURCE
!!
  subroutine int_one_gsumv(variable,len)

    implicit none

    ! Passed variables
    integer :: len
    integer, dimension(len) :: variable

    ! Local variables
    integer, allocatable, dimension(:) :: temp
    integer :: ierr, stat, tmpsize, j

    t1 = MPI_wtime()
    if(len<=0) call cq_abort('int_gsumv: length of vector passed as zero')
    ! Establish size of temporary variable and allocate
    tmpsize = min(len,int_gsum_buf)
    allocate(temp(tmpsize),STAT=stat)
    if(stat/=0) call cq_abort('int_gsumv: error allocating temp')
    ! Loop over chunks
    do j=1,len,int_gsum_buf
       if(j+int_gsum_buf-1>len) tmpsize = len-j+1 ! Catch last bit
       ! Copy data into temp
       temp(1:tmpsize) = variable(j:j+tmpsize-1)
       call MPI_allreduce(temp,variable(j),tmpsize,MPI_Integer,&
            MPI_Sum,MPI_COMM_WORLD,ierr)
       if(ierr/=MPI_Success) call cq_abort('int_gsumv: Problem with allreduce')
    enddo
    ! Tidy up
    deallocate(temp,STAT=stat)
    if(stat/=0) call cq_abort('int_gsumv: error deallocating temp')
    t2 = MPI_wtime()
    timings(4) = timings(4) + t2 - t1
    return
  end subroutine int_one_gsumv
!!***

!!****f* GenComms/double_one_gsumv *
!!
!!  NAME 
!!   double_one_gsumv - global sum on an integer vector
!!  USAGE
!!   double_one_gsumv(variable,len)
!!  PURPOSE
!!   Performs a global sum on a specified real(double) 1D vector - in
!!   other words, sum variable over all nodes and ensure
!!   that all processors get the results
!!
!!   Since we don't know in advance how big this will be, we'll
!!   have a parameter in the module header int_gsum_buf which 
!!   will allow the sum to be done piecewise (for memory conservation)
!!
!!   The shenanigans with temp are necessary as there's no in-place
!!   global sum in MPI
!!  INPUTS
!!   integer :: len
!!   real(double), dimension(len) :: variable
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   05/06/2001 dave
!!  MODIFICATION HISTORY
!!   06/06/2001 dave
!!    Declared j
!!    Corrected minor errors in headers above
!!   05/07/2001 dave
!!    Corrected header to refer to double not int !
!!  SOURCE
!!
  subroutine double_one_gsumv(variable,len)

    implicit none

    ! Passed variables
    integer :: len
    real(double), dimension(len) :: variable

    ! Local variables
    real(double), allocatable, dimension(:) :: temp
    integer :: ierr, stat, tmpsize, j

    t1 = MPI_wtime()
    if(len<=0) call cq_abort('int_gsumv: length of vector passed as zero')
    ! Establish size of temporary variable and allocate
    tmpsize = min(len,double_gsum_buf)
    allocate(temp(tmpsize),STAT=stat)
    if(stat/=0) call cq_abort('int_gsumv: error allocating temp')
    ! Loop over chunks
    do j=1,len,double_gsum_buf
       if(j+double_gsum_buf-1>len) tmpsize = len-j+1 ! Catch last bit
       ! Copy data into temp
       temp(1:tmpsize) = variable(j:j+tmpsize-1)
       call MPI_allreduce(temp,variable(j),tmpsize,MPI_double_precision,&
            MPI_Sum,MPI_COMM_WORLD,ierr)
       if(ierr/=MPI_Success) call cq_abort('int_gsumv: Problem with allreduce')
    enddo
    ! Tidy up
    deallocate(temp,STAT=stat)
    if(stat/=0) call cq_abort('double_one_gsumv: error deallocating temp')
    t2 = MPI_wtime()
    timings(4) = timings(4) + t2 - t1
    return
  end subroutine double_one_gsumv
!!***

!!****f* GenComms/double_two_gsumv *
!!
!!  NAME 
!!   double_two_gsumv - global sum on a real(double) vector
!!  USAGE
!!   double_two_gsumv(variable, len1, len2)
!!  PURPOSE
!!   Performs a global sum on a specified real(double) 2D vector- in
!!   other words, sum vector over all processors and ensure
!!   that all processors get the results
!!
!!   Since we don't know in advance how big this will be, we'll
!!   have a parameter in the module header double_gsum_buf which 
!!   will allow the sum to be done piecewise (for memory conservation)
!!
!!   The shenanigans with temp are necessary as there's no in-place
!!   global sum in MPI
!!
!!   N.B. This routine (double_two_gsumv) assumes that len1 is smaller
!!   than double_gsum_buf - if not, things might get a little hairy.
!!  INPUTS
!!   integer :: len1, len2
!!   real(double), dimension(len1,len2) :: variable
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   05/06/2001 dave
!!  MODIFICATION HISTORY
!!   06/06/2001 dave
!!    Properly implemented
!!   08/06/2001 dave
!!    Fixed copy so that temp = variable (not the other way around !)
!!   05/12/2014 lat
!!    Fixed mod_gsum rounded to 0 
!!  SOURCE
!!
  subroutine double_two_gsumv(variable, len1, len2)

    use datatypes

    implicit none

    ! Passed variables
    integer :: len1, len2
    real(double), dimension(len1,len2) :: variable

    ! Local variables
    real(double), allocatable, dimension(:,:) :: temp
    integer :: ierr, stat, tmpsize, j, mod_dgsum_buf, mod_dgsum

    ! Make the buffersize an integer multiple of len1
    mod_dgsum = floor(double_gsum_buf/real(len1,double)) ! Round down
    ! mod_dgsum must be at least equal to 1
    if (mod_dgsum == 0) mod_dgsum = 1
    !print*, 'mod_dgsum', mod_dgsum, 'double_gsum_buf', double_gsum_buf, 'real(len1,double)', real(len1,double)
    mod_dgsum_buf = mod_dgsum*len1
    t1 = MPI_wtime()
    if(len1<=0.OR.len2<=0) &
         call cq_abort('double_gsumv: length of vector passed as zero')
    ! Establish size of temporary variable and allocate
    tmpsize = min(len2,mod_dgsum)
    allocate(temp(len1,tmpsize),STAT=stat)
    if(stat/=0) call cq_abort('double_gsumv: error allocating temp')
    ! Loop over chunks of len2
    do j=1,len2,mod_dgsum
       if(j+mod_dgsum-1>len2) tmpsize = len2-j+1 ! Catch last bit
       ! Copy data into temp
       call dcopy(len1*tmpsize,variable(1:len1,j:j+tmpsize-1),1,temp,1)
       call MPI_allreduce(temp,variable(1,j),len1*tmpsize,&
            MPI_Double_precision, MPI_Sum,MPI_COMM_WORLD,ierr)
       if(ierr/=MPI_Success) &
            call cq_abort('double_gsumv: Problem with allreduce')
    enddo
    ! Tidy up
    deallocate(temp,STAT=stat)
    if(stat/=0) call cq_abort('double_gsumv: error deallocating temp')
    t2 = MPI_wtime()
    timings(5) = timings(5) + t2 - t1
    return
  end subroutine double_two_gsumv
!!***

!!****f* GenComms/double_three_gsumv *
!!
!!  NAME 
!!   double_three_gsumv - global sum on a real(double) vector
!!  USAGE
!!   double_three_gsumv(variable, len1, len2, len3)
!!  PURPOSE
!!   Performs a global sum on a specified real(double) 3D vector- in
!!   other words, sum vector over all processors and ensure
!!   that all processors get the results
!!
!!   Since we don't know in advance how big this will be, we'll
!!   have a parameter in the module header double_gsum_buf which 
!!   will allow the sum to be done piecewise (for memory conservation)
!!
!!   The shenanigans with temp are necessary as there's no in-place
!!   global sum in MPI
!!
!!   N.B. This routine (double_two_gsumv) assumes that len1 is smaller
!!   than double_gsum_buf - if not, things might get a little hairy.
!!  INPUTS
!!   integer :: len1, len2, len3
!!   real(double), dimension(len1,len2,len3) :: variable
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   14:44, 27/07/2004 drb 
!!  MODIFICATION HISTORY
!! 
!!  SOURCE
!!
  subroutine double_three_gsumv(variable, len1, len2, len3)

    use datatypes

    implicit none

    ! Passed variables
    integer :: len1, len2, len3
    real(double), dimension(len1,len2,len3) :: variable

    ! Local variables
    real(double), allocatable, dimension(:,:,:) :: temp
    integer :: ierr, stat

    allocate(temp(len1,len2,len3),STAT=stat)
    if(stat/=0) call cq_abort('double_gsumv: error allocating temp')
    t1 = MPI_wtime()
    call dcopy(len1*len2*len3,variable,1,temp,1)
    call MPI_allreduce(temp,variable,len1*len2*len3,&
         MPI_Double_precision, MPI_Sum,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_Success) &
         call cq_abort('double_gsumv: Problem with allreduce')
    ! Tidy up
    deallocate(temp,STAT=stat)
    if(stat/=0) call cq_abort('double_gsumv: error deallocating temp')
    t2 = MPI_wtime()
    timings(5) = timings(5) + t2 - t1
    return
  end subroutine double_three_gsumv
!!***

!!****f* GenComms/double_four_gsumv *
!!
!!  NAME 
!!   double_four_gsumv - global sum on a real(double) vector
!!  USAGE
!!   double_four_gsumv(variable, len1, len2, len3, len4)
!!  PURPOSE
!!   Performs a global sum on a specified real(double) 3D vector- in
!!   other words, sum vector over all processors and ensure
!!   that all processors get the results
!!
!!   Since we don't know in advance how big this will be, we'll
!!   have a parameter in the module header double_gsum_buf which 
!!   will allow the sum to be done piecewise (for memory conservation)
!!
!!   The shenanigans with temp are necessary as there's no in-place
!!   global sum in MPI
!!
!!   N.B. This routine (double_two_gsumv) assumes that len1 is smaller
!!   than double_gsum_buf - if not, things might get a little hairy.
!!  INPUTS
!!   integer :: len1, len2, len3
!!   real(double), dimension(len1,len2,len3) :: variable
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2018/10/22 15:03 dave
!!  MODIFICATION HISTORY
!! 
!!  SOURCE
!!
  subroutine double_four_gsumv(variable, len1, len2, len3, len4)

    use datatypes

    implicit none

    ! Passed variables
    integer :: len1, len2, len3, len4
    real(double), dimension(len1,len2,len3, len4) :: variable

    ! Local variables
    real(double), allocatable, dimension(:,:,:,:) :: temp
    integer :: ierr, stat

    allocate(temp(len1,len2,len3,len4),STAT=stat)
    if(stat/=0) call cq_abort('double_four_gsumv: error allocating temp')
    t1 = MPI_wtime()
    call dcopy(len1*len2*len3*len4,variable,1,temp,1)
    call MPI_allreduce(temp,variable,len1*len2*len3*len4,&
         MPI_Double_precision, MPI_Sum,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_Success) &
         call cq_abort('double_gsumv: Problem with allreduce')
    ! Tidy up
    deallocate(temp,STAT=stat)
    if(stat/=0) call cq_abort('double_gsumv: error deallocating temp')
    t2 = MPI_wtime()
    timings(5) = timings(5) + t2 - t1
    return
  end subroutine double_four_gsumv
!!***

!!****f* GenComms/dcplx_one_gsumv *
!!
!!  NAME 
!!   dcplx_one_gsumv - global sum on a complex(double_cplx) 1D vector
!!  USAGE
!!   dcplx_one_gsumv(variable, len)
!!  PURPOSE
!!   Performs a global sum on a specified 1D complex(double_cplx) vector -
!!   in other words, sum vector over all processors and ensure
!!   that all processors get the results
!!
!!   Since we don't know in advance how big this will be, we'll
!!   have a parameter in the module header double_gsum_buf which 
!!   will allow the sum to be done piecewise (for memory conservation)
!!
!!   The shenanigans with temp are necessary as there's no in-place
!!   global sum in MPI
!!  INPUTS
!!   integer :: len
!!   complex(double_cplx), dimension(len) :: variable
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   05/06/2001 dave
!!  MODIFICATION HISTORY
!!   06/06/2001 dave
!!    Implemented properly
!!  SOURCE
!!
  subroutine dcplx_one_gsumv(variable, len)

    use datatypes

    implicit none

    ! Passed variables
    integer :: len
    complex(double_cplx), dimension(len) :: variable

    ! Local variables
    complex(double_cplx), allocatable, dimension(:) :: temp
    integer :: ierr, stat, tmpsize, j

    t1 = MPI_wtime()
    if(len<=0) call cq_abort('dcplx_gsumv: length of vector passed as zero')
    ! Establish size of temporary variable and allocate
    tmpsize = min(len,double_cplx_gsum_buf)
    allocate(temp(tmpsize),STAT=stat)
    if(stat/=0) call cq_abort('dcplx_gsumv: error allocating temp')
    ! Loop over chunks
    do j=1,len,double_cplx_gsum_buf
       if(j+double_cplx_gsum_buf-1>len) tmpsize = len-j+1 ! Catch last bit
       ! Copy data into temp
       temp(1:tmpsize) = variable(j:j+tmpsize-1)
       call MPI_allreduce(temp,variable(j),tmpsize,MPI_Double_complex,&
            MPI_Sum,MPI_COMM_WORLD,ierr)
       if(ierr/=MPI_Success) &
            call cq_abort('dcplx_gsumv: Problem with allreduce')
    enddo
    ! Tidy up
    deallocate(temp,STAT=stat)
    if(stat/=0) call cq_abort('dcplx_gsumv: error deallocating temp')
    t2 = MPI_wtime()
    timings(6) = timings(6) + t2 - t1
    return
  end subroutine dcplx_one_gsumv
!!***

!!****f* GenComms/double_gmax *
!!
!!  NAME 
!!   double_gmax - global max on a real(double)
!!  USAGE
!!   double_gmax(variable)
!!  PURPOSE
!!   Performs a global max on a specified real(double) - in
!!   other words, find the max value of the variable over all 
!!   processors and ensure
!!   that all processors get the results
!!  INPUTS
!!   real(double) :: variable
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   27/06/2001
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine double_gmax(variable)

    use datatypes

    implicit none

    ! Passed variables
    real(double) :: variable

    ! Local variables
    integer :: ierr

    real(double) :: temp

    t1 = MPI_wtime()
    temp = variable
    call MPI_allreduce(temp, variable, 1, MPI_Double_precision, &
         MPI_max, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_Success) then
       call cq_abort('double_gmax: Problem with allreduce')
    endif
    t2 = MPI_wtime()
    timings(2) = timings(2) + t2 - t1
    return
  end subroutine double_gmax
!!***

!!****f* GenComms/integer_gmax *
!!
!!  NAME 
!!   integer_gmax - global max on an integer
!!  USAGE
!!   integer_gmax(variable)
!!  PURPOSE
!!   Performs a global max on a specified integer - in
!!   other words, find the max value of the variable over all 
!!   processors and ensure
!!   that all processors get the results
!!  INPUTS
!!   integer :: variable
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2006/08/09
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine integer_gmax(variable)

    use datatypes

    implicit none

    ! Passed variables
    integer :: variable

    ! Local variables
    integer :: ierr

    integer :: temp

    t1 = MPI_wtime()
    temp = variable
    call MPI_allreduce(temp, variable, 1, MPI_Integer, &
         MPI_max, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_Success) then
       call cq_abort('integer_gmax: Problem with allreduce')
    endif
    t2 = MPI_wtime()
    timings(2) = timings(2) + t2 - t1
    return
  end subroutine integer_gmax
!!***

!!****f* GenComms/double_gmin *
!!
!!  NAME 
  !!   double_gmin - global min on a real(double)
!!  USAGE
!!   double_gmin(variable)
!!  PURPOSE
!!   Performs a global min on a specified real(double) - in
!!   other words, find the min value of the variable over all 
!!   processors and ensure
!!   that all processors get the results
!!  INPUTS
!!   real(double) :: variable
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   27/06/2001
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine double_gmin(variable)

    use datatypes

    implicit none

    ! Passed variables
    real(double) :: variable

    ! Local variables
    integer :: ierr

    real(double) :: temp

    t1 = MPI_wtime()
    temp = variable
    call MPI_allreduce(temp, variable, 1, MPI_Double_precision, &
         MPI_min, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_Success) then
       call cq_abort('double_gmin: Problem with allreduce')
    endif
    t2 = MPI_wtime()
    timings(2) = timings(2) + t2 - t1
    return
  end subroutine double_gmin
!!***

!!****f* GenComms/dcplx_exchv *
!!
!!  NAME 
!!   dcplx_exchv
!!  USAGE
!! 
!!  PURPOSE
!!   Complex global exchange routine.  Need to distribute chunks of 
!!   sendarray to all processors, to be put into recvarray.  The elements
!!   from psnd(j) to psnd(j+1)-1 should be sent to processor j and 
!!   prcv(j) to prcv(j+1)-1 in recv array will be received from processor
!!   j.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   06/06/2001
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine dcplx_exchv(sendarray, index_send, &
       recvarray, index_recv, size)

    use datatypes

    implicit none

    ! Passed variables
    integer :: size
    ! Original integer, dimension(numprocs+1) :: index_send,index_recv
    integer, dimension(:) :: index_send,index_recv
    complex(double_cplx), dimension(size) :: sendarray, recvarray

    ! Local variables
    integer :: i, ierror
    integer :: sndcounts(numprocs), rcvcounts(numprocs)
    integer :: snddisps(numprocs), rcvdisps(numprocs)

    t1 = MPI_wtime()

    ! Set up arguments required for MPI_alltoallv. Note that the
    ! send and receive displacement arguments to this routine are
    ! offsets from the start of the arrays, starting at zero (not 
    ! array indices starting at 1).
    do i = 1, numprocs
       sndcounts(i) = index_send(i+1) - index_send(i)
       rcvcounts(i) = index_recv(i+1) - index_recv(i)
       snddisps(i) = index_send(i) - 1
       rcvdisps(i) = index_recv(i) - 1
    enddo

    call MPI_alltoallv( sendarray, sndcounts, snddisps, MPI_double_complex, &
         recvarray, rcvcounts, rcvdisps, MPI_double_complex, &
         MPI_comm_world, ierror )
    if ( ierror /= MPI_success ) &
         call cq_abort('dcplx_exch: error in alltoallv')
    t2 = MPI_wtime()
    timings(13) = timings(13) + (t2-t1)
    return
  end subroutine dcplx_exchv
!!***

!!****f* GenComms/int_exch *
!!
!!  NAME 
!!   int_exch
!!  USAGE
!!   int_exch(sendarray, recvarray, count)
!!  PURPOSE
!!   Performs a global exchange on an integer array without 
!!   flexibility for where data comes from or goes to
!! (for which exchv is needed).  Sends count elements of
!!   sendarray to every processor and places in recvarray
!!  INPUTS
!!   integer :: count
!!   integer, dimension(nodes*count) :: sendarray, recvarray
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine int_exch(sendarray, recvarray, count)

    ! Passed variables
    integer :: count
    integer, dimension(numprocs*count) :: sendarray, recvarray

    ! Local variables
    integer :: ierr

    ! Send count elements of sendarray to every processor and
    ! recv count elements into recvarray from every processor
    call MPI_alltoall(sendarray,count,MPI_Integer, &
         recvarray,count,MPI_Integer, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_Success) call cq_abort('integer_exch: error in MPI_alltoall')
    return
  end subroutine int_exch
!!***

!!****f* GenComms/int_copy *
!!
!!  NAME 
!!   int_copy
!!  USAGE
!!   int_copy(vector)
!!  PURPOSE
!!   Copies the integer vec from root processor to
!!   all other processors via MPI bcast
!!  INPUTS
!!   integer :: vec
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   22/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine int_copy(vec)

    implicit none

    ! Passed variables
    integer :: vec

    ! Local variables
    integer :: ierror

    call MPI_bcast( vec, 1, MPI_integer, root, MPI_comm_world, ierror )
    if (ierror/=MPI_success) call cq_abort('int_copy: error in MPI_bcast')
    return
  end subroutine int_copy
!!***

!!****f* GenComms/int_one_copy *
!!
!!  NAME 
!!   int_one_copy
!!  USAGE
!!   int_one_copy(vector, length)
!!  PURPOSE
!!   Copies the 1D integer vector vec from root processor to
!!   all other processors via MPI bcast
!!  INPUTS
!!   integer :: length
!!   integer, dimension(length) :: vec
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine int_one_copy(vec, length)

    implicit none

    ! Passed variables
    integer :: length
    integer, dimension(length) :: vec

    ! Local variables
    integer :: ierror

    if (length<=0) call cq_abort('int_one_copy: invalid vector length')
    call MPI_bcast( vec, length, MPI_integer, root, MPI_comm_world, ierror )
    if (ierror/=MPI_success) call cq_abort('int_one_copy: error in MPI_bcast')
    return
  end subroutine int_one_copy
!!***

!!****f* GenComms/int_two_copy *
!!
!!  NAME 
!!   int_two_copy
!!  USAGE
!!   int_two_copy(vector, len1, len2)
!!  PURPOSE
!!   Copies the 2D integer vector vec from root processor to
!!   all other processors via MPI bcast
!!  INPUTS
!!   integer :: len1, len2
!!   integer, dimension(len1,len2) :: vec
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine int_two_copy(vec, len1, len2)

    implicit none

    ! Passed variables
    integer :: len1, len2
    integer, dimension(len1,len2) :: vec

    ! Local variables
    integer :: ierror

    if (len1<=0) call cq_abort('int_two_copy: invalid vector len1')
    if (len2<=0) call cq_abort('int_two_copy: invalid vector len2')
    call MPI_bcast( vec, len1*len2, MPI_integer, root, MPI_comm_world, ierror )
    if (ierror/=MPI_success) call cq_abort('int_two_copy: error in MPI_bcast')
    return
  end subroutine int_two_copy
!!***

!!****f* GenComms/double_copy *
!!
!!  NAME 
!!   double_copy
!!  USAGE
!!   double_copy(vector)
!!  PURPOSE
!!   Copies the real(double) vec from root processor to
!!   all other processors via MPI bcast
!!  INPUTS
!!   real(double) :: vec
!!  USES
!!   datatypes
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   22/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine double_copy(vec)

    use datatypes

    implicit none

    ! Passed variables
    real(double) :: vec

    ! Local variables
    integer :: ierror

    call MPI_bcast( vec, 1, MPI_Double_precision, root, &
         MPI_comm_world, ierror )
    if (ierror/=MPI_success) &
         call cq_abort('double_copy: error in MPI_bcast')
    return
  end subroutine double_copy
!!***

!!****f* GenComms/double_one_copy *
!!
!!  NAME 
!!   double_one_copy
!!  USAGE
!!   double_one_copy(vector, length)
!!  PURPOSE
!!   Copies the 1D real(double) vector vec from root processor to
!!   all other processors via MPI bcast
!!  INPUTS
!!   integer :: length
!!   real(double), dimension(length) :: vec
!!  USES
!!   datatypes
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine double_one_copy(vec, length)

    use datatypes

    implicit none

    ! Passed variables
    integer :: length
    real(double), dimension(length) :: vec

    ! Local variables
    integer :: ierror

    if (length<=0) call cq_abort('double_one_copy: invalid vector length')
    call MPI_bcast( vec, length, MPI_Double_precision, root, &
         MPI_comm_world, ierror )
    if (ierror/=MPI_success) &
         call cq_abort('double_one_copy: error in MPI_bcast')
    return
  end subroutine double_one_copy
!!***

!!****f* GenComms/double_two_copy *
!!
!!  NAME 
!!   double_two_copy
!!  USAGE
!!   double_two_copy(vector, len1, len2)
!!  PURPOSE
!!   Copies the 2D real(double) vector vec from root processor to
!!   all other processors via MPI bcast
!!  INPUTS
!!   integer :: len1, len2
!!   real(double), dimension(len1,len2) :: vec
!!  USES
!!   datatypes
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine double_two_copy(vec, len1, len2)

    use datatypes

    implicit none

    ! Passed variables
    integer :: len1, len2
    real(double), dimension(len1, len2) :: vec

    ! Local variables
    integer :: ierror

    if (len1<=0) call cq_abort('double_two_copy: invalid vector len1')
    if (len2<=0) call cq_abort('double_two_copy: invalid vector len2')
    call MPI_bcast( vec, len1*len2, MPI_Double_precision, root, &
         MPI_comm_world, ierror )
    if (ierror/=MPI_success) &
         call cq_abort('double_two_copy: error in MPI_bcast')
    return
  end subroutine double_two_copy
!!***

!!****f* GenComms/double_three_copy *
!!
!!  NAME 
!!   double_three_copy
!!  USAGE
!!   double_three_copy(vector, len1, len2, len3)
!!  PURPOSE
!!   Copies the 3D real(double) vector vec from root processor to
!!   all other processors via MPI bcast
!!  INPUTS
!!   integer :: len1, len2, len3
!!   real(double), dimension(length) :: vec
!!  USES
!!   datatypes
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine double_three_copy(vec, len1, len2, len3)

    use datatypes

    implicit none

    ! Passed variables
    integer :: len1, len2, len3
    real(double), dimension(len1, len2, len3) :: vec

    ! Local variables
    integer :: ierror

    if (len1<=0) call cq_abort('double_three_copy: invalid vector len1')
    if (len2<=0) call cq_abort('double_three_copy: invalid vector len2')
    if (len3<=0) call cq_abort('double_three_copy: invalid vector len3')
    call MPI_bcast( vec, len1*len2*len3, MPI_Double_precision, root, &
         MPI_comm_world, ierror )
    if (ierror/=MPI_success) &
         call cq_abort('double_three_copy: error in MPI_bcast')
    return
  end subroutine double_three_copy
!!***


  !!****f* GenComms/double_four_copy *
  !!
  !!  NAME 
  !!   double_four_copy
  !!  USAGE
  !!   call double_four_copy(vector, len1, len2, len3)
  !!  PURPOSE
  !!   Copies the 4D real(double) vector vec from root processor to
  !!   all other processors via MPI bcast
  !!  INPUTS
  !!   integer :: len1, len2, len3, len4
  !!   real(double), dimension(length) :: vec
  !!  USES
  !!   datatypes
  !!  AUTHOR
  !!   L.Tong
  !!  CREATION DATE
  !!   2012/04/03 L.Tong
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine double_four_copy(vec, len1, len2, len3, len4)

    use datatypes

    implicit none

    ! Passed variables
    integer :: len1, len2, len3, len4
    real(double), dimension(len1,len2,len3,len4) :: vec

    ! Local variables
    integer :: ierror

    if (len1 <= 0) call cq_abort('double_three_copy: invalid vector len1')
    if (len2 <= 0) call cq_abort('double_three_copy: invalid vector len2')
    if (len3 <= 0) call cq_abort('double_three_copy: invalid vector len3')
    if (len4 <= 0) call cq_abort('double_three_copy: invalid vector len4')
    call MPI_bcast(vec, len1*len2*len3*len4, MPI_Double_precision, &
                   root, MPI_comm_world, ierror)
    if (ierror /= MPI_success) &
         call cq_abort('double_four_copy: error in MPI_bcast')
    return
  end subroutine double_four_copy
  !!***


!!****f* GenComms/char_one_copy *
!!
!!  NAME 
!!   char_one_copy
!!  USAGE
!!   char_one_copy(vector, length)
!!  PURPOSE
!!   Copies the 1D character vector vec from root processor to
!!   all other processors via MPI bcast
!!  INPUTS
!!   integer :: length
!!   character(len=length) :: vec
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   22/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine char_one_copy(vec, length)

    implicit none

    ! Passed variables
    integer :: length
    character(len=length) :: vec

    ! Local variables
    integer :: ierror

    if (length<=0) call cq_abort('int_one_copy: invalid vector length')
    call MPI_bcast( vec, length, MPI_character, root, MPI_comm_world, ierror )
    if (ierror/=MPI_success) call cq_abort('int_one_copy: error in MPI_bcast')
    return
  end subroutine char_one_copy
!!***

!!****f* GenComms/logical_copy *
!!
!!  NAME 
!!   logical_copy
!!  USAGE
!!   logical_copy(vector)
!!  PURPOSE
!!   Copies the logical variable vec from root processor to
!!   all other processors via MPI bcast
!!  INPUTS
!!   integer :: length
!!   logical, dimension(length) :: vec
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   04/07/2001
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine logical_copy(vec)

    implicit none

    ! Passed variables
    logical :: vec

    ! Local variables
    integer :: ierror

    call MPI_bcast( vec, 1, MPI_logical, root, MPI_comm_world, ierror )
    if (ierror/=MPI_success) call cq_abort('logical_copy: error in MPI_bcast')
    return
  end subroutine logical_copy
!!***

!!****f* GenComms/logical_one_copy *
!!
!!  NAME 
!!   logical_one_copy
!!  USAGE
!!   logical_one_copy(vector, length)
!!  PURPOSE
!!   Copies the 1D logical vector vec from root processor to
!!   all other processors via MPI bcast
!!  INPUTS
!!   integer :: length
!!   logical, dimension(length) :: vec
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   22/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine logical_one_copy(vec, length)

    implicit none

    ! Passed variables
    integer :: length
    logical, dimension(length) :: vec

    ! Local variables
    integer :: ierror

    if (length<=0) call cq_abort('logical_one_copy: invalid vector length')
    call MPI_bcast( vec, length, MPI_logical, root, MPI_comm_world, ierror )
    if (ierror/=MPI_success) call cq_abort('logical_one_copy: error in MPI_bcast')
    return
  end subroutine logical_one_copy
!!***

!!****f* GenComms/logical_two_copy *
!!
!!  NAME 
!!   logical_two
!!  USAGE
!!   logical_two(vector, length1, length2)
!!  PURPOSE
!!   Copies the 2D logical vector vec from root processor to
!!   all other processors via MPI bcast
!!  INPUTS
!!   integer :: length
!!   logical, dimension(length1, length2) :: vec
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   22/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine logical_two_copy(vec, length1, length2)

    implicit none

    ! Passed variables
    integer :: length1, length2
    logical, dimension(length1, length2) :: vec

    ! Local variables
    integer :: ierror

    if (length1<=0) call cq_abort('logical_one_copy: invalid vector length')
    call MPI_bcast( vec, length1*length2, MPI_logical, root, MPI_comm_world, ierror )
    if (ierror/=MPI_success) call cq_abort('logical_one_copy: error in MPI_bcast')
    return
  end subroutine logical_two_copy
!!***

!!****f* GenComms/integer_copy_diff *
!!
!!  NAME 
!!   integer_copy_diff
!!  USAGE
!!   integer_copy_diff(target, source, node)
!!  PURPOSE
!!   Copies the integer target from processor node to
!!   all other processors via MPI bcast
!!  INPUTS
!!   integer :: node
!!   integer :: target
!!   integer :: source
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine integer_copy_diff(target, source, node)

    implicit none

    ! Passed variables
    integer :: node
    integer :: target
    integer :: source

    ! Local variables
    integer :: ierror

    if(node==inode) target = source
    call MPI_bcast( target, 1, MPI_integer, node-1, MPI_comm_world, ierror )
    if (ierror/=MPI_success) call cq_abort('int_one_copy: error in MPI_bcast')
    return
  end subroutine integer_copy_diff
!!***

!!****f* GenComms/int_one_copy_diff *
!!
!!  NAME 
!!   int_one_copy_diff
!!  USAGE
!!   int_one_copy_diff(vec_target, vec_source, node, length)
!!  PURPOSE
!!   Copies the 1D integer vector vec from processor node to
!!   all other processors via MPI bcast
!!  INPUTS
!!   integer :: length, node
!!   integer, dimension(length) :: vec_target
!!   integer, dimension(length) :: vec_source
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine int_one_copy_diff(vec_target, vec_source, node, length)

    implicit none

    ! Passed variables
    integer :: length, node
    integer, dimension(length) :: vec_target
    integer, dimension(length) :: vec_source

    ! Local variables
    integer :: ierror

    if (length<=0) call cq_abort('int_one_copy_diff: invalid vector length')
    if(node==inode) vec_target = vec_source
    call MPI_bcast( vec_target, length, MPI_integer, node-1, &
         MPI_comm_world, ierror )
    if (ierror/=MPI_success) call cq_abort('int_one_copy: error in MPI_bcast')
    return
  end subroutine int_one_copy_diff
!!***

!!****f* GenComms/my_barrier *
!!
!!  NAME 
!!   my_barrier
!!  USAGE
!! 
!!  PURPOSE
!!   Performs a barrier call across processors
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07/06/2001 dave
!!  MODIFICATION HISTORY
!!   2008/07/24 ast
!!    - Added timers
!!   2015/06/08 lat
!!    - Added experimental backtrace
!!  SOURCE
!!
  subroutine my_barrier()

    use global_module, ONLY: IPRINT_TIME_THRES3
    use timer_module,  ONLY: start_timer,stop_print_timer,stop_timer,WITH_LEVEL
    use timer_module,  ONLY: start_backtrace,stop_backtrace,cq_timer
    
    implicit none

    integer        :: ierr
    type(cq_timer) :: tmr_l_tmp1
    type(cq_timer) :: backtrace_timer

!****lat<$
    call start_backtrace(t=backtrace_timer,who='my_barrier',where=9)
!****lat>$

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    call MPI_barrier(MPI_COMM_WORLD, ierr)
    call stop_print_timer(tmr_l_tmp1,"my_barrier",IPRINT_TIME_THRES3)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='my_barrier')
!****lat>$

    if(ierr/=MPI_success) call cq_abort('my_barrier: problem with MPI_barrier')

    return
  end subroutine my_barrier
!!***

!!****f* GenComms/mtmini *
!!
!!  NAME 
!!   mtmini - Initialises mtime counter
!!  USAGE
!!   mtmini
!!  PURPOSE
!!   Initialises mtime counter
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07/06/2001 dave
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine mtmini()

    implicit none

    ! Record the start time
    secnd0 = MPI_wtime()

    return
  end subroutine mtmini
!!***

!!****f* GenComms/mtime *
!!
!!  NAME 
!!   mtime - timer
!!  USAGE
!!   mtime()
!!  PURPOSE
!!   Records elapsed time in milli-seconds using the
!!   MPI standard timing function
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07/06/2001 dave
!!  MODIFICATION HISTORY
!!   19/07/2001 dave
!!    Changed to be real
!!  SOURCE
!!
  real(double) function mtime()
    implicit none

    mtime = ( 1000.0_double * (MPI_wtime()-secnd0) )

    return
  end function mtime
!!***
end module GenComms
