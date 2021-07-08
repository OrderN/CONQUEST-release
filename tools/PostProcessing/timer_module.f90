! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id:$
! ------------------------------------------------------------------------------
! Module timer_module
! ------------------------------------------------------------------------------
! Code area 9: General
! ------------------------------------------------------------------------------

!!****h* Conquest/timer_module *
!!  NAME
!!   timer_module
!!
!!  USES
!!
!!  PURPOSE
!!   Timing functions
!!
!!  AUTHOR
!!   A.S.Torralba
!!
!!  CREATION DATE
!!   23/04/2008
!!
!!  MODIFICATION HISTORY
!!   2008/09/08 07:46 dave
!!    Changed variable types
!!   2008/09/08 12:11 dave & ast
!!    Added TimingOn flag
!!   2014/10/03 lat
!!    Added new variables in cq_timer for futur extensions 
!!    Added new timer subroutines
!! 
!!  TODO
!!
!!  SOURCE
!!
module timer_module

  use mpi
  use datatypes
  use global_module, only: zero
  use global_module, only: iprint, iprint_time

  implicit none

  ! Parameters
  logical, parameter :: TIME_ACCUMULATE_YES = .true.
  logical, parameter :: TIME_ACCUMULATE_NO  = .false.
  logical, parameter :: WITH_LEVEL          = .true.


  ! Flag to test for timers
  logical            :: TimingOn
  logical            :: TimersWriteOut ! This switches on output ON LOCAL PROCESSOR

  ! Flag to switch on/off backtracing
  logical            :: BackTraceOn

  ! Node information
  integer, save      :: mynode
  integer, save      :: myionode
  integer            :: lun_tmr
  integer            :: lun_user = 6

  ! Printing tolerance (the minimum time for a timer to be printed)
  real(double), save :: time_threshold

  ! Level indicator (basically, the number of timers started and not stopped)
  ! This is mainly intended for private use
  integer            :: cq_timer_level = 0

  
  type cq_timer
      logical :: have_ini         ! Signals whether an initial time mark has been obtained yet
      logical :: errors           ! True if errors, e.g. stop_timer without starting it
      logical :: first_use = .true.
      integer :: level
      integer :: calls
      real(double)      :: t_ini
      real(double)      :: t_end
      real(double)      :: t_tot
      integer           :: t_area  
      integer           :: t_level
      character(len=48)               :: t_name 
      character(len=1), dimension(12) :: t_index ! level visual indicator
  end type cq_timer
!!****

contains

! ------------------------------------------------------------------------------
! Subroutine init_timing_system
! ------------------------------------------------------------------------------

!!****f* H_matrix_module/init_timing_system *
!!
!!  NAME
!!   init_timing_system
!!  USAGE
!!
!!  PURPOSE
!!   To set information about the node 
!!    (direct access to GenComms creates circular dependence conflict)
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   A.S.Torralba
!!  CREATION DATE
!!   2008/07/24
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine init_timing_system(node)

    implicit none

    ! Passed parameters
    integer :: node

    mynode   = node 
    myionode = 1

  end subroutine init_timing_system
!!****


! ------------------------------------------------------------------------------
! Subroutine initialise_timer
! ------------------------------------------------------------------------------

!!****f* H_matrix_module/initialise_timer *
!!
!!  NAME 
!!   initialise_timer
!!  USAGE
!! 
!!  PURPOSE
!!   Initialises the timer
!!   NOTE: The user of timers NEEDS NOT call this explicitly, because
!!         start_timer will call it if it hasn't been
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   A.S.Torralba
!!  CREATION DATE
!!   23/04/08 
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine initialise_timer(t)

    implicit none

    ! Passed variables
    type(cq_timer)   :: t   ! The timer to be used

    t%first_use = .false.
    t%have_ini=.false.
    t%t_ini=0.0
    t%t_end=0.0
    t%t_tot=0.0
    t%errors=.false.

    return
  end subroutine initialise_timer
!!****

! ------------------------------------------------------------------------------
! Subroutine start_timer
! ------------------------------------------------------------------------------

!!****f* H_matrix_module/start_timer *
!!
!!  NAME 
!!   start_timer
!!  USAGE
!! 
!!  PURPOSE
!!   Gets initial time of the timer. This operation can (should) be reapeated after
!!   a final-time subroutine, such as stop_timer or stop_acummulate_timer, is called
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   A.S.Torralba
!!  CREATION DATE
!!   23/04/08 
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine start_timer(t,l)

    implicit none


    ! Passed variables
    type(cq_timer) :: t      ! The timer to be used
    logical,optional :: l    ! Do we want to assign a level? Ignored if already initialised

    if(.NOT.TimingOn) return
    if(t%first_use) then
       call initialise_timer(t)
    end if

    ! This subroutine is not, in principle, called again without calling a final-time routine,
    !   so, if it is, issue a warning
    if(t%have_ini.AND.TimersWriteOut.AND.iprint >= 3) then
       write(unit=lun_tmr,fmt='("start_timer: Warning: Second call in a row at node ", i3)') mynode
    else 
       t%have_ini=.true.   ! The timer is now started
       if(PRESENT(l)) then
         if(l.eqv.WITH_LEVEL) then
            t%level=cq_timer_level+1
            cq_timer_level=cq_timer_level+1
         endif
       else
         t%level=-1        ! This means the timer has no level (use mainly for standard total timers)
       endif
    end if

    ! In any case, get the time; for proper timing, this should be the last line in the subroutine
    !   so DON'T write any code after this
    t%t_ini = MPI_WTIME()

    return
  end subroutine start_timer
!!****

! ------------------------------------------------------------------------------
! Subroutine stop_timer
! ------------------------------------------------------------------------------

!!****f* H_matrix_module/stop_timer *
!!
!!  NAME 
!!   stop_timer
!!  USAGE
!! 
!!  PURPOSE
!!   Gets the time at the moment of the call, calculates the difference with the initial
!!   time, which MUST have been stored with start_timer, and stores the result in
!!   t_tot. It replaces its value, so if what you want is to get a total time, 
!!   call stop_accumulate_timer instead
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   A.S.Torralba
!!  CREATION DATE
!!   23/04/08 
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine stop_timer(t,accumulate)

    implicit none

    ! Passed variables
    type(cq_timer) :: t                ! The timer to be used
    logical, optional :: accumulate    ! Do we want to accumulate or restart

    ! This check is officially allowed ! 
    if(.NOT.TimingOn) return
    ! First things first: Get the time
    !   Don't EVER think of doing anything before this line
    !   You don't want to overestimate the time, do you?
    t%t_end = MPI_WTIME()

    ! Once that is done, check that:
    !  1) We have actually initialised the timer (if not, do it, but complain)
    !  2) An initial time mark is actually available (t%have_ini==.true.)
    ! If any of those are false, complain 
    if(.not.t%first_use .and. (t%have_ini .eqv. .true.)) then
      if(present(accumulate)) then
        if(accumulate.eqv.TIME_ACCUMULATE_NO) then  ! Don't accumulate times
          t%t_tot=t%t_end-t%t_ini
        else
          t%t_tot=t%t_tot+(t%t_end-t%t_ini) 
        end if
      else  
        t%t_tot=t%t_tot+(t%t_end-t%t_ini) 
      end if
      t%have_ini=.false.                                   ! But now we need to restart the timer again
      if(t%level > 0) then
         cq_timer_level = cq_timer_level - 1               ! And we go up one level in the hierarchy
      endif                                                !   so we have to ask for a level again, if needed
    else
      if(t%first_use)  call initialise_timer(t)
      t%errors=.true.
      if(TimersWriteOut.AND.iprint >= 3) then
         write(unit=lun_tmr,&
         fmt='("stop_timer: Error : Tried to calculate time difference without initial time mark &
        &at node", i3)') mynode
      end if
    end if

    return
  end subroutine stop_timer
!!****

! ------------------------------------------------------------------------------
! Subroutine print_timer
! ------------------------------------------------------------------------------

!!****f* timer_module/print_timer *
!!
!!  NAME 
!!   print_timer
!!  USAGE
!! 
!!  PURPOSE
!!   Prints the total time stored in a timer
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   A.S.Torralba
!!  CREATION DATE
!!   23/04/08 
!!  MODIFICATION HISTORY
!!   2014/01/18 lat
!!   - Added optional argument for printing
!!
!!  SOURCE
!!
  subroutine print_timer(t,m,lun)

    implicit none
    
    ! Passed variables
    type(cq_timer)    :: t            ! The timer to be used
    character(len=*)  :: m            ! Label for the printout
    character(len=40) :: message      ! Label for the printout
    integer, optional :: lun          ! To print somewhere else
    integer           :: lun_sav
    
    ! write(unit=lun_tmr,fmt='("Time resolution = ", e18.10e2)') MPI_WTICK()

    if( .not. TimingOn ) return
    if(  present(lun)  ) lun_sav = lun_tmr
    if(  present(lun)  ) lun_tmr = lun
    
    if( TimersWriteOut.AND.t%t_tot > time_threshold )  then
       message = trim(m)
       if( t%level >= 0 ) then             ! Print the level if assigned
          write(unit=lun_tmr,fmt='("Timing: Level ",i3," - Proc ",i6,": &
              &Time spent in ", a," = ", f12.5," s")') &
               t%level,mynode,message, t%t_tot
       else
          if( present(lun) ) lun_tmr = lun
          write(unit=lun_tmr,fmt='("Timing:                Proc ",i6,": &
              &Time spent in ", a," = ", f12.5," s")') &
               mynode, message, t%t_tot
       endif
    endif
    
    if( present(lun) ) lun_tmr = lun_sav
    
    return
  end subroutine print_timer
  !****
  
!!$  subroutine print_timer(t,m)
!!$    
!!$    implicit none
!!$
!!$    ! Passed variables
!!$    type(cq_timer) :: t               ! The timer to be used
!!$    character(len=*) :: m             ! Label for the printout
!!$    ! character(len=50) :: message      ! Label for the printout
!!$
!!$!    write(unit=lun_tmr,fmt='("Time resolution = ", e18.10e2)') MPI_WTICK()
!!$
!!$    if(.NOT.TimingOn) return
!!$    if(TimersWriteOut.AND.t%t_tot > time_threshold)  then
!!$      if(t%level >= 0) then             ! Print the level if assigned
!!$         write(unit=lun_tmr,fmt='("Timing: Level ",i3," - &
!!$              &Proc ",i6,": Time spent in ", a," = ", f12.5," s")') &
!!$              t%level, mynode, trim(m), t%t_tot
!!$      else
!!$         write(unit=lun_tmr,fmt='("Timing: Proc ",i6,": &
!!$              &Time spent in ", a," = ", f12.5," s")') mynode, trim(m), t%t_tot
!!$      endif
!!$   endif
!!$   
!!$   return
!!$  end subroutine print_timer
!!***

! ------------------------------------------------------------------------------
! Subroutine stop_print_timer
! ------------------------------------------------------------------------------

!!****f* timer_module/stop_print_timer *
!!
!!  NAME
!!   stop_print_timer
!!  USAGE
!!
!!  PURPOSE
!!   Stops the timer WITHOUT accumulation and prints the total time stored
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   A.S.Torralba
!!  CREATION DATE
!!   17/07/08
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine stop_print_timer(t,m,i,a)

    implicit none

    ! Passed variables
    type(cq_timer)   :: t             ! The timer to be used
    character(len=*) :: m             ! Label for the printout
    logical,optional :: a             ! Accumulate time (Default = Don't acummulate)
    integer :: i                      ! iprint_time level to be used

    ! Local variables
    logical :: accumulate

    if(.NOT.TimingOn) return
    accumulate = TIME_ACCUMULATE_NO

    if( present(a) ) accumulate = a

    call stop_timer(t, accumulate) 

    if( TimersWriteOut .and. iprint_time >= i ) then
       call print_timer(t,m)
    endif

    return
  end subroutine stop_print_timer
!!****


!****************************************************************************************
!
!  **<lat>** below new backtracing subroutines called at the begining/end of important 
!            subroutines mainly for debugging purpose ; also interesting to find out where you are
!            in the code during execution (to be improved)
!
!****************************************************************************************

!!****f* timer_module/init_backtrace *
!!
!!  NAME 
!!   initialise_timer_new
!!
!!  USAGE
!! 
!!  PURPOSE
!!
!!  INPUTS 
!! 
!!  USES
!! 
!!  AUTHOR
!!   A.S.Torralba/L.Truflandier
!!
!!  CREATION DATE
!!   2014/10/03
!!
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine initialise_backtrace(t,name,area,state)

    implicit none

    ! Passed variables
    type(cq_timer)              :: t    ! The timer to be used
    character(len=*),  optional :: name 
    integer,           optional :: area
    integer,           optional :: state

    ! Local variables
    character(len=1 ), parameter:: symb = '*'
    character(len=48)           :: string
    integer                     :: tmp

    t%first_use = .false.
    t%have_ini  = .false.
    t%errors    = .false.
    t%calls     =  0
    t%t_ini     =  zero
    t%t_end     =  zero
    t%t_tot     =  zero

    if ( present(name) ) then
       if ( mod(len(name),2) == 0 ) then
          tmp      = (48-2-len(name))/2
          string   = repeat(symb,tmp)//' '//name//' '//repeat(symb,tmp)
          t%t_name = string
       else
          tmp      = (48-2-len(name)-1)/2
          string   = repeat(symb,tmp)//' '//name//' '//repeat(symb,tmp+1)
          t%t_name = string
       end if
    else       
       string   = repeat(symb,19)//' unknown '//repeat(symb,20)
       t%t_name = string
    end if

    if ( present(area) ) then
       t%t_area = area
    else
       t%t_area = 0
    end if

    t%t_index = ' '
    !print*, 'state', state
    if ( present(state) ) then
       if ( state==0 ) then
          t%t_level = state
          t%t_index = '-'
          !
       elseif ( state>0 ) then
          t%t_level        = state
          t%t_index(state) =  '*'
          !
       end if
    else
       t%t_level = -1
       t%t_index = ' '
       !
    end if

    return
  end subroutine initialise_backtrace
!!****

!!****f* timer_module/start_backtrace *
!!
!!  NAME 
!!   start_backtrace
!!
!!  USAGE
!! 
!!  PURPOSE
!!
!!  INPUTS
!! 
!!  USES
!! 
!!  AUTHOR
!!   A.S.Torralba/L.A.Truflandier
!!
!!  CREATION DATE
!!   2014/10/03
!!
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine start_backtrace(t,l,who,where,level,echo)

    implicit none

    ! Passed variables
    type(cq_timer)             :: t     ! The timer to be used
    logical,          optional :: l     ! Do we want to assign a level ? 
                                        ! Ignored if already initialised
    character(len=*), optional :: who   ! what subroutine   ?
    integer,          optional :: where ! ... in which area ?
    integer,          optional :: level ! where in the code's tree ? 
    logical,          optional :: echo  ! print in std

    ! Local variables
    integer                    :: i
   
    if( .not. BackTraceOn) return

    if( t%first_use ) then       
       if ( present(level) ) then
          call initialise_backtrace(t,who,where,level)
       else
          call initialise_backtrace(t,who,where, -1  )
       end if
    end if
    
    ! This subroutine is not, in principle, called again without calling a 
    ! final-time routine, so, if it is, issue a warning
    if( t%have_ini .and. iprint >= 3 ) then
       write(unit=lun_user,fmt='("start_backtrace: Warning: Second call in a row for ",a, &
           &"at node ", i6)') who, mynode
       if( present(who) .and. present(echo) .and. mynode==myionode ) then
          write(lun_user,3) t%t_area, t%t_name, (t%t_index(i), i=1,12)
       end if
    else
       t%have_ini = .true.   ! The timer is now started
       if( present(l) ) then
          if( l.eqv.WITH_LEVEL ) then
             t%level       = cq_timer_level + 1
             cq_timer_level= cq_timer_level + 1
          endif
       else
          t%level = - 1 ! This means the timer has no level (use mainly for standard total timers)
       endif
       if( present(who) .and. present(echo) .and. mynode==myionode ) then
          write(lun_user,2) t%t_area, t%t_name, (t%t_index(i), i=1,12)
       end if
    end if

    ! In any case, get the time; for proper timing, this should be the last line in the subroutine
    !   so DON'T write any code after this
    t%t_ini = MPI_WTIME()

    return

2   format( '[Area:',i2,']',1x,'Backtrace  STAR:',1x,a48,1x,'[',12a,']' )
3   format( '[Area:',i2,']',1x,'Backtrace  WARN:',1x,a48,1x,'[',12a,']' )


  end subroutine start_backtrace
!!***

!!****f* timer_module/stop_backtrace *
!!
!!  NAME 
!!   stop_backtrace
!!
!!  USAGE
!! 
!!  PURPOSE
!!
!!  INPUTS 
!! 
!!  USES
!! 
!!  AUTHOR
!!   A.S.Torralba/L.A.Truflandier
!!
!!  CREATION DATE
!!   2014/10/03
!!
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine stop_backtrace(t,accumulate,who,echo)

    implicit none

    ! Passed variables
    type(cq_timer)    :: t           ! The timer to be used
    logical, optional :: accumulate  ! Do we want to accumulate or restart
    character(len=*), optional :: who
    logical,          optional :: echo

    ! Local variables
    integer           :: i

    ! This check is officially allowed ! 
    if( .not. BackTraceOn) return

    ! First things first: Get the time
    !   Don't EVER think of doing anything before this line
    !   You don't want to overestimate the time, do you?
    t%t_end = MPI_WTIME()

    ! Once that is done, check that:
    !  1) We have actually initialised the timer (if not, do it, but complain)
    !  2) An initial time mark is actually available (t%have_ini==.true.)
    !     If any of those are false, complain 
    if( .not. t%first_use .and. (t%have_ini .eqv. .true.) ) then
       if( present(accumulate) ) then
          if( accumulate.eqv.TIME_ACCUMULATE_NO ) then  ! Don't accumulate times
             t%t_tot=t%t_end-t%t_ini
          else
             t%t_tot=t%t_tot+(t%t_end-t%t_ini) 
          end if
          !if ( present(who) ) then
          !   write(lun_user,1) (t%t_index(i), i=1,6), t%t_name, t%t_area
          !end if
       else  
          t%t_tot=t%t_tot+(t%t_end-t%t_ini) 
          if ( present(who) .and. present(echo) .and. mynode==myionode ) then 
             write(lun_user,1) & 
                   t%t_area, t%t_tot, t%t_name, (t%t_index(i), i=1,12)
          end if
          t%have_ini = .false.  ! But now we need to restart the timer again
       end if
       if( t%level > 0 ) then
          cq_timer_level = cq_timer_level - 1 ! And we go up one level in the hierarchy
       endif ! so we have to ask for a level again, if needed
    else
       if( t%first_use ) call initialise_backtrace(t)
       t%errors = .true.
       if( iprint >= 3 ) then
          write(unit=lun_user,&
               fmt='("stop_backtrace: Warning: Tracing of ",a," not possible &
              &at node", i6)') who, mynode
       end if
       if ( present(who) .and. present(echo) .and. mynode==myionode ) then        
          write(lun_user,2) t%t_area, t%t_tot, t%t_name, (t%t_index(i), i=1,12)
       end if
    end if
    
    return
    
1   format( /'[Area:',i2,']',1x,'time = ',1x,f8.3,1x,'s'/ &
                          9x,1x,'Backtrace  STOP:',1x,a48,1x,'[',12a,']' )
2   format( /'[Area:',i2,']',1x,'time = ',1x,f8.3,1x,'s'/ &
                          9x,1x,'Backtrace  WARN:',1x,a48,1x,'[',12a,']' )
  end subroutine stop_backtrace
!!****

end module timer_module
