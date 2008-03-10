! $Id$
! -----------------------------------------------------------
! Module PosTan
! -----------------------------------------------------------
! Code area 9: general
! -----------------------------------------------------------

!!****h* Conquest/PosTan
!!  NAME
!!   PosTan
!!  PURPOSE
!!   Calculates the coefficients used to relate an energy 
!!   tolerance (i.e. dE) to a residual by a process of least
!!   squares fitting, subject to the fit being below all points.
!!
!!   Conquest makes extensive use of residual minimisation 
!!   techniques, and we need a way to relate a tolerance on the
!!   energy (which is easily specified by the user) to the
!!   magnitude of the residual -- hence this routine.
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   15/11/00
!!  MODIFICATION HISTORY
!!   22/05/2001 dave
!!    ROBODoc header
!!   21/06/2001 dave
!!    Added RCS Id and Log tags and removed stop statements.
!!    Created fit_coeff; bug fixes
!!   2008/02/01 17:48 dave
!!    Changes for output to file not stdout
!!  SOURCE
!!
module PosTan

  use datatypes
  use global_module, ONLY: io_lun

  implicit none

  integer, parameter :: max_iters = 50
  real(double), parameter :: eprec = 1.0e-12_double
  real(double), parameter :: cscale = 0.01_double

  real(double) :: PulayE(max_iters),PulayR(max_iters)
  real(double) :: SCE(max_iters),SCR(max_iters)
  real(double) :: SupFnE(max_iters),SupFnR(max_iters)
  real(double) :: PulayC, PulayBeta
  real(double) :: SCC,SCBeta
  real(double) :: SupFnC, SupFnBeta

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"

!!***

contains

!!****f* PosTan/fit_coeff *
!!
!!  NAME 
!!   fit_coeff - fits coefficients
!!  USAGE
!!   fit_coeff(C, beta, E, R, n)
!!  PURPOSE
!!   Fits the coefficients C and beta (such that the tolerance applied
!!   to the residual is C*dE^beta where dE is the energy tolerance) to
!!   the arrays E and R
!!  INPUTS
!!   integer :: n - size of arrays
!!   real(double), dimension(n) :: E, R - energy and residual at each iter
!!   real(double) :: C, beta - fitted coefficients
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   14/06/2001 dave
!!  MODIFICATION HISTORY
!!   2008/03/03 18:34 dave
!!    Removed dlog, dabs
!!  SOURCE
!!
  subroutine fit_coeff(C, beta, E, R, n)

    use datatypes
    use numbers
    use GenComms, ONLY: inode, ionode
    use global_module, ONLY: iprint_gen

    implicit none

    ! Passed variables
    integer, intent(IN) :: n
    real(double), dimension(n), intent(IN) :: E, R
    real(double), intent(OUT) :: C, beta

    ! Local variables
    integer :: nkeep, ndelta, i
    real(double) :: LastE
    real(double), dimension(n) :: SCE, SCR ! Automatic arrays

    if(n<=2) then
       nkeep = 0
    else
       ndelta = n-1
       LastE = E(n)
       do i=1,ndelta
          SCE(i) = abs(E(i)-LastE)
       enddo
       nkeep=0
       do i=1,ndelta
          if(SCE(i)>eprec*abs(LastE))then
             nkeep = nkeep+1
             SCR(i) = log(R(i))
             SCE(i) = log(SCE(i))
          endif
       enddo
    endif
    if(nkeep<2) then
       if(inode==ionode.AND.iprint_gen>1) write(io_lun,*) '  No energy deviations kept'
       C = zero
       beta = zero
    else
       if(inode==ionode.AND.iprint_gen>1) then
          write(io_lun,*) '  Logs of residual and energy deviations'
          do i=1,nkeep
             write(io_lun,1) i,SCR(i),SCE(i)
          enddo
       endif
       call pos_tan(n,nkeep,SCE,SCR,C,beta, inode-1)
       C = cscale*dexp(C)
    endif
1   format(2x,i4,2f15.8)
    return
  end subroutine fit_coeff
!!***

!!****f* PosTan/pos_tan *
!!
!!  NAME 
!!   pos_tan
!!  USAGE
!! 
!!  PURPOSE
!!   This subroutine performs the fitting described above. 
!!  INPUTS
!!   integer :: mxnpt, npt   (max. number of points, number of points)
!!   real(double) :: xpy, ypt (energy and residual at each iteration)
!!   real(double) :: cconst, bslope (fitted parameters)
!!  USES
!!   datatypes, global_module, GenComms
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   15/11/00
!!  MODIFICATION HISTORY
!!  22/05/2001 dave
!!   ROBODoc header
!!  SOURCE
!!
  subroutine pos_tan(mxnpt,npt,xpt,ypt,cconst,bslope,myid)

    ! Module usage
    use datatypes
    use global_module, ONLY: iprint_gen
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer, intent(in) :: mxnpt,npt, myid
    real(double) :: cconst, bslope
    real(double) :: xpt(mxnpt),ypt(mxnpt)

    ! Local variables
    ! Automatic arrays
    real(double) :: slopmin(mxnpt),slopmax(mxnpt)
    integer :: icand(mxnpt)
    integer :: indx(mxnpt)
    ! Other
    integer :: i,ir,il,n,nr,nptmin,io,imin,ip
    real(double) :: slope,alpha,beta,gamma
    real(double) :: bmin, sumsq, bminmin,sumsq1

    ! --- test array bounds ---------------------------------------------
    if((npt.lt.2).or.(npt.gt.mxnpt)) then
       call cq_abort('pos_tan: no. of points out of range',npt)
    endif
    ! --- make index for ordering of points from right to left ----------
    call picksort(mxnpt,npt,xpt,indx)
    if(myid==0.AND.iprint_gen>1) then 
       write(io_lun,15)
15     format(/'sequence nos. of points, going from right to left:'/)
       do i=1,npt
          write(io_lun,16) i,indx(i)
16        format(5x,i5,2x,i5)
       enddo
    endif
    ! --- test candidacy of points --------------------------------------
    ! --- right-most point: test slopes to pts on left ------------------
    ir=indx(1)
    do nr=2,npt
       il=indx(nr)
       slope=(ypt(ir)-ypt(il))/(xpt(ir)-xpt(il))
       if(nr.eq.2) then
          slopmin(1)=slope
       else
          if(slope.gt.slopmin(1)) slopmin(1)=slope
       endif
    enddo
    if(slopmin(1).lt.0.0_double) slopmin(1)=0.0_double
    icand(1)=1
    ! --- other points --------------------------------------------------
    do n=2,npt
       icand(n)=0
       il=indx(n)
       ! --- test slopes to pts on right -----------------------------------
       do nr=1,n-1
          ir=indx(nr)
          slope=(ypt(ir)-ypt(il))/(xpt(ir)-xpt(il))
          if(nr.eq.1) then
             slopmax(n)=slope
          else
             if(slope.lt.slopmax(n)) slopmax(n)=slope
          endif
       enddo
       ! --- test slopes to pts on left -----------------------------------
       if(slopmax(n).ge.0.) then
          if(n.eq.npt) then
             icand(n)=1
             slopmin(n)=0.
          else
             ir=indx(n)
             do nr=n+1,npt
                il=indx(nr)
                slope=(ypt(ir)-ypt(il))/(xpt(ir)-xpt(il))
                if(nr.eq.n+1) then
                   slopmin(n)=slope
                else
                   if(slope.gt.slopmin(n)) slopmin(n)=slope
                endif
             enddo
             if(slopmin(n).lt.0.) slopmin(n)=0.
             if(slopmin(n).le.slopmax(n)) icand(n)=1 
          endif
       endif
    enddo
    if(myid==0.AND.iprint_gen>1) then
       write(io_lun,21)
21     format(/'surviving candidates:'/)
       do n=1,npt
          write(io_lun,22) n,icand(n)
22        format(i5,2x,i3)
       enddo
       write(io_lun,23)
23     format(/'min. and max. slopes:'/)
       n=1
       write(io_lun,24) n,slopmin(n)
24     format(i5,2x,e15.6)
       do n=2,npt
          if(icand(n).eq.1) then
             write(io_lun,25) n,slopmin(n),slopmax(n)
25           format(i5,2x,2e15.6)
          endif
       enddo
    endif
    ! --- find point giving smallest sum of squares -----------------------
    ! --- right-most point ------------------------------------------------
    alpha=0.0_double
    beta=0.0_double
    gamma=0.0_double
    ir=indx(1)
    do nr=2,npt
       il=indx(nr)
       alpha=alpha+(xpt(il)-xpt(ir))**2
       beta=beta+(ypt(il)-ypt(ir))*(xpt(il)-xpt(ir))
       gamma=gamma+(ypt(il)-ypt(ir))**2
    enddo
    bmin=beta/alpha
    if(bmin.lt.slopmin(1)) bmin=slopmin(1)
    sumsq=alpha*bmin*bmin-2.*beta*bmin+gamma
    nptmin=1
    bminmin=bmin
    ! --- other points ---------------------------------------------------
    do n=2,npt
       if(icand(n).eq.1) then
          ip=indx(n)
          alpha=0.0_double
          beta=0.0_double
          gamma=0.0_double
          do nr=1,n-1
             io=indx(nr)
             alpha=alpha+(xpt(io)-xpt(ip))**2
             beta=beta+(ypt(io)-ypt(ip))*(xpt(io)-xpt(ip))
             gamma=gamma+(ypt(io)-ypt(ip))**2
          enddo
          if(n.lt.npt) then
             do nr=n+1,npt
                io=indx(nr)
                alpha=alpha+(xpt(io)-xpt(ip))**2
                beta=beta+(ypt(io)-ypt(ip))*(xpt(io)-xpt(ip))
                gamma=gamma+(ypt(io)-ypt(ip))**2
             enddo
          endif
          bmin=beta/alpha
          if(bmin.lt.slopmin(n)) bmin=slopmin(n)
          if(bmin.gt.slopmax(n)) bmin=slopmax(n)
          sumsq1=alpha*bmin*bmin-2.*beta*bmin+gamma
          if(sumsq1.lt.sumsq) then
             sumsq=sumsq1
             nptmin=n
             bminmin=bmin
          endif
       endif
    enddo
    if(myid==0.AND.iprint_gen>1) write(io_lun,31) nptmin,bminmin,sumsq
31  format(/'point with smallest sum of squares:',i5/ &
         'slope with smallest sum of squares:',e15.6/ &
         'sum of squares:',e15.6)
    imin=indx(nptmin)
    cconst=ypt(imin)-bminmin*xpt(imin)
    bslope=bminmin
    return 
  end subroutine pos_tan
!!***

!!****f* PosTan/picksort *
!!
!!  NAME 
!!   picksort
!!  USAGE
!! 
!!  PURPOSE
!!   Sorts elements of an array into descending order - makes
!!   an index array indx so that arr(indx(i)) is correctly 
!!   ordered.
!!   Taken from Numerical Recipes
!!  INPUTS
!!   integer :: n, nmax - points, array size
!!   integer, dimension(nmax) :: indx - index giving sorted array
!!   real(double), dimension(nmax) :: arr - array to be sorted
!!  USES
!!   datatypes, GenComms
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   15/11/00
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Changed stop to cq_abort and added use GenComms
!!    Added to documentation a little
!!  SOURCE
!!
  subroutine picksort(nmax,n,arr,indx)

    use datatypes
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer :: n,nmax
    integer, dimension(nmax) :: indx
    real(double), dimension(nmax) :: arr

    ! Local variables
    integer :: i,j,ifind
    real(double) :: a

    if((n<2).OR.(n>nmax)) then
       call cq_abort('picksort: no. of elements out of range',n)
    endif
    ! Initiate index array
    do i=1,n
       indx(i)=i
    enddo
    do j=2,n ! Now sort array
       a=arr(j)
       ifind=0
       i=j
       do while((ifind==0).AND.(i>1))
          i=i-1
          if(arr(indx(i)).lt.a) then
             indx(i+1)=indx(i)
          else
             indx(i+1)=j
             ifind=1
          endif
       enddo ! End do while(ifind==0.AND.i>1)
       if(ifind==0) then
          indx(1)=j
       endif
    enddo ! End do j=2,n
    return
  end subroutine picksort
!!***
end module PosTan
