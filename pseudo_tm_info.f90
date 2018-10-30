! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module pseudo_tm_info
! ------------------------------------------------------------------------------
! Code area 10: Pseudopotentials
! ------------------------------------------------------------------------------

!!****h* Conquest/pseudo_tm_info 
!!
!!  NAME
!!   pseudo_tm_info
!!  PURPOSE
!!    Sets up and stores the information of pseudopotentials for 
!!   each species of atom.
!!    For using TM's pseudopotetials, 
!!
!!  USES
!!   datatypes, GenComms, numbers, spline_module, global_module
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   18/06/2002
!!  MODIFICATION HISTORY
!!   2008/02/06 08:32 dave
!!    Changed for output to file not stdout
!!   2008/05/25
!!    Added timers
!!   2011/03/30 19:07 M.Arita
!!    Added statements for P.C.C.
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2015/11/09 17:21 dave with TM, NW (Mizuho)
!!    Adding new members to pseudo_info derived type for neutral atom
!!   2017/11/27 15:52 dave with TM, NW
!!    NA projectors added to pseudo_info type
!!  SOURCE
!!
module pseudo_tm_info

  use datatypes,              only: double
  use GenComms,               only: cq_abort
  use global_module,          only: area_pseudo, io_lun, flag_pcc_global
  use timer_module,           only: start_timer,stop_timer
  use timer_stdclocks_module, only: tmr_std_allocation

  implicit none

  save
  private
  public :: rad_func, pseudo_info, pseudo, setup_pseudo_info, &
           rad_alloc, rad_dealloc, loc_pot, loc_chg

  ! Choose between representing local part of pseudopotential as a local potential or charge density
  integer :: loc_pot = 1
  integer :: loc_chg = 2

  !type rad_func is same as in TM
  !  see radial.f of SIESTA
  type rad_func
     integer      :: n               ! number of grids along radial direction
     real(double) :: cutoff          ! cutoff radius
     real(double) :: delta           ! step
     real(double), pointer :: f(:)   ! Actual data
     real(double), pointer :: d2(:)  ! Second derivative
     ! Reciprocal space for integrals
     integer :: k_length
     real(double) :: k_delta
     real(double), pointer, dimension(:) :: k_table
  end type rad_func

  type pseudo_info
     character(len=80) :: filename   ! name of *.ion file
     integer :: lmax                 ! maximum of the angular momentum
     integer :: n_pjnl               ! number of projector functions
     logical :: flag_pcc             ! flag for partial core correction
     integer :: z                    ! atomic weight
     real(double) :: zval            ! number of valence electrons
     real(double) :: alpha           ! exponent of gaussian for long range term
     ! of local pseudopotential (might not be used)
     real(double) :: prefac          ! (alpha/pi)**1.5 stored for efficiency
     integer :: tm_loc_pot
     type(rad_func) :: vlocal            ! local potential
     type(rad_func) :: chlocal            ! local charge
     ! for Neutral atom potential
     type(rad_func) :: chna               ! neutral atom charge
     type(rad_func) :: vna                ! neutral atom potential
     real(double)   :: eshift             ! shift of vchlocal-zval/r
     
     type(rad_func) :: chpcc              ! partial core charge
     type(rad_func), pointer :: pjnl(:)   ! size = n_pjnl, projector functions
     integer, pointer :: pjnl_l(:)        ! size = n_pjnl, angular momentum
     integer, pointer :: pjnl_n(:)        ! size = n_pjnl, index for projectors
     real(double), pointer :: pjnl_ekb(:) ! size = n_pjnl, <phi_ln|dV_l|phi_ln>
     
     ! Neutral atom projectors
     integer :: n_pjna               ! number of NA projector functions
     type(rad_func), pointer :: pjna(:)   ! size = n_pjna, projector functions
     integer, pointer :: pjna_l(:)        ! size = n_pjna, angular momentum
     integer, pointer :: pjna_n(:)        ! size = n_pjna, index for projectors
     real(double), pointer :: pjna_ekb(:) ! size = n_pjna, <phi_ln|V_na|phi_ln>
     
  end type pseudo_info

  type(pseudo_info), allocatable :: pseudo(:)

  ! RCS tag for object file identification
  character(len=80), private :: RCSid = "$Id$"
!!***

contains

!---------------------------------------------------------------
! sbrt setup_pseudo_info
!---------------------------------------------------------------

!!****f* pseudo_tm_info/setup_pseudo_info *
!!
!!  NAME 
!!   setup_pseudo_info
!!  USAGE
!! 
!!  PURPOSE
!!   allocates pseudo(:), whose size is nspecies.
!!   sets up the name of files *.ion and reads the files.
!!  INPUTS
!!   nspecies, species_label
!!  USES
!! 
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   18/06/2002
!!  MODIFICATION HISTORY
!!   2008/05/25
!!    Added timers
!!   2011/03/30 M.Arita
!!    Added statements for P.C.C.
!!   2011/09/16 11:00 dave
!!    Changed to use n_species and species_label from species_module, and changed atomicrad to atomicnum
!!   2014/10/12 12:21 lat
!!    Added possibility of user defined filename for PAO/pseudo *ion
!!   2017/02/17 15:59 dave
!!    Changes to allow reading of Hamann's optimised norm-conserving pseudopotentials (remove ABINIT)
!!   2018/06/19 20:30 nakata
!!    Add IF with pseudo_type for calling init_rad when ghost atoms are used
!!  SOURCE
!!
  subroutine setup_pseudo_info

    use datatypes
    use numbers,        ONLY: zero
    use pao_format,     ONLY: pao
    use species_module, ONLY: n_species, species_label, species_file, species_from_files
    use species_module, ONLY: npao_species, nsf_species, type_species, charge
    use global_module,  ONLY: iprint_pseudo
    use dimens,         ONLY: RadiusSupport, atomicnum
    use GenComms,       ONLY: inode, ionode, cq_abort, gcopy
    use pseudopotential_common, ONLY: pseudo_type, SIESTA, ABINIT

    implicit none

    !local
    integer :: stat, ispecies, l, zeta
    real(double) :: cutoff
    character(len=80) :: filename
    integer :: ii

    if(allocated(pseudo)) then
       if(iprint_pseudo>2.AND.inode==ionode) write(io_lun,fmt='(10x," setup_pseudo_info is skipped because it is already called")')
    else
       call start_timer(tmr_std_allocation)
       allocate(pseudo(n_species),STAT=stat)
       if(stat /= 0) call cq_abort ('allocating pseudo in setup_pseudo_info',stat)
       allocate(pao(n_species),STAT=stat)
       if(stat /= 0) call cq_abort ('allocating pao in setup_pseudo_info',stat)
       call stop_timer(tmr_std_allocation)
       flag_pcc_global = .false.
       do ispecies=1,n_species
          if(pseudo_type==SIESTA) then
             pseudo(ispecies)%tm_loc_pot = loc_chg
          else if(pseudo_type==ABINIT) then
             pseudo(ispecies)%tm_loc_pot = loc_pot
          else
             call cq_abort("Error in pseudopotential type: ",pseudo_type)
          end if
          !**<lat>** 2014/10/12
          if ( species_from_files ) then
             write(filename,'(a,a)') trim(species_file(ispecies))
          else
             write(filename,'(a,a)') trim(species_label(ispecies)),".ion"
          end if
          if(iprint_pseudo>3.AND.inode==ionode) &
               write(io_lun,fmt='(10x,"ispecies = ",i5," file = ",a)') ispecies, filename

          pseudo(ispecies)%filename = filename
          call read_ion_ascii_tmp(pseudo(ispecies),pao(ispecies))

          npao_species(ispecies) = pao(ispecies)%count
          atomicnum(ispecies)    = pseudo(ispecies)%z
          if(charge(ispecies)/=pseudo(ispecies)%zval) then
             if(inode==ionode) &
                  write(io_lun,fmt='(/,2x,"WARNING: Mismatch between valence charge in input and pseudopotential: ",2f7.2,/)') &
                  charge(ispecies),pseudo(ispecies)%zval
          end if
          ! For P.C.C.
          if (pseudo(ispecies)%flag_pcc) flag_pcc_global = .true.
          !For Ghost atoms
          if(type_species(ispecies) < 0) then
            pseudo(ispecies)%zval = zero
           if(pseudo(ispecies)%n_pjnl > 0) then
            do ii=1, pseudo(ispecies)%n_pjnl
             call init_rad(pseudo(ispecies)%pjnl(ii))
            enddo
           endif
            if(pseudo_type==ABINIT) call init_rad(pseudo(ispecies)%vlocal)
            if(pseudo_type==SIESTA) call init_rad(pseudo(ispecies)%chlocal)
            if(pseudo(ispecies)%flag_pcc) call init_rad(pseudo(ispecies)%chpcc)
            pseudo(ispecies)%n_pjnl      = 0
            pseudo(ispecies)%flag_pcc    = .false.
            pseudo(ispecies)%pjnl_l(:)   = 0
            pseudo(ispecies)%pjnl_n(:)   = 0
            pseudo(ispecies)%pjnl_ekb(:) = zero
            ! pseudo(ispecies)%alpha
            ! pseudo(ispecies)%prefac
          endif

          if(npao_species(ispecies)<nsf_species(ispecies)) &
               call cq_abort("Error ! Less PAOs than SFs.  Decrease NumberOfSupports: ", &
               npao_species(ispecies),nsf_species(ispecies))
          cutoff = zero
          do l=0,pao(ispecies)%greatest_angmom
             if(pao(ispecies)%angmom(l)%n_zeta_in_angmom>0) then
                do zeta = 1,pao(ispecies)%angmom(l)%n_zeta_in_angmom
                   cutoff = max(cutoff,pao(ispecies)%angmom(l)%zeta(zeta)%cutoff)
                end do
             end if
          end do
          if(cutoff>RadiusSupport(ispecies)) then
             if(inode==ionode.AND.iprint_pseudo>0) &
                  write(io_lun,fmt='(10x,"Warning ! Species ",i3," support radius less than PAO radius ",f8.3)') ispecies,cutoff
             !RadiusSupport(ispecies) = cutoff
          endif
       enddo
    endif
    call gcopy(flag_pcc_global)
    if (iprint_pseudo>0.AND.inode==ionode .AND.flag_pcc_global) &
         write (io_lun,fmt='(10x,a)') "P.C.C. is taken into account."
    return
  end subroutine setup_pseudo_info
!!***

!---------------------------------------------------------------
! sbrt alloc_pseudo_info
!---------------------------------------------------------------

!!****f* pseudo_tm_info/alloc_pseudo_info *
!!
!!  NAME 
!!   alloc_pseudo_info
!!  USAGE
!! 
!!  PURPOSE
!!   Allocates memory for ps_info type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   May 2003
!!  MODIFICATION HISTORY
!!   2008/05/25
!!    Added timers
!!  SOURCE
!!
  subroutine alloc_pseudo_info(ps_info,n)

    use memory_module, ONLY: reg_alloc_mem, type_int, type_dbl    

    implicit none

    integer, intent(in) :: n
    type(pseudo_info),intent(out):: ps_info

    integer :: stat

    ps_info%n_pjnl = n
    call start_timer(tmr_std_allocation)
    allocate(ps_info%pjnl(n), STAT=stat)  
    if(stat /= 0) call cq_abort('alloc pjnl in alloc_pseudo_info',stat)
    allocate(ps_info%pjnl_l(n), ps_info%pjnl_n(n), STAT=stat)  
    if(stat /= 0) call cq_abort('alloc pjnl_l&n in alloc_pseudo_info',stat)
    call reg_alloc_mem(area_pseudo,n,type_int)
    allocate(ps_info%pjnl_ekb(n), STAT=stat)  
    if(stat /= 0) call cq_abort('alloc pjnl_ekb in alloc_pseudo_info',stat)
    call reg_alloc_mem(area_pseudo,n,type_dbl)
    call stop_timer(tmr_std_allocation)

    return
  end subroutine alloc_pseudo_info
!!***

!---------------------------------------------------------------
! sbrt rad_alloc : from TM
!---------------------------------------------------------------

!!****f* pseudo_tm_info/rad_alloc *
!!
!!  NAME 
!!   rad_alloc
!!  USAGE
!! 
!!  PURPOSE
!!   Sets the 'size' n of the arrays and allocates f and d2.
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!   2008/05/25
!!    Added timers
!!  SOURCE
!!
  subroutine rad_alloc(func,n)

    use memory_module, ONLY: reg_alloc_mem, type_dbl    

    implicit none

    type(rad_func), intent(out)    :: func
    integer, intent(in)        :: n
    integer :: stat

    func%n = n
    call start_timer(tmr_std_allocation)
    allocate(func%f(n),func%d2(n),STAT=stat)
    if(stat /= 0) call cq_abort('allocating radial func in rad_alloc',stat)
    call reg_alloc_mem(area_pseudo,2*n,type_dbl)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine rad_alloc
!!***

!---------------------------------------------------------------
! sbrt rad_dealloc 
!---------------------------------------------------------------

!!****f* pseudo_tm_info/rad_dealloc *
!!
!!  NAME 
!!   rad_dealloc
!!  USAGE
!! 
!!  PURPOSE
!!   Deallocates f and d2
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!   2008/05/25
!!    Added timers
!!  SOURCE
!!
  subroutine rad_dealloc(func)

    use memory_module, ONLY: reg_dealloc_mem, type_dbl    

    implicit none

    type(rad_func), intent(inout)    :: func

    integer :: stat

    call start_timer(tmr_std_allocation)
    deallocate(func%f,func%d2,STAT=stat)
    if(stat /= 0) call cq_abort('deallocating radial func in rad_alloc',stat)
    call reg_dealloc_mem(area_pseudo,2*func%n,type_dbl)
    nullify(func%f, func%d2)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine rad_dealloc
!!***

!---------------------------------------------------------------
! sbrt init_rad: from TM
!---------------------------------------------------------------

!!****f* pseudo_tm_info/init_rad *
!!
!!  NAME
!!   init_rad
!!  USAGE
!!
!!  PURPOSE
!!   initilise all members of rad_func
!!
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!
!!  MODIFICATION HISTORY
!!    2012/06/17 L.Tong
!!    - Changed func to intent(inout) from intent(out). If set to out,
!!      then all members of func upon entering the subroutine
!!      becomes undefined! This means f and d2 becomes disassociated
!!      and non-allocated.
!!  SOURCE
!!
  subroutine init_rad(func)

    use numbers, only: zero
    implicit none

    type(rad_func), intent(inout)    :: func

    func%n = 0
    func%cutoff = zero
    func%delta = zero
    func%f(:) = zero
    func%d2(:) = zero
    return
  end subroutine init_rad
!!***



!---------------------------------------------------------------
! sbrt radial_read_ascii : from SIESTA, with small changes
!---------------------------------------------------------------

!!****f* pseudo_tm_info/radial_read_ascii *
!!
!!  NAME 
!!   radial_read_ascii
!!  USAGE
!! 
!!  PURPOSE
!!   Reads radial table
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!   10:17, 2004/01/12 dave
!!    Changed formatting for reads to read new Siesta (1.3) pseudos
!!   2017/11/10 14:24 dave
!!    Changed yp1 and ypn to use first derivative (simple FD) instead of zero
!!   2018/03/08 09:54 dave
!!    Added consistency check between cutoff and step size/number of points
!!  SOURCE
!!
  subroutine radial_read_ascii(op,lun)

    use numbers, ONLY: BIG, zero, RD_ERR
    use spline_module, ONLY: spline
    use global_module, ONLY: iprint_pseudo
    use GenComms, ONLY: myid

    implicit none

    type(rad_func),intent(out) :: op
    integer,intent(in)         :: lun

    integer :: j, npts
    real(double) :: dummy, r0, r1
    real(double) :: yp1, ypn
    real(double) :: delta, cutoff

    !ori read(lun,'(i4,2g25.15)') npts, op%delta, op%cutoff
!    read(lun,'(i4,2g25.15)') npts, delta, cutoff
    read(lun,*) npts, delta, cutoff
    ! DRB 2018/03/08 Adding consistency check between cutoff and delta
    if(abs(cutoff - delta*real(npts-1,double))>RD_ERR) then
       if(myid==0) write(io_lun,fmt='(/4x,"Warning: cutoff and step inconsistent, cutoff will be adjusted: ",2f8.3/)') &
            cutoff,delta*real(npts-1,double)
       cutoff=delta*real(npts-1,double)
    end if
    if(myid==0.AND.iprint_pseudo>3) write(io_lun,fmt='(10x,"Radius: ",f15.10)') cutoff
    op%delta = delta
    op%cutoff= cutoff
    call rad_alloc(op,npts)
    ! Test step size in radial grid is consistent
    read(lun,*) r0, op%f(1)
    read(lun,*) r1, op%f(2)
    if(abs(r1-r0-delta)>RD_ERR) then
       if(myid==0) write(io_lun,fmt='(/4x,"Warning: radial grid and step inconsistent ! ",2f8.3/)') r1-r0,delta
       ! DRB 2018/03/08 10:55
       ! I can see an argument to abort here - if the step size and grid are inconsistent we may have problems
    end if
    do j=3,npts
!       read(lun,'(2g25.15)') dummy, op%f(j)
       read(lun,*) dummy, op%f(j)
    enddo
    !ori call rad_setup_d2(op)
    ! conquest version 
    !  yp1= BIG * 1.1
    !  ypn= BIG * 1.1
    yp1= (op%f(2)-op%f(1))/delta!zero
    ypn= (op%f(npts)-op%f(npts-1))/delta!zero
    call spline(op%n, op%delta, op%f, yp1, ypn, op%d2)
    return
  end subroutine radial_read_ascii
!!***

!---------------------------------------------------------------
! sbrt read_ion_ascii_tmp
!---------------------------------------------------------------

!!****f* pseudo_tm_info/read_ion_ascii_tmp *
!!
!!  NAME 
!!   read_ion_ascii_tmp
!!  USAGE
!! 
!!  PURPOSE
!!   this subroutine only reads KB projectors,
!!   chlocal and partial core charge
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!   13:04, 31/10/2003 drb 
!!    Added routine to create lmax for pseudos
!!   11:00, 13/02/2006 drb 
!!    Changed line length
!!   2006/10/17 08:18 dave
!!    Rewrote to read once on inode and broadcast
!!   2007/02/12 08:21 dave
!!    Added changes for TM potentials from ABINIT
!!   2007/08/16 15:36 dave
!!    More tweaks for local potential (particularly short/long range division)
!!   09:30, 27/11/2007 drb 
!!    Added calculation of prefac on all processors for local potential
!!   2008/05/25
!!    Added timers
!!   2008/09/01 08:21 dave
!!    Added io_ routines from input_module
!!   2013/07/05 dave
!!    Added reading and copying of z, charge state of ion
!!    Changed to read ALL radial functions for a given l (now picks up "semi-core" PAOs)
!!   2016/01/28 16:46 dave
!!    Updated module name to ion_electrostatic
!!   2016/02/09 08:17 dave
!!    Changed to use erfc from functions module
!!   2016/01/09 20:00 nakata
!!    Added pao_info%angmom%prncpl
!!   2017/02/24 10:06 dave & tsuyoshi
!!    Bug fix for read_header_tmp: removed formatting for zval
!!   2017/03/13 dave
!!    Repeated bug fix for other if branch in read_header_tmp (credit Jac van Driel for finding issue)
!!   2017/03/14 dave
!!    Changed string comparisons to use leqi from input_module (direct comparison breaks on Cray)
!!   2017/03/23 dave
!!    Added storage of spacing of PAO table
!!   2017/11/10 14:22 dave
!!    Bug fix: added Ry->Ha conversion for Siesta VNA d2 table (as well as f table)
!!   2018/01/22 14:39 JST dave
!!    Added test for lmax_pao/lmax_ps to find maximum PAO/PP angular momentum
!!   2018/10/30 11:02 dave
!!    Added semicore flag for each zeta
!!  SOURCE
!!
  subroutine read_ion_ascii_tmp(ps_info,pao_info)

    use numbers, ONLY: six, half, zero, pi, one
    use global_module, ONLY: numprocs, iprint_pseudo
    use GenComms, ONLY: inode, ionode, gcopy
    use pao_format, ONLY: species_pao
    use spline_module, ONLY : spline
    use memory_module, ONLY: reg_alloc_mem, type_dbl    
    use functions, ONLY: erfc
    use input_module, ONLY: io_assign, io_close
    use pseudopotential_common, ONLY: pseudo_type, SIESTA, ABINIT
    use maxima_module, only: lmax_pao, lmax_ps

    implicit none

    type(pseudo_info),intent(inout) :: ps_info
    type(species_pao),intent(inout) :: pao_info
    type(rad_func), dimension(:), allocatable :: dummy_rada
    type(rad_func) :: dummy_rad

    character(len=80) :: filename
    integer :: i, lun , i1, i2, i3, i4, z
    real(double) :: dummy, a, r
    integer :: n_orbnl, n_pjnl
    real(double) :: zval, yp1, ypn, erfarg, tmpv
    real(double), parameter :: ln10 = 2.302585092994_double

    integer :: iproc, lmax, maxz, alls, nzeta, l, count,tzl
    real(double), allocatable :: thispop(:)
    integer, allocatable :: thisl(:), thisn(:), thisz(:), zl(:), indexlz(:,:)

    if(inode==ionode) then
       filename=ps_info%filename
       call io_assign(lun)
       open(lun,file=filename,status='old',form='formatted')
       rewind(lun)

       call read_header_tmp(n_orbnl,lmax,n_pjnl, zval, z, lun)
       call alloc_pseudo_info(ps_info, n_pjnl)
       call start_timer(tmr_std_allocation)
       allocate(dummy_rada(n_orbnl),thisl(n_orbnl),thisn(n_orbnl),thisz(n_orbnl),thispop(n_orbnl),zl(0:lmax))
       call stop_timer(tmr_std_allocation)

       ps_info%z = z
       ps_info%zval = zval

       read(lun,*)
       zl = 0
       maxz = 0
       if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(10x,"Reading PAOs")')
       do i=1,n_orbnl
          read(lun,*) i1,i2,i3,i4, dummy
          thisl(i)=i1
          thisn(i)=i2
          !thisz(i)=i3
          thispop(i)=dummy
          zl(thisl(i)) = zl(thisl(i))+1
          thisz(i)=zl(thisl(i))
          !if(thisz(i)>zl(thisl(i))) zl(thisl(i))=thisz(i)
          maxz = max(maxz,thisz(i))
          if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(10x,"l: ",i3," z: ",i3)') i1,i3
          call radial_read_ascii(dummy_rada(i),lun)
       enddo
       allocate(indexlz(maxz,0:lmax))
       indexlz = 0
       do i=1,n_orbnl
          indexlz(thisz(i),thisl(i)) = i
          if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(10x,"indexlz: ",3i5)') thisz(i),thisl(i),i
       end do
       ! Now store data
       pao_info%greatest_angmom = lmax
       call start_timer(tmr_std_allocation)
       allocate(pao_info%angmom(0:lmax),STAT=alls)
       call stop_timer(tmr_std_allocation)
       count = 0
       if(alls/=0) call cq_abort('Failed to allocate PAOs')
       if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(10x,"Storing PAOs lmax: ",i5)') lmax
       do l=0,lmax
          pao_info%angmom(l)%n_zeta_in_angmom = zl(l)
          if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(10x,"l, zl: ",2i5)') l,zl(l)
          if(zl(l)>0) then
             call start_timer(tmr_std_allocation)
             allocate(pao_info%angmom(l)%zeta(zl(l)),pao_info%angmom(l)%prncpl(zl(l)), &
                  pao_info%angmom(l)%occ(zl(l)),pao_info%angmom(l)%semicore(zl(l)),STAT=alls)
             pao_info%angmom(l)%semicore(:) = 0
             if(alls/=0) call cq_abort('Failed to allocate PAOs zeta')
             call stop_timer(tmr_std_allocation)
             count = count + zl(l)*(2*l+1)
             do nzeta = 1,zl(l)
                i = indexlz(nzeta,l)
                if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(10x,"i,z,l: ",3i5)') i,nzeta,l
                pao_info%angmom(l)%zeta(nzeta)%length = dummy_rada(i)%n
                pao_info%angmom(l)%zeta(nzeta)%cutoff = dummy_rada(i)%cutoff
                pao_info%angmom(l)%zeta(nzeta)%delta = dummy_rada(i)%delta
                pao_info%angmom(l)%prncpl(nzeta) = thisn(i)
                pao_info%angmom(l)%occ(nzeta) = thispop(i)
                if(pao_info%angmom(l)%zeta(nzeta)%length>=1) then
                   call start_timer(tmr_std_allocation)
                   allocate(pao_info%angmom(l)%zeta(nzeta)%table(pao_info%angmom(l)%zeta(nzeta)%length),STAT=alls)
                   if(alls/=0) call cq_abort('Failed to allocate PAOs zeta')
                   call reg_alloc_mem(area_pseudo,pao_info%angmom(l)%zeta(nzeta)%length,type_dbl)
                   call stop_timer(tmr_std_allocation)
                   pao_info%angmom(l)%zeta(nzeta)%table(1:dummy_rada(i)%n) = dummy_rada(i)%f(1:dummy_rada(i)%n)
                   call start_timer(tmr_std_allocation)
                   allocate(pao_info%angmom(l)%zeta(nzeta)%table2(pao_info%angmom(l)%zeta(nzeta)%length),STAT=alls)
                   if(alls/=0) call cq_abort('Failed to allocate PAOs zeta')
                   call reg_alloc_mem(area_pseudo,pao_info%angmom(l)%zeta(nzeta)%length,type_dbl)
                   call stop_timer(tmr_std_allocation)
                else
                   call cq_abort('PAO with zero length: ',l,nzeta)
                end if
                yp1 = (pao_info%angmom(l)%zeta(nzeta)%table(2)-pao_info%angmom(l)%zeta(nzeta)%table(1))/ &
                     dummy_rada(i)%delta
                ypn = (pao_info%angmom(l)%zeta(nzeta)%table(dummy_rada(i)%n)-&
                     pao_info%angmom(l)%zeta(nzeta)%table(dummy_rada(i)%n-1))/ &
                     dummy_rada(i)%delta
                call spline(dummy_rada(i)%n,dummy_rada(i)%delta,pao_info%angmom(l)%zeta(nzeta)%table, &
                     yp1,ypn,pao_info%angmom(l)%zeta(nzeta)%table2)
             end do
          end if
       end do
       pao_info%count = count
       deallocate(indexlz)
       do i=1,n_orbnl
          call rad_dealloc(dummy_rada(i))
       end do
       deallocate(dummy_rada,thisl,thisn,thisz,thispop,zl)

       ! KBs
       if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(10x,"Reading KB projectors ")')
       read(lun,*)
       ps_info%lmax = 0
       do i=1,ps_info%n_pjnl
!          read(lun,'(2i3,f22.16)') &
          read(lun,*) &
               ps_info%pjnl_l(i), ps_info%pjnl_n(i), ps_info%pjnl_ekb(i)
          !!   25/Jul/2002 TM : Rydberg units -> Hartree units
          ps_info%pjnl_ekb(i) = half* ps_info%pjnl_ekb(i)
          if(iprint_pseudo>1.AND.inode==ionode) write(io_lun,fmt='(10x,"PSEUDO_TM: i, l,n, ekb = ",2i3,f22.16)') &
               ps_info%pjnl_l(i), ps_info%pjnl_n(i),ps_info%pjnl_ekb(i)
          if(ps_info%pjnl_l(i)>ps_info%lmax) ps_info%lmax = ps_info%pjnl_l(i)
          call radial_read_ascii(ps_info%pjnl(i),lun)
       enddo
       !Vlocal
       if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(10x,"Reading Vna ")')
       read(lun,*)
       ! Read and store local potential; it may be local part of pseudo, or neutral atom
       ! DRB 2017/02/21 switch here: after vna, get vlocal for Hamann, chlocal for siesta
       ! DRB 2017/02/21 This is the NA potential
       if(pseudo_type==SIESTA) then
          call radial_read_ascii(ps_info%vna,lun)
          ps_info%vna%f  = half*ps_info%vna%f  ! Siesta works in Ry
          ps_info%vna%d2 = half*ps_info%vna%d2 ! Siesta works in Ry
          if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(10x,"Reading Chlocal ")')
          read(lun,*)
          call radial_read_ascii(ps_info%chlocal,lun)
       else if(pseudo_type==ABINIT) then
          call radial_read_ascii(ps_info%vna,lun)
          if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(10x,"Reading Vlocal ")')
          read(lun,*)
          call radial_read_ascii(ps_info%vlocal,lun)
       end if
       !sets up gaussian which is used to calculate G=0 term
       !  exp( - alpha * cutoff **2 ) = 10 **(-six) 
       if(ps_info%tm_loc_pot==loc_pot) then
          ps_info%alpha = six * ln10 / (ps_info%vlocal%cutoff ** 2)
          ps_info%prefac = (ps_info%alpha/pi)**1.5_double
          ! Remove long-range part from short-range part
          a = sqrt(ps_info%alpha)
          !r = ps_info%vlocal%delta
          r = 1.0e-7_double
          erfarg = one - erfc(a*r)
          !erfarg = derf(a*r)
          tmpv = ps_info%vlocal%f(1)
          ps_info%vlocal%f(1) = ps_info%vlocal%f(1) + ps_info%zval * erfarg/r
          !write(52,*) r,ps_info%vlocal%f(1),tmpv,ps_info%zval * erfarg/r!ps_info%zval * derf(a*r)/r
          do i=2,ps_info%vlocal%n
             r = ps_info%vlocal%delta*real(i-1,double)
             erfarg = one - erfc(a*r)
             tmpv = ps_info%vlocal%f(i)
             ps_info%vlocal%f(i) = ps_info%vlocal%f(i) + ps_info%zval * erfarg/r
             !write(52,*) r,ps_info%vlocal%f(i), tmpv,ps_info%zval * erfarg/r!,ps_info%zval * derf(a*r)/r
          end do
          call spline(ps_info%vlocal%n, ps_info%vlocal%delta, ps_info%vlocal%f, zero, zero, ps_info%vlocal%d2)
       else
          ps_info%alpha = six * ln10 / (ps_info%chlocal%cutoff ** 2)
       endif
       !Core
       ! 2017/03/17 dave
       ! I'm removing this because we read the PCC flag from the pseudopotential header
       !ps_info%flag_pcc = .false.
       !! I *really* don't like this, but we'll stay with it for now
       !read(lun,*,end=9999)
       if(ps_info%flag_pcc) then
          read(lun,*)
          if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(10x,"Reading pcc ")')
          call radial_read_ascii(ps_info%chpcc,lun)
       end if
       !ps_info%flag_pcc = .true.

!9999   continue
       call io_close(lun)
    endif !  (inode == ionode) then
    ! Now broadcast the information
    call gcopy(lmax   )
    if(lmax>lmax_pao) lmax_pao = lmax    
    call gcopy(n_pjnl )
    call gcopy(zval)
    call gcopy(z)
    if(inode/=ionode) then
       call alloc_pseudo_info(ps_info, n_pjnl)
       ps_info%zval = zval
       ps_info%z = z
       pao_info%greatest_angmom = lmax
       allocate(pao_info%angmom(0:lmax),STAT=alls)
    end if
    if(numprocs>1) then
       count = 0
       if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,*) 'PAOs'
       do l=0,lmax
          if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(10x," l is ",i5)') l
          call gcopy(pao_info%angmom(l)%n_zeta_in_angmom)
          tzl = pao_info%angmom(l)%n_zeta_in_angmom
          if(tzl>0) then
             if(inode/=ionode) then
                call start_timer(tmr_std_allocation)
                allocate(pao_info%angmom(l)%zeta(tzl),pao_info%angmom(l)%prncpl(tzl),pao_info%angmom(l)%occ(tzl),STAT=alls)
                if(alls/=0) call cq_abort('Failed to allocate PAOs zeta')
                call stop_timer(tmr_std_allocation)
             end if
             count = count + tzl*(2*l+1)
             do nzeta = 1,tzl
                call gcopy(pao_info%angmom(l)%zeta(nzeta)%length)
                call gcopy(pao_info%angmom(l)%zeta(nzeta)%cutoff)
                call gcopy(pao_info%angmom(l)%zeta(nzeta)%delta)
                call gcopy(pao_info%angmom(l)%prncpl(nzeta))
                call gcopy(pao_info%angmom(l)%occ(nzeta))
                if(inode/=ionode) then
                   if(pao_info%angmom(l)%zeta(nzeta)%length>=1) then
                      call start_timer(tmr_std_allocation)
                      allocate(pao_info%angmom(l)%zeta(nzeta)%table(pao_info%angmom(l)%zeta(nzeta)%length),STAT=alls)
                      if(alls/=0) call cq_abort('Failed to allocate PAOs zeta',pao_info%angmom(l)%zeta(nzeta)%length)
                      call reg_alloc_mem(area_pseudo,pao_info%angmom(l)%zeta(nzeta)%length,type_dbl)

                      allocate(pao_info%angmom(l)%zeta(nzeta)%table2(pao_info%angmom(l)%zeta(nzeta)%length),STAT=alls)
                      if(alls/=0) call cq_abort('Failed to allocate PAOs zeta ',pao_info%angmom(l)%zeta(nzeta)%length)
                      call reg_alloc_mem(area_pseudo,pao_info%angmom(l)%zeta(nzeta)%length,type_dbl)
                      call stop_timer(tmr_std_allocation)
                   end if
                end if
                call gcopy(pao_info%angmom(l)%zeta(nzeta)%table,pao_info%angmom(l)%zeta(nzeta)%length)
                call gcopy(pao_info%angmom(l)%zeta(nzeta)%table2,pao_info%angmom(l)%zeta(nzeta)%length)
             end do
          end if
       end do
       pao_info%count = count
       if(inode/=ionode) ps_info%lmax = 0
       if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,*) 'PPs'
       do i=1,ps_info%n_pjnl
          if(iprint_pseudo>3.AND.inode==ionode) write(io_lun,fmt='(2x,"Projector number: ",i3)') i
          call gcopy(ps_info%pjnl_l(i))
          call gcopy(ps_info%pjnl_n(i))
          call gcopy(ps_info%pjnl_ekb(i))
          if(inode/=ionode.AND.ps_info%pjnl_l(i)>ps_info%lmax) ps_info%lmax = ps_info%pjnl_l(i)
          call gcopy(ps_info%pjnl(i)%delta)
          call gcopy(ps_info%pjnl(i)%cutoff)
          call gcopy(ps_info%pjnl(i)%n)
          if(inode/=ionode) call rad_alloc(ps_info%pjnl(i),ps_info%pjnl(i)%n)
          call gcopy(ps_info%pjnl(i)%f,ps_info%pjnl(i)%n)
          call gcopy(ps_info%pjnl(i)%d2,ps_info%pjnl(i)%n)
       end do
       if(ps_info%lmax>lmax_ps) lmax_ps = ps_info%lmax
       if(pseudo_type==SIESTA) then
          !Chlocal
          call gcopy(ps_info%tm_loc_pot)
          call gcopy(ps_info%chlocal%n)
          call gcopy(ps_info%chlocal%cutoff)
          call gcopy(ps_info%chlocal%delta)
          if(inode/=ionode) call rad_alloc(ps_info%chlocal,ps_info%chlocal%n)
          call gcopy(ps_info%chlocal%f,ps_info%chlocal%n)
          call gcopy(ps_info%chlocal%d2,ps_info%chlocal%n)
          ! Vna
          call gcopy(ps_info%vna%n)
          call gcopy(ps_info%vna%cutoff)
          call gcopy(ps_info%vna%delta)
          if(inode/=ionode) call rad_alloc(ps_info%vna,ps_info%vna%n)
          call gcopy(ps_info%vna%f,ps_info%vna%n)
          call gcopy(ps_info%vna%d2,ps_info%vna%n)
       else if(pseudo_type==ABINIT) then
          !Vlocal
          call gcopy(ps_info%tm_loc_pot)
          call gcopy(ps_info%vlocal%n)
          call gcopy(ps_info%vlocal%cutoff)
          call gcopy(ps_info%vlocal%delta)
          if(inode/=ionode) call rad_alloc(ps_info%vlocal,ps_info%vlocal%n)
          call gcopy(ps_info%vlocal%f,ps_info%vlocal%n)
          call gcopy(ps_info%vlocal%d2,ps_info%vlocal%n)
          ! Vna
          call gcopy(ps_info%vna%n)
          call gcopy(ps_info%vna%cutoff)
          call gcopy(ps_info%vna%delta)
          if(inode/=ionode) call rad_alloc(ps_info%vna,ps_info%vna%n)
          call gcopy(ps_info%vna%f,ps_info%vna%n)
          call gcopy(ps_info%vna%d2,ps_info%vna%n)
       end if
       if(inode/=ionode) then
          if(ps_info%tm_loc_pot==loc_pot) then
             ps_info%alpha = six * ln10 / (ps_info%vlocal%cutoff ** 2)
             ps_info%prefac = (ps_info%alpha/pi)**1.5_double
          else
             ps_info%alpha = six * ln10 / (ps_info%chlocal%cutoff ** 2)
          end if
       end if
       call gcopy(ps_info%flag_pcc)
       if(ps_info%flag_pcc) then
          call gcopy(ps_info%chpcc%n)
          call gcopy(ps_info%chpcc%cutoff)
          call gcopy(ps_info%chpcc%delta)
          if(inode/=ionode) call rad_alloc(ps_info%chpcc,ps_info%chpcc%n)
          call gcopy(ps_info%chpcc%f,ps_info%chpcc%n)
          call gcopy(ps_info%chpcc%d2,ps_info%chpcc%n)
       end if
    end if ! numprocs>1
  contains

    subroutine read_header_tmp(n_orbnl, lmax_basis, n_pjnl, zval, z, unit)

      use global_module, ONLY: iprint_pseudo
      use input_module, ONLY: leqi

      implicit none

      integer, intent(in)         :: unit
      integer, intent(out) :: n_orbnl, n_pjnl, z
      real(double), intent(out) :: zval

      character(len=78) :: line, trim_line

      character(len=2)  :: symbol
      character(len=20) :: label
      integer :: lmax_basis, lmax_projs, zval_int
      real(double) :: mass, self_energy

      ! For judging P.C.C.
      character(len=2) :: symbol2, xc
      character(len=3) :: rel
      character(len=4) :: pcc

      ! Judge if P.C.C. is considered
      read(unit,'(a)') line
      trim_line = trim(line)
      if (leqi(trim_line(1:10),'<preamble>')) then
         read(unit,'(a)') line
         trim_line = trim(line)
         do while(.NOT.leqi(trim_line(1:11),'</preamble>'))
            read(unit,'(a)') line
            trim_line = trim(line)
            if (leqi(trim_line(1:24),'<pseudopotential_header>')) then
              read (unit, '(1x,a2,1x,a2,1x,a3,1x,a4)') symbol2, xc, rel, pcc
              if (leqi(pcc,'pcec')) then
                ps_info%flag_pcc = .true.
              else !if (pcc .NE. 'pcec') then
                ps_info%flag_pcc = .false.
              endif
            endif
         end do
      endif

      if(iprint_pseudo>1) then
         read(unit,'(a2)') symbol
         write(io_lun,fmt='(10x,"symbol = ",a2)') symbol
         read(unit,'(a20)') label
         write(io_lun,fmt='(10x,"label = ",a20)') label
         read(unit,*) z
         write(io_lun,fmt='(10x,"z = ",i5)') z
         read(unit,*) zval
         write(io_lun,fmt='(10x,"zval_int, zval = ",f12.5)') zval
         read(unit,*) mass
         write(io_lun,fmt='(10x,"mass = ",g25.15)') mass
         read(unit,*) self_energy
         write(io_lun,fmt='(10x,"self_energy = ",g25.15)') self_energy
         read(unit,'(2i4)') lmax_basis, n_orbnl
         write(io_lun,fmt='(10x,"lmax_basis, n_orbnl = ",2i4)') lmax_basis, n_orbnl
         read(unit,'(2i4)') lmax_projs, n_pjnl
         write(io_lun,fmt='(10x,"lmax_projs, n_pjnl = ",2i4)') lmax_projs, n_pjnl
      else
         read(unit,'(a2)') symbol
         read(unit,'(a20)') label
         read(unit,*) z
         read(unit,*) zval!_int
         !zval=zval_int
         read(unit,'(g25.15)') mass
         read(unit,'(g25.15)') self_energy
         read(unit,'(2i4)') lmax_basis, n_orbnl
         read(unit,'(2i4)') lmax_projs, n_pjnl
      end if
    end subroutine read_header_tmp
  end subroutine read_ion_ascii_tmp

!!***
end module pseudo_tm_info
