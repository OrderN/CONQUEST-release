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
!!
!!  SOURCE
!!
module pseudo_tm_info

  use datatypes, ONLY: double
  use GenComms, ONLY: cq_abort
  use global_module, ONLY: area_pseudo

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
  end type rad_func

  type pseudo_info
     character(len=20) :: filename   ! name of *.ion file
     integer :: lmax                 ! maximum of the angular momentum
     integer :: n_pjnl               ! number of projector functions
     logical :: flag_pcc             ! flag for partial core correction
     real(double) :: zval            ! number of valence electrons
     real(double) :: alpha           ! exponent of gaussian for long range term
     ! of local pseudopotential (might not be used)
     real(double) :: prefac          ! (alpha/pi)**1.5 stored for efficiency
     integer :: tm_loc_pot
     type(rad_func) :: vlocal            ! local potential
     type(rad_func) :: chlocal            ! local charge
     type(rad_func) :: chpcc              ! partial core charge
     type(rad_func), pointer :: pjnl(:)   ! size = n_pjnl, projector functions
     integer, pointer :: pjnl_l(:)        ! size = n_pjnl, angular momentum
     integer, pointer :: pjnl_n(:)        ! size = n_pjnl, index for projectors
     real(double), pointer :: pjnl_ekb(:) ! size = n_pjnl, <phi_ln|dV_l|phi_ln>
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
!!
!!  SOURCE
!!
  subroutine setup_pseudo_info(nspecies, species_label)

    use datatypes
    use numbers, ONLY: zero
    use pao_format, ONLY: pao
    use species_module, ONLY: npao_species, nsf_species, type_species
    use global_module, ONLY: iprint_pseudo
    use dimens, ONLY: RadiusSupport
    use GenComms, ONLY: inode, ionode, cq_abort
    use pseudopotential_common, ONLY: pseudo_type, SIESTA, ABINIT

    implicit none

    !dummy arguments
    integer, intent(in) :: nspecies
    character(len=10),intent(in) :: species_label(nspecies)
    !local
    integer :: stat, ispecies, l, zeta
    real(double) :: cutoff
    character(len=20) :: filename
    integer :: ii

    if(allocated(pseudo)) then
       if(iprint_pseudo>2) write(*,fmt='(10x," setup_pseudo_info is skipped because it is already called")')
    else
       allocate(pseudo(nspecies),STAT=stat)
       if(stat /= 0) call cq_abort ('allocating pseudo in setup_pseudo_info',stat)
       allocate(pao(nspecies),STAT=stat)
       if(stat /= 0) call cq_abort ('allocating pao in setup_pseudo_info',stat)
       do ispecies=1,nspecies
          if(pseudo_type==SIESTA) then
             pseudo(ispecies)%tm_loc_pot = loc_chg
          else if(pseudo_type==ABINIT) then
             pseudo(ispecies)%tm_loc_pot = loc_pot
          else
             call cq_abort("Error in pseudopotential type: ",pseudo_type)
          end if
          write(filename,'(a,a)') trim(species_label(ispecies)), ".ion"
          if(iprint_pseudo>3) write(*,fmt='(10x,"ispecies = ",i5," file = ",a)') ispecies,filename
          pseudo(ispecies)%filename = filename
          call read_ion_ascii_tmp(pseudo(ispecies),pao(ispecies))
          npao_species(ispecies) = pao(ispecies)%count
          !For Ghost atoms
          if(type_species(ispecies) < 0) then
            pseudo(ispecies)%zval = zero
           if(pseudo(ispecies)%n_pjnl > 0) then
            do ii=1, pseudo(ispecies)%n_pjnl
             call init_rad(pseudo(ispecies)%pjnl(ii))
            enddo
           endif
            call init_rad(pseudo(ispecies)%vlocal)
            call init_rad(pseudo(ispecies)%chlocal)
            if(pseudo(ispecies)%flag_pcc) call init_rad(pseudo(ispecies)%chpcc)
            pseudo(ispecies)%n_pjnl = 0
            pseudo(ispecies)%flag_pcc = .false.
            pseudo(ispecies)%pjnl_l(:) = 0
            pseudo(ispecies)%pjnl_n(:) = 0
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
                  write(*,fmt='(10x,"Warning ! Species ",i3," support radius less than PAO radius ",f8.3)') ispecies,cutoff
             !RadiusSupport(ispecies) = cutoff
          endif
       enddo
    endif
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
!!
!!  SOURCE
!!
  subroutine alloc_pseudo_info(ps_info,n)

    use memory_module, ONLY: reg_alloc_mem, type_int, type_dbl    

    implicit none

    integer, intent(in) :: n
    type(pseudo_info),intent(out):: ps_info

    integer :: stat

    ps_info%n_pjnl = n
    allocate(ps_info%pjnl(n), STAT=stat)  
    if(stat /= 0) call cq_abort('alloc pjnl in alloc_pseudo_info',stat)
    allocate(ps_info%pjnl_l(n), ps_info%pjnl_n(n), STAT=stat)  
    if(stat /= 0) call cq_abort('alloc pjnl_l&n in alloc_pseudo_info',stat)
    call reg_alloc_mem(area_pseudo,n,type_int)
    allocate(ps_info%pjnl_ekb(n), STAT=stat)  
    if(stat /= 0) call cq_abort('alloc pjnl_ekb in alloc_pseudo_info',stat)
    call reg_alloc_mem(area_pseudo,n,type_dbl)

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
!!
!!  SOURCE
!!
  subroutine rad_alloc(func,n)

    use memory_module, ONLY: reg_alloc_mem, type_dbl    

    implicit none

    type(rad_func), intent(out)    :: func
    integer, intent(in)        :: n
    integer :: stat

    func%n = n
    allocate(func%f(n),func%d2(n),STAT=stat)
    if(stat /= 0) call cq_abort('allocating radial func in rad_alloc',stat)
    call reg_alloc_mem(area_pseudo,2*n,type_dbl)
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
!!
!!  SOURCE
!!
  subroutine rad_dealloc(func)

    use memory_module, ONLY: reg_dealloc_mem, type_dbl    

    implicit none

    type(rad_func), intent(inout)    :: func

    integer :: stat

    deallocate(func%f,func%d2,STAT=stat)
    if(stat /= 0) call cq_abort('deallocating radial func in rad_alloc',stat)
    call reg_dealloc_mem(area_pseudo,2*func%n,type_dbl)
    nullify(func%f, func%d2)
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
!!
!!  SOURCE
!!
  subroutine init_rad(func)

    use numbers, only: zero
    implicit none

    type(rad_func), intent(out)    :: func

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
!!  SOURCE
!!
  subroutine radial_read_ascii(op,lun)

    use numbers, ONLY: BIG, zero
    use spline_module, ONLY: spline
    use global_module, ONLY: iprint_pseudo
    use GenComms, ONLY: myid

    implicit none

    type(rad_func),intent(out) :: op
    integer,intent(in)         :: lun

    integer :: j, npts
    real(double) :: dummy
    real(double) :: yp1, ypn
    real(double) :: delta, cutoff

    !ori read(lun,'(i4,2g25.15)') npts, op%delta, op%cutoff
!    read(lun,'(i4,2g25.15)') npts, delta, cutoff
    read(lun,*) npts, delta, cutoff
    if(myid==0.AND.iprint_pseudo>3) write(*,fmt='(10x,"Radius: ",f15.10)') cutoff
    op%delta = delta
    op%cutoff= cutoff
    call rad_alloc(op,npts)
    do j=1,npts
!       read(lun,'(2g25.15)') dummy, op%f(j)
       read(lun,*) dummy, op%f(j)
    enddo
    !ori call rad_setup_d2(op)
    ! conquest version 
    !  yp1= BIG * 1.1
    !  ypn= BIG * 1.1
    yp1= zero
    ypn= zero
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
!!  SOURCE
!!
  subroutine read_ion_ascii_tmp(ps_info,pao_info)

    use numbers, ONLY: six, half, zero, pi, one
    use global_module, ONLY: numprocs, iprint_pseudo
    use GenComms, ONLY: inode, ionode, gcopy
    use pao_format, ONLY: species_pao
    use spline_module, ONLY : spline
    use memory_module, ONLY: reg_alloc_mem, type_dbl    
    use ewald_module, ONLY: erfc

    implicit none

    type(pseudo_info),intent(inout) :: ps_info
    type(species_pao),intent(inout) :: pao_info
    type(rad_func), dimension(:), allocatable :: dummy_rada
    type(rad_func) :: dummy_rad

    character(len=20) :: filename
    integer :: i, lun , i1, i2, i3, i4
    real(double) :: dummy, a, r
    integer :: n_orbnl, n_pjnl
    real(double) :: zval, yp1, ypn, erfarg, tmpv
    real(double), parameter :: ln10 = 2.302585092994_double

    integer :: iproc, lmax, maxz, alls, nzeta, l, count,tzl
    real(double), allocatable :: thispop(:)
    integer, allocatable :: thisl(:), thisz(:), zl(:), indexlz(:,:)

    if(inode==ionode) then
       filename=ps_info%filename
       call io_assign(lun)
       open(lun,file=filename,status='old',form='formatted')
       rewind(lun)

       call read_header_tmp(n_orbnl,lmax,n_pjnl, zval, lun)
       call alloc_pseudo_info(ps_info, n_pjnl)
       allocate(dummy_rada(n_orbnl),thisl(n_orbnl),thisz(n_orbnl),thispop(n_orbnl),zl(0:lmax))

       ps_info%zval = zval

       read(lun,*)
       zl = 0
       maxz = 0
       if(iprint_pseudo>3) write(*,fmt='(10x,"Reading PAOs")')
       do i=1,n_orbnl
          read(lun,*) i1,i2,i3,i4, dummy
          thisl(i)=i1
          thisz(i)=i3
          thispop(i)=dummy
          if(thisz(i)>zl(thisl(i))) zl(thisl(i))=thisz(i)
          maxz = max(maxz,thisz(i))
          if(iprint_pseudo>3) write(*,fmt='(10x,"l: ",i3," z: ",i3)') i1,i3
          call radial_read_ascii(dummy_rada(i),lun)
       enddo
       allocate(indexlz(maxz,0:lmax))
       indexlz = 0
       do i=1,n_orbnl
          indexlz(thisz(i),thisl(i)) = i
          if(iprint_pseudo>3) write(*,fmt='(10x,"indexlz: ",3i5)') thisz(i),thisl(i),i
       end do
       ! Now store data
       pao_info%greatest_angmom = lmax
       allocate(pao_info%angmom(0:lmax),STAT=alls)
       count = 0
       if(alls/=0) call cq_abort('Failed to allocate PAOs')
       if(iprint_pseudo>3) write(*,fmt='(10x,"Storing PAOs lmax: ",i5)') lmax
       do l=0,lmax
          pao_info%angmom(l)%n_zeta_in_angmom = zl(l)
          if(iprint_pseudo>3) write(*,fmt='(10x,"l, zl: ",2i5)') l,zl(l)
          if(zl(l)>0) then
             allocate(pao_info%angmom(l)%zeta(zl(l)),pao_info%angmom(l)%occ(zl(l)),STAT=alls)
             if(alls/=0) call cq_abort('Failed to allocate PAOs zeta')
             count = count + zl(l)*(2*l+1)
             do nzeta = 1,zl(l)
                i = indexlz(nzeta,l)
                if(iprint_pseudo>3) write(*,fmt='(10x,"i,z,l: ",3i5)') i,nzeta,l
                pao_info%angmom(l)%zeta(nzeta)%length = dummy_rada(i)%n
                pao_info%angmom(l)%zeta(nzeta)%cutoff = dummy_rada(i)%cutoff
                pao_info%angmom(l)%occ(nzeta) = thispop(i)
                if(pao_info%angmom(l)%zeta(nzeta)%length>=1) then
                   allocate(pao_info%angmom(l)%zeta(nzeta)%table(pao_info%angmom(l)%zeta(nzeta)%length),STAT=alls)
                   if(alls/=0) call cq_abort('Failed to allocate PAOs zeta')
                   call reg_alloc_mem(area_pseudo,pao_info%angmom(l)%zeta(nzeta)%length,type_dbl)
                   pao_info%angmom(l)%zeta(nzeta)%table(1:dummy_rada(i)%n) = dummy_rada(i)%f(1:dummy_rada(i)%n)
                   allocate(pao_info%angmom(l)%zeta(nzeta)%table2(pao_info%angmom(l)%zeta(nzeta)%length),STAT=alls)
                   if(alls/=0) call cq_abort('Failed to allocate PAOs zeta')
                   call reg_alloc_mem(area_pseudo,pao_info%angmom(l)%zeta(nzeta)%length,type_dbl)
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
       deallocate(dummy_rada,thisl,thisz,thispop,zl)

       ! KBs
       if(iprint_pseudo>3) write(*,fmt='(10x,"Reading KB projectors ")')
       read(lun,*)
       ps_info%lmax = 0
       do i=1,ps_info%n_pjnl
          read(lun,'(2i3,f22.16)') &
               ps_info%pjnl_l(i), ps_info%pjnl_n(i), ps_info%pjnl_ekb(i)
          !!   25/Jul/2002 TM : Rydberg units -> Hartree units
          ps_info%pjnl_ekb(i) = half* ps_info%pjnl_ekb(i)
          if(iprint_pseudo>1) write(*,fmt='(10x,"PSEUDO_TM: i, l,n, ekb = ",2i3,f22.16)') &
               ps_info%pjnl_l(i), ps_info%pjnl_n(i),ps_info%pjnl_ekb(i)
          if(ps_info%pjnl_l(i)>ps_info%lmax) ps_info%lmax = ps_info%pjnl_l(i)
          call radial_read_ascii(ps_info%pjnl(i),lun)
       enddo
       !Vlocal
       if(iprint_pseudo>3) write(*,fmt='(10x,"Reading Vlocal ")')
       read(lun,*)
       ! Read and store local potential; it may be local part of pseudo, or neutral atom
       call radial_read_ascii(ps_info%vlocal,lun)
       !call radial_read_ascii(dummy_rad,lun)
       !call rad_dealloc(dummy_rad)
       !Chlocal
       if(iprint_pseudo>3) write(*,fmt='(10x,"Reading Chlocal ")')
       read(lun,*)
       call radial_read_ascii(ps_info%chlocal,lun)
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
       ps_info%flag_pcc = .false.
       ! I *really* don't like this, but we'll stay with it for now
       read(lun,*,end=9999)
       if(iprint_pseudo>3) write(*,fmt='(10x,"Reading pcc ")')
       call radial_read_ascii(ps_info%chpcc,lun)
       ps_info%flag_pcc = .true.

9999   continue
       call io_close(lun)
    endif !  (inode == ionode) then
    ! Now broadcast the information
    call gcopy(lmax   )
    call gcopy(n_pjnl )
    call gcopy(zval)
    if(inode/=ionode) then
       call alloc_pseudo_info(ps_info, n_pjnl)
       ps_info%zval = zval
       pao_info%greatest_angmom = lmax
       allocate(pao_info%angmom(0:lmax),STAT=alls)
    end if
    if(numprocs>1) then
       count = 0
       do l=0,lmax
          if(iprint_pseudo>3) write(*,fmt='(10x," l is ",i5)') l
          call gcopy(pao_info%angmom(l)%n_zeta_in_angmom)
          tzl = pao_info%angmom(l)%n_zeta_in_angmom
          if(tzl>0) then
             if(inode/=ionode) then
                allocate(pao_info%angmom(l)%zeta(tzl),pao_info%angmom(l)%occ(tzl),STAT=alls)
                if(alls/=0) call cq_abort('Failed to allocate PAOs zeta')
             end if
             count = count + tzl*(2*l+1)
             do nzeta = 1,tzl
                call gcopy(pao_info%angmom(l)%zeta(nzeta)%length)
                call gcopy(pao_info%angmom(l)%zeta(nzeta)%cutoff)
                call gcopy(pao_info%angmom(l)%occ(nzeta))
                if(inode/=ionode) then
                   if(pao_info%angmom(l)%zeta(nzeta)%length>=1) then
                      allocate(pao_info%angmom(l)%zeta(nzeta)%table(pao_info%angmom(l)%zeta(nzeta)%length),STAT=alls)
                      if(alls/=0) call cq_abort('Failed to allocate PAOs zeta',pao_info%angmom(l)%zeta(nzeta)%length)
                      call reg_alloc_mem(area_pseudo,pao_info%angmom(l)%zeta(nzeta)%length,type_dbl)

                      allocate(pao_info%angmom(l)%zeta(nzeta)%table2(pao_info%angmom(l)%zeta(nzeta)%length),STAT=alls)
                      if(alls/=0) call cq_abort('Failed to allocate PAOs zeta ',pao_info%angmom(l)%zeta(nzeta)%length)
                      call reg_alloc_mem(area_pseudo,pao_info%angmom(l)%zeta(nzeta)%length,type_dbl)
                   end if
                end if
                call gcopy(pao_info%angmom(l)%zeta(nzeta)%table,pao_info%angmom(l)%zeta(nzeta)%length)
                call gcopy(pao_info%angmom(l)%zeta(nzeta)%table2,pao_info%angmom(l)%zeta(nzeta)%length)
             end do
          end if
       end do
       pao_info%count = count
       if(inode/=ionode) ps_info%lmax = 0
       do i=1,ps_info%n_pjnl
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
       !Vlocal
       call gcopy(ps_info%tm_loc_pot)
       call gcopy(ps_info%vlocal%n)
       call gcopy(ps_info%vlocal%cutoff)
       call gcopy(ps_info%vlocal%delta)
       if(inode/=ionode) call rad_alloc(ps_info%vlocal,ps_info%vlocal%n)
       call gcopy(ps_info%vlocal%f,ps_info%vlocal%n)
       call gcopy(ps_info%vlocal%d2,ps_info%vlocal%n)
       !Chlocal
       call gcopy(ps_info%chlocal%n)
       call gcopy(ps_info%chlocal%cutoff)
       call gcopy(ps_info%chlocal%delta)
       if(inode/=ionode) call rad_alloc(ps_info%chlocal,ps_info%chlocal%n)
       call gcopy(ps_info%chlocal%f,ps_info%chlocal%n)
       call gcopy(ps_info%chlocal%d2,ps_info%chlocal%n)
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

    subroutine read_header_tmp(n_orbnl, lmax_basis, n_pjnl, zval, unit)

      use global_module, ONLY: iprint_pseudo

      implicit none

      integer, intent(in)         :: unit
      integer, intent(out) :: n_orbnl, n_pjnl
      real(double), intent(out) :: zval

      character(len=78) :: line

      character(len=2)  :: symbol
      character(len=20) :: label
      integer :: z, lmax_basis, lmax_projs, zval_int
      real(double) :: mass, self_energy

      read(unit,'(a)') line
      if (trim(line) .eq. '<preamble>') then
         read(unit,'(a)') line
         do while(trim(line) .ne. '</preamble>')
            read(unit,'(a)') line
         end do
      endif

      if(iprint_pseudo>1) then
         read(unit,'(a2)') symbol
         write(*,fmt='(10x,"symbol = ",a2)') symbol
         read(unit,'(a20)') label
         write(*,fmt='(10x,"label = ",a20)') label
         read(unit,'(i5)') z
         write(*,fmt='(10x,"z = ",i5)') z
         read(unit,'(i5)') zval_int
         zval=real(zval_int,double)
         write(*,fmt='(10x,"zval_int, zval = ",i5,f12.5)') zval_int, zval
         read(unit,'(g25.15)') mass
         write(*,fmt='(10x,"mass = ",g25.15)') mass
         read(unit,'(g25.15)') self_energy
         write(*,fmt='(10x,"self_energy = ",g25.15)') self_energy
         read(unit,'(2i4)') lmax_basis, n_orbnl
         write(*,fmt='(10x,"lmax_basis, n_orbnl = ",2i4)') lmax_basis, n_orbnl
         read(unit,'(2i4)') lmax_projs, n_pjnl
         write(*,fmt='(10x,"lmax_projs, n_pjnl = ",2i4)') lmax_projs, n_pjnl
      else
         read(unit,'(a2)') symbol
         read(unit,'(a20)') label
         read(unit,'(i5)') z
         read(unit,'(i5)') zval_int
         zval=zval_int
         read(unit,'(g25.15)') mass
         read(unit,'(g25.15)') self_energy
         read(unit,'(2i4)') lmax_basis, n_orbnl
         read(unit,'(2i4)') lmax_projs, n_pjnl
      end if
    end subroutine read_header_tmp
  end subroutine read_ion_ascii_tmp

!!***
end module pseudo_tm_info
