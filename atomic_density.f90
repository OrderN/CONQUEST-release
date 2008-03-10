! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module atomic_density
! ------------------------------------------------------------------------------
! Code area 5: Self consistency
! ------------------------------------------------------------------------------

!!****h* Conquest/atomic_density *
!!  NAME
!!   atomic_density
!!  PURPOSE
!!   Holds data structures for storing atomic densities and routines for reading 
!!   or building these densities
!!  USES
!!
!!  AUTHOR
!!   D.R.Bowler and M.J.Gillan
!!  CREATION DATE
!!   24/09/2002 
!!  MODIFICATION HISTORY
!!   12:15, 25/09/2002 mjg & drb 
!!    Added rcut_dens - maximum atomic density table cutoff
!!   17:25, 2003/03/18 dave
!!    Added second derivative to radial density type and splining routine
!!   15:32, 08/04/2003 drb 
!!    Added debugging statements
!!   17:28, 2003/06/10 tm
!!    spline problem fixed
!!   2008/02/04 08:21 dave
!!    Changed for output to file not stdout
!!  SOURCE
!!
module atomic_density

  use datatypes
  use global_module, ONLY: io_lun

  implicit none
  save
  
  type radial_density
     integer :: length
     real(double) :: cutoff
     real(double), pointer, dimension(:) :: table
     real(double), pointer, dimension(:) :: d2_table
  end type radial_density

  type(radial_density), allocatable, dimension(:) :: atomic_density_table

  ! Maximum cutoff atomic on charge density tables
  real(double), allocatable, dimension(:) :: rcut_dens
  logical :: flag_atomic_density_from_pao
  character(len=80) :: read_atomic_density_file
  character(len=10) :: atomic_density_method

  ! RCS tag for object file identification
  character(len=80), private :: RCSid = "$Id$"
!!***

contains

! -----------------------------------------------------------
! Subroutine read_atomic_density
! -----------------------------------------------------------

!!****f* read_density/read_atomic_density *
!!
!!  NAME 
!!   read_atomic_density
!!  USAGE
!!   read_atomic_density(inode, ionode)
!!  PURPOSE
!!   Reads tables of atomic electron density
!!  INPUTS
!!   inode
!!   ionode
!!  USES
!!   datatypes, global_module, GenComms, numbers
!!  AUTHOR
!!   L.K.Dash
!!  CREATION DATE
!!   14/05/02
!!  MODIFICATION HISTORY
!!   18/06/2002 lkd
!!    changed "open/close" file  statements to io_assign/io_close
!!   17/07/2002 lkd
!!    moved to new density_read_module, changed to new radial_density derived type
!!   10/08/2002 MJG:
!!    Minor renaming
!!   11/08/2002 MJG:
!!    new scheme for broadcasting data
!!   11:23, 24/09/2002 mjg & drb 
!!    Minor reformatting
!!   11:53, 24/09/2002 mjg & drb 
!!    Included in atomic_density
!!   12:14, 25/09/2002 mjg & drb 
!!    Added rcut_dens to keep track of maximum cutoff on atomic charge density
!!   2006/09/20 17:19 dave
!!    Moved read of file to initial_read
!!  SOURCE
!!
  subroutine read_atomic_density(inode,ionode,n_species)

    use datatypes
    use global_module, ONLY: iprint_SC, area_SC
    use GenComms, ONLY: gcopy, cq_abort
    use numbers
    use memory_module, ONLY: reg_alloc_mem, type_dbl

    implicit none

    ! Passed variables
    integer, intent(in) :: inode, ionode, n_species

    ! Local variables
    character(len=80) :: read_atomic_density_file
    integer :: i, nt, n_sp ! do loop variable
    integer :: lun ! local unit number for io_assign
    integer :: ios ! status indicator for opening unit
    integer :: alls ! status indicator for allocations
    real(double) :: r_dummy

    if(inode == ionode) then
       if(iprint_SC >= 1) then
          write(unit=io_lun,fmt='(//10x,50("+")/10x,"read_atomic_density: &
               &no. of species for reading atomic density:",i5/10x,50("+"))') n_species
       end if

       if(n_species < 1) call cq_abort('read_atomic_density: &
            &no. of species must be positive',n_species)

       call io_assign(lun)
       if(iprint_SC >= 2) then
          write(unit=io_lun,fmt='(/10x," read_atomic_density: io_assign unit no. lun:",i3)') lun
       end if

       if(iprint_SC >= 2) then
          write(unit=io_lun,fmt='(/10x," read_atomic_density: name of input file:",a80)') &
               &read_atomic_density_file
       end if

       open(unit=lun,file=read_atomic_density_file,status='old',iostat=ios)
       if(ios /= 0) call cq_abort('read_atomic_density: failed to open input file')
    end if

    ! allocate memory for atomic density tables
    if(allocated(atomic_density_table)) then
       do i=1,size(atomic_density_table)
          deallocate(atomic_density_table(i)%table)
       end do
       deallocate(atomic_density_table)
    end if
    allocate(atomic_density_table(n_species),STAT=alls)
    if (alls /= 0) call cq_abort('read_atomic_density: error allocating atomic_density_table ',n_species)

    allocate(rcut_dens(n_species),STAT=alls)
    if (alls /= 0) call cq_abort('read_atomic_density: error allocating rcut_dens ',n_species)
    call reg_alloc_mem(area_SC, n_species, type_dbl)
    rcut_dens = 0.0_double
    do n_sp = 1, n_species

       if(iprint_SC >= 2) write(unit=io_lun,fmt='(/10x," reading atomic density data for species no:",i3)') n_sp
       if(inode == ionode) then
          read(lun,fmt=*) atomic_density_table(n_sp)%length
          read(lun,fmt=*) atomic_density_table(n_sp)%cutoff
          if(iprint_SC >= 2) then
             write(unit=io_lun,fmt='(/10x," table length:",i5)') atomic_density_table(n_sp)%length
             write(unit=io_lun,fmt='(/10x," radial cut-off distance:",f12.6)') atomic_density_table(n_sp)%cutoff
          end if
       end if

       call gcopy(atomic_density_table(n_sp)%length)
       call gcopy(atomic_density_table(n_sp)%cutoff)
       rcut_dens(n_sp) = atomic_density_table(n_sp)%cutoff
       !if(atomic_density_table(n_sp)%cutoff>rcut_dens) rcut_dens=atomic_density_table(n_sp)%cutoff

       allocate(atomic_density_table(n_sp)%table(atomic_density_table(n_sp)%length),stat = alls)
       if(alls /= 0) call cq_abort('read_atomic_density: &
            &failed to allocate atomic_density_table(n_sp)%table')
       call reg_alloc_mem(area_SC,atomic_density_table(n_sp)%length, type_dbl)

       if(atomic_density_table(n_sp)%length > 0) then

          if(inode == ionode) then
             do nt = 1, atomic_density_table(n_sp)%length
                read(unit=lun,fmt=*) r_dummy, atomic_density_table(n_sp)%table(nt)
             end do
          end if

          call gcopy(atomic_density_table(n_sp)%table,atomic_density_table(n_sp)%length)
       end if

    end do

    if(inode == ionode) then
       call io_close(lun)
       if(iprint_SC>2) then
          do n_sp = 1,n_species
            write(io_lun,fmt='(10x,"Atomic density cutoff for species ",i4," : ",f15.8)') n_sp,rcut_dens(n_sp)
         end do
      end if
   end if
  end subroutine read_atomic_density
!!***
  
! -----------------------------------------------------------
! Subroutine make_atomic_density_from_paos
! -----------------------------------------------------------

!!****f* atomic_density/make_atomic_density_from_paos *
!!
!!  NAME 
!!   make_atomic_density_from_paos
!!  USAGE
!! 
!!  PURPOSE
!!   Makes atomic densities from PAOs read in 
!!  INPUTS
!! 
!! 
!!  USES
!!   datatypes, GenComms, global_module, numbers, pao_format
!!  AUTHOR
!!   M.J.Gillan and D.R.Bowler
!!  CREATION DATE
!!   Summer 2002
!!  MODIFICATION HISTORY
!!   11:53, 24/09/2002 mjg & drb 
!!    Incorporated into atomic_density
!!   12:35, 25/09/2002 mjg & drb 
!!    Added rcut_dens to keep track of maximum cutoff on atomic charge density
!!   13:52, 29/07/2003 drb 
!!    Changed iprint level at which 2001 point table printed out...
!!   2007/11/16 10:19 dave
!!    Changed linear interpolation to spline interpolation to fix forces problem
!!   2008/03/03 18:40 dave
!!    Changed float to real()
!!  SOURCE
!!
  subroutine make_atomic_density_from_paos(inode,ionode,n_species)

    use datatypes
    use GenComms, ONLY : cq_abort, gcopy
    use global_module, ONLY : iprint_SC, area_SC
    use numbers, ONLY : zero, one, four, very_small, pi
    use pao_format
    use memory_module, ONLY: reg_alloc_mem, type_dbl
    use spline_module, ONLY: splint

    implicit none

    real(double), parameter :: one_over_four_pi = one/(four*pi)
    integer, intent(in) :: inode,ionode, n_species
    integer :: alls, i, lun, nt, n_am, n_sp, n_zeta
    integer, parameter :: default_atomic_density_length = 2001
    real(double) :: alpha, cutoff, density_deltar, pao_deltar, r, rn_am, val
    logical :: range_flag

    if(allocated(atomic_density_table)) then
       do i=1,size(atomic_density_table)
          deallocate(atomic_density_table(i)%table)
       end do
       deallocate(atomic_density_table)
    end if
    allocate(atomic_density_table(n_species),STAT = alls)
    if(alls /= 0) call cq_abort('make_atomic_density_from_paos: error allocating atomic_density_table ',n_species)

    allocate(rcut_dens(n_species),STAT=alls)
    call reg_alloc_mem(area_SC, n_species, type_dbl)
    if (alls /= 0) call cq_abort('make_atomic_density_from_paos: error allocating rcut_dens ',n_species)
    rcut_dens = 0.0_double
    do n_sp = 1, n_species
       ! By default, use max PAO cut-off radius as density cut-off
       cutoff = zero
       do n_am = 0, pao(n_sp)%greatest_angmom
          if(pao(n_sp)%angmom(n_am)%n_zeta_in_angmom > 0) then
             do n_zeta = 1, pao(n_sp)%angmom(n_am)%n_zeta_in_angmom
                cutoff = max(cutoff,pao(n_sp)%angmom(n_am)%zeta(n_zeta)%cutoff)
             end do
          end if
       end do
       atomic_density_table(n_sp)%cutoff = cutoff
       rcut_dens(n_sp)=atomic_density_table(n_sp)%cutoff
       !if(atomic_density_table(n_sp)%cutoff>rcut_dens) rcut_dens=atomic_density_table(n_sp)%cutoff
       if((inode == ionode).and.(iprint_SC >= 2)) &
            &write(unit=io_lun,fmt='(/10x," radial cut-off to be used is taken to be max PAO cut-off radius for &
            &this species:",f12.6)') atomic_density_table(n_sp)%cutoff     
       ! By default, use the parameter given above for length
       atomic_density_table(n_sp)%length = default_atomic_density_length
       ! Check for sensible length
       if(atomic_density_table(n_sp)%length < 2) &
            &call cq_abort('make_atomic_density_from_paos: table length must be >= 2',&
            &atomic_density_table(n_sp)%length)
       ! Allocate space for atomic density
       allocate(atomic_density_table(n_sp)%table(atomic_density_table(n_sp)%length),STAT = alls)
       if(alls /= 0) call cq_abort('make_atomic_density_from_paos: &
            &failed to allocate atomic_density_table(n_sp)%table()')
       call reg_alloc_mem(area_SC,atomic_density_table(n_sp)%length, type_dbl)

       ! initiate to zero table of atomic density for current species and calculate table spacing
       do nt = 1, atomic_density_table(n_sp)%length
          atomic_density_table(n_sp)%table(nt) = zero
       end do
       ! Find spacing of table
       density_deltar = atomic_density_table(n_sp)%cutoff/&
            &real(atomic_density_table(n_sp)%length-1,double)
       ! Write out info and check angular momentum
       if((inode == ionode).and.(iprint_SC >= 2)) &
            write(unit=io_lun,fmt='(/10x," greatest ang. mom. for making density from PAOs:",i3)') &
            &pao(n_sp)%greatest_angmom
       if(pao(n_sp)%greatest_angmom < 0) &
            &call cq_abort('make_atomic_density_from_paos: greatest ang. mom. cannot be negative')
       ! Loop over angular momenta
       do n_am = 0, pao(n_sp)%greatest_angmom
          ! Check for zeta
          if(pao(n_sp)%angmom(n_am)%n_zeta_in_angmom > 0) then
             ! Loop over zetas
             do n_zeta = 1, pao(n_sp)%angmom(n_am)%n_zeta_in_angmom
                if(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length > 1) then
                   pao_deltar = pao(n_sp)%angmom(n_am)%zeta(n_zeta)%cutoff/&
                        &real(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length-1,double)
                   do nt = 1, atomic_density_table(n_sp)%length
                      r = (nt-1)*density_deltar
                      i = one + very_small + r/pao_deltar
                      if(i+1 <= pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length) then
                         if(n_am /=0) then
                           rn_am = r**n_am
                         else
                           rn_am = one
                         endif
                         call splint(pao_deltar,pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(:), &
                              pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table2(:), &
                              pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length, &
                              r,val,range_flag)
                         atomic_density_table(n_sp)%table(nt) = atomic_density_table(n_sp)%table(nt) + &
                              &one_over_four_pi * pao(n_sp)%angmom(n_am)%occ(n_zeta) * &
                              &(rn_am * val )**2
                         ! Linear interpolation
                         !alpha = one + r/pao_deltar - i
                         !atomic_density_table(n_sp)%table(nt) = atomic_density_table(n_sp)%table(nt) + &
                         !     &one_over_four_pi * pao(n_sp)%angmom(n_am)%occ(n_zeta) * &
                         !     &(rn_am * ((one-alpha)*pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(i) + &
                         !     &alpha*pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(i+1)))**2
                      end if ! if(i+1<=pao(...)%length
                   end do ! do nt = atomic_density_table()%length
                end if ! if(pao(...)%length > 1
             end do ! do n_zeta = pao(...)%n_zeta_in_angmom
          end if ! pao(...)%n_zeta_in_angmom > 0
       end do ! n_am = pao(...)%greatest_angmom
       do nt = 1, atomic_density_table(n_sp)%length
          r = (nt-1)*density_deltar
          if(inode==ionode.AND.iprint_SC>3) write(io_lun,fmt='(10x,"Radial table: ",i5,2f15.8)') &
               nt,r,atomic_density_table(n_sp)%table(nt)
       end do
    end do ! n_sp = n_species
    do n_sp = 1,n_species
       if(inode == ionode.AND.iprint_SC>2) &
            write(io_lun,fmt='(10x,"Atomic density cutoff for species ",i4," : ",f15.8)') n_sp,rcut_dens(n_sp)
    end do
  end subroutine make_atomic_density_from_paos
!!***

! -----------------------------------------------------------
! Subroutine spline_atomic_density
! -----------------------------------------------------------

!!****f* atomic_density/spline_atomic_density *
!!
!!  NAME 
!!   spline_atomic_density
!!  USAGE
!! 
!!  PURPOSE
!!   Build spline tables for the radial tables of atomic densities
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   17:18, 2003/03/18 dave
!!  MODIFICATION HISTORY
!!   17:27, 2003/06/10 tm
!!    Fixed over-run problem with splining
!!  SOURCE
!!
  subroutine spline_atomic_density(n_species)

    use datatypes
    use numbers
    use spline_module, ONLY: spline,splint,dsplint
    use GenComms, ONLY: cq_abort, inode, ionode
    use global_module, ONLY: iprint_SC, area_SC
    use memory_module, ONLY: reg_alloc_mem, type_dbl

    implicit none

    ! Passed variables
    integer :: n_species
    ! Local variables

    integer :: i, j, n, n_l, stat

    real(double) :: d_end, d_origin, delta_r, r, local_density, derivative
    logical :: range_flag

    ! loop over species and do the interpolation
    do n=1, n_species
       allocate(atomic_density_table(n)%d2_table(atomic_density_table(n)%length), STAT=stat)
       if(stat/=0) call cq_abort('spline_atomic_density: error allocating d2_table ! ',stat)
       call reg_alloc_mem(area_SC,atomic_density_table(n)%length, type_dbl)
       ! do the splining for the table
       delta_r = atomic_density_table(n)%cutoff/real(atomic_density_table(n)%length-1,double)
       d_origin = (atomic_density_table(n)%table(2)- atomic_density_table(n)%table(1))/delta_r
       d_end = (atomic_density_table(n)%table(atomic_density_table(n)%length)- &
            atomic_density_table(n)%table(atomic_density_table(n)%length-1))/delta_r
       call spline( atomic_density_table(n)%length, delta_r, atomic_density_table(n)%table(:),  &
            d_origin, d_end, atomic_density_table(n)%d2_table(:) )
       if(inode==ionode.AND.iprint_SC>3) then
          write(io_lun,fmt='(10x,"Atomic density for species ",i5)') n
          do i=1,atomic_density_table(n)%length
             r = real(i-1,double)*delta_r
             if(r < atomic_density_table(n)%cutoff) then   !TM
!                call splint(delta_r,atomic_density_table(n)%table(:), & 
!                     atomic_density_table(n)%d2_table(:), &
!                     atomic_density_table(n)%length, & 
!                     r,local_density,range_flag)
                call dsplint(delta_r,atomic_density_table(n)%table(:), & 
                     atomic_density_table(n)%d2_table(:), &
                     atomic_density_table(n)%length, & 
                     r,local_density,derivative,range_flag)
                write(io_lun,fmt='(10x,3f20.12)') r,local_density,derivative
             else                                                    !TM
                write(io_lun,fmt='(10x,"r in spline_atomic_density > cutoff  ",2f20.12)') &
                     r, atomic_density_table(n)%cutoff
                write(io_lun,fmt='(10x,3f20.12)') r, &                                   !TM
                     atomic_density_table(n)%table(atomic_density_table(n)%length), & !TM
                     atomic_density_table(n)%d2_table(atomic_density_table(n)%length) !TM
             endif                                                   !TM
          end do
       end if
    end do
    return
  end subroutine spline_atomic_density
!!***

end module atomic_density
