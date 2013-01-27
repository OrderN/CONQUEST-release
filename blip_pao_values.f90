! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module blip_pao_values
! ------------------------------------------------------------------------------
! Code area 11: Basis operations
! ------------------------------------------------------------------------------

!!****h* Conquest/blip_pao_values *
!!  NAME
!!   blip_pao_values
!!  PURPOSE
!!   Sole purpose is to contain the four routines b_value,
!!   f0_value, f1_value and f2_value. The purpose of these
!!   routines is to calculate contributions to the vector
!!   of overlap integrals of blip functions with PAO's.
!!   Purposes of the four routines are described in more
!!   detail below.
!!  USES
!!
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!   12/08/2002 mjg
!!    i)  All pao tables changed so that enumeration starts from one not zero
!!    ii) All radial tables changed so that the expectation is that
!!        they are divided by r**l (in line with Siesta radial tables
!!        and spherical harmonic conventions)
!!   2008/06/01 ast
!!    Added timers
!!  SOURCE
module blip_pao_values

  use global_module, only: area_basis

  implicit none
  save

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: &
       RCSid = "$Id$"
!!***

contains

! -----------------------------------------------------------
! Subroutine b_value
! -----------------------------------------------------------

!!****f* blip_pao_values/b_value *
!!
!!  NAME 
!!   b_value
!!  USAGE
!! 
!!  PURPOSE
!!   Calculates values of the one-dimensional B-spline on the
!!   integration grid used to calculate the overlap integral
!!   of blip functions and PAO's. The calculated values are
!!   stored in the array bv(). The table of values is stored
!!   in four pieces, specified by the four values of
!!   second array dimension of bv(). The four pieces correspond
!!   to the four blip-grid intervals covered by each B-spline.
!! 
!!  INPUTS
!!   nu_int: number of integration-grid intervals in each
!!     blip-grid interval. (NB: this integration grid is
!!     for private use of the PAO-to-blip conversion, and
!!     has no connection with the integration grid used
!!     elsewhere in Conquest.)
!! 
!!  USES
!! 
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!   2008/03/03 18:42 dave
!!    Changed float to real
!!  SOURCE
  subroutine b_value(nu_int,bv)
    use datatypes
    use numbers
    implicit none
    integer, intent(in) :: nu_int
    real(double), intent(out), dimension(:,:) :: bv
    integer :: nu
    real(double) :: deltax, x, f1, f2

    deltax = one/real(nu_int,double)

    do nu = 1, nu_int-1
       x = nu*deltax
       f1 = one - three_halves*x*x + three_quarters*x*x*x
       f2 = quarter*x*x*x
       bv(nu,1) = f2
       bv(nu_int-nu,2) = f1
       bv(nu,3) = f1
       bv(nu_int-nu,4) = f2
    end do
    bv(nu_int,1) = quarter
    bv(nu_int,2) = one
    bv(nu_int,3) = quarter
    bv(nu_int,4) = zero

    return

  end subroutine b_value
!!*** 

! -----------------------------------------------------------
! Subroutine f0_value
! -----------------------------------------------------------
!!****f* blip_pao_values/f0_value *
!!
!!  NAME 
!!   f0_value
!!  USAGE
!! 
!!  PURPOSE
!!   Evaluates pao for a  given value of position
!! (x,y,z), which is supplied in argument list. In production
!!   use, pao is obtained by interpolation in the appropriate
!!   radial table of pao values. For test purposes, there is
!!   also provision in the sbrt for inserting code to evaluate
!!   pao represented as a superposition of blips, or in other ways.
!!   Production use is obtained when the flag i_test is set equal
!!   to zero. To use the test version, flag i_test must be re-set
!!   to unity and the routine must be recompiled. At the time
!!   that this routine is first put into Conquest, code for
!!   doing tests has been removed. However, this notional
!!   facility has been left in place in order to enable
!!   test code to be put in should it be needed.
!!
!!  INPUTS
!!   sym_type: specifies symmetry type in current use.
!!   n_sp: species number
!!   n_am: angular momentum number (s = 0, p = 1, d = 2)
!!   n_zeta: zeta number
!!   x, y, z: cartesian components of position where PAO is
!!   to be evaluated (atomic units of length).
!!   fv0: output calculated value of PAO.
!!
!!  USES
!! 
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!   12/08/2002 mjg
!!    i) All pao tables changed so that enumeration starts from one not zero
!!    ii) All radial tables changed so that the expectation is that they are divided by
!!    r**l (in line with Siesta radial tables and spherical harmonic conventions)
!!  SOURCE
  subroutine f0_value(sym_type,n_sp,n_am,n_zeta,x,y,z,fv0)
    use datatypes
    use numbers
    use pao_format
    use GenComms, ONLY: cq_abort

    implicit none

    character(len=*), intent(in) :: sym_type
    integer, intent(in) :: n_sp, n_am, n_zeta
    real(double), intent(in) ::  x, y, z
    real(double), intent(out) :: fv0
    integer :: i_test, j0
    real(double) :: alpha, ang_fac, deltar, r
    real(double) :: a1g_norm, t1u_norm, t2g_norm, eg_norm

    ! for routine use, i_test = 0; for test, i_test = 1
    i_test = 0

    select case(i_test)

    ! code for routine use -----------------------------------------
    case(0)

       ! angular normalisation constants
       a1g_norm = sqrt(one/(four*pi))
       t1u_norm = sqrt(three/(four*pi))
       t2g_norm = sqrt(fifteen/(four*pi))
       eg_norm = sqrt(fifteen/(sixteen*pi))

       ! prepare to consult appropriate table
       deltar = pao(n_sp)%angmom(n_am)%zeta(n_zeta)%cutoff/&
            &(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length - 1)

       r = sqrt(x*x + y*y + z*z)
       j0 = 1 + int(r/deltar)
       if(j0+1 > pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length) then
          fv0 = zero
       else
          select case(sym_type)
          case("a1g")
             ang_fac = a1g_norm
          case("t1u")
             ang_fac = t1u_norm*z
          case("t2g")
             ang_fac = t2g_norm*x*y
          case("eg")
             ang_fac = eg_norm*(x*x - y*y)
          case default
             call cq_abort('f0_value: unrecognised symmetry type')
          end select
          alpha = one + r/deltar - j0
          fv0 = ang_fac*((one - alpha)*pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(j0) + &
               &alpha*pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(j0+1))
       end if

       ! code for test use ----------------------------------------------
    case(1)
       call cq_abort('f0_value: entered dummy part of code')

    end select

    return

  end subroutine f0_value
!!*** 

! -----------------------------------------------------------
! Subroutine f1_value
! -----------------------------------------------------------

!!****f* blip_pao_values/f1_value *
!!
!!  NAME 
!!   f1_value
!!  USAGE
!! 
!!  PURPOSE
!!   Evaluates PAO for a sequence of x-values for
!!   given constant values of y and z. The number of x-values
!!   in the sequence is nu_int, the number of integration-grid
!!   intervals in the blip-grid interval. The value of the PAO
!!   for each triplet (x,y,z) is obtained by calling the routine
!!   f0_value
!! 
!!  INPUTS
!!   sym_type: specifies symmetry type required.
!!   n_sp: species number.
!!   n_am: angular momentum number (s = 0, p = 1, d = 2).
!!   n_zeta: zeta number.
!!   nu_int: number of integration-grid intervals in each
!!   blip-grid interval.
!!   delta_ig: length of integration-grid interval (atomic units).
!!   x0: initial x-point in the sequence (atomic units).
!!   y, z: fixed y- and z-values (atomic units).
!!  USES
!! 
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!
!!  SOURCE
  subroutine f1_value(sym_type,n_sp,n_am,n_zeta,nu_int,&
       &delta_ig,x0,y,z,fv1)
    use datatypes
    implicit none

    character(len=*), intent(in) :: sym_type
    integer, intent(in) :: n_sp, n_am, n_zeta, nu_int
    real(double), intent(in) :: delta_ig, x0, y, z
    real(double), intent(out), dimension(:) :: fv1
    integer ::  n
    real(double) :: fv0, x

    ! loop over x values
    do n = 1, nu_int
       x = x0 + n*delta_ig
       call f0_value(sym_type,n_sp,n_am,n_zeta,x,y,z,fv0)
       fv1(n) = fv0
    end do

    return

  end subroutine f1_value
!!*** 

! -----------------------------------------------------------
! Subroutine f2_value
! -----------------------------------------------------------

!!****f* blip_pao_values/f2_value *
!!
!!  NAME 
!!   f2_value
!!  USAGE
!! 
!!  PURPOSE
!!   For fixed z, and for a sequence of y values
!!   covering one blip interval, evaluates all scalar products
!!   over x.
!! 
!!  INPUTS
!!   Essentially the same as for sbrt f0_value (see above).
!! 
!!  USES
!! 
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!   2008/03/03 18:42 dave
!!    Changed float to real
!!  SOURCE
  subroutine f2_value(sym_type,n_sp,n_am,n_zeta,nu_int,&
       &delta_ig,n_blip,bv,y0,z,fv2)
    use datatypes
    use numbers
    use GenComms, only: cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    implicit none

    character(len=*), intent(in) :: sym_type
    integer, intent(in) :: nu_int, n_blip, n_sp, n_am, n_zeta
    real(double), intent(in) :: delta_ig, y0, z
    real(double), intent(in), dimension(:,:) :: bv
    real(double), intent(out), dimension(:,:) :: fv2
    integer :: i, i_store, j, j_add, j_store, kb_min, kb_max, &
         m, n, n_b_int_in_diam, n_b_int_in_radius, stat
    real(double) :: b_int, deltax, sum, x0, y
    real(double), dimension(:), allocatable :: fv1
    real(double), dimension(4,4) :: store

    allocate(fv1(nu_int), STAT=stat)
    if (stat /= 0) call cq_abort("f2_value: Error alloc mem: ", nu_int)
    call reg_alloc_mem(area_basis, nu_int, type_dbl)

    ! check that n_blip is odd
    if(2*(n_blip/2) == n_blip) then
       call cq_abort('f2_value: n_blip must be odd',n_blip)
    end if

    ! constants
    n_b_int_in_diam = n_blip + 3
    n_b_int_in_radius = n_b_int_in_diam/2
    deltax = one/real(nu_int,double)
    b_int = nu_int*delta_ig

    ! compute and assemble contributions to scalar products
    do m = 1, nu_int
       y = y0 + m*delta_ig
       do i = 1, n_b_int_in_diam
          x0 = b_int*(i-1-n_b_int_in_radius)
          call f1_value(sym_type,n_sp,n_am,n_zeta,nu_int,&
               &delta_ig,x0,y,z,fv1)
          kb_min = 4 - min(n_b_int_in_diam-i,3)
          kb_max = min(i,4)
          i_store = 1 + mod(i-1,4)
          do j = kb_min, kb_max
             j_store = 1 + mod(i_store-j+4,4)
             sum = zero
             do n = 1, nu_int
                sum = sum + bv(n,j)*fv1(n)
             end do
             store(i_store,j_store) = deltax*sum
          end do
          if(i >= 4) then
             j_add = 1 + mod(i_store,4)
             fv2(m,i-3) = store(1,j_add) + store(2,j_add) + &
                  & store(3,j_add) + store(4,j_add)
          endif
       end do
    end do

    deallocate(fv1, STAT=stat)
    if (stat /= 0) call cq_abort("f2_value: Error dealloc mem")
    call reg_dealloc_mem(area_basis, nu_int, type_dbl)

    return

  end subroutine f2_value
!!*** 

! -----------------------------------------------------------
! Subroutine f3_value
! -----------------------------------------------------------

!!****f* blip_pao_values/f3_value *
!!
!!  NAME 
!!   f3_value
!!  USAGE
!! 
!!  PURPOSE
!!   For a sequence of z-values, evaluates all scalar products over y.
!! 
!!  INPUTS
!!    Essentially the same as for sbrt's f1_value and f2_value
!! 
!!  USES
!! 
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!   2008/03/03 18:42 dave
!!    Changed float to real
!!   2008/06/01 ast
!!    Added timers
!!  SOURCE
  subroutine f3_value(sym_type, n_sp, n_am, n_zeta, nu_int, delta_ig, &
                      n_blip, bv, z0, fv3)
    use datatypes
    use numbers
    use GenComms,               only: cq_abort
    use timer_stdclocks_module, only: start_timer, stop_timer, &
                                      tmr_std_basis
    use memory_module,          only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    character(len=*), intent(in) :: sym_type
    integer, intent(in) :: nu_int, n_blip, n_sp, n_am, n_zeta
    real(double), intent(in) :: delta_ig, z0
    real(double), intent(in), dimension(:,:) :: bv
    real(double), intent(out), dimension(:,:,:) :: fv3
    integer :: i, i_store, j, j_add, j_store, k, kb_min, kb_max, &
         n, n_b_int_in_radius, n_b_int_in_diam, p
    real(double) :: b_int, deltax, sum, y0, z
    real(double), dimension(:,:),   allocatable :: fv2
    real(double), dimension(:,:,:), allocatable :: store
    integer :: stat

    call start_timer(tmr_std_basis)

    allocate(fv2(nu_int,n_blip), store(4,4,n_blip), STAT=stat)
    if (stat /= 0) call cq_abort("f3_value: Error alloc mem: ", nu_int, n_blip)
    call reg_alloc_mem(area_basis, (nu_int+16)*n_blip, type_dbl)

    ! check that n_blip is odd
    if(2*(n_blip/2) == n_blip) then
       call cq_abort('f3_value: n_blip must be odd',n_blip)
    end if

    ! constants
    n_b_int_in_diam = n_blip + 3
    n_b_int_in_radius = n_b_int_in_diam/2
    deltax = one/real(nu_int,double)
    b_int = nu_int*delta_ig

    ! compute and assemble contributions to scalar products
    do k = 1, nu_int
       z = z0 + k*delta_ig
       do i = 1, n_b_int_in_diam
          y0 = b_int*(i-1-n_b_int_in_radius)
          call f2_value(sym_type,n_sp,n_am,n_zeta,nu_int,&
               &delta_ig,n_blip,bv,y0,z,fv2)
          kb_min = 4 - min(n_b_int_in_diam-i,3)
          kb_max = min(i,4)
          i_store = 1 + mod(i-1,4)
          do j = kb_min, kb_max
             j_store = 1 + mod(i_store-j+4,4)
             do p = 1, n_blip
                sum = zero
                do n = 1, nu_int
                   sum = sum + bv(n,j)*fv2(n,p)
                end do
                store(i_store,j_store,p) = deltax*sum
             end do
          end do
          if(i >= 4) then
             j_add = 1 + mod(i_store,4)
             do p = 1, n_blip
                fv3(k,p,i-3) = store(1,j_add,p) + store(2,j_add,p) + &
                     & store(3,j_add,p) + store(4,j_add,p)
             end do
          end if
       end do
    end do

    deallocate(fv2, store, STAT=stat)
    if (stat /= 0) call cq_abort("f3_value: Error dealloc mem")
    call reg_dealloc_mem(area_basis, (nu_int+16)*n_blip, type_dbl)

    call stop_timer(tmr_std_basis)

    return

  end subroutine f3_value
!!*** 

end module blip_pao_values
