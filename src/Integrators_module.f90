! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module Integrators_module
! ------------------------------------------------------------------------------
! Code area : 
! ------------------------------------------------------------------------------

!!****h* Conquest/Integrators_module
!! NAME
!!   Integrators_module
!! PURPOSE
!!   Evolves particle positions and velocities in MD loop
!! AUTHOR
!!   M.Arita
!! CREATION DATE
!!   2014/02/03
!! MODIFICATION HISTORY
!!   2015/06/19 15:25 dave
!!    Including FIRE from SA and COR with modifications
!!   2017/10/20 09:20 dave
!!    Moved fire parameters here from global, and added user parameter fire_max_step
!!   2017/11/10 15:16 dave
!!    Added fire_N_below_thresh
!!   2017/11/14 tsuyoshi
!!    Removed atom_coord_diff
!! SOURCE
!!
module Integrators

  use datatypes
  implicit none

  ! FIRE relaxation method
  ! Parameters adjusted after this many steps if there is slow convergence
  integer :: fire_N_max, fire_N_min, fire_N_below_thresh
  real(double) :: fire_alpha0, fire_f_inc, fire_f_dec, fire_f_alpha, fire_max_step
  
contains

  !****f* Integrators/vVerlet_r_dt
  ! PURPOSE
  !   Evolve particle positions by dt via velocity Verlet
  !    x(t+dt) = x(t) + dt*v(t) + (dt*dt/2)*a(t)
  !            = x(t) + dt*v(t+dt/2)
  ! USAGE
  !   call vVerlet_r_dt(dt,v,flag_movable)
  ! INPUTS
  !   double dt: MD time step
  !   double v : velocities at t+dt/2
  !   logical flag_movable: flag to tell if atoms move
  ! AUTHOR
  !   M.Arita
  ! CREATION DATE
  !   2014/02/03
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine vVerlet_r_dt(dt,v,flag_movable)
   ! Module usage
   use global_module, ONLY: ni_in_cell,id_glob,x_atom_cell,y_atom_cell,z_atom_cell, &
                            iprint_MD, io_lun
   use GenComms, only: myid
   use species_module, ONLY: species,mass
   use move_atoms, ONLY: fac

   implicit none
   ! passed variables
   real(double),intent(in) :: dt
   real(double),intent(in) :: v(3,ni_in_cell)
   logical,intent(in) :: flag_movable(3*ni_in_cell)
   ! local variables
   integer :: atom,speca,gatom,k,ibeg_atom
   real(double) :: massa
   logical :: flagx,flagy,flagz
   real(double) :: dx, dy, dz

   if (myid==0 .and. iprint_MD>1) write(io_lun,*) "Welcome to vVerlet_r_dt"

   ibeg_atom = 1
   do atom = 1, ni_in_cell
     speca = species(atom)
     massa = mass(speca)*fac
     gatom = id_glob(atom)
     !do k = 1, 3
     !  if (.NOT.flag_movable(ibeg_atom+k-1)) cycle
     !enddo
     flagx = flag_movable(ibeg_atom)
     flagy = flag_movable(ibeg_atom+1)
     flagz = flag_movable(ibeg_atom+2)
     ibeg_atom = ibeg_atom + 3
     ! X
     if (flagx) then
       dx = dt*v(1,atom)
       x_atom_cell(atom) = x_atom_cell(atom) + dx
     endif
     ! Y
     if (flagy) then
       dy = dt*v(2,atom)
       y_atom_cell(atom) = y_atom_cell(atom) + dy
     endif
     ! Z
     if (flagz) then
       dz = dt*v(3,atom)
       z_atom_cell(atom) = z_atom_cell(atom) + dz
     endif
   enddo

   return
  end subroutine vVerlet_r_dt
  !*****

  !!****f* Integrators/vVerlet_r_dt
  !!  PURPOSE
  !!    Evolve particle velocities by dt/2 via velocity Verlet
  !!     v(t+dt/2) = v(t) + (dt/2)*a(t)
  !!    When quenched-MD applies, calculate inner product v*f
  !!    If v*f < 0, set v = 0
  !! 
  !!    Note that the initial velocity is defined as v(0), NOT v(-dt/2)
  !!    any longer as in velocityVerlet at move_atoms.module
  !!  USAGE
  !!    call vVerlet_v_dthalf(dt,v,f,flag_movable,second_call)
  !!  INPUTS
  !!    double dt           : MD time step
  !!    double v            : velocities at t+dt/2
  !!    double force        : forces at t
  !!    logical flag_movable: flag to tell if atoms are movable
  !!    logical second_call : tell this is the 2nd call
  !!  OUTPUT
  !!    real(double), v: half-a-step evolved particle velocities
  !!  AUTHOR
  !!    M.Arita
  !!  CREATION DATE
  !!    2014/02/03
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !! 
  subroutine vVerlet_v_dthalf(dt,v,f,flag_movable,second_call)
   ! Module usage
   use numbers, ONLY: half,zero, one
   use global_module, ONLY: ni_in_cell,id_glob,flag_quench_MD, flag_MDdebug, io_lun, iprint_MD
   use species_module, oNLY: species,mass
   use move_atoms, ONLY: fac
   use GenComms, only: inode,ionode

   implicit none
   ! passed variables
   real(double),intent(in) :: dt
   real(double),dimension(3,ni_in_cell),intent(in)    :: f
   real(double),dimension(3,ni_in_cell),intent(inout) :: v
   logical,dimension(3*ni_in_cell),intent(in)         :: flag_movable
   logical,optional :: second_call
   ! local variables
   integer :: atom,speca,gatom,k,ibeg_atom
   real(double) :: vf,massa,ff,scale_fac
   logical :: flag_first

   flag_first=.true.
   if(present(second_call)) flag_first = .false.
   if (inode==ionode .and. flag_MDdebug) write(io_lun,*) "Welcome to vVerlet_v_dthalf"
   ibeg_atom=1
   ! For quenched-MD
   ! 2018/07/16 DRB It seems that zeroing velocities after either call is more efficient
   if(flag_quench_MD.and.flag_first) then
      ! Alternative implementation: update v
      do atom = 1, ni_in_cell
         speca = species(atom)
         massa = fac ! MDMin routine sets all masses = 1
         gatom = id_glob(atom)
         do k = 1, 3
            !if (.NOT.flag_movable(ibeg_atom+k-1)) cycle
            v(k,atom) = v(k,atom) + half*dt*f(k,gatom)/massa
         enddo
         ibeg_atom = ibeg_atom + 3
      enddo
      ! Accumulate system-wide F.v and set v=0 if F.v<0
      vf = zero
      ff = zero
      do atom = 1, ni_in_cell
         gatom = id_glob(atom)
         do k = 1, 3
            !if (.NOT.flag_movable(ibeg_atom+k-1)) cycle
            vf = vf + v(k,atom)*f(k,gatom)
            ff = ff + f(k,gatom)*f(k,gatom)
         enddo
         ibeg_atom = ibeg_atom + 3
      enddo
      scale_fac = vf/ff
      ibeg_atom=1
      ! Now zero or update
      if(vf<zero) then
         if(inode==ionode.AND.iprint_MD>2) write(io_lun,fmt='(2x,"In quenched MD, zeroing velocities")')
         v(:,:) = zero
      else
         do atom = 1, ni_in_cell
            speca = species(atom)
            massa = fac ! MDMin routine sets all masses = 1
            gatom = id_glob(atom)
            do k = 1, 3
               !if (.NOT.flag_movable(ibeg_atom+k-1)) cycle
               ! Original
               !v(k,atom) = v(k,atom) + half*dt*f(k,gatom)/massa
               ! Alternative ! Hack by DRB - choose proportion of force direction to add
               v(k,atom) = (one - fire_alpha0)*v(k,atom) + fire_alpha0*f(k,gatom)*scale_fac
            enddo
            ibeg_atom = ibeg_atom + 3
         enddo
      end if
      ! Alternative implementation: update v
      do atom = 1, ni_in_cell
         speca = species(atom)
         massa = fac ! MDMin routine sets all masses = 1
         gatom = id_glob(atom)
         do k = 1, 3
            !if (.NOT.flag_movable(ibeg_atom+k-1)) cycle
            v(k,atom) = v(k,atom) + half*dt*f(k,gatom)/massa
         enddo
         ibeg_atom = ibeg_atom + 3
      enddo
   ! Ordinary MD update
   else if(.NOT.flag_quench_MD) then
      do atom = 1, ni_in_cell
         speca = species(atom)
         massa = mass(speca)*fac
         if(flag_quench_MD) massa = fac
         gatom = id_glob(atom)
         do k = 1, 3
            v(k,atom) = v(k,atom) + half*dt*f(k,gatom)/massa
            if (.NOT.flag_movable(ibeg_atom+k-1)) v(k,atom) = zero
         enddo
         ibeg_atom = ibeg_atom + 3
      enddo
    !! Use global labelling
    !   do gatom = 1, ni_in_cell
    !     speca = species_glob(gatom)
    !     massa = mass(speca)*fac
    !     do k = 1,3
    !      atom_vels(k,gatom) = atom_vels(k,gatom) + half*dt*f(k,gatom)/massa
    !      if(.NOT.flag_move_atom(k,gatom)) atom_vels(k,gatom) = zero
    !     enddo
    !   enddo
   endif

   return
  end subroutine vVerlet_v_dthalf
  !!***
  
!!****f* Integrators/fire_qMD
!! PURPOSE
!!   performs structure optimisation via (modified) FIRE quenched MD
!!   Phys. Rev. Lett. 97, 170201 (2006).
!!   1. velocities are updated according to:
!!      v -> v + dt * f (m = 1 for all atomic species)
!!   2. the product P = v * f is calculated:
!!      if P > 0 , v -> v'
!!      else v -> dt * f  
!!      MD time step & other paremeters of the method are
!!      modified accordingly
!!   3. r -> r + dt * v'
!!   4. I have noticed a slow convergence for large systems
!!      ( > 10^4 atoms) and introduced an extra control parameter 
!!      for adjusting velocities, MD time step & other FIRE parameters:
!!      when convergence is slow P(t+dt)/(P(t)*dt) oscillates about a constant 
!!      value (dt is usually const. thus P almost does not change). If P
!!      does not change for fire_N_max steps, velocities & FIRE parameters are 
!!      changed as for the P <= 0 case. 
!!      fire_N_max is set by 'AtomMove.FireNMaxSlowQMD', which is 10 steps by default
!! USAGE
!!   call fire_qMD(dt,v,f,flag_movable,iter)
!!   to use FIRE set 'AtomMove.FIRE = T' in Conquest_input
!! INPUTS
!!   double dt           : MD time step
!!   double v            : velocities at t+dt/2
!!   double force        : forces at t
!!   logical flag_movable: flag to tell if atoms are movable
!!   integer iter        : MD step
!! OUTPUT
!!   double dt           : changed MD time step
!!   real(double), v     : changed particle velocities
!! AUTHOR
!!   S. Arapan
!! CREATION DATE
!!   2014/07/31
!! MODIFICATION HISTORY
!!   2015/02/04          : modifying in order to allow continuation of previous FIRE run
!!   2015/06/19 15:27 dave
!!    Adapted to use global module for parameters read from input
!!   2017/10/20 09:21 dave
!!    Removed fire_step_max argument (now module variable and user input fire_max_step)
!! SOURCE
!!
  subroutine fire_qMD(dt,v,f,flag_movable,iter,fire_N,fire_N2,fire_P0,fire_alpha)

    ! Module usage

    use numbers, ONLY: half,zero
    use global_module, ONLY: ni_in_cell,id_glob,x_atom_cell,y_atom_cell,z_atom_cell, &
         flag_move_atom, io_lun,  iprint_MD
    use GenComms, ONLY: myid
    use io_module, ONLY: write_fire

    implicit none

    ! passed variables
    real(double),intent(inout) :: dt, fire_P0, fire_alpha
    real(double),dimension(3,ni_in_cell),intent(in)    :: f
    real(double),dimension(3,ni_in_cell),intent(inout) :: v
    logical,dimension(3*ni_in_cell),intent(in)         :: flag_movable
    integer :: iter
    integer,intent(inout) :: fire_N, fire_N2
    
    ! local variables
    integer :: atom,speca,gatom,k,ibeg_atom
    integer :: lun
    real(double) :: vf,massa
    logical :: flagx,flagy,flagz
    real(double) :: dx, dy, dz

    ! modified FIRE quenched MD
    real(double) :: fire_P, fire_norm_F, fire_norm_v, fire_r
    real(double) :: fire_r1, fire_r2

    fire_r1 = 0.75_double / fire_max_step
    fire_r2 = 1.25_double / fire_max_step

    ! update velocities
    ! for FIRE quench MD all atomic masses are set to unity
    do atom = 1, ni_in_cell
       gatom = id_glob(atom)
       do k = 1, 3
          v(k,atom) = v(k,atom) + dt * f(k,gatom)
       end do
    end do


    ! FIRE quenched MD
    fire_P = 0.0_double
    fire_norm_F = 0.0_double
    fire_norm_v = 0.0_double
    do atom = 1, ni_in_cell
       gatom = id_glob(atom)
       do k = 1, 3
          fire_P = fire_P + v(k,atom) * f(k,gatom)
          fire_norm_F = fire_norm_F + f(k,gatom)**2
          fire_norm_v = fire_norm_v + v(k,atom)**2
       end do
    end do
    fire_norm_F = sqrt(fire_norm_F)
    fire_norm_v = sqrt(fire_norm_v)
    fire_norm_v = fire_norm_v / fire_norm_F
    fire_r = fire_P / (dt * fire_P0)
    fire_P0 = fire_P

    ! Power, F.v, greater than zero means we're still going down hill
    if (fire_P > 0.0_double) then
       fire_N = fire_N + 1
       ! Check for slow convergence
       if (fire_r > fire_r1 .AND. fire_r < fire_r2) then
          fire_N2 = fire_N2 + 1
       else
          fire_N2 = 0
       end if
       ! Update velocities
       do atom = 1, ni_in_cell
          gatom = id_glob(atom)
          do k = 1, 3
             v(k,atom) = (1 - fire_alpha) * v(k,atom) + &
                  fire_alpha * fire_norm_v * f(k,gatom)
          end do
       end do
       ! Update FIRE parameters
       if (fire_N > fire_N_min) then
          ! Slow convergence
          if (fire_N2 > fire_N_max) then
             dt = dt * fire_f_dec ! 2015/01/29 SA
             do atom = 1, ni_in_cell
                gatom = id_glob(atom)
                do k = 1, 3
                   ! This is the normal FIRE approach
                   !v(k,atom) = 0.0_double
                   v(k,atom) = dt * f(k,gatom) ! 2014/08/01 SA
                end do
             end do
             !dt = dt * fire_f_dec
             fire_alpha = fire_alpha0
             fire_N = 1
             fire_N2 = 0
          else
             fire_alpha = fire_alpha * fire_f_alpha
             dt = min(dt * fire_f_inc,fire_max_step)
          end if
       end if
    else
       ! Decrease timestep and reset alpha
       dt = dt * fire_f_dec ! 2015/01/29 SA
       do atom = 1, ni_in_cell
          gatom = id_glob(atom)
          do k = 1, 3
             ! Normally FIRE sets v to zero if P<0
             v(k,atom) = 0.0_double
             !v(k,atom) = dt * f(k,gatom) ! 2014/08/01 SA
          end do
       end do
       !dt = dt * fire_f_dec
       fire_alpha = fire_alpha0
       fire_N = 1
    end if

    ! Output to file
    if (myid == 0) then
       if(iprint_MD>1) write (io_lun,11) iter, fire_N, fire_N2, fire_P, dt, fire_alpha, fire_norm_v
       call write_fire(fire_N, fire_N2, fire_P, dt, fire_alpha)
    end if

    ! update positions 
    do atom = 1, ni_in_cell
       gatom = id_glob(atom)
       flagx = flag_move_atom(1,gatom)
       flagy = flag_move_atom(2,gatom)
       flagz = flag_move_atom(3,gatom)
       ! X
       if (flagx) then
          dx = dt * v(1,atom)
          x_atom_cell(atom) = x_atom_cell(atom) + dx
       endif
       ! Y
       if (flagy) then
          dy = dt * v(2,atom)
          y_atom_cell(atom) = y_atom_cell(atom) + dy
       endif
       ! Z
       if (flagz) then
          dz = dt * v(3,atom)
          z_atom_cell(atom) = z_atom_cell(atom) + dz
       endif
    end do

11  format(4x,'iter ',i3, 2x, 'N= ', i3, 2x, 'N2= ', i3, 2x, 'P= ', f12.8,&
         2x,'dt= ', f9.6, 2x, 'alpha= ' ,f9.6, 2x, '|v|/|F|= ',f9.6)
12  format(i4,i4,f22.15,f22.15,f22.15)
    return

  end subroutine fire_qMD

  
end module Integrators
!*****
