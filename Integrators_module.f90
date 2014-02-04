! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module Integrators_module
! ------------------------------------------------------------------------------
! Code area : 
! ------------------------------------------------------------------------------

!****h* Conquest/Integrators_module
! NAME
!   Integrators_module
! PURPOSE
!   Evolves particle positions and velocities in MD loop
! AUTHOR
!   M.Arita
! CREATION DATE
!   2014/02/03
! MODIFICATION HISTORY
! SOURCE
!
module Integrators

  use datatypes
  implicit none
  character(80),save,private :: RCSid = "$Id$"
  
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
                            atom_coord_diff,flag_move_atom
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

   ibeg_atom = 1
   do atom = 1, ni_in_cell
     speca = species(atom)
     massa = mass(speca)*fac
     gatom = id_glob(atom)
     !do k = 1, 3
     !  if (.NOT.flag_movable(ibeg_atom+k-1)) cycle
     !  atom_coord_diff(k,gatom) = dt*v(k,atom)
     !enddo
     flagx = flag_movable(ibeg_atom)
     flagy = flag_movable(ibeg_atom+1)
     flagz = flag_movable(ibeg_atom+2)
     ibeg_atom = ibeg_atom + 3
     ! X
     if (flagx) then
       atom_coord_diff(1,gatom) = dt*v(1,atom)
       x_atom_cell(atom) = x_atom_cell(atom) + atom_coord_diff(1,gatom)
     endif
     ! Y
     if (flagy) then
       atom_coord_diff(2,gatom) = dt*v(2,atom)
       y_atom_cell(atom) = y_atom_cell(atom) + atom_coord_diff(2,gatom)
     endif
     ! Z
     if (flagz) then
       atom_coord_diff(3,gatom) = dt*v(3,atom)
       z_atom_cell(atom) = z_atom_cell(atom) + atom_coord_diff(3,gatom)
     endif
   enddo

   return
  end subroutine vVerlet_r_dt
  !*****

  !****f* Integrators/vVerlet_r_dt
  ! PURPOSE
  !   Evolve particle velocities by dt/2 via velocity Verlet
  !    v(t+dt/2) = v(t) + (dt/2)*a(t)
  !   When quenched-MD applies, calculate inner product v*f
  !   If v*f < 0, set v = 0
  !
  !   Note that the initial velocity is defined as v(0), NOT v(-dt/2)
  !   any longer as in velocityVerlet at move_atoms.module
  ! USAGE
  !   call vVerlet_v_dthalf(dt,v,f,flag_movable,second_call)
  ! INPUTS
  !   double dt           : MD time step
  !   double v            : velocities at t+dt/2
  !   double force        : forces at t
  !   logical flag_movable: flag to tell if atoms are movable
  !   logical second_call : tell this is the 2nd call
  ! OUTPUT
  !   real(double), v: half-a-step evolved particle velocities
  ! AUTHOR
  !   M.Arita
  ! CREATION DATE
  !   2014/02/03
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine vVerlet_v_dthalf(dt,v,f,flag_movable,second_call)
   ! Module usage
   use numbers, ONLY: half,zero
   use global_module, ONLY: ni_in_cell,id_glob,flag_quench_MD
   use species_module, oNLY: species,mass
   use move_atoms, ONLY: fac

   implicit none
   ! passed variables
   real(double),intent(in) :: dt
   real(double),dimension(3,ni_in_cell),intent(in)    :: f
   real(double),dimension(3,ni_in_cell),intent(inout) :: v
   logical,dimension(3*ni_in_cell),intent(in)         :: flag_movable
   logical,optional :: second_call
   ! local variables
   integer :: atom,speca,gatom,k,ibeg_atom
   real(double) :: vf,massa
   logical :: flagx,flagy,flagz

   ibeg_atom=1
   ! for quenched-MD
   if (present(second_call) .AND. flag_quench_MD) then
     do atom = 1, ni_in_cell
       speca = species(atom)
       massa = mass(speca)*fac
       gatom = id_glob(atom)
       do k = 1, 3
         if (.NOT.flag_movable(ibeg_atom+k-1)) cycle
         vf = v(k,atom)+f(k,gatom)
         if (vf.LT.0) v(k,atom) = zero
         v(k,atom) = v(k,atom) + half*dt*f(k,gatom)/massa
       enddo
     enddo
   ! otherwise
   else
     do atom = 1, ni_in_cell
       speca = species(atom)
       massa = mass(speca)*fac
       gatom = id_glob(atom)
       do k = 1, 3
         if (.NOT.flag_movable(ibeg_atom+k-1)) cycle
         v(k,atom) = v(k,atom) + half*dt*f(k,gatom)/massa
       enddo
       ibeg_atom = ibeg_atom + 3
!      flagx = flag_move_atom(1,gatom)
!      flagy = flag_move_atom(2,gatom)
!      flagz = flag_move_atom(3,gatom)
!      ! X
!      if (flagx) v(1,atom) = v(1,atom) + dt*half*f(1,gatom)/massa
!      ! Y
!      if (flagy) v(2,atom) = v(2,atom) + dt*half*f(2,gatom)/massa
!      ! Z
!      if (flagz) v(3,atom) = v(3,atom) + dt*half*f(3,gatom)/massa
     enddo
   endif

   return
  end subroutine vVerlet_v_dthalf
  
end module Integrators
!*****
