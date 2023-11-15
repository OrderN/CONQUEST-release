! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module constraint_module
! ------------------------------------------------------------------------------
! Code area : 
! ------------------------------------------------------------------------------

!****h* Conquest/constraint_module
! NAME
!  constraint_module
! PURPOSE
!  Constrain atoms during MD simulations
! AUTHOR
!  M.Arita
! CREATION DATE
!  2014/02/04
! MODIFICATION HISTORY
! SOURCE
!
module constraint_module

  use datatypes
  use auxiliary_types
  implicit none

  ! public
  logical,public :: flag_RigidBonds
  integer,public :: maxiterSHAKE,maxiterRATTLE
  integer,allocatable,public :: n_bond(:)
  real(double),public :: SHAKE_tol,RATTLE_tol,const_range
  type(group_aux),target,public :: constraints
  ! private
  real(double),allocatable,private :: dij2(:),vec_rji(:)
  integer,allocatable,private :: ibeg_dij2(:),ibeg_vec_rji(:)

contains

  !****f* constraint_module/ready_constraint
  ! NAME
  !  ready_constraint
  ! PURPOSE
  !  Calculates the fixed bond lengths and interatomic
  !  distances beforehand
  ! USAGE
  !  call ready_constraint()
  ! INPUTS
  ! OUTPUT
  ! AUTHOR
  !  M.Arita
  ! CREATION DATE
  !  2014/02/04
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine ready_constraint()
   ! Module usage
   use datatypes
   use numbers, ONLY: zero
   use global_module, ONLY: atom_coord,rcellx,rcelly,rcellz
   use auxiliary_types, ONLY: group_aux
   use GenComms, ONLY: cq_abort

   implicit none
   ! local variables
   integer :: i,j,k,l,ll,atom,bsize,jsize,stat
   integer :: nbond,gatom,max_nbond
   real(double) :: cutoff2,rx,ry,rz,rx2,ry2,rz2
   !real(double) :: x(5),y(5),z(5) ! Suppose group composed of 5 atoms at most
   real(double),allocatable :: x(:),y(:),z(:)
   type(group_aux),pointer :: SHAKE_RATTLE

   !% NOTE: I (michi) first thought the arrays in this subroutine should be
   !%       stored by ionode and broadcasted at the end. But it may take
   !%       too much time for communication depending on the system.
   !%       So all processors do the following operations.

   SHAKE_RATTLE => constraints   
   ! Allocations
   max_nbond = maxval(SHAKE_RATTLE%n_atom_in_grp)
   allocate (x(max_nbond),y(max_nbond),z(max_nbond), STAT=stat)
   if (stat.NE.0) call cq_abort('Error allocating x,y,z at ready_constraint: ', &
                                max_nbond)
   bsize = 0
   jsize = 0
   do i = 1, SHAKE_RATTLE%n_grp
     ! Count bonds
     if (SHAKE_RATTLE%grp_name(i)(1:5) .EQ. 'water') then
       n_bond(i) = 3
     else
       n_bond(i) = SHAKE_RATTLE%n_atom_in_grp(i) - 1
     endif
     bsize = bsize + n_bond(i)*SHAKE_RATTLE%n_subgrp(i)
     jsize = jsize + SHAKE_RATTLE%n_subgrp(i)
   enddo
   allocate (dij2(bsize), STAT=stat)
   if (stat.NE.0) call cq_abort('Error allocating dij2: ',bsize)
   allocate (ibeg_vec_rji(bsize), STAT=stat)
   if (stat.NE.0) call cq_abort('Error allocating ibeg_vec_rji: ',bsize)
   allocate (vec_rji(3*bsize), STAT=stat)
   if (stat.NE.0) call cq_abort('Error allocating vec_rji: ',3*bsize)
   allocate (ibeg_dij2(jsize), STAT=stat)
   if (stat.NE.0) call cq_abort('Error allocating ibeg_dji2: ',jsize)
   ! Reckon addresses
   ibeg_dij2(1) = 1
   l = 0
   nbond = n_bond(1)
   do i = 1, SHAKE_RATTLE%n_grp
     do j = 1, SHAKE_RATTLE%n_subgrp(i)
       l = l + 1
       if (l.NE.1) then
         ibeg_dij2(l) = ibeg_dij2(l-1) + nbond
       else
         cycle
       endif
       nbond = n_bond(i)
     enddo ! j
   enddo ! i
   ll = 0
   ibeg_vec_rji(1) = 1
   do i = 1, SHAKE_RATTLE%n_grp
     do j = 1, SHAKE_RATTLE%n_subgrp(i)
       do k = 1, n_bond(i)
         ll = ll + 1
         if (ll.NE.1) then
           ibeg_vec_rji(ll) = ibeg_vec_rji(ll-1) + 3
         else
           cycle
         endif
       enddo ! k
     enddo ! j
   enddo ! i
   cutoff2 = const_range*const_range ! in bohr
   l = 0
   ll = 0
   do i = 1, SHAKE_RATTLE%n_grp
     do j = 1, SHAKE_RATTLE%n_subgrp(i)
       l = l + 1
       do atom = 1, SHAKE_RATTLE%n_atom_in_grp(i)
         gatom = SHAKE_RATTLE%glob_atom(SHAKE_RATTLE%iatom_beg(l)+atom-1)
         x(atom) = atom_coord(1,gatom)
         y(atom) = atom_coord(2,gatom)
         z(atom) = atom_coord(3,gatom)
       enddo ! atom
       ! Calculate interatomic distances and wrap atoms back in the cell if necessary.
       ! The following do loop applies except group 'water***'
       ! NOTE: In most cases, n_bond = n_atom_in_grp - 1
       do atom = 1, SHAKE_RATTLE%n_atom_in_grp(i)-1
         ll = ll + 1
         ! x(1),y(1),z(1) must always be reference (mutually-associated atom)
         rx = x(atom+1)-x(1)
         ry = y(atom+1)-y(1)
         rz = z(atom+1)-z(1)
         rx2 = rx*rx
         ry2 = ry*ry
         rz2 = rz*rz
         ! X
         if (rx2 .GT. cutoff2) then
           if (rx.GT.zero) then
             x(atom+1) = x(atom+1) - rcellx
           elseif (rx.LT.zero) then
             x(atom+1) = x(atom+1) + rcellx
           endif
         endif
         ! Y
         if (ry2 .GT. cutoff2) then
           if (ry.GT.zero) then
             y(atom+1) = y(atom+1) - rcelly
           elseif (ry.LT.zero) then
             y(atom+1) = y(atom+1) + rcelly
           endif
         endif
         ! Z
         if (rz2 .GT. cutoff2) then
           if (rz.GT.zero) then
             z(atom+1) = z(atom+1) - rcellz
           elseif (rz.LT.zero) then
             z(atom+1) = z(atom+1) + rcellz
           endif
         endif
         ! Calculate and store fixed lengths dij2 and interatomic distances vec_rji
         ! NOTE: Cannot reuse rx2,ry2 and rz2 b/c some of x,y,z may be wrapped back
         dij2(ibeg_dij2(l)+atom-1) = (x(atom+1)-x(1))*(x(atom+1)-x(1)) + &
                                     (y(atom+1)-y(1))*(y(atom+1)-y(1)) + &
                                     (z(atom+1)-z(1))*(z(atom+1)-z(1))
         vec_rji(ibeg_vec_rji(ll)  ) = x(atom+1) - x(1)
         vec_rji(ibeg_vec_rji(ll)+1) = y(atom+1) - y(1)
         vec_rji(ibeg_vec_rji(ll)+2) = z(atom+1) - z(1)
       enddo ! atom
       ! Calculate H-H length for group 'water***'
       if (SHAKE_RATTLE%grp_name(i)(1:5) .EQ. 'water') then
         ll = ll + 1
         dij2(ibeg_dij2(l)+2) = (x(3)-x(2))*(x(3)-x(2)) + &
                                (y(3)-y(2))*(y(3)-y(2)) + &
                                (z(3)-z(2))*(z(3)-z(2))
         vec_rji(ibeg_vec_rji(ll)  ) = x(3) - x(2)
         vec_rji(ibeg_vec_rji(ll)+1) = y(3) - y(2)
         vec_rji(ibeg_vec_rji(ll)+2) = z(3) - z(2)
         !vec_rji
       endif
     enddo ! j
   enddo ! i
   ! Deallocations
   deallocate (x,y,z, STAT=stat)
   if (stat.NE.0) call cq_abort('Error deallocating x,y,z at ready_constraint')

   return
  end subroutine ready_constraint
  
  !****f* constraint_module/correct_atomic_position
  ! NAME
  !  correct_atomic_position
  ! PURPOSE
  !  Selector subroutine for constraining atomic positions
  ! USAGE
  !  call correct_atomic_position(v,dt)
  ! INPUTS
  !  real(double), v : atomic velocities
  !  real(double), dt: time step
  ! OUTPUT
  !  real(double), v: atomic velocities
  ! AUTHOR
  !  M.Arita
  ! CREATION DATE
  !  2014/02/04
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine correct_atomic_position(v,dt)
    ! Module usage
    use datatypes
    use global_module, ONLY: ni_in_cell,io_lun
    use GenComms, ONLY: myid

    implicit none
    ! passed variables
    real(double),intent(inout) :: v(3,ni_in_cell)
    real(double),intent(in)    :: dt
    ! local variables
    integer :: i,l,ll,lll

    l = 0
    ll = 0
    lll = 0
    do i = 1, constraints%n_grp
      if (constraints%grp_name(i)(1:5) .NE. 'water') then
        call Do_nSHAKE(i,l,ll,lll,v,dt)
      else
!%%!    if (flag_SETTLE) then
!%%!      call Do_rSETTLE()
!%%!    else
          call Do_waterSHAKE(i,l,ll,lll,v,dt)
!%%!    endif
      endif
    enddo ! i

    !if (myid.EQ.0) write (io_lun,*) "Finished correcting positions"

    return
  end subroutine correct_atomic_position

  !****f* constraint_module/Do_nSHAKE
  ! NAME
  !  Do_nSHAKE
  ! PURPOSE
  !  Correct atomic positions by the SHAKE algorithm
  ! USAGE
  !  call Do_nSHAKE(grpid,l,ll,lll,v,dt)
  ! INPUTS
  !  integer,      grp: group id
  !  integer,      l  : subgroup id
  !  integer,      ll : bond id for ibeg_vec_rji
  !  integer,      lll: bond id for ibeg_vec_rji (storing)
  !  real(double), v  : atomic velocities
  !  real(double), dt : time step
  ! OUTPUTS
  !  integer,      l  : subgroup id
  !  integer,      ll : bond id for ibeg_vec_rji
  !  integer,      lll: bond id for ibeg_vec_rji (stroing)
  !  real(double), v  : atomic velocities
  ! AUTHOR
  !  M.Arita
  ! CREATION DATE
  !  2014/02/04
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine Do_nSHAKE(grpid,l,ll,lll,v,dt)
    ! Module usage
    use datatypes
    use numbers, ONLY: zero,half,one
    use global_module, ONLY: ni_in_cell,id_glob_inv,x_atom_cell,y_atom_cell,z_atom_cell, &
                             rcellx,rcelly,rcellz,flag_LmatrixReuse,atom_coord_diff    , &
                             id_glob,iprint_MD,io_lun
    use auxiliary_types, ONLY: group_aux
    use GenComms, ONLY: cq_abort,myid
    use species_module, ONLY: species,mass

    implicit none
    ! passed variables
    integer,intent(in)    :: grpid
    integer,intent(inout) :: l,ll,lll
    real(double),intent(inout) :: v(3,ni_in_cell)
    real(double),intent(in)    :: dt
    ! local variables
    integer :: j,iter,nsize,stat,gatom_i,patom_i,spec_i,gatom_j,spec_j,bond,k
    integer,allocatable :: patom_j(:)
    real(double) :: mass_i,inv_mass_i,x_ji_old,y_ji_old,z_ji_old,cutoff2,inv_dt,  &
                    numerator,denominator,gamma,xac_lc_i,yac_lc_i,zac_lc_i,rij2,  &
                    residual,xi_old,yi_old,zi_old,xj_old,yj_old,zj_old
    real(double),dimension(3) :: delta,delta_i,diff_i,v_lc_i
    real(double),allocatable :: mass_j(:),inv_mass_j(:)
    real(double),allocatable :: xi_new(:),yi_new(:),zi_new(:), &
                                xj_new(:),yj_new(:),zj_new(:), &
                                delta_j(:,:),diff_j(:,:),v_lc_j(:,:)
    real(double),allocatable :: xac_lc_j(:),yac_lc_j(:),zac_lc_j(:)
    logical :: done,fail
    type(group_aux),pointer :: SHAKE

    SHAKE => constraints
    cutoff2 = const_range*const_range ! in bohr
    inv_dt = one/dt

    ! Allocations
    nsize  = n_bond(grpid)
    allocate (patom_j(nsize),mass_j(nsize),inv_mass_j(nsize), &
              xac_lc_j(nsize),yac_lc_j(nsize),zac_lc_j(nsize), STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating arrays for j: ', nsize)
    allocate (xi_new(nsize),yi_new(nsize),zi_new(nsize), &
              xj_new(nsize),yj_new(nsize),zj_new(nsize), &
              STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating arrays for new positions: ',nsize)
    allocate (delta_j(3,nsize),diff_j(3,nsize),v_lc_j(3,nsize), STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating arrays for j: ',3,nsize)

    do j = 1, SHAKE%n_subgrp(grpid)
      l = l + 1
      iter = 0
      done = .false.
      ! Stores information of reference "i"
      gatom_i = SHAKE%glob_atom(SHAKE%iatom_beg(l))
      patom_i = id_glob_inv(gatom_i)
      spec_i  = species(patom_i)
      mass_i  = mass(spec_i)
      inv_mass_i = one/mass_i
      xac_lc_i = x_atom_cell(patom_i)
      yac_lc_i = y_atom_cell(patom_i)
      zac_lc_i = z_atom_cell(patom_i)
      if (flag_LmatrixReuse) &
        diff_i(1:3) = atom_coord_diff(1:3,gatom_i)
      v_lc_i(1:3) = v(1:3,patom_i)
      ! Stores information of "j"
      do bond = 1, nsize
        gatom_j = SHAKE%glob_atom(SHAKE%iatom_beg(l)+bond)
        patom_j(bond) = id_glob_inv(gatom_j)
        spec_j  = species(patom_j(bond))
        mass_j(bond)     = mass(spec_j)
        inv_mass_j(bond) = one/mass_j(bond)
        xac_lc_j(bond) = x_atom_cell(patom_j(bond))
        yac_lc_j(bond) = y_atom_cell(patom_j(bond))
        zac_lc_j(bond) = z_atom_cell(patom_j(bond))
        if (flag_LmatrixReuse) &
          diff_j(1:3,bond) = atom_coord_diff(1:3,gatom_j)
        v_lc_j(1:3,bond) = v(1:3,patom_j(bond))
      enddo ! bond
      !-----------!
      ! MAIN LOOP !
      !-----------!
      do while (.NOT.done)
        iter = iter + 1
        fail = .false.
        do bond = 1, nsize
          ll = ll + 1
          xi_old = xac_lc_i
          yi_old = yac_lc_i
          zi_old = zac_lc_i
          xj_old = xac_lc_j(bond)
          yj_old = yac_lc_j(bond)
          zj_old = zac_lc_j(bond)
          ! r^{old}_{j}(t+dt) - r^{old}_{i}(t+dt)
          x_ji_old = xj_old - xi_old
          y_ji_old = yj_old - yi_old
          z_ji_old = zj_old - zi_old
          ! Wrap "j" in/out if necessary
          ! X
          if (x_ji_old*x_ji_old .GT. cutoff2) then
            if (x_ji_old.GT.zero) then
              xj_old = xac_lc_j(bond) - rcellx
            elseif (x_ji_old.LT.zero) then
              xj_old = xac_lc_j(bond) + rcellx
            endif
            x_ji_old = xj_old - xi_old
          endif
          ! Y
          if (y_ji_old*y_ji_old .GT. cutoff2) then
            if (y_ji_old.GT.zero) then
              yj_old = yac_lc_j(bond) - rcelly
            elseif (y_ji_old.LT.zero) then
              yj_old = yac_lc_j(bond) + rcelly
            endif
            y_ji_old = yj_old - yi_old
          endif
          ! Z
          if (z_ji_old*z_ji_old .GT. cutoff2) then
            if (z_ji_old.GT.zero) then
              zj_old = zac_lc_j(bond) - rcellz
            elseif (z_ji_old.LT.zero) then
              zj_old = zac_lc_j(bond) + rcellz
            endif
            z_ji_old = zj_old - zi_old
          endif
          ! Calculates gamma
          numerator = x_ji_old*x_ji_old + y_ji_old*y_ji_old + z_ji_old*z_ji_old - &
                      dij2(ibeg_dij2(l)+bond-1)
          numerator = mass_i*mass_j(bond)*numerator          
          denominator = vec_rji(ibeg_vec_rji(ll)  )*x_ji_old + &
                        vec_rji(ibeg_vec_rji(ll)+1)*y_ji_old + &
                        vec_rji(ibeg_vec_rji(ll)+2)*z_ji_old
          denominator = (mass_i+mass_j(bond)) * denominator
          gamma = half*numerator/denominator
          ! Corrects atomic positions
          do k = 1, 3
            delta(k)        = - gamma*vec_rji(ibeg_vec_rji(ll)+k-1)
            delta_i(k)      = - inv_mass_i*delta(k)
            delta_j(k,bond) = inv_mass_j(bond)*delta(k)
          enddo
          ! "i" (reference)
          xi_new(bond) = xi_old   + delta_i(1)
          yi_new(bond) = yi_old   + delta_i(2)
          zi_new(bond) = zi_old   + delta_i(3)
          xac_lc_i     = xac_lc_i + delta_i(1)
          yac_lc_i     = yac_lc_i + delta_i(2)
          zac_lc_i     = zac_lc_i + delta_i(3)
          ! "j"
          xj_new(bond)   = xj_old         + delta_j(1,bond)
          yj_new(bond)   = yj_old         + delta_j(2,bond)
          zj_new(bond)   = zj_old         + delta_j(3,bond)
          xac_lc_j(bond) = xac_lc_j(bond) + delta_j(1,bond)
          yac_lc_j(bond) = yac_lc_j(bond) + delta_j(2,bond)
          zac_lc_j(bond) = zac_lc_j(bond) + delta_j(3,bond)
          ! Corrects atomic displacements
          if (flag_LmatrixReuse) then
            diff_i(1:3) = diff_i(1:3) + delta_i(1:3)                ! "i" (reference)
            diff_j(1:3,bond) = diff_j(1:3,bond) + delta_j(1:3,bond) ! "j"
          endif
          ! Corrects velocities
          v_lc_i(1:3) = v_lc_i(1:3) + delta_i(1:3)*inv_dt                ! "i" (reference)
          v_lc_j(1:3,bond) = v_lc_j(1:3,bond) + delta_j(1:3,bond)*inv_dt ! "j"
        enddo ! bond
        ! Convergence check
        do bond = 1, nsize
          rij2 = (xj_new(bond)-xi_new(bond))*(xj_new(bond)-xi_new(bond)) + &
                 (yj_new(bond)-yi_new(bond))*(yj_new(bond)-yi_new(bond)) + &
                 (zj_new(bond)-zi_new(bond))*(zj_new(bond)-zi_new(bond))
          residual = rij2 - dij2(ibeg_dij2(l)+bond-1)
          residual = abs(residual)
          if (residual .GT. SHAKE_tol) then
            ll = ll - nsize
            fail = .true.
            exit
          endif
        enddo
        if (.NOT.fail) done = .true.
        if (iter.GT.maxiterSHAKE) &
          call cq_abort('Error: Too many SHAKE iterations !')
      enddo ! END of MAIN LOOP
     ! Stores the corrected interatomic distances, positions, displacements
     ! and velocities
     do bond = 1, nsize
       lll = lll + 1
       vec_rji(ibeg_vec_rji(lll)  ) = xj_new(bond) - xi_new(bond)
       vec_rji(ibeg_vec_rji(lll)+1) = yj_new(bond) - yi_new(bond)
       vec_rji(ibeg_vec_rji(lll)+2) = zj_new(bond) - zi_new(bond)
       ! Puts back for "j"
       x_atom_cell(patom_j(bond)) = xac_lc_j(bond)
       y_atom_cell(patom_j(bond)) = yac_lc_j(bond)
       z_atom_cell(patom_j(bond)) = zac_lc_j(bond)
       if (flag_LmatrixReuse) &
         atom_coord_diff(1:3,id_glob(patom_j(bond))) = diff_j(1:3,bond)
       v(1:3,patom_j(bond)) = v_lc_j(1:3,bond)
     enddo
     ! Puts back for "i" (reference)
     x_atom_cell(patom_i) = xac_lc_i
     y_atom_cell(patom_i) = yac_lc_i
     z_atom_cell(patom_i) = zac_lc_i
     if (flag_LmatrixReuse) &
       atom_coord_diff(1:3,gatom_i) = diff_i(1:3)
     v(1:3,patom_i) = v_lc_i(1:3)
     if (iprint_MD.GE.4 .AND. myid.EQ.0) &
       write (io_lun,fmt='(10x,a,4i5)') "grp, subgrp, SHAKE iter, residual: ", grpid,j,iter,residual
    enddo ! j

    ! Deallocations
    deallocate (v_lc_j,diff_j,delta_j, STAT=stat)
    if (stat.NE.0) call cq_abort('Error deallocating arrays for j: ',3,nsize)
    deallocate (zj_new,yj_new,xj_new,zi_new,yi_new,xi_new, STAT=stat)
    if (stat.NE.0) call cq_abort('Error deallocating arrays for new positions: ',nsize)
    deallocate (zac_lc_j,yac_lc_j,xac_lc_j,inv_mass_j,mass_j,patom_j, STAT=stat)
    if (stat.NE.0) call cq_abort('Error deallocating arrays for j: ',nsize)

    return
  end subroutine Do_nSHAKE

  !****f* constraint_module/Do_waterSHAKE
  ! NAME
  !  Do_waterSHAKE
  ! PURPOSE
  !  Correct water geometries by the SHAKE algorithm
  ! USAGE
  !  call Do_waterSHAKE(grpid,l,ll,lll,v,dt)
  ! INPUTS
  !  integer,      grp: group id
  !  integer,      l  : subgroup id
  !  integer,      ll : bond id for ibeg_vec_rji
  !  integer,      lll: bond id for ibeg_vec_rji (storing)
  !  real(double), v  : atomic velocities
  !  real(double), dt : time step
  ! OUTPUTS
  !  integer,      l  : subgroup id
  !  integer,      ll : bond id for ibeg_vec_rji
  !  integer,      lll: bond id for ibeg_vec_rji (stroing)
  !  real(double), v  : atomic velocities
  ! AUTHOR
  !  M.Arita
  ! CREATION DATE
  !  2014/01/21
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine Do_waterSHAKE(grpid,l,ll,lll,v,dt)
    ! Module usage
    use datatypes
    use numbers, ONLY: zero,half,one
    use global_module, ONLY: ni_in_cell,id_glob_inv,x_atom_cell,y_atom_cell,z_atom_cell   , &
                             flag_LmatrixReuse,atom_coord_diff,io_lun,rcellx,rcelly,rcellz, &
                             id_glob,iprint_MD
    use auxiliary_types, ONLY: group_aux
    use species_module, ONLY: species,mass
    use GenComms, ONLY: cq_abort,myid

    implicit none
    ! passed variables
    integer,intent(in) :: grpid
    integer,intent(inout) :: l,ll,lll
    real(double),intent(inout) :: v(3,ni_in_cell)
    real(double),intent(in) :: dt
    ! local variables
    integer :: j,iter,gatom_i,patom_i,spec_i,OHbond,gatom_j,spec_j,k,bond
    integer,dimension(2) :: patom_j
    real(double) :: cutoff2,inv_dt,mass_i,inv_mass_i,xac_lc_i,yac_lc_i,zac_lc_i, &
                    xi_old,yi_old,zi_old,xj_old,yj_old,zj_old,x_ji_old,y_ji_old, &
                    z_ji_old,numerator,denominator,gamma,rij2,residual
    real(double),dimension(2) :: mass_j,inv_mass_j,xac_lc_j,yac_lc_j,zac_lc_j
    real(double),dimension(3) :: diff_i,v_lc_i,delta_i,delta,xi_new,yi_new,zi_new, &
                                 xj_new,yj_new,zj_new
    real(double),dimension(3,2) :: diff_j,v_lc_j,delta_j
    logical :: done,fail
    type(group_aux),pointer :: SHAKE

    SHAKE => constraints
    cutoff2 = const_range*const_range
    inv_dt = one/dt

    do j = 1, SHAKE%n_subgrp(grpid)
      l = l + 1
      iter = 0
      done = .false.
      ! Stores information of reference "i"
      gatom_i = SHAKE%glob_atom(SHAKE%iatom_beg(l))
      patom_i = id_glob_inv(gatom_i)
      spec_i  = species(patom_i)
      mass_i  = mass(spec_i)
      inv_mass_i = one/mass_i
      xac_lc_i = x_atom_cell(patom_i)
      yac_lc_i = y_atom_cell(patom_i)
      zac_lc_i = z_atom_cell(patom_i)
      if (flag_LmatrixReuse) &
        diff_i(1:3) = atom_coord_diff(1:3,gatom_i)
      v_lc_i(1:3) = v(1:3,patom_i)
      ! Stores information of "j" or water hydrogen
      do OHbond = 1, 2
        gatom_j = SHAKE%glob_atom(SHAKE%iatom_beg(l)+OHbond)
        patom_j(OHbond) = id_glob_inv(gatom_j)
        spec_j  = species(patom_j(OHbond))
        mass_j(OHbond) = mass(spec_j)
        inv_mass_j(OHbond) = one/mass_j(OHbond)
        xac_lc_j(OHbond) = x_atom_cell(patom_j(OHbond))
        yac_lc_j(OHbond) = y_atom_cell(patom_j(OHbond))
        zac_lc_j(OHbond) = z_atom_cell(patom_j(OHbond))
        if (flag_LmatrixReuse) &
          diff_j(1:3,OHbond) = atom_coord_diff(1:3,gatom_j)
        v_lc_j(1:3,OHbond) = v(1:3,patom_j(OHbond))
      enddo ! OHbond
      !-----------!
      ! MAIN LOOP !
      !-----------!
      do while (.NOT.done)
        iter = iter + 1
        fail = .false.
        do OHbond = 1, 2
          ll = ll + 1
          xi_old = xac_lc_i
          yi_old = yac_lc_i
          zi_old = zac_lc_i
          xj_old = xac_lc_j(OHbond)
          yj_old = yac_lc_j(OHbond)
          zj_old = zac_lc_j(OHbond)
           ! r^{old}_{j}(t+dt) - r^{old}_{i}(t+dt)
          x_ji_old = xj_old - xi_old
          y_ji_old = yj_old - yi_old
          z_ji_old = zj_old - zi_old
          ! Wrap "j" in/out if necessary
          ! X
          if (x_ji_old*x_ji_old .GT. cutoff2) then
            if (x_ji_old.GT.zero) then
              xj_old = xac_lc_j(OHbond) - rcellx
            elseif (x_ji_old.LT.zero) then
              xj_old = xac_lc_j(OHbond) + rcellx
            endif
            x_ji_old = xj_old - xi_old
          endif
          ! Y
          if (y_ji_old*y_ji_old .GT. cutoff2) then
            if (y_ji_old.GT.zero) then
              yj_old = yac_lc_j(OHbond) - rcelly
            elseif (y_ji_old.LT.zero) then
              yj_old = yac_lc_j(OHbond) + rcelly
            endif
            y_ji_old = yj_old - yi_old
          endif
          ! Z
          if (z_ji_old*z_ji_old .GT. cutoff2) then
            if (z_ji_old.GT.zero) then
              zj_old = zac_lc_j(OHbond) - rcellz
            elseif (z_ji_old.LT.zero) then
              zj_old = zac_lc_j(OHbond) + rcellz
            endif
            z_ji_old = zj_old - zi_old
          endif
          ! Calculates gamma
          numerator = x_ji_old*x_ji_old + y_ji_old*y_ji_old + z_ji_old*z_ji_old - &
                      dij2(ibeg_dij2(l)+OHbond-1)
          numerator = mass_i*mass_j(OHbond)*numerator          
          denominator = vec_rji(ibeg_vec_rji(ll)  )*x_ji_old + &
                        vec_rji(ibeg_vec_rji(ll)+1)*y_ji_old + &
                        vec_rji(ibeg_vec_rji(ll)+2)*z_ji_old
          denominator = (mass_i+mass_j(OHbond)) * denominator
          gamma = half*numerator/denominator
          ! Corrects atomic positions
          do k = 1, 3
            delta(k)        = - gamma*vec_rji(ibeg_vec_rji(ll)+k-1)
            delta_i(k)      = - inv_mass_i*delta(k)
            delta_j(k,OHbond) = inv_mass_j(OHbond)*delta(k)
          enddo
          ! "i" (reference) or water oxygen
          xi_new(OHbond) = xi_old   + delta_i(1)
          yi_new(OHbond) = yi_old   + delta_i(2)
          zi_new(OHbond) = zi_old   + delta_i(3)
          xac_lc_i       = xac_lc_i + delta_i(1)
          yac_lc_i       = yac_lc_i + delta_i(2)
          zac_lc_i       = zac_lc_i + delta_i(3)
          ! "j" or water hydrogen
          xj_new(OHbond)   = xj_old           + delta_j(1,OHbond)
          yj_new(OHbond)   = yj_old           + delta_j(2,OHbond)
          zj_new(OHbond)   = zj_old           + delta_j(3,OHbond)
          xac_lc_j(OHbond) = xac_lc_j(OHbond) + delta_j(1,OHbond)
          yac_lc_j(OHbond) = yac_lc_j(OHbond) + delta_j(2,OHbond)
          zac_lc_j(OHbond) = zac_lc_j(OHbond) + delta_j(3,OHbond)
          ! Corrects atomic displacements
          if (flag_LmatrixReuse) then
            diff_i(1:3) = diff_i(1:3) + delta_i(1:3)                      ! "i" (reference)
            diff_j(1:3,OHbond) = diff_j(1:3,OHbond) + delta_j(1:3,OHbond) ! "j"
          endif
          ! Corrects velocities
          v_lc_i(1:3) = v_lc_i(1:3) + delta_i(1:3)*inv_dt                      ! "i" (reference)
          v_lc_j(1:3,OHbond) = v_lc_j(1:3,OHbond) + delta_j(1:3,OHbond)*inv_dt ! "j"
        enddo ! OHbond
        ! Now constrains H-H length. Think of H1 as a reference.
        ll = ll + 1
        xi_old = xac_lc_j(1) ! x of H1
        yi_old = yac_lc_j(1) ! y of H1
        zi_old = zac_lc_j(1) ! z of H1
        xj_old = xac_lc_j(2) ! x of H2
        yj_old = yac_lc_j(2) ! y of H2
        zj_old = zac_lc_j(2) ! z of H2
        ! r^{old}_{H2}(t+dt) - r^{old}_{H1}(t+dt)
        x_ji_old = xj_old - xi_old
        y_ji_old = yj_old - yi_old
        z_ji_old = zj_old - zi_old
        ! Wrap "H2" in/out if necessary
        ! X
        if (x_ji_old*x_ji_old .GT. cutoff2) then
          if (x_ji_old.GT.zero) then
            xj_old = xac_lc_j(2) - rcellx
          elseif (x_ji_old.LT.zero) then
            xj_old = xac_lc_j(2) + rcellx
          endif
          x_ji_old = xj_old - xi_old
        endif
        ! Y
        if (y_ji_old*y_ji_old .GT. cutoff2) then
          if (y_ji_old.GT.zero) then
            yj_old = yac_lc_j(2) - rcelly
          elseif (y_ji_old.LT.zero) then
            yj_old = yac_lc_j(2) + rcelly
          endif
          y_ji_old = yj_old - yi_old
        endif
        ! Z
        if (z_ji_old*z_ji_old .GT. cutoff2) then
          if (z_ji_old.GT.zero) then
            zj_old = zac_lc_j(2) - rcellz
          elseif (y_ji_old.LT.zero) then
            zj_old = zac_lc_j(2) + rcellz
          endif
          z_ji_old = zj_old - zi_old
        endif
        ! Calculates gamma
        numerator = x_ji_old*x_ji_old + y_ji_old*y_ji_old + z_ji_old*z_ji_old - &
                    dij2(ibeg_dij2(l)+2)
        numerator = mass_j(1)*mass_j(2)*numerator          
        denominator = vec_rji(ibeg_vec_rji(ll)  )*x_ji_old + &
                      vec_rji(ibeg_vec_rji(ll)+1)*y_ji_old + &
                      vec_rji(ibeg_vec_rji(ll)+2)*z_ji_old
        denominator = (mass_j(1)+mass_j(2)) * denominator
        gamma = half*numerator/denominator
        ! Corrects atomic positions
        do k = 1, 3
          delta(k)     = - gamma*vec_rji(ibeg_vec_rji(ll)+k-1)
          delta_j(k,1) = - inv_mass_j(1)*delta(k)
          delta_j(k,2) =   inv_mass_j(2)*delta(k)
        enddo
        ! "i" (reference)
        xi_new(3)   = xi_old      + delta_j(1,1)
        yi_new(3)   = yi_old      + delta_j(2,1)
        zi_new(3)   = zi_old      + delta_j(3,1)
        xac_lc_j(1) = xac_lc_j(1) + delta_j(1,1)
        yac_lc_j(1) = yac_lc_j(1) + delta_j(2,1)
        zac_lc_j(1) = zac_lc_j(1) + delta_j(3,1)
        ! "j"
        xj_new(3)   = xj_old      + delta_j(1,2)
        yj_new(3)   = yj_old      + delta_j(2,2)
        zj_new(3)   = zj_old      + delta_j(3,2)
        xac_lc_j(2) = xac_lc_j(2) + delta_j(1,2)
        yac_lc_j(2) = yac_lc_j(2) + delta_j(2,2)
        zac_lc_j(2) = zac_lc_j(2) + delta_j(3,2)
        ! Corrects atomic displacements
        if (flag_LmatrixReuse) then
          diff_j(1:3,1) = diff_j(1:3,1) + delta_j(1:3,1) ! "H1" (reference)
          diff_j(1:3,2) = diff_j(1:3,2) + delta_j(1:3,2) ! "H2"
        endif
        ! Corrects velocities
        v_lc_j(1:3,1) = v_lc_j(1:3,1) + delta_j(1:3,1)*inv_dt ! "H1" (reference)
        v_lc_j(1:3,2) = v_lc_j(1:3,2) + delta_j(1:3,2)*inv_dt ! "H2"
        ! Convergence check
        do bond = 1, 3
          rij2 = (xj_new(bond)-xi_new(bond))*(xj_new(bond)-xi_new(bond)) + &
                 (yj_new(bond)-yi_new(bond))*(yj_new(bond)-yi_new(bond)) + &
                 (zj_new(bond)-zi_new(bond))*(zj_new(bond)-zi_new(bond))
          residual = rij2 - dij2(ibeg_dij2(l)+bond-1)
          residual = abs(residual)
          if (residual .GT. SHAKE_tol) then
            ll = ll - 3
            fail = .true.
            exit
          endif
        enddo ! bond
        if (.NOT.fail) done = .true.
        if (iter.GT.maxiterSHAKE) &
          call cq_abort('Error: Too many SHAKE iterations !')
      enddo ! END of MAIN LOOP
      ! Stores the corrected interatomic distances, positions, displacements
      ! and velocities
      do bond = 1, 3
        lll = lll + 1
        vec_rji(ibeg_vec_rji(lll)  ) = xj_new(bond) - xi_new(bond)
        vec_rji(ibeg_vec_rji(lll)+1) = yj_new(bond) - yi_new(bond)
        vec_rji(ibeg_vec_rji(lll)+2) = zj_new(bond) - zi_new(bond)
      enddo ! bond
      do OHbond = 1, 2
        ! Puts back for "j" (water hydrogen)
        x_atom_cell(patom_j(OHbond)) = xac_lc_j(OHbond)
        y_atom_cell(patom_j(OHbond)) = yac_lc_j(OHbond)
        z_atom_cell(patom_j(OHbond)) = zac_lc_j(OHbond)
        if (flag_LmatrixReuse) &
          atom_coord_diff(1:3,id_glob(patom_j(OHbond))) = diff_j(1:3,OHbond)
        v(1:3,patom_j(OHbond)) = v_lc_j(1:3,OHbond)
      enddo ! OHbond
      ! Puts back for "i" (water oxygen)
      x_atom_cell(patom_i) = xac_lc_i
      y_atom_cell(patom_i) = yac_lc_i
      z_atom_cell(patom_i) = zac_lc_i
      if (flag_LmatrixReuse) &
        atom_coord_diff(1:3,gatom_i) = diff_i(1:3)
      v(1:3,patom_i) = v_lc_i(1:3)
      if (iprint_MD.GE.4 .AND. myid.EQ.0) &
        write (io_lun,fmt='(10x,a,4i5)') "grp, subgrp, SHAKE iter, residual: ", grpid,j,iter,residual
    enddo ! j

    return
  end subroutine Do_waterSHAKE

  !****f* constraint_module/correct_atomic_velocity
  ! NAME
  !  correct_atomic_velocity
  ! PURPOSE
  !  Selector subroutine for constraining atomic velocities
  ! USAGE
  !  call correct_atomic_velocity(v,dt)
  ! INPUT
  !  real(double), v: atomic velocities
  ! OUTPUT
  !  real(double), v: atomic velocities
  ! AUTHOR
  !  M.Arita
  ! CREATION DATE
  !  2014/01/20
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine correct_atomic_velocity(v)
    ! Module usage
    use datatypes
    use global_module, ONLY: ni_in_cell,io_lun
    use GenComms, ONLY: myid

    implicit none
    ! passed variables
    real(double),intent(inout) :: v(3,ni_in_cell)
    ! local variables
    integer :: i,l,ll,lll

    l = 0
    ll = 0
    lll = 0
    do i = 1, constraints%n_grp
      if (constraints%grp_name(i)(1:5) .NE. 'water') then
        call Do_nRATTLE(i,l,ll,lll,v)
      else
!%%!    if (flag_SETTLE) then
!%%!      call Do_vSETTLE()
!%%!    else
          call Do_waterRATTLE(i,l,ll,lll,v)
!%%!    endif
      endif
    enddo ! i

    !if (myid.EQ.0) write (io_lun,*) "Finished correcting velocities"

    return
  end subroutine correct_atomic_velocity

  !****f* constraint_module/Do_nRATTLE
  ! NAME
  !  Do_nRATTLE
  ! PURPOSE
  !  Correct atomic velocities by the RATTLE algorithm
  ! USAGE
  !  call Do_nRATTLE(grpid,l,ll,lll,v)
  ! INPUTS
  !  integer,      grp: group id
  !  integer,      l  : subgroup id
  !  integer,      ll : bond id for ibeg_vec_rji
  !  integer,      lll: bond id for ibeg_vec_rji
  !  real(double), v  : atomic velocities
  ! OUTPUTS
  !  integer,      l  : subgroup id
  !  integer,      ll : bond id for ibeg_vec_rji
  !  integer,      lll: bond id for ibeg_vec_rji
  !  real(double), v  : atomic velocities
  ! AUTHOR
  !  M.Arita
  ! CREATION DATE
  !  2014/01/20
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine Do_nRATTLE(grpid,l,ll,lll,v)
    ! Module usage
    use datatypes
    use numbers, ONLY: one
    use global_module, ONLY: ni_in_cell,id_glob_inv,iprint_MD,io_lun
    use auxiliary_types, ONLY: group_aux
    use GenComms, ONLY: cq_abort,myid
    use species_module, ONLY: species,mass

    implicit none
    ! passed variables
    integer,intent(in) :: grpid
    integer,intent(inout) :: l,ll,lll
    real(double),intent(inout) :: v(3,ni_in_cell)
    ! local variables
    integer :: nsize,stat,j,iter,gatom_i,patom_i,spec_i,bond,gatom_j,spec_j,k
    integer,allocatable :: patom_j(:)
    real(double) :: mass_i,inv_mass_i,vxi_lc,vyi_lc,vzi_lc,vji_x_old,vji_y_old, &
                    vji_z_old,numerator,denominator,eta,residual
    real(double),dimension(3) :: delta(3),delta_i(3)
    real(double),allocatable :: mass_j(:),inv_mass_j(:)
    real(double),allocatable :: vxj_lc(:),vyj_lc(:),vzj_lc(:), &
                                vji_x_new(:),vji_y_new(:),vji_z_new(:)
    real(double),allocatable :: delta_j(:,:)
    logical :: done,fail
    type(group_aux),pointer :: RATTLE

    RATTLE => constraints

    ! Allocations
    nsize = n_bond(grpid)
    allocate (patom_j(nsize),mass_j(nsize),inv_mass_j(nsize), STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating arrays for j:', nsize)
    allocate (vxj_lc(nsize),vyj_lc(nsize),vzj_lc(nsize), STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating velocities of j: ',nsize)
    allocate (vji_x_new(nsize),vji_y_new(nsize),vji_z_new(nsize), STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating relative velocities: ',nsize)
    allocate (delta_j(3,nsize), STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating delta_j: ',3,nsize)

    do j = 1, RATTLE%n_subgrp(grpid)
      l = l + 1
      iter = 0
      done = .false.
      ! Stores information of reference "i"
      gatom_i = RATTLE%glob_atom(RATTLE%iatom_beg(l))
      patom_i = id_glob_inv(gatom_i)
      spec_i  = species(patom_i)
      mass_i  = mass(spec_i)
      inv_mass_i = one/mass_i
      vxi_lc = v(1,patom_i)
      vyi_lc = v(2,patom_i)
      vzi_lc = v(3,patom_i)
      ! Stores information of "j"
      do bond = 1, nsize
        gatom_j = RATTLE%glob_atom(RATTLE%iatom_beg(l)+bond)
        patom_j(bond) = id_glob_inv(gatom_j)
        spec_j  = species(patom_j(bond))
        mass_j(bond) = mass(spec_j)
        inv_mass_j(bond) = one/mass_j(bond)
        vxj_lc(bond) = v(1,patom_j(bond))
        vyj_lc(bond) = v(2,patom_j(bond))
        vzj_lc(bond) = v(3,patom_j(bond))
      enddo ! bond
      !-----------!
      ! MAIN LOOP !
      !-----------!
      do while (.NOT.done)
        iter = iter + 1
        fail = .false.
        do bond = 1, nsize
          ll = ll + 1
          ! old relative velocities
          vji_x_old = vxj_lc(bond) - vxi_lc
          vji_y_old = vyj_lc(bond) - vyi_lc
          vji_z_old = vzj_lc(bond) - vzi_lc
          ! Calculates eta
          numerator = vec_rji(ibeg_vec_rji(ll)  )*vji_x_old + &
                      vec_rji(ibeg_vec_rji(ll)+1)*vji_y_old + &
                      vec_rji(ibeg_vec_rji(ll)+2)*vji_z_old
          numerator = mass_i*mass_j(bond)*numerator
          denominator = (mass_i+mass_j(bond)) * dij2(ibeg_dij2(l)+bond-1)
          eta = numerator/denominator
          ! Corrects atomic velocities
          do k = 1, 3
            delta(k)        = - eta*vec_rji(ibeg_vec_rji(ll)+k-1)
            delta_i(k)      = - inv_mass_i*delta(k)
            delta_j(k,bond) =   inv_mass_j(bond)*delta(k)
          enddo
          ! "i" (reference)
          vxi_lc = vxi_lc + delta_i(1)
          vyi_lc = vyi_lc + delta_i(2)
          vzi_lc = vzi_lc + delta_i(3)
          ! "j"
          vxj_lc(bond) = vxj_lc(bond) + delta_j(1,bond)
          vyj_lc(bond) = vyj_lc(bond) + delta_j(2,bond)
          vzj_lc(bond) = vzj_lc(bond) + delta_j(3,bond)
          ! Stores new relative velocities
          vji_x_new(bond) = vxj_lc(bond) - vxi_lc
          vji_y_new(bond) = vyj_lc(bond) - vyi_lc
          vji_z_new(bond) = vzj_lc(bond) - vzi_lc
        enddo ! bond
        ! Convergence check
        do bond = 1, nsize
          lll = lll + 1
          residual = vec_rji(ibeg_vec_rji(lll)  )*vji_x_new(bond) + &
                     vec_rji(ibeg_vec_rji(lll)+1)*vji_y_new(bond) + &
                     vec_rji(ibeg_vec_rji(lll)+2)*vji_z_new(bond)
          residual = abs(residual)
          if (residual.GT.RATTLE_tol) then
            ll = ll - nsize
            lll = lll - bond
            fail = .true.
            exit
          endif
        enddo ! bond
        if (.NOT.fail) done = .true.
        if (iter.GT.maxiterRATTLE) &
          call cq_abort('Error: Too many RATTLE iterations !')
      enddo ! END of MAIN LOOP
      ! Stores the corrected velocities
      ! Puts back for "j"
      do bond = 1, nsize
        v(1,patom_j(bond)) = vxj_lc(bond)
        v(2,patom_j(bond)) = vyj_lc(bond)
        v(3,patom_j(bond)) = vzj_lc(bond)
      enddo ! bond
      ! Puts back for "i"
      v(1,patom_i) = vxi_lc
      v(2,patom_i) = vyi_lc
      v(3,patom_i) = vzi_lc

      if (iprint_MD.GE.4 .AND. myid.EQ.0) &
        write (io_lun,fmt='(10x,a,4i5)') "grp, subgrp, RATTLE iter, residual:",grpid,j,iter,residual
    enddo ! j

    ! Deallocations
    deallocate (delta_j, STAT=stat)
    if (stat.NE.0) call cq_abort('Error deallocating delta_j: ',nsize)
    deallocate (vji_z_new,vji_y_new,vji_x_new, STAT=stat)
    if (stat.NE.0) call cq_abort('Error deallocating relative velocities: ',nsize)
    deallocate (vzj_lc,vyj_lc,vxj_lc, STAT=stat)
    if (stat.NE.0) call cq_abort('Error deallocating velocities of j: ',nsize)
    deallocate (inv_mass_j,mass_j,patom_j, STAT=stat)
    if (stat.NE.0) call cq_abort('Error deallocating arrays for j', nsize)

    return
  end subroutine Do_nRATTLE

  !****f* constraint_module/Do_waterRATTLE
  ! NAME
  !  Do_waterRATTLE
  ! PURPOSE
  !  Correct atomic velocities by the RATTLE algorithm for
  !  rigid water models
  ! USAGE
  !  call Do_waterRATTLE(grpid,l,ll,lll,v)
  ! INPUTS
  !  integer,      grp: group id
  !  integer,      l  : subgroup id
  !  integer,      ll : bond id for ibeg_vec_rji
  !  integer,      lll: bond id for ibeg_vec_rji
  !  real(double), v  : atomic velocities
  ! OUTPUTS
  !  integer,      l  : subgroup id
  !  integer,      ll : bond id for ibeg_vec_rji
  !  integer,      lll: bond id for ibeg_vec_rji
  !  real(double), v  : atomic velocities
  ! AUTHOR
  !  M.Arita
  ! CREATION DATE
  !  2014/01/20
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine Do_waterRATTLE(grpid,l,ll,lll,v)
    ! Module usage
    use datatypes
    use numbers, ONLY: one
    use global_module, ONLY: ni_in_cell,id_glob_inv,iprint_MD,io_lun
    use auxiliary_types, ONLY: group_aux
    use species_module, ONLY: species,mass
    use GenComms, ONLY: cq_abort,myid

    implicit none
    ! passed variables
    integer,intent(in) :: grpid
    integer,intent(inout) :: l,ll,lll
    real(double),intent(inout) :: v(3,ni_in_cell)
    ! local variables
    integer :: j,iter,gatom_i,patom_i,spec_i,OHbond,gatom_j,spec_j,k,bond
    integer,dimension(2) :: patom_j
    real(double) :: mass_i,inv_mass_i,vxi_lc,vyi_lc,vzi_lc,vji_x_old,vji_y_old,vji_z_old, &
                    numerator,denominator,eta,vxi_new,vyi_new,vzi_new,residual
    real(double),dimension(2) :: mass_j,inv_mass_j,vxj_lc,vyj_lc,vzj_lc
    real(double),dimension(3) :: delta,delta_i,vji_x_new,vji_y_new,vji_z_new
    real(double),dimension(3,2) :: delta_j
    logical :: done,fail
    type(group_aux),pointer :: RATTLE
    
    RATTLE => constraints

    do j = 1, RATTLE%n_subgrp(grpid)
      l = l + 1
      iter = 0
      done = .false.
      ! Stores information of reference "i" or water oxygen
      gatom_i = RATTLE%glob_atom(RATTLE%iatom_beg(l))
      patom_i = id_glob_inv(gatom_i)
      spec_i  = species(patom_i)
      mass_i  = mass(spec_i)
      inv_mass_i = one/mass_i
      vxi_lc = v(1,patom_i)
      vyi_lc = v(2,patom_i)
      vzi_lc = v(3,patom_i)
      ! Stores information of "j"
      do OHbond = 1, 2
        gatom_j = RATTLE%glob_atom(RATTLE%iatom_beg(l)+OHbond)
        patom_j(OHbond) = id_glob_inv(gatom_j)
        spec_j  = species(patom_j(OHbond))
        mass_j(OHbond)  = mass(spec_j)
        inv_mass_j(OHbond) = one/mass_j(OHbond)
        vxj_lc(OHbond) = v(1,patom_j(OHbond))
        vyj_lc(OHbond) = v(2,patom_j(OHbond))
        vzj_lc(OHbond) = v(3,patom_j(OHbond))
      enddo ! OHbond
      !-----------!
      ! MAIN LOOP !
      !-----------!
      do while (.NOT.done)
        iter = iter + 1
        fail = .false.
        do OHbond = 1, 2
          ll = ll + 1
          ! old relative velocities
          vji_x_old = vxj_lc(OHbond) - vxi_lc
          vji_y_old = vyj_lc(OHbond) - vyi_lc
          vji_z_old = vzj_lc(OHbond) - vzi_lc
          ! Calculates eta
          numerator = vec_rji(ibeg_vec_rji(ll)  )*vji_x_old + &
                      vec_rji(ibeg_vec_rji(ll)+1)*vji_y_old + &
                      vec_rji(ibeg_vec_rji(ll)+2)*vji_z_old
          numerator = mass_i*mass_j(OHbond)*numerator
          denominator = (mass_i+mass_j(OHbond)) * dij2(ibeg_dij2(l)+OHbond-1)
          eta = numerator/denominator
          ! Corrects atomic velocities
          do k = 1, 3
            delta(k)          = - eta*vec_rji(ibeg_vec_rji(ll)+k-1)
            delta_i(k)        = - inv_mass_i*delta(k)
            delta_j(k,OHbond) =   inv_mass_j(OHbond)*delta(k)
          enddo
          ! "i" (reference or water oxygen)
          vxi_lc = vxi_lc + delta_i(1)
          vyi_lc = vyi_lc + delta_i(2)
          vzi_lc = vzi_lc + delta_i(3)
          ! "j"
          vxj_lc(OHbond) = vxj_lc(OHbond) + delta_j(1,OHbond)
          vyj_lc(OHbond) = vyj_lc(OHbond) + delta_j(2,OHbond)
          vzj_lc(OHbond) = vzj_lc(OHbond) + delta_j(3,OHbond)
          ! Stores new relative velocities
          vji_x_new(OHbond) = vxj_lc(OHbond) - vxi_lc
          vji_y_new(OHbond) = vyj_lc(OHbond) - vyi_lc
          vji_z_new(OHbond) = vzj_lc(OHbond) - vzi_lc
        enddo ! OHbond
        ! Now constrains H-H length. Think of H1 as a reference
        ll = ll + 1
        ! old relative velocities
        vji_x_old = vxj_lc(2) - vxj_lc(1)
        vji_y_old = vyj_lc(2) - vyj_lc(1)
        vji_z_old = vzj_lc(2) - vzj_lc(1)
        ! Calculates eta
        numerator = vec_rji(ibeg_vec_rji(ll)  )*vji_x_old + &
                    vec_rji(ibeg_vec_rji(ll)+1)*vji_y_old + &
                    vec_rji(ibeg_vec_rji(ll)+2)*vji_z_old
        numerator = mass_j(1)*mass_j(2)*numerator
        denominator = (mass_j(1)+mass_j(2)) * dij2(ibeg_dij2(l)+2)
        eta = numerator/denominator
        ! Corrects atomic velocities
        do k = 1, 3
          delta(k)     = - eta*vec_rji(ibeg_vec_rji(ll)+k-1)
          delta_j(k,1) = - inv_mass_j(1)*delta(k)
          delta_j(k,2) =   inv_mass_j(2)*delta(k)
        enddo
        ! "i" (reference or H1)
        vxj_lc(1) = vxj_lc(1) + delta_j(1,1)
        vyj_lc(1) = vyj_lc(1) + delta_j(2,1)
        vzj_lc(1) = vzj_lc(1) + delta_j(3,1)
        ! "j" or H2
        vxj_lc(2) = vxj_lc(2) + delta_j(1,2)
        vyj_lc(2) = vyj_lc(2) + delta_j(2,2)
        vzj_lc(2) = vzj_lc(2) + delta_j(3,2)
        ! Stores new relative velocities
        vji_x_new(3) = vxj_lc(2) - vxj_lc(1)
        vji_y_new(3) = vyj_lc(2) - vyj_lc(1)
        vji_z_new(3) = vzj_lc(2) - vzj_lc(1)
        ! Convergence check
        do bond = 1, 3
          lll = lll + 1
          residual = vec_rji(ibeg_vec_rji(lll)  )*vji_x_new(bond) + &
                     vec_rji(ibeg_vec_rji(lll)+1)*vji_y_new(bond) + &
                     vec_rji(ibeg_vec_rji(lll)+2)*vji_z_new(bond)
          residual = abs(residual)
          if (residual.GT.RATTLE_tol) then
            ll = ll - 3
            lll = lll - bond
            fail = .true.
            exit
          endif
        enddo ! bond
        if (.NOT.fail) done = .true.
        if (iter.GT.maxiterRATTLE) &
          call cq_abort('Error: Too many RATTLE iterations !')
      enddo ! END of MAIN LOOP
      ! Stores the corrected velocities
      ! Puts back for "j" or water hydrogen
      do OHbond = 1, 2
        v(1,patom_j(OHbond)) = vxj_lc(OHbond)
        v(2,patom_j(OHbond)) = vyj_lc(OHbond)
        v(3,patom_j(OHbond)) = vzj_lc(OHbond)
      enddo ! OHbond
      ! Puts back for "i" or water oxygen
      v(1,patom_i) = vxi_lc
      v(2,patom_i) = vyi_lc
      v(3,patom_i) = vzi_lc
      if (iprint_MD.GE.4 .AND. myid.EQ.0) &
        write (io_lun,fmt='(10x,a,4i5)') "grp, subgrp, RATTLE iter, residual:",grpid,j,iter,residual
    enddo ! j

    return
  end subroutine Do_waterRATTLE

end module constraint_module
!*****
