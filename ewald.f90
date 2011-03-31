! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module ewald
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/ewald_module *
!!  NAME
!!   ewald_module
!!  PURPOSE
!!   Collects variables and subroutines related to ewald calculation
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   14/05/01
!!  MODIFICATION HISTORY
!!   15/05/2001 dave
!!    Small changes to fix bugs
!!   08/06/2001 dave
!!    Added RCS Id and Log tags and changed structure_factor to use
!!    GenComms and gsumv
!!    Changed gsumv to gsum in accordance with GenComms
!!    Added my_barrier to use GenComms
!!   17/06/2002 dave
!!    Added RCS static object, improved headers
!!   14:29, 26/02/2003 drb 
!!    Rewrote to use x_atom_cell (etc) instead of rionx
!!   14:35, 26/02/2003 drb 
!!    Added save so that gamma is preserved
!!   13:36, 22/09/2003 drb 
!!    Major debugging: added accurate erfc function, fixed ordering problem with forces and atomic coordinates, 
!!    and tidied
!!   12:41, 2003/10/22 drb & mjg
!!    Checking in all Mike's new code which uses partitions instead of supercells for realspace sum and has a 
!!    well-defined (and understood) definition of gamma
!!   13:00, 2003/10/22 drb & mjg
!!    Made mike's variables dynamically allocated
!!   2005/10/06 09:55 dave
!!    Bug fix in old ewald routines - was checking rij2 not r2_min in set_ewald
!!   2005/10/06 10:15 dave
!!    Further bug fix to partition_distance: various variables were not recalculated if pd_init not true
!!   13:02, 11/10/2005 drb 
!!    Small changes: check for log(0) in erfc and wrap lines to less than 132
!!   2008/02/04 17:23 dave
!!    Changes for output to file not stdout
!!   2008/07/18 ast
!!    Added timers
!!  SOURCE
!!
module ewald_module

  use datatypes
  use global_module, ONLY: io_lun
  use GenComms, ONLY: cq_abort

  implicit none
  save

  integer :: n_g_vectors, n_supercells
  integer :: number_of_g_vectors, n_ewald_partition_neighbours
  real(double) :: gamma
  real(double) :: ewald_gamma, ewald_real_cutoff, ewald_recip_cutoff
  real(double) :: ewald_real_cell_volume
  real(double) :: ewald_gaussian_self_energy, ewald_net_charge_energy
  real(double), allocatable, dimension(:) :: gx, gy, gz, g_factor, structure2
  complex(double_cplx), allocatable, dimension(:) :: cstructure
  real(double), allocatable, dimension(:) :: supercell_vec_x, supercell_vec_y, supercell_vec_z
  integer, dimension(:,:), allocatable :: ewald_partition_neighbour_list
  real(double), dimension(:), allocatable :: ewald_g_vector_x, ewald_g_vector_y, ewald_g_vector_z, &
       &ewald_g_factor, struc_fac_r, struc_fac_i

  real(double), allocatable, dimension(:,:) :: ewald_force
  real(double) :: ewald_energy
  real(double) :: ewald_accuracy 


  logical :: flag_old_ewald

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id$"
  !!***

contains

  !!****f* ewald_module/ewald *
  !!
  !!  NAME 
  !!   ewald
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Performs ewald summation for the Coulomb energy and finds
  !!   forces due to this
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   E.H.Hernandez
  !!  CREATION DATE
  !!   20/02/96
  !!  MODIFICATION HISTORY
  !!   14/05/01 DRB
  !!    Added ROBODoc headers and integrated into ewald module
  !!   08/06/2001 dave
  !!    Changed dgsum to gsum
  !!    Added my_barrier to use GenComms
  !!   14:47, 14/02/2003 drb 
  !!    Changed to use x_atom_cell not rionx
  !!   2008/02/06 08:00 dave & veronika
  !!    Subtle bug fix: species_glob not species
  !!  SOURCE
  !!
  subroutine ewald( )

    use datatypes
    use numbers
    use units
    use atoms
    use dimens
    use species_module, ONLY: charge, species
    use GenComms, ONLY: gsum, my_barrier, inode, ionode
    use global_module, ONLY: atom_coord, ni_in_cell, iprint_gen, area_general, &
                             species_glob, IPRINT_TIME_THRES3
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module, ONLY: cq_timer,start_timer,stop_print_timer,WITH_LEVEL

    implicit none

    ! Local variables

    integer :: i, j, k, n, stat

    real(double) :: argument, cosine, ThisErfc, e_real, e_reciprocal, &
         e_constant, exponential, force, g_dot_r, &
         q_i, q_j, rij, rij2, rxij, ryij, rzij,&
         sine, sum, sum_x, sum_y, sum_z, total_charge
    real(double), allocatable, dimension(:,:) :: c1_force, c2_force
    type(cq_timer) :: tmr_l_tmp1

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    allocate(c1_force(3,ni_in_cell),c2_force(3,ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating c1_force in ewald: ",ni_in_cell,stat)
    call reg_alloc_mem(area_general,6*ni_in_cell,type_dbl)
    ! initialise the force array
    c1_force = zero
    c2_force = zero
    ewald_force = zero

    ! first we'll calculate the  real-space part
    sum = zero

    ! first we loop over pairs of atoms in the simulation cell, then
    ! we loop over pairs involving an atom in the central cell and
    ! an atom in a supercell. These cases are treated differently to 
    ! avoid the self-interaction of atom in the simulation cell
    do k=1, n_my_atoms
       i = index_my_atoms(k)
       q_i = charge(species_glob(i))
       sum_x = zero
       sum_y = zero
       sum_z = zero
       do j=1, ni_in_cell
          q_j = charge(species_glob(j))
          if ( j/=i ) then
             rxij = ( atom_coord(1,i) - atom_coord(1,j) )
             ryij = ( atom_coord(2,i) - atom_coord(2,j) )
             rzij = ( atom_coord(3,i) - atom_coord(3,j) )
             rij2 = rxij * rxij + ryij * ryij + rzij * rzij
             rij = sqrt( rij2 )
             argument = gamma * rij 
             ThisErfc = erfc( argument )
             exponential = exp(-gamma * gamma * rij2)
             sum = sum + q_i * q_j * ThisErfc / rij
             force = q_j * (ThisErfc / rij &
                  + two * gamma * exponential / sqrt_pi) / rij2
             ! Changed 14:30, 20/08/2003 drb from + two... to - two...
             sum_x = sum_x + force * rxij
             sum_y = sum_y + force * ryij
             sum_z = sum_z + force * rzij
          end if
       end do
       c1_force(1,i) = sum_x
       c1_force(2,i) = sum_y
       c1_force(3,i) = sum_z
    end do
    ! now neighbouring cells
    call my_barrier()
    do k=1, n_my_atoms
       i = index_my_atoms(k) 
       q_i = charge(species_glob(i))
       sum_x = zero
       sum_y = zero
       sum_z = zero
       do j=1, ni_in_cell
          q_j = charge(species_glob(j))
          do n=1, n_supercells
             rxij = ( atom_coord(1,i) - atom_coord(1,j)) + supercell_vec_x(n) * r_super_x
             ryij = ( atom_coord(2,i) - atom_coord(2,j)) + supercell_vec_y(n) * r_super_y
             rzij = ( atom_coord(3,i) - atom_coord(3,j)) + supercell_vec_z(n) * r_super_z
             rij2 = rxij * rxij + ryij * ryij + rzij * rzij
             rij = sqrt( rij2 )
             argument = gamma * rij 
             ThisErfc = erfc( argument )
             exponential = exp(-gamma * gamma * rij2)
             sum = sum + q_i * q_j * ThisErfc / rij
             force = q_j * (ThisErfc / rij&
                  + two * gamma * exponential / sqrt_pi) / rij2 
             ! Changed 14:30, 20/08/2003 drb from + two... to - two...
             sum_x = sum_x + force * rxij
             sum_y = sum_y + force * ryij
             sum_z = sum_z + force * rzij
          end do
       end do
       c1_force(1,i) = (c1_force(1,i) + sum_x) * q_i 
       c1_force(2,i) = (c1_force(2,i) + sum_y) * q_i 
       c1_force(3,i) = (c1_force(3,i) + sum_z) * q_i 
    end do
    e_real = sum / two
    call gsum(e_real)
    !    if(inode==ionode) write(io_lun,*) 'Real energy: ',e_real
    ! now the reciprocal part of the sum
    ! first calculate the structure factor
    call my_barrier()
    call structure_factor( )

    sum = zero
    do n=1, n_g_vectors
       sum = sum + structure2(n) * g_factor(n)
    end do
    e_reciprocal = two * pi * sum / volume
    !    if(inode==ionode) write(io_lun,*) 'Recip energy: ',e_reciprocal

    ! and calculate the contribution to the forces from the reciprocal 
    ! part of the sum
    call my_barrier()
    do k=1, n_my_atoms
       i = index_my_atoms(k)
       q_i = charge(species_glob(i))
       sum_x = zero
       sum_y = zero
       sum_z = zero
       do n=1, n_g_vectors
          g_dot_r = gx(n) * atom_coord(1,i) + &
               gy(n) * atom_coord(2,i) + &
               gz(n) * atom_coord(3,i) 
          cosine = cos( g_dot_r )
          sine = sin( g_dot_r )
          force = ( aimag( cstructure(n) ) * cosine +&
               real( cstructure(n), double ) * sine ) * g_factor(n)
          sum_x = sum_x + force * gx(n)
          sum_y = sum_y + force * gy(n)
          sum_z = sum_z + force * gz(n)
       end do
       c2_force(1,i) = four * pi * q_i * sum_x / volume
       c2_force(2,i) = four * pi * q_i * sum_y / volume
       c2_force(3,i) = four * pi * q_i * sum_z / volume
    end do

    ! finally there is a constant to be added 
    e_constant = zero
    total_charge = zero
    do i=1, ni_in_cell
       e_constant = e_constant + charge(species(i)) * charge(species(i))
       total_charge = total_charge + charge(species(i))
    end do
    e_constant = -gamma * e_constant / sqrt_pi - &
         pi * total_charge * total_charge / &
         ( two * gamma * gamma * volume )

    ewald_energy = e_real + e_reciprocal + e_constant
    !ewald_energy = e_real + e_reciprocal! + e_constant

    if (inode==ionode.AND.iprint_gen>1) &
         write(io_lun,1) en_conv*ewald_energy,en_units(energy_units), en_conv*e_real,en_units(energy_units), &
         en_conv*e_reciprocal,en_units(energy_units), en_conv*e_constant,en_units(energy_units)

    ! add the real-part and reciprocal-part forces
    call gsum(c1_force,3,ni_in_cell)
    call gsum(c2_force,3,ni_in_cell)
    do i=1, ni_in_cell
       ewald_force(1,i) = c1_force(1,i) + c2_force(1,i)
       ewald_force(2,i) = c1_force(2,i) + c2_force(2,i)
       ewald_force(3,i) = c1_force(3,i) + c2_force(3,i)
       !       if(inode==ionode) write(io_lun,fmt='(2x,"Real force: ",3f20.12)') c1_force(1,i),c1_force(2,i),c1_force(3,i)
       !       if(inode==ionode) write(io_lun,fmt='(2x,"Recip force: ",3f20.12)') c2_force(1,i),c2_force(2,i),c2_force(3,i)
    end do
    deallocate(c2_force,c1_force,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating c1_force in ewald: ",ni_in_cell,stat)
    call reg_dealloc_mem(area_general,6*ni_in_cell,type_dbl)
    call stop_print_timer(tmr_l_tmp1,"ewald",IPRINT_TIME_THRES3)
1   format(/20x,'E_coulomb = ',f25.17,' ',a2,/,&
         20x,'E_real = ',f25.17,' ',a2,/,&
         20x,'E_recipr. = ',f25.17,' ',a2,/,&
         20x,'E_const. = ',f25.17,' ',a2,/)
    return
  end subroutine ewald
  !!***

  ! -----------------------------------------------------------
  ! Subroutine set_ewald
  ! -----------------------------------------------------------

  !!****f* ewald_module/set_ewald *
  !!
  !!  NAME 
  !!   set_ewald
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Sets up parameters for Ewald sum.  Works out G reciprocal lattice
  !!   vectors of the simulation cell (assumes cubic/orthorhombic shape).
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   E.H.Hernandez
  !!  CREATION DATE
  !!   05/02/96
  !!  MODIFICATION HISTORY
  !!   14/05/01 DRB
  !!    Added ROBODoc header and incorporated into ewald_module
  !!   2005/10/06 09:56 dave
  !!    Bug fix: was checking rij2 not r2_min
  !!   2008/03/03 18:44 dave
  !!    Changed float to real
  !!  TODO 
  !!   Update for non-orthorhombic cells and to use partitions DRB 14/05/01
  !!  SOURCE
  !!
  subroutine set_ewald(inode,ionode)

    use datatypes
    use numbers
    use atoms
    use dimens
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: x_atom_cell, y_atom_cell, z_atom_cell, ni_in_cell, &
         atom_coord, iprint_gen, area_general
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Shared variables
    integer :: inode, ionode

    ! Local variables
    integer :: i, j, k, l, m, n_shell, n_vectors, n, nx, ny, nz, stat
    integer :: i_corner_x( 8 ), i_corner_y( 8 ), i_corner_z( 8 ), &
         j_corner_x( 8 ), j_corner_y( 8 ), j_corner_z( 8 )

    real(double) :: a_min, argument, b_min, b_x, b_y, b_z,  &
         g2, g2_min, rij2, rij, rxij, ryij, rzij,  &
         r2_min, test, vec_x, vec_y, vec_z
    !    real(double), parameter :: epsilon = 1.0d-5
    real(double) :: epsilon
    real(double), parameter :: huge = 1.0e30_double

    logical :: contributing_shell

    ! +++
    if(inode == ionode.AND.iprint_gen>1) write(io_lun,fmt='(//"+++ entering original Hernandez set_ewald +++"//)')
    ! +++

    ! first set up the reciprocal lattice vectors
    b_x = two * pi / r_super_x
    b_y = two * pi / r_super_y
    b_z = two * pi / r_super_z

    ! the best Ewald parameter gamma is given by the relation 
    !         gamma^2 = b_min / ( 2 * a_min )
    ! where a_min and b_min are the smallest direct and reciprocal 
    ! lattice vectors respectively, so the first thing is to find them
    a_min = min( r_super_x, min( r_super_y, r_super_z ))
    b_min = min( b_x, min( b_y, b_z ))
    gamma = sqrt( b_min / ( two * a_min ))

    ! set the eight corners of the central box
    n = 1
    do k=1,2
       do j=1,2
          do i=1,2
             i_corner_x(n) = i - 1 
             i_corner_y(n) = j - 1 
             i_corner_z(n) = k - 1 
             n = n + 1
          end do
       end do
    end do

    ! Count supercells
    n_supercells = 0
    n_shell = 0
    contributing_shell = .true.
    do while ( contributing_shell )
       n_shell = n_shell + 1
       contributing_shell = .false.
       do nz=-n_shell, n_shell
          do ny=-n_shell, n_shell
             do nx=-n_shell, n_shell
                ! test if this box is on the surface of the 'sphere' of boxes
                if ( (abs(nz)==n_shell ) .or. &
                     (abs(ny)==n_shell ) .or. &
                     (abs(nx)==n_shell ) ) then
                   ! find out the minimum distance between atoms in this box 
                   ! and atoms in the central simulation box (nx = ny = nz = 0)
                   ! the minimum distance is larger than the distance between 
                   ! the nearest corners of both cells, so we will take that 
                   ! distance set the eight corners of the outer cell
                   n = 1
                   do k=1,2
                      do j=1,2
                         do i=1,2
                            j_corner_x(n) = nx + i - 1 
                            j_corner_y(n) = ny + j - 1 
                            j_corner_z(n) = nz + k - 1 
                            n = n + 1
                         end do
                      end do
                   end do
                   r2_min = huge
                   do i=1, 8
                      do j=1, 8
                         rxij = real( i_corner_x(i) - j_corner_x(j), double )
                         ryij = real( i_corner_y(i) - j_corner_y(j), double )
                         rzij = real( i_corner_z(i) - j_corner_z(j), double )
                         rij2 = rxij * rxij * r_super_x_squared + &
                              ryij * ryij * r_super_y_squared + &
                              rzij * rzij * r_super_z_squared
                         !                         if ( ( rij2<=r2_min ) .and. &
                         !                              ( rij2>zero ) ) r2_min = rij2
                         if (rij2<=r2_min) r2_min = rij2
                      end do
                   end do
                   if(r2_min>zero) then
                      rij = sqrt( r2_min )
                      argument = gamma * rij 
                      test = real( ni_in_cell * ni_in_cell, double ) * &
                           erfc( argument ) / rij
                   else
                      test = zero
                   end if
                   ! if this box contributes, keep its components in a list
                   if ( test>=ewald_accuracy.OR.r2_min==zero ) then
                      contributing_shell = .true.
                      n_supercells = n_supercells + 1
                   end if ! test>=ewald_accuracy
                end if ! abs(nx)==nshell etc
             end do ! nz
          end do ! ny 
       end do ! nx
    end do ! while contributing shell
    allocate(supercell_vec_x(n_supercells), supercell_vec_y(n_supercells), &
         supercell_vec_z(n_supercells), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating supercells in set_ewald: ",n_supercells, stat)
    call reg_alloc_mem(area_general,3*n_supercells,type_dbl)
    ! now find out which supercells contribute to the real-space sum 
    n_supercells = 0
    n_shell = 0
    contributing_shell = .true.
    do while ( contributing_shell )
       n_shell = n_shell + 1
       contributing_shell = .false.
       do nz=-n_shell, n_shell
          do ny=-n_shell, n_shell
             do nx=-n_shell, n_shell
                ! test if this box is on the surface of the 'sphere' of boxes
                if ( (abs(nz)==n_shell ) .or. &
                     (abs(ny)==n_shell ) .or. &
                     (abs(nx)==n_shell ) ) then
                   ! find out the minimum distance between atoms in this box 
                   ! and atoms in the central simulation box (nx = ny = nz = 0)
                   ! the minimum distance is larger than the distance between 
                   ! the nearest corners of both cells, so we will take that 
                   ! distance set the eight corners of the outer cell
                   n = 1
                   do k=1,2
                      do j=1,2
                         do i=1,2
                            j_corner_x(n) = nx + i - 1 
                            j_corner_y(n) = ny + j - 1 
                            j_corner_z(n) = nz + k - 1 
                            n = n + 1
                         end do
                      end do
                   end do
                   r2_min = huge
                   do i=1, 8
                      do j=1, 8
                         rxij = real( i_corner_x(i) - j_corner_x(j), double )
                         ryij = real( i_corner_y(i) - j_corner_y(j), double )
                         rzij = real( i_corner_z(i) - j_corner_z(j), double )
                         rij2 = rxij * rxij * r_super_x_squared + &
                              ryij * ryij * r_super_y_squared + &
                              rzij * rzij * r_super_z_squared
                         !                         if ( ( rij2<=r2_min ) .and. &
                         !                              ( rij2>zero ) ) r2_min = rij2
                         if (rij2<=r2_min) r2_min = rij2
                      end do
                   end do
                   if(r2_min>zero) then
                      rij = sqrt( r2_min )
                      argument = gamma * rij 
                      test = real( ni_in_cell * ni_in_cell, double ) * &
                           erfc( argument ) / rij
                   else
                      test = zero
                   end if
                   ! if this box contributes, keep its components in a list
                   if ( test>=ewald_accuracy.OR.r2_min==zero ) then
                      contributing_shell = .true.
                      n_supercells = n_supercells + 1
                      supercell_vec_x(n_supercells) = real( nx, double )
                      supercell_vec_y(n_supercells) = real( ny, double )
                      supercell_vec_z(n_supercells) = real( nz, double )
                   end if ! test>=ewald_accuracy
                end if ! abs(nx)==nshell etc
             end do ! nz
          end do ! ny 
       end do ! nx
    end do ! while contributing shell
    ! now find out which is the maximum reciprocal lattice vector modulus
    ! that we need to include in the sum to include terms bigger than ewald_accuracy
    n_shell = 1
    contributing_shell = .true.
    do while ( contributing_shell )
       contributing_shell = .false.
       g2_min = huge
       n_vectors = 0
       do nz=-n_shell, n_shell
          do ny=-n_shell, n_shell
             do nx=-n_shell, n_shell
                ! test if the tip of this vector is on the surface of the 
                ! 'sphere' of G vectors
                if ( (abs(nz)==n_shell ) .or. &
                     (abs(ny)==n_shell ) .or. &
                     (abs(nx)==n_shell ) ) then
                   ! find out the modulus of this vector and see if it 
                   ! contributes
                   n_vectors = n_vectors + 1
                   vec_x = real(nx,double) * b_x
                   vec_y = real(ny,double) * b_y
                   vec_z = real(nz,double) * b_z
                   g2 = vec_x*vec_x + vec_y*vec_y + vec_z*vec_z
                   if ( g2<=g2_min ) g2_min = g2
                end if ! abs(nx)==n_shell etc
             end do ! nx
          end do ! ny
       end do ! nz
       ! if this shell of G vectors has a contribution, then include it
       test = real( n_vectors, double ) *  &
            exp( -g2_min / (four * gamma * gamma) ) / g2_min
       if ( test>=ewald_accuracy ) then
          contributing_shell = .true.
          n_shell = n_shell + 1
       else
          n_shell = n_shell - 1
       end if
    end do ! contributing shell
    ! Count g-vectors
    n_g_vectors = 0
    do i=1, n_shell
       do m=-i, i
          do l=-i, i
             do k=-i, i
                if ( (abs(m)==i) .or.  &
                     (abs(l)==i) .or. &
                     (abs(k)==i) ) then
                   n_g_vectors = n_g_vectors + 1
                end if
             end do ! k=-i,i
          end do ! l=-i,i
       end do ! m=-i,i
    end do ! i=1,n_shell
    allocate(gx(n_g_vectors), gy(n_g_vectors), gz(n_g_vectors), g_factor(n_g_vectors), structure2(n_g_vectors), &
         cstructure(n_g_vectors), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating gx, gy etc in set_ewald: ",n_g_vectors,stat)
    call reg_alloc_mem(area_general,6*n_g_vectors,type_dbl)
    ! now use this information to construct the actual G vectors
    n_g_vectors = 0
    do i=1, n_shell
       do m=-i, i
          do l=-i, i
             do k=-i, i
                if ( (abs(m)==i) .or.  &
                     (abs(l)==i) .or. &
                     (abs(k)==i) ) then
                   n_g_vectors = n_g_vectors + 1
                   gx(n_g_vectors) = real( k, double ) * b_x
                   gy(n_g_vectors) = real( l, double ) * b_y
                   gz(n_g_vectors) = real( m, double ) * b_z
                   g2 =  real( k * k, double) * b_x * b_x + &
                        &real( l * l, double) * b_y * b_y + &
                        &real( m * m, double) * b_z * b_z
                   argument = - g2 / ( four * gamma * gamma )
                   g_factor(n_g_vectors) = exp( argument ) / g2
                end if
             end do ! k=-i,i
          end do ! l=-i,i
       end do ! m=-i,i
    end do ! i=1,n_shell
    allocate(ewald_force(3,ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating ewald_force in set_ewald: ",ni_in_cell,stat)
    call reg_alloc_mem(area_general,3*ni_in_cell,type_dbl)
    ! write out some information
    if ( inode .eq. ionode .AND.iprint_gen>0)  &
         write(io_lun,1) ewald_accuracy, gamma, n_supercells, n_g_vectors
1   format(/10x,'PrecisionQ required of the Ewald summation = ',g15.7, &
         /10x,'Ewald parameter gamma = ',f10.5, &
         /10x,'With this choice of ewald_accuracy and gamma we have that', &
         /10x,i5,' direct space supercells and ', &
         /10x,i5,' reciprocal lattice vectors ', &
         /10x,'need to be included in the direct and reciprocal parts' &
         /10x,'of the Ewald sumation respectively.'/)
    ! all done so return
    return
  end subroutine set_ewald
  !!***

  ! -----------------------------------------------------------
  ! Function obsolete_erfcc
  ! -----------------------------------------------------------

  !!****f* ewald_module/obsolete_erfcc *
  !!
  !!  NAME 
  !!   obsolete_erfcc
  !!  USAGE
  !!   obsolete_erfcc(x)
  !!  PURPOSE
  !!   Finds the complementary error function with fractional error
  !!   less than 1.2e-7 everywhere - taken from Numerical Recipes.
  !!   A fractional error of 1.2e-7 is inadequate for Conquest, and
  !!   this routine is kept here only as a museum exhibit. It is
  !!   not intended to be called.
  !!  INPUTS
  !!   real(double) :: x - point at which to evaluate function
  !!  USES
  !! 
  !!  AUTHOR
  !!   E.H.Hernandez
  !!  CREATION DATE
  !!   Maybe around 1996?
  !!  MODIFICATION HISTORY
  !!   14/05/01 DRB
  !!    Incorporated into ewald_module, and updated to F90 format
  !!    ROBODoc header added
  !!  TODO
  !! 
  !!  SOURCE
  !!
  real(double) function obsolete_erfcc(x)

    use datatypes
    use numbers

    implicit none

    real(double) :: x
    real(double) :: argument, t, z

    z = abs(x)
    t = one / ( one + half * z )

    argument = - z * z -1.26551223_double + t * ( 1.00002368_double + &
         t * ( 0.37409196_double + t * ( 0.09678418_double + &
         t * (-0.18628806_double + t * ( 0.27886807_double + &
         t * (-1.13520398_double + t * ( 1.48851587_double + &
         t * (-0.82215223_double + t * 0.17087277_double))))))))

    obsolete_erfcc = t * exp( argument )

    if (x < zero) obsolete_erfcc = two - obsolete_erfcc

    return
  end function obsolete_erfcc
  !!***

  ! -----------------------------------------------------------
  ! Function erfc
  ! -----------------------------------------------------------

  !!****f* ewald_module/erfc *
  !!
  !!  NAME 
  !!   erfc
  !!  USAGE
  !!   erfc(x)
  !!  PURPOSE
  !!   Calculated the complementary error function to rounding-error
  !!   accuracy
  !!  INPUTS
  !!   real(double) :: argument of complementary error function
  !!  USES
  !! 
  !!  AUTHOR
  !!   Iain Chapman
  !!  CREATION DATE
  !!   Some time in 2001
  !!  MODIFICATION HISTORY
  !!   14/05/01 DRB
  !!    Put into Conquest by MJG 10/2003, various tidying done at this time
  !!  TODO
  !! 
  !!  SOURCE
  !!
  real(double) function erfc(x)
    use datatypes
    use numbers, ONLY: very_small, one
    use GenComms, ONLY: cq_abort

    real(double), parameter :: erfc_delta = 1.0e-12_double, erfc_gln = 0.5723649429247447e0_double, &
         &erfc_fpmax =1.e30_double
    integer, parameter:: erfc_iterations = 10000

    real(double), intent(in)::x

    ! Local variables
    real(double) :: x2
    real(double) :: ap, sum, del
    real(double) :: an, b,c,d,h

    integer :: i

    ! This expects x^2...
    x2 = x*x
    !if(x2<very_small) call cq_abort("Error in ewald: x2 small")
    if(x<very_small) then
       erfc = one
       return
    end if
    !if x2 less than (1.0 + 0.5) squared
    if (x2 < 2.25_double) then
       ap = 0.5_double
       sum = 2.0_double
       del = sum
       do i = 1, erfc_iterations
          ap = ap + 1.0_double
          del = del * x2 /ap
          sum = sum + del
          if (abs(del) < abs(sum) * erfc_delta) exit
       end do
       erfc = 1.0_double - sum * exp(-x2 + 0.5_double * log(x2) - erfc_gln)

    else
       b = x2 + 0.5_double
       c = erfc_fpmax
       d = 1.0_double / b
       sum = d

       do i = 1, erfc_iterations
          an = - i * (i - 0.5_double)
          b = b + 2.0_double
          d = an * d + b
          !  if (abs(d) < fpmin) d = fpmin
          c = b + an / c
          !  if (abs(c) < fpmin) c = fpmin
          d = 1.0_double / d
          del = d * c
          sum = sum * del         
          if (abs(del - 1.0_double) < erfc_delta) exit
       end do
       erfc = sum * exp(-x2 + 0.5_double * log(x2) - erfc_gln)
    end if
    return

  end function erfc
  !!***

  !!****f* ewald_module/structure_factor *
  !!
  !!  NAME 
  !!   structure_factor
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Calculates the structure factor and its squared modulus, needed
  !!    for the reciprocal lattice part of the Ewald summation
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   E.H.Hernandez
  !!  CREATION DATE
  !!   20/02/96 
  !!  MODIFICATION HISTORY
  !!   14/05/01 DRB
  !!    Added ROBODoc header, stripped shared variables and updated to
  !!    F90 format
  !!   08/06/2001 dave
  !!    Changed zgsum to gsumv and used GenComms
  !!    Changed gsumv to gsum
  !!     2008/02/06 08:00 dave & veronika
  !!      Subtle bug fix: species_glob not species
  !!  SOURCE
  !!
  subroutine structure_factor( )

    use numbers
    use atoms
    use species_module
    use dimens
    use GenComms, ONLY: gsum
    use global_module, ONLY: atom_coord, species_glob

    implicit none

    ! Local variables
    integer i, k, n
    real(double) :: g_dot_r, q_i
    complex(double_cplx) :: csum

    ! calculate the charge-weighted structure factor 
    cstructure = (zero, zero)
    do n=1, n_g_vectors
       csum = (zero, zero)
       do k=1, n_my_atoms
          i = index_my_atoms(k)
          q_i = charge(species_glob(i))
          g_dot_r = gx(n) * atom_coord(1,i) + &
               gy(n) * atom_coord(2,i) + &
               gz(n) * atom_coord(3,i) 
          csum = csum + q_i * exp( cmplx(zero, -one, double_cplx) * g_dot_r )
       end do
       cstructure(n) = csum
    end do
    call gsum( cstructure, n_g_vectors )
    ! now get the squared modulus of the structure factor
    do n=1, n_g_vectors
       structure2(n) = cstructure(n) * conjg( cstructure(n) )
    end do
    return
  end subroutine structure_factor
  !!***

  ! =============================================================================
  !   New routines under development by Mike, starting 04/09/03
  ! =============================================================================

! Bug fix, 2005/10/06 10:13 dave
!  It was assumed that a1_dot_a1 etc and b1_dot_b1 etc would remain defined after exit
!  from routine; I added an explicit recalculation if pd_init .false.
  ! ==========================================================
  subroutine partition_distance(pd_init,aa,beta,rr,distance)

    use numbers
    use GenComms, ONLY: cq_abort, myid

    implicit none

    logical, intent(in) :: pd_init
    real(double), intent(in), optional, dimension(3) :: rr
    real(double), intent(in), dimension(3,3) :: aa
    real(double), intent(inout), dimension(3,3) :: beta
    real(double), intent(out), optional :: distance

    integer :: i, n
    integer :: isig1, isig2, isig3
    real(double) :: alpha1, alpha2, alpha3
    real(double) :: asig1, asig2, asig3, a1_dot_a1, a2_dot_a2, a3_dot_a3, &
         &a2_dot_a3, a3_dot_a1, a1_dot_a2, b1_dot_b1, b2_dot_b2, b3_dot_b3, &
         &b2_dot_b3, b3_dot_b1, b1_dot_b2, dx, dy, dz, p1, p2, p3, x1, x2, x3, &
         &y1, y2, y3, z1, z2, z3
    real(double), dimension(3) :: edg_a, edg_b, q_vec


    ! --- initialisation for given cell vectors: do only if pd_int = .true.
    if(pd_init) then
       ! --- construct required quantities from cell vectors -----

       a1_dot_a1 = aa(1,1)*aa(1,1) + aa(1,2) * aa(1,2) + aa(1,3)*aa(1,3)
       a2_dot_a2 = aa(2,1)*aa(2,1) + aa(2,2) * aa(2,2) + aa(2,3)*aa(2,3)
       a3_dot_a3 = aa(3,1)*aa(3,1) + aa(3,2) * aa(3,2) + aa(3,3)*aa(3,3)
       a2_dot_a3 = aa(2,1)*aa(3,1) + aa(2,2) * aa(3,2) + aa(2,3)*aa(3,3)
       a3_dot_a1 = aa(3,1)*aa(1,1) + aa(3,2) * aa(1,2) + aa(3,3)*aa(1,3)
       a1_dot_a2 = aa(1,1)*aa(2,1) + aa(1,2) * aa(2,2) + aa(1,3)*aa(2,3)

       beta(1,1) = aa(2,2)*aa(3,3) - aa(2,3)*aa(3,2)
       beta(1,2) = aa(2,3)*aa(3,1) - aa(2,1)*aa(3,3)
       beta(1,3) = aa(2,1)*aa(3,2) - aa(2,2)*aa(3,1)
       beta(2,1) = aa(3,2)*aa(1,3) - aa(3,3)*aa(1,2)
       beta(2,2) = aa(3,3)*aa(1,1) - aa(3,1)*aa(1,3)
       beta(2,3) = aa(3,1)*aa(1,2) - aa(3,2)*aa(1,1)
       beta(3,1) = aa(1,2)*aa(2,3) - aa(1,3)*aa(2,2)
       beta(3,2) = aa(1,3)*aa(2,1) - aa(1,1)*aa(2,3)
       beta(3,3) = aa(1,1)*aa(2,2) - aa(1,2)*aa(2,1)


       alpha1 = (aa(1,1)*beta(1,1) + aa(1,2)*beta(1,2) + aa(1,3)*beta(1,3)) / &
            &(beta(1,1)*beta(1,1) + beta(1,2)*beta(1,2) + beta(1,3)*beta(1,3))
       alpha2 = (aa(2,1)*beta(2,1) + aa(2,2)*beta(2,2) + aa(2,3)*beta(2,3)) / &
            &(beta(2,1)*beta(2,1) + beta(2,2)*beta(2,2) + beta(2,3)*beta(2,3))
       alpha3 = (aa(3,1)*beta(3,1) + aa(3,2)*beta(3,2) + aa(3,3)*beta(3,3)) / &
            &(beta(3,1)*beta(3,1) + beta(3,2)*beta(3,2) + beta(3,3)*beta(3,3))

       beta(1,1) = alpha1*beta(1,1)
       beta(1,2) = alpha1*beta(1,2)
       beta(1,3) = alpha1*beta(1,3)
       beta(2,1) = alpha2*beta(2,1)
       beta(2,2) = alpha2*beta(2,2)
       beta(2,3) = alpha2*beta(2,3)
       beta(3,1) = alpha3*beta(3,1)
       beta(3,2) = alpha3*beta(3,2)
       beta(3,3) = alpha3*beta(3,3)

       b1_dot_b1 = beta(1,1)*beta(1,1) + beta(1,2)*beta(1,2) + beta(1,3)*beta(1,3)
       b2_dot_b2 = beta(2,1)*beta(2,1) + beta(2,2)*beta(2,2) + beta(2,3)*beta(2,3)
       b3_dot_b3 = beta(3,1)*beta(3,1) + beta(3,2)*beta(3,2) + beta(3,3)*beta(3,3)
       b2_dot_b3 = beta(2,1)*beta(3,1) + beta(2,2)*beta(3,2) + beta(2,3)*beta(3,3)
       b3_dot_b1 = beta(3,1)*beta(1,1) + beta(3,2)*beta(1,2) + beta(3,3)*beta(1,3)
       b1_dot_b2 = beta(1,1)*beta(2,1) + beta(1,2)*beta(2,2) + beta(1,3)*beta(2,3)
    else
       a1_dot_a1 = aa(1,1)*aa(1,1) + aa(1,2) * aa(1,2) + aa(1,3)*aa(1,3)
       a2_dot_a2 = aa(2,1)*aa(2,1) + aa(2,2) * aa(2,2) + aa(2,3)*aa(2,3)
       a3_dot_a3 = aa(3,1)*aa(3,1) + aa(3,2) * aa(3,2) + aa(3,3)*aa(3,3)
       a2_dot_a3 = aa(2,1)*aa(3,1) + aa(2,2) * aa(3,2) + aa(2,3)*aa(3,3)
       a3_dot_a1 = aa(3,1)*aa(1,1) + aa(3,2) * aa(1,2) + aa(3,3)*aa(1,3)
       a1_dot_a2 = aa(1,1)*aa(2,1) + aa(1,2) * aa(2,2) + aa(1,3)*aa(2,3)
       b1_dot_b1 = beta(1,1)*beta(1,1) + beta(1,2)*beta(1,2) + beta(1,3)*beta(1,3)
       b2_dot_b2 = beta(2,1)*beta(2,1) + beta(2,2)*beta(2,2) + beta(2,3)*beta(2,3)
       b3_dot_b3 = beta(3,1)*beta(3,1) + beta(3,2)*beta(3,2) + beta(3,3)*beta(3,3)
       b2_dot_b3 = beta(2,1)*beta(3,1) + beta(2,2)*beta(3,2) + beta(2,3)*beta(3,3)
       b3_dot_b1 = beta(3,1)*beta(1,1) + beta(3,2)*beta(1,2) + beta(3,3)*beta(1,3)
       b1_dot_b2 = beta(1,1)*beta(2,1) + beta(1,2)*beta(2,2) + beta(1,3)*beta(2,3)
    endif
    if(pd_init) return

    ! --- components of supplied position ----------------------
    x1 = (beta(1,1)*rr(1) + beta(1,2)*rr(2) + beta(1,3)*rr(3))/b1_dot_b1
    x2 = (beta(2,1)*rr(1) + beta(2,2)*rr(2) + beta(2,3)*rr(3))/b2_dot_b2
    x3 = (beta(3,1)*rr(1) + beta(3,2)*rr(2) + beta(3,3)*rr(3))/b3_dot_b3

    if(abs(x1)>1e5.AND.myid==0) then
       write(io_lun,*) 'Error: betas, rs: ',beta(1,:),rr(1)
       write(io_lun,*) 'Error: betas, rs: ',beta(2,:),rr(2)
       write(io_lun,*) 'Error: betas, rs: ',beta(3,:),rr(3)
    end if
    y1 = (aa(1,1)*rr(1) + aa(1,2)*rr(2) + aa(1,3)*rr(3))/b1_dot_b1
    y2 = (aa(2,1)*rr(1) + aa(2,2)*rr(2) + aa(2,3)*rr(3))/b2_dot_b2
    y3 = (aa(3,1)*rr(1) + aa(3,2)*rr(2) + aa(3,3)*rr(3))/b3_dot_b3

    z2 = two*x2 - (two*x3-sign(min(abs(two*x3),one),x3))*b2_dot_b3/b2_dot_b2
    z3 = two*x3 - (two*x2-sign(min(abs(two*x2),one),x2))*b2_dot_b3/b3_dot_b3
    p1 = (two/a1_dot_a1) * (y1*b1_dot_b1 - &
         &half*a1_dot_a2*sign(min(abs(z2),one),z2) -&
         &half*a3_dot_a1*sign(min(abs(z3),one),z3))
    z3 = two*x3 - (two*x1-sign(min(abs(two*x1),one),x1))*b3_dot_b1/b3_dot_b3
    z1 = two*x1 - (two*x3-sign(min(abs(two*x3),one),x3))*b3_dot_b1/b1_dot_b1
    p2 = (two/a2_dot_a2) * (y2*b2_dot_b2 - &
         &half*a2_dot_a3*sign(min(abs(z3),one),z3) -&
         &half*a1_dot_a2*sign(min(abs(z1),one),z1))
    z1 = two*x1 - (two*x2-sign(min(abs(two*x2),one),x2))*b1_dot_b2/b1_dot_b1
    z2 = two*x2 - (two*x1-sign(min(abs(two*x1),one),x1))*b1_dot_b2/b2_dot_b2
    p3 = (two/a3_dot_a3) * (y3*b3_dot_b3 - &
         &half*a3_dot_a1*sign(min(abs(z1),one),z1) -&
         &half*a2_dot_a3*sign(min(abs(z2),one),z2))

    isig1 = sign(min(1,floor(abs(p1))),nint(p1))
    isig2 = sign(min(1,floor(abs(p2))),nint(p2))
    isig3 = sign(min(1,floor(abs(p3))),nint(p3))

      !write(io_lun,fmt='(/" info for supplied position:"/)')
      !write(io_lun,fmt='(" x1, x2, x3:",3f12.6)') x1, x2, x3
      !write(io_lun,fmt='(" y1, y2, y3:",3f12.6)') y1, y2, y3
      !write(io_lun,fmt='(" p1, p2, p3:",3f12.6)') p1, p2, p3
      !write(io_lun,fmt='(/" signatures of position:",3i3)') isig1, isig2, isig3

    ! --- calculate shortest distance --------------------------
    asig1 = real(isig1,double)
    asig2 = real(isig2,double)
    asig3 = real(isig3,double)

    select case(abs(isig1)+abs(isig2)+abs(isig3))
       ! --- vertex: no zero signatures
    case(3)
       dx = rr(1) - half*(asig1*aa(1,1) + asig2*aa(2,1) + asig3*aa(3,1))
       dy = rr(2) - half*(asig1*aa(1,2) + asig2*aa(2,2) + asig3*aa(3,2))
       dz = rr(3) - half*(asig1*aa(1,3) + asig2*aa(2,3) + asig3*aa(3,3))
       distance = sqrt(dx*dx+dy*dy+dz*dz)
       ! --- edge: only one zero signature
    case(2)
       edg_a(1) = half*(asig1*aa(1,1) + asig2*aa(2,1) + asig3*aa(3,1)) - rr(1)
       edg_a(2) = half*(asig1*aa(1,2) + asig2*aa(2,2) + asig3*aa(3,2)) - rr(2)
       edg_a(3) = half*(asig1*aa(1,3) + asig2*aa(2,3) + asig3*aa(3,3)) - rr(3)
       edg_b(1) = (one-abs(asig1))*aa(1,1) + (one-abs(asig2))*aa(2,1) + &
            &(one-abs(asig3))*aa(3,1)
       edg_b(2) = (one-abs(asig1))*aa(1,2) + (one-abs(asig2))*aa(2,2) + &
            &(one-abs(asig3))*aa(3,2)
       edg_b(3) = (one-abs(asig1))*aa(1,3) + (one-abs(asig2))*aa(2,3) + &
            &(one-abs(asig3))*aa(3,3)
       !ORIdistance = sqrt(edg_a(1)*edg_a(1) + edg_a(2)*edg_a(2) + edg_a(3)*edg_a(3) - &
       !ORI     &((edg_a(1)*edg_b(1) + edg_a(2)*edg_b(2) + edg_a(3)*edg_b(3))**2) / &
       !ORI     &(edg_b(1)*edg_b(1) + edg_b(2)*edg_b(2) + edg_b(3)*edg_b(3)))
       distance = edg_a(1)*edg_a(1) + edg_a(2)*edg_a(2) + edg_a(3)*edg_a(3) - &
            &((edg_a(1)*edg_b(1) + edg_a(2)*edg_b(2) + edg_a(3)*edg_b(3))**2) / &
            &(edg_b(1)*edg_b(1) + edg_b(2)*edg_b(2) + edg_b(3)*edg_b(3))
       if(distance > very_small) then 
           distance = sqrt(distance)
       else
           distance = zero
       endif

       ! --- face: only one non-zero signature
    case(1)
       q_vec(1) = half*(asig1*beta(1,1) + asig2*beta(2,1) + asig3*beta(3,1))
       q_vec(2) = half*(asig1*beta(1,2) + asig2*beta(2,2) + asig3*beta(3,2))
       q_vec(3) = half*(asig1*beta(1,3) + asig2*beta(2,3) + asig3*beta(3,3))
       distance = abs(q_vec(1)*q_vec(1) + q_vec(2)*q_vec(2) + q_vec(3)*q_vec(3) - &
            &q_vec(1)*rr(1) - q_vec(2)*rr(2) - q_vec(3)*rr(3)) / &
            &sqrt(q_vec(1)*q_vec(1) + q_vec(2)*q_vec(2) + q_vec(3)*q_vec(3))
       ! --- body: all zero signatures
    case(0)
       distance = zero
       ! --- case default included as a check on coding
    case default
       write(io_lun,*) 'sigs: ',isig1,isig2,isig3
       call cq_abort('partition_distance: default found ',abs(isig1)+abs(isig2)+abs(isig3))
    end select

  end subroutine partition_distance
  ! ==========================================================

!!   2009/07/08 16:58 dave
!!    Removed double loop over all atoms to find lambda as it's not used

  ! ==========================================================
  subroutine mikes_set_ewald(inode,ionode)

    use construct_module, ONLY : init_cover
    use cover_module, ONLY : ewald_CS, make_cs
    use datatypes
    use dimens
    use GenComms, ONLY : cq_abort, my_barrier
    use group_module, ONLY : parts
    use global_module, ONLY : x_atom_cell, y_atom_cell, z_atom_cell, iprint_gen, ni_in_cell, area_general
    use numbers
    use primary_module, ONLY : bundle
    use species_module, ONLY : charge, species
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl


    implicit none

    integer, intent(in) :: inode, ionode

    integer :: i, ind_part, ip, m_limit_1, m_limit_2, m_limit_3, m1, m2, m3, n, nc, nccx, nccy, nccz, ni,&
         &np, nppx, nppy, nppz, np1, np2, rec_lim_1, rec_lim_2, rec_lim_3, run_mode, stat
    real(double) :: ewald_recip_cutoff_squared, mean_number_density, mean_interparticle_distance, &
         &mean_square_charge, partition_radius_squared, partition_diameter, total_charge
    real(double) :: argument, const_phi_0, const_phi_1, const_zeta_0, const_zeta_1, distance, &
         &dummy1, dummy2, dummy3, g1, g2, g3, g_squared, lambda, rsq, sep, vol_prefac, v1, v2, v3
    real(double), dimension(3) :: rr, vert
!    real(double), dimension(3,3) :: abs_sep, rel_sep
! 16:52, 2003/10/22 dave Changed to allow compilation under ifc
    real(double), dimension(3) :: abs_sep, rel_sep
    real(double), dimension(3,3) :: real_cell_vec, recip_cell_vec, part_cell_vec, part_cell_dual

    if(inode == ionode.AND.iprint_gen>1) write(io_lun,fmt='(/8x," edge lengths of orthorhombic supercell:", &
         &3f12.6," a.u.")') r_super_x, r_super_y, r_super_z

    ! --- setting primitive translation vectors of cell ----------------------------
    !   At time of writing, Conquest is restricted to orthorhombic cells. However
    !   the present Ewald routine is written for general cells. The variables
    !   r_super_x, r_super_y, r_super_z are the edge lengths of the orthorhombic
    !   cell. These are used to make the general real-space primitive translation vectors, which in
    !   turn are used to make the reciprocal-space vectors. Here, real_cell_vec(n,i) and
    !   recip_cell_vec(n,i) are the Cartesian components (labelled by i) of the real- and reciprocal-space
    !   vectors (labelled by n).
    real_cell_vec(1,1) = r_super_x
    real_cell_vec(1,2) = zero
    real_cell_vec(1,3) = zero
    real_cell_vec(2,1) = zero
    real_cell_vec(2,2) = r_super_y
    real_cell_vec(2,3) = zero
    real_cell_vec(3,1) = zero
    real_cell_vec(3,2) = zero
    real_cell_vec(3,3) = r_super_z
    if(inode == ionode.AND.iprint_gen>1) then
       write(io_lun,fmt='(/8x," mikes_set_ewald: cartesian components of real-cell &
            &primitive translations vectors:"/)')
       do n = 1, 3
          write(io_lun,fmt='(8x," vector no:",i3,":",3f12.6)') n, (real_cell_vec(n,i), i = 1, 3)
       end do
    end if

    ! --- volume of real-space unit cell ------------------------------
    ewald_real_cell_volume = &
         &real_cell_vec(1,1)*(real_cell_vec(2,2)*real_cell_vec(3,3) - &
         &real_cell_vec(2,3)*real_cell_vec(3,2)) + &
         &real_cell_vec(1,2)*(real_cell_vec(2,3)*real_cell_vec(3,1) - &
         &real_cell_vec(2,1)*real_cell_vec(3,3)) + &
         &real_cell_vec(1,3)*(real_cell_vec(2,1)*real_cell_vec(3,2) - &
         &real_cell_vec(2,2)*real_cell_vec(3,1))
    ewald_real_cell_volume = abs(ewald_real_cell_volume)
    if(inode==ionode.AND.iprint_gen>1) &
         write(io_lun,fmt='(/8x," volume of real-space cell:",f19.6)') ewald_real_cell_volume

    ! --- quantities needed for calculation of Ewald gamma, and real_space and recip-space cutoffs
    total_charge = zero
    mean_square_charge = zero
    do np = 1, ni_in_cell
       total_charge = total_charge + charge(species(np))
       mean_square_charge = mean_square_charge + charge(species(np))*charge(species(np))
    enddo
    mean_square_charge = mean_square_charge/real(ni_in_cell,double)
    mean_number_density = real(ni_in_cell,double)/ewald_real_cell_volume
    mean_interparticle_distance = (three/(four*pi*mean_number_density))**one_third
    ! Rescale ewald_accuracy to a dimensionless number
    ewald_accuracy = ewald_accuracy*mean_interparticle_distance/mean_square_charge
    ! --- miscellaneous constants used for making Ewald gamma and the real-space and recip-space cutoffs
    const_zeta_0 = (pi**five_sixths)/(six**one_third)
    const_zeta_1 = (two**two_thirds)*(pi**five_sixths)/(three**one_third)
    argument = min(exp(-one),&
         &(const_zeta_0*ewald_accuracy)/(real(ni_in_cell,double)**one_third))
    const_phi_0 = sqrt(-log(argument))
    argument = min(exp(-one),&
         &(const_zeta_1*ewald_accuracy)/(real(ni_in_cell,double)**five_sixths))
    const_phi_1 = sqrt(-log(argument))
    if(inode == ionode.AND.iprint_gen>1) then
       write(unit=io_lun,fmt='(/8x," const_zeta_0:",e20.12," const_zeta_1:",e20.12/&
            &8x," const_phi_0:",e20.12," const_phi_1:",e20.12/&
            &8x," total charge:",e20.12," mean-square charge:",e20.12/&
            &8x," mean number density:",e20.12," mean interparticle distance:",e20.12)') &
            &const_zeta_0, const_zeta_1, const_phi_0, const_phi_1, total_charge, mean_square_charge, &
            &mean_number_density, mean_interparticle_distance
    endif



    ! *************************************************************************************************
    !                          SET-UP INFORMATION FOR EWALD RECIPROCAL-SPACE SUM
    ! *************************************************************************************************

    ! --- reciprocal-lattice vectors ----------------------------------
    vol_prefac = 2.0_double*pi/ewald_real_cell_volume  

    recip_cell_vec(1,1) = vol_prefac*& 
         &(real_cell_vec(2,2)*real_cell_vec(3,3) - real_cell_vec(2,3)*real_cell_vec(3,2))
    recip_cell_vec(1,2) = vol_prefac*& 
         &(real_cell_vec(2,3)*real_cell_vec(3,1) - real_cell_vec(2,1)*real_cell_vec(3,3))
    recip_cell_vec(1,3) = vol_prefac*& 
         &(real_cell_vec(2,1)*real_cell_vec(3,2) - real_cell_vec(2,2)*real_cell_vec(3,1))

    recip_cell_vec(2,1) = vol_prefac*& 
         &(real_cell_vec(3,2)*real_cell_vec(1,3) - real_cell_vec(3,3)*real_cell_vec(1,2))
    recip_cell_vec(2,2) = vol_prefac*& 
         &(real_cell_vec(3,3)*real_cell_vec(1,1) - real_cell_vec(3,1)*real_cell_vec(1,3))
    recip_cell_vec(2,3) = vol_prefac*& 
         &(real_cell_vec(3,1)*real_cell_vec(1,2) - real_cell_vec(3,2)*real_cell_vec(1,1))

    recip_cell_vec(3,1) = vol_prefac*& 
         &(real_cell_vec(1,2)*real_cell_vec(2,3) - real_cell_vec(1,3)*real_cell_vec(2,2))
    recip_cell_vec(3,2) = vol_prefac*& 
         &(real_cell_vec(1,3)*real_cell_vec(2,1) - real_cell_vec(1,1)*real_cell_vec(2,3))
    recip_cell_vec(3,3) = vol_prefac*& 
         &(real_cell_vec(1,1)*real_cell_vec(2,2) - real_cell_vec(1,2)*real_cell_vec(2,1))

    if(inode==ionode.AND.iprint_gen>1) then
       write(unit=io_lun,fmt='(/8x," cartesian components of reciprocal-space &
            &primitive vectors:"/)')
       do n = 1, 3
          write(unit=io_lun,fmt='(8x," basis vector no:",i3,":",3f12.6)') &
               &n, (recip_cell_vec(n,i), i = 1, 3)
       enddo
    end if

    ! ++++++ temporary code to calculate lower bound on inter-particle separation
    ! --- radius of sphere inscribed in real-space cell --------------
    dummy1 = pi/sqrt(recip_cell_vec(1,1)**2 + recip_cell_vec(1,2)**2 + &
         &recip_cell_vec(1,3)**2)
    dummy2 = pi/sqrt(recip_cell_vec(2,1)**2 + recip_cell_vec(2,2)**2 + &
         &recip_cell_vec(2,3)**2)
    dummy3 = pi/sqrt(recip_cell_vec(3,1)**2 + recip_cell_vec(3,2)**2 + &
         &recip_cell_vec(3,3)**2)
    lambda = min(dummy1,dummy2,dummy3)
    if(inode==ionode.AND.iprint_gen>1) write(unit=io_lun,fmt='(/8x," radius of sphere inscribed in real-space cell:",&
         &f12.6)') lambda

    if(inode == ionode.AND.iprint_gen>1) write(unit=io_lun,fmt='(/8x," no. of atoms in cell:",i5)') ni_in_cell

    if(.false.) then ! ni_in_cell > 1) then
       do np1 = 1, ni_in_cell-1
          do np2 = np1+1, ni_in_cell
             abs_sep(1) = x_atom_cell(np1) - x_atom_cell(np2)
             abs_sep(2) = y_atom_cell(np1) - y_atom_cell(np2)
             abs_sep(3) = z_atom_cell(np1) - z_atom_cell(np2)
             do n = 1, 3
                rel_sep(n) = zero
                do i = 1, 3
                   rel_sep(n) = rel_sep(n) + recip_cell_vec(n,i)*abs_sep(i)
                enddo
                rel_sep(n) = rel_sep(n)/(two*pi)
                rel_sep(n) = rel_sep(n) - anint(rel_sep(n))
             enddo
             do i = 1, 3
                abs_sep(i) = zero
                do n = 1, 3
                   abs_sep(i) = abs_sep(i) + rel_sep(n)*real_cell_vec(n,i)
                enddo
                sep = sqrt(abs_sep(1)**2 + abs_sep(2)**2 + abs_sep(3)**2)
                lambda = min(sep,lambda)
             enddo
          enddo
       enddo
    endif
    if(inode == ionode.AND.iprint_gen>1) &
         write(unit=io_lun,fmt='(/8x," lower bound on inter-particle separation:",f12.6)') lambda

    run_mode = 0
    select case(run_mode)

    case(0)
       ewald_gamma = (pi*(mean_number_density**two_thirds)/real(ni_in_cell,double)**one_third)*&
            &(const_phi_0/const_phi_1)
       argument = min(exp(-one),&
            &(ewald_accuracy*ewald_gamma)/(two*sqrt(pi)*mean_interparticle_distance*mean_number_density))
       ewald_real_cutoff = (one/sqrt(ewald_gamma))*sqrt(-log(argument))
       argument = min(exp(-one),&
            &(pi*ewald_accuracy)/(mean_interparticle_distance*real(ni_in_cell,double)*sqrt(ewald_gamma)))
       ewald_recip_cutoff = two*sqrt(ewald_gamma)*sqrt(-log(argument))

    case(1)
       ewald_gamma = (one/lambda**2)*log(one/ewald_accuracy) 
       ewald_real_cutoff = lambda
       ewald_recip_cutoff = sqrt(four*ewald_gamma*log(one/ewald_accuracy))

    end select

    if(inode == ionode.AND.iprint_gen>1) then
       write(unit=io_lun,fmt='(/8x," run_mode:",i3/" Ewald gamma:",e20.12/&
            &" Ewald real_space cutoff distance:",e20.12/" Ewald recip-space cutoff wavevector:",e20.12)') &
            &run_mode, ewald_gamma, ewald_real_cutoff, ewald_recip_cutoff
    endif

    ! --- Gaussian self-energy and net charge contributions to the Ewald energy
    ewald_gaussian_self_energy = -sqrt(ewald_gamma/pi)*real(ni_in_cell,double)*mean_square_charge
    ewald_net_charge_energy = -(pi/(two*ewald_real_cell_volume*ewald_gamma))*total_charge*total_charge
    if(inode == ionode.AND.iprint_gen>1) &
         write(unit=io_lun,fmt='(/8x," mikes_set_ewald: Ewald Gaussian self-energy:",e20.12,&
         &" Ewald net-charge energy:",e20.12)') ewald_gaussian_self_energy, ewald_net_charge_energy

    ! --- limits for recip-space covering set
    rec_lim_1 = int(ewald_recip_cutoff*&
         &sqrt(real_cell_vec(1,1)**2 + real_cell_vec(1,2)**2 + real_cell_vec(1,3)**2)/(two*pi))
    rec_lim_2 = int(ewald_recip_cutoff*&
         &sqrt(real_cell_vec(2,1)**2 + real_cell_vec(2,2)**2 + real_cell_vec(2,3)**2)/(two*pi))
    rec_lim_3 = int(ewald_recip_cutoff*&
         &sqrt(real_cell_vec(3,1)**2 + real_cell_vec(3,2)**2 + real_cell_vec(3,3)**2)/(two*pi))
    if(inode == ionode.AND.iprint_gen>1) &
         write(unit=io_lun,fmt='(/8x," integer limits for recip-space covering set:",3i5)') &
         &rec_lim_1, rec_lim_2, rec_lim_3

    ! --- determine and store wavevectors to be included in recip-space sum
    ewald_recip_cutoff_squared = ewald_recip_cutoff*ewald_recip_cutoff
    number_of_g_vectors = 0
    do m1 = -rec_lim_1, rec_lim_1
       do m2 = -rec_lim_2, rec_lim_2
          do m3 = -rec_lim_3, rec_lim_3
             if(m1*m1+m2*m2+m3*m3 /= 0) then
                g1 = m1*recip_cell_vec(1,1) + m2*recip_cell_vec(2,1) + m3*recip_cell_vec(3,1)
                g2 = m1*recip_cell_vec(1,2) + m2*recip_cell_vec(2,2) + m3*recip_cell_vec(3,2)
                g3 = m1*recip_cell_vec(1,3) + m2*recip_cell_vec(2,3) + m3*recip_cell_vec(3,3)
                if(g1*g1 + g2*g2 + g3*g3 < ewald_recip_cutoff_squared) then
                   number_of_g_vectors = number_of_g_vectors + 1
                endif
             endif
          enddo
       enddo
    enddo
    ! Now we have number_of_g_vectors
    ! Allocate ewald_g_vector, g_factor
    allocate(ewald_g_vector_x(number_of_g_vectors),ewald_g_vector_y(number_of_g_vectors), &
         ewald_g_vector_z(number_of_g_vectors),ewald_g_factor(number_of_g_vectors),STAT=stat)
    if(stat/=0) call cq_abort("set_ewald: error allocating ewald_g_vector ",number_of_g_vectors,stat)
    call reg_alloc_mem(area_general,4*number_of_g_vectors,type_dbl)
    ! Zero because we use it to count
    number_of_g_vectors = 0
    do m1 = -rec_lim_1, rec_lim_1
       do m2 = -rec_lim_2, rec_lim_2
          do m3 = -rec_lim_3, rec_lim_3
             if(m1*m1+m2*m2+m3*m3 /= 0) then
                g1 = m1*recip_cell_vec(1,1) + m2*recip_cell_vec(2,1) + m3*recip_cell_vec(3,1)
                g2 = m1*recip_cell_vec(1,2) + m2*recip_cell_vec(2,2) + m3*recip_cell_vec(3,2)
                g3 = m1*recip_cell_vec(1,3) + m2*recip_cell_vec(2,3) + m3*recip_cell_vec(3,3)
                if(g1*g1 + g2*g2 + g3*g3 < ewald_recip_cutoff_squared) then
                   number_of_g_vectors = number_of_g_vectors + 1
                   g_squared = g1*g1 + g2*g2 + g3*g3
                   ewald_g_vector_x(number_of_g_vectors) = g1
                   ewald_g_vector_y(number_of_g_vectors) = g2
                   ewald_g_vector_z(number_of_g_vectors) = g3
                   ewald_g_factor(number_of_g_vectors) = exp(-g_squared/(four*ewald_gamma))/g_squared
                endif
             endif
          enddo
       enddo
    enddo
    !    if(inode == ionode) then
    !       write(unit=io_lun,fmt='(/" no. of non-zero recips in Ewald sum:",i10)') number_of_g_vectors
    !       write(unit=io_lun,fmt='(//" recip vectors and recip factor:"/)')
    !       do n = 1, number_of_g_vectors
    !          write(unit=io_lun,fmt='(3x,i5,3x,3e15.6,3x,e15.6)') n, &
    !               &ewald_g_vector_x(n), ewald_g_vector_y(n), ewald_g_vector_z(n), ewald_g_factor(n)
    !       enddo
    !    endif

    ! *************************************************************************************************
    !                          SET-UP INFORMATION FOR EWALD REAL-SPACE SUM
    ! *************************************************************************************************

    ! +++ test printing of complementary error function
    if(inode == ionode.AND.iprint_gen>=6) then
       write(unit=io_lun,fmt='(/8x," +++ test printing of complementary error function calculated by 2 routines"/)')
       do n = 0, 100
          dummy1 = real(n,double) * 0.1_double
          dummy2 = obsolete_erfcc(dummy1)
          dummy3 = erfc(dummy1)
          write(unit=io_lun,fmt='(3x,i5,1x,f12.6,2x,2e20.12)') n, dummy1, dummy2, dummy3
       enddo
    endif
    ! +++

    ! --- primitive translation vectors for partition
    if(inode == ionode.AND.iprint_gen>1) write(unit=io_lun,fmt='(/8x," nos of partitions in cell edges:",3i5)') &
         &parts%ngcellx, parts%ngcelly, parts%ngcellz
    do i = 1, 3
       part_cell_vec(1,i) = real_cell_vec(1,i) / real(parts%ngcellx,double)
       part_cell_vec(2,i) = real_cell_vec(2,i) / real(parts%ngcelly,double)
       part_cell_vec(3,i) = real_cell_vec(3,i) / real(parts%ngcellz,double)
    enddo
    if(inode == ionode.AND.iprint_gen>1) then
       write(unit=io_lun,fmt='(/8x," partition cell vectors:"/)')
       do n = 1, 3
          write(unit=io_lun,fmt='(8x," partition cell vector no:",i5,3x,3f12.6)') n, (part_cell_vec(n,i), i = 1, 3)
       enddo
    endif

    ! --- diameter of sphere circumscribed round partition unit cell -
    partition_radius_squared = zero
    do m1 = 0, 1
       v1 = real(m1,double) - half
       do m2 = 0, 1
          v2 = real(m2,double) - half
          do m3 = 0, 1
             v3 = real(m3,double) - half
             vert(1) = v1*part_cell_vec(1,1) + v2*part_cell_vec(2,1) + v3*part_cell_vec(3,1)
             vert(2) = v1*part_cell_vec(1,2) + v2*part_cell_vec(2,2) + v3*part_cell_vec(3,2)
             vert(3) = v1*part_cell_vec(1,3) + v2*part_cell_vec(2,3) + v3*part_cell_vec(3,3)
             rsq = vert(1)*vert(1) + vert(2)*vert(2) + vert(3)*vert(3)
             if(rsq > partition_radius_squared) partition_radius_squared = rsq
          enddo
       enddo
    enddo
    partition_diameter = two*sqrt(partition_radius_squared)
    if(inode == ionode.AND.iprint_gen>1) write(unit=io_lun,fmt='(/8x," diameter of sphere circumscribed round &
         &partition cell:",f12.6)') partition_diameter 

    ! --- make Ewald covering set
    call make_cs(inode-1,ewald_real_cutoff,ewald_CS,parts,bundle,ni_in_cell, x_atom_cell,y_atom_cell,z_atom_cell)
    ! +++
    if(inode == ionode.AND.iprint_gen>1) then
       write(unit=io_lun,fmt='(/8x,"+++ ng_cover:",i10)') ewald_CS%ng_cover
       write(unit=io_lun,fmt='(8x,"+++ ncoverx, y, z:",3i8)') ewald_CS%ncoverx, ewald_CS%ncovery, ewald_CS%ncoverz
       write(unit=io_lun,fmt='(8x,"+++ nspanlx, y, z:",3i8)') ewald_CS%nspanlx, ewald_CS%nspanly, ewald_CS%nspanlz
       write(unit=io_lun,fmt='(8x,"+++ nx_origin, y, z:",3i8)') ewald_CS%nx_origin, ewald_CS%ny_origin, ewald_CS%nz_origin
    endif
    ! +++

    ! --- initiate partition_distance routine
    call partition_distance(.true.,part_cell_vec,part_cell_dual)

!    ! --- loop over primary-set partitions of present node
!    if(bundle%groups_on_node > mx_nponn) call cq_abort('mikes_set_ewald: mx_nponn to small',&
!         &bundle%groups_on_node)

    ! Work out n_ewald_partition_neighbours for first partition in group
    ip = 1

    ! --- Cartesian integer indices of current primary-set partition with respect
    !     to left of Ewald covering set (to be more precise: for a hypothetical primary-set 
    !     partition in bottom left-hand corner of covering set, all the following
    !     indices nppx, nppy, nppz would be equal to one).
    nppx = 1 + bundle%idisp_primx(ip) + ewald_CS%nspanlx
    nppy = 1 + bundle%idisp_primy(ip) + ewald_CS%nspanly
    nppz = 1 + bundle%idisp_primz(ip) + ewald_CS%nspanlz

    n_ewald_partition_neighbours = 0
    ! --- loop over partitions in Ewald covering set in node-periodic grouped order
    do nc = 1, ewald_CS%ng_cover
       ! --- for current partition nc in Ewald covering set, get the Cartesian-composite index
       ind_part = ewald_CS%lab_cover(nc)
       ! --- unpack the Cartesian-composite index to obtain the Cartesian integer indices of current
       !     partition nc in Ewald covering set
       nccx = 1 + (ind_part - 1)/(ewald_CS%ncovery * ewald_CS%ncoverz)
       nccy = 1 + (ind_part - 1 - (nccx-1) * ewald_CS%ncovery * ewald_CS%ncoverz)/ewald_CS%ncoverz
       nccz = ind_part - (nccx-1) * ewald_CS%ncovery * ewald_CS%ncoverz - (nccy-1) * ewald_CS%ncoverz
       ! --- integer triplet specifying offset of covering-set partition from primary-set partition
       m1 = nccx - nppx
       m2 = nccy - nppy
       m3 = nccz - nppz
       ! --- exclude from list if all components of offset are zero
       if(m1*m1 + m2*m2 + m3*m3 /= 0) then
          ! --- calculate distance (in a.u.) between centre of primary-set partition and the point
          !     mid-way between the primary- and covering-set partitions
          v1 = real(m1,double)
          v2 = real(m2,double)
          v3 = real(m3,double)
          rr(1) = half*(v1*part_cell_vec(1,1) + v2*part_cell_vec(2,1) + v3*part_cell_vec(3,1))
          rr(2) = half*(v1*part_cell_vec(1,2) + v2*part_cell_vec(2,2) + v3*part_cell_vec(3,2))
          rr(3) = half*(v1*part_cell_vec(1,3) + v2*part_cell_vec(2,3) + v3*part_cell_vec(3,3))
          call partition_distance(.false.,part_cell_vec,part_cell_dual,rr,distance)
          distance = two*distance
          if(distance < ewald_real_cutoff) then
             n_ewald_partition_neighbours = n_ewald_partition_neighbours + 1
          endif
       endif
    enddo
    ! Allocate memory
    allocate(ewald_partition_neighbour_list(bundle%groups_on_node,n_ewald_partition_neighbours), STAT=stat)
    if(stat/=0) call cq_abort("set_ewald: error allocating partition_neighbour_list ", &
         bundle%groups_on_node,n_ewald_partition_neighbours)
    call reg_alloc_mem(area_general,bundle%groups_on_node*n_ewald_partition_neighbours,type_dbl)

    do ip = 1, bundle%groups_on_node

       ! --- Cartesian integer indices of current primary-set partition with respect
       !     to left of Ewald covering set (to be more precise: for a hypothetical primary-set 
       !     partition in bottom left-hand corner of covering set, all the following
       !     indices nppx, nppy, nppz would be equal to one).
       nppx = 1 + bundle%idisp_primx(ip) + ewald_CS%nspanlx
       nppy = 1 + bundle%idisp_primy(ip) + ewald_CS%nspanly
       nppz = 1 + bundle%idisp_primz(ip) + ewald_CS%nspanlz

       n_ewald_partition_neighbours = 0
       ! --- loop over partitions in Ewald covering set in node-periodic grouped order
       do nc = 1, ewald_CS%ng_cover
          ! --- for current partition nc in Ewald covering set, get the Cartesian-composite index
          ind_part = ewald_CS%lab_cover(nc)
          ! --- unpack the Cartesian-composite index to obtain the Cartesian integer indices of current
          !     partition nc in Ewald covering set
          nccx = 1 + (ind_part - 1)/(ewald_CS%ncovery * ewald_CS%ncoverz)
          nccy = 1 + (ind_part - 1 - (nccx-1) * ewald_CS%ncovery * ewald_CS%ncoverz)/ewald_CS%ncoverz
          nccz = ind_part - (nccx-1) * ewald_CS%ncovery * ewald_CS%ncoverz - (nccy-1) * ewald_CS%ncoverz
          ! --- integer triplet specifying offset of covering-set partition from primary-set partition
          m1 = nccx - nppx
          m2 = nccy - nppy
          m3 = nccz - nppz
          ! --- exclude from list if all components of offset are zero
          if(m1*m1 + m2*m2 + m3*m3 /= 0) then
             ! --- calculate distance (in a.u.) between centre of primary-set partition and the point
             !     mid-way between the primary- and covering-set partitions
             v1 = real(m1,double)
             v2 = real(m2,double)
             v3 = real(m3,double)
             rr(1) = half*(v1*part_cell_vec(1,1) + v2*part_cell_vec(2,1) + v3*part_cell_vec(3,1))
             rr(2) = half*(v1*part_cell_vec(1,2) + v2*part_cell_vec(2,2) + v3*part_cell_vec(3,2))
             rr(3) = half*(v1*part_cell_vec(1,3) + v2*part_cell_vec(2,3) + v3*part_cell_vec(3,3))
             call partition_distance(.false.,part_cell_vec,part_cell_dual,rr,distance)
             distance = two*distance
             if(distance < ewald_real_cutoff) then
                n_ewald_partition_neighbours = n_ewald_partition_neighbours + 1
                ! --- for current primary-set partition ip, record the PG index nc of latest covering-set partition within cutoff
                ewald_partition_neighbour_list(ip,n_ewald_partition_neighbours) = nc
             endif
          endif
       enddo
       ! +++
       !if(inode == 2) then
       !   write(unit=io_lun,fmt='(" +++ for node no",i3," partition ",i3," no of Ewald neighbour partitions:",i5)') &
       !        &inode, ip, n_ewald_partition_neighbours
       !   write(unit=io_lun,fmt='(" no. of atoms in current partition:",i5)') bundle%nm_nodgroup(ip)
       !   if(bundle%nm_nodgroup(ip) > 0) then
       !      do ni = 1, bundle%nm_nodgroup(ip)
       !         i = bundle%ig_prim(bundle%nm_nodbeg(ip)+ni-1)
       !         write(unit=io_lun,fmt='(" +++ global index of atom no.",i5," in current partition:",i5)') ni, i
       !      enddo
       !   endif
       !endif
       ! +++
    enddo

    allocate(ewald_force(3,ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating ewald_force in set_ewald: ",ni_in_cell,stat)
    call reg_alloc_mem(area_general,3*ni_in_cell,type_dbl)
    call my_barrier()

    return

  end subroutine mikes_set_ewald
  ! =============================================================================

  ! ==============================================================================
  subroutine mikes_structure_factor(inode,ionode)

    use atoms, ONLY : index_my_atoms, n_my_atoms
    use GenComms, ONLY : gsum
    use global_module, ONLY : x_atom_cell, y_atom_cell, z_atom_cell, iprint_gen
    use numbers, ONLY : zero
    use primary_module, ONLY : bundle
    use species_module, ONLY : charge, species

    implicit none

    integer, intent(in) :: inode, ionode

    integer :: i, ip, k, n, ni
    real(double) :: g_dot_r, q_i, x, y, z


    ! >>>>>>
    !do ip = 1, bundle%groups_on_node
    !   if(inode /= ionode.AND.iprint_gen>=6) write(unit=io_lun,fmt='(/" node:",i3," primary-set partition no:",i3," contains",i3,&
    !        &" atoms")') inode, ip, bundle%nm_nodgroup(ip)
    !   if(bundle%nm_nodgroup(ip) > 0) then
    !      do ni = 1, bundle%nm_nodgroup(ip)
    !         x = bundle%xprim(bundle%nm_nodbeg(ip)+ni-1)
    !         y = bundle%yprim(bundle%nm_nodbeg(ip)+ni-1)
    !         z = bundle%zprim(bundle%nm_nodbeg(ip)+ni-1)
    !         i = bundle%ig_prim(bundle%nm_nodbeg(ip)+ni-1)
    !         !q_i = charge(species(i))
    !         q_i = charge(bundle%species(bundle%nm_nodbeg(ip)+ni-1))
    !         if(inode /= ionode.AND.iprint_gen>=6) write(unit=io_lun,fmt='(2x,2i3,2x,3e15.6,2x,f12.6)') ni, i, x, y, z, q_i
    !      enddo
    !   endif
    !enddo
    ! <<<<<<

    do n = 1, number_of_g_vectors
       struc_fac_r(n) = zero
       struc_fac_i(n) = zero
       do ip = 1, bundle%groups_on_node
          if(bundle%nm_nodgroup(ip) > 0) then
             do ni = 1, bundle%nm_nodgroup(ip)
                x = bundle%xprim(bundle%nm_nodbeg(ip)+ni-1)
                y = bundle%yprim(bundle%nm_nodbeg(ip)+ni-1)
                z = bundle%zprim(bundle%nm_nodbeg(ip)+ni-1)
                i = bundle%ig_prim(bundle%nm_nodbeg(ip)+ni-1)
                !q_i = charge(species(i))
                q_i = charge(bundle%species(bundle%nm_nodbeg(ip)+ni-1))
                g_dot_r = ewald_g_vector_x(n)*x + ewald_g_vector_y(n)*y + ewald_g_vector_z(n)*z
                struc_fac_r(n) = struc_fac_r(n) + q_i * cos(g_dot_r)
                struc_fac_i(n) = struc_fac_i(n) + q_i * sin(g_dot_r)
             enddo
          endif
       enddo
    enddo

    call gsum(struc_fac_r,number_of_g_vectors)
    call gsum(struc_fac_i,number_of_g_vectors)

    ! +++ temporary printing
    if(inode == ionode.AND.iprint_gen>=6) then
       write(unit=io_lun,fmt='(/"+++ structure factors:"/)')
       do n = 1, number_of_g_vectors
          write(unit=io_lun,fmt='(i10,2x,3f12.6,2x,2e15.6)') n, &
               &ewald_g_vector_x(n), ewald_g_vector_y(n), ewald_g_vector_z(n), &
               &struc_fac_r(n), struc_fac_i(n)
       enddo
    endif
    ! +++

  end subroutine mikes_structure_factor
  ! ==============================================================================

!!   10:42, 13/02/2006 drb 
!!    Corrected line length error
  ! =============================================================================
  subroutine mikes_ewald()

    use cover_module, ONLY : ewald_CS
    use datatypes
    use GenComms, ONLY : gsum, cq_abort, inode, ionode
    use global_module, ONLY : id_glob, iprint_gen, species_glob, ni_in_cell, area_general, IPRINT_TIME_THRES3
    use group_module, ONLY : parts
    use maxima_module, ONLY : maxatomsproc
    use numbers
    use primary_module, ONLY : bundle
    use species_module, ONLY : charge, species
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module, ONLY: cq_timer,start_timer,stop_print_timer,WITH_LEVEL

    implicit none

    integer :: i, ip, ig_atom_beg, j, n, nc, ni, nj, nn, stat
    real(double) :: arg_1, arg_2, coeff_1, dummy, ewald_real_energy, ewald_real_sum_inter, ewald_real_sum_intra, &
         &ewald_real_sum_ip, ewald_recip_energy, ewald_total_energy, g_dot_r, q_i, q_j, rij, rijx, rijy, rijz, rij_squared, x, y, z
    real(double), allocatable, dimension(:) :: ewald_recip_force_x, ewald_recip_force_y, ewald_recip_force_z, &
         &ewald_intra_force_x, ewald_intra_force_y, ewald_intra_force_z, &
         &ewald_inter_force_x, ewald_inter_force_y, ewald_inter_force_z
    real(double), allocatable, dimension(:) :: ewald_total_force_x, ewald_total_force_y, ewald_total_force_z
    type(cq_timer) :: tmr_l_tmp1

    ! --- Ewald reciprocal-space energy and forces  --------------------------------
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    allocate(ewald_recip_force_x(maxatomsproc),ewald_recip_force_y(maxatomsproc), &
         ewald_recip_force_z(maxatomsproc), ewald_intra_force_x(maxatomsproc), &
         ewald_intra_force_y(maxatomsproc), ewald_intra_force_z(maxatomsproc), &
         &ewald_inter_force_x(maxatomsproc), ewald_inter_force_y(maxatomsproc), &
         ewald_inter_force_z(maxatomsproc),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating ewald forces(1): ",maxatomsproc,stat)
    call reg_alloc_mem(area_general,9*maxatomsproc,type_dbl)
    allocate(ewald_total_force_x(ni_in_cell),ewald_total_force_y(ni_in_cell), &
         ewald_total_force_z(ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating ewald forces(2): ",ni_in_cell,stat)
    call reg_alloc_mem(area_general,3*ni_in_cell,type_dbl)
    ! ... initialise energy and forces to zero
    ewald_recip_energy = zero
    do ip = 1, bundle%groups_on_node
       if(bundle%nm_nodgroup(ip) > 0) then
          do ni = 1, bundle%nm_nodgroup(ip)
             ewald_recip_force_x(bundle%nm_nodbeg(ip)+ni-1) = zero
             ewald_recip_force_y(bundle%nm_nodbeg(ip)+ni-1) = zero
             ewald_recip_force_z(bundle%nm_nodbeg(ip)+ni-1) = zero
          enddo
       endif
    enddo
    ! Allocate memory for and call structure_factor
    allocate(struc_fac_r(number_of_g_vectors),struc_fac_i(number_of_g_vectors),STAT=stat)
    if(stat/=0) call cq_abort("ewald: error allocating struc_fac ",number_of_g_vectors,stat)
    call reg_alloc_mem(area_general,2*number_of_g_vectors,type_dbl)
    call mikes_structure_factor(inode,ionode)
    ! ... loop over reciprocal lattice vectors
    do n = 1, number_of_g_vectors
       ewald_recip_energy = ewald_recip_energy + &
            &(struc_fac_r(n)*struc_fac_r(n) + struc_fac_i(n)*struc_fac_i(n)) * &
            &ewald_g_factor(n)
       do ip = 1, bundle%groups_on_node
          if(bundle%nm_nodgroup(ip) > 0) then
             do ni = 1, bundle%nm_nodgroup(ip)
                i = bundle%ig_prim(bundle%nm_nodbeg(ip)+ni-1)
                !q_i = charge(species(i))
                q_i = charge(bundle%species(bundle%nm_nodbeg(ip)+ni-1))
                x = bundle%xprim(bundle%nm_nodbeg(ip)+ni-1)
                y = bundle%yprim(bundle%nm_nodbeg(ip)+ni-1)
                z = bundle%zprim(bundle%nm_nodbeg(ip)+ni-1)
                g_dot_r = ewald_g_vector_x(n)*x + ewald_g_vector_y(n)*y + ewald_g_vector_z(n)*z
                ewald_recip_force_x(bundle%nm_nodbeg(ip)+ni-1) = ewald_recip_force_x(bundle%nm_nodbeg(ip)+ni-1) + &
                     &q_i * ewald_g_vector_x(n) * ewald_g_factor(n) * (cos(g_dot_r)*struc_fac_i(n) - sin(g_dot_r)*struc_fac_r(n))
                ewald_recip_force_y(bundle%nm_nodbeg(ip)+ni-1) = ewald_recip_force_y(bundle%nm_nodbeg(ip)+ni-1) + &
                     &q_i * ewald_g_vector_y(n) * ewald_g_factor(n) * (cos(g_dot_r)*struc_fac_i(n) - sin(g_dot_r)*struc_fac_r(n))
                ewald_recip_force_z(bundle%nm_nodbeg(ip)+ni-1) = ewald_recip_force_z(bundle%nm_nodbeg(ip)+ni-1) + &
                     &q_i * ewald_g_vector_z(n) * ewald_g_factor(n) * (cos(g_dot_r)*struc_fac_i(n) - sin(g_dot_r)*struc_fac_r(n))
             enddo
          endif
       enddo
    enddo
    deallocate(struc_fac_r,struc_fac_i,STAT=stat)
    if(stat/=0) call cq_abort("ewald: error deallocating struc_fac",stat)
    call reg_dealloc_mem(area_general,2*number_of_g_vectors,type_dbl)
    ewald_recip_energy = two*pi*ewald_recip_energy/ewald_real_cell_volume
    do ip = 1, bundle%groups_on_node
       if(bundle%nm_nodgroup(ip) > 0) then
          do ni = 1, bundle%nm_nodgroup(ip)
             ewald_recip_force_x(bundle%nm_nodbeg(ip)+ni-1) = &
                  &-(four*pi/ewald_real_cell_volume)*ewald_recip_force_x(bundle%nm_nodbeg(ip)+ni-1)
             ewald_recip_force_y(bundle%nm_nodbeg(ip)+ni-1) = &
                  &-(four*pi/ewald_real_cell_volume)*ewald_recip_force_y(bundle%nm_nodbeg(ip)+ni-1)
             ewald_recip_force_z(bundle%nm_nodbeg(ip)+ni-1) = &
                  &-(four*pi/ewald_real_cell_volume)*ewald_recip_force_z(bundle%nm_nodbeg(ip)+ni-1)
          enddo
       endif
    enddo

    !if(inode /= ionode.AND.iprint_gen>=6) then
    !   write(unit=io_lun,fmt='(/" Ewald recip-space info from node:",i3)') inode
    !   write(unit=io_lun,fmt='(/" Ewald recip-space energy:",e20.12)') &
    !        &ewald_recip_energy
    !   write(unit=io_lun,fmt='(/" Ewald recip-space forces:"/)')
    !   do ip = 1, bundle%groups_on_node
    !      if(bundle%nm_nodgroup(ip) > 0) then
    !         do ni = 1, bundle%nm_nodgroup(ip)
    !            i = bundle%ig_prim(bundle%nm_nodbeg(ip)+ni-1)
    !            write(unit=io_lun,fmt='(2x,3i5,2x,3e20.12)') ip, ni, i, &
    !                 &ewald_recip_force_x(bundle%nm_nodbeg(ip)+ni-1), ewald_recip_force_y(bundle%nm_nodbeg(ip)+ni-1), &
    !                 &ewald_recip_force_z(bundle%nm_nodbeg(ip)+ni-1)
    !         enddo
    !      endif
    !   enddo
    !endif

    ! --- Ewald real-space sum -----------------------------------------------------

    ! --- intra-partition contribution ---------------------------------------------
    coeff_1 = two*sqrt(ewald_gamma/pi)
    ! ... initialise energy and forces to zero
    ewald_real_sum_intra = zero
    do ip = 1, bundle%groups_on_node
       if(bundle%nm_nodgroup(ip) > 0) then
          do ni = 1, bundle%nm_nodgroup(ip)
             ewald_intra_force_x(bundle%nm_nodbeg(ip)+ni-1) = zero
             ewald_intra_force_y(bundle%nm_nodbeg(ip)+ni-1) = zero
             ewald_intra_force_z(bundle%nm_nodbeg(ip)+ni-1) = zero
          enddo
       endif
    enddo

    ! --- loop over primary-set partitions
    do ip = 1, bundle%groups_on_node
       ! --- loop over distinct pairs of atoms in current primary-set partition
       ewald_real_sum_ip = zero
       if(bundle%nm_nodgroup(ip) > 1) then
          do ni = 1, bundle%nm_nodgroup(ip) - 1
             ! ... i is global label of atom, q_i is its charge
             i = bundle%ig_prim(bundle%nm_nodbeg(ip)+ni-1)
             !q_i = charge(species(i))
             q_i = charge(bundle%species(bundle%nm_nodbeg(ip)+ni-1))
             do nj = ni+1, bundle%nm_nodgroup(ip)
                ! ... global label and charge, as before
                j = bundle%ig_prim(bundle%nm_nodbeg(ip)+nj-1)
                !q_j = charge(species(j))
                q_j = charge(bundle%species(bundle%nm_nodbeg(ip)+nj-1))
                ! ... Cartesian components of vector separation of the two atoms
                rijx = bundle%xprim(bundle%nm_nodbeg(ip)+ni-1) - bundle%xprim(bundle%nm_nodbeg(ip)+nj-1)
                rijy = bundle%yprim(bundle%nm_nodbeg(ip)+ni-1) - bundle%yprim(bundle%nm_nodbeg(ip)+nj-1)
                rijz = bundle%zprim(bundle%nm_nodbeg(ip)+ni-1) - bundle%zprim(bundle%nm_nodbeg(ip)+nj-1)
                rij_squared = rijx*rijx + rijy*rijy + rijz*rijz
                rij = sqrt(rij_squared)
                if(rij < ewald_real_cutoff) then
                   arg_1 = ewald_gamma*rij_squared
                   arg_2 = sqrt(arg_1)
                   ewald_real_sum_ip = ewald_real_sum_ip + &
                        &q_i * q_j * erfc(arg_2) / rij
                   dummy = q_i * q_j * (erfc(arg_2)/rij + coeff_1*exp(-arg_1))/rij_squared
                   ewald_intra_force_x(bundle%nm_nodbeg(ip)+ni-1) = ewald_intra_force_x(bundle%nm_nodbeg(ip)+ni-1) + dummy*rijx
                   ewald_intra_force_x(bundle%nm_nodbeg(ip)+nj-1) = ewald_intra_force_x(bundle%nm_nodbeg(ip)+nj-1) - dummy*rijx
                   ewald_intra_force_y(bundle%nm_nodbeg(ip)+ni-1) = ewald_intra_force_y(bundle%nm_nodbeg(ip)+ni-1) + dummy*rijy
                   ewald_intra_force_y(bundle%nm_nodbeg(ip)+nj-1) = ewald_intra_force_y(bundle%nm_nodbeg(ip)+nj-1) - dummy*rijy
                   ewald_intra_force_z(bundle%nm_nodbeg(ip)+ni-1) = ewald_intra_force_z(bundle%nm_nodbeg(ip)+ni-1) + dummy*rijz
                   ewald_intra_force_z(bundle%nm_nodbeg(ip)+nj-1) = ewald_intra_force_z(bundle%nm_nodbeg(ip)+nj-1) - dummy*rijz
                endif
             enddo
          enddo
       endif
       if(inode == ionode.AND.iprint_gen>=4) &
            write(unit=io_lun,fmt='(/" >>> mikes_ewald: node:",i3:" real_sum_intra for partition no:",i3,&
            &" is:",e20.12)') inode, ip, ewald_real_sum_ip
       ewald_real_sum_intra = ewald_real_sum_intra + ewald_real_sum_ip
    enddo

    call gsum(ewald_real_sum_intra)

    if(inode == ionode.AND.iprint_gen>=6) then
       write(unit=io_lun,fmt='(/" Ewald real-space self info for node:",i3/)') inode 
       write(unit=io_lun,fmt='(/" self-partition part of real_space Ewald:",e20.12)') ewald_real_sum_intra
       write(unit=io_lun,fmt='(/" self-partition part of real-space Ewald forces:"/)')
       do ip = 1, bundle%groups_on_node
          if(bundle%nm_nodgroup(ip) > 0) then
             do ni = 1, bundle%nm_nodgroup(ip)
                i = bundle%ig_prim(bundle%nm_nodbeg(ip)+ni-1)
                write(unit=io_lun,fmt='(2x,3i5,2x,3e20.12)') ip, ni, i, &
                     &ewald_intra_force_x(bundle%nm_nodbeg(ip)+ni-1), ewald_intra_force_y(bundle%nm_nodbeg(ip)+ni-1), &
                     &ewald_intra_force_z(bundle%nm_nodbeg(ip)+ni-1)
             enddo
          endif
       enddo
    endif

    ! --- inter-partition contribution ---------------------------------------------

    ! ... initialise energy and forces to zero
    ewald_real_sum_inter = zero
    do ip = 1, bundle%groups_on_node
       if(bundle%nm_nodgroup(ip) > 0) then
          do ni = 1, bundle%nm_nodgroup(ip)
             ewald_inter_force_x(bundle%nm_nodbeg(ip)+ni-1) = zero
             ewald_inter_force_y(bundle%nm_nodbeg(ip)+ni-1) = zero
             ewald_inter_force_z(bundle%nm_nodbeg(ip)+ni-1) = zero
          enddo
       endif
    enddo

    ! --- loop over primary-set partitions
    do ip = 1, bundle%groups_on_node
       ! --- loop over atoms in current primary-set partition
       ewald_real_sum_ip = zero
       if(bundle%nm_nodgroup(ip) > 0) then
          do ni = 1, bundle%nm_nodgroup(ip)
             ! ... global label and charge of current primary-set atom
             i = bundle%ig_prim(bundle%nm_nodbeg(ip)+ni-1)
             !q_i = charge(species(i))
             q_i = charge(bundle%species(bundle%nm_nodbeg(ip)+ni-1))
             ! --- loop over Ewald neighbour-list of partitions in Ewald covering set
             do nn = 1, n_ewald_partition_neighbours
                ! ... get PG label of current ECS partition
                nc = ewald_partition_neighbour_list(ip,nn)
                ig_atom_beg = parts%icell_beg(ewald_CS%lab_cell(nc))
                ! --- loop over atoms in ECS partition
                if(ewald_CS%n_ing_cover(nc) > 0) then
                   do nj = 1, ewald_CS%n_ing_cover(nc)
                      ! ... global label and charge of current atom in ECS partition
                      !q_j = charge(species(id_glob(ig_atom_beg+nj-1)))
                      q_j = charge(species_glob(id_glob(ig_atom_beg+nj-1)))
                      ! ... Cartesian components of vector separation of the two atoms
                      rijx = bundle%xprim(bundle%nm_nodbeg(ip)+ni-1) - ewald_CS%xcover(ewald_CS%icover_ibeg(nc)+nj-1)
                      rijy = bundle%yprim(bundle%nm_nodbeg(ip)+ni-1) - ewald_CS%ycover(ewald_CS%icover_ibeg(nc)+nj-1)
                      rijz = bundle%zprim(bundle%nm_nodbeg(ip)+ni-1) - ewald_CS%zcover(ewald_CS%icover_ibeg(nc)+nj-1)
                      rij_squared = rijx*rijx + rijy*rijy + rijz*rijz
                      rij = sqrt(rij_squared)
                      if(rij < ewald_real_cutoff) then
                         arg_1 = ewald_gamma*rij_squared
                         arg_2 = sqrt(arg_1)
                         ewald_real_sum_ip = ewald_real_sum_ip + &
                              &q_i * q_j * erfc(arg_2) / rij
                         dummy = q_i * q_j * (erfc(arg_2)/rij + coeff_1*exp(-arg_1))/rij_squared
                         ewald_inter_force_x(bundle%nm_nodbeg(ip)+ni-1) = &
                              ewald_inter_force_x(bundle%nm_nodbeg(ip)+ni-1) + dummy*rijx
                         ewald_inter_force_y(bundle%nm_nodbeg(ip)+ni-1) = &
                              ewald_inter_force_y(bundle%nm_nodbeg(ip)+ni-1) + dummy*rijy
                         ewald_inter_force_z(bundle%nm_nodbeg(ip)+ni-1) = &
                              ewald_inter_force_z(bundle%nm_nodbeg(ip)+ni-1) + dummy*rijz
                      endif
                   enddo
                endif
             enddo
          enddo
       endif
       if(inode == ionode.AND.iprint_gen>=4) &
            write(unit=io_lun,fmt='(/" >>> mikes_ewald: node:",i3:" real_sum_inter for partition no:",i3,&
            &" is:",e20.12)') inode, ip, ewald_real_sum_ip
       ewald_real_sum_inter = ewald_real_sum_inter + ewald_real_sum_ip
    enddo

    call gsum(ewald_real_sum_inter)

    if(inode == ionode.AND.iprint_gen>=6) then
       write(unit=io_lun,fmt='(/" Ewald real-space other info for node:",i3/)') inode 
       write(unit=io_lun,fmt='(/" other-partition part of real_space Ewald:",e20.12)') ewald_real_sum_inter
       write(unit=io_lun,fmt='(/" other-partition part of real-space Ewald forces:"/)')
       do ip = 1, bundle%groups_on_node
          if(bundle%nm_nodgroup(ip) > 0) then
             do ni = 1, bundle%nm_nodgroup(ip)
                i = bundle%ig_prim(bundle%nm_nodbeg(ip)+ni-1)
                write(unit=io_lun,fmt='(2x,3i5,2x,3e20.12)') ip, ni, i, &
                     &ewald_inter_force_x(bundle%nm_nodbeg(ip)+ni-1), ewald_inter_force_y(bundle%nm_nodbeg(ip)+ni-1), &
                     &ewald_inter_force_z(bundle%nm_nodbeg(ip)+ni-1)
             enddo
          endif
       enddo
    endif

    ewald_real_energy = ewald_real_sum_intra + half*ewald_real_sum_inter
    if(inode == ionode.AND.iprint_gen>0) write(unit=io_lun,fmt='(/8x," ++++++ Ewald real-space energy:",e20.12)') &
         &ewald_real_energy

    ! --- add all contributions to total Ewald energy and forces
    ewald_total_energy = ewald_recip_energy + ewald_real_energy + ewald_gaussian_self_energy + ewald_net_charge_energy
    ewald_energy = ewald_total_energy
    do i = 1, ni_in_cell
       ewald_total_force_x(i) = zero
       ewald_total_force_y(i) = zero
       ewald_total_force_z(i) = zero
    enddo
    do ip = 1, bundle%groups_on_node
       if(bundle%nm_nodgroup(ip) > 0) then
          do ni = 1, bundle%nm_nodgroup(ip)
             i = bundle%ig_prim(bundle%nm_nodbeg(ip)+ni-1)
             ewald_total_force_x(i) = ewald_recip_force_x(bundle%nm_nodbeg(ip)+ni-1) + &
                  &ewald_intra_force_x(bundle%nm_nodbeg(ip)+ni-1) + ewald_inter_force_x(bundle%nm_nodbeg(ip)+ni-1)
             ewald_total_force_y(i) = ewald_recip_force_y(bundle%nm_nodbeg(ip)+ni-1) + &
                  &ewald_intra_force_y(bundle%nm_nodbeg(ip)+ni-1) + ewald_inter_force_y(bundle%nm_nodbeg(ip)+ni-1)
             ewald_total_force_z(i) = ewald_recip_force_z(bundle%nm_nodbeg(ip)+ni-1) + &
                  &ewald_intra_force_z(bundle%nm_nodbeg(ip)+ni-1) + ewald_inter_force_z(bundle%nm_nodbeg(ip)+ni-1)
          enddo
       endif
    enddo

    call gsum(ewald_total_force_x,ni_in_cell)
    call gsum(ewald_total_force_y,ni_in_cell)
    call gsum(ewald_total_force_z,ni_in_cell)
    ! Added DRB & MJG 2003/10/22 to make force available to rest of code
    do i=1,ni_in_cell
       ewald_force(1,i) = ewald_total_force_x(i)
       ewald_force(2,i) = ewald_total_force_y(i)
       ewald_force(3,i) = ewald_total_force_z(i)
    end do
    if(inode == ionode.AND.iprint_gen>=2) then
       !write(unit=io_lun,fmt='(/" Ewald total energy and forces from node:",i3)') inode 
       write(unit=io_lun,fmt='(/8x," Ewald total energy:",e20.12)') ewald_total_energy
       if(inode == ionode.AND.iprint_gen>=4) then
          write(unit=io_lun,fmt='(/8x," Ewald total forces:"/)')
          do i = 1, ni_in_cell
             write(unit=io_lun,fmt='(10x,i5,2x,3e20.12)') i, ewald_total_force_x(i), ewald_total_force_y(i), ewald_total_force_z(i)
          enddo
       end if
    endif
    deallocate(ewald_total_force_z,ewald_total_force_y, ewald_total_force_x,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating ewald forces(2): ",ni_in_cell,stat)
    call reg_dealloc_mem(area_general,3*ni_in_cell,type_dbl)
    deallocate(ewald_recip_force_x,ewald_recip_force_y, ewald_recip_force_z, ewald_intra_force_x, &
         ewald_intra_force_y, ewald_intra_force_z, ewald_inter_force_x, ewald_inter_force_y, &
         ewald_inter_force_z,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating ewald forces: ",maxatomsproc)
    call reg_dealloc_mem(area_general,9*maxatomsproc,type_dbl)
    call stop_print_timer(tmr_l_tmp1,"mikes_ewald",IPRINT_TIME_THRES3)

  end subroutine mikes_ewald
  ! =============================================================================

end module ewald_module
