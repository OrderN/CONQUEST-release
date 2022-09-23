! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module ion_electrostatic
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/ion_electrostatic *
!!  NAME
!!   ion_electrostatic
!!  PURPOSE
!!   Collects variables and subroutines related to electrostatic interactions between ions
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
!!   2015/05/01 13:24 dave and sym
!!    Adding stress
!!   2015/11/24 08:21 dave with TM and N. Watanabe of Mizuho
!!    Introducing neutral atom variables and routines, and removing old ewald
!!   2015/11/26 15:20 dave
!!    Changed standard (ewald or other method) names to ion_interaction_ instead of ewald_
!!   2015/12/21 17:25 dave
!!    Changed n_ewald_partition_neighbours to n_partition_neighbours (use for both ewald and neutral atom)
!!    and ewald_partition_list
!!   2016/01/28 16:35 dave
!!    Tidying and simplifying the neutral-atom screened ion interactions, and changing module name
!!    from ewald_module to ion_electrostatic (file name will be changed to ion_electrostatic_module.f90)
!!   2016/02/09 08:10 dave
!!    Moved erfc to functions module; also energies to energy.module
!!   2019/04/09 zamaan
!!    Off-diagonal stress tensor contributions
!!  SOURCE
!!
module ion_electrostatic

  use datatypes
  use global_module, ONLY: io_lun
  use GenComms,      ONLY: cq_abort

  implicit none
  save

  integer :: n_g_vectors, n_supercells
  integer :: number_of_g_vectors, n_partition_neighbours
  real(double) :: gamma
  real(double) :: ewald_gamma, ewald_real_cutoff, ewald_recip_cutoff
  real(double) :: ewald_real_cell_volume
  real(double) :: ewald_gaussian_self_energy, ewald_net_charge_energy
  real(double), allocatable, dimension(:) :: gx, gy, gz, g_factor, structure2
  complex(double_cplx), allocatable, dimension(:) :: cstructure
  real(double), allocatable, dimension(:) :: supercell_vec_x, supercell_vec_y, supercell_vec_z
  integer, dimension(:,:), allocatable :: partition_neighbour_list
  real(double), dimension(:), allocatable :: ewald_g_vector_x, ewald_g_vector_y, ewald_g_vector_z, &
       &ewald_g_factor, struc_fac_r, struc_fac_i


  ! Standard pseudopotential (not neutral atom)
  real(double), allocatable, dimension(:,:) :: ion_interaction_force
  real(double), dimension(3,3) :: ion_interaction_stress
  real(double), dimension(3,3) :: ewald_gaussian_self_stress, &
                                  ewald_recip_stress, ewald_intra_stress, &
                                  ewald_inter_stress
  real(double) :: ewald_accuracy, ion_ion_cutoff

  ! For neutral atom potential
  real(double), allocatable, dimension(:,:) :: screened_ion_force
  real(double), dimension(3,3) :: screened_ion_stress  
  
  !!***

contains

!!****f* ion_electrostatic/partition_distance *
!!
!!  NAME 
!!   partition_distance
!!  USAGE
!!   
!!  PURPOSE
!!   Finds the distance between a primary set partition and a covering set partition
!!   (though for reasons that are not clear, we actually find half this distance)
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   M. J. Gillan
!!  CREATION DATE
!!   2003
!!  MODIFICATION HISTORY
!!   2005/10/06 10:13 dave
!!    It was assumed that a1_dot_a1 etc and b1_dot_b1 etc would remain defined after exit
!!    from routine; I added an explicit recalculation if pd_init .false.
!!  SOURCE
!!  
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
       if(distance > RD_ERR) then 
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

    return
  end subroutine partition_distance
!!***

!!****f* ion_electrostatic/set_ewald *
!!
!!  NAME 
!!   set_ewald
!!  USAGE
!!   
!!  PURPOSE
!!   Sets up all the arrays and routines needed for Ewald sum calculations
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   M. J. Gillan
!!  CREATION DATE
!!   2003
!!  MODIFICATION HISTORY
!!   2009/07/08 16:58 dave
!!    Removed double loop over all atoms to find lambda as it's not used
!!   2015/05/01 13:26 dave and sym
!!    Adding stress calculation  
!!   2015/11/24 08:28 dave
!!    Removing call to obsolete_erfcc (only in a redundant test loop)
!!   2016/01/29 14:57 dave
!!    Name change: only one ewald routine left, so removed prefix mikes_
!!  SOURCE
!!  
  subroutine set_ewald(inode,ionode)

    use cover_module, ONLY : ion_ion_CS, make_cs
    use datatypes
    use dimens
    use GenComms, ONLY : cq_abort, my_barrier
    use group_module, ONLY : parts
    use global_module, ONLY : x_atom_cell, y_atom_cell, z_atom_cell, &
                              iprint_gen, ni_in_cell, area_general, flag_stress
    use numbers
    use primary_module, ONLY : bundle
    use species_module, ONLY : charge, species
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl


    implicit none

    integer, intent(in) :: inode, ionode

    integer :: i, ind_part, ip, m_limit_1, m_limit_2, m_limit_3, m1, &
         m2, m3, n, nc, nccx, nccy, nccz, ni, np, nppx, nppy, nppz, &
         np1, np2, rec_lim_1, rec_lim_2, rec_lim_3, run_mode, stat, direction
    real(double) :: ewald_recip_cutoff_squared, mean_number_density, &
         mean_interparticle_distance, mean_square_charge, &
         partition_radius_squared, partition_diameter, total_charge
    real(double) :: argument, const_phi_0, const_phi_1, const_zeta_0, &
         const_zeta_1, distance, dummy1, dummy2, dummy3, g1, g2, g3, &
         g_squared, lambda, rsq, sep, vol_prefac, v1, v2, v3
    real(double) :: gaussian_self_stress
    real(double), dimension(3) :: rr, vert
!    real(double), dimension(3,3) :: abs_sep, rel_sep
! 16:52, 2003/10/22 dave Changed to allow compilation under ifc
    real(double), dimension(3) :: abs_sep, rel_sep
    real(double), dimension(3,3) :: real_cell_vec, recip_cell_vec, &
         part_cell_vec, part_cell_dual

    if(inode == ionode.AND.iprint_gen>1) &
         write(io_lun,fmt='(/8x," edge lengths of orthorhombic supercell:", &
         &3f12.6," a.u.")') r_super_x, r_super_y, r_super_z

    ! ------------ setting primitive translation vectors of cell ------------
    !   At time of writing, Conquest is restricted to orthorhombic
    !   cells. However the present Ewald routine is written for
    !   general cells. The variables r_super_x, r_super_y, r_super_z
    !   are the edge lengths of the orthorhombic cell. These are used
    !   to make the general real-space primitive translation vectors,
    !   which in turn are used to make the reciprocal-space
    !   vectors. Here, real_cell_vec(n,i) and recip_cell_vec(n,i) are
    !   the Cartesian components (labelled by i) of the real- and
    !   reciprocal-space vectors (labelled by n).
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
       write(io_lun,fmt='(/8x," set_ewald: cartesian components &
            &of real-cell primitive translations vectors:"/)')
       do n = 1, 3
          write(io_lun,fmt='(8x," vector no:",i3,":",3f12.6)') &
               n, (real_cell_vec(n,i), i = 1, 3)
       end do
    end if

    ! --- volume of real-space unit cell ------------------------------
    ewald_real_cell_volume = &
         real_cell_vec(1,1)*(real_cell_vec(2,2)*real_cell_vec(3,3) - &
         real_cell_vec(2,3)*real_cell_vec(3,2)) + &
         real_cell_vec(1,2)*(real_cell_vec(2,3)*real_cell_vec(3,1) - &
         real_cell_vec(2,1)*real_cell_vec(3,3)) + &
         real_cell_vec(1,3)*(real_cell_vec(2,1)*real_cell_vec(3,2) - &
         real_cell_vec(2,2)*real_cell_vec(3,1))
    ewald_real_cell_volume = abs(ewald_real_cell_volume)
    if(inode==ionode.AND.iprint_gen>1) &
         write(io_lun,fmt='(/8x," volume of real-space cell:",f19.6)')&
         ewald_real_cell_volume

    ! --- quantities needed for calculation of Ewald gamma, and real_space and recip-space cutoffs
    total_charge = zero
    mean_square_charge = zero
    do np = 1, ni_in_cell
       total_charge = total_charge + charge(species(np))
       mean_square_charge = mean_square_charge + &
            charge(species(np))*charge(species(np))
    enddo
    mean_square_charge = mean_square_charge/real(ni_in_cell,double)
    mean_number_density = real(ni_in_cell,double)/ewald_real_cell_volume
    mean_interparticle_distance = (three/(four*pi*mean_number_density))**one_third
    ! Rescale ewald_accuracy to a dimensionless number
    ewald_accuracy = ewald_accuracy*mean_interparticle_distance/mean_square_charge
    ! --- miscellaneous constants used for making Ewald gamma and
    ! the real-space and recip-space cutoffs
    const_zeta_0 = (pi**five_sixths)/(six**one_third)
    const_zeta_1 = (two**two_thirds)*(pi**five_sixths)/(three**one_third)
    argument = min(exp(-one),&
         &(const_zeta_0*ewald_accuracy)/(real(ni_in_cell,double)**one_third))
    const_phi_0 = sqrt(-log(argument))
    argument = min(exp(-one),&
         &(const_zeta_1*ewald_accuracy)/(real(ni_in_cell,double)**five_sixths))
    const_phi_1 = sqrt(-log(argument))
    if(inode == ionode.AND.iprint_gen>1) then
       write (unit=io_lun, &
            fmt='(/8x," const_zeta_0:",e20.12," const_zeta_1:",e20.12/&
            &8x," const_phi_0:",e20.12," const_phi_1:",e20.12/&
            &8x," total charge:",e20.12," mean-square charge:",e20.12/&
            &8x," mean number density:",e20.12," mean interparticle &
            &distance:",e20.12)') &
            const_zeta_0, const_zeta_1, const_phi_0, const_phi_1, &
            total_charge, mean_square_charge, mean_number_density, &
            mean_interparticle_distance
    endif

    ! *************************************************************************
    !              SET-UP INFORMATION FOR EWALD RECIPROCAL-SPACE SUM
    ! *************************************************************************

    ! --- reciprocal-lattice vectors -------------------
    vol_prefac = 2.0_double*pi/ewald_real_cell_volume  

    recip_cell_vec(1,1) = vol_prefac * (real_cell_vec(2,2) * &
         real_cell_vec(3,3) - real_cell_vec(2,3)*real_cell_vec(3,2))
    recip_cell_vec(1,2) = vol_prefac * (real_cell_vec(2,3) * &
         real_cell_vec(3,1) - real_cell_vec(2,1)*real_cell_vec(3,3))
    recip_cell_vec(1,3) = vol_prefac * (real_cell_vec(2,1) * &
         real_cell_vec(3,2) - real_cell_vec(2,2)*real_cell_vec(3,1))

    recip_cell_vec(2,1) = vol_prefac * (real_cell_vec(3,2) * &
         real_cell_vec(1,3) - real_cell_vec(3,3)*real_cell_vec(1,2))
    recip_cell_vec(2,2) = vol_prefac * (real_cell_vec(3,3) * &
         real_cell_vec(1,1) - real_cell_vec(3,1)*real_cell_vec(1,3))
    recip_cell_vec(2,3) = vol_prefac * (real_cell_vec(3,1) * &
         real_cell_vec(1,2) - real_cell_vec(3,2)*real_cell_vec(1,1))

    recip_cell_vec(3,1) = vol_prefac * (real_cell_vec(1,2) * &
         real_cell_vec(2,3) - real_cell_vec(1,3)*real_cell_vec(2,2))
    recip_cell_vec(3,2) = vol_prefac * (real_cell_vec(1,3) * &
         real_cell_vec(2,1) - real_cell_vec(1,1)*real_cell_vec(2,3))
    recip_cell_vec(3,3) = vol_prefac * (real_cell_vec(1,1) * &
         real_cell_vec(2,2) - real_cell_vec(1,2)*real_cell_vec(2,1))

    if(inode==ionode.AND.iprint_gen>1) then
       write(unit=io_lun,fmt='(/8x," cartesian components of &
            &reciprocal-space primitive vectors:"/)')
       do n = 1, 3
          write(unit=io_lun,fmt='(8x," basis vector no:",i3,":",3f12.6)') &
               n, (recip_cell_vec(n,i), i = 1, 3)
       enddo
    end if

    ! ++++++ temporary code to calculate lower bound on inter-particle separation
    ! --- radius of sphere inscribed in real-space cell --------------
    dummy1 = pi/sqrt(recip_cell_vec(1,1)**2 + recip_cell_vec(1,2)**2 +&
         recip_cell_vec(1,3)**2)
    dummy2 = pi/sqrt(recip_cell_vec(2,1)**2 + recip_cell_vec(2,2)**2 +&
         recip_cell_vec(2,3)**2)
    dummy3 = pi/sqrt(recip_cell_vec(3,1)**2 + recip_cell_vec(3,2)**2 +&
         recip_cell_vec(3,3)**2)
    lambda = min(dummy1,dummy2,dummy3)
    if(inode==ionode.AND.iprint_gen>1) &
         write(unit=io_lun,fmt='(/8x," radius of sphere inscribed in &
         &real-space cell:",f12.6)') lambda

    if(inode == ionode.AND.iprint_gen>1) &
         write(unit=io_lun,fmt='(/8x," no. of atoms in cell:",i5)') &
         ni_in_cell

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
         write(unit=io_lun,fmt='(/8x," lower bound on inter-particle &
         &separation:",f12.6)') lambda

    run_mode = 0
    select case(run_mode)

    case(0)
       ewald_gamma = (pi*(mean_number_density**two_thirds)/&
            real(ni_in_cell,double)**one_third)* (const_phi_0/&
            const_phi_1)
       argument = min(exp(-one), (ewald_accuracy*ewald_gamma)/(two*&
            sqrt(pi)*mean_interparticle_distance*mean_number_density))
       ewald_real_cutoff = (one/sqrt(ewald_gamma))*sqrt(-&
            log(argument))
       argument = min(exp(-one), (pi*ewald_accuracy)/&
            (mean_interparticle_distance*real(ni_in_cell,double)*&
            sqrt(ewald_gamma)))
       ewald_recip_cutoff = two*sqrt(ewald_gamma)*sqrt(-log(argument))

    case(1)
       ewald_gamma = (one/lambda**2)*log(one/ewald_accuracy) 
       ewald_real_cutoff = lambda
       ewald_recip_cutoff = sqrt(four*ewald_gamma*log(one/&
            ewald_accuracy))

    end select

    if(inode == ionode.AND.iprint_gen>1) then
       write(unit=io_lun,fmt='(/8x," run_mode:",i3/" Ewald gamma:",&
            &e20.12/" Ewald real_space cutoff distance:",e20.12/" &
            &Ewald recip-space cutoff wavevector:",e20.12)') run_mode,&
            ewald_gamma, ewald_real_cutoff, ewald_recip_cutoff
    endif

    ! --- Gaussian self-energy and net charge contributions to the Ewald energy
    ewald_gaussian_self_energy = -sqrt(ewald_gamma/pi)*&
         real(ni_in_cell,double)*mean_square_charge
    ewald_net_charge_energy = -(pi/(two*ewald_real_cell_volume*&
         ewald_gamma))*total_charge*total_charge
    if(inode == ionode.AND.iprint_gen>1) &
         write(unit=io_lun,fmt='(/8x," set_ewald: Ewald &
         &Gaussian self-energy:",e20.12," Ewald net-charge energy:",&
         &e20.12)') ewald_gaussian_self_energy, &
         ewald_net_charge_energy

    ! --- SYM 2014/10/16 20:20
    ! --- Stress due to Gaussian self-energy
    ! I'm assuming ewald_gaussian_self_stress is only a diagonal
    ! contribution - zamaan
    ewald_gaussian_self_stress = zero
    if (flag_stress) then
      gaussian_self_stress = pi * total_charge * total_charge/(two * ewald_real_cell_volume * ewald_gamma)
      do direction = 1, 3
         ewald_gaussian_self_stress(direction,direction) = gaussian_self_stress
      enddo
    end if

    ! --- limits for recip-space covering set
    rec_lim_1 = int(ewald_recip_cutoff* sqrt(real_cell_vec(1,1)**2 + &
         real_cell_vec(1,2)**2 + real_cell_vec(1,3)**2)/(two*pi))
    rec_lim_2 = int(ewald_recip_cutoff* sqrt(real_cell_vec(2,1)**2 + &
         real_cell_vec(2,2)**2 + real_cell_vec(2,3)**2)/(two*pi))
    rec_lim_3 = int(ewald_recip_cutoff* sqrt(real_cell_vec(3,1)**2 + &
         real_cell_vec(3,2)**2 + real_cell_vec(3,3)**2)/(two*pi))
    if(inode == ionode.AND.iprint_gen>1) write(unit=io_lun,fmt='(/8x,&
         &" integer limits for recip-space covering set:",3i5)') &
         rec_lim_1, rec_lim_2, rec_lim_3

    ! --- determine and store wavevectors to be included in recip-space sum
    ewald_recip_cutoff_squared = ewald_recip_cutoff*ewald_recip_cutoff
    number_of_g_vectors = 0
    do m1 = -rec_lim_1, rec_lim_1
       do m2 = -rec_lim_2, rec_lim_2
          do m3 = -rec_lim_3, rec_lim_3
             if(m1*m1+m2*m2+m3*m3 /= 0) then
                g1 = m1*recip_cell_vec(1,1) + m2*recip_cell_vec(2,1) +&
                     m3*recip_cell_vec(3,1)
                g2 = m1*recip_cell_vec(1,2) + m2*recip_cell_vec(2,2) +&
                     m3*recip_cell_vec(3,2)
                g3 = m1*recip_cell_vec(1,3) + m2*recip_cell_vec(2,3) +&
                     m3*recip_cell_vec(3,3)
                if(g1*g1 + g2*g2 + g3*g3 < ewald_recip_cutoff_squared) then
                   number_of_g_vectors = number_of_g_vectors + 1
                endif
             endif
          enddo
       enddo
    enddo
    ! Now we have number_of_g_vectors
    ! Allocate ewald_g_vector, g_factor
    allocate(ewald_g_vector_x(number_of_g_vectors),&
         ewald_g_vector_y(number_of_g_vectors), &
         ewald_g_vector_z(number_of_g_vectors),&
         ewald_g_factor(number_of_g_vectors),STAT=stat)
    ewald_g_vector_x = zero
    ewald_g_vector_y = zero
    ewald_g_vector_z = zero
    ewald_g_factor = zero
    if(stat/=0) &
         call cq_abort("set_ewald: error allocating ewald_g_vector ",&
         number_of_g_vectors,stat)
    call reg_alloc_mem(area_general,4*number_of_g_vectors,type_dbl)
    ! Zero because we use it to count
    number_of_g_vectors = 0
    do m1 = -rec_lim_1, rec_lim_1
       do m2 = -rec_lim_2, rec_lim_2
          do m3 = -rec_lim_3, rec_lim_3
             if(m1*m1+m2*m2+m3*m3 /= 0) then
                g1 = m1*recip_cell_vec(1,1) + m2*recip_cell_vec(2,1) +&
                     m3*recip_cell_vec(3,1)
                g2 = m1*recip_cell_vec(1,2) + m2*recip_cell_vec(2,2) +&
                     m3*recip_cell_vec(3,2)
                g3 = m1*recip_cell_vec(1,3) + m2*recip_cell_vec(2,3) +&
                     m3*recip_cell_vec(3,3)
                if(g1*g1 + g2*g2 + g3*g3 < ewald_recip_cutoff_squared) then
                   number_of_g_vectors = number_of_g_vectors + 1
                   g_squared = g1*g1 + g2*g2 + g3*g3
                   ewald_g_vector_x(number_of_g_vectors) = g1
                   ewald_g_vector_y(number_of_g_vectors) = g2
                   ewald_g_vector_z(number_of_g_vectors) = g3
                   ewald_g_factor(number_of_g_vectors) = &
                        exp(-g_squared/(four*ewald_gamma))/g_squared
                endif
             endif
          enddo
       enddo
    enddo
    ! ********************************************************************
    !              SET-UP INFORMATION FOR EWALD REAL-SPACE SUM
    ! ********************************************************************

    ! --- primitive translation vectors for partition
    if(inode == ionode.AND.iprint_gen>1) &
         write(unit=io_lun,fmt='(/8x," nos of partitions in cell edges:",3i5)') &
         parts%ngcellx, parts%ngcelly, parts%ngcellz
    do i = 1, 3
       part_cell_vec(1,i) = real_cell_vec(1,i) / real(parts%ngcellx,double)
       part_cell_vec(2,i) = real_cell_vec(2,i) / real(parts%ngcelly,double)
       part_cell_vec(3,i) = real_cell_vec(3,i) / real(parts%ngcellz,double)
    enddo
    if(inode == ionode.AND.iprint_gen>1) then
       write(unit=io_lun,fmt='(/8x," partition cell vectors:"/)')
       do n = 1, 3
          write(unit=io_lun,&
               fmt='(8x," partition cell vector no:",i5,3x,3f12.6)') &
               n, (part_cell_vec(n,i), i = 1, 3)
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
             vert(1) = v1*part_cell_vec(1,1) + v2*part_cell_vec(2,1) +&
                  v3*part_cell_vec(3,1)
             vert(2) = v1*part_cell_vec(1,2) + v2*part_cell_vec(2,2) +&
                  v3*part_cell_vec(3,2)
             vert(3) = v1*part_cell_vec(1,3) + v2*part_cell_vec(2,3) +&
                  v3*part_cell_vec(3,3)
             rsq = vert(1)*vert(1) + vert(2)*vert(2) + vert(3)*vert(3)
             if(rsq > partition_radius_squared) partition_radius_squared = rsq
          enddo
       enddo
    enddo
    partition_diameter = two*sqrt(partition_radius_squared)
    if(inode == ionode.AND.iprint_gen>1) &
         write(unit=io_lun,fmt='(/8x," diameter of sphere circumscribed round &
         &partition cell:",f12.6)') partition_diameter 

    ! --- make Ewald covering set
    call make_cs(inode-1,ewald_real_cutoff,ion_ion_CS,parts,bundle,&
         ni_in_cell, x_atom_cell,y_atom_cell,z_atom_cell)
    ! +++
    if(inode == ionode.AND.iprint_gen>1) then
       write(unit=io_lun,fmt='(/8x,"+++ ng_cover:",i10)') &
            ion_ion_CS%ng_cover
       write(unit=io_lun,fmt='(8x,"+++ ncoverx, y, z:",3i8)') &
            ion_ion_CS%ncoverx, ion_ion_CS%ncovery, ion_ion_CS%ncoverz
       write(unit=io_lun,fmt='(8x,"+++ nspanlx, y, z:",3i8)') &
            ion_ion_CS%nspanlx, ion_ion_CS%nspanly, ion_ion_CS%nspanlz
       write(unit=io_lun,fmt='(8x,"+++ nx_origin, y, z:",3i8)') &
            ion_ion_CS%nx_origin, ion_ion_CS%ny_origin, ion_ion_CS%nz_origin
    endif
    ! +++

    ! --- initiate partition_distance routine
    call partition_distance(.true.,part_cell_vec,part_cell_dual)

    ! --- loop over primary-set partitions of present node
    ! Work out n_partition_neighbours for first partition in group
    ip = 1

    ! --- Cartesian integer indices of current primary-set partition with respect
    !     to left of Ewald covering set (to be more precise: for a hypothetical primary-set 
    !     partition in bottom left-hand corner of covering set, all the following
    !     indices nppx, nppy, nppz would be equal to one).
    nppx = 1 + bundle%idisp_primx(ip) + ion_ion_CS%nspanlx
    nppy = 1 + bundle%idisp_primy(ip) + ion_ion_CS%nspanly
    nppz = 1 + bundle%idisp_primz(ip) + ion_ion_CS%nspanlz

    n_partition_neighbours = 0
    ! --- loop over partitions in Ewald covering set in node-periodic grouped order
    do nc = 1, ion_ion_CS%ng_cover
       ! --- for current partition nc in Ewald covering set, get the Cartesian-composite index
       ind_part = ion_ion_CS%lab_cover(nc)
       ! --- unpack the Cartesian-composite index to obtain the Cartesian integer indices of current
       !     partition nc in Ewald covering set
       nccx = 1 + (ind_part - 1)/(ion_ion_CS%ncovery * ion_ion_CS%ncoverz)
       nccy = 1 + (ind_part - 1 - (nccx-1) * ion_ion_CS%ncovery * &
            ion_ion_CS%ncoverz)/ion_ion_CS%ncoverz
       nccz = ind_part - (nccx-1) * ion_ion_CS%ncovery * ion_ion_CS%ncoverz - &
            (nccy-1) * ion_ion_CS%ncoverz
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
          rr(1) = half*(v1*part_cell_vec(1,1) + v2*part_cell_vec(2,1) &
               + v3*part_cell_vec(3,1))
          rr(2) = half*(v1*part_cell_vec(1,2) + v2*part_cell_vec(2,2) &
               + v3*part_cell_vec(3,2))
          rr(3) = half*(v1*part_cell_vec(1,3) + v2*part_cell_vec(2,3) &
               + v3*part_cell_vec(3,3))
          call partition_distance(.false.,part_cell_vec,part_cell_dual,rr,distance)
          distance = two*distance
          if(distance < ewald_real_cutoff) then
             n_partition_neighbours = n_partition_neighbours + 1
          endif
       endif
    enddo
    ! Allocate memory
    allocate(partition_neighbour_list(bundle%groups_on_node,&
         n_partition_neighbours), STAT=stat)
    if(stat/=0) &
         call cq_abort("set_ewald: error allocating partition_neighbour_list ", &
         bundle%groups_on_node,n_partition_neighbours)
    call reg_alloc_mem(area_general,bundle%groups_on_node*&
         n_partition_neighbours,type_dbl)

    do ip = 1, bundle%groups_on_node

       ! --- Cartesian integer indices of current primary-set
       !     partition with respect to left of Ewald covering set (to
       !     be more precise: for a hypothetical primary-set partition
       !     in bottom left-hand corner of covering set, all the
       !     following indices nppx, nppy, nppz would be equal to
       !     one).
       nppx = 1 + bundle%idisp_primx(ip) + ion_ion_CS%nspanlx
       nppy = 1 + bundle%idisp_primy(ip) + ion_ion_CS%nspanly
       nppz = 1 + bundle%idisp_primz(ip) + ion_ion_CS%nspanlz

       n_partition_neighbours = 0
       ! --- loop over partitions in Ewald covering set in
       !      node-periodic grouped order
       do nc = 1, ion_ion_CS%ng_cover
          ! --- for current partition nc in Ewald covering set, get
          !     the Cartesian-composite index
          ind_part = ion_ion_CS%lab_cover(nc)
          ! --- unpack the Cartesian-composite index to obtain the
          !     Cartesian integer indices of current partition nc in
          !     Ewald covering set
          nccx = 1 + (ind_part - 1)/(ion_ion_CS%ncovery * ion_ion_CS%ncoverz)
          nccy = 1 + (ind_part - 1 - (nccx-1) * ion_ion_CS%ncovery * &
               ion_ion_CS%ncoverz)/ion_ion_CS%ncoverz
          nccz = ind_part - (nccx-1) * ion_ion_CS%ncovery * ion_ion_CS%ncoverz - &
               (nccy-1) * ion_ion_CS%ncoverz
          ! --- integer triplet specifying offset of covering-set
          !     partition from primary-set partition
          m1 = nccx - nppx
          m2 = nccy - nppy
          m3 = nccz - nppz
          ! --- exclude from list if all components of offset are zero
          if(m1*m1 + m2*m2 + m3*m3 /= 0) then
             ! --- calculate distance (in a.u.) between centre of
             !     primary-set partition and the point mid-way between
             !     the primary- and covering-set partitions
             v1 = real(m1,double)
             v2 = real(m2,double)
             v3 = real(m3,double)
             rr(1) = half*(v1*part_cell_vec(1,1) + v2*part_cell_vec(2,1) + &
                  v3*part_cell_vec(3,1))
             rr(2) = half*(v1*part_cell_vec(1,2) + v2*part_cell_vec(2,2) + &
                  v3*part_cell_vec(3,2))
             rr(3) = half*(v1*part_cell_vec(1,3) + v2*part_cell_vec(2,3) + &
                  v3*part_cell_vec(3,3))
             call partition_distance(.false.,part_cell_vec,&
                  part_cell_dual,rr,distance)
             distance = two*distance
             if(distance < ewald_real_cutoff) then
                n_partition_neighbours = n_partition_neighbours + 1
                ! --- for current primary-set partition ip, record the
                !     PG index nc of latest covering-set partition
                !     within cutoff
                partition_neighbour_list(ip,&
                     n_partition_neighbours) = nc
             endif
          endif
       enddo
    enddo

    allocate(ion_interaction_force(3,ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating ion_interaction_force in set_ewald: ",ni_in_cell,stat)
    call reg_alloc_mem(area_general,3*ni_in_cell,type_dbl)
    ion_interaction_force = zero
    call my_barrier()

    return

  end subroutine set_ewald
!!***

!!****f* ion_electrostatic/structure_factor *
!!
!!  NAME 
!!   structure_factor
!!  USAGE
!!   
!!  PURPOSE
!!   Calculates structure factors
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   M. J. Gillan
!!  CREATION DATE
!!   2003
!!  MODIFICATION HISTORY
!!   2016/01/29 14:58 dave
!!    Rename to remove prefix mikes_
!!  SOURCE
!!  
  subroutine structure_factor(inode,ionode)

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
    return
  end subroutine structure_factor
!!***

!!****f* ion_electrostatic/ewald *
!!
!!  NAME 
!!   mikes_ewald
!!  USAGE
!!   
!!  PURPOSE
!!   Calculates the ion-ion interaction energy, forces and stresses using the Ewald sum
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   M. J. Gillan
!!  CREATION DATE
!!   2003
!!  MODIFICATION HISTORY
!!   10:42, 13/02/2006 drb 
!!    Corrected line length error
!!   2015/05/01 13:29 dave and sym
!!     Added stress
!!   2016/01/29 14:58 dave
!!     Removed prefix mikes_
!!   2016/02/09 08:10 dave  
!!    Added use statement for erfc
!!   2019/04/09 zamaan
!!    Modifications to compute off-diagonal stress tensor elements
!!   2019/05/08 zamaan
!!    Added atomic stress contributions
!!   2020/09/08 17:18 dave
!!    Bug fix for reciprocal-space component of stress
!!  SOURCE
!!  
  subroutine ewald()

    use cover_module,   ONLY: ion_ion_CS
    use datatypes
    use GenComms,       ONLY: gsum, cq_abort, inode, ionode
    use global_module,  ONLY: id_glob, iprint_gen, species_glob, ni_in_cell, &
                              area_general, IPRINT_TIME_THRES3, &
                              flag_full_stress, flag_stress, &
                              flag_atomic_stress, atomic_stress, &
                              non_atomic_stress
    use group_module,   ONLY: parts
    use maxima_module,  ONLY: maxatomsproc
    use numbers
    use primary_module, ONLY: bundle
    use species_module, ONLY: charge, species
    use memory_module,  ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module,   ONLY: cq_timer,start_timer,stop_timer, &
                              stop_print_timer,WITH_LEVEL
    use timer_module,   ONLY: start_backtrace,stop_backtrace
    use functions,      ONLY: erfc_cq
    use energy,         ONLY: ion_interaction_energy
    
    implicit none

    integer :: i, ip, ig_atom_beg, j, n, nc, ni, nj, nn, stat, direction, dir1, dir2, iatom, jatom
    real(double) :: arg_1, arg_2, coeff_1, dummy, ewald_real_energy, ewald_real_sum_inter, ewald_real_sum_intra, &
         &ewald_real_sum_ip, ewald_recip_energy, ewald_total_energy, g_dot_r, q_i, q_j, rij, rijx, rijy, rijz, &
         rij_squared, x, y, z, one_over_g_squared, delta_ab, one_over_four_gamma
    real(double), dimension(3) :: egv_n, rij_vec
    real(double), allocatable, dimension(:) :: ewald_recip_force_x, ewald_recip_force_y, ewald_recip_force_z, &
         &ewald_intra_force_x, ewald_intra_force_y, ewald_intra_force_z, &
         &ewald_inter_force_x, ewald_inter_force_y, ewald_inter_force_z
    real(double), allocatable, dimension(:) :: ewald_total_force_x, ewald_total_force_y, ewald_total_force_z
    type(cq_timer) :: tmr_l_tmp1
    type(cq_timer) :: backtrace_timer

!****lat<$
    call start_backtrace(t=backtrace_timer,who='ewald',where=9,level=2)
!****lat>$

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
    ! Zero stress
    ewald_recip_stress = zero
    ewald_intra_stress = zero
    ewald_inter_stress = zero
    ion_interaction_stress = zero
    ! Allocate memory for and call structure_factor
    allocate(struc_fac_r(number_of_g_vectors),struc_fac_i(number_of_g_vectors),STAT=stat)
    if(stat/=0) call cq_abort("ewald: error allocating struc_fac ",number_of_g_vectors,stat)
    call reg_alloc_mem(area_general,2*number_of_g_vectors,type_dbl)
    call structure_factor(inode,ionode)
    ! ... loop over reciprocal lattice vectors
    one_over_four_gamma = one/(four*ewald_gamma)
    do n = 1, number_of_g_vectors
       ewald_recip_energy = ewald_recip_energy + &
            &(struc_fac_r(n)*struc_fac_r(n) + struc_fac_i(n)*struc_fac_i(n)) * &
            &ewald_g_factor(n)
       ! Stress added SYM 2014/10/16 21:38
       ! Here we calculate the stress from the reciprocal space part of the Ewald sum
       ! I'm creating a new 3-vector to make the double loop more 
       ! compact - zamaan
       egv_n(1) = ewald_g_vector_x(n)
       egv_n(2) = ewald_g_vector_y(n)
       egv_n(3) = ewald_g_vector_z(n)

       one_over_g_squared = &
         one/(egv_n(1)*egv_n(1) + egv_n(2)*egv_n(2) + egv_n(3)*egv_n(3))
       if (flag_stress) then
          do dir1=1,3
           if (flag_full_stress) then
              do dir2=1,3
                 if(dir2==dir1) then
                    delta_ab = one
                 else
                    delta_ab = zero
                 end if
                 ewald_recip_stress(dir1,dir2) = ewald_recip_stress(dir1,dir2)+ &
                      ewald_g_factor(n) * (struc_fac_r(n)*struc_fac_r(n) + &
                      struc_fac_i(n)*struc_fac_i(n)) * &
                      (two*egv_n(dir1)*egv_n(dir2) * &
                      (one_over_four_gamma + one_over_g_squared) - delta_ab)
             end do
           else
             ewald_recip_stress(dir1,dir1) = ewald_recip_stress(dir1,dir1) + &
                  ewald_g_factor(n) * (struc_fac_r(n)*struc_fac_r(n) + &
                  struc_fac_i(n)*struc_fac_i(n)) * &
                  (two*egv_n(dir1)*egv_n(dir1) * &
                  (one_over_four_gamma + one_over_g_squared) - one)
           end if
         end do
       end if

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
    ! Correctly scale stress
    if (flag_stress) then
       do dir1=1,3 
        do dir2=1,3
          ewald_recip_stress(dir1,dir2) = ewald_recip_stress(dir1,dir2) * &
                                          two * pi/(ewald_real_cell_volume)
          end do
      end do
    end if
    ! 2015/05/01 13:35 dave
    !  I think that the reciprocal space part doesn't need this sum, but I should check
    !call gsum(ewald_recip_stress,3)

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
                q_j = charge(bundle%species(bundle%nm_nodbeg(ip)+nj-1))
                ! ... Cartesian components of vector separation of the two atoms
                iatom = bundle%nm_nodbeg(ip)+ni-1
                jatom = bundle%nm_nodbeg(ip)+nj-1
                rij_vec(1) = bundle%xprim(iatom) - bundle%xprim(jatom)
                rij_vec(2) = bundle%yprim(iatom) - bundle%yprim(jatom)
                rij_vec(3) = bundle%zprim(iatom) - bundle%zprim(jatom)
                rij_squared = rij_vec(1)*rij_vec(1) + rij_vec(2)*rij_vec(2) + &
                              rij_vec(3)*rij_vec(3)
                rij = sqrt(rij_squared)
                if(rij < ewald_real_cutoff) then
                   arg_1 = ewald_gamma*rij_squared
                   arg_2 = sqrt(arg_1)
                   ewald_real_sum_ip = ewald_real_sum_ip + &
                     q_i * q_j * erfc_cq(arg_2) / rij
                   dummy = q_i * q_j * (erfc_cq(arg_2)/rij + coeff_1*exp(-arg_1))/rij_squared
                   ewald_intra_force_x(iatom) = &
                     ewald_intra_force_x(iatom) + dummy*rij_vec(1)
                   ewald_intra_force_x(jatom) = &
                     ewald_intra_force_x(jatom) - dummy*rij_vec(1)
                   ewald_intra_force_y(iatom) = &
                     ewald_intra_force_y(iatom) + dummy*rij_vec(2)
                   ewald_intra_force_y(jatom) = &
                     ewald_intra_force_y(jatom) - dummy*rij_vec(2)
                   ewald_intra_force_z(iatom) = &
                     ewald_intra_force_z(iatom) + dummy*rij_vec(3)
                   ewald_intra_force_z(jatom) = &
                     ewald_intra_force_z(jatom) - dummy*rij_vec(3)
                   ! --- Edited SYM 2014/10/16 14:23 Ewald stress
                   if (flag_stress) then
                     do dir1=1,3
                       if (flag_full_stress) then
                         do dir2=1,3
                           ewald_intra_stress(dir1,dir2) = &
                             ewald_intra_stress(dir1,dir2) - &
                             (dummy * rij_vec(dir1) * rij_vec(dir2))
                           if (flag_atomic_stress) then
                             ! how is this contributions split between the 
                             ! two atoms? - zamaan
                             atomic_stress(dir1,dir2,iatom) = &
                               atomic_stress(dir1,dir2,iatom) - &
                               (dummy * rij_vec(dir1) * rij_vec(dir2)) * half
                           end if
                         end do ! dir2
                       else
                         ewald_intra_stress(dir1,dir1) = &
                           ewald_intra_stress(dir1,dir1) - &
                           (dummy * rij_vec(dir1) * rij_vec(dir1))
                       end if
                     end do ! dir1
                   end if
                endif
             enddo
          enddo
       endif
       if(inode == ionode.AND.iprint_gen>=4) &
            write(unit=io_lun,fmt='(/" >>> ewald: node:",i3:" real_sum_intra for partition no:",i3,&
            &" is:",e20.12)') inode, ip, ewald_real_sum_ip
       ewald_real_sum_intra = ewald_real_sum_intra + ewald_real_sum_ip
    enddo

    call gsum(ewald_real_sum_intra)
    if (flag_stress) call gsum(ewald_intra_stress,3,3)
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
             iatom = bundle%nm_nodbeg(ip)+ni-1
             i = bundle%ig_prim(iatom)
             q_i = charge(bundle%species(iatom))
             ! --- loop over Ewald neighbour-list of partitions in Ewald covering set
             do nn = 1, n_partition_neighbours
                ! ... get PG label of current ECS partition
                nc = partition_neighbour_list(ip,nn)
                ig_atom_beg = parts%icell_beg(ion_ion_CS%lab_cell(nc))
                ! --- loop over atoms in ECS partition
                if(ion_ion_CS%n_ing_cover(nc) > 0) then
                   do nj = 1, ion_ion_CS%n_ing_cover(nc)
                      ! ... global label and charge of current atom in ECS partition
                      q_j = charge(species_glob(id_glob(ig_atom_beg+nj-1)))
                      ! ... Cartesian components of vector separation of the two atoms
                      rij_vec(1) = bundle%xprim(iatom) - &
                        ion_ion_CS%xcover(ion_ion_CS%icover_ibeg(nc)+nj-1)
                      rij_vec(2) = bundle%yprim(iatom) - &
                        ion_ion_CS%ycover(ion_ion_CS%icover_ibeg(nc)+nj-1)
                      rij_vec(3) = bundle%zprim(iatom) - &
                        ion_ion_CS%zcover(ion_ion_CS%icover_ibeg(nc)+nj-1)
                      rij_squared = rij_vec(1)*rij_vec(1) + &
                        rij_vec(2)*rij_vec(2) + rij_vec(3)*rij_vec(3)
                      rij = sqrt(rij_squared)
                      if(rij < ewald_real_cutoff) then
                         arg_1 = ewald_gamma*rij_squared
                         arg_2 = sqrt(arg_1)
                         ewald_real_sum_ip = ewald_real_sum_ip + &
                              &q_i * q_j * erfc_cq(arg_2) / rij
                         dummy = q_i * q_j * (erfc_cq(arg_2)/rij + coeff_1*exp(-arg_1))/rij_squared
                         ewald_inter_force_x(iatom) = &
                              ewald_inter_force_x(iatom) + dummy*rij_vec(1)
                         ewald_inter_force_y(iatom) = &
                              ewald_inter_force_y(iatom) + dummy*rij_vec(2)
                         ewald_inter_force_z(iatom) = &
                              ewald_inter_force_z(iatom) + dummy*rij_vec(3)
                         if (flag_stress) then
                           do dir1=1,3
                             if (flag_full_stress) then
                               do dir2=1,3
                                 ewald_inter_stress(dir1,dir2) = &
                                   ewald_inter_stress(dir1,dir2) - &
                                   dummy * rij_vec(dir1) * rij_vec(dir2)
                                 if (flag_atomic_stress) then
                                   atomic_stress(dir1,dir2,iatom) = &
                                     atomic_stress(dir1,dir2,iatom) - &
                                     dummy * rij_vec(dir1) * &
                                     rij_vec(dir2) * half
                                 end if
                               end do ! dir2
                             else
                               ewald_inter_stress(dir1,dir1) = &
                                 ewald_inter_stress(dir1,dir1) - &
                                 dummy * rij_vec(dir1) * rij_vec(dir1)
                             end if
                           end do ! dir1
                         end if
                      endif
                   enddo
                endif
             enddo
          enddo
       endif
       if(inode == ionode.AND.iprint_gen>=4) &
            write(unit=io_lun,fmt='(/" >>> ewald: node:",i3:" real_sum_inter for partition no:",i3,&
            &" is:",e20.12)') inode, ip, ewald_real_sum_ip
       ewald_real_sum_inter = ewald_real_sum_inter + ewald_real_sum_ip
    enddo
    if (flag_stress) then
      do dir1=1,3
         do dir2=1,3
            ewald_inter_stress(dir1,dir2) = half*ewald_inter_stress(dir1,dir2)
         end do
      end do
    end if

    call gsum(ewald_real_sum_inter)
    if (flag_stress) call gsum(ewald_inter_stress,3,3)
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
    ion_interaction_energy = ewald_total_energy
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
    ! SYM 2014/10/22 14:34: Summ all the stress contributions
    if (flag_stress) then
       do dir1=1,3
         do dir2=1,3
            ion_interaction_stress(dir1,dir2) = &
              ewald_intra_stress(dir1,dir2) + ewald_inter_stress(dir1,dir2) + &
              ewald_recip_stress(dir1,dir2) + &
              ewald_gaussian_self_stress(dir1,dir2)
            if (flag_atomic_stress) then
              non_atomic_stress(dir1,dir2) = &
                non_atomic_stress(dir1,dir2) + &
                ewald_recip_stress(dir1,dir2) + &
                ewald_gaussian_self_stress(dir1,dir2)
            end if
          end do
      end do
    end if
    ! Added DRB & MJG 2003/10/22 to make force available to rest of code
    do i=1,ni_in_cell
       ion_interaction_force(1,i) = ewald_total_force_x(i)
       ion_interaction_force(2,i) = ewald_total_force_y(i)
       ion_interaction_force(3,i) = ewald_total_force_z(i)
    end do
    if(inode == ionode.AND.iprint_gen>=2) then
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
    call stop_print_timer(tmr_l_tmp1,"ewald",IPRINT_TIME_THRES3)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='ewald')
!****lat>$

    return
  end subroutine ewald
!!***

!!****f* ion_electrostatic/setup_screened_ion_interaction *
!!
!!  NAME 
!!   setup_screened_ion_interaction
!!  USAGE
!!   
!!  PURPOSE
!!   Sets up the necessary arrays and covering set for calculating ion-ion interactions
!!   within the neutral atom method (using same arrays and loops as Ewald sum but shorter range)
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2015/12/09
!!  MODIFICATION HISTORY
!!   2016/01/14 to 2016/01/28
!!    Removed the old-style ewald approach and replaced with simple loop over neighbours
!!    which vastly simplifies the calculation and the code
!!   2016/02/04 17:27 dave
!!    Bug fix: the cutoff needed to be doubled
!!   2018/11/13 TM
!!    Bug fix for ghost atoms
!!  SOURCE
!!  
  subroutine setup_screened_ion_interaction

    use datatypes
    use GenComms, ONLY : cq_abort, my_barrier, inode, ionode
    use group_module, ONLY : parts
    use global_module, ONLY : x_atom_cell, y_atom_cell, z_atom_cell, &
         iprint_gen, ni_in_cell, area_general
    use numbers
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use cover_module, ONLY : ion_ion_CS, make_cs
    use primary_module, ONLY : bundle
    use species_module, only: n_species, type_species ! for Neutral atom potential
    use atomic_density, only: atomic_density_table ! for Neutral atom potential

    implicit none

    integer :: stat, isp

    ! Find the cutoff distance
    ion_ion_cutoff = zero
    do isp=1, n_species
      if(type_species(isp) < 0) cycle
       if( ion_ion_cutoff<atomic_density_table(isp)%cutoff ) then
          ion_ion_cutoff=atomic_density_table(isp)%cutoff
       end if
    end do
    ! The cutoff needs to be twice the radius so that atoms with overlapping charges interact
    ! It might be possible to shrink this to reduce the time at the cost of reduced accuracy,
    ! but the time required should be small
    ion_ion_cutoff = two*ion_ion_cutoff
    ! --- make covering set
    call make_cs(inode-1,ion_ion_cutoff,ion_ion_CS,parts,bundle,&
         ni_in_cell, x_atom_cell,y_atom_cell,z_atom_cell)
    ! +++ Output if necessary
    if(inode == ionode.AND.iprint_gen>2) then
       write(unit=io_lun,fmt='(/8x,"+++ ng_cover:",i10)') &
            ion_ion_CS%ng_cover
       write(unit=io_lun,fmt='(8x,"+++ ncoverx, y, z:",3i8)') &
            ion_ion_CS%ncoverx, ion_ion_CS%ncovery, ion_ion_CS%ncoverz
       write(unit=io_lun,fmt='(8x,"+++ nspanlx, y, z:",3i8)') &
            ion_ion_CS%nspanlx, ion_ion_CS%nspanly, ion_ion_CS%nspanlz
       write(unit=io_lun,fmt='(8x,"+++ nx_origin, y, z:",3i8)') &
            ion_ion_CS%nx_origin, ion_ion_CS%ny_origin, ion_ion_CS%nz_origin
    endif
    ! And allocate forces
    allocate(screened_ion_force(3,ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating screened_ion_force in setup_ion_interaction: ",ni_in_cell,stat)
    call reg_alloc_mem(area_general,3*ni_in_cell,type_dbl)
    screened_ion_force = zero
    call my_barrier()

  end subroutine setup_screened_ion_interaction
!!***

!!****f* ion_electrostatic/screened_ion_interaction *
!!
!!  NAME 
!!   screened_ion_interaction
!!  USAGE
!!   
!!  PURPOSE
!!   Calculates the ion-ion interaction within the neutral atom approach
!!   (based on the standard Ewald routine for looping over atom pairs)
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   D. R. Bowler with TM and NW (Mizuho)
!!  CREATION DATE
!!   2015/12/22
!!  MODIFICATION HISTORY
!!   2016/01/14 to 2016/01/28
!!    Removed the old-style ewald approach and replaced with simple loop over neighbours
!!   2016/02/02 11:06 dave
!!    Changed force reference - now saves to correct atoms
!!   2017/03/23 dave
!!    Changed cutoff criterion on interactions to use individual species size
!!   2018/11/13 TM
!!    Bug fix for ghost atoms
!!   2019/05/08 zamaan
!!    Added atomic stress contributions
!!  SOURCE
!!  
  subroutine screened_ion_interaction

    use cover_module, ONLY : ion_ion_CS
    use datatypes
    use GenComms, ONLY : gsum, cq_abort, inode, ionode
    use global_module, ONLY : id_glob, iprint_gen, species_glob, ni_in_cell, &
                              area_general, IPRINT_TIME_THRES3, &
                              flag_full_stress, flag_stress, &
                              flag_atomic_stress, atomic_stress, min_layer
    use group_module, ONLY : parts
    use maxima_module, ONLY : maxatomsproc
    use numbers
    use primary_module, ONLY : bundle
    use species_module, ONLY : charge, species, type_species
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module, ONLY: cq_timer,start_timer,stop_print_timer,WITH_LEVEL
    use energy, ONLY: local_ps_energy, screened_ion_interaction_energy
    use atomic_density, only: atomic_density_table ! for Neutral atom potential
    use io_module, only: return_prefix
    use units, only: energy_units, en_units
    
    implicit none

    integer :: i, ip, ig_atom_beg, j, n, nc, ni, nj, nn, stat, direction, spi, spj, dir1, dir2
    real(double), dimension(3) :: rij_vec, goverlap
    real(double) :: q_i, q_j, rij, rijx, rijy, rijz, rij_squared
    real(double) :: screenedE_sum_intra, screenedE_sum_inter
    real(double) :: screenedE_sum_self
    real(double) :: overlap, dummy
    real(double) :: goverlap_x, goverlap_y, goverlap_z
    character(len=12) :: subname = "Screened ion int: "
    character(len=120) :: prefix

    prefix = return_prefix(subname, min_layer)
    if(iprint_gen>4) write(io_lun,fmt='(4x,a)') trim(prefix)//" Starting"
    screened_ion_force   = zero
    screened_ion_stress  = zero
    screened_ion_interaction_energy = zero
    if(iprint_gen>4) write(io_lun,fmt='(4x,a)') trim(prefix)//" Interactions"
    ! --- loop over primary-set partitions
    do ip = 1, bundle%groups_on_node
       if(bundle%nm_nodgroup(ip) > 0) then
          do ni = 1, bundle%nm_nodgroup(ip)
             i = bundle%ig_prim(bundle%nm_nodbeg(ip)+ni-1)
             spi = bundle%species(bundle%nm_nodbeg(ip)+ni-1)
             if(type_species(spi) < 0) cycle
             q_i = charge(spi)
             do nc = 1, ion_ion_CS%ng_cover
                ig_atom_beg = parts%icell_beg(ion_ion_CS%lab_cell(nc))

                if(ion_ion_CS%n_ing_cover(nc) > 0) then
                   do nj = 1, ion_ion_CS%n_ing_cover(nc)
                      spj = species_glob(id_glob(ig_atom_beg+nj-1))
                      if(type_species(spj) < 0) cycle
                      q_j = charge(spj)
                      rij_vec(1) = bundle%xprim(bundle%nm_nodbeg(ip)+ni-1) - &
                        ion_ion_CS%xcover(ion_ion_CS%icover_ibeg(nc)+nj-1)
                      rij_vec(2) = bundle%yprim(bundle%nm_nodbeg(ip)+ni-1) - &
                        ion_ion_CS%ycover(ion_ion_CS%icover_ibeg(nc)+nj-1)
                      rij_vec(3) = bundle%zprim(bundle%nm_nodbeg(ip)+ni-1) - &
                        ion_ion_CS%zcover(ion_ion_CS%icover_ibeg(nc)+nj-1)
                      rij_squared = rij_vec(1)*rij_vec(1) + &
                        rij_vec(2)*rij_vec(2) + rij_vec(3)*rij_vec(3)
                      rij = sqrt(rij_squared)
                      if(rij<(atomic_density_table(spi)%cutoff + atomic_density_table(spj)%cutoff)) then 
                         if(rij>very_small) then
                            ! Find electrostatic energy due to atomic densities i and j
                            call calc_overlap( overlap, spi, spj,rij )
                            ! and the differential
                            call calc_goverlap( goverlap(1), goverlap(2), &
                              goverlap(3), spi, spj, rij_vec(1), rij_vec(2), &
                              rij_vec(3), rij )
                            ! Energy is the difference between atomic densities and point charges
                            screened_ion_interaction_energy = screened_ion_interaction_energy &
                                 + half*(q_i*q_j/rij - overlap)
                            ! Forces and stresses
                            dummy = q_i * q_j /(rij_squared*rij)

                            ! I'm really just guessing here - zamaan
                            ! NB the factor of half here comes from summing over all i and j
                            do dir1=1,3
                              screened_ion_force(dir1,i) = &
                                screened_ion_force(dir1,i) + &
                                dummy*rij_vec(dir1) + goverlap(dir1)
                              if (flag_stress) then
                                if (flag_full_stress) then
                                  do dir2=1,3
                                    screened_ion_stress(dir1,dir2) = &
                                      screened_ion_stress(dir1,dir2) - &
                                      half * &
                                      (dummy*rij_vec(dir1) + goverlap(dir1))* &
                                      rij_vec(dir2)
                                      if (flag_atomic_stress) then
                                        atomic_stress(dir1,dir2,i) = &
                                          atomic_stress(dir1,dir2,i) - &
                                          half * (dummy*rij_vec(dir1) + &
                                          goverlap(dir1)) * rij_vec(dir2)
                                      end if
                                  end do
                                else
                                  screened_ion_stress(dir1,dir1) = &
                                    screened_ion_stress(dir1,dir1) - &
                                    half * &
                                    (dummy*rij_vec(dir1) + goverlap(dir1)) * &
                                    rij_vec(dir1)
                                end if
                              end if
                            end do
                         else ! i=j we just need -0.5*overlap
                            ! Find electrostatic energy due to atomic densities i and j
                            call calc_overlap( overlap, spi, spj,rij )
                            screened_ion_interaction_energy = screened_ion_interaction_energy - half*overlap
                         end if
                      end if
                   enddo
                endif
             enddo
          enddo
       endif
    enddo
    call gsum(screened_ion_interaction_energy)
    call gsum(screened_ion_force,3,ni_in_cell)
    if (flag_stress) call gsum(screened_ion_stress,3,3)
    if(inode == ionode.AND.iprint_gen>1) &
         write(unit=io_lun,fmt='(/4x,a,e20.12,a2/)') &
         trim(prefix)//" Energy: ", screened_ion_interaction_energy,en_units(energy_units)

  end subroutine screened_ion_interaction
!!***
  
!!****f* ion_electrostatic/calc_overlap *
!!
!!  NAME 
!!   calc_overlap
!!  USAGE
!!   
!!  PURPOSE
!!   Calculates the electrostatic energy for atomic charge density
!!   for species 1 in the Hartree potential from atomic density for species 2
!!  INPUTS
!!
!!   
!!  USES
!!   
!!  AUTHOR
!!   NW (Mizuho) with TM and DRB
!!  CREATION DATE
!!   Mid 2014
!!  MODIFICATION HISTORY
!!  SOURCE
!!  
  subroutine calc_overlap( overlap, isp1, isp2, distance )

    use numbers
    use atomic_density, ONLY: atomic_density_table
    use functions, ONLY: j0

    implicit none

    real(double), intent(out) :: overlap
    integer, intent(in) :: isp1, isp2
    real(double), intent(in) :: distance

    integer :: ik
    real(double) :: dk, k, value

    overlap = zero
    dk = atomic_density_table(isp1)%k_delta
    do ik=0, atomic_density_table(isp1)%k_length
       k = real(ik,double)*dk
       value = j0(k*distance) &
            * atomic_density_table(isp1)%k_table(ik) &
            * atomic_density_table(isp2)%k_table(ik)

       if( ik==0 ) then
          overlap = overlap + dk*value*half
       else
          overlap = overlap + dk*value
       end if
    end do
    overlap = overlap*four*pi*four*pi

    return
  end subroutine calc_overlap
!!***

!!****f* ion_electrostatic/calc_goverlap *
!!
!!  NAME 
!!   calc_goverlap
!!  USAGE
!!   
!!  PURPOSE
!!   Calculates the differential of the electrostatic energy for atomic charge density
!!   for species 1 in the Hartree potential from atomic density for species 2
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   NW (Mizuho) with TM and DRB
!!  CREATION DATE
!!   Mid 2014
!!  MODIFICATION HISTORY
!!  
!!  SOURCE
!!  
  subroutine calc_goverlap( &
       goverlap_x, goverlap_y, goverlap_z, isp1, isp2, &
       distance_x, distance_y, distance_z, distance )
    
    use numbers
    use atomic_density, ONLY: atomic_density_table
    use functions, ONLY: j1x

    implicit none

    real(double), intent(out) :: goverlap_x
    real(double), intent(out) :: goverlap_y
    real(double), intent(out) :: goverlap_z
    integer, intent(in) :: isp1, isp2
    real(double), intent(in) :: distance_x
    real(double), intent(in) :: distance_y
    real(double), intent(in) :: distance_z
    real(double), intent(in) :: distance

    integer :: ik
    real(double) :: dk, k, doverlap

    goverlap_x = zero
    goverlap_y = zero
    goverlap_z = zero

    dk = atomic_density_table(isp1)%k_delta
    do ik=0, atomic_density_table(isp1)%k_length
       k = real(ik,double)*dk
       doverlap = minus_one*k*k*j1x(k*distance) &
            * atomic_density_table(isp1)%k_table(ik) &
            * atomic_density_table(isp2)%k_table(ik)

       if( ik==0 ) then
          goverlap_x = goverlap_x + dk*doverlap*distance_x*half
          goverlap_y = goverlap_y + dk*doverlap*distance_y*half
          goverlap_z = goverlap_z + dk*doverlap*distance_z*half
       else
          goverlap_x = goverlap_x + dk*doverlap*distance_x
          goverlap_y = goverlap_y + dk*doverlap*distance_y
          goverlap_z = goverlap_z + dk*doverlap*distance_z
       end if
    end do
    goverlap_x = goverlap_x*four*pi*four*pi
    goverlap_y = goverlap_y*four*pi*four*pi
    goverlap_z = goverlap_z*four*pi*four*pi

    return
  end subroutine calc_goverlap
!!***

end module ion_electrostatic
