! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_solve
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_solve *
!!
!!NAME
!! cq_ut_psp2pao_solve
!!PURPOSE
!! Equations are solved here
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_numeric
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 14/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_solve

  implicit none

contains

! -----------------------------------------------------------------------------
! Subroutine solve_eigenstates
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_solve/solve_eigenstates *
!!
!!NAME
!! solve_eigenstates
!!USAGE
!!
!!PURPOSE
!! Get the eigenstates
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 14/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine solve_eigenstates

     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types

     implicit none

     integer :: orbital, i
     integer :: ini_orbital, ini_v_nonlocal
     integer :: no_loop
     integer :: shift
     real(double) :: integ_diff      ! Difference of integral, for convergence
     real(double) :: norm_old, norm_diff, diff

     ! For the moment, UNPOLARISED CASE
     ini_orbital = 0     ! This would be gl_no_orbitals if e.g. spin = dn
     ini_v_nonlocal = 0  ! This would be gl_no_orbitals if e.g. spin = dn

     do orbital=1,gl_no_orbitals
        no_loop = 0
        integ_diff = 1.0_double

        do while(integ_diff > 5.0e-13_double .and. no_loop < 1000)
           norm_old = gl_norm(ini_orbital + orbital)
           call solve_orbital(orbital)
           norm_diff = gl_norm(ini_orbital + orbital) - norm_old 

           integ_diff = 0.0_double
           do i=1,gl_psp_in%psp_comp
              shift = ini_v_nonlocal + (orbital-1) * gl_psp_in%psp_comp + i
              diff = gl_v_nonlocal_int_out(shift) - gl_v_nonlocal_int_in(shift)
              gl_v_nonlocal_int_in(shift) = gl_v_nonlocal_int_in(shift) &
                                          + gl_orbital_mix * diff
              integ_diff = integ_diff + diff * diff
           end do
           integ_diff = sqrt(integ_diff)
           gl_norm(ini_orbital + orbital) = norm_old + gl_orbital_mix * norm_diff

           no_loop = no_loop + 1
        end do
     end do

  end subroutine solve_eigenstates
!!***


! -----------------------------------------------------------------------------
! Subroutine solve_orbital
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_solve/solve_orbital *
!!
!!NAME
!! solve_orbital
!!USAGE
!!
!!PURPOSE
!! Solve Schroedinger and refine eigenvalues
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_numeric
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 14/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine solve_orbital(orbital)

     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_numeric

     implicit none

     !! Passed variables
     integer :: orbital

     !! Local variables
     integer :: i, j
     integer :: shift, l1
     integer :: ini_ul, ini_orbital, ini_v_nonlocal
     integer :: nodes_expected, nodes_count, nodes_last
     real(double) :: eigen_min, eigen_mid, eigen_max
     real(double) :: step, ul_extreme1, ul_extreme2, ul_extreme_diff
     real(double) :: result, norm

     !! Constant parameters
     real(double), parameter :: four_pi = 4.0_double * 3.1415926535897932_double

     ini_ul = 0           ! If spin dn allowed: gl_no_orbitals * gl_points_mesh
     ini_orbital = 0      ! If spin dn allowed: gl_no_orbitals
     ini_v_nonlocal = 0   ! If spin dn allowed: gl_no_orbitals * gl_psp_in%psp_comp

     do i=1, gl_points_mesh
        gl_tmp1(i) = 0.0_double
     end do

     ! Initialise non-local pseudopotential
     do j=1, gl_psp_in%psp_comp
        if( gl_psp_in%l_nonlocal(j) == gl_orbitals(orbital)%l ) then
          do i=1, gl_points_mesh
             shift = ini_v_nonlocal + (orbital - 1) * gl_psp_in%psp_comp + j

             ! Note that the factor 4.0 (2.0 from projector and 2.0 from integral) makes units consistent
             ! This is the nonlocal potential. gl_v_nonlocal is a projector 
             ! (includes inverse of Kleinman-Bylander integral)
             gl_tmp1(i) = gl_tmp1(i) &
                        + 4.0_double * (gl_r(i)**3) &
                        * gl_v_nonlocal_int_in(shift) &
                        * gl_v_nonlocal(i + (j-1) * gl_points_mesh) &
                        * gl_psp_in%sign_nonlocal(j)
          end do
        end if
     end do


     ! Initialise wavefunction
     shift = ini_ul + (orbital-1) * gl_points_mesh
     l1 = gl_orbitals(orbital)%l + 1
     do i=1,6
        gl_ul(shift + i) = gl_norm(ini_orbital + orbital) * ( gl_r(i)**l1 )
        gl_tmp2(i) = l1 * gl_ul(shift + i)          ! Zeta
        gl_tmp3(i) = l1 * l1 * gl_ul(shift + i)     ! Derivative of zeta
     end do

     ! Restrict eigenvalue range
     nodes_expected = gl_orbitals(orbital)%n - gl_orbitals(orbital)%l
     nodes_count = 0

     eigen_max = 0.0_double
     call solve_schroedinger(eigen_max, orbital, nodes_last, nodes_count)
     do while (nodes_count < nodes_expected)
        eigen_max =  eigen_max + 1.0_double
        call solve_schroedinger(eigen_max, orbital, nodes_last, nodes_count)
     end do

     eigen_min = eigen_max - 1.0_double
     call solve_schroedinger(eigen_min, orbital, nodes_last, nodes_count)
     do while (nodes_count > nodes_expected - 1.0_double)
        eigen_min =  eigen_min - 1.0_double
        call solve_schroedinger(eigen_min, orbital, nodes_last, nodes_count)
     end do

     eigen_mid = 0.5_double * (eigen_max + eigen_min)

     do while ( (eigen_max - eigen_min) > 1.0e-15_double * abs(eigen_max))
        call solve_schroedinger(eigen_mid, orbital, nodes_last, nodes_count)

        if(nodes_count >= nodes_expected) then
           eigen_max = eigen_mid
        else
           eigen_min = eigen_mid
        end if

        eigen_mid = 0.5_double * (eigen_max + eigen_min)
     end do

     !! Make sure that the last node is at cutoff
     call solve_schroedinger(eigen_mid, orbital, nodes_last, nodes_count)
     if(nodes_count < nodes_expected) then
        eigen_mid = eigen_max
        call solve_schroedinger(eigen_mid, orbital, nodes_last, nodes_count)
     end if


     ! If still not enough nodes, warn about precission
     ! If number of nodes is correct, make sure that 
     !   last node is at last grid point
     if(nodes_count /= nodes_expected) then
        write(*,*) 'WARNING: Expected number of nodes = ', nodes_expected
        write(*,*) '         Obtained number of nodes = ', nodes_count
        write(*,*) '         Last node at point = ', nodes_last
        write(*,*) '         This may be due to lack of precission'
     else
        i=0
        shift = ini_ul + (orbital-1) * gl_points_mesh
        ul_extreme1 = gl_ul(shift + gl_points_mesh)
        step = 1.0e-10_double * abs(eigen_mid)

        do while (abs(ul_extreme1) > 1.0e-12_double .and. i < 20)
           call solve_schroedinger(eigen_mid+step, orbital, nodes_last, nodes_count)
           ul_extreme2 = gl_ul(shift + gl_points_mesh)
           ul_extreme_diff = ul_extreme2 - ul_extreme1

           if(abs(ul_extreme_diff) < 0.1_double) then
              eigen_mid = eigen_mid - step*ul_extreme1/ul_extreme_diff
              call solve_schroedinger(eigen_mid,orbital,nodes_last,nodes_count)
              ul_extreme1 = gl_ul(shift + gl_points_mesh)
           else
              ul_extreme1 = 0.0_double
           end if
           i = i + 1
        end do
     end if

     ! Tail
     call solve_schroedinger(eigen_mid,orbital,nodes_last,nodes_count)
     do i=nodes_last, gl_points_mesh
        gl_ul(shift + i) = 0.0_double
     end do

     ! Back to Hartrees
     gl_eigenvalues( ini_orbital + orbital ) = eigen_mid/2.0_double

     ! Normalisation
     do i=1, gl_points_mesh
        ! Beware! Can use this because solve_schroedinger not called anymore
        gl_tmp1(i) = four_pi * gl_ul(shift + i) * gl_ul(shift + i)
     end do
     call finite_differentiation(gl_tmp1, gl_dtmp1)
     call integral_cubic(gl_r, gl_tmp1, gl_dtmp1, 1, gl_points_mesh, result)
     norm = 1.0_double / sqrt(result)
     do i=1, gl_points_mesh
        ! Back to hartrees
        gl_ul(shift + i) = gl_ul(shift +i) * norm
     end do
     gl_norm(ini_orbital + orbital) = gl_norm(ini_orbital + orbital) * norm

     ! Integrals of nonlocal pseudopotential
     do j=0, gl_psp_in%psp_comp-1
        if( gl_psp_in%l_nonlocal(j+1) == gl_orbitals(orbital)%l ) then
           do i=1, gl_points_mesh
              gl_tmp1(i) = gl_r(i) * gl_v_nonlocal( i + j*gl_points_mesh ) &
                         * gl_ul( shift + i )
           end do
           call simpson_integral_prod(gl_r, gl_tmp1, gl_points_mesh, gl_dnu, &
                gl_v_nonlocal_int_out( ini_v_nonlocal + &
                (orbital-1) * gl_psp_in%psp_comp + j + 1) )
        else
           gl_v_nonlocal_int_out(ini_v_nonlocal + &
               (orbital-1) * gl_psp_in%psp_comp + j + 1) =  0.0_double
        end if
     end do

  end subroutine solve_orbital
!!***

! -----------------------------------------------------------------------------
! Subroutine solve_schroedinger
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_solve/solve_schroedinger *
!!
!!NAME
!! solve_schroedinger
!!USAGE
!!
!!PURPOSE
!! Actuaally solve Schroedinger and count number of nodes
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 15/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine solve_schroedinger(eigen, orbital, nodes_last, nodes_count)

     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types

     implicit none

     !! Passed variables
     real(double) :: eigen
     integer :: orbital
     integer :: nodes_last
     integer :: nodes_count

     !! Local variables
     integer :: i, shift
     integer :: ini_v      ! If spin dn allowed: gl_points_mesh
     integer :: ini_ul
     integer :: nodes_expected

     real(double) :: l

     ! Local temporals
     real(double) :: lc_tmp1, lc_tmp2, lc_dtmp2, lc_tmp3, lc_dtmp3

     !! Parameters
     real(double), parameter :: k19 = 19.0_double
     real(double), parameter :: k27 = 27.0_double
     real(double), parameter :: k106 = 106.0_double
     real(double), parameter :: k251 = 251.0_double
     real(double), parameter :: k264 = 264.0_double
     real(double), parameter :: k475 = 475.0_double
     real(double), parameter :: k502 = 502.0_double
     real(double), parameter :: k646 = 646.0_double
     real(double), parameter :: k720 = 720.0_double
     real(double), parameter :: k1274 = 1274.0_double
     real(double), parameter :: k1901 = 1901.0_double
     real(double), parameter :: k2616 = 2616.0_double
     real(double), parameter :: k2774 = 2774.0_double


     ini_v = 0      ! If spin dn allowed: gl_points_mesh
     ini_ul = 0     ! If spin dn allowed: gl_no_orbitals * gl_points_mesh

     shift = ini_ul + (orbital - 1) * gl_points_mesh

     l = gl_orbitals(orbital)%l

     do i=7,gl_points_mesh
        lc_tmp1 = l * ( l + 1.0_double ) &
                + gl_r(i) * gl_r(i) * ( 2.0_double * gl_v( i + ini_v ) - eigen)
        lc_tmp2 = gl_ul(shift + i - 1) &
                + (gl_dnu / k720) &
                * (k1901 * gl_tmp2(i-1) &
                 - k2774 * gl_tmp2(i-2) &
                 + k2616 * gl_tmp2(i-3) &
                 - k1274 * gl_tmp2(i-4) &
                 + k251 * gl_tmp2(i-5))
        lc_dtmp2 = gl_tmp2(i-1) &
                 + (gl_dnu / k720) &
                 * (k1901 * gl_tmp3(i-1) &
                  - k2774 * gl_tmp3(i-2) &
                  + k2616 * gl_tmp3(i-3) &
                  - k1274 * gl_tmp3(i-4) &
                  + k251 * gl_tmp3(i-5))
        gl_tmp3(i) = lc_dtmp2 + lc_tmp1 * lc_tmp2 + gl_tmp1(i)


        lc_tmp3 = gl_ul(shift + i - 1) &
                + (gl_dnu / k720) &
                * (k251 * lc_dtmp2 &
                 + k646 * gl_tmp2(i-1) &
                 - k264 * gl_tmp2(i-2) &
                 + k106 * gl_tmp2(i-3) &
                 - k19 * gl_tmp2(i-4)) 
        lc_dtmp3 = gl_tmp2(i-1) &
                + (gl_dnu / k720) &
                * (k251 * gl_tmp3(i) &
                 + k646 * gl_tmp3(i-1) &
                 - k264 * gl_tmp3(i-2) &
                 + k106 * gl_tmp3(i-3) &
                 - k19 * gl_tmp3(i-4)) 

       gl_ul(shift + i) = (k475 * lc_tmp3 + k27 * lc_tmp2) / k502
       gl_tmp2(i) = (k475 * lc_dtmp3 + k27 * lc_dtmp2) / k502
       gl_tmp3(i) = gl_tmp2(i) + lc_tmp1 * gl_ul(shift + i) + gl_tmp1(i)
    end do

    ! Count nodes
    nodes_expected = gl_orbitals(orbital)%n - gl_orbitals(orbital)%l
    nodes_last = gl_points_mesh + 1
    nodes_count = 0

    do i=2,gl_points_mesh-1
       if ( (gl_ul(shift + i) > 0.0_double .and. gl_ul(shift + i - 1) < 0.0_double) &
       .or. (gl_ul(shift + i) < 0.0_double .and. gl_ul(shift + i - 1) > 0.0_double) &
       .or. (gl_ul(shift + i) == 0.0_double)) then
           nodes_count = nodes_count + 1
           if(nodes_count  == nodes_expected) then
              nodes_last = i
           end if
       end if
    end do

  end subroutine solve_schroedinger
!!***

end module cq_ut_psp2pao_solve


