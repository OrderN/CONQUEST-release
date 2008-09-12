! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_basis
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_basis *
!!
!!NAME
!! cq_ut_psp2pao_basis
!!PURPOSE
!! Calculation of the basis
!!USES
!!  fdf
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_psp
!!  cq_ut_psp2pao_numeric
!!  cq_ut_psp2pao_hartree
!!  cq_ut_psp2pao_xc
!!  cq_ut_psp2pao_alloc
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 08/05/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_basis

  implicit none

contains


! -----------------------------------------------------------------------------
! Subroutine get_radial_function
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_basis/get_radial_function *
!!
!!NAME
!! get_radial_function
!!USAGE
!!
!!PURPOSE
!! Divide Ul by the radius
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 15/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine get_radial_function

    use cq_ut_psp2pao_global

    implicit none

    integer :: i,j,shift

     do i=1, gl_no_orbitals
        do  j=1, gl_points_mesh
            shift = j+(i-1)*gl_points_mesh 
            gl_ul(shift) = gl_ul(shift)/gl_r(j)
        end do
     end do

  end subroutine get_radial_function
!!***


! -----------------------------------------------------------------------------
! Subroutine fill_wf_set
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_basis/fill_wf_set *
!!
!!NAME
!! fill_wf_set
!!USAGE
!!
!!PURPOSE
!! Store relevant wavefunction information for use in basis construction
!!INPUTS
!!
!!USES
!!  fdf
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_psp
!!  cq_ut_psp2pao_numeric
!!  cq_ut_psp2pao_hartree
!!  cq_ut_psp2pao_xc
!!  cq_ut_psp2pao_alloc
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 15/08/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine fill_wf_set(no_set)

!     use input_module, ONLY: io_assign, io_close
     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_psp
     use cq_ut_psp2pao_numeric
     use cq_ut_psp2pao_hartree
     use cq_ut_psp2pao_xc
     use cq_ut_psp2pao_alloc

     implicit none

     ! Passed variables
     integer :: no_set

     ! Local variables

     integer :: i, j, stat

     gl_wf_set(no_set)%no_orbitals = gl_no_orbitals
     gl_wf_set(no_set)%points_mesh = gl_points_mesh
     gl_wf_set(no_set)%partial_core = gl_psp_in%partial_core
     gl_wf_set(no_set)%cutoff_radius = gl_cutoff_radius
     gl_wf_set(no_set)%e_total = gl_e_total
     gl_wf_set(no_set)%core_charge = gl_psp_in%valence_elec
     gl_wf_set(no_set)%dnu = gl_dnu
     gl_wf_set(no_set)%psp_comp = gl_psp_in%psp_comp

     call allocate_wavefunction_set(no_set, gl_wf_set(no_set)%no_orbitals, &
                                    gl_wf_set(no_set)%points_mesh, &
                                    gl_wf_set(no_set)%psp_comp, &
                                    gl_wf_set(no_set)%partial_core)

!     allocate(gl_wf_set(no_set)%r(gl_wf_set(no_set)%points_mesh), stat=STAT)
!     allocate(gl_wf_set(no_set)%rho(2*gl_wf_set(no_set)%points_mesh), stat=STAT)
!     if(gl_wf_set(no_set)%partial_core) then
!        allocate(gl_wf_set(no_set)%rhopc(gl_wf_set(no_set)%points_mesh), stat=STAT)
!        allocate(gl_wf_set(no_set)%rhopc2(gl_wf_set(no_set)%points_mesh), stat=STAT)
!     end if

     do i=1, gl_wf_set(no_set)%points_mesh
        gl_wf_set(no_set)%r(i) = gl_r(i)
        gl_wf_set(no_set)%rho(i) = gl_rho(i)
        gl_wf_set(no_set)%rho(i + gl_wf_set(no_set)%points_mesh) = gl_rho(i + gl_points_mesh)
        if(gl_wf_set(no_set)%partial_core) then
           gl_wf_set(no_set)%rhopc(i) = gl_rhopc(i)
        end if
     end do


!     allocate(gl_wf_set(no_set)%l_nonlocal(gl_wf_set(no_set)%psp_comp), stat=STAT)
!     allocate(gl_wf_set(no_set)%sign_nonlocal(gl_wf_set(no_set)%psp_comp), stat=STAT)
!     allocate(gl_wf_set(no_set)%v_local(gl_wf_set(no_set)%points_mesh), stat=STAT)
!     allocate(gl_wf_set(no_set)%v_nonlocal(gl_wf_set(no_set)%points_mesh*gl_wf_set(no_set)%psp_comp), stat=STAT)
!
!     ! Second derivatives for splines
!     allocate(gl_wf_set(no_set)%v_local2(gl_wf_set(no_set)%points_mesh), stat=STAT)
!     allocate(gl_wf_set(no_set)%v_nonlocal2(gl_wf_set(no_set)%points_mesh*gl_wf_set(no_set)%psp_comp), stat=STAT)

     do i=1, gl_wf_set(no_set)%psp_comp
        gl_wf_set(no_set)%l_nonlocal(i) = gl_psp_in%l_nonlocal(i)
        gl_wf_set(no_set)%sign_nonlocal(i) = gl_psp_in%sign_nonlocal(i)
        gl_wf_set(no_set)%v_nonlocal_cutoff(i) = gl_psp_in%v_nonlocal_cutoff(i)
     end do

     gl_wf_set(no_set)%v_local_cutoff = gl_psp_in%v_local_cutoff

     do i=1, gl_wf_set(no_set)%points_mesh
        gl_wf_set(no_set)%v_local(i) = gl_v_nuclear(i)
        do j=1, gl_wf_set(no_set)%psp_comp
           gl_wf_set(no_set)%v_nonlocal(i+(j-1)*gl_wf_set(no_set)%points_mesh) = &
                                          gl_v_nonlocal(i+(j-1)*gl_wf_set(no_set)%points_mesh)
        end do
     end do

     ! Orbitals
!     allocate(gl_wf_set(no_set)%orb_n(gl_wf_set(no_set)%no_orbitals), stat=STAT)
!     allocate(gl_wf_set(no_set)%orb_l(gl_wf_set(no_set)%no_orbitals), stat=STAT)
!     allocate(gl_wf_set(no_set)%orb_occ(gl_wf_set(no_set)%no_orbitals), stat=STAT)
!     allocate(gl_wf_set(no_set)%orb_keep(gl_wf_set(no_set)%no_orbitals), stat=STAT)
!     allocate(gl_wf_set(no_set)%orb_cutoff_radius(gl_wf_set(no_set)%no_orbitals), stat=STAT)
!     allocate(gl_wf_set(no_set)%orb_eigenvalues(gl_wf_set(no_set)%no_orbitals), stat=STAT)
!     allocate(gl_wf_set(no_set)%orb_ul(gl_wf_set(no_set)%no_orbitals, gl_wf_set(no_set)%points_mesh), stat=STAT)
!     allocate(gl_wf_set(no_set)%orb_ul2(gl_wf_set(no_set)%no_orbitals, gl_wf_set(no_set)%points_mesh), stat=STAT)

     do i=1, gl_wf_set(no_set)%no_orbitals
        gl_wf_set(no_set)%orb_l(i) = gl_orbitals(i)%l
        gl_wf_set(no_set)%orb_n(i) = gl_orbitals(i)%n
        gl_wf_set(no_set)%orb_cutoff_radius(i) = gl_cutoff_radius
        gl_wf_set(no_set)%orb_occ(i) = gl_occ(i)
        gl_wf_set(no_set)%orb_keep(i) = gl_orbitals(i)%keep
        gl_wf_set(no_set)%orb_zeta(i) = gl_orbitals(i)%zeta
        gl_wf_set(no_set)%orb_eigenvalues(i) = gl_eigenvalues(i)
        do  j=1, gl_wf_set(no_set)%points_mesh
            gl_wf_set(no_set)%orb_ul(i, j) = gl_ul(j+(i-1)*gl_wf_set(no_set)%points_mesh)
        end do
     end do

     ! Splines
     call spline(gl_wf_set(no_set)%r, gl_wf_set(no_set)%v_local, &
                 gl_wf_set(no_set)%points_mesh, gl_wf_set(no_set)%v_local2)
     do i=0, gl_wf_set(no_set)%psp_comp-1
        call spline(gl_wf_set(no_set)%r, &
             gl_wf_set(no_set)%v_nonlocal(i*gl_wf_set(no_set)%points_mesh+1:(i+1)*gl_wf_set(no_set)%points_mesh), &
             gl_wf_set(no_set)%points_mesh, &
             gl_wf_set(no_set)%v_nonlocal2(i*gl_wf_set(no_set)%points_mesh+1:(i+1)*gl_wf_set(no_set)%points_mesh))
     end do
     call spline(gl_wf_set(no_set)%r, gl_wf_set(no_set)%v_nonlocal, &
                 gl_wf_set(no_set)%points_mesh, gl_wf_set(no_set)%v_nonlocal2)
     do i=1, gl_wf_set(no_set)%no_orbitals
        call spline(gl_wf_set(no_set)%r, gl_wf_set(no_set)%orb_ul(i,:), &
                    gl_wf_set(no_set)%points_mesh, gl_wf_set(no_set)%orb_ul2(i,:))
     end do
     if(gl_wf_set(no_set)%partial_core) then
        call spline(gl_wf_set(no_set)%r, gl_wf_set(no_set)%rhopc, &
                    gl_wf_set(no_set)%points_mesh, gl_wf_set(no_set)%rhopc2)
     end if

  end subroutine fill_wf_set
!!***

! -----------------------------------------------------------------------------
! Subroutine calculate_basis
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_basis/calculate_basis *
!!
!!NAME
!! calculate_basis
!!USAGE
!!
!!PURPOSE
!! Calculates a basis from a set of orbitals
!!INPUTS
!!  fdf
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_psp
!!  cq_ut_psp2pao_numeric
!!  cq_ut_psp2pao_hartree
!!  cq_ut_psp2pao_xc
!!  cq_ut_psp2pao_alloc
!!USES
!! 
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 08/05/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine calculate_basis(no_sets)

!     use input_module, ONLY: io_assign, io_close
     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_psp
     use cq_ut_psp2pao_numeric
     use cq_ut_psp2pao_hartree
     use cq_ut_psp2pao_xc
     use cq_ut_psp2pao_alloc

     implicit none

     ! Passed variables
     integer :: no_sets

     ! Local variables
     integer :: i, j, k, stat
     integer :: max_set, orb_no
     real(double) :: r_min, r_max


     ! Count number of states in the basis and find maximum and minimum radii
     max_set = -1
     r_min = 1.0e100_double;
     r_max = 0.0_double;
     do i=1,no_sets
        do j=1, gl_wf_set(i)%no_orbitals
           if(gl_wf_set(i)%orb_keep(j)) then
              gl_basis%no_orbitals = gl_basis%no_orbitals + 1
           end if
        end do
        if(gl_basis%no_orbitals > 0) then
           if(r_min > gl_wf_set(i)%r(1)) then
              r_min = gl_wf_set(i)%r(1)
           end if
           if(r_max < gl_wf_set(i)%r(gl_wf_set(i)%points_mesh)) then
              max_set = i
              r_max = gl_wf_set(i)%r(gl_wf_set(i)%points_mesh)
           end if
        end if
     end do


     gl_basis%points_mesh = gl_wf_set(max_set)%points_mesh
     gl_basis%partial_core = gl_wf_set(1)%partial_core
     gl_basis%e_total = gl_wf_set(1)%e_total
     gl_basis%cutoff_radius = gl_wf_set(max_set)%cutoff_radius
     gl_basis%core_charge = gl_wf_set(1)%core_charge
     gl_basis%dnu = gl_wf_set(max_set)%dnu
     gl_basis%psp_comp = gl_wf_set(1)%psp_comp

     call allocate_basis(gl_basis%no_orbitals, gl_basis%points_mesh, gl_basis%psp_comp, gl_basis%partial_core)

     do i=1, gl_basis%points_mesh
        gl_basis%r(i) = gl_wf_set(max_set)%r(i)
        gl_basis%rho(i) = gl_wf_set(1)%rho(i)
        gl_basis%rho(i + gl_basis%points_mesh) = gl_wf_set(1)%rho(i + gl_wf_set(1)%points_mesh)
        if(gl_wf_set(1)%partial_core) then
           gl_basis%rhopc(i) = gl_wf_set(1)%rhopc(i)
        end if
     end do

     do i=1, gl_basis%psp_comp
        gl_basis%l_nonlocal(i) = gl_wf_set(1)%l_nonlocal(i)
        gl_basis%sign_nonlocal(i) = gl_wf_set(1)%sign_nonlocal(i)
        gl_basis%v_nonlocal_cutoff(i) = gl_wf_set(1)%v_nonlocal_cutoff(i)
     end do

     gl_basis%v_local_cutoff = gl_wf_set(1)%v_local_cutoff

     do i=1, gl_basis%points_mesh
        if(gl_basis%r(i) > gl_wf_set(1)%r(gl_wf_set(1)%points_mesh)) then
           gl_basis%v_local(i) = -gl_basis%core_charge/gl_basis%r(i)
           do j=1, gl_basis%psp_comp
              gl_basis%v_nonlocal(i+(j-1)*gl_basis%points_mesh) = 0.0_double 
           end do
        else
           call spline_interpolation(gl_wf_set(1)%r, &
                                     gl_wf_set(1)%v_local, &
                                     gl_wf_set(1)%points_mesh, &
                                     gl_wf_set(1)%v_local2, &
                                     gl_basis%r(i), &
                                     gl_basis%v_local(i))
           do j=1, gl_basis%psp_comp
              call spline_interpolation(gl_wf_set(1)%r, &
                   gl_wf_set(1)%v_nonlocal((j-1)*gl_wf_set(1)%points_mesh+1:j*gl_wf_set(1)%points_mesh), &
                   gl_wf_set(1)%points_mesh, &
                   gl_wf_set(1)%v_nonlocal2((j-1)*gl_wf_set(1)%points_mesh+1:j*gl_wf_set(1)%points_mesh), &
                   gl_basis%r(i), &
                   gl_basis%v_nonlocal((j-1)*gl_basis%points_mesh + i))
           end do
        end if
     end do

     orb_no = 1
     do i=1,no_sets
       do j=1, gl_wf_set(i)%no_orbitals
          if(gl_wf_set(i)%orb_keep(j)) then
             gl_basis%orb_n(orb_no) = gl_wf_set(i)%orb_n(j)
             gl_basis%orb_l(orb_no) = gl_wf_set(i)%orb_l(j)
             gl_basis%orb_keep(orb_no) = gl_wf_set(i)%orb_keep(j)        ! Not really necessary
             gl_basis%orb_zeta(orb_no) = gl_wf_set(i)%orb_zeta(j)
             gl_basis%orb_cutoff_radius(orb_no) = gl_wf_set(i)%orb_cutoff_radius(j)

             ! Only the first set is occupied
             if(i==1) then
                gl_basis%orb_occ(orb_no) = gl_wf_set(i)%orb_occ(j)
             else
                gl_basis%orb_occ(orb_no) = 0.0_double
             end if
             gl_basis%orb_eigenvalues(orb_no) = gl_wf_set(i)%orb_eigenvalues(j)
             do  k=1, gl_basis%points_mesh
                 call spline_interpolation(gl_wf_set(i)%r, &
                                           gl_wf_set(i)%orb_ul(j,:), &
                                           gl_wf_set(i)%points_mesh, &
                                           gl_wf_set(i)%orb_ul2(j,:), &
                                           gl_basis%r(k), gl_basis%orb_ul(orb_no, k))
             end do
             orb_no = orb_no + 1
          end if
       end do
     end do

  end subroutine calculate_basis
!!***

end module cq_ut_psp2pao_basis

