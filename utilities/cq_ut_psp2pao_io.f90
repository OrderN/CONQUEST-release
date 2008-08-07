! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_io
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_io *
!!
!!NAME
!! cq_ut_psp2pao_io
!!PURPOSE
!! Input/output
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_alloc
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 03/04/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_io

  implicit none

contains

! -----------------------------------------------------------------------------
! Subroutine write_wavefunctions
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_io/write_wavefunctions *
!!
!!NAME
!! write_wavefunctions
!!USAGE
!!
!!PURPOSE
!! Write the wavefunctions to file
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_alloc
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 03/04/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine write_wavefunctions(filename)

     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_alloc

     implicit none

     !! Passed variables

     character(len=*) :: filename	

     !! Local variables

     integer :: i, j, lun

     character(len=200) :: format1

     ! Create file for writing

     call io_assign(lun)
     open(unit=lun,file=filename,status='replace')


     ! Write fields to file

     write(lun, fmt='(i5,i10,l5,2e18.10e2,i5,e18.10e2,i5)') &
                    gl_no_orbitals, gl_points_mesh, &
                    gl_psp_in%partial_core, gl_e_total, &
                    gl_cutoff_radius, &
                    gl_psp_in%valence_elec, gl_dnu, gl_psp_in%psp_comp

     do i=1, gl_no_orbitals
        write(lun, fmt='(2i5,2e18.10e2,l2)') gl_orbitals(i)%n, gl_orbitals(i)%l, &
                                             gl_occ(i)+gl_occ(i+gl_no_orbitals), &
                                             gl_eigenvalues(i), gl_orbitals(i)%keep
     end do

     write(lun, fmt=*) gl_psp_in%l_nonlocal(1:gl_psp_in%psp_comp)
     write(lun, fmt=*) gl_psp_in%sign_nonlocal(1:gl_psp_in%psp_comp)
 
     do i=1,gl_points_mesh
        if(gl_psp_in%partial_core) then
          write(format1, '(a,i3,a)') '(',4+gl_no_orbitals+gl_psp_in%psp_comp,'e18.10e2)'
          write(lun, fmt=format1) gl_r(i), &
                    (gl_ul(i+j*gl_points_mesh)/gl_r(i), j=0,gl_no_orbitals-1), &
                    gl_rho(i)+gl_rho(i+gl_points_mesh), gl_rhopc(i), &
                    gl_v_nuclear(i), &
                    (gl_v_nonlocal(i+j*gl_points_mesh),j=0,gl_psp_in%psp_comp-1)
        else
          write(format1, '(a,i3,a)') '(',3+gl_no_orbitals+gl_psp_in%psp_comp,'e18.10e2)'
          write(lun, fmt=format1) gl_r(i), &
                    (gl_ul(i+j*gl_points_mesh)/gl_r(i), j=0,gl_no_orbitals-1), &
                    gl_rho(i)+gl_rho(i+gl_points_mesh), &
                    gl_v_nuclear(i), &
                    (gl_v_nonlocal(i+j*gl_points_mesh),j=0,gl_psp_in%psp_comp-1)
        end if
     end do

     ! Close file

     call io_close(lun)

  end subroutine write_wavefunctions
!!***


! -----------------------------------------------------------------------------
! Subroutine read_plato_wavefunctions
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_io/read_plato_wavefunctions *
!!
!!NAME
!! read_plato_wavefunctions
!!USAGE
!!
!!PURPOSE
!! Read wavefunctions from a Plato file (for testing purposes only)
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_alloc
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 08/05/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine read_plato_wavefunctions(filename)

     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_alloc

     implicit none

     !! Passed variables

     character(len=*) :: filename	

     !! Local variables

     integer :: i, j, lun

     integer :: stat   ! Needed only for reallocation of gl_v_nonlocal

     ! Open file for reading

     call io_assign(lun)
     open(unit=lun,file=filename,status='old')

     ! Read fields from file

     read(lun, fmt=*) &
                    gl_no_orbitals, gl_points_mesh, &
                    gl_psp_in%partial_core, gl_e_total, &
                    gl_cutoff_radius, &
                    gl_psp_in%valence_elec, gl_dnu


     ! NOTE that here we don't know gl_psp_in%psp_comp yet
     !      We bypass the problem by assigning a value (1)
     !      And then deallocating and reallocating, once
     !      we now the real one
     !      This is needed because this format is from Plato, 
     !      not from Conquest.

     gl_psp_in%psp_comp = 1   ! Only until it's read from the file
     call allocate_writable_globals(gl_no_orbitals, gl_points_mesh, &
                                    gl_psp_in%psp_comp, gl_psp_in%partial_core)


     do i=1, gl_no_orbitals
        read(lun, fmt=*) gl_orbitals(i)%n, gl_orbitals(i)%l, &
                         gl_occ(i), gl_eigenvalues(i)
     end do

     read(lun, fmt=*) gl_psp_in%psp_comp, &
                      gl_psp_in%l_nonlocal(1:gl_psp_in%psp_comp)

     ! The following conditional is not ellegant
     ! See comment above
     if(gl_psp_in%psp_comp > 0) then
        nullify(gl_v_nonlocal)
        allocate(gl_v_nonlocal(gl_psp_in%psp_comp * gl_points_mesh), STAT=stat)
     end if


     read(lun, fmt=*) gl_psp_in%sign_nonlocal(1:gl_psp_in%psp_comp)

     do i=1,gl_points_mesh
        if(gl_psp_in%partial_core) then
          read(lun, fmt=*) gl_r(i), &
                    (gl_ul(i+j*gl_points_mesh), j=0,gl_no_orbitals-1), &
                    gl_rho(i), gl_rhopc(i), &
                    gl_v_nuclear(i), &
                    (gl_v_nonlocal(i+j*gl_points_mesh), j=0,gl_psp_in%psp_comp-1)
        else
          read(lun, fmt=*) gl_r(i), &
                    (gl_ul(i+j*gl_points_mesh), j=0,gl_no_orbitals-1), &
                    gl_rho(i), &
                    gl_v_nuclear(i), &
                    (gl_v_nonlocal(i+j*gl_points_mesh), j=0,gl_psp_in%psp_comp-1)
        end if

        ! This is because we restrict ourselves to the non-polarised case
        gl_rho(i) = gl_rho(i) / 2.0_double
        gl_rho(i+gl_points_mesh) = gl_rho(i)
     end do

     ! Close file

     call io_close(lun)

  end subroutine read_plato_wavefunctions
!!***


! -----------------------------------------------------------------------------
! Subroutine read_wavefunctions
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_io/read_wavefunctions *
!!
!!NAME
!! read_wavefunctions
!!USAGE
!!
!!PURPOSE
!! Read wavefunctions from a file
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_alloc
!!
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 08/05/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine read_wavefunctions(filename)

     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_alloc

     implicit none

     !! Passed variables

     character(len=*) :: filename	

     !! Local variables

     integer :: i, j, lun
     integer :: stat   ! Needed only for reallocation of gl_v_nonlocal
     real(double) :: tmpcore

     ! Open file for reading

     call io_assign(lun)
     open(unit=lun,file=filename,status='old')

     ! Read fields from file

     read(lun, fmt=*) &
                    gl_no_orbitals, gl_points_mesh, &
                    gl_psp_in%partial_core, gl_e_total, &
                    gl_cutoff_radius, &
                    gl_psp_in%valence_elec, gl_dnu, gl_psp_in%psp_comp

     call allocate_writable_globals(gl_no_orbitals, gl_points_mesh, &
                                    gl_psp_in%psp_comp, gl_psp_in%partial_core)

     if(gl_psp_in%partial_core) then
        tmpcore = 1.0_double
     else
        tmpcore = 0.0_double
     end if

     call allocate_psp_globals(gl_points_mesh, gl_psp_in%psp_comp, tmpcore)

     do i=1, gl_no_orbitals
        read(lun, fmt=*) gl_orbitals(i)%n, gl_orbitals(i)%l, &
                         gl_occ(i), gl_eigenvalues(i), gl_orbitals(i)%keep
        ! For the moment, NON-POLARISED CASE is assumed
        ! For the polarised case, we should read two occupancies 
        !   for up and down, assign them to the global variable 
        !   and assign the nonlocal psp integrals
        ! Half the occupancy read because unpolarised
        gl_occ(i) = 0.5_double * gl_occ(i)
        gl_occ(i + gl_no_orbitals) = gl_occ(i)
     end do

     read(lun, fmt=*) gl_psp_in%l_nonlocal(1:gl_psp_in%psp_comp)
     read(lun, fmt=*) gl_psp_in%sign_nonlocal(1:gl_psp_in%psp_comp)

     do i=1,gl_points_mesh
        if(gl_psp_in%partial_core) then
          read(lun, fmt=*) gl_r(i), &
                    (gl_ul(i+j*gl_points_mesh), j=0,gl_no_orbitals-1), &
                    gl_rho(i), gl_rhopc(i), &
                    gl_v_nuclear(i), &
                    (gl_v_nonlocal(i+j*gl_points_mesh), j=0,gl_psp_in%psp_comp-1)
        else
          read(lun, fmt=*) gl_r(i), &
                    (gl_ul(i+j*gl_points_mesh), j=0,gl_no_orbitals-1), &
                    gl_rho(i), &
                    gl_v_nuclear(i), &
                    (gl_v_nonlocal(i+j*gl_points_mesh), j=0,gl_psp_in%psp_comp-1)
        end if

        ! This is because we restrict ourselves to the non-polarised case
        gl_rho(i) = gl_rho(i) / 2.0_double
        gl_rho(i+gl_points_mesh) = gl_rho(i)
     end do

     ! Close file

     call io_close(lun)

  end subroutine read_wavefunctions
!!***


! -----------------------------------------------------------------------------
! Subroutine write_basis_plato
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_io/write_basis_plato *
!!
!!NAME
!! write_basis
!!USAGE
!!
!!PURPOSE
!! Write the basis to file (in Plato style)
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 17/05/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine write_basis_plato(filename)

     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types

     implicit none

     !! Passed variables

     character(len=*) :: filename	

     !! Local variables

     integer :: i, j, lun

     character(len=200) :: format1


     ! Create file for writing

     call io_assign(lun)
     open(unit=lun,file=filename,status='replace')


     ! Write fields to file

     write(lun, fmt='(i5,i10,l5,5e18.10e2)') &
                    gl_basis%no_orbitals, gl_basis%points_mesh, &
                    gl_basis%partial_core, &
                    gl_basis%e_total, &
                    gl_basis%cutoff_radius, &
                    gl_basis%core_charge, gl_basis%dnu

     do i=1, gl_basis%no_orbitals
        write(lun, fmt='(2i5,3e18.10e2)') gl_basis%orb_n(i), &
                                          gl_basis%orb_l(i), &
                                          gl_basis%orb_occ(i), &
                                          gl_basis%orb_eigenvalues(i), &
                                          gl_basis%orb_cutoff_radius(i)
     end do

     write(lun, fmt=*) gl_basis%psp_comp, &
                       gl_basis%l_nonlocal(1:gl_basis%psp_comp)
     write(lun, fmt=*) gl_basis%sign_nonlocal(1:gl_basis%psp_comp)
 
     do i=1,gl_basis%points_mesh
        if(gl_basis%partial_core) then
          write(format1, '(a,i3,a)') '(',4+gl_basis%no_orbitals+gl_basis%psp_comp,'e18.10e2)'
          write(lun, fmt=format1) gl_basis%r(i), &
                    (gl_basis%orb_ul(j, i), j=1,gl_basis%no_orbitals), &
                    gl_basis%rho(i)+gl_basis%rho(i+gl_basis%points_mesh), &
                    gl_basis%rhopc(i), &
                    gl_basis%v_local(i), &
                    (gl_basis%v_nonlocal(i+j*gl_basis%points_mesh), &
                    j=0,gl_basis%psp_comp-1)
        else
          write(format1, '(a,i3,a)') '(',3+gl_basis%no_orbitals+gl_basis%psp_comp,'e18.10e2)'
          write(lun, fmt=format1) gl_basis%r(i), &
                    (gl_basis%orb_ul(j, i), j=1,gl_basis%no_orbitals), &
                    gl_basis%rho(i)+gl_basis%rho(i+gl_basis%points_mesh), &
                    gl_basis%v_local(i), &
                    (gl_basis%v_nonlocal(i+j*gl_basis%points_mesh), &
                    j=0,gl_basis%psp_comp-1)
        end if
     end do

     call io_close(lun)

  end subroutine write_basis_plato
!!***

! -----------------------------------------------------------------------------
! Subroutine write_basis
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_io/write_basis *
!!
!!NAME
!! write_basis
!!USAGE
!!
!!PURPOSE
!! Write the basis to file (in Conquest -pseudo-siesta- style)
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 17/05/2006
!!MODIFICATION HISTORY
!! 20/09/2007 Clean code and fix wrong infinite for r=0 and l>1
!!SOURCE
!!
  subroutine write_basis(filename)

     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_numeric

     implicit none

     !! Passed variables

     character(len=*) :: filename	

     !! Parameters
     real(double), parameter :: four_pi = 4.0_double * 3.1415926535897932_double

     !! Local variables

     integer :: i, j, lun
     integer :: max_l_orb, max_l_psp
integer :: max_point


     character(len=200) :: format1
     real(double) :: splineder(gl_basis%points_mesh)  ! Second derivatives for splines     
     real(double) :: delta, splineval
real(double), dimension(3) :: primer, primer_der

     ! Create file for writing

     call io_assign(lun)
     open(unit=lun,file=filename,status='replace')


     ! Write fields to file

     ! Preamble, empty (to be ignored by main code)
     write(lun, '(a)') '<preamble>'
     write(lun, '(a)') '</preamble>'

     ! Header
     !*ast* TO DO Fix symbols, mass, self-energy
     write(lun, fmt='(a2,23x,"# Symbol")') 'XX'
     write(lun, fmt='(a20,5x,"# Label")') 'Label               '
     write(lun, fmt='(i5,20x,"# Atomic number")')  0
     write(lun, fmt='(i5,20x,"# Valence charge")')  int(gl_basis%core_charge)
     write(lun, fmt='(g25.15,"# Mass")')  0.000000_double
     write(lun, fmt='(g25.15,"# Self energy")')  0.000000_double
     max_l_orb=0
     max_l_psp=0
     do i=1, gl_basis%no_orbitals
        if(gl_basis%orb_l(i) > max_l_orb)   max_l_orb=gl_basis%orb_l(i)
     end do
     do i=1, gl_basis%psp_comp
        if(gl_basis%l_nonlocal(i) > max_l_psp)   max_l_psp=gl_basis%l_nonlocal(i)
     end do

     write(lun, fmt='(2i4,17x,"# Lmax for basis, no. of nl orbitals")')   max_l_orb, gl_basis%no_orbitals
     write(lun, fmt='(2i4,17x,"# Lmax for projectors, no. of nl KB projectors")')   max_l_psp, gl_basis%psp_comp

     write(lun, fmt='(a)')   "#PAOs:_______________________________"
     do i=1, gl_basis%no_orbitals
        write(lun, fmt='(4i5,f10.6,3x,"#orbital l, n, z, is_polarized, population")') &
             gl_basis%orb_l(i), gl_basis%orb_n(i), gl_basis%orb_zeta(i), 0, 2.0_double*gl_basis%orb_occ(i)
        delta = gl_basis%orb_cutoff_radius(i)/(gl_basis%points_mesh-1)

        do j=1, gl_basis%points_mesh
           if(gl_basis%r(j) < 1e-300_double) then
              gl_basis%orb_ul(i, j) = 0.0_double
           else
              gl_basis%orb_ul(i, j) = ((four_pi)**0.5)*gl_basis%orb_ul(i, j)/(gl_basis%r(j)**(gl_basis%orb_l(i)))
           end if
        end do

        call spline(gl_basis%r, gl_basis%orb_ul(i, 1:gl_basis%points_mesh), &
                    gl_basis%points_mesh, splineder)
        write(lun, fmt='(i4,2g25.15,3x,"# npts, delta, cutoff")')  &
              gl_basis%points_mesh, delta, gl_basis%orb_cutoff_radius(i)
        do j=1, gl_basis%points_mesh
           call spline_interpolation(gl_basis%r, gl_basis%orb_ul(i, 1:gl_basis%points_mesh), &
                                     gl_basis%points_mesh, splineder, (j-1)*delta, splineval)
           write(lun, fmt='(2g25.15)')  (j-1)*delta, splineval !gl_basis%r(j), gl_basis%orb_ul(i, j)
        end do
     end do

     write(lun, fmt='(a)')   "#KBs:________________________________"
     do i=1, gl_basis%psp_comp
        ! Since the energy is included in the projectors, only the sign 
        ! has to be given here as the scaling factor.
        ! This could be changed if needed.
        ! n is just the number in sequence of zetas (for the moment, just one)
        write(lun, fmt='(2i3,x,g23.16,3x,"#kb l, n, Reference energy")') &
             gl_basis%l_nonlocal(i), 1, real(gl_basis%sign_nonlocal(i))
        delta = gl_basis%v_nonlocal_cutoff(i)/(gl_basis%points_mesh-1)

        ! Do proper scaling for Conquest
        do j=1, gl_basis%points_mesh
           if(gl_basis%r(j) < 1e-300_double) then
              gl_basis%v_nonlocal((i-1)*gl_basis%points_mesh+j) = 0.0_double
           else
              gl_basis%v_nonlocal((i-1)*gl_basis%points_mesh+j) = 2.0_double &
                                                          * gl_basis%v_nonlocal((i-1)*gl_basis%points_mesh + j) &
                                                          / (gl_basis%r(j)**(gl_basis%l_nonlocal(i)))
           end if
        end do

        call spline(gl_basis%r, gl_basis%v_nonlocal(1+(i-1)*gl_basis%points_mesh:i*gl_basis%points_mesh), &
                    gl_basis%points_mesh, splineder)
        write(lun, fmt='(i4,2g25.15,3x,"# npts, delta, cutoff")')  &
              gl_basis%points_mesh, delta, gl_basis%v_nonlocal_cutoff(1)

!if(gl_basis%points_mesh < 4) then
!  max_point=gl_basis%points_mesh
!else
!  max_point=4
!end if
!
!do j=2,max_point
!   call spline_interpolation(gl_basis%r, &
!                             gl_basis%v_nonlocal(1+(i-1)*gl_basis%points_mesh:i*gl_basis%points_mesh), &
!                             gl_basis%points_mesh, splineder, (j-1)*delta, primer(j-1))
!end do
!call spline(gl_basis%r(2:4), primer, 3, primer_der)
!call spline_interpolation(gl_basis%r(2:4), primer, 3, primer_der, 0.0_double, splineval)
!!write(lun, fmt='(2g25.15)')  0.0_double, splineval


        ! First point of the projector, extrapolated linearly to prevent  
        !   non-sensical values that happen for l>0 on the original logarithmic scale 
        call spline_interpolation(gl_basis%r, &
                                  gl_basis%v_nonlocal(1+(i-1)*gl_basis%points_mesh:i*gl_basis%points_mesh), &
                                  gl_basis%points_mesh, splineder, 1*delta, primer(1))
        call spline_interpolation(gl_basis%r, &
                                  gl_basis%v_nonlocal(1+(i-1)*gl_basis%points_mesh:i*gl_basis%points_mesh), &
                                  gl_basis%points_mesh, splineder, 2*delta, primer(2))
        !Linear extrapolation to 0.0
        write(lun, fmt='(2g25.15)')  0.0_double, (primer(1)*2-primer(2))

        do j=2, gl_basis%points_mesh
           call spline_interpolation(gl_basis%r, &
                                     gl_basis%v_nonlocal(1+(i-1)*gl_basis%points_mesh:i*gl_basis%points_mesh), &
                                     gl_basis%points_mesh, splineder, (j-1)*delta, splineval)
           write(lun, fmt='(2g25.15)')  (j-1)*delta, splineval
        end do
     end do

     ! This potential is the local part
     ! In the siesta version, this is Vna. Here, it is the local potential
     write(lun, fmt='(a)')   "#Vlocal:________________________________"
     delta = gl_basis%v_local_cutoff/(gl_basis%points_mesh-1)
     write(lun, fmt='(i4,2g25.15,3x,"# npts, delta, cutoff")')  &
              gl_basis%points_mesh, delta, gl_basis%v_local_cutoff
     call spline(gl_basis%r, gl_basis%v_local, gl_basis%points_mesh, splineder)
     do i=1, gl_basis%points_mesh
        call spline_interpolation(gl_basis%r, gl_basis%v_local, &
                                  gl_basis%points_mesh, splineder, (i-1)*delta, splineval)
        write(lun, fmt='(2g25.15)')  (i-1)*delta, splineval
     end do

     write(lun, fmt='(a)')   "#Chlocal:____________________________"
     delta = gl_basis%orb_cutoff_radius(1)/(gl_basis%points_mesh-1)
     do i=1, gl_basis%points_mesh
        gl_basis%rho(i) = gl_basis%rho(i)+gl_basis%rho(i+gl_basis%points_mesh)
     end do
     call spline(gl_basis%r, gl_basis%rho, gl_basis%points_mesh, splineder)
     write(lun, fmt='(i4,2g25.15,3x,"# npts, delta, cutoff")')  &
              gl_basis%points_mesh, delta, gl_basis%orb_cutoff_radius(1)

     do i=1, gl_basis%points_mesh
        call spline_interpolation(gl_basis%r, gl_basis%rho, &
                                  gl_basis%points_mesh, splineder, (i-1)*delta, splineval)
        write(lun, fmt='(2g25.15)')  (i-1)*delta, -splineval
     end do

     if(gl_basis%partial_core) then
       write(lun, fmt='(a)') "#Partial core:_______________________"
       delta = gl_basis%orb_cutoff_radius(1)/(gl_basis%points_mesh-1)
       call spline(gl_basis%r, gl_basis%rhopc, gl_basis%points_mesh, splineder)
       write(lun, fmt='(i4,2g25.15,3x,"# npts, delta, cutoff")')  &
                gl_basis%points_mesh, delta, gl_basis%orb_cutoff_radius(1)
       do i=1, gl_basis%points_mesh
          call spline_interpolation(gl_basis%r, gl_basis%rhopc, &
                                    gl_basis%points_mesh, splineder, (i-1)*delta, splineval)
          write(lun, fmt='(2g25.15)')  (i-1)*delta, splineval
       end do
     end if

     call io_close(lun)

  end subroutine write_basis
!!***


end module cq_ut_psp2pao_io





