! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_initialise
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_initialise *
!!
!!NAME
!! cq_ut_psp2pao_initialise
!!PURPOSE
!! Initialisation  of the job
!!USES
!!  fdf
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_psp
!!  cq_ut_psp2pao_numeric
!!  cq_ut_psp2pao_hartree
!!  cq_ut_psp2pao_xc
!!  cq_ut_psp2pao_xc
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 23/02/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_initialise

  implicit none

contains

! -----------------------------------------------------------------------------
! Subroutine initialise
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_initialise/initialise *
!!
!!NAME
!! initialise
!!USAGE
!!
!!PURPOSE
!! Reads the instructions for the job from 'input_psp2pao' and the
!! pseudopotential from the specified file
!!INPUTS
!!
!!USES
!!  fdf
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_psp
!!  cq_ut_psp2pao_numeric
!!  cq_ut_psp2pao_hartree
!!  cq_ut_psp2pao_alloc
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 23/02/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine initialise(wf_set_no)

     use fdf
     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_psp
     use cq_ut_psp2pao_numeric
     use cq_ut_psp2pao_hartree
     use cq_ut_psp2pao_xc
     use cq_ut_psp2pao_alloc

     implicit none

     ! Passed variables

     integer :: wf_set_no

     ! Local variables

     character(len=200) :: label
     character(len=20) :: def
     character(len=200) :: pseudopotential_file

     integer :: i, j, stat, fdf_file

     real(double) :: no_elec, norm_factor, r_ini, r_ini2, alpha

     ! Constant parameters
     real(double), parameter :: rho_alpha = 1.0_double
     real(double), parameter :: rho_norm_exp = 1.5_double
     real(double), parameter :: pi = 3.1415926535897932_double
     real(double), parameter :: small = 1.0e-13_double

     def = ' '

     !!! Start reading job specification from 'input_psp2pao'

     !! Read  name of pseudopotential file
     write(label,'(a,i0)') 'Pseudopotential',wf_set_no
     pseudopotential_file = fdf_string(trim(label), def)


     !! Read grid information

     ! Read number of grid points (mesh)
     write(label,'(a,i0)') 'GridPointsMesh',wf_set_no
     gl_points_mesh = fdf_integer(trim(label), 1)
     gl_points_mesh = 2 * (gl_points_mesh / 2) + 1

     ! Read number of grid points (pseudopotential)
     ! This MUST agree with the pseudopotential file
     write(label,'(a,i0)') 'GridPointsPsp',wf_set_no
     gl_points_psp=fdf_integer(trim(label), 1)

     write(label,'(a,i0)') 'GridRMin',wf_set_no
     gl_grid_rmin=fdf_double(trim(label), 0.0_double)
     write(label,'(a,i0)') 'GridRMax',wf_set_no
     gl_grid_rmax=fdf_double(trim(label), 1.0_double)

     ! Read number of orbitals
     write(label,'(a,i0)') 'NumberOrbitals',wf_set_no
     gl_no_orbitals = fdf_integer(trim(label), 1)

     ! Read cutoff radius
     write(label,'(a,i0)') 'CutoffRadius',wf_set_no
     gl_cutoff_radius = fdf_double(trim(label), 10.0_double)

     ! Read local cutoff tolerance - MUST be the same for all wavefunction sets
     gl_local_cutoff_tolerance = fdf_integer('LocalCutoffTolerance', 6)


     !! Read pseudopotential file

     call read_pseudopotential(pseudopotential_file)

     !! Allocate global arrays to be used for output
     call allocate_writable_globals(gl_no_orbitals, gl_points_mesh, &
                                    gl_psp_in%psp_comp, gl_psp_in%partial_core)

     ! Read XC flag
     !  0 = PW92 LDA (Default)
     !  1 = PBE GGA
     write(label,'(a,i0)') 'XCFlag',wf_set_no
     gl_xcflag = fdf_integer(trim(label), 0)

     if(gl_xcflag < 0 .or. gl_xcflag > 1) then
        write(*,*) 'gl_xcflag = ',gl_xcflag,' not permitted'
        stop
     end if

     !! Read orbital information

     no_elec = 0.0_double

     call cq_allocate(gl_norm, 2 * gl_no_orbitals, 'initialise', 'gl_norm')   ! Normalisation
     do i=1, 2*gl_no_orbitals
        gl_norm(i) = 1.0_double
     end do

     write(label,'(a,i0)') 'OrbitalInformation',wf_set_no
     if(fdf_block(trim(label),fdf_file)) then
       ! Read n, l and occupancy for each orbital
       ! IMPORTANT If spin is allowed, two separate occupancies should be read (NOT IMPLEMENTED YET)
       do i=1,gl_no_orbitals
          read(fdf_file,*) gl_orbitals(i)%n,gl_orbitals(i)%l,gl_orbitals(i)%occ,gl_orbitals(i)%keep

          ! Set occupancies

          ! For the moment, NON-POLARISED CASE is assumed
          ! For the polarised case, we should read two occupancies 
          !   for up and down, assign them to the global variable 
          !   and assign the nonlocal psp integrals
          ! Half the occupancy read because unpolarised
          gl_occ(i) = 0.5_double * gl_orbitals(i)%occ
          gl_occ(i + gl_no_orbitals) = 0.5_double * gl_orbitals(i)%occ

          ! Check that the quantum numbers are consistent
          if(gl_orbitals(i)%n <= gl_orbitals(i)%l) then
             write(*,*) 'Incompatible quantum numbers: n = ',gl_orbitals(i)%n, ' l = ', gl_orbitals(i)%l
             stop
          end if

          ! Check that there are not too many electrons
          ! NOTE: When spin is implemented, the first 2.0 should be 1.0
          if(gl_orbitals(i)%occ > 2.0_double * (2.0_double * gl_orbitals(i)%l + 1.0_double) ) then
             write(*,*) 'Two many electrons (',gl_orbitals(i)%occ,') in orbital: n = ',gl_orbitals(i)%n, ' l = ', gl_orbitals(i)%l
             stop
          end if

          no_elec = no_elec + gl_orbitals(i)%occ
       end do
     end if


     !! Read SCF instructions

     ! Read maximum number of cycles
     write(label,'(a,i0)') 'MaxSCFCycles',wf_set_no
     gl_max_scf_cycles = fdf_integer(trim(label), 100)

     ! Read potential mixing
     write(label,'(a,i0)') 'PotentialMix',wf_set_no
     gl_potential_mix = fdf_double(trim(label), 0.7_double)

     ! Read orbital mixing
     write(label,'(a,i0)') 'OrbitalMix',wf_set_no
     gl_orbital_mix = fdf_double(trim(label), 0.7_double)


     !! Scaling
     gl_dnu = log(gl_cutoff_radius/small) / (gl_points_mesh - 1.0_double)
     alpha = exp(gl_dnu)

     if(gl_psp_in%psp_comp > 0) then
        ! Integrals
        call cq_allocate(gl_v_nonlocal_int_in, 2 * gl_psp_in%psp_comp * gl_no_orbitals, &
                         'initialise', 'gl_v_nonlocal_int_in')
        call cq_allocate(gl_v_nonlocal_int_out, 2 * gl_psp_in%psp_comp * gl_no_orbitals, &
                         'initialise', 'gl_v_nonlocal_int_out')
     end if

     call cq_allocate(gl_v_ext, gl_points_mesh, 'initialise','gl_v_ext')
     call cq_allocate(gl_v_hartree, gl_points_mesh, 'initialise','gl_v_hartree')
     call cq_allocate(gl_v_xc, 2 * gl_points_mesh, 'initialise','gl_v_xc')
     call cq_allocate(gl_v, 2 * gl_points_mesh, 'initialise','gl_v')

     ! Temporal vectors and derivatives
     call cq_allocate(gl_tmp1, gl_points_mesh, 'initialise','gl_tmp1')
     call cq_allocate(gl_tmp2, gl_points_mesh, 'initialise','gl_tmp2')
     call cq_allocate(gl_tmp3, gl_points_mesh, 'initialise','gl_tmp3')
     call cq_allocate(gl_dtmp1, gl_points_mesh, 'initialise','gl_dtmp1')
     call cq_allocate(gl_dtmp2, gl_points_mesh, 'initialise','gl_dtmp2')
     call cq_allocate(gl_dtmp3, gl_points_mesh, 'initialise','gl_dtmp3')

     !! Initialisations

     norm_factor = 2.0_double * no_elec * rho_alpha &
                 * ((rho_alpha / pi)**rho_norm_exp) / 3.0_double
     r_ini = small

     do i=1,gl_points_mesh
        gl_r(i) = r_ini

        ! In principle, spin polarised
        r_ini2 = r_ini * r_ini
        gl_rho(i) = 0.5_double*norm_factor*r_ini2*exp(-rho_alpha*r_ini2)
        gl_rho(gl_points_mesh+i) = gl_rho(i)

        if(r_ini > gl_psp_in%r_max) then
            gl_v_nuclear(i) = -gl_psp_in%valence_elec / r_ini

            do j=0,gl_psp_in%psp_comp-1
               gl_v_nonlocal(i + j * gl_points_mesh) = 0.0_double
            end do

            if(gl_psp_in%partial_core) then
               gl_rhopc(i) = 0.0_double
            end if
        else
            ! Result in gl_v_nuclear
            call spline_interpolation(gl_psp_in%r, gl_psp_in%v_nuclear, &
                                      gl_points_psp, gl_psp_in%v_nuclear2, &
                                      r_ini, gl_v_nuclear(i))

            do j=0,gl_psp_in%psp_comp-1

             ! Result in gl_v_nonlocal
             call spline_interpolation(gl_psp_in%r, &
                  gl_psp_in%v_nonlocal(j*gl_points_psp+1:(j+1)*gl_points_psp), &
                  gl_points_psp, &
                  gl_psp_in%v_nonlocal2(j*gl_points_psp+1:(j+1)*gl_points_psp), &
                  r_ini, &
                  gl_v_nonlocal(i+j*gl_points_mesh))
            end do

            if(gl_psp_in%partial_core) then

             ! Result in gl_rhopc
             call spline_interpolation(gl_psp_in%r, &
                  gl_psp_in%rhopc, gl_points_psp, gl_psp_in%rhopc2, &
                  r_ini, gl_rhopc(i))
            end if
        end if

        !*ast* IMP Alternative external potentials should be assigned here
        !*ast*     if (when) their parameters are read from file
        !*ast*     For the moment, just set to zero

        gl_v_ext(i) = 0.0_double

        r_ini = r_ini * alpha
     end do

     ! Integrals of non-local pseudopotential
     if(gl_psp_in%psp_comp > 0) then
        do i=1, 2 * gl_psp_in%psp_comp * gl_no_orbitals
           gl_v_nonlocal_int_in(i) = 0.0_double
           gl_v_nonlocal_int_out(i) = 0.0_double
        end do
     end if

     ! Hartree and exchange-correlation potentials
     call get_v_hartree
     call get_v_xc

     do i=1, gl_points_mesh
        gl_v(i) = gl_v_nuclear(i) + gl_v_hartree(i) + gl_v_xc(i) + gl_v_ext(i)
        gl_v(i + gl_points_mesh) = gl_v(i)
     end do

  end subroutine initialise
!!***

end module cq_ut_psp2pao_initialise

