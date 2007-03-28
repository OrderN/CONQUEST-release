! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! cq_ut_psp2pao: PAO-basis generator from pseudopotentials
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao *
!!
!!NAME
!! cq_ut_psp2pao
!!PURPOSE
!! Main routine of the PAO basis generator
!!USES
!!  fdf
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_initialise
!!  cq_ut_psp2pao_scf
!!  cq_ut_psp2pao_density
!!  cq_ut_psp2pao_io
!!  cq_ut_psp2pao_basis
!!  cq_ut_psp2pao_orbital
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_alloc
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 23/02/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
program cq_ut_psp2pao

   use fdf
   use cq_ut_psp2pao_types
   use cq_ut_psp2pao_initialise
   use cq_ut_psp2pao_scf, ONLY : scf_loop
   use cq_ut_psp2pao_density
   use cq_ut_psp2pao_io
   use cq_ut_psp2pao_basis
   use cq_ut_psp2pao_orbital
   use cq_ut_psp2pao_global
   use cq_ut_psp2pao_alloc

   implicit none

   integer :: i, j
   integer :: no_wf_sets
   integer :: no_set
   logical :: read_option
   character(len=200) :: label
   character(len=200) :: wf_file
   real(double) :: sigma

   call fdf_init('input_psp2pao', 'fdf.out')

   no_wf_sets = fdf_integer('NumberWavefunctionSets', 1)
   if(no_wf_sets > 0) then
     call cq_allocate(gl_wf_set, no_wf_sets, 'cq_ut_psp2pao_main', 'gl_wf_set')
   else
     write(*,*) 'The number of wave functions sets must be a positive integer'
     stop
   end  if

   do no_set=1, no_wf_sets

      ! 'ReadFile' = 0 for a calculation of the wavefunctions
      !            = 1 for reading the wavefunctions from a file 
      write(label,'(a,i0)') 'ReadFile',no_set
      read_option = fdf_boolean(trim(label), .false.)

      if(.not.read_option) then
        call initialise(no_set)
        call scf_loop

        write(wf_file,'(a,i0)') 'wavefunctions.txt.',no_set
        write(label,'(a,i0)') 'WavefunctionFile',no_set
        wf_file = fdf_string(trim(label), wf_file)
        write(*,*)
        write(*,*) 'Writing wavefunction file: ', wf_file
        call write_wavefunctions(wf_file)

        call get_radial_function
        call fill_wf_set(no_set)
        call deallocate_writable_globals(gl_psp_in%psp_comp, gl_psp_in%partial_core)
        call deallocate_nonwritable_globals
        call deallocate_psp_globals(gl_psp_in%partial_core)
      else
        write(wf_file,'(a,i0)') 'wavefunctions.txt.',no_set
        write(label,'(a,i0)') 'WavefunctionFile',no_set
        wf_file = fdf_string(trim(label), wf_file)
        write(*,*)
        write(*,*) 'Reading wavefunction set ', no_set, ' from file ', trim(wf_file)
        call read_wavefunctions(wf_file)
        call fill_wf_set(no_set)
        call deallocate_writable_globals(gl_psp_in%psp_comp, gl_psp_in%partial_core)
        call deallocate_psp_globals(gl_psp_in%partial_core)
      end if
   end do

   call calculate_basis(no_wf_sets)

   if(fdf_boolean('SmoothBasis', .false.)) then
      sigma = fdf_double('SmoothSigma', 1.5_double)
      call smooth_basis(sigma)
   end if

   if(fdf_boolean('OrthogonaliseBasis', .true.)) then
      call orthogonalise
   end if

   call get_rho_basis

   wf_file = fdf_string("BasisFile", "basis.txt")
   write(*,*)
   write(*,*) 'Writing basis file: ', wf_file
   call write_basis(wf_file)

   call fdf_shutdown

   call deallocate_wavefunction_sets(no_wf_sets)
   call deallocate_basis(gl_basis%partial_core)

end program cq_ut_psp2pao
!!***
