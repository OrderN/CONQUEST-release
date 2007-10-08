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
! Subroutine initialise_orbital_list
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_initialise/initialise_orbital_list *
!!
!!NAME
!! initialise_orbital_list
!!USAGE
!!
!!*ast*
!!
!!PURPOSE
!! Reads selected parts of the input file and pseudopotentials
!! to generate the list of orbitals that will be calculated
!!INPUTS
!!
!!USES
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 26/06/2007
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine initialise_orbital_list(no_wf_sets)

     use fdf
     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_alloc

     ! Passed variables

     integer :: no_wf_sets

     ! Local variables
     character(len=15) :: zeta_method
     character(len=200) :: label
     character(len=200) :: pseudopotential_file
     character(len=200) :: wf_file
     character(len=250) :: line

     logical :: is_multiple_zeta, is_polarised
     integer stat, i, j, k, l, lun, fdf_file
     type(wf_set_description), dimension(:), allocatable :: wf_descr
     real(double) :: ion_charge, cutoff
     real(double) :: ion_charge_def, cutoff_def
     real(double) :: tmp_charge, junk_double
     integer :: no_orbitals
     integer :: tmp_elec, tmp_n, tmp_l, min, shell_elec, core_elec
     type(ang_momentum_list) :: l_table, l_table_kept, l_tmp
     type(ang_momentum_list) :: l_expected, l_achieved
     logical :: read_option
     integer :: general_maximum_zeta, general_maximum_pol_zeta
     integer :: polarisation_level, polarisation_levels_achieved
     integer :: zeta_levels, maximum_zeta, tmp_zeta
     logical :: found     

     ! Constants
     real(double), parameter :: huge_value = 1.0e30_double

     zeta_method=fdf_string('MultipleZetaMethod', 'eigenstates')
     if(zeta_method == 'eigenstates') then
       no_wf_sets = 1
       is_multiple_zeta = .false.
       ion_charge_def = 0.0_double
       cutoff_def = 10.0_double
     else if(zeta_method == 'ions') then
       no_wf_sets = 2 
       is_multiple_zeta = .true.
       ion_charge_def = huge_value    ! Unrealistic default value
       cutoff_def = 10.0_double
     else if(zeta_method == 'cutoffs') then
       no_wf_sets = 2 
       is_multiple_zeta = .true.
       ion_charge_def = 0.0_double
       cutoff_def = huge_value    ! Unrealistic default value
     else
          write(*,*) "Invalid MultipleZetaMethod: '",trim(zeta_method),"'"
          write(*,*) "Valid values are 'eigenstates', 'ions' and 'cutoffs'"
          stop
     end if

     write(*,'(3a)') "Using method '",trim(zeta_method),"'"

     ! Decide number of wavefunction sets from method and user input
     no_wf_sets=fdf_integer('NumberWaveFunctionSets',no_wf_sets)

     if(no_wf_sets <= 0) then
       write(*,*) 'The number of wave functions sets must be a positive integer'
       stop
     else if((zeta_method == 'ions' .or. zeta_method == 'cutoffs') &
             .and. no_wf_sets < 2) then
       write(*,'(3a)') "Error: The number of wave functions sets for method '",&
                       trim(zeta_method),"' must be at least 2"
       stop
     end  if
    
     write(*,'(a,i4)') 'Number of wavefunction sets in the basis = ',no_wf_sets
     write(*,*) ''

     ! Store basic information for each wavefunction set

     allocate(wf_descr(no_wf_sets), STAT=stat)

     gl_orb_list%size_alloc=0

     ! Initialise l table (to annotate how many orbitals per l channel)
     do i=0,l_table_max-1
       l_table%no_l(i)=0
       l_table%occ_l(i)=0
       l_table_kept%no_l(i)=0
       l_table_kept%occ_l(i)=0
       l_expected%no_l(i)=0
       l_expected%occ_l(i)=0
       l_achieved%no_l(i)=0
       l_achieved%occ_l(i)=0
     end do

     do i=1,no_wf_sets
       write (*,'(a,i4)') "WAVEFUNTION SET", i

       ! Initialise the temporal l table, depending on the method
       if(zeta_method == 'eigenvalues') then
         do j=0,l_table_max-1
           l_tmp%no_l(j) = l_table%no_l(j)
         end do
       else
         do j=0,l_table_max-1
           l_tmp%no_l(j) = 0
         end do
       end if

       ! Open the pseudopotential file to get some data
       write(label,'(a,i0)') 'Pseudopotential',i
       pseudopotential_file = fdf_string(trim(label), '')
       if(pseudopotential_file == '') then
          write(*,'(a,i4,a)') 'Error: A name for the pseudopotential file ',i, &
                              ' is required'
          stop
       end if

       write (*,'(2a)') "---------Pseudopotential : ", trim(pseudopotential_file)

       call io_assign(lun)
       open(unit=lun,file=pseudopotential_file,status='old')

       ! Read first line and discard
       read(lun,fmt=*, iostat=stat) line

       ! Read atom number and valence electrons
       read(lun,fmt=*, iostat=stat) wf_descr(i)%zatom, wf_descr(i)%zion
       call io_close(lun)

       ! If the user wants to read the wavefunction set from file,
       !   do it and continue with next set
       write(label,'(a,i0)') 'ReadFile',i
       read_option = fdf_boolean(trim(label), .false.)
       if(read_option) then
         write(wf_file,'(a,i0)') 'wf.dat.',i
         write(label,'(a,i0)') 'WavefunctionFile',i
         wf_file = fdf_string(trim(label), wf_file)

         call io_assign(lun)
         open(unit=lun,file=wf_file,status='old')

         ! Read the number of orbitals from first line
         read(lun, fmt=*) no_orbitals
         tmp_charge = 0.0_double
         do j=1, no_orbitals
            call allocate_orbital_in_table()
            read(lun, fmt=*)  gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%n, &
                              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%l, &
                              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%occ, &
                              junk_double, &
                              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%keep
            gl_orb_list%table(gl_orb_list%size_alloc)%p%wf_set = i
            tmp_charge = tmp_charge &
                       + gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%occ
            ! This will be set properly later
            gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%zeta=0

            gl_orb_list%table(gl_orb_list%size_alloc)%p%read=.true.

            ! Update l table
            l = gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%l
            l_table%no_l(l) = l_table%no_l(l) + 1
            ! Only annotate the number of orbitals occupied, not the no. of electrons
            if(gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%occ .gt. 0.0) then
              l_table%occ_l(l) = l_table%occ_l(l) + 1
            end if

            ! Update l table of kept orbitals
            if(gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%keep .eqv. .true.) then
              l_table_kept%no_l(l) = l_table_kept%no_l(l) + 1
              ! Only annotate the number of orbitals occupied, not the no. of electrons
              if(gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%occ .gt. 0.0) then
                l_table_kept%occ_l(l) = l_table_kept%occ_l(l) + 1
              end if
            end if
         end do

         ! Close file
         call io_close(lun)

         ion_charge = wf_descr(i)%zion - tmp_charge

         write(*, '(2a)') '---------Read from file : ',trim(wf_file)
         write(*, '(a, f8.3)') '---------IonCharge = ', ion_charge
         cycle            ! Go to the next wavefunction set
       end if

       ! Try to get the ion charge either from IonCharge
       !   or, with priority, from OrbitalInformation, if it exists

       write(label,'(a,i0)') 'NumberOrbitals',i
       no_orbitals = fdf_integer(trim(label), 1)
 
       write(label,'(a,i0)') 'OrbitalInformation',i
       if(fdf_block(trim(label),fdf_file)) then
         ! Read n, l and occupancy for each orbital
         ! IMPORTANT If spin is allowed, two separate occupancies should be read 
         !           (NOT IMPLEMENTED YET)
         tmp_charge = 0.0_double
         do j=1,no_orbitals
            call allocate_orbital_in_table()

            read(fdf_file,*) gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%n, &
                gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%l, &
                gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%occ, &
                gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%keep
            gl_orb_list%table(gl_orb_list%size_alloc)%p%wf_set = i

            ! This will be corrected later
            gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%zeta=0

            gl_orb_list%table(gl_orb_list%size_alloc)%p%read=.false.

            tmp_charge = tmp_charge &
                       + gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%occ

            l = gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%l
            l_table%no_l(l) = l_table%no_l(l) + 1
            ! Only annotate the number of orbitals occupied, not the no. of electrons
            if(gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%occ .gt. 0.0) then
              l_table%occ_l(l) = l_table%occ_l(l) + 1
            end if
            ! If the orbital is kept, annotate
            if(gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%keep .eqv. .true.) then
              l_table_kept%no_l(l) = l_table_kept%no_l(l) + 1
              ! Only annotate the number of orbitals occupied, not the no. of electrons
              if(gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%occ .gt. 0.0) then
                l_table_kept%occ_l(l) = l_table_kept%occ_l(l) + 1
              end if
            end if
         end do
         ion_charge = wf_descr(i)%zion - tmp_charge

         ! If IonCharge is defined, warn that it will be ignored
         write(label,'(a,i0)') 'IonCharge',i
         if(fdf_double(trim(label), huge_value) /= huge_value) then
           write (*,*) "WARNING: '",trim(label),&
                       "' ignored; using OrbitalInformation instead"
         end if
         write(*, '(a)') '---------Read from OrbitalInformation'
         write(*, '(a, f8.3)') '---------IonCharge = ', ion_charge
       else 
         ! First, try to get the ion charge from splicit tag
         write(label,'(a,i0)') 'IonCharge',i
         ion_charge=fdf_double(trim(label), ion_charge_def)

         if(ion_charge == huge_value .and. zeta_method == 'ions') then
           ! The 'ions' method requires IonCharge to be defined
           write(*,'(a,i4,a)') &
                'Error: IonCharge is required for wavefunction set ',&
                i," in method 'ions'"
           stop
       end if

     ! Calculate standard electron configuration from atomic number
     tmp_elec = wf_descr(i)%zatom - ion_charge
     tmp_n=1
     tmp_l=0

     if((wf_descr(i)%zion - ion_charge) .lt. 0) then
        write(*,'(a,i4,a)') &
             'Error: The ion charge (wavefunction set ',i,') has exhausted'
        write(*,'(a)') '       the valence electrons and goes into the core'
        stop
     end if

     core_elec = wf_descr(i)%zatom-wf_descr(i)%zion
     do while(tmp_elec > 0)
       ! The Afbau principle imposes iterating twice over the ang. momentum
       do j=0,1
         k=tmp_n+j          ! This is the n value
         do l=tmp_l,0,-1    ! This is the l value
            shell_elec=2+4*l
            if(shell_elec .le. tmp_elec) then
              min=shell_elec
            else
              min=tmp_elec
            end if
            k=k+1
            tmp_elec  = tmp_elec  - shell_elec
            core_elec = core_elec - min

            ! If all core electron have been calculated
            !   the rest are "valence" electrons and their
            !   orbitals are part of the basis
            ! Although it shouldn't be frequent, the following
            !   allows that only part of a shell is valence
            if(core_elec .lt. 0) then
              call allocate_orbital_in_table()


              
              ! In this part of the calculation, only occupied orbitals
              !   are assigned

              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%n = l + l_tmp%no_l(l) + 1
              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%l = l
              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%occ = -core_elec
              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%keep = .true.
              gl_orb_list%table(gl_orb_list%size_alloc)%p%wf_set = i

              ! This will be corrected later
              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%zeta=0

              gl_orb_list%table(gl_orb_list%size_alloc)%p%read=.false.

              ! Update l tables (in this case, all orbitals are kept)
              l_table%no_l(l) = l_table%no_l(l) + 1
              l_table_kept%no_l(l) = l_table_kept%no_l(l) + 1
              ! The following, just in case. Not sure I'll use it
              if(core_elec /= 0) then
                l_table%occ_l(l) = l_table%occ_l(l) + 1
                l_table_kept%occ_l(l) = l_table_kept%occ_l(l) + 1
              end if

              ! Make sure that the next shell starts anew
              core_elec = 0.0_double
            end if
            if(tmp_elec .le. 0) then
              exit
            end if
         end do
         if(tmp_elec .le. 0) then
           exit
         end if
       end do
       ! Move to next shell
       tmp_n=tmp_n+1
       tmp_l=tmp_l+1
     end do
         write(*, '(a)') '---------Calculated using the Afbau principle'
         write(*, '(a, f8.3)') '---------IonCharge = ', ion_charge
       end if
       
       write(label,'(a,i0)') 'CutoffRadius',i
       cutoff = fdf_double(trim(label), cutoff_def)


       ! The 'cutoffs' method requires CutoffRadius to be defined
       ! IMPORTANT: These cutoffs are only informative 
       !           (i.e. no need to store them for later)
       if(zeta_method == 'cutoffs' .and. cutoff == huge_value) then
         write(*,'(a,i4,a)') &
              'Error: CutoffRadius is required for wavefunction set ',&
              i," in method 'cutoffs'"
        stop
       end if
     end do

     ! Up to now, the user's desires as to the zeta and polarisation levels
     !   have been ignored and the orbital list has been built using other
     !   specified tags
     ! Now, it is time to check what the zetas, etc. of the list are and,
     !   if they are not enough, to complete the list with empty orbitals.
     !   These will be assigned to the first wf set, so if the user is not
     !   happy, he/she will have to override this behaviour explicitly

     ! The structure l_table contains information on the zeta/polarisation
     !   per l channel. In a different table, we store the user's
     !   requirements. Use MaximumZetaPerChannel if available, if not
     !   MaximumZeta, and complete with PolarisationLevel

     ! For the 'ions' and 'cutoff' methods, warn the user if 
     !   MultipleZetaBasis is false
     is_multiple_zeta=fdf_boolean('MultipleZetaBasis',is_multiple_zeta)
     if(.not.is_multiple_zeta .and. &
       (zeta_method == 'ions' .or. zeta_method == 'cutoffs')) then
         write(*,*) "WARNING: The '",trim(zeta_method),&
                    "' method usually implies a multiple zeta basis,"
         write(*,*) "         but the input specifies 'MultipleZetaBasis F'"
         write(*,*) "         The user-defined value might be ignored"
     end if

     ! Should the basis be polarised?
     is_polarised=fdf_boolean('PolarisedBasis', .false.)

     ! Default values for the zeta and polarisation levels
     if(is_multiple_zeta .eqv. .true.) then
       general_maximum_zeta = 2
     else
       general_maximum_zeta = 1
     end if

     if(is_polarised .eqv. .true.) then
       polarisation_level = 1
     else
       polarisation_level = 0
     end if
     
     ! User-specified values for zetas
     general_maximum_zeta=fdf_integer('MaximumZeta', general_maximum_zeta)
     general_maximum_pol_zeta=fdf_integer('MaximumPolZeta', 1)
     polarisation_level=fdf_integer('PolarisationLevel', polarisation_level)

     ! Determine the number of l values 
     !   from which there is at least one orbital occupied (and keep is T)
     zeta_levels=0           ! Occupied, as opposed to polarisation levels
     polarisation_levels_achieved = 0
     do i=0,l_table_max-1
       if( (l_table_kept%no_l(i) .ge. 1) .and. (l_table_kept%occ_l(i) .ge. 1)) then
         zeta_levels = zeta_levels + 1
         ! Store the zeta value for that l channel
         l_expected%no_l(i) = general_maximum_zeta 
         l_expected%occ_l(i) = 1                       ! To indicate "zeta level"
       else
         ! Not a "zeta channel"; it could be a "polarisation channel"; check
         l_expected%occ_l(i) = 0     ! To signal a polarisation channel
         if( (l_table_kept%no_l(i) .ge. 1) .and. (l_table_kept%occ_l(i) .eq. 0)) then
           l_expected%no_l(i) = general_maximum_pol_zeta
           if(is_polarised .eqv. .false.) then
             write(*,'(a)') 'WARNING: The basis will be polarised, although not requested'
             is_polarised = .true.
           end if
           polarisation_levels_achieved = polarisation_levels_achieved + 1
         else
           l_expected%no_l(i) = 0    ! This channel doesn't count at all!!!
         end if
       end if
     end do

     ! Create table of expected zeta values for each occupied l
     if(fdf_block('MaximumZetaPerChannel',fdf_file)) then
       write(*,'(a,i2,a)') 'Expected exactly ',zeta_levels,' l channels in block MaximumZetaPerChannel'
       do i=1,zeta_levels
          read(fdf_file,*) l, maximum_zeta
          if( (l_table_kept%no_l(l) .le. 0) .or. (l_table_kept%occ_l(l) .le. 0) ) then
            write(*,'(a,i2,a)') "WARNING: In 'MaximumZetaPerChannel', channel l = ", l, &
                                " user-defined zeta value will be ignored"
            write(*,'(a)'), "         There are no occupied orbitals for that channel" 
          else 
            if(l_table_kept%no_l(l) .le. maximum_zeta) then
              l_expected%no_l(l)=maximum_zeta
            end if
          end if
       end do    
     end if

     ! Create table of expected zeta values for polarisations orbitals
     if(fdf_block('MaximumZetaPerPolChannel',fdf_file) .and. polarisation_level .ge. 1) then
       write(*,'(a,i2,a)') 'Expected exactly ',polarisation_level,' l channels in block MaximumZetaPerPolChannel'
       do i=1,polarisation_level
          read(fdf_file,*) l, maximum_zeta
          if( (l_table_kept%no_l(l) .ge. 1) .or. (l_table_kept%occ_l(l) .ge. 1) ) then
            write(*,'(a,i2,a)') "WARNING: In 'MaximumZetaPerPolChannel', channel l = ", l, &
                                " user-defined zeta value will be ignored"
            write(*,'(a)'), "         This is not a polarisation channel" 
          else 
            if(l_table_kept%no_l(l) .le. maximum_zeta) then
              l_expected%no_l(l)=maximum_zeta
            end if
          end if
       end do    
     end if

     ! Complete orbital table with polarisation levels
     l=0
     do while(polarisation_level .gt. polarisation_levels_achieved)
       if(l_table%no_l(l) == 0 .or. l_table%occ_l(l) == 0) then          ! This is the next polarisation level
         call allocate_orbital_in_table()
         gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%n = l+1         ! Only one polarisation level per l
         gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%l = l
         gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%occ = 0         ! Empty orbital
         gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%keep = .true.
         gl_orb_list%table(gl_orb_list%size_alloc)%p%wf_set = 1          ! Extra orbitals, always set 1
         gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%zeta=0          ! To be assigned later
         gl_orb_list%table(gl_orb_list%size_alloc)%p%read=.false.

         ! Update l tables
         l_table%no_l(l) = l_table%no_l(l) + 1                           ! This would work in general, for the first zeta
         l_table%occ_l(l) = 0
         l_table_kept%no_l(l) = l_table_kept%no_l(l) + 1
         l_table_kept%occ_l(l) = 0

         ! Assign default polarisation zeta value, if it hasn't been
         if(l_expected%no_l(l) .eq. 0) then
           l_expected%no_l(l) = general_maximum_pol_zeta
         end if
     
         polarisation_levels_achieved = polarisation_levels_achieved + 1
       end if
       l = l + 1
     end do

     ! Complete orbital table with missing zetas
     do i=0,l_table_max-1
       if(l_table_kept%no_l(i) .gt. l_expected%no_l(i)) then
         ! The requested zeta level for this l has been surpased
         write(*,'(a,i2,a,i2,a,i2)')  "WARNING: Channel l = ", i, " has ",l_table_kept%no_l(i),&
                                      " zetas, while the request was ",l_expected%no_l(i)
       else if(l_table_kept%no_l(i) .lt. l_expected%no_l(i)) then
         ! Add extra orbitals to the list until completion of requested zeta
         k=i                     ! This works as the n value
         do while(l_table_kept%no_l(i) .lt. l_expected%no_l(i))
           ! Very primitive search in the list of orbitals
           found = .false.
           k = k + 1             ! Go to next n
           do j=1,gl_orb_list%size_alloc
             ! If the n value is in the list and is kept, go to next value
             ! If it is in the list but not kept, change keep to T
             ! If it is not in the list, add it
             if(gl_orb_list%table(j)%p%orb%l == i) then   ! First, check if it's the right l
               if(gl_orb_list%table(j)%p%orb%n == k) then ! Value of n found
                 found = .true.                           ! Pair (n,l) found
                 tmp_zeta = tmp_zeta + 1
                 if(gl_orb_list%table(j)%p%orb%keep .eqv. .true.) then
                   exit                                   ! Go to next value of n
                 else                  ! The orbital will be calculated, so keep it
                   gl_orb_list%table(j)%p%orb%keep = .true.
                   l_table_kept%no_l(i) = l_table_kept%no_l(i) + 1
                   if(gl_orb_list%table(j)%p%orb%occ .gt. 0.0) then
                     l_table_kept%occ_l(i) = l_table_kept%occ_l(i) + 1
                   end if
                   exit
                 end if
               end if
             end if
           end do
           if(found .eqv. .false.) then     ! Pair (n,l) not found; add
              call allocate_orbital_in_table()
              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%n = k
              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%l = i
              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%occ = 0         ! Empty orbital
              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%keep = .true.
              gl_orb_list%table(gl_orb_list%size_alloc)%p%wf_set = 1          ! Extra orbitals, always set 1
              gl_orb_list%table(gl_orb_list%size_alloc)%p%orb%zeta=0          ! To be assigned later
              gl_orb_list%table(gl_orb_list%size_alloc)%p%read=.false.

              ! Update l tables
              l_table%no_l(i) = l_table%no_l(i) + 1
              l_table_kept%no_l(i) = l_table_kept%no_l(i) + 1
           end if
         end do
       end if
     end do

     ! Finally, assing zeta values
     do i=1,gl_orb_list%size_alloc
       l=gl_orb_list%table(i)%p%orb%l
       l_achieved%no_l(l) = l_achieved%no_l(l) + 1
       gl_orb_list%table(i)%p%orb%zeta = l_achieved%no_l(l)
     end do
     
     print *,""
     write(*,'(a)') 'Summary of orbitals'
     print *,""
     write(*,'(a)') ' Orb #   n   l   occupancy   kept   zeta   atomic set'
     write(*,'(a)') '------- --- --- ----------- ------ ------ ------------'
     do i=1,gl_orb_list%size_alloc
       write(*,'("    ",i2,"  ",i2,"  ",i2,"     ",f6.4,"     ",l2,"     ",i2,"        ",i2)') &
                i, gl_orb_list%table(i)%p%orb%n, &
                gl_orb_list%table(i)%p%orb%l, &
                gl_orb_list%table(i)%p%orb%occ, &
                gl_orb_list%table(i)%p%orb%keep, &
                gl_orb_list%table(i)%p%orb%zeta, &
                gl_orb_list%table(i)%p%wf_set
     end do

     if(allocated(wf_descr))  deallocate(wf_descr)

  end subroutine initialise_orbital_list
 
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
!! 26/06/2007 Slight changes to accomodate the new list of orbitals
!!            This makes user's life easier
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

     ! Read cutoff radius
     write(label,'(a,i0)') 'CutoffRadius',wf_set_no
     gl_cutoff_radius = fdf_double(trim(label), 10.0_double)

     ! Read local cutoff tolerance - MUST be the same for all wavefunction sets
     gl_local_cutoff_tolerance = fdf_integer('LocalCutoffTolerance', 6)

     ! Get the number of orbitals in this wavefunction set from the total list
     gl_no_orbitals=0
     do i=1,gl_orb_list%size_alloc
       if(gl_orb_list%table(i)%p%wf_set == wf_set_no) then
         gl_no_orbitals = gl_no_orbitals + 1
       end if
     end do

     !! Read pseudopotential file

     call read_pseudopotential(pseudopotential_file, wf_set_no)

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

       ! Get info for the orbitals from total table
       i=0
       do j=1,gl_orb_list%size_alloc
         if(gl_orb_list%table(j)%p%wf_set == wf_set_no) then
            i = i + 1  
 
            gl_orbitals(i)%n=gl_orb_list%table(j)%p%orb%n
            gl_orbitals(i)%l=gl_orb_list%table(j)%p%orb%l
            gl_orbitals(i)%occ=gl_orb_list%table(j)%p%orb%occ
            gl_orbitals(i)%keep=gl_orb_list%table(j)%p%orb%keep
            gl_orbitals(i)%zeta=gl_orb_list%table(j)%p%orb%zeta

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
               write(*,*) 'Two many electrons (',gl_orbitals(i)%occ,') in orbital: n = ',&
                          gl_orbitals(i)%n, ' l = ', gl_orbitals(i)%l
               stop
            end if


            no_elec = no_elec + gl_orbitals(i)%occ
         end if
       end do


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

