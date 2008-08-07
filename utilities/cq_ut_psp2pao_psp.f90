! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_psp
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_psp *
!!
!!NAME
!! cq_ut_psp2pao_psp
!!PURPOSE
!! Reads and handles the pseudopotential
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_numeric
!!  cq_ut_psp2pao_alloc
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 23/02/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_psp

   implicit none

contains

! -----------------------------------------------------------------------------
! Subroutine read_pseudopotential
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_psp/read_pseudopotential *
!!
!!NAME
!! read_pseudopotential
!!USAGE
!!
!!PURPOSE
!! Reads the pseudopotential from the specified file
!!   and builds tables for later use
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_numeric
!!  cq_ut_psp2pao_alloc
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 23/02/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
   subroutine read_pseudopotential(filename,wf_no)

     use cq_ut_psp2pao_global, ONLY: gl_grid_rmin, gl_grid_rmax, gl_points_psp, &
                                     gl_psp_in, gl_local_cutoff_tolerance, gl_cutoff_radius
     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_numeric, ONLY: spline, spline_interpolation, &
                                      simpson_integral_prod, very_small, zero, big
     use cq_ut_psp2pao_alloc

     implicit none

     real(double), parameter :: four_pi = 4.0_double * 3.1415926535897932_double

     ! Passed variables

     character(len=*) :: filename
     integer :: wf_no

     ! Local variables

     integer:: lun, i, j, k

     character(len=250) :: line
     integer :: stat

     integer :: junk_int
     real(double) :: junk_real, z, zat
     real(double) :: r_incr, r
     real(double), dimension(:), allocatable :: tmp, tmp1, tmp2, tmp3, tmp4, tmp5
     real(double) :: tmp_v, tmp_u
     real(double) :: tail, tail_nl
     real(double) :: min_cutoff
!*ast* integer :: tmp_elec, tmp_n, tmp_l, min, shell_elec

     integer :: psp_comp
     integer :: local_comp         ! Local component of the psp
     real(double) :: partial_core  ! > 0.0 if there is a partial core
     real(double) :: valence_elec  ! Number of valence electrons

     type(psp_in_table_fhi), dimension(:), pointer :: readin_psp


     ! Constant parameters

     real(double), parameter :: maxtail = 1.0e-10_double


     !! Read the pseudopotential from a file

     call io_assign(lun)
     open(unit=lun,file=filename,status='old')

     ! NOTE: Fritz-Haber-Institut Troullier-Martins-type Abinit type 6 psp assumed

     ! Read first line and discard
     read(lun,fmt=*, iostat=stat) line
     !read(lun,fmt=*, iostat=stat) line

     ! Read atom number and valence electrons
     read(lun,fmt=*, iostat=stat) zat, z, junk_int

     ! Read the local component of the pseudopotential
     read(lun, fmt=*) junk_int, junk_int, &
                      junk_int, local_comp

     ! Change to indexes starting at 1
     local_comp = local_comp + 1

     ! Read whether there is a partial core
     read(lun, fmt=*) junk_real, partial_core

     ! Read the next 3 lines and discard
     do i=1,3
        read(lun,fmt=*, iostat=stat) line
     end do

     ! Read the number of valence electrons and the number of psp components
     read(lun, fmt=*, iostat=stat) valence_elec, psp_comp

     ! Discard the next 10 lines
     do i=1,10
        read(lun,fmt=*, iostat=stat) line
     end do

     call cq_allocate(readin_psp, psp_comp, 'read_pseudopotential', 'readin_psp')

     ! Allocate memory for the pseudopotential

     ! The actual pseudopotential is read in here
     do i=1,psp_comp
        read(lun,fmt=*, iostat=stat)  readin_psp(i)%nmesh, &
                                      readin_psp(i)%amesh
        readin_psp(i)%lmesh = log(readin_psp(i)%amesh)

        allocate(readin_psp(i)%r(readin_psp(i)%nmesh), STAT=stat)
        allocate(readin_psp(i)%u(readin_psp(i)%nmesh), STAT=stat)
        allocate(readin_psp(i)%u2(readin_psp(i)%nmesh), STAT=stat)
        allocate(readin_psp(i)%v(readin_psp(i)%nmesh), STAT=stat)
        allocate(readin_psp(i)%v2(readin_psp(i)%nmesh), STAT=stat)

        do j=1, readin_psp(i)%nmesh
           read(lun,fmt=*, iostat=stat) &
               junk_int, readin_psp(i)%r(j), readin_psp(i)%u(j), readin_psp(i)%v(j)
        end do

     end do

     ! If there is a partial core, read it
     if(partial_core > 0.0_double) then
   
        ! Allocate partial core
        allocate(readin_psp(1)%rpc(readin_psp(1)%nmesh))
        allocate(readin_psp(1)%rhopc(readin_psp(1)%nmesh))
        allocate(readin_psp(1)%rhopc2(readin_psp(1)%nmesh))
   
        do i=1,readin_psp(1)%nmesh
           read(lun, fmt=*) readin_psp(1)%rpc(i), &
                            readin_psp(1)%rhopc(i)
           readin_psp(1)%rhopc(i) = readin_psp(1)%rhopc(i) / four_pi
        end do
   
        call spline(readin_psp(1)%rpc, readin_psp(1)%rhopc, &
                    readin_psp(1)%nmesh, readin_psp(1)%rhopc2)
     end if

     call io_close(lun)

     !! Kleinman-Bylander transformation

     do i=1,psp_comp
        if(i /= local_comp) then
            allocate(tmp(readin_psp(i)%nmesh), STAT=stat)
            min_cutoff = big
            do j=1,readin_psp(i)%nmesh
               readin_psp(i)%v(j) = readin_psp(i)%v(j) &
                                  - readin_psp(local_comp)%v(j)
               tmp(j) = readin_psp(i)%u(j)*readin_psp(i)%v(j)*readin_psp(i)%u(j)
            
               ! Determine the cutoff for the nonlocal projectors
               !  from the difference v_nonlocal - v_local
               ! In principle, it should be zero, but it's not good 
               !  to check idendities of reals
               if(abs(readin_psp(i)%v(j)) < very_small) then
                  if( min_cutoff > readin_psp(i)%r(j) ) then
                      min_cutoff = readin_psp(i)%r(j)

!                     if(j-1 > 1) then
!                        readin_psp(i)%cutoff = readin_psp(i)%r(j-1)
!                     else
!                        readin_psp(i)%cutoff = zero
!                     end if
                  end if
               else
                  min_cutoff = big      ! This prevents to take accidental crossings as the cutoff
               end if

            end do

            ! Fix cutoff for the (nonlocal) component
            if(min_cutoff < big) then
               readin_psp(i)%cutoff = min_cutoff
            else
               readin_psp(i)%cutoff = readin_psp(i)%r(readin_psp(i)%nmesh)
               write(*,*) 'WARNING: Cutoff tolerance not reached for the nonlocal component', i 
               write(*,*) '         Cutoff set to PAOs cutoff' 
            end if
            if(readin_psp(i)%cutoff > gl_cutoff_radius) then
               write(*,*) 'ERROR: Cutoff for nonlocal component', i, '(', readin_psp(i)%cutoff,')'
               write(*,*) '       is larger than the PAOs cutoff (', gl_cutoff_radius,')'
               write(*,*) '       CutoffRadius should be increased'
               stop
            end if
            

            ! The result of the integral is returned in readin_psp%kb_integral(i)
            call simpson_integral_prod (readin_psp(i)%r, tmp, &
                                        readin_psp(i)%nmesh, readin_psp(i)%lmesh, &
                                        readin_psp(i)%kb_integral)

            do j=1,readin_psp(i)%nmesh
               readin_psp(i)%u(j) = readin_psp(i)%u(j) / readin_psp(i)%r(j)
            end do

            if(allocated(tmp)) deallocate(tmp)
        else        ! Find cutoff for the local component
            min_cutoff = big
            do j=1,readin_psp(i)%nmesh
               if(readin_psp(i)%r(j) > zero) then
                  if( abs( readin_psp(i)%v(j) + ( valence_elec / readin_psp(i)%r(j) ) ) &
                    < 10.0**(-gl_local_cutoff_tolerance) ) then
                      if( min_cutoff > readin_psp(i)%r(j) ) then
                          min_cutoff = readin_psp(i)%r(j)
                      end if
                  else
                    min_cutoff = big
                  end if
               end if
            end do
            if(min_cutoff < big) then
               readin_psp(i)%cutoff = min_cutoff
            else
               readin_psp(i)%cutoff = readin_psp(i)%r(readin_psp(i)%nmesh)
               write(*,*) 'WARNING: Cutoff tolerance not reached for the local component' 
               write(*,*) '         Local cutoff set to PAOs cutoff' 
            end if
            if(readin_psp(i)%cutoff > gl_cutoff_radius) then
               write(*,*) 'ERROR: Cutoff for local component (', readin_psp(i)%cutoff,')'
               write(*,*) '       is larger than the PAOs cutoff (', gl_cutoff_radius,')'
               write(*,*) '       CutoffRadius should be increased (or LocalCutoffTolerance reduced)'
               stop
            end if
        end if

        call spline(readin_psp(i)%r, readin_psp(i)%u, &
                    readin_psp(i)%nmesh, readin_psp(i)%u2)
!        call spline(readin_psp(i)%r, readin_psp(i)%v, &
!                    readin_psp(i)%nmesh, readin_psp(i)%v2)
        call spline(readin_psp(i)%r, readin_psp(i)%v, &
                    readin_psp(i)%nmesh, readin_psp(i)%v2,YP1=0.0_double, &
                    YPN=z/(readin_psp(i)%r(readin_psp(i)%nmesh)*readin_psp(i)%r(readin_psp(i)%nmesh)))
     end do


     if(wf_no == 1) then
       ! Inform the user about the cutoffs
       write (*,*)
       write (*,'(a)') 'The following pseudopotential cutoffs will be assigned'
       write (*,*)
       do i=1,psp_comp
          if(i /= local_comp) then
             write (*,'(a,i3,a,f9.6)') 'Component', i, ' (nonlocal) Cutoff = ', readin_psp(i)%cutoff
          else
             write (*,'(a,i3,a,f9.6)') 'Component', i, ' (local)    Cutoff = ', readin_psp(i)%cutoff
          end if
       end do
       write (*,*)
     end if


     !! Build PSP tables
     r_incr = (gl_grid_rmax/gl_grid_rmin)**(1.0_double/(gl_points_psp-1.0_double))

     gl_psp_in%valence_elec = valence_elec
     gl_psp_in%psp_comp = psp_comp - 1

     call allocate_psp_globals(gl_points_psp, psp_comp, partial_core)

     r = gl_grid_rmin
     do i=1,gl_points_psp
        gl_psp_in%r(i) = r

        ! Result stored in gl_psp_in%v_nuclear(i)
        call spline_interpolation(readin_psp(local_comp)%r,&
                                  readin_psp(local_comp)%v, &
                                  readin_psp(local_comp)%nmesh, &
                                  readin_psp(local_comp)%v2, &
                                  r, gl_psp_in%v_nuclear(i))
        r = r * r_incr
     end do

     ! Result stored in gl_psp_in%v_nuclear2
     call spline(gl_psp_in%r, gl_psp_in%v_nuclear, &
                 gl_points_psp, gl_psp_in%v_nuclear2)

     if(gl_psp_in%r(gl_points_psp) &
        > readin_psp(local_comp)%r(readin_psp(local_comp)%nmesh)) then
       write(*,*) 'Grid range should be reduced: Larger than pseudopotential range'
       stop
     end if

     j=1
     do i=1,psp_comp
        if(i /= local_comp) then
           gl_psp_in%l_nonlocal(j) = i-1
           if(readin_psp(i)%kb_integral < 0.0_double) then
              gl_psp_in%sign_nonlocal(j) = -1
           else
              gl_psp_in%sign_nonlocal(j) = 1
           end if

           allocate(tmp1(readin_psp(i)%nmesh), stat=STAT)
           allocate(tmp2(readin_psp(i)%nmesh), stat=STAT)
           allocate(tmp3(readin_psp(i)%nmesh), stat=STAT)
           allocate(tmp4(readin_psp(i)%nmesh), stat=STAT)
           allocate(tmp5(readin_psp(i)%nmesh), stat=STAT)
           do k=1,readin_psp(i)%nmesh
            tmp1(k)=readin_psp(i)%r(k)
            tmp2(k)=readin_psp(i)%v(k)
            tmp3(k)=readin_psp(i)%v2(k)
            tmp4(k)=readin_psp(i)%u(k)
            tmp5(k)=readin_psp(i)%u2(k)
           end do

           do k=1,gl_points_psp
              ! Result in tmp_v and tmp_u
              call spline_interpolation(tmp1, tmp2, &
                                        readin_psp(i)%nmesh, &
                                        tmp3, &
                                        gl_psp_in%r(k), &
                                        tmp_v)
              call spline_interpolation(tmp1, tmp4, &
                                        readin_psp(i)%nmesh, &
                                        tmp5, &
                                        gl_psp_in%r(k), &
                                        tmp_u)
   
              ! VERY IMPORTANT. The factor 2.0 is needed for consistency of the units
              !                 It seems that the formula works in Rydbergs,
              !                 but I don't know why.
              !                 The potential is in Hartrees, though
      
              gl_psp_in%v_nonlocal((j-1)*gl_points_psp + k) = &
                    tmp_u * tmp_v / &
                    abs(2.0_double * readin_psp(i)%kb_integral)**0.5_double
           end do

           ! Set nonlocal cutoff   
           gl_psp_in%v_nonlocal_cutoff(j) = readin_psp(i)%cutoff

           j = j + 1
   
           if(allocated(tmp1)) deallocate(tmp1)
           if(allocated(tmp2)) deallocate(tmp2)
           if(allocated(tmp3)) deallocate(tmp3)
           if(allocated(tmp4)) deallocate(tmp4)
           if(allocated(tmp5)) deallocate(tmp5)
        else
           ! Set local cutoff
           gl_psp_in%v_local_cutoff = readin_psp(i)%cutoff
        end if
     end do

     do i=0,gl_psp_in%psp_comp-1
        call spline(gl_psp_in%r, &
                 gl_psp_in%v_nonlocal(i*gl_points_psp+1:(i+1)*gl_points_psp), &
                 gl_points_psp, &
                 gl_psp_in%v_nonlocal2(i*gl_points_psp+1:(i+1)*gl_points_psp))
     end do
   
     if(partial_core > 0.0_double) then
        gl_psp_in%partial_core = .true.
   
        do i=1,gl_points_psp
            ! Result in gl_psp_in.rhopc(i)
            call spline_interpolation(readin_psp(1)%rpc, readin_psp(1)%rhopc, &
                                      readin_psp(1)%nmesh, readin_psp(1)%rhopc2, &
                                      gl_psp_in%r(i), gl_psp_in%rhopc(i))
        end do
   
        call spline(gl_psp_in%r, gl_psp_in%rhopc, gl_points_psp, gl_psp_in%rhopc2)
     end if
   
!     tail = 2.0_double*abs(gl_psp_in%v_nuclear(gl_points_psp) &
!              * gl_psp_in%r(gl_points_psp) &
!              + 1.0_double * gl_psp_in%valence_elec)
!   
!     if(tail > maxtail) then
!        write(*,*) 'WARNING: Tail error is too large = ', tail
!        write(*,*) '         Try increasing the cutoff of the pseudopotential'
!     end if
!   
!     tail_nl = 0.0_double
!     do i=1,gl_psp_in%psp_comp
!        if(abs(gl_psp_in%v_nonlocal(i*gl_points_psp)) > tail_nl) then
!           tail_nl = abs(gl_psp_in%v_nonlocal(i*gl_points_psp))
!        end if
!     end do
!   
!     tail_nl = 2.0_double * tail_nl
!   
!     if(tail_nl > maxtail) then
!        write(*,*) 'WARNING: The tail of at least one nonlocal potential is too large = ', tail_nl
!        write(*,*) '         Try increasing the cutoff'
!     end if
   
     gl_psp_in%r_max = gl_psp_in%r(gl_points_psp)

     do i=1,psp_comp
       if(associated(readin_psp(i)%r))    deallocate(readin_psp(i)%r)
       if(associated(readin_psp(i)%u))    deallocate(readin_psp(i)%u)
       if(associated(readin_psp(i)%u2))   deallocate(readin_psp(i)%u2)
       if(associated(readin_psp(i)%v))    deallocate(readin_psp(i)%v)
       if(associated(readin_psp(i)%v2))   deallocate(readin_psp(i)%v2)
     end do  

     if(associated(readin_psp))   deallocate(readin_psp)

   end subroutine read_pseudopotential
!!***

end module cq_ut_psp2pao_psp


