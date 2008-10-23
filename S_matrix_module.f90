! $Id$
! -----------------------------------------------------------
! Module S_matrix_module
! -----------------------------------------------------------
! Code area 3: Operators
! -----------------------------------------------------------

!!****h* Conquest/S_matrix_module
!!  NAME
!!   S_matrix_module
!!  PURPOSE
!!   Collects together the routines needed to get the S matrix
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/05/2001
!!  MODIFICATION HISTORY
!!   16:03, 04/02/2003 drb 
!!    Created get_onsite_S (mainly taken from kinetic energy)
!!   14:26, 26/02/2003 drb 
!!    Corrected get_S_matrix
!!   07:52, 2003/09/22 dave
!!    Added flags to choose between blips and PAOs as basis set
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2008/02/01 17:53 dave
!!    Changes for output to file not stdout
!!   2008/05/22 ast
!!    Added timers
!!  SOURCE
!!
module S_matrix_module

  use datatypes
  use global_module, ONLY: io_lun
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_smatrix

  implicit none

  real(double) :: InvSTolerance

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
!!***

contains

! -----------------------------------------------------------
! Subroutine get_S_matrix
! -----------------------------------------------------------

!!****f* S_matrix_module/get_S_matrix *
!!
!!  NAME 
!!   get_S_matrix
!!  USAGE
!! 
!!  PURPOSE
!!   Gets a new S matrix by doing blip_to_support, integrating
!!   and updating inverse S matrix
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/05/2001
!!  MODIFICATION HISTORY
!!   08:35, 2003/02/05 dave
!!    Added lines and call to get_onsite_S to replace onsite elements with analytical values
!!   14:25, 26/02/2003 drb
!!    Removed - breaking code for some reason
!!   07:53, 2003/09/22 dave
!!    Added choice between blips and PAOs
!!   08:18, 2003/10/01 dave
!!    A little more for PAOs
!!   12:52, 31/10/2003 drb 
!!    Changed call to assemble
!!   08:31, 2004/07/23 dave
!!    Added call to allow calculation of derivative of S wrt PAO coefficient
!!   2008/05/22 ast
!!    Added timer
!!  TODO
!!    Find out why on-site evaluation gave problems
!!  SOURCE
!!
  subroutine get_S_matrix(inode, ionode)

    use datatypes
    use global_module, only: iprint_ops, flag_basis_set, blips, PAOs, flag_vary_basis, &
                             ni_in_cell, IPRINT_TIME_THRES1
    use matrix_data, ONLY: Srange
    use mult_module, ONLY: matS, matdS
    use set_bucket_module, ONLY: rem_bucket
    use calc_matrix_elements_module, ONLY: get_matrix_elements_new
    use blip_grid_transform_module, ONLY: blip_to_support_new
    use primary_module , ONLY : bundle
    use build_PAO_matrices, ONLY: assemble_2
    use species_module, ONLY: n_species
    use PAO_grid_transform_module, ONLY: PAO_to_grid
    use functions_on_grid, ONLY: supportfns
    use io_module, ONLY: dump_matrix
    use timer_module, ONLY: cq_timer, start_timer, stop_print_timer, WITH_LEVEL

    implicit none

    ! Passed variables
    integer :: inode, ionode
    
    ! Local variables
    integer :: np, ni, iprim
    type(cq_timer) :: tmr_l_tmp1

    call start_timer(tmr_std_smatrix)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    ! Project support functions onto grid
    if(flag_basis_set==blips) then
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Doing blip-to-support ',supportfns
       call blip_to_support_new(inode-1, supportfns)

       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Doing integration ',supportfns
       ! Integrate
       call get_matrix_elements_new(inode-1,rem_bucket(1),matS,supportfns,supportfns)
       ! Do the onsite elements analytically
       ! iprim=0
       !    if(flag_basis_set==blips) then
       !       do np=1,bundle%groups_on_node
       !          if(bundle%nm_nodgroup(np) > 0) then
       !             do ni=1,bundle%nm_nodgroup(np)
       !                iprim=iprim+1
       !                call get_onsite_S(data_blip(1,1,iprim), data_S(1,1,mat(np,Srange)%onsite(ni)), &
       !                     inode, ionode,nsf_species(bundle%speices(iprim), MAX_N_BLIPS)
       !             end do
       !          end if
       !       end do
       !    end if
    else if(flag_basis_set==PAOs) then
       ! Get S matrix with assemble
       if(flag_vary_basis) then
          if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Calling assemble_2 for S, dS: ',matS,matdS
          call assemble_2(Srange, matS,1,matdS)
          !call dump_matrix("NS",matS,inode)
          !call dump_matrix("NdS",matdS,inode)
       else
          if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Calling assemble_2 for S: ',matS
          call assemble_2(Srange, matS,1)
       end if
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'PAO to grid ',supportfns
       ! Also generate support with a call to PAO_to_grid
       call PAO_to_grid(inode-1,supportfns)
    end if
    !call dump_matrix("NS",matS,inode)    
    ! get the new InvS matrix
    call  Iter_Hott_InvS( iprint_ops, 100, 0.0001_double,ni_in_cell, inode, ionode)
    call stop_print_timer(tmr_l_tmp1,"get_S_matrix",IPRINT_TIME_THRES1)
    call stop_timer(tmr_std_smatrix)
    return
  end subroutine get_S_matrix
!!***

! -----------------------------------------------------------
! Subroutine Iter_Hott_InvS
! -----------------------------------------------------------

!!****f* S_matrix_module/Iter_Hott_InvS *
!!
!!  NAME 
!!   Iter_Hott_InvS
!!  USAGE
!! 
!!  PURPOSE
!!   This subroutine finds the inverse of S by Hotelling's
!!   method - see Numerical Recipes (ed 2), pp49-50.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   16/10/98 DRB
!!  MODIFICATION HISTORY
!!   04/05/00 by D.R.Bowler to use new matrix mults
!!   22/05/2001 dave
!!    Added ROBODoc header, indented
!!   25/05/2001 dave
!!    Included in S_matrix_module, relocated HotInvSmm
!!   11/06/2001 dave
!!    Changed to use GenComms
!!   14:54, 28/08/2003 drb 
!!    Added a logical test for diagonalisation (which doesn't need S^-1)
!!   2004/10/28 drb
!!    Added Told so that we can revert to the best T at exit, and changed omega definition 
!!    so that it's divided by total number of orbitals; added user-set criterion on omega 
!!    tolerance
!!  SOURCE
!!
  subroutine Iter_Hott_InvS(output_level, n_L_iterations, tolerance,n_atoms,&
       inode, ionode)

    use datatypes
    use numbers
    use global_module, ONLY: IPRINT_TIME_THRES1
    use matrix_data, ONLY: Trange, TSrange, mat, Srange
    use mult_module, ONLY: allocate_temp_matrix, free_temp_matrix, store_matrix_value, matrix_scale, matrix_sum, &
         matT, matS, return_matrix_value
    use primary_module, only: bundle
    use GenComms, ONLY: gsum, my_barrier
    use DiagModule, ONLY: diagon
    use species_module, ONLY: nsf_species, species
    use timer_module, ONLY: cq_timer,start_timer,stop_print_timer,WITH_LEVEL

    implicit none

    ! Passed variables
    integer :: output_level, inode, ionode,n_atoms, n_L_iterations
    real(double) ::  tolerance

    ! Local variables
    integer :: n_iterations, nn, np, nb, ist, ip, i,j,nsf1,nsf2
    integer :: matI, matT1, matTM, matTold
    real(double) :: step, tot, eps, x, omega, oldomega, deltaomega, n_orbs
    type(cq_timer) :: tmr_l_tmp1

    matI = allocate_temp_matrix(TSrange,0)
    matT1 = allocate_temp_matrix(Trange,0)
    matTold = allocate_temp_matrix(Trange,0)
    matTM = allocate_temp_matrix(TSrange,0)

    if(diagon) then ! Don't need S^-1 for diagonalisation
       call matrix_scale(zero,matT)
       ip = 1
       nb = 1
       do np = 1,bundle%groups_on_node
          if(bundle%nm_nodgroup(np)>0) then
             do i=1,bundle%nm_nodgroup(np)
                do j = 1, nsf_species(bundle%species(ip))
                   call store_matrix_value(matT,np,i,ip,nb,j,j,one,1)
                enddo
                ip = ip+1
             enddo
          end if ! bundle%nm_nodgroup(np)>0
       enddo
    else
       call start_timer(tmr_l_tmp1,WITH_LEVEL)
       n_orbs = zero
       do i=1,n_atoms
          n_orbs = n_orbs + real(nsf_species(species(i)),double)
       end do
       ! First construct the identity
       if (inode.eq.ionode.and.output_level>=2) write(io_lun,*) 'Zeroing data'
       call matrix_scale(zero,matT1)
       call matrix_scale(zero,matI)
       call matrix_scale(zero,matT)
       call matrix_scale(zero,matTold)
       call matrix_scale(zero,matTM)
       if (inode.eq.ionode.and.output_level>=2) write(io_lun,*) 'Creating I'
       ip = 1
       nb = 1
       do np = 1,bundle%groups_on_node
          if(bundle%nm_nodgroup(np)>0) then
             do i=1,bundle%nm_nodgroup(np)
                do j = 1, nsf_species(bundle%species(ip))
                   call store_matrix_value(matT,np,i,ip,nb,j,j,one,1)
                   call store_matrix_value(matI,np,i,ip,nb,j,j,one,1)
                enddo
                ip = ip+1
             enddo
          end if ! bundle%nm_nodgroup(np)>0
       enddo
       call my_barrier


       ! Construct the initial guess for T as epsilon.S^T, where epsilon is
       ! given as 1/(sum_jk S^2_jk)
       tot = 0.0_double
       ip = 1
       do np = 1,bundle%groups_on_node
          if(bundle%nm_nodgroup(np)>0) then
             do i=1,bundle%nm_nodgroup(np)
                do nb=1,mat(np,Srange)%n_nab(i)
                   ist = mat(np,Srange)%i_acc(i)+nb-1
                   do nsf1 = 1,mat(np,Srange)%ndimi(i)
                      do nsf2 = 1,mat(np,Srange)%ndimj(ist)
                         eps = return_matrix_value(matS,np,i,ip,nb,nsf1,nsf2)
                         tot = tot + eps*eps
                      enddo
                   enddo
                enddo
                ip = ip+1
             enddo
          end if ! bundle%nm_nodgroup(np)>0
       enddo
       call gsum(tot)
       eps = 1.0_double/(tot)
       if(output_level>1.and.inode==ionode) write(io_lun,*) 'Eps, tot: ',eps,tot
       call matrix_scale(zero,matT)
       call matrix_sum(zero,matT,eps,matS)
       call stop_print_timer(tmr_l_tmp1,"inverse S preliminaries",IPRINT_TIME_THRES1)
       ! and evaluate the current value of the functional and its gradient
       deltaomega = zero
       oldomega = zero
       if (inode==ionode.and.output_level>=2) write(io_lun,*) 'Starting loop'
       do n_iterations=1,n_L_iterations
          call start_timer(tmr_l_tmp1,WITH_LEVEL)
          if (inode==ionode.and.output_level>=2) &
               write(io_lun,2) n_iterations
          deltaomega = deltaomega * half
          ! check for convergence
          if(n_iterations<3.or.abs(deltaomega)>tolerance) then
             call HotInvS_mm( matI, matS, matT, matT1, matTM, omega,n_iterations)
             deltaomega = omega - oldomega
             if(inode==ionode.and.output_level>=1) then
                write(io_lun,*) 'Omega is ',omega/n_orbs
                if(omega>zero) write(io_lun,*) 'R is ',sqrt(omega)/n_orbs
                write(io_lun,*) 'deltaomega is ',n_iterations,deltaomega
             endif
             if ( omega>oldomega.and.oldomega/=0.0_double) then
                if(inode==ionode) write(io_lun,*) 'Truncation error reached !'
                call matrix_sum(zero,matT,one,matTold)
                call stop_print_timer(tmr_l_tmp1,"an inverse S iteration",IPRINT_TIME_THRES1)
                exit
             endif
             oldomega = omega 
             call matrix_sum(zero,matTold,one,matT)
             call matrix_sum(zero,matT,one,matT1)
          else
             call stop_print_timer(tmr_l_tmp1,"an inverse S iteration",IPRINT_TIME_THRES1)
             exit  ! Leave the do loop
          endif
          call stop_print_timer(tmr_l_tmp1,"an inverse S iteration",IPRINT_TIME_THRES1)
       end do
       ! If this isn't a good guess, then reset to I
       if((omega/n_orbs)>InvSTolerance) then
          if(inode==ionode) write(io_lun,*) 'Setting InvS to I'
          call matrix_scale(zero,matT)
          ip = 1
          nb = 1
          do np = 1,bundle%groups_on_node
             if(bundle%nm_nodgroup(np)>0) then
                do i=1,bundle%nm_nodgroup(np)
                   do j = 1, nsf_species(bundle%species(ip))
                      call store_matrix_value(matT,np,i,ip,nb,j,j,one,1)
                   enddo
                   ip = ip+1
                enddo
             end if ! bundle%nm_nodgroup(np)>0
          enddo
       endif
    end if! End if NOT diagon
    call free_temp_matrix(matTM)
    call free_temp_matrix(matTold)
    call free_temp_matrix(matT1)
    call free_temp_matrix(matI)

1   format(20x,'Starting functional value: ',f15.7,' a.u.')
2   format(/,20x,'Conjugate Gradients InvS iteration:',i5)
3   format(/,20x,'Functional value reached after ',i5,' InvS iterations: ',&
         /,20x,' Omega: ', f15.7, ' DeltaOmega: ', f15.7)
4   format('InvS is ',4i5,f15.7)
5   format('T0S,A is ',4i5,2f15.7)
    return

  end subroutine Iter_Hott_InvS
!!***

!!****f* S_matrix_module/HotInvS_mm *
!!
!!  NAME 
!! 
!!  USAGE
!! 
!!  PURPOSE
!!   To perform all the required matrix multiplications to produce
!!   an inverse S matrix, which will give "tensorially correct" 
!!   gradients for the density matrix
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   15/10/98
!!  MODIFICATION HISTORY
!!   11/06/2001 dave
!!    Changed to use GenComms
!!   12:59, 16/11/2004 dave & tsuyoshi
!!    Removed TS_TS_T multiplication and replaced omega with Frobenius norm (basically same !)
!!   12:12, 03/02/2006 drb 
!   Reworked for new scheme
!!  SOURCE
!!
  subroutine HotInvS_mm( matI, matS, matT0, matT1, matA, omega,n)

    use datatypes
    use numbers
    use mult_module, ONLY: matrix_scale, matrix_product, matrix_sum, free_temp_matrix, allocate_temp_matrix, &
         mult, T_S_TS, TS_T_T, matrix_product_trace
    use matrix_data, ONLY: TSrange, Trange
    use GenComms, ONLY: gsum, inode, my_barrier

    implicit none

    ! Passed variables
    integer :: matI, matS, matT0, matT1, matA, n
    real( double ) :: omega

    !     Local Variables
    integer :: matT0S, matGrad

    matT0S = allocate_temp_matrix(TSrange,0)
    matGrad = allocate_temp_matrix(Trange,0)

    call matrix_scale(zero,matT0S)
    call matrix_scale(zero,matGrad)
    ! Create T0.S
    call matrix_product(matT0, matS, matT0S, mult( T_S_TS ) )
    ! Now A=I-TS, as a diagnostic
    call my_barrier()
    call matrix_sum(zero,matA,one,matI)
    call matrix_sum(one,matA,-one,matT0S)
    ! Create T0S.T0
    call matrix_product(matT0S, matT0, matGrad, mult( TS_T_T ) )
    ! T1 = 2T0 - Grad
    call matrix_sum(zero,matT1,two,matT0)
    call matrix_sum(one,matT1,-one,matGrad)
    omega = matrix_product_trace(matA,matA)
    call free_temp_matrix(matGrad)
    call free_temp_matrix(matT0S)

  end subroutine HotInvS_mm
!!***

! -----------------------------------------------------------
! Subroutine get_onsite_S
! -----------------------------------------------------------

!!****f* H_matrix_module/get_onsite_S *
!!
!!  NAME 
!!   get_onsite_S
!!  USAGE
!! 
!!  PURPOSE
!!   This routine evalutes the overlay matrix elements for the block
!!   diagonal. This can be done analytically relatively quickly, because the
!!   support functions are represented on the _same_ blip grid, as they are
!!   associated with a single atom. 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler/C.M.Goringe
!!  CREATION DATE
!!   16:03, 04/02/2003 drb 
!!  MODIFICATION HISTORY
!! 
!!  SOURCE
!!
  subroutine get_onsite_S(blip_co, resulting_submatrix, &
       inode, ionode,this_nsf, spec)

    use datatypes
    use numbers
    use GenBlas, ONLY: axpy, copy, scal, gemm
    use blip, ONLY: blip_info, BlipArraySize, OneArraySize, FullArraySize, SupportGridSpacing
    use support_spec_format, ONLY: support_function
    use GenComms, ONLY: cq_abort

    implicit none

    ! Shared Variables
    integer :: this_nsf, spec, inode, ionode

    type(support_function) :: blip_co
    real(double) :: resulting_submatrix(this_nsf,this_nsf)
    ! Local Variables
    integer, parameter :: MAX_D = 3

    real(double) ::  FAC(0:MAX_D), D2FAC(0:MAX_D)

    real(double), allocatable, dimension(:) :: work1, work2, work3, work4, work5, work6

    integer :: dx, dy, dz, offset, l, at, nsf1, stat

    allocate(work1(FullArraySize(spec)*this_nsf),work2(FullArraySize(spec)*this_nsf),work3(FullArraySize(spec)*this_nsf), &
         work4(FullArraySize(spec)*this_nsf),work5(FullArraySize(spec)*this_nsf),work6(FullArraySize(spec)*this_nsf), &
         STAT=stat)
    if(stat/=0) call cq_abort("Error allocating arrays for onsite S elements: ",FullArraySize(spec),this_nsf)
    ! first, we copy the blip functions for this atom onto a cubic grid;
    ! we make this grid 'too big' in order to have a fast routine below.

    FAC(0) = 151.0_double/140.0_double
    FAC(1) = 1191.0_double/2240.0_double
    FAC(2) = 3.0_double/56.0_double
    FAC(3) = 1.0_double/2240.0_double
    D2FAC(0) = 3.0_double/2.0_double
    D2FAC(1) = -9.0_double/32.0_double
    D2FAC(2) = -18.0_double/40.0_double
    D2FAC(3) = -3.0_double/160.0_double

    work1 = zero
    offset = BlipArraySize(spec)+1

    do dx = -BlipArraySize(spec), BlipArraySize(spec)
       do dy = -BlipArraySize(spec), BlipArraySize(spec)
          do dz = -BlipArraySize(spec), BlipArraySize(spec)
             l = blip_info(spec)%blip_number(dx,dy,dz)
             if (l.ne.0) then
                at = (((dz+offset)*OneArraySize(spec) + (dy+offset))*OneArraySize(spec) + &
                     (dx+offset)) * this_nsf
                do nsf1 = 1,this_nsf
                   work1(nsf1+at) = blip_co%supp_func(nsf1)%coefficients(l)
                enddo
             end if
          end do
       end do
    end do

    ! now, for each direction in turn, we need to apply a number of
    ! 'spreading operations'. Do z first... put blip(z) in 2, and
    ! del2blip(z) in 3

    call copy(FullArraySize(spec)*this_nsf,work1,1,work2,1)
    call scal(FullArraySize(spec)*this_nsf,FAC(0),work2,1)
    do dz = 1, MAX_D
       offset = dz * OneArraySize(spec) * OneArraySize(spec) * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dz), &
            work1(1:), 1, work2(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dz), &
            work1(1+offset:), 1, work2(1:), 1 )
    end do

!    call copy(FullArraySize(spec)*this_nsf,work1,1,work3,1)
!    call scal(FullArraySize(spec)*this_nsf,FAC(0),work3,1)
!    do dz = 1, MAX_D
!       offset = dz * OneArraySize(spec) * OneArraySize(spec) * this_nsf
!       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dz), &
!            work1(1:), 1, work3(1+offset:), 1 )
!       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dz), &
!            work1(1+offset:), 1, work3(1:), 1 )
!    end do


    ! now do y : put blip(y).blip(z) in 4,
    ! blip(y).del2blip(z) + del2blip(y).blip(z)  in 5

    call copy(FullArraySize(spec)*this_nsf,work2,1,work4,1)
    call scal(FullArraySize(spec)*this_nsf,FAC(0),work4,1)
    do dy = 1, MAX_D
       offset = dy * OneArraySize(spec) * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dy), &
            work2(1:), 1, work4(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dy), &
            work2(1+offset:), 1, work4(1:), 1 )
    end do

!    call copy(FullArraySize(spec)*this_nsf,work2,1,work5,1)
!    call scal(FullArraySize(spec)*this_nsf,FAC(0),work5,1)
!    do dy = 1, MAX_D
!       offset = dy * OneArraySize(spec) * this_nsf
!       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dy), &
!            work2(1:), 1, work5(1+offset:), 1 )
!       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dy), &
!            work2(1+offset:), 1, work5(1:), 1 )
!    end do

!    call axpy(FullArraySize(spec)*this_nsf,FAC(0),work3,1,work5,1)
!    do dy = 1, MAX_D
!       offset = dy * OneArraySize(spec) * this_nsf
!       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dy), &
!            work3(1:), 1, work5(1+offset:), 1 )
!       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dy), &
!            work3(1+offset:), 1, work5(1:), 1 )
!    end do
!
!    ! and x - put it all into 6
!
!    call copy(FullArraySize(spec)*this_nsf,work5,1,work6,1)
!    call scal(FullArraySize(spec)*this_nsf,FAC(0),work6,1)
!    do dx = 1, MAX_D
!       offset = dx * this_nsf
!       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dx), &
!            work5(1:), 1, work6(1+offset:), 1 )
!       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dx), &
!            work5(1+offset:), 1, work6(1:), 1 )
!    end do
    work6 = zero
    call axpy(FullArraySize(spec)*this_nsf,FAC(0),work4,1,work6,1)
    do dx = 1, MAX_D
       offset = dx * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dx), &
            work4(1:), 1, work6(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dx), &
            work4(1+offset:), 1, work6(1:), 1 )
    end do

    ! and now get the matrix elements by multiplication...

    call gemm('n','t',this_nsf,this_nsf,OneArraySize(spec)*OneArraySize(spec)*OneArraySize(spec), &
         one,work1,this_nsf,work6,this_nsf,zero,resulting_submatrix,this_nsf )
    call scal(this_nsf*this_nsf,SupportGridSpacing(spec),resulting_submatrix,1)
    deallocate(work1,work2,work3, work4,work5,work6, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating arrays for onsite S blip elements: ",FullArraySize(spec),this_nsf)
    return
  end subroutine get_onsite_S
!!***

!%%!  subroutine basic_get_onsite_S(blip_co, resulting_submatrix, &
!%%!       inode, ionode,this_nsf, MAX_N_BLIPS)
!%%!
!%%!    use datatypes
!%%!    use numbers
!%%!    use GenBlas
!%%!    use blip, ONLY: blip_info(spec)%blip_number, BlipArraySize(spec)
!%%!    use dimens, ONLY: support_grid_spacing
!%%!
!%%!    implicit none
!%%!
!%%!    ! Passed Variables
!%%!    integer :: this_nsf, MAX_N_BLIPS, inode, ionode
!%%!    real(double) :: blip_co(this_nsf,MAX_N_BLIPS), resulting_submatrix(this_nsf,this_nsf)
!%%!
!%%!    ! Local variables
!%%!    real(double) :: FAC(-3:3)
!%%!    real(double) :: work1(this_nsf,-BlipArraySize(spec)-3:BlipArraySize(spec)+3,&
!%%!         -BlipArraySize(spec)-3:BlipArraySize(spec)+3,-BlipArraySize(spec)-3:BlipArraySize(spec)+3), &
!%%!         work2(this_nsf,-BlipArraySize(spec)-3:BlipArraySize(spec)+3,&
!%%!         -BlipArraySize(spec)-3:BlipArraySize(spec)+3,-BlipArraySize(spec)-3:BlipArraySize(spec)+3)
!%%!    integer :: dx,dy,dz,nx,nx1,l,l1,nsf1,ny,ny1,nz,nz1,nsf2
!%%!
!%%!    FAC(-3) = 1.0_double/2240.0_double
!%%!    FAC(-2) = 3.0_double/56.0_double
!%%!    FAC(-1) = 1191.0_double/2240.0_double
!%%!    FAC(0) = 151.0_double/140.0_double
!%%!    FAC(1) = 1191.0_double/2240.0_double
!%%!    FAC(2) = 3.0_double/56.0_double
!%%!    FAC(3) = 1.0_double/2240.0_double
!%%!
!%%!    ! Start by convolving blip coefficients with integrals
!%%!
!%%!    ! x first
!%%!    work1 = zero
!%%!    do dx = -BlipArraySize(spec), BlipArraySize(spec)
!%%!       do dy = -BlipArraySize(spec), BlipArraySize(spec)
!%%!          do dz = -BlipArraySize(spec), BlipArraySize(spec)
!%%!             do nx = -3,3
!%%!                nx1 = nx+dx
!%%!                l = blip_info(spec)%blip_number(dx,dy,dz)
!%%!                if (l/=0) then
!%%!                   do nsf1 = 1,this_nsf
!%%!                      work1(nsf1,dz,dy,nx1) = work1(nsf1,dz,dy,nx1)+FAC(nx)*blip_co(nsf1,l)
!%%!                   enddo
!%%!                end if
!%%!             end do
!%%!          end do
!%%!       end do
!%%!    end do
!%%!    ! Now y
!%%!    work2 = zero
!%%!    do dx = -BlipArraySize(spec)-3, BlipArraySize(spec)+3
!%%!       do dy = -BlipArraySize(spec), BlipArraySize(spec)
!%%!          do dz = -BlipArraySize(spec), BlipArraySize(spec)
!%%!             do ny = -3,3
!%%!                ny1 = ny+dy
!%%!                do nsf1 = 1,this_nsf
!%%!                   work2(nsf1,dz,ny1,dx) = work2(nsf1,dz,ny1,dx)+FAC(ny)*work1(nsf1,dz,dy,dx)
!%%!                enddo
!%%!             end do
!%%!          end do
!%%!       end do
!%%!    end do
!%%!    ! Finally z
!%%!    work1 = zero
!%%!    do dx = -BlipArraySize(spec)-3, BlipArraySize(spec)+3
!%%!       do dy = -BlipArraySize(spec)-3, BlipArraySize(spec)+3
!%%!          do dz = -BlipArraySize(spec), BlipArraySize(spec)
!%%!             do nz = -3,3
!%%!                nz1 = nz+dz
!%%!                do nsf1 = 1,this_nsf
!%%!                   work1(nsf1,nz1,dy,dx) = work1(nsf1,nz1,dy,dx)+FAC(nz)*work2(nsf1,dz,dy,dx)
!%%!                enddo
!%%!             end do
!%%!          end do
!%%!       end do
!%%!    end do
!%%!    ! Now work1 holds a complete convolution of one set of blip coefficients with integrals (blipcon in Cquest)
!%%!    ! So sum over blip values
!%%!    resulting_submatrix = zero
!%%!    do dx = -BlipArraySize(spec), BlipArraySize(spec)
!%%!       do dy = -BlipArraySize(spec), BlipArraySize(spec)
!%%!          do dz = -BlipArraySize(spec), BlipArraySize(spec)
!%%!             l = blip_info(spec)%blip_number(dx,dy,dz)
!%%!             if (l/=0) then
!%%!                do nsf1 = 1,this_nsf
!%%!                   do nsf2 = 1,this_nsf
!%%!                      resulting_submatrix(nsf1,nsf2) = resulting_submatrix(nsf1,nsf2) + &
!%%!                           work1(nsf1,dz,dy,dx)*blip_co(nsf2,l)
!%%!                   end do
!%%!                end do
!%%!             end if
!%%!          end do
!%%!       end do
!%%!    end do
!%%!    call scal(this_nsf*this_nsf,support_grid_spacing,resulting_submatrix,1)
!%%!    return
!%%!  end subroutine basic_get_onsite_S

end module S_matrix_module
