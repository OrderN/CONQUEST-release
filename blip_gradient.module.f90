! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: blip_gradient.module.f90,v 1.4.2.3 2006/03/31 12:09:07 drb Exp $
! ------------------------------------------------------------------------------
! Module blip_gradient
! ------------------------------------------------------------------------------
! Code area 6: Minimising energy
! ------------------------------------------------------------------------------

!!****h* Conquest/blip_gradient *
!!  NAME
!!   blip_gradient
!!  PURPOSE
!!   Contains all subroutines relating to finding the gradient of energy etc with respect to blips
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   08:25, 2003/04/03 dave
!!  MODIFICATION HISTORY
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!  SOURCE
!!
module blip_gradient

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id: blip_gradient.module.f90,v 1.4.2.3 2006/03/31 12:09:07 drb Exp $"
!!***

contains

! -----------------------------------------------------------
! Subroutine get_blip_gradient
! -----------------------------------------------------------

!!****f* blip_gradient/get_blip_gradient *
!!
!!  NAME 
!!   get_blip_gradient
!!  USAGE
!! 
!!  PURPOSE
!!   Evaluates the gradient of energy wrt blip coefficients
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   C.M.Goringe
!!  CREATION DATE
!!   19/10/95
!!  MODIFICATION HISTORY
!!   14/08/2000 TM
!!    Added new act_on_vector
!!   15/01/2001 TM
!!    Added new blip_to_grad
!!   15/05/2001 dave
!!    Added ROBODoc header, indented, changed act_on_vector
!!   16/05/2001 dave
!!    Simplified get_onsite_KE call
!!   17/05/2001 dave
!!    Shortened blip_to_grad call and inverse calls
!!   23/05/2001 dave
!!    Shortened subroutine call
!!   08:10, 31/07/2002 dave
!!    Changed call to get_support_gradient not to pass data_M12 - it's now used from matrix_data
!!   08:26, 2003/04/03 dave
!!    Included in blip_gradient
!!   08:34, 2003/04/03 dave
!!    Tidied use of GenBlas
!!   2006/09/19 08:22 dave
!!    Use temporary support fn on grid
!!  SOURCE
!!
  subroutine get_blip_gradient(inode, ionode)

    use datatypes
    use numbers
    use dimens, ONLY: grid_point_volume
    use primary_module, ONLY : bundle
    use matrix_data, ONLY : mat, Hrange
    use mult_module, ONLY: matK, scale_matrix_value, return_matrix_block_pos
    use GenBlas, ONLY: scal, copy
    use pseudopotential_data, ONLY: non_local
    use support_spec_format, ONLY: supports_on_atom, support_gradient
    use species_module, ONLY: nsf_species
    use set_bucket_module, ONLY: rem_bucket
    use calc_matrix_elements_module, ONLY: act_on_vectors_new
    use blip_grid_transform_module, ONLY: blip_to_grad_new,&
         inverse_blip_transform_new,&
         inverse_blip_to_grad_new
    use global_module, ONLY: WhichPulay, BothPulay, sf
    use functions_on_grid, ONLY: H_on_supportfns, gridfunctions, fn_on_grid, &
         allocate_temp_fn_on_grid, free_temp_fn_on_grid

    implicit none

    ! Passed variables
    integer :: inode, ionode

    !     local variables

    real(double) ::  ke, ke2, tmp
    real(double), allocatable, dimension(:,:) :: this_data_K

    integer :: i, direction, np, nn, n1, n2, this_nsf, tmp_fn


    !
    ! first, get the gradient of energy wrt the support functions.
    ! this does NOT include the change in energy due to the change in the
    ! KE matrix, but it DOES include the change in KE due to the change
    ! in S matrix.
    !
    ! a note on units... get_support_gradient and get_non_local_gradient
    ! both return values in units of 'per grid point'. The blip inverse
    ! transform assumes units of 'per au', so we rescale by multiplying
    ! by the grid point volume just before the transform. Similar arguments
    ! apply to the kinetic energy.
    !
    ! a note on signs... both the above mentioned routines give -Grad,
    ! ie the gradient downwards.
    !
    !
    ! We have support and h_on_support in support and workspace_support,
    ! so use workspace2_support as working area...
    ! This returns the support_gradient in workspace_support
    WhichPulay = BothPulay
    call get_support_gradient( inode, ionode)

    ! now we need to accumulate the 'type 1' gradient of the non-local energy
    ! (see notes, 3/1/97)
    if (non_local) call get_non_local_gradient( H_on_supportfns, inode, ionode)

    ! now we need to transform this into the blip basis
    ! first, apply the scaling for grid size,

    call scal( gridfunctions(H_on_supportfns)%size, grid_point_volume, gridfunctions(H_on_supportfns)%griddata, 1 )

    ! and then 
    call inverse_blip_transform_new(inode-1, H_on_supportfns, support_gradient, bundle%n_prim)

    ! Kinetic Energy; first, the change in onsite T matrix elements,
    ! and as we do so, clear the diagonal blocks of data K
    i = 1
    this_nsf = bundle%species(1)
    allocate(this_data_K(this_nsf,this_nsf))
    do np = 1, bundle%groups_on_node
       do nn = 1,bundle%nm_nodgroup(np)
          if(nsf_species(bundle%species(i))/=this_nsf) then
             deallocate(this_data_K)
             this_nsf = nsf_species(bundle%species(i))
             allocate(this_data_K(this_nsf,this_nsf))
          end if
          this_data_K = zero
          call return_matrix_block_pos(matK,mat(np,Hrange)%onsite(nn),this_data_K,this_nsf*this_nsf)
          call get_onsite_KE_gradient( supports_on_atom(i),this_data_K, support_gradient(i), this_nsf, bundle%species(i))
          do n1 = 1,this_nsf
             do n2 = 1,this_nsf
                call scale_matrix_value(matK,np,nn,i,0,n1,n2,zero,1)
             end do
          end do
          i = i+1
       enddo
    end do
    deallocate(this_data_K)
    ! Now, for the offsite part, done one the integration grid.
    ! to do the KE bit we need to blip_to_grad transform, act with K,
    ! and inverse_blip_grad_transform back (this routine ACCUMULATES onto
    ! data_dblip)
    !
    ! because we now do the onsite KE analytically, we need to remove
    ! the onsite parts of this term. This can be done by removing the
    ! onsite parts of the K matrix (because KE = Tr[KT], this is
    ! equivalent to removing the onsite parts of T, and somewhat easier!
    !
    !
    ! first, for x: do the blip_grad transform into workspace_support 
    tmp_fn = allocate_temp_fn_on_grid(sf)
    do direction = 1, 3 
       call blip_to_grad_new(inode-1, direction, H_on_supportfns)
       ! now act with K - result into workspace2_support
       gridfunctions(tmp_fn)%griddata = zero
       call act_on_vectors_new(inode-1,rem_bucket(3),matK, tmp_fn,H_on_supportfns)

       ! apply the appropriate scaling factor
       call scal( gridfunctions(tmp_fn)%size, -two*grid_point_volume, &
            gridfunctions(tmp_fn)%griddata, 1 )
       call inverse_blip_to_grad_new(inode-1,direction, tmp_fn, support_gradient,bundle%n_prim)
    end do
    call free_temp_fn_on_grid(tmp_fn)
    return
  end subroutine get_blip_gradient
!!***

! -----------------------------------------------------------
! Subroutine get_electron_gradient
! -----------------------------------------------------------

!!****f* blip_gradient/get_electron_gradient *
!!
!!  NAME 
!!   get_electron_gradient
!!  USAGE
!! 
!!  PURPOSE
!!   This subroutine obtains the gradient of the electron
!!   number with respect to the support functions on
!!   the grid. The gradient with respect to
!! |phi_i> at grid point n is given by
!!   24 sum_j (LSLij  - LSLSLij) |phi_j(n)>
!! 
!!   M4 = 24 (LSL - LSLSL) 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   E.H.Hernandez
!!  CREATION DATE
!!   01/03/95
!!  MODIFICATION HISTORY
!!   12/9/95 EH to conform to new data storage
!!   8/5/96 CG this routine becomes trivial...
!!   7/7/98 DRB to use act_on_vector
!!   19/9/2000 TM
!!    Added new act_on_vector calls
!!   15/05/2001 dave
!!    Converted to F90, added ROBODoc and updated act_on_vector call
!!   23/05/2001 dave
!!    Shortened subroutine call
!!   24/05/2001 dave
!!    Bug fix - added use maxima_module
!!   08:27, 2003/04/03 dave
!!    Included in blip_gradient
!!   08:34, 2003/04/03 dave
!!    Tidied use of GenBlas
!!  SOURCE
!!
  subroutine get_electron_gradient(support, electron_gradient, inode, ionode)

    use datatypes
    use numbers
    use GenBlas, ONLY: scal
    use dimens, ONLY: grid_point_volume
    use set_bucket_module, ONLY: rem_bucket
    use calc_matrix_elements_module, ONLY: act_on_vectors_new
    use blip_grid_transform_module, ONLY: inverse_blip_transform_new
    use mult_module, ONLY: matM4
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid
    use primary_module, ONLY : bundle
    use support_spec_format, ONLY: support_elec_gradient

    implicit none

    ! Passed variables
    integer :: inode, ionode

    integer :: support, electron_gradient

    gridfunctions(electron_gradient)%griddata = zero
    call act_on_vectors_new(inode-1,rem_bucket(1),matM4,electron_gradient,support)
    call scal( gridfunctions(electron_gradient)%size, -one * grid_point_volume, &
         gridfunctions(electron_gradient)%griddata, 1 )
    call inverse_blip_transform_new(inode-1,electron_gradient, support_elec_gradient, bundle%n_prim)
    return
  end subroutine get_electron_gradient
!!***

! -----------------------------------------------------------
! Subroutine get_support_gradient
! -----------------------------------------------------------

!!****f* blip_gradient/get_support_gradient *
!!
!!  NAME 
!!   get_support_gradient
!!  USAGE
!! 
!!  PURPOSE
!!   This subroutine obtains the gradient of the total
!!   energy with respect to the support functions on
!!   the grid. The gradient of E with respect to
!! |phi_i> at grid point n is given by
!!
!!    dE/d|phi_i(n)> = 4 sum_j 
!!            (      Kij H +          <- Type I
!!                   M1ij - M2ij    )   <- Type II
!!   |phi_j(n)>
!! 
!!    where H is the Hamiltonian operator (not its 
!!    matrix representation), which must be made to 
!!    act on |phi_j>, and the matrices M1 and M2 are
!!    defined as
!!    M1 = 3 LHL 
!!    M2 = 2 LSLHL + LHLSL
!!    respectively.
!!
!!    Notes 4/6/94, section 11.2 (pp24-34)
!!    The Hamiltonian matrix here is the complete matrix, but the Hamiltonian
!!    operator is (H_local - T); the non-local type I derivative and the
!!    kinetic type I derivative have to be done seperately. 
!!    See Notes 24/2/96 pp9-10 and Conquest paper IV for KE discussion,
!!        Notes 31/1/97 pp1-5 for Non-Local energy
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   E.H.Hernandez
!!  CREATION DATE
!!   01/03/95
!!  MODIFICATION HISTORY
!!   12/09/95 EH to conform to new data storage
!!   02/05/96 EH to speed up and take advantage of
!!             block data distribution.
!!   03/05/96 CG to simplify above changes
!!   19/9/2000 TM
!!    Added new act_on_vector calls
!!   15/05/2001 dave
!!    F90, ROBODoc and act_on_vector sorted
!!   21/05/2001 dave
!!    Shortened variables passed
!!   31/07/2002 dave
!!    Changed to use data_M12 from matrix_data instead of passing it
!!   08:27, 2003/04/03 dave
!!    Included in blip_gradient
!!   08:34, 2003/04/03 dave
!!    Tidied use of GenBlas
!!  SOURCE
!!
  subroutine get_support_gradient( inode, ionode)

    use datatypes
    use numbers
    use mult_module, ONLY: matM12, matK
    use GenBlas, ONLY: axpy, copy, scal
    use set_bucket_module,     ONLY: rem_bucket
    use calc_matrix_elements_module, ONLY: act_on_vectors_new
    use global_module, ONLY: flag_basis_set, PAOs, blips, WhichPulay, PhiPulay, BothPulay, SPulay, sf
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid, supportfns, H_on_supportfns, &
         allocate_temp_fn_on_grid, free_temp_fn_on_grid

    implicit none

    !     Shared variables
    integer :: inode, ionode

    ! Local variables
    integer :: tmp_fn

    tmp_fn = allocate_temp_fn_on_grid(sf)
    ! first, act on workspace_support (which holds h_on_support on entry)
    ! with the K matrix) to get type I variation
    gridfunctions(tmp_fn)%griddata = zero
    if(WhichPulay==PhiPulay.OR.WhichPulay==BothPulay) then
       call act_on_vectors_new(inode-1,rem_bucket(3),matK,tmp_fn, H_on_supportfns)

    ! store the result in workspace_support (which will hold support_gradient
    ! on exit)

       call copy( gridfunctions(H_on_supportfns)%size, gridfunctions(tmp_fn)%griddata,1, &
            gridfunctions(H_on_supportfns)%griddata,1)
    else
       gridfunctions(H_on_supportfns)%griddata = zero
    end if
    ! now act with M12 on support to get type II variation
    gridfunctions(tmp_fn)%griddata = zero
    if(flag_basis_set==blips.AND.(WhichPulay==SPulay.OR.WhichPulay==BothPulay)) then
       call act_on_vectors_new (inode-1,rem_bucket(1),matM12,tmp_fn,supportfns)
       call axpy( gridfunctions(H_on_supportfns)%size, one,gridfunctions(tmp_fn)%griddata,1, &
            gridfunctions(H_on_supportfns)%griddata, 1)
    end if
    call scal( gridfunctions(H_on_supportfns)%size, minus_four, gridfunctions(H_on_supportfns)%griddata, 1 )
    call free_temp_fn_on_grid(tmp_fn)
    return
  end subroutine get_support_gradient
!!***

! -----------------------------------------------------------
! Subroutine get_non_local_gradient
! -----------------------------------------------------------

!!****f* blip_gradient/get_non_local_gradient *
!!
!!  NAME 
!!   get_non_local_gradient
!!  USAGE
!! 
!!  PURPOSE
!!   accumulates the gradient of energy due to type 1
!!   variations of the non-local energy onto the grid
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   C.M.Goringe
!!  CREATION DATE
!!   24/01/97
!!  MODIFICATION HISTORY
!!   07/07/98 by Dave Bowler to use act_on_vector
!!   15/05/2001 dave
!!    Changed act_on_vector, added ROBODoc header, indented
!!   23/05/2001 dave
!!    Shortened subroutine call
!!   15:42, 04/02/2003 drb 
!!    Changed to use temp_TCS (problem with transposes if NSF/=NCF)
!!   08:15, 2003/04/04 dave
!!    Included into blip_gradient
!!  SOURCE
!!
  subroutine get_non_local_gradient(support_gradient, inode, ionode)

    use datatypes
    Use mult_module, ONLY: H_SP_SP, SP_trans, mult, matrix_product, matrix_scale, matrix_transpose, &
         matSC, matCS, matU, matK
    use numbers, ONLY: zero, minus_four
    use set_bucket_module, ONLY: rem_bucket
    use calc_matrix_elements_module, ONLY: act_on_vectors_new
    use functions_on_grid, ONLY: pseudofns

    implicit none

    ! Passed variables
    integer :: inode, ionode
    integer :: support_gradient

    ! First of all, find U (=K.SC)
    call matrix_transpose(matSC,matCS)
    call matrix_product(matK,matCS,matU,mult(H_SP_SP))
    call matrix_scale(minus_four,matU)

    ! Finally, act on the non-local projectors with the scaled U
    call act_on_vectors_new(inode-1,rem_bucket(2),matU,support_gradient,pseudofns)

    return
  end subroutine get_non_local_gradient
!!***

!!****f* blip_gradient/get_onsite_KE_gradient *
!!
!!  NAME 
!!   get_onsite_KE_gradient
!!  USAGE
!! 
!!  PURPOSE
!!   This routine accumulates onto the blip gradient the derivative in
!!   KE due to the change in T wrt the blip coefficients
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   C.M.Goringe
!!  CREATION DATE
!!   18/11/96 
!!  MODIFICATION HISTORY
!!   16/05/2001 dave
!!    F90, ROBODoc and:
!!    REMOVED DEPENDENCE ON NSF=4
!!   16/05/2001 dave
!!    Simplified call
!!   08:17, 2003/04/04 dave
!!    Included into blip_gradient
!!  SOURCE
!!
  subroutine get_onsite_KE_gradient( this_data_blip,this_data_K,this_blip_grad, this_nsf, spec)

    use datatypes
    use numbers
    use GenBlas, ONLY: axpy, copy, scal
    use blip, ONLY: blip_info, BlipArraySize, OneArraySize, FullArraySize, SupportGridSpacing
    use support_spec_format, ONLY: support_function
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_minE
    use memory_module, ONLY: reg_alloc_mem, type_dbl, reg_dealloc_mem

    implicit none

    ! Passed Variables
    integer :: this_nsf, spec

    real(double) :: this_data_K(this_nsf*this_nsf) ! N.B. Here NSF is for ONE site only so is OK
    type(support_function) :: this_data_blip, this_blip_grad

    ! Local Variables

    integer, parameter :: MAX_D=3

    real(double) ::  FAC(0:MAX_D), D2FAC(0:MAX_D)

    real(double), allocatable, dimension(:) :: work1, work2, work3, work4, work5, work6

    integer :: dx, dy, dz, offset, l, at, nsf1, nsf2, stat

    allocate(work1(FullArraySize(spec)*this_nsf),work2(FullArraySize(spec)*this_nsf),work3(FullArraySize(spec)*this_nsf), &
         work4(FullArraySize(spec)*this_nsf),work5(FullArraySize(spec)*this_nsf),work6(FullArraySize(spec)*this_nsf), &
         STAT=stat)
    if(stat/=0) call cq_abort("Error allocating arrays for onsite KE blip grad: ",FullArraySize(spec),this_nsf)
    call reg_alloc_mem(area_minE,6*FullArraySize(spec)*this_nsf,type_dbl)
    ! first, we copy the blip functions for this atom onto a cubic grid;
    ! we make this grid 'too big' in order to have a fast routine below.
    !
    ! To get gradients, we start by acting with the K matrix upon the
    ! coefficients

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
                   work1(nsf1+at) = zero
                   do nsf2 = 1,this_nsf
                      work1(nsf1+at) = work1(nsf1+at) + &
                           this_data_K((nsf2-1)*this_nsf+nsf1)*this_data_blip%supp_func(nsf2)%coefficients(l)
                   enddo
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

    call copy(FullArraySize(spec)*this_nsf,work1,1,work3,1)
    call scal(FullArraySize(spec)*this_nsf,D2FAC(0),work3,1)
    do dz = 1, MAX_D
       offset = dz * OneArraySize(spec) * OneArraySize(spec) * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), D2FAC(dz), &
            work1(1:), 1, work3(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), D2FAC(dz), &
            work1(1+offset:), 1, work3(1:), 1 )
    end do

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

    call copy(FullArraySize(spec)*this_nsf,work2,1,work5,1)
    call scal(FullArraySize(spec)*this_nsf,D2FAC(0),work5,1)
    do dy = 1, MAX_D
       offset = dy * OneArraySize(spec) * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), D2FAC(dy), &
            work2(1:), 1, work5(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), D2FAC(dy), &
            work2(1+offset:), 1, work5(1:), 1 )
    end do

    call axpy(FullArraySize(spec)*this_nsf,FAC(0),work3,1,work5,1)
    do dy = 1, MAX_D
       offset = dy * OneArraySize(spec) * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dy), &
            work3(1:), 1, work5(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dy), &
            work3(1+offset:), 1, work5(1:), 1 )
    end do

    ! and x - put it all into 6

    call copy(FullArraySize(spec)*this_nsf,work5,1,work6,1)
    call scal(FullArraySize(spec)*this_nsf,FAC(0),work6,1)
    do dx = 1, MAX_D
       offset = dx * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dx), &
            work5(1:), 1, work6(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dx), &
            work5(1+offset:), 1, work6(1:), 1 )
    end do

    call axpy(FullArraySize(spec)*this_nsf,D2FAC(0),work4,1,work6,1)
    do dx = 1, MAX_D
       offset = dx * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), D2FAC(dx), &
            work4(1:), 1, work6(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), D2FAC(dx), &
            work4(1+offset:), 1, work6(1:), 1 )
    end do

    ! now accumulate work6 onto the gradient

    call scal(FullArraySize(spec)*this_nsf, two*SupportGridSpacing(spec),work6,1)

    offset = BlipArraySize(spec) + 1
    do dx = -BlipArraySize(spec), BlipArraySize(spec)
       do dy = -BlipArraySize(spec), BlipArraySize(spec)
          do dz = -BlipArraySize(spec), BlipArraySize(spec)
             l = blip_info(spec)%blip_number(dx,dy,dz)
             if (l.ne.0) then
                at = (((dz+offset)*OneArraySize(spec) + (dy+offset))*OneArraySize(spec) + &
                     (dx+offset)) * this_nsf
                do nsf1 = 1,this_nsf
                   this_blip_grad%supp_func(nsf1)%coefficients(l) = &
                        this_blip_grad%supp_func(nsf1)%coefficients(l) - work6(nsf1+at)
                enddo
             end if
          end do
       end do
    end do

    deallocate(work1,work2,work3, work4,work5,work6, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating arrays for onsite KE blip grad: ",FullArraySize(spec),this_nsf)
    call reg_dealloc_mem(area_minE,6*FullArraySize(spec)*this_nsf,type_dbl)
    return
  end subroutine get_onsite_KE_gradient
!!***

end module blip_gradient
