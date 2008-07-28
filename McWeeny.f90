! $Id$
! -----------------------------------------------------------
! Module McWeeny
! -----------------------------------------------------------
! Code area 4: density matrix
! -----------------------------------------------------------

!!****h* Conquest/McWeeny
!!  NAME
!!   McWeeny
!!  PURPOSE
!!   To handle the McWeeny/Manolopoulos initialisation of density matrix
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   21/10/98
!!  MODIFICATION HISTORY
!!   11/05/01 DRB
!!    ROBODoc headers and better output
!!   08/06/2001 dave
!!    Added RCS Id and Log tags; used GenComms throughout
!!    Changed my_barrier to come from GenComms
!!   11:09, 04/02/2003 drb 
!!    Major change: stability fix in InitMcW to ensure mu lies between Hamiltonian limits
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2006/11/21 16:53 dave
!!    Changed iprint statements
!!   2008/02/01 17:48 dave
!!    Changes for output to file not stdout
!!***
module McWeeny

  use global_module, ONLY: io_lun

  implicit none

contains

!!****f* McWeeny/McWMin *
!!
!!  NAME 
!!   McWMin
!!  USAGE
!! 
!!  PURPOSE
!!   Performs iterative minimisation using McWeeny algorithm
!!  INPUTS
!! 
!! 
!!  USES
!!   datatypes, numbers, logicals, maxima_module, matrix_data
!!   mult_module, multiply_module, GenBlas
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   21/10/98
!!  MODIFICATION HISTORY
!!   11/05/01 DRB
!!    ROBODoc headers and better output also nsf from maxima_module
!!   08/06/2001 dave
!!    Tidied and indented
!!   10:38, 16/05/2007 drb 
!!    Moved matrix_sum to stop error at end
!!  SOURCE
!!
  subroutine McWMin( n_L_iterations, deltaE)

    use datatypes
    use numbers
    use matrix_data, ONLY: Lrange
    use mult_module, ONLY: allocate_temp_matrix, matrix_sum, free_temp_matrix, matL, matH, matS, matT, symmetrise_L
    use global_module, ONLY: iprint_DM
    use GenComms, ONLY: inode, ionode

    implicit none

    ! Passed variables
    integer :: n_L_iterations
    real(double) :: deltaE

    ! Local variables
    integer :: n_iterations
    integer :: matRhoNew, mat_oldL

    real(double) :: cn, oldE, omega1, c_old

    oldE = 1.0e30_double
    mat_oldL = allocate_temp_matrix(Lrange,0)
    matRhoNew = allocate_temp_matrix(Lrange,0)
    call matrix_sum(zero,mat_oldL,one,matL)
    do n_iterations=1,n_L_iterations
       if (inode==ionode.and.iprint_DM>=1) write(io_lun,2) n_iterations

       ! find new rho and energy
       call McW_matrix_multiply( matH, matS, matL, matT, matRhoNew, omega1, cn  )
       ! If new energy is lower, copy new rho into old rho and repeat
       if(omega1.GT.oldE.OR.cn.GT.one.OR.cn.LT.zero) then
          if(inode==ionode.AND.iprint_DM>0) write(io_lun,4) c_old, oldE, cn,omega1
4         format(2x,'Rounding error found: ',4f15.6)
          exit
       endif
       deltaE = omega1 - oldE
       call matrix_sum(zero,matL,one,matRhoNew)
       oldE = omega1
       c_old = cn
       if(inode==ionode.AND.iprint_DM>=0) write(io_lun,5) oldE,c_old
5      format(3x,'Energy: ',f15.6,' Mid-point: ',f15.6)
       call symmetrise_L()
    end do
    if(n_iterations<3) call matrix_sum(zero,matL,one,mat_oldL)
    call free_temp_matrix(matRhoNew)
    call free_temp_matrix(mat_oldL)
    if (inode.eq.ionode.and.iprint_DM.ge.1) write(io_lun,3) n_iterations, omega1
2   format(/,20x,'McWeeny L iteration:',i5)
3   format(/,20x,'Functional value reached after ',i5,' L iterations: ', &
         /,20x,' Omega: ', f15.7)
    return
  end subroutine McWMin
!!***

!!****f* McWeeny/InitMcW *
!!
!!  NAME 
!!   InitMcW 
!!  USAGE
!! 
!!  PURPOSE
!!   Initialises L to be a function of H with eigenvalues between 0 and 1
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   End 10/98
!!  MODIFICATION HISTORY
!!   11/05/01 DRB
!!    ROBODoc headers and better output also nsf from maxima_module
!!   08/06/2001 dave
!!    Tidying and added GenComms for gsum
!!    Added my_barrier to use GenComms statement
!!   15:25, 10/12/2002 drb 
!!    Bug fix - changed location of output for mubar and added gmin and gmax calls on hmin and hmax
!!   11:08, 04/02/2003 drb 
!!    Implemented quick fix for bracketing mu with Hamiltonian limits (required for stability)
!!   12:30, 2004/06/09 dave
!!    Fixed assumption that T and S have same range
!!   2004/10/29 drb
!!    Removed gsum on SX: MAJOR bug found by TM which was messing up initialisation; 
!!    also removed loops to change hmin and hmax (unnecessary, I think)
!!  SOURCE
!!
  subroutine InitMcW(number_of_bands,inode,ionode)

    use datatypes
    use numbers
    use global_module, ONLY: iprint_DM, ni_in_cell
    use species_module, ONLY: species, nsf_species
    use mult_module, ONLY: allocate_temp_matrix, matH, matS, matT, matL, mult, &
         matrix_trace, matrix_product_trace, matrix_product, T_S_TS, free_temp_matrix
    use matrix_data, ONLY: Lrange, Srange, TSrange
    use GenBlas
    use GenComms, ONLY: gsum, my_barrier, gmin, gmax

    implicit none

    ! Passed variables
    real(double) ::  number_of_bands
    integer :: inode, ionode

    ! Local variables
    real(double) :: n_e, n_o, hmax, hmin, SX, SXHX
    real(double) :: A, mu1, mu2, mubar, lambda
    integer :: matXHX, mat_temp, matTS
    integer :: length, i

    if(inode==ionode.AND.iprint_DM>1) write(io_lun,1)
1   format(1x,'Welcome to InitMcW')
    ! We must first initialise rho
    n_e = number_of_bands
    !n_o = 2.0_double * number_of_bands
    n_o = zero
    do i=1,ni_in_cell
       n_o = n_o + nsf_species(species(i))
    end do
    if(inode==ionode.AND.iprint_DM>1) write(io_lun,2) n_e,n_o
2   format(2x,'Electrons: ',f15.6,' Orbitals: ',f15.6)
    matXHX = allocate_temp_matrix(Lrange,0)
    mat_temp = allocate_temp_matrix(Srange,0)
    matTS = allocate_temp_matrix(TSrange,0)
    call McWXHX(matH, matT, matXHX, mat_temp)
    call matrix_product(matT,matS,matTS, mult( T_S_TS ))
    SX = matrix_trace(matTS)
    write(io_lun,*) 'Trace of temp: ',matrix_trace(mat_temp)
    SXHX = matrix_product_trace(matS,mat_temp)
    call my_barrier()
    if(inode==ionode.AND.iprint_DM>1) write(io_lun,*) 'SX, SXHX are ',SX,SXHX
    if(inode==ionode.AND.iprint_DM>1) write(io_lun,*) 'SX/n_o is ',SX/n_o
    ! Get Gershgorin limits on H
    call GetHLimits(hmin,hmax)
    call gmin(hmin)
    call gmax(hmax)
    if(inode==ionode.AND.iprint_DM>1) write(io_lun,3) hmin,hmax
3   format(2x,'Minimum and maximum limits on H are ',2f15.6)
    A = (1.0_double - SX/n_o) 
    mu1 = SXHX/n_o + hmax*A
    ! DRB 2004/09/16 Fixing typo ?
    mu2 = (hmin*A*n_e/(n_o-n_e) - SXHX/n_o)/(n_o*A/(n_o-n_e) - one)
    if(inode==ionode) write(io_lun,*) 'Mu1,2: ',mu1,mu2
    if((mu1<hmax.AND.mu1>hmin).AND.(abs(n_e/(hmax-mu1)).LT.abs((n_o - n_e)/(mu2-hmin)))) then
       mubar = mu1
       lambda = n_e/(hmax-mu1)
    else
       mubar = mu2
       lambda = (n_o - n_e)/(mu2-hmin)
    endif
    if(inode==ionode.AND.iprint_DM>1) write(io_lun,4) mubar
4   format(2x,'Mubar is ',f15.6)
    if(inode==ionode.AND.iprint_DM>1) write(io_lun,*) 'lambda and mu are ',lambda,mubar
    call McWRho0(matL, matT, matXHX, mubar, lambda, n_e, n_o)
    call free_temp_matrix(matTS)
    call free_temp_matrix(mat_temp)
    call free_temp_matrix(matXHX)
    return
  end subroutine InitMcW
!!***

!!****f* McWeeny/McW_matrix_multiply *
!!
!!  NAME 
!!   McW_matrix_multiply
!!  USAGE
!! 
!!  PURPOSE
!!   Performs one iteration of McWeeny iterative scheme
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   21/10/98
!!  MODIFICATION HISTORY
!!   11/05/01 DRB
!!    ROBODoc headers and better output also nsf from maxima_module
!!   08/06/2001 dave
!!    Added GenComms for gsum
!!  SOURCE
!!
  subroutine McW_matrix_multiply( matH, matS, matL, matInvS,matRhoNew, energy, c  )

    use datatypes
    use numbers
    use primary_module, ONLY : bundle
    use matrix_data, ONLY: LSLrange, Lrange, Hrange, Srange
    use mult_module, ONLY: matrix_product, matrix_sum, matrix_product_trace, matrix_scale, matrix_transpose, &
         allocate_temp_matrix, free_temp_matrix, matLS, matSL, L_S_LS, LS_L_LSL, LSL_SL_L, mult
    use global_module, ONLY: iprint_DM
    use GenComms, ONLY: gsum, cq_abort, inode, ionode, my_barrier

    implicit none

    ! Passed variables
    integer :: matH, matS, matL, matInvS, matRhoNew
    real( double ) energy, c

    ! Local Variables
    real( double ) c1,c2,c3,cn
    integer :: matLSL, matLSLSL, mat_top, mat_bottom, matL2

    call matrix_product(matL, matS, matLS, mult(L_S_LS))
    ! Generate LSL, LSLSL
    matLSL = allocate_temp_matrix(LSLrange,0)
    matLSLSL = allocate_temp_matrix(Lrange,0)
    mat_top = allocate_temp_matrix(Srange,0)
    mat_bottom = allocate_temp_matrix(Srange,0)
    call my_barrier
    call matrix_product(matLS,matL,matLSL,mult( LS_L_LSL))
    call matrix_product(matLSL,matLS,matLSLSL,mult( LSL_SL_L))
    !****************************************************************
    ! Shorten L, LSL and LSLSL to range S and take traces (can do it
    ! with sdot because S is symmetric)
    !****************************************************************
    call matrix_sum(zero,mat_top,one,matLSL)
    call matrix_sum(one,mat_top,-one,matLSLSL)
    c1 = matrix_product_trace(matS,mat_top)
    if(inode==ionode.AND.iprint_DM>=2) write(io_lun,*) 'S.top is ',c1
    call matrix_sum(zero,mat_bottom,one,matL)
    c1 = matrix_product_trace(matS,mat_bottom)
    if(inode==ionode.AND.iprint_DM>=2) write(io_lun,*) 'N_e is ',c1
    call matrix_sum(one,mat_bottom,-one,matLSL)
    c1 = matrix_product_trace(matS,mat_top)
    c2 = matrix_product_trace(matS,mat_bottom)
    if(c2/=zero) then
       cn = c1/c2
    else
       call cq_abort('McW_matrix_multiply: c2 is zero')
    endif
    c = cn
    if(inode.eq.ionode.AND.iprint_DM>=2) write(io_lun,*) 'c, c1,c2 are ',cn,c1,c2
    ! Shorten LSL and LSLSL to range L
    call matrix_scale(zero,matRhoNew)
    call matrix_sum(zero,matRhoNew,-one,matLSLSL)
    call matrix_sum(one,matRhoNew,one+cn,matLSL)
    if(cn<0.5_double) then
       call matrix_sum(one,matRhoNew,(one - two*cn),matL)
       call matrix_scale((one/(one - cn)),matRhoNew)
    else
       call matrix_scale(one/cn,matRhoNew)
    endif
    ! Idempotency checks
    !%%! data_L3 = zero
    !%%! call axpy(NSF*NSF*mat(1,Lrange)%length,one,dat_L,1,data_L3,1)
    !%%! ! LSLSL -> 2LSLSL - 3LSL
    !%%! call matrix_add(2.0_double, data_LSLSL, mat(:,Lrange), &
    !%%!      -3.0_double, data_LSL, mat(:,LSLrange), bundle, inode)
    !%%! ! Lnew -> Lnew+(2LSLSL-3LSL) [= Lnew - McW(Lold)]
    !%%! call axpy(NSF*NSF*mat(1,Lrange)%length,one,data_LSLSL,1,data_L3,1)
    !%%! call mult_wrap(inode-1,data_L3,dat_S,data_LS,mult(L_S_LS))
    !%%! call mult_wrap(inode-1,dat_S,data_SL,data_XHX,mult(S_LS_L))
    !%%! cn = dot(NSF*NSF*mat(1,Lrange)%length,data_SLS2, 1, data_L3, 1)
    !%%! call gsum(cn)
    !%%! if(inode.eq.ionode.AND.iprint_DM>=2) write(io_lun,*) 'Idempot monitor is ',cn/256.0_double 
    ! Shorten rhonew to range H into temp
    call free_temp_matrix(mat_bottom)
    mat_bottom = allocate_temp_matrix(Hrange,0)
    call matrix_sum(zero,mat_bottom,one,matRhoNew)
    energy = two*matrix_product_trace(mat_bottom,matH)
    if(inode==ionode.AND.iprint_DM>=3) write(io_lun,*) 'energy is ',energy
    call matrix_sum(zero,mat_top,one,matRhoNew)
    c1 = matrix_product_trace(mat_top,matS)
    if(inode==ionode.AND.iprint_DM>=3) write(io_lun,*) 'N_e(2) is ',c1
    call free_temp_matrix(mat_bottom)
    call free_temp_matrix(mat_top)
    call free_temp_matrix(matLSLSL)
    call free_temp_matrix(matLSL)
  End Subroutine McW_matrix_multiply
!!***

! sbrt McWXHX: Pre- and post- multiplies H by S^-1, and returns the result 
! to both L range (bab) and S range (tm)
  Subroutine McWXHX(matA, matB, matBAB, mat_temp)

    use datatypes
    use numbers
    use matrix_data, ONLY: THrange
    use mult_module, ONLY: matrix_product, matrix_sum, T_H_TH, TH_T_L, &
         allocate_temp_matrix, mult, free_temp_matrix, matrix_trace
    use GenComms, ONLY : inode, ionode
    use global_module, ONLY: iprint_DM

    implicit none

    ! Passed variables
    integer :: matA, matB, matBAB, mat_temp

    ! Local variable
    integer :: matBA
    real(double) :: xx

    matBA = allocate_temp_matrix(THrange,0)
    call matrix_product(matB, matA, matBA, mult( T_H_TH ))
    if(iprint_DM >= 2) xx=matrix_trace(matBA)
    if(inode.eq.ionode.AND.iprint_DM>=2) write(io_lun,*) 'Trace of BA: ',xx
    call matrix_product(matBA, matB, matBAB, mult( TH_T_L ))
    if(iprint_DM >= 2) xx=matrix_trace(matBAB)
    if(inode.eq.ionode.AND.iprint_DM>=2) write(io_lun,*) 'Trace of BAB: ',xx
    call matrix_sum(zero, mat_temp, one, matBAB)
    if(iprint_DM >= 2) xx=matrix_trace(mat_temp)
    if(inode.eq.ionode.AND.iprint_DM>=2) write(io_lun,*) 'Trace of BAB: ',xx
    call free_temp_matrix(matBA)

  End Subroutine McWXHX

  Subroutine McWRho0(matA, matB, matC, m, l, n_e, n_o)

    use datatypes
    use numbers
    use GenComms, ONLY : inode, ionode
    use mult_module, ONLY : matrix_sum
    use global_module, ONLY: iprint_DM

    implicit none

    integer :: matA, matB, matC
    Real( double ) :: m, l, n_e, n_o, tmp

    tmp = (l*m+n_e)/n_o
    if(inode.eq.ionode.AND.iprint_DM>=2) write(io_lun,*) 'l,m,n_e,n_o and tmp are ',l,m,n_e,n_o,tmp
    call matrix_sum(zero,matA, tmp, matB)
    tmp = -l/n_o
    if(inode.eq.ionode.AND.iprint_DM>=2) write(io_lun,*) 'tmp is ',tmp
    call matrix_sum(one, matA, tmp, matC)

  End Subroutine McWRho0

! -----------------------------------------------------------------------
! Finds the Gershgorin limits on the eigenvalues of the Hamiltonian
! -----------------------------------------------------------------------
  subroutine GetHLimits(MinH,MaxH)

    use datatypes
    use matrix_module, ONLY: matrix, matrix_halo
    use matrix_data, ONLY: mat, Hrange, halo
    use cover_module, ONLY: BCS_parts
    use primary_module, ONLY: bundle
    use mult_module, ONLY: matH, return_matrix_value, matrix_pos
    use numbers
    use GenComms,ONLY: cq_abort
    use maxima_module, ONLY: maxnsf

    implicit none

    ! Passed variables
    real(double) :: MinH, MaxH
    
    ! Local variables
    integer :: np, nat, nb, ist, i, j, gcspart, wheremat, iprim, nsfi, nsfj, stat
    real(double), allocatable, dimension(:) :: HMin, HMax
    real(double) :: Hval

    MinH = 1.0e30_double
    MaxH = -1.0e30_double
    iprim = 0
    allocate(Hmin(maxnsf),Hmax(maxnsf),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating in GetHLimits: ",stat,nsfi)
    do np = 1, bundle%groups_on_node
       if(bundle%nm_nodgroup(np)>0) then
          do nat = 1,bundle%nm_nodgroup(np)
             iprim=iprim+1
             nsfi = mat(np,Hrange)%ndimi(nat)             
             Hmin = zero
             Hmax = zero
             do nb = 1, mat(np,Hrange)%n_nab(nat)
                ist = mat(np,Hrange)%i_acc(nat)+nb-1
                nsfj = mat(np,Hrange)%ndimj(ist)
                gcspart = BCS_parts%icover_ibeg(mat(np,Hrange)%i_part(ist))+mat(np,Hrange)%i_seq(ist)-1
                !wheremat = matrix_pos(matH,iprim,halo(Hrange)%i_halo(gcspart),1,1)
                wheremat = matrix_pos(matH,iprim,halo(Hrange)%i_halo(gcspart))
                if(wheremat==mat(np,Hrange)%onsite(nat)) then
                   do i = 1,nsfj
                      do j=1,nsfi
                         Hval = return_matrix_value(matH,np,nat,iprim,nb,j,i)
                         !if(i==j) then
                         HMin(i) = HMin(i) + abs(Hval)
                         !else
                         !   HMin(i) = HMin(i) - abs(return_matrix_value(matH,np,nat,iprim,nb,j,i))
                         !end if
                         HMax(i) = HMax(i) + abs(Hval)
                      enddo
                   enddo
                else
                   do i = 1,nsfj
                      do j=1,nsfi
                         Hval = return_matrix_value(matH,np,nat,iprim,nb,j,i)
                         HMin(i) = HMin(i) - abs(Hval)
                         HMax(i) = HMax(i) + abs(Hval)
                      enddo
                   enddo
                endif
             enddo
             do i=1,nsfi
                if(HMin(i)<MinH) MinH = HMin(i)
                if(HMax(i)>MaxH) MaxH = HMax(i)
             enddo
          end do ! nat = bundle%nm_nodgroup
       endif
    enddo
    deallocate(Hmin,Hmax)
    return
  end subroutine GetHLimits

end module McWeeny
