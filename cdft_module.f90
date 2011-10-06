! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cdft_module
! ------------------------------------------------------------------------------
! Code area 6: energy minimisation
! ------------------------------------------------------------------------------

!!***h* Conquest/cdft_module *
!!  NAME
!!   cdft_module
!!  CREATION DATE
!!   2009/02/01 A. M. P. Sena and D. R. Bowler
!!  MODIFICATION HISTORY
!!   2011/07/21 16:41 dave
!!    Preparing for inclusion in main trunk
!!  SOURCE
!!
module cdft_module

  use datatypes
  use cdft_data
  use global_module, ONLY: io_lun, iprint_SC
  use GenComms, ONLY: inode, ionode

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"

!!***

contains

!!****f* cdt_module/init_cdft *
!!
!!  NAME 
!!   init_cdft
!!  USAGE
!!   
!!  PURPOSE
!!   Initialises cDFT
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   DRB
!!  CREATION DATE
!!   2011/08
!!  MODIFICATION HISTORY
!!  
!!  SOURCE
!!  
  subroutine init_cdft
    
    use datatypes
    use numbers
    use global_module, ONLY: ni_in_cell, flag_cdft_atom
    use mult_module, ONLY: allocate_temp_matrix
    use matrix_data, ONLY: Hrange
    use density_module, ONLY: bwgrid
    use maxima_module, ONLY: maxngrid
    use GenComms, ONLY: cq_abort

    implicit none

    integer :: i, j, stat

    if(cDFT_NumberAtomGroups>2) then
       call cq_abort("Maximum of two groups permitted for cDFT (for charge difference)",cDFT_NumberAtomGroups)
    else if(cDFT_NumberAtomGroups>1.AND.cDFT_Type==cDFT_Fix_Charge) then
       call cq_abort("Maximum of one group permitted for cDFT (for charge fixed)",cDFT_NumberAtomGroups)
    end if
    allocate(matWc(cDFT_NumberAtomGroups),cDFT_Vc(cDFT_NumberAtomGroups),cDFT_W(cDFT_NumberAtomGroups),STAT=stat)
    if(stat/=0) call cq_abort("Failure to allocate matWc and cDFT_Vc: ",cDFT_NumberAtomGroups)
    cDFT_Vc = zero
    do i=1,cDFT_NumberAtomGroups
       matWc(i) = allocate_temp_matrix(Hrange,0)
    end do
    allocate(flag_cdft_atom(ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Failure to allocate flag_cdft_atom: ",ni_in_cell)
    flag_cdft_atom = 0
    allocate(bwgrid(maxngrid,cDFT_NumberAtomGroups),STAT=stat)
    do i=1,cDFT_NumberAtomGroups
       if(inode==ionode.AND.iprint_SC>=0)  write(io_lun,fmt='(6x,"cDFT Atom Group ",i4," Target: ",f12.5)') i,cDFT_Target(i)
       do j=1,cDFT_NAtoms(i)
          flag_cdft_atom(cDFT_AtomList(i)%Numbers(j)) = i
          if(inode==ionode.AND.iprint_SC>1) write(io_lun,fmt='(4x,"Atom ",2i8)') j,cDFT_AtomList(i)%Numbers(j)
       enddo
    end do
    matHzero = allocate_temp_matrix(Hrange,0)
    return
  end subroutine init_cdft
!!***

!!****f* cdt_module/make_weights *
!!
!!  NAME 
!!   make_weights
!!  USAGE
!!   
!!  PURPOSE
!!   This subroutine abstracts the process of making the weight matrix
!!   However it's only ever called when the support functions change
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   DRB
!!  CREATION DATE
!!   2011/08
!!  MODIFICATION HISTORY
!!  
!!  SOURCE
!!  
  subroutine make_weights

    use numbers
    use mult_module, ONLY: matrix_sum, matrix_scale
    use density_module, ONLY: build_Becke_weight_matrix
    use io_module, ONLY: dump_matrix

    implicit none

    call build_Becke_weight_matrix(matWc,cDFT_NumberAtomGroups)
    call dump_matrix("matW1",matWc(1),inode)
    if(cDFT_Type==cDFT_Fix_ChargeDifference) then
       call dump_matrix("matW2",matWc(2),inode)
       if(inode==ionode.AND.iprint_SC>2) write(io_lun,fmt='(6x,"Adjusting Wc for charge difference")') 
       call matrix_sum(one,matWc(1),-one,matWc(2))
       call matrix_scale(zero,matWc(2))
       call dump_matrix("matWD",matWc(1),inode)
    end if
    
  end subroutine make_weights
!!***

!!****f* cdft_module/cdft_min *
!!
!!  NAME 
!!   cdft_min
!!  USAGE
!!   
!!  PURPOSE
!!   Finds the correct value of the constraining potential, Vc, to minimise the cDFT constraint
!!   Uses a simple root finder from Numerical Reciples
!!
!!   N.B. As it stands, this really only works for ONE constraint (charge difference between two
!!   groups counts as one).  It should be possible to implement multiple constraints in principle,
!!   but fairly extensive tests suggest that this is hard without gradients (i.e. we would require
!!   more than the simple root finder used here).
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   DRB
!!  CREATION DATE
!!   2011/08
!!  MODIFICATION HISTORY
!!  
!!  SOURCE
!!  
  subroutine cdft_min(reset_L, fixed_potential, vary_mu, &
       n_CG_L_iterations, number_of_bands, tolerance, mu, total_energy)

    use datatypes
    use numbers
    use global_module, ONLY: iprint_SC, sf,iprint_SC,io_lun
    use GenComms, ONLY: cq_abort
    use energy, ONLY: cdft_energy, get_energy
    use density_module, ONLY: get_cdft_constraint

    implicit none

    ! Passed Variables
    logical :: reset_L, fixed_potential, vary_mu
    integer :: n_CG_L_iterations
    real(double) :: number_of_bands, tolerance, mu, total_energy

    ! Local variables
    real(double),dimension(:), allocatable :: ax,bx,cx,fa,fb,fc, ax_lo, cx_hi, fa_lo,fc_hi
    real(double) :: gold,glimit,tiny,m, P, S, T
    real(double) :: dum,fu,q,r,u,ulim, inc
    integer :: i, isum, iaim, stat, ngroups, n_iter
    logical :: done
    logical, dimension(:), allocatable :: notdone

    if(cDFT_Type==cDFT_Fix_ChargeDifference) then
       ngroups = 1
    else if(cDFT_NumberAtomGroups>1) then
       call cq_abort("You can only have one cDFT group (unless doing charge differences)")
    else
       ngroups = cDFT_NumberAtomGroups
    end if
    allocate(ax(ngroups),bx(ngroups),cx(ngroups),fa(ngroups),fb(ngroups),fc(ngroups), &
         ax_lo(ngroups),cx_hi(ngroups),fa_lo(ngroups),fc_hi(ngroups),notdone(ngroups),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating ax etc in cdft_min: ",ngroups)
    isum = 0
    iaim = 2**(ngroups)-1 ! If we need more than 31 constraints, this may overflow...
    notdone = .true.
    do i=1,ngroups
       ax(i)=cDFT_Vc(i)
    end do
    reset_L = .true.
    call evaluate_cdft_function(reset_L, fixed_potential, vary_mu, &
         n_CG_L_iterations, number_of_bands, tolerance, mu, total_energy)
    do i=1,ngroups
       fa(i) = cDFT_W(i)
       if(abs(fa(i))<cDFT_Tolerance) then
          isum = isum + 2**(i-1)
          notdone(i) = .false.
       end if
    end do
    if(isum==iaim) then
       if(inode==ionode) then 
          if(iprint_SC>1) write(io_lun,fmt='(4x,"cDFT Minimum found: ")') 
          if(iprint_SC>2) then
             do i=1,ngroups
                write(io_lun,fmt='(4x,2f12.5)') ax(i),fa(i)
             end do
          end if
       end if
       deallocate(ax,bx,cx,fa,fb,fc, ax_lo,cx_hi,fa_lo,fc_hi,notdone)
       return
    end if
    do i=1,ngroups
       if(notdone(i)) then
          if(fa(i)<zero) then
             cDFT_Vc(i) = cDFT_Vc(i)-0.03675_double*abs(cDFT_Target(i))
          else
             cDFT_Vc(i) = cDFT_Vc(i)+0.03675_double*abs(cDFT_Target(i))
          end if
          cx(i)=cDFT_Vc(i)
       end if
    end do
    reset_L = .true.
    call evaluate_cdft_function(reset_L, fixed_potential, vary_mu, &
         n_CG_L_iterations, number_of_bands, tolerance, mu, total_energy)
    do i=1,ngroups
       if(notdone(i)) then
          fc(i) = cDFT_W(i)
          if(abs(fc(i))<cDFT_Tolerance) then
             isum = isum + 2**(i-1)
             notdone(i) = .false.
          end if
       end if
    end do
    if(isum==iaim) then
       if(inode==ionode) then 
          if(iprint_SC>1) write(io_lun,fmt='(4x,"cDFT Minimum found: ")') 
          if(iprint_SC>2) then
             do i=1,ngroups
                write(io_lun,fmt='(4x,2f12.5)') cx(i),fc(i)
             end do
          end if
       end if
       deallocate(ax,bx,cx,fa,fb,fc, ax_lo,cx_hi,fa_lo,fc_hi,notdone)
       return
    end if
    ! Here I need to add a real bracketing routine
    done = .false.
    if(fc(1)*fa(1)<zero) then
       done = .true.
       if(fc(1)<=fa(1))then!Switch a and b so we can go downhill
          dum=ax(1)
          ax(1)=cx(1)
          cx(1)=dum
          dum=fc(1)
          fc(1)=fa(1)
          fa(1)=dum
       end if
    end if
    n_iter = 0
    do while(.NOT.done.AND.n_iter<=cDFT_MaxIterations) 
       do i=1,ngroups
          if(notdone(i)) then
             if(fc(i)<=fa(i))then!Switch a and b so we can go downhill
                dum=ax(i)
                ax(i)=cx(i)
                cx(i)=dum
                dum=fc(i)
                fc(i)=fa(i)
                fa(i)=dum
             end if
             if(inode==ionode.AND.iprint_SC>2) write(io_lun,fmt='(6x,"cDFT: ax, bx, fa, fb: ",4f12.5)') ax(i),cx(i),fa(i),fc(i)
             m = (fc(i)-fa(i))/(cx(i)-ax(i))
             bx(i) = -(fc(i)-m*cx(i))/m
             if(inode==ionode.AND.iprint_SC>2) write(io_lun,fmt='(6x,"cDFT: a,b, prediction: ",6f12.5)') &
                  ax(i),cx(i),fa(i),fc(i),m,bx(i)
             if(abs(bx(i)-cx(i))>four*abs(cx(i)-ax(i))) then
                bx(i) = ax(i) + four*(cx(i)-ax(i))
                if(inode==ionode.AND.iprint_SC>2) write(io_lun,fmt='(6x,"cDFT: bx too large; adjusted to: ",f12.5)') bx(i)
             end if
             cDFT_Vc(i) = bx(i)
          end if
       end do
       reset_L = .true.
       call evaluate_cdft_function(reset_L, fixed_potential, vary_mu, &
            n_CG_L_iterations, number_of_bands, tolerance, mu, total_energy)
       do i=1,ngroups
          if(notdone(i)) then
             fb(i) = cDFT_W(i)
             if(abs(fb(i))<cDFT_Tolerance) then
                isum = isum + 2**(i-1)
                notdone(i) = .false.
             end if
          end if
       end do
       fa(1) = fc(1)
       ax(1) = cx(1)
       fc(1) = fb(1)
       cx(1) = bx(1)
       if(fb(1)*fa(1)<zero.OR.isum==iaim) then
          done = .true.
       end if
       n_iter = n_iter+1
    end do
    if(n_iter>=cDFT_MaxIterations) call cq_abort("Exceeded cDFT.MaxIterations: ",n_iter)
    if(inode==ionode.AND.iprint_SC>2) write(io_lun,fmt='(6x,"cDFT done bracketing: ",4f12.5)') ax(1),cx(1),fa(1),fc(1)
    if(isum==iaim) then
       if(inode==ionode) then 
          if(iprint_SC>1) write(io_lun,fmt='(4x,"cDFT Minimum found: ")') 
          if(iprint_SC>2) then
             do i=1,ngroups
                write(io_lun,fmt='(4x,2f12.5)') bx(i),fb(i)
             end do
          end if
       end if
       deallocate(ax,bx,cx,fa,fb,fc, ax_lo,cx_hi,fa_lo,fc_hi,notdone)
       return
    end if
    n_iter = 0
    do while(isum/=iaim.AND.n_iter<=cDFT_MaxIterations)!abs(fb)>cDFT_Tolerance)
       bx(1) = half*(ax(1)+cx(1))
       cDFT_Vc(1) = bx(1)
       reset_L = .true.
       call evaluate_cdft_function(reset_L, fixed_potential, vary_mu, &
            n_CG_L_iterations, number_of_bands, tolerance, mu, total_energy)
       do i=1,ngroups
          if(notdone(i)) then
             fb(i) = cDFT_W(i)
             if(abs(fb(i))<cDFT_Tolerance) then
                isum = isum + 2**(i-1)
                notdone(i) = .false.
             end if
          end if
       end do
       if(inode==ionode.AND.iprint_SC>2) write(io_lun,fmt='(6x,"cDFT: bisection: ",2f12.5)') bx(1),fb(1)
       if(isum==iaim) exit
       S = sqrt(fb(1)*fb(1)-fa(1)*fc(1))
       if(fa(1)-fc(1)<zero) then
          T = bx(1) - (bx(1)-ax(1))*fb(1)/S
       else
          T = bx(1) + (bx(1)-ax(1))*fb(1)/S
       end if
       cDFT_Vc(1) = T
       reset_L = .true.
       call evaluate_cdft_function(reset_L, fixed_potential, vary_mu, &
            n_CG_L_iterations, number_of_bands, tolerance, mu, total_energy)
       do i=1,ngroups
          if(notdone(i)) then
             fc_hi(i) = cDFT_W(i)
             if(abs(fc_hi(i))<cDFT_Tolerance) then
                isum = isum + 2**(i-1)
                notdone(i) = .false.
             end if
          end if
       end do
       if(inode==ionode.AND.iprint_SC>2) write(io_lun,fmt='(6x,"cDFT: prediction: ",2f12.5)') T,fc_hi(1)
       if(isum==iaim) then
          bx(1) = T
          fb(1) = fc_hi(1)
          exit
       end if
       if(fc_hi(1)*fb(1)<zero) then
          ax(1)=bx(1)
          fa(1)=fb(1)
          cx(1)=T
          fc(1)=fc_hi(1)
       else if(fc_hi(1)*fa(1)<zero) then
          cx(1)=T
          fc(1)=fc_hi(1)
       else if(fc_hi(1)*fc(1)<zero) then
          ax(1)=T
          fa(1)=fc_hi(1)
       else
          call cq_abort("Error !")
       end if
       n_iter = n_iter+1
    end do
    if(n_iter>=cDFT_MaxIterations) call cq_abort("Exceeded cDFT.MaxIterations: ",n_iter)
    if(inode==ionode) then 
       if(iprint_SC>1) write(io_lun,fmt='(4x,"cDFT Minimum found: ")') 
       if(iprint_SC>2) then
          do i=1,ngroups
             write(io_lun,fmt='(4x,2f12.5)') bx(i),fb(i)
          end do
       end if
    end if
    deallocate(ax,bx,cx,fa,fb,fc, ax_lo,cx_hi,fa_lo,fc_hi,notdone)
    return
  end subroutine cdft_min
!!***

!!****f* cdft_module/evaluate_cdft_function *
!!
!!  NAME 
!!   evaluate_cdft_function
!!  USAGE
!!   
!!  PURPOSE
!!   Evaluates the cDFT energy functional by finding ground state DM after applying potential
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   DRB
!!  CREATION DATE
!!   2011/08
!!  MODIFICATION HISTORY
!!  
!!  SOURCE
!!  
   subroutine evaluate_cdft_function( reset_L, fixed_potential, vary_mu, &
        n_CG_L_iterations, number_of_bands, tolerance, mu, total_energy)

     use datatypes
     use numbers
     use logicals
     use mult_module, ONLY: LNV_matrix_multiply, matH, matrix_sum
     use DMMin, ONLY: FindMinDM
     use global_module, ONLY: iprint_SC, sf,iprint_SC,io_lun
     use DiagModule, ONLY: diagon
     use energy, ONLY: get_energy
     use GenComms, ONLY: inode, ionode
     use density_module, ONLY: get_cdft_constraint

     implicit none

     ! Passed Variables
     logical ::  reset_L, fixed_potential, vary_mu

     integer :: n_CG_L_iterations

     real(double) :: number_of_bands, tolerance, mu, total_energy

     ! Local variables
     real(double) :: electrons, total_energy_1, start_BE, new_BE, Ltol
     integer :: i, temp_supp_fn

     ! Change the Hamiltonian
     call matrix_sum(zero,matH,one,matHzero)
     do i=1,cDFT_NumberAtomGroups
        call matrix_sum(one,matH,cDFT_Vc(i),matWc(i))
        if(inode==ionode.AND.iprint_SC>2) write(io_lun,fmt='(4x,"Group ",i4," Vc ",f12.5)') i,cDFT_Vc(i)
     end do
     ! Find minimum density matrix
     Ltol = tolerance
     !reset_L = .true.
     call FindMinDM(n_CG_L_iterations, number_of_bands, vary_mu, Ltol, mu, inode, ionode, reset_L, .false.)
     ! If we're using O(N), we only have L, and we need K - if diagonalisation, we have K
     if(.NOT.diagon) call LNV_matrix_multiply(electrons, total_energy,  &
          doK, dontM1, dontM2, dontM3, dontM4, dontphi, dontE,0,0,0,0)
     ! Get total energy
     call matrix_sum(zero,matH,one,matHzero)
     call get_energy(total_energy)
     call get_cdft_constraint
     return
   end subroutine evaluate_cdft_function
!!***

end module cdft_module


