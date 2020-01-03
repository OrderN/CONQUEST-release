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
!!   2012/04/21 L.Tong
!!   - Moved get_cdft_constraint from density module to here
!!   2019/10/24 11:52 dave
!!    Changed function calls to FindMinDM
!!  SOURCE
!!
module cdft_module

  use datatypes
  use cdft_data
  use global_module, only: io_lun, iprint_SC
  use GenComms, only: inode, ionode

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
  !!   2011/12/10 L.Tong
  !!     Added spin polarisation: initialisation for matHzero_dn
  !!   2012/03/18 L.Tong
  !!     Changed spin implementation
  !!  SOURCE
  !!  
  subroutine init_cdft
    
    use datatypes
    use numbers
    use global_module,  only: ni_in_cell, flag_cdft_atom, &
                              nspin
    use mult_module,    only: allocate_temp_matrix
    use matrix_data,    only: Hrange
    use density_module, only: bwgrid
    use maxima_module,  only: maxngrid
    use GenComms,       only: cq_abort

    implicit none

    integer :: i, j, stat, spin

    if (cDFT_NumberAtomGroups > 2) then
       call cq_abort("Maximum of two groups permitted for cDFT (for &
                      &charge difference)", cDFT_NumberAtomGroups)
    else if (cDFT_NumberAtomGroups > 1 .AND. &
         cDFT_Type == cDFT_Fix_Charge) then
       call cq_abort("Maximum of one group permitted for cDFT (for &
                      &charge fixed)", cDFT_NumberAtomGroups)
    end if
    allocate (matWc(cDFT_NumberAtomGroups), &
         cDFT_Vc(cDFT_NumberAtomGroups), &
         cDFT_W(cDFT_NumberAtomGroups), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Failure to allocate matWc and cDFT_Vc: ", &
                       cDFT_NumberAtomGroups)
    cDFT_Vc = zero
    do i = 1, cDFT_NumberAtomGroups
       matWc(i) = allocate_temp_matrix(Hrange, 0)
    end do
    allocate(flag_cdft_atom(ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Failure to allocate flag_cdft_atom: ", &
                       ni_in_cell)
    flag_cdft_atom = 0
    allocate(bwgrid(maxngrid, cDFT_NumberAtomGroups), STAT=stat)
    do i = 1, cDFT_NumberAtomGroups
       if (inode == ionode .AND. iprint_SC >= 0) &
            write (io_lun, fmt='(6x,"cDFT Atom Group ",i4," Target: ",f12.5)') &
                  i, cDFT_Target(i)
       do j = 1, cDFT_NAtoms(i)
          flag_cdft_atom(cDFT_AtomList(i)%Numbers(j)) = i
          if (inode == ionode .AND. iprint_SC > 1) &
               write (io_lun, fmt='(4x,"Atom ",2i8)') &
                      j, cDFT_AtomList(i)%Numbers(j)
       enddo
    end do
    allocate(matHzero(nspin), STAT=stat)
    if (stat /= 0) call cq_abort('init_cdft: failed to allocate matHzero', stat)
    do spin = 1, nspin
       matHzero(spin) = allocate_temp_matrix(Hrange, 0)
    end do

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
  !!   2013/03/06 17:00 dave
  !!   - Removing matrix dumping
  !!  SOURCE
  !!  
  subroutine make_weights

    use numbers
    use mult_module, only: matrix_sum, matrix_scale
    use density_module, only: build_Becke_weight_matrix
    use io_module, only: dump_matrix

    implicit none

    call build_Becke_weight_matrix(matWc,cDFT_NumberAtomGroups)
    !call dump_matrix("matW1",matWc(1),inode)
    if(cDFT_Type==cDFT_Fix_ChargeDifference) then
       !call dump_matrix("matW2",matWc(2),inode)
       if(inode==ionode.AND.iprint_SC>2) write(io_lun,fmt='(6x,"Adjusting Wc for charge difference")') 
       call matrix_sum(one,matWc(1),-one,matWc(2))
       call matrix_scale(zero,matWc(2))
       !call dump_matrix("matWD",matWc(1),inode)
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
  !!   Finds the correct value of the constraining potential, Vc, to
  !!   minimise the cDFT constraint Uses a simple root finder from
  !!   Numerical Reciples
  !!
  !!   N.B. As it stands, this really only works for ONE constraint
  !!   (charge difference between two groups counts as one).  It should
  !!   be possible to implement multiple constraints in principle, but
  !!   fairly extensive tests suggest that this is hard without
  !!   gradients (i.e. we would require more than the simple root finder
  !!   used here).
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
  !!   2011/12/10 L.Tong
  !!     Removed redundant parameter number_of_bands
  !!   2012/03/18 L.Tong
  !!   - Removed redundant input parameter real(double) mu
  !!   2013/03/06 17:00 dave
  !!   - Normalise step to number of constrained atoms
  !!   2016/08/08 15:30 nakata
  !!   - Removed unused sf in global_module
  !!  SOURCE
  !!  
  subroutine cdft_min(reset_L, fixed_potential, vary_mu, &
                      n_CG_L_iterations, tolerance, total_energy)

    use datatypes
    use numbers
    use global_module,  only: iprint_SC, iprint_SC,io_lun
    use GenComms,       only: cq_abort
    use energy,         only: cdft_energy, get_energy

    implicit none

    ! Passed Variables
    logical      :: reset_L, fixed_potential, vary_mu
    integer      :: n_CG_L_iterations
    real(double) :: tolerance, total_energy

    ! Local variables
    real(double) :: gold, glimit, tiny, m, P, S, T
    real(double) :: dum, fu, q, r, u, ulim, inc
    integer      :: i, isum, iaim, stat, ngroups, n_iter
    logical      :: done
    logical,      dimension(:), allocatable :: notdone
    real(double), dimension(:), allocatable :: ax, bx, cx, fa, fb, fc,&
                                               ax_lo, cx_hi, fa_lo,fc_hi

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
         n_CG_L_iterations, tolerance, total_energy)
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
             cDFT_Vc(i) = cDFT_Vc(i)-0.03675_double*abs(cDFT_Target(i))/real(cDFT_NAtoms(i),double)
          else
             cDFT_Vc(i) = cDFT_Vc(i)+0.03675_double*abs(cDFT_Target(i))/real(cDFT_NAtoms(i),double)
          end if
          cx(i)=cDFT_Vc(i)
       end if
    end do
    reset_L = .true.
    call evaluate_cdft_function(reset_L, fixed_potential, vary_mu, &
         n_CG_L_iterations, tolerance, total_energy)
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
            n_CG_L_iterations, tolerance, total_energy)
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
            n_CG_L_iterations, tolerance, total_energy)
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
            n_CG_L_iterations, tolerance, total_energy)
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
  !!  PURPOSE
  !!   Evaluates the cDFT energy functional by finding ground state DM
  !!   after applying potential
  !!  INPUTS
  !!  USES
  !!  AUTHOR
  !!   DRB
  !!  CREATION DATE
  !!   2011/08
  !!  MODIFICATION HISTORY
  !!   2011/12/09 L.Tong
  !!     - Removed redundant parameter number_of_bands
  !!     - Added spin polarisation: we assume constraint potantial is
  !!       the same for both spin components (components in the direct
  !!       sum)
  !!   2012/03/18 L.Tong
  !!     - Rewrite for change in spin implementation
  !!     - Removed redundant input parameter real(double) mu
  !!   2016/08/08 nakata
  !!     - Removed unused sf in global_module
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!  SOURCE
  !!  
  subroutine evaluate_cdft_function(reset_L, fixed_potential,   &
                                    vary_mu, n_CG_L_iterations, &
                                    tolerance, total_energy)

     use datatypes
     use numbers
     use logicals
     use mult_module,    only: LNV_matrix_multiply, matH, matrix_sum
     use DMMin,          only: FindMinDM
     use global_module,  only: iprint_SC, iprint_SC, io_lun, &
                               nspin, spin_factor, flag_diagonalisation
     !use DiagModule,     only: diagon
     use energy,         only: get_energy
     use GenComms,       only: inode, ionode

     implicit none

     ! Passed Variables
     logical ::  reset_L, fixed_potential, vary_mu
     integer :: n_CG_L_iterations
     real(double) :: tolerance, total_energy

     ! Local variables
     real(double), dimension(nspin) :: electrons, tmp_energy
     real(double) :: start_BE, new_BE, Ltol
     integer :: i, temp_supp_fn, spin

     ! Change the Hamiltonian
     do spin = 1, nspin
        call matrix_sum(zero, matH(spin), one, matHzero(spin))
     end do
     
     do i = 1, cDFT_NumberAtomGroups
        do spin = 1, nspin
           call matrix_sum(one, matH(spin), cDFT_Vc(i), matWc(i))
        end do
        if (inode == ionode .AND. iprint_SC > 2) &
             write (io_lun,fmt='(4x,"Group ",i4," Vc ",f12.5)') &
                   i, cDFT_Vc(i)
     end do
     ! Find minimum density matrix
     Ltol = tolerance
     !reset_L = .true.
     call FindMinDM(n_CG_L_iterations, vary_mu, Ltol, &
                    reset_L, .false.)
     ! If we're using O(N), we only have L, and we need K - if
     ! diagonalisation, we have K
     if (.not. flag_diagonalisation) then
        ! electrons and tmp_energy are used as dump for electrons and
        ! energies calculated from LNV_matrix_multipy
        ! total_energy is calculated from get_energy below
        call LNV_matrix_multiply(electrons, tmp_energy, doK, dontM1, &
                                 dontM2, dontM3, dontM4, dontphi, dontE)
     end if
     ! Get total energy
     do spin = 1, nspin
        call matrix_sum(zero, matH(spin), one, matHzero(spin))
     end do
     call get_energy(total_energy)
     call get_cdft_constraint
     return
   end subroutine evaluate_cdft_function
   !!***


   !!****f* cdft_module/get_cdft_constraint *
   !!
   !!  NAME 
   !!   get_cdft_constraint
   !!  USAGE
   !!   
   !!  PURPOSE
   !!   Gets constraint for cDFT
   !!  INPUTS
   !!   
   !!  USES
   !!   
   !!  AUTHOR
   !!   DRB and Alex Sena
   !!  CREATION DATE
   !!   2009
   !!  MODIFICATION HISTORY
   !!   2011/08 and 2011/09
   !!    Incorporated into new trunk
   !!   2011/12/10 L.Tong
   !!    Removed redundant dependence on matHzero from cdft_data module
   !!   2012/03/13 L.Tong
   !!    Added spin polarisation
   !!   2013/03/06 17:01 dave
   !!   - Added output to show charge on group
   !!  SOURCE
   !!  
   subroutine get_cdft_constraint

     use numbers
     use mult_module,   only: matH, matK, matrix_sum,               &
                              matrix_product_trace
     use cdft_data,     only: cDFT_Type, cDFT_Fix_Charge,           &
                              cDFT_Fix_ChargeDifference,            &
                              cDFT_NumberAtomGroups, cDFT_W, matWc, &
                              cDFT_Target, cDFT_Vc
     use energy,        only: cdft_energy
     use GenComms,      only: inode, ionode
     use global_module, only: nspin, spin_factor

     implicit none
     
     integer :: i, spin
     real(double) :: group_charge

     cdft_energy = zero
     if (cDFT_Type == cDFT_Fix_Charge .or. &
          cDFT_Type == cDFT_Fix_ChargeDifference) then
        do i = 1, cDFT_NumberAtomGroups
           cDFT_W(i) = zero
           do spin = 1, nspin
              group_charge = spin_factor * &
                          matrix_product_trace(matK(spin), matWc(i))
              cDFT_W(i) = cDFT_W(i) + (group_charge - cDFT_Target(i))
              if(inode == ionode .and. iprint_SC > 2) write(io_lun,fmt='(4x,"Group ",i3,"Charge: ",f20.12)') i,group_charge
           end do
           cdft_energy = cdft_energy + cDFT_Vc(i) * cDFT_W(i)
           if (inode == ionode .and. iprint_SC > 2) &
                write (io_lun, fmt='(4x,"Group ",i3,"Vc, W, E: ",3f20.12)') &
                i, cDFT_Vc(i), cDFT_W(i), cDFT_Vc(i) * cDFT_W(i)
        end do
     end if
     ! Other cDFT types go above
     return
   end subroutine get_cdft_constraint
   !!***


end module cdft_module


