module schro

  use datatypes
  use pseudo_atom_info
  use global_module, ONLY: iprint
  
  implicit none

contains

  subroutine make_paos(i_species)

    use datatypes
    use numbers
    use pseudo_tm_info, ONLY: pseudo
    use mesh, ONLY: nmesh, make_mesh, rr
    use global_module, ONLY: iprint
    use GenComms, ONLY: cq_abort
    use radial_xc, ONLY: get_vxc
    use input_module, ONLY: io_assign, io_close
    use periodic_table, ONLY: pte
    use units, ONLY: HaToeV
    
    implicit none

    ! Passed variables
    integer :: i_species

    ! Local variables
    real(double), allocatable, dimension(:) :: vha, vxc, atomic_rho, vha_conf
    
    write(*,fmt='(/"Generating PAOs for ",a2/)') pte(nint(pseudo(i_species)%z))
    !
    ! Create mesh
    !
    nmesh = local_and_vkb%ngrid
    call make_mesh(i_species)
    !
    ! Allocate and zero
    !
    allocate(vha(nmesh), vxc(nmesh), vha_conf(nmesh), atomic_rho(nmesh))
    atomic_rho = zero
    vha = zero
    vxc = zero
    ! 
    ! Potentials: Hartree and XC; add and remove PCC charge if used
    !
    call radial_hartree(nmesh,local_and_vkb%charge,vha)
    if(pseudo(i_species)%flag_pcc) local_and_vkb%charge = local_and_vkb%charge + local_and_vkb%pcc
    call get_vxc(nmesh,rr,local_and_vkb%charge,vxc)
    if(pseudo(i_species)%flag_pcc) local_and_vkb%charge = local_and_vkb%charge - local_and_vkb%pcc
    !
    ! Find energies of valence states without confinement
    !
    call find_unconfined_valence_states(i_species,vha,vxc)
    !
    ! Find default radii or radii from cutoffs
    !
    call set_radii(i_species,vha,vxc)
    !
    ! Solve for PAOs
    !
    write(*,fmt='(/2x,"Solving for PAOs")')
    if(paos%flag_zetas==1.and.iprint>2) write(*,fmt='(2x,"Using split-norm approach for multiple zetas")')
    call solve_for_occupied_paos(i_species,vha,vxc, atomic_rho)
    if(paos%n_shells>val%n_occ) &
         call solve_for_polarisation(i_species,vha,vxc)
    !
    ! Interpolation
    ! 
    write(*,fmt='(/2x,"Interpolating onto regular mesh")')
    !
    ! Interpolate PAOs
    !
    call interpolate_paos(i_species)
    !
    ! Interpolate potentials
    !
    write(*,fmt='(/2x,"Calculating potentials")')
    ! Build the neutral atom potential as sum of local and pseudo-atomic hartree potentials
    call radial_hartree(nmesh,atomic_rho,vha_conf)
    vha_conf = vha_conf + local_and_vkb%local
    call interpolate_potentials(i_species,vha_conf)
    deallocate(vha,vxc,vha_conf,atomic_rho)
    write(*,fmt='(/2x,"Finished ",a2)') pte(nint(pseudo(i_species)%z))
  end subroutine make_paos

  ! Solve for the unconfined (atomic) pseudo-functions
  subroutine find_unconfined_valence_states(i_species,vha,vxc)

    use datatypes
    use numbers
    use mesh, ONLY: rr, nmesh, drdi_squared, rr_squared,drdi
    use radial_xc, ONLY: get_vxc
    use read, ONLY: ps_format, hgh, alpha_scf, max_scf_iters
    use pseudo_tm_info, ONLY: pseudo
    use GenComms, ONLY: cq_abort
    
    implicit none

    ! Passed variables
    integer :: i_species
    real(double), dimension(nmesh) :: vha, vxc

    ! Local variables
    integer :: i_shell, ell, en, i, iter, maxiter
    real(double) :: radius_large, large_energy, resid, total_charge, check
    real(double), allocatable, dimension(:) :: psi, newcharge

    allocate(psi(nmesh),newcharge(nmesh))
    psi = zero
    newcharge = zero
    maxiter = max_scf_iters
    iter = 0
    if(iprint>2) write(*,fmt='(/2x,"Finding unconfined energies for valence states")')
    resid = one
    do while(resid>1e-12_double.and.iter<maxiter)
       iter = iter+1
       newcharge = zero
       do i_shell = 1, val%n_occ
          if(iprint>2) write(*,fmt='(2x,"n=",i2," l=",i2)') val%n(i_shell), val%l(i_shell)
          radius_large = rr(nmesh)
          ell = val%l(i_shell)
          en = val%npao(i_shell)
          large_energy = val%en_ps(i_shell)
          if(ps_format==hgh) then
             if(resid<0.01_double) then
                large_energy = val%en_ps(i_shell)
             else
                large_energy = zero!half*val%en_ps(i_shell)
             end if
          end if
          !if(abs(large_energy)<RD_ERR.and.i_shell>1) large_energy = val%en_pao(i_shell-1)
          call find_eigenstate_and_energy_vkb(i_species,en,ell,radius_large, psi,large_energy,vha,vxc)
          val%en_pao(i_shell) = large_energy
          ! Accumulate output charge
          if(ell==0) then
             newcharge = newcharge + val%occ(i_shell)*psi*psi
          else if(ell==1) then
             newcharge = newcharge + val%occ(i_shell)*psi*psi*rr*rr
          else if(ell==2) then
             newcharge = newcharge + val%occ(i_shell)*psi*psi*rr*rr*rr*rr
          end if
       end do
       ! Find residual
       resid = zero
       check = zero
       do i=1,nmesh
          resid = resid + drdi(i)*rr_squared(i)*(local_and_vkb%charge(i)-newcharge(i))**2
          check = check + drdi(i)*rr_squared(i)*local_and_vkb%charge(i)
       end do
       local_and_vkb%charge = alpha_scf*newcharge + (one-alpha_scf)*local_and_vkb%charge
       ! We can use these lines to write out the charge for future solvers
       !open(unit=70,file='charge_out.dat',position="append")
       !do i=1,nmesh
       !   write(70,*) rr(i), newcharge(i)!local_and_vkb%charge(i)
       !end do
       !close(70)
       ! Integrate
       total_charge = zero
       check = zero
       do i=1,nmesh
          resid = resid + drdi(i)*rr_squared(i)*(local_and_vkb%charge(i)-newcharge(i))**2
          total_charge = total_charge + drdi(i)*rr_squared(i)*newcharge(i)
          check = check + drdi(i)*rr_squared(i)*local_and_vkb%charge(i)
       end do
       ! This rescales charge to have full valence value, but can be unstable
       !if(abs(check-total_charge)>RD_ERR) local_and_vkb%charge = local_and_vkb%charge*total_charge/check
       if(iprint>2) write(*,fmt='("Iteration ",i4," residual ",e14.6)') iter,resid
       if(iprint>0) then
          write(*,fmt='(/2x,"Unconfined valence state energies (Ha)")')
          write(*,fmt='(2x,"  n  l         AE energy        PAO energy")')
          do i_shell = 1, val%n_occ
             ell = val%l(i_shell)
             en = val%n(i_shell)
             write(*,fmt='(2x,2i3,2f18.10)') en, ell, val%en_ps(i_shell), &
                  val%en_pao(i_shell)
             if(ps_format==hgh) val%en_ps(i_shell) = val%en_pao(i_shell)
          end do
       else if(ps_format==hgh) then
          do i_shell = 1, val%n_occ
             val%en_ps(i_shell) = val%en_pao(i_shell)
          end do
       end if
       call radial_hartree(nmesh,local_and_vkb%charge,vha)
       if(pseudo(i_species)%flag_pcc) local_and_vkb%charge = local_and_vkb%charge + local_and_vkb%pcc
       call get_vxc(nmesh,rr,local_and_vkb%charge,vxc)
       if(pseudo(i_species)%flag_pcc) local_and_vkb%charge = local_and_vkb%charge - local_and_vkb%pcc
    end do
    if(iter>maxiter) call cq_abort("Exceeded iterations in SCF")
    ! We can use these lines to write out the charge for future solvers
    open(unit=70,file='charge_out.dat')
    do i=1,nmesh
       write(70,*) rr(i), newcharge(i)!local_and_vkb%charge(i)
    end do
    close(70)
    return
  end subroutine find_unconfined_valence_states

  ! Set radii of PAOs following user settings
  subroutine set_radii(i_species,vha,vxc)

    use datatypes
    use mesh, ONLY: rr, delta_r_reg, nmesh
    use units, ONLY: HaToeV
    
    implicit none

    ! Passed variables
    integer :: i_species
    real(double), dimension(nmesh) :: vha, vxc

    ! Local variables
    integer :: i_shell, ell, en, zeta, grid_size
    
    write(*,fmt='(/2x,"Setting cutoff radii"/)')
    if(flag_default_cutoffs) then ! Work out radii
       if(paos%flag_cutoff==pao_cutoff_energies.OR.paos%flag_cutoff==pao_cutoff_default) then
          write(*,fmt='(4x,"Default cutoffs, shells share energy shifts")')
       else if(paos%flag_cutoff==pao_cutoff_radii) then
          write(*,fmt='(4x,"Default cutoffs, averaging radii over shells")')
       end if
       call find_default_cutoffs(i_species,vha,vxc)
    else if(paos%flag_cutoff==pao_cutoff_energies) then ! Work out radii from energy shifts
       write(*,fmt='(4x,"User-specified energy shifts")')
       write(*,fmt='(4x,"  n  l  z  delta E (Ha)  delta E (eV)")')
       do i_shell = 1,val%n_occ
          ell = paos%l(i_shell)
          en = paos%n(i_shell)
          do zeta = 1,paos%nzeta(i_shell)
             write(*,fmt='(4x,3i3,1x,2f13.8)') en, ell, zeta, paos%energy(zeta,i_shell), &
                  paos%energy(zeta,i_shell)*HaToeV
             en = paos%npao(i_shell)
             call find_radius_from_energy(i_species,en, ell, paos%cutoff(zeta,i_shell), &
                  val%en_ps(i_shell) + paos%energy(zeta,i_shell), vha, vxc, .false.)
             ! Set cutoff to the first regular mesh point BEYOND the actual cutoff
             grid_size = 1+floor(paos%cutoff(zeta,i_shell)/delta_r_reg)
             paos%cutoff(zeta,i_shell) = delta_r_reg*real(grid_size,double)
          end do
       end do
       do i_shell = val%n_occ+1,paos%n_shells
          if(i_shell>val%n_occ+1) write(*,fmt='(2x,"Dangerous to set unoccupied shell radii via energy !")')
          if(paos%inner(i_shell)>0) then
             do zeta = 1,paos%nzeta(i_shell)
                paos%cutoff(zeta,i_shell) = paos%cutoff(zeta,paos%inner(i_shell))
             end do
          else
             do zeta = 1,paos%nzeta(i_shell)
                paos%cutoff(zeta,i_shell) = paos%cutoff(zeta,i_shell-1)
             end do
          end if
       end do
    else ! Radii
       write(*,fmt='(4x,"User-specified radii")')          
    end if
    write(*,fmt='(/4x,"Cutoff radii for PAOs")')
    write(*,fmt='(4x,"  n  l  z      R (bohr)")')
    do i_shell = 1,paos%n_shells
       ell = paos%l(i_shell)
       en = paos%n(i_shell)
       if(paos%flag_perturb_polarise.AND.i_shell==paos%n_shells) then
          if(en<3) en = en+1
          do zeta = 1,paos%nzeta(i_shell)
             write(*,fmt='(4x,3i3,f12.4)') en,ell,zeta,paos%cutoff(zeta,i_shell)
          end do
       else
          do zeta = 1,paos%nzeta(i_shell)
             write(*,fmt='(4x,3i3,f12.4)') en,ell,zeta,paos%cutoff(zeta,i_shell)
          end do
       end if
    end do
    return
  end subroutine set_radii

  ! Find the occupied (non-polarisation) PAOs
  !
  ! Changes
  !
  ! 2019/07/15 11:33 dave
  ! 2021/09/27 16:31 dave
  !  Restore orthogonalisation to semi-core states and correctly
  !  normalise resulting functions.
  ! 2021/09/28 14:38 dave
  !  Remove orthogonalisation (!) because it breaks Ba perturbative
  !  polarisation (see below)
  subroutine solve_for_occupied_paos(i_species,vha,vxc,atomic_density)

    use datatypes
    use numbers, ONLY: zero, RD_ERR
    use mesh, ONLY: rr, delta_r_reg, drdi, nmesh
    use units, ONLY: HaToeV
    
    implicit none

    ! Passed variables
    integer :: i_species
    real(double), dimension(nmesh) :: vha, vxc, atomic_density

    ! Local variables
    integer :: i, i_shell, ell, en, zeta
    real(double) :: large_energy, dot_p
    real(double), allocatable, dimension(:) :: psi

    allocate(psi(nmesh))
    psi = zero
    do i_shell = 1,val%n_occ
       ell = paos%l(i_shell)
       en = paos%npao(i_shell)
       do zeta = 1,paos%nzeta(i_shell)
          allocate(paos%psi(zeta,i_shell)%f(nmesh))
          paos%psi(zeta,i_shell)%f = zero
          if(iprint>2) write(*,fmt='(2x,"Species ",i2," n=",i2," l=",i2," zeta=",i2, " Rc=",f4.1," bohr")') &
               i_species, en, ell, zeta, paos%cutoff(zeta,i_shell)
          if(zeta>1.AND.paos%flag_zetas==1) then
             if(iprint>2) write(*,fmt='(2x,"Using split-norm approach for multiple zetas")')
             call find_split_norm(en,ell,paos%cutoff(zeta,i_shell),&
                  paos%psi(zeta,i_shell)%f,paos%psi(1,i_shell)%f)
          else
             large_energy = val%en_ps(i_shell) + paos%energy(zeta,i_shell)
             call find_eigenstate_and_energy_vkb(i_species,en,ell,paos%cutoff(zeta,i_shell),&
                  paos%psi(zeta,i_shell)%f,large_energy, vha,vxc,&
                  paos%width(i_shell),paos%prefac(i_shell))
             if(iprint>2) write(*,fmt='(2x,"Final energy shift required: ",f10.5," Ha")') &
                  large_energy - val%en_ps(i_shell)
             paos%energy(zeta,i_shell) = large_energy
             ! Accumulate atomic density
             if(zeta==1.AND.i_shell<=val%n_shells) then
                psi = paos%psi(zeta,i_shell)%f
                ! Scale by r^l
                if(ell>0) then
                   do i=1,ell
                      psi = psi * rr
                   end do
                end if
                ! Accumulate atomic charge density
                atomic_density = atomic_density + val%occ(i_shell)*psi*psi
                !do i=1,nmesh
                !   write(50+ell,*) rr(i),psi(i)*psi(i)*val%occ(i_shell)
                !end do
                !flush(50+ell)
             end if
             ! These lines orthogonalise to semi-core states where necessary
             ! They are left for completeness, but I found that they can cause
             ! problems with the perturbative polarisation for Ba (basically the
             ! orthogonalisation means that the 6s shell which we perturb is not
             ! negative near r=0 which breaks the sign-change detection in the
             ! perturbative polarisation routines and means they fail).  If we can
             ! sort this out, we could restore this, but I think that it's not
             ! very important.  Dave Bowler, 2021/09/28 14:38
             !if(paos%inner(i_shell)>0) then
             !   ! Dot product of two
             !   dot_p = zero
             !   do i=1,nmesh
             !      dot_p = dot_p + rr(i)**(2*ell+2)*paos%psi(zeta,i_shell)%f(i)* &
             !           paos%psi(1,val%inner(i_shell))%f(i)*drdi(i)
             !   end do
             !   if(iprint>2) write(*,fmt='(2x,"Orthogonalising to semi-core; overlap is ",f10.5)') dot_p
             !   ! Orthogonalise
             !   paos%psi(zeta,i_shell)%f = paos%psi(zeta,i_shell)%f - &
             !        dot_p * paos%psi(1,val%inner(i_shell))%f
             !   ! Normalise
             !   dot_p = zero
             !   do i=1,nmesh
             !      dot_p = dot_p + rr(i)**(2*ell+2)*paos%psi(zeta,i_shell)%f(i)* &
             !           paos%psi(zeta,i_shell)%f(i)*drdi(i)
             !   end do
             !   if(iprint>2) write(*,fmt='(2x,"Normalising: ",f10.5)') dot_p
             !   paos%psi(zeta,i_shell)%f = paos%psi(zeta,i_shell)%f/sqrt(dot_p)
             !end if ! Inner shell orthogonalisation
          end if ! Split-norm or confined state
       end do ! zeta = 1, paos%nzeta
    end do ! i_shell = 1,paos%n_shells
    deallocate(psi)
    return
  end subroutine solve_for_occupied_paos

  ! Solve for polarisation (unoccupied) shells
  !
  ! Changes
  !
  ! 2019/07/15 11:33 dave
  ! 2021/09/27 16:31 dave
  !  Restore orthogonalisation to semi-core states and correctly
  !  normalise resulting functions.  Particularly important for
  !  elements such as Ga or Ge with semi-core d shell and l=2
  !  polarisation shell
  subroutine solve_for_polarisation(i_species,vha,vxc)

    use datatypes
    use numbers, ONLY: zero, one
    use mesh, ONLY: rr, delta_r_reg, drdi, nmesh
    use units, ONLY: HaToeV
    
    implicit none

    ! Passed variables
    integer :: i_species
    real(double), dimension(nmesh) :: vha, vxc

    ! Local variables
    integer :: i, i_shell, ell, en, zeta, i_end
    real(double) :: large_energy, dot_p

    !
    ! If perturbative polarisation is selected, it is the last shell
    !
    if(paos%flag_perturb_polarise) then
       i_end = paos%n_shells-1
    else
       i_end = paos%n_shells
    end if
    !
    ! Empty shells, solved for confined atom
    !
    do i_shell = val%n_occ+1,i_end
       ell = paos%l(i_shell)
       en = paos%npao(i_shell)
       do zeta = 1,paos%nzeta(i_shell)
          allocate(paos%psi(zeta,i_shell)%f(nmesh))
          paos%psi(zeta,i_shell)%f = zero
          if(iprint>2) write(*,fmt='(2x,"Species ",i2," n=",i2," l=",i2," zeta=",i2, " Rc=",f4.1," bohr pol")') &
               i_species, en, ell, zeta, paos%cutoff(zeta,i_shell)!, &
          if(zeta>1.AND.paos%flag_zetas==1) then
             if(iprint>2) write(*,fmt='(2x,"Using split-norm approach for multiple zetas")')
             call find_split_norm(en,ell,paos%cutoff(zeta,i_shell),&
                  paos%psi(zeta,i_shell)%f,paos%psi(1,i_shell)%f)
          else
             paos%energy(zeta,i_shell) = one
             call find_eigenstate_and_energy_vkb(i_species,en,ell,paos%cutoff(zeta,i_shell),&
                  paos%psi(zeta,i_shell)%f,paos%energy(zeta,i_shell), &
                  vha,vxc,paos%width(i_shell),paos%prefac(i_shell))
             ! Orthogonalise to semi-core state
             if(paos%inner(i_shell)>0) then
                ! Dot product of two
                dot_p = zero
                do i=1,nmesh
                   dot_p = dot_p + rr(i)**(2*ell+2)*paos%psi(zeta,i_shell)%f(i)* &
                        paos%psi(1,paos%inner(i_shell))%f(i)*drdi(i)
                end do
                if(iprint>2) write(*,fmt='(2x,"Orthogonalising to semi-core; overlap is ",f10.5)') dot_p
                ! Orthogonalise
                paos%psi(zeta,i_shell)%f = paos%psi(zeta,i_shell)%f - &
                     dot_p * paos%psi(1,paos%inner(i_shell))%f
                ! Normalise
                dot_p = zero
                do i=1,nmesh
                   dot_p = dot_p + rr(i)**(2*ell+2)*paos%psi(zeta,i_shell)%f(i)* &
                        paos%psi(zeta,i_shell)%f(i)*drdi(i)
                end do
                if(iprint>2) write(*,fmt='(2x,"Normalising: ",f10.5)') dot_p
                paos%psi(zeta,i_shell)%f = paos%psi(zeta,i_shell)%f/sqrt(dot_p)
             end if
          end if
       end do
    end do
    !
    ! If perturbative polarisation is selected, it is the last shell
    !
    if(paos%flag_perturb_polarise) then       
       i_shell = paos%n_shells
       ! NB we pass l-1 here as the angular momentum of the shell being perturbed
       ell = paos%l(i_shell)-1
       ! Note that npao is set to the value of n for the shell being perturbed
       en = paos%npao(i_shell)
       do zeta = 1,paos%nzeta(i_shell)
          allocate(paos%psi(zeta,i_shell)%f(nmesh))
          paos%psi(zeta,i_shell)%f = zero
          if(iprint>2) then
             write(*,fmt='(2x,"Perturbative polarisation")')
             write(*,fmt='(2x,"Species ",i2," n=",i2," l=",i2," zeta=",i2, " Rc=",f4.1," bohr")') &
                  i_species, en, ell, zeta, paos%cutoff(zeta,i_shell)
             write(*,fmt='(4x,"Prefactor scaled by ",f5.1)') paos%pol_pf
             write(*,fmt='(4x,"Energy: ",f12.5)') paos%energy(zeta,paos%polarised_shell)
             write(*,fmt='(4x,"Polarising shell: ",i2)') paos%polarised_shell
          end if
          if(zeta>1.AND.paos%flag_zetas==1) then
             call find_split_norm(en,ell+1,paos%cutoff(zeta,i_shell),&
                  paos%psi(zeta,i_shell)%f,paos%psi(1,i_shell)%f)
          else
             ! Set radius of polarisation PAO
             paos%cutoff(zeta,i_shell) = paos%cutoff(zeta,paos%polarised_shell)
             call find_polarisation(i_species,en,ell,paos%cutoff(zeta,i_shell),&
                 paos%psi(zeta,paos%polarised_shell)%f,paos%psi(zeta,i_shell)%f,&
                 paos%energy(zeta,paos%polarised_shell),vha,vxc,paos%pol_pf)
             ! Orthogonalise to semi-core state
             if(paos%inner(i_shell)>0) then
                ! Dot product of two
                ! NB r^(2l+4) is because each psi needs to be scaled by r^(l+1)
                ! to remove the Siesta normalisation, and then the integral needs
                ! another r^2.  Also below when normalising.
                dot_p = zero
                do i=1,nmesh
                   dot_p = dot_p + rr(i)**(2*ell+4)*paos%psi(zeta,i_shell)%f(i)* &
                        paos%psi(1,paos%inner(i_shell))%f(i)*drdi(i)
                end do
                if(iprint>2) write(*,fmt='(2x,"Orthogonalising to semi-core; overlap is ",f10.5)') dot_p
                ! Orthgonalise
                paos%psi(zeta,i_shell)%f = paos%psi(zeta,i_shell)%f - &
                     dot_p * paos%psi(1,paos%inner(i_shell))%f
                ! Normalise
                dot_p = zero
                do i=1,nmesh
                   dot_p = dot_p + rr(i)**(2*ell+4)*paos%psi(zeta,i_shell)%f(i)* &
                        paos%psi(zeta,i_shell)%f(i)*drdi(i)
                end do
                if(iprint>2) write(*,fmt='(2x,"Normalising: ",f10.5)') dot_p
                paos%psi(zeta,i_shell)%f = paos%psi(zeta,i_shell)%f/sqrt(dot_p)
             end if
          end if
       end do
    end if ! paos%flag_perturb_polarise
    return
  end subroutine solve_for_polarisation

  ! Interpolate PAOs onto regular grid
  subroutine interpolate_paos(i_species)

    use datatypes
    use numbers
    use mesh, ONLY: nmesh, rr, delta_r_reg, interpolate, make_mesh_reg, convert_r_to_i
    use pseudo_tm_info, ONLY: pseudo
    use write, ONLY: write_pao_plot

    implicit none

    ! Passed variables
    integer :: i_species

    ! Local variables
    integer :: i_shell, en, ell, zeta, nmesh_pao, nrc
 
    do i_shell = 1,paos%n_shells
       ell = paos%l(i_shell)
       en = paos%n(i_shell)
       do zeta = 1,paos%nzeta(i_shell)
          call convert_r_to_i(paos%cutoff(zeta,i_shell),nrc)
          paos%cutoff(zeta,i_shell) = rr(nrc-1)
          nmesh_pao = floor(paos%cutoff(zeta,i_shell)/delta_r_reg) + 1
          allocate(paos%psi_reg(zeta,i_shell)%f(nmesh_pao), paos%psi_reg(zeta,i_shell)%x(nmesh_pao))
          call make_mesh_reg(paos%psi_reg(zeta,i_shell)%x,nmesh_pao,paos%cutoff(zeta,i_shell))
          paos%psi_reg(zeta,i_shell)%delta = paos%cutoff(zeta,i_shell)/real(nmesh_pao-1,double)
          paos%psi_reg(zeta,i_shell)%n = nmesh_pao
          ! Interpolate
          call interpolate(paos%psi_reg(zeta,i_shell)%x,paos%psi_reg(zeta,i_shell)%f,nmesh_pao,&
               rr(1:nrc-1),paos%psi(zeta,i_shell)%f(1:nrc-1),nrc-1,zero)
          if(flag_plot_output) call write_pao_plot(nint(pseudo(i_species)%z),paos%psi_reg(zeta,i_shell)%x, &
               paos%psi_reg(zeta,i_shell)%f, nmesh_pao,"PAO", en,ell,zeta)
       end do
    end do
  end subroutine interpolate_paos
  
  ! Interpolate potentials and charges onto regular grid
  subroutine interpolate_potentials(i_species,vha_conf)

    use datatypes
    use numbers
    use mesh, ONLY: nmesh, rr, delta_r_reg, interpolate, make_mesh_reg, convert_r_to_i
    use pseudo_tm_info, ONLY: pseudo, rad_alloc
    use write, ONLY: write_pao_plot
    use read, ONLY: ps_format, oncvpsp

    implicit none

    ! Passed variables
    integer :: i_species
    real(double), dimension(nmesh) :: vha_conf

    ! Local variables
    integer :: i_shell, en, ell, zeta, nmesh_pot, i, j, k, nrc, istart
    real(double) :: max_cutoff
    real(double) :: zz
    real(double), allocatable, dimension(:) :: x_reg

    !
    ! Kleinman-Bylander projectors
    !
    j = 0
    if(iprint>1) write(*,fmt='(/4x,"VKB projectors")')
    do ell = 0, pseudo(i_species)%lmax
       do i=1,local_and_vkb%n_proj(ell)
          j = j+1
          ! Scale projector by r**(l+1)
          !if(ps_format==oncvpsp) then
          do k=0,ell
             local_and_vkb%projector(:,i,ell) = local_and_vkb%projector(:,i,ell)/rr
          end do
          !end if
          ! Find actual cutoff: two successive points with magnitude less than RD_ERR
          ! We may want to start this somewhere r=0.1 to avoid errors
          nrc = local_and_vkb%ngrid_vkb
          call convert_r_to_i(0.1_double,istart)
          do k = istart, local_and_vkb%ngrid_vkb-1
             if(abs(local_and_vkb%projector(k,i,ell))<RD_ERR.AND.&
                  abs(local_and_vkb%projector(k+1,i,ell))<RD_ERR) then
                nrc = k
                exit
             end if
          end do
          ! Create regular mesh
          pseudo(i_species)%pjnl(j)%cutoff = rr(nrc)
          nmesh_pot = floor(rr(nrc)/delta_r_reg) + 1 ! +2 in original: why?!
          allocate(x_reg(nmesh_pot))
          call make_mesh_reg(x_reg,nmesh_pot,rr(nrc))
          call rad_alloc(pseudo(i_species)%pjnl(j),nmesh_pot)
          pseudo(i_species)%pjnl(j)%delta = rr(nrc)/real(nmesh_pot-1,double)
          call interpolate(x_reg,pseudo(i_species)%pjnl(j)%f,nmesh_pot, &
               rr(1:nrc),local_and_vkb%projector(1:nrc,i,ell),nrc,zero)
          if(flag_plot_output) call write_pao_plot(nint(pseudo(i_species)%z),x_reg, &
               pseudo(i_species)%pjnl(j)%f,nmesh_pot,"KB",ell=ell,zeta=i)
          deallocate(x_reg)
       end do
    end do
    !
    ! Partial core charge
    !
    if(pseudo(i_species)%flag_pcc) then
       if(iprint>1) write(*,fmt='(/4x,"Partial core correction")')
       ! Find core charge cutoff
       do i=1,nmesh
          ! Find cutoff radius: two successive points less than 1e-8 
          if(local_and_vkb%pcc(i)<RD_ERR.AND.i<nmesh) then
             if(local_and_vkb%pcc(i+1)<RD_ERR) then
                pseudo(i_species)%chpcc%cutoff = rr(i)
                exit
             end if
          end if
       end do
       nmesh_pot = floor(pseudo(i_species)%chpcc%cutoff/delta_r_reg) + 1 ! 2 in original?
       allocate(x_reg(nmesh_pot))
       call rad_alloc(pseudo(i_species)%chpcc,nmesh_pot)
       pseudo(i_species)%chpcc%delta = pseudo(i_species)%chpcc%cutoff/real(nmesh_pot-1,double)
       call make_mesh_reg(x_reg,nmesh_pot,pseudo(i_species)%chpcc%cutoff)
       call interpolate(x_reg,pseudo(i_species)%chpcc%f,nmesh_pot, &
            rr,local_and_vkb%pcc,local_and_vkb%ngrid,zero)
       if(flag_plot_output) call write_pao_plot(nint(pseudo(i_species)%z),x_reg, &
            pseudo(i_species)%chpcc%f,nmesh_pot,"PCC")
       deallocate(x_reg)
    end if
    !
    ! Local - the radius is taken as largest PAO for compatibility with NA
    !
    if(iprint>1) write(*,fmt='(/4x,"Local potential")')
    max_cutoff = maxval(paos%cutoff)
    nmesh_pot = max_cutoff/delta_r_reg + 1
    allocate(x_reg(nmesh_pot))
    call rad_alloc(pseudo(i_species)%vlocal,nmesh_pot)
    pseudo(i_species)%vlocal%cutoff = max_cutoff
    pseudo(i_species)%vlocal%delta = max_cutoff/real(nmesh_pot-1,double)
    call make_mesh_reg(x_reg,nmesh_pot,pseudo(i_species)%vlocal%cutoff)
    !
    ! Interpolate using exact value of local potential at cutoff: -Z/r
    !
    zz = pseudo(i_species)%z - pseudo(i_species)%zcore
    call interpolate(x_reg,pseudo(i_species)%vlocal%f,nmesh_pot, &
         !rr,local_and_vkb%local,local_and_vkb%ngrid,-pseudo(i_species)%zval/pseudo(i_species)%vlocal%cutoff)
         rr,local_and_vkb%local,local_and_vkb%ngrid,-zz/pseudo(i_species)%vlocal%cutoff)
    if(flag_plot_output) call write_pao_plot(nint(pseudo(i_species)%z),x_reg, &
         pseudo(i_species)%vlocal%f,nmesh_pot,"Vlocal")
    !
    ! Neutral atom potential - same mesh as local
    !
    if(iprint>1) write(*,fmt='(/4x,"Neutral atom potential")')
    call rad_alloc( pseudo(i_species)%vna, nmesh_pot )
    pseudo(i_species)%vna%delta = pseudo(i_species)%vlocal%delta
    pseudo(i_species)%vna%cutoff = pseudo(i_species)%vlocal%cutoff
    call interpolate(x_reg,pseudo(i_species)%vna%f,pseudo(i_species)%vna%n,&
         rr,vha_conf,nmesh,zero)
    if(flag_plot_output) call write_pao_plot(nint(pseudo(i_species)%z),x_reg, &
         pseudo(i_species)%vna%f,nmesh_pot,"VNA")
    deallocate(x_reg)
    return
  end subroutine interpolate_potentials
  
  ! For the default basis, find the cutoff radii for all shells (based on energies)
  subroutine find_default_cutoffs(i_species,vha,vxc)

    use datatypes
    use numbers
    use mesh, ONLY: nmesh, rr, delta_r_reg, convert_r_to_i
    use units, ONLY: HaToeV
    
    implicit none

    ! Passed variables
    integer :: i_species
    real(double), dimension(nmesh) :: vha,vxc

    ! Local variables
    integer :: i_shell, grid_size, n_shells_average
    real(double), dimension(:), allocatable :: large_cutoff, small_cutoff
    real(double) :: energy, average_large, average_small, average_mid
    real(double), dimension(:), allocatable :: psi

    ! Find the cutoffs for the valence shells first
    allocate(large_cutoff(val%n_occ),small_cutoff(val%n_occ))
    large_cutoff = zero
    small_cutoff = zero
    write(*,fmt='(/4x,"Default energy shifts")')
    write(*,fmt='(4x,"  n  l   delta E (Ha) delta E (eV)")')
    ! Loop over valence states, find large/small cutoffs
    do i_shell = 1, val%n_occ !paos%n_shells-1 
       if(iprint>3) write(*,*) '# Finding radius for ',paos%npao(i_shell), paos%l(i_shell), &
            val%en_ps(i_shell)+deltaE_large_radius
       call find_radius_from_energy(i_species,paos%npao(i_shell), paos%l(i_shell), &
            large_cutoff(i_shell), val%en_ps(i_shell)+deltaE_large_radius, vha, vxc, .false.)
       ! Round to grid step
       if(val%semicore(i_shell)==0) then
          write(*,fmt='(4x,2i3,2f13.8," (large radius)")') paos%n(i_shell), paos%l(i_shell), &
               deltaE_large_radius, deltaE_large_radius*HaToeV
          write(*,fmt='(4x,2i3,2f13.8," (small radius)")') paos%n(i_shell), paos%l(i_shell), &
               deltaE_small_radius, deltaE_small_radius*HaToeV
          call find_radius_from_energy(i_species,paos%npao(i_shell), paos%l(i_shell), &
               small_cutoff(i_shell), val%en_ps(i_shell)+deltaE_small_radius, vha, vxc, .false.)
       else
          write(*,fmt='(4x,2i3,2f13.8," (only radius)")') paos%n(i_shell), paos%l(i_shell), &
               deltaE_large_radius, deltaE_large_radius*HaToeV
          small_cutoff(i_shell) = large_cutoff(i_shell)
       end if
       write(*,*) '# Radii: ',large_cutoff(i_shell),small_cutoff(i_shell)
    end do
    ! Create cutoffs based on defaults chosen by user
    if(paos%flag_cutoff==pao_cutoff_energies.OR.paos%flag_cutoff==pao_cutoff_default) then ! Same energy for all l/n shells
       do i_shell = 1, val%n_occ !paos%n_shells-1
          ! NB Semi-core orbitals have large = small
          paos%cutoff(:,i_shell) = zero
          paos%cutoff(1,i_shell) = large_cutoff(i_shell)
          if(paos%nzeta(i_shell)==2) then
             paos%cutoff(2,i_shell) = small_cutoff(i_shell)
          else if(paos%nzeta(i_shell)==3) then
             paos%cutoff(2,i_shell) = half*(large_cutoff(i_shell) + small_cutoff(i_shell))
             ! Locate nearest logarithmic grid point
             call convert_r_to_i(paos%cutoff(2,i_shell),grid_size)
             paos%cutoff(2,i_shell) = rr(grid_size-1)
             paos%cutoff(3,i_shell) = small_cutoff(i_shell)
          end if
       end do
       ! Set polarisation radii
       if(paos%n_shells>val%n_occ) then ! Polarisation
          paos%cutoff(1,paos%n_shells)=paos%cutoff(1,paos%n_shells-1)
          if(paos%nzeta(i_shell)==2) then
             paos%cutoff(2,paos%n_shells)=paos%cutoff(2,paos%n_shells-1)
          else if(paos%nzeta(i_shell)==3) then
             paos%cutoff(2,paos%n_shells)=paos%cutoff(2,paos%n_shells-1)
             paos%cutoff(3,paos%n_shells)=paos%cutoff(3,paos%n_shells-1)
          end if
       end if
    else if(paos%flag_cutoff==pao_cutoff_radii) then ! Same radius for all l/n shells
       ! Sum over radii of different shells, excluding polarisation
       average_large = zero
       average_small = zero
       n_shells_average = 0
       do i_shell = 1, val%n_occ!paos%n_shells-1
          if(val%semicore(i_shell)/=1) then
             average_large = average_large + large_cutoff(i_shell)
             average_small = average_small + small_cutoff(i_shell)
             n_shells_average = n_shells_average + 1
          end if
       end do
       ! Find averages for large and small
       average_large = average_large/real(n_shells_average, double)
       average_small = average_small/real(n_shells_average, double)
       ! Round to grid step
       call convert_r_to_i(average_large,grid_size)
       average_large = rr(grid_size-1)
       call convert_r_to_i(average_small,grid_size)
       average_small = rr(grid_size-1)
       call convert_r_to_i(half*(average_large+average_small),grid_size)
       average_mid = rr(grid_size-1)
       ! Set radii for all shells
       do i_shell = 1, val%n_occ
          if(val%semicore(i_shell)==1) then
             paos%cutoff(:,i_shell) = zero
             paos%cutoff(1,i_shell) = large_cutoff(i_shell)
          else
             paos%cutoff(:,i_shell) = zero
             paos%cutoff(1,i_shell) = average_large
             if(paos%nzeta(i_shell)==2) then
                paos%cutoff(2,i_shell) = average_small
             else if(paos%nzeta(i_shell)==3) then
                paos%cutoff(3,i_shell) = average_small
                ! Now average large and small
                paos%cutoff(2,i_shell) = average_mid !half*(average_large + average_small)
             end if
          end if
       end do
       do i_shell = val%n_occ + 1, paos%n_shells
          paos%cutoff(:,i_shell) = zero
          paos%cutoff(1,i_shell) = average_large
          if(paos%nzeta(i_shell)==2) then
             paos%cutoff(2,i_shell) = average_small
          else if(paos%nzeta(i_shell)==3) then
             paos%cutoff(3,i_shell) = average_small
             ! Now average large and small
             paos%cutoff(2,i_shell) = average_mid !half*(average_large + average_small)
          end if
       end do
    end if
    deallocate(large_cutoff,small_cutoff)
    return
  end subroutine find_default_cutoffs
  
  ! Given an energy shift, integrate outwards to find the corresponding radius
  subroutine find_radius_from_energy(i_species,en,ell,Rc,energy,vha,vxc,flag_use_semilocal)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: rr, rr_squared, nmesh, alpha, make_mesh, beta, drdi, drdi_squared
    use pseudo_tm_info, ONLY: pseudo
    
    implicit none

    integer :: i_species, ell, en
    real(double) :: energy, Rc
    real(double), dimension(nmesh) :: vha,vxc

    ! Local variables
    real(double) :: g_temp, dy_L, dy_R
    real(double), dimension(:), allocatable :: f, potential, psi
    integer :: classical_tp, i, n_crossings, n_nodes, n_loop, loop, nmax
    integer :: j, n_kink, n_nodes_lower, n_nodes_upper
    real(double) :: fac, norm, d_energy, r_lower, r_upper, df_cusp, cusp_psi, tol
    real(double) :: l_l_plus_one, alpha_sq_over_four, xkap
    real(double) :: gin, gout, gsgin, gsgout, xin, xout
    logical :: flag_use_semilocal

    if(iprint>5) write(*,*) '# Finding radius for energy: ',energy, en, ell
    l_l_plus_one = real(ell*(ell+1),double)
    alpha_sq_over_four = alpha*alpha/four
    n_nodes = en - ell - 1 
    allocate(f(nmesh),potential(nmesh),psi(nmesh))
    f = zero
    potential = zero
    psi = zero
    nmax = nmesh
    n_crossings = 0
    if(flag_use_semilocal) then
       write(*,*) '# Using semi-local potential'
       do i=1,nmesh
          potential(i) = local_and_vkb%semilocal_potential(i,ell) + vha(i) + vxc(i)
          g_temp = (drdi_squared(i)*(two*(energy - potential(i))-l_l_plus_one/rr_squared(i)) - alpha_sq_over_four)/twelve
          f(i) = one + g_temp
       end do
       psi(1) = one
       psi(2) = ( (twelve - ten*f(1)) * psi(1) )/f(2) ! Implicit psi(0) = zero
       dy_R = psi(2) - psi(1)
       do i=2,nmax-1
          dy_L = dy_R
          psi(i+1) = ( (twelve - ten*f(i)) * psi(i) - f(i-1)*psi(i-1) )/f(i+1)
          dy_R = psi(i+1) - psi(i)
          if(psi(i)*psi(i+1)<zero.OR.(abs(psi(i+1))<1e-16_double)) then ! Crossing the x-axis
             !write(*,*) '# Node ! ',rr(i),psi(i),psi(i+1)
             n_crossings = n_crossings + 1
             if(n_crossings>=n_nodes+1) exit
          end if
       end do
    else
       !write(*,*) '# Using VKB potentials'
       do i=1,nmesh
          potential(i) = local_and_vkb%local(i) + vha(i) + vxc(i)
          g_temp = (drdi_squared(i)*(two*(energy - potential(i))-l_l_plus_one/rr_squared(i)) - alpha_sq_over_four)/twelve
          f(i) = one + g_temp
       end do
       n_kink = nmax
       call integrate_vkb_outwards(i_species,n_kink,ell,psi,f,n_crossings,n_nodes) ! We want to integrate to max before final node
       if(n_crossings>=n_nodes+1) then
          write(*,fmt='(2x,"Found cutoff inside the VKB projector cutoff: this suggests too great an energy shift")')
          write(*,fmt='(2x,"Try reducing Atom.dE_small_radius from its default of 2 eV")')
          write(*,fmt='(2x,"Alternatively this shell should be semi-core - adjust General.SemicoreEnergy")')
          call cq_abort("Aborting")
       else
          do i=n_kink,nmax-1
             psi(i+1) = ( (twelve - ten*f(i)) * psi(i) - f(i-1)*psi(i-1) )/f(i+1)
             if(psi(i)*psi(i+1)<zero.OR.(abs(psi(i+1))<1e-16_double)) then ! Crossing the x-axis
                !write(*,*) '# Node ! ',rr(i),psi(i),psi(i+1)
                n_crossings = n_crossings + 1
                if(n_crossings>=n_nodes+1) exit
             end if
          end do
       end if
       !write(*,*) '# Crossings, nodes: ',n_crossings, n_nodes
       if(n_crossings<n_nodes+1) call cq_abort("Failed to find confined state - please check input")
    end if
    ! Find radius by integrating outwards
    ! Interpolate between i+1 and i
    if(i>=nmax) i=nmax-1
    Rc = rr(i) - psi(i)*(rr(i+1)-rr(i))/(psi(i+1)-psi(i))
    !write(*,*) 'ri, ri+1 and interp are: ',rr(i), rr(i+1),psi(i),psi(i+1), - psi(i)*(rr(i+1)-rr(i))/(psi(i+1)-psi(i))
    if(iprint>5) write(*,*) '# Found radius ',Rc
    deallocate(f,potential,psi)
    return
  end subroutine find_radius_from_energy
    
  ! Rc has been specified - find the energy that gives this cutoff
  ! VKB version
  ! Use local potential for homogeneous equation, and solve inhomogeneous
  ! equations for projectors (then combine solutions)
  subroutine find_eigenstate_and_energy_vkb(i_species,en,ell,Rc,psi,energy,vha,vxc,width,prefac)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: rr, rr_squared, nmesh, alpha, make_mesh, beta, drdi, drdi_squared, convert_r_to_i
    use pseudo_tm_info, ONLY: pseudo
    use read, ONLY: e_step, max_solver_iters
    
    implicit none

    integer :: i_species, ell, en
    real(double) :: energy, Rc
    real(double), dimension(nmesh) :: psi,vha,vxc
    real(double), OPTIONAL :: width, prefac

    ! Local variables
    real(double) :: g_temp, dy_L, dy_R
    real(double), dimension(:), allocatable :: f, potential
    integer :: classical_tp, i, n_crossings, n_nodes, n_loop, loop, nmax, n_kink, n_nodes_lower, n_nodes_upper, n_kink_vkb
    real(double) :: l_half_sq, dx_sq_over_twelve, fac, norm, d_energy, e_lower, e_upper, df_cusp, cusp_psi, tol
    real(double) :: delta_energy_bracket, zval, l_l_plus_one, alpha_sq_over_four, xkap
    real(double) :: gin, gout, gsgin, gsgout, xin, xout, exponent, delta
    logical :: flag_find_radius = .false.
    logical :: flag_confine

    if(iprint>2) write(*,fmt='(2x,"Entering find_eigenstate_and_energy_vkb")')
    flag_confine = .false.
    if(present(width).AND.present(prefac)) then
       if(abs(prefac)>RD_ERR) flag_confine = .true.
    end if
    l_half_sq = real(ell,double) + half
    l_half_sq = l_half_sq*l_half_sq
    l_l_plus_one = real(ell*(ell+1),double)
    zval = pseudo(i_species)%zval ! Added
    !zval = pseudo(i_species)%z  (z or zval?)
    !write(*,*) '# zval is ',zval
    dx_sq_over_twelve = alpha*alpha/twelve
    alpha_sq_over_four = alpha*alpha/four
    n_nodes = en - ell - 1 
    allocate(f(nmesh),potential(nmesh))
    if(abs(energy)<RD_ERR) then
       e_lower = -two*zval!-zval*zval/real(en*en,double)
    else if(energy<zero) then
       e_lower = energy*1.2_double
    else ! Unbound (polarisation) state
       e_lower = -half!zero
    end if
    ! Energy bounds - allow for unbound states
    e_upper = half!five ! One failed to find the tightest PAO for O
    do i=1,nmesh
       potential(i) = local_and_vkb%local(i) + vha(i) + vxc(i)  ! Half when using Siesta
       !g_temp = l_l_plus_one/(rr_squared(i)) + potential(i)
       !if(g_temp<e_lower) e_lower = g_temp
    end do
    ! Now set number of loops and maximum radius
    n_loop = max_solver_iters
    call convert_r_to_i(Rc,nmax)
    nmax = nmax - 1
    ! NEW !
    Rc = rr(nmax)
    ! NEW !
    if(abs(energy)<RD_ERR) then
       energy = zero!half*(e_lower+e_upper)
    end if
    tol = 1.0e-6_double
    if(flag_confine) then
       delta = 0.01_double
       do i=1,nmax
          ! Locally rc = Rc  + delta, ri = Rc - width
          if(rr(i)>(Rc + delta - width)) then
             exponent = (width + delta)/(rr(i) - Rc + width)
             potential(i) = potential(i) + prefac*exp(-exponent)/(Rc + delta - rr(i))
          end if
       end do
    end if
    if(iprint>2) write(*,fmt='(2x,"For this eigenstate, we require ",i2," nodes")') n_nodes
    ! Now loop to find the correct energy
    do loop = 1,n_loop
       if(iprint>4) write(*,fmt='(2x,"Loop ",i3," with energy brackets of ",3f18.10)') loop, e_lower,energy,e_upper
       n_kink = nmax
       do i=1,nmax
          g_temp = (drdi_squared(i)*(two*(energy - potential(i))-l_l_plus_one/rr_squared(i)) - alpha_sq_over_four)/twelve
          f(i) = one + g_temp
       end do
       psi = zero
       ! Outward
       ! We will need homogeneous and n_proj inhomogeneous
       n_crossings = 0
       call integrate_vkb_outwards(i_species,n_kink,ell,psi,f,n_crossings,n_nodes) ! We want to integrate to max before final node
       n_kink_vkb = n_kink
       if(iprint>4) write(*,fmt='(2x,"Kink is at ",f18.10," with ",i2," crossings")') rr(n_kink),n_crossings
       ! If we haven't found enough nodes, we need to try further
       if(n_crossings/=n_nodes) then
          !write(*,*) 'Found ',n_crossings,' crossings so far; continuing ',n_kink_vkb
          n_kink = nmax-n_kink_vkb+1
          call numerov(n_kink_vkb-1,n_kink,nmax,psi,rr,f,1,n_crossings,n_nodes,xkap,0)
          !write(*,*) 'Left numerov with kink, crossings: ',n_kink,n_crossings
       end if
       if(iprint>4) write(*,fmt='(2x,"Kink is at ",f18.10)') rr(n_kink)
       xout = psi(n_kink)*f(n_kink) - psi(n_kink-1)*f(n_kink-1)
       gout = psi(n_kink)
       gsgout = psi(n_kink)*psi(n_kink)*drdi_squared(n_kink)
       if(n_kink == nmax) then
          if(iprint>4) write(*,fmt='(2x,"No kink found - adjusting lower bound")')
          e_lower = energy
          if(e_upper-e_lower>e_step) then
             energy = energy + e_step
          else
             energy = half*(e_lower+e_upper)
          end if
          cycle
       end if
       if(n_crossings /= n_nodes) then
          if ( n_crossings > n_nodes ) then
             e_upper = energy
             if(e_upper-e_lower>e_step) then
                energy = energy - e_step
             else
                energy = half * ( e_upper + e_lower )
             end if
          else
             e_lower = energy
             if(e_upper-e_lower>e_step) then
                energy = energy + e_step
             else
                energy = half * ( e_upper + e_lower )
             end if
          end if
          !energy = half * ( e_upper + e_lower )
          if(iprint>4) write(*,fmt='(2x,"Nodes found: ",i3," but required: ",i3)') n_crossings, n_nodes
          cycle
       end if       
       ! Used for matching
       fac = psi(n_kink)
       ! Inward
       i = nmax
       xkap = sqrt(two*(potential(i) - energy)+l_l_plus_one/rr_squared(i))
       psi(nmax) = zero!exp(-xkap)
       psi(nmax-1) = one !exp(-xkap*(rr(nmax) - rr(nmax-1)))
       call numerov(n_kink+1,nmax-n_kink,nmax,psi,rr,f,-1,n_crossings,n_nodes,xkap,0)
       if(iprint>5) write(*,fmt='(2x,"Values of psi to L and R of kink point",2f18.10)') fac,psi(n_kink)
       fac = fac/psi(n_kink)
       psi(n_kink:nmax) = psi(n_kink:nmax)*fac
       xin = psi(n_kink)*f(n_kink)- psi(n_kink+1)*f(n_kink+1)
       gin = psi(n_kink)
       gsgin = psi(n_kink)*psi(n_kink)*drdi_squared(n_kink)       
       ! Remember that psi is y in numerov - don't forget factor of root(r)
       ! Normalise
       norm = zero
       !call integrate_simpson(psi*psi,drdi_squared,nmax,norm)
       do i=1,nmax
          norm = norm + psi(i)*psi(i)*drdi_squared(i)
       end do
       psi = psi/sqrt(norm)
       ! It is possible to improve the search by predicting the necessary shift in energy
       ! The lines below feature various ways of doing this but at the moment they aren't
       ! needed as the search is very fast.  Kept as useful diagnostics.
       ! All the lines below are redundant, I think
       !%%! ! [yL(i-1) + yR(i+1) - (14-12f(i))y(i)]/dx > 0 means energy too high
       !%%! dy_L = (psi(n_kink-1)*sqrt(drdi(n_kink-1))/rr(n_kink-1)-psi(n_kink-2)*sqrt(drdi(n_kink-2))/rr(n_kink-2)) &
       !%%!      /(rr(n_kink-1)-rr(n_kink-2))
       !%%! dy_R = (psi(n_kink+2)*sqrt(drdi(n_kink+2))/rr(n_kink+2)-psi(n_kink+1)*sqrt(drdi(n_kink+1))/rr(n_kink+1)) &
       !%%!      /(rr(n_kink+2)-rr(n_kink+1))
       !%%! if(iprint>5) then
       !%%!    write(*,fmt='("Gradient of psi on L: ",f18.10)') dy_L
       !%%!    write(*,fmt='("Gradient of psi on L (alt): ",f18.10)') (psi(n_kink-1)-psi(n_kink-2))/(rr(n_kink-1)-rr(n_kink-2))
       !%%!    write(*,fmt='("Gradient of psi on R: ",f18.10)') dy_R
       !%%!    write(*,fmt='("Gradient of psi on R (alt): ",f18.10)') (psi(n_kink+2)-psi(n_kink+1))/(rr(n_kink+2)-rr(n_kink+1))
       !%%! end if
       !%%! d_energy = dy_R - dy_L
       !%%! if(iprint>4) write(*,fmt='("Difference in gradients: ",f18.10)') d_energy
       i = n_kink
       ! Perturbation theory to improve energy
       cusp_psi = ( psi(i-1)*f(i-1) + f(i+1)*psi(i+1) + ten*f(i)*psi(i)) / twelve
       df_cusp = f(i)*( psi(i)/cusp_psi - one )
       d_energy = (df_cusp*twelve * cusp_psi * cusp_psi)/drdi(i)
       if(iprint>5) write(*,fmt='(2x,"Energies (low, mid, high): ",3f18.10)') e_lower,energy,e_upper
       if(iprint>4) write(*,fmt='(2x,"Energy shift      : ",f18.10)') d_energy
       ! Alternate approach
       !cusp_psi = xin + xout
       !cusp_psi = cusp_psi + twelve*(one-f(n_kink))*gout
       !cusp_psi = cusp_psi*gout/(gsgin + gsgout)
       !d_energy = cusp_psi
       if(iprint>4) write(*,fmt='(2x,"Energy shift (alt): ",f18.10)') cusp_psi
       if(iprint>5) write(*,fmt='(2x,"Number of nodes: ",i4)') n_crossings
       if ( n_crossings /= n_nodes) then
          if ( n_crossings > n_nodes ) then
             e_upper = energy
             if(e_upper-e_lower>e_step) then
                energy = energy - e_step
             else
                energy = half * ( e_upper + e_lower )
             end if
          else
             e_lower = energy
             if(e_upper-e_lower>e_step) then
                energy = energy + e_step
             else
                energy = half * ( e_upper + e_lower )
             end if
          end if
          !energy = half * ( e_upper + e_lower )
          cycle
       end if       
       if(d_energy>zero) then
          e_lower = energy
       else
          e_upper = energy
       end if
       if(energy+d_energy<e_upper.AND.energy+d_energy>e_lower) then
          if(abs(d_energy)<e_step) then
             energy = energy + d_energy
          else
             if(d_energy>zero) then
                energy = energy + e_step
             else
                energy = energy - e_step
             end if
          end if
       else
          energy = half*(e_lower + e_upper)
       end if
       if(abs(d_energy)<tol) exit
    end do
    if(loop>=n_loop.AND.abs(d_energy)>tol) call cq_abort("Error: failed to find energy for n,l: ",en,ell)
    if(iprint>2) write(*,fmt='(2x,"Final energy found: ",f11.6," Ha")') energy
    ! Rescale - remove factor of sqrt r
    do i=1,nmax
       psi(i) = psi(i)*sqrt(drdi(i))/rr(i)
    end do
    ! Adjust so that psi approaches zero from above for large r
    if(psi(nmax - 5)<zero) psi = -psi
    if(ell>0) then
       do loop = 1,ell
          psi = psi/rr
       end do
    end if
    deallocate(f,potential)
  end subroutine find_eigenstate_and_energy_vkb

  ! Fit a simple polynomial (a + b*r*r)*r**(l+1) up to Rc and then match psi_match
  ! from that point on.  We then subtract this function off psi_match to generate a
  ! more compressed orbital, but one whose gradient goes to zero at Rc.
  !
  ! Follows the ideas of Siesta (e.g. JPCM 14, 2745 2002)
  subroutine find_split_norm(en,ell,Rc,psi,psi_match)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: rr, rr_squared, nmesh, alpha, make_mesh, beta, drdi, drdi_squared, convert_r_to_i
    use pseudo_tm_info, ONLY: pseudo
    
    implicit none

    ! Passed variables
    integer :: ell, en
    real(double) :: Rc
    real(double), dimension(nmesh) :: psi,psi_match

    ! Local variables
    integer :: nrc, i
    real(double) :: psi_val, dpsi_val, a, b, c, d, rmatch, norm, rl, dpsip, dpsim, d2psi
    
    ! Find matching point
    call convert_r_to_i(Rc,nrc)
    nrc = nrc-1
    rmatch = rr(nrc)
    ! Values of psi and dpsi at Rc
    psi_val = rr(nrc)**ell*psi_match(nrc)
    dpsi_val = half*(rr(nrc+1)**ell*psi_match(nrc+1) - rr(nrc-1)**ell*psi_match(nrc-1)) / drdi(nrc)
    ! These are used for higher-order polynomial matching
    !dpsip = half*(rr(nrc+2)**ell*psi_match(nrc+2) - rr(nrc)**ell*psi_match(nrc)) / drdi(nrc+1)
    !dpsim = half*(rr(nrc)**ell*psi_match(nrc) - rr(nrc-2)**ell*psi_match(nrc-2)) / drdi(nrc-1)
    !d2psi = half*(dpsip - dpsim) / drdi(nrc)
    ! Coefficients of polynomial
    rl = one
    do i=1,ell
       rl = rl*rmatch
    end do
    ! Original, simple quadratic a - b*r*r, matching psi and dpsi at rc
    b = half*(dpsi_val*rmatch - real(ell,double)*psi_val)/(rl*rmatch*rmatch)
    a = psi_val/rl - b*rmatch*rmatch
    if(iprint>3) write(*,*) "Split norm coefficients, radius: ",a,b,rmatch
    ! -------------------------------------------------------------
    ! General quadratic a + b*r + c*r*r (added matching psi at r=0)
    ! -------------------------------------------------------------
    !if(ell>0) then
    !   a = zero
    !else
    !   a = psi_match(1)
    !end if
    !c = (dpsi_val*rmatch + rl*a - real(ell+1,double)*psi_val)/(rl*rmatch*rmatch)
    !b = (dpsi_val - real(ell,double) * psi_val/rmatch - two*c*rl*rmatch)/rl
    !if(iprint>3) write(*,*) "Split norm coefficients, radius: ",a,b,c,rmatch
    ! -------------------------------------------------------------
    ! Alternate general quadratic a + b*r + c*r*r (added matching d2psi at rc)
    ! -------------------------------------------------------------
    !c = half*(d2psi*rmatch*rmatch - real(ell,double) * (dpsi_val * rmatch - psi_val) )/(rl*rmatch*rmatch)
    !b = (dpsi_val - real(ell,double) * psi_val/rmatch - two*c*rl*rmatch)/rl
    !a = psi_val/rl - c * rmatch * rmatch - b * rmatch
    !if(iprint>3) write(*,*) "Split norm coefficients, radius: ",a,b,c,rmatch
    ! -------------------------------------------------------------
    ! Cubic a + b*r + c*r*r + d*r*r*r (added d2psi = 0 at rc)
    ! -------------------------------------------------------------
    !if(ell>0) then
    !   a = zero
    !else
    !   a = psi_match(1)
    !end if
    !d = ( (psi_val - a) - dpsi_val * rmatch) / (rmatch * rmatch * rmatch)
    !c = -three * d * rmatch
    !b = dpsi_val - two * c * rmatch - three * d * rmatch * rmatch 
    !if(iprint>3) write(*,*) "Split norm coefficients, radius: ",a,b,c,d,rmatch
    ! -------------------------------------------------------------
    ! Cubic a + b*r + c*r*r + d*r*r*r (added d2psi matches zeta at rc)
    ! -------------------------------------------------------------
    !if(ell>0) then
    !   a = zero
    !else
    !   a = psi_match(1)
    !end if
    !d = half*( d2psi * rmatch * rmatch - real(2*(ell+1),double) * dpsi_val * rmatch &
    !     + real((ell+1)*(ell+2),double) * psi_val - two * a) / (rl * rmatch * rmatch * rmatch)
    !c = dpsi_val/(rl*rmatch) - real(ell+1,double) * psi_val/(rl*rmatch*rmatch) - &
    !     two * d * rmatch + a/(rmatch * rmatch)
    !b = psi_val/(rl*rmatch) - c * rmatch - d * rmatch * rmatch  - a/rmatch
    !if(iprint>3) write(*,*) "Split norm coefficients, radius: ",a,b,c,d,rmatch
    ! Now create the difference between smooth and original
    do i=1,nrc
       ! Cubic
       !psi(i) = psi_match(i) - (a + b*rr(i) + c*rr(i)*rr(i) + d*rr(i)*rr(i)*rr(i))
       ! General quadratic
       !psi(i) = psi_match(i) - (a + b*rr(i) + c*rr(i)*rr(i))
       ! Original quadratic
       psi(i) = psi_match(i) - (a + b*rr(i)*rr(i))
    end do
    ! Normalise
    norm = zero
    do i=1,nrc
       !norm = norm + rr(i)**(2*ell+2)*psi(i)*psi(i)*drdi_squared(i)
       norm = norm + rr(i)**(2*ell+2)*psi(i)*psi(i)*drdi(i)
    end do
    psi = psi/sqrt(norm)
  end subroutine find_split_norm

  ! Perform outward integration when using Vanderbilt-Kleinman-Bylander projectors
  ! We solve for the homogeneous equation (local potential) and then nproj inhomogeneous
  ! equations (one for each projector) which are then combined
  ! This serves as a layer between the eigenstate finder and numerov
  subroutine integrate_vkb_outwards(i_species,n_kink,ell,psi,f,n_crossings,n_nodes_inside)

    use datatypes
    use numbers
    use mesh, ONLY: rr, nmesh, drdi, drdi_squared, convert_r_to_i
    use pseudo_tm_info, ONLY: pseudo
    use GenComms, ONLY: cq_abort
    
    implicit none

    ! Passed variables
    integer :: i_species, n_kink, ell, n_crossings, n_nodes_inside
    real(double), dimension(nmesh) :: psi, f
    
    ! Local variables
    integer :: i_proj, j_proj, i, info, nproj_acc, j
    integer, dimension(4) :: ipiv ! Larger than needed
    real(double) :: integral, xkap, r_l_p_1, dy, dy_last
    real(double), allocatable, dimension(:) :: pot_vector, psi_h, s, integrand
    real(double), allocatable, dimension(:,:) :: pot_matrix, psi_inh

    !if(local_and_vkb%n_proj(ell)>0) then
       allocate(pot_matrix(local_and_vkb%n_proj(ell),local_and_vkb%n_proj(ell)),pot_vector(local_and_vkb%n_proj(ell)))
       pot_matrix = zero
       pot_vector = zero
       allocate(psi_h(nmesh),psi_inh(nmesh,local_and_vkb%n_proj(ell)))
       psi_h = zero
       psi_inh = zero
    !else
    !   allocate(psi_h(nmesh))
    !   psi_h = zero
    !end if
    if(ell == 0) then
       nproj_acc = 0
    else
       nproj_acc = sum(local_and_vkb%n_proj(0:ell-1))
    end if
    ! Homogeneous - standard numerov - but only to projector radius
    psi_h = zero
    n_kink = local_and_vkb%ngrid_vkb
    if(n_kink<5)     call convert_r_to_i(two,n_kink)
    if(iprint>5) write(*,fmt='("In integrate_vkb, outer limit is ",f18.10)') rr(n_kink)
    allocate(s(n_kink),integrand(n_kink))
    s = zero
    integrand = zero
    ! Initial values of the homogeneous solution from series expansion
    xkap = one
    psi_h(1) = rr(1)**real(ell+1,double)/sqrt(drdi(1))
    psi_h(2) = rr(2)**real(ell+1,double)/sqrt(drdi(2))
    r_l_p_1 = rr(1)**real(ell+1,double)
    call numerov(1,n_kink,nmesh,psi_h,rr,f,1,n_crossings,n_nodes_inside,xkap,n_kink,s) ! We want to integrate to max before final node
    ! Integral between solution and VKB projectors
    do i_proj = 1,local_and_vkb%n_proj(ell)
       do j=1,n_kink
          integrand(j) = psi_h(j)*sqrt(drdi(j))* &
               local_and_vkb%projector(j,i_proj,ell)
       end do
       call integrate_simpson(integrand,drdi(1:n_kink),n_kink,integral)
       if(iprint>6) write(*,fmt='("In integrate_vkb, homogeneous integral is ",f18.10)') integral
       pot_vector(i_proj) = integral*pseudo(i_species)%pjnl_ekb(nproj_acc + i_proj)
    end do
    !stop
    ! Inhomogeneous - add VKB projector as source term
    do i_proj = 1,local_and_vkb%n_proj(ell)
       s(1:n_kink) = two*local_and_vkb%projector(1:n_kink,i_proj,ell)*drdi(1:n_kink)*sqrt(drdi(1:n_kink))
       ! Series expansion based on projector (FIND APPROPRIATE REFERENCE !)
       xkap = two*local_and_vkb%projector(1,i_proj,ell)/(r_l_p_1*(six + four*real(ell,double)))
       psi_inh(1,i_proj) = xkap*rr(1)**real(ell+3,double)/sqrt(drdi(1))
       psi_inh(2,i_proj) = xkap*rr(2)**real(ell+3,double)/sqrt(drdi(2))
       call numerov(1,n_kink,nmesh,psi_inh(:,i_proj),rr,f,1,n_crossings,n_nodes_inside, &
            xkap,n_kink,s) 
       ! Integrals between this solution and VKB projectors
       do j_proj = 1,local_and_vkb%n_proj(ell)
          call integrate_simpson(psi_inh(1:n_kink,i_proj)*(sqrt(drdi(1:n_kink)))* & ! /rr(1:n_kink)
               local_and_vkb%projector(1:n_kink,j_proj,ell), drdi,n_kink,integral)
          if(iprint>6) write(*,fmt='("In integrate_vkb, inhomogeneous integral is ",f18.10)') integral
          pot_matrix(j_proj,i_proj) = -integral*pseudo(i_species)%pjnl_ekb(nproj_acc + j_proj)
       end do
       pot_matrix(i_proj,i_proj) = one + pot_matrix(i_proj,i_proj)
    end do
    !write(*,*) '# Pot mat and vec: ',pot_matrix, pot_vector
    ! Invert matrix
    if(local_and_vkb%n_proj(ell)>0) then
       call dgesv(local_and_vkb%n_proj(ell), 1, pot_matrix, local_and_vkb%n_proj(ell), ipiv, &
            pot_vector, local_and_vkb%n_proj(ell), info)
       if(info/=0) call cq_abort("Error from dgesv called in integrate_vkb_outward: ",info)
       if(iprint>6) write(*,fmt='("In integrate_vkb, coefficients are ",3f18.10)') pot_vector
    else
       pot_vector = zero
    end if
    ! Construct total outward wavefunction
    psi(1:n_kink) = psi_h(1:n_kink)
    do i_proj=1,local_and_vkb%n_proj(ell)
       psi(1:n_kink) = psi(1:n_kink) + pot_vector(i_proj)*psi_inh(1:n_kink,i_proj)
    end do
    ! Search for nodes
    n_crossings = 0
    dy = psi(2) - psi(1)
    do i=2,n_kink
       if(rr(i)>half.AND.psi(i)*psi(i-1)<zero) n_crossings = n_crossings+1
    end do
    !if(local_and_vkb%n_proj(ell)>0) then
       deallocate(pot_matrix, pot_vector, s, psi_h, psi_inh, integrand)
    !else
    !   deallocate(s, psi_h, integrand)
    !end if
  end subroutine integrate_vkb_outwards
  
  subroutine find_polarisation(i_species,en,ell,Rc,psi_l,psi_pol,energy,vha,vxc,pf_sign)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: rr, rr_squared, nmesh, alpha, make_mesh, beta, drdi, drdi_squared, convert_r_to_i
    use pseudo_tm_info, ONLY: pseudo
    
    implicit none

    integer :: i_species, ell, en
    real(double) :: energy, Rc, pf_sign
    real(double), dimension(nmesh) :: psi_l, psi_pol, vha, vxc

    ! Local variables
    real(double) :: g_temp, dy_L, dy_R
    real(double), dimension(:), allocatable :: f, g, s, potential
    integer :: classical_tp, i, n_crossings, n_nodes, n_loop, loop, nmax, n_kink, n_nodes_lower, n_nodes_upper
    real(double) :: l_half_sq, dx_sq_over_twelve, fac, norm, d_energy, e_lower, e_upper, df_cusp, cusp_psi, tol
    real(double) :: delta_energy_bracket, zval, l_l_plus_one, alpha_sq_over_four
    real(double) :: prefac, pf_low, pf_high

    do loop = 1,ell
       psi_l = psi_l*rr
    end do
    pf_low = 1.0e-3_double
    pf_high = 10000.0_double
    prefac = 1000.0_double ! Arbitrary initial value
    l_half_sq = real(ell+1,double) + half
    l_half_sq = l_half_sq*l_half_sq
    l_l_plus_one = real((ell+2)*(ell+1),double)
    ! NB The values of n and l we have here are for the perturbed function, so this is correct
    n_nodes = en - ell - 1
    if(iprint>2) write(*,fmt='(2x,"For this polarisation function, we require ",i2," nodes")') n_nodes
    ! This is needed when we are perturbing a valence state with a semi-core state
    if(n_nodes>0 .and. pf_sign>zero .and. psi_l(1)<zero) then
       pf_sign = -one
       if(iprint>5) write(*,fmt='(4x,"Changing sign of function as required")')
    end if
    zval = pseudo(i_species)%z
    dx_sq_over_twelve = alpha*alpha/twelve
    alpha_sq_over_four = alpha*alpha/four
    allocate(f(nmesh),g(nmesh), s(nmesh), potential(nmesh))
    !call radial_hartree(nmesh,semilocal%charge,vha)
    !call vxc_pz_ca(nmesh,rr,semilocal%charge,vxc)
    e_lower = zero!-zval*zval/real(en*en,double)
    ! It is both correct and better to use V_{l+1} in semi-local but we allow V_{l} both for
    ! compatibility with old Siesta pseudopotentials and for cases where we do not have semi-local
    ! potentials for l+1
    if((ell+1>pseudo(i_species)%lmax).OR.flag_use_Vl) then
       if(iprint>2) then
          if(ell+1>pseudo(i_species)%lmax) then
             write(*,*) 'lmax is ',pseudo(i_species)%lmax,' so perturbing using l not l+1'
          else
             write(*,*) 'Using V_{l} not V_{l+1} for perturbation'
          end if
       end if
       !l_l_plus_one = real((ell)*(ell+1),double)
       do i=1,nmesh
          potential(i) = local_and_vkb%semilocal_potential(i,ell) + vha(i) + vxc(i)
       end do
    else 
       do i=1,nmesh
          potential(i) = local_and_vkb%semilocal_potential(i,ell+1) + vha(i) + vxc(i)
       end do
    end if
    n_loop = 100
    tol = 1.0e-8_double
    ! Energy bounds - allow for unbound states
    nmax = nmesh ! Adjust later to confine
    ! Test
    call convert_r_to_i(Rc,nmax)
    nmax = nmax - 1
    ! Numerov
    do i=1,nmax
       g(i) = (drdi_squared(i)*(two*(energy - potential(i))-l_l_plus_one/rr_squared(i)) - alpha_sq_over_four)/twelve
       f(i) = one + g(i)
       s(i) = -two*rr(i)*rr(i)*psi_l(i) ! Need r*phi here
    end do
    ! Iterate to find correct initial boundary condition
    do loop = 1,n_loop
       ! Small radius solution - dV/dr is zero near centre
       ! Use prefac to set psi(1) and gradient
       psi_pol(1) = pf_sign*prefac*rr(1)**(ell+2)!*(one - zval*r(1)/(real(ell+1,double)))
       psi_pol(1) = psi_pol(1)/sqrt(drdi(1))
       psi_pol(2) = pf_sign*prefac*rr(2)**(ell+2)!*(one - zval*r(2)/(real(ell+1,double)))
       psi_pol(2) = psi_pol(2)/sqrt(drdi(2))
       ! We need crossings counted - if there are crossings, then increase, else decrease
       n_crossings = 0
       do i=2,nmax-1
          psi_pol(i+1) = (two*psi_pol(i)*(one-five*g(i)) - psi_pol(i-1)*f(i-1) + (s(i+1)+ten*s(i)+s(i-1))/twelve)/f(i+1)
          if(psi_pol(i+1)*psi_pol(i)<zero) n_crossings = n_crossings + 1
       end do
       if(iprint>5) write(*,fmt='("Crossings found: ",i3)') n_crossings
       if(n_crossings > n_nodes) then ! Increase
          pf_low = prefac
          prefac = half*(pf_low + pf_high)
          if(iprint>5) write(*,fmt='("More crossings than nodes required: ",2i3)') n_crossings, n_nodes
          if(iprint>5) write(*,fmt='("Polarisation loop ",i3," brackets and prefactor: ",3f13.5)')loop,pf_low,prefac,pf_high
          !do i=1,nmax
          !   write(13,*) rr(i),psi_pol(i),-rr(i)*rr(i)*psi_l(i)
          !end do
          !write(13,*) '&'
          cycle
       else
          pf_high = prefac
          prefac = half*(pf_low+pf_high)
          !do i=1,nmax
          !   write(13,*) rr(i),psi_pol(i),-rr(i)*rr(i)*psi_l(i)
          !end do
          !write(13,*) '&'
          if(iprint>5) write(*,fmt='("Polarisation loop ",i3," brackets and prefactor: ",3f13.5)') loop,pf_low,prefac,pf_high
       end if
       ! Normalise
       norm = zero
       do i=1,nmax
          norm = norm + psi_pol(i)*psi_pol(i)*drdi_squared(i)
       end do
       psi_pol = psi_pol/sqrt(norm)
       if(pf_high - pf_low<1e-3_double.AND.abs(psi_pol(nmax))<1e-3_double) then
          !if(iprint>4) write(*,*) '# Found prefac'
          exit
       end if
    end do
    if(loop>=n_loop) then
       if(iprint>5) then
          do i=1,nmax
             write(50,*) rr(i),psi_pol(i),-rr(i)*rr(i)*psi_l(i)
          end do
          call flush(50)
       end if
       call cq_abort("ERROR in perturbative polarisation routines - prefactor not found")
    else
       if(iprint>2) write(*,fmt='("Finished perturbative polarisation search with prefactor ",f13.5)') prefac
    end if
    !do i=1,nmax
    !   write(12,*) rr(i),psi_l(i),psi_pol(i)!,psi_pol(i)*sqrt(drdi(i))/rr(i)
    !end do
    !write(12,*) '&'
    ! Rescale - remove factor of sqrt r
    do i=1,nmax
       psi_pol(i) = psi_pol(i)*sqrt(drdi(i))/rr(i)
       !write(13,*) rr(i),(two*rr_squared(i)*potential(i)+l_half_sq)*alpha*alpha, &
       !     alpha*alpha*rr_squared(i)
    end do
    !do i=1,nmesh
    !   write(25+ell+1,*) rr(i),psi_pol(i)
    !end do
    !write(25+ell+1,*) '&'    
    do loop = 1,ell+1
       psi_pol = psi_pol/rr
       if(loop<ell+1) psi_l = psi_l/rr
    end do
    deallocate(f, g, s, potential)
  end subroutine find_polarisation

  ! I've added the possibility to include a source term, for the VKB inhomogeneous solutions
  ! In this case we won't count nodes (done in the calling routine) so we'll need to make sure
  ! that the limits are correctly set
  subroutine numerov(n_st,n_pts,n_tot,y,r,f,direction,n_crossings,n_nodes_inside,prefac,n_source,source)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: drdi
    
    implicit none

    ! Passed variables
    integer :: n_st, n_pts, direction, n_tot
    integer :: n_crossings, n_nodes_inside
    real(double) :: energy, prefac
    real(double), dimension(n_tot) :: y, r, f
    integer :: n_source
    real(double), OPTIONAL, dimension(n_source) :: source
    
    ! Local variables
    integer :: i, j, n_end, n_stop_rl, en
    real(double) :: l_half_sq, dx_sq_over_twelve, g_temp, dy, dy_last, expon
    !real(double), dimension(:), allocatable :: y
    logical :: flag_source

    flag_source = .false.
    if(present(source)) flag_source = .true.
    n_end = n_st + n_pts - 1
    if(direction==1) then
       ! Set first two points - seies expansion
       !expon = real(ell+1,double)
       !y(1) = prefac*r(1)**expon!*(one - z*r(1)/(real(ell+1,double)))
       !y(2) = prefac*r(2)**expon!*(one - z*r(2)/(real(ell+1,double)))
       !y(1) = y(1)/sqrt(drdi(1))
       !y(2) = y(2)/sqrt(drdi(2))
       !%%! if(ell==0) then
       !%%!    expon = real(ell+1,double)
       !%%!    !y(1) = r(1)**expon*(one - z*r(1)/(real(ell+1,double)) + two*(z/real(en,double))**1.5_double*r(1)*r(1))
       !%%!    y(1) = r(1)**expon*(one - z*r(1)/(real(ell+1,double)) + (200.0_double/real(en,double))*r(1)*r(1))
       !%%!    y(1) = y(1)/sqrt(drdi(1))
       !%%!    y(2) = r(2)**expon*(one - z*r(2)/(real(ell+1,double)) + (200.0_double/real(en,double))*r(2)*r(2))
       !%%!    y(2) = y(2)/sqrt(drdi(2))
       !%%! else
       !%%!    y(1) = r(1)**expon*(one - z*r(1)/(real(ell+1,double)))
       !%%!    y(1) = y(1)/sqrt(drdi(1))
       !%%!    y(2) = r(2)**expon*(one - z*r(2)/(real(ell+1,double)))
       !%%!    y(2) = y(2)/sqrt(drdi(2))
       !%%! end if
       !y(1) = f(1)/(twelve - ten*f(1))
       !dy = y(1)
       !y(2) = ( (twelve-ten*f(1)) * y(1) )/f(2)
       !y(2) = one
       !y(3) = r(3)**(ell+1)*(one - z*r(3)/(real(ell+1,double)))
       !y(3) = y(3)/sqrt(drdi(3))
       !y(4) = r(4)**(ell+1)*(one - z*r(4)/(real(ell+1,double)))
       !y(4) = y(4)/sqrt(drdi(4))
       !dy = y(4) - y(3)
       dy = y(n_st+1) - y(n_st)
       ! Solve for y
       if(flag_source) then
          do i=n_st+1,n_end-1
             dy_last = dy
             y(i+1) = ( (twelve - ten*f(i)) * y(i) - f(i-1)*y(i-1) + &
                  (source(i+1) + ten*source(i) + source(i-1))/twelve)/f(i+1)
             dy = y(i+1) - y(i)
             !if(y(i)*y(i+1)<zero.OR.(abs(y(i+1))<1e-16_double)) then ! Crossing the x-axis
             !   n_crossings = n_crossings + 1
             !end if
             !if(dy*dy_last<zero.AND.n_crossings>=n_nodes_inside) then
             !   n_pts = i ! Maximum
             !   exit
             !end if
          end do
       else
          do i=n_st+1,n_end-1
             dy_last = dy
             y(i+1) = ( (twelve - ten*f(i)) * y(i) - f(i-1)*y(i-1) )/f(i+1)
             dy = y(i+1) - y(i)
             if(y(i)*y(i+1)<zero.OR.(abs(y(i+1))<1e-16_double)) then ! Crossing the x-axis
                n_crossings = n_crossings + 1
             end if
             n_pts = i ! Maximum
             if(dy*dy_last<zero.AND.n_crossings>=n_nodes_inside) then
                n_pts = i ! Maximum
                exit
             end if
          end do
       end if
       !n_pts = i ! Point of maximum
    else if(direction==-1) then
       !write(*,*) '# Inward: ',n_end,n_st
       y(n_end) = zero!alpha
       y(n_end-1) = 0.001_double!one ! Arbitrary - we scale later 
       !
       ! inward integration 
       !
       do i = n_end-1,n_st,-1
          y(i-1) = ( (twelve - ten*f(i)) * y(i) - f(i+1)*y(i+1) )/f(i-1)
          if(y(i)*y(i-1)<zero.OR.(abs(y(i-1))<1e-16_double)) then
             n_crossings = n_crossings + 1
          end if
          ! This rescaling may be needed - but not at the moment
          if (y(i-1) > 1.0e7_double) then
             y(i-1:n_end) = y(i-1:n_end)/y(i-1)
          end if
       end do
    end if
  end subroutine numerov

  ! Also uses numerov - may be possible to combine in future
  subroutine radial_hartree(n_tot,rho,vha)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: alpha, drdi, rr, rr_squared, sqrt_drdi
    
    implicit none

    ! Passed variables
    integer :: n_tot
    real(double), dimension(n_tot) :: rho,vha
   
    ! Local variables
    integer :: i
    real(double) :: dx_sq_over_twelve, qtot, V0, f, g, prefac
    real(double), dimension(:), allocatable :: y, s

    ! Set constants
    dx_sq_over_twelve = alpha*alpha/twelve
    allocate(y(n_tot),s(n_tot))
    y = zero
    ! Set first two points
    ! Point 1: integrate outwards 4.pi.r on log mesh
    qtot = zero
    V0 = zero
    vha = zero
    call integrate_simpson(rho*rr_squared,drdi,n_tot,qtot)
    call integrate_simpson(rho*rr,drdi,n_tot,V0)
    do i=1,n_tot
    !%%!    ! The charge density is actually 4.pi.rho for Hamann
    !%%!    ! The Siesta density is given as 4.pi.r^2.rho but read_module divides by r^2.
    !%%!    ! NB \int dr -> \int (dr/di) di with (dr/di) found in mesh_module.f90
    !%%!    ! Total charge
    !%%!    qtot = qtot + rho(i)*rr_squared(i)*drdi(i)
    !%%!    ! Potential at r=0
    !%%!    V0 = V0 + rho(i)*drdi(i)*rr(i)
    !%%!    ! array for solution, below
       s(i) = -rr(i)*rho(i)*drdi(i)*sqrt_drdi(i)
    end do
    if(iprint>4) write(*,fmt='("In Hartree, total valence charge is ",f12.5)') qtot
    ! Numerov
    g = -alpha*alpha/four ! alpha is d/dr(dr/di) for both meshes
    f = one + g/twelve
    ! Small radius solution - dV/dr is zero near centre
    y(1) = V0*rr(1)/sqrt_drdi(1)
    y(2) = V0*rr(2)/sqrt_drdi(2)
    do i=2,n_tot-1
       y(i+1) = (two*y(i)*(one-five*g/twelve) - y(i-1)*f + (s(i+1)+ten*s(i)+s(i-1))/twelve)/f
    end do
    V0 = (Qtot - sqrt_drdi(n_tot)*y(n_tot))/rr(n_tot) ! Solution of homogeneous equation so that we reach asymptote Z/r
    !write(*,*) '# Correction shift: ',V0
    do i=1,n_tot
       vha(i) = y(i)*sqrt_drdi(i)/rr(i) + V0
       !write(25,*) rr(i),vha(i), rho(i)
    end do
    deallocate(y,s)
  end subroutine radial_hartree

  subroutine integrate_simpson(f,metric,np,total)

    use numbers
    
    implicit none

    ! Passed variables
    integer :: np
    real(double), dimension(np) :: f, metric
    real(double) :: total

    ! Local variables
    integer :: n, i

    ! Even or odd ?
    if(mod(np,2)==1) then
       n = np
    else
       n = np - 3 ! Fix at end
    end if

    ! Integrate with simpson's rule
    total = zero
    ! Four thirds terms
    do i=2,n-1,2
       total = total + f(i)*metric(i)
    end do
    total = total * two
    ! Two thirds terms
    do i=3,n-2,2
       total = total + f(i)*metric(i)
    end do
    total = total * two
    ! End points
    total = total + f(1)*metric(1)
    total = total + f(n)*metric(n)
    total = total*third
    if(mod(np,2)==0) then ! Correct final segment
       i=np-3
       total = total + half*three_quarters*(f(i)*metric(i) + f(i+3)*metric(i+3) + &
            three*(f(i+1)*metric(i+1) + f(i+2)*metric(i+2)) )
    end if
  end subroutine integrate_simpson
end module schro
