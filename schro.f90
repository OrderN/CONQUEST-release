module schro

  use datatypes
  use pao_info
  use global_module, ONLY: iprint
  
  implicit none

!  type(potential), allocatable, dimension(:) :: semilocal
  type(potential_vkb), allocatable, dimension(:) :: local_and_vkb

  integer :: pseudo_type

  type(valence_info), allocatable, dimension(:) :: val
  logical :: flag_user_specified = .false.
  integer, dimension(0:6) :: n_proj ! Number of KB projectors for each l (declare as 6 to avoid alloc)
  logical :: flag_run_debug = .false.
  logical :: flag_plot_output

  integer :: n_debug_run, l_debug_run
  real(double) :: E_debug_run
  
  real(double), allocatable, dimension(:) :: deltaE_large_radius! = 0.00073498_double
  real(double), allocatable, dimension(:) :: deltaE_small_radius! = 0.073498_double
  
  ! Useful parameters to improve code readability
  integer, parameter :: pao_cutoff_energies = 1
  integer, parameter :: pao_cutoff_radii = 2
  integer, parameter :: pao_cutoff_default = 3 ! This will be set properly at some point
contains

  subroutine make_paos

    use datatypes
    use numbers
    use species_module, ONLY: n_species
    use pseudo_tm_info, ONLY: pseudo, rad_alloc
    use mesh, ONLY: nmesh, make_mesh, rr, nmesh_reg, rmesh_reg, interpolate, make_mesh_reg, rr_squared, delta_r_reg
    use global_module, ONLY: flag_pcc_global, iprint
    use GenComms, ONLY: cq_abort
    use radial_xc, ONLY: get_vxc
    use input_module, ONLY: io_assign, io_close
    use write, ONLY: write_pao_plot, pte
    
    implicit none

    integer :: i_species, ell, en, i, j, k, zeta, n_shells, i_shell, grid_size, n_val, lun
    real(double), allocatable, dimension(:) :: psi_reg, x_reg, vha, vxc, vha_reg, psi, atomic_rho, &
         small_cutoff, large_cutoff, cutoffs, vha_conf
    real(double) :: energy, max_cutoff, small_energy, large_energy
    real(double) :: radius_small, radius_med, radius_large, sum_small, sum_large
    real(double), dimension(3) :: scaling
    !logical :: flag_plot = .true.
    
    ! The default cutoff radii are based on energies which scale by ten
    scaling(1) = one
    do zeta = 2,3
       scaling(zeta) = ten*scaling(zeta-1)
    end do
    do i_species = 1,n_species
       write(*,fmt='(/"Generating PAOs for ",a2/)') pte(pseudo(i_species)%z)
       nmesh = local_and_vkb(i_species)%ngrid
       call make_mesh(i_species)
       !allocate(psi_reg(nmesh_reg))
       allocate(vha(nmesh),vxc(nmesh),vha_conf(nmesh))
       allocate(psi(nmesh),atomic_rho(nmesh))
       ! We need max cutoff to work out largest PAO radius (for neutral atom)
       max_cutoff = zero
       atomic_rho = zero
       ! Find VHa and VXC for atom to screen the semi-local potentials
       call radial_hartree(nmesh,local_and_vkb(i_species)%charge,vha,i_species)
       if(flag_pcc_global) & 
            local_and_vkb(i_species)%charge = local_and_vkb(i_species)%charge + local_and_vkb(i_species)%pcc
       call get_vxc(nmesh,rr,local_and_vkb(i_species)%charge,val(i_species)%functional,vxc)
       if(flag_run_debug) then
          call io_assign(lun)
          open(unit=lun,file=pte(pseudo(i_species)%z)//"_AtomicVHaVXC.dat")
          do i=1,nmesh
             write(lun,fmt='(3e18.10)') rr(i),vha(i),vxc(i)
          end do
          call io_close(lun)
       end if
       ! Polarisation will be done separately at the end
       n_shells = val(i_species)%n_occ
       ! Find energy of valence states without confinement
       if(iprint>2) write(*,fmt='(2x,"Finding unconfined energies for valence states")')
       do i_shell = 1, val(i_species)%n_occ
          if(iprint>2) write(*,fmt='(2x,"n=",i2," l=",i2)') val(i_species)%n(i_shell), val(i_species)%l(i_shell)
          radius_large = rr(nmesh)
          ell = val(i_species)%l(i_shell)
          en = val(i_species)%npao(i_shell)
          large_energy = val(i_species)%en_ps(i_shell)
          call find_eigenstate_and_energy_vkb(i_species,en,ell,radius_large, psi,large_energy,vha,vxc)
          val(i_species)%en_pao(i_shell) = large_energy
       end do
       if(iprint>0) then
          write(*,fmt='(2x,"Unconfined valence state energies")')
          write(*,fmt='(2x,"  n  l         AE energy        PAO energy")')
          do i_shell = 1, val(i_species)%n_occ
             ell = val(i_species)%l(i_shell)
             en = val(i_species)%n(i_shell)
             write(*,fmt='(2x,2i3,2f18.10)') en, ell, val(i_species)%en_ps(i_shell), &
                  val(i_species)%en_pao(i_shell)
          end do
       end if
       ! Find default radii or radii from cutoffs
       write(*,fmt='(/2x,"Setting cutoff radii"/)')
       if(flag_default_cutoffs) then ! Work out radii
          call find_default_cutoffs(i_species,vha,vxc)
       else if(paos(i_species)%flag_cutoff==pao_cutoff_energies) then ! Work out radii from energy shifts
          do i_shell = 1,val(i_species)%n_occ
             ell = val(i_species)%l(i_shell)
             en = val(i_species)%npao(i_shell)
             do zeta = 1,paos(i_species)%nzeta(i_shell)
                !if(iprint>2) write(*,*) '# shell, n, l, zeta: ',i_shell, ell, en, zeta
                call find_radius_from_energy(i_species,paos(i_species)%npao(i_shell), paos(i_species)%l(i_shell), &
                     paos(i_species)%cutoff(zeta,i_shell), val(i_species)%en_ps(i_shell) + paos(i_species)%energy(zeta,i_shell), &
                     vha, vxc, .false.)
                ! Now ensure that we have exact cutoff - effectively round up (2+ not 1+)
                grid_size = 2+floor(paos(i_species)%cutoff(zeta,i_shell)/delta_r_reg)
                paos(i_species)%cutoff(zeta,i_shell) = delta_r_reg*real(grid_size-1,double)
                !if(iprint>2) write(*,*) '# Setting cutoff to ',paos(i_species)%cutoff(zeta,i_shell)
             end do
          end do
          do i_shell = val(i_species)%n_occ+1,paos(i_species)%n_shells
             if(i_shell>val(i_species)%n_occ+1) write(*,fmt='(2x,"Dangerous to set unoccupied shell radii via energy !")')
             if(paos(i_species)%inner_shell>0) then
                do zeta = 1,paos(i_species)%nzeta(i_shell)
                   paos(i_species)%cutoff(zeta,i_shell) = paos(i_species)%cutoff(zeta,paos(i_species)%inner_shell)
                end do
             else
                do zeta = 1,paos(i_species)%nzeta(i_shell)
                   paos(i_species)%cutoff(zeta,i_shell) = paos(i_species)%cutoff(zeta,i_shell-1)
                end do
             end if
          end do
       end if
       write(*,fmt='(4x,"Cutoff radii for PAOs")')
       write(*,fmt='(4x,"  n  l  z        R")')
       do i_shell = 1,paos(i_species)%n_shells
          ell = paos(i_species)%l(i_shell)
          en = paos(i_species)%n(i_shell)
          if(paos(i_species)%flag_perturb_polarise.AND.i_shell==paos(i_species)%n_shells) then
             if(en<3) en = en+1
             do zeta = 1,paos(i_species)%nzeta(i_shell)
                write(*,fmt='(4x,3i3,f18.10)') en,ell+1,zeta,paos(i_species)%cutoff(zeta,i_shell)
             end do
          else
             do zeta = 1,paos(i_species)%nzeta(i_shell)
                write(*,fmt='(4x,3i3,f18.10)') en,ell,zeta,paos(i_species)%cutoff(zeta,i_shell)
             end do
          end if
       end do
       if(flag_run_debug) then
          if(n_debug_run>0) then
             call find_eigenstate_and_radius_vkb(i_species,n_debug_run,l_debug_run,radius_med,&
                  psi,E_debug_run,vha,vxc,.false.)
          end if
       end if
       ! Solve for PAOs
       write(*,fmt='(/2x,"Solving for PAOs")')
       do i_shell = 1,paos(i_species)%n_shells
          ell = paos(i_species)%l(i_shell)
          en = paos(i_species)%npao(i_shell)
          do zeta = 1,paos(i_species)%nzeta(i_shell)
             allocate(paos(i_species)%psi(zeta,i_shell)%f(nmesh))
             !write(*,*) '# Occ for shell: ',paos(i_species)%occ(i_shell)
             if(paos(i_species)%flag_perturb_polarise.AND.i_shell==paos(i_species)%n_shells) then ! Polarisation functions
                ell = paos(i_species)%l(i_shell)!-1 !i_shell-1)
                en = paos(i_species)%npao(i_shell)!i_shell-1)
                if(iprint>2) write(*,fmt='(2x,"Perturbative polarisation")')
                if(iprint>2) write(*,fmt='(2x,"Species ",i2," n=",i2," l=",i2," zeta=",i2, "Rc=",f4.1)') &
                     i_species, en, ell, zeta, paos(i_species)%cutoff(zeta,i_shell)
                call find_polarisation(i_species,en,ell,paos(i_species)%cutoff(zeta,i_shell-1),&
                     paos(i_species)%psi(zeta,i_shell-1)%f,paos(i_species)%psi(zeta,i_shell)%f,&
                     paos(i_species)%energy(zeta,i_shell-1),vha,vxc)
             else if(i_shell>val(i_species)%n_occ) then  ! Polarisation shells
             !else if(i_shell==paos(i_species)%n_shells) then
                if(iprint>2) write(*,fmt='(2x,"Species ",i2," n=",i2," l=",i2," zeta=",i2, "Rc=",f4.1," pol ",f4.1)') &
                     i_species, en, ell, zeta, paos(i_species)%cutoff(zeta,i_shell), &
                     paos(i_species)%energy(1,val(i_species)%n_occ)
                !if(iprint>2) write(*,*) '# Species, n, l, zeta, cutoff: ',i_species, en, ell, zeta, &
                !     paos(i_species)%cutoff(zeta,i_shell)," pol ",paos(i_species)%energy(1,val(i_species)%n_occ)
                paos(i_species)%energy(zeta,i_shell) = paos(i_species)%energy(1,val(i_species)%n_occ)
                call find_eigenstate_and_energy_vkb(i_species,en,ell,paos(i_species)%cutoff(zeta,i_shell),&
                     paos(i_species)%psi(zeta,i_shell)%f,paos(i_species)%energy(zeta,i_shell), &
                     vha,vxc) 
             else
                if(iprint>2) write(*,fmt='(2x,"Species ",i2," n=",i2," l=",i2," zeta=",i2, "Rc=",f4.1)') &
                     i_species, en, ell, zeta, paos(i_species)%cutoff(zeta,i_shell)
                !if(iprint>2) write(*,*) '# Species, n, l, zeta, cutoff: ',i_species, en, ell, zeta, &
                !     paos(i_species)%cutoff(zeta,i_shell)
                large_energy = val(i_species)%en_ps(i_shell) + paos(i_species)%energy(zeta,i_shell)
                call find_eigenstate_and_energy_vkb(i_species,en,ell,paos(i_species)%cutoff(zeta,i_shell),&
                     paos(i_species)%psi(zeta,i_shell)%f,large_energy, vha,vxc)
                if(iprint>2) write(*,*) '# Energy shift required: ',large_energy - val(i_species)%en_ps(i_shell)
                paos(i_species)%energy(zeta,i_shell) = large_energy
             end if
             if(paos(i_species)%cutoff(zeta,i_shell)>max_cutoff) max_cutoff = paos(i_species)%cutoff(zeta,i_shell)
          end do
       end do
       write(*,fmt='(/2x,"Interpolating onto regular mesh")')
       do i_shell = 1,paos(i_species)%n_shells
          ell = paos(i_species)%l(i_shell)
          en = paos(i_species)%n(i_shell)
          if(paos(i_species)%flag_perturb_polarise.AND.i_shell==paos(i_species)%n_shells) then
             ell = ell + 1
             if(en<3) en = en + 1
          end if
          do zeta = 1,paos(i_species)%nzeta(i_shell)
             nmesh_reg = paos(i_species)%cutoff(zeta,i_shell)/delta_r_reg + 1
             allocate(paos(i_species)%psi_reg(zeta,i_shell)%f(nmesh_reg))
             allocate(paos(i_species)%psi_reg(zeta,i_shell)%x(nmesh_reg))
             !if(iprint>2) write(*,*) '# regular mesh ',paos(i_species)%cutoff(zeta,i_shell)
             call make_mesh_reg(paos(i_species)%cutoff(zeta,i_shell))
             paos(i_species)%psi_reg(zeta,i_shell)%x = rmesh_reg
             paos(i_species)%psi_reg(zeta,i_shell)%delta = delta_r_reg
             paos(i_species)%psi_reg(zeta,i_shell)%n = nmesh_reg
             ! Use this for atomic density
             psi = zero
             psi = paos(i_species)%psi(zeta,i_shell)%f
             if(ell>0) then
                do i=1,ell
                   !write(*,*) '# Scaling psi by r for density'
                   psi = psi * rr
                end do
             end if
             ! Allocate regular mesh
             if(zeta==1.AND.i_shell<=val(i_species)%n_shells) atomic_rho = atomic_rho + val(i_species)%occ(i_shell)*psi*psi
             ! Interpolate
             !if(iprint>2) write(*,*) '# interpolate ',nmesh_reg, nmesh
             call interpolate(rmesh_reg,paos(i_species)%psi_reg(zeta,i_shell)%f,nmesh_reg,&
                  paos(i_species)%psi(zeta,i_shell)%f,nmesh,zero)
             if(flag_plot_output) call write_pao_plot(pseudo(i_species)%z,rmesh_reg,paos(i_species)%psi_reg(zeta,i_shell)%f, &
                  nmesh_reg, en,ell,zeta)
          end do
       end do
       ! Interpolate the logarithmic grid pseudopotential terms
       write(*,fmt='(/2x,"Calculating potentials")')
       !if(iprint>2) write(*,*) '# Moving to VKB projectors'
       ! VKB projectors
       ! Round up VKB cutoff
       grid_size = floor(local_and_vkb(i_species)%r_vkb/delta_r_reg) + 2
       local_and_vkb(i_species)%r_vkb = (grid_size-1)*delta_r_reg
       call make_mesh_reg(local_and_vkb(i_species)%r_vkb)
       !if(iprint>2) write(*,*) '# Cutoff and size: ',local_and_vkb(i_species)%r_vkb,local_and_vkb(i_species)%ngrid_vkb
       j = 0
       if(iprint>1) write(*,fmt='(/4x,"VKB projectors")')
       do ell = 0, pseudo(i_species)%lmax
          !if(iprint>2) write(*,*) '# l: ',ell
          do i=1,n_proj(ell)
             !if(iprint>2) write(*,*) '# projector number: ',i
             j = j+1
             call rad_alloc(pseudo(i_species)%pjnl(j),grid_size)
             pseudo(i_species)%pjnl(j)%delta = delta_r_reg
             ! Scale by r**(l+1)
             do k=0,ell
                local_and_vkb(i_species)%projector(:,i,ell) = local_and_vkb(i_species)%projector(:,i,ell)/rr
             end do
             ! Grid point, projector, ell
             call interpolate(rmesh_reg,pseudo(i_species)%pjnl(j)%f,grid_size, &
                  local_and_vkb(i_species)%projector(:,i,ell),local_and_vkb(i_species)%ngrid_vkb,zero)
             pseudo(i_species)%pjnl(j)%cutoff = local_and_vkb(i_species)%r_vkb
          end do
       end do
       ! Partial core charge
       if(pseudo(i_species)%flag_pcc) then
          if(iprint>1) write(*,fmt='(/4x,"Partial core correction")')
          !if(iprint>2) write(*,*) '# PCC'
          pseudo(i_species)%chpcc%delta = delta_r_reg
          ! Find core charge cutoff
          do i=1,nmesh
             ! Find cutoff radius: two successive points less than 1e-8 
             if(local_and_vkb(i_species)%pcc(i)<RD_ERR.AND.i<nmesh) then
                if(local_and_vkb(i_species)%pcc(i+1)<RD_ERR) then
                   pseudo(i_species)%chpcc%cutoff = rr(i)
                   exit
                end if
             end if
          end do
          grid_size = floor(pseudo(i_species)%chpcc%cutoff/pseudo(i_species)%chpcc%delta) + 2
          pseudo(i_species)%chpcc%cutoff = real(grid_size-1,double)*delta_r_reg
          !if(iprint>3) write(*,*) '# PCC cutoff ',pseudo(i_species)%chpcc%cutoff
          call rad_alloc(pseudo(i_species)%chpcc,grid_size)
          call make_mesh_reg(pseudo(i_species)%chpcc%cutoff)!,grid_size)
          call interpolate(rmesh_reg,pseudo(i_species)%chpcc%f,grid_size, &
               local_and_vkb(i_species)%pcc,local_and_vkb(i_species)%ngrid,zero)
       end if
       ! Local - the radius is taken as largest PAO for compatibility with NA
       if(iprint>1) write(*,fmt='(/4x,"Local potential")')
       !if(iprint>2) write(*,*) '# Local'
       grid_size = max_cutoff/delta_r_reg + 1
       call rad_alloc(pseudo(i_species)%vlocal,grid_size)
       pseudo(i_species)%vlocal%delta = delta_r_reg
       pseudo(i_species)%vlocal%cutoff = (grid_size-1)*delta_r_reg
       call make_mesh_reg(pseudo(i_species)%vlocal%cutoff)!,grid_size)
       call interpolate(rmesh_reg,pseudo(i_species)%vlocal%f,grid_size, &
            local_and_vkb(i_species)%local,local_and_vkb(i_species)%ngrid)
       ! Build the neutral atom potential as sum of local and pseudo-atomic hartree potentials
       if(iprint>1) write(*,fmt='(/4x,"Neutral atom potential")')
       !if(iprint>2) write(*,*) '# NA'
       ! Find pseudo-atomic hartree potential from pseudo-atomic density
       call radial_hartree(nmesh,atomic_rho,vha_conf,i_species)
       do i=1,nmesh
          !write(60,*) rr(i),vha_conf(i),local_and_vkb(i_species)%local(i) + vha_conf(i), atomic_rho(i) 
          vha_conf(i) = vha_conf(i) + local_and_vkb(i_species)%local(i)
       end do
       !if(iprint>2) write(*,*) '# Building VNA'
       ! Allocate space and assign variables - use same delta as for local potential
       !grid_size = max_cutoff/pseudo(i_species)%vlocal%delta + 1
       call rad_alloc( pseudo(i_species)%vna, grid_size )
       pseudo(i_species)%vna%delta = delta_r_reg
       pseudo(i_species)%vna%cutoff = max_cutoff
       ! Interpolate onto regular mesh
       !allocate(vha_reg(pseudo(i_species)%vna%n))
       !call make_mesh_reg(pseudo(i_species)%vna%cutoff,pseudo(i_species)%vna%n)
       call interpolate(rmesh_reg,pseudo(i_species)%vna%f,pseudo(i_species)%vna%n,vha_conf,nmesh,zero)
       if(iprint>1) write(*,fmt='(/4x,"XC potential")')
       vxc = zero
       energy = zero
       call get_vxc(nmesh,rr,atomic_rho,val(i_species)%functional,vxc,energy)
       if(iprint>2) write(*,*) '# XC energy w/o core: ',energy
       !do i=1,nmesh
       !   write(70,*) rr(i),vxc(i),atomic_rho(i)
       !end do
       vxc = zero
       energy = zero
       if(flag_pcc_global) & 
            atomic_rho = atomic_rho + local_and_vkb(i_species)%pcc
       call get_vxc(nmesh,rr,atomic_rho,val(i_species)%functional,vxc,energy)
       if(iprint>2) write(*,*) '# XC energy with core: ',energy
       !do i=1,nmesh
       !   write(71,*) rr(i),vxc(i),atomic_rho(i)
       !end do
       !%%! ! The local potential is normally defined to a specific cut-off, beyond which it goes
       !%%! ! as zval/r, whereas the pseudo-atomic density extends to the PAO range, so we need to
       !%%! ! add the long-range asymptote of the local potential to make the neutral atom potential
       !%%! if(pseudo(i_species)%vlocal%n<pseudo(i_species)%vna%n) then
       !%%!    pseudo(i_species)%vna%f(1:pseudo(i_species)%vlocal%n) = vha_reg(1:pseudo(i_species)%vlocal%n) + &
       !%%!         pseudo(i_species)%vlocal%f(1:pseudo(i_species)%vlocal%n)
       !%%!    do i=pseudo(i_species)%vlocal%n+1,pseudo(i_species)%vna%n
       !%%!       pseudo(i_species)%vna%f(i) = vha_reg(i) - &
       !%%!            pseudo(i_species)%zval/(pseudo(i_species)%vna%delta*real(i-1,double))
       !%%!    end do
       !%%! else
       !%%!    pseudo(i_species)%vna%f(1:pseudo(i_species)%vna%n) = vha_reg(1:pseudo(i_species)%vna%n) + &
       !%%!         pseudo(i_species)%vlocal%f(1:pseudo(i_species)%vna%n)
       !%%! end if
       deallocate(vha,vxc,vha_conf)
       deallocate(psi,atomic_rho)
       write(*,fmt='(/2x,"Finished ",a2)') pte(pseudo(i_species)%z)
    end do
  end subroutine make_paos

  ! For the default basis, find the cutoff radii for all shells (based on energies)
  subroutine find_default_cutoffs(species,vha,vxc)

    use datatypes
    use numbers
    use mesh, ONLY: nmesh, rr
    
    implicit none

    ! Passed variables
    integer :: species
    real(double), dimension(nmesh) :: vha,vxc

    ! Local variables
    integer :: i_shell
    real(double), dimension(:), allocatable :: large_cutoff, small_cutoff
    real(double) :: energy, average_large, average_small
    real(double), dimension(:), allocatable :: psi

    ! Find the cutoffs for the valence shells first
    !allocate(psi(nmesh))
    !psi = zero
    allocate(large_cutoff(val(species)%n_occ),small_cutoff(val(species)%n_occ))
    large_cutoff = zero
    small_cutoff = zero
    ! Loop over valence states, find large/small cutoffs
    do i_shell = 1, val(species)%n_occ !paos(species)%n_shells-1 
       if(iprint>3) write(*,*) '# Finding radius for ',paos(species)%npao(i_shell), paos(species)%l(i_shell)
       call find_radius_from_energy(species,paos(species)%npao(i_shell), paos(species)%l(i_shell), &
            large_cutoff(i_shell), val(species)%en_ps(i_shell)+deltaE_large_radius(species), vha, vxc, .false.)
       if(val(species)%semicore(i_shell)==0) then
          call find_radius_from_energy(species,paos(species)%npao(i_shell), paos(species)%l(i_shell), &
               small_cutoff(i_shell), val(species)%en_ps(i_shell)+deltaE_small_radius(species), vha, vxc, .false.)
          !write(*,*) '# Cutoffs: ',large_cutoff, small_cutoff
       else
          !write(*,*) '# Cutoff: ',large_cutoff
          small_cutoff(i_shell) = large_cutoff(i_shell)
       end if
    end do
    ! Create cutoffs based on defaults chosen by user
    if(paos(species)%flag_cutoff==pao_cutoff_energies.OR.paos(species)%flag_cutoff==pao_cutoff_default) then ! Same energy for all l/n shells
       do i_shell = 1, val(species)%n_occ !paos(species)%n_shells-1
          ! NB Semi-core orbitals have large = small
          paos(species)%cutoff(:,i_shell) = zero
          paos(species)%cutoff(1,i_shell) = large_cutoff(i_shell)
          if(paos(species)%nzeta(i_shell)==2) then
             paos(species)%cutoff(2,i_shell) = small_cutoff(i_shell)
          else if(paos(species)%nzeta(i_shell)==3) then
             paos(species)%cutoff(2,i_shell) = half*(large_cutoff(i_shell) + small_cutoff(i_shell))
             paos(species)%cutoff(3,i_shell) = small_cutoff(i_shell)
          end if
          !write(*,*) '# Cutoffs: ',paos(species)%cutoff(:,i_shell)
       end do
       ! Set polarisation radii
       if(paos(species)%n_shells>val(species)%n_occ) then ! Polarisation
          paos(species)%cutoff(1,paos(species)%n_shells)=paos(species)%cutoff(1,paos(species)%n_shells-1)
          if(paos(species)%nzeta(i_shell)==2) then
             paos(species)%cutoff(2,paos(species)%n_shells)=paos(species)%cutoff(2,paos(species)%n_shells-1)
          else if(paos(species)%nzeta(i_shell)==3) then
             paos(species)%cutoff(2,paos(species)%n_shells)=paos(species)%cutoff(2,paos(species)%n_shells-1)
             paos(species)%cutoff(3,paos(species)%n_shells)=paos(species)%cutoff(3,paos(species)%n_shells-1)
          end if
          !write(*,*) '# Cutoffs: ',paos(species)%cutoff(:,paos(species)%n_shells)
       end if
    else if(paos(species)%flag_cutoff==pao_cutoff_radii) then ! Same radius for all l/n shells
       ! Sum over radii of different shells, excluding polarisation
       average_large = zero
       average_small = zero
       do i_shell = 1, val(species)%n_occ!paos(species)%n_shells-1
          average_large = average_large + large_cutoff(i_shell)
          average_small = average_small + small_cutoff(i_shell)
       end do
       ! Find averages for large and small
       average_large = average_large/real(paos(species)%n_shells-1, double)
       average_small = average_small/real(paos(species)%n_shells-1, double)
       ! Set radii for all shells
       paos(species)%cutoff(:,i_shell) = zero
       paos(species)%cutoff(1,:) = average_large
       do i_shell = 1, paos(species)%n_shells
          if(paos(species)%nzeta(i_shell)==2) then
             paos(species)%cutoff(2,:) = average_small
          else if(paos(species)%nzeta(i_shell)==3) then
             paos(species)%cutoff(3,:) = average_small
             ! Now average large and small
             paos(species)%cutoff(2,:) = half*(average_large + average_small)
          end if
          !write(*,*) '# Cutoffs: ',paos(species)%cutoff(:,i_shell)
       end do
    end if
    !deallocate(psi)
    return
  end subroutine find_default_cutoffs
  
  ! Given an energy shift, integrate outwards to find the corresponding radius
  subroutine find_radius_from_energy(species,en,ell,Rc,energy,vha,vxc,flag_use_semilocal)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: rr, rr_squared, nmesh, alpha, make_mesh, beta, sqrt_rr, drdi, drdi_squared, convert_r_to_i
    use pseudo_tm_info, ONLY: pseudo
    
    implicit none

    integer :: species, ell, en
    real(double) :: energy, Rc
    real(double), dimension(nmesh) :: vha,vxc

    ! Local variables
    real(double) :: g_temp, dy_L, dy_R
    real(double), dimension(:), allocatable :: f, potential, psi
    integer :: classical_tp, i, n_crossings, n_nodes, n_loop, loop, nmax, n_kink, n_nodes_lower, n_nodes_upper
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
       !write(*,*) '# Using semi-local potential'
       do i=1,nmesh
          potential(i) = local_and_vkb(species)%semilocal_potential(i,ell) + vha(i) + vxc(i)
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
          potential(i) = local_and_vkb(species)%local(i) + vha(i) + vxc(i)
          g_temp = (drdi_squared(i)*(two*(energy - potential(i))-l_l_plus_one/rr_squared(i)) - alpha_sq_over_four)/twelve
          f(i) = one + g_temp
       end do
       n_kink = nmax
       call integrate_vkb_outwards(species,n_kink,ell,psi,f,n_crossings,n_nodes) ! We want to integrate to max before final node
       do i=n_kink,nmax-1
          psi(i+1) = ( (twelve - ten*f(i)) * psi(i) - f(i-1)*psi(i-1) )/f(i+1)
          if(psi(i)*psi(i+1)<zero.OR.(abs(psi(i+1))<1e-16_double)) then ! Crossing the x-axis
             !write(*,*) '# Node ! ',rr(i),psi(i),psi(i+1)
             n_crossings = n_crossings + 1
             if(n_crossings>=n_nodes+1) exit
          end if
       end do
    end if
    ! Find radius by integrating outwards
    ! Interpolate between i+1 and i
    if(i>=nmax) i=nmax-1
    Rc = rr(i+1)! - psi(i)*(rr(i+1)-rr(i))/(psi(i+1)-psi(i))
    if(iprint>5) write(*,*) '# Found radius ',Rc
    return
  end subroutine find_radius_from_energy
  
  subroutine find_valence_states_semilocal

    use datatypes
    use numbers
    use species_module, ONLY: n_species
    use pseudo_tm_info, ONLY: pseudo, rad_alloc
    use mesh, ONLY: nmesh, make_mesh, rr, nmesh_reg, rmesh_reg, interpolate, make_mesh_reg, rr_squared
    use global_module, ONLY: flag_pcc_global
    use radial_xc, ONLY: get_vxc
    
    implicit none

    integer :: i_species, ell, en, i, zeta, n_shells, i_shell, grid_size, n_val
    real(double), allocatable, dimension(:) :: psi_reg, x_reg, vha, vxc, vha_reg, psi, atomic_rho, &
         small_cutoff, large_cutoff, cutoffs
    real(double) :: energy, max_cutoff, small_energy, large_energy
    real(double) :: radius_small, radius_med, radius_large, sum_small, sum_large
    real(double), dimension(3) :: scaling
    !logical :: flag_plot = .true.
    
    do i_species = 1,n_species
       nmesh = local_and_vkb(i_species)%ngrid
       call make_mesh(i_species)
       allocate(psi_reg(nmesh_reg))
       allocate(vha(nmesh),vxc(nmesh))
       allocate(psi(nmesh),atomic_rho(nmesh))
       ! We need max cutoff to work out largest PAO radius (for neutral atom)
       max_cutoff = zero
       atomic_rho = zero
       call radial_hartree(nmesh,local_and_vkb(i_species)%charge,vha,i_species)
       if(flag_pcc_global) & 
            local_and_vkb(i_species)%charge = local_and_vkb(i_species)%charge + local_and_vkb(i_species)%pcc
       call get_vxc(nmesh,rr,local_and_vkb(i_species)%charge,val(i_species)%functional,vxc)
       if(flag_pcc_global) & 
            local_and_vkb(i_species)%charge = local_and_vkb(i_species)%charge - local_and_vkb(i_species)%pcc
       do i=1,nmesh
          write(35,*) rr(i),vha(i),vxc(i)
       end do
       write(35,*) '&'
       ! Polarisation will be done separately at the end
       n_shells = val(i_species)%n_occ
       do i_shell = 1,n_shells
          ell = val(i_species)%l(i_shell)
          en = val(i_species)%npao(i_shell)
          energy = val(i_species)%en_ps(i_shell)
          zeta = 1
          allocate(paos(i_species)%psi(zeta,i_shell)%f(nmesh))
          write(*,*) '# Species, n, l, zeta, cutoff: ',i_species, en, ell, zeta, paos(i_species)%cutoff(zeta,i_shell)
          write(*,*) '# Occ for shell: ',val(i_species)%occ(i_shell)
          call find_eigenstate_and_energy(i_species,en,ell,paos(i_species)%cutoff(zeta,i_shell),&
               paos(i_species)%psi(zeta,i_shell)%f,energy,vha,vxc)
          paos(i_species)%energy(zeta,i_shell) = energy
       end do
    end do
  end subroutine find_valence_states_semilocal

  subroutine find_valence_states_vkb

    use datatypes
    use numbers
    use species_module, ONLY: n_species
    use pseudo_tm_info, ONLY: pseudo, rad_alloc
    use mesh, ONLY: nmesh, make_mesh, rr, nmesh_reg, rmesh_reg, interpolate, make_mesh_reg, rr_squared
    use global_module, ONLY: flag_pcc_global
    use radial_xc, ONLY: get_vxc
    
    implicit none

    integer :: i_species, ell, en, i, zeta, n_shells, i_shell, grid_size, n_val
    real(double), allocatable, dimension(:) :: psi_reg, x_reg, vha, vxc, vha_reg, psi, atomic_rho, &
         small_cutoff, large_cutoff, cutoffs
    real(double) :: energy, max_cutoff, small_energy, large_energy
    real(double) :: radius_small, radius_med, radius_large, sum_small, sum_large
    real(double), dimension(3) :: scaling
    !logical :: flag_plot = .true.
    
    do i_species = 1,n_species
       nmesh = local_and_vkb(i_species)%ngrid
       call make_mesh(i_species)
       allocate(psi_reg(nmesh_reg))
       allocate(vha(nmesh),vxc(nmesh))
       allocate(psi(nmesh),atomic_rho(nmesh))
       ! We need max cutoff to work out largest PAO radius (for neutral atom)
       max_cutoff = zero
       atomic_rho = zero
       call radial_hartree(nmesh,local_and_vkb(i_species)%charge,vha,i_species)
       if(flag_pcc_global) & 
            local_and_vkb(i_species)%charge = local_and_vkb(i_species)%charge + local_and_vkb(i_species)%pcc
       call get_vxc(nmesh,rr,local_and_vkb(i_species)%charge,val(i_species)%functional,vxc)
       do i=1,nmesh
          write(35,*) rr(i),vha(i),vxc(i)
       end do
       ! Polarisation will be done separately at the end
       n_shells = val(i_species)%n_occ
       do i_shell = 1,n_shells
          ell = val(i_species)%l(i_shell)
          en = val(i_species)%npao(i_shell)
          zeta = 1
          allocate(paos(i_species)%psi(zeta,i_shell)%f(nmesh))
          write(*,*) '# Species, n, l, zeta, cutoff: ',i_species, en, ell, zeta, paos(i_species)%cutoff(zeta,i_shell)
          write(*,*) '# Occ for shell: ',val(i_species)%occ(i_shell)
          call find_eigenstate_and_radius_vkb(i_species,en,ell,paos(i_species)%cutoff(zeta,i_shell),&
               paos(i_species)%psi(zeta,i_shell)%f,paos(i_species)%energy(zeta,i_shell),vha,vxc,.false.)
       end do
    end do
  end subroutine find_valence_states_vkb
  
  ! Rc has been specified - find the energy that gives this cutoff
  ! VKB version
  ! Use local potential for homogeneous equation, and solve inhomogeneous
  ! equations for projectors (then combine solutions)
  subroutine find_eigenstate_and_energy_vkb(species,en,ell,Rc,psi,energy,vha,vxc)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: rr, rr_squared, nmesh, alpha, make_mesh, beta, sqrt_rr, drdi, drdi_squared, convert_r_to_i
    use pseudo_tm_info, ONLY: pseudo
    
    implicit none

    integer :: species, ell, en
    real(double) :: energy, Rc
    real(double), dimension(nmesh) :: psi,vha,vxc

    ! Local variables
    real(double) :: g_temp, dy_L, dy_R
    real(double), dimension(:), allocatable :: f, potential
    integer :: classical_tp, i, n_crossings, n_nodes, n_loop, loop, nmax, n_kink, n_nodes_lower, n_nodes_upper, n_kink_vkb
    real(double) :: l_half_sq, dx_sq_over_twelve, fac, norm, d_energy, e_lower, e_upper, df_cusp, cusp_psi, tol
    real(double) :: delta_energy_bracket, zval, l_l_plus_one, alpha_sq_over_four, xkap
    real(double) :: gin, gout, gsgin, gsgout, xin, xout
    logical :: flag_find_radius = .false.

    l_half_sq = real(ell,double) + half
    l_half_sq = l_half_sq*l_half_sq
    l_l_plus_one = real(ell*(ell+1),double)
    zval = pseudo(species)%zval ! Added
    !write(*,*) '# zval is ',zval
    dx_sq_over_twelve = alpha*alpha/twelve
    alpha_sq_over_four = alpha*alpha/four
    n_nodes = en - ell - 1 
    allocate(f(nmesh),potential(nmesh))
    if(abs(energy)<RD_ERR) then
       e_lower = -zval!-zval*zval/real(en*en,double)
    else if(energy<zero) then
       e_lower = energy*1.2_double
    else ! Unbound (polarisation) state
       e_lower = zero
    end if
    ! Energy bounds - allow for unbound states
    e_upper = five ! One failed to find the tightest PAO for O
    do i=1,nmesh
       potential(i) = local_and_vkb(species)%local(i) + vha(i) + vxc(i)  ! Half when using Siesta
       !g_temp = l_l_plus_one/(rr_squared(i)) + potential(i)
       !if(g_temp<e_lower) e_lower = g_temp
    end do
    ! Now set number of loops and maximum radius
    n_loop = 100
    call convert_r_to_i(Rc,nmax)
    nmax = nmax - 1
    if(abs(energy)<RD_ERR) then
       energy = half*(e_lower+e_upper)
    end if
    tol = 1.0e-8_double
    !write(*,*) "# Nmax is ",nmax
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
       call integrate_vkb_outwards(species,n_kink,ell,psi,f,n_crossings,n_nodes) ! We want to integrate to max before final node
       n_kink_vkb = n_kink
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
          energy = half*(e_lower+e_upper)
          cycle
       end if
       if(n_crossings /= n_nodes) then
          if ( n_crossings > n_nodes ) then
             e_upper = energy
          else
             e_lower = energy
          end if
          energy = half * ( e_upper + e_lower )
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
       cusp_psi = xin + xout
       cusp_psi = cusp_psi + twelve*(one-f(n_kink))*gout
       cusp_psi = cusp_psi*gout/(gsgin + gsgout)
       d_energy = cusp_psi
       if(iprint>4) write(*,fmt='(2x,"Energy shift (alt): ",f18.10)') cusp_psi
       if(iprint>5) write(*,fmt='(2x,"Number of nodes: ",i4)') n_crossings
       if ( n_crossings /= n_nodes) then
          if ( n_crossings > n_nodes ) then
             e_upper = energy
          else
             e_lower = energy
          end if
          energy = half * ( e_upper + e_lower )
          cycle
       end if       
       if(d_energy>zero) then
          e_lower = energy
       else
          e_upper = energy
       end if
       if(energy+d_energy<e_upper.AND.energy+d_energy>e_lower) then
          energy = energy + d_energy
       else
          energy = half*(e_lower + e_upper)
       end if
       if(abs(d_energy)<tol) exit
    end do
    if(loop>=100.AND.abs(d_energy)>tol) call cq_abort("Error: failed to find energy for n,l: ",en,ell)
    if(iprint>2) write(*,fmt='(2x,"Final energy found: ",f18.10)') energy
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
  
  ! Energy shift has been specified - find the radius that gives this shift
  ! VKB version
  ! Use local potential for homogeneous equation, and solve inhomogeneous
  ! equations for projectors (then combine solutions)
  subroutine find_eigenstate_and_radius_vkb(species,en,ell,Rc,psi,energy,vha,vxc,flag_use_semilocal)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: rr, rr_squared, nmesh, alpha, make_mesh, beta, sqrt_rr, drdi, drdi_squared, convert_r_to_i
    use pseudo_tm_info, ONLY: pseudo
    
    implicit none

    integer :: species, ell, en
    real(double) :: energy, Rc
    real(double), dimension(nmesh) :: psi,vha,vxc
    logical :: flag_use_semilocal

    ! Local variables
    real(double) :: g_temp, dy_L, dy_R
    real(double), dimension(:), allocatable :: f, potential
    integer :: classical_tp, i, n_crossings, n_nodes, n_loop, loop, nmax, n_kink, n_nodes_lower, n_nodes_upper
    real(double) :: l_half_sq, dx_sq_over_twelve, fac, norm, d_energy, e_lower, e_upper, df_cusp, cusp_psi, tol
    real(double) :: delta_energy_bracket, zval, l_l_plus_one, alpha_sq_over_four, xkap
    real(double) :: gin, gout, gsgin, gsgout, xin, xout
    logical :: flag_find_radius = .false.

    nmax = nmesh
    l_l_plus_one = real(ell*(ell+1),double)
    alpha_sq_over_four = alpha*alpha/four
    n_nodes = en - ell - 1 
    allocate(f(nmesh),potential(nmesh))
    tol = 1.0e-8_double
    !write(*,*) "# Required number of nodes is ",n_nodes
    n_crossings = 0
    psi = zero
    if(flag_use_semilocal) then
       !write(*,*) '# Using semi-local potential'
       ! Find kink roughly
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
             if(n_crossings>=n_nodes) exit
          end if
       end do
       ! Interpolate between i+1 and i
       if(i>=nmax) i=nmax-1
       Rc = rr(i+1)! - psi(i)*(rr(i+1)-rr(i))/(psi(i+1)-psi(i))
       !write(*,*) '# Found cutoff ',Rc
       nmax = i+1
       psi = zero
       n_kink = nmax
       xkap = two*local_and_vkb(species)%semilocal_potential(1,ell)/(rr(1)**real(ell+1,double)*(six + four*real(ell)))
       psi(1) = rr(1)**real(ell+1,double)/sqrt(drdi(1))
       psi(2) = rr(2)**real(ell+1,double)/sqrt(drdi(2))
       call numerov(1,n_kink,nmesh,psi,rr,f,1,n_crossings,n_nodes,xkap,0) ! Integrate out
    else
       !write(*,*) '# Using VKB potentials'
       do i=1,nmesh
          potential(i) = local_and_vkb(species)%local(i) + vha(i) + vxc(i)
          g_temp = (drdi_squared(i)*(two*(energy - potential(i))-l_l_plus_one/rr_squared(i)) - alpha_sq_over_four)/twelve
          f(i) = one + g_temp
       end do
       n_kink = nmax
       call integrate_vkb_outwards(species,n_kink,ell,psi,f,n_crossings,n_nodes) ! We want to integrate to max before final node
       !write(*,*) '# Kink, crossings: ',rr(n_kink),n_crossings
       do i=n_kink,nmax-1
          psi(i+1) = ( (twelve - ten*f(i)) * psi(i) - f(i-1)*psi(i-1) )/f(i+1)
          if(psi(i)*psi(i+1)<zero.OR.(abs(psi(i+1))<1e-16_double)) then ! Crossing the x-axis
             !write(*,*) '# Node ! ',rr(i),psi(i),psi(i+1)
             n_crossings = n_crossings + 1
             if(n_crossings>=n_nodes+1) exit
          end if
       end do
       ! Interpolate between i+1 and i
       if(i>=nmax) i=nmax-1
       Rc = rr(i+1)! - psi(i)*(rr(i+1)-rr(i))/(psi(i+1)-psi(i))
       !write(*,*) '# Found cutoff ',Rc
       nmax = i+1
       psi = zero
       n_kink = nmax
       call integrate_vkb_outwards(species,n_kink,ell,psi,f,n_crossings,n_nodes) ! We want to integrate to max before final node
    end if
    fac = psi(n_kink)
    psi(nmax) = zero!exp(-xkap)
    psi(nmax-1) = one !exp(-xkap*(rr(nmax) - rr(nmax-1)))
    call numerov(n_kink+1,nmax-n_kink,nmax,psi,rr,f,-1,n_crossings,n_nodes,xkap,0)
    fac = fac/psi(n_kink)
    psi(n_kink:nmax) = psi(n_kink:nmax)*fac
    ! Normalise
    norm = zero
    !call integrate_simpson(psi*psi,drdi_squared,nmax,norm)
    do i=1,nmax
       norm = norm + psi(i)*psi(i)*drdi_squared(i)
    end do
    psi = psi/sqrt(norm)
    ! Rescale - remove factor of sqrt r
    do i=1,nmax
       psi(i) = psi(i)*sqrt(drdi(i))/rr(i)
       !write(13,*) rr(i),(two*rr_squared(i)*potential(i)+l_half_sq)*alpha*alpha, &
       !     alpha*alpha*rr_squared(i)
    end do
    if(psi(nmax - 5)<zero) psi = -psi
    if(flag_run_debug) then
       write(25+ell,*) '# VKB find radius'
       do i=1,nmesh
          write(25+ell,*) rr(i),psi(i)
       end do
       write(25+ell,*) '&'
    end if
    if(ell>0) then
       do loop = 1,ell
          psi = psi/rr
       end do
    end if
    deallocate(f,potential)
  end subroutine find_eigenstate_and_radius_vkb
  
  ! Perform outward integration when using Vanderbilt-Kleinman-Bylander projectors
  ! We solve for the homogeneous equation (local potential) and then nproj inhomogeneous
  ! equations (one for each projector) which are then combined
  ! This serves as a layer between the eigenstate finder and numerov
  subroutine integrate_vkb_outwards(species,n_kink,ell,psi,f,n_crossings,n_nodes_inside)

    use datatypes
    use numbers
    use mesh, ONLY: rr, nmesh, drdi, drdi_squared
    use pseudo_tm_info, ONLY: pseudo
    use GenComms, ONLY: cq_abort
    
    implicit none

    ! Passed variables
    integer :: species, n_kink, ell, n_crossings, n_nodes_inside
    real(double), dimension(nmesh) :: psi, f
    
    ! Local variables
    integer :: i_proj, j_proj, i, info, nproj_acc, j
    integer, dimension(4) :: ipiv ! Larger than needed
    real(double) :: integral, xkap, r_l_p_1, dy, dy_last
    real(double), allocatable, dimension(:) :: pot_vector, psi_h, s, integrand
    real(double), allocatable, dimension(:,:) :: pot_matrix, psi_inh

    allocate(pot_matrix(n_proj(ell),n_proj(ell)),pot_vector(n_proj(ell)))
    pot_matrix = zero
    pot_vector = zero
    allocate(psi_h(nmesh),psi_inh(nmesh,n_proj(ell)))
    psi_h = zero
    psi_inh = zero
    if(ell == 0) then
       nproj_acc = 0
    else
       nproj_acc = sum(n_proj(0:ell-1))
    end if
    ! Homogeneous - standard numerov - but only to projector radius
    psi_h = zero
    n_kink = local_and_vkb(species)%ngrid_vkb
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
    do i_proj = 1,n_proj(ell)
       do j=1,n_kink
          integrand(j) = psi_h(j)*sqrt(drdi(j))* &
               local_and_vkb(species)%projector(j,i_proj,ell)
       end do
       call integrate_simpson(integrand,drdi(1:n_kink),n_kink,integral)
       if(iprint>6) write(*,fmt='("In integrate_vkb, homogeneous integral is ",f18.10)') integral
       pot_vector(i_proj) = integral*pseudo(species)%pjnl_ekb(nproj_acc + i_proj)
    end do
    ! Inhomogeneous - add VKB projector as source term
    do i_proj = 1,n_proj(ell)
       s(1:n_kink) = two*local_and_vkb(species)%projector(1:n_kink,i_proj,ell)*drdi(1:n_kink)*sqrt(drdi(1:n_kink))
       ! Series expansion based on projector (FIND APPROPRIATE REFERENCE !)
       xkap = two*local_and_vkb(species)%projector(1,i_proj,ell)/(r_l_p_1*(six + four*real(ell,double)))
       psi_inh(1,i_proj) = xkap*rr(1)**real(ell+3,double)/sqrt(drdi(1))
       psi_inh(2,i_proj) = xkap*rr(2)**real(ell+3,double)/sqrt(drdi(2))
       call numerov(1,n_kink,nmesh,psi_inh(:,i_proj),rr,f,1,n_crossings,n_nodes_inside, &
            xkap,n_kink,s) 
       ! Integrals between this solution and VKB projectors
       do j_proj = 1,n_proj(ell)
          call integrate_simpson(psi_inh(1:n_kink,i_proj)*(sqrt(drdi(1:n_kink)))* & ! /rr(1:n_kink)
               local_and_vkb(species)%projector(1:n_kink,j_proj,ell), drdi,n_kink,integral)
          if(iprint>6) write(*,fmt='("In integrate_vkb, inhomogeneous integral is ",f18.10)') integral
          pot_matrix(j_proj,i_proj) = -integral*pseudo(species)%pjnl_ekb(nproj_acc + j_proj)
       end do
       pot_matrix(i_proj,i_proj) = one + pot_matrix(i_proj,i_proj)
    end do
    !write(*,*) '# Pot mat and vec: ',pot_matrix, pot_vector
    ! Invert matrix
    call dgesv(n_proj(ell), 1, pot_matrix, n_proj(ell), ipiv, pot_vector, n_proj(ell), info)
    if(info/=0) call cq_abort("Error from dgesv called in integrate_vkb_outward: ",info)
    if(iprint>6) write(*,fmt='("In integrate_vkb, coefficients are ",3f18.10)') pot_vector
    ! Construct total outward wavefunction
    psi(1:n_kink) = psi_h(1:n_kink)
    do i_proj=1,n_proj(ell)
       psi(1:n_kink) = psi(1:n_kink) + pot_vector(i_proj)*psi_inh(1:n_kink,i_proj)
    end do
    ! Search for nodes
    n_crossings = 0
    dy = psi(2) - psi(1)
    do i=2,n_kink
       if(rr(i)>half.AND.psi(i)*psi(i-1)<zero) n_crossings = n_crossings+1
    end do
    deallocate(pot_matrix, pot_vector, s, psi_h, psi_inh)
  end subroutine integrate_vkb_outwards
  
  ! Energy has been specified - find the cutoff that gives this energy
  subroutine find_eigenstate_and_radius(species,en,ell,Rc,psi,energy,vha,vxc)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: rr, rr_squared, nmesh, alpha, make_mesh, beta, sqrt_rr, drdi, drdi_squared, convert_r_to_i
    use pseudo_tm_info, ONLY: pseudo
    
    implicit none

    integer :: species, ell, en
    real(double) :: energy, Rc
    real(double), dimension(nmesh) :: psi,vha,vxc

    ! Local variables
    real(double) :: g_temp, dy_L, dy_R
    real(double), dimension(:), allocatable :: f, potential
    integer :: classical_tp, i, n_crossings, n_nodes, n_loop, loop, nmax, n_kink, n_nodes_lower, n_nodes_upper
    real(double) :: l_half_sq, dx_sq_over_twelve, fac, norm, d_energy, r_lower, r_upper, df_cusp, cusp_psi, tol
    real(double) :: delta_energy_bracket, zval, l_l_plus_one, alpha_sq_over_four, xkap
    real(double) :: gin, gout, gsgin, gsgout, xin, xout
    logical :: flag_find_radius = .false.

    if(iprint>6) write(*,fmt='("SL Finding cutoff radius that gives energy ",f8.3)') energy
    l_half_sq = real(ell,double) + half
    l_half_sq = l_half_sq*l_half_sq
    l_l_plus_one = real(ell*(ell+1),double)
    zval = pseudo(species)%zval ! Added
    dx_sq_over_twelve = alpha*alpha/twelve
    alpha_sq_over_four = alpha*alpha/four
    n_nodes = en - ell - 1 
    allocate(f(nmesh),potential(nmesh))
    do i=1,nmesh
       potential(i) = local_and_vkb(species)%semilocal_potential(i,ell) + vha(i) + vxc(i)
    end do
    ! Find radius by integrating outwards
    nmax = nmesh
    n_crossings = 0
    do i=1,nmax
       g_temp = (drdi_squared(i)*(two*(energy - potential(i))-l_l_plus_one/rr_squared(i)) - alpha_sq_over_four)/twelve!*dx_sq_over_twelve
       f(i) = one + g_temp
    end do
    psi = zero
    psi(1) = one
    psi(2) = ( (twelve - ten*f(1)) * psi(1) )/f(2) ! Implicit psi(0) = zero
    dy_R = psi(2) - psi(1)
    do i=2,nmax-1
       dy_L = dy_R
       psi(i+1) = ( (twelve - ten*f(i)) * psi(i) - f(i-1)*psi(i-1) )/f(i+1)
       dy_R = psi(i+1) - psi(i)
       if(psi(i)*psi(i+1)<zero.OR.(abs(psi(i+1))<1e-16_double)) then ! Crossing the x-axis
          n_crossings = n_crossings + 1
          if(n_crossings>=n_nodes) exit
       end if
    end do
    ! Interpolate between i+1 and i
    if(i>=nmax) i=nmax-1
    Rc = rr(i+1)! - psi(i)*(rr(i+1)-rr(i))/(psi(i+1)-psi(i))
    nmax = i+1
    psi = zero
    n_kink = nmax
    xkap = two*local_and_vkb(species)%semilocal_potential(1,ell)/(rr(1)**real(ell+1,double)*(six + four*real(ell)))
    psi(1) = rr(1)**real(ell+1,double)/sqrt(drdi(1))
    psi(2) = rr(2)**real(ell+1,double)/sqrt(drdi(2))
    call numerov(1,n_kink,nmesh,psi,rr,f,1,n_crossings,n_nodes,xkap,0) ! We want to integrate to max before final node
    fac = psi(n_kink)
    xkap = sqrt(two*(potential(i) - energy)+l_l_plus_one/rr_squared(i))
    psi(nmax) = zero!exp(-xkap)
    psi(nmax-1) = one !exp(-xkap*(rr(nmax) - rr(nmax-1)))
    call numerov(n_kink+1,nmax-n_kink,nmax,psi,rr,f,-1,n_crossings,n_nodes,xkap,0)
    fac = fac/psi(n_kink)
    psi(n_kink:nmax) = psi(n_kink:nmax)*fac
    ! Normalise
    norm = zero
    !call integrate_simpson(psi*psi,drdi_squared,nmax,norm)
    do i=1,nmax
       norm = norm + psi(i)*psi(i)*drdi_squared(i)
    end do
    psi = psi/sqrt(norm)
    ! Rescale - remove factor of sqrt r
    do i=1,nmax
       psi(i) = psi(i)*sqrt(drdi(i))/rr(i)
       !write(13,*) rr(i),(two*rr_squared(i)*potential(i)+l_half_sq)*alpha*alpha, &
       !     alpha*alpha*rr_squared(i)
    end do
    if(psi(nmax - 5)<zero) psi = -psi
    if(ell>0) then
       do loop = 1,ell
          psi = psi/rr
       end do
    end if
    deallocate(f,potential)
    return
  end subroutine find_eigenstate_and_radius

  ! Rc has been specified - find the energy that gives this cutoff
  subroutine find_eigenstate_and_energy(species,en,ell,Rc,psi,energy,vha,vxc)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: rr, rr_squared, nmesh, alpha, make_mesh, beta, sqrt_rr, drdi, drdi_squared, convert_r_to_i
    use pseudo_tm_info, ONLY: pseudo
    
    implicit none

    integer :: species, ell, en
    real(double) :: energy, Rc
    real(double), dimension(nmesh) :: psi,vha,vxc

    ! Local variables
    real(double) :: g_temp, dy_L, dy_R
    real(double), dimension(:), allocatable :: f, potential
    integer :: classical_tp, i, n_crossings, n_nodes, n_loop, loop, nmax, n_kink, n_nodes_lower, n_nodes_upper
    real(double) :: l_half_sq, dx_sq_over_twelve, fac, norm, d_energy, e_lower, e_upper, df_cusp, cusp_psi, tol
    real(double) :: delta_energy_bracket, zval, l_l_plus_one, alpha_sq_over_four, xkap
    real(double) :: gin, gout, gsgin, gsgout, xin, xout
    logical :: flag_find_radius = .false.

    l_half_sq = real(ell,double) + half
    l_half_sq = l_half_sq*l_half_sq
    l_l_plus_one = real(ell*(ell+1),double)
    zval = pseudo(species)%zval ! Added
    write(*,*) '# zval is ',zval
    dx_sq_over_twelve = alpha*alpha/twelve
    alpha_sq_over_four = alpha*alpha/four
    n_nodes = en - ell - 1 
    allocate(f(nmesh),potential(nmesh))
    if(abs(energy)<RD_ERR) then
       e_lower = -zval!-zval*zval/real(en*en,double)
    else
       e_lower = energy*1.2_double
    end if
    ! Energy bounds - allow for unbound states
    e_upper = one
    do i=1,nmesh
       potential(i) = local_and_vkb(species)%semilocal_potential(i,ell) + vha(i) + vxc(i)
       !g_temp = l_l_plus_one/(rr_squared(i)) + potential(i)
       !if(g_temp<e_lower) e_lower = g_temp
    end do
    ! Now set number of loops and maximum radius
    n_loop = 100
    call convert_r_to_i(Rc,nmax)
    nmax = nmax - 1
    if(abs(energy)<RD_ERR) then
       energy = half*(e_lower+e_upper)
    end if
    tol = 1.0e-8_double
    write(*,*) "# Nmax is ",nmax
    write(*,*) "# Required number of nodes is ",n_nodes
    ! Now loop to find the correct energy
    do loop = 1,n_loop
       write(*,*) "# loop ",loop
       write(*,*) '# Energies: ',e_lower,energy,e_upper
       n_kink = nmax
       do i=1,nmax
          g_temp = (drdi_squared(i)*(two*(energy - potential(i))-l_l_plus_one/rr_squared(i)) - alpha_sq_over_four)/twelve!*dx_sq_over_twelve
          !if(i>1) then
          !   if(g_temp*(f(i-1)-one)<zero) n_kink = i
          !end if
          f(i) = one + g_temp
       end do
       !write(*,*) '# Kink at ',n_kink,rr(n_kink)
       psi = zero
       ! Outward
       n_crossings = 0
       xkap = two*local_and_vkb(species)%semilocal_potential(1,ell)/(rr(1)**real(ell+1,double)*(six + four*real(ell)))
       psi(1) = rr(1)**real(ell+1,double)/sqrt(drdi(1))
       psi(2) = rr(2)**real(ell+1,double)/sqrt(drdi(2))
       call numerov(1,n_kink,nmax,psi,rr,f,1,n_crossings,n_nodes,xkap,0) ! We want to integrate to max before final node
       write(*,*) '# Kink at ',n_kink,rr(n_kink)
       xout = psi(n_kink)*f(n_kink) - psi(n_kink-1)*f(n_kink-1)
       gout = psi(n_kink)
       gsgout = psi(n_kink)*psi(n_kink)*drdi_squared(n_kink)
       if(n_kink == nmax) then
          write(*,*) '#No kink'
          e_lower = energy
          energy = half*(e_lower+e_upper)
          cycle
       end if
       if(n_crossings /= n_nodes) then
          if ( n_crossings > n_nodes ) then
             e_upper = energy
          else
             e_lower = energy
          end if
          energy = half * ( e_upper + e_lower )
          write(*,*) '# Nodes wrong: ',n_crossings, n_nodes
          !write(8,*) '# Nodes wrong: ',n_crossings, n_nodes
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
       write(*,*) '#psi at ctp: ',fac,psi(n_kink)
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
       !if(en==3.AND.ell==2) then
       !   do i=1,nmax
       !      write(12,*) rr(i),psi(i)*sqrt(drdi(i))/rr(i)
       !   end do
       !   write(12,*) '&'
       !end if
       ! Cusp condition to adjust 
       ! [yL(i-1) + yR(i+1) - (14-12f(i))y(i)]/dx > 0 means energy too high
       dy_L = (psi(n_kink-1)*sqrt(drdi(n_kink-1))/rr(n_kink-1)-psi(n_kink-2)*sqrt(drdi(n_kink-2))/rr(n_kink-2)) &
            /(rr(n_kink-1)-rr(n_kink-2))
       dy_R = (psi(n_kink+2)*sqrt(drdi(n_kink+2))/rr(n_kink+2)-psi(n_kink+1)*sqrt(drdi(n_kink+1))/rr(n_kink+1)) &
            /(rr(n_kink+2)-rr(n_kink+1))
       write(*,*) '# dy_L: ',dy_L
       write(*,*) '# dy_L alt: ',(psi(n_kink-1)-psi(n_kink-2))/(rr(n_kink-1)-rr(n_kink-2))
       write(*,*) '# dy_R: ',dy_R
       write(*,*) '# dy_R alt: ',(psi(n_kink+2)-psi(n_kink+1))/(rr(n_kink+2)-rr(n_kink+1))
       !write(*,*) '# Energy indicator: ',(psi(i-1) + psi(i+1) - psi(i)*(14.0_double-twelve*f(i)))/alpha
       !d_energy = (psi(i-1) + psi(i+1) - psi(i)*(14.0_double-twelve*f(i)))/alpha
       d_energy = dy_R - dy_L
       !if(dy_R>zero) &
       !     d_energy = dy_R
       write(*,*) '# Energy indicator: ',d_energy
       !dy_L = (psi(n_kink-1)-psi(n_kink-2))/(rr(n_kink-1)-rr(n_kink-2))
       !dy_R = (psi(n_kink+2)-psi(n_kink+1))/(rr(n_kink+2)-rr(n_kink+1))
       !d_energy = dy_R - dy_L
       i = n_kink
       ! Perturbation theory
       cusp_psi = ( psi(i-1)*f(i-1) + f(i+1)*psi(i+1) + ten*f(i)*psi(i)) / twelve
       df_cusp = f(i)*( psi(i)/cusp_psi - one )
       !%%!!d_energy = half*df_cusp/dx_sq_over_twelve * cusp_psi * cusp_psi * alpha
       !write(*,*)  '# Alt dE: ',half*df_cusp*twelve * cusp_psi * cusp_psi * (drdi(i)/rr(i))*(drdi(i)/rr(i))
       !write(*,*)  '# Alt dE: ',(two*psi(n_kink)*f(n_kink) - psi(n_kink-1)*f(n_kink-1) - psi(n_kink+1)*f(n_kink+1) + &
       !     twelve*(one - f(n_kink))*psi(n_kink))*psi(n_kink)/(psi(n_kink)*psi(n_kink)*drdi_squared(n_kink))
       !write(*,*) '# Alt dE(2): ',(df_cusp*twelve * cusp_psi * cusp_psi)/drdi(i) !* (drdi(i)/rr(i))!*(drdi(i)/rr(i))
       d_energy = (df_cusp*twelve * cusp_psi * cusp_psi)/drdi(i)
       write(*,*) '# Energies: ',e_lower,energy,e_upper
       write(*,*) '# dE: ',d_energy
       cusp_psi = xin + xout!two*psi(n_kink)*f(n_kink) - psi(n_kink-1)*f(n_kink-1) - psi(n_kink+1)*f(n_kink+1) ! XIN + XOUT
       ! ARW lines
       cusp_psi = cusp_psi + twelve*(one-f(n_kink))*gout!psi(n_kink)
       cusp_psi = cusp_psi*gout/(gsgin + gsgout)!psi(n_kink)/(psi(n_kink)*psi(n_kink)*drdi_squared(n_kink))
       ! Hamann lines
       !cusp_psi = half*cusp_psi*psi(n_kink)/drdi(n_kink)
       d_energy = cusp_psi
       write(*,*) '# Alt dE: ',cusp_psi
       if ( n_crossings == n_nodes) then
          write(*,*) '# Crossings: ',n_crossings
       else
          write(*,*) '# Crossings: ',n_crossings
          if ( n_crossings > n_nodes ) then
             e_upper = energy
          else
             e_lower = energy
          end if
          energy = half * ( e_upper + e_lower )
          write(*,*) '# Nodes wrong: ',n_crossings, n_nodes
          !write(8,*) '# Nodes wrong: ',n_crossings, n_nodes
          cycle
       end if       
       !%%!if(psi(n_kink)>0) then ! Gradients will be negative
       !%%!   if(dy_L<dy_R) then ! Deconfine
       !%%!      e_upper = energy
       !%%!   else
       !%%!      e_lower = energy
       !%%!   end if
       !%%!else
       !%%!   if(dy_L<dy_R) then ! Confine
       !%%!      e_lower = energy
       !%%!   else
       !%%!      e_upper = energy
       !%%!   end if
       !%%!end if
       if(d_energy>zero) then
          !%%! if(d_energy>(e_upper - energy)) d_energy = half*(e_upper-energy)
          e_lower = energy
          !e_upper = energy
       else
          !%%! if(d_energy<(e_lower - energy)) d_energy = half*(e_lower-energy)
          e_upper = energy
          !e_lower = energy
       end if
       !energy = half*(e_lower+e_upper)
       if(energy+d_energy<e_upper.AND.energy+d_energy>e_lower) then
          energy = energy + d_energy
       else
          energy = half*(e_lower + e_upper)
       end if
       if(abs(d_energy)<tol) exit
    end do
    write(*,*) '# Final Energy: ',energy
    ! Rescale - remove factor of sqrt r
    do i=1,nmax
       psi(i) = psi(i)*sqrt(drdi(i))/rr(i)
       !write(13,*) rr(i),(two*rr_squared(i)*potential(i)+l_half_sq)*alpha*alpha, &
       !     alpha*alpha*rr_squared(i)
    end do
    if(psi(nmax - 5)<zero) psi = -psi
    write(25+ell,*) '# Semilocal find energy'
    do i=1,nmesh
       write(25+ell,*) rr(i),psi(i)
    end do
    write(25+ell,*) '&'
    if(ell>0) then
       do loop = 1,ell
          psi = psi/rr
       end do
    end if
    deallocate(f,potential)
  end subroutine find_eigenstate_and_energy
  
  subroutine find_polarisation(species,en,ell,Rc,psi_l,psi_pol,energy,vha,vxc)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: rr, rr_squared, nmesh, alpha, make_mesh, beta, sqrt_rr, drdi, drdi_squared, convert_r_to_i
    use pseudo_tm_info, ONLY: pseudo
    
    implicit none

    integer :: species, ell, en
    real(double) :: energy, Rc
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
    pf_low = one
    pf_high = 10000.0_double
    prefac = 1000.0_double ! Arbitrary initial value
    l_half_sq = real(ell+1,double) + half
    l_half_sq = l_half_sq*l_half_sq
    l_l_plus_one = real((ell+2)*(ell+1),double)
    n_nodes = en - ell - 1 
    zval = pseudo(species)%z
    dx_sq_over_twelve = alpha*alpha/twelve
    alpha_sq_over_four = alpha*alpha/four
    allocate(f(nmesh),g(nmesh), s(nmesh), potential(nmesh))
    !call radial_hartree(nmesh,semilocal(species)%charge,vha,species)
    !call vxc_pz_ca(nmesh,rr,semilocal(species)%charge,vxc)
    e_lower = zero!-zval*zval/real(en*en,double)
    do i=1,nmesh
       potential(i) = local_and_vkb(species)%semilocal_potential(i,ell) + vha(i) + vxc(i)
    end do
    n_loop = 100
    tol = 1.0e-8_double
    ! Energy bounds - allow for unbound states
    nmax = nmesh ! Adjust later to confine
    ! Test
    !Rc = 8.0_double
    call convert_r_to_i(Rc,nmax)
    nmax = nmax - 1
    !write(*,*) "# Nmax is ",nmax
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
       psi_pol(1) = prefac*rr(1)**(ell+2)!*(one - zval*r(1)/(real(ell+1,double)))
       psi_pol(1) = psi_pol(1)/sqrt(drdi(1))
       psi_pol(2) = prefac*rr(2)**(ell+2)!*(one - zval*r(2)/(real(ell+1,double)))
       psi_pol(2) = psi_pol(2)/sqrt(drdi(2))
       ! We need crossings counted - if there are crossings, then increase, else decrease
       n_crossings = 0
       do i=2,nmax-1
          psi_pol(i+1) = (two*psi_pol(i)*(one-five*g(i)) - psi_pol(i-1)*f(i-1) + (s(i+1)+ten*s(i)+s(i-1))/twelve)/f(i+1)
          if(psi_pol(i+1)*psi_pol(i)<zero) n_crossings = n_crossings + 1
       end do
       if(n_crossings > n_nodes) then ! Increase
          pf_low = prefac
          prefac = half*(pf_low + pf_high)
          if(iprint>5) write(*,fmt='("More crossings than nodes required: ",2i3)') n_crossings, n_nodes
          if(iprint>5) write(*,fmt='("Polarisation loop ",i3," brackets and prefactor: ",3f13.5)')loop,pf_low,prefac,pf_high
          !do i=1,nmax
          !   write(13,*) rr(i),psi_pol(i)
          !end do
          !write(13,*) '&'
          cycle
       else
          pf_high = prefac
          prefac = half*(pf_low+pf_high)
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
    deallocate(f,potential)
  end subroutine find_polarisation

  subroutine set_pao_specification

    use species_module, ONLY: n_species
    use GenComms, ONLY: cq_abort
    
    implicit none

    integer :: i_species, i_shell

    if(flag_user_specified) then
       ! Check compatibility with valence shells
       do i_species = 1,n_species
          paos(i_species)%lmax = 0
          if(paos(i_species)%flag_perturb_polarise) then
             !if(paos(i_species)%n_shells/=val(i_species)%n_occ+1) & 
             !     call cq_abort("Incompatible PSP and PAO valence shells: ",paos(i_species)%n_shells, n_valence_shells(i_species))
             do i_shell = 1,val(i_species)%n_occ 
                if(paos(i_species)%l(i_shell)/=val(i_species)%l(i_shell)) &
                     call cq_abort("Incompatible l values for shell: ",paos(i_species)%l(i_shell),val(i_species)%l(i_shell))
                if(paos(i_species)%l(i_shell)>paos(i_species)%lmax) paos(i_species)%lmax = paos(i_species)%l(i_shell)
             end do
          else
             !if(paos(i_species)%n_shells/=val(i_species)%n_occ) &
             !     call cq_abort("Incompatible PSP and PAO valence shells: ",paos(i_species)%n_shells, n_valence_shells(i_species))
             do i_shell = 1,val(i_species)%n_occ
                if(paos(i_species)%l(i_shell)/=val(i_species)%l(i_shell)) &
                     call cq_abort("Incompatible l values for shell: ",paos(i_species)%l(i_shell),val(i_species)%l(i_shell))
                if(paos(i_species)%l(i_shell)>paos(i_species)%lmax) paos(i_species)%lmax = paos(i_species)%l(i_shell)
             end do
          end if
       end do
    else
       do i_species = 1,n_species
          paos(i_species)%lmax = 0
          allocate(paos(i_species)%nzeta(paos(i_species)%n_shells),paos(i_species)%l(paos(i_species)%n_shells), &
               paos(i_species)%n(paos(i_species)%n_shells))
          do i_shell = 1,paos(i_species)%n_shells
             paos(i_species)%l(i_shell)=val(i_species)%l(i_shell)
             paos(i_species)%n(i_shell)=val(i_species)%n(i_shell)
             if(val(i_species)%semicore(i_shell)==1) then ! Semi-core
                paos(i_species)%nzeta(i_shell)=1
             else
                paos(i_species)%nzeta(i_shell)=3
             end if
             if(paos(i_species)%l(i_shell)>paos(i_species)%lmax) paos(i_species)%lmax = paos(i_species)%l(i_shell)
          end do
       end do
    end if
  end subroutine set_pao_specification

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
  subroutine radial_hartree(n_tot,rho,vha,species)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort
    use mesh, ONLY: alpha, drdi, rr, rr_squared, sqrt_rr, sqrt_drdi
    
    implicit none

    ! Passed variables
    integer :: n_tot, species
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
    deallocate(y)
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
