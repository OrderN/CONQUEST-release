! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module energy
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/energy
!!  NAME
!!   energy
!!  PURPOSE
!!   Stores various energy variables and routines
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   17:17, 2003/03/10
!!  MODIFICATION HISTORY
!!   13:21, 22/09/2003 drb
!!   - Added band_energy
!!   10:09, 13/02/2006 drb
!!   - Removed all explicit references to data_ variables and rewrote in terms of new
!!     matrix routines
!!   2007/08/14 17:31 dave
!!   - Added entropy, free energy output
!!   2008/02/04 17:16 dave
!!   - Changes for output to file not stdout
!!   2012/04/16 L.Tong
!!   - Added vdW_xc_energy
!!   2014/01/18 lat
!!   - Added exx_energy for "exact" exchange (Hartree-Fock)
!!   - Added   x_energy for   DFT   exchange
!!   2015/06/08 lat
!!    - Added experimental backtrace
!!   2015/11/10 17:14 dave with TM and NW (Mizuho)
!!    - Adding variables for neutral atom (and generalising hartree_energy to hartree_energy_total_rho)
!!   2015/11/25 08:30 dave with TM and NW (Mizuho)
!!    - Adding neutral atom energy expressions
!!   2016/03/02 17:19 dave
!!    - Moved ion-ion energies into this module, removed hartree_energy_atom_rho
!!   2021/07/23 16:46 dave
!!    Changes to make the DFT energy correct, in particular creating hartree_energy_drho_input which
!!    stores the value of hartree_energy_drho made from the input charge density when we have check_DFT T
!!   2021/07/30 12:15 dave
!!    Remove check_DFT
!!   2022/06/09 12:29 dave
!!    Moved disp_energy here from DFT-D2 module
!!  SOURCE
!!
module energy

  use datatypes
  use global_module, only: io_lun
  use numbers,       only: zero
  use timer_module,  only: start_timer, stop_timer, cq_timer
  use timer_module,  only: start_backtrace, stop_backtrace

  !use DiagModule,             only: nkp
  use ScalapackFormat,        only: matrix_size

  
  implicit none

  ! Area identification
  integer, parameter, private :: area = 9

  real(double) ::  hartree_energy_total_rho
  real(double) :: local_ps_energy
  real(double) ::       xc_energy

  real(double) ::       nl_energy
  real(double) ::     band_energy
  real(double) ::  kinetic_energy
  real(double) ::     cdft_energy
  real(double) ::      exx_energy
  real(double) ::        x_energy
  real(double) ::     disp_energy
  real(double) :: delta_E_hartree
  real(double) :: delta_E_xc
  real(double) :: entropy = zero
  real(double) :: ion_interaction_energy
  
  ! Neutral atom potential
  real(double) :: hartree_energy_drho          ! Self-energy for drho = rho - rho_atom
  real(double) :: hartree_energy_drho_input    ! Self-energy for drho with input charge
  real(double) :: hartree_energy_drho_atom_rho ! Energy for rho_atom in drho potential
  real(double) :: screened_ion_interaction_energy
  
  ! To avoid cyclic dependancy with DiagModule, the local variables here record information needed
  ! from DiagModule
  logical :: flag_check_Diag = .false.
  integer :: SmearingType, MPOrder

!!*** energy

contains

  !!****f* energy/get_energy
  !!
  !!  NAME
  !!   get_energy
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Builds and writes out total energy from 2Tr[KH] and appropriate corrections
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   17:20, 2003/03/10
  !!  MODIFICATION HISTORY
  !!   10:09, 13/02/2006 drb
  !!    Removed all explicit references to data_ variables and rewrote
  !!    in terms of new matrix routines
  !!   2007/08/16 15:31 dave
  !!    Changed output format for energy, added Free Energy for smeared
  !!    electron distribution and Gillan E-0.5*TS kT=0 extrapolation
  !!    (see J. Phys.:Condens. Matter 1, 689 (1989)
  !!   2011/08/29 L.Tong
  !!    Added spin polarisation
  !!   2011/08/02 13:52 dave
  !!    Changes for cDFT
  !!   2012/03/18 L.Tong
  !!    Rewrite for changed in spin implementation
  !!   2012/05/29 L.Tong
  !!   - Changed myid == 0 to inode == ionode for consistency
  !!   - Added output for total spin polarisation for spin polarised
  !!     calculations
  !!   2014/05/29 LAT
  !!   - Added exchange energy calculation
  !!   2015/11/24 08:38 dave
  !!    Adjusted name of hartree_energy to hartree_energy_total_rho for neutral atom implementation
  !!   2015/11/25 08:33 dave with TM and NW (Mizuho)
  !!    Adding neutral atom energy expressions
  !!   2016/03/02 17:20 dave
  !!    Tidying up output for neutral atom
  !!   2018/05/24 19:00 nakata
  !!    Changed matKE, matNL and matNA to be spin_SF dependent
  !!   2021/07/26 11:13 dave
  !!    Removed output of hartree_energy_drho_atom_rho (now included
  !!    in delta_E_hartree)
  !!   2021/07/28 10:55 dave
  !!    Change behaviour to print Harris etc always, and DFT only if
  !!    printDFT = T
  !!  SOURCE
  !!
  subroutine get_energy(total_energy, printDFT, level)

    use datatypes
    use numbers
    use units
    use mult_module,            only: matrix_product_trace, matH,     &
                                      matK, matKE, matNL, matX, matNA
    use GenComms,               only: inode, ionode, cq_warn
    use global_module,          only: iprint_gen, nspin, spin_factor, &
                                      flag_SpinDependentSF,           &
                                      flag_dft_d2,                    &
                                      flag_SCconverged_D2,            &
                                      flag_self_consistent,           &
                                      flag_perform_cdft,              &
                                      flag_vdWDFT,                    &
                                      flag_exx, exx_alpha,            &
                                      flag_neutral_atom, min_layer
    use pseudopotential_common, only: core_correction, &
                                      flag_neutral_atom_projector
    use density_module,         only: electron_number
    use exx_module, only: matK_Xrange
    

    implicit none

    ! Passed variables
    real(double) :: total_energy
    logical, intent(in), optional :: printDFT
    integer, intent(in), optional :: level
    
    ! Local variables
    character(len=80) :: sub_name = "get_energy"
    integer        :: spin, spin_SF
    logical        :: print_Harris, print_DFT
    real(double)   :: total_energy2
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level

    ! electron number information
    real(double), dimension(nspin) :: electrons
    real(double) :: electrons_tot

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='get_energy',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$

    print_DFT    = .false.
    print_Harris = .true.
    if (present(printDFT)) then
       print_DFT = printDFT
    else
       if (iprint_gen >= 2) print_DFT = .true.
    end if

    ! Find energies
    nl_energy      = zero
    band_energy    = zero
    kinetic_energy = zero
    if(flag_neutral_atom_projector) local_ps_energy     = zero
    exx_energy     = zero
    spin_SF = 1
    do spin = 1, nspin
       if (flag_SpinDependentSF) spin_SF = spin
       nl_energy   = nl_energy   + spin_factor * &
                     matrix_product_trace(matK(spin), matNL(spin_SF))
       if(flag_neutral_atom_projector) local_ps_energy = local_ps_energy      &
                        + spin_factor*matrix_product_trace(matK(spin), matNA(spin_SF))
       band_energy = band_energy + spin_factor * &
                     matrix_product_trace(matK(spin), matH(spin))
       ! note that matKE is < phi_i | - grad^2 | phi_j >
       kinetic_energy = kinetic_energy + spin_factor * half * &
                        matrix_product_trace(matK(spin), matKE(spin_SF))
       exx_energy = exx_energy - spin_factor * half * exx_alpha * &
                    matrix_product_trace(matK_Xrange(spin), matX(spin))
    end do

    ! Find exx energy
    !exx_energy = zero
    !if (flag_exx) then
    !   do spin = 1, nspin
    !      exx_energy = exx_energy -                        &
    !           half * spin_factor * exx_alpha *            &
    !           matrix_product_trace(matK(spin), matX(spin))
    !   end do
    !end if

    ! Find total pure DFT energy
    if(flag_neutral_atom) then
       ! NB I have changed this to use delta_E_hartree = -hartree_energy_drho -hartree_energy_drho_atom_rho
       ! which is calculated in H_matrix_module so that we can calculate hartree_energy_drho for the output
       ! DRB 2021/07/23 16:33
       total_energy = band_energy + delta_E_hartree + delta_E_xc + &
            screened_ion_interaction_energy 
    else
       total_energy = band_energy  + delta_E_hartree + delta_E_xc + &
            ion_interaction_energy + core_correction
    end if
    ! Add contribution from constrained (cDFT)
    if (flag_perform_cdft) total_energy = total_energy + cdft_energy

    ! Add contribution from dispersion (DFT-D2)
    if (flag_dft_d2)       total_energy = total_energy + disp_energy

    ! Add contribution from exact-exchange (EXX)
    !if (flag_exx)          total_energy = total_energy + exx_energy

    ! Write out data
    if (inode == ionode) then
       if(print_Harris) then
          !if(iprint_gen>=1) write(io_lun,2) electrons
          if (iprint_gen + min_layer >= 3) then
             if(flag_neutral_atom) then
                write (io_lun, 1) en_conv*band_energy,     en_units(energy_units)
                write (io_lun,33) en_conv*hartree_energy_drho,  en_units(energy_units)
                !write (io_lun,36) en_conv*hartree_energy_drho_atom_rho,  en_units(energy_units)
                write (io_lun, 4) en_conv*(xc_energy+exx_energy), en_units(energy_units)
                write (io_lun,23) en_conv*x_energy,               en_units(energy_units)
                write (io_lun,24) en_conv*(xc_energy-x_energy),   en_units(energy_units)
                write (io_lun,22) en_conv*exx_energy,      en_units(energy_units)
                write (io_lun,35) en_conv*local_ps_energy, en_units(energy_units)
                write (io_lun, 7) en_conv*nl_energy,       en_units(energy_units)
                write (io_lun, 8) en_conv*kinetic_energy,  en_units(energy_units)
                write (io_lun,39) en_conv*screened_ion_interaction_energy,    en_units(energy_units)
                write (io_lun,11) en_conv*delta_E_hartree, en_units(energy_units)
                write (io_lun,12) en_conv*delta_E_xc,      en_units(energy_units)
             else
                write (io_lun, 1) en_conv*band_energy,     en_units(energy_units)
                write (io_lun, 3) en_conv*hartree_energy_total_rho,  en_units(energy_units)
                write (io_lun, 4) en_conv*(xc_energy+exx_energy), en_units(energy_units)
                write (io_lun,23) en_conv*x_energy,               en_units(energy_units)
                write (io_lun,24) en_conv*(xc_energy-x_energy),   en_units(energy_units)
                write (io_lun,22) en_conv*exx_energy,      en_units(energy_units)
                write (io_lun, 5) en_conv*local_ps_energy, en_units(energy_units)
                write (io_lun, 6) en_conv*core_correction, en_units(energy_units)
                write (io_lun, 7) en_conv*nl_energy,       en_units(energy_units)
                write (io_lun, 8) en_conv*kinetic_energy,  en_units(energy_units)
                write (io_lun, 9) en_conv*ion_interaction_energy,    en_units(energy_units)
                write (io_lun,11) en_conv*delta_E_hartree, en_units(energy_units)
                write (io_lun,12) en_conv*delta_E_xc,      en_units(energy_units)
             end if
             if (flag_perform_cdft) write (io_lun,&
                      '(10x,"cDFT Energy, 2Tr[K.W]            : ",f25.15," ",a2)')&
                     en_conv*cdft_energy, en_units(energy_units)

             if (flag_dft_d2) write (io_lun,17) en_conv*disp_energy, en_units(energy_units)
          end if

          if (abs(entropy) >= RD_ERR) then
             if (iprint_gen + min_layer >= 1) &
                  write(io_lun,10) en_conv*total_energy, en_units(energy_units)
             if (flag_check_Diag) then
                select case (SmearingType)
                case (0) ! Fermi smearing
                   if (entropy < zero) &
                        call cq_warn(sub_name,'Calculated entropy is less than zero; something is wrong ', entropy)
                   if (iprint_gen + min_layer >= 2) &
                        write (io_lun,14) en_conv*(total_energy-half*entropy), &
                                          en_units(energy_units)
                case (1) ! Methfessel-Paxton smearing
                   if (iprint_gen + min_layer >= 2)                                     &
                        write (io_lun,16)                                   &
                              en_conv * (total_energy -                     &
                                         (real(MPOrder+1,double) /          &
                                          real(MPOrder+2,double))*entropy), &
                              en_units(energy_units)
                end select
             else
                if (iprint_gen + min_layer >= 2) &
                     write (io_lun,14) en_conv*(total_energy-half*entropy), &
                                       en_units(energy_units)
             end if
             if (iprint_gen + min_layer >= 1) &
                  write (io_lun,15) en_conv*(total_energy-entropy), &
                                    en_units(energy_units)
          else
             if (iprint_gen + min_layer >= 1) &
                  write (io_lun,10) en_conv*total_energy, en_units(energy_units)
             if (iprint_gen + min_layer >= 2) &
                  write (io_lun, '(10x,"(TS=0 as O(N) or entropic &
                                  &contribution is negligible)")')
          end if
       end if ! (print_Harris)
    end if
    ! Check on validity of band energy
    if(print_DFT) then
       if(flag_neutral_atom) then
          total_energy2 = hartree_energy_drho  + &
               xc_energy       + &     
               exx_energy      + &     
               local_ps_energy + &
               nl_energy       + &
               kinetic_energy  + &
               screened_ion_interaction_energy     
       else
          total_energy2 = hartree_energy_total_rho  + &
               xc_energy       + &     
               exx_energy      + &     
               local_ps_energy + &
               nl_energy       + &
               kinetic_energy  + &
               core_correction + &
               ion_interaction_energy     
       end if
       if (flag_perform_cdft) total_energy2 = total_energy2 + cdft_energy
       if (flag_dft_d2)       total_energy2 = total_energy2 + disp_energy

       if (inode == ionode .and. iprint_gen + min_layer>=2) then
          write(io_lun,13) en_conv*total_energy2, en_units(energy_units)
       end if
    end if

    ! print electron number and spin polarisation information
    if (iprint_gen + min_layer >= 0) then
       call electron_number(electrons)
       if (inode == ionode) then
          if (iprint_gen + min_layer >= 2) then
             electrons_tot = electrons(1) + electrons(nspin)
             write (io_lun,18) electrons_tot
             if (nspin == 2) then
                write (io_lun,19) electrons(1)
                write (io_lun,20) electrons(2)
             end if
             if (nspin == 2) &
                  write (io_lun,21) electrons(1) - electrons(2)
          end if
       end if
    end if

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_energy',echo=.true.)
!****lat>$

    return

1   format(10x,'Band Energy, 2Tr[K.H]            : ',f25.15,' ',a2)
2   format(10x,'Electron Count                   : ',f25.15,' ',a2)
3   format(10x,'Hartree Energy (rho)             : ',f25.15,' ',a2)
4   format(10x,'X-C     Energy                   : ',f25.15,' ',a2)
5   format(10x,'Local PsePot  Energy             : ',f25.15,' ',a2)
6   format(10x,'Core Correction Energy           : ',f25.15,' ',a2)
7   format(10x,'NonLocal PsePot Energy           : ',f25.15,' ',a2)
8   format(10x,'Kinetic Energy                   : ',f25.15,' ',a2)
9   format(10x,'Ion-Ion Interaction Energy       : ',f25.15,' ',a2)
10  format(10x,'Harris-Foulkes Energy            : ',f25.15,' ',a2)
14  format(10x,'GroundState Energy (E-(1/2)TS)   : ',f25.15,' ',a2)
16  format(10x,'GroundState Energy (kT --> 0)    : ',f25.15,' ',a2)
15  format(10x,'Free Energy (E-TS)               : ',f25.15,' ',a2)
13  format(10x,'DFT Total Energy                 : ',f25.15,' ',a2)
11  format(10x,'Ha Correction                    : ',f25.15,' ',a2)
12  format(10x,'XC Correction                    : ',f25.15,' ',a2)
17  format(10x,'dispersion (DFT-D2)              : ',f25.15,' ',a2)
18  format(10x,'Total number of electrons        : ',f25.15)
19  format(10x,'Number of electrons spin up      : ',f25.15)
20  format(10x,'Number of electrons spin down    : ',f25.15)
21  format(10x,'Spin polarisation (NeUP - NeDN)  : ',f25.15)
22  format(10x,'EXX Energy, -Tr[K.X]             : ',f25.15,' ',a2)
23  format(10x,'X only  Energy                   : ',f25.15,' ',a2)
24  format(10x,'C only  Energy                   : ',f25.15,' ',a2)
33  format(10x,'Hartree Energy (delta rho)       : ',f25.15,' ',a2)
35  format(10x,'Neutral Atom  Energy             : ',f25.15,' ',a2)
!36  format(10x,'Hartree Energy (drho-atom rho)   : ',f25.15,' ',a2)
39  format(10x,'Screened Ion-Ion Energy          : ',f25.15,' ',a2)

  end subroutine get_energy
  !!*** get_energy


  !!****f* energy/final_energy
  !!
  !!  NAME
  !!   final_energy
  !!  USAGE
  !!
  !!  PURPOSE
  !!  Writes out energy contributions at the very last step
  !!   of the calculation 
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler/T.Miyazaki/L.Truflandier
  !!  CREATION DATE
  !!   2014/09/24
  !!  MODIFICATION HISTORY
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2015/11/30 17:13 dave
  !!    - Added neutral atom output
  !!   2018/11/13 17:30 nakata
  !!    Changed matS, matKE, matNL and matNA to be spin_SF dependent
  !!    Activated "electrons_tot2" calculation
  !!   2019/12/03 08:09 dave
  !!    Removed broken code for electrons via Tr[KS]
  !!   2021/10/28 17:13 lionel
  !!    Added 'DFT total energy' printing in ASE output
  !!    nkp must to be passed to avoid circular dependency
  !!   2023/01/10 18:44 lionel
  !!    Secure ASE printing when using ordern
  !!  SOURCE
  !!
  subroutine final_energy(nkp,level)

    use datatypes
    use numbers
    use units
    use GenComms,               only: inode, ionode, cq_warn, cq_abort
    use mult_module,            only: matrix_product_trace, matH,     &
                                      matrix_product_trace_length,    &
                                      matrix_trace,                   &
                                      matK, matKE, matNL, matX, matS, matNA

    use global_module,          only: iprint_gen, nspin, spin_factor, &
                                      flag_SpinDependentSF,           &
                                      flag_dft_d2,                    &
                                      flag_SCconverged_D2,            &
                                      flag_self_consistent,           &
                                      flag_perform_cdft,              &
                                      flag_vdWDFT,                    &
                                      flag_exx, exx_alpha,            &
                                      flag_neutral_atom, min_layer,   &
                                      flag_fix_spin_population,       &
                                      io_ase, write_ase, ase_file,    &
                                      flag_diagonalisation

    use density_module,         only: electron_number
    use pseudopotential_common, only: core_correction, &
                                      flag_neutral_atom_projector
    use species_module,         only: n_species
    use input_module,           only: io_close
    use exx_module, only: matK_Xrange
    
    implicit none

    ! Passed variables
    integer, optional   :: level
    integer, intent(in) :: nkp
    
    ! Local variables
    character(len=80) :: sub_name = "final_energy"
    integer        :: spin, spin_SF
    real(double)   :: total_energy1
    real(double)   :: total_energy2
    real(double)   :: one_electron_energy
    real(double)   :: potential_energy
    real(double)   :: virial
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level

    ! electron number information
    real(double), dimension(nspin) :: electrons
    real(double)                   :: electrons_tot, electrons_tot2
    integer :: i, stat, counter
    character(len=80) :: tmp

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='final_energy',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$

    ! Initialise energies
    nl_energy           = zero
    band_energy         = zero
    kinetic_energy      = zero
    if(flag_neutral_atom_projector) local_ps_energy     = zero
    exx_energy          = zero
    one_electron_energy = zero
    potential_energy    = zero
    total_energy1       = zero
    total_energy2       = zero

    spin_SF = 1

    ! Nonlocal pseudop, kinetic and band energies
    do spin = 1, nspin
       if (flag_SpinDependentSF) spin_SF = spin
       ! 2*Tr[K NL]
       nl_energy      = nl_energy      &
                        + spin_factor*matrix_product_trace(matK(spin), matNL(spin_SF))
       if(flag_neutral_atom_projector) local_ps_energy = local_ps_energy      &
                        + spin_factor*matrix_product_trace(matK(spin), matNA(spin_SF))
       ! 2*Tr[K KE] with KE = - < grad**2 >
       kinetic_energy = kinetic_energy &
                        + spin_factor*half*matrix_product_trace(matK(spin), matKE(spin_SF))
       ! 2*Tr[K H]
       band_energy    = band_energy    &
                        + spin_factor*matrix_product_trace(matK(spin), matH(spin))
       ! -alpha*Tr[K X]
       exx_energy     = exx_energy     &
                        - spin_factor*half*exx_alpha*matrix_product_trace(matK_Xrange(spin), matX(spin))
    end do

    ! Find total pure DFT energy
    if(flag_neutral_atom) then
       ! NB I have changed this to use delta_E_hartree = -hartree_energy_drho -hartree_energy_drho_atom_rho
       ! which is calculated in H_matrix_module so that we can calculate hartree_energy_drho for the output
       ! DRB 2021/07/23 16:33
       total_energy1 = band_energy     + &
            delta_E_hartree + &
            delta_E_xc + &
            screened_ion_interaction_energy
    else
       total_energy1 = band_energy     + &
            delta_E_hartree + &
            delta_E_xc      + &
            ion_interaction_energy    + &
            core_correction
    end if

    ! Add contribution from constrained (cDFT)
    if (flag_perform_cdft) total_energy1 = total_energy1 + cdft_energy

    ! Add contribution from dispersion (DFT-D2)
    if (flag_dft_d2)       total_energy1 = total_energy1 + disp_energy

    !Write out data
    !...
    !
    ! print electron number and spin polarisation information
    
    call electron_number(electrons)
   
    !if (inode == ionode) write(io_lun,*) 'electrons_tot2 start'
    !electrons_tot2 = matrix_product_trace_length(matK(1),matS(1))
    !if (inode == ionode) write(io_lun,*) 'electrons_tot2 stop'

    if (inode == ionode) then
       electrons_tot = electrons(1) + electrons(nspin)
       !
       if(iprint_gen + min_layer<2 .and. iprint_gen + min_layer>-1) then
          if (nspin == 1) then
             write (io_lun,fmt='(6x,"| Number of electrons      = ",f16.6)') electrons_tot
          else if(nspin == 2) then
             write (io_lun,fmt='(6x,"| Number of electrons (u/d)= ",2f16.6)') electrons(1),electrons(nspin)
          end if
       else if (iprint_gen + min_layer ==2) then
          if (nspin == 1) then
             write (io_lun,fmt='(6x,"| Number of electrons      = ",f16.6)') electrons_tot
          else if(nspin == 2) then
             write (io_lun,fmt='(6x,"| Number of electrons (u/d)= ",2f16.6)') electrons(1),electrons(nspin)
          end if
          write (io_lun, 6) en_conv *    band_energy,  en_units(energy_units)
          if(flag_neutral_atom) then
             write (io_lun,68) en_conv * screened_ion_interaction_energy,  en_units(energy_units)
          else
             write (io_lun, 8) en_conv *   ion_interaction_energy,  en_units(energy_units)
          end if
          if(flag_neutral_atom) then
             write (io_lun,11) en_conv*(- hartree_energy_drho - hartree_energy_drho_atom_rho), &
                  en_units(energy_units)
          else
             write (io_lun,11) en_conv* delta_E_hartree, en_units(energy_units)
          end if
          write (io_lun,12) en_conv* delta_E_xc,      en_units(energy_units)
       else if (iprint_gen + min_layer > 2) then
          write (io_lun, *) 
          !write (io_lun, *) 
          !write (io_lun, 1) 
          !
          if (nspin == 1) then
             write (io_lun,fmt='(6x,"| Number of electrons      = ",f25.15)') electrons_tot
          else if (nspin == 2) then
             write (io_lun,19) electrons(1)
             write (io_lun,20) electrons(2)
             write (io_lun,21) electrons(1) - electrons(2)
             !
          end if
          !
          write (io_lun, 6) en_conv *    band_energy,  en_units(energy_units)
          if(flag_neutral_atom) then
             write (io_lun,67) en_conv * hartree_energy_drho,  en_units(energy_units)
             write (io_lun,68) en_conv * screened_ion_interaction_energy,  en_units(energy_units)
          else
             write (io_lun, 7) en_conv * hartree_energy_total_rho,  en_units(energy_units)
             write (io_lun, 8) en_conv *   ion_interaction_energy,  en_units(energy_units)
          end if
          write (io_lun, 9) en_conv * kinetic_energy,  en_units(energy_units)
          write (io_lun, 2)
          write (io_lun,30) en_conv* (xc_energy + exx_energy),  en_units(energy_units)
          write (io_lun,31) en_conv*   x_energy,                en_units(energy_units)
          write (io_lun,32) en_conv* (xc_energy - x_energy),    en_units(energy_units)
          write (io_lun,33) en_conv* exx_energy,  en_units(energy_units)
          write (io_lun, 2)
          if(flag_neutral_atom) then
             write (io_lun,60) en_conv* (local_ps_energy + &
                                         nl_energy    ),  en_units(energy_units)
             write (io_lun,62) en_conv* local_ps_energy,  en_units(energy_units)
          else
             write (io_lun,40) en_conv* (core_correction + &
                                         local_ps_energy + &
                                         nl_energy    ),  en_units(energy_units)
             write (io_lun,41) en_conv* core_correction,  en_units(energy_units) 
             write (io_lun,42) en_conv* local_ps_energy,  en_units(energy_units)
          end if
          write (io_lun,43) en_conv*       nl_energy,  en_units(energy_units)
          write (io_lun, 2)
          write (io_lun,11) en_conv* delta_E_hartree, en_units(energy_units)
          write (io_lun,12) en_conv* delta_E_xc,      en_units(energy_units)
          if (flag_dft_d2) &
               write (io_lun,17) en_conv*disp_energy,en_units(energy_units)
          if (flag_perform_cdft) &          
               write (io_lun,18) en_conv*cdft_energy,en_units(energy_units)
          write (io_lun, 2)
       end if
    end if
    
    if ( inode == ionode ) then
       !
       if (abs(entropy) >= RD_ERR) then
          
          !if (iprint_gen >= 0) &
          !     write(io_lun,10) en_conv*total_energy1, en_units(energy_units)
       
          if (flag_check_Diag) then
             !
             select case (SmearingType)
             case (0) ! Fermi smearing
                if (entropy < zero) &
                     call cq_warn(sub_name, 'Calculated entropy is less than zero; something is wrong ', entropy)
                !
                if (iprint_gen + min_layer >= 1) &
                     write (io_lun,14) en_conv*(total_energy1-half*entropy), &
                     en_units(energy_units)
                !
                !
             case (1) ! Methfessel-Paxton smearing
                if (iprint_gen + min_layer >= 1)                             &
                     write (io_lun,16) en_conv * (total_energy1 - &
                     (real(MPOrder+1,double) /                   &
                     real(MPOrder+2,double))*entropy),          &
                     en_units(energy_units)
                !
                !
             end select
             !
          else
             if (iprint_gen + min_layer >= 1) &
                  write (io_lun,14) en_conv*(total_energy1-half*entropy), &
                  en_units(energy_units)
          end if
          !
          !
          if (iprint_gen + min_layer >= 1) &
               write (io_lun,15) en_conv*(total_energy1-entropy), &
               en_units(energy_units)
       else
          !if (iprint_gen >= 0) &
          !     write (io_lun,10) en_conv*total_energy1, en_units(energy_units)
          !if (iprint_gen >= 0) &
          !     write (io_lun, '(10x,"(TS=0 as O(N) or entropic &
          !     &contribution is negligible)")')
       end if
    end if
    
    ! Check on validity of band energy
    if(flag_neutral_atom) then
       total_energy2 = hartree_energy_drho  + &
            xc_energy       + &     
            exx_energy      + &
            local_ps_energy + &
            nl_energy       + &
            kinetic_energy  + &
            screened_ion_interaction_energy
    else
       total_energy2 = hartree_energy_total_rho  + &
            xc_energy       + &     
            exx_energy      + &
            local_ps_energy + &
            nl_energy       + &
            kinetic_energy  + &
            core_correction + &
            ion_interaction_energy     
    end if
    if (flag_perform_cdft) total_energy2 = total_energy2 + cdft_energy
    if (flag_dft_d2)       total_energy2 = total_energy2 + disp_energy

    ! One-electron energy
    one_electron_energy = local_ps_energy + &
                          nl_energy       + &
                          kinetic_energy  

    ! Potential energy
    potential_energy = local_ps_energy + &
                       nl_energy       + &
                       hartree_energy_total_rho

    if (inode == ionode) then
       if (iprint_gen + min_layer >= 0) then
          write(io_lun,10) en_conv*total_energy1, en_units(energy_units)
          if(iprint_gen + min_layer>0) then
             
             write(io_lun,13)  en_conv*total_energy2, en_units(energy_units)             
             if(iprint_gen + min_layer>1) &
                  write(io_lun,22) en_conv*(total_energy1 - total_energy2), &
                  en_units(energy_units) 
          end if
       end if
       !
       ! BEGIN %%%% ASE printing %%%%
       !
       if ( write_ase ) then
          !
          open(io_ase,file=ase_file, status='old', action='readwrite', iostat=stat, position='rewind')
          !
          if (stat .ne. 0) call cq_abort('ASE/energy error opening file !')
          !
          if ( flag_diagonalisation ) then
             if ( nspin == 2 ) then
                counter = nkp*3 + (nspin+1)*nkp + nspin*nkp*(matrix_size/3) + 1 + 2 
                if ( mod(matrix_size,3) > 0 ) counter = counter + nspin*nkp
                
             else
                counter = nkp*3 + nkp*(matrix_size/3) + 1 + 1
                if ( mod(matrix_size,3) > 0 ) counter = counter + nkp
             
             end if
             counter = counter + 7 + n_species + nkp + 1 + 1
             
          else
             counter =  7 + n_species
          end if
          
          do i = 1, counter
             read (io_ase,*)
          end do
          !
          write(io_ase,133) en_conv*total_energy2, en_units(energy_units)       
          !
          close(io_ase)
          !
          ! END %%%% ASE printing %%%%
          !
       end if
    end if
    
    !call electron_number(electrons)
    !
    !electrons_tot2 = zero
    !do spin = 1, nspin
    !   if (flag_SpinDependentSF) spin_SF = spin
    !   electrons_tot2 = electrons_tot2 + &
    !        spin_factor * matrix_product_trace_length(matK(spin),matS(spin_SF))
    !end do
    !
    !if (inode == ionode) then
    !   electrons_tot = electrons(1) + electrons(nspin)
    !   if(iprint_gen + min_layer>-1) then
    !      write (io_lun,241) electrons_tot
    !      if (nspin == 2 .and. (.not. flag_fix_spin_population)) &
    !           write (io_lun,21) electrons(1) - electrons(2)
    !   elseif (iprint_gen + min_layer > 2) then
    !      write (io_lun,23) 
    !      write (io_lun,24) electrons_tot
    !      write (io_lun,25) electrons_tot2
    !      !write (io_lun,26) one_electron_energy
    !      !write (io_lun,27) potential_energy
    !      !write (io_lun,28) kinetic_energy
    !      write (io_lun,29) (total_energy2 - kinetic_energy)/kinetic_energy 
    !      if (exx_alpha > zero .and. exx_alpha < one) then             
    !         write (io_lun,50) x_energy/(one - exx_alpha)
    !         write (io_lun,51) exx_energy/exx_alpha
    !      end if
    !      !write (io_lun, 1)
    !      !write (io_lun, *) 
    !      write (io_lun, *) 
    !   end if
    !end if


!****lat<$
    call stop_backtrace(t=backtrace_timer,who='final_energy',echo=.true.)
!****lat>$

    return


1   format(6x, '****===============================', &
               '===================================', &
               '===========****')

2   format(6x, ' ')


6   format(6x,'|* band energy as 2Tr[K.H] = ',f25.15,' ',a2)
7   format(6x,'|  hartree energy (rho)    = ',f25.15,' ',a2)
8   format(6x,'|  ion-ion energy          = ',f25.15,' ',a2)
67  format(6x,'|  hartree energy (drho)   = ',f25.15,' ',a2)
68  format(6x,'|  screened ion-ion energy = ',f25.15,' ',a2)
9   format(6x,'|  kinetic energy          = ',f25.15,' ',a2)

30  format(6x,'|* xc total energy         = ',f25.15,' ',a2)
31  format(6x,'|    DFT exchange          = ',f25.15,' ',a2)
32  format(6x,'|    DFT correlation       = ',f25.15,' ',a2)
33  format(6x,'|    EXX contribution      = ',f25.15,' ',a2)


40  format(6x,'|* pseudopotential energy  = ',f25.15,' ',a2)
41  format(6x,'|    core correction       = ',f25.15,' ',a2)
42  format(6x,'|    local contribution    = ',f25.15,' ',a2)
43  format(6x,'|    nonlocal contribution = ',f25.15,' ',a2)
60  format(6x,'|* pseudo/NA energy        = ',f25.15,' ',a2)
62  format(6x,'|    NA contribution       = ',f25.15,' ',a2)

11  format(6x,'|  Ha correction           = ',f25.15,' ',a2)
12  format(6x,'|  XC correction           = ',f25.15,' ',a2)

10  format(6x,'|* Harris-Foulkes energy   = ',f25.15,' ',a2)
13  format(6x,'|* DFT total energy        = ',f25.15,' ',a2)
133 format(/4x,'DFT total energy        = ',f25.15,' ',a2)
    
22  format(6x,'|  estimated accuracy      = ',f25.15,' ',a2)

14  format(6x,'| GS Energy as E-(1/2)TS   = ',f25.15,' ',a2)
15  format(6x,'| Free Energy as E-TS      = ',f25.15,' ',a2)
16  format(6x,'| GS Energy with kT -> 0   = ',f25.15,' ',a2)
17  format(6x,'| Dispersion (DFT-D2)      = ',f25.15,' ',a2)
18  format(6x,'| cDFT Energy as 2Tr[K.W]  = ',f25.15,' ',a2)
19  format(6x,'| Number of e- spin up     = ',f25.15)
20  format(6x,'| Number of e- spin down   = ',f25.15)
21  format(6x,'| Spin pol. as (up - down) = ',f25.15)

23  format(6x,'|* check for accuracy      = ')
24  format(6x,'|  number of e- num. int.  = ',f25.15)
241 format(6x,'|  number of electrons     = ',f25.15)
25  format(6x,'|  number of e- as 2Tr[KS] = ',f25.15)
26  format(6x,'|  one-electron energy     = ',f25.15,' ',a2)
27  format(6x,'|  potential energy V      = ',f25.15,' ',a2)
28  format(6x,'|  kinetic energy T        = ',f25.15,' ',a2)
29  format(6x,'|* virial V/T              = ',f25.15,' ',a2)
50  format(6x,'|  rescaled DFT exchange   = ',f25.15,' ',a2)
51  format(6x,'|  rescaled exact exchange = ',f25.15,' ',a2)


  end subroutine final_energy
  !!*** get_energy


end module energy
