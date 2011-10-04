! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!    Added band_energy
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2007/08/14 17:31 dave
!!    Added entropy, free energy output
!!   2008/02/04 17:16 dave
!!    Changes for output to file not stdout
!!  SOURCE
!!
module energy

  use datatypes
  use global_module, ONLY: io_lun
  use numbers, ONLY: zero

  implicit none

  real(double) :: hartree_energy
  real(double) :: local_ps_energy
  real(double) :: xc_energy
  real(double) :: nl_energy
  real(double) :: band_energy
  real(double) :: kinetic_energy
  real(double) :: delta_E_hartree
  real(double) :: delta_E_xc
  real(double) :: entropy=zero

  logical :: flag_check_DFT=.false.

  ! To avoid cyclic dependancy with DiagModule, the local variables here record information needed 
  ! from DiagModule
  logical :: flag_check_Diag=.false.
  integer :: SmearingType, MPOrder

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
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
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2007/08/16 15:31 dave
!!    Changed output format for energy, added Free Energy for smeared electron distribution
!!    and Gillan E-0.5*TS kT=0 extrapolation (see J. Phys.:Condens. Matter 1, 689 (1989)
!!  SOURCE
!!
  subroutine get_energy(total_energy,printDFT)

    use datatypes
    use numbers
    use units
    use mult_module, ONLY: matrix_product_trace, matH, matK, matKE, matNL
    use GenComms, ONLY: myid
    use global_module, ONLY: iprint_gen, flag_dft_d2, flag_SCconverged_D2, flag_self_consistent
    use ewald_module, ONLY: ewald_energy
    use pseudopotential_common, ONLY: core_correction  
    use DFT_D2, ONLY: disp_energy

    implicit none

    ! Passed variables
    real(double) :: total_energy

    ! Local variables
    real(double) :: total_energy2
    ! check DFT energy mode  : by TM Nov2007
    logical, intent(in), OPTIONAL :: printDFT
    logical :: print_Harris, print_DFT

 ! For Now, 
 ! If printDFT does not exist (as in the previous version), 
 !  DFT energy will be printed out if iprint_gen >= 2
 !
 ! If it exists,
 !  printDFT = .true.  -> prints out only DFT energy
 !  printDFT = .false. -> prints out Energies except DFT energy
 !  
    print_DFT=.false.
    print_Harris = .true.
    if(PRESENT(printDFT)) then
       if(printDFT) print_Harris = .false.
       print_DFT=printDFT
    else
       if(iprint_gen>=2) print_DFT=.true.
    endif
    
    ! Find energies
    nl_energy = two*matrix_product_trace(matK,matNL)
    kinetic_energy = matrix_product_trace(matK,matKE)
    band_energy = two*matrix_product_trace(matK,matH)
    total_energy = band_energy + delta_E_hartree + delta_E_xc + ewald_energy + core_correction
    ! for DFT-D2
    if (flag_dft_d2) total_energy = total_energy + disp_energy
    ! Write out data
    if(myid==0) then
       if(print_Harris) then
          !if(iprint_gen>=1) write(io_lun,2) electrons
          if(iprint_gen>=1) write(io_lun,1) en_conv*band_energy, en_units(energy_units)
          if(iprint_gen>=1) write(io_lun,3) en_conv*hartree_energy, en_units(energy_units)
          if(iprint_gen>=1) write(io_lun,4) en_conv*xc_energy, en_units(energy_units)
          if(iprint_gen>=1) write(io_lun,5) en_conv*local_ps_energy, en_units(energy_units)
          if(iprint_gen>=1) write(io_lun,6) en_conv*core_correction, en_units(energy_units)
          if(iprint_gen>=1) write(io_lun,7) en_conv*nl_energy, en_units(energy_units)
          if(iprint_gen>=1) write(io_lun,8) en_conv*kinetic_energy, en_units(energy_units)
          if(iprint_gen>=1) write(io_lun,9) en_conv*ewald_energy, en_units(energy_units)
          if(iprint_gen>=1) write(io_lun,11) en_conv*delta_E_hartree, en_units(energy_units)
          if(iprint_gen>=1) write(io_lun,12) en_conv*delta_E_xc, en_units(energy_units)
          if (flag_dft_d2) then
            if (iprint_gen>=1) write(io_lun,17) en_conv*disp_energy, en_units(energy_units)
          endif
          if(abs(entropy) >= very_small) then
             if(iprint_gen>=0) write(io_lun,10) en_conv*total_energy, en_units(energy_units)
             if(flag_check_Diag) then
                select case (SmearingType)
                case (0) ! Fermi smearing
                   if(entropy < zero) write(io_lun,*) ' WARNING !!!!    entropy < 0??? ', entropy
                   if(iprint_gen>=0) write(io_lun,14) en_conv*(total_energy- half*entropy), en_units(energy_units)
                case (1) ! Methfessel-Paxton smearing
                   if(iprint_gen>=0) write(io_lun,16) &
                        en_conv*(total_energy - (real(MPOrder+1,double)/real(MPOrder+2,double))*entropy), &
                        en_units(energy_units)
                end select
             else 
                if(iprint_gen>=0) write(io_lun,14) en_conv*(total_energy- half*entropy), en_units(energy_units)
             end if
             if(iprint_gen>=1) write(io_lun,15) en_conv*(total_energy - entropy), en_units(energy_units)
          else
             if(iprint_gen>=0) write(io_lun,10) en_conv*total_energy, en_units(energy_units)
             if(iprint_gen>=0) write(io_lun,*) '  (TS=0 as O(N) or entropic contribution is negligible) '
          endif
       endif ! (print_Harris)
    end if
    ! Check on validity of band energy
    if(print_DFT) then
       total_energy2 = hartree_energy + xc_energy + local_ps_energy + &
            nl_energy + kinetic_energy + core_correction + ewald_energy
       if (flag_dft_d2) then
         total_energy2 = hartree_energy + xc_energy + local_ps_energy + &
              nl_energy + kinetic_energy + core_correction + ewald_energy + disp_energy
       endif
       if(myid==0) write(io_lun,13) en_conv*total_energy2, en_units(energy_units)
    end if
    return
1   format(10x,'Band Energy, 2Tr[K.H]            : ',f25.15,' ',a2)
2   format(10x,'Electron Count                   : ',f25.15,' ',a2)
3   format(10x,'Hartree Energy                   : ',f25.15,' ',a2)
4   format(10x,'X-C     Energy                   : ',f25.15,' ',a2)
5   format(10x,'Local PsePot  Energy             : ',f25.15,' ',a2)
6   format(10x,'Core Correction Energy           : ',f25.15,' ',a2)
7   format(10x,'NonLocal PsePot Energy           : ',f25.15,' ',a2)
8   format(10x,'Kinetic Energy                   : ',f25.15,' ',a2)
9   format(10x,'Ewald Energy                     : ',f25.15,' ',a2)
10  format(10x,'Harris-Foulkes Energy            : ',f25.15,' ',a2)
14  format(10x,'GroundState Energy (E-(1/2)TS)   : ',f25.15,' ',a2)
16  format(10x,'GroundState Energy (kT --> 0)    : ',f25.15,' ',a2)
15  format(10x,'Free Energy (E-TS)               : ',f25.15,' ',a2)
13  format(10x,'DFT Total Energy                 : ',f25.15,' ',a2)
11  format(10x,'Ha Correction                    : ',f25.15,' ',a2)
12  format(10x,'XC Correction                    : ',f25.15,' ',a2)
17  format(10x,'dispersion (DFT-D2)              : ',f25.15,' ',a2)
  end subroutine get_energy
!!*** get_energy

end module energy
