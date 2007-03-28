! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: energy.module.f90,v 1.4.2.2 2006/03/07 07:36:42 drb Exp $
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
!!  SOURCE
!!
module energy

  use datatypes

  implicit none

  real(double) :: hartree_energy
  real(double) :: local_ps_energy
  real(double) :: xc_energy
  real(double) :: nl_energy
  real(double) :: band_energy
  real(double) :: kinetic_energy
  real(double) :: delta_E_hartree
  real(double) :: delta_E_xc

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id: energy.module.f90,v 1.4.2.2 2006/03/07 07:36:42 drb Exp $"
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
!!  SOURCE
!!
  subroutine get_energy(total_energy)

    use datatypes
    use numbers
    use units
    use mult_module, ONLY: matrix_product_trace, matH, matK, matKE, matNL
    use GenComms, ONLY: myid
    use global_module, ONLY: iprint_gen
    use ewald_module, ONLY: ewald_energy
    use pseudopotential_common, ONLY: core_correction    

    implicit none

    ! Passed variables
    real(double) :: total_energy

    ! Local variables
    real(double) :: total_energy2

    ! Find energies
    nl_energy = two*matrix_product_trace(matK,matNL)
    kinetic_energy = matrix_product_trace(matK,matKE)
    band_energy = two*matrix_product_trace(matK,matH)
    total_energy = band_energy + delta_E_hartree + delta_E_xc + ewald_energy + core_correction
    ! Write out data
    if(myid==0) then
       !if(iprint_gen>=1) write(*,2) electrons
       if(iprint_gen>=1) write(*,1) en_conv*band_energy, en_units(energy_units)
       if(iprint_gen>=1) write(*,3) en_conv*hartree_energy, en_units(energy_units)
       if(iprint_gen>=1) write(*,4) en_conv*xc_energy, en_units(energy_units)
       if(iprint_gen>=1) write(*,5) en_conv*local_ps_energy, en_units(energy_units)
       if(iprint_gen>=1) write(*,6) en_conv*core_correction, en_units(energy_units)
       if(iprint_gen>=1) write(*,7) en_conv*nl_energy, en_units(energy_units)
       if(iprint_gen>=1) write(*,8) en_conv*kinetic_energy, en_units(energy_units)
       if(iprint_gen>=1) write(*,9) en_conv*ewald_energy, en_units(energy_units)
       if(iprint_gen>=1) write(*,11) en_conv*delta_E_hartree, en_units(energy_units)
       if(iprint_gen>=1) write(*,12) en_conv*delta_E_xc, en_units(energy_units)
       if(iprint_gen>=0) write(*,10) en_conv*total_energy, en_units(energy_units)
    end if
    ! Check on validity of band energy
    if(iprint_gen>=2) then
       total_energy2 = hartree_energy + xc_energy + local_ps_energy + &
            nl_energy + kinetic_energy + core_correction + ewald_energy
       if(myid==0) write(*,13) en_conv*total_energy2, en_units(energy_units)
    end if
    return
1   format(20x,'Band Energy, 2Tr[K.H]  : ',f25.15,' ',a2)
2   format(20x,'Electron Count         : ',f25.15,' ',a2)
3   format(20x,'Hartree Energy         : ',f25.15,' ',a2)
4   format(20x,'X-C     Energy         : ',f25.15,' ',a2)
5   format(20x,'Local PsePot  Energy   : ',f25.15,' ',a2)
6   format(20x,'Core Correction Energy : ',f25.15,' ',a2)
7   format(20x,'NonLocal PsePot Energy : ',f25.15,' ',a2)
8   format(20x,'Kinetic Energy         : ',f25.15,' ',a2)
9   format(20x,'Ewald Energy           : ',f25.15,' ',a2)
10  format(20x,'Harris-Foulkes Energy  : ',f25.15,' ',a2)
13  format(20x,'DFT Total Energy       : ',f25.15,' ',a2)
11  format(20x,'Ha Correction          : ',f25.15,' ',a2)
12  format(20x,'XC Correction          : ',f25.15,' ',a2)
  end subroutine get_energy
!!*** get_energy

end module energy
