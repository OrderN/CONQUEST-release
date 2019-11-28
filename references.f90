!!****h* Conquest/references *
!!  NAME
!!   references
!!  PURPOSE
!!   Literal data for bibliography generation, plus wrapper for citing all 
!!   relevant references.
!!  AUTHOR
!!   Zamaan Raza
!!  CREATION DATE
!!   2019/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
module references

  use biblio

  implicit none

  contains

    !!****f* referencds/get_bib_db *
    !!
    !!  NAME 
    !!   compile_biblio
    !!  PURPOSE
    !!   Generate a bibliography database
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/07/04
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    subroutine get_bib_db(bib)

      ! Passed variables
      type(type_bibliography), intent(inout) :: bib

      call bib%add_ref( &
"Arita2014", &
"M. Arita and D. R. Bowler and T. Miyazaki", &
"Stable and efficient linear scaling first-principles molecular dyanmamics for 10000+ atoms", &
"J. Chem. Theor. Comput.", &
10, &
5419, &
2014, &
"10.1021/ct500847y", &
"CONQUEST XL-BOMD implementation")

      call bib%add_ref( &
"Bowler2001", &
"D. R. Bowler and T. Miyazaki and M. J. Gillan", &
"Parallel sparse matrix multiplication for linear scaling electronic structure calculations", &
"Comp. Phys. Commun.", &
137, &
255, &
2001, &
"10.1016/S0010-4655(01)00164-3", &
"Order(N) matrix multiplication")

      call bib%add_ref( &
"Nakata2014", &
"A. Nakata and D. R. Bowler and T. Miyazaki", &
"Efficient calculations with multisite local orbitals in a large scale DFT code CONQUEST", &
"J. Chem. Theor. Comput.", &
10, &
4813, &
2014, &
"10.1021/ct5004934", &
"CONQUEST multisite support functions")

      call bib%add_ref( &
"Gillan2007", &
"M. J. Gillan and D. R. Bowler and A. S. Torralba and T. Miyazaki", &
"Order-N first-principles calculations with the CONQUEST code", &
"Comp. Phys. Commun.", &
177, &
14, &
2007, &
"10.1016/j.cpc.2007.02.075", &
"CONQUEST order(N) method")

      call bib%add_ref( &
"Bowler2012", &
"D. R. Bowler and T. Miyazaki", &
"O(N) methods in electronic structure calculations", &
"Rep. Prog. Phys.", &
75, &
036503, &
2012, &
"10.1088/0034-4885/75/3/036503", &
"Order(N) review paper")

      call bib%add_ref( &
"Hirakawa2017", &
"T. Hirakawa and T. Suzuki and D. R. Bowler and T. Miyazaki", &
"Canonical-ensemble extended Lagrangian Born-Oppenheimer molecular dynamics for the linear scaling density functional theory", &
"J. Phys. Condens. Matter", &
29, &
405901, &
2017, &
"10.1088/1361-648X/aa810d", &
"CONQUEST canonical ensemble MD with XL-BOMD")

      call bib%add_ref( &
"Bowler2010", &
"D. R. Bowler and T. Miyazaki", &
"Calculations for millions of atoms with density functional theory: linear scaling shows its potential", &
"J. Phys.: Condens. Matter", &
22, &
074207, &
2010, &
"10.1088/0953-8984/22/7/074207", &
"CONQUEST order(N) paper")

      call bib%add_ref( &
"Niklasson2006", &
"A. M. N. Niklasson and C. J. Tymczak and M. Challacombe", &
"Time-reversible Born-Oppenheimer molecular dynamics", &
"Phys. Rev. Lett.", &
97, &
123001, &
2006, &
"10.1103/PhysRevLett.97.123001", &
"Original extended-Lagrangian Born-Oppenheimer MD paper")

      call bib%add_ref( &
"Niklasson2009", &
"A. M. N. Niklasson and P. Steneteg and A. Odell and N. Bock and M. Challacombe &
 &and C. J. Tymczak and E. Holmstrom and G. Zhang and V. Weber", &
"Extended Lagrangian Born-Oppenheimer molecular dynamics with dissipation", &
"J. Chem. Phys.", &
130, &
214109, &
2009, &
"10.1063/1.3148075", &
"XL-BOMD with dissipation")

      call bib%add_ref( &
"Torralba2008", &
"A. S. Torralba and M. Todorovic and V. Brazdova and R. Choudhury and T. Miyazaki and M. J. Gillan and D. R. Bowler", &
"Pseudo-atomic orbitals as basis sets for O(N) DFT code CONQUEST", &
"J. Phys.: Condens. Matter", &
20, &
294206, &
2008, &
"10.1088/0953-8984/20/29/294206", &
"PAO basis sets for CONQUEST")

      call bib%add_ref( &
"Bowler2019", &
"D. R. Bowler and J. S. Baker and J. T. L. Poulton and S. Y. Mujahed and J. Lin and S. Yadav and Z. Raza and T. Miyazaki", &
"Highly accurate local basis sets for large-scale DFT calculations in CONQUEST", &
"Japan. J. Appl. Phys. ", &
58, &
100503, &
2019, &
"10.7567/1347-4065/ab45af", &
"Details of new CONQUEST PAOs from basis generation tool, using Hamann-type pseudopotentials")

      call bib%add_ref( &
"Martyna1994", &
"G. J. Martyna and D. J. Tobias and M. L. Klein", &
"Constant pressure molecular dynamics algorithsm", &
"J. Chem. Phys.", &
101, &
4177, &
1994, &
"10.1063/1.467468", &
"Equations of motion for NPT molecular dynamics")

      call bib%add_ref( &
"Martyna1996", &
"G. J. Martyna and M. E. Tuckerman, D. J. Tobias and M. L. Klein", &
"Explicit reversible integrators for extended systems dynamics", &
"Mol. Phys.", &
87, &
1117, &
1996, &
"10.1080/002689799163235", &
"Time-reversible integrators for extended system MD")

      call bib%add_ref( &
"Shinoda2004", &
"W. Shinoda and M. Shiga and M. Mikami", &
"Rapid estimation of elastic constants by molecular dynamics simulation under constant stress", &
"Phys. Rev. B", &
69, &
134103, &
2004, &
"10.1103/PhysRevB.69.134103", &
"Liouvillian splitting for NPT MD used in CONQUEST")

      call bib%add_ref( &
"Hamann2013", &
"D. R. Hamann", &
"Optimized norm-conserving Vanderbilt pseudopotentials", &
"Phys. Rev. B", &
88, &
085117, &
2013, &
"10.1103/PhysRevB.88.085117", &
"Hamann-type pseudopotentials generated via ONCVPSP")

      call bib%add_ref( &
"Miyazakij2004", &
"T. Miyazaki and D. R. Bowler and R. Choudhury and M. J. Gillan", &
"Atomic force algorithms in density functional theory electronic-structure techniques based on local orbitals", &
"J. Chem. Phys.", &
121, &
6186, &
2004, &
"10.1063/1.1787832", &
"Force calculations in CONQUEST")

      call bib%add_ref( &
"Carbogno2017", &
"C. Carbogno and R. Ramprasad and M. Scheffler", &
"Ab Initio Green-Kubo approach for thermal conductivity of solids", &
"Phys. Rev. Lett.", &
118, &
175901, &
2017, &
"10.1103/PhysRevLett.118.175901", &
"Method for Green-Kubo thermoal conductivity calculations for solids")

      call bib%add_ref( &
"Pfrommer1997", &
"B. G. Pfrommer and M. Cote and S. G. Louie and M. L. Cohen", &
"Relaxation of crystals with the quasi-Newton method", &
"J. Comput. Phys.", &
131, &
233, &
1997, &
"10.1006/jcph.1996.5612", &
"Method for simultaneous cell and geometry optimisation")

    end subroutine get_bib_db
    !!***

    !!****f* references/compile_biblio *
    !!
    !!  NAME 
    !!   compile_biblio
    !!  PURPOSE
    !!   Check the relevant input flags, and compile a bibliography
    !!  AUTHOR
    !!   Zamaan Raza
    !!  CREATION DATE
    !!   2019/07/04
    !!  MODIFICATION HISTORY
    !!
    !!  SOURCE
    !!
    subroutine compile_biblio

      use global_module, only: flag_diagonalisation, flag_Multisite, &
                               flag_XLBOMD, flag_basis_set, PAOs, &
                               flag_heat_flux, optcell_method
      use pseudopotential_common, only: pseudo_type, ABINIT
      use input_module,  only: leqi
      use control,       only: md_ensemble
      use md_control,    only: md_thermo_type

      ! local variables
      type(type_bibliography) :: bib
      type(type_reference)    :: reference

      call bib%init_bib
      call get_bib_db(bib)

      if (inode==ionode) then
        write(io_lun,*)
        write(io_lun,'(2x,a)') "BIBLIOGRAPHY: If you publish results obtained &
                               &with CONQUEST, please consider citing the &
                               &following:"
        write(io_lun,*)
      end if

      ! Cite publications for *any* CONQUEST calculation
      call bib%cite("Gillan2007")
      call bib%cite("Bowler2010")

      if (inode==ionode) then
        write(io_lun,'(2x,a)') "The following papers detail methodology used &
                               &in this CONQUEST calculation:"
        write(io_lun,*)
      end if

      ! Cite method-specific publications
      if (flag_Multisite)               call bib%cite("Nakata2014")
      if (flag_basis_set == PAOs)       call bib%cite("Torralba2008")
      if (pseudo_type == ABINIT)        call bib%cite("Hamann2013")
      if (pseudo_type == ABINIT)        call bib%cite("Bowler2019")
      if (leqi(md_thermo_type, 'nhc'))  call bib%cite("Martyna1996")
      if (leqi(md_thermo_type, 'nhc'))  call bib%cite("Hirakawa2017")
      if (leqi(md_ensemble, 'npt'))     call bib%cite("Martyna1994")
      if (leqi(md_ensemble, 'npt'))     call bib%cite("Shinoda2004")
      if (flag_XLBOMD)                  call bib%cite("Arita2014")
      if (flag_XLBOMD)                  call bib%cite("Niklasson2006")
      if (flag_XLBOMD)                  call bib%cite("Niklasson2009")
      if (flag_heat_flux)               call bib%cite("Carbogno2017")
      if (optcell_method == 3)          call bib%cite("Pfrommer1997")

      call bib%close_bib

    end subroutine compile_biblio
    !!***

end module references
