module references

  use biblio

  implicit none

  contains

    subroutine biblio_data(bib)

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

    end subroutine biblio_data
    !!***

    !!****f* biblio/compile_biblio *
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
                               flag_XLBOMD
      use input_module,  only: leqi
      use control,       only: md_ensemble
      use md_control,    only: md_thermo_type

      ! local variables
      type(type_bibliography) :: bib
      type(type_reference)    :: reference

      call bib%init_bib
      call biblio_data(bib)

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
      if (leqi(md_thermo_type, 'nhc'))  call bib%cite("Hirakawa2017")
      if (flag_XLBOMD)                  call bib%cite("Arita2014")
      if (flag_Multisite)               call bib%cite("Nakata2014")

      write(io_lun,*)
      call bib%close_bib

    end subroutine compile_biblio
    !!***

end module references
