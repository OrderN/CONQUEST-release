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
    !!   2020/01/03 10:11 dave
    !!    Arrange references into Conquest-specific and other
    !!    and order by date
    !!   2020/01/03 14:09 dave
    !!    Add further papers, including XC functional references
    !!  SOURCE
    !!
    subroutine get_bib_db(bib)

      ! Passed variables
      type(type_bibliography), intent(inout) :: bib

      ! Conquest-specific papers
      call bib%add_ref( &
           "Hernandez1997", &
           "{Hern\'andez, E. and Gillan, M. J. and Goringe, C. M.}", &
           "{Basis functions for linear-scaling first-principles calculations}", &
           "{Phys. Rev. B}", &
           55, &
           13485, &
           1997, &
           "10.1103/PhysRevB.55.13485", &
           "Blip functions in CONQUEST")

      call bib%add_ref( &
           "Bowler:1998lo", &
           "{Bowler, D R and Gillan, M J}", &
           "{Length-scale ill conditioning in linear-scaling DFT}", &
           "{Comp. Phys. Commun.}", &
           112, &
           103, &
           1998, &
           "10.1016/S0010-4655(98)00061-7", &
           "Optimisation of blips in CONQUEST")

      call bib%add_ref( &
           "Bowler:1999if", &
           "{Bowler, D R and Gillan, M J}", &
           "{Density matrices in O(N) electronic structure calculations: &
           &theory and applications}", &
           "{Comp. Phys. Commun.}", &
           120, &
           95, &
           1999, &
           "10.1016/S0010-4655(99)00221-0", &
           "CONQUEST O(N) implementation")

      call bib%add_ref( &
           "Bowler2001", &
           "{D. R. Bowler and T. Miyazaki and M. J. Gillan}", &
           "{Parallel sparse matrix multiplication for linear scaling &
           &electronic structure calculations}", &
           "{Comp. Phys. Commun.}", &
           137, &
           255, &
           2001, &
           "10.1016/S0010-4655(01)00164-3", &
           "O(N) matrix multiplication implementation")

      
      call bib%add_ref(&
           "Bowler2002pt", &
           "{Bowler, D R and Miyazaki, T and Gillan, M J}", &
           "{Recent progress in linear scaling $\backslash$textit{\{}ab initio{\}} &
           &electronic structure techniques}", &
           "{J. Phys. Condens. Matter}", &
           14, &
           2781, &
           2002, &
           "10.1088/0953-8984/14/11/303",&
           "Overview of CONQUEST methodology")

      call bib%add_ref( &
           "Miyazaki2004", &
           "{T. Miyazaki and D. R. Bowler and R. Choudhury and M. J. Gillan}", &
           "{Atomic force algorithms in density functional theory &
           &electronic-structure techniques based on local orbitals}", &
           "{J. Chem. Phys.}", &
           121, &
           6186, &
           2004, &
           "10.1063/1.1787832", &
           "Force calculations in CONQUEST")

      call bib%add_ref( &
           "Bowler:2006xr", &
           "{Bowler, D R and Choudhury, R and Gillan, M J and Miyazaki, T}", &
           "{Recent progress with large-scale ab initio calculations: &
           &the $\backslash$textsc{\{}CONQUEST{\}} code}", &
           "{phys. stat. sol. b}", &
           243, &
           989, &
           2006, &
           "10.1002/pssb.200541386", &
           "Overview of CONQUEST methodology, including diagonalisation")

      call bib%add_ref( &
           "Torralba2008", &
           "{A. S. Torralba and M. Todorovic and V. Brazdova and R. Choudhury &
           &and T. Miyazaki and M. J. Gillan and D. R. Bowler}", &
           "{Pseudo-atomic orbitals as basis sets for O(N) DFT code CONQUEST}", &
           "{J. Phys.: Condens. Matter}", &
           20, &
           294206, &
           2008, &
           "10.1088/0953-8984/20/29/294206", &
           "PAO basis sets for CONQUEST")

      call bib%add_ref( &
           "Bowler2010", &
           "{D. R. Bowler and T. Miyazaki}", &
           "{Calculations for millions of atoms with density functional theory: &
           &linear scaling shows its potential}", &
           "{J. Phys.: Condens. Matter}", &
           22, &
           074207, &
           2010, &
           "10.1088/0953-8984/22/7/074207", &
           "Further details on O(N) methodology in CONQUEST")

      call bib%add_ref( &
           "Sena:2011fc", &
           "{Sena, Alex M P and Miyazaki, Tsuyoshi and Bowler, David R}", &
           "{Linear Scaling Constrained Density Functional Theory in CONQUEST}", &
           "{J. Chem. Theory Comput.}", &
           4, &
           884, &
           2011, &
           "10.1021/ct100601n", &
           "Constrained DFT in CONQUEST")

      call bib%add_ref( &
           "Terranova:2013fk", &
           "{Terranova, U and Bowler, D R}", &
           "{$\Delta$ Self-Consistent Field Method for Natural Anthocyanidin Dyes}", &
           "{J. Chem. Theory Comput.}", &
           9, &
           3181, &
           2013, &
           "10.1021/ct400356k", &
           "Implementation of delta-SCF method in CONQUEST")

      call bib%add_ref( &
           "Arita2014", &
           "{M. Arita and D. R. Bowler and T. Miyazaki}", &
           "{Stable and efficient linear scaling first-principles molecular dynamics &
           &for 10000+ atoms}", &
           "{J. Chem. Theor. Comput.}", &
           10, &
           5419, &
           2014, &
           "10.1021/ct500847y", &
           "CONQUEST XL-BOMD implementation")

      call bib%add_ref( &
           "Nakata2014", &
           "{A. Nakata and D. R. Bowler and T. Miyazaki}", &
           "{Efficient calculations with multisite local orbitals &
           &in a large scale DFT code CONQUEST}", &
           "{J. Chem. Theor. Comput.}", &
           10, &
           4813, &
           2014, &
           "10.1021/ct5004934", &
           "CONQUEST multisite support functions")


      call bib%add_ref( &
           "Nakata2015", &
           "{Nakata, Ayako and Bowler, David and Miyazaki, Tsuyoshi}", &
           "{Optimized multi-site local orbitals in the large-scale DFT program CONQUEST}", &
           "{Phys. Chem. Chem. Phys.}", &
           17, &
           31427, &
           2015, &
           "10.1039/C5CP00934K", &
           "Optimisation of CONQUEST multisite support functions")

      call bib%add_ref( &
           "Hirakawa2017", &
           "{T. Hirakawa and T. Suzuki and D. R. Bowler and T. Miyazaki}", &
           "{Canonical-ensemble extended Lagrangian Born-Oppenheimer molecular dynamics &
           &for the linear scaling density functional theory}", &
           "{J. Phys. Condens. Matter}", &
           29, &
           405901, &
           2017, &
           "10.1088/1361-648X/aa810d", &
           "CONQUEST canonical ensemble MD with XL-BOMD")

      call bib%add_ref( &
           "Bowler2019", &
           "{D. R. Bowler and J. S. Baker and J. T. L. Poulton and S. Y. Mujahed &
           &and J. Lin and S. Yadav and Z. Raza and T. Miyazaki}", &
           "{Highly accurate local basis sets for large-scale DFT calculations in CONQUEST}", &
           "{Japan. J. Appl. Phys.}", &
           58, &
           100503, &
           2019, &
           "10.7567/1347-4065/ab45af", &
           "Details of new CONQUEST PAOs from basis generation tool, &
           &using Hamann-type pseudopotentials")

      ! Other papers: not Conquest specific
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
           "Pfrommer1997", &
           "B. G. Pfrommer and M. Cote and S. G. Louie and M. L. Cohen", &
           "Relaxation of crystals with the quasi-Newton method", &
           "J. Comput. Phys.", &
           131, &
           233, &
           1997, &
           "10.1006/jcph.1996.5612", &
           "Method for simultaneous cell and geometry optimisation")

      call bib%add_ref( &
           "Shinoda2004", &
           "W. Shinoda and M. Shiga and M. Mikami", &
           "Rapid estimation of elastic constants by molecular dynamics simulation&
           &under constant stress", &
           "Phys. Rev. B", &
           69, &
           134103, &
           2004, &
           "10.1103/PhysRevB.69.134103", &
           "Liouvillian splitting for NPT MD used in CONQUEST")

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
           "Carbogno2017", &
           "C. Carbogno and R. Ramprasad and M. Scheffler", &
           "Ab Initio Green-Kubo approach for thermal conductivity of solids", &
           "Phys. Rev. Lett.", &
           118, &
           175901, &
           2017, &
           "10.1103/PhysRevLett.118.175901", &
           "Method for Green-Kubo thermoal conductivity calculations for solids")

      ! Functionals
      call bib%add_ref( &
           "Perdew1981", &
           "{Perdew, J. P. and Zunger, Alex}", &
           "{Self-interaction correction to density-functional approximations &
           &for many-electron systems}", &
           "{Phys. Rev. B}", &
           23, &
           5048, &
           1981, &
           "10.1103/PhysRevB.23.5048", &
           "PZ81 LDA parameterisation")
      call bib%add_ref( &
           "Perdew1992", &
           "{Perdew, John P. and Wang, Yue}", &
           "{Accurate and simple analytic representation of the electron-gas &
           &correlation energy}", &
           "{Phys. Rev. B}", &
           45, &
           13244, &
           1992, &
           "10.1103/PhysRevB.45.13244", &
           "PW92 LDA parameterisation")
      call bib%add_ref( &
           "Perdew1996", &
           "{Perdew, John P. and Burke, Kieron and Ernzerhof, Matthias}", &
           "{Generalized Gradient Approximation Made Simple}", &
           "{Phys. Rev. Lett.}", &
           77, &
           3865, &
           1996, &
           "10.1103/PhysRevLett.77.3865", &
           "PBE GGA XC functional")
      call bib%add_ref( &
           "Zhang1998", &
           "{Comment on ``Generalized Gradient Approximation Made Simple''}", &
           "{Zhang, Yingkai and Yang, Weitao}", &
           "{Phys. Rev. Lett.}", &
           80, &
           890, &
           1998, &
           "10.1103/PhysRevLett.80.890", &
           "PBE GGA functional revised by Zhang & Yang 1998 (revPBE)")
      call bib%add_ref( &
           "Hammer1999", &
           "{Improved adsorption energetics within density-functional theory using &
           &revised Perdew-Burke-Ernzerhof functionals}", &
           "{Hammer, B. and Hansen, L. B. and N\o{}rskov, J. K.}", &
           "{Phys. Rev. B}", &
           59, &
           7413, &
           1999, &
           "10.1103/PhysRevB.59.7413", &
           "PBE GGA functional revised by Hammer et al (rPBE)")
      call bib%add_ref( &
           "Wu2006", &
           "{More accurate generalized gradient approximation for solids}", &
           "{Wu, Zhigang and Cohen, R. E.}", &
           "{Phys. Rev. B}", &
           73, &
           235116, &
           2006, &
           "10.1103/PhysRevB.73.235116", &
           "Wu-Cohen reformulation of PBE GGA functional")
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
    !!   2020/03/24 14:22 dave
    !!    Changed to write reference key only if iprint<2
    !!
    !!  SOURCE
    !!
    subroutine compile_biblio

      use global_module, only: flag_diagonalisation, flag_Multisite, &
                               flag_XLBOMD, flag_basis_set, PAOs, &
                               flag_heat_flux, optcell_method, &
                               flag_perform_cDFT, flag_DeltaSCF, runtype, iprint_init
      use pseudopotential_common, only: pseudo_type, ABINIT
      use input_module,  only: leqi
      use control,       only: md_ensemble
      use md_control,    only: md_thermo_type
      use XC,            only: flag_functional_type, functional_lda_pz81, &
           functional_lda_gth96, functional_lda_pw92, functional_gga_pbe96, &
           functional_gga_pbe96_rev98, functional_gga_pbe96_r99, functional_gga_pbe96_wc, &
           write_xc_refs

      ! local variables
      type(type_bibliography) :: bib
      type(type_reference)    :: reference

      call bib%init_bib
      call get_bib_db(bib)

      if (inode==ionode) then
         if(iprint_init>1) then
            write(io_lun,*)
            write(io_lun,'(4x,a)') "BIBLIOGRAPHY: If you publish results obtained &
                 &with CONQUEST, please consider citing the &
                 &following:"
            write(io_lun,*)
         else
            write(io_lun,fmt='(/4x,a/)') "BIBLIOGRAPHY: Please consider citing the following &
                 &references in the conquest.bib file"
         end if
      end if

      ! Cite publications for *any* CONQUEST calculation
      call bib%cite("Bowler2002pt",punc=", ",pre="CONQUEST: ")
      call bib%cite("Miyazaki2004",punc=", ")
      call bib%cite("Bowler:2006xr") ! Replace with 2020 JCP when submitted
      write(io_lun,*)

      !if (inode==ionode) then
      !  write(io_lun,'(/4x,a)') "The following papers detail methodology used &
      !                         &in this CONQUEST calculation:"
      !  write(io_lun,*)
      !end if

      ! Cite method-specific publications
      ! Basis sets
      if (flag_basis_set == PAOs) then
         call bib%cite("Bowler2019",punc=", ",pre="Basis: ")
         if (flag_Multisite) then
            call bib%cite("Nakata2014",punc=", ")
            call bib%cite("Nakata2015")
         end if
      else
         call bib%cite("Hernandez1997",punc=", ",pre="Basis: ")
         call bib%cite("Bowler:1998lo")
      end if
      write(io_lun,*)
      ! Finding DM
      if (.NOT.flag_diagonalisation) then ! O(N)
         call bib%cite("Bowler:1999if",punc=", ",pre="DM: ")
         call bib%cite("Bowler2001",punc=", ")
         call bib%cite("Bowler2010")
      else
         call bib%cite("Bowler:2006xr",pre="DM: ")
      end if
      write(io_lun,*)
      if(flag_perform_cDFT)             call bib%cite("Sena:2011fc")
      if(flag_DeltaSCF)                 call bib%cite("Terranova:2013fk")
      ! Pseudopotentials
      if (pseudo_type == ABINIT) then
         call bib%cite("Hamann2013",punc=", ",pre="Pseudopotentials: ")
         call bib%cite("Bowler2019")
      end if
      write(io_lun,*)
      ! MD
      if (leqi(runtype,'md')) then
         call bib%cite("Arita2014",punc=", ",pre="MD: ")
         if (flag_XLBOMD)                  call bib%cite("Niklasson2006",punc=", ")
         if (flag_XLBOMD)                  call bib%cite("Niklasson2009",punc=", ")
         if (leqi(md_ensemble, 'npt')) then
            if (leqi(md_thermo_type, 'nhc')) then
               call bib%cite("Martyna1996",punc=", ")
               call bib%cite("Hirakawa2017",punc=", ")
            end if
            call bib%cite("Martyna1994",punc=", ")
            call bib%cite("Shinoda2004")
         else
            if (leqi(md_thermo_type, 'nhc')) then
               call bib%cite("Martyna1996",punc=", ")
               call bib%cite("Hirakawa2017")
            end if
         end if
         !if (flag_heat_flux)               call bib%cite("Carbogno2017")
         write(io_lun,*)
      end if
      ! Cell optimisation
      if (optcell_method == 3) then
         call bib%cite("Pfrommer1997",pre="Cell optimisation: ")
         write(io_lun,*)
      end if
      ! Functionals
      if(flag_functional_type>0) then
         !if (inode==ionode) then
         !   write(io_lun,'(/4x,a)') "The following paper details the functional used &
         !        &in this CONQUEST calculation:"
         !   write(io_lun,*)
         !end if
         if(flag_functional_type == functional_lda_pz81) &
              call bib%cite("Perdew1981",pre="XC functional: ")
         if(flag_functional_type == functional_lda_pw92) &
              call bib%cite("Perdew1992",pre="XC functional: ")
         if(flag_functional_type == functional_gga_pbe96) &
              call bib%cite("Perdew1996",pre="XC functional: ")
         if(flag_functional_type == functional_gga_pbe96_rev98) &
              call bib%cite("Zhang1998",pre="XC functional: ")
         if(flag_functional_type == functional_gga_pbe96_r99) &
              call bib%cite("Hammer1999",pre="XC functional: ")
         if(flag_functional_type == functional_gga_pbe96_wc) &
              call bib%cite("Wu2006",pre="XC functional: ")
         write(io_lun,*)
      else ! LibXC
         call write_xc_refs
      end if
      
      call bib%close_bib

    end subroutine compile_biblio
    !!***

end module references
