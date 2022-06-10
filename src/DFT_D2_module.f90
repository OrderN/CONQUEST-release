! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id: DFT_D2_module.f90 85 2011-06-09 14:46:30Z M. Arita $
! ------------------------------------------------------------------------------
! Module DFT_D2
!-------------------------------------------------------------------------------
! Code area 9: general
!-------------------------------------------------------------------------------

!!****h* Conquest_DFT_D2 *
!! NAME
!!  DFT_D2
!! PURPOSE
!!  Stores various useful dimensions and a subroutine
!! USES
!!  
!! AUTHOR
!!  M. Arita
!! CREATION DATE
!!  11/06/09 M. Arita
!! MODIFICATION HISTORY
!!  2022/06/09 12:27 dave
!!   Moved disp_energy to energy module, partition_distance taken
!!   from ion_electrostatic
!! SOURCE
module DFT_D2

  use datatypes
  use energy, only: disp_energy

  implicit none

  save

  real(double), allocatable :: disp_force(:, :)
  real(double), allocatable :: atm_disp_coeff(:), vdW_rad(:)
  real(double), dimension(3,3) :: disp_stress

contains


! --------------------------------------------------------------
! Subroutine dispersion_D2
! --------------------------------------------------------------

!!****f* DFT_D2/dispersion_D2 *
!!
!! NAME
!!  dispersion_D2
!!  USAGE
!!
!!  PURPOSE
!!    Compute the dispersion interaction and the atomic fores
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!    M. Arita
!!  CREATION DATE
!!    11/06/10
!!  MODIFICATION
!!    2022/06/08 16:11 dave
!!     Adding stress calculation
!!  SOURCE
!!
  subroutine dispersion_D2

    use cover_module, ONLY: D2_CS, make_CS
    use datatypes
    use units
    use dimens, ONLY: atomicnum, r_dft_d2, r_super_x, r_super_y, r_super_z
    use GenComms, ONLY: cq_abort, gcopy, inode, ionode, gsum, &
         my_barrier
    use group_module, ONLY: parts
    use global_module, ONLY: io_lun, numprocs, &
         ni_in_cell, species_glob, id_glob, flag_only_dispersion, &
         flag_stress, flag_full_stress
    use numbers
    use primary_module, ONLY: bundle
    use io_module, ONLY: get_file_name
    use species_module, ONLY: n_species, species_label
    use input_module, ONLY: io_assign, io_close
    use XC, ONLY: s_6
    use ion_electrostatic, ONLY: partition_distance

    implicit none

    ! local variables
    integer :: i, j, ip, ia, jp, ja, ia_glob
    integer :: stat, Z_i, Z_j
    integer :: nppx, nppy, nppz, ind_part, nccx, nccy, nccz, m1, m2, &
         m3, m, ig_atom_beg
    integer, parameter :: d = 20.0_double
    integer, allocatable :: neighbour_part(:)
    real(double) :: real_cell_vec(3, 3), part_cell_vec(3, 3), v1,&
         v2, v3, rr(3), part_cell_dual(3, 3), distance, x_ia, y_ia, &
         z_ia, x_ja, y_ja, z_ja, rijx, rijy, rijz, rij, rij2, rij6, &
         power, f_damp, denominator, emp_energy, emp_force, &
         unit_vec_x, unit_vec_y, unit_vec_z
    real(double) :: C_i, C_j, C_ij, r_vdW_i, r_vdW_j, R_r

    ! Initialisation of the dispersion energy and the atomic forces
    disp_energy = 0.0_double
    disp_force  = 0.0_double
    if(flag_stress) disp_stress = zero

    allocate ( neighbour_part(D2_CS%ng_cover), STAT=stat )
    if (stat .NE. 0) call cq_abort("Eror allocating neighbour_list: ", D2_CS%ng_cover)

    real_cell_vec(1, 1) = r_super_x
    real_cell_vec(1, 2) = zero
    real_cell_vec(1, 3) = zero
    real_cell_vec(2, 1) = zero
    real_cell_vec(2, 2) = r_super_y
    real_cell_vec(2, 3) = zero
    real_cell_vec(3, 1) = zero
    real_cell_vec(3, 2) = zero
    real_cell_vec(3, 3) = r_super_z

    ! primitive translation vectors for partitions
    do i = 1, 3
      part_cell_vec(1, i) = real_cell_vec(1, i) / real(parts%ngcellx, double)
      part_cell_vec(2, i) = real_cell_vec(2, i) / real(parts%ngcelly, double)
      part_cell_vec(3, i) = real_cell_vec(3, i) / real(parts%ngcellz, double)
    enddo

    ! Initialise the partition distance routine
    call partition_distance(.true., part_cell_vec, part_cell_dual)

    call my_barrier()

    ! Kernel of the DFT-D2 scheme
    do ip = 1, bundle%groups_on_node
      ! The address of the curret partition "ip"
      nppx = 1 + bundle%idisp_primx(ip) + D2_CS%nspanlx
      nppy = 1 + bundle%idisp_primy(ip) + D2_CS%nspanly
      nppz = 1 + bundle%idisp_primz(ip) + D2_CS%nspanlz
      ! Initialisation of the neighbour_list
      ! 0 -> not neighbour, 1 -> neighbour
      neighbour_part = 0
      do jp = 1, D2_CS%ng_cover
        ind_part = D2_CS%lab_cover(jp) ! PARTITION INDEX (CC)
        ! The address of the partition "ind_part" in CS
        nccx = 1 + (ind_part - 1) / (D2_CS%ncovery * D2_CS%ncoverz)
        nccy = 1 + (ind_part - 1 - (nccx - 1) * D2_CS%ncovery * D2_CS%ncoverz) / D2_CS%ncoverz
        nccz = ind_part - (nccx - 1) * D2_CS%ncovery * D2_CS%ncoverz - (nccy - 1) * D2_CS%ncoverz
        ! Calculate the interpartition distance between "ip" and "jp"
        m1 = nccx - nppx
        m2 = nccy - nppy
        m3 = nccz - nppz
        m = m1*m1 + m2*m2 + m3*m3
        if (m .NE. 0) then ! IF "jp" IS NOT "ip"
          v1 = real(m1, double)
          v2 = real(m2, double)
          v3 = real(m3, double)
          rr(1) = half * (v1*part_cell_vec(1, 1) + v2*part_cell_vec(2, 1) + v3*part_cell_vec(3, 1))
          rr(2) = half * (v1*part_cell_vec(1, 2) + v2*part_cell_vec(2, 2) + v3*part_cell_vec(3, 2))
          rr(3) = half * (v1*part_cell_vec(1, 3) + v2*part_cell_vec(2, 3) + v3*part_cell_vec(3, 3))
          call partition_distance(.false., part_cell_vec, part_cell_dual, rr, distance)
          distance = distance * two ! INTERPARTITION DISTANCE
          ! See which partitions are neighbours of the partition "ip"
          if (distance .LT. r_dft_d2) then
            neighbour_part(jp) = 1
          endif !(if, distance)
        elseif (m .EQ. 0) then
          neighbour_part(jp) = 1
        endif !(m)
      enddo !(jp, D2_CS%ng_cover)
      if (bundle%nm_nodgroup(ip) .GT. 0) then
        do ia = 1, bundle%nm_nodgroup(ip) ! LOOP OVER ATOMS IN PARTITION "ip"
          ! Coordinates of the atom "ia"
          x_ia = bundle%xprim(bundle%nm_nodbeg(ip) + ia - 1)
          y_ia = bundle%yprim(bundle%nm_nodbeg(ip) + ia - 1)
          z_ia = bundle%zprim(bundle%nm_nodbeg(ip) + ia - 1)
          do jp = 1, D2_CS%ng_cover ! LOOP OVER PARTITIONS OF THE CS
            ! Choose neighbour partitions
            !OLD if (neighbour_part(jp) .EQ. 1) then
            if ( (neighbour_part(jp) .EQ. 1) .AND. (D2_CS%n_ing_cover(jp) .GT. 0) ) then
              do ja = 1, D2_CS%n_ing_cover(jp)
                x_ja = D2_CS%xcover(D2_CS%icover_ibeg(jp) + ja - 1)
                y_ja = D2_CS%ycover(D2_CS%icover_ibeg(jp) + ja - 1)
                z_ja = D2_CS%zcover(D2_CS%icover_ibeg(jp) + ja - 1)
                rijx = x_ia - x_ja
                rijy = y_ia - y_ja
                rijz = z_ia - z_ja
                rij2 = rijx*rijx + rijy*rijy + rijz*rijz
                ! If the interatomic distance is smaller than the cutoff, culculate the dispersion,
                ! but exclude the situation in which "ia" = "ja".
                if ( (rij2 .GT. RD_ERR) .AND. (rij2 .LE. r_dft_d2*r_dft_d2) ) then
                  ig_atom_beg = parts%icell_beg(D2_CS%lab_cell(jp))
                  ! Calculate the molecular dispersion coefficient and vdW radius
                  Z_i = atomicnum(bundle%species(bundle%nm_nodbeg(ip) + ia - 1))
                  Z_j = atomicnum(species_glob(id_glob(ig_atom_beg+ja-1)))
                  C_i = atm_disp_coeff(Z_i)
                  C_j = atm_disp_coeff(Z_j)
                  ! Following J Comput Chem 27, 1787
                  C_ij = sqrt(C_i * C_j)
                  r_vdW_i = vdW_rad(Z_i)
                  r_vdW_j = vdW_rad(Z_j)
                  R_r = r_vdW_i + r_vdW_j
                  ! Compute the damping function
                  rij = sqrt(rij2)
                  power = rij / R_r - 1.0_double
                  power = -d * power
                  denominator = 1.0_double + exp(power)
                  f_damp = 1.0_double / denominator
                  rij6 = rij2 * rij2 * rij2
                  ! Calculate the dispersion energy
                  emp_energy = C_ij * f_damp / rij6
                  disp_energy = disp_energy + emp_energy
                  ! Calculate the dispersion force
                  unit_vec_x = rijx / rij
                  unit_vec_y = rijy / rij
                  unit_vec_z = rijz / rij
                  !OLD emp_force = 6.0_double / rij - d / R_r * f_damp * exp(power)
                  emp_force = d / R_r * f_damp * exp(power) - 6.0_double / rij
                  ia_glob = bundle%ig_prim(bundle%nm_nodbeg(ip) + ia - 1)
                  disp_force(1, ia_glob) = disp_force(1, ia_glob) + emp_energy * emp_force * unit_vec_x
                  disp_force(2, ia_glob) = disp_force(2, ia_glob) + emp_energy * emp_force * unit_vec_y
                  disp_force(3, ia_glob) = disp_force(3, ia_glob) + emp_energy * emp_force * unit_vec_z
                  if(flag_stress) then
                     if(flag_full_stress) then
                        ! NB force direction is set by unit_vec so tracks first index
                        disp_stress(1,1) = disp_stress(1,1) + emp_energy * emp_force &
                             * rijx * unit_vec_x
                        disp_stress(1,2) = disp_stress(1,2) + emp_energy * emp_force &
                             * rijy * unit_vec_x
                        disp_stress(1,3) = disp_stress(1,3) + emp_energy * emp_force &
                             * rijz * unit_vec_x
                        disp_stress(2,1) = disp_stress(2,1) + emp_energy * emp_force &
                             * rijx * unit_vec_y
                        disp_stress(2,2) = disp_stress(2,2) + emp_energy * emp_force &
                             * rijy * unit_vec_y
                        disp_stress(2,3) = disp_stress(2,3) + emp_energy * emp_force &
                             * rijz * unit_vec_y
                        disp_stress(3,1) = disp_stress(3,1) + emp_energy * emp_force &
                             * rijx * unit_vec_z
                        disp_stress(3,2) = disp_stress(3,2) + emp_energy * emp_force &
                             * rijy * unit_vec_z
                        disp_stress(3,3) = disp_stress(3,3) + emp_energy * emp_force &
                             * rijz * unit_vec_z
                     else
                        disp_stress(1,1) = disp_stress(1,1) + emp_energy * emp_force &
                             * rijx * unit_vec_x
                        disp_stress(2,2) = disp_stress(2,2) + emp_energy * emp_force &
                             * rijy * unit_vec_y
                        disp_stress(3,3) = disp_stress(3,3) + emp_energy * emp_force &
                             * rijz * unit_vec_z
                     end if
                  end if
                  !OLD disp_force(1, ia) = disp_force(1, ia) + emp_energy * emp_force * unit_vec_x
                  !OLD disp_force(2, ia) = disp_force(2, ia) + emp_energy * emp_force * unit_vec_y
                  !OLD disp_force(3, ia) = disp_force(3, ia) + emp_energy * emp_force * unit_vec_z
                endif !(if, cutoff)
              enddo !(do, ja)
            else
              cycle
            endif !(if, neighbour)
          enddo !(jp, D2_CS%ng_cover)
        enddo !(ia, bundle%nm_groups_on_node)
      endif !(bundle%nm_nodgroup(ip) > 0)
    enddo !(ip)

    call gsum(disp_energy)
    call gsum(disp_force, 3, ni_in_cell)
    disp_energy = - disp_energy * s_6 * half ! CANCEL THE DOUBLE-COUNTING
    disp_force = s_6 * disp_force
    if(flag_stress) then
       call gsum(disp_stress,3,3)
       disp_stress = -half * s_6 * disp_stress  ! Double counting again
    end if

    ! Works only when you want to get the dispersion
    if (flag_only_dispersion) then
      if (inode == ionode) then
        write (io_lun, '(/a)') "***** Dispersion Information: DFT-D2 *****"
        write (io_lun,16) en_conv*disp_energy, en_units(energy_units)
        write (io_lun,fmt='(/,20x,"Disp. Force on atoms (",a2,"/",a2,")"/)') en_units(energy_units), d_units(dist_units)
        write (io_lun,fmt='("Atom        X               Y               Z")')
        do i = 1, ni_in_cell
            write (io_lun, 109) i, (for_conv*disp_force(j,i),j=1,3)
         enddo
         write(io_lun,'(/8x,"Dispersion stress: ",3f15.8)') disp_stress
      endif
    endif

16  format('E_disp:',f25.15,' ',a2)
101 format('Force on atom ', i9)
109 format(i10, 3f15.10)

    deallocate ( neighbour_part, STAT=stat )
    if (stat .NE. 0) call cq_abort("Eror deallocating neighbour_part: ", D2_CS%ng_cover, stat)

    return
  end subroutine dispersion_D2
!!***

! --------------------------------------------------------------
! Subroutine set_para_D2
! --------------------------------------------------------------

!!****f* DFT_D2/set_para_D2 *
!!
!!  NAME
!!    set_para_D2
!!  USAGE
!!
!!  PURPOSE
!!    Reads and stores the necessary parameters in the DFT-D2 scheme
!!    The coefficients are from J. Comp. Chem. 27, 15, 2006 up to Xe
!!    and from J. Chem. Phys, 132 (2010), 154104 (DFT-D3) thereafter
!!    Change of units via translation parameters
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!    M. Arita
!!  CREATION DATE
!!    11/06/09
!!  MODIFICATION HISTORY
!!    2011/10/06 14:01 dave
!!     Bug fixes
!!    2022/06/09 08:31 dave
!!     Hard-code parameters (not read from file)
!!  SOURCE
!!
  subroutine set_para_D2

    use datatypes
    use numbers
    use units, ONLY: AngToBohr, eVToJ, HaToeV, AvogadroC
    use GenComms, ONLY: cq_abort, gcopy, inode, ionode, my_barrier
    use global_module, ONLY: io_lun, numprocs, ni_in_cell
    use io_module, ONLY: get_file_name
    use species_module, ONLY: n_species, species_label
    use input_module, ONLY: io_assign, io_close

    implicit none

    ! local variables
    integer :: i, j, stat, ios, lun
    real(double) :: coeff, radius
    real(double) :: trans1 ! [Jnm^6/mol]-->[Ha bohr^6]
    real(double) :: trans2 ! [ang]-->[bohr]
    character(len=20) dummy
    character(len=2)  kind
    integer, parameter :: n_elements = 86

    ! debug variables
    integer :: lun_para
    character(len=30) :: filename

    trans1 = (ten*AngToBohr)**6/(eVToJ*HaToeV*AvogadroC)
    trans2 = AngToBohr
    ! Allocate storage
    allocate( atm_disp_coeff(n_elements), vdW_rad(n_elements), STAT=stat )
    if (stat .NE. 0) call cq_abort("Error allocating DFT-D2 parameters: ",n_elements)
    ! Set coefficients
    ! Atomic dispersion coefficients, C_i
    atm_disp_coeff = zero
    atm_disp_coeff(1) = 0.14  * trans1   ! H
    atm_disp_coeff(2) = 0.08  * trans1

    atm_disp_coeff(3 ) = 1.61  * trans1  ! Li
    atm_disp_coeff(4 ) = 1.61  * trans1
    atm_disp_coeff(5 ) = 3.13  * trans1
    atm_disp_coeff(6 ) = 1.75  * trans1
    atm_disp_coeff(7 ) = 1.23  * trans1
    atm_disp_coeff(8 ) = 0.70  * trans1
    atm_disp_coeff(9 ) = 0.75  * trans1   
    atm_disp_coeff(10) = 0.63  * trans1
    
    atm_disp_coeff(11) = 5.71  * trans1  ! Na
    atm_disp_coeff(12) = 5.71  * trans1
    atm_disp_coeff(13) = 10.79 * trans1
    atm_disp_coeff(14) = 9.23  * trans1
    atm_disp_coeff(15) = 7.84  * trans1
    atm_disp_coeff(16) = 5.57  * trans1
    atm_disp_coeff(17) = 5.07  * trans1
    atm_disp_coeff(18) = 4.61  * trans1
    
    atm_disp_coeff(19) = 10.80 * trans1  ! K
    atm_disp_coeff(20) = 10.80 * trans1
    atm_disp_coeff(21) = 10.80 * trans1
    atm_disp_coeff(22) = 10.80 * trans1
    atm_disp_coeff(23) = 10.80 * trans1
    atm_disp_coeff(24) = 10.80 * trans1
    atm_disp_coeff(25) = 10.80 * trans1
    atm_disp_coeff(26) = 10.80 * trans1
    atm_disp_coeff(27) = 10.80 * trans1
    atm_disp_coeff(28) = 10.80 * trans1
    atm_disp_coeff(29) = 10.80 * trans1
    atm_disp_coeff(30) = 10.80 * trans1
    atm_disp_coeff(31) = 16.99 * trans1
    atm_disp_coeff(32) = 17.10 * trans1
    atm_disp_coeff(33) = 16.37 * trans1
    atm_disp_coeff(34) = 12.64 * trans1
    atm_disp_coeff(35) = 12.47 * trans1
    atm_disp_coeff(36) = 12.01 * trans1
    
    atm_disp_coeff(37) = 24.67 * trans1  ! Rb
    atm_disp_coeff(38) = 24.67 * trans1
    atm_disp_coeff(39) = 24.67 * trans1
    atm_disp_coeff(40) = 24.67 * trans1
    atm_disp_coeff(41) = 24.67 * trans1
    atm_disp_coeff(42) = 24.67 * trans1
    atm_disp_coeff(43) = 24.67 * trans1
    atm_disp_coeff(44) = 24.67 * trans1
    atm_disp_coeff(45) = 24.67 * trans1
    atm_disp_coeff(46) = 24.67 * trans1
    atm_disp_coeff(47) = 24.67 * trans1
    atm_disp_coeff(48) = 24.67 * trans1
    atm_disp_coeff(49) = 37.32 * trans1
    atm_disp_coeff(50) = 38.71 * trans1
    atm_disp_coeff(51) = 38.44 * trans1
    atm_disp_coeff(52) = 31.74 * trans1
    atm_disp_coeff(53) = 31.50 * trans1
    atm_disp_coeff(54) = 29.99 * trans1
    
    atm_disp_coeff(55) = 315.28* trans1  ! Cs
    atm_disp_coeff(56) = 226.99* trans1
    atm_disp_coeff(57) = 176.25* trans1
    atm_disp_coeff(58) = 140.68* trans1
    atm_disp_coeff(59) = 140.68* trans1
    atm_disp_coeff(60) = 140.68* trans1
    atm_disp_coeff(61) = 140.68* trans1
    atm_disp_coeff(62) = 140.68* trans1
    atm_disp_coeff(63) = 140.68* trans1
    atm_disp_coeff(64) = 140.68* trans1
    atm_disp_coeff(65) = 140.68* trans1
    atm_disp_coeff(66) = 140.68* trans1
    atm_disp_coeff(67) = 140.68* trans1
    atm_disp_coeff(68) = 140.68* trans1
    atm_disp_coeff(69) = 140.68* trans1
    atm_disp_coeff(70) = 140.68* trans1
    atm_disp_coeff(71) = 140.68* trans1
    atm_disp_coeff(72) = 105.11* trans1
    atm_disp_coeff(73) = 81.24 * trans1
    atm_disp_coeff(74) = 81.24 * trans1
    atm_disp_coeff(75) = 81.24 * trans1
    atm_disp_coeff(76) = 81.24 * trans1
    atm_disp_coeff(77) = 81.24 * trans1
    atm_disp_coeff(78) = 81.24 * trans1
    atm_disp_coeff(79) = 81.24 * trans1
    atm_disp_coeff(80) = 57.36 * trans1
    atm_disp_coeff(81) = 57.25 * trans1
    atm_disp_coeff(82) = 63.16 * trans1
    atm_disp_coeff(83) = 63.54 * trans1
    atm_disp_coeff(84) = 55.28 * trans1
    atm_disp_coeff(85) = 57.17 * trans1
    atm_disp_coeff(86) = 56.64 * trans1
    ! atomic vdW radii
    vdW_rad = zero
    vdW_rad(1 ) = 1.001 * trans2   ! H
    vdW_rad(2 ) = 1.012 * trans2
    vdW_rad(3 ) = 0.825 * trans2  ! Li
    vdW_rad(4 ) = 1.408 * trans2
    vdW_rad(5 ) = 1.485 * trans2
    vdW_rad(6 ) = 1.452 * trans2
    vdW_rad(7 ) = 1.397 * trans2
    vdW_rad(8 ) = 1.342 * trans2
    vdW_rad(9 ) = 1.287 * trans2   
    vdW_rad(10) = 1.243 * trans2
    vdW_rad(11) = 1.144 * trans2  ! Na
    vdW_rad(12) = 1.364 * trans2
    vdW_rad(13) = 1.639 * trans2
    vdW_rad(14) = 1.716 * trans2
    vdW_rad(15) = 1.705 * trans2
    vdW_rad(16) = 1.683 * trans2
    vdW_rad(17) = 1.639 * trans2
    vdW_rad(18) = 1.595 * trans2
    vdW_rad(19) = 1.485 * trans2  ! K
    vdW_rad(20) = 1.474 * trans2
    vdW_rad(21) = 1.562 * trans2
    vdW_rad(22) = 1.562 * trans2
    vdW_rad(23) = 1.562 * trans2
    vdW_rad(24) = 1.562 * trans2
    vdW_rad(25) = 1.562 * trans2
    vdW_rad(26) = 1.562 * trans2
    vdW_rad(27) = 1.562 * trans2
    vdW_rad(28) = 1.562 * trans2
    vdW_rad(29) = 1.562 * trans2
    vdW_rad(30) = 1.562 * trans2
    vdW_rad(31) = 1.650 * trans2
    vdW_rad(32) = 1.727 * trans2
    vdW_rad(33) = 1.760 * trans2
    vdW_rad(34) = 1.771 * trans2
    vdW_rad(35) = 1.749 * trans2
    vdW_rad(36) = 1.727 * trans2
    vdW_rad(37) = 1.628 * trans2  ! Rb
    vdW_rad(38) = 1.606 * trans2
    vdW_rad(39) = 1.639 * trans2
    vdW_rad(40) = 1.639 * trans2
    vdW_rad(41) = 1.639 * trans2
    vdW_rad(42) = 1.639 * trans2
    vdW_rad(43) = 1.639 * trans2
    vdW_rad(44) = 1.639 * trans2
    vdW_rad(45) = 1.639 * trans2
    vdW_rad(46) = 1.639 * trans2
    vdW_rad(47) = 1.639 * trans2
    vdW_rad(48) = 1.639 * trans2
    vdW_rad(49) = 1.672 * trans2
    vdW_rad(50) = 1.804 * trans2
    vdW_rad(51) = 1.881 * trans2
    vdW_rad(52) = 1.892 * trans2
    vdW_rad(53) = 1.892 * trans2
    vdW_rad(54) = 1.881 * trans2
    vdW_rad(55) = 1.802 * trans2  ! Cs
    vdW_rad(56) = 1.762 * trans2
    vdW_rad(57) = 1.720 * trans2
    vdW_rad(58) = 1.753 * trans2
    vdW_rad(59) = 1.753 * trans2
    vdW_rad(60) = 1.753 * trans2
    vdW_rad(61) = 1.753 * trans2
    vdW_rad(62) = 1.753 * trans2
    vdW_rad(63) = 1.753 * trans2
    vdW_rad(64) = 1.753 * trans2
    vdW_rad(65) = 1.753 * trans2
    vdW_rad(66) = 1.753 * trans2
    vdW_rad(67) = 1.753 * trans2
    vdW_rad(68) = 1.753 * trans2
    vdW_rad(69) = 1.753 * trans2
    vdW_rad(70) = 1.753 * trans2
    vdW_rad(71) = 1.753 * trans2
    vdW_rad(72) = 1.788 * trans2
    vdW_rad(73) = 1.772 * trans2
    vdW_rad(74) = 1.772 * trans2
    vdW_rad(75) = 1.772 * trans2
    vdW_rad(76) = 1.772 * trans2
    vdW_rad(77) = 1.772 * trans2
    vdW_rad(78) = 1.772 * trans2
    vdW_rad(79) = 1.772 * trans2
    vdW_rad(80) = 1.758 * trans2
    vdW_rad(81) = 1.989 * trans2
    vdW_rad(82) = 1.944 * trans2
    vdW_rad(83) = 1.898 * trans2
    vdW_rad(84) = 2.005 * trans2
    vdW_rad(85) = 1.991 * trans2
    vdW_rad(86) = 1.924 * trans2
    
    allocate ( disp_force(3, ni_in_cell), STAT = stat )
    if (stat .NE. 0) call cq_abort("Error allocating disp_force: ", ni_in_cell)

    return
  end subroutine set_para_D2
!!***

end module DFT_D2
