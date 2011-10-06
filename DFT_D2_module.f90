! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!
!! SOURCE
module DFT_D2

  use datatypes

  implicit none
  
  ! RCS tag for object file identification
  character(len=80), private :: RCSid = "$Id: DFT_D2_module.f90 85 2011-06-09 14:46:30Z M. Arita $"

  real(double) :: disp_energy
  real(double), allocatable :: disp_force(:, :)
  real(double), allocatable :: atm_disp_coeff(:), vdW_rad(:)

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
!!  SOURCE
!!
  subroutine dispersion_D2

    use cover_module, ONLY: D2_CS, make_CS
    use datatypes
    use units
    use dimens
    use GenComms, ONLY: cq_abort, gcopy, inode, ionode, gsum, my_barrier
    use group_module, ONLY: parts
    use global_module, ONLY: io_lun, numprocs, flag_functional_type, ni_in_cell, species_glob, &
                             id_glob, flag_only_dispersion
    use numbers
    use primary_module, ONLY: bundle
    use io_module, ONLY: get_file_name
    use species_module, ONLY: n_species, species_label
    use input_module, ONLY: io_assign, io_close
    implicit none

    ! local variables
    integer :: i, j, ip, ia, jp, ja, ia_glob
    integer :: stat
    integer :: nppx, nppy, nppz, ind_part, nccx, nccy, nccz, m1, m2, m3, m, ig_atom_beg
    integer, parameter :: d = 20.0_double
    integer, allocatable :: neighbour_part(:)
    real(double) :: s_6, real_cell_vec(3, 3), part_cell_vec(3, 3), v1, v2, v3, rr(3), part_cell_dual(3, 3), &
                    distance, x_ia, y_ia, z_ia, x_ja, y_ja, z_ja, rijx, rijy, rijz, rij, rij2, rij6, power, &
                    f_damp, denominator, emp_energy, emp_force, unit_vec_x, unit_vec_y, unit_vec_z
    real(double) :: C_i, C_j, C_ij, r_vdW_i, r_vdW_j, R_r

    ! Determine the s_6 value
    if (flag_functional_type .EQ. 102) then
      ! for revPBE
      s_6 = 1.25_double
    elseif (flag_functional_type .EQ. 101) then
      ! for PBE
      s_6 = 0.75_double
    else
      call cq_abort("Error determining s_6: ")
    endif

    ! Initialisation of the dispersion energy and the atomic forces
    disp_energy = 0.0_double
    disp_force  = 0.0_double

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
                if ( (rij2 .GT. very_small) .AND. (rij2 .LE. r_dft_d2*r_dft_d2) ) then
                  ig_atom_beg = parts%icell_beg(D2_CS%lab_cell(jp))
                  ! Calculate the molecular dispersion coefficient and vdW radius
                  C_i = atm_disp_coeff(bundle%species(bundle%nm_nodbeg(ip) + ia - 1))
                  C_j = atm_disp_coeff(species_glob(id_glob(ig_atom_beg+ja-1)))
                  C_ij = sqrt(C_i * C_j)
                  r_vdW_i = vdW_rad(bundle%species(bundle%nm_nodbeg(ip) + ia - 1))
                  r_vdW_j = vdW_rad(species_glob(id_glob(ig_atom_beg+ja-1)))
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

    ! Works only when you want to get the dispersion
    if (flag_only_dispersion) then
      if (inode == ionode) then
        write (io_lun, '(/a)') "***** Dispersion Information: DFT-D2 *****"
        write (io_lun,16) en_conv*disp_energy, en_units(energy_units)
        write (io_lun,fmt='(/,20x,"Disp. Force on atoms (",a2"/",a2,")"/)') en_units(energy_units), d_units(dist_units)
        write (io_lun,fmt='("Atom        X               Y               Z")')
        do i = 1, ni_in_cell
            write (io_lun, 109) i, (for_conv*disp_force(j,i),j=1,3)
        enddo
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
! Subroutine read_para_D2
! --------------------------------------------------------------

!!****f* DFT_D2/read_para_D2 *
!!
!!  NAME
!!    read_para_D2
!!  USAGE
!!
!!  PURPOSE
!!    Reads and stores the necessary parameters in the DFT-D2 scheme
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
!!  SOURCE
!!
  subroutine read_para_D2

    use datatypes
    use GenComms, ONLY: cq_abort, gcopy, inode, ionode, my_barrier
    use global_module, ONLY: io_lun, numprocs, ni_in_cell
    use io_module, ONLY: get_file_name
    use species_module, ONLY: n_species, species_label
    use input_module, ONLY: io_assign, io_close
    implicit none

    ! local variables
    integer :: i, j, stat, ios, lun
    real(double) :: coeff, radius
    real(double), parameter :: trans1 = 17.3452771_double ! [Jnm^6/mol]-->[Ha bohr^6]
    real(double), parameter :: trans2 = 1.889726125_double ! [ang]-->[bohr]
    character(len=20) dummy
    character(len=2)  kind

    ! debug variables
    integer :: lun_para
    character(len=30) :: filename

    allocate( atm_disp_coeff(n_species), vdW_rad(n_species), STAT=stat )
    if (stat .NE. 0) call cq_abort("Error allocating DFT-D2 parameters: ",n_species)
    ! Fetch and store the parameters
    if (inode==ionode) then
       call io_assign(lun)
       open (unit=lun, file='para_D2.dat', status='old', action = 'read', &
            iostat=ios)
       if (ios>0) call cq_abort("Reading para_D2.dat: file error")
       do i = 1, n_species
          read  (lun, *) dummy
          do while ( trim(dummy)/="</parameters>" )
             read (lun, '(a)')  dummy
          enddo !while
          do j = 1, 55
             read  (lun, *) kind, coeff, radius
             if ( trim(kind) .EQ. species_label(i) ) then
                atm_disp_coeff(i) = coeff * trans1 ![Ha bohr^6]
                vdW_rad(i) = radius * trans2          ![bohr]
                exit
             endif
          enddo !(do, j)
          rewind (lun)
       enddo !(i, n_species)
       call io_close(lun)
    endif

    call gcopy(atm_disp_coeff, n_species)
    call gcopy(vdW_rad, n_species)

    allocate ( disp_force(3, ni_in_cell), STAT = stat )
    if (stat .NE. 0) call cq_abort("Error allocating disp_force: ", ni_in_cell)

    return
  end subroutine read_para_D2
!!***

  subroutine partition_distance(pd_init,aa,beta,rr,distance)

    use numbers
    use global_module, ONLY: io_lun
    use GenComms, ONLY: cq_abort, myid

    implicit none

    logical, intent(in) :: pd_init
    real(double), intent(in), optional, dimension(3) :: rr
    real(double), intent(in), dimension(3,3) :: aa
    real(double), intent(inout), dimension(3,3) :: beta
    real(double), intent(out), optional :: distance

    integer :: i, n
    integer :: isig1, isig2, isig3
    real(double) :: alpha1, alpha2, alpha3
    real(double) :: asig1, asig2, asig3, a1_dot_a1, a2_dot_a2, a3_dot_a3, &
         &a2_dot_a3, a3_dot_a1, a1_dot_a2, b1_dot_b1, b2_dot_b2, b3_dot_b3, &
         &b2_dot_b3, b3_dot_b1, b1_dot_b2, dx, dy, dz, p1, p2, p3, x1, x2, x3, &
         &y1, y2, y3, z1, z2, z3
    real(double), dimension(3) :: edg_a, edg_b, q_vec


    ! --- initialisation for given cell vectors: do only if pd_int = .true.
    if(pd_init) then
       ! --- construct required quantities from cell vectors -----

       a1_dot_a1 = aa(1,1)*aa(1,1) + aa(1,2) * aa(1,2) + aa(1,3)*aa(1,3)
       a2_dot_a2 = aa(2,1)*aa(2,1) + aa(2,2) * aa(2,2) + aa(2,3)*aa(2,3)
       a3_dot_a3 = aa(3,1)*aa(3,1) + aa(3,2) * aa(3,2) + aa(3,3)*aa(3,3)
       a2_dot_a3 = aa(2,1)*aa(3,1) + aa(2,2) * aa(3,2) + aa(2,3)*aa(3,3)
       a3_dot_a1 = aa(3,1)*aa(1,1) + aa(3,2) * aa(1,2) + aa(3,3)*aa(1,3)
       a1_dot_a2 = aa(1,1)*aa(2,1) + aa(1,2) * aa(2,2) + aa(1,3)*aa(2,3)

       beta(1,1) = aa(2,2)*aa(3,3) - aa(2,3)*aa(3,2)
       beta(1,2) = aa(2,3)*aa(3,1) - aa(2,1)*aa(3,3)
       beta(1,3) = aa(2,1)*aa(3,2) - aa(2,2)*aa(3,1)
       beta(2,1) = aa(3,2)*aa(1,3) - aa(3,3)*aa(1,2)
       beta(2,2) = aa(3,3)*aa(1,1) - aa(3,1)*aa(1,3)
       beta(2,3) = aa(3,1)*aa(1,2) - aa(3,2)*aa(1,1)
       beta(3,1) = aa(1,2)*aa(2,3) - aa(1,3)*aa(2,2)
       beta(3,2) = aa(1,3)*aa(2,1) - aa(1,1)*aa(2,3)
       beta(3,3) = aa(1,1)*aa(2,2) - aa(1,2)*aa(2,1)

       alpha1 = (aa(1,1)*beta(1,1) + aa(1,2)*beta(1,2) + aa(1,3)*beta(1,3)) / &
            &(beta(1,1)*beta(1,1) + beta(1,2)*beta(1,2) + beta(1,3)*beta(1,3))
       alpha2 = (aa(2,1)*beta(2,1) + aa(2,2)*beta(2,2) + aa(2,3)*beta(2,3)) / &
            &(beta(2,1)*beta(2,1) + beta(2,2)*beta(2,2) + beta(2,3)*beta(2,3))
       alpha3 = (aa(3,1)*beta(3,1) + aa(3,2)*beta(3,2) + aa(3,3)*beta(3,3)) / &
            &(beta(3,1)*beta(3,1) + beta(3,2)*beta(3,2) + beta(3,3)*beta(3,3))

       beta(1,1) = alpha1*beta(1,1)
       beta(1,2) = alpha1*beta(1,2)
       beta(1,3) = alpha1*beta(1,3)
       beta(2,1) = alpha2*beta(2,1)
       beta(2,2) = alpha2*beta(2,2)
       beta(2,3) = alpha2*beta(2,3)
       beta(3,1) = alpha3*beta(3,1)
       beta(3,2) = alpha3*beta(3,2)
       beta(3,3) = alpha3*beta(3,3)

       b1_dot_b1 = beta(1,1)*beta(1,1) + beta(1,2)*beta(1,2) + beta(1,3)*beta(1,3)
       b2_dot_b2 = beta(2,1)*beta(2,1) + beta(2,2)*beta(2,2) + beta(2,3)*beta(2,3)
       b3_dot_b3 = beta(3,1)*beta(3,1) + beta(3,2)*beta(3,2) + beta(3,3)*beta(3,3)
       b2_dot_b3 = beta(2,1)*beta(3,1) + beta(2,2)*beta(3,2) + beta(2,3)*beta(3,3)
       b3_dot_b1 = beta(3,1)*beta(1,1) + beta(3,2)*beta(1,2) + beta(3,3)*beta(1,3)
       b1_dot_b2 = beta(1,1)*beta(2,1) + beta(1,2)*beta(2,2) + beta(1,3)*beta(2,3)
    else
       a1_dot_a1 = aa(1,1)*aa(1,1) + aa(1,2) * aa(1,2) + aa(1,3)*aa(1,3)
       a2_dot_a2 = aa(2,1)*aa(2,1) + aa(2,2) * aa(2,2) + aa(2,3)*aa(2,3)
       a3_dot_a3 = aa(3,1)*aa(3,1) + aa(3,2) * aa(3,2) + aa(3,3)*aa(3,3)
       a2_dot_a3 = aa(2,1)*aa(3,1) + aa(2,2) * aa(3,2) + aa(2,3)*aa(3,3)
       a3_dot_a1 = aa(3,1)*aa(1,1) + aa(3,2) * aa(1,2) + aa(3,3)*aa(1,3)
       a1_dot_a2 = aa(1,1)*aa(2,1) + aa(1,2) * aa(2,2) + aa(1,3)*aa(2,3)
       b1_dot_b1 = beta(1,1)*beta(1,1) + beta(1,2)*beta(1,2) + beta(1,3)*beta(1,3)
       b2_dot_b2 = beta(2,1)*beta(2,1) + beta(2,2)*beta(2,2) + beta(2,3)*beta(2,3)
       b3_dot_b3 = beta(3,1)*beta(3,1) + beta(3,2)*beta(3,2) + beta(3,3)*beta(3,3)
       b2_dot_b3 = beta(2,1)*beta(3,1) + beta(2,2)*beta(3,2) + beta(2,3)*beta(3,3)
       b3_dot_b1 = beta(3,1)*beta(1,1) + beta(3,2)*beta(1,2) + beta(3,3)*beta(1,3)
       b1_dot_b2 = beta(1,1)*beta(2,1) + beta(1,2)*beta(2,2) + beta(1,3)*beta(2,3)
    endif
    if(pd_init) return

    ! --- components of supplied position ----------------------
    x1 = (beta(1,1)*rr(1) + beta(1,2)*rr(2) + beta(1,3)*rr(3))/b1_dot_b1
    x2 = (beta(2,1)*rr(1) + beta(2,2)*rr(2) + beta(2,3)*rr(3))/b2_dot_b2
    x3 = (beta(3,1)*rr(1) + beta(3,2)*rr(2) + beta(3,3)*rr(3))/b3_dot_b3

    if(abs(x1)>1e5.AND.myid==0) then
       write(io_lun,*) 'Error: betas, rs: ',beta(1,:),rr(1)
       write(io_lun,*) 'Error: betas, rs: ',beta(2,:),rr(2)
       write(io_lun,*) 'Error: betas, rs: ',beta(3,:),rr(3)
    end if
    y1 = (aa(1,1)*rr(1) + aa(1,2)*rr(2) + aa(1,3)*rr(3))/b1_dot_b1
    y2 = (aa(2,1)*rr(1) + aa(2,2)*rr(2) + aa(2,3)*rr(3))/b2_dot_b2
    y3 = (aa(3,1)*rr(1) + aa(3,2)*rr(2) + aa(3,3)*rr(3))/b3_dot_b3

    z2 = two*x2 - (two*x3-sign(min(abs(two*x3),one),x3))*b2_dot_b3/b2_dot_b2
    z3 = two*x3 - (two*x2-sign(min(abs(two*x2),one),x2))*b2_dot_b3/b3_dot_b3
    p1 = (two/a1_dot_a1) * (y1*b1_dot_b1 - &
         &half*a1_dot_a2*sign(min(abs(z2),one),z2) -&
         &half*a3_dot_a1*sign(min(abs(z3),one),z3))
    z3 = two*x3 - (two*x1-sign(min(abs(two*x1),one),x1))*b3_dot_b1/b3_dot_b3
    z1 = two*x1 - (two*x3-sign(min(abs(two*x3),one),x3))*b3_dot_b1/b1_dot_b1
    p2 = (two/a2_dot_a2) * (y2*b2_dot_b2 - &
         &half*a2_dot_a3*sign(min(abs(z3),one),z3) -&
         &half*a1_dot_a2*sign(min(abs(z1),one),z1))
    z1 = two*x1 - (two*x2-sign(min(abs(two*x2),one),x2))*b1_dot_b2/b1_dot_b1
    z2 = two*x2 - (two*x1-sign(min(abs(two*x1),one),x1))*b1_dot_b2/b2_dot_b2
    p3 = (two/a3_dot_a3) * (y3*b3_dot_b3 - &
         &half*a3_dot_a1*sign(min(abs(z1),one),z1) -&
         &half*a2_dot_a3*sign(min(abs(z2),one),z2))

    isig1 = sign(min(1,floor(abs(p1))),nint(p1))
    isig2 = sign(min(1,floor(abs(p2))),nint(p2))
    isig3 = sign(min(1,floor(abs(p3))),nint(p3))

      !write(io_lun,fmt='(/" info for supplied position:"/)')
      !write(io_lun,fmt='(" x1, x2, x3:",3f12.6)') x1, x2, x3
      !write(io_lun,fmt='(" y1, y2, y3:",3f12.6)') y1, y2, y3
      !write(io_lun,fmt='(" p1, p2, p3:",3f12.6)') p1, p2, p3
      !write(io_lun,fmt='(/" signatures of position:",3i3)') isig1, isig2, isig3

    ! --- calculate shortest distance --------------------------
    asig1 = real(isig1,double)
    asig2 = real(isig2,double)
    asig3 = real(isig3,double)

    select case(abs(isig1)+abs(isig2)+abs(isig3))
       ! --- vertex: no zero signatures
    case(3)
       dx = rr(1) - half*(asig1*aa(1,1) + asig2*aa(2,1) + asig3*aa(3,1))
       dy = rr(2) - half*(asig1*aa(1,2) + asig2*aa(2,2) + asig3*aa(3,2))
       dz = rr(3) - half*(asig1*aa(1,3) + asig2*aa(2,3) + asig3*aa(3,3))
       distance = sqrt(dx*dx+dy*dy+dz*dz)
       ! --- edge: only one zero signature
    case(2)
       edg_a(1) = half*(asig1*aa(1,1) + asig2*aa(2,1) + asig3*aa(3,1)) - rr(1)
       edg_a(2) = half*(asig1*aa(1,2) + asig2*aa(2,2) + asig3*aa(3,2)) - rr(2)
       edg_a(3) = half*(asig1*aa(1,3) + asig2*aa(2,3) + asig3*aa(3,3)) - rr(3)
       edg_b(1) = (one-abs(asig1))*aa(1,1) + (one-abs(asig2))*aa(2,1) + &
            &(one-abs(asig3))*aa(3,1)
       edg_b(2) = (one-abs(asig1))*aa(1,2) + (one-abs(asig2))*aa(2,2) + &
            &(one-abs(asig3))*aa(3,2)
       edg_b(3) = (one-abs(asig1))*aa(1,3) + (one-abs(asig2))*aa(2,3) + &
            &(one-abs(asig3))*aa(3,3)
       !ORIdistance = sqrt(edg_a(1)*edg_a(1) + edg_a(2)*edg_a(2) + edg_a(3)*edg_a(3) - &
       !ORI     &((edg_a(1)*edg_b(1) + edg_a(2)*edg_b(2) + edg_a(3)*edg_b(3))**2) / &
       !ORI     &(edg_b(1)*edg_b(1) + edg_b(2)*edg_b(2) + edg_b(3)*edg_b(3)))
       distance = edg_a(1)*edg_a(1) + edg_a(2)*edg_a(2) + edg_a(3)*edg_a(3) - &
            &((edg_a(1)*edg_b(1) + edg_a(2)*edg_b(2) + edg_a(3)*edg_b(3))**2) / &
            &(edg_b(1)*edg_b(1) + edg_b(2)*edg_b(2) + edg_b(3)*edg_b(3))
       if(distance > very_small) then 
           distance = sqrt(distance)
       else
           distance = zero
       endif

       ! --- face: only one non-zero signature
    case(1)
       q_vec(1) = half*(asig1*beta(1,1) + asig2*beta(2,1) + asig3*beta(3,1))
       q_vec(2) = half*(asig1*beta(1,2) + asig2*beta(2,2) + asig3*beta(3,2))
       q_vec(3) = half*(asig1*beta(1,3) + asig2*beta(2,3) + asig3*beta(3,3))
       distance = abs(q_vec(1)*q_vec(1) + q_vec(2)*q_vec(2) + q_vec(3)*q_vec(3) - &
            &q_vec(1)*rr(1) - q_vec(2)*rr(2) - q_vec(3)*rr(3)) / &
            &sqrt(q_vec(1)*q_vec(1) + q_vec(2)*q_vec(2) + q_vec(3)*q_vec(3))
       ! --- body: all zero signatures
    case(0)
       distance = zero
       ! --- case default included as a check on coding
    case default
       write(io_lun,*) 'sigs: ',isig1,isig2,isig3
       call cq_abort('partition_distance: default found ',abs(isig1)+abs(isig2)+abs(isig3))
    end select

  end subroutine partition_distance




end module DFT_D2
