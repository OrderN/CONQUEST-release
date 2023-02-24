! $Id$
! -----------------------------------------------------------
! Module polarisation
! -----------------------------------------------------------
! Code area 3: Operators
! -----------------------------------------------------------

!!****h* Conquest/polarisation
!!  NAME
!!   polarisation
!!  PURPOSE
!!   Calculates polarisation (using Resta approach among others)
!!  AUTHOR
!!   D.R.Bowler (based on work by K. Shenton)
!!  CREATION DATE
!!   2023/01/12
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
module polarisation

  use datatypes

  implicit none

  real(double) :: Pel_gamma

!!***

contains

  ! -----------------------------------------------------------
  ! Subroutine 
  ! -----------------------------------------------------------

  !!****f* polarisation/get_polarisation *
  !!
  !!  NAME 
  !!   get_polarisation
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Calculates polarisation
  !!  INPUTS
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   2023/01/12
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  subroutine get_polarisation()

    use numbers
    use global_module, only: polS, ne_spin_in_cell, nspin, atomf, &
         mat_polX_re, mat_polX_im, mat_polX_re_atomf, mat_polX_im_atomf, io_lun, flag_do_pol_calc
    use GenComms, only: cq_abort, cq_warn, inode, ionode
    use S_matrix_module, only: get_r_on_atomfns
    use functions_on_grid, only: atomfns, allocate_temp_fn_on_grid, &
         free_temp_fn_on_grid, gridfunctions
    use calc_matrix_elements_module, only: get_matrix_elements_new
    use set_bucket_module,           only: rem_bucket, atomf_atomf_rem
    use DiagModule, only: FindEvals
    use matrix_data, only: Srange, aSa_range
    use mult_module, only: S_trans, allocate_temp_matrix

    implicit none

    ! Passed variables

    ! Local variables
    integer, dimension(:), allocatable :: ipiv
    integer, dimension(2) :: number_of_bands
    integer :: direction, flag_func, info
    integer :: stat, ispin, tmp_fn, size, i
    character(len=20) :: subname = "get_polarisation: "
    complex(double_cplx) :: detS
    real(double) :: Pion_x, Pion_y, Pion_z

    number_of_bands = int(ne_spin_in_cell)
    direction = 1
    tmp_fn = allocate_temp_fn_on_grid(atomfns)
    gridfunctions(tmp_fn)%griddata = zero
    mat_polX_re_atomf = allocate_temp_matrix(aSa_range, 0, atomf, atomf)
    mat_polX_im_atomf = allocate_temp_matrix(aSa_range, 0, atomf, atomf)    
    flag_do_pol_calc = .true.
    do ispin=1,nspin
       size = number_of_bands(ispin)
       ! Get electronic contribution
       ! Allocate ipiv and polS
       allocate(polS(size,size), STAT=stat)
       allocate(ipiv(size),STAT=stat)
       ! Calculate polX matrix
       ! get_r_on_atomfns: real
       flag_func = 1
       call get_r_on_atomfns(direction,flag_func,atomfns,tmp_fn)
       call get_matrix_elements_new(inode-1, rem_bucket(atomf_atomf_rem),&
                                       mat_polX_re_atomf, atomfns, tmp_fn)
       ! get_r_on_atomfns: imag
       flag_func = 2
       call get_r_on_atomfns(direction,flag_func,atomfns,tmp_fn)
       call get_matrix_elements_new(inode-1, rem_bucket(atomf_atomf_rem),&
            mat_polX_im_atomf, atomfns, tmp_fn)
       mat_polX_re = mat_polX_re_atomf
       mat_polX_im = mat_polX_im_atomf
       ! Resta: call to diagonalisation to get polS matrix
       call FindEvals(ne_spin_in_cell)
       ! Find determinant: call dgetrf to decompose polS into P.U.L
       call zgetrf(size,size,polS,size,ipiv,info)
       ! Take product of diagonals of polS (on output) to get determinant
       detS = cmplx(one,zero,double_cplx)
       do i=1,size
           detS = detS * polS(i,i)
           ! Permutation: if ipiv(i)/=i, scale by -1
           if(ipiv(i)/=i) detS = -detS
        end do
        write(io_lun,fmt='(4x,"detS real ",e20.12," detS imag ",e20.12)') real(detS), aimag(detS)
       ! Calculate electronic P: Im ln polS is just finding phase
       Pel_gamma = -atan2(aimag(detS),real(detS))/pi
       write(io_lun,fmt='(4x,"Pel is ",e20.12)') Pel_gamma
       deallocate(ipiv)
       deallocate(polS)
    end do
    ! Get ionic contribution
    call get_P_ionic(Pion_x, Pion_y, Pion_z)
    ! Output
    ! Do we also test Tr[K.polX] for different forms of polX? (x, Resta and other)
    ! Free space
    call free_temp_fn_on_grid(tmp_fn)
    flag_do_pol_calc = .false.
    return
  end subroutine get_polarisation
  !!***

  ! -----------------------------------------------------------
  ! Subroutine get_P_ionic
  ! -----------------------------------------------------------
  
  !!****f* polarisation_module/get_P_ionic *
  !!
  !!  NAME
  !!   get_P_ionic
  !!  USAGE
  !!
  !!  PURPOSE
  !!   P_ionic is given by:
  !!   $\mathbf{P}_{ion} = \frac{e}{V_{cell}} \sum_I {Z^{ion}_I \mathbf{R}_I}$
  !!
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   J. Kane Shenton and Dave Bowler
  !!  CREATION DATE
  !!   22/04/2015
  !!  MODIFICATION HISTORY
  !!   2023/01/13 08:31 dave
  !!    Incorporated into new repository
  !!  SOURCE
  !!
  subroutine get_P_ionic(Pion_x,Pion_y,Pion_z)

    use datatypes
    use numbers
    use global_module,               only: atom_coord, ni_in_cell, &
                                           nspin, area_ops, &
                                           iprint_gen, iprint_ops, io_lun
    use primary_module,              only: bundle, domain
    use species_module,              only: charge, species
    use dimens,                      only: r_super_x, r_super_y, r_super_z
    use GenComms,                    only: gsum,my_barrier, inode, ionode
  !     use io_module,                   only: dump_locps
    implicit none

    ! Passed variables
    real(double), intent(out):: Pion_x,Pion_y,Pion_z

    ! Local variables
    integer     :: ipoint, i, ni
    integer     :: atom
    real(double):: q_i, x, y, z
    real(double):: testphase,wrappedtestphase

    Pion_x = zero
    Pion_y = zero
    Pion_z = zero
    atom = 0
    do ipoint = 1, bundle%groups_on_node
      if(bundle%nm_nodgroup(ipoint) > 0) then
         do ni = 1, bundle%nm_nodgroup(ipoint)
            atom = atom + 1
            i = bundle%ig_prim(bundle%nm_nodbeg(ipoint)+ni-1)
!           --- We want the positions of atom i in fractional coordinates ---
!           --- (hence we divide by r_super_x)                            ---
            x = atom_coord(1,i) !/ r_super_x
            y = atom_coord(2,i) !/ r_super_y
            z = atom_coord(3,i) !/ r_super_z
            q_i = charge(bundle%species(bundle%nm_nodbeg(ipoint)+ni-1))

            ! write (io_lun, *) "q_i: ",q_i
            Pion_x = Pion_x + (x * q_i)
            Pion_y = Pion_y + (y * q_i)
            Pion_z = Pion_z + (z * q_i)
         enddo
      endif
    enddo
    ! We wait for all cores to finish
    call my_barrier()
    ! Sum contributions from each core
    call gsum(Pion_x)
    call gsum(Pion_y)
    call gsum(Pion_z)
    !call wrapphase(Pion_x, one, Pion_x)
    !call wrapphase(Pion_y, one, Pion_y)
    !call wrapphase(Pion_z, one, Pion_z)
    if(inode==ionode) then
       write (io_lun, *) "Pion_x: ",Pion_x
       write (io_lun, *) "Pion_y: ",Pion_y
       write (io_lun, *) "Pion_z: ",Pion_z
    end if
    return
  end subroutine get_P_ionic
  !!***

end module polarisation
