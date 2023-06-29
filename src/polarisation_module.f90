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

  real(double), dimension(3) :: Pel_gamma

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
  !!
  !!   Note that we use a slightly confusing approach to the polarisation direction.
  !!   We only need the actual direction, which is found using i_pol_dir(direction),
  !!   to pass to the subroutine get_r_on_atomfns and for output (ionic contribution
  !!   and cell vector).
  !!
  !!   With one direction only i_pol_dir_st = i_pol_dir_end (=1) and the direction is found
  !!   in i_pol_dir(1).  For all three directions we go from 1 to 3 as usual.  The Resta
  !!   S matrix is found for all three directions at once in FindEvals in this case for
  !!   computational efficiency.
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
    use global_module, only: polS, ne_spin_in_cell, atomf, iprint, &
         mat_polX_re, mat_polX_im, mat_polX_re_atomf, mat_polX_im_atomf, io_lun, &
         flag_do_pol_calc, i_pol_dir, i_pol_dir_st, i_pol_dir_end, nspin, spin_factor
    use GenComms, only: cq_abort, cq_warn, inode, ionode, gsum
    use S_matrix_module, only: get_r_on_atomfns
    use functions_on_grid, only: atomfns, allocate_temp_fn_on_grid, &
         free_temp_fn_on_grid, gridfunctions
    use calc_matrix_elements_module, only: get_matrix_elements_new
    use set_bucket_module,           only: rem_bucket, atomf_atomf_rem
    use DiagModule, only: FindEvals
    use matrix_data, only: Srange, aSa_range
    use mult_module, only: S_trans, allocate_temp_matrix, free_temp_matrix, &
         matrix_product_trace, matK, matrix_sum
    use dimens,                      only: r_super_x, r_super_y, r_super_z
    use units, only: eVToJ, BohrToAng

    implicit none

    ! Passed variables

    ! Local variables
    integer, dimension(:), allocatable :: ipiv
    integer :: direction, flag_func, info, matK_rangeS
    integer :: stat, tmp_fn, size, i, spin
    character(len=20) :: subname = "get_polarisation: "
    complex(double_cplx) :: detS
    real(double), dimension(3) :: Pion, cell_vec
    real(double) :: cell_vol, matRe, matIm

    if(inode==ionode) then
       if(iprint>1) then
          write (io_lun, "(2/, 6x, 56('='))")
          write (io_lun, fmt="(2/,20x, 'Bulk Polarization')")
          write (io_lun, fmt="(20x, 18('='),/)")
       else
       end if
    end if
    cell_vec(1) = r_super_x
    cell_vec(2) = r_super_y
    cell_vec(3) = r_super_z
    cell_vol = r_super_x * r_super_y * r_super_z
    ! Electronic contribution
    tmp_fn = allocate_temp_fn_on_grid(atomfns)
    gridfunctions(tmp_fn)%griddata = zero
    do direction = i_pol_dir_st, i_pol_dir_end
       mat_polX_re_atomf(direction) = allocate_temp_matrix(aSa_range, 0, atomf, atomf)
       mat_polX_im_atomf(direction) = allocate_temp_matrix(aSa_range, 0, atomf, atomf)
    end do
    flag_do_pol_calc = .true.
    Pel_gamma = zero
    ! Get ionic contribution
    call get_P_ionic(Pion)
    size = maxval(ne_spin_in_cell)
    ! Allocate ipiv and polS
    allocate(polS(size,size,i_pol_dir_end,nspin), STAT=stat)
    allocate(ipiv(size),STAT=stat)
    polS = zero
    ! ---------------------------------------
    ! Experimental approaches to polarisation
    ! ---------------------------------------
    ! Simple: x = <phi|x|phi>, P = Tr[Kx]/Vol
    if(inode==ionode) write(io_lun,fmt='(/4x,"Experimental approaches to polarisation")')
    if(inode==ionode) write(io_lun,fmt='(/4x,"Use P=Tr[Kx]/Vol with x=<phi|x|phi>")')
    do direction = i_pol_dir_st, i_pol_dir_end
       ! This is not something that should work, particularly as x here is relative to ionic
       ! positions; it's here as a place holder for now
       flag_func = 0 ! Just x/y/z
       call get_r_on_atomfns(i_pol_dir(direction),flag_func,atomfns,tmp_fn)
       call get_matrix_elements_new(inode-1, rem_bucket(atomf_atomf_rem),&
            mat_polX_re_atomf(direction), atomfns, tmp_fn)
       mat_polX_re(direction) = mat_polX_re_atomf(direction)
    end do
    matK_rangeS = allocate_temp_matrix(aSa_range, 0, atomf, atomf)
    do direction = i_pol_dir_st, i_pol_dir_end
       do spin=1,nspin
          ! K to S range
          call matrix_sum(zero, matK_rangeS, one, matK(spin))
          ! Evaluate Tr[KX]
          matRe = matrix_product_trace(matK_rangeS,mat_polX_re(direction))
          if(inode==ionode) then
             write(io_lun,fmt='(4x,"Direction ",i2," Tr[Kx] is ",2e20.12)') &
                  i_pol_dir(direction),matRe
             Pel_gamma(direction) = matRe/cell_vec(i_pol_dir(direction))
             write(io_lun,fmt='(4x,"Total polarisation from Tr[Kx] is   ",e20.12," e / Bohr^2")') &
                  (Pel_gamma(direction) + Pion(i_pol_dir(direction))) &
                  * cell_vec(i_pol_dir(direction))/cell_vol
          end if
       end do
    end do
    call free_temp_matrix(matK_rangeS)
    Pel_gamma = zero
    ! Define X = <phi|exp(i.2pi.x/L)|phi>
    if(inode==ionode) write(io_lun,fmt='(/4x,"Use P=(-1/pi) Im ln S, S=Tr[KX] with X=<phi|exp(i.2pi.x/L|phi>")')
    do direction = i_pol_dir_st, i_pol_dir_end
       ! Calculate polX matrix: <phi_{i\alpha}|exp(i.2pi.r/L|phi_{j\beta}>
       ! This is in *fractional* coordinates so quantum of polarisation is 1
       flag_func = 1 ! real part
       call get_r_on_atomfns(i_pol_dir(direction),flag_func,atomfns,tmp_fn)
       call get_matrix_elements_new(inode-1, rem_bucket(atomf_atomf_rem),&
            mat_polX_re_atomf(direction), atomfns, tmp_fn)
       ! get_r_on_atomfns: imag
       flag_func = 2 ! imag part
       call get_r_on_atomfns(i_pol_dir(direction),flag_func,atomfns,tmp_fn)
       call get_matrix_elements_new(inode-1, rem_bucket(atomf_atomf_rem),&
            mat_polX_im_atomf(direction), atomfns, tmp_fn)
       ! This needs to be a proper MSSF transform
       mat_polX_re(direction) = mat_polX_re_atomf(direction)
       mat_polX_im(direction) = mat_polX_im_atomf(direction)
    end do
    matK_rangeS = allocate_temp_matrix(aSa_range, 0, atomf, atomf)
    do direction = i_pol_dir_st, i_pol_dir_end
       do spin=1,nspin
          call matrix_sum(zero, matK_rangeS, one, matK(spin))
          ! Real and imaginary parts
          matRe = matrix_product_trace(matK_rangeS,mat_polX_re(direction))
          matIm = matrix_product_trace(matK_rangeS,mat_polX_im(direction))
          Pel_gamma(direction) = Pel_gamma(direction) - spin_factor*atan2(matIm,matRe)/twopi
          if(inode==ionode) then
             write(io_lun,fmt='(4x,"Direction ",i2," Tr[KX] is ",2e20.12)') &
                  i_pol_dir(direction),matRe, matIm
             write(io_lun,fmt='(4x,"Direction ",i2," Pel from Tr[KX] is ",e20.12)') &
                  i_pol_dir(direction),Pel_gamma(direction)
             write(io_lun,fmt='(4x,"Total polarisation from Tr[KX] is   ",e20.12," e / Bohr^2")') &
                  (Pel_gamma(direction) + Pion(i_pol_dir(direction))) &
                  * cell_vec(i_pol_dir(direction))/cell_vol
          end if
       end do
    end do
    call free_temp_matrix(matK_rangeS)
    ! Now the Resta approach
    if(inode==ionode) write(io_lun,fmt='(/4x,"Use Resta approach")')
    Pel_gamma = zero
    call FindEvals(ne_spin_in_cell)
    call gsum(polS,size,size,i_pol_dir_end,nspin)
    ! Construct polarisation from determinant of polS matrix
    do direction = i_pol_dir_st, i_pol_dir_end
       do spin=1,nspin
          ! Find determinant: call zgetrf to decompose polS into P.U.L
          call zgetrf(size,size,polS(:,:,direction,spin),size,ipiv,info)
          if(info>0) write(io_lun,fmt='(4x,"Problem with zgetrf; info=",i4)') spin,info
          ! Take product of diagonals of polS (on output) to get determinant
          detS = cmplx(one,zero,double_cplx)
          do i=1,size
             detS = detS * polS(i,i,direction,spin)
             ! Permutation: if ipiv(i)/=i, scale by -1
             if(ipiv(i)/=i) detS = -detS
          end do
          if(inode==ionode .and. iprint>2) &
               write(io_lun,fmt='(/4x,"detS real ",e20.12," detS imag ",e20.12/)') &
               real(detS), aimag(detS)
          ! Calculate electronic P: Im ln polS is just finding phase NB 1/pi includes factor of 2 for spin
          Pel_gamma(direction) = Pel_gamma(direction) - spin_factor*atan2(aimag(detS),real(detS))/twopi
          if(inode==ionode .and. iprint>2) &
               write(io_lun,fmt='(4x,"Direction ",i2," Pel is ",e20.12)') &
               i_pol_dir(direction),Pel_gamma(direction)
       end do
    end do
    deallocate(ipiv)
    deallocate(polS)
    ! Output - include quantum of polarisation
    ! The quantum is \frac{e}{V_{cell}} \mathbf{R} for lattice vector R
    if(inode==ionode) then
       if(iprint>1) then
          ! NB we always calculate all three directions for ionic, but potentially limit electronic
          ! so the indexing is different - this is correct
          do direction=i_pol_dir_st, i_pol_dir_end
             if(inode==ionode .and. iprint>2) &
                  write(io_lun,fmt='(4x,"Direction ",i2," Pion is ",e20.12/)') &
                  i_pol_dir(direction),Pion(i_pol_dir(direction))
             write(io_lun,fmt='(4x,"Direction: ",i2)') i_pol_dir(direction)
             write(io_lun,fmt='(4x,"Total polarisation:      ",e20.12," e / Bohr^2")') &
                  (Pel_gamma(direction) + Pion(i_pol_dir(direction))) &
                  * cell_vec(i_pol_dir(direction))/cell_vol
             write(io_lun,fmt='(4x,"Quantum of polarisation: ",e20.12," e / Bohr^2")') &
                  cell_vec(i_pol_dir(direction))/cell_vol
             write(io_lun,fmt='(4x,"Total polarisation:      ",e20.12," C / m^2")') &
                  (Pel_gamma(direction) + Pion(i_pol_dir(direction))) &
                  * cell_vec(i_pol_dir(direction))*eVToJ/(cell_vol*BohrToAng*BohrToAng*1e-20_double)
          end do
       else
          do direction=i_pol_dir_st, i_pol_dir_end
             write(io_lun,fmt='(4x,"Direction: ",i2)') i_pol_dir(direction)
             write(io_lun,fmt='(4x,"Total polarisation:      ",e20.12," e / Bohr^2")') &
                  (Pel_gamma(direction) + Pion(i_pol_dir(direction))) &
                  * cell_vec(i_pol_dir(direction))/cell_vol
             write(io_lun,fmt='(4x,"Quantum of polarisation: ",e20.12," e / Bohr^2")') &
                  cell_vec(i_pol_dir(direction))/cell_vol
          end do
       end if
    end if
    ! Free space
    do direction = i_pol_dir_end, i_pol_dir_st, -1
       call free_temp_matrix(mat_polX_im_atomf(direction))
       call free_temp_matrix(mat_polX_re_atomf(direction))
    end do
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
  subroutine get_P_ionic(Pion)

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
    real(double), dimension(3) :: Pion

    ! Local variables
    integer     :: ipoint, i, ni
    integer     :: atom
    real(double):: q_i, x, y, z
    real(double):: testphase,wrappedtestphase

    Pion = zero
    atom = 0
    do ipoint = 1, bundle%groups_on_node
      if(bundle%nm_nodgroup(ipoint) > 0) then
         do ni = 1, bundle%nm_nodgroup(ipoint)
            atom = atom + 1
            i = bundle%ig_prim(bundle%nm_nodbeg(ipoint)+ni-1)
!           --- We want the positions of atom i in fractional coordinates ---
!           --- (hence we divide by r_super_x)                            ---
            x = atom_coord(1,i) / r_super_x
            y = atom_coord(2,i) / r_super_y
            z = atom_coord(3,i) / r_super_z
            q_i = charge(bundle%species(bundle%nm_nodbeg(ipoint)+ni-1))

            ! write (io_lun, *) "q_i: ",q_i
            Pion(1) = Pion(1) + (x * q_i)
            Pion(2) = Pion(2) + (y * q_i)
            Pion(3) = Pion(3) + (z * q_i)
         enddo
      endif
    enddo
    ! We wait for all cores to finish
    call my_barrier()
    ! Sum contributions from each core
    call gsum(Pion,3)
    ! We have found P in fractional coordinates, so quantum is 1
    Pion(1) = Pion(1) - floor(Pion(1))
    Pion(2) = Pion(2) - floor(Pion(2))
    Pion(3) = Pion(3) - floor(Pion(3))
    return
  end subroutine get_P_ionic
  !!***

end module polarisation
