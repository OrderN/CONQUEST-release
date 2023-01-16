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

    use global_module, only: polS, ne_spin_in_cell, nspin
    use GenComms, only: cq_abort, cq_warn

    implicit none

    ! Passed variables

    ! Local variables
    integer, dimension(:), allocatable :: number_of_bands
    integer :: stat, ispin
    character(len=20) :: subname = "get_polarisation: "

    number_of_bands = ne_spin_in_cell
    do ispin=1,nspin
       ! Get electronic contribution
       ! Allocate polS
       allocate(polS(number_of_bands(ispin), number_of_bands(ispin)), stat=stat)
       if(stat/=0) call cq_abort("Error allocating polS matrix ",number_of_bands(ispin))
       ! Calculate polX matrix
       ! Resta: call to diagonalisation to get polS matrix
       ! Find determinant
       ! Calculate electronic P
       deallocate(polS)
    end do
    ! Get ionic contribution
    ! Output
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
    use GenComms,                    only: gsum,my_barrier
  !     use io_module,                   only: dump_locps
    implicit none

    ! Passed variables
    real(double), intent(out):: Pion_x,Pion_y,Pion_z

    ! Local variables
    integer     :: ipoint, i, ni
    integer     :: atom
    real(double):: q_i, x, y, z
    real(double):: pion_atom_x(ni_in_cell)
    real(double):: pion_atom_y(ni_in_cell)
    real(double):: pion_atom_z(ni_in_cell)
    integer     :: ion_mod(ni_in_cell)
    logical     :: ion_odd
    real(double):: testphase,wrappedtestphase

    pion_atom_x = zero
    pion_atom_y = zero
    pion_atom_z = zero
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
            pion_atom_x(atom) = pion_atom_x(atom) + (x * q_i)
            pion_atom_y(atom) = pion_atom_y(atom) + (y * q_i)
            pion_atom_z(atom) = pion_atom_z(atom) + (z * q_i)
         enddo
      endif
    enddo
    Pion_x = sum(pion_atom_x(1:atom))
    Pion_y = sum(pion_atom_y(1:atom))
    Pion_z = sum(pion_atom_z(1:atom))
    ! We wait for all cores to finish
    call my_barrier()
    ! Sum contributions from each core
    call gsum(Pion_x)
    call gsum(Pion_y)
    call gsum(Pion_z)
    !call wrapphase(Pion_x, one, Pion_x)
    !call wrapphase(Pion_y, one, Pion_y)
    !call wrapphase(Pion_z, one, Pion_z)
    write (io_lun, *) "Pion_x: ",Pion_x
    write (io_lun, *) "Pion_y: ",Pion_y
    write (io_lun, *) "Pion_z: ",Pion_z
    return
  end subroutine get_P_ionic
  !!***

end module polarisation
