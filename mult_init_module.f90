! $Id: mult_init_module.f90,v 1.2 2002/03/12 14:56:38 drb Exp $
! -----------------------------------------------------------
! Module mult_init_module
! -----------------------------------------------------------
! Code area 2: matrices
! -----------------------------------------------------------

!!****h* Conquest/mult_init_module *
!!  NAME
!!   mult_init_module
!!  PURPOSE
!!   Routines to initialise matrix multiplications.
!!   Contains mult_ini and check_mm
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!   27/01/00 by D.R.Bowler to use new communications
!!   19/02/00 by DRB to implement transpose initialisation
!!   02/03/00 by DRB to finish transpose work
!!   12/04/00 by DRB to add itran_addr to trans_remote
!!   09/06/00 by DRB to adapt for Conquest
!!   20/06/2001 dave
!!    Added ROBODoc headers, RCS Id and Log tags and cq_abort
!!***
module mult_init_module

contains

!!****f* mult_init_module/mult_ini *
!!
!!  NAME 
!!   mult_ini
!!  USAGE
!! 
!!  PURPOSE
!!   Initiates a multiplication given by the
!!   matrices in a_b_c.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!   20/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!  SOURCE
!!
  subroutine mult_ini(a_b_c,ind,myid,n_prim,parts)

    ! Module usage
    use GenComms, ONLY: cq_abort
    use matrix_module
    use basic_types
    use matrix_comms_module

    ! Passed variables
    integer :: n_prim,myid,mx_part
    type(matrix_mult)::a_b_c
    type(group_set) :: parts
    integer, dimension(:), target :: ind

    ! Local variables
    integer :: ierr, indlen

    indlen = size(ind)
    a_b_c%bindex => ind(1:indlen)
    ! Check for mult_type
    if(a_b_c%mult_type/=1.AND.a_b_c%mult_type/=2) then
       call cq_abort('mult_ini: mult_type must be 1 or 2: ',a_b_c%mult_type)
    endif
    ! Check that the dimension are all correct
    if(n_prim.gt.0)then
       call check_mm( a_b_c%ahalo%np_in_halo,&
            a_b_c%ahalo%lab_hcell,a_b_c%ahalo%j_beg,a_b_c%ahalo%nh_part,&
            a_b_c%ahalo%j_seq,a_b_c%ahalo%i_halo,&
            a_b_c%ahalo%mx_part, &
            parts%mx_gcell,a_b_c%ahalo%mx_halo,parts%mx_mem_grp)
       call check_mm( a_b_c%chalo%np_in_halo,&
            a_b_c%chalo%lab_hcell,a_b_c%chalo%j_beg,a_b_c%chalo%nh_part,&
            a_b_c%chalo%j_seq,a_b_c%chalo%i_halo,&
            a_b_c%chalo%mx_part, &
            parts%mx_gcell,a_b_c%chalo%mx_halo,parts%mx_mem_grp)
    endif
    call init_mult_comms(parts,a_b_c,myid)
    return
  end subroutine mult_ini
!!***

!!****f* mult_init_module/check_mm *
!!
!!  NAME 
!!   check_mm
!!  USAGE
!! 
!!  PURPOSE
!!   Checks limits for a matrix multiplication
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!   20/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!  SOURCE
!!
  subroutine check_mm(np_in_ahalo, lab_ahcell,kbeg,nah_part,kseq,k_halo,&
       mx_apart,mx_pcell, mx_ahalo,mx_part)

    ! Module usage
    use datatypes
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer :: mx_apart,mx_pcell
    integer :: mx_ahalo,mx_part
    integer :: np_in_ahalo
    integer :: lab_ahcell(mx_apart),kbeg(mx_apart),nah_part(mx_apart)
    integer :: kseq(mx_ahalo),k_halo(mx_ahalo)

    ! Local variables
    integer :: kpart,k

    if((np_in_ahalo<=0).or.(np_in_ahalo>mx_apart)) then
       call cq_abort('check_mm: np_in_halo wrong')
    endif
    do kpart=1,np_in_ahalo
       if(lab_ahcell(kpart)>mx_pcell) then
          call cq_abort('check_mm: lab_ahcell wrong')
       endif
       if(kbeg(kpart)+nah_part(kpart)-1>mx_ahalo) then
          call cq_abort('check_mm: kseq arg wrong')
       endif
       if(nah_part(kpart)<=0) then
          call cq_abort('check_mm: nahpart < 1')
       endif
       do k=1,nah_part(kpart)
          if(kseq(kbeg(kpart)+k-1)>mx_part) then
             call cq_abort('check_mm: ibind_rem wrong')
          endif
          if(k_halo(kbeg(kpart)+k-1)>mx_ahalo) then
             call cq_abort('check_mm: nahnab arg wrong')
          endif
       enddo
    enddo
    return
  end subroutine check_mm
!!***
end module mult_init_module
