module dynpro_sub

  use kinds,      only : DP
  use constants,  only : bohr => BOHR_RADIUS_ANGS

  implicit none

  save

contains

  subroutine dynpro_sub_file(xyz_)

    use dynpro_inp, only : nat, label, ityp
    use dynpro_inp, only : filesub_out, sub_nsp  

    implicit none

    real(DP), intent(in)  :: xyz_(3,nat)

    integer :: i, j

    !ALLOCATE(xyz_xml(3,nat_xml))

    !xyz_xml = MATMUL(ht_xml,stau_xml)

    open(911, file=filesub_out, status='unknown', position='rewind')

    write(911,*) nat
    write(911,*) sub_nsp 

    !DO i = 1, nat 
    !   WRITE(911,'(a2,x,3f12.4,2x,i2,x,3i)')  &
    !        label(i), (xyz_xml(j,i)*bohr, j=1,3), ityp(i), i     
    !END DO

    do i = 1, nat 
       write(911,*)  &
            label(i), (xyz_(j,i), j=1,3), ityp(i), i     
    end do

    !DEALLOCATE(xyz_xml)

    close(911)

  end subroutine dynpro_sub_file


  subroutine dynpro_sub_atom()

    use dynpro_inp, only : ldebug        
    use dynpro_inp, only : nat, nsp, na, label, ityp, atm, label, atomic_number, &
         tmp_atm, syst_ID_spe, syst_ID_lab, atm_xml, ct_atom, sph_atom
    use dynpro_inp, only : iteration
    use dynpro_inp, only : ityp_old, label_old
    use dynpro_inp, only : perm_xyz, xyz_sub, perm_xyz_tmp, perm_n
    use dynpro_inp, only : filesub_in, sub_nsp
    use dynpro_inp, only : xyz, vel
    
    use dynpro_mod, only : permute, very_tricky, tricky_sort

    
    implicit none

    integer :: i, j

    !ALLOCATE(ityp_old(nat))
    !ALLOCATE(label_old(nat))

    ityp_old  = ityp
    !nsp_old   = nsp
    label_old = label

    deallocate(ityp,label,atm,tmp_atm,na,atomic_number,syst_ID_spe,syst_ID_lab)

    allocate(ityp(nat))
    allocate(label(nat))
    allocate(xyz_sub(3,nat))


    open(912, file=filesub_in, status='old', position='rewind')

    read(912,*)
    read(912,*) sub_nsp

    !DO i = 1, nat
    !   READ(912,'(a2,x,3f12.4,2x,i2,x,3i)') &
    !        label(i), (xyz_sub(j,i), j=1,3), ityp(i)     
    !END DO

    do i = 1, nat
       read(912,*) &
            label(i), (xyz_sub(j,i), j=1,3), ityp(i)
       if (ldebug) print*, label(i), (xyz_sub(j,i), j=1,3), ityp(i)
    end do


    close(912)
    !---------!

    ityp_old  = ityp
    !label_old = label

    allocate(atm(nsp+sub_nsp))
    allocate(tmp_atm(nsp+sub_nsp))
    allocate(na(nsp+sub_nsp))
    allocate(atomic_number(nsp+sub_nsp))
    allocate(syst_ID_spe(nsp+sub_nsp,2))
    allocate(syst_ID_lab(nsp+sub_nsp))

    atm(1:nsp) = atm_xml(1:nsp)

    if (ldebug) print*, 'call very_tricky()', perm_n

    call VERY_TRICKY(sub_nsp, nsp+sub_nsp, nat, atm, ityp, label_old, label, &
         syst_ID_spe, syst_ID_lab, perm_xyz_tmp, perm_n)

    allocate(perm_xyz(perm_n))

    perm_xyz(:) = perm_xyz_tmp(1:perm_n)

    nsp = nsp + sub_nsp

    if (ldebug) write(*,*) ' call permute()'
    if (ldebug) write(*,*) perm_xyz(:) 
    if (ldebug) write(*,*) syst_ID_spe
    if (ldebug) write(*,*) syst_ID_lab
    if (ldebug) write(*,*) label
    if (ldebug) write(*,*)
    if (ldebug) write(*,*) nsp
    if (ldebug) write(*,*)
    if (ldebug) write(*,*) ityp_old
    if (ldebug) write(*,*) 
    if (ldebug) write(*,*) ityp
    
    !CALL PERMUTE(nat, iteration, xyz, perm_xyz, perm_n)
    !CALL PERMUTE(nat, iteration, vel, perm_xyz, perm_n)
    !CALL PERMUTE(nat, iteration, for, perm_xyz, perm_n)

    !nsp = nsp + sub_nsp

    call tricky_sort(nsp, nat, iteration, xyz, vel, ityp_old, ityp, syst_ID_spe, &
         label_old, syst_ID_lab, .false.)

    ct_atom  = ct_atom  - perm_n
    sph_atom = sph_atom - perm_n

    deallocate(xyz_sub, perm_xyz, label_old, ityp_old)


    !print*, 'I am here'


  end subroutine dynpro_sub_atom

end module dynpro_sub
