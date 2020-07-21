module process

  implicit none
  
contains

  subroutine assign_blocks

    use datatypes
    use local, ONLY: block_store, nprocs, block_size_x, block_size_y, block_size_z, &
         stm_z_min, stm_z_max, stm_x_min, stm_x_max, stm_y_min, stm_y_max
    
    implicit none

    integer :: proc, iblock, ig1, ind_group, block_x, block_y, block_z, nblock
    real(double) :: rbx, rby, rbz

    write(*,*) 'Range for x: ',stm_x_min, stm_x_max
    write(*,*) 'Range for y: ',stm_y_min, stm_y_max
    write(*,*) 'Range for z: ',stm_z_min, stm_z_max
    nblock = 0
    do proc=1,nprocs
       block_store(proc)%active = 0
       do iblock=1, block_store(proc)%num_blocks
          ! Create starting point for block using ind_block
          nblock = nblock + 1
          ! Check min area - if the RHS of block is in LHS of area, keep it
          rbx = block_size_x*real(block_store(proc)%nx(iblock),double)
          rby = block_size_y*real(block_store(proc)%ny(iblock),double)
          rbz = block_size_z*real(block_store(proc)%nz(iblock),double)
          if(rbx>=stm_x_min.AND.rby>=stm_y_min.AND.rbz>=stm_z_min) &
               block_store(proc)%active(iblock) = block_store(proc)%active(iblock) + 1
          !write(*,*) 'RHS: ',rbx,rby,rbz
          ! Now LHS of block and RHS of area
          rbx = block_size_x*real(block_store(proc)%nx(iblock)-1,double)
          rby = block_size_y*real(block_store(proc)%ny(iblock)-1,double)
          rbz = block_size_z*real(block_store(proc)%nz(iblock)-1,double) 
          if(rbx<=stm_x_max.AND.rby<=stm_y_max.AND.rbz<=stm_z_max) &
               block_store(proc)%active(iblock) = block_store(proc)%active(iblock) + 1
          !write(*,*) 'RHS: ',rbx,rby,rbz
          !write(*,*) 'Active: ',iblock!block_store(proc)%active(iblock)
          if(block_store(proc)%active(iblock)==2) then
             !write(*,*) 'Active: ',iblock!block_store(proc)%active(iblock)
             block_store(proc)%active(iblock)=1
          else
             block_store(proc)%active(iblock)=0
          end if
       end do
    end do
    write(*,*) "Blocks assigned"
    return
  end subroutine assign_blocks

  subroutine read_acc_charge

    use datatypes
    use numbers
    use local, ONLY: nprocs, n_bands_active, nkp, block_store, efermi, wtk, stm_broad, &
         nxmin, nymin, nzmin, current, nptsx, nptsy, nptsz, &
         eigenvalues, band_no, stm_bias
    use block_module, only: n_pts_in_block, in_block_x,in_block_y,in_block_z
    use io_module, ONLY: get_file_name

    implicit none

    integer :: proc, band, nk, iblock, point, ptx, pty, ptz, npx, npy, npz, idum1, idum2
    real(double) :: weight, rbx, rby, rbz, sq, test
    real(double), dimension(:), allocatable :: local_grid
    character(len=50) :: ci, filename

    allocate(current(nptsx,nptsy,nptsz))
    current = zero
    allocate(local_grid(n_pts_in_block))
    local_grid = zero
    do proc = 1, nprocs
       !write(*,*) 'Proc: ',proc
       ! Loop over active bands
       do band = 1, n_bands_active
          write(ci,'("Band",I0.6,"Kwfden")') band_no(band)
          call get_file_name(ci,nprocs,proc,filename)
          !write(*,*) 'Band, file: ',band, filename
          ! Open file
          open(unit=17,file=filename)
          read(17,*) idum1,idum2 ! nk, num points
          ! Read kpoints, grid points
          do nk = 1, nkp
             read(17,*) rbx,rby,rbz,sq ! kpt and energy
             !write(*,*) 'Kpoint, energy: ',nk,eigenvalues(nk,band)
             ! Read k-point and energy
             ! If energy lies in window, calculate weight
             if(eigenvalues(nk,band)>-1e20_double) then
                if(stm_bias<0) then
                   if(eigenvalues(nk,band)>efermi) then
                      sq = (eigenvalues(nk,band)-efermi)/stm_broad
                      weight = wtk(nk)/(exp(sq)+one)
                   else if(eigenvalues(nk,band)<(efermi+stm_bias)) then
                      sq = (eigenvalues(nk,band)-efermi-stm_bias)/stm_broad
                      weight = wtk(nk)/(exp(sq)+one)
                   else
                      weight = wtk(nk)
                   end if
                else
                   if(eigenvalues(nk,band)<efermi) then
                      sq = (eigenvalues(nk,band)-efermi)/stm_broad
                      weight = wtk(nk)/(exp(sq)+one)
                   else if(eigenvalues(nk,band)>(efermi+stm_bias)) then
                      sq = (eigenvalues(nk,band)-efermi-stm_bias)/stm_broad
                      weight = wtk(nk)/(exp(sq)+one)
                   else
                      weight = wtk(nk)
                   end if
                end if
                write(*,*) 'Weight: ',weight,eigenvalues(nk,band),nk,band
                call read_domain(17,local_grid,proc)
             else
                call read_domain(17,local_grid,proc,1)
             end if ! Eigenvalue in energy window
          end do ! nkp
          close(unit=17)
       end do ! band
    end do ! proc
    return
  end subroutine read_acc_charge

  subroutine read_write_charge

    use datatypes
    use numbers
    use local, ONLY: nprocs, n_bands_active, nkp, block_store, efermi, wtk, stm_broad, &
         nxmin, nymin, nzmin, current, nptsx, nptsy, nptsz, &
         eigenvalues, band_no, stm_bias, charge_stub, n_bands_active, band_no, nkp, flag_by_kpoint, &
         flag_output, dx, cube
    use block_module, only: n_pts_in_block, in_block_x,in_block_y,in_block_z
    use io_module, ONLY: get_file_name
    use output, ONLY: write_dx_density, write_cube, write_dx_coords
    use global_module, only : nspin

    implicit none

    integer :: proc, band, nk, iblock, point, ptx, pty, ptz, npx, npy, npz, idum1, idum2, kp, ispin
    real(double) :: weight, rbx, rby, rbz, sq, test
    real(double), dimension(:), allocatable :: local_grid
    character(len=50) :: filename, ci

    real(double), allocatable, dimension(:,:,:) :: tmp_current

    write(*,*) 'Processing charge'
    allocate(current(nptsx,nptsy,nptsz),tmp_current(nptsx,nptsy,nptsz))
    allocate(local_grid(n_pts_in_block))
    do ispin=1,nspin
       if(flag_by_kpoint) then ! Separate bands by k-point
          do band=1,n_bands_active
             tmp_current = zero
             ! Open files
             ! In future we might need to move this inside the loop inside k-points
             do proc = 1, nprocs
                ! Loop over active bands
                if(nspin==1) then
                   write(ci,'("Band",I0.6,"Kwfden")') band_no(band)
                else
                   if(ispin==1) write(ci,'("Band",I0.6,"Kwfden_up")') band_no(band)
                   if(ispin==2) write(ci,'("Band",I0.6,"Kwfden_dn")') band_no(band)
                end if
                call get_file_name(ci,nprocs,proc,filename)
                !call get_file_name(charge_stub,nprocs,proc,filename)
                ! Open file
                open(unit=17+proc,file=filename)
                read(17+proc,*) idum1,idum2 ! nk, num points
             end do
             do kp = 1,nkp
                current = zero
                do proc = 1, nprocs
                   read(17+proc,*) rbx,rby,rbz,sq ! kpt and energy
                   call read_domain(17+proc,local_grid,proc)
                end do ! proc
                write(ci,'("Band",I0.6,"Kwfden",I0.3)') band_no(band),kp
                if(flag_output==dx) then
                   call write_dx_density(ci)
                else
                   call write_cube(ci)
                end if
                ! Accumulate charge
                tmp_current = tmp_current + current
             end do ! kp
             current = tmp_current ! Total charge
             write(ci,'("Band",I0.6,"wfdenTot")') band_no(band)
             if(flag_output==dx) then
                call write_dx_density(ci)
             else
                call write_cube(ci)
             end if
             do proc = 1, nprocs
                close(unit=17+proc)
             end do
          end do ! bands = 1, n_bands_active
       else ! Not all k-points
          if(n_bands_active==0) then ! Just charge
             current = zero
             do proc = 1, nprocs
                ! Loop over active bands
                !write(ci,'("chden")') 
                !call get_file_name(ci,nprocs,proc,filename)
                if(nspin==1) then
                   ci = TRIM(charge_stub)
                   call get_file_name(ci,nprocs,proc,filename)
                else
                   if(ispin==1) ci = TRIM(charge_stub)//"_up"
                   if(ispin==2) ci = TRIM(charge_stub)//"_dn"
                   call get_file_name(ci,nprocs,proc,filename)
                end if
                ! Open file
                open(unit=17,file=filename)
                call read_domain(17,local_grid,proc)
                close(unit=17)
             end do ! proc
             if(flag_output==dx) then
                call write_dx_density(ci)
             else
                call write_cube(ci)
             end if
          else
             do band=1,n_bands_active
                current = zero
                do proc = 1, nprocs
                   ! Loop over active bands
                   if(nspin==1) then
                      write(ci,'("Band",I0.6,"wfden")') band_no(band)
                   else
                      if(ispin==1) write(ci,'("Band",I0.6,"wfden_up")') band_no(band)
                      if(ispin==2) write(ci,'("Band",I0.6,"wfden_dn")') band_no(band)
                   end if
                   call get_file_name(ci,nprocs,proc,filename)
                   !call get_file_name(charge_stub,nprocs,proc,filename)
                   ! Open file
                   open(unit=17,file=filename)
                   call read_domain(17,local_grid,proc)
                   close(unit=17)
                end do ! proc
                if(flag_output==dx) then
                   call write_dx_density(ci)
                else
                   call write_cube(ci)
                end if
             end do ! bands
          end if
       end if
    end do
    if(flag_output==dx) call write_dx_coords(charge_stub)
    return
  end subroutine read_write_charge

  subroutine read_domain(lun,local_grid,proc,dummy)

    use datatypes
    use numbers
    use local, ONLY: block_store, nxmin, nymin, nzmin, current, nptsx, nptsy, nptsz
    use block_module, only: n_pts_in_block, in_block_x,in_block_y,in_block_z
    
    implicit none

    ! Passed
    integer :: lun, proc
    real(double), dimension(n_pts_in_block) :: local_grid
    integer, OPTIONAL :: dummy

    ! Local
    integer :: iblock, point, ptx, pty, ptz, npx, npy, npz
    logical :: is_dummy

    is_dummy = .false.
    if(PRESENT(dummy)) then
       if(dummy==1) is_dummy = .true.
    end if
    do iblock=1,block_store(proc)%num_blocks
       local_grid = zero
       ! Read block
       do point = 1,n_pts_in_block
          read(lun,*) local_grid(point)
       end do
       if(.NOT.is_dummy) then
          point = 0
          if(block_store(proc)%active(iblock)==1) then
             do ptz = 1, in_block_z
                npz = in_block_z*(block_store(proc)%nz(iblock)-1) + ptz - nzmin
                if(npz>0.AND.npz<=nptsz) then
                   do pty = 1, in_block_y
                      npy = in_block_y*(block_store(proc)%ny(iblock)-1) + pty - nymin
                      if(npy>0.AND.npy<=nptsy) then
                         do ptx = 1, in_block_x
                            npx = in_block_x*(block_store(proc)%nx(iblock)-1) + ptx - nxmin
                            point = point + 1
                            if(npx>0.AND.npx<=nptsx) then
                               current(npx,npy,npz) = current(npx,npy,npz) + local_grid(point)
                            end if
                         end do
                      else
                         do ptx = 1, in_block_x
                            point = point + 1
                         end do
                      end if
                   end do
                else
                   do pty = 1, in_block_y
                      do ptx = 1, in_block_x
                         point = point + 1
                      end do
                   end do
                end if
             end do
          end if
       end if ! not is_dummy
    end do ! iblock
  end subroutine read_domain
  
end module process
