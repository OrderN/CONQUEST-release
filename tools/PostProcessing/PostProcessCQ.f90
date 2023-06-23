! Read and convert Conquest charge density files into files that can be read either by
! OpenDX or Cube files (other formats ?)
program ConvertCharge

  use GenComms, ONLY: myid, numprocs, init_comms, end_comms
  use read
  use process
  use local
  use output
  use memory_module,     only: init_reg_mem
  use stm
  
  implicit none

  call init_comms(myid,numprocs)
  call init_reg_mem

  call write_banner
  ! Load input and read
  call read_input
  ! Now select job
  if(i_job==1) then ! Convert coordinates
     if(coord_format==1) then
        write(*,fmt='(2x,"Converting coordinates to XYZ format"/)')
        call write_xyz
     else if(coord_format==2) then
        write(*,fmt='(2x,"Converting coordinates to CASTEP cell format"/)')
        call write_cell
     else if(coord_format==3) then
        write(*,fmt='(2x,"Converting coordinates to XSF format"/)')
        call write_xsf
     end if
  else if(i_job==2) then ! Convert chden to cube
     write(*,fmt='(2x,"Converting charge density to cube format"/)')
     ! Read blocks
     call read_block_input
     ! Calculate active processes/blocks
     call assign_blocks
     call process_charge
  else if(i_job==3) then ! Create band-resolved charge densities
     write(*,fmt='(2x,"Creating band-resolved charge densities in cube format"/)')
     ! Read blocks to get grid
     call read_block_input
     call process_bands
  else if(i_job==4) then ! STM via TH sum over band densities
     write(*,fmt='(2x,"Creating Tersoff-Hamman simulated STM image in cube format"/)')
     ! Read blocks to get grid
     call read_block_input
     call tersoff_hamann
  else if(i_job==5) then ! STM via Paz-Soler TH
     write(*,*) "Paz-Soler simulation is still under development and not available by default"
     if(.false.) then
        ! Read blocks
        call read_block_input
        ! Calculate active processes/blocks
        call assign_blocks
        call make_stm_current
     end if
  else if(i_job==6) then ! DOS output
     write(*,fmt='(2x,"Creating density of states (DOS)"/)')
     call process_dos
  else if(i_job==7) then ! pDOS output
     write(*,fmt='(2x,"Creating projected density of states (pDOS); also generating DOS"/)')
     call process_dos
     call process_pdos
  else if(i_job==8) then ! DOS output
     write(*,fmt='(2x,"Creating band structure"/)')
     call process_band_structure
  end if
  ! Read eigenvalues and calculate weight for bands by kpt
  !if(flag_only_charge) then
  !   if(flag_range.and.n_bands_active==0) call read_eigenvalues_bands
  !   call read_write_charge
  !else
  !   call read_eigenvalues_stm
  !   ! STM output
  !   call make_stm_current
  !   call end_comms
  !   return
  !   ! Read charge and accumulate
  !   call read_acc_charge
  !   ! Write output files
  !   !if(flag_output==dx) then
  !   !   call write_dx_coords
  !   !   call write_dx_density
  !   !else
  !   !   call write_cube
  !   !end if
  !end if
  call end_comms
end program ConvertCharge
  
!   
!   integer, parameter :: prec = selected_real_kind(12,50)
! 
!   real(prec), dimension(:,:,:), allocatable :: density
!   real(prec) :: test, sx, sy, sz
!   integer, dimension(:), allocatable :: first_block
!   integer :: nnodes, blocksx, blocksy, blocksz, pointsx, pointsy, pointsz
!   integer :: n_blocks, n_per_node, n_first, n_points, stat, nnd, locx, locy, locz
!   integer :: nfb_x, nfb_y, nfb_z, nc_x, nc_y, nc_z, n_add_x, n_add_y, n_add_z, nx, ny, nz
!   integer :: nhund, ntens, nunit, n1, icount, inode, ptx, pty, ptz, blocks_on_proc, isx, isy, isz, samp
!   integer :: block_x, block_y, block_z,iblock,ind_group,i,proc,blocksx1,blocksy1,blocksz1,nnodes1,ig1,proc1,idum
!   logical :: flag_read_blocks
!   character(len=11) :: digitstr = "01234567890"
!   character(len=50) :: filename
!   character(len=50) :: stub
! 
!   open(15,file="convert_charge.dat")
!   read(15,*) nnodes
!   read(15,*) blocksx, blocksy, blocksz
!   read(15,*) pointsx, pointsy, pointsz
!   read(15,*) flag_read_blocks
!   read(15,*) sx,sy,sz
!   read(15,*) samp
!   read(15,*) stub
!   isx = sx*blocksx*pointsx
!   isy = sy*blocksy*pointsy
!   isz = sz*blocksz*pointsz
!   write(*,*) 'Blocks: ',blocksx, blocksy, blocksz
!   write(*,*) 'Points: ',pointsx, pointsy, pointsz
!   write(*,*) 'Shift: ',isx, isy,isz
!   write(*,*) 'Sampling: ',samp
!   n_points = pointsx*blocksx*pointsy*blocksy*pointsz*blocksz
!   write(*,*) 'Total points: ',n_points
!   if(flag_read_blocks) then
!      write(*,*) 'Reading make_blk.dat for blocks'
!   else
!      write(*,*) 'Old-style blocks'
!   end if
!   allocate(density(pointsx*blocksx,pointsy*blocksy,pointsz*blocksz),STAT=stat)
!   if(stat/=0) write(*,*) 'Error allocating memory'
!   ! At this point we fork into old- and new-style blocks
!   if(flag_read_blocks) then
!      ! Open make_blk.dat
!      open(unit=17,file='make_blk.dat')
!      read(17,*) blocksx1,blocksy1,blocksz1
!      if(blocksx1/=blocksx.OR.blocksy1/=blocksy.OR.blocksz1/=blocksz) then
!         write(*,*) 'Block specification error: ',blocksx,blocksx1,blocksy,blocksy1,blocksz,blocksz1
!         stop
!      end if
!      read(17,*) nnodes1
!      if(nnodes1/=nnodes) then
!         write(*,*) 'Processors error: ',nnodes,nnodes1
!         stop
!      end if
!      ! Loop over processors
!      do proc = 1,nnodes
!         ! Read from make_blk
!         read(17,*) proc1,blocks_on_proc,idum
!         write(*,*) 'Processor, blocks: ',proc,blocks_on_proc
!         ! For each processor, open the file
!         nhund = aint(real(proc/100))+1
!         n1 = proc - (100*nhund) + 100
!         ntens = aint(real(n1/10))+1
!         nunit = n1 - 10*ntens + 11
!         filename = trim(stub)//'.'//digitstr(Nhund:Nhund)//digitstr(Ntens:Ntens)//digitstr(Nunit:Nunit)
!         !filename = 'chden.'//digitstr(Nhund:Nhund)//digitstr(Ntens:Ntens)//digitstr(Nunit:Nunit)
!         open(unit=16,file=filename)
!         ! Now, loop over grid points on processor and read in charge
!         do iblock=1, blocks_on_proc
!            ! Define ind_group: does this come from make_blk.dat ?
!            read(17,*) ig1,ind_group
!            ! Create starting point for block using ind_block
!            block_x = 1+(ind_group-1)/(blocksy*blocksz)
!            block_y = 1+(ind_group-1-(block_x-1)*blocksy*blocksz)/blocksz
!            block_z = ind_group - (block_x-1)*blocksy*blocksz - (block_y-1)*blocksz
!            do ptz = 1,pointsz
!               do pty = 1,pointsy
!                  do ptx = 1,pointsx
!                     locx = (block_x-1)*pointsx + ptx + isx
!                     if(locx<1) locx = locx + blocksx*pointsx
!                     if(locx>blocksx*pointsx) locx = locx - blocksx*pointsx
!                     locy = (block_y-1)*pointsy + pty + isy
!                     if(locy<1) locy = locy + blocksy*pointsy
!                     if(locy>blocksy*pointsy) locy = locy - blocksy*pointsy
!                     locz = (block_z-1)*pointsz + ptz + isz
!                     if(locz<1) locz = locz + blocksz*pointsz
!                     if(locz>blocksz*pointsz) locz = locz - blocksz*pointsz
!                     read(16,*) test
!                     density(locx,locy,locz) = test
!                  end do
!               end do
!            end do
!         end do ! iblock=domain%groups_on_node
!         ! Close file
!         close(16)
!      end do
!   else
!      allocate(first_block(nnodes))
!      n_per_node = blocksx*blocksy*blocksz/nnodes
!      n_first = blocksx*blocksy*blocksz - n_per_node*nnodes
!      if (n_first.ne.0) then
!         write(*,*) 'the first ', n_first, ' nodes each have ', &
!              n_per_node + 1, ' blocks;'
!         write(*,*) 'the remaining ', nnodes - n_first,          &
!              ' nodes each have ', n_per_node, ' blocks.'
!      else
!         write(*,*) 'All nodes have ', n_per_node, ' blocks.'
!      end if
!      do nnd=1, nnodes
!         if(nnd == 1) then
!            first_block(nnd)=1
!         else
!            first_block(nnd)= first_block(nnd-1)+n_per_node
!            if(nnd<n_first) first_block(nnd) = first_block(nnd)+1
!         endif
!      enddo
!      nfb_x = 1
!      nfb_y = 1
!      nfb_z = 1
!      nc_x = 1
!      nc_y = 1
!      nc_z = 1
!      n_add_z = 0
!      icount = 0
!      inode = 1
!      nhund = aint(real(inode/100))+1
!      n1 = inode - (100*nhund) + 100
!      ntens = aint(real(n1/10))+1
!      nunit = n1 - 10*ntens + 11
!      filename = trim(stub)//'.'//digitstr(Nhund:Nhund)//digitstr(Ntens:Ntens)//digitstr(Nunit:Nunit)
!      !filename = 'chden.'//digitstr(Nhund:Nhund)//digitstr(Ntens:Ntens)//digitstr(Nunit:Nunit)
!      open(unit=16,file=filename)
!      !read(16,*) test
!      !write(*,*) 'Test is : ',test
!      do nz = 1, blocksz
!         nc_z = nc_z + n_add_z
!         n_add_z = nfb_z
!         n_add_y = 0
!         do ny = 1, blocksy
!            nc_y = nc_y + n_add_y
!            n_add_y = nfb_y
! 
!            n_add_x = 0
!            do nx = 1, blocksx
!               nc_x = nc_x + n_add_x
!               n_add_x = nfb_x
!               icount=icount+1
!               !write(*,*) 'Block xyz is ',nc_x, nc_y, nc_z
!               ! Check - is this a new file ?
!               if(inode<nnodes.AND.icount==first_block(inode+1)) then
!                  inode = inode+1
!                  write(*,*) 'Moving to a new node: ',inode
!                  close(16)
!                  nhund = aint(real(inode/100))+1
!                  n1 = inode - (100*nhund) + 100
!                  ntens = aint(real(n1/10))+1
!                  nunit = n1 - 10*ntens + 11
!                  filename = 'chden.'//digitstr(Nhund:Nhund)//digitstr(Ntens:Ntens)//digitstr(Nunit:Nunit)
!                  open(unit=16,file=filename)
!                  !read(16,*) test
!                  !write(*,*) 'Test is : ',test
!               end if
!               do ptz = 1,pointsz
!                  do pty = 1,pointsy
!                     do ptx = 1,pointsx
!                        locx = (nc_x-1)*pointsx + ptx + isx
!                        if(locx<1) locx = locx + blocksx*pointsx
!                        if(locx>blocksx*pointsx) locx = locx - blocksx*pointsx
!                        locy = (nc_y-1)*pointsy + pty + isy
!                        if(locy<1) locy = locy + blocksy*pointsy
!                        if(locy>blocksy*pointsy) locy = locy - blocksy*pointsy
!                        locz = (nc_z-1)*pointsz + ptz + isz
!                        if(locz<1) locz = locz + blocksz*pointsz
!                        if(locz>blocksz*pointsz) locz = locz - blocksz*pointsz
!                        read(16,*) test
!                        density(locx,locy,locz) = test
!                     end do
!                  end do
!               end do
!               ! Loop over block grid points
!            enddo
!            nfb_x = -nfb_x
! 
!         enddo
!         nfb_y = -nfb_y
! 
!      enddo
!      close(16)
!   end if ! if flag_read_blocks else
!   open(unit=16,file="Charge.dat")
!   write(16,fmt='(2x,3i6)') pointsx*blocksx/samp,pointsy*blocksy/samp,pointsz*blocksz/samp
!   icount = 0
!   do ptz = 1,pointsz*blocksz,samp
!      do pty = 1,pointsy*blocksy,samp
!         do ptx = 1,pointsx*blocksx,samp
!            icount = icount+1
!            if(mod(icount,5)==0) then
!               write(16,fmt='(e18.9)') density(ptx,pty,ptz)
!            else
!               write(16,fmt='(e18.9)',advance='no') density(ptx,pty,ptz)
!            end if
!         end do
!      end do
!   end do
!   close(16)
!   if(.NOT.flag_read_blocks) deallocate(first_block)
!   deallocate(density)
! end program ConvertCharge
