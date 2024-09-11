! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id: exx_module.f90 XX year-mm-dd XX:XX:XXZ Lionel $
! -----------------------------------------------------------
! Module exx_module
! -----------------------------------------------------------
! Code area 13: EXX
! -----------------------------------------------------------

!!***h* Conquest/exx_module *
!!  NAME
!!   exx_module
!! 
!!  PURPOSE
!!   Holds routines called in exx_kernel_module
!!   to compute exact exchange
!!
!!  USES
!!   GenComms, dimens, datatypes, numbers, units 
!!   group_module, cover_module, global_module
!!   species_module, primary_module, matrix_data
!!   atomic_density, support_spec_format
!!   fft_interface_module, exx_types, exx_io 
!!   DiagModule
!!
!!  AUTHOR
!!   L.A.Truflandier (lat)
!!
!!  CREATION DATE
!!   2011/02/11
!!
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
module exx_module
  
  use datatypes
  use GenComms,                  only: my_barrier, cq_abort, mtime
  use GenComms,                  only: inode, ionode, myid, root
  use timer_module,              only: start_backtrace, stop_backtrace, cq_timer
  !**<lat>** ISF Poisson solver Will be available in the forthcoming version 
  !use Poisson_Solver,            only: PSolver, createKernel, gequad
  use global_module,             only: io_lun, exx_cutoff
  use exx_types,                 only: reckernel_3d, fftwrho3d, exx_debug
  use exx_io
  
  use fft_interface_module,      only: fft3_exec_wrapper, fft3_init_wrapper
  
  implicit none 
  
  ! Area identification
  integer, parameter, private :: area = 13
  
  !!***
  
contains
  !
  subroutine get_X_params(level)
    
    use exx_types, only: exx_psolver,exx_pscheme,exx_pscheme_default,p_cutoff,p_factor,p_omega
    use exx_types, only: pulay_factor,pulay_radius, magic_number,ewald_alpha,isf_order
    use exx_types, only: extent,ngrid,r_int,grid_spacing,volume,edge,screen

    use exx_types, only: exx_scheme, exx_mem, exx_screen, exx_alloc
    use exx_types, only: exx_cartesian, exx_overlap, exx_radius, exx_screen_pao, exx_hgrid
    use exx_types, only: exx_gto, exx_debug
    use exx_types, only: tmr_std_exx_setup, exx_store_eris

    use atomic_density, only: atomic_density_table
    use species_module, only: n_species

    use global_module,  only: iprint_exx, exx_hgrid_medium, exx_hgrid_coarse, exx_alpha
    use numbers,        only: zero, very_small, one, two, three
    use units,          only: BohrToAng
    use dimens,         only: r_super_x, r_super_y, r_super_z, &
                              n_grid_x, n_grid_y, n_grid_z,    &
                              r_h, r_nl, r_h   
    implicit none

    !real(double), optional, intent(in) :: exx_scf_hgrid
    integer, optional :: level

    character(100) :: solver, scheme
    character(100) :: phil_scheme, eri_scheme, alloc_scheme, pao_scheme
    character(100) ::  mem_scheme, scr_scheme
    character(100) :: filename
    real(double)   :: threshold = 1.0d-4
    real(double)   :: r_shift   = 0.0_double
    real(double)   :: gs(3)     = 0.0_double
    real(double)   :: r_max, r_maxs, gs_min
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level
    integer        :: i,j,k,l,unit
   
    !
!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='get_X_params',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$
    !
    call start_timer(tmr_std_exx_setup) 
    !
    ! **<lat>** This routine need to be re-written !
    !
    exx_screen      = .false.
    exx_screen_pao  = .false.
    exx_store_eris  = .false.
    !exx_cutoff      = zero              ! do not touch  
    !
    ! Find out the finest grid spacing from input
    ! Input grid spacing from input (Bohr unit)
    gs(1) = r_super_x/n_grid_x
    gs(2) = r_super_y/n_grid_y
    gs(3) = r_super_z/n_grid_z
    !
    ! Find out the finest grid spacing from input
    ! we use the same grid spacing for Ox,Oy,Oz   
    gs_min = minval(gs)
    !
    ! Setup grid spacing of EXX
    if ( exx_hgrid < very_small ) then
       grid_spacing = gs_min
    else if ( exx_hgrid >= zero ) then
       grid_spacing = exx_hgrid       
    else
       call cq_abort('EXX: unrecognised grid_spacing ',exx_hgrid)
       !
    end if
    !
    !
    ! **<lat>** Need to be adapted to GPFA ; not too difficult
    ! Find r_max for integration
    r_max = zero
    do i = 1, n_species
       r_max = max(r_max,atomic_density_table(i)%cutoff)
    end do
    !
    ! Setup EXX r_int
    if (exx_radius > very_small) then
       r_int = exx_radius
    else
       r_int = real(ceiling(r_max),double)
    end if
    !
    ! Number of grid points alongthe ]O,X] segment
    extent = ceiling(r_int/grid_spacing)    
    !
    ! Re-compute the real grid spacing with respect to extent
    grid_spacing = r_int/real(extent,double)
    !
    ! Number of grid points along [-X,X] segment
    ngrid  = 2*extent+1
    edge   = two*r_int 
    volume = edge**3    
    !
    ! Setup Poisson solver
    if (exx_psolver == 'fftw') then 
       solver = 'reciprocal space'
       scheme = trim(exx_psolver)     
    else 
       exx_psolver = 'isf'
       solver   = 'real space'
       scheme   = trim(exx_psolver)
       if (isf_order <= 0) then
          isf_order = 8
       end if
    end if
    !
    ! Below for output purpose
    if (exx_scheme==1) then
       eri_scheme = "3center reduction integrals"
    else if ( exx_scheme==2 ) then
       eri_scheme = "4center eris"
    else if ( exx_scheme==3 ) then
       eri_scheme = "4center eris + storage"
    else
       eri_scheme  = "none"
       exx_psolver = "none"
       solver = 'none'
       scheme = trim(exx_psolver) 
    end if
    !
    !if (exx_phil) then
    !   phil_scheme = 'storage'
    !else
    !   phil_scheme = 'on-the-fly'
    !end if
    !
    if (exx_alloc) then
       alloc_scheme = 'on-the-fly'
    else
       alloc_scheme = 'once'
    end if
    !
    if (exx_cartesian) then
       pao_scheme = 'cartesian'
    else
       pao_scheme = 'spherical'
    end if
    !
    !if (exx_mem == 1) then
    !   mem_scheme = 'high'
    !else if (exx_mem == 2) then
    !   mem_scheme = 'low'
    !else if (exx_mem == 3) then
    !   mem_scheme = 'medium'
    !else
    !   mem_scheme = 'low'
    !   exx_mem    = 2
    !end if
    !
    if (exx_screen) then
       if (exx_cutoff > very_small) then
          scr_scheme     = 'user'
          exx_screen_pao = .false.
          screen         = exx_cutoff
       else
          exx_screen_pao =.true.
          scr_scheme = 'pao'
       end if
    else
       exx_screen_pao = .false.
       scr_scheme = 'n/a'
       screen     = 100.0d0
    end if
    !
    if (exx_scheme==2) then
       phil_scheme  = 'on-the-fly'
       alloc_scheme = 'on-the-fly'
       !mem_scheme   = 'high'
    end if
    !
    if ( inode == ionode .and. iprint_exx > 3 ) then
       write(io_lun,2) ('Entering in the EXX module')
       write(io_lun,41) eri_scheme
       !write(io_lun,42) phil_scheme
       write(io_lun,43) alloc_scheme
       !write(io_lun,48) mem_scheme
       write(io_lun,44) exx_screen, scr_scheme, screen, screen*BohrToAng
       
       write(io_lun,46) exx_overlap
       write(io_lun,47) pao_scheme
       write(io_lun,20)
       write(io_lun,21) r_max, r_max*BohrToAng
       write(io_lun,22) r_int, r_int*BohrToAng
       write(io_lun,23) grid_spacing, grid_spacing*BohrToAng
       write(io_lun,24) ngrid
       write(io_lun,25) 
       write(io_lun,26) edge, edge*BohrToAng
       write(io_lun,27) volume, volume*BohrToAng**3
       write(io_lun,28) ngrid**3

    end if

    if (exx_psolver == 'fftw') then       
       scheme = trim(scheme)//'/'//trim(exx_pscheme)
       poisson_fftw: select case(exx_pscheme)          
          
       case('default')
          scheme = trim(scheme)//trim(exx_pscheme_default)
          if ( inode == ionode .and. iprint_exx > 2 ) then
             write(io_lun,30) solver, scheme
             write(io_lun,31) ('G=0 component neglected... warning: inaccurate !')
          end if
       case('ewald')
          if (ewald_alpha < very_small) then
             ewald_alpha = real(3.0,double)
          end if
          if ( inode == ionode .and. iprint_exx > 2 ) then
             write(io_lun,30) solver, scheme
             write(io_lun,32) ewald_alpha
          end if
          
       case('pulay')
          if (pulay_radius < very_small) then
             !pulay_radius = r_maxs*sqrt(three)
             pulay_radius = r_int*sqrt(three)
             !pulay_radius = exx_radius

          end if
          if (pulay_factor < very_small) then
             pulay_factor = one
          end if 
          pulay_radius = pulay_factor*pulay_radius
          if ( inode == ionode .and. iprint_exx > 2 ) then
             write(io_lun,30) solver, scheme
             write(io_lun,33) pulay_factor
             write(io_lun,34) pulay_radius
          end if
       case('yukawa')
          if (p_omega < very_small) then
             p_omega = magic_number
             p_omega = -log(threshold*r_int)/r_int
          end if
          if ( inode == ionode .and. iprint_exx > 2 ) then
             write(io_lun,30) solver, scheme
             write(io_lun,35) p_omega
          end if
       case('gauss')
          if ( inode == ionode  .and. iprint_exx > 2 ) then
             write(io_lun,30) solver, scheme
          end if
       case default
          scheme = trim(exx_pscheme_default)
          if ( inode == ionode  .and. iprint_exx > 2 ) then
             write(io_lun,30) solver, scheme
             write(io_lun,31) ('WARNING: G=0 component neglected !')
          end if
          
       end select poisson_fftw
       
    else if (exx_psolver == 'isf') then
       if ( inode == ionode  .and. iprint_exx > 2 ) then
          write(io_lun,30) solver, scheme
          write(io_lun,36) isf_order
       end if
    end if
   
!****lat<$
    if ( inode == ionode .and. iprint_exx > 3 ) then
       write (io_lun,50) &
            grid_spacing,     &
            r_int,            &
            extent,           &
            trim(exx_psolver),&
            exx_alpha 
       
    end if

50 format (' EXX: gs = ',f6.4,', rc = ',f7.4,', extent = ',i4,', psolv = ',a,', alpha = ', f4.2)
!****lat>$
    !
    !
    call stop_timer(tmr_std_exx_setup,.true.)
    !
!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_X_matrix',echo=.true.)
!****lat>$
    !
    
    return

    !! EXX Formats start here
1   format(/1x,104a/)
2   format(4x,34a,8x,/)
    
    !! Grid settings
20  format(/25x,'Grid Settings: Cubic')  
21  format( 29x,'r_max_pao: ',    f12.4,' a0,',f12.4,' Ang')  
22  format( 25x,'r_integration: ',f12.4,' a0,',f12.4,' Ang')  
23  format( 26x,'grid_spacing: ', f12.8,' a0,',f12.8,' Ang')  
24  format( 25x,'n_grid_points: ',i4)  
    
25  format(/25x,'Local FFT Box: Cubic')  
26  format( 27x,'edge_length: ',f12.4,' a0,',f12.4,' Ang')
27  format( 32x,'volume: ',     f12.4,' a0,',f12.4,' Ang')
28  format( 23x,'n^3_grid_points: ',i12)  
    
    !! Poisson solver
30  format(/20x,'EXX Poisson Solver: ',a20,/32x,'Scheme: ',a20)  
    ! Warning
31  format( 24x,a34) 
    ! Ewald scheme
32  format( 33x,'alpha: ',f8.6,' a0^{-1/2}')
    ! Pulay scheme
33  format( 27x,'cutoff_fact: ',f4.1, ' a0')
34  format( 28x,'cutoff_pot: ',f4.1, ' a0')
    ! Yukawa scheme
35  format( 33x,'omega: ',f8.6, ' a0^{-1}')
    ! ISF Psolver
36  format( 31x,'n_order: ',i3)

    !! User settings
40  format(/25x,'Exact exchange settings:')  
41  format( 32x,'method: ',a30)
42  format( 34x,'phil: ',a30)
43  format( 28x,'allocation: ',a30)
44  format( 29x,'screening: ',l1,', method: ',a4, ', cutoff: ',f12.4,' a0,',f12.4,' Ang')

46  format( 27x,'overlap_box: ',l1)
47  format( 27x,'pao_on_grid: ',a30)
48  format( 24x,'memory_storage: ',a30)
    
  end subroutine get_X_params
  !
  !
  subroutine get_halodat(hl,kl,ind,part_cover,part,which,verbose,unit)

    use numbers,        only: three
    use group_module,   only: parts
    use cover_module,   only: BCS_parts
    use global_module,  only: id_glob,atom_coord, species_glob
    use species_module, only: species_label, nsf_species
    use atomic_density, only: atomic_density_table
    use numbers,        only: zero, one, two, twopi, pi
    !
    use matrix_data, only: mat, Hrange, Srange
    use exx_types,   only: neigh_atomic_data
    use exx_types,   only: tmr_std_exx_fetch
    use pao_format,  only: pao
    !
    use species_module, only: nsf_species
    !
    implicit none
    !
    type(neigh_atomic_data), intent(inout) :: hl
    type(neigh_atomic_data), intent(inout) :: kl
    integer,          intent(in)           :: ind
    integer,          intent(in)           :: part_cover
    integer,          intent(in)           :: part
    logical,          intent(in)           :: verbose
    character(len=1), intent(in)           :: which
    !
    real(double)      :: xyz_Ang(3)      
    character(len=20) :: filename
    integer           :: l1, acz1, m1, count, max_nsup
    integer, optional :: unit
    !
    call start_timer(tmr_std_exx_fetch)
    !
    max_nsup = maxval(nsf_species) 
    !
    ! Get u(v) data
    hl%npc = part_cover
    hl%nic = ind
    !
    hl%global_part = parts%icell_beg(part)
    hl%global_num  = id_glob(parts%icell_beg(part)+ind-1)
    !
    !if (which == 'k') then
       !hl%spec  = BCS_parts%spec_cover(BCS_parts%icover_ibeg(hl%npc)+hl%nic-1)
    hl%spec  = species_glob(hl%global_num)
    !else
    !   hl%spec = BCS_parts%spec_cover(part_cover+ind-1)
    !end if
    !
    !hl%spec     = BCS_parts%spec_cover(BCS_parts%icover_ibeg(hl%npc)+hl%nic-1)
    !
    hl%name = species_label(hl%spec)               
    hl%radi = atomic_density_table(hl%spec)%cutoff
    hl%nsup = nsf_species(hl%spec)
    !
    ! Calculate R_u
    hl%xyz_hl(1) = atom_coord(1,hl%global_num)
    hl%xyz_hl(2) = atom_coord(2,hl%global_num)
    hl%xyz_hl(3) = atom_coord(3,hl%global_num)
    !
    hl%xyz_cv(1) = BCS_parts%xcover(part_cover+ind-1)
    hl%xyz_cv(2) = BCS_parts%ycover(part_cover+ind-1)
    hl%xyz_cv(3) = BCS_parts%zcover(part_cover+ind-1)
    !
    ! Calculate R_iu
    hl%xyz(1) = - hl%xyz_cv(1) + kl%xyz_hl(1)
    hl%xyz(2) = - hl%xyz_cv(2) + kl%xyz_hl(2)
    hl%xyz(3) = - hl%xyz_cv(3) + kl%xyz_hl(3)
    hl%r      = sqrt(dot_product(hl%xyz,hl%xyz))
    !
    ! Calculate D_iu
    hl%d      = sqrt(three)*(hl%radi + kl%radi)
    !
    xyz_Ang(1) = hl%xyz_hl(1)
    xyz_Ang(2) = hl%xyz_hl(2)
    xyz_Ang(3) = hl%xyz_hl(3)
    !
    if ( exx_debug ) then
       if (which == 'k') then
          write(unit,'(I8,6X,A9,2X,A,I8,A2,I3,A,6X,A2,1X,I8,1X,4F12.4,3X,I3,2I8,1X,F7.3)')     &
               part,'{k\gamma}','{',hl%global_num,'\ ',hl%nsup,'}', hl%name, hl%global_num, &
               xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), zero, hl%spec, ind, hl%nsup, hl%radi
       else if  (which == 'l') then
          write(unit,'(I8,6X,A9,2X,A,I8,A2,I3,A,6X,A2,1X,I8,1X,4F12.4,3X,I3,2I8,1X,F7.3)')     &
               part,'{l\delta}','{',hl%global_num,'\ ',hl%nsup,'}', hl%name, hl%global_num, &
               xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), zero, hl%spec, ind, hl%nsup, hl%radi
       else if  (which == 'j') then
          write(unit,'(I8,6X,A9,2X,A,I8,A2,I3,A,6X,A2,1X,I8,1X,4F12.4,3X,I3,2I8,1X,F7.3)')     &
               part,'{j\beta}','{',hl%global_num,'\ ',hl%nsup,'}', hl%name, hl%global_num, &
               xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), zero, hl%spec, ind, hl%nsup, hl%radi
       end if

       
    end if    
    call stop_timer(tmr_std_exx_fetch,.true.)

    return
  end subroutine get_halodat

!!   2016/07/29 18:30 nakata
!!    Renamed supports_on_atom -> blips_on_atom
!!   2016/12/29 18:30 nakata
!!    Removed no longer used blips_on_atom and flag_one_to_one
  subroutine get_iprimdat(ia,hl,ind,iprim,part,verbose,unit)
    
    use units,          only: BohrToAng
    use numbers,        only: very_small, zero, three
    use primary_module, only: bundle
    use species_module, only: species_label,nsf_species
    use atomic_density, only: atomic_density_table
    !
    use matrix_data, only: mat, Hrange, Srange
    use exx_types, only: prim_atomic_data, neigh_atomic_data
    use exx_types, only: tmr_std_exx_fetch
    use group_module,only: parts 
    use GenComms,    only: inode
    !
    implicit none
    !
    type(prim_atomic_data), intent(inout) :: ia
    type(neigh_atomic_data),intent(inout) :: hl
    integer, intent(in)                   :: ind
    integer, intent(in)                   :: iprim
    integer, intent(in)                   :: part
    logical, intent(in)                   :: verbose
    !
    real(double)      :: xyz_Ang(3)      
    character(len=20) :: filename
    integer, optional :: unit
    !
    call start_timer(tmr_std_exx_fetch)
    ! Get u(v) data
    ia%pr  = ind
    ia%num = bundle%nm_nodbeg(part) + ia%pr - 1 
    ia%ip  = bundle%ig_prim(iprim)             

    ia%spec = bundle%species(ia%num)
    ia%name = species_label(ia%spec)             
    ia%radi = atomic_density_table(ia%spec)%cutoff
    ia%labcell = ia%ip 
    !if (flag_one_to_one) then
    !   ia%nsup = blips_on_atom(ia%spec)%nsuppfuncs
    !else
    !   ia%nsup = blips_on_atom(ia%ip)%nsuppfuncs
    !end if
    ia%nsup = nsf_species(ia%spec)

    ! Calculate R_u    
    ia%xyz_ip(1) = bundle%xprim(ia%num)
    ia%xyz_ip(2) = bundle%yprim(ia%num)
    ia%xyz_ip(3) = bundle%zprim(ia%num)

    ! Calculate R_iu
    ia%xyz(1) = - ia%xyz_ip(1) + hl%xyz_hl(1)
    ia%xyz(2) = - ia%xyz_ip(2) + hl%xyz_hl(2)
    ia%xyz(3) = - ia%xyz_ip(3) + hl%xyz_hl(3)
    ia%r      = sqrt(dot_product(ia%xyz,ia%xyz))
    !
    ! Calculate D_iu
    ia%d      = sqrt(three)*(ia%radi + hl%radi)
    !    
    xyz_Ang(1) = ia%xyz_ip(1)
    xyz_Ang(2) = ia%xyz_ip(2)
    xyz_Ang(3) = ia%xyz_ip(3)
    !
    if ( exx_debug ) then
       write(unit,'(I8,2X,A9,6X,A,I8,A2,I3,A,6X,A2,1X,I8,1X,4F12.4,3X,I3,2I8,1X,F7.3,F12.8)')           &
            parts%ngnode(parts%inode_beg(inode)+part-1),'{i\alpha}','{',ia%labcell,'\ ',ia%nsup,'}', &
            ia%name, ia%labcell, xyz_Ang(1), xyz_Ang(2), xyz_Ang(3), zero, ia%spec, ind, ia%nsup, zero
    end if
    call stop_timer(tmr_std_exx_fetch,.true.)

    return
  end subroutine get_iprimdat
  !!
  !!***
!!$  subroutine initialise_exx( scheme )
!!$    !
!!$    use timer_stdclocks_module, only: tmr_std_exx
!!$    use global_module, only: iprint_exx
!!$    use timer_module,  only: start_timer, stop_timer    
!!$    use species_module,only: nsf_species
!!$    use mult_module,   only: S_X_SX, mult 
!!$    
!!$    use exx_types,     only: extent, exx_psolver, exx_pscheme, r_int, eris
!!$    use exx_types,     only: ewald_alpha, ewald_charge, ewald_rho, ewald_pot
!!$    use exx_types,     only: p_omega, p_ngauss, p_gauss, w_gauss, pulay_radius
!!$    use exx_poisson,   only: exx_scal_rho_3d, exx_ewald_rho, exx_ewald_pot
!!$    
!!$    use exx_types,     only: tmr_std_exx_setup, unit_eri_debug, unit_exx_debug
!!$    use exx_types,     only: unit_memory_write, unit_timers_write
!!$    use exx_memory,    only: exx_mem_alloc
!!$    
!!$    integer, intent(in) :: scheme
!!$            
!!$    character(len=20) :: filename1, filename2, filename3, filename4, filename5, filename6
!!$    integer :: maxsuppfuncs, stat
!!$
!!$    maxsuppfuncs = maxval(nsf_species)
!!$    !
!!$    call start_timer(tmr_std_exx)   
!!$    !
!!$    call start_timer(tmr_std_exx_setup)    
!!$    !
!!$    if ( iprint_exx > 3 ) then
!!$       
!!$       call io_assign(unit_timers_write)
!!$       call get_file_name('exx_timers',numprocs,inode,filename3)
!!$       open(unit_timers_write,file=filename3)
!!$       !
!!$       !
!!$       call io_assign(unit_memory_write)
!!$       call get_file_name('exx_memory',numprocs,inode,filename4)
!!$       open(unit_memory_write,file=filename4)
!!$       
!!$    end if
!!$    
!!$    if ( exx_debug ) then
!!$       
!!$       call io_assign(unit_exx_debug)
!!$       call get_file_name('exx_debug',numprocs,inode,filename5)
!!$       open(unit_exx_debug,file=filename5)       
!!$       
!!$       call io_assign(unit_eri_debug)
!!$       call get_file_name('eri_debug',numprocs,inode,filename6)
!!$       open(unit_eri_debug,file=filename6)       
!!$
!!$       !call exx_write_head(unit_exx_debug,inode,bundle%groups_on_node) 
!!$       !
!!$       !call exx_global_write()
!!$       !
!!$    end if
!!$    !
!!$    maxsuppfuncs = maxval(nsf_species)
!!$    !
!!$    if ( scheme > 0 ) then
!!$       call exx_mem_alloc(extent,0,0,'work_3d' ,'alloc')
!!$    end if
!!$    !
!!$    !DRB! This appears to set up the different Poisson solvers
!!$    if (exx_psolver == 'fftw')     then
!!$       call exx_mem_alloc(extent,0,0,'fftw_3d','alloc')  
!!$       call exx_mem_alloc(extent,0,0,'reckernel_3d','alloc')
!!$
!!$       poisson_fftw: select case(exx_pscheme)          
!!$
!!$       case('default')
!!$          call exx_scal_rho_3d(inode,extent,r_int,exx_pscheme,pulay_radius, &
!!$               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)
!!$
!!$       case('ewald')
!!$          call exx_mem_alloc(extent,0,0,'ewald_3d','alloc')
!!$          call exx_ewald_rho(ewald_rho,extent,ewald_alpha,r_int)
!!$          call exx_ewald_pot(ewald_pot,extent,ewald_alpha,r_int)          
!!$          call exx_scal_rho_3d(inode,extent,r_int,exx_pscheme,pulay_radius, &
!!$               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)
!!$
!!$       case('pulay')
!!$          call exx_scal_rho_3d(inode,extent,r_int,exx_pscheme,pulay_radius, &
!!$               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)
!!$
!!$       case('yukawa')
!!$          call exx_scal_rho_3d(inode,extent,r_int,exx_pscheme,pulay_radius, &
!!$               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)
!!$
!!$       case('gauss')
!!$          call cq_abort('EXX: Gaussian representation if 1/r for solving &
!!$               &the Poisson equation &
!!$               &is currently under testing...')
!!$          !call createBeylkin(p_gauss,w_gauss,r_int)     
!!$          !call exx_scal_rho_3d(inode,extent,r_int,exx_pscheme,pulay_radius, &
!!$          !     p_omega,p_ngauss,p_gauss,w_gauss)
!!$
!!$       case default
!!$          call exx_scal_rho_3d(inode,extent,r_int,exx_pscheme,pulay_radius, &
!!$               p_omega,p_ngauss,p_gauss,w_gauss,reckernel_3d)
!!$
!!$       end select poisson_fftw
!!$
!!$    else if (exx_psolver == 'isf') then
!!$       call cq_abort('EXX: ISF Poisson solver is not available yet ; &
!!$            &under optimisation...')
!!$       !call exx_mem_alloc(extent,0,0,'isf_rho','alloc')
!!$       !call createKernel('F',ngrid,ngrid,ngrid,grid_spacing,grid_spacing,grid_spacing,isf_order,&
!!$       !     0,1,kernel)       
!!$
!!$    end if
!!$    call stop_timer(tmr_std_exx_setup,.true.)    
!!$    
!!$    if ( scheme > 0 ) then
!!$       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_i','alloc')       
!!$       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_j','alloc')
!!$       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_k','alloc')
!!$       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_l','alloc')
!!$       !
!!$       if ( scheme == 1 ) then
!!$          call exx_mem_alloc(extent,maxsuppfuncs,0,'Phy_k','alloc')
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'Ome_kj','alloc')       
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_kj','alloc')       
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_kj','alloc')
!!$          !
!!$       else if ( scheme == 2 ) then       
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_ki','alloc')       
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_lj','alloc')
!!$          !
!!$       else if ( scheme == 3 ) then       
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_ki','alloc')       
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_lj','alloc')
!!$          !
!!$          !if ( niter == 1 ) then
!!$          allocate( eris( mult(S_X_SX)%ahalo%np_in_halo ), STAT=stat)
!!$          if(stat/=0) call cq_abort('Error allocating memory to eris/exx !',stat)
!!$          !end if
!!$          !for now allocated here ; should not ; better to have this before
!!$          ! SCF calculation ; note that here it is never deallocated!
!!$          !eris_size = int( sqrt(dble(size( mat_p(matX(  exxspin  ))%matrix ))) )
!!$          !print*, 'inode', inode, inode, size( mat_p(matX(  exxspin  ))%matrix), eris_size, &
!!$          !     & mat_p(matX(  exxspin  ))%length, mat_p(matX(  exxspin  ))%sf1_type, mat_p(matX(  exxspin  ))%sf2_type
!!$          !call exx_mem_alloc(eris_size**4,0,0,'eris','alloc')
!!$          !
!!$          !
!!$       end if
!!$       !
!!$       !
!!$    end if
!!$    
!!$  end subroutine initialise_exx
!!$  !
!!$  !
!!$  subroutine finalise_exx( scheme )
!!$    !
!!$    use global_module, only: iprint_exx, area_exx
!!$    use timer_stdclocks_module, only: tmr_std_exx
!!$    use timer_module,  only: start_timer, stop_timer, print_timer
!!$    use species_module,only: nsf_species
!!$    use memory_module, only: write_mem_use
!!$    
!!$    use exx_types,     only: extent, exx_psolver ,kernel, exx_pscheme
!!$    use exx_types,     only: tmr_std_exx_setup, unit_eri_debug, unit_exx_debug
!!$    use exx_types,     only: tmr_std_exx_evalpao, tmr_std_exx_poisson, tmr_std_exx_matmult
!!$    use exx_types,     only: tmr_std_exx_accumul, tmr_std_exx_allocat, tmr_std_exx_dealloc
!!$    use exx_types,     only: tmr_std_exx_fetch, tmr_std_exx_barrier, tmr_std_exx_kernel
!!$    use exx_types,     only: tmr_std_exx_kernel
!!$    use exx_types,     only: exx_total_time
!!$    
!!$    use exx_types,     only: unit_memory_write, unit_eri_debug, unit_exx_debug, unit_timers_write
!!$    use exx_memory,    only: exx_mem_alloc
!!$    !
!!$    integer, intent(in) :: scheme
!!$    integer :: maxsuppfuncs
!!$
!!$    maxsuppfuncs = maxval(nsf_species)
!!$    
!!$    if ( scheme > 0 ) then
!!$       call exx_mem_alloc(extent,0,0,'work_3d' ,'dealloc')
!!$       !
!!$       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_i','dealloc')       
!!$       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_j','dealloc')
!!$       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_k','dealloc')    
!!$       call exx_mem_alloc(extent,maxsuppfuncs,0,'phi_l','dealloc')
!!$       !
!!$       if ( scheme == 1 ) then
!!$          call exx_mem_alloc(extent,maxsuppfuncs,0,'Phy_k','dealloc')
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'Ome_kj','dealloc')   
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_kj','dealloc')       
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_kj','dealloc')
!!$          !
!!$       else if( scheme == 2 ) then
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_ki','dealloc')       
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_lj','dealloc')
!!$          !
!!$       else if( scheme == 3 ) then
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'rho_ki','dealloc')       
!!$          call exx_mem_alloc(extent,maxsuppfuncs,maxsuppfuncs,'vhf_lj','dealloc')
!!$          !
!!$       end if
!!$       !
!!$    end if
!!$    !
!!$    !
!!$    if (exx_psolver == 'fftw')  then
!!$       call exx_mem_alloc(extent,0,0,'fftw_3d',     'dealloc')  
!!$       call exx_mem_alloc(extent,0,0,'reckernel_3d','dealloc')
!!$       select case(exx_pscheme)
!!$       case('default')         
!!$       case('gauss')          
!!$       case('pulay')
!!$       case('ewald')
!!$          call exx_mem_alloc(extent,0,0,'ewald_3d','dealloc')
!!$       case default
!!$       end select
!!$
!!$    else if (exx_psolver == 'isf') then
!!$       call exx_mem_alloc(extent,0,0,'isf_rho','dealloc')
!!$       deallocate(kernel)
!!$    end if
!!$    !
!!$    if (inode == ionode) call write_mem_use(unit_memory_write,area_exx)
!!$    !
!!$    call start_timer(tmr_std_exx_barrier)
!!$    call my_barrier()
!!$    call stop_timer(tmr_std_exx_barrier,.true.)
!!$    !
!!$    call stop_timer(tmr_std_exx,.true.)
!!$    !    
!!$    ! Timers and elapse time
!!$    exx_total_time = &
!!$         tmr_std_exx_evalpao%t_tot + &
!!$         tmr_std_exx_poisson%t_tot + &
!!$         tmr_std_exx_matmult%t_tot + &
!!$         tmr_std_exx_accumul%t_tot + &
!!$         tmr_std_exx_allocat%t_tot + &
!!$         tmr_std_exx_dealloc%t_tot + &
!!$         tmr_std_exx_setup%t_tot   + &
!!$         tmr_std_exx_write%t_tot   + &
!!$         tmr_std_exx_fetch%t_tot   + &
!!$         tmr_std_exx_barrier%t_tot
!!$
!!$    if ( iprint_exx > 3 ) then
!!$       
!!$       call print_timer(tmr_std_exx_setup,  "exx_setup   time:", unit_timers_write)    
!!$       call print_timer(tmr_std_exx_write,  "exx_write   time:", unit_timers_write)    
!!$       !call print_timer(tmr_std_exx_kernel, "exx_kernel  time:", unit_timers_write)    
!!$       call print_timer(tmr_std_exx_fetch,  "exx_fetch   time:", unit_timers_write)    
!!$       call print_timer(tmr_std_exx_evalpao,"exx_evalpao time:", unit_timers_write)    
!!$       call print_timer(tmr_std_exx_poisson,"exx_poisson time:", unit_timers_write)
!!$       call print_timer(tmr_std_exx_matmult,"exx_matmult time:", unit_timers_write)    
!!$       call print_timer(tmr_std_exx_accumul,"exx_accumul time:", unit_timers_write)    
!!$       call print_timer(tmr_std_exx_allocat,"exx_allocat time:", unit_timers_write)    
!!$       call print_timer(tmr_std_exx_dealloc,"exx_dealloc time:", unit_timers_write)    
!!$       call print_timer(tmr_std_exx_barrier,"exx_barrier time:", unit_timers_write)    
!!$       write(unit_timers_write,*)
!!$       call print_timer(tmr_std_exx_kernel, "exx_kernel  time:", unit_timers_write)    
!!$       call print_timer(tmr_std_exx,        "exx_total   time:", unit_timers_write)    
!!$
!!$       write(unit=unit_timers_write,fmt='("Timing: Proc ",i6,": Time spent in ", a50, " = ", &
!!$            &f12.5," s")') inode, 'get_X_matrix',  tmr_std_exx%t_tot
!!$
!!$       write(unit=unit_timers_write,fmt='("Timing: Proc ",i6,": Time spent in ", a50, " = ", &
!!$            &f12.5," s")') inode, 'timer calls',  tmr_std_exx%t_tot-exx_total_time 
!!$       !call io_close(unit_matrix_write)
!!$
!!$       !call io_close(unit_output_write)
!!$       !call io_close(unit_screen_write)
!!$       !
!!$       !call io_close(unit_eri_debug)
!!$       !call io_close(unit_exx_debug)
!!$       call io_close(unit_memory_write)
!!$       call io_close(unit_timers_write)
!!$       !
!!$       !
!!$    end if
!!$
!!$    if ( exx_debug ) then
!!$       call io_close(unit_eri_debug)
!!$       call io_close(unit_exx_debug)       
!!$    end if
!!$    
!!$  end subroutine finalise_exx
  !
  end module exx_module
