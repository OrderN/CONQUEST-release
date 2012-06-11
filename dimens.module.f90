! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module dimens
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/dimens *
!!  NAME
!!   dimens
!!  PURPOSE
!!   Stores various useful dimensions and a subroutine
!!  USES
!!   GenComms, datatypes, global_module, matrix_data, numbers, species_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   08/05/01 DRB
!!  MODIFICATION HISTORY
!!   10/05/01 drb
!!     Added blip_width, four_on_blip_width, fobw2, fobw3
!!   11/06/2001 dave
!!    Added RCS Id and Log tags and altered argument list 
!!    to set_dimensions
!!   20/06/2001 dave
!!    Used cq_abort instead of stop
!!   18/03/2002 dave
!!    Added header and static object for object file id
!!   17:28, 2003/06/10 tm
!!    Flag to allow maximum numbers not to be checked
!!   2006/10/19 16:49 dave
!!    Added automatic grid density finding, and GridCutoff variable
!!   2008/02/04 17:14 dave
!!    Changes for output to file not stdout
!!   12:15, 14/02/2008 drb 
!!    Added variables for buffer around primary and covering sets
!!   2011/09/16 11:04 dave
!!    Changed variable atomicrad to atomicnum
!!  SOURCE
module dimens

  use datatypes
  use global_module, ONLY: io_lun

  implicit none
  save

  ! RCS tag for object file identification 
  character(len=80), private :: RCSid = "$Id$"

  real(double) :: r_super_x, r_super_y, r_super_z, volume
  real(double) :: r_super_x_squared, r_super_y_squared, r_super_z_squared
  real(double) :: r_s, r_h, r_c, r_core_squared, r_dft_d2
  real(double) :: grid_point_volume, one_over_grid_point_volume
  real(double) :: support_grid_volume
  real(double) :: GridCutoff, min_blip_sp
  real(double) :: AtomMove_buffer ! Buffer around primary sets and covering sets
  logical :: flag_buffer_old
  !real(double) :: support_grid_spacing, support_grid_volume
  !real(double) ::  blip_width, four_on_blip_width, fobw2, fobw3

  real(double), allocatable, dimension(:) :: RadiusSupport, &
                                             NonLocalFactor, InvSRange
  integer, allocatable, dimension(:) :: atomicnum

  integer :: n_grid_x, n_grid_y, n_grid_z, n_my_grid_points

  real(double) :: x_grid, y_grid, z_grid

  integer, parameter :: ngrids = 105
  !All grid spacings: multiples of 3, 4, 5
  integer, dimension(ngrids), parameter :: grid_spacings = (/   &
       4,    8,    12,   16,   20,   24,   32,   36,   40,   48,  &
       60,   64,   72,   80,   96,   100,  108,  120,  128,  144, &
       160,  180,  192,  200,  216,  240,  256,  288,  300,  320, &
       324,  360,  384,  400,  432,  480,  500,  512,  540,  576, &
       600,  640,  648,  720,  768,  800,  864,  900,  960,  972, &
       1000, 1024, 1080, 1152, 1200, 1280, 1296, 1440, 1500,&
       1536, 1600, 1620, 1728, 1800, 1920, 1944, 2000, 2048, 2160,&
       2304, 2400, 2500, 2560, 2592, 2700, 2880, 3000, 3072, 3200,&
       3240, 3456, 3600, 3840, 3888, 4000, 4096, 4320, 4500,&
       4608, 4800, 5000, 5120, 5184, 5400, 5760, 6000, 6144, 6400,&
       6480, 6912, 7200, 7500, 7680, 8000, 8192 /)

!!***

contains

! -----------------------------------------------------------
! Subroutine set_dimensions
! -----------------------------------------------------------

!!****f* dimens/set_dimensions *
!!
!!  NAME 
!!   set_dimensions
!!  USAGE
!! 
!!  PURPOSE
!!   Sets up various parameters and dimensions 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   EHH, CMG, D.R.Bowler
!!  CREATION DATE
!!   5/5/95
!!  MODIFICATION HISTORY
!!   3/10/95 CMG 
!!   23/1/96 EHH/CMG Added dimension setting
!!   08/05/01 DRB changed format to ROBODoc and simplified output
!!   10/05/01 DRB added matrix names
!!   11/06/2001 dave
!!    Changed to pass data from pseudopotential_data
!!   20/06/2001 dave
!!    Used cq_abort instead of stop
!!   15:57, 04/02/2003 drb 
!!    Changed scan over core radii to be for ALL cases: we need this
!!    for local pseudos too
!!   12:17, 2004/06/09 dave
!!    Added scan for separate InvSrange to allow larger radii.
!!   2007/04/17 09:42 dave
!!    Removed rcut_dens
!!   2008/03/03 18:44 dave
!!    Changed float to real
!!   2008/03/12 06:13 dave
!!    Added loop to find longest range matrix
!!   2009/07/08 16:40 dave
!!    Zero r_t
!!   2011/09/29 M. Arita
!!    Set the cutoff for DFT-D2
!!  SOURCE
!!
  subroutine set_dimensions(inode, ionode,HNL_fac,non_local, n_species, non_local_species, core_radius)

    use datatypes
    use numbers
    use global_module, ONLY: iprint_init, flag_basis_set, blips, runtype, flag_dft_d2
    use matrix_data
    use GenComms, ONLY: cq_abort
    use pseudopotential_common, ONLY: pseudo_type, OLDPS, SIESTA, ABINIT
    use block_module, ONLY: in_block_x, in_block_y, in_block_z, n_pts_in_block, n_blocks
    use input_module, ONLY: leqi

    implicit none

    ! Passed variables
    logical :: non_local, non_local_species(:)
    integer :: inode, ionode, n_species

    real(double) :: HNL_fac, core_radius(:)

    ! Local variables
    integer :: max_extent, n, stat
    real(double) :: r_core, r_t, rcutmax

    !n_my_grid_points = n_pts_in_block * n_blocks    
    ! decide if the flag non_local needs to be set to true
    if(pseudo_type == OLDPS) then
       non_local = .false.
       do n=1, n_species
          if(non_local_species(n)) non_local = .true.
       end do
    elseif(pseudo_type == SIESTA.OR.pseudo_type==ABINIT) then
       non_local = .true.
    else
       call cq_abort('ERROR : pseudo_type in set_dimension',pseudo_type)
    endif

    ! we need to decide which is the largest core radius
    r_core = zero
    r_h = zero
    r_t = zero
    !    if(non_local) then
    do n=1, n_species
       r_core = max(r_core,core_radius(n))
       r_h = max(r_h,RadiusSupport(n))
       r_t = max(r_t,InvSRange(n))
    end do
    !    else
    !       do n=1, n_species
    !          r_core = max(r_core,core_radius(n))
    !       end do
    !    end if
    if(non_local.and.(inode==ionode).and.iprint_init>0) then
       write(io_lun,2) r_core
2      format(8x,'This calculation includes non-local pseudopotentials,'/&
            9x,'with a maximum core radius of ',f15.8)
    end if
    if (r_core>r_h) then
       call cq_abort('set_dimens: r_core > r_support')
    end if

    ! then, obtain the blip width
    !four_on_blip_width = four / blip_width
    !fobw2 = four_on_blip_width*four_on_blip_width
    !fobw3 = four_on_blip_width*fobw2

    ! Set range of S matrix
    r_s = r_h
    if(two*r_s>r_c) then
       if(inode==ionode) write(io_lun,fmt='(8x,"WARNING ! S range greater than L !")')
       !r_s = r_c
       !r_h = r_c 
    endif
    if(.NOT.leqi(runtype,'static')) then
       if(flag_buffer_old) then
          r_s = 1.1_double * r_s
          r_c = 1.1_double * r_c
          r_h = 1.1_double * r_h
          r_core = 1.1_double * r_core
          if (flag_dft_d2) r_dft_d2 = 1.1_double * r_dft_d2 ! for DFT-D2
       else
          r_s = AtomMove_buffer +  r_s
          r_c = AtomMove_buffer +  r_c
          r_h = AtomMove_buffer +  r_h
          r_core = AtomMove_buffer +  r_core
          if (flag_dft_d2) r_dft_d2 = AtomMove_buffer + r_dft_d2 ! for DFT-D2
       end if
     endif

    ! Set other ranges
    r_core_squared = r_core * r_core
    ! Set matrix ranges (from matrix_data)
    rcut(Srange)   = (two*r_s)
    rcut(Trange)   = (two*r_t)
    rcut(Lrange)   = (r_c)
    rcut(Hrange)   = (two*(r_h+HNL_fac*r_core))
    rcut(SPrange)  = (r_core + r_h)
    rcut(LSrange)  = (rcut(Lrange) + rcut(Srange))
    rcut(LHrange)  = (rcut(Lrange) + rcut(Hrange))
    rcut(LSLrange) = (two*rcut(Lrange) + rcut(Srange))
    rcut(SLSrange) = (two*rcut(Srange) + rcut(Lrange))
    ! New 12:16, 2004/06/09 dave
    !rcut(Trange)   = (rcut(Srange))
    rcut(TSrange)  = (rcut(Trange)+rcut(Srange))
    rcut(THrange)  = (rcut(Trange)+rcut(Hrange))
    rcut(TLrange)  = (rcut(Trange)+rcut(Lrange))
    rcut(PSrange)  = rcut(SPrange)
    rcut(LTrrange) = rcut(Lrange)
    rcut(SLrange) = rcut(LSrange)
    rcut(TTrrange) = rcut(Trange)
    rcut(dSrange) = rcut(Srange)
    rcut(dHrange) = rcut(Hrange)
    rcut(PAOPrange) = rcut(SPrange)
    rcut(HLrange) = rcut(LHrange)
    rcutmax = zero
    do n=1,mx_matrices
       if(rcut(n)>rcutmax) then
          max_range = n
          rcutmax = rcut(n)
       end if
    enddo
    ! Seemingly trivial, but may be quite useful - matrix names
    mat_name(Srange)   = "S  "
    mat_name(Lrange)   = "L  "
    mat_name(Hrange)   = "H  "
    mat_name(SPrange)  = "SP "
    mat_name(LSrange)  = "LS "
    mat_name(LHrange)  = "LH "
    mat_name(LSLrange) = "LSL"
    mat_name(SLSrange) = "SLS"
    mat_name(Trange)   = "T  "
    mat_name(TSrange)  = "TS "
    mat_name(THrange)  = "TH "
    mat_name(TLrange)  = "TL "
    if(inode==ionode.AND.iprint_init>1) then
       do n=1,12!mx_matrices
          write(io_lun,1) mat_name(n),rcut(n)
       enddo
    endif
1   format(8x,'Matrix ',a3,' has range ',f15.8)
    ! Various useful parameters
    r_super_x_squared = r_super_x * r_super_x
    r_super_y_squared = r_super_y * r_super_y
    r_super_z_squared = r_super_z * r_super_z
    volume = r_super_x * r_super_y * r_super_z
    grid_point_volume = volume/(n_grid_x*n_grid_y*n_grid_z)
    one_over_grid_point_volume = one / grid_point_volume
    !support_grid_volume = support_grid_spacing**3
    x_grid = one / real( n_grid_x, double)
    y_grid = one / real( n_grid_y, double)
    z_grid = one / real( n_grid_z, double)
    ! Grid points in a block
    n_pts_in_block = in_block_x * in_block_y * in_block_z
    return
  end subroutine set_dimensions
!!***

  subroutine find_grid

    use datatypes
    use numbers, ONLY: pi, very_small,two, half
    use global_module, ONLY: iprint_init, flag_basis_set, blips
    use GenComms, ONLY: inode,ionode,cq_abort
    use species_module, ONLY: n_species
    
    implicit none
    
    ! Local variables
    real(double) :: sqKE, blipKE, blip_sp
    integer :: mingrid, i

    if(n_grid_x>0.AND.n_grid_y>0.AND.n_grid_z>0) then ! We dont need to find grids
       if(iprint_init>1.AND.inode==ionode) &
            write(io_lun,fmt='(2x,"User specified grid dimensions: ",3i5)') n_grid_x,n_grid_y,n_grid_z
       return
    else
       if(GridCutoff<very_small) call cq_abort("Grid cutoff too small: ",GridCutoff)
       if(flag_basis_set==blips) then ! Test for fine enough grid
          blip_sp = half*min_blip_sp
          blipKE = half*(pi/blip_sp)*(pi/blip_sp)
          if(GridCutoff<blipKE) then
             if(iprint_init>0.AND.inode==ionode) &
                  write(io_lun,fmt='(2x,"Warning ! Adjusted grid cutoff to ",f8.3)') blipKE
             GridCutoff = blipKE
          end if
       end if
       sqKE = sqrt(two*GridCutoff) ! Use KE = 0.5*G*G
       ! For now, we assume orthorhombic cells
       ! x first
       mingrid = int(sqKE*r_super_x/pi)
       do i=1,ngrids
          if(grid_spacings(i)>mingrid) then
             n_grid_x = grid_spacings(i)
             exit
          end if
       end do
       ! y
       mingrid = int(sqKE*r_super_y/pi)
       do i=1,ngrids
          if(grid_spacings(i)>mingrid) then
             n_grid_y = grid_spacings(i)
             exit
          end if
       end do
       ! z
       mingrid = int(sqKE*r_super_z/pi)
       do i=1,ngrids
          if(grid_spacings(i)>mingrid) then
             n_grid_z = grid_spacings(i)
             exit
          end if
       end do
    endif

  end subroutine find_grid

end module dimens
