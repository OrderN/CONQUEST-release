! -*- mode: F90; mode: font-lock -*-
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
!!   2013/07/01 M.Arita
!!    Commented out some tricks used in MD & structure optimisation
!!   2016/09/16 16:00 nakata
!!    Added variables for PAO-based matrices and multi-site SFs
!!  SOURCE
module dimens

  use datatypes
  use global_module, only: io_lun
  use timer_module,  only: start_timer,     stop_timer, cq_timer
  use timer_module,  only: start_backtrace, stop_backtrace

  implicit none
  save

  real(double) :: r_super_x, r_super_y, r_super_z, volume
  real(double) :: r_super_x_squared, r_super_y_squared, r_super_z_squared
  real(double) :: r_s, r_h, r_c, r_nl, r_core_squared, r_dft_d2, r_exx
  real(double) :: r_s_atomf, r_h_atomf, r_MS, r_LD
  real(double) :: grid_point_volume, one_over_grid_point_volume
  real(double) :: support_grid_volume
  real(double) :: GridCutoff, min_blip_sp
  real(double) :: AtomMove_buffer ! Buffer around primary sets and covering sets
  logical :: flag_buffer_old
  !real(double) :: support_grid_spacing, support_grid_volume
  !real(double) ::  blip_width, four_on_blip_width, fobw2, fobw3

  real(double), allocatable, dimension(:) :: RadiusSupport, RadiusAtomf, &
                                             RadiusMS, RadiusLD,         &
                                             NonLocalFactor, InvSRange
  integer, allocatable, dimension(:) :: atomicnum

  integer :: n_grid_x, n_grid_y, n_grid_z, n_my_grid_points

  real(double) :: x_grid, y_grid, z_grid

  integer, parameter :: ngrids = 105
  !All grid spacings: multiples of 3, 4, 5
  integer, dimension(ngrids), parameter :: grid_spacings = (/     &
       4,    8,    12,   16,   20,   24,   32,   36,   40,   48,  &
       60,   64,   72,   80,   96,   100,  108,  120,  128,  144, &
       160,  180,  192,  200,  216,  240,  256,  288,  300,  320, &
       324,  360,  384,  400,  432,  480,  500,  512,  540,  576, &
       600,  640,  648,  720,  768,  800,  864,  900,  960,  972, &
       1000, 1024, 1080, 1152, 1200, 1280, 1296, 1440, 1500,      &
       1536, 1600, 1620, 1728, 1800, 1920, 1944, 2000, 2048, 2160,&
       2304, 2400, 2500, 2560, 2592, 2700, 2880, 3000, 3072, 3200,&
       3240, 3456, 3600, 3840, 3888, 4000, 4096, 4320, 4500,      &
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
!!   2013/07/03 M.Arita
!!    Introduced flag_MDold to choose the way of member updates.
!!    Now the code does not expand cutoff ranges automatically by default.
!!   2014/01/17 lat 
!!    Added r_nl: fix unconsistency between HNL_fac/NonLocalFactor
!!   2014/01/17 lat
!!    Added r_exx and assigned to rcut(Xrange)
!!   2016/09/16 16:00 nakata
!!    Set ranges of atomf-based matrices
!!   2016/11/02 dave
!!    Added flag so that we only warn Lrange<Srange if O(N) is used
!!   2017/03/08 15:00 nakata
!!    Removed dSrange, dHrange and PAOPrange which are no longer used
!!   2017/02/23 dave
!!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
!!   2017/12/05 11:00 dave
!!    Adding NA projector ranges
!!   2018/05/15 12:21 dave
!!    Bug fix: mataNA range was wrong with MSSF (didn't correctly use atom function ranges)
!!   2019/11/18 tsuyoshi
!!    flag_MDold was removed.
!!   2019/12/02 15:17 dave 
!!    Added checks to round RadiusAtomf and RadiusSupport to safe value (including grid points)
!!   2034/09/16 14:18 lionel 
!!    Check consistency of Xrange wrt r_exx read from input 
!!  SOURCE
!!
  subroutine set_dimensions(inode, ionode,HNL_fac,non_local, n_species, non_local_species, core_radius)

    use datatypes
    use numbers
    use matrix_data
    use GenComms,      only: cq_abort, cq_warn
    use global_module, only: iprint_init, atomf, sf, flag_Multisite, flag_exx, flag_diagonalisation
    use block_module,  only: in_block_x, in_block_y, in_block_z, n_pts_in_block
    use pseudopotential_common, only: pseudo_type, OLDPS, SIESTA, ABINIT, flag_neutral_atom_projector

    implicit none

    ! Passed variables
    logical :: non_local, non_local_species(:)
    integer :: inode, ionode, n_species

    real(double) :: HNL_fac, core_radius(:)

    ! Local variables
    character(len=80) :: sub_name = "set_dimensions"
    type(cq_timer) :: backtrace_timer
    integer        :: n, mx_matrices_tmp
    real(double)   :: r_core, r_t, rcutmax, max_grid

!****lat<$
    call start_backtrace(t=backtrace_timer,who='set_dimensions',where=9,level=2)
!****lat>$

    ! Various grid parameters
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
    ! Calculate largest grid spacing in any direction (local to routine for now)
    max_grid = max(r_super_x*x_grid,r_super_y*y_grid,r_super_z*z_grid)
    ! Grid points in a block
    n_pts_in_block = in_block_x * in_block_y * in_block_z

    ! Set indices for atomf-based-matrix ranges
    if (atomf==sf) then
       aSa_range = Srange
       aHa_range = Hrange
       if(flag_neutral_atom_projector) then
          aNArange = 20
          NAarange = 21
          mx_matrices_tmp = 21
       else
          mx_matrices_tmp = 19
       end if
    else
       aSa_range       = 20
       aHa_range       = 21
       STr_range       = 22
       HTr_range       = 23
       aSs_range       = 24
       aHs_range       = 25
       sSa_range       = 26
       sHa_range       = 27
       SFcoeff_range   = 28
       SFcoeffTr_range = 29
       LD_range        = 30
       if(flag_neutral_atom_projector) then
          aNArange        = 31
          NAarange        = 32
          mx_matrices_tmp = mx_matrices ! = 30
       else
          mx_matrices_tmp = 30
       end if
    endif

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
    r_core    = zero
    r_h       = zero
    r_t       = zero
    r_nl      = zero
    r_h_atomf = zero
    r_MS      = zero
    r_LD      = zero
    do n=1, n_species
       ! Round the atom function radius to the nearest multiple of the largest grid spacing
       ! (add one to account for displacement relative to grid point)
       RadiusAtomf(n) = max_grid*(floor(RadiusAtomf(n)/max_grid)+two)
       RadiusSupport(n) = max(RadiusSupport(n), RadiusAtomf(n))
       r_core    = max(r_core,core_radius(n))
       r_h       = max(r_h,RadiusSupport(n))
       r_t       = max(r_t,InvSRange(n))
       r_nl      = max(r_nl,NonLocalFactor(n)*core_radius(n))
       r_h_atomf = max(r_h_atomf,RadiusAtomf(n))
       r_MS      = max(r_MS,RadiusMS(n))
       r_LD      = max(r_LD,RadiusLD(n))
    end do
    if(non_local.and.(inode==ionode).and.iprint_init>0) then
       write(io_lun,2) r_core
2      format(4x,'This calculation includes non-local pseudopotentials,'/&
            4x,'with a maximum core radius of ',f15.8)
    end if
    if (r_core-r_h>1e-4_double) then
       call cq_abort('set_dimens: r_core > r_support',r_core,r_h)
    end if

    ! then, obtain the blip width
    !four_on_blip_width = four / blip_width
    !fobw2 = four_on_blip_width*four_on_blip_width
    !fobw3 = four_on_blip_width*fobw2

    ! Set range of S matrix
    r_s       = r_h
    r_s_atomf = r_h_atomf
    if(two*r_s>r_c.AND.(.NOT.flag_diagonalisation)) then ! Only warn if using O(N)
       call cq_warn(sub_name, "S range greater than L !",two*r_s,r_c)
    endif
    ! Set ranges for contracted SFs
    if (flag_Multisite) then
       r_s = r_s_atomf + r_MS
       r_h = r_h_atomf + r_MS
       if (r_t.lt.r_s) r_t = r_s
    endif

    ! Set other ranges
    r_core_squared = r_core * r_core
    ! Set matrix ranges (from matrix_data)
    if (HNL_fac <= 0.1_double) then
       rcut(Hrange) = (two*(r_h + r_nl)) 
    else
       rcut(Hrange) = (two*(r_h+HNL_fac*r_core))
    end if

    ! **< lat >** don't touch
    if(flag_exx) then
       if ( r_exx <= zero ) then
          rcut(Xrange)  = rcut(Hrange)
          rcut(SXrange) = rcut(Hrange)
       else
          rcut(Xrange)   = two*r_exx
          rcut(SXrange)  = two*r_exx
       endif
    else
       rcut(Xrange)   = two
       rcut(SXrange)  = two
    end if
    
    ! New 16:25, 2014/08/29  **< lat >**
    !rcut(Hrange)   = (two*(r_h + r_nl))
    !rcut(Hrange)   = (two*(r_h+HNL_fac*r_core))
    rcut(Srange)   = (two*r_s)
    rcut(Trange)   = (two*r_t)
    rcut(Lrange)   = (r_c)
    rcut(APrange)  = (r_core + r_h_atomf)
    rcut(LSrange)  = (rcut(Lrange) + rcut(Srange))
    rcut(LHrange)  = (rcut(Lrange) + rcut(Hrange))
    rcut(LSLrange) = (two*rcut(Lrange) + rcut(Srange))
    rcut(SLSrange) = (two*rcut(Srange) + rcut(Lrange))
    ! New 12:16, 2004/06/09 dave
    !rcut(Trange)   = (rcut(Srange))
    rcut(TSrange)  = (rcut(Trange)+rcut(Srange))
    rcut(THrange)  = (rcut(Trange)+rcut(Hrange))
    rcut(TLrange)  = (rcut(Trange)+rcut(Lrange))
    rcut(PArange)  = rcut(APrange)
    rcut(LTrrange) = rcut(Lrange)
    rcut(SLrange)  = rcut(LSrange)
    rcut(TTrrange) = rcut(Trange)
    rcut(HLrange)  = rcut(LHrange)
    ! for atomf-based matrices
    if (atomf.ne.sf) then
       rcut(aSa_range)       = two*r_s_atomf
       rcut(aHa_range)       = rcut(Hrange) - two*r_MS   ! = two*(r_h_atomf+HNL_fac*r_core)
       rcut(STr_range)       = rcut(Srange)
       rcut(HTr_range)       = rcut(Hrange)
       rcut(aSs_range)       = r_s_atomf + r_s
       rcut(aHs_range)       = rcut(Hrange) - r_MS       ! = (r_h_atomf+HNL_fac*r_core)+(r_h+HNL_fac*r_core)
       rcut(sSa_range)       = rcut(aSs_range)
       rcut(sHa_range)       = rcut(aHs_range)
       rcut(SFcoeff_range)   = r_MS
       rcut(SFcoeffTr_range) = r_MS
       rcut(LD_range)        = r_LD
       if (abs(r_MS)<very_small) then
          rcut(SFcoeff_range)   = 0.001_double
          rcut(SFcoeffTr_range) = 0.001_double
       endif
       if (abs(r_LD)<very_small) rcut(LD_range) = 0.001_double
    endif
    if(flag_neutral_atom_projector) then
       rcut(aNArange)   = r_s_atomf + r_h_atomf
       rcut(NAarange)   = rcut(aNArange)
       mat_name(aNArange) = "aNA"
       mat_name(NAarange) = "NAa"
    end if
    rcutmax = zero
    do n=1,mx_matrices_tmp
       if(rcut(n)>rcutmax) then
          max_range = n
          rcutmax = rcut(n)
       end if
    enddo
    ! Seemingly trivial, but may be quite useful - matrix names
    ! agreed... 2014/08/29 LAT
    mat_name(Srange)   = "S  "
    mat_name(Lrange)   = "L  "
    mat_name(Hrange)   = "H  "
    mat_name(APrange)  = "AP "
    mat_name(LSrange)  = "LS "
    mat_name(LHrange)  = "LH "
    mat_name(LSLrange) = "LSL"
    mat_name(SLSrange) = "SLS"
    mat_name(Trange)   = "T  "
    mat_name(TSrange)  = "TS "
    mat_name(THrange)  = "TH "
    mat_name(TLrange)  = "TL "
    mat_name(PArange)  = "PA "
    mat_name(LTrrange) = "LT "
    mat_name(SLrange)  = "SL "
    mat_name(TTrrange) = "TT "
    mat_name(HLrange)  = "HL "
    mat_name(Xrange)   = "X  "
    mat_name(SXrange)  = "SX "
    if (atomf.ne.sf) then
       mat_name(aSa_range)       = "aSa"
       mat_name(aHa_range)       = "aHa"
       mat_name(STr_range)       = "St"
       mat_name(HTr_range)       = "Ht"
       mat_name(aSs_range)       = "aSs"
       mat_name(aHs_range)       = "aHs"
       mat_name(sSa_range)       = "sSa"
       mat_name(sHa_range)       = "sHa"
       mat_name(SFcoeff_range)   = "MS"
       mat_name(SFcoeffTr_range) = "MSt"
       mat_name(LD_range)        = "LD"
    endif
    if(inode==ionode.AND.iprint_init>1) then
       do n=1,mx_matrices_tmp
          write(io_lun,1) mat_name(n),rcut(n)
       enddo
    endif
1   format(8x,'Matrix ',a3,' has range ',f15.8)
!****lat<$
    call stop_backtrace(t=backtrace_timer,who='set_dimensions')
!****lat>$

    return
  end subroutine set_dimensions
!!***

  subroutine find_grid

    use datatypes
    use numbers,        only: pi, very_small,two, half
    use global_module,  only: iprint_init, flag_basis_set, blips
    use GenComms,       only: inode,ionode,cq_abort, cq_warn
    
    implicit none
    
    ! Local variables
    character(len=80) :: sub_name="find_grid"
    type(cq_timer) :: backtrace_timer
    real(double)   :: sqKE, blipKE, blip_sp
    integer        :: mingrid, i

!****lat<$
    call start_backtrace(t=backtrace_timer,who='find_grid',where=9,level=2)
!****lat>$

    if(n_grid_x>0.AND.n_grid_y>0.AND.n_grid_z>0) then ! We dont need to find grids
       if(iprint_init>1.AND.inode==ionode) &
            write(io_lun,fmt='(8x,"User specified grid dimensions: ",3i5)') n_grid_x,n_grid_y,n_grid_z
       call stop_backtrace(t=backtrace_timer,who='find_grid')
       return
    else
       if(GridCutoff<very_small) call cq_abort("Grid cutoff too small: ",GridCutoff)
       if(flag_basis_set==blips) then ! Test for fine enough grid
          blip_sp = half*min_blip_sp
          blipKE = half*(pi/blip_sp)*(pi/blip_sp)
          if(GridCutoff<blipKE) then
             if(iprint_init>0) call cq_warn(sub_name, "Adjusted grid cutoff to ",blipKE)
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

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='find_grid')
!****lat>$

    return
  end subroutine find_grid

end module dimens
