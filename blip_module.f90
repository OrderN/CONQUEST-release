! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module blip
! ------------------------------------------------------------------------------
! Code area 11: Basis functions
! ------------------------------------------------------------------------------

!!****h* Conquest/blip *
!!  NAME
!!   blip_module
!!  PURPOSE
!!   Stores various variables relating to blips
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   16/05/01
!!  MODIFICATION HISTORY
!!   28/05/2001 dave
!!    Added ROBODoc header
!!   18/03/2002 dave
!!    Added RCS id and log tags in header and static RCS id for object file
!!   04/07/2002 mike
!!   Added arrays region_single and region_double, defined as follows:
!!    the location of a blip function on the blip-grid is specified
!!    by the integer triplet (nx,ny,nz), with (0,0,0) at the centre
!!    of the support region. All the integers are in the range
!!  -bliparraysize <= nx,ny,nz <= bliparraysize. Now fix the
!!    magnitude of one of the integers to be na, with 0 <= na <= bliparraysize.
!!    Then the maximum magnitude of either of the other two integers
!!    is region_single(na). With one of the integers fixed at na
!!    with 0 <= na <= bliparraysize, and another of the integers fixed
!!    at nb, with 0 <= nb <= region_single(na), the maximum magnitude
!!    of the third integer is given by region_double(na,nb). 
!!   11:30, 12/11/2004 dave 
!!    Changed to use proper variables from maxima_module
!!   2006/07/12 16:56 dave
!!    Changing variables to dynamically allocatable
!!   2006/07/20 08:33 dave
!!    Removed use common, tidied - all variables now dynamically allocated, basis functions
!!    stored in derived types in support_spec_format
!!   2006/10/20 14:01 dave
!!    Included make_pre and set_blip_index
!!   2006/11/01 17:04 dave
!!    Added alpha, beta and gauss2blip
!!   2008/02/04 08:28 dave
!!    Changed for output to file not stdout
!!   2008/05/28 ast
!!    Added timers
!!   2011/11/15 07:59 dave
!!    Updating the derived type for blips
!!  SOURCE
!!
module blip

  use datatypes
  use global_module, ONLY: io_lun
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_basis,tmr_std_allocation

  implicit none

  save

  ! RCS tag for object file identification 
  character(len=80), private :: RCSid = "$Id$"

  real(double) :: alpha, beta

  type blip_data
     integer, pointer, dimension(:,:) :: blip_location ! dim( 3, n_blips)
     integer, pointer, dimension(:,:,:) :: blip_number ! dim(-bliparraysize:bliparraysize x3 )
     integer :: BlipArraySize
     integer :: NBlipsRegion
     integer :: OneArraySize, FullArraySize
     integer :: Extent
     real(double) :: SupportGridSpacing, RegionRadius
     real(double) :: BlipWidth
     real(double) :: FourOnBlipWidth
     ! These last two are temporary and will be removed once G2B is working
     integer, pointer, dimension(:) :: region_single   ! (0:bliparraysize)
     integer, pointer, dimension(:,:) :: region_double ! (0:bliparraysize,0:bliparraysize)
  end type blip_data

  type blip_precon
     integer :: size
     real(double), pointer, dimension(:,:) :: coeffs
  end type blip_precon

  type(blip_data), allocatable, dimension(:) :: blip_info
  type(blip_precon), allocatable, dimension(:) :: PreCond !(n_blips_region, n_blips_region )
  character(len=10) :: init_blip_flag
!!***

contains

!!****f* blip/make_pre *
!!
!!  NAME 
!!   make_pre
!!  USAGE
!! 
!!  PURPOSE
!!   Calculates the preconditioning matrix for length-scale
!!    ill conditioning
!!  INPUTS
!! 
!! 
!!  USES
!!   datatypes, numbers, global_module, blip, dimens, GenComms, GenBlas
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   27/2/98
!!  MODIFICATION HISTORY
!!   10/05/01 DRB
!!    Converted to ROBODoc format and F90
!!   28/05/2001 dave
!!    Stripped subroutine calls
!!   08/06/2001 dave
!!    Added RCS Id and Log tags and GenComms for my_barrier
!!    Misc tidying
!!   11/06/2001 dave
!!    Changed dpot calls to generic ones through GenBlas
!!   25/06/2002 drb
!!    Tweaked headers, added USES line
!!   15:40, 2003/06/09 dave
!!    Added flag to switch preconditioning on/off
!!   2006/09/13 08:06 dave
!!    Changed to use species-dependent parameters from blip
!!   2006/10/20 13:57 dave
!!    Added memory registration and fixed small bug (stat not passed to allocate)
!!   2006/10/20 14:00 dave
!!    Included into blip
!!   2008/05/28 ast
!!    Added timers
!!  SOURCE
!!
  subroutine make_pre(inode, ionode) 

    use datatypes
    use numbers
    use global_module, ONLY: iprint_basis, area_basis, flag_precondition_blips
    use GenComms, ONLY: my_barrier, cq_abort
    use GenBlas, ONLY: potrf, potri
    use species_module, ONLY: n_species
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    integer :: inode, ionode

    ! Local variables
    integer :: i, N, n_blip1, n_blip2, nx, n_blip3,  &
         nx1, nx2, ny, ny1, ny2, nz, nz1, nz2

    real(double) :: p0(0:3),y0(0:3),p2(0:3),asq, sum
    real(double), allocatable, dimension(:,:)  :: Kmat
    !real(double) :: Kmat(n_blips_region, n_blips_region)
    real(double), parameter :: k0=1.3229425_double
    integer :: n_blips,info,stat, spec

    call start_timer(tmr_std_basis)
    if(flag_precondition_blips) then
       stat=0
       call start_timer(tmr_std_allocation)
       allocate(PreCond(n_species),STAT=stat)
       if(stat/=0) call cq_abort('ERROR allocating PreCond: ',n_species)
       call stop_timer(tmr_std_allocation)
       do spec = 1,n_species
          call start_timer(tmr_std_allocation)
          allocate(PreCond(spec)%coeffs(blip_info(spec)%NBlipsRegion,blip_info(spec)%NBlipsRegion), &
               Kmat(blip_info(spec)%NBlipsRegion,blip_info(spec)%NBlipsRegion),STAT=stat)
          if(stat/=0) call cq_abort('ERROR allocating PreCond and Kmat: ',blip_info(spec)%NBlipsRegion)
          call reg_alloc_mem(area_basis,2*blip_info(spec)%NBlipsRegion*blip_info(spec)%NBlipsRegion,type_dbl)
          call stop_timer(tmr_std_allocation)
          PreCond(spec)%size = blip_info(spec)%NBlipsRegion
          PreCond(spec)%coeffs = zero
          Kmat = zero
          if(inode==ionode.AND.iprint_basis>0) write(io_lun,fmt='(6x,"Applying preconditioning")')
          asq = one/(k0*k0*blip_info(spec)%SupportGridSpacing*blip_info(spec)%SupportGridSpacing)
          ! Set up the coefficients for the blip multiples
          p0(0)=151.0_double/140.0_double
          p0(1)=1191.0_double/2240.0_double
          p0(2)=3.0_double/56.0_double
          p0(3)=1.0_double/2240.0_double
          y0(0)=3.0_double/2.0_double
          y0(1)=-9.0_double/32.0_double
          y0(2)=-9.0_double/20.0_double
          y0(3)=-3.0_double/160.0_double
          p2(0)=y0(0)/p0(0)
          p2(1)=y0(1)/p0(1)
          p2(2)=y0(2)/p0(2)
          p2(3)=y0(3)/p0(3)
          ! Build the approximate Hessian
          call my_barrier()
          do n_blip1 = 1, blip_info(spec)%NBlipsRegion
             do n_blip2 = 1, blip_info(spec)%NBlipsRegion

                nx1 = blip_info(spec)%blip_location( 1, n_blip1 )
                ny1 = blip_info(spec)%blip_location( 2, n_blip1 )
                nz1 = blip_info(spec)%blip_location( 3, n_blip1 )

                nx2 = blip_info(spec)%blip_location( 1, n_blip2 )
                ny2 = blip_info(spec)%blip_location( 2, n_blip2 )
                nz2 = blip_info(spec)%blip_location( 3, n_blip2 )

                nx = iabs(nx1-nx2)
                ny = iabs(ny1-ny2)
                nz = iabs(nz1-nz2)
                if((nx.le.3).AND.(ny.LE.3).AND.(nz.LE.3)) then
                   Kmat(n_blip1,n_blip2) = p0(nx)*p0(ny)*p0(nz)* &
                        (one+asq*(p2(nx)+p2(ny)+p2(nz)))
                else
                   Kmat(n_blip1,n_blip2) = zero
                endif
             end do
          end do
          ! Invert
          ! Cholesky decomposition
          call my_barrier()
          info = 0
          do n_blip1 = 1,blip_info(spec)%NBlipsRegion
             do n_blip2 = 1,blip_info(spec)%NBlipsRegion
                PreCond(spec)%coeffs(n_blip1,n_blip2) = Kmat(n_blip1,n_blip2)
             enddo
          enddo
          call potrf('U',blip_info(spec)%NBlipsRegion,Kmat,blip_info(spec)%NBlipsRegion,info)
          if(info/=0) call cq_abort("Problem with potrf in preconditioning: ",info)
          ! Inversion
          call my_barrier()
          info = 0
          call potri('U',blip_info(spec)%NBlipsRegion,Kmat,blip_info(spec)%NBlipsRegion,info)
          if(info/=0) call cq_abort("Problem with potri in preconditioning: ",info)
          ! Now make the upper triangle a full matrix
          call my_barrier()
          do n_blip1 = 1, blip_info(spec)%NBlipsRegion
             do n_blip2 = n_blip1, blip_info(spec)%NBlipsRegion
                Kmat(n_blip2,n_blip1) = Kmat(n_blip1,n_blip2)
             enddo
          enddo
          do n_blip1 = 1, blip_info(spec)%NBlipsRegion
             do n_blip2 = 1, blip_info(spec)%NBlipsRegion
                PreCond(spec)%coeffs(n_blip1,n_blip2) = Kmat(n_blip1,n_blip2)
             enddo
          enddo
          deallocate(Kmat)
          call reg_dealloc_mem(area_basis,blip_info(spec)%NBlipsRegion*blip_info(spec)%NBlipsRegion,type_dbl)
       end do
    else
       if(inode==ionode.AND.iprint_basis>0) &
            write(io_lun,fmt='(6x,"Not applying preconditioning")') 
       !do n_blip1 = 1, blip_info(spec)%NBlipsRegion
       !   do n_blip2 = 1, blip_info(spec)%NBlipsRegion
       !      if(n_blip1.EQ.n_blip2) then
       !         PreCond(spec)%coeffs(n_blip1,n_blip2) = one
       !      else
       !         PreCond(spec)%coeffs(n_blip1,n_blip2) = zero
       !      endif
       !   enddo
       !enddo
    endif
    call stop_timer(tmr_std_basis)
    return
  end subroutine make_pre
!!***

!!****f* blip/set_blip_index *
!!
!!  NAME 
!!   set_blip_index
!!  USAGE
!! 
!!  PURPOSE
!!   Creates variables and index arrays characterising the set of
!!   blip functions wholly contained within the support region.
!!   Specifically:
!!   bliparraysize: the smallest integer such that for all blip
!!     functions wholly contained within the support region, these
!!     blip functions being centred on blip-grid sites (nx,ny,nz),
!!     the absolute values of nx, ny and nz never exceed bliparraysize.
!!   n_blips_region: total number of blip functions wholly contained
!!     within the support region.
!!   blip_number(nx,ny,nz): gives the sequence number of blip function
!!     at site (nx,ny,nz) in the enumeration in which nx varies
!!     most rapidly and nz least rapidly.
!!   blip_location(i,n): for blip function having sequence number n,
!    gives its blip-grid site (nx,ny,nz), with components nx, ny and
!!     nz corresponding to index i = 1, 2, 3 respectively.
!!   region_single(na), region_double(na,nb): defined as follows:
!!    the location of a blip function on the blip-grid is specified
!!    by the integer triplet (nx,ny,nz), with (0,0,0) at the centre
!!    of the support region. All the integers are in the range
!!  -bliparraysize <= nx,ny,nz <= bliparraysize. Now fix the
!!    magnitude of one of the integers to be na, with 0 <= na <= bliparraysize.
!!    Then the maximum magnitude of either of the other two integers
!!    is region_single(na). With one of the integers fixed at na
!!    with 0 <= na <= bliparraysize, and another of the integers fixed
!!    at nb, with 0 <= nb <= region_single(na), the maximum magnitude
!!    of the third integer is given by region_double(na,nb). 
!!  INPUTS
!!   inode: id number of node running executable
!!   ionode: id number of input/output node
!! 
!!  USES
!! 
!!  AUTHOR
!!   E.H.Hernandez
!!  CREATION DATE
!!   14/04/95
!!  MODIFICATION HISTORY
!!   16/05/2001 dave
!!    Converted to F90, added ROBODoc header and reduced number
!!    of passed variables
!!   28/05/2001 dave
!!    Removed r_H as passed variables and included from dimens
!!   20/06/2001 dave
!!    Added cq_abort and RCS Id and Log tags
!!   25/06/2002 mike
!!    Code simplified, protection against array-bound overflow
!!    improved, comments added
!!   29/06/2002 mike
!!    Code added to create arrays region_single and region_double 
!!   2006/10/20 14:00 dave
!!    Included into blip
!!   2008/05/28 ast
!!    Added timers
!!   2011/11/17 10:38 dave
!!    Changes for new blip data (removing allocation)
!!  SOURCE
!!
  subroutine set_blip_index(inode,ionode)

    use datatypes
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: iprint_init, area_basis
    use numbers
    use support_spec_format, ONLY : supports_on_atom, support_gradient, support_elec_gradient, &
         coefficient_array,grad_coeff_array, elec_grad_coeff_array, coeff_array_size,&
         allocate_supp_coeff_array, associate_supp_coeff_array
    use species_module, ONLY: nsf_species, n_species
    use primary_module, ONLY: bundle
    use dimens, ONLY: RadiusSupport
    use memory_module, ONLY: reg_alloc_mem, type_dbl, type_int

    implicit none

    ! Variables passed in argument list
    integer, intent(in) :: inode, ionode

    ! Local variables
    integer :: i, n_blips, na, nb, nx, ny, nz, stat, size, this_nsf, spec
    real(double) :: a2, b2, dx, dy, dz, r2_over_b2


    call start_timer(tmr_std_basis)
    if((inode == ionode).and.(iprint_init >= 2)) then
       write(unit=io_lun,fmt='(//10x,65("*")/25x,"REPORT FROM SET_BLIP_INDEX"/&
            &10x,65("*"))')
    end if

    do spec = 1, n_species
       ! It is required that the full width of the B-spline blip-function
       ! be exactly four times the blip-grid spacing. This condition
       ! is checked here, and the calculation stops if the condition
       ! is not satisfied:
       !blip_info(spec)%BlipWidth = four*blip_info(spec)%SupportGridSpacing
       if(inode==ionode.AND.iprint_init>2) write(io_lun,fmt='(10x,"Blip width: ",f20.12)') blip_info(spec)%BlipWidth
       blip_info(spec)%FourOnBlipWidth = four/blip_info(spec)%BlipWidth
       !if(abs(blip_width - four*support_grid_spacing) > very_small) then
       !   call cq_abort('set_blip_index: blip width must be exactly four &
       !        &times blip-grid spacing')
       !end if

       ! check that the relation between region radius and blip-grid
       ! spacing is such that at least one blip function is
       ! wholly contained within the region (RadiusSupport(spec) is support-region radius,
       ! support_grid_spacing is blip-grid spacing):
       r2_over_b2 = (RadiusSupport(spec)*RadiusSupport(spec))/ &
            (blip_info(spec)%SupportGridSpacing*blip_info(spec)%SupportGridSpacing)
       if(r2_over_b2 < 12.0_double) then
          call cq_abort('set_blip_grid: support region too small to hold&
               & one blip function')
       end if

       ! For the meaning of the integer bliparraysize, see header. Here,
       ! bliparraysize is calculated by formula, printed, and
       ! tested against the relevant array bound:
       blip_info(spec)%BlipArraySize = sqrt(r2_over_b2 - 8.0_double) - 2.0_double
       if((inode == ionode).and.(iprint_init >= 2)) then
          write(unit=io_lun,fmt='(/10x," set_blip_index: bliparraysize:",i5/10x&
               &"(This is the smallest integer such that for all blip functions"/10x&
               &" wholly contained within the support region, these blip"/10x&
               &" functions being centred on sites (nx,ny,nz), the absolute"/10x&
               &" values of nx, ny and nz never exceed this integer)")') &
               &blip_info(spec)%BlipArraySize
       end if
       blip_info(spec)%OneArraySize = 1 + 2*(blip_info(spec)%BlipArraySize+4)
       blip_info(spec)%FullArraySize = blip_info(spec)%OneArraySize*blip_info(spec)%OneArraySize*&
            blip_info(spec)%OneArraySize
       !if(bliparraysize > MAXBAS) then
       !   call cq_abort('set_blip_index: MAXBAS too small',bliparraysize,MAXBAS)
       !end if
       call start_timer(tmr_std_allocation)
       allocate(blip_info(spec)%region_single(0:blip_info(spec)%BlipArraySize),&
            blip_info(spec)%region_double(0:blip_info(spec)%BlipArraySize,0:blip_info(spec)%BlipArraySize),&
            STAT=stat)
       if(stat/=0) call cq_abort('Error allocating region_single and region_double: ',&
            blip_info(spec)%BlipArraySize)
       call reg_alloc_mem(area_basis,(blip_info(spec)%BlipArraySize+1)*(blip_info(spec)%BlipArraySize+2),type_int)
       call stop_timer(tmr_std_allocation)
       ! Calculate the index arrays region_single, region_double
       do na = 0, blip_info(spec)%BlipArraySize
          a2 = (na+2)*(na+2)
          blip_info(spec)%region_single(na) = sqrt(r2_over_b2 - a2 - 4.0_double) - 2.0_double
          do nb = 0, blip_info(spec)%region_single(na)
             b2 = (nb+2)*(nb+2) + a2
             blip_info(spec)%region_double(na,nb) = sqrt(r2_over_b2 - b2) - 2.0_double
          end do
       end do

       n_blips = 0
       do nz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
          dz = 2 + abs(nz)
          do ny = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             dy = 2 + abs(ny)
             do nx = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
                dx = 2 + abs(nx)
                if (dx*dx + dy*dy + dz*dz <= r2_over_b2) n_blips = n_blips + 1
             end do
          end do
       end do
       ! For meaning of the array blip_info(nx,ny,nz), see description
       ! at head of routine. The entire array
       ! is initialised to zero here to maintain compatibility with
       ! the original version of the code.

       call start_timer(tmr_std_allocation)
       allocate(blip_info(spec)%blip_number(-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,&
            -blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,&
            -blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize),STAT=stat)
       if(stat/=0) call cq_abort('Error allocating blip_number: ',blip_info(spec)%BlipArraySize)
       call reg_alloc_mem(area_basis,(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%BlipArraySize+1)*&
            (2*blip_info(spec)%BlipArraySize+1),type_int)
       call stop_timer(tmr_std_allocation)
       blip_info(spec)%blip_number = 0

       call start_timer(tmr_std_allocation)
       allocate(blip_info(spec)%blip_location(3,n_blips),STAT=stat)
       if(stat/=0) call cq_abort('Error allocating blip_number and blip_location: ',n_blips)
       call reg_alloc_mem(area_basis,3*n_blips,type_int)
       call stop_timer(tmr_std_allocation)
       blip_info(spec)%NBlipsRegion = n_blips
       n_blips = 0
       do nz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
          dz = 2 + abs(nz)
          do ny = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             dy = 2 + abs(ny)
             do nx = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
                dx = 2 + abs(nx)
                if (dx*dx + dy*dy + dz*dz <= r2_over_b2) then
                   n_blips = n_blips + 1
                   blip_info(spec)%blip_number(nx,ny,nz) = n_blips
                   blip_info(spec)%blip_location(1,n_blips) = nx
                   blip_info(spec)%blip_location(2,n_blips) = ny
                   blip_info(spec)%blip_location(3,n_blips) = nz
                end if
             end do
          end do
       end do
       if((inode == ionode).and.(iprint_init >= 2)) then
          write(unit=io_lun,fmt='(/10x," set_blip_index: total number of blip &
               &functions wholly contained within support region:",i5)') &
               &blip_info(spec)%NBlipsRegion
       end if

       if((inode == ionode).and.(iprint_init > 2)) then
          write(unit=io_lun,fmt='(/10x," set_blip_index: index array region_single:"/)')
          do na = 0, blip_info(spec)%BlipArraySize
             write(unit=io_lun,fmt='(10x,i5,3x,i5)') na, blip_info(spec)%region_single(na)
          end do
          write(unit=io_lun,fmt='(/10x," set_blip_index: index array region_double:"/)')
          do na = 0, blip_info(spec)%BlipArraySize
             do nb = 0, blip_info(spec)%region_single(na)
                write(unit=io_lun,fmt='(10x,2i5,3x,i5)') na, nb, blip_info(spec)%region_double(na,nb)
             end do
          end do
          write(unit=io_lun,fmt='(/10x," set_blip_index: index array blip_number:"/)')
          do nz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             do ny = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
                do nx = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
                   write(unit=io_lun,fmt='(10x,3i5,3x,i5)') nx, ny, nz, &
                        &blip_info(spec)%blip_number(nx,ny,nz)
                end do
             end do
          end do
          write(unit=io_lun,fmt='(/10x," set_blip_index: index array blip_location:"/)')
          do i = 1, blip_info(spec)%NBlipsRegion
             write(unit=io_lun,fmt='(10x,i5,3x,3i5)') i, blip_info(spec)%blip_location(1,i), &
                  &blip_info(spec)%blip_location(2,i), blip_info(spec)%blip_location(3,i)
          end do
       end if
    end do ! Loop over species
    ! Set up data storage for atoms on this processor
    size = 0
    call start_timer(tmr_std_allocation)
    allocate(supports_on_atom(bundle%n_prim),support_gradient(bundle%n_prim),support_elec_gradient(bundle%n_prim))
    call stop_timer(tmr_std_allocation)
    do i=1,bundle%n_prim
       ! Check on species
       spec = bundle%species(i)
       this_nsf = nsf_species(spec)
       size = size + this_nsf*blip_info(spec)%NBlipsRegion
       supports_on_atom(i)%nsuppfuncs = this_nsf
       call start_timer(tmr_std_allocation)
       allocate(supports_on_atom(i)%supp_func(this_nsf),support_gradient(i)%supp_func(this_nsf),&
            support_elec_gradient(i)%supp_func(this_nsf))
       call stop_timer(tmr_std_allocation)
       supports_on_atom(i)%supp_func(:)%ncoeffs = blip_info(spec)%NBlipsRegion
       support_gradient(i)%nsuppfuncs = this_nsf
       support_gradient(i)%supp_func(:)%ncoeffs = blip_info(spec)%NBlipsRegion
       support_elec_gradient(i)%nsuppfuncs = this_nsf
       support_elec_gradient(i)%supp_func(:)%ncoeffs = blip_info(spec)%NBlipsRegion
    end do
    coeff_array_size = size
    call allocate_supp_coeff_array(size)
    call associate_supp_coeff_array(supports_on_atom,bundle%n_prim,coefficient_array,size)
    call associate_supp_coeff_array(support_gradient,bundle%n_prim,grad_coeff_array,size)
    call associate_supp_coeff_array(support_elec_gradient,bundle%n_prim,elec_grad_coeff_array,size)

    call stop_timer(tmr_std_basis)
    return

  end subroutine set_blip_index
!!***

!!****f* blip/gauss2blip *
!!
!!  NAME 
!!   gauss2blip
!!  USAGE
!! 
!!  PURPOSE
!!   Sets blip coefficients as values of a Gaussian function
!!   at the distance of each blip from the atom position. Can
!!   only be used if there are four supports on each atom,
!!   in which case the support have the form of teh s-function
!!   and three p-functions. The Gaussian exponents for s- and
!!   p-functions are alpha and beta (atomic units). This way
!!   of initiating blip coefficients is a relic of the very
!!   early history of Conquest, and its use is strongly
!!   discouraged, unless you know what you are doing. It
!!   continues to be included in Conquest mainly in order
!!   to maintain compatibility with earlier version.
!!
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   E.H.Hernandez
!!  CREATION DATE
!!   11/2/95
!!  MODIFICATION HISTORY
!!   10/05/01 DRB
!!    Added ROBODoc header and updated to F90
!!   28/05/2001 dave
!!    Stripped subroutine call
!!   13/06/02 MJG
!!    Name of routine changed from set_support to gauss2blip
!!    and comments added. In addition, the number of
!!    supports on each atom is checked, and the code aborts
!!    if this number is not equal to 4. fortran coding tidied up.
!!   15:48, 04/02/2003 drb 
!!    Changed output so that it's only seen for iprint>2 and changed format statement for grid spacing
!!   2006/09/13 07:49 dave
!!    Changed loop to be over primary set and added species dependence
!!   2006/11/01 17:06 dave
!!    Included in blip
!!   2008/03/03 18:41 dave
!!    Changed float to real (probably unnecessary)
!!   2008/05/28 ast
!!    Added timers
!!   2008/09/17 10:01 dave
!!    Added nsf=1 case 
!!  SOURCE
!!
  subroutine gauss2blip

    use datatypes
    use numbers
    use global_module, ONLY: iprint_basis
    use GenComms, ONLY: cq_abort, inode, ionode
    use support_spec_format, ONLY: coefficient_array, supports_on_atom
    use primary_module, ONLY: bundle

    implicit none

    ! Local variables
    integer :: i, n_blip, nx, ny, nz, spec
    real(double) :: gauss_1, gauss_2, r2, x, y, z

    ! variables accessed by use:
    ! -------------------------------------------
    !   variable                  module
    ! -------------------------------------------
    !   blip_info(spec)%blip_location:  blip
    !   data_blip:                      blip
    !   iprint:                         global module
    !   blip_info(spec)%NBlipsRegion:             blip
    !   blip_info(spec)%SupportGridSpacing:       blip
    ! -------------------------------------------

    ! since the initialisation method provided in this routine
    ! can only be used if there are 4 support functions on
    ! every atom, this condition is checked. It is a fatal
    ! error if the condition is not satisfied, and the
    ! calculation is aborted.

    call start_timer(tmr_std_basis)
    if((inode == ionode).and.(iprint_basis >= 0)) then
       write(unit=io_lun,fmt='(/10x," gauss2blip: sbrt entered")')
    end if

    !if((inode == ionode) .AND. (iprint_basis >= 0)) then
    !   write(unit=io_lun,fmt='(10x," gauss2blip: support-grid spacing:",f20.12)') &
    !        &blip_info(spec)%SupportGridSpacing
    !end if

    ! data_blip is the array holding the blip coeffs for the support
    ! support functions of atoms in the primary set of my processor. T
    ! maintain consistency with earlier versions of the code, the
    ! entire array is initialised to zero before values are
    ! assigned to the blip coeffs.

    coefficient_array = zero

    ! loop over all atoms in my primary set, and over all blips in
    ! support region of each atom. At present, all species of atoms
    ! have the same support region and the same blip grid.
    do i = 1, bundle%n_prim
       spec = bundle%species(i)
       do n_blip = 1, blip_info(spec)%NBlipsRegion
          nx = blip_info(spec)%blip_location( 1, n_blip )
          ny = blip_info(spec)%blip_location( 2, n_blip )
          nz = blip_info(spec)%blip_location( 3, n_blip )
          x = real( nx, double ) * blip_info(spec)%SupportGridSpacing
          y = real( ny, double ) * blip_info(spec)%SupportGridSpacing
          z = real( nz, double ) * blip_info(spec)%SupportGridSpacing
          r2 = x * x + y * y + z * z
          gauss_1 = exp( - alpha * r2)
          gauss_2 = exp( - beta * r2)
          if(supports_on_atom(i)%nsuppfuncs>=4) then
             supports_on_atom(i)%supp_func(1)%coefficients(n_blip) = gauss_1
             supports_on_atom(i)%supp_func(2)%coefficients(n_blip) = x * gauss_2
             supports_on_atom(i)%supp_func(3)%coefficients(n_blip) = y * gauss_2
             supports_on_atom(i)%supp_func(4)%coefficients(n_blip) = z * gauss_2
             if(iprint_basis>2.AND.inode==ionode) then
                write(io_lun,fmt='(10x,"Blip Values: ",4i4,f20.12)') inode, i, 1, n_blip, &
                     supports_on_atom(i)%supp_func(1)%coefficients(n_blip)
                write(io_lun,fmt='(10x,"Blip Values: ",4i4,f20.12)') inode, i, 2, n_blip, &
                     supports_on_atom(i)%supp_func(2)%coefficients(n_blip)
                write(io_lun,fmt='(10x,"Blip Values: ",4i4,f20.12)') inode, i, 3, n_blip, &
                     supports_on_atom(i)%supp_func(3)%coefficients(n_blip)
                write(io_lun,fmt='(10x,"Blip Values: ",4i4,f20.12)') inode, i, 4, n_blip, &
                     supports_on_atom(i)%supp_func(4)%coefficients(n_blip)
             end if
          else if (supports_on_atom(i)%nsuppfuncs==1) then
             supports_on_atom(i)%supp_func(1)%coefficients(n_blip) = gauss_1
             if(iprint_basis>2.AND.inode==ionode) then
                write(io_lun,fmt='(10x,"Blip Values: ",4i4,f20.12)') inode, i, 1, n_blip, &
                     supports_on_atom(i)%supp_func(1)%coefficients(n_blip)
             end if
          else
             call cq_abort("With this version of Conquest, gaussian initialisation needs 1 or 4 support functions ",&
                  supports_on_atom(i)%nsuppfuncs)
          end if
       end do
    end do
    call stop_timer(tmr_std_basis)
    return
  end subroutine gauss2blip
!!***
  
end module blip
