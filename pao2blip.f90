! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module pao2blip
! ------------------------------------------------------------------------------
! Code area 11: Basis operations
! ------------------------------------------------------------------------------

!!****h* Conquest/pao2blip *
!!  NAME
!!   pao2blip
!!  PURPOSE
!!   The sole purpose of the present module is to contain
!!   the sbrt make_blips_from_paos, whose purpose is described
!!   below.
!!
!!  USES
!!
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22 June 2002
!!  MODIFICATION HISTORY
!!   12/08/2002 mjg
!!    Changed names of members of pao derived type to reflect reorganisation 
!!    for dynamic allocation of memory
!!   05/09/2002 mjg & drb 
!!    Removed read for PAO data and moved to separate file
!!   2006/09/13 08:27 dave
!!    Removed maxima from common
!!   2008/02/06 08:09 dave
!!    Changed for output to file not stdout
!!  SOURCE
module pao2blip

  use global_module, ONLY: io_lun

  implicit none
  save

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id$"

!!***

contains

! -------------------------------------------------------------
! subroutine make_blips_from_paos
! -------------------------------------------------------------

!!****f* pao2blip/make_blips_from_paos *
!!
!!  NAME 
!!   make_blips_from_paos
!!  USAGE
!! 
!!  PURPOSE
!!   This is the top routine of a suite of routines whose purpose
!!   is to calculate the blip coefficients for a set of support
!!   functions. Each of the support functions is specified
!!   by its representation in terms of pseudo-atomic orbitals
!! (PAO's), and the job of the present suite of routines is to
!!   determine the blip coefficients that give the best least-
!!   squares fit to this PAO specification. The PAO data itself
!!   is read by a call from the sbrt make_blips_from_paos to
!!   the routine read_pao, and the specification of the supports
!!   in terms of PAO's is read by a call from make_blips_to_paos
!!   to the routine read_support. The sbrt make_blips_from_paos
!!   is contained in the present module pao2blip. The sbrt's
!!   read_pao and read_support are contained in the modules
!!   read_pao_info and read_support_spec respectively. 
!!
!!  INPUTS
!!   inode: identification of node running the executable
!!   ionode: identification number of input/ouput node
!!   n_species: no. of species of atom in the simulated system
!!   r_h: cut-off radius of the support region
!!   support_grid_spacing: spacing of grid on which blip-functions sit
!!  USES
!!
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22 June 2002
!!  MODIFICATION HISTORY
!!   12/08/2002 mjg
!!    Changed names of members of pao derived type to reflect reorganisation 
!!    for dynamic allocation of memory
!!   14/08/2002 mjg
!!    Call to read_pao_info removed, because reading of pao's
!!    is now done elsewhere by the routine get_ionic_data
!!  SOURCE
  subroutine make_blips_from_paos(inode,ionode,n_species)

    use blip, ONLY: NBlipsRegion, SupportGridSpacing, BlipArraySize, blip_info
    use blip_pao_values
    use datatypes
    use dimens, ONLY : RadiusSupport
    use GenComms, ONLY : cq_abort, gcopy
    use global_module, ONLY : iprint_basis
    use linear_equation_solver
    use numbers
    use pao_format
    use primary_module, ONLY : bundle
    use read_support_spec, ONLY: read_support, support_info
    use symmetry
    use species_module, ONLY: nsf_species
    use support_spec_format, ONLY: supports_on_atom
    use maxima_module, ONLY: maxnsf, max_blip_nu_int

    implicit none

    integer, intent(in) :: inode, ionode, n_species
    integer :: i, i_stop, n, na, nb, nx, ny, nz, ns, n_ac, n_acz, n_am, &
         &n_s, n_sp, n_sup, n_zeta, min_ac, max_ac, l, m, count
    integer :: n_blip, nu_int, stat
    integer, allocatable, dimension(:) :: n_cube2sphere, n_b_half
    real(double) :: c, deltax, r2_over_b2, x
    real(double), dimension(max_blip_nu_int,4) :: bv
    real(double), allocatable, dimension(:) :: scal_prod_sym
    real(double), allocatable, dimension(:) :: blip_contrib

    type blip_store
       real(double), pointer, dimension(:,:) :: coeff
    end type blip_store
    type(blip_store), allocatable, dimension(:) :: blip_coeff

    type scan_pao
       integer, dimension(:,:), pointer :: n_support
       integer, dimension(:,:,:), pointer :: which_support
       real(double), dimension(:,:,:), pointer :: coeff
    end type scan_pao

    type(scan_pao), dimension(:), allocatable :: this_scan

    if((inode == ionode).and.(iprint_basis >= 2)) then
       write(unit=io_lun,fmt='(//" make_blips_from_paos: sbrt entered")')
    end if

    ! number of species
    if((inode == ionode).and.(iprint_basis >= 1)) then
       write(unit=io_lun,fmt='(//" make_paos_from_blips:&
            & no. of species:",i5)') n_species
    end if
    allocate(blip_coeff(n_species),n_b_half(n_species), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating blip_coeff in pao2blip: ",n_sp)

    ! N.B. This should be species-dependent
    ! Message to the future code developer:
    ! All the routines for obtaining blip coefficients from PAO data
    ! are constructed on the assumption that the region radii and
    ! blip-grid intervals are generally different for different
    ! species. However, at the time these routines are being
    ! put into Conquest, the latter requires that region radii and
    ! blip-grid spacings are the same for all species. Quantities
    ! such as RadiusSupport (the support-region radius), support_grid_spacing
    ! (the blip-grid spacing), etc., are therefore just single
    ! variables, independent of species number. The generalisation
    ! is straightforward: just make these quantities into
    ! arrays RadiusSupport(n_species), etc... and use them inside the
    ! loops over species number that occur in the present
    ! routine. To reflect the intended level of generality,
    ! some of the loops over species number in this routine
    ! are actually redundant, but they are written like this
    ! to ease the job of the future developer.... 

    nu_int = max_blip_nu_int

    ! loop over species
    do n_sp = 1, n_species
       allocate(blip_coeff(n_sp)%coeff(NBlipsRegion(n_sp), maxnsf),STAT=stat)
       if(stat/=0) call cq_abort("Error allocating blip_coeff in pao2blip: ",NBlipsRegion(n_sp),maxnsf)
       if((inode == ionode).and.(iprint_basis >= 1)) then
          write(unit=io_lun,fmt='(//" make_blips_from_paos: region and&
               & blip data for species no.",&
               &i3,":"/)') n_sp
          write(unit=io_lun,fmt='(" region radius:",f12.6)') RadiusSupport(n_sp)
          write(unit=io_lun,fmt='(" length of blip-grid interval:",f12.6)') &
               &SupportGridSpacing(n_sp)
          write(unit=io_lun,fmt='(" no. of integration-grid intervals in",&
               &" blip-grid interval:",i5)') nu_int
       end if
    end do

    ! check input data and make constants for later use
    do n_sp = 1, n_species
       if((inode == ionode).and.(iprint_basis >= 2)) then
          write(unit=io_lun,fmt='(/" make_blips_from_paos: checking input",&
               &" for species no:",i3)') n_sp
       end if
       ! check that there are enough blip intervals in region radius
       r2_over_b2 = (RadiusSupport(n_sp)*RadiusSupport(n_sp))/&
            &(SupportGridSpacing(n_sp)*SupportGridSpacing(n_sp))
       if(r2_over_b2 <= 12.0_double) then
          call cq_abort('make_blips_from_paos: support region too&
               & small to hold one blip')
       end if
       n_b_half(n_sp) = sqrt(r2_over_b2 - 8.0_double) - 2.0_double + very_small
       n_blip = 1 + 2*n_b_half(n_sp)
       if((inode == ionode).and.(iprint_basis >= 1)) then
          write(unit=io_lun,fmt='(" no. of entire blips in region diameter:",&
               &i5)') n_blip
       end if
       if(n_b_half(n_sp) > BlipArraySize(n_sp)) then
          call cq_abort('make_blips_from_paos: value of n_b_half&
               & gives too many blips',n_b_half(n_sp),BlipArraySize(n_sp))
       end if
    end do

    ! Here, in previous versions, there was a call to
    ! read_pao_info. Now removed, because read of pao's is
    ! now done earlier in the routine get_ionic_data

    ! read specification of all support functions in terms of PAO's
    call read_support(inode,ionode,n_species)

    ! check that specification of support functions is compatible
    ! with PAO data. NOTE: at the time of putting this routine
    ! into Conquest (3/07/02), the number of support functions for
    ! each species is hard-wired into the code as the parameter NSF.
    ! At the earliest possible opportunity, this should be replaced
    ! by a variable set in input, which could be called
    ! no_support_in_species(n_sp).
!** Change this to loop over atoms if it remains **! 
    do n_sp = 1, n_species
       do n_s = 1, nsf_species(n_sp)
          do n_acz = 1, support_info(n_sp,n_s)%no_of_acz
             n_ac = support_info(n_sp,n_s)%acz_info(n_acz)%which_ac
             if(n_ac == 1) then
                n_am = 0
             else if((n_ac >= 2) .and. (n_ac <= 4)) then
                n_am = 1
             else if((n_ac >= 5) .and. (n_ac <= 9)) then
                n_am = 2
             else
                if((inode == ionode).and.(iprint_basis >= 1)) then
                   write(unit=io_lun,fmt='(//" n_sp:",i3," n_s:",&
                        &i3," n_acz:",i3)') n_sp, n_s, n_acz
                end if
                call cq_abort('make_blips_from_paos: n_ac out of range',n_ac)
             end if
             if(n_am > pao(n_sp)%greatest_angmom) then
                if((inode == ionode).and.(iprint_basis >= 1)) then
                   write(unit=io_lun,fmt='(//" n_sp:",i3," n_s:",i3," n_acz:",&
                        &i3)') n_sp, n_s, n_acz
                end if
                call cq_abort('make_blips_from_paos: ang. mom. requested&
                     &by support spec too big')
             end if
             n_zeta = support_info(n_sp,n_s)%acz_info(n_acz)%which_zeta
             if(n_zeta > pao(n_sp)%angmom(n_am)%n_zeta_in_angmom) then
                if((inode == ionode).and.(iprint_basis >= 1)) then
                   write(unit=io_lun,fmt='(//" n_am:",i3," n_sp:",i3," n_s:",i3,&
                        &" n_acz:",i3)') n_am, n_sp, n_s, n_acz
                end if
                call cq_abort('make_blips_from_paos: support spec requested&
                     & too many zetas')
             end if
          end do
       end do
    end do

    ! -------------------------------------------------------------
    !   loop over species
    ! -------------------------------------------------------------
    allocate(this_scan(n_species),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating this_scan in pao2blip ",n_species)
    do n_sp = 1, n_species
       n_zeta = 0
       min_ac = 1
       do n_am = 0, pao(n_sp)%greatest_angmom
          max_ac = min_ac + 2*n_am
          min_ac = max_ac + 1
          if(pao(n_sp)%angmom(n_am)%n_zeta_in_angmom>n_zeta) n_zeta = pao(n_sp)%angmom(n_am)%n_zeta_in_angmom
       end do
       allocate(this_scan(n_sp)%n_support(max_ac,n_zeta),STAT=stat)
       allocate(this_scan(n_sp)%which_support(max_ac,n_zeta,maxnsf),STAT=stat)
       allocate(this_scan(n_sp)%coeff(max_ac,n_zeta,maxnsf),STAT=stat)
       ! initialise scan table to zero
       min_ac = 1
       do n_am = 0, pao(n_sp)%greatest_angmom
          max_ac = min_ac + 2*n_am
          do n_ac = min_ac, max_ac
             if(pao(n_sp)%angmom(n_am)%n_zeta_in_angmom > 0) then
                do n_zeta = 1, pao(n_sp)%angmom(n_am)%n_zeta_in_angmom
                   this_scan(n_sp)%n_support(n_ac,n_zeta) = 0
                end do
             end if
          end do
          min_ac = max_ac + 1
       end do

       ! build up scan table
       ! loop over support functions for current species
       do n_s = 1, nsf_species(n_sp)
          do n_acz = 1, support_info(n_sp,n_s)%no_of_acz
             n_ac = support_info(n_sp,n_s)%acz_info(n_acz)%which_ac
             n_zeta = support_info(n_sp,n_s)%acz_info(n_acz)%which_zeta
             this_scan(n_sp)%n_support(n_ac,n_zeta) = &
                  this_scan(n_sp)%n_support(n_ac,n_zeta) + 1
             n = this_scan(n_sp)%n_support(n_ac,n_zeta)
             this_scan(n_sp)%which_support(n_ac,n_zeta,n) = n_s
             this_scan(n_sp)%coeff(n_ac,n_zeta,n) = &
                  support_info(n_sp,n_s)%acz_info(n_acz)%coeff
          end do
          !%%!count = 1
          !%%!!   loop over angular components and over zeta's
          !%%!do l=0,support_info(n_sp)%lmax
          !%%!   n_ac = l+1
          !%%!   do n_zeta = 1,support_info(n_sp)%naczs(l)
          !%%!      do m=-l,l
          !%%!         this_scan(n_sp)%n_support(n_ac,n_zeta) = &
          !%%!              this_scan(n_sp)%n_support(n_ac,n_zeta) + 1
          !%%!         n = this_scan(n_sp)%n_support(n_ac,n_zeta)
          !%%!         this_scan(n_sp)%which_support(n_ac,n_zeta,n) = n_s
          !%%!         this_scan(n_sp)%coeff(n_ac,n_zeta,n) = &
          !%%!              support_info(n_sp)%supp_func(n_s)%coefficients(count)
          !%%!         count = count+1
          !%%!      end do
          !%%!   end do
          !%%!end do
       end do
       ! Temporary output
       !min_ac = 1
       !do n_s = 1,nsf_species(n_sp)
       !   do n_am = 0, pao(n_sp)%greatest_angmom
       !      max_ac = min_ac + 2*n_am
       !      do n_ac = min_ac, max_ac
       !         if(pao(n_sp)%angmom(n_am)%n_zeta_in_angmom > 0) then
       !            do n_zeta = 1, pao(n_sp)%angmom(n_am)%n_zeta_in_angmom
       !               write(io_lun,*) n_am,n_ac,n_zeta,this_scan(n_sp)%n_support(n_ac,n_zeta),this_scan(n_sp)%coeff(n_ac,n_zeta,n)
       !            end do
       !         end if
       !      end do
       !      min_ac = max_ac + 1
       !   end do
       !end do

       ! fetch values of blip function on integration grid
       call b_value(nu_int,bv)
       if((inode == ionode).and.(iprint_basis >= 2)) then
          write(unit=io_lun,fmt='(/" values of 4 sectors of blip function",&
               &" on integration grid:"/)')
          deltax = one/float(nu_int)
          do n = 1, nu_int
             x = n*deltax  
             write(unit=io_lun,fmt='(i3,2x,f12.6,2x,4e15.6)') n, x, &
                  &(bv(n,i), i = 1, 4)
          end do
       end if

       ! the variable n_blips_region is in module blips
       if((inode == ionode).and.(iprint_basis >= 1)) then
          write(unit=io_lun,fmt='(/" make_blips_from_paos: total no. of &
               &blips wholly within region:", &
               & i5)') NBlipsRegion(n_sp)
       end if
       ! initialise blip coefficients to zero
       do ns = 1, nsf_species(n_sp)
          do n = 1, NBlipsRegion(n_sp)
             blip_coeff(n_sp)%coeff(n,ns) = 0.0_double
          end do
       end do

       allocate(n_cube2sphere((BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+1)/2),STAT=stat)
       if(stat/=0) call cq_abort("Error allocating n_cube2sphere: ",(BlipArraySize(n_sp)+1))
       allocate(scal_prod_sym((BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+1)/2),STAT=stat)
       if(stat/=0) call cq_abort("Error allocating scal_prod_sym: ",(BlipArraySize(n_sp)+1))
       allocate(blip_contrib(NBlipsRegion(n_sp)), STAT=stat)
       if(stat/=0) call cq_abort("Error allocating blip_contrib in pao2blip: ",NBlipsRegion(n_sp))
       ! ... loop over angular momenta ...............................
       do n_am = 0, pao(n_sp)%greatest_angmom
          if(pao(n_sp)%angmom(n_am)%n_zeta_in_angmom > 0) then

             ! . . loop over zeta's . . . . . . . . . . . . . . . . . . . .
             do n_zeta = 1, pao(n_sp)%angmom(n_am)%n_zeta_in_angmom

                select case(n_am)

                   ! ----------------------------------------------
                   !   angular momentum = 0
                   ! ----------------------------------------------
                case(0)
                   if(this_scan(n_sp)%n_support(1,n_zeta) > 0) then
                      call blips_in_star(inode,ionode,"a1g",&
                           &n_sp,n_am,n_zeta,nu_int,&
                           &n_b_half(n_sp),SupportGridSpacing(n_sp),&
                           &blip_info(n_sp)%region_single,blip_info(n_sp)%region_double,&
                           &bv,n_cube2sphere,scal_prod_sym)
                      n_ac = 1
                      call star2sphere("a1g",n_ac,n_b_half(n_sp),&
                           &NBlipsRegion(n_sp),blip_info(n_sp)%blip_location,&
                           &n_cube2sphere,scal_prod_sym,&
                           &blip_contrib)
                      do ns = 1, this_scan(n_sp)%n_support(n_ac,n_zeta)
                         n_sup = this_scan(n_sp)%which_support(n_ac,n_zeta,ns)
                         c = this_scan(n_sp)%coeff(n_ac,n_zeta,ns)
                         do n = 1, NBlipsRegion(n_sp)
                            blip_coeff(n_sp)%coeff(n,n_sup) = &
                                 &blip_coeff(n_sp)%coeff(n,n_sup) + &
                                 &c*blip_contrib(n)
                         end do
                      end do
                   end if

                   ! -----------------------------------------------
                   !   angular momentum = 1
                   ! -----------------------------------------------
                case(1)
                   if((this_scan(n_sp)%n_support(2,n_zeta) > 0) .or. &
                        (this_scan(n_sp)%n_support(3,n_zeta) > 0) .or. &
                        (this_scan(n_sp)%n_support(4,n_zeta) > 0)) then
                      call blips_in_star(inode,ionode,"t1u",&
                           &n_sp,n_am,n_zeta,nu_int,&
                           &n_b_half(n_sp),SupportGridSpacing(n_sp),&
                           &blip_info(n_sp)%region_single,blip_info(n_sp)%region_double,&
                           &bv,n_cube2sphere,scal_prod_sym)
                      do n_ac = 2, 4
                         if(this_scan(n_sp)%n_support(n_ac,n_zeta) > 0) then
                            call star2sphere("t1u",n_ac,n_b_half(n_sp),&
                                 &NBlipsRegion(n_sp),blip_info(n_sp)%blip_location,&
                                 &n_cube2sphere,scal_prod_sym,&
                                 &blip_contrib)
                            do ns = 1, this_scan(n_sp)%n_support(n_ac,n_zeta)
                               n_sup = this_scan(n_sp)%which_support(n_ac,n_zeta,ns)
                               c = this_scan(n_sp)%coeff(n_ac,n_zeta,ns)
                               do n = 1, NBlipsRegion(n_sp)
                                  blip_coeff(n_sp)%coeff(n,n_sup) = &
                                       &blip_coeff(n_sp)%coeff(n,n_sup) + &
                                       &c*blip_contrib(n)
                               end do
                            end do
                         end if
                      end do
                   end if

                   ! ----------------------------------------------------
                   !   angular momentum = 2
                   ! ----------------------------------------------------
                case(2)
                   if((this_scan(n_sp)%n_support(5,n_zeta) > 0) .or. &
                        (this_scan(n_sp)%n_support(6,n_zeta) > 0) .or. &
                        (this_scan(n_sp)%n_support(7,n_zeta) > 0)) then
                      call blips_in_star(inode,ionode,"t2g",&
                           &n_sp,n_am,n_zeta,nu_int,&
                           &n_b_half(n_sp),SupportGridSpacing(n_sp),&
                           &blip_info(n_sp)%region_single,blip_info(n_sp)%region_double,&
                           &bv,n_cube2sphere,scal_prod_sym)
                      do n_ac = 5, 7
                         if(this_scan(n_sp)%n_support(n_ac,n_zeta) > 0) then
                            call star2sphere("t2g",n_ac,n_b_half(n_sp),&
                                 &NBlipsRegion(n_sp),blip_info(n_sp)%blip_location,&
                                 &n_cube2sphere,scal_prod_sym,&
                                 &blip_contrib)
                            do ns = 1, this_scan(n_sp)%n_support(n_ac,n_zeta)
                               n_sup = this_scan(n_sp)%which_support(n_ac,n_zeta,ns)
                               c = this_scan(n_sp)%coeff(n_ac,n_zeta,ns)
                               do n = 1, NBlipsRegion(n_sp)
                                  blip_coeff(n_sp)%coeff(n,n_sup) = &
                                       &blip_coeff(n_sp)%coeff(n,n_sup) + &
                                       &c*blip_contrib(n)
                               end do
                            end do
                         end if
                      end do
                   end if
                   if((this_scan(n_sp)%n_support(8,n_zeta) > 0) .or. &
                        (this_scan(n_sp)%n_support(9,n_zeta) > 0)) then
                      call blips_in_star(inode,ionode,"eg",&
                           &n_sp,n_am,n_zeta,nu_int,&
                           &n_b_half(n_sp),SupportGridSpacing(n_sp),&
                           &blip_info(n_sp)%region_single,blip_info(n_sp)%region_double,&
                           &bv,n_cube2sphere,scal_prod_sym)
                      do n_ac = 8, 9
                         if(this_scan(n_sp)%n_support(n_ac,n_zeta) > 0) then
                            call star2sphere("eg",n_ac,n_b_half(n_sp),&
                                 &NBlipsRegion(n_sp),blip_info(n_sp)%blip_location,&
                                 &n_cube2sphere,scal_prod_sym,&
                                 &blip_contrib)
                            do ns = 1, this_scan(n_sp)%n_support(n_ac,n_zeta)
                               n_sup = this_scan(n_sp)%which_support(n_ac,n_zeta,ns)
                               c = this_scan(n_sp)%coeff(n_ac,n_zeta,ns)
                               do n = 1, NBlipsRegion(n_sp)
                                  blip_coeff(n_sp)%coeff(n,n_sup) = &
                                       &blip_coeff(n_sp)%coeff(n,n_sup) + &
                                       &c*blip_contrib(n)
                               end do
                            end do
                         end if
                      end do
                   end if

                end select

             end do
             ! . . end loop over zeta's . . . . . . . . . . . . . . . . . .

          end if
       end do
       ! ... end loop over angular momenta ...........................
       deallocate(n_cube2sphere,STAT=stat)
       if(stat/=0) call cq_abort("Error deallocating n_cube2sphere: ",(BlipArraySize(n_sp)+1))
       deallocate(scal_prod_sym,STAT=stat)
       if(stat/=0) call cq_abort("Error deallocating scal_prod_sym: ",(BlipArraySize(n_sp)+1))
       deallocate(blip_contrib, STAT=stat)
       if(stat/=0) call cq_abort("Error deallocating blip_contrib in pao2blip: ",NBlipsRegion(n_sp))

       ! for current species, loop over support functions, and for
       ! each support function determine deviation of blip-fn repn from
       ! original pao repn.
       do n_sup = 1, nsf_species(n_sp)
          call deviation(inode,ionode,n_sp,n_sup,nu_int,&
               &n_b_half(n_sp),RadiusSupport(n_sp),SupportGridSpacing(n_sp),&
               &blip_info(n_sp)%region_single,blip_info(n_sp)%region_double,&
               &blip_info(n_sp)%blip_number,BlipArraySize(n_sp),bv,blip_coeff(n_sp)%coeff)
       end do

    end do ! n_sp
    ! -------------------------------------------------------------
    !   end loop over species
    ! -------------------------------------------------------------

    ! loop over atoms i in primary set of my node, look up species n_sp of
    ! each atom, and put the blip coeffs blip_coeff(n_blip,n_sup,n_sp)
    ! for each support-function n_sup of that species into array 
    ! data_blip(n_sup,n_blip,i).
    do i = 1, bundle%n_prim
       n_sp = bundle%species(i)
       do n_sup = 1, nsf_species(n_sp)
          ! loop over blip functions in the support region of this species
          do n = 1, NBlipsRegion(n_sp)
             supports_on_atom(i)%supp_func(n_sup)%coefficients(n) = blip_coeff(n_sp)%coeff(n,n_sup)
          end do
       end do
    end do
    do n_sp=1,n_species
       deallocate(this_scan(n_sp)%n_support,this_scan(n_sp)%which_support,this_scan(n_sp)%coeff)
       deallocate(blip_coeff(n_sp)%coeff)
    end do
    deallocate(blip_coeff,this_scan,n_b_half)
  end subroutine make_blips_from_paos
!!*** 

! -----------------------------------------------------------
! Subroutine blips_in_star
! -----------------------------------------------------------

!!****f* find_blip_coeffs/blips_in_star *
!!
!!  NAME 
!!   blips_in_star
!!  USAGE
!! 
!!  PURPOSE
!!   Manages the entire calculation of blip-coefficient representation
!!   of a single pseudo-atomic orbital (PAO) for requested
!!   species number n_sp, angular momentum n_am (s = 0, p = 1, d = 2)
!!   and zeta number (n_zeta = 1, 2,...). The main steps in this
!!   calculation correspond to the main subroutines called, as follows:
!! (1)  call star_coeffs: do analysis of stars of blip-grid points
!!   for the given symmetry, the latter being specified by the 
!!   string-valued parameter sym_type;
!! (2) call symmat: construct the block of the blip-overlap matrix
!!   belonging to the requested symmetry type; 
!! (3) call f3_value: together with the block of surrounding code, this
!!   calculates the vector of overlap integrals between the blip functions
!!   and the specified PAO;
!! (4) call lin_solv: solves the resulting set of linear equations
!!   A.x = b, with A the symmetry-adapted block of the blip-overlap
!!   matrix, b the vector of blip-PAO overlaps, and x the vector
!!   of optimal blip coeffs.
!! 
!!  INPUTS
!!   inode: id of node running the executable
!!   ionode: id of input/output node
!!   sym_type: string-valued parameter specifying the symmetry
!!     type of the PAO whose optimal blip representation is to
!!     be found: acceptable values are: "a1g" (totally symmetric
!!     representation, corresponding to ang. mom. = 0);
!!   "t1u" (ang. mom. = 1, x, y, or z); "t2g" (ang. mom. = 2, yz, zx, 
!!     or xy); "eg" (ang. mom. = 2, x^2 - y^2, or 3z^2 - r^2);
!!   n_sp: species number
!!   n_am: angular momentum (s: 0; p: 1, d: 2). Note that the specification
!!     is already made through the sym_type parameter, so n_am is
!!     in practice redundant. But in principle, the two parameters
!!     serve slightly different purposes: sym_type specifies the
!!     angular symmetry, while n_am governs which table of PAO's
!!     is selected.
!!   n_zeta: zeta number
!!   nu_int: number of integration-grid intervals in blip-grid interval.
!!   (NB: this integration grid is purely for internal purposes of
!!     the PAO-to-blip calculation, and has no connection with the
!!     integration grid used elsewhere in Conquest.)
!!   n_b_half: with blip-grid site specified by (nx,ny,nz), each of the
!!     integers lies in the range [-n_b_half,n_b_half] for blips
!!     wholly contained within the support region.
!!   n_int: the blip-grid interval.
!!   region_bound_single, region_bound_double: arrays specifying the 
!!     maximum values of cartesian indices for blips wholly contained 
!!     within the support region. (For a detailed description of these 
!!     index arrays, see the header to the routine "sphere".)
!!   bv: precalculated table of values of the B-spline.
!!   n_cube2sphere: index array relating sequence number of a blip-site
!!     in the support-region cube to the corresponding sequence number
!!     in the inscribed support-region sphere.
!!   scal_prod_sym: array holding the vector of blip coefficients.   
!! 
!!  USES
!! 
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!
!!
!!  SOURCE
  subroutine blips_in_star(inode,ionode,sym_type,&
       &n_sp,n_am,n_zeta,nu_int,n_b_half,b_int,region_bound_single,&
       &region_bound_double,bv,n_cube2sphere,scal_prod_sym)
    use blip_pao_values
    use blip, ONLY: BlipArraySize
    use maxima_module, ONLY : max_blip_nu_int
    use datatypes
    use GenComms, ONLY: cq_abort
    use global_module, ONLY : iprint_basis
    use linear_equation_solver
    use numbers
    use symmetry

    implicit none

    character(len=*), intent(in) :: sym_type
    integer, intent(in) :: inode, ionode, n_sp, &
         &n_am, n_zeta, nu_int, n_b_half
    integer, intent(in), dimension(0:) :: region_bound_single
    integer, intent(in), dimension(0:,0:) :: region_bound_double
    real(double), intent(in) :: b_int
    real(double), intent(in), dimension(:,:) :: bv
    integer, intent(out), dimension(:) :: n_cube2sphere
    real(double), intent(out), dimension(:) :: scal_prod_sym
    integer :: i, i_store, j, j_add, j_store, kb_min, kb_max, &
         &m, mult, n, n1, n2, n3, n_b_int_in_radius, n_b_int_in_diam, &
         &n_blip, n_blip_squared, n_star_in_cube, n_star_in_sphere
    integer, dimension((BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+2)/2) :: xsite, ysite, zsite
    real(double) :: delta_ig, deltax, sum, z0
    real(double), dimension((BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+2)/2) :: coeff
    real(double), dimension(4,4,2*BlipArraySize(n_sp)+1,2*BlipArraySize(n_sp)+1) :: store
    real(double), dimension((2*BlipArraySize(n_sp)+1)*(2*BlipArraySize(n_sp)+1)*(2*BlipArraySize(n_sp)+1)) :: scal_prod
    real(double), dimension((BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+2)/2,&
         (BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+2)/2) :: amat
    real(double), dimension(max_blip_nu_int,2*BlipArraySize(n_sp)+1,2*BlipArraySize(n_sp)+1) :: fv3

    ! constants
    n_blip = 1 + 2*n_b_half
    n_blip_squared = n_blip*n_blip
    n_b_int_in_radius = n_b_half + 2
    n_b_int_in_diam = 2*n_b_int_in_radius
    delta_ig = b_int/float(nu_int)
    deltax = one/float(nu_int)

    ! check enough blips for requested symmetry
    if(n_b_half < 0) then
       call cq_abort('blips_in_star: n_b_half cannot be negative',&
            &n_b_half)
    end if
    if((sym_type /= "a1g") .and. (n_b_half < 1)) then
       call cq_abort('blips_in_star: value of n_b_half gives too few&
            & blips for requested symmetry',n_b_half)
    end if
    if(n_b_half > BlipArraySize(n_sp)) then
       call cq_abort('blips_in_star: n_b_half too big',n_b_half)
    end if

    ! determine and check number of stars in cube
    select case(sym_type)
    case("a1g")
       n_star_in_cube = ((n_b_half+1)*(n_b_half+2)*(n_b_half+3))/6
    case("t1u")
       n_star_in_cube = ((n_b_half+1)*(n_b_half+2)*n_b_half)/2
    case("t2g")
       n_star_in_cube = ((n_b_half+1)*(n_b_half+1)*n_b_half)/2
    case("eg")
       n_star_in_cube = ((n_b_half+1)*(n_b_half+1)*n_b_half)/2
    case default
       if(inode == ionode) then
          call cq_abort('blips_in_star: unrecognised symmetry type')
       end if
    end select
    if((inode == ionode).and.(iprint_basis >= 2)) then
       write(unit=io_lun,fmt='(/" number of stars in cubic region:",i5)') &
            & n_star_in_cube
    end if
    if(n_star_in_cube > (BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+2)/2) then
       call cq_abort('blips_in_star: too many stars',n_star_in_cube,&
            &(BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+1)*(BlipArraySize(n_sp)+2)/2)
    end if

    ! do analysis of stars of blip-grid points
    call star_coeffs(sym_type,n_b_half,region_bound_single,&
         &region_bound_double,n_star_in_sphere, &
         n_cube2sphere,coeff,xsite,ysite,zsite)
    if((inode == ionode).and.(iprint_basis >= 2)) then
       write(unit=io_lun,fmt='(/" summary of info about stars:"/ &
            & " star no.",3x,"x",4x,"y",4x,"z",7x,"coeff",9x,"mult")')
       do n = 1, n_star_in_sphere
          mult = very_small + one/(coeff(n)*coeff(n))
          write(unit=io_lun,fmt='(1x,i4,3x,3i5,3x,e15.6,3x,i4)') n, &
               xsite(n), ysite(n), zsite(n), coeff(n), mult
       end do
    end if

    ! make block of blip-overlap matrix having requested symmetry
    call symmat(sym_type,n_b_half,n_star_in_sphere,n_cube2sphere, &
         xsite,ysite,zsite,coeff,amat)
    if((inode == ionode).and.(iprint_basis >= 3)) then
       write(unit=io_lun,fmt='(/" elements of block of blip-overlap", &
            & " matrix having requested symmetry:"/)')
       do m = 1, n_star_in_sphere
          do n = 1, n_star_in_sphere
             write(unit=io_lun,fmt='(3x,2i5,3x,e15.6)') m, n, amat(m,n)
          end do
       end do
    end if

    ! compute and assemble contributions to scalar product
    do i = 1, n_b_int_in_diam
       z0 = b_int*(i-1-n_b_int_in_radius)
       call f3_value(sym_type,n_sp,n_am,n_zeta,nu_int,&
            &delta_ig,n_blip,bv,z0,fv3)
       kb_min = 4 - min(n_b_int_in_diam-i,3)
       kb_max = min(i,4)
       i_store = 1 + mod(i-1,4)
       do j = kb_min, kb_max
          j_store = 1 + mod(i_store-j+4,4)
          do n1 = 1, n_blip
             do n2 = 1, n_blip
                sum = zero
                do n = 1, nu_int
                   sum = sum + bv(n,j)*fv3(n,n1,n2)
                end do
                store(i_store,j_store,n1,n2) = deltax*sum
             end do
          end do
       end do
       if(i >= 4) then
          n3 = i - 3
          j_add = 1 + mod(i_store,4)
          do n1 = 1, n_blip
             do n2 = 1, n_blip
                n = (n1-1)*n_blip_squared + (n2-1)*n_blip + n3
                scal_prod(n) = store(1,j_add,n1,n2) + store(2,j_add,n1,n2) + &
                     store(3,j_add,n1,n2) + store(4,j_add,n1,n2)
             end do
          end do
       end if
    end do

    ! make star components of scalar product
    do n = 1, n_star_in_sphere
       n1 = xsite(n) + n_b_int_in_radius - 1
       n2 = ysite(n) + n_b_int_in_radius - 1
       n3 = zsite(n) + n_b_int_in_radius - 1
       m = (n1-1)*n_blip_squared + (n2-1)*n_blip + n3
       scal_prod_sym(n) = scal_prod(m)/coeff(n)
    end do
    if((inode == ionode).and.(iprint_basis >= 2)) then
       write(unit=io_lun,fmt='(/" symmetric scalar product:"/)')
       do n = 1, n_star_in_sphere
          write(unit=io_lun,fmt='(3x,i5,3x,e15.6)') n, scal_prod_sym(n)
       end do
    end if

    ! solve linear equations to get blip coefficients
    call linsolv(amat,n_star_in_sphere,scal_prod_sym)
    ! renormalise to get blip coefficients
    do n = 1, n_star_in_sphere
       scal_prod_sym(n) = coeff(n)*scal_prod_sym(n)
    end do
    if((inode == ionode).and.(iprint_basis >= 2)) then
       write(unit=io_lun,fmt='(/" blip coefficients:"/)')
       do n = 1, n_star_in_sphere
          write(unit=io_lun,fmt='(1x,i5,3x,e15.6)') n, scal_prod_sym(n)
       end do
    end if

    return

  end subroutine blips_in_star
!!*** 

! -----------------------------------------------------------
! Subroutine deviation
! -----------------------------------------------------------

!!****f* diagnostics/deviation *
!!
!!  NAME 
!!   deviation
!!  USAGE
!! 
!!  PURPOSE
!!   Calculates a measure of the difference between the blip
!!   representation of the support function and the original
!!   PAO representation. This measure is the rms difference
!!   between the two functions, normalised by dividing
!!   by the norm of the support functions in its
!!   PAO representation.
!! 
!!  INPUTS
!!   inode: id of node running the executable
!!   ionode: id of input/output node
!!   n_sp: id number of the atomic species
!!   n_sup: id number of the support function
!!   nu_int: number of integration-grid intervals in the
!!     blip-grid interval
!!   n_b_half: with blip-grid sites identified by the triplet
!!     of integers (nx,ny,nz), with (0,0,0) being the site in
!!     centre of the support region, n_b_half is the maximum
!!     magnitude of nx, ny or nz for blip-functions contained
!!     wholly within the support region
!!   reg_rad: support-region radius for current species.
!!   b_int: blip-grid spacing for current species.
!!   region_bound_single, region_bound_double: index arrays 
!!     specifying max. cartesian indices
!!     for blip functions contained wholly within support region.
!!   (For a detailed description of these index arrays, see header
!!     to the routine "sphere".)
!!   blip_number: index array: blip_number(n1,n2,n3) gives the
!!     sequence number of a blip function. Here, n1, n2 and n3
!!     run over range 1 to 1+2*n_b_half
!!   bv: pre-calculated values of the B-spline
!!   blip: array holding values of blip coefficients 
!! 
!!  USES
!! 
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!   12/08/2002
!!    i) All pao tables changed so that enumeration starts from one not zero
!!    ii) All radial tables changed so that the expectation is that they are divided by
!!    r**l (in line with Siesta radial tables and spherical harmonic conventions)
!!  SOURCE
  subroutine deviation(inode,ionode,n_sp,n_sup,nu_int,&
       &n_b_half,reg_rad,b_int,region_bound_single,region_bound_double,&
       &blip_number,bas,bv,blip)
    use maxima_module, ONLY : max_blip_nu_int
    use datatypes
    use GenComms, ONLY: cq_abort
    use global_module, ONLY : iprint_basis
    use numbers
    use pao_format
    use read_support_spec, ONLY: support_info

    implicit none

    integer, intent(in) :: inode, ionode, n_sp, n_sup, nu_int, n_b_half,bas
    integer, intent(in), dimension(0:bas) :: region_bound_single
    integer, intent(in), dimension(0:bas,0:bas) :: region_bound_double
    integer, intent(in), dimension(-bas:bas,-bas:bas, -bas:bas) :: blip_number
    real(double), intent(in) :: reg_rad, b_int
    real(double), intent(in), dimension(:,:) :: bv
    real(double), intent(in), dimension(:,:) :: blip
    integer :: i, j, k, j0, ired, jred, kred, n, nx, ny, nz, &
         &n1, n2, n3, abs_nx, abs_ny, n_ac, n_acz, n_am, n_zeta, &
         &n_int_grid, n_blip_add, n_int_add, n_int_points, &
         &nb_beg0, nb_beg, nb_end0, nb_end
    real(double) :: a1g_norm, t1u_norm, t2g_norm, eg_a_norm, eg_b_norm
    real(double) :: acc_dphi_squared, acc_phi_squared, alpha, ang_fac, &
         &coeff, delta_ig, deltar, dphi_rms, phi_rms, phi_blip, phi_pao, &
         &r, sum, x, y, z
    real(double), dimension(2*bas+1,2*bas+1,(2*bas+11)*max_blip_nu_int) :: blip1
    real(double), dimension(2*bas+1,(2*bas+11)*max_blip_nu_int, &
         (2*bas+11)*max_blip_nu_int) :: blip2

    a1g_norm = sqrt(one/(four*pi))
    t1u_norm = sqrt(three/(four*pi))
    t2g_norm = sqrt(fifteen/(four*pi))
    eg_a_norm = sqrt(fifteen/(sixteen*pi))
    eg_b_norm = sqrt(five/(sixteen*pi))

    delta_ig = b_int/float(nu_int)
    n_int_grid = int(very_small + reg_rad/delta_ig)
    k = 2*n_int_grid + 1
    if(k > (2*bas+11)*max_blip_nu_int) then
       call cq_abort('deviation: too many integration-grid pts',&
            &k,(2*bas+11)*max_blip_nu_int)
    end if
    n_blip_add = 1 + n_int_grid/nu_int
    n_int_add = n_blip_add*nu_int

    !   blip -> grid z-transform --------------------------------------
    do nx = -n_b_half, n_b_half
       n1 = nx + 1 + n_b_half
       abs_nx = abs(nx)
       do ny = -region_bound_single(abs_nx), region_bound_single(abs_nx)
          n2 = ny + 1 + n_b_half
          abs_ny = abs(ny)
          do k = -n_int_grid, n_int_grid
             nb_beg0 = (k+n_int_add-1)/nu_int - n_blip_add - 1
             nb_end0 = nb_beg0 + 3
             nb_beg = max(-region_bound_double(abs_nx,abs_ny), nb_beg0)
             nb_end = min(region_bound_double(abs_nx,abs_ny), nb_end0)
             sum = zero
             if(nb_beg <= nb_end) then
                kred = 1 + mod(k+n_int_add-1,nu_int)
                do nz = nb_beg, nb_end
                   n3 = nz + 1 + n_b_half
                   n = blip_number(nx,ny,nz)
                   sum = sum + bv(kred,4+nb_beg0-nz)*blip(n,n_sup)
                end do
             end if
             blip1(n1,n2,k+1+n_int_grid) = sum
          end do
       end do
    end do

    !   blip -> grid y-transform --------------------------------------
    do k = -n_int_grid, n_int_grid
       do nx = -n_b_half, n_b_half
          n1 = nx + 1 + n_b_half
          abs_nx = abs(nx)
          do j = -n_int_grid, n_int_grid
             nb_beg0 = (j+n_int_add-1)/nu_int - n_blip_add - 1
             nb_end0 = nb_beg0 + 3
             nb_beg = max(-region_bound_single(abs_nx), nb_beg0)
             nb_end = min(region_bound_single(abs_nx), nb_end0)
             sum = zero
             if(nb_beg <= nb_end) then
                jred = 1 + mod(j+n_int_add-1,nu_int)
                do ny = nb_beg, nb_end
                   n2 = ny + 1 + n_b_half
                   sum = sum + bv(jred,4+nb_beg0-ny)*blip1(n1,n2,k+1+n_int_grid)
                end do
             end if
             blip2(n1,j+1+n_int_grid,k+1+n_int_grid) = sum
          end do
       end do
    end do

    ! accumulator for mean square deviation
    acc_phi_squared = zero
    acc_dphi_squared = zero
    n_int_points = 0

    !   blip -> grid z-transform
    do j = -n_int_grid, n_int_grid
       y = j*delta_ig
       do k = -n_int_grid, n_int_grid
          z = k*delta_ig
          do i = -n_int_grid, n_int_grid
             x = i*delta_ig

             ! integration grid-point must be inside region
             if(x*x + y*y + z*z <= reg_rad*reg_rad) then
                n_int_points = n_int_points + 1
                r = sqrt(x*x + y*y + z*z)

                ! calculate value of support function in PAO representation
                phi_pao = zero
                !%%!count = 1
                !%%!!   loop over angular components and over zeta's
                !%%!do n_am=0,support_info(n_sp)%lmax
                !%%!   n_ac = l+1
                !%%!   do n_zeta = 1,support_info(n_sp)%naczs(l)
                !%%!      do m=-l,l
                !%%!         coeff = support_info(n_sp)%supp_func(n_sup)%coefficients(count)
                do n_acz = 1, support_info(n_sp,n_sup)%no_of_acz
                   n_ac = support_info(n_sp,n_sup)%acz_info(n_acz)%which_ac
                   n_zeta = support_info(n_sp,n_sup)%acz_info(n_acz)%which_zeta
                   coeff = support_info(n_sp,n_sup)%acz_info(n_acz)%coeff
                   if(n_ac == 1) then
                      n_am = 0
                   else if((n_ac >= 2) .and. (n_ac <= 4)) then
                      n_am = 1
                   else if((n_ac >= 5) .and. (n_ac <= 9)) then
                      n_am = 2
                   else
                      call cq_abort('deviation: n_ac out of range',n_ac)
                   end if
                   deltar = pao(n_sp)%angmom(n_am)%zeta(n_zeta)%cutoff/&
                        &(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length - 1)
                   j0 = 1 + int(r/deltar)
                   if(j0+1 <= pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length) then
                      select case(n_ac)
                      case(1)
                         ang_fac = a1g_norm
                      case(2)
                         ang_fac = t1u_norm*x
                      case(3)
                         ang_fac = t1u_norm*y
                      case(4)
                         ang_fac = t1u_norm*z
                      case(5)
                         ang_fac = t2g_norm*y*z
                      case(6)
                         ang_fac = t2g_norm*z*x
                      case(7)
                         ang_fac = t2g_norm*x*y
                      case(8)
                         ang_fac = eg_a_norm*(x*x - y*y)
                      case(9)
                         ang_fac = eg_b_norm*(three*z*z - r*r)
                      end select
                      alpha = one + r/deltar - j0
                      phi_pao = phi_pao + &
                           &((one-alpha)*pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(j0) + &
                           &alpha*pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(j0+1))*ang_fac*coeff
                   end if
                end do
          
                acc_phi_squared = acc_phi_squared + phi_pao*phi_pao
                nb_beg0 = (i+n_int_add-1)/nu_int - n_blip_add - 1
                nb_end0 = nb_beg0 + 3
                nb_beg = max(-n_b_half, nb_beg0)
                nb_end = min(n_b_half, nb_end0)
                phi_blip = zero
                if(nb_beg <= nb_end) then
                   ired = 1 + mod(i+n_int_add-1,nu_int)
                   do nx = nb_beg, nb_end
                      n1 = nx + 1 + n_b_half
                      phi_blip = phi_blip + bv(ired,4+nb_beg0-nx)*&
                           &blip2(n1,j+1+n_int_grid,k+1+n_int_grid)
                   end do
                end if
                acc_dphi_squared = acc_dphi_squared + &
                     &(phi_blip - phi_pao)*(phi_blip - phi_pao)

             end if
             
          end do
       end do
    end do

    ! calculate r.m.s. deviation and normalise it
    phi_rms = sqrt(acc_phi_squared/float(n_int_points))
    dphi_rms = sqrt(acc_dphi_squared/float(n_int_points))
    if((inode == ionode).and.(iprint_basis >= 1)) then
       write(unit=io_lun,fmt='(/" species:",i3," support fn no:",i3,&
            &", normalised rms deviation:",e15.6)') n_sp, n_sup, &
            &dphi_rms/phi_rms
    end if

    return

  end subroutine deviation
!!*** 

end module pao2blip
