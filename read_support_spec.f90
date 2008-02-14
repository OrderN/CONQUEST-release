! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module read_support_spec
! ------------------------------------------------------------------------------
! Code area 11: Basis operations
! ------------------------------------------------------------------------------

!!****h* Conquest/read_support_spec *
!!  NAME
!!   read_support_spec
!!  PURPOSE
!!   Sole purpose of the module is to contain the sbrt read_support
!!  USES
!!
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!   2006/06/12 08:10 dave
!!    Added per-species PAO coefficient information derived type
!!   2008/02/06 08:10 dave
!!    Changed for output to file not stdout
!!  SOURCE
module read_support_spec

  use datatypes
  use global_module, ONLY: io_lun

  implicit none

  save

  type acz_spec
     integer :: which_ac
     integer :: which_zeta
     real(double) :: coeff
  end type acz_spec

  type support_spec
     integer :: no_of_acz
     type(acz_spec), pointer, dimension(:) :: acz_info
  end type support_spec

  type(support_spec), allocatable, dimension(:,:) :: support_info

  character(len=80) :: support_spec_file
  logical :: flag_read_support_spec
    
  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id$"

!!***
  
contains

! -----------------------------------------------------------
! Subroutine read_support
! -----------------------------------------------------------

!!****f* read_support_spec/read_support *
!!
!!  NAME 
!!   read_support
!!  USAGE
!! 
!!  PURPOSE
!!   Reads all data specifying the initial support functions of
!!   all atomic species in terms of pseudo-atomic orbitals.
!! 
!!  INPUTS
!!   inode: id of node running the executable
!!   ionode: id of input/output node
!!   n_species: number of atomic species in the simulated system
!! 
!!  USES
!! 
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!  22/6/02
!!  MODIFICATION HISTORY
!!
!!  SOURCE
  subroutine read_support(inode,ionode,n_species)

    use datatypes
    use GenComms, ONLY: cq_abort, gcopy, my_barrier
    use global_module, ONLY : iprint_basis
    use species_module, ONLY: nsf_species, npao_species
    use maxima_module, ONLY: maxnsf
    use numbers, ONLY: one

    implicit none

    integer, intent(in) :: inode, ionode, n_species
    integer :: ios, lun, n_acz, n_sp, n_s, info, l, n_coeff, i, maxacz
    integer, allocatable, dimension(:) :: ibuff1, ibuff2 ! n_species_max*NSF*max_acz_in_spec
    real(double), allocatable, dimension(:) :: rbuff

    allocate(support_info(n_species,maxnsf),STAT=info)
    if(info/=0) call cq_abort("read_support_spec: support info: ",n_species)
    if(inode == ionode) then
       if(iprint_basis >= 3) then
          write(unit=io_lun,fmt='(//1x,60("+")/10x,"no of species for reading",&
               &" support spec:",i5/1x,60("+"))') n_species
       end if
       call io_assign(lun)
       if(iprint_basis >= 3) then
          write(unit=io_lun,fmt='(/" read_support_spec: input file unit no:",&
               &i3)') lun
       end if
       if(iprint_basis >= 3) then
          write(unit=io_lun,fmt='(/" read_support_spec: name of input file:",&
               &a80)') support_spec_file
       end if
       if(flag_read_support_spec) then
          open(unit=lun,file=support_spec_file,status='old',iostat=ios)
          if(ios /= 0) then
             call cq_abort('read_support_spec: input file not found')
          end if
          do n_sp = 1, n_species
             if(iprint_basis >= 1) then
                write(unit=io_lun,fmt='(/1x,60("*")/5x,"species no",i5/1x,60("*")/)') n_sp
                write(unit=io_lun,fmt='(/" no. of support functions:",i5)') nsf_species(n_sp)
             end if

             do n_s = 1, nsf_species(n_sp)
                if(iprint_basis >= 1) then
                   write(unit=io_lun,fmt='(/1x,40("-")/8x," support function no:",i5/&
                        &1x,40("-"))') n_s
                end if

                read(unit=lun,fmt=*) support_info(n_sp,n_s)%no_of_acz
                if(iprint_basis >= 1) then
                   write(unit=io_lun,fmt='(" no. of acz to specify support&
                        & function:",i5)') support_info(n_sp,n_s)%no_of_acz
                end if
                allocate(support_info(n_sp,n_s)%acz_info(support_info(n_sp,n_s)%no_of_acz))
                do n_acz = 1, support_info(n_sp,n_s)%no_of_acz
                   read(unit=lun,fmt=*) &
                        support_info(n_sp,n_s)%acz_info(n_acz)%which_ac, &
                        support_info(n_sp,n_s)%acz_info(n_acz)%which_zeta, &
                        support_info(n_sp,n_s)%acz_info(n_acz)%coeff
                   if(iprint_basis >= 1) then
                      write(unit=io_lun,fmt='(" acz no:",i5," which ac:",i5,& 
                           &" which zeta:",i5," coeff:",e15.6)') n_acz, &
                           &support_info(n_sp,n_s)%acz_info(n_acz)%which_ac, &
                           &support_info(n_sp,n_s)%acz_info(n_acz)%which_zeta, &
                           &support_info(n_sp,n_s)%acz_info(n_acz)%coeff
                   end if
                end do ! n_acz
             end do ! n_s = nsf_species
          end do ! n_sp = n_species
          call io_close(lun)
       else
          do n_sp = 1, n_species
             if(iprint_basis >= 1) then
                write(unit=io_lun,fmt='(/1x,60("*")/5x,"species no",i5/1x,60("*")/)') n_sp
                write(unit=io_lun,fmt='(/" no. of support functions:",i5)') nsf_species(n_sp)
             end if
             if(npao_species(n_sp)/=nsf_species(n_sp)) then
                write(io_lun,fmt='(10x,"If there are more PAOs than SFs you must specify how they map using support.dat")')
                call cq_abort("PAO-to-blip map not specified.")
             end if
             do n_s = 1, nsf_species(n_sp)
                if(iprint_basis >= 1) then
                   write(unit=io_lun,fmt='(/1x,40("-")/8x," support function no:",i5/&
                        &1x,40("-"))') n_s
                end if

                support_info(n_sp,n_s)%no_of_acz = 1
                if(iprint_basis >= 1) then
                   write(unit=io_lun,fmt='(" no. of acz to specify support&
                        & function:",i5)') support_info(n_sp,n_s)%no_of_acz
                end if
                allocate(support_info(n_sp,n_s)%acz_info(support_info(n_sp,n_s)%no_of_acz))
                support_info(n_sp,n_s)%acz_info(1)%which_ac = n_s
                support_info(n_sp,n_s)%acz_info(1)%which_zeta = 1
                support_info(n_sp,n_s)%acz_info(1)%coeff = one
                if(iprint_basis >= 1) then
                   write(unit=io_lun,fmt='(" acz no:",i5," which ac:",i5,& 
                        &" which zeta:",i5," coeff:",e15.6)') 1, &
                        &support_info(n_sp,n_s)%acz_info(1)%which_ac, &
                        &support_info(n_sp,n_s)%acz_info(1)%which_zeta, &
                        &support_info(n_sp,n_s)%acz_info(1)%coeff
                end if
             end do ! n_s = nsf_species
          end do ! n_sp = n_species
       end if
    end if ! inode==ionode

    ! distribute to all nodes all the data read above
    allocate(ibuff1(n_species*maxnsf))
    if(inode == ionode) then
       do n_sp = 1, n_species
          do n_s = 1, nsf_species(n_sp)
             ibuff1((n_sp-1)*maxnsf + n_s) = support_info(n_sp,n_s)%no_of_acz
          end do
       end do
    end if

    call gcopy(ibuff1,n_species*maxnsf)

    maxacz = 0
    do n_sp = 1, n_species
       do n_s = 1, nsf_species(n_sp)
          support_info(n_sp,n_s)%no_of_acz = ibuff1((n_sp-1)*maxnsf + n_s)
          if(support_info(n_sp,n_s)%no_of_acz>maxacz) maxacz = support_info(n_sp,n_s)%no_of_acz
          if(inode/=ionode) allocate(support_info(n_sp,n_s)%acz_info(support_info(n_sp,n_s)%no_of_acz))
       end do
    end do
    deallocate(ibuff1)
    if((inode == ionode).and.(iprint_basis >= 2)) then
       write(unit=io_lun,fmt='(/" read_support_spec:&
            & support_info(n_sp,n_s)%no_of_acz distributed")')
    end if

    allocate(ibuff1(n_species*maxnsf*maxacz),ibuff2(n_species*maxnsf*maxacz),rbuff(n_species*maxnsf*maxacz))
    if(inode == ionode) then
       do n_sp = 1, n_species
          do n_s = 1, nsf_species(n_sp)
             do n_acz = 1, support_info(n_sp,n_s)%no_of_acz
                ibuff1((n_sp-1)*maxnsf*maxacz +&
                     &(n_s-1)*maxacz + n_acz) = &
                     &support_info(n_sp,n_s)%acz_info(n_acz)%which_ac
                ibuff2((n_sp-1)*maxnsf*maxacz +&
                     &(n_s-1)*maxacz + n_acz) = &
                     &support_info(n_sp,n_s)%acz_info(n_acz)%which_zeta
                rbuff((n_sp-1)*maxnsf*maxacz +&
                     &(n_s-1)*maxacz + n_acz) = &
                     &support_info(n_sp,n_s)%acz_info(n_acz)%coeff
             end do
          end do
       end do
    end if

    call gcopy(ibuff1,n_species*maxnsf*maxacz)
    call gcopy(ibuff2,n_species*maxnsf*maxacz)
    call gcopy(rbuff,n_species*maxnsf*maxacz)

    do n_sp = 1, n_species
       do n_s = 1, nsf_species(n_sp)
          do n_acz = 1, support_info(n_sp,n_s)%no_of_acz
             support_info(n_sp,n_s)%acz_info(n_acz)%which_ac = &
                  &ibuff1((n_sp-1)*maxnsf*maxacz +&
                  &(n_s-1)*maxacz + n_acz)
             support_info(n_sp,n_s)%acz_info(n_acz)%which_zeta = &
                  &ibuff2((n_sp-1)*maxnsf*maxacz +&
                  &(n_s-1)*maxacz + n_acz)
             support_info(n_sp,n_s)%acz_info(n_acz)%coeff = &
                  &rbuff((n_sp-1)*maxnsf*maxacz +&
                  &(n_s-1)*maxacz + n_acz)
          end do
       end do
    end do
    deallocate(rbuff,ibuff2,ibuff1)
    if((inode == ionode).and.(iprint_basis >= 2)) then
       write(unit=io_lun,fmt='(/" read_support_spec: support_info(n_sp,n_s)%&
            &acz_info(n_acz)%which_ac distributed")')
       write(unit=io_lun,fmt='(/" read_support_spec: support_info(n_sp,n_s)%&
            &acz_info(n_acz)%which_zeta distributed")')
       write(unit=io_lun,fmt='(/" read_support_spec: support_info(n_sp,n_s)%&
            &acz_info(n_acz)%coeff distributed")')
    end if
    if(inode==ionode) then
       close(unit=lun)
    end if

    !do n_sp = n_species,1,-1
    !   do n_s = nsf_species(n_sp),1,-1
    !      deallocate(support_info(n_sp,n_s)%acz_info)
    !   end do
    !end do
    !deallocate(support_info,STAT=info)
    call my_barrier()

  end subroutine read_support
!!***
  
end module read_support_spec
