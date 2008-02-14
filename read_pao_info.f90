! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module read_pao_info
! ------------------------------------------------------------------------------
! Code area 1: Initialisation
! ------------------------------------------------------------------------------

!!****h* Conquest/read_pao_info *
!!  NAME
!!   read_pao_info
!!  PURPOSE
!!   Sole purpose is to contain the sbrt read_pao
!!  USES
!!   common, datatypes, fdf, GenComms, global_module, pao_format
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!   12/8/2002 mjg
!!    Changed to do dynamic allocation of memory while reading is in progress
!!   14:04, 2003/12/19 dave & rc
!!    Changed reading to follow new derived types: occupation numbers outside table
!!   2008/02/06 08:09 dave
!!    Changed for output to file not stdout
!!  SOURCE
module read_pao_info

  use global_module, ONLY: io_lun

  implicit none

  save

  character(len=80) :: pao_info_file
  integer :: pao_norm_flag


  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id$"

!!***

contains

! -----------------------------------------------------------
! Subroutine read_pao
! -----------------------------------------------------------

!!****f* read_pao_info/read_pao *
!!
!!  NAME 
!!   read_pao
!!  USAGE
!!   read_pao(inode,ionode,n_species)
!!  PURPOSE
!!   Reads all data for pseudo-atomic orbitals (PAO's). For the
!!   format used to read and store the date, see the header
!!   of the module pao_format.
!!  INPUTS
!!   inode: id of node running the executable
!!   ionode: id of input/output node
!!   n_species: number of atomic species in simulated system
!!  USES
!!   common, datatypes, fdf, GenComms, global_module, pao_format
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!   12/08/2002 mjg
!!    i) Changed enumeration of tables to start from one (not zero)
!!    ii) Changed to use dynamic allocation of memory and allocate during input
!!   15/08/2002 mjg
!!    added code to do optional normalisation of PAO's
!!   11:55, 25/09/2002 mjg & drb 
!!    Changed default normalisation option to 1 - i.e. normalisation of PAOs is done by default
!!   08:19, 2004/07/27 dave
!!    Added lines to work out maximum value for npao
!!  SOURCE
!!
  subroutine read_pao(inode,ionode,n_species)
    use datatypes
    use GenComms, ONLY: cq_abort, gcopy
    use global_module, ONLY : iprint_init, area_init
    use numbers, ONLY : zero, one, two, three, four
    use spline_module, ONLY : spline
    use pao_format
    use memory_module, ONLY: reg_alloc_mem, type_dbl

    implicit none

    integer, intent(in) :: inode, ionode, n_species
    integer :: alls, ios, limit, lun, n_sp, n_am, n_zeta, nt
    integer :: count
    real(double) :: deltar, norm_factor, r, sum, yp1, ypn

    if(inode == ionode) then

       if(iprint_init >= 0) then
          write(unit=io_lun,fmt='(//1x,60("+")/10x,"read_pao_info: no. of species for",&
               &" reading pao data:",i5/1x,60("+"))') n_species
       end if
       if(n_species < 1) then
          call cq_abort('read_pao: no. of species must be positive',&
               &n_species)
       end if

       call io_assign(lun)
       if(iprint_init >= 0) then
          write(unit=io_lun,fmt='(/" read_pao_info: io_assign unit no. lun:",&
               &i3)') lun
       end if
       if(iprint_init >= 0) then
          write(unit=io_lun,fmt='(/" read_pao_info: name of input file:",a80)') &
               &pao_info_file
       end if

       open(unit=lun,file=pao_info_file,status='old',iostat=ios)
       if(ios /= 0) call cq_abort('read_pao: failed to open input file' )

    end if

    allocate(pao(n_species), stat = alls)
    if(alls /= 0) call cq_abort('read_pao: failed to allocate pao' )

    do n_sp = 1, n_species

       if(inode == ionode) then      
          read(unit=lun,fmt=*) pao(n_sp)%greatest_angmom
          if(iprint_init >= 0) then
             write(unit=io_lun,fmt='(//" species number:",i5/1x,60("-"))') n_sp
             write(unit=io_lun,fmt='(" greatest angular momentum:",i3)') &
                  &pao(n_sp)%greatest_angmom
          end if
          if(pao(n_sp)%greatest_angmom < 0) &
               &call cq_abort('read_pao: greatest ang. mom. cannot be negative',pao(n_sp)%greatest_angmom)
       end if

       call gcopy(pao(n_sp)%greatest_angmom)

       allocate(pao(n_sp)%angmom(0:pao(n_sp)%greatest_angmom), stat = alls)
       if(alls /= 0) call cq_abort('read_pao: failed to allocate pao%angmom' )

       count = 0
       do n_am = 0, pao(n_sp)%greatest_angmom

          if(inode == ionode) then
             read(unit=lun, fmt=*) pao(n_sp)%angmom(n_am)%n_zeta_in_angmom
             if(iprint_init >= 0) then
                write(unit=io_lun,fmt='(//" angular momentum:",i5/1x,60("."))') n_am
                write(unit=io_lun,fmt='(/" no. of zetas in this angular momentum:",i3)') &
                     &pao(n_sp)%angmom(n_am)%n_zeta_in_angmom
             end if
          end if

          call gcopy(pao(n_sp)%angmom(n_am)%n_zeta_in_angmom)

          if(pao(n_sp)%angmom(n_am)%n_zeta_in_angmom > 0) then
             allocate(pao(n_sp)%angmom(n_am)%zeta(pao(n_sp)%angmom(n_am)%n_zeta_in_angmom), stat = alls)
             if(alls /= 0) call cq_abort('read_pao: failed to allocate pao%angmom%zeta' )
             !RC adding occ pointer here
             allocate(pao(n_sp)%angmom(n_am)%occ(pao(n_sp)%angmom(n_am)%n_zeta_in_angmom), stat = alls)
             if(alls /= 0) call cq_abort('read_pao: failed to allocate pao%angmom%occ' )

             do n_zeta = 1, pao(n_sp)%angmom(n_am)%n_zeta_in_angmom

                count = count + 2*n_am + 1
                if(inode == ionode) then
                   read(unit=lun,fmt=*) pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length
                   read(unit=lun,fmt=*) pao(n_sp)%angmom(n_am)%zeta(n_zeta)%cutoff
                   read(unit=lun,fmt=*) pao(n_sp)%angmom(n_am)%occ(n_zeta)
                   if(iprint_init >= 0) then
                      write(unit=io_lun,fmt='(/" zeta no:",i5)') n_zeta
                      write(unit=io_lun,fmt='(" table length for this ang. mom. and zeta no:",i5)') &
                           pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length
                      write(unit=io_lun,fmt='(" table cut-off radius for this ang. mom. and zeta no:",f12.6)') &
                           pao(n_sp)%angmom(n_am)%zeta(n_zeta)%cutoff
                      write(unit=io_lun,fmt='(" occupancy for this ang. mom. and zeta no:",f12.6)') &
                           pao(n_sp)%angmom(n_am)%occ(n_zeta)
                   end if
                end if

                call gcopy(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length)
                call gcopy(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%cutoff)
                call gcopy(pao(n_sp)%angmom(n_am)%occ(n_zeta))

                if(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length >= 1) then
                   allocate(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(&
                        &pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length), stat = alls)
                   if(alls /= 0) call cq_abort('read_pao: failed to allocate pao%angmom%zeta%table' )
                   call reg_alloc_mem(area_init, pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length, type_dbl)
                   if(inode == ionode) then
                      do nt = 1, pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length
                         read(unit=lun,fmt=*) r, pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(nt)
                      end do
                   end if
                   
                   call gcopy(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table,&
                        &pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length)
                   
                end if

             end do ! n_zeta_in_angmom
          end if
       end do ! greatest_angmom
!       if(count > npao) npao = count
    end do ! n_species
    if(inode == ionode) call io_close(lun)
    
    !RC loop to create tables of second derivatives
    do n_sp = 1, n_species
       do n_am = 0, pao(n_sp)%greatest_angmom
          do n_zeta = 1, pao(n_sp)%angmom(n_am)%n_zeta_in_angmom

             deltar = pao(n_sp)%angmom(n_am)%zeta(n_zeta)%cutoff / &
                  &float(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length - 1)
             limit = pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length

             allocate(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table2(limit),STAT=alls)
             if(alls /= 0) call cq_abort('read_pao: failed to allocate pao%angmom%zeta%table2' )
             call reg_alloc_mem(area_init, limit, type_dbl)

             yp1 = (pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(2)-&
                  &pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(1))/deltar
             ypn = (pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(limit)-&
                  &pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(limit-1))/deltar
             

             call spline(limit,deltar,pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table,&
                  yp1,ypn,pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table2)             
          enddo
       enddo
    enddo
    if(iprint_init>3.AND.inode==ionode) write(io_lun,*) 'finished making tables of PAO second derivatives'
    
    ! Optional normalisation of PAO's
    if((inode == ionode).and.(iprint_init >= 1)) &
         &write(unit=io_lun,fmt='(/" read_pao: pao_norm_flag:",i3)') pao_norm_flag

    if(pao_norm_flag == 1) then

       do n_sp = 1, n_species
          do n_am = 0, pao(n_sp)%greatest_angmom
             if(pao(n_sp)%angmom(n_am)%n_zeta_in_angmom > 0) then
                do n_zeta = 1, pao(n_sp)%angmom(n_am)%n_zeta_in_angmom
                   if(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length < 3) then
                      if((inode == ionode).and.(iprint_init > 0)) then
                         write(unit=io_lun,fmt='(//" read_pao: n_sp:",i3," n_am:",i3," n_zeta:",i3,":"/&
                              &" PAO normalisation requested, but cannot be done because table length < 3")') &
                              &n_sp, n_am, n_zeta
                      end if
                   else
                      deltar = pao(n_sp)%angmom(n_am)%zeta(n_zeta)%cutoff / &
                           &float(pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length - 1)
                      sum = zero
                      limit = (pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length - 1) / 2
                      limit = 1 + 2*limit
                      do nt = 2, limit-1, 2
                         r = (nt-1)*deltar
                         sum = sum + four*((r**(1+n_am))*pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(nt))**2
                      end do
                      if(limit-2 >= 3) then
                         do nt = 3, limit-2, 2
                            r = (nt-1)*deltar
                            sum = sum + two*((r**(1+n_am))*pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(nt))**2
                         end do
                      end if
                      r = (limit-1)*deltar
                      sum = sum + ((r**(1+n_am))*pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(nt))**2
                      sum = sum*deltar/three
                      norm_factor = sqrt(sum)
                      do nt = 1, pao(n_sp)%angmom(n_am)%zeta(n_zeta)%length
                         pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(nt) = norm_factor*&
                              &pao(n_sp)%angmom(n_am)%zeta(n_zeta)%table(nt)
                      end do
                   end if
                end do
             end if
          end do
       end do

       if((inode == ionode).and.(iprint_init >= 1)) &
            &write(unit=io_lun,fmt='(/" read_pao: normalisation of PAOs done")')

    end if

    return

  end subroutine read_pao
!!*** 

end module read_pao_info
