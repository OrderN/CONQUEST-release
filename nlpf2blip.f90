module nlpf2blip

  use datatypes
  use numbers
  use global_module,          only: io_lun
  use timer_stdclocks_module, only: start_timer, stop_timer, &
                                    tmr_std_basis, tmr_std_allocation
  use support_spec_format,    only: support_function

  implicit none
  save

  type proj_function
     integer :: ncoeffs
     real(double), pointer, dimension(:) :: coefficients ! Dimension npaos (sum over l of  acz(l)*(2l+1)
  end type proj_function

  type projector_function
     integer :: nsuppfuncs
     type(proj_function), pointer, dimension(:) :: proj_func ! coefficients for a sf, dimension nsuppfuncs
  end type projector_function

  !type(projector_function), allocatable, dimension(:), target :: nlpf_on_atom ! Dimension mx_atoms (flag above)
  type(support_function), allocatable, dimension(:), target :: nlpf_on_atom
  real(double), allocatable, dimension(:), target :: nlpf_array


  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id$"

contains

  subroutine make_blips_from_nlpfs

    use datatypes
    use numbers
    use blip!, only: blip_data, blip_info, blip_FFT_size, blip_FFT_off
    use dimens, only: RadiusSupport
    use GenComms, only: cq_abort, gcopy, myid, my_barrier, inode, ionode
    use global_module, only: iprint_basis
    use pseudo_tm_info
    use species_module, only: nlpf_species, n_species
    use blip_grid_transform_module, only: do_local_grid_to_blip, &
                                          do_local_blip_to_grid
    use io_module, only: grab_blip_coeffs
    
    implicit none

    integer :: i, i_stop, n, na, nb, nx, ny, nz, ns, n_ac, n_acz, &
               n_am, n_s, n_sp, n_sup, n_zeta, min_ac, max_ac, l, m, &
               count, marker, size
    integer :: n_blip, nu_int, stat, max_nlpf, the_l, nsf_send
    integer :: iblock, ipart, ia, no_of_ib_ia, pos, ip, lenx, nl, &
               mval, j, off, offx, offy, offz
    integer, allocatable, dimension(:) :: n_cube2sphere, n_b_half
    real(double) :: c, deltax, r2_over_b2, x, y, z, r, step, rr, a, b,&
                    d, r1, r2, r3, r4, nl_potential, val, xg, yg, zg, &
                    ang_fac
    real(double), allocatable, dimension(:,:) :: local_blip
    real(double), allocatable, dimension(:,:,:,:) :: local_grid
    real(double) :: a1g_norm, t1u_norm, t2g_norm, eg_a_norm, eg_b_norm

    call start_timer(tmr_std_basis)
    if((inode == ionode).and.(iprint_basis >= 2)) then
       write(unit=io_lun,fmt='(//" make_blips_from_nlpfs: sbrt entered")')
    end if
    a1g_norm = sqrt(one/(four*pi))
    t1u_norm = sqrt(three/(four*pi))
    t2g_norm = sqrt(fifteen/(four*pi))
    eg_a_norm = sqrt(fifteen/(sixteen*pi))
    eg_b_norm = sqrt(five/(sixteen*pi))
    ! number of species
    if((inode == ionode).and.(iprint_basis >= 1)) then
       write(unit=io_lun,fmt='(//" make_nlpfs_from_blips:&
            & no. of species:",i5)') n_species
    end if
    ! -------------------------------------------------------------
    !   loop over species
    ! -------------------------------------------------------------
    ! Storage for NLPF blip coefficients
    size=0
    allocate(nlpf_on_atom(n_species))
    do i = 1, n_species
       size=size+nlpf_species(i)*blip_info(i)%NBlipsRegion
       allocate(nlpf_on_atom(i)%supp_func(nlpf_species(i)))
    end do
    allocate(nlpf_array(size))
    marker=0
    do i = 1, n_species
       n_sp = i
       do n_sup = 1, nlpf_species(n_sp)
          count = blip_info(n_sp)%NBlipsRegion
          nlpf_on_atom(i)%supp_func(n_sup)%coefficients => nlpf_array(marker+1:marker+count)
          marker = marker + count          
       end do
    end do
    marker = 0
    do n_sp = 1, n_species
       xg = blip_info(n_sp)%SupportGridSpacing
       yg = blip_info(n_sp)%SupportGridSpacing
       zg = blip_info(n_sp)%SupportGridSpacing
       lenx = blip_FFT_size!2*(blip_info(n_sp)%BlipArraySize+4)+1
       allocate(local_grid(lenx,lenx,lenx,nlpf_species(n_sp)))
       local_grid = zero
       ! Create NLPF on grid
       off = blip_FFT_off+1!blip_info(n_sp)%BlipArraySize+4+1
       do nz=-blip_info(n_sp)%BlipArraySize,blip_info(n_sp)%BlipArraySize
          !if(nz<0) then
          !   offz = blip_FFT_size+1
          !else
          !   offz = 1
          !end if
          do ny=-blip_info(n_sp)%BlipArraySize,blip_info(n_sp)%BlipArraySize
             !if(ny<0) then
             !   offy = blip_FFT_size+1
             !else
             !   offy = 1
             !end if
             do nx=-blip_info(n_sp)%BlipArraySize,blip_info(n_sp)%BlipArraySize
                !if(nx<0) then
                !   offx = blip_FFT_size+1
                !else
                !   offx = 1
                !end if
                x = xg*real(nx,double)
                y = yg*real(ny,double)
                z = zg*real(nz,double)
                r = sqrt(x*x+y*y+z*z)
                if(r>very_small) then
                   x=x/r
                   y=y/r
                   z=z/r
                else
                   x = zero
                   y = zero
                   z = zero
                end if
                !write(70,*) nx,ny,nz,r,nlpf_species(n_sp),pseudo(n_sp)%pjnl(1)%cutoff
                n_sup = 0
                do nl= 1, pseudo(n_sp)%n_pjnl
                   the_l = pseudo(n_sp)%pjnl_l(nl)  
                   !do n_sup = 1,nlpf_species(n_sp)
                   !if(n_sup>4) then 
                   !   nl = 3
                   !   mval = n_sup - 7
                   !else if(n_sup>1) then
                   !   nl = 2
                   !   mval = n_sup - 3
                   !else
                   !   nl = 1
                   !   mval = 0
                   !end if
                   do mval=-the_l,the_l
                   n_sup = n_sup+1
                   if(r<pseudo(n_sp)%pjnl(nl)%cutoff) then
                      step = pseudo(n_sp)%pjnl(nl)%delta
                      j = aint(r/step)+1
                      if(j<pseudo(n_sp)%pjnl(nl)%n) then
                         rr = real(j,double)*step
                         a = ( rr - r ) / step
                         b = one - a
                         c = a * ( a * a - one ) * step * step / six
                         d = b * ( b * b - one ) * step * step / six
                         
                         r1=pseudo(n_sp)%pjnl(nl)%f(j)
                         r2=pseudo(n_sp)%pjnl(nl)%f(j+1)
                         r3=pseudo(n_sp)%pjnl(nl)%d2(j)
                         r4=pseudo(n_sp)%pjnl(nl)%d2(j+1)
                         
                         nl_potential =  a * r1 + b * r2 + c * r3 + d * r4
                         select case(the_l)
                         case(0)
                            ang_fac = a1g_norm
                         case(1)
                            select case(mval)
                            case(-1)
                               ang_fac = t1u_norm*y*r
                            case(0)
                               ang_fac = t1u_norm*z*r
                            case(1)
                               ang_fac = t1u_norm*x*r
                            end select
                         case(2)
                            select case(mval)
                            case(-2)
                               ang_fac = t2g_norm*x*y*r*r
                            case(-1)
                               ang_fac = t2g_norm*y*z*r*r
                            case(0)
                               ang_fac = eg_b_norm*(three*z*z - one)*r*r
                            case(1)
                               ang_fac = t2g_norm*z*x*r*r
                            case(2)
                               ang_fac = eg_a_norm*(x*x - y*y)*r*r
                            end select
                         end select
                         local_grid(nx+off,ny+off,nz+off,n_sup) = nl_potential*ang_fac ! val
                         !write(110+n_sup,*) nl_potential*ang_fac!nx+off,ny+off,nz+off,nl_potential,ang_fac
                      end if ! j<pseudo(n_sp)%pjnl(nl)%n
                   end if ! r<pseudo(n_sp)%pjnl(nl)%cutoff
                   end do
                end do ! n_sup
             end do ! nx
          end do ! ny
       end do ! nz
       !write(*,*) "Calling G2B, size: ",lenx
       call do_local_grid_to_blip(blip_info(n_sp), nlpf_species(n_sp), &
            nlpf_array(marker+1:marker+nlpf_species(n_sp)*blip_info(n_sp)%NBlipsRegion), &
            local_grid,lenx)
       marker = marker+nlpf_species(n_sp)*blip_info(n_sp)%NBlipsRegion
       deallocate(local_grid)
    end do ! n_sp
    !call grab_blip_coeffs(nlpf_array,size,inode,"nlpfs")
    ! -------------------------------------------------------------
    !   end loop over species
    ! -------------------------------------------------------------
  end subroutine make_blips_from_nlpfs
!!*** 

  subroutine get_SP(blip_co, nlpf_co, matSP, ip, j_in_halo, dx, dy, &
                    dz, speci, specj)

    use datatypes
    use numbers
    use GenBlas, only: axpy, copy, scal, gemm
    use blip, only: blip_info
    use support_spec_format, only: support_function
    use GenComms, only: cq_abort, mtime, inode, ionode
    use mult_module, only: store_matrix_value, scale_matrix_value, &
                           return_matrix_value_pos, &
                           store_matrix_value_pos, matrix_pos
    use species_module, only: nsf_species, nlpf_species

    implicit none

    ! Shared Variables
    integer :: this_nsf, this_nlpf, speci, specj, np, nn, ip, matSP, j_in_halo

    type(support_function) :: blip_co
    type(support_function) :: nlpf_co
    !type(projector_function) :: nlpf_co
    real(double) :: dx, dy, dz

    ! Local Variables
    integer, parameter :: MAX_D = 4

    real(double) ::  FAC(-MAX_D:MAX_D,3)

    real(double), allocatable, dimension(:) :: work1, work2, work3, &
                                               work4, work5, work6
    real(double), allocatable, dimension(:,:) :: temp
    real(double) :: tmp, facx, facy, facz, factot, t0, t1, mat_val

    integer :: ix, iy, iz, offset, l, at, nsf1, stat, i1, i2, m, x,y,&
               z, xmin, xmax, ymin, ymax, zmin, zmax, idx, idy, idz, &
               wheremat

    this_nsf = nsf_species(speci)
    this_nlpf = nlpf_species(specj)
    allocate(temp(this_nlpf,this_nsf))
    idx = floor(abs(dx/blip_info(speci)%SupportGridSpacing))
    if(dx<zero) then
       idx=-idx
       if(abs(dx-real(idx,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idx=idx-1
    else
       if(abs(dx-real(idx,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idx=idx+1
    end if
    idy = floor(abs(dy/blip_info(speci)%SupportGridSpacing))
    if(dy<zero) then
       idy=-idy
       if(abs(dy-real(idy,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idy=idy-1
    else
       if(abs(dy-real(idy,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idy=idy+1
    end if
    idz = floor(abs(dz/blip_info(speci)%SupportGridSpacing))
    if(dz<zero) then
       idz=-idz
       if(abs(dz-real(idz,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idz=idz-1
    else
       if(abs(dz-real(idz,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idz=idz+1
    end if
    call do_blip_integrals(FAC,dx/blip_info(speci)%SupportGridSpacing-real(idx,double), &
         dy/blip_info(speci)%SupportGridSpacing-real(idy,double),dz/blip_info(speci)%SupportGridSpacing-real(idz,double))
    allocate(work1(blip_info(speci)%FullArraySize*this_nsf),work2(blip_info(speci)%FullArraySize*this_nsf), &
         work4(blip_info(speci)%FullArraySize*this_nsf),work6(blip_info(speci)%FullArraySize*this_nsf), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating arrays for onsite S elements: ",blip_info(speci)%FullArraySize,this_nsf)
    ! first, we copy the blip functions for this atom onto a cubic grid;
    ! we make this grid 'too big' in order to have a fast routine below.

    work1 = zero
    offset = blip_info(speci)%BlipArraySize+1+3

    
    do ix = -blip_info(speci)%BlipArraySize, blip_info(speci)%BlipArraySize
       do iy = -blip_info(speci)%BlipArraySize, blip_info(speci)%BlipArraySize
          do iz = -blip_info(speci)%BlipArraySize, blip_info(speci)%BlipArraySize
             l = blip_info(speci)%blip_number(ix,iy,iz)
             if (l.ne.0) then
                at = (((iz+offset)*blip_info(speci)%OneArraySize + (iy+offset))*blip_info(speci)%OneArraySize + &
                     (ix+offset)) * this_nsf
                do nsf1 = 1,this_nsf
                   work1(nsf1+at) = blip_co%supp_func(nsf1)%coefficients(l)
                enddo
             end if
          end do
       end do
    end do
    
    ! now, for each direction in turn, we need to apply a number of
    ! 'spreading operations'. Do z first... put blip(z) in 2, and
    ! del2blip(z) in 3
    
    call copy(blip_info(speci)%FullArraySize*this_nsf,work1,1,work2,1)
    call scal(blip_info(speci)%FullArraySize*this_nsf,FAC(0,3),work2,1)
    do iz = 1, MAX_D
       offset = iz * blip_info(speci)%OneArraySize * blip_info(speci)%OneArraySize * this_nsf
       call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), FAC(iz,3), &
            work1(1:), 1, work2(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), FAC(-iz,3), &
            work1(1+offset:), 1, work2(1:), 1 )
    end do
    
    ! now do y : put blip(y).blip(z) in 4,
    ! blip(y).del2blip(z) + del2blip(y).blip(z)  in 5
    
    call copy(blip_info(speci)%FullArraySize*this_nsf,work2,1,work4,1)
    call scal(blip_info(speci)%FullArraySize*this_nsf,FAC(0,2),work4,1)
    do iy = 1, MAX_D
       offset = iy * blip_info(speci)%OneArraySize * this_nsf
       call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), FAC(iy,2), &
            work2(1:), 1, work4(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), FAC(-iy,2), &
            work2(1+offset:), 1, work4(1:), 1 )
    end do

    ! Now x
    work6 = zero
    call axpy(blip_info(speci)%FullArraySize*this_nsf,FAC(0,1),work4,1,work6,1)
    do ix = 1, MAX_D
       offset = ix * this_nsf
       call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), FAC(ix,1), &
            work4(1:), 1, work6(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), FAC(-ix,1), &
            work4(1+offset:), 1, work6(1:), 1 )
    end do
    
    ! and now get the matrix elements by multiplication...
    
    temp = zero
    ! Now set work1 to Projector coeffs
    deallocate(work1)
    allocate(work1(blip_info(specj)%FullArraySize*this_nlpf))
    work1 = zero
    offset = blip_info(specj)%BlipArraySize+1+3
    
    xmin = -blip_info(specj)%BlipArraySize
    xmax =  blip_info(specj)%BlipArraySize
    ymin = -blip_info(specj)%BlipArraySize
    ymax =  blip_info(specj)%BlipArraySize
    zmin = -blip_info(specj)%BlipArraySize
    zmax =  blip_info(specj)%BlipArraySize
    if(idx<0) then 
       xmin = xmin -3
    else if(idx>0) then 
       xmax = xmax +3
    end if
    if(idy<0) then 
       ymin = ymin -3
    else if(idy>0) then 
       ymax = ymax +3
    end if
    if(idz<0) then 
       zmin = zmin -3
    else if(idz>0) then 
       zmax = zmax +3
    end if
    do ix = xmin,xmax
       do iy = ymin,ymax
          do iz = zmin,zmax
             if(abs(ix-idx)<=blip_info(speci)%BlipArraySize.AND.abs(iy-idy)<=blip_info(speci)%BlipArraySize.AND. &
                  abs(iz-idz)<=blip_info(speci)%BlipArraySize) then
                l = blip_info(specj)%blip_number(ix-idx,iy-idy,iz-idz)
             else
                l = 0
             end if
             if (l.ne.0) then
                at = (((iz+offset)*blip_info(specj)%OneArraySize + (iy+offset))*blip_info(specj)%OneArraySize + &
                     (ix+offset)) * this_nlpf
                do nsf1 = 1,this_nlpf
                   work1(nsf1+at) = nlpf_co%supp_func(nsf1)%coefficients(l)
                enddo
             end if
          end do
       end do
    end do
    call gemm('n','t',this_nlpf,this_nsf,&
         blip_info(specj)%OneArraySize*blip_info(specj)%OneArraySize*blip_info(specj)%OneArraySize, &
         one,work1,this_nlpf,work6,this_nsf,zero,temp,this_nlpf )
    tmp = zero
    do i1 = 1,this_nsf
       do i2=1,this_nlpf
          mat_val = blip_info(specj)%SupportGridSpacing*blip_info(specj)%SupportGridSpacing* &
               blip_info(specj)%SupportGridSpacing*temp(i2,i1)
          !write(60,*) mat_val
          !write(io_lun,fmt='(2x,"Spread SP: ",2i4,f20.12)') i1,i2, mat_val
          wheremat = matrix_pos(matSP,ip,j_in_halo,i1,i2)
          !write(io_lun,fmt='(2x,"Grid SP  : ",2i4,f20.12)') i1,i2, return_matrix_value_pos(matSP,wheremat)
          !if(abs(return_matrix_value_pos(matSP,wheremat))>very_small) &
          !     write(30+ip,fmt='(f20.12)') mat_val - return_matrix_value_pos(matSP,wheremat)
          !if(abs(mat_val - return_matrix_value_pos(matSP,wheremat))>tmp) &
          !     tmp = abs(mat_val - return_matrix_value_pos(matSP,wheremat))
          !call scale_matrix_value(matSP,np,nn,ip,nabj,i1,i2,zero)
          call store_matrix_value_pos(matSP,wheremat,mat_val)
          !if(abs(temp(i2,i1))>1e-6_double) then
          !call scale_matrix_value(matSP,np,nn,ip,0,i1,i2,zero,1)
          !call store_matrix_value(matSP,np,nn,ip,0,i1,i2, &
          !     SupportGridSpacing(spec)*SupportGridSpacing(spec)*SupportGridSpacing(spec)*temp(i2,i1),1)
          !end if
       end do
    end do
    !write(io_lun,fmt='(2x,"Maximum SP deviation: ",f20.12)') tmp
    deallocate(work1,work2, work4,work6, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating arrays for onsite S blip elements: ",blip_info(speci)%FullArraySize,this_nsf)
    return
  end subroutine get_SP


  subroutine get_dSP(blip_co, nlpf_co, matSP, ip, j_in_halo, dx, dy, &
                     dz, speci, specj,dir)

    use datatypes
    use numbers
    use GenBlas, only: axpy, copy, scal, gemm
    use blip, only: blip_info
    use support_spec_format, only: support_function
    use GenComms, only: cq_abort, mtime, inode, ionode
    use mult_module, only: store_matrix_value, scale_matrix_value, &
                           return_matrix_value_pos, &
                           store_matrix_value_pos, matrix_pos
    use species_module, only: nsf_species, nlpf_species

    implicit none

    ! Shared Variables
    integer :: this_nsf, this_nlpf, speci, specj, np, nn, ip, matSP, &
               j_in_halo, dir

    type(support_function) :: blip_co
    type(support_function) :: nlpf_co
    !type(projector_function) :: nlpf_co
    real(double) :: dx, dy, dz

    ! Local Variables
    integer, parameter :: MAX_D = 4

    real(double) ::  FAC(-MAX_D:MAX_D,3), DFAC(-MAX_D:MAX_D,3)

    real(double), allocatable, dimension(:) :: work1, work2, work3, &
                                               work4, work5, work6
    real(double), allocatable, dimension(:,:) :: temp
    real(double) :: tmp, facx, facy, facz, factot, t0, t1, mat_val

    integer :: ix, iy, iz, offset, l, at, nsf1, stat, i1, i2, m, x,y,&
               z, xmin, xmax, ymin, ymax, zmin, zmax, idx, idy, idz, &
               wheremat

    this_nsf = nsf_species(speci)
    this_nlpf = nlpf_species(specj)
    allocate(temp(this_nlpf,this_nsf))
    !FAC(0) = 151.0_double/140.0_double
    !FAC(1) = 1191.0_double/2240.0_double
    !FAC(2) = 3.0_double/56.0_double
    !FAC(3) = 1.0_double/2240.0_double
    !DFAC(0) = zero
    !DFAC(1) = 49.0_double/64.0_double
    !DFAC(2) = 7.0_double/40.0_double
    !DFAC(3) = 1.0_double/320.0_double
    idx = floor(abs(dx/blip_info(speci)%SupportGridSpacing))
    if(dx<zero) then
       idx=-idx
       if(abs(dx - real(idx,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idx = idx-1
    else
       if(abs(dx - real(idx,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idx = idx+1
    end if
    idy = floor(abs(dy/blip_info(speci)%SupportGridSpacing))
    if(dy<zero) then
       idy=-idy
       if(abs(dy - real(idy,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idy = idy-1
    else
       if(abs(dy - real(idy,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idy = idy+1
    end if
    idz = floor(abs(dz/blip_info(speci)%SupportGridSpacing))
    if(dz<zero) then
       idz=-idz
       if(abs(dz - real(idz,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idz = idz-1
    else
       if(abs(dz - real(idz,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idz = idz+1
    end if
    call do_blip_integrals(FAC,dx/blip_info(speci)%SupportGridSpacing-real(idx,double), &
         dy/blip_info(speci)%SupportGridSpacing-real(idy,double),dz/blip_info(speci)%SupportGridSpacing-real(idz,double))
    call do_dblip_integrals(DFAC,dx/blip_info(speci)%SupportGridSpacing-real(idx,double), &
         dy/blip_info(speci)%SupportGridSpacing-real(idy,double),dz/blip_info(speci)%SupportGridSpacing-real(idz,double))
    !write(io_lun,fmt='(2x,"Offset: ",3f12.7," in BG ",3i2)') dx, dy, dz, idx, idy, idz
    allocate(work1(blip_info(speci)%FullArraySize*this_nsf),work2(blip_info(speci)%FullArraySize*this_nsf), &
         work4(blip_info(speci)%FullArraySize*this_nsf),work6(blip_info(speci)%FullArraySize*this_nsf), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating arrays for onsite S elements: ",blip_info(speci)%FullArraySize,this_nsf)
    ! first, we copy the blip functions for this atom onto a cubic grid;
    ! we make this grid 'too big' in order to have a fast routine below.

    work1 = zero
    offset = blip_info(speci)%BlipArraySize+1+3

    
    do ix = -blip_info(speci)%BlipArraySize, blip_info(speci)%BlipArraySize
       do iy = -blip_info(speci)%BlipArraySize, blip_info(speci)%BlipArraySize
          do iz = -blip_info(speci)%BlipArraySize, blip_info(speci)%BlipArraySize
             l = blip_info(speci)%blip_number(ix,iy,iz)
             if (l.ne.0) then
                at = (((iz+offset)*blip_info(speci)%OneArraySize + (iy+offset))*blip_info(speci)%OneArraySize + &
                     (ix+offset)) * this_nsf
                do nsf1 = 1,this_nsf
                   work1(nsf1+at) = blip_co%supp_func(nsf1)%coefficients(l)
                enddo
             end if
          end do
       end do
    end do
    
    ! now, for each direction in turn, we need to apply a number of
    ! 'spreading operations'. Do z first... put blip(z) in 2, and
    ! del2blip(z) in 3
    if(dir==3) then
       call copy(blip_info(speci)%FullArraySize*this_nsf,work1,1,work2,1)
       call scal(blip_info(speci)%FullArraySize*this_nsf,DFAC(0,3),work2,1)
       do iz = 1, MAX_D
          offset = iz * blip_info(speci)%OneArraySize * blip_info(speci)%OneArraySize * this_nsf
          call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), DFAC(iz,3), &
               work1(1:), 1, work2(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), DFAC(-iz,3), &
               work1(1+offset:), 1, work2(1:), 1 )
       end do
    else
       call copy(blip_info(speci)%FullArraySize*this_nsf,work1,1,work2,1)
       call scal(blip_info(speci)%FullArraySize*this_nsf,FAC(0,3),work2,1)
       do iz = 1, MAX_D
          offset = iz * blip_info(speci)%OneArraySize * blip_info(speci)%OneArraySize * this_nsf
          call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), FAC(iz,3), &
               work1(1:), 1, work2(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), FAC(-iz,3), &
               work1(1+offset:), 1, work2(1:), 1 )
       end do
    end if

    ! now do y : put blip(y).blip(z) in 4,
    ! blip(y).del2blip(z) + del2blip(y).blip(z)  in 5
    if(dir==2) then
       call copy(blip_info(speci)%FullArraySize*this_nsf,work2,1,work4,1)
       call scal(blip_info(speci)%FullArraySize*this_nsf,DFAC(0,2),work4,1)
       do iy = 1, MAX_D
          offset = iy * blip_info(speci)%OneArraySize * this_nsf
          call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), DFAC(iy,2), &
               work2(1:), 1, work4(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsf-offset),  DFAC(-iy,2), &
               work2(1+offset:), 1, work4(1:), 1 )
       end do
    else
       call copy(blip_info(speci)%FullArraySize*this_nsf,work2,1,work4,1)
       call scal(blip_info(speci)%FullArraySize*this_nsf,FAC(0,2),work4,1)
       do iy = 1, MAX_D
          offset = iy * blip_info(speci)%OneArraySize * this_nsf
          call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), FAC(iy,2), &
               work2(1:), 1, work4(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), FAC(-iy,2), &
               work2(1+offset:), 1, work4(1:), 1 )
       end do
    end if

    ! x
    work6 = zero
    if(dir==1) then
       call axpy(blip_info(speci)%FullArraySize*this_nsf,DFAC(0,1),work4,1,work6,1)
       do ix = 1, MAX_D
          offset = ix * this_nsf
          call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), DFAC(ix,1), &
               work4(1:), 1, work6(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsf-offset),  DFAC(-ix,1), &
               work4(1+offset:), 1, work6(1:), 1 )
       end do
    else
       call axpy(blip_info(speci)%FullArraySize*this_nsf,FAC(0,2),work4,1,work6,1)
       do ix = 1, MAX_D
          offset = ix * this_nsf
          call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), FAC(ix,1), &
               work4(1:), 1, work6(1+offset:), 1 )
          call axpy((blip_info(speci)%FullArraySize*this_nsf-offset), FAC(-ix,1), &
               work4(1+offset:), 1, work6(1:), 1 )
       end do
    end if
    ! and now get the matrix elements by multiplication...
    
    temp = zero
    ! Now set work1 to Projector coeffs
    deallocate(work1)
    allocate(work1(blip_info(specj)%FullArraySize*this_nlpf))
    work1 = zero
    offset = blip_info(specj)%BlipArraySize+1+3
    
    xmin = -blip_info(specj)%BlipArraySize
    xmax =  blip_info(specj)%BlipArraySize
    ymin = -blip_info(specj)%BlipArraySize
    ymax =  blip_info(specj)%BlipArraySize
    zmin = -blip_info(specj)%BlipArraySize
    zmax =  blip_info(specj)%BlipArraySize
    if(idx<0) then 
       xmin = xmin -3
    else if(idx>0) then 
       xmax = xmax +3
    end if
    if(idy<0) then 
       ymin = ymin -3
    else if(idy>0) then 
       ymax = ymax +3
    end if
    if(idz<0) then 
       zmin = zmin -3
    else if(idz>0) then 
       zmax = zmax +3
    end if
    do ix = xmin,xmax
       do iy = ymin,ymax
          do iz = zmin,zmax
             if(abs(ix-idx)<=blip_info(speci)%BlipArraySize.AND.abs(iy-idy)<=blip_info(speci)%BlipArraySize.AND. &
                  abs(iz-idz)<=blip_info(speci)%BlipArraySize) then
                l = blip_info(specj)%blip_number(ix-idx,iy-idy,iz-idz)
             else
                l = 0
             endif
             if (l.ne.0) then
                at = (((iz+offset)*blip_info(specj)%OneArraySize + (iy+offset))*blip_info(specj)%OneArraySize + &
                     (ix+offset)) * this_nlpf
                do nsf1 = 1,this_nlpf
                   work1(nsf1+at) = nlpf_co%supp_func(nsf1)%coefficients(l)
                enddo
             end if
          end do
       end do
    end do
    call gemm('n','t',this_nlpf,this_nsf,&
         blip_info(specj)%OneArraySize*blip_info(specj)%OneArraySize*blip_info(specj)%OneArraySize, &
         one,work1,this_nlpf,work6,this_nsf,zero,temp,this_nlpf )
    tmp = zero
    do i1 = 1,this_nsf
       do i2=1,this_nlpf
          mat_val = blip_info(specj)%SupportGridSpacing*blip_info(specj)%SupportGridSpacing*temp(i2,i1)
          wheremat = matrix_pos(matSP,ip,j_in_halo,i1,i2)
          !write(io_lun,fmt='(2x,"Spread dSP: ",2i4,2f20.12)') i1,i2, mat_val, return_matrix_value_pos(matSP,wheremat)
          !write(io_lun,fmt='(2x,"Grid dSP  : ",2i4,f20.12)') i1,i2, return_matrix_value_pos(matSP,wheremat)
          !if(abs(return_matrix_value_pos(matSP,wheremat))>very_small) &
          !     write(30+ip,fmt='(f20.12)') mat_val - return_matrix_value_pos(matSP,wheremat)
          !if(abs(mat_val - return_matrix_value_pos(matSP,wheremat))>tmp) &
          !     tmp = abs(mat_val - return_matrix_value_pos(matSP,wheremat))
          call store_matrix_value_pos(matSP,wheremat,mat_val)
          !if(abs(temp(i2,i1))>1e-6_double) then
          !call scale_matrix_value(matSP,np,nn,ip,0,i1,i2,zero,1)
          !call store_matrix_value(matSP,np,nn,ip,0,i1,i2, &
          !     SupportGridSpacing(spec)*SupportGridSpacing(spec)*SupportGridSpacing(spec)*temp(i2,i1),1)
          !end if
       end do
    end do
    !write(io_lun,fmt='(2x,"Maximum SP deviation: ",f20.12)') tmp
    deallocate(work1,work2, work4,work6, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating arrays for onsite S blip elements: ",blip_info(specj)%FullArraySize,this_nsf)
    return
  end subroutine get_dSP


  !!****f* 
  !! PURPOSE
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   David Bowler
  !! CREATION DATE 
  !!   2012/02/27
  !! MODIFICATION HISTORY
  !! SOURCE
  !!  
  subroutine get_blipP(blip_grad, nlpf_co, dataU, ip, j_in_halo, dx, &
                       dy, dz, speci, specj)

    use datatypes
    use numbers
    use GenBlas,             only: axpy, copy, scal, gemm
    use blip,                only: blip_info
    use support_spec_format, only: support_function
    use GenComms,            only: cq_abort, mtime, inode, ionode
    use mult_module,         only: store_matrix_value,      &
                                   scale_matrix_value,      &
                                   return_matrix_value_pos, &
                                   store_matrix_value_pos, matrix_pos
    use species_module,      only: nsf_species, nlpf_species

    implicit none

    ! Shared Variables
    integer :: this_nsf, this_nlpf, speci, specj, np, nn, ip, matSP, j_in_halo

    type(support_function) :: blip_grad
    type(support_function) :: nlpf_co
    real(double), dimension(:,:) :: dataU
    !type(projector_function) :: nlpf_co
    real(double) :: dx, dy, dz

    ! Local Variables
    integer, parameter :: MAX_D = 4

    real(double) ::  FAC(-MAX_D:MAX_D,3)

    real(double), allocatable, dimension(:) :: work1, work2, work3, &
                                               work4, work5, work6
    real(double), allocatable, dimension(:,:) :: temp
    real(double) :: tmp, facx, facy, facz, factot, t0, t1, mat_val

    integer :: ix, iy, iz, offset, l, at, nsf1, stat, i1, i2, m, x,y,&
               z, xmin, xmax, ymin, ymax, zmin, zmax, idx, idy, idz, &
               wheremat, nlpf1

    this_nsf = nsf_species(speci)
    this_nlpf = nlpf_species(specj)
    allocate(temp(this_nlpf,this_nsf))
    !FAC(0) = 151.0_double/140.0_double
    !FAC(1) = 1191.0_double/2240.0_double
    !FAC(2) = 3.0_double/56.0_double
    !FAC(3) = 1.0_double/2240.0_double
    idx = floor(abs(dx/blip_info(speci)%SupportGridSpacing))
    if(dx<zero) then
       idx=-idx
       if(abs(dx - real(idx,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idx = idx-1
    else
       if(abs(dx - real(idx,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idx = idx+1
    end if
    idy = floor(abs(dy/blip_info(speci)%SupportGridSpacing))
    if(dy<zero) then
       idy=-idy
       if(abs(dy - real(idy,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idy = idy-1
    else
       if(abs(dy - real(idy,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idy = idy+1
    end if
    idz = floor(abs(dz/blip_info(speci)%SupportGridSpacing))
    if(dz<zero) then
       idz=-idz
       if(abs(dz - real(idz,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idz = idz-1
    else
       if(abs(dz - real(idz,double)*blip_info(speci)%SupportGridSpacing)>half*blip_info(speci)%SupportGridSpacing) idz = idz+1
    end if
    !write(io_lun,fmt='(2x,"Offset: ",3f12.7," in BG ",3i2)') dx, dy, dz, idx, idy, idz
    call do_blip_integrals(FAC,dx/blip_info(speci)%SupportGridSpacing-real(idx,double), &
         dy/blip_info(speci)%SupportGridSpacing-real(idy,double),dz/blip_info(speci)%SupportGridSpacing-real(idz,double))
    allocate(work1(blip_info(speci)%FullArraySize*this_nlpf),work2(blip_info(speci)%FullArraySize*this_nlpf), &
         work4(blip_info(speci)%FullArraySize*this_nlpf),work6(blip_info(speci)%FullArraySize*this_nlpf), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating arrays for onsite S elements: ",blip_info(speci)%FullArraySize,this_nlpf)
    ! first, we copy the blip functions for this atom onto a cubic grid;
    ! we make this grid 'too big' in order to have a fast routine below.

    work1 = zero
    offset = blip_info(speci)%BlipArraySize+1+3

    
    xmin = -blip_info(specj)%BlipArraySize
    xmax =  blip_info(specj)%BlipArraySize
    ymin = -blip_info(specj)%BlipArraySize
    ymax =  blip_info(specj)%BlipArraySize
    zmin = -blip_info(specj)%BlipArraySize
    zmax =  blip_info(specj)%BlipArraySize
    if(idx<0) then 
       xmin = xmin -3
    else if(idx>0) then 
       xmax = xmax +3
    end if
    if(idy<0) then 
       ymin = ymin -3
    else if(idy>0) then 
       ymax = ymax +3
    end if
    if(idz<0) then 
       zmin = zmin -3
    else if(idz>0) then 
       zmax = zmax +3
    end if
    do ix = xmin,xmax
       do iy = ymin,ymax
          do iz = zmin,zmax
             if(abs(ix-idx)<=blip_info(speci)%BlipArraySize.AND.abs(iy-idy)<=blip_info(speci)%BlipArraySize.AND. &
                  abs(iz-idz)<=blip_info(speci)%BlipArraySize) then
                l = blip_info(specj)%blip_number(ix-idx,iy-idy,iz-idz)
             else
                l=0
             end if
             !l = blip_info(speci)%blip_number(ix,iy,iz)
             if (l.ne.0) then
                at = (((iz+offset)*blip_info(speci)%OneArraySize + (iy+offset))*blip_info(speci)%OneArraySize + &
                     (ix+offset)) * this_nlpf
                do nlpf1 = 1,this_nlpf
                   !do nsf1 =1,this_nsf
                   !   work1(nlpf1+at) + &
                   !        dataU(nsf1,nlpf1)*nlpf_co%supp_func(nlpf1)%coefficients(l)
                   !end do
                   work1(nlpf1+at) = nlpf_co%supp_func(nlpf1)%coefficients(l)
                enddo
             end if
          end do
       end do
    end do
    
    ! now, for each direction in turn, we need to apply a number of
    ! 'spreading operations'. Do z first... put blip(z) in 2, and
    ! del2blip(z) in 3
    
    call copy(blip_info(speci)%FullArraySize*this_nlpf,work1,1,work2,1)
    call scal(blip_info(speci)%FullArraySize*this_nlpf,FAC(0,3),work2,1)
    do iz = 1, MAX_D
       offset = iz * blip_info(speci)%OneArraySize * blip_info(speci)%OneArraySize * this_nlpf
       call axpy((blip_info(speci)%FullArraySize*this_nlpf-offset), FAC(-iz,3), &
            work1(1:), 1, work2(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nlpf-offset), FAC(iz,3), &
            work1(1+offset:), 1, work2(1:), 1 )
    end do
    
    ! now do y : put blip(y).blip(z) in 4,
    ! blip(y).del2blip(z) + del2blip(y).blip(z)  in 5
    
    call copy(blip_info(speci)%FullArraySize*this_nlpf,work2,1,work4,1)
    call scal(blip_info(speci)%FullArraySize*this_nlpf,FAC(0,2),work4,1)
    do iy = 1, MAX_D
       offset = iy * blip_info(speci)%OneArraySize * this_nlpf
       call axpy((blip_info(speci)%FullArraySize*this_nlpf-offset), FAC(-iy,2), &
            work2(1:), 1, work4(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nlpf-offset), FAC(iy,2), &
            work2(1+offset:), 1, work4(1:), 1 )
    end do
    
    work6 = zero
    call axpy(blip_info(speci)%FullArraySize*this_nlpf,FAC(0,1),work4,1,work6,1)
    do ix = 1, MAX_D
       offset = ix * this_nlpf
       call axpy((blip_info(speci)%FullArraySize*this_nlpf-offset), FAC(-ix,1), &
            work4(1:), 1, work6(1+offset:), 1 )
       call axpy((blip_info(speci)%FullArraySize*this_nlpf-offset), FAC(ix,1), &
            work4(1+offset:), 1, work6(1:), 1 )
    end do
    
    ! and now get the matrix elements by multiplication...
    
    call scal(blip_info(specj)%FullArraySize*this_nlpf, &
         blip_info(specj)%SupportGridSpacing*blip_info(specj)%SupportGridSpacing*blip_info(specj)%SupportGridSpacing,&
         work6,1)
    offset = blip_info(specj)%BlipArraySize+1+3
    
    xmin =  -blip_info(speci)%BlipArraySize!max(-blip_info(speci)%BlipArraySize+abs(idx),-blip_info(specj)%BlipArraySize)
    xmax =  blip_info(speci)%BlipArraySize!min( blip_info(speci)%BlipArraySize-abs(idx), blip_info(specj)%BlipArraySize)
    ymin =  -blip_info(speci)%BlipArraySize!max(-blip_info(speci)%BlipArraySize+abs(idy),-blip_info(specj)%BlipArraySize)
    ymax =  blip_info(speci)%BlipArraySize!min( blip_info(speci)%BlipArraySize-abs(idy), blip_info(specj)%BlipArraySize)
    zmin =  -blip_info(speci)%BlipArraySize!max(-blip_info(speci)%BlipArraySize+abs(idz),-blip_info(specj)%BlipArraySize)
    zmax =  blip_info(speci)%BlipArraySize!min( blip_info(speci)%BlipArraySize-abs(idz), blip_info(specj)%BlipArraySize)
    !write(55,*) 'Limits: ',xmin,xmax,ymin,ymax,zmin,zmax,idx,idy,idz
    do ix = xmin,xmax
       do iy = ymin,ymax
          do iz = zmin,zmax
             !write(55,*) ix,iy,iz,ix-idx,iy-idy,iz-idz
             !l = blip_info(specj)%blip_number(ix-idx,iy-idy,iz-idz)
             m = blip_info(specj)%blip_number(ix,iy,iz)
             if (m/=0) then
                at = (((iz+offset)*blip_info(specj)%OneArraySize + (iy+offset))*blip_info(specj)%OneArraySize + &
                     (ix+offset)) * this_nlpf
                do nlpf1 =1,this_nlpf
                   !write(34,*) work6(nlpf1+at)
                   do nsf1 = 1,this_nsf
                      blip_grad%supp_func(nsf1)%coefficients(m) = &
                           blip_grad%supp_func(nsf1)%coefficients(m) &
                           + dataU(nsf1,nlpf1)*work6(nlpf1+at)
                   end do
                enddo
             end if
          end do
       end do
    end do
    deallocate(work1,work2, work4,work6, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating arrays for onsite S blip elements: ",blip_info(specj)%FullArraySize,this_nsf)
    !allocate(work1(this_nlpf*NBlipsRegion(speci)))
    !do l=1,NBlipsRegion(speci)
    !   x = blip_info(speci)%blip_location(1,l)
    !   y = blip_info(speci)%blip_location(2,l)
    !   z = blip_info(speci)%blip_location(3,l)
    !   xmin = max(-blip_info(speci)%BlipArraySize,x-3)
    !   xmax = min( blip_info(speci)%BlipArraySize,x+3)
    !   ymin = max(-blip_info(speci)%BlipArraySize,y-3)
    !   ymax = min( blip_info(speci)%BlipArraySize,y+3)
    !   zmin = max(-blip_info(speci)%BlipArraySize,z-3)
    !   zmax = min( blip_info(speci)%BlipArraySize,z+3)
    !   do ix = xmin,xmax
    !      facx = FAC(abs(x-ix))
    !      do iy = ymin,ymax
    !         facy = FAC(abs(y-iy))
    !         do iz = zmin,zmax
    !            facz = FAC(abs(z-iz))
    !            if((abs(ix-idx)<=blip_info(speci)%BlipArraySize).AND.(abs(iy-idy)<=blip_info(speci)%BlipArraySize) &
    !                 .AND.(abs(iz-idz)<=blip_info(speci)%BlipArraySize)) then
    !               m = blip_info(speci)%blip_number(ix-idx,iy-idy,iz-idz)
    !            else
    !               m=0
    !            end if
    !            if(m/=0) then
    !               factot = facx*facy*facz
    !               do i2=1,this_nlpf
    !                  work1(m+(i2-1)*NBlipsRegion(speci)) = &
    !                       work1(m+(i2-1)*NBlipsRegion(speci)) + factot*nlpf_co%supp_func(1)%coefficients(l)
    !               end do
    !            end if
    !         end do
    !      end do
    !   end do
    !end do
    !do l=1,this_nsf
    !   do m=1,NBlipsRegion(speci)
    !      factot=zero
    !      do i1=1,this_nlpf
    !         factot=factot+dataU(l,i1)*work1(m+(i1-1)*NBlipsRegion(speci))
    !      end do
    !      write(71,*) factot*blip_info(specj)%SupportGridSpacing*blip_info(specj)%SupportGridSpacing*blip_info(specj)%SupportGridSpacing, &
    !           blip_grad%supp_func(l)%coefficients(m)
    !   end do
    !end do
    !call flush(71)
    !t1 = mtime()
    !if(inode==ionode) write(io_lun,*) 'Loop time: ',t1-t0
    !do i1 = 1,this_nsf
    !   do i2=1,this_nlpf
    !      tmp = return_matrix_value(matSP,np,nn,ip,0,i1,i2,1)
    !      write(io_lun,fmt='(2x,"Grid SP: ",2i4,f20.12)') i1,i2,tmp
    !      write(io_lun,fmt='(2x,"Blip SP: ",2i4,f20.12)') i1,i2, &
    !           SupportGridSpacing(spec)*SupportGridSpacing(spec)*SupportGridSpacing(spec)*temp(i2,i1)
    !      call scale_matrix_value(matSP,np,nn,ip,0,i1,i2,zero,1)
    !      call store_matrix_value(matSP,np,nn,ip,0,i1,i2, &
    !           SupportGridSpacing(spec)*SupportGridSpacing(spec)*SupportGridSpacing(spec)*temp(i2,i1),1)
    !   end do
    !end do

    return
  end subroutine get_blipP
  !!*****

  subroutine do_blip_integrals(FAC, dx, dy, dz)

    use datatypes
    use numbers

    real(double), dimension(-4:4,3) :: FAC
    real(double) :: dx, dy, dz, delta, sum, x, xmin, xmax

    integer :: i, n, nbg

    delta = 0.001_double
    ! x 
    do nbg = -4,4
       xmin = max(-two,-two+dx+real(nbg,double))
       xmax = min(two,two+dx+real(nbg,double))
       n = floor((xmax-xmin)/delta) + 1
       sum = zero
       do i=1,n
          x = xmin + real(i-1,double)*delta
          sum = sum + blipf(x)*blipf(x-dx-real(nbg,double))
       end do
       FAC(nbg,1) = delta*sum
    end do
    ! y
    do nbg = -4,4
       xmin = max(-two,-two+dy+real(nbg,double))
       xmax = min(two,two+dy+real(nbg,double))
       n = floor((xmax-xmin)/delta) + 1
       sum = zero
       do i=1,n
          x = xmin + real(i-1,double)*delta
          sum = sum + blipf(x)*blipf(x-dy-real(nbg,double))
       end do
       FAC(nbg,2) = delta*sum
    end do
    ! z
    do nbg = -4,4
       xmin = max(-two,-two+dz+real(nbg,double))
       xmax = min(two,two+dz+real(nbg,double))
       n = floor((xmax-xmin)/delta) + 1
       sum = zero
       do i=1,n
          x = xmin + real(i-1,double)*delta
          sum = sum + blipf(x)*blipf(x-dz-real(nbg,double))
       end do
       FAC(nbg,3) = delta*sum
    end do
    return
  end subroutine do_blip_integrals


  subroutine do_dblip_integrals(FAC, dx, dy, dz)

    use datatypes
    use numbers

    real(double), dimension(-4:4,3) :: FAC
    real(double) :: dx, dy, dz, delta, sum, x, xmin, xmax

    integer :: i, n, nbg

    delta = 0.001_double
    ! x 
    do nbg = -4,4
       xmin = max(-two,-two+dx+real(nbg,double))
       xmax = min(two,two+dx+real(nbg,double))
       n = floor((xmax-xmin)/delta) + 1
       sum = zero
       do i=1,n
          x = xmin + real(i-1,double)*delta
          sum = sum + dblipf(x)*blipf(x-dx-real(nbg,double))
       end do
       FAC(nbg,1) = delta*sum
    end do
    ! y
    do nbg = -4,4
       xmin = max(-two,-two+dy+real(nbg,double))
       xmax = min(two,two+dy+real(nbg,double))
       n = floor((xmax-xmin)/delta) + 1
       sum = zero
       do i=1,n
          x = xmin + real(i-1,double)*delta
          sum = sum + dblipf(x)*blipf(x-dy-real(nbg,double))
       end do
       FAC(nbg,2) = delta*sum
    end do
    ! z
    do nbg = -4,4
       xmin = max(-two,-two+dz+real(nbg,double))
       xmax = min(two,two+dz+real(nbg,double))
       n = floor((xmax-xmin)/delta) + 1
       sum = zero
       do i=1,n
          x = xmin + real(i-1,double)*delta
          sum = sum + dblipf(x)*blipf(x-dz-real(nbg,double))
       end do
       FAC(nbg,3) = delta*sum
    end do
    return
  end subroutine do_dblip_integrals


  subroutine do_d2blip_integrals(FAC, dx, dy, dz)

    use datatypes
    use numbers

    real(double), dimension(-4:4,3) :: FAC
    real(double) :: dx, dy, dz, delta, sum, x, xmin, xmax

    integer :: i, n, nbg

    delta = 0.001_double
    ! x 
    do nbg = -4,4
       xmin = max(-two,-two+dx+real(nbg,double))
       xmax = min(two,two+dx+real(nbg,double))
       n = floor((xmax-xmin)/delta) + 1
       sum = zero
       do i=1,n
          x = xmin + real(i-1,double)*delta
          sum = sum + d2blipf(x)*blipf(x-dx-real(nbg,double))
       end do
       FAC(nbg,1) = delta*sum
    end do
    ! y
    do nbg = -4,4
       xmin = max(-two,-two+dy+real(nbg,double))
       xmax = min(two,two+dy+real(nbg,double))
       n = floor((xmax-xmin)/delta) + 1
       sum = zero
       do i=1,n
          x = xmin + real(i-1,double)*delta
          sum = sum + d2blipf(x)*blipf(x-dy-real(nbg,double))
       end do
       FAC(nbg,2) = delta*sum
    end do
    ! z
    do nbg = -4,4
       xmin = max(-two,-two+dz+real(nbg,double))
       xmax = min(two,two+dz+real(nbg,double))
       n = floor((xmax-xmin)/delta) + 1
       sum = zero
       do i=1,n
          x = xmin + real(i-1,double)*delta
          sum = sum + d2blipf(x)*blipf(x-dz-real(nbg,double))
       end do
       FAC(nbg,3) = delta*sum
    end do
    return
  end subroutine do_d2blip_integrals


! NB The minus sign on delta*sum is there because we should really differentiate
! the blip on the atom but we can't because this becomes discontinuous
  subroutine do_d3blip_integrals(FAC, dx, dy, dz)

    use datatypes
    use numbers

    real(double), dimension(-4:4,3) :: FAC
    real(double) :: dx, dy, dz, delta, sum, x, xmin, xmax

    integer :: i, n, nbg

    delta = 0.001_double
    ! x 
    do nbg = -4,4
       xmin = max(-two,-two+dx+real(nbg,double))
       xmax = min(two,two+dx+real(nbg,double))
       n = floor((xmax-xmin)/delta) + 1
       sum = zero
       do i=1,n
          x = xmin + real(i-1,double)*delta
          sum = sum + d2blipf(x)*dblipf(x-dx-real(nbg,double))
       end do
       FAC(nbg,1) = -delta*sum
    end do
    ! y
    do nbg = -4,4
       xmin = max(-two,-two+dy+real(nbg,double))
       xmax = min(two,two+dy+real(nbg,double))
       n = floor((xmax-xmin)/delta) + 1
       sum = zero
       do i=1,n
          x = xmin + real(i-1,double)*delta
          sum = sum + d2blipf(x)*dblipf(x-dy-real(nbg,double))
       end do
       FAC(nbg,2) = -delta*sum
    end do
    ! z
    do nbg = -4,4
       xmin = max(-two,-two+dz+real(nbg,double))
       xmax = min(two,two+dz+real(nbg,double))
       n = floor((xmax-xmin)/delta) + 1
       sum = zero
       do i=1,n
          x = xmin + real(i-1,double)*delta
          sum = sum + d2blipf(x)*dblipf(x-dz-real(nbg,double))
       end do
       FAC(nbg,3) = -delta*sum
    end do
    return
  end subroutine do_d3blip_integrals


  real(double) function blipf(x)

    use datatypes
    use numbers

    implicit none

    real(double) :: x, y

    y = x
    if(x<0) y = -y
    if(y>two) then
       blipf = zero
    else if(y<one) then
       blipf = one-three_halves*y*y+three_quarters*y*y*y
    else
       blipf = quarter*(two-y)*(two-y)*(two-y)
    end if
  end function blipf


  real(double) function dblipf(x)

    use datatypes
    use numbers

    implicit none

    real(double) :: x, y

    y = x
    if(x<0) y = -y
    if(y>two) then
       dblipf = zero
    else if(y<one) then
       dblipf = -three*y+nine_quarters*y*y
    else
       dblipf = -three_quarters*(two-y)*(two-y)
    end if
    if(x<0) dblipf = -dblipf
  end function dblipf


  real(double) function d2blipf(x)

    use datatypes
    use numbers

    implicit none

    real(double) :: x, y

    y = x
    if(x<0) y = -y
    if(y>two) then
       d2blipf = zero
    else if(y<one) then
       d2blipf = three*(three_halves*y-one) ! -three+three*three_halves*y
    else
       d2blipf = three_halves*(two-y)
    end if
  end function d2blipf

end module nlpf2blip
