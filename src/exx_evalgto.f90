module exx_evalgto

  use datatypes
  use numbers,                ONLY: zero, one, two, three, four, five, six, fifteen, sixteen
  use numbers,                ONLY: half, quarter, very_small, pi, twopi, fourpi
  implicit none

contains


  subroutine exx_gto_on_grid_new_new(inode,atom,spec,extent,xyz,nsuppfuncs,phi_on_grid,r_int,rst)

    use timer_module,           ONLY: start_timer, stop_timer
    use exx_types,      ONLY: unit_matrix_write, tmr_std_exx_evalgto
    use exx_types,      ONLY: exx_overlap, exx_cartesian, exx_gto
    use GenComms,               ONLY: cq_abort

    use angular_coeff_routines, ONLY: evaluate_pao
    use dimens,                 ONLY: r_h   
    use gto_format_new,         ONLY: gto

    implicit none

    ! << Passed variables >>
    integer, intent(in)      :: extent   ! number of grid points [for 1 dimension]
    integer, intent(in)      :: inode
    integer, intent(in)      :: atom     
    integer, intent(in)      :: spec     
    integer, intent(in)      :: nsuppfuncs 

    real(double), intent(in) :: xyz(3)
    real(double), intent(in) :: rst(3)
    real(double), intent(in) :: r_int

    real(double), dimension(2*extent+1,2*extent+1,2*extent+1,nsuppfuncs), &
         intent(inout) :: phi_on_grid

    real(double),  parameter :: a1g_norm  = sqrt(one/(four*pi))
    real(double),  parameter :: t1u_norm  = sqrt(three/(four*pi))

    ! << Local variables >>
    !real(double), dimension(nsuppfuncs) :: sfsum
    integer             :: nx, ny, nz 
    integer             :: mx, my, mz   
    integer             :: px, py, pz  
    integer             :: ex, ey, ez  

    real(double)        :: grid_spacing 
    real(double)        :: x, y, z, r
    real(double)        :: int, n, rest
    real(double)        :: xyz_delta(3)

    integer             :: count1, ierr, stat
    integer             :: i, j, p, l1, m1, acz1, ngrid
    integer             :: ijk_delta(3), ijk(3), kji(3), njk(3), nji(3)

    real(double)        :: gto_rad, gto_cart, gto_val
    real(double)        :: a, d, factor

    ngrid        = 2*extent+1
    grid_spacing = r_int/real(extent,double)                   
    phi_on_grid  = zero

    mx = -extent ; px = extent
    my = mx      ; py = px 
    mz = my      ; pz = py

    call start_timer(tmr_std_exx_evalgto)    
    ! Define the overlap box
    overlap_box: if (exx_overlap) then
       do i = 1, 3
          n    = xyz(i)/grid_spacing
          rest = n - real(aint(n))

          if (abs(rest) >= half) then
             int = ceiling(abs(n))
             int = sign(int,n)
             xyz_delta(i) = int - n

          else if (abs(rest) < half) then
             int = floor(abs(n))
             int = sign(int,n)
             xyz_delta(i) = int - n

          else if ((rest == zero)) then
             int = n
             xyz_delta(i) = zero

          end if

          ijk_delta(i) = int
          xyz_delta(i) = xyz_delta(i)*grid_spacing 
       end do

       do i = 1, 3       
          if (ijk_delta(i) > 0) then
             ijk(i) = ijk_delta(i) + 1                                 
             njk(i) = ngrid 
             kji(i) = 1
             nji(i) = ngrid - ijk_delta(i)

          else if (ijk_delta(i) < 0) then
             ijk(i) = 1
             njk(i) = ngrid + ijk_delta(i)
             kji(i) = -ijk_delta(i) + 1
             nji(i) = ngrid 

          else
             ijk(i) = 1
             njk(i) = ngrid
             kji(i) = 1
             nji(i) = ngrid

          end if
       end do
       mx = mx +kji(1)-1
       my = my +kji(2)-1
       mz = mz +kji(3)-1
       px = px -ijk(1)+1
       py = py -ijk(2)+1
       pz = pz -ijk(3)+1
    end if overlap_box
    
    do i = 1, gto( spec )%nsf_tot

       ex = gto( spec )%sf( i )%nx
       ey = gto( spec )%sf( i )%ny
       ez = gto( spec )%sf( i )%nz
       
       grid_x_loop: do nx = mx, px
          x = xyz(1) + real(nx,double)*grid_spacing + rst(1)
          
          grid_y_loop: do ny = my, py
             y = xyz(2) + real(ny,double)*grid_spacing + rst(2)
             
             grid_z_loop: do nz = mz, pz
                z = xyz(3) + real(nz,double)*grid_spacing + rst(3)
                
                r = sqrt(x*x+y*y+z*z)  
                
                gto_cart = (x**real(ex))*(y**real(ey))*(z**real(ez))
                
                gto_rad = zero
                primitive: do p = 1, gto( spec )%sf( i )%ngto

                   a = gto( spec )%sf( i )%a( p )

                   gto_rad = gto_rad + exp(-a*r**2)
                   !gto_rad = gto_rad + d*exp(-a*r**2)
                   
                end do primitive
                
                gto_val = gto_rad * gto_cart !* factor 
                phi_on_grid(nx+extent+1,ny+extent+1,nz+extent+1,i) = gto_val
                
             end do grid_z_loop
          end do grid_y_loop
       end do grid_x_loop

    end do
    call stop_timer(tmr_std_exx_evalgto,.true.)

    return 
  end subroutine exx_gto_on_grid_new_new


  
  subroutine exx_gto_on_grid_new(inode,atom,spec,extent,xyz,nsuppfuncs,phi_on_grid,r_int,rst)

    use timer_module,           ONLY: start_timer, stop_timer
    use exx_types,      ONLY: unit_matrix_write, tmr_std_exx_evalgto
    use exx_types,      ONLY: exx_overlap, exx_cartesian, exx_gto
    use GenComms,               ONLY: cq_abort

    use angular_coeff_routines, ONLY: evaluate_pao
    use dimens,                 ONLY: r_h   
    use gto_format_new,         ONLY: gto

    implicit none

    ! << Passed variables >>
    integer, intent(in)      :: extent   ! number of grid points [for 1 dimension]
    integer, intent(in)      :: inode
    integer, intent(in)      :: atom     
    integer, intent(in)      :: spec     
    integer, intent(in)      :: nsuppfuncs 

    real(double), intent(in) :: xyz(3)
    real(double), intent(in) :: rst(3)
    real(double), intent(in) :: r_int

    real(double), dimension(2*extent+1,2*extent+1,2*extent+1,nsuppfuncs), &
         intent(inout) :: phi_on_grid

    real(double),  parameter :: a1g_norm  = sqrt(one/(four*pi))
    real(double),  parameter :: t1u_norm  = sqrt(three/(four*pi))

    ! << Local variables >>
    !real(double), dimension(nsuppfuncs) :: sfsum
    integer             :: nx, ny, nz 
    integer             :: mx, my, mz   
    integer             :: px, py, pz  
    integer             :: ex, ey, ez  

    real(double)        :: grid_spacing 
    real(double)        :: x, y, z, r
    real(double)        :: int, n, rest
    real(double)        :: xyz_delta(3)

    integer             :: count1, ierr, stat
    integer             :: i, j, p, l1, m1, acz1, ngrid
    integer             :: ijk_delta(3), ijk(3), kji(3), njk(3), nji(3)

    real(double)        :: gto_rad, gto_cart, gto_val
    real(double)        :: a, d, factor

    ngrid        = 2*extent+1
    grid_spacing = r_int/real(extent,double)                   
    phi_on_grid  = zero

    mx = -extent ; px = extent
    my = mx      ; py = px 
    mz = my      ; pz = py

    call start_timer(tmr_std_exx_evalgto)    
    ! Define the overlap box
    overlap_box: if (exx_overlap) then
       do i = 1, 3
          n    = xyz(i)/grid_spacing
          rest = n - real(aint(n))

          if (abs(rest) >= half) then
             int = ceiling(abs(n))
             int = sign(int,n)
             xyz_delta(i) = int - n

          else if (abs(rest) < half) then
             int = floor(abs(n))
             int = sign(int,n)
             xyz_delta(i) = int - n

          else if ((rest == zero)) then
             int = n
             xyz_delta(i) = zero

          end if

          ijk_delta(i) = int
          xyz_delta(i) = xyz_delta(i)*grid_spacing 
       end do

       do i = 1, 3       
          if (ijk_delta(i) > 0) then
             ijk(i) = ijk_delta(i) + 1                                 
             njk(i) = ngrid 
             kji(i) = 1
             nji(i) = ngrid - ijk_delta(i)

          else if (ijk_delta(i) < 0) then
             ijk(i) = 1
             njk(i) = ngrid + ijk_delta(i)
             kji(i) = -ijk_delta(i) + 1
             nji(i) = ngrid 

          else
             ijk(i) = 1
             njk(i) = ngrid
             kji(i) = 1
             nji(i) = ngrid

          end if
       end do
       mx = mx +kji(1)-1
       my = my +kji(2)-1
       mz = mz +kji(3)-1
       px = px -ijk(1)+1
       py = py -ijk(2)+1
       pz = pz -ijk(3)+1
    end if overlap_box
    
    count1 = 1
    angu_loop: do l1 = 0, gto(spec)%greatest_angmom

       if (l1 == 0) factor = a1g_norm
       if (l1 == 1) factor = t1u_norm 

       zeta_loop: do acz1 = 1, gto(spec)%angmom(l1)%n_zeta_in_angmom

          magn_loop: do m1 = 1, gto(spec)%angmom(l1)%nf 

             ex = gto(spec)%angmom(l1)%nx(m1)
             ey = gto(spec)%angmom(l1)%ny(m1)
             ez = gto(spec)%angmom(l1)%nz(m1)
             
             grid_x_loop: do nx = mx, px
                x = xyz(1) + real(nx,double)*grid_spacing + rst(1)

                grid_y_loop: do ny = my, py
                   y = xyz(2) + real(ny,double)*grid_spacing + rst(2)

                   grid_z_loop: do nz = mz, pz
                      z = xyz(3) + real(nz,double)*grid_spacing + rst(3)

                      r = sqrt(x*x+y*y+z*z)  

                      gto_cart = (x**real(ex))*(y**real(ey))*(z**real(ez))

                      gto_rad = zero
                      primitive: do p = 1, gto(spec)%angmom(l1)%zeta(acz1)%ngto 

                         a = gto(spec)%angmom(l1)%zeta(acz1)%a(p)
                         d = gto(spec)%angmom(l1)%zeta(acz1)%d(p)
                         
                         !gto_rad = gto_rad + exp(-a*r**2)
                         gto_rad = gto_rad + d*exp(-a*r**2)

                      end do primitive

                      gto_val = gto_rad * gto_cart * factor 
                      phi_on_grid(nx+extent+1,ny+extent+1,nz+extent+1,count1) = gto_val

                   end do grid_z_loop
                end do grid_y_loop
             end do grid_x_loop

             count1 = count1 + 1
          end do magn_loop
       end do zeta_loop
    end do angu_loop

call stop_timer(tmr_std_exx_evalgto,.true.)

return 
end subroutine exx_gto_on_grid_new

  
  subroutine exx_gto_on_grid(inode,atom,spec,extent,xyz,nsuppfuncs,phi_on_grid,r_int,rst)

    use timer_module,           ONLY: start_timer, stop_timer
    use exx_types,      ONLY: unit_matrix_write, tmr_std_exx_evalgto
    use exx_types,      ONLY: exx_overlap, exx_cartesian, exx_gto
    use GenComms,               ONLY: cq_abort

    use angular_coeff_routines, ONLY: evaluate_pao
    use dimens,                 ONLY: r_h   
    use gto_format,             ONLY: gto

    implicit none

    ! << Passed variables >>
    integer, intent(in)      :: extent   ! number of grid points [for 1 dimension]
    integer, intent(in)      :: inode
    integer, intent(in)      :: atom     
    integer, intent(in)      :: spec     
    integer, intent(in)      :: nsuppfuncs 

    real(double), intent(in) :: xyz(3)
    real(double), intent(in) :: rst(3)
    real(double), intent(in) :: r_int

    real(double), dimension(nsuppfuncs,2*extent+1,2*extent+1,2*extent+1), &
         intent(inout) :: phi_on_grid

    ! << Local variables >>
    !real(double), dimension(nsuppfuncs) :: sfsum
    integer             :: nx, ny, nz 
    integer             :: mx, my, mz   
    integer             :: px, py, pz  
    integer             :: ex, ey, ez  

    real(double)        :: grid_spacing 
    real(double)        :: x, y, z, r
    real(double)        :: int, n, rest
    real(double)        :: xyz_delta(3)

    integer             :: count1, nsf1
    integer             :: ierr, stat
    integer             :: i, j, k, kk, ngrid   
    integer             :: ijk_delta(3), ijk(3), kji(3), njk(3), nji(3)

    integer             :: gto_lmn, index
    real(double)        :: gto_val, gto_sph, gto_rad, gto_cart
    real(double)        :: ak, dk, Nk
    character(len=16)   :: nt
    
    ngrid        = 2*extent+1
    grid_spacing = r_int/real(extent,double)                   
    phi_on_grid  = zero

    mx = -extent ; px = extent
    my = mx      ; py = px 
    mz = my      ; pz = py

    call start_timer(tmr_std_exx_evalgto)    
    ! Define the overlap box
    overlap_box: if (exx_overlap) then
       do i = 1, 3
          n    = xyz(i)/grid_spacing
          rest = n - real(aint(n))

          if (abs(rest) >= half) then
             int = ceiling(abs(n))
             int = sign(int,n)
             xyz_delta(i) = int - n

          else if (abs(rest) < half) then
             int = floor(abs(n))
             int = sign(int,n)
             xyz_delta(i) = int - n

          else if ((rest == zero)) then
             int = n
             xyz_delta(i) = zero

          end if

          ijk_delta(i) = int
          xyz_delta(i) = xyz_delta(i)*grid_spacing 
       end do

       do i = 1, 3       
          if (ijk_delta(i) > 0) then
             ijk(i) = ijk_delta(i) + 1                                 
             njk(i) = ngrid 
             kji(i) = 1
             nji(i) = ngrid - ijk_delta(i)

          else if (ijk_delta(i) < 0) then
             ijk(i) = 1
             njk(i) = ngrid + ijk_delta(i)
             kji(i) = -ijk_delta(i) + 1
             nji(i) = ngrid 

          else
             ijk(i) = 1
             njk(i) = ngrid
             kji(i) = 1
             nji(i) = ngrid

          end if
       end do
       mx = mx +kji(1)-1
       my = my +kji(2)-1
       mz = mz +kji(3)-1
       px = px -ijk(1)+1
       py = py -ijk(2)+1
       pz = pz -ijk(3)+1
    end if overlap_box
!!$    
    count1 = 1
    !print*,
    shell_loop: do nsf1 = 1, gto( spec )%nshell
       
       do index = 1, gto( spec )%shell( nsf1 )%nf
          
          ex = gto( spec )%shell( nsf1 )%nx( index )
          ey = gto( spec )%shell( nsf1 )%ny( index )
          ez = gto( spec )%shell( nsf1 )%nz( index )
          nt = gto( spec )%shell( nsf1 )%nt( index )

          !print*, gto( spec )%label, 'ex ey ez', ex, ey, ez, trim(nt), &
          !     (gto( spec )%shell( nsf1 )%ak( k ), k = 1,  gto( spec )%shell( nsf1 )%n)
          !
          grid_x_loop: do nx = mx, px
             x =  xyz(1) + real(nx,double)*grid_spacing + rst(1)
             
             grid_y_loop: do ny = my, py
                y = xyz(2) + real(ny,double)*grid_spacing + rst(2)
                
                grid_z_loop: do nz = mz, pz
                   z = xyz(3) + real(nz,double)*grid_spacing + rst(3)
                   
                   r = sqrt(x*x+y*y+z*z)                                
                   
                   gto_cart = (x**real(ex))*(y**real(ey))*(z**real(ez))
                   
                   gto_rad  = zero
                   do k = 1, gto( spec )%shell( nsf1 )%n                                      
                      ak = gto( spec )%shell( nsf1 )%ak( k )
                      dk = gto( spec )%shell( nsf1 )%dk( k )
                      Nk = gto( spec )%shell( nsf1 )%Nk( k, index )

                      gto_rad = gto_rad + exp(-ak*r**2)!*Nk 

                   end do
                   
                   
                   gto_val = gto_rad * gto_cart ! / ( sqrt(4*pi) ) 

                   phi_on_grid(count1,&
                        nx+extent+1,  &
                        ny+extent+1,  &
                        nz+extent+1) = gto_val
                   
                end do grid_z_loop
             end do grid_y_loop
          end do grid_x_loop

          count1 = count1 + 1
       end do
    end do shell_loop
    
    call stop_timer(tmr_std_exx_evalgto,.true.)

    return 
  end subroutine exx_gto_on_grid
  !
  subroutine exx_gto_on_grid_prim(inode,atom,spec,extent,xyz,phi_on_grid,r_int,rst,ex,ey,ez,ak)

    use timer_module,           ONLY: start_timer, stop_timer
    use exx_types,      ONLY: unit_matrix_write, tmr_std_exx_evalgto
    use exx_types,      ONLY: exx_overlap, exx_cartesian, exx_gto
    use GenComms,               ONLY: cq_abort 

    use angular_coeff_routines, ONLY: evaluate_pao
    use dimens,                 ONLY: r_h   
    use gto_format,             ONLY: gto

    implicit none

    ! << Passed variables >>
    integer, intent(in)      :: extent   ! number of grid points [for 1 dimension]
    integer, intent(in)      :: inode
    integer, intent(in)      :: atom     
    integer, intent(in)      :: spec     
    !integer, intent(in)      :: nsuppfuncs 

    real(double), intent(in) :: xyz(3)
    real(double), intent(in) :: rst(3)
    real(double), intent(in) :: r_int

    real(double), dimension(2*extent+1,2*extent+1,2*extent+1), &
         intent(inout) :: phi_on_grid

    ! << Local variables >>
    !real(double), dimension(nsuppfuncs) :: sfsum
    integer             :: nx, ny, nz 
    integer             :: mx, my, mz   
    integer             :: px, py, pz  
    integer             :: ex, ey, ez  

    real(double)        :: grid_spacing 
    real(double)        :: x, y, z, r
    real(double)        :: int, n, rest
    real(double)        :: xyz_delta(3)

    integer             :: count1, nsf1
    integer             :: ierr, stat
    integer             :: i, j, k, kk, ngrid   
    integer             :: ijk_delta(3), ijk(3), kji(3), njk(3), nji(3)

    integer             :: gto_lmn, index, l1
    real(double)        :: gto_val, gto_sph, gto_rad, gto_cart
    real(double)        :: ak, dk, Nk, factor

    real(double),  parameter :: a1g_norm  = sqrt(one/(four*pi))
    real(double),  parameter :: t1u_norm  = sqrt(three/(four*pi))

    
    ngrid        = 2*extent+1
    grid_spacing = r_int/real(extent,double)                   
    phi_on_grid  = zero

    mx = -extent ; px = extent
    my = mx      ; py = px 
    mz = my      ; pz = py

    call start_timer(tmr_std_exx_evalgto)
    !
    factor = 1.0d0
    !l1 = nx + ny + nz
    !if (l1 == 0) factor = a1g_norm 
    !if (l1 == 1) factor = t1u_norm    
    !
    grid_x_loop: do nx = mx, px
       x =  xyz(1) + real(nx,double)*grid_spacing + rst(1)
       
       grid_y_loop: do ny = my, py
          y = xyz(2) + real(ny,double)*grid_spacing + rst(2)
          
          grid_z_loop: do nz = mz, pz
             z = xyz(3) + real(nz,double)*grid_spacing + rst(3)
             
             r = sqrt(x*x+y*y+z*z)                                
             
             gto_cart = (x**real(ex))*(y**real(ey))*(z**real(ez))
             
             gto_rad = exp(-ak*r**2)
                
             gto_val = gto_rad * gto_cart                   
             
             phi_on_grid(&
                  nx+extent+1,  &
                  ny+extent+1,  &
                  nz+extent+1) = gto_val * factor
             
          end do grid_z_loop
       end do grid_y_loop
    end do grid_x_loop

    call stop_timer(tmr_std_exx_evalgto,.true.)

    return 
  end subroutine exx_gto_on_grid_prim
  !
  
end module exx_evalgto
