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

    use angular_coeff_routines, ONLY: evaluate_pao, re_cart_norm
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
    
    do i = 1, gto( spec )%nsf_gto_tot

       ex = gto( spec )%sf( i )%nx
       ey = gto( spec )%sf( i )%ny
       ez = gto( spec )%sf( i )%nz

       factor = re_cart_norm( l1,m1-l1-1 ) !call re_cart_norm(l1,m1-l1-1,factor)
       
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

                   !gto_rad = gto_rad + exp(-a*r**2)
                   gto_rad = gto_rad + d*exp(-a*r**2)
                   
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

    use angular_coeff_routines, ONLY: evaluate_pao, re_cart_norm
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
    integer             :: i, j, p, l1, m1, acz1, n1, ngrid
    integer             :: ijk_delta(3), ijk(3), kji(3), njk(3), nji(3)

    real(double)        :: gto_rad, gto_cart, gto_val
    real(double)        :: a, d, factor, c

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

       !if (l1 == 0) factor = a1g_norm
       !if (l1 == 1) factor = t1u_norm 

       zeta_loop: do acz1 = 1, gto(spec)%angmom(l1)%n_zeta_in_angmom

          !magn_loop: do m1 = 1, gto(spec)%angmom(l1)%nf_gto

          magn_loop: do m1 = -l1, +l1

             factor = re_cart_norm(l1,m1)! call re_cart_norm(l1,m1-l1-1,factor)

             !sph_loop: do n1 = 1, gto(nsp)%angmom(l1)%transform_sph(m1)%size

             !ex = gto(spec)%angmom(l1)%nx(m1)
             !ey = gto(spec)%angmom(l1)%ny(m1)
             !ez = gto(spec)%angmom(l1)%nz(m1)

             !factor = re_cart_norm(l1,m1)! call re_cart_norm(l1,m1-l1-1,factor)

             grid_x_loop: do nx = mx, px
                x = xyz(1) + real(nx,double)*grid_spacing + rst(1)

                grid_y_loop: do ny = my, py
                   y = xyz(2) + real(ny,double)*grid_spacing + rst(2)

                   grid_z_loop: do nz = mz, pz
                      z = xyz(3) + real(nz,double)*grid_spacing + rst(3)
                      
                      r = sqrt(x*x+y*y+z*z)  

                      gto_cart = zero
                      sph_loop: do n1 = 1, gto(spec)%angmom(l1)%transform_sph(m1)%size

                         ex = gto(spec)%angmom(l1)%transform_sph(m1)%nx(n1)
                         ey = gto(spec)%angmom(l1)%transform_sph(m1)%ny(n1)
                         ez = gto(spec)%angmom(l1)%transform_sph(m1)%nz(n1)
                         c  = gto(spec)%angmom(l1)%transform_sph(m1)%c( n1)

                         gto_cart = gto_cart + c * x**real(ex) * y**real(ey) * z**real(ez)

                      end do sph_loop
                      !
                      gto_rad = zero 
                      primitive: do p = 1, gto(spec)%angmom(l1)%zeta(acz1)%ngto 

                         a = gto(spec)%angmom(l1)%zeta(acz1)%a(p)
                         d = gto(spec)%angmom(l1)%zeta(acz1)%d(p)

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
!!$  subroutine re_cart_norm(l,m,norm)
!!$    
!!$    use datatypes
!!$    use numbers
!!$    use GenComms,ONLY: cq_abort
!!$
!!$    implicit none
!!$    !
!!$    !
!!$    integer,      intent(in) :: l, m
!!$    real(double), intent(out):: norm        
!!$    !
!!$    real(double),  parameter :: a1g_norm  = sqrt(one/(four*pi))
!!$    real(double),  parameter :: t1u_norm  = sqrt(three/(four*pi))
!!$    real(double),  parameter :: t2g_norm  = sqrt(fifteen/(four*pi))
!!$    real(double),  parameter :: eg_a_norm = sqrt(fifteen/(sixteen*pi)) 
!!$    real(double),  parameter :: eg_b_norm = sqrt(five/(sixteen*pi))     
!!$    !
!!$    real(double),  parameter :: f0_norm = sqrt(seven/(sixteen*pi))
!!$    real(double),  parameter :: f1_norm = sqrt( 21.0_double/(sixteen*two*pi))
!!$    real(double),  parameter :: f2_norm = sqrt(105.0_double/(sixteen*pi)) 
!!$    real(double),  parameter :: f3_norm = sqrt( 35.0_double/(sixteen*two*pi)) 
!!$    !
!!$    integer :: i
!!$    !
!!$    if(l == 0) then !s-type function
!!$       !
!!$       norm = a1g_norm
!!$       !
!!$    else if(l == 1) then !p-type function
!!$       !
!!$       norm = t1u_norm
!!$       !
!!$    else if(l == 2) then !d-type function
!!$       !
!!$       select case( m )
!!$
!!$       case( 3 ) ! d_{x2-y2}
!!$          norm = 1.0d0
!!$       case( 2 ) ! d_{x2-y2}
!!$          norm = eg_a_norm !(x*x-y*y) 
!!$       case( 1 ) ! d_{xz}
!!$          norm = t2g_norm !(x*z)         
!!$       case( 0 ) ! d_{z2}
!!$          norm = eg_b_norm!(3*z*z-r*r)
!!$       case(-1 ) ! d_{yz}
!!$          norm = t2g_norm !(y*z)
!!$       case(-2 ) ! d_{xy}             
!!$          norm =-t2g_norm !(x*y) ! take care phase factor
!!$       case default
!!$          call cq_abort('re_cart_norm/problem with (l,m) =',l,m)
!!$       end select
!!$       !
!!$    else if(l == 3) then !f-type function
!!$       !
!!$       select case( m )
!!$       case( 3 ) ! f_{x3-xy2}
!!$          norm = f3_norm !*(x*x - 3*y*y)*x
!!$       case( 2 ) ! f_{zx2-zy2}
!!$          norm = f2_norm !*(x*x - y*y)*z
!!$       case( 1 ) ! f_{xz2}
!!$          norm = f1_norm !*(5*z*z - r*r)*x
!!$       case( 0 ) ! f_{z3}
!!$          norm = f0_norm !*(5*z*z - 3*r*r)*z
!!$       case(-1 ) ! f_{yz2}             
!!$          norm = f1_norm !(5*z*z - r*r)*y
!!$       case(-2 ) ! f_{xyz}             
!!$          norm = f2_norm !x*y*z
!!$       case(-3 ) ! f_{3yx2-y3}             
!!$          norm = f3_norm !(3*x*x - y*y)*y
!!$       case default
!!$          call cq_abort('re_cart_norm/problem with (l,m) =',l,m)
!!$       end select
!!$       !
!!$    else if( l > 3) then
!!$       call cq_abort('re_cart_norm/not implemented for l > 3')
!!$    else                               
!!$       call cq_abort('re_cart_norm/problem with l =',l)
!!$    end if
!!$    
!!$  end subroutine re_cart_norm
  
end module exx_evalgto
