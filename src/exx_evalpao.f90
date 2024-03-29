module exx_evalpao

    use datatypes
    use GenComms,     only: cq_abort
    use timer_module, only: start_timer, stop_timer

    implicit none

contains

!!   2016/07/29 18:30 nakata
!!    Renamed supports_on_atom -> blips_on_atom
!!   2016/08/08 15:30 nakata
!!    Removed unused supportfns
!!   2016/12/29 21:30 nakata
!!    Used pao information instead of blips_on_atom, and 
!!    removed support_spec_format (blips_on_atom), which is no longer needed
  subroutine exx_phi_on_grid(inode,atom,spec,extent,xyz,nsuppfuncs,phi_on_grid,r_int,rst)

    use numbers,      only: zero, one, two, three, four, five, six, fifteen, sixteen
    use numbers,      only: half, quarter, very_small, pi, twopi, fourpi
    use exx_types,    only: unit_matrix_write, tmr_std_exx_evalpao
    use exx_types,    only: exx_overlap, exx_cartesian

    use angular_coeff_routines, only: evaluate_pao
    use dimens,                 only: r_h   
    use pao_format,             only: pao

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

    ! << Local variables >>
    real(double), dimension(nsuppfuncs) :: sfsum

    integer             :: nx, ny, nz 
    integer             :: mx, my, mz   
    integer             :: px, py, pz   
    integer             :: l1, m1, acz

    real(double)        :: grid_spacing 
    real(double)        :: x, y, z, r
    real(double)        :: int, n, rest
    real(double)        :: xyz_delta(3)

    integer             :: count1, nsf1
    integer             :: ierr, stat
    integer             :: i, j, ngrid   
    integer             :: ijk_delta(3), ijk(3), kji(3), njk(3), nji(3)
    integer             :: n1, n2, npts

    real(double)        :: xj1, xj2, a, b, c, d, del_r, alpha
    real(double)        :: a_table, b_table, c_table, d_table

    ! << Called subroutine variables >>
    real(double)        :: pao_val, y_val    

    real(double),  parameter :: a1g_norm  = sqrt(one/(four*pi))
    real(double),  parameter :: t1u_norm  = sqrt(three/(four*pi))
    real(double),  parameter :: t2g_norm  = sqrt(fifteen/(four*pi))
    real(double),  parameter :: eg_a_norm = sqrt(fifteen/(sixteen*pi)) 
    real(double),  parameter :: eg_b_norm = sqrt(five/(sixteen*pi))     

    ngrid        = 2*extent+1
    grid_spacing = r_int/real(extent,double)                   
    phi_on_grid  = zero
    sfsum        = zero

    mx = -extent ; px = extent
    my = mx      ; py = px 
    mz = my      ; pz = py

    call start_timer(tmr_std_exx_evalpao)    

    !Define the overlap box
    overlap_box: if (exx_overlap) then
       do i = 1, 3
          n    = xyz(i)/ grid_spacing
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
          !write(*,'(5F12.6,4I6)') xyz(i), n, int, rest, &
          !     xyz_delta(i), ijk(i), njk(i), kji(i), nji(i)
          !write(*,'(4I6)') ijk(i), njk(i), kji(i), nji(i)

       end do
       mx = mx +kji(1)-1
       my = my +kji(2)-1
       mz = mz +kji(3)-1
       px = px -ijk(1)+1
       py = py -ijk(2)+1
       pz = pz -ijk(3)+1
    end if overlap_box

    grid_x_loop: do nx = mx, px
       x = xyz(1) + real(nx,double)*grid_spacing + rst(1)

       grid_y_loop: do ny = my, py
          y = xyz(2) + real(ny,double)*grid_spacing + rst(2)

          grid_z_loop: do nz = mz, pz
             z = xyz(3) + real(nz,double)*grid_spacing + rst(3)

             !norm = sqrt((x-xyz(1))**2+(y-xyz(2))**2+(z-xyz(3))**2)
             !if (norm <= r_h) then

             r = sqrt(x*x+y*y+z*z)             
             !if(r < very_small) then
             !   r = zero
             !end if

             !print*, '1 cycle start'

             count1  = 1
             sfsum   = zero
             angu_loop: do l1 = 0, pao(spec)%greatest_angmom

                zeta_loop: do acz = 1, pao(spec)%angmom(l1)%n_zeta_in_angmom

                   magn_loop: do m1 = -l1, l1                      

                      pao_val = zero
                      y_val   = zero             
                      rec_ylm: if (.not.exx_cartesian) then                      
                         call evaluate_pao(0,spec,l1,acz,m1,x,y,z,pao_val)
                      else
                         npts  = pao(spec)%angmom(l1)%zeta(acz)%length
                         del_r = (pao(spec)%angmom(l1)%zeta(acz)%cutoff/&
                              (pao(spec)%angmom(l1)%zeta(acz)%length-1))

                         n1 = floor(r/del_r)
                         n2 = n1 + 1

                         if(n2 > npts-1) then
                            pao_val = 0.0_double

                         else
                            xj1 = n1  * del_r
                            xj2 = xj1 + del_r

                            a = (xj2 - r)/del_r
                            b = 1.0_double - a
                            c = a*(a*a - 1.0_double)*del_r*del_r/6.0_double
                            d = b*(b*b - 1.0_double)*del_r*del_r/6.0_double

                            a_table = pao(spec)%angmom(l1)%zeta(acz)%table(n1+1)
                            b_table = pao(spec)%angmom(l1)%zeta(acz)%table(n2+1)
                            c_table = pao(spec)%angmom(l1)%zeta(acz)%table2(n1+1)
                            d_table = pao(spec)%angmom(l1)%zeta(acz)%table2(n2+1)

                            pao_val = a*a_table + b*b_table + c*c_table + d*d_table

                            if (l1 == 0) then !s-type function
                               y_val = a1g_norm

                            else if(l1 == 1) then !p-type function
                               r_p: if (r > zero) then                               
                                  select case(m1)
                                  case( 1) !p_x
                                     y_val = t1u_norm*x
                                  case( 0) !p_z
                                     y_val = t1u_norm*z 
                                  case(-1) !p_y
                                     y_val = t1u_norm*y 
                                  case default
                                     call cq_abort('exx_phi_on_grid/Unrecognized magn. m ',m1)
                                  end select
                               else
                                  y_val = zero
                               end if r_p
                            else if(l1 == 2) then !d-type function                                 
                               r_d: if (r > zero) then
                                  select case(m1)
                                  case( 2) ! d_{x2-y2}
                                     y_val =  eg_a_norm*(x*x-y*y) 
                                  case( 1) ! d_{xz}
                                     y_val =  t2g_norm*(x*z)         
                                  case( 0) ! d_{z2}
                                     y_val =  eg_b_norm*(3.0d0*z*z-r*r)
                                  case(-1) ! d_{yz}
                                     y_val =  t2g_norm*(y*z)
                                  case(-2) ! d_{yz}             
                                     y_val = -t2g_norm*(x*y) ! Take care phase factor
                                  case default
                                     call cq_abort('exx_phi_on_grid/Unrecognized mag. m ',m1)
                                  end select
                               else
                                  y_val = zero
                               end if r_d
                            else                               
                               call cq_abort('exx_phi_on_grid/Unrecognized angul. l ',l1)
                            end if

                            pao_val = pao_val*y_val
                         end if
                      end if rec_ylm

                      ! Put pao_val directly into phi_on_grid
                      ! (only for primitive PAOs and not for blips)
                      phi_on_grid(nx+extent+1,ny+extent+1,nz+extent+1,count1) = pao_val

                      count1 = count1 + 1
                   end do magn_loop
                end do zeta_loop
             end do angu_loop

          end do grid_z_loop
       end do grid_y_loop
    end do grid_x_loop

    call stop_timer(tmr_std_exx_evalpao,.true.)

    return 
  end subroutine exx_phi_on_grid

end module exx_evalpao
