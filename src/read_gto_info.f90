module read_gto_info

  use datatypes
  use numbers,       ONLY: pi,two, one, zero
  use global_module, ONLY: io_lun
  use timer_stdclocks_module, &
       ONLY: start_timer,stop_timer,tmr_std_initialisation

  ! Modifications in:
  ! 1. created gto_format.f90 (cf. pao_format.f90)
  ! 2. pseudo_tm_info.f90/ allocate derived type gto(n_species)
  ! 3. initial_read.f90/ allocate gto_file(n_species)
  ! 4. ionic_data.f90/ call read_gto(...)

  implicit none

  character(len=1), dimension(0:5), parameter :: l_name_lowercase = (/'s','p','d','f','g','h'/)
  character(len=1), dimension(0:5), parameter :: l_name_uppercase = (/'S','P','D','F','G','H'/)
  integer, parameter :: max_angmom = 5

contains

  subroutine read_gto_new(inode,ionode,n_species)

    use numbers       , ONLY: pi, four, three, one
    use gto_format_new, ONLY: gto
    use species_module, ONLY: gto_file, species_label
    use GenComms,       ONLY: cq_abort, gcopy, my_barrier
    use exx_erigto,     ONLY: norm_gto, overlap_gto

    implicit none

    ! Arguments
    integer, intent(in)  :: inode, ionode, n_species

    ! Local variables 
    integer :: ios, match, tmp_int, gto_unit, nsp
    character(len=256) :: match_atom, tmp
    character(len=8)   :: charx, chary, charz

    integer :: i, j, k, l, lambda
    integer :: x, y, z, index, count
    integer :: l1, m1, n1, l2, m2, n2
    real(double) :: x1, y1, z1    
    real(double) :: x2, y2, z2    
    real(double) :: test, norm, S_12, alpha1, alpha2

    integer, dimension(0:max_angmom) :: l_nzeta 

    character(len=2) :: tmp_nl
    character(len=1) :: tmp_l, tmp_n   
    integer          :: tmp_nG, tmp_kind, tmp_n_int
    integer          :: l_value, l_max, acz1, p
    real(double)     :: tmp_occ
    real(double), dimension(0:1), parameter :: l_cart_norm = (/ sqrt(one/(four*pi)), sqrt(three/(four*pi))/)
    
    !
    ! gto_file(nsp) allocated in initial_read.module
    ! gto(nsp)      allocated in pseudo_tm_info
    !

    gto_unit = 1000

    call start_timer(tmr_std_initialisation)

    species_loop: do nsp = 1, n_species

       l_nzeta = 0
       !
       if(inode == ionode) then   
          !
          ! setup gto filename, eg for H atom we get H.gto
          gto_file(nsp)   = trim(species_label(nsp))//trim('.gto')
          gto(nsp)%label  = trim(species_label(nsp))
          !
          open(unit=gto_unit,file=gto_file(nsp),status="old",action="read", &
               iostat=ios,position="rewind")  
          if(ios /= 0) call cq_abort('read_gto: failed to open input file')
          !
          match = 0
          gto(nsp)%nsf_tot = 0
          !
          ! search for atom name
          do while (match == 0)
             read(unit=gto_unit,fmt=*,iostat=ios) match_atom
             match_atom = trim(match_atom)
             if (ios < 0) call cq_abort('read_gto: end error GTO data file', gto_unit)        
             if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)
             if (match_atom == gto(nsp)%label) then
                match = 1
             end if
          end do
          !
          ! get total number of zeta
          read(gto_unit,fmt=*,iostat=ios) gto(nsp)%n_zeta_tot
          if (ios < 0) call cq_abort('read_gto: end error GTO data file',gto_unit)        
          if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)
          !
          write(*,*)
          write(*,fmt='("#### read GTO for ",A2,"/n_zeta_tot =",I3," ####")') gto(nsp)%label, gto(nsp)%n_zeta_tot
          !
       end if
       !
       if(inode == ionode) then   
          !
          ! search for greatest angular momentum among zetas
          do i = 1, gto(nsp)%n_zeta_tot
             !
             read(unit=gto_unit,fmt=*,iostat=ios) tmp_nl, tmp_nG, tmp_occ, tmp_kind
             if (ios < 0) call cq_abort('read_gto: end error GTO data file', gto_unit)        
             if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)
             !
             ! convert angular momentum name to integer
             tmp_n = tmp_nl(1:1)
             tmp_l = tmp_nl(2:2)
             !
             call angmom_convert(tmp_l,l_value)
             print*, 'inode', inode, 'nsp', nsp, 'n_zeta', i,'tmp_nl', tmp_nl, 'l_value', l_value
             !
             l_max = 0
             ! get max angular momentum value
             if ( l_value > l_max ) then
                l_max = l_value
             end if
             !
             ! count and store number of zeta for each l value
             do l = 0, max_angmom
                if ( l_value == l ) then
                   l_nzeta(l) = l_nzeta(l) + 1
                end if
             end do
             !
             ! pass through GTO's parameters
             do j = 1, tmp_nG
                read(unit=gto_unit,fmt=*,iostat=ios)                
                if (ios < 0) call cq_abort('read_gto: end error GTO data file', gto_unit)        
                if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)
                !
             end do
             !
          end do
          !
          gto(nsp)%greatest_angmom = l_max
          !if (inode == ionode) then
          !print*,'inode', inode, 'nsp', nsp, 'gto(nsp)%greatest_angmom', gto(nsp)%greatest_angmom
          !print*,'inode', inode, 'nsp', nsp, 'l_nzeta', l_nzeta
          !end if
          !
          if ( sum(l_nzeta) /= gto(nsp)%n_zeta_tot ) &
               call cq_abort('read_gto: problem with number of zeta', gto_unit)
          !
       end if
       ! Broadcast to all the proc.
       call gcopy( gto(nsp)%label, 2)
       call gcopy( gto(nsp)%greatest_angmom )
       call gcopy( gto(nsp)%n_zeta_tot      )
       !
       print*,'inode', inode, 'nsp', nsp, 'gto(nsp)%greatest_angmom', gto(nsp)%greatest_angmom
       print*,'inode', inode, 'nsp', nsp, 'gto(nsp)%label', gto(nsp)%label       
       !
       allocate( gto(nsp)%angmom( 0:gto(nsp)%greatest_angmom ) )
       !
       if(inode == ionode) then
          do i = 0, gto(nsp)%greatest_angmom
             gto(nsp)%angmom(i)%n_zeta_in_angmom = l_nzeta(i)   
             gto(nsp)%angmom(i)%l_value = i                             
             gto(nsp)%angmom(i)%l_name  = l_name_uppercase(i)
             !print*,'inode', inode, 'nsp', nsp, 'angmom', i
             !print*,'inode', inode, '   gto(nsp)%angmom(i)%n_zeta_in_angmom', gto(nsp)%angmom(i)%n_zeta_in_angmom
             !print*,'inode', inode, '   gto(nsp)%angmom(i)%l_value', gto(nsp)%angmom(i)%l_value
             !print*,'inode', inode, '   gto(nsp)%angmom(i)%l_name ', gto(nsp)%angmom(i)%l_name
             
          end do
       end if
       !
       do i = 0, gto(nsp)%greatest_angmom
          call gcopy( gto(nsp)%angmom(i)%n_zeta_in_angmom )
          call gcopy( gto(nsp)%angmom(i)%l_value  )
          call gcopy( gto(nsp)%angmom(i)%l_name, 1)
          !
          print*,'inode', inode, 'nsp', nsp, 'angmom', i
          print*,'inode', inode,'nsp', nsp, 'i', i,  gto(nsp)%angmom(i)%n_zeta_in_angmom
          print*,'inode', inode, '   gto(nsp)%angmom(i)%l_value', gto(nsp)%angmom(i)%l_value
          print*,'inode', inode, '   gto(nsp)%angmom(i)%l_name ', gto(nsp)%angmom(i)%l_name
          !
          allocate( gto(nsp)%angmom(i)%zeta( gto(nsp)%angmom(i)%n_zeta_in_angmom ) )
          !
       end do

       if ( inode == ionode ) then
          !
          ! rewind and perform reading to grad GTO parameters
          rewind(gto_unit,iostat=ios)
          if (ios > 0) call cq_abort('read_gto: rewind error',gto_unit)
          !
          ! skip first 2 lines
          read(gto_unit,fmt=*,iostat=ios)
          read(gto_unit,fmt=*,iostat=ios)
       end if
       !
       do l = 0, gto(nsp)%greatest_angmom
          !
          do i = 1, gto(nsp)%angmom(l)%n_zeta_in_angmom
             !
             if ( inode == ionode ) then
                !
                read(unit=gto_unit,fmt=*,iostat=ios) tmp_nl, tmp_nG, tmp_occ, tmp_kind
                if (ios < 0) call cq_abort('read_gto: end error GTO data file', gto_unit)
                if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)
                !
                gto(nsp)%angmom(l)%zeta(i)%ngto = tmp_nG
                !gto(nsp)%angmom(l)%occ  = tmp_occ
                !gto(nsp)%angmom(l)%kind = tmp_kind
                !
                !tmp_n = tmp_nl(1:1)
                !read(tmp,*) tmp_n_int
                !gto(nsp)%angmom(l)%n = tmp_n_int
                                
                !
             end if
             !
             call gcopy(gto(nsp)%angmom(l)%zeta(i)%ngto)
             !call gcopy(gto(nsp)%angmom(l)%occ )
             !call gcopy(gto(nsp)%angmom(l)%kind)
             !call gcopy(gto(nsp)%angmom(l)%n   )             
             !
             print*, 'inode', inode, 'tmp_nG', tmp_nG, gto(nsp)%angmom(l)%zeta(i)%ngto

             allocate( gto(nsp)%angmom(l)%zeta(i)%a ( gto(nsp)%angmom(l)%zeta(i)%ngto) )
             gto(nsp)%angmom(l)%zeta(i)%a = zero
             allocate( gto(nsp)%angmom(l)%zeta(i)%d ( gto(nsp)%angmom(l)%zeta(i)%ngto) )
             gto(nsp)%angmom(l)%zeta(i)%d = zero
             allocate( gto(nsp)%angmom(l)%zeta(i)%c ( gto(nsp)%angmom(l)%zeta(i)%ngto) )
             gto(nsp)%angmom(l)%zeta(i)%c = zero
             !
             !if ( inode == ionode ) then
             print*,'inode', inode, 'nsp', nsp, 'l', l, 'i', i, 'gto(nsp)%angmom(l)%zeta(i)%a', gto(nsp)%angmom(l)%zeta(i)%a
             !end if

             !
             if ( inode == ionode ) then
                !             
                do j = 1, gto(nsp)%angmom(l)%zeta(i)%ngto
                   !
                   read(unit=gto_unit,fmt=*,iostat=ios)  &
                        gto(nsp)%angmom(l)%zeta(i)%a(j), &
                        gto(nsp)%angmom(l)%zeta(i)%d(j), &
                        gto(nsp)%angmom(l)%zeta(i)%c(j)
                   if (ios < 0) call cq_abort('read_gto: end error GTO data file', gto_unit)        
                   if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)
                   print*,'inode', inode, 'G', j, 'a', gto(nsp)%angmom(l)%zeta(i)%a(j), 'd', gto(nsp)%angmom(l)%zeta(i)%d(j), &
                        'c', gto(nsp)%angmom(l)%zeta(i)%c(j)
                   !
                end do
                !
             end if
             !
             !if ( inode == ionode ) then
             !print*,'inode', inode, 'gto(nsp)%angmom(l)%zeta(i)%a', gto(nsp)%angmom(l)%zeta(i)%a, gto(nsp)%angmom(l)%zeta(i)%ngto
             !end if
             print*,'inode', inode, 'gto(nsp)%angmom(l)%zeta(i)%ngto', gto(nsp)%angmom(l)%zeta(i)%ngto
             !print*, 'gto(nsp)%angmom(l)%zeta(i)%a', gto(nsp)%angmom(l)%zeta(i)%a, gto(nsp)%angmom(l)%zeta(i)%ngto
             !print*, 'gto(nsp)%angmom(l)%zeta(i)%d', gto(nsp)%angmom(l)%zeta(i)%d, gto(nsp)%angmom(l)%zeta(i)%ngto
             !print*, 'gto(nsp)%angmom(l)%zeta(i)%C', gto(nsp)%angmom(l)%zeta(i)%c, gto(nsp)%angmom(l)%zeta(i)%ngto
             
             !do j = 1, gto(nsp)%angmom(l)%zeta(i)%ngto
             call gcopy( gto(nsp)%angmom(l)%zeta(i)%a, gto(nsp)%angmom(l)%zeta(i)%ngto )
             call gcopy( gto(nsp)%angmom(l)%zeta(i)%d, gto(nsp)%angmom(l)%zeta(i)%ngto )
             call gcopy( gto(nsp)%angmom(l)%zeta(i)%c, gto(nsp)%angmom(l)%zeta(i)%ngto )
             !end do
             !
             print*,'gcopy inode', inode, 'gto(nsp)%angmom(l)%zeta(i)%a', gto(nsp)%angmom(l)%zeta(i)%a,&
                  gto(nsp)%angmom(l)%zeta(i)%ngto
             !
             !print*,  gto(nsp)%angmom(l)%zeta(i)%a
             !print*,  gto(nsp)%angmom(l)%zeta(i)%d
             !print*,  gto(nsp)%angmom(l)%zeta(i)%c
                          
          end do
          !
          if ( inode == ionode ) then
             gto(nsp)%angmom(l)%nf = (l+1)*(l+2)/2
             !
             do i = 1, gto(nsp)%angmom(l)%n_zeta_in_angmom             
                gto(nsp)%nsf_tot = gto(nsp)%nsf_tot + gto(nsp)%angmom(l)%nf
             end do
             !
          end if
          !
          call gcopy( gto(nsp)%angmom(l)%nf )
          print*,'inode', inode, 'nsp', nsp, 'l', l, 'gto(nsp)%angmom(l)%nf', gto(nsp)%angmom(l)%nf
          !          !
          allocate( gto(nsp)%angmom(l)%nx( (l+1)*(l+2)/2 ) )
          allocate( gto(nsp)%angmom(l)%ny( (l+1)*(l+2)/2 ) )
          allocate( gto(nsp)%angmom(l)%nz( (l+1)*(l+2)/2 ) )
          allocate( gto(nsp)%angmom(l)%nt( (l+1)*(l+2)/2 ) )
          allocate( gto(nsp)%angmom(l)%norm( (l+1)*(l+2)/2 ) )
          !
          if ( inode == ionode ) then
             !
             call cart_gauss_function( inode, ionode, &
                  gto(nsp)%angmom(l)%l_value, &
                  gto(nsp)%angmom(l)%l_name,  &
                  0.0d0,&
                  gto(nsp)%angmom(l)%nx, &
                  gto(nsp)%angmom(l)%ny, &
                  gto(nsp)%angmom(l)%nz, &
                  gto(nsp)%angmom(l)%nt )
             !
          end if
          !
          call gcopy( gto(nsp)%angmom(l)%nx, gto(nsp)%angmom(l)%nf )
          call gcopy( gto(nsp)%angmom(l)%ny, gto(nsp)%angmom(l)%nf )
          call gcopy( gto(nsp)%angmom(l)%nz, gto(nsp)%angmom(l)%nf )
          print*,'inode', inode, 'nsp', nsp, 'l', l,'gto(nsp)%angmom(l)%nx',  gto(nsp)%angmom(l)%nx
          print*,'inode', inode, 'nsp', nsp, 'l', l,'gto(nsp)%angmom(l)%ny',  gto(nsp)%angmom(l)%ny
          print*,'inode', inode, 'nsp', nsp, 'l', l,'gto(nsp)%angmom(l)%nz',  gto(nsp)%angmom(l)%nz
          
          do j = 1,  gto(nsp)%angmom(l)%nf
             call gcopy( gto(nsp)%angmom(l)%nt(j), 8 )
          end do
          !
       end do
       !
       call gcopy( gto(nsp)%nsf_tot )
       !
       print*,'inode', inode, 'nsf_tot',  gto(nsp)%nsf_tot
       !
       allocate( gto(nsp)%sf( gto(nsp)%nsf_tot ) )
       allocate( gto(nsp)%sf( gto(nsp)%nsf_tot ) )
       allocate( gto(nsp)%sf( gto(nsp)%nsf_tot ) )
       !
       count = 1
       angu_loop: do l1 = 0, gto(nsp)%greatest_angmom
          
          zeta_loop: do acz1 = 1, gto(nsp)%angmom(l1)%n_zeta_in_angmom
             
             magn_loop: do m1 = 1, gto(nsp)%angmom(l1)%nf

                if ( inode == ionode ) then
                   gto(nsp)%sf( count )%nx   = gto(nsp)%angmom(l1)%nx( m1 )
                   gto(nsp)%sf( count )%ny   = gto(nsp)%angmom(l1)%ny( m1 )
                   gto(nsp)%sf( count )%nz   = gto(nsp)%angmom(l1)%nz( m1 )
                   gto(nsp)%sf( count )%nt   = gto(nsp)%angmom(l1)%nt( m1 )
                   gto(nsp)%sf( count )%ngto = gto(nsp)%angmom(l1)%zeta(acz1)%ngto
                   gto(nsp)%sf( count )%norm = l_cart_norm(l1)
                end if

                call gcopy( gto(nsp)%sf( count )%nx )
                call gcopy( gto(nsp)%sf( count )%ny )
                call gcopy( gto(nsp)%sf( count )%nz )
                call gcopy( gto(nsp)%sf( count )%nt, 8)
                call gcopy( gto(nsp)%sf( count )%ngto)
                call gcopy( gto(nsp)%sf( count )%norm)

                print*,'inode', inode, 'l1, acz1, m1', l1, acz1, m1, 'ex, ey, ez', &
                     gto(nsp)%sf( count )%nx, gto(nsp)%sf( count )%ny, gto(nsp)%sf( count )%nz, &
                     'nt', gto(nsp)%sf( count )%nt, 'norm', gto(nsp)%sf( count )%norm, 'count', count, &
                     'ngto', gto(nsp)%sf( count )%ngto
                                
                !if (.not. allocated( gto(nsp)%sf( count )%a ) ) &
                allocate( gto(nsp)%sf( count )%a( gto(nsp)%sf( count )%ngto ) ) 
                !if (.not. allocated( gto(nsp)%sf( count )%d ) ) &
                allocate( gto(nsp)%sf( count )%d( gto(nsp)%sf( count )%ngto ) )
                !if (.not. allocated( gto(nsp)%sf( count )%c ) ) &
                allocate( gto(nsp)%sf( count )%c( gto(nsp)%sf( count )%ngto ) )
                
                primitive: do p = 1, gto(nsp)%angmom(l1)%zeta(acz1)%ngto 

                   if ( inode == ionode ) then
                      
                      gto(nsp)%sf( count )%a( p ) = gto(nsp)%angmom(l1)%zeta(acz1)%a( p )
                      gto(nsp)%sf( count )%d( p ) = gto(nsp)%angmom(l1)%zeta(acz1)%d( p )
                      gto(nsp)%sf( count )%c( p ) = gto(nsp)%angmom(l1)%zeta(acz1)%c( p )
                                            
                   end if

                   !print*,'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%a( p )',  &
                   !     gto(nsp)%sf( count )%a( p )
                   !print*,'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%d( p )',  &
                   !     gto(nsp)%sf( count )%d( p )
                   !print*,'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%c( p )',  &
                   !     gto(nsp)%sf( count )%c( p )

                   
                end do primitive
                call gcopy(  gto(nsp)%sf( count )%a,  gto(nsp)%angmom(l1)%zeta(acz1)%ngto )
                call gcopy(  gto(nsp)%sf( count )%d,  gto(nsp)%angmom(l1)%zeta(acz1)%ngto )
                call gcopy(  gto(nsp)%sf( count )%c,  gto(nsp)%angmom(l1)%zeta(acz1)%ngto )

                print*,'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%a( p )',  &
                     gto(nsp)%sf( count )%a
                print*,'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%d( p )',  &
                     gto(nsp)%sf( count )%d
                print*,'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%c( p )',  &
                     gto(nsp)%sf( count )%c

                
                count = count + 1
             end do magn_loop
          end do zeta_loop
       end do angu_loop

       print*,'inode', inode, 'fin'
          
    end do species_loop

    
    call stop_timer(tmr_std_initialisation)

    return
  end subroutine read_gto_new
  !
  subroutine cart_gauss_function( inode, ionode, lambda, name, a, nx, ny, nz, nt )
    
    use exx_erigto,     ONLY: norm_gto
    
    implicit none

    integer,          intent(in) :: lambda
    integer,          intent(in) :: inode, ionode
    character(len=1), intent(in) :: name
    real(double),     intent(in) :: a
    
    integer, intent(out), dimension( (lambda+1)*(lambda+2)/2 ) :: nx, ny, nz
    !integer, intent(out), dimension( (lambda+1)*(lambda+2)/2 ) :: norm
    character(len=8), intent(out), dimension( (lambda+1)*(lambda+2)/2 ) :: nt

    
    character(len=8) :: charx, chary, charz
    integer :: x, y, z, k, index
    
    if(inode == ionode) then
       
       index = 1
       do x = 0, lambda
          do z = 0, lambda - x
             y = lambda - (x + z)
             !
             if(inode == ionode) then
                ! Exponents for Cartesian Gaussian functions
                nx(index) = x !x
                ny(index) = y !y
                nz(index) = z !z
                !
             end if
             !
             charx = '' ; chary = '' ; charz = '' 
             if (x > 0) then
                do  k = 1, x
                   charx = trim(charx)//'x' 
                end do
             end if
             if (y > 0) then
                do  k = 1, y
                   chary = trim(chary)//'y' 
                end do
             end if
             if (z > 0) then
                do  k = 1, z
                   charz = trim(charz)//'z' 
                end do
             end if
             !
             nt(index) = trim(name)//trim(charx)//trim(chary)//trim(charz)
             !
             !             
             !norm(index)= norm_gto( &
             !     a,     & ! alpha
             !     nx(index), & ! l
             !     ny(index), & ! m
             !     nz(index))   ! n
                !
             write(*,'(6x,"(nx, ny, nz) = (",I2,x,I2,x,I2,") =",x,A6,"/","norm = ",F14.8,I8)') &
                  nx(index), &
                  ny(index), &
                  nz(index), &
                  nt(index)
             !     norm(index), index
             
             !                      
             index = index + 1
             !
          end do
       end do
    end if
    !
  end subroutine cart_gauss_function

  
  subroutine angmom_convert(l_name,l_value)

    use GenComms, ONLY: cq_abort

    implicit none

    character(len=1), intent(inout) :: l_name    
    integer, intent(out) :: l_value
    !
    integer :: i, l, lambda
    !
    ! if necessary convert lowercase to uppercase
    do l = 0, size(l_name_lowercase)-1
       if ( l_name == l_name_lowercase(l)) then
          l_name  = l_name_uppercase(l)
       end if
    end do
    !
    l_value = -1
    do l = 0, size(l_name_uppercase)-1
       if ( l_name == l_name_uppercase(l)) then
          l_value = l
       end if
    end do
    !
    if ( l_value == -1 ) then
       call cq_abort('read_gto: error angular momentum convert')
    end if
    !
    !select case(l_name)             
    !case('S')
    !   lambda = 0
    !case('P')
    !   lambda = 1
    !case('D')
    !   lambda = 2
    !case('F')
    !   lambda = 3
    !case('G')
    !   lambda = 4
    !case('H')
    !   lambda = 5
    !case default
    !   call cq_abort('read_gto: error angular momentum') 
    !end select
    !l_value = lambda

    return
  end subroutine angmom_convert


  subroutine read_gto(inode,ionode,n_species)

    use gto_format,     ONLY: gto
    use species_module, ONLY: gto_file, species_label
    use GenComms,       ONLY: cq_abort, gcopy, my_barrier
    use exx_erigto,     ONLY: norm_gto, overlap_gto

    implicit none

    ! Arguments
    integer, intent(in)  :: inode, ionode, n_species

    ! Local variables 
    integer :: ios, match, tmp_int, gto_unit, nsp
    character(len=256) :: match_atom
    character(len=8)   :: charx, chary, charz

    integer :: i, j, k, l, lambda
    integer :: x, y, z, index, count
    integer :: l1, m1, n1, l2, m2, n2
    real(double) :: x1, y1, z1    
    real(double) :: x2, y2, z2    
    real(double) :: test, norm, S_12, alpha1, alpha2

    !
    ! gto_file(nsp) allocated in initial_read.module
    ! gto(nsp)      allocated in pseudo_tm_info
    !

    gto_unit = 1000

    call start_timer(tmr_std_initialisation)

    species_loop: do nsp = 1, n_species

       if(inode == ionode) then   

          gto_file(nsp)   = trim(species_label(nsp))//trim('.gto')
          gto(nsp)%label  = trim(species_label(nsp))
          !
          open(unit=gto_unit,file=gto_file(nsp),status="old",action="read", &
               iostat=ios,position="rewind")  
          if(ios /= 0) call cq_abort('read_gto: failed to open input file')
          !
          !call gcopy(gto(nsp)%label,2)     
          !
          match = 0
          count = 1
          !
          do while (match == 0)
             read(unit=gto_unit,fmt=*,iostat=ios) match_atom
             match_atom = trim(match_atom)
             if (ios < 0) call cq_abort('read_gto: end error GTO data file', gto_unit)        
             if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)
             if (match_atom == gto(nsp)%label) then
                match = 1
             end if
          end do
          !
          read(gto_unit,fmt=*,iostat=ios) gto(nsp)%nshell
          if (ios < 0) call cq_abort('read_gto: end error GTO data file',gto_unit)        
          if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)
          !
          write(*,*)
          write(*,fmt='("#### read GTO for ",A2,"/nshell =",I3," ####")') gto(nsp)%label, gto(nsp)%nshell
          !
          !write(*,fmt='(2x,"nshell = ",I3)') gto(nsp)%nshell
          !
       end if
       !
       call gcopy(gto(nsp)%label,2)
       call gcopy(gto(nsp)%nshell)
       !
       allocate(gto(nsp)%shell(gto(nsp)%nshell))
       !print*, 'node', inode, allocated(gto(nsp)%shell), size(gto(nsp)%shell)
       !  
       !
       gto(nsp)%nprim = 0
       gto(nsp)%nang  = 0
       !
       shell_loop: do i = 1, gto(nsp)%nshell
          !
          if(inode == ionode) then
             !
             read(gto_unit,fmt=*,iostat=ios) &
                  gto(nsp)%shell(i)%t, & ! type s,p,sp,d,f...
                  gto(nsp)%shell(i)%n, & ! number of primitive
                  gto(nsp)%shell(i)%f    ! scaling factor
             if (ios < 0) call cq_abort('read_gto: end error GTO data file', gto_unit)        
             if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)             
             !
             write(*,'(">>> shell ",I3)') i
             write(*,'(2x,"type           =",3x,A2)'   ) gto(nsp)%shell(i)%t
             write(*,'(2x,"scaling factor =",3x,F8.6)' ) gto(nsp)%shell(i)%f

             select case(gto(nsp)%shell(i)%t)
             case('S')
                lambda = 0
                !gto(nsp)%shell(i)%nn = 0
                !gto(nsp)%shell(i)%nf = 1
             case('P')
                lambda = 1
                !gto(nsp)%shell(i)%nn = 1
                !gto(nsp)%shell(i)%nf = 3
             case('D')
                lambda = 2
                !gto(nsp)%shell(i)%nn = 2
                !gto(nsp)%shell(i)%nf = 6
             case('F')
                lambda = 3
                !gto(nsp)%shell(i)%nn = 3
                !gto(nsp)%shell(i)%nf = 10
             case('G')
                lambda = 4
                !gto(nsp)%shell(i)%nn = 4
                !gto(nsp)%shell(i)%nf = 15
             case('H')
                lambda = 5
                !gto(nsp)%shell(i)%nn = 5
                !gto(nsp)%shell(i)%nf = 24
             case default
                call cq_abort('read_gto: error angular momentum',gto(nsp)%shell(i)%nn) 
             end select
             !
             gto(nsp)%shell(i)%nn     = lambda
             gto(nsp)%shell(i)%nf     = (lambda+1)*(lambda+2)/2
             gto(nsp)%nang            = gto(nsp)%nang  + gto(nsp)%shell(i)%nf
             gto(nsp)%nprim           = gto(nsp)%nprim + gto(nsp)%shell(i)%n*gto(nsp)%shell(i)%nf
             !
             write(*,'(2x,"angular momentum     l =",I3)') gto(nsp)%shell(i)%nn  
             write(*,'(2x,"number of values for m =",I3)') gto(nsp)%shell(i)%nf              
             write(*,'(2x,"number of primitives   =",I3)') gto(nsp)%shell(i)%n              
             !
             !
          end if
          !
          call gcopy(gto(nsp)%shell(i)%t,2)     
          call gcopy(gto(nsp)%shell(i)%n  )     
          call gcopy(gto(nsp)%shell(i)%f  )     
          call gcopy(gto(nsp)%shell(i)%nn )     
          call gcopy(gto(nsp)%shell(i)%nf )
          !
          allocate(gto(nsp)%shell(i)%ak( gto(nsp)%shell(i)%n ))
          allocate(gto(nsp)%shell(i)%dk( gto(nsp)%shell(i)%n ))                 
          allocate(gto(nsp)%shell(i)%Nk( gto(nsp)%shell(i)%n,gto(nsp)%shell(i)%nf ))
          !
          allocate(gto(nsp)%shell(i)%nx( gto(nsp)%shell(i)%nf ))
          allocate(gto(nsp)%shell(i)%ny( gto(nsp)%shell(i)%nf ))                 
          allocate(gto(nsp)%shell(i)%nz( gto(nsp)%shell(i)%nf ))                 
          allocate(gto(nsp)%shell(i)%nt( gto(nsp)%shell(i)%nf ))
          allocate(gto(nsp)%shell(i)%norm(  gto(nsp)%shell(i)%nf ))  
          !
          prim_loop: do j = 1, gto(nsp)%shell(i)%n 
             !
             if(inode == ionode) then

                read(gto_unit,fmt=*,iostat=ios) &
                     gto(nsp)%shell(i)%ak(j),   &
                     gto(nsp)%shell(i)%dk(j)
                if (ios < 0) call cq_abort('read_gto: end error GTO data file', gto_unit)        
                if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)

                write(*,'(2x,"primitive",x,"(a",I3,", d",I3,") = (",F14.8,x,F14.8,")")') j, j,  &
                     gto(nsp)%shell(i)%ak(j), &
                     gto(nsp)%shell(i)%dk(j)                

                !write(*,'(6x,"(a",I3,", d",I3,") = (",F14.8,x,F14.8,")")') j, j, &
                !gto(nsp)%shell(i)%ak(j), &
                !gto(nsp)%shell(i)%dk(j)
             end if
             !
             call gcopy(gto(nsp)%shell(i)%ak(j))
             call gcopy(gto(nsp)%shell(i)%dk(j))
             !
             index  = 1
             lambda = gto(nsp)%shell(i)%nn 
             do x = 0, lambda
                do y = 0, lambda - x
                   z = lambda - (x + y)
                   !
                   if(inode == ionode) then
                      ! Exponents for Cartesian Gaussian functions
                      gto(nsp)%shell(i)%nx(index) = x!x
                      gto(nsp)%shell(i)%ny(index) = y!y
                      gto(nsp)%shell(i)%nz(index) = z!z
                      !
                   end if
                   !
                   charx = '' ; chary = '' ; charz = '' 
                   if (x > 0) then
                      do  k = 1, x
                         charx = trim(charx)//'x' 
                      end do
                   end if
                   if (y > 0) then
                      do  k = 1, y
                         chary = trim(chary)//'y' 
                      end do
                   end if
                   if (z > 0) then
                      do  k = 1, z
                         charz = trim(charz)//'z' 
                      end do
                   end if
                   !
                   if(inode == ionode) then
                      !
                      gto(nsp)%shell(i)%nt(index) = &
                           trim(gto(nsp)%shell(i)%t)//trim(charx)//trim(chary)//trim(charz)
                      !
                   end if
                   !
                   call gcopy(gto(nsp)%shell(i)%nt(index),16)
                   !
                   if(inode == ionode) then                   
                      ! Compute the norm of each Cartesian components

                      gto(nsp)%shell(i)%Nk(j,index) = norm_gto( &
                           gto(nsp)%shell(i)%ak(j),     & ! alpha
                           gto(nsp)%shell(i)%nx(index), & ! l
                           gto(nsp)%shell(i)%ny(index), & ! m
                           gto(nsp)%shell(i)%nz(index))   ! n
                      !
                      write(*,'(6x,"(nx, ny, nz) = (",I2,x,I2,x,I2,") =",x,A6,"/","norm = ",F14.8,I8)') &
                           gto(nsp)%shell(i)%nx(index), &
                           gto(nsp)%shell(i)%ny(index), &
                           gto(nsp)%shell(i)%nz(index), &
                           gto(nsp)%shell(i)%nt(index), &
                           gto(nsp)%shell(i)%Nk(j,index), index

                   end if
                   !
                   call gcopy(gto(nsp)%shell(i)%Nk(j,index))
                   !                      
                   index = index + 1
                   count = count + 1
                   !
                end do
             end do
             !
          end do prim_loop
          !
          if(inode == ionode) then
             !
             gto(nsp)%shell(i)%norm(:) = zero
             !
             do j = 1, gto(nsp)%shell(i)%n
                !
                do k = 1, gto(nsp)%shell(i)%n

                   do l = 1, gto(nsp)%shell(i)%nf

                      alpha1 = gto(nsp)%shell(i)%ak(j)
                      l1     = gto(nsp)%shell(i)%nx(l)
                      m1     = gto(nsp)%shell(i)%ny(l)
                      n1     = gto(nsp)%shell(i)%nz(l)
                      x1 = zero ; y1 = zero ; z1 = zero

                      alpha2 = gto(nsp)%shell(i)%ak(k)
                      l2     = gto(nsp)%shell(i)%nx(l)
                      m2     = gto(nsp)%shell(i)%ny(l)
                      n2     = gto(nsp)%shell(i)%nz(l)
                      x2 = zero ; y2 = zero ; z2 = zero

                      S_12   = overlap_gto(alpha1,l1,m1,n1,x1,y1,z1, &
                           alpha2,l2,m2,n2,x2,y2,z2)

                      !print*, j, k, l1, m1, n1, l2, m2, n2, alpha1, alpha2, l, S_12

                      gto(nsp)%shell(i)%norm(l) = gto(nsp)%shell(i)%norm(l) + &
                           gto(nsp)%shell(i)%ak(j) * gto(nsp)%shell(i)%ak(k) * S_12 
                      !
                   end do
                   !
                end do
                !
             end do
             !
          end if
          !
          if(inode == ionode) then
             !
             write(*,'(2x,"number of angular funtions = ",I3)') gto(nsp)%shell(i)%nf
             !
             do l = 1, gto(nsp)%shell(i)%nf
                write(*,'(2x,"norm =",F12.6)') gto(nsp)%shell(i)%norm(l)
             end do
             !
          end if
          !
          call gcopy( gto(nsp)%shell(i)%norm, gto(nsp)%shell(i)%nf)
          !
          !
       end do shell_loop

       if(inode == ionode) then
          !
          write(*,'(">>> total number of angular functions =",I10)') gto(nsp)%nang
          write(*,'(">>> total number of primitives        =",I10)') gto(nsp)%nprim
          !
          close(unit=gto_unit,iostat=ios)  
          if(ios /= 0) call cq_abort('read_gto: failed to close input file',gto_unit)
       end if

       call gcopy(gto(nsp)%nang )
       call gcopy(gto(nsp)%nprim)

    end do species_loop

    call stop_timer(tmr_std_initialisation)

    return
  end subroutine read_gto
  !  
end module read_gto_info
