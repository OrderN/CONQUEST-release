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

    use angular_coeff_routines, ONLY: re_cart_norm
    use numbers,        ONLY: pi, four, three, one
    use gto_format_new, ONLY: gto
    !use exx_evalgto,    ONLY: re_cart_norm
    use species_module, ONLY: gto_file, species_label
    use GenComms,       ONLY: cq_abort, gcopy, my_barrier
    use exx_erigto,     ONLY: norm_gto, overlap_gto

    implicit none

    ! Arguments
    integer, intent(in)  :: inode, ionode, n_species

    ! Local variables 
    integer :: ios, match, tmp_int, gto_unit, nsp
    character(len=256) :: match_atom, tmp
    !character(len=8)   :: charx, chary, charz

    integer :: i, j, k, l, lambda
    integer :: x, y, z, index, count, count2
    integer :: l1, m1, n1, l2, m2, n2
    
    real(double) :: x1, y1, z1    
    real(double) :: x2, y2, z2    
    real(double) :: test, norm, S_12, alpha1, alpha2

    integer, dimension(0:max_angmom) :: l_nzeta 

    character(len=2) :: tmp_nl
    character(len=1) :: tmp_l, tmp_n   
    integer          :: tmp_nG, tmp_kind, tmp_n_int
    integer          :: l_value, l_max, acz1, p, sph_size
    integer          :: sph_nx(5), sph_ny(5), sph_nz(5)
    real(double)     :: tmp_occ, sph_coef(5)
    character(len=12), dimension(0:2) :: zeta_kind = (/'regular     ','polarisation','other       '/)
    character(len=12):: sph_name

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
          write(io_lun,fmt='(10x,"######## read GTO for ",A2,"/n_zeta_tot =",I3," ########")') &
               gto(nsp)%label, gto(nsp)%n_zeta_tot
          !
       end if
       !
       if(inode == ionode) then   
          !
          gto(nsp)%nsf_gto_tot = 0
          gto(nsp)%nsf_sph_tot = 0          
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
             !write(*,*) 'inode', inode, 'nsp', nsp, 'n_zeta', i,'tmp_nl', tmp_nl, 'l_value', l_value
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
          !          
          !if ( sum(l_nzeta) /= gto(nsp)%n_zeta_tot ) &
          !     call cq_abort('read_gto: problem with number of zeta', gto_unit)
          !
       end if
       ! Broadcast to all the proc.
       call gcopy( gto(nsp)%label, 2)
       call gcopy( gto(nsp)%greatest_angmom )
       call gcopy( gto(nsp)%n_zeta_tot      )
       !
       !write(*,*) 'inode', inode, 'nsp', nsp, 'gto(nsp)%greatest_angmom', gto(nsp)%greatest_angmom
       !write(*,*) 'inode', inode, 'nsp', nsp, 'gto(nsp)%label', gto(nsp)%label
       if ( inode == ionode) write(io_lun,fmt='(10x,"greatest angular momentum =",i4)')&
            gto(nsp)%greatest_angmom       
       !
       allocate( gto(nsp)%angmom( 0:gto(nsp)%greatest_angmom ) )
       !
       if(inode == ionode) then
          do i = 0, gto(nsp)%greatest_angmom
             gto(nsp)%angmom(i)%n_zeta_in_angmom = l_nzeta(i)   
             gto(nsp)%angmom(i)%l_value = i                             
             gto(nsp)%angmom(i)%l_name  = l_name_uppercase(i)
             !write(*,*) 'inode', inode, 'nsp', nsp, 'angmom', i
             !write(*,*) 'inode', inode, '   gto(nsp)%angmom(i)%n_zeta_in_angmom', gto(nsp)%angmom(i)%n_zeta_in_angmom
             !write(*,*) 'inode', inode, '   gto(nsp)%angmom(i)%l_value', gto(nsp)%angmom(i)%l_value
             !write(*,*) 'inode', inode, '   gto(nsp)%angmom(i)%l_name ', gto(nsp)%angmom(i)%l_name
             
          end do
       end if
       !
       do i = 0, gto(nsp)%greatest_angmom
          !
          call gcopy( gto(nsp)%angmom(i)%n_zeta_in_angmom )
          call gcopy( gto(nsp)%angmom(i)%l_value  )
          call gcopy( gto(nsp)%angmom(i)%l_name, 1)
          !
          !write(*,*) 'inode', inode, 'nsp', nsp, 'angmom', i
          !write(*,*) 'inode', inode,'nsp', nsp, 'i', i,  gto(nsp)%angmom(i)%n_zeta_in_angmom
          !write(*,*) 'inode', inode, '   gto(nsp)%angmom(i)%l_value', gto(nsp)%angmom(i)%l_value
          !write(*,*) 'inode', inode, '   gto(nsp)%angmom(i)%l_name ', gto(nsp)%angmom(i)%l_name
          !
          allocate( gto(nsp)%angmom(i)%zeta( gto(nsp)%angmom(i)%n_zeta_in_angmom ) )
          !
       end do

       if ( inode == ionode ) then
          !
          ! rewind and perform reading to grab GTO parameters
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
                gto(nsp)%angmom(l)%zeta(i)%occ  = tmp_occ
                gto(nsp)%angmom(l)%zeta(i)%kind = tmp_kind
                
                tmp_n = tmp_nl(1:1) ; read(tmp_n,*) tmp_n_int
                gto(nsp)%angmom(l)%zeta(i)%n = tmp_n_int                                
                !
             end if
             !
             call gcopy(gto(nsp)%angmom(l)%zeta(i)%ngto)
             call gcopy(gto(nsp)%angmom(l)%zeta(i)%occ )
             call gcopy(gto(nsp)%angmom(l)%zeta(i)%kind)
             call gcopy(gto(nsp)%angmom(l)%zeta(i)%n   )             
             !
             if ( inode == ionode) then
                write(io_lun,fmt='(10x,">> zeta:",x,i2," / state",x,i2,A)') &
                  i, gto(nsp)%angmom(l)%zeta(i)%n, gto(nsp)%angmom(l)%l_name
                write(io_lun,fmt='(10x,"occupation =",x,f4.2)') &
                     gto(nsp)%angmom(l)%zeta(i)%occ
                write(io_lun,fmt='(10x,"kind =",x,a12)') &
                     zeta_kind (gto(nsp)%angmom(l)%zeta(i)%kind )
             end if
             
             allocate( gto(nsp)%angmom(l)%zeta(i)%a ( gto(nsp)%angmom(l)%zeta(i)%ngto) )
             gto(nsp)%angmom(l)%zeta(i)%a = zero
             allocate( gto(nsp)%angmom(l)%zeta(i)%d ( gto(nsp)%angmom(l)%zeta(i)%ngto) )
             gto(nsp)%angmom(l)%zeta(i)%d = zero
             allocate( gto(nsp)%angmom(l)%zeta(i)%c ( gto(nsp)%angmom(l)%zeta(i)%ngto) )
             gto(nsp)%angmom(l)%zeta(i)%c = zero
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
                   !
                   write(io_lun,fmt='(10x,"GTO primitive",x,i3,":",x,"a =",x,f16.10,x,"d =",x,f16.10)') &
                        j, gto(nsp)%angmom(l)%zeta(i)%a(j), gto(nsp)%angmom(l)%zeta(i)%d(j)
                   !
                end do
                !
             end if
             !
             call gcopy( gto(nsp)%angmom(l)%zeta(i)%a, gto(nsp)%angmom(l)%zeta(i)%ngto )
             call gcopy( gto(nsp)%angmom(l)%zeta(i)%d, gto(nsp)%angmom(l)%zeta(i)%ngto )
             call gcopy( gto(nsp)%angmom(l)%zeta(i)%c, gto(nsp)%angmom(l)%zeta(i)%ngto )
             !                          
          end do
          !
          if ( inode == ionode ) then
             ! given angular momentum, compute the number of cartesian Gaussian functions
             gto(nsp)%angmom(l)%nf_gto = (l+1)*(l+2)/2
             ! given angular momentum, accumulate the number of cartesian Gaussian functions
             do i = 1, gto(nsp)%angmom(l)%n_zeta_in_angmom             
                gto(nsp)%nsf_gto_tot = gto(nsp)%nsf_gto_tot + gto(nsp)%angmom(l)%nf_gto
             end do
             !
             ! given angular momentum, compute the number of spherical harmonic functions
             gto(nsp)%angmom(l)%nf_sph = 2*l + 1
             ! given angular momentum, accumulate the number of spherical harmonic functions
             do i = 1, gto(nsp)%angmom(l)%n_zeta_in_angmom             
                gto(nsp)%nsf_sph_tot = gto(nsp)%nsf_sph_tot + gto(nsp)%angmom(l)%nf_sph
             end do
             !
          end if
          !
          call gcopy( gto(nsp)%angmom(l)%nf_gto )
          call gcopy( gto(nsp)%angmom(l)%nf_sph )
          !
          if (inode == ionode) write(io_lun,fmt='(10x,"number of cartesian Gaussian functions =",x,i4)') &
               gto(nsp)%angmom(l)%nf_gto
          !
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
          call gcopy( gto(nsp)%angmom(l)%nx, gto(nsp)%angmom(l)%nf_gto )
          call gcopy( gto(nsp)%angmom(l)%ny, gto(nsp)%angmom(l)%nf_gto )
          call gcopy( gto(nsp)%angmom(l)%nz, gto(nsp)%angmom(l)%nf_gto )
          !
          do j = 1,  gto(nsp)%angmom(l)%nf_gto
             call gcopy( gto(nsp)%angmom(l)%nt(j), 12 )
          end do
          !
       end do
       !
       if (inode == ionode) write(io_lun,fmt='(10x,">> size of the cartesian gaussian basis set =", i4)') &
            gto(nsp)%nsf_gto_tot
       !
       call gcopy( gto(nsp)%nsf_gto_tot )
       call gcopy( gto(nsp)%nsf_sph_tot )
       !
       do l1 = 0, gto(nsp)%greatest_angmom
                         
          allocate( gto(nsp)%angmom(l1)%transform_sph(-l1:+l1) )
          !if (inode == ionode) write(io_lun,*) l1
          if (inode == ionode) write(io_lun,fmt='(10x,"number of cartesian spherical functions =",x,i4)') &
               gto(nsp)%angmom(l1)%nf_sph
          
          do m1 = -l1, +l1
             
             if ( inode == ionode ) then
             
                call sph_gauss_function( inode, ionode, l1, m1, sph_size, sph_nx, sph_ny, sph_nz, sph_coef, sph_name )             
                gto(nsp)%angmom(l1)%transform_sph(m1)%size = sph_size

                !write(*,*) 'nsp, l1, m1, coef 0', nsp, l1, m1, sph_coef
                
             end if
             !
             call gcopy( gto(nsp)%angmom(l1)%transform_sph(m1)%size )
             !
             allocate( gto(nsp)%angmom(l1)%transform_sph(m1)%nx( gto(nsp)%angmom(l1)%transform_sph(m1)%size ) )
             allocate( gto(nsp)%angmom(l1)%transform_sph(m1)%ny( gto(nsp)%angmom(l1)%transform_sph(m1)%size ) )
             allocate( gto(nsp)%angmom(l1)%transform_sph(m1)%nz( gto(nsp)%angmom(l1)%transform_sph(m1)%size ) )
             allocate( gto(nsp)%angmom(l1)%transform_sph(m1)%c ( gto(nsp)%angmom(l1)%transform_sph(m1)%size ) )
             !
             if ( inode == ionode ) then
                !                
                gto(nsp)%angmom(l1)%transform_sph(m1)%nt = sph_name
                !
                gto(nsp)%angmom(l1)%transform_sph(m1)%nx = sph_nx ( 1:gto(nsp)%angmom(l1)%transform_sph(m1)%size )
                gto(nsp)%angmom(l1)%transform_sph(m1)%ny = sph_ny ( 1:gto(nsp)%angmom(l1)%transform_sph(m1)%size )
                gto(nsp)%angmom(l1)%transform_sph(m1)%nz = sph_nz ( 1:gto(nsp)%angmom(l1)%transform_sph(m1)%size )
                gto(nsp)%angmom(l1)%transform_sph(m1)%c  = sph_coef(1:gto(nsp)%angmom(l1)%transform_sph(m1)%size )
                !write(*,*) 'nsp, l1, m1, coef 1', nsp, l1, m1, sph_coef(1:gto(nsp)%angmom(l1)%transform_sph(m1)%size)

                !write(*,*) 'ionode', inode, gto(nsp)%angmom(l1)%transform_sph(m1)%nx, gto(nsp)%angmom(l1)%transform_sph(m1)%ny,&
                !     gto(nsp)%angmom(l1)%transform_sph(m1)%nz, gto(nsp)%angmom(l1)%transform_sph(m1)%c
             end if
             !
             call gcopy( gto(nsp)%angmom(l1)%transform_sph(m1)%nt, 12 )
             !
             call gcopy( gto(nsp)%angmom(l1)%transform_sph(m1)%nx, gto(nsp)%angmom(l1)%transform_sph(m1)%size )
             call gcopy( gto(nsp)%angmom(l1)%transform_sph(m1)%ny, gto(nsp)%angmom(l1)%transform_sph(m1)%size )
             call gcopy( gto(nsp)%angmom(l1)%transform_sph(m1)%nz, gto(nsp)%angmom(l1)%transform_sph(m1)%size )
             !if (inode == ionode) write(*,*) 'nsp, l1, m1, coef 2', nsp, l1, m1, gto(nsp)%angmom(l1)%transform_sph(m1)%c
             call gcopy( gto(nsp)%angmom(l1)%transform_sph(m1)%c,  gto(nsp)%angmom(l1)%transform_sph(m1)%size )

          end do
       end do
       !
       if (inode == ionode) write(io_lun,fmt='(10x,">> size of the cartesian spherical basis set =", i4)') &
            gto(nsp)%nsf_sph_tot
       !
       !#############################
       ! sf derived type for Conquest
       !
       allocate( gto(nsp)%sf( gto(nsp)%nsf_sph_tot ) )
       !
       count = 1
       angu_loop: do l1 = 0, gto(nsp)%greatest_angmom
          
          zeta_loop: do acz1 = 1, gto(nsp)%angmom(l1)%n_zeta_in_angmom
             
             magn_loop: do m1 = -l1, l1

                if ( inode == ionode ) gto(nsp)%sf(count)%sph_size = gto(nsp)%angmom(l1)%transform_sph(m1)%size

                call gcopy( gto(nsp)%sf( count )%sph_size )
                
                allocate(  gto(nsp)%sf(count)%sph_nx  ( gto(nsp)%angmom(l1)%transform_sph(m1)%size ) )
                allocate(  gto(nsp)%sf(count)%sph_ny  ( gto(nsp)%angmom(l1)%transform_sph(m1)%size ) )
                allocate(  gto(nsp)%sf(count)%sph_nz  ( gto(nsp)%angmom(l1)%transform_sph(m1)%size ) )
                allocate(  gto(nsp)%sf(count)%sph_c   ( gto(nsp)%angmom(l1)%transform_sph(m1)%size ) )
                
                sph_loop: do n1 = 1, gto(nsp)%angmom(l1)%transform_sph(m1)%size
                   
                   if ( inode == ionode ) then
                      gto(nsp)%sf( count )%sph_nx(n1) = gto(nsp)%angmom(l1)%transform_sph(m1)%nx( n1 )
                      gto(nsp)%sf( count )%sph_ny(n1) = gto(nsp)%angmom(l1)%transform_sph(m1)%ny( n1 )
                      gto(nsp)%sf( count )%sph_nz(n1) = gto(nsp)%angmom(l1)%transform_sph(m1)%nz( n1 )
                      gto(nsp)%sf( count )%sph_c (n1) = gto(nsp)%angmom(l1)%transform_sph(m1)%c ( n1 )

                      !write(*,*) 'l1, acz1, m1, n1', l1, acz1, m1, n1, 'nx',  gto(nsp)%sf( count )%sph_nx(n1), &
                      !     'ny',  gto(nsp)%sf( count )%sph_ny(n1), 'nz',  gto(nsp)%sf( count )%sph_nz(n1), inode
                      
                   end if
                                      
                   call gcopy( gto(nsp)%sf( count )%sph_nx(n1) )
                   call gcopy( gto(nsp)%sf( count )%sph_ny(n1) )
                   call gcopy( gto(nsp)%sf( count )%sph_nz(n1) )
                   call gcopy( gto(nsp)%sf( count )%sph_c (n1) )

                end do sph_loop

                !if ( inode == 1 ) write(*,*) inode, 'l1, acz1, m1', l1, acz1, m1, 'nx', gto(nsp)%sf( count )%sph_nx 
                !if ( inode == 1 ) write(*,*) inode, 'l1, acz1, m1', l1, acz1, m1, 'ny', gto(nsp)%sf( count )%sph_ny 
                !if ( inode == 1 ) write(*,*) inode, 'l1, acz1, m1', l1, acz1, m1, 'nz', gto(nsp)%sf( count )%sph_nz 

                !if ( inode == 2 ) write(*,*) inode, 'l1, acz1, m1', l1, acz1, m1, 'nx', gto(nsp)%sf( count )%sph_nx 
                !if ( inode == 2 ) write(*,*) inode, 'l1, acz1, m1', l1, acz1, m1, 'ny', gto(nsp)%sf( count )%sph_ny 
                !if ( inode == 2 ) write(*,*) inode, 'l1, acz1, m1', l1, acz1, m1, 'nz', gto(nsp)%sf( count )%sph_nz 

                gto(nsp)%sf( count )%nt   = gto(nsp)%angmom(l1)%transform_sph(m1)%nt
                gto(nsp)%sf( count )%ngto = gto(nsp)%angmom(l1)%zeta(acz1)%ngto
                gto(nsp)%sf( count )%norm = re_cart_norm(l1,m1)                   
                
                call gcopy( gto(nsp)%sf( count )%nt, 12 )
                call gcopy( gto(nsp)%sf( count )%ngto )
                call gcopy( gto(nsp)%sf( count )%norm )
                
                !write(*,*) 'inode', inode, 'l1, acz1, m1', l1, acz1, m1, 'ex, ey, ez', &
                !     gto(nsp)%sf( count )%sph_nx, gto(nsp)%sf( count )%sph_ny, gto(nsp)%sf( count )%sph_nz, &
                !     'nt', gto(nsp)%sf( count )%nt, 'norm', gto(nsp)%sf( count )%norm, 'count', count, &
                !     'ngto', gto(nsp)%sf( count )%ngto
                                
                allocate( gto(nsp)%sf( count )%a( gto(nsp)%sf( count )%ngto ) ) 
                allocate( gto(nsp)%sf( count )%d( gto(nsp)%sf( count )%ngto ) )
                allocate( gto(nsp)%sf( count )%c( gto(nsp)%sf( count )%ngto ) )
                
                primitive: do p = 1, gto(nsp)%angmom(l1)%zeta(acz1)%ngto 

                   if ( inode == ionode ) then
                      
                      gto(nsp)%sf( count )%a( p ) = gto(nsp)%angmom(l1)%zeta(acz1)%a( p )
                      gto(nsp)%sf( count )%d( p ) = gto(nsp)%angmom(l1)%zeta(acz1)%d( p )
                      gto(nsp)%sf( count )%c( p ) = gto(nsp)%angmom(l1)%zeta(acz1)%c( p )
                                            
                   end if
                   
                   !write(*,*) 'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%a( p )',  &
                   !     gto(nsp)%sf( count )%a( p )
                   !write(*,*) 'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%d( p )',  &
                   !     gto(nsp)%sf( count )%d( p )
                   !write(*,*) 'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%c( p )',  &
                   !     gto(nsp)%sf( count )%c( p )

                   
                end do primitive
                call gcopy(  gto(nsp)%sf( count )%a,  gto(nsp)%angmom(l1)%zeta(acz1)%ngto )
                call gcopy(  gto(nsp)%sf( count )%d,  gto(nsp)%angmom(l1)%zeta(acz1)%ngto )
                call gcopy(  gto(nsp)%sf( count )%c,  gto(nsp)%angmom(l1)%zeta(acz1)%ngto )

                !write(*,*) 'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%a( p )',  &
                !     gto(nsp)%sf( count )%a
                !write(*,*) 'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%d( p )',  &
                !     gto(nsp)%sf( count )%d
                !write(*,*) 'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%c( p )',  &
                !     gto(nsp)%sf( count )%c

                count = count + 1
             end do magn_loop
          end do zeta_loop
       end do angu_loop


!!$       allocate( gto(nsp)%sf( gto(nsp)%nsf_gto_tot ) )
!!$       !
!!$       count = 1
!!$       angu_loop: do l1 = 0, gto(nsp)%greatest_angmom
!!$          
!!$          zeta_loop: do acz1 = 1, gto(nsp)%angmom(l1)%n_zeta_in_angmom
!!$             
!!$             magn_loop: do m1 = 1, gto(nsp)%angmom(l1)%nf_gto
!!$
!!$                
!!$                
!!$                if ( inode == ionode ) then
!!$                   gto(nsp)%sf( count )%nx   = gto(nsp)%angmom(l1)%nx( m1 )
!!$                   gto(nsp)%sf( count )%ny   = gto(nsp)%angmom(l1)%ny( m1 )
!!$                   gto(nsp)%sf( count )%nz   = gto(nsp)%angmom(l1)%nz( m1 )
!!$                   gto(nsp)%sf( count )%nt   = gto(nsp)%angmom(l1)%nt( m1 )
!!$                   gto(nsp)%sf( count )%ngto = gto(nsp)%angmom(l1)%zeta(acz1)%ngto                   
!!$                   !gto(nsp)%sf( count )%norm = l_cart_norm(l1)
!!$                   gto(nsp)%sf( count )%norm = re_cart_norm(l1,m1-l1-1)                   
!!$                end if
!!$
!!$                call gcopy( gto(nsp)%sf( count )%nx )
!!$                call gcopy( gto(nsp)%sf( count )%ny )
!!$                call gcopy( gto(nsp)%sf( count )%nz )
!!$                call gcopy( gto(nsp)%sf( count )%nt, 12)
!!$                call gcopy( gto(nsp)%sf( count )%ngto)
!!$                call gcopy( gto(nsp)%sf( count )%norm)
!!$
!!$                print*,'inode', inode, 'l1, acz1, m1', l1, acz1, m1, 'ex, ey, ez', &
!!$                     gto(nsp)%sf( count )%nx, gto(nsp)%sf( count )%ny, gto(nsp)%sf( count )%nz, &
!!$                     'nt', gto(nsp)%sf( count )%nt, 'norm', gto(nsp)%sf( count )%norm, 'count', count, &
!!$                     'ngto', gto(nsp)%sf( count )%ngto
!!$                                
!!$                !if (.not. allocated( gto(nsp)%sf( count )%a ) ) &
!!$                allocate( gto(nsp)%sf( count )%a( gto(nsp)%sf( count )%ngto ) ) 
!!$                !if (.not. allocated( gto(nsp)%sf( count )%d ) ) &
!!$                allocate( gto(nsp)%sf( count )%d( gto(nsp)%sf( count )%ngto ) )
!!$                !if (.not. allocated( gto(nsp)%sf( count )%c ) ) &
!!$                allocate( gto(nsp)%sf( count )%c( gto(nsp)%sf( count )%ngto ) )
!!$                
!!$                primitive: do p = 1, gto(nsp)%angmom(l1)%zeta(acz1)%ngto 
!!$
!!$                   if ( inode == ionode ) then
!!$                      
!!$                      gto(nsp)%sf( count )%a( p ) = gto(nsp)%angmom(l1)%zeta(acz1)%a( p )
!!$                      gto(nsp)%sf( count )%d( p ) = gto(nsp)%angmom(l1)%zeta(acz1)%d( p )
!!$                      gto(nsp)%sf( count )%c( p ) = gto(nsp)%angmom(l1)%zeta(acz1)%c( p )
!!$                                            
!!$                   end if
!!$                   
!!$                   !print*,'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%a( p )',  &
!!$                   !     gto(nsp)%sf( count )%a( p )
!!$                   !print*,'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%d( p )',  &
!!$                   !     gto(nsp)%sf( count )%d( p )
!!$                   !print*,'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%c( p )',  &
!!$                   !     gto(nsp)%sf( count )%c( p )
!!$
!!$                   
!!$                end do primitive
!!$                call gcopy(  gto(nsp)%sf( count )%a,  gto(nsp)%angmom(l1)%zeta(acz1)%ngto )
!!$                call gcopy(  gto(nsp)%sf( count )%d,  gto(nsp)%angmom(l1)%zeta(acz1)%ngto )
!!$                call gcopy(  gto(nsp)%sf( count )%c,  gto(nsp)%angmom(l1)%zeta(acz1)%ngto )
!!$
!!$                print*,'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%a( p )',  &
!!$                     gto(nsp)%sf( count )%a
!!$                print*,'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%d( p )',  &
!!$                     gto(nsp)%sf( count )%d
!!$                print*,'inode', inode, 'nsp', nsp, 'count', count, 'p', p, ' gto(nsp)%sf( count )%c( p )',  &
!!$                     gto(nsp)%sf( count )%c
!!$
!!$                
!!$                count = count + 1
!!$             end do magn_loop
!!$          end do zeta_loop
!!$       end do angu_loop
!!$
!!$       
!!$       deallocate( gto(nsp)%sf )
!!$
!!$       do i = 1, count-1
!!$          deallocate( gto(nsp)%sf( i )%a ) 
!!$          deallocate( gto(nsp)%sf( i )%d )
!!$          deallocate( gto(nsp)%sf( i )%c )
!!$       end do
!!$       
!!$       count2 = 1
!!$       angu_loop2: do l1 = 0, gto(nsp)%greatest_angmom
!!$          
!!$          zeta_loop2: do acz1 = 1, gto(nsp)%angmom(l1)%n_zeta_in_angmom
!!$             
!!$             magn_loop2: do m1 = -l1, +l1
!!$
!!$                sph_loop2: do n1 = 1, gto(nsp)%angmom(l1)%transform_sph(m1)%size
!!$
!!$                   if (inode == ionode ) then
!!$
!!$                      print*, 'l1, acz1, m1, n1', l1, acz1, m1, n1
!!$                   end if
!!$
!!$                   
!!$
!!$                end do sph_loop2
!!$
!!$                
!!$                count2 = count2 + 1
!!$
!!$             end do magn_loop2
!!$          end do zeta_loop2
!!$       end do angu_loop2
!!$       
       !write(*,*)'inode', inode, 'fin'
          
    end do species_loop

    
    call stop_timer(tmr_std_initialisation)
    
    return
  end subroutine read_gto_new
  !
  !  
  !
  subroutine sph_gauss_function( inode, ionode, l, m, size, nx, ny, nz, c, name )

    use GenComms, only: cq_abort

    implicit none

    integer,          intent(in)  :: inode, ionode
    integer,          intent(in)  :: l ! angular momentum
    integer,          intent(in)  :: m ! magnetic moment

    character(len=12),intent(out) :: name
    integer,          intent(out) :: size
    integer,          intent(out), dimension( 5 ) :: nx, ny, nz
    real(double),     intent(out), dimension( 5 ) :: c

    size = 0; nx = 0; ny = 0; nz = 0; c = 0.0d0
    
    if(inode == ionode) then
       !
       select case(l)             
       case(0)
          !
          select case(m)             
          case(0) ! s
             size  = 1
             nx(1) = 0; ny(1) = 0; nz(1) = 0
             c (1) = 1.0d0
             name  = 'S'
             !
             write(io_lun,'(12x,"(nx, ny, nz) = (",I2,x,I2,x,I2,") =",x,A12)') &
                  nx(1), ny(1), nz(1), name           
             !
          case default
             call cq_abort('sph_gauss_function: error angular magnetic moment',l,m)
          end select
          !
       case(1)
          !
          select case(m) 
          case(-1) ! py
             size  = 1
             nx(1) = 0; ny(1) = 1; nz(1) = 0
             c (1) = 1.0d0
             name  = 'P_{y}'
             write(io_lun,'(12x,"(nx, ny, nz) = (",I2,x,I2,x,I2,") =",x,A12)') &
                  nx(1), ny(1), nz(1), name
             !
          case( 0) ! pz
             size  = 1
             nx(1) = 0; ny(1) = 0; nz(1) = 1
             c (1) = 1.0d0
             name  = 'P_{z}'
             write(io_lun,'(12x,"(nx, ny, nz) = (",I2,x,I2,x,I2,") =",x,A12)') &
                  nx(1), ny(1), nz(1), name
             !
          case(+1) ! px
             size  = 1
             nx(1) = 1; ny(1) = 0; nz(1) = 0
             c (1) = 1.0d0
             name  = 'P_{x}'
             write(io_lun,'(12x,"(nx, ny, nz) = (",I2,x,I2,x,I2,") =",x,A12)') &
                  nx(1), ny(1), nz(1), name             
             !
          case default
             call cq_abort('sph_gauss_function: error angular magnetic moment',l,m)
          end select
          !
       case(2)
          !
          select case(m) 
          case(-2) ! d_{xy} = (x*y)
             size  = 1
             nx(1) = 1; ny(1) = 1; nz(1) = 0
             c (1) =-1.0d0 ! phase factor
             name  = 'D_{xy}'
             write(io_lun,'(12x,"(nx, ny, nz) = (",I2,x,I2,x,I2,") =",x,A12)') &
                  nx(1), ny(1), nz(1), name
             !
          case(-1) ! d_{yz} = (y*z)
             size  = 1
             nx(1) = 0; ny(1) = 1; nz(1) = 1
             c (1) = 1.0d0
             name  = 'D_{yz}'
             write(io_lun,'(12x,"(nx, ny, nz) = (",I2,x,I2,x,I2,") =",x,A12)') &
                  nx(1), ny(1), nz(1), name
             !
          case( 0) ! d_{z2} = (3*z*z-r*r) =  (2*z*z-x*x-y*y)
             size  = 3
             nx(1) = 0; ny(1) = 0; nz(1) = 2
             nx(2) = 2; ny(2) = 0; nz(2) = 0
             nx(3) = 0; ny(3) = 2; nz(3) = 0
             c (1) = 2.0d0
             c (2) =-1.0d0
             c (3) =-1.0d0
             name  = 'D_{3z2-r2}'             
             write(io_lun,'(12x,  "(nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4)') &
                  nx(1), ny(1), nz(1), c(1)
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4)') &
                  nx(2), ny(2), nz(2), c(2)
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4,x,"=",x,A12)') &
                  nx(3), ny(3), nz(3), c(3), name
             !
          case(+1) ! d_{xz} = (x*z)
             size  = 1
             nx(1) = 1; ny(1) = 0; nz(1) = 1
             c (1) = 1.0d0
             name  = 'D_{xz}'             
             write(io_lun,'(12x,"(nx, ny, nz) = (",I2,x,I2,x,I2,") =",x,A12)') &
                  nx(1), ny(1), nz(1), name
             !
          case(+2) ! d_{x2-y2} = (x*x-y*y)
             size  = 2
             nx(1) = 2; ny(1) = 0; nz(1) = 0
             nx(2) = 0; ny(2) = 2; nz(2) = 0
             c (1) = 1.0d0
             c (2) =-1.0d0
             name  = 'D_{x2-y2}' 
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4)') &
                  nx(1), ny(1), nz(1), c(1)
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4,x,"=",x,A12)') &
                  nx(2), ny(2), nz(2), c(2), name
             !
          case default
             call cq_abort('sph_gauss_function: error angular magnetic moment',l,m)
          end select
          !
       case(3)
          !
          select case(m) 
          case(-3) ! f_{y(3x2-y2)} = (3x*x*y-y*y*y)
             size  = 2
             nx(1) = 2; ny(1) = 1; nz(1) = 0
             nx(2) = 0; ny(2) = 3; nz(2) = 0             
             c (1) = 3.0d0 
             c (2) =-1.0d0 
             name  = 'F_{y(3x2-y2)}'
             write(io_lun,'(12x,  "(nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4)') &
                  nx(1), ny(1), nz(1), c(1)
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4,x,"=",x,A12)') &
                  nx(3), ny(3), nz(3), c(3), name
             !
          case(-2) ! f_{xyz} = (x*y*z)
             size  = 1
             nx(1) = 1; ny(1) = 1; nz(1) = 1
             c (1) = 1.0d0
             name  = 'D_{xyz}'             
             write(io_lun,'(12x,"(nx, ny, nz) = (",I2,x,I2,x,I2,") =",x,A12)') &
                  nx(1), ny(1), nz(1), name
             !
          case(-1) ! f_{yz2} = (4*z*z*y-y*x*x-y*y*y)
             size  = 3
             nx(1) = 0; ny(1) = 1; nz(1) = 2
             nx(2) = 2; ny(2) = 1; nz(2) = 0
             nx(3) = 0; ny(3) = 3; nz(3) = 0
             c (1) = 4.0d0
             c (2) =-1.0d0
             c (3) =-1.0d0
             name  = 'F_{yz2}'
             write(io_lun,'(12x,  "(nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4)') &
                  nx(1), ny(1), nz(1), c(1)
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4)') &
                  nx(2), ny(2), nz(2), c(2)
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4,x,"=",x,A12)') &
                  nx(3), ny(3), nz(3), c(3), name
             !
          case( 0) ! f_{z3} = (2*z*z*z-3*x*x*z-3*y*y*z) 
             size  = 3
             nx(1) = 0; ny(1) = 0; nz(1) = 3
             nx(2) = 2; ny(2) = 0; nz(2) = 1
             nx(3) = 0; ny(3) = 2; nz(3) = 1
             c (1) = 2.0d0
             c (2) =-3.0d0
             c (3) =-3.0d0
             name  = 'F_{z3}'             
             write(io_lun,'(12x,  "(nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4)') &
                  nx(1), ny(1), nz(1), c(1)
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4)') &
                  nx(2), ny(2), nz(2), c(2)
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4,x,"=",x,A12)') &
                  nx(3), ny(3), nz(3), c(3), name
             !
          case(+1) ! f_{xz2} = (4*x*z*z-x*x*x-x*y*y)
             size  = 3
             nx(1) = 1; ny(1) = 0; nz(1) = 2
             nx(2) = 3; ny(2) = 0; nz(2) = 0
             nx(3) = 1; ny(3) = 2; nz(3) = 0
             c (1) = 4.0d0
             c (2) =-1.0d0
             c (3) =-1.0d0
             name  = 'F_{xz2}'
             write(io_lun,'(12x,  "(nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4)') &
                  nx(1), ny(1), nz(1), c(1)
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4)') &
                  nx(2), ny(2), nz(2), c(2)
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4,x,"=",x,A12)') &
                  nx(3), ny(3), nz(3), c(3), name
                 !
          case(+2) ! f_{z(x2-y2)} = (x*x*z-y*y*z)
             size  = 2
             nx(1) = 2; ny(1) = 0; nz(1) = 1
             nx(2) = 0; ny(2) = 2; nz(2) = 1
             c (1) = 1.0d0
             c (2) =-1.0d0
             name  = 'F_{z(x2-y2)}'
             write(io_lun,'(12x,  "(nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4)') &
                  nx(1), ny(1), nz(1), c(1)
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4,x,"=",x,A12)') &
                  nx(2), ny(2), nz(2), c(2), name
             !
          case(+3) ! f_{x(x2-3y2)} = (x*x*x-3*x*y*y)
             size  = 2
             nx(1) = 3; ny(1) = 0; nz(1) = 0
             nx(2) = 1; ny(2) = 2; nz(2) = 0
             c (1) = 1.0d0
             c (2) =-3.0d0
             name  = 'F_{x(x2-3y2}'
             write(io_lun,'(12x,  "(nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4)') &
                  nx(1), ny(1), nz(1), c(1)
             write(io_lun,'(10x,"+ (nx, ny, nz) = (",I2,x,I2,x,I2,")",x,"x",x,f8.4,x,"=",x,A12)') &
                  nx(2), ny(2), nz(2), c(2), name

             !
          case default
             call cq_abort('sph_gauss_function: error angular magnetic moment',l,m)
          end select
          !
       case default
          call cq_abort('sph_gauss_function: error angular momentum',l) 
       end select
       !
    end if
    !
  end subroutine sph_gauss_function
  !
  !
  !
  subroutine cart_gauss_function( inode, ionode, lambda, name, a, nx, ny, nz, nt )

    implicit none

    integer,          intent(in) :: lambda
    integer,          intent(in) :: inode, ionode
    character(len=1), intent(in) :: name
    real(double),     intent(in) :: a

    integer,           intent(out), dimension( (lambda+1)*(lambda+2)/2 ) :: nx, ny, nz
    character(len=12), intent(out), dimension( (lambda+1)*(lambda+2)/2 ) :: nt
    character(len=12) :: charx, chary, charz
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
             nt(index) = trim(name)//'_{'//trim(charx)//trim(chary)//trim(charz)//'}'
             !
             !norm(index)= norm_gto( &
             !     a,     & ! alpha
             !     nx(index), & ! l
             !     ny(index), & ! m
             !     nz(index))   ! n
             !
             write(io_lun,'(12x,"(nx, ny, nz) = (",I2,x,I2,x,I2,") =",x,A12)') &
                  nx(index), &
                  ny(index), &
                  nz(index), &
                  nt(index)
             !                      
             index = index + 1
             !
          end do
       end do
    end if
    !
  end subroutine cart_gauss_function
  !
  ! 
  !
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


!!$  subroutine read_gto(inode,ionode,n_species)
!!$
!!$    use gto_format,     ONLY: gto
!!$    use species_module, ONLY: gto_file, species_label
!!$    use GenComms,       ONLY: cq_abort, gcopy, my_barrier
!!$    use exx_erigto,     ONLY: norm_gto, overlap_gto
!!$
!!$    implicit none
!!$
!!$    ! Arguments
!!$    integer, intent(in)  :: inode, ionode, n_species
!!$
!!$    ! Local variables 
!!$    integer :: ios, match, tmp_int, gto_unit, nsp
!!$    character(len=256) :: match_atom
!!$    character(len=8)   :: charx, chary, charz
!!$
!!$    integer :: i, j, k, l, lambda
!!$    integer :: x, y, z, index, count
!!$    integer :: l1, m1, n1, l2, m2, n2
!!$    real(double) :: x1, y1, z1    
!!$    real(double) :: x2, y2, z2    
!!$    real(double) :: test, norm, S_12, alpha1, alpha2
!!$
!!$    !
!!$    ! gto_file(nsp) allocated in initial_read.module
!!$    ! gto(nsp)      allocated in pseudo_tm_info
!!$    !
!!$
!!$    gto_unit = 1000
!!$
!!$    call start_timer(tmr_std_initialisation)
!!$
!!$    species_loop: do nsp = 1, n_species
!!$
!!$       if(inode == ionode) then   
!!$
!!$          gto_file(nsp)   = trim(species_label(nsp))//trim('.gto')
!!$          gto(nsp)%label  = trim(species_label(nsp))
!!$          !
!!$          open(unit=gto_unit,file=gto_file(nsp),status="old",action="read", &
!!$               iostat=ios,position="rewind")  
!!$          if(ios /= 0) call cq_abort('read_gto: failed to open input file')
!!$          !
!!$          !call gcopy(gto(nsp)%label,2)     
!!$          !
!!$          match = 0
!!$          count = 1
!!$          !
!!$          do while (match == 0)
!!$             read(unit=gto_unit,fmt=*,iostat=ios) match_atom
!!$             match_atom = trim(match_atom)
!!$             if (ios < 0) call cq_abort('read_gto: end error GTO data file', gto_unit)        
!!$             if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)
!!$             if (match_atom == gto(nsp)%label) then
!!$                match = 1
!!$             end if
!!$          end do
!!$          !
!!$          read(gto_unit,fmt=*,iostat=ios) gto(nsp)%nshell
!!$          if (ios < 0) call cq_abort('read_gto: end error GTO data file',gto_unit)        
!!$          if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)
!!$          !
!!$          write(*,*)
!!$          write(*,fmt='("#### read GTO for ",A2,"/nshell =",I3," ####")') gto(nsp)%label, gto(nsp)%nshell
!!$          !
!!$          !write(*,fmt='(2x,"nshell = ",I3)') gto(nsp)%nshell
!!$          !
!!$       end if
!!$       !
!!$       call gcopy(gto(nsp)%label,2)
!!$       call gcopy(gto(nsp)%nshell)
!!$       !
!!$       allocate(gto(nsp)%shell(gto(nsp)%nshell))
!!$       !print*, 'node', inode, allocated(gto(nsp)%shell), size(gto(nsp)%shell)
!!$       !  
!!$       !
!!$       gto(nsp)%nprim = 0
!!$       gto(nsp)%nang  = 0
!!$       !
!!$       shell_loop: do i = 1, gto(nsp)%nshell
!!$          !
!!$          if(inode == ionode) then
!!$             !
!!$             read(gto_unit,fmt=*,iostat=ios) &
!!$                  gto(nsp)%shell(i)%t, & ! type s,p,sp,d,f...
!!$                  gto(nsp)%shell(i)%n, & ! number of primitive
!!$                  gto(nsp)%shell(i)%f    ! scaling factor
!!$             if (ios < 0) call cq_abort('read_gto: end error GTO data file', gto_unit)        
!!$             if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)             
!!$             !
!!$             write(*,'(">>> shell ",I3)') i
!!$             write(*,'(2x,"type           =",3x,A2)'   ) gto(nsp)%shell(i)%t
!!$             write(*,'(2x,"scaling factor =",3x,F8.6)' ) gto(nsp)%shell(i)%f
!!$
!!$             select case(gto(nsp)%shell(i)%t)
!!$             case('S')
!!$                lambda = 0
!!$                !gto(nsp)%shell(i)%nn = 0
!!$                !gto(nsp)%shell(i)%nf = 1
!!$             case('P')
!!$                lambda = 1
!!$                !gto(nsp)%shell(i)%nn = 1
!!$                !gto(nsp)%shell(i)%nf = 3
!!$             case('D')
!!$                lambda = 2
!!$                !gto(nsp)%shell(i)%nn = 2
!!$                !gto(nsp)%shell(i)%nf = 6
!!$             case('F')
!!$                lambda = 3
!!$                !gto(nsp)%shell(i)%nn = 3
!!$                !gto(nsp)%shell(i)%nf = 10
!!$             case('G')
!!$                lambda = 4
!!$                !gto(nsp)%shell(i)%nn = 4
!!$                !gto(nsp)%shell(i)%nf = 15
!!$             case('H')
!!$                lambda = 5
!!$                !gto(nsp)%shell(i)%nn = 5
!!$                !gto(nsp)%shell(i)%nf = 24
!!$             case default
!!$                call cq_abort('read_gto: error angular momentum',gto(nsp)%shell(i)%nn) 
!!$             end select
!!$             !
!!$             gto(nsp)%shell(i)%nn     = lambda
!!$             gto(nsp)%shell(i)%nf     = (lambda+1)*(lambda+2)/2
!!$             gto(nsp)%nang            = gto(nsp)%nang  + gto(nsp)%shell(i)%nf
!!$             gto(nsp)%nprim           = gto(nsp)%nprim + gto(nsp)%shell(i)%n*gto(nsp)%shell(i)%nf
!!$             !
!!$             write(*,'(2x,"angular momentum     l =",I3)') gto(nsp)%shell(i)%nn  
!!$             write(*,'(2x,"number of values for m =",I3)') gto(nsp)%shell(i)%nf              
!!$             write(*,'(2x,"number of primitives   =",I3)') gto(nsp)%shell(i)%n              
!!$             !
!!$             !
!!$          end if
!!$          !
!!$          call gcopy(gto(nsp)%shell(i)%t,2)     
!!$          call gcopy(gto(nsp)%shell(i)%n  )     
!!$          call gcopy(gto(nsp)%shell(i)%f  )     
!!$          call gcopy(gto(nsp)%shell(i)%nn )     
!!$          call gcopy(gto(nsp)%shell(i)%nf )
!!$          !
!!$          allocate(gto(nsp)%shell(i)%ak( gto(nsp)%shell(i)%n ))
!!$          allocate(gto(nsp)%shell(i)%dk( gto(nsp)%shell(i)%n ))                 
!!$          allocate(gto(nsp)%shell(i)%Nk( gto(nsp)%shell(i)%n,gto(nsp)%shell(i)%nf ))
!!$          !
!!$          allocate(gto(nsp)%shell(i)%nx( gto(nsp)%shell(i)%nf ))
!!$          allocate(gto(nsp)%shell(i)%ny( gto(nsp)%shell(i)%nf ))                 
!!$          allocate(gto(nsp)%shell(i)%nz( gto(nsp)%shell(i)%nf ))                 
!!$          allocate(gto(nsp)%shell(i)%nt( gto(nsp)%shell(i)%nf ))
!!$          allocate(gto(nsp)%shell(i)%norm(  gto(nsp)%shell(i)%nf ))  
!!$          !
!!$          prim_loop: do j = 1, gto(nsp)%shell(i)%n 
!!$             !
!!$             if(inode == ionode) then
!!$
!!$                read(gto_unit,fmt=*,iostat=ios) &
!!$                     gto(nsp)%shell(i)%ak(j),   &
!!$                     gto(nsp)%shell(i)%dk(j)
!!$                if (ios < 0) call cq_abort('read_gto: end error GTO data file', gto_unit)        
!!$                if (ios > 0) call cq_abort('read_gto: read error GTO data file',gto_unit)
!!$
!!$                write(*,'(2x,"primitive",x,"(a",I3,", d",I3,") = (",F14.8,x,F14.8,")")') j, j,  &
!!$                     gto(nsp)%shell(i)%ak(j), &
!!$                     gto(nsp)%shell(i)%dk(j)                
!!$
!!$                !write(*,'(6x,"(a",I3,", d",I3,") = (",F14.8,x,F14.8,")")') j, j, &
!!$                !gto(nsp)%shell(i)%ak(j), &
!!$                !gto(nsp)%shell(i)%dk(j)
!!$             end if
!!$             !
!!$             call gcopy(gto(nsp)%shell(i)%ak(j))
!!$             call gcopy(gto(nsp)%shell(i)%dk(j))
!!$             !
!!$             index  = 1
!!$             lambda = gto(nsp)%shell(i)%nn 
!!$             do x = 0, lambda
!!$                do y = 0, lambda - x
!!$                   z = lambda - (x + y)
!!$                   !
!!$                   if(inode == ionode) then
!!$                      ! Exponents for Cartesian Gaussian functions
!!$                      gto(nsp)%shell(i)%nx(index) = x!x
!!$                      gto(nsp)%shell(i)%ny(index) = y!y
!!$                      gto(nsp)%shell(i)%nz(index) = z!z
!!$                      !
!!$                   end if
!!$                   !
!!$                   charx = '' ; chary = '' ; charz = '' 
!!$                   if (x > 0) then
!!$                      do  k = 1, x
!!$                         charx = trim(charx)//'x' 
!!$                      end do
!!$                   end if
!!$                   if (y > 0) then
!!$                      do  k = 1, y
!!$                         chary = trim(chary)//'y' 
!!$                      end do
!!$                   end if
!!$                   if (z > 0) then
!!$                      do  k = 1, z
!!$                         charz = trim(charz)//'z' 
!!$                      end do
!!$                   end if
!!$                   !
!!$                   if(inode == ionode) then
!!$                      !
!!$                      gto(nsp)%shell(i)%nt(index) = &
!!$                           trim(gto(nsp)%shell(i)%t)//trim(charx)//trim(chary)//trim(charz)
!!$                      !
!!$                   end if
!!$                   !
!!$                   call gcopy(gto(nsp)%shell(i)%nt(index),16)
!!$                   !
!!$                   if(inode == ionode) then                   
!!$                      ! Compute the norm of each Cartesian components
!!$
!!$                      gto(nsp)%shell(i)%Nk(j,index) = norm_gto( &
!!$                           gto(nsp)%shell(i)%ak(j),     & ! alpha
!!$                           gto(nsp)%shell(i)%nx(index), & ! l
!!$                           gto(nsp)%shell(i)%ny(index), & ! m
!!$                           gto(nsp)%shell(i)%nz(index))   ! n
!!$                      !
!!$                      write(*,'(6x,"(nx, ny, nz) = (",I2,x,I2,x,I2,") =",x,A6,"/","norm = ",F14.8,I8)') &
!!$                           gto(nsp)%shell(i)%nx(index), &
!!$                           gto(nsp)%shell(i)%ny(index), &
!!$                           gto(nsp)%shell(i)%nz(index), &
!!$                           gto(nsp)%shell(i)%nt(index), &
!!$                           gto(nsp)%shell(i)%Nk(j,index), index
!!$
!!$                   end if
!!$                   !
!!$                   call gcopy(gto(nsp)%shell(i)%Nk(j,index))
!!$                   !                      
!!$                   index = index + 1
!!$                   count = count + 1
!!$                   !
!!$                end do
!!$             end do
!!$             !
!!$          end do prim_loop
!!$          !
!!$          if(inode == ionode) then
!!$             !
!!$             gto(nsp)%shell(i)%norm(:) = zero
!!$             !
!!$             do j = 1, gto(nsp)%shell(i)%n
!!$                !
!!$                do k = 1, gto(nsp)%shell(i)%n
!!$
!!$                   do l = 1, gto(nsp)%shell(i)%nf
!!$
!!$                      alpha1 = gto(nsp)%shell(i)%ak(j)
!!$                      l1     = gto(nsp)%shell(i)%nx(l)
!!$                      m1     = gto(nsp)%shell(i)%ny(l)
!!$                      n1     = gto(nsp)%shell(i)%nz(l)
!!$                      x1 = zero ; y1 = zero ; z1 = zero
!!$
!!$                      alpha2 = gto(nsp)%shell(i)%ak(k)
!!$                      l2     = gto(nsp)%shell(i)%nx(l)
!!$                      m2     = gto(nsp)%shell(i)%ny(l)
!!$                      n2     = gto(nsp)%shell(i)%nz(l)
!!$                      x2 = zero ; y2 = zero ; z2 = zero
!!$
!!$                      S_12   = overlap_gto(alpha1,l1,m1,n1,x1,y1,z1, &
!!$                           alpha2,l2,m2,n2,x2,y2,z2)
!!$
!!$                      !print*, j, k, l1, m1, n1, l2, m2, n2, alpha1, alpha2, l, S_12
!!$
!!$                      gto(nsp)%shell(i)%norm(l) = gto(nsp)%shell(i)%norm(l) + &
!!$                           gto(nsp)%shell(i)%ak(j) * gto(nsp)%shell(i)%ak(k) * S_12 
!!$                      !
!!$                   end do
!!$                   !
!!$                end do
!!$                !
!!$             end do
!!$             !
!!$          end if
!!$          !
!!$          if(inode == ionode) then
!!$             !
!!$             write(*,'(2x,"number of angular funtions = ",I3)') gto(nsp)%shell(i)%nf
!!$             !
!!$             do l = 1, gto(nsp)%shell(i)%nf
!!$                write(*,'(2x,"norm =",F12.6)') gto(nsp)%shell(i)%norm(l)
!!$             end do
!!$             !
!!$          end if
!!$          !
!!$          call gcopy( gto(nsp)%shell(i)%norm, gto(nsp)%shell(i)%nf)
!!$          !
!!$          !
!!$       end do shell_loop
!!$
!!$       if(inode == ionode) then
!!$          !
!!$          write(*,'(">>> total number of angular functions =",I10)') gto(nsp)%nang
!!$          write(*,'(">>> total number of primitives        =",I10)') gto(nsp)%nprim
!!$          !
!!$          close(unit=gto_unit,iostat=ios)  
!!$          if(ios /= 0) call cq_abort('read_gto: failed to close input file',gto_unit)
!!$       end if
!!$
!!$       call gcopy(gto(nsp)%nang )
!!$       call gcopy(gto(nsp)%nprim)
!!$
!!$    end do species_loop
!!$
!!$    call stop_timer(tmr_std_initialisation)
!!$
!!$    return
!!$  end subroutine read_gto
  !  
end module read_gto_info
