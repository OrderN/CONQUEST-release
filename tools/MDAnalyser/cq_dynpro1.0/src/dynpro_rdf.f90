module dynpro_rdf

  use kinds, only      : DP
  use constants,  only : pi, eps12
  use constants,  only : ps   => AU_PS
  use constants,  only : c    => C_SI
  use constants,  only : kb   => K_BOLTZMANN_SI
  use constants,  only : hp   => H_PLANCK_SI
  use constants,  only : bohr => BOHR_RADIUS_ANGS
  use constants,  only : ry   => AUTOEV

  use dynpro_cell, only : center_cell
  use dynpro_math, only : vec_product
  
  implicit none

contains

  subroutine dNab (nsp, natoms,  nframes, nr, rmin, rmax,   & 
       ityp, syst_ID_spe, cel, abc,  &
       rd_atom, rd_count, dN, event, V)
    !==-----------------------------------------------------------------------==
    !    dN_{ab}(r) = "number of particle pair probability" 
    !==-----------------------------------------------------------------------==
    !
    ! given 2 definite species a and b separated by a distance ranging
    ! in [r, r+dr], dN_{ab} is number of pair a-b in this range
    !
    !==-----------------------------------------------------------------------==

    implicit none

    integer,  intent(in)        :: nsp, natoms, rd_atom, rd_count, nr
    integer,  intent(in)        :: ityp(natoms), syst_ID_spe(nsp,2)

    integer,  intent(in)        :: nframes
    real(DP), intent(in)        :: rmin, rmax

    real(DP), intent(in)        :: cel(3,3*nframes)
    real(DP), intent(in)        :: abc(3,natoms*nframes)

    integer,  intent(out)       :: dN(nframes,rd_count,nsp,nr)
    integer,  intent(out)       :: event(nframes,rd_count,nsp)
    real(DP)                    :: vol(nframes), V

    integer                     :: i, j, k, l, m, n, o, p, q
    integer                     :: atom_ID_nbr(rd_count)
    real(DP)                    :: rad, r, rstep
    real(DP)                    :: atom_ID_abc(3,natoms)
    real(DP)                    :: atom_ID_tmp(3,rd_count)
    real(DP)                    :: cell_ID_tmp(3,3)
    real(DP)                    :: a(3), b(3), c(3), ab(3)
    real(DP)                    :: ct_abc(3,natoms), ct_xyz(3,natoms)

    o     = 0
    p     = 0
    event = 0
    dN    = 0
    r     = 0.0d0
    a     = 0.0d0    
    b     = 0.0d0
    c     = 0.0d0
    ab    = 0.0d0
    rstep = (rmax - rmin) / dble(nr - 1)
    !
    j = 1
    do i = 1, natoms     
       if (ityp(i) == real(rd_atom)) then
          atom_ID_nbr(j) = i 
          j = j + 1
       end if
    end do
    !
    frames_loop: do q = 1, nframes
       !
       atom_ID_abc(:,:) = abc(:,(q+p):(q*natoms))
       cell_ID_tmp(:,:) = cel(:,(q+o):(q*3))
       a(:)             = cel(1,(q+o):(q*3))
       b(:)             = cel(2,(q+o):(q*3))
       c(:)             = cel(3,(q+o):(q*3))
       !
       call vec_product(a, b, ab)
       vol(q)           = DOT_product(ab,c)
       !
       rd_atom_loop: do i = 1, rd_count
          m  = 1      
          !
          call center_cell(natoms, 1, atom_ID_nbr(i), cell_ID_tmp, &
               atom_ID_abc, ct_abc, ct_xyz)
          !
          j = 1
          do k = 1, natoms     
             if (ityp(k) == real(rd_atom)) then
                atom_ID_tmp(1:3,j) = ct_xyz(1:3,k)       
                j = j + 1
             end if
          end do
          !          
          species_loop: do j = 1, nsp     
             !
             ext_atom_loop: do k = 1, syst_ID_spe(j,2)           
                !
                r = sqrt((atom_ID_tmp(1,i) - ct_xyz(1,m))**2 &                
                     +   (atom_ID_tmp(2,i) - ct_xyz(2,m))**2 &
                     +   (atom_ID_tmp(3,i) - ct_xyz(3,m))**2)
                !
                m   = m + 1
                rad = rmin
                !
                rad_point_loop: do l = 1, nr   
                   !                   
                   if (r <= rad .and. r /= 0.0d0) then
                      dN(q,i,j,l)  = dN(q,i,j,l)  + 1      
                      event(q,i,j) = event(q,i,j) + 1
                      exit 
                   end if
                   !
                   rad = rad + rstep
                   !

                end do rad_point_loop
                !
             end do ext_atom_loop
             !
          end do species_loop
          !
       end do rd_atom_loop
       !
       o = o + 2
       p = p + (natoms - 1)
       !
    end do frames_loop
    !
    V = sum(vol)/nframes
    !
    return
  end subroutine DNAB


  subroutine RDF (dN, event, count, rd_atom, rd_count, rd_lab, &
       natoms, nsp, nr, rmin, rmax, rint, V, syst_ID_spe, syst_ID_lab, &
       filerdf, lau_unit)
    
    !==-----------------------------------------------------------------------==
    !    RDF = g_{ab}(r) = "pair correlation function" (Radial Distribution)
    !==-----------------------------------------------------------------------==
    !
    !  g_{ab}(r) = (V/N_{ab))<rho_{ab}(r)>
    !
    !  where
    !        <rho_{ab}(r)> : the local density of (pair) particle
    !                      = SUM_{a,b,a/=b}(1/(4*pi*r**2*dr))*dN_{ab}
    !                           
    !        V             : volume of the cell
    !        N_{ab}        : number of pair of particles
    !                      = Na*Nb      if a /= b
    !                      = Na(Na - 1) if a == b 
    !  Remark:
    !        (N_{ab}/V)*int[dr g_{ab}(r)*4*pi*r**2] 
    !                      = Na*Nb       if a /= b
    !                      = Na*(Na - 1) if a == b
    !==-----------------------------------------------------------------------==

    implicit none

    logical,  intent(in)        :: lau_unit    
    integer,  intent(in)        :: natoms, nsp, nr, count, rd_atom, rd_count
    integer,  intent(in)        :: dN(count,rd_count,nsp,nr)
    integer,  intent(in)        :: event(count,rd_count,nsp)
    integer,  intent(in)        :: syst_ID_spe(nsp,2)
    real(DP), intent(in)        :: rmin, rmax
    real(DP), intent(inout)     :: rint(100)
    real(DP), intent(in)        :: V

    character(256), intent(in)  :: filerdf
    character(2),   intent(in)  :: syst_ID_lab(nsp), rd_lab

    character(4)                :: r_unit
    character(256)              :: rdfout(nsp),    int_g_out(nsp)
    real(DP)                    :: g_int_max(nsp), int_g_max(nsp), int_g(nsp)
    real(DP)                    :: g_tmp(nsp,nr) 
    real(DP)                    :: rad, rstep, Nab, scal
    !real(DP)                    :: g_sum(nsp)  

    integer                     :: nr_int(100)
    integer                     :: event_tmp(count,rd_count,nsp)
    integer                     :: i, j, k, l, m, n

    !
    ! Setup output and output file names
    do i = 1, nsp
       rdfout(i)    = trim(filerdf)//'.'//trim(rd_lab)//'_'//trim(syst_ID_lab(i))//'.dat'
       int_g_out(i) = trim(rd_lab)//'_'//trim(syst_ID_lab(i))
       !print*, ' syst_ID_lab', syst_ID_lab(i), syst_ID_spe(i,1),  syst_ID_spe(i,2) 
       !print*, rdfout(i),  int_g_out(i)
    end do
    !
    ! Compute integration step
    rstep = (rmax - rmin) / dble(nr - 1)
    !
    ! Open RDF output files
    if(lau_unit) then
       r_unit = 'Bohr'
    else
       r_unit = 'Ang.'
    end if
    !
    do m = 1, nsp
       open((1000+m), file=rdfout(m), status='unknown')
       write(1000+m,'(a5,a4,a1,10x,a20,4x,a21)') '# r [',r_unit,']','g_ab(r) [Normalized]','g_int(r) [Normalized]'
    end do
    !
    ! Setup nr_int for integration
    nr_int = 0
    do m = 1, nsp
       ! 
       if ( rint(m) < rmin ) then
          rint(m) = rmin
       end if
       !
       nr_int(m) = ceiling( (rint(m) - rmin)/rstep )
       !
       print*, 'nsp/r_int/nr_int:', m, rint(m), nr_int(m), (rint(m) - rmin)/rstep 
       !
    end do
    !
    g_int_max = 0.0d0
    int_g_max = 0.0d0
    int_g = 0.0d0
    !
    g_tmp = 0.0d0
    !
    Nab   = 0.0d0
    !
    scal  = 1.0d0
    !
    ! Compute RDF for each species
    species_loop: do i = 1, nsp
       !
       rad = rmin
       Nab = rd_count*syst_ID_spe(i,2)
       !
       rad_point_loop: do j = 1, nr
          !
          g_tmp(i,j)   = sum(dble(dN(:,:,i,j)))/(count*4*pi*rad**2*rstep*Nab/V)
          !
          g_int_max(i) = g_int_max(i) + (4*pi*rad**2)*rstep*Nab/V*g_tmp(i,j)
          !
          write((1000+i),'(f10.6,1X,2(f24.8))') rad, g_tmp(i,j), g_int_max(i) / rd_count !syst_ID_spe(i,2)
          !
          rad = rad + rstep 
          !
          if(j == nr_int(i)) then
             int_g(i) = g_int_max(i) / rd_count !/syst_ID_spe(i,2)
             !
          end if
          !
       end do rad_point_loop
       !
       do m = 1, nsp
          int_g_max(m) = g_int_max(m) / syst_ID_spe(m,2)
       end do
       !
       !g_sum(i) = sum(g_tmp(i,:)) 
       !
    end do species_loop
    !
    write(*,*) '====== RDF parameters ====================================>>'
    write(*,'(a22,2x,a10)')   'rd_atom', rd_lab
    write(*,'(a22,2x,i10)')   'rd_count', rd_count
    write(*,'(a22,2x,i10)')   'nb. of species', nsp
    write(*,'(a22,2x,i10)')   'nb. of atoms', natoms
    write(*,'(a22,2x,i10)')
    write(*,'(a22,2x,i10)')   'nb. radial point', nr
    write(*,'(a22,2x,f10.4)') 'r_min (Ang.)', rmin
    write(*,'(a22,2x,f10.4)') 'r_max (Ang.)', rmax
    write(*,'(a22,2x,f10.4)') 'r_step (Ang.)', rstep
    write(*,'(a22,2x,i10)')   'nframes', count
    write(*,'(a22,2x,i10)')
    write(*,'(a22)')         'output file(s)'
    !
    do i = 1, nsp
       write(*,'(20x,a20)') rdfout(i)
    end do
    !
    write(*,'(a22)')         'int[g_{ab}(r)dr]'
    !
    do i = 1, nsp
       write(*,'(20x,a5,2x,f8.2,2x,f12.5)') int_g_out(i), int_g_max(i), int_g(i)
    end do
    !
    write(*,*) '==========================================================<<'
    !
    do m = 1, nsp
       close(1000+m)
    end do
    !
    return

  end subroutine RDF

end module dynpro_rdf
