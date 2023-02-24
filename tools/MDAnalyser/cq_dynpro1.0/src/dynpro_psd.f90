module dynpro_psd

  use kinds, only      : DP
  use constants,  only : pi, eps12
  use constants,  only : ps   => AU_PS
  use constants,  only : c    => C_SI
  use constants,  only : kb   => K_BOLTZMANN_SI
  use constants,  only : hp   => H_PLANCK_SI
  use constants,  only : bohr => BOHR_RADIUS_ANGS
  use constants,  only : ry   => AUTOEV

  use dynpro_fft, only : EVLMEM, MEMCOF, SPCTRM
  
  implicit none

contains

  subroutine PSD(type, nsp, lspe, iteration, syst_ID_lab, syst_ID_spe, &
       filevac, loverlap, n_seg, win, unit_vac, m_cof, time_step, & 
       interval_i_m, iprint, unit, file, qcfactor, temp)

    implicit none

    integer,        intent(in)        :: iprint, unit, nsp, type
    real(DP),       intent(in)        :: time_step
    integer,        intent(in)        :: iteration, interval_i_m
    integer,        intent(in)        :: syst_ID_spe(nsp,2)
    character(2),   intent(in)        :: syst_ID_lab(nsp)

    integer,        intent(in)        :: n_seg, unit_vac ! SPC
    character(256), intent(in)        :: file, filevac
    logical,        intent(in)        :: loverlap, lspe

    integer,        intent(inout)     :: m_cof           ! MEM

    integer,        intent(in)        :: qcfactor        ! QCF
    real(DP),       intent(in)        :: temp

    real(DP),       allocatable       :: C_tot(:)

    character(18)                     :: qcf_name
    character(256)                    :: memfile, spcfile, vacfile
    character(256)                    :: vacfilespe(nsp), spcfilespe(nsp)

    real(DP),        allocatable      :: J(:),  omega(:), Q(:) 
    real(DP),        allocatable      :: w1(:), w2(:), wm(:), cof(:), p(:), fdt(:)
    real(DP)                          :: a0, dt_mem         
    real(DP)                          :: dt_spc

    character(11)                     :: win_name

    real(DP)                          :: time_conv, freq_conv, ps_to_cm, norm ! conversion

    integer                           :: i, ii, m, win, mod_m, ios

    ! 1 cm^-1 = 2.99792458E10  s^-1

    ps_to_cm  = c*eps12*100.0d0
    time_conv = time_step*iprint*interval_i_m*ps ! ps
    freq_conv = (1.0d0/time_conv)/ps_to_cm       ! 1/dt 

    print*, '   ps_to_cm :', ps_to_cm
    print*, '   time_conv:', time_conv
    print*, '   freq_conv:', freq_conv
    print*, '   time_step:', time_step
    print*, '      iprint:', iprint
    print*, 'interval_i_m:', interval_i_m
    print*, '          ps:', ps


    !print*, type,  nsp, lspe, iteration, syst_ID_lab, syst_ID_spe, &
    !   filevac, loverlap, n_seg, win, unit_vac, m_cof, time_step,  & 
    !   interval_i_m, iprint, unit, file, qcfactor, temp

    !SELECT CASE(type)
    !CASE(0)
    !   vacfile = TRIM(filevac)//'1.dat'
    !CASE(1)   
    vacfile = trim(filevac)//'_raw.dat'

    !END SELECT
    memfile = trim(file)//'_mem.dat'
    spcfile = trim(file)//'_spc.dat'

    !--------------------------------------------------------------------

    print*, 'M_COF= ', m_cof
    if (m_cof == 0) then
       m_cof  = iteration/10
    else
       if (m_cof >= iteration) then
          write(*,*)'Error -PSD- check the m_cof and nframes values !'
          write(*,*)'m_cof =', m_cof, 'sp_n =', iteration
          stop
       end if
    end if

    allocate(fdt(0:m_cof))
    allocate(J(m_cof))
    allocate(omega(m_cof))
    allocate(Q(m_cof))
    allocate(w1(iteration))
    allocate(w2(iteration))  
    allocate(wm(m_cof))
    allocate(cof(m_cof))
    !
    allocate(C_tot(iteration))
    !
    open(unit_vac, action="read", iostat=ios, file=vacfile, status='old', position='rewind')
    !
    if (ios /= 0) then
       print *,"Error -DynPro PSD- opening files:  ", vacfile
       stop
    end if
    !
    print*, ios
    do i = 1, iteration
       read(unit_vac,'(f20.10)') C_tot(i)
    end do
    !
    close(unit_vac)
    !
    call MEMCOF(C_tot, iteration, m_cof, a0, cof, w1, w2, wm)

    dt_mem    =  1.d0/(m_cof*2.d0)  
    fdt(0)    = -dt_mem

    do i = 1, m_cof
       fdt(i) = fdt(i-1) + dt_mem
       J(i)   = EVLMEM(fdt(i),cof,m_cof,a0)
    end do

    norm = sum(J(2:m_cof))

    do i = 2, m_cof
       J(i) = J(i)/norm
    end do

    open(unit, file=memfile, status='unknown')
    write(unit,'(a21,4x,a21)') '#   WaveNumber [cm-1]','Intensity [arb. unit]'
    do i = 2, m_cof
       write(unit,'(f20.6,4x,e20.10)') fdt(i)*freq_conv, J(i)
    end do
    close(unit)

    deallocate(C_tot)

    deallocate(J)
    deallocate(fdt)
    deallocate(omega)
    deallocate(Q)
    deallocate(w1)
    deallocate(w2)
    deallocate(wm)
    deallocate(cof)

    !---------------------------------------------------------------------

    if (.not. loverlap) then
       m     = iteration/(n_seg*4)
       mod_m = mod(iteration, n_seg*4)
    else
       m = iteration/(2*n_seg + 1)
       mod_m = mod(iteration, 2*n_seg+1)
    end if

    allocate(J(m))
    allocate(fdt(0:m))
    allocate(omega(m))
    allocate(Q(m))
    allocate(p(m))
    allocate(w1(4*m))
    allocate(w2(m))

    dt_spc  = 1.d0/(m*2.d0)  
    fdt(0)  = -dt_spc

    do i = 1, m
       fdt(i) = fdt(i-1) + dt_spc
    end do

    omega(1:m) = fdt(1:m)*(1.0d0/(time_conv*eps12))

    p       =  0.0d0
    w1      =  0.0d0
    w2      =  0.0d0
    J       =  0.0d0
    norm    =  0.0d0

    call SPCTRM (p,m,n_seg,loverlap,w1,w2,unit_vac,vacfile,win,win_name)  
    call QCF    (omega,qcfactor,temp,m,qcf_name,Q)

    norm = sum(p(2:m)*Q(2:m))

    do i = 2, m
       J(i) = p(i)*Q(i)/norm
    end do

    open(unit, file=spcfile, status='unknown')
    write(unit,'(a21,4x,a21)') '#   WaveNumber [cm-1]','Intensity [arb. unit]'
    do i = 2, m
       write(unit,'(f20.6,4x,e20.10)') fdt(i)*freq_conv, J(i)
    end do
    close(unit)


    if (lspe) then

       do i = 1, nsp
          vacfilespe(i)  = trim(filevac)//'_raw.'//trim(syst_ID_lab(i))//'.dat'
          spcfilespe(i)  = trim(file)    //'.'//trim(syst_ID_lab(i))//'.dat'
       end do

       do ii = 1, nsp

          p       =  0.0d0
          w1      =  0.0d0
          w2      =  0.0d0
          J       =  0.0d0
          norm    =  0.0d0

          call SPCTRM(p,m,n_seg,loverlap,w1,w2,unit_vac,vacfilespe(ii),win,win_name)

          norm = sum(p(2:m)*Q(2:m))/syst_ID_spe(ii,2)*sum(syst_ID_spe(:,2))

          do i = 2, m
             J(i) = p(i)*Q(i)/norm
          end do

          open(3000+ii, file=spcfilespe(ii), status='unknown')
          write(3000+ii,'(a21,4x,a21)') '#   WaveNumber [cm-1]','Intensity [arb. unit]'
          do i = 2, m
             write(3000+ii,'(f20.6,4x,e20.10)') fdt(i)*freq_conv, J(i)
          end do
          close(3000+ii)

       end do

    end if


    deallocate(J)
    deallocate(fdt)
    deallocate(omega)
    deallocate(Q)
    deallocate(p,w1,w2)

    write(*,*) '==== PSD parameters ======================================>>'  
    write(*,*)
    write(*,'(a29)')                     'Power Spectrum Estimation'
    write(*,'(a31)')                   '>>=========================<<'
    write(*,*)
    write(*,'(a32)')                      '1.Maximum Entropy Method (MEM)'
    write(*,'(a22,2x,i10)')               'MEM order', m_cof
    write(*,'(a22,2x,f10.6,a6,f12.6,a5)') 'fdt',       dt_mem/time_conv,' ps =>', dt_mem*freq_conv, ' cm-1'
    write(*,'(a24)')                      'output file(s)'
    write(*,'(20x,a20)')                   memfile
    write(*,*)
    write(*,'(a18)')                      '2.Data Windowing'
    write(*,'(a22,2x,a11)')               'window function', win_name 
    write(*,'(a22,2x,f10.6,a6,f12.6,a5)') 'fdt',        dt_spc/time_conv,' ps =>', dt_spc*freq_conv, ' cm-1'
    write(*,'(a22,2x,l)')                 'overlap',    loverlap
    write(*,'(a22,2x,l)')                 'speciation', lspe
    write(*,'(a22,2x,i10)')               'n_seg',      n_seg
    write(*,'(a22,2x,i10)')               'm',          m
    write(*,'(a22,2x,i10)')               'mod_m',      mod_m 

    write(*,'(a24)')                      'output file(s)'
    write(*,'(20x,a20)')                   spcfile
    if (nsp > 1 .and. lspe) then
       do i = 1, nsp
          write(*,'(20x,a20)')                spcfilespe(i)
       end do
    end if
    write(*,*)
    write(*,'(a29)')                      '3.Quantum Correction Factor'
    write(*,'(a22,2x,f12.1,2a)')          'Temperature', temp,' K'
    write(*,'(a22,2x,a18)')               'QC-factor', qcf_name
    write(*,*) '==========================================================<<'


    return

  end subroutine PSD



  subroutine QCF(omega,k,T,n,qcf_name,q)

    !==-----------------------------------------------------------------------==
    !    QCF = q(omega) = "Quantum Correction Factor"
    !==-----------------------------------------------------------------------==
    !   
    !   see for example J.Chem.Phys.(2004) Vol. 12, p. 6412  and ref. therein
    !   for the definition of the various QCF or the work of A. Kohlmeyer
    !
    !==-----------------------------------------------------------------------==

    implicit none

    integer,           intent(in)  :: k, n
    real(DP),          intent(in)  :: omega(n), T
    real(DP)                       :: beta, hbar 

    character(18),     intent(out) :: qcf_name
    real(DP),          intent(out) :: q(n)

    integer                        :: i

    hbar = hp/(2.0d0*pi)
    beta = 1.0d0/(kb*T)

    q = 0.0d0

    select case(k)
    case(1)
       do i = 2, n
          q(i) = 2.0d0/(1.0d0 + exp(-beta*hbar*omega(i)))
       end do
       qcf_name = 'standard          '
    case(2)
       do i = 2, n
          write(*,*) omega(i)
          q(i) = beta*hbar*omega(i)/(1.0d0 - exp(-beta*hbar*omega(i)))

          write(*,*) q(i)
       end do
       qcf_name = 'harmonic          '
    case(3)
       do i = 2, n
          q(i) = exp(beta*hbar*omega(i)/2.0d0)
       end do
       qcf_name = 'Schofield         '
    case(4)
       do i = 2, n
          q(i) = (((beta*hbar*omega(i))/(1.0d0 - exp(-beta*hbar*omega(i))))**0.5d0)*exp(beta*hbar*omega(i)/4.0d0)
       end do
       qcf_name = 'harmonic-Schofield'
    case DEFAULT
       q = 1.0d0
       qcf_name = 'none              '

    end select

  end subroutine QCF



  subroutine VACF(nat, iteration, vel, nsp, syst_ID_spe, syst_ID_lab, C_nsp, C_tot, file, unit, &
       time_step, interval_i_m, iprint)

    !==-----------------------------------------------------------------------==
    !    VACF = C(t) = "velocity auto-correlation function" 
    !==-----------------------------------------------------------------------==
    !
    !  C(t) = 1/N*SUM^nt_i SUM^n_j SUM^3_k v_jk(t_i)v_jk(t_i + t)
    !
    !  where
    !        N             : normalization factor
    !                      = SUM^nt_i SUM^n_j SUM^3_k v_jk(t_i)v_jk(t_i)
    !        nt            : number of time origins
    !        n             : number of particles                   
    !        v_jk          : velocity component (k) of the particle (j)
    !       
    !  Remark:
    !        the code also give the C(t) of each atomic species
    !
    !==-----------------------------------------------------------------------==

    implicit none

    integer,            intent(in)             :: nat, iteration, nsp, unit
    real(DP),           intent(in)             :: vel(3,nat*iteration)
    integer,            intent(in)             :: syst_ID_spe(nsp,2)
    character(2),       intent(in)             :: syst_ID_lab(nsp)

    real(DP),           intent(in)             :: time_step
    integer,            intent(in)             :: interval_i_m, iprint

    real(DP),           intent(out)            :: C_nsp(nsp,iteration), C_tot(iteration)   
    character(256),     intent(inout)          :: file   


    real(DP)                         :: N_nsp(nsp), C_tot2(iteration)
    real(DP)                         :: N_tot
    real(DP)                         :: time_conv
    character(256)                   :: vacfile1, vacfile2, vacfile3, vacout1(nsp), vacout2(nsp)

    integer                          :: i, j, k, l, m, t, tt, t0, time, na
    integer                          :: unit_vac1, unit_vac2, unit_vac3

    vacfile1 = trim(file)//'_raw.dat'
    !vacfile2 = TRIM(file)//'2.dat'
    vacfile3 = trim(file)//'.dat'

    unit_vac1 = unit
    !unit_vac2 = unit + 1
    unit_vac3 = unit + 2

    N_tot = 0.0d0  
    C_tot = 0.0d0
    C_tot2= 0.0d0

    time_conv = time_step*iprint*interval_i_m*ps

    do t0 = 1, iteration*nat, nat        
       do k = 0, nat-1        
          N_tot = N_tot + vel(1,t0+k)*vel(1,t0+k) + vel(2,t0+k)*vel(2,t0+k) + vel(3,t0+k)*vel(3,t0+k)        
       end do
    end do

    k = 0
    t = 1
    time_loop: do time = 1, iteration-1     
       origin_loop: do t0 = 1, (iteration-k)*nat, nat                
          t = t0 + k*nat                  
          atom_loop: do na = 0, nat-1           
             C_tot(time) = C_tot(time) + vel(1,t0+na)*vel(1,t+na) + vel(2,t0+na)*vel(2,t+na) + vel(3,t0+na)*vel(3,t+na)           
          end do atom_loop
       end do origin_loop
       k           = k + 1
       C_tot(time) = C_tot(time)/N_tot     
    end do time_loop

    !-----------------------------------------!


    if (nsp > 1) then

       N_nsp = 0.d0
       k     = 0

       do t0 = 1, iteration*nat, nat            
          k = 0    
          do m = 1, nsp       
             do l = 1, syst_ID_spe(m,2)         
                N_nsp(m) = N_nsp(m) + vel(1,t0+k)*vel(1,t0+k) + vel(2,t0+k)*vel(2,t0+k) + vel(3,t0+k)*vel(3,t0+k)                   
                k = k + 1           
             end do
          end do
       end do

       C_nsp = 0.0d0

       k = 0
       t = 1

       time_loop2: do time = 1, iteration-1     
          origin_loop2: do t0 = 1, (iteration-k)*nat, nat                
             t  = t0 + k*nat           
             na = 0
             nsp_loop2: do m = 1, nsp                 
                atom_loop2:  do l = 1, syst_ID_spe(m,2)                  
                   C_nsp(m,time) = C_nsp(m,time) + vel(1,t0+na)*vel(1,t+na) + vel(2,t0+na)*vel(2,t+na) + vel(3,t0+na)*vel(3,t+na)              
                   na = na + 1             
                end do atom_loop2
             end do nsp_loop2
          end do origin_loop2

          k             = k + 1
          C_nsp(:,time) = C_nsp(:,time)/N_nsp(:) ! *syst_ID_spe(:,2)/nat

       end do time_loop2

       do i = 1, nsp
          C_tot2(1:iteration) = C_tot2(1:iteration) + C_nsp(i,1:iteration)*N_nsp(i) ! /syst_ID_spe(i,2)*nat
       end do

       C_tot2 = C_tot2/sum(N_nsp(1:nsp))

    end if

    do i = 1, nsp
       vacout1(i) = trim(file)//'_raw.'//trim(syst_ID_lab(i))//'.dat'
       vacout2(i) = trim(file)//'.'//trim(syst_ID_lab(i))//'.dat'
    end do


    open(unit_vac1, file=vacfile1, status='unknown')
    do i = 1, iteration
       write(unit_vac1,'(f20.10)') C_tot(i) 
    end do
    close(unit_vac1)


    !OPEN(unit_vac2, file=vacfile2, status='unknown')
    !WRITE(unit_vac2,'(1a,7x,36a)') '#','time [ps]           c(t) [arb. unit]'
    !DO i = 1, iteration     
    !   WRITE(unit_vac2,'(f20.10,f20.10)') (i-1)*time_conv, C_tot(i)  
    !END DO
    !CLOSE(unit_vac2)


    open(unit_vac3, file=vacfile3, status='unknown')
    write(unit_vac3,'(1a,7x,37a)') '#','time [ps] normalized c(t) [arb. unit]'
    do i = 1, iteration     
       write(unit_vac3,'(f20.10,f20.10)') (i-1)*time_conv, C_tot2(i)
    end do
    close(unit_vac3)

    if (nsp > 1) then
       do m = 1, nsp

          open((3000+m), file=vacout1(m), status='unknown')
          do i = 1, iteration
             write(3000+m,'(f20.10)') C_nsp(m,i)     
          end do
          close(3000+m)

          open((2000+m), file=vacout2(m), status='unknown')
          write(2000+m,'(1a,7x,37a)') '#','time [ps] normalized c(t) [arb. unit]'
          do i = 1, iteration
             write(2000+m,'(f20.10,f20.10)') (i-1)*time_conv, C_nsp(m,i)     
          end do
          close(2000+m)

       end do
    end if

    write(*,*) '==== VACF parameters =====================================>>'  
    write(*,'(a24)')            '(total) output file(s)'
    write(*,'(10x,a9, x,a20,x,a17) ') 'raw data:', vacfile1,  '(step index in x)'
    !WRITE(*,'(10x,a9, x,a20,x,a14) ') 'full set:', vacfile2,  '(time ps in x)'
    write(*,'(2x, a17,x,a20,x,a14) ') 'sum over species:', vacfile3, '(time ps in x)'
    write(*,*)
    if (nsp > 1) then
       write(*,'(a24)')         'species output file(s)'
       do i = 1, nsp
          write(*,'(20x,a20)') vacout1(i)
          write(*,'(20x,a20)') vacout2(i)
       end do
    end if
    write(*,'(a27)') 't in ps / c(t) normalized'
    write(*,*) '==========================================================<<'


    return
  end subroutine VACF




end module dynpro_psd
