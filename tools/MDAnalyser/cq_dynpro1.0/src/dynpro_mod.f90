MODULE dynpro_mod

  USE kinds, ONLY      : DP
  USE constants,  ONLY : pi, eps12
  USE constants,  ONLY : ps   => AU_PS
  USE constants,  ONLY : c    => C_SI
  USE constants,  ONLY : kb   => K_BOLTZMANN_SI
  USE constants,  ONLY : hp   => H_PLANCK_SI
  USE constants,  ONLY : bohr => BOHR_RADIUS_ANGS
  USE constants,  ONLY : ry   => AUTOEV

  USE dynpro_cell, ONLY : center_cell
  USE dynpro_inp,  ONLY : ldebug
  
  
IMPLICIT NONE

CONTAINS


!!$SUBROUTINE PSD(type, nsp, lspe, iteration, syst_ID_lab, syst_ID_spe, &
!!$     filevac, loverlap, n_seg, win, unit_vac, m_cof, time_step, & 
!!$     interval_i_m, iprint, unit, file, qcfactor, temp)
!!$
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER,        INTENT(in)        :: iprint, unit, nsp, type
!!$  REAL(DP),       INTENT(in)        :: time_step
!!$  INTEGER,        INTENT(in)        :: iteration, interval_i_m
!!$  INTEGER,        INTENT(in)        :: syst_ID_spe(nsp,2)
!!$  CHARACTER(2),   INTENT(in)        :: syst_ID_lab(nsp)
!!$
!!$  INTEGER,        INTENT(in)        :: n_seg, unit_vac ! SPC
!!$  CHARACTER(256), INTENT(in)        :: file, filevac
!!$  LOGICAL,        INTENT(in)        :: loverlap, lspe
!!$
!!$  INTEGER,        INTENT(inout)     :: m_cof           ! MEM
!!$
!!$  INTEGER,        INTENT(in)        :: qcfactor        ! QCF
!!$  REAL(DP),       INTENT(in)        :: temp
!!$
!!$  REAL(DP),       ALLOCATABLE       :: C_tot(:)
!!$
!!$  CHARACTER(18)                     :: qcf_name
!!$  CHARACTER(256)                    :: memfile, spcfile, vacfile
!!$  CHARACTER(256)                    :: vacfilespe(nsp), spcfilespe(nsp)
!!$
!!$  REAL(DP),        ALLOCATABLE      :: J(:),  omega(:), Q(:) 
!!$  REAL(DP),        ALLOCATABLE      :: w1(:), w2(:), wm(:), cof(:), p(:), fdt(:)
!!$  REAL(DP)                          :: a0, dt_mem         
!!$  REAL(DP)                          :: dt_spc
!!$
!!$  CHARACTER(11)                     :: win_name
!!$
!!$  REAL(DP)                          :: time_conv, freq_conv, ps_to_cm, norm ! conversion
!!$
!!$  INTEGER                           :: i, ii, m, win, mod_m, ios
!!$ 
!!$  ! 1 cm^-1 = 2.99792458E10  s^-1
!!$  
!!$  ps_to_cm  = c*eps12*100.0d0
!!$  time_conv = time_step*iprint*interval_i_m*ps ! ps
!!$  freq_conv = (1.0d0/time_conv)/ps_to_cm       ! 1/dt 
!!$
!!$  print*, '   ps_to_cm :', ps_to_cm
!!$  print*, '   time_conv:', time_conv
!!$  print*, '   time_conv:', freq_conv
!!$
!!$  !print*, type,  nsp, lspe, iteration, syst_ID_lab, syst_ID_spe, &
!!$  !   filevac, loverlap, n_seg, win, unit_vac, m_cof, time_step,  & 
!!$  !   interval_i_m, iprint, unit, file, qcfactor, temp
!!$  
!!$  !SELECT CASE(type)
!!$  !CASE(0)
!!$  !   vacfile = TRIM(filevac)//'1.dat'
!!$  !CASE(1)   
!!$  vacfile = TRIM(filevac)//'_raw.dat'
!!$
!!$  !END SELECT
!!$  memfile = TRIM(file)//'_mem.dat'
!!$  spcfile = TRIM(file)//'_spc.dat'
!!$ 
!!$  !--------------------------------------------------------------------
!!$
!!$  IF (m_cof == 0) THEN
!!$     m_cof  = iteration/10
!!$  ELSE
!!$     IF (m_cof >= iteration) THEN
!!$        WRITE(*,*)'Error -PSD- check the m_cof and nframes values !'
!!$        WRITE(*,*)'m_cof =', m_cof, 'sp_n =', iteration
!!$        STOP
!!$     END IF
!!$  END IF
!!$
!!$  ALLOCATE(fdt(0:m_cof))
!!$  ALLOCATE(J(m_cof))
!!$  ALLOCATE(omega(m_cof))
!!$  ALLOCATE(Q(m_cof))
!!$  ALLOCATE(w1(iteration))
!!$  ALLOCATE(w2(iteration))  
!!$  ALLOCATE(wm(m_cof))
!!$  ALLOCATE(cof(m_cof))
!!$  !
!!$  ALLOCATE(C_tot(iteration))
!!$  !
!!$  open(unit_vac, action="read", iostat=ios, file=vacfile, status='old', position='rewind')
!!$  !
!!$  IF (ios /= 0) THEN
!!$     PRINT *,"Error -DynPro PSD- opening files:  ", vacfile
!!$     STOP
!!$  END IF
!!$  !
!!$  print*, ios
!!$  DO i = 1, iteration
!!$     READ(unit_vac,'(f20.10)') C_tot(i)
!!$  END DO
!!$  !
!!$  close(unit_vac)
!!$  !
!!$  CALL MEMCOF(C_tot, iteration, m_cof, a0, cof, w1, w2, wm)
!!$
!!$  dt_mem    =  1.d0/(m_cof*2.d0)  
!!$  fdt(0)    = -dt_mem
!!$  
!!$  DO i = 1, m_cof
!!$     fdt(i) = fdt(i-1) + dt_mem
!!$     J(i)   = EVLMEM(fdt(i),cof,m_cof,a0)
!!$  END DO
!!$
!!$  norm = SUM(J(2:m_cof))
!!$  
!!$  DO i = 2, m_cof
!!$     J(i) = J(i)/norm
!!$  END DO
!!$
!!$  OPEN(unit, file=memfile, status='unknown')
!!$  WRITE(unit,'(a21,4x,a21)') '#   WaveNumber [cm-1]','Intensity [arb. unit]'
!!$  DO i = 2, m_cof
!!$     WRITE(unit,'(f20.6,4x,e20.10)') fdt(i)*freq_conv, J(i)
!!$  END DO
!!$  CLOSE(unit)
!!$
!!$  DEALLOCATE(C_tot)
!!$  
!!$  DEALLOCATE(J)
!!$  DEALLOCATE(fdt)
!!$  DEALLOCATE(omega)
!!$  DEALLOCATE(Q)
!!$  DEALLOCATE(w1)
!!$  DEALLOCATE(w2)
!!$  DEALLOCATE(wm)
!!$  DEALLOCATE(cof)
!!$
!!$  !---------------------------------------------------------------------
!!$
!!$  IF (.not. loverlap) THEN
!!$     m     = iteration/(n_seg*4)
!!$     mod_m = mod(iteration, n_seg*4)
!!$  ELSE
!!$     m = iteration/(2*n_seg + 1)
!!$     mod_m = mod(iteration, 2*n_seg+1)
!!$  END IF
!!$
!!$  ALLOCATE(J(m))
!!$  ALLOCATE(fdt(0:m))
!!$  ALLOCATE(omega(m))
!!$  ALLOCATE(Q(m))
!!$  ALLOCATE(p(m))
!!$  ALLOCATE(w1(4*m))
!!$  ALLOCATE(w2(m))
!!$
!!$  dt_spc  = 1.d0/(m*2.d0)  
!!$  fdt(0)  = -dt_spc
!!$ 
!!$  DO i = 1, m
!!$     fdt(i) = fdt(i-1) + dt_spc
!!$  END DO
!!$
!!$  omega(1:m) = fdt(1:m)*(1.0d0/(time_conv*eps12))
!!$
!!$  p       =  0.0d0
!!$  w1      =  0.0d0
!!$  w2      =  0.0d0
!!$  J       =  0.0d0
!!$  norm    =  0.0d0
!!$
!!$  CALL SPCTRM (p,m,n_seg,loverlap,w1,w2,unit_vac,vacfile,win,win_name)  
!!$  CALL QCF    (omega,qcfactor,temp,m,qcf_name,Q)
!!$  
!!$  norm = SUM(p(2:m)*Q(2:m))
!!$
!!$  DO i = 2, m
!!$     J(i) = p(i)*Q(i)/norm
!!$  END DO
!!$
!!$  OPEN(unit, file=spcfile, status='unknown')
!!$  WRITE(unit,'(a21,4x,a21)') '#   WaveNumber [cm-1]','Intensity [arb. unit]'
!!$  DO i = 2, m
!!$     WRITE(unit,'(f20.6,4x,e20.10)') fdt(i)*freq_conv, J(i)
!!$  END DO
!!$  CLOSE(unit)
!!$
!!$
!!$  IF (lspe) THEN
!!$
!!$     DO i = 1, nsp
!!$        vacfilespe(i)  = TRIM(filevac)//'_raw.'//TRIM(syst_ID_lab(i))//'.dat'
!!$        spcfilespe(i)  = TRIM(file)    //'.'//TRIM(syst_ID_lab(i))//'.dat'
!!$     END DO
!!$
!!$     DO ii = 1, nsp
!!$        
!!$        p       =  0.0d0
!!$        w1      =  0.0d0
!!$        w2      =  0.0d0
!!$        J       =  0.0d0
!!$        norm    =  0.0d0
!!$
!!$        CALL SPCTRM(p,m,n_seg,loverlap,w1,w2,unit_vac,vacfilespe(ii),win,win_name)
!!$     
!!$        norm = SUM(p(2:m)*Q(2:m))/syst_ID_spe(ii,2)*SUM(syst_ID_spe(:,2))
!!$        
!!$        DO i = 2, m
!!$           J(i) = p(i)*Q(i)/norm
!!$        END DO
!!$        
!!$        OPEN(3000+ii, file=spcfilespe(ii), status='unknown')
!!$        WRITE(3000+ii,'(a21,4x,a21)') '#   WaveNumber [cm-1]','Intensity [arb. unit]'
!!$        DO i = 2, m
!!$           WRITE(3000+ii,'(f20.6,4x,e20.10)') fdt(i)*freq_conv, J(i)
!!$        END DO
!!$        CLOSE(3000+ii)
!!$        
!!$     END DO
!!$     
!!$  END IF
!!$
!!$
!!$  DEALLOCATE(J)
!!$  DEALLOCATE(fdt)
!!$  DEALLOCATE(omega)
!!$  DEALLOCATE(Q)
!!$  DEALLOCATE(p,w1,w2)
!!$
!!$  WRITE(*,*) '==== PSD parameters ======================================>>'  
!!$  WRITE(*,*)
!!$  WRITE(*,'(a29)')                     'Power Spectrum Estimation'
!!$  WRITE(*,'(a31)')                   '>>=========================<<'
!!$  WRITE(*,*)
!!$  WRITE(*,'(a32)')                      '1.Maximum Entropy Method (MEM)'
!!$  WRITE(*,'(a22,2x,i10)')               'MEM order', m_cof
!!$  WRITE(*,'(a22,2x,f10.6,a6,f12.6,a5)') 'fdt',       dt_mem/time_conv,' ps =>', dt_mem*freq_conv, ' cm-1'
!!$  WRITE(*,'(a24)')                      'output file(s)'
!!$  WRITE(*,'(20x,a20)')                   memfile
!!$  WRITE(*,*)
!!$  WRITE(*,'(a18)')                      '2.Data Windowing'
!!$  WRITE(*,'(a22,2x,a11)')               'window function', win_name 
!!$  WRITE(*,'(a22,2x,f10.6,a6,f12.6,a5)') 'fdt',        dt_spc/time_conv,' ps =>', dt_spc*freq_conv, ' cm-1'
!!$  WRITE(*,'(a22,2x,l)')                 'overlap',    loverlap
!!$  WRITE(*,'(a22,2x,l)')                 'speciation', lspe
!!$  WRITE(*,'(a22,2x,i10)')               'n_seg',      n_seg
!!$  WRITE(*,'(a22,2x,i10)')               'm',          m
!!$  WRITE(*,'(a22,2x,i10)')               'mod_m',      mod_m 
!!$  
!!$  WRITE(*,'(a24)')                      'output file(s)'
!!$  WRITE(*,'(20x,a20)')                   spcfile
!!$  IF (nsp > 1 .and. lspe) THEN
!!$     DO i = 1, nsp
!!$        WRITE(*,'(20x,a20)')                spcfilespe(i)
!!$     END DO
!!$  END IF
!!$  WRITE(*,*)
!!$  WRITE(*,'(a29)')                      '3.Quantum Correction Factor'
!!$  WRITE(*,'(a22,2x,f12.1,2a)')          'Temperature', temp,' K'
!!$  WRITE(*,'(a22,2x,a18)')               'QC-factor', qcf_name
!!$  WRITE(*,*) '==========================================================<<'
!!$
!!$
!!$RETURN
!!$
!!$END SUBROUTINE PSD
!!$
!!$
!!$
!!$SUBROUTINE QCF(omega,k,T,n,qcf_name,q)
!!$  
!!$  !==-----------------------------------------------------------------------==
!!$  !    QCF = q(omega) = "Quantum Correction Factor"
!!$  !==-----------------------------------------------------------------------==
!!$  !   
!!$  !   see for example J.Chem.Phys.(2004) Vol. 12, p. 6412  and ref. therein
!!$  !   for the definition of the various QCF or the work of A. Kohlmeyer
!!$  !
!!$  !==-----------------------------------------------------------------------==
!!$
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER,           INTENT(in)  :: k, n
!!$  REAL(DP),          INTENT(in)  :: omega(n), T
!!$  REAL(DP)                       :: beta, hbar 
!!$
!!$  CHARACTER(18),     INTENT(out) :: qcf_name
!!$  REAL(DP),          INTENT(out) :: q(n)
!!$
!!$  INTEGER                        :: i
!!$
!!$  hbar = hp/(2.0d0*pi)
!!$  beta = 1.0d0/(kb*T)
!!$
!!$  q = 0.0d0
!!$
!!$  SELECT CASE(k)
!!$  CASE(1)
!!$     DO i = 2, n
!!$        q(i) = 2.0d0/(1.0d0 + exp(-beta*hbar*omega(i)))
!!$     END DO
!!$     qcf_name = 'standard          '
!!$  CASE(2)
!!$     DO i = 2, n
!!$        write(*,*) omega(i)
!!$        q(i) = beta*hbar*omega(i)/(1.0d0 - exp(-beta*hbar*omega(i)))
!!$        
!!$        write(*,*) q(i)
!!$     END DO
!!$     qcf_name = 'harmonic          '
!!$  CASE(3)
!!$     DO i = 2, n
!!$        q(i) = exp(beta*hbar*omega(i)/2.0d0)
!!$     END DO
!!$     qcf_name = 'Schofield         '
!!$  CASE(4)
!!$     DO i = 2, n
!!$        q(i) = (((beta*hbar*omega(i))/(1.0d0 - exp(-beta*hbar*omega(i))))**0.5d0)*exp(beta*hbar*omega(i)/4.0d0)
!!$     END DO
!!$     qcf_name = 'harmonic-Schofield'
!!$  CASE DEFAULT
!!$     q = 1.0d0
!!$     qcf_name = 'none              '
!!$     
!!$  END SELECT  
!!$
!!$END SUBROUTINE QCF
!!$
!!$
!!$
!!$SUBROUTINE VACF(nat, iteration, vel, nsp, syst_ID_spe, syst_ID_lab, C_nsp, C_tot, file, unit, &
!!$     time_step, interval_i_m, iprint)
!!$  
!!$  !==-----------------------------------------------------------------------==
!!$  !    VACF = C(t) = "velocity auto-correlation function" 
!!$  !==-----------------------------------------------------------------------==
!!$  !
!!$  !  C(t) = 1/N*SUM^nt_i SUM^n_j SUM^3_k v_jk(t_i)v_jk(t_i + t)
!!$  !
!!$  !  where
!!$  !        N             : normalization factor
!!$  !                      = SUM^nt_i SUM^n_j SUM^3_k v_jk(t_i)v_jk(t_i)
!!$  !        nt            : number of time origins
!!$  !        n             : number of particles                   
!!$  !        v_jk          : velocity component (k) of the particle (j)
!!$  !       
!!$  !  Remark:
!!$  !        the code also give the C(t) of each atomic species
!!$  !
!!$  !==-----------------------------------------------------------------------==
!!$
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER,            INTENT(in)             :: nat, iteration, nsp, unit
!!$  REAL(DP),           INTENT(in)             :: vel(3,nat*iteration)
!!$  INTEGER,            INTENT(in)             :: syst_ID_spe(nsp,2)
!!$  CHARACTER(2),       INTENT(in)             :: syst_ID_lab(nsp)
!!$
!!$  REAL(DP),           INTENT(in)             :: time_step
!!$  INTEGER,            INTENT(in)             :: interval_i_m, iprint
!!$  
!!$  REAL(DP),           INTENT(out)            :: C_nsp(nsp,iteration), C_tot(iteration)   
!!$  CHARACTER(256),     INTENT(inout)          :: file   
!!$
!!$
!!$  REAL(DP)                         :: N_nsp(nsp), C_tot2(iteration)
!!$  REAL(DP)                         :: N_tot
!!$  REAL(DP)                         :: time_conv
!!$  CHARACTER(256)                   :: vacfile1, vacfile2, vacfile3, vacout1(nsp), vacout2(nsp)
!!$
!!$  INTEGER                          :: i, j, k, l, m, t, tt, t0, time, na
!!$  INTEGER                          :: unit_vac1, unit_vac2, unit_vac3
!!$  
!!$  vacfile1 = TRIM(file)//'_raw.dat'
!!$  !vacfile2 = TRIM(file)//'2.dat'
!!$  vacfile3 = TRIM(file)//'.dat'
!!$
!!$  unit_vac1 = unit
!!$  !unit_vac2 = unit + 1
!!$  unit_vac3 = unit + 2
!!$
!!$  N_tot = 0.0d0  
!!$  C_tot = 0.0d0
!!$  C_tot2= 0.0d0
!!$
!!$  time_conv = time_step*iprint*interval_i_m*ps
!!$
!!$  DO t0 = 1, iteration*nat, nat        
!!$     DO k = 0, nat-1        
!!$        N_tot = N_tot + vel(1,t0+k)*vel(1,t0+k) + vel(2,t0+k)*vel(2,t0+k) + vel(3,t0+k)*vel(3,t0+k)        
!!$     END DO
!!$  END DO
!!$    
!!$  k = 0
!!$  t = 1
!!$  time_loop: DO time = 1, iteration-1     
!!$     origin_loop: DO t0 = 1, (iteration-k)*nat, nat                
!!$        t = t0 + k*nat                  
!!$        atom_loop: DO na = 0, nat-1           
!!$           C_tot(time) = C_tot(time) + vel(1,t0+na)*vel(1,t+na) + vel(2,t0+na)*vel(2,t+na) + vel(3,t0+na)*vel(3,t+na)           
!!$        END DO atom_loop        
!!$     END DO origin_loop     
!!$     k           = k + 1
!!$     C_tot(time) = C_tot(time)/N_tot     
!!$  END DO time_loop
!!$  
!!$  !-----------------------------------------!
!!$  
!!$  
!!$  IF (nsp > 1) THEN
!!$     
!!$     N_nsp = 0.d0
!!$     k     = 0
!!$     
!!$     DO t0 = 1, iteration*nat, nat            
!!$        k = 0    
!!$        DO m = 1, nsp       
!!$           DO l = 1, syst_ID_spe(m,2)         
!!$              N_nsp(m) = N_nsp(m) + vel(1,t0+k)*vel(1,t0+k) + vel(2,t0+k)*vel(2,t0+k) + vel(3,t0+k)*vel(3,t0+k)                   
!!$              k = k + 1           
!!$           END DO
!!$        END DO
!!$     END DO
!!$     
!!$     C_nsp = 0.0d0
!!$     
!!$     k = 0
!!$     t = 1
!!$     
!!$     time_loop2: DO time = 1, iteration-1     
!!$        origin_loop2: DO t0 = 1, (iteration-k)*nat, nat                
!!$           t  = t0 + k*nat           
!!$           na = 0
!!$           nsp_loop2: DO m = 1, nsp                 
!!$              atom_loop2:  DO l = 1, syst_ID_spe(m,2)                  
!!$                 C_nsp(m,time) = C_nsp(m,time) + vel(1,t0+na)*vel(1,t+na) + vel(2,t0+na)*vel(2,t+na) + vel(3,t0+na)*vel(3,t+na)              
!!$                 na = na + 1             
!!$              END DO atom_loop2
!!$           END DO nsp_loop2
!!$        END DO origin_loop2
!!$        
!!$        k             = k + 1
!!$        C_nsp(:,time) = C_nsp(:,time)/N_nsp(:) ! *syst_ID_spe(:,2)/nat
!!$        
!!$     END DO time_loop2
!!$     
!!$     DO i = 1, nsp
!!$        C_tot2(1:iteration) = C_tot2(1:iteration) + C_nsp(i,1:iteration)*N_nsp(i) ! /syst_ID_spe(i,2)*nat
!!$     END DO
!!$     
!!$     C_tot2 = C_tot2/SUM(N_nsp(1:nsp))
!!$     
!!$  END IF
!!$
!!$  DO i = 1, nsp
!!$     vacout1(i) = TRIM(file)//'_raw.'//TRIM(syst_ID_lab(i))//'.dat'
!!$     vacout2(i) = TRIM(file)//'.'//TRIM(syst_ID_lab(i))//'.dat'
!!$  END DO
!!$  
!!$
!!$  OPEN(unit_vac1, file=vacfile1, status='unknown')
!!$  DO i = 1, iteration
!!$     WRITE(unit_vac1,'(f20.10)') C_tot(i) 
!!$  END DO
!!$  CLOSE(unit_vac1)
!!$  
!!$
!!$  !OPEN(unit_vac2, file=vacfile2, status='unknown')
!!$  !WRITE(unit_vac2,'(1a,7x,36a)') '#','time [ps]           c(t) [arb. unit]'
!!$  !DO i = 1, iteration     
!!$  !   WRITE(unit_vac2,'(f20.10,f20.10)') (i-1)*time_conv, C_tot(i)  
!!$  !END DO
!!$  !CLOSE(unit_vac2)
!!$
!!$
!!$  OPEN(unit_vac3, file=vacfile3, status='unknown')
!!$  WRITE(unit_vac3,'(1a,7x,37a)') '#','time [ps] normalized c(t) [arb. unit]'
!!$  DO i = 1, iteration     
!!$     WRITE(unit_vac3,'(f20.10,f20.10)') (i-1)*time_conv, C_tot2(i)
!!$  END DO
!!$  CLOSE(unit_vac3)
!!$  
!!$  IF (nsp > 1) THEN
!!$     DO m = 1, nsp
!!$ 
!!$        OPEN((3000+m), file=vacout1(m), status='unknown')
!!$        DO i = 1, iteration
!!$           WRITE(3000+m,'(f20.10)') C_nsp(m,i)     
!!$        END DO
!!$        CLOSE(3000+m)
!!$        
!!$        OPEN((2000+m), file=vacout2(m), status='unknown')
!!$        WRITE(2000+m,'(1a,7x,37a)') '#','time [ps] normalized c(t) [arb. unit]'
!!$        DO i = 1, iteration
!!$           WRITE(2000+m,'(f20.10,f20.10)') (i-1)*time_conv, C_nsp(m,i)     
!!$        END DO
!!$        CLOSE(2000+m)
!!$        
!!$     END DO
!!$  END IF
!!$
!!$  WRITE(*,*) '==== VACF parameters =====================================>>'  
!!$  WRITE(*,'(a24)')            '(total) output file(s)'
!!$  WRITE(*,'(10x,a9, x,a20,x,a17) ') 'raw data:', vacfile1,  '(step index in x)'
!!$  !WRITE(*,'(10x,a9, x,a20,x,a14) ') 'full set:', vacfile2,  '(time ps in x)'
!!$  WRITE(*,'(2x, a17,x,a20,x,a14) ') 'sum over species:', vacfile3, '(time ps in x)'
!!$  WRITE(*,*)
!!$  IF (nsp > 1) THEN
!!$     WRITE(*,'(a24)')         'species output file(s)'
!!$     DO i = 1, nsp
!!$        WRITE(*,'(20x,a20)') vacout1(i)
!!$        WRITE(*,'(20x,a20)') vacout2(i)
!!$     END DO
!!$  END IF
!!$  WRITE(*,'(a27)') 't in ps / c(t) normalized'
!!$  WRITE(*,*) '==========================================================<<'
!!$
!!$
!!$  RETURN
!!$END SUBROUTINE VACF


!!$FUNCTION EVLMEM(FDT,COF,M,PM)
!!$  
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER                      :: M
!!$  REAL(DP)                     :: PM,FDT,COF(M)
!!$  REAL(DP)                     :: WR,WI,WPR,WPI,WTEMP,THETA
!!$  REAL(DP)                     :: SUMI,SUMR,EVLMEM
!!$  INTEGER                      :: I
!!$  
!!$  THETA=6.28318530717959D0*FDT
!!$  WPR=DCOS(THETA)
!!$  WPI=DSIN(THETA)
!!$  WR=1.D0
!!$  WI=0.D0
!!$  SUMR=1.d0
!!$  SUMI=0.d0
!!$  DO 11 I=1,M
!!$     WTEMP=WR
!!$     WR=WR*WPR-WI*WPI
!!$     WI=WI*WPR+WTEMP*WPI
!!$     SUMR=SUMR-COF(I)*SNGL(WR)
!!$     SUMI=SUMI-COF(I)*SNGL(WI)
!!$11 END DO
!!$  EVLMEM=PM/(SUMR**2+SUMI**2)
!!$  RETURN
!!$
!!$END FUNCTION EVLMEM


!!$SUBROUTINE MEMCOF(DATA,N,M,PM,COF,WK1,WK2,WKM)
!!$  
!!$  !==-----------------------------------------------------------------------==
!!$  !    MEMCOF = "Maximum Entropy Method" 
!!$  !==-----------------------------------------------------------------------==
!!$  !
!!$  !  From: Numerical Recipies in Fortran 77, vol. 1, p. 561
!!$  !
!!$  ! " ...given a real vector data(1...n), ad given m, this routine returns 
!!$  !      vector cof(1...m) with cof(j) = a_j, and a scalar pm=a0, which are
!!$  !      the coefficient for Maximum Entropy Method spectral estimation..."
!!$  ! 
!!$  !  Remark:
!!$  !      enforced variable declaration have been considered here
!!$  !
!!$  !==-----------------------------------------------------------------------==
!!$
!!$
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER,           INTENT(in)     :: N,M
!!$  REAL(DP),          INTENT(out)    :: PM,COF(M)
!!$  REAL(DP),          INTENT(inout)  :: DATA(N),WK1(N),WK2(N),WKM(M)    
!!$ 
!!$  INTEGER                           :: I,J,K
!!$  REAL(DP)                          :: P,PNEUM,DENOM
!!$    
!!$  P=0.d0 
!!$  DO 11 J=1,N
!!$     P=P+DATA(J)**2
!!$11 END DO
!!$  PM=P/N
!!$  WK1(1)=DATA(1)
!!$  WK2(N-1)=DATA(N)
!!$  DO 12 J=2,N-1
!!$     WK1(J)=DATA(J)
!!$     WK2(J-1)=DATA(J)
!!$12 END DO
!!$  DO 17 K=1,M
!!$     PNEUM=0.d0
!!$     DENOM=0.d0
!!$     DO 13 J=1,N-K
!!$        PNEUM=PNEUM+WK1(J)*WK2(J)
!!$        DENOM=DENOM+WK1(J)**2+WK2(J)**2
!!$13   END DO
!!$     COF(K)=2.d0*PNEUM/DENOM
!!$     PM=PM*(1.d0-COF(K)**2)
!!$     IF(K.NE.1)THEN
!!$        DO 14 I=1,K-1
!!$           COF(I)=WKM(I)-COF(K)*WKM(K-I)
!!$14      END DO
!!$     ENDIF
!!$     IF(K.EQ.M)RETURN
!!$     DO 15 I=1,K
!!$        WKM(I)=COF(I)
!!$15   END DO
!!$     DO 16 J=1,N-K-1
!!$        WK1(J)=WK1(J)-WKM(K)*WK2(J)
!!$        WK2(J)=WK2(J+1)-WKM(K)*WK1(J+1)
!!$16   END DO
!!$17 END DO
!!$
!!$  !PAUSE 'never get here'
!!$  
!!$END SUBROUTINE MEMCOF


SUBROUTINE sampling(matrix_in, n, m, o, p, start, end, delta, iprint, &
     matrix_out)

  IMPLICIT NONE

  INTEGER,  INTENT(in)         :: n, m, o, p, start, delta, end, iprint
  REAL(DP), INTENT(in)         :: matrix_in(n,m*o)
  REAL(DP), INTENT(inout)      :: matrix_out(n,m*p)
  INTEGER                      :: i, j, k

  k = 1
  DO i = start-1, end-1, delta
     DO j = 1, m
        matrix_out(:,k) = matrix_in(:,i*m+j)        
        k = k + 1        
     END DO
 END DO

END SUBROUTINE sampling


SUBROUTINE tricky(nsp, nat, atm, ityp, label, syst_ID_spe, syst_ID_lab)

  IMPLICIT NONE

  INTEGER,          INTENT(in)     :: nsp, nat
  CHARACTER(len=3), INTENT(in)     :: atm(nsp)
  INTEGER,          INTENT(inout)  :: ityp(nat)
  CHARACTER(len=3), INTENT(inout)  :: label(nat)
  INTEGER,          INTENT(inout)  :: syst_ID_spe(nsp,2)
  CHARACTER(len=2), INTENT(inout)  :: syst_ID_lab(nsp)


  INTEGER                         :: i, j, k
  INTEGER                         :: atomic_number(nsp), na(nsp)
  CHARACTER(len=3)                :: tmp_atm(nsp)
  CHARACTER(len=2)                :: pte(103)

  !==---- Periodic Table of the Elements  --------------------------------

  DATA pte /  "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", &
              "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", &
              "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
              "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr", &
              "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "St", &
              "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
              "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
              "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
              "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
              "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", &
              "Md", "No", "Lr"/

  
  tmp_atm = ADJUSTL(atm) 


  !print*, '################ TRICKY ################'
  !print*, tmp_atm
  !print*,
  !print*, ityp
  
  DO i = 1, nsp
     k = 1
     DO WHILE (TRIM(tmp_atm(i)) /= TRIM(pte(k)))
       k = k + 1
     END DO
     atomic_number(i) = k
  END DO

  if (ldebug) print*, 'atomic_number'
  if (ldebug) print*, atomic_number
  
  na = 0
  DO i = 1, nat
     na(ityp(i)) = na(ityp(i)) + 1
  END DO

  if (ldebug) print*, 'na'
  if (ldebug) print*, na
  
  k = 0
  DO i = 1, nsp
     DO j = 1, na(i)
        k = k + 1
        ityp(k)  = atomic_number(i)
        label(k) = pte(atomic_number(i))
     END DO
  END DO  

  if (ldebug) print*, 'ityp'
  if (ldebug) print*, ityp
  if (ldebug) print*, 'label'
  if (ldebug) print*, label
  
  DO i = 1, nsp
     DO j = 1, 2
        syst_ID_spe(i,1) = atomic_number(i)
        syst_ID_spe(i,2) = na(i)
        syst_ID_lab(i)   = pte(atomic_number(i))
     END DO
  END DO

  if (ldebug)  print*,  syst_ID_spe
  if (ldebug)  print*,  syst_ID_lab
  
  !print*, '################ ###### ################'

  
END SUBROUTINE tricky


SUBROUTINE tricky_sort(nsp, nat, iteration, xyz, vel, ityp_xml, ityp, syst_ID_spe, &
     label_xml, syst_ID_lab, species)

  IMPLICIT NONE

  INTEGER,    INTENT(in)    :: nsp, nat, iteration
  INTEGER,    INTENT(in)    :: ityp(nat)
  INTEGER,    INTENT(in)    :: syst_ID_spe(nsp,2)
  LOGICAL,    INTENT(in)    :: species
  CHARACTER(len=3), INTENT(in)  :: label_xml(nat)
  CHARACTER(len=2), INTENT(in)  :: syst_ID_lab(nsp)
  
  INTEGER,    INTENT(inout) :: ityp_xml(nat)  
  REAL(DP),   INTENT(inout) :: xyz(3,nat*iteration)
  REAL(DP),   INTENT(inout) :: vel(3,nat*iteration)

  INTEGER                   :: i, j, k, n, m, pos, count, a, b, c, d
  INTEGER                   :: index_int(nat)
  REAL(DP)                  :: index_real(nat)
  
  REAL(DP)                  :: xyz_tmp(3,nat)
  REAL(DP)                  :: vel_tmp(3,nat)
  REAL(DP)                  :: test_tmp(3,nat) 
  
  if (ldebug) write(*,*)'################ START TRICKY SORT DEBUG ################'
  !
  if (ldebug) write(*,*) 'tricky sort setup'
  if (ldebug) write(*,*)
  if (ldebug) write(*,*) 'ityp_xml'  
  if (ldebug) write(*,*) ityp_xml
  if (ldebug) write(*,*)
  if (ldebug) write(*,*) 'label_xml'  
  if (ldebug) write(*,*) label_xml
  if (ldebug) write(*,*)
  if (ldebug) write(*,*) 'syst_ID_lab'  
  if (ldebug) write(*,*) syst_ID_lab
  if (ldebug) write(*,*) 
  if (ldebug) write(*,*) 'syst_ID_spe'  
  if (ldebug) write(*,*) syst_ID_spe
  !
  if (ldebug) write(*,*) 'ityp_xml befor loop'
  if (ldebug) write(*,*) ityp_xml
  !
  if ( species ) then
     !
     !do i = 1, nat
     !   do j = 1, nsp       
     !      if ( ityp_xml(i) == j ) then
     !         ityp_xml(i) = syst_ID_spe(j,1)
     !      end if
     !   end do
     !end do
     !
     do i = 1, nat
        do j = 1, nsp       
           !
           !print*, label_xml(i),  syst_ID_lab(j)
           if ( label_xml(i) == syst_ID_lab(j) ) then
              ityp_xml(i) = syst_ID_spe(j,1)
              !
           end if
           !
        end do
     end do
     !     
  end if
  !
  if (ldebug) print*, 'ityp_xml after loop'
  if (ldebug) print*, ityp_xml
  !
  index_int = 0
  pos       = 1
  !
  do j = 1, nsp        
     do i = 1, nat
        !
        if ( ityp_xml(i) == syst_ID_spe(j,1) ) then
           index_int(i) = pos            
           pos = pos + 1
           !
        end if
        !
     end do
  end do
  !
  index_real = real(index_int)
  !
  if (ldebug) print*, 'index_int'
  if (ldebug) print*, index_int
  !
  do i = 1, iteration

     a = nat * (i - 1) + 1
     b = nat *  i

     xyz_tmp = xyz(:, a:b)
     vel_tmp = vel(:, a:b)

     test_tmp   = xyz_tmp
     index_real = real(index_int)
  
     call sort_vect(index_real, nat, xyz_tmp, vel_tmp)

     !print*, 'iteration', i
     !print*, index_real
     !do j = 1, nat
     !   print*, xyz_tmp(1,j), test_tmp(1,j)
        !print*, index_real
       
     !end do

     xyz(:, a:b) = xyz_tmp
     vel(:, a:b) = vel_tmp
          
  end do

  if (ldebug) print*, '################ END TRICKY SORT DEBUG ################'
  
END SUBROUTINE tricky_sort



SUBROUTINE very_tricky(nsp_n, nsp, nat, atm, ityp, label_old, label, syst_ID_spe, &
     syst_ID_lab, sub_xyz, k)

  IMPLICIT NONE

  INTEGER,          INTENT(in)     :: nsp, nat, nsp_n
  CHARACTER(len=3), INTENT(inout)  :: atm(nsp)
  INTEGER,          INTENT(inout)  :: ityp(nat)
  CHARACTER(len=3), INTENT(inout)  :: label(nat)
  CHARACTER(len=3), INTENT(inout)  :: label_old(nat)
  INTEGER,          INTENT(inout)  :: syst_ID_spe(nsp,2)
  CHARACTER(len=2), INTENT(inout)  :: syst_ID_lab(nsp)
 

  INTEGER                          :: sub_row(100), n_sub, sub_xyz1(100), sub_xyz2(100)
  INTEGER                          :: sub_xyz(100)
  INTEGER                          :: ityp_perm(nat), ityp_1, ityp_2
  CHARACTER(len=3)                 :: sub_lab(100), tmp_lab(100)


  CHARACTER(len=3)                 :: label_perm(nat)
  CHARACTER(len=2)                 :: pte(103)
  
  INTEGER                          :: i, j, k, k1, k2, l, tmp


  !==---- Periodic Table of the Elements  --------------------------------

  DATA pte /  "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", &
              "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", &
              "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
              "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr", &
              "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "St", &
              "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
              "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
              "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
              "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
              "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", &
              "Md", "No", "Lr"/

  
  k          = 1
  l          = 1
  DO i = 1, nat
     IF (label(i) /= label_old(i)) THEN
        sub_row(k) = i
        sub_lab(k) = label(i)
        k = k + 1
     ELSE        
        label_perm(l) = label(i)
        ityp_perm(l) = ityp(i)
        l = l + 1
     END IF
  END DO
  

  n_sub   = k - 1
  tmp_lab = sub_lab

  sub_xyz1(1)      =  sub_row(1)        
  atm(nsp-nsp_n+1) =  sub_lab(1)

  k1 = 2
  k2 = 1
  DO i = 2, n_sub
     IF (tmp_lab(1) == sub_lab(i)) THEN 
        sub_xyz1(k1) = sub_row(i)
        k1 = k1 + 1
     ELSE
        IF (nsp_n > 1) THEN
           atm(nsp)   =  sub_lab(i)
           sub_xyz2(k2) = sub_row(i)           
           k2 = k2 + 1
        END IF
     END IF
     
  END DO

  DO i = 1, 103
     IF(ADJUSTL(atm(nsp-nsp_n+1)) == pte(i)) THEN
        ityp_1 = i
     END IF
     IF(ADJUSTL(atm(nsp)) == pte(i)) THEN
        ityp_2 = i
     END IF
  END DO

  k1 = k1 - 1
  k2 = k2 - 1
  l  = l 

  if (ldebug) print*, 'label_perm'
  if (ldebug) print*, label_perm
  
  DO i = 1, k1
     label_perm(l-1+i) = atm(nsp-nsp_n+1)
     ityp_perm(l-1+i)  = ityp_1
  END DO
  
  IF (nsp_n > 1) THEN
     DO i = 1, k2
        label_perm(l-1+k1+i) = atm(nsp)
        ityp_perm(l-1+k1+i)  = ityp_2
     END DO
  END IF

  tmp = 1
  DO i = 1, nat
     ityp(i) = tmp
     IF(ityp_perm(i+1) /= ityp_perm(i)) THEN
        tmp = tmp + 1
     END IF
  END DO

  CALL tricky(nsp, nat, atm, ityp, label_perm, syst_ID_spe, syst_ID_lab)

  k = k1 + k2
  
  sub_xyz(1:k1) = sub_xyz1(1:k1)

  IF (nsp_n > 1) THEN
     sub_xyz(k1+1:k) = sub_xyz2(:)
  END IF

  label = label_perm

END SUBROUTINE very_tricky



SUBROUTINE permute(n, m, matrix, perm, p)
  
  IMPLICIT NONE

  INTEGER,  INTENT(in)         :: n, m, p
  INTEGER,  INTENT(in)         :: perm(p)

  REAL(DP), INTENT(inout)      :: matrix(3,n*m)
  REAL(DP)                     :: matrix_tmp(3,n*m)

  REAL(DP)                     :: matrix_tmp1(3,(n-p)*m)
  REAL(DP)                     :: matrix_tmp2(3,(p)*m)

  INTEGER                      :: i, j, k, l, o, q
  

  q = 0

  DO j = 1, m
 
     k = 1
     l = 1
     o = 1

     DO i = 1, n

        !print*, 'here j,i', j, i
        
        IF (i /= perm(o)) THEN
           matrix_tmp1(:,l) = matrix(:,i+q)
           l = l + 1
        ELSE
           matrix_tmp2(:,k) = matrix(:,i+q)
           k = k + 1
           o = o + 1
        END IF
        
        matrix_tmp(:,1+q:q+n-p)   =  matrix_tmp1(:,:)
        matrix_tmp(:,n-p+1+q:q+n) =  matrix_tmp2(:,:)

        !print*
        
     END DO

     q = q + n
         
  END DO
  
  matrix = matrix_tmp


END SUBROUTINE permute



SUBROUTINE check_sampling(sp_start, sp_end, evp, iteration, iprint, &
     sp_n, n, m, interval_i_m, time_step, time_dt)
  
  IMPLICIT NONE

  INTEGER,  INTENT(in)         :: iteration, iprint
  REAL(DP), INTENT(in)         :: evp(8,iteration), time_step
  INTEGER,  INTENT(inout)      :: sp_start, sp_end, sp_n
  REAL(DP), INTENT(out)        :: time_dt
  INTEGER,  INTENT(out)        :: n, m, interval_i_m
  INTEGER                      :: nconf, interval_m, interval_i, sp_end_m
  REAL(DP)                     :: interval_r, time, nyquist, ps_to_cm

  ! 1 cm^-1 = 2.99792458E10  s^-1
  !         = 2.99792458E-2 ps^-1
  ps_to_cm  = 0.0299792458d0

  interval_r = 0.0d0
  interval_i = 0 

  IF (sp_start == 0) THEN
     sp_start = INT(evp(1,1))          ! start sampling 
  END IF
  IF (sp_end == 0) THEN
     sp_end   = INT(evp(1,iteration))  ! end sampling
  END IF
 
  print*, 'SAMPLING:', iteration
 
  n = 1
  DO WHILE (sp_start /= INT(evp(1,n)))
     n = n + 1
     IF (n == iteration) THEN
        WRITE(*,*)'Error -sampling- check the sp_start/sp_end values !'
        WRITE(*,*)'sp_start =',n
        STOP
     END IF
  END DO
  m = n
  DO WHILE (sp_end   /= INT(evp(1,m)))
     m = m + 1
     IF (m > iteration) THEN
        WRITE(*,*)'Error -sampling- check the sp_start/sp_end values !'
        WRITE(*,*)'sp_end =',m
        STOP
     END IF
  END DO

  nconf = (sp_end - sp_start)/iprint + 1

  IF (sp_n == 0) THEN
     sp_n = nconf
  END IF

  IF (sp_n > nconf) THEN
     WRITE(*,*)'Error -sampling- the number of frame (sp_n) is too high !'
     WRITE(*,*)'Maximum number of frame =',nconf
     STOP
  ELSE
     IF (sp_n == 1) THEN
        WRITE(*,*)'Error -sampling- sp_n must be > 1 !'
        WRITE(*,*)'sp_n =',sp_n
        STOP          
     ELSE
        IF (sp_n <= nconf) THEN
           interval_r = DBLE(sp_end - sp_start)/(DBLE(sp_n)-1)
           interval_i = INT(interval_r)
        END IF
     END IF
  END IF

  interval_m    = MOD(interval_i,iprint)
  interval_i_m  = (interval_i - interval_m)/iprint
  sp_end_m      = sp_start + interval_i_m*(sp_n-1)*iprint
  m             = n + interval_i_m*(sp_n-1)
  time          = (nconf-1)*time_step*ps*iprint 
  time_dt       = time_step*iprint*interval_i_m
  nyquist       = 1.0d0/(2.0d0*time_dt*ps)

  WRITE(*,*) '====== Sampling parameters ===============================>>'
  WRITE(*,'(a22,2x,f10.1,a6,f12.6,a3)')  'CP time step', time_step,' au => ', time_step*ps, ' ps'
  WRITE(*,'(a22,2x,i10)')                'sp_n (nframes)', sp_n
  WRITE(*,'(a22,2x,i10)')                'iprint', iprint
  WRITE(*,*)
  WRITE(*,'(a22,2x,2i10)')               'sp_start', sp_start, n
  WRITE(*,'(a22,2x,i10)')                'sp_end', sp_end
  WRITE(*,'(a22,2x,i10)')                'nconf max allowed', nconf
  WRITE(*,'(a22,2x,f10.6,a3)')           'time allowed', time, ' ps'
  WRITE(*,'(a22,2x,i10)')                'mod_n', interval_m
  WRITE(*,'(a22,2x,2i10)')               'interval_n', interval_i, interval_i_m
  WRITE(*,'(a22,2x,2i10)')               'sp_end_n', sp_end_m, m
  WRITE(*,*)
  WRITE(*,'(a22,2x,f10.1,a7,f12.6,a3)') 'time sampling interval',      time_dt, ' au => ', time_dt*ps,          ' ps'
  WRITE(*,'(a22,2x,f10.1,a7,f12.6,a5)') 'Nyquist freq. ',               nyquist, ' THz=> ', nyquist*ps_to_cm,    ' cm-1'
  WRITE(*,'(a22,2x,f10.1,a7,f12.6,a3)') 'total sampled time',  time_dt*(sp_n-1), ' au => ', time_dt*(sp_n-1)*ps, ' ps'
  WRITE(*,*) '==========================================================<<'
           
  RETURN
  
END SUBROUTINE check_sampling


SUBROUTINE title(version)

  IMPLICIT NONE
 
  
  CHARACTER(len=5),   INTENT(in)  :: version

  WRITE(*,*) '==========================================================>>'
  WRITE(*,*) '| DYNAMICS | STRUCTURAL PROPERTIES >>'
  WRITE(*,*) '==========================================================<<'


  RETURN
END SUBROUTINE title


SUBROUTINE sort(nsp, natoms, nframes, ct_atom, syst_ID_spe, &
     cel, abc, vel, sort_xyz, sort_vel)

  IMPLICIT NONE

  INTEGER,  INTENT(in)        :: nsp, ct_atom, natoms, nframes
  INTEGER,  INTENT(in)        :: syst_ID_spe(nsp,2)
  
  REAL(DP), INTENT(in)        :: cel(3,3*nframes)
  REAL(DP), INTENT(in)        :: abc(3,natoms*nframes)
  REAL(DP), INTENT(in)        :: vel(3,natoms*nframes)
  
  REAL(DP), INTENT(out)       :: sort_xyz(3,natoms*nframes)
  REAL(DP), INTENT(out)       :: sort_vel(3,natoms*nframes)
  

  INTEGER                     :: i, j, k, l, m, n, o, p, q, cc, dd

  REAL(DP)                    :: atom_ID_abc(3,natoms)
  REAL(DP)                    :: velo_ID_tmp(3,natoms)
  REAL(DP)                    :: atom_ID_tmp(3)
  REAL(DP)                    :: cell_ID_tmp(3,3)
  
  
  REAL(DP)                    :: r(natoms)
  REAL(DP)                    :: ct_abc(3,natoms), ct_xyz(3,natoms)
  

  p = 0
  o = 0
  r = 0.0d0
  frames_loop1: DO q = 1, nframes
     
     atom_ID_abc(:,:) = abc(:,(q+p):(q*natoms))

     velo_ID_tmp(:,:) = vel(:,(q+p):(q*natoms))

     cell_ID_tmp(:,:) = cel(:,(q+o):(q*3))
     
     
     CALL center_cell(natoms, 1, ct_atom, cell_ID_tmp, &
          atom_ID_abc, ct_abc, ct_xyz)
     
     atom_ID_tmp(1:3) = ct_xyz(1:3,ct_atom)            
     
     cc = 0
     dd = 0
     m  = 1

     species_loop1: DO j = 1, nsp     
        
        ext_atom_loop1: DO k = 1, syst_ID_spe(j,2)           
           
           r(m) = sqrt((atom_ID_tmp(1) - ct_xyz(1,m))**2 &                
                +        (atom_ID_tmp(2) - ct_xyz(2,m))**2 &
                +        (atom_ID_tmp(3) - ct_xyz(3,m))**2)
           
           m    = m + 1
           
        END DO ext_atom_loop1
        
        dd = dd + syst_ID_spe(j,2)
        
        IF (syst_ID_spe(j,1) == 8 .or. syst_ID_spe(j,1) == 1) THEN

           CALL sort_vect(r(1+cc:dd), syst_ID_spe(j,2), ct_xyz(1:3,1+cc:dd), velo_ID_tmp(1:3,1+cc:dd))

        END IF

        cc = cc + syst_ID_spe(j,2)

     END DO species_loop1
   
     sort_xyz(:,(q+p):(q*natoms)) =  ct_xyz(:,:)
     sort_vel(:,(q+p):(q*natoms)) =  velo_ID_tmp(:,:)
     
     o = o + 2
     p = p + (natoms - 1)
    
  END DO frames_loop1
  
  RETURN
  
END SUBROUTINE sort

SUBROUTINE class(nsp, natoms, nframes, sph_atom, syst_ID_spe, syst_ID_lab, &
     sort_xyz, sort_vel, sph_atn)
  
  IMPLICIT NONE

  INTEGER,  INTENT(in)           :: nsp, natoms, nframes
  INTEGER,  INTENT(inout)        :: sph_atom
  INTEGER,  INTENT(inout)        :: syst_ID_spe(nsp,2)
  REAL(DP), INTENT(inout)        :: sort_xyz(3,natoms*nframes)  
  CHARACTER(len=2),INTENT(inout) :: syst_ID_lab(nsp) 
  INTEGER,  INTENT(in)           :: sph_atn(100)

  REAL(DP), INTENT(inout)        :: sort_vel(3,natoms*nframes)
  

  REAL(DP)                       :: class_xyz(3,natoms*nframes)

  REAL(DP)                       :: class_vel(3,natoms*nframes)
   

  REAL(DP)                       :: sph_atom_xyz(3)
  INTEGER                        :: class_syst_ID_spe(nsp,2)
  CHARACTER(len=2)               :: class_syst_ID_lab(nsp) 

  INTEGER                        :: i, j, k, l, m, n, o, p
    
  sph_atom_xyz(:) = sort_xyz(:,sph_atom)

  m = 0
  DO i = 1, nsp
     p = 0
     l = 1
     DO WHILE (sph_atn(i) /= syst_ID_spe(l,1)) 
        l = l + 1
     END DO

     class_syst_ID_lab(i)   = syst_ID_lab(l)
     class_syst_ID_spe(i,1) = syst_ID_spe(l,1)
     class_syst_ID_spe(i,2) = syst_ID_spe(l,2)

     IF (l == 1) THEN
        n = 0
     ELSE
        n = SUM(syst_ID_spe(1:l-1,2))
     END IF
     DO j = 1, nframes
        DO k = 1, syst_ID_spe(l,2)

           class_xyz(:,k+p+m) = sort_xyz(:,k+p+n)         
           class_vel(:,k+p+m) = sort_vel(:,k+p+n)

        END DO
        p = p + natoms
     END DO
     m = m + syst_ID_spe(l,2)
  END DO

  syst_ID_spe = class_syst_ID_spe
  syst_ID_lab = class_syst_ID_lab
  sort_xyz    = class_xyz
  sort_vel    = class_vel

  k        = 1
  sph_atom = 1
  DO WHILE (sph_atom_xyz(1) /= sort_xyz(1,k) .and. sph_atom_xyz(2) /= sort_xyz(1,k) &
       .and. sph_atom_xyz(3) /= sort_xyz(1,k))
     sph_atom = sph_atom + 1
     k = k + 1
  END DO 


END SUBROUTINE class


SUBROUTINE sort_vect(v, n, xyz, vel)

  IMPLICIT NONE

  INTEGER,  INTENT(in)        :: n
  REAL(DP), INTENT(inout)     :: v(n), xyz(3,n), vel(3,n)
  
  INTEGER                     :: i, j, rgmax
  REAL(DP)                    :: xyz_m(3), vel_m(3) 
  REAL(DP)                    :: max

  ! max <= min !

    DO i = 1, n
     max      = v(i)
     rgmax    = i
     xyz_m(:) = xyz(:,i)
     vel_m(:) = vel(:,i)
     DO j = i + 1, n
        IF (max > v(j)) THEN
           max      = v(j)
           xyz_m(:) = xyz(:,j)
           vel_m(:) = vel(:,j)
           rgmax    = j
        END IF        
     END DO
     v(rgmax)     = v(i)
     v(i)         = max
     xyz(:,rgmax) = xyz(:,i)
     vel(:,rgmax) = vel(:,i)
     xyz(:,i)     = xyz_m(:)
     vel(:,i)     = vel_m(:)
  END DO


  RETURN
END SUBROUTINE sort_vect

SUBROUTINE moment(DATA,N,AVE,ADEV,SDEV,VAR)
  
  IMPLICIT NONE
  
  INTEGER,  INTENT(in)     :: N
  REAL(DP), INTENT(in)     :: DATA(N)
  
  REAL(DP), INTENT(out)    :: AVE, ADEV, SDEV, VAR
  REAL(DP)                 :: SKEW, CURT
  
  INTEGER                  :: I, J
  REAL(DP)                 :: S, P
  
  IF(N.LE.1) THEN
     WRITE(*,*)  'STOP: N must be at least 2'
     STOP
  END IF

  S = 0.0D0
  P = 0.0D0     
  DO 11 J = 1, N
     S = S + DATA(J)
11 END DO
  
  AVE  = S/N
  ADEV = 0.0D0
  VAR  = 0.0D0
  SKEW = 0.0D0
  CURT = 0.0D0
  
  DO 12 J = 1, N
     S    = DATA(J) - AVE
     ADEV = ADEV + ABS(S)
     P   = S*S
     VAR = VAR + P
     P   = P*S
     SKEW= SKEW + P
     P   = P*S
     CURT= CURT + P
12 END DO
  
  ADEV = ADEV/N
  VAR  = VAR/(N - 1)
  SDEV = SQRT(VAR)
  IF(VAR.NE.0.) THEN
     SKEW = SKEW/(N*SDEV**3)
     CURT = CURT/(N*VAR**2) - 3.0D0
  ELSE
     WRITE(*,*) 'no skew or kurtosis when zero variance'
     SKEW = 0.0d0
     CURT = 0.0d0
  ENDIF
  RETURN
END SUBROUTINE moment



!!$SUBROUTINE dNab (nsp, natoms,  nframes, nr, rmin, rmax,   & 
!!$     ityp, syst_ID_spe, cel, abc,  &
!!$     rd_atom, rd_count, dN, event, V)
!!$!==-----------------------------------------------------------------------==
!!$!    dN_{ab}(r) = "number of particle pair probability" 
!!$!==-----------------------------------------------------------------------==
!!$!
!!$! given 2 definite species a and b separated by a distance ranging
!!$! in [r, r+dr], dN_{ab} is number of pair a-b in this range
!!$!
!!$!==-----------------------------------------------------------------------==
!!$
!!$  IMPLICIT NONE
!!$ 
!!$  INTEGER,  INTENT(in)        :: nsp, natoms, rd_atom, rd_count, nr
!!$  INTEGER,  INTENT(in)        :: ityp(natoms), syst_ID_spe(nsp,2)
!!$
!!$  INTEGER,  INTENT(in)        :: nframes
!!$  REAL(DP), INTENT(in)        :: rmin, rmax
!!$ 
!!$  REAL(DP), INTENT(in)        :: cel(3,3*nframes)
!!$  REAL(DP), INTENT(in)        :: abc(3,natoms*nframes)
!!$
!!$  INTEGER,  INTENT(out)       :: dN(nframes,rd_count,nsp,nr)
!!$  INTEGER,  INTENT(out)       :: event(nframes,rd_count,nsp)
!!$  REAL(DP)                    :: vol(nframes), V
!!$
!!$  INTEGER                     :: i, j, k, l, m, n, o, p, q
!!$  INTEGER                     :: atom_ID_nbr(rd_count)
!!$  REAL(DP)                    :: rad, r, rstep
!!$  REAL(DP)                    :: atom_ID_abc(3,natoms)
!!$  REAL(DP)                    :: atom_ID_tmp(3,rd_count)
!!$  REAL(DP)                    :: cell_ID_tmp(3,3)
!!$  REAL(DP)                    :: a(3), b(3), c(3), ab(3)
!!$  REAL(DP)                    :: ct_abc(3,natoms), ct_xyz(3,natoms)
!!$
!!$  o     = 0
!!$  p     = 0
!!$  event = 0
!!$  dN    = 0
!!$  r     = 0.0d0
!!$  a     = 0.0d0    
!!$  b     = 0.0d0
!!$  c     = 0.0d0
!!$  ab    = 0.0d0
!!$  rstep = (rmax - rmin) / DBLE(nr - 1)
!!$ 
!!$  j = 1
!!$  DO i = 1, natoms     
!!$     IF (ityp(i) == REAL(rd_atom)) THEN
!!$        atom_ID_nbr(j) = i 
!!$        j = j + 1
!!$     END IF
!!$  END DO
!!$ 
!!$  frames_loop: DO q = 1, nframes
!!$     
!!$     atom_ID_abc(:,:) = abc(:,(q+p):(q*natoms))
!!$     cell_ID_tmp(:,:) = cel(:,(q+o):(q*3))
!!$     a(:)             = cel(1,(q+o):(q*3))
!!$     b(:)             = cel(2,(q+o):(q*3))
!!$     c(:)             = cel(3,(q+o):(q*3))
!!$     call vec_product(a, b, ab)
!!$     vol(q)           = DOT_PRODUCT(ab,c)
!!$     
!!$     rd_atom_loop: DO i = 1, rd_count
!!$        m  = 1      
!!$        
!!$        CALL center_cell(natoms, 1, atom_ID_nbr(i), cell_ID_tmp, &
!!$             atom_ID_abc, ct_abc, ct_xyz)
!!$                
!!$        j = 1
!!$        DO k = 1, natoms     
!!$           IF (ityp(k) == REAL(rd_atom)) THEN
!!$              atom_ID_tmp(1:3,j) = ct_xyz(1:3,k)       
!!$              j = j + 1
!!$           END IF
!!$        END DO
!!$        
!!$        species_loop: DO j = 1, nsp     
!!$           
!!$           ext_atom_loop: DO k = 1, syst_ID_spe(j,2)           
!!$               
!!$              r = sqrt((atom_ID_tmp(1,i) - ct_xyz(1,m))**2 &                
!!$                   +   (atom_ID_tmp(2,i) - ct_xyz(2,m))**2 &
!!$                   +   (atom_ID_tmp(3,i) - ct_xyz(3,m))**2)
!!$              
!!$              m   = m + 1
!!$              rad = rmin
!!$              
!!$              rad_point_loop: DO l = 1, nr   
!!$                 
!!$                 IF (r <= rad .and. r /= 0.0d0) THEN
!!$                    dN(q,i,j,l)  = dN(q,i,j,l)  + 1      
!!$                    event(q,i,j) = event(q,i,j) + 1
!!$                    EXIT 
!!$                 END IF
!!$                 rad = rad + rstep
!!$                 
!!$              END DO rad_point_loop
!!$              
!!$           END DO ext_atom_loop
!!$           
!!$        END DO species_loop
!!$        
!!$     END DO rd_atom_loop
!!$     
!!$     o = o + 2
!!$     p = p + (natoms - 1)
!!$
!!$  END DO frames_loop
!!$
!!$  V = SUM(vol)/nframes  
!!$
!!$  RETURN
!!$END SUBROUTINE DNAB


!!$SUBROUTINE RDF (dN, event, count, rd_atom, rd_count, rd_lab, &
!!$     natoms, nsp, nr, rmin, rmax, V, & 
!!$     syst_ID_spe, syst_ID_lab, filerdf, lau_unit)
!!$!==-----------------------------------------------------------------------==
!!$!    RDF = g_{ab}(r) = "pair correlation function" (Radial Distribution)
!!$!==-----------------------------------------------------------------------==
!!$!
!!$!  g_{ab}(r) = (V/N_{ab))<rho_{ab}(r)>
!!$!
!!$!  where
!!$!        <rho_{ab}(r)> : the local density of (pair) particle
!!$!                      = SUM_{a,b,a/=b}(1/(4*pi*r**2*dr))*dN_{ab}
!!$!                           
!!$!        V             : volume of the cell
!!$!        N_{ab}        : number of pair of particles
!!$!                      = Na*Nb      if a /= b
!!$!                      = Na(Na - 1) if a == b 
!!$!  Remark:
!!$!        (N_{ab}/V)*int[dr g_{ab}(r)*4*pi*r**2] 
!!$!                      = Na*Nb       if a /= b
!!$!                      = Na*(Na - 1) if a == b
!!$!==-----------------------------------------------------------------------==
!!$  
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER,  INTENT(in)        :: natoms, nsp, nr, count, rd_atom, rd_count
!!$  INTEGER,  INTENT(in)        :: dN(count,rd_count,nsp,nr)
!!$  INTEGER,  INTENT(in)        :: event(count,rd_count,nsp)
!!$  INTEGER,  INTENT(in)        :: syst_ID_spe(nsp,2)
!!$  REAL(DP), INTENT(in)        :: rmin, rmax
!!$  REAL(DP), INTENT(in)        :: V
!!$  LOGICAL,    INTENT(in)  :: lau_unit
!!$  CHARACTER(256), INTENT(in)  :: filerdf
!!$  CHARACTER(2),   INTENT(in)  :: syst_ID_lab(nsp), rd_lab
!!$
!!$  CHARACTER(256)              :: rdfout(nsp), int_g_out(nsp)
!!$  CHARACTER(4)                ::  r_unit
!!$  REAL(DP)                    :: g_int(nsp), int_g(nsp)  
!!$  REAL(DP)                    :: g_tmp(nsp,nr)
!!$  REAL(DP)                    :: rad, rstep, Nab, scal
!!$  REAL(DP)                    :: g_sum(nsp)  
!!$ 
!!$  INTEGER                     :: event_tmp(count,rd_count,nsp)
!!$  INTEGER                     :: i, j, k, l, m, n
!!$
!!$
!!$  IF(lau_unit) THEN
!!$     r_unit = 'Bohr'
!!$  ELSE
!!$     r_unit = 'Ang.'
!!$  END IF
!!$  
!!$  DO i = 1, nsp
!!$     rdfout(i) = TRIM(filerdf)//'.'//TRIM(rd_lab)//'_'//TRIM(syst_ID_lab(i))//'.dat'
!!$     int_g_out(i) = TRIM(rd_lab)//'_'//TRIM(syst_ID_lab(i))
!!$     print*, ' syst_ID_lab', syst_ID_lab(i), syst_ID_spe(i,1),  syst_ID_spe(i,2) 
!!$     print*, rdfout(i),  int_g_out(i)
!!$  END DO
!!$
!!$  rstep = (rmax - rmin) / DBLE(nr - 1)
!!$  
!!$  DO m = 1, nsp
!!$     OPEN((1000+m), file=rdfout(m), status='unknown')
!!$     WRITE(1000+m,'(a5,a4,a1,10x,a20,4x,a21)') '# r [',r_unit,']','g_ab(r) [Normalized]','g_int(r) [Normalized]'
!!$  END DO
!!$
!!$  g_int = 0.0d0
!!$  g_tmp = 0.0d0
!!$  Nab   = 0.0d0
!!$
!!$  scal  = 1.0d0
!!$
!!$  species_loop: DO i = 1, nsp 
!!$     rad = rmin
!!$     Nab = rd_count*syst_ID_spe(i,2)
!!$     rad_point_loop: DO j = 1, nr
!!$        
!!$        g_tmp(i,j) = SUM(DBLE(dN(:,:,i,j)))/(count*4*pi*rad**2*rstep*Nab/V)
!!$       
!!$        g_int(i)   = g_int(i) + (4*pi*rad**2)*rstep*Nab/V*g_tmp(i,j)
!!$        
!!$        WRITE((1000+i),'(f10.6,1X,2(f24.8))') rad, g_tmp(i,j), g_int(i)/syst_ID_spe(i,2)
!!$        
!!$        rad = rad + rstep 
!!$        
!!$        IF(j == nr) THEN
!!$           int_g(i) = g_int(i)/syst_ID_spe(i,2)
!!$        END IF
!!$
!!$     END DO rad_point_loop
!!$     
!!$     g_sum(i) = SUM(g_tmp(i,:)) 
!!$
!!$  END DO species_loop
!!$  
!!$  WRITE(*,*) '====== RDF parameters ====================================>>'
!!$  WRITE(*,'(a22,2x,a10)')   'rd_atom', rd_lab
!!$  WRITE(*,'(a22,2x,i10)')   'rd_count', rd_count
!!$  WRITE(*,'(a22,2x,i10)')   'nb. of species', nsp
!!$  WRITE(*,'(a22,2x,i10)')   'nb. of atoms', natoms
!!$  WRITE(*,'(a22,2x,i10)')
!!$  WRITE(*,'(a22,2x,i10)')   'nb. radial point', nr
!!$  WRITE(*,'(a22,2x,f10.4)') 'r_min (Ang.)', rmin
!!$  WRITE(*,'(a22,2x,f10.4)') 'r_max (Ang.)', rmax
!!$  WRITE(*,'(a22,2x,f10.4)') 'r_step (Ang.)', rstep
!!$  WRITE(*,'(a22,2x,i10)')   'nframes', count
!!$  WRITE(*,'(a22,2x,i10)')
!!$  WRITE(*,'(a22)')         'output file(s)'
!!$  DO i = 1, nsp
!!$     WRITE(*,'(20x,a20)') rdfout(i)
!!$  END DO
!!$  WRITE(*,'(a22)')         'int[g_{ab}(r)dr]'
!!$  DO i = 1, nsp
!!$     WRITE(*,'(20x,a5,2x,f8.2)') int_g_out(i), int_g(i)
!!$  END DO
!!$  WRITE(*,*) '==========================================================<<'
!!$
!!$  DO m = 1, nsp
!!$     CLOSE(1000+m)
!!$  END DO
!!$  
!!$  RETURN
!!$  
!!$END SUBROUTINE RDF


!!$SUBROUTINE vec_product(a1, a2, a3)
!!$  
!!$  IMPLICIT NONE
!!$  
!!$  REAL(DP), INTENT(in)         :: a1(3), a2(3)
!!$  REAL(DP), INTENT(out)        :: a3(3)
!!$
!!$  a3(1) = a1(2)*a2(3) - a1(3)*a2(2)
!!$  a3(2) = a1(3)*a2(1) - a1(1)*a2(3)
!!$  a3(3) = a1(1)*a2(2) - a1(2)*a2(1)
!!$
!!$  RETURN
!!$
!!$END SUBROUTINE VEC_PRODUCT

!!$SUBROUTINE xyz_to_abc(n, m, xyz, abc, cel)
!!$  
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,  INTENT(in)         :: n, m
!!$  REAL(DP), INTENT(in)         :: xyz(3,n*m)
!!$  REAL(DP), INTENT(in)         :: cel(3,3*m)
!!$  REAL(DP), INTENT(inout)      :: abc(3,n*m)
!!$
!!$  REAL(DP)                     :: at(3,3), at_inv(3,3), sigma(3,n)
!!$  INTEGER                      :: i, j, k, l, q
!!$
!!$  k = 0
!!$  l = 0
!!$  DO i = 1, m
!!$     at(:,:) = cel(1:3,(i+k):(i*3))
!!$     CALL inverse(at, at_inv)
!!$     abc(:,(i+l):(i*n)) = MATMUL(at_inv(:,:), xyz(:,(i+l):(i*n)))
!!$     k = k + 2
!!$     l = l + (n - 1)
!!$  END DO
!!$
!!$END SUBROUTINE XYZ_TO_ABC
!!$
!!$SUBROUTINE center_cell(n, m, ct_atom, cel, abc, ct_abc, ct_xyz)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,  INTENT(in)  :: n, m, ct_atom
!!$  REAL(DP), INTENT(in)  :: abc(3,n*m), cel(3,3*m)
!!$  REAL(DP), INTENT(out) :: ct_abc(3,n*m) , ct_xyz(3,n*m)
!!$
!!$  REAL(DP)      :: sigma(3, n), sigma_catom(3, n)
!!$  INTEGER       :: i, j, k, l, q
!!$
!!$  k = 0
!!$  l = 0
!!$
!!$  DO q = 1, m
!!$
!!$     sigma(:,:)  = abc(:,(q+l):(q*n))
!!$
!!$     sigma_catom = 0.0d0
!!$
!!$     DO i = 1, 3
!!$        sigma_catom(i,:) =  sigma(i,:) - sigma(i,ct_atom)
!!$        DO j = 1, n
!!$           IF (sigma_catom(i,j) < -0.5d0) THEN
!!$              sigma_catom(i,j) = sigma_catom(i,j) + 1.0d0
!!$           ELSE
!!$              IF (sigma_catom(i,j) > 0.5d0) THEN
!!$                 sigma_catom(i,j) = sigma_catom(i,j) - 1.0d0
!!$              ENDIF
!!$           END IF
!!$        END DO
!!$     END DO
!!$
!!$     ct_abc(:,(q+l):(q*n)) = sigma_catom(:,1:n)
!!$     ct_xyz(:,(q+l):(q*n)) = MATMUL(cel(:,(k+q):(q*3)), &
!!$          ct_abc(:,(q+l):(q*n)))
!!$
!!$     k = k + 2
!!$     l = l + (n - 1)
!!$
!!$  END DO
!!$
!!$  RETURN
!!$END SUBROUTINE center_cell


SUBROUTINE write_xyz(n, m, label, xyz, evp, file, unit)

  IMPLICIT NONE

  INTEGER,           INTENT(in)   :: n, m, unit
  REAL(DP),          INTENT(in)   :: xyz(3,n*m), evp(8,m)
  CHARACTER(len=3),  INTENT(in)   :: label(n)
  CHARACTER(len=256),INTENT(in)   :: file

  INTEGER                         :: i, j, k

  OPEN(unit, file=file, status='unknown')
    
  DO k = 1, m
     
     WRITE(unit,'(i10)') n
     WRITE(unit,'("frame",1X,I8,1X,"count",1X,f8.0)') k, evp(1,k)

     !print*,
     
     DO i = 1, n
        WRITE (unit,'(a2,3x,3f15.9)') TRIM(label(i)), &
             (xyz(j, n*(k - 1)+i), j = 1,3)

        !WRITE (*,'(a2,3x,3f15.9)') TRIM(label(i)), &
        !     (xyz(j, n*(k - 1)+i), j = 1,3)
        
     END DO

  END DO

  CLOSE(unit)

  RETURN
END SUBROUTINE write_xyz

SUBROUTINE write_xml(time_xml, iteration_xml, nat_xml, nsp_xml)

  IMPLICIT NONE

  
  REAL(DP), INTENT(in)          :: time_xml      
  INTEGER,  INTENT(in)          :: iteration_xml
  INTEGER,  INTENT(in)          :: nat_xml, nsp_xml

  WRITE(*,*) '==== Dynamic parameters (XML file) =======================>>'
  WRITE(*,'(a22,2x,i10)')      'iteration', iteration_xml
  WRITE(*,'(a22,2x,f10.6,a3)') 'time', time_xml, ' ps'
  WRITE(*,'(a22,2x,i10)')      'nb. of atoms', nat_xml
  WRITE(*,'(a22,2x,i10)')      'nb. of species', nsp_xml     
  WRITE(*,*) '==========================================================<<'
  

  RETURN
END SUBROUTINE write_xml

!!$SUBROUTINE inverse(at, atinv)
!!$
!!$  ! compute inverse of 3*3 matrix
!!$
!!$  IMPLICIT NONE
!!$
!!$  REAL(DP), INTENT(in)  :: at(3, 3)
!!$  REAL(DP), INTENT(out) :: atinv(3, 3)
!!$
!!$  REAL(DP) :: det
!!$
!!$  atinv(1, 1) = at(2, 2) * at(3, 3) - at(2, 3) * at(3, 2)
!!$  atinv(2, 1) = at(2, 3) * at(3, 1) - at(2, 1) * at(3, 3)
!!$  atinv(3, 1) = at(2, 1) * at(3, 2) - at(2, 2) * at(3, 1)
!!$  atinv(1, 2) = at(1, 3) * at(3, 2) - at(1, 2) * at(3, 3)
!!$  atinv(2, 2) = at(1, 1) * at(3, 3) - at(1, 3) * at(3, 1)
!!$  atinv(3, 2) = at(1, 2) * at(3, 1) - at(1, 1) * at(3, 2)
!!$  atinv(1, 3) = at(1, 2) * at(2, 3) - at(1, 3) * at(2, 2)
!!$  atinv(2, 3) = at(1, 3) * at(2, 1) - at(1, 1) * at(2, 3)
!!$  atinv(3, 3) = at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1)
!!$
!!$  det = at(1, 1) * atinv(1, 1) + at(1, 2) * atinv(2, 1) + &
!!$        at(1, 3) * atinv(3, 1)
!!$  atinv(:,:) = atinv(:,:) / det
!!$
!!$  RETURN
!!$END SUBROUTINE inverse



SUBROUTINE check_file(filename, n, m, iteration)

  IMPLICIT NONE

  CHARACTER (len=256), INTENT(IN)    :: filename
  INTEGER,             INTENT(IN)    :: n, m
  INTEGER,             INTENT(OUT)   :: iteration
  INTEGER                            :: iunit, ios, mod, rest
  INTEGER                            :: i, count

  CHARACTER(len=256)                 :: line
  
  iunit = 5
  
  OPEN (UNIT = iunit, FILE = filename, STATUS = "OLD", ACTION = "READ", &
       IOSTAT = ios, POSITION = "REWIND" )
  IF (ios /= 0) THEN
     PRINT *,"Error -DynPro read- opening file ", filename
     STOP
  END IF
  
  count = -1
  DO WHILE (ios == 0)
     READ (iunit, FMT='(A)', IOSTAT = ios) line     
     count = count + 1
  END DO
  
  rest = MODULO(count,(n + m))
  
  IF (rest /= 0) THEN
     PRINT *,"Error -DynPro read- ending file ", filename
     STOP     
  END IF

  iteration = count / (n + m)
 
END SUBROUTINE check_file


SUBROUTINE read_file(filename, unit, o, n, natoms, iteration, matrix)

  IMPLICIT NONE

  CHARACTER (len=256), INTENT(in)    :: filename
  INTEGER,             INTENT(in)    :: o, n, natoms, iteration, unit
  
  REAL(DP),            INTENT(inout) :: matrix(n,natoms*iteration)
 

  INTEGER                            :: i, j, k, l, ios
  
  OPEN(unit, file=filename, status='old', position='rewind')
 
  k   = 1
  ios = 0

  write(*,'(3x,a15,x,a80)') 'file_name read:', filename

  DO i = 1, iteration
     IF (o > 0) THEN
        DO j = 1, o
           READ(unit,*)
        END DO
     END IF
        DO j = 1, natoms
           READ(unit,*) (matrix(l,k), l = 1,n)
           !PRINT*,  (matrix(l,k), l = 1,n)
           !print*, i
           IF (ios /= 0) THEN
              PRINT *,"Error -DynPro read- file: probably not a real number ", filename
              STOP
           END IF
           k = k + 1
        END DO
  END DO
  
  CLOSE(unit)
 
END SUBROUTINE read_file





!!$SUBROUTINE SPCTRM(P,M,K,OVRLAP,W1,W2,unit_vac,filevac,win,win_name)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,           INTENT(in)     :: K,M,unit_vac
!!$  CHARACTER(len=256),INTENT(inout)  :: filevac
!!$  REAL(DP),          INTENT(inout)  :: P(M),W1(4*M),W2(m)
!!$  LOGICAL,           INTENT(in)     :: OVRLAP
!!$  
!!$  INTEGER                           :: I,J,J2,JOFF,JOFFN,KK,M4,M43,M44,MM
!!$  REAL(DP)                          :: DEN,FACM,FACP,SUMW,W
!!$  REAL(DP),          ALLOCATABLE    :: WINDOW(:)
!!$  
!!$  INTEGER,           INTENT(in)     :: win
!!$  CHARACTER(11),     INTENT(out)    :: win_name
!!$ 
!!$  MM   = M  + M
!!$  M4   = MM + MM
!!$  M44  = M4 + 4
!!$  M43  = M4 + 3
!!$  DEN  = 0.d0
!!$  FACM = M  - 0.5d0
!!$  FACP = 1.d0/(M + 0.5d0)
!!$  SUMW = 0.d0
!!$
!!$  ALLOCATE(WINDOW(MM))
!!$
!!$  SELECT CASE(win)
!!$
!!$  CASE(1)
!!$     DO J = 1, MM
!!$        WINDOW(J) = 1.d0 - ABS(((J - 1) - FACM)*FACP)
!!$     END DO
!!$     win_name = 'Bartlett   '
!!$
!!$  CASE(2)
!!$     WINDOW = 1.0d0      
!!$     win_name = 'rectangular'
!!$
!!$  CASE(3)    
!!$     DO J = 1, MM
!!$        WINDOW(J) = 1.d0 - (((J - 1) - FACM)*FACP)**2
!!$     END DO
!!$     win_name = 'Welch      '
!!$
!!$  CASE(4)   
!!$     DO J = 1, MM
!!$        WINDOW(J) = 0.5d0*(1.d0 - cos(pi*j/FACM))
!!$     END DO
!!$     win_name = 'Hann       '    
!!$
!!$  CASE DEFAULT
!!$     WINDOW = 1.0d0      
!!$     win_name = 'rectangular'
!!$      
!!$  END SELECT
!!$
!!$  open(unit_vac, file=filevac, status='unknown', position='rewind')
!!$
!!$  DO 11 J=1,MM
!!$     SUMW=SUMW+WINDOW(J)**2
!!$11 END DO
!!$  DO 12 J=1,M
!!$     P(J)=0.d0
!!$12 END DO
!!$ 
!!$  IF(OVRLAP)THEN
!!$     READ(unit_vac,'(f20.10)') (W2(J),J=1,M)
!!$  ENDIF
!!$  
!!$  DO 18 KK=1,K
!!$     DO 15 JOFF=-1,0,1
!!$        IF (OVRLAP) THEN
!!$           DO 13 J=1,M
!!$              W1(JOFF+J+J)=W2(J)
!!$13         END DO
!!$           READ(unit_vac,'(f20.10)') (W2(J),J=1,M)
!!$           JOFFN=JOFF+MM
!!$           DO 14 J=1,M
!!$              W1(JOFFN+J+J)=W2(J)
!!$14         END DO
!!$        ELSE
!!$           READ(unit_vac,'(f20.10)') (W1(J),J=JOFF+2,M4,2)
!!$        ENDIF
!!$15   END DO
!!$     DO 16 J=1,MM
!!$        J2=J+J
!!$        W=WINDOW(J)
!!$        W1(J2)=W1(J2)*W
!!$        W1(J2-1)=W1(J2-1)*W
!!$16   END DO
!!$     CALL FOUR1(W1,MM,1)
!!$     P(1)=P(1)+W1(1)**2+W1(2)**2
!!$     DO 17 J=2,M
!!$        J2=J+J
!!$        P(J)=P(J)+W1(J2)**2+W1(J2-1)**2+W1(M44-J2)**2+W1(M43-J2)**2
!!$17   END DO
!!$     DEN=DEN+SUMW
!!$18 END DO
!!$
!!$  DEN=M4*DEN
!!$
!!$  DO 19 J=1,M
!!$     P(J)=P(J)/DEN
!!$19 END DO
!!$
!!$  close(unit_vac)
!!$
!!$  DEALLOCATE(WINDOW)
!!$
!!$  RETURN
!!$END SUBROUTINE SPCTRM
!!$
!!$SUBROUTINE FOUR1(DATA,NN,ISIGN)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,   INTENT(in)    :: NN,ISIGN
!!$  REAL(DP),  INTENT(inout) :: DATA(*)
!!$
!!$  REAL(DP)                 :: WR,WI,WPR,WPI,WTEMP
!!$  REAL(DP)                 :: TEMPR,TEMPI,THETA
!!$  INTEGER                  :: N,M,MMAX,ISTEP
!!$  INTEGER                  :: I,J
!!$
!!$  N=2*NN
!!$  J=1
!!$
!!$  DO 11 I=1,N,2
!!$     IF(J.GT.I)THEN
!!$        TEMPR=DATA(J)
!!$        TEMPI=DATA(J+1)
!!$        DATA(J)=DATA(I)
!!$        DATA(J+1)=DATA(I+1)
!!$        DATA(I)=TEMPR
!!$        DATA(I+1)=TEMPI
!!$     ENDIF
!!$     M=N/2
!!$1    IF ((M.GE.2).AND.(J.GT.M)) THEN
!!$        J=J-M
!!$        M=M/2
!!$        GO TO 1
!!$     ENDIF
!!$     J=J+M
!!$11 END DO
!!$  MMAX=2
!!$2 IF (N.GT.MMAX) THEN
!!$     ISTEP=2*MMAX
!!$     THETA=6.28318530717959D0/(ISIGN*MMAX)
!!$     WPR=-2.D0*DSIN(0.5D0*THETA)**2
!!$     WPI=DSIN(THETA)
!!$     WR=1.D0
!!$     WI=0.D0
!!$     DO 13 M=1,MMAX,2
!!$        DO 12 I=M,N,ISTEP
!!$           J=I+MMAX
!!$           TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
!!$           TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
!!$           DATA(J)=DATA(I)-TEMPR
!!$           DATA(J+1)=DATA(I+1)-TEMPI
!!$           DATA(I)=DATA(I)+TEMPR
!!$           DATA(I+1)=DATA(I+1)+TEMPI
!!$12      END DO
!!$        WTEMP=WR
!!$        WR=WR*WPR-WI*WPI+WR
!!$        WI=WI*WPR+WTEMP*WPI+WI
!!$13   END DO
!!$     MMAX=ISTEP
!!$     GO TO 2
!!$  ENDIF
!!$
!!$  RETURN
!!$        
!!$END SUBROUTINE FOUR1
!!$
!!$  
!!$SUBROUTINE REALFT(DATA,N,ISIGN)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,   INTENT(in)    :: N,ISIGN
!!$  REAL(DP),  INTENT(inout) :: DATA(2*N)
!!$
!!$  REAL(DP)                 :: WR,WRS,WI,WIS,WPR,WPI,WTEMP
!!$  REAL(DP)                 :: THETA,C1,C2 
!!$  REAL(DP)                 :: H1R,H1I,H2R,H2I
!!$  INTEGER                  :: I,J,I1,I2,I3,I4,N2P3
!!$
!!$  THETA=6.28318530717959D0/2.0D0/DBLE(N)
!!$  C1=0.5
!!$
!!$  IF (ISIGN.EQ.1) THEN
!!$     C2=-0.5
!!$     CALL FOUR1(DATA,N,+1)
!!$  ELSE
!!$     C2=0.5
!!$     THETA=-THETA
!!$  ENDIF
!!$  WPR=-2.0D0*DSIN(0.5D0*THETA)**2
!!$  WPI=DSIN(THETA)
!!$  WR=1.0D0+WPR
!!$  WI=WPI
!!$  N2P3=2*N+3
!!$  DO 11 I=2,N/2+1
!!$     I1=2*I-1
!!$     I2=I1+1
!!$     I3=N2P3-I2
!!$     I4=I3+1
!!$     WRS=SNGL(WR)
!!$     WIS=SNGL(WI)
!!$     H1R=C1*(DATA(I1)+DATA(I3))
!!$     H1I=C1*(DATA(I2)-DATA(I4))
!!$     H2R=-C2*(DATA(I2)+DATA(I4))
!!$     H2I=C2*(DATA(I1)-DATA(I3))
!!$     DATA(I1)=H1R+WRS*H2R-WIS*H2I
!!$     DATA(I2)=H1I+WRS*H2I+WIS*H2R
!!$     DATA(I3)=H1R-WRS*H2R+WIS*H2I
!!$     DATA(I4)=-H1I+WRS*H2I+WIS*H2R
!!$     WTEMP=WR
!!$     WR=WR*WPR-WI*WPI+WR
!!$     WI=WI*WPR+WTEMP*WPI+WI
!!$11 END DO
!!$  IF (ISIGN.EQ.1) THEN
!!$     H1R=DATA(1)
!!$     DATA(1)=H1R+DATA(2)
!!$     DATA(2)=H1R-DATA(2)
!!$  ELSE
!!$     H1R=DATA(1)
!!$     DATA(1)=C1*(H1R+DATA(2))
!!$     DATA(2)=C1*(H1R-DATA(2))
!!$     CALL FOUR1(DATA,N,-1)
!!$  ENDIF
!!$  
!!$  RETURN
!!$END SUBROUTINE REALFT


SUBROUTINE write_evp(evp, n, time_step, iprint, lgrace, filener)

  IMPLICIT NONE

  REAL(DP),           INTENT(in)    :: time_step
  INTEGER,            INTENT(in)    :: iprint, n
  LOGICAL,            INTENT(in)    :: lgrace
  CHARACTER(len=256), INTENT(inout) :: filener
  REAL(DP),           INTENT(inout) :: evp(8,n)

  INTEGER                           :: grace
  INTEGER                           :: i, j, k

  REAL(DP)                          :: max_e, min_e
  REAL(DP)                          :: max_ks, min_ks, max_t, min_t
  REAL(DP)                          :: max_tot, min_tot, max_cp, min_cp
  REAL(DP)                          :: De, Dks, Dt, Dtot, Dcp, max_time

  grace    = 30

  evp(1,:) =  (evp(1,:)-evp(1,1))*time_step*ps


  OPEN(grace, file=filener, status='unknown')

  max_time= MAXVAL(evp(1,:))
  max_e   = MAXVAL(evp(2,:)) ; min_e   = MINVAL(evp(2,:)) ; De   = ABS(max_e   - min_e)/2
  max_t   = MAXVAL(evp(4,:)) ; min_t   = MINVAL(evp(4,:)) ; Dt   = ABS(max_t   - min_t)/2
  max_ks  = MAXVAL(evp(5,:)) ; min_ks  = MINVAL(evp(5,:)) ; Dks  = ABS(max_ks  - min_ks)/2
  max_tot = MAXVAL(evp(7,:)) ; min_tot = MINVAL(evp(7,:)) ; Dtot = ABS(max_tot - min_tot)/2
  max_cp  = MAXVAL(evp(8,:)) ; min_cp  = MINVAL(evp(8,:)) ; Dcp  = ABS(max_cp  - min_cp)/2
  
  !==---------------------------------------------------------------------==
  !   .evp file contains the following quantities :
  !   ---------------------------------------------
  !   NFI    [int]          - step index
  !   EKINC  [HARTREE A.U.] - kinetic energy of the fictitious electronic dynamics
  !   TEMPH  [K]            - Temperature of the fictitious cell dynamics
  !   TEMP   [K]            - Ionic temperature
  !   ETOT   [HARTREE A.U.] - Scf total energy (Kohn-Sham hamiltonian)
  !   ENTHAL [HARTREE A.U.] - Enthalpy ( ETOT + P * V )
  !   ECONS  [HARTREE A.U.] - Enthalpy + kinetic energy of ions and cell
  !   ECONT  [HARTREE A.U.] - Constant of motion for the CP lagrangian
  !==---------------------------------------------------------------------==
  !  1 / 2 / 4 / 5 / 7/ 8 


  IF (.not. lgrace) THEN
     
     WRITE(grace,'(6(a,4X))') '# time(ps)','e-kin(au)','T(K)','Eks(au)','Etot(au)','CP(au)' 
     
     DO k = 1, n
        WRITE(grace,'(6(f0.8,4X))') evp(1,k),evp(2,k),evp(4,k),evp(5,k),evp(7,k),evp(8,k)  
     END DO
     
  ELSE
     
     ! Xmgrace Format
write(grace,*) '# Grace project file                                            '    
write(grace,*) '#								'
write(grace,*) '@version 50122							'
write(grace,*) '@page size 792, 612						'
write(grace,*) '@page scroll 5% 						'
write(grace,*) '@page inout 5%							'
write(grace,*) '@link page off							'
write(grace,*) '@map font 0 to "Times-Roman", "Times-Roman"			'
write(grace,*) '@map font 1 to "Times-Italic", "Times-Italic"			'
write(grace,*) '@map font 2 to "Times-Bold", "Times-Bold"			'
write(grace,*) '@map font 3 to "Times-BoldItalic", "Times-BoldItalic"		'
write(grace,*) '@map font 4 to "Helvetica", "Helvetica" 			'
write(grace,*) '@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique" 	'
write(grace,*) '@map font 6 to "Helvetica-Bold", "Helvetica-Bold"		'
write(grace,*) '@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique" '
write(grace,*) '@map font 8 to "Courier", "Courier"				'
write(grace,*) '@map font 9 to "Courier-Oblique", "Courier-Oblique"		'
write(grace,*) '@map font 10 to "Courier-Bold", "Courier-Bold"			'
write(grace,*) '@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"	'
write(grace,*) '@map font 12 to "Symbol", "Symbol"				'
write(grace,*) '@map font 13 to "ZapfDingbats", "ZapfDingbats"			'
write(grace,*) '@map color 0 to (255, 255, 255), "white"			'
write(grace,*) '@map color 1 to (0, 0, 0), "black"				'
write(grace,*) '@map color 2 to (255, 0, 0), "red"				'
write(grace,*) '@map color 3 to (0, 255, 0), "green"				'
write(grace,*) '@map color 4 to (0, 0, 255), "blue"				'
write(grace,*) '@map color 5 to (255, 255, 0), "yellow" 			'
write(grace,*) '@map color 6 to (188, 143, 143), "brown"			'
write(grace,*) '@map color 7 to (220, 220, 220), "grey" 			'
write(grace,*) '@map color 8 to (148, 0, 211), "violet" 			'
write(grace,*) '@map color 9 to (0, 255, 255), "cyan"				'
write(grace,*) '@map color 10 to (255, 0, 255), "magenta"			'
write(grace,*) '@map color 11 to (255, 165, 0), "orange"			'
write(grace,*) '@map color 12 to (114, 33, 188), "indigo"			'
write(grace,*) '@map color 13 to (103, 7, 72), "maroon" 			'
write(grace,*) '@map color 14 to (64, 224, 208), "turquoise"			'
write(grace,*) '@map color 15 to (0, 139, 0), "green4"				'
write(grace,*) '@reference date 0						'
write(grace,*) '@date wrap off							'
write(grace,*) '@date wrap year 1950						'
write(grace,*) '@default linewidth 1.0						'
write(grace,*) '@default linestyle 1						'
write(grace,*) '@default color 1						'
write(grace,*) '@default pattern 1						'
write(grace,*) '@default font 0 						'
write(grace,*) '@default char size 1.000000					'
write(grace,*) '@default symbol size 1.000000					'
write(grace,*) '@default sformat "%.8g" 					'
write(grace,*) '@background color 0						'
write(grace,*) '@page background fill on					'
write(grace,*) '@timestamp off							'
write(grace,*) '@timestamp 0.03, 0.03						'
write(grace,*) '@timestamp color 1						'
write(grace,*) '@timestamp rot 0						'
write(grace,*) '@timestamp font 0						'
write(grace,*) '@timestamp char size 1.000000					'
write(grace,*) '@timestamp def "Mon May 25 03:32:54 2009"			'
write(grace,*) '@r0 off 							'
write(grace,*) '@link r0 to g0							'
write(grace,*) '@r0 type above							'
write(grace,*) '@r0 linestyle 1 						'
write(grace,*) '@r0 linewidth 1.0						'
write(grace,*) '@r0 color 1							'
write(grace,*) '@r0 line 0, 0, 0, 0						'
write(grace,*) '@r1 off 							'
write(grace,*) '@link r1 to g0							'
write(grace,*) '@r1 type above							'
write(grace,*) '@r1 linestyle 1 						'
write(grace,*) '@r1 linewidth 1.0						'
write(grace,*) '@r1 color 1							'
write(grace,*) '@r1 line 0, 0, 0, 0						'
write(grace,*) '@r2 off 							'
write(grace,*) '@link r2 to g0							'
write(grace,*) '@r2 type above							'
write(grace,*) '@r2 linestyle 1 						'
write(grace,*) '@r2 linewidth 1.0						'
write(grace,*) '@r2 color 1							'
write(grace,*) '@r2 line 0, 0, 0, 0						'
write(grace,*) '@r3 off 							'
write(grace,*) '@link r3 to g0							'
write(grace,*) '@r3 type above							'
write(grace,*) '@r3 linestyle 1 						'
write(grace,*) '@r3 linewidth 1.0						'
write(grace,*) '@r3 color 1							'
write(grace,*) '@r3 line 0, 0, 0, 0						'
write(grace,*) '@r4 off 							'
write(grace,*) '@link r4 to g0							'
write(grace,*) '@r4 type above							'
write(grace,*) '@r4 linestyle 1 						'
write(grace,*) '@r4 linewidth 1.0						'
write(grace,*) '@r4 color 1							'
write(grace,*) '@r4 line 0, 0, 0, 0						'
write(grace,*) '@g0 on								'
write(grace,*) '@g0 hidden false						'
write(grace,*) '@g0 type XY							'
write(grace,*) '@g0 stacked false						'
write(grace,*) '@g0 bar hgap 0.000000						'
write(grace,*) '@g0 fixedpoint off						'
write(grace,*) '@g0 fixedpoint type 0						'
write(grace,*) '@g0 fixedpoint xy 0.000000, 0.000000				'
write(grace,*) '@g0 fixedpoint format general general				'
write(grace,*) '@g0 fixedpoint prec 6, 6					'
write(grace,*) '@with g0							'
write(grace,'(a14,f0.6,a,f0.2,a,f0.6)')' @    world 0,', min_e-De,',',max_time,',',max_e+De
write(grace,*) '@    stack world 0, 0, 0, 0					'
write(grace,*) '@    znorm 1							'
write(grace,*) '@    view 0.150000, 0.531818, 0.602279, 0.850000		'
write(grace,*) '@    title ""							'
write(grace,*) '@    title font 0						'
write(grace,*) '@    title size 1.500000					'
write(grace,*) '@    title color 1						'
write(grace,*) '@    subtitle ""						'
write(grace,*) '@    subtitle font 0						'
write(grace,*) '@    subtitle size 1.000000					'
write(grace,*) '@    subtitle color 1						'
write(grace,*) '@    xaxes scale Normal 					'
write(grace,*) '@    yaxes scale Normal 					'
write(grace,*) '@    xaxes invert off						'
write(grace,*) '@    yaxes invert off						'
write(grace,*) '@    xaxis  on							'
write(grace,*) '@    xaxis  type zero false					'
write(grace,*) '@    xaxis  offset 0.000000 , 0.000000				'
write(grace,*) '@    xaxis  bar on						'
write(grace,*) '@    xaxis  bar color 1 					'
write(grace,*) '@    xaxis  bar linestyle 1					'
write(grace,*) '@    xaxis  bar linewidth 1.0					'
write(grace,*) '@    xaxis  label ""						'
write(grace,*) '@    xaxis  label layout para					'
write(grace,*) '@    xaxis  label place auto					'
write(grace,*) '@    xaxis  label char size 1.000000				'
write(grace,*) '@    xaxis  label font 0					'
write(grace,*) '@    xaxis  label color 1					'
write(grace,*) '@    xaxis  label place normal					'
write(grace,*) '@    xaxis  tick on						'
write(grace,*) '@    xaxis  tick major 0.2					'
write(grace,*) '@    xaxis  tick minor ticks 4					'
write(grace,*) '@    xaxis  tick default 6					'
write(grace,*) '@    xaxis  tick place rounded true				'
write(grace,*) '@    xaxis  tick in						'
write(grace,*) '@    xaxis  tick major size 0.75				'
write(grace,*) '@    xaxis  tick major color 1					'
write(grace,*) '@    xaxis  tick major linewidth 1.0				'
write(grace,*) '@    xaxis  tick major linestyle 1				'
write(grace,*) '@    xaxis  tick major grid off 				'
write(grace,*) '@    xaxis  tick minor color 1					'
write(grace,*) '@    xaxis  tick minor linewidth 1.0				'
write(grace,*) '@    xaxis  tick minor linestyle 1				'
write(grace,*) '@    xaxis  tick minor grid off 				'
write(grace,*) '@    xaxis  tick minor size 0.250000				'
write(grace,*) '@    xaxis  ticklabel on					'
write(grace,*) '@    xaxis  ticklabel format general				'
write(grace,*) '@    xaxis  ticklabel prec 5					'
write(grace,*) '@    xaxis  ticklabel formula ""				'
write(grace,*) '@    xaxis  ticklabel append "" 				'
write(grace,*) '@    xaxis  ticklabel prepend ""				'
write(grace,*) '@    xaxis  ticklabel angle 0					'
write(grace,*) '@    xaxis  ticklabel skip 0					'
write(grace,*) '@    xaxis  ticklabel stagger 0 				'
write(grace,*) '@    xaxis  ticklabel place normal				'
write(grace,*) '@    xaxis  ticklabel offset auto				'
write(grace,*) '@    xaxis  ticklabel offset 0.000000 , 0.010000		'
write(grace,*) '@    xaxis  ticklabel start type auto				'
write(grace,*) '@    xaxis  ticklabel start 0.000000				'
write(grace,*) '@    xaxis  ticklabel stop type auto				'
write(grace,*) '@    xaxis  ticklabel stop 0.000000				'
write(grace,*) '@    xaxis  ticklabel char size 0.750000			'
write(grace,*) '@    xaxis  ticklabel font 0					'
write(grace,*) '@    xaxis  ticklabel color 1					'
write(grace,*) '@    xaxis  tick place both					'
write(grace,*) '@    xaxis  tick spec type none 				'
write(grace,*) '@    yaxis  on							'
write(grace,*) '@    yaxis  type zero false					'
write(grace,*) '@    yaxis  offset 0.000000 , 0.000000				'
write(grace,*) '@    yaxis  bar on						'
write(grace,*) '@    yaxis  bar color 1 					'
write(grace,*) '@    yaxis  bar linestyle 1					'
write(grace,*) '@    yaxis  bar linewidth 1.0					'
write(grace,*) '@    yaxis  label "e-kin [a.u.]"				'
write(grace,*) '@    yaxis  label layout para					'
write(grace,*) '@    yaxis  label place auto					'
write(grace,*) '@    yaxis  label char size 1.000000				'
write(grace,*) '@    yaxis  label font 0					'
write(grace,*) '@    yaxis  label color 1					'
write(grace,*) '@    yaxis  label place normal					'
write(grace,*) '@    yaxis  tick on						'
write(grace,*) '@    yaxis  tick major',De,'  				'
write(grace,*) '@    yaxis  tick minor ticks 4					'
write(grace,*) '@    yaxis  tick default 6					'
write(grace,*) '@    yaxis  tick place rounded true				'
write(grace,*) '@    yaxis  tick in						'
write(grace,*) '@    yaxis  tick major size 0.750000				'
write(grace,*) '@    yaxis  tick major color 1					'
write(grace,*) '@    yaxis  tick major linewidth 1.0				'
write(grace,*) '@    yaxis  tick major linestyle 1				'
write(grace,*) '@    yaxis  tick major grid off 				'
write(grace,*) '@    yaxis  tick minor color 1					'
write(grace,*) '@    yaxis  tick minor linewidth 1.0				'
write(grace,*) '@    yaxis  tick minor linestyle 1				'
write(grace,*) '@    yaxis  tick minor grid off 				'
write(grace,*) '@    yaxis  tick minor size 0.250000				'
write(grace,*) '@    yaxis  ticklabel on					'
write(grace,*) '@    yaxis  ticklabel format general				'
write(grace,*) '@    yaxis  ticklabel prec 3					'
write(grace,*) '@    yaxis  ticklabel formula ""				'
write(grace,*) '@    yaxis  ticklabel append "" 				'
write(grace,*) '@    yaxis  ticklabel prepend ""				'
write(grace,*) '@    yaxis  ticklabel angle 0					'
write(grace,*) '@    yaxis  ticklabel skip 0					'
write(grace,*) '@    yaxis  ticklabel stagger 0 				'
write(grace,*) '@    yaxis  ticklabel place normal				'
write(grace,*) '@    yaxis  ticklabel offset auto				'
write(grace,*) '@    yaxis  ticklabel offset 0.000000 , 0.010000		'
write(grace,*) '@    yaxis  ticklabel start type auto				'
write(grace,*) '@    yaxis  ticklabel start 0.000000				'
write(grace,*) '@    yaxis  ticklabel stop type auto				'
write(grace,*) '@    yaxis  ticklabel stop 0.000000				'
write(grace,*) '@    yaxis  ticklabel char size 0.750000			'
write(grace,*) '@    yaxis  ticklabel font 0					'
write(grace,*) '@    yaxis  ticklabel color 1					'
write(grace,*) '@    yaxis  tick place both					'
write(grace,*) '@    yaxis  tick spec type none 				'
write(grace,*) '@    altxaxis  off						'
write(grace,*) '@    altyaxis  off						'
write(grace,*) '@    legend on							'
write(grace,*) '@    legend loctype view					'
write(grace,*) '@    legend 0.85, 0.8						'
write(grace,*) '@    legend box color 1 					'
write(grace,*) '@    legend box pattern 1					'
write(grace,*) '@    legend box linewidth 1.0					'
write(grace,*) '@    legend box linestyle 1					'
write(grace,*) '@    legend box fill color 0					'
write(grace,*) '@    legend box fill pattern 1					'
write(grace,*) '@    legend font 0						'
write(grace,*) '@    legend char size 1.000000					'
write(grace,*) '@    legend color 1						'
write(grace,*) '@    legend length 4						'
write(grace,*) '@    legend vgap 1						'
write(grace,*) '@    legend hgap 1						'
write(grace,*) '@    legend invert false					'
write(grace,*) '@    frame type 0						'
write(grace,*) '@    frame linestyle 1						'
write(grace,*) '@    frame linewidth 1.0					'
write(grace,*) '@    frame color 1						'
write(grace,*) '@    frame pattern 1						'
write(grace,*) '@    frame background color 0					'
write(grace,*) '@    frame background pattern 0 				'
write(grace,*) '@    s0 hidden false						'
write(grace,*) '@    s0 type xy 						'
write(grace,*) '@    s0 symbol 0						'
write(grace,*) '@    s0 symbol size 1.000000					'
write(grace,*) '@    s0 symbol color 1						'
write(grace,*) '@    s0 symbol pattern 1					'
write(grace,*) '@    s0 symbol fill color 1					'
write(grace,*) '@    s0 symbol fill pattern 0					'
write(grace,*) '@    s0 symbol linewidth 1.0					'
write(grace,*) '@    s0 symbol linestyle 1					'
write(grace,*) '@    s0 symbol char 65						'
write(grace,*) '@    s0 symbol char font 0					'
write(grace,*) '@    s0 symbol skip 0						'
write(grace,*) '@    s0 line type 1						'
write(grace,*) '@    s0 line linestyle 1					'
write(grace,*) '@    s0 line linewidth 1.0					'
write(grace,*) '@    s0 line color 1						'
write(grace,*) '@    s0 line pattern 1						'
write(grace,*) '@    s0 baseline type 0 					'
write(grace,*) '@    s0 baseline off						'
write(grace,*) '@    s0 dropline off						'
write(grace,*) '@    s0 fill type 0						'
write(grace,*) '@    s0 fill rule 0						'
write(grace,*) '@    s0 fill color 1						'
write(grace,*) '@    s0 fill pattern 1						'
write(grace,*) '@    s0 avalue off						'
write(grace,*) '@    s0 avalue type 2						'
write(grace,*) '@    s0 avalue char size 1.000000				'
write(grace,*) '@    s0 avalue font 0						'
write(grace,*) '@    s0 avalue color 1						'
write(grace,*) '@    s0 avalue rot 0						'
write(grace,*) '@    s0 avalue format general					'
write(grace,*) '@    s0 avalue prec 3						'
write(grace,*) '@    s0 avalue prepend ""					'
write(grace,*) '@    s0 avalue append ""					'
write(grace,*) '@    s0 avalue offset 0.000000 , 0.000000			'
write(grace,*) '@    s0 errorbar on						'
write(grace,*) '@    s0 errorbar place both					'
write(grace,*) '@    s0 errorbar color 1					'
write(grace,*) '@    s0 errorbar pattern 1					'
write(grace,*) '@    s0 errorbar size 1.000000					'
write(grace,*) '@    s0 errorbar linewidth 1.0					'
write(grace,*) '@    s0 errorbar linestyle 1					'
write(grace,*) '@    s0 errorbar riser linewidth 1.0				'
write(grace,*) '@    s0 errorbar riser linestyle 1				'
write(grace,*) '@    s0 errorbar riser clip off 				'
write(grace,*) '@    s0 errorbar riser clip length 0.100000			'
write(grace,*) '@    s0 comment "Cols 1:2"					'
write(grace,*) '@    s0 legend  ""						'
write(grace,*) '@g1 on								'
write(grace,*) '@g1 hidden false						'
write(grace,*) '@g1 type XY							'
write(grace,*) '@g1 stacked false						'
write(grace,*) '@g1 bar hgap 0.000000						'
write(grace,*) '@g1 fixedpoint off						'
write(grace,*) '@g1 fixedpoint type 0						'
write(grace,*) '@g1 fixedpoint xy 0.000000, 0.000000				'
write(grace,*) '@g1 fixedpoint format general general				'
write(grace,*) '@g1 fixedpoint prec 6, 6					'
write(grace,*) '@with g1							'
write(grace,'(a14,f0.0,a,f0.2,a,f0.0)')' @    world 0,', min_t-Dt,',',max_time,',',max_t+Dt
write(grace,*) '@    stack world 0, 0, 0, 0					'
write(grace,*) '@    znorm 1							'
write(grace,*) '@    view 0.692735, 0.531818, 1.145014, 0.850000		'
write(grace,*) '@    title ""							'
write(grace,*) '@    title font 0						'
write(grace,*) '@    title size 1.500000					'
write(grace,*) '@    title color 1						'
write(grace,*) '@    subtitle ""						'
write(grace,*) '@    subtitle font 0						'
write(grace,*) '@    subtitle size 1.000000					'
write(grace,*) '@    subtitle color 1						'
write(grace,*) '@    xaxes scale Normal 					'
write(grace,*) '@    yaxes scale Normal 					'
write(grace,*) '@    xaxes invert off						'
write(grace,*) '@    yaxes invert off						'
write(grace,*) '@    xaxis  on							'
write(grace,*) '@    xaxis  type zero false					'
write(grace,*) '@    xaxis  offset 0.000000 , 0.000000				'
write(grace,*) '@    xaxis  bar on						'
write(grace,*) '@    xaxis  bar color 1 					'
write(grace,*) '@    xaxis  bar linestyle 1					'
write(grace,*) '@    xaxis  bar linewidth 1.0					'
write(grace,*) '@    xaxis  label ""						'
write(grace,*) '@    xaxis  label layout para					'
write(grace,*) '@    xaxis  label place auto					'
write(grace,*) '@    xaxis  label char size 1.000000				'
write(grace,*) '@    xaxis  label font 0					'
write(grace,*) '@    xaxis  label color 1					'
write(grace,*) '@    xaxis  label place normal					'
write(grace,*) '@    xaxis  tick on						'
write(grace,*) '@    xaxis  tick major 0.2             				'
write(grace,*) '@    xaxis  tick minor ticks 3					'
write(grace,*) '@    xaxis  tick default 6					'
write(grace,*) '@    xaxis  tick place rounded true				'
write(grace,*) '@    xaxis  tick in						'
write(grace,*) '@    xaxis  tick major size 0.75		              	'
write(grace,*) '@    xaxis  tick major color 1					'
write(grace,*) '@    xaxis  tick major linewidth 1.0				'
write(grace,*) '@    xaxis  tick major linestyle 1				'
write(grace,*) '@    xaxis  tick major grid off 				'
write(grace,*) '@    xaxis  tick minor color 1					'
write(grace,*) '@    xaxis  tick minor linewidth 1.0				'
write(grace,*) '@    xaxis  tick minor linestyle 1				'
write(grace,*) '@    xaxis  tick minor grid off 				'
write(grace,*) '@    xaxis  tick minor size 0.260000				'
write(grace,*) '@    xaxis  ticklabel on					'
write(grace,*) '@    xaxis  ticklabel format general				'
write(grace,*) '@    xaxis  ticklabel prec 5					'
write(grace,*) '@    xaxis  ticklabel formula ""				'
write(grace,*) '@    xaxis  ticklabel append "" 				'
write(grace,*) '@    xaxis  ticklabel prepend ""				'
write(grace,*) '@    xaxis  ticklabel angle 0					'
write(grace,*) '@    xaxis  ticklabel skip 0					'
write(grace,*) '@    xaxis  ticklabel stagger 0 				'
write(grace,*) '@    xaxis  ticklabel place normal				'
write(grace,*) '@    xaxis  ticklabel offset auto				'
write(grace,*) '@    xaxis  ticklabel offset 0.000000 , 0.010000		'
write(grace,*) '@    xaxis  ticklabel start type auto				'
write(grace,*) '@    xaxis  ticklabel start 0.000000				'
write(grace,*) '@    xaxis  ticklabel stop type auto				'
write(grace,*) '@    xaxis  ticklabel stop 0.000000				'
write(grace,*) '@    xaxis  ticklabel char size 0.750000			'
write(grace,*) '@    xaxis  ticklabel font 0					'
write(grace,*) '@    xaxis  ticklabel color 1					'
write(grace,*) '@    xaxis  tick place both					'
write(grace,*) '@    xaxis  tick spec type none 				'
write(grace,*) '@    yaxis  on							'
write(grace,*) '@    yaxis  type zero false					'
write(grace,*) '@    yaxis  offset 0.000000 , 0.000000				'
write(grace,*) '@    yaxis  bar on						'
write(grace,*) '@    yaxis  bar color 1 					'
write(grace,*) '@    yaxis  bar linestyle 1					'
write(grace,*) '@    yaxis  bar linewidth 1.0					'
write(grace,*) '@    yaxis  label "T [K]"					'
write(grace,*) '@    yaxis  label layout para					'
write(grace,*) '@    yaxis  label place auto					'
write(grace,*) '@    yaxis  label char size 1.000000				'
write(grace,*) '@    yaxis  label font 0					'
write(grace,*) '@    yaxis  label color 1					'
write(grace,*) '@    yaxis  label place opposite					'
write(grace,*) '@    yaxis  tick on						'
write(grace,*) '@    yaxis  tick major', 50,'					'
write(grace,*) '@    yaxis  tick minor ticks 3					'
write(grace,*) '@    yaxis  tick default 6					'
write(grace,*) '@    yaxis  tick place rounded true				'
write(grace,*) '@    yaxis  tick in						'
write(grace,*) '@    yaxis  tick major size 0.750000				'
write(grace,*) '@    yaxis  tick major color 1					'
write(grace,*) '@    yaxis  tick major linewidth 1.0				'
write(grace,*) '@    yaxis  tick major linestyle 1				'
write(grace,*) '@    yaxis  tick major grid off 				'
write(grace,*) '@    yaxis  tick minor color 1					'
write(grace,*) '@    yaxis  tick minor linewidth 1.0				'
write(grace,*) '@    yaxis  tick minor linestyle 1				'
write(grace,*) '@    yaxis  tick minor grid off 				'
write(grace,*) '@    yaxis  tick minor size 0.250000				'
write(grace,*) '@    yaxis  ticklabel on					'
write(grace,*) '@    yaxis  ticklabel format general				'
write(grace,*) '@    yaxis  ticklabel prec 5					'
write(grace,*) '@    yaxis  ticklabel formula ""				'
write(grace,*) '@    yaxis  ticklabel append "" 				'
write(grace,*) '@    yaxis  ticklabel prepend ""				'
write(grace,*) '@    yaxis  ticklabel angle 0					'
write(grace,*) '@    yaxis  ticklabel skip 0					'
write(grace,*) '@    yaxis  ticklabel stagger 0 				'
write(grace,*) '@    yaxis  ticklabel place normal				'
write(grace,*) '@    yaxis  ticklabel offset auto				'
write(grace,*) '@    yaxis  ticklabel offset 0.000000 , 0.010000		'
write(grace,*) '@    yaxis  ticklabel start type auto				'
write(grace,*) '@    yaxis  ticklabel start 0.000000				'
write(grace,*) '@    yaxis  ticklabel stop type auto				'
write(grace,*) '@    yaxis  ticklabel stop 0.000000				'
write(grace,*) '@    yaxis  ticklabel char size 0.750000			'
write(grace,*) '@    yaxis  ticklabel font 0					'
write(grace,*) '@    yaxis  ticklabel color 1					'
write(grace,*) '@    yaxis  tick place both					'
write(grace,*) '@    yaxis  tick spec type none 				'
write(grace,*) '@    altxaxis  off						'
write(grace,*) '@    altyaxis  off						'
write(grace,*) '@    legend on							'
write(grace,*) '@    legend loctype view					'
write(grace,*) '@    legend 0.5, 0.8						'
write(grace,*) '@    legend box color 1 					'
write(grace,*) '@    legend box pattern 1					'
write(grace,*) '@    legend box linewidth 1.0					'
write(grace,*) '@    legend box linestyle 1					'
write(grace,*) '@    legend box fill color 0					'
write(grace,*) '@    legend box fill pattern 1					'
write(grace,*) '@    legend font 0						'
write(grace,*) '@    legend char size 1.000000					'
write(grace,*) '@    legend color 1						'
write(grace,*) '@    legend length 4						'
write(grace,*) '@    legend vgap 1						'
write(grace,*) '@    legend hgap 1						'
write(grace,*) '@    legend invert false					'
write(grace,*) '@    frame type 0						'
write(grace,*) '@    frame linestyle 1						'
write(grace,*) '@    frame linewidth 1.0					'
write(grace,*) '@    frame color 1						'
write(grace,*) '@    frame pattern 1						'
write(grace,*) '@    frame background color 0					'
write(grace,*) '@    frame background pattern 0 				'
write(grace,*) '@    s0 hidden false						'
write(grace,*) '@    s0 type xy 						'
write(grace,*) '@    s0 symbol 0						'
write(grace,*) '@    s0 symbol size 1.000000					'
write(grace,*) '@    s0 symbol color 4						'
write(grace,*) '@    s0 symbol pattern 1					'
write(grace,*) '@    s0 symbol fill color 4					'
write(grace,*) '@    s0 symbol fill pattern 0					'
write(grace,*) '@    s0 symbol linewidth 1.0					'
write(grace,*) '@    s0 symbol linestyle 1					'
write(grace,*) '@    s0 symbol char 65						'
write(grace,*) '@    s0 symbol char font 0					'
write(grace,*) '@    s0 symbol skip 0						'
write(grace,*) '@    s0 line type 1						'
write(grace,*) '@    s0 line linestyle 1					'
write(grace,*) '@    s0 line linewidth 1.0					'
write(grace,*) '@    s0 line color 4						'
write(grace,*) '@    s0 line pattern 1						'
write(grace,*) '@    s0 baseline type 0 					'
write(grace,*) '@    s0 baseline off						'
write(grace,*) '@    s0 dropline off						'
write(grace,*) '@    s0 fill type 0						'
write(grace,*) '@    s0 fill rule 0						'
write(grace,*) '@    s0 fill color 1						'
write(grace,*) '@    s0 fill pattern 1						'
write(grace,*) '@    s0 avalue off						'
write(grace,*) '@    s0 avalue type 2						'
write(grace,*) '@    s0 avalue char size 1.000000				'
write(grace,*) '@    s0 avalue font 0						'
write(grace,*) '@    s0 avalue color 1						'
write(grace,*) '@    s0 avalue rot 0						'
write(grace,*) '@    s0 avalue format general					'
write(grace,*) '@    s0 avalue prec 3						'
write(grace,*) '@    s0 avalue prepend ""					'
write(grace,*) '@    s0 avalue append ""					'
write(grace,*) '@    s0 avalue offset 0.000000 , 0.000000			'
write(grace,*) '@    s0 errorbar on						'
write(grace,*) '@    s0 errorbar place both					'
write(grace,*) '@    s0 errorbar color 4					'
write(grace,*) '@    s0 errorbar pattern 1					'
write(grace,*) '@    s0 errorbar size 1.000000					'
write(grace,*) '@    s0 errorbar linewidth 1.0					'
write(grace,*) '@    s0 errorbar linestyle 1					'
write(grace,*) '@    s0 errorbar riser linewidth 1.0				'
write(grace,*) '@    s0 errorbar riser linestyle 1				'
write(grace,*) '@    s0 errorbar riser clip off 				'
write(grace,*) '@    s0 errorbar riser clip length 0.100000			'
write(grace,*) '@    s0 comment "Cols 1:4"					'
write(grace,*) '@    s0 legend  ""						'
write(grace,*) '@g2 on								'
write(grace,*) '@g2 hidden false						'
write(grace,*) '@g2 type XY							'
write(grace,*) '@g2 stacked false						'
write(grace,*) '@g2 bar hgap 0.000000						'
write(grace,*) '@g2 fixedpoint off						'
write(grace,*) '@g2 fixedpoint type 0						'
write(grace,*) '@g2 fixedpoint xy 0.000000, 0.000000				'
write(grace,*) '@g2 fixedpoint format general general				'
write(grace,*) '@g2 fixedpoint prec 6, 6					'
write(grace,*) '@with g2							'
write(grace,'(a14,f0.6,a,f0.2,a,f0.6)')' @    world 0,', min_ks-Dks,',',max_time,',',max_ks+Dks
write(grace,*) '@    stack world 0, 0, 0, 0					'
write(grace,*) '@    znorm 1							'
write(grace,*) '@    view 0.150000, 0.150000, 0.602279, 0.468182		'
write(grace,*) '@    title ""							'
write(grace,*) '@    title font 0						'
write(grace,*) '@    title size 1.500000					'
write(grace,*) '@    title color 1						'
write(grace,*) '@    subtitle ""						'
write(grace,*) '@    subtitle font 0						'
write(grace,*) '@    subtitle size 1.000000					'
write(grace,*) '@    subtitle color 1						'
write(grace,*) '@    xaxes scale Normal 					'
write(grace,*) '@    yaxes scale Normal 					'
write(grace,*) '@    xaxes invert off						'
write(grace,*) '@    yaxes invert off						'
write(grace,*) '@    xaxis  on							'
write(grace,*) '@    xaxis  type zero false					'
write(grace,*) '@    xaxis  offset 0.000000 , 0.000000				'
write(grace,*) '@    xaxis  bar on						'
write(grace,*) '@    xaxis  bar color 1 					'
write(grace,*) '@    xaxis  bar linestyle 1					'
write(grace,*) '@    xaxis  bar linewidth 1.0					'
write(grace,*) '@    xaxis  label "time [ps]"					'
write(grace,*) '@    xaxis  label layout para					'
write(grace,*) '@    xaxis  label place auto					'
write(grace,*) '@    xaxis  label char size 1.000000				'
write(grace,*) '@    xaxis  label font 0					'
write(grace,*) '@    xaxis  label color 1					'
write(grace,*) '@    xaxis  label place normal					'
write(grace,*) '@    xaxis  tick on						'
write(grace,*) '@    xaxis  tick major 0.2	 				'
write(grace,*) '@    xaxis  tick minor ticks 4					'
write(grace,*) '@    xaxis  tick default 6					'
write(grace,*) '@    xaxis  tick place rounded true				'
write(grace,*) '@    xaxis  tick in						'
write(grace,*) '@    xaxis  tick major size 0.750000				'
write(grace,*) '@    xaxis  tick major color 1					'
write(grace,*) '@    xaxis  tick major linewidth 1.0				'
write(grace,*) '@    xaxis  tick major linestyle 1				'
write(grace,*) '@    xaxis  tick major grid off 				'
write(grace,*) '@    xaxis  tick minor color 1					'
write(grace,*) '@    xaxis  tick minor linewidth 1.0				'
write(grace,*) '@    xaxis  tick minor linestyle 1				'
write(grace,*) '@    xaxis  tick minor grid off 				'
write(grace,*) '@    xaxis  tick minor size 0.250000				'
write(grace,*) '@    xaxis  ticklabel on					'
write(grace,*) '@    xaxis  ticklabel format general				'
write(grace,*) '@    xaxis  ticklabel prec 5					'
write(grace,*) '@    xaxis  ticklabel formula ""				'
write(grace,*) '@    xaxis  ticklabel append "" 				'
write(grace,*) '@    xaxis  ticklabel prepend ""				'
write(grace,*) '@    xaxis  ticklabel angle 0					'
write(grace,*) '@    xaxis  ticklabel skip 0					'
write(grace,*) '@    xaxis  ticklabel stagger 0 				'
write(grace,*) '@    xaxis  ticklabel place normal				'
write(grace,*) '@    xaxis  ticklabel offset auto				'
write(grace,*) '@    xaxis  ticklabel offset 0.000000 , 0.010000		'
write(grace,*) '@    xaxis  ticklabel start type auto				'
write(grace,*) '@    xaxis  ticklabel start 0.000000				'
write(grace,*) '@    xaxis  ticklabel stop type auto				'
write(grace,*) '@    xaxis  ticklabel stop 0.000000				'
write(grace,*) '@    xaxis  ticklabel char size 0.750000			'
write(grace,*) '@    xaxis  ticklabel font 0					'
write(grace,*) '@    xaxis  ticklabel color 1					'
write(grace,*) '@    xaxis  tick place both					'
write(grace,*) '@    xaxis  tick spec type none 				'
write(grace,*) '@    yaxis  on							'
write(grace,*) '@    yaxis  type zero false					'
write(grace,*) '@    yaxis  offset 0.000000 , 0.000000				'
write(grace,*) '@    yaxis  bar on						'
write(grace,*) '@    yaxis  bar color 1 					'
write(grace,*) '@    yaxis  bar linestyle 1					'
write(grace,*) '@    yaxis  bar linewidth 1.0					'
write(grace,*) '@    yaxis  label "KS energy [a.u.]"				'
write(grace,*) '@    yaxis  label layout para					'
write(grace,*) '@    yaxis  label place auto					'
write(grace,*) '@    yaxis  label char size 1.000000				'
write(grace,*) '@    yaxis  label font 0					'
write(grace,*) '@    yaxis  label color 1					'
write(grace,*) '@    yaxis  label place normal					'
write(grace,*) '@    yaxis  tick on						'
write(grace,*) '@    yaxis  tick major',Dks,'					'
write(grace,*) '@    yaxis  tick minor ticks 4					'
write(grace,*) '@    yaxis  tick default 6					'
write(grace,*) '@    yaxis  tick place rounded true				'
write(grace,*) '@    yaxis  tick in						'
write(grace,*) '@    yaxis  tick major size 0.750000				'
write(grace,*) '@    yaxis  tick major color 1					'
write(grace,*) '@    yaxis  tick major linewidth 1.0				'
write(grace,*) '@    yaxis  tick major linestyle 1				'
write(grace,*) '@    yaxis  tick major grid off 				'
write(grace,*) '@    yaxis  tick minor color 1					'
write(grace,*) '@    yaxis  tick minor linewidth 1.0				'
write(grace,*) '@    yaxis  tick minor linestyle 1				'
write(grace,*) '@    yaxis  tick minor grid off 				'
write(grace,*) '@    yaxis  tick minor size 0.250000				'
write(grace,*) '@    yaxis  ticklabel on					'
write(grace,*) '@    yaxis  ticklabel format general				'
write(grace,*) '@    yaxis  ticklabel prec 5					'
write(grace,*) '@    yaxis  ticklabel formula ""				'
write(grace,*) '@    yaxis  ticklabel append "" 				'
write(grace,*) '@    yaxis  ticklabel prepend ""				'
write(grace,*) '@    yaxis  ticklabel angle 0					'
write(grace,*) '@    yaxis  ticklabel skip 0					'
write(grace,*) '@    yaxis  ticklabel stagger 0 				'
write(grace,*) '@    yaxis  ticklabel place normal				'
write(grace,*) '@    yaxis  ticklabel offset auto				'
write(grace,*) '@    yaxis  ticklabel offset 0.000000 , 0.010000		'
write(grace,*) '@    yaxis  ticklabel start type auto				'
write(grace,*) '@    yaxis  ticklabel start 0.000000				'
write(grace,*) '@    yaxis  ticklabel stop type auto				'
write(grace,*) '@    yaxis  ticklabel stop 0.000000				'
write(grace,*) '@    yaxis  ticklabel char size 0.750000			'
write(grace,*) '@    yaxis  ticklabel font 0					'
write(grace,*) '@    yaxis  ticklabel color 1					'
write(grace,*) '@    yaxis  tick place both					'
write(grace,*) '@    yaxis  tick spec type none 				'
write(grace,*) '@    altxaxis  off						'
write(grace,*) '@    altyaxis  off						'
write(grace,*) '@    legend on							'
write(grace,*) '@    legend loctype view					'
write(grace,*) '@    legend 0.5, 0.8						'
write(grace,*) '@    legend box color 1 					'
write(grace,*) '@    legend box pattern 1					'
write(grace,*) '@    legend box linewidth 1.0					'
write(grace,*) '@    legend box linestyle 1					'
write(grace,*) '@    legend box fill color 0					'
write(grace,*) '@    legend box fill pattern 1					'
write(grace,*) '@    legend font 0						'
write(grace,*) '@    legend char size 1.000000					'
write(grace,*) '@    legend color 1						'
write(grace,*) '@    legend length 4						'
write(grace,*) '@    legend vgap 1						'
write(grace,*) '@    legend hgap 1						'
write(grace,*) '@    legend invert false					'
write(grace,*) '@    frame type 0						'
write(grace,*) '@    frame linestyle 1						'
write(grace,*) '@    frame linewidth 1.0					'
write(grace,*) '@    frame color 1						'
write(grace,*) '@    frame pattern 1						'
write(grace,*) '@    frame background color 0					'
write(grace,*) '@    frame background pattern 0 				'
write(grace,*) '@    s0 hidden false						'
write(grace,*) '@    s0 type xy 						'
write(grace,*) '@    s0 symbol 0						'
write(grace,*) '@    s0 symbol size 1.000000					'
write(grace,*) '@    s0 symbol color 1						'
write(grace,*) '@    s0 symbol pattern 1					'
write(grace,*) '@    s0 symbol fill color 1					'
write(grace,*) '@    s0 symbol fill pattern 0					'
write(grace,*) '@    s0 symbol linewidth 1.0					'
write(grace,*) '@    s0 symbol linestyle 1					'
write(grace,*) '@    s0 symbol char 65						'
write(grace,*) '@    s0 symbol char font 0					'
write(grace,*) '@    s0 symbol skip 0						'
write(grace,*) '@    s0 line type 1						'
write(grace,*) '@    s0 line linestyle 1					'
write(grace,*) '@    s0 line linewidth 1.0					'
write(grace,*) '@    s0 line color 1						'
write(grace,*) '@    s0 line pattern 1						'
write(grace,*) '@    s0 baseline type 0 					'
write(grace,*) '@    s0 baseline off						'
write(grace,*) '@    s0 dropline off						'
write(grace,*) '@    s0 fill type 0						'
write(grace,*) '@    s0 fill rule 0						'
write(grace,*) '@    s0 fill color 1						'
write(grace,*) '@    s0 fill pattern 1						'
write(grace,*) '@    s0 avalue off						'
write(grace,*) '@    s0 avalue type 2						'
write(grace,*) '@    s0 avalue char size 1.000000				'
write(grace,*) '@    s0 avalue font 0						'
write(grace,*) '@    s0 avalue color 1						'
write(grace,*) '@    s0 avalue rot 0						'
write(grace,*) '@    s0 avalue format general					'
write(grace,*) '@    s0 avalue prec 3						'
write(grace,*) '@    s0 avalue prepend ""					'
write(grace,*) '@    s0 avalue append ""					'
write(grace,*) '@    s0 avalue offset 0.000000 , 0.000000			'
write(grace,*) '@    s0 errorbar on						'
write(grace,*) '@    s0 errorbar place both					'
write(grace,*) '@    s0 errorbar color 1					'
write(grace,*) '@    s0 errorbar pattern 1					'
write(grace,*) '@    s0 errorbar size 1.000000					'
write(grace,*) '@    s0 errorbar linewidth 1.0					'
write(grace,*) '@    s0 errorbar linestyle 1					'
write(grace,*) '@    s0 errorbar riser linewidth 1.0				'
write(grace,*) '@    s0 errorbar riser linestyle 1				'
write(grace,*) '@    s0 errorbar riser clip off 				'
write(grace,*) '@    s0 errorbar riser clip length 0.100000			'
write(grace,*) '@    s0 comment "Cols 1:5"					'
write(grace,*) '@    s0 legend  ""						'
write(grace,*) '@g3 on								'
write(grace,*) '@g3 hidden false						'
write(grace,*) '@g3 type XY							'
write(grace,*) '@g3 stacked false						'
write(grace,*) '@g3 bar hgap 0.000000						'
write(grace,*) '@g3 fixedpoint off						'
write(grace,*) '@g3 fixedpoint type 0						'
write(grace,*) '@g3 fixedpoint xy 0.000000, 0.000000				'
write(grace,*) '@g3 fixedpoint format general general				'
write(grace,*) '@g3 fixedpoint prec 6, 6					'
write(grace,*) '@with g3							'
write(grace,'(a14,f0.6,a,f0.2,a,f0.6)')' @    world 0,', min_tot-Dtot,',',max_time,',',max_tot+Dtot
write(grace,*) '@    stack world 0, 0, 0, 0					'
write(grace,*) '@    znorm 1							'
write(grace,*) '@    view 0.692735, 0.150000, 1.145014, 0.468182		'
write(grace,*) '@    title ""							'
write(grace,*) '@    title font 0						'
write(grace,*) '@    title size 1.500000					'
write(grace,*) '@    title color 1						'
write(grace,*) '@    subtitle ""						'
write(grace,*) '@    subtitle font 0						'
write(grace,*) '@    subtitle size 1.000000					'
write(grace,*) '@    subtitle color 1						'
write(grace,*) '@    xaxes scale Normal 					'
write(grace,*) '@    yaxes scale Normal 					'
write(grace,*) '@    xaxes invert off						'
write(grace,*) '@    yaxes invert off						'
write(grace,*) '@    xaxis  on							'
write(grace,*) '@    xaxis  type zero false					'
write(grace,*) '@    xaxis  offset 0.000000 , 0.000000				'
write(grace,*) '@    xaxis  bar on						'
write(grace,*) '@    xaxis  bar color 1 					'
write(grace,*) '@    xaxis  bar linestyle 1					'
write(grace,*) '@    xaxis  bar linewidth 1.0					'
write(grace,*) '@    xaxis  label "time [ps]"					'
write(grace,*) '@    xaxis  label layout para					'
write(grace,*) '@    xaxis  label place auto					'
write(grace,*) '@    xaxis  label char size 1.000000				'
write(grace,*) '@    xaxis  label font 0					'
write(grace,*) '@    xaxis  label color 1					'
write(grace,*) '@    xaxis  label place normal					'
write(grace,*) '@    xaxis  tick on						'
write(grace,*) '@    xaxis  tick major 0.2 					'
write(grace,*) '@    xaxis  tick minor ticks 4					'
write(grace,*) '@    xaxis  tick default 6					'
write(grace,*) '@    xaxis  tick place rounded true				'
write(grace,*) '@    xaxis  tick in						'
write(grace,*) '@    xaxis  tick major size 0.750000				'
write(grace,*) '@    xaxis  tick major color 1					'
write(grace,*) '@    xaxis  tick major linewidth 1.0				'
write(grace,*) '@    xaxis  tick major linestyle 1				'
write(grace,*) '@    xaxis  tick major grid off 				'
write(grace,*) '@    xaxis  tick minor color 1					'
write(grace,*) '@    xaxis  tick minor linewidth 1.0				'
write(grace,*) '@    xaxis  tick minor linestyle 1				'
write(grace,*) '@    xaxis  tick minor grid off 				'
write(grace,*) '@    xaxis  tick minor size 0.250000				'
write(grace,*) '@    xaxis  ticklabel on					'
write(grace,*) '@    xaxis  ticklabel format general				'
write(grace,*) '@    xaxis  ticklabel prec 5					'
write(grace,*) '@    xaxis  ticklabel formula ""				'
write(grace,*) '@    xaxis  ticklabel append "" 				'
write(grace,*) '@    xaxis  ticklabel prepend ""				'
write(grace,*) '@    xaxis  ticklabel angle 0					'
write(grace,*) '@    xaxis  ticklabel skip 0					'
write(grace,*) '@    xaxis  ticklabel stagger 0 				'
write(grace,*) '@    xaxis  ticklabel place normal				'
write(grace,*) '@    xaxis  ticklabel offset auto				'
write(grace,*) '@    xaxis  ticklabel offset 0.000000 , 0.010000		'
write(grace,*) '@    xaxis  ticklabel start type auto				'
write(grace,*) '@    xaxis  ticklabel start 0.000000				'
write(grace,*) '@    xaxis  ticklabel stop type auto				'
write(grace,*) '@    xaxis  ticklabel stop 0.000000				'
write(grace,*) '@    xaxis  ticklabel char size 0.750000			'
write(grace,*) '@    xaxis  ticklabel font 0					'
write(grace,*) '@    xaxis  ticklabel color 1					'
write(grace,*) '@    xaxis  tick place both					'
write(grace,*) '@    xaxis  tick spec type none 				'
write(grace,*) '@    yaxis  on							'
write(grace,*) '@    yaxis  type zero false					'
write(grace,*) '@    yaxis  offset 0.000000 , 0.000000				'
write(grace,*) '@    yaxis  bar on						'
write(grace,*) '@    yaxis  bar color 1 					'
write(grace,*) '@    yaxis  bar linestyle 1					'
write(grace,*) '@    yaxis  bar linewidth 1.0					'
write(grace,*) '@    yaxis  label "Etot [au]"					'
write(grace,*) '@    yaxis  label layout para					'
write(grace,*) '@    yaxis  label place auto					'
write(grace,*) '@    yaxis  label char size 1.000000				'
write(grace,*) '@    yaxis  label font 0					'
write(grace,*) '@    yaxis  label color 1					'
write(grace,*) '@    yaxis  label place opposite					'
write(grace,*) '@    yaxis  tick on						'
write(grace,*) '@    yaxis  tick major',Dtot,'					'
write(grace,*) '@    yaxis  tick minor ticks 4					'
write(grace,*) '@    yaxis  tick default 6					'
write(grace,*) '@    yaxis  tick place rounded true				'
write(grace,*) '@    yaxis  tick in						'
write(grace,*) '@    yaxis  tick major size 0.750000				'
write(grace,*) '@    yaxis  tick major color 1					'
write(grace,*) '@    yaxis  tick major linewidth 1.0				'
write(grace,*) '@    yaxis  tick major linestyle 1				'
write(grace,*) '@    yaxis  tick major grid off 				'
write(grace,*) '@    yaxis  tick minor color 1					'
write(grace,*) '@    yaxis  tick minor linewidth 1.0				'
write(grace,*) '@    yaxis  tick minor linestyle 1				'
write(grace,*) '@    yaxis  tick minor grid off 				'
write(grace,*) '@    yaxis  tick minor size 0.250000				'
write(grace,*) '@    yaxis  ticklabel on					'
write(grace,*) '@    yaxis  ticklabel format general				'
write(grace,*) '@    yaxis  ticklabel prec 5					'
write(grace,*) '@    yaxis  ticklabel formula ""				'
write(grace,*) '@    yaxis  ticklabel append "" 				'
write(grace,*) '@    yaxis  ticklabel prepend ""				'
write(grace,*) '@    yaxis  ticklabel angle 0					'
write(grace,*) '@    yaxis  ticklabel skip 0					'
write(grace,*) '@    yaxis  ticklabel stagger 0 				'
write(grace,*) '@    yaxis  ticklabel place normal				'
write(grace,*) '@    yaxis  ticklabel offset auto				'
write(grace,*) '@    yaxis  ticklabel offset 0.000000 , 0.010000		'
write(grace,*) '@    yaxis  ticklabel start type auto				'
write(grace,*) '@    yaxis  ticklabel start 0.000000				'
write(grace,*) '@    yaxis  ticklabel stop type auto				'
write(grace,*) '@    yaxis  ticklabel stop 0.000000				'
write(grace,*) '@    yaxis  ticklabel char size 0.750000			'
write(grace,*) '@    yaxis  ticklabel font 0					'
write(grace,*) '@    yaxis  ticklabel color 1					'
write(grace,*) '@    yaxis  tick place both					'
write(grace,*) '@    yaxis  tick spec type none 				'
write(grace,*) '@    altxaxis  off						'
write(grace,*) '@    altyaxis  off						'
write(grace,*) '@    legend on							'
write(grace,*) '@    legend loctype view					'
write(grace,*) '@    legend 0.5, 0.8						'
write(grace,*) '@    legend box color 1 					'
write(grace,*) '@    legend box pattern 1					'
write(grace,*) '@    legend box linewidth 1.0					'
write(grace,*) '@    legend box linestyle 1					'
write(grace,*) '@    legend box fill color 0					'
write(grace,*) '@    legend box fill pattern 1					'
write(grace,*) '@    legend font 0						'
write(grace,*) '@    legend char size 1.000000					'
write(grace,*) '@    legend color 1						'
write(grace,*) '@    legend length 4						'
write(grace,*) '@    legend vgap 1						'
write(grace,*) '@    legend hgap 1						'
write(grace,*) '@    legend invert false					'
write(grace,*) '@    frame type 0						'
write(grace,*) '@    frame linestyle 1						'
write(grace,*) '@    frame linewidth 1.0					'
write(grace,*) '@    frame color 1						'
write(grace,*) '@    frame pattern 1						'
write(grace,*) '@    frame background color 0					'
write(grace,*) '@    frame background pattern 0 				'
write(grace,*) '@    s0 hidden false						'
write(grace,*) '@    s0 type xy 						'
write(grace,*) '@    s0 symbol 0						'
write(grace,*) '@    s0 symbol size 1.000000					'
write(grace,*) '@    s0 symbol color 1						'
write(grace,*) '@    s0 symbol pattern 1					'
write(grace,*) '@    s0 symbol fill color 1					'
write(grace,*) '@    s0 symbol fill pattern 0					'
write(grace,*) '@    s0 symbol linewidth 1.0					'
write(grace,*) '@    s0 symbol linestyle 1					'
write(grace,*) '@    s0 symbol char 65						'
write(grace,*) '@    s0 symbol char font 0					'
write(grace,*) '@    s0 symbol skip 0						'
write(grace,*) '@    s0 line type 1						'
write(grace,*) '@    s0 line linestyle 1					'
write(grace,*) '@    s0 line linewidth 1.0					'
write(grace,*) '@    s0 line color 1						'
write(grace,*) '@    s0 line pattern 1						'
write(grace,*) '@    s0 baseline type 0 					'
write(grace,*) '@    s0 baseline off						'
write(grace,*) '@    s0 dropline off						'
write(grace,*) '@    s0 fill type 0						'
write(grace,*) '@    s0 fill rule 0						'
write(grace,*) '@    s0 fill color 1						'
write(grace,*) '@    s0 fill pattern 1						'
write(grace,*) '@    s0 avalue off						'
write(grace,*) '@    s0 avalue type 2						'
write(grace,*) '@    s0 avalue char size 1.000000				'
write(grace,*) '@    s0 avalue font 0						'
write(grace,*) '@    s0 avalue color 1						'
write(grace,*) '@    s0 avalue rot 0						'
write(grace,*) '@    s0 avalue format general					'
write(grace,*) '@    s0 avalue prec 3						'
write(grace,*) '@    s0 avalue prepend ""					'
write(grace,*) '@    s0 avalue append ""					'
write(grace,*) '@    s0 avalue offset 0.000000 , 0.000000			'
write(grace,*) '@    s0 errorbar on						'
write(grace,*) '@    s0 errorbar place both					'
write(grace,*) '@    s0 errorbar color 1					'
write(grace,*) '@    s0 errorbar pattern 1					'
write(grace,*) '@    s0 errorbar size 1.000000					'
write(grace,*) '@    s0 errorbar linewidth 1.0					'
write(grace,*) '@    s0 errorbar linestyle 1					'
write(grace,*) '@    s0 errorbar riser linewidth 1.0				'
write(grace,*) '@    s0 errorbar riser linestyle 1				'
write(grace,*) '@    s0 errorbar riser clip off 				'
write(grace,*) '@    s0 errorbar riser clip length 0.100000			'
write(grace,*) '@    s0 comment "Cols 1:8"					'
write(grace,*) '@    s0 legend  ""						'
write(grace,*) '@    s1 hidden false						'
write(grace,*) '@    s1 type xy 						'
write(grace,*) '@    s1 symbol 0						'
write(grace,*) '@    s1 symbol size 1.000000					'
write(grace,*) '@    s1 symbol color 2						'
write(grace,*) '@    s1 symbol pattern 1					'
write(grace,*) '@    s1 symbol fill color 2					'
write(grace,*) '@    s1 symbol fill pattern 0					'
write(grace,*) '@    s1 symbol linewidth 1.0					'
write(grace,*) '@    s1 symbol linestyle 1					'
write(grace,*) '@    s1 symbol char 65						'
write(grace,*) '@    s1 symbol char font 0					'
write(grace,*) '@    s1 symbol skip 0						'
write(grace,*) '@    s1 line type 1						'
write(grace,*) '@    s1 line linestyle 1					'
write(grace,*) '@    s1 line linewidth 1.0					'
write(grace,*) '@    s1 line color 2						'
write(grace,*) '@    s1 line pattern 1						'
write(grace,*) '@    s1 baseline type 0 					'
write(grace,*) '@    s1 baseline off						'
write(grace,*) '@    s1 dropline off						'
write(grace,*) '@    s1 fill type 0						'
write(grace,*) '@    s1 fill rule 0						'
write(grace,*) '@    s1 fill color 1						'
write(grace,*) '@    s1 fill pattern 1						'
write(grace,*) '@    s1 avalue off						'
write(grace,*) '@    s1 avalue type 2						'
write(grace,*) '@    s1 avalue char size 1.000000				'
write(grace,*) '@    s1 avalue font 0						'
write(grace,*) '@    s1 avalue color 1						'
write(grace,*) '@    s1 avalue rot 0						'
write(grace,*) '@    s1 avalue format general					'
write(grace,*) '@    s1 avalue prec 3						'
write(grace,*) '@    s1 avalue prepend ""					'
write(grace,*) '@    s1 avalue append ""					'
write(grace,*) '@    s1 avalue offset 0.000000 , 0.000000			'
write(grace,*) '@    s1 errorbar on						'
write(grace,*) '@    s1 errorbar place both					'
write(grace,*) '@    s1 errorbar color 2					'
write(grace,*) '@    s1 errorbar pattern 1					'
write(grace,*) '@    s1 errorbar size 1.000000					'
write(grace,*) '@    s1 errorbar linewidth 1.0					'
write(grace,*) '@    s1 errorbar linestyle 1					'
write(grace,*) '@    s1 errorbar riser linewidth 1.0				'
write(grace,*) '@    s1 errorbar riser linestyle 1				'
write(grace,*) '@    s1 errorbar riser clip off 				'
write(grace,*) '@    s1 errorbar riser clip length 0.100000			'
write(grace,*) '@    s1 comment "Cols 1:7"					'
write(grace,*) '@    s1 legend  ""						'
write(grace,*) '@target G0.S0							'
write(grace,*) '@type xy							'
!
!
DO i = 1, n
   WRITE(grace,'(2(f0.8,1X))') evp(1,i), evp(2,i)
END DO
!
write(grace,*) '&             '
write(grace,*) '@target G1.S0 '
write(grace,*) '@type xy      '
!
!
DO i = 1, n
   WRITE(grace,'(2(f0.8,1X))') evp(1,i), evp(4,i)
END DO
!
write(grace,*) '&             '
write(grace,*) '@target G2.S0 '
write(grace,*) '@type xy      '
!
!
DO i = 1, n
   WRITE(grace,'(2(f0.8,1X))') evp(1,i), evp(5,i)
END DO
!
write(grace,*) '&             '
write(grace,*) '@target G3.S0 '
write(grace,*) '@type xy      '
!
!
DO i = 1, n   
   WRITE(grace,'(2(f0.8,1X))') evp(1,i), evp(7,i)
END DO
!
write(grace,*) '&             '
write(grace,*) '@target G3.S1 '
write(grace,*) '@type xy      '
!
!
DO i = 1, n   
   WRITE(grace,'(2(f0.8,1X))') evp(1,i), evp(8,i)
END DO

END IF

 CLOSE(grace)
 
  RETURN
END SUBROUTINE write_evp

SUBROUTINE read_dip(filer,n)
  ! read dipole values from ASCII file
  
  CHARACTER (len=256), INTENT(in)    :: filer
  INTEGER,             INTENT(in)    :: n
  
  CHARACTER (len=256)                :: filew

  INTEGER                            :: unitr, unitw
  INTEGER                            :: i, j 
  REAL(DP)                           :: v(3)
  
  unitr = 10
  unitw = 11
  
  OPEN(unitr, file=filer, status='old', position='rewind')
  OPEN(unitw, status='unknown', position='rewind')
  
  DO i = 1, n
     READ(unitr, *) (v(j), j =1,3)
     WRITE(unitw,'(4(f16.8))') (v(j), j =1,3), SQRT(DOT_PRODUCT(v,v))
  END DO
  
  CLOSE(unitr)
  CLOSE(unitw)
  
END SUBROUTINE read_dip
 

SUBROUTINE tcf(filein,fileout,n,m,time_step,interval,iprint,sph_natoms)

  INTEGER,             INTENT(in)    :: n, m, iprint, interval, sph_natoms
  REAL(DP),            INTENT(in)    :: time_step
  CHARACTER(len=256),  INTENT(in)    :: filein, fileout

  CHARACTER(len=256)                 :: filetcf(sph_natoms)

  REAL(DP)                           :: time_dt, tscal

  REAL(DP)                           :: C(n), Cxx(n), Cyy(n), Czz(n)
  REAL(DP)                           :: Norm, Normxx, Normyy, Normzz
  
  INTEGER                            :: iteration, t0, t, time
  INTEGER                            :: i, j, k, p, count
  INTEGER                            :: unitr, unitw

  INTEGER                            :: num(n)
  REAL(DP)                           :: v(n,m), q(n)   
  CHARACTER(len=2)                   :: name(n)

  unitr    = 10
  unitw    = 20

  time_dt = time_step*iprint*interval
  tscal   = time_dt*ps

  DO p = 1, sph_natoms
     
     OPEN(unitr, file=filein, status='old', position='rewind')

     DO i = 1, p-1
        READ(unitr,*)
     END DO

     count = 0
     num   = 0
     q     = 0.0d0
     v     = 0.0d0
     name  = 'XX'

     DO i = 1, n
        
        READ(unitr,*) num(i), name(i), q(i), (v(i,j), j = 1,m)

     !   WRITE(*,'(i2,x,A2,x,4(f16.8))')  num(i), name(i), q(i), (v(i,j), j = 1,m)        

        count = count + 1

        IF (count < n) THEN
           DO j = 1, sph_natoms-1
              READ(unitr,*)
           END DO
        END IF
     END DO

     CLOSE(unitr)

     IF (p < 10) THEN
        filetcf(p) = TRIM(fileout)//'.'//TRIM(name(1))//ACHAR(48+p)//'.dat'
     ELSE IF  (10 <= p .and. p < 20) THEN
        filetcf(p) = TRIM(fileout)//'.'//TRIM(name(1))//ACHAR(48+1)//ACHAR(48+p-10)//'.dat'
     ELSE IF  (20 <= p .and. p < 30) THEN
        filetcf(p) = TRIM(fileout)//'.'//TRIM(name(1))//ACHAR(48+2)//ACHAR(48+p-20)//'.dat'
     ELSE IF  (30 <= p .and. p < 40) THEN
        filetcf(p) = TRIM(fileout)//'.'//TRIM(name(1))//ACHAR(48+3)//ACHAR(48+p-30)//'.dat'
     ELSE IF  (40 <= p .and. p < 50) THEN
        filetcf(p) = TRIM(fileout)//'.'//TRIM(name(1))//ACHAR(48+4)//ACHAR(48+p-40)//'.dat'
     END IF     
  
     iteration = n
     C    = 0.0D0 ; Cxx = 0.0D0   ; Cyy = 0.0D0    ; Czz = 0.0D0
     Norm = 0.0D0 ; Normxx =0.0D0 ; Normyy = 0.0D0 ; Normzz = 0.0D0
     
     DO t0 = 1, iteration
        
        Norm   = Norm   + dot_product(v(t0,:),v(t0,:))
        Normxx = Normxx + v(t0,1)*v(t0,1)
        Normyy = Normyy + v(t0,2)*v(t0,2)
        Normzz = Normzz + v(t0,3)*v(t0,3)
        
     END DO
     
     k  = 0
     t  = 1
     t0 = 1
     
     OPEN(unitw, file=filetcf(p), status='unknown', position='rewind')

     time_loop: DO time = 1, iteration
        origin_loop: DO t0 = 1, iteration-k                
           
           t = t0 + k
           
           C(time)   = C(time)   + v(t0,1)*v(t,1) + v(t0,2)*v(t,2) +  v(t0,3)*v(t,3)
           Cxx(time) = Cxx(time) + v(t0,1)*v(t,1)
           Cyy(time) = Cyy(time) + v(t0,2)*v(t,2)
           Czz(time) = Czz(time) + v(t0,3)*v(t,3)
           
        END DO origin_loop
        
        k         = k + 1
        C(time)   = C(time)/Norm
        Cxx(time) = Cxx(time)/Normxx
        Cyy(time) = Cyy(time)/Normyy
        Czz(time) = Czz(time)/Normzz
        
        WRITE(unitw,'(f10.4,7f10.4)') (time-1)*tscal, C(time), Cxx(time), Cyy(time), Czz(time), &
              (v(time,j), j = 1,m)
        
     END DO time_loop
     
     CLOSE(unitw)


  END DO
END SUBROUTINE tcf


SUBROUTINE tcffull(filein,fileout,n,m,time_step,interval,iprint)
  
  INTEGER,             INTENT(in)    :: n, m, iprint, interval
  REAL(DP),            INTENT(in)    :: time_step
  CHARACTER(len=256),  INTENT(in)    :: filein, fileout
  
  CHARACTER(len=256)                 :: filetcf
  
  REAL(DP)                           :: time_dt, tscal
  
  REAL(DP)                           :: C(n), Cxx(n), Cyy(n), Czz(n)
  REAL(DP)                           :: Norm, Normxx, Normyy, Normzz
  
  INTEGER                            :: iteration, t0, t, time
  INTEGER                            :: i, j, k, p, count
  INTEGER                            :: unitr, unitw
  
  REAL(DP)                           :: v(n,m), q(n)   
  
  
  unitr    = 10
  unitw    = 20
  
  time_dt = time_step*iprint*interval
  tscal   = time_dt*ps
  
  OPEN(unitr, file=filein, status='old', position='rewind')
  
  DO i = 1, n
     
     READ(unitr,*) (v(i,j), j = 1,m)
     
     !   WRITE(*,'(i2,x,A2,x,4(f16.8))')  num(i), name(i), q(i), (v(i,j), j = 1,m)        
     
  END DO
  
  CLOSE(unitr)
  
  iteration = n
  C    = 0.0D0 ; Cxx = 0.0D0   ; Cyy = 0.0D0    ; Czz = 0.0D0
  Norm = 0.0D0 ; Normxx =0.0D0 ; Normyy = 0.0D0 ; Normzz = 0.0D0
  
  DO t0 = 1, iteration
     
     Norm   = Norm   + dot_product(v(t0,:),v(t0,:))
     Normxx = Normxx + v(t0,1)*v(t0,1)
     Normyy = Normyy + v(t0,2)*v(t0,2)
     Normzz = Normzz + v(t0,3)*v(t0,3)
     
  END DO
  
  k  = 0
  t  = 1
  t0 = 1
  
  OPEN(unitw, file=fileout, status='unknown', position='rewind')
  
  time_loop: DO time = 1, iteration
     origin_loop: DO t0 = 1, iteration-k                
        
        t = t0 + k
        
        C(time)   = C(time)   + v(t0,1)*v(t,1) + v(t0,2)*v(t,2) +  v(t0,3)*v(t,3)
        Cxx(time) = Cxx(time) + v(t0,1)*v(t,1)
        Cyy(time) = Cyy(time) + v(t0,2)*v(t,2)
        Czz(time) = Czz(time) + v(t0,3)*v(t,3)
        
     END DO origin_loop
     
     k         = k + 1
     C(time)   = C(time)/Norm
     Cxx(time) = Cxx(time)/Normxx
     Cyy(time) = Cyy(time)/Normyy
     Czz(time) = Czz(time)/Normzz
     
     WRITE(unitw,'(f10.4,7f10.4)') (time-1)*tscal, C(time), Cxx(time), Cyy(time), Czz(time), &
          (v(time,j), j = 1,m)
     
  END DO time_loop
  
  CLOSE(unitw)
  
END SUBROUTINE tcffull

END MODULE dynpro_mod
