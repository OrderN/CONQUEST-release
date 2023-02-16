!MIT License
!
!Copyright (c) 2019 DynPro: Dynamical Properties
!
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
!
!#include "f_defs.h"
!
PROGRAM dynpro_postproc
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : bohr => BOHR_RADIUS_ANGS
  USE constants,     ONLY : ry   => AUTOEV
  USE constants,     ONLY : ps   => AU_PS
  USE io_files,      ONLY : outdir
  !
  !USE io_global,     ONLY : io_global_start
  !USE mp_global,     ONLY : mp_global_start, mp_global_end
  !USE mp_global,     ONLY : mp_startup, mp_global_end
  !USE mp,            ONLY : mp_start
  !
  !USE dynpro_mod
  !USE dynpro_rax,   ONLY : tensor, tensor_acf
  !
  !USE dynpro_inp
  !
  USE dynpro_inp, ONLY : ldebug
  !
  ! (old version) xml variables need to be clean
  !
  USE dynpro_inp, ONLY : ityp_xml, label_xml, atm_xml, n_xml, nat_xml, nsp_xml
  !
  USE dynpro_inp, ONLY : ityp, ityp_old, label, label_old, nsp
  USE dynpro_inp, ONLY : atm, tmp_atm, na, nat, syst_ID_spe, syst_ID_lab, atomic_number
  !
  USE dynpro_inp, ONLY : lau_unit, lgrace, lspe
  !
  ! Conquest input variables/routines
  !
  USE dynpro_inp, ONLY : lcq_coordforce, nat_cqinp, nsp_cqinp
  USE dynpro_inp, ONLY : lcq_mdframes
  
  USE dynpro_inp, ONLY : ityp_cqcoord, label_cqinp, atm_cqinp, nat_cqinp
  USE dynpro_cqc, ONLY : dynpro_read_cqinp, dynpro_read_cqCoordForce, dynpro_read_cqmdframes
  !
  ! Centering cell variables/routines
  !
  USE dynpro_inp,  ONLY : ct_atom, lcenter
  USE dynpro_cell, ONLY : center_cell, xyz_to_abc
  !
  ! Time/MD snapshot variables
  !
  USE dynpro_inp, ONLY : time_dt, time_step, iprint, iteration
  !
  ! Sampling variables
  !
  USE dynpro_inp, ONLY : lsample, sp_start, sp_end, sp_n, end, start, interval
  !
  ! RDF variables/routines
  !
  USE dynpro_inp, ONLY : lrdf, rd_atom, rmax, rmin, rint, nr
  USE dynpro_rdf, ONLY : dNab, RDF
  !
  ! VACF variables
  !
  USE dynpro_inp, ONLY : lvacf, unit_vac
  USE dynpro_inp, ONLY : C_nsp, C_tot  
  USE dynpro_psd, ONLY : VACF

  !
  ! PSD/FFT variables/routines
  !
  USE dynpro_inp, ONLY : lpsd, unit_psd
  USE dynpro_inp, ONLY : loverlap, m_cof, n_seg, win, qcfactor, temp
  USE dynpro_psd, ONLY : PSD
  USE dynpro_math,ONLY : check_power2
  !
  ! SUB variables/routines
  !
  USE dynpro_inp, ONLY : lsub
  USE dynpro_sub, ONLY : dynpro_sub_file, dynpro_sub_atom
  !
  ! Sphere variables/routines
  !    
  USE dynpro_inp, ONLY : lsphere  
  USE dynpro_sph, ONLY : dynpro_sphere
  !
  ! Output files
  !
  USE dynpro_inp, ONLY : filevel, fileevp, filefor, filecel, filepos
  USE dynpro_inp, ONLY : filemol, filener, filepsd, filevac, filerdf
  !
  USE dynpro_inp, ONLY : xyz, abc, vel, for, cel, evp
  !
  USE dynpro_mod, ONLY : check_file, check_sampling, read_file
  USE dynpro_mod, ONLY : sampling, tricky, tricky_sort, write_evp, write_xyz
  !
  USE dynpro_set, ONLY : dynpro_read_input !, dynpro_read_xml
  !
  IMPLICIT NONE

  INTEGER                       :: cunit, punit, funit, vunit, eunit

  !==---- CENTER cell --------==
  REAL(DP),        ALLOCATABLE  :: ct_abc(:,:), ct_xyz(:,:)

  !==---- RDF calculation ----==
  INTEGER,         ALLOCATABLE  :: event(:,:,:)
  INTEGER,         ALLOCATABLE  :: dN(:,:,:,:)
  INTEGER                       :: rd_count
  CHARACTER(len=2)              :: rd_lab    
  REAL(DP)                      :: vol
  
  !------ Sampling -----------==
  REAL(DP),        ALLOCATABLE  :: sp_xyz(:,:), sp_vel(:,:), xyz_(:,:)
  REAL(DP),        ALLOCATABLE  :: sp_for(:,:), sp_cel(:,:), sp_evp(:,:)

  INTEGER                       :: i, j, k, l, m, n, o, p, tmp_xml

  
  !==---------------------------------------------------------------------==

#if defined __INTEL
  ! ... Intel compilers v .ge.8 allocate a lot of stack space
  ! ... Stack limit is often small, thus causing SIGSEGV and crash
  CALL remove_stack_limit ( )
#endif  
      
  !  see cprstart.f90 for the meaning of the following 4 calls
  !CALL mp_start()
  !CALL mp_env( nproc, mpime, world )
  !CALL mp_global_start( root, mpime, world, nproc )
  !CALL io_global_start( mpime, root )
  !CALL get_env( 'ESPRESSO_TMPDIR', outdir )

  !  initialize mpi
  !CALL mp_start( nproc, mpime, world )
  !CALL mp_global_start( root, mpime, world, nproc )
  !CALL io_global_start( mpime, root )
  !  initialize mpi 
  !CALL mp_startup  ( )                                                             

  CALL get_env( 'ESPRESSO_TMPDIR', outdir )

  IF ( TRIM( outdir ) == ' ' ) outdir = './tmp/'

  !==---------------------------------------------------------------------==
  if (ldebug) write(*,*) 'call dynpro_read_cqinp() '
  CALL dynpro_read_cqinp( )
  
  if (ldebug) write(*,*) 'call dynpro_read_input() '  
  CALL dynpro_read_input( )
  
  if (lcq_coordforce) then
     if (ldebug) write(*,*) 'call dynpro_read_cqCoordForce() '
     CALL dynpro_read_cqCoordForce( )
  else if (lcq_mdframes) then
     if (ldebug) write(*,*) 'call dynpro_read_cqmdframes() '
     CALL dynpro_read_cqmdframes( )
  end if

  if (ldebug) print*, 'call check_file() '
  CALL check_file(filecel,3      ,1,i)
  CALL check_file(filepos,nat_xml,1,j)
  CALL check_file(filevel,nat_xml,1,k)
  CALL check_file(filefor,nat_xml,1,l)
  CALL check_file(fileevp,1      ,0,m)

  tmp_xml = MODULO((i+j+k+l+m), 5)

  n_xml   = i

  nat_xml = nat_cqinp
  nsp_xml = nsp_cqinp

  ALLOCATE(ityp_xml(nat_xml))
  ALLOCATE(label_xml(nat_xml))
  ALLOCATE(atm_xml(nsp_xml))
  
  ityp_xml  = ityp_cqcoord
  label_xml = label_cqinp
  atm_xml   = atm_cqinp

  if (ldebug) print*, '    n steps = ',     n_xml  
  if (ldebug) print*, '    n atoms = ',   nat_xml
  if (ldebug) print*, '    species = ',   nsp_xml
  if (ldebug) print*, '    atoms   = ',   atm_xml
  if (ldebug) print*, '    label   = ', label_xml
  if (ldebug) print*, '    ityp    = ',  ityp_xml
  
  !==---------------------------------------------------------------------==

  !<<< variables settings >>>
  iteration = n_xml         
  nat       = nat_xml       
  nsp       = nsp_xml       
  
  ALLOCATE(xyz(3,nat*iteration))
  ALLOCATE(abc(3,nat*iteration))
  ALLOCATE(vel(3,nat*iteration))
  ALLOCATE(for(3,nat*iteration))
  ALLOCATE(cel(3,3*iteration))
  ALLOCATE(evp(8,iteration))  
  ALLOCATE(ct_abc(3,nat*iteration))
  ALLOCATE(ct_xyz(3,nat*iteration))

  if (ldebug) print*, 'call dynpro_files( ) '  
  punit = 10 
  CALL read_file(filepos, punit, 1, 3, nat, iteration, xyz)
  vunit = 11
  CALL read_file(filevel, vunit, 1, 3, nat, iteration, vel)
  funit = 12 
  CALL read_file(filefor, funit, 1, 3, nat, iteration, for)
  cunit = 13 
  CALL read_file(filecel, cunit, 1, 3, 3,  iteration, cel)
  eunit = 14
  CALL read_file(fileevp, eunit, 0, 8, 1,  iteration, evp)

  !==---------------------------------------------------------------------==
  !<<< sampling function >>>
  IF (lsample) THEN
     
     if (ldebug) print*, 'call check_sampling( ) '  
     CALL check_sampling(sp_start, sp_end, evp, iteration, iprint, &
          sp_n, start, end, interval, time_step, time_dt)

     ALLOCATE(sp_xyz(3,nat*iteration))
     ALLOCATE(sp_vel(3,nat*iteration))
     ALLOCATE(sp_for(3,nat*iteration))
     ALLOCATE(sp_cel(3,3*iteration))
     ALLOCATE(sp_evp(8,iteration))  
     
     sp_xyz = xyz
     sp_vel = vel
     sp_for = for
     sp_cel = cel 
     sp_evp = evp
     
     DEALLOCATE(xyz, abc, vel, for, cel, evp)
     DEALLOCATE(ct_abc,ct_xyz)     
 
     ALLOCATE(xyz(3,nat*sp_n))
     ALLOCATE(abc(3,nat*sp_n))
     ALLOCATE(vel(3,nat*sp_n))
     ALLOCATE(for(3,nat*sp_n))
     ALLOCATE(cel(3,3*sp_n))
     ALLOCATE(evp(8,sp_n))
     ALLOCATE(ct_abc(3,nat*sp_n))
     ALLOCATE(ct_xyz(3,nat*sp_n))

     if (ldebug) print*, 'call sampling( ) '  
     CALL sampling(sp_xyz, 3, nat, iteration, sp_n, start, end, & 
          interval, iprint, xyz)
     CALL sampling(sp_cel, 3, 3, iteration, sp_n, start, end, & 
          interval, iprint, cel)
     CALL sampling(sp_vel, 3, nat, iteration, sp_n, start, end, & 
          interval, iprint, vel)
     CALL sampling(sp_for, 3, nat, iteration, sp_n, start, end, & 
          interval, iprint, for)
     CALL sampling(sp_evp, 8, 1, iteration, sp_n, start, end, & 
          interval, iprint, evp)

     iteration = sp_n
     
     DEALLOCATE(sp_xyz)
     DEALLOCATE(sp_vel)
     DEALLOCATE(sp_for)
     DEALLOCATE(sp_cel)
     DEALLOCATE(sp_evp)  

     ELSE
        interval = 1
  END IF   

  !==---------------------------------------------------------------------==
  ! Add distance / angle analysis here

  !==---------------------------------------------------------------------==

  !<<< variables settings >>>
  IF (rmax == 0.0d0) THEN
     rmax = cel(1,1)
  END IF

  IF( .not. lau_unit) THEN
     cel  = cel*bohr 
     xyz  = xyz*bohr
     for  = for*ry/bohr
     rmax = rmax*bohr
     rmin = rmin*bohr
     !rint = rint*bohr
     vel  = vel*bohr/ps
  END IF

  ALLOCATE(ityp(nat))
  ALLOCATE(ityp_old(nat))

  ALLOCATE(label(nat))
  ALLOCATE(label_old(nat))

 
  ALLOCATE(atm(nsp))
  ALLOCATE(tmp_atm(nsp))
  ALLOCATE(na(nsp))
  ALLOCATE(atomic_number(nsp))
  ALLOCATE(syst_ID_spe(nsp,2))
  ALLOCATE(syst_ID_lab(nsp)) 

  ityp    = ityp_xml
  atm     = atm_xml   

  if (ldebug) print*, 'dynpro: atm'
  if (ldebug) print*, atm
  !
  if (ldebug) print*, 'dynpro: ityp'
  if (ldebug) print*, ityp
  !
  print*, 'call tricky( ) '
  !
  CALL tricky(nsp, nat, atm, ityp, label, syst_ID_spe, syst_ID_lab)

  if (ldebug) print*, 'dynpro before tricky_sort: label_xml'
  if (ldebug) print*, label_xml
  !
  if (ldebug) print*, 'dynpro  before tricky_sort: ityp_xml'
  if (ldebug) print*, ityp_xml
  !
  if (ldebug) print*, 'call tricky_sort( ) '
  !
  CALL tricky_sort(nsp, nat, iteration, xyz, vel, ityp_xml, ityp, syst_ID_spe, &
       label_xml, syst_ID_lab, .true.)

  if (ldebug) print*, 'dynpro: label'
  if (ldebug) print*, label
  !
  if (ldebug) print*, 'dynpro: ityp'
  if (ldebug) print*, ityp
  
  ALLOCATE(xyz_(3,nat))

  xyz_(:,:) = xyz(:,1:nat)
  
  CALL dynpro_sub_file(xyz_)

  DEALLOCATE(xyz_)

  IF (lsub) THEN
     if (ldebug) print*, 'call dynpro_sub_atom( ) '  
     CALL dynpro_sub_atom()
     
  END IF
  
  !==---------------------------------------------------------------------==  
  print*, 'call xyz_to_abc( ) '  
  CALL xyz_to_abc(nat, iteration, xyz, abc, cel)

  IF (lcenter) THEN
     if (ldebug) print*, 'entering CENTER '            
     if (ldebug) print*, '     call center_cell( ) '       
     CALL center_cell(nat, iteration, ct_atom, cel, abc, ct_abc, ct_xyz)
     if (ldebug) print*, '     call write_xyz( ) '       
     CALL write_xyz(nat, iteration, label, ct_xyz, evp, filemol, 21)
     if (ldebug) print*, 'leaving CENTER '   
  ELSE
     if (ldebug) print*, ' call write_xyz( ) '       
     CALL write_xyz(nat, iteration, label, xyz, evp, filemol, 21)
     if (ldebug) print*, ' leaving CENTER '     
  END IF
  
  if (ldebug) print*, 'call write_evp( ) '       
  CALL write_evp(evp, iteration, time_step, iprint, lgrace, filener)

  !==---------------------------------------------------------------------==

  IF (lrdf) THEN
     
     if (ldebug) print*, 'entering RDF ',  rmin, rmax, rint, nr      

     rd_count = 0
     k        = 1
     DO i = 1, nsp     
        IF (rd_atom == syst_ID_spe(i,1)) THEN
           rd_count = syst_ID_spe(i,2)
           rd_lab   = syst_ID_lab(i)
        END IF
     END DO
     
     ALLOCATE(dN(iteration,rd_count,nsp,nr))
     ALLOCATE(event(iteration,rd_count,nsp))

     if (ldebug) print*, '   call dNab '   
     CALL dNab (nsp, nat,  iteration, nr, rmin, rmax,  & 
          ityp, syst_ID_spe, cel, abc, &
          rd_atom, rd_count, dN, event, vol)
     if (ldebug) print*, '   call RDF '        
     CALL RDF (dN, event, iteration, rd_atom, rd_count, rd_lab, &
          nat, nsp, nr, rmin, rmax, rint, vol, & 
          syst_ID_spe, syst_ID_lab, filerdf,  lau_unit)
     DEALLOCATE(dN)
     DEALLOCATE(event)
     if (ldebug) print*, 'leaving RDF '     
  END IF


  !==---------------------------------------------------------------------==

  IF (lvacf) THEN
     if (ldebug) print*, 'entering VACF ' 
     unit_vac = 777
     unit_psd = 888      
     !
     ALLOCATE(C_nsp(nsp,iteration))
     ALLOCATE(C_tot(iteration))  
     !
     if (ldebug) print*, '   call VACF ' 
     !
     CALL VACF(nat, iteration, vel, nsp, syst_ID_spe, syst_ID_lab, &
          C_nsp, C_tot, filevac, unit_vac, time_step, interval, iprint)
     !
     DEALLOCATE(C_nsp)
     DEALLOCATE(C_tot)
     !
     if (ldebug) print*, 'leaving VACF ' 
     !
  END IF
  !
  IF (lpsd) THEN
     !
     if (ldebug) print*, '   call PSD '
     if (ldebug) print*, '      call check_power2'
     !
     CALL check_power2(iteration)
     !
     write(*,*) m_cof
     CALL PSD(0, nsp, lspe, iteration, syst_ID_lab, syst_ID_spe, &
          filevac, loverlap, n_seg, win, unit_vac, m_cof, &
          time_step, interval, iprint, unit_psd, filepsd, qcfactor, temp)
     !  
  END IF
  !    
  !DEALLOCATE(C_nsp)
  !DEALLOCATE(C_tot)
  !
 
!!$  IF (lpsd) THEN
!!$     print*, '   call PSD '
!!$     
!!$     print*, '      call check_power2'
!!$     
!!$     CALL check_power2(iteration)
!!$     
!!$     CALL PSD(0, C_tot, nsp, lspe, iteration, syst_ID_lab, syst_ID_spe, &
!!$          filevac, loverlap, n_seg, win, unit_vac, m_cof, &
!!$          time_step, interval, iprint, unit_psd, filepsd, qcfactor, temp)
!!$
!!$  END IF

  !==---------------------------------------------------------------------==

  IF (lsphere) THEN
     if (ldebug) print*, 'entering SPHERE ' 
     CALL dynpro_sphere()
     if (ldebug) print*, 'leaving SPHERE ' 
  END IF
  
  !==---------------------------------------------------------------------==
  
  !IF (lnmr .and. lsample) THEN
  !   CALL tensor()
  !   CALL tensor_acf()
  !END IF

  !==---------------------------------------------------------------------==
  DEALLOCATE(label_xml)

  DEALLOCATE(ityp_xml)
  DEALLOCATE(atm_xml)
  
  DEALLOCATE(xyz)
  DEALLOCATE(abc)
  DEALLOCATE(vel)
  DEALLOCATE(for)
  DEALLOCATE(cel)
  DEALLOCATE(evp)  
  DEALLOCATE(ct_abc)
  DEALLOCATE(ct_xyz)

  DEALLOCATE(ityp)
  DEALLOCATE(label)
  DEALLOCATE(atm)
  DEALLOCATE(tmp_atm)
  DEALLOCATE(na)
  DEALLOCATE(atomic_number)
  DEALLOCATE(syst_ID_spe)
  DEALLOCATE(syst_ID_lab)   

  !CALL mp_global_end()

END PROGRAM dynpro_postproc



