MODULE dynpro_sph
  
  USE kinds,        ONLY : DP
  
  USE dynpro_inp,   ONLY : ios
  USE dynpro_inp,   ONLY : sort_xyz, sort_vel
  USE dynpro_inp,   ONLY : sort_xyz, sort_vel
  USE dynpro_inp,   ONLY : cel, abc, vel

  USE dynpro_inp,   ONLY : langle, langle_spef, ldihed_spef, ldip, ldist, ldist_spef, &
       ang_grid, ang_start, ang_stop, ang_test
  USE dynpro_inp,   ONLY : angle_atom1, angle_atom2, angle_atom3, dist_atom

  
  USE dynpro_inp,   ONLY : sph_spe, sph_vel, sph_xyz, sph_label  
  USE dynpro_inp,   ONLY : sph_m, sph_n, sph_natoms, sph_spe, sph_atom, sph_q, sph_p, sph_atn, sph_rcut

  USE dynpro_inp,   ONLY : syst_ID_lab, syst_ID_spe, label
  USE dynpro_inp,   ONLY : interval, iteration, nat, nsp, time_step, iprint

  USE dynpro_inp,   ONLY : evp, fileclas, filedip, filedip_tcf
  USE dynpro_inp,   ONLY : filemol, fileout, filevac, filepsd, filesort, filesph, filesph_pc, filesph_tmp
  
  USE dynpro_mod,   ONLY : sort, class, write_xyz
  USE dynpro_mod,   ONLY : tcffull
  !
  ! RDF variables/routines
  !
  USE dynpro_inp, ONLY : lrdf, rd_atom, rmax, rmin, nr
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
  USE dynpro_inp, ONLY : lpsd, lvacf, loverlap, lspe

  IMPLICIT NONE
  
  SAVE  
CONTAINS
  
  SUBROUTINE dynpro_sphere   
    
    IMPLICIT NONE

    CHARACTER(256)  :: vacfile
    INTEGER         :: i, j , k
    
    DO i = 1, nsp
       IF (sph_atn(i) == 0) THEN
           WRITE(*,*)'Error -sphere- no default values for this function !'
           WRITE(*,*)'Check the keywords sph_atom/sph_atn'
           STOP
        END IF
     END DO
     
     ALLOCATE(sort_xyz(3,nat*iteration))
     ALLOCATE(sort_vel(3,nat*iteration))
     
     CALL sort(nsp, nat, iteration, sph_atom, syst_ID_spe, &
          cel, abc, vel, sort_xyz, sort_vel)
     
     CALL write_xyz(nat, iteration, label, sort_xyz, evp,  &
          filesort, 22)
     
     ALLOCATE(sph_spe(nsp))
     
     sph_spe(1) = sph_n
     sph_spe(2) = sph_spe(1)*2
     
     DO i = 3, nsp
        k = 1
        DO WHILE (sph_atn(i) /= syst_ID_spe(k,1))
           k = k + 1
        END DO
        sph_spe(i) = syst_ID_spe(k,2)
     END DO
     
     sph_natoms = SUM(sph_spe)

     ALLOCATE(sph_xyz(3,sph_natoms*iteration))
     ALLOCATE(sph_vel(3,sph_natoms*iteration))
     ALLOCATE(sph_label(sph_natoms))
     
     CALL class(nsp, nat, iteration, sph_atom, syst_ID_spe, syst_ID_lab, &
          sort_xyz, sort_vel, sph_atn)
     
     CALL write_xyz(nat, iteration, label, sort_xyz, evp,  &
          fileclas, 23)

     k = 1
     DO i = 1, nsp
        DO j = 1, syst_ID_spe(i,2)
           label(k) = syst_ID_lab(i)
           k = k + 1
        END DO
     END DO

     CALL sphere(sort_xyz, sort_vel, sph_spe, sph_natoms, nsp, nat, iteration, &
          syst_ID_spe, sph_xyz, sph_vel, sph_rcut, syst_ID_lab, sph_label, &
          sph_atom, rmax, sph_q)     

     CALL write_xyz(sph_natoms, iteration, sph_label, sph_xyz, evp,  &
          filesph, 25)

     IF (sph_m > 0 .and. sph_n > 1 .and. sph_m < sph_n) THEN
        CALL sphere_m(sph_label, sph_xyz, sph_n, sph_m, sph_natoms, iteration, evp, filesph_tmp, sph_q) 
     END IF

      IF (sph_p > 0 .and. sph_n > 1 .and. sph_m < sph_n) THEN
        CALL sphere_p(nsp, sph_spe, sph_label, sph_xyz, sph_n, sph_p, sph_natoms, iteration, evp, filesph_pc, sph_q) 
     END IF

     IF (sph_n .eq. 0 .and. langle .eqv. .true.) THEN
        CALL ANGLE1(filesph, 417, nsp, iteration, sph_spe, sph_atn, sph_natoms, sph_label, ang_grid, ang_start, ang_stop, ang_test)
     END IF
     
     IF (sph_n .eq. 0 .and. ldist .eqv. .true.) THEN
        CALL DISTANCE1(filesph, 26, nsp, iteration, sph_spe, sph_atn, sph_natoms, sph_label)
     END IF
     
     IF (sph_n .eq. 0 .and. langle_spef .eqv. .true.) THEN

        ! CISPLATIN
        ! analysis 1
        !
        angle_atom1 = 4
        angle_atom2 = 8
        CALL ANGLE2(filesph, 26, iteration, sph_natoms, sph_label, angle_atom1, angle_atom2)
        angle_atom1 = 2
        angle_atom2 = 8
        CALL ANGLE2(filesph, 26, iteration, sph_natoms, sph_label, angle_atom1, angle_atom2)
        angle_atom1 = 2
        angle_atom2 = 4
        CALL ANGLE2(filesph, 26, iteration, sph_natoms, sph_label, angle_atom1, angle_atom2)
        angle_atom1 = 3
        angle_atom2 = 5
        CALL ANGLE2(filesph, 26, iteration, sph_natoms, sph_label, angle_atom1, angle_atom2)
        angle_atom1 = 7
        angle_atom2 = 8
        CALL ANGLE2(filesph, 26, iteration, sph_natoms, sph_label, angle_atom1, angle_atom2)
        !
        angle_atom1 = 8
        angle_atom2 = 7
        angle_atom3 = 1
        !
        CALL ANGLE3(filesph, 26, iteration, sph_natoms, sph_label, angle_atom1, angle_atom2, angle_atom3)

        angle_atom1 = 8
        angle_atom2 = 6
        angle_atom3 = 1
        !
        CALL ANGLE3(filesph, 26, iteration, sph_natoms, sph_label, angle_atom1, angle_atom2, angle_atom3)        
        !
        ! analysis 2
        !
        !angle_atom1 = 3
        !angle_atom2 = 2
        !CALL ANGLE2(filesph, 26, iteration, sph_natoms, sph_label, angle_atom1, angle_atom2)
  
     END IF

     IF (sph_n .eq. 0 .and. ldihed_spef .eqv. .true.) THEN

        ! CISPLATIN
        ! analysis 1
        !

        !dihed_atom1 = 1
        !dihed_atom2 = 2
        !dihed_atom3 = 3
        !dihed_atom4 = 4

        CALL DIHEDRAL(filesph, 26, iteration, sph_natoms, sph_label, 2, 1, 8, 6)
        CALL DIHEDRAL(filesph, 26, iteration, sph_natoms, sph_label, 2, 1, 8, 7)
  
     END IF
     
     IF (sph_n .eq. 0 .and. ldist_spef .eqv. .true.) THEN
        ! analysis 1
        !
        dist_atom = 2
        CALL DISTANCE2(filesph, 26, iteration, sph_natoms, sph_label, dist_atom)
        dist_atom = 3
        CALL DISTANCE2(filesph, 26, iteration, sph_natoms, sph_label, dist_atom)
        dist_atom = 4
        CALL DISTANCE2(filesph, 26, iteration, sph_natoms, sph_label, dist_atom)
        dist_atom = 5
        CALL DISTANCE2(filesph, 26, iteration, sph_natoms, sph_label, dist_atom)
        dist_atom = 6
        CALL DISTANCE2(filesph, 26, iteration, sph_natoms, sph_label, dist_atom)
        dist_atom = 7
        CALL DISTANCE2(filesph, 26, iteration, sph_natoms, sph_label, dist_atom)
        dist_atom = 8
        CALL DISTANCE2(filesph, 26, iteration, sph_natoms, sph_label, dist_atom)
        !
        ! analysis 2
        !
        !dist_atom = 2
        !CALL DISTANCE2(filesph, 26, iteration, sph_natoms, sph_label, dist_atom)
        !dist_atom = 3
        !CALL DISTANCE2(filesph, 26, iteration, sph_natoms, sph_label, dist_atom)
        
     END IF
     !
     unit_vac = 778
     filevac  = TRIM(fileout)//'.vac_sph'
     vacfile = TRIM(filevac)//'_raw.dat'
     !
     IF (lvacf) THEN
        !
        ALLOCATE(C_tot(iteration))  
        ALLOCATE(C_nsp(nsp,iteration))
        !
        CALL VACF(sph_natoms, iteration, sph_vel, 1, syst_ID_spe, syst_ID_lab, C_nsp, C_tot, filevac, unit_vac, &
             time_step, interval, iprint)
        !
        DEALLOCATE(C_tot)  
        DEALLOCATE(C_nsp)
        !
     END IF
     !
     write(*,*) 
     !       
     IF (lpsd) THEN
        !
        ALLOCATE(C_tot(iteration))  
        !
        open(unit_vac, action="read", iostat=ios, file=vacfile, status='old', position='rewind')
        !
        IF (ios /= 0) THEN
           PRINT *,"Error -DynPro SPH- opening file:  ", vacfile
           STOP
        END IF
        !
        DO i = 1, iteration
           READ(unit_vac,'(f20.10)') C_tot(i)
           !print*, C_tot(i)
        END DO
        !
        close(unit_vac)
        !     
        unit_psd = 889
        filepsd  = TRIM(fileout)//'.psd_sph'
        lspe     = .false.
           
        CALL PSD(0, 1, lspe, iteration, syst_ID_lab, syst_ID_spe, filevac, loverlap, n_seg, win, unit_vac, m_cof, &
             time_step, interval, iprint, unit_psd, filepsd, qcfactor, temp)
        !
        DEALLOCATE(C_tot)         
        !
     END IF
     !   
     !DEALLOCATE(C_tot)  
     !DEALLOCATE(C_nsp)
     !  
     !END IF
     !
     CALL  write_pointcharge(nat,sph_natoms,iteration,filemol,filesph)
     !
     IF (ldip) THEN
        ! 
        write(*,*) syst_ID_lab
        write(*,*) syst_ID_spe        
        write(*,*) sph_natoms
        write(*,*) sph_spe
        write(*,*) sph_label
        !
        ! CALL read_dip(filedip,iteration)
        ! CALL tcf(filespe,filedip_tcf,iteration,3,time_step,interval,iprint,sph_natoms)
        !
        CALL tcffull(filedip,filedip_tcf,iteration,3,time_step,interval,iprint)
        ! 
     END IF
     !
     DEALLOCATE(sort_xyz)
     DEALLOCATE(sort_vel)
     DEALLOCATE(sph_spe)
     DEALLOCATE(sph_xyz)
     DEALLOCATE(sph_label)
     !
   END SUBROUTINE dynpro_sphere
   
   !==---------------------------------------------------------------------=
   
   SUBROUTINE ANGLE1(filename, unit, nsp, n, sph_spe, sph_atn, sph_natoms, sph_label, grid, ang_start, ang_stop, ang_test)
    
     USE constants,  ONLY : pi, eps12
 
     IMPLICIT NONE
     
     INTEGER,        INTENT(in)   :: unit, nsp, n, sph_natoms
     CHARACTER(256), INTENT(in)   :: filename
     CHARACTER(3),   INTENT(in)   :: sph_label(sph_natoms)
     
     INTEGER,        INTENT(in)   :: sph_atn(100)
     INTEGER,        INTENT(in)   :: sph_spe(nsp)
     
     INTEGER,        INTENT(in)   :: grid
     REAL(DP),       INTENT(inout):: ang_start, ang_stop, ang_test
     
     
     REAL(DP)                :: xyz_d(3), g_a(nsp-3,0:grid), g_m(nsp-3)
     REAL(DP)                :: cos_abs(0:grid), mean(nsp-3,n)
     REAL(DP),   ALLOCATABLE :: xyz(:,:)
     REAL(DP)                :: cosine
     REAL(DP)                :: d1, d2, d3, ang, step, ang_mean
     
     INTEGER                 :: unitsph, unitspha, unitsphax
     INTEGER                 :: i, j, k, l, m, o, count
     
     CHARACTER(256)          :: filesph, filespha
     CHARACTER(256)          :: filesphax(nsp-3)         
     CHARACTER(2)            :: atom 
     CHARACTER(6)            :: ang_label(nsp-3)
     
     filesph   = filename
     filespha  = TRIM(filename) // '.angl'
     
     l = 1
     DO k = 4, nsp
        l = l + sph_spe(k)
        ang_label(k-3) = 'Pt_'//TRIM(sph_label(l))
        filesphax(k-3) = TRIM(filename)//'.'//TRIM(ang_label(k-3))// '.angl'
     END DO
     
     unitsph  = unit
     unitspha = unit + 2
     
     g_a    = 0.0d0
     cosine = 0.0d0
     mean   = 0.0d0
     
     ang_start = cos(ang_start*pi/180.0d0)
     ang_stop  = cos(ang_stop*pi/180.0d0)
     ang_test  = cos(ang_test*pi/180.0d0)
     
     step = -(ang_start - ang_stop)/grid
     
     OPEN(unitsph,   file=filesph,  status='old', position='rewind')
     OPEN(unitspha,  file=filespha, status='unknown')
     
     DO k = 1, nsp-3
        unitsphax = unitspha + k
        !print*, unitsphax
        OPEN(unitsphax,  file=filesphax(k), status='unknown')
     END DO
     
     DO m = 1, n
        
        READ(unitsph,*)   
        READ(unitsph,*)
        READ(unitsph,'(a2,3x,3f15.9)') atom, (xyz_d(i), i=1,3)
        
        DO k = 4, nsp
           
           ALLOCATE(xyz(sph_spe(k),3))
           
           DO l = 1, sph_spe(k)
              READ(unitsph,'(a2,3x,3f15.9)') atom, (xyz(l,i),i=1,3)           
           END DO
           
           DO i = 1, sph_spe(k)
              
              DO j = i+1, sph_spe(k)
                 
                 d1 = sqrt((xyz_d(1) - xyz(i,1))**2 &                
                      +   (xyz_d(2)  - xyz(i,2))**2 &
                      +   (xyz_d(3)  - xyz(i,3))**2)
                 
                 d2 = sqrt((xyz_d(1) - xyz(j,1))**2 &                
                      +   (xyz_d(2)  - xyz(j,2))**2 &
                      +   (xyz_d(3)  - xyz(j,3))**2)
                 
                 d3 = dot_product(xyz(i,:),xyz(j,:))
                 
                 cosine = d3/(d1*d2)
                 
                 IF (cosine > ang_test) THEN
                    mean(k-3,m) = mean(k-3,m) + acos(cosine)*180.0d0/pi
                 END IF
                 
                 ang   = ang_start
                 count = 0
                 DO WHILE (cosine > ang)
                    ang   = ang   + step
                    count = count + 1
                 END DO
                 
                 g_a(k-3,count) = g_a(k-3,count) + 1.0d0
                 
              END DO
              
           END DO
           
           DEALLOCATE(xyz)
           
        END DO
        
     END DO
     
     DO k = 1, nsp-3
        DO m = 1, n
           WRITE(unitspha+k,*) mean(k,m) ! /4.0d0
        END DO
     END DO
     
     DO i = 0, grid
        cos_abs(i) = acos(ang_start+(i)*step)*180.0d0/pi
     END DO
     
     DO k = 4, nsp
        g_a(k-3,:)  = g_a(k-3,:)/(sph_spe(k)/2*n)
     END DO
     
     WRITE(unitspha,'(a1,2x,a3,2x,10a10)')   '#', 'cos', (ang_label(i), i=1,(nsp-3))
     
     DO i = 1, grid
        WRITE(unitspha,'(10f10.6)') cos_abs(i), (g_a(k,i), k=1,nsp-3)
     END DO
     
     CLOSE(unitsph)
     CLOSE(unitspha)
     DO k = 1, n
        CLOSE(unitspha+k)
     END DO
     
     RETURN
   END SUBROUTINE ANGLE1

   !==---------------------------------------------------------------------=

   SUBROUTINE ANGLE2(filename, unit, n, sph_natoms, sph_label, atom1, atom2)
     
     USE constants,  ONLY : pi, eps12

     IMPLICIT NONE
     
     INTEGER,        INTENT(in)   :: unit, n, sph_natoms
     CHARACTER(256), INTENT(in)   :: filename
     CHARACTER(3),   INTENT(in)   :: sph_label(sph_natoms)
     INTEGER,        INTENT(in)   :: atom1, atom2
     
     REAL(DP)                :: xyz_d(3)
     REAL(DP),   ALLOCATABLE :: xyz(:,:)
     REAL(DP)                :: cosine, deg
     REAL(DP)                :: d1, d2, d3
     
     INTEGER                 :: unitsph, unitspha, unitspha2
     INTEGER                 :: i, j, k, l, m
     
     CHARACTER(256)          :: filesph, filespha2
     CHARACTER(2)            :: atom 
     CHARACTER(8)            :: tag1, tag2, tagx

     write(tagx,'(I0)') 1
     write(tag1,'(I0)') atom1 
     write(tag2,'(I0)') atom2 

     filesph    = filename
     filespha2  = trim(filename) &
          //'.'//adjustl(trim(sph_label(atom1)))//adjustl(trim(tag1)) &
          //'_'//adjustl(trim(sph_label(  1  )))//adjustl(trim(tagx)) &          
          //'_'//adjustl(trim(sph_label(atom2)))//adjustl(trim(tag2)) &
          //'.ang'
     
     unitsph   = unit
     unitspha2 = unit + 1
     
     cosine = 0.0d0
     
     ALLOCATE(xyz(sph_natoms-1,3))
     
     OPEN(unitsph,    file=filesph,  status='old', position='rewind')
     OPEN(unitspha2,  file=filespha2, status='unknown')
     
     frame_loop: DO m = 1, n
        
        READ(unitsph,*)   
        READ(unitsph,*)
        READ(unitsph,'(a2,3x,3f15.9)') atom, (xyz_d(i), i=1,3)
        !print*, atom, (xyz_d(i), i=1,3)    
        DO j = 1, sph_natoms - 1
           READ(unitsph,'(a2,3x,3f15.9)') atom, (xyz(j,i),i=1,3)       
        !print*, atom, (xyz(j,i), i=1,3)
        END DO
         
        i = atom1-1
        j = atom2-1
        !print*, atom1, atom2
        d1 = sqrt((xyz_d(1) - xyz(i,1))**2 &                
             +   (xyz_d(2)  - xyz(i,2))**2 &
             +   (xyz_d(3)  - xyz(i,3))**2)
        
        d2 = sqrt((xyz_d(1) - xyz(j,1))**2 &                
             +   (xyz_d(2)  - xyz(j,2))**2 &
             +   (xyz_d(3)  - xyz(j,3))**2)
        
        d3 = dot_product(xyz_d(:)-xyz(i,:), xyz_d(:)-xyz(j,:))
        
        cosine = d3/(d1*d2)
        deg    = acos(cosine)*180.0d0/pi	
        
        ! IF (deg < 120.0d0) THEN	
        write(unitspha2,'(f6.2)') deg
        ! END IF  
        
     END DO frame_loop
     
     DEALLOCATE(xyz)
     
     
     CLOSE(unitsph)
     CLOSE(unitspha2)
     
     RETURN
   END SUBROUTINE ANGLE2

   SUBROUTINE ANGLE3(filename, unit, n, sph_natoms, sph_label, atom1, atom2, atom3)
     
     USE constants,  ONLY : pi, eps12

     IMPLICIT NONE
     
     INTEGER,        INTENT(in)   :: unit, n, sph_natoms
     CHARACTER(256), INTENT(in)   :: filename
     CHARACTER(3),   INTENT(in)   :: sph_label(sph_natoms)
     INTEGER,        INTENT(in)   :: atom1, atom2, atom3
     
     REAL(DP)                :: xyz_d(3)
     REAL(DP),   ALLOCATABLE :: xyz(:,:)
     REAL(DP)                :: cosine, deg
     REAL(DP)                :: d1, d2, d3
     
     INTEGER                 :: unitsph, unitspha, unitspha2
     INTEGER                 :: i, j, k, l, m
     
     CHARACTER(256)          :: filesph, filespha2
     CHARACTER(2)            :: atom 
     CHARACTER(8)            :: tag1, tag2, tag3

     write(tag1,'(I0)') atom1
     write(tag2,'(I0)') atom2 
     write(tag3,'(I0)') atom3 

     filesph    = filename
     filespha2  = trim(filename) &
          //'.'//adjustl(trim(sph_label(atom1)))//adjustl(trim(tag1)) &
          //'_'//adjustl(trim(sph_label(atom2)))//adjustl(trim(tag2)) &          
          //'_'//adjustl(trim(sph_label(atom3)))//adjustl(trim(tag3)) &
          //'.ang'
     
     unitsph   = unit
     unitspha2 = unit + 1
     
     cosine = 0.0d0
     
     ALLOCATE(xyz(sph_natoms,3))
     
     OPEN(unitsph,    file=filesph,  status='old', position='rewind')
     OPEN(unitspha2,  file=filespha2, status='unknown')
     
     frame_loop: DO m = 1, n
        
        READ(unitsph,*)   
        READ(unitsph,*)
        !READ(unitsph,'(a2,3x,3f15.9)') atom, (xyz_d(i), i=1,3)
        !print*, atom, (xyz_d(i), i=1,3)    
        DO j = 1, sph_natoms
           READ(unitsph,'(a2,3x,3f15.9)') atom, (xyz(j,i),i=1,3)       
        !print*, atom, (xyz(j,i), i=1,3)
        END DO
         
        !i = atom1
        !j = atom2
        
        !print*, atom1, atom2
        d1 = sqrt((xyz(atom1,1) - xyz(atom2,1))**2 &                
             +   ( xyz(atom1,2) - xyz(atom2,2))**2 &
             +   ( xyz(atom1,3) - xyz(atom2,3))**2)
        
        d2 = sqrt((xyz(atom1,1) - xyz(atom3,1))**2 &                
             +   ( xyz(atom1,2) - xyz(atom3,2))**2 &
             +   ( xyz(atom1,3) - xyz(atom3,3))**2)
        
        d3 = dot_product(xyz(atom1,:) - xyz(atom2,:), xyz(atom1,:) - xyz(atom3,:))
        
        cosine = d3/(d1*d2)
        deg    = acos(cosine)*180.0d0/pi	
        
        ! IF (deg < 120.0d0) THEN	
        write(unitspha2,'(f6.2)') deg
        ! END IF  
        
     END DO frame_loop
     
     DEALLOCATE(xyz)
     
     
     CLOSE(unitsph)
     CLOSE(unitspha2)
     
     RETURN
   END SUBROUTINE ANGLE3
   
   SUBROUTINE DIHEDRAL(filename, unit, n, sph_natoms, sph_label, atom1, atom2, atom3, atom4)
     !http://www.bio.net/bionet/mm/xtal-log/1997-February/002899.html
     USE constants,  ONLY : pi, eps12

     IMPLICIT NONE
     
     INTEGER,        INTENT(in)   :: unit, n, sph_natoms
     CHARACTER(256), INTENT(in)   :: filename
     CHARACTER(3),   INTENT(in)   :: sph_label(sph_natoms)
     INTEGER,        INTENT(in)   :: atom1, atom2, atom3, atom4
     
     REAL(DP)                :: xyz_d(3)
     REAL(DP),   ALLOCATABLE :: xyz(:,:)
     REAL(DP)                :: cosine, deg
     REAL(DP)                :: d12, d23, d34, d13, d14, d24
     REAL(DP)                :: P, Q
     
     INTEGER                 :: unitsph, unitspha, unitspha2
     INTEGER                 :: i, j, k, l, m
     
     CHARACTER(256)          :: filesph, filespha2
     CHARACTER(2)            :: atom 
     CHARACTER(8)            :: tag1, tag2, tag3, tag4

     write(tag1,'(I0)') atom1 
     write(tag2,'(I0)') atom2 
     write(tag3,'(I0)') atom3
     write(tag4,'(I0)') atom4

     filesph    = filename
     filespha2  = trim(filename) &
          //'.'//adjustl(trim(sph_label(atom1)))//adjustl(trim(tag1)) &
          //'_'//adjustl(trim(sph_label(atom2)))//adjustl(trim(tag2)) &          
          //'_'//adjustl(trim(sph_label(atom3)))//adjustl(trim(tag3)) &
          //'_'//adjustl(trim(sph_label(atom4)))//adjustl(trim(tag4)) &          
          //'.dih'
     
     unitsph   = unit
     unitspha2 = unit + 1
     
     cosine = 0.0d0
     
     ALLOCATE(xyz(sph_natoms,3))
     
     OPEN(unitsph,    file=filesph,   status='old', position='rewind')
     OPEN(unitspha2,  file=filespha2, status='unknown')
     
     frame_loop: DO m = 1, n
        
        READ(unitsph,*)   
        READ(unitsph,*)
        !READ(unitsph,'(a2,3x,3f15.9)') atom, (xyz_d(i), i=1,3)
        !
        !print*, atom, (xyz_d(i), i=1,3)
        !
        DO j = 1, sph_natoms
           READ(unitsph,'(a2,3x,3f15.9)') atom, (xyz(j,i),i=1,3)
           !
           !print*, atom, (xyz(j,i), i=1,3)
           !
        END DO
        !
        !print*, atom1, atom2, atom3, atom4
        !
        d12 = sqrt((xyz(atom1,1) - xyz(atom2,1))**2 &                
             +   (  xyz(atom1,2) - xyz(atom2,2))**2 &
             +   (  xyz(atom1,3) - xyz(atom2,3))**2)
        !                   
        d23 = sqrt((xyz(atom2,1) - xyz(atom3,1))**2 &                
             +   (  xyz(atom2,2) - xyz(atom3,2))**2 &
             +   (  xyz(atom2,3) - xyz(atom3,3))**2)
        !                                    
        d13 = sqrt((xyz(atom1,1) - xyz(atom3,1))**2 &                
             +   (  xyz(atom1,2) - xyz(atom3,2))**2 &
             +   (  xyz(atom1,3) - xyz(atom3,3))**2)
        !                                    
        d34 = sqrt((xyz(atom3,1) - xyz(atom4,1))**2 &                
             +   (  xyz(atom3,2) - xyz(atom4,2))**2 &
             +   (  xyz(atom3,3) - xyz(atom4,3))**2)
        !                                    
        d13 = sqrt((xyz(atom1,1) - xyz(atom3,1))**2 &                
             +   (  xyz(atom1,2) - xyz(atom3,2))**2 &
             +   (  xyz(atom1,3) - xyz(atom3,3))**2)
        !                                    
        d14 = sqrt((xyz(atom1,1) - xyz(atom4,1))**2 &                
             +   (  xyz(atom1,2) - xyz(atom4,2))**2 &
             +   (  xyz(atom1,3) - xyz(atom4,3))**2)
        !                                    
        d24 = sqrt((xyz(atom2,1) - xyz(atom4,1))**2 &                
             +   (  xyz(atom2,2) - xyz(atom4,2))**2 &
             +   (  xyz(atom2,3) - xyz(atom4,3))**2)

        !d12 = 2.38d0
	!d23 = 1.481d0
	!d34 = 1.487d0
	!d13 = 3.564d0
	!d14 = 3.619d0
	!d24 = 2.408d0
        !dihedral = 13.2 deg.
        
        !
        P =  d12**2 * ( d23**2 + d34**2 - d24**2 ) + &
             d23**2 * (-d23**2 + d34**2 + d24**2 ) + &
             d13**2 * ( d23**2 - d34**2 + d24**2 ) - &
             2 * d23**2 * d14**2
        !
        Q =  (d12 + d23 + d13) * ( d12 + d23 - d13)  * &
             (d12 - d23 + d13) * (-d12 + d23 + d13 ) * &
             (d23 + d34 + d24) * ( d23 + d34 - d24 ) * &
             (d23 - d34 + d24) * (-d23 + d34 + d24 ) 
        
        
        !d3 = dot_product(xyz_d(:)-xyz(i,:), xyz_d(:)-xyz(j,:))
        
        cosine = P / sqrt(Q)
        !print*, cosine
        deg    = acos(cosine) * 180.0d0 / pi	
        
        ! IF (deg < 120.0d0) THEN	
        write(unitspha2,'(f6.2)') deg
        ! END IF  
        
     END DO frame_loop
     
     DEALLOCATE(xyz)
     
     
     CLOSE(unitsph)
     CLOSE(unitspha2)
     
     RETURN
   END SUBROUTINE DIHEDRAL


   
   !==---------------------------------------------------------------------=
   
   SUBROUTINE DISTANCE1(filename, unit, nsp, n, sph_spe, sph_atn, sph_natoms, sph_label)
     
     IMPLICIT NONE
     
     INTEGER,        INTENT(in)   :: unit, nsp, n, sph_natoms
     CHARACTER(256), INTENT(in)   :: filename
     CHARACTER(3),   INTENT(in)   :: sph_label(sph_natoms)
     
     INTEGER,        INTENT(in)   :: sph_atn(100)
     INTEGER,        INTENT(in)   :: sph_spe(nsp)
     
     
     REAL(DP)                :: dist(nsp-3,n), mean(nsp-3)
     REAL(DP)                :: xyz(3), xyz_d(3), r, d
     
     INTEGER                 :: unitsph, unitsphd
     INTEGER                 :: i, j, k, l
     
     CHARACTER(256)          :: filesph, filesphd
     CHARACTER(2)            :: atom 
     CHARACTER(6)            :: dis_label(nsp-3)
     
     
     filesph  = filename
     filesphd = TRIM(filename) // '.dist'
     
     unitsph  = unit
     unitsphd = unit + 2
     
     OPEN(unitsph,  file=filesph,  status='old', position='rewind')
     OPEN(unitsphd, file=filesphd, status='unknown')
     
     DO j = 1, n
        
        READ(unitsph,*)   
        READ(unitsph,*)
        READ(unitsph,'(a2,3x,3f15.9)') atom, (xyz_d(i), i=1,3)
        
        DO k = 4, nsp
           
           d = 0.0
           r = 0.0
           
           
           DO l = 1, sph_spe(k)
              
              READ(unitsph,'(a2,3x,3f15.9)') atom, (xyz(i),i=1,3)
              
              r = sqrt((xyz_d(1) - xyz(1))**2 &                
                   +   (xyz_d(2) - xyz(2))**2 &
                   +   (xyz_d(3) - xyz(3))**2)
              
              d = d + r
              
           END DO
           
           dist(k-3,j) = d / sph_spe(k)
           
        END DO
        
     END DO
     
     l = 1
     DO i = 4, nsp
        l              = l + sph_spe(i)
        mean(i-3)      = SUM(dist(i-3,:)) / n
        dis_label(i-3) = 'Pt_'//TRIM(sph_label(l))
     END DO
     
     WRITE(unitsphd,'(a1,2x,a5,2x,10a10)') '#', 'frame',(dis_label(i), i=1,(nsp-3))
     WRITE(unitsphd,'(a1,3x,a4,2x,10f10.6)') '#', 'mean' , (mean(i), i=1,(nsp-3))
     
     DO i = 1, n
        WRITE(unitsphd,'(i10,2x,10f10.6)') i, (dist(k,i), k=1,(nsp-3))
     END DO
     
     CLOSE(unitsph)
     CLOSE(unitsphd)
     
     RETURN
   END SUBROUTINE DISTANCE1
   
   !==---------------------------------------------------------------------=
   
   SUBROUTINE DISTANCE2(filename, unit, n, sph_natoms, sph_label, atom1)
     
     IMPLICIT NONE
     
     INTEGER,        INTENT(in)   :: unit, n, sph_natoms, atom1
     CHARACTER(256), INTENT(in)   :: filename
     CHARACTER(3),   INTENT(in)   :: sph_label(sph_natoms)
     
     REAL(DP)                :: xyz(sph_natoms-1,3), xyz_d(3), r
     
     INTEGER                 :: unitsph, unitsphd2
     INTEGER                 :: i, j, k, l
     
     CHARACTER(256)          :: filesph, filesphd2
     CHARACTER(2)            :: atom
     CHARACTER(8)            :: tagx, tagy     
     
     
     filesph  = filename
     !filesphd2 = TRIM(filename) // '.2dist'

     write(tagx,'(I0)') 1
     write(tagy,'(I0)') atom1

     filesph    = filename
     filesphd2  = trim(filename) &
          //'.'//adjustl(trim(sph_label( 1   )))//adjustl(trim(tagx)) &
          //'_'//adjustl(trim(sph_label(atom1)))//adjustl(trim(tagy)) &
          //'.dist'
     
     unitsph   = unit
     unitsphd2 = unit + 1
     
     OPEN(unitsph,   file=filesph,   status='old', position='rewind')
     OPEN(unitsphd2, file=filesphd2, status='unknown')
     
     l = atom1 - 1
     
     DO k = 1, n
        
        READ(unitsph,*)   
        READ(unitsph,*)
        READ(unitsph,'(a2,3x,3f15.9)') atom, (xyz_d(i), i=1,3)
        
        DO j = 1, sph_natoms - 1
           READ(unitsph,'(a2,3x,3f15.9)') atom, (xyz(j,i),i=1,3)       
        END DO
        
        r = sqrt((xyz_d(1) - xyz(l,1))**2 &                
             +   (xyz_d(2) - xyz(l,2))**2 &
             +   (xyz_d(3) - xyz(l,3))**2)
        
        WRITE(unitsphd2,'(f6.4)') r
        
     END DO
       
     CLOSE(unitsph)
     CLOSE(unitsphd2)
     
     RETURN
   END SUBROUTINE DISTANCE2

   !==---------------------------------------------------------------------==
   
   SUBROUTINE sphere(sort_xyz, sort_vel, sph_spe, sph_natoms, nsp, natoms, nframes, &
        syst_ID_spe, sph_xyz, sph_vel, d, syst_ID_lab, sph_label, sph_atom, rmax, sph_q)
     
     
     IMPLICIT NONE
     
     INTEGER,  INTENT(in)        :: nsp, natoms, nframes, sph_natoms, sph_q
     
     REAL(DP), INTENT(in)        :: d, rmax
     REAL(DP), INTENT(in)        :: sort_xyz(3,natoms*nframes)  
     
     
     REAL(DP), INTENT(in)        :: sort_vel(3,natoms*nframes)  
     REAL(DP), INTENT(inout)     :: sph_vel(3,sph_natoms*nframes)
     REAL(DP)                    :: vel(3,natoms) 
     REAL(DP), ALLOCATABLE       :: vel_h(:,:)
     
     
     INTEGER,  INTENT(in)        :: sph_spe(nsp), sph_atom
     INTEGER,  INTENT(in)        :: syst_ID_spe(nsp,2)
     CHARACTER(len=2),INTENT(in) :: syst_ID_lab(nsp) 
     REAL(DP), INTENT(inout)     :: sph_xyz(3,sph_natoms*nframes)
     CHARACTER(len=3),INTENT(out):: sph_label(sph_natoms) 
     
     REAL(DP), ALLOCATABLE       :: sph_h(:,:)
     REAL(DP)                    :: xyz(3,natoms), tmp_vel(3), tmp_xyz(3)  
     REAL(DP)                    :: atom_ID_tmp(3), do_h, ghost_xyz(3)
     REAL(DP)                    :: rOH, rXO, rprobe
     
     CHARACTER(len=6)            :: solvent_name
     CHARACTER(len=6)            :: test
     CHARACTER(len=256)          :: filesphlog, filesphadf, toto
     
     INTEGER                     :: charge(nframes)
     INTEGER                     :: h, o, n_o, n_h, nbh
     INTEGER                     :: n, q, m, p, k, i, j, jj, pp, ii, iii, ty
     INTEGER                     :: count, ghost_count, nbo_count, exc_count
     
     ! solvent description (here H2O) !
     ! should also work for other di-hetero-atomic solvent !
     o    = sph_spe(1)            ! O atom
     h    = sph_spe(2)            ! H atom
     n_o  = sph_spe(1)            ! nb. of O atoms, i.e. H2O molecules
     n_h  = sph_spe(2)            ! nb. of H atoms
     nbh  = 2                     ! nb. of H bond to O
     
     ghost_xyz(:) = 0.123456789d0 ! ghost atom coordinate
     rprobe = rmax/2.0d0 - d      ! max for the probing radius 
     
     filesphlog = 'sph.log'       ! log file for the -sphere- routine
     filesphadf = 'sph.adf'       ! adf xyz file for aiMD averaging
     
     n  = syst_ID_spe(1,2) + syst_ID_spe(2,2)
     
     ii = 1
     DO i = 1, nsp
        DO j = 1, sph_spe(i)
           sph_label(ii) = syst_ID_lab(i)
           ii = ii + 1
        END DO
     END DO
     
     test = ACHAR(48+nbh)
     
     solvent_name = TRIM(syst_ID_lab(2))//TRIM(test)//TRIM(syst_ID_lab(1))
     
     ALLOCATE(sph_h(3,n_h))
     ALLOCATE(vel_h(3,n_h))
     
     p           = 0
     pp          = 0
     
     charge      = sph_q
     nbo_count   = 0
     exc_count   = 0
     ghost_count = 0
     
     OPEN(777, file=filesphlog, status='unknown')
     
     WRITE(777,*)' | logfile for the sphere shell reconstruction > '
     
     frames_loop1: DO q = 1, nframes
        
        xyz(:,:)                               = sort_xyz(:,(q+p):(q*natoms))
        sph_xyz(:,(q+pp):(q+pp+n_o))           = xyz(:,1:n_o)
        sph_xyz(:,(q+pp+n_o+n_h):q*sph_natoms) = xyz(:,n+1:natoms)
        vel(:,:)                               = sort_vel(:,(q+p):(q*natoms))
        sph_vel(:,(q+pp):(q+pp+n_o))           = vel(:,1:n_o)
        sph_vel(:,(q+pp+n_o+n_h):q*sph_natoms) = vel(:,n+1:natoms)
        
        ii    = 1    
        sph_h = 0.0d0
        vel_h = 0.0d0
        O_atom_loop1: DO j = 1, n_o
           
           atom_ID_tmp(1:3) = xyz(1:3,j)
           
           m   = syst_ID_spe(1,2) + 1
           i   = 1
           iii = 1
           H_atom_loop1: DO WHILE(i < nbh + 1)
              
              rOH = sqrt((atom_ID_tmp(1) - xyz(1,m))**2 &                
                   +      (atom_ID_tmp(2) - xyz(2,m))**2 &
                   +      (atom_ID_tmp(3) - xyz(3,m))**2)
              
              rXO = sqrt((atom_ID_tmp(1))**2 &                
                   +      (atom_ID_tmp(2))**2 &
                   +      (atom_ID_tmp(3))**2)
              
              IF (rXO >= rprobe ) THEN
                 
                 WRITE(777,*)'Error -sphere- You are out of the min. periodic image constraint !'
                 WRITE(777,*)'frame', q, 'atom',j, 'r =', rXO
                 WRITE(777,*)'Try to adjust the cutoff-radius or/and the number of solvent molecules'
                 STOP
                 
              END IF
              
              IF (rOH <= d ) THEN
                 
                 tmp_xyz = xyz(:,m)
                 tmp_vel = vel(:,m)
                 
              DO jj = 1, ii                 
                 IF (xyz(1,m) == sph_h(1,jj) .and. xyz(2,m) == sph_h(2,jj) .and. xyz(3,m) == sph_h(3,jj)) THEN
                    WRITE(777,*) '  H exchange -sphere- frame', q, 'atom', j
                    tmp_xyz   = ghost_xyz
                    tmp_vel   =  0.0d0
                    charge(q) = charge(q) - 1
                    
                    exc_count = exc_count + 1
                    
                 END IF
              END DO
              
              sph_h(:,ii) = tmp_xyz(:)
              vel_h(:,ii) = tmp_vel(:)
              
              i   = i  + 1
              ii  = ii + 1          
              
           END IF
           
           iii = iii + 1
           
           IF (iii > syst_ID_spe(2,2)) THEN              
              sph_h(:,ii) = ghost_xyz(:)
              vel_h(:,ii) = 0.0d0
              charge(q)   = charge(q) - 1
              

              IF (abs(charge(q) - sph_q)  == nbh) THEN
                 WRITE(777,*) 'Warning -sphere- isolated ', TRIM(syst_ID_lab(1)),' found !'
                 WRITE(777,*) 'frame',q, 'atom',j

                 nbo_count = nbo_count + 1

              END IF

              WRITE(777,*)'  Ghost atom -sphere- frame', q, 'atom',j, 'charge', charge(q)

              ghost_count = ghost_count + 1

              i  = i  + 1
              ii = ii + 1        
              
           END IF
           
           m   = m + 1
           
        END DO H_atom_loop1
        
     END DO O_atom_loop1
     
     sph_xyz(:,(q+pp+n_o):(q+pp+n_o+n_h-1)) = sph_h(:,1:n_h)
     sph_vel(:,(q+pp+n_o):(q+pp+n_o+n_h-1)) = vel_h(:,1:n_h)
     
     p  = p  + (natoms - 1)
     pp = pp + (sph_natoms - 1)
     
  END DO frames_loop1
  
  CALL  write_adf(sph_natoms, nframes, sph_label, sph_xyz, 888, charge, filesphadf)
  
  close(777)
  
  WRITE(*,*) '==== Sphere parameters ===================================>>'
  WRITE(*,'(a22,2x,i10)')      'sph_atom', sph_atom 
  WRITE(*,'(a22,2x,a6)')       'solvent molecule', solvent_name
  WRITE(*,'(a22,2x,i10)')      'shell size', n_o
  WRITE(*,'(a22,2x,f6.2)')     'distance (Ang.)', d
  WRITE(*,'(a22,2x,f6.2)')     'r_probe (Ang.)', rprobe
  WRITE(*,'(a22,2x,i10)')      'nb. of ghost frame', ghost_count
  WRITE(*,'(a22,2x,i10)')      'nb. of NBO', nbo_count
  WRITE(*,'(a22,2x,i10)')      'nb. of H exchange', exc_count
  DO i = 1, nsp
     WRITE(*,'(a22,2x,a2,i8)') 'atom species, number', syst_ID_lab(i), sph_spe(i)
  END DO
  WRITE(*,*) '==========================================================<<'
  
  DEALLOCATE(sph_h)
  DEALLOCATE(vel_h)
  
END SUBROUTINE sphere

!==---------------------------------------------------------------------==

SUBROUTINE sphere_m(sph_label, sph_xyz, sph_n, sph_m, sph_natoms, iteration, evp, &
     filesph_tmp, sph_q)

  IMPLICIT NONE

  INTEGER,            INTENT(in)   :: sph_n, sph_m, sph_natoms, iteration, sph_q

  REAL(DP),           INTENT(in)   :: sph_xyz(3,sph_natoms*iteration)
  REAL(DP),           INTENT(in)   :: evp(iteration)
  REAL(DP),           ALLOCATABLE  :: sph_xyz_tmp(:,:) 

  CHARACTER(len=3),   INTENT(in)   :: sph_label(sph_natoms)
  CHARACTER(len=3),   ALLOCATABLE  :: sph_label_tmp(:)

  CHARACTER(len=256), INTENT(in)   :: filesph_tmp

  CHARACTER(len=256)               :: filesphadf

  INTEGER                          :: scal1, scal2, n_o, n_h, i
  INTEGER                          :: sph_natoms_tmp
  INTEGER                          :: charge(iteration)

  charge = sph_q

  filesphadf = 'sph.adf_m'

  sph_natoms_tmp = sph_natoms - sph_m*3

  n_o =  sph_n - sph_m
  n_h = (sph_n - sph_m)*2
  
  ALLOCATE(sph_xyz_tmp(3,sph_natoms_tmp*iteration))
  ALLOCATE(sph_label_tmp(sph_natoms_tmp)) 
  
  sph_label_tmp(1:n_o)                    = sph_label(sph_m+1:sph_n)
  sph_label_tmp(n_o+1:n_h+n_o)            = sph_label(sph_n+sph_m*2+1:sph_n*2)
  sph_label_tmp(n_h+n_o+1:sph_natoms_tmp) = sph_label(sph_n*3+1:sph_natoms)        
  
  scal1 = 0
  scal2 = 0
  DO i = 1, iteration

     sph_xyz_tmp(:, 1+scal1:scal1+n_o)                     = sph_xyz(:, scal2+sph_m+1:scal2+sph_n)
     sph_xyz_tmp(:, scal1+n_o+1:scal1+n_h+n_o)             = sph_xyz(:, scal2+sph_n+sph_m*2+1:scal2+sph_n*2)
     sph_xyz_tmp(:, scal1+n_h+n_o+1:scal1+sph_natoms_tmp)  = sph_xyz(:, scal2+sph_n*3+1:scal2+sph_natoms)  
     
     scal1 = scal1 + sph_natoms_tmp
     scal2 = scal2 + sph_natoms
     
  END DO
  
  CALL write_adf(sph_natoms_tmp, iteration, sph_label_tmp, sph_xyz_tmp, 889, charge, filesphadf)

  CALL write_xyz(sph_natoms_tmp, iteration, sph_label_tmp, sph_xyz_tmp, evp,  &
       filesph_tmp, 27)
  
  DEALLOCATE(sph_xyz_tmp)
  DEALLOCATE(sph_label_tmp)

  RETURN

END SUBROUTINE sphere_m

!==---------------------------------------------------------------------== 
 
SUBROUTINE sphere_p(nsp, sph_spe, sph_label, sph_xyz, sph_n, sph_p, sph_natoms, iteration, evp, &
     filesph_tmp,  sph_q)

  IMPLICIT NONE

  INTEGER,            INTENT(in)   :: sph_n, sph_p, sph_natoms, iteration, sph_q, nsp

  INTEGER,            INTENT(in)   :: sph_spe(nsp)

  REAL(DP),           INTENT(in)   :: sph_xyz(3,sph_natoms*iteration)
  REAL(DP),           INTENT(in)   :: evp(iteration)
  REAL(DP),           ALLOCATABLE  :: sph_xyz_tmp(:,:) 



  CHARACTER(len=3),   INTENT(in)   :: sph_label(sph_natoms)
  CHARACTER(len=3),   ALLOCATABLE  :: sph_label_tmp(:)

  CHARACTER(len=256), INTENT(in)   :: filesph_tmp

  CHARACTER(len=256)               :: filesphadf

  INTEGER                          :: scal1, scal2, n_o, n_h, i
  INTEGER                          :: sph_natoms_tmp
  INTEGER                          :: charge(iteration)

  charge = sph_q

  filesphadf = 'sph.adf_p'

  sph_natoms_tmp = sph_natoms - sph_p*3 - SUM(sph_spe(3:nsp))

  n_o =  sph_n - sph_p
  n_h = (sph_n - sph_p)*2
  
  ALLOCATE(sph_xyz_tmp(3,sph_natoms_tmp*iteration))
  ALLOCATE(sph_label_tmp(sph_natoms_tmp)) 
  
  sph_label_tmp(1:n_o)                    = sph_label(sph_p+1:sph_n)
  sph_label_tmp(n_o+1:n_h+n_o)            = sph_label(sph_n+sph_p*2+1:sph_n*2)
  !sph_label_tmp(n_h+n_o+1:sph_natoms_tmp) = sph_label(sph_n*3+1:sph_natoms)        
  
  scal1 = 0
  scal2 = 0
  DO i = 1, iteration

     sph_xyz_tmp(:, 1+scal1:scal1+n_o)                     = sph_xyz(:, scal2+sph_p+1:scal2+sph_n)
     sph_xyz_tmp(:, scal1+n_o+1:scal1+n_h+n_o)             = sph_xyz(:, scal2+sph_n+sph_p*2+1:scal2+sph_n*2)
     !sph_xyz_tmp(:, scal1+n_h+n_o+1:scal1+sph_natoms_tmp)  = sph_xyz(:, scal2+sph_n*3+1:scal2+sph_natoms)  
     
     scal1 = scal1 + sph_natoms_tmp
     scal2 = scal2 + sph_natoms
     
  END DO
  
  CALL write_adf(sph_natoms_tmp, iteration, sph_label_tmp, sph_xyz_tmp, 889, charge, filesphadf)

  CALL write_xyz(sph_natoms_tmp, iteration, sph_label_tmp, sph_xyz_tmp, evp,  &
       filesph_tmp, 27)
  
  DEALLOCATE(sph_xyz_tmp)
  DEALLOCATE(sph_label_tmp)

  RETURN

END SUBROUTINE sphere_p

!==---------------------------------------------------------------------==

SUBROUTINE write_adf(n, m, label, xyz, unit, q, file)
  
  IMPLICIT NONE
  
  INTEGER,           INTENT(in)   :: n, m, unit
  INTEGER,           INTENT(in)   :: q(m)
  REAL(DP),          INTENT(in)   :: xyz(3,n*m)
  CHARACTER(len=3),  INTENT(in)   :: label(n)
  CHARACTER(len=256),INTENT(in)   :: file

  INTEGER                         :: i, j, k

  OPEN(unit, file=file, status='unknown')
    
  DO k = 1, m
     
     WRITE(unit,'(i10)') n
     WRITE(unit,'("frame",1X,I8,1X,"charge",1X,I8)') k, q(k)
     
     DO i = 1, n
        IF (xyz(1, n*(k - 1)+i) == 0.123456789d0) THEN
           WRITE (unit,'(a2)') '::'
        ELSE        
           WRITE (unit,'(a2,3x,3f15.9)') TRIM(label(i)), &
                (xyz(j, n*(k - 1)+i), j = 1,3)
        END IF
     END DO
     
  END DO
  
  CLOSE(unit)


  RETURN
END SUBROUTINE write_adf

!==---------------------------------------------------------------------==

SUBROUTINE write_pointcharge(nat,sphnat,nframes,filemol,filesph)

  IMPLICIT NONE

  INTEGER,           INTENT(in)    :: nat, nframes, sphnat

  CHARACTER(len=256),INTENT(in)    :: filemol, filesph
  CHARACTER(len=3)                 :: label(nat), labelsph(sphnat)
  CHARACTER(len=3), ALLOCATABLE    :: labelpc(:)

  CHARACTER(len=256)               :: filepchxyz, filepchadf, filepchmol

  REAL(DP)                         :: xyz(3,nat),  xyzsph(3,sphnat)
  REAL(DP)                         :: xyzsum1, xyzsum2
  REAL(DP), ALLOCATABLE            :: xyzpc(:,:)

  REAL(DP)                         :: qO, qH, qX

  LOGICAL                          :: pointcharge(nat)

  INTEGER                          :: i, j, k, l, count, pcnat, pcnatadf, ghost
  INTEGER                          :: unit1, unit2, unit3, unit4, unit5

  unit1 = 100
  unit2 = 101
  unit3 = 102
  unit4 = 103
  unit5 = 104

  filepchxyz = 'pch.xyz' ! complementary of the sphere output
  filepchadf = 'pch.adf' ! charge file for ADF
  filepchmol = 'pch.mol' ! charge file for vizu

  pcnatadf = nat - sphnat + 1

  ! TIP3P point charge for H2O
  ! from J.Chem.Phys. 79, 926 (1983)
  qO = - 0.8340d0
  qH =   0.4170d0
  qX =   0.0000d0


  OPEN(unit1, file=filemol, status='unknown',position='rewind')
  OPEN(unit2, file=filesph, status='unknown',position='rewind')
  OPEN(unit3, file=filepchxyz, status='unknown')
  OPEN(unit4, file=filepchadf, status='unknown')
  OPEN(unit5, file=filepchmol, status='unknown')

  DO i = 1, nframes
     
     READ(unit1,*) 
     READ(unit1,*) 
     READ(unit2,*) 
     READ(unit2,*) 

     DO j = 1, nat
        READ(unit1,'(a2,3x,3f15.9)') label(j), (xyz(k,j), k=1,3)                    
     END DO
     
     DO j = 1, sphnat
        READ(unit2,'(a2,3x,3f15.9)')  labelsph(j), (xyzsph(k,j), k=1,3)
     END DO
     
     pointcharge = .true.
     count       = 0
         
     DO j = 1, nat
        DO k = 1, sphnat           
           IF (xyz(1,j) == xyzsph(1,k)  .and. xyz(2,j) == xyzsph(2,k) &
                .and. xyz(3,j) == xyzsph(3,k)) THEN 
              pointcharge(j) = .false.
              count = count + 1
           END IF
        END DO
     END DO
     
     pcnat = nat - count

     WRITE(unit5,*) nat
     WRITE(unit5,*) i

     DO  j = 1, sphnat
        IF (xyzsph(1,j) /= 0.123456789d0 .and. xyzsph(2,j) /= 0.123456789d0 &
             .and. xyzsph(3,j) /= 0.123456789d0) THEN 
           
           WRITE(unit5,'(a2,3x,3f15.9)') labelsph(j), (xyzsph(k,j), k=1,3)
        END IF
     END DO     

     DO k = 1, sphnat
        IF (xyzsph(1,k) == 0.123456789d0 .and. xyzsph(2,k) == 0.123456789d0 &
             .and. xyzsph(3,k) == 0.123456789d0) THEN 
           xyzsph(:,k) = 0.0d0
        END IF
     END DO

     ALLOCATE(labelpc(pcnat))
     ALLOCATE(xyzpc(3,pcnat))

     xyzpc       = 0.0d0
     labelpc     = 'XX'
     count       = 1

     DO j = 1, nat
        IF (pointcharge(j) .eqv. .true.) THEN
           labelpc(count)  = label(j)
           xyzpc(:,count)  = xyz(:,j)             
           count = count + 1
        END IF
     END DO
     
     xyzsum1 = SUM(xyzpc) + SUM(xyzsph)
     xyzsum2 = SUM(xyz)

     WRITE(unit3,*) pcnat
     WRITE(unit3,*) xyzsum1, xyzsum2

     IF (pcnat < pcnatadf) THEN

        WRITE(unit4,*) pcnat + 1
        WRITE(unit4,*) i
     ELSE
        WRITE(unit4,*) pcnat
        WRITE(unit4,*) i
     END IF

     DO  j = 1, pcnat            
        WRITE(unit3,'(a2,3x,3f15.9)') labelpc(j), (xyzpc(k,j), k=1,3)
        WRITE(unit5,'(a2,3x,3f15.9)') 'XX', (xyzpc(k,j), k=1,3)

        IF (labelpc(j) == ' H') THEN      
           WRITE(unit4,'(3f15.9,3x,f7.4)') (xyzpc(k,j), k=1,3), qH
        ELSE
           IF  (labelpc(j) == ' O') THEN
              WRITE(unit4,'(3f15.9,3x,f7.4)') (xyzpc(k,j), k=1,3), qO
           ELSE
              WRITE(unit4,'(3f15.9,3x,f7.4)') (xyzpc(k,j), k=1,3), qX
           END IF
        END IF

     END DO

     IF (pcnat < pcnatadf) THEN
        ghost = pcnatadf - pcnat
        DO j = 1, ghost
            WRITE(unit4,'(a2)') '::'
        END DO        
     END IF

     
     DEALLOCATE(labelpc)
     DEALLOCATE(xyzpc)
          
  END DO
  
  CLOSE(unit1)
  CLOSE(unit2)
  CLOSE(unit3)
  CLOSE(unit4)
  CLOSE(unit5)  

  RETURN
END SUBROUTINE write_pointcharge


END MODULE dynpro_sph
