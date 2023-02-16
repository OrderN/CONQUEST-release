MODULE dynpro_cqc

  IMPLICIT NONE

CONTAINS
  
  SUBROUTINE dynpro_read_cqinp( )

    USE dynpro_inp, ONLY : nat_cqinp, nsp_cqinp, atm_cqinp, ityp_cqcoord, label_cqinp, cel_cqinp
    
    IMPLICIT NONE

    REAL(KIND=8),     ALLOCATABLE :: abc(:,:)
    REAL(KIND=8),     ALLOCATABLE :: ijk(:,:)
    REAL(KIND=8),     ALLOCATABLE :: xyz(:,:)
    REAL(KIND=8),     ALLOCATABLE :: for(:,:)
    REAL(KIND=8),     ALLOCATABLE :: nrg(:)
    REAL(KIND=8),     ALLOCATABLE :: maxf(:)
    REAL(KIND=8),     ALLOCATABLE :: rmsf(:)
    REAL(KIND=8),     ALLOCATABLE :: resf(:)
    INTEGER,          ALLOCATABLE :: ind(:)

    INTEGER,          ALLOCATABLE :: species_ind(:)
    REAL(KIND=8),     ALLOCATABLE :: species_mass(:)
    CHARACTER(len=2), ALLOCATABLE :: species_name(:)

    INTEGER                       :: species_number
    LOGICAL                       :: frac_coord
    
    CHARACTER(len=16)       :: conq_updatom_dat = 'Conquest_coord'
    CHARACTER(len=16)       :: conq_updatom_for = 'UpdatedAtoms.for'
    CHARACTER(len=16)       :: conq_updatom_mld = 'UpdatedAtoms.mld'
    CHARACTER(len=16)       :: conq_updatom_xyz = 'UpdatedAtoms.xyz'
    
    CHARACTER(len=14)       :: conq_in = 'Conquest_input'
    CHARACTER(len=256)      :: sent_conq    = ''
    CHARACTER(len=27)       :: sent_species = '%block ChemicalSpeciesLabel'
    CHARACTER(len=9)        :: sent_end     = '%endblock'
    CHARACTER(len=23)       :: sent_general_species  = 'General.NumberOfSpecies'
    CHARACTER(len=23)       :: sent_general_species_ = 'general.numberofspecies'    
    CHARACTER(len=25)       :: sent_general_fracoor  = 'IO.FractionalAtomicCoords'
    CHARACTER(len=25)       :: sent_general_fracoor_ = 'io.fractionalatomiccoords'
    
    CHARACTER(len=256)      :: str ! To convert integer into strings

    
    CHARACTER(len=12)       :: conq_out     = 'Conquest_out'

    CHARACTER(len=5)       :: tmp_char

    CHARACTER(len=15), PARAMETER ::  sent_mld_fmt         = '[Molden Format]'
    CHARACTER(len=9),  PARAMETER ::  sent_mld_geoconv     = '[GEOCONV]'
    CHARACTER(len=6),  PARAMETER ::  sent_mld_energy      = 'energy'
    CHARACTER(len=9),  PARAMETER ::  sent_mld_maxf        = 'max-force'
    CHARACTER(len=9),  PARAMETER ::  sent_mld_rmsf        = 'rms-force'
    CHARACTER(len=9),  PARAMETER ::  sent_mld_resf        = 'rms-step'
    CHARACTER(len=16), PARAMETER ::  sent_mld_geoxyz      = '[GEOMETRIES] XYZ'
    CHARACTER(len=9),  PARAMETER ::  sent_mld_scf         = 'scf done:'
    CHARACTER(len=8),  PARAMETER ::  sent_mld_forces      = '[FORCES]'
    CHARACTER(len=5),  PARAMETER ::  sent_mld_point       = 'point'

    REAL(KIND=8), PARAMETER :: a_0 = 0.52917720859d0

    INTEGER                 :: lines_number, atoms_number, steps_number, lines_count
    INTEGER                 :: i, j, k, l, n, io
    INTEGER                 :: tmp_int
    REAL(KIND=8)            :: tmp_real(2)

    OPEN(10,file=conq_in,position='rewind')

    !==============================================================>>
    ! Get fractional coordinates for 'Conquest_coord' file
    DO WHILE( sent_conq(1:25) /= sent_general_fracoor .and. sent_conq(1:25) /= sent_general_fracoor_ )
       !READ(10,'(A25)') sent_conq
       READ(10,'(A)') sent_conq
       sent_conq = TRIM(sent_conq)
       !write(*,*) sent_conq
    END DO

    BACKSPACE(10)

    open(1005, file="tmp.txt") ! Temporary file to convert integer into strings.
    write(1005,*) INDEX(sent_conq,"T")-2
    close(1005)

    open(1005, file="tmp.txt")
    read(1005,*) str
    !write(*,*) str
    close(1005)
    !write(*,'(A'//trim(str)//',X,A)') sent_conq,"-"
    READ(10,'(A'//trim(str)//',X,L)') sent_conq, frac_coord
    !READ(10,'(A32,X,L)') sent_conq, frac_coord
    !==============================================================<<

    REWIND(10)

    !==============================================================>>
    ! Get the number of species from the 'Conquest_input' file
    DO WHILE(sent_conq(1:23) /= sent_general_species .and. sent_conq(1:23) /= sent_general_species_ )
       !READ(10,'(A23)') sent_conq
       READ(10,'(A)') sent_conq
       sent_conq = TRIM(sent_conq)
       !write(*,*) sent_conq
    END DO

    BACKSPACE(10)
    
    open(1005, file="tmp2.txt") ! Temporary file to convert integer into strings.
    write(1005,*) len(TRIM(sent_conq))-2
    close(1005)

    open(1005, file="tmp2.txt")
    read(1005,*) str
    close(1005)
    READ(10,'(A'//trim(str)//',X,I4)') sent_conq, species_number
    !READ(10,'(A23,X,I4)') sent_conq, species_number
    ALLOCATE(species_ind(species_number))
    ALLOCATE(species_mass(species_number))
    ALLOCATE(species_name(species_number))
    !==============================================================<<

    REWIND(10)
    
    !==============================================================>>
    ! Get the number of species from the 'Conquest_input' file
    DO WHILE(sent_conq /= sent_species)
       READ(10,'(A64)') sent_conq
       sent_conq = TRIM(sent_conq)
    END DO

    WRITE(*,*) 'species_number:', species_number
    n = 1
    DO WHILE(n /= species_number + 1)
       READ(10,*) species_ind(n), species_mass(n), species_name(n)
       WRITE(*,*) species_ind(n), species_mass(n), species_name(n)
       n = n + 1
    END DO
    !==============================================================<<

    CLOSE(10)

    OPEN(10,file=conq_updatom_dat,position='rewind')  

    !==============================================================>>
    ! Get the number of lines from the 'UpdatedAtoms.dat' file
    n = 0
    DO
       READ(10, *, iostat=io) 
       IF( io < 0 ) THEN
          EXIT
       ELSE   
          n = n + 1
       END IF
    END DO
    lines_number = n - 1
    !==============================================================<<
    CLOSE(10)

    OPEN(10,file=conq_updatom_dat,position='rewind')    

    !==============================================================>>
    ! Get the number of atoms from the 'UpdatedAtoms.dat' file
    DO i = 1, 4
       READ(10,*)
    END DO
    READ(10,'(I12)') atoms_number
    !==============================================================<<

    CLOSE(10)

    !==============================================================>>
    !Calualte the number of points
    steps_number = lines_number/(5 + atoms_number) + 1
    lines_count  = 4 + atoms_number
    !==============================================================<<

    ALLOCATE(abc(3,3))
    ALLOCATE(ijk(atoms_number,3))
    ALLOCATE(xyz(atoms_number,3))
    ALLOCATE(ind(atoms_number))
    !ALLOCATE(for(atoms_number,3))
    !ALLOCATE(maxf(steps_number))
    !ALLOCATE(rmsf(steps_number))
    !ALLOCATE(resf(steps_number))
    !ALLOCATE(nrg(steps_number))

    ALLOCATE(ityp_cqcoord( atoms_number   ))
    ALLOCATE(atm_cqinp   ( species_number ))
    ALLOCATE(label_cqinp ( atoms_number   ))
 

    OPEN(10,file=conq_updatom_dat,position='rewind')   
    !OPEN(20,file=conq_updatom_for,position='rewind') 
    OPEN(30,file=conq_updatom_xyz) ! xyz          file
    !OPEN(40,file=conq_updatom_mld) ! molden (mld) file

    !W=mld=========================================>>
    !WRITE(40,*) sent_mld_fmt
    !WRITE(40,*) sent_mld_geoxyz
    !W=mld=========================================<<

    WRITE(*,*) 'atoms_number:', atoms_number
    !WRITE(*,*) 'steps_number:', steps_number

    DO k = 1, steps_number 

       READ(10,*) !tmp_int, nrg(k)
       DO i = 1, 3
          READ(10,*) (abc(i,j), j=1,3)
       END DO
       !abc = abc*a_0
       READ(10,*)

       !W=mld=xyz========================================>>
       WRITE(30,'(I12)') atoms_number
       WRITE(30,*) sent_mld_point, k
       !WRITE(40,'(I12)') atoms_number
       !WRITE(40,*) sent_mld_point, k
       !W=mld=xyz========================================<<

       DO i = 1, atoms_number
          READ(10,*) (ijk(i,j), j=1,3), ind(i)
       END DO

       IF (MINVAL(ijk) < 0.0d0) THEN
          ijk = ijk + 0.5d0
       END IF

       DO i = 1, atoms_number

          DO j = 1, 3
             IF (ijk(i,j) > 1.0d0) THEN     
                ijk(i,j) = -1.0d0 + ijk(i,j)
             ELSE 
                IF (ijk(i,j) < -1.0d0) THEN     
                   ijk(i,j) =  1.0d0 + ijk(i,j)
                END IF
             END IF
          END DO

          IF ( frac_coord ) THEN
             xyz(i,:) = MATMUL(abc,ijk(i,:))

          ELSE
             xyz(i,:) = ijk(i,:)
             
          END IF
          
          !xyz(i,:) = MATMUL(abc,ijk(i,:))
          !xyz = ijk

          l = 1
          DO WHILE (ind(i) /= species_ind(l))
             l = l + 1
          END DO
          !W=mld=xyz========================================>>
          WRITE(30,'(A2,3F12.6)') species_name(l), (xyz(i,j)*a_0, j=1, 3)
          !WRITE(40,'(A2,3F12.6)') species_name(l), (xyz(i,j), j=1, 3)
          !W=mld=xyz========================================<<

          label_cqinp(i) = trim(species_name(l))
          
       END DO

    END DO

    cel_cqinp = abc
    
    
    !W==========================================>>
    !WRITE(40,*) sent_mld_forces
    !W==========================================<<

    !DO k = 1, steps_number 

       !READ(20,*)
       !DO i = 1, 3
       !   READ(20,*)
       !END DO
       !READ(20,*)

       !DO i = 1, atoms_number
       !   READ(20,*) (for(i,j), j=1,3)
       !END DO

       !rmsf(k) = SQRT(SUM(for**2))/atoms_number
       !maxf(k) = MAXVAL(ABS(for))

       !W==========================================>>
       !WRITE(40,*) sent_mld_point, k
       !WRITE(40,'(I12)') atoms_number
       !W==========================================<<    

       !DO i = 1, atoms_number
       !   WRITE(40,'(3F12.6)') (for(i,j), j=1, 3)
       !END DO
    !END DO

    !W==========================================>>
    !WRITE(40,*) sent_mld_geoconv
    !WRITE(40,*) sent_mld_energy
    !W==========================================<<
    !DO k = 1, steps_number 
       !WRITE(40,'(F14.8)') nrg(k)
    !END DO

    !W==========================================>>
    !WRITE(40,*) sent_mld_maxf
    !W==========================================<<
    !DO k = 1, steps_number 
    !   WRITE(40,'(F14.8)') maxf(k)
    !END DO

    !W==========================================>>
    !WRITE(40,*) sent_mld_rmsf
    !W==========================================<<
    !DO k = 1, steps_number 
    !   WRITE(40,'(F14.8)') rmsf(k)
    !END DO

    CLOSE(10)
    CLOSE(20)
    CLOSE(30)
    !CLOSE(40)

    nat_cqinp    = atoms_number
    nsp_cqinp    = species_number
    ityp_cqcoord = ind
    atm_cqinp    = species_name

    
    !print*, species_name
    !print*,
    !print*, label_cqinp

    
    DEALLOCATE(species_ind)
    DEALLOCATE(species_mass)
    DEALLOCATE(species_name)  

    DEALLOCATE(abc)
    DEALLOCATE(ijk)
    DEALLOCATE(xyz)
    !DEALLOCATE(for)
    !DEALLOCATE(nrg)
    !DEALLOCATE(maxf)
    !DEALLOCATE(rmsf)
    !DEALLOCATE(resf)
    DEALLOCATE(ind)

  END SUBROUTINE dynpro_read_cqinp
  !
  !
  !
  SUBROUTINE dynpro_read_cqCoordForce( )

    USE dynpro_inp, ONLY : nat_cqinp, nsp_cqinp, atm_cqinp, ityp_cqcoord, label_cqinp, filecqd, evp_cqinp, cel_cqinp
    USE dynpro_inp, ONLY : filepos, filepos, filecel, filevel, fileevp, filefor, nat

    
    IMPLICIT NONE

    REAL(KIND=8)  :: xyz_cq(3)
    REAL(KIND=8)  :: for_cq(3)
    REAL(KIND=8)  :: vel_cq(3)
    REAL(KIND=8)  :: cell_cq(3,3)
       
    REAL(KIND=8)  :: energy_    
    INTEGER       :: lines_number, atoms_number, steps_number, lines_count
    INTEGER       :: i, j, k, l, n, io, index_, index__, k_print

    OPEN(10,file=filecqd,position='rewind')  

    !==============================================================>>
    ! Get the number of lines from the 'CoordFile.dat' file
    n = 0
    DO
       READ(10, *, iostat=io) 
       IF( io < 0 ) THEN
          EXIT
       ELSE   
          n = n + 1
       END IF
    END DO

    lines_number = n 
    print*, 'nb. lines:', lines_number
    !==============================================================<<

    CLOSE(10)

    !==============================================================>>
    !Calualte the number of points
    atoms_number = nat_cqinp
    steps_number = lines_number/(atoms_number + 1)
    print*, 'nb. steps:', steps_number
    !==============================================================<<

    WRITE(*,*) 'atoms_number:', atoms_number
    WRITE(*,*) 'steps_number:', steps_number
 
    OPEN(10,file=filecqd)!,position='unknown')   
    OPEN(20,file=filepos)!,position='unknown')   
    OPEN(30,file=filefor)!,position='unknown')   
    OPEN(40,file=filevel)!,position='unknown')   
    OPEN(50,file=filecel)!,position='unknown')   
    OPEN(60,file=fileevp)!,position='unknown')   

    evp_cqinp = 0.0d0

    !k_print = 0 
    DO k = 1, steps_number 
       
       read(10,*) index_, energy_
       !print*, index_, energy_
       evp_cqinp(4) = energy_
       
       index_ = k      
 
       write(20,'(I12,X,F20.8)') index_, energy_
       write(30,'(I12,X,F20.8)') index_, energy_
       write(40,'(I12,X,F20.8)') index_, energy_
       write(50,'(I12,X,F20.8)') index_, energy_

       write(60,'(I12,10F20.8)') index_, (evp_cqinp(j), j=1,10)

       do n = 1, 3
          write(50,'(3F12.8)') (cel_cqinp(n,j), j=1,3)
       end do
       
       DO i = 1, atoms_number
          
          read(10,*) index__ , (xyz_cq(j), j=1,3), (for_cq(j), j=1,3), (vel_cq(j), j=1,3)

          write(20,'(3F20.8)') (xyz_cq(j), j=1,3)
          write(30,'(3F20.8)') (for_cq(j), j=1,3)
          write(40,'(3F20.8)') (vel_cq(j), j=1,3)
          
       END DO

      !k_print = k_print + nprint
       
    END DO
    
    CLOSE(10)
    CLOSE(20)
    CLOSE(30)
    CLOSE(40)
    CLOSE(50)
    CLOSE(60)

    RETURN
  END SUBROUTINE dynpro_read_cqCoordForce


  SUBROUTINE dynpro_read_cqmdframes( )

    USE dynpro_inp, ONLY : nat_cqinp, nsp_cqinp, atm_cqinp, ityp_cqcoord, label_cqinp, filecqd, evp_cqinp, cel_cqinp
    USE dynpro_inp, ONLY : filepos, filepos, filecel, filevel, fileevp, filefor, nat

    
    IMPLICIT NONE

    REAL(KIND=8)  :: xyz_cq(3)
    REAL(KIND=8)  :: for_cq(3)
    REAL(KIND=8)  :: vel_cq(3)
    REAL(KIND=8)  :: cell_cq(3,3)
    REAL(KIND=8)  :: stress_cq(3,3)
       
    REAL(KIND=8)  :: energy_    
    INTEGER       :: lines_number, atoms_number, steps_number, lines_count
    INTEGER       :: i, j, k, l, n, io, index_, index__, k_print
    CHARACTER(len=5)  :: tmp 

    INTEGER :: tmp1, tmp2
    
    
    OPEN(10,file=filecqd,position='rewind')  

    write(*,*) 'extract MD data from file ', filecqd
    
    !==============================================================>>
    ! Get the number of lines from the 'md.frames' file
    n = 0
    DO
       READ(10, *, iostat=io) 
       IF( io < 0 ) THEN
          EXIT
       ELSE   
          n = n + 1
       END IF
    END DO

    lines_number = n 
    write(*,*) 'nb. lines:', lines_number
    !==============================================================<<

    CLOSE(10)

    !==============================================================>>
    !Calulate the number of points
    atoms_number = nat_cqinp
    steps_number = lines_number/(atoms_number*3 + 18)
    print*, 'nb. steps:', steps_number
    !==============================================================<<

    WRITE(*,*) 'atoms_number:', atoms_number
    WRITE(*,*) 'steps_number:', steps_number
 
    OPEN(10,file=filecqd)!,position='unknown')   
    OPEN(20,file=filepos)!,position='unknown')   
    OPEN(30,file=filefor)!,position='unknown')   
    OPEN(40,file=filevel)!,position='unknown')   
    OPEN(50,file=filecel)!,position='unknown')   
    OPEN(60,file=fileevp)!,position='unknown')   

    evp_cqinp = 0.0d0

    !k_print = 0 
    DO k = 1, steps_number
       
       read(10,*) tmp, index_
       !print*, index_, energy_
       evp_cqinp(4) = 0.0d0

       index_ = k      
 
       write(20,'(I12,X,F20.8)') index_, energy_
       write(30,'(I12,X,F20.8)') index_, energy_
       write(40,'(I12,X,F20.8)') index_, energy_
       write(50,'(I12,X,F20.8)') index_, energy_
       write(60,'(I12,10F20.8)') index_, (evp_cqinp(j), j=1,10)

       do n = 1, 3
          write(50,'(3F12.8)') (cel_cqinp(n,j), j=1,3)
       end do

       read(10,*)

       read(10,*) (cell_cq(j,1), j=1,3) 
       read(10,*) (cell_cq(j,2), j=1,3)
       read(10,*) (cell_cq(j,3), j=1,3)

       !write(*,*) (cell_cq(j,1), j=1,3) 
       !write(*,*) (cell_cq(j,2), j=1,3)
       !write(*,*) (cell_cq(j,3), j=1,3)

       
       read(10,*)
       read(10,*)

       read(10,*) (stress_cq(j,1), j=1,3)
       read(10,*) (stress_cq(j,2), j=1,3)
       read(10,*) (stress_cq(j,3), j=1,3)

       read(10,*)
       read(10,*)

       DO i = 1, atoms_number          
          read(10,*) tmp1, tmp2, (xyz_cq(j), j=1,3)          
          write(20,'(3F20.8)') (xyz_cq(j), j=1,3)          
       END DO
       
       read(10,*)
       read(10,*)
       
       DO i = 1, atoms_number          
          read(10,*) tmp1, tmp2, (vel_cq(j), j=1,3)          
          write(30,'(3F20.8)') (vel_cq(j), j=1,3)          
       END DO
       
       read(10,*)
       read(10,*)
       
       DO i = 1, atoms_number          
          read(10,*) tmp1, tmp2, (for_cq(j), j=1,3)          
          write(40,'(3F20.8)') (for_cq(j), j=1,3)          
       END DO

       read(10,*)
       read(10,*)
      !k_print = k_print + nprint
       
    END DO
    
    CLOSE(10)
    CLOSE(20)
    CLOSE(30)
    CLOSE(40)
    CLOSE(50)
    CLOSE(60)

    RETURN
  END SUBROUTINE dynpro_read_cqmdframes

  
END MODULE dynpro_cqc
