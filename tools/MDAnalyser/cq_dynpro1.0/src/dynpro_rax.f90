MODULE dynpro_rax

  USE kinds,      ONLY : DP
  
  IMPLICIT NONE

CONTAINS
  
  SUBROUTINE tensor()

    USE dynpro_inp,  ONLY : fileeig, filevec, filetens, fileran
    USE dynpro_inp,  ONLY : fileR, fileRiso, fileRant, fileRsym,  fileRsph
    USE dynpro_inp,  ONLY : lrandom, verbose_nmr
    USE dynpro_inp,  ONLY : iteration, time_dt
    USE dynpro_ran,  ONLY : random

    IMPLICIT NONE    

    CHARACTER(len=256)         :: fileRsphIm, fileRsphRe
    CHARACTER(len=256)         :: fileR_kq(0:2,0:2)

    REAL(DP),    ALLOCATABLE   :: U(:,:), U_t(:,:)
    REAL(DP),    ALLOCATABLE   :: E(:,:)
    REAL(DP),    ALLOCATABLE   :: eig(:)
    REAL(DP),    ALLOCATABLE   :: tmp(:,:)
    REAL(DP),    ALLOCATABLE   :: Unity(:,:)

    REAL(DP),    ALLOCATABLE   :: R(:,:)
    REAL(DP),    ALLOCATABLE   :: Riso(:,:)
    REAL(DP),    ALLOCATABLE   :: Rant(:,:)
    REAL(DP),    ALLOCATABLE   :: Rsym(:,:)

    COMPLEX(DP)                :: R_kq(0:2,0:2)
    COMPLEX                    :: R_kq_coef(0:2,0:2)
    CHARACTER(len=10)          :: R_kq_name(0:2,0:2)

    REAL(DP),    ALLOCATABLE   :: work(:,:)
    COMPLEX(DP)                :: workc(0:2,0:2)

    REAL(DP)                  :: Lambda, Delta, eta

    INTEGER    ::  i, j, k, l, n, m, q

    n = iteration
    m = 3

    IF (verbose_nmr < 2) THEN
       !CALL errore(' Dynpro/acf ',' verbose_nmr >= 2 ', verbose_nmr)
    END IF

    ALLOCATE(E(m,m), U(m,m), U_t(m,m), eig(m), Unity(m,m))
    ALLOCATE(R(m,m), Riso(m,m), Rant(m,m), Rsym(m,m)) 
    ALLOCATE(work(m,m))

    E     = 0.0d0 ! Eigenvalues  matrix 
    eig   = 0.0d0 ! Eigenvalues  vector
    U     = 0.0d0 ! Eigenvectors matrix
    U_t   = 0.0d0 ! ... transposed
    work  = 0.0d0 ! Working matrix (real)
    Unity = 0.0d0 ! Unitary matrix

    R     = 0.0d0 ! Full tensor  (Cartesian frame)
    Riso  = 0.0d0 ! Isotropic     component
    Rant  = 0.0d0 ! Antisymmetric component
    Rsym  = 0.0d0 ! Symmetric     component

    R_kq  = 0.0d0 ! Matrix containing the irreducible (spherical) components

    fileRsphIm = TRIM(fileRsph) // 'Im'
    fileRsphRe = TRIM(fileRsph) // 'Re'

    DO i = 1, m
       Unity(i,i) = 1.0d0
    END DO

    OPEN(7,  file=filetens, status='replace', position='rewind')
    OPEN(8,  file=fileR,    status='replace', position='rewind')
    OPEN(10, file=fileRiso, status='replace', position='rewind')
    OPEN(11, file=fileRant, status='replace', position='rewind')
    OPEN(12, file=fileRsym, status='replace', position='rewind')
    OPEN(13, file=fileRsph, status='replace', position='rewind')
    OPEN(14, file=fileRsphRe, status='replace', position='rewind')
    OPEN(15, file=fileRsphIm, status='replace', position='rewind')

    CALL tensor_char(0,fileR_kq,R_kq_coef,R_kq_name)

    l = 0
    DO k = 0, m - 1
       DO q = 0 , m - 1
          OPEN(30+l,file=fileR_kq(k,q), status='replace', position='rewind')
          l = l + 1
       END DO
    END DO

    IF (lrandom) THEN
       CALL random()
       OPEN(4, file=fileran, status='old', position='rewind')
    ELSE
       OPEN(5, file=filevec, status='old', position='rewind')
       OPEN(6, file=fileeig, status='old', position='rewind')
    END IF


    DO i = 1, n

       !==---------------------------------------------------------------------==
       
       IF (lrandom) THEN          
          DO j = 1, m
             READ(4,'(3F14.6)') (R(j,k), k=1,m)
          END DO
       ELSE
          DO j = 1, m
             READ(5,*) (U(j,k), k = 1,m)
          END DO
          READ(6,*)  (E(j,j), j = 1,m)
          
          U_t   = TRANSPOSE(U)
          work  = MATMUL(U,E)
          R     = MATMUL(work,U_t)         

          U_t   = 0.0d0
          U     = 0.0d0
          work  = 0.0d0
       END IF

       !==---------------------------------------------------------------------==

       CALL cartesian(m,R,Riso,Rant,Rsym)
              
       work = Riso + Rsym  

       !CALL CDIAGH(3, work, 3, eig, U)
       eig = 0d0
       
       DO j = 1, m
          E(j,j) = eig(j)
       END DO

       CALL spherical(m,R,R_kq)

       !==---------------------------------------------------------------------==
          
       WRITE(7,'(7X,A6,37X,A4,39X,A4)') 'Tensor', 'Riso', 'Rant'
       DO j = 1, m         
          WRITE(7,'(3F14.6,X,3F14.6,X,3F14.6)') &
               (R(j,k), k=1,m) , (Riso(j,k), k=1,m) , (Rant(j,k), k=1,m)
       END DO
       
       WRITE(7,'(7X,A4,39X,A11,32X,A12)') 'Rsym', 'Eigenvalues', 'Eigenvectors'
       DO j =1, m
          WRITE(7,'(3F14.6,X,3F14.6,X,3F14.6)') &
               (Rsym(j,k), k=1,m), (E(j,k), k=1,m), (U(j,k), k=1,m)
       END DO
       WRITE(7,'(3(A43))') ('<<=======================================>>', k=1,3)
                 
       work = Riso + Rant + Rsym

       DO j = 1, m
          WRITE(8, '(3F14.6)') (work(j,k), k=1,m)
          WRITE(10,'(3F14.6)') (Riso(j,k), k=1,m)
          WRITE(11,'(3F14.6)') (Rant(j,k), k=1,m)
          WRITE(12,'(3F14.6)') (Rsym(j,k), k=1,m)
       END DO
       DO k = 0, m - 1
          WRITE(13,'(3F14.6,3F14.6)') (REAL(R_kq(k,q)), q = 0, m-1), (AIMAG(R_kq(k,q)), q = 0, m-1) 
          WRITE(14,'(3F14.6)') ( REAL(R_kq(k,q)), q = 0, m-1)
          WRITE(15,'(3F14.6)') (AIMAG(R_kq(k,q)), q = 0, m-1)
       END DO
 
       l = 0
       DO k = 0, m - 1
          DO q = 0, m - 1
             WRITE(30+l,'(F14.6,F14.6)') REAL(R_kq(k,q)), AIMAG(R_kq(k,q))
             l = l + 1
          END DO
       END DO

    END DO

    CALL tensor_char(1,fileR_kq,R_kq_coef,R_kq_name)

    IF (lrandom) THEN
       CLOSE(4)
    ELSE
       CLOSE(5)
       CLOSE(6)
    END IF

    DO i = 7,15
       CLOSE(i)
    END DO


    DEALLOCATE(E, U, U_t)
    DEALLOCATE(R, Riso, Rant, Rsym)
    DEALLOCATE(Unity,eig)
    DEALLOCATE(work)
    
  END SUBROUTINE tensor
  
  SUBROUTINE spherical(m,R,R_kq)
    
    IMPLICIT NONE    

    INTEGER,      INTENT(in)   :: m
    REAL(DP),     INTENT(in)   :: R(m,m)
    COMPLEX(DP),  INTENT(out)  :: R_kq(0:2,0:2)

    COMPLEX(DP)                :: R_0
    COMPLEX(DP)                :: R_1(-1:1)
    COMPLEX(DP)                :: R_2(-2:2)
    REAL(DP)                   :: Lambda

    REAL(DP)   :: r2(-1:1) 
    REAL(DP)   :: r3(-1:1) 
    REAL(DP)   :: r6(-1:1) 
    REAL(DP)   :: s2(-1:1) 

    INTEGER               :: i, j, k

     r2(-1:1) = (/ -1.0d0/sqrt(2.0d0), 0.0d0, 1.0d0/sqrt(2.0d0) /)
     r3(-1:1) = (/ -1.0d0/sqrt(3.0d0), 0.0d0, 1.0d0/sqrt(3.0d0) /)
     r6(-1:1) = (/ -1.0d0/sqrt(6.0d0), 0.0d0, 1.0d0/sqrt(6.0d0) /)
     s2(-1:1) = (/ -1.0d0/2.0d0      , 0.0d0, 1.0d0/2.0d0       /)



    Lambda = SUM(R, mask=reshape(source = (/ ((j==k,j=1,m),k=1,m) /), &
         shape = shape(R)))

       R_0     = r3(-1)*CMPLX( Lambda         ,   0.0d0           )

       R_1( 0) = r2(-1)*CMPLX( 0.0d0          ,   R(1,2) - R(2,1) )
       R_1( 1) = s2(-1)*CMPLX( R(3,1) - R(1,3),   R(3,2) - R(2,3) )
       R_1(-1) = s2(-1)*CMPLX( R(3,1) - R(1,3), - R(3,2) + R(2,3) )
       
       R_2( 0) = r6( 1)*CMPLX(3*R(3,3)- Lambda,   0.0d0           )
       R_2( 1) = s2(-1)*CMPLX( R(1,3) + R(3,1),   R(2,3) + R(3,2) )
       R_2(-1) = s2( 1)*CMPLX( R(1,3) + R(3,1), - R(2,3) - R(3,2) )
       R_2( 2) = s2( 1)*CMPLX( R(1,1) - R(2,2),   R(1,2) + R(2,1) )
       R_2(-2) = s2( 1)*CMPLX( R(1,1) - R(2,2), - R(1,2) - R(2,1) )

       R_kq(0, 0) = R_0       ;  R_kq(0, 1) = R_1( 1) ;  R_kq(0, 2) = R_2( 1)
       R_kq(1, 0) = R_1(-1)   ;  R_kq(1, 1) = R_1( 0) ;  R_kq(1, 2) = R_2( 2) 
       R_kq(2, 0) = R_2(-1)   ;  R_kq(2, 1) = R_2(-2) ;  R_kq(2, 2) = R_2( 0)

  END SUBROUTINE spherical

  SUBROUTINE tensor_char(in,fileR_kq,R_kq_coef,R_kq_name)

    USE dynpro_inp,  ONLY : fileR

    IMPLICIT NONE
    
    INTEGER,            INTENT(in)  :: in
    CHARACTER(len=256), INTENT(out) :: fileR_kq(0:2,0:2)
    CHARACTER(len=10),  INTENT(out) :: R_kq_name(0:2,0:2)
    COMPLEX,            INTENT(out) :: R_kq_coef(0:2,0:2)
    CHARACTER(len=3)                :: char_Re, char_Im

    INTEGER                         :: tmp_Re, tmp_Im
    INTEGER                         :: i, j, k, q, l, m

    SELECT CASE (in)
       
    CASE(0)
       
       m = 3
       
       R_kq_coef(0, 0) = CMPLX(0,0)  ; R_kq_coef(0, 1) = CMPLX(1,1)  ; R_kq_coef(0, 2) = CMPLX(2,1)
       R_kq_coef(1, 0) = CMPLX(1,-1) ; R_kq_coef(1, 1) = CMPLX(1,0)  ; R_kq_coef(1, 2) = CMPLX(2,2)
       R_kq_coef(2, 0) = CMPLX(2,-1) ; R_kq_coef(2, 1) = CMPLX(2,-2) ; R_kq_coef(2, 2) = CMPLX(2,0)
       
       DO k = 0, m - 1
          DO q = 0, m -1
             tmp_Re  = REAL(R_kq_coef(k,q)) ; tmp_Im = AIMAG(R_kq_coef(k,q))
             char_Re = CHAR(48+tmp_Re)
             char_Im = CHAR(48+ABS(tmp_Im))
             IF (tmp_Im < 0) THEN
                char_Im = '-' // TRIM(char_Im)
             ELSE
                char_Im = '+' // TRIM(char_Im)
             END IF
             fileR_kq(k,q) = TRIM(fileR) // TRIM(char_Re) // &
                  TRIM(char_Im)
             R_kq_name(k,q) = 'R' // '^{' // TRIM(char_Re)// '}' // '_{' // TRIM(char_Im) // '}'
          END DO
       END DO
    
    CASE(1)

       DO i = 30, 39
          CLOSE(i)
       END DO

    END SELECT

  END SUBROUTINE tensor_char

  SUBROUTINE cartesian(m,R,Riso,Rant,Rsym)
    
    IMPLICIT NONE    

    INTEGER,      INTENT(in)   :: m
    REAL(DP),     INTENT(in)   :: R(m,m)
    REAL(DP),     INTENT(out)  :: Riso(m,m)
    REAL(DP),     INTENT(out)  :: Rant(m,m)
    REAL(DP),     INTENT(out)  :: Rsym(m,m)

    REAL(DP)                   :: Unity(m,m)
    REAL(DP)                   :: Lambda

    INTEGER                    :: i, j, k

    Lambda = SUM(R, mask=reshape(source = (/ ((j==k,j=1,m),k=1,m) /), &
         shape = shape(R)))

    Unity= 0.0d0
    DO i = 1, m
       Unity(i,i) = 1.0d0
    END DO
  
    Riso   = (1.d0/3.d0)*(Lambda*Unity)
    Rant   = (1.d0/2.d0)*(R - TRANSPOSE(R))
    Rsym   = (1.d0/2.d0)*(R + TRANSPOSE(R)) - Riso 


  END SUBROUTINE cartesian

  SUBROUTINE tensor_stat(file,n,m,name,V,type)

    USE dynpro_mod,  ONLY : moment

    IMPLICIT NONE

    CHARACTER(len=256),  INTENT(in)    :: file
    CHARACTER(len=7),    INTENT(in)    :: name
    INTEGER,             INTENT(in)    :: type
    INTEGER,             INTENT(in)    :: n, m
    COMPLEX(DP),         INTENT(out)   :: V(m,m)
    
    CHARACTER(len=256)                 :: file_Im, file_Re


    REAL(DP)                           :: Re_U(m,m), Im_U(m,m)     
    REAL(DP)                           :: Re_R(n,m,m), Im_R(n,m,m)
    REAL(DP)                           :: Re_Uadev(m,m), Re_Usdev(m,m), Re_Uvar(m,m)
    REAL(DP)                           :: Im_Uadev(m,m), Im_Usdev(m,m), Im_Uvar(m,m)
    REAL(DP)                           :: work(n)

    INTEGER                            :: i, j, k, unit

    Re_R     = 0.0D0 ;  Im_R     = 0.0D0
    Re_U     = 0.0D0 ;  Im_U     = 0.0D0
    Re_Uadev = 0.0D0 ;  Im_Uadev = 0.0D0
    Re_Usdev = 0.0D0 ;  Im_Usdev = 0.0D0
    Re_Uvar  = 0.0D0 ;  Im_Uvar  = 0.0D0

    V     = 0.0D0
    work  = 0.0d0

    file_Im = TRIM(file) // 'Im'
    file_Re = TRIM(file) // 'Re'
    
    SELECT CASE(type)

    CASE (0)
       
       OPEN(10, file=file, status='old', position='rewind')
       DO i = 1, n
          DO j = 1, m
             READ(10,'(3F14.6)') (Re_R(i,j,k), k = 1,m)
          END DO
       END DO
       CLOSE(10)

       DO k = 1, m
          DO j = 1, m
             CALL moment(Re_R(:,j,k),n,Re_U(j,k),Re_Uadev(j,k),Re_Usdev(j,k),Re_Uvar(j,k))
          END DO
       END DO
       
       V = CMPLX( Re_U, 0.0d0 )

       WRITE(*,*) '======= Tensor Statistic  ================================>>'
       WRITE(*,'(2X,7A)') name
       DO j = 1, m
          WRITE(*,'(3F14.6)') (Re_U(j,k), k = 1,m)
       END DO
       WRITE(*,'(2X,8A)') '+/- sdev'
       DO j = 1, m
          WRITE(*,'(3F14.6)') (Re_Usdev(j,k), k = 1,m)
       END DO
       WRITE(*,*) '==========================================================<<'
       
    CASE (1)

       OPEN(10, file=file_Re, status='old', position='rewind')
       OPEN(11, file=file_Im, status='old', position='rewind')
       DO i = 1, n
          DO j = 1, m
             READ(10,'(3F14.6)') (Re_R(i,j,k), k = 1, m)
             READ(11,'(3F14.6)') (Im_R(i,j,k), k = 1, m)
          END DO
       END DO
       CLOSE(10)
       CLOSE(11)

       DO k = 1, m
          DO j = 1, m
             CALL moment(Re_R(:,j,k),n,Re_U(j,k),Re_Uadev(j,k), &
                  Re_Usdev(j,k),Re_Uvar(j,k))
             CALL moment(Im_R(:,j,k),n,Im_U(j,k),Im_Uadev(j,k), &
                  Im_Usdev(j,k),Im_Uvar(j,k))
          END DO
       END DO

       V = CMPLX( Re_U, Im_U )

       WRITE(*,*) '======= Tensor Statistic  ================================>>'
       WRITE(*,'(2X,7A,7A)') name, '_{Real}'
       DO j = 1, m
          WRITE(*,'(3F14.6)') (Re_U(j,k), k = 1,m)
       END DO
       WRITE(*,'(2X,8A)') '+/- sdev'
       DO j = 1, m
          WRITE(*,'(3F14.6)') (Re_Usdev(j,k), k = 1,m)
       END DO
       WRITE(*,'(2X,7A,5A)') name, '_{Im}'
       DO j = 1, m
          WRITE(*,'(3F14.6)') (Im_U(j,k), k = 1,m)
       END DO
       WRITE(*,'(2X,8A)') '+/- sdev'
       DO j = 1, m
          WRITE(*,'(3F14.6)') (Im_Usdev(j,k), k = 1,m)
       END DO

       WRITE(*,*) '==========================================================<<'
       
    CASE DEFAULT   
       !CALL errore(' - Dynpro/tensor_stat - ',' no default here ',type)
    END SELECT


  END SUBROUTINE tensor_stat


  SUBROUTINE tensor_acf()

    USE dynpro_inp,  ONLY : iteration, time_dt
    USE dynpro_inp,  ONLY : fileR, fileRiso, fileRant, fileRsym
    USE dynpro_inp,  ONLY : fileRsph
    USE dynpro_inp,  ONLY : lstat, verbose_nmr, cor_function

    IMPLICIT NONE

    CHARACTER(len=256)                 :: fileacf_in(16)
    CHARACTER(len=256)                 :: fileacf_out(16)
    CHARACTER(len=7)                   :: R_name(16)
    CHARACTER(len=10)                  :: R_kq_name(0:2,0:2)
    CHARACTER(len=256)                 :: fileR_kq(0:2,0:2)

    REAL(DP),            ALLOCATABLE   :: ReR_t0(:,:,:)
    REAL(DP),            ALLOCATABLE   :: ImR_t0(:,:,:)
    
    COMPLEX(DP),         ALLOCATABLE   :: R_t0(:,:,:)
    COMPLEX(DP),         ALLOCATABLE   :: R_t(:,:,:)
    COMPLEX(DP),         ALLOCATABLE   :: C(:), U(:,:)

    COMPLEX                            :: R_kq_coef(0:2,0:2)
    COMPLEX(DP)                        :: Norm(16), tmp_Norm
    REAL(DP)                           :: tmp_Re, tmp_Im

    INTEGER    :: type
    INTEGER    :: time, t0, t
    INTEGER    :: i, ii, j, k, l, n, m, q, frame

    !==---------------------------------------------------------------------==
    IF (verbose_nmr >= 10) THEN
       verbose_nmr = 16
       CALL tensor_char(0,fileR_kq,R_kq_coef,R_kq_name)
    END IF
    !==---------------------------------------------------------------------==
    n     = iteration
    m     = 3
    ALLOCATE(R_t0(n,m,m), R_t(n,m,m), C(n))
    ALLOCATE(ReR_t0(n,m,m), ImR_t0(n,m,m))
    ALLOCATE(U(m,m))

    R_t0 = 0.0d0
    R_t  = 0.0d0
    Norm = 0.0d0
    C    = 0.0d0
    U    = 0.0d0
    !==---------------------------------------------------------------------==
    R_name(1) = 'R^{car}'
    R_name(2) = 'R^{sph}'
    R_name(3) = 'R^{iso}'
    R_name(4) = 'R^{ant}'
    R_name(5) = 'R^{sym}'
    R_name(6) = 'R^{Re}'
    R_name(7) = 'R^{Im}'

    fileacf_in(1) = TRIM(fileR)
    fileacf_in(2) = TRIM(fileRsph)
    fileacf_in(3) = TRIM(fileRiso)
    fileacf_in(4) = TRIM(fileRant)
    fileacf_in(5) = TRIM(fileRsym)
    fileacf_in(6) = TRIM(fileRsph) // 'Re'
    fileacf_in(7) = TRIM(fileRsph) // 'Im'

    l = 0
    DO k = 0, m - 1 
       DO q = 0, m - 1
          fileacf_in (8+l) = TRIM(fileR_kq(k,q))
          fileacf_out(8+l) = TRIM(fileR_kq(k,q)) // '.acf'    
          R_name(8+l)      = TRIM(R_kq_name(k,q))
          l = l + 1
       END DO
    END DO

    fileacf_out(1) = TRIM(fileR)    // '.acf'
    fileacf_out(2) = TRIM(fileRsph) // '.acf'
    fileacf_out(3) = TRIM(fileRiso) // '.acf'
    fileacf_out(4) = TRIM(fileRant) // '.acf'
    fileacf_out(5) = TRIM(fileRsym) // '.acf'
    fileacf_out(6) = TRIM(fileacf_in(6)) // '.acf'
    fileacf_out(7) = TRIM(fileacf_in(7)) // '.acf'    
    !==---------------------------------------------------------------------== 
    
    DO ii = 1, verbose_nmr

       IF (lstat) THEN
          IF (ii == 2) THEN
             CALL tensor_stat(fileacf_in(ii),n,m,R_name(ii),U,1)
          ELSE IF (ii < 2 .or. ii < 8) THEN
             CALL tensor_stat(fileacf_in(ii),n,m,R_name(ii),U,0)
          END IF
          
          cor_function = TRIM(ADJUSTL(cor_function))

          SELECT CASE(cor_function)             
          CASE('default')
             U = 0.0d0
          CASE('acf')
             U = 0.0d0
          CASE('acv')             
          CASE DEFAULT   
             !CALL errore(' - Dynpro/tensor_acf - ',' unrecognized ', cor_function)
          END SELECT

       END IF
  
       R_t0     = 0.0d0
       R_t      = 0.0d0
       ReR_t0   = 0.0d0
       ImR_t0   = 0.0d0
       C        = 0.0d0
       tmp_Re   = 0.0d0
       tmp_Im   = 0.0d0

       SELECT CASE(ii)
          
       CASE DEFAULT
          
          OPEN(5, file=fileacf_in(ii), status='old', position='rewind')
          DO t0 = 1, n
             DO j = 1, m
                READ(5,'(3F14.6)') (ReR_t0(t0,j,k), k = 1,m)
             END DO
          END DO
          CLOSE(5)

       CASE(2)
          
          OPEN(5, file=fileacf_in(ii),   status='old', position='rewind')
          DO t0 = 1, n
             DO j = 1, m
                READ(5,'(3F14.6, 3F14.6)') (ReR_t0(t0,j,k), k=1,m), (ImR_t0(t0,j,k), k=1,m)
             END DO
          END DO
          CLOSE(5)
          
       CASE(6)
          
          OPEN(5, file=fileacf_in(ii), status='old', position='rewind')
          DO t0 = 1, n
             DO j = 1, m
                READ(5,'(3F14.6)') (ReR_t0(t0,j,k), k = 1,m)
             END DO
          END DO
          CLOSE(5)
          
       CASE(7)
          
          OPEN(5, file=fileacf_in(ii), status='old', position='rewind')
          DO t0 = 1, n
             DO j = 1, m
                READ(5,'(3F14.6)') (ImR_t0(t0,j,k), k = 1,m)
             END DO
          END DO
          CLOSE(5)

       CASE(8,9,10,11,12,13,14,15,16)

          OPEN(5, file=fileacf_in(ii), status='old', position='rewind')
          DO t0 = 1, n
             READ(5,'(3F14.6, 3F14.6)') tmp_Re, tmp_Im
             DO j = 1, m
                ReR_t0(t0,j,j) = tmp_Re
                ImR_t0(t0,j,j) = tmp_Im
             END DO
          END DO
          CLOSE(5)          
   
       END SELECT

       R_t0 = CMPLX( ReR_t0, ImR_t0 )

       CALL acf_order_n(ii,m,n,fileacf_out(ii), R_t0, R_t, R_kq_coef, U, tmp_Norm)

       Norm(ii) = tmp_Norm

       CALL tensor_psd(n,C,fileacf_out(ii), fileacf_out(ii)) 
       
    END DO
    !==---------------------------------------------------------------------==
    WRITE(*,*) '======= Tensor Norm  =====================================>>'
    DO ii = 1, 2
       WRITE(*,'(2X,A7,2X,2F20.6,A)') &
            R_name(ii), REAL(Norm(ii)), AIMAG(Norm(ii)),'i'
    END DO
    IF (verbose_nmr >= 5) THEN
       WRITE(*,*)
       DO ii = 3, 5
          WRITE(*,'(2X,A7,2X,2F20.6,A)') &
               R_name(ii), REAL(Norm(ii)),AIMAG(Norm(ii)),'i'
       ENDDO
       WRITE(*,'(2X,A7,2X,2F20.6,A)') 'Sum', &
            SUM(REAL(Norm(3:5))), SUM(AIMAG(Norm(2:4))),'i' 
    ELSE  IF (verbose_nmr > 5) THEN
       DO ii = 6, verbose_nmr
          WRITE(*,'(2X,A7,2X,2F20.6,A)') &
               R_name(ii), REAL(Norm(ii)),AIMAG(Norm(ii)),'i'
       END DO
    END IF
    WRITE(*,*) '==========================================================<<'
    !==---------------------------------------------------------------------==

    DEALLOCATE(R_t0,R_t,C)
    DEALLOCATE(ReR_t0,ImR_t0,U)

    
  END SUBROUTINE tensor_acf


  SUBROUTINE tensor_psd(n,C,filein,fileout)

    USE dynpro_mod,  ONLY : PSD
    USE dynpro_inp,  ONLY : iteration, time_step, iprint, interval
    USE dynpro_inp,  ONLY : syst_ID_lab, syst_ID_spe

    IMPLICIT NONE

    INTEGER,            INTENT(in)     :: n
    COMPLEX(DP),        INTENT(in)     :: C(n)
    CHARACTER(len=256), INTENT(in)     :: filein, fileout
    REAL(DP)                           :: C_Re(n)


    INTEGER :: m_cof = 512

    C_Re = REAL(C)    

    CALL PSD(1, C_Re, 1, .false., n, syst_ID_lab, syst_ID_spe, & 
         filein, .false., 1, 1, 30, m_cof, time_step, & 
         interval, iprint, 40, fileout, 1, 0.0d0)

  END SUBROUTINE tensor_psd


  SUBROUTINE acf_order_n(debug,m,n,file, R_t0, R_t, R_kq_coef, U, Norm)
    
    USE dynpro_inp,  ONLY : iteration
    USE dynpro_inp,  ONLY : cor_order

    IMPLICIT NONE
    
    INTEGER,            INTENT(in)     :: debug, m, n                      
    CHARACTER(len=256), INTENT(in)     :: file
    
    COMPLEX(DP),        INTENT(in)     :: R_t0(n,m,m)
    COMPLEX(DP),        INTENT(in)     :: R_t(n,m,m)
    COMPLEX(DP),        INTENT(in)     :: U(m,m)
    COMPLEX,            INTENT(in)     :: R_kq_coef(0:2,0:2)
    COMPLEX(DP),        INTENT(out)    :: Norm

    COMPLEX(DP)                        :: work_t0(n,m,m)
    COMPLEX(DP)                        :: C(n)

    INTEGER    :: time, t0, tau, t
    INTEGER    :: i, j, k, q, p, order

    CHARACTER(len=256)                 :: file_fft

    Norm    = 0.0d0
    C       = 0.0d0
    order   = cor_order
    work_t0 = 0.0d0

    file_fft= TRIM(file)//'1.dat'

    DO t0 = 1, n
       DO k = 1, m
          DO q = 1, m
             !IF (debug == 2) THEN
             !   Norm = Norm + REAL((-1)**(REAL(R_kq_coef(k-1,q-1)) + &
             !        AIMAG(R_kq_coef(k-1,q-1))))*(R_t0(t0,k,q) - U(k,q))   * &
             !        (R_t0(t0,q,k) - U(q,k))
             !ELSE
             Norm = Norm + (R_t0(t0,k,q) - U(k,q))*CONJG(R_t0(t0,k,q) - U(k,q))
             !work_t0(t0,k,q)    = -(R_t0(t0,k,q) - U(k,q))*CMPLX(0.0d0,1.0d0)
             !DO p = 1, order - 1
             !   work_t0(t0,k,q) = work_t0(t0,k,q)*CONJG(R_t0(t0,k,q) - U(k,q))*CMPLX(0.0d0,1.0d0)                
             !END DO
             !Norm               = Norm + work_t0(t0,k,q)
             !END IF                          
          END DO
       END DO
    END DO
    
    tau  = 0
    t    = 1
    t0   = 1   
      
    OPEN(6, file=file, status='replace', position='rewind')
    OPEN(7, file=file_fft, status='replace', position='rewind')

    time_loop: DO time = 1, n
       origin_loop: DO t0 = 1, n - tau
          
          t = t0 + tau
          
          DO k = 1, m
             DO q = 1, m                   
                !IF (debug == 2) THEN
                !   C(time) = C(time) + REAL((-1)**(REAL(R_kq_coef(k-1,q-1)) + &
                !        AIMAG(R_kq_coef(k-1,q-1))))*(R_t0(t0,k,q) - U(k,q)) * &
                !        (R_t0(t,q,k) - U(q,k))
                !ELSE
                C(time) = C(time) + (R_t0(t0,k,q) - U(k,q))*CONJG(R_t0(t,k,q) - U(k,q))
                !END IF
             END DO
          END DO
          
       END DO origin_loop
       
       tau       = tau + 1
       C(time)   = C(time)/Norm
       
       WRITE(6,*)  time-1,   REAL(C(time)), AIMAG(C(time))
       WRITE(7,'(f20.10)')   REAL(C(time))

       
    END DO time_loop

    CLOSE(6)
    CLOSE(7)

  END SUBROUTINE acf_order_n
  

END MODULE dynpro_rax
