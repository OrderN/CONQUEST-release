MODULE dynpro_ran

  USE kinds,      ONLY : DP

  IMPLICIT NONE

CONTAINS

  SUBROUTINE random()
    
    USE dynpro_inp,  ONLY : fileran, filetoto
    USE dynpro_inp,  ONLY : iteration
    
    
    IMPLICIT NONE
    
    REAL(DP),   ALLOCATABLE   :: H(:,:),  Q(:,:,:)
    INTEGER,    ALLOCATABLE   :: SEED(:)
    REAL(DP)                  :: x, y

    INTEGER                   :: i, j, k, n, m, p, length

    length = 10000
    
    n =  iteration
    m =  3
    p = 10
    y = 0.0d0

    ALLOCATE(H(m,m))
    ALLOCATE(SEED(p))
    ALLOCATE(Q(n,m,m))


    OPEN(20, file=fileran, status='unknown', position='rewind')

    OPEN(21, file=filetoto, status='unknown', position='rewind')


    H = 0.0d0
    Q = 0.0d0
    x = 0.0d0

    !DO i = 1,n 
    !   
    !   READ(21,*) (H(k,k), k=1,m)      
    !   DO j =1, m
    !      WRITE(20,'(3F14.6)') (1000.0d0*H(j,k), k=1,m)
    !   END DO
       
    !END DO
    
    !DO i = 1,n 
       
    !   DO j =1, m
    !      READ(21,*) (H(j,k), k=1,m)      
    !      WRITE(20,'(3F14.6)') (H(j,k), k=1,m)
    !   END DO
       
    !END DO
    

    !CALL RANDOM_SEED(SIZE=length)
    
    DO i = 1, n
       DO k = 1, 50000000
          CALL RANDOM_NUMBER(H)
       END DO
       DO j = 1, m
          WRITE(20,'(3F14.6)') (H(j,k), k=1,m)
       END DO
    END DO
    
    
    
    

    !DO i = 1, n
    !   DO j = 1, m
    !      DO k = 1, m
    !         H(j,k) = cos(H(j,k)+x) 
    !         WRITE(20,'(3F14.6)') H(j,k)
    !      END DO
    !   END DO
    !   x = x + 0.5
    !END DO

    CLOSE(20)
    CLOSE(21)

    DEALLOCATE(H, SEED)

  END SUBROUTINE random
  
END MODULE dynpro_ran
