C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2014, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@cs.tsukuba.ac.jp
C
C
C     FACTORIZATION ROUTINE
C
C     FORTRAN77 SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE FACTOR(N,IP)
      DIMENSION IP(*)
C
      IP(1)=0
      IP(2)=0
      IP(3)=0
      N2=N
      IF (MOD(N,2) .NE. 0 .AND. MOD(N,3) .NE. 0 .AND.
     1    MOD(N,5) .NE. 0) RETURN
   10 IF (N2 .LE. 1) RETURN
      IF (MOD(N2,2) .EQ. 0) THEN
        IP(1)=IP(1)+1
        N2=N2/2
        GO TO 10
      ELSE IF (MOD(N2,3) .EQ. 0) THEN
        IP(2)=IP(2)+1
        N2=N2/3
        GO TO 10
      ELSE IF (MOD(N2,5) .EQ. 0) THEN
        IP(3)=IP(3)+1
        N2=N2/5
        GO TO 10
      END IF
      RETURN
      END
      SUBROUTINE GETNXNY(N,NX,NY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IP(3),LNX(3),LNY(3)
C
      ISQRTN=IDINT(DSQRT(DBLE(N)))
      CALL FACTOR(N,IP)
      DO 10 I=1,3
        LNX(I)=0
   10 CONTINUE
      IRES=ISQRTN
      DO 40 K=0,(IP(3)+1)/2
        DO 30 J=0,(IP(2)+1)/2
          DO 20 I=0,(IP(1)+1)/2
            NX=(2**I)*(3**J)*(5**K)
            IF (NX .LE. ISQRTN) THEN
              IRES2=ISQRTN-NX
              IF (IRES2 .LT. IRES) THEN
                LNX(1)=I
                LNX(2)=J
                LNX(3)=K
                IRES=IRES2
              END IF
            END IF
   20     CONTINUE
   30   CONTINUE
   40 CONTINUE
      DO 50 I=1,3
        LNY(I)=IP(I)-LNX(I)
   50 CONTINUE
      NX=(2**LNX(1))*(3**LNX(2))*(5**LNX(3))
      NY=(2**LNY(1))*(3**LNY(2))*(5**LNY(3))
      RETURN
      END
