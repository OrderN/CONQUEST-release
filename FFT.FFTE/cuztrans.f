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
C     TRANSPOSE ROUTINE (FOR NVIDIA GPUS)
C
C     CUDA FORTRAN SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE ZTRANS(A_d,B_d,NX,NY)
      use cudafor
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      complex(8),device :: A_d(NX,*),B_d(NY,*),C_d(NB,NB)
C
      DO 60 II=1,NX,NB
        DO 50 JJ=1,NY,NB
!$cuf kernel do(2) <<<*,*>>>
          DO 20 I=II,MIN0(II+NB-1,NX)
            DO 10 J=JJ,MIN0(JJ+NB-1,NY)
              C_d(J-JJ+1,I-II+1)=A_d(I,J)
   10       CONTINUE
   20     CONTINUE
!$cuf kernel do(2) <<<*,*>>>
          DO 40 I=II,MIN0(II+NB-1,NX)
            DO 30 J=JJ,MIN0(JJ+NB-1,NY)
              B_d(J,I)=C_d(J-JJ+1,I-II+1)
   30       CONTINUE
   40     CONTINUE
   50   CONTINUE
   60 CONTINUE
      RETURN
      END
      SUBROUTINE MZTRANS(A_d,B_d,NS,NX,NY)
      use cudafor
      IMPLICIT REAL*8 (A-H,O-Z)
      complex(8),device :: A_d(NS,NX,*),B_d(NS,NY,*)
C
!$cuf kernel do(3) <<<*,*>>>
      DO 30 I=1,NX
        DO 20 J=1,NY
          DO 10 K=1,NS
            B_d(K,J,I)=A_d(K,I,J)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
