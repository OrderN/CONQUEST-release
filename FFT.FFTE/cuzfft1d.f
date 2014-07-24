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
C     1-D COMPLEX FFT ROUTINE (FOR NVIDIA GPUS)
C
C     CUDA FORTRAN SOURCE PROGRAM
C
C     CALL ZFFT1D(A,N,IOPT,B)
C
C     A(N) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16)
C     B(N*2) IS WORK/COEFFICIENT VECTOR (COMPLEX*16)
C     N IS THE LENGTH OF THE TRANSFORMS (INTEGER*4)
C       -----------------------------------
C         N = (2**IP) * (3**IQ) * (5**IR)
C       -----------------------------------
C     IOPT = 0 CREATE AN FFT PLAN (INTEGER*4)
C          = -1 FOR FORWARD TRANSFORM
C          = +1 FOR INVERSE TRANSFORM
C          = 3 DESTROY THE FFT PLAN
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE ZFFT1D(A,N,IOPT,B)
      use cudafor
      use cufft
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*)
      complex(8),device,allocatable :: A_d(:),B_d(:)
      INTEGER*4 PLANX,PLANY
      SAVE PLANX,PLANY
C
      CALL GETNXNY(N,NX,NY)
C
      IF (IOPT .EQ. 0) THEN
        istat=cufftPlan1D(PLANX,NX,CUFFT_Z2Z,NY)
        istat=cufftPlan1D(PLANY,NY,CUFFT_Z2Z,NX)
        RETURN
      END IF
C
      IF (IOPT .EQ. 3) THEN
        istat=cufftDestroy(PLANX)
        istat=cufftDestroy(PLANY)
        RETURN
      END IF
C
      ALLOCATE(A_d(N),B_d(N))
      A_d=A(1:N)
C
      IF (IOPT .EQ. 1) THEN
!$cuf kernel do <<<*,*>>>
        DO 10 I=1,N
          A_d(I)=DCONJG(A_d(I))
   10   CONTINUE
      END IF
C
      CALL ZTRANS(A_d,B_d,NX,NY)
      istat=cufftExecZ2Z(PLANY,B_d,B_d,CUFFT_FORWARD)
      CALL ZTRANSMUL(B_d,A_d,NY,NX)
      istat=cufftExecZ2Z(PLANX,A_d,A_d,CUFFT_FORWARD)
      CALL ZTRANS(A_d,B_d,NX,NY)
C
      IF (IOPT .EQ. 1) THEN
        DN=1.0D0/DBLE(N)
!$cuf kernel do <<<*,*>>>
        DO 20 I=1,N
          B_d(I)=DCONJG(B_d(I))*DN
   20   CONTINUE
      END IF
C
      A(1:N)=B_d
      DEALLOCATE(A_d,B_d)
      RETURN
      END
      SUBROUTINE ZTRANSMUL(A_d,B_d,NX,NY)
      use cudafor
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      complex(8),device :: A_d(NX,*),B_d(NY,*),C_d(NB,NB)
C
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(NX)*DBLE(NY))
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
              TEMP=PX*DBLE(I-1)*DBLE(J-1)
              B_d(J,I)=C_d(J-JJ+1,I-II+1)*DCMPLX(DCOS(TEMP),DSIN(TEMP))
   30       CONTINUE
   40     CONTINUE
   50   CONTINUE
   60 CONTINUE
      RETURN
      END
