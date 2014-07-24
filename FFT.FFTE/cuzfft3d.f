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
C     3-D COMPLEX FFT ROUTINE (FOR NVIDIA GPUS)
C
C     CUDA FORTRAN SOURCE PROGRAM
C
C     CALL ZFFT3D(A,NX,NY,NZ,IOPT)
C
C     A(NX,NY,NZ) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16)
C     NX IS THE LENGTH OF THE TRANSFORMS IN THE X-DIRECTION (INTEGER*4)
C     NY IS THE LENGTH OF THE TRANSFORMS IN THE Y-DIRECTION (INTEGER*4)
C     NZ IS THE LENGTH OF THE TRANSFORMS IN THE Z-DIRECTION (INTEGER*4)
C       ------------------------------------
C         NX = (2**IP) * (3**IQ) * (5**IR)
C         NY = (2**JP) * (3**JQ) * (5**JR)
C         NZ = (2**KP) * (3**KQ) * (5**KR)
C       ------------------------------------
C     IOPT = 0 CREATE AN FFT PLAN (INTEGER*4)
C          = -1 FOR FORWARD TRANSFORM
C          = +1 FOR INVERSE TRANSFORM
C          = 3 DESTROY THE FFT PLAN
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE ZFFT3D(A,NX,NY,NZ,IOPT)
      use cudafor
      use cufft
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*)
      complex(8),device,allocatable :: A_d(:),B_d(:)
      INTEGER*4 PLANX,PLANY,PLANZ
      SAVE PLANX,PLANY,PLANZ
C
      IF (IOPT .EQ. 0) THEN
        istat=cufftPlan1D(PLANX,NX,CUFFT_Z2Z,NY*NZ)
        istat=cufftPlan1D(PLANY,NY,CUFFT_Z2Z,NZ*NX)
        istat=cufftPlan1D(PLANZ,NZ,CUFFT_Z2Z,NX*NY)
        RETURN
      END IF
C
      IF (IOPT .EQ. 3) THEN
        istat=cufftDestroy(PLANX)
        istat=cufftDestroy(PLANY)
        istat=cufftDestroy(PLANZ)
        RETURN
      END IF
C
      ALLOCATE(A_d(NX*NY*NZ),B_d(NX*NY*NZ))
      A_d=A(1:NX*NY*NZ)
C
      IF (IOPT .EQ. 1) THEN
!$cuf kernel do <<<*,*>>>
        DO 10 I=1,NX*NY*NZ
          A_d(I)=DCONJG(A_d(I))
   10   CONTINUE
      END IF
C
      CALL ZTRANS(A_d,B_d,NX*NY,NZ)
      istat=cufftExecZ2Z(PLANZ,B_d,B_d,CUFFT_FORWARD)
      CALL ZTRANS(B_d,A_d,NZ*NX,NY)
      istat=cufftExecZ2Z(PLANY,A_d,A_d,CUFFT_FORWARD)
      CALL ZTRANS(A_d,B_d,NY*NZ,NX)
      istat=cufftExecZ2Z(PLANX,B_d,A_d,CUFFT_FORWARD)
C
      IF (IOPT .EQ. 1) THEN
        DN=1.0D0/(DBLE(NX)*DBLE(NY)*DBLE(NZ))
!$cuf kernel do <<<*,*>>>
        DO 20 I=1,NX*NY*NZ
          A_d(I)=DCONJG(A_d(I))*DN
   20   CONTINUE
      END IF
C
      A(1:NX*NY*NZ)=A_d
      DEALLOCATE(A_d,B_d)
      RETURN
      END
