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
C     2-D COMPLEX FFT ROUTINE (FOR NVIDIA GPUS)
C
C     CUDA FORTRAN SOURCE PROGRAM
C
C     CALL ZFFT2D(A,NX,NY,IOPT)
C
C     A(NX,NY) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16)
C     NX IS THE LENGTH OF THE TRANSFORMS IN THE X-DIRECTION (INTEGER*4)
C     NY IS THE LENGTH OF THE TRANSFORMS IN THE Y-DIRECTION (INTEGER*4)
C       ------------------------------------
C         NX = (2**IP) * (3**IQ) * (5**IR)
C         NY = (2**JP) * (3**JQ) * (5**JR)
C       ------------------------------------
C     IOPT = 0 CREATE AN FFT PLAN (INTEGER*4)
C          = -1 FOR FORWARD TRANSFORM
C          = +1 FOR INVERSE TRANSFORM
C          = 3 DESTROY THE FFT PLAN
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE ZFFT2D(A,NX,NY,IOPT)
      use cudafor
      use cufft
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*)
      complex(8),device,allocatable :: A_d(:),B_d(:)
      INTEGER*4 PLANX,PLANY
      SAVE PLANX,PLANY
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
      ALLOCATE(A_d(NX*NY),B_d(NX*NY))
      A_d=A(1:NX*NY)
C
      IF (IOPT .EQ. 1) THEN
!$cuf kernel do <<<*,*>>>
        DO 10 I=1,NX*NY
          A_d(I)=DCONJG(A_d(I))
   10   CONTINUE
      END IF
C
      CALL ZTRANS(A_d,B_d,NX,NY)
      istat=cufftExecZ2Z(PLANY,B_d,B_d,CUFFT_FORWARD)
      CALL ZTRANS(B_d,A_d,NY,NX)
      istat=cufftExecZ2Z(PLANX,A_d,A_d,CUFFT_FORWARD)
C
      IF (IOPT .EQ. 1) THEN
        DN=1.0D0/(DBLE(NX)*DBLE(NY))
!$cuf kernel do <<<*,*>>>
        DO 20 I=1,NX*NY
          A_d(I)=DCONJG(A_d(I))*DN
   20   CONTINUE
      END IF
C
      A(1:NX*NY)=A_d
      DEALLOCATE(A_d,B_d)
      RETURN
      END
