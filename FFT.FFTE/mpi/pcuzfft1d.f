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
C     PARALLEL 1-D COMPLEX FFT ROUTINE (FOR NVIDIA GPUS)
C
C     CUDA FORTRAN + MPI SOURCE PROGRAM
C
C     CALL PZFFT1D(A,B,W,N,ICOMM,ME,NPU,IOPT)
C
C     W(N/NPU) IS COEFFICIENT VECTOR (COMPLEX*16)
C     N IS THE LENGTH OF THE TRANSFORMS (INTEGER*8)
C       -----------------------------------
C         N = (2**IP) * (3**IQ) * (5**IR)
C       -----------------------------------
C     ICOMM IS THE COMMUNICATOR (INTEGER*4)
C     ME IS THE RANK (INTEGER*4)
C     NPU IS THE NUMBER OF PROCESSORS (INTEGER*4)
C     IOPT = 0 CREATE AN FFT PLAN (INTEGER*4)
C     IOPT = -1 FOR FORWARD TRANSFORM WHERE
C              A(N/NPU) IS COMPLEX INPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE A(BLOCK)
C              B(N/NPU) IS COMPLEX OUTPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE B(BLOCK)
C     IOPT = +1 FOR INVERSE TRANSFORM WHERE
C              A(N/NPU) IS COMPLEX INPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE A(BLOCK)
C              B(N/NPU) IS COMPLEX OUTPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE B(BLOCK)
C     IOPT = 3 DESTROY THE FFT PLAN
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE PZFFT1D(A,B,W,N,ICOMM,ME,NPU,IOPT)
      use cudafor
      use cufft
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'mpif.h'
      COMPLEX*16 A(*),B(*),W(*)
      complex(8),device,allocatable :: A_d(:),B_d(:)
      INTEGER*8 N
      INTEGER*4 PLANX,PLANY
      SAVE PLANX,PLANY
C
      CALL PGETNXNY(N,NX,NY,NPU)
C
      NNX=NX/NPU
      NNY=NY/NPU
      NN=NX*NNY
C
      IF (IOPT .EQ. 0) THEN
        istat=cufftPlan1D(PLANX,NX,CUFFT_Z2Z,NNY)
        istat=cufftPlan1D(PLANY,NY,CUFFT_Z2Z,NNX)
        RETURN
      END IF
C
      IF (IOPT .EQ. 3) THEN
        istat=cufftDestroy(PLANX)
        istat=cufftDestroy(PLANY)
        RETURN
      END IF
C
      ALLOCATE(A_d(NN),B_d(NN))
      A_d=A(1:NN)
C
      IF (IOPT .EQ. 1) THEN
!$cuf kernel do <<<*,*>>>
        DO 10 I=1,NN
          A_d(I)=DCONJG(A_d(I))
   10   CONTINUE
      END IF
C
      CALL MZTRANS(A_d,B_d,NNX,NPU,NNY)
      CALL MPI_ALLTOALL(B_d,NN/NPU,MPI_DOUBLE_COMPLEX,A_d,NN/NPU,
     1                  MPI_DOUBLE_COMPLEX,ICOMM,IERR)
      CALL ZTRANS(A_d,B_d,NNX,NY)
      istat=cufftExecZ2Z(PLANY,B_d,B_d,CUFFT_FORWARD)
      CALL PMZTRANSMUL(B_d,A_d,NNY,NPU,NNX,ME,NPU)
      CALL MPI_ALLTOALL(A_d,NN/NPU,MPI_DOUBLE_COMPLEX,B_d,NN/NPU,
     1                  MPI_DOUBLE_COMPLEX,ICOMM,IERR)
      CALL ZTRANS(B_d,A_d,NNY,NX)
      istat=cufftExecZ2Z(PLANX,A_d,A_d,CUFFT_FORWARD)
      CALL MZTRANS(A_d,B_d,NNX,NPU,NNY)
      CALL MPI_ALLTOALL(B_d,NN/NPU,MPI_DOUBLE_COMPLEX,A_d,NN/NPU,
     1                  MPI_DOUBLE_COMPLEX,ICOMM,IERR)
      CALL ZTRANS(A_d,B_d,NNX,NY)
C
      IF (IOPT .EQ. 1) THEN
        DN=1.0D0/DBLE(N)
!$cuf kernel do <<<*,*>>>
        DO 20 I=1,NN
          B_d(I)=DCONJG(B_d(I))*DN
   20   CONTINUE
      END IF
C
      B(1:NN)=B_d
      DEALLOCATE(A_d,B_d)
      RETURN
      END
      SUBROUTINE PMZTRANSMUL(A_d,B_d,NS,NX,NY,ME,NPU)
      use cudafor
      IMPLICIT REAL*8 (A-H,O-Z)
      complex(8),device :: A_d(NS,NX,*),B_d(NS,NY,*)
C
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(NS)*DBLE(NX)*DBLE(NY)*DBLE(NPU))
!$cuf kernel do(3) <<<*,*>>>
      DO 30 I=1,NX
        DO 20 J=1,NY
          DO 10 K=1,NS
            TEMP=PX*DBLE((K-1)+(I-1)*NS)*(DBLE(J-1)+DBLE(ME)*DBLE(NY))
            B_d(K,J,I)=A_d(K,I,J)*DCMPLX(DCOS(TEMP),DSIN(TEMP))
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
