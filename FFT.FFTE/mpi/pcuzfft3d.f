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
C     PARALLEL 3-D COMPLEX FFT ROUTINE (FOR NVIDIA GPUS)
C
C     CUDA FORTRAN + MPI SOURCE PROGRAM
C
C     CALL PZFFT3D(A,B,NX,NY,NZ,ICOMM,NPU,IOPT)
C
C     NX IS THE LENGTH OF THE TRANSFORMS IN THE X-DIRECTION (INTEGER*4)
C     NY IS THE LENGTH OF THE TRANSFORMS IN THE Y-DIRECTION (INTEGER*4)
C     NZ IS THE LENGTH OF THE TRANSFORMS IN THE Z-DIRECTION (INTEGER*4)
C       ------------------------------------
C         NX = (2**IP) * (3**IQ) * (5**IR)
C         NY = (2**JP) * (3**JQ) * (5**JR)
C         NZ = (2**KP) * (3**KQ) * (5**KR)
C       ------------------------------------
C     ICOMM IS THE COMMUNICATOR (INTEGER*4)
C     NPU IS THE NUMBER OF PROCESSORS (INTEGER*4)
C     IOPT = 0 CREATE AN FFT PLAN (INTEGER*4)
C     IOPT = -1 FOR FORWARD TRANSFORM WHERE
C              A(NX,NY,NZ/NPU) IS COMPLEX INPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE A(*,*,BLOCK)
C              B(NX,NY,NZ/NPU) IS COMPLEX OUTPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE B(*,*,BLOCK)
C     IOPT = +1 FOR INVERSE TRANSFORM WHERE
C              A(NX,NY,NZ/NPU) IS COMPLEX INPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE A(*,*,BLOCK)
C              B(NX,NY,NZ/NPU) IS COMPLEX OUTPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE B(*,*,BLOCK)
C     IOPT = -2 FOR FORWARD TRANSFORM WHERE
C              A(NX,NY,NZ/NPU) IS COMPLEX INPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE A(*,*,BLOCK)
C              B(NX/NPU,NY,NZ) IS COMPLEX OUTPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE B(BLOCK,*,*)
C     IOPT = +2 FOR INVERSE TRANSFORM WHERE
C              A(NX/NPU,NY,NZ) IS COMPLEX INPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE A(BLOCK,*,*)
C              B(NX,NY,NZ/NPU) IS COMPLEX OUTPUT VECTOR (COMPLEX*16)
C!HPF$ DISTRIBUTE B(*,*,BLOCK)
C     IOPT = 3 DESTROY THE FFT PLAN
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE PZFFT3D(A,B,NX,NY,NZ,ICOMM,NPU,IOPT)
      use cudafor
      use cufft
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*)
      complex(8),device,allocatable :: A_d(:),B_d(:)
      INTEGER*4 PLANX,PLANY,PLANZ
      SAVE PLANX,PLANY,PLANZ
C
      NNX=NX/NPU
      NNZ=NZ/NPU
      NN=NX*NY*NNZ
C
      IF (IOPT .EQ. 0) THEN
        istat=cufftPlan1D(PLANX,NX,CUFFT_Z2Z,NY*NNZ)
        istat=cufftPlan1D(PLANY,NY,CUFFT_Z2Z,NX*NNZ)
        istat=cufftPlan1D(PLANZ,NZ,CUFFT_Z2Z,NNX*NY)
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
      ALLOCATE(A_d(NN),B_d(NN))
      A_d=A(1:NN)
C
      IF (IOPT .EQ. 1 .OR. IOPT .EQ. 2) THEN
!$cuf kernel do <<<*,*>>>
        DO 10 I=1,NN
          A_d(I)=DCONJG(A_d(I))
   10   CONTINUE
      END IF
C
      IF (IOPT .EQ. -1 .OR. IOPT .EQ. -2) THEN
        CALL PZFFT3DF(A_d,B_d,NX,NY,NZ,ICOMM,NPU,IOPT,PLANX,PLANY,PLANZ)
      ELSE
        CALL PZFFT3DB(A_d,B_d,NX,NY,NZ,ICOMM,NPU,IOPT,PLANX,PLANY,PLANZ)
      END IF
C
      IF (IOPT .EQ. 1 .OR. IOPT .EQ. 2) THEN
        DN=1.0D0/(DBLE(NX)*DBLE(NY)*DBLE(NZ))
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
      SUBROUTINE PZFFT3DF(A_d,B_d,NX,NY,NZ,ICOMM,NPU,IOPT,PLANX,PLANY,
     1                    PLANZ)
      use cudafor
      use cufft
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'mpif.h'
      complex(8),device :: A_d(*),B_d(*)
      INTEGER*4 PLANX,PLANY,PLANZ
C
      NNX=NX/NPU
      NNZ=NZ/NPU
      NN=NX*NY*NNZ
C
      istat=cufftExecZ2Z(PLANX,A_d,A_d,CUFFT_FORWARD)
      CALL ZTRANS2(A_d,B_d,NX,NY,NNZ)
      istat=cufftExecZ2Z(PLANY,B_d,B_d,CUFFT_FORWARD)
      CALL ZTRANS3(B_d,A_d,NY,NNX,NPU,NNZ)
      CALL MPI_ALLTOALL(A_d,NN/NPU,MPI_DOUBLE_COMPLEX,B_d,NN/NPU,
     1                  MPI_DOUBLE_COMPLEX,ICOMM,IERR)
      CALL ZTRANS(B_d,A_d,NNX*NY,NZ)
      istat=cufftExecZ2Z(PLANZ,A_d,A_d,CUFFT_FORWARD)
      CALL ZTRANS(A_d,B_d,NZ,NNX*NY)
      IF (IOPT .EQ. -2) RETURN
      CALL MPI_ALLTOALL(B_d,NN/NPU,MPI_DOUBLE_COMPLEX,A_d,NN/NPU,
     1                  MPI_DOUBLE_COMPLEX,ICOMM,IERR)
      CALL MZTRANS(A_d,B_d,NNX,NY*NNZ,NPU)
      RETURN
      END
      SUBROUTINE PZFFT3DB(A_d,B_d,NX,NY,NZ,ICOMM,NPU,IOPT,PLANX,PLANY,
     1                    PLANZ)
      use cudafor
      use cufft
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'mpif.h'
      complex(8),device :: A_d(*),B_d(*)
      INTEGER*4 PLANX,PLANY,PLANZ
C
      NNX=NX/NPU
      NNZ=NZ/NPU
      NN=NX*NY*NNZ
C
      IF (IOPT .EQ. 1) THEN
        CALL MZTRANS(A_d,B_d,NNX,NPU,NY*NNZ)
        CALL MPI_ALLTOALL(B_d,NN/NPU,MPI_DOUBLE_COMPLEX,A_d,NN/NPU,
     1                    MPI_DOUBLE_COMPLEX,ICOMM,IERR)
      END IF
C
      CALL ZTRANS(A_d,B_d,NNX*NY,NZ)
      istat=cufftExecZ2Z(PLANZ,B_d,B_d,CUFFT_FORWARD)
      CALL ZTRANS(B_d,A_d,NZ,NNX*NY)
      CALL MPI_ALLTOALL(A_d,NN/NPU,MPI_DOUBLE_COMPLEX,B_d,NN/NPU,
     1                  MPI_DOUBLE_COMPLEX,ICOMM,IERR)
      CALL ZTRANS3(B_d,A_d,NNX,NY,NNZ,NPU)
      istat=cufftExecZ2Z(PLANY,A_d,A_d,CUFFT_FORWARD)
      CALL ZTRANS2(A_d,B_d,NY,NX,NNZ)
      istat=cufftExecZ2Z(PLANX,B_d,B_d,CUFFT_FORWARD)
      RETURN
      END
      SUBROUTINE ZTRANS2(A_d,B_d,NX,NY,NZ)
      use cudafor
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      complex(8),device :: A_d(NX,NY,*),B_d(NY,NX,*),C_d(NB,NB)
C
      DO 70 K=1,NZ
        DO 60 II=1,NX,NB
          DO 50 JJ=1,NY,NB
!$cuf kernel do(2) <<<*,*>>>
            DO 20 I=II,MIN0(II+NB-1,NX)
              DO 10 J=JJ,MIN0(JJ+NB-1,NY)
                C_d(J-JJ+1,I-II+1)=A_d(I,J,K)
   10         CONTINUE
   20       CONTINUE
!$cuf kernel do(2) <<<*,*>>>
            DO 40 I=II,MIN0(II+NB-1,NX)
              DO 30 J=JJ,MIN0(JJ+NB-1,NY)
                B_d(J,I,K)=C_d(J-JJ+1,I-II+1)
   30         CONTINUE
   40       CONTINUE
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE
      RETURN
      END
      SUBROUTINE ZTRANS3(A_d,B_d,NX,NY,NZ,NW)
      use cudafor
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.h'
      complex(8),device :: A_d(NX,NY,NZ,*),B_d(NY,NX,NW,*),C_d(NB,NB)
C
      DO 80 K=1,NZ
        DO 70 L=1,NW
          DO 60 II=1,NX,NB
            DO 50 JJ=1,NY,NB
!$cuf kernel do(2) <<<*,*>>>
              DO 20 I=II,MIN0(II+NB-1,NX)
                DO 10 J=JJ,MIN0(JJ+NB-1,NY)
                  C_d(J-JJ+1,I-II+1)=A_d(I,J,K,L)
   10           CONTINUE
   20         CONTINUE
!$cuf kernel do(2) <<<*,*>>>
              DO 40 I=II,MIN0(II+NB-1,NX)
                DO 30 J=JJ,MIN0(JJ+NB-1,NY)
                  B_d(J,I,L,K)=C_d(J-JJ+1,I-II+1)
   30           CONTINUE
   40         CONTINUE
   50       CONTINUE
   60     CONTINUE
   70   CONTINUE
   80 CONTINUE
      RETURN
      END
