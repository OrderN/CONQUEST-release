module dynpro_fft

  use kinds, only      : DP
  use constants,  only : pi, eps12
  use constants,  only : ps   => AU_PS
  use constants,  only : c    => C_SI
  use constants,  only : kb   => K_BOLTZMANN_SI
  use constants,  only : hp   => H_PLANCK_SI
  use constants,  only : bohr => BOHR_RADIUS_ANGS
  use constants,  only : ry   => AUTOEV

  implicit none

contains


  function EVLMEM(FDT,COF,M,PM)

    implicit none

    integer                      :: M
    real(DP)                     :: PM,FDT,COF(M)
    real(DP)                     :: WR,WI,WPR,WPI,WTEMP,THETA
    real(DP)                     :: SUMI,SUMR,EVLMEM
    integer                      :: I

    THETA=6.28318530717959D0*FDT
    WPR=DCOS(THETA)
    WPI=DSIN(THETA)
    WR=1.D0
    WI=0.D0
    SUMR=1.d0
    SUMI=0.d0
    do 11 I=1,M
       WTEMP=WR
       WR=WR*WPR-WI*WPI
       WI=WI*WPR+WTEMP*WPI
       SUMR=SUMR-COF(I)*SNGL(WR)
       SUMI=SUMI-COF(I)*SNGL(WI)
11  end do
    EVLMEM=PM/(SUMR**2+SUMI**2)
    return

  end function EVLMEM


  subroutine MEMCOF(data,N,M,PM,COF,WK1,WK2,WKM)

    !==-----------------------------------------------------------------------==
    !    MEMCOF = "Maximum Entropy Method" 
    !==-----------------------------------------------------------------------==
    !
    !  From: Numerical Recipies in Fortran 77, vol. 1, p. 561
    !
    ! " ...given a real vector data(1...n), ad given m, this routine returns 
    !      vector cof(1...m) with cof(j) = a_j, and a scalar pm=a0, which are
    !      the coefficient for Maximum Entropy Method spectral estimation..."
    ! 
    !  Remark:
    !      enforced variable declaration have been considered here
    !
    !==-----------------------------------------------------------------------==


    implicit none

    integer,           intent(in)     :: N,M
    real(DP),          intent(out)    :: PM,COF(M)
    real(DP),          intent(inout)  :: data(N),WK1(N),WK2(N),WKM(M)    

    integer                           :: I,J,K
    real(DP)                          :: P,PNEUM,DENOM

    P=0.d0 
    do 11 J=1,N
       P=P+data(J)**2
11  end do
    PM=P/N
    WK1(1)=data(1)
    WK2(N-1)=data(N)
    do 12 J=2,N-1
       WK1(J)=data(J)
       WK2(J-1)=data(J)
12  end do
    do 17 K=1,M
       PNEUM=0.d0
       DENOM=0.d0
       do 13 J=1,N-K
          PNEUM=PNEUM+WK1(J)*WK2(J)
          DENOM=DENOM+WK1(J)**2+WK2(J)**2
13     end do
       COF(K)=2.d0*PNEUM/DENOM
       PM=PM*(1.d0-COF(K)**2)
       if(K.ne.1)then
          do 14 I=1,K-1
             COF(I)=WKM(I)-COF(K)*WKM(K-I)
14        end do
       endif
       if(K.eq.M)return
       do 15 I=1,K
          WKM(I)=COF(I)
15     end do
       do 16 J=1,N-K-1
          WK1(J)=WK1(J)-WKM(K)*WK2(J)
          WK2(J)=WK2(J+1)-WKM(K)*WK1(J+1)
16     end do
17  end do

    !PAUSE 'never get here'

  end subroutine MEMCOF

  subroutine SPCTRM(P,M,K,OVRLAP,W1,W2,unit_vac,filevac,win,win_name)

    implicit none

    integer,           intent(in)     :: K,M,unit_vac
    character(len=256),intent(inout)  :: filevac
    real(DP),          intent(inout)  :: P(M),W1(4*M),W2(m)
    logical,           intent(in)     :: OVRLAP

    integer                           :: I,J,J2,JOFF,JOFFN,KK,M4,M43,M44,MM
    real(DP)                          :: DEN,FACM,FACP,SUMW,W
    real(DP),          allocatable    :: WINDOW(:)

    integer,           intent(in)     :: win
    character(11),     intent(out)    :: win_name

    MM   = M  + M
    M4   = MM + MM
    M44  = M4 + 4
    M43  = M4 + 3
    DEN  = 0.d0
    FACM = M  - 0.5d0
    FACP = 1.d0/(M + 0.5d0)
    SUMW = 0.d0

    allocate(WINDOW(MM))

    select case(win)

    case(1)
       do J = 1, MM
          WINDOW(J) = 1.d0 - abs(((J - 1) - FACM)*FACP)
       end do
       win_name = 'Bartlett   '

    case(2)
       WINDOW = 1.0d0      
       win_name = 'rectangular'

    case(3)    
       do J = 1, MM
          WINDOW(J) = 1.d0 - (((J - 1) - FACM)*FACP)**2
       end do
       win_name = 'Welch      '

    case(4)   
       do J = 1, MM
          WINDOW(J) = 0.5d0*(1.d0 - cos(pi*j/FACM))
       end do
       win_name = 'Hann       '    

    case DEFAULT
       WINDOW = 1.0d0      
       win_name = 'rectangular'

    end select

    open(unit_vac, file=filevac, status='unknown', position='rewind')

    do 11 J=1,MM
       SUMW=SUMW+WINDOW(J)**2
11  end do
    do 12 J=1,M
       P(J)=0.d0
12  end do

    if(OVRLAP)then
       read(unit_vac,'(f20.10)') (W2(J),J=1,M)
    endif

    do 18 KK=1,K
       do 15 JOFF=-1,0,1
          if (OVRLAP) then
             do 13 J=1,M
                W1(JOFF+J+J)=W2(J)
13           end do
             read(unit_vac,'(f20.10)') (W2(J),J=1,M)
             JOFFN=JOFF+MM
             do 14 J=1,M
                W1(JOFFN+J+J)=W2(J)
14           end do
          else
             read(unit_vac,'(f20.10)') (W1(J),J=JOFF+2,M4,2)
          endif
15     end do
       do 16 J=1,MM
          J2=J+J
          W=WINDOW(J)
          W1(J2)=W1(J2)*W
          W1(J2-1)=W1(J2-1)*W
16     end do
       call FOUR1(W1,MM,1)
       P(1)=P(1)+W1(1)**2+W1(2)**2
       do 17 J=2,M
          J2=J+J
          P(J)=P(J)+W1(J2)**2+W1(J2-1)**2+W1(M44-J2)**2+W1(M43-J2)**2
17     end do
       DEN=DEN+SUMW
18  end do

    DEN=M4*DEN

    do 19 J=1,M
       P(J)=P(J)/DEN
19  end do

    close(unit_vac)

    deallocate(WINDOW)

    return
  end subroutine SPCTRM

  subroutine FOUR1(data,NN,ISIGN)

    implicit none

    integer,   intent(in)    :: NN,ISIGN
    real(DP),  intent(inout) :: data(*)

    real(DP)                 :: WR,WI,WPR,WPI,WTEMP
    real(DP)                 :: TEMPR,TEMPI,THETA
    integer                  :: N,M,MMAX,ISTEP
    integer                  :: I,J

    N=2*NN
    J=1

    do 11 I=1,N,2
       if(J.gt.I)then
          TEMPR=data(J)
          TEMPI=data(J+1)
          data(J)=data(I)
          data(J+1)=data(I+1)
          data(I)=TEMPR
          data(I+1)=TEMPI
       endif
       M=N/2
1      if ((M.ge.2).and.(J.gt.M)) then
          J=J-M
          M=M/2
          GO TO 1
       endif
       J=J+M
11  end do
    MMAX=2
2   if (N.gt.MMAX) then
       ISTEP=2*MMAX
       THETA=6.28318530717959D0/(ISIGN*MMAX)
       WPR=-2.D0*DSIN(0.5D0*THETA)**2
       WPI=DSIN(THETA)
       WR=1.D0
       WI=0.D0
       do 13 M=1,MMAX,2
          do 12 I=M,N,ISTEP
             J=I+MMAX
             TEMPR=SNGL(WR)*data(J)-SNGL(WI)*data(J+1)
             TEMPI=SNGL(WR)*data(J+1)+SNGL(WI)*data(J)
             data(J)=data(I)-TEMPR
             data(J+1)=data(I+1)-TEMPI
             data(I)=data(I)+TEMPR
             data(I+1)=data(I+1)+TEMPI
12        end do
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13     end do
       MMAX=ISTEP
       GO TO 2
    endif

    return

  end subroutine FOUR1


  subroutine REALFT(data,N,ISIGN)

    implicit none

    integer,   intent(in)    :: N,ISIGN
    real(DP),  intent(inout) :: data(2*N)

    real(DP)                 :: WR,WRS,WI,WIS,WPR,WPI,WTEMP
    real(DP)                 :: THETA,C1,C2 
    real(DP)                 :: H1R,H1I,H2R,H2I
    integer                  :: I,J,I1,I2,I3,I4,N2P3

    THETA=6.28318530717959D0/2.0D0/dble(N)
    C1=0.5

    if (ISIGN.eq.1) then
       C2=-0.5
       call FOUR1(data,N,+1)
    else
       C2=0.5
       THETA=-THETA
    endif
    WPR=-2.0D0*DSIN(0.5D0*THETA)**2
    WPI=DSIN(THETA)
    WR=1.0D0+WPR
    WI=WPI
    N2P3=2*N+3
    do 11 I=2,N/2+1
       I1=2*I-1
       I2=I1+1
       I3=N2P3-I2
       I4=I3+1
       WRS=SNGL(WR)
       WIS=SNGL(WI)
       H1R=C1*(data(I1)+data(I3))
       H1I=C1*(data(I2)-data(I4))
       H2R=-C2*(data(I2)+data(I4))
       H2I=C2*(data(I1)-data(I3))
       data(I1)=H1R+WRS*H2R-WIS*H2I
       data(I2)=H1I+WRS*H2I+WIS*H2R
       data(I3)=H1R-WRS*H2R+WIS*H2I
       data(I4)=-H1I+WRS*H2I+WIS*H2R
       WTEMP=WR
       WR=WR*WPR-WI*WPI+WR
       WI=WI*WPR+WTEMP*WPI+WI
11  end do
    if (ISIGN.eq.1) then
       H1R=data(1)
       data(1)=H1R+data(2)
       data(2)=H1R-data(2)
    else
       H1R=data(1)
       data(1)=C1*(H1R+data(2))
       data(2)=C1*(H1R-data(2))
       call FOUR1(data,N,-1)
    endif

    return
  end subroutine REALFT


end module dynpro_fft
