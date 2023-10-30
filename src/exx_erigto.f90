module exx_erigto

  use datatypes
  use numbers,  ONLY: pi,zero,one,two,four,half,three_halves


  USE ISO_C_BINDING, ONLY: C_DOUBLE, C_F_POINTER, C_F_PROCPOINTER, C_NULL_PTR
  !use libint_f,  only: libint_t, libint2_static_init, libint2_static_cleanup,&
  !     libint2_build, libint2_max_am_eri, libint2_cleanup_eri, &
  !     compute_eri_f, libint2_init_eri, print_eri

  
  implicit none
  
  !TYPE(libint_t), DIMENSION(1) :: erieval

  !INTEGER, PARAMETER :: dp = C_DOUBLE

  
  ! The equations herein are based upon
  ! 'Gaussian Expansion Methods for Molecular Orbitals' H. Taketa,
  ! S. Huzinaga, and K. O-ohata. H. Phys. Soc. Japan, 21, 2313, 1966.
  !
contains
  ! Functions in this module:
  ! eri_gto:     Compute the coulomb repulsion between four primitive gaussian functions
  ! overlap_gto: Compute the overlap between two primitive gaussians
  !
  function eri_gto_hoh( &
       xa, ya, za, norma, la, ma, na, alphaa, &
       xb, yb, zb, normb, lb, mb, nb, alphab, &
       xc, yc, zc, normc, lc, mc, nc, alphac, &
       xd, yd, zd, normd, ld, md, nd, alphad)
    ! use functions:   dist2, product_center_1D, Fgamma
    ! use subroutines: B_array
    implicit none
    ! Funtion
    real(double) :: eri_gto_hoh

    ! Parameters
    integer, parameter :: MAXLEN=25

    ! Arguments
    real(double), intent(in) :: xa,ya,za,norma,alphaa,xb,yb,zb,normb,alphab
    real(double), intent(in) :: xc,yc,zc,normc,alphac,xd,yd,zd,normd,alphad
    integer,      intent(in) :: la,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd

    ! Local variables
    real(double) :: rab2,rcd2,rpq2
    real(double) :: xp,yp,zp,xq,yq,zq 
    real(double) :: gamma1,gamma2,delta,sum,xval

    real(double) ::  Bx(0:MAXLEN-1), By(0:MAXLEN-1), Bz(0:MAXLEN-1)
    integer      ::  i,j,k


    Bx = zero ; By = zero ; Bz = zero

    rab2 = dist2(xa,ya,za,xb,yb,zb)
    rcd2 = dist2(xc,yc,zc,xd,yd,zd)

    xp = product_center_1D(alphaa,xa,alphab,xb)
    yp = product_center_1D(alphaa,ya,alphab,yb)
    zp = product_center_1D(alphaa,za,alphab,zb)
    xq = product_center_1D(alphac,xc,alphad,xd)
    yq = product_center_1D(alphac,yc,alphad,yd)
    zq = product_center_1D(alphac,zc,alphad,zd)

    rpq2 = dist2(xp,yp,zp,xq,yq,zq)

    gamma1 = alphaa+alphab
    gamma2 = alphac+alphad

    delta = (one/gamma1+one/gamma2)/four

    call B_array(la+lb+lc+ld+1,Bx, &
         la,lb,lc,ld,xp,xa,xb,xq,xc,xd,gamma1,gamma2,delta)
    call B_array(ma+mb+mc+md+1,By, &
         ma,mb,mc,md,yp,ya,yb,yq,yc,yd,gamma1,gamma2,delta)
    call B_array(na+nb+nc+nd+1,Bz, &
         na,nb,nc,nd,zp,za,zb,zq,zc,zd,gamma1,gamma2,delta)

    sum = zero
    do i = 0,la+lb+lc+ld
       do j = 0, ma+mb+mc+md
          do k = 0, na+nb+nc+nd
             xval = float(i+j+k) ! cast to double
             sum  = sum + &
                  Bx(i)*By(j)*Bz(k) &
                  *Fgamma(xval,0.25d0*rpq2/delta)
             !print*, Bx(i), By(j), Bz(k),Fgamma(xval,0.25d0*rpq2/delta)
          enddo
       enddo
    enddo

    !print*, 'sum =', sum

    eri_gto_hoh = two*pi**2.5d0 &
         /(gamma1*gamma2*sqrt(gamma1+gamma2)) &
         *exp(-alphaa*alphab*rab2/gamma1) &
         *exp(-alphac*alphad*rcd2/gamma2) &
         *sum*norma*normb*normc*normd     

    return
  end function eri_gto_hoh
  !
  function overlap_gto(alpha1,l1,m1,n1,x1,y1,z1, &
       alpha2,l2,m2,n2,x2,y2,z2)
    ! Note: this function does multiply by the normalization constant    
    ! Use functions: product_center_1D,dist2,overlap1D
    implicit none
    ! Function
    real(double) :: overlap_gto
    ! Arguments
    real(double), intent(in) :: alpha1,x1,y1,z1,alpha2,x2,y2,z2
    integer,      intent(in) :: l1,m1,n1,l2,m2,n2
    ! Local variables
    real(double) :: rab2,gamma,xp,yp,zp,factor,wx,wy,wz

    rab2 = dist2(x1,y1,z1,x2,y2,z2)

    gamma = alpha1+alpha2
    xp = product_center_1D(alpha1,x1,alpha2,x2)
    yp = product_center_1D(alpha1,y1,alpha2,y2)
    zp = product_center_1D(alpha1,z1,alpha2,z2)

    factor = ((pi/gamma)**1.5d0)*exp(-alpha1*alpha2*rab2/gamma)

    wx = overlap1D(l1,l2,xp-x1,xp-x2,gamma)
    wy = overlap1D(m1,m2,yp-y1,yp-y2,gamma)
    wz = overlap1D(n1,n2,zp-z1,zp-z2,gamma)
    overlap_gto = factor*wx*wy*wz;
    return
  end function overlap_gto
  !
  function norm_gto(alpha,l,m,n)
    ! Compute the norm of any Cartesian Gaussian function
    implicit none
    ! Function
    real(double) :: norm_gto
    ! Arguments
    real(double), intent(in) :: alpha
    integer,      intent(in) :: l,m,n                
    ! Local variables
    real(double) :: numo,deno

    numo = two**(two*(l+m+n)+1.5d0) &
         *alpha**((l+m+n)+1.5d0)
    deno = fact2(2*l-1) &
         *fact2(2*m-1)  &
         *fact2(2*n-1)  &         
         *pi**1.5d0
    norm_gto = sqrt(numo/deno)
    return
  end function norm_gto
  !
  function overlap1D(l1,l2,xap,xbp,gamma)
    ! Use functions: binomial_prefactor,fact2
    implicit none
    ! Function
    real(double) :: overlap1D
    ! Arguments
    integer,      intent(in) :: l1,l2
    real(double), intent(in) :: xap,xbp,gamma
    ! Local variables
    integer      :: i
    real(double) :: sum
    sum = zero

    do i=0, int(half*(l1+l2))
       sum = sum + &
            binomial_prefactor(2*i,l1,l2,xap,xbp)* &
            fact2(2*i-1)/((2*gamma)**i)
    enddo
    overlap1D = sum
    return 
  end function overlap1D
  !
  function dist2(xa,ya,za,xb,yb,zb)

    implicit none
    ! Function
    real(double) :: dist2
    ! Arguments
    real(double), intent(in) :: xa,ya,za,xb,yb,zb
    ! Local variables
    real(double) :: dx,dy,dz

    dx = xa-xb
    dy = ya-yb
    dz = za-zb
    dist2 = dx*dx+dy*dy+dz*dz
    return
  end function dist2
  !
  function product_center_1D(alphaa,xa,alphab,xb)

    implicit none
    real(double) :: product_center_1D    
    real(double), intent(in) :: alphaa,xa,alphab,xb

    product_center_1D = (alphaa*xa+alphab*xb)/(alphaa+alphab)
    return 
  end function product_center_1D
  !    
  function Fgamma(a,x)
    ! use function: gammap
    implicit none
    real(double) :: a,x
    real(double) :: Fgamma,val

    if (x .lt. 0.00001d0) x = 0.000001d0    
    val = gammp(a+half,x)
    Fgamma = half*(x**(-a-half))*val
    return
  end function Fgamma
  !  
  subroutine B_array(n,barry,l1,l2,l3,l4, &
       p,a,b,q,c,d,g1,g2,delta)
    ! use functions: barry, b_term
    implicit none
    ! Arguments
    integer,      intent(in)    :: n,l1,l2,l3,l4
    real(double), intent(in)    :: p,a,b,q,c,d,g1,g2,delta
    real(double), intent(inout) :: barry(0:n-1)
    ! Local variables
    integer                     ::i,i1,i2,j1,j2,k,index

    barry = zero
    !do i = 0, n-1
    !   barry(i) = zero
    !enddo

    do i1=0, l1+l2
       do i2=0, l3+l4
          do j1=0, i1/2
             do j2=0, i2/2
                do k=0, (i1+i2)/2-j1-j2
                   index = i1+i2-2*(j1+j2)-k
                   barry(index) = barry(index)              &
                        + b_term(i1,i2,j1,j2,k,l1,l2,l3,l4, &
                        p,a,b,q,c,d,g1,g2,delta)
                enddo
             enddo
          enddo
       enddo
    enddo
    return 
  end subroutine B_array
  ! 
  function b_term(i1,i2,j1,j2,k,l1,l2,l3,l4, &
       p,a,b,q,c,d,g1,g2,delta)
    ! use functions: fB, fact_ratio2
    implicit none    
    ! Function
    real(double) :: b_term
    ! Arguments
    integer,      intent(in) :: i1,i2,j1,j2,k,l1,l2,l3,l4
    real(double), intent(in) :: p,a,b,q,c,d,g1,g2,delta

    !real*8 fB
    !integer fact_ratio2

    b_term = fB(i1,l1,l2,p,a,b,j1,g1)             &
         *((-1)**i2)*fB(i2,l3,l4,q,c,d,j2,g2)       &
         *((-1)**k)*fact_ratio2(i1+i2-2*(j1+j2),k)  &
         *(q-p)**(i1+i2-2*(j1+j2)-2*k)            &
         /(delta)**(i1+i2-2*(j1+j2)-k)            

    return
  end function b_term
  ! 
  function fB(i,l1,l2,px,ax,bx,ir,g)
    ! use functions: binomial_prefactor, Bfunc
    implicit none
    real(double) :: fB
    real(double), intent(in) :: px,ax,bx,g
    integer,      intent(in) :: i,l1,l2,ir

    fB = binomial_prefactor(i,l1,l2,px-ax,px-bx)*Bfunc(i,ir,g)
    return
  end function fB
  ! 
  function Bfunc(i,ir,g)
    ! use function: fact_ratio2
    implicit none
    real(double) :: Bfunc
    real(double) :: g
    integer      :: i,ir

    Bfunc = fact_ratio2(i,ir)*((4*g)**(ir-i))
    return
  end function Bfunc
  !
  function binomial_prefactor(s,ia,ib,xpa,xpb)
    ! use functions: binomial
    implicit none
    ! Function
    real(double) :: binomial_prefactor
    ! Arguments
    integer,      intent(in) :: s,ia,ib
    real(double), intent(in) :: xpa,xpb
    ! Local variables
    integer      :: j
    real(double) :: sum

    sum = zero
    do j=0,s
       if ((s-ia .le. j) .and. (j .le. ib))             &
            sum = sum + binomial(ia,s-j)*binomial(ib,j) &
            *(xpa**(ia-s+j))*(xpb**(ib-j))
    enddo
    binomial_prefactor = sum
    return
  end function binomial_prefactor
  !
  function binomial(i,j)
    ! use function: fact
    implicit none      
    integer             :: binomial
    integer, intent(in) :: i,j

    binomial = fact(i)/fact(j)/fact(i-j)      

    return
  end function binomial
  !
  function fact(n)

    implicit none
    integer :: i,n,fact

    fact=1
    if (n .le. 1) return
    do i=1,n
       fact = fact*i
    enddo
    return
  end function fact
  !
  function fact2(n)
    ! Compute double factorial n!!
    implicit none
    integer :: fact2,i,n

    fact2=1
    if (n .le. 1) return
    do i = 1,n,2
       fact2 = fact2*i
    enddo
    return
  end function fact2
  !
  function fact_ratio2(i,j)
    ! use function: fact
    implicit none
    integer, intent(in) :: i,j
    integer :: fact_ratio2

    fact_ratio2 = fact(i)/fact(j)/fact(i-2*j)
    return
  end function fact_ratio2
  !
  ! The following are from Numerical Recipes
  ! gammp is hacked a bit to return exp(gln)*gamma
  function gammp(a,x)
    ! use subroutines: gcf, gser
    implicit none
    real(double), intent(in) :: a,x
    real(double) :: gammp
    real(double) :: gammcf,gamser,gln

    if(x.lt.0..or.a.le.0.) stop 'exx_erigto: error, bad arguments in gammp'
    if(x.lt.a+one)then
       call gser(gamser,a,x,gln)
       gammp=exp(gln)*gamser
    else
       call gcf(gammcf,a,x,gln)
       gammp=exp(gln)*(one-gammcf)
    endif
    return
  end function gammp
  !
  subroutine gser(gamser,a,x,gln)

    implicit none
    ! Arguments
    real(double), intent(in)  :: a,x
    real(double), intent(out) :: gamser,gln
    ! Local parameters
    real(double),parameter    :: EPS=3.e-9
    integer,     parameter    :: ITMAX=100
    ! Local variables
    integer                   :: n
    real(double)              :: sum,del,ap

    gln=gammln(a)
    if(x.le.0.)then
       if(x.lt.0.) stop 'exx_erigto: error, x < 0 in gser'
       gamser=zero
       return
    endif
    ap =a
    sum=one/a
    del=sum
    do 11 n=1,ITMAX
       ap =ap+one
       del=del*x/ap
       sum=sum+del
       if(abs(del).lt.abs(sum)*EPS) goto 1
11  end do
    stop 'exx_erigto: error, a too large, ITMAX too small in gser'
1   gamser=sum*exp(-x+a*log(x)-gln)
    return
  end subroutine gser
  !
  function gammln(xx)

    implicit none
    ! Function
    real(double) :: gammln
    ! Arguments
    real(double), intent(in) :: xx
    ! Local variables
    integer      :: j
    real(double) :: ser,stp,tmp,x,y,cof(6)
    SAVE cof,stp
    DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
         24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
         -.5395239384953d-5,2.5066282746310005d0/

    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do 11 j=1,6
       y=y+one
       ser=ser+cof(j)/y
11  end do
    gammln=tmp+log(stp*ser/x)
    return
  end function gammln
  !
  subroutine gcf(gammcf,a,x,gln)
    ! use function gammln
    implicit none
    ! Arguments
    real(double), intent(in)  :: a,x
    real(double), intent(out) :: gln,gammcf
    ! Local parameters
    real(double),parameter    :: EPS=3.e-9
    real(double),parameter    :: FPMIN=1.e-30
    integer,     parameter    :: ITMAX=100

    ! Local variables
    integer                   :: i
    real(double)              :: an,b,c,d,del,h

    gln=gammln(a)
    !gln=dlgamma(a)
    b=x+one-a
    c=one/FPMIN
    d=one/b
    h=d
    do 11 i=1,ITMAX
       an=-i*(i-a)
       b=b+two
       d=an*d+b
       if(abs(d).lt.FPMIN)d=FPMIN
       c=b+an/c
       if(abs(c).lt.FPMIN)c=FPMIN
       d=one/d
       del=d*c
       h=h*del
       if(abs(del-one).lt.EPS) goto 1
11  end do
    stop 'exx_erigto: error, a too large, ITMAX too small in gcf'
1   gammcf=exp(-x+a*log(x)-gln)*h
    return
  end subroutine gcf
  !
  function dlgamma(x)
    implicit real(double) (a - h, o - z)
    dimension a(0 : 21), b(0 : 97), c(0 : 64), d(0 : 6)
    integer :: i, k
    data (a(i), i = 0, 10) / &
         0.00009967270908702825d0, -0.00019831672170162227d0, &
         -0.00117085315349625822d0, 0.00722012810948319552d0, &
         -0.00962213009367802970d0, -0.04219772092994235254d0, &
         0.16653861065243609743d0, -0.04200263501129018037d0, &
         -0.65587807152061930091d0, 0.57721566490153514421d0, &
         0.99999999999999999764d0 / 
    data (a(i), i = 11, 21) / &
         0.00004672097259011420d0, -0.00006812300803992063d0, &
         -0.00132531159076610073d0, 0.00733521178107202770d0, &
         -0.00968095666383935949d0, -0.04217642811873540280d0, &
         0.16653313644244428256d0, -0.04200165481709274859d0, &
         -0.65587818792782740945d0, 0.57721567315209190522d0, &
         0.99999999973565236061d0 / 
    data (b(i), i = 0, 13) / &
         -0.00000000004587497028d0, 0.00000000019023633960d0, &
         0.00000000086377323367d0, 0.00000001155136788610d0, &
         -0.00000002556403058605d0, -0.00000015236723372486d0,& 
         -0.00000316805106385740d0, 0.00000122903704923381d0, &
         0.00002334372474572637d0, 0.00111544038088797696d0, &
         0.00344717051723468982d0, 0.03198287045148788384d0, &
         -0.32705333652955399526d0, 0.40120442440953927615d0 / 
    data (b(i), i = 14, 27) / &
         &    -0.00000000005184290387d0, -0.00000000083355121068d0, &
         &    -0.00000000256167239813d0, 0.00000001455875381397d0, &
         &    0.00000013512178394703d0, 0.00000029898826810905d0, &
         &    -0.00000358107254612779d0, -0.00002445260816156224d0,& 
         &    -0.00004417127762011821d0, 0.00112859455189416567d0, &
         &    0.00804694454346728197d0, 0.04919775747126691372d0, &
         &    -0.24818372840948854178d0, 0.11071780856646862561d0 / 
    data (b(i), i = 28, 41) / &
         &    0.00000000030279161576d0, 0.00000000160742167357d0,  &
         &    -0.00000000405596009522d0, -0.00000005089259920266d0,& 
         &    -0.00000002029496209743d0, 0.00000135130272477793d0, &
         &    0.00000391430041115376d0, -0.00002871505678061895d0, &
         &    -0.00023052137536922035d0, 0.00045534656385400747d0, &
         &    0.01153444585593040046d0, 0.07924014651650476036d0,  &
         &    -0.12152192626936502982d0, -0.07916438300260539592d0 / 
    data (b(i), i = 42, 55) / &
         &    -0.00000000050919149580d0, -0.00000000115274986907d0,&
         &    0.00000001237873512188d0, 0.00000002937383549209d0,  &
         &    -0.00000030621450667958d0, -0.00000077409414949954d0,& 
         &    0.00000816753874325579d0, 0.00002412433382517375d0,  &
         &    -0.00026061217606063700d0, -0.00091000087658659231d0,& 
         &    0.01068093850598380797d0, 0.11395654404408482305d0,  &
         &    0.07209569059984075595d0, -0.10971041451764266684d0 / 
    data (b(i), i = 56, 69) / &
         &    0.00000000040119897187d0, -0.00000000013224526679d0, &
         &    -0.00000001002723190355d0, 0.00000002569249716518d0, &
         &    0.00000020336011868466d0, -0.00000118097682726060d0, &
         &    -0.00000300660303810663d0, 0.00004402212897757763d0, &
         &    -0.00001462405876235375d0, -0.00164873795596001280d0, &
         &    0.00513927520866443706d0, 0.13843580753590579416d0,   &
         &    0.32730190978254056722d0, 0.08588339725978624973d0 / 
    data (b(i), i = 70, 83) / &
         &    -0.00000000015413428348d0, 0.00000000064905779353d0, &
         &    0.00000000160702811151d0, -0.00000002655645793815d0, &
         &    0.00000007619544277956d0, 0.00000047604380765353d0,  &
         &    -0.00000490748870866195d0, 0.00000821513040821212d0, &
         &    0.00014804944070262948d0, -0.00122152255762163238d0, &
         &    -0.00087425289205498532d0, 0.14438703699657968310d0, &
         &    0.61315889733595543766d0, 0.55513708159976477557d0 / 
    data (b(i), i = 84, 97) / &
         &    0.00000000001049740243d0, -0.00000000025832017855d0,  &
         &    0.00000000139591845075d0, -0.00000000021177278325d0,  &
         &    -0.00000005082950464905d0, 0.00000037801785193343d0,  &
         &    -0.00000073982266659145d0, -0.00001088918441519888d0, &
         &    0.00012491810452478905d0, -0.00049171790705139895d0,  &
         &    -0.00425707089448266460d0, 0.13595080378472757216d0,  &
         &    0.89518356003149514744d0, 1.31073912535196238583d0 / 
    data (c(i), i = 0, 12) / &
         &    0.0000000116333640008d0, -0.0000000833156123568d0, &
         &    0.0000003832869977018d0, -0.0000015814047847688d0, &
         &    0.0000065010672324100d0, -0.0000274514060128677d0, &
         &    0.0001209015360925566d0, -0.0005666333178228163d0, &
         &    0.0029294103665559733d0, -0.0180340086069185819d0, &
         &    0.1651788780501166204d0, 1.1031566406452431944d0,  &
         &    1.2009736023470742248d0 / 
    data (c(i), i = 13, 25) / &
         &    0.0000000013842760642d0, -0.0000000069417501176d0, &
         &    0.0000000342976459827d0, -0.0000001785317236779d0, &
         &    0.0000009525947257118d0, -0.0000052483007560905d0, &
         &    0.0000302364659535708d0, -0.0001858396115473822d0, &
         &    0.0012634378559425382d0, -0.0102594702201954322d0, &
         &    0.1243625515195050218d0, 1.3888709263595291174d0, &
         &    2.4537365708424422209d0 / 
    data (c(i), i = 26, 38) / &
         &    0.0000000001298977078d0, -0.0000000008029574890d0, &
         &    0.0000000049454846150d0, -0.0000000317563534834d0, &
         &    0.0000002092136698089d0, -0.0000014252023958462d0, &
         &    0.0000101652510114008d0, -0.0000774550502862323d0, &
         &    0.0006537746948291078d0, -0.0066014912535521830d0, &
         &    0.0996711934948138193d0, 1.6110931485817511402d0, &
         &    3.9578139676187162939d0 / 
    data (c(i), i = 39, 51) / &
         &    0.0000000000183995642d0, -0.0000000001353537034d0, &
         &    0.0000000009984676809d0, -0.0000000076346363974d0, &
         &    0.0000000599311464148d0, -0.0000004868554120177d0, &
         &    0.0000041441957716669d0, -0.0000377160856623282d0, &
         &    0.0003805693126824884d0, -0.0045979851178130194d0, &
         &    0.0831422678749791178d0, 1.7929113303999329439d0, &
         &    5.6625620598571415285d0 / 
    data (c(i), i = 52, 64) / &
         &    0.0000000000034858778d0, -0.0000000000297587783d0, &
         &    0.0000000002557677575d0, -0.0000000022705728282d0, &
         &    0.0000000207024992450d0, -0.0000001954426390917d0, &
         &    0.0000019343161886722d0, -0.0000204790249102570d0, &
         &    0.0002405181940241215d0, -0.0033842087561074799d0, &
         &    0.0713079483483518997d0, 1.9467574842460867884d0, &
         &    7.5343642367587329552d0 / 
    data (d(i), i = 0, 6) / &
         &    -0.00163312359200500807d0, 0.00083644533703385956d0, &
         &    -0.00059518947575728181d0, 0.00079365057505415415d0, &
         &    -0.00277777777735463043d0, 0.08333333333333309869d0, &
         &    0.91893853320467274178d0 / 
    w = x
    if (x .lt. 0) w = 1 - x
    if (w .lt. 0.5d0) then
       k = 0
       if (w .ge. 0.25d0) k = 11
       y = ((((((((((a(k) * w + a(k + 1)) * w + &
            a(k + 2)) * w + a(k + 3)) * w + a(k + 4)) * w + &
            a(k + 5)) * w + a(k + 6)) * w + a(k + 7)) * w + &
            a(k + 8)) * w + a(k + 9)) * w + a(k + 10)) * w  
       y = -log(y)
    else if (w .lt. 3.5d0) then
       t = w - 4.5d0 / (w + 0.5d0)
       k = int(t + 4)
       t = t - (k - 3.5d0)
       k = k * 14
       y = ((((((((((((b(k) * t + b(k + 1)) * t + &
            b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + &
            b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + &
            b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + &
            b(k + 11)) * t + b(k + 12)) * t + b(k + 13) 
    else if (w .lt. 8) then
       k = (int(w)) - 3
       t = w - (k + 3.5d0)
       k = k * 13
       y = (((((((((((c(k) * t + c(k + 1)) * t + &
            c(k + 2)) * t + c(k + 3)) * t + c(k + 4)) * t + &
            c(k + 5)) * t + c(k + 6)) * t + c(k + 7)) * t + &
            c(k + 8)) * t + c(k + 9)) * t + c(k + 10)) * t + &
            c(k + 11)) * t + c(k + 12)
    else
       v = 1 / w
       t = v * v
       y = (((((d(0) * t + d(1)) * t + d(2)) * t + &
            d(3)) * t + d(4)) * t + d(5)) * v + d(6)
       y = y + ((w - 0.5d0) * log(w) - w)
    end if
    if (x .lt. 0) y = log(pi / sin(pi * x)) - y
    dlgamma = y
  end function dlgamma
  !
  !
  !
  subroutine compute_eri_hoh(i,j,k,l,i_nsp,j_nsp,k_nsp,l_nsp,i_xyz,j_xyz,k_xyz,l_xyz,&
        i_nt, j_nt, k_nt, l_nt, eri_gto)

    use gto_format_new, only: gto
        
    implicit none

    integer,      intent(in)  :: i,j,k,l
    integer,      intent(in)  :: i_nsp,j_nsp,k_nsp,l_nsp
    real(double), intent(in)  :: i_xyz(3),j_xyz(3),k_xyz(3),l_xyz(3)
    real(double), intent(out) :: eri_gto
    character(len=8), intent(out) :: i_nt, j_nt, k_nt, l_nt

    integer          :: i_nx, j_nx, k_nx, l_nx
    integer          :: i_ny, j_ny, k_ny, l_ny
    integer          :: i_nz, j_nz, k_nz, l_nz
    integer          :: i_n,  j_n,  k_n,  l_n
    
    integer          :: ia_gto, jb_gto, kg_gto, ld_gto
    real(double)     :: ai, aj, ak, al, di, dj, dk, dl 
    real(double)     :: i_norm, j_norm, k_norm, l_norm
    real(double)     :: ci, cj, ck, cl
    
    real(double), dimension(5) :: F

    i_nt = trim(gto( i_nsp )%sf( i )%nt)
    j_nt = trim(gto( j_nsp )%sf( j )%nt)
    k_nt = trim(gto( k_nsp )%sf( k )%nt)
    l_nt = trim(gto( l_nsp )%sf( l )%nt)

    spherical_ia: do i_n = 1,  gto( i_nsp )%sf( i )%sph_size
       spherical_jb: do j_n = 1,  gto( j_nsp )%sf( j )%sph_size
          spherical_kg: do k_n = 1,  gto( k_nsp )%sf( k )%sph_size
             spherical_ld: do l_n = 1,  gto( l_nsp )%sf( l )%sph_size

                i_nx = gto( i_nsp )%sf( i )%sph_nx( i_n )
                i_ny = gto( i_nsp )%sf( i )%sph_ny( i_n )
                i_nz = gto( i_nsp )%sf( i )%sph_nz( i_n )
                ci   = gto( i_nsp )%sf( i )%sph_c ( i_n )


                j_nx = gto( j_nsp )%sf( j )%sph_nx( j_n )
                j_ny = gto( j_nsp )%sf( j )%sph_ny( j_n ) 
                j_nz = gto( j_nsp )%sf( j )%sph_nz( j_n )
                cj   = gto( j_nsp )%sf( j )%sph_c ( j_n )

                k_nx = gto( k_nsp )%sf( k )%sph_nx( k_n )
                k_ny = gto( k_nsp )%sf( k )%sph_ny( k_n )                                        
                k_nz = gto( k_nsp )%sf( k )%sph_nz( k_n )
                ck   = gto( k_nsp )%sf( k )%sph_c ( k_n )

                l_nx = gto( l_nsp )%sf( l )%sph_nx( l_n )
                l_ny = gto( l_nsp )%sf( l )%sph_ny( l_n )
                l_nz = gto( l_nsp )%sf( l )%sph_nz( l_n )
                cl   = gto( l_nsp )%sf( l )%sph_c ( l_n )
                
                prim_ld: do ld_gto = 1, gto( l_nsp )%sf( l )%ngto
                   prim_kg: do kg_gto = 1, gto( k_nsp )%sf( k )%ngto
                      prim_jb: do jb_gto = 1, gto( j_nsp )%sf( j )%ngto
                         prim_ia: do ia_gto = 1, gto( i_nsp )%sf( i )%ngto

                            ai = gto( i_nsp )%sf( i )%a( ia_gto )
                            aj = gto( j_nsp )%sf( j )%a( jb_gto )
                            ak = gto( k_nsp )%sf( k )%a( kg_gto )
                            al = gto( l_nsp )%sf( l )%a( ld_gto )

                            di = gto( i_nsp )%sf( i )%d( ia_gto )
                            dj = gto( j_nsp )%sf( j )%d( jb_gto )
                            dk = gto( k_nsp )%sf( k )%d( kg_gto )
                            dl = gto( l_nsp )%sf( l )%d( ld_gto )

                            i_norm  = gto( i_nsp )%sf( i )%norm
                            j_norm  = gto( j_nsp )%sf( j )%norm
                            k_norm  = gto( k_nsp )%sf( k )%norm
                            l_norm  = gto( l_nsp )%sf( l )%norm

                            !i_norm = 1.0d0
                            !j_norm = 1.0d0
                            !k_norm = 1.0d0
                            !l_norm = 1.0d0
                            
                            eri_gto = eri_gto + (ci*cj*ck*cl) * (di*dj*dk*dl) * &
                                 eri_gto_hoh( &
                                 i_xyz(1), i_xyz(2), i_xyz(3), i_norm, i_nx, i_ny, i_nz, ai, &
                                 k_xyz(1), k_xyz(2), k_xyz(3), k_norm, k_nx, k_ny, k_nz, ak, &
                                 l_xyz(1), l_xyz(2), l_xyz(3), l_norm, l_nx, l_ny, l_nz, al, &
                                 j_xyz(1), j_xyz(2), j_xyz(3), j_norm, j_nx, j_ny, j_nz, aj)


                            !call compute_eri(1, ai, ia%xyz_ip, &
                            !     1, aj, jb%xyz_cv, &
                            !     1, ak, kg%xyz_cv, &
                            !     1, al, ld%xyz_cv)

                            !F = [0.41608906765397796d0,  0.044889937015574935d0, &
                            !     0.013706554295511562d0, 0.0063780699489852013d0, &
                            !     0.39523364424416996d0]
                            
                            !F(1) = 1.0d0

                            
                            !print*,  eri_gto_hoh( &
                            !     i_xyz(1), i_xyz(2), i_xyz(3), 1.0d0, i_nx, i_ny, i_nz, ai, &
                            !     k_xyz(1), k_xyz(2), k_xyz(3), 1.0d0, k_nx, k_ny, k_nz, ak, &
                            !     l_xyz(1), l_xyz(2), l_xyz(3), 1.0d0, l_nx, l_ny, l_nz, al, &
                            !     j_xyz(1), j_xyz(2), j_xyz(3), 1.0d0, j_nx, j_ny, j_nz, aj)


                            !CALL libint2_static_init()
                            !CALL libint2_init_eri(erieval, 0, C_NULL_PTR)
                            !CALL compute_eri_f(1, 0, &
                            !     0, [1.0d0], [ai], i_xyz, &
                            !     0, [1.0d0], [ak], k_xyz, &
                            !     0, [1.0d0], [al], l_xyz, &
                            !     0, [1.0d0], [aj], j_xyz, &
                            !     F, erieval)
                            !CALL print_eri(0, 0, 0, 0, 0, erieval)
                            !print*, erieval
                            !CALL libint2_cleanup_eri(erieval)

!!$                                  allocate(phi_on_grid_ia(2*extent+1,2*extent+1,2*extent+1))
!!$                                  phi_on_grid_ia = zero
!!$                                  allocate(phi_on_grid_jb(2*extent+1,2*extent+1,2*extent+1))
!!$                                  phi_on_grid_jb = zero
!!$                                  allocate(phi_on_grid_kg(2*extent+1,2*extent+1,2*extent+1))
!!$                                  phi_on_grid_kg = zero
!!$                                  allocate(phi_on_grid_ld(2*extent+1,2*extent+1,2*extent+1))
!!$                                  phi_on_grid_ld = zero
!!$                                  allocate(work_out_3d_(2*extent+1,2*extent+1,2*extent+1))
!!$                                  work_out_3d_ = zero
!!$                                  allocate(work_in_3d_(2*extent+1,2*extent+1,2*extent+1))
!!$                                  work_in_3d_ = zero
!!$
!!$
!!$                                  call exx_gto_on_grid_prim(inode,ia%ip, &
!!$                                       i_nsp, &
!!$                                       extent,  ia%xyz,       &
!!$                                       phi_on_grid_ia,r_int,xyz_zero,i_nx,i_ny,i_nz,ai)
!!$
!!$                                  call exx_gto_on_grid_prim(inode,jb%global_num,&
!!$                                       j_nsp,  &
!!$                                       extent, jb%xyz,&
!!$                                       phi_on_grid_jb,r_int,xyz_zero,j_nx,j_ny,j_nz,aj)
!!$
!!$                                  call exx_gto_on_grid_prim(inode,kg%global_num,&
!!$                                       k_nsp,&
!!$                                       extent,kg%xyz,&
!!$                                       phi_on_grid_kg,r_int,xyz_zero,k_nx,k_ny,k_nz,ak)
!!$
!!$                                  call exx_gto_on_grid_prim(inode,ld%global_num,&
!!$                                       l_nsp,&
!!$                                       extent,ld%xyz,&
!!$                                       phi_on_grid_ld,r_int,xyz_zero,l_nx,l_ny,l_nz,al)
!!$
!!$
!!$                                  work_in_3d_(:,:,:) = phi_on_grid_ld(:,:,:) * &
!!$                                       phi_on_grid_jb(:,:,:)
!!$
!!$                                  !call hf_v_on_grid(inode,extent,work_in_3d,&
!!$                                  !     work_out_3d,r_int, &
!!$                                  !     exx_psolver,p_scheme,pulay_radius,p_omega,&
!!$                                  !     p_ngauss,p_gauss, &
!!$                                  !     w_gauss)
!!$
!!$                                  call exx_v_on_grid(inode,extent,work_in_3d_,work_out_3d_,&
!!$                                       r_int,   &
!!$                                       exx_psolver,exx_pscheme,pulay_radius,p_omega,&
!!$                                       p_ngauss,p_gauss,&
!!$                                       w_gauss,fftwrho3d,reckernel_3d)
!!$
!!$                                  work_in_3d_ = phi_on_grid_ia(:,:,:) * phi_on_grid_kg(:,:,:)
!!$
!!$                                  !eri_pao = zero
!!$                                  do ii = 1, 2*extent + 1
!!$                                     do jj = 1, 2*extent + 1
!!$                                        do kk = 1, 2*extent + 1       
!!$                                           eri_pao = eri_pao + &
!!$                                                work_in_3d_ (ii,jj,kk) * &
!!$                                                work_out_3d_(ii,jj,kk) * dv 
!!$                                                !i_norm * j_norm * k_norm * l_norm
!!$                                        end do
!!$                                     end do
!!$                                  end do
!!$
!!$                                  !eri_pao = eri_pao 
!!$                                  
!!$                                  deallocate(phi_on_grid_ia)
!!$                                  deallocate(phi_on_grid_jb)
!!$                                  deallocate(phi_on_grid_kg) 
!!$                                  deallocate(phi_on_grid_ld) 
!!$                                  deallocate(work_out_3d_) 
!!$                                  deallocate(work_in_3d_) 



                         end do prim_ia
                      end do prim_jb
                   end do prim_kg
                end do prim_ld
                !

             end do spherical_ld
          end do spherical_kg
       end do spherical_jb
    end do spherical_ia

  end subroutine compute_eri_hoh
  !
end module exx_erigto
