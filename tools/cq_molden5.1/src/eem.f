      subroutine eem(iop,iasel,istat)
c
c     program for the calculation of eem parameters
c     patrick bultinck, jrf, 2001
c
c     J. Phys. Chem. A (2002)
c
      implicit double precision(a-h,o-z)
      parameter (mxeat=300)
      dimension var(mxeat,2),ipntr(mxeat)

      istat = 0
      numat = 0
      call valdis(var,ipntr,numat,iop,iasel,istat)
      if (istat.eq.0) call eemcalc(var,ipntr,numat)

   
      return

      end

      subroutine valdid(var,ipntr,numat,iop,iasel,istat,
     &                  ianz,iresid,qat)

      implicit double precision(a-h,o-z), integer (i-n)
      parameter (mxeat=300)
      parameter (maxtyp=19)
      parameter (mxtyp4=4*maxtyp)
      parameter (mxel=100)
      parameter (mexcl=10)

      common /eempar/param(maxtyp,2),pparam(maxtyp,2),
     &               peparm(10,2),esffprm(19,2),ipt(mxel,4)
      common /athlp/ iatoms, mxnat
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /metexc/ qexcl(mexcl),ianexc(mexcl)
      dimension ianz(*),iresid(*),qat(*),var(mxeat,2),ipntr(mxeat)

c eem NPA

      data ((param(i,k),k=1,2),i=1,maxtyp) /
     & -0.00298,0.89729,0.54541,0.00000,0.93882,0.00000,
     &  0.00000,0.10570,0.00000,0.61634,0.35970,0.33639,
     &  0.53938,0.37574,1.04219,0.72138,1.44195,1.57956,
     &  0.28351,0.00000,0.00000,0.00000,0.00000,0.00000,
     &  0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,
     &  0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,
     &  0.00000,0.00000/
c eem Mull + P en CL estimated
      data ((pparam(i,k),k=1,2),i=1,maxtyp) /
     &  0.20606,0.65971,0.00000,0.00000,0.00000,0.75331,
     &  0.00000,0.00000,0.00000,0.00000,0.36237,0.32966,
     &  0.49279,0.34519,0.73013,0.54428,0.72052,0.72664,
     &  0.78058,0.75467,0.00000,0.00000,0.00000,0.00000,
     &  0.00000,0.00000,0.00000,0.00000,0.36200,0.32000,
     &  0.00000,0.00000,0.41200,0.36600,0.00000,0.00000,
     &  0.00000,0.00000/
c
c Pesp, eem parameters
c        khi       eta      atom    MMtype                      Definition
c    0.68903   0.39888         1         H                 common hydrogen
c    0.71998   0.34395         6         C                   common carbon
c    0.97104   0.46763         7         N                 common nitrogen
c    0.95419   0.43623         8         O                   common oxigen
c    0.78155   0.24617        16         S                  common sulphur
c    0.67920   0.17446        12        Mg                common magnesium
c    0.98387   0.13229        15         P               common phosporous
c    0.86408   0.40514         9         F                 common fluorine
c    0.76957   0.17213        17        Cl                 common chlorine
c    0.77709   0.15474        35        Br                  common bromine
c    0.31627   0.90661         1       HAc                   acidic proton
c    0.58475   0.54282         1       Hpo         H bound to electrophile
c    0.61286   0.82238         1       Hpm    H bound to weak electrophile
c    0.67301   0.72481         1        HC        H bound to alkane carbon
c    0.68654   0.67974         1       HC3                       H in -CH3
c    0.54822   0.44969         6        CO                 carbonil carbon
c    0.79906   0.33085         6       CAr                 aromatic carbon
c    0.75858   0.28429         6      CArp       polarized aromatic carbon
c    0.72259   0.28054         6        CT                      sp3 carbon
c    0.84720   0.70094         6       CT3                  carbon in -CH3
c    0.78794   0.30423         6        CX                      sp2 carbon
c    0.85476   0.38483         6        CZ                       sp carbon
c    0.91123   0.40572         7        NT                    sp3 nitrogen
c    0.48744   0.43578         7        Np     quaterner positive nitrogen
c    0.80779   0.28151         7        NX                    sp2 nitrogen
c    0.98950   0.48874         7        NZ                     sp nitrogen
c    1.15452   0.50429         7       NAr               aromatic nitrogen
c    0.77898   0.35097         7      NArp      aromatic-positive nitrogen
c    0.71947   0.31954         7       NCO                  amide nitrogen
c    0.96737   0.43334         8        OH                       OH oxigen
c    0.65563   0.29674         8        Op                 positive oxigen
c    1.10775   0.41545         8        On                      oxidanione
c    0.81599   0.35269         8       OAr                 aromatic oxigen
c    0.82099   0.26674         8       Oox                      oxo-oxigen
c    0.84823   0.33950         8        OS                 aetheric oxigen
c    0.69661   0.25867         8      OAcH                acidic OH oxigen
c    0.77413   0.25308        16       SAc                sulphate sulphur
c    0.88299   0.41097        16        ST                     sp3 sulphur
c    0.84700   0.26207        16        SX                     sp2 sulphur
c
c currently only common types implemented HCNOF Mg P S Cl Br
c
      data ((peparm(i,k),k=1,2),i=1,10) /
     &  0.68903,0.39888,
     &  0.71998,0.34395,
     &  0.97104,0.46763,
     &  0.95419,0.43623,
     &  0.86408,0.40514,
     &  0.67920,0.17446,
     &  0.98387,0.13229,
     &  0.78155,0.24617,
     &  0.76957,0.17213,
     &  0.77709,0.15474/

c from ESFF parameterisation
c
c 1,2   H,H+
c 3-5   Csp3, Csp2, Csp
c 6-8   Nsp3, Nsp2, Nsp
c 9-11  Osp3, Osp2, Osp
c 12    F
c 13-14 Psp3, Psp2
c 15-16 Ssp3, Ssp2
c 17    Cl
c 18    Br
c 19    I
c
c hybridisation specific currently not used, only first atom 

      data ((esffprm(i,k),k=1,2),i=1,19) /
     &  0.233,0.440,0.252,0.440,
     &  0.2904,0.382,0.316,0.384,0.366,0.389,
     &  0.392,0.454,0.419,0.454,0.490,0.464,
     &  0.480,0.518,0.527,0.519,0.625,0.529,
     &  0.431,0.582,
     &  0.306,0.309,0.323,0.311,
     &  0.365,0.346,0.399,0.349,
     &  0.328,0.374,
     &  0.300,0.343,
     &  0.358,0.304/

c eem NPA supported elements HCNOF

      data (ipt(i,1),i=1,mxel) /1,0,
     1        0,0,0,6,7,8,9,0,
     2        0,0,0,0,0,0,0,0,
     3        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     4        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     5        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     3        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     4        0,0,0,0,0/

c eem Mull supported elements HCNOF + est. P Cl

      data (ipt(i,2),i=1,mxel) /1,0,
     1        0,0,0,6,7,8,9,0,
     2        0,0,0,0,15,0,17,0,
     3        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     4        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     5        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     3        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     4        0,0,0,0,0/

c eem PESP supported elements HCNOF Mg P S Cl Br

      data (ipt(i,3),i=1,mxel) /1,0,
     1        0,0,0,2,3,4,5,0,
     2        0,6,0,0,7,8,9,0,
     3        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,
     4        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     5        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     3        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     4        0,0,0,0,0/


c eem ESFF supported elements HCNOF  P S Cl Br I
c plus hybridisation specific, currently not used
 
      data (ipt(i,4),i=1,mxel) /1,0,
     1        0,0,0,3,6,9,12,0,
     2        0,0,0,0,13,15,17,0,
     3        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,18,0,
     4        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0,
     5        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     3        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     4        0,0,0,0,0/

      istat = 0

      icnt = 0
      ioverf = 0

      do i=1,iatoms
          if (iasel.gt.0.or.(iasel.lt.-3.and.iasel.eq.iresid(i)) )
     &    then

            iexcl = 0
            do k=1,mexcl
               if (ianz(i).eq.ianexc(k)) then
                  iexcl = 1
                  qat(i) = qexcl(k)
                  print*,' '
                  print*,'Excluded metal center ',ianz(i)
                  print*,'Set formal charge: ',qat(i)
                  print*,' '
               endif
            end do

            if (ianz(i).eq.100) iexcl = 1

            if (iexcl.eq.0) then
             icnt = icnt + 1
             if (icnt.lt.mxeat) then
                ipntr(icnt) = i
                do j=1,2
                   ipatm = ipt(ianz(i),iop)
                   if (ipatm.eq.0) then
                      call inferr(
     &                 'no parameters for element '//elemnt(ianz(i)),0)
                      ihasq = 0
                      istat = 1
                      return
                   endif
                   if (iop.eq.1) then
c NPA
                      var(icnt,j) = param(ipatm,j)
                   elseif (iop.eq.2) then
c Mull
                      var(icnt,j) = pparam(ipatm,j)
                   elseif (iop.eq.3) then
c PESP
                      var(icnt,j) = peparm(ipatm,j)
                   elseif (iop.eq.4) then
c ESFF
                      var(icnt,j) = esffprm(ipatm,j)
                   endif
                end do
c                print*,i,' ',ianz(i),(var(icnt,j),j=1,2)
             else
                icnt = mxeat
                ioverf = 1
             endif

            endif
          endif
      end do

      if (ioverf.eq.1) then
          call inferr(
     &     'increase parameter mxeat in subroutine valdis',0)
      endif

      numat = icnt

      return
      end

      subroutine eemcald(var,ipntr,numat,ianz,coo,qat)
c
      implicit double precision(a-h,o-z), integer (i-n)
      parameter (mxeat=300)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /totchg/ itot
      common /athlp/ iatoms, mxnat
      dimension det(2)
      dimension coo(3,*),qat(*),ianz(*)
      dimension var(mxeat,2),ipntr(mxeat),work(mxeat+1,mxeat+1)
      dimension x(mxeat+1,mxeat+1),y(mxeat+1), ipvt(mxeat+1),q(mxeat)

      do i=1,(numat-1)
         x(i,i)       = 2*var(i,2)
         x(i,numat+1) = -1
         x(numat+1,i) = 1

         do j=i+1,numat
            x(i,j) = ( coo(1,ipntr(i)) - coo(1,ipntr(j)) )**2 +
     &               ( coo(2,ipntr(i)) - coo(2,ipntr(j)) )**2 +
     &               ( coo(3,ipntr(i)) - coo(3,ipntr(j)) )**2 
            x(i,j) = 1.0d0/dsqrt(x(i,j))
            x(j,i) = x(i,j)
         end do

      end do
      x(numat,numat)     = 2*var(numat,2)
      x(numat,numat+1)   = -1
      x(numat+1,numat)   = 1
      x(numat+1,numat+1) = 0

      y(numat+1)       = dble(itot)

      do i=1,numat
         y(i) = -var(i,1)
      end do

      call dgefa(x,mxeat+1,numat+1,ipvt,info)
      call dgedi(x,mxeat+1,(numat+1),ipvt,det,work,01)

      call mtmul(x,y,q,numat)

      write (*,'(a)') ' '
      write (*,'(a)') 'EEM charges'
      write (*,'(a)') ' '
      qtot = 0.0d0
      do i=1,(numat)
          j = ipntr(i)
          qat(j) = q(i)
          qtot = qtot + q(i)
          write (*,'(i5,1x,i3,1x,f10.5)') j,ianz(j),qat(j)
      end do
      write (*,'(a)') ' '
      write (*,'(a,f10.3)') 'Sum of EEM charges = ',qtot
      ihasq = 1

      return
      end

c
      subroutine mtmul(a,y,b,numat)
c
      parameter (mxeat=300)
c
      integer i,j,numat
      double precision a(mxeat+1,mxeat+1),y(mxeat+1),b(mxeat+1)
c
      do i=1,numat+1

         b(i) = 0.0

         do j=1,numat+1
            b(i) = b(i) + a(i,j)*y(j)
         end do

      end do
c
      return
      end
c
      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      implicit double precision(a-h,o-z)
      dimension a(lda,*),det(2),work(*),ipvt(*)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. abs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran abs,mod
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d+00
         det(2) = 0.0d+00
         ten = 10.0d+00
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d+00) go to 60
   10       if (abs(det(1)) .ge. 1.0d+00) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d+00
            go to 10
   20       continue
   30       if (abs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d+00
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) return
         do k = 1, n
            a(k,k) = 1.0d+00/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d+00
               call daxpy(k,t,a(1,k),1,a(1,j),1)
            end do
   90       continue
         end do
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) return
         do kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d+00
            end do
            do j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
            end do
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
         end do

      return
      end

      subroutine dgefa(a,lda,n,ipvt,info)
      implicit double precision(a-h,o-z)
      dimension a(lda,*),ipvt(*)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1

      if (nm1.lt.1) goto 70

      do k=1,nm1
         kp1 = k + 1

c        find l = pivot index

         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l

c        zero pivot implies this column already triangularized

         if (a(l,k).eq.0.0d0) then
            info = k
         else

c           interchange if necessary

            if (l.ne.k) then
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
            endif

c           compute multipliers

            t = -1.0d+00/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)

c           row elimination with column indexing

            do j = kp1, n

               t = a(l,j)

               if (l.ne.k) then
                  a(l,j) = a(k,j)
                  a(k,j) = t
               endif

               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
            end do

         endif

      end do

   70 continue

      ipvt(n) = n
      if (a(n,n) .eq. 0.0d+00) info = n

      return
      end

      double precision function dasum(n,dx,incx)

c     takes the sum of the absolute values.

      double precision dx(1),dtemp
      integer i,incx,m,mp1,n,nincx

      dasum = 0.0d+00
      dtemp = 0.0d+00

      if (n.le.0) return

      if (incx.ne.1) then

         nincx = n*incx

         do i=1,nincx,incx
           dtemp = dtemp + abs(dx(i))
         end do

         dasum = dtemp

      else

         m = mod(n,6)

         if (m.ne.0) then

            do i = 1,m
              dtemp = dtemp + abs(dx(i))
            end do

            if (n.lt.6) then
               dasum = dtemp
               return
            endif

         endif

         mp1 = m + 1

         do i = mp1,n,6
           dtemp = dtemp + abs(dx(i)) + abs(dx(i + 1)) + abs(dx(i + 2))
     &     + abs(dx(i + 3)) + abs(dx(i + 4)) + abs(dx(i + 5))
         end do

         dasum = dtemp

      endif

      return
      end

      subroutine daxpy(n,da,dx,incx,dy,incy)
      implicit double precision(a-h,o-z)
      dimension dx(1),dy(1)
c
c     constant times a vector plus a vector.
c           dy(i) = dy(i) + da * dx(i)
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      if (n.le.0) return
      if (da .eq. 0.0d0) return

      if (incx.eq.1.and.incy.eq.1) then
         m = mod(n,4)
         if (m.ne.0) then

            do i = 1,m
              dy(i) = dy(i) + da*dx(i)
            end do

            if (n.lt.4) return
         endif

         mp1 = m + 1

         do i = mp1,n,4
           dy(i) = dy(i) + da*dx(i)
           dy(i + 1) = dy(i + 1) + da*dx(i + 1)
           dy(i + 2) = dy(i + 2) + da*dx(i + 2)
           dy(i + 3) = dy(i + 3) + da*dx(i + 3)
         end do

      else

c        code for unequal increments or equal increments
c        not equal to 1

         ix = 1
         iy = 1

         if (incx.lt.0) ix = (-n+1)*incx + 1
         if (incy.lt.0) iy = (-n+1)*incy + 1

         do i = 1,n
           dy(iy) = dy(iy) + da*dx(ix)
           ix = ix + incx
           iy = iy + incy
         end do

      endif

      return
      end

      subroutine dcopy(n,dx,incx,dy,incy)
      implicit double precision(a-h,o-z)
      dimension dx(*),dy(*)
c
c     copies a vector.
c           dy(i) <== dx(i)
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end

      double precision function ddot(n,dx,incx,dy,incy)
      implicit double precision(a-h,o-z)
      dimension dx(1),dy(1)
c
c     forms the dot product of two vectors.
c           dot = dx(i) * dy(i)
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      ddot = 0.0d+00
      dtemp = 0.0d+00
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end

      double precision function dnrm2 ( n, dx, incx)
      integer          next
      double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d+00, 1.0d+00/
c
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  sqrt(u/eps)  over all known machines.
c         cuthi = minimum of  sqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d+19
c     data cutlo, cuthi / 8.232d-11,  1.304d+19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d+19 /
c
      j=0
      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( abs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( dx(i) .eq. zero) go to 200
      if( abs(dx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = abs(dx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( abs(dx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( abs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = abs(dx(i))
         go to 200
c
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/n
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j =i,nn,incx
      if(abs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = sqrt( sum )
      go to 300
c
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      dnrm2 = xmax * sqrt(sum)
  300 continue
      return
      end

      subroutine drot(n,dx,incx,dy,incy,c,s)
      implicit double precision(a-h,o-z)
      dimension dx(1),dy(1)
c
c     applies a plane rotation.
c           dx(i) =  c*dx(i) + s*dy(i)
c           dy(i) = -s*dx(i) + c*dy(i)
c     jack dongarra, linpack, 3/11/78.
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
      end do

      return
c
c       code for both increments equal to 1
c
   20 do i = 1,n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
      end do

      return
      end

      subroutine drotg(da,db,c,s)
c
c     construct givens plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
      double precision da,db,c,s,roe,scale,r,z
      double precision zero, one
c
      parameter (zero=0.0d+00, one=1.0d+00)
c
c-----------------------------------------------------------------------
c
c
      roe = db
      if( abs(da) .gt. abs(db) ) roe = da
      scale = abs(da) + abs(db)
      if ( scale .eq. zero ) then
         c = one
         s = zero
         r = zero
      else
         r = scale*sqrt((da/scale)**2 + (db/scale)**2)
         r = sign(one,roe)*r
         c = da/r
         s = db/r
      endif
      z = one
      if( abs(da) .gt. abs(db) ) z = s
      if( abs(db) .ge. abs(da) .and. c .ne. zero ) z = one/c
      da = r
      db = z

      return
      end

      subroutine dscal(n,da,dx,incx)
      implicit double precision(a-h,o-z)
      dimension dx(1)
c
c     scales a vector by a constant.
c           dx(i) = da * dx(i)
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx

      do i = 1,nincx,incx
        dx(i) = da*dx(i)
      end do

      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if ( m .eq. 0 ) go to 40
      do i = 1,m
        dx(i) = da*dx(i)
      end do
      if( n .lt. 5 ) return

   40 mp1 = m + 1
      do i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
      end do

      return
      end

      subroutine dswap(n,dx,incx,dy,incy)
      implicit double precision(a-h,o-z)
      dimension dx(1),dy(1)
c
c     interchanges two vectors.
c           dx(i) <==> dy(i)
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
      end do

      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end

      integer function idamax(n,dx,incx)
      implicit double precision(a-h,o-z)
      dimension dx(1)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      rmax = abs(dx(1))
      ix = ix + incx

      do i = 2,n
         if (abs(dx(ix)).gt.rmax) then
            idamax = i
            rmax = abs(dx(ix))
         endif
         ix = ix + incx
      end do

      return
c
c        code for increment equal to 1
c
   20 rmax = abs(dx(1))
      do i = 2,n
         if(abs(dx(i)).gt.rmax) then
            idamax = i
            rmax = abs(dx(i))
         endif
      end do

      return
      end

      subroutine dgemv(forma,m,n,alpha,a,lda,x,incx,beta,y,incy)
      implicit double precision(a-h,o-z)
      character*1 forma
      dimension a(lda,*),x(*),y(*)
      parameter (zero=0.0d+00, one=1.0d+00)
c
c        clone of -dgemv- written by mike schmidt
c
      locy = 1
      if(forma.eq.'t') go to 200
c
c                  y = alpha * a * x + beta * y
c
      if (alpha.eq.one  .and.  beta.eq.zero) then
         do i=1,m
            y(locy) =       ddot(n,a(i,1),lda,x,incx)
            locy = locy+incy
         end do
      else
         do i=1,m
            y(locy) = alpha*ddot(n,a(i,1),lda,x,incx) + beta*y(locy)
            locy = locy+incy
         end do
      end if
      return
c
c                  y = alpha * a-transpose * x + beta * y
c
  200 continue

      if(alpha.eq.one  .and.  beta.eq.zero) then
         do i=1,n
            y(locy) =       ddot(m,a(1,i),1,x,incx)
            locy = locy+incy
         end do
      else
         do i=1,n
            y(locy) = alpha*ddot(m,a(1,i),1,x,incx) + beta*y(locy)
            locy = locy+incy
         end do
      end if

      return
      end

      subroutine filnam(filn,nnn)
      character name*70,filn*(*)

      i1 = 0
      do i = 1,nnn
         if (filn(i:i).ne.' ') then
            i1 = i1 + 1
            name(i1:i1) = filn(i:i)
         endif
      end do

      do i = 1,nnn
         filn(i:i) = ' '
      end do

      do i = 1,i1
         filn(i:i) = name(i:i)
      end do

      return
      end

      subroutine distchd(coo)
c
c
      implicit double precision(a-h,o-z)
      common /athlp/ iatoms, mxnat
      dimension coo(3,*)

      tol = 0.1d0
      do i=1,iatoms
         do j=i+1,iatoms
             dd = 0.0d0
             do k=1,3
                d1 = coo(k,i) - coo(k,j)
                d2 = d1*d1
                dd = dd + d2
             end do 
             if (dd.lt.tol) print*,'close ',i,' ',j
         end do
      end do

      return
      end

      subroutine sigini
      implicit double precision(a-h,o-z)
      parameter (mxsigm=16)
      parameter (mxmol2=41)
      common /sigma/ siga(mxsigm),sigb(mxsigm),sigc(mxsigm),
     &               sigd(mxsigm),impmol2(mxmol2)
c      data sigs  /"H ","C3","C2","C1","N3","N2","N1","O3","O2","F ",
c     &       "Cl","Br","I ","S3","P ","DU"/     
      data siga /7.17,7.98,8.79,10.39,11.54,12.87,15.68,14.18,
     &       17.07,14.66, 11.00, 10.08,9.90,10.14,8.9,0.0/
      data sigb /6.24,9.18,9.32,9.45,10.82,11.15,11.70,12.92,
     &           13.79,13.85,9.69,8.47,7.96,9.13,8.24,0.0/
      data sigc /-0.56,1.88,1.51,0.73,1.36,0.85,-0.27,1.39,
     &           0.47,2.31,1.35,1.16,0.96,1.38,0.96,0.0/
      data impmol2 /16,16,16,16,2,3,4,3,3,5,6,7,6,6,6,5,
     &           8,9,9,8,8,14,14,14,14,15,1,1,1,10,11,12,13,
     &           16,16,16,16,16,16,16,16/

c sybyl atom types to gasteiger types
c1 any
c2 hal
c3 het
c4 hev
c5 C.3     2
c6 C.2     3
c7 C.1     4
c8 C.ar    3
c9 C.cat   ?
c10 N.3    5
c11 N.2    6
c12 N.1    7
c13 N.ar   6
c14 N.am   6
c15 N.pl3  6
c16 N.4    5
c17 O.3    8
c18 O.2    9
c19 O.co2  9
c20 O.spc  8
c21 O.t3p  8
c22 S.3    14
c23 S.2    14
c24 S.O    14
c25 S.O2   14
c26 P.3    15
c27 H      1
c28 H.spc  1
c29 H.t3p  1
c30 F      10
c31 Cl     11
c32 Br     12
c33 I      13
c34 Si
c35 Lp
c36 Du
c37 Na
c38 K
c39 Ca
c40 Li
c41 Al

c unknown is DU 16

      do i=1,mxsigm
         sigd(i) = siga(i) + sigb(i) + sigc(i)
      end do
      sigd(15) = 1.0d0

      return
      end
                              
      subroutine inigad(qtot,qat,ianz,iconn,ityp)
      implicit double precision(a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (mxsigm=16)
      parameter (mxmol2=41)
      common /sigma/ siga(mxsigm),sigb(mxsigm),sigc(mxsigm),
     &               sigd(mxsigm),impmol2(mxmol2)
      common /athlp/ iatoms, mxnat
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      integer*2 ityp
      common /types/ iff
      dimension qat(*),ianz(*),iconn(mxcon+1,*),ityp(*)

      ihasq = 0
      iff = 5
      call dotyp(0)
      qtot = 0.0d0
  
      do i=1,iatoms
         qat(i) = 0.0d0
         ia = ianz(i)

         if (ia.eq.6) then
           if (ityp(i).eq.9) qat(i) = 1.0d0
         elseif (ia.eq.7) then
           ibnd = 0
           do j=1,iconn(1,i)
              if (iconn(1+j,i).gt.0) ibnd = ibnd + 1
           end do
           if (ibnd.eq.4) qat(i) = 1.0d0
         elseif (ia.eq.8) then
c carboxyl
           if (ityp(i).eq.19) qat(i) = -0.5d0
         endif
         qtot = qtot + qat(i)
      end do

      return
      end

      subroutine clqgas(ishoh)
      implicit double precision(a-h,o-z)
      logical dozme
      common /getpnt/irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &               iconv,ircus,dozme
      parameter (mxheta=150)
      character*3 hetz
      common /clfstr/ ihashz,ihetq(mxheta),ihqset(mxheta),ihhadd(mxheta)
     &                ,hetz(mxheta)

      istat = 0

      if (ipdbon.eq.1) then

         call numhet(nhmol)
         do i=4,nhmol
            if (i.ne.ishoh) then
               print*,' '
               if (ihashz.eq.1) then
                  print*,'HETATM residue ',hetz(i+1)
               else
                  print*,'HETATM residue ',(i+1-4)
               endif
               print*,' '
               call calgas(-i,istat)
            endif
         end do

        print*,' '
        print*,'WARNING: Total charges of the different HETATM residues'
        print*,'         assumed zero. If this is NOT correct, assign '
        print*,'         HETATM charges individually, by clicking with'
        print*,'         2nd mouse button on the residue and select   '
        print*,'         Calculate Charges'
        print*,' '

      else
         call calgas(1,istat)
      endif

      if (istat.eq.1) call messg(14)

      return
      end

      subroutine calgad(iasel,istat,qat,ianz,iconn,iresid,ityp)
      implicit double precision(a-h,o-z)
      parameter (numat1=20000)
      parameter (mxcon=10)
      parameter (mxsigm=16)
      parameter (mxmol2=41)
      common /sigma/ siga(mxsigm),sigb(mxsigm),sigc(mxsigm),
     &               sigd(mxsigm),impmol2(mxmol2)
      common /athlp/ iatoms, mxnat
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      integer*2 ityp
      common /types/ iff
      dimension zz(numat1)
      dimension qat(*),ianz(*),iconn(mxcon+1,*),ityp(*),iresid(*)

      call sigini
  
      ihasq = 0
      istat = 0
      iff = 5
      call dotyp(0)
  
      do i=1,iatoms
         if (iasel.gt.0.or.(iasel.lt.-3.and.iasel.eq.iresid(i)) )
     &   then

            it = ityp(i)

            if (it.gt.0.and.it.le.mxmol2) then
                l = impmol2(it)
            else
                print*,'Gasteiger: Untyped mol2 atom type.'
                istat = 1
                return
            endif

            if (i.le.numat1) then
                zz(i) = siga(l)
            else
                print*,'array for gasteiger calculation to small'
                istat = 1
                return
            endif

            qat(i) = 0.0d0
            ia = ianz(i)

            if (ia.eq.6) then

              if (ityp(i).eq.9) qat(i) = 1.0d0

            elseif (ia.eq.7) then

              ibnd = 0
              do j=1,iconn(1,i)
                 if (iconn(1+j,i).gt.0) ibnd = ibnd + 1
              end do
              if (ibnd.eq.4) qat(i) = 1.0d0

c              if (N3+) qat(i) = 1.0d0

            elseif (ia.eq.8) then
c carboxyl
              if (ityp(i).eq.19) qat(i) = -0.5d0

            endif
         endif
      end do

      fac = 1.0d0
      icnt = 0

      do while (.true.)

          fac = fac*0.5d0
          sd1 = 0.0d0

          do i=1,iatoms

           if (iasel.gt.0.or.(iasel.lt.-3.and.iasel.eq.iresid(i)) )
     &     then

             l = impmol2(ityp(i))

             if (sigd(l).ne.1.0d0) then

                qt = qat(i)
                do j=1,iconn(1,i)
                   if (iconn(1+j,i).gt.0) then
                      jj = iabs(iconn(1+j,i))
                      ll = impmol2(ityp(jj))
                      if (sigd(ll).ne.1.0d0) then
                         sd2 = sigd(ll)
                         if (i.le.numat1.and.jj.le.numat1) then
                            if (zz(jj).gt.zz(i)) sd2 = sigd(l)
                         endif
                         if (ianz(i).eq.1.or.ianz(jj).eq.1) 
     &                      sd2 = 20.02d0
                         if (i.le.numat1.and.jj.le.numat1) then
                            qat(i) = qat(i) + (zz(jj) - zz(i))*fac/sd2
                         endif
                      endif
                   endif
                end do

                qt = dabs(qat(i) - qt)
                if (qt.gt.sd1) sd1 = qt

             endif

           endif
          end do

          if (sd1.ge.1.0d-3) then

             do i=1,iatoms

              if (iasel.gt.0.or.(iasel.lt.-3.and.iasel.eq.iresid(i)) )
     &        then
                l = impmol2(ityp(i))
                zz(i) = siga(l) + sigb(l)*qat(i) + sigc(l)*qat(i)*qat(i)
              endif

             end do

          endif 

          icnt = icnt + 1

          if (.not.(sd1.gt.1.0d-3.and.icnt.le.5)) goto 100
      end do

100   continue

      write (*,'(a)') ' '
      write (*,'(a)') 'Gasteiger charges'
      write (*,'(a)') ' '
      qtot = 0.0d0
      do i=1,iatoms
          if (iasel.gt.0.or.(iasel.lt.-3.and.iasel.eq.iresid(i)) )
     &    then
             qtot = qtot + qat(i)
             write (*,'(i5,1x,i3,1x,f10.5)') i,ianz(i),qat(i)
          endif
      end do
      write (*,'(a)') ' '
      write (*,'(a,f10.3)') 'Sum of Gasteiger charges = ',qtot

      ihasq = 1
      return
      end

      subroutine clqeem(iop,ishoh)
      implicit double precision(a-h,o-z)
      logical dozme
      common /getpnt/irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &               iconv,ircus,dozme
      parameter (mxheta=150)
      character*3 hetz
      common /clfstr/ ihashz,ihetq(mxheta),ihqset(mxheta),ihhadd(mxheta)
     &                ,hetz(mxheta)

      istat = 0

      if (ipdbon.eq.1) then

         call numhet(nhmol)
         do i=4,nhmol
            if (i.ne.ishoh) then
               print*,' '
               if (ihashz.eq.1) then
                  print*,'HETATM residue ',hetz(i+1)
               else
                  print*,'HETATM residue ',(i+1-4)
               endif
               print*,' '
               call eem(iop,-i,istat)
            endif
         end do

        print*,' '
        print*,'WARNING: Total charges of the different HETATM residues'
        print*,'         assumed zero. If this is NOT correct, assign '
        print*,'         HETATM charges individually, by clicking with'
        print*,'         2nd mouse button on the residue and select   '
        print*,'         Calculate Charges'
        print*,' '

      else
         call eem(iop,1,istat)
      endif

      if (istat.eq.1) call messg(14)

      return
      end

      subroutine clceem(ishoh)
      implicit double precision(a-h,o-z)
      logical dozme
      common /getpnt/irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &               iconv,ircus,dozme
      common /types/ iff
      istat = 0

      if (iff.eq.7) then
         if (ipdbon.eq.1) then
            call numhet(nhmol)
            do i=4,nhmol
               if (i.ne.ishoh) then
                  call eem(3,-i,istat)
               endif
            end do
         else
            call eem(3,1,istat)
         endif
      endif

      print*,' '
      print*,'EEM: Electronegativity Equalization Method:'
      print*,' '
      print*,' Patrick Bultinck et al'
      print*,' J. Phys. Chem. A (2002)'
      print*,' '

      if (istat.eq.1) call messg(14)

      return
      end

