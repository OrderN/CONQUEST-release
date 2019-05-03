      subroutine mpolefit(idip,nesp,esp,connl,dx,dy,dz,iz,dmachg,
     &                    keywrd,ichadd)

c Routine for fitting of atomic monopoles, dipoles and quadrupoles to
c the electrostatic potential. Multipoles are obtained in spherical 
c tensor form, in coordinate system of the molecule, 
c and are given in atomic units.
c
c       Q_00
c       Q_10  Q_1c  Q_1s
c       Q_20  Q_21c Q_21s Q_22c Q_22s
c
c Fitting strategy is based on:
c Hinsen and Roux, J. Comp. Chem. 1997, 18, 368
c using Singular Value Decomposition 
c constraints are imposed by elimination
c
c Wijnand Mooij, may 1998

c
c     process keywords ATPOL MPTOL MPEQUIV
c
c ATPOL sets the use of charge, dipole and quadrupole on each atom
c Default is M+D+Q, except for hydrogen M+D
c
c     ATPOL=(atnum1/icode,atnum2-atnum10/icode,...)
c
c     icode = +1 (if monopole) +2 (if dipole) +4 (if quadrupole)
c     
c  SV tolerance (default 1D-5)
c
c     MPTOL=value
c
c  Equivalences can be set for groups of atoms.
c  Note that equivalences are set in variable numbers, not in atom number;
c  these will differ when dipoles, quadrupoles are used.
c  examp: if atom 1=M+D+Q, then charge on atom 2 has variable number 10
c
c     MPEQUIV=(var1,var2,var3/var4,var5/...)
c

      implicit double precision (a-h,o-z)
      logical dmachg
      dimension esp(*), connl(3,*)
      parameter (numatm=2000)
      parameter (mxel=100)
      parameter (mesp = 30000)
      parameter (np = 200)
      parameter (tolc=1d-10)
      parameter (tol=1d-5)
      parameter (lfmax=2)
      parameter (lfmaxf=(lfmax+1)*(lfmax+1))     
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line
      common /curlin/ line
      character*(*) keywrd
      character*137 keyhlp,tstr,strt
      logical keyr
      dimension u(mesp,np),v(np,np),w(np)
      dimension b(np,np),c(np),p(np,np),binvc(np),esp1(mesp)
      dimension scratch(np)
      dimension pvec(3),temp(5)
      dimension ict(50)
      logical monopole,dipole,quadrupole
      common /mpfit/ q(lfmaxf,numatm),
     &     monopole(numatm),dipole(numatm),quadrupole(numatm)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon

      write(iun3,*) ' ---------------- '
      write(iun3,*) '     MPOLEFIT'   
      write(iun3,*) ' ---------------- '

      do i=1,numatm
         monopole(i)   = .true.
         dipole(i)     = .false.
         quadrupole(i) = .false.
      end do

      do i=natoms+1,numatm
         nat(i) = 99
      end do

      elemnt(99) = ' X'

      call setmp(keywrd)

c Check for user-input SV tolerance

      if (.not.keyr(keywrd,'MPTOL',toll)) then
         toll = tol
      endif
      write(iun3,*) 'MPTOL =',toll

      rt3 = dsqrt(3.0d0)
c
c     conversion factor for debye to atomic units
c
      cf = 5.2917715d-11*1.601917d-19 / 3.33564d-30
      toang = 0.52917706d0
     
c     calculate total charge on the molecule

      iz = 0
      do i=1,natoms
          iz = iz + nat(i)
      end do
      iz = iz - nelecs
      if (iz.ne.0) write(iun3,*) 'Charge of molecule = ',iz

      if (ichadd.ne.0) then

         call chadd(natoms)

         do j=1,natoms
            write(iun3,'(i4,3f10.5,i4,a)') j,(xyz(m,j)*toang,m=1,3),
     &           nat(j),elemnt(nat(j))
         end do

      endif


c Set the number of variables

      nvar = 0
      do i=1,natoms
         if (monopole(i)) nvar = nvar + 1
         if (dipole(i)) nvar = nvar + 3
         if (quadrupole(i)) nvar = nvar + 5
      end do
      
      if (nvar.gt.np) then
         write(iun3,*) 'Too many parameters in mpolefit'
         stop
      endif
      if (nesp.gt.mesp) then
         write(iun3,*) 'Too many ESP points in mpolefit'
         stop
      endif

c Set up the constraint relation b*q=c

      ivar = 0
      do j=1,natoms

         if (monopole(j)) then
            ivar = ivar + 1
            b(1,ivar) = 1.d0
            if (idip.eq.1) then
               b(2,ivar) = xyz(1,j)
               b(3,ivar) = xyz(2,j)
               b(4,ivar) = xyz(3,j)
            endif
         endif

         if (dipole(j)) then
            do k=1,3
               b(1,ivar+k) = 0.d0
            end do
            if (idip.eq.1) then

c              X-component  of dipole

               b(2,ivar+1) = 0.d0
               b(2,ivar+2) = 1.d0
               b(2,ivar+3) = 0.d0

c              Y-component  of dipole

               b(3,ivar+1) = 0.d0
               b(3,ivar+2) = 0.d0
               b(3,ivar+3) = 1.d0

c              Z-component  of dipole

               b(4,ivar+1) = 1.d0
               b(4,ivar+2) = 0.d0
               b(4,ivar+3) = 0.d0
            endif
            ivar = ivar + 3
         endif

         if (quadrupole(j)) then
            do k=1,5
               b(1,ivar+k) = 0.d0
            end do
            if (idip.eq.1) then
               do k=1,5
                  do i=2,4
                     b(i,ivar+k) = 0.d0
                  end do
               end do
            endif
            ivar = ivar + 5
         endif

      end do
      
      c(1) = dfloat(iz)
      ncon = 1

      if (idip.eq.1) then
         c(2) = dx/cf
         c(3) = dy/cf
         c(4) = dz/cf
         ncon = 4
      endif

c     Augment the constraint matrix to square

      do i=ncon+1,nvar
         do j=1,nvar      
            b(i,j) = 0.d0
         end do
         c(i) = 0.d0
      end do


c Set equivalences in constraint matrix

      i = index(keywrd,'MPEQUIV')
      if (i.ne.0) then
         i1 = index(keywrd(i+1:),'(')
         if (i1.ne.0) then
           i2 = index(keywrd(i+1:),')')
           if (i2.ne.0) then
              keyhlp = keywrd(i+i1+1:i+i2-1)
              l = i2 - i1 - 1
              call spatrm(keyhlp,l)
              call setlin(keyhlp,47)
              keyhlp = line(1:l)
              do while (nxtwrd(tstr,nstr,itype,rtype).eq.1)
                  keyhlp = line(1:l)
                  call setlin(tstr,44)
                  n = 0
                  do while (nxtwrd(strt,nstr,itype,rtype).eq.2)
                      n = n + 1
                      ict(n) = itype
                  end do
                  do k=2,n
                     ncon = ncon + 1
                     b(ncon,ict(1)) = 1.d0
                     b(ncon,ict(k)) = -1.d0
c                     print*,'b(',ncon,',',ict(1),') = 1.0d0',
c  &                         ' b(',ncon,',',ict(k),') = -1.0d0'
                  end do
                  call setlin(keyhlp,0)
              end do
           endif
         endif
      endif

 13   format(6f9.4)

      write(iun3,*) 'Number of variables :       ',nvar
      write(iun3,*) 'Number of constraints :     ',ncon

      if (ncon.gt.nvar) then
         write(iun3,*) 'More constraints than variables in mpolefit'
         stop
      endif


c Copy b to u. B is needed later on.

      do i=1,nvar
         do j=1,nvar
            u(j,i) = b(j,i)
         end do
      end do

      call svd(mesp,nvar,nvar,u,w,.true.,u,.true.,v,ierr,scratch,np)

      wmax = 0.d0
      do j=1,nvar
         if (w(j).gt.wmax) wmax = w(j)
      end do     
       

      thresh = tolc*wmax
      if (thresh.lt.1d-10) thresh = 1d-10

      nsvd = 0
      do j=1,nvar
         if (w(j).lt.thresh) then
              nsvd = nsvd + 1
              w(j) = 0.d0
         endif
      end do
      write(iun3,*) 'Rank of constraint matrix : ',nvar-nsvd
    
c  Calculate  1/w(i)*Utranspose, in p            

      do i=1,nvar
         if (w(i).ne.0.d0) then
            do j=1,nvar
               p(i,j) = u(j,i) / w(i)
            end do
         else
            do j=1,nvar
               p(i,j) = 0.d0
            end do
         endif
      end do
      
c Inverse of B, in u

      do i=1,nvar
         do j=1,nvar
            u(i,j) = 0.d0
            do k=1,nvar
               u(i,j) = u(i,j) + v(i,k)*p(k,j)
            end do
         end do
      end do

c     Calculate p=1-binv*b

      do i=1,nvar
         do j=1,nvar
            p(i,j) = 0.d0
            do k=1,nvar
               p(i,j) = p(i,j) - u(i,k)*b(k,j)
            end do
         end do
         p(i,i) = p(i,i) + 1.d0
      end do

c Calculate binv*c

      do i=1,nvar
         binvc(i) = 0.0d0
         do j=1,nvar
            binvc(i) = binvc(i) + u(i,j)*c(j)
         end do
      end do


      do i=1,nesp
         do j=1,nvar
            u(i,j) = 0.d0
         end do
         esp1(i) = esp(i)
      end do

c     Set up design matrix A`=AP en p`=p-A*Binv*C

      do i=1,nesp
         ivar = 0
         do j=1,natoms

            if (monopole(j)) then
               ivar = ivar + 1
               pvec(1) = connl(1,i) - xyz(1,j)
               pvec(2) = connl(2,i) - xyz(2,j)
               pvec(3) = connl(3,i) - xyz(3,j)
               qp1 = pvec(1)*pvec(1)
               qp2 = pvec(2)*pvec(2)
               qp3 = pvec(3)*pvec(3)
               r2 = qp1 + qp2 + qp3
               r = dsqrt(r2)
               temp1 = 1.d0 / r
               do k=1,nvar
                    u(i,k) = u(i,k) + temp1*p(ivar,k)
               end do
               esp1(i) = esp1(i) - temp1*binvc(ivar)
            endif

            if (dipole(j)) then
               pvec(1) = connl(1,i) - xyz(1,j)
               pvec(2) = connl(2,i) - xyz(2,j)
               pvec(3) = connl(3,i) - xyz(3,j)
               qp1 = pvec(1)*pvec(1)
               qp2 = pvec(2)*pvec(2)
               qp3 = pvec(3)*pvec(3)
               r2 = qp1 + qp2 + qp3
               r = dsqrt(r2)
               r3= r2*r
               temp(1) = pvec(3) / r3
               temp(2) = pvec(1) / r3
               temp(3) = pvec(2) / r3
               do kk=1,3
                  do k=1,nvar
                     u(i,k) = u(i,k) + temp(kk)*p(ivar+kk,k)
                  end do
                  esp1(i) = esp1(i) - temp(kk)*binvc(ivar+kk)
               end do
               ivar = ivar + 3
            endif

            if (quadrupole(j)) then
               pvec(1) = connl(1,i) - xyz(1,j)
               pvec(2) = connl(2,i) - xyz(2,j)
               pvec(3) = connl(3,i) - xyz(3,j)
               qp1 = pvec(1)*pvec(1)
               qp2 = pvec(2)*pvec(2)
               qp3 = pvec(3)*pvec(3)
               r2 = qp1 + qp2 + qp3
               r = dsqrt(r2)
               r5 = r2*r2*r
               temp(1) = 0.5d0 * (3.0d0*qp3 - r2)/ r5
               temp(2) = rt3*pvec(1)*pvec(3)/r5
               temp(3) = rt3*pvec(2)*pvec(3)/r5
               temp(4) = 0.5d0*rt3*(qp1 - qp2)/r5
               temp(5) = rt3*pvec(1)*pvec(2)/r5
               do kk=1,5
                  do k=1,nvar
                     u(i,k) = u(i,k) + temp(kk)*p(ivar+kk,k)
                  end do
                  esp1(i) = esp1(i) - temp(kk)*binvc(ivar+kk)
               end do
               ivar = ivar + 5
            endif
         end do
      end do

      call svd(mesp,nesp,nvar,u,w,.true.,u,.true.,v,ierr,scratch,np)

      wmax = 0.d0
      do j=1,nvar
         if (w(j).gt.wmax) wmax = w(j)
      end do     
      thresh = toll*wmax
      if (thresh.lt.1d-10) thresh = 1d-10

      nsvd = 0
      do j=1,nvar
         if (w(j).lt.thresh) then
             w(j) = 0.d0
             nsvd = nsvd + 1
         endif
      end do

      write(iun3,*) 'Rank of design matrix AP :  ',nvar-nsvd

c Solve for p`

      do i=1,nvar
         s = 0.d0
         if (w(i).ne.0.d0) then
             do j=1,nesp
                s = s + u(j,i)*esp1(j)
             end do
             s = s / w(i)
         endif
         scratch(i) = s
      end do
      do i=1,nvar
         s = 0.d0
         do j=1,nvar
            s = s + v(i,j)*scratch(j)
         end do
         c(i) = s
      end do

      
c Solution obeying constraints from q=BinvC + p*c
c stored in w

      do i=1,nvar
         w(i) = binvc(i)
         do j=1,nvar
           w(i) = w(i) + p(i,j)*c(j)
         end do
      end do

c     Store fitted multipoles  

      ivar = 0
      do i=1, natoms

         if (monopole(i)) then
            ivar = ivar + 1
            q(1,i) = w(ivar)
         else 
            q(1,i) = 0.d0
         endif

         if (dipole(i)) then
            do j=1,3
               ivar = ivar + 1
               q(1+j,i) = w(ivar)
            end do
         else
            do j=1,3
               q(1+j,i) = 0.d0
            end do
         endif

         if (quadrupole(i)) then
            do j=1,5
               ivar = ivar + 1
               q(4+j,i) = w(ivar)
            end do
         else
            do j=1,5
               q(4+j,i) = 0.d0
            end do
         endif

      end do
            
c
c     calculate root mean square fits and relative root mean square fits
c
      rms = 0.d0
      rrms = 0.d0
      do i=1,nesp
         espc = 0.d0
         call calc2(connl(1,i),connl(2,i),connl(3,i),espc)
         rms = rms + (espc-esp(i))**2
         rrms = rrms + esp(i)**2
      end do

      rms = dsqrt(rms/nesp)
      rrms = rms / dsqrt(rrms/nesp)
      rms = rms*627.51d0

      write(iun3,'(5x,''ATOM NO.    TYPE'')')
      write(iun3,*)' '

      qtot = 0.d0
      do i=1,natoms
         qtot = qtot + q(1,i)
         write(iun3,'(7x,i2,9x,a2,1x,f9.5)')i,elemnt(nat(i)),q(1,i)
         write(iun3,'(21x,3f9.5)') (q(j,i),j=2,4)
         write(iun3,'(21x,5f9.5)') (q(j,i),j=5,9)
      end do

      write(iun3,*)' '
      write(iun3,*) 'THE TOTAL CHARGE IS:       ',qtot
      write(iun3,*)' '

      write(iun3,*) 'THE NUMBER OF POINTS IS:   ',nesp
      write(iun3,*) 'THE RMS DEVIATION IS:      ',rms
      write(iun3,*) 'THE RRMS DEVIATION IS:     ',rrms

      write(iun3,*)' '
      write(iun3,*) 
     &     'DIPOLE MOMENT EVALUATED FROM Monopoles and Dipoles (debye)'
      write(iun3,'(12x,'' X        Y        Z       TOTAL'')')
      write(iun3,*)' '

      dipx = 0
      dipy = 0
      dipz = 0

      do i=1,natoms
         dipx = dipx + xyz(1,i)*q(1,i)
         dipy = dipy + xyz(2,i)*q(1,i)
         dipz = dipz + xyz(3,i)*q(1,i)
      end do

      dip = dsqrt(dipx**2+dipy**2+dipz**2)

      write(iun3,'(8x,4f9.4,4x,"Charge")')
     &              dipx*cf,dipy*cf,dipz*cf,dip*cf
      dipx = 0
      dipy = 0
      dipz = 0

      do i=1,natoms
         dipx = dipx + q(3,i)
         dipy = dipy + q(4,i)
         dipz = dipz + q(2,i)
      end do

      dip = dsqrt(dipx**2+dipy**2+dipz**2)
      write(iun3,'(8x,4f9.4,4x,"Dipole")')
     &                dipx*cf,dipy*cf,dipz*cf,dip*cf

      dipx = 0
      dipy = 0
      dipz = 0

      do i=1,natoms
         dipx = dipx + xyz(1,i)*q(1,i)
         dipy = dipy + xyz(2,i)*q(1,i)
         dipz = dipz + xyz(3,i)*q(1,i)
         dipx = dipx + q(3,i)
         dipy = dipy + q(4,i)
         dipz = dipz + q(2,i)
      end do

      dip = dsqrt(dipx**2+dipy**2+dipz**2)
      write(iun3,'(8x,4f9.4,4x,"Total/")') 
     &               dipx*cf,dipy*cf,dipz*cf,dip*cf

      dipo(1) = dipx
      dipo(2) = dipy
      dipo(3) = dipz
      ihsdp = 3

      call cpypol

      open(unit=46,form='formatted',file='esp.xyz',
     &    status='unknown',err=200)

      write(46,'(i5)') natoms
      write(46,'(a)') 'Molden MPOLEFIT fitted charges'
      do i=1,natoms
          write(46,'(a2,1x,3(f9.6,1x),f9.6)')
     &        elemnt(nat(i)),(xyz(j,i)*0.52917706d0,j=1,3),q(1,i)
      end do
      close(46)

      return

  200 write(iun3,*) 'Couldnt write XYZ file esp.xyz'
      end

    
      subroutine setmp(keywrd)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (lfmax=2)
      parameter (lfmaxf=(lfmax+1)*(lfmax+1))     
      logical monopole,dipole,quadrupole
      common /mpfit/ q(lfmaxf,numatm),
     &     monopole(numatm),dipole(numatm),quadrupole(numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      character*(*) keywrd
      dimension poles(numatm)

      indx = index(keywrd,'ATPOL')

      if (indx.ne.0) then
         do i=1,natoms
            monopole(i) = .false.
            dipole(i) = .false.
            quadrupole(i) = .false.
            poles(i) = 7
            if (nat(i).eq.1) poles(i) = 3
         end do
         call occin(indx,poles,numatm)
         do i=1,natoms
            j = int(poles(i))
            if (j.ge.4) then
               quadrupole(i) = .true.
               j = j - 4
            endif
            if (j.ge.2) then
               dipole(i) = .true.
               j = j - 2
            endif
            if (j.eq.1) monopole(i) = .true.
         end do
      else
         do i=1,natoms
            monopole(i) = .true.
            dipole(i) = .true.
            quadrupole(i) = .true.
            if (nat(i).eq.1) then
               quadrupole(i)=.false.
           endif
         end do
      endif

      return
      end

      subroutine calc2(x,y,z,pot)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (lfmax=2)
      parameter (lfmaxf=(lfmax+1)*(lfmax+1))
      parameter (numatm=2000)
      logical monopole,dipole,quadrupole
      common /mpfit/ q(lfmaxf,numatm),
     &     monopole(numatm),dipole(numatm),quadrupole(numatm)
      common /coord/ xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)

      dimension pvec(3)
c in this procedure the calculation order is inversed
c we use
c 
c   a      b      c      d       
c ---- + ---- + ---- + ---- = ((((d /x**2+c) /x**2)+b) /x**2) +a) /x
c   x    x**3   x**5   x**7   
c
      rt3=dsqrt(3.0d0)

      pot=0.0d0 
      do 10 i=1,natoms
         pvec(1)=x-xyz(1,i)
         pvec(2)=y-xyz(2,i)
         pvec(3)=z-xyz(3,i)
         qp1=pvec(1)*pvec(1)
         qp2=pvec(2)*pvec(2)
         qp3=pvec(3)*pvec(3)
         r2= qp1+qp2+qp3
         r=dsqrt(r2)
      pot = pot 
     & +(((
     & (q(5,i)*0.5d0*(3.0d0*qp3-r2)+
     & q(6,i)*rt3*pvec(1)*pvec(3)+
     & q(7,i)*rt3*pvec(2)*pvec(3)+
     & q(8,i)*0.5d0*rt3*(qp1-qp2)+
     & q(9,i)*rt3*pvec(1)*pvec(2)))/r2
     & +(q(2,i)*pvec(3)+q(3,i)*pvec(1)+q(4,i)*pvec(2)))/r2
     & +q(1,i))/r
10    continue

      return
      end 

      subroutine cpypol
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (numatm=2000)
      parameter (lfmax=2)
      parameter (lfmaxf=(lfmax+1)*(lfmax+1))
      parameter (mxsite=300)
      parameter (lmax=4)
      parameter (lmaxf=(lmax+1)*(lmax+1))
      logical monopole,dipole,quadrupole
      common /mpfit/ q(lfmaxf,numatm),
     &     monopole(numatm),dipole(numatm),quadrupole(numatm)
      character*8   ctag
      common /multip/ qmom(lmaxf,mxsite),car(3,mxsite),ctag(mxsite),
     &                nsites
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)

      do i=1,natoms
         do j=1,lmaxf
             qmom(j,i) = 0.0d0
         end do
         do j=1,lfmaxf
             qmom(j,i) = q(j,i)
         end do
         do j=1,3
             car(j,i) = xyz(j,i)
         end do
      end do

      nsites = natoms

      return
      end


      SUBROUTINE SVD (NM, M, N, A, W, MATU, U, MATV, V, IERR, RV1, nn)
C
C***PURPOSE  Perform the singular value decomposition of a rectangular
C            matrix.
C***LIBRARY   SLATEC
C
C     This subroutine is a translation of the ALGOL procedure SVD,
C     NUM. MATH. 14, 403-420(1970) by Golub and Reinsch.
C     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
C
C     This subroutine determines the singular value decomposition
C          T
C     A=USV  of a REAL M by N rectangular matrix.  Householder
C     bidiagonalization and a variant of the QR algorithm are used.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A, U, V, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C          Note that NM must be at least as large as the maximum
C          of M and N.
C
C        M is the number of rows of A and U.
C
C        N is the number of columns of A and U and the order of V.
C
C        A contains the rectangular input matrix to be decomposed.  A is
C          a two-dimensional REAL array, dimensioned A(NM,N).
C
C        MATU should be set to .TRUE. if the U matrix in the
C          decomposition is desired, and to .FALSE. otherwise.
C          MATU is a LOGICAL variable.
C
C        MATV should be set to .TRUE. if the V matrix in the
C          decomposition is desired, and to .FALSE. otherwise.
C          MATV is a LOGICAL variable.
C
C     On Output
C
C        A is unaltered (unless overwritten by U or V).
C
C        W contains the N (non-negative) singular values of A (the
C          diagonal elements of S).  They are unordered.  If an
C          error exit is made, the singular values should be correct
C          for indices IERR+1, IERR+2, ..., N.  W is a one-dimensional
C          REAL array, dimensioned W(N).
C
C        U contains the matrix U (orthogonal column vectors) of the
C          decomposition if MATU has been set to .TRUE.  Otherwise,
C          U is used as a temporary array.  U may coincide with A.
C          If an error exit is made, the columns of U corresponding
C          to indices of correct singular values should be correct.
C          U is a two-dimensional REAL array, dimensioned U(NM,N).
C
C        V contains the matrix V (orthogonal) of the decomposition if
C          MATV has been set to .TRUE.  Otherwise, V is not referenced.
C          V may also coincide with A if U does not.  If an error
C          exit is made, the columns of V corresponding to indices of
C          correct singular values should be correct.  V is a two-
C          dimensional REAL array, dimensioned V(NM,N).
c+++ Changed: V can not coincide with A, as for memory reasons v is limited
c to V(nn,nn). 

C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          K          if the K-th singular value has not been
C                     determined after 30 iterations.
C
C        RV1 is a one-dimensional REAL array used for temporary storage,
C          dimensioned RV1(N).
c+++Changed: NN must be set to the row dimension of the two-dimensional
C          array parameters V, as declared in the calling
C          program dimension statement.  NN is an INTEGER variable.
C          Note that NN must be at least as large as the maximum of N.

c
c

C

C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
      INTEGER I,J,K,L,M,N,II,I1,KK,K1,LL,L1,MN,NM,ITS,IERR,nn
      double precision A(NM,*),W(*),U(NM,*),V(nn,*),RV1(*)
      double precision C,F,G,H,S,X,Y,Z,SCALE,S1
      LOGICAL MATU,MATV
C
      IERR = 0

C
      DO 100 I = 1, M
C
         DO 100 J = 1, N
            U(I,J) = A(I,J)
  100 CONTINUE
C     .......... HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM ..........
      G = 0.0d0
      SCALE = 0.0d0
      S1 = 0.0d0
C
      DO 300 I = 1, N
         L = I + 1
         RV1(I) = SCALE * G
         G = 0.0d0
         S = 0.0d0
         SCALE = 0.0d0
         IF (I .GT. M) GO TO 210
C
         DO 120 K = I, M
  120    SCALE = SCALE + ABS(U(K,I))
C
         IF (SCALE .EQ. 0.0d0) GO TO 210
C
         DO 130 K = I, M
            U(K,I) = U(K,I) / SCALE
            S = S + U(K,I)**2
  130    CONTINUE
C
         F = U(I,I)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         U(I,I) = F - G
         IF (I .EQ. N) GO TO 190
C
         DO 150 J = L, N
            S = 0.0d0
C
            DO 140 K = I, M
  140       S = S + U(K,I) * U(K,J)
C
            F = S / H
C
            DO 150 K = I, M
               U(K,J) = U(K,J) + F * U(K,I)
  150    CONTINUE
C
  190    DO 200 K = I, M
  200    U(K,I) = SCALE * U(K,I)
C
  210    W(I) = SCALE * G
         G = 0.0d0
         S = 0.0d0
         SCALE = 0.0d0
         IF (I .GT. M .OR. I .EQ. N) GO TO 290
C
         DO 220 K = L, N
  220    SCALE = SCALE + ABS(U(I,K))
C
         IF (SCALE .EQ. 0.0d0) GO TO 290
C
         DO 230 K = L, N
            U(I,K) = U(I,K) / SCALE
            S = S + U(I,K)**2
  230    CONTINUE
C
         F = U(I,L)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         U(I,L) = F - G
C
         DO 240 K = L, N
  240    RV1(K) = U(I,K) / H
C
         IF (I .EQ. M) GO TO 270
C
         DO 260 J = L, M
            S = 0.0d0
C
            DO 250 K = L, N
  250       S = S + U(J,K) * U(I,K)
C
            DO 260 K = L, N
               U(J,K) = U(J,K) + S * RV1(K)
  260    CONTINUE
C
  270    DO 280 K = L, N
  280    U(I,K) = SCALE * U(I,K)
C
  290    S1 = MAX(S1,ABS(W(I))+ABS(RV1(I)))
  300 CONTINUE
C     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS ..........
      IF (.NOT. MATV) GO TO 410
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 400 II = 1, N
         I = N + 1 - II
         IF (I .EQ. N) GO TO 390
         IF (G .EQ. 0.0d0) GO TO 360
C
         DO 320 J = L, N
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
  320    V(J,I) = (U(I,J) / U(I,L)) / G
C
         DO 350 J = L, N
            S = 0.0d0
C
            DO 340 K = L, N
  340       S = S + U(I,K) * V(K,J)
C
            DO 350 K = L, N
               V(K,J) = V(K,J) + S * V(K,I)
  350    CONTINUE
C
  360    DO 380 J = L, N
            V(I,J) = 0.0d0
            V(J,I) = 0.0d0
  380    CONTINUE
C
  390    V(I,I) = 1.0d0
         G = RV1(I)
         L = I
  400 CONTINUE
C     .......... ACCUMULATION OF LEFT-HAND TRANSFORMATIONS ..........
  410 IF (.NOT. MATU) GO TO 510
C     ..........FOR I=MIN(M,N) STEP -1 UNTIL 1 DO -- ..........
      MN = N
      IF (M .LT. N) MN = M
C
      DO 500 II = 1, MN
         I = MN + 1 - II
         L = I + 1
         G = W(I)
         IF (I .EQ. N) GO TO 430
C
         DO 420 J = L, N
  420    U(I,J) = 0.0d0
C
  430    IF (G .EQ. 0.0d0) GO TO 475
         IF (I .EQ. MN) GO TO 460
C
         DO 450 J = L, N
            S = 0.0d0
C
            DO 440 K = L, M
  440       S = S + U(K,I) * U(K,J)
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            F = (S / U(I,I)) / G
C
            DO 450 K = I, M
               U(K,J) = U(K,J) + F * U(K,I)
  450    CONTINUE
C
  460    DO 470 J = I, M
  470    U(J,I) = U(J,I) / G
C
         GO TO 490
C
  475    DO 480 J = I, M
  480    U(J,I) = 0.0d0
C
  490    U(I,I) = U(I,I) + 1.0d0
  500 CONTINUE
C     .......... DIAGONALIZATION OF THE BIDIAGONAL FORM ..........
  510 CONTINUE
C     .......... FOR K=N STEP -1 UNTIL 1 DO -- ..........
      DO 700 KK = 1, N
         K1 = N - KK
         K = K1 + 1
         ITS = 0
C     .......... TEST FOR SPLITTING.
C                FOR L=K STEP -1 UNTIL 1 DO -- ..........
  520    DO 530 LL = 1, K
            L1 = K - LL
            L = L1 + 1
            IF (S1 + ABS(RV1(L)) .EQ. S1) GO TO 565
C     .......... RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
            IF (S1 + ABS(W(L1)) .EQ. S1) GO TO 540
  530    CONTINUE
C     .......... CANCELLATION OF RV1(L) IF L GREATER THAN 1 ..........
  540    C = 0.0d0
         S = 1.0d0
C
         DO 560 I = L, K
            F = S * RV1(I)
            RV1(I) = C * RV1(I)
            IF (S1 + ABS(F) .EQ. S1) GO TO 565
            G = W(I)
            H = dsqrt(f*f+g*g)
            W(I) = H
            C = G / H
            S = -F / H
            IF (.NOT. MATU) GO TO 560
C
            DO 550 J = 1, M
               Y = U(J,L1)
               Z = U(J,I)
               U(J,L1) = Y * C + Z * S
               U(J,I) = -Y * S + Z * C
  550       CONTINUE
C
  560    CONTINUE
C     .......... TEST FOR CONVERGENCE ..........
  565    Z = W(K)
         IF (L .EQ. K) GO TO 650
C     .......... SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
         IF (ITS .EQ. 30) GO TO 1000
         ITS = ITS + 1
         X = W(L)
         Y = W(K1)
         G = RV1(K1)
         H = RV1(K)
         F = 0.5d0 * (((G + Z) / H) * ((G - Z) / Y) + Y / H - H / Y)
         G = dsqrt(F*F+1.0d0)
         F = X - (Z / X) * Z + (H / X) * (Y / (F + SIGN(G,F)) - H)
C     .......... NEXT QR TRANSFORMATION ..........
         C = 1.0d0
         S = 1.0d0
C
         DO 600 I1 = L, K1
            I = I1 + 1
            G = RV1(I)
            Y = W(I)
            H = S * G
            G = C * G
            Z = dsqrt(F*f+H*h)
            RV1(I1) = Z
            C = F / Z
            S = H / Z
            F = X * C + G * S
            G = -X * S + G * C
            H = Y * S
            Y = Y * C
            IF (.NOT. MATV) GO TO 575
C
            DO 570 J = 1, N
               X = V(J,I1)
               Z = V(J,I)
               V(J,I1) = X * C + Z * S
               V(J,I) = -X * S + Z * C
  570       CONTINUE
C
  575       Z = dsqrt(F*f+H*h)
            W(I1) = Z
C     .......... ROTATION CAN BE ARBITRARY IF Z IS ZERO ..........
            IF (Z .EQ. 0.0d0) GO TO 580
            C = F / Z
            S = H / Z
  580       F = C * G + S * Y
            X = -S * G + C * Y
            IF (.NOT. MATU) GO TO 600
C
            DO 590 J = 1, M
               Y = U(J,I1)
               Z = U(J,I)
               U(J,I1) = Y * C + Z * S
               U(J,I) = -Y * S + Z * C
  590       CONTINUE
C
  600    CONTINUE
C
         RV1(L) = 0.0d0
         RV1(K) = F
         W(K) = X
         GO TO 520
C     .......... CONVERGENCE ..........
  650    IF (Z .GE. 0.0d0) GO TO 700
C     .......... W(K) IS MADE NON-NEGATIVE ..........
         W(K) = -Z
         IF (.NOT. MATV) GO TO 700
C
         DO 690 J = 1, N
  690    V(J,K) = -V(J,K)
C
  700 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO A
C                SINGULAR VALUE AFTER 30 ITERATIONS ..........
 1000 IERR = K
 1001 RETURN
      END
