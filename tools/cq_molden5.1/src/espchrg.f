      subroutine espchrg(valc,nvalc,doiso,dmachg)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /coord/ xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      character keywrd*320,keyori*320
      common /keywrd/ keywrd,keyori
      parameter (max3d=61)
      parameter (mesp = (max3d*max3d*max3d+max3d)/4)
      parameter (mdum = (max3d*max3d*max3d+max3d) - mesp*4)
      common /spa3d/  connl(3,mesp), esp(mesp), dum(mdum)
      logical keyr,keyi,dmachg,doiso
      dimension valc(*)

c
c     KEYWORDS:                                              Default
c
c     ESPCH          triggers calculation esp charges
c     DMACH          triggers calculation dmaesp charges
c     MPFIT          triggers calculation Multipole Fit
c     NUMSURF=n      the number of connolly surfaces           (4)
c     CONNSC=n.n     initial scalefactor vdwaals radii         (1.4)
c     CONNINCR=n.n   increment scalefactor for next surface    (0.2)
c     PTDEN=n.n      Density of points per Unit Area           (3.0)
c
c     DIPX=n.n       specifies dipole moment to be fitted to
c     DIPY=n.n       specifies dipole moment to be fitted to
c     DIPZ=n.n       specifies dipole moment to be fitted to
c
c     The default is not to fit the dipole moment
c
c     When constructing atomspheres use as principle axis :
c     (if you want to reproduce symmetry)
c
c     AXIS-X,AXIS-Y,AXIS-Z                                    AXIS-Z
c
c     After BESLER,MERZ,KOLLMAN J. COMPUT. CHEM.
c
C     Added keyword (B. P. van Eijck, 2006):
C     CHADD     read non-atomic charge centers from file "molden.def"

      if (.not.keyr(keywrd,'CONNSC',scale)) then
          scale = 1.4d0
      endif

      if (.not.keyr(keywrd,'CONNINCR',scincr)) then
          scincr = 0.20d0
      endif

      if (keyi(keywrd,'NUMSURF',nsurf)) then
          if (doiso) nvalc = nsurf
      else
          nsurf = 4
      endif

c     Density of points per Unit Area

      if (.not.keyr(keywrd,'PTDEN',den)) then
          den = 3.0d0
      endif
      
c     Fit to dipole as well ?

      idip = 0

c     dx,dy,dz not used when idip = 0


      if (keyr(keywrd,'DIPX',dx).and.keyr(keywrd,'DIPY',dy).and.
     &    keyr(keywrd,'DIPZ',dz)) then
          idip = 1
      else
          dx = 0.0d0
          dy = 0.0d0
          dz = 0.0d0
      endif

      iaxis = 3
      if (index(keywrd,'AXIS-X').ne.0) then
          iaxis = 1      
      endif
      if (index(keywrd,'AXIS-Y').ne.0) then
          iaxis = 2      
      endif
      if (index(keywrd,'AXIS-Z').ne.0) then
          iaxis = 3      
      endif
c
c     now calculate the surface points
c
      nesp = 0
      if (doiso) then
         print*,'======================================'
         print*,'charges from isodensity surfaces'
         print*,'Number of surface(s): ',nvalc
         print*,'Surface contour value(s): ',(valc(j),j=1,nvalc)
         print*,'======================================'
         call isoden(valc,nvalc,scincr,nesp,-1)
      else
         do i = 1,nsurf
            call connol(scale,den,nesp,iaxis)
            scale = scale + scincr
         end do
      endif
c
c     next calculate the esp at the points calculated by connol
c
      do i = 1,nesp
         if (dmachg) then
            call calc(connl(1,i),connl(2,i),connl(3,i),esp(i))
         else
            call espot(connl(1,i),connl(2,i),connl(3,i),esp(i),0)
         endif
      end do

      ichadd = 0

      if (index(keywrd,'CHADD').ne.0) then
c
c        Read non-atomic charge sites from file CHADD
c
         print*,' Reading non-atomic charge sites'
         open(48,file='molden.def',status='OLD',err=200)
         ichadd = 1
         goto 210
  200    continue
         stop ' ***** No file MOLDEN.DEF present. Continue?'
  210    continue
      endif

      ichrg = 0
      if (index(keywrd,'MPFIT').ne.0) then
         call mpolefit(idip,nesp,esp,connl,dx,dy,dz,ichrg,
     &                 dmachg,keywrd,ichadd)
      else
         if (dmachg) then
            idmach = 1
         else
            idmach = 0
         endif
         call espfit(idip,nesp,esp,connl,dx,dy,dz,ichrg,idmach,ichadd)
      endif

      if (index(keywrd,'ARESP').ne.0) call aresp(nesp,ichrg)

      return
      end

      subroutine espfid(idip,nesp,esp,connl,dx,dy,dz,iz,idmach,ichadd,
     &                  qat)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension al((numatm+4)**2),a(numatm,numatm),b(numatm+4)
      dimension esp(*),connl(3,*),qat(*)

c
c     conversion factor for debye to atomic units
c
      cf = 5.2917715d-11*1.601917d-19/3.33564d-30

c     calculate total charge on the molecule

      iz = 0
      do i=1,natoms
          iz = iz + nat(i)
      end do
      iz = iz - nelecs

      if (iz.ne.0) write(iun3,*) 'Charge of molecule = ',iz
c      print*,'charge nuclei ',iz,' nelec ',nelecs

      rms = 0.0d0
      rrms = 0.0d0

      do i=natoms+1,numatm
         nat(i) = 99
      end do

      elemnt(99) = ' X'

c     Read additional charge sites from file CHADD

      if (ichadd.ne.0) then

         toang = 0.52917706d0
         call chadd(natoms)

         do j=1,natoms
            write(iun3,'(i4,3f10.5,i4,a)') j,(xyz(m,j)*toang,m=1,3),
     &                                     nat(j),elemnt(nat(j))
         end do

      endif

c
c     the following sets up the linear equation a*q=b
c     set up the a(j,k) array
c
      do j=1,natoms
         b(j) = 0.d0
      end do
      do i=1,nesp
         do j=1,natoms
            rij = dsqrt((xyz(1,j)-connl(1,i))**2 +
     &      (xyz(2,j)-connl(2,i))**2 + (xyz(3,j)-connl(3,i))**2)
            b(j) = b(j) + esp(i)*1.d0/rij
         end do
      end do

      do k=1,natoms
         do j=1,natoms
            a(j,k) = 0.0d0
            do i=1,nesp
               rik = dsqrt((xyz(1,k)-connl(1,i))**2 + 
     &         (xyz(2,k)-connl(2,i))**2 + (xyz(3,k)-connl(3,i))**2)
               rij = dsqrt((xyz(1,j)-connl(1,i))**2 + 
     &         (xyz(2,j)-connl(2,i))**2 + (xyz(3,j)-connl(3,i))**2)
               a(j,k) = a(j,k) + 1.d0/rik*1.d0/rij
            end do
         end do
         a(natoms+1,k) = 1.d0
         a(k,natoms+1) = 1.d0
         a(natoms+1,natoms+1) = 0.d0
         if(idip .eq. 1) then
            a(natoms+2,k) = xyz(1,k)
            a(k,natoms+2) = xyz(1,k)
            a(natoms+2,natoms+2) = 0.d0
            a(natoms+3,k) = xyz(2,k)
            a(k,natoms+3) = xyz(2,k)
            a(natoms+3,natoms+3) = 0.d0
            a(natoms+4,k) = xyz(3,k)
            a(k,natoms+4) = xyz(3,k)
            a(natoms+4,natoms+4) = 0.d0
         endif
      end do
      b(natoms+1) = dfloat(iz)
      b(natoms+2) = dx/cf
      b(natoms+3) = dy/cf
      b(natoms+4) = dz/cf


c
c     insert charge and dipolar (if desired) constraints
c
      if (idip.eq.1) then
         ndim = natoms + 4
      else
         ndim = natoms + 1
      endif

      l = 0
      do i=1,ndim
         do j=1,ndim
            l = l + 1
            al(l) = a(i,j)
         end do
      end do
      call matinv(al,ndim,det)
      l = 0
      do i=1,ndim
         do j=1,ndim
            l = l + 1
            a(i,j) = al(l)
         end do
      end do
      do i=1,ndim
         qat(i) = 0.d0
         do j=1,ndim
            qat(i) = qat(i) + a(i,j)*b(j)
         end do
      end do
c
c     calculate root mean square fits and relative root mean square fits
c
      do i=1,nesp
         espc = 0.d0
         do j=1,natoms
            rij=dsqrt((xyz(1,j)-connl(1,i))**2+(xyz(2,j)-connl(2,i))**2
     &      +(xyz(3,j)-connl(3,i))**2)
            espc = espc + qat(j)/rij
         end do
         rms = rms + (espc-esp(i))**2
         rrms = rrms + esp(i)**2
      end do
      rms = dsqrt(rms/nesp)
      rrms = rms / dsqrt(rrms/nesp)
      rms = rms*627.51d0

      write(iun3,'(15x,''ATOM NO.    TYPE    CHARGE'')')
      write(iun3,*)' '
      qtot = 0.d0
      do i=1,natoms
         qtot = qtot + qat(i)
         write(iun3,'(17x,i2,9x,a2,1x,f10.4)')i,elemnt(nat(i)),qat(i)
      end do

      write(iun3,*)' '
      write(iun3,*) 'THE TOTAL CHARGE IS:       ',qtot
      write(iun3,*)' '

      write(iun3,*) 'THE NUMBER OF POINTS IS:   ',nesp
      write(iun3,*) 'THE RMS DEVIATION IS (kcal/mol):      ',rms
      write(iun3,*) 'THE RRMS DEVIATION IS:     ',rrms

c      if (iz.ne.0) return

      write(iun3,*)' '
      write(iun3,*) 
     & 'DIPOLE MOMENT EVALUATED FROM POINT CHARGES (debye)'
      write(iun3,'(12x,'' X        Y        Z       TOTAL'')')
      write(iun3,*)' '

      dipx = 0
      dipy = 0
      dipz = 0
      do i=1,natoms
         dipx = dipx + xyz(1,i)*qat(i)
         dipy = dipy + xyz(2,i)*qat(i)
         dipz = dipz + xyz(3,i)*qat(i)
      end do
      dip = dsqrt(dipx**2+dipy**2+dipz**2)

      write(iun3,'(8x,4f9.4)') dipx*cf,dipy*cf,dipz*cf,dip*cf

      ihasq = 1

      jmode = 0
      if (idmach.eq.1) jmode = 1
      call wrxyz(jmode)

      return
      end

      subroutine connol(scale,dens,nesp,iaxis)
      implicit double precision (a-h,o-z)
      parameter (mxptn=5000)
      parameter (numatm=2000)
      parameter (mxel=100)
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      parameter (max3d=61)
      parameter (mesp = (max3d*max3d*max3d+max3d)/4)
      parameter (mdum = (max3d*max3d*max3d+max3d) - mesp*4)
      common /spa3d/  connl(3,mesp), esp(mesp), dum(mdum)
      common /esprad/ vander(mxel)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension rad(numatm)
      dimension con(3,mxptn)
      dimension inbr(200),cnbr(3,200),rnbr(200)
      logical mnbr(200)
      dimension ci(3), temp0(3)
      dimension cw(3)
      logical collid

      pi = 4.d0*datan(1.d0)
      bohr = 0.529177249d+00
      iop = 1
      den = dens
      do i=1,natoms
         rad(i) = vander(nat(i))*scale/bohr
         if (rad(i) .lt. 0.01d0) then
            write(iun3,*) 'zero van der waals radius for atom ',nat(i)
         endif
      end do
c
c     big loop for each atom
c
      do iatom = 1, natoms
c
c     transfer values from large arrays to iatom variables
c
         ri = rad(iatom)
         do k = 1,3
            ci(k) = xyz(k,iatom)
         end do
c
c     gather the neighboring atoms of iatom
c
         nnbr = 0
         do jatom = 1, natoms
            if (iatom .eq. jatom) goto 60
            d2 = dist2(ci,xyz(1,jatom))
            if (d2 .ge. (ri+rad(jatom)) ** 2) goto 60
c
c     we have a new neighbor
c     transfer atom coordinates, radius and surface request number
c
            nnbr = nnbr + 1
            if (nnbr .gt. 200)then
               write(iun3,*) 'too many neighbors:',nnbr
               stop
            endif
            inbr(nnbr) = jatom
            do k = 1,3
               cnbr(k,nnbr) = xyz(k,jatom)
            end do
            rnbr(nnbr) = rad(jatom)
60          continue
         end do
c
c     contact surface
c
         ncon = (4 * pi * ri ** 2) * den
         if (ncon .gt. mxptn) then
             write(iun3,*) ncon,' points requested'
             ncon = mxptn
             write(iun3,*) 'restrain grid to ',ncon,' pts /atom'
         endif
c
c     this call may decrease ncon somewhat
c
         if ( ncon .eq. 0) then
            write(iun3,*) 'vector length of zero in connol'
            stop
         endif
         call mksph1(con,ncon,iaxis)
c
c     contact probe placement loop
c
         do i = 1,ncon
            do k = 1,3
               cw(k) = ci(k) + ri * con(k,i)
            end do
c
c     check for collision with neighboring atoms
c
            if (collid(cw,cnbr,rnbr,mnbr,nnbr,1,
     &      jnbr,knbr)) goto 100
            do kk=1,3
               temp0(kk) = ci(kk)+ri*con(kk,i)
            end do
c
c     store point in connl and increment nesp
c
            nesp = nesp + 1
            if (nesp .gt. mesp) then
               write(iun3,*) 'to many points generated in connol'
               write(iun3,*) ' reduce nsurf, scale, den, or scincr'
               stop
            endif
            connl(1,nesp) = temp0(1)
            connl(2,nesp) = temp0(2)
            connl(3,nesp) = temp0(3)
  100       continue
         end do
  110    continue
      end do

      return
      end

      logical function collid(cw,cnbr,rnbr,mnbr,nnbr,ishape,
     &                        jnbr,knbr)
      implicit double precision (a-h,o-z)
      dimension cw(3)
      dimension cnbr(3,200)
      dimension rnbr(200)
      logical mnbr(200)
      if (nnbr .le. 0) goto 20
c
c     check whether probe is too close to any neighbor
c
      do 10 i = 1, nnbr
         if (ishape .gt. 1 .and. i .eq. jnbr) goto 10
         if (ishape .eq. 3 .and. (i .eq. knbr .or. .not. mnbr(i)))
     1   goto 10
         sumrad = rnbr(i)
         vect1 = dabs(cw(1) - cnbr(1,i))
         if (vect1 .ge. sumrad) goto 10
         vect2 = dabs(cw(2) - cnbr(2,i))
         if (vect2 .ge. sumrad) goto 10
         vect3 = dabs(cw(3) - cnbr(3,i))
         if (vect3 .ge. sumrad) goto 10
         sr2 = sumrad * sumrad
         dd2 = vect1*vect1 + vect2*vect2 + vect3*vect3
         if (dd2 .lt. sr2) goto 30
   10 continue
   20 continue
      collid = .false.
      goto 40
   30 continue
      collid = .true.
   40 continue
      return
      end

      subroutine mksph1(u,n,iaxis)
      implicit double precision (a-h,o-z)
      dimension u(3,n)

      pi = 4.d0*datan(1.d0)
      nequat = dsqrt(n * pi)
      nvert = nequat/2
      nu = 0
      do i = 1,nvert+1
         fi = (pi * (i-1)) / nvert
         z = dcos(fi)
         xy = dsin(fi)
         nhor = nequat * xy
         if (nhor .lt. 1) nhor = 1
         do j = 1,nhor
            fj = (2.d0 * pi * (j-1)) / nhor
            x = dcos(fj) * xy
            y = dsin(fj) * xy
            if (nu .ge. n) goto 30
            nu = nu + 1
            if (iaxis.eq.1) then
               u(1,nu) = z
               u(2,nu) = y
               u(3,nu) = x
            elseif (iaxis.eq.2) then
               u(1,nu) = x
               u(2,nu) = z
               u(3,nu) = y
            elseif (iaxis.eq.3) then
               u(1,nu) = x
               u(2,nu) = y
               u(3,nu) = z
            endif
         end do
      end do
   30 continue

      n = nu

      return
      end

c      subroutine mksph2(c,norder,npts,t)
c      implicit double precision (a-h,o-z)
c      dimension c(3,*)
c
cc     Extended octahedron
c
c      pi = 4.d0*datan(1.d0)
c      r2 = dsqrt(2.0d0)/2.0d0
c      r3 = dsqrt(3.0d0)/3.0d0
c      r6 = dsqrt(6.0d0)/6.0d0
c
c      npts = 6
c
c      c(1,1) = -r6
c      c(2,1) = r2
c      c(3,1) = r3
c
c      c(1,2) = -r6
c      c(2,2) = -r2
c      c(3,2) = r3
c
c      c(1,3) = 2.0d0*r6
c      c(2,3) = 0.0d0
c      c(3,3) = r3
c
c      c(1,4) = r6
c      c(2,4) = r2
c      c(3,4) = -r3
c
c      c(1,5) = r6
c      c(2,5) = -r2
c      c(3,5) = -r3
c
c      c(1,6) = -2.0d0*r6
c      c(2,6) = 0.0d0
c      c(3,6) = -r3
c
cc      t = dsqrt(2.0d0)
cc      return
c
c      do n=1,norder-1
c         t = dsqrt(2.0d0*(1.0d0 - dcos(pi/2.0d0**n)))
c         t2 = t*1.40d0
c         nn = npts
c         do i=1,nn
c            do j=i+1,nn
c               if (dist2(c(1,i),c(1,j))-t*t.lt.0.01d0) then
cc                if (dabs(dist2(c(1,i),c(1,j))-t2*t2).lt.0.01d0) then
cc               if (dist2(c(1,i),c(1,j)).lt.1.01d0) then
c                  npts = npts + 1
c                  call newpt(c(1,i),c(1,j),c(1,npts))
c               endif
c            end do
c         end do
c      end do
c
c      t = dsqrt(2.0d0*(1.0d0 - dcos(pi/2.0d0**norder)))
c
c      return
c      end
c
c      subroutine mksph3(cc,npts,t)
c      implicit double precision (a-h,o-z)
c      dimension cc(3,*)
c
cc     Extended icosahedron
c
c      pi = 4.d0*datan(1.d0)
c      r5 = dsqrt(5.0d0)
c      a = dsqrt(10.0d0+2.0d0*r5)/2.0d0
c      b = dsqrt(50.0d0+10.0d0*r5)/10.0d0
c      c = dsqrt(50.0d0-10.0d0*r5)/10.0d0
c      d = dsqrt(25.0d0+10.0d0*r5)/5.0d0
c      e = (1.0d0+r5)/2.0d0
c      d0 = 0.0d0
c      d1 = 1.0d0
c
c      npts = 12
c
c      cc(1,1) = d0
c      cc(2,1) = d0
c      cc(3,1) = a
c
c      cc(1,2) = 2.0d0*b
c      cc(2,2) = d0
c      cc(3,2) = b
c
c      cc(1,3) = c
c      cc(2,3) = e
c      cc(3,3) = b
c
c      cc(1,4) = -d
c      cc(2,4) = d1
c      cc(3,4) = b
c
c      cc(1,5) = -d
c      cc(2,5) = -d1
c      cc(3,5) = b
c
c      cc(1,6) = c
c      cc(2,6) = -e
c      cc(3,6) = b
c
c      cc(1,7) = d
c      cc(2,7) = d1
c      cc(3,7) = -b
c
c      cc(1,8) = -c
c      cc(2,8) = e
c      cc(3,8) = -b
c
c      cc(1,9) = -2.0d0*b
c      cc(2,9) = d0
c      cc(3,9) = -b
c
c      cc(1,10) = -c
c      cc(2,10) = -e
c      cc(3,10) = -b
c
c      cc(1,11) = d
c      cc(2,11) = -d1
c      cc(3,11) = -b
c
c      cc(1,12) = d0
c      cc(2,12) = d0
c      cc(3,12) = -a
c
c      do i=1,npts
c         do j=1,3
c            cc(j,i) = cc(j,i)/(b*r5)
c         end do
c      end do
c
c      t = 2.0d0/(b*r5)
c
c         nn = npts
c         do i=1,nn
c            do j=i+1,nn
c               if (dabs(dist2(cc(1,i),cc(1,j))-t*t).lt.0.01d0) then
c                  npts = npts + 1
c                  call newpt(cc(1,i),cc(1,j),cc(1,npts))
c               endif
c            end do
c         end do
c
c      t = 0.765d0
c
c      return
c      end

      subroutine mksph4(c,rad,grdw,cn,npts)
      implicit double precision (a-h,o-z)
      dimension c(3),cn(3,*)

      npts = 0
      rad2 = rad*rad

      imxz = (c(3)+rad)/grdw
      imnz = (c(3)-rad)/grdw

      do i=imnz,imxz
          z = i*grdw
          zz = dabs(z-c(3))
          r12 = rad2 - zz*zz
          if (r12.gt.0.0d0) then
            r1 = dsqrt(r12)
            imxx = (c(1)+r1)/grdw
            imnx = (c(1)-r1)/grdw
            imxy = (c(2)+r1)/grdw
            imny = (c(2)-r1)/grdw
            do j=imnx,imxx
                x = j*grdw
                xx = x-c(1)
                yy = r1*r1 - xx*xx
                if (yy.gt.0.0d0) then
                   yy = dsqrt(yy)
                   npts = npts + 1
                   cn(1,npts) = x
                   cn(2,npts) = c(2) + yy
                   cn(3,npts) = z
                   npts = npts + 1
                   cn(1,npts) = x
                   cn(2,npts) = c(2) - yy
                   cn(3,npts) = z
                endif
            end do
            do j=imny,imxy
                y = j*grdw
                yy = y-c(2)
                xx = r1*r1 - yy*yy
                if (xx.gt.0.0d0) then
                   xx = dsqrt(xx)
                   npts = npts + 1
                   cn(1,npts) = c(1) + xx
                   cn(2,npts) = y
                   cn(3,npts) = z
                   npts = npts + 1
                   cn(1,npts) = c(1) - xx
                   cn(2,npts) = y
                   cn(3,npts) = z
                endif
            end do
          endif
      end do

      imxx = (c(1)+rad)/grdw
      imnx = (c(1)-rad)/grdw

      imxy = (c(2)+rad)/grdw
      imny = (c(2)-rad)/grdw


      do i=imnx,imxx
          x = i*grdw
          xx = x-c(1)
          do j=imny,imxy
             y = j*grdw
             yy = y-c(2)
             r1 = xx*xx+yy*yy
             if (r1.lt.rad2) then
                zz = dsqrt(rad2-r1)
                npts = npts + 1
                cn(1,npts) = x
                cn(2,npts) = y
                cn(3,npts) = c(3) + zz
                npts = npts + 1
                cn(1,npts) = x
                cn(2,npts) = y
                cn(3,npts) = c(3) - zz
             endif
          end do
      end do

      return
      end

      subroutine mkscod(c,rad,grdw,iptr,npts,coo,iconn)
      implicit double precision (a-h,o-z)
      parameter (mxipts=500)
      parameter (mxcon=10)
      common /athlp/ iatoms, mxnat
      dimension c(3),ic(mxipts),v1(3),v2(3),d(4),ii(4)
      dimension coo(3,*),iconn(mxcon+1,*)

      gr2 = grdw*grdw*2.000001d0
      tol = 0.001d0

      do n=1,3

        imx = dint((c(n)+rad)/grdw)
        imn = dint((c(n)-rad)/grdw)
  
        do i=imn,imx

            z = dble(i)*grdw
            nc = 0

            do j=1,npts
               if (dabs(coo(n,iptr+j)-z).lt.tol) then
                  if (nc.lt.mxipts) then
                     nc = nc + 1
                     ic(nc) = iptr+j
                  endif
               endif
            end do

            do j=1,nc

               jj = ic(j)

               do l=1,4
                  d(l)  = 10000.0d0
                  ii(l) = 0
               end do

               do k=1,nc

                  if (j.ne.k) then

                    dd2 = dist2(coo(1,ic(j)),coo(1,ic(k)))

                    il = 0
                    if (dd2.lt.d(1)) then
                       il = 1
                    elseif (dd2.lt.d(2)) then
                       il = 2
                    elseif (dd2.lt.d(3)) then
                       il = 3
                    elseif (dd2.lt.d(4)) then
                       il = 4
                    endif

                    if (il.ne.0) then

                       do l=4,1+il,-1
                          ii(l) = ii(l-1)
                          d(l)  = d(l-1)
                       end do

                       ii(il) = k
                       d (il) = dd2

                    endif

                  endif

               end do

               if (ii(1).ne.0.and.d(1).le.gr2) then
                  iconn(1,jj) = iconn(1,jj) + 1
                  iconn(iconn(1,jj)+1,jj) = ic(ii(1))
               endif

               if (ii(1).ne.0.and.ii(2).ne.0) then

                  it = 0

                  do l=1,3
                     v1(l) = coo(l,jj) - coo(l,ic(ii(1)))
                     v2(l) = coo(l,jj) - coo(l,ic(ii(2)))
                  end do

                  v = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

                  if (v.lt.0.0d0) then

                     it = ii(2)
                     dt = d(2)

                  elseif (ii(3).ne.0) then

                     do l=1,3
                        v2(l) = coo(l,jj) - coo(l,ic(ii(3)))
                     end do

                     v = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

                     if (v.lt.0.0d0) then
                        it = ii(3)
                        dt = d(3)
                     else
                        it = ii(4)
                        dt = d(4)
                     endif

                  endif

                  if ((it.gt.0.and.it.le.mxipts).and.dt.le.gr2) then

                     if (jj.gt.0.and.jj.lt.mxnat) then

                        ijj = iconn(1,jj)

                        if (ijj.lt.mxcon) then
                           iconn(1,jj) = ijj + 1
                           iconn(iconn(1,jj)+1,jj) = ic(it)
                        endif

                     endif

                  endif

               endif
            end do
        end do
      end do

      return
      end

      subroutine newpt(c1,c2,c3)
      implicit double precision (a-h,o-z)
      dimension c1(3),c2(3),c3(3)

      do i=1,3
         c3(i) = c1(i) + c2(i)
      end do

      cl = vlen(c3)

      do i=1,3
         c3(i) = c3(i)/cl
      end do

      return
      end

      double precision function dist2(a,b)
c
c     determine distances between neighboring atoms
c
      implicit double precision (a-h,o-z)
      dimension a(3)
      dimension b(3)

      d1 = a(1)-b(1)
      d2 = a(2)-b(2)
      d3 = a(3)-b(3)

      dist2 = d1*d1 + d2*d2 + d3*d3
c      dist2 = (a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2

      return
      end

      subroutine matinv(a,n,d)
      implicit double precision (a-h,o-z)
      dimension a(*)
c
c     inverts a general square matrix
c     
      dimension l(500), m(500)

      tol = 1.d-8

      d=1.d0
      nk=-n
      do 180 k=1,n
         nk=nk+n
         l(k)=k
         m(k)=k
         kk=nk+k
         biga=a(kk)
         do 20 j=k,n
            iz=n*(j-1)
            do 20 i=k,n
               ij=iz+i
               if (dabs(biga)-dabs(a(ij))) 10,20,20
   10          biga=a(ij)
               l(k)=i
               m(k)=j
   20    continue
         j=l(k)
         if (j-k) 50,50,30
   30    ki=k-n
         do 40 i=1,n
            ki=ki+n
            holo=-a(ki)
            ji=ki-k+j
            a(ki)=a(ji)
   40    a(ji)=holo
   50    i=m(k)
         if (i-k) 80,80,60
   60    jp=n*(i-1)
         do 70 j=1,n
            jk=nk+j
            ji=jp+j
            holo=-a(jk)
            a(jk)=a(ji)
   70    a(ji)=holo
   80    if (dabs(biga)-tol) 90,100,100
   90    d=0.d0
         return
  100    do 120 i=1,n
            if (i-k) 110,120,110
  110       ik=nk+i
            a(ik)=a(ik)/(-biga)
  120    continue
         do 150 i=1,n
            ik=nk+i
            ij=i-n
            do 150 j=1,n
               ij=ij+n
               if (i-k) 130,150,130
  130          if (j-k) 140,150,140
  140          kj=ij-i+k
               a(ij)=a(ik)*a(kj)+a(ij)
  150    continue
         kj=k-n
         do 170 j=1,n
            kj=kj+n
            if (j-k) 160,170,160
  160       a(kj)=a(kj)/biga
  170    continue
         d=min(d*biga,1.d10)
         a(kk)=1.d0/biga
  180 continue
      k=n
  190 k=k-1
      if (k) 260,260,200
  200 i=l(k)
      if (i-k) 230,230,210
  210 jq=n*(k-1)
      jr=n*(i-1)
      do 220 j=1,n
         jk=jq+j
         holo=a(jk)
         ji=jr+j
         a(jk)=-a(ji)
  220 a(ji)=holo
  230 j=m(k)
      if (j-k) 190,190,240
  240 ki=k-n
      do 250 i=1,n
         ki=ki+n
         holo=a(ki)
         ji=ki+j-k
         a(ki)=-a(ji)
  250 a(ji)=holo
      goto 190
  260 return
c
      end

      double precision function rad(i,ianz)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      common /esprad/ vander(mxel)
      dimension ianz(*)

      bohr = 0.529177249d+00
      rad = vander(ianz(i))/bohr

      return
      end

      subroutine connld(dens,idomap,isp,
     &                  coo,ianz,iaton,iatclr,iresid,
     &                  iconn,isurf)
      implicit double precision (a-h,o-z)
      parameter (mxptn=5000)
      parameter (mxcon=10)
      parameter (mxel=100)
      parameter (max3d=61)
      parameter (mxesp=max3d*max3d*max3d+max3d) 
      parameter (mxplev=5)
      common /athlp/  iatoms, mxnat
      common /surf/   natorg,nosncd
      common /spa3d/  esp(mxesp)

      common /potlev/ grdw, plevel(mxplev),ipcol(mxplev+1),nplev
      common /esprad/ vander(mxel)
      integer srfmap, srfloc
      common /hlpsrf/ npts1,npts2,npts3,srfmap,srfloc,ifogl,itsrf
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension con(3,mxptn)
      dimension inbr(200),cnbr(3,200),rnbr(200)
      logical mnbr(200)
      logical collid
      dimension ci(3), temp0(3)
      dimension cw(3), jt(2,3),dmn(2,3),ic(3)
      dimension coo(3,*),ianz(*),iaton(*),iatclr(*),iresid(*),
     &          iconn(mxcon+1,*),isurf(*)

      data plevel /-0.1d0,-0.05d0,0.0d0,0.05d0,0.1d0/
      data ipcol /7,5,4,2,9,1/

      itsrf = 0
      nplev = 5
      call setcod(0.05d0)

      call parptr(18,esp,fdum,idum)

      idebug = 1

      pi = 4.d0*datan(1.d0)


      bohr = 0.529177249d+00
      den = dens
      t = 0.0d0

      if (isp.eq.5) then
          nesp = iatoms
          iatoms = natorg
          goto 200
      endif

      if (natorg.eq.0) then
          natorg = iatoms
      else
          iatoms = natorg
      endif
      nesp = iatoms

c      if (isp.eq.3) call mksph3(con,ncon,t)
c
c     big loop for each atom
c
      call inferr(
     &          'Starting generation of surface',0)
      do iatom = 1, iatoms

         if (ianz(iatom).ge.99.or.isurf(iatom).lt.1) goto 110
c
c     transfer values from large arrays to iatom variables
c
         ri = rad(iatom,ianz)
         do k = 1,3
            ci(k) = coo(k,iatom)
         end do
c
c     gather the neighboring atoms of iatom
c
         nnbr = 0
         do jatom = 1, iatoms
            if (ianz(jatom).ge.99.or.isurf(jatom).lt.1) goto 60
            if (iatom .eq. jatom) goto 60
            d2 = dist2(ci,coo(1,jatom))
            if (d2 .ge. (ri+rad(jatom,ianz)) ** 2) goto 60
c
c     we have a new neighbor
c     transfer atom coordinates, radius and surface request number
c
            nnbr = nnbr + 1
            if (nnbr .gt. 200)then
               write(iun3,*) 'too many neighbors:',nnbr
               stop
            endif
            inbr(nnbr) = jatom
            do k = 1,3
               cnbr(k,nnbr) = coo(k,jatom)
            end do
            rnbr(nnbr) = rad(jatom,ianz)
60          continue
         end do
c
c     contact surface
c

         if (isp.eq.1) then
            ncon = (4 * pi * ri ** 2) * den
            if (ncon .gt. mxptn) ncon = mxptn
         elseif (isp.eq.2.or.isp.eq.4) then
            ncon = 1
         endif

c
c     this call may decrease ncon somewhat
c
         if ( ncon .eq. 0) then
            write(iun3,*) 'vector length of zero in connol'
            stop
         endif

         if (isp.eq.1) then
             call mksph1(con,ncon,3)
c         elseif (isp.eq.2) then
c             call mksph2(con,2,ncon,t)
         elseif (isp.eq.4) then
             call mksph4(ci,ri,grdw,con,ncon)
         endif
c
c     contact probe placement loop
c
         nold = nesp
         do i = 1,ncon
            do k = 1,3
               if (isp.eq.4) then
                  cw(k) = con(k,i)
               else
                  cw(k) = ci(k) + ri * con(k,i)
               endif
            end do
c
c     check for collision with neighboring atoms
c
            if (collid(cw,cnbr,rnbr,mnbr,nnbr,1,
     &      jnbr,knbr)) goto 100
            do j=iatoms+1,nold
               if (dist2(coo(1,j),cw).lt.1.0d-5) goto 100
            end do
            do kk=1,3
               if (isp.eq.4) then
                  temp0(kk) = con(kk,i)
               else
                  temp0(kk) = ci(kk)+ri*con(kk,i)
               endif
            end do
c
c     store point in connl and increment nesp
c
            nesp = nesp + 1
            if (nesp .gt. mxnat) then
               call inferr(
     &          'to many points, surface incomplete !',0)
               nesp = nesp - 1
               goto 200
            endif
            coo(1,nesp) = temp0(1)
            coo(2,nesp) = temp0(2)
            coo(3,nesp) = temp0(3)
            ianz(nesp) = 100
            iaton(nesp) = 1
            iresid(nesp) = iatom
            iatclr(nesp) = iatclr(iatom)
            isurf(nesp) = 0
            iconn(1,nesp) = 0
  100       continue
         end do
         if (isp.eq.4) then
            npts = nesp - nold
            call mkscon(ci,ri,grdw,nold,npts)
            goto 110
         endif

         if (isp.eq.3) then
            tl = t*ri
            thresh = tl*tl
         elseif (isp.eq.2) then
            tl = t*ri*1.4d0
            thresh = tl*tl
         else
            thresh = 2.5d0
         endif

         do i=nold+1,nesp
            do j=i+1,nesp
               if (dist2(coo(1,i),coo(1,j))-thresh.le.0.05d0) then
                if (iconn(1,i).lt.mxcon.and.iconn(1,j).lt.mxcon) then
                  nc = iconn(1,i)
                  if (nc.lt.mxcon) then
                     iconn(1,i) = nc + 1
                     iconn(nc+2,i) = j
                  endif
                  nc = iconn(1,j)
                  if (nc.lt.mxcon) then
                     iconn(1,j) = nc + 1
                     iconn(nc+2,j) = i
                  endif
                endif
               endif
            end do
         end do

  110    continue
      end do

  200 continue


      if (idomap.ne.0) then
         emin = 10000.0d0
         emax = -emin
         do i=iatoms+1,nesp
             if (idomap.eq.1) then
                call espot(coo(1,i),coo(2,i),coo(3,i),esp(i),0)
             elseif (idomap.eq.2) then
                call calc(coo(1,i),coo(2,i),coo(3,i),esp(i))
             elseif (idomap.eq.3) then
                call clmon(coo(1,i),esp(i))
             endif
             if (esp(i).lt.emin) emin = esp(i)
             if (esp(i).gt.emax) emax = esp(i)
         end do
         if (idebug.eq.1) then
             print*,' '
             print*,'emin = ',emin,' emax = ',emax
             print*,' '
         endif
         call clrcod(iatoms,nesp,idebug)
      endif

      if (isp.eq.5) goto 2222

      if (isp.eq.4) then
         thresh = grdw*grdw*2.000001d0
      else
         thresh = 2.5d0
      endif

      do i=iatoms+1,nesp
         jlast = 0
         if (isp.eq.4) then
            do l=1,3
               dmn(1,l) = 10000.0d0
               dmn(2,l) = 10000.0d0
               jt(1,l) = 0
               jt(2,l) = 0
            end do
         endif

         icn = 0
         iatom = iresid(i)

         do j=iatoms+1,nesp

               jatom = iresid(j)

               if (jatom.ne.jlast) then
                   jlast = jatom
                   icn = 0
                   if (iatom.ne.jatom) then
                      d2 = dist2(coo(1,iatom),coo(1,jatom))
                      if (d2.le.(rad(iatom,ianz)+rad(jatom,ianz))**2) 
     &                   icn = 1
                   endif
               endif

               if (icn.eq.1) then
                  if (isp.eq.4) then
                     ii = 0
                     do l=1,3
                        if (coo(l,i).eq.coo(l,j)) then
                           ii = l
                        endif
                     end do
                  else
                     ii = 1
                  endif
                  dd = dist2(coo(1,i),coo(1,j))
                  if (dd.lt.thresh.and.ii.ne.0.and.dd.ne.0.0d0) then
                     if (isp.eq.4) then
                        if (dd.lt.dmn(1,ii)) then
                           dmn(2,ii) = dmn(1,ii)
                           jt(2,ii) = jt(1,ii)
                           dmn(1,ii) = dd
                           jt(1,ii) = j
                        elseif (dd.lt.dmn(2,ii)) then
                           dmn(2,ii) = dd
                           jt(2,ii) = j
                        endif
                     else
                        nc = iconn(1,i)
                        if (nc.lt.mxcon) then
                           iconn(1,i) = nc + 1
                           iconn(nc+2,i) = j
                        endif
                        nc = iconn(1,j)
                        if (nc.lt.mxcon) then
                           iconn(1,j) = nc + 1
                           iconn(nc+2,j) = i
                        endif
                     endif
                  endif
               endif
         end do
         if (isp.eq.4) then
            do l=1,3
               ic(l) = 0
            end do
            nc = iconn(1,i)
            do j=1,nc
               do l=1,3
                  if (coo(l,i).eq.coo(l,iconn(1+j,i))) then
                     ic(l) = ic(l) + 1
                  endif
               end do
            end do
            do l=1,3
              ii = 0
              if (jt(1,l).ne.0.and.ic(l).le.1) then
                 ii = 1
                 if (jt(2,l).ne.0.and.ic(l).eq.0) ii = 2
              endif
              do j=1,ii
                 nc = iconn(1,i)
                 if (nc.lt.mxcon) then
                    iconn(1,i) = nc + 1
                    iconn(nc+2,i) = jt(j,l)
                 endif
              end do
            end do
         endif
      end do

2222  continue

      do i=iatoms+1,nesp
         iresid(i) = -4
      end do

      iatoms = nesp

      call inferr(
     &          'Surface generation completed !',0)

      return
      end

      subroutine clrsrd(isurf)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      common /surf/ natorg,nosncd
      integer srfmap, srfloc
      common /hlpsrf/ npts1,npts2,npts3,srfmap,srfloc,ifogl,itsrf
      dimension isurf(*)

      if (natorg.ne.0) iatoms = natorg

      do i=1,iatoms
          isurf(i) = 0
      end do
      
      if (ifogl.eq.1) call srfclr
         
      return
      end

      subroutine csrft(isurf)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      common /surf/ natorg,nosncd
      dimension isurf(*)

      if (natorg.ne.0) iatoms = natorg

      do i=1,iatoms
          isurf(i) = 0
      end do
      
      return
      end

      subroutine allsrd(idocol,idomap,idocal,isp,
     &                  isurf,iaton,ianz,iatclr)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      common /surf/ natorg,nosncd
      integer srfmap, srfloc
      common /hlpsrf/ npts1,npts2,npts3,srfmap,srfloc,ifogl,itsrf
      dimension isurf(*),iaton(*),ianz(*),iatclr(*)

      do i=1,iatoms
          if (isurf(i).eq.0.and.iaton(i).gt.0
     &        .and.ianz(i).lt.99) then
             isurf(i) = 1
             if (idocol.eq.1) iatclr(i) = 15
          endif
      end do
      
      if (ifogl.eq.1) then
         srfloc = 0
         call asurf(idomap,idocal)
      else
         call connlp(1.0d0,idomap,isp)
      endif

      return
      end

      subroutine alasrd(iresid,isurf,iams,ihets)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      dimension iresid(*),isurf(*),iams(*),ihets(*)

      do j=1,iatoms
         if (iresid(j).gt.0) then
            if (iams(iresid(j)).eq.1) then
               isurf(j) = 1
            endif
         elseif (iresid(j).le.-4) then
            if (ihets(iabs(iresid(j))+1).eq.1) then
               isurf(j) = 1
            endif
         endif
      end do

      return
      end

      subroutine clrcdd(natorg,natoms,idebug,iatclr)
      implicit double precision (a-h,o-z)
      parameter (max3d=61)
      parameter (mxesp=max3d*max3d*max3d+max3d) 
      parameter (mxplev=5)
      common /potlev/ grdw, plevel(mxplev),ipcol(mxplev+1),nplev
      common /spa3d/  esp(mxesp)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension iatclr(*), iespst(mxplev+1)

      do i=1,nplev+1
         iespst(i) = 0
      end do

      do i=natorg+1,natoms
          icolor = ipcol(1)
          ilev = 1
          do j=1,nplev
             if (esp(i).gt.plevel(j)) then
                icolor = ipcol(j+1)
                ilev = j+1
             endif
          end do
          iespst(ilev) = iespst(ilev) + 1
          iatclr(i) = icolor
      end do

      if (idebug.eq.1) then
          write(iun3,'(''-----------------------------------'')')
          write(iun3,'('' 1              esp < '',f10.5,'' '',i5)')
     &          plevel(1),iespst(1)
          do i=1,nplev-1
           write(iun3,'(i2,'' '',f10.5,'' < esp < '',
     &     f10.5,'' '',i5)') i+1,plevel(i),plevel(i+1),iespst(i+1)
          end do
          i = nplev+1
          write(iun3,'(i2,''              esp > '',f10.5,'' '',i5)')
     &       i,plevel(i-1),iespst(i)
          write(iun3,'(''-----------------------------------'')')
      endif

      return
      end

      subroutine setcod(spacin)
      implicit double precision (a-h,o-z)
      parameter (mxplev=5)
      common /potlev/ grdw, plevel(mxplev),ipcol(mxplev+1),nplev
      common /mapcls/ valcol(5),colmap(3,5)

      valcol(1) = -2.0d0*spacin
      valcol(2) = -1.0d0*spacin
      valcol(3) = 0.0d0
      valcol(4) = spacin
      valcol(5) = 2.0d0*spacin

c      print*,'valcol=',(valcol(i),i=1,5)

      nlev = nplev/2
      if (nplev-nlev*2.ne.1) return

      plevel(nlev+1) = 0.0d0
      do i=1,nlev
          plevel(i) = -(nlev+1-i)*spacin
          plevel(nlev+1+i) = i*spacin
      end do

      return
      end

      subroutine aresp(nesp,ichrg)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)
      parameter (max3d=61)
      parameter (mesp = (max3d*max3d*max3d+max3d)/4)
      parameter (mdum = (max3d*max3d*max3d+max3d) - mesp*4)
      common /spa3d/  connl(3,mesp), esp(mesp), dum(mdum)

      open(unit=50,form='formatted',file='resp.in',status='unknown',
     &     err=1000)
      write(50,'(a)') 'Molden generated aresp inputfile'
      write(50,'(a)') ' &cntrl nmol=1, ihfree=1'
      write(50,'(a)') ' &end'
      write(50,'(a)') '1.0'
      write(50,'(a)') 'Subtitle'
      write(50,'(2i5)') ichrg,natoms
      izero = 0
      do i=1,natoms
          write(50,'(2i5)') nat(i),izero
      end do
      write(50,'(a)') ' '
      close(50)

      open(unit=50,form='formatted',file='esp.in',status='unknown',
     &     err=1000)
      write(50,'(2i5)') natoms,nesp
      do i=1,natoms
          write(50,'(17x,3E16.7)') (xyz(j,i),j=1,3)
      end do
      do i=1,nesp
          write(50,'(1x,4E16.7)') esp(i),(connl(j,i),j=1,3)
      end do
      write(50,'(a)') ' '
      close(50)

      return
1000  call inferr('Could open file',0)
      return
      end

      subroutine wrxyd(jmode,q)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxel=100)
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension q(*)

      toang  = 0.52917706d0
      
      open(unit=46,form='formatted',file='esp.xyz',
     &    status='unknown',err=2000)

      write(46,'(i5)') natoms
      if (jmode.eq.1) then
         write(46,'(a)') 'Molden DMA esp fitted charges'
      elseif (jmode.eq.2) then
         write(46,'(a)') 'Molden Mulliken charges'
      else
         write(46,'(a)') 'Molden esp fitted charges'
      endif

      do i=1,natoms
          write(46,'(a2,1x,3(f12.6,1x),f9.6)') 
     &        elemnt(nat(i)),(xyz(j,i)*toang,j=1,3),q(i)
      end do

      close(46)

      write(iun3,*) 'Wrote XYZ file esp.xyz'
      return

2000  write(iun3,*) 'Couldnt write XYZ file esp.xyz'
      return
      end

      subroutine chadd(natoms)
c
c     Reads information from file 48 and adds virtual coordinates
c
      parameter (numatm=2000)
      implicit double precision (a-h,o-z)
      common /coord / xyz(3,numatm)
      character*137 line,str,tstr
      common /curlin/ line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      integer getlin
      double precision xyz
      dimension rba(3),rca(3),r(3)

      iunt = iun2
      iun2 = 48

      do while (getlin(0).eq.1)

c        Add a new virtual atom. Note: distances are atomic units!

         p = 0.0d0
         q = 0.0d0
         s = 0.0d0

         natoms = natoms + 1

         tstr = line
         ivatyp = 0

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            ia = itype
         else
            ivatyp = -1
            goto 200
         endif

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            ib = itype
         else
            ivatyp = -1
            goto 200
         endif

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.3) then
            p = rtype
         else
            ivatyp = -1
            goto 200
         endif

         do m=1,3
            rba(m) = xyz(m,ib) - xyz(m,ia)
         end do

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            ic = itype
            ivatyp = 1
         elseif (ktype.eq.0) then
            goto 200
         else
            ivatyp = -1
            goto 200
         endif

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.3) then
            q = rtype
         elseif (ktype.eq.0) then
            goto 200
         else
            ivatyp = -1
            goto 200
         endif

         do m=1,3
            rca(m) = xyz(m,ic) - xyz(m,ia)
         end do

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.3) then
            s = rtype
            ivatyp = 2
         elseif (ktype.eq.0) then
            goto 200
         else
            ivatyp = -1
            goto 200
         endif

200      if (ivatyp.eq.-1) then

            print*,' Incorrect line on unit :'
            print*,tstr
            stop

         elseif (ivatyp.eq.0) then

c           Virtual atom in a bond

            do m=1,3
               xyz(m,natoms) = xyz(m,ia) + p*rba(m)
            end do

            write(iun3,20) natoms,(xyz(m,natoms),m=1,3),ia,ib,p

         elseif (ivatyp.eq.1) then

c           Virtual atom in a plane

            do m=1,3
               xyz(m,natoms) = xyz(m,ia) + p*rba(m) + q*rca(m)
            end do

            write(iun3,20) natoms,(xyz(m,natoms),m=1,3),ia,ib,p,ic,q

         elseif (ivatyp.eq.2) then

c           Virtual atom in three dimensions

            s = s*0.52917706d0

c           s krijgt (tijdelijk?) de dimensie A^(-1)

            r(1) = rba(2)*rca(3) - rba(3)*rca(2)
            r(2) = rba(3)*rca(1) - rba(1)*rca(3)
            r(3) = rba(1)*rca(2) - rba(2)*rca(1)

            d2 = 0

            do m=1,3
               d = p*rba(m) + q*rca(m) + s*r(m)
               xyz(m,natoms) = xyz(m,ia) + d
               d2 = d2 + d*d
            end do

            write(iun3,20) natoms,(xyz(m,natoms),m=1,3),ia,ib,p,ic,q,s,
     &           sqrt(d2)*0.52917706d0

c               print*,(rba(m),m=1,3),sqrt(rba(1)**2+rba(2)**2+rba(3)**2)
c               print*,(rca(m),m=1,3),sqrt(rca(1)**2+rca(2)**2+rca(3)**2)
c               print*,(r(m),m=1,3),sqrt(r(1)**2+r(2)**2+r(3)**2)

         endif

      end do

      iun2 = iunt

   20 format(' New atom',i4,3f10.5,2i4,f8.3,i4,2f8.3,f12.3)

      return
      end


