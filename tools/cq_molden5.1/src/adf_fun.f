      subroutine adffud(x,y,z,psi,
     &                  stoalfa,stobnorm,istos,naorbs)
c     THIS iS REALLY adffun
c
c     data description
c
c     Integer for SLater functions in ADF
c     istos(5,maxao)      (1) = p orbital is centered on atom (1,p)
c                         (2) = x exponent of p orbtal
c                         (3) = y exponent of p orbtal
c                         (4) = z exponent of p orbtal
c                         (5) = r exponent of p orbtal
c
c     Double for SLater functions in ADF
c     dsladf(2,maxao)     (1) = alfa axponent of p orbital
c                         (2) = normalzzation factor for p orbital
c
c     Double Atomic to Molecular ADF coefficients
c     damadf(maxmo,maxao) LCAO coefficients
c
c
c     xyz(3,numatm)       Atomic coordinates as defined il molden source.
c
c     Number Of SLater functions
c     nosf                number of total slater function defined.
c
c     Number Of Molecular Orbitals
c     nomo                number of defined molecular orbitals
c                         In molden is norbs.
c
c     Number Of total AToms
c     noat                number of total atoms. In molden is natoms.
c
c     Comment:
c             These variables must be seved in a common block and this can
c             be a dedidcated ADF common block or one of the just used
c             common block. Using an old common block we may use the same
c             variables or we can define a new sets of variables with the
c             same memory space. In this last case the portability of source
c             code must be evaluated.
c             Now this function make a loop over all defined molecular
c             orbitals. Code for select the MO to be use in a grid must be
c             write.
c

c
c     declarations
c
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /orbhlp/ mxorb,iuhf,ispd
      common/coord/ xyz(3,numatm)
      common /moldat/ natoms,norbs,nelecs,nat(numatm) 
      real stoalfa,stobnorm
      dimension istos(5,*),stoalfa(*),stobnorm(*)
      dimension psi(*)
c
c     local declarations
c
c     reset psi
c
      do i=1,mxorb
         psi(i) = 0.0d0
      end do


c
c     loop over atomic orbitals
c
         do iao = 1,naorbs
c
c     detect distance from right atom
c
               xk = x - xyz(1,istos(1,iao))
               yk = y - xyz(2,istos(1,iao))
               zk = z - xyz(3,istos(1,iao))
               r2 = xk*xk + yk*yk + zk*zk
               r2 = r2 + 1.d-10
               r1 = dsqrt(r2)


c
c     evaluate function
c
               adffunct = adfpsi(r1,xk,yk,zk,
     &              istos(2,iao),istos(3,iao),
     &              istos(4,iao),istos(5,iao),
     &              stoalfa(iao),stobnorm(iao))
               psi(iao) = adffunct


c            endif
         
         end do

      return
      end

      double precision function adfpsi(r,x,y,z,ia,ib,ic,id,da,dn)
c     
c     variables descriptions:
c     r         atom-point distance
c     x,y,z     point coordinates
c     ia,ib,ic  x,y,z exponents
c     id        r exponent
c     da,dn     alfa and normalizzation factor
c
c     function:
c     psi = dn * x^ia * y^ib * z^ic * r^id * exp(da*r)
c     
      parameter (small=1.0d-10, big=30., zero=0.0d0)
      integer ia,ib,ic,id
      real da,dn
      double precision r,x,y,z,rho,desp
      integer l,n
c
c
      l = ia + ib + ic
      n = l + id + 1
c
      if (dabs(x).lt.small) x = zero
      if (dabs(y).lt.small) y = zero
      if (dabs(z).lt.small) z = zero

c
c     the point is on the nucleus      
c
      if (r.lt.small) then
         if (l.eq.0.and.n.eq.1) then
            adfpsi = 1.0d0
         else
            adfpsi = zero
         endif
         return
      endif
c
c     general case
c
      rho = da * r
      if (rho.lt.big) then
         adfpsi = dn * desp(x,ia) * desp(y,ib) * desp(z,ic) * 
     1        desp(r,id) * dexp(-rho)
      else
         adfpsi = zero
      endif

      return
      end


      double precision function desp(x,a)
      integer a
      double precision x

      if (a.eq.0) then
         desp = 1.0d0
      else 
         desp = x**a
      endif

      return
      end

