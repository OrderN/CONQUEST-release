      subroutine ryspol(n,x,t,w)
      implicit double precision (a-h,o-z)
      dimension t(*),w(*)
      dimension f(10,10),c(10,10),ra(20),dr(20),ap(20)
      dimension rap(20),fmt(20),a(100)

      if (n.le.0) stop 'ryspol n illegal'

      call fmtgen(fmt,x,2*n+1,istatus)
      if (istatus.ne.0) stop 'ryspol in fmtgen'

      do j=1,n+1
         do i=1,n+1
            f(i,j) = fmt(i+j-1)
         end do
      end do

c start with identity matrix

      do i = 1,n+1
         do j = 1,n+1
            c(i,j) = 0.0d0
         end do
         c(i,i) = 1.0d0
      end do
 
      c(1,1) = 1.0d0 / dsqrt(f(1,1))
      temp = 1.0d0 / (dsqrt(f(1,1)*(f(1,1)*f(2,2) - f(2,1)*f(2,1))))
      c(1,2) =  -f(2,1)*temp
      c(2,2) =   f(1,1)*temp

      if (n+1.gt.100) stop 'ryspol nbasis.gt.100'
      if (n+1.le.1) goto 100
 
      do k=1,n + 1
         if (k.ne.1) then
            do i=1,n + 1
               a(i) = 0.0d0
               do j=1,n+1
                  a(i) = a(i) + c(j,k)*f(j,i)
               end do
            end do
            do j = 1,k - 1
               b = 0.0d0
               do i=1,n+1
                  b = b + c(i,j)*a(i)
               end do
               do i=1,n+1
                  c(i,k) = c(i,k) - b*c(i,j)
               end do
            end do
         endif
         tt = 0.0d0
         do i=1,n+1
            t1 = 0.0d0
            do j=1,n+1
               t1 = t1 + c(j,k)*f(j,i)
            end do
            tt = tt + c(i,k)*t1
         end do
         tt = 1.0d0 / dsqrt(tt)
         do i=1,n+1
            c(i,k) = c(i,k)*tt
         end do

      end do

100   continue

      if (n.le.1) then
         t(1) = dsqrt( - c(1,2)/c(2,2))
         w(1) = f(1,1)
         goto 300
      endif

      k = 2*n+1
      do i = 1,n+1
         ra(k) = c(i,n+1)
         k = k - 2
      end do
      ndeg = 2*n

      call rysrot(ndeg,ra,dr,rap)
      call fdpol(ndeg,nfd,a,ap)

      do i = 1,n
         root = dr(i)
         do j = 1,50
            f1 = fpol(ndeg,root,a)
            f2 = fpol(nfd,root,ap)
            if (dabs(f2).lt.1.0d-12) goto 200
            if (dabs(f1/(f2*root)).le.1.0d-12) goto 200
            root = root - f1 / f2
         end do
200      t(i) = root
      end do

      do i = 1,n
         w(i) = 0.0d0
         root = t(i)
         do j = 1,n
            sum = 0.0d0
            do k = 1,j
               km1 = k - 1
               sum = sum + (c(k,j)*(root**(2*km1)))
            end do
            w(i) = w(i) + sum*sum
         end do
         w(i) = 1.0d0/w(i)
      end do

300   continue

      do i = 1,n
         t(i) = t(i)*t(i)
      end do

      return
      end
      double precision function fpol(npol,x,f)
      implicit double precision (a-h,o-z)
      dimension f(*)
      fpol = 0.0d0
      if (npol.eq.0) fpol = f(1)
      if (npol.lt.1) return
      sum = f(npol+1)
      xx = x
      do i = 1, npol
         sum = sum + f(npol-i+1) * xx
         xx = xx*x
      end do
      fpol = sum

      return
      end
      subroutine fdpol(npow,nwpow,f,g)
      implicit double precision (a-h,o-z)
      dimension f(*), g(*)
 
      nwpow = npow - 1
      if (npow.lt.1) return
      do i = 1, npow
          g(i) = f(i) * dfloat(npow-i+1)
      end do

      return
      end
