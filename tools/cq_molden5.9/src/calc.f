      subroutine calc(x,y,z,pot)
      implicit double precision (a-h,o-z), integer ( i-n)
      parameter (mxsite=300)
      parameter (lmax=4)
      parameter (lmaxf=(lmax+1)*(lmax+1))
      character*8   ctag
      common /multip/ q(lmaxf,mxsite),car(3,mxsite),ctag(mxsite),
     &                nsites
      common /extchg/ exchg(3,3),iexchg(3),nexchg
      dimension pvec(3)
c in this procedure the calculation order is inversed
c we use
c 
c   a      b      c      d       
c ---- + ---- + ---- + ---- = ((((d /x**2+c) /x**2)+b) /x**2) +a) /x
c   x    x**3   x**5   x**7   
c
      rt3  =        dsqrt(3.0d0)
      rt6  = 0.25d0*dsqrt(6.0d0)
      rt10 = 0.25d0*dsqrt(10.0d0)
      rt15 = 0.5d0 *dsqrt(15.0d0)

      pot = 0.0d0 
      do i=1,nsites

         pvec(1) = x-car(1,i)
         pvec(2) = y-car(2,i)
         pvec(3) = z-car(3,i)
         qp1     = pvec(1)*pvec(1)
         qp2     = pvec(2)*pvec(2)
         qp3     = pvec(3)*pvec(3)
         r2      = qp1+qp2+qp3
         r       = dsqrt(r2)
         pot     = pot 
     & +((((
     &(q(10,i)*0.5d0*(5.0d0*pvec(3)*qp3-3.0d0*r2*pvec(3))+
     & q(11,i)*rt6*(5.0d0*pvec(1)*qp3-r2*pvec(1))+
     & q(12,i)*rt6*(5.0d0*pvec(2)*qp3-r2*pvec(2))+
     & q(13,i)*rt15*pvec(3)*(qp1-qp2)+
     & q(14,i)*2.0d0*rt15*pvec(1)*pvec(2)*pvec(3)+
     & q(15,i)*rt10*pvec(1)*(qp1-3.0d0*qp2)+
     & q(16,i)*rt10*(-pvec(2))*(qp2-3.0d0*qp1)))/r2
     & +(q(5,i)*0.5d0*(3.0d0*qp3-r2)+
     & q(6,i)*rt3*pvec(1)*pvec(3)+
     & q(7,i)*rt3*pvec(2)*pvec(3)+
     & q(8,i)*0.5d0*rt3*(qp1-qp2)+
     & q(9,i)*rt3*pvec(1)*pvec(2)))/r2
     & +(q(2,i)*pvec(3)+q(3,i)*pvec(1)+q(4,i)*pvec(2)))/r2
     & +q(1,i))/r
      end do


      if (nexchg.ne.0) then
          do i=1,nexchg
             xt = exchg(1,i) - x
             yt = exchg(2,i) - y
             zt = exchg(3,i) - z
             rsq = xt*xt + yt*yt + zt*zt
             if (rsq.ge.1.0d-8) then
                if (iexchg(i).eq.1) then
                   pot = pot + dsqrt(1.0d0 / rsq)
                else
                   pot = pot - dsqrt(1.0d0 / rsq)
                endif
             endif
          end do
      endif

      return
      end 

      subroutine clmod(car,pot,coo,qat)
      implicit double precision (a-h,o-z), integer( i-n)
      common /athlp/ iatoms, mxnat
      dimension car(3),coo(3,*),qat(*)

      pot = 0.0d0
      do i=1,iatoms
         r = dsqrt(dist2(car,coo(1,i)))
         pot = pot + qat(i)/r
      end do

      return
      end 

      subroutine clmond(car,pot,idoloc,coo,qat,icont,ncont)
      implicit double precision (a-h,o-z), integer( i-n)
      common /athlp/ iatoms, mxnat
      dimension car(3),coo(3,*),qat(*),icont(*)

      cutoff = 900.0d0
      pot = 0.0d0

      if (idoloc.eq.1) then
         do l=1,ncont
            i = icont(l)
            if (i.gt.0.and.i.le.mxnat) then
               r2 = dist2(car,coo(1,i))
               pot = pot + qat(i)/r2
            endif
         end do
      else
         do i=1,iatoms
            r2 = dist2(car,coo(1,i))
            r = dsqrt(r2)
            pot = pot + qat(i)/r
         end do
      endif

      return
      end 

