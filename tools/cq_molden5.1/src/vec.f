      subroutine vec(small,ohoh,u,c,j,k)
c
c
      implicit double precision (a-h,o-z)
      logical ohoh
c
      dimension c(*),r(3),u(3)
c
      data dzero/0.0d0/
c
      jtemp=(j-1)*3
      ktemp=(k-1)*3
      r2=dzero
      do 10 i=1,3
      r(i)=c(i+jtemp)-c(i+ktemp)
      r2=r2+r(i)*r(i)
 10   continue
      r2=dsqrt(r2)
      ohoh=r2.lt.small
      if(ohoh)return
      do 20 i=1,3
 20   u(i)=r(i)/r2
      return
c
      end
