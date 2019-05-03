      subroutine crpsin(a,b,c,d)
c
c       calculates cross product:  (b-a) x (c-a) = d
c                                   ---     ---    -
      implicit real (a-h,o-z)
      dimension a(3),b(3),c(3),d(3)
      dimension v1(3),v2(3)
      do 10 i=1,3
         v1(i)=b(i)-a(i)
         v2(i)=c(i)-a(i)
10    continue
c
      d(1)=v2(2)*v1(3)-v2(3)*v1(2)
      d(2)=v2(3)*v1(1)-v2(1)*v1(3)
      d(3)=v2(1)*v1(2)-v2(2)*v1(1)
      dlen=sqrt(d(1)**2+d(2)**2+d(3)**2)
      do 20 i=1,3
20       d(i)=d(i)/dlen
      return
      end
