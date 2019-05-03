csg      subroutine crpsin(a,b,c,d)
csgc
csgc       calculates cross product:  (b-a) x (c-a) = d
csgc                                   ---     ---    -
csg      implicit real (a-h,o-z)
csg      dimension a(3),b(3),c(3),d(3)
csg      dimension v1(3),v2(3)
csg      do 10 i=1,3
csg         v1(i)=b(i)-a(i)
csg         v2(i)=c(i)-a(i)
csg10    continue
csgc
csg      d(1)=v2(2)*v1(3)-v2(3)*v1(2)
csg      d(2)=v2(3)*v1(1)-v2(1)*v1(3)
csg      d(3)=v2(1)*v1(2)-v2(2)*v1(1)
csg      dlen=sqrt(d(1)**2+d(2)**2+d(3)**2)
csg      do 20 i=1,3
csg20       d(i)=d(i)/dlen
csg      return
csg      end
csg      subroutine sgdisp
csg#include <gl/fgl.h>
csg#include <gl/fdevice.h>
csg      integer*4 iang,jang,kang
csg      dimension amat(15),alt(10),al(12),bmat(15),umat(4,4)
csg      data amat/ambien,.1,.1,.1,diffus,.200,.369,.065,
csg     &          specul,.5,.5,.5,shinin,30.0,lmnull/
csg      data bmat/ambien,.0,.0,.0,diffus,.0,.0,.0,
csg     &          specul,.0,.0,.0,shinin,.0,lmnull/
csg      data umat/1.0,0.0,0.0,0.0,
csg     &          0.0,1.0,0.0,0.0,
csg     &          0.0,0.0,1.0,0.0,
csg     &          0.0,0.0,0.0,1.0/
csg      data alt/lcolor,.8,1.0,.8,positi,2.,2.,2.,1.,lmnull/
csg      data al/ambien,.4,.1,.1,localv,0,twosid,1,attenu,.5,.1,lmnull/
csgc      data al/ambien,.1,.1,.1,localv,0.,twosid,0.,lmnull/
csg      call prefpo(0,1500,0,1200)
csg      call winope("molden",6)
csg      call zbuffe(.true.)
csg      call double
csg      call rgbmod
csg      call gconfi
csg      call mmode(mviewi)
csg      call loadma(umat)
csg      call ortho(-40.0,40.0,-40.0,40.0,-40.0,40.0)
csg      call lmdef(defmat,1,0,amat)
csg      call lmdef(deflig,2,0,alt)
csg      call lmdef(deflmo,3,0,al)
csg      call lmbind(materi,1)
csgc      call lmbind(backma,1)
csg      call lmbind(light0,2)
csg      call lmbind(lmodel,3)
csg      iang=0
csg      jang=0
csg      kang=0
csg      do 100 i=0,1000
csg         call pushma 
csg         call rotate(iang,'x')
csg         call rotate(jang,'y')
csg         call rotate(kang,'z')
csg         iang=iang+10
csg         jang=jang+13
csg         if(iang+jang.gt.3000)kang=kang+17
csg         if(iang.gt.3600) iang=iang-3600
csg         if(jang.gt.3600) jang=jang-3600
csg         if(kang.gt.3600) kang=kang-3600
csg         call cpack($0)
csg         call clear
csg         call zclear
csg         call drawsf
csg         call swapbu
csg         call popmat
csg100   continue
csg      return
csg      end
csg      subroutine drawsf
csg#include <gl/fgl.h>
csg#include <gl/fdevice.h>
csg      common /sgcon/ vec1(3,12800),vec2(3,12800),vec3(3,12800),
csg     &          vnorm(3,12800),npol
csg      do 100 i=1,npol
csg      call bgnpol
csg         call n3f(vnorm(1,i))
csg         call v3f(vec1(1,i))
csg         call v3f(vec2(1,i))
csg         call v3f(vec3(1,i))
csg      call endpol
csg100   continue
csg      return
csg      end
