      subroutine silld(const,npts,dens)
      implicit real (a-h,o-z)

c THIS IS REALLY silly

      double precision dens,const
      common /sgcon/ vec1(3,12800),vec2(3,12800),vec3(3,12800),
     &          vnorm(3,12800),npol
      dimension dens(*)

c     a grid of n*n points has n-1*n-1 squares
c     and twice as much triangular polygons
c     so 2*(npts-1)**2
c     npts is 80 at maximum so lets make it 12800
      ij=0
      noff=npts/2+1
      rpts=float(npts-1)
      do 100 i=1,npts-1
      do 100 j=1,npts-1
         ij=ij+1
c        first triangle
         vec1(1,ij)=float(j-noff)
         vec1(2,ij)=float(i-noff)
         vec1(3,ij)=sngl(dens(j+(i-1)*npts)*const)*rpts
         vec2(1,ij)=float(j+1-noff)
         vec2(2,ij)=float(i-noff)
         vec2(3,ij)=sngl(dens(j+1+(i-1)*npts)*const)*rpts
         vec3(1,ij)=float(j-noff)
         vec3(2,ij)=float(i+1-noff)
         vec3(3,ij)=sngl(dens(j+i*npts)*const)*rpts
         call crpsin(vec1(1,ij),vec2(1,ij),vec3(1,ij),vnorm(1,ij))
c        second triangle
         ij=ij+1
         vec1(1,ij)=float(j+1-noff)
         vec1(2,ij)=float(i-noff)
         vec1(3,ij)=sngl(dens(j+1+(i-1)*npts)*const)*rpts
         vec2(1,ij)=float(j+1-noff)
         vec2(2,ij)=float(i+1-noff)
         vec2(3,ij)=sngl(dens(j+1+i*npts)*const)*rpts
         vec3(1,ij)=float(j-noff)
         vec3(2,ij)=float(i+1-noff)
         vec3(3,ij)=sngl(dens(j+i*npts)*const)*rpts
         call crpsin(vec1(1,ij),vec2(1,ij),vec3(1,ij),vnorm(1,ij))
100   continue
      npol=ij
      call sgdisp
      return
      end
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
      subroutine sgdisp
#include <gl/fgl.h>
#include <gl/fdevice.h>
      integer*4 iang,jang,kang
      dimension amat(15),alt(10),al(12),bmat(15),umat(4,4)
      data amat/ambien,.1,.1,.1,diffus,.200,.369,.065,
     &          specul,.5,.5,.5,shinin,30.0,lmnull/
      data bmat/ambien,.0,.0,.0,diffus,.0,.0,.0,
     &          specul,.0,.0,.0,shinin,.0,lmnull/
      data umat/1.0,0.0,0.0,0.0,
     &          0.0,1.0,0.0,0.0,
     &          0.0,0.0,1.0,0.0,
     &          0.0,0.0,0.0,1.0/
      data alt/lcolor,.8,1.0,.8,positi,2.,2.,2.,1.,lmnull/
      data al/ambien,.4,.1,.1,localv,0,twosid,1,attenu,.5,.1,lmnull/
c      data al/ambien,.1,.1,.1,localv,0.,twosid,0.,lmnull/
      call prefpo(0,1500,0,1200)
      call winope("molden",6)
      call zbuffe(.true.)
      call double
      call rgbmod
      call gconfi
      call mmode(mviewi)
      call loadma(umat)
      call ortho(-40.0,40.0,-40.0,40.0,-40.0,40.0)
      call lmdef(defmat,1,0,amat)
      call lmdef(deflig,2,0,alt)
      call lmdef(deflmo,3,0,al)
      call lmbind(materi,1)
c      call lmbind(backma,1)
      call lmbind(light0,2)
      call lmbind(lmodel,3)
      iang=0
      jang=0
      kang=0
      do 100 i=0,1000
         call pushma 
         call rotate(iang,'x')
         call rotate(jang,'y')
         call rotate(kang,'z')
         iang=iang+10
         jang=jang+13
         if(iang+jang.gt.3000)kang=kang+17
         if(iang.gt.3600) iang=iang-3600
         if(jang.gt.3600) jang=jang-3600
         if(kang.gt.3600) kang=kang-3600
         call cpack($0)
         call clear
         call zclear
         call drawsf
         call swapbu
         call popmat
100   continue
      return
      end
      subroutine drawsf
#include <gl/fgl.h>
#include <gl/fdevice.h>
      common /sgcon/ vec1(3,12800),vec2(3,12800),vec3(3,12800),
     &          vnorm(3,12800),npol
      do 100 i=1,npol
      call bgnpol
         call n3f(vnorm(1,i))
         call v3f(vec1(1,i))
         call v3f(vec2(1,i))
         call v3f(vec3(1,i))
      call endpol
100   continue
      return
      end
