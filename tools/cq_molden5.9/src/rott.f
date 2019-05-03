      subroutine rotd(x,y,z,xc,yc,zc,itran,
     &                rx,ry,rz,t)
      implicit double precision (a-h,o-z)
      dimension rx(3),ry(3),rz(3),t(3)

      if (itran.eq.1) then
         xc = (x-t(1))*rx(1)+(y-t(2))*rx(2)+(z-t(3))*rx(3)
         yc = (x-t(1))*ry(1)+(y-t(2))*ry(2)+(z-t(3))*ry(3)
         zc = (x-t(1))*rz(1)+(y-t(2))*rz(2)+(z-t(3))*rz(3)
      else
         xc = x*rx(1)+y*rx(2)+z*rx(3)
         yc = x*ry(1)+y*ry(2)+z*ry(3)
         zc = x*rz(1)+y*rz(2)+z*rz(3)
      endif

      return
      end

      subroutine rottd(x,y,z,xc,yc,zc,itran,
     &                 rx,ry,rz,t)
      implicit double precision (a-h,o-z)
      real x,y,z
      dimension rx(3),ry(3),rz(3),t(3)

      if (itran.eq.1) then
         xc = (x-t(1))*rx(1)+(y-t(2))*rx(2)+(z-t(3))*rx(3)
         yc = (x-t(1))*ry(1)+(y-t(2))*ry(2)+(z-t(3))*ry(3)
         zc = (x-t(1))*rz(1)+(y-t(2))*rz(2)+(z-t(3))*rz(3)
      else
         xc = x*rx(1)+y*rx(2)+z*rx(3)
         yc = x*ry(1)+y*ry(2)+z*ry(3)
         zc = x*rz(1)+y*rz(2)+z*rz(3)
      endif

      return
      end

      subroutine mktrd(inct,incp,
     &                 xv,yv,zv,pincr,
     &                 scal,scali,smag)
      implicit double precision (a-h,o-z)
      parameter (ssincr=1.05d0)
      parameter (xincr=0.10d0)
      parameter (yincr=0.10d0)
      parameter (rincr=1.0d0)
      common /cllmat/rx(3),ry(3),rz(3),t(3),tz(3),tzorg(3),itz
      common /coars/ icoars
      integer dolabs,fancy,persp,shade,atcol,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt
      dimension vec(3)

      celinc = 0.001d0*0.52917706d0

      if (irtcel.ne.0) then
         do i=1,3
            vec(i) = 0.0d0
         end do
      endif

      idir = 0
      if (inct.eq.290) then
         if (irtcel.ne.0) then
            idir = incp
            iax = 3
         else
            if (persp.eq.1) then
               if (incp.eq.1)  zv = zv + pincr
               if (incp.eq.-1) zv = zv - pincr
            else
               if (incp.eq.1)  smag = smag*ssincr
               if (incp.eq.-1) smag = smag/ssincr
               smag = max(smag,0.01d0)
               scal = scali*2.4d0*smag
            endif
         endif
      elseif (inct.eq.415) then
         if (irtcel.ne.0) then
            idir = 1
            if (incp.lt.0) idir = -1
            iax = 2
         else
            yv = yv + yincr*(5**icoars)*incp
         endif
      elseif (inct.eq.416) then
         if (irtcel.ne.0) then
            idir = 1
            if (incp.lt.0) idir = -1
            iax = 1
         else
            xv = xv + xincr*(5**icoars)*incp
         endif
      elseif (inct.eq.417) then
         if (irtcel.ne.0) then
            idir = 1
            if (incp.lt.0) idir = -1
            iax = 3
         else
            if (persp.eq.1) then
               zv = zv + pincr*incp
            else
               ssc = abs(incp)*ssincr
               if (incp.gt.0) smag = smag*ssc
               if (incp.lt.0) smag = smag/ssc
               smag = max(smag,0.01d0)
               smag = min(smag,1000.0d0)
               scal = scali*2.4d0*smag
            endif
         endif
      elseif (inct.eq.420) then
         if (irtcel.ne.0) then
            idir = -1
            iax = 2
         else
            yv = yv - yincr*(5**icoars)
         endif
      elseif (inct.eq.430) then
         if (irtcel.ne.0) then
            idir = 1
            iax = 2
         else
            yv = yv + yincr*(5**icoars)
         endif
      elseif (inct.eq.440) then
         if (irtcel.ne.0) then
            idir = -1
            iax = 1
         else
            xv = xv - xincr*(5**icoars)
         endif
      elseif (inct.eq.450) then
         if (irtcel.ne.0) then
            idir = 1
            iax = 1
         else
            xv = xv + xincr*(5**icoars)
         endif
      else 
         goto 100
      endif

      if (irtcel.ne.0) then
         if (itz.eq.1.and.iax.eq.3) then
            do i=1,3
               vec(i) =  idir*tz(i)*(10**icoars)
            end do
         else
            vec(iax) =  idir*celinc*(10**icoars)
            if ((inct.eq.415.or.inct.eq.416.or.inct.eq.417)
     &          .and.irtcel.eq.2) then
               do i=1,3
                  vec(i) =  celinc*(10**icoars)*incp*vecs(iax,i)
               end do
            endif
         endif
         if (irtcel.eq.1) then
            call cllrot(vec,0,ifd)
         else
            call alnrot(vec,0)
         endif
      endif

      goto 200

100   continue
      if (abs(incp).eq.1.or.abs(inct).eq.421.or.inct.eq.422) then
         if (abs(inct).eq.421) then
            theinc = -incp
            if (inct.eq.421) then
               inct = -3
            else
               inct = -2
            endif
         elseif(inct.eq.422) then
            theinc = -incp
            inct = -1
         else
            if (icoars.eq.2) then
               theinc = incp*rincr*45.0d0
            elseif (icoars.eq.1) then
               theinc = incp*rincr*5.0d0
            else
               theinc = incp*rincr
            endif
         endif

         if (irtcel.gt.0) then
            if (irtcel.eq.1) then
               call clmat(inct,theinc,1)
               call cllrot(vec,1,ifd)
            else
               if (ialtyp.eq.1) then
                  call clmat(inct,theinc,0)
               else
                  call clmat(inct,theinc,1)
               endif
               call alnrot(vec,1)
            endif
         else
            call xyzrot(inct,theinc)
         endif
      endif

200   if (irtcel.gt.1) then
         if (ialtyp.eq.1) then
            call totpmf(tpmf)
c            print*,'totpmf=',tpmf
            call upsco()
         endif
      endif

      return
      end

      subroutine mtind3(rx,ry,rz)
      implicit double precision (a-h,o-z)
      dimension al(9),alt(9)
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt
      dimension rx(3),ry(3),rz(3)
   
      do i=1,3
         do j=1,3
            if (i.eq.j) then
               vecs(i,j) = 1.0d0
            else
               vecs(i,j) = 0.0d0
            endif
         end do
      end do

      do j=1,3
         al(j)   = rx(j)
         al(3+j) = ry(j)
         al(6+j) = rz(j)
         alt(j)   = rx(j)
         alt(3+j) = ry(j)
         alt(6+j) = rz(j)
      end do

      call matinv(al,3,det)

      do j=1,3
         rx(j) = al(j)
         ry(j) = al(3+j)
         rz(j) = al(6+j)
      end do

      call rott(vecs(1,1),vecs(1,2),vecs(1,3),xc,yc,zc,0)
      vecs(1,1) = xc
      vecs(1,2) = yc
      vecs(1,3) = zc
      call rott(vecs(2,1),vecs(2,2),vecs(2,3),xc,yc,zc,0)
      vecs(2,1) = xc
      vecs(2,2) = yc
      vecs(2,3) = zc
      call rott(vecs(3,1),vecs(3,2),vecs(3,3),xc,yc,zc,0)
      vecs(3,1) = xc
      vecs(3,2) = yc
      vecs(3,3) = zc

c restore original rotation

      do j=1,3
         rx(j) = alt(j)
         ry(j) = alt(3+j)
         rz(j) = alt(6+j)
      end do

      return
      end

      subroutine xyzrod(inct,theang,rx,ry,rz)
      implicit double precision (a-h,o-z)
      dimension rx(3),ry(3),rz(3)

      todeg = 45.0d0 / datan(1.0d0)

      sa = dsin(theang/todeg)
      ca = dcos(theang/todeg)

      if (inct.eq.-3) then
         do i=1,3
            y = ry(i)
            z = rz(i)
            ry(i) = ca*y + sa*z
            rz(i) = ca*z - sa*y
         end do
      elseif (inct.eq.-2) then
         do i=1,3
            x = rx(i)
            z = rz(i)
            rx(i) = ca*x + sa*z
            rz(i) = ca*z - sa*x
         end do
      elseif (inct.eq.-1) then
         do i=1,3
            x = rx(i)
            y = ry(i)
            rx(i) = ca*x - sa*y
            ry(i) = ca*y + sa*x
         end do
      endif

      return
      end

      subroutine inirod(rx,ry,rz,t)
      implicit double precision (a-h,o-z)
      dimension rx(3),ry(3),rz(3),t(3)

      do i=1,3
         rx(i) = 0.0d0
         ry(i) = 0.0d0
         rz(i) = 0.0d0
         t(i)  = 0.0d0
      end do

      rx(1) = 1.0d0
      ry(2) = 1.0d0
      rz(3) = 1.0d0

      return
      end

      subroutine strot(rx,ry,rz)
      implicit double precision (a-h,o-z)
      dimension rx(3),ry(3),rz(3)

      do i=1,3
         rx(i) = 0.0d0
         ry(i) = 0.0d0
         rz(i) = 0.0d0
      end do

      rx(1) = 1.0d0
      ry(3) = 1.0d0
      rz(2) = -1.0d0

      return
      end

      subroutine setord(iatom,t,coo)
      implicit double precision (a-h,o-z)
      integer dolabs,fancy,persp,shade,atcol,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      dimension t(3),coo(3,*)

      if (irtcel.eq.2) then
         call alnsorg(iatom,coo)
      else
         do i=1,3
            t(i) = coo(i,iatom)
         end do
      endif

c      call doscal

      call setxyv

      return
      end

      subroutine setxyd(xv,yv)
      implicit double precision (a-h,o-z)

      xv = 0.0d0
      yv = 0.0d0

      return
      end

      subroutine prtrot(rx,ry,rz,t)
      implicit double precision (a-h,o-z)
      dimension rx(3),ry(3),rz(3),t(3)

      print*,(rx(i),i=1,3)
      print*,(ry(i),i=1,3)
      print*,(rz(i),i=1,3)

      print*,(t(i),i=1,3)

      return
      end

      subroutine rarbxi
      implicit double precision (a-h,o-z)
      dimension v(3)
      common /rligx/ ry(3,3),rz(3,3),rzi(3,3),ryi(3,3)
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt

      todeg = 45.0d0 / datan(1.0d0)
      do i=1,3
         v(i) = vecs(1,i)
      end do

c Rotation about an arbitrary axis v (going through the origin)

      pi = 4.d0*datan(1.d0)

c calculate rho, theta, phi (Spherical coordinates)

      rho = vlen(v)

      if (v(1).gt.0.0d0) then
         theta = datan(v(2)/v(1))
      elseif (v(1).lt.0.0d0) then
         theta = pi + datan(v(2)/v(1))
      elseif (v(1).eq.0.0d0) then
         if (v(2).ge.0.0d0) then
            theta = pi / 2.0d0
         else
            theta = 3.0d0 * pi / 2.0d0
         endif
      endif
      
      phi = dacos(v(3)/rho)
      
      cs = dcos(theta)
      si = dsin(theta)

      rz(1,1) = cs
      rz(1,2) = si
      rz(1,3) = 0.0d0

      rz(2,1) = -si
      rz(2,2) = cs
      rz(2,3) = 0.0d0

      rz(3,1) = 0.0d0
      rz(3,2) = 0.0d0
      rz(3,3) = 1.0d0

      rzi(1,1) = cs
      rzi(1,2) = -si
      rzi(1,3) = 0.0d0

      rzi(2,1) = si
      rzi(2,2) = cs
      rzi(2,3) = 0.0d0

      rzi(3,1) = 0.0d0
      rzi(3,2) = 0.0d0
      rzi(3,3) = 1.0d0

      cs1 = dcos(phi)
      si1 = dsin(phi)

      ry(1,1) = cs1
      ry(1,2) = 0.0d0
      ry(1,3) = -si1

      ry(2,1) = 0.0d0
      ry(2,2) = 1.0d0
      ry(2,3) = 0.0d0

      ry(3,1) = si1
      ry(3,2) = 0.0d0
      ry(3,3) = cs1

      ryi(1,1) = cs1
      ryi(1,2) = 0.0d0
      ryi(1,3) = si1

      ryi(2,1) = 0.0d0
      ryi(2,2) = 1.0d0
      ryi(2,3) = 0.0d0

      ryi(3,1) = -si1
      ryi(3,2) = 0.0d0
      ryi(3,3) = cs1

c      print*,'matrix ry'
c      call prtmat(ry)
c      print*,'matrix rz'
c      call prtmat(rz)
c      print*,'matrix ryi'
c      call prtmat(ryi)
c      print*,'matrix rzi'
c      call prtmat(rzi)
      return

      end

      subroutine rarbx(alfat)
      implicit double precision (a-h,o-z)
      dimension r(3,3),rv(3,3),rout(3,3),rtmp(3,3),rcpy(3,3)
      common /rligx/ ry(3,3),rz(3,3),rzi(3,3),ryi(3,3)
      common /cllmat/rxl(3),ryl(3),rzl(3),t(3),tz(3),tzorg(3),itz

      do i=1,3
         rcpy(i,1) = rxl(i)
         rcpy(i,2) = ryl(i)
         rcpy(i,3) = rzl(i)
      end do

c Rotation about an arbitrary axis v (going through the origin)
c alpha is the angle over which is rotated along axis v

      alfa = -alfat

      ca = dcos(alfa)
      sa = dsin(alfa)

      rv(1,1) = ca
      rv(1,2) = sa
      rv(1,3) = 0.0d0

      rv(2,1) = -sa
      rv(2,2) = ca
      rv(2,3) = 0.0d0

      rv(3,1) = 0.0d0
      rv(3,2) = 0.0d0
      rv(3,3) = 1.0d0

      call matmult(rzi,ryi,rtmp)
      call matmult(rtmp,rv,rout)
      call matmult(rout,ry,rtmp)
      call matmult(rtmp,rz,r)
      call matmult(rcpy,r,rtmp)

      do i=1,3
         rxl(i) = rtmp(i,1)
         ryl(i) = rtmp(i,2)
         rzl(i) = rtmp(i,3)
      end do

      return
      end

      subroutine rarbyi
      implicit double precision (a-h,o-z)
      dimension v(3)
      common /rligy/ ry(3,3),rz(3,3),rzi(3,3),ryi(3,3)
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt

      do i=1,3
         v(i) = vecs(2,i)
      end do

c Rotation about an arbitrary axis v (going through the origin)

      pi = 4.d0*datan(1.d0)

c calculate rho, theta, phi (Spherical coordinates)

      rho = vlen(v)

      if (v(1).gt.0.0d0) then
         theta = datan(v(2)/v(1))
      elseif (v(1).lt.0.0d0) then
         theta = pi + datan(v(2)/v(1))
      elseif (v(1).eq.0.0d0) then
         if (v(2).ge.0.0d0) then
            theta = pi / 2.0d0
         else
            theta = 3.0d0 * pi / 2.0d0
         endif
      endif
      
      phi = dacos(v(3)/rho)
      
      cs = dcos(theta)
      si = dsin(theta)

      rz(1,1) = cs
      rz(1,2) = si
      rz(1,3) = 0.0d0

      rz(2,1) = -si
      rz(2,2) = cs
      rz(2,3) = 0.0d0

      rz(3,1) = 0.0d0
      rz(3,2) = 0.0d0
      rz(3,3) = 1.0d0

      rzi(1,1) = cs
      rzi(1,2) = -si
      rzi(1,3) = 0.0d0

      rzi(2,1) = si
      rzi(2,2) = cs
      rzi(2,3) = 0.0d0

      rzi(3,1) = 0.0d0
      rzi(3,2) = 0.0d0
      rzi(3,3) = 1.0d0

      cs1 = dcos(phi)
      si1 = dsin(phi)

      ry(1,1) = cs1
      ry(1,2) = 0.0d0
      ry(1,3) = -si1

      ry(2,1) = 0.0d0
      ry(2,2) = 1.0d0
      ry(2,3) = 0.0d0

      ry(3,1) = si1
      ry(3,2) = 0.0d0
      ry(3,3) = cs1

      ryi(1,1) = cs1
      ryi(1,2) = 0.0d0
      ryi(1,3) = si1

      ryi(2,1) = 0.0d0
      ryi(2,2) = 1.0d0
      ryi(2,3) = 0.0d0

      ryi(3,1) = -si1
      ryi(3,2) = 0.0d0
      ryi(3,3) = cs1

c      print*,'y matrix ry'
c      call prtmat(ry)
c      print*,'y matrix rz'
c      call prtmat(rz)
c      print*,'y matrix ryi'
c      call prtmat(ryi)
c      print*,'y matrix rzi'
c      call prtmat(rzi)

      return
      end

      subroutine rarby(alfa)
      implicit double precision (a-h,o-z)
      dimension r(3,3),rv(3,3),rout(3,3),rtmp(3,3),rcpy(3,3)
      common /rligy/ ry(3,3),rz(3,3),rzi(3,3),ryi(3,3)
      common /cllmat/rxl(3),ryl(3),rzl(3),t(3),tz(3),tzorg(3),itz

      do i=1,3
         rcpy(i,1) = rxl(i)
         rcpy(i,2) = ryl(i)
         rcpy(i,3) = rzl(i)
      end do


c Rotation about an arbitrary axis v (going through the origin)
c alpha is the angle over which is rotated along axis v

      ca = dcos(alfa)
      sa = dsin(alfa)

      rv(1,1) = ca
      rv(1,2) = sa
      rv(1,3) = 0.0d0

      rv(2,1) = -sa
      rv(2,2) = ca
      rv(2,3) = 0.0d0

      rv(3,1) = 0.0d0
      rv(3,2) = 0.0d0
      rv(3,3) = 1.0d0

      call matmult(rzi,ryi,rtmp)
      call matmult(rtmp,rv,rout)
      call matmult(rout,ry,rtmp)
      call matmult(rtmp,rz,r)
      call matmult(rcpy,r,rtmp)

      do i=1,3
         rxl(i) = rtmp(i,1)
         ryl(i) = rtmp(i,2)
         rzl(i) = rtmp(i,3)
      end do

      return
      end

      subroutine matmult(a,b,c)
      implicit double precision (a-h,o-z)
      dimension a(3,3),b(3,3),c(3,3)

      do i=1,3
         do j=1,3
            c(i,j)=0.0d0
            do k=1,3
               c(i,j) = c(i,j)+a(i,k)*b(k,j)
            end do
         end do
      end do

      return
      end

      subroutine rarbzi
      implicit double precision (a-h,o-z)
      dimension v(3)
      common /rligz/ ry(3,3),rz(3,3),rzi(3,3),ryi(3,3)
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt

      todeg = 45.0d0 / datan(1.0d0)
      do i=1,3
         v(i) = vecs(3,i)
      end do

c Rotation about an arbitrary axis v (going through the origin)

      pi = 4.d0*datan(1.d0)

c calculate rho, theta, phi (Spherical coordinates)

      rho = vlen(v)

      if (v(1).gt.0.0d0) then
         theta = datan(v(2)/v(1))
      elseif (v(1).lt.0.0d0) then
         theta = pi + datan(v(2)/v(1))
      elseif (v(1).eq.0.0d0) then
         if (v(2).ge.0.0d0) then
            theta = pi / 2.0d0
         else
            theta = 3.0d0 * pi / 2.0d0
         endif
      endif
      
      phi = dacos(v(3)/rho)
      
      cs = dcos(theta)
      si = dsin(theta)

      rz(1,1) = cs
      rz(1,2) = si
      rz(1,3) = 0.0d0

      rz(2,1) = -si
      rz(2,2) = cs
      rz(2,3) = 0.0d0

      rz(3,1) = 0.0d0
      rz(3,2) = 0.0d0
      rz(3,3) = 1.0d0

      rzi(1,1) = cs
      rzi(1,2) = -si
      rzi(1,3) = 0.0d0

      rzi(2,1) = si
      rzi(2,2) = cs
      rzi(2,3) = 0.0d0

      rzi(3,1) = 0.0d0
      rzi(3,2) = 0.0d0
      rzi(3,3) = 1.0d0

      cs1 = dcos(phi)
      si1 = dsin(phi)

      ry(1,1) = cs1
      ry(1,2) = 0.0d0
      ry(1,3) = -si1

      ry(2,1) = 0.0d0
      ry(2,2) = 1.0d0
      ry(2,3) = 0.0d0

      ry(3,1) = si1
      ry(3,2) = 0.0d0
      ry(3,3) = cs1

      ryi(1,1) = cs1
      ryi(1,2) = 0.0d0
      ryi(1,3) = si1

      ryi(2,1) = 0.0d0
      ryi(2,2) = 1.0d0
      ryi(2,3) = 0.0d0

      ryi(3,1) = -si1
      ryi(3,2) = 0.0d0
      ryi(3,3) = cs1

c      print*,'matrix ry'
c      call prtmat(ry)
c      print*,'matrix rz'
c      call prtmat(rz)
c      print*,'matrix ryi'
c      call prtmat(ryi)
c      print*,'matrix rzi'
c      call prtmat(rzi)
      return

      end

      subroutine rarbz(alfa)
      implicit double precision (a-h,o-z)
      dimension r(3,3),rv(3,3),rout(3,3),rtmp(3,3),rcpy(3,3)
      common /rligz/ ry(3,3),rz(3,3),rzi(3,3),ryi(3,3)
      common /cllmat/rxl(3),ryl(3),rzl(3),t(3),tz(3),tzorg(3),itz

      do i=1,3
         rcpy(i,1) = rxl(i)
         rcpy(i,2) = ryl(i)
         rcpy(i,3) = rzl(i)
      end do

c Rotation about an arbitrary axis v (going through the origin)
c alpha is the angle over which is rotated along axis v


      ca = dcos(alfa)
      sa = dsin(alfa)

      rv(1,1) = ca
      rv(1,2) = sa
      rv(1,3) = 0.0d0

      rv(2,1) = -sa
      rv(2,2) = ca
      rv(2,3) = 0.0d0

      rv(3,1) = 0.0d0
      rv(3,2) = 0.0d0
      rv(3,3) = 1.0d0

      call matmult(rzi,ryi,rtmp)
      call matmult(rtmp,rv,rout)
      call matmult(rout,ry,rtmp)
      call matmult(rtmp,rz,r)
      call matmult(rcpy,r,rtmp)

      do i=1,3
         rxl(i) = rtmp(i,1)
         ryl(i) = rtmp(i,2)
         rzl(i) = rtmp(i,3)
      end do

      return
      end

      subroutine prtmat(a)
      implicit double precision (a-h,o-z)
      dimension a(3,3)

      do i=1,3
         write(6,'(3f12.4)') (a(i,j),j=1,3)
      end do

      return
      end

      subroutine rotmt(rot,t1,t2,t3)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265d0)
      parameter (rad=pi/180.0d0)
      dimension rot(3,3)

      tp = t1
      tm = t3
      tt1 = (tp + tm)/2.0d0
      tt3 = (tp - tm)/2.0d0
 
      s1 = dsin(rad*tt1)
      s2 = dsin(rad*t2)
      s3 = dsin(rad*tt3)
      c1 = dcos(rad*tt1)
      c2 = dcos(rad*t2)
      c3 = dcos(rad*tt3)

      rot(1,1) = -s1*c2*s3 + c1*c3
      rot(1,2) =  c1*c2*s3 + s1*c3                   
      rot(1,3) =  s2*s3
      rot(2,1) = -s1*c2*c3 - c1*s3
      rot(2,2) =  c1*c2*c3 - s1*s3
      rot(2,3) =  s2*c3
      rot(3,1) =  s1*s2
      rot(3,2) = -c1*s2
      rot(3,3) =  c2

      return
      end

      subroutine supimp(coo1,coo2,iatoms1,iatoms2,
     &                  iamin1,iamin2,icalf1,icalf2,ncalf1,ncalf2,
     &                  isal1,isal2)
      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      logical ok
      dimension coo1(3,*),coo2(3,*),iamin1(*),iamin2(*),
     &          icalf1(6,*),icalf2(6,*),isal1(*),isal2(*)
      dimension ca1(3,numcal),ca2(3,numcal),
     &          x(3,numcal), y(3,numcal),
     &          iaatp1(numcal),iaatp2u(numcal)
      dimension ang(3),trans(3),rot(3,3),trans2(3)

      call fident(coo1,coo2,iamin1,iamin2,icalf1,icalf2,
     &            ncalf1,ncalf2,isal1,isal2,x,y,neq,ok)

      if (.not.ok) then

         print*,'identity too low switching to slower algorithm'

         call getca(coo1,coo2,iamin1,iamin2,icalf1,icalf2,
     &              ncalf1,ncalf2,ca1,ca2,iaatp1,iaatp2u)

c angscan is the slow part

         call angscan(ca1,ncalf1,ca2,ncalf2,ang,trans,iaatp1,iaatp2u)

         call getca(coo1,coo2,iamin1,iamin2,icalf1,icalf2,
     &              ncalf1,ncalf2,ca1,ca2,iaatp1,iaatp2u)

         call getequ(ca1,ncalf1,ca2,ncalf2,ang,trans,x,y,neq)

      endif

      call kabsch(x,y,neq,rot,trans2,ier)

      call rotcoo(coo2,iatoms2,trans2,rot)

      return
      end

      subroutine rotcoo(coo,iatoms,trans,rot)
      implicit double precision (a-h,o-z)
      dimension coo(3,*),rot(3,3),trans(3),x(3)

      do i=1,iatoms

         do j=1,3
            x(j) = rot(j,1)*coo(1,i) + 
     &             rot(j,2)*coo(2,i) + 
     &             rot(j,3)*coo(3,i) + trans(j)/0.52917706d0
         end do

         do j=1,3
            coo(j,i) = x(j)
         end do

      end do

      return
      end

      subroutine getca(coo1,coo2,iamin1,iamin2,icalf1,icalf2,
     &                 ncalf1,ncalf2,ca1,ca2,iaatp1,iaatp2u)
      implicit double precision (a-h,o-z)
      integer aatyp
      dimension coo1(3,*),coo2(3,*),iamin1(*),iamin2(*),
     &          icalf1(6,*),icalf2(6,*)
      dimension ca1(3,*),ca2(3,*),iaatp1(*),iaatp2u(*)
      dimension aatyp(42)

c aatyp:
c
c   1 hydrophilic 
c   2 hydrophobic 
c   3 large aromatic
c   4 nucleotide
      
      data aatyp/2,2,1,2,1,
     &           2,2,2,1,1,
     &           2,1,1,1,2,
     &           1,3,3,3,3,
     &           1,1,3,4,4,
     &           4,4,4,4,4,
     &           4,4,4,4,4,
     &           4,4,4,4,4,
     &           4,4/

      do i=1,ncalf1
         do j=1,3
            ca1(j,i) = coo1(j,icalf1(1,i))*0.52917706d0
         end do
         iaatp1(i) = aatyp(iamin1(i))
      end do

      do i=1,ncalf2
         do j=1,3
            ca2(j,i) = coo2(j,icalf2(1,i))*0.52917706d0
         end do
         iaatp2u(i) = aatyp(iamin2(i))
      end do

      return
      end

      subroutine getequ(ca1,ncalf1,ca2,ncalf2,ang,trans,x,y,neq)
      implicit double precision (a-h,o-z)
      dimension ca1(3,*),ca2(3,*),x(3,*),y(3,*),
     &          rot(3,3),ang(3),trans(3),t(3)

      cutoff = 1.9d0
      neq  = 0
      dsum = 0.0d0      

      call rotmt(rot,ang(1),ang(2),ang(3))

c transform coordinates according to best solution

      do j=1,ncalf2

        do i=1,3
           t(i) = rot(i,1)*ca2(1,j) + rot(i,2)*ca2(2,j) + 
     &            rot(i,3)*ca2(3,j) + trans(i)
        end do

c find closest equivalent in reference. 

        dmin = 10000

        do k=1,ncalf1

          d = (ca1(1,k)-t(1))**2 + (ca1(2,k)-t(2))**2 + 
     &        (ca1(3,k)-t(3))**2

          if (d.lt.dmin.and.d.lt.cutoff**2) then
            dmin  = d
            ksave = k
          endif

        end do

        if (dmin.lt.cutoff**2) then 

          dsum = dsum + dmin
          neq  = neq  + 1

          do k=1,3
            x(k,neq) = ca2(k,j)
            y(k,neq) = ca1(k,ksave)
          end do

        endif

      end do

      return
      end

      subroutine stosol(over,corr,np,nm,qaqa,
     &        n2cc,n2n2,phip,phi2,phim,
     &        trans,transb,qapr,ccpr,n2pr,
     &        n2n2m,corrp,n2prp,
     &        ccprp,ang)

      implicit double precision (a-h,o-z)
      dimension ang(3),trans(3),transb(3)

      if (np.gt.n2cc) then
        n2cc    = ncntp
      endif

      if (np+nm.gt.n2n2+n2n2m) then
        n2n2    = np
        n2n2m   = nm
      endif

      if (np*corr.gt.n2pr*ccpr) then
        qapr    = over
        ccpr    = corr
        n2pr    = np
      endif

      if (np*corrp.gt.n2prp*ccprp) then
        ccprp     = corrp
        n2prp     = np
        ang(1)  = phip
        ang(2)  = phi2
        ang(3)  = phim

        transb(1) = trans(1)
        transb(2) = trans(2)
        transb(3) = trans(3)
      endif

      if (over.gt.qaqa) then
        qaqa   = over
      endif

      return
      end

      subroutine calcor(ca1,cat,ncalf1,ncalf2,trans,n,np,nm,
     &                  corr,corrp)
      implicit double precision (a-h,o-z)
      dimension ca1(3,*),cat(3,*),trans(3)

      j1tot = 0
      j2tot = 0
      j1j2t = 0
      j1j1t = 0
      j2j2t = 0
      i1tot = 0
      i2tot = 0
      i1i2t = 0
      i1i1t = 0
      i2i2t = 0
      n     = 0
      np    = 0
      nm    = 0
      ilast = -1
      naddp = 1
      naddm = 1

      do k=1,ncalf1

        dstmin = 10000.0d0

c look for nearest equivalent atom

        do i=1,ncalf2

          dist = 0.0d0

          do j=1,3
             dist = dist + (cat(j,i) + trans(j) - ca1(j,k))**2 
          end do

          if (dist.lt.dstmin) then
            dstmin = dist
            near   = i
          endif

        end do

        if (dstmin.lt.(2.5d0)**2) then

          n     = n + 1
          i1    = k
          i2    = near
          i1tot = i1tot + i1        
          i2tot = i2tot + i2
          i1i2t = i1i2t + i1*i2
          i1i1t = i1i1t + i1*i1
          i2i2t = i2i2t + i2*i2

          if (near-ilast.eq.1) then

            naddp = naddp + 1

            if (naddp.eq.2.or.naddp.eq.1) then 

              np = np + naddp 

              do j=0,naddp-1
                i1    = k     - j
                i2    = near  - j
                j1tot = j1tot + i1
                j2tot = j2tot + i2
                j1j2t = j1j2t + i1*i2
                j1j1t = j1j1t + i1*i1
                j2j2t = j2j2t + i2*i2
              end do

              naddp = 0

            endif

            naddm = 1

          elseif (ilast-near.eq.1) then

            naddm = naddm + 1

            if (naddm.eq.2.or.naddm.eq.1) then
              nm = nm + naddm  
              naddm = 0
            endif

            naddp = 1

          else

            naddp = 1
            naddm = 1

          endif

          ilast = near

        else

          ilast = -1
          naddp =  1
          naddm =  1

        endif

      end do

      if (n.gt.0.and.n*i1i1t.gt.i1tot**2.and.n*i2i2t.gt.i2tot**2) then        

        corr =  (dble(i1i2t) - dble(i1tot*i2tot)/n) /
     &     dsqrt(dble(i1i1t) - dble(i1tot*i1tot)/n) /
     &     dsqrt(dble(i2i2t) - dble(i2tot*i2tot)/n)

      else

        corr = 0.0d0

      endif

      if (np.gt.0.and.np*j1j1t.gt.j1tot**2.and.
     &    np*j2j2t.gt.j2tot**2) then        

        corrp =  (dble(j1j2t) - dble(j1tot*j2tot)/np) /
     &      dsqrt(dble(j1j1t) - dble(j1tot*j1tot)/np) /
     &      dsqrt(dble(j2j2t) - dble(j2tot*j2tot)/np)

      else

        corrp = 0.0d0

      endif

      return
      end

      subroutine sort(x,n,iaatp)
      implicit double precision (a-h,o-z)
      parameter (np=850)
      parameter (numcal=50000)
      dimension x(3,*),tmp(np),indx(numcal),iaatp(*)

      if (n.gt.np) return 

      call v3cpy2(tmp,x,n,1)

      call rqsrt(n,tmp,indx)

      do i=1,n
         x(1,i) = tmp(indx(i))
      end do

      do i=1,n
         tmp(i) = x(2,indx(i))
      end do

      call v3cpy1(x,tmp,n,2)

      do i=1,n
         tmp(i) = x(3,indx(i))
      end do

      call v3cpy1(x,tmp,n,3)

      do i=1,n
         tmp(i) = iaatp(indx(i))
      end do

      do i=1,n
         iaatp(i) = nint(tmp(i))
      end do

      return
      end

      subroutine dotran(ca1,ncalf1,ca2,ncalf2,trans,best,maxd,
     &                  iaatp1,iaatp2)
      implicit double precision (a-h,o-z)
      parameter (np=850)
      dimension grid(-maxd:maxd,-maxd:maxd,-maxd:maxd)
      dimension ca1(3,*),ca2(3,*),trans(3),pam(3,3),
     &          ica2(3,np),iaatp1(*),iaatp2(*)
      data pam/1.0d0,0.3d0,0.6d0,
     &         0.3d0,1.0d0,0.6d0,
     &         0.6d0,0.6d0,1.0d0/

      if (ncalf2.gt.np) stop 'dotran: np underdim.'

      xcut = dble(min(ncalf1,ncalf2))/100.0d0
      
      maxd1 = (maxd-1)
      maxd2 = maxd1*maxd1
      best  = 0.0d0

      do i1=-maxd,maxd
         do i2=-maxd,maxd
            do i3=-maxd,maxd
               grid(i3,i2,i1) = 0.0d0
            end do
         end do
      end do

      do i=1,ncalf2
         do j=1,3
            ica2(j,i) = nint(ca2(j,i))
         end do
      end do

      jstart = 0

      do i=1,ncalf1

        ica1x = nint(ca1(1,i))
        ica1y = nint(ca1(2,i))
        ica1z = nint(ca1(3,i))

        do j=jstart+1,ncalf2
          i1 = ica1x - ica2(1,j)

          if (i1.gt.maxd1) then
            jstart = j
          else
             if (i1.lt.-maxd1) goto 100

             i2 = ica1y - ica2(2,j)

             if (abs(i2).le.maxd1) then
                i3 = ica1z - ica2(3,j)

                if (abs(i3).le.maxd1) then
                   grid(i3,i2,i1) = grid(i3,i2,i1) + 
     &                              pam(iaatp1(i),iaatp2(j))
                endif

             endif

          endif

        end do

100     continue
      end do

c     search within a sphere around the origin

      do i1=-maxd+1,maxd-1

         i = int(dsqrt(dble(maxd2-i1*i1)))

         do i2=-i,i

            j = int(dsqrt(dble(maxd2-i1*i1-i2*i2)))

            do i3=-j,j

              if (grid(i3,i2,i1).gt.xcut) then

                over = 0.0d0

                do k=-1,1

                   do l=-1,1

                      do m=-1,1
                         over = over + grid(m+i3,l+i2,k+i1)
                      end do

                   end do

                end do

                if (over.gt.best) then

                  best = over

                  call vsetr(trans,3,0.0d0)

                  do k=-1,1
                     do l=-1,1
                        do m=-1,1
                          trans(1) = trans(1) + 
     &                               grid(m+i3,l+i2,k+i1)*(k+i1)
                          trans(2) = trans(2) + 
     &                               grid(m+i3,l+i2,k+i1)*(l+i2)
                          trans(3) = trans(3) + 
     &                               grid(m+i3,l+i2,k+i1)*(m+i3)
                        end do
                     end do
                  end do

                  call vscal(trans,3,1.0d0/over)

                endif

              endif

            end do
         end do
      end do

      return
      end

      subroutine angscan(ca1,ncalf1,ca2,ncalf2,ang,transb,
     &                   iaatp1,iaatp2u)
      implicit double precision (a-h,o-z)
      parameter (rad=3.14159265d0/180.0d0)
      parameter (numcal=50000,maxd=20)
      dimension ca1(3,*),ca2(3,*),rot(3,3),trans(3),
     &  cen1(3),cen2(3),transb(3),ang(3)
      dimension iaatp1(*),iaatp2u(*)
      dimension cat(3,numcal),ca1u(3,numcal),catu(3,numcal),
     &          iaatp2(numcal)

      delt2 = 10.0d0

      rmom1 = 0.0d0
      rmom2 = 0.0d0

      call cntvc2(cen1,ca1,ncalf1)
      call cntvc2(cen2,ca2,ncalf2)

      do i=1,ncalf1
         do j=1,3
            ca1(j,i) = ca1(j,i) - cen1(j)
            rmom1 = rmom1 + ca1(j,i)**2
            ca1u(j,i) = ca1(j,i)
         end do
      end do

      do i=1,ncalf2
         do j=1,3
            ca2(j,i) = ca2(j,i) - cen2(j)
            rmom2 = rmom2 + ca2(j,i)**2
         end do
      end do

      if (ncalf1.gt.0) rmom1 = dsqrt(rmom1/ncalf1)
      if (ncalf2.gt.0) rmom2 = dsqrt(rmom2/ncalf2)

      maxd1 = nint(9.5d0 + dabs(rmom1-rmom2))

      if (maxd1.gt.maxd) maxd1 = maxd

      ccprp    = 0.0d0
      ccpr     = 0.0d0
      qaqa     = 0.0d0
      qapr     = 0.0d0
      totover2 = 0.0d0
      totover  = 0.0d0

      n2cc   = 0
      n2n2   = 0
      n2pr   = 0
      n2n2m  = 0
      n2prp  = 0
      nang   = 0

      call sort(ca1,ncalf1,iaatp1)

      phi2 = -delt2

      do iphi2=0,int(180.01d0/delt2)

        phi2 = phi2 + delt2

        if (phi2.lt.179.9d0) then
          deltp = delt2 / dcos(phi2/2.0d0*rad)
        else
          deltp = 720.01d0
        endif

        phip = -deltp

        do iphip=0,int(719.99d0/deltp)

          phip = phip + deltp

          if (iphi2.ne.0) then

            deltm = delt2/dsin(phi2/2.0d0*rad)

          else

            if (phip.gt.359.99d0) goto 100
            deltm = 360.02d0

          endif

          phim = -deltm

          do iphim=0,int(360.01d0/deltm)

            phim = phim + deltm

            if (phi2.gt.179.99d0.and.phim.gt.359.99d0) goto 100

            call rotmt(rot,phip,phi2,phim)

            do i=1,ncalf2

               iaatp2(i) = iaatp2u(i)

               do j=1,3
                  cat(j,i) = rot(j,1)*ca2(1,i) + rot(j,2)*ca2(2,i) + 
     &                       rot(j,3)*ca2(3,i)       
                  catu(j,i) = cat(j,i)
               end do

            end do

            call sort(cat,ncalf2,iaatp2)

            call dotran(ca1,ncalf1,cat,ncalf2,trans,over,maxd1,iaatp1,
     &                  iaatp2)

            nang     = nang + 1
            totover  = totover + over
            totover2 = totover2 + over*over
            corr     = 0.0d0
            corrp    = 0.0d0
            np       = 0
            nm       = 0

            overc = totover/nang + 2.0d0*
     &              dsqrt(dabs(totover2/nang - (totover/nang)**2))

            if (over.gt.qaqa.or.over.gt.overc)
     &          call calcor(ca1u,catu,ncalf1,ncalf2,trans,n,np,nm,
     &                      corr,corrp)        

c store solution

            call stosol(over,corr,np,nm,qaqa,n2cc,n2n2,
     &        phip,phi2,phim,
     &        trans,transb,qapr,ccpr,n2pr,n2n2m,corrp,n2prp,
     &        ccprp,ang)

c fine scan 

            overc = totover/nang + 3.0d0*
     &              dsqrt(dabs(totover2/nang - (totover/nang)**2))

            if (corrp*np.gt.0.9d0*(n2prp*ccprp).or.
     &          np+nm.gt.0.9d0*(n2n2+n2n2m).or.over.gt.overc) then

             do i1=-1,1
                do i2=-1,1
                   do i3=-1,1

                      if (.not.(i1.eq.0.and.i2.eq.0.and.i3.eq.0)) then

                          call rotmt(rot,phip + i1*deltp/2.0d0,
     &                                   phi2 + i2*delt2/2.0d0,
     &                                   phim + i3*deltm/2.0d0)
                          do i=1,ncalf2

                             iaatp2(i) = iaatp2u(i)

                             do j=1,3
                                cat(j,i) = rot(j,1)*ca2(1,i) + 
     &                                     rot(j,2)*ca2(2,i) + 
     &                                     rot(j,3)*ca2(3,i)      
                                catu(j,i) = cat(j,i)
                             end do

                          end do

                          call sort(cat,ncalf2,iaatp2)

                          call dotran(ca1,ncalf1,cat,ncalf2,trans,over,
     &                                maxd1,iaatp1,iaatp2)

                          call calcor(ca1u,catu,ncalf1,ncalf2,trans,
     &                                n,np,nm,corr,corrp)        

                          call stosol(over,corr,np,nm,qaqa,
     &                                n2cc,n2n2,
     &                                phip+i1*deltp/2.0d0,
     &                                phi2+i2*delt2/2.0d0,
     &                                phim+i3*deltm/2.0d0,
     &                                trans,transb,qapr,ccpr,n2pr,n2n2m,
     &                                corrp,n2prp,ccprp,ang)

                      endif
                   end do
                end do
             end do

            endif

          end do
        end do
100     continue
      end do

      totover = totover / nang
      totover2 = dsqrt(dabs(totover2/nang - totover**2))

      call rotmt(rot,ang(1),ang(2),ang(3)) 

      do j=1,3
         cat(j,1) = rot(j,1)*cen2(1) + rot(j,2)*cen2(2) + 
     &              rot(j,3)*cen2(3)
         transb(j) = transb(j) + cen1(j) - cat(j,1)
      end do

      return
      end         

      subroutine kabsch(x,y,n,u,t,ier)
      implicit double precision (a-h,o-z)
      double precision mu(3)
      dimension x(3,*),y(3,*),u(3,3),t(3),ip(9),r(3,3),
     &          a(3,3),b(3,3),rr(6),cen1(3),cen2(3)
      data ip /1,2,4,  2,3,5,  4,5,6/

      ier = -1

      if (n.lt.2) then

         do i=1,3
            do j=1,3
               if (i.eq.j) then
                  u(i,j) = 1.0d0
               else
                  u(i,j) = 0.0d0
               endif
            end do
         end do

         call vsetr(t,3,0.0d0)

         return
      endif

      call cntvc2(cen1,x,n)
      call cntvc2(cen2,y,n)

      do i=1,3

         do j=1,3
            r(i,j) = 0.0d0
         end do

      end do

      xyavl = 0.0d0

c R is the covariance matrix

      do m=1,n
         do i=1,3

            xyavl = xyavl + 
     &              ((x(i,m) - cen1(i))**2 + (y(i,m) - cen2(i))**2)

            do j=1,3
               r(i,j) = r(i,j) + (y(i,m) - cen2(i)) * (x(j,m) - cen1(j))
            end do

         end do
      end do

c     determinant of r(i,j)

      det = r(1,1)*(r(2,2)*r(3,3)-r(2,3)*r(3,2))
     &     -r(1,2)*(r(2,1)*r(3,3)-r(2,3)*r(3,1))
     &     +r(1,3)*(r(2,1)*r(3,2)-r(2,2)*r(3,1))

      dett = det

c     RR = transposed(R)*R

      m = 0
      do j=1,3
         do i=1,j
            m = m + 1
            rr(m) = r(1,i)*r(1,j) + r(2,i)*r(2,j) + r(3,i)*r(3,j)
         end do
      end do

      tracerr = (rr(1) + rr(3) + rr(6)) / 3.0d0
      cof = (rr(3)*rr(6) - rr(5)*rr(5) + rr(1)*rr(6) - rr(4)*rr(4) + 
     &       rr(1)*rr(3) - rr(2)*rr(2) ) / 3.0d0
      det = det*det

c     reduce cubic to standard form y**3-3hy+2g=0 by putting x=y-trace

      d = tracerr*tracerr
      h = d - cof
      g = tracerr*(cof*1.5d0 - d) - det*0.5d0

c     eigenvalues RR: mu(3), eigenvectors of RR: a(3,3) (A)

      if (h.le.d*1.0d-9) then

c     three roots the same

         mu(1) = tracerr
         mu(2) = tracerr
         mu(3) = tracerr

         do i=1,3
            do j=1,3
   
               if (i.eq.j) then
                  a(i,j) = 1.0d0
               else
                  a(i,j) = 0.0d0
               endif

            end do
         end do

      else

         sqrth = dsqrt(h)
         d = -g / (h*sqrth)
   
         if (d.gt. 0.9999999d0) then

c        two roots the same, two different 

            mu(1) = tracerr + sqrth + sqrth
            mu(2) = tracerr - sqrth
            mu(3) = mu(2)

            call root2(a,rr,mu,ip,1,3,1,2)

         else

            if (d.lt.-0.9999999d0) then

c           two roots the same, two different 

               mu(1) = tracerr + sqrth
               mu(2) = mu(1)
               mu(3) = tracerr - sqrth - sqrth

               call root2(a,rr,mu,ip,3,1,2,3)

            else

c           three different roots

               d   = dacos(d)/3.0d0
               cphi = sqrth*dcos(d)
               sphi = sqrth*dsin(d)*dsqrt(3.0d0)

               mu(1)  = tracerr + cphi + cphi
               mu(2)  = tracerr - cphi + sphi
               mu(3)  = tracerr - cphi - sphi

               call root3(a,rr,mu,ip)

            endif

         endif

      endif

c B = R*A

      do l=1,2

         d = 0.0d0

         do i=1,3
            b(i,l) = r(i,1)*a(1,l) + r(i,2)*a(2,l) + r(i,3)*a(3,l)
            d = d + b(i,l)**2
         end do

         d = dsqrt(d)

         do i=1,3
            b(i,l) = b(i,l) / d
         end do

      end do

      call axa(b,3,1,2)

c     rotation matrix U = B*A (U = At*R*A )

      do i=1,3
         do j=1,3
            u(i,j) = b(i,1)*a(j,1) + b(i,2)*a(j,2) + b(i,3)*a(j,3)
         end do
      end do

c     translation vector: center y - rotated center x

      do i=1,3
         t(i) = cen2(i) - 
     &         (u(i,1)*cen1(1) + u(i,2)*cen1(2) + u(i,3)*cen1(3))
      end do

c     rms error

      if (mu(3).lt.0.0d0) mu(3) = 0.0d0
      if (mu(2).lt.0.0d0) mu(2) = 0.0d0
      if (mu(1).lt.0.0d0) mu(1) = 0.0d0

      d = dsqrt(mu(3))
      if (dett.lt.0.0) d = -d
      d = d + dsqrt(mu(2)) + dsqrt(mu(1))
      rms = dabs(xyavl-2.0d0*d)
      print*,'rms unnorm ',rms
      rms = dsqrt(rms/n)

      print*,'rms ',rms,' (angstroms)'

      ier = 0

      return
      end

c     eigenvectors three distinct roots

      subroutine root3(a,rr,e,ip)
      implicit double precision (a-h,o-z)
      dimension a(3,3),ss(6),rr(6),e(3),ip(9)

      do l=1,2

         ss(1) = (e(l)-rr(3))*(e(l)-rr(6)) - rr(5)*rr(5)
         ss(2) = (e(l)-rr(6))*      rr(2)  + rr(4)*rr(5)
         ss(3) = (e(l)-rr(1))*(e(l)-rr(6)) - rr(4)*rr(4)
         ss(4) = (e(l)-rr(3))*      rr(4)  + rr(2)*rr(5)
         ss(5) = (e(l)-rr(1))*      rr(5)  + rr(2)*rr(4)
         ss(6) = (e(l)-rr(1))*(e(l)-rr(3)) - rr(2)*rr(2)

         if ( dabs(ss(1)).ge.dabs(ss(3)) ) then
            if ( dabs(ss(1)).lt.dabs(ss(6)) ) then
               j = 3
            else
               j = 1
            endif
         else 
            if ( dabs(ss(3)).ge.dabs(ss(6)) ) then
               j = 2
            else
               j = 3
            endif
         endif

         sum = 0.0d0
         j = 3*(j-1)

         do i=1,3
            a(i,l) = ss(ip(i+j))
            sum = sum + ss(ip(i+j))*ss(ip(i+j))
         end do
   
         sum = dsqrt(sum)

         do i=1,3
            a(i,l) = a(i,l) / sum
         end do

      end do

      call axa(a,3,1,2)

      return
      end

c     eigenvectors two distinct roots
c     h=trace+g

      subroutine root2(a,rr,e,ip,m,m1,m2,m3)
      implicit double precision (a-h,o-z)
      dimension a(3,3),rr(6),e(3),ip(9)

      p = 0.0d0

      h = e(2)

      do i=1,3

         k = (i*i + i)/2
         d = dabs(rr(k) - h)

         if (d.ge.p) then
            j = i
            p = d
         endif

      end do

      p = 0.0d0
      d = 0.0d0
      l = 3*(j-1)

      do i=1,3

         k = ip(i+l)
         a(i,2) = 1.0d0

         if (i.ne.j) then
            a(i,m) = rr(k)
            p = p - a(i,m)
         else
            a(i,m) = rr(k) - h
         endif

         d = d + a(i,m)**2
      end do

      a(j,2) = p / a(j,m)

      d = dsqrt(d)
      p = dsqrt(a(1,2)**2 + a(2,2)**2 + a(3,2)**2)

      do i=1,3
         a(i,2) = a(i,2) / p
         a(i,m) = a(i,m) / d
      end do

      call axa(a,m1,m2,m3)

      return
      end

      subroutine axa(a,m1,m2,m3)
      implicit double precision (a-h,o-z)
      dimension a(3,3)

      a(1,m1) = a(2,m2)*a(3,m3) - a(2,m3)*a(3,m2)
      a(2,m1) = a(3,m2)*a(1,m3) - a(3,m3)*a(1,m2)
      a(3,m1) = a(1,m2)*a(2,m3) - a(1,m3)*a(2,m2)

      return
      end

      logical function chksec(isal,ncalf)
      implicit double precision (a-h,o-z)
      dimension isal(*)

      chksec = .false.

      icoil = 0
      ihelsh = 0

      do i=1,ncalf
         is = isal(i)
         if (is.le.1) ihelsh = ihelsh + 1
         if (is.eq.3) icoil = icoil + 1
      end do

      ncalft = ihelsh + icoil
      if (dble(ihelsh)/dble(ncalft).gt.0.5d0) chksec = .true.

      return
      end

      subroutine fident(coo1,coo2,iamin1,iamin2,icalf1,icalf2,
     &                  ncalf1,ncalf2,isal1,isal2,x,y,neq,ok)
      implicit double precision (a-h,o-z)
      parameter (mxhits = 1000)
      logical ohit,ok,chksec,oscnd
      dimension coo1(3,*),coo2(3,*),iamin1(*),iamin2(*),
     &          icalf1(6,*),icalf2(6,*),isal1(*),isal2(*)
      dimension ihit(mxhits,2), x(3,*),y(3,*)

      oscnd = .false.
      if (chksec(isal1,ncalf1).and.chksec(isal2,ncalf2)) oscnd = .true.

      neql = 5
      nhit = 0

      ind1 = 1
      ind2 = 1

      do while(.true.)

         if (ohit(iamin1,iamin2,ncalf1,ncalf2,ind1,ind2,neql,
     &            isal1,isal2,oscnd)) then

            nhit = nhit + 1

            if (nhit.le.mxhits) then
               ihit(nhit,1) = ind1
               ihit(nhit,2) = ind2
            endif

            ind1 = ind1 + neql
            ind2 = ind2 + neql

            if (ind1.gt.ncalf1-neql) goto 10

            if (ind2.gt.ncalf2-neql) then
               ind1 = ind1 + 1
               if (nhit.gt.0) then
                  ind2 = ihit(nhit,2)
               else
                  ind2 = 1
               endif
            endif

         else

            ind2 = ind2 + 1

            if (ind2.gt.ncalf2-neql) then
               ind1 = ind1 + 1
               if (nhit.gt.0) then
                  ind2 = ihit(nhit,2)
               else
                  ind2 = 1
               endif
            endif

         endif 

      end do

10    if (nhit*neql.gt.(ncalf1+ncalf2)/6.0d0) then

         ok = .true.
         neq = 0

         do i=1,nhit

            i1 = ihit(i,1)
            i2 = ihit(i,2)

            do j=1,neql

               if (i1.le.ncalf1.and.i2.le.ncalf2) then

                  neq = neq + 1

                  do k=1,3
                     y(k,neq) = coo1(k,icalf1(1,i1))*0.52917706d0
                     x(k,neq) = coo2(k,icalf2(1,i2))*0.52917706d0
                  end do

                  i1 = i1 + 1
                  i2 = i2 + 1

               else
                  return
               endif

            end do

         end do

      else
         ok = .false.
      endif

      return
      end

      logical function ohit(iamin1,iamin2,ncalf1,ncalf2,ind1,ind2,neql,
     &                      isal1,isal2,oscnd)
      implicit double precision (a-h,o-z)
      logical oscnd
      dimension iamin1(*),iamin2(*),isal1(*),isal2(*)

      ohit = .true.
      do i=1,neql
         i1 = ind1+i-1
         i2 = ind2+i-1

         if (i1.le.ncalf1.and.i2.le.ncalf2) then
            if (oscnd) then
               if (iamin1(i1).ne.iamin2(i2).or.
     &             (isal1(i1).gt.1.or.isal2(i2).gt.1)) then
                  ohit = .false.
                  return
               endif
            else
               if (iamin1(i1).ne.iamin2(i2)) then
                  ohit = .false.
                  return
               endif
            endif
         endif

      end do

      return
      end

