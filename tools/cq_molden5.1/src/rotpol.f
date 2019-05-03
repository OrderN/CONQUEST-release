c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     "rtmat" find the rotation matrix that converts from the local
c     coordinate system at each multipole site to the global system
c
c     Modified version for MOLDEN
c     Distributed with the permission of Prof. Jay Ponder
c
      subroutine rtmat(i,j,k,itype,a)
      implicit double precision (a-h,o-z)
      parameter (mxsite=300)
      parameter (lmax=4)
      parameter (lmaxf=(lmax+1)*(lmax+1))
      character*8 ctag
      common /multip/ qmom(lmaxf,mxsite),car(3,mxsite),ctag(mxsite),
     &                nsites
      dimension a(3,3)
c
c
c     rotation matrix elements for z- and x-axes, first z then x
c
      if (itype.eq.1) then
         dx = car(1,j) - car(1,i)
         dy = car(2,j) - car(2,i)
         dz = car(3,j) - car(3,i)
         r = dsqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = car(1,k) - car(1,i)
         dy = car(2,k) - car(2,i)
         dz = car(3,k) - car(3,i)
         dotxz = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         dx = dx - dotxz*a(1,3)
         dy = dy - dotxz*a(2,3)
         dz = dz - dotxz*a(3,3)
         r = dsqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     rotation matrix elements for z- and x-axes, bisector method
c
      else

         dx = car(1,j) - car(1,i)
         dy = car(2,j) - car(2,i)
         dz = car(3,j) - car(3,i)
         r = dsqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = car(1,k) - car(1,i)
         dy = car(2,k) - car(2,i)
         dz = car(3,k) - car(3,i)
         r = dsqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = dsqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dotxz = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
         dx = dx2 - dotxz*a(1,3)
         dy = dy2 - dotxz*a(2,3)
         dz = dz2 - dotxz*a(3,3)
         r = dsqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
      endif
c
c     finally, find rotation matrix elements for the y-axis
c
      a(1,2) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
      a(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
      a(3,2) = a(2,1)*a(1,3) - a(1,1)*a(2,3)

      return
      end
      
c
c     "rotpole" computes the atomic multipole values in the global
c     coordinate frame by applying a rotation matrix to a set of
c     locally defined multipoles
c
c     Modified version for MOLDEN
c      
      subroutine rotpole (a,qrot)
      implicit double precision (a-h,o-z)
      dimension p2(3,3),r2(3,3),a(3,3),rpole(13),qrot(*)
c
c
c     monopoles have the same value in any coordinate frame
c
      rpole(1) = qrot(1)
c
c     rotate the dipoles
c
      do i= 2,4
         rpole(i) = 0.0d0
         do j= 2,4
            rpole(i) = rpole(i) + qrot(j)*a(i-1,j-1)
         end do
      end do
c
c     rotate the quadrupoles
c
      k = 5
      do i=1,3
         do j=1,3
            p2(i,j) = qrot(k)
            r2(i,j) = 0.0d0
            k = k + 1
         end do
      end do

      do i=1,3
         do j=1,3
            if (j.lt.i) then
               r2(i,j) = r2(j,i)
            else
               do k=1,3
                  do m=1,3
                     r2(i,j) = r2(i,j) + a(i,k)*a(j,m)*p2(k,m)
                  end do
               end do
            end if
         end do
      end do

      k = 5
      do i=1,3
         do j=1,3
            rpole(k) = r2(i,j)
            k = k + 1
         end do
      end do

c     Overwrite global multipoles in qrot  

      do i=1,13
         qrot(i) = rpole(i)
      end do

      return
      end
      
