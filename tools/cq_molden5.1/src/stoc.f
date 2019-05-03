      subroutine stoc(maxnz,nz,ipart,imn,imx,ianz,iz,bl,alph,bet,
     $           ottest,natoms,ian,c,cz,imap,alpha,beta,
     $           oerror,oconv,odum)
c
c     imap(zmatline) = xyzline
c
      implicit double precision (a-h,p-z), logical (o)
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      character*80 errmsg
      dimension ianz(*),iz(4,*),bl(*),alpha(*),beta(*),ian(*),
     $          c(*),cz(*),alph(*),bet(*),imap(*)
      dimension u1(3),u2(3),u3(3),u4(3),vj(3),vp(3),v3(3)
      dimension t1(3),t2(3),t3(3),s1(3),s2(3),s3(3),s4(3)
      logical door
c
      data dzero/0.0d0/,done/1.0d0/,two/2.0d0/
      data tenm5/1.0d-5/,tenm6/1.0d-6/
      data f180/180.0d0/,four/4.0d0/
      data tenm10/1.0d-10/
      data tetdat/109.471d0/, toldat/0.001d0/, three/3.d0/
c
 2001 format(1x,i3,' tetrahedral angles replaced')
c
      otest(ccc) =  dabs(ccc-tettst) .lt. tettol
      errmsg = 'incipient floating point error detected !'

      door = .true.
      id = 0
      do i=1,3
         if (cz((i-1)*3+1).eq.0.0d0.and.cz((i-1)*3+2).eq.0.0d0
     &       .and.cz((i-1)*3+3).eq.0.0d0) id = id +1
      end do
      if (id.ge.2) door = .false.
c
c
c     check for potential overflow.
c
      if (oconv) then
         conver = 1.0d0/0.52917706d0
      else
         conver = 1.0d0
      endif

      if (ipart.ne.0) then
         il = imn
         if (il.lt.2) il = 2
         ih = imx
      else
         il = 2
         ih = nz
      endif
      
      if (nz.gt.maxnz) then
         oerror = .true.
         call zmterr('z-matrix cards is greater than the maximum !',
     &               -1,0,1)
         return
      endif
      oerror = .false.
      if (ih .lt. 2) goto 14

      do i=il,ih
         if (i.gt.3.and.iabs(iz(4,i)).gt.1) then
            oerror = .true.
            call zmterr('invalid beta angle type (z4) !',i,-4,1)
            return
         endif
         if (i.eq.2.and.iz(1,2).eq.1) goto 9
         if (i.eq.3.and.iz(1,3).lt.i.and.iz(2,3).lt.i.and.
     &       iz(1,3).gt.0.and.iz(2,3).gt.0) goto 9
         if (i.gt.3.and.iz(1,i).lt.i.and.iz(2,i).lt.i.and.iz(3,i).lt.i
     &       .and.iz(1,i).gt.0.and.iz(2,i).gt.0.and.iz(3,i).gt.0) goto 9
            oerror = .true.
            l = 0
            if (iz(1,i).ge.i) l = -1
            if (iz(2,i).ge.i) l = -2
            if (iz(3,i).ge.i) l = -3
            call zmterr('reference made to an undefined center !',i,l,1)
            return
    9    continue
         if (i.eq.3.and.iz(1,i).eq.iz(2,i)) then
             oerror = .true.
             call zmterr(
     &      'multiple references to a center on the same card !',i,0,1)
             return
         endif
         if (i.gt.3.and.(iz(1,i).eq.iz(2,i).or.iz(1,i).eq.iz(3,i).or.
     &        iz(2,i).eq.iz(3,i))) then
            oerror = .true.
            call zmterr(
     &     'multiple references to a center on the same card !',i,0,1)
            return
         endif
      end do
   14 continue

      pi = four * datan(done)
      torad = pi/f180
c
c     set up for laundering tetrahedral angles.
c     this feature is only invoked when test=.true..
c
      tetang = dacos(-done/three)
      tettst = tetdat * torad
      tettol = toldat * torad
c
c     zero temporary coordinate array cz
c
      if (door) then
         do l=1,3
            t1(l) = cz(l)
            t2(l) = cz(3+l)
            t3(l) = cz(6+l)
         end do
      endif

      nz3 = 3*nz
      if (ipart.eq.0) call vclr(cz,1,nz3)
c
c     move angles to local arrays and optionally test for
c     tetrahedral angles
c     test alpha for out of range 0 to 180 degrees
c     test for negative bond lengths.
c
      numtet = 0

      do i=il,ih

         alpha(i) = alph(i) * torad
         beta(i)  = bet(i) * torad

         if (bl(i).le.dzero .and. i.ne.1) then
            oerror = .true.
            call zmterr('negative bond length !',i,1,1)
            return
         endif

         if (i.gt.2 .and..not.(alpha(i).ge.dzero.and.alpha(i).le.pi))
     &   then
            oerror = .true.
            call zmterr(
     &     'angle alpha is outside the valid range of 0 to 180 !',i,2,1)
            return
         endif
         if (ottest .and. i.gt.3) then
            if (otest(alpha(i))) then
               alpha(i) = tetang
               alph(i)  = tetdat
               numtet = numtet + 1
            endif
            if (otest(beta(i))) then
               beta(i) = tetang
               bet(i)  = tetdat
               numtet = numtet + 1
            endif
            if (iz(4,i).ne.0) then
               if (.not.(beta(i).ge.dzero .and. beta(i).le.pi))
     &         then
                  oerror = .true.
                  call zmterr(
     &    'Dihedral is outside the valid range of 0 to 180 !',i,3,1)
                  return
               endif
            endif
         endif

      end do

c      if( (numtet.ne.0) .and. (iun3.ne.0) ) write(iun3,2001)numtet
c
c     z-coordinate, atom 2.
c
      if (door) then
         do l=1,3
            cz(l) = t1(l)
            s1(l) = t2(l) - t1(l)
         end do
         call vsc1(s1,1.0d0,1.0d-4)
         do l=1,3
            cz(3+l) = t1(l)+s1(l)*bl(2)*conver
         end do
      else
         cz(6) = bl(2)*conver
      endif

      if (ih.lt.3) goto 260
c
c     x-coordinate, center 3.
c
      if (door) then
         if (iz(1,3).eq.1) then
            do l=1,3
               s2(l) = t3(l) - t1(l)
            end do
         else
            do l=1,3
               s2(l) = t3(l) - t2(l)
               s1(l) = -s1(l)
            end do
         endif

         call vsc1(s2,1.0d0,1.0d-4)
         call crprod(s1,s2,s3)
         call crprod(s1,s3,s4)
         call vsc1(s4,1.0d0,1.0d-4)

         do l=1,3
            s3(l) = -dsin(alpha(3))*s4(l) + dcos(alpha(3))*s1(l)
         end do

         call vsc1(s3,bl(3)*conver,1.0d-4)
         if (iz(1,3).eq.1) then
            do l=1,3
               cz(6+l) = cz(l) + s3(l)
            end do
         else
            do l=1,3
               cz(6+l) = cz(3+l) + s3(l)
            end do
         endif
      else
         cz(7) = bl(3)*conver* dsin(alpha(3))

         if (iz(1,3).eq.1) then
c
c     z-coordinate, center 3.
c
            cz(9) = bl(3)*conver* dcos(alpha(3))
         else
c
c     z-coordinate on center 3 as a function of z-coordinate, center 2.
c
            cz(9) = cz(6)-bl(3)*conver* dcos(alpha(3))
         endif
      endif
c
c     beware of linear molecule.
c
      if (ih.lt.4) goto 260

      if (door) then 
         if (ipart.ne.0) then
            i = imn
            if (i.lt.4) i = 4
         else
            i = 4
         endif
         goto 90
      endif

      if (ipart.ne.0) then
         il1 = imn
         if (il1.lt.4) il1 = 4
         ih = imx
      else
         il1 = 4
         ih = nz
      endif

      do i=il1,ih
         ind3 = (i-1)*3
         if (dabs(cz(1+ind3-3)).ge.tenm5) goto 90
         cz(1+ind3) = bl(i)*conver* dsin(alpha(i))
         itemp = (iz(1,i)-1)*3
         jtemp = (iz(2,i)-1)*3
         cz(3+ind3) = cz(3+itemp)-bl(i)*conver* dcos(alpha(i))*
     $                dsign(done,cz(3+itemp)-cz(3+jtemp))
      end do

   90 k=i
      if (k.gt.ih) goto 260
      do j=k,ih
         jnd3 = (j-1)*3
         dcaj = dcos(alpha(j))
         dsaj = dsin(alpha(j))
         dcbj = dcos(beta(j))
         dsbj = dsin(beta(j))
         if (iz(4,j).eq.0) then
   
            call vec(tenm6,ovec,u1,cz,iz(2,j),iz(3,j))
            if (ovec) then
               oerror = .true.
               call zmterr(errmsg,i,0,1)
               return
            endif
   
            call vec(tenm6,ovec,u2,cz,iz(1,j),iz(2,j))
            if (ovec) then
               oerror = .true.
               call zmterr(errmsg,i,0,1)
               return
            endif
   
            call crprod(u1,u2,vp)
   
            arg = done - (u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))**2
            if (arg .lt. dzero) then
               oerror = .true.
               call zmterr(errmsg,i,0,1)
               return
            endif
   
            r = dsqrt(arg)
            if (r .lt. tenm6) then
               oerror = .true.
               call zmterr(errmsg,i,0,1)
               return
            endif
   
            do i=1,3
               u3(i) = vp(i)/r
            end do
            call crprod(u3,u2,u4)
            do i=1,3
               vj(i)      = bl(j)*conver*(-u2(i)*dcaj+u4(i)*dsaj*dcbj
     &                      +u3(i)*dsaj*dsbj)
               itemp      = (iz(1,j)-1)*3
               cz(i+jnd3) = vj(i)+cz(i+itemp)
            end do
   
         else
   
            if (iabs(iz(4,j))-1.eq.0) then
   
               call vec(tenm6,ovec,u1,cz,iz(1,j),iz(3,j))
   
               if (ovec) then
                   oerror = .true.
                   call zmterr(errmsg,i,0,1)
                   return
               endif
   
               call vec(tenm6,ovec,u2,cz,iz(2,j),iz(1,j))
               if (ovec) then
                  oerror = .true.
                  call zmterr(errmsg,i,0,1)
                  return
               endif
   
               czeta=-(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))
               denom = done - czeta ** 2
   
               if ( dabs(denom) .lt. tenm6) then
                  oerror = .true.
                  call zmterr(errmsg,i,0,1)
                  return
               endif
   
               aj    = (-dcbj+czeta*dcaj)/denom
               bj    = (dcaj-czeta*dcbj)/denom
               r     = dzero
               gamma = pi/two
               if (dabs(czeta) .ge. tenm6) then
                  if (czeta .lt. dzero) r = pi
                  if (denom .lt. dzero) then
                     oerror = .true.
                     call zmterr(errmsg,i,0,1)
                     return
                  endif
                  gamma = datan(dsqrt(denom)/czeta)+r
               endif
   
               dj = dzero
   
               if (dabs(gamma+alpha(j)+beta(j)-two*pi)
     &             -tenm6.ge.0.0d0) then
                  arg = (done + aj*dcbj - bj*dcaj)  /  denom
                  if (arg .lt. dzero) then
                     oerror = .true.
                     call zmterr(errmsg,i,0,1)
                     return
                  endif
                  dj =  dfloat(iz(4,j)) * dsqrt(arg)
               endif
   
               call crprod(u1,u2,v3)
               do i=1,3
                  u3(i)      = aj*u1(i)+bj*u2(i)+dj*v3(i)
                  vj(i)      = bl(j)*conver*u3(i)
                  itemp      = (iz(1,j)-1)*3
                  cz(i+jnd3) = vj(i)+cz(i+itemp)
               end do
   
            else
   
               call vec(tenm6,ovec,u1,cz,iz(1,j),iz(3,j))
               if (ovec) then
                  oerror = .true.
                  call zmterr(errmsg,i,0,1)
                  return
               endif
   
               call vec(tenm6,ovec,u2,cz,iz(2,j),iz(1,j))
               if (ovec) then
                  oerror = .true.
                  call zmterr(errmsg,i,0,1)
                  return
               endif
   
               czeta=-(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))
               call crprod(u1,u2,v3)
               v3mag = dsqrt(v3(1)*v3(1)+v3(2)*v3(2)+v3(3)*v3(3))
               denom = done - czeta**2
   
               if (dabs(denom) .le. tenm6) then
                  oerror = .true.
                  call zmterr(errmsg,i,0,1)
                  return
               endif
   
               aj = v3mag*dcbj / denom
               arg = (done-dcaj*dcaj-aj*dcbj*v3mag) / denom
               if (arg .lt. dzero) then
                  oerror = .true.
                  call zmterr(errmsg,i,0,1)
                  return
               endif
   
               bj = dsqrt(arg)
   
               if (iz(4,j)-2.ne.0) bj = -bj
   
               dj = bj*czeta+dcaj
               do i=1,3
                  u3(i)      = bj*u1(i)+dj*u2(i)+aj*v3(i)
                  vj(i)      = bl(j)*conver*u3(i)
                  itemp      = (iz(1,j)-1)*3
                  cz(i+jnd3) = vj(i)+cz(i+itemp)
               end do
   
            endif
   
         endif

      end do
c
c     eliminate dummy atoms.  dummy atoms are characterized by
c     negative atomic numbers.  ghost atoms have zero atomic
c     numbers.  ghost atoms are not eliminated.
c
  260 continue

      if (ipart.ne.0) then
         iaind = (il-1)*3
         do i=il,ih
            imap(i)     = i
            ian(i) = ianz(i)
            c(1+iaind)  = cz(1+iaind)
            c(2+iaind)  = cz(2+iaind)
            c(3+iaind)  = cz(3+iaind)
            iaind       = iaind+3
         end do
         return
      endif

      natoms=0
      iaind = 0
      naind = 0
      do i=1,nz
         imap(i) = 0
         if (.not.(.not.odum.and.ianz(i).eq.99)) then
            natoms      = natoms+1
            imap(i)     = natoms
            ian(natoms) = ianz(i)
            c(1+naind)  = cz(1+iaind)
            c(2+naind)  = cz(2+iaind)
            c(3+naind)  = cz(3+iaind)
            naind       = naind+3
         endif
         iaind       = iaind+3
      end do

      if (ipart.ne.0) return
c
c     check for superposed atoms
c
      do i=nz,1,-1
         do j=1,nz
            if (i.ne.j.and.imap(i).ne.0.and.imap(j).ne.0) then
               if (dist2(c((imap(i)-1)*3+1),c((imap(j)-1)*3+1))
     &             .lt.tenm6) then
                   call zmterr('Superimposed atoms !',i,0,0)
                   print*,'Zmat to xyz: Superimposed atoms:'
                   print*,'atoms ',j,',',i,' ',ian(j),' ',ian(i)
                   print*,'mapatoms ',imap(j),',',imap(i),
     &                     ian(imap(j)),ian(imap(i))
                   goto 285
               endif
            endif
         end do
      end do
285   continue
c
c     'tidy' up the coordinates.
c
      nat3 = 3*natoms
      do i=1,nat3
         if (dabs(c(i)).le.tenm10) then
            c(i) = dzero
         endif
      end do

      return
      end

      subroutine stocc(maxnz,nz,iz,bl,alph,bet,
     $           alpha,beta,ierror)

c     this routine checks the number of error that would be made by stoc

      implicit double precision (a-h,p-z), logical (o)
      dimension iz(4,*),bl(*),alpha(*),beta(*),
     $          alph(*),bet(*)
c
      data dzero/0.0d0/,done/1.0d0/
      data f180/180.0d0/,four/4.0d0/
      data tetdat/109.471d0/, toldat/0.001d0/, three/3.d0/
c
      otest(ccc) =  dabs(ccc-tettst) .lt. tettol
c
      ierror = 0
      if (nz.gt.maxnz) ierror = ierror + 1
      if (nz .lt. 2) goto 14
      do i=2,nz
         if (i.gt.3.and.iabs(iz(4,i)).gt.1) ierror = ierror + 1
         if (i.eq.2.and.iz(1,2).eq.1) goto 9
         if (i.eq.3.and.iz(1,3).lt.i.and.iz(2,3).lt.i.and.
     &       iz(1,3).gt.0.and.iz(2,3).gt.0) goto 9
         if (i.gt.3.and.iz(1,i).lt.i.and.iz(2,i).lt.i.and.iz(3,i).lt.i
     &       .and.iz(1,i).gt.0.and.iz(2,i).gt.0.and.iz(3,i).gt.0) goto 9
            ierror = ierror + 1
    9    continue
         if (i.eq.3.and.iz(1,i).eq.iz(2,i)) then
            ierror = ierror + 1
         endif
         if (i.gt.3.and.(iz(1,i).eq.iz(2,i).or.iz(1,i).eq.iz(3,i).or.
     &        iz(2,i).eq.iz(3,i))) then
            ierror = ierror + 1
         endif
      end do
   14 continue
      pi=four* datan(done)
      torad=pi/f180
c
      tetang = dacos(-done/three)
      tettst = tetdat * torad
      tettol = toldat * torad
c
      numtet = 0
      do 20 i=1,nz
         alpha(i) = alph(i) * torad
         beta(i)  = bet(i) * torad
         if (bl(i).gt.dzero .or. i.eq.1) goto 15
            ierror = ierror + 1
   15    continue
         if (i.le.2 .or. (alpha(i).ge.dzero.and.alpha(i).le.pi)) goto 16
            ierror = ierror + 1
   16    continue
         if (.not. otest(alpha(i))) goto 17
            alpha(i) = tetang
            alph(i)  = tetdat
            numtet = numtet + 1
   17    if (.not. otest(beta(i))) goto 18
            beta(i) = tetang
            bet(i)  = tetdat
            numtet = numtet + 1
   18    continue
         if (iz(4,i).eq.0 .or. i.le.3) goto 20
         if (beta(i).ge.dzero .and. beta(i).le.pi) goto 20
            ierror = ierror + 1
   20 continue

      return
      end
