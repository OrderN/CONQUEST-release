      subroutine rotmod(ipoint,ifav,coo)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (maxfat=1000)
      real fxyz,fc
      common /forcom/fxyz(3,maxfat),fc(3,maxfat) 
      common /athlp/ iatoms, mxnat

      common/rotxyz/amat(3,3),bmat(3,3),cmat(3,3),rmat(3,3),todeg
      parameter (maxpt=1000)
      common /savcor/ angs(3,maxpt),idon(maxpt)
      common /gamori/imap(numat1),cstand(3,numat1),istand(numat1),
     &               msucc
      dimension angle(3),cori(3,numat1),ftemp(3,numat1),coo(3,*)

      if (ipoint .gt. maxpt) return

      todeg=45.0d0/datan(1.0d0)
      delold=1.0d5

      amat(1,3)=0.0d0
      amat(2,3)=0.0d0
      amat(3,1)=0.0d0
      amat(3,2)=0.0d0
      amat(3,3)=1.0d0

      bmat(1,2)=0.0d0
      bmat(2,1)=0.0d0
      bmat(2,2)=1.0d0
      bmat(2,3)=0.0d0
      bmat(3,2)=0.0d0

      cmat(1,1)=1.0d0
      cmat(1,2)=0.0d0
      cmat(1,3)=0.0d0
      cmat(2,1)=0.0d0
      cmat(3,1)=0.0d0

c     put forces in correct order

      if (ifav.eq.1) then

         do i=1,iatoms
            do j=1,3
               ftemp(j,i) = fxyz(j,imap(i))
            end do 
         end do
   
         do i=1,iatoms
            do j=1,3
               fxyz(j,i) = ftemp(j,i)
            end do 
         end do

      endif

c     check if point already done

      if (idon(ipoint).eq.1) then

          do j=1,3
             angle(j) = angs(j,ipoint)
          end do

          goto 100

      endif

      if (ifav.eq.0) then

          do j=1,3
             angle(j) = angs(j,ipoint-1)
             angs(j,ipoint) = angs(j,ipoint-1)
          end do

          idon(ipoint) = 1
          goto 100

      endif

c----- start of medium search ----

      alco = angs(1,1)
      beco = angs(2,1)
      gaco = angs(3,1)

      do ia=1,21

         angle(1) = alco + dfloat(ia-1) - 10.0d0
         call rota(angle(1))

         do ib=1,21

            angle(2) = beco + dfloat(ib-1) - 10.0d0
            call rotb(angle(2))

            do ic=1,21

               angle(3) = gaco + dfloat(ic-1) - 10.0d0
               call rotc(angle(3))
               delta = rmomen(cori)

               if (delta.lt.delold) then
                 delold = delta
                 alopt  = angle(1)
                 beopt  = angle(2)
                 gaopt  = angle(3)
               endif

            end do

         end do

      end do

c----- end of medium search ------

c----- start of fine search ----

      alco = alopt
      beco = beopt
      gaco = gaopt

      do ia=1,21

         angle(1) = alco + dfloat(ia-1)*0.1d0 - 1.0d0
         call rota(angle(1))

         do ib=1,21

            angle(2) = beco + dfloat(ib-1)*0.1d0 - 1.0d0
            call rotb(angle(2))

            do ic=1,21

               angle(3) = gaco + dfloat(ic-1)*0.1d0 - 1.0d0
               call rotc(angle(3))
               delta=rmomen(cori)

               if (delta.lt.delold) then
                 delold = delta
                 alopt  = angle(1)
                 beopt  = angle(2)
                 gaopt  = angle(3)
               endif

            end do

         end do

      end do

c----- end of fine search ------
c     final angles

      angle(1)=alopt
      angle(2)=beopt
      angle(3)=gaopt

      idon(ipoint) = 1

      do j=1,3
         angs(j,ipoint) = angle(j)
      end do

100   continue

      call rota(angle(1))
      call rotb(angle(2))
      call rotc(angle(3))
      call rotcor(cori)

      do i=1,iatoms
         do j=1,3
            coo(j,i) = cori(j,i)
         end do
      end do

      return
      end
