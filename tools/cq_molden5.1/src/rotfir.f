      subroutine rotfid(ioxyz,ianz,coo)
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (maxpt=1000)
      parameter (bignum=1.0d12)
      parameter (small=1.0d-4)
      common /athlp/ iatoms, mxnat
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      common /rotxyz/amat(3,3),bmat(3,3),cmat(3,3),rmat(3,3),todeg
      common /gamori/imap(numat1),cstand(3,numat1),istand(numat1),
     &               msucc
      common /savcor/ angs(3,maxpt),idon(maxpt)
      character line*80
      dimension angfir(3)
      dimension vect1(3),vect2(3),vect3(3)
      dimension fdist(numat1),cdist(numat1),cori(3,numat1),
     &          cttt(3,numat1),coo(3,*),ianz(*)
 
      tol = small
      msucc = 1
cd     write(iun3,'(a)')'enter subroutine rotfir'

c     read in coords prior to orient. + forces first point

      rewind iun2
      ip = 1
      call gampoi(ip,istat,ioxyz)
      call docct
      
c     read in standard orientation

      rewind iun2
      call search(line,'molecular geometry',istat)
      if (istat.eq.0) then
         call inferr('No molecular geometry found !',1)
         goto 1000
      endif
      call search(line,'atom   atomic',istat)
      call readel(line,3)

      katoms = 0
5        read(iun2,'(a)',err=1000) line
         if (line(10:14).eq.'*****') goto 10
         if (line(11:70).eq.' ') goto 5
         katoms = katoms + 1
         read(line,'(22x,f5.1,3(f12.7,3x))',err=1000)
     &    charge,cstand(1,katoms),cstand(2,katoms),cstand(3,katoms)
         istand(katoms) = charge
         goto 5
10    continue
      
      if (katoms .ne. iatoms) then
         call inferr('rotfir: #atoms orientation prior != standard',1)
         goto 1000
      endif

c     establish the mapping between z-mat and standard orientation

      do i=1,iatoms
         if (i.le.numat1) then
            fdist(i) = vlen(coo(1,i))
            cdist(i) = vlen(cstand(1,i))
         endif
      end do

      if ((ianz(1).eq.istand(1)).and.(fdist(1).gt.tol).and.
     &   ((fdist(1)-cdist(1))**2.lt.tol)) then
         ist = 1
      elseif ((ianz(2).eq.istand(2)).and.
     &   ((fdist(1)-cdist(1))**2.lt.tol)) then
         ist = 2
      else
         call inferr('rotfir: mapping scheme doesnt work',1)
         goto 1000
      endif

c     rotate atom 1 or 2 to standard orientation

      call crprod(coo(1,ist),cstand(1,ist),vect1)

      if (vlen(vect1).gt.small) then

         call impsc(coo(1,ist),cstand(1,ist),co1)
         call crprod(vect1,coo(1,ist),vect2)
         call impsc(cstand(1,ist),vect2,si1)
         call tessa(vect1,co1,si1,fdist,coo,cori,iatoms,tol)

      else

         do i=1,iatoms
            do j=1,3
              cori(j,i) = coo(j,i)
            end do
         end do

      endif

c     now rotate every other atom of the prior orientation to 
c     the next atom of the standard orientation (atom 2 or 3)
c     and see which rotation gives the smallest error 

      delta = bignum
      do n=1,iatoms-ist
         call crprod(cstand(1,ist),cstand(1,ist+n),vect2)
         call crprod(vect2,cstand(1,ist),vect1)
         if (vlen(vect2).gt.small) goto 30
      end do

30    continue
      
      inum = 1

      do i=ist+1,iatoms
         if((ianz(i).eq.istand(ist+n)).and.
     &      ((fdist(i)-cdist(ist+n))**2.lt.tol)) then

c           establish rotation angle

            call crprod(cstand(1,ist),cori(1,i),vect2)
            if(vlen(vect2).gt.small) then
               call crprod(vect2,cstand(1,ist),vect3)
               call impsc(vect1,vect2,si2)
               call impsc(vect3,vect1,co2)

c           the rotatation part

            call tessa(cstand(1,ist),co2,si2,fdist,cori,cttt,iatoms,tol)
            
            else
                do k=1,iatoms
                   do j=1,3
                     cttt(j,k) = cori(j,k)
                   end do
                end do
            endif
c           the mapping part

            call caldis(dis,cttt,fdist,cdist,tol,ianz)

            if (dis.lt.delta) then
               delta = dis
               inum = i
            endif

         endif
      end do
       
c     activate rotation and mapping

      call crprod(cstand(1,ist),cori(1,inum),vect2)
      if(vlen(vect2).gt.small) then
         call crprod(vect2,cstand(1,ist),vect3)
         call impsc(vect1,vect2,si2)
         call impsc(vect3,vect1,co2)
         call tessa(cstand(1,ist),co2,si2,fdist,cori,cttt,iatoms,tol)
      else
         do k=1,iatoms
            do j=1,3
              cttt(j,k) = cori(j,k)
            end do
         end do
      endif
      call caldis(dis,cttt,fdist,cdist,tol,ianz)

c     put standard orientation in correct order

      do i=1,iatoms
         do j=1,3
            cori(j,i) = cstand(j,imap(i))
         end do 
      end do
  
      do i=1,iatoms
         do j=1,3
            cstand(j,i) = cori(j,i)
         end do 
      end do

c     end of mapping start finding rotation angles first point

      todeg     = 45.0d0 / datan(1.0d0)
      delold    = 1.0d5

      amat(1,3) = 0.0d0
      amat(2,3) = 0.0d0
      amat(3,1) = 0.0d0
      amat(3,2) = 0.0d0
      amat(3,3) = 1.0d0

      bmat(1,2) = 0.0d0
      bmat(2,1) = 0.0d0
      bmat(2,2) = 1.0d0
      bmat(2,3) = 0.0d0
      bmat(3,2) = 0.0d0

      cmat(1,1) = 1.0d0
      cmat(1,2) = 0.0d0
      cmat(1,3) = 0.0d0
      cmat(2,1) = 0.0d0
      cmat(3,1) = 0.0d0


c----- start of coarse search ----

      do ia=0,360,10
         angfir(1)=dfloat(ia)
         call rota(angfir(1))
         do ib=0,360,10
            angfir(2)=dfloat(ib)
            call rotb(angfir(2))
            do ic=0,360,10
               angfir(3)=dfloat(ic)
               call rotc(angfir(3))
               delta=distot()
               if(delta.lt.delold)then
                 delold=delta
                 alopt=angfir(1)
                 beopt=angfir(2)
                 gaopt=angfir(3)
               endif
            end do
         end do
      end do

c----- end of coarse search ------

c----- start of medium search ----

      alco=alopt
      beco=beopt
      gaco=gaopt
      do ia=1,21
         angfir(1)=alco+dfloat(ia-1)-10.0d0
         call rota(angfir(1))
         do ib=1,21
            angfir(2)=beco+dfloat(ib-1)-10.0d0
            call rotb(angfir(2))
            do ic=1,21
               angfir(3)=gaco+dfloat(ic-1)-10.0d0
               call rotc(angfir(3))
               delta=distot()
               if(delta.lt.delold)then
                 delold=delta
                 alopt=angfir(1)
                 beopt=angfir(2)
                 gaopt=angfir(3)
               endif
            end do
         end do
      end do

c----- end of medium search ------

c----- start of fine search ----

      alco=alopt
      beco=beopt
      gaco=gaopt
      do ia=1,21
         angfir(1)=alco+dfloat(ia-1)*0.1d0-1.0d0
         call rota(angfir(1))
         do ib=1,21
            angfir(2)=beco+dfloat(ib-1)*0.1d0-1.0d0
            call rotb(angfir(2))
            do ic=1,21
               angfir(3)=gaco+dfloat(ic-1)*0.1d0-1.0d0
               call rotc(angfir(3))
               delta=distot()
               if(delta.lt.delold)then
                 delold=delta
                 alopt=angfir(1)
                 beopt=angfir(2)
                 gaopt=angfir(3)
               endif
            end do
         end do
      end do

c----- end of fine search ------
c     final angles

      angfir(1) = alopt
      angfir(2) = beopt
      angfir(3) = gaopt
      angs(1,1) = alopt
      angs(2,1) = beopt
      angs(3,1) = gaopt
      idon(1)   = 1

cd     write(iun3,'(a)')'leave subroutine rotfir'
      return
1000  msucc = 0
      call inferr('ERROR in rotation to Standard Orientation!',1)
cd     write(iun3,'(a)')'leave subroutine rotfir'
      return
      end
