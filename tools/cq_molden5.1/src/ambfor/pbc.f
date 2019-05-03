      subroutine reddis(vr)
      implicit real (a-h,o-z), integer (i-n)
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      dimension vr(3)
      
      
      
      do i=1,3
         vr(i) = mod(vr(i),abc(i))
         if (vr(i).gt.abc2(i)) then
            vr(i) = vr(i) - abc(i)
         endif
         if (vr(i).lt.-abc2(i)) then
            vr(i) = vr(i) + abc(i)
         endif
      end do

      return
      end

      subroutine rddisr(vr)
      implicit real (a-h,o-z), integer (i-n)
      real abc,abc2,angles
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      dimension vr(3)
      
      do i=1,3
         do while (vr(i).gt.real(abc2(i)))
            vr(i) = vr(i) - real(abc(i))
         end do
         do while (vr(i).lt.-real(abc2(i)))
            vr(i) = vr(i) + real(abc(i))
         end do
      end do

      return
      end

      subroutine appbnd(coo,ityp)
      implicit real (a-h,o-z), integer (i-n)
      integer*2 ityp
      logical box,cell,fast,outbox
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /athlp/  iatoms, mxnat
      dimension coo(3,*),cw(3,3),vec(3),ityp(*)
      
c put water molecules which are entirely out of the box
c back in the box

      do i=1,iatoms

         if (ityp(i)  .eq.649.and.
     &       ityp(i+1).eq.650.and.
     &       ityp(i+2).eq.650) then

            outbox = .true.

            do j=1,3
               vec(j) = 0.0e0
            end do

            do k=1,3

               l = 0

               do j=1,3
                  cw(j,k) = coo(j,i+k-1)
                  if (cw(j,k).gt.abc2(j)) then
                     l = l + 1
                     if (k.eq.1) vec(j) = -abc(j)
                  endif
                  if (cw(j,k).lt.-abc2(j)) then
                     l = l + 1
                     if (k.eq.1) vec(j) = abc(j)
                  endif
               end do

               if (l.eq.0) outbox = .false.

            end do

c water out of the box, put it back

            if (outbox) then
               do k=1,3
                  do j=1,3
                     coo(j,i+k-1) = cw(j,k) + vec(j)
                  end do
               end do
            endif

         endif
      end do

      return
      end

      subroutine makbod(coo)
      implicit real (a-h,o-z), integer (i-n)
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      dimension coo(3,*),vec(3)

c     get largest diameter protein and add a default
c     set abc and abc2
c     set protein in center of the box

      offs = 14.0e0 

      call docnt(vec,coo)

      do i=1,3
         abc(i) = 2.0e0*vec(i) + offs
         abc2(i) = 0.5e0*abc(i)
      end do

      return
      end

      subroutine allbox
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxion=2000)
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /h2oer/ numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      common /athlp/  iatoms, mxnat
      dimension nbox(3)

      watd = 18.620e0
      nwater = 648/3

      do i=1,3
         t = abc(i)/watd
         nbox(i) = int(t) + 1
      end do

      newat = nbox(1)*nbox(2)*nbox(3)*nwater*3

      call allcoo(newat)
      if (nion.gt.0) call allq

      return
      end

      subroutine filbod(water,coo,iconn,iresid,ityp,
     &                  nac,iac,nad,
     &                  nbnd,ibnd,bl,bk,
     &                  nang,iang,ango,ak,q,iopt,iwtpr)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxac=3*mxcon)
      parameter (numres=50000)
      parameter (mxion=2000)
      logical box,cell,fast,chkwat,chkbox,chkion
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /athlp/  iatoms, mxnat
      common /residu/  qtot,ihsres,
     &                 ires(numres),ibeg(numres),iend(numres)
      common /h2oer/   numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      integer*2 ityp
      dimension water(3,*), coo(3,*), nbox(3), doff(3)
      dimension ityp(*),iconn(mxcon+1,*),iresid(*)
      dimension nac(*),nad(*),iac(mxac,*)
      dimension ibnd(2,*), bl(*), bk(*)
      dimension iang(3,*),ango(*),ak(*)
      dimension q(*),iopt(*),iwtpr(*),cwat(3,3),ibox(3)

      watd = 18.620e0
      nwater = 648/3

c 18.620 ang box
c     take prefilled water box and tile the real box with this
c     sub box
c     remove water that overlap with protein

      do i=1,3
         t = abc(i)/watd
         nbox(i) = int(t) + 1
      end do

      ishoh = ires(ihsres)
      if (ishoh.gt.0) then
          ishoh = -4
      else
          ishoh = ishoh - 1
      endif

      natnow = iatoms
      ioff = iatoms
      numwat = 0
      nwat = 0

      ibox(1) = 0

      do i=1,nbox(1)
         if (i.eq.nbox(1)) ibox(1) = 1
         ibox(2) = 0

         do j=1,nbox(2)
            if (j.eq.nbox(2)) ibox(2) = 1
            ibox(3) = 0

            do k=1,nbox(3)
               if (k.eq.nbox(3)) ibox(3) = 1

               doff(1) = -abc2(1) + 0.5e0*watd + watd*dble(i-1)
               doff(2) = -abc2(2) + 0.5e0*watd + watd*dble(j-1)
               doff(3) = -abc2(3) + 0.5e0*watd + watd*dble(k-1)

               chkbox = .false.
               if (i.eq.nbox(1).or.j.eq.nbox(2).or.k.eq.nbox(3))
     &             chkbox = .true.

               do n=1,nwater

                  in = (n-1)*3

                  do m=1,3
                     cwat(m,1) = water(m,in+1) + doff(m)
                     cwat(m,2) = water(m,in+2) + doff(m)
                     cwat(m,3) = water(m,in+3) + doff(m)
                  end do

                  if (chkwat(cwat,coo,ityp,iopt,ioff,ibox,
     &                       iwt,chkbox)) then


                     nwat = nwat + 1

                     if (.not.chkion(nwat)) then

c ADD WATER

                        numwat = numwat + 1
                        iwtpr(numwat) = iwt

                        if (ihsres.lt.numres) then
                           ihsres = ihsres + 1
                        endif

                        ires(ihsres) = ishoh
                        ibeg(ihsres) = ioff+1
                        iend(ihsres) = ioff+3
                        iresid(ioff+1) = ishoh
                        iresid(ioff+2) = ishoh
                        iresid(ioff+3) = ishoh
                        ishoh = ishoh - 1

                        do m=1,3
                           coo(m,ioff+1) = cwat(m,1)
                           coo(m,ioff+2) = cwat(m,2)
                           coo(m,ioff+3) = cwat(m,3)
                        end do
   
                        ityp(ioff+1) = 649
                        ityp(ioff+2) = 650
                        ityp(ioff+3) = 650

                        iwtpr(numwat) = 0
   
                        iconn(1,ioff+1) = 2
                        iconn(2,ioff+1) = ioff+2
                        iconn(3,ioff+1) = ioff+3
   
                        iconn(1,ioff+2) = 1
                        iconn(2,ioff+2) = ioff+1
   
                        iconn(1,ioff+3) = 1
                        iconn(2,ioff+3) = ioff+1

                        ibnd(1,nbnd+1) = ioff+1
                        ibnd(2,nbnd+1) = ioff+2
                        bl(nbnd+1)     = 0.9572e0
                        bk(nbnd+1)     = 553.0e0

                        ibnd(1,nbnd+2) = ioff+1
                        ibnd(2,nbnd+2) = ioff+3
                        bl(nbnd+2)     = 0.9572e0
                        bk(nbnd+2)     = 553.0e0
                        nbnd = nbnd + 2

                        nang = nang + 1
                        iang(1,nang) = ioff+2
                        iang(2,nang) = ioff+1
                        iang(3,nang) = ioff+3
                        ango(nang)   = 104.52e0
                        ak(nang)     = 100.00e0

                        q(ioff+1)    = -0.8340e0
                        q(ioff+2)    = 0.4170e0
                        q(ioff+3)    = 0.4170e0
   
                        if (i.eq.1.or.j.eq.1.or.k.eq.1) then
                           iopt(ioff+1) = 1
                           iopt(ioff+2) = 1
                           iopt(ioff+3) = 1
                        else
                           iopt(ioff+1) = 0
                           iopt(ioff+2) = 0
                           iopt(ioff+3) = 0
                        endif
   
                        nac(ioff+1) = 0
                        nac(ioff+2) = 1
                        nac(ioff+3) = 1
                        iac(1,ioff+2) = ioff+3
                        iac(1,ioff+3) = ioff+2
                        nad(ioff+1) = 0
                        nad(ioff+2) = 0
                        nad(ioff+3) = 0

                        ioff = ioff + 3

                     endif

                  endif

               end do
               
            end do
         end do
      end do

      do i=iatoms,ioff
         iopt(i) = 1
      end do

      iatoms = ioff

      return
      end

      logical function chkion(iwat)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxion=2000)
      common /h2oer/  numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat

      chkion = .false.
      if (nion.eq.0) return

      do i=1,nion
         if (iwat.eq.iw(i)) then
            chkion = .true.
            return
         endif 
      end do

      return
      end

      subroutine addions(coo,iconn,iresid,ityp,q,iopt)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (numres=50000)
      parameter (mxion=2000)
      common /athlp/  iatoms, mxnat
      common /residu/  qtot,ihsres,
     &                 ires(numres),ibeg(numres),iend(numres)
      common /h2oer/   numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      integer*2 ityp
      dimension coo(3,*), ityp(*),iconn(mxcon+1,*),iresid(*)
      dimension q(*),iopt(*)

c ADD ION

      ishoh = ires(ihsres)
      if (ishoh.gt.0) then
          ishoh = -4
      else
          ishoh = ishoh - 1
      endif

      ioff = iatoms

      do i=1,nion
         ihsres = ihsres + 1
         ires(ihsres) = ishoh
         ibeg(ihsres) = ioff+1
         iend(ihsres) = ioff+1
         iresid(ioff+1) = ishoh

         ishoh = ishoh - 1

         iatptr = iatoms + (iw(i)-1)*3 + 1

         do j=1,3
            coo(j,ioff+1) = coo(j,iatptr)
         end do

         if (iontyp.eq.1) then
            ityp(ioff+1) = 659
            q(ioff+1)    = -1.0e0
         else
            ityp(ioff+1) = 652
            q(ioff+1)    = +1.0e0
         endif

         iconn(1,ioff+1) = 0
         iopt(ioff+1) = 1

         ioff = ioff + 1

      end do

      iatoms = ioff

      return
      end

      subroutine cntwad(water,coo,iconn,iresid,ityp,q,iopt,qpot)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxac=3*mxcon)
      parameter (numres=50000)
      parameter (mxion=2000)
      logical doind,exclu,box,cell,fast,chkwat,chkbox
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /athlp/  iatoms, mxnat
      common /residu/  qtot,ihsres,
     &                 ires(numres),ibeg(numres),iend(numres)
      common /h2oer/   numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      integer*2 ityp
      dimension water(3,*), coo(3,*), nbox(3), doff(3)
      dimension ityp(*),iconn(mxcon+1,*),iresid(*)
      dimension q(*),iopt(*),cwat(3,3),c1(3),c2(3),ibox(3),qpot(*)

      watd = 18.620e0
      nwater = 648/3

c 18.620 ang box
c     take prefilled water box and tile the real box with this
c     sub box
c     remove water that overlap with protein

      do i=1,3
         t = abc(i)/watd
         nbox(i) = int(t) + 1
      end do

      ioff = iatoms
      nitmp = 0

      ibox(1) = 0

      do i=1,nbox(1)
         if (i.eq.nbox(1)) ibox(1) = 1
         ibox(2) = 0

         do j=1,nbox(2)
            if (j.eq.nbox(2)) ibox(2) = 1
            ibox(3) = 0

            do k=1,nbox(3)
               if (k.eq.nbox(3)) ibox(3) = 1

               doff(1) = -abc2(1) + 0.5e0*watd + watd*dble(i-1)
               doff(2) = -abc2(2) + 0.5e0*watd + watd*dble(j-1)
               doff(3) = -abc2(3) + 0.5e0*watd + watd*dble(k-1)

               chkbox = .false.
               if (i.eq.nbox(1).or.j.eq.nbox(2).or.k.eq.nbox(3))
     &             chkbox = .true.

               do n=1,nwater

                  in = (n-1)*3

                  do m=1,3
                     cwat(m,1) = water(m,in+1) + doff(m)
                     cwat(m,2) = water(m,in+2) + doff(m)
                     cwat(m,3) = water(m,in+3) + doff(m)
                  end do

                  if (chkwat(cwat,coo,ityp,iopt,ioff,ibox,
     &                       iwt,chkbox)) then

                     numwat = numwat + 1

                     if (i.eq.1.or.j.eq.1.or.k.eq.1) then
                           iopt(ioff+1) = 1
                           iopt(ioff+2) = 1
                           iopt(ioff+3) = 1
                     else
                           iopt(ioff+1) = 0
                           iopt(ioff+2) = 0
                           iopt(ioff+3) = 0
                     endif

                     ityp(ioff+1) = 649
                     ityp(ioff+2) = 650
                     ityp(ioff+3) = 650

                     do m=1,3
                        coo(m,ioff+1) = cwat(m,1)
                        coo(m,ioff+2) = cwat(m,2)
                        coo(m,ioff+3) = cwat(m,3)
                     end do

                     if (ionpl.eq.1) then
                        call clmond(cwat(1,1),pot,coo,q)
                        qpot(numwat) = pot
                     endif

                     ioff = ioff + 3

                  endif

               end do
               
            end do
         end do
      end do

      if (ionpl.eq.1) then

         nitmp = 0

         do while (nitmp.lt.nion)

            doind = .false.
            potl = 1.0e10
            iptr = 0

            do i=1,numwat

               pot = qpot(i)

               exclu = .false.
               do j=1,nitmp
                  if (iw(j).eq.i) exclu = .true.
               end do

               if (pot.lt.potl.and..not.exclu) then
                  doind = .true.
                  potl = pot
                  iptr = i
               endif

            end do

            if (doind) then
               nitmp = nitmp + 1
               iw(nitmp) = iptr

               do i=1,numwat
                  iatptr = iatoms + (i-1)*3 + 1
                  do j=1,3
                     c1(j) = coo(j,iatptr)
                  end do

                  jatptr = iatoms + (iptr-1)*3 + 1
                  do j=1,3
                     c2(j) = coo(j,jatptr)
                  end do
                  r2 = dist2(c1,c2)
                  r = sqrt(r2)
                  qpot(i) = qpot(i) + 1.0e0/r

               end do

            endif

         end do

      else

         call crionp

      endif

      call addions(coo,iconn,iresid,ityp,q,iopt)

      return
      end

      logical function chkvec(vec,cwat,coo,iopt,ityp,ioff)
      implicit real (a-h,p-w),integer (i-n),logical (o)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      integer ambcls
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxcls=49)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      integer amb2gf
      common /gftyp/  gfvdw(2,mxgff),amb2gf(mxamb)
      common /athlp/  iatoms, mxnat
      integer*2 ityp
      dimension v(3),vec(3),cwat(3,3),coo(3,*),ityp(*)
      dimension iopt(*)

      chkvec = .true.

      vdwat = 2.7683e0

      do i=iatoms,ioff
            
         if (iopt(i).eq.1) then

            i1 = int(ityp(i))
            if (i1.gt.0) then
               il = ambcls(i1)
               vdwr = ambvdwr(il)
            else
               il = iabs(i1)
               vdwr = gfvdw(1,il)
            endif


               do k=1,3
                  v(k) = coo(k,i) + vec(k) - cwat(k,1)
               end do

               rab2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
               dmaxsq = vdwr + vdwat
               dmaxsq = dmaxsq * dmaxsq

               if (rab2.lt.dmaxsq) then
                  chkvec = .false.
                  return
               endif


         endif

      end do

      return
      end

      logical function chkwat(cwat,coo,ityp,iopt,ioff,ibox,iwrpr,chkbox)
      implicit real (a-h,p-w),integer (i-n),logical (o)
      parameter (mxamb=1590)
      integer ambcls
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxcls=49)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      parameter (mxgff=72)
      integer amb2gf
      common /gftyp/  gfvdw(2,mxgff),amb2gf(mxamb)
      integer*2 ityp
      logical chkbox,chkvec
      common /athlp/  iatoms, mxnat
      logical box,cell,fast,owat
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      dimension coo(3,*),ityp(*),cwat(3,3),v(3),vdwat(3),ibox(3)
      dimension vec(3),iopt(*)

      chkwat = .false.


      owat = .true.
      iwrpr = 0

      if (chkbox) then

         do i=1,3
            do j=1,3
               if (cwat(j,i).gt.abc2(j)) owat = .false.
            end do
         end do

         if (.not.owat) return

         do i=1,3

            do j=1,3
               vec(j) = 0.0e0
            end do

            if (ibox(i).eq.1) then
               vec(i) = abc(i)
               if (.not.chkvec(vec,cwat,coo,iopt,ityp,ioff)) return
            endif

         end do

c         chkwat = .true.
c         return

      endif

c check overlap with protein atoms 

      owat = .true.

      vdwat(1) = 2.7683e0

      do i=1,iatoms
         
         i1 = int(ityp(i))
         if (i1.gt.0) then
            il = ambcls(i1)
            vdwr = ambvdwr(il)
         elseif (i1.le.0) then
            i1 = iabs(i1)
            vdwr = gfvdw(1,i1)
         endif

         do k=1,3
            v(k) = coo(k,i) - cwat(k,1)
         end do
         rab2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
         dmaxsq = vdwr + vdwat(1)
         dmaxsq = dmaxsq * dmaxsq
         if (rab2.lt.dmaxsq) owat = .false.
         if (rab2.lt.81.0e0) iwrpr = 1

      end do

      if (owat) chkwat = .true.

      return
      end

      subroutine docnt(vecm,coo)
      implicit real (a-h,p-w),integer (i-n),logical (o)
      common /athlp/ iatoms, mxnat
      dimension vec(3),vecm(3)
      dimension coo(3,*)

      do i=1,3
         vecm(i) = 0.0e0
      end do

      call cntvec(vec,coo,iatoms)

      do i=1,iatoms
         do j=1,3
            coo(j,i) = coo(j,i) - vec(j)
            da = abs(coo(j,i))
            if (da.gt.vecm(j)) vecm(j) = da
         end do
      end do

      return
      end

      subroutine cntvec(vec,coo,iatoms)
      implicit real (a-h,p-w),integer (i-n),logical (o)
      dimension vec(3), coo(3,*)

      do i=1,3
         vec(i) = 0.0e0
      end do

      if (iatoms.le.0) return

      do i=1,iatoms
         do j=1,3
            vec(j) = vec(j) + coo(j,i)
         end do
      end do
 
      do j=1,3
         vec(j) = vec(j) / dble(iatoms)
      end do

      return
      end

      subroutine setop(a,b,c,alpha,beta,gamma)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /cell/ xa,ya,yb,za,zb,zc
      real rxa,rya,ryb,rza,rzb,rzc
      common /cellr/ rxa,rya,ryb,rza,rzb,rzc

      torad = atan(1.0e0) / 45.0e0
      alpha = alpha*torad
      beta = beta*torad
      gamma = gamma*torad

      ca = cos(alpha)
      cb = cos(beta)
      cc = cos(gamma)
      sc = sin(gamma)
      sa = sin(alpha)
      cvol = a*b*c*
     &   sqrt(1-ca**2-cb**2-cc**2+2.0e0*ca*cb*cc)
      xa = cvol / (b*c*sa)
      ya = (a*(cc-(ca*cb)))/sa
      yb = (b*sa)
      za = (a*cb)
      zb = (b*ca)
      zc = c

      rxa = real(xa)
      rya = real(ya)
      ryb = real(yb)
      rza = real(za)
      rzb = real(zb)
      rzc = real(zc)

c      print*,' '
c      write(*,'(a,f9.3)') ' Cell Volume = ',cvol

      return
      end

      subroutine getcel(a,b,c,alpha,beta,gamma)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /cell/ xa,ya,yb,za,zb,zc

      torad = atan(1.0e0) / 45.0e0

      a = sqrt(xa**2+ya**2+za**2)
      b = sqrt(yb**2+zb**2)
      c = zc

      sa = yb/b
      ca = zb/b
      cb = za/a
      cc = ca*cb + ya*sa/a

      alpha = acos(ca) / torad
      beta  = acos(cb) / torad
      gamma = acos(cc) / torad

      return
      end

      subroutine fr2crt(vec)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /cell/ xa,ya,yb,za,zb,zc
      dimension vec(3)

      vec(3) = za*vec(1) + zb*vec(2) + zc*vec(3)
      vec(2) = ya*vec(1) + yb*vec(2)
      vec(1) = xa*vec(1)

      return
      end

      subroutine crt2fr(veci,veco)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /cell/ xa,ya,yb,za,zb,zc
      dimension veci(3),veco(3)

      veco(1) =  veci(1) / xa
      veco(2) = (veci(2) - ya*veco(1)) / yb
      veco(3) = (veci(3) - za*veco(1) - zb*veco(2)) / zc 

      return
      end

      subroutine fr2crr(vec)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /cellr/ xa,ya,yb,za,zb,zc
      dimension vec(3)

      vec(3) = za*vec(1) + zb*vec(2) + zc*vec(3)
      vec(2) = ya*vec(1) + yb*vec(2)
      vec(1) = xa*vec(1)

      return
      end

      subroutine crr2fr(veci,veco)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /cellr/ xa,ya,yb,za,zb,zc
      dimension veci(3),veco(3)

      veco(1) =  veci(1) / xa
      veco(2) = (veci(2) - ya*veco(1)) / yb
      veco(3) = (veci(3) - za*veco(1) - zb*veco(2)) / zc

      return
      end

      subroutine redfr(xyz)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      dimension xyz(3),fr(3)

      call crt2fr(xyz,fr)
      
      do i=1,3
         frt = 1.0e0
         if (fr(i).lt.0.0e0) then
             fr(i) = abs(fr(i))
             frt = -1.0e0
         endif

         fr(i) = mod(fr(i),1.0e0) 
         if (fr(i) .gt. 0.5e0) fr(i) = fr(i) - 1.0e0

         xyz(i) = frt*fr(i)
      end do

      call fr2crt(xyz)

      return
      end

      subroutine rddfrr(xyz)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      real abc,abc2,angles
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /debug/ ideb
      dimension xyz(3),fr(3)


      call crr2fr(xyz,fr)
      

      do i=1,3
         frt = 1.0e0
         if (fr(i).lt.0.0e0) then
             fr(i) = abs(fr(i))
             frt = -1.0e0
         endif

         fr(i) = amod(fr(i),1.0e0) 
         if (fr(i) .gt. 0.5e0) fr(i) = fr(i) - 1.0e0

         xyz(i) = frt*fr(i)
      end do


      call fr2crr(xyz)

      return
      end

      subroutine adcldr(cellder,ded,de,fct,lex)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /cell/ xa,ya,yb,za,zb,zc
      common /expnd/ addfr(3,125),ie(5),izero
      dimension ded(3),cellder(6),addf(3)

      do j=1,3
         addf(j) = addfr(j,lex)
      end do

      cellder(1) = cellder(1) + addf(1)*ded(1)*de*fct

      cellder(2) = cellder(2) + addf(1)*ded(2)*de*fct
      cellder(3) = cellder(3) + addf(2)*ded(2)*de*fct

      cellder(4) = cellder(4) + addf(1)*ded(3)*de*fct
      cellder(5) = cellder(5) + addf(2)*ded(3)*de*fct
      cellder(6) = cellder(6) + addf(3)*ded(3)*de*fct

      return
      end

      subroutine prcldr(cellder,par)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /athlp/ iatoms, mxnat
      dimension cellder(6),par(3,*)

c      do i=1,6
c         print*,'cellder (',i,')=',cellder(i)
c      end do

      print*,'val ',par(1,iatoms+1),' cellder (1)=',cellder(1)
      print*,'val ',par(2,iatoms+1),' cellder (2)=',cellder(2)
      print*,'val ',par(3,iatoms+1),' cellder (3)=',cellder(3)
      print*,'val ',par(1,iatoms+2),' cellder (4)=',cellder(4)
      print*,'val ',par(2,iatoms+2),' cellder (5)=',cellder(5)
      print*,'val ',par(3,iatoms+2),' cellder (6)=',cellder(6)

      return
      end

      subroutine expfr(in,cxyz,exyz,lex,nexpnd,oje)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /expnd/ addfr(3,125),ie(5),izero
      logical oje(125)
      dimension lex(125),cxyz(3,125),exyz(3,125,*)

      l = 0

      do i=1,125

         if (i.ne.izero) then
            l = l + 1
            do j=1,3
               cxyz(j,l) = exyz(j,i,in) - exyz(j,izero,in)
            end do
            lex(l) = i
            oje(l) = .false.
         endif

      end do
      
      return
      end

      subroutine setfrac(coo,frac,exyz)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /athlp/  iatoms, mxnat
      common /expnd/ addfr(3,125),ie(5),izero
      dimension coo(3,*),frac(3,*),exyz(3,125,*),cxyz(3,125),dfr(3)

      do j=1,3
         dfr(j) = 0.0e0
      end do

      do i=1,iatoms

         do j=1,3
            if (abs(coo(j,i)).lt.0.001e0) coo(j,i) = 0.0e0
         end do

         call crt2fr(coo(1,i),frac(1,i))

         do j=1,3
            if (frac(j,i).lt.dfr(j)) dfr(j) = frac(j,i)
         end do

      end do

c huh ???

      do i=1,5
         do j=1,3
            if ((dfr(j).ge.dble(-1*i)).and.(dfr(j).lt.dble(-1*(i-1)))) 
     &          dfr(j) = dble(-1*i)
         end do
      end do

      do i=1,iatoms

          do j=1,3
             frac(j,i) = frac(j,i) - dfr(j)
          end do

          do k=1,125
             do j=1,3
                exyz(j,k,i) = frac(j,i) + addfr(j,k)
             end do
             call fr2crt(exyz(1,k,i))
          end do

      end do

      return
      end

      subroutine xpfr(in,kn,cxyz,exyz,lex,nexpnd,oje)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /expnd/ addfr(3,125),ie(5),izero
      logical oje(125)
      dimension cxyz(3,125),exyz(3,125,*),lex(125)

      do i=1,nexpnd+1
         oje(i) = .false.
         if (i.eq.izero) oje(i) = .true.
         do j=1,3
            cxyz(j,i) = exyz(j,i,in) - exyz(j,izero,kn)
         end do
         lex(i) = i
      end do
      
      return
      end

      subroutine initxp
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /expnd/ addfr(3,125),ie(5),izero
      data (ie(j),j=1,5)  / -2,-1,0,1,2/

      izero = 0

      l = 0
      do i=1,5
         do j=1,5
            do k=1,5
              l = l + 1
              addfr(1,l) = dble(ie(i))
              addfr(2,l) = dble(ie(j))
              addfr(3,l) = dble(ie(k))
              if (ie(i).eq.0.and.ie(j).eq.0.and.ie(k).eq.0) izero = l
            end do
         end do
      end do

      return
      end

      subroutine var2cl(x,n)
c copy cell variables to cell common
      implicit real (a-h,p-z),integer (i-n),logical (o)
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /cell/ xa,ya,yb,za,zb,zc
      real rxa,rya,ryb,rza,rzb,rzc
      common /cellr/ rxa,rya,ryb,rza,rzb,rzc
      dimension x(*)

      n3 = n*3

      xa = x(n3+1)
      ya = x(n3+2)
      yb = x(n3+3)
      za = x(n3+4)
      zb = x(n3+5)
      zc = x(n3+6)

      abc(1) = sqrt(xa**2+ya**2+za**2)
      abc(2) = sqrt(yb**2+zb**2)
      abc(3) = zc

      do i=1,3
         abc2(i) = abc(i)*0.5e0
      end do

      rxa = real(x(n3+1))
      rya = real(x(n3+2))
      ryb = real(x(n3+3))
      rza = real(x(n3+4))
      rzb = real(x(n3+5))
      rzc = real(x(n3+6))

      return
      end

      subroutine grd2var(cellder,f,n)
c copy cell variables to cell common
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /cell/ xa,ya,yb,za,zb,zc
      dimension cellder(6),f(*)

c      print*,'cellval ',xa,ya,yb,za,zb,zc

      n3 = n*3

      do i=1,6
         f(n3+i) = -cellder(i)
      end do

      return
      end

      subroutine uinner(uin,cellder,coo,f,q,n)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /cell/ xa,ya,yb,za,zb,zc
      dimension cellder(6),pole(3),coo(3,*),f(3,*),q(*)

      pi = 3.141592654e0
      pi43 = 4.0e0*pi/3.0e0
      econv = 332.05382e0
      fin =  pi43*econv/(xa*yb*zc)

      do j=1,3
         pole(j) = 0.0e0
      end do

      do i=1,n
         do j=1,3
            pole(j) = pole(j) + coo(j,i)*q(i)
         end do
      end do

      poll = 0.0e0
      do j=1,3
         poll = poll + pole(j)*pole(j)
      end do

      uin = -0.5e0*fin*poll

      cellder(1) = -uin/xa
      cellder(3) = -uin/yb
      cellder(6) = -uin/zc

      do i=1,n
         do j=1,3
            f(j,i) = f(j,i) - fin*pole(j)*q(i)
         end do
      end do

      return
      end

      subroutine crionp
      parameter (mxion=2000)
      implicit real (a-h,p-z),integer (i-n),logical (o)
      common /h2oer/   numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      logical doit

      nitmp = 0

      do while (nitmp.lt.nion)

         doit = .true.
         irnd =  numwat*random()

         do i=1,nitmp
            if (iw(i).eq.irnd) doit = .false.
         end do

         if (doit) then
           nitmp = nitmp + 1
           iw(nitmp) = irnd
         endif

      end do

      return
      end

      real function random()
      parameter (mplier=16807,modlus=2147483647,mobymp=127773,
     &           momdmp=2836)
      common  /seed/jseed,ifrst,nextn
c     mseed comes from alloc.c
      integer hvlue, lvlue, testv, nextn, mseed
 
      if (ifrst.eq.0) then
        jseed = mseed()
        nextn = jseed
        ifrst = 1
      endif
 
      hvlue = nextn / mobymp
      lvlue = mod(nextn, mobymp)
      testv = mplier*lvlue - momdmp*hvlue

      if (testv.gt.0) then
        nextn = testv
      else
        nextn = testv + modlus
      endif

      random = dble(nextn)/dble(modlus)
 
      return
      end

