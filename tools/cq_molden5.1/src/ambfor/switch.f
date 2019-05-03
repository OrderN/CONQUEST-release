      subroutine inisw
      implicit real (a-h,o-z)
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw

      cutof2 = cutoff*cutoff
      cuton  = cutoff - 1.0e0

      cuton2 = cuton*cuton
      tmp    = cutof2 - cuton2
      conof3 = tmp*tmp*tmp
      isw = 1


      return
      end

      real function swvdw(rij2)
      implicit real (a-h,o-z)
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw

      if (rij2.ge.cutof2) then
          swvdw = 0.0e0
      elseif (rij2.le.cuton2) then
          swvdw = 1.0e0
      else
          rd1 = cutof2 - rij2
          rd2 = cutof2 + 2.0e0*rij2 - 3.0e0*cuton2
          swvdw = (rd1*rd1*rd2)/conof3
      endif

      return
      end

      real function swdvdw(rij2)
c 
c derivative of VDW switch function
c
      implicit real (a-h,o-z)
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw

      if (rij2.ge.cutof2) then
          swdvdw = 0.0e0
      elseif (rij2.le.cuton2) then
          swdvdw = 1.0e0
      else
          rd1 = cutof2 - rij2
c          rd2 = rij2 - cuton2
c          swdvdw = 12.0e0*rd1*rd2/conof3
          rd2 = cutof2 + 2.0e0*rij2 - 3.0e0*cuton2
          rij = sqrt(rij2)
          swdvdw = 4.0e0*rij*rd1*(rd1-rd2)/conof3
      endif

      return
      end

c swchg not used, instead same function as swvdw used for charge cutoff

      subroutine swchg(rij2,s,ds)
      implicit real (a-h,o-z)
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw


      if (rij2.ge.cutof2) then
         s = 0.0e0
         ds = 0.0e0
      else
         frac = rij2/cutof2
         frac2 = frac*frac
         s = 1.0e0 - frac
         s = s*s
         ds = (frac2 - frac)*4.0e0
c        ds = ds/rij
c        we put the rij of de = de / r also in here
c           de = e*ds - e*s/r , de = de / r => 
c                de = e*ds/rij - e*s/rij2
c now       de = e*ds - e*s
         s = s/rij2
         ds = ds/rij2
      endif
          
      return
      end

      subroutine bldlsd(coo,iresid,nlst,lst,istat)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxneib=200)
      parameter (numres=50000)
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw
      common /athlp/  iatoms, mxnat
      common /residu/  qtot,ihsres,
     &                 ires(numres),ibeg(numres),iend(numres)
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      logical doit
      dimension coo(3,*),iresid(*),nlst(*),lst(*),vr(3),loclst(mxneib)

      if (istat.eq.0) then
         usesw = .false.
         print*,'Not enough memory to use cutoff/switch'
         return
      endif

      cutf1 = cutoff + 1.0e0
      cutf2 = cutf1*cutf1

c Build nonbonded interaction list, iatoms*mxneib
c each atom can have a max of mxneib neighbours within cutoff of 8 angs
c now not an atom is stored, but a residue

      do i=1,iatoms
         nlst(i) = 0
         do j=i+1,iatoms

            if (box) then
               do k=1,3
                  vr(k) = coo(k,i) - coo(k,j)
               end do
               call reddis(vr)
               r2 = vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3)
               doit = (r2.lt.cutf2)
            else
               doit = (dist2(coo(1,i),coo(1,j)).lt.cutf2)
            endif

            if (doit) then
               if (nlst(i).lt.mxneib) then
                  ido = 1
                  do k=1,nlst(i)
                     if (iresid(j).eq.lst((i-1)*mxneib+k)) ido = 0
                  end do
                  if (ido.eq.1) then
                     nlst(i) = nlst(i) + 1
                     lst((i-1)*mxneib+nlst(i)) = iresid(j)
                  endif
               else
                  print*,"neighbour list full for atom ",i
                  print*,"increase mxneib "
               endif
            endif
         end do
      end do

      do i=1,ihsres
         if (ires(i).gt.0) then
            nloc = nlst(ibeg(i))
            do j=1,nloc
               loclst(j) = lst((ibeg(i)-1)*mxneib+j)
            end do
            do j=ibeg(i)+1,iend(i)
               do k=1,nlst(j)
                  ir = lst((j-1)*mxneib+k)
                  ido = 1
                  do l=1,nloc
                     if (ir.eq.loclst(l)) ido = 0
                  end do
                  if (ido.eq.1) then
                     if (nloc.lt.mxneib) then
                        nloc = nloc + 1
                        loclst(nloc) = ir
                     else
                        print*,"neighbour list full for atom ",ibeg(i)
                        print*,"increase mxneib "
                     endif
                  endif
               end do
            end do
            do j=ibeg(i),iend(i)
               nlst(j) = nloc
               do l=1,nloc
                  lst((j-1)*mxneib+l) = loclst(l)
               end do
            end do
         endif
      end do

      return
      end

      logical function resinc(iatom,nlst,lst,ires)
      implicit real (a-h,o-z)
      parameter (mxneib=200)
      dimension lst(*)


      resinc = .false.
      do i=1,nlst
         if (lst((iatom-1)*mxneib + i).eq.ires) then
             resinc = .true.
             return
         endif
      end do

      return
      end

      real function dist2(a,b)
c
c     determine distances between neighboring atoms
c
      implicit real (a-h,o-z)
      dimension a(3)
      dimension b(3)

      d1 = a(1)-b(1)
      d2 = a(2)-b(2)
      d3 = a(3)-b(3)

      dist2 = d1*d1 + d2*d2 + d3*d3

      return
      end

      subroutine premul
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcls=50)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      common /prem/ rs(mxcls,mxcls),es(mxcls,mxcls),iprem

      do i=1,mxcls
         do j=1,mxcls
            rst = ambvdwr(i) + ambvdwr(j)
            rs(i,j) = rst*rst*rst*rst*rst*rst
            es(i,j) = sqrt(ambvdwe(i)*ambvdwe(j))
         end do
      end do

      call allerf
      call alltab

      iprem = 1

      return
      end

      subroutine prmulr
      implicit real(a-h,o-z), integer (i-n)
      parameter (mxcls=50)
      real ambvdwr,ambvdwe
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      common /premr/ rs(mxcls,mxcls),es(mxcls,mxcls),iprem

      do i=1,mxcls
         do j=1,mxcls
            rst = real(ambvdwr(i) + ambvdwr(j))
            rs(i,j) = rst*rst*rst*rst*rst*rst
            es(i,j) = sqrt(real(ambvdwe(i)*ambvdwe(j)))
         end do
      end do

      iprem = 1

      call alltab

      return
      end

      subroutine alltab
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxion=2000)
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /h2oer/ numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      common /atab/ ialtab

      if (ialtab.eq.1) return

      call watcnt
      call watlst(niwat)
      call stutab

      ialtab = 1

      return
      end

      subroutine stutad(woo,whh,woh,dewoo,dewhh,dewoh)
      implicit real (a-h,o-z), integer (i-n)
      dimension woo(*),whh(*),woh(*),dewoo(*),dewhh(*),dewoh(*)

      rcut  = 9.0e0
      rcinv = 1.0e0 / rcut
      rc2inv =  rcinv * rcinv

      econv = 332.05382e0
      qo    = -0.8340e0
      qh    = 0.4170e0
      vro   = 1.7683e0
      veo   = 0.1520e0
      vrh   = 0.0001e0
      veh   = 0.0e0

      poo   = econv*qo*qo
      rst   = 2.0e0*vro
      rst2  = rst*rst
      pr6oo = rst2*rst2*rst2
      epsoo = veo

      phh   = econv*qh*qh
      pr6hh = 0.0e0
      epshh = 0.0e0

      poh   = econv*qo*qh
      rst   = vro + vrh
      rst2  = rst*rst
      pr6oh = rst2*rst2*rst2
      epsoh = sqrt(veo*veh)

      nsize = 9000

      do i=1,nsize

         r     = dble(i) / 1000.0e0
         rinv  = 1.0e0 / r
         r2inv = rinv*rinv
         r6inv = r2inv*r2inv*r2inv

         p6    = pr6oo*r6inv
         p12   = p6*p6
         eq    = poo * (rinv - rcinv + (r-rcut)*rc2inv)
         ev    = epsoo*(p12-2.0e0*p6)

         woo(i) = eq + ev

c         if (i.lt.8000) print*,i,' woo ',woo(i),
c     &                  ' dwoo ',woo(i-1)-woo(i)
         de = epsoo * (p12 - p6) * (-12.0e0)
         dewoo(i) =  (de-eq)*r2inv + poo*rinv*rc2inv

         eq    = phh * (rinv - rcinv + (r-rcut)*rc2inv)

         whh(i) = eq

         dewhh(i) =  -eq*r2inv + phh*rinv*rc2inv

         p6    = pr6oh*r6inv
         p12   = p6*p6
         eq    = poh * (rinv - rcinv + (r-rcut)*rc2inv)
         ev    = epsoh*(p12-2.0e0*p6)

         woh(i) = eq + ev

         de = epsoh * (p12 - p6) * (-12.0e0)
         dewoh(i) =  (de-eq)*r2inv + poh*rinv*rc2inv

      end do

      do i=9000,10000
         woo(i) = 0.0e0
         woh(i) = 0.0e0
         whh(i) = 0.0e0
         dewoo(i) = 0.0e0
         dewoh(i) = 0.0e0
         dewhh(i) = 0.0e0
      end do

      return
      end

      subroutine etutad(nsize,woo,whh,woh,dewoo,dewhh,dewoh)
      implicit real (a-h,o-z), integer (i-n)
      dimension woo(*),whh(*),woh(*),dewoo(*),dewhh(*),dewoh(*)

      econv = 332.05382e0
      qo    = -0.8340e0
      qh    = 0.4170e0
      vro   = 1.7683e0
      veo   = 0.1520e0
      vrh   = 0.0001e0
      veh   = 0.0e0

      poo   = econv*qo*qo
      rst   = 2.0e0*vro
      rst2  = rst*rst
      pr6oo = rst2*rst2*rst2
      epsoo = veo

      phh   = econv*qh*qh
      pr6hh = 0.0e0
      epshh = 0.0e0

      poh   = econv*qo*qh
      rst   = vro + vrh
      rst2  = rst*rst
      pr6oh = rst2*rst2*rst2
      epsoh = sqrt(veo*veh)

      do i=1,nsize

         r     = real(i) / 1000.0e0
         rinv  = 1.0e0 / r
         r2inv = rinv*rinv
         r6inv = r2inv*r2inv*r2inv

         p6    = pr6oo*r6inv
         p12   = p6*p6
         eq    = poo * rinv
         ev    = epsoo*(p12-2.0e0*p6)

         woo(i) = eq + ev

         de = epsoo * (p12 - p6) * (-12.0e0)
         dewoo(i) =  (de-eq)*r2inv

         eq    = phh * rinv

         whh(i) = eq

         dewhh(i) =  -eq*r2inv

         p6    = pr6oh*r6inv
         p12   = p6*p6
         eq    = poh * rinv
         ev    = epsoh*(p12-2.0e0*p6)

         woh(i) = eq + ev

         de = epsoh * (p12 - p6) * (-12.0e0)
         dewoh(i) =  (de-eq)*r2inv

      end do

      return
      end

      subroutine allerd(aprerf,aprexp)
      implicit real (a-h,o-z), integer (i-n)
      real aprerf,aprexp,alaf,pi,spi,f
      dimension aprerf(*),aprexp(*)

      alfa   = 0.3e0
      pi     = 4.e0*atan(1.e0)
      spi    = sqrt(pi)
      
      do i=1,1800000

         r     = dble(i) / 1000000.0e0
         aprerf(i) = real(erfc(r))

      end do

      do i=1,3240000

         f     = - float(i) / 1000000.0e0
         aprexp(i) = 2.0*alfa*exp(f)/spi

      end do

      return
      end

c      subroutine allerd(aprerf,aprexp)
c      implicit real (a-h,o-z), integer (i-n)
c      real aprerf,aprexp,alaf,pi,spi,f
c      dimension aprerf(*),aprexp(*)
c
c      alfa   = 0.2e0
c      pi     = 4.e0*atan(1.e0)
c      spi    = sqrt(pi)
c      
c      do i=1,3000000
c
c         r     = dble(i) / 1000000.0e0
c         aprerf(i) = real(erfc(r))
c
c      end do
c
c      do i=1,9000000
c
c         f     = - float(i) / 1000000.0e0
c         aprexp(i) = 2.0*alfa*exp(f)/spi
c
c      end do
c
c      return
c      end

      subroutine apperfd(r,fc,aprerf)
      implicit real (a-h,o-z), integer (i-n)
      real aprerf,r1000,x,w
      dimension aprerf(*)

      thou = 1000000.0e0

      if (r.le.10.0e0) then
          r1000 = real(r)*thou
          ir = int(r1000)
          x = (r1000 - float(ir))
          fc = dble((1.0e0-x)*aprerf(ir) + x*aprerf(ir+1))
      endif

      return
      end

