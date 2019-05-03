      subroutine qvdwd(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
     &                  q,iopt,ityp,fx,fy,fz,
     &                  cx,cy,cz,vdwe,vdwr,potq,potv,
     &                  eps,pre6,iag,ikl,ftmp)
      implicit real(a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxac=3*mxcon)
      parameter (mxad=9*mxcon)
      parameter (mxneib=200)
      parameter (mxion=2000)
cmpi      include 'mpif.h'
cmpi      integer chunk
cmpi      integer status(mpi_status_size)
      real q,coo
      integer taskid,offset
      common /mpih/ nproc,taskid
      common /athlp/  iatoms, mxnat
      integer ambcls
      parameter (mxamb=1590)
      real ambchg
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxcls=50)
      real ambvdwr,ambvdwe
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      common /premr/ rs(mxcls,mxcls),es(mxcls,mxcls),iprem
      common /debug/ ideb
      real gfvdw
      parameter (mxgff=72)
      integer amb2gf
      common /gftyp/  gfvdw(2,mxgff),amb2gf(mxamb)
      logical box, cell, fast, dompi, dovdw
      real abc,abc2,angles
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /h2oer/ numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      integer*2 ityp
      dimension vr(3),ded(3)
      dimension coo(3,*),q(*)
      dimension cx(*),cy(*),cz(*),ityp(*),iconn(mxcon+1,*),
     &          iresid(*),nac(*),nad(*),iac(mxac,*),iad(mxad,*),
     &          fx(*),fy(*),fz(*),potq(*),potv(*),iopt(*),ftmp(*)
      dimension vdwe(*),vdwr(*),iag(*),ikl(*),eps(*),pre6(*),
     &          qbck(13*mxcon),vbck(13*mxcon),ibck(13*mxcon)

cmpi      call mpi_comm_rank(mpi_comm_world, taskid, ierr)


      if (iprem.ne.1) call prmulr

cmpi      dompi = .true.
      econv = 332.05382e0
      v14scq = 1.0e0 / 1.2e0
      v14scv = 0.5e0
      dzero = 0.0e0
      cutvdw = 7.86e0
      dovdw = .true.

      iats  = iatoms
      if (fast.and.box) iats  = natnow

      ec = dzero
      ev = dzero

      do i=1,iatoms
         cx(i) = real(coo(1,i))
         cy(i) = real(coo(2,i))
         cz(i) = real(coo(3,i))
         fx(i) = dzero
         fy(i) = dzero
         fz(i) = dzero

         if (iopt(i).eq.1) then
            potq(i) = econv*real(q(i))
            potv(i) = 1.0e0
         else
            potq(i) = dzero
            potv(i) = dzero
         endif

         i1 = int(ityp(i))
         iag(i) = i1
         if (i1.gt.0) then
            il = ambcls(i1)
            vdwr(i) = real(ambvdwr(il))
            vdwe(i) = real(ambvdwe(il))
            ikl(i) = il
         elseif (i1.le.0) then
            i1 = iabs(i1)
            vdwr(i) = real(gfvdw(1,i1))
            vdwe(i) = real(gfvdw(2,i1))
         endif

      end do

      do i=1,iats

         nbck = 0
         vdwr1 = vdwr(i)
         vdwe1 = vdwe(i)
         i1 = iag(i)
         il = ikl(i)
         qi = real(q(i))

         do j=1,iconn(1,i)
            jj = iconn(1+j,i)
            if (jj.gt.0) then
               nbck = nbck + 1
               ibck(nbck) = jj
               qbck(nbck) = potq(jj)
               vbck(nbck) = potv(jj)
               potq(jj) = dzero
               potv(jj) = dzero
            endif
         end do

         do j=1,nac(i)
            jj = iac(j,i)
            nbck = nbck + 1
            ibck(nbck) = jj
            qbck(nbck) = potq(jj)
            vbck(nbck) = potv(jj)
            potq(jj) = dzero
            potv(jj) = dzero
         end do

         do j=1,nad(i)
            jj = iad(j,i)
            nbck = nbck + 1
            ibck(nbck) = jj
            qbck(nbck) = potq(jj)
            vbck(nbck) = potv(jj)
            potq(jj) = potq(jj)*v14scq
            potv(jj) = potv(jj)*v14scv
         end do

         dovdw = (vdwe1.ne.0.0e0)

         if (dovdw) then

            do k=i+1,iatoms

               vdwe2 = vdwe(k)

               if (vdwe2.ne.dzero) then
   
                  vdwr2 = vdwr(k)
                  i2 = iag(k)
                  kl = ikl(k)

                  if (i1.gt.0.and.i2.gt.0) then
                     pre6(k)  = rs(il,kl)
                     eps(k)   = es(il,kl)
                  else
                     rsum     = vdwr1 + vdwr2
                     rsum3    = rsum*rsum*rsum
                     eps(k)   = sqrt(vdwe1 * vdwe2)
                     pre6(k)  = rsum3*rsum3
                  endif

               else
                  pre6(k)  = 0.0e0
                  eps(k)   = 0.0e0
               endif

            end do
         endif

cmpi         chunk = ((iatoms-i)/nproc)
cmpi         ichkh = chunk/4
cmpi         chunk = ichkh*4
cmpi         if (chunk.lt.8) dompi = .false.
cmpi
cmpi         if (dompi) then
cmpi            if (taskid.eq.0) then
cmpi
cmpi
cmpi               offset = chunk + i + 1
cmpi
cmpi               do it=1, nproc - 1
cmpi                  call mpi_send(offset, 1, mpi_integer, it, 5,
cmpi     &              mpi_comm_world, ierr)
cmpi                  offset = offset + chunk
cmpi               end do
cmpi
cmpi               offset = i + 1
cmpi
cmpi            else
cmpi               call mpi_recv(offset, 1, mpi_integer, 0, 5,
cmpi     &            mpi_comm_world, status, ierr)
cmpi            endif
cmpi
cmpi            iend = offset + chunk - 1
cmpi            if (taskid.eq.nproc-1) iend = iatoms
cmpi
cmpi         else
cmpi
cmpi            if (taskid.eq.0) then
               offset = i + 1
               iend = iatoms
cmpi            else
cmpi               offset = 1
cmpi               iend = 0
cmpi            endif
cmpi
cmpi         endif


         evt = 0.0e0

         do k=offset,iend

               vr(1) = cx(i) - cx(k)
               vr(2) = cy(i) - cy(k)
               vr(3) = cz(i) - cz(k)

               if (box) call rddisr(vr)

               r2 = vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3)

               r = sqrt(r2)

               rinv = 1.0e0 / r
               r2inv = rinv*rinv

               e = qi * potq(k) * rinv

               de = -e * r2inv

               do j=1,3
                  ded(j) = de * vr(j)
               end do

               ec = ec + e

               epsm = eps(k)

               if (dovdw.and.epsm.ne.dzero.and.r.lt.cutvdw) then

                  r6inv = r2inv*r2inv*r2inv
                  p6 = pre6(k) * r6inv

                  epsm = epsm * potv(k)
                  p12  = p6 * p6

                  e    = epsm * (p12 - 2.0e0*p6)
                  de   = epsm * (p12 - p6) * (-12.0e0)

                  de   = de * r2inv

                  do j=1,3
                     ded(j) = ded(j) + de * vr(j)
                  end do

c                  ev   = ev + e
                  evt   = evt + e

               endif

               fx(i) = fx(i) + ded(1)
               fy(i) = fy(i) + ded(2)
               fz(i) = fz(i) + ded(3)

               fx(k) = fx(k) - ded(1)
               fy(k) = fy(k) - ded(2)
               fz(k) = fz(k) - ded(3)

         end do

         ev = ev + evt

         do k=1,nbck
            jj = ibck(k)
            potq(jj) = qbck(k)
            potv(jj) = vbck(k)
         end do


c         call mpi_barrier(mpi_comm_world, ierr)

      end do

cmpi      if (taskid.eq.0) then
cmpi
cmpi         do it=1, nproc - 1
cmpi            itmp = it
cmpi
cmpi            call mpi_recv(ftmp(1), iatoms, mpi_real,
cmpi     &        itmp, 3, mpi_comm_world, status, ierr)
cmpi            do k=1,iatoms
cmpi               fx(k) = fx(k) + ftmp(k)
cmpi            end do
cmpi
cmpi            call mpi_recv(ftmp(1), iatoms, mpi_real,
cmpi     &        itmp, 4, mpi_comm_world, status, ierr)
cmpi            do k=1,iatoms
cmpi               fy(k) = fy(k) + ftmp(k)
cmpi            end do
cmpi
cmpi            call mpi_recv(ftmp(1), iatoms, mpi_real,
cmpi     &        itmp, 5, mpi_comm_world, status, ierr)
cmpi            do k=1,iatoms
cmpi               fz(k) = fz(k) + ftmp(k)
cmpi            end do
cmpi
cmpi         end do
cmpi
cmpi         call mpi_reduce(ec, ectmp, 1, mpi_real, mpi_sum, 0,
cmpi     &                mpi_comm_world, ierr)
cmpi         call mpi_reduce(ev, evtmp, 1, mpi_real, mpi_sum, 0,
cmpi     &                mpi_comm_world, ierr)
cmpi
cmpi         ec = ectmp
cmpi         ev = evtmp
cmpi
cmpi         do it=1, nproc - 1
cmpi            call mpi_send(ectmp, 1, mpi_real, it, 6,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(evtmp, 1, mpi_real, it, 7,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(fx(1), iatoms, mpi_real, it, 8,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(fy(1), iatoms, mpi_real, it, 9,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(fz(1), iatoms, mpi_real, it, 10,
cmpi     &              mpi_comm_world, ierr)
cmpi         end do
cmpi
cmpi      else
cmpi
cmpi         call mpi_send(fx(1), iatoms, mpi_real, 0, 3,
cmpi     &              mpi_comm_world, ierr)
cmpi         call mpi_send(fy(1), iatoms, mpi_real, 0, 4,
cmpi     &              mpi_comm_world, ierr)
cmpi         call mpi_send(fz(1), iatoms, mpi_real, 0, 5,
cmpi     &              mpi_comm_world, ierr)
cmpi
cmpi         call mpi_reduce(ec, ectmp, 1, mpi_real, mpi_sum, 0,
cmpi     &              mpi_comm_world, ierr)
cmpi         call mpi_reduce(ev, evtmp, 1, mpi_real, mpi_sum, 0,
cmpi     &              mpi_comm_world, ierr)
cmpi
cmpi         call mpi_recv(ectmp, 1, mpi_real, 0, 6,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(evtmp, 1, mpi_real, 0, 7,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(fx(1), iatoms, mpi_real, 0, 8,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(fy(1), iatoms, mpi_real, 0, 9,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(fz(1), iatoms, mpi_real, 0, 10,
cmpi     &              mpi_comm_world, status, ierr)
cmpi
cmpi         ec = ectmp
cmpi         ev = evtmp
cmpi
cmpi      endif

      return
      end

      subroutine qvdwdd(ec,ev,ftmp,nac,iac,nad,iad,iresid,coo,iconn,
     &                  q,iopt,ityp,
     &                  vdwe,vdwr,potq,potv,
     &                  eps,pre6,iag,ikl,fftmp)
      implicit real(a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxac=3*mxcon)
      parameter (mxad=9*mxcon)
      parameter (mxneib=200)
      parameter (mxion=2000)
cmpi      include 'mpif.h'
cmpi      integer chunk
cmpi      integer status(mpi_status_size)
      integer taskid,offset
      common /mpih/ nproc,taskid
      common /athlp/  iatoms, mxnat
      integer ambcls
      parameter (mxamb=1590)
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxcls=50)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      common /prem/ rs(mxcls,mxcls),es(mxcls,mxcls),iprem
      parameter (mxgff=72)
      integer amb2gf
      common /gftyp/  gfvdw(2,mxgff),amb2gf(mxamb)
      logical box, cell, fast, dompi, dovdw
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /h2oer/ numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      integer*2 ityp
      dimension vr(3),ded(3)
      dimension coo(3,*),q(*)
      dimension ityp(*),iconn(mxcon+1,*),
     &          iresid(*),nac(*),nad(*),iac(mxac,*),iad(mxad,*),
     &          potq(*),potv(*),iopt(*),ftmp(3,*),fftmp(3,*)
      dimension vdwe(*),vdwr(*),iag(*),ikl(*),eps(*),pre6(*),
     &          qbck(13*mxcon),vbck(13*mxcon),ibck(13*mxcon)

cmpi      call mpi_comm_rank(mpi_comm_world, taskid, ierr)


      if (iprem.ne.1) call premul

      dompi = .true.
      econv = 332.05382e0
      v14scq = 1.0e0 / 1.2e0
      v14scv = 0.5e0
      dzero = 0.0e0
      cutvdw = 7.86e0
      dovdw = .true.

      ec = 0.0e0
      ev = 0.0e0

      iat3 = iatoms*3

      iats  = iatoms
      if (fast.and.box) iats  = natnow

      do i=1,iatoms

         do j=1,3
            ftmp(j,i) = dzero
         end do

         if (iopt(i).eq.1) then
            potq(i) = econv*q(i)
            potv(i) = 1.0e0
         else
            potq(i) = dzero
            potv(i) = dzero
         endif

         i1 = int(ityp(i))
         iag(i) = i1
         if (i1.gt.0) then
            il = ambcls(i1)
            vdwr(i) = ambvdwr(il)
            vdwe(i) = ambvdwe(il)
            ikl(i) = il
         elseif (i1.le.0) then
            i1 = iabs(i1)
            vdwr(i) = gfvdw(1,i1)
            vdwe(i) = gfvdw(2,i1)
         endif

      end do

      do i=1,iats

         nbck = 0
         vdwr1 = vdwr(i)
         vdwe1 = vdwe(i)
         i1 = iag(i)
         il = ikl(i)
         qi = q(i)

         do j=1,iconn(1,i)
            jj = iconn(1+j,i)
            if (jj.gt.0) then
               nbck = nbck + 1
               ibck(nbck) = jj
               qbck(nbck) = potq(jj)
               vbck(nbck) = potv(jj)
               potq(jj) = dzero
               potv(jj) = dzero
            endif
         end do

         do j=1,nac(i)
            jj = iac(j,i)
            nbck = nbck + 1
            ibck(nbck) = jj
            qbck(nbck) = potq(jj)
            vbck(nbck) = potv(jj)
            potq(jj) = dzero
            potv(jj) = dzero
         end do

         do j=1,nad(i)
            jj = iad(j,i)
            nbck = nbck + 1
            ibck(nbck) = jj
            qbck(nbck) = potq(jj)
            vbck(nbck) = potv(jj)
            potq(jj) = potq(jj)*v14scq
            potv(jj) = potv(jj)*v14scv
         end do

         dovdw = (vdwe1.ne.0.0e0)

         if (dovdw) then

            do k=i+1,iatoms

               vdwe2 = vdwe(k)

               if (vdwe2.ne.dzero) then
   
                  vdwr2 = vdwr(k)
                  i2 = iag(k)
                  kl = ikl(k)

                  if (i1.gt.0.and.i2.gt.0) then
                     pre6(k)  = rs(il,kl)
                     eps(k)   = es(il,kl)
                  else
                     rsum     = vdwr1 + vdwr2
                     rsum3    = rsum*rsum*rsum
                     eps(k)   = sqrt(vdwe1 * vdwe2)
                     pre6(k)  = rsum3*rsum3
                  endif

               else

                  pre6(k)  = dzero
                  eps(k)   = dzero

               endif
            end do
         endif

cmpi         chunk = ((iatoms-i)/nproc)
cmpi         ichkh = chunk/4
cmpi         chunk = ichkh*4
cmpi         if (chunk.lt.8) dompi = .false.
cmpi
cmpi         if (dompi) then
cmpi            if (taskid.eq.0) then
cmpi
cmpi
cmpi               offset = chunk + i + 1
cmpi
cmpi               do it=1, nproc - 1
cmpi                  call mpi_send(offset, 1, mpi_integer, it, 5,
cmpi     &              mpi_comm_world, ierr)
cmpi                  offset = offset + chunk
cmpi               end do
cmpi
cmpi               offset = i + 1
cmpi
cmpi            else
cmpi               call mpi_recv(offset, 1, mpi_integer, 0, 5,
cmpi     &            mpi_comm_world, status, ierr)
cmpi            endif
cmpi
cmpi            iend = offset + chunk - 1
cmpi            if (taskid.eq.nproc-1) iend = iatoms
cmpi
cmpi         else
cmpi
cmpi            if (taskid.eq.0) then
               offset = i + 1
               iend = iatoms
cmpi            else
cmpi               offset = 1
cmpi               iend = 0
cmpi            endif
cmpi
cmpi         endif


         do k=offset,iend

               do j=1,3
                  vr(j) = coo(j,i) - coo(j,k)
               end do

               if (box) call reddis(vr)

               r2 = vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3)

               r = sqrt(r2)

               rinv = 1.0e0 / r
               r2inv = rinv*rinv

               e = qi * potq(k) * rinv

               de = -e * r2inv

               do j=1,3
                  ded(j) = de * vr(j)
               end do

               ec = ec + e

               epsm = eps(k)

               if (dovdw.and.epsm.ne.dzero.and.r.lt.cutvdw) then

                  r6inv = r2inv*r2inv*r2inv
                  p6 = pre6(k) * r6inv

                  epsm = epsm * potv(k)
                  p12  = p6 * p6

                  e    = epsm * (p12 - 2.0e0*p6)
                  de   = epsm * (p12 - p6) * (-12.0e0)

                  de   = de * r2inv

                  do j=1,3
                     ded(j) = ded(j) + de * vr(j)
                  end do

                  ev   = ev + e

               endif

               do j=1,3
                  ftmp(j,i) = ftmp(j,i) + ded(j)
                  ftmp(j,k) = ftmp(j,k) - ded(j)
               end do

         end do

         do k=1,nbck
            jj = ibck(k)
            potq(jj) = qbck(k)
            potv(jj) = vbck(k)
         end do

      end do

cmpi      if (taskid.eq.0) then
cmpi
cmpi         do it=1, nproc - 1
cmpi            itmp = it
cmpi
cmpi            call mpi_recv(fftmp(1,1), iat3, mpi_real,
cmpi     &        itmp, 3, mpi_comm_world, status, ierr)
cmpi            do k=1,iatoms
cmpi               do j=1,3
cmpi                  ftmp(j,k) = ftmp(j,k) + fftmp(j,k)
cmpi               end do
cmpi            end do
cmpi
cmpi         end do
cmpi
cmpi         call mpi_reduce(ec, ectmp, 1, mpi_real, mpi_sum,
cmpi     &                0, mpi_comm_world, ierr)
cmpi         call mpi_reduce(ev, evtmp, 1, mpi_real, mpi_sum,
cmpi     &                0, mpi_comm_world, ierr)
cmpi
cmpi         ec = ectmp
cmpi         ev = evtmp
cmpi
cmpi         do it=1, nproc - 1
cmpi            call mpi_send(ectmp, 1, mpi_real, it, 6,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(evtmp, 1, mpi_real, it, 7,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(ftmp(1,1), iat3, 
cmpi     &              mpi_real, it, 8,mpi_comm_world, ierr)
cmpi         end do
cmpi
cmpi      else
cmpi
cmpi         call mpi_send(ftmp(1,1), iat3, mpi_real, 0, 3,
cmpi     &              mpi_comm_world, ierr)
cmpi
cmpi         call mpi_reduce(ec, ectmp, 1, mpi_real, mpi_sum,
cmpi     &              0, mpi_comm_world, ierr)
cmpi         call mpi_reduce(ev, evtmp, 1, mpi_real, mpi_sum, 
cmpi     &              0, mpi_comm_world, ierr)
cmpi
cmpi         call mpi_recv(ectmp, 1, mpi_real, 0, 6,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(evtmp, 1, mpi_real, 0, 7,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(ftmp(1,1), iat3, mpi_real, 
cmpi     &              0, 8, mpi_comm_world, status, ierr)
cmpi
cmpi         ec = ectmp
cmpi         ev = evtmp
cmpi
cmpi      endif

      return
      end

      subroutine qvd(ec,ev,nac,iac,nad,iad,iresid,coo,iconn,
     &                  q,forces,cellder,iopt,ityp,
     &                  vdwe,vdwr,potq,potv,
     &                  eps,pre6,iag,ikl,ftmp,fftmp,frac,exyz)
      implicit real(a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxac=3*mxcon)
      parameter (mxad=9*mxcon)
      parameter (mxneib=200)
cmpi      include 'mpif.h'
cmpi      integer chunk
cmpi      integer status(mpi_status_size)
      integer taskid,offset
      common /mpih/ nproc,taskid
      common /athlp/  iatoms, mxnat
      integer ambcls
      parameter (mxamb=1590)
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxcls=50)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      common /prem/ rs(mxcls,mxcls),es(mxcls,mxcls),iprem
      parameter (mxgff=72)
      integer amb2gf
      common /gftyp/  gfvdw(2,mxgff),amb2gf(mxamb)
      logical box, cell, fast, dompi, lcon
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /expnd/ addfr(3,125),ie(5),izero
      logical oje(125)
      integer*2 ityp
      dimension vr(3),ded(3),cellder(*),cxyz(3,125),lex(125)
      dimension coo(3,*),forces(3,*),q(*),frac(3,*),exyz(3,125,*)
      dimension ityp(*),iconn(mxcon+1,*),
     &          iresid(*),nac(*),nad(*),iac(mxac,*),iad(mxad,*),
     &          potq(*),potv(*),iopt(*),ftmp(3,*),fftmp(3,*),
     &          cder(6)
      dimension vdwe(*),vdwr(*),iag(*),ikl(*),eps(*),pre6(*),
     &          qbck(13*mxcon),vbck(13*mxcon),ibck(13*mxcon)

c difference with other routines: second loop runs from 1,iatoms
c instead of i,iatoms
c
c Secondly, routine could be faster if we go back to i,iatoms
c           then we have to correctly incoporate the second
c           expansion
c
cmpi      call mpi_comm_rank(mpi_comm_world, taskid, ierr)


      nexpnd = 124

cmpi      dompi = .true.
cmpi      chunk = iatoms/nproc
cmpi      if (chunk.lt.5) dompi = .false.

      if (iprem.ne.1) call premul

      call setfrac(coo,frac,exyz)

      do i=1,6
         cellder(i) = 0.0e0
      end do



      econv = 332.05382e0
      v14scq = 1.0e0 / 1.2e0
      v14scv = 0.5e0

      ec = 0.0e0
      ev = 0.0e0
      iat3 = iatoms*3

      do i=1,iatoms

         do j=1,3
            ftmp(j,i) = 0.0e0
         end do

         if (iopt(i).eq.1) then
            potq(i) = econv*q(i)
            potv(i) = 1.0e0
         else
            potq(i) = 0.0e0
            potv(i) = 0.0e0
         endif

         i1 = int(ityp(i))
         iag(i) = i1
         if (i1.gt.0) then
            il = ambcls(i1)
            vdwr(i) = ambvdwr(il)
            vdwe(i) = ambvdwe(il)
            ikl(i) = il
         elseif (i1.le.0) then
            i1 = iabs(i1)
            vdwr(i) = gfvdw(1,i1)
            vdwe(i) = gfvdw(2,i1)
         endif

      end do

      do i=1,iatoms

         nbck = 0
         vdwr1 = vdwr(i)
         vdwe1 = vdwe(i)
         i1 = iag(i)
         il = ikl(i)
         qi = q(i)

         do j=1,iconn(1,i)
            jj = iconn(1+j,i)
            if (jj.gt.0) then
               nbck = nbck + 1
               ibck(nbck) = jj
               qbck(nbck) = potq(jj)
               vbck(nbck) = potv(jj)
               potq(jj) = 0.0e0
               potv(jj) = 0.0e0
            endif
         end do

         do j=1,nac(i)
            jj = iac(j,i)
            nbck = nbck + 1
            ibck(nbck) = jj
            qbck(nbck) = potq(jj)
            vbck(nbck) = potv(jj)
            potq(jj) = 0.0e0
            potv(jj) = 0.0e0
         end do

         do j=1,nad(i)
            jj = iad(j,i)
            nbck = nbck + 1
            ibck(nbck) = jj
            qbck(nbck) = potq(jj)
            vbck(nbck) = potv(jj)
            potq(jj) = potq(jj)*v14scq
            potv(jj) = potv(jj)*v14scv
         end do

         do k=1,iatoms

            if (vdwe1.ne.0.0e0) then
               vdwe2 = vdwe(k)

               if (vdwe2.ne.0.0e0) then
   
                  vdwr2 = vdwr(k)
                  i2 = iag(k)
                  kl = ikl(k)

                  if (i1.gt.0.and.i2.gt.0) then
                     pre6(k)  = rs(il,kl)
                     eps(k)   = es(il,kl)
                  else
                     rsum     = vdwr1 + vdwr2
                     rsum3    = rsum*rsum*rsum
                     eps(k)   = sqrt(vdwe1 * vdwe2)
                     pre6(k)  = rsum3*rsum3
                  endif

               else

                  pre6(k)  = 0.0e0
                  eps(k)   = 0.0e0

               endif
            else
               pre6(k)  = 0.0e0
               eps(k)   = 0.0e0
            endif
         end do

cmpi         if (dompi) then
cmpi            ibeg = chunk*taskid + 1
cmpi            iend = chunk*(taskid+1)
cmpi            if (iend.gt.iatoms) iend = iatoms
cmpi         else
            ibeg = 1
            iend = iatoms
cmpi         endif

         do k=ibeg,iend

               epsm = eps(k)*potv(k)
               ptq  = qi*potq(k)
               pr6  = pre6(k)

               ist = 1
               if (i.eq.k) then
                  call expfr(i,cxyz,exyz,lex,nexpnd,oje)
                  lst = nexpnd
               else
                  call xpfr(i,k,cxyz,exyz,lex,nexpnd,oje)
                  lst = 1 + nexpnd
               endif

               do l=ist,lst

                  fct = 0.5e0

                  r2 = cxyz(1,l)*cxyz(1,l) + 
     &                 cxyz(2,l)*cxyz(2,l) + 
     &                 cxyz(3,l)*cxyz(3,l)

                  r = sqrt(r2)

                  if (oje(l).and.i.ne.k) then
                     e = ptq / r
                  else
                     e = econv*q(i)*q(k) / r
                  endif

                  de = -e / r2

                  call adcldr(cellder,cxyz(1,l),de,fct,lex(l))

                  do j=1,3
                     ded(j) = de * fct * cxyz(j,l)
                  end do

                  ec = ec + fct*e

                  p6 = pr6 / (r2*r2*r2)

                  p12  = p6 * p6

                
                  if (oje(l).and.i.ne.k) then
                     epsm = eps(k)*potv(k)
                  else
                     epsm = eps(k)
                  endif

                  e    = epsm * (p12 - 2.0e0*p6)
                  de   = epsm * (p12 - p6) * (-12.0e0)

                  de   = de / r2

                  do j=1,3
                     ded(j) = ded(j) + de * fct * cxyz(j,l)
                  end do


                  ev = ev + fct*e

                  do j=1,3
                     ftmp(j,i) = ftmp(j,i) + ded(j)
                     ftmp(j,k) = ftmp(j,k) - ded(j)
                  end do

                  call adcldr(cellder,cxyz(1,l),de,fct,lex(l))

               end do

         end do

         do k=1,nbck
            jj = ibck(k)
            potq(jj) = qbck(k)
            potv(jj) = vbck(k)
         end do

      end do

cmpi      if (taskid.eq.0) then
cmpi
cmpi         do it=1, nproc - 1
cmpi            itmp = it
cmpi
cmpi            call mpi_recv(fftmp, iat3, mpi_real,
cmpi     &        itmp, 3, mpi_comm_world, status, ierr)
cmpi
cmpi            do i=1,iatoms
cmpi               do j=1,3
cmpi                  ftmp(j,i) = ftmp(j,i) + fftmp(j,i)
cmpi               end do
cmpi            end do
cmpi
cmpi            call mpi_recv(cder, 6, mpi_real,
cmpi     &        itmp, 5, mpi_comm_world, status, ierr)
cmpi
cmpi            do j=1,6
cmpi               cellder(j) = cellder(j) + cder(j)
cmpi            end do
cmpi
cmpi         end do
cmpi
cmpi         call mpi_reduce(ec, ectmp, 1, mpi_real, mpi_sum,
cmpi     &                0, mpi_comm_world, ierr)
cmpi         call mpi_reduce(ev, evtmp, 1, mpi_real, mpi_sum,
cmpi     &                0, mpi_comm_world, ierr)
cmpi
cmpi         ec = ectmp
cmpi         ev = evtmp
cmpi
cmpi         do it=1, nproc - 1
cmpi            call mpi_send(ectmp, 1, mpi_real, it, 6,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(evtmp, 1, mpi_real, it, 7,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(ftmp, iat3, 
cmpi     &              mpi_real, it, 8,mpi_comm_world, ierr)
cmpi            call mpi_send(cellder, 6, 
cmpi     &              mpi_real, it, 10,mpi_comm_world, ierr)
cmpi         end do
cmpi
cmpi      else
cmpi
cmpi         call mpi_send(ftmp, iat3, mpi_real, 0, 3,
cmpi     &              mpi_comm_world, ierr)
cmpi
cmpi         call mpi_send(cellder, 6, mpi_real, 0, 5,
cmpi     &              mpi_comm_world, ierr)
cmpi
cmpi         call mpi_reduce(ec, ectmp, 1, mpi_real, mpi_sum,
cmpi     &              0, mpi_comm_world, ierr)
cmpi         call mpi_reduce(ev, evtmp, 1, mpi_real, mpi_sum, 
cmpi     &              0, mpi_comm_world, ierr)
cmpi
cmpi         call mpi_recv(ectmp, 1, mpi_real, 0, 6,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(evtmp, 1, mpi_real, 0, 7,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(ftmp, iat3, mpi_real, 
cmpi     &              0, 8, mpi_comm_world, status, ierr)
cmpi         call mpi_recv(cellder, 6, mpi_real, 
cmpi     &              0, 10, mpi_comm_world, status, ierr)
cmpi
cmpi         ec = ectmp
cmpi         ev = evtmp
cmpi
cmpi      endif

      do i=1,iatoms
         do j=1,3
            forces(j,i) = forces(j,i) + ftmp(j,i)
         end do
      end do

c      call prcldr(cellder,coo)

      return
      end

      subroutine qenvdd(ec,ev,nac,iac,nad,iad,iresid,
     &                  coo,iconn,q,forces,iopt,
     &                  nlst,lst,ityp,
     &                  idoq,idov,vdwe,vdwr,potq,potv,
     &                  eps,pre6,iag,ikl,ftmp,fftmp)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxac=3*mxcon)
      parameter (mxad=9*mxcon)
      parameter (mxneib=200)
      parameter (mxion=2000)
cmpi      include 'mpif.h'
cmpi      integer chunk
cmpi      integer status(mpi_status_size)
      integer taskid,offset
      common /mpih/ nproc,taskid
      common /athlp/  iatoms, mxnat
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw
      integer ambcls
      parameter (mxamb=1590)
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxcls=50)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      common /prem/ rs(mxcls,mxcls),es(mxcls,mxcls),iprem
      parameter (mxgff=72)
      integer amb2gf
      common /gftyp/  gfvdw(2,mxgff),amb2gf(mxamb)
      logical box,cell,fast,resinc,doit,dompi
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /h2oer/ numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      integer*2 ityp
      dimension vr(3),ded(3)
      dimension coo(3,*),q(*),ityp(*),iconn(mxcon+1,*),iresid(*),
     &          nac(*),nad(*),iac(mxac,*),iad(mxad,*),
     &          forces(3,*),iopt(*),nlst(*),lst(*)
      dimension vdwe(*),vdwr(*),potq(*),potv(*),eps(*),pre6(*),
     &          iag(*),ikl(*),idoq(*),idov(*),ftmp(3,*),fftmp(3,*),
     &          qbck(13*mxcon),vbck(13*mxcon),ibck(13*mxcon)

cmpi      call mpi_comm_rank(mpi_comm_world, taskid, ierr)

      if (iprem.ne.1) call premul

      dompi = .true.
      econv = 332.05382e0
      v14scq = 1.0e0 / 1.2e0
      v14scv = 0.5e0

      ec = 0.0e0
      ev = 0.0e0
      iat3 = iatoms*3
      iats  = iatoms
      if (fast.and.box) iats  = natnow

      do i=1,iatoms

         do j=1,3
            ftmp(j,i) = 0.0e0
         end do

         if (iopt(i).eq.1) then
            potq(i) = econv*q(i)
            potv(i) = 1.0e0
         else
            potq(i) = 0.0e0
            potv(i) = 0.0e0
         endif

         i1 = int(ityp(i))
         iag(i) = i1
         if (i1.gt.0) then
            il = ambcls(i1)
            vdwr(i) = ambvdwr(il)
            vdwe(i) = ambvdwe(il)
            ikl(i) = il
         elseif (i1.le.0) then
            i1 = iabs(i1)
            vdwr(i) = gfvdw(1,i1)
            vdwe(i) = gfvdw(2,i1)
         endif

      end do

      do i=1,iats

         nbck = 0
         vdwr1 = vdwr(i)
         vdwe1 = vdwe(i)
         i1 = iag(i)
         il = ikl(i)
         qi = q(i)

         do j=1,iconn(1,i)
            jj = iconn(1+j,i)
            if (jj.gt.0) then
               nbck = nbck + 1
               ibck(nbck) = jj
               qbck(nbck) = potq(jj)
               vbck(nbck) = potv(jj)
               potq(jj) = 0.0e0
               potv(jj) = 0.0e0
            endif
         end do

         do j=1,nac(i)
            jj = iac(j,i)
            nbck = nbck + 1
            ibck(nbck) = jj
            qbck(nbck) = potq(jj)
            vbck(nbck) = potv(jj)
            potq(jj) = 0.0e0
            potv(jj) = 0.0e0
         end do

         do j=1,nad(i)
            jj = iad(j,i)
            nbck = nbck + 1
            ibck(nbck) = jj
            qbck(nbck) = potq(jj)
            vbck(nbck) = potv(jj)
            potq(jj) = potq(jj)*v14scq
            potv(jj) = potv(jj)*v14scv
         end do

         do k=i+1,iatoms

            idoq(k) = 0
            idov(k) = 0

            if (iopt(i).eq.1.or.iopt(k).eq.1) then
               idoq(k) = 1
               idov(k) = 1
            endif

            if (vdwe1.ne.0.0e0) then
               vdwe2 = vdwe(k)

               if (vdwe2.ne.0.0e0) then
   
                  vdwr2 = vdwr(k)
                  i2 = iag(k)
                  kl = ikl(k)

                  if (i1.gt.0.and.i2.gt.0) then
                     pre6(k)  = rs(il,kl)
                     eps(k)   = es(il,kl)
                  else
                     rsum     = vdwr1 + vdwr2
                     rsum3    = rsum*rsum*rsum
                     eps(k)   = sqrt(vdwe1 * vdwe2)
                     pre6(k)  = rsum3*rsum3
                  endif

               else

                  idov(k) = 0

               endif
            else
               idov(k) = 0
            endif

            if (.not.resinc(i,nlst(i),lst,iresid(k))) then
               idoq(k) = 0
               idov(k) = 0
            endif

         end do

cmpi         chunk = ((iatoms-i)/nproc)
cmpi         ichkh = chunk/4
cmpi         chunk = ichkh*4
cmpi         if (chunk.lt.8) dompi = .false.
cmpi
cmpi         if (dompi) then
cmpi            if (taskid.eq.0) then
cmpi
cmpi
cmpi               offset = chunk + i + 1
cmpi
cmpi               do it=1, nproc - 1
cmpi                  call mpi_send(offset, 1, mpi_integer, it, 5,
cmpi     &              mpi_comm_world, ierr)
cmpi                  offset = offset + chunk
cmpi               end do
cmpi
cmpi               offset = i + 1
cmpi
cmpi            else
cmpi               call mpi_recv(offset, 1, mpi_integer, 0, 5,
cmpi     &            mpi_comm_world, status, ierr)
cmpi            endif
cmpi
cmpi            iend = offset + chunk - 1
cmpi            if (taskid.eq.nproc-1) iend = iatoms
cmpi
cmpi         else
cmpi
cmpi            if (taskid.eq.0) then
               offset = i + 1
               iend = iatoms
cmpi            else
cmpi               offset = 1
cmpi               iend = 0
cmpi            endif
cmpi
cmpi         endif

         do k=offset,iend

           do j=1,3
              ded(j) = 0.0e0
           end do

           if (idoq(k).eq.1.or.idov(k).eq.1) then
               do j=1,3
                  vr(j) = coo(j,i) - coo(j,k)
               end do
               r2 = vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3)
               if (box) call reddis(vr)
               r = sqrt(r2)
           endif

           if (idoq(k).eq.1) then


               if (r2.le.cutof2) then


                  e = qi * potq(k) / r

c since we do: de = de / r, and de*vr(j), we might aswell
c put 1/r in the dw

                  call swchg(r2,s,ds)

                  de = e*(ds - s)
                  e = e * s * r2

                  do j=1,3
                     ded(j) = de * vr(j)
                  end do

                  ec = ec + e

              endif

           end if

           if (idov(k).eq.1) then


               p6 = pre6(k) / (r2*r2*r2)
               epsm = eps(k)
               epsm = epsm * potv(k)


               p12  = p6 * p6

               e    = epsm * (p12 - 2.0e0*p6)
               de   = epsm * (p12 - p6) * (-12.0e0/r)

               if (r2.gt.cuton2) then
                  sw = swvdw(r2)
                  dsw = swdvdw(r2)
                  desw = de*sw
                  edsw = e*dsw
                  de = e*dsw + de*sw
                  e = e * sw
               endif

               de   = de / r

               do j=1,3
                  ded(j) = ded(j) + de * vr(j)
               end do

               ev   = ev + e
   
           end if

           do j=1,3
              ftmp(j,i) = ftmp(j,i) + ded(j)
              ftmp(j,k) = ftmp(j,k) - ded(j)
           end do

         end do

         do k=1,nbck
            jj = ibck(k)
            potq(jj) = qbck(k)
            potv(jj) = vbck(k)
         end do

      end do

cmpi      if (taskid.eq.0) then
cmpi
cmpi         do it=1, nproc - 1
cmpi            itmp = it
cmpi
cmpi            call mpi_recv(fftmp(1,1), iat3, mpi_real,
cmpi     &        itmp, 3, mpi_comm_world, status, ierr)
cmpi            do k=1,iatoms
cmpi               do j=1,3
cmpi                  ftmp(j,k) = ftmp(j,k) + fftmp(j,k)
cmpi               end do
cmpi            end do
cmpi
cmpi         end do
cmpi
cmpi         call mpi_reduce(ec, ectmp, 1, mpi_real, mpi_sum,
cmpi     &                0, mpi_comm_world, ierr)
cmpi         call mpi_reduce(ev, evtmp, 1, mpi_real, mpi_sum,
cmpi     &                0, mpi_comm_world, ierr)
cmpi
cmpi         ec = ectmp
cmpi         ev = evtmp
cmpi
cmpi         do it=1, nproc - 1
cmpi            call mpi_send(ectmp, 1, mpi_real, it, 6,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(evtmp, 1, mpi_real, it, 7,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(ftmp(1,1), iat3, 
cmpi     &              mpi_real, it, 8,mpi_comm_world, ierr)
cmpi         end do
cmpi
cmpi      else
cmpi
cmpi         call mpi_send(ftmp(1,1), iat3, mpi_real, 0, 3,
cmpi     &              mpi_comm_world, ierr)
cmpi
cmpi         call mpi_reduce(ec, ectmp, 1, mpi_real, mpi_sum,
cmpi     &              0, mpi_comm_world, ierr)
cmpi         call mpi_reduce(ev, evtmp, 1, mpi_real, mpi_sum, 
cmpi     &              0, mpi_comm_world, ierr)
cmpi
cmpi         call mpi_recv(ectmp, 1, mpi_real, 0, 6,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(evtmp, 1, mpi_real, 0, 7,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(ftmp(1,1), iat3, mpi_real, 
cmpi     &              0, 8, mpi_comm_world, status, ierr)
cmpi
cmpi         ec = ectmp
cmpi         ev = evtmp
cmpi
cmpi      endif

      do i=1,iatoms
         do j=1,3
            forces(j,i) = forces(j,i) + ftmp(j,i)
         end do
      end do

      return
      end

      subroutine qvnod(ec,ev,nac,iac,nad,iad,iresid,
     &                  coo,iconn,q,forces,iopt,ityp,iwtpr,
     &                  vdwe,vdwr,potq,potv,
     &                  eps,pre6,iag,ikl,ftmp,fftmp)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxac=3*mxcon)
      parameter (mxad=9*mxcon)
      parameter (mxneib=200)
      parameter (mxion=2000)
cmpi      include 'mpif.h'
cmpi      integer chunk
cmpi      integer status(mpi_status_size)
      integer taskid,offset
      common /mpih/ nproc,taskid
      common /athlp/  iatoms, mxnat
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw
      integer ambcls
      parameter (mxamb=1590)
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxcls=50)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      common /prem/ rs(mxcls,mxcls),es(mxcls,mxcls),iprem
      parameter (mxgff=72)
      integer amb2gf
      common /gftyp/  gfvdw(2,mxgff),amb2gf(mxamb)
      logical box,cell,fast,resinc,doit,dompi
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /h2oer/ numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      logical opt
      integer*2 ityp
      dimension vr(3),ded(3)
      dimension coo(3,*),q(*),ityp(*),iconn(mxcon+1,*),iresid(*),
     &          nac(*),nad(*),iac(mxac,*),iad(mxad,*),
     &          forces(3,*),iopt(*),iwtpr(*)
      dimension vdwe(*),vdwr(*),potq(*),potv(*),eps(*),pre6(*),
     &          iag(*),ikl(*),ftmp(3,*),fftmp(3,*),
     &          qbck(13*mxcon),vbck(13*mxcon),ibck(13*mxcon)

cmpi      call mpi_comm_rank(mpi_comm_world, taskid, ierr)

      if (iprem.ne.1) call premul

      dompi = .true.
      econv = 332.05382e0
      v14scq = 1.0e0 / 1.2e0
      v14scv = 0.5e0
      pi     = 4.d0*datan(1.d0)
      spi    = sqrt(pi)
      alfa   = 0.3e0
      alfa2  = alfa*alfa
      rc     = 9.0e0
      rc2    = rc*rc
      erfrc  = erfc(alfa*rc) / rc

c     max afstand rc = 9.0, alfa = 0.2:  max alfa*r = 1.8
c     approx array 1,8 miljoen ( 1.8 * 1.000000).
c     max -alfa**2*rc**2 = -3.24, approx exp array  3240000

c     max afstand rc = 9.0, alfa = 0.3:  max alfa*r = 2.7
c     approx array 2,7 miljoen ( 2.7 * 1.000000).
c     max -alfa**2*rc**2 = -7.29, approx exp array  7290000


c     max afstand rc = 15.0, alfa = 0.2:  max alfa*r = 3.0
c     approx array 3 miljoen ( 3 * 1.000000).
c
c     max -alfa**2*rc**2 = -9.0, approx exp array 9 miljoen
c     

      erfrc2 = (erfrc  + 
     &         2.0e0*alfa*exp(-1.0e0*alfa2*rc2)/spi)/rc

      ec = 0.0e0
      ev = 0.0e0
      iat3 = iatoms*3

      iats  = iatoms
      if (fast.and.box) iats  = natnow


      do i=1,iatoms

         do j=1,3
            ftmp(j,i) = 0.0e0
         end do

         if (iopt(i).eq.1) then
            potq(i) = econv*q(i)
            potv(i) = 1.0e0
         else
            potq(i) = 0.0e0
            potv(i) = 0.0e0
         endif

         i1 = int(ityp(i))
         iag(i) = i1
         if (i1.gt.0) then
            il = ambcls(i1)
            vdwr(i) = ambvdwr(il)
            vdwe(i) = ambvdwe(il)
            ikl(i) = il
         elseif (i1.le.0) then
            i1 = iabs(i1)
            vdwr(i) = gfvdw(1,i1)
            vdwe(i) = gfvdw(2,i1)
         endif

      end do

      do i=1,iats

         nbck = 0
         vdwr1 = vdwr(i)
         vdwe1 = vdwe(i)
         i1 = iag(i)
         il = ikl(i)
         qi = q(i)

         do j=1,iconn(1,i)
            jj = iconn(1+j,i)
            if (jj.gt.0) then
               nbck = nbck + 1
               ibck(nbck) = jj
               qbck(nbck) = potq(jj)
               vbck(nbck) = potv(jj)
               potq(jj) = 0.0e0
               potv(jj) = 0.0e0
            endif
         end do

         do j=1,nac(i)
            jj = iac(j,i)
            nbck = nbck + 1
            ibck(nbck) = jj
            qbck(nbck) = potq(jj)
            vbck(nbck) = potv(jj)
            potq(jj) = 0.0e0
            potv(jj) = 0.0e0
         end do

         do j=1,nad(i)
            jj = iad(j,i)
            nbck = nbck + 1
            ibck(nbck) = jj
            qbck(nbck) = potq(jj)
            vbck(nbck) = potv(jj)
            potq(jj) = potq(jj)*v14scq
            potv(jj) = potv(jj)*v14scv
         end do

cmpi         chunk = ((iats-i)/nproc)
cmpi         ichkh = chunk/4
cmpi         chunk = ichkh*4
cmpi         if (chunk.lt.8) dompi = .false.
cmpi
cmpi         if (dompi) then
cmpi            if (taskid.eq.0) then
cmpi
cmpi
cmpi               offset = chunk + i + 1
cmpi
cmpi               do it=1, nproc - 1
cmpi                  call mpi_send(offset, 1, mpi_integer, it, 5,
cmpi     &              mpi_comm_world, ierr)
cmpi                  offset = offset + chunk
cmpi               end do
cmpi
cmpi               offset = i + 1
cmpi
cmpi            else
cmpi               call mpi_recv(offset, 1, mpi_integer, 0, 5,
cmpi     &            mpi_comm_world, status, ierr)
cmpi            endif
cmpi
cmpi            iend = offset + chunk - 1
cmpi            if (taskid.eq.nproc-1) iend = iats
cmpi
cmpi         else
cmpi
cmpi            if (taskid.eq.0) then
               offset = i + 1
               iend = iats
cmpi            else
cmpi               offset = 1
cmpi               iend = 0
cmpi            endif
cmpi
cmpi         endif

         do k=offset,iend

           do j=1,3
              ded(j) = 0.0e0
           end do

           do j=1,3
              vr(j) = coo(j,i) - coo(j,k)
           end do
           if (box) call reddis(vr)
           r2 = vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3)
           r = sqrt(r2)

           if (r2.le.rc2) then

               rinv = 1.0e0/r
               rinv2 = rinv*rinv
               opt = (iopt(i).eq.1.or.iopt(k).eq.1)


               if (opt) then

                  call apperfc(alfa*r,alerfc,alexp)
c                  alexp1 = 2.0e0*alfa*exp(-1.0e0*alfa2*r2)/spi


                  e = qi * potq(k) * 
     &                     (alerfc*rinv - erfrc + erfrc2*(r-rc))

c                  etmp =  qi *potq(k) * rinv

c                  print*,i,' ',k,' eqc=',e,' etmp=',etmp,
c     &                   alerfc,erfrc,erfrc2*(r-rc),' r-rc=',r-rc

                  de = qi * potq(k) * 
     &                     (alerfc*rinv2 + alexp*rinv - erfrc2)*rinv

                  do j=1,3
                     ded(j) = -de * vr(j)
                  end do

                  ec = ec + e

               endif


               if (opt.and.vdwe1.ne.0.0e0.and.vdwe(k).ne.0.0e0) then

                  i2 = iag(k)
                  kl = ikl(k)

                  if (i1.gt.0.and.i2.gt.0) then
                     pre6(k)  = rs(il,kl)
                     eps(k)   = es(il,kl)
                  else
                     rsum     = vdwr1 + vdwr(k)
                     rsum3    = rsum*rsum*rsum
                     eps(k)   = sqrt(vdwe1 * vdwe(k))
                     pre6(k)  = rsum3*rsum3
                  endif

                  p6 = pre6(k) * rinv2 * rinv2 * rinv2
                  epsm = eps(k)
                  epsm = epsm * potv(k)


                  p12  = p6 * p6

                  e    = epsm * (p12 - 2.0e0*p6)
                  de   = -12.0e0*rinv2*epsm*(p12 - p6)

                  do j=1,3
                     ded(j) = ded(j) + de * vr(j)
                  end do

                  ev   = ev + e
   
               end if

           end if

           do j=1,3
              ftmp(j,i) = ftmp(j,i) + ded(j)
              ftmp(j,k) = ftmp(j,k) - ded(j)
           end do

         end do

         do k=1,nbck
            jj = ibck(k)
            potq(jj) = qbck(k)
            potv(jj) = vbck(k)
         end do

      end do

cmpi      if (taskid.eq.0) then
cmpi
cmpi         do it=1, nproc - 1
cmpi            itmp = it
cmpi
cmpi            call mpi_recv(fftmp(1,1), iat3, mpi_real,
cmpi     &        itmp, 3, mpi_comm_world, status, ierr)
cmpi            do k=1,iatoms
cmpi               do j=1,3
cmpi                  ftmp(j,k) = ftmp(j,k) + fftmp(j,k)
cmpi               end do
cmpi            end do
cmpi
cmpi         end do
cmpi
cmpi         call mpi_reduce(ec, ectmp, 1, mpi_real, mpi_sum,
cmpi     &                0, mpi_comm_world, ierr)
cmpi         call mpi_reduce(ev, evtmp, 1, mpi_real, mpi_sum,
cmpi     &                0, mpi_comm_world, ierr)
cmpi
cmpi         ec = ectmp
cmpi         ev = evtmp
cmpi
cmpi         do it=1, nproc - 1
cmpi            call mpi_send(ectmp, 1, mpi_real, it, 6,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(evtmp, 1, mpi_real, it, 7,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(ftmp(1,1), iat3, 
cmpi     &              mpi_real, it, 8,mpi_comm_world, ierr)
cmpi         end do
cmpi
cmpi      else
cmpi
cmpi         call mpi_send(ftmp(1,1), iat3, mpi_real, 0, 3,
cmpi     &              mpi_comm_world, ierr)
cmpi
cmpi         call mpi_reduce(ec, ectmp, 1, mpi_real, mpi_sum,
cmpi     &              0, mpi_comm_world, ierr)
cmpi         call mpi_reduce(ev, evtmp, 1, mpi_real, mpi_sum, 
cmpi     &              0, mpi_comm_world, ierr)
cmpi
cmpi         call mpi_recv(ectmp, 1, mpi_real, 0, 6,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(evtmp, 1, mpi_real, 0, 7,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(ftmp(1,1), iat3, mpi_real, 
cmpi     &              0, 8, mpi_comm_world, status, ierr)
cmpi
cmpi         ec = ectmp
cmpi         ev = evtmp
cmpi
cmpi      endif


      ecbck =  ec
      evbck =  ev

      do i=1,iatoms
         do j=1,3
            forces(j,i) = forces(j,i) + ftmp(j,i)
            ftmp(j,i) = 0.0e0
         end do
      end do

      if (.not.(fast.and.box)) return
c         call mpi_barrier(mpi_comm_world, ierr)


      ec = 0.0e0
      ev = 0.0e0

      do i=1,numwat

         if (iwtpr(i).eq.1) then
            ii0 = natnow + (i-1)*3

            do l=1,3
         
               ii = ii0 + l
               vdwr1 = vdwr(ii)
               vdwe1 = vdwe(ii)
               i1 = iag(ii)
               il = ikl(ii)
               qi = q(ii)

cmpi               chunk = (iats/nproc)
cmpi               ichkh = chunk/4
cmpi               chunk = ichkh*4
cmpi               if (chunk.lt.8) dompi = .false.
cmpi
cmpi               if (dompi) then
cmpi                  if (taskid.eq.0) then
cmpi
cmpi                     offset = chunk + 1
cmpi
cmpi                     do it=1, nproc - 1
cmpi                        call mpi_send(offset, 1, mpi_integer, it,17,
cmpi     &                    mpi_comm_world, ierr)
cmpi                        offset = offset + chunk
cmpi                     end do
cmpi
cmpi                     offset = 1
cmpi
cmpi                  else
cmpi                     call mpi_recv(offset, 1, mpi_integer, 0,17,
cmpi     &                  mpi_comm_world, status, ierr)
cmpi                  endif
cmpi
cmpi                  iend = offset + chunk - 1
cmpi                  if (taskid.eq.nproc-1) iend = iats
cmpi
cmpi               else
cmpi
cmpi                  if (taskid.eq.0) then
                     offset = 1
                     iend = iats
cmpi                  else
cmpi                     offset = 1
cmpi                     iend = 0
cmpi                  endif
cmpi
cmpi               endif

               do k=offset,iend

                 do j=1,3
                    ded(j) = 0.0e0
                 end do

                 do j=1,3
                    vr(j) = coo(j,ii) - coo(j,k)
                 end do
                 if (box) call reddis(vr)
                 r2 = vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3)
                 r = sqrt(r2)

                 if (r2.le.rc2) then

                     rinv = 1.0e0/r
                     rinv2 = rinv*rinv
                     opt = (iopt(ii).eq.1.or.iopt(k).eq.1)

                     if (opt) then

                        call apperfc(alfa*r,alerfc,alexp)


                        e = qi * potq(k) * 
     &                           (alerfc*rinv - erfrc + erfrc2*(r-rc))


                        de = qi * potq(k) * 
     &                      (alerfc*rinv2 + alexp*rinv - erfrc2)*rinv

                        do j=1,3
                           ded(j) = -de * vr(j)
                        end do

                        ec = ec + e

                     endif


                     if (opt.and.vdwe1.ne.0.0e0
     &                  .and.vdwe(k).ne.0.0e0) then

                        i2 = iag(k)
                        kl = ikl(k)

                        if (i1.gt.0.and.i2.gt.0) then
                           pre6(k)  = rs(il,kl)
                           eps(k)   = es(il,kl)
                        else
                           rsum     = vdwr1 + vdwr(k)
                           rsum3    = rsum*rsum*rsum
                           eps(k)   = sqrt(vdwe1 * vdwe(k))
                           pre6(k)  = rsum3*rsum3
                        endif

                        p6 = pre6(k) * rinv2 * rinv2 * rinv2
                        epsm = eps(k)
                        epsm = epsm * potv(k)


                        p12  = p6 * p6

                        e    = epsm * (p12 - 2.0e0*p6)
                        de   = -12.0e0*rinv2*epsm*(p12 - p6)

                        do j=1,3
                           ded(j) = ded(j) + de * vr(j)
                        end do

                        ev   = ev + e
   
                     endif

                 endif

                 do j=1,3
                    ftmp(j,ii) = ftmp(j,ii) + ded(j)
                    ftmp(j,k) = ftmp(j,k) - ded(j)
                 end do

               end do

            end do
         endif
      end do

cmpi      if (taskid.eq.0) then
cmpi
cmpi         do it=1, nproc - 1
cmpi            itmp = it
cmpi
cmpi            call mpi_recv(fftmp(1,1), iat3, mpi_real,
cmpi     &        itmp,18, mpi_comm_world, status, ierr)
cmpi            do k=1,iatoms
cmpi               do j=1,3
cmpi                  ftmp(j,k) = ftmp(j,k) + fftmp(j,k)
cmpi               end do
cmpi            end do
cmpi
cmpi         end do
cmpi
cmpi         call mpi_reduce(ec, ectmp, 1, mpi_real, mpi_sum,
cmpi     &                0, mpi_comm_world, ierr)
cmpi         call mpi_reduce(ev, evtmp, 1, mpi_real, mpi_sum,
cmpi     &                0, mpi_comm_world, ierr)
cmpi
cmpi         ec = ectmp
cmpi         ev = evtmp
cmpi
cmpi         do it=1, nproc - 1
cmpi            call mpi_send(ectmp, 1, mpi_real, it,19,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(evtmp, 1, mpi_real, it,20,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(ftmp(1,1), iat3, 
cmpi     &              mpi_real, it,21,mpi_comm_world, ierr)
cmpi         end do
cmpi
cmpi      else
cmpi
cmpi         call mpi_send(ftmp(1,1), iat3, mpi_real, 0,18,
cmpi     &              mpi_comm_world, ierr)
cmpi
cmpi         call mpi_reduce(ec, ectmp, 1, mpi_real, mpi_sum,
cmpi     &              0, mpi_comm_world, ierr)
cmpi         call mpi_reduce(ev, evtmp, 1, mpi_real, mpi_sum, 
cmpi     &              0, mpi_comm_world, ierr)
cmpi
cmpi         call mpi_recv(ectmp, 1, mpi_real, 0,19,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(evtmp, 1, mpi_real, 0,20,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(ftmp(1,1), iat3, mpi_real, 
cmpi     &              0,21, mpi_comm_world, status, ierr)
cmpi
cmpi         ec = ectmp
cmpi         ev = evtmp
cmpi
cmpi      endif

      ec = ec + ecbck
      ev = ev + evbck

      do i=1,iatoms
         do j=1,3
            forces(j,i) = forces(j,i) + ftmp(j,i)
         end do
      end do

      return
      end

      subroutine watcnd(coo)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxion=2000)
      common /h2oer/ numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      logical dowat
      dimension coo(3,*),vr(3)

      rcut = 7.86e0
      rcut2 = rcut*rcut
      zero = 0.0e0
      niwat = 0
      iaddr = natnow

      do i=1,numwat

         jaddr = iaddr + 3

         do j=i+1,numwat

           kaddr = iaddr + 1
           laddr = jaddr + 1

           do m=1,3
              vr(m) = coo(m,kaddr) - coo(m,laddr)
           end do
           
           call reddis(vr)

           r2 = zero
           do m=1,3
              r2 = r2 + vr(m)*vr(m)
           end do
           
           if (r2.lt.rcut2) niwat = niwat + 1

           jaddr = jaddr + 3

         end do

         iaddr =  iaddr + 3

      end do

      return
      end

      subroutine watlsd(coo,iwtpr,listw)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxion=2000)
      common /h2oer/ numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      parameter (mxproc=1000)
      common /mpiwat/ iprb(mxproc),ipre(mxproc),iniwat
      integer taskid,offset
      common /mpih/ nproc,taskid
      logical dowat
      dimension coo(3,*),vr(3),listw(2,*),iwtpr(*)

      rcut = 7.86e0
      rcut2 = rcut*rcut
      zero = 0.0e0
      niwat = 0
      iaddr = natnow

      do i=1,numwat

         kaddr = iaddr + 1
         jaddr = iaddr + 3

         do j=i+1,numwat

           laddr = jaddr + 1

           do m=1,3
              vr(m) = coo(m,kaddr) - coo(m,laddr)
           end do
           
           call reddis(vr)

           r2 = zero
           do m=1,3
              r2 = r2 + vr(m)*vr(m)
           end do
           
           if (r2.lt.rcut2) then
               niwat = niwat + 1
               listw(1,niwat) = iaddr + 1
               listw(2,niwat) = jaddr + 1
           endif

           jaddr = jaddr + 3


         end do

         iwtpr(i) = 0

         do j=1,natnow

            do k=1,3
               vr(k) = coo(k,iaddr+1) - coo(k,j)
            end do

            call reddis(vr)

            rab2 = vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3)
            if (rab2.lt.103.0e0) iwtpr(i) = 1

         end do

         iaddr =  iaddr + 3

      end do

      iniwat = 0

c      print*,"Updated water list"

      return
      end

      subroutine fstwad(ew,coo,ftmp,q,fftmp,listw,
     &                  woo,whh,woh,dewoo,dewhh,dewoh)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxac=3*mxcon)
      common /athlp/  iatoms, mxnat
      parameter (mxion=2000)
      common /h2oer/ numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      parameter (mxcls=50)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      integer ambcls
      parameter (mxamb=1590)
cmpi      include 'mpif.h'
cmpi      integer chunk
cmpi      integer status(mpi_status_size)
      logical dompi
      integer taskid
      common /mpih/ nproc,taskid
      parameter (mxproc=1000)
      common /mpiwat/ iprb(mxproc),ipre(mxproc),iniwat
      logical box, cell, fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      dimension coo(3,*),vr(3),q(*),listw(2,*),
     &          ftmp(3,*),fftmp(3,*)
      dimension woo(*),whh(*),woh(*),dewoo(*),dewhh(*),dewoh(*)

c
c Adaptation from:
c
c Cristopher J. Fennell and J. Daniel Gezelter,
c Journal of Chemical Physics. 124, 234104 (2006)
c
      dompi = .false.

cmpi      call mpi_comm_rank(mpi_comm_world, taskid, ierr)
cmpi      dompi = .true.
cmpi      call iniprc

      ew   = 0.0e0
      zero = 0.0e0
      thou = 1000.0e0

      rcut = 7.86e0
      rcut2 = rcut*rcut
      rcinv = 1.0e0 / rcut
      rc2inv = rcinv*rcinv

c interpolation Fapprox = (1-x)*F(i) + x*F(i+1)
c where i is: i = int(R/1000)
c where x is remainder: x = R - int(R/1000)

      iat3 = iatoms*3

      if (dompi) then
         ibeg = iprb(taskid+1)
         iend = ipre(taskid+1)
         iaddr = natnow + (ibeg-1)*3
      else
         ibeg = 1
         iend = niwat
      endif

      do i=ibeg,iend
       
            kaddr = listw(1,i)
            jaddr = listw(2,i)
            laddr = jaddr

c    O w1 <->  O w2
       
            do k=1,3
               vr(k) = coo(k,kaddr) - coo(k,laddr)
            end do
           
            call reddis(vr)

            r2 = zero
            do k=1,3
               r2 = r2 + vr(k)*vr(k)
            end do
           

            r = sqrt(r2) 
            if (r.le.8.999e0) then
               r1000 = r*thou
               ir = int(r1000)
               x = (r1000 - dble(ir))
               w = (1.0e0-x)*woo(ir) + x*woo(ir+1)

               ew = ew + w
               de = dewoo(ir)
       
               do k=1,3
                  ftmp(k,kaddr) = ftmp(k,kaddr) + de*vr(k)
                  ftmp(k,laddr) = ftmp(k,laddr) - de*vr(k)
               end do
       

c    O w1 <->  H1 w2
       
               laddr = laddr + 1
               
               do k=1,3
                  vr(k) = coo(k,kaddr) - coo(k,laddr)
               end do
           
               call reddis(vr)

               r2 = zero
               do k=1,3
                  r2 = r2 + vr(k)*vr(k)
               end do
           

               r = sqrt(r2) 
               if (r.le.8.999e0) then
                  r1000 = r*thou
                  ir = int(r1000)
                  x = (r1000 - dble(ir))
                  w = (1.0e0-x)*woh(ir) + x*woh(ir+1)
                  ew = ew + w
                  de = dewoh(ir)
       
                  do k=1,3
                     ftmp(k,kaddr) = ftmp(k,kaddr) + de*vr(k)
                     ftmp(k,laddr) = ftmp(k,laddr) - de*vr(k)
                  end do
               endif
       

c    O w1 <->  H2 w2
       
               laddr = laddr + 1
            
               do k=1,3
                  vr(k) = coo(k,kaddr) - coo(k,laddr)
               end do
           
               call reddis(vr)

               r2 = zero
               do k=1,3
                  r2 = r2 + vr(k)*vr(k)
               end do
           

               r = sqrt(r2) 
               if (r.le.8.999e0) then
                  r1000 = r*thou
                  ir = int(r1000)
                  x = (r1000 - dble(ir))
                  w = (1.0e0-x)*woh(ir) + x*woh(ir+1)
                  ew = ew + w
                  de = dewoh(ir)
       
                  do k=1,3
                     ftmp(k,kaddr) = ftmp(k,kaddr) + de*vr(k)
                     ftmp(k,laddr) = ftmp(k,laddr) - de*vr(k)
                  end do
               endif
       

c    H1 w1 <->  O w2
       
               kaddr = kaddr + 1
             
               laddr = jaddr
       
               do k=1,3
                  vr(k) = coo(k,kaddr) - coo(k,laddr)
               end do
           
               call reddis(vr)

               r2 = zero
               do k=1,3
                  r2 = r2 + vr(k)*vr(k)
               end do
           

               r = sqrt(r2) 
               if (r.le.8.999e0) then
                  r1000 = r*thou
                  ir = int(r1000)
                  x = (r1000 - dble(ir))
                  w = (1.0e0-x)*woh(ir) + x*woh(ir+1)
                  ew = ew + w
                  de = dewoh(ir)
       
                  do k=1,3
                     ftmp(k,kaddr) = ftmp(k,kaddr) + de*vr(k)
                     ftmp(k,laddr) = ftmp(k,laddr) - de*vr(k)
                  end do
               endif
       

c    H1 w1 <->  H1 w2
       
               laddr = laddr + 1
            
               do k=1,3
                  vr(k) = coo(k,kaddr) - coo(k,laddr)
               end do
           
               call reddis(vr)

               r2 = zero
               do k=1,3
                  r2 = r2 + vr(k)*vr(k)
               end do
           
               r = sqrt(r2) 
               if (r.le.8.999e0) then
                  r1000 = r*thou
                  ir = int(r1000)
                  x = (r1000 - dble(ir))
                  w = (1.0e0-x)*whh(ir) + x*whh(ir+1)
                  ew = ew + w
                  de = dewhh(ir)
       
                  do k=1,3
                     ftmp(k,kaddr) = ftmp(k,kaddr) + de*vr(k)
                     ftmp(k,laddr) = ftmp(k,laddr) - de*vr(k)
                  end do
               endif

       
c    H1 w1 <->  H2 w2
       
               laddr = laddr + 1
            
               do k=1,3
                  vr(k) = coo(k,kaddr) - coo(k,laddr)
               end do
           
               call reddis(vr)

               r2 = zero
               do k=1,3
                  r2 = r2 + vr(k)*vr(k)
               end do
           
               r = sqrt(r2) 
               if (r.le.8.999e0) then
                  r1000 = r*thou
                  ir = int(r1000)
                  x = (r1000 - dble(ir))
                  w = (1.0e0-x)*whh(ir) + x*whh(ir+1)
                  ew = ew + w
                  de = dewhh(ir)
       
                  do k=1,3
                     ftmp(k,kaddr) = ftmp(k,kaddr) + de*vr(k)
                     ftmp(k,laddr) = ftmp(k,laddr) - de*vr(k)
                  end do
               endif
       

c    H2 w1 <->  O w2
       
               kaddr = kaddr + 1
                
               laddr = jaddr
       
               do k=1,3
                  vr(k) = coo(k,kaddr) - coo(k,laddr)
               end do
           
               call reddis(vr)

               r2 = zero
               do k=1,3
                  r2 = r2 + vr(k)*vr(k)
               end do
           
               r = sqrt(r2) 
               if (r.le.8.999e0) then
                  r1000 = r*thou
                  ir = int(r1000)
                  x = (r1000 - dble(ir))
                  w = (1.0e0-x)*woh(ir) + x*woh(ir+1)
                  ew = ew + w
                  de = dewoh(ir)
       
                  do k=1,3
                     ftmp(k,kaddr) = ftmp(k,kaddr) + de*vr(k)
                     ftmp(k,laddr) = ftmp(k,laddr) - de*vr(k)
                  end do
               endif
       

c    H2 w1 <->  H1 w2
       
               laddr = laddr + 1
           
               do k=1,3
                  vr(k) = coo(k,kaddr) - coo(k,laddr)
               end do
           
               call reddis(vr)

               r2 = zero
               do k=1,3
                  r2 = r2 + vr(k)*vr(k)
               end do
              
               r = sqrt(r2) 
               if (r.le.8.999e0) then
                  r1000 = r*thou
                  ir = int(r1000)
                  x = (r1000 - dble(ir))
                  w = (1.0e0-x)*whh(ir) + x*whh(ir+1)
                  ew = ew + w
                  de = dewhh(ir)
       
                  do k=1,3
                     ftmp(k,kaddr) = ftmp(k,kaddr) + de*vr(k)
                     ftmp(k,laddr) = ftmp(k,laddr) - de*vr(k)
                  end do
               endif
       

c    H2 w1 <->  H2 w2
       
               laddr = laddr + 1
               
               do k=1,3
                  vr(k) = coo(k,kaddr) - coo(k,laddr)
               end do
              
               call reddis(vr)

               r2 = zero
               do k=1,3
                  r2 = r2 + vr(k)*vr(k)
               end do
           
               r = sqrt(r2) 
               if (r.le.8.999e0) then
                  r1000 = r*thou
                  ir = int(r1000)
                  x = (r1000 - dble(ir))
                  w = (1.0e0-x)*whh(ir) + x*whh(ir+1)
                  ew = ew + w
                  de = dewhh(ir)
       
                  do k=1,3
                     ftmp(k,kaddr) = ftmp(k,kaddr) + de*vr(k)
                     ftmp(k,laddr) = ftmp(k,laddr) - de*vr(k)
                  end do
               endif
       
            endif

      end do
           
cmpi      if (taskid.eq.0) then
cmpi
cmpi         do it=1, nproc - 1
cmpi            itmp = it
cmpi
cmpi            call mpi_recv(fftmp(1,1), iat3, mpi_real,
cmpi     &        itmp, 10, mpi_comm_world, status, ierr)
cmpi            do k=1,iatoms
cmpi               do j=1,3
cmpi                  ftmp(j,k) = ftmp(j,k) + fftmp(j,k)
cmpi               end do
cmpi            end do
cmpi
cmpi         end do
cmpi
cmpi         call mpi_reduce(ew, ewtmp, 1, mpi_real, mpi_sum,
cmpi     &                0, mpi_comm_world, ierr)
cmpi
cmpi         ew = ewtmp
cmpi
cmpi         do it=1, nproc - 1
cmpi            call mpi_send(ewtmp, 1, mpi_real, it, 13,
cmpi     &              mpi_comm_world, ierr)
cmpi            call mpi_send(ftmp(1,1), iat3, 
cmpi     &              mpi_real, it, 14,mpi_comm_world, ierr)
cmpi         end do
cmpi
cmpi      else
cmpi
cmpi         call mpi_send(ftmp(1,1), iat3, mpi_real, 0, 10,
cmpi     &              mpi_comm_world, ierr)
cmpi
cmpi         call mpi_reduce(ew, ewtmp, 1, mpi_real, mpi_sum,
cmpi     &              0, mpi_comm_world, ierr)
cmpi
cmpi         call mpi_recv(ewtmp, 1, mpi_real, 0, 13,
cmpi     &              mpi_comm_world, status, ierr)
cmpi         call mpi_recv(ftmp(1,1), iat3, mpi_real, 
cmpi     &              0, 14, mpi_comm_world, status, ierr)
cmpi
cmpi         ew = ewtmp
cmpi
cmpi      endif

      return
      end

      subroutine iniprc
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxion=2000)
      common /h2oer/ numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      parameter (mxproc=1000)
      common /mpiwat/ iprb(mxproc),ipre(mxproc),iniwat
      integer taskid,chunk
      common /mpih/ nproc,taskid

      if (iniwat.eq.0) then
         chunk = niwat /nproc
         do i=1,nproc
            iprb(i) = (i-1)*chunk + 1
            ipre(i) = i*chunk
         end do
         ipre(nproc) = niwat
         iniwat = 1
      endif
           
      return
      end

