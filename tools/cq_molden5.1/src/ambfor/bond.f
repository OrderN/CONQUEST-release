      subroutine gtbpar(istat,i1,i2,bl,bk)
      implicit real (a-h,o-z), integer (i-n)
      integer*2 i1, i2
      parameter (mxbnd=86)
      integer bndcon
      common /bndpar/ bnd(2,mxbnd),bndcon(2,mxbnd)
      parameter (mxamb=1590)
      integer ambcls
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxcls=50)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      logical doamb, dogaff
      parameter (mxgff=72)
      parameter (mxgbnd=747)
      integer gfbcon
      common /gfbpar/ gfbnd(2,mxgbnd),gfbcon(2,mxgbnd)
      integer amb2gf
      common /gftyp/  gfvdw(2,mxgff),amb2gf(mxamb)
      parameter (mxcor=42)
      common /corr/   iatcor(2,mxcor)
      logical docora, docorb, doinit
      
      istat = 0
      doamb = .false.
      dogaff = .false.
      docora = .true.
      docorb = .false.
      doinit = .true.

      if (i1.gt.0.and.i2.gt.0) then
         doamb = .true.
      elseif (i1.lt.0.or.i2.lt.0) then
         dogaff = .true.
      endif

      if (doamb) then
         iac = ambcls(i1)
         ibc = ambcls(i2)

         do k=1,mxbnd 
            if (
     &       (iac.eq.bndcon(1,k).and.ibc.eq.bndcon(2,k))
     &          .or.
     &       (ibc.eq.bndcon(1,k).and.iac.eq.bndcon(2,k))
     &      ) then
               istat = 1
               bl = bnd(2,k)
               bk = bnd(1,k)
            endif
         end do

         if (istat.eq.0.and.(i1.ge.1254.and.i2.ge.1254)) then
            iac = amb2gf(i1)
            ibc = amb2gf(i2)
            dogaff = .true.
            doinit = .false.
         endif

      endif

      if (dogaff) then

         if (doinit) then
            if (i1.gt.0) then
               iac = mapagf(ambcls(i1))
            else
               iac = abs(i1)
               if (iac.eq.72) iac = 6
            endif
            if (i2.gt.0) then
               ibc = mapagf(ambcls(i2))
            else
               ibc = abs(i2)
               if (ibc.eq.72) ibc = 6
            endif
         endif

10       continue

         do k=1,mxgbnd 
            if (
     &       (iac.eq.gfbcon(1,k).and.ibc.eq.gfbcon(2,k))
     &          .or.
     &       (ibc.eq.gfbcon(1,k).and.iac.eq.gfbcon(2,k))
     &      ) then
               istat = 1
               bl = gfbnd(2,k)
               bk = gfbnd(1,k)
            endif
         end do

         if (istat.ne.1.and.docora) then

            do i=1,mxcor
               if (iac.eq.iatcor(1,i)) iac = iatcor(2,i)
            end do

            docora = .false.
            docorb = .true.
            goto 10
         endif

         if (istat.ne.1.and.docorb) then

            do i=1,mxcor
               if (ibc.eq.iatcor(1,i)) ibc = iatcor(2,i)
            end do

            docorb = .false.
            goto 10
         endif

      endif

      return
      end

      subroutine bndard(istat,nbnd,ibnd,bl,bk,maxbnd,iconn,ityp,iopt)
c
c     create bond list from the connection array and assign parameters
c
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      common /athlp/  iatoms, mxnat
      parameter (mxbnd=86)
      integer bndcon
      common /bndpar/ bnd(2,mxbnd),bndcon(2,mxbnd)
      parameter (mxamb=1590)
      integer ambcls
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      integer*2 ityp
      dimension iconn(mxcon+1,*),ityp(*),iopt(*),
     &          ibnd(2,*),bl(*),bk(*)

      nbnd = 0
      istat = 1

      do i=1,iatoms
         do j=1,iconn(1,i)
            jj = iconn(j+1,i)
            if (jj.gt.0) then
               if (i.lt.jj.and.(iopt(i).eq.1.or.iopt(jj).eq.1)) then
                  if (nbnd.lt.maxbnd) then
                     nbnd = nbnd + 1
                     ibnd(1,nbnd) = i
                     ibnd(2,nbnd) = jj
                     call gtbpar(istat,
     &                    ityp(i),ityp(jj),bl(nbnd),bk(nbnd))
                     if (istat.eq.0) write(iun5,*) 
     &                       'Missing Bond parameter: ',
     &                               abs(ityp(i)),abs(ityp(jj))
                  else
                     istat = 0
                  endif
               endif
            endif
         end do
      end do

      return
      end

      subroutine bond(eb,nbnd,ibnd,bl,bk,coo,forces)
      implicit real (a-h,o-z), integer (i-n)
      common /athlp/  iatoms, mxnat
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast

      dimension ibnd(2,*), bl(*), bk(*), coo(3,*), forces(3,*)
      dimension vab(3)

      eb = 0.0e0

      do i = 1, nbnd
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         bopt = bl(i)
         fk = bk(i)

            do j=1,3
               vab(j) = coo(j,ia) - coo(j,ib)
            end do

            if (box) call reddis(vab)

            rab = vlen(vab)
            dt = rab - bopt

            dt2 = dt * dt
            e = fk * dt2 
            deddt = 2.0e0 * fk * dt

            de = deddt / rab

            call vscal(vab,3,de)

            eb = eb + e
          
            do j=1,3
               forces(j,ia) = forces(j,ia) + vab(j)
               forces(j,ib) = forces(j,ib) - vab(j)
            end do

      end do

      return
      end
