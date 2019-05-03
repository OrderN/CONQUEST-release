      subroutine gtapar(istat,swap,tryabc,i1,i2,i3,al,ak)
      implicit real (a-h,o-z), integer (i-n)
      integer*2 i1,i2,i3
      parameter (mxamb=1590)
      integer ambcls,angcon,gfacon
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxcls=50)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      parameter (mxang=199)
      common /angpar/ ang(2,mxang),angcon(3,mxang)
      parameter (mxgff=72)
      parameter (mxgang=3506)
      common /gfapar/ gfang(2,mxgang),angint(2,mxgff),
     &                gfacon(3,mxgang)
      integer amb2gf
      common /gftyp/  gfvdw(2,mxgff),amb2gf(mxamb)
      character*3 ambstr
      character*2 gffstr
      common /ffstr/  ambstr(mxamb), gffstr(mxgff)
      parameter (mxcor=42)
      common /corr/   iatcor(2,mxcor)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      logical swap, tryabc, doamb, dogaff
      logical docora, docorb, docorc, docrab, docrbc, docrac,
     &        docabc, doinit
      
      istat = 0
      doamb = .false.
      dogaff = .false.

      docora = .true.
      docorb = .true.
      docorc = .true.
      docrab = .true.
      docrbc = .true.
      docrac = .true.
      docabc = .true.
      doinit = .true.

      if (i1.gt.0 .and. i2.gt.0 .and. i3.gt.0) then
         doamb = .true.
      elseif (i1.lt.0 .or. i2.lt.0 .or.i3.lt.0) then
         dogaff = .true.
      endif

      if (doamb) then
         if (swap) then
            iac = ambcls(i1)
            ibc = ambcls(i2)
            icc = ambcls(i3)
         else
            iac = ambcls(i3)
            ibc = ambcls(i2)
            icc = ambcls(i1)
         endif

         if (iac.gt.icc) then
            ict = iac
            iac = icc
            icc = ict
         endif

         do l=1,mxang
            if ( iac.eq.angcon(1,l).and.ibc.eq.angcon(2,l)
     &    .and.  icc.eq.angcon(3,l) ) then
               istat = 1
               ak = ang(1,l)
               al = ang(2,l)
            endif
         end do

         if (istat.eq.0.and.
     &      (i1.ge.1254.and.i2.ge.1254.and.i3.ge.1254)) then
            if (swap) then
               iac = amb2gf(i1)
               ibc = amb2gf(i2)
               icc = amb2gf(i3)
            else
               iac = amb2gf(i3)
               ibc = amb2gf(i2)
               icc = amb2gf(i1)
            endif

            dogaff = .true.
            doinit = .false.
         endif

      endif

      if (dogaff) then

         if (doinit) then

            if (swap) then
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
               if (i3.gt.0) then
                  icc = mapagf(ambcls(i3))
               else
                  icc = abs(i3)
                  if (icc.eq.72) icc = 6
               endif
            else
               if (i3.gt.0) then
                  iac = mapagf(ambcls(i3))
               else
                  iac = abs(i3)
                  if (iac.eq.72) iac = 6
               endif
               if (i2.gt.0) then
                  ibc = mapagf(ambcls(i2))
               else
                  ibc = abs(i2)
                  if (ibc.eq.72) ibc = 6
               endif
               if (i1.gt.0) then
                  icc = mapagf(ambcls(i1))
               else
                  icc = abs(i1)
                  if (icc.eq.72) icc = 6
               endif
            endif
         endif

         if (iac.gt.icc) then
            ict = iac
            iac = icc
            icc = ict
         endif

         iact = iac
         ibct = ibc
         icct = icc

         irun = 0
10       continue

         irun = irun + 1
         do l=1,mxgang
            if ( iac.eq.gfacon(1,l).and.ibc.eq.gfacon(2,l)
     &    .and.  icc.eq.gfacon(3,l) ) then
               istat = 1
               ak = gfang(1,l)
               al = gfang(2,l)
            endif
         end do

         if (istat.ne.1.and.docorb) then

            iac = iact
            ibc = ibct
            icc = icct
            ibcs = 0

            do i=1,mxcor
               if (ibc.eq.iatcor(1,i)) ibc = iatcor(2,i)
            end do

            if (ibc.ne.ibct) ibcs = ibc

            docorb = .false.
            if (ibcs.ne.0) goto 10
         endif

         if (istat.ne.1.and.docora) then

            iac = iact
            ibc = ibct
            icc = icct
            iacs = 0

            do i=1,mxcor
               if (iac.eq.iatcor(1,i)) iac = iatcor(2,i)
            end do

            if (iac.ne.iact) iacs = iac

            if (iac.gt.icc) then
               ict = iac
               iac = icc
               icc = ict
            endif


            docora = .false.
            if (iacs.ne.0) goto 10
         endif

         if (istat.ne.1.and.docorc) then

            iac = iact
            ibc = ibct
            icc = icct
            iccs = 0

            do i=1,mxcor
               if (icc.eq.iatcor(1,i)) icc = iatcor(2,i)
            end do

            if (icc.ne.icct) iccs = icc

            if (iac.gt.icc) then
               ict = iac
               iac = icc
               icc = ict
            endif

            docorc = .false.
            if (iccs.ne.0) goto 10

         endif

         if (istat.ne.1.and.docrab) then

            docrab = .false.
            if (iacs.ne.0.and.ibcs.ne.0) then
               iac = iacs
               ibc = ibcs
               icc = icct

               if (iac.gt.icc) then
                  ict = iac
                  iac = icc
                  icc = ict
               endif

               goto 10
            endif

         endif

         if (istat.ne.1.and.docrbc) then

            docrbc = .false.

            if (ibcs.ne.0.and.iccs.ne.0) then
               iac = iact
               ibc = ibcs
               icc = iccs

               if (iac.gt.icc) then
                  ict = iac
                  iac = icc
                  icc = ict
               endif

               goto 10
            endif

         endif

         if (istat.ne.1.and.docrac) then

            docrac = .false.

            if (iacs.ne.0.and.iccs.ne.0) then
               iac = iacs
               ibc = ibct
               icc = iccs

               if (iac.gt.icc) then
                  ict = iac
                  iac = icc
                  icc = ict
               endif

               goto 10
            endif

         endif

         if (istat.ne.1.and.docabc) then

            docabc = .false.

            if (iacs.ne.0.and.ibcs.ne.0.and.iccs.ne.0) then
               iac = iacs
               ibc = ibcs
               icc = iccs

               if (iac.gt.icc) then
                  ict = iac
                  iac = icc
                  icc = ict
               endif

               goto 10
            endif

         endif

         iac = iact
         ibc = ibct
         icc = icct

         if (istat.eq.0.and.tryabc) then
            call parabc(istat,i1,i2,i3,al,ak)
         endif

      endif

      if (istat.eq.0) then

         if (doamb) then
            write(iun5,'(a,3(i4,1x),a,3(a3,1x))') 
     &         'amber type ', i1,i2,i3,
     &         ' label ',ambstr(i1),ambstr(i2),ambstr(i3)
         endif

         if (dogaff) then
            if (swap) then
               write(iun5,'(a,3(i4,1x),a,3(a3,1x),a,3(i5,1x))') 
     &         'type ',i3,i2,i1,
     &         ' label ',gffstr(iac),gffstr(ibc),gffstr(icc),
     &         ' gaff type ',iac,ibc,icc
            else
               write(iun5,'(a,3(i4,1x),a,3(a3,1x),a,3(i5,1x))') 
     &         'type ',i1,i2,i3,
     &         ' label ',gffstr(iac),gffstr(ibc),gffstr(icc),
     &         ' gaff type ',iac,ibc,icc
            endif
         endif
      endif

      return
      end

      subroutine parabc(istat,ia,ib,ic,angabc,fabc)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxgff=72)
      parameter (mxgbnd=747)
      integer gfbcon,gfacon
      common /gfbpar/ gfbnd(2,mxgbnd),gfbcon(2,mxgbnd)
      parameter (mxgang=3506)
      common /gfapar/ gfang(2,mxgang),angint(2,mxgff),
     &                gfacon(3,mxgang)
      integer*2 ia,ib,ic

      todeg = 45.0e0 / atan(1.0e0)
      istat = 0

      call gtapar(istat,.true.,.false.,ia,ib,ia,angaba,ak)
      if (istat.eq.0) return
      call gtapar(istat,.true.,.false.,ic,ib,ic,angcbc,ak)
      if (istat.eq.0) return

      angabc = (angaba + angcbc) / 2.0e0
      sanabc = sqrt(angabc / todeg)

      call gtbpar(istat,ia,ib,rab,bk)
      if (istat.eq.0) return
      call gtbpar(istat,ib,ic,rbc,bk)
      if (istat.eq.0) return

c
c formula gaff publication
c
c      d = (rab - rbc)**2 / (rab + rbc)**2
c
c as in parmchk:
c
      d = (rab-rbc)*(rab-rbc)
      d = d / (rab+rbc)*(rab+rbc)

      aa = angint(2,abs(ia))
      ab = angint(1,abs(ib))
      ac = angint(2,abs(ic))
c
c formula gaff publication
c
c      fabc = todeg*todeg*143.9e0*aa*ab*ac*
c     &       exp(-2.0e0*d) / ( (rab+rbc)*angabc*angabc )

c
c as in parmchk:
c
      fabc = 143.9e0*aa*ab*ac*
     &       exp(-2.0e0*d) / ( (rab+rbc)*sanabc )
      return
      end

      subroutine angard(istat,nang,iang,ango,ak,maxang,iconn,ityp,iopt)
c
c     create angle list from the connection array and assign parameters
c
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      common /athlp/  iatoms, mxnat
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      integer*2 ityp
      dimension iconn(mxcon+1,*),ityp(*),iopt(*),
     &          iang(3,*),ango(*),ak(*)

      nang = 0
      istat = 1

      do i=1,iatoms

         do j=1,iconn(1,i) - 1
            jj = iconn(j+1,i)
            if (jj.gt.0) then

               do k=j+1,iconn(1,i)
                kk = iconn(k+1,i)
                if (kk.gt.0.and.kk.ne.jj.and.
     &          (iopt(i).eq.1.or.iopt(jj).eq.1.or.iopt(kk).eq.1)) then

                  if (nang.lt.maxang) then
                     nang = nang + 1
                     if (kk.ge.jj) then
                        iang(1,nang) = jj
                        iang(2,nang) = i
                        iang(3,nang) = kk
                     else
                        iang(1,nang) = kk
                        iang(2,nang) = i
                        iang(3,nang) = jj
                     endif

                     call gtapar(istat,
     &                    kk.ge.jj,.true.,ityp(jj),ityp(i),ityp(kk),
     &                    ango(nang),ak(nang))
                     if (istat.eq.0) then
                          write(iun5,*) 'angle parameter not found: ',
     &                    abs(ityp(jj)),abs(ityp(i)),abs(ityp(kk)),
     &                    ' atnr ',jj,i,kk

                     endif

                  else

                     istat = 0

                  endif

                endif
               end do

            endif

         end do

      end do

      return
      end

      subroutine angle(ea,nang,iang,ango,ak,coo,forces)
      implicit real (a-h,o-z), integer (i-n)
      common /athlp/  iatoms, mxnat
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      dimension iang(3,*),ango(*),ak(*),coo(3,*),forces(3,*)
      dimension va(3),vb(3),vc(3),vab(3),vcb(3),vp(3),
     &          deda(3),dedb(3),dedc(3)

c initialize

      ea = 0.0e0
      todeg = 45.0e0 / atan(1.0e0)
      angu = 1.0e0 / (todeg*todeg)

      do i=1,nang

         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)

         aopt = ango(i)
         fk   = ak(i)

            do j=1,3
               va(j) = coo(j,ia)
               vb(j) = coo(j,ib)
               vc(j) = coo(j,ic)
            end do

            call vsub(va,vb,vab,3)
            call vsub(vc,vb,vcb,3)

            if (box) then
               call reddis(vab)
               call reddis(vcb)
            endif

            rab2 = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
            rcb2 = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)

            if (rab2.ne.0.0e0 .and. rcb2.ne.0.0e0) then
                  call crprod(vcb,vab,vp)
                  rp = vlen(vp)
                  rp = max(rp,0.000001e0)
                  call impsc(vab,vcb,cosa)
                  cosa = dble(min(1.0e0,max(-1.0e0,cosa)))
                  angl = todeg * acos(cosa)

                  dt  = angl - aopt
                  dt2 = dt * dt
                  e   = angu * fk * dt2
                  deddt = 2.d0 * fk * dt * angu * todeg

                  call crprod(vab,vp,deda)
                  rmul1 = -deddt / (rab2*rp)
                  call vscal(deda,3,rmul1)

                  call crprod(vcb,vp,dedc)
                  rmul2 =  deddt / (rcb2*rp)
                  call vscal(dedc,3,rmul2)

                  call vadd(deda,dedc,dedb,3)
                  call vscal(dedb,3,-1.0e0)

                  ea = ea + e
                  do j=1,3
                     forces(j,ia) = forces(j,ia) + deda(j)
                     forces(j,ib) = forces(j,ib) + dedb(j)
                     forces(j,ic) = forces(j,ic) + dedc(j)
                  end do
            end if
      end do

      return
      end
