      subroutine pmfasd(iopt,idochg,
     &                  ianz,iresid,iconn,ityp,
     &                  ncalf,iamino)
      implicit double precision (a-h,o-z)
      parameter (mxcon=10)
      parameter (numcal=50000)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxliga=200)
      parameter (mxewin=100)
      parameter (mxpmfl=5)
      common /athlp/ iatoms, mxnat
      common /pmfres/ ewin(mxewin),npmfs
      real pmf,dfire
      common /pmfpar/ ipmfrs(20,20),lpmfrs(20,20),ipmfl(20),
     &                pmf(16,26,60),dfire(13,13,20)
c
c pmlig   contains the total pmf score of this ligand atom
c
c ipmlig  contains index of ligand atoms with pmf types
c mpmlig  the length of pmlig
c ipmres  contains index to residues with pmf types
c npmres  the number of residues
c
c ityp    is used to contain the protein and ligand pmf types
c 
c all residues have pmf assigned but only those close enough to
c the ligand have there indexes kept in ipmres
c
      common /pmfrk/ pmlig(mxliga),ipmlig(mxliga),npmlig,
     &               ipmres(numcal),npmres
      common /pmflvl/ ipmfm,ipmfh,pmflev(mxpmfl),levcol(mxpmfl)
      integer*2 ityp
      common /types/ iff
      common /align/ vecs(3,3),nscnd,iscst,ialtyp,iocnt
      common /cllmat/ rdum(12),tz(3),tzorg(3),itz
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension ianz(*),iresid(*),ityp(*),iconn(mxcon+1,*)
      dimension iamino(*)
      data levcol /1,6,3,5,5/

      toang = 0.52917706d0

c iopt = 0, initialize, assign protein and ligand pmf types
c iopt = 1, only reassign ligand pmf types, and assign which
c           residues interact

      if (iopt.eq.1) then
         do i=iscst+1,iscst+nscnd
            ityp(i) = 0
         end do
      else
         do i=1,iatoms
            ityp(i) = 0
         end do
      endif

      npmlig = 0

      if (iopt.eq.1) goto 100

      call inipmf

      iff = 8


      do i=1,ncalf
         call getpdb(i,ipdb,ihpdb)
         iam = iamino(i)
         if (iam.ge.1.and.iam.le.20) then
            do j=1,mxsym
               if (ipdb(j).gt.0) then
c Backbone atoms 1-4
                   if (j.eq.1) then
                      ityp(ipdb(j)) = 6
                   elseif (j.eq.2) then
                      ityp(ipdb(j)) = 2
                   elseif (j.eq.3) then
                      ityp(ipdb(j)) = 2
                   elseif (j.eq.4) then
                      ityp(ipdb(j)) = 11
                   else
c SideChain atoms
                      do k=1,ipmfl(iam)
                         if (j.eq.ipmfrs(iam,(k-1)*2+1)) then
                            ityp(ipdb(j)) = ipmfrs(iam,(k-1)*2+2)
                            goto 10
                         endif
                      end do
10                    continue
                   endif
               endif
            end do
            do j=1,mxhsym*3
               if (ihpdb(j).gt.0) then
                   ityp(ihpdb(j)) = 16
               endif
            end do
         endif
      end do

100   continue

c Ligand types

      do i=1,iatoms
         if (iresid(i).le.-4.and.ityp(i).eq.0) then
c exclude waters
            if (.not.(ianz(i).eq.8.and.iconn(1,i).eq.0)) then
               call ipmtyp(ipmt,i,ianz(i),idochg)
               ityp(i) = ipmt
               if (ityp(i).ne.0) then
                  if (npmlig.lt.mxliga) then
                     npmlig = npmlig + 1
                     ipmlig(npmlig) = i
                  else
                     print*,'array to hold ligand atoms full'
                  endif
               endif
            endif
         endif
      end do

      call updres

      if (iopt.eq.1) return

      npmfs = 0

      call parptr(115,ewin,fdum,idum)
      call parptr(116,fdum,fdum,npmfs)
      call parptr(123,fdum,fdum,ipmfm)
      call parptr(124,fdum,fdum,levcol)
      call parptr(125,fdum,fdum,ipmfh)
      call crsco()

      return
      end

      subroutine dfiasd(ityp,ncalf,iamino)
      implicit double precision (a-h,o-z)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxdres=40)
      parameter (mxrot=81)
      common /athlp/ iatoms, mxnat
      real pmf,dfire
      common /pmfpar/ ipmfrs(20,20),lpmfrs(20,20),ipmfl(20),
     &                pmf(16,26,60),dfire(13,13,20)
      integer*2 ityp
      common /types/ iff
      common /dfiwrk/estat(mxdres,mxrot),ndrs,idrs(mxdres),lrot(mxdres),
     &               lpmres(mxdres,mxdres),nlpmrs(mxdres),
     &               kpmres(mxdres,mxdres),mlpmrs(mxdres)
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension ityp(*),iamino(*)

      toang = 0.52917706d0

      do i=1,iatoms
         ityp(i) = 0
      end do

      call inipmf

      iff = 8

      do i=1,ncalf
         call getpdb(i,ipdb,ihpdb)
         iam = iamino(i)
         if (iam.ge.1.and.iam.le.20) then
            do j=1,mxsym
               if (ipdb(j).gt.0) then
c Backbone atoms 1-4
                   if (j.eq.1) then
                      ityp(ipdb(j)) = 10
                   elseif (j.eq.2) then
                      ityp(ipdb(j)) = 2
                   elseif (j.eq.3) then
                      ityp(ipdb(j)) = 4
                   elseif (j.eq.4) then
                      ityp(ipdb(j)) = 6
                   else
c SideChain atoms
                      do k=1,ipmfl(iam)
                         if (j.eq.lpmfrs(iam,(k-1)*2+1)) then
                            ityp(ipdb(j)) = lpmfrs(iam,(k-1)*2+2)
                            goto 10
                         endif
                      end do

10                    continue
                   endif
               endif
            end do

            do j=1,mxhsym*3
               if (ihpdb(j).gt.0) then
                   ityp(ihpdb(j)) = 1
               endif
            end do

         endif
      end do

      call dfires

      call parptr(153,estat,fdum,idum)

      return
      end

      subroutine dfipar
      implicit double precision (a-h,o-z)
      parameter (mxdres=40)
      parameter (mxrot=81)
      common /dfiwrk/estat(mxdres,mxrot),ndrs,idrs(mxdres),lrot(mxdres),
     &               lpmres(mxdres,mxdres),nlpmrs(mxdres),
     &               kpmres(mxdres,mxdres),mlpmrs(mxdres)

      call parptr(153,estat,fdum,idum)

      return
      end

      subroutine updred(coo,ncalf,icalf)
c
c update the array ipmres holding the residues which are close
c enough for an pmf interaction with the ligand
c since the ligand changes position, this has to be updated
c every once in a while
c
      implicit double precision (a-h,o-z)
      parameter (mxliga=200)
      parameter (numcal=50000)
      common /athlp/ iatoms, mxnat
      common /pmfrk/ pmlig(mxliga),ipmlig(mxliga),npmlig,
     &               ipmres(numcal),npmres
      common /cllmat/ rdum(12),tz(3),tzorg(3),itz
      dimension vec(3)
      dimension coo(3,*),icalf(6,*)

      do i=1,3
         vec(i) = 0.0d0
      end do

      do i=1,npmlig
         k = ipmlig(i)
         do j=1,3
            vec(j) = vec(j) + coo(j,k)
         end do
      end do

c center of ligand atoms

      do j=1,3
         vec(j) = vec(j) / dble(npmlig)
      end do

c estimate radius of ligand

      erad = 0.0d0
      do i=1,npmlig
         d2 = dist2(coo(1,ipmlig(i)),vec)
         if (d2.gt.erad) erad = d2
      end do

      erad = dsqrt(erad) + 19.0d0
      erad = erad * erad
c max radius amino acid estim. 10 Angs

      npmres = 0
      do i=1,ncalf
         d2 = dist2(coo(1,icalf(1,i)),vec)
         if (d2.lt.erad) then
            if (npmres.lt.numcal) then
               npmres = npmres + 1
               ipmres(npmres) = i
            endif
         endif
      end do

c copy ligand translation to tzorg

      do i=1,3
         tzorg(i) = tz(i)
      end do

      print*,'updating interacting residues'

      return
      end

      subroutine dfired(coo,icalf,ncalf,iamino)
c
c for each residue for which dfire score has to be calculated idrs(1..ndrs)
c update the array lpmres holding the residues which are close
c enough for an pmf interaction with the residue of choice
c
      implicit double precision (a-h,o-z)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxres=42)
      parameter (mxdres=40)
      parameter (mxrot=81)
      common /athlp/ iatoms, mxnat
c lpmres holds flexible residues within reach of flexible residue l
c kpmres holds static   residues within reach of flexible residue l
      common /dfiwrk/estat(mxdres,mxrot),ndrs,idrs(mxdres),lrot(mxdres),
     &               lpmres(mxdres,mxdres),nlpmrs(mxdres),
     &               kpmres(mxdres,mxdres),mlpmrs(mxdres)
      common /dfisto/ rres(mxres),
     &                lpdb(mxdres,mxsym),lhpdb(mxdres,mxhsym*3)
      common /doh/ idoh
      data rres /0.0,1.6,3.4,3.0,3.5,4.7,3.5,6.4,3.8,4.6,4.8,7.2,
     &           5.0,5.9,3.0,8.3,5.8,6.2,6.9,7.9,22*0.0/
      dimension ipdb(mxsym),ihpdb(mxhsym*3),vec(3),dmaxa(mxdres)
      dimension coo(3,*),icalf(6,*),iamino(*)

      toang = 0.52917706d0
      iprtms = 0

      do l=1,ndrs

         jin = idrs(l)

         call getpdb(jin,ipdb,ihpdb)

         do i=1,mxsym
            lpdb(l,i) = ipdb(i)
         end do

         do i=1,mxhsym*3
            lhpdb(l,i) = ihpdb(i)
         end do

c center of residue

         do j=1,3
            vec(j) = coo(j,icalf(1,jin))
         end do

c max radius amino acid estim. 10 Angs, max dist two res = 20, * 2 for a.u.

c         erad = 20.0d0*2.0d0


         nlpmrs(l) = 0
         mlpmrs(l) = 0
         dmax = 0.0d0
         nmax = 0
         do i=1,ncalf
            if (i.ne.jin) then

c               d2 = dist2(coo(1,icalf(1,i)),vec)
c               erad = (rres(iamino(jin)) + rres(iamino(i))) / toang
c               erad = erad * erad
c
c               if (d2.lt.erad) then

               d2 = dsqrt(dist2(coo(1,icalf(1,i)),vec))
               erad = (rres(iamino(jin)) + rres(iamino(i)))
               
               if (d2*toang - erad .lt. 14.5d0) then

                  nostat = 0
                  do j=1,ndrs
                     if (i.eq.idrs(j).and.jin.ne.idrs(j)) nostat = j
                  end do

                  if (nostat.ne.0) then
c flexible residues
                     if (nlpmrs(l).lt.mxdres) then
                        nlpmrs(l) = nlpmrs(l) + 1
                        lpmres(l,nlpmrs(l)) = nostat
                     endif
                  else
c non flexible residues
                     if (mlpmrs(l).lt.mxdres) then
                        mlpmrs(l) = mlpmrs(l) + 1
                        kpmres(l,mlpmrs(l)) = i
                        dmaxa(mlpmrs(l)) = d2
                        if (d2.gt.dmax) dmax = d2
                     else

c dmaxa list is not sorted, dmax indicates element with highest
c distance, this one will be replaced with a new smaller distance
c in the array to hold static residues within contact

                        if (d2.lt.dmax) then
                           dmaxt = d2
                           do m=1,mxdres
                              if (dmaxa(m).eq.dmax) then
                                 kpmres(l,m) = i
                                 dmaxa(m) = d2
                              elseif (dmaxa(m).gt.dmaxt) then
                                 dmaxt = dmaxa(m)
                              endif
                           end do
                           dmax = dmaxt
                        endif
                        if (iprtms.eq.0) then
                        print*,'Interactionlist array to small,'// 
     &                         ' keeping only closest residues'
                            iprtms = 1
                        endif
                     endif
                  endif

               endif
            endif
         end do

      end do

      return
      end

      subroutine dfirot(l)
c
c use this function to calculate estat:
c
c rotamer with all non flexible sidechains and rotamer with backbone part
c flexible sidechains and internal rotamer energy (fixed or self pmf ?)
c
c rotamer must be specified via lrot(l) and new geometry must have been
c generated
c
      implicit double precision (a-h,o-z)
      parameter (mxdres=40)
      parameter (mxrot=81)
      common /dfiwrk/estat(mxdres,mxrot),ndrs,idrs(mxdres),lrot(mxdres),
     &               lpmres(mxdres,mxdres),nlpmrs(mxdres),
     &               kpmres(mxdres,mxdres),mlpmrs(mxdres)


      jin = idrs(l)
      irot = lrot(l)

      estat(l,irot) = 0.0d0

c contrib backbone flexible side-chains

c      print*,'irot=',irot
      do i=1,nlpmrs(l)
         k = idrs(lpmres(l,i))
         call twodfib(twod,jin,k)
         estat(l,irot) = estat(l,irot) + twod
      end do

c contrib backbone + side chain static side-chains

      do i=1,mlpmrs(l)
         k = kpmres(l,i)
	 call twodfib(e1,jin,k)
	 call twodfi(e2,jin,k)
         estat(l,irot) = estat(l,irot) + e1 + e2
c         print*,'est+ ',jin,',',k,' ',e1,' ',e2
      end do

c self contribution

      call onedfi(e3,jin)
      estat(l,irot) = estat(l,irot) + e3
c      print*,'est+self ',jin,',',e3

      return
      end

      subroutine dfiflx(etot)
c
c use this function to calculate all interactions between flexible
c residues (except with backbone flex residues)
c
      implicit double precision (a-h,o-z)
      parameter (mxdres=40)
      parameter (mxrot=81)
      common /dfiwrk/estat(mxdres,mxrot),ndrs,idrs(mxdres),lrot(mxdres),
     &               lpmres(mxdres,mxdres),nlpmrs(mxdres),
     &               kpmres(mxdres,mxdres),mlpmrs(mxdres)


      etot = 0.0d0

      do i=1,ndrs
         do j=i+1,ndrs
	    call twodfi(e2,idrs(i),idrs(j))
            etot = etot + e2
c         print*,'est+ ',i,'=',idrs(i),j,idrs(j),' e=',e2
         end do
      end do


      return
      end

      double precision function pmfsco(ptyp,ltyp,d)
c
c calculates pmfscore of a protein atom - ligand atom couple
c
      implicit double precision (a-h,o-z)
      integer*2 ptyp,ltyp
      integer ptype
      real pmf,dfire
      common /pmfpar/ ipmfrs(20,20),lpmfrs(20,20),ipmfl(20),
     &                pmf(16,26,60),dfire(13,13,20)

      pmfsco = 0.0d0

      ptype = int(ptyp) 
      ltype = int(ltyp) 

      cutoff = 9.0d0
      if ((ptype.ge.1.and.ptype.le.4).and.ltype.eq.27) cutoff = 6.0d0
c actually Cl (=27) has no data

      if (d.gt.cutoff) return
      if (ptype.lt.1.or.ptype.gt.16) return
      if (ltype.lt.1.or.ltype.gt.26) return


      it2 = int(d/0.2d0)
      it1 = it2 - 1
      it3 = it2 + 1

      if (it1.lt.0) then
         it1 = 0
      elseif (it1.gt.59) then
         it1 = 59
      endif

      if (it2.lt.0) then
         it2 = 0
      elseif (it2.gt.59) then
         it2 = 59
      endif

      if (it3.lt.0) then
         it3 = 0
      elseif (it3.gt.59) then
         it3 = 59
      endif

      it1 = it1 + 1
      it2 = it2 + 1
      it3 = it3 + 1

c      pmfsco = pmf(ptype,ltype,it2)

c smooth version

      pmfsco = 0.0d0
      pmfsco = pmfsco + 0.2d0*(dble(pmf(ptype,ltype,it1))+3.0d0)
      pmfsco = pmfsco + 0.6d0*(dble(pmf(ptype,ltype,it2))+3.0d0)
      pmfsco = pmfsco + 0.2d0*(dble(pmf(ptype,ltype,it3))+3.0d0)

      return
      end

      double precision function dfisco(ptyp,ltyp,d)
c
c calculates dfire score of a protein atom - protein atom couple
c
      implicit double precision (a-h,o-z)
      integer*2 ptyp,ltyp
      integer ptype
      real pmf,dfire
      common /pmfpar/ ipmfrs(20,20),lpmfrs(20,20),ipmfl(20),
     &                pmf(16,26,60),dfire(13,13,20)

      dfisco = 0.0d0

      ptype = int(ptyp) 
      ltype = int(ltyp) 

c d aangeleverd in angstrom ?

      cutoff = 15.0d0

      if (d.gt.cutoff) return
      if (ptype.lt.1.or.ptype.gt.13) return
      if (ltype.lt.1.or.ltype.gt.13) return


      if (d.lt.2.0d0) then
c 1
c 0-2
         it = 1
      elseif (d.lt.8.0d0) then
c  2     3      4     5     6     7     8     9    10     11   12    13
c 2-2.5,2.5-3,3-3.5,3.5-4,4-4.5,4.5-5,5-5.5,5.5-6,6-6.5,6.5-7,7-7.5,7.5-8
         it = int((d-2.0d0)/0.5d0) + 2
      elseif (d.lt.15.d0) then
c  14  15    16    17    18    19    20
c 8-9,9-10,10-11,11-12,12-13,13-14,14-15
         it = int(d-8.0d0) + 14
      else
         return
      endif

      dfisco = dble(dfire(ptype,ltype,it))
      
      return
      it2 = int(d/0.5d0)
      it1 = it2 - 1
      it3 = it2 + 1

      if (it1.lt.0) then
         it1 = 0
      elseif (it1.gt.19) then
         it1 = 19
      endif

      if (it2.lt.0) then
         it2 = 0
      elseif (it2.gt.19) then
         it2 = 19
      endif

      if (it3.lt.0) then
         it3 = 0
      elseif (it3.gt.19) then
         it3 = 19
      endif

      it1 = it1 + 1
      it2 = it2 + 1
      it3 = it3 + 1

c      dfisco = dfire(ptype,ltype,it2)

c smooth version

      dfisco = 0.0d0
      dfisco = dfisco + 0.2d0*(dble(dfire(ptype,ltype,it1)))
      dfisco = dfisco + 0.6d0*(dble(dfire(ptype,ltype,it2)))
      dfisco = dfisco + 0.2d0*(dble(dfire(ptype,ltype,it3)))

      return
      end

      subroutine totpmd(totpmf,
     &                coo,ianz,iaton,iatclr,ityp)
      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      parameter (mxchai=50)
      parameter (mxheta=150)
      parameter (mxres=42)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxliga=200)
      parameter (mxewin=100)
      parameter (mxpmfl=5)
      common /pmfrk/ pmlig(mxliga),ipmlig(mxliga),npmlig,
     &               ipmres(numcal),npmres
      integer*2 ityp,iltp,iptp
      common /pmfres/ ewin(mxewin),npmfs
      common /pmflvl/ ipmfm,ipmfh,pmflev(mxpmfl),levcol(mxpmfl)
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension coo(3,*),ianz(*),iaton(*),iatclr(*),ityp(*)


      toang = 0.52917706d0

      do i=1,npmlig
         pmlig(i) = 0.0d0
         ilat = ipmlig(i)
         iatclr(ilat) = 8
      end do

      pmfmn = 1.0d10
      pmfmx = -1.0d10
      dstmn = 1.0d10
      dstmx = -1.0d10
      do i=1,npmres
         call getpdb(ipmres(i),ipdb,ihpdb)
         do j=1,mxsym
            if (ipdb(j).ne.0) then
               ipat = ipdb(j)
               iptp = ityp(ipat)
               iatclr(ipat) = 9
               iaton(ipat) = 1
               do k=1,npmlig
                  ilat = ipmlig(k)
                  iaton(ilat) = 1
                  if (.not.(ipmfh.eq.0.and.ianz(ilat).eq.1)) then
                     iltp = ityp(ilat)
                     d2 = dist2(coo(1,ipat),coo(1,ilat))
                     d2 = dsqrt(d2)*toang
                     if (d2.gt.dstmx) dstmx = d2
                     if (d2.lt.dstmn) dstmn = d2
                     pmf = pmfsco(iptp,iltp,d2)
                     if (pmf.gt.pmfmx) pmfmx = pmf
                     if (pmf.lt.pmfmn) pmfmn = pmf
                     if (pmf.gt.2.0d0) then
                        iatclr(ipat) = 4
                        iatclr(ilat) = 4
                     endif
                     pmlig(k) = pmlig(k) + pmf
                  endif
               end do
            endif
         end do

         if (ipmfh.eq.1) then

            do j=1,mxhsym*3
               if (ihpdb(j).ne.0) then
                  ipat = ihpdb(j)
                  iatclr(ipat) = 9
                  iaton(ipat) = 1
                  iptp = ityp(ipat)
                  do k=1,npmlig
                     ilat = ipmlig(k)
                     iaton(ilat) = 1
                     iltp = ityp(ilat)
                     d2 = dist2(coo(1,ipat),coo(1,ilat))
                     d2 = dsqrt(d2)*toang
                     pmf = pmfsco(iptp,iltp,d2)
                     if (pmf.gt.pmfmx) pmfmx = pmf
                     if (pmf.lt.pmfmn) pmfmn = pmf
                     if (pmf.gt.2.0d0) then
                        iatclr(ipat) = 4
                        iatclr(ilat) = 4
                     endif
                     pmlig(k) = pmlig(k) + pmf
                  end do
               endif
            end do
         endif

      end do


      totpmf = 0.0d0
      pmfmin = 1.0d10
      pmfmax = -1.0d10
      nligef = 0
      do i=1,npmlig
         if (pmlig(i).gt.pmfmax) pmfmax = pmlig(i)
         if (pmlig(i).lt.pmfmin) pmfmin = pmlig(i)
         totpmf = totpmf + pmlig(i)
         if (pmlig(i).ne.0.0d0) nligef = nligef + 1
      end do

      if (ipmfm.eq.1) then
         pmfavg = totpmf / dble(nligef)
         pmflev(1) = pmfmax
         pmflev(3) = pmfavg
         pmflev(5) = pmfmin
         pmflev(2) = (pmflev(1) + pmflev(3))/ 2.0d0
         pmflev(4) = (pmflev(5) + pmflev(3))/ 2.0d0
         do i=1,npmlig
            pm = pmlig(i)
            ilcol = 1
            if (pm.le.pmflev(1).and.pm.gt.pmflev(2)) then
                ilcol = levcol(1)
            elseif (pm.le.pmflev(2).and.pm.gt.pmflev(3)) then
                ilcol = levcol(2)
            elseif (pm.le.pmflev(3).and.pm.gt.pmflev(4)) then
                ilcol = levcol(3)
            elseif (pm.le.pmflev(4).and.pm.ge.pmflev(5)) then
                ilcol = levcol(4)
            endif
            if (.not.(ianz(ipmlig(i)).eq.1.and.ipmfh.eq.0)) then
               iatclr(ipmlig(i)) = ilcol
            end if
         end do
      endif

      if (npmfs.lt.mxewin) then
         npmfs = npmfs + 1
         ewin(npmfs) = totpmf
      else
         do i=1,mxewin-1
            ewin(i) = ewin(i+1)
         end do
         ewin(mxewin) = totpmf
      endif
      
      return
      end

      double precision function totdfi(idum)
      implicit double precision (a-h,o-z)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxdres=40)
      parameter (mxrot=81)
      common /dfiwrk/estat(mxdres,mxrot),ndrs,idrs(mxdres),lrot(mxdres),
     &               lpmres(mxdres,mxdres),nlpmrs(mxdres),
     &               kpmres(mxdres,mxdres),mlpmrs(mxdres)

      toang = 0.52917706d0

      totdfi = 0.0d0

c totdfi = static part + selfpart + rotres1*rotres2

      do i=1,ndrs

c static part + self part

         totdfi = totdfi + estat(i,lrot(i))

c rotres1*rotres2

         do j=1,nlpmrs(i)
            call twodfi(twod,idrs(i),idrs(lpmres(i,j)))
            totdfi = totdfi + twod
         end do

      end do

      return
      end

      subroutine twodfd(twodfi,ires1,ires2,coo,ityp)
c
c pmf score of the two sidechains without the backbone scores
c
      implicit double precision (a-h,o-z)
      parameter (mxres=42)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxdres=40)
      integer*2 ityp,iltp,iptp
      common /dfisto/ rres(mxres),
     &                lpdb(mxdres,mxsym),lhpdb(mxdres,mxhsym*3)
      common /doh/ idoh
      dimension ipdb(mxsym),ipdb2(mxsym),
     &          ihpdb(mxhsym*3),ihpdb2(mxhsym*3)
      dimension coo(3,*),ityp(*)

      twodfi = 0.0d0

      toang = 0.52917706d0

      ifl1 = iflex(ires1)
      ifl2 = iflex(ires2)

      if (ifl1.eq.0) call getpdb(ires1,ipdb,ihpdb)
      if (ifl2.eq.0) call getpdb(ires2,ipdb2,ihpdb2)

      do j=5,mxsym
         
         if (ifl1.eq.0) then
             ipat = ipdb(j)
         else 
             ipat = lpdb(ifl1,j)
         endif

         if (ipat.ne.0) then

            iptp = ityp(ipat)

            do k=5,mxsym

               if (ifl2.eq.0) then
                  ilat = ipdb2(k)
               else 
                  ilat = lpdb(ifl2,k)
               endif

               if (ilat.ne.0) then
                  iltp = ityp(ilat)
                  d2 = dist2(coo(1,ipat),coo(1,ilat))
                  d2 = dsqrt(d2)*toang
                  dfi = dfisco(iptp,iltp,d2)
c                  if (k.eq.10.and.j.eq.24) 
c     &   print*,'d2=',d2,' dfi=',dfi,' iptp=',iptp,' iltp=',iltp
                  twodfi = twodfi + dfi
               endif
            end do
         endif
      end do


c also do H....H contacts ?

      if (idoh.eq.1) then

         do j=4,mxhsym*3
         
            if (ifl1.eq.0) then
                ipat = ihpdb(j)
            else 
                ipat = lhpdb(ifl1,j)
            endif

            if (ipat.ne.0) then

               iptp = ityp(ipat)

               do k=4,mxhsym*3

                  if (ifl2.eq.0) then
                     ilat = ihpdb2(k)
                  else 
                     ilat = lhpdb(ifl2,k)
                  endif

                  if (ilat.ne.0) then
                     iltp = ityp(ilat)
                     d2 = dist2(coo(1,ipat),coo(1,ilat))
                     d2 = dsqrt(d2)*toang
                     dfi = dfisco(iptp,iltp,d2)
                     twodfi = twodfi + dfi
                  endif
               end do
            endif
         end do

      endif

c      print*,'twodfi  ',ires1,' ',ires2,' =',twodfi

      return
      end

      subroutine twodfid(twodfib,ires1,ires2,coo,ityp)
c
c pmf score of sidechain ires1 with the backbone of ires2
c
      implicit double precision (a-h,o-z)
      parameter (mxres=42)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxdres=40)
      integer*2 ityp,iltp,iptp
      common /dfisto/ rres(mxres),
     &                lpdb(mxdres,mxsym),lhpdb(mxdres,mxhsym*3)
      common /doh/ idoh
      dimension ipdb(mxsym),ipdb2(mxsym),
     &          ihpdb(mxhsym*3),ihpdb2(mxhsym*3)
      dimension coo(3,*),ityp(*)

      twodfib = 0.0d0

      toang = 0.52917706d0

      ifl1 = iflex(ires1)
      ifl2 = iflex(ires2)

      if (ifl1.eq.0) call getpdb(ires1,ipdb,ihpdb)
      if (ifl2.eq.0) call getpdb(ires2,ipdb2,ihpdb2)

      do j=4,mxsym

         if (ifl1.eq.0) then
            ipat = ipdb(j)
         else 
            ipat = lpdb(ifl1,j)
         endif

         if (ipat.ne.0) then

            iptp = ityp(ipat)

            do k=1,4

               if (ifl2.eq.0) then
                  ilat = ipdb2(k)
               else 
                  ilat = lpdb(ifl2,k)
               endif

               if (ilat.ne.0) then
                  iltp = ityp(ilat)
                  d2 = dist2(coo(1,ipat),coo(1,ilat))
                  d2 = dsqrt(d2)*toang
                  dfi = dfisco(iptp,iltp,d2)
c                  if (k.eq.2.and.j.eq.13) 
c     &   print*,'d2=',d2,' dfi=',dfi,' iptp=',iptp,' iltp=',iltp,
c     &          ' ires1=',ires1,' ires2=',ires2
                  twodfib = twodfib + dfi
               endif
            end do
         endif
      end do

c also do H....H contacts ?

      if (idoh.eq.1) then

         do j=5,mxhsym*3

            if (ifl1.eq.0) then
               ipat = ihpdb(j)
            else 
               ipat = lhpdb(ifl1,j)
            endif

            if (ipat.ne.0) then

               iptp = ityp(ipat)

               do k=1,4

                  if (ifl2.eq.0) then
                     ilat = ihpdb2(k)
                  else 
                     ilat = lhpdb(ifl2,k)
                  endif

                  if (ilat.ne.0) then
                     iltp = ityp(ilat)
                     d2 = dist2(coo(1,ipat),coo(1,ilat))
                     d2 = dsqrt(d2)*toang
                     dfi = dfisco(iptp,iltp,d2)
                     twodfib = twodfib + dfi
                  endif

               end do

            endif
         end do

      endif

c      print*,'twodfib ',ires1,' ',ires2,' =',twodfib

      return
      end

      integer function iflex(ires)
      implicit double precision (a-h,o-z)
      parameter (mxdres=40)
      parameter (mxrot=81)
      common /dfiwrk/estat(mxdres,mxrot),ndrs,idrs(mxdres),lrot(mxdres),
     &               lpmres(mxdres,mxdres),nlpmrs(mxdres),
     &               kpmres(mxdres,mxdres),mlpmrs(mxdres)

      iflex = 0

      do i=1,ndrs
         if (ires.eq.idrs(i)) iflex = i
      end do

      return
      end

      subroutine onedfd(onedfi,ires1,coo,ityp)
c self contribution of residue with it self
      implicit double precision (a-h,o-z)
      parameter (mxres=42)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxdres=40)
      integer*2 ityp,iltp,iptp
      common /dfisto/ rres(mxres),
     &                lpdb(mxdres,mxsym),lhpdb(mxdres,mxhsym*3)
      logical o34
      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension coo(3,*),ityp(*)

      onedfi = 0.0d0

      toang = 0.52917706d0

      ifl1 = iflex(ires1)
      if (ifl1.eq.0) call getpdb(ires1,ipdb,ihpdb)

      do j=6,mxsym

         if (ifl1.eq.0) then
            ipat = ipdb(j)
         else 
            ipat = lpdb(ifl1,j)
         endif

         if (ipat.ne.0) then
            iptp = ityp(ipat)
            do k=1,4

               if (ifl1.eq.0) then
                  ilat = ipdb(k)
               else 
                  ilat = lpdb(ifl1,k)
               endif

               if (ilat.ne.0) then
                  if (o34(k,j)) then
                     iltp = ityp(ilat)
                     d2 = dist2(coo(1,ipat),coo(1,ilat))
                     d2 = dsqrt(d2)*toang
                     dfi = dfisco(iptp,iltp,d2)
                     onedfi = onedfi + dfi
                  endif
               endif
            end do
         endif
      end do

      return
      end

      logical function o34(k,j)
      implicit double precision (a-h,o-z)

      o34 = .true.

      if (k.ne.2) return

      if (j.eq.6)  goto 100
      if (j.eq.7)  goto 100
      if (j.eq.8)  goto 100
      if (j.eq.31) goto 100
      if (j.eq.32) goto 100
      if (j.eq.37) goto 100

      return

100   o34 = .false.
      return
      end 

      subroutine pmfind(iatm,coo,ianz,ityp)
      implicit double precision (a-h,o-z)
      parameter (numcal=50000)
      parameter (mxsym=103)
      parameter (mxhsym=64)
      parameter (mxliga=200)
      parameter (mxewin=100)
      parameter (mxmx=20)
      parameter (mxppmf=16)
      parameter (mxlpmf=26)
      parameter (mxmol2=41)
      parameter (mxmm3=164)
      parameter (mxmsf=235)
      parameter (mxamb=1590)
      parameter (mxamo=201)
      parameter (mxchtp=136)
      parameter (mxpmfl=5)
      common /pmfrk/ pmlig(mxliga),ipmlig(mxliga),npmlig,
     &               ipmres(numcal),npmres
      integer*2 ityp,iltp,iptp
      common /pmfres/ ewin(mxewin),npmfs
      common /pmfint/ pmfmx(mxmx),pmfmn(mxmx),pmflga,
     &                ipmfmx(mxmx),ipmfmn(mxmx),npmfmx
      character*2 ppmf, lpmf
      character*5 mol2
      character*19 mm3
      character*20 chmtnk
      character*4 chmsf
      character*20 ambstr
      character*20 amostr
      common /ftypes/ihasl(11),mol2(mxmol2),mm3(mxmm3),chmtnk(mxchtp),
     &               chmsf(mxmsf),ambstr(mxamb),amostr(mxamo),
     &               ppmf(mxppmf),lpmf(mxlpmf)
      common /pmflvl/ ipmfm,ipmfh,pmflev(mxpmfl),levcol(mxpmfl)

      dimension ipdb(mxsym),ihpdb(mxhsym*3)
      dimension dstmx(mxmx),dstmn(mxmx)
      dimension coo(3,*),ianz(*),ityp(*)

      toang = 0.52917706d0

      npmfmx = 10
      pmftot = 0.0d0
      ilig = 0
      do i=1,npmlig
         if (iatm.eq.ipmlig(i)) ilig = i
      end do

      if (ilig.eq.0) then
         print*,'not a pmf ligand atom'
         return
      endif

      do i=1,npmfmx
         pmfmn(i) = 1.0d10
         pmfmx(i) = -1.0d10
         dstmn(i) = 1.0d10
         dstmx(i) = -1.0d10
      end do 
      ddstmn = 1.0d10
      ddstmx = -1.0d10

      do i=1,npmres
         call getpdb(ipmres(i),ipdb,ihpdb)
         do j=1,mxsym
            if (ipdb(j).ne.0) then
               ipat = ipdb(j)
               iptp = ityp(ipat)
               ilat = ipmlig(ilig)
               if (.not.(ipmfh.eq.0.and.ianz(ilat).eq.1)) then
                  iltp = ityp(ilat)
                  d2 = dist2(coo(1,ipat),coo(1,ilat))
                  d2 = dsqrt(d2)*toang
                  if (d2.gt.ddstmx) ddstmx = d2
                  if (d2.lt.ddstmn) ddstmn = d2
                  pmf = pmfsco(iptp,iltp,d2)
                  do k=1,npmfmx
                     if (pmf.gt.pmfmx(k)) then
                        do l=npmfmx,k,-1
                           if (l-1.gt.0) then
                              pmfmx(l) = pmfmx(l-1)
                              ipmfmx(l) = ipmfmx(l-1)
                              dstmx(l) = dstmx(l-1)
                           endif
                        end do
                        pmfmx(k) = pmf
                        ipmfmx(k) = ipat
                        dstmx(k) = d2
                        goto 10
                     endif
                     if (pmf.lt.pmfmn(k)) then
                        do l=npmfmx,k,-1
                           if (l-1.gt.0) then
                              pmfmn(l) = pmfmn(l-1)
                              ipmfmn(l) = ipmfmn(l-1)
                              dstmn(l) = dstmn(l-1)
                           endif
                        end do
                        pmfmn(k) = pmf
                        ipmfmn(k) = ipat
                        dstmn(k) = d2
                        goto 10
                     endif
                  end do
10                continue
                  pmftot = pmftot + pmf
               endif
            endif
         end do

         if (ipmfh.eq.1) then
            do j=1,mxhsym*3
               if (ihpdb(j).ne.0) then
                  ipat = ihpdb(j)
                  iptp = ityp(ipat)
                  ilat = ipmlig(ilig)
                  iltp = ityp(ilat)
                  d2 = dist2(coo(1,ipat),coo(1,ilat))
                  d2 = dsqrt(d2)*toang
                  if (d2.gt.ddstmx) ddstmx = d2
                  if (d2.lt.ddstmn) ddstmn = d2
                  pmf = pmfsco(iptp,iltp,d2)
                  do k=1,npmfmx
                     if (pmf.gt.pmfmx(k)) then
                        do l=npmfmx,k,-1
                           if (l-1.gt.0) then
                              pmfmx(l) = pmfmx(l-1)
                              ipmfmx(l) = ipmfmx(l-1)
                              dstmx(l) = dstmx(l-1)
                           endif
                        end do
                        pmfmx(k) = pmf
                        ipmfmx(k) = ipat
                        dstmx(k) = d2
                        goto 20
                     endif
                     if (pmf.lt.pmfmn(k)) then
                        do l=npmfmx,k,-1
                           if (l-1.gt.0) then
                              pmfmn(l) = pmfmn(l-1)
                              ipmfmn(l) = ipmfmn(l-1)
                              dstmn(l) = dstmn(l-1)
                           endif
                        end do
                        pmfmn(k) = pmf
                        ipmfmn(k) = ipat
                        dstmn(k) = d2
                        goto 20
                     endif
                  end do
20                continue
                  pmftot = pmftot + pmf
               endif
            end do
         endif

      end do

      pmflga = pmftot
c      do i=1,npmfmx
c         print*,'pmfmn=',pmfmn(i),' pmfmx=',pmfmx(i)
c         print*,'ipmfmn=',ipmfmn(i),' ipmfmx=',ipmfmx(i)
c         print*,'dstmn=',dstmn(i),' dstmx=',dstmx(i)
c      end do
c      print*,'ddstmn=',ddstmn,' ddstmx=',ddstmx

      call parptr(117,pmfmn,pmfmx,npmfmx)
      call parptr(118,fdum,fdum,ipmfmn)
      call parptr(119,fdum,fdum,ipmfmx)
      call parptr(120,pmflga,fdum,idum)
c to xwindow and show the pmfinf window
      call cpmf

      return
      end

      subroutine ipmtyd(ipmtyp,iat,ian,idochg,ianz,iconn,qat)
      implicit double precision (a-h,p-z),integer (i-n),logical (o)
      parameter (mxcon=10)
      dimension icnn(mxcon)
      dimension ianz(*),iconn(mxcon+1,*),qat(*)

      ibnds = 0
      io = 0
      in = 0
      ic = 0
      ih = 0
      ihet = 0
      ico = 0
      icn = 0
      ipmtyp = -1
      do i=1,iconn(1,iat)
         if (iconn(i+1,iat).gt.0) then
            ibnds = ibnds + 1
            icnn(ibnds) = iconn(i+1,iat)
            ia = ianz(icnn(ibnds))
            if (ia.eq.1) ih = ih + 1
            if (ia.eq.6) ic = ic + 1
            if (ia.eq.7) then
               in = in + 1
               if (qat(icnn(ibnds)).ge.0.01d0) icn = icn + 1
            endif
            if (ia.eq.8) then
               io = io + 1
               if (qat(icnn(ibnds)).le.-0.01d0) ico = ico + 1
            endif
c HET are N,O,F, Si,P,S,Cl, Br, I
            if ((ia.ge.7.and.ia.le.9).or.(ia.ge.14.and.ia.le.17).or.
     &         ia.eq.35.or.ia.eq.53) ihet = ihet + 1
         endif
      end do

      call ispn(ispt,iat,irng,idochg,1)
      
      if (ian.eq.1) then

c     Hydrogen -> HL

          ipmtyp = 25

      elseif (ian.eq.6) then

c     Carbon

          if (ispt.eq.2) then
c            C.1 -> C0
             ipmtyp = 34
          elseif (ispt.eq.7) then
c            C.cat -> CN
             ipmtyp = 8
          elseif (idochg.eq.1.and.ico.gt.0) then
c            C attached to neg O  -> CO
             ipmtyp = 7
          elseif (idochg.eq.1.and.icn.gt.0) then
c            C attached to pos N  -> CN
             ipmtyp = 8
          else
             if (irng.gt.0.or.ispt.eq.6) then
c            aromatic ring -> cF or cP
                 if (ihet.eq.0) then
                    ipmtyp = 3
                 else
                    ipmtyp = 4
                 endif
             elseif (ispt.eq.4) then
c            C.3 -> CF or CP
                 if (ihet.eq.0) then
                    ipmtyp = 1
                 else
                    ipmtyp = 2
                 endif
             elseif (ispt.eq.3) then
c            C.2 -> C3 or CW
                 if (ihet.eq.0) then
                    ipmtyp = 5
                 else
                    ipmtyp = 6
                 endif
             endif
          endif
      elseif (ian.eq.7) then

c Nitrogen

          if (ispt.eq.6) then
c            N.ar -> NR
             ipmtyp = 13
          elseif (ispt.eq.2) then
c            N.1 -> N0
             ipmtyp = 14
          elseif (qat(iat).ge.0.01d0) then
c            positively charged nitrogen -> NC
             ipmtyp = 9
          else

             if (irng.eq.-1) then

c                non aromatic ring -> NR, NS, NP

                 if (ih.ne.0) then
                    ipmtyp = 13
                 elseif (ihet.gt.0) then
                    ipmtyp = 15
                 else
                    ipmtyp = 10
                 endif

             elseif (irng.gt.0) then

c                aromatic ring -> NR

                 ipmtyp = 13
             else

c                no ring -> ND, NS, NP, NA

                 if (ih.ne.0) then
                    ipmtyp = 12
                 elseif (ihet.gt.0) then
                    ipmtyp = 15
                 elseif (ispt.eq.9.or.ispt.eq.8) then
                    ipmtyp = 10
                 else
                    ipmtyp = 11
                 endif

             endif
          endif

      elseif (ian.eq.8) then

c Oxygen

          if (ispt.eq.4) then

c         O.3
             if (irng.gt.0) then
c               O in aromatic ring -> OR
                ipmtyp = 19
             elseif (ih.gt.0) then
c               O attached to H -> OD
                ipmtyp = 21
             elseif (ihet.gt.0) then
c               O attached to non C, H -> OS
                ipmtyp = 20
             else
c               other aromatic O -> OE
                ipmtyp = 18
             endif

          elseif (ispt.eq.3) then

c         O.2 -> OA
             ipmtyp = 17

          elseif (ispt.eq.10) then

c         O.co2 -> OC
             ipmtyp = 16

          endif

      elseif (ian.eq.16) then

c Sulphur

          if (ispt.eq.4) then
c         S.3  -> SD, SA
             if (ibnds.ge.1) then
                ipmtyp = 24
             else
                ipmtyp = 23
             endif
          elseif (ispt.eq.3.or.ispt.eq.13.or.ispt.eq.14) then

c         S.2, S.o, S.o2  -> SD
             ipmtyp = 24
          endif

      elseif (ian.eq.15) then

c Phosphor -> P

          ipmtyp = 22
      elseif (ian.eq.9) then
c Fluor -> F
          ipmtyp = 26
      endif

c types not assigned: Zn, Mn, Mg, Fe, V, Cl, Br, C0 there is no data

      return
      end

c Protein types:
c
c 1  CF
c 2  CP
c 3  cF
c 4  cP
c 5  CO
c 6  CN
c 7  NC
c 8  ND
c 9  NR
c 10 OC
c 11 OA
c 12 OD
c 13 OW
c 14 SA
c 15 SD
c 16 HH
c
c Ligand types:
c
c CF	1
c CP	2
c cF	3
c cP	4
c C3	5
c CW	6
c CO	7
c CN	8
c NC	9
c NP	10
c NA	11
c ND	12
c NR	13
c N0	14 is actually also no data, but will leave it in
c NS	15
c OC	16
c OA	17
c OE	18
c OR	19
c OS	20
c OD	21
c P 	22
c SA	23
c SD 	24
c HL 	25
c Zn	------ no data
c Cl	------ no data
c Mn	------ no data
c Mg	------ no data
c F 	26
c Fe	------ no data
c Br	------ no data
c V 	------ no data
c C0	------ no data
c
      subroutine inipmf
      implicit double precision (a-h,o-z)
      real pmf,dfire
      common /pmfpar/ ipmfrs(20,20),lpmfrs(20,20),ipmfl(20),
     &                pmf(16,26,60),dfire(13,13,20)
c      dimension ipmfrs(nres,10*2)
c     for all
c     1,6, 3,2, 4,11, 2,2
      data (ipmfl(i),i=1,20) /
     &  0,1,2,2,3,4,3,4,4,4,4,5,5,5,3,7,6,7,8,10 /
c     gly
      data (ipmfrs(1,i),i=1,20)
     &    /0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     ala
      data (ipmfrs(2,i),i=1,20)
     &    /5,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     ser
      data (ipmfrs(3,i),i=1,20)
     &    /5,2, 31,12, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     cys
      data (ipmfrs(4,i),i=1,20)
     &    /5,2, 37,14, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     thr
      data (ipmfrs(5,i),i=1,20)
     &    /5,2, 32,12, 8,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     ile
      data (ipmfrs(6,i),i=1,20)
     &    /5,1, 8,1, 7,1, 10,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     val
      data (ipmfrs(7,i),i=1,20)
     &    /5,1, 7,1, 8,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     met
      data (ipmfrs(8,i),i=1,20)
     &    /5,1, 6,2, 36,14, 12,2, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     asp
      data (ipmfrs(9,i),i=1,20)
     &    /5,1, 6,5, 29,10, 30,10, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     asn
      data (ipmfrs(10,i),i=1,20)
     &    /5,1, 6,2, 29,11, 21,8, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     leu
      data (ipmfrs(11,i),i=1,20)
     &    /5,1, 6,1, 10,1, 11,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     lys
      data (ipmfrs(12,i),i=1,20)
     &    /5,1, 6,1, 9,1, 12,6, 27,7, 0,0, 0,0, 0,0, 0,0, 0,0/
c     glu
      data (ipmfrs(13,i),i=1,20)
     &    /5,1, 6,1, 9,5, 34,10, 35,10, 0,0, 0,0, 0,0, 0,0, 0,0/
c     gln
      data (ipmfrs(14,i),i=1,20)
     &    /5,1, 6,1, 9,2, 34,11, 24,8, 0,0, 0,0, 0,0, 0,0, 0,0/
c     pro
      data (ipmfrs(15,i),i=1,20)
     &    /5,1, 6,1, 9,2, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     arg
      data (ipmfrs(16,i),i=1,20)
     &    /5,1, 6,1, 9,6, 22,7, 17,6, 25,7, 26,7, 0,0, 0,0, 0,0/
c     his
      data (ipmfrs(17,i),i=1,20)
     &    /5,1, 6,4, 20,9, 11,4, 13,4, 24,9, 0,0, 0,0, 0,0, 0,0/
c     phe
      data (ipmfrs(18,i),i=1,20)
     &    /5,1, 6,3, 10,3, 11,3, 13,3, 14,3, 17,3, 0,0, 0,0, 0,0/
c     tyr
      data (ipmfrs(19,i),i=1,20)
     &    /5,1, 6,3, 10,3, 11,3, 13,3, 14,3, 17,4, 33,12, 0,0, 0,0/
c     trp
      data (ipmfrs(20,i),i=1,20)
     &    /5,1, 6,3, 11,3, 14,4, 15,3, 10,4, 23,8, 18,3, 19,3, 16,3/

c dfire res-res types

c  1 H
c  2 C3
c  3 Car
c  4 C2
c  5 Ccat
c  6 O2
c  7 O3
c  8 Oco2
c  9 Npl3
c 10 Nam
c 11 N4
c 12 N2
c 13 S3

c check met aleks atoom types op gln (N -> N.am) 
c                                his (C -> C.2)
c                                trp (C 5 meb ring  -> C.2)

c     gly
      data (lpmfrs(1,i),i=1,20)
     &    /0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     ala
      data (lpmfrs(2,i),i=1,20)
     &    /5,2, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     ser
      data (lpmfrs(3,i),i=1,20)
     &    /5,2, 31,7, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     cys
      data (lpmfrs(4,i),i=1,20)
     &    /5,2, 37,13, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     thr
      data (lpmfrs(5,i),i=1,20)
     &    /5,2, 32,7, 8,2, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     ile
      data (lpmfrs(6,i),i=1,20)
     &    /5,2, 8,2, 7,2, 10,2, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     val
      data (lpmfrs(7,i),i=1,20)
     &    /5,2, 7,2, 8,2, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     met
      data (lpmfrs(8,i),i=1,20)
     &    /5,2, 6,2, 36,13, 12,2, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     asp
      data (lpmfrs(9,i),i=1,20)
     &    /5,2, 6,4, 29,8, 30,8, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     asn N 
      data (lpmfrs(10,i),i=1,20)
     &    /5,2, 6,4, 29,6, 21,10, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     leu
      data (lpmfrs(11,i),i=1,20)
     &    /5,2, 6,2, 10,2, 11,2, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     lys
      data (lpmfrs(12,i),i=1,20)
     &    /5,2, 6,2, 9,2, 12,2, 27,11, 0,0, 0,0, 0,0, 0,0, 0,0/
c     glu
      data (lpmfrs(13,i),i=1,20)
     &    /5,2, 6,2, 9,4, 34,8, 35,8, 0,0, 0,0, 0,0, 0,0, 0,0/
c     gln N -> Nam, niet Npl3
      data (lpmfrs(14,i),i=1,20)
     &    /5,2, 6,2, 9,4, 34,6, 24,10, 0,0, 0,0, 0,0, 0,0, 0,0/
c     pro 
      data (lpmfrs(15,i),i=1,20)
     &    /5,2, 6,2, 9,2, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0/
c     arg
      data (lpmfrs(16,i),i=1,20)
     &    /5,2, 6,2, 9,2, 22,9, 17,5, 25,9, 26,9, 0,0, 0,0, 0,0/
c     his (gekozen voor 2 x N.pl3, moet zijn N2 en N.pl3
      data (lpmfrs(17,i),i=1,20)
     &    /5,2, 6,4, 20,9, 11,4, 13,4, 24,9, 0,0, 0,0, 0,0, 0,0/
c     phe
      data (lpmfrs(18,i),i=1,20)
     &    /5,2, 6,3, 10,3, 11,3, 13,3, 14,3, 17,3, 0,0, 0,0, 0,0/
c     tyr
      data (lpmfrs(19,i),i=1,20)
     &    /5,2, 6,3, 10,3, 11,3, 13,3, 14,3, 17,3, 33,7, 0,0, 0,0/
c     trp 
      data (lpmfrs(20,i),i=1,20)
     &    /5,2, 6,4, 11,3, 14,3, 15,3, 10,4, 23,9, 18,3, 19,3, 16,3/

c add 3.0 to values
c PAIR CF CF
      data (pmf( 1, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.4182,  0.0000,  0.0000,  0.0000,  0.1469, -0.1548,
     & -0.2771, -0.4697, -0.8331, -1.4843, -2.0998, -2.4846,
     & -2.9138, -3.1265, -3.2178, -3.2463, -3.2239, -3.1897,
     & -3.2106, -3.1748, -3.1702, -3.1140, -3.0983, -3.0767,
     & -3.0541, -3.0123, -3.0034, -2.9908, -2.9956, -2.9857,
     & -3.0123, -2.9863, -3.0115, -2.9874, -2.9887, -2.9859,
     & -2.9795, -3.0014, -3.0235, -3.0060, -3.0097, -2.9917,
     & -3.0033, -3.0133, -3.0125, -3.0168, -2.9957, -3.0117,
     & -3.0074, -2.9873, -2.9799, -3.0042, -3.0021, -2.9979/
c PAIR CF CP
      data (pmf( 1, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000, -0.9953,  0.0000, -0.4644,
     &  0.0000, -0.0991,  0.0000,  0.0000,  0.0000, -0.2434,
     &  0.0000,  0.6674, -0.3628, -1.2658, -1.6809, -2.2643,
     & -2.6577, -2.8107, -2.9706, -3.0261, -3.0529, -3.0725,
     & -3.0449, -3.0365, -3.0741, -3.0850, -3.0682, -3.0597,
     & -3.0352, -3.0221, -3.0192, -3.0114, -2.9856, -2.9577,
     & -2.9696, -2.9907, -2.9730, -2.9980, -3.0149, -3.0068,
     & -3.0268, -3.0099, -3.0174, -3.0189, -3.0081, -3.0177,
     & -3.0173, -3.0334, -3.0193, -3.0297, -3.0239, -3.0317,
     & -3.0257, -3.0203, -3.0038, -3.0066, -3.0060, -2.9962/
c PAIR CF cF
      data (pmf( 1, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.1056,  0.0000,
     &  0.0000,  0.0428, -0.9067, -1.5722, -2.2672, -2.7808,
     & -3.0469, -3.1765, -3.3100, -3.3019, -3.2622, -3.2745,
     & -3.3228, -3.2565, -3.2366, -3.2283, -3.1358, -3.1166,
     & -3.1237, -3.0785, -3.0405, -3.0040, -2.9999, -2.9565,
     & -2.9493, -2.9616, -2.9570, -2.9617, -2.9665, -2.9693,
     & -2.9475, -2.9557, -2.9764, -2.9835, -2.9716, -2.9787,
     & -3.0103, -3.0116, -3.0315, -3.0236, -2.9974, -3.0110,
     & -3.0039, -3.0037, -2.9829, -2.9888, -2.9778, -2.9932/
c PAIR CF cP
      data (pmf( 1, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.4159,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.1817,  0.0000, -0.3413, -1.2607, -2.1536, -2.7910,
     & -3.1483, -3.2931, -3.2954, -3.2462, -3.2762, -3.2722,
     & -3.2340, -3.1960, -3.1565, -3.1467, -3.0886, -3.0881,
     & -3.0438, -2.9829, -2.9707, -2.9817, -2.9783, -2.9734,
     & -2.9883, -2.9833, -2.9922, -2.9856, -2.9790, -2.9734,
     & -2.9920, -2.9924, -3.0296, -3.0135, -3.0193, -3.0185,
     & -3.0148, -3.0254, -3.0152, -3.0098, -2.9875, -3.0068,
     & -3.0016, -2.9934, -2.9737, -2.9840, -2.9902, -2.9805/
c PAIR CF C3
      data (pmf( 1, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -0.5923,
     & -0.4628,  0.0000, -0.2299, -1.4829, -2.0116, -2.5612,
     & -2.8455, -3.1346, -3.1110, -3.1054, -3.2195, -3.2172,
     & -3.2037, -3.2426, -3.1386, -3.1414, -3.0929, -3.0881,
     & -3.0754, -3.0568, -3.0660, -3.0528, -3.0232, -3.0151,
     & -3.0048, -3.0286, -2.9815, -3.0153, -3.0218, -3.0346,
     & -3.0497, -3.0396, -3.0151, -3.0329, -3.0262, -2.9924,
     & -2.9792, -2.9768, -2.9987, -3.0107, -3.0163, -2.9904,
     & -2.9691, -2.9765, -2.9890, -2.9946, -2.9842, -2.9748/
c PAIR CF CW
      data (pmf( 1, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.9606,  0.0000, -0.6820,  0.0000,  0.0000,
     & -0.3205, -0.6031, -0.4853, -0.6141, -1.5337, -2.3555,
     & -2.5491, -2.9369, -3.0076, -3.0860, -3.0879, -3.1033,
     & -3.1438, -3.1690, -3.1303, -3.0670, -3.1007, -3.1035,
     & -3.0564, -3.0384, -3.0067, -3.0254, -3.0164, -3.0300,
     & -3.0382, -3.0180, -3.0285, -3.0439, -3.0358, -3.0152,
     & -3.0104, -2.9968, -3.0181, -3.0340, -3.0108, -3.0244,
     & -3.0210, -3.0388, -3.0034, -3.0018, -3.0097, -3.0119,
     & -2.9850, -2.9914, -2.9823, -2.9959, -2.9943, -2.9715/
c PAIR CF CO
      data (pmf( 1, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.2470, -1.9527, -2.2711,
     & -2.7301, -3.0463, -3.1917, -3.1261, -3.1362, -3.0725,
     & -3.1178, -3.1590, -3.1156, -3.0394, -3.0048, -3.0430,
     & -2.9574, -3.0436, -3.0611, -3.0385, -3.0131, -3.0080,
     & -3.0180, -3.0016, -3.0078, -3.0474, -3.0552, -3.0380,
     & -2.9934, -2.9762, -2.9736, -2.9998, -3.0206, -3.0092,
     & -2.9966, -3.0173, -2.9932, -2.9952, -3.0524, -3.0336,
     & -3.0125, -3.0053, -2.9960, -2.9902, -3.0001, -2.9782/
c PAIR CF CN
      data (pmf( 1, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -0.5743, -1.2889, -2.1108, -2.7217,
     & -3.1036, -3.2827, -3.3797, -3.3774, -3.3752, -3.3376,
     & -3.3196, -3.2587, -3.2271, -3.1774, -3.1113, -3.0286,
     & -2.9934, -3.0025, -2.9268, -2.9147, -2.9502, -2.9247,
     & -2.9423, -2.9422, -2.9992, -2.9649, -2.9557, -2.9982,
     & -2.9944, -3.0230, -3.0574, -3.0221, -3.0144, -3.0007,
     & -3.0166, -3.0194, -3.0101, -3.0304, -3.0138, -2.9801,
     & -3.0005, -3.0053, -2.9853, -2.9743, -2.9422, -2.9532/
c PAIR CF NC
      data (pmf( 1, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.3634, -1.3158, -1.5194, -1.9440, -2.5676, -2.9990,
     & -3.1968, -3.4225, -3.4216, -3.3826, -3.2922, -3.3717,
     & -3.3833, -3.3536, -3.2475, -3.1711, -3.0568, -3.0325,
     & -2.9666, -2.8981, -2.9680, -2.9274, -2.8979, -2.8745,
     & -2.9096, -2.9443, -2.9840, -3.0238, -3.0546, -3.0355,
     & -3.0176, -3.0124, -3.0588, -3.0404, -3.0678, -3.0008,
     & -2.9935, -3.0019, -2.9838, -2.9958, -2.9903, -2.9689,
     & -2.9614, -2.9534, -2.9640, -2.9729, -2.9654, -2.9700/
c PAIR CF NP
      data (pmf( 1,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.9704,  0.0000,  0.0000, -1.1700, -2.0310, -2.4918,
     & -2.9777, -3.2232, -3.3025, -3.2368, -3.2785, -3.2859,
     & -3.2448, -3.2067, -3.1435, -3.1049, -3.0594, -3.0752,
     & -3.0347, -2.9727, -3.0213, -3.0030, -2.9864, -2.9827,
     & -2.9798, -3.0166, -3.0307, -2.9954, -2.9929, -2.9808,
     & -2.9938, -3.0211, -3.0122, -3.0127, -3.0050, -3.0002,
     & -3.0382, -3.0109, -3.0248, -3.0220, -2.9752, -2.9661,
     & -2.9816, -2.9944, -3.0206, -2.9919, -2.9929, -2.9685/
c PAIR CF NA
      data (pmf( 1,11,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.4673,
     & -2.6107, -2.5199, -2.4334, -2.9992, -2.7698, -3.0574,
     & -3.0268, -3.3206, -3.1535, -3.1117, -3.0713, -3.1547,
     & -2.9248, -3.0623, -3.0557, -2.9741, -2.8947, -3.0674,
     & -2.7765, -3.0231, -2.8271, -3.0209, -3.0921, -3.0168,
     & -3.0187, -2.9078, -2.9911, -2.9321, -3.0720, -3.0122,
     & -3.0880, -3.0333, -3.1061, -2.9845, -3.0420, -3.0262,
     & -3.0944, -2.9066, -3.0346, -3.0196, -2.9073, -3.0334/
c PAIR CF ND
      data (pmf( 1,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.8390,  0.0000, -1.2469, -2.0469, -3.0755,
     & -2.7681, -3.1180, -3.0191, -3.2182, -3.0480, -3.0599,
     & -3.0806, -3.1576, -3.1626, -3.1187, -3.1135, -3.1783,
     & -3.1241, -3.2435, -3.0081, -2.9928, -2.9771, -3.0674,
     & -3.0839, -3.0607, -2.9996, -2.9985, -3.0626, -3.0920,
     & -3.0832, -3.0279, -3.0284, -3.0033, -3.0553, -3.0487,
     & -3.0078, -3.0232, -3.0588, -3.0195, -3.0084, -2.9430,
     & -2.9391, -2.9660, -2.9765, -2.9408, -2.9438, -2.9211/
c PAIR CF NR
      data (pmf( 1,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.0530, -0.8900, -1.7098, -2.5138, -3.0828,
     & -3.3595, -3.4094, -3.3672, -3.3235, -3.3174, -3.3111,
     & -3.3685, -3.2217, -3.1506, -3.1381, -3.1017, -3.0838,
     & -3.0633, -2.9899, -3.0117, -2.9964, -2.9852, -2.9962,
     & -2.9872, -2.9934, -2.9510, -3.0020, -2.9948, -3.0303,
     & -3.0164, -2.9846, -3.0065, -3.0311, -3.0366, -3.0429,
     & -3.0218, -3.0141, -3.0140, -2.9938, -2.9930, -2.9720,
     & -2.9585, -2.9499, -2.9623, -2.9686, -2.9572, -2.9448/
c PAIR CF N0
      data (pmf( 1,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CF NS
      data (pmf( 1,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.0762, -0.9707, -1.8211, -2.4496, -2.7480,
     & -2.8931, -3.0604, -2.9882, -3.0618, -2.9095, -3.0947,
     & -3.0120, -3.1472, -3.0464, -3.2730, -3.1571, -3.0628,
     & -3.1555, -2.9513, -3.1185, -2.9992, -3.0320, -2.9785,
     & -2.9634, -3.0563, -3.0466, -2.9736, -3.0177, -3.0322,
     & -3.0045, -2.9664, -3.0245, -2.9334, -2.9886, -2.9666,
     & -3.0125, -3.0522, -3.0361, -2.9915, -2.9339, -3.0500,
     & -3.0282, -2.9551, -3.0204, -2.9903, -2.9740, -3.0258/
c PAIR CF OC
      data (pmf( 1,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.0048, -0.5501, -1.2205, -2.1112, -2.6975, -2.9850,
     & -3.0945, -3.0991, -3.0705, -3.0057, -3.0476, -3.1420,
     & -3.1426, -3.0465, -3.0562, -3.0549, -3.0272, -3.0253,
     & -3.0152, -3.0254, -3.0016, -3.0099, -3.0018, -3.0238,
     & -3.0378, -3.0235, -3.0486, -3.0206, -2.9990, -3.0222,
     & -3.0029, -2.9901, -3.0179, -3.0074, -3.0128, -3.0043,
     & -2.9937, -3.0028, -2.9920, -3.0001, -3.0045, -2.9912,
     & -3.0041, -2.9985, -3.0137, -3.0052, -3.0028, -3.0181/
c PAIR CF OA
      data (pmf( 1,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.7402,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.9424, -0.6686, -1.7469, -2.0778, -2.7653, -3.1684,
     & -3.3470, -3.3024, -3.2210, -3.1495, -3.0624, -3.1063,
     & -3.1471, -3.0613, -3.0764, -3.0654, -3.0441, -2.9958,
     & -3.0308, -3.0087, -3.0615, -3.0205, -3.0137, -3.0070,
     & -3.0357, -3.0576, -3.0573, -3.0772, -3.0454, -3.0437,
     & -3.0220, -3.0136, -3.0062, -3.0052, -2.9672, -3.0071,
     & -2.9880, -2.9899, -2.9680, -2.9640, -3.0277, -2.9655,
     & -2.9733, -2.9894, -2.9853, -3.0030, -2.9766, -3.0164/
c PAIR CF OE
      data (pmf( 1,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.1088, -0.9699, -1.4893, -2.2989, -2.7311, -3.0663,
     & -3.0538, -3.2342, -3.1753, -3.1069, -3.1592, -3.1134,
     & -3.0859, -3.0206, -3.0163, -2.9768, -3.0675, -3.0752,
     & -3.0359, -3.0049, -2.9287, -2.9384, -2.9752, -2.9559,
     & -2.9523, -3.0032, -3.0054, -2.9873, -2.9965, -3.0282,
     & -3.0326, -3.0224, -3.0100, -3.0199, -3.0185, -3.0328,
     & -3.0279, -3.0133, -3.0238, -3.0175, -3.0358, -3.0025,
     & -2.9958, -2.9730, -3.0136, -3.0249, -3.0214, -2.9901/
c PAIR CF OR
      data (pmf( 1,19,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.9285, -2.8185, -3.1119,
     & -3.2604, -3.4528, -3.5815, -3.1618, -3.2523, -2.9167,
     & -2.9249, -2.6126, -2.6202, -2.6695, -2.8484, -2.9792,
     & -3.0442, -3.0438, -2.9919, -2.9596, -3.0648, -3.0383,
     & -3.0090, -3.0399, -3.1034, -3.0591, -3.1003, -3.0255,
     & -3.0006, -2.9717, -3.0966, -3.0193, -3.0697, -3.0399,
     & -2.9794, -3.0320, -3.0732, -2.9065, -2.9265, -2.9624,
     & -3.0244, -2.9164, -2.9425, -3.0123, -3.0515, -3.0371/
c PAIR CF OS
      data (pmf( 1,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.3817,
     & -0.4567, -0.9170, -1.4981, -2.1881, -2.6744, -2.9534,
     & -3.0142, -2.9639, -2.8829, -2.8860, -2.9033, -2.9175,
     & -2.9223, -2.8754, -2.8855, -2.9478, -2.9335, -2.9347,
     & -2.9169, -2.9113, -2.9123, -2.9284, -2.9325, -2.9454,
     & -2.9821, -3.0005, -3.0218, -3.0298, -3.0475, -3.0104,
     & -3.0215, -3.0373, -3.0383, -3.0433, -3.0380, -3.0323,
     & -3.0366, -3.0294, -3.0267, -3.0332, -3.0265, -3.0189,
     & -3.0320, -3.0297, -3.0428, -3.0383, -3.0340, -3.0274/
c PAIR CF OD
      data (pmf( 1,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -0.1557,
     & -0.9256, -0.5772, -1.5154, -2.0731, -2.6417, -2.9305,
     & -3.0529, -2.9848, -2.9536, -3.0086, -3.0521, -3.0936,
     & -3.1478, -3.0901, -3.0306, -3.0386, -2.9727, -2.9654,
     & -2.9500, -2.9677, -2.9698, -2.9989, -3.0076, -2.9797,
     & -2.9900, -3.0018, -2.9954, -3.0153, -2.9730, -2.9786,
     & -3.0012, -2.9972, -3.0072, -3.0135, -2.9931, -3.0183,
     & -3.0208, -3.0435, -3.0134, -3.0093, -3.0314, -3.0258,
     & -3.0165, -3.0215, -3.0316, -3.0312, -3.0224, -3.0069/
c PAIR CF P
      data (pmf( 1,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -0.6380, -0.7760, -1.3908, -1.8828,
     & -2.3209, -2.7908, -3.0346, -3.1076, -3.1699, -3.0624,
     & -2.9623, -2.9700, -3.0119, -2.9497, -2.8787, -2.9312,
     & -2.9254, -2.7875, -2.8371, -2.8811, -2.9436, -2.8612,
     & -2.9661, -3.0192, -2.9963, -3.0300, -3.0563, -3.1192,
     & -3.0770, -3.0736, -3.0840, -3.0919, -3.0539, -3.0404,
     & -3.0374, -3.0094, -3.0052, -3.0306, -3.0236, -3.0658,
     & -3.0084, -3.0328, -3.0677, -3.0340, -3.0392, -3.0167/
c PAIR CF SA
      data (pmf( 1,23,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.1108, -2.0145, -2.3266,
     & -2.2349, -3.2048, -3.1690, -3.1741, -3.2069, -3.4165,
     & -3.4414, -3.2605, -3.3920, -3.2793, -3.0886, -3.0206,
     & -2.8051, -3.0620, -2.7959, -2.7285, -2.8276, -2.8110,
     & -2.8990, -3.1486, -2.9308, -2.9093, -2.9375, -3.1028,
     & -3.0701, -3.0297, -2.9289, -3.0615, -3.1170, -3.0843,
     & -3.1595, -2.9932, -2.9103, -2.9515, -2.9971, -3.1245,
     & -2.9344, -2.9007, -2.9743, -3.0008, -2.9387, -3.0088/
c PAIR CF SD
      data (pmf( 1,24,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.3425,
     & -2.2571, -2.5848, -2.7465,  0.0000, -2.7707, -2.5326,
     & -1.8186, -3.2212, -3.2511, -2.7812, -3.2796, -3.2633,
     & -2.9927, -2.7925, -2.7483, -2.9038, -3.2546, -3.3721,
     & -3.0635, -3.3815, -3.1778, -2.6951, -3.0719, -3.0189,
     & -3.0098, -3.2465, -3.1913, -2.7264, -3.0300, -2.9507,
     & -2.9430, -2.7299, -3.0354, -2.9520, -2.7907, -2.8065,
     & -3.0597, -2.9968, -2.9506, -3.0135, -3.1044, -3.0751/
c PAIR CF HL
      data (pmf( 1,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -0.6478,
     & -1.4834, -2.1952, -2.5268, -2.9554, -2.9843, -3.0399,
     & -3.0988, -3.1175, -3.1310, -3.1405, -3.1196, -3.1176,
     & -3.0506, -3.0737, -3.0440, -3.0819, -3.0647, -3.0246,
     & -3.0214, -3.0260, -3.0018, -3.0141, -3.0092, -3.0415,
     & -3.0372, -2.9945, -3.0146, -3.0185, -3.0301, -3.0358,
     & -3.0203, -3.0165, -3.0194, -2.9912, -2.9954, -2.9857,
     & -2.9967, -3.0315, -3.0107, -3.0156, -3.0111, -3.0117,
     & -3.0138, -3.0100, -3.0073, -3.0092, -3.0021, -3.0034/
c PAIR CF F
      data (pmf( 1,26,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.7680, -2.6228, -2.9422, -3.1665,
     & -3.3766, -3.4279, -3.3166, -3.1670, -3.1887, -3.2441,
     & -3.4442, -3.4555, -3.2106, -3.1805, -3.1313, -2.9985,
     & -3.0391, -3.2376, -2.9905, -2.9389, -3.0510, -2.8634,
     & -3.0132, -2.8279, -2.8682, -2.9091, -2.9300, -2.9702,
     & -2.9667, -3.0787, -2.9136, -2.9513, -2.9147, -2.9818,
     & -2.9760, -2.9932, -2.9797, -3.0011, -3.0400, -3.0024,
     & -3.0019, -2.9216, -2.9630, -3.0504, -3.0048, -2.9874/
c PAIR CP CF
      data (pmf( 2, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.3320, -0.8333,  0.0000,  0.0000,  0.0000, -0.7384,
     & -1.2445, -0.7006, -0.8199, -1.2794, -1.5455, -2.1967,
     & -2.4329, -2.6419, -2.7068, -2.7956, -2.8591, -2.8591,
     & -2.8943, -2.9015, -2.9228, -2.9218, -2.9332, -2.9609,
     & -2.9502, -2.9374, -2.9455, -2.9560, -2.9450, -2.9566,
     & -2.9778, -3.0046, -3.0158, -3.0316, -3.0345, -3.0476,
     & -3.0566, -3.0442, -3.0624, -3.0492, -3.0471, -3.0583,
     & -3.0586, -3.0620, -3.0513, -3.0500, -3.0323, -3.0499,
     & -3.0346, -3.0192, -3.0295, -3.0179, -2.9982, -2.9935/
c PAIR CP CP
      data (pmf( 2, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.4684,  0.0000,  0.0000,  0.4449,  0.0000, -0.5567,
     & -1.2003, -0.5806, -0.8716, -1.3455, -1.7822, -2.2224,
     & -2.3574, -2.6167, -2.7732, -2.8810, -2.9915, -3.0038,
     & -2.9955, -3.0141, -2.9721, -2.9395, -2.9273, -2.9539,
     & -2.9518, -2.9685, -2.9388, -2.9387, -2.9478, -2.9711,
     & -2.9891, -2.9870, -3.0023, -2.9913, -2.9992, -3.0264,
     & -3.0434, -3.0709, -3.0625, -3.0540, -3.0632, -3.0517,
     & -3.0487, -3.0615, -3.0580, -3.0488, -3.0348, -3.0221,
     & -3.0292, -3.0239, -3.0179, -3.0049, -2.9880, -2.9717/
c PAIR CP cF
      data (pmf( 2, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.6944, -0.4687,  0.0000, -1.5184,  0.3369, -0.3738,
     & -1.8071, -1.0845, -1.0331, -1.7825, -1.8892, -2.3263,
     & -2.5400, -2.7117, -2.8195, -2.8461, -2.8655, -2.8560,
     & -2.8761, -2.8841, -2.9190, -2.9234, -2.9349, -2.9509,
     & -2.9726, -2.9350, -2.9542, -2.9588, -2.9663, -2.9674,
     & -2.9832, -3.0099, -3.0147, -3.0258, -3.0241, -3.0348,
     & -3.0430, -3.0432, -3.0548, -3.0462, -3.0556, -3.0636,
     & -3.0633, -3.0573, -3.0652, -3.0502, -3.0514, -3.0351,
     & -3.0335, -3.0188, -3.0109, -2.9931, -2.9938, -2.9876/
c PAIR CP cP
      data (pmf( 2, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.0109, -0.1942, -0.0467,  0.0000,  0.0000, -0.0993,
     & -0.9095, -0.4790, -0.2512, -1.4081, -1.7381, -2.4093,
     & -2.6199, -2.8616, -2.9394, -2.9932, -3.0397, -3.0521,
     & -3.0300, -3.0156, -2.9940, -2.9621, -2.9935, -3.0123,
     & -3.0250, -3.0208, -3.0155, -3.0054, -3.0073, -3.0142,
     & -3.0068, -3.0273, -3.0199, -3.0280, -3.0581, -3.0512,
     & -3.0398, -3.0521, -3.0580, -3.0715, -3.0492, -3.0205,
     & -3.0204, -3.0288, -3.0312, -3.0298, -3.0197, -3.0042,
     & -2.9836, -2.9810, -2.9925, -2.9905, -2.9907, -2.9632/
c PAIR CP C3
      data (pmf( 2, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.4654, -1.2969,  0.0000, -0.6093,  0.0000, -0.3840,
     & -1.0724, -1.0837, -0.4306, -1.5877, -1.7826, -2.3663,
     & -2.4874, -2.5933, -2.8305, -2.8727, -2.9219, -2.8681,
     & -2.8752, -2.9480, -2.9725, -2.9516, -2.8994, -2.9245,
     & -2.9057, -2.9245, -2.8834, -2.9365, -2.9019, -2.9629,
     & -2.9704, -2.9987, -3.0211, -2.9951, -3.0056, -3.0271,
     & -3.0628, -3.0706, -3.0652, -3.0863, -3.0735, -3.0674,
     & -3.0533, -3.0718, -3.0611, -3.0516, -3.0287, -3.0319,
     & -3.0276, -3.0256, -3.0257, -3.0133, -3.0067, -3.0071/
c PAIR CP CW
      data (pmf( 2, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.7063,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.8841, -0.3488, -1.2368, -1.3094, -1.6338, -2.1561,
     & -2.3329, -2.5773, -2.6996, -2.8556, -3.0259, -3.0717,
     & -3.0312, -2.9868, -2.9599, -2.9326, -2.9222, -2.9021,
     & -2.9570, -2.9515, -2.9281, -2.9109, -2.9386, -2.9450,
     & -2.9611, -2.9689, -2.9558, -3.0019, -2.9886, -3.0308,
     & -3.0247, -3.0327, -3.0457, -3.0452, -3.0694, -3.0489,
     & -3.0258, -3.0706, -3.0691, -3.0557, -3.0583, -3.0555,
     & -3.0427, -3.0277, -3.0188, -3.0289, -2.9970, -3.0075/
c PAIR CP CO
      data (pmf( 2, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.7897,  0.0000, -0.7478, -1.2998, -1.8892, -1.9486,
     & -2.2776, -2.6669, -2.8939, -2.9150, -2.9598, -2.9316,
     & -2.9233, -2.8921, -2.8143, -2.7803, -2.8836, -2.8294,
     & -2.8972, -2.9198, -2.9770, -2.9665, -2.9481, -2.9769,
     & -3.0050, -3.0101, -3.0533, -3.0422, -3.0306, -3.0104,
     & -3.0187, -3.0493, -3.0656, -3.0731, -3.0865, -3.0706,
     & -3.0320, -3.0509, -3.0408, -3.0473, -3.0381, -3.0406,
     & -3.0213, -3.0284, -3.0153, -3.0084, -3.0057, -3.0307/
c PAIR CP CN
      data (pmf( 2, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.9778, -0.8093,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.1788, -0.7093, -0.7690, -2.1439, -2.0928, -2.3089,
     & -2.5383, -2.8292, -2.9763, -3.0392, -3.1819, -3.2324,
     & -3.1924, -3.1140, -3.0738, -3.0151, -3.0405, -3.0360,
     & -3.0368, -3.0678, -3.1150, -3.0759, -3.0275, -2.9923,
     & -2.9760, -3.0244, -3.0407, -3.0194, -3.0328, -3.0348,
     & -3.0509, -3.0526, -3.0368, -3.0306, -3.0589, -3.0397,
     & -3.0406, -3.0504, -3.0072, -3.0153, -2.9933, -2.9848,
     & -2.9905, -2.9823, -2.9880, -2.9445, -2.9378, -2.9373/
c PAIR CP NC
      data (pmf( 2, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.7935, -0.0407, -1.0909, -1.7954, -2.2792, -2.8364,
     & -3.1000, -3.2574, -3.2568, -3.1426, -3.0255, -3.0692,
     & -3.0856, -3.0705, -3.0796, -3.1309, -3.0874, -3.0369,
     & -3.0239, -3.0241, -2.9483, -2.9605, -3.0100, -2.9980,
     & -3.0358, -3.0148, -3.0068, -3.0261, -3.0292, -3.0514,
     & -3.0456, -3.0713, -3.0900, -3.1006, -3.0528, -3.0521,
     & -3.0544, -3.0069, -3.0319, -3.0221, -3.0186, -2.9800,
     & -2.9734, -2.9937, -2.9630, -2.9341, -2.9086, -2.9105/
c PAIR CP NP
      data (pmf( 2,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.9779, -0.8291, -1.3406, -1.8284, -2.1512,
     & -2.2378, -2.5247, -2.7307, -2.9174, -2.8926, -2.9987,
     & -3.0631, -3.0564, -3.1175, -3.0528, -3.0934, -3.0970,
     & -3.0818, -3.0328, -3.0000, -2.9774, -2.9829, -2.9694,
     & -3.0012, -3.0528, -3.0588, -3.0762, -3.0404, -3.0516,
     & -3.0133, -2.9910, -3.0087, -2.9844, -3.0407, -3.0319,
     & -3.0319, -3.0230, -3.0339, -3.0390, -3.0182, -3.0529,
     & -3.0225, -3.0005, -2.9614, -2.9684, -2.9756, -2.9742/
c PAIR CP NA
      data (pmf( 2,11,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.3486,  0.0000,  0.0000, -2.0067,  0.0000, -1.8050,
     & -2.1182, -2.5679, -1.9409,  0.0000, -2.0165, -2.8041,
     & -2.6312, -2.8166, -2.9002, -2.6800, -2.9089, -2.9394,
     & -2.7700, -2.8828, -2.9340, -3.0103, -2.9213, -3.1272,
     & -3.0067, -2.9995, -2.9734, -3.0110, -2.9704, -2.9555,
     & -3.0940, -2.9916, -2.9904, -3.0984, -3.1726, -3.0804,
     & -2.9215, -2.9720, -2.9688, -3.1325, -3.0409, -2.9788,
     & -3.1335, -3.0319, -3.0262, -3.0277, -3.0825, -3.0657/
c PAIR CP ND
      data (pmf( 2,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.1800, -1.5496, -2.2375,
     & -3.0090, -3.0501, -3.0362, -2.8312, -2.9052, -2.9080,
     & -2.9338, -2.9688, -2.9177, -2.9961, -2.9459, -2.9760,
     & -2.8961, -2.8991, -2.9311, -2.8710, -2.8703, -2.8950,
     & -2.9285, -2.9966, -2.9655, -3.0004, -2.9636, -3.0183,
     & -3.0314, -3.0519, -3.0668, -2.9706, -3.0268, -3.0379,
     & -3.0832, -3.0753, -3.1089, -3.0609, -3.0209, -3.0621,
     & -3.0507, -3.0408, -3.0249, -3.0066, -3.0434, -3.0362/
c PAIR CP NR
      data (pmf( 2,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0259, -0.6715, -1.1558, -1.7658, -2.3985, -2.9555,
     & -3.1918, -3.2028, -3.0992, -2.9733, -3.0407, -3.0866,
     & -3.0579, -3.1290, -3.0946, -3.0173, -3.0150, -3.0460,
     & -3.0343, -3.0276, -3.0227, -3.0241, -3.0462, -3.0107,
     & -3.0423, -3.0325, -3.0501, -3.0586, -3.0762, -3.0492,
     & -3.0431, -3.0391, -3.0372, -3.0649, -3.0462, -3.0394,
     & -3.0003, -2.9994, -3.0092, -2.9931, -3.0025, -2.9902,
     & -2.9821, -2.9747, -2.9823, -2.9810, -2.9493, -2.9343/
c PAIR CP N0
      data (pmf( 2,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CP NS
      data (pmf( 2,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.5103,  0.0000,  0.0000, -1.8934,  0.0000,  0.0000,
     & -2.0468, -1.1949, -1.6300, -1.7294, -2.2250, -2.6140,
     & -2.5656, -2.8440, -3.0147, -2.9291, -3.0282, -2.9764,
     & -2.8624, -2.9152, -2.8795, -2.7850, -2.9091, -2.9042,
     & -2.8070, -2.7742, -2.8209, -2.8714, -2.9597, -3.0321,
     & -3.0559, -2.8870, -2.9186, -3.0059, -3.1298, -2.9615,
     & -3.0555, -3.0106, -3.0554, -3.0710, -3.1420, -3.0674,
     & -3.0205, -3.0938, -3.0300, -3.1085, -3.0809, -3.0505,
     & -2.9802, -3.0602, -2.9866, -3.0274, -3.0059, -2.9646/
c PAIR CP OC
      data (pmf( 2,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -0.7292,  0.0000,  0.0000,
     & -0.6946,  0.3578, -1.5528, -2.2127, -2.6353, -2.9201,
     & -2.9772, -2.8858, -2.8453, -2.8072, -2.7887, -2.7983,
     & -2.8197, -2.8686, -2.8480, -2.9269, -2.9211, -2.9241,
     & -2.9628, -2.9451, -2.9332, -2.9790, -2.9808, -2.9883,
     & -3.0243, -3.0159, -3.0320, -3.0435, -3.0296, -3.0298,
     & -3.0546, -3.0576, -3.0277, -3.0162, -3.0169, -3.0312,
     & -3.0337, -3.0403, -3.0221, -3.0276, -3.0406, -3.0192,
     & -3.0284, -3.0316, -3.0479, -3.0231, -3.0203, -3.0158/
c PAIR CP OA
      data (pmf( 2,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0032,
     & -0.5331, -1.0078, -1.5074, -2.5903, -2.9020, -3.2445,
     & -3.3136, -3.1891, -3.0925, -2.9007, -3.0202, -2.9444,
     & -2.8513, -2.8520, -2.9101, -2.8929, -2.8799, -2.9224,
     & -2.9465, -2.9507, -2.9257, -2.8714, -2.9273, -2.9614,
     & -2.9817, -2.9975, -3.0160, -3.0497, -3.0407, -3.0564,
     & -3.0683, -3.0610, -3.0281, -3.0486, -3.0469, -3.0266,
     & -3.0153, -3.0190, -3.0148, -3.0203, -3.0126, -3.0308,
     & -3.0381, -3.0278, -3.0087, -2.9892, -3.0198, -3.0124/
c PAIR CP OE
      data (pmf( 2,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.2585,  0.0000,  0.0000, -0.8113,  0.0000, -0.5860,
     & -1.3849, -0.9446, -1.2248, -2.1619, -2.7032, -2.8986,
     & -2.8915, -2.9373, -2.8586, -2.9252, -2.9098, -2.8567,
     & -2.9634, -2.9564, -2.9443, -2.9260, -2.8965, -2.9166,
     & -2.9864, -2.9928, -2.9773, -2.9715, -2.9872, -2.9680,
     & -2.9932, -2.9775, -3.0378, -3.0274, -2.9849, -3.0848,
     & -3.1002, -3.0632, -3.0731, -3.0506, -3.0343, -3.0281,
     & -3.0551, -3.0546, -3.0262, -3.0185, -3.0335, -3.0251,
     & -3.0361, -2.9923, -2.9860, -2.9778, -2.9997, -2.9824/
c PAIR CP OR
      data (pmf( 2,19,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.7258, -1.5815, -1.8572, -2.7904, -3.2345, -3.4783,
     & -3.3803, -3.0736, -3.0824, -3.3020, -3.1813, -2.9700,
     & -2.8589, -2.9471, -3.0935, -3.2363, -3.1017, -3.1588,
     & -3.0237, -3.0284, -2.9641, -3.0648, -3.0726, -2.9828,
     & -2.8984, -2.9852, -3.0056, -3.0678, -3.0380, -3.1228,
     & -3.0727, -3.1573, -3.1479, -3.0043, -2.9805, -2.9302,
     & -2.9912, -3.0793, -3.0084, -2.9690, -2.9394, -2.8921,
     & -2.9636, -2.9827, -2.9217, -2.9529, -2.9813, -2.9686/
c PAIR CP OS
      data (pmf( 2,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.0432, -0.2836,  0.0000,  0.0000,  0.0000,  0.2203,
     & -0.6743, -1.1428, -1.7078, -2.3720, -2.9358, -3.2252,
     & -3.3507, -3.2601, -3.0879, -3.0094, -2.9808, -2.9985,
     & -3.0284, -3.0426, -3.0891, -3.1170, -3.1176, -3.0903,
     & -3.0502, -3.0262, -2.9826, -2.9619, -2.9861, -3.0019,
     & -3.0372, -3.0174, -3.0202, -3.0420, -3.0628, -3.0434,
     & -3.0240, -3.0288, -3.0067, -3.0065, -3.0050, -3.0040,
     & -3.0155, -3.0167, -3.0143, -2.9976, -3.0036, -2.9859,
     & -2.9988, -2.9807, -2.9801, -2.9776, -2.9756, -2.9766/
c PAIR CP OD
      data (pmf( 2,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.1010, -0.3918,  0.0000,  0.0000,  0.4138,  0.1121,
     & -0.8748, -0.9895, -1.3264, -2.2389, -2.7215, -2.9882,
     & -3.0752, -3.0599, -2.9158, -2.8196, -2.8413, -2.9020,
     & -2.9108, -2.9063, -2.9376, -2.9875, -2.9479, -2.9290,
     & -2.9473, -2.9646, -2.9748, -2.9787, -2.9765, -2.9743,
     & -2.9739, -3.0282, -3.0437, -3.0360, -3.0224, -3.0069,
     & -3.0266, -3.0423, -3.0600, -3.0654, -3.0639, -3.0480,
     & -3.0348, -3.0469, -3.0360, -3.0303, -3.0205, -3.0159,
     & -3.0188, -3.0218, -3.0113, -2.9920, -2.9749, -3.0000/
c PAIR CP P
      data (pmf( 2,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.2051, -0.4988, -0.8000, -0.8304, -1.5226, -1.9240,
     & -2.4316, -2.8903, -3.1429, -3.2949, -3.3556, -3.4367,
     & -3.4356, -3.2526, -3.1563, -2.9587, -2.9825, -3.0166,
     & -2.9940, -3.0058, -2.9600, -3.0542, -3.0986, -3.0659,
     & -3.0327, -3.0195, -3.0316, -3.0126, -2.9668, -2.9443,
     & -2.9497, -2.9985, -3.0372, -3.0511, -3.0685, -3.0731,
     & -3.0349, -2.9731, -3.0007, -3.0154, -3.0075, -3.0277,
     & -2.9963, -2.9713, -2.9857, -2.9961, -2.9812, -2.9922/
c PAIR CP SA
      data (pmf( 2,23,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.4829, -2.9148,  0.0000, -2.5848, -2.0790,
     & -1.9873, -2.4407, -2.4652, -2.8421, -2.9255, -3.0572,
     & -3.0123, -3.0535, -2.9696, -2.8620, -2.8803, -2.9047,
     & -2.8417, -2.8144, -2.7875, -2.8636, -2.9787, -3.1195,
     & -3.0105, -3.0003, -3.0988, -3.1025, -3.0432, -3.0656,
     & -3.1503, -3.1513, -3.1922, -3.2041, -3.0163, -3.0333,
     & -3.1761, -3.0688, -3.1195, -3.0514, -3.0735, -3.0787,
     & -2.8916, -2.9133, -3.0148, -2.8506, -2.8259, -2.8901/
c PAIR CP SD
      data (pmf( 2,24,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -3.2271,  0.0000,
     &  0.0000,  0.0000, -3.1944, -2.4098, -2.0993, -2.9805,
     & -2.7830, -2.9588, -2.0798, -2.5976, -2.6704, -2.8614,
     & -3.1161, -2.9952, -2.8326, -3.1103, -3.2387, -2.8358,
     & -2.9641, -3.0888, -3.1546, -2.8202, -2.8501, -3.0835,
     & -2.6313, -2.9705, -2.7684, -3.1121, -3.1143, -3.0392,
     & -2.9926, -3.1629, -3.1071, -2.8746, -3.1070, -3.3337,
     & -3.0007, -2.8147, -3.1344, -3.1600, -2.8420, -3.0031/
c PAIR CP HL
      data (pmf( 2,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -0.0149,
     & -1.2595, -1.9564, -2.3506, -2.6858, -2.8838, -2.8631,
     & -2.7240, -2.8219, -2.9425, -2.9466, -2.8898, -2.8783,
     & -2.8929, -2.8576, -2.9194, -2.9113, -2.9300, -2.9513,
     & -2.9772, -2.9714, -2.9917, -2.9822, -2.9799, -2.9641,
     & -3.0088, -3.0340, -3.0506, -3.0583, -3.0316, -3.0493,
     & -3.0521, -3.0809, -3.0577, -3.0502, -3.0568, -3.0439,
     & -3.0558, -3.0661, -3.0567, -3.0325, -3.0288, -3.0230,
     & -3.0309, -3.0117, -3.0166, -2.9968, -3.0024, -2.9857/
c PAIR CP F
      data (pmf( 2,26,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.5194, -2.2425, -1.9833, -2.4732,
     & -2.8992, -2.8189, -2.5621, -2.6689, -2.4573, -2.8532,
     & -2.9001, -2.7401, -2.8566, -2.9112, -2.9835, -2.6628,
     & -2.7804, -2.8801, -3.0262, -2.9955, -2.9222, -2.9312,
     & -3.0419, -3.0637, -3.1228, -3.0625, -2.9056, -2.9883,
     & -3.1080, -3.1930, -3.1872, -3.0383, -3.1304, -3.0882,
     & -3.0644, -2.9726, -3.0033, -3.0744, -3.0744, -2.9992,
     & -3.0246, -2.9932, -2.9407, -2.9476, -2.9771, -3.0272/
c PAIR cF CF
      data (pmf( 3, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.6324, -1.0305, -1.9172, -2.3771, -2.8995,
     & -3.2239, -3.5377, -3.4919, -3.4525, -3.4366, -3.3859,
     & -3.4105, -3.3423, -3.3116, -3.2769, -3.2455, -3.2299,
     & -3.2246, -3.1975, -3.1369, -3.1117, -3.1258, -3.0940,
     & -3.0805, -3.0385, -3.0049, -3.0053, -3.0248, -3.0012,
     & -2.9932, -2.9839, -2.9711, -3.0035, -2.9920, -2.9693,
     & -2.9925, -2.9411, -2.9389, -2.9548, -2.9088, -2.9343,
     & -2.9324, -2.9248, -2.9027, -2.8540, -2.8952, -2.8924/
c PAIR cF CP
      data (pmf( 3, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.1480, -0.7912, -0.9822, -2.0264, -2.8396,
     & -3.3618, -3.5107, -3.4947, -3.4463, -3.4465, -3.4330,
     & -3.4075, -3.3816, -3.3503, -3.2598, -3.2205, -3.2945,
     & -3.2502, -3.2310, -3.1852, -3.1366, -3.1301, -3.1201,
     & -3.0986, -3.0759, -3.0457, -3.0032, -2.9804, -2.9859,
     & -2.9728, -2.9652, -2.9214, -2.9422, -2.9279, -2.9467,
     & -2.9041, -2.9110, -2.9160, -2.9102, -2.9196, -2.9235,
     & -2.9338, -2.9254, -2.9107, -2.9295, -2.9449, -2.9582/
c PAIR cF cF
      data (pmf( 3, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -0.5785, -1.8739, -2.6577, -3.1947,
     & -3.4825, -3.5813, -3.4996, -3.4822, -3.4592, -3.4424,
     & -3.4333, -3.3720, -3.3404, -3.2847, -3.2636, -3.2686,
     & -3.2873, -3.2426, -3.1835, -3.1353, -3.1566, -3.1127,
     & -3.1007, -3.0867, -3.0407, -2.9875, -2.9809, -2.9881,
     & -2.9572, -2.9178, -2.9282, -2.9104, -2.9210, -2.9079,
     & -2.9149, -2.9375, -2.9324, -2.9330, -2.9041, -2.9129,
     & -2.9296, -2.9306, -2.9446, -2.9304, -2.9085, -2.8961/
c PAIR cF cP
      data (pmf( 3, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.0367, -1.9812, -2.6420, -3.3289,
     & -3.6271, -3.6062, -3.4922, -3.4900, -3.4096, -3.3565,
     & -3.3363, -3.2367, -3.2239, -3.2402, -3.1746, -3.1344,
     & -3.1211, -3.1270, -3.1024, -3.0693, -3.0605, -2.9837,
     & -3.0051, -2.9941, -2.9871, -2.9417, -2.9997, -2.9389,
     & -2.9510, -2.9972, -2.9978, -2.9475, -2.9446, -2.9784,
     & -2.9722, -2.9681, -2.9674, -2.9578, -2.9734, -2.9588,
     & -2.9664, -2.9685, -2.9869, -2.9625, -2.9560, -2.9332/
c PAIR cF C3
      data (pmf( 3, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.7672, -0.6547, -1.8455, -2.6823, -3.0189,
     & -3.2982, -3.3801, -3.3133, -3.4038, -3.4110, -3.4457,
     & -3.3336, -3.4246, -3.3434, -3.3010, -3.2279, -3.2568,
     & -3.2689, -3.2574, -3.1854, -3.1754, -3.2092, -3.1467,
     & -3.1613, -3.0925, -3.0368, -3.0339, -3.0170, -2.9945,
     & -2.9498, -2.9769, -2.9452, -3.0223, -3.0074, -3.0258,
     & -2.9440, -2.9258, -2.9254, -2.9479, -2.9423, -2.8940,
     & -2.8922, -2.8842, -2.8714, -2.8816, -2.9188, -2.8634/
c PAIR cF CW
      data (pmf( 3, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.7335, -1.0246, -0.5053, -2.2513, -2.8181,
     & -3.2680, -3.3778, -3.3563, -3.3320, -3.3493, -3.3877,
     & -3.3632, -3.3301, -3.3628, -3.2689, -3.3135, -3.2875,
     & -3.2629, -3.2227, -3.2295, -3.1969, -3.2087, -3.1706,
     & -3.0935, -3.0696, -2.9943, -2.9695, -3.0336, -2.9536,
     & -2.9647, -2.9384, -2.9485, -2.9487, -2.9238, -2.8910,
     & -2.9542, -2.9376, -2.9116, -2.9493, -2.8933, -2.9488,
     & -2.9387, -2.9134, -2.9434, -2.9112, -2.9479, -2.9765/
c PAIR cF CO
      data (pmf( 3, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.7143,  0.0000,  0.0000,  0.0000,  0.0000, -1.4508,
     &  0.0000,  0.0000, -0.7127, -1.5661, -1.8843, -2.7777,
     & -3.4543, -3.4886, -3.4155, -3.3483, -3.3868, -3.3503,
     & -3.4464, -3.3727, -3.4117, -3.4105, -3.3622, -3.2713,
     & -3.3115, -3.3049, -3.2319, -3.1827, -3.2238, -3.1517,
     & -3.1101, -3.0630, -3.0333, -2.9829, -3.0449, -3.0130,
     & -2.9891, -2.9692, -2.8051, -2.8951, -2.9452, -2.9228,
     & -2.9656, -2.8413, -2.7991, -2.8749, -2.8171, -2.8549,
     & -2.8999, -2.8869, -2.9215, -2.9348, -2.9689, -2.9412/
c PAIR cF CN
      data (pmf( 3, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.1977,  0.0000,  0.0000,  0.0000,  0.0000, -1.7645,
     &  0.0000, -0.8720, -1.4102, -1.8855, -2.6755, -3.2488,
     & -3.5312, -3.3636, -3.3301, -3.3587, -3.4551, -3.4372,
     & -3.3159, -3.2639, -3.2140, -3.1651, -3.2213, -3.1894,
     & -3.2064, -3.1551, -3.0445, -3.1150, -3.1280, -3.0576,
     & -3.0192, -3.1220, -3.0332, -2.9907, -3.0214, -3.0436,
     & -2.9773, -2.9428, -2.9924, -2.9613, -2.9800, -2.9834,
     & -2.9652, -2.9946, -2.9930, -2.9576, -2.9474, -2.9246,
     & -2.9504, -2.9383, -2.9214, -2.9252, -2.8767, -2.8883/
c PAIR cF NC
      data (pmf( 3, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.2464, -0.7396, -1.9440, -2.6914, -3.0432,
     & -3.3046, -3.3497, -3.2420, -3.1692, -3.2040, -3.2656,
     & -3.2711, -3.2296, -3.2484, -3.2151, -3.1590, -3.2118,
     & -3.2097, -3.1628, -3.1723, -3.0283, -3.0526, -3.0648,
     & -3.0976, -3.0901, -3.0762, -3.0942, -3.0853, -3.0608,
     & -3.0164, -2.9913, -3.0361, -3.0380, -2.9963, -2.9662,
     & -3.0138, -2.9493, -2.9451, -2.9942, -2.9231, -2.9326,
     & -2.9006, -2.8863, -2.9286, -2.9056, -2.9216, -2.9043/
c PAIR cF NP
      data (pmf( 3,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.9185, -1.9515, -2.2338, -3.0171,
     & -3.4446, -3.3158, -3.3784, -3.3253, -3.3102, -3.2934,
     & -3.3814, -3.2828, -3.3154, -3.3607, -3.2136, -3.2742,
     & -3.1935, -3.1825, -3.0589, -3.0407, -3.0199, -2.9758,
     & -2.9325, -2.9658, -2.9886, -2.9620, -2.8959, -3.0203,
     & -3.0016, -2.9552, -2.9659, -2.9719, -2.9686, -2.9461,
     & -3.0126, -2.9890, -2.9906, -3.0081, -2.9794, -2.9039,
     & -2.9555, -2.9808, -2.9772, -3.0059, -2.9525, -2.9467/
c PAIR cF NA
      data (pmf( 3,11,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR cF ND
      data (pmf( 3,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.2188, -3.1288, -3.1401, -3.3155,
     & -3.5199, -3.5808, -3.7018, -3.4085, -3.5051, -3.4706,
     & -3.2515, -3.4506, -3.4321, -3.2623, -3.1265, -3.1959,
     & -3.2769, -3.1931, -3.2427, -3.1789, -3.0328, -3.1099,
     & -3.0493, -3.1001, -3.0974, -2.9807, -2.9726, -3.0170,
     & -3.1327, -2.9973, -3.0435, -2.9965, -2.8941, -3.0056,
     & -2.9263, -2.9194, -2.9033, -2.8712, -2.7551, -2.8971,
     & -2.9257, -2.9069, -2.9013, -2.9877, -2.8643, -2.9560/
c PAIR cF NR
      data (pmf( 3,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.6518, -1.1874, -2.2834, -2.8653, -3.3122,
     & -3.4029, -3.4991, -3.5045, -3.4374, -3.4446, -3.3406,
     & -3.3062, -3.2486, -3.2358, -3.1842, -3.1804, -3.1521,
     & -3.1481, -3.0783, -3.0921, -3.0599, -2.9643, -3.0302,
     & -2.9818, -2.9421, -2.9356, -2.9380, -2.9541, -2.9925,
     & -2.9889, -2.9755, -2.9884, -3.0009, -2.9974, -3.0204,
     & -2.9617, -2.9250, -3.0002, -3.0088, -3.0037, -2.9848,
     & -2.9441, -2.9631, -2.9718, -2.9667, -2.9621, -2.9178/
c PAIR cF N0
      data (pmf( 3,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR cF NS
      data (pmf( 3,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.3491, -1.8982, -2.5150, -3.4171,
     & -3.3504, -3.4148, -3.3667, -3.4280, -3.3944, -3.1272,
     & -3.2357, -3.3101, -3.0684, -3.0359, -3.0836, -3.1579,
     & -3.1558, -3.0734, -3.0770, -3.0714, -3.0232, -3.0409,
     & -2.9753, -3.0403, -3.0791, -2.9007, -2.8950, -2.8645,
     & -3.0248, -3.1461, -2.9491, -2.9602, -2.9335, -2.9636,
     & -2.9965, -2.8929, -3.0023, -2.9866, -2.9895, -3.0314,
     & -2.9827, -2.9982, -2.9536, -2.9632, -2.9495, -2.9655/
c PAIR cF OC
      data (pmf( 3,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.3143,  0.0000,  0.0000,  0.0000,  0.0000, -0.6419,
     &  0.0000, -1.0765, -1.3893, -2.6673, -3.0889, -3.3440,
     & -3.2232, -3.2402, -3.2348, -3.2744, -3.3439, -3.4653,
     & -3.4329, -3.3058, -3.2252, -3.2590, -3.3118, -3.3497,
     & -3.2290, -3.1896, -3.1600, -3.1762, -3.1081, -3.1040,
     & -3.0963, -3.0599, -3.0532, -2.9948, -3.0268, -3.0186,
     & -2.9989, -2.9949, -2.9594, -2.9377, -2.9170, -2.9119,
     & -2.9370, -2.9033, -2.8926, -2.9251, -2.8883, -2.8982,
     & -2.9192, -2.9121, -2.8833, -2.9412, -2.9526, -2.9585/
c PAIR cF OA
      data (pmf( 3,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.9243, -1.4753, -2.4730, -2.9636, -3.1438,
     & -3.4101, -3.4678, -3.2845, -3.3861, -3.2786, -3.3824,
     & -3.2666, -3.2944, -3.2926, -3.2264, -3.2179, -3.2216,
     & -3.2550, -3.1626, -3.0673, -3.0695, -3.1218, -3.1382,
     & -3.1203, -3.1357, -3.0406, -3.0432, -3.0558, -3.0065,
     & -2.9617, -2.9551, -2.9853, -2.9450, -2.9585, -2.9205,
     & -2.9488, -2.9573, -2.9213, -2.9324, -2.8834, -2.9574,
     & -2.9269, -2.9179, -2.9386, -2.9525, -2.9365, -2.9546/
c PAIR cF OE
      data (pmf( 3,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.0586, -1.3388, -2.2759, -2.8738, -3.0061,
     & -3.1705, -3.2627, -3.3287, -3.3258, -3.4120, -3.4071,
     & -3.4415, -3.4550, -3.4307, -3.4073, -3.3667, -3.3600,
     & -3.2982, -3.1841, -3.1409, -3.0796, -3.1651, -3.1823,
     & -3.0951, -3.0136, -2.9790, -3.0374, -3.0249, -3.0235,
     & -3.0042, -2.9588, -2.9674, -2.9278, -2.9248, -2.9080,
     & -2.9793, -2.9777, -2.9704, -2.8957, -2.8808, -2.9200,
     & -2.9303, -2.9226, -2.9194, -2.9060, -2.9157, -2.8598/
c PAIR cF OR
      data (pmf( 3,19,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -2.8664, -3.0522,
     & -3.2212, -3.2509, -2.9441, -3.1708, -2.6702, -3.1976,
     & -3.1148, -2.9619, -2.8863, -2.9497, -3.1657, -3.1691,
     & -3.1695, -2.6626, -2.4570, -2.6031, -2.5615, -2.7380,
     & -2.6262, -2.8678, -2.7969, -2.8703, -2.8906, -3.0048,
     & -3.1812, -3.0589, -3.1946, -3.0199, -3.0090, -3.1144,
     & -3.2165, -3.0638, -3.0757, -3.0400, -3.0399, -3.0176,
     & -3.1421, -2.8721, -3.0407, -3.0939, -2.9694, -2.9284/
c PAIR cF OS
      data (pmf( 3,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.8763, -1.4278, -1.8576, -2.4467, -2.9691,
     & -2.9303, -2.8220, -2.7082, -2.9241, -2.8993, -2.9468,
     & -3.1154, -2.9299, -2.8769, -2.9153, -2.9958, -3.0344,
     & -2.9807, -2.9394, -2.9152, -2.9399, -2.9227, -2.9750,
     & -3.0039, -2.9642, -2.9502, -2.9834, -2.9800, -3.0046,
     & -2.9460, -2.9782, -2.9843, -2.9664, -2.9849, -3.0110,
     & -2.9990, -3.0109, -3.0303, -3.0175, -3.0492, -3.0393,
     & -3.1011, -3.0671, -3.0653, -3.0474, -3.0610, -3.0849/
c PAIR cF OD
      data (pmf( 3,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.4436, -1.2910, -2.3345, -2.5897, -2.9610,
     & -3.1227, -3.2541, -3.3045, -3.2960, -3.3301, -3.3753,
     & -3.4194, -3.3792, -3.3709, -3.3330, -3.2856, -3.2823,
     & -3.2402, -3.1964, -3.1696, -3.1396, -3.1384, -3.0626,
     & -3.0951, -3.0823, -3.0714, -3.0113, -3.0052, -2.9997,
     & -2.9375, -2.9625, -2.9521, -2.9719, -2.8949, -2.9300,
     & -2.9102, -2.9131, -2.8971, -2.8979, -2.9279, -2.9347,
     & -2.9362, -2.9236, -2.9323, -2.9511, -2.9597, -2.9931/
c PAIR cF P
      data (pmf( 3,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -0.8349, -1.3876, -1.5949,
     & -2.3524, -2.4159, -2.5669, -2.8445, -2.9187, -2.9860,
     & -2.8812, -2.8604, -2.7347, -2.9023, -2.9887, -3.0485,
     & -3.0404, -3.0016, -2.9372, -2.9578, -2.9919, -3.0908,
     & -2.9499, -2.9072, -2.9481, -3.0347, -2.9756, -2.9950,
     & -2.9525, -3.0347, -2.9595, -2.9907, -2.9571, -3.0233,
     & -3.0607, -3.0406, -2.9465, -2.9634, -3.0376, -3.0678,
     & -3.0652, -3.0434, -3.0560, -3.1298, -3.1213, -3.1593/
c PAIR cF SA
      data (pmf( 3,23,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR cF SD
      data (pmf( 3,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR cF HL
      data (pmf( 3,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -0.9585, -1.6692,
     & -1.3854, -2.7118, -3.1670, -3.3460, -3.3490, -3.4280,
     & -3.2948, -3.3526, -3.2806, -3.3045, -3.3119, -3.2100,
     & -3.2772, -3.2510, -3.2247, -3.2012, -3.2597, -3.2435,
     & -3.1956, -3.2269, -3.1115, -3.1346, -3.0732, -3.1616,
     & -3.0719, -3.0774, -3.0735, -3.0173, -3.0681, -3.0608,
     & -3.0528, -2.9947, -2.9723, -2.9822, -2.9887, -2.9730,
     & -2.9247, -2.9261, -2.9316, -2.9531, -2.8919, -2.9015,
     & -2.9163, -2.9336, -2.9957, -2.9287, -2.9127, -2.9010/
c PAIR cF F
      data (pmf( 3,26,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.8534, -1.9460, -3.3706,
     & -3.2390, -3.4298, -3.0821, -3.2480, -3.5153, -3.5317,
     & -3.3976, -3.2572, -3.0837, -3.3087, -3.3468, -3.3373,
     & -3.2919, -3.1844, -3.1924, -3.2516, -3.2236, -3.2352,
     & -3.1294, -3.2165, -3.0889, -3.0560, -3.0900, -3.1418,
     & -3.0453, -3.0883, -3.0900, -2.9763, -3.0207, -2.9056,
     & -2.9754, -2.8544, -2.8553, -2.7791, -2.9329, -2.7757,
     & -2.8399, -2.8515, -2.8637, -2.7744, -2.8481, -2.7555/
c PAIR cP CF
      data (pmf( 4, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -0.9646,  0.0000,
     & -1.1494, -1.6888, -2.4522, -2.8299, -2.9220, -2.8830,
     & -3.2209, -3.3184, -3.4254, -3.4609, -3.4728, -3.4301,
     & -3.3028, -3.2656, -3.3054, -3.3237, -3.3242, -3.2376,
     & -3.2489, -3.1930, -3.1600, -3.1926, -3.1456, -3.1049,
     & -3.1199, -3.1472, -3.1208, -3.0841, -3.0393, -3.0155,
     & -2.9888, -3.0532, -3.0221, -2.9836, -2.9777, -2.9354,
     & -2.9342, -2.9156, -2.8749, -2.8948, -2.8838, -2.9111,
     & -2.9071, -2.8887, -2.9081, -2.8831, -2.8250, -2.8739/
c PAIR cP CP
      data (pmf( 4, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.9379, -1.3522, -2.4462, -2.9476, -2.8437, -3.0732,
     & -3.3204, -3.4443, -3.5749, -3.5997, -3.5243, -3.5515,
     & -3.5865, -3.4996, -3.4509, -3.3752, -3.3439, -3.3043,
     & -3.3127, -3.2236, -3.1908, -3.1660, -3.1064, -3.1363,
     & -3.0803, -3.0523, -3.1007, -3.0396, -3.0501, -3.0074,
     & -2.9661, -2.9510, -2.9596, -2.9967, -2.9232, -2.9010,
     & -2.9422, -2.9194, -2.9280, -2.9224, -2.8885, -2.8683,
     & -2.8435, -2.8370, -2.8448, -2.8330, -2.8594, -2.8595/
c PAIR cP cF
      data (pmf( 4, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.7590, -1.0529, -2.7554, -2.9980, -3.2719, -3.2761,
     & -3.3890, -3.3748, -3.4703, -3.3970, -3.3055, -3.3263,
     & -3.3790, -3.2984, -3.3783, -3.4423, -3.3613, -3.2660,
     & -3.3083, -3.2689, -3.1848, -3.2005, -3.1611, -3.0672,
     & -3.0069, -3.0784, -3.0633, -2.9641, -2.9901, -3.0475,
     & -2.9925, -3.0254, -2.9937, -2.9970, -2.9620, -2.9429,
     & -2.8903, -2.8817, -2.9278, -2.9331, -2.8881, -2.8985,
     & -2.8767, -2.8958, -2.9233, -2.8774, -2.8800, -2.8919/
c PAIR cP cP
      data (pmf( 4, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.9824,  0.0000, -2.5567, -3.0058, -3.0697, -3.3552,
     & -3.5189, -3.4622, -3.5396, -3.5694, -3.4692, -3.4437,
     & -3.4511, -3.3285, -3.3168, -3.3192, -3.3042, -3.3097,
     & -3.3267, -3.2208, -3.1338, -3.1151, -3.1354, -3.1485,
     & -3.1630, -3.0154, -2.9859, -3.0432, -3.0380, -2.9853,
     & -2.9967, -2.9242, -2.9707, -3.0177, -2.9414, -2.9347,
     & -2.9552, -2.9303, -2.8951, -2.8540, -2.8463, -2.8970,
     & -2.9083, -2.8690, -2.8793, -2.9436, -2.8977, -2.8849/
c PAIR cP C3
      data (pmf( 4, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.7078,
     & -1.9873, -1.8669, -1.3455, -2.4667, -2.6065, -3.1723,
     & -3.4717, -3.5081, -3.6473, -3.4944, -3.3768, -3.4443,
     & -3.4759, -3.5445, -3.4591, -3.3642, -3.3164, -3.2880,
     & -3.2116, -3.1747, -3.0875, -3.1592, -3.1770, -3.0601,
     & -3.0145, -2.9927, -3.0353, -3.0549, -2.9867, -2.9543,
     & -2.9615, -2.9172, -2.9035, -2.8110, -2.9341, -2.9491,
     & -2.9602, -2.9885, -3.0167, -2.9359, -2.9780, -2.9088,
     & -2.8869, -3.0121, -2.8661, -2.9437, -2.8576, -2.8338/
c PAIR cP CW
      data (pmf( 4, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.5530,
     &  0.0000, -1.9387, -2.4690, -2.9373, -2.6956, -3.0283,
     & -3.3142, -3.2912, -3.5301, -3.6121, -3.5614, -3.5179,
     & -3.5554, -3.4371, -3.4166, -3.3651, -3.2790, -3.3459,
     & -3.2493, -3.1639, -3.1617, -3.1931, -3.1040, -3.1404,
     & -3.1651, -3.0348, -3.0350, -2.9444, -3.0492, -3.0290,
     & -3.0442, -2.9526, -2.9603, -2.9569, -2.9969, -3.0042,
     & -2.8672, -2.9353, -2.9966, -2.9160, -2.9310, -2.9069,
     & -2.7883, -2.8640, -2.8887, -2.9015, -2.8193, -2.8362/
c PAIR cP CO
      data (pmf( 4, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.4022, -1.2928, -2.8621, -3.1011, -3.1341, -3.1971,
     & -3.5841, -3.5845, -3.7211, -3.8383, -3.8137, -3.6552,
     & -3.4899, -3.4518, -3.3533, -3.4314, -3.3596, -3.3033,
     & -3.2905, -3.1837, -3.1104, -3.0849, -2.9868, -2.9469,
     & -3.0063, -2.9308, -2.9977, -2.9605, -3.0016, -2.8881,
     & -2.9099, -2.8907, -3.0261, -2.9051, -2.8779, -2.9505,
     & -2.9980, -2.8574, -2.9178, -2.9185, -2.8176, -2.8393,
     & -2.9225, -2.8312, -2.8907, -2.9253, -2.8359, -2.8759/
c PAIR cP CN
      data (pmf( 4, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.2066,
     &  0.0000, -1.9624, -1.8523, -2.5669, -2.6011, -3.2439,
     & -3.3630, -3.4831, -3.4483, -3.3694, -3.4163, -3.2973,
     & -3.4094, -3.4730, -3.3539, -3.3664, -3.2406, -3.1951,
     & -3.1252, -3.0507, -3.0736, -3.1021, -3.1066, -3.0963,
     & -3.1594, -3.1320, -3.0917, -2.9887, -2.9658, -2.9926,
     & -3.0539, -3.0423, -3.0787, -2.9784, -2.9024, -2.8999,
     & -2.9308, -2.8187, -2.8662, -2.9377, -2.9325, -2.9718,
     & -2.9550, -2.9158, -2.8488, -2.8908, -2.9242, -2.9235/
c PAIR cP NC
      data (pmf( 4, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.6212, -1.9255, -2.3683, -2.5537, -2.8759, -3.2823,
     & -3.3522, -3.6300, -3.4902, -3.4174, -3.3971, -3.3013,
     & -3.3534, -3.3576, -3.2062, -3.0883, -3.1951, -3.2726,
     & -3.2749, -3.1529, -3.0197, -2.9922, -3.0715, -3.1595,
     & -3.0653, -2.9137, -2.9344, -3.0781, -3.0778, -3.0463,
     & -3.0008, -2.9282, -2.9418, -3.0297, -2.9999, -2.9824,
     & -3.0163, -3.0631, -2.9993, -2.9545, -2.9685, -2.9501,
     & -2.8627, -2.8540, -2.8557, -2.8855, -2.8870, -2.8968/
c PAIR cP NP
      data (pmf( 4,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.5372, -2.1612, -2.7737, -3.0466,
     & -2.8775, -3.2544, -3.2532, -3.2828, -3.1434, -3.1531,
     & -3.4050, -3.3780, -3.4273, -3.3039, -3.1080, -3.0932,
     & -2.9157, -3.2831, -3.1558, -3.2056, -3.0572, -2.9794,
     & -3.1664, -3.1547, -3.1311, -3.0920, -3.1761, -3.1755,
     & -3.0286, -2.9629, -2.9637, -2.9637, -2.9656, -3.0300,
     & -2.9173, -2.9625, -2.9240, -2.8061, -2.9019, -2.9472,
     & -2.9030, -2.8950, -2.9014, -2.9747, -2.9241, -2.9177/
c PAIR cP NA
      data (pmf( 4,11,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR cP ND
      data (pmf( 4,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.6450, -3.2681, -3.2380, -3.5517,
     & -3.2814, -3.5279, -3.5068, -3.4167, -3.2673, -3.4945,
     & -3.3065, -3.4536, -3.1826, -3.3429, -3.4262, -3.3036,
     & -3.3214, -3.0695, -2.9683, -2.9820, -3.2511, -3.1233,
     & -2.9805, -2.8642, -3.0681, -3.1148, -3.1210, -3.1193,
     & -2.9535, -2.8568, -3.2057, -2.9620, -2.9207, -2.9763,
     & -2.9134, -2.9203, -3.0055, -2.8816, -2.7894, -2.9486,
     & -2.9908, -2.9988, -2.8670, -2.9727, -2.8893, -2.7403/
c PAIR cP NR
      data (pmf( 4,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.6074, -3.4272, -3.4080, -3.4047,
     & -3.5070, -3.2776, -3.4190, -3.5129, -3.5703, -3.4121,
     & -3.3943, -3.4664, -3.2331, -3.1434, -3.1131, -3.2172,
     & -3.2417, -3.1568, -3.1613, -3.0004, -3.0071, -3.0645,
     & -3.0361, -3.0392, -2.9930, -3.0185, -3.0079, -3.0042,
     & -3.0316, -2.9591, -2.9678, -2.9690, -2.9205, -2.9849,
     & -3.0007, -3.0172, -2.9786, -2.9797, -2.8948, -2.8245,
     & -2.8997, -2.8773, -2.9056, -2.9879, -2.9505, -2.8842/
c PAIR cP N0
      data (pmf( 4,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR cP NS
      data (pmf( 4,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.0770, -2.9138, -3.1551, -3.5558, -3.1794, -3.1302,
     & -3.4917, -3.7436, -3.8525, -3.8995, -3.7009, -3.7116,
     & -3.7024, -3.4153, -3.4447, -3.6830, -3.2884, -3.4730,
     & -3.4109, -3.2667, -3.1273, -3.3740, -3.1094, -3.1284,
     & -2.9594, -2.9566, -2.8732, -3.1477, -3.1062, -3.1487,
     & -3.0952, -3.0365, -2.9883, -2.9237, -2.8859, -2.8871,
     & -2.6809, -2.9215, -2.7502, -2.7497, -2.7789, -2.7669,
     & -2.4437, -2.6520, -2.7317, -2.7929, -2.6809, -2.8281/
c PAIR cP OC
      data (pmf( 4,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.1257,
     & -1.8330, -1.8618, -3.0620, -3.3872, -3.7142, -3.7569,
     & -3.7312, -3.6893, -3.5502, -3.5420, -3.5037, -3.4201,
     & -3.4306, -3.4140, -3.4014, -3.4076, -3.2648, -3.1756,
     & -3.1562, -3.1587, -3.0440, -3.0149, -3.0949, -3.0241,
     & -2.9405, -2.9739, -3.0049, -2.9004, -2.9143, -2.9070,
     & -3.0068, -2.9910, -2.9286, -2.9186, -2.9408, -2.9509,
     & -2.9639, -2.9335, -2.9204, -2.9303, -2.9809, -2.8981,
     & -2.9597, -2.9144, -2.8307, -2.7794, -2.8820, -2.8823/
c PAIR cP OA
      data (pmf( 4,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -1.4370,  0.0000,
     & -1.6268, -1.5228, -2.8666, -3.6244, -3.7915, -3.7311,
     & -3.7119, -3.4323, -3.4325, -3.3654, -3.3741, -3.3470,
     & -3.4574, -3.4683, -3.3472, -3.1578, -3.0461, -3.0587,
     & -3.1682, -3.0965, -3.1714, -3.0757, -3.1174, -3.1285,
     & -3.1155, -3.1755, -3.2318, -3.0442, -3.0905, -3.0749,
     & -2.9333, -2.9647, -2.9541, -3.0683, -3.0238, -3.0147,
     & -3.0289, -2.9317, -2.9088, -2.8884, -2.9687, -2.8848,
     & -2.7738, -2.8832, -2.8118, -2.7871, -2.8311, -2.7703/
c PAIR cP OE
      data (pmf( 4,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.8019, -1.6630, -2.6823, -2.6411, -2.9734, -3.1297,
     & -3.0826, -3.3584, -3.3567, -3.3153, -3.3724, -3.4718,
     & -3.3311, -3.4927, -3.5068, -3.3899, -3.3054, -3.3162,
     & -3.3479, -3.2915, -3.3245, -3.2250, -3.2320, -3.1877,
     & -3.1100, -3.0646, -3.0054, -2.9718, -2.9820, -2.9834,
     & -2.9732, -2.9295, -2.9325, -3.0340, -3.0387, -2.9768,
     & -2.8525, -2.8779, -2.9594, -2.8576, -2.8950, -2.9241,
     & -2.8675, -2.9264, -2.8483, -2.8375, -2.8771, -2.8879/
c PAIR cP OR
      data (pmf( 4,19,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR cP OS
      data (pmf( 4,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.6976, -1.5438, -2.6438, -3.1059, -3.4076, -3.4910,
     & -3.5177, -3.4203, -3.2318, -3.2757, -3.2664, -3.3357,
     & -3.3380, -3.4081, -3.2982, -3.2459, -3.2548, -3.2018,
     & -3.1240, -2.9743, -2.9097, -2.9718, -3.0587, -3.0170,
     & -2.9585, -2.9431, -3.0001, -3.0426, -2.9839, -2.9908,
     & -2.9223, -2.9514, -2.9257, -2.9837, -2.9424, -2.9942,
     & -2.9865, -2.9899, -3.0243, -2.9684, -2.9712, -2.9861,
     & -2.9655, -2.9619, -2.9592, -2.9292, -2.9413, -2.9026/
c PAIR cP OD
      data (pmf( 4,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.4899, -1.9212, -2.6657, -3.1158, -3.3569, -3.6241,
     & -3.6095, -3.5647, -3.4299, -3.5186, -3.4330, -3.4733,
     & -3.4570, -3.4157, -3.4568, -3.4157, -3.2938, -3.2797,
     & -3.3022, -3.2056, -3.1969, -3.1796, -3.1394, -3.0353,
     & -3.0895, -3.1235, -3.0962, -3.0852, -3.0735, -3.0244,
     & -2.9501, -2.9813, -3.0366, -2.9427, -2.9546, -3.0026,
     & -2.9691, -2.9059, -2.8962, -2.8814, -2.8875, -2.8528,
     & -2.8158, -2.8183, -2.8045, -2.8362, -2.8357, -2.8116/
c PAIR cP P
      data (pmf( 4,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -1.9028,  0.0000,
     &  0.0000,  0.0000, -2.8072, -2.8607, -2.9235, -2.5164,
     & -2.5820, -3.0286, -3.3925, -3.7538, -3.7070, -3.7777,
     & -3.4632, -3.1917, -3.0865, -3.1183, -3.3172, -3.0902,
     & -3.0487, -3.0853, -3.0271, -2.7492, -2.9485, -2.9663,
     & -2.8485, -2.8388, -2.9319, -3.0120, -3.0560, -2.9956,
     & -2.9814, -2.9825, -2.9669, -3.0133, -2.9437, -2.9264,
     & -3.0514, -3.0821, -2.9592, -3.0569, -2.9994, -2.9416,
     & -2.9386, -2.9933, -2.9278, -2.9589, -2.9117, -3.0568/
c PAIR cP SA
      data (pmf( 4,23,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR cP SD
      data (pmf( 4,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR cP HL
      data (pmf( 4,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.4856,
     & -2.3212, -2.7317, -2.9036, -3.3694, -3.2815, -3.4114,
     & -3.1763, -3.3244, -3.5176, -3.3217, -3.3286, -3.2589,
     & -3.1968, -3.1930, -3.1525, -3.2192, -3.1474, -3.1080,
     & -3.1505, -3.1531, -3.1254, -3.1113, -3.0960, -2.9879,
     & -2.9655, -3.0585, -3.0712, -3.0866, -3.0598, -3.0373,
     & -3.0344, -3.0754, -2.9957, -3.0528, -2.9867, -2.9596,
     & -2.9285, -2.9512, -3.0002, -2.9510, -2.9192, -2.8632,
     & -2.9033, -2.9749, -2.9549, -2.9509, -2.9306, -2.9646/
c PAIR cP F
      data (pmf( 4,26,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CO CF
      data (pmf( 5, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.5001,
     & -2.0314, -2.6330, -2.7705, -2.6748, -2.9736, -2.8845,
     & -2.8199, -3.1581, -3.2337, -3.3714, -3.3817, -3.2938,
     & -3.0847, -3.1138, -3.1145, -3.1535, -3.1276, -3.0700,
     & -3.1834, -3.1396, -2.9073, -3.0818, -3.0773, -3.0769,
     & -3.0752, -3.0270, -3.1586, -3.0512, -3.0602, -3.0498,
     & -3.0910, -3.0723, -3.0658, -3.0676, -3.0441, -2.9719,
     & -2.9111, -2.8147, -2.9059, -2.9658, -2.9872, -3.0013,
     & -2.9312, -2.8809, -2.9089, -2.9725, -2.9785, -2.9430/
c PAIR CO CP
      data (pmf( 5, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.5473,
     & -2.0593, -2.8203, -3.1535, -3.5283, -3.6870, -3.4253,
     & -3.5132, -3.7308, -3.8068, -3.7444, -3.7106, -3.5671,
     & -3.5461, -3.4333, -3.4528, -3.3396, -3.3458, -3.2331,
     & -3.2104, -3.2099, -3.2404, -3.1817, -3.1286, -3.0246,
     & -2.9452, -2.9827, -3.0404, -2.9944, -2.9602, -2.9562,
     & -2.9650, -2.9052, -2.9534, -2.9725, -2.9613, -2.9030,
     & -2.9138, -2.8746, -2.9275, -2.9667, -2.8944, -2.8572,
     & -2.8619, -2.8705, -2.8234, -2.8394, -2.8546, -2.8311/
c PAIR CO cF
      data (pmf( 5, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.8491, -2.8203, -2.9070, -3.3733, -3.1098,
     & -2.9602, -2.9571, -2.8331, -2.8908, -2.9414, -2.8386,
     & -2.9217, -3.0092, -3.0481, -3.0747, -3.1589, -3.1703,
     & -2.9967, -2.9669, -3.1190, -3.0885, -3.1063, -3.0221,
     & -3.0305, -3.1360, -3.0639, -3.2072, -3.1896, -3.1209,
     & -3.0503, -3.0945, -3.0885, -3.1153, -3.0420, -3.0246,
     & -3.0839, -3.0388, -3.0901, -3.0804, -2.9529, -2.9338,
     & -2.8868, -2.8936, -2.8662, -2.8497, -2.8513, -2.7976/
c PAIR CO cP
      data (pmf( 5, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -1.7838,  0.0000,
     & -1.9462, -2.4652, -2.5147, -2.6413, -3.1324, -2.8367,
     & -3.0216, -3.0356, -3.3841, -3.4303, -3.6778, -3.5743,
     & -3.2488, -3.1869, -3.1826, -3.2390, -3.4222, -3.3943,
     & -3.3155, -3.1739, -2.9120, -2.7720, -2.8111, -2.8357,
     & -2.9944, -2.9393, -2.9610, -3.1353, -3.1513, -3.0471,
     & -2.9825, -2.9633, -3.0693, -3.0764, -2.9899, -2.8880,
     & -2.8841, -2.9445, -2.9845, -3.0407, -2.9912, -2.9273,
     & -2.9346, -2.9547, -2.9677, -2.8757, -2.9108, -2.8627/
c PAIR CO C3
      data (pmf( 5, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.6359, -3.0641, -2.9585, -3.1363, -2.9590,
     & -2.7747, -2.7770, -2.6926, -3.2602, -3.3494, -3.2930,
     & -3.1316, -3.1683, -3.2077, -3.1714, -3.1359, -3.1435,
     & -2.9940, -2.7622, -2.8539, -2.8606, -3.0442, -2.8433,
     & -2.7764, -2.9080, -2.8477, -3.0207, -3.0040, -3.0952,
     & -3.0824, -2.9240, -3.0768, -2.9616, -2.9731, -2.9998,
     & -3.0484, -3.0563, -2.9544, -3.0330, -2.9432, -3.0675,
     & -3.0983, -3.0160, -2.9707, -3.0199, -2.9378, -3.0835/
c PAIR CO CW
      data (pmf( 5, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.6318, -2.5055, -3.2965, -3.4790, -3.2927, -3.1240,
     & -3.4796, -3.5192, -3.6478, -3.5943, -3.7461, -3.6492,
     & -3.5270, -3.3079, -3.3267, -3.3346, -3.1696, -3.3110,
     & -3.3918, -3.3491, -3.1147, -3.0587, -3.2020, -3.0790,
     & -3.0087, -2.9567, -3.1561, -3.1191, -2.9890, -3.1314,
     & -3.0120, -2.9532, -3.0335, -2.9301, -2.8922, -2.8405,
     & -2.9841, -2.9719, -2.9393, -2.9356, -2.9135, -2.8132,
     & -2.8436, -2.8309, -2.8756, -2.7655, -2.8062, -2.8349/
c PAIR CO CO
      data (pmf( 5, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.5937, -3.0700, -3.1055, -3.5313, -3.4408,
     & -2.9459, -3.2479, -3.4175, -3.5292, -3.4919, -3.2162,
     & -3.0606, -3.3232, -3.4699, -3.1537, -3.1174, -3.1253,
     & -3.0148, -3.1496, -3.0527, -3.0557, -2.9440, -3.0620,
     & -3.0905, -2.9831, -2.9485, -2.9005, -3.0143, -3.0949,
     & -3.1222, -3.0687, -3.0257, -2.9077, -2.9506, -2.8648,
     & -2.9023, -2.9019, -3.0261, -3.1082, -2.9385, -2.9551,
     & -2.9554, -3.0415, -2.9671, -2.8189, -2.8492, -2.8758/
c PAIR CO CN
      data (pmf( 5, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -2.3307, -2.2234,
     &  0.0000, -2.7971, -2.6870, -2.4139, -2.3165, -2.5215,
     & -3.2179, -3.5849, -4.0341, -3.5857, -3.6784, -3.6668,
     & -3.4456, -3.3278, -3.1228, -3.1208, -3.0348, -3.0251,
     & -3.0568, -3.2255, -3.1458, -3.1039, -3.2298, -3.1196,
     & -3.0962, -3.1093, -3.0445, -3.0811, -3.0268, -2.9663,
     & -2.7947, -2.8152, -2.8644, -3.0527, -2.9814, -2.8818,
     & -2.7668, -2.7767, -2.8553, -2.9025, -2.9030, -3.0798,
     & -2.9849, -2.9240, -2.8177, -2.8281, -3.0052, -2.9825/
c PAIR CO NC
      data (pmf( 5, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.0692, -2.3735, -2.8163, -3.4106, -3.6748, -4.3050,
     & -4.2391, -4.0291, -3.4232, -3.0549, -3.2430, -3.0716,
     & -3.1820, -3.0846, -3.1228, -2.8674, -2.9242, -3.1278,
     & -3.1154, -3.1322, -3.1612, -2.6834, -2.8421, -3.0188,
     & -2.9370, -3.2166, -3.2467, -3.1435, -3.1914, -3.0853,
     & -3.1858, -2.8464, -2.7864, -2.8413, -3.0531, -3.0711,
     & -3.0575, -3.0406, -2.8639, -2.7395, -2.7525, -2.6322,
     & -2.6581, -2.8214, -2.8762, -2.9036, -2.9857, -2.9592/
c PAIR CO NP
      data (pmf( 5,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -3.0625,  0.0000, -2.3417, -2.6139, -2.7265, -3.2387,
     & -3.4146, -2.6500, -2.8881, -3.4556, -4.0051, -3.8127,
     & -3.5738, -3.2871, -3.1262, -3.2002, -3.3085, -3.0735,
     & -2.9047, -2.9467, -3.1559, -2.9994, -3.1138, -2.9949,
     & -3.0560, -3.1312, -3.0237, -3.1388, -3.1577, -3.0865,
     & -3.0128, -2.7172, -2.9344, -2.9441, -3.0016, -3.0487,
     & -2.9986, -2.9971, -2.9087, -2.8360, -2.8644, -2.9099,
     & -2.9105, -2.9367, -2.9843, -2.9341, -2.8787, -2.9050/
c PAIR CO NA
      data (pmf( 5,11,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CO ND
      data (pmf( 5,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.5504, -2.4262, -2.9584, -3.4976, -4.0746,
     & -3.7489, -3.2396, -3.2995, -3.1622, -2.7556, -3.2285,
     & -3.0786, -3.0712, -3.5249, -3.5272, -3.4330, -3.4406,
     & -3.2678, -3.0075, -3.2494, -2.9517, -2.8331, -2.7943,
     & -3.2713, -3.3376, -2.9738, -3.0553, -2.7769, -2.8723,
     & -2.9765, -3.0470, -3.1937, -3.0459, -3.0838, -3.0053,
     & -3.1590, -2.8799, -2.9564, -2.7945, -2.9195, -2.4741,
     & -3.0607, -2.8715, -2.6785, -2.9178, -3.0372, -2.7826/
c PAIR CO NR
      data (pmf( 5,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.0794,
     &  0.0000, -2.4775, -2.8648, -2.6682, -3.4560, -4.1611,
     & -4.0125, -3.6178, -3.2658, -3.3724, -3.0853, -3.1257,
     & -3.1496, -3.0548, -3.1759, -3.4814, -3.3278, -3.0499,
     & -3.1020, -3.1742, -3.2874, -3.1353, -3.1011, -3.2425,
     & -3.0907, -3.1277, -2.9435, -2.9344, -3.1217, -3.1431,
     & -3.1353, -2.9163, -2.8839, -2.8718, -2.9855, -2.9913,
     & -2.9291, -2.9569, -2.9237, -2.8793, -2.7624, -2.8054,
     & -2.7804, -2.8886, -2.8675, -2.9012, -2.9843, -2.9547/
c PAIR CO N0
      data (pmf( 5,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CO NS
      data (pmf( 5,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -3.1572, -3.7908, -3.3615, -3.3996, -4.2322,
     & -3.7351, -3.1347, -3.4637, -3.2556, -3.5106, -2.6994,
     & -3.2270, -3.5693, -3.3173, -3.4147, -2.5071, -2.9259,
     & -3.2035, -3.0204, -2.9370, -3.0756, -2.8578, -2.8212,
     & -3.2352, -3.2739, -2.7569, -3.2414, -3.1765, -3.0100,
     & -3.2680, -3.2774, -2.9440, -2.7182, -2.6613, -2.3640,
     & -2.7945, -3.1556, -3.1766, -2.7438, -2.6237, -2.6506,
     & -2.9061, -3.0154, -2.8324, -2.7461, -3.1434, -2.8017/
c PAIR CO OC
      data (pmf( 5,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.8006,
     &  0.0000, -3.1004, -2.9057, -3.2998, -3.7784, -3.7460,
     & -3.2651, -2.9967, -3.2855, -3.3273, -3.0241, -3.1712,
     & -2.7585, -3.1248, -3.1183, -3.1534, -3.0365, -3.2242,
     & -3.2055, -3.0108, -3.1585, -3.0768, -2.9794, -3.0919,
     & -3.0333, -2.9975, -2.9806, -3.0825, -3.0884, -3.1983,
     & -3.1244, -3.0134, -3.0937, -2.9572, -2.9170, -3.0009,
     & -2.9578, -3.0164, -3.0110, -2.9892, -2.9462, -2.9839,
     & -2.8824, -2.8174, -2.8646, -2.8899, -2.8717, -2.9319/
c PAIR CO OA
      data (pmf( 5,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.8401, -1.7361, -2.6960, -3.0136, -3.4059, -3.7759,
     & -3.6792, -3.4202, -3.4760, -3.5539, -3.5360, -3.6012,
     & -3.4550, -3.1535, -3.1676, -3.0066, -3.1524, -3.1388,
     & -3.0303, -2.9423, -3.0368, -2.9707, -3.2238, -3.0235,
     & -3.1412, -3.1295, -3.0814, -3.0684, -3.0357, -3.1233,
     & -3.0684, -2.9084, -2.9876, -3.0771, -2.9416, -2.7951,
     & -2.7507, -2.8142, -2.9465, -2.8009, -2.9080, -2.8942,
     & -2.8861, -2.8994, -3.0550, -2.8922, -2.9608, -3.0070/
c PAIR CO OE
      data (pmf( 5,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.4016,
     & -2.6599, -3.0616, -3.2102, -3.2780, -3.3486, -3.3714,
     & -3.1607, -3.0579, -2.9608, -3.2397, -3.5264, -3.3848,
     & -3.4242, -3.1992, -3.3219, -3.3401, -3.3164, -3.3391,
     & -3.3904, -3.4409, -3.2780, -3.3020, -3.0792, -2.9819,
     & -3.0005, -3.0146, -2.8674, -3.0069, -2.8995, -2.9664,
     & -2.9490, -2.9928, -2.9305, -2.9880, -3.0925, -2.9416,
     & -2.9586, -3.0178, -2.8663, -2.8079, -2.8369, -2.8926,
     & -2.7928, -2.9067, -2.9806, -3.0353, -2.9322, -2.8870/
c PAIR CO OR
      data (pmf( 5,19,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CO OS
      data (pmf( 5,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.6116,
     & -2.2395, -2.5451, -3.2338, -3.1963, -3.3566, -3.3465,
     & -3.0201, -3.1191, -3.1702, -3.0721, -3.1195, -3.0799,
     & -3.1561, -2.9539, -3.1095, -3.1113, -3.1436, -3.0264,
     & -3.0746, -3.1534, -3.1170, -3.1013, -3.0700, -3.0728,
     & -3.1107, -3.1644, -3.1134, -3.0943, -3.1127, -3.1352,
     & -3.1353, -3.0791, -3.1244, -3.1064, -3.0168, -2.9555,
     & -2.9312, -3.0320, -2.9913, -2.9452, -2.9244, -2.8597,
     & -2.8650, -2.9131, -2.9823, -2.8927, -2.7863, -2.8500/
c PAIR CO OD
      data (pmf( 5,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.5412, -2.6946, -3.2957, -4.1679, -4.3252, -4.3293,
     & -4.0353, -3.7609, -3.5712, -3.3064, -3.3114, -3.0932,
     & -3.2204, -3.3649, -3.2711, -3.3756, -3.4868, -3.3520,
     & -3.2828, -3.2247, -3.0189, -3.0306, -3.0745, -2.9788,
     & -3.0122, -2.9779, -3.0127, -3.0445, -3.0347, -3.0068,
     & -3.0155, -3.0162, -2.9658, -2.8604, -2.8232, -2.8806,
     & -2.8906, -2.8524, -2.8927, -2.8379, -2.8156, -2.8886,
     & -2.8944, -2.8661, -2.8542, -2.8220, -2.7914, -2.7955/
c PAIR CO P
      data (pmf( 5,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.9168, -2.9571, -2.8559, -2.8392, -1.9272,
     & -2.4862, -2.9801, -3.1856, -3.1728, -3.0103, -3.2719,
     & -3.2664, -2.9105, -3.3027, -3.1589, -3.0631, -3.0869,
     & -3.2001, -3.1244, -3.0721, -3.1027, -3.1222, -3.2372,
     & -3.0455, -3.0731, -3.2227, -3.1751, -3.1889, -3.1998,
     & -2.9537, -3.1889, -3.1375, -3.0643, -2.9465, -3.0864,
     & -3.0440, -3.0416, -3.1191, -3.0207, -2.9083, -2.8824,
     & -3.0159, -2.7750, -2.9237, -2.7713, -2.7322, -2.7491/
c PAIR CO SA
      data (pmf( 5,23,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CO SD
      data (pmf( 5,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CO HL
      data (pmf( 5,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.4182,  0.0000, -3.0804,
     & -3.1833, -3.8229, -3.5791, -3.6305, -3.6895, -3.6298,
     & -3.5061, -3.4651, -3.4683, -3.2653, -3.1498, -3.2690,
     & -3.3320, -3.2477, -3.2053, -3.2519, -3.3115, -3.2379,
     & -3.2989, -3.1765, -3.2187, -3.2807, -3.2521, -3.2182,
     & -3.1852, -3.1671, -3.0072, -3.1073, -3.0868, -3.0065,
     & -3.0297, -2.9991, -2.9270, -2.9127, -2.8898, -2.9045,
     & -2.9397, -2.9359, -2.9523, -2.9274, -2.9118, -2.8626,
     & -2.8978, -2.8515, -2.8641, -2.7921, -2.8472, -2.7737/
c PAIR CO F
      data (pmf( 5,26,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CN CF
      data (pmf( 6, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.4581,  0.0000, -1.2476, -1.9699, -2.2886, -2.9594,
     & -3.2735, -3.1156, -3.1561, -3.1885, -3.2128, -3.1767,
     & -3.1290, -3.0709, -3.2377, -3.3981, -3.2809, -3.1913,
     & -3.3252, -3.3573, -3.1420, -3.1539, -3.0050, -2.9868,
     & -3.0831, -3.0464, -3.1417, -3.0090, -3.0112, -3.0310,
     & -3.1239, -3.0740, -3.0381, -3.0028, -2.9552, -3.0393,
     & -3.0405, -2.9634, -3.0150, -2.9321, -2.8598, -2.9487,
     & -2.9652, -2.8743, -2.9040, -2.8807, -2.8747, -2.8456/
c PAIR CN CP
      data (pmf( 6, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.9695, -1.8590, -2.2958, -2.5588,
     & -2.8892, -3.2576, -3.3472, -3.4075, -3.5715, -3.6860,
     & -3.6520, -3.5532, -3.4231, -3.4542, -3.4215, -3.3265,
     & -3.3131, -3.2176, -3.1884, -3.2580, -3.2356, -3.1268,
     & -3.1190, -3.0584, -3.0944, -3.1079, -3.1050, -3.0234,
     & -3.0041, -2.9478, -2.9465, -2.8968, -2.9304, -2.8771,
     & -2.8336, -2.8901, -2.9161, -2.8934, -2.8952, -2.8650,
     & -2.8491, -2.8462, -2.8586, -2.8616, -2.8301, -2.8327/
c PAIR CN cF
      data (pmf( 6, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.5443, -1.4292, -1.7305, -2.0383, -1.5339, -2.9501,
     & -3.2555, -3.3202, -3.2094, -3.1392, -3.2996, -3.2714,
     & -3.3573, -3.2237, -3.2506, -3.3194, -3.2742, -3.2246,
     & -3.3487, -3.3091, -3.0647, -3.0134, -3.1070, -3.0659,
     & -3.0794, -3.1938, -3.2529, -3.1005, -2.9697, -2.8201,
     & -2.9438, -3.0391, -2.9282, -2.9299, -2.9865, -3.0532,
     & -2.9625, -3.0183, -2.9767, -2.9128, -2.9274, -2.9422,
     & -2.9652, -2.9178, -2.9292, -2.8471, -2.8722, -2.9046/
c PAIR CN cP
      data (pmf( 6, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.5647,  0.0000, -1.7243, -2.3508, -2.2450, -3.1196,
     & -3.5156, -3.6210, -3.7093, -3.5779, -3.5516, -3.5108,
     & -3.4460, -3.2587, -3.2737, -3.3183, -3.1955, -3.0859,
     & -2.9820, -3.1018, -3.1193, -3.0747, -3.0766, -2.9974,
     & -2.9918, -3.0805, -3.0152, -2.9228, -3.0054, -3.0412,
     & -3.0034, -3.0370, -3.0747, -3.0679, -3.0173, -2.9435,
     & -2.9116, -2.9782, -2.9830, -3.0079, -2.9648, -2.9195,
     & -2.8713, -2.9096, -2.9270, -2.8476, -2.8681, -2.8782/
c PAIR CN C3
      data (pmf( 6, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.8631,
     &  0.0000, -2.6132, -2.5007, -2.6343, -1.8866, -2.1970,
     & -2.9911, -3.0209, -3.3180, -3.3403, -3.4749, -3.3283,
     & -3.1964, -3.2692, -3.4018, -3.4653, -3.3010, -3.1339,
     & -3.2051, -3.1906, -3.1462, -3.1328, -3.2459, -3.0740,
     & -3.1235, -3.1245, -3.0318, -3.0901, -3.0278, -3.0789,
     & -2.9781, -2.9955, -3.0711, -2.9902, -2.9376, -2.8974,
     & -2.8773, -2.9787, -2.8749, -2.9430, -2.7975, -2.8343,
     & -2.8712, -3.0915, -2.9354, -3.0292, -2.8551, -2.8868/
c PAIR CN CW
      data (pmf( 6, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.0739,  0.0000,  0.0000, -2.5373, -2.0245, -2.5661,
     & -3.0801, -3.3001, -3.3846, -3.1821, -3.7268, -3.7608,
     & -3.6458, -3.6737, -3.4721, -3.5509, -3.3046, -3.2427,
     & -3.2371, -3.2409, -3.2154, -3.2698, -3.1848, -2.9653,
     & -3.0942, -3.0568, -3.1342, -3.1465, -3.1281, -3.0701,
     & -2.9448, -2.8679, -3.0181, -3.0070, -2.9181, -2.7876,
     & -2.8822, -2.7993, -2.8947, -2.8704, -2.9286, -2.9739,
     & -2.7950, -2.7985, -2.8924, -2.8743, -2.8246, -2.8155/
c PAIR CN CO
      data (pmf( 6, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.0036,  0.0000, -1.7919,  0.0000, -2.6622, -2.8730,
     & -3.1001, -3.9466, -4.1714, -3.9564, -4.0279, -4.0487,
     & -3.8989, -3.3628, -3.3791, -3.4811, -3.2330, -3.1694,
     & -3.1419, -3.2351, -3.3356, -3.4134, -3.3842, -3.1271,
     & -3.0912, -3.0977, -3.1541, -2.9697, -3.0434, -2.8939,
     & -2.8499, -2.8953, -3.0171, -2.9515, -2.7561, -2.7026,
     & -2.7927, -2.7153, -2.5545, -2.6543, -2.6600, -2.6012,
     & -2.7185, -2.7540, -2.7461, -2.8257, -2.8814, -2.7779/
c PAIR CN CN
      data (pmf( 6, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.0873,  0.0000, -2.2828, -1.7765, -3.1462,
     & -3.3563, -3.2983, -3.5336, -3.3064, -3.5349, -3.5430,
     & -3.6082, -3.5487, -3.4104, -3.4275, -3.1902, -3.0723,
     & -2.8836, -2.8880, -2.9162, -2.9841, -2.9678, -2.9088,
     & -2.9359, -3.0258, -3.0596, -2.8805, -2.9652, -2.9106,
     & -2.9252, -2.9890, -2.8805, -2.8081, -2.8318, -2.8454,
     & -2.9684, -2.9636, -3.0510, -3.0189, -2.9431, -3.0746,
     & -3.0597, -3.0592, -3.0077, -2.8808, -2.9596, -2.9071/
c PAIR CN NC
      data (pmf( 6, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.9843, -2.8419, -2.8628, -3.2308,
     & -3.4012, -3.5987, -3.6542, -3.6765, -3.4807, -3.3812,
     & -3.3481, -3.3445, -3.2979, -3.2517, -3.2810, -3.1731,
     & -3.1146, -3.1779, -3.1004, -2.9840, -2.7456, -2.8555,
     & -3.0042, -2.9998, -2.9851, -2.9520, -2.8367, -2.8977,
     & -2.9297, -2.8907, -3.0142, -2.9658, -2.8802, -2.8052,
     & -2.8031, -2.8800, -2.9650, -2.9836, -3.0136, -2.9051,
     & -2.9203, -3.0894, -3.0458, -3.1038, -2.9773, -2.8938/
c PAIR CN NP
      data (pmf( 6,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.2712, -2.5536, -2.9582,
     & -3.3965, -3.7734, -3.5571, -3.5004, -3.4639, -3.4269,
     & -3.1508, -3.3683, -3.3029, -3.3763, -3.1918, -3.0533,
     & -2.9872, -2.9988, -2.7038, -3.1179, -3.0159, -3.0986,
     & -3.1990, -3.1838, -3.1552, -3.2398, -3.2661, -3.1469,
     & -2.9765, -2.9383, -3.0738, -3.0458, -3.0731, -2.9799,
     & -2.9763, -2.9104, -3.0387, -3.0418, -2.9306, -2.7572,
     & -2.8384, -2.8571, -2.9135, -2.7834, -2.8001, -2.7903/
c PAIR CN NA
      data (pmf( 6,11,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CN ND
      data (pmf( 6,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.3462, -2.8855, -3.0735,
     & -2.0174, -2.8666, -3.2369, -2.8802, -2.9399, -2.3771,
     & -3.0765, -3.6346, -3.6123, -3.5976, -3.3760, -3.3668,
     & -3.1990, -3.2828, -3.1802, -3.1889, -3.4582, -3.4271,
     & -3.1105, -3.0598, -3.1363, -3.1513, -2.9541, -2.9083,
     & -3.0887, -3.0743, -2.9198, -2.8465, -3.0380, -3.0182,
     & -2.9502, -2.9328, -2.9996, -2.9530, -2.7668, -2.9471,
     & -3.0651, -2.8355, -2.7497, -2.8089, -2.6734, -2.8047/
c PAIR CN NR
      data (pmf( 6,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.7224, -2.0257, -2.8135, -3.3599,
     & -3.4955, -3.7592, -3.7622, -3.5649, -3.4021, -3.3524,
     & -3.2777, -3.3594, -3.3825, -3.2894, -3.0728, -3.0398,
     & -3.0455, -2.9770, -3.0476, -2.9375, -2.9515, -2.8392,
     & -2.8952, -2.7581, -2.8973, -2.9353, -2.9308, -3.0297,
     & -3.1380, -2.9640, -2.8526, -2.9673, -2.9788, -2.9065,
     & -2.8733, -3.0085, -2.9916, -3.0554, -3.0152, -3.0388,
     & -2.9456, -2.9880, -2.9795, -2.9629, -2.9975, -2.9467/
c PAIR CN N0
      data (pmf( 6,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CN NS
      data (pmf( 6,15,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CN OC
      data (pmf( 6,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.5791,  0.0000, -2.1980, -3.0861, -3.8530, -4.2579,
     & -4.3791, -3.9972, -3.8996, -3.5437, -3.4508, -3.5002,
     & -3.1921, -3.3546, -3.4122, -3.6688, -3.5153, -3.4039,
     & -3.1149, -3.1235, -3.1488, -3.1228, -3.1445, -3.1180,
     & -3.1416, -3.1152, -3.0677, -3.1423, -3.0723, -3.0385,
     & -2.9729, -2.9216, -2.9383, -2.8200, -2.7995, -2.6549,
     & -2.7011, -2.7143, -2.7096, -2.6308, -2.7632, -2.8553,
     & -2.7864, -2.6930, -2.7369, -2.8134, -2.7631, -2.8583/
c PAIR CN OA
      data (pmf( 6,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.9470,  0.0000, -2.6953, -3.6204, -3.8880, -3.6823,
     & -3.5047, -3.4991, -3.6739, -3.6775, -2.9301, -3.0538,
     & -3.1460, -3.2707, -3.3896, -3.3908, -3.3631, -3.1194,
     & -2.9856, -3.0576, -3.0268, -3.1954, -3.0659, -2.9128,
     & -3.1663, -3.0945, -3.2551, -2.9712, -2.9442, -2.9450,
     & -2.8460, -2.9480, -2.8373, -2.9316, -2.9125, -2.8892,
     & -2.8123, -2.7477, -2.8885, -2.8326, -3.0344, -3.0230,
     & -2.9419, -2.9217, -2.9750, -2.9747, -2.9568, -2.9823/
c PAIR CN OE
      data (pmf( 6,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.4021, -3.4785, -3.8446,
     & -3.9137, -3.7284, -3.5691, -3.6709, -3.2575, -2.6968,
     & -2.8678, -2.9201, -3.1778, -3.2954, -3.4121, -3.3460,
     & -3.5003, -3.4714, -3.4079, -3.4318, -3.2152, -3.1507,
     & -3.0927, -3.1874, -3.1603, -3.0735, -3.1089, -3.0304,
     & -2.9123, -2.8559, -2.8653, -2.9872, -2.9786, -2.8841,
     & -2.9697, -2.8823, -2.8946, -2.8839, -2.9447, -2.9148,
     & -2.8244, -2.7511, -2.7876, -2.7277, -2.8576, -2.7941/
c PAIR CN OR
      data (pmf( 6,19,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CN OS
      data (pmf( 6,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.7550, -1.9530, -2.6440, -3.2317, -3.8652, -4.0482,
     & -4.0822, -3.9863, -3.8120, -3.5970, -3.5346, -3.5609,
     & -3.4596, -3.5391, -3.5475, -3.5218, -3.5750, -3.4717,
     & -3.3552, -3.2611, -3.0607, -3.0384, -3.0697, -3.0163,
     & -2.9672, -3.0612, -3.0939, -2.9581, -3.0776, -3.0301,
     & -2.9747, -2.8822, -2.8974, -2.9077, -2.8384, -2.8060,
     & -2.8714, -2.8994, -2.8449, -2.8059, -2.7650, -2.8179,
     & -2.7569, -2.7391, -2.7102, -2.7443, -2.7694, -2.7349/
c PAIR CN OD
      data (pmf( 6,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.4461,
     &  0.0000, -1.8676, -2.5833, -2.9871, -3.4200, -3.7191,
     & -3.8113, -3.7876, -3.6006, -3.5662, -3.4499, -3.3958,
     & -3.4079, -3.3386, -3.2434, -3.3361, -3.3482, -3.3035,
     & -3.2273, -3.1879, -3.2682, -3.1464, -3.1580, -3.1251,
     & -3.1202, -3.1037, -3.1558, -3.1348, -3.0721, -3.1011,
     & -3.0510, -2.9038, -2.9501, -2.8783, -2.8928, -2.9063,
     & -2.9281, -2.9368, -2.8948, -2.8643, -2.9165, -2.8931,
     & -2.7976, -2.7882, -2.7972, -2.8785, -2.8330, -2.7802/
c PAIR CN P
      data (pmf( 6,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -2.6661, -2.6630,
     & -3.3480, -3.5024, -4.2796, -4.2908, -4.1973, -4.0965,
     & -3.9612, -3.7902, -3.5043, -3.2146, -3.0025, -2.9170,
     & -2.8573, -3.0555, -3.2839, -3.3077, -3.2530, -3.1938,
     & -3.0358, -3.2048, -3.0764, -2.9691, -2.9578, -2.8951,
     & -2.8187, -2.8114, -2.9484, -3.0461, -2.9302, -2.9538,
     & -2.7880, -2.8523, -2.7869, -2.7674, -2.7340, -2.6676,
     & -2.7440, -2.7476, -2.6335, -2.7058, -2.6646, -2.8767/
c PAIR CN SA
      data (pmf( 6,23,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CN SD
      data (pmf( 6,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR CN HL
      data (pmf( 6,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.2156, -2.3488, -2.9600, -3.1441, -3.3649, -3.5225,
     & -3.2538, -3.4507, -3.4435, -3.4885, -3.5219, -3.3280,
     & -3.2610, -3.2494, -3.2585, -3.4453, -3.2371, -3.0397,
     & -3.1826, -3.1903, -3.2228, -3.1244, -3.1368, -3.1372,
     & -3.1972, -3.1166, -3.0243, -2.9534, -3.0517, -3.0901,
     & -3.0588, -2.9504, -2.9732, -3.0215, -3.0324, -2.9918,
     & -2.9546, -2.9955, -2.9590, -2.8905, -2.9557, -2.8795,
     & -2.9727, -2.9265, -2.8392, -2.8345, -2.7505, -2.8563/
c PAIR CN F
      data (pmf( 6,26,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR NC CF
      data (pmf( 7, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.8083,  0.0000,  0.0000,  0.0000, -2.2540,
     & -1.5965, -2.4935, -2.7733, -2.9075, -3.1109, -3.2725,
     & -3.3800, -3.3517, -3.2993, -3.3506, -3.3122, -3.2748,
     & -3.3242, -3.4296, -3.3886, -3.2451, -3.2796, -3.3382,
     & -3.3586, -3.1980, -3.0925, -3.1199, -3.0533, -3.0858,
     & -3.1104, -3.1602, -3.1124, -3.0788, -3.0447, -3.0286,
     & -3.0350, -3.0064, -2.9676, -2.9767, -2.9842, -2.9285,
     & -2.9472, -2.8941, -2.8744, -2.9635, -2.9345, -2.8743,
     & -2.9139, -2.8853, -2.8652, -2.9075, -2.8775, -2.8250/
c PAIR NC CP
      data (pmf( 7, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.4382,
     & -1.3020, -2.7112, -2.7950, -2.9734, -3.1274, -3.4342,
     & -3.8265, -3.7067, -3.7967, -3.5765, -3.4688, -3.5475,
     & -3.5900, -3.5715, -3.4648, -3.3933, -3.3523, -3.3149,
     & -3.2621, -3.2453, -3.2398, -3.1967, -3.1858, -3.1627,
     & -3.1085, -3.1182, -3.1170, -3.0909, -3.0520, -3.0040,
     & -2.9655, -2.9550, -2.9231, -2.8877, -2.9214, -2.8528,
     & -2.8580, -2.8179, -2.8873, -2.8482, -2.8942, -2.8464,
     & -2.8351, -2.8622, -2.8233, -2.8108, -2.7823, -2.8444/
c PAIR NC cF
      data (pmf( 7, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.2823,  0.0000,  0.0000,  0.0000, -1.3695,
     & -2.1953, -2.5454, -2.8193, -2.8442, -3.0854, -3.2687,
     & -3.2921, -3.3985, -3.4027, -3.2940, -3.2988, -3.2833,
     & -3.3290, -3.3602, -3.3093, -3.3909, -3.2958, -3.2572,
     & -3.1640, -3.0790, -3.2299, -3.1655, -3.2344, -3.2261,
     & -3.1862, -3.0730, -2.9871, -3.0367, -3.0614, -3.0279,
     & -3.0276, -2.9931, -2.9838, -2.9572, -2.9452, -2.9237,
     & -2.9153, -2.9423, -2.9156, -2.9783, -2.9317, -2.8607,
     & -2.9238, -2.9451, -2.8906, -2.8303, -2.8703, -2.8739/
c PAIR NC cP
      data (pmf( 7, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.9323,  0.0000,  0.0000,  0.0000, -1.8374,
     & -1.9373, -2.8653, -2.8360, -2.6947, -3.2271, -3.5809,
     & -3.4607, -3.6449, -3.6385, -3.4672, -3.4815, -3.4220,
     & -3.3857, -3.3561, -3.2887, -3.2282, -3.1541, -3.1289,
     & -3.1227, -3.1435, -3.1594, -3.1586, -3.1289, -3.1052,
     & -3.0427, -3.0391, -3.0620, -3.1002, -3.0950, -2.9495,
     & -2.9742, -3.0115, -2.9910, -3.0371, -3.0176, -2.9959,
     & -2.9357, -2.9595, -2.9627, -2.9278, -2.8852, -2.8998,
     & -2.8263, -2.8759, -2.9480, -2.8812, -2.8899, -2.8125/
c PAIR NC C3
      data (pmf( 7, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -3.0894,  0.0000,  0.0000,  0.0000, -2.1766,
     & -2.8650, -2.7446, -3.1106, -2.5265, -3.4054, -3.3605,
     & -3.5787, -3.5612, -3.1665, -3.4096, -3.2292, -3.3485,
     & -3.3672, -3.4540, -3.3163, -3.2230, -3.3270, -3.3729,
     & -3.3365, -3.1657, -3.1721, -3.0977, -3.1950, -3.1515,
     & -3.0491, -3.1483, -3.1464, -3.0635, -3.0620, -2.9980,
     & -2.8703, -2.9059, -2.9462, -2.9760, -2.9939, -2.9986,
     & -2.8746, -2.9593, -2.9694, -2.9440, -2.9003, -2.9657,
     & -2.8227, -2.7865, -2.8473, -2.9452, -2.8754, -2.8923/
c PAIR NC CW
      data (pmf( 7, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.4392,  0.0000,  0.0000, -2.4515,  0.0000,
     & -2.2081, -2.4907, -3.1120, -2.9578, -3.3474, -3.3956,
     & -3.7407, -3.8482, -3.9085, -3.6357, -3.5616, -3.4912,
     & -3.5718, -3.4218, -3.5148, -3.2895, -3.3199, -3.3673,
     & -3.3351, -3.2728, -3.3364, -3.1379, -3.2244, -3.1215,
     & -3.1370, -3.1351, -3.1347, -3.1133, -3.0502, -3.0049,
     & -2.9989, -3.0651, -2.9487, -2.9350, -3.0009, -2.8448,
     & -2.9211, -2.7848, -2.8645, -2.8531, -2.8474, -2.8533,
     & -2.7801, -2.7629, -2.7944, -2.7337, -2.8052, -2.7852/
c PAIR NC CO
      data (pmf( 7, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.3403,  0.0000,  0.0000,  0.0000, -1.8364,
     & -1.7189, -2.0185, -2.8036, -3.3334, -3.9687, -4.3733,
     & -4.4381, -4.1515, -3.9480, -3.7580, -3.6559, -3.4948,
     & -3.4501, -3.2549, -3.2941, -3.4125, -3.3658, -3.4056,
     & -3.4029, -3.3387, -3.3229, -3.3291, -3.2050, -3.1000,
     & -3.1602, -3.1954, -3.1961, -3.1130, -2.9534, -3.0375,
     & -2.9875, -2.8768, -2.6698, -2.7137, -2.7171, -2.8337,
     & -2.7798, -2.7921, -2.7301, -2.7434, -2.7254, -2.6377,
     & -2.6910, -2.6490, -2.6631, -2.5939, -2.6434, -2.7826/
c PAIR NC CN
      data (pmf( 7, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.9360, -2.4665, -2.5262, -2.6621, -3.0645, -3.3103,
     & -3.3646, -3.5643, -3.5039, -3.3795, -3.5798, -3.2849,
     & -3.2247, -3.3596, -3.4390, -3.3133, -3.3469, -3.2674,
     & -3.2999, -3.0430, -3.0073, -3.0897, -2.9828, -2.9716,
     & -3.0126, -3.0301, -3.0760, -2.9093, -2.8187, -2.8939,
     & -2.7980, -2.8959, -2.8484, -2.9669, -2.8890, -2.9653,
     & -2.9652, -2.9273, -2.9658, -2.9238, -3.0898, -3.0426,
     & -3.0479, -2.9065, -2.9196, -2.8807, -2.9662, -2.9418/
c PAIR NC NC
      data (pmf( 7, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.0191,
     & -2.7245, -2.9500, -3.1188, -3.1692, -3.2307, -3.4820,
     & -3.3602, -3.3714, -3.6033, -3.6055, -3.5262, -3.4496,
     & -3.3471, -3.3485, -3.4402, -3.3315, -3.1876, -3.1147,
     & -3.1738, -2.9979, -3.0555, -3.0505, -2.9519, -3.0954,
     & -3.1059, -2.9687, -3.0139, -2.9575, -2.9008, -2.8373,
     & -2.8192, -2.8942, -2.8326, -2.8159, -2.9224, -2.9298,
     & -2.8688, -2.9525, -2.9203, -2.9826, -2.9866, -2.9360,
     & -2.9425, -2.9835, -2.9970, -2.9789, -2.9611, -3.0157/
c PAIR NC NP
      data (pmf( 7,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -3.0741, -2.9374, -3.2759, -3.1399,
     & -3.5864, -3.5093, -3.5050, -3.6368, -3.4453, -3.5525,
     & -3.5850, -3.2282, -3.2242, -3.0816, -2.9084, -3.0194,
     & -2.9376, -2.9301, -3.2545, -3.3228, -3.1414, -3.1148,
     & -3.2141, -3.2699, -3.2231, -3.1344, -3.0937, -3.0620,
     & -3.0907, -3.0407, -3.0539, -3.0475, -3.0272, -3.0139,
     & -2.9072, -2.8116, -2.8375, -2.8783, -2.9050, -2.8827,
     & -2.9470, -2.9903, -2.8959, -2.7711, -2.7744, -2.7061/
c PAIR NC NA
      data (pmf( 7,11,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR NC ND
      data (pmf( 7,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.8589,  0.0000, -3.2493, -2.4851, -3.0243, -3.0807,
     & -2.8954, -3.0055, -3.2574, -3.8045, -3.8363, -3.2114,
     & -3.5168, -3.5572, -3.6724, -3.7195, -3.5291, -3.2894,
     & -3.2866, -3.2386, -3.1999, -3.1173, -3.0775, -3.2084,
     & -3.2047, -3.1796, -3.2153, -3.2141, -3.1175, -3.1633,
     & -2.8940, -3.0549, -3.0256, -2.9309, -2.8566, -2.9237,
     & -2.9627, -2.9210, -2.9115, -2.7256, -2.7883, -2.8637,
     & -2.9169, -2.8173, -2.7654, -2.7447, -2.6516, -2.6774/
c PAIR NC NR
      data (pmf( 7,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.8231,
     & -2.3416, -2.9877, -3.3595, -3.4352, -3.4315, -3.3585,
     & -3.3532, -3.2498, -3.2644, -3.3029, -3.5343, -3.4514,
     & -3.4648, -3.2949, -3.3453, -3.2975, -3.1456, -3.1057,
     & -3.1626, -3.0164, -3.0144, -2.9782, -2.9984, -2.9413,
     & -3.0359, -2.9539, -2.9280, -3.0136, -2.9702, -3.0320,
     & -3.0348, -2.9653, -2.9623, -2.9801, -3.0735, -3.0079,
     & -3.0236, -2.9985, -3.0319, -3.0123, -2.9939, -2.9046,
     & -2.9283, -2.8829, -2.9656, -2.9255, -2.8953, -2.8475/
c PAIR NC N0
      data (pmf( 7,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR NC NS
      data (pmf( 7,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -3.3812, -2.6203, -2.9237, -3.2335, -3.2717, -3.7425,
     & -3.0906, -3.6142, -3.6488, -3.2592, -3.2916, -2.9016,
     & -3.5080, -3.4414, -3.2907, -3.3403, -3.2258, -2.9809,
     & -3.4359, -3.4999, -3.5360, -3.0854, -2.8444, -2.9599,
     & -2.8058, -3.0099, -2.7986, -2.9859, -3.0309, -2.7049,
     & -3.0216, -2.8001, -2.9476, -2.5600, -2.9423, -2.6691,
     & -2.8248, -2.9105, -3.0595, -2.8328, -2.9539, -2.9660,
     & -2.8069, -2.7985, -3.0815, -3.1230, -3.0067, -2.8929/
c PAIR NC OC
      data (pmf( 7,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.9088,  0.0000,  0.0000,  0.0000, -2.7012,
     & -3.4852, -4.3763, -4.4545, -4.1872, -4.0546, -4.0282,
     & -3.8752, -3.6068, -3.5775, -3.5695, -3.5274, -3.5025,
     & -3.5897, -3.5039, -3.3623, -3.3532, -3.2836, -3.3198,
     & -3.2560, -3.1879, -3.1497, -3.1976, -3.2468, -3.2507,
     & -3.2024, -3.0716, -3.1541, -3.1497, -3.0442, -2.9517,
     & -3.0025, -2.9248, -2.8360, -2.9390, -2.8640, -2.7190,
     & -2.6836, -2.6809, -2.6496, -2.6907, -2.7310, -2.7650,
     & -2.7317, -2.7534, -2.7089, -2.7326, -2.6974, -2.6350/
c PAIR NC OA
      data (pmf( 7,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.4339,
     & -2.7310, -4.0143, -4.2013, -4.1271, -3.2611, -3.1030,
     & -3.1433, -3.1325, -3.2875, -3.2967, -3.2420, -3.1416,
     & -3.3418, -3.3248, -3.2063, -3.2439, -3.2691, -3.1852,
     & -3.1162, -3.2138, -3.2926, -3.2532, -3.3333, -3.2466,
     & -3.0617, -2.8989, -2.9110, -2.9485, -2.9553, -2.9796,
     & -2.9936, -2.9656, -2.9394, -2.9765, -2.9910, -2.8256,
     & -2.8753, -2.8054, -2.8139, -2.9744, -2.9122, -2.8963,
     & -2.8510, -2.8729, -2.8518, -2.9490, -2.9820, -2.9739/
c PAIR NC OE
      data (pmf( 7,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.7628, -3.7432, -4.2219, -4.1226, -3.7122, -3.4424,
     & -3.2877, -3.4162, -3.4213, -3.0834, -3.1811, -3.3617,
     & -3.4121, -3.6445, -3.5643, -3.4801, -3.3709, -3.4120,
     & -3.4051, -3.3823, -3.2804, -3.3453, -3.2978, -3.1909,
     & -3.2577, -3.1249, -3.1020, -3.0277, -2.9674, -3.0105,
     & -2.9425, -2.9078, -2.8981, -2.9742, -2.9188, -2.8855,
     & -2.9690, -2.9629, -2.9516, -2.8570, -2.7942, -2.7625,
     & -2.7438, -2.7646, -2.7578, -2.7483, -2.7909, -2.7481/
c PAIR NC OR
      data (pmf( 7,19,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR NC OS
      data (pmf( 7,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.1860, -1.8859, -2.5178,
     & -3.5953, -4.0600, -4.2380, -3.9802, -3.7955, -3.7797,
     & -3.7242, -3.7520, -3.7300, -3.6252, -3.6141, -3.6901,
     & -3.7234, -3.8244, -3.6319, -3.4112, -3.2715, -3.1534,
     & -3.0726, -3.1623, -3.0562, -3.1326, -3.2142, -3.1625,
     & -3.1271, -2.9706, -2.9843, -2.9658, -2.9132, -2.9482,
     & -3.0194, -2.8991, -2.9374, -2.8914, -2.8709, -2.8313,
     & -2.8193, -2.8001, -2.7596, -2.7765, -2.7992, -2.7554,
     & -2.7396, -2.7573, -2.7792, -2.7487, -2.7654, -2.7162/
c PAIR NC OD
      data (pmf( 7,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.6892,  0.0000,  0.0000,  0.0000, -1.1854,
     & -3.0546, -3.7882, -4.0700, -3.8467, -3.6140, -3.5828,
     & -3.5889, -3.5148, -3.4680, -3.4801, -3.3953, -3.3880,
     & -3.4664, -3.4226, -3.3900, -3.3812, -3.3631, -3.3302,
     & -3.2170, -3.2150, -3.2252, -3.0961, -3.1557, -3.1783,
     & -3.1972, -3.0957, -3.0515, -3.0151, -3.0689, -3.0441,
     & -3.0299, -2.9609, -2.9612, -2.9590, -2.9605, -2.8657,
     & -2.8193, -2.8099, -2.7949, -2.8635, -2.9035, -2.8394,
     & -2.8437, -2.8386, -2.8570, -2.7851, -2.8393, -2.8213/
c PAIR NC P
      data (pmf( 7,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.6630, -2.9137, -3.4200, -3.7541, -4.0741,
     & -4.4469, -4.5029, -4.1809, -3.8963, -3.7107, -3.4403,
     & -3.4451, -3.3729, -3.4006, -3.6185, -3.4789, -3.3559,
     & -3.3092, -3.2873, -3.2233, -3.1839, -3.0641, -2.9329,
     & -3.0309, -2.9632, -3.0563, -3.0453, -3.0276, -3.0005,
     & -3.0662, -3.0311, -3.0167, -2.8460, -2.7980, -2.8162,
     & -2.8255, -2.7058, -2.7630, -2.7324, -2.8303, -2.7083,
     & -2.7284, -2.6668, -2.6488, -2.7212, -2.7461, -2.6117/
c PAIR NC SA
      data (pmf( 7,23,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR NC SD
      data (pmf( 7,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR NC HL
      data (pmf( 7,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.1308,  0.0000,  0.0000, -2.6838, -3.0417,
     & -2.9277, -3.1442, -3.1828, -3.2594, -3.4692, -3.4728,
     & -3.4600, -3.6075, -3.4446, -3.3371, -3.3181, -3.3801,
     & -3.3763, -3.4687, -3.3738, -3.3026, -3.3649, -3.3091,
     & -3.2079, -3.2689, -3.2333, -3.1758, -3.1168, -3.1431,
     & -3.1979, -3.1201, -3.0455, -3.1110, -2.9986, -3.0908,
     & -3.0315, -3.0074, -2.9412, -2.9521, -2.9465, -2.9773,
     & -2.8940, -2.8904, -3.0139, -2.9982, -2.9593, -2.8463,
     & -2.8459, -2.8365, -2.7600, -2.8262, -2.8075, -2.7611/
c PAIR NC F
      data (pmf( 7,26,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR ND CF
      data (pmf( 8, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.6522, -0.4265, -0.2790,  0.0000, -0.6780, -1.0707,
     & -1.0326, -1.2540, -1.2310, -1.5801, -1.9407, -2.3146,
     & -2.5474, -2.6551, -2.6774, -2.7046, -2.7645, -2.7925,
     & -2.8153, -2.8543, -2.9000, -2.9036, -2.9622, -2.9604,
     & -2.9317, -2.9695, -2.9937, -3.0004, -2.9820, -2.9889,
     & -3.0208, -3.0158, -3.0162, -3.0169, -3.0355, -3.0212,
     & -3.0500, -3.0384, -3.0682, -3.0592, -3.0618, -3.0628,
     & -3.0692, -3.0551, -3.0440, -3.0524, -3.0614, -3.0370,
     & -3.0295, -3.0406, -3.0094, -2.9995, -2.9821, -2.9771/
c PAIR ND CP
      data (pmf( 8, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.6275,  0.0000,  0.0000,  0.0466, -0.2444, -1.3640,
     & -1.1897, -1.3406, -1.8105, -1.9302, -1.9773, -2.4221,
     & -2.7681, -2.9497, -2.9374, -2.9325, -3.0115, -2.9731,
     & -2.9867, -3.0207, -3.0476, -2.9912, -2.9747, -3.0009,
     & -2.9499, -3.0008, -2.9829, -3.0057, -2.9687, -2.9924,
     & -2.9960, -3.0092, -3.0400, -3.0208, -3.0300, -3.0231,
     & -3.0146, -3.0373, -3.0429, -3.0385, -3.0422, -3.0568,
     & -3.0582, -3.0525, -3.0471, -3.0349, -3.0212, -3.0258,
     & -3.0131, -2.9924, -2.9815, -2.9870, -2.9807, -2.9812/
c PAIR ND cF
      data (pmf( 8, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.3423, -0.4684,  0.0000, -0.1898, -0.8897, -1.9905,
     & -1.4766, -1.2827, -1.6247, -1.5798, -1.9766, -2.3296,
     & -2.5714, -2.6369, -2.6750, -2.6885, -2.7606, -2.8391,
     & -2.9588, -2.9591, -2.9138, -2.8907, -2.9079, -2.8950,
     & -2.9057, -2.9592, -2.9823, -3.0313, -3.0097, -3.0101,
     & -3.0027, -3.0109, -3.0272, -3.0623, -3.0424, -3.0555,
     & -3.0497, -3.0540, -3.0769, -3.0850, -3.0586, -3.0815,
     & -3.0726, -3.0566, -3.0392, -3.0401, -3.0162, -3.0221,
     & -3.0028, -3.0150, -3.0111, -2.9982, -2.9846, -2.9517/
c PAIR ND cP
      data (pmf( 8, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.1811,  0.0000,  0.0000,  0.0000, -0.2070, -0.5087,
     & -1.0175, -1.2461, -1.5073, -1.6961, -2.2271, -2.5777,
     & -2.9209, -3.0788, -2.9917, -2.8772, -2.9402, -3.0054,
     & -2.9801, -2.9513, -2.9384, -2.9573, -2.9842, -3.0541,
     & -3.0532, -3.0258, -3.0270, -3.0088, -3.0049, -2.9919,
     & -3.0379, -3.0471, -3.0912, -3.0895, -3.0714, -3.0567,
     & -3.0497, -3.0480, -3.0497, -3.0049, -3.0563, -3.0404,
     & -3.0462, -3.0434, -3.0266, -3.0196, -3.0088, -2.9796,
     & -2.9907, -2.9752, -2.9883, -2.9792, -2.9686, -2.9424/
c PAIR ND C3
      data (pmf( 8, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.6656,
     &  0.0000,  0.0000,  0.0000,  0.0000, -0.9038, -1.2055,
     & -1.4849, -1.6947, -1.9912, -1.8419, -2.0253, -2.5363,
     & -2.5179, -2.6164, -2.8085, -2.7816, -2.8201, -2.8620,
     & -2.8724, -2.8790, -2.8726, -2.7874, -2.9142, -2.9314,
     & -2.9135, -2.9016, -2.9062, -2.9603, -2.9558, -2.9860,
     & -2.9982, -3.0125, -3.0142, -3.0396, -3.0475, -3.0498,
     & -3.0612, -3.0584, -3.0711, -3.0400, -3.0860, -3.0713,
     & -3.0563, -3.0537, -3.0573, -3.0557, -3.0565, -3.0181,
     & -3.0267, -3.0269, -3.0049, -3.0252, -3.0145, -2.9825/
c PAIR ND CW
      data (pmf( 8, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.6899,  0.0000,  0.0000, -0.8338,  0.0000, -1.6656,
     & -0.8813, -1.2955, -1.8640, -1.7893, -1.8101, -2.6643,
     & -2.7670, -2.9379, -3.0211, -2.8210, -2.8393, -2.9771,
     & -2.9486, -2.9469, -3.0367, -3.0221, -3.0068, -2.9818,
     & -2.9822, -2.9625, -2.9608, -2.9518, -2.9440, -2.9644,
     & -2.9781, -2.9871, -3.0115, -3.0128, -3.0044, -3.0097,
     & -3.0095, -3.0127, -3.0659, -3.0382, -3.0314, -3.0534,
     & -3.0709, -3.0458, -3.0754, -3.0488, -3.0526, -3.0368,
     & -3.0113, -3.0132, -3.0017, -3.0192, -2.9747, -2.9977/
c PAIR ND CO
      data (pmf( 8, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.3249,  0.0000,  0.0000,  0.0000,  0.0000, -0.6524,
     & -1.3528, -0.8345, -1.7380, -1.8629, -2.3690, -2.9355,
     & -2.8920, -3.0839, -2.9347, -2.8393, -2.7225, -2.7879,
     & -2.9289, -2.8888, -2.9452, -3.0909, -2.9922, -2.8770,
     & -2.8794, -2.9447, -2.9708, -2.9607, -3.0045, -3.0267,
     & -2.9866, -2.9512, -2.9740, -3.0237, -3.1087, -3.1048,
     & -3.0461, -3.0092, -3.0107, -3.0153, -3.0153, -3.0168,
     & -3.0300, -3.0411, -3.0775, -3.0552, -3.0337, -3.0497,
     & -3.0377, -3.0075, -3.0107, -2.9928, -2.9993, -2.9958/
c PAIR ND CN
      data (pmf( 8, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.1276,
     & -0.5921, -1.1226, -1.9214, -1.7744, -2.1503, -2.6672,
     & -3.0930, -3.1301, -3.0276, -3.0284, -2.9448, -2.9855,
     & -3.0198, -3.0029, -3.0653, -3.1357, -3.0556, -3.1103,
     & -3.1472, -3.0419, -3.0681, -3.0797, -3.0738, -3.0062,
     & -3.0051, -3.0721, -3.0786, -3.0560, -3.0388, -3.0358,
     & -3.0628, -3.0385, -3.0334, -3.0435, -3.0236, -3.0219,
     & -3.0089, -3.0316, -3.0691, -3.0044, -3.0113, -3.0103,
     & -2.9765, -2.9520, -2.9276, -2.9321, -2.9459, -2.9276/
c PAIR ND NC
      data (pmf( 8, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.7515,  0.0000, -1.0265, -0.8954,  0.0000,  0.0000,
     & -0.9665, -1.1011, -2.1788, -2.3936, -2.4845, -2.8284,
     & -2.9528, -3.1150, -3.0772, -3.1429, -3.1975, -3.0574,
     & -3.0790, -3.2110, -3.1495, -3.0493, -3.0527, -3.0153,
     & -3.0116, -3.0389, -3.0356, -3.0287, -3.0170, -3.0073,
     & -3.0081, -3.0445, -3.0823, -3.0785, -3.0718, -3.0402,
     & -3.0502, -3.0663, -3.0271, -3.0187, -3.0294, -3.0191,
     & -3.0300, -3.0506, -3.0239, -3.0220, -3.0194, -3.0076,
     & -2.9910, -2.9626, -2.9868, -2.9203, -2.9028, -2.8877/
c PAIR ND NP
      data (pmf( 8,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.7300,
     &  0.0000, -0.9771, -1.7779, -1.8397, -2.2703, -2.3950,
     & -2.2258, -2.4587, -2.5228, -2.7442, -2.9248, -3.0697,
     & -3.2023, -3.1091, -3.0260, -3.0677, -3.0578, -3.0807,
     & -3.0363, -3.0003, -2.9597, -3.0562, -3.0243, -3.0603,
     & -3.0672, -3.0202, -3.0069, -3.0536, -3.0521, -3.0249,
     & -3.0488, -3.0636, -3.0756, -3.0441, -3.0373, -3.0642,
     & -3.0467, -3.0507, -3.0328, -2.9906, -2.9768, -3.0099,
     & -3.0232, -2.9998, -2.9608, -2.9602, -2.9520, -2.9098/
c PAIR ND NA
      data (pmf( 8,11,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -3.5727,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.6477,  0.0000,  0.0000, -2.3270,  0.0000,
     &  0.0000,  0.0000, -2.3616, -1.8702, -2.9372, -2.1204,
     & -2.8640, -3.0791, -3.0105, -2.7877, -2.6600, -2.7301,
     & -3.1184, -2.8207, -2.8499, -3.0571, -3.1221, -2.6561,
     & -2.7891, -2.8089, -3.0302, -2.9427, -3.0609, -2.9992,
     & -3.1752, -2.9216, -2.8434, -2.8494, -3.0615, -3.0349,
     & -3.1535, -3.1854, -3.1607, -3.0352, -3.1173, -3.1738,
     & -2.9641, -2.9514, -3.0303, -3.0817, -3.0044, -3.0276/
c PAIR ND ND
      data (pmf( 8,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.1386, -1.0051, -1.9380, -1.9917, -2.2537, -2.9280,
     & -2.9675, -2.7611, -2.8262, -2.8121, -3.0065, -3.0798,
     & -3.1383, -3.2509, -2.8075, -2.9523, -2.9920, -2.9616,
     & -2.9284, -2.9664, -2.9869, -2.9154, -2.8980, -2.9304,
     & -2.9010, -2.9178, -3.0124, -2.9663, -3.0822, -2.9895,
     & -3.0295, -3.0077, -3.0509, -3.0210, -3.0011, -3.0460,
     & -3.0234, -3.0613, -3.0415, -3.0350, -2.9733, -3.0525,
     & -3.0643, -3.0270, -3.1078, -3.0337, -2.9999, -3.0173/
c PAIR ND NR
      data (pmf( 8,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -0.5150,
     & -1.7439, -2.0886, -3.0507, -2.9829, -2.8455, -2.8925,
     & -2.9419, -2.8291, -2.9498, -3.1208, -3.1628, -3.0988,
     & -3.0330, -3.0294, -3.0121, -3.0854, -3.0781, -3.0832,
     & -3.0043, -3.0255, -2.9880, -3.0123, -3.0489, -3.1363,
     & -3.0703, -3.0471, -3.0418, -3.0907, -3.0201, -3.0193,
     & -3.0154, -3.0756, -3.0640, -3.0465, -3.0618, -3.0512,
     & -3.0726, -3.0480, -3.0191, -2.9727, -2.9716, -2.9691,
     & -2.9613, -2.9419, -2.9546, -2.9479, -2.9391, -2.9354/
c PAIR ND N0
      data (pmf( 8,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR ND NS
      data (pmf( 8,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.5000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.5936, -1.4881, -2.4461, -2.6237, -2.7381,
     & -2.8818, -2.9055, -2.9572, -2.8505, -2.8437, -2.9320,
     & -3.0000, -2.8811, -2.9474, -2.8892, -2.9422, -2.8373,
     & -2.9953, -2.9671, -2.9244, -2.7979, -2.9379, -2.9881,
     & -2.9957, -2.9479, -2.9666, -3.1188, -2.9817, -3.1444,
     & -3.0562, -2.9874, -3.0305, -3.0502, -3.0388, -3.0408,
     & -3.0598, -3.0561, -3.0327, -3.1092, -2.9373, -3.0281,
     & -3.0362, -3.0928, -3.0562, -2.9870, -2.9593, -2.9838/
c PAIR ND OC
      data (pmf( 8,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.9247,  0.0000,  0.0000,  0.0000,  0.0000, -0.2522,
     & -1.6986, -2.7438, -3.1194, -2.9986, -2.6612, -2.6107,
     & -2.5468, -2.6947, -2.8064, -2.8411, -2.9616, -2.9522,
     & -3.1382, -3.0501, -2.9326, -2.8827, -2.9000, -2.8993,
     & -2.8910, -2.9383, -2.9516, -3.0278, -3.0265, -3.0452,
     & -3.0406, -3.0380, -3.0344, -3.0447, -3.0216, -3.0034,
     & -2.9745, -3.0012, -3.0283, -3.0815, -3.0455, -3.0173,
     & -2.9985, -3.0030, -3.0489, -3.0321, -3.0057, -3.0360,
     & -3.0328, -3.0281, -3.0155, -3.0014, -2.9992, -3.0024/
c PAIR ND OA
      data (pmf( 8,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.3540,
     & -2.1676, -3.1207, -3.5304, -3.2785, -3.0934, -2.9337,
     & -2.9393, -2.6443, -2.7961, -2.9031, -3.0457, -3.0975,
     & -2.9983, -2.9551, -2.9324, -2.8835, -2.8825, -2.8519,
     & -2.8939, -2.8394, -2.9808, -3.0218, -3.0562, -3.0206,
     & -3.0263, -3.0146, -3.0248, -3.0141, -3.0101, -3.0125,
     & -3.0115, -3.0431, -3.0164, -3.0607, -3.0307, -3.0374,
     & -3.0430, -3.0467, -3.0033, -3.0041, -3.0435, -3.0365,
     & -3.0187, -3.0117, -3.0354, -2.9893, -2.9871, -2.9513/
c PAIR ND OE
      data (pmf( 8,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -1.0903, -1.8009,
     & -1.7818, -2.2912, -2.8362, -2.9679, -2.6747, -2.6785,
     & -2.4706, -2.7493, -2.8587, -2.9607, -2.9435, -3.0157,
     & -3.0779, -2.9925, -2.9490, -2.9298, -3.0187, -3.0584,
     & -3.0248, -3.0646, -3.0869, -3.0415, -2.9860, -3.0025,
     & -2.9983, -3.0257, -2.9647, -3.0044, -3.0555, -3.0699,
     & -3.0243, -3.0084, -2.9990, -3.0234, -3.0429, -3.0922,
     & -3.0699, -3.0625, -3.0239, -3.0249, -3.0251, -2.9923,
     & -2.9700, -2.9740, -2.9927, -2.9957, -2.9732, -2.9723/
c PAIR ND OR
      data (pmf( 8,19,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.9597,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.1301,  0.0000,  0.0000, -2.6782, -2.7609, -2.8558,
     & -2.6853, -2.5138, -3.0897, -3.6386, -3.4705, -3.3363,
     & -3.2435, -3.2306, -3.1726, -2.6137, -2.8182, -2.9182,
     & -2.9569, -2.9249, -3.1007, -3.1267, -3.1789, -3.2774,
     & -3.0606, -3.0803, -3.0342, -3.0722, -3.1251, -2.9768,
     & -2.9257, -2.9603, -2.9792, -2.9422, -2.9504, -3.1029,
     & -3.0564, -2.9978, -2.9506, -2.9190, -2.9964, -2.9447,
     & -2.9457, -3.0258, -3.0177, -2.9688, -2.9674, -2.8888/
c PAIR ND OS
      data (pmf( 8,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.8651,  0.0000,  0.0000, -1.1571, -0.7089, -1.6587,
     & -2.0881, -3.2017, -3.4956, -3.2867, -3.0210, -3.0401,
     & -3.0078, -3.0947, -3.1314, -3.1368, -3.1694, -3.1978,
     & -3.1913, -3.2042, -3.1226, -3.0157, -2.9621, -2.9288,
     & -2.9451, -2.9602, -2.9894, -3.0405, -3.0622, -3.0150,
     & -3.0247, -3.0188, -3.0065, -2.9948, -3.0090, -3.0213,
     & -3.0367, -3.0472, -3.0574, -3.0389, -3.0287, -2.9824,
     & -2.9910, -3.0173, -2.9951, -3.0064, -2.9867, -2.9871,
     & -2.9805, -2.9646, -2.9684, -2.9562, -2.9693, -2.9658/
c PAIR ND OD
      data (pmf( 8,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.3813,  0.0000, -0.1026,  0.0154, -1.0255,
     & -1.6099, -2.5910, -3.0547, -3.0529, -2.9475, -2.8607,
     & -2.8343, -2.9009, -2.9197, -2.8882, -2.9886, -3.0061,
     & -3.0208, -3.0016, -2.9313, -2.9221, -2.9446, -2.9867,
     & -2.9997, -2.9872, -2.9779, -3.0085, -3.0110, -3.0268,
     & -3.0130, -3.0421, -3.0507, -3.0420, -3.0368, -3.0161,
     & -3.0260, -3.0434, -3.0692, -3.0441, -3.0169, -3.0559,
     & -3.0547, -3.0514, -3.0197, -3.0144, -3.0181, -3.0170,
     & -3.0127, -2.9966, -2.9859, -2.9710, -2.9601, -2.9566/
c PAIR ND P
      data (pmf( 8,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.4110,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -0.6147, -1.4491, -1.2096, -1.5174, -2.2210, -3.1345,
     & -3.5423, -3.6498, -3.6093, -3.5137, -3.1682, -3.0103,
     & -2.9635, -2.9593, -2.9834, -3.1449, -3.1386, -3.1348,
     & -3.0750, -2.9743, -2.9415, -2.9174, -2.9403, -2.9227,
     & -2.9369, -2.9908, -3.0389, -3.0172, -3.0589, -3.1158,
     & -3.0808, -3.0180, -3.0196, -3.0000, -3.0107, -3.0036,
     & -3.0035, -2.9809, -3.0343, -2.9898, -2.9804, -2.9930,
     & -2.9885, -2.9645, -3.0158, -3.0164, -2.9764, -2.9725/
c PAIR ND SA
      data (pmf( 8,23,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -3.2335, -3.5351,
     & -2.5923,  0.0000, -2.3674,  0.0000, -2.1691, -2.4812,
     & -2.6287, -2.8429, -3.1066, -2.8790, -2.8830, -2.7324,
     & -2.5725, -2.8075, -2.7436, -2.9912, -2.6712, -2.8808,
     & -2.8350, -3.0634, -2.6731, -2.5464, -2.7327, -3.0036,
     & -3.0996, -3.0785, -2.9669, -3.0638, -3.0428, -3.1664,
     & -3.2575, -3.1580, -3.2209, -3.0160, -3.1270, -3.1774,
     & -3.1797, -3.1120, -3.0494, -3.0174, -2.9999, -3.0852,
     & -2.8622, -2.8682, -2.9735, -2.9528, -2.9253, -2.8980/
c PAIR ND SD
      data (pmf( 8,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR ND HL
      data (pmf( 8,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.0873,  0.0000,  0.0000, -1.0491, -1.6702, -1.7113,
     & -2.0062, -2.0699, -2.4382, -2.6508, -2.6985, -2.9307,
     & -2.9127, -2.8073, -2.9475, -2.9287, -2.9552, -2.8824,
     & -2.9060, -2.9357, -2.9759, -2.9373, -2.9846, -2.9714,
     & -2.9235, -2.9863, -2.9844, -3.0678, -3.0288, -3.0083,
     & -3.0013, -3.0036, -3.0172, -3.0506, -3.0883, -3.0644,
     & -3.0548, -3.0824, -3.0589, -3.0880, -3.0470, -3.0448,
     & -3.0257, -3.0541, -3.0270, -3.0345, -3.0208, -3.0232,
     & -3.0121, -3.0173, -2.9785, -2.9714, -3.0031, -2.9890/
c PAIR ND F
      data (pmf( 8,26,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.5864, -2.6613, -2.4020, -2.0740,
     & -2.3987, -2.4500, -2.2417, -2.9075, -2.8354, -2.8068,
     & -2.7402, -2.2264, -2.5298, -2.9417, -2.8402, -2.9119,
     & -2.8665, -2.6818, -2.9072, -2.9200, -2.9893, -2.9183,
     & -3.0603, -3.0575, -3.0383, -3.1868, -3.1430, -3.0428,
     & -3.1631, -3.0620, -3.0676, -3.0922, -3.1358, -3.0436,
     & -3.0559, -2.9974, -3.0016, -3.0438, -2.9618, -3.0090,
     & -3.1188, -3.1144, -3.0707, -2.8965, -2.9371, -2.9471/
c PAIR NR CF
      data (pmf( 9, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.0475,  0.0000, -3.1274, -3.5954, -3.1006,
     & -2.0749, -2.3751, -2.6309, -2.6776, -2.9054, -3.1322,
     & -3.3462, -3.3234, -3.5238, -3.5413, -3.4206, -3.4201,
     & -3.4295, -3.4816, -3.2996, -3.4527, -3.4812, -3.2943,
     & -3.4488, -3.3009, -3.1701, -3.2404, -3.2905, -3.1925,
     & -3.2182, -3.1542, -3.1827, -3.1235, -3.1110, -3.0716,
     & -3.1589, -2.9367, -2.9498, -2.8670, -2.9127, -2.8871,
     & -2.7589, -2.8582, -2.8365, -2.9405, -2.8495, -2.8429,
     & -2.8940, -2.7155, -2.7323, -2.7896, -2.7619, -2.7665/
c PAIR NR CP
      data (pmf( 9, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.5361, -3.7401, -3.2835,
     & -3.1241, -3.0646, -2.5378, -3.2660, -3.4893, -3.5146,
     & -3.5488, -3.5350, -3.6410, -3.6509, -3.6640, -3.5719,
     & -3.5204, -3.5257, -3.4419, -3.4230, -3.3810, -3.3855,
     & -3.3881, -3.3277, -3.3393, -3.2495, -3.1832, -3.1177,
     & -3.1373, -3.2252, -3.1365, -3.1686, -3.1964, -3.0935,
     & -3.0146, -3.0108, -3.0234, -2.9879, -2.8749, -2.8389,
     & -2.8782, -2.8576, -2.8028, -2.7991, -2.7969, -2.7450,
     & -2.7768, -2.7134, -2.7288, -2.7174, -2.6836, -2.7172/
c PAIR NR cF
      data (pmf( 9, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -3.2232, -3.6690, -3.4632,
     &  0.0000, -2.2186, -2.7077, -2.8277, -2.6098, -2.7046,
     & -3.1531, -3.0250, -3.2430, -3.3914, -3.4408, -3.6015,
     & -3.5651, -3.4956, -3.3675, -3.4431, -3.2361, -3.3823,
     & -3.3411, -3.2944, -3.3938, -3.3007, -3.1668, -3.1051,
     & -3.0533, -3.0159, -3.1155, -3.1825, -3.0488, -2.9066,
     & -2.8189, -2.9380, -2.8936, -2.8493, -2.8325, -2.8318,
     & -2.8741, -2.9212, -2.8790, -2.8610, -2.8845, -2.9015,
     & -2.9852, -2.8957, -2.8834, -2.9355, -2.8704, -2.8528/
c PAIR NR cP
      data (pmf( 9, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.1843,  0.0000, -3.6405, -3.7502,
     & -2.7458, -2.5091, -2.9053, -2.9983, -2.8186, -3.6835,
     & -3.8834, -3.7689, -3.6539, -3.7782, -3.5595, -3.6045,
     & -3.6191, -3.5817, -3.3414, -3.4745, -3.3315, -3.3442,
     & -3.4382, -3.2088, -3.1964, -3.1861, -3.1443, -3.0745,
     & -3.1722, -2.9643, -3.1207, -3.0917, -3.0297, -3.0379,
     & -2.9755, -2.9075, -2.8465, -2.8674, -2.8535, -2.8924,
     & -2.7941, -2.8394, -2.7620, -2.8676, -2.7952, -2.8107,
     & -2.8259, -2.8333, -2.8725, -2.8778, -2.8507, -2.8598/
c PAIR NR C3
      data (pmf( 9, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.8573,  0.0000,  0.0000, -3.2786, -2.7623,
     &  0.0000, -2.5124, -2.3999, -2.2944, -2.6038, -3.2851,
     & -3.2656, -3.9111, -3.8560, -3.7048, -3.5796, -3.5578,
     & -3.5718, -3.6380, -3.3953, -3.3645, -3.4502, -3.5102,
     & -3.4915, -3.3291, -3.1892, -3.1540, -2.9629, -3.1704,
     & -3.0387, -2.9683, -2.9142, -2.9822, -2.9491, -2.7975,
     & -3.0525, -2.9210, -2.8662, -2.8643, -2.8039, -2.8705,
     & -2.8610, -2.7100, -2.7975, -2.9318, -2.9041, -3.0804,
     & -2.9415, -2.9907, -2.7954, -2.8012, -2.7361, -2.8950/
c PAIR NR CW
      data (pmf( 9, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -3.2500, -3.3712, -2.6157,
     & -2.7188, -2.5924, -2.9745, -3.3518, -3.3175, -3.6966,
     & -3.9756, -3.7158, -3.5038, -3.3799, -3.5558, -3.5156,
     & -3.4383, -3.5221, -3.4973, -3.4271, -3.5153, -3.2817,
     & -3.3917, -3.2790, -3.2374, -3.2487, -3.1719, -3.1568,
     & -3.3104, -3.1261, -3.1309, -3.1577, -3.0229, -3.0742,
     & -3.0921, -3.1163, -2.8679, -2.9898, -2.9148, -2.9814,
     & -2.9381, -2.8926, -2.9265, -2.6692, -2.6803, -2.7931,
     & -2.7246, -2.8230, -2.7866, -2.7628, -2.7425, -2.5454/
c PAIR NR CO
      data (pmf( 9, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -3.8233, -3.7160,
     & -2.7334, -2.6240, -2.8231, -3.1921, -3.7127, -3.9042,
     & -3.9730, -3.6335, -3.5693, -3.7474, -3.5449, -3.6131,
     & -3.6861, -3.5066, -3.6545, -3.5502, -3.5156, -3.2147,
     & -3.2320, -3.0975, -2.8845, -3.0773, -3.0981, -2.8862,
     & -3.0174, -3.0286, -2.9310, -3.0184, -3.0862, -3.1252,
     & -3.0738, -2.9212, -2.7209, -2.8178, -2.9910, -2.7273,
     & -3.0214, -3.0218, -3.0256, -2.9209, -2.7647, -2.6164,
     & -2.6610, -2.5618, -2.6142, -2.8809, -2.8795, -2.8328/
c PAIR NR CN
      data (pmf( 9, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.7059, -2.5959, -3.1407, -3.2130, -3.7398,
     & -3.6736, -3.6340, -3.5736, -3.2217, -3.3230, -3.3438,
     & -3.4427, -3.4513, -3.3562, -3.2767, -3.1991, -3.2820,
     & -3.2988, -3.1272, -3.1391, -3.2863, -3.0721, -3.1185,
     & -3.0823, -3.0003, -3.1760, -3.3767, -2.8503, -2.9373,
     & -2.9675, -2.9748, -3.0685, -3.0101, -2.9902, -2.8809,
     & -2.9662, -2.7444, -2.6887, -2.7450, -2.7572, -2.7523,
     & -3.0162, -3.1154, -2.9274, -2.9010, -2.7673, -2.7333/
c PAIR NR NC
      data (pmf( 9, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.4679,
     & -2.3553, -2.2507, -3.3009, -2.8788, -3.5719, -3.5913,
     & -3.7245, -3.8030, -3.8078, -3.3107, -3.0164, -2.9925,
     & -3.2696, -3.4382, -3.2340, -3.3180, -3.3624, -3.1633,
     & -3.3333, -3.1693, -2.9939, -2.9513, -3.2005, -2.7322,
     & -3.1402, -2.8602, -3.0122, -3.1205, -2.9826, -3.0662,
     & -2.9536, -2.9754, -3.0885, -3.1413, -2.8712, -2.9125,
     & -2.9041, -3.0486, -2.9215, -2.9358, -2.8142, -2.9111,
     & -2.9761, -2.9184, -2.9654, -2.6937, -2.8497, -2.7956/
c PAIR NR NP
      data (pmf( 9,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -3.8580,  0.0000,
     &  0.0000, -3.1675, -4.0246,  0.0000, -2.7554, -2.8586,
     & -3.4435, -3.0257, -2.9170, -2.8149, -2.7162, -3.4158,
     & -3.3999, -3.2928, -2.9496, -3.2224, -3.2058, -3.0081,
     & -3.4335, -3.4904, -3.1376, -3.3097, -3.2178, -3.2726,
     & -3.1863, -3.2364, -3.3069, -3.1495, -3.2561, -3.1631,
     & -2.9428, -3.2046, -2.9856, -2.9669, -2.8056, -2.6779,
     & -2.7039, -2.6912, -2.6665, -3.0325, -2.9329, -2.8340,
     & -3.0458, -2.8994, -2.7908, -2.8439, -2.8954, -2.7824/
c PAIR NR NA
      data (pmf( 9,11,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR NR ND
      data (pmf( 9,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -3.8405, -2.9153,
     & -3.5887, -3.9958, -4.1847, -3.8631, -3.5150, -3.4711,
     & -3.6253, -2.6161, -3.1693, -3.1977, -3.5935, -3.8153,
     & -3.4956, -3.3278, -3.1163, -3.3320, -3.7831, -3.1716,
     & -3.1423, -3.2019, -2.9623, -2.8688, -2.7502, -3.1478,
     & -3.3123, -3.0155, -3.5696, -2.9724, -2.8948, -3.0959,
     & -2.9265, -3.2284, -3.3955, -3.1090, -2.7405, -2.9071,
     & -2.9198, -3.0540, -2.8590, -2.5646, -2.6381, -2.6415,
     & -2.9964, -2.7489, -2.7653, -2.6567, -2.8165, -2.6668/
c PAIR NR NR
      data (pmf( 9,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -3.0877, -2.7412,
     & -2.8508, -3.7874, -4.2657, -4.1227, -3.4696, -3.4299,
     & -3.3337, -3.2122, -3.5791, -3.6366, -3.6588, -3.7874,
     & -3.4023, -3.3192, -3.1216, -3.2580, -3.1391, -3.2949,
     & -3.2472, -3.2242, -3.1859, -3.1029, -3.0795, -2.9493,
     & -3.0635, -3.0196, -3.2205, -3.0175, -3.0218, -2.9718,
     & -2.9762, -3.0932, -3.0080, -2.8893, -2.9772, -2.6692,
     & -2.7542, -2.7227, -2.7761, -2.8207, -3.0254, -2.8996,
     & -2.9500, -2.8690, -2.8673, -2.9698, -2.8631, -2.8024/
c PAIR NR N0
      data (pmf( 9,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR NR NS
      data (pmf( 9,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -3.9897, -4.2185, -3.8720,
     &  0.0000, -3.2291, -3.1236, -4.0121, -4.2050, -4.0649,
     & -3.9516, -3.8415, -4.0840, -4.0378, -3.5537, -3.7009,
     & -3.7079, -3.7835, -4.0291, -3.4771, -3.5954, -3.2884,
     & -3.2856, -3.3530, -3.1572, -2.9948, -3.4415, -2.9468,
     & -2.8286, -3.0521, -2.8669, -3.0741, -3.2367, -2.9787,
     & -3.2105, -2.9014, -2.9722, -2.8282, -2.6318, -2.7373,
     & -2.6495, -2.6024, -2.5058, -2.0981, -2.1635, -2.0975,
     & -2.7388, -2.7165, -2.6651, -2.7572, -2.7619, -2.4446/
c PAIR NR OC
      data (pmf( 9,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -3.3754, -3.5695,
     & -3.4286, -3.9420, -3.8348, -3.7770, -3.7126, -3.5522,
     & -3.7074, -3.5296, -3.6577, -3.6009, -3.6737, -3.6741,
     & -3.5436, -3.4821, -3.2810, -3.3941, -3.3329, -3.2256,
     & -3.2069, -3.1885, -2.9690, -3.0016, -2.9212, -2.9323,
     & -2.9407, -3.0387, -3.0123, -3.0086, -2.9688, -2.8596,
     & -2.9691, -2.9053, -2.9992, -3.1450, -3.0840, -3.0316,
     & -2.9592, -2.8782, -2.8283, -2.9416, -2.7867, -2.6970,
     & -2.8258, -2.7886, -2.7703, -2.7059, -2.8643, -2.7706/
c PAIR NR OA
      data (pmf( 9,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -2.9595, -3.4491,
     & -3.3372, -3.7407, -3.7472, -3.9723, -3.9186, -3.7837,
     & -3.6248, -3.6032, -3.7976, -3.7645, -3.5093, -3.6523,
     & -3.3318, -3.4646, -3.3575, -3.3860, -3.0897, -3.0678,
     & -3.1126, -3.0412, -3.0413, -2.9954, -3.2717, -3.1275,
     & -3.2278, -3.0648, -3.0445, -3.3083, -3.0546, -3.0227,
     & -2.9637, -2.8364, -2.9659, -3.1298, -3.0642, -2.8942,
     & -2.8153, -2.9365, -2.9008, -2.8833, -2.9464, -2.7763,
     & -2.6477, -2.8577, -2.6190, -2.7460, -2.7327, -2.7197/
c PAIR NR OE
      data (pmf( 9,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.9391, -3.4693, -3.1228,
     & -2.5631, -2.4243, -2.7045, -2.8237, -3.2900, -3.4125,
     & -3.0640, -3.2476, -3.0519, -3.0585, -3.2206, -3.2408,
     & -3.4442, -3.2839, -3.5719, -3.6523, -3.5802, -3.4912,
     & -3.2315, -3.3010, -3.3510, -3.3404, -3.2156, -3.2941,
     & -3.3066, -3.1608, -3.2310, -3.2095, -3.0189, -3.0675,
     & -2.9630, -3.0574, -2.9339, -3.0260, -2.9008, -2.9411,
     & -2.8478, -2.7787, -3.0211, -2.9924, -2.9117, -2.8145,
     & -2.5347, -2.7035, -2.7714, -2.6950, -2.6616, -2.7557/
c PAIR NR OR
      data (pmf( 9,19,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR NR OS
      data (pmf( 9,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.8045, -2.6230, -3.1532, -3.1854,
     & -3.2361, -3.7649, -3.5325, -3.6315, -3.6535, -3.4454,
     & -3.3251, -3.3018, -3.3230, -3.6613, -3.4304, -3.6909,
     & -3.7065, -3.4645, -3.3935, -3.3372, -3.1244, -3.2501,
     & -3.2627, -3.0963, -3.0143, -3.1324, -3.1469, -3.1028,
     & -3.1608, -3.0377, -3.0319, -3.0120, -2.7929, -2.9202,
     & -3.0612, -2.8338, -2.9181, -2.8961, -2.9609, -2.7983,
     & -2.9302, -3.0260, -2.8992, -2.7773, -2.8400, -2.9396,
     & -2.8972, -2.9034, -2.8465, -2.9082, -2.8524, -2.7681/
c PAIR NR OD
      data (pmf( 9,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.1421,  0.0000,  0.0000, -2.3431, -3.2127, -3.1412,
     & -3.5432, -3.9944, -3.9554, -3.7848, -3.6408, -3.5847,
     & -3.4357, -3.2224, -3.2600, -3.5301, -3.7855, -3.7870,
     & -3.7068, -3.5070, -3.2859, -3.3422, -3.2571, -3.3073,
     & -3.3408, -3.3138, -3.1743, -3.1162, -3.0442, -3.1706,
     & -3.1896, -3.1949, -3.2094, -3.0271, -3.1639, -3.1123,
     & -3.1319, -3.1457, -3.0240, -3.0121, -2.9528, -2.8120,
     & -2.8502, -2.7445, -2.7879, -2.7853, -2.7587, -2.8869,
     & -2.7484, -2.7297, -2.6914, -2.6964, -2.6371, -2.6019/
c PAIR NR P
      data (pmf( 9,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -3.2849, -3.4683, -3.5596,
     &  0.0000, -2.8206, -2.8826, -2.3724, -3.0948, -3.7606,
     & -4.1004, -3.6975, -3.6651, -3.4484, -3.1630, -3.2816,
     & -3.3873, -3.0093, -3.4808, -3.4448, -3.6001, -3.5312,
     & -3.3332, -3.0301, -2.9266, -3.1428, -2.8737, -3.0227,
     & -3.0555, -2.8669, -2.7372, -2.8422, -2.8922, -2.9148,
     & -3.0647, -3.1665, -3.1007, -2.8998, -2.9818, -2.8826,
     & -2.8772, -2.9040, -2.6906, -2.8099, -2.9386, -2.8449,
     & -2.9198, -2.9641, -3.0869, -2.9134, -2.9020, -2.8775/
c PAIR NR SA
      data (pmf( 9,23,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR NR SD
      data (pmf( 9,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR NR HL
      data (pmf( 9,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -3.0444, -3.3461,
     & -3.1245, -3.3653, -3.1966, -3.3899, -3.1472, -3.2568,
     & -3.2870, -3.1708, -3.3574, -3.6040, -3.2114, -3.2686,
     & -3.3060, -3.3437, -3.2103, -3.2730, -3.1076, -3.2160,
     & -2.9203, -3.1042, -3.0600, -3.3504, -3.0457, -3.2548,
     & -3.2046, -3.1463, -3.0881, -3.0993, -3.1589, -3.1088,
     & -2.6931, -3.0654, -3.0291, -3.1410, -2.9519, -2.9481,
     & -2.9612, -2.9737, -2.8375, -2.9072, -2.9573, -2.9005,
     & -2.9725, -2.8382, -2.8623, -2.9100, -2.7864, -2.8478/
c PAIR NR F
      data (pmf( 9,26,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OC CF
      data (pmf(10, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.9178, -2.9074, -3.2091,
     & -2.8016, -2.0172, -2.8122, -2.8482, -3.2494, -3.4831,
     & -3.3562, -3.2037, -3.2067, -3.1822, -3.1909, -3.3044,
     & -3.2800, -3.1496, -3.1851, -3.1723, -3.1367, -3.0681,
     & -3.0414, -3.1317, -3.1356, -3.1221, -3.1193, -3.0924,
     & -3.0944, -3.1407, -3.0884, -3.0409, -3.0620, -3.0968,
     & -3.0836, -3.0525, -3.0228, -2.9953, -2.9939, -2.9830,
     & -3.0077, -2.9578, -2.8969, -2.9259, -2.9105, -2.8501,
     & -2.8862, -2.9042, -2.9159, -2.9371, -2.9637, -2.9156/
c PAIR OC CP
      data (pmf(10, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.2271, -2.0292, -3.1100, -3.7828, -3.7671,
     & -3.5664, -2.9672, -2.8146, -3.3899, -3.8713, -4.0698,
     & -3.9760, -3.8433, -3.7161, -3.6966, -3.5670, -3.5610,
     & -3.5306, -3.4637, -3.4516, -3.3778, -3.3041, -3.2184,
     & -3.2545, -3.1530, -3.1356, -3.1209, -3.0487, -3.0942,
     & -3.0146, -3.0101, -2.9764, -2.9455, -2.9346, -2.9679,
     & -2.9546, -2.9507, -2.9601, -2.9679, -2.9321, -2.9106,
     & -2.8920, -2.8767, -2.9332, -2.9042, -2.8782, -2.8238,
     & -2.8372, -2.8414, -2.8354, -2.8312, -2.8217, -2.8187/
c PAIR OC cF
      data (pmf(10, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.3010, -3.0297, -3.3978,
     & -3.4562, -2.6836, -2.7243, -2.7591, -3.0981, -3.0174,
     & -3.1169, -2.8919, -2.9910, -3.1597, -3.1187, -3.1271,
     & -3.1255, -3.1440, -2.9814, -3.0297, -3.1121, -3.0806,
     & -3.1278, -3.1028, -3.0891, -3.1176, -3.1675, -3.1348,
     & -3.0898, -3.1235, -3.1505, -3.1508, -3.1499, -3.1053,
     & -3.1690, -3.0988, -3.1030, -3.0905, -3.0587, -3.0091,
     & -3.0171, -3.0160, -2.9494, -2.9701, -2.9329, -2.9327,
     & -2.8798, -2.9408, -2.8935, -2.7959, -2.8114, -2.7741/
c PAIR OC cP
      data (pmf(10, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.6333, -2.3201, -3.0558, -3.1269,
     & -3.1636, -2.0656, -2.3543, -2.8733, -3.2625, -3.6578,
     & -3.6209, -3.4508, -3.3339, -3.2825, -3.3590, -3.4676,
     & -3.5072, -3.4601, -3.3146, -3.3517, -3.3494, -3.3289,
     & -3.1741, -3.0373, -2.9527, -2.9560, -2.9971, -3.0769,
     & -3.0953, -3.0959, -3.0110, -3.0122, -3.0313, -2.9370,
     & -2.9967, -3.0062, -3.0203, -2.9976, -2.9775, -2.9487,
     & -2.9142, -2.9488, -2.9262, -2.8513, -2.9154, -2.9594,
     & -2.9402, -2.9402, -2.9498, -2.8714, -2.8478, -2.9252/
c PAIR OC C3
      data (pmf(10, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.7125, -3.0034, -3.7141,
     & -3.0968,  0.0000, -2.3641, -2.6674, -2.9769, -3.2085,
     & -3.3995, -3.1748, -3.0209, -3.4188, -3.2688, -3.3972,
     & -3.2783, -3.2258, -3.0372, -3.1817, -3.1388, -3.0673,
     & -2.8999, -2.8728, -2.9747, -3.0486, -2.9726, -2.9216,
     & -3.0310, -2.8952, -2.9089, -3.0712, -3.0829, -2.9530,
     & -2.9376, -2.8894, -2.9250, -3.0790, -3.0140, -3.0752,
     & -3.0115, -2.9906, -3.0503, -3.0411, -2.8981, -3.0265,
     & -2.9895, -2.9548, -3.0509, -3.0059, -2.9522, -2.9656/
c PAIR OC CW
      data (pmf(10, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.2103, -2.4718, -3.5675, -3.3708, -3.3780,
     & -3.3665, -3.1155, -3.1760, -3.3133, -3.6432, -3.8540,
     & -3.8842, -3.6441, -3.6179, -3.5699, -3.5649, -3.5118,
     & -3.4895, -3.5199, -3.4304, -3.3468, -3.2793, -3.3415,
     & -3.2563, -3.1299, -3.1102, -3.1529, -3.1200, -3.2294,
     & -3.1194, -3.0516, -3.1247, -3.0156, -3.0057, -3.0153,
     & -2.9620, -2.9831, -2.9055, -2.9813, -2.9256, -2.9335,
     & -2.8598, -2.9218, -2.9108, -2.8769, -2.8424, -2.8754,
     & -2.8637, -2.8418, -2.8689, -2.7448, -2.7809, -2.8545/
c PAIR OC CO
      data (pmf(10, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.3957, -3.0956, -3.3592,
     & -3.5431, -3.0916, -3.1998, -3.2807, -3.6097, -3.5993,
     & -3.6162, -3.3928, -3.3985, -3.3971, -3.4095, -3.3327,
     & -3.4333, -3.3834, -3.3159, -3.2238, -3.0804, -3.1285,
     & -3.0733, -3.0701, -3.1087, -3.0362, -2.9470, -3.0369,
     & -3.0428, -3.0830, -3.1008, -2.9654, -3.0124, -2.9720,
     & -3.0415, -2.9706, -3.0113, -2.9143, -2.8581, -3.0199,
     & -3.0020, -2.9811, -2.9523, -2.9217, -2.8595, -2.9123,
     & -2.9886, -2.9713, -3.0124, -2.8860, -2.8126, -2.7684/
c PAIR OC CN
      data (pmf(10, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -3.0056, -2.5862, -2.6487,
     & -3.2613, -2.2347, -2.8913, -3.0090, -3.4948, -3.8905,
     & -4.1617, -3.8954, -3.5865, -3.4025, -3.1880, -3.2282,
     & -3.2033, -3.2489, -3.2130, -3.3919, -3.1679, -3.2283,
     & -3.1539, -3.1976, -3.1659, -3.2014, -3.2046, -3.3643,
     & -3.1405, -2.9807, -2.9191, -2.9093, -2.9344, -2.7907,
     & -3.0139, -2.8653, -2.9812, -2.8993, -2.8632, -2.8529,
     & -2.9461, -2.9145, -2.8856, -2.8525, -2.7615, -2.9965,
     & -2.8832, -2.9264, -2.8727, -2.8847, -2.9736, -2.8786/
c PAIR OC NC
      data (pmf(10, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.7823, -2.9526, -2.2940, -3.0741,
     & -3.3705, -4.1748, -4.3254, -4.0819, -3.6898, -3.8030,
     & -3.7402, -3.6323, -3.2760, -3.3088, -3.3444, -3.1991,
     & -3.2432, -3.1409, -3.1853, -3.1116, -3.0879, -3.0475,
     & -2.9470, -3.0487, -3.0787, -3.0928, -3.1624, -3.1070,
     & -3.1447, -3.1197, -3.1019, -3.0945, -3.0431, -3.1043,
     & -3.0026, -3.0144, -2.9405, -2.9095, -2.9449, -2.9207,
     & -2.8348, -2.8309, -2.8121, -2.8472, -2.8843, -2.9154,
     & -2.9259, -2.8543, -2.8811, -2.8033, -2.8525, -2.9128/
c PAIR OC NP
      data (pmf(10,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.9335, -2.7860,  0.0000, -3.1851, -3.3792,
     & -2.8970, -2.7339, -3.5728, -2.8574, -3.1869, -3.2092,
     & -3.4569, -3.6492, -3.6765, -3.7198, -3.9101, -3.6499,
     & -3.2004, -3.1927, -3.2336, -3.3623, -3.2396, -3.0642,
     & -3.2391, -2.9716, -3.0452, -3.0457, -3.2325, -3.2383,
     & -3.1333, -3.0360, -3.1167, -2.9385, -3.0444, -3.0324,
     & -2.9645, -2.9926, -2.9288, -2.8272, -2.9024, -2.9265,
     & -2.9906, -2.9165, -2.9882, -3.0266, -2.9107, -2.8386,
     & -2.8346, -2.8326, -2.8537, -2.9785, -2.9402, -2.8445/
c PAIR OC NA
      data (pmf(10,11,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OC ND
      data (pmf(10,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.6375, -2.9285, -3.2302,
     & -3.5641, -3.7700, -3.9091, -3.6902, -3.6088, -3.5940,
     & -3.2705, -3.1007, -3.1374, -3.2021, -3.3823, -3.5127,
     & -3.4894, -3.3554, -3.4312, -3.2748, -3.2641, -3.2213,
     & -3.3784, -3.3075, -3.1370, -3.1390, -2.9490, -3.0332,
     & -2.9597, -2.9065, -3.2350, -3.0533, -3.1542, -2.9916,
     & -3.0885, -3.0267, -3.0476, -2.9816, -2.8459, -2.8481,
     & -2.9477, -2.9412, -2.8765, -3.1030, -2.9296, -2.6715,
     & -2.7648, -2.6596, -2.9074, -3.0751, -2.8960, -2.9015/
c PAIR OC NR
      data (pmf(10,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.4490, -2.3179, -2.4391, -3.0422,
     & -3.2593, -4.0547, -3.9063, -3.5465, -3.5786, -3.6622,
     & -3.5725, -3.4680, -3.4264, -3.5332, -3.4048, -3.3870,
     & -3.5366, -3.3303, -3.0606, -3.1418, -3.2281, -3.2453,
     & -3.2228, -3.2628, -3.1679, -3.3592, -3.2129, -3.1851,
     & -3.1620, -3.1232, -3.0354, -2.9999, -2.9245, -2.9377,
     & -2.8706, -3.0193, -3.0319, -2.8792, -2.8972, -2.8899,
     & -2.8705, -2.8840, -2.8941, -2.8400, -2.8694, -2.8427,
     & -2.9302, -2.9005, -2.8732, -2.8462, -2.9108, -2.8841/
c PAIR OC N0
      data (pmf(10,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OC NS
      data (pmf(10,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -2.6969, -3.8165,
     & -4.1737, -3.9534, -3.7162, -3.7868, -3.9074, -3.4307,
     & -3.4334, -3.4284, -3.7573, -3.3077, -3.3572, -3.6031,
     & -3.3419, -2.9888, -3.2755, -3.3234, -3.2953, -3.1253,
     & -3.0182, -2.9523, -3.0134, -3.0263, -3.0846, -2.9858,
     & -3.2213, -3.2636, -3.1415, -3.0770, -2.8648, -2.9055,
     & -3.0281, -3.1198, -3.1853, -3.0210, -3.1169, -2.7642,
     & -2.8650, -2.6010, -2.6894, -2.7602, -2.9591, -2.8739,
     & -2.8265, -2.7401, -2.8003, -2.8124, -2.9332, -2.9613/
c PAIR OC OC
      data (pmf(10,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.7426, -2.0205, -2.7204, -3.6578,
     & -3.7822, -3.6959, -3.4301, -3.4020, -3.3266, -3.4448,
     & -3.2260, -3.3083, -3.1942, -3.1851, -3.3423, -3.2259,
     & -3.2665, -3.2323, -3.3321, -3.2348, -3.1334, -3.0803,
     & -2.9786, -3.1383, -3.0250, -3.0816, -3.0802, -3.0407,
     & -3.1713, -3.1319, -3.1472, -3.0555, -2.9780, -3.0196,
     & -2.9714, -2.9645, -2.9747, -2.9983, -3.0517, -3.0049,
     & -2.9772, -2.9329, -2.9900, -2.9258, -2.8928, -2.8572,
     & -2.9332, -2.9037, -2.8532, -2.8845, -2.8718, -2.9000/
c PAIR OC OA
      data (pmf(10,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.7717, -1.6537, -2.9050,
     & -3.3337, -3.5070, -3.6683, -3.5063, -3.6031, -3.6875,
     & -3.6009, -3.4372, -3.3519, -3.2901, -3.4441, -3.3012,
     & -3.5019, -3.4521, -3.3856, -3.0759, -3.1045, -3.2463,
     & -3.1670, -3.1321, -3.0026, -3.1003, -3.2478, -3.1813,
     & -3.0983, -2.9727, -2.9806, -3.0698, -3.0570, -2.9721,
     & -2.9993, -2.9986, -2.9857, -2.8912, -2.8790, -2.9328,
     & -2.8433, -2.8214, -2.8989, -2.8758, -2.8986, -3.0348,
     & -2.8799, -2.9129, -2.9188, -2.9052, -2.9114, -2.8981/
c PAIR OC OE
      data (pmf(10,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.4946, -2.3471, -2.8641, -3.2460, -3.7869,
     & -3.5772, -3.1672, -3.2439, -3.1577, -3.3059, -3.3997,
     & -3.4675, -3.5740, -3.2102, -3.3534, -3.4168, -3.4623,
     & -3.4666, -3.3699, -3.5062, -3.5006, -3.4394, -3.3163,
     & -3.3180, -3.2909, -3.1531, -3.1332, -3.1212, -3.1404,
     & -3.0253, -2.9354, -2.9890, -2.9974, -3.0013, -3.0106,
     & -2.9400, -2.9390, -2.9827, -2.9530, -2.9338, -2.9264,
     & -2.9318, -2.8811, -2.9075, -2.8837, -3.0051, -2.8679,
     & -2.9771, -2.9061, -2.8312, -2.8571, -2.8867, -2.8491/
c PAIR OC OR
      data (pmf(10,19,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OC OS
      data (pmf(10,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.2992, -2.2088, -2.5776, -2.7525, -3.2848,
     & -3.6403, -3.2223, -3.0840, -3.2129, -3.0755, -3.0565,
     & -3.0731, -3.0831, -3.0977, -3.2497, -3.2745, -3.2043,
     & -3.1525, -3.1692, -3.2146, -3.1478, -3.1763, -3.1505,
     & -3.0577, -3.1605, -3.1785, -3.1843, -3.1299, -3.1479,
     & -3.1130, -3.1308, -3.1324, -3.1575, -3.1227, -3.0969,
     & -3.0638, -3.0504, -3.0436, -3.0012, -3.0291, -3.0392,
     & -3.0201, -2.9436, -2.9344, -2.9137, -2.9152, -2.9006,
     & -2.9198, -2.8721, -2.8613, -2.8316, -2.8287, -2.7679/
c PAIR OC OD
      data (pmf(10,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.1364, -2.8232, -3.6779, -3.6932,
     & -4.1891, -4.4604, -4.2747, -4.0662, -3.9020, -3.7906,
     & -3.6552, -3.4910, -3.5278, -3.5367, -3.5093, -3.5455,
     & -3.5160, -3.5015, -3.3982, -3.2141, -3.2129, -3.2154,
     & -3.1731, -3.1764, -3.1906, -3.1146, -3.1415, -3.1283,
     & -3.0251, -2.9608, -2.9651, -2.9596, -2.9450, -2.9216,
     & -2.9464, -2.9266, -2.9628, -2.9818, -2.8886, -2.8656,
     & -2.8504, -2.8036, -2.8170, -2.8228, -2.8859, -2.8422,
     & -2.8352, -2.8539, -2.8493, -2.8137, -2.8184, -2.8246/
c PAIR OC P
      data (pmf(10,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.9729, -3.0030, -2.6972,
     & -2.9202, -2.6567, -2.4580, -2.4477, -2.9094, -2.8473,
     & -3.2797, -3.2672, -3.0166, -3.2497, -3.0883, -3.0939,
     & -3.3134, -3.3809, -3.2755, -3.1596, -3.1999, -3.0157,
     & -3.1896, -3.1877, -3.2366, -3.2843, -3.2645, -3.1293,
     & -3.2352, -3.1224, -3.0699, -3.1095, -3.1309, -3.1443,
     & -3.0662, -3.1421, -3.1263, -3.0369, -3.0270, -3.0436,
     & -2.9650, -2.9854, -3.0426, -3.0379, -2.9989, -2.8540,
     & -2.8805, -2.8650, -2.6925, -2.7398, -2.7405, -2.7913/
c PAIR OC SA
      data (pmf(10,23,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OC SD
      data (pmf(10,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OC HL
      data (pmf(10,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.5311, -3.0317, -3.5488, -3.7579, -3.6871,
     & -3.7095, -3.8006, -3.6227, -3.6837, -3.4641, -3.5738,
     & -3.6376, -3.4398, -3.3053, -3.4112, -3.3179, -3.3765,
     & -3.3229, -3.3388, -3.2823, -3.2139, -3.3133, -3.2624,
     & -3.2084, -3.2395, -3.2845, -3.2551, -3.2529, -3.1981,
     & -3.1108, -3.0831, -3.1216, -3.0931, -3.0638, -2.9398,
     & -2.9582, -2.9159, -2.9278, -3.0417, -2.9915, -2.9303,
     & -2.9858, -2.9026, -2.8948, -2.8834, -2.8580, -2.8979,
     & -2.9038, -2.7986, -2.8311, -2.7663, -2.7934, -2.7701/
c PAIR OC F
      data (pmf(10,26,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OA CF
      data (pmf(11, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -0.3327,  0.0000, -0.7317, -2.2352,
     & -2.1183, -1.2639, -1.5934, -1.9787, -2.5556, -2.8043,
     & -2.9182, -2.8313, -2.7506, -2.6910, -2.7552, -2.7232,
     & -2.6404, -2.7249, -2.7574, -2.7644, -2.8041, -2.8309,
     & -2.8252, -2.8797, -2.8797, -2.9162, -2.9575, -3.0225,
     & -3.0407, -3.0712, -3.0673, -3.0581, -3.0566, -3.0641,
     & -3.0572, -3.0616, -3.0489, -3.0561, -3.0485, -3.0393,
     & -3.0443, -3.0164, -3.0427, -3.0431, -3.0479, -3.0468,
     & -3.0283, -3.0383, -3.0451, -3.0300, -3.0161, -3.0118/
c PAIR OA CP
      data (pmf(11, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.1208, -1.3035, -1.2512, -1.9130, -2.5449,
     & -2.4674, -1.7823, -1.9814, -2.4515, -2.8087, -2.9140,
     & -2.8256, -2.7850, -2.7961, -2.7963, -2.7792, -2.7338,
     & -2.7452, -2.7756, -2.8237, -2.7815, -2.8498, -2.8552,
     & -2.9348, -2.9119, -2.9118, -2.9333, -2.9704, -2.9927,
     & -3.0117, -3.0323, -3.0145, -3.0401, -3.0225, -3.0329,
     & -3.0326, -3.0079, -3.0293, -3.0422, -3.0779, -3.0754,
     & -3.0769, -3.0606, -3.0653, -3.0558, -3.0386, -3.0353,
     & -3.0209, -3.0102, -3.0065, -3.0090, -3.0191, -3.0117/
c PAIR OA cF
      data (pmf(11, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -1.7284, -2.6446,
     & -2.2422, -1.5521, -1.6836, -2.3382, -2.5915, -2.5891,
     & -2.7108, -2.8198, -2.7870, -2.7374, -2.7352, -2.7450,
     & -2.7310, -2.8174, -2.7850, -2.7975, -2.8163, -2.8534,
     & -2.8592, -2.8601, -2.8812, -2.9057, -2.9607, -2.9585,
     & -2.9884, -3.0233, -3.0463, -3.0348, -3.0459, -3.0604,
     & -3.0594, -3.0599, -3.0467, -3.0759, -3.0577, -3.0665,
     & -3.0625, -3.0404, -3.0454, -3.0228, -3.0325, -3.0412,
     & -3.0377, -3.0396, -3.0299, -3.0224, -3.0173, -3.0260/
c PAIR OA cP
      data (pmf(11, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -0.3961, -1.6928, -2.0459,
     & -1.7990, -1.4594, -2.3054, -2.7168, -2.7696, -2.7634,
     & -2.8038, -2.6762, -2.7324, -2.7365, -2.7007, -2.7187,
     & -2.8514, -2.9174, -2.9373, -2.9461, -2.9905, -2.9858,
     & -2.9559, -2.9257, -2.9411, -2.9529, -3.0153, -3.0498,
     & -3.0951, -3.0827, -3.0767, -3.0720, -3.0692, -3.0696,
     & -3.0246, -3.0409, -3.0259, -3.0245, -3.0533, -3.0527,
     & -3.0492, -3.0384, -3.0380, -3.0461, -2.9991, -3.0041,
     & -3.0100, -2.9966, -2.9915, -2.9861, -2.9863, -2.9846/
c PAIR OA C3
      data (pmf(11, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -1.3647, -2.5201,
     & -2.0153, -0.5986, -1.5432, -1.8938, -2.4115, -2.6543,
     & -2.6520, -2.7519, -2.5160, -2.7340, -2.9013, -2.8812,
     & -2.6872, -2.7813, -2.6760, -2.7460, -2.8006, -2.8934,
     & -2.8678, -2.8343, -2.7813, -2.8660, -2.8901, -2.9537,
     & -2.9792, -3.0168, -3.0478, -3.0831, -3.0673, -3.0331,
     & -3.0526, -3.0725, -3.0415, -3.0983, -3.0682, -3.0652,
     & -3.0667, -3.0614, -3.0453, -3.0591, -3.0518, -3.0239,
     & -3.0462, -3.0129, -3.0287, -3.0371, -3.0492, -3.0247/
c PAIR OA CW
      data (pmf(11, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.3072, -0.7802, -2.2300,
     & -1.7636, -1.5585, -1.9284, -2.1689, -2.6783, -3.0257,
     & -2.9014, -3.0073, -2.9645, -2.7850, -2.7938, -2.8017,
     & -2.7605, -2.7273, -2.7514, -2.8037, -2.7515, -2.8968,
     & -2.8961, -2.8602, -2.8827, -2.8738, -2.9502, -2.9680,
     & -2.9682, -3.0016, -2.9997, -3.0670, -3.0013, -3.0511,
     & -3.0030, -3.0183, -3.0226, -3.0338, -3.0487, -3.0477,
     & -3.0574, -3.0550, -3.0540, -3.0320, -3.0620, -3.0658,
     & -3.0521, -3.0781, -3.0517, -3.0395, -3.0256, -2.9917/
c PAIR OA CO
      data (pmf(11, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -1.4847, -2.2425,
     & -2.0778, -1.7292, -1.3496, -1.8609, -2.4721, -2.4424,
     & -2.3791, -2.5897, -2.5510, -2.5365, -2.5336, -2.6532,
     & -2.7462, -2.8096, -2.8709, -2.8703, -2.8809, -2.8633,
     & -2.8993, -2.8678, -2.9284, -2.8843, -2.9410, -2.9415,
     & -2.9674, -2.9859, -2.9953, -3.0172, -3.0574, -3.0529,
     & -3.0427, -3.0692, -3.0696, -3.0550, -3.0408, -3.0487,
     & -3.1166, -3.0879, -3.0666, -3.0327, -3.0314, -3.0236,
     & -3.0143, -3.0177, -3.0486, -3.0533, -3.0486, -3.0302/
c PAIR OA CN
      data (pmf(11, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.1351,  0.0000, -1.5342, -1.9268,
     & -1.7093, -1.5916, -1.7830, -2.4976, -2.7865, -3.0082,
     & -3.2445, -3.3562, -3.1099, -2.8515, -2.6468, -2.7944,
     & -2.7455, -2.7943, -2.9050, -2.9699, -2.9993, -3.0118,
     & -3.0233, -2.9606, -3.0278, -3.0877, -3.0338, -3.0669,
     & -3.0145, -3.0665, -3.0483, -3.0537, -3.1339, -3.0930,
     & -3.0765, -3.0454, -3.0239, -3.0053, -3.0309, -3.0066,
     & -2.9948, -3.0017, -2.9895, -2.9833, -3.0154, -2.9876,
     & -2.9929, -2.9766, -2.9871, -2.9685, -3.0079, -3.0198/
c PAIR OA NC
      data (pmf(11, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.6005,  0.0000, -2.1418,
     & -2.6950, -3.4379, -3.6973, -3.4401, -2.9385, -2.8066,
     & -2.6165, -2.6807, -2.6342, -2.5194, -2.6903, -2.8323,
     & -2.9036, -2.9183, -3.0081, -2.9787, -3.0602, -2.9476,
     & -2.9188, -2.9663, -2.9655, -2.9997, -2.9919, -3.0405,
     & -3.0716, -3.0817, -3.1056, -3.0988, -3.0800, -3.1125,
     & -3.0584, -3.0606, -3.0352, -3.0206, -3.0042, -2.9823,
     & -2.9890, -2.9860, -2.9369, -3.0224, -3.0241, -3.0092,
     & -2.9937, -2.9859, -2.9803, -3.0090, -3.0051, -2.9838/
c PAIR OA NP
      data (pmf(11,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.7491, -2.0270, -2.1482, -2.0409,
     & -1.2119, -1.6970, -2.9806, -2.8081, -2.7510, -2.8542,
     & -2.8819, -3.0785, -2.9233, -2.8371, -2.6235, -2.5790,
     & -2.7472, -2.8526, -2.8227, -2.8546, -2.8490, -2.9192,
     & -2.9889, -3.0270, -3.0853, -3.0810, -3.1262, -3.1026,
     & -3.0409, -2.9590, -3.0391, -3.0472, -3.0685, -3.0575,
     & -2.9761, -2.9495, -3.0199, -3.0508, -3.0091, -3.0622,
     & -3.0736, -3.0554, -3.0423, -3.0434, -2.9917, -2.9872,
     & -3.0008, -3.0047, -2.9845, -2.9922, -2.9961, -3.0097/
c PAIR OA NA
      data (pmf(11,11,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.5962, -2.0964, -2.8278, -2.5757, -3.0732, -3.2348,
     & -3.2921, -2.6802, -2.6116, -2.4360, -2.3705, -2.7311,
     & -2.8828, -2.5071, -2.5319, -2.7898, -2.8249, -2.5304,
     & -2.7842, -2.7478, -3.0875, -3.1792, -3.2007, -3.1943,
     & -3.0842, -3.1324, -3.1282, -3.1122, -3.0851, -2.9976,
     & -3.0088, -3.1504, -3.0686, -3.0445, -3.0925, -2.9749,
     & -2.8989, -2.8771, -3.0229, -3.1389, -3.0466, -3.0779/
c PAIR OA ND
      data (pmf(11,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.9996,
     & -2.5033, -3.1876, -3.7225, -3.2488, -2.8048, -2.2467,
     & -2.6187, -2.8759, -2.6945, -2.6819, -2.5362, -2.7976,
     & -3.0197, -2.6675, -2.9948, -2.8124, -2.8295, -2.7952,
     & -2.7827, -2.8649, -2.8850, -2.9010, -2.9839, -2.9022,
     & -2.9906, -2.9214, -3.1445, -2.9554, -3.0269, -2.9783,
     & -3.0402, -3.0324, -3.0052, -3.0219, -3.0860, -3.0328,
     & -3.1142, -3.0408, -3.0203, -3.0566, -2.9848, -3.0133,
     & -3.0675, -3.0677, -3.0470, -3.0977, -3.0170, -3.0289/
c PAIR OA NR
      data (pmf(11,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.4524, -1.0951, -1.6360,
     & -2.0470, -2.2510, -2.6016, -2.1734, -2.6109, -2.9393,
     & -3.1072, -3.0482, -3.0994, -2.9718, -2.7950, -2.8604,
     & -3.0007, -2.9711, -2.9054, -2.8372, -2.9051, -2.9448,
     & -3.0121, -3.0235, -3.0070, -3.0615, -3.0941, -3.1019,
     & -3.0912, -3.1144, -3.0962, -3.0870, -3.0949, -3.0855,
     & -3.0405, -3.0243, -3.0043, -2.9919, -2.9679, -2.9805,
     & -3.0099, -2.9977, -2.9975, -3.0538, -3.0274, -3.0234,
     & -3.0080, -2.9838, -2.9899, -2.9780, -2.9486, -2.9756/
c PAIR OA N0
      data (pmf(11,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OA NS
      data (pmf(11,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -2.2420, -2.6346,
     & -2.6617, -2.2022, -2.7829, -2.4627, -2.5522, -2.6578,
     & -2.3304, -2.4858, -2.5086, -2.8292, -2.8222, -2.6049,
     & -2.6809, -2.7219, -2.9260, -2.6701, -2.8458, -2.8874,
     & -2.7649, -2.9250, -2.7696, -2.8623, -3.0628, -3.0627,
     & -2.9975, -2.9554, -3.1043, -3.0084, -2.9849, -3.0293,
     & -2.9856, -3.0422, -3.0333, -3.0814, -3.0522, -3.0605,
     & -3.0831, -3.0682, -3.1434, -3.0739, -3.0697, -3.0407,
     & -3.0208, -3.0203, -3.0716, -3.0110, -2.9782, -2.9862/
c PAIR OA OC
      data (pmf(11,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.8307,  0.0000, -0.5521, -1.4912, -2.2018,
     & -2.1821, -2.1000, -1.8922, -2.2109, -2.3312, -2.4802,
     & -2.4967, -2.7025, -2.6852, -2.6914, -2.6604, -2.7753,
     & -2.8938, -2.9136, -2.9121, -2.8379, -2.7946, -2.8782,
     & -2.9423, -2.9571, -2.9599, -2.8797, -2.9270, -2.9463,
     & -2.9780, -3.0008, -3.0198, -3.0352, -2.9974, -3.0508,
     & -3.0569, -3.0422, -3.0814, -3.0555, -3.0338, -3.0685,
     & -3.0594, -3.0548, -3.0483, -2.9992, -3.0142, -3.0538,
     & -3.0600, -3.0419, -3.0258, -3.0475, -3.0137, -3.0327/
c PAIR OA OA
      data (pmf(11,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.6905,
     & -2.1192, -2.3781, -2.1165, -2.4548, -2.7827, -3.0223,
     & -3.0385, -2.9936, -2.8823, -2.8685, -2.9169, -3.0886,
     & -3.1978, -3.1590, -3.0146, -2.9030, -2.7845, -2.8254,
     & -2.7900, -2.8436, -2.8593, -2.8787, -2.9350, -3.0085,
     & -2.9907, -2.9948, -2.9998, -3.0121, -3.0221, -3.0235,
     & -3.0288, -3.0370, -3.0259, -3.0517, -3.0321, -3.0755,
     & -3.0186, -3.0067, -3.0278, -3.0214, -3.0216, -3.0715,
     & -3.0336, -3.0296, -3.0446, -3.0153, -3.0127, -2.9995/
c PAIR OA OE
      data (pmf(11,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.4104, -1.9275, -1.1613, -2.3504,
     & -2.3694, -2.0608, -1.9942, -2.4817, -2.6749, -2.8099,
     & -2.9870, -3.0228, -2.9680, -2.9590, -2.9252, -2.8909,
     & -3.0036, -2.8806, -2.8047, -2.8815, -2.9728, -2.9037,
     & -2.8668, -2.7982, -2.8647, -2.8913, -2.9052, -2.9102,
     & -2.9633, -3.0177, -3.0084, -3.0667, -3.1178, -3.0825,
     & -3.0652, -3.0640, -3.0526, -3.0455, -3.0673, -3.0595,
     & -3.0432, -3.0323, -3.0014, -3.0111, -2.9945, -3.0409,
     & -3.0300, -3.0119, -3.0208, -3.0172, -3.0188, -3.0070/
c PAIR OA OR
      data (pmf(11,19,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.7990, -2.3314, -3.3914,
     & -3.3816, -3.3233, -2.7270, -3.0919, -3.3239, -2.9977,
     & -3.2043, -3.0880, -2.9726, -3.1493, -3.1027, -3.0516,
     & -2.9098, -2.8909, -2.9153, -2.9554, -2.8444, -2.8587,
     & -2.9403, -2.9252, -3.0177, -3.0540, -3.2570, -3.1871,
     & -3.0076, -3.0131, -3.0747, -2.9508, -3.0146, -3.0958,
     & -3.0803, -3.0322, -3.0338, -3.0250, -3.1074, -2.9767,
     & -3.0048, -2.9455, -2.8359, -2.9085, -2.9709, -2.9509/
c PAIR OA OS
      data (pmf(11,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.3550, -0.2075, -0.4854, -1.5154, -2.3049,
     & -2.1844, -1.8954, -1.9355, -2.0633, -2.1782, -2.2342,
     & -2.4204, -2.4905, -2.6667, -2.8245, -2.8780, -2.9629,
     & -3.1203, -3.0648, -3.0158, -2.9929, -3.0172, -3.0481,
     & -3.0866, -3.0722, -3.0093, -3.0324, -3.0685, -3.0976,
     & -3.0880, -3.0712, -3.0247, -3.0111, -2.9976, -3.0112,
     & -3.0008, -2.9707, -2.9886, -3.0311, -3.0011, -3.0305,
     & -3.0156, -3.0360, -3.0010, -3.0132, -3.0266, -3.0220,
     & -3.0084, -2.9895, -2.9989, -3.0231, -3.0118, -3.0000/
c PAIR OA OD
      data (pmf(11,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -0.4514, -0.9521, -1.1224, -1.8223, -2.1240,
     & -2.4784, -2.9520, -2.7706, -2.5929, -2.6752, -2.6021,
     & -2.6233, -2.7069, -2.7044, -2.8512, -2.8518, -2.9106,
     & -2.9263, -2.9057, -2.8398, -2.8505, -2.8914, -2.8416,
     & -2.8915, -2.9598, -2.9249, -2.9636, -3.0058, -3.0094,
     & -3.0418, -3.0441, -3.0244, -3.0415, -3.0409, -3.0264,
     & -3.0314, -3.0290, -3.0462, -3.0304, -3.0236, -3.0495,
     & -3.0423, -3.0472, -3.0462, -3.0458, -3.0053, -3.0282,
     & -3.0378, -3.0308, -3.0152, -3.0177, -3.0059, -3.0002/
c PAIR OA P
      data (pmf(11,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.5731,  0.0000, -1.8646, -1.8649,
     & -0.6840,  0.0000, -0.8699, -1.5078, -1.4123, -1.5850,
     & -1.8485, -2.0409, -2.1993, -2.4796, -2.4801, -2.9047,
     & -2.9424, -3.0034, -2.9743, -3.0543, -3.2191, -3.2631,
     & -3.3260, -3.1592, -3.0818, -3.0273, -3.0606, -3.0406,
     & -3.0829, -2.9610, -2.9500, -2.9870, -3.0009, -3.0099,
     & -2.9999, -3.0422, -3.0306, -3.0342, -2.9548, -2.9632,
     & -2.9298, -2.9689, -3.0228, -3.0326, -3.0141, -3.0507,
     & -3.0291, -3.0055, -3.0531, -3.0828, -3.0611, -3.0076/
c PAIR OA SA
      data (pmf(11,23,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.3113, -2.6240,  0.0000,
     &  0.0000, -1.9393, -3.0043, -2.7265, -2.5200, -2.8572,
     & -2.7094, -2.5521, -2.7274, -2.7803, -2.8595, -3.1997,
     & -2.9147, -2.9886, -3.0437, -3.1767, -3.0074, -2.9708,
     & -2.9558, -3.2175, -3.0291, -2.9615, -2.8747, -3.0158,
     & -3.1895, -3.1071, -2.9500, -2.9606, -3.0942, -3.1129,
     & -3.0605, -2.9883, -3.0713, -2.9094, -3.0218, -2.9913,
     & -3.0544, -3.0883, -3.0404, -3.0531, -2.9342, -3.0192/
c PAIR OA SD
      data (pmf(11,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OA HL
      data (pmf(11,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.3023, -2.8356, -2.7781, -2.8548,
     & -2.5569, -2.4789, -2.6757, -2.5741, -2.6664, -2.7394,
     & -2.7665, -2.6861, -2.6862, -2.7580, -2.7142, -2.8237,
     & -2.8880, -2.8340, -2.8872, -2.8594, -2.8923, -2.9123,
     & -2.9635, -2.9183, -2.9529, -2.9910, -2.9810, -3.0411,
     & -3.0022, -3.0221, -3.0333, -3.0315, -3.0289, -3.0624,
     & -3.0529, -3.0503, -3.0879, -3.0324, -3.0314, -3.0531,
     & -3.0575, -3.0390, -3.0156, -3.0485, -3.0571, -3.0528,
     & -3.0601, -3.0276, -3.0339, -3.0282, -3.0155, -2.9652/
c PAIR OA F
      data (pmf(11,26,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.9044, -2.2239, -1.7261,
     & -2.2900, -2.6186, -2.7117, -2.2294, -1.9875, -2.5665,
     & -2.6182, -2.7752, -2.6345, -2.7214, -2.8790, -2.5077,
     & -2.7637, -2.8466, -2.7433, -3.1980, -2.9295, -2.9794,
     & -3.0389, -3.1022, -3.1374, -3.0173, -3.1353, -3.1361,
     & -3.0246, -3.1231, -3.0832, -3.0788, -3.1924, -3.0528,
     & -3.0158, -3.0682, -3.0532, -2.9637, -2.9848, -2.9165,
     & -3.1286, -3.0204, -3.0235, -3.0261, -3.0282, -2.9664/
c PAIR OD CF
      data (pmf(12, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.8992,
     & -2.6360, -2.2142, -2.3516, -2.8909, -3.0994, -3.2733,
     & -3.3099, -3.4369, -3.3717, -3.2072, -3.2542, -3.2376,
     & -3.2445, -3.1115, -3.1000, -3.1280, -3.0840, -3.0508,
     & -3.1185, -3.1197, -3.1301, -3.0879, -3.0801, -3.0658,
     & -3.0422, -3.0442, -3.0225, -3.0142, -2.9886, -2.9997,
     & -2.9838, -2.9681, -3.0203, -3.0225, -2.9782, -2.9329,
     & -2.9833, -2.9161, -2.9433, -2.9521, -2.9796, -2.9626,
     & -2.9778, -3.0059, -2.9947, -2.9651, -2.9551, -2.9608/
c PAIR OD CP
      data (pmf(12, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.3302, -1.2332, -1.7502, -2.2804, -1.6947,
     & -2.5905, -2.1712, -2.1679, -2.8113, -3.2909, -3.5117,
     & -3.5060, -3.4724, -3.3799, -3.3572, -3.3123, -3.3888,
     & -3.3569, -3.2648, -3.2436, -3.1351, -3.1427, -3.1454,
     & -3.1011, -3.0719, -3.0382, -2.9959, -2.9695, -2.9954,
     & -3.0035, -2.9928, -2.9718, -2.9483, -3.0447, -3.0272,
     & -2.9777, -2.9802, -2.9785, -3.0498, -3.0568, -2.9482,
     & -2.9746, -2.9438, -2.9684, -2.9472, -2.9971, -3.0233,
     & -2.9695, -2.9264, -2.9267, -2.9847, -2.9294, -2.9288/
c PAIR OD cF
      data (pmf(12, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.0972, -2.9157, -2.4963, -2.5074,
     & -2.6405, -2.8605, -3.0424, -3.0895, -3.2239, -3.5097,
     & -3.4864, -3.4684, -3.4772, -3.5152, -3.3745, -3.4496,
     & -3.3431, -3.2408, -3.2100, -3.1347, -3.1709, -3.1771,
     & -3.1618, -3.0435, -3.0790, -3.0234, -3.0926, -3.0927,
     & -3.0996, -3.1133, -3.0929, -3.0373, -3.0093, -3.0641,
     & -2.9832, -2.9977, -2.9328, -2.9479, -2.9461, -2.9259,
     & -2.9406, -2.9078, -2.9045, -2.9193, -2.9182, -2.9151,
     & -2.9210, -2.9662, -2.9224, -2.9536, -2.9618, -2.9442/
c PAIR OD cP
      data (pmf(12, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.7508, -2.4212, -2.7686, -2.5811, -2.5433,
     & -2.8751, -2.6501, -2.8028, -3.0172, -3.3510, -3.6024,
     & -3.5200, -3.4394, -3.4120, -3.4144, -3.3476, -3.3060,
     & -3.3024, -3.2464, -3.3492, -3.2972, -3.2266, -3.1013,
     & -3.1416, -3.1094, -3.0456, -3.0928, -3.0059, -3.1196,
     & -3.0129, -3.0614, -3.0910, -3.0636, -3.0550, -3.0767,
     & -3.0519, -3.0249, -3.0474, -2.9419, -2.9719, -2.9837,
     & -2.8764, -2.9357, -2.9518, -2.9432, -2.8746, -2.9230,
     & -2.9580, -2.9086, -2.8866, -2.8928, -2.9282, -2.9243/
c PAIR OD C3
      data (pmf(12, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.4613,  0.0000,  0.0000, -2.4737, -2.7753,
     & -2.9760, -2.6571, -2.8914, -3.3394, -3.3767, -3.4057,
     & -3.3745, -3.5152, -3.3675, -3.2870, -3.3177, -3.3573,
     & -3.2758, -3.0378, -3.1291, -3.1606, -3.1206, -3.2159,
     & -3.2205, -3.0212, -3.0834, -3.0825, -3.0610, -3.0484,
     & -3.0935, -3.1979, -3.0606, -3.0695, -3.0203, -2.8692,
     & -2.9967, -3.0395, -3.0767, -3.0090, -3.0761, -2.9670,
     & -2.9528, -2.9622, -2.8678, -2.8293, -2.8949, -2.9370,
     & -2.8730, -2.9094, -2.9860, -2.9866, -2.9790, -2.9542/
c PAIR OD CW
      data (pmf(12, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.2769,  0.0000,  0.0000, -2.2893,  0.0000,
     & -2.4548, -2.7374, -2.8588, -3.2046, -3.4500, -3.6270,
     & -3.4609, -3.5312, -3.3932, -3.4735, -3.4546, -3.2457,
     & -3.1684, -3.2459, -3.0800, -3.0666, -3.1026, -3.2114,
     & -3.1003, -3.1074, -3.0106, -2.8910, -2.9186, -3.0263,
     & -3.0183, -2.9143, -3.0065, -3.0417, -2.9938, -3.1160,
     & -2.9232, -2.9560, -3.0355, -2.9823, -2.9337, -3.0032,
     & -2.9456, -2.9995, -3.0396, -3.0346, -2.9845, -3.0614,
     & -2.9960, -2.9968, -2.9990, -2.8984, -2.8792, -2.8448/
c PAIR OD CO
      data (pmf(12, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.7860, -3.2336, -2.7067, -2.5994,
     & -2.0729, -2.6117, -3.1576, -3.2312, -3.7283, -4.0690,
     & -3.9470, -3.5137, -3.2488, -3.2700, -3.1562, -3.1267,
     & -3.1066, -3.1261, -3.1578, -3.0618, -3.2716, -3.2117,
     & -3.2873, -3.1369, -3.0321, -2.9645, -3.0048, -3.0034,
     & -2.9589, -3.0630, -3.0619, -3.1028, -3.0008, -3.0044,
     & -2.9022, -2.8180, -2.8540, -2.9229, -2.9346, -2.7952,
     & -3.0156, -3.0339, -3.0270, -3.0457, -2.9744, -2.9880,
     & -2.9399, -2.9930, -2.8432, -2.8686, -2.9139, -2.8710/
c PAIR OD CN
      data (pmf(12, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.7470, -2.3898, -3.6108,
     & -3.3447, -2.7774, -2.7461, -2.7744, -3.0557, -3.1667,
     & -3.4733, -3.5422, -3.6843, -3.5994, -3.5623, -3.3798,
     & -3.3718, -3.2068, -3.2013, -3.0770, -3.2305, -3.2227,
     & -3.1832, -3.0577, -3.1577, -2.9883, -2.9346, -3.0915,
     & -2.9530, -3.1595, -3.1839, -3.0359, -3.0417, -3.1061,
     & -3.1226, -2.9643, -3.0290, -3.0306, -3.1049, -2.9120,
     & -2.9013, -2.8997, -2.8976, -2.8388, -2.8497, -2.9004,
     & -2.8589, -2.9039, -2.9723, -2.7565, -2.8928, -2.9062/
c PAIR OD NC
      data (pmf(12, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.3308,  0.0000,  0.0000,  0.0000, -1.8269,
     & -2.7715, -3.1668, -3.7446, -3.7735, -3.5792, -3.7170,
     & -3.4072, -3.1621, -3.2245, -3.2166, -3.0532, -3.0786,
     & -3.1353, -3.2063, -3.4161, -3.4000, -3.1879, -3.1286,
     & -3.0768, -3.0807, -3.0427, -3.1143, -3.1767, -3.0548,
     & -3.0287, -2.9055, -2.9986, -3.0298, -3.0520, -3.2039,
     & -3.1651, -3.1793, -3.1387, -3.0129, -2.9530, -2.9042,
     & -2.8283, -2.9456, -3.0129, -2.9355, -2.8984, -2.7753,
     & -2.7685, -2.8897, -2.8686, -2.9012, -2.9098, -2.8608/
c PAIR OD NP
      data (pmf(12,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.4729,
     & -2.2921, -2.9469, -3.2071, -3.4005, -3.4843, -3.4308,
     & -3.5531, -3.6925, -3.4493, -3.3741, -3.4352, -3.2675,
     & -3.4031, -3.3990, -3.3456, -3.2202, -3.2182, -3.3040,
     & -3.3143, -3.1744, -2.9664, -3.0189, -3.0654, -3.0177,
     & -3.0340, -3.0671, -2.9574, -2.9149, -3.0776, -3.0527,
     & -2.9151, -2.8165, -2.8717, -2.7991, -2.9200, -2.9311,
     & -2.8731, -2.9123, -2.9484, -3.0511, -3.0445, -2.9220,
     & -3.0186, -2.9180, -2.8620, -2.9983, -2.9761, -2.9814/
c PAIR OD NA
      data (pmf(12,11,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OD ND
      data (pmf(12,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -2.5875,  0.0000,
     & -2.3358, -2.6111, -2.7262, -3.2583, -2.9102, -3.3374,
     & -3.0480, -3.3003, -3.0737, -3.0851, -3.2460, -3.2010,
     & -3.1013, -2.9622, -3.1350, -2.7460, -2.9312, -3.0580,
     & -2.8310, -2.6393, -2.7027, -3.1138, -2.9633, -3.3193,
     & -3.0906, -3.0523, -3.0658, -2.9533, -2.7427, -3.0849,
     & -2.9728, -2.9351, -2.9860, -3.2194, -3.2285, -3.1977,
     & -2.8840, -2.9739, -2.7922, -3.1064, -3.0339, -2.9888,
     & -2.9326, -3.1893, -2.8763, -2.9785, -2.8902, -2.9353/
c PAIR OD NR
      data (pmf(12,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -1.8156, -2.9352,
     & -2.8055, -3.0559, -3.2447, -3.0362, -3.2209, -3.2974,
     & -3.2621, -3.3904, -3.4256, -3.5061, -3.3624, -3.3659,
     & -3.2619, -3.2050, -3.1920, -3.2364, -3.2300, -3.0916,
     & -2.9831, -3.0770, -3.2322, -3.1829, -3.0773, -3.1225,
     & -3.1197, -3.0556, -3.0557, -3.1039, -3.0882, -3.1518,
     & -3.0232, -3.0232, -2.9948, -2.9894, -3.0801, -3.0546,
     & -2.9354, -2.8665, -2.9063, -2.8743, -2.8897, -2.8921,
     & -2.8396, -2.8932, -2.9078, -2.8696, -2.8909, -2.9669/
c PAIR OD N0
      data (pmf(12,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OD NS
      data (pmf(12,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -3.2579, -2.7178,  0.0000,  0.0000,
     & -3.0195, -3.0765, -3.7508, -3.2808, -3.3752, -3.2260,
     & -3.3364, -3.2088, -3.0816, -3.2762, -3.3138, -3.5970,
     & -3.2671, -3.2426, -3.0417, -3.0732, -2.8353, -3.1020,
     & -2.9426, -3.0634, -3.1618, -2.8935, -2.9254, -2.9509,
     & -3.0228, -3.1844, -2.9534, -3.0331, -2.8652, -3.0295,
     & -3.0606, -3.1303, -3.2053, -2.8961, -2.9895, -2.8432,
     & -2.8363, -3.1417, -2.8808, -3.0325, -3.0836, -3.0020,
     & -2.8253, -3.1660, -2.8108, -2.5816, -3.0711, -2.8098/
c PAIR OD OC
      data (pmf(12,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.7861, -3.0639, -3.2171, -3.5636,
     & -3.7063, -4.1192, -3.6841, -3.3784, -3.4655, -3.4447,
     & -3.1385, -3.2463, -3.3136, -3.3233, -3.2972, -3.3552,
     & -3.3748, -3.2410, -3.2033, -3.1555, -3.0332, -2.9752,
     & -2.9612, -3.0394, -3.0261, -3.0679, -3.0222, -2.9428,
     & -3.0476, -3.0753, -3.0364, -2.9614, -2.9339, -2.8274,
     & -2.8661, -2.8790, -3.0079, -3.0235, -3.0341, -3.0439,
     & -2.9860, -2.9750, -3.0036, -3.0032, -2.9235, -2.9073,
     & -2.8670, -2.8820, -2.9423, -2.9185, -2.9600, -2.9739/
c PAIR OD OA
      data (pmf(12,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.8795, -2.8055, -3.0964, -3.2284,
     & -3.6022, -4.0301, -3.5061, -3.3238, -3.1130, -3.2128,
     & -3.4618, -3.4893, -3.3742, -3.1417, -3.1223, -2.9540,
     & -2.8678, -3.1283, -2.7908, -2.7838, -2.9422, -3.0184,
     & -3.0410, -3.0931, -3.0553, -3.0741, -3.0601, -3.1600,
     & -3.0863, -3.0910, -3.0799, -3.0525, -3.0317, -3.1384,
     & -2.9721, -2.9011, -2.9274, -3.0229, -2.9965, -2.9486,
     & -2.8420, -2.9342, -3.0577, -2.9020, -2.9369, -2.9608,
     & -2.9244, -2.9535, -2.9746, -3.0375, -3.0226, -2.9391/
c PAIR OD OE
      data (pmf(12,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.6091, -2.8706, -2.3305, -2.6215, -2.5142,
     & -3.1026, -2.4639, -2.7441, -3.3632, -3.5174, -3.0094,
     & -3.4745, -3.3309, -3.2745, -2.9814, -3.2602, -3.4121,
     & -3.4474, -3.4194, -3.4260, -3.3677, -3.2251, -3.2721,
     & -3.0891, -3.1190, -3.0661, -2.9687, -2.9866, -2.9510,
     & -3.0786, -3.1727, -3.0613, -3.1068, -2.9710, -2.9480,
     & -3.0401, -2.9793, -3.0275, -3.0389, -2.9404, -2.9230,
     & -2.8431, -2.9765, -2.9159, -2.9380, -2.9585, -3.0051,
     & -2.9422, -2.9461, -2.9535, -2.9282, -2.8455, -2.9499/
c PAIR OD OR
      data (pmf(12,19,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OD OS
      data (pmf(12,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.6052, -2.2132, -1.8967, -2.8265,
     & -3.8516, -4.0780, -3.7695, -3.2575, -3.3647, -3.2827,
     & -3.3623, -3.3367, -3.4149, -3.4666, -3.5266, -3.5045,
     & -3.5174, -3.3796, -3.1632, -3.0112, -2.9867, -2.8708,
     & -2.9379, -3.0890, -3.1265, -3.0236, -3.0158, -3.0480,
     & -2.9895, -2.9677, -3.0120, -2.9678, -3.0227, -2.9998,
     & -3.0331, -2.9713, -2.9899, -2.9826, -2.9618, -2.9847,
     & -2.9732, -2.9306, -2.9340, -2.8684, -2.9759, -2.9514,
     & -2.9193, -2.9392, -2.9533, -2.8986, -2.8791, -2.8835/
c PAIR OD OD
      data (pmf(12,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.3689,  0.0000, -1.5288, -2.3710,
     & -3.1913, -3.5955, -3.5180, -3.1802, -3.1647, -3.1449,
     & -3.1064, -3.3831, -3.3395, -3.4032, -3.2732, -3.2396,
     & -3.2847, -3.2395, -3.1364, -3.1063, -3.0340, -3.0812,
     & -2.9916, -3.0637, -2.9989, -2.9756, -3.0385, -3.0333,
     & -2.9957, -3.0209, -2.9272, -3.0086, -3.0157, -3.0043,
     & -3.0806, -3.0162, -3.0145, -3.0323, -2.9901, -2.9461,
     & -2.9156, -2.9720, -2.9362, -2.9923, -2.9807, -2.9209,
     & -2.9546, -3.0015, -3.0003, -2.9985, -2.9180, -2.9771/
c PAIR OD P
      data (pmf(12,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.1549, -2.4327, -1.9058, -2.4467,
     & -2.3228, -1.9684, -1.4517, -2.7089, -3.5630, -3.9552,
     & -4.1705, -4.1129, -3.6405, -3.4358, -3.3465, -3.2370,
     & -3.4275, -3.3034, -3.3168, -3.1563, -3.2636, -3.0773,
     & -2.9187, -2.8119, -2.9460, -2.9027, -3.0164, -3.0543,
     & -2.8562, -2.8186, -3.0936, -3.0882, -3.0947, -3.0922,
     & -2.9843, -3.0993, -3.0916, -2.9044, -2.9851, -2.8483,
     & -2.8629, -2.9159, -2.9641, -2.9398, -2.8703, -2.9446,
     & -2.9192, -2.9309, -2.8793, -2.9258, -2.9823, -2.9201/
c PAIR OD SA
      data (pmf(12,23,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OD SD
      data (pmf(12,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OD HL
      data (pmf(12,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.9256, -3.2092, -3.0350, -3.1670,
     & -3.3303, -3.1991, -3.0175, -3.2016, -3.3149, -3.2390,
     & -3.0867, -3.1688, -3.1240, -3.1801, -3.2793, -3.2961,
     & -3.2237, -3.2257, -3.1438, -3.1180, -3.0695, -2.9542,
     & -3.0356, -2.9705, -3.1090, -3.0261, -3.0310, -3.1224,
     & -3.1235, -3.0723, -3.0772, -3.0842, -2.9737, -3.0893,
     & -3.0163, -3.0138, -2.9953, -3.0613, -3.0258, -3.0433,
     & -3.0362, -2.9816, -2.9453, -3.0026, -2.9479, -2.8839,
     & -2.9755, -2.9607, -2.9713, -2.9704, -2.9338, -2.8684/
c PAIR OD F
      data (pmf(12,26,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OW CF
      data (pmf(13, 1,i),i=1,60) /
     &  0.0000, -2.8844,  0.0000,  0.0000,  0.0000, -1.3711,
     &  0.0000, -1.0058, -1.2673, -1.6768, -2.2808, -2.3770,
     & -2.4827, -2.5479, -2.7059, -2.7780, -3.0687, -3.3524,
     & -3.5160, -3.5270, -3.3890, -3.3165, -3.2927, -3.2762,
     & -3.2765, -3.2339, -3.1650, -3.1313, -3.0658, -3.0588,
     & -3.0611, -3.0492, -3.0274, -3.0228, -3.0127, -3.0186,
     & -3.0220, -3.0181, -3.0218, -3.0058, -2.9800, -2.9994,
     & -2.9905, -3.0025, -3.0101, -3.0091, -3.0207, -3.0093,
     & -3.0195, -2.9777, -2.9796, -2.9925, -2.9723, -2.9533,
     & -2.9275, -2.9405, -2.9605, -2.9320, -2.9601, -2.9494/
c PAIR OW CP
      data (pmf(13, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000, -1.2803, -1.4533,
     &  0.0000, -2.7051, -1.8901, -2.2755, -2.7269, -2.8033,
     & -2.6670, -2.7565, -2.6800, -2.9696, -3.3021, -3.6054,
     & -3.7440, -3.7172, -3.6227, -3.5394, -3.4901, -3.4360,
     & -3.4008, -3.3951, -3.3192, -3.2652, -3.1932, -3.1673,
     & -3.1501, -3.1462, -3.1381, -3.1128, -3.0986, -3.0903,
     & -3.0873, -3.0623, -3.0824, -3.0508, -3.0281, -3.0169,
     & -3.0225, -2.9956, -2.9878, -2.9553, -2.9843, -2.9699,
     & -2.9624, -2.9450, -2.9474, -2.9187, -2.9116, -2.8912,
     & -2.8964, -2.8958, -2.8767, -2.8368, -2.8335, -2.8225/
c PAIR OW cF
      data (pmf(13, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -0.9013, -1.8273, -2.3894, -2.8984,
     & -2.8931, -2.7217, -2.7216, -2.8260, -3.0722, -3.2271,
     & -3.3095, -3.3762, -3.3679, -3.2710, -3.2874, -3.3076,
     & -3.2525, -3.1997, -3.1770, -3.1432, -3.1327, -3.1344,
     & -3.1132, -3.0591, -3.0469, -3.0478, -3.0694, -3.0541,
     & -3.0291, -3.0333, -3.0604, -3.0695, -3.0612, -3.0427,
     & -3.0459, -3.0165, -3.0305, -3.0436, -3.0112, -2.9947,
     & -2.9752, -2.9749, -2.9724, -2.9168, -2.9241, -2.9419,
     & -2.9153, -2.9296, -2.9536, -2.9594, -2.9129, -2.9525/
c PAIR OW cP
      data (pmf(13, 4,i),i=1,60) /
     & -4.2337, -3.0856,  0.0000,  0.0000, -1.8083,  0.0000,
     & -2.0238, -2.5656, -2.1167, -2.1553, -2.4082, -2.9359,
     & -2.7965, -2.8097, -2.7202, -2.9533, -3.2041, -3.4062,
     & -3.5697, -3.5520, -3.3285, -3.2772, -3.2196, -3.2291,
     & -3.2360, -3.2658, -3.2198, -3.1793, -3.2079, -3.1827,
     & -3.1449, -3.0654, -3.0147, -3.0627, -3.0562, -3.1036,
     & -3.0886, -3.0983, -3.0802, -3.1050, -3.1170, -3.0699,
     & -3.0478, -3.0231, -3.0303, -2.9958, -3.0075, -2.9748,
     & -2.9476, -2.9571, -2.9390, -2.9374, -2.9430, -2.9199,
     & -2.9127, -2.9129, -2.9166, -2.9108, -2.8925, -2.8803/
c PAIR OW C3
      data (pmf(13, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -2.3800, -3.0118,
     & -2.8386, -2.7182, -2.6494, -2.6921, -2.7624, -3.0347,
     & -3.3568, -3.3238, -3.3987, -3.3489, -3.3519, -3.2522,
     & -3.2643, -3.3054, -3.1481, -3.1320, -3.0865, -3.0117,
     & -3.0317, -3.0413, -3.0388, -3.0262, -3.0551, -3.0898,
     & -3.0573, -2.9580, -2.9984, -3.0445, -3.0199, -2.9648,
     & -2.9210, -2.8969, -2.9721, -3.0000, -2.9769, -2.9530,
     & -2.9986, -2.9585, -3.0245, -3.0386, -2.9923, -2.9918,
     & -3.0075, -2.9100, -2.9732, -2.9850, -2.9862, -3.0188/
c PAIR OW CW
      data (pmf(13, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.1423, -2.0112, -2.0629, -2.6510,
     & -2.7069, -2.8197, -2.6565, -2.8522, -3.1914, -3.5004,
     & -3.7521, -3.7441, -3.6052, -3.4379, -3.4034, -3.3638,
     & -3.3700, -3.3691, -3.2292, -3.1533, -3.2267, -3.1564,
     & -3.1236, -3.0203, -3.0832, -3.0892, -3.1263, -3.0902,
     & -3.0757, -3.1345, -3.0883, -3.0918, -3.0333, -3.0331,
     & -3.0000, -3.0012, -2.9684, -2.9697, -2.9701, -2.9090,
     & -2.9673, -2.9842, -2.9532, -2.8995, -2.9792, -2.9117,
     & -2.9160, -2.8859, -2.9251, -2.8803, -2.8715, -2.8857/
c PAIR OW CO
      data (pmf(13, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.9810,
     &  0.0000,  0.0000, -1.8773, -1.7462, -2.5156, -2.4083,
     & -2.4605, -2.6525, -2.8275, -3.1310, -3.4493, -3.8319,
     & -3.7982, -3.6641, -3.4336, -3.3723, -3.1418, -3.1279,
     & -3.1826, -3.1099, -3.1183, -3.2012, -3.0730, -3.0753,
     & -3.0907, -3.0266, -3.0384, -3.0667, -3.0671, -3.0421,
     & -3.0607, -3.1675, -3.0285, -3.0669, -3.0629, -3.0215,
     & -3.0214, -3.0110, -3.0207, -2.9565, -2.9405, -2.9517,
     & -2.9442, -2.9113, -2.9182, -2.9490, -2.9439, -3.0147,
     & -2.9838, -2.9535, -2.8724, -2.8983, -2.9400, -2.8431/
c PAIR OW CN
      data (pmf(13, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.0763, -2.4858, -2.7145, -3.0481,
     & -2.4112, -2.7025, -2.4699, -2.8400, -3.1908, -3.4939,
     & -3.7001, -3.7554, -3.5783, -3.5020, -3.3658, -3.2769,
     & -3.2667, -3.2255, -3.1107, -3.1531, -3.1281, -3.0114,
     & -3.0181, -3.0155, -3.0020, -2.9862, -3.0025, -2.9378,
     & -3.0459, -3.0161, -3.0422, -3.0536, -3.0850, -3.0403,
     & -3.0599, -3.0459, -3.0280, -2.9654, -2.9742, -2.9557,
     & -2.9816, -2.9690, -2.9952, -2.9536, -2.9901, -2.9483,
     & -2.9272, -2.9047, -2.9023, -2.8990, -2.9319, -2.8987/
c PAIR OW NC
      data (pmf(13, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.0190,  0.0000, -2.3105, -2.5500,
     & -3.1671, -3.5345, -3.7482, -3.7230, -3.5224, -3.5298,
     & -3.6648, -3.5418, -3.4374, -3.2577, -3.1561, -3.1065,
     & -3.1394, -3.1340, -3.0087, -2.9894, -3.0355, -2.8986,
     & -3.0686, -3.0457, -3.0202, -2.9910, -3.0675, -3.0631,
     & -3.0596, -3.1064, -3.0385, -3.0304, -3.0095, -3.0632,
     & -3.0197, -3.0209, -2.9959, -2.9620, -2.9453, -2.9609,
     & -2.8986, -2.9540, -2.9718, -2.9197, -2.9707, -2.9415,
     & -2.9552, -2.9880, -2.9557, -2.8966, -2.9578, -2.9392/
c PAIR OW NP
      data (pmf(13,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.9344,  0.0000, -2.2094, -3.0278, -2.7782, -2.5012,
     & -2.6217, -3.2765, -3.1566, -2.9608, -3.2112, -3.2690,
     & -3.3336, -3.3471, -3.5131, -3.5588, -3.4701, -3.5350,
     & -3.4997, -3.3073, -3.2438, -3.1121, -3.1174, -3.1403,
     & -3.1293, -3.0363, -3.1313, -3.1950, -3.2683, -3.2343,
     & -3.1080, -3.1745, -3.1256, -3.0468, -3.0489, -3.0160,
     & -2.9921, -3.0234, -2.9825, -3.0107, -2.9864, -2.9695,
     & -2.8959, -2.9328, -2.9498, -2.9004, -2.9367, -2.8889,
     & -2.8698, -2.8976, -2.8782, -2.8946, -2.8961, -2.8372/
c PAIR OW NA
      data (pmf(13,11,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OW ND
      data (pmf(13,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.5022, -2.3842, -1.8680,
     & -2.5414, -3.1877, -3.6119, -3.3758, -3.2669, -3.2417,
     & -2.9655, -3.2404, -3.2886, -3.2337, -3.2887, -3.4436,
     & -3.4267, -3.3998, -3.3722, -3.1395, -3.1248, -3.0779,
     & -3.0826, -3.0384, -3.0400, -3.1651, -3.0529, -3.0338,
     & -3.1108, -3.1858, -3.0811, -3.0207, -3.0393, -3.0375,
     & -3.0030, -2.9575, -2.9929, -3.0193, -2.9975, -2.9585,
     & -3.0329, -2.9578, -2.9115, -2.8553, -2.8542, -2.9336,
     & -2.9662, -3.0012, -2.9449, -2.9253, -3.0057, -2.9265/
c PAIR OW NR
      data (pmf(13,13,i),i=1,60) /
     & -4.6522,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.0345, -1.4780, -2.1648, -2.0468, -2.8891,
     & -3.1382, -3.5198, -3.6398, -3.3523, -3.2405, -3.1324,
     & -3.3040, -3.3412, -3.4012, -3.3665, -3.2756, -3.2965,
     & -3.2820, -3.2746, -3.1573, -3.0532, -3.0861, -3.0851,
     & -3.1501, -3.1214, -3.0581, -3.0148, -2.9887, -3.0788,
     & -3.1081, -3.0812, -3.0535, -3.0556, -3.0059, -3.0514,
     & -3.0470, -3.0497, -3.0616, -3.0742, -3.0310, -2.9956,
     & -2.9601, -2.9169, -2.9229, -2.9110, -2.9853, -2.9314,
     & -2.8986, -2.8887, -2.9275, -2.8971, -2.9441, -2.9229/
c PAIR OW N0
      data (pmf(13,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OW NS
      data (pmf(13,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.3917, -2.2606, -2.9605, -3.3317,
     & -3.3802, -3.4730, -3.6619, -3.3340, -3.6494, -3.1475,
     & -3.4986, -3.4293, -3.3774, -3.4775, -3.2907, -3.2658,
     & -3.3635, -3.2242, -3.2065, -3.0969, -2.9531, -3.1194,
     & -3.1265, -2.9653, -3.0258, -3.0398, -2.9393, -3.2286,
     & -3.1303, -3.1229, -2.9962, -3.0098, -3.1322, -3.0942,
     & -3.1965, -3.0459, -2.9029, -2.9522, -2.9566, -2.9527,
     & -2.9123, -2.8875, -2.8519, -2.7968, -2.9410, -2.8026,
     & -2.8466, -3.0228, -2.9606, -2.9220, -2.9513, -2.9168/
c PAIR OW OC
      data (pmf(13,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.7870,  0.0000,  0.0000, -1.3398, -2.5501, -3.2760,
     & -3.5818, -3.8413, -3.8138, -3.5115, -3.3753, -3.1117,
     & -3.2327, -3.2179, -3.1962, -3.3240, -3.3185, -3.3200,
     & -3.2523, -3.1752, -3.0734, -3.0563, -3.0066, -3.0342,
     & -3.0466, -3.0176, -3.0529, -3.0673, -3.0815, -3.1019,
     & -3.1182, -3.0327, -3.0747, -3.0205, -3.0239, -3.0743,
     & -2.9694, -3.0274, -2.9739, -2.9925, -2.9265, -2.9468,
     & -2.9753, -2.9666, -2.9543, -2.9493, -2.9239, -2.9316,
     & -2.9253, -2.9493, -2.9305, -2.9486, -2.9229, -2.9114/
c PAIR OW OA
      data (pmf(13,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.4768,
     & -1.6320, -1.4634, -1.3159, -2.0027, -1.8847, -2.5573,
     & -3.3195, -3.9285, -3.9088, -3.6254, -3.3845, -3.2674,
     & -3.2780, -3.2030, -3.1773, -3.1381, -3.1841, -3.2384,
     & -3.2628, -3.1659, -3.1193, -3.0401, -3.0721, -3.0593,
     & -3.0181, -3.0555, -3.0835, -3.1541, -3.1174, -3.0561,
     & -3.1527, -3.1586, -3.0880, -3.1155, -3.0937, -3.0509,
     & -3.0214, -2.9974, -2.9928, -3.0015, -2.9827, -2.9398,
     & -2.9132, -2.9154, -2.9008, -2.8933, -2.9092, -2.9367,
     & -2.9476, -2.9197, -2.9206, -2.9378, -2.9418, -2.9439/
c PAIR OW OE
      data (pmf(13,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.7466,  0.0000, -2.6960, -2.8554,
     & -3.0862, -3.3563, -3.5617, -3.5234, -3.4248, -3.3783,
     & -3.2507, -3.3932, -3.4408, -3.3950, -3.4470, -3.4959,
     & -3.5157, -3.4486, -3.4052, -3.4154, -3.4180, -3.3299,
     & -3.2784, -3.2058, -3.0945, -3.0294, -3.0467, -3.0541,
     & -3.0996, -3.0563, -3.0429, -3.0652, -3.0618, -3.0397,
     & -3.0056, -2.9696, -3.0192, -3.0018, -2.9592, -2.9294,
     & -2.9533, -2.9258, -2.9441, -2.9241, -2.8853, -2.8918,
     & -2.9003, -2.8975, -2.8460, -2.8681, -2.8198, -2.8317/
c PAIR OW OR
      data (pmf(13,19,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -2.7194, -2.9841,  0.0000, -2.7269, -2.6111, -2.0866,
     & -2.6264, -3.5478, -3.1946, -2.8366, -3.1947, -3.8061,
     & -3.5886, -3.5533, -3.4339, -3.4290, -3.5392, -3.4579,
     & -3.3427, -3.1659, -3.2116, -3.1007, -3.1117, -2.9611,
     & -3.0191, -2.7473, -3.1182, -3.0647, -3.0213, -2.9473,
     & -2.9152, -2.8945, -3.0390, -3.0993, -2.9771, -2.8559,
     & -2.9784, -2.8504, -2.9146, -2.8194, -2.8597, -2.7984,
     & -2.8663, -2.9890, -2.9514, -2.9526, -2.9585, -2.8132/
c PAIR OW OS
      data (pmf(13,20,i),i=1,60) /
     &  0.0000, -2.6520, -2.0628,  0.0000,  0.0000, -1.1387,
     & -0.9420,  0.0000, -1.5755, -1.4444, -2.4745, -2.9563,
     & -3.6615, -4.0025, -3.7788, -3.5076, -3.3124, -3.3414,
     & -3.4338, -3.4179, -3.4366, -3.5244, -3.4173, -3.5026,
     & -3.4336, -3.4395, -3.3022, -3.2162, -3.0913, -3.0895,
     & -3.0767, -3.0988, -3.0895, -3.1363, -3.1151, -3.1540,
     & -3.1393, -3.0959, -3.0464, -3.0360, -3.0175, -3.0192,
     & -3.0108, -2.9842, -3.0022, -2.9946, -2.9755, -2.9621,
     & -2.9374, -2.8989, -2.8817, -2.8808, -2.8721, -2.8677,
     & -2.8787, -2.8676, -2.8770, -2.8467, -2.8514, -2.8160/
c PAIR OW OD
      data (pmf(13,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000, -1.7143,  0.0000,  0.0000,
     &  0.0000, -1.2270, -1.6201, -2.0528, -2.4475, -2.8123,
     & -3.4001, -3.9475, -3.8711, -3.6099, -3.4722, -3.4227,
     & -3.3086, -3.3310, -3.3691, -3.4212, -3.4210, -3.4060,
     & -3.3861, -3.3091, -3.2904, -3.2193, -3.1517, -3.1241,
     & -3.1415, -3.1112, -3.1143, -3.0845, -3.1060, -3.1027,
     & -3.1122, -3.0853, -3.0684, -3.0475, -3.0298, -3.0268,
     & -3.0067, -3.0053, -2.9994, -3.0070, -2.9815, -2.9454,
     & -2.9323, -2.9194, -2.9260, -2.9138, -2.9250, -2.9025,
     & -2.8829, -2.8842, -2.8991, -2.8842, -2.8624, -2.8379/
c PAIR OW P
      data (pmf(13,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.1031, -1.3370, -2.2868,
     & -2.9296, -2.5477, -2.5545, -2.8959, -3.4320, -3.7797,
     & -4.2006, -4.1962, -3.8419, -3.5564, -3.4355, -3.2953,
     & -3.2343, -3.1937, -3.2486, -3.1885, -3.2480, -3.2368,
     & -3.2494, -3.1208, -3.1217, -3.0677, -3.0170, -3.0556,
     & -3.0150, -3.0806, -3.1246, -3.1405, -3.0949, -3.1117,
     & -3.0412, -3.0004, -2.9880, -2.8767, -2.9124, -2.8608,
     & -2.9750, -2.9459, -2.8957, -2.8474, -2.8784, -2.9150,
     & -2.8732, -2.8860, -2.8421, -2.8353, -2.8003, -2.8520/
c PAIR OW SA
      data (pmf(13,23,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OW SD
      data (pmf(13,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR OW HL
      data (pmf(13,25,i),i=1,60) /
     & -4.5803, -3.4322,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.7222, -1.5536, -2.8722, -3.3198, -3.2018, -3.1640,
     & -3.2439, -3.3926, -3.5793, -3.5227, -3.3982, -3.4837,
     & -3.4498, -3.2072, -3.3302, -3.3719, -3.3260, -3.2443,
     & -3.2836, -3.2449, -3.2539, -3.2435, -3.2090, -3.2832,
     & -3.2407, -3.1108, -3.1079, -3.1639, -3.0666, -3.1363,
     & -3.0960, -3.0950, -3.0806, -3.1045, -3.0393, -3.0789,
     & -3.0338, -2.9376, -2.9202, -2.9709, -2.9606, -2.9582,
     & -2.9636, -2.9816, -2.9232, -2.9229, -2.9292, -2.9422,
     & -2.8820, -2.8878, -2.8965, -2.9297, -2.8671, -2.8627/
c PAIR OW F
      data (pmf(13,26,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.5920, -3.1392, -3.3458, -2.7157, -3.3660,
     & -3.4300, -3.2014, -3.4428, -2.9605, -3.1275, -2.9888,
     & -2.9917, -3.2289, -2.9841, -3.1437, -2.6754, -2.9995,
     & -2.8697, -2.7785, -2.7824, -3.1020, -2.7457, -2.9751,
     & -2.9650, -2.6768, -2.9884, -2.8072, -2.9011, -3.1384,
     & -2.8824, -3.0772, -3.0490, -3.1266, -2.9649, -3.0227,
     & -2.9559, -3.0096, -2.9599, -3.0081, -3.0677, -2.8886,
     & -3.0136, -3.0378, -3.0349, -2.9961, -3.0848, -3.1019/
c PAIR SA CF
      data (pmf(14, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.2752, -2.4188,  0.0000, -2.8877,
     & -2.8031, -3.4523, -3.6082, -3.7012, -3.1901, -3.2801,
     & -3.3633, -2.8265, -2.8853, -3.1989, -3.2253, -3.0570,
     & -3.0476, -3.1190, -2.8733, -2.8574, -2.8439, -3.0455,
     & -2.8281, -2.9722, -3.0680, -2.8425, -2.9749, -2.8927,
     & -2.9179, -2.9449, -3.0993, -2.8363, -2.7523, -2.8073,
     & -2.9272, -2.8914, -2.9814, -3.2605, -3.0086, -2.9715,
     & -2.9896, -3.0637, -3.0027, -3.0838, -3.0472, -3.0407/
c PAIR SA CP
      data (pmf(14, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.9943, -1.8679,  0.0000,  0.0000, -1.9446, -2.7268,
     & -2.7459, -2.9087, -2.8532, -3.0496, -2.7751, -3.0601,
     & -3.1753, -2.9810, -2.8829, -3.2627, -3.1822, -3.0149,
     & -2.9836, -2.8610, -2.7273, -2.9221, -2.9913, -2.9469,
     & -3.0504, -3.0316, -3.0869, -2.8313, -2.7986, -3.0544,
     & -3.0366, -3.1462, -2.9939, -3.0193, -2.9119, -2.8675,
     & -2.9475, -3.0681, -3.2513, -3.0195, -3.1099, -3.1136,
     & -2.9623, -2.8760, -3.0426, -3.1147, -2.9298, -2.9497/
c PAIR SA cF
      data (pmf(14, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -2.2580, -2.5698,
     & -2.8079, -3.0384, -3.1944, -3.1148, -3.1544, -3.0227,
     & -3.1489, -3.0126, -3.2376, -2.9019, -2.9218, -2.8916,
     & -2.7166, -2.9012, -2.5212, -2.7945, -2.8868, -2.8864,
     & -2.8490, -3.0343, -3.0616, -2.9555, -3.0593, -2.9395,
     & -2.8500, -2.9293, -3.1449, -3.0302, -3.0139, -2.9816,
     & -3.0463, -3.0484, -2.9213, -2.9527, -2.9652, -2.9824,
     & -3.0248, -3.0151, -3.0014, -3.0498, -3.2129, -3.2047/
c PAIR SA cP
      data (pmf(14, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.9066, -2.6187, -2.6426,
     & -3.3287, -3.5760, -3.6503, -3.5636, -3.3129, -3.1957,
     & -3.1360, -3.0433, -3.0076, -2.9706, -2.8860, -2.9007,
     & -2.7176, -2.5676, -2.7610, -2.6693, -2.5606, -2.7857,
     & -2.8462, -2.9103, -2.9678, -2.6659, -2.9165, -2.9512,
     & -2.8466, -2.9768, -2.7512, -2.9678, -2.9794, -2.9802,
     & -3.0373, -3.0426, -3.1164, -3.0893, -3.1631, -3.0910,
     & -3.0515, -3.1231, -3.0060, -3.0474, -3.1188, -3.0825/
c PAIR SA C3
      data (pmf(14, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.5443,  0.0000, -2.7552,
     & -3.3100, -3.8359, -3.6776, -3.6709, -3.4321, -3.2002,
     & -3.5916, -3.1117, -3.1421, -3.2659, -2.9610, -3.0820,
     & -2.6953, -2.5800, -2.6672, -3.0028, -2.8619, -3.2854,
     & -3.0253, -2.8175, -2.6320, -3.0278, -3.0147, -2.9824,
     & -3.3185, -2.9381, -3.1423, -2.6516, -2.7724, -2.7230,
     & -2.7608, -2.9804, -2.9808, -2.9681, -2.7832, -2.9547,
     & -3.0405, -3.1304, -3.0499, -2.9923, -3.0998, -3.0949/
c PAIR SA CW
      data (pmf(14, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.2340,
     & -2.5425, -3.0956, -3.5052, -3.4601, -3.0739, -2.8448,
     & -3.1327, -2.6025, -3.1796, -3.1716, -2.9725, -2.8850,
     & -2.5978, -3.2375, -3.0453, -2.9999, -3.2226, -2.8976,
     & -2.8840, -2.9542, -2.8555, -2.7522, -2.7665, -2.6870,
     & -2.9416, -3.1646, -3.0624, -3.0872, -2.7942, -2.6623,
     & -3.0988, -3.3185, -3.0400, -2.9603, -2.9160, -2.9946,
     & -3.1010, -3.0685, -2.9666, -3.0600, -3.0266, -3.0985/
c PAIR SA CO
      data (pmf(14, 7,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SA CN
      data (pmf(14, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -2.3416,
     & -2.8985, -3.1131, -3.2289, -2.8198, -3.6302, -3.5883,
     & -3.0071, -3.7188, -2.8734, -1.9913, -2.5776, -2.7713,
     & -2.7272, -2.9858, -2.9021, -2.5995, -2.4919, -2.1774,
     & -2.4880, -2.5145, -2.7570, -2.7936, -2.8237, -3.1054,
     & -2.7312, -2.7312, -3.0654, -3.0044, -2.8452, -3.1360,
     & -2.8701, -3.1648, -3.1074, -2.7591, -3.1727, -3.1017,
     & -3.2155, -3.0351, -3.2390, -3.1657, -3.2329, -2.8107/
c PAIR SA NC
      data (pmf(14, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -2.4483, -2.3601,
     & -3.5727, -3.2537,  0.0000, -3.1047, -2.7951, -3.3244,
     & -2.9927, -3.2952, -3.3573, -3.3333, -3.3713, -2.8532,
     & -2.1595, -2.1152, -2.7196, -2.1986, -2.9793, -2.9430,
     & -3.0393, -3.0329, -2.8381, -2.9929, -2.9606, -3.0276,
     & -3.1548, -2.6796, -2.9598, -3.0679, -3.0405, -2.8774,
     & -2.7563, -2.8878, -2.9010, -3.2856, -3.1825, -3.1810,
     & -2.9826, -2.8839, -2.8930, -3.1116, -3.1651, -3.0091/
c PAIR SA NP
      data (pmf(14,10,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SA NA
      data (pmf(14,11,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SA ND
      data (pmf(14,12,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SA NR
      data (pmf(14,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.5554, -2.4428,  0.0000, -2.8858,  0.0000,
     & -3.3359, -3.4144, -3.6577, -3.5462, -3.5963, -3.5414,
     & -3.2817, -3.3227, -3.0016, -2.8927, -2.8270, -2.6356,
     & -2.6866, -2.5931, -2.3776, -2.6726, -2.8214, -2.8588,
     & -2.8449, -2.7609, -2.8947, -2.5840, -2.8816, -2.9574,
     & -3.0520, -2.6768, -2.7593, -3.0689, -3.1539, -2.8878,
     & -3.0927, -2.9911, -3.0429, -2.8851, -2.9086, -3.1896,
     & -3.0722, -3.0539, -3.0820, -3.0248, -3.0468, -3.2641/
c PAIR SA N0
      data (pmf(14,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SA NS
      data (pmf(14,15,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SA OC
      data (pmf(14,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.7394,  0.0000, -2.1551,  0.0000,
     &  0.0000, -2.3229,  0.0000,  0.0000, -2.8502, -2.2840,
     & -2.6292, -3.1083, -2.0994, -2.0412, -2.5636, -2.5173,
     & -2.5419, -2.2592, -2.7947, -2.3437, -2.4373, -2.8092,
     & -2.7731, -2.9143, -3.1259, -3.0000, -2.9683, -3.0012,
     & -2.8721, -3.0120, -3.1153, -3.0429, -3.1678, -2.9888,
     & -3.0870, -3.1156, -3.1161, -3.2858, -3.1153, -3.2337,
     & -2.9989, -3.0061, -3.0316, -3.0812, -3.1471, -3.1213/
c PAIR SA OA
      data (pmf(14,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.3838,  0.0000, -2.2064,
     & -2.5288, -2.6856, -3.5566, -3.1105, -3.4073, -3.0348,
     & -3.1362, -3.1177, -2.9558, -3.2158, -3.2420, -2.7266,
     & -2.8326, -2.7853, -2.5833, -2.7360, -2.8457, -2.7386,
     & -3.2572, -3.1101, -3.0915, -3.2565, -3.3027, -2.8618,
     & -2.8903, -2.9139, -3.0234, -2.9372, -2.8454, -2.8344,
     & -2.9121, -3.1211, -2.9621, -2.9489, -3.0178, -3.0051,
     & -2.9337, -2.9759, -2.9548, -3.1546, -3.1149, -3.0817/
c PAIR SA OE
      data (pmf(14,18,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SA OR
      data (pmf(14,19,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SA OS
      data (pmf(14,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -1.9008,  0.0000, -1.7244, -2.2871,
     &  0.0000, -1.4803, -2.7649, -2.6320, -2.4143, -2.3476,
     & -2.4939, -2.6713, -2.4800, -2.9615, -2.9238, -3.0700,
     & -3.0246, -2.7274, -2.5676, -2.4426, -2.4044, -2.7612,
     & -2.9591, -2.8701, -2.9278, -2.7561, -2.7240, -2.9819,
     & -2.7714, -2.7891, -2.9486, -3.0251, -2.8732, -2.9710,
     & -3.0706, -3.1339, -2.9117, -2.9832, -3.0951, -3.0802,
     & -3.0762, -3.2369, -3.2557, -3.2398, -3.2194, -3.2066/
c PAIR SA OD
      data (pmf(14,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -1.7650, -1.6728,
     & -2.7334, -3.1379, -3.1901, -2.9038, -3.2980, -3.1869,
     & -2.8686, -3.3366, -3.0079, -3.1423, -3.1714, -2.8711,
     & -2.6538, -2.5001, -2.6526, -2.7048, -3.1337, -3.3039,
     & -3.1610, -3.2780, -2.9243, -2.8903, -3.0529, -2.8839,
     & -3.1301, -3.2942, -3.1757, -2.8531, -2.9330, -3.0534,
     & -3.1335, -3.1518, -3.1201, -3.0246, -2.8608, -2.9402,
     & -3.0159, -2.9915, -2.9741, -2.8522, -2.7832, -2.9477/
c PAIR SA P
      data (pmf(14,22,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SA SA
      data (pmf(14,23,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SA SD
      data (pmf(14,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SA HL
      data (pmf(14,25,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SA F
      data (pmf(14,26,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD CF
      data (pmf(15, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.7707, -2.6232, -3.6402, -3.8874, -3.6816,
     & -2.9678, -2.4501, -1.9394, -1.8438, -2.1625, -3.0792,
     & -3.5459, -3.5799, -3.1597, -3.0861, -3.0440, -3.0821,
     & -3.2032, -3.3218, -3.2598, -3.0499, -3.0388, -3.3611,
     & -3.3242, -3.4439, -3.0781, -3.1697, -3.2879, -3.0399,
     & -3.0217, -3.0453, -3.1279, -2.8877, -3.0682, -3.0621,
     & -2.9189, -3.0595, -3.1542, -2.9848, -2.9268, -3.0518,
     & -3.0103, -3.0180, -2.9080, -2.9094, -3.0439, -2.9287,
     & -2.8021, -2.9377, -2.8690, -2.7452, -2.7408, -2.8597/
c PAIR SD CP
      data (pmf(15, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.9567, -3.2346, -4.0661, -4.2719,
     & -3.2820, -2.3377, -2.7605, -3.1152, -2.4145, -2.9574,
     & -3.3143, -3.2426, -2.9762, -3.1498, -3.2161, -3.1210,
     & -3.4827, -3.3643, -3.1343, -3.0541, -2.8977, -2.9061,
     & -3.1857, -3.0802, -2.9913, -2.9615, -3.0522, -3.1626,
     & -2.9991, -2.8644, -2.7030, -3.0091, -2.9080, -2.8709,
     & -2.9531, -2.9679, -2.8254, -2.7778, -2.9447, -3.1424,
     & -2.9463, -2.9500, -3.0257, -3.0069, -3.0252, -2.9196,
     & -2.9563, -2.9060, -3.0255, -3.1369, -3.0420, -3.0640/
c PAIR SD cF
      data (pmf(15, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -3.3765, -4.5211, -3.5601,
     & -3.1591, -2.7426, -2.3958, -3.7848, -3.7770, -4.1185,
     & -3.9208, -4.0031, -3.6355, -3.2498, -3.5308, -3.6772,
     & -3.7095, -3.7029, -3.4415, -3.3927, -3.3007, -3.2696,
     & -3.4265, -3.0763, -3.2803, -2.9750, -3.1748, -2.8025,
     & -2.7522, -3.0240, -2.8861, -3.0058, -3.0827, -3.1961,
     & -3.2172, -2.8100, -2.7488, -2.8561, -2.8551, -2.8218,
     & -2.8511, -2.8199, -2.7769, -2.7647, -2.6694, -2.9020,
     & -2.8438, -2.8855, -2.8387, -2.8003, -2.8064, -2.8477/
c PAIR SD cP
      data (pmf(15, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -3.4993, -3.8681, -4.2379, -3.9095,
     & -3.0036,  0.0000, -2.7543, -3.7610, -4.0641, -4.1111,
     & -4.2016, -4.1965, -3.9720, -3.5507, -3.9046, -3.4837,
     & -3.2229, -3.0299, -3.5533, -3.4905, -3.2502, -3.2798,
     & -3.2503, -3.1672, -2.9704, -2.8878, -3.1340, -2.9180,
     & -2.9429, -3.2076, -3.1030, -2.8934, -3.0172, -2.8470,
     & -3.0492, -3.0069, -2.9544, -3.0343, -2.7833, -2.5819,
     & -2.7812, -2.8093, -3.0135, -2.8913, -2.7735, -2.8080,
     & -2.9296, -2.8135, -2.7925, -2.9162, -2.6887, -2.7264/
c PAIR SD C3
      data (pmf(15, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -3.4583,  0.0000, -3.8574, -3.5109,
     & -2.9724, -3.5002, -2.7395, -3.0429, -2.9434, -3.6628,
     & -3.5694, -3.9592, -4.2672, -4.3606, -3.9951, -3.7551,
     & -3.6250, -3.1451, -3.1845, -3.0650, -3.1416, -3.6164,
     & -3.4546, -2.9998, -2.9554, -2.7137, -2.6743, -2.7435,
     & -2.9772, -2.9391, -2.9884, -2.9270, -3.0844, -3.0099,
     & -3.0961, -3.0647, -3.2190, -3.1360, -2.8167, -3.0046,
     & -2.8704, -2.6873, -2.9420, -2.4915, -2.7653, -2.5862,
     & -2.9055, -2.8502, -2.6122, -2.8234, -2.8593, -2.7209/
c PAIR SD CW
      data (pmf(15, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -3.8509, -4.0343, -4.2739,
     & -3.6592,  0.0000,  0.0000, -3.3046, -3.3324, -3.4244,
     & -3.4027, -4.0030, -3.5038, -3.2632, -2.7695, -2.7962,
     & -3.3989, -3.2235, -3.1524, -2.9597, -3.0148, -3.0011,
     & -3.0490, -3.1101, -3.3794, -2.9930, -2.9272, -3.1454,
     & -3.0532, -2.9765, -2.7770, -2.9262, -3.3402, -3.1002,
     & -3.1867, -3.1110, -2.8071, -2.9798, -3.0485, -2.9866,
     & -3.0573, -3.1158, -3.2086, -3.1086, -3.0919, -2.6607,
     & -2.8938, -2.7286, -2.6231, -2.8366, -2.6397, -2.8635/
c PAIR SD CO
      data (pmf(15, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -3.4652,  0.0000, -3.1865, -3.0685, -2.9613,
     &  0.0000,  0.0000, -2.6321, -2.9449, -2.4452, -3.1726,
     & -3.0866, -3.1364, -2.7569,  0.0000, -2.9106, -2.2998,
     & -3.2382, -3.4958, -2.7531, -3.2703, -3.0412, -3.0421,
     & -3.2514, -3.0381, -2.8632, -2.9529, -3.0545, -3.3434,
     & -3.0683, -2.6753, -3.2007, -3.3880, -3.3691, -3.0452,
     & -3.3325, -3.3934, -3.1204, -2.9688, -2.8099, -2.7579,
     & -2.4385, -2.7063, -2.7763, -2.9987, -3.0158, -3.0921,
     & -3.0925, -2.9185, -2.8227, -2.9042, -2.7317, -2.8503/
c PAIR SD CN
      data (pmf(15, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -3.5579, -3.8194,  0.0000, -3.5703, -5.3000,
     & -3.7454,  0.0000, -2.6998, -2.5964, -2.4990, -3.7612,
     & -3.8684, -2.6335, -3.1990, -2.7111, -2.8038, -2.9692,
     & -3.6375, -3.5103, -3.1736, -3.4933, -3.3777, -3.3034,
     & -3.2291, -3.0468, -2.5541, -3.0000, -2.9237, -2.7549,
     & -3.0489, -3.2975, -2.8180, -2.4457, -2.7530, -2.9609,
     & -3.0997, -2.9470, -2.8951, -2.9745, -2.8393, -2.8790,
     & -3.2506, -3.2150, -2.7363, -2.9438, -2.8522, -2.5907,
     & -2.8756, -3.0173, -3.0081, -2.8846, -2.9049, -2.6050/
c PAIR SD NC
      data (pmf(15, 9,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD NP
      data (pmf(15,10,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD NA
      data (pmf(15,11,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD ND
      data (pmf(15,12,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD NR
      data (pmf(15,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -3.3215,  0.0000, -3.0724,  0.0000,
     &  0.0000, -3.5329, -4.3699, -4.3959, -4.3417, -3.7101,
     & -4.0491, -3.8132, -3.9656, -3.6436, -3.1531, -3.0056,
     & -2.1130, -3.0985, -3.1993, -3.2915, -2.7881, -3.2042,
     & -3.3776, -3.5193, -3.2069, -3.2031, -3.0073, -2.9943,
     & -2.6800, -2.8691, -2.7467, -2.9628, -2.8427, -2.7868,
     & -2.9252, -3.0549, -2.9356, -2.8382, -2.7529, -2.8369,
     & -2.9181, -3.0303, -2.9471, -2.9230, -2.9570, -2.8762,
     & -3.0770, -3.0085, -2.7705, -2.7627, -2.8585, -2.8060/
c PAIR SD N0
      data (pmf(15,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD NS
      data (pmf(15,15,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD OC
      data (pmf(15,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -3.6101, -4.0890, -3.3849,
     &  0.0000, -2.3534,  0.0000, -3.5255, -2.0819, -2.4064,
     & -2.7351, -3.4766, -3.0642, -2.3461, -3.0950, -2.3805,
     & -2.7952, -2.7959, -3.2280, -3.1697, -3.1134, -3.3315,
     & -3.2353, -2.8963, -2.8531, -2.8111, -3.0861, -3.1300,
     & -3.0300, -2.9942, -3.0808, -2.9596, -2.9590, -3.2646,
     & -3.4647, -3.1957, -3.1407, -3.2128, -3.1209, -2.9915,
     & -2.9948, -2.8283, -2.7786, -2.8380, -3.0190, -2.9389,
     & -2.8232, -2.9869, -2.9804, -2.9377, -2.7839, -2.6619/
c PAIR SD OA
      data (pmf(15,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -3.6413, -3.5233, -2.7678,
     & -3.3041,  0.0000,  0.0000, -3.4205, -2.9254, -3.4824,
     & -3.6972, -3.3755, -3.4954, -3.8699, -3.7551, -3.5417,
     & -3.3213, -3.2893, -3.2881, -3.5211, -3.4045, -3.2921,
     & -3.3056, -2.9887, -2.7614, -3.2033, -3.1066, -3.1629,
     & -3.4414, -3.3289, -3.5207, -3.2361, -2.9702, -2.7976,
     & -3.2788, -2.4846, -2.5513, -2.7213, -2.6496, -2.9456,
     & -2.9184, -2.8234, -2.8807, -3.0867, -2.8681, -2.6959,
     & -2.7331, -2.7857, -2.7647, -2.8024, -2.7708, -2.8153/
c PAIR SD OE
      data (pmf(15,18,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD OR
      data (pmf(15,19,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD OS
      data (pmf(15,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -3.3022, -2.5065, -3.0236, -3.9290, -3.7072,
     & -2.8568, -1.9355,  0.0000, -3.2142, -2.9588, -3.4007,
     & -3.0935, -2.8844, -3.1944, -3.1729, -2.8021, -2.8748,
     & -3.0802, -2.9562, -2.8019, -2.9190, -2.8425, -2.8337,
     & -2.7095, -3.1439, -3.2973, -2.9444, -2.8947, -2.9534,
     & -3.0488, -2.9195, -3.2700, -3.2275, -3.0008, -3.1690,
     & -3.1845, -3.0780, -2.9502, -3.1158, -3.1182, -3.1035,
     & -3.0396, -3.0103, -2.9545, -2.7561, -2.9416, -2.8792,
     & -2.9103, -2.9134, -3.0710, -2.9845, -2.9195, -2.8353/
c PAIR SD OD
      data (pmf(15,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -3.2430, -3.2948, -2.3696,
     & -3.0701, -2.5518,  0.0000, -3.2409, -3.1503, -3.3188,
     & -3.3459, -3.1888, -3.2791, -2.8505, -2.7194, -3.0855,
     & -2.8919, -3.5075, -3.6530, -2.9920, -2.8865, -2.9600,
     & -3.2531, -3.0932, -2.5128, -2.4361, -2.8527, -2.9839,
     & -3.1624, -3.3457, -2.9764, -2.8498, -2.7702, -2.6874,
     & -2.6284, -2.7528, -2.9490, -2.9887, -3.0334, -3.0757,
     & -2.9563, -2.8353, -2.8508, -2.9166, -2.9396, -3.0418,
     & -3.0403, -2.9434, -3.0594, -3.2331, -3.1442, -3.0091/
c PAIR SD P
      data (pmf(15,22,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD SA
      data (pmf(15,23,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD SD
      data (pmf(15,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD HL
      data (pmf(15,25,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR SD F
      data (pmf(15,26,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR HH CF
      data (pmf(16, 1,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.2709,
     &  0.0000,  0.0000,  0.0000, -1.4450, -1.3270, -1.2197,
     & -1.7982, -2.2209, -2.5094, -2.7132, -2.8717, -2.9949,
     & -3.0394, -3.0558, -3.1098, -3.1398, -3.1300, -3.1769,
     & -3.1154, -3.1064, -3.1124, -3.0063, -3.0338, -3.0055,
     & -3.0061, -2.9811, -2.9885, -3.0734, -3.0267, -3.0479,
     & -2.9978, -2.9874, -3.0249, -2.9907, -3.0224, -3.0500,
     & -3.0152, -3.0757, -3.0482, -3.0285, -3.0344, -3.0419,
     & -3.0184, -2.9705, -3.0209, -3.0181, -3.0130, -3.0045,
     & -2.9723, -2.9878, -2.9687, -2.9924, -2.9492, -2.9543/
c PAIR HH CP
      data (pmf(16, 2,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -1.3453, -1.2273, -1.8866,
     & -2.4161, -2.9825, -3.1932, -3.1790, -3.1822, -3.2195,
     & -3.3162, -3.3427, -3.3982, -3.3377, -3.2544, -3.2048,
     & -3.2085, -3.1740, -3.1672, -3.1713, -3.1510, -3.1100,
     & -3.1432, -3.1232, -3.1035, -3.0396, -3.0558, -3.0258,
     & -3.0568, -3.0295, -3.0079, -3.0230, -3.0155, -3.0380,
     & -3.0403, -3.0227, -3.0202, -3.0196, -2.9875, -2.9651,
     & -2.9953, -2.9748, -2.9682, -2.9730, -2.9789, -2.9661,
     & -2.9652, -2.9654, -2.9294, -2.9076, -2.9208, -2.9245/
c PAIR HH cF
      data (pmf(16, 3,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.3741,  0.0000, -1.8760,  0.0000, -1.7585, -2.1164,
     & -2.6040, -2.6393, -2.7996, -3.0219, -3.1030, -3.1087,
     & -3.0398, -3.2393, -3.2531, -3.2686, -3.2434, -3.2037,
     & -3.1789, -3.1802, -3.0980, -3.1535, -3.0515, -3.0911,
     & -3.0643, -3.0663, -3.0602, -3.0339, -3.0165, -3.0022,
     & -3.0160, -3.0454, -3.0220, -3.0520, -3.0621, -3.0543,
     & -3.0590, -3.0272, -3.0136, -2.9979, -2.9973, -2.9883,
     & -2.9993, -2.9913, -3.0272, -3.0068, -2.9834, -2.9873,
     & -2.9524, -2.9427, -2.9446, -2.9510, -2.9422, -2.9421/
c PAIR HH cP
      data (pmf(16, 4,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.9032,
     &  0.0000, -1.5379,  0.0000, -2.4074, -2.1985, -2.4488,
     & -2.4081, -2.8799, -3.1520, -3.3323, -3.3123, -3.2676,
     & -3.1539, -3.2820, -3.3136, -3.1775, -3.1893, -3.2074,
     & -3.1945, -3.1940, -3.2163, -3.1711, -3.1438, -3.1145,
     & -3.1553, -3.1019, -3.0543, -3.0862, -3.1069, -3.1059,
     & -3.0872, -3.0543, -3.0783, -3.0691, -3.0373, -3.0398,
     & -3.0537, -3.0484, -2.9895, -2.9977, -2.9907, -3.0021,
     & -3.0066, -2.9935, -3.0091, -2.9544, -2.9339, -2.9681,
     & -2.9145, -2.8984, -2.9279, -2.9141, -2.9112, -2.8978/
c PAIR HH C3
      data (pmf(16, 5,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.9633,  0.0000, -1.6472, -1.5161, -1.8071, -2.5177,
     & -1.9792, -2.7462, -3.1475, -3.1426, -2.9425, -3.1398,
     & -3.0464, -3.1716, -3.1317, -3.0404, -3.1457, -2.9546,
     & -3.1102, -3.1132, -3.1531, -3.1085, -3.1166, -2.9558,
     & -2.9999, -3.1011, -3.0451, -3.0746, -3.0905, -3.0862,
     & -3.0372, -2.9859, -3.0166, -3.0339, -2.9964, -2.9813,
     & -3.0012, -3.0781, -3.0288, -3.0181, -2.9915, -2.9141,
     & -2.9816, -3.0881, -3.0339, -2.9839, -2.9713, -2.9561,
     & -2.9518, -3.0315, -2.9797, -3.0070, -2.9689, -2.9727/
c PAIR HH CW
      data (pmf(16, 6,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.0787, -1.5222, -2.0393, -1.6821, -2.4622,
     & -2.6655, -2.8478, -3.3087, -3.5285, -3.3604, -3.1637,
     & -3.2881, -3.2707, -3.2261, -3.2160, -3.1878, -3.1707,
     & -3.0687, -3.0193, -3.1370, -3.1828, -3.1037, -3.0241,
     & -3.0912, -3.1120, -3.1134, -3.0291, -3.0709, -3.0627,
     & -3.0882, -2.9986, -2.9187, -2.9602, -2.9977, -3.0032,
     & -3.0229, -3.0342, -2.9934, -3.0266, -3.0666, -2.9840,
     & -3.0076, -2.9803, -3.0346, -3.0164, -3.0035, -3.0044,
     & -2.9615, -2.9642, -2.9761, -2.9437, -2.9342, -2.8865/
c PAIR HH CO
      data (pmf(16, 7,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.2352, -1.7083, -3.1580,
     & -3.4495, -3.5793, -3.6253, -3.5611, -3.3309, -3.3798,
     & -3.1275, -3.2015, -3.2204, -3.3287, -3.3119, -3.2111,
     & -3.1977, -3.0725, -3.0036, -3.0087, -3.0588, -3.0582,
     & -3.0392, -3.2207, -3.2028, -3.1635, -3.0830, -3.0097,
     & -3.0665, -2.9866, -3.0072, -2.9954, -3.0216, -2.9694,
     & -2.8541, -2.9546, -3.0008, -3.0143, -3.0257, -3.0122,
     & -2.9655, -2.9455, -3.0225, -2.9251, -2.9275, -2.9558,
     & -2.9481, -2.9881, -2.9177, -2.9767, -2.9395, -2.9304/
c PAIR HH CN
      data (pmf(16, 8,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.6910, -2.7991, -2.9825,  0.0000,
     & -2.2081, -2.7386, -2.9299, -3.0251, -3.1669, -3.4204,
     & -3.4021, -3.3802, -3.2045, -3.3302, -3.2732, -3.1894,
     & -3.2110, -3.2227, -3.1105, -3.1828, -3.1460, -3.1246,
     & -2.9029, -3.1253, -3.0616, -3.0701, -3.2381, -3.0044,
     & -3.1223, -3.1392, -2.9694, -3.0370, -3.0293, -3.1038,
     & -2.9871, -3.0650, -2.9936, -2.9444, -2.9887, -2.9165,
     & -2.9915, -3.0094, -3.0053, -2.9331, -2.9959, -3.0019,
     & -2.9155, -2.9400, -2.9189, -2.9003, -2.9133, -2.8639/
c PAIR HH NC
      data (pmf(16, 9,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.9020, -2.3751, -3.2736,
     & -2.8034, -3.0002, -3.1009, -3.1891, -3.4628, -3.3030,
     & -3.2738, -3.3326, -3.2962, -3.2606, -3.2690, -2.9569,
     & -3.2059, -3.3199, -3.3126, -3.1785, -3.0912, -3.0210,
     & -2.9431, -3.0950, -3.1793, -3.1022, -3.0442, -3.0132,
     & -3.1643, -3.1215, -3.0473, -3.1197, -3.0663, -3.0054,
     & -2.9109, -3.0179, -2.9657, -2.9238, -2.9600, -2.9424,
     & -2.9522, -2.9855, -2.9554, -3.0204, -3.0975, -2.9964,
     & -2.8836, -2.9075, -2.9322, -2.9245, -2.8454, -2.8927/
c PAIR HH NP
      data (pmf(16,10,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.5350, -2.4170, -2.7187,
     & -2.5379, -2.6140, -2.4653, -2.6299, -3.0198, -3.1230,
     & -3.3183, -3.3523, -3.2997, -3.0279, -3.1122, -3.2268,
     & -3.3264, -3.0728, -3.2381, -3.2576, -3.1415, -3.1054,
     & -3.1599, -3.1809, -3.0923, -3.0616, -3.1593, -3.1375,
     & -3.0313, -3.1265, -3.0910, -2.9569, -3.0943, -3.0601,
     & -3.0769, -3.0486, -3.0449, -2.9237, -3.0463, -3.0133,
     & -2.9382, -2.9717, -2.9894, -2.9652, -2.9828, -2.9249,
     & -2.9386, -2.9330, -2.9317, -2.9238, -2.8825, -2.8786/
c PAIR HH NA
      data (pmf(16,11,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR HH ND
      data (pmf(16,12,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.9369,  0.0000, -1.6583, -2.8988, -2.2509,
     & -2.5154, -2.7120, -2.8271, -3.0452, -2.9520, -3.4036,
     & -3.1333, -3.3003, -3.3996, -3.2077, -3.1483, -3.1783,
     & -3.2481, -3.1686, -3.0710, -3.0244, -3.0631, -3.1169,
     & -3.0566, -3.0548, -2.9739, -3.0450, -3.0826, -3.1192,
     & -3.0214, -3.0328, -3.0477, -3.0916, -3.0032, -3.0339,
     & -3.0480, -2.9802, -3.0338, -2.9730, -3.0067, -3.0382,
     & -2.9877, -3.0292, -2.9751, -3.0298, -3.0092, -2.9579,
     & -2.9752, -2.9611, -2.9802, -3.0053, -2.9247, -2.9171/
c PAIR HH NR
      data (pmf(16,13,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.7258, -3.1353, -3.1248, -3.1492,
     & -2.9789, -3.0976, -3.0125, -2.9581, -3.0812, -2.9794,
     & -3.1279, -3.0796, -3.2816, -3.1174, -3.2374, -3.0200,
     & -3.1531, -3.1343, -3.1834, -3.1957, -3.0028, -3.0928,
     & -3.0768, -3.1065, -3.1170, -3.0521, -3.0357, -3.1395,
     & -3.0493, -3.0862, -3.0570, -3.0353, -3.0737, -3.0730,
     & -3.0319, -3.0636, -3.0250, -3.0476, -3.0250, -3.0376,
     & -3.0665, -2.9684, -2.9527, -2.8690, -2.9641, -2.9908,
     & -2.9540, -2.9052, -2.9395, -2.9043, -2.9265, -2.9217/
c PAIR HH N0
      data (pmf(16,14,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR HH NS
      data (pmf(16,15,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000, -2.7143, -3.7444, -3.3069,
     & -3.0160, -3.3122, -2.9675, -3.1985, -3.1838, -3.0121,
     & -3.5433, -3.3635, -3.4870, -3.5422, -3.3354, -3.0599,
     & -2.9168, -2.9240, -3.0633, -2.9215, -3.2105, -3.0047,
     & -2.9798, -3.2010, -3.1708, -3.0362, -2.8881, -3.1244,
     & -3.0069, -3.1983, -3.0084, -2.9421, -3.1656, -3.0442,
     & -3.0571, -3.0021, -2.9477, -3.0204, -3.0380, -2.8676,
     & -3.0720, -2.7499, -2.9511, -3.0575, -3.0264, -3.0262,
     & -2.9823, -2.7921, -2.9501, -2.9182, -2.7869, -2.9712/
c PAIR HH OC
      data (pmf(16,16,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.6904, -2.9014, -3.8562, -3.7657, -3.5063,
     & -3.3347, -3.3479, -3.2725, -3.3260, -3.3431, -3.3399,
     & -3.2992, -3.2850, -3.3854, -3.1170, -3.0829, -3.1424,
     & -3.0382, -3.1680, -3.1798, -3.0594, -3.1383, -3.1355,
     & -3.0906, -3.1443, -3.0704, -3.0953, -3.1236, -3.0722,
     & -3.1143, -3.0377, -2.9750, -2.9673, -2.9822, -2.9970,
     & -2.9969, -2.9156, -3.0212, -3.0103, -2.9758, -3.0386,
     & -2.9796, -2.9789, -2.9714, -2.9844, -2.9081, -2.9564,
     & -2.9553, -2.9355, -2.9106, -2.8851, -2.9449, -2.9081/
c PAIR HH OA
      data (pmf(16,17,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     & -1.7063, -1.9467, -2.7487, -3.8445, -3.7623, -3.5319,
     & -3.0681, -3.0768, -3.0632, -2.9718, -3.1432, -3.0300,
     & -3.2188, -3.1076, -2.9502, -2.9439, -2.9994, -2.9820,
     & -2.9957, -3.1146, -3.0251, -3.1200, -2.9239, -3.0474,
     & -3.0909, -3.0540, -3.1635, -3.0911, -2.9771, -2.9462,
     & -2.9016, -3.0273, -3.0608, -3.0744, -3.1216, -3.0488,
     & -3.0633, -3.0417, -3.0553, -2.9982, -2.9475, -3.0178,
     & -3.0395, -3.0590, -2.9877, -2.9251, -2.9235, -2.9470,
     & -2.9977, -2.9379, -2.9674, -2.9777, -2.9604, -2.9643/
c PAIR HH OE
      data (pmf(16,18,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -2.2617, -3.0637, -3.4964, -3.5009, -3.2239,
     & -3.1204, -3.2054, -3.1632, -3.4654, -3.2375, -3.3102,
     & -3.1400, -3.1532, -3.3023, -3.1745, -3.3598, -3.3854,
     & -3.3772, -3.3209, -3.3951, -3.2704, -3.1734, -3.1572,
     & -3.2359, -3.1401, -3.1254, -3.1464, -3.1132, -3.1107,
     & -3.1955, -3.1051, -3.0580, -3.0744, -3.0534, -3.0638,
     & -2.9999, -2.9193, -2.9917, -2.9835, -2.9408, -2.9087,
     & -2.8865, -2.9675, -2.9519, -2.9937, -2.9737, -2.9022,
     & -2.8979, -2.9062, -2.8626, -2.9282, -2.9182, -2.8879/
c PAIR HH OR
      data (pmf(16,19,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR HH OS
      data (pmf(16,20,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.5706,
     & -1.7829, -2.7187, -3.7321, -3.9638, -3.8213, -3.4417,
     & -3.4099, -3.5170, -3.4573, -3.5306, -3.6371, -3.6315,
     & -3.5313, -3.5061, -3.5396, -3.4750, -3.3807, -3.2400,
     & -3.1801, -3.2355, -3.1434, -3.1363, -3.1564, -3.1204,
     & -3.0419, -3.0472, -3.0082, -3.0263, -3.0523, -2.9782,
     & -2.9997, -3.0612, -2.9949, -3.0302, -2.9923, -2.9621,
     & -2.9689, -2.9935, -2.9771, -2.9768, -2.9226, -2.9449,
     & -2.9172, -2.9330, -2.9161, -2.9685, -2.9139, -2.9104,
     & -2.9291, -2.9379, -2.9635, -2.9360, -2.8981, -2.8917/
c PAIR HH OD
      data (pmf(16,21,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000, -1.6090, -2.7241, -3.4195, -3.3962, -3.2356,
     & -3.2717, -3.1746, -3.3156, -3.2618, -3.2573, -3.3219,
     & -3.2245, -3.2096, -3.2790, -3.2512, -3.2254, -3.2029,
     & -3.1694, -3.1849, -3.1688, -3.1640, -3.1500, -3.1244,
     & -3.0736, -3.0413, -3.1169, -3.0715, -3.0827, -3.1030,
     & -3.0857, -3.0706, -3.0633, -3.0652, -3.0368, -2.9809,
     & -3.0117, -3.0041, -2.9208, -2.9739, -2.9645, -2.9888,
     & -2.9960, -3.0054, -2.9709, -2.9775, -2.9741, -2.9683,
     & -2.9528, -2.9201, -2.9269, -2.9025, -2.9242, -2.9263/
c PAIR HH P
      data (pmf(16,22,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000, -2.0488, -1.9177, -1.7997, -2.5103,
     & -3.5553, -3.9766, -4.1896, -4.2313, -3.9484, -3.7227,
     & -3.6205, -3.5613, -3.4929, -3.5323, -3.4263, -3.5080,
     & -3.1873, -3.0820, -3.1442, -3.0779, -3.0854, -3.1343,
     & -3.1637, -2.9912, -3.0293, -3.0686, -2.9942, -2.9482,
     & -3.0619, -2.9839, -3.0381, -3.1412, -2.9736, -3.0203,
     & -2.9412, -2.9657, -2.9133, -2.8745, -2.8193, -2.9282,
     & -2.8665, -2.9709, -2.9495, -2.9203, -2.9221, -2.9423,
     & -3.0248, -2.9478, -3.0032, -2.9166, -2.8976, -2.8737/
c PAIR HH SA
      data (pmf(16,23,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR HH SD
      data (pmf(16,24,i),i=1,60) /
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000,
     & -3.0000, -3.0000, -3.0000, -3.0000, -3.0000, -3.0000/
c PAIR HH HL
      data (pmf(16,25,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000, -0.7589, -1.4039,
     & -1.8586, -2.0653, -2.2329, -2.3610, -2.6839, -2.9824,
     & -3.0767, -3.1760, -3.1906, -3.1557, -3.1802, -3.1521,
     & -3.1790, -3.1571, -3.1899, -3.1777, -3.1834, -3.1791,
     & -3.1757, -3.1583, -3.1348, -3.1277, -3.1169, -3.1056,
     & -3.0883, -3.0871, -3.0775, -3.0704, -3.0704, -3.0534,
     & -3.0584, -3.0781, -3.0737, -3.0542, -3.0350, -3.0491,
     & -3.0594, -3.0267, -3.0106, -3.0088, -3.0201, -3.0073,
     & -2.9916, -2.9891, -2.9832, -2.9883, -2.9654, -2.9809,
     & -2.9544, -2.9503, -2.9466, -2.9354, -2.9296, -2.9408/
c PAIR HH F
      data (pmf(16,26,i),i=1,60) /
     &  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
     &  0.0000,  0.0000,  0.0000,  0.0000, -3.0234, -3.3251,
     & -3.2091, -3.1013, -2.5914,  0.0000, -2.8160, -3.2678,
     & -2.8821, -3.2108, -3.1341, -2.9532, -2.9886, -3.1586,
     & -2.9437, -3.0279, -2.8177, -2.9668, -2.8453, -3.0592,
     & -3.0138, -3.1180, -2.7841, -3.0326, -2.6475, -3.2384,
     & -3.1377, -2.9469, -2.9405, -3.0558, -2.8754, -2.8714,
     & -2.9390, -2.8364, -3.0031, -2.9174, -3.0947, -2.9747,
     & -3.0420, -3.0169, -2.9156, -3.1522, -3.0604, -2.9627,
     & -2.9409, -3.0902, -3.1171, -2.9624, -3.0077, -3.0277/

c Pair H H
      data (dfire(1,1,i),i=1,20) /
     &  7494.35, 0.279, -0.010, 0.013, 0.007, 0.003,
     &  0.002, 0.001, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair H C3
      data (dfire(1,2,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair H Car
      data (dfire(1,3,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair H C2
      data (dfire(1,4,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair H Ccat
      data (dfire(1,5,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair H O2
      data (dfire(1,6,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair H O3
      data (dfire(1,7,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair H Oco2
      data (dfire(1,8,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair H Npl3
      data (dfire(1,9,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair H Nam
      data (dfire(1,10,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair H N4
      data (dfire(1,11,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair H N2
      data (dfire(1,12,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair H S3
      data (dfire(1,13,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair C3 H
      data (dfire(2,1,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair C3 C3                                                              
      data (dfire(2,2,i),i=1,20) /
     &  6.570, 5.587, 2.007, 1.984, 0.110, 0.079,
     & -0.276,-0.253,-0.275,-0.238,-0.187,-0.107,
     & -0.049,-0.081,-0.128,-0.121,-0.076,-0.032,
     & -0.011, 0.000/
c Pair C3 Car                                                             
      data (dfire(2,3,i),i=1,20) /
     &  6.437, 7.516, 4.764, 1.591, 0.183,-0.056,
     & -0.151,-0.178,-0.156,-0.165,-0.150,-0.135,
     & -0.141,-0.172,-0.172,-0.143,-0.104,-0.062,
     & -0.031, 0.000/
c Pair C3 C2                                                              
      data (dfire(2,4,i),i=1,20) /
     &  6.889, 2.391, 1.895, 0.955, 0.165,-0.181,
     & -0.316,-0.221,-0.184,-0.046,-0.127,-0.138,
     & -0.080,-0.078,-0.096,-0.088,-0.060,-0.030,
     & -0.014, 0.000/
c Pair C3 Ccat                                                            
      data (dfire(2,5,i),i=1,20) /
     &  5.919, 4.792, 4.985, 1.769, 0.523, 0.246,
     &  0.107, 0.023, 0.027, 0.058, 0.090, 0.015,
     & -0.035,-0.075,-0.108,-0.074,-0.027,-0.010,
     & -0.002, 0.000/
c Pair C3 O2                                                              
      data (dfire(2,6,i),i=1,20) /
     &  6.358, 4.100,-0.162, 0.205,-0.283,-0.260,
     & -0.093,-0.077, 0.036,-0.119,-0.157,-0.087,
     & -0.190,-0.167,-0.067,-0.059,-0.083,-0.067,
     & -0.030, 0.000/
c Pair C3 O3                                                              
      data (dfire(2,7,i),i=1,20) /
     &  5.619, 5.182, 2.284, 0.578, 0.034,-0.035,
     & -0.151,-0.132,-0.119,-0.124,-0.112,-0.095,
     & -0.095,-0.087,-0.075,-0.090,-0.070,-0.034,
     & -0.011, 0.000/
c Pair C3 Oco2                                                            
      data (dfire(2,8,i),i=1,20) /
     &  5.762, 4.366, 2.192, 0.552, 0.093, 0.126,
     &  0.063, 0.017, 0.010,-0.038,-0.086,-0.084,
     & -0.082,-0.088,-0.061,-0.039,-0.042,-0.033,
     & -0.020, 0.000/
c Pair C3 Np13                                                            
      data (dfire(2,9,i),i=1,20) /
     &  6.735, 5.608, 3.027, 1.250, 0.309, 0.165,
     &  0.057, 0.019, 0.027, 0.031,-0.009,-0.035,
     & -0.054,-0.099,-0.102,-0.078,-0.053,-0.028,
     & -0.015, 0.000/
c Pair C3 Nam                                                             
      data (dfire(2,10,i),i=1,20) /
     &  6.376, 5.336, 1.767, 0.155, 0.523,-0.057,
     & -0.251,-0.277,-0.272,-0.143,-0.095,-0.104,
     & -0.029,-0.053,-0.105,-0.096,-0.056,-0.018,
     & -0.010, 0.000/
c Pair C3 N4                                                              
      data (dfire(2,11,i),i=1,20) /
     &  4.243, 4.219, 2.199, 1.144, 0.820, 0.474,
     &  0.242, 0.140, 0.109, 0.098, 0.077, 0.048,
     &  0.021,-0.045,-0.090,-0.061,-0.038,-0.004,
     & -0.005, 0.000/
c Pair C3 N2                                                              
      data (dfire(2,12,i),i=1,20) /
     &  5.574, 5.550, 3.331, 1.210, 0.233, 0.125,
     & -0.003,-0.066,-0.033, 0.013, 0.029, 0.004,
     & -0.047,-0.142,-0.115,-0.074,-0.040,-0.030,
     & -0.020, 0.000/
c Pair C3 S3                                                              
      data (dfire(2,13,i),i=1,20) /
     &  4.378, 4.766, 1.896, 0.773, 0.077,-0.124,
     & -0.170,-0.230,-0.186,-0.220,-0.187,-0.175,
     & -0.169,-0.197,-0.175,-0.137,-0.105,-0.061,
     & -0.028, 0.000/    
c Pair Car H                                                             
      data (dfire(3,1,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair Car C3                                                             
      data (dfire(3,2,i),i=1,20) /
     &  6.437, 7.516, 4.764, 1.591, 0.183,-0.056,
     & -0.151,-0.178,-0.156,-0.165,-0.150,-0.135,
     & -0.141,-0.172,-0.172,-0.143,-0.104,-0.062,
     & -0.031, 0.000/
c Pair Car Car                                                            
      data (dfire(3,3,i),i=1,20) /
     &  6.633, 6.609, 4.873, 1.359,-0.019,-0.172,
     & -0.258,-0.270,-0.279,-0.325,-0.239,-0.174,
     & -0.139,-0.135,-0.140,-0.140,-0.122,-0.093,
     & -0.054, 0.000/
c Pair Car C2                                                             
      data (dfire(3,4,i),i=1,20) /
     &  7.029, 7.005, 4.722, 1.424, 0.359, 0.067,
     & -0.025,-0.077,-0.125,-0.147,-0.207,-0.215,
     & -0.206,-0.193,-0.166,-0.122,-0.091,-0.047,
     & -0.023, 0.000/
c Pair Car Ccat                                                           
      data (dfire(3,5,i),i=1,20) /
     &  5.006, 3.879, 4.072, 1.088, 0.092, 0.030,
     &  0.015, 0.076, 0.094, 0.132, 0.176, 0.183,
     &  0.195, 0.144, 0.062, 0.008,-0.021,-0.038,
     & -0.028, 0.000/
c Pair Car O2                                                             
      data (dfire(3,6,i),i=1,20) /
     &  6.918, 5.137, 3.075, 0.751, 0.189, 0.054,
     & -0.037,-0.034,-0.122,-0.149,-0.183,-0.223,
     & -0.235,-0.227,-0.160,-0.126,-0.088,-0.058,
     & -0.025, 0.000/
c Pair Car O3                                                             
      data (dfire(3,7,i),i=1,20) /
     &  5.804, 5.780, 3.084, 0.810, 0.305, 0.217,
     &  0.072, 0.058,-0.028,-0.062,-0.051,-0.054,
     & -0.066,-0.097,-0.109,-0.110,-0.097,-0.068,
     & -0.032, 0.000/
c Pair Car Oco2                                                           
      data (dfire(3,8,i),i=1,20) /
     &  5.946, 4.406, 2.860, 0.842, 0.451, 0.392,
     &  0.216, 0.208, 0.095, 0.120, 0.145, 0.127,
     &  0.100, 0.050,-0.015,-0.059,-0.069,-0.045,
     & -0.025, 0.000/
c Pair Car Np13                                                           
      data (dfire(3,9,i),i=1,20) /
     &  5.825, 4.697, 3.048, 0.862, 0.056, 0.047,
     &  0.014, 0.060, 0.037, 0.068, 0.104, 0.097,
     &  0.094, 0.056, 0.014,-0.030,-0.047,-0.050,
     & -0.033, 0.000/
c Pair Car Nam                                                            
      data (dfire(3,10,i),i=1,20) /
     &  6.922, 6.898, 3.464, 1.191, 0.383, 0.282,
     &  0.032,-0.064,-0.178,-0.237,-0.267,-0.259,
     & -0.220,-0.203,-0.174,-0.130,-0.081,-0.046,
     & -0.021, 0.000/
c Pair Car N4                                                             
      data (dfire(3,11,i),i=1,20) /
     &  5.102, 3.975, 2.858, 1.304, 0.635, 0.365,
     &  0.205, 0.160, 0.179, 0.180, 0.203, 0.222,
     &  0.257, 0.236, 0.157, 0.084, 0.015,-0.022,
     & -0.027, 0.000/
c Pair Car N2                                                             
      data (dfire(3,12,i),i=1,20) /
     &  4.675, 4.651, 3.328, 0.974, 0.068,-0.027,
     & -0.096, 0.009, 0.038, 0.083, 0.081, 0.104,
     &  0.056,-0.026,-0.059,-0.062,-0.048,-0.052,
     & -0.030, 0.000/
c Pair Car S3                                                             
      data (dfire(3,13,i),i=1,20) /
     &  3.415, 4.907, 3.037, 1.148,-0.118,-0.284,
     & -0.353,-0.397,-0.292,-0.265,-0.182,-0.117,
     & -0.117,-0.162,-0.202,-0.199,-0.170,-0.117,
     & -0.061, 0.000/
c Pair C2 H                                                              
      data (dfire(4,1,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair C2 C3                                                              
      data (dfire(4,2,i),i=1,20) /
     &  6.889, 2.391, 1.895, 0.955, 0.165,-0.181,
     & -0.316,-0.221,-0.184,-0.046,-0.127,-0.138,
     & -0.080,-0.078,-0.096,-0.088,-0.060,-0.030,
     & -0.014, 0.000/
c Pair C2 Car                                                             
      data (dfire(4,3,i),i=1,20) /
     &  7.029, 7.005, 4.722, 1.424, 0.359, 0.067,
     & -0.025,-0.077,-0.125,-0.147,-0.207,-0.215,
     & -0.206,-0.193,-0.166,-0.122,-0.091,-0.047,
     & -0.023, 0.000/
c Pair C2 C2                                                              
      data (dfire(4,4,i),i=1,20) /
     &  7.479, 5.939, 0.111,-0.572, 0.087, 0.131,
     & -0.476,-0.108,-0.074,-0.144,-0.061, 0.030,
     &  0.004,-0.048,-0.100,-0.086,-0.055,-0.044,
     & -0.029, 0.000/
c Pair C2 Ccat                                                            
      data (dfire(4,5,i),i=1,20) /
     &  5.369, 5.345, 4.021, 2.007, 0.582,-0.056,
     & -0.154, 0.024, 0.081,-0.008,-0.143,-0.215,
     & -0.185,-0.168,-0.148,-0.087,-0.040,-0.001,
     &  0.001, 0.000/
c Pair C2 O2                                                              
      data (dfire(4,6,i),i=1,20) /
     &  7.359, 4.005, 0.044,-0.190,-0.616,-0.386,
     &  0.027, 0.085,-0.098,-0.065,-0.039,-0.164,
     & -0.069,-0.074,-0.068,-0.088,-0.081,-0.077,
     & -0.045, 0.000/
c Pair C2 O3                                                              
      data (dfire(4,7,i),i=1,20) /
     &  5.107, 5.083, 2.172, 0.377,-0.216,-0.334,
     & -0.113,-0.076,-0.133,-0.066,-0.126,-0.094,
     & -0.120,-0.098,-0.065,-0.046,-0.043,-0.029,
     & -0.012, 0.000/
c Pair C2 Oco2                                                            
      data (dfire(4,8,i),i=1,20) /
     &  6.348, 5.220, 2.227, 0.879, 0.088, 0.043,
     &  0.049,-0.077,-0.092,-0.086,-0.085,-0.135,
     & -0.138,-0.079,-0.064,-0.027,-0.007,-0.012,
     & -0.007, 0.000/
c Pair C2 Np13                                                            
      data (dfire(4,9,i),i=1,20) /
     &  6.198, 6.174, 2.645, 0.867, 0.016, 0.105,
     &  0.094,-0.028,-0.095,-0.095,-0.123,-0.158,
     & -0.185,-0.169,-0.130,-0.088,-0.049,-0.019,
     & -0.010, 0.000/
c Pair C2 Nam                                                             
      data (dfire(4,10,i),i=1,20) /
     &  6.283, 5.191, 1.805,-0.123,-0.168,-0.091,
     & -0.099,-0.065,-0.290,-0.109, 0.011, 0.044,
     &  0.005,-0.029,-0.053,-0.099,-0.054,-0.035,
     & -0.015, 0.000/
c Pair C2 N4                                                              
      data (dfire(4,11,i),i=1,20) /
     &  4.346, 4.322, 2.086, 0.645, 0.035, 0.199,
     &  0.260, 0.144, 0.039, 0.016,-0.020,-0.079,
     & -0.165,-0.162,-0.137,-0.086,-0.042,-0.006,
     &  0.005, 0.000/
c Pair C2 N2                                                              
      data (dfire(4,12,i),i=1,20) /
     &  5.049, 5.025, 3.701, 0.973, 0.210, 0.072,
     & -0.133,-0.143,-0.128,-0.022,-0.160,-0.148,
     & -0.159,-0.169,-0.107,-0.070,-0.049,-0.022,
     & -0.022, 0.000/
c Pair C2 S3                                                              
      data (dfire(4,13,i),i=1,20) /
     &  4.280, 5.359, 3.490, 1.252, 0.233,-0.145,
     & -0.190,-0.219,-0.214,-0.171,-0.184,-0.297,
     & -0.236,-0.198,-0.156,-0.105,-0.075,-0.052,
     & -0.030, 0.000/
c Pair Ccat H                                                            
      data (dfire(5,1,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair Ccat C3                                                            
      data (dfire(5,2,i),i=1,20) /
     &  5.919, 4.792, 4.985, 1.769, 0.523, 0.246,
     &  0.107, 0.023, 0.027, 0.058, 0.090, 0.015,
     & -0.035,-0.075,-0.108,-0.074,-0.027,-0.010,
     & -0.002, 0.000/
c Pair Ccat Car                                                           
      data (dfire(5,3,i),i=1,20) /
     &  5.006, 3.879, 4.072, 1.088, 0.092, 0.030,
     &  0.015, 0.076, 0.094, 0.132, 0.176, 0.183,
     &  0.195, 0.144, 0.062, 0.008,-0.021,-0.038,
     & -0.028, 0.000/
c Pair Ccat C2                                                            
      data (dfire(5,4,i),i=1,20) /
     &  5.369, 5.345, 4.021, 2.007, 0.582,-0.056,
     & -0.154, 0.024, 0.081,-0.008,-0.143,-0.215,
     & -0.185,-0.168,-0.148,-0.087,-0.040,-0.001,
     &  0.001, 0.000/
c Pair Ccat Ccat                                                          
      data (dfire(5,5,i),i=1,20) /
     &  3.384, 3.360, 3.552, 0.888, 0.003, 0.173,
     &  0.029, 0.010, 0.025, 0.029, 0.136, 0.053,
     & -0.008, 0.004, 0.013,-0.070,-0.053,-0.076,
     & -0.067, 0.000/
c Pair Ccat O2                                                            
      data (dfire(5,6,i),i=1,20) /
     &  5.248, 5.224, 2.372, 0.414,-0.059, 0.091,
     &  0.233, 0.124, 0.064, 0.012,-0.085,-0.202,
     & -0.279,-0.205,-0.110,-0.082,-0.069,-0.034,
     &  0.003, 0.000/
c Pair Ccat O3                                                            
      data (dfire(5,7,i),i=1,20) /
     &  4.115, 4.091, 1.941, 0.508,-0.175,-0.085,
     &  0.041,-0.052,-0.040, 0.033, 0.022,-0.023,
     & -0.009,-0.051,-0.065,-0.055,-0.048,-0.029,
     &  0.013, 0.000/
c Pair Ccat Oco2                                                          
      data (dfire(5,8,i),i=1,20) /
     &  4.253, 2.471, 1.106,-0.498,-1.092,-0.478,
     & -0.251,-0.325,-0.307,-0.155,-0.062,-0.095,
     & -0.075,-0.094,-0.096,-0.096,-0.082,-0.046,
     & -0.011, 0.000/
c Pair Ccat Np13                                                          
      data (dfire(5,9,i),i=1,20) /
     &  4.215, 4.191, 2.213, 0.787, 0.017, 0.029,
     &  0.015, 0.068, 0.069, 0.088, 0.056, 0.055,
     &  0.023, 0.020, 0.003,-0.022,-0.047,-0.044,
     & -0.021, 0.000/
c Pair Ccat Nam                                                           
      data (dfire(5,10,i),i=1,20) /
     &  5.268, 5.244, 3.508, 1.898, 1.023, 0.533,
     &  0.297, 0.037,-0.083,-0.146,-0.180,-0.187,
     & -0.177,-0.163,-0.171,-0.076, 0.008, 0.011,
     &  0.004, 0.000/
c Pair Ccat N4                                                            
      data (dfire(5,11,i),i=1,20) /
     &  3.409, 3.385, 3.578, 1.326, 0.589, 0.353,
     &  0.215, 0.134, 0.082, 0.192, 0.133, 0.074,
     &  0.051, 0.076, 0.055, 0.014, 0.008, 0.009,
     &  0.007, 0.000/
c Pair Ccat N2                                                            
      data (dfire(5,12,i),i=1,20) /
     &  2.954, 2.930, 3.122, 0.425,-0.188,-0.143,
     & -0.052,-0.098,-0.164, 0.040, 0.137,-0.003,
     & -0.022,-0.092,-0.102,-0.097,-0.096,-0.105,
     & -0.096, 0.000/
c Pair Ccat S3                                                            
      data (dfire(5,13,i),i=1,20) /
     &  3.360, 3.336, 3.529, 1.760, 0.264, 0.218,
     &  0.128, 0.146, 0.276, 0.309, 0.290, 0.258,
     &  0.125, 0.069,-0.008, 0.000,-0.132,-0.095,
     & -0.024, 0.000/
c Pair O2 H                                                              
      data (dfire(6,1,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair O2 C3                                                              
      data (dfire(6,2,i),i=1,20) /
     &  6.358, 4.100,-0.162, 0.205,-0.283,-0.260,
     & -0.093,-0.077, 0.036,-0.119,-0.157,-0.087,
     & -0.190,-0.167,-0.067,-0.059,-0.083,-0.067,
     & -0.030, 0.000/
c Pair O2 Car                                                             
      data (dfire(6,3,i),i=1,20) /
     &  6.918, 5.137, 3.075, 0.751, 0.189, 0.054,
     & -0.037,-0.034,-0.122,-0.149,-0.183,-0.223,
     & -0.235,-0.227,-0.160,-0.126,-0.088,-0.058,
     & -0.025, 0.000/
c Pair O2 C2                                                              
      data (dfire(6,4,i),i=1,20) /
     &  7.359, 4.005, 0.044,-0.190,-0.616,-0.386,
     &  0.027, 0.085,-0.098,-0.065,-0.039,-0.164,
     & -0.069,-0.074,-0.068,-0.088,-0.081,-0.077,
     & -0.045, 0.000/
c Pair O2 Ccat                                                            
      data (dfire(6,5,i),i=1,20) /
     &  5.248, 5.224, 2.372, 0.414,-0.059, 0.091,
     &  0.233, 0.124, 0.064, 0.012,-0.085,-0.202,
     & -0.279,-0.205,-0.110,-0.082,-0.069,-0.034,
     &  0.003, 0.000/
c Pair O2 O2                                                              
      data (dfire(6,6,i),i=1,20) /
     &  5.339, 2.990, 1.443,-0.382,-0.347,-0.115,
     & -0.574,-0.225, 0.206,-0.063, 0.008,-0.016,
     &  0.033,-0.085,-0.077,-0.078,-0.058,-0.055,
     & -0.024, 0.000/
c Pair O2 O3                                                              
      data (dfire(6,7,i),i=1,20) /
     &  6.108, 1.943, 0.067, 0.300, 0.147,-0.048,
     & -0.353,-0.111, 0.020,-0.115,-0.079,-0.112,
     & -0.147,-0.131,-0.044,-0.032,-0.038,-0.043,
     & -0.011, 0.000/
c Pair O2 Oco2                                                            
      data (dfire(6,8,i),i=1,20) /
     &  6.245, 2.836, 1.506, 0.699, 0.370, 0.193,
     & -0.097,-0.079,-0.025,-0.090,-0.129,-0.118,
     & -0.070,-0.101,-0.045,-0.030,-0.010, 0.002,
     & -0.001, 0.000/
c Pair O2 Np13                                                            
      data (dfire(6,9,i),i=1,20) /
     &  4.985, 2.512, 0.290, 0.382, 0.339, 0.195,
     &  0.037, 0.003, 0.031,-0.027,-0.128,-0.188,
     & -0.194,-0.171,-0.116,-0.078,-0.053,-0.028,
     & -0.004, 0.000/
c Pair O2 Nam                                                             
      data (dfire(6,10,i),i=1,20) /
     &  5.736, 2.963,-0.231,-0.518, 0.042, 0.164,
     &  0.061, 0.077,-0.071,-0.092,-0.226,-0.050,
     & -0.111,-0.023,-0.090,-0.071,-0.104,-0.088,
     & -0.053, 0.000/
c Pair O2 N4                                                              
      data (dfire(6,11,i),i=1,20) /
     &  3.813, 2.176, 0.080, 0.282, 0.405, 0.383,
     &  0.230, 0.152, 0.134, 0.058,-0.057,-0.158,
     & -0.177,-0.163,-0.128,-0.089,-0.071,-0.032,
     &  0.001, 0.000/
c Pair O2 N2                                                              
      data (dfire(6,12,i),i=1,20) /
     &  4.939, 3.399, 0.665, 0.581, 0.229, 0.053,
     & -0.034, 0.079, 0.055,-0.062,-0.228,-0.186,
     & -0.172,-0.134,-0.100,-0.087,-0.037,-0.031,
     & -0.003, 0.000/
c Pair O2 S3                                                              
      data (dfire(6,13,i),i=1,20) /
     &  5.272, 5.248, 2.352, 0.476, 0.023,-0.097,
     & -0.241,-0.260,-0.101,-0.129,-0.211,-0.204,
     & -0.277,-0.253,-0.153,-0.085,-0.076,-0.064,
     & -0.032, 0.000/
c Pair O3 H                                                              
      data (dfire(7,1,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair O3 C3                                                              
      data (dfire(7,2,i),i=1,20) /
     &  5.619, 5.182, 2.284, 0.578, 0.034,-0.035,
     & -0.151,-0.132,-0.119,-0.124,-0.112,-0.095,
     & -0.095,-0.087,-0.075,-0.090,-0.070,-0.034,
     & -0.011, 0.000/
c Pair O3 Car                                                             
      data (dfire(7,3,i),i=1,20) /
     &  5.804, 5.780, 3.084, 0.810, 0.305, 0.217,
     &  0.072, 0.058,-0.028,-0.062,-0.051,-0.054,
     & -0.066,-0.097,-0.109,-0.110,-0.097,-0.068,
     & -0.032, 0.000/
c Pair O3 C2                                                              
      data (dfire(7,4,i),i=1,20) /
     &  5.107, 5.083, 2.172, 0.377,-0.216,-0.334,
     & -0.113,-0.076,-0.133,-0.066,-0.126,-0.094,
     & -0.120,-0.098,-0.065,-0.046,-0.043,-0.029,
     & -0.012, 0.000/
c Pair O3 Ccat                                                            
      data (dfire(7,5,i),i=1,20) /
     &  4.115, 4.091, 1.941, 0.508,-0.175,-0.085,
     &  0.041,-0.052,-0.040, 0.033, 0.022,-0.023,
     & -0.009,-0.051,-0.065,-0.055,-0.048,-0.029,
     &  0.013, 0.000/
c Pair O3 O2                                                              
      data (dfire(7,6,i),i=1,20) /
     &  6.108, 1.943, 0.067, 0.300, 0.147,-0.048,
     & -0.353,-0.111, 0.020,-0.115,-0.079,-0.112,
     & -0.147,-0.131,-0.044,-0.032,-0.038,-0.043,
     & -0.011, 0.000/
c Pair O3 O3                                                              
      data (dfire(7,7,i),i=1,20) /
     &  4.962, 1.770,-0.117, 0.168, 0.228,-0.051,
     & -0.161,-0.080,-0.052,-0.163,-0.110,-0.101,
     & -0.086,-0.096,-0.087,-0.094,-0.069,-0.058,
     & -0.015, 0.000/
c Pair O3 Oco2                                                            
      data (dfire(7,8,i),i=1,20) /
     &  3.948, 0.727,-0.469, 0.000,-0.044,-0.081,
     & -0.206, 0.004, 0.007,-0.031,-0.059,-0.094,
     & -0.079,-0.088,-0.058,-0.052,-0.041,-0.045,
     & -0.024, 0.000/
c Pair O3 Np13                                                            
      data (dfire(7,9,i),i=1,20) /
     &  4.930, 2.275, 0.194, 0.171, 0.153, 0.008,
     & -0.070,-0.072, 0.022,-0.018,-0.045,-0.063,
     & -0.063,-0.072,-0.083,-0.079,-0.065,-0.050,
     & -0.028, 0.000/
c Pair O3 Nam                                                             
      data (dfire(7,10,i),i=1,20) /
     &  6.112, 2.850, 0.496,-0.135, 0.080, 0.046,
     & -0.209,-0.120,-0.148,-0.199,-0.191,-0.088,
     & -0.007,-0.072,-0.075,-0.062,-0.030,-0.015,
     & -0.015, 0.000/
c Pair O3 N4                                                              
      data (dfire(7,11,i),i=1,20) /
     &  4.189, 2.236, 0.178, 0.123, 0.144, 0.026,
     & -0.043,-0.011, 0.001,-0.012,-0.028,-0.021,
     & -0.005,-0.018,-0.064,-0.036,-0.025,-0.013,
     & -0.012, 0.000/
c Pair O3 N2                                                              
      data (dfire(7,12,i),i=1,20) /
     &  3.783, 1.830,-0.047, 0.152,-0.025,-0.047,
     & -0.244, 0.053, 0.080, 0.047,-0.062,-0.018,
     & -0.050,-0.074,-0.114,-0.057,-0.058,-0.056,
     & -0.018, 0.000/
c Pair O3 S3                                                              
      data (dfire(7,13,i),i=1,20) /
     &  4.118, 4.094, 2.224, 0.779, 0.148, 0.056,
     & -0.082,-0.026,-0.026,-0.060,-0.109,-0.127,
     & -0.161,-0.168,-0.171,-0.135,-0.120,-0.079,
     & -0.040, 0.000/
c Pair Oco2 H                                                            
      data (dfire(8,1,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair Oco2 C3                                                            
      data (dfire(8,2,i),i=1,20) /
     &  5.762, 4.366, 2.192, 0.552, 0.093, 0.126,
     &  0.063, 0.017, 0.010,-0.038,-0.086,-0.084,
     & -0.082,-0.088,-0.061,-0.039,-0.042,-0.033,
     & -0.020, 0.000/
c Pair Oco2 Car                                                           
      data (dfire(8,3,i),i=1,20) /
     &  5.946, 4.406, 2.860, 0.842, 0.451, 0.392,
     &  0.216, 0.208, 0.095, 0.120, 0.145, 0.127,
     &  0.100, 0.050,-0.015,-0.059,-0.069,-0.045,
     & -0.025, 0.000/
c Pair Oco2 C2                                                            
      data (dfire(8,4,i),i=1,20) /
     &  6.348, 5.220, 2.227, 0.879, 0.088, 0.043,
     &  0.049,-0.077,-0.092,-0.086,-0.085,-0.135,
     & -0.138,-0.079,-0.064,-0.027,-0.007,-0.012,
     & -0.007, 0.000/
c Pair Oco2 Ccat                                                          
      data (dfire(8,5,i),i=1,20) /
     &  4.253, 2.471, 1.106,-0.498,-1.092,-0.478,
     & -0.251,-0.325,-0.307,-0.155,-0.062,-0.095,
     & -0.075,-0.094,-0.096,-0.096,-0.082,-0.046,
     & -0.011, 0.000/
c Pair Oco2 O2                                                            
      data (dfire(8,6,i),i=1,20) /
     &  6.245, 2.836, 1.506, 0.699, 0.370, 0.193,
     & -0.097,-0.079,-0.025,-0.090,-0.129,-0.118,
     & -0.070,-0.101,-0.045,-0.030,-0.010, 0.002,
     & -0.001, 0.000/
c Pair Oco2 O3                                                            
      data (dfire(8,7,i),i=1,20) /
     &  3.948, 0.727,-0.469, 0.000,-0.044,-0.081,
     & -0.206, 0.004, 0.007,-0.031,-0.059,-0.094,
     & -0.079,-0.088,-0.058,-0.052,-0.041,-0.045,
     & -0.024, 0.000/
c Pair Oco2 Oco2                                                          
      data (dfire(8,8,i),i=1,20) /
     &  5.235, 1.669, 0.749, 0.478, 0.416, 0.202,
     &  0.042, 0.069, 0.008,-0.029,-0.082,-0.070,
     & -0.067,-0.088,-0.104,-0.069,-0.044,-0.032,
     & -0.007, 0.000/
c Pair Oco2 Np13                                                          
      data (dfire(8,9,i),i=1,20) /
     &  3.554, 1.182,-0.681,-0.458,-0.388,-0.269,
     & -0.527,-0.357,-0.127,-0.131,-0.138,-0.110,
     & -0.084,-0.084,-0.102,-0.093,-0.084,-0.053,
     & -0.029, 0.000/
c Pair Oco2 Nam                                                           
      data (dfire(8,10,i),i=1,20) /
     &  5.139, 2.822, 0.266, 0.273, 0.471, 0.316,
     &  0.019,-0.134,-0.106,-0.088,-0.167,-0.150,
     & -0.107,-0.085,-0.038,-0.015,-0.010,-0.006,
     & -0.003, 0.000/
c Pair Oco2 N4                                                            
      data (dfire(8,11,i),i=1,20) /
     &  2.852, 0.898,-0.787,-0.545,-0.425,-0.413,
     & -0.444,-0.294,-0.189,-0.174,-0.143,-0.129,
     & -0.125,-0.099,-0.094,-0.076,-0.076,-0.046,
     & -0.036, 0.000/
c Pair Oco2 N2                                                            
      data (dfire(8,12,i),i=1,20) /
     &  3.887, 2.105,-0.361,-0.098,-0.040,-0.228,
     & -0.469,-0.209,-0.091,-0.076,-0.046,-0.035,
     & -0.016,-0.045,-0.067,-0.110,-0.110,-0.055,
     & -0.026, 0.000/
c Pair Oco2 S3                                                            
      data (dfire(8,13,i),i=1,20) /
     &  4.287, 4.263, 2.114, 1.182, 0.623, 0.413,
     &  0.287, 0.266, 0.250, 0.163, 0.191, 0.131,
     &  0.095, 0.022,-0.067,-0.077,-0.066,-0.051,
     & -0.017, 0.000/
c Pair Np13 H                                                            
      data (dfire(9,1,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair Np13 C3                                                            
      data (dfire(9,2,i),i=1,20) /
     &  6.735, 5.608, 3.027, 1.250, 0.309, 0.165,
     &  0.057, 0.019, 0.027, 0.031,-0.009,-0.035,
     & -0.054,-0.099,-0.102,-0.078,-0.053,-0.028,
     & -0.015, 0.000/
c Pair Np13 Car                                                           
      data (dfire(9,3,i),i=1,20) /
     &  5.825, 4.697, 3.048, 0.862, 0.056, 0.047,
     &  0.014, 0.060, 0.037, 0.068, 0.104, 0.097,
     &  0.094, 0.056, 0.014,-0.030,-0.047,-0.050,
     & -0.033, 0.000/
c Pair Np13 C2                                                            
      data (dfire(9,4,i),i=1,20) /
     &  6.198, 6.174, 2.645, 0.867, 0.016, 0.105,
     &  0.094,-0.028,-0.095,-0.095,-0.123,-0.158,
     & -0.185,-0.169,-0.130,-0.088,-0.049,-0.019,
     & -0.010, 0.000/
c Pair Np13 Ccat                                                          
      data (dfire(9,5,i),i=1,20) /
     &  4.215, 4.191, 2.213, 0.787, 0.017, 0.029,
     &  0.015, 0.068, 0.069, 0.088, 0.056, 0.055,
     &  0.023, 0.020, 0.003,-0.022,-0.047,-0.044,
     & -0.021, 0.000/
c Pair Np13 O2                                                            
      data (dfire(9,6,i),i=1,20) /
     &  4.985, 2.512, 0.290, 0.382, 0.339, 0.195,
     &  0.037, 0.003, 0.031,-0.027,-0.128,-0.188,
     & -0.194,-0.171,-0.116,-0.078,-0.053,-0.028,
     & -0.004, 0.000/
c Pair Np13 O3                                                            
      data (dfire(9,7,i),i=1,20) /
     &  4.930, 2.275, 0.194, 0.171, 0.153, 0.008,
     & -0.070,-0.072, 0.022,-0.018,-0.045,-0.063,
     & -0.063,-0.072,-0.083,-0.079,-0.065,-0.050,
     & -0.028, 0.000/
c Pair Np13 Oco2                                                          
      data (dfire(9,8,i),i=1,20) /
     &  3.554, 1.182,-0.681,-0.458,-0.388,-0.269,
     & -0.527,-0.357,-0.127,-0.131,-0.138,-0.110,
     & -0.084,-0.084,-0.102,-0.093,-0.084,-0.053,
     & -0.029, 0.000/
c Pair Np13 Np13                                                          
      data (dfire(9,9,i),i=1,20) /
     &  5.022, 4.998, 1.216, 0.423, 0.055, 0.023,
     & -0.038, 0.050, 0.055, 0.038, 0.029, 0.040,
     &  0.018,-0.020,-0.020,-0.037,-0.060,-0.048,
     & -0.029, 0.000/
c Pair Np13 Nam                                                           
      data (dfire(9,10,i),i=1,20) /
     &  6.094, 4.141, 1.868, 1.207, 0.628, 0.344,
     &  0.103, 0.019,-0.083,-0.132,-0.191,-0.195,
     & -0.173,-0.163,-0.138,-0.080,-0.031,-0.004,
     & -0.007, 0.000/
c Pair Np13 N4                                                            
      data (dfire(9,11,i),i=1,20) /
     &  4.221, 3.094, 1.221, 0.941, 0.507, 0.324,
     &  0.121, 0.161, 0.134, 0.112, 0.112, 0.127,
     &  0.089, 0.067, 0.038, 0.009,-0.003,-0.013,
     &  0.000, 0.000/
c Pair Np13 N2                                                            
      data (dfire(9,12,i),i=1,20) /
     &  3.824, 3.800, 0.883, 0.306,-0.117,-0.119,
     & -0.230,-0.096, 0.002,-0.130, 0.070,-0.011,
     & -0.034,-0.080,-0.079,-0.069,-0.077,-0.074,
     & -0.032, 0.000/
c Pair Np13 S3                                                            
      data (dfire(9,13,i),i=1,20) /
     &  4.173, 4.149, 2.584, 0.655, 0.113, 0.214,
     &  0.080, 0.061, 0.119, 0.186, 0.130, 0.074,
     &  0.064, 0.005,-0.039,-0.080,-0.099,-0.089,
     & -0.043, 0.000/
c Pair Nam H                                                             
      data (dfire(10,1,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair Nam C3                                                             
      data (dfire(10,2,i),i=1,20) /
     &  6.376, 5.336, 1.767, 0.155, 0.523,-0.057,
     & -0.251,-0.277,-0.272,-0.143,-0.095,-0.104,
     & -0.029,-0.053,-0.105,-0.096,-0.056,-0.018,
     & -0.010, 0.000/
c Pair Nam Car                                                            
      data (dfire(10,3,i),i=1,20) /
     &  6.922, 6.898, 3.464, 1.191, 0.383, 0.282,
     &  0.032,-0.064,-0.178,-0.237,-0.267,-0.259,
     & -0.220,-0.203,-0.174,-0.130,-0.081,-0.046,
     & -0.021, 0.000/
c Pair Nam C2                                                             
      data (dfire(10,4,i),i=1,20) /
     &  6.283, 5.191, 1.805,-0.123,-0.168,-0.091,
     & -0.099,-0.065,-0.290,-0.109, 0.011, 0.044,
     &  0.005,-0.029,-0.053,-0.099,-0.054,-0.035,
     & -0.015, 0.000/
c Pair Nam Ccat                                                           
      data (dfire(10,5,i),i=1,20) /
     &  5.268, 5.244, 3.508, 1.898, 1.023, 0.533,
     &  0.297, 0.037,-0.083,-0.146,-0.180,-0.187,
     & -0.177,-0.163,-0.171,-0.076, 0.008, 0.011,
     &  0.004, 0.000/
c Pair Nam O2                                                             
      data (dfire(10,6,i),i=1,20) /
     &  5.736, 2.963,-0.231,-0.518, 0.042, 0.164,
     &  0.061, 0.077,-0.071,-0.092,-0.226,-0.050,
     & -0.111,-0.023,-0.090,-0.071,-0.104,-0.088,
     & -0.053, 0.000/
c Pair Nam O3                                                             
      data (dfire(10,7,i),i=1,20) /
     &  6.112, 2.850, 0.496,-0.135, 0.080, 0.046,
     & -0.209,-0.120,-0.148,-0.199,-0.191,-0.088,
     & -0.007,-0.072,-0.075,-0.062,-0.030,-0.015,
     & -0.015, 0.000/
c Pair Nam Oco2                                                           
      data (dfire(10,8,i),i=1,20) /
     &  5.139, 2.822, 0.266, 0.273, 0.471, 0.316,
     &  0.019,-0.134,-0.106,-0.088,-0.167,-0.150,
     & -0.107,-0.085,-0.038,-0.015,-0.010,-0.006,
     & -0.003, 0.000/
c Pair Nam Np13                                                           
      data (dfire(10,9,i),i=1,20) /
     &  6.094, 4.141, 1.868, 1.207, 0.628, 0.344,
     &  0.103, 0.019,-0.083,-0.132,-0.191,-0.195,
     & -0.173,-0.163,-0.138,-0.080,-0.031,-0.004,
     & -0.007, 0.000/
c Pair Nam Nam                                                            
      data (dfire(10,10,i),i=1,20) /
     &  5.338, 2.927,-0.787, 0.024,-0.185,-0.251,
     & -0.318,-0.055,-0.029,-0.184,-0.015, 0.107,
     &  0.053,-0.016,-0.140,-0.091,-0.070,-0.067,
     & -0.033, 0.000/
c Pair Nam N4                                                             
      data (dfire(10,11,i),i=1,20) /
     &  5.337, 3.142, 0.855, 1.117, 0.758, 0.530,
     &  0.190, 0.194, 0.106, 0.026,-0.047,-0.158,
     & -0.168,-0.178,-0.172,-0.104,-0.009,-0.003,
     & -0.004, 0.000/
c Pair Nam N2                                                             
      data (dfire(10,12,i),i=1,20) /
     &  4.937, 4.913, 1.956, 1.475, 0.705, 0.272,
     &  0.016,-0.051,-0.093,-0.166,-0.304,-0.134,
     & -0.142,-0.182,-0.096,-0.077,-0.034,-0.034,
     & -0.017, 0.000/
c Pair Nam S3                                                             
      data (dfire(10,13,i),i=1,20) /
     &  5.282, 5.258, 2.775, 0.838, 0.198, 0.224,
     & -0.152,-0.217,-0.278,-0.302,-0.331,-0.304,
     & -0.217,-0.148,-0.142,-0.125,-0.068,-0.029,
     & -0.019, 0.000/
c Pair N4 H                                                              
      data (dfire(11,1,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair N4 C3                                                              
      data (dfire(11,2,i),i=1,20) /
     &  4.243, 4.219, 2.199, 1.144, 0.820, 0.474,
     &  0.242, 0.140, 0.109, 0.098, 0.077, 0.048,
     &  0.021,-0.045,-0.090,-0.061,-0.038,-0.004,
     & -0.005, 0.000/
c Pair N4 Car                                                             
      data (dfire(11,3,i),i=1,20) /
     &  5.102, 3.975, 2.858, 1.304, 0.635, 0.365,
     &  0.205, 0.160, 0.179, 0.180, 0.203, 0.222,
     &  0.257, 0.236, 0.157, 0.084, 0.015,-0.022,
     & -0.027, 0.000/
c Pair N4 C2                                                              
      data (dfire(11,4,i),i=1,20) /
     &  4.346, 4.322, 2.086, 0.645, 0.035, 0.199,
     &  0.260, 0.144, 0.039, 0.016,-0.020,-0.079,
     & -0.165,-0.162,-0.137,-0.086,-0.042,-0.006,
     &  0.005, 0.000/
c Pair N4 Ccat                                                            
      data (dfire(11,5,i),i=1,20) /
     &  3.409, 3.385, 3.578, 1.326, 0.589, 0.353,
     &  0.215, 0.134, 0.082, 0.192, 0.133, 0.074,
     &  0.051, 0.076, 0.055, 0.014, 0.008, 0.009,
     &  0.007, 0.000/
c Pair N4 O2                                                              
      data (dfire(11,6,i),i=1,20) /
     &  3.813, 2.176, 0.080, 0.282, 0.405, 0.383,
     &  0.230, 0.152, 0.134, 0.058,-0.057,-0.158,
     & -0.177,-0.163,-0.128,-0.089,-0.071,-0.032,
     &  0.001, 0.000/
c Pair N4 O3                                                              
      data (dfire(11,7,i),i=1,20) /
     &  4.189, 2.236, 0.178, 0.123, 0.144, 0.026,
     & -0.043,-0.011, 0.001,-0.012,-0.028,-0.021,
     & -0.005,-0.018,-0.064,-0.036,-0.025,-0.013,
     & -0.012, 0.000/
c Pair N4 Oco2                                                            
      data (dfire(11,8,i),i=1,20) /
     &  2.852, 0.898,-0.787,-0.545,-0.425,-0.413,
     & -0.444,-0.294,-0.189,-0.174,-0.143,-0.129,
     & -0.125,-0.099,-0.094,-0.076,-0.076,-0.046,
     & -0.036, 0.000/
c Pair N4 Np13                                                            
      data (dfire(11,9,i),i=1,20) /
     &  4.221, 3.094, 1.221, 0.941, 0.507, 0.324,
     &  0.121, 0.161, 0.134, 0.112, 0.112, 0.127,
     &  0.089, 0.067, 0.038, 0.009,-0.003,-0.013,
     &  0.000, 0.000/
c Pair N4 Nam                                                             
      data (dfire(11,10,i),i=1,20) /
     &  5.337, 3.142, 0.855, 1.117, 0.758, 0.530,
     &  0.190, 0.194, 0.106, 0.026,-0.047,-0.158,
     & -0.168,-0.178,-0.172,-0.104,-0.009,-0.003,
     & -0.004, 0.000/
c Pair N4 N4                                                              
      data (dfire(11,11,i),i=1,20) /
     &  3.634, 2.094, 0.977, 1.075, 0.716, 0.552,
     &  0.136, 0.134, 0.274, 0.088, 0.047,-0.021,
     & -0.024, 0.033, 0.050,-0.008,-0.005,-0.037,
     & -0.017, 0.000/
c Pair N4 N2                                                              
      data (dfire(11,12,i),i=1,20) /
     &  3.052, 3.028, 1.291, 0.797, 0.369, 0.254,
     & -0.072, 0.063, 0.253, 0.088, 0.102, 0.195,
     &  0.216, 0.118, 0.035, 0.032, 0.016,-0.025,
     & -0.047, 0.000/
c Pair N4 S3                                                              
      data (dfire(11,13,i),i=1,20) /
     &  3.447, 3.423, 1.858, 1.364, 0.764, 0.481,
     &  0.293, 0.251, 0.348, 0.369, 0.417, 0.265,
     &  0.320, 0.200, 0.134, 0.022,-0.027,-0.017,
     & -0.003, 0.000/
c Pair N2 H                                                              
      data (dfire(12,1,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair N2 C3                                                              
      data (dfire(12,2,i),i=1,20) /
     &  5.574, 5.550, 3.331, 1.210, 0.233, 0.125,
     & -0.003,-0.066,-0.033, 0.013, 0.029, 0.004,
     & -0.047,-0.142,-0.115,-0.074,-0.040,-0.030,
     & -0.020, 0.000/
c Pair N2 Car                                                             
      data (dfire(12,3,i),i=1,20) /
     &  4.675, 4.651, 3.328, 0.974, 0.068,-0.027,
     & -0.096, 0.009, 0.038, 0.083, 0.081, 0.104,
     &  0.056,-0.026,-0.059,-0.062,-0.048,-0.052,
     & -0.030, 0.000/
c Pair N2 C2                                                              
      data (dfire(12,4,i),i=1,20) /
     &  5.049, 5.025, 3.701, 0.973, 0.210, 0.072,
     & -0.133,-0.143,-0.128,-0.022,-0.160,-0.148,
     & -0.159,-0.169,-0.107,-0.070,-0.049,-0.022,
     & -0.022, 0.000/
c Pair N2 Ccat                                                            
      data (dfire(12,5,i),i=1,20) /
     &  2.954, 2.930, 3.122, 0.425,-0.188,-0.143,
     & -0.052,-0.098,-0.164, 0.040, 0.137,-0.003,
     & -0.022,-0.092,-0.102,-0.097,-0.096,-0.105,
     & -0.096, 0.000/
c Pair N2 O2                                                              
      data (dfire(12,6,i),i=1,20) /
     &  4.939, 3.399, 0.665, 0.581, 0.229, 0.053,
     & -0.034, 0.079, 0.055,-0.062,-0.228,-0.186,
     & -0.172,-0.134,-0.100,-0.087,-0.037,-0.031,
     & -0.003, 0.000/
c Pair N2 O3                                                              
      data (dfire(12,7,i),i=1,20) /
     &  3.783, 1.830,-0.047, 0.152,-0.025,-0.047,
     & -0.244, 0.053, 0.080, 0.047,-0.062,-0.018,
     & -0.050,-0.074,-0.114,-0.057,-0.058,-0.056,
     & -0.018, 0.000/
c Pair N2 Oco2                                                            
      data (dfire(12,8,i),i=1,20) /
     &  3.887, 2.105,-0.361,-0.098,-0.040,-0.228,
     & -0.469,-0.209,-0.091,-0.076,-0.046,-0.035,
     & -0.016,-0.045,-0.067,-0.110,-0.110,-0.055,
     & -0.026, 0.000/
c Pair N2 Np13                                                            
      data (dfire(12,9,i),i=1,20) /
     &  3.824, 3.800, 0.883, 0.306,-0.117,-0.119,
     & -0.230,-0.096, 0.002,-0.130, 0.070,-0.011,
     & -0.034,-0.080,-0.079,-0.069,-0.077,-0.074,
     & -0.032, 0.000/
c Pair N2 Nam                                                             
      data (dfire(12,10,i),i=1,20) /
     &  4.937, 4.913, 1.956, 1.475, 0.705, 0.272,
     &  0.016,-0.051,-0.093,-0.166,-0.304,-0.134,
     & -0.142,-0.182,-0.096,-0.077,-0.034,-0.034,
     & -0.017, 0.000/
c Pair N2 N4                                                              
      data (dfire(12,11,i),i=1,20) /
     &  3.052, 3.028, 1.291, 0.797, 0.369, 0.254,
     & -0.072, 0.063, 0.253, 0.088, 0.102, 0.195,
     &  0.216, 0.118, 0.035, 0.032, 0.016,-0.025,
     & -0.047, 0.000/
c Pair N2 N2                                                              
      data (dfire(12,12,i),i=1,20) /
     &  2.682, 2.658,-0.146,-0.657,-0.787,-0.677,
     & -0.409,-0.129,-0.324,-0.119,-0.045,-0.054,
     & -0.159,-0.159,-0.123,-0.127,-0.099,-0.151,
     & -0.107, 0.000/
c Pair N2 S3                                                              
      data (dfire(12,13,i),i=1,20) /
     &  3.037, 3.013, 1.689, 0.575,-0.019,-0.307,
     & -0.190,-0.139,-0.032,-0.217,-0.137,-0.124,
     & -0.069,-0.098,-0.105,-0.134,-0.114,-0.047,
     &  0.000, 0.000/
c Pair S3 H                                                              
      data (dfire(13,1,i),i=1,20) /
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  0.000, 0.000/
c Pair S3 C3                                                              
      data (dfire(13,2,i),i=1,20) /
     &  4.378, 4.766, 1.896, 0.773, 0.077,-0.124,
     & -0.170,-0.230,-0.186,-0.220,-0.187,-0.175,
     & -0.169,-0.197,-0.175,-0.137,-0.105,-0.061,
     & -0.028, 0.000/
c Pair S3 Car                                                             
      data (dfire(13,3,i),i=1,20) /
     &  3.415, 4.907, 3.037, 1.148,-0.118,-0.284,
     & -0.353,-0.397,-0.292,-0.265,-0.182,-0.117,
     & -0.117,-0.162,-0.202,-0.199,-0.170,-0.117,
     & -0.061, 0.000/
c Pair S3 C2                                                              
      data (dfire(13,4,i),i=1,20) /
     &  4.280, 5.359, 3.490, 1.252, 0.233,-0.145,
     & -0.190,-0.219,-0.214,-0.171,-0.184,-0.297,
     & -0.236,-0.198,-0.156,-0.105,-0.075,-0.052,
     & -0.030, 0.000/
c Pair S3 Ccat                                                            
      data (dfire(13,5,i),i=1,20) /
     &  3.360, 3.336, 3.529, 1.760, 0.264, 0.218,
     &  0.128, 0.146, 0.276, 0.309, 0.290, 0.258,
     &  0.125, 0.069,-0.008, 0.000,-0.132,-0.095,
     & -0.024, 0.000/
c Pair S3 O2                                                              
      data (dfire(13,6,i),i=1,20) /
     &  5.272, 5.248, 2.352, 0.476, 0.023,-0.097,
     & -0.241,-0.260,-0.101,-0.129,-0.211,-0.204,
     & -0.277,-0.253,-0.153,-0.085,-0.076,-0.064,
     & -0.032, 0.000/
c Pair S3 O3                                                              
      data (dfire(13,7,i),i=1,20) /
     &  4.118, 4.094, 2.224, 0.779, 0.148, 0.056,
     & -0.082,-0.026,-0.026,-0.060,-0.109,-0.127,
     & -0.161,-0.168,-0.171,-0.135,-0.120,-0.079,
     & -0.040, 0.000/
c Pair S3 Oco2                                                            
      data (dfire(13,8,i),i=1,20) /
     &  4.287, 4.263, 2.114, 1.182, 0.623, 0.413,
     &  0.287, 0.266, 0.250, 0.163, 0.191, 0.131,
     &  0.095, 0.022,-0.067,-0.077,-0.066,-0.051,
     & -0.017, 0.000/
c Pair S3 Np13                                                            
      data (dfire(13,9,i),i=1,20) /
     &  4.173, 4.149, 2.584, 0.655, 0.113, 0.214,
     &  0.080, 0.061, 0.119, 0.186, 0.130, 0.074,
     &  0.064, 0.005,-0.039,-0.080,-0.099,-0.089,
     & -0.043, 0.000/
c Pair S3 Nam                                                             
      data (dfire(13,10,i),i=1,20) /
     &  5.282, 5.258, 2.775, 0.838, 0.198, 0.224,
     & -0.152,-0.217,-0.278,-0.302,-0.331,-0.304,
     & -0.217,-0.148,-0.142,-0.125,-0.068,-0.029,
     & -0.019, 0.000/
c Pair S3 N4                                                              
      data (dfire(13,11,i),i=1,20) /
     &  3.447, 3.423, 1.858, 1.364, 0.764, 0.481,
     &  0.293, 0.251, 0.348, 0.369, 0.417, 0.265,
     &  0.320, 0.200, 0.134, 0.022,-0.027,-0.017,
     & -0.003, 0.000/
c Pair S3 N2                                                              
      data (dfire(13,12,i),i=1,20) /
     &  3.037, 3.013, 1.689, 0.575,-0.019,-0.307,
     & -0.190,-0.139,-0.032,-0.217,-0.137,-0.124,
     & -0.069,-0.098,-0.105,-0.134,-0.114,-0.047,
     &  0.000, 0.000/
c Pair S3 S3                                                              
      data (dfire(13,13,i),i=1,20) /
     & -0.087,-1.852, 0.713, 0.368,-0.612,-0.504,
     & -0.607,-0.535,-0.210,-0.195,-0.236,-0.202,
     & -0.299,-0.247,-0.260,-0.212,-0.195,-0.121,
     & -0.044, 0.000/

      return
      end

