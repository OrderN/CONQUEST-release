      subroutine mulddd(vdwr,moddma,idomul,idodma,p,qat,iconn)

c THIS IS REALLY muldma

      implicit double precision(a-h,o-z)
c
c     Mulliken Charges and/or Distributed Multipole Analysis
c
c
c max rank hexadecapole, l=4
c
c explicit energy expressions only available upto hexadecapole
c d functions contribute upto l=4
c we will be missing some higher contributions from f-functions
c
c total number of functions is (lmax+1)**2
c
c for l you have k=-l,-l+1,..,0,l-1,l
c for k.ne.0 the Qlk (or Rlk) are complex:
c
c Ql, k = 1/sqrt(2)(-1)**k {Qlkc + i Qlks}
c Ql,-k = 1/sqrt(2) {Qlkc - i Qlks}
c
c same thing for the regular solid harmonics Rlk, both expressions
c are used shifting overlap density contributions to atomic sites
c
c an s-p overlap kan have a maximum rank of 1(=0+1) but the shifted
c multipole doesnt converge at the same rank, but much higher
c
c moddma
c        = 0 stadard dma; all overlap density shifted to atoms
c        = 1 sites consist of atoms + points halfway bonds
c        = 2 sites consist of atoms + centers of gaussian products
c

      parameter (mxcon=10)
      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      parameter (mxel=100)
      parameter (mxsite=300)
      parameter (lmax=4)
      parameter (lmaxf=(lmax+1)*(lmax+1))
      parameter (mxint=225)
      logical debug,dohb
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /gauhlp/ d6to5(6,5),f10to7(10,7),g15to9(15,9),
     &                plmn(35),lmn(3,35)
      common /orbhlp/ mxorb,iuhf,ispd
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /pseudo/ ipseud,ivale(numatm)
      common /coord / xyz(3,numatm)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      character*8   ctag,ctagd
      common /multip/ qmom(lmaxf,mxsite),car(3,mxsite),ctag(mxsite),
     &                nsites
      common /setp/   pgau(3),p1(3),p2(3),gam,pigam2,fac
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /athlp/  iatoms, mxnat

      dimension co1(3),co2(3),ca(35),cb(35),s2(mxint,34)
      dimension p(*),q(lmaxf), vdwr(*), temp(3),dp(3),dpn(3),ltmp(3)
      dimension shlp(mxint),shlpd(mxint,3)
      dimension qat(*),iconn(mxcon+1,*)

      if (isgau.ne.1) return

      debug = .false.
      dohb = .true.
      bohr = 0.529177249d+00
      rt3 = dsqrt(3.0d0)
      rt5 = dsqrt(5.0d0)
      rt8 = dsqrt(8.0d0)
      rt35 = dsqrt(35.0d0)
      cf = 5.2917715d-11*1.601917d-19/3.33564d-30
      ctagd = '        '

      call denini
      call setcon

      if (idomul.eq.1) then
         do i=1,3
            dp(i) = 0.0d0
            dpn(i) = 0.0d0
         end do
         do i=1,natoms
            if (ipseud.eq.1) then
               qat(i) = dble(ivale(i))
            else
               if (nat(i).eq.99) then
                  qat(i) = 0.0d0
               else
                  qat(i) = dble(nat(i))
               endif
            endif
            do j=1,3
               dpn(j) = dpn(j) + qat(i)*xyz(j,i)
            end do
         end do
         ihasq = 1
      endif

      if (idodma.eq.1) then
         nsites = natoms
         if (nsites.gt.mxsite) then
            nsites = mxsite
            write(iun3,*) 'exceded max. # of sites, increase mxsite'
         endif
         do i=1,natoms
            if (ipseud.eq.1) then
               qmom(1,i) = dble(ivale(i))
            else
               qmom(1,i) = dble(nat(i))
            endif
            ctag(i) = ctagd
            ctag(i) = elemnt(nat(i))
            do j=2,lmaxf
               qmom(j,i) = 0.0d0
            end do
            do j=1,3
               car(j,i) = xyz(j,i)
            end do
         end do

         if (moddma.eq.1) then
            call xyzcoo(1,0,0)
            call doconn
            if (dohb) call dohcon(0)
            do i=1,natoms
                do k=1,iconn(1,i)
                   j = abs(iconn(1+k,i))
c               do j=i+1,natoms
c                   ttot = 0.0d0
c                   do k=1,3
c                       temp(k) = xyz(k,j) - xyz(k,i)
c                       ttot = ttot + temp(k)*temp(k)
c                   end do
c                   dvdw = (vdwr(nat(i))+vdwr(nat(j))) / bohr
c                   dvdw2 = dvdw*dvdw
c                   if (ttot.lt.dvdw2) then
                   if (j.gt.i) then
                       if (nsites.lt.mxsite) then
                          nsites = nsites + 1
                          ctag(nsites) = ctagd
                          ctag(nsites)(1:2) = elemnt(nat(i))
                          ctag(nsites)(3:4) = ' -'
                          ctag(nsites)(5:6) = elemnt(nat(j))
                          do l=1,3
                             temp(l) = xyz(l,j) - xyz(l,i)
                             car(l,nsites) = xyz(l,i) + temp(l)/2.0d0
                          end do
                          do l=1,lmaxf
                             qmom(l,nsites) = 0.0d0
                          end do
                        else
                          write(iun3,*) 'exceded max. # of sites'
                        endif
                   endif
               end do
            end do
         endif
      endif
c
c     loop over shell pairs.
c

      iao = 0
      do ishell = 1, nshell
        isha = shella(ishell)
        isht = shellt(ishell)
        ishc = shellc(ishell)
        ishn = shelln(ishell)
        iat = jan(ishell)
        co1(1) = gx(ishell)
        co1(2) = gy(ishell)
        co1(3) = gz(ishell)
        call qtype(isht,ishc,ibeg,iend)
        jao = 0
        iendt = iend

        do jshell = 1,ishell
c        do jshell = 1,nshell
          jsha = shella(jshell)
          jsht = shellt(jshell)
          jshc = shellc(jshell)
          jshn = shelln(jshell)
          jat = jan(jshell)
          co2(1) = gx(jshell)
          co2(2) = gy(jshell)
          co2(3) = gz(jshell)
          call qtype(jsht,jshc,jbeg,jend)

          nrank = isht+jsht
          if (nrank.gt.lmax) nrank = lmax

          s1or2 = 2.0d0
          if (ishell.eq.jshell) s1or2 = 1.0d0

          jendt = jend

c===================================================
c         i loop over primitive gaussians

          do igauss=1,ishn
             ex1 = exx(isha+igauss-1)
             call fcij(isht,isha,isha+igauss-1,shladf(ishell),ca)

c            j loop over primitive gaussians

             
             jdim = jshn
             if (ishell.eq.jshell) jdim = igauss

             do jgauss=1,jdim
c             do jgauss=1,jshn
                iend = iendt
                jend = jendt
                ex2 = exx(jsha+jgauss-1)
                call fcij(jsht,jsha,jsha+jgauss-1,shladf(jshell),cb)
                call psetup(ex1,co1,ex2,co2)
                f1or2 = 1.0d0
                if (ishell.eq.jshell.and.igauss.ne.jgauss) f1or2 = 2.0d0

c
c loop over basis functions
c
              do i=1,lmaxf
                q(i) = 0.0d0
              end do

              kk = 0
              do ifun=ibeg,iend
                 do jfun=jbeg,jend
                   kk = kk + 1

                   call sint(lmn(1,ifun),lmn(1,jfun),s)
                   shlp(kk) = ca(ifun)*cb(jfun)*s

                   do id=1,3
                      do idd=1,3
                         ltmp(idd) = lmn(idd,ifun)
                         if (idd.eq.id) ltmp(idd) = ltmp(idd) + 1
                      end do
                      call sint(ltmp,lmn(1,jfun),s0)
                      shlpd(kk,id) = ca(ifun)*cb(jfun)*s0
                   end do

                   if (idodma.eq.1) then
                      nn = 1
                      if (nrank.ge.1) nn = 3
                      if (nrank.ge.2) nn = 9
                      if (nrank.ge.3) nn = 19
                      if (nrank.ge.4) nn = 34

                      do i=1,nn
                         call 
     &                dsint(lmn(1,ifun),lmn(1,jfun),lmn(1,1+i),ss)
                         s2(kk,i) = ca(ifun)*cb(jfun)*ss
                      end do

                   endif

                 end do
              end do

              call purdf(isht,jsht,ibeg,jbeg,iend,jend,shlp)
              iend = iendt
              jend = jendt
              call purdf(isht,jsht,ibeg,jbeg,iend,jend,shlpd(1,1))
              iend = iendt
              jend = jendt
              call purdf(isht,jsht,ibeg,jbeg,iend,jend,shlpd(1,2))
              iend = iendt
              jend = jendt
              call purdf(isht,jsht,ibeg,jbeg,iend,jend,shlpd(1,3))

              if (idodma.eq.1) then
                 nn = 1
                 if (nrank.ge.1) nn = 3
                 if (nrank.ge.2) nn = 9
                 if (nrank.ge.3) nn = 19
                 if (nrank.ge.4) nn = 34
                 do i=1,nn
                    iend = iendt
                    jend = jendt
                    call purdf(isht,jsht,ibeg,jbeg,iend,jend,s2(1,i))
                 end do
              endif

              kk = 0
              iiao = iao
              do ifun=ibeg,iend
               iiao = iiao + 1
               jjao = jao
              do jfun=jbeg,jend
               jjao = jjao + 1
                kk = kk + 1

                fcom = ca(ifun)*cb(jfun)*p((jjao-1)*mxorb+iiao)
     &                 *s1or2*f1or2
                fcomh = p((jjao-1)*mxorb+iiao)
     &                 *s1or2*f1or2

c q0
                call sint(lmn(1,ifun),lmn(1,jfun),s)
                q(1) = q(1) - fcomh*shlp(kk)

                do id=1,3
                   dp(id) = dp(id) - 
     &               (shlpd(kk,id) + co1(id)*shlp(kk))*fcomh
                end do

                if (idodma.eq.1) then

                   if (nrank.ge.1) then
c q10,q11c,q11s
                      q(2) = q(2) - fcomh*s2(kk,3)
   
                      q(3) = q(3) - fcomh*s2(kk,1)

                      q(4) = q(4) - fcomh*s2(kk,2)
                   endif


                   if (nrank.ge.2) then
c q20,q21c,q21s,q22c,q22s

                      q(5) = q(5) - 
     &                       fcomh*(s2(kk,6)-0.5d0*(s2(kk,4)+s2(kk,5)))
                      q(6) = q(6) - rt3*fcomh*s2(kk,8)
                      q(7) = q(7) - rt3*fcomh*s2(kk,9)
                      q(8) = q(8) - 0.5d0*rt3*fcomh*(s2(kk,4)-s2(kk,5))
                      q(9) = q(9) - rt3*fcomh*s2(kk,7)
                   endif

                   if (nrank.ge.3) then
c q30,q31c,q31s,q32c,q32s,q33c,q33s

                      q(10) = q(10) - fcomh*
     &                        (s2(kk,12) - 1.5d0*(s2(kk,15)+ s2(kk,18)))
                      q(11) = q(11) - fcomh*rt3/rt8*
     &                        (4.0d0*s2(kk,16) - s2(kk,10) - s2(kk,13))
                      q(12) = q(12) - fcomh*rt3/rt8*
     &                        (4.0d0*s2(kk,17) - s2(kk,14) - s2(kk,11))
                      q(13) = q(13) - fcomh*0.5d0*rt3*rt5*
     &                        (s2(kk,15) - s2(kk,18))
                      q(14) = q(14) - fcomh*rt3*rt5*s2(kk,19)
                      q(15) = q(15) - fcomh*rt5/rt8*
     &                        (      s2(kk,10) - 3.0d0*s2(kk,13))
                      q(16) = q(16) - fcomh*rt5/rt8*
     &                        (3.0d0*s2(kk,14) -       s2(kk,11))
                   endif

                   if (nrank.ge.4) then
c q40,q41c,q41s,q42c,q42s,q43c,q43s,q44c,q44s

                      q(17) = q(17) - fcomh*0.125d0*
     &                     (8.0d0*s2(kk,22)+3.0d0*(s2(kk,20)+s2(kk,21))+
     &                     6.0d0*s2(kk,29)-24.0d0*(s2(kk,30)+s2(kk,31)))
                      q(18) = q(18) - fcomh*rt5/rt8*
     &                     (4.0d0*s2(kk,27)-3.0d0*(s2(kk,24)+s2(kk,33)))
                      q(19) = q(19) - fcomh*rt5/rt8*
     &                     (4.0d0*s2(kk,28)-3.0d0*(s2(kk,26)+s2(kk,32)))
                      q(20) = q(20) - fcomh*0.25d0*rt5*
     &                     (s2(kk,21)-s2(kk,20)+
     &                      6.0d0*(s2(kk,30)-s2(kk,31)))
                      q(21) = q(21) - fcomh*0.5d0*rt5*
     &                       (6.0d0*s2(kk,34)-s2(kk,23)-s2(kk,25))
                      q(22) = q(22) - fcomh*rt35/rt8*
     &                       (s2(kk,24)-3.0d0*s2(kk,33))
                      q(23) = q(23) - fcomh*rt35/rt8*
     &                       (3.0d0*s2(kk,32)-s2(kk,26))
                      q(24) = q(24) - fcomh*0.125d0*rt35*
     &                       (s2(kk,20)+s2(kk,21)-6.0d0*s2(kk,29))
                      q(25) = q(25) - fcomh*0.5d0*rt35*
     &                       (s2(kk,23)-s2(kk,25))
                   endif
                end if

               end do
               end do
c--- end loop basis functions

               if (debug) then
                  write(iun3,'(6i4)') iat,jat,
     &              ishell,jshell,isha+igauss-1,jsha+jgauss-1
                  write(iun3,'(f10.6)') q(1)
                  if (idodma.eq.1) then
                    write(iun3,'(3f10.6)') (q(i),i=2,4)
                    write(iun3,'(5f10.6)') (q(i),i=5,9)
                    write(iun3,'(7f10.6)') (q(i),i=10,16)
                    write(iun3,'(9f10.6)') (q(i),i=17,25)
                  endif
               endif

               if (idomul.eq.1) then
                  if (iat.eq.jat) then
                      qat(iat) = qat(iat) + q(1)
                  else
                      qat(iat) = qat(iat) + 0.5d0*q(1)
                      qat(jat) = qat(jat) + 0.5d0*q(1)
                  endif
               endif
               if (idodma.eq.1) then
                  if (iat.eq.jat) then
c
c atomic contribution
c
                      do i=1,(nrank+1)**2
                         qmom(i,iat) = qmom(i,iat) + q(i)
                      end do
                  else
c
c partition and shift overlap contribution
c
                      if (moddma.le.1.or.
     &                   (moddma.eq.2.and.nrank.eq.0)) then
                         call partq(q,nrank,pgau)
                      else
                        ttot = 0.0d0
                        do k=1,3
                            temp(k) = xyz(k,jat) - xyz(k,iat)
                            ttot = ttot + temp(k)*temp(k)
                        end do
                        dvdw = (vdwr(nat(iat))+vdwr(nat(jat))) / bohr
                        dvdw2 = dvdw*dvdw
                        if (ttot.lt.dvdw2) then
                          if (nsites.lt.mxsite) then
c
c Create New site at P
c
                            nsites = nsites + 1
                            ctag(nsites) = ctagd
                            ctag(nsites)(1:2) = elemnt(nat(iat))
                            ctag(nsites)(3:4) = ' -'
                            ctag(nsites)(5:6) = elemnt(nat(jat))
                            do k=1,3
                               car(k,nsites) = pgau(k)
                            end do
                            do k=1,lmaxf
                               qmom(k,nsites) = q(k)
                            end do
                          else
                            if (nsites.ge.mxsite) then
            write(iun3,*) 'exceded max. # of sites, increase mxsite'
                            endif
                            call partq(q,nrank,pgau)
                          endif
                        else
                            call partq(q,nrank,pgau)
                        endif
                      endif
                  endif
               endif
             end do
          end do
c--- end loops primitives

c===================================================

        jao = jao + (jend-jbeg+1)
        jend = jendt
        end do
        iao = iao + (iend-ibeg+1)
        iend = iendt
      end do
c---- end loop shells

      if (idomul.eq.1) then
         write(iun3,'(a)') ' '
         write(iun3,'(24x,a)') 'Mulliken Charges'
         write(iun3,'(24x,a)') '================'
         write(iun3,'(a)') ' '
         qtot = 0.0d0
         do i=1,natoms
            write(iun3,'(14x,a2,a,f8.4)') elemnt(nat(i)),' ',qat(i)
c            write(iun3,'(14x,a2,a,f10.6)') elemnt(nat(i)),' ',qat(i)
            qtot = qtot + qat(i)
         end do

         write(iun3,'(a)') ' '
         write(iun3,'(a,f8.4)') 'Sum of Mulliken Charges: ',qtot

         if (.true.) then
         write(iun3,'(a)') ' '
         write(iun3,'(24x,a)') 'Dipole moment (Debye)'
         write(iun3,'(24x,a)') '======================='
         write(iun3,'(17x,a)') 
     &            '      X           Y           Z          Scalar'
         write(iun3,'(a)') ' '
         write(iun3,'(a,4(f11.4,1x))') 'Nuclear          ',
     &       (dpn(i)*cf,i=1,3),vlen(dpn)*cf
         write(iun3,'(a,4(f11.4,1x))') 'Electronic       ',
     &       (dp(i)*cf,i=1,3),vlen(dp)*cf
         do i=1,3
            dipo(i) = (dp(i) + dpn(i))
         end do
         ihsdp = 1
         do i=1,3
            dp(i) = (dp(i) + dpn(i))*cf
         end do
         write(iun3,'(a,4(f11.4,1x))') 'Total            ',
     &     (dp(i),i=1,3),vlen(dp)
         endif
         write(iun3,'(a)') ' '
         
      endif

      if (idodma.eq.1) then
      write(iun3,*) ' '
         if (moddma.eq.0) then
            write(iun3,*) 
     &    'stadard dma; all overlap density shifted to atoms'
         elseif (moddma.eq.1) then
            write(iun3,*) 
     &    'DMA sites consist of atoms + points halfway bonds'
         else
            write(iun3,*) 
     &    'DMA sites consist of atoms + centers of gaussian products'
         endif
         write(iun3,'(a)') ' '
         write(iun3,'(24x,a)') '     DMA        '
         write(iun3,'(24x,a)') '================'
         write(iun3,'(a)') ' '
         write(iun3,'(a,i4)') 'Number of sites ',nsites
         write(iun3,'(a)') ' '
         do i=1,nsites
            write(iun3,'(14x,a8,a,3f10.4)') ctag(i),' ',
     &                                     (car(j,i),j=1,3)
            write(iun3,'(f10.6)') qmom(1,i)
            write(iun3,'(3f10.6)') (qmom(j,i),j=2,4)
            write(iun3,'(5f10.6)') (qmom(j,i),j=5,9)
            write(iun3,'(7f10.6)') (qmom(j,i),j=10,16)
            write(iun3,'(9f10.6)') (qmom(j,i),j=17,25)
         end do
      endif

      return
      end
   
      subroutine stind(p,qat)
c THIS IS REALLY stint
      implicit double precision(a-h,o-z)
      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      parameter (mxel=100)
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /gauhlp/ d6to5(6,5),f10to7(10,7),g15to9(15,9),
     &                plmn(35),lmn(3,35)
      common /orbhlp/ mxorb,iuhf,ispd
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /pseudo / ipseud,ivale(numatm)
      common /coord / xyz(3,numatm)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      character*2 elemnt
      common /elem/   elemnt(mxel)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      dimension co1(3),co2(3),ca(35),cb(35)
      dimension p(*),qat(*)

      call denini

c
c     loop over shell pairs.
c

      do i=1,natoms
         if (ipseud.eq.1) then
            qat(i) = dble(ivale(i))
         else
            qat(i) = dble(nat(i))
         endif
      end do

      iao = 0
      do ishell = 1, nshell
        isha = shella(ishell)
        isht = shellt(ishell)
        ishc = shellc(ishell)
        ishn = shelln(ishell)
        iat = jan(ishell)
        co1(1) = gx(ishell)
        co1(2) = gy(ishell)
        co1(3) = gz(ishell)
        call qtype(isht,ishc,ibeg,iend)

        do ifun=ibeg,iend
          iao = iao + 1
          jao = 0
c        do jshell = 1,ishell
        do jshell = 1,nshell
          jsha = shella(jshell)
          jsht = shellt(jshell)
          jshc = shellc(jshell)
          jshn = shelln(jshell)
          jat = jan(jshell)
          co2(1) = gx(jshell)
          co2(2) = gy(jshell)
          co2(3) = gz(jshell)
          call qtype(jsht,jshc,jbeg,jend)

          do jfun=jbeg,jend
          jao = jao + 1
          s12 = 0.0d0
          t12 = 0.0d0

c===================================================
c         i loop over primitive gaussians

          do igauss=1,ishn
             ex1 = exx(isha+igauss-1)
             call fcij(isht,isha,isha+igauss-1,shladf(ishell),ca)

c            j loop over primitive gaussians

             do jgauss=1,jshn
                ex2 = exx(jsha+jgauss-1)
                call fcij(jsht,jsha,jsha+jgauss-1,shladf(jshell),cb)
                call psetup(ex1,co1,ex2,co2)

                call sint(lmn(1,ifun),lmn(1,jfun),s)
                s12 = s12 + ca(ifun)*cb(jfun)*s

c
c Kinetic energy integral
                call tint(lmn(1,ifun),lmn(1,jfun),t,ex2)
                t12 = t12 + ca(ifun)*cb(jfun)*t
             end do
          end do

         if (iat.eq.jat) then
             qat(iat) = qat(iat) - p((jao-1)*mxorb+iao)*s12
         else
             qat(iat) = qat(iat) - 0.5d0*p((jao-1)*mxorb+iao)*s12
             qat(jat) = qat(jat) - 0.5d0*p((jao-1)*mxorb+iao)*s12
         endif
        write(*,'(2(a,i3,a,i3,a,f16.8))') 's(',ifun,',',jfun,') = ',s12,
     &            't(',ifun,',',jfun,') = ',t12
c===================================================
          end do
          end do

        end do
      end do

c
c Nuclear attraction integrals, we have these, they are the same integrals
c We need for the electrostatic potential
c We sum over the atoms and multiply by atom number.
c
c      do i=1,natoms
c          print*,'atom ',i,' ========V=========='
c          call espot(xyz(1,i),xyz(2,i),xyz(3,i),epje,1)
c          print*,'==================='
c      end do

      return
      end
   
      subroutine qtype(itype,icons,ibeg,iend)
      implicit double precision(a-h,o-z)

      ibeg = 1
      if (itype.eq.4) then
          ibeg = 21
          iend = 35
      elseif (itype.eq.3) then
          ibeg = 11
          iend = 20
      elseif (itype.eq.2) then
          iend = 10
          if (icons.eq.2) ibeg = 5
      elseif (itype.eq.1) then
          iend = 4
          if (icons.eq.1) ibeg = 2
      elseif (itype.eq.0) then
          iend = 1
      endif

      return
      end

      subroutine tint(lmn1,lmn2,t12,alpha2)
      implicit double precision(a-h,o-z)
      dimension lmn1(3),lmn2(3)
      dimension lmnb(3)

      t12 = 0.0d0

      do i=1,3

c a2(2l2+1)<0|0>

          call sint(lmn1,lmn2,s12)
          t12 = t12 + alpha2*(2.0d0*lmn2(i)+1.0d0)*s12

c -2a2**2<0|+2>

          do j=1,3
             lmnb(j) = lmn2(j)
          end do
          lmnb(i) = lmnb(i) + 2
          call sint(lmn1,lmnb,s12)
          t12 = t12 - 2.0d0*alpha2*alpha2*s12

c -0.5l2(l2-1)<0|-2>

          if (lmn2(i).gt.1) then
             do j=1,3
                lmnb(j) = lmn2(j)
             end do
             lmnb(i) = lmnb(i) - 2
             call sint(lmn1,lmnb,s12)
             t12 = t12 - 0.5d0*lmn2(i)*(lmn2(i)-1)*s12
          endif

      end do
       
      return
      end

      subroutine dsint(lmn1,lmn2,ipinc,s12)
      implicit double precision(a-h,o-z)
      common /setp/ p(3),p1(3),p2(3),gam,pigam2,fac
      dimension lmn1(3),lmn2(3),ipinc(3)
      dimension xyzi(3)

      do i=1,3
         xyzi(i) = 0.0d0
      end do

      do j=1,3
         do i=0,lmn1(j)+lmn2(j)
            ip = i + ipinc(j)
            if (ip-(ip/2)*2.eq.0) 
     &      xyzi(j) = xyzi(j) + (fkxyz(p1,p2,lmn1(j),lmn2(j),i,j)
     &                *ifac2(ip-1)*pigam2/(2.0d0*gam)**(dble(ip)/2.0d0))
         end do
      end do

      s12 = fac*xyzi(1)*xyzi(2)*xyzi(3)

      return
      end

      subroutine sint(lmn1,lmn2,s12)
      implicit double precision(a-h,o-z)
      common /setp/ p(3),p1(3),p2(3),gam,pigam2,fac
      dimension lmn1(3),lmn2(3)
      dimension xyzi(3)

      do i=1,3
         xyzi(i) = 0.0d0
      end do

      do j=1,3
         do i=0,(lmn1(j)+lmn2(j))/2
            xyzi(j) = xyzi(j) + (fkxyz(p1,p2,lmn1(j),lmn2(j),2*i,j)
     &                *ifac2(2*i-1)*pigam2/(2.0d0*gam)**i)
         end do
      end do

      s12 = fac*xyzi(1)*xyzi(2)*xyzi(3)

      return
      end

      subroutine psetup(alpha1,c1,alpha2,c2)
      implicit double precision(a-h,o-z)
      common /setp/ p(3),p1(3),p2(3),gam,pigam2,fac
      dimension c1(3),c2(3)

      pi = 4.d0*datan(1.d0)
      gam = alpha1 + alpha2

      d2 = 0.0d0
      do i=1,3
         c12 = c1(i) - c2(i)
         d2 = d2 + c12*c12
         p(i) = (alpha1*c1(i) + alpha2*c2(i)) / gam
         p1(i) = p(i) - c1(i)
         p2(i) = p(i) - c2(i)
      end do
      
      fac = dexp(-alpha1*alpha2*d2/gam)

      pigam2 = dsqrt(pi/gam)

      return
      end

      double precision function fkxyz(pa,pb,l1,l2,k,iaxis)
      implicit double precision(a-h,o-z)
      dimension pa(3),pb(3)

      fkxyz = 0.0d0

      do i=0,l1
         do j=0,l2
            if (i+j.eq.k) then
               prod = fioverj(l1,i)*fioverj(l2,j)
               if (l1-i.gt.0) prod = prod*(pa(iaxis)**(l1-i))
               if (l2-j.gt.0) prod = prod*(pb(iaxis)**(l2-j))
               fkxyz = fkxyz + prod
            endif
         end do
      end do

      return
      end

      integer function ifac(n)
       
      ifac = 1
      do i=1,n
        ifac = ifac * i
      end do

      return
      end
          
      integer function ifac2(n)
       
      ifac2 = 1
      do i=1,n,2
        ifac2 = ifac2 * i
      end do

      return
      end
          
      double precision function fioverj(i,j)

      fioverj = dble(ifac(i)) / dble((ifac(j) * ifac(i-j)))

      return
      end

      double precision function drfmol(arg)
      implicit double precision(a-h,o-z)
      common /merrf/ p,a(5)
      data p/0.3275911d0/
      data a/0.254829592d0,-0.284496736d0,1.421413741d0,
     &      -1.453152027d0,1.061405429d0/
      
c      check this
      t = 1.0d0/(1.0d0 + p*arg)
      tn = t
      poly = a(1)*tn
      do i=2,5
         tn = tn*t
         poly = poly + a(i)*tn
      end do
      
      drfmol = 1.0d0 -poly*dexp(-arg*arg) 

      return
      end

      subroutine partq(q,nrank,p)
      implicit double precision (a-h, o-z)
      parameter (mxsite=300)
      parameter (lmax=4)
      parameter (lmaxf=(lmax+1)*(lmax+1))
      parameter (rtol=1.0d-4)
      character*8   ctag
      common /multip/ qmom(lmaxf,mxsite),car(3,mxsite),ctag(mxsite),
     &                nsites
      dimension q(*), p(3), r2(mxsite), iquiv(mxsite)

      lfun = (nrank+1)*(nrank+1)
c
c Partition and shift overlap multipoles to nearest center(s)
c
      r2min = 100000.0d0
      ihit = 0

      do i=1,nsites
         r2(i) = (p(1)-car(1,i))**2 + (p(2)-car(2,i))**2 +
     &           (p(3)-car(3,i))**2
         if (r2(i).lt.r2min) then
             r2min = r2(i)
             ihit = i
         endif
      end do

      if (ihit.ne.0) then
         
c Check for other centers with the same distance to P

         nequiv = 1
         iquiv(nequiv) = ihit
         do i=1,nsites
            if (i.ne.ihit.and.dabs(r2(i)-r2min).lt.rtol) then
                nequiv = nequiv + 1
                iquiv(nequiv) = i
            endif
         end do
         
c Partition multipoles over equivalent centers

         if (nequiv.gt.1) then
            do i=1,lfun
                q(i) = q(i)/dble(nequiv)
            end do
         endif

c Shift multipoles to equivalent centers

         do i=1,nequiv
            k = iquiv(i)
            call shiftq(q,0,nrank,p(1)-car(1,k),p(2)-car(2,k),
     &            p(3)-car(3,k),qmom(1,k),lmax)
         end do
      endif

      return
      end

      subroutine shiftq(q,lomin,lomax,x,y,z,qn,lnmax)
      implicit double precision (a-h, o-z)
      parameter (lmax=4)
      parameter (lmaxf=(lmax+1)*(lmax+1))
 
c
c shift multipoles q of ranks lmin to lmax at P 
c (center of overlap gaussian) to new point S (qn)
c
c          l  k  -               - (1/2)
c          -- -- | (l+m)   (l-m) |
c Qlm(S) = \  \  | |   | * |   | |  Qkq(P) * Rl-k,m-1(S-P)
c          /  /  | (k+q)   (k-q) |
c          -- -- -               -
c         k=0 q=-k
c
c x,y,z = (S-P)
c
c for k.ne.0 the Qlk (or Rlk) are complex:
c
c Ql, k = 1/sqrt(2)(-1)**k {Qlkc + i Qlks}
c Ql,-k = 1/sqrt(2) {Qlkc - i Qlks}
c
 
      dimension q(*), qn(*), r(lmaxf)

      call solidh(x,y,z,lnmax,r)

c  ln is rank of target multipole

      do ln = lomin,lnmax
         kn = ln*ln+1
         mn = 0

c  l is rank of source multipoles

100      do l = lomin,min(ln,lomax)
            lr = ln - l
            l2 = l*l
            lr2 = lr*lr
            do m = max(-l,mn-lr),min(l,mn+lr)
                if (m.eq.0) then
                  qc =  q(l2+1)
                  qs =  0.0d0
                else
                  qc =  q(l2+2*iabs(m))
                  qs =  q(l2+2*iabs(m)+1)
                  if (m.lt.0) qs = -qs
                endif

                mr = mn - m

                if (mr.eq.0) then
                  rc =  r(lr2+1)
                  rs =  0.0d0
                else
                  rc =  r(lr2+2*iabs(mr))
                  rs =  r(lr2+2*iabs(mr)+1)
                  if (mr.lt.0) rs = -rs
                endif

                fac = 0.5d0*qrsgn(mn)/(qrsgn(m)*qrsgn(mr))
                root = dsqrt(fioverj(ln-mn,l-m)*fioverj(ln+mn,l+m))

                qn(kn) = qn(kn) + root*fac*(qc*rc-qs*rs)
                if (mn.gt.0)
     &              qn(kn+1) = qn(kn+1) + root*fac*(qc*rs+qs*rc)
            end do
         end do
 
         mn = mn + 1
         kn = ln*ln+2*mn
         if (mn.le.ln) goto 100
      end do

      return
      end

      double precision function qrsgn(m)
      implicit double precision (a-h, o-z)
      parameter (hrt2=0.7071067811865475244d0)

      if (m.eq.0) then
         qrsgn = 0.5d0
      elseif (m.lt.0) then
         qrsgn = hrt2
      else
         if (mod(m,2).eq.0) then
            qrsgn = hrt2
         else
            qrsgn = -hrt2
         endif
      endif

      return
      end

      subroutine solidh (x,y,z,j,r)
      implicit double precision (a-h, o-z)
 
c  computes regular solid harmonics r**k ckq(theta,phi) for ranks k up
c  to j, if j > =  0;
c  or irregular solid harmonics r**(-k-1) ckq(theta,phi) for ranks k up
c  to |j|, if j < 0.
 
      dimension r(*)
 
c  locations in r are used as follows:
 
c        1    2    3    4    5    6    7    8    9   10   11  ...
c  kq  =  00   10   11c  11s  20   21c  21s  22c  22s  30   31c ...
 
c  r(k,0) is real and is left in location k**2 + 1.
c  r(k,mc) and r(k,ms) are sqrt(2) times the real and imaginary parts
c  respectively of the complex solid harmonic r(k,-m)*  =  (-1)**m r(k,m),
c  and are left in locations k**2 + 2m and k**2 + 2m + 1 respectively.
c
c  S.l. Price et al, Molecular Physics, 1984, 52, 987: Table 1
c
c  R00  = 1
c 
c  R10  = z
c  R11c = x
c  R11s = y
c
c  R20  = 0.5(3z**2-r**2)
c  R21c = dsqrt(3)xz
c  R21s = dsqrt(3)yz
c  R22c = 0.5*dsqrt(3)(x**2-y**2)
c  R22s = dsqrt(3)xy
c
c  etc ..
c

      l = abs(j)

      rr = x*x + y*y + z*z
 
      if (j.lt.0) then
c  irregular
        rr = 1.0/rr
        rfx = x*rr
        rfy = y*rr
        rfz = z*rr
        r(1) = dsqrt(rr)
        r(2) = rfz*r(1)
        r(3) = rfx*r(1)
        r(4) = rfy*r(1)
 
      else
c  regular
        r(1) = 1.0
        r(2) = z
        r(3) = x
        r(4) = y
        rfz = z
        rfx = x
        rfy = y
      endif
 
c  remaining values are found using recursion formulae, relating
c  the new set n to the current set k and the previous set p.
 
      k = 1
10    n  = k + 1
      ln = n*n + 1
      lk = k*k + 1
      lp = (k-1)**2 + 1
      a2kp1 = k+k + 1
 
c  obtain r(k+1,0) from r(k,0)*r(1,0) and r(k-1,0)
 
      r(ln) = (a2kp1*r(lk)*rfz-k*rr*r(lp))/(k+1)
 
      m = 1
      ln = ln + 1
      lk = lk + 1
      lp = lp + 1
      if (k.eq.1) goto 40
 
c  obtain r(k+1,m) from r(k,m)*r(1,0) and r(k-1,m)
 
20    r(ln) = (a2kp1*r(lk)*rfz-dsqrt(dble((k+m)*(k-m)))*rr*r(lp))
     &                /dsqrt(dble((n+m)*(n-m)))
      r(ln+1) = (a2kp1*r(lk+1)*rfz-dsqrt(dble((k+m)*(k-m)))*rr*r(lp+1))
     &                /dsqrt(dble((n+m)*(n-m)))
      m  = m + 1
      ln = ln + 2
      lk = lk + 2
      lp = lp + 2
      if (m.lt.k) goto 20
 
c  obtain r(k+1,k) from r(k,k)*r(1,0)
 
40    r(ln) = dsqrt(dble(n+k))*r(lk)*rfz
      r(ln+1) = dsqrt(dble(n+k))*r(lk+1)*rfz
      ln = ln + 2
 
c  obtain r(k+1,k+1) from r(k,k)*r(1,1)
 
      s = dsqrt(dble(n+k))/dsqrt(dble(n+n))
      r(ln) = s*(rfx*r(lk)-rfy*r(lk+1))
      r(ln+1) = s*(rfx*r(lk+1)+rfy*r(lk))
 
      k = k + 1
      if (k.lt.l) goto 10
      return
 
      end

      subroutine caldid(q,coo)
      implicit double precision (a-h, o-z)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /athlp/ iatoms, mxnat
      dimension q(*),coo(3,*)

      do i=1,3
         dipo(i) = 0.0d0
      end do

      do i=1,iatoms
         dipo(1) = dipo(1) + coo(1,i)*q(i)
         dipo(2) = dipo(2) + coo(2,i)*q(i)
         dipo(3) = dipo(3) + coo(3,i)*q(i)
      end do
 
      ihsdp = 2

      return
 
      end

      subroutine prtdip
      implicit double precision (a-h, o-z)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      cf = 5.2917715d-11*1.601917d-19/3.33564d-30

      write(iun3,'(a)') '   '
      if (ihsdp.eq.1) then
         write(iun3,'(a)') 
     &         'Dipole moment calculated from Wavefunction:'
      elseif (ihsdp.eq.2) then
         write(iun3,'(a)') 
     &         'Dipole moment calculated from partial charges:'
      elseif (ihsdp.eq.3) then
         write(iun3,'(a)') 
     &         'Dipole moment evaluated  from Monopoles and Dipoles:'
      endif
      write(iun3,'(a)') '   '

      write(iun3,'(a)') ' '
      write(iun3,'(24x,a)') 'Dipole moment (Debye)'
      write(iun3,'(24x,a)') '======================='
      write(iun3,'(17x,a)') 
     &         '      X           Y           Z          Scalar'
      write(iun3,'(a)') ' '
      write(iun3,'(a,4(f11.4,1x))') '                 ',
     &    (dipo(i)*cf,i=1,3),vlen(dipo)*cf

      return
      end
