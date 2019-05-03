      subroutine gtipar(istat,iat,iti,trsi1,trsi2,ityp)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxamb=1590)
      parameter (mxitor=42)
      integer*2 ityp
      integer ambcls
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxcls=50)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      common /itrpar/ tori1(2,mxitor),tori2(2,mxitor),
     &                itrcon(4,mxitor)
      parameter (mxgitr=42)
      integer gficon
      common /gfipar/ gftri1(2,mxgitr),gftri2(2,mxgitr),
     &                gficon(4,mxgitr)

      logical doamb, dogaff
      save isq
      dimension iat(4),icls(4),iti(4),trsi1(4),trsi2(4),isq(3,6),
     &          it(4),itb(4),ityp(*)
      data isq / 1,2,4, 1,4,2, 2,1,4, 4,1,2, 2,4,1, 4,2,1/

      todeg = 45.0e0 / atan(1.0e0)
      wght6 = dble(1)/dble(6)
      doamb = .false.
      dogaff = .false.
      istat = 0

      i1 = int(ityp(iat(1)))
      i2 = int(ityp(iat(2)))
      i3 = int(ityp(iat(3)))
      i4 = int(ityp(iat(4)))

      if (i1.gt.0 .and. i2.gt.0 .and. i3.gt.0 .and. i4.gt.0) then
         doamb = .true.
      elseif (i1.lt.0 .or. i2.lt.0 .or. i3.lt.0 .or. i4.lt.0) then
         dogaff = .true.
      endif

      if (doamb) then

            do j=1,4
               icls(j) = ambcls(ityp(iat(j)))
            end do
               
            do j = 1, mxitor

               if (icls(3).eq.itrcon(3,j)) then

                  do k=1,6

                     if (icls(isq(1,k)).eq.itrcon(1,j).and.
     &                   icls(isq(2,k)).eq.itrcon(2,j).and.
     &                   icls(isq(3,k)).eq.itrcon(4,j)) then


                           istat = 1
                           iti(1) = iat(isq(1,k))
                           iti(2) = iat(isq(2,k))
                           iti(3) = iat(3)
                           iti(4) = iat(isq(3,k))
                 
                           trsi1(1) = tori1(1,j)
                           trsi1(2) = tori1(2,j)
                           angl = tori1(2,j) / todeg
                           trsi1(3) = cos(angl)
                           trsi1(4) = sin(angl)

                           trsi2(1) = tori2(1,j)
                           trsi2(2) = tori2(2,j)
                           angl = tori2(2,j) / todeg
                           trsi2(3) = cos(angl)
                           trsi2(4) = sin(angl)

                     end if

                  end do

               end if

            end do

      endif

      if (dogaff) then

          do j=1,4
             if (ityp(iat(j)).gt.0) then
                itb(j) = mapagf(ambcls(ityp(iat(j))))
             else
                itb(j) = iabs(int(ityp(iat(j))))
             endif
          end do

          do j=1,4
             it(j) = itb(j)
          end do

          call gtgfi(istat1,isq,iat,it,iti,trsi1,trsi2,ityp)
          if (istat1.ne.0) then
              istat = 1
              return
          endif

c try one X, in all permutations

          do j=1,4
             it(j) = itb(j)
          end do
          it(1) = 1
          call gtgfi(istat1,isq,iat,it,iti,trsi1,trsi2,ityp)
          if (istat1.ne.0) then
              istat = 1
              return
          endif

          do j=1,4
             it(j) = itb(j)
          end do
          it(2) = 1
          call gtgfi(istat1,isq,iat,it,iti,trsi1,trsi2,ityp)
          if (istat1.ne.0) then
              istat = 1
              return
          endif

          do j=1,4
             it(j) = itb(j)
          end do
          it(4) = 1
          call gtgfi(istat1,isq,iat,it,iti,trsi1,trsi2,ityp)
          if (istat1.ne.0) then
              istat = 1
              return
          endif

c try two X's, in all permutations

          do j=1,4
             it(j) = itb(j)
          end do
          it(1) = 1
          it(2) = 1
          call gtgfi(istat1,isq,iat,it,iti,trsi1,trsi2,ityp)
          if (istat1.ne.0) then
              istat = 1
              return
          endif

          do j=1,4
             it(j) = itb(j)
          end do
          it(1) = 1
          it(4) = 1
          call gtgfi(istat1,isq,iat,it,iti,trsi1,trsi2,ityp)
          if (istat1.ne.0) then
              istat = 1
              return
          endif

          do j=1,4
             it(j) = itb(j)
          end do
          it(2) = 1
          it(4) = 1
          call gtgfi(istat1,isq,iat,it,iti,trsi1,trsi2,ityp)
          if (istat1.ne.0) then
              istat = 1
              return
          endif

c try three X's

          do j=1,4
             it(j) = itb(j)
          end do

          it(1) = 1
          it(2) = 1
          it(4) = 1

          call gtgfi(istat1,isq,iat,it,iti,trsi1,trsi2,ityp)
          if (istat1.ne.0) then
              istat = 1
              return
          endif

      endif

      return
      end

      subroutine gtgfi(istat,isq,iat,it,iti,trsi1,trsi2,ityp)
      implicit real (a-h,o-z), integer (i-n)
      integer*2 ityp
      parameter (mxgitr=42)
      integer gficon
      common /gfipar/ gftri1(2,mxgitr),gftri2(2,mxgitr),
     &                gficon(4,mxgitr)

      dimension iat(4),icls(4),iti(4),trsi1(4),trsi2(4),isq(3,6),
     &          it(4),ityp(*)

      todeg = 45.0e0 / atan(1.0e0)
      istat = 0

      do j = 1, mxgitr

         if (it(3).eq.gficon(3,j)) then

            do k=1,6

               if (it(isq(1,k)).eq.gficon(1,j).and.
     &             it(isq(2,k)).eq.gficon(2,j).and.
     &             it(isq(3,k)).eq.gficon(4,j)) then


                     istat = 1
                     iti(1) = iat(isq(1,k))
                     iti(2) = iat(isq(2,k))
                     iti(3) = iat(3)
                     iti(4) = iat(isq(3,k))
                 
                     trsi1(1) = gftri1(1,j) 
                     trsi1(2) = gftri1(2,j)
                     angl = gftri1(2,j) / todeg
                     trsi1(3) = cos(angl)
                     trsi1(4) = sin(angl)

                     trsi2(1) = gftri2(1,j)
                     trsi2(2) = gftri2(2,j)
                     angl = gftri2(2,j) / todeg
                     trsi2(3) = cos(angl)
                     trsi2(4) = sin(angl)

               end if

            end do

         end if

      end do

      return
      end

      subroutine itrard(istat,mxitrs,
     &                  nti,iti,trsi1,trsi2,
     &                  iconn,ityp,iopt)

      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)

      common /athlp/  iatoms, mxnat
      integer*2 ityp
      dimension iconn(mxcon+1,*),ityp(*),iopt(*),
     &          iti(4,*),trsi1(4,*),trsi2(4,*),
     &          iat(4),icnn(mxcon)
      dimension itit(4),trsi1t(4),trsi2t(4)
      
      istat = 1
      nti = 0

      do i=1,iatoms

         ncnn = 0
         do j=1,iconn(1,i)
            jj = iconn(1+j,i)
            if (jj.gt.0) then
               ncnn = ncnn + 1
               icnn(ncnn) = jj
            endif
         end do
         
         if (ncnn.eq.3) then

            iat(1)  = icnn(1)
            iat(2)  = icnn(2)
            iat(3)  = i
            iat(4)  = icnn(3)

            ido = 0
            do j=1,4
               if (iopt(iabs(iat(j))).eq.1) ido = 1
            end do

            call gtipar(istat,iat,itit,trsi1t,trsi2t,ityp)

            if (istat.eq.1.and.ido.eq.1) then

               nti = nti + 1
               if (nti.lt.mxitrs) then
                  do j=1,4
                     iti(j,nti) = itit(j)
                     trsi1(j,nti) = trsi1t(j)
                     trsi2(j,nti) = trsi2t(j)
                  end do
               else

                  istat = 0
                  return

               endif

            endif

         end if

      end do

      return
      end 

      subroutine gttpar(istat,iat,it,trs1,trs2,trs3,trs4,ityp)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      integer*2 ityp
      integer ambcls
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxcls=50)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      parameter (mxtor=351)
      integer torcon
      common /torpar/ tor1(2,mxtor),tor2(2,mxtor),
     &                tor3(2,mxtor),tor4(2,mxtor),
     &                torcon(4,mxtor)
      parameter (mxgtor=638)
      integer gftcon
      common /gftpar/ gftor1(2,mxgtor),gftor2(2,mxgtor),
     &                gftor3(2,mxgtor),gftor4(2,mxgtor),
     &                gftcon(4,mxgtor)
      integer amb2gf
      common /gftyp/  gfvdw(2,mxgff),amb2gf(mxamb)
      character*3 ambstr
      character*2 gffstr
      common /ffstr/  ambstr(mxamb), gffstr(mxgff)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5

      logical doamb, dogaff, first
      dimension iat(4),it(4),trs1(4),trs2(4),trs3(4),trs4(4),
     &          ityp(*)

      todeg = 45.0e0 / atan(1.0e0)
      doamb = .false.
      dogaff = .false.
      first = .true.
      istat = 0

      i1 = int(ityp(iat(1)))
      i2 = int(ityp(iat(2)))
      i3 = int(ityp(iat(3)))
      i4 = int(ityp(iat(4)))

      if (i1.gt.0 .and. i2.gt.0 .and. i3.gt.0 .and. i4.gt.0) then
         doamb = .true.
      elseif (i1.lt.0 .or. i2.lt.0 .or. i3.lt.0 .or. i4.lt.0) then
         dogaff = .true.
      endif

      it(1) = iat(1)
      it(2) = iat(2)
      it(3) = iat(3)
      it(4) = iat(4)

      if (doamb) then
          ita = ambcls(i1)
          itb = ambcls(i2)
          itc = ambcls(i3)
          itd = ambcls(i4)
      endif
      
      if (dogaff) then

         if (i1.gt.0) then
            ita = mapagf(ambcls(i1))
         else
            ita = iabs(i1)
            if (ita.eq.72) ita = 6
         endif
         if (i2.gt.0) then
            itb = mapagf(ambcls(i2))
         else
            itb = iabs(i2)
            if (itb.eq.72) itb = 6
         endif
         if (i3.gt.0) then
            itc = mapagf(ambcls(i3))
         else
            itc = iabs(i3)
            if (itc.eq.72) itc = 6
         endif
         if (i4.gt.0) then
            itd = mapagf(ambcls(i4))
         else
            itd = iabs(i4)
            if (itd.eq.72) itd = 6
         endif
      endif

      if (itb.lt.itc) then
         it1 = ita
         it2 = itb
         it3 = itc
         it4 = itd
      elseif (itc.lt.itb) then
         it1 = itd
         it2 = itc
         it3 = itb
         it4 = ita
      elseif (ita.le.itd) then
         it1 = ita
         it2 = itb
         it3 = itc
         it4 = itd
      elseif (itd.lt.ita) then
         it1 = itd
         it2 = itc
         it3 = itb
         it4 = ita
      endif

      if (doamb) then

10        continue

          do l=1,mxtor
             if (torcon(1,l).eq.it1.and.
     &           torcon(2,l).eq.it2.and.
     &           torcon(3,l).eq.it3.and.
     &           torcon(4,l).eq.it4) then

                 istat = 1
                 trs1(1) = tor1(1,l)
                 trs1(2) = tor1(2,l)
                 angl = tor1(2,l) / todeg
                 trs1(3) = cos(angl)
                 trs1(4) = sin(angl)

                 trs2(1) = tor2(1,l)
                 trs2(2) = tor2(2,l)
                 angl = tor2(2,l) / todeg
                 trs2(3) = cos(angl)
                 trs2(4) = sin(angl)

                 trs3(1) = tor3(1,l)
                 trs3(2) = tor3(2,l)
                 angl = tor3(2,l) / todeg
                 trs3(3) = cos(angl)
                 trs3(4) = sin(angl)

                 trs4(1) = tor4(1,l)
                 trs4(2) = tor4(2,l)
                 angl = tor4(2,l) / todeg
                 trs4(3) = cos(angl)
                 trs4(4) = sin(angl)

             endif
          end do

          if (istat.eq.0.and.it2.eq.it3.and.first) then
             itmp = it1
             it1 = it4
             it4 = itmp
             first = .false.
             goto 10
          endif

          if (istat.eq.0.and.(i1.ge.1254.and.i2.ge.1254)
     &        .and.(i3.ge.1254.and.i3.ge.1254)) then
             ita = amb2gf(i1)
             itb = amb2gf(i2)
             itc = amb2gf(i3)
             itd = amb2gf(i4)

             if (itb.lt.itc) then
                it1 = ita
                it2 = itb
                it3 = itc
                it4 = itd
             elseif (itc.lt.itb) then
                it1 = itd
                it2 = itc
                it3 = itb
                it4 = ita
             elseif (ita.le.itd) then
                it1 = ita
                it2 = itb
                it3 = itc
                it4 = itd
             elseif (itd.lt.ita) then
                it1 = itd
                it2 = itc
                it3 = itb
                it4 = ita
             endif

             first = .true.
             dogaff = .true.
          endif

      endif

      if (dogaff) then

100       continue

          do l=1,mxgtor
             if (gftcon(1,l).eq.it1.and.
     &           gftcon(2,l).eq.it2.and.
     &           gftcon(3,l).eq.it3.and.
     &           gftcon(4,l).eq.it4) then

                 istat = 1
                 trs1(1) = gftor1(1,l)
                 trs1(2) = gftor1(2,l)
                 angl = gftor1(2,l) / todeg
                 trs1(3) = cos(angl)
                 trs1(4) = sin(angl)

                 trs2(1) = gftor2(1,l)
                 trs2(2) = gftor2(2,l)
                 angl = gftor2(2,l) / todeg
                 trs2(3) = cos(angl)
                 trs2(4) = sin(angl)

                 trs3(1) = gftor3(1,l)
                 trs3(2) = gftor3(2,l)
                 angl = gftor3(2,l) / todeg
                 trs3(3) = cos(angl)
                 trs3(4) = sin(angl)

                 trs4(1) = gftor4(1,l)
                 trs4(2) = gftor4(2,l)
                 angl = gftor4(2,l) / todeg
                 trs4(3) = cos(angl)
                 trs4(4) = sin(angl)

             endif
          end do

          if (istat.eq.0.and.first) then
             it1 = 1
             it4 = 1
             first = .false.
             goto 100
          endif
      endif

      if (istat.eq.0) then
         if (doamb) then
            write(iun5,'(a,4(i4,1x),a,4(a3,1x),a,4(i4,1x))') 'atom ',
     &         it(1),it(2),it(3),it(4),' type ',
     &         ambstr(i1),ambstr(i2),ambstr(i3),ambstr(i4),
     &         ' class ',ita,itb,itc,itd
         endif
         if (dogaff) then
            write(iun5,'(4(i4,1x),4(a2,1x))') 
     &         it(1),it(2),it(3),it(4),
     &         gffstr(ita),gffstr(itb),gffstr(itc),gffstr(itd)
         endif
      endif

      return
      end

      subroutine torard(istat,maxtors,
     &                  nbnd,ibnd,
     &                  nt,it,trs1,trs2,trs3,trs4,
     &                  iconn,ityp,iopt)

      implicit real (a-h,o-z), integer (i-n)

      parameter (mxcon=10)
      parameter (mxamb=1590)
      integer ambcls
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxtor=351)
      integer torcon
      common /torpar/ tor1(2,mxtor),tor2(2,mxtor),
     &                tor3(2,mxtor),tor4(2,mxtor),
     &                torcon(4,mxtor)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      integer*2 ityp
      dimension ibnd(2,*),iconn(mxcon+1,*),ityp(*),iopt(*),
     &          it(4,*),trs1(4,*),trs2(4,*),trs3(4,*),trs4(4,*),
     &          itt(4),trs1t(4),trs2t(4),trs3t(4),trs4t(4),iat(4)

      todeg = 45.0e0 / atan(1.0e0)
      istat = 1
      nt = 0

      do i=1,nbnd
         iat(2) = ibnd(1,i)
         iat(3) = ibnd(2,i)

         do j=1,iconn(1,iat(2))
            iat(1) = iconn(1+j,iat(2))
            if (iat(1).gt.0.and.iat(1).ne.iat(3)) then

               do k=1,iconn(1,iat(3))
                  iat(4) = iconn(1+k,iat(3))

                  if (iat(4).gt.0.and.
     &               (iat(4).ne.iat(2).and.iat(4).ne.iat(1))) then

                     ido = 0
                     do l=1,4
                        if (iopt(iabs(iat(l))).eq.1) ido = 1
                     end do

                     call gttpar(istat,
     &                           iat,itt,trs1t,trs2t,trs3t,trs4t,ityp)

                     if (istat.eq.1.and.ido.eq.1) then
                        nt = nt + 1
                        if (nt.lt.maxtors) then
                           do l=1,4
                              it(l,nt)   = itt(l)
                              trs1(l,nt) = trs1t(l)
                              trs2(l,nt) = trs2t(l)
                              trs3(l,nt) = trs3t(l)
                              trs4(l,nt) = trs4t(l)
                           end do

                        else

                           istat = 0
                           return

                        endif
                     else
                        write(iun5,*) 'torsion parameter not found'
                     endif

                  endif

               end do

            endif

         end do

      end do

      return
      end

      subroutine prttod(istat,maxtors,
     &                  nt,it,trs1,trs2,trs3,trs4,
     &                  iconn,ityp,iopt)

      implicit real (a-h,o-z), integer (i-n)

      parameter (mxcon=10)
      parameter (mxamb=1590)
      integer ambcls
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      integer*2 ityp
      dimension iconn(mxcon+1,*),ityp(*),iopt(*),
     &          it(4,*),trs1(4,*),trs2(4,*),trs3(4,*),trs4(4,*)

      print*,''
      print*,'Torsion parameters:'
      print*,''

      do i=1,nt
         print*,'it ',(it(l,i),l=1,4)
         print*,' trs1 ',(trs1(l,i),l=1,4)
         print*,' trs2 ',(trs2(l,i),l=1,4)
         print*,' trs3 ',(trs3(l,i),l=1,4)
         print*,' trs4 ',(trs4(l,i),l=1,4)
      end do

      print*,''

      return
      end

      subroutine pritod(istat,mxitrs,
     &                  nti,iti,trsi1,trsi2,
     &                  iconn,ityp,iopt)

      implicit real (a-h,o-z), integer (i-n)

      parameter (mxcon=10)
      parameter (mxamb=1590)
      integer ambcls
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      integer*2 ityp
      dimension iconn(mxcon+1,*),ityp(*),iopt(*),
     &          iti(4,*),trsi1(4,*),trsi2(4,*)

      print*,''
      print*,'Impromper Torsion parameters:'
      print*,''

      do i=1,nti
         print*,'iti ',(iti(l,i),l=1,4)
         print*,' trsi1 ',(trsi1(l,i),l=1,4)
         print*,' trsi2 ',(trsi2(l,i),l=1,4)
      end do

      print*,''

      return
      end

      subroutine tors(et,idoimp,
     &                nt,it,trs1,trs2,trs3,trs4,
     &                nti,iti,trsi1,trsi2,
     &                coo,forces)

      implicit real (a-h,o-z), integer (i-n)

      common /athlp/  iatoms, mxnat

      dimension cosna(4),sinna(4),v(4),sa(4),ca(4),phi(4),dphi(4)
      dimension va(3),vb(3),vc(3),vd(3),vba(3),vcb(3),vdc(3)
      dimension vca(3),vdb(3),vt(3),vu(3),vtu(3),dedt(3),dedu(3)
      dimension deda(3),dedb(3),dedc(3),dedd(3),vtmp1(3),vtmp2(3)

      dimension coo(3,*),forces(3,*)
      dimension iti(4,*),trsi1(4,*),trsi2(4,*)
      dimension it(4,*),trs1(4,*),trs2(4,*),trs3(4,*),trs4(4,*)

      todeg = 45.0e0 / atan(1.0e0)
      et = 0.0e0

      if (idoimp.eq.1) then
         ntor = nti
         nterms = 2
      else
         ntor = nt
         nterms = 4
      endif

      do i =1,ntor
            if (idoimp.eq.1) then
               ia = iti(1,i)
               ib = iti(2,i)
               ic = iti(3,i)
               id = iti(4,i)
            else
               ia = it(1,i)
               ib = it(2,i)
               ic = it(3,i)
               id = it(4,i)
            endif

c coordinates of atoms a,b,c,d
            do j=1,3
              va(j) = coo(j,ia)
              vb(j) = coo(j,ib)
              vc(j) = coo(j,ic)
              vd(j) = coo(j,id)
            end do

c vectors ba, cb, dc, ca, db

            call vsub(vb,va,vba,3)
            call vsub(vc,vb,vcb,3)
            call vsub(vd,vc,vdc,3)
            call vsub(vc,va,vca,3)
            call vsub(vd,vb,vdb,3)

c t outvector ba x cb, u outvector cb x dc

            call crprod(vba,vcb,vt)
            call crprod(vcb,vdc,vu)

c tu outvector t x u

            call crprod(vt,vu,vtu)


            rt2 = vt(1)*vt(1) + vt(2)*vt(2) + vt(3)*vt(3)
            ru2 = vu(1)*vu(1) + vu(2)*vu(2) + vu(3)*vu(3)

            rtru = vlen(vu)*vlen(vt)

            if (rtru .ne. 0.0e0) then

               rcb = vlen(vcb)

c improduct t.u /|t||u| , cosine of angle  between t and u

               call impsc(vt,vu,cosna(1))

c improduct cb.tu /|cb||t||u|, this is mainly for the sign

               call timpsc(vcb,vtu,sinna(1))
               sinna(1) = sinna(1) / (rcb*rtru)

c torsional parameters 

               if (idoimp.eq.1) then

                  v (1) = trsi1(1,i)
                  ca(1) = trsi1(3,i)
                  sa(1) = trsi1(4,i)

                  v (2) = trsi2(1,i)
                  ca(2) = trsi2(3,i)
                  sa(2) = trsi2(4,i)

               else

                  v (1) = trs1(1,i)
                  ca(1) = trs1(3,i)
                  sa(1) = trs1(4,i)

                  v (2) = trs2(1,i)
                  ca(2) = trs2(3,i)
                  sa(2) = trs2(4,i)

                  v (3) = trs3(1,i)
                  ca(3) = trs3(3,i)
                  sa(3) = trs3(4,i)

                  v (4) = trs4(1,i)
                  ca(4) = trs4(3,i)
                  sa(4) = trs4(4,i)

               endif

c knowing    cos(1*angle) and sin(1*angle)
c calculate  cos(2*angle) and sin(2*angle) etc.
c
c cos(2*angle) = cos(angle+angle) = 
c                cos(angle)*cos(angle) - sin(angle)*sin(angle)
c
c cos(2*angle) = cos(angle)*cos(angle) - sin(angle)*sin(angle)
c sin(2*angle) = cos(angle)*sin(angle) + cos(angle)*sin(angle)

               do j=1,nterms-1
                  cosna(j+1) = cosna(1)*cosna(j) - sinna(1)*sinna(j)
                  sinna(j+1) = cosna(1)*sinna(j) + sinna(1)*cosna(j)
               end do

c Etors = sum n=1,nterms Vn*[1 + cos(n*angle-SHIFTANGLEn]
c parameters Vn and SHIFTANGLEn

c cos(alpha-beta) = cos(alpha)*cos(beta) + sin(alpha)*sin(beta)
c alpha=n*angle, beta = SHIFTANGLEn, 
c cos(SHIFTANGLEn) and sin(SHIFTANGLEn) stored in ca(n) and sa(n)

               do j=1,nterms
                  phi(j)  = 1.0e0 + (cosna(j)*ca(j) + sinna(j)*sa(j))
               end do

               do j=1,nterms
                  dphi(j) = dble(j)*(cosna(j)*sa(j) - sinna(j)*ca(j))
               end do

               e      = 0.0e0
               detor  = 0.0e0

               do j=1,nterms
                   e     = e     + v(j)*phi(j)
                   detor = detor + v(j)*dphi(j)
               end do

c detor -> u and t: du, dt
c detor * (t x cb) / |t| |cb|

               call crprod(vt,vcb,dedt)
               s = detor / (rt2*rcb)
               call vscal(dedt,3,s)

c detor * (u x cb) / |u| |cb|

               call crprod(vu,vcb,dedu)
               s = -detor / (ru2*rcb)
               call vscal(dedu,3,s)

c outvector cb x dt


               call crprod(vcb,dedt,deda)
               call vscal(deda,3,-1.0e0)

c outvector ca x dt + outvector dc x du

               call crprod(vca,dedt,vtmp1)
               call crprod(vdc,dedu,vtmp2)
               call vsub(vtmp1,vtmp2,dedb,3)

c outvector ba x dt + outvector db x du

               call crprod(vba,dedt,vtmp1)
               call crprod(vdb,dedu,vtmp2)
               call vsub(vtmp2,vtmp1,dedc,3)

c outvector cb x du

               call crprod(vcb,dedu,dedd)
               call vscal(dedd,3,-1.0e0)

c               if (idoimp.eq.0) 
c     &         write(6,'(i4,a,4(i3,a),f9.3,a,12f9.3)') 
c     &               i,' ',ia,' ',ib,' ',ic,' ',id,' ',e,' ded ',
c     &               (deda(j),j=1,3),(dedb(j),j=1,3),
c     &               (dedc(j),j=1,3),(dedd(j),j=1,3)
               

               et = et + e

               do j=1,3
                  forces(j,ia) = forces(j,ia) + deda(j)
                  forces(j,ib) = forces(j,ib) + dedb(j)
                  forces(j,ic) = forces(j,ic) + dedc(j)
                  forces(j,id) = forces(j,id) + dedd(j)
               end do

            end if
      end do

      return
      end
