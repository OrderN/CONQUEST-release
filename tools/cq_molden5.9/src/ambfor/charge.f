      subroutine asschd(idebug,ityp,iconn,q)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxamb=1590)
      parameter (mxcon=10)
      integer ambcls
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      common /athlp/  iatoms, mxnat
      logical opfil,gargpl,osingl,dolbfgs,oqscal
      common /opts/ ideb,imon,iarc,ilog,iout,osingl,dolbfgs,oqscal
      integer*2 ityp
      dimension q(*),ityp(*),iconn(mxcon+1,*)


      do i=1,iatoms
          if (ityp(i).gt.0.and.ityp(i).le.mxamb) then
             q(i) = ambchg(ityp(i))
          elseif (ityp(i).eq.0) then
             q(i) = 0.0e0
          endif
      end do

      call ckncys(q,ityp,iconn)

      call discq(q)

      call sumq(idebug,q)

      if (oqscal) call qscal(q)

      return
      end

      subroutine ckncys(q,ityp,iconn)
c  check for negative cysteine
      implicit real (a-h,o-z), integer ( i-n)
      parameter (numres=50000)
      parameter (mxcon=10)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      common /residu/  qtot,ihsres,
     &                 ires(numres),ibeg(numres),iend(numres)
      common /cyschg/ cysneg(9)
      integer*2 ityp
      dimension q(*),ityp(*),iconn(mxcon+1,*)

      do i=1,ihsres
         ineg = 0
         do j=ibeg(i),iend(i)
            if (ityp(j).eq.85) then
               if (iconn(1,j).eq.1) then
                  if (ityp(iconn(2,j)).eq.83) ineg = 1
               endif
            endif
         end do

         if (ineg.eq.1) then
            do j=ibeg(i),iend(i)
               jj = ityp(j) - 76
               if (jj.ge.1.and.jj.le.9) then
                  q(j) = cysneg(jj)
               endif
            end do
         endif

      end do

      return
      end

      subroutine sumq(idebug,q)
      implicit real (a-h,o-z), integer ( i-n)
      parameter (numres=50000)
      parameter (mxion=2000)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      common /residu/  qtot,ihsres,
     &                 ires(numres),ibeg(numres),iend(numres)
      common /h2oer/   numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      character*5 strion
      common /ions/    strion(2)
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      dimension q(*)
      data strion /' Cl- ',' Na+ '/

      qtot = 0.0e0
      write(iun5,'(a)') ' '
      do i=1,ihsres
         sum = 0.0e0
c         print*,'ibeg ',ibeg(i),' iend ',iend(i)
         do j=ibeg(i),iend(i)
            sum = sum + q(j)
c            print*,'   ',q(j)
         end do
         qtot = qtot + sum
         if (idebug.eq.1) 
     &      write(iun5,'(a,i5,a,f8.3)') 'res ',i,' tot charge ',sum
      end do

      nion = abs(nint(qtot))
      if (qtot.lt.0.0e0) iontyp = 2

      if (nion.ne.0) then
         write(iun5,'(a,f7.3)') ' Total charge box ',qtot
         write(iun5,'(a)') ' '
         if (.not.cell) then
            write(iun5,'(a,i4,a,a)') ' Adding ',nion,strion(iontyp),
     &                     ' ions, to neutralise box'
         endif
      endif

      write(iun5,'(a)') ' '

      return
      end

      subroutine discq(q)
      implicit real (a-h,o-z), integer ( i-n)
      parameter (numres=50000)
      common /rdwr/ iun1,iun2,iun3,iun4,iun5
      common /residu/  qtot,ihsres,
     &                 ires(numres),ibeg(numres),iend(numres)
      logical corr
      dimension q(*)

      write(iun5,'(a)') ' '


      do i=1,ihsres
         n = iend(i) - ibeg(i) + 1
         corr = .false.
         if (n.eq.1) then

            qt = q(ibeg(i))
            rq = nint(qt)
            if (qt-rq.ne.0.0e0) then
               write(iun5,'(a,i5,a,f9.3)') ' Single atom ',ibeg(i),
     &                           ' has non integer charge: ',qt
            endif

         else if (ires(i).lt.0) then

            sum = 0.0e0
            do j=ibeg(i),iend(i)
               sum = sum + q(j)
            end do
            dis = nint(sum)
            if (abs(dis-sum).gt.0.001e0) corr = .true.

            if (corr) then
               dq = (dis - sum)/dble(n)
               sumn = 0.0e0
               do j=ibeg(i),iend(i)
                  q(j) = q(j) + dq
                  sumn = sumn + q(j)
               end do
               write(iun5,'(a,i5)') 
     &            ' Corrected total charge of residue: ',ires(i)
               write(iun5,'(a,f9.3,a,f9.3)') 
     &            '  Sum charges: old ',sum,' new ',sumn
            endif

         endif
      end do

      write(iun5,'(a)') ' '

      return
      end

      subroutine charge(ec,nac,iac,nad,iad,
     &                  coo,iconn,q,forces,potscl,iopt)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxac=3*mxcon)
      parameter (mxad=9*mxcon)
      common /athlp/  iatoms, mxnat
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      dimension vr(3),ded(3)
      dimension coo(3,*),q(*),iconn(mxcon+1,*),
     &          nac(*),nad(*),iac(mxac,*),iad(mxad,*),
     &          forces(3,*),potscl(*),iopt(*)

      econv = 332.05382e0
      v14sc = 1.0e0 / 1.2e0
      ec = 0.0e0

      do i=1,iatoms

         do j=1,iatoms
            potscl(j) = 1.0e0
         end do

         do j=1,iconn(1,i)
            jj = iconn(1+j,i)
            if (jj.gt.0) then
               potscl(jj) = 0.0e0
            endif
         end do

         do j=1,nac(i)
            potscl(iac(j,i)) = 0.0e0
         end do

         do j=1,nad(i)
            potscl(iad(j,i)) = v14sc
         end do

         do k=i+1,iatoms

            if (potscl(k).ne.0.0e0.and.
     &          (iopt(i).eq.1.or.iopt(k).eq.1) ) then

               do j=1,3
                  vr(j) = coo(j,i) - coo(j,k)
               end do

               if (box) call reddis(vr)

               r2 = vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3)

               if (r2.le.cutof2) then
                  r = sqrt(r2)

                  e = econv * q(i) * q(k) * potscl(k) / r

                  de = -e / r
                  de = de / r

                  do j=1,3
                     ded(j) = de * vr(j)
                  end do

                  ec = ec + e

                  do j=1,3
                     forces(j,i) = forces(j,i) + ded(j)
                     forces(j,k) = forces(j,k) - ded(j)
                  end do

               end if
            end if
         end do
      end do

      return
      end

      subroutine clmond(car,pot,coo,qat)
      implicit real (a-h,o-z), integer( i-n)
      common /athlp/ iatoms, mxnat
      dimension car(3),coo(3,*),qat(*)

      pot = 0.0e0

      do i=1,iatoms
         r2 = dist2(car,coo(1,i))
         r = sqrt(r2)
         pot = pot + qat(i)/r
      end do

      return
      end 

      subroutine qscal(q)
      implicit real (a-h,o-z), integer( i-n)
      common /passq/ scalq
      common /athlp/ iatoms, mxnat
      dimension q(*)

      print*,'Old charges:\\'
      do i=1,iatoms
         print*,i,' ',q(i)
      end do
      print*,'\\'

      print*,'New charges:\\'
      do i=1,iatoms
         q(i) = q(i)*scalq
         print*,i,' ',q(i)
      end do
      print*,'\\'

      return
      end 

