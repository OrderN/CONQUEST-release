      subroutine vdwlj(ev,nac,iac,nad,iad,coo,iconn,ityp,forces,potscl,
     &                 iopt)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxac=3*mxcon)
      parameter (mxad=9*mxcon)
      parameter (mxamb=1590)
      integer ambcls
      common /typpar/ ambchg(mxamb),ambcls(mxamb)
      parameter (mxcls=50)
      common /clspar/ ambvdwr(mxcls),ambvdwe(mxcls),mapagf(mxcls)
      parameter (mxgff=72)
      integer amb2gf
      common /gftyp/  gfvdw(2,mxgff),amb2gf(mxamb)
      logical usecut,usesw
      common /limits/ cutoff, cutof2,cuton, cuton2,conof3,usecut,usesw
      common /athlp/  iatoms, mxnat
      logical box,cell,fast
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      integer*2 ityp
      dimension ded(3),vr(3)
      dimension coo(3,*),ityp(*),iconn(mxcon+1,*),
     &          forces(3,*),potscl(*),iopt(*),
     &          nac(*),nad(*),iac(mxac,*),iad(mxad,*)

      ev = 0.0e0
      v14sc = 0.5e0

      do i=1,iatoms

         i1 = int(ityp(i))
         if (i1.gt.0) then
            il = ambcls(i1)
            vdwr1 = ambvdwr(il)
            vdwe1 = ambvdwe(il)
         elseif (i1.le.0) then
            i1 = iabs(i1)
            vdwr1 = gfvdw(1,i1)
            vdwe1 = gfvdw(2,i1)
         endif

         do j=i+1,iatoms
            potscl(j) = 1.0e0
         end do

         do j=1,iconn(1,i)
            jj = iconn(j+1,i)
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

         if (vdwe1.ne.0.0e0) then

           do k=i+1,iatoms

            if (potscl(k).ne.0.0e0.and.
     &          (iopt(i).eq.1.or.iopt(k).eq.1) ) then

               i2 = int(ityp(k))
               if (i2.gt.0) then
                  kl = ambcls(i2)
                  vdwr2 = ambvdwr(kl)
                  vdwe2 = ambvdwe(kl)
               elseif (i2.le.0) then
                  i2 = iabs(i2)
                  vdwr2 = gfvdw(1,i2)
                  vdwe2 = gfvdw(2,i2)
               endif

               if (vdwe2.ne.0.0e0) then

                  do j=1,3
                     vr(j) = coo(j,i) - coo(j,k)
                  end do
 
                  if (box) call reddis(vr)

                  rv2 = vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3)

                  if (rv2.le.cutof2) then

c          [ (Rmin)**12       (Rmin)**6 ]
c e = eps  [ (----)     - 2.0 (----)    ]
c          [ ( r  )           ( r  )    ]


c alternatively we could precalculate these vdwr(il,kl)
c                                           vdwe(il,kl)
                     rsum = vdwr1 + vdwr2
                     epsm = sqrt(vdwe1 * vdwe2)

                     epsm = epsm * potscl(k)
                     rv   = sqrt(rv2)
		     rs2  = rsum*rsum
		     rs3  = rs2*rsum
                     p6   = (rs3*rs3) / (rv2*rv2*rv2)
                     p12  = p6 * p6

                     e    = epsm * (p12 - 2.0e0*p6)
                     de   = epsm * (p12 - p6) * (-12.0e0/rv)

                     de   = de / rv

                     do j=1,3
                        ded(j) = de * vr(j)
                     end do

                     ev   = ev + e
   
                     do j=1,3
                        forces(j,i) = forces(j,i) + ded(j)
                        forces(j,k) = forces(j,k) - ded(j)
                     end do

                  end if

               end if
            end if
           end do
         end if
      end do

      return
      end

