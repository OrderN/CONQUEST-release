      subroutine prprsc(rs,sc,iconn,ityp)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      parameter (mxamb=1590)
      parameter (mxgff=72)
      integer*2 ityp
      common /athlp/  iatoms, mxnat
      common /ambgff/ iambat(mxamb),igffat(mxgff)
      dimension iconn(mxcon+1,*),ityp(*)
      dimension rs(*),sc(*)

      do i=1,iatoms
            itypi = ityp(i)
            if (itypi.lt.0) then
               ia = igffat(itypi)
            else
               ia = iambat(itypi)
            endif
            n = iconn(1,i)
            if (ia.eq.1) then
               rs(i) = 1.25e0
               if (n.ne.0) then
                  k = iconn(2,i)
                  if (ia.eq.7)  rs(i) = 1.15e0
                  if (ia.eq.8)  rs(i) = 1.05e0
               endif
            else if (ia.eq.3) then
               rs(i) = 1.432e0
            else if (ia.eq.6) then
               rs(i) = 1.90e0
               if (n.ne.3) rs(i) = 1.875e0
               if (n.ne.2) rs(i) = 1.825e0
            else if (ia.eq.7) then
               rs(i) = 1.7063e0
               if (n.eq.4)  rs(i) = 1.625e0
               if (n.eq.1)  rs(i) = 1.60e0
            else if (ia.eq.8) then
               rs(i) = 1.535e0
               if (n.eq.1) rs(i) = 1.48e0
            else if (ia.eq.9) then
               rs(i) = 1.47e0
            else if (ia.eq.10) then
               rs(i) = 1.39e0
            else if (ia.eq.11) then
               rs(i) = 1.992e0
            else if (ia.eq.12) then
               rs(i) = 1.70e0
            else if (ia.eq.14) then
               rs(i) = 1.80e0
            else if (ia.eq.15) then
               rs(i) = 1.87e0
            else if (ia.eq.16) then
               rs(i) = 1.775e0
            else if (ia.eq.17) then
               rs(i) = 1.735e0
            else if (ia.eq.18) then
               rs(i) = 1.70e0
            else if (ia.eq.19) then
               rs(i) = 2.123e0
            else if (ia.eq.20) then
               rs(i) = 1.817e0
            else if (ia.eq.35) then
               rs(i) = 1.90e0
            else if (ia.eq.36) then
               rs(i) = 1.812e0
            else if (ia.eq.37) then
               rs(i) = 2.26e0
            else if (ia.eq.53) then
               rs(i) = 2.10e0
            else if (ia.eq.54) then
               rs(i) = 1.967e0
            else if (ia.eq.55) then
               rs(i) = 2.507e0
            else if (ia.eq.56) then
               rs(i) = 2.188e0
            else
               rs(i) = 2.0e0
            end if
      end do

      do i=1,iatoms

          sc(i) = 0.80e0
          itypi = ityp(i)

          if (itypi.lt.0) then
             ia = igffat(itypi)
          else
             ia = iambat(itypi)
          endif

          if (ia.eq.1)   sc(i) = 0.85e0
          if (ia.eq.6)   sc(i) = 0.72e0
          if (ia.eq.7)   sc(i) = 0.79e0
          if (ia.eq.8)   sc(i) = 0.85e0
          if (ia.eq.9)   sc(i) = 0.88e0
          if (ia.eq.15)  sc(i) = 0.86e0
          if (ia.eq.16)  sc(i) = 0.96e0
          if (ia.eq.26)  sc(i) = 0.88e0
      end do

      return
      end

      subroutine born(coo,rs,brad,sc,iconn,ityp)
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      common /athlp/  iatoms, mxnat
      integer*2 ityp
      dimension coo(3,*),iconn(mxcon+1,*),ityp(*)
      dimension rs(*),brad(*),sc(*),vr(3)

      call prprsc(rs,sc,iconn,ityp)

      do i=1,iatoms
         rs(i) = rs(i) - 0.09e0
      end do

c calculate born radii, HTC scheme

      do i=1,iatoms

            ri = rs(i)
            sum = 1.0e0 / ri

            do k=1,iatoms

               if (i.ne.k) then

                  do j=1,3
                     vr(j) = coo(j,k) - coo(j,i)
                  end do

                  r2 = vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3)

                  r = sqrt(r2)
                  rk = rs(k)
                  sk = rk*sc(k)
                  sk2 = sk*sk

                  if (ri.lt.r+sk) then
                     dl = 1.0e0 / dble(max(ri,r-sk))
                     du = 1.0e0 / (r+sk)
                     dl2 = dl * dl
                     du2 = du * du
                     term = dl - du + 0.25e0*r*(du2-dl2)
     &                         + (0.5e0/r)*log(du/dl)
     &                         + (0.25e0*sk2/r)*(dl2-du2)
                     sum = sum - 0.5e0*term
                  end if

               end if

            end do

            brad(i) = 1.0e0 / sum
            brad(i) = dble(max(ri,brad(i)))

      end do

      return
      end

      subroutine esolv(coo,forces,drb,rs,brad,sc,q,
     &                 iconn,ityp)
c calculate born solvation
      implicit real (a-h,o-z), integer (i-n)
      parameter (mxcon=10)
      common /athlp/  iatoms, mxnat
      integer*2 ityp
      dimension coo(3,*),iconn(mxcon+1,*),ityp(*)
      dimension forces(3,*),q(*),drb(*)
      dimension rs(*),brad(*),sc(*)

      es = 0.0e0
      do i = 1, iatoms
         drb(i) = 0.0e0
      end do

c calculate rs, brad and sc

      call born(coo,rs,brad,sc,iconn,ityp)

      term = 16.0e0*datan(1.d0)

      do i = 1, iatoms

c do we have to use + 0.09 here, doffset ??

          ri = rs(i) + 0.09e0
          rb = brad(i)
          if (rb .ne. 0.0e0) then
             e = 0.0054e0 * term * (ri+1.4e0)**2 * (ri/rb)**6
             es = es + e
             drb(i) = drb(i) - 6.0e0*e/rb
          end if
      end do

      call egbsa(coo,forces,drb,q,brad)

      call bself(coo,forces,drb,rs,brad,sc)

      return
      end

      subroutine egbsa(coo,forces,drb,q,brad)
      implicit real (a-h,o-z), integer (i-n)
      common /athlp/  iatoms, mxnat
      dimension coo(3,*),forces(3,*),q(*),drb(*),vr(3),ded(3)
      dimension brad(*)

      dwater = 78.3e0
      f = -332.05382e0 * (1.0e0 - 1.0e0/dwater)

      do i=1,iatoms

         fi = f * q(i)
         rbi = brad(i)

         do k=i,iatoms

               do j=1,3
                  vr(j) = coo(j,i) - coo(j,k)
               end do
               r2 = vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3)

               r = sqrt(r2)
               rbk = brad(k)
               fik = fi * q(k)
               rb2 = rbi * rbk
               ex = exp(-0.25e0*r2/rb2)
               fgb2 = r2 + rb2*ex
               fgb = sqrt(fgb2)
               e = fik / fgb
               de = -e * (r-0.25e0*r*ex) / fgb2
               derb = -e * ex*(0.5e0+0.125e0*r2/rb2) / fgb2

c should this happen ??

               if (i.eq.k) then

                  e = 0.5e0 * e
                  es = es + e
                  drbi = derb * rbk
                  drb(i) = drb(i) + drbi

               else

                  es = es + e
                  de = de / r

                  do j=1,3
                     ded(j) = de * vr(j)
                     forces(j,i) = forces(j,i) + ded(j)
                     forces(j,k) = forces(j,k) - ded(j)
                  end do

                  drbi = derb * rbk
                  drbk = derb * rbi
                  drb(i) = drb(i) + drbi
                  drb(k) = drb(k) + drbk

               end if

         end do
      end do

      return
      end

      subroutine bself(coo,forces,drb,rs,brad,sc)
      implicit real (a-h,o-z), integer (i-n)
      common /athlp/  iatoms, mxnat
      dimension coo(3,*),forces(3,*),drb(*),ded(3),vr(3)
      dimension rs(*),brad(*),sc(*)

      do i =1,iatoms

         ri = rs(i)
         rb2 = brad(i) * brad(i)

         do k =1,iatoms

               if (k.ne.i) then

                  do j=1,3
                     vr(j) = coo(j,k) - coo(j,i)
                  end do

                  rk = rs(k)
                  sk = rk * sc(k)
                  sk2 = sk * sk
                  r2 = vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3)
                  r = sqrt(r2)

                  if (ri.lt.r+sk) then

                     dl = 1.0e0 / dble(max(ri,r-sk))
                     du = 1.0e0 / (r+sk)
                     dl2 = dl * dl
                     du2 = du * du
                     dl3 = dl * dl2
                     du3 = du * du2

                     fac = 1.0e0
                     if (ri.ge.r-sk) fac = 0.0e0

                     t1 = 0.5e0*dl2 + 0.25e0*sk2*dl3/r
     &                       - 0.25e0*(dl/r+dl3*r)
                     t2 = -0.5e0*du2 - 0.25e0*sk2*du3/r
     &                       + 0.25e0*(du/r+du3*r)
                     t3 = 0.125e0*(1.0e0+sk2/r2)*(dl2-du2)
     &                       + 0.25e0*log(du/dl)/r2
                     de = drb(i) * rb2 * (fac*t1+t2+t3) / r

                     do j=1,3
                        ded(j) = de * vr(j)
                     end do

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
