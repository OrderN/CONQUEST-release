      subroutine grdcad(dens,npts1,npts2,iprnt,space,
     &                  p,paa,psi,grd,hess)

c THIS IS REALLY grdcal

c     routine calculates the z-values of the plot plane
c     output array dens and step
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      common /elmcom/ vdwr(mxel),vrad(mxel),icol(mxel)
      logical calcu,valenc,bonds,ovrlap,atomic,doori,
     &  ostep,fine,oscal,dolap,purden,onden
      integer space
      parameter (numatm=2000)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /plspec/ step,stepf,scale,scalf,cut,pmin,pmax,one,
     &                oscal,ostep,fine
      common /eulx/   ca,cb,sa,sb,cc,sc
      common /coord / xyz(3,numatm)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /orbhlp/ mxorb,iuhf,ispd

      common /option/ ipsi,idebug,valenc,bonds,ovrlap,atomic,doori,
     &                dolap
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      dimension dens(*),p(*),paa(*),psi(*),grd(3,*),hess(6,*)
      dimension hes(6),g(3),car(3)

      toang = 0.52917706d0

      purden = onden(idum)

      call curs(1)
      call parrat

      if (iprnt.eq.1) call inferr('Calculation of Grid',0)
      cutoff = 1.0d-7
c      if (.false.) call cmpgau(paa,psi,grd,hess,ipsi)

      radius1 = 0.5d0 * r(1)
      radius2 = 0.5d0 * r(2)
      rf1     = r(1) / dble(npts1-1)
      rf2     = r(2) / dble(npts2-1)
      pmax    = -100000.d0
      pmin    =  100000.d0
      ij      = 0

      do ic=1,npts1

          x2 = -radius1+dble(ic-1)*rf1

          do jc=1,npts2

              ij =  ij + 1
              y2 = -radius2+dble(jc-1)*rf2
              x  = v1(1)*x2 + v2(1)*y2 + px
              y  = v1(2)*x2 + v2(2)*y2 + py
              z  = v1(3)*x2 + v2(3)*y2 + pz

              if (molpot.or.elpot.or.chpot) then

c  Molpot = calculate electrostatic potential using Multipoles 
c  Elpot  = calculate ab initio electrostatic potential 

                 calcu=.true.

                 do l=1,natoms
                     xk = x - xyz(1,l)
                     yk = y - xyz(2,l)
                     zk = z - xyz(3,l)
                     r2 = xk*xk + yk*yk + zk*zk
                     r2 = r2 + 1.d-10
                     r1 = dsqrt(r2)

                     if (space.eq.1) then
                        if (r1.lt.0.1d0) calcu = .false.
                     else
                        if (r1.lt.vdwr(nat(l))/toang) calcu = .false.
                     endif

                 end do

                 if (calcu) then
                   if (chpot) then
                      car(1) = x
                      car(2) = y
                      car(3) = z
                      call clmon(car,dens(ij))
                   else
                      if (molpot) call calc(x,y,z,dens(ij))
                      if (elpot) call espot(x,y,z,dens(ij),0)
                   endif
                 else
                   dens(ij) = 0.0d0
                 endif

              else

                 if (iftyp.eq.1) then
                    call slater(x,y,z,psi)
                 elseif (iftyp.eq.5.and.isgau.eq.0) then
                    call adffun(x,y,z,psi)
                 else
                    if (ipsi.eq.0.and.dolap) then
                       call denhes(x,y,z,psi,grd,hess,norbs,dolap)
                    else
                       call gaussian(x,y,z,psi,norbs,1,ic,jc)
                    endif
                 endif

                 sum = 0.0d0
                 if (ipsi.eq.0)  then
                     if (dolap) then
                        call calhes(psi,grd,hess,den,g,hes)
                        sum = hes(1) + hes(2) + hes(3)
                        if (dabs(sum).lt.cutoff) sum = 0.0d0
                        if (dabs(sum).gt.20.0d0) then
                            if (sum.gt.0.0d0) then
                               sum = 20.0d0
                            else
                               sum = -20.0d0
                            endif
                        endif

                     else

                        if (purden) then
                           call denfst(sum,psi)
                        else
                           do i=1,norbs
                               sum = sum - psi(i)*psi(i)
     &                                *p((i-1)*mxorb+i)*0.5d0
                               do j=1,i
                                   sum = sum + psi(i)*psi(j)
     &                                *p((j-1)*mxorb+i)
                               end do
                           end do
                           sum = sum + sum
                        endif

                     endif

                 else
              

                     do i=1,norbs
                         sum = sum + paa(i)*psi(i)
                     end do

                 endif
                 dens(ij) = sum * one
              endif

              if (dens(ij).gt.pmax) pmax = dens(ij)
              if (dens(ij).lt.pmin) pmin = dens(ij)

           end do
      end do

c
c------------- end grid points loop -----------------------
c

      call curs(0)
      if (iprnt.eq.1) then
      write(iun3,'(//''***** GRID POINT DENSITY/INTENSITY *****''//)')
      write(iun3,'('' MAXIMUM DENSITY/INTENSITY = '',f13.5)')pmax
      write(iun3,'('' MINIMUM DENSITY/INTENSITY = '',f13.5)')pmin
      write(iun3,'(/''ON A TOTAL OF '',i3,''*'',i3,'' POINTS'')')npts1,
     & npts2
      endif
      pmax = max(pmax,-pmin)
      if (iprnt.eq.1) then
         if (pmax.lt.1.d-5) write(iun3,2222)
         if (pmax.gt.1.d 5) write(iun3,2223)
      endif
      if (pmax.lt.1.d-5) pmax = 0.0d0

      call parstp

2222  format(' DENSITY OF PLOT IS VERY LOW-',
     &       ' SUGGEST YOU CHECK YOUR DATA')
2223  format(' DENSITY OF PLOT IS VERY HIGH-',
     &       ' SUGGEST YOU CHECK YOUR DATA')

      return
      end

      subroutine cmpgau(paa,psi,grd,hess,ipsi)
c     test: compare with the gaussian program
c# hf/6-31g** 6d gfinput iop(6/7=3) cube=orbitals \
c  iop(6/49=-1) iop(6/33=2)
c  oh = 0.942774, hoh = 106.053
cgridfile
c  -80    0.000000    0.000000    0.000000
c   -1    1.000000    0.000000    0.000000
c    3    0.000000    1.000000    0.000000
c    3    0.000000    0.000000    1.000000
c1	(orbital number)
c or
c# hf/6-31g** 6d gfinput iop(6/7=3) cube=density iop(6/49=-1) iop(6/23=2)
c iop(6/24=2)
c The last iop is to stop gaussian to disregard core orbitals
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      dimension xyz(3,9),psi(*),grd(3,*),hess(6,*),
     &          hes(6),g(3),paa(*)
      data ((xyz(j,i),j=1,3),i=1,9) /
     & 0.000000, 0.000000, 0.000000,
     & 0.000000, 0.000000, 1.000000,
     & 0.000000, 0.000000, 2.000000,
     & 0.000000, 1.000000, 0.000000,
     & 0.000000, 1.000000, 1.000000,
     & 0.000000, 1.000000, 2.000000,
     & 0.000000, 2.000000, 0.000000,
     & 0.000000, 2.000000, 1.000000,
     & 0.000000, 2.000000, 2.000000/


      print*,'density and its gradient'
      print*,' '
      do i=1,9
        call denhes(xyz(1,i),xyz(2,i),xyz(3,i),
     &              psi,grd,hess,norbs,.true.)
c         call gaussian(xyz(1,i),xyz(2,i),xyz(3,i),psi,norbs)
         if (ipsi.eq.0)  then
             call calhes(psi,grd,hess,den,g,hes)
             write(*,'(i3,1(1pe13.5))') i,den
             write(*,'(i3,3(1pe13.5))') i,(g(l),l=1,3)
             write(*,'(i3,6(1pe13.5))') i,(hes(l),l=1,6)
         else
            sum = 0.0d0
            do j=1,norbs
               sum = sum + paa(j)*psi(j)
            end do
            print*,i,sum
         endif
      end do

      return
      end

      subroutine grdcd(npts1,npts2,iprnt,space,dens)
      implicit double precision (a-h,o-z)
      integer space

c THIS IS REALLY grdcl

      dimension dens(*)

      call preczo
      call precal(npts1,npts2)
      call grdcal(dens,npts1,npts2,iprnt,space)

      return
      end

      subroutine resedd(iedlog)
      implicit double precision (a-h,o-z)

c THIS IS REALLY resedl

      common /grdhlp/ mx3d,mx3d2
      dimension iedlog(*)

      do i=1,mx3d2
         iedlog(i) = 1
      end do

      return
      end

      subroutine precal(npts1,npts2)
      implicit double precision (a-h,o-z)
      parameter (numpre=500)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common/prefa/ prexpx(numpre,numcex),prexpy(numpre,numcex),
     &              prexpz(numcex)
      dimension d1(3),v1t(3),v2t(3)


      radius1 = 0.5d0 * r(1)
      radius2 = 0.5d0 * r(2)
      rf1     = radius1 / dble(npts1-1) * 2.d0
      rf2     = radius2 / dble(npts2-1) * 2.d0

      vl = vlen(v1)
      do i=1,3
         v1t(i) = v1(i) / vl
      end do

      vl = vlen(v2)
      do i=1,3
         v2t(i) = v2(i) / vl
      end do

      ipnt = 0
      do ishell=1,nshell

          d1(1) = px - gx(ishell)
          d1(2) = py - gy(ishell)
          d1(3) = pz - gz(ishell)

          call timpsc(d1,v1t,dx)
          call timpsc(d1,v2t,dy)

          do igauss=1,shelln(ishell)

             iex = shella(ishell) + igauss - 1
             ipnth = ipnt + igauss
             ex  = exx(iex)

             do i=1,npts1

                x2 = -radius1 + dble(i-1)*rf1
   
                xk = x2 + dx

                prexpx(i,ipnth) = dexp(-ex*xk*xk)

             end do

             do j=1,npts2

                y2 = -radius2 + dble(j-1)*rf2
   
                yk = y2 + dy

                prexpy(j,ipnth) = dexp(-ex*yk*yk)

             end do

          end do
          ipnt = ipnt + shelln(ishell)
      end do

      return
      end

      subroutine precz(rz,nptsz,k)
      implicit double precision (a-h,o-z)
      parameter (numpre=500)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common/prefa/ prexpx(numpre,numcex),prexpy(numpre,numcex),
     &              prexpz(numcex)
      dimension d1(3),v3(3)

      radius3 = 0.5d0 * rz
      rf3     = radius3 / dble(nptsz-1) * 2.d0

      v3(1) = cx
      v3(2) = cy
      v3(3) = cz

      vl = vlen(v3)
      do i=1,3
         v3(i) = v3(i) / vl
      end do

      ipnt = 0
      do ishell=1,nshell

          d1(1) = px - gx(ishell)
          d1(2) = py - gy(ishell)
          d1(3) = pz - gz(ishell)

          call timpsc(d1,v3,dz)

          do igauss=1,shelln(ishell)

             iex = shella(ishell) + igauss - 1
             ipnth = ipnt + igauss
             ex  = exx(iex)


             z2 = radius3 - dble(k-1)*rf3
   
             zk = z2 + dz

             prexpz(ipnth) = dexp(-ex*zk*zk)


          end do
          ipnt = ipnt + shelln(ishell)
      end do

      return
      end

      subroutine preczo
      implicit double precision (a-h,o-z)
      parameter (numpre=500)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common/prefa/ prexpx(numpre,numcex),prexpy(numpre,numcex),
     &              prexpz(numcex)
      dimension d1(3),v3(3)

      v3(1) = cx
      v3(2) = cy
      v3(3) = cz

      vl = vlen(v3)
      do i=1,3
         v3(i) = v3(i) / vl
      end do

      ipnt = 0
      do ishell=1,nshell

          d1(1) = px - gx(ishell)
          d1(2) = py - gy(ishell)
          d1(3) = pz - gz(ishell)

          call timpsc(d1,v3,zk)

          do igauss=1,shelln(ishell)

             iex = shella(ishell) + igauss - 1
             ipnth = ipnt + igauss
             ex  = exx(iex)

             prexpz(ipnth) = dexp(-ex*zk*zk)


          end do

          ipnt = ipnt + shelln(ishell)

      end do

      return
      end

      subroutine denfsd(sum,psi,vectrs,vectrb,focc,focb)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /orbhlp/ mxorb,iuhf,ispd
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      dimension psi(*),vectrs(*),vectrb(*),focc(*),focb(*)

      sum = 0.0d0
      do k=1,norbs
         if (focc(k).gt.0.00001d0) then
            summ = 0.0d0
            do i=1,norbs
                summ = summ + 
     &                 vectrs((k-1)*mxorb+i)*psi(i)
            end do
            sum = sum + focc(k)*summ*summ
         endif
      end do

      if (iuhf.eq.1) then
         do k=1,norbs
            if (focb(k).gt.0.0d0) then
               summ = 0.0d0
               do i=1,norbs
                   summ = summ + 
     &                    vectrb((k-1)*mxorb+i)*psi(i)
               end do
               if (ispd.eq.0) then
                   sum = sum + focb(k)*summ*summ
               else
                   sum = sum - focb(k)*summ*summ
               endif
            endif
         end do
      endif

      return
      end
