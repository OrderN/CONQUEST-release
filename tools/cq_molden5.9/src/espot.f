      subroutine espod(x,y,z,epot,idebug,
     &                 p)
c THIS IS REALLY espot
      implicit double precision(a-h,o-z)
      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      parameter (numtmp=4000)
      logical oeerst
      integer shella,shelln,shellt,shellc,shladf,aos
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /bounds/ ilow(5,5),iupp(5)
      common /indxe/  jx(35),jy(35),jz(35),ix(20),iy(20),iz(20)
      common /strt/   oeerst
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /pseudo / ipseud,ivale(numatm)
      common /coord / xyz(3,numatm)
      common /orbhlp/ mxorb,iuhf,ispd
      common /extchg/ exchg(3,3),iexchg(3),nexchg
      dimension tt(4),zz(4),pre(45),ca(35),cb(35),
     &          coeffx(192),coeffy(192),coeffz(192)
      dimension h(numtmp),f(numtmp),h2x(9),h2y(9),h2z(9),
     &          h3x(16),h3y(16),h3z(16)
      dimension p(*)
      data jx/1,2,1,1,3,1,1,2,2,1,4,1,1,2,3,3,2,1,1,2,
     &           5,1,1,4,4,2,1,2,1,3,3,1,3,2,2/
      data jy/1,1,2,1,1,3,1,2,1,2,1,4,1,3,2,1,1,2,3,2,
     &           1,5,1,2,1,4,4,1,2,3,1,3,2,3,2/
      data jz/1,1,1,2,1,1,3,1,2,2,1,1,4,1,1,2,3,3,2,2,
     &           1,1,5,1,2,1,2,4,4,1,3,3,2,2,3/
      data ix/0,4,0,0,8,0,0,4,4,0,12,0,0,4,8,8,4,0,0,4/
      data iy/0,0,4,0,0,8,0,4,0,4,0,12,0,8,4,0,0,4,8,4/
      data iz/0,0,0,4,0,0,8,0,4,4,0,0,12,0,0,4,8,8,4,4/
c iupp 1  4  10 20 35
c      s  3p 6d 10f 15g
c 
c ilow 0    1       2        3      4    <- shellt
c---------------------------------------
c 0    1(s) 1 (sp)  1 (spd)  1      1
c 1    1    2 (p)   5       11     21
c 2    1    1       5 (d)   11     21
c 3    1    1       5       11 (f) 21
c 4    1    1       5       11     21 (g)
c---------------------------------------
c ^
c | shellc
c
      data iupp/1,4,10,20,35/
      data ilow/5*1,1,2,5,11,21,1,1,5,11,21,1,1,5,11,21,
     &          1,1,5,11,21/

      if (oeerst) then
         call epint
         call setcon
         oeerst = .false.
      endif

      epot = 0.0d0
      pi = 4.d0*datan(1.d0)
c
c     loop over shell pairs.
c
      do ishell = 1, nshell
        xa = gx(ishell)
        ya = gy(ishell)
        za = gz(ishell)
        ifrst = shella(ishell)
        ilast = ifrst + shelln(ishell) - 1
        itype = shellt(ishell)
        itype1 = itype + 1
        istart = ilow(itype1,shellc(ishell)+1)
        iend = iupp(itype1)
        iendt = iend

        do jshell = 1,ishell
          xb = gx(jshell)
          yb = gy(jshell)
          zb = gz(jshell)
          jfsrt = shella(jshell)
          jlast = jfsrt + shelln(jshell) - 1
          jtype = shellt(jshell)
          jtype1 = jtype + 1
          jstart = ilow(jtype1,shellc(jshell)+1)
          jend = iupp(jtype1)
          ijtyp = itype1 + jtype1 - 1
          ndim = (iend-istart+1)*(jend-jstart+1)
          iminj = iabs(ishell - jshell)
          n = (itype + jtype) / 2 + 1
          xt = xb - xa
          yt = yb - ya
          zt = zb - za
          rsq = xt*xt + yt*yt + zt*zt
          if (ndim.gt.numtmp) 
     &       print*,'WARNING: exceding array limit in espot'
          do ii=1,ndim
             h(ii) = 0.0d0
          end do
          jendt = jend
c
c         loop over primitive gaussians.
c
          do igauss = ifrst, ilast
            ei = exx(igauss)
            call fcij(itype,ifrst,igauss,shladf(ishell),ca)
            do jgauss = jfsrt,jlast
              iend = iendt
              jend = jendt
              ej = exx(jgauss)
              call fcij(jtype,jfsrt,jgauss,shladf(jshell),cb)
              ep = ei + ej
              exptmp = ej*ei*rsq / ep
              if (exptmp.ge.600.0d0) goto 100
              expp = 2.0d0*dexp(-exptmp)
              tx = (ei*xa + ej*xb) / ep
              ty = (ei*ya + ej*yb) / ep
              tz = (ei*za + ej*zb) / ep
              xta = tx - xa
              yta = ty - ya
              zta = tz - za
              xtb = tx - xb
              ytb = ty - yb
              ztb = tz - zb
              call coeffs(coeffx,xta,xtb,itype1,jtype1)
              call coeffs(coeffy,yta,ytb,itype1,jtype1)
              call coeffs(coeffz,zta,ztb,itype1,jtype1)
              if (ndim.gt.numtmp) 
     &           print*,'WARNING: exceding array limit in espot'
              do ii=1,ndim
                  f(ii) = 0.0d0
              end do
              xct = x - tx
              yct = y - ty
              zct = z - tz
              expont = ep * (xct * xct + yct * yct + zct * zct)
              call rys(n,expont,tt,zz)
              epp = 0.5d0 / ep
              call calct(pre,epp,ijtyp)
              do ii = 1, n
                fac1 = (ep + ep)*tt(ii)
                zf = pi * expp * zz(ii) / ep
                call twocen(h2x,xct,1.d0,pre,fac1,ijtyp)
                call twocen(h2y,yct,1.d0,pre,fac1,ijtyp)
                call twocen(h2z,zct,zf,pre,fac1,ijtyp)
                call thrcen(h3x,h2x,coeffx,itype1,jtype1)
                call thrcen(h3y,h2y,coeffy,itype1,jtype1)
                call thrcen(h3z,h2z,coeffz,itype1,jtype1)

                k = 0
                do i = istart, iend
                  do j = jstart, jend
                    k = k + 1
                    if (k.gt.numtmp) 
     &                print*,'WARNING: exceding array limit in espot'
                    f(k) = f(k)-
     &           (h3x(jx(j)+ix(i))*h3y(jy(j)+iy(i))*h3z(jz(j)+iz(i)))
                  end do
                end do
              end do

              k = 0
              do i = istart, iend
                do j = jstart,jend
                  k = k + 1
                  if (k.gt.numtmp) 
     &              print*,'WARNING: exceding array limit in espot'
                  h(k) = h(k) + f(k)*ca(i)*cb(j)
                end do
              end do
              kend = k

100           continue
              end do
            end do

            call purdf(itype,jtype,istart,jstart,iend,jend,h)
c
c Density matrix contribution
c

            ii = 0
            ist = aos(ishell) - 1
            jst = aos(jshell) - 1
            do i = istart,iend
               do j = jstart,jend
                  ii = ii + 1
                  if (iminj.eq.0) then
                     a1or2 = 1.0d0
                     if (i-j.lt.0) then
                        pt = p((ist+i-1)*mxorb + (jst+j))
                     else
                        pt = p((jst+j-1)*mxorb + (ist+i))
                     endif
                  else
                     a1or2 = 2.0d0
                     pt = p((jst+j-1)*mxorb + (ist+i))
                  endif
                  phelp = pt*a1or2
                  epot = epot + phelp*h(ii)
            if (idebug.eq.1) print*,'v(',ist+i,',',jst+j,')=',h(ii)
               end do
            end do

            jend = jendt
            iend = iendt
        end do
      end do

c sum in nuclear contribution

c
      do ii = 1 ,natoms
          xt = xyz(1,ii) - x
          yt = xyz(2,ii) - y
          zt = xyz(3,ii) - z
          rsq = xt*xt + yt*yt + zt*zt
          if (rsq.ge.1.0d-8) then
             if (ipseud.eq.1) then
                epot = epot + dfloat(ivale(ii))*dsqrt(1.0d0 / rsq)
             else
                epot = epot + dfloat(nat(ii))*dsqrt(1.0d0 / rsq)
             endif
          endif
      end do

      if (nexchg.ne.0) then
          do ii=1,nexchg
             xt = exchg(1,ii) - x
             yt = exchg(2,ii) - y
             zt = exchg(3,ii) - z
             rsq = xt*xt + yt*yt + zt*zt
             if (rsq.ge.1.0d-8) then
                if (iexchg(ii).eq.1) then
                   epot = epot + dsqrt(1.0d0 / rsq)
                else
                   epot = epot - dsqrt(1.0d0 / rsq)
                endif
             endif
          end do
      endif

      return
      end

      subroutine espgrdd(npts1,npts2,npts3,idebug,
     &                   denn,p,xden,yden,zden)
c THIS IS REALLY espgrd
      implicit double precision(a-h,o-z)
      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      parameter (nmcex2=numcex*numcex)
      parameter (numtmp=4000)
      character monstr*100
      character*2 gstr
      logical oeerst
      integer shella,shelln,shellt,shellc,shladf,aos
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /bounds/ ilow(5,5),iupp(5)
      common /indxe/  jx(35),jy(35),jz(35),ix(20),iy(20),iz(20)
      common /strt/   oeerst
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /pseudo / ipseud,ivale(numatm)
      common /coord / xyz(3,numatm)
      common /orbhlp/ mxorb,iuhf,ispd
      common /extchg/ exchg(3,3),iexchg(3),nexchg
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /grdhlp/ mx3d,mx3d2
      dimension tt(4),zz(4),pre(45),ca(35),cb(35),
     &          coeffx(192),coeffy(192),coeffz(192)
      dimension h(numtmp),f(numtmp),h2x(9),h2y(9),h2z(9),
     &          h3x(16),h3y(16),h3z(16),v3(3)
      dimension p(*),denn(*),xden(*),yden(*),zden(*)

      if (oeerst) then
         call epint
         call setcon
         oeerst = .false.
         do i=1,192
            coeffx(i) = 0.d0
         end do
      endif


      call pregrd(npts1,npts2,npts3,xden,yden,zden)

      do kc=1,npts3
          do ic=1,npts1
             do jc=1,npts2
                iadrs = (kc-1)*mx3d2 + ij
                denn(iadrs) = 0.0d0
             end do
          end do
      end do

      pi = 4.d0*datan(1.d0)
c
c     loop over shell pairs.
c
      nn = 0
      iup = 0
      do ishell = 1, nshell

        xa = gx(ishell)
        ya = gy(ishell)
        za = gz(ishell)

        ifrst = shella(ishell)
        ilast = ifrst + shelln(ishell) - 1
        itype = shellt(ishell)
        itype1 = itype + 1
        istart = ilow(itype1,shellc(ishell)+1)
        iend = iupp(itype1)
        iendt = iend

        do jshell = 1,ishell

          xb = gx(jshell)
          yb = gy(jshell)
          zb = gz(jshell)

          jfsrt = shella(jshell)
          jlast = jfsrt + shelln(jshell) - 1
          jtype = shellt(jshell)
          jtype1 = jtype + 1
          jstart = ilow(jtype1,shellc(jshell)+1)
          jend = iupp(jtype1)
          ijtyp = itype1 + jtype1 - 1
          ndim = (iend-istart+1)*(jend-jstart+1)
          iminj = iabs(ishell - jshell)

          n = (itype + jtype) / 2 + 1

          xt = xb - xa
          yt = yb - ya
          zt = zb - za

          rsq = xt*xt + yt*yt + zt*zt

          if (ndim.gt.numtmp) 
     &       print*,'WARNING: exceding array limit in espot'

c deze jendt moet niet in binnenste loop ??
          jendt = jend

          do kc=1,npts3

             monstr = 'Progress Monitor [1-' // gstr(npts3)//'] ' //
     &                 gstr(kc)
             ij = 0

             do ic=1,npts1
   
                do jc=1,npts2


                ij =  ij + 1
                iadrs = (kc-1)*mx3d2 + ij

                x  = xden(iadrs)
                y  = yden(iadrs)
                z  = zden(iadrs)


                do ii=1,ndim
                   h(ii) = 0.0d0
                end do

 
c         loop over primitive gaussians.
 
                do igauss = ifrst, ilast

                  ei = exx(igauss)
                  call fcij(itype,ifrst,igauss,shladf(ishell),ca)

                  do jgauss = jfsrt,jlast

                    iend = iendt
                    jend = jendt

                    ej = exx(jgauss)
                    call fcij(jtype,jfsrt,jgauss,shladf(jshell),cb)

                    ep = ei + ej
                    epr = 1.0d0 / ep

                    exptmp = ej*ei*rsq*epr

                    if (exptmp.ge.600.0d0) goto 100

                    expp = 2.0d0*dexp(-exptmp)
                    tx = (ei*xa + ej*xb) * epr
                    ty = (ei*ya + ej*yb) * epr
                    tz = (ei*za + ej*zb) * epr

                    xta = tx - xa
                    yta = ty - ya
                    zta = tz - za

                    xtb = tx - xb
                    ytb = ty - yb
                    ztb = tz - zb

                    call coeffs(coeffx,xta,xtb,itype1,jtype1)
                    call coeffs(coeffy,yta,ytb,itype1,jtype1)
                    call coeffs(coeffz,zta,ztb,itype1,jtype1)

                    if (ndim.gt.numtmp) 
     &                 print*,'WARNING: exceding array limit in espot'

                    do ii=1,ndim
                        f(ii) = 0.0d0
                    end do

                    xct = x - tx
                    yct = y - ty
                    zct = z - tz

                    expont = ep * (xct * xct + yct * yct + zct * zct)

                    call rys(n,expont,tt,zz)
                    epp = 0.5d0 * epr
                    call calct(pre,epp,ijtyp)
                    ep2 = ep + ep
                    pee = pi*expp*epr

                    do ii = 1, n
                      fac1 = ep2 * tt(ii)
                      zf = pee * zz(ii)
                      call twocen(h2x,xct,1.d0,pre,fac1,ijtyp)
                      call twocen(h2y,yct,1.d0,pre,fac1,ijtyp)
                      call twocen(h2z,zct,zf,pre,fac1,ijtyp)
                      call thrcen(h3x,h2x,coeffx,itype1,jtype1)
                      call thrcen(h3y,h2y,coeffy,itype1,jtype1)
                      call thrcen(h3z,h2z,coeffz,itype1,jtype1)

                      k = 0
                      do i = istart, iend
                        do j = jstart, jend
                          k = k + 1
                          if (k.gt.numtmp) print*,
     &                      'WARNING: exceding array limit in espot'
                          f(k) = f(k) - (h3x(jx(j)+ix(i))*
     &                                   h3y(jy(j)+iy(i))*
     &                                   h3z(jz(j)+iz(i)))

                        end do
                      end do

                    end do

                    k = 0
                    do i = istart, iend
                      do j = jstart,jend
                        k = k + 1

                        if (k.gt.numtmp) print*,
     &                    'WARNING: exceding array limit in espot'

                        h(k) = h(k) + f(k)*ca(i)*cb(j)

                      end do
                    end do


100                 continue

                  end do
                end do

                call purdf(itype,jtype,istart,jstart,iend,jend,h)
c
c Density matrix contribution
c

                ii = 0
                ist = aos(ishell) - 1
                jst = aos(jshell) - 1

                do i = istart,iend
                   do j = jstart,jend
                      ii = ii + 1
                      if (iminj.eq.0) then
                         a1or2 = 1.0d0
                         if (i-j.lt.0) then
                            pt = p((ist+i-1)*mxorb + (jst+j))
                         else
                            pt = p((jst+j-1)*mxorb + (ist+i))
                         endif
                      else
                         a1or2 = 2.0d0
                         pt = p((jst+j-1)*mxorb + (ist+i))
                      endif
                      phelp = pt*a1or2
                      denn(iadrs) =  denn(iadrs) + phelp*h(ii)

                      if (idebug.eq.1) 
     &                   print*,'v(',ist+i,',',jst+j,')=',h(ii)

                   end do
                end do

                jend = jendt
                iend = iendt

c  end do's kc,ic,jc

                end do
              end do
            end do

c  end do's ishell, jshell

        end do
      end do

c sum in nuclear contribution

      do kc=1,npts3
          ij = 0
          do ic=1,npts1
             do jc=1,npts2

                ij =  ij + 1
                iadrs = (kc-1)*mx3d2 + ij

                x  = xden(iadrs)
                y  = yden(iadrs)
                z  = zden(iadrs)

                do ii = 1,natoms

                    xt = xyz(1,ii) - x
                    yt = xyz(2,ii) - y
                    zt = xyz(3,ii) - z

                    rsq = xt*xt + yt*yt + zt*zt

                    if (rsq.ge.1.0d-8) then
                       if (ipseud.eq.1) then
                          denn(iadrs) = denn(iadrs) + 
     &                        dfloat(ivale(ii))*dsqrt(1.0d0 / rsq)
                       else
                          denn(iadrs) = denn(iadrs) + 
     &                        dfloat(nat(ii))*dsqrt(1.0d0 / rsq)
                       endif
                    endif
                end do

                if (nexchg.ne.0) then
                    do ii=1,nexchg
                       xt = exchg(1,ii) - x
                       yt = exchg(2,ii) - y
                       zt = exchg(3,ii) - z
                       rsq = xt*xt + yt*yt + zt*zt
                       if (rsq.ge.1.0d-8) then
                          if (iexchg(ii).eq.1) then
                             denn(iadrs) = denn(iadrs) + 
     &                                       dsqrt(1.0d0 / rsq)
                          else
                             denn(iadrs) = denn(iadrs) - 
     &                                       dsqrt(1.0d0 / rsq)
                          endif
                       endif
                    end do
                endif

             end do
          end do
      end do

      return
      end

      subroutine purdf(itype,jtype,istart,jstart,iend,jend,h)
c
c only iend and jend are really nescessary to return
c
c iupp 1  4  10 20 35
c      s  3p 6d 10f 15g
c 
c ilow 0    1       2        3      4    <- shellt
c---------------------------------------
c 0    1(s) 1 (sp)  1 (spd)  1      1
c 1    1    2 (p)   5       11     21
c 2    1    1       5 (d)   11     21
c 3    1    1       5       11 (f) 21
c 4    1    1       5       11     21 (g)
c---------------------------------------
c ^
c | shellc
c
c question does h come filled in shell form ?, eg d is really spd
c with start at 5 (instead of 1) (we treat d such but not f?)
c
c xxxx,yyyy,zzzz *1
c xxxy ...       *d7
c xxyy ...       *d5*d7/d3
c xxyz ...       *d5*d7
c
c via fcij ca and cb have these, and therfor h() has them 
c
      implicit double precision (a-h,o-z)
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /intcon/ d3,d5,d7,r1,r2,r3,r4,r3ov2,z1,z2,z3,
     &                g1,g2,g3,g4,g5,g6
      logical debug
      dimension h(*),inc(14)
      debug = .false.
c
c     convert  gaussians to pure angular functions.
c     (6D -> 5D, 10F -> 7F, 15G -> 9G)
c

      b1 = 1.0d0/d7
      b2 = d3/(d5*d7)
      b3 = 1.0d0/(d5*d7)

      idim = iend - istart + 1
      jdim = jend - jstart + 1

      if (ido5d.eq.1.and.jtype.eq.2) then
c
c     The row side of the matrix: pure d
c     

         i1 = 5 - jstart + 1

c     when d   shell jstart = 5 , i1 = 1, jend = 10, jdim = 6
c     when spd shell jstart = 1 , i1 = 5, jend = 10, jdim = 10
c     single d shell is stored in H() as:
c
c     d1 d2 d3 d4 d5 d6 NOT s px py pz d1 d2 etc

         do i=1,idim
            dz2 = h(i1+2) - 0.5d0*(h(i1) + h(i1+1))
            dx2y2 = r3ov2*(h(i1) - h(i1+1))
            h(i1  ) = dz2
            h(i1+1) = h(i1+4)
            h(i1+2) = h(i1+5)
            h(i1+4) = h(i1+3)
            h(i1+3) = dx2y2
            i1 = i1 + jdim
         end do

      endif
 
      if (ido7f.eq.1.and.jtype.eq.3) then
c
c     The row side af the matrix: pure f
c     
c     when f  shell jstart = 11 , jend = 20, jdim = 10
c
         i1 = 0
         do i=1,idim
 
            f0  = h(i1+3) - r2*(h(i1+6) + h(i1+9))
            f1p = r4*(z1*h(i1+7) - h(i1+1) - z2*h(i1+4))
            f1m = r4*(z1*h(i1+8) - h(i1+2) - z2*h(i1+5))
            f2p = r3ov2*(h(i1+6) - h(i1+9))
            f2m = h(i1+10)
            f3p = r1*(h(i1+1) - z3*h(i1+4))
            f3m = r1*(z3*h(i1+5) - h(i1+2))
 
            h(i1+1) = f0
            h(i1+2) = f1p
            h(i1+3) = f1m
            h(i1+4) = f2p
            h(i1+5) = f2m
            h(i1+6) = f3p
            h(i1+7) = f3m
 
            i1 = i1 + jdim

         end do
      endif

      if (ido9g.eq.1.and.jtype.eq.4) then
c
c     The row side af the matrix: pure g
c     
c     when g  shell jstart = 21 , jend = 35, jdim = 15
c
         i1 = 0
         do i=1,idim
 
          g0  = 
     &       0.375d0*h(i1+1)  + 0.375d0*  h(i1+2) +           h(i1+3) 
     &     -3.0d0*b2*h(i1+11) - 3.0d0*b2*h(i1+12) + 0.75d0*b2*h(i1+10)

          g1p = 
     &  g1*(4.0d0*b1*h(i1+8)  - 3.0d0*b3*h(i1+14) - 3.0d0*b1* h(i1+5))

          g1m = 
     &  g1*(4.0d0*b1*h(i1+9)  - 3.0d0*b3*h(i1+13) - 3.0d0*b1* h(i1+7))

          g2p = 
     &  g2*(6.0d0*b2*h(i1+11) - 6.0d0*b2*h(i1+12) -           h(i1+1) 
     &    +          h(i1+2))

          g2m = 
     &  g3*(6.0d0*b3*h(i1+15) - b1*      h(i1+4)  - b1*       h(i1+6))

          g3p = 
     &  g4*(   b1*   h(i1+5)  - 3.0d0*b3*h(i1+14))

          g3m = 
     &  g4*(3.0d0*b3*h(i1+13) -  b1*     h(i1+7))
 
          g4p = 
     &  g5*(         h(i1+1)  - 6.0d0*b2*h(i1+10) +           h(i1+2))

          g4m = 
     &  g6*b1*(      h(i1+4)  -          h(i1+6))

            h(i1+1) = g0
            h(i1+2) = g1p
            h(i1+3) = g1m
            h(i1+4) = g2p
            h(i1+5) = g2m
            h(i1+6) = g3p
            h(i1+7) = g3m
            h(i1+8) = g4p
            h(i1+9) = g4m
 
            i1 = i1 + jdim

         end do
      endif
c
c     the row size of the matrix has changed, so get rid of
c     the superflous functions.
c
      if ((ido5d.eq.1.and.jtype.eq.2).or.
     &    (ido7f.eq.1.and.jtype.eq.3).or.
     &    (ido9g.eq.1.and.jtype.eq.4)) then

         if (jtype.eq.2) then
            jendp = 9
         else if (jtype.eq.3) then
c     when f   shell jstart = 11 , jend = 17, jpure = jdim = 7, jend = 17
            jendp = 17
         else if (jtype.eq.4) then
c     when g   shell jstart = 21 , jendp = 29, jpure = jdim = 9, jend = 29
            jendp = 29
         endif
         if (debug) print*,'row.old.new ',itype,jtype,jend,jendp
         jrpure = jendp - jstart + 1
 
         i1 = 0
         i2 = 0
 
         do i=1,idim
            do j=1,jrpure
               h(i2+j) = h(i1+j)
            end do
            i1 = i1 + jdim
            i2 = i2 + jrpure
         end do
         jend = jendp
         jdim = jrpure

      endif
c
c     transformation at row side complete, start of column side
c

      if (ido5d.eq.1.and.itype.eq.2) then
c
c     The column side of the matrix: pure d
c
         i1 = (5-istart)*jdim + 1
         do i=1,5
            inc(i) = i*jdim
         end do
 
         do j=1,jdim
            dz2 = h(i1+inc(2)) - 0.5d0*(h(i1) + h(i1+inc(1)))
            dx2y2 = r3ov2*(h(i1) - h(i1+inc(1)))
            h(i1       ) = dz2
            h(i1+inc(1)) = h(i1+inc(4))
            h(i1+inc(2)) = h(i1+inc(5))
            h(i1+inc(4)) = h(i1+inc(3))
            h(i1+inc(3)) = dx2y2
            i1 = i1 + 1
         end do

      endif
 
      if (ido7f.eq.1.and.itype.eq.3) then
c
c     The column side af the matrix: pure F
c
         i1 = 1
         do i=1,9
            inc(i) = i*jdim
         end do
 
         do j=1,jdim
 
            f0  = h(i1+inc(2)) - r2*(h(i1+inc(5)) + h(i1+inc(8)))
            f1p = r4*(z1*h(i1+inc(6)) - h(i1) - z2*h(i1+inc(3)))
            f1m = r4*(z1*h(i1+inc(7)) - h(i1+inc(1))-z2*h(i1+inc(4)))
            f2p = r3ov2*(h(i1+inc(5)) - h(i1+inc(8)))
            f2m = h(i1+inc(9))
            f3p = r1*(h(i1) - z3*h(i1+inc(3)))
            f3m = r1*(z3*h(i1+inc(4)) - h(i1+inc(1)))
    
            h(i1       ) = f0
            h(i1+inc(1)) = f1p
            h(i1+inc(2)) = f1m
            h(i1+inc(3)) = f2p
            h(i1+inc(4)) = f2m
            h(i1+inc(5)) = f3p
            h(i1+inc(6)) = f3m
    
            i1 = i1 + 1
         end do

      endif
 
      if (ido9g.eq.1.and.itype.eq.4) then
c
c     The column side af the matrix: pure G
c
         i1 = 1
         do i=1,14
            inc(i) = i*jdim
         end do
 
         do j=1,jdim
 
          g0  = 
     &          0.375d0*h(i1)         + 0.375d0*  h(i1+inc(1))  +
     &                  h(i1+inc(2))  - 3.0d0*b2* h(i1+inc(10)) - 
     &         3.0d0*b2*h(i1+inc(11)) + 0.75d0*b2*h(i1+inc(9))

          g1p = 
     &     g1*(4.0d0*b1*h(i1+inc(7))  - 3.0d0*b3* h(i1+inc(13)) - 
     &         3.0d0*b1*h(i1+inc(4)))

          g1m = 
     &     g1*(4.0d0*b1*h(i1+inc(8))  - 3.0d0*b3* h(i1+inc(12)) - 
     &         3.0d0*b1*h(i1+inc(6)))

          g2p = 
     &     g2*(6.0d0*b2*h(i1+inc(10)) - 6.0d0*b2* h(i1+inc(11)) -
     &                  h(i1)         +           h(i1+inc(1)))

          g2m = 
     &     g3*(6.0d0*b3*h(i1+inc(14)) - b1*       h(i1+inc(3))  -
     &               b1*h(i1+inc(5)))

          g3p = 
     &     g4*(b1*      h(i1+inc(4))  - 3.0d0*b3* h(i1+inc(13)))

          g3m = 
     &     g4*(3.0d0*b3*h(i1+inc(12)) - b1*       h(i1+inc(6)))
 
          g4p = 
     &      g5*(        h(i1)         - 6.0d0*b2* h(i1+inc(9))  + 
     &                  h(i1+inc(1)))

          g4m = 
     &      g6*b1*(     h(i1+inc(3))  -           h(i1+inc(5)) )

            h(i1       ) = g0
            h(i1+inc(1)) = g1p
            h(i1+inc(2)) = g1m
            h(i1+inc(3)) = g2p
            h(i1+inc(4)) = g2m
            h(i1+inc(5)) = g3p
            h(i1+inc(6)) = g3m
            h(i1+inc(7)) = g4p
            h(i1+inc(8)) = g4m
 
            i1 = i1 + 1
         end do

      endif
      if ((ido5d.eq.1.and.itype.eq.2) .or.
     &    (ido7f.eq.1.and.itype.eq.3) .or.
     &    (ido9g.eq.1.and.itype.eq.4)) then

         if (debug) print*,'colm.old ',itype,jtype,iend
         if (itype.eq.2) then
            iend = 9
         else if (itype.eq.3) then
            iend = 17
         else if (itype.eq.4) then
            iend = 29
         endif
         if (debug) print*,'colm.new ',itype,jtype,iend

      endif
 
      return
      end

      subroutine setcon
      implicit double precision (a-h,o-z)
      common /intcon/ d3,d5,d7,r1,r2,r3,r4,r3ov2,z1,z2,z3,
     &                g1,g2,g3,g4,g5,g6
      d3 = dsqrt(3.0d0)
      d5 = dsqrt(5.0d0)
      d7 = dsqrt(7.0d0)
      z1 = 4.0d0/d5
      z2 = 1.0d0/d5
      z3 = 3.0d0/d5
      r1 = 0.5d0*dsqrt(5.0d0/2.0d0)
      r2 = 1.5d0/d5
      r3ov2 = 0.5d0*dsqrt(3.0d0)
      r3 = r3ov2
      r4 = 0.5d0*dsqrt(3.0d0/2.0d0)

      g1 = dsqrt(5.0d0/8.0d0)
      g2 = dsqrt(5.0d0/16.0d0)
      g3 = dsqrt(5.0d0/4.0d0)
      g4 = dsqrt(35.0d0/8.0d0)
      g5 = dsqrt(35.0d0/64.0d0)
      g6 = dsqrt(35.0d0/4.0d0)

      return
      end

      subroutine dipold(ntt,nor,dmao,focc,focb,vectrs,vectrb,p)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)

c     dipole integrals 

      logical debug
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /bounds/ ilow(5,5),iupp(5)
      common /orbhlp/ mxorb,iuhf,ispd
      common /coord / xyz(3,numatm)

      integer shella,shelln,shellt,shellc,shladf,aos
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp

      common /indxe/  jx(35),jy(35),jz(35),ix(20),iy(20),iz(20)
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      double precision imprd
      integer genaos

      dimension focc(*),focb(*),vectrs(*),vectrb(*),nor(*)
      dimension ca(35),cb(35),sx(36),sy(36),sz(36),dd(3,100),
     &    ccx(330),ccy(330),ccz(330),p(*),dmao(ntt,3),
     &    d12(100),dipon(3)

      lind(i,j) = ((max0(i,j)*(max0(i,j)-1)/2)+min0(i,j))

c iupp 1  4  10 20 35
c      s  3p 6d 10f 15g
c 
c ilow 0    1       2        3      4    <- shellt
c---------------------------------------
c 0    1(s) 1 (sp)  1 (spd)  1      1
c 1    1    2 (p)   5       11     21
c 2    1    1       5 (d)   11     21
c 3    1    1       5       11 (f) 21
c 4    1    1       5       11     21 (g)
c---------------------------------------
c ^
c | shellc
c

      debug = .false.
      cf = 2.5417463d0
      idum = genaos(.false.)

      pi = 4.d0*datan(1.d0)

      do i=1,3
         dipo(i) = 0.d0
         dipon(i) = 0.d0
      end do

      do i=1,natoms
         do j=1,3
            dipon(j) = dipon(j) + dble(nat(i))*xyz(j,i)
         end do
      end do

      call denini

      do ishell=1,nshell

         xi = gx(ishell)
         yi = gy(ishell)
         zi = gz(ishell)

         ifrst = shella(ishell)
         ilast = ifrst + shelln(ishell) - 1
         itype = shellt(ishell)

         itype1 = itype + 1
         lidim = itype1 + 2

         istart = ilow(itype1,shellc(ishell)+1)
         iend   = iupp(itype1)
         iendt  = iend

         do jshell=1,ishell

            xj = gx(jshell)
            yj = gy(jshell)
            zj = gz(jshell)

            jfrst = shella(jshell)
            jlast = jfrst + shelln(jshell) - 1
            jtype = shellt(jshell)

            jtype1 = jtype+1
            ljdim  = jtype+2

            jstart = ilow(jtype1,shellc(jshell)+1)
            jend   = iupp(jtype1)
            jendt  = jend
            lpdim  = lidim + ljdim

            xt = xi - xj
            yt = yi - yj
            zt = zi - zj

            r2 = xt*xt + yt*yt + zt*zt

            do l=1,100
               do k=1,3
                  dd(k,l) = 0.d0
               end do
            end do
 
            irp = iendt - istart + 1
            jrp = jendt - jstart + 1
            lentqp = irp*jrp

            if (ishell.eq.jshell) lentqp = (irp*(irp+1))/2

            do igauss=ifrst,ilast

               ei = exx(igauss)
               twoei = ei + ei

               call fcij(itype,ifrst,igauss,shladf(ishell),ca)
 
               do jgauss=jfrst,jlast

                  ej = exx(jgauss)
                  twoej = ej + ej

                  call fcij(jtype,jfrst,jgauss,shladf(jshell),cb)
    
                  ep = ei + ej

                  tx = (ei*xi + ej*xj) / ep
                  ty = (ei*yi + ej*yj) / ep
                  tz = (ei*zi + ej*zj) / ep

                  xit = tx - xi
                  yit = ty - yi
                  zit = tz - zi

                  xjt = tx - xj
                  yjt = ty - yj
                  zjt = tz - zj

                  expont = (ej*ei*r2) / ep
                  pexp   = 0.d0

                  if (expont.lt.600.d0) pexp = dexp(-expont)
 
                  call dipint(lidim,ljdim,lpdim,
     &                  xit,xjt,ccx, yit,yjt,ccy, zit,zjt,ccz,
     &                  ei,ej,pi,pexp,sx,sy,sz)
   
                  call dipao((ishell.eq.jshell),istart,iend,jstart,jend,
     &                 lidim,ljdim,ca,cb,
     &                 sx,sy,sz,xi,yi,zi,dd)

               end do
            end do

            do l=1,3
               call purdf(itype,jtype,istart,jstart,iend,jend,dd(l,1))
            end do

            ii = 0
            ist = aos(ishell) - 1
            jst = aos(jshell) - 1

            do i=istart,iend
               jnd = jend
               if (ishell.eq.jshell) jnd = i

               do j=jstart,jnd
                  ii = ii + 1
                  d12(ii) = p((ist+i-1)*mxorb + (jst+j))
                  if (ist+i.ne.jst+j) d12(ii) = 2.0d0*d12(ii)
                  ll = lind(ist+i,jst+j)
                  do l=1,3
                     dmao(ll,l) =  dd(l,ii)
                  end do
               end do

            end do

            do l=1,3
               dipo(l) = dipo(l) - imprd(lentqp,d12,dd,l)
            end do

         end do
      end do
 
      if (debug) then
         write(6,'(a)') 'Dipole moment vector (Debye)'
         write(6,'(a,3f10.4)') 'elec dipole ',(cf*dipo(i),i=1,3)
         write(6,'(a,3f10.4)') 'nuc  dipole ',(cf*dipon(i),i=1,3)

         dipt = 0.d0
         do l=1,3
            dipo(l) = dipo(l) + dipon(l)
            dipt = dipt + dipo(l)*dipo(l)
         end do
         dipt = dsqrt(dipt)
         write(6,'(a,3f10.4)') 'tot  dipole ',(cf*dipo(i),i=1,3)
         write(6,'(a,f10.4)') 'total dipole scalar ',cf*dipt
      endif
   
      call setnor(nocc,norbs,nor,focc)
      call boys(nocc,nor,dmao,vectrs)

      if (iuhf.ne.0) then
          call setnor(nocc,norbs,nor,focb)
          call boys(nocc,nor,dmao,vectrb)
      endif

      return
      end

      subroutine dipint(l1max,l2max,l12mxp,
     &                  ax,bx,ccx, ay,by,ccy, az,bz,ccz,
     &                  as,bs,pi,pexp,sx,sy,sz)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension ccx(l12mxp,l2max,l1max)
      dimension ccy(l12mxp,l2max,l1max)
      dimension ccz(l12mxp,l2max,l1max)
      dimension sx(l2max,l1max),sy(l2max,l1max),sz(l2max,l1max)
      dimension s1(21)
 
      do l1=1,l1max
         do l2=1,l2max
            do l12=1,l12mxp
               ccx(l12,l2,l1) = 0.d0
               ccy(l12,l2,l1) = 0.d0
               ccz(l12,l2,l1) = 0.d0
            end do
         end do
      end do

      ccx(2,1,1) = 1.d0
      ccy(2,1,1) = 1.d0
      ccz(2,1,1) = 1.d0

      do l1=1,l1max
         do l2=1,l2max
            iwmax = l1 + l2
            do iw=2,iwmax

               if (l1.gt.1) then

                   ccx(iw,l2,l1) = 
     &                  ax*ccx(iw,l2,l1-1) + ccx(iw-1,l2,l1-1)
                   ccy(iw,l2,l1) = 
     &                  ay*ccy(iw,l2,l1-1) + ccy(iw-1,l2,l1-1)
                   ccz(iw,l2,l1) = 
     &                  az*ccz(iw,l2,l1-1) + ccz(iw-1,l2,l1-1)

               elseif (l2.gt.1) then

                   ccx(iw,l2,l1) = 
     &                  bx*ccx(iw,l2-1,l1) + ccx(iw-1,l2-1,l1)
                   ccy(iw,l2,l1) =  
     &                  by*ccy(iw,l2-1,l1) + ccy(iw-1,l2-1,l1)
                   ccz(iw,l2,l1) = 
     &                  bz*ccz(iw,l2-1,l1) + ccz(iw-1,l2-1,l1)

               endif

            end do
         end do
      end do

      temp = 1.d0/(2.d0*(as+bs))
      s1(1) = dsqrt(pi/(as+bs))

      lim = l1max + l2max - 1

      if (lim.ge.2) then

         fact = 1.d0
   
         do i=2,lim
            if (mod(i,2).eq.0) then
               s1(i) = 0.d0
            else
               s1(i) = s1(i-2)*temp*fact
               fact   = fact + 2.d0
            endif
         end do

      endif

      do l1=1,l1max
         do l2=1,l2max

            x = 0.d0
            y = 0.d0
            z = 0.d0
            iwmax = l1 + l2 - 1

            do iw=1,iwmax
               x = x + ccx(iw+1,l2,l1)*s1(iw)
               y = y + ccy(iw+1,l2,l1)*s1(iw)
               z = z + ccz(iw+1,l2,l1)*s1(iw)
            end do

            sx(l2,l1) = x
            sy(l2,l1) = y
            sz(l2,l1) = z*pexp

         end do
      end do

      return
      end

      subroutine dipao(ieqj,istart,iend,jstart,jend,
     &                 l1max,l2max,ca,cb,sx,sy,sz,xa,ya,za,dd)
      implicit double precision (a-h,o-z), integer (i-n)
      logical ieqj
      common /indxe/  lx(35),ly(35),lz(35),kx(20),ky(20),kz(20)
      dimension ca(*),cb(*),dd(3,*)
      dimension sx(l2max,l1max),sy(l2max,l1max),sz(l2max,l1max)
      dimension rxyz(3)

      ii = 0

      do i=istart,iend

         ix = lx(i)
         iy = ly(i)
         iz = lz(i)
 
         jnd = jend

         if (ieqj) jnd = i

         do j=jstart,jnd

            jx = lx(j)
            jy = ly(j)
            jz = lz(j)

            ii = ii + 1
            coef = ca(i)*cb(j)

            rxyz(1) = (xa*sx(jx,ix) + sx(jx,ix+1))
     &                *sy(jy,iy)*sz(jz,iz)
            rxyz(2) = (ya*sy(jy,iy) + sy(jy,iy+1))
     &                *sx(jx,ix)*sz(jz,iz)
            rxyz(3) = (za*sz(jz,iz) + sz(jz,iz+1))
     &                *sx(jx,ix)*sy(jy,iy)

            do l=1,3
               dd(l,ii) = dd(l,ii) + coef*rxyz(l)
            end do

         end do
      end do

      return
      end

      double precision function imprd(n,a,b,l)
      implicit double precision (a-h,o-z), integer (i-n)
c     improduct of a and b.
      dimension a(*), b(3,*)
 
      imprd = 0.d0

      if (n.lt.1) return

      do i=1,n
          imprd = imprd + a(i)*b(l,i)
      end do

      return
      end

      subroutine boyd(norb,nor,dmao,cc,nbasis,ntt,cl,
     &                iord,iir,rij,qpix,qpjx)
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numatm=2000)
      logical debug
 
C     Boys Localization.  After QCPE program 354
C
C     NOrb   ... Number of orbitals to localize.
C     CC     ... Input and output orbitals (NBasis,).
C     CL     ... Scratch matrix (NBasis,NOrb).

      common /orbhlp/ mxorb,iuhf,ispd
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      dimension cc(*), cl(nbasis,*), dmao(ntt,3), nor(nbasis,nbasis),
     &          iord(*), iir(*), rij(ntt,3), qpix(*), qpjx(*)

      lind(i,j) = ((max0(i,j)*(max0(i,j)-1)/2)+min0(i,j))
      kind(j,i) = ((i-1)*mxorb + j)

c essential that in above definition i and j have swapped places !!!


      debug = .false.
      cf = 2.5417463d0

      do itry=1,3

          do i=1,norb
              ii = i

              do j=1,i
                  jj = j
                  sumx = 0.0d0
                  sumy = 0.0d0
                  sumz = 0.0d0

                  do k=1,nbasis
                      do l=1,nbasis
                          cpr = cc(kind(k,ii))*cc(kind(l,jj))
                          sumx = sumx + cpr*dmao(lind(k,l),1)
                          sumy = sumy + cpr*dmao(lind(k,l),2)
                          sumz = sumz + cpr*dmao(lind(k,l),3)
                      end do
                  end do

                  rij(lind(i,j),1) = sumx
                  rij(lind(i,j),2) = sumy
                  rij(lind(i,j),3) = sumz

              end do

          end do


      call caltra(norb,ntt,rij,tracia)

      if (debug) write(6,'(a,f15.8)') 
     &          ' Initial TraceA=',tracia*cf*cf

      do i=1,nbasis
          do j=1,nbasis
             cl(i,j) = 0.0d0
          end do
      end do

      do i=1,norb
          cl(i,i) = 1.0d0
      end do

      iter = 0
      shift = datan(1.0d0)

      change = 1.0d0

       do while (iter.lt.150.and.change.gt.1.d-6)

         if (change.ge.1.d-10) then
            change = 0.0d0
   
            iter = iter + 1

            do i=1,norb
                iir(i) = i
            end do

            nnn = norb

            do i=1,norb
                iii = irnumb(nnn)
                iord(i) = iir(iii)
                iir(iii) = iir(nnn)
                nnn = nnn - 1
            end do

         endif

         do iii=1,norb
   
           i  = iord(iii)
           if (nor(i,i).eq.0) then

             li = lind(i,i)
             
             jm = 1
             rm = 0.0d0
             tm = 0.0d0
             sm = 0.0d0
             cm = 1.0d0

             do j=1,norb

               if (i.ne.j.and.nor(i,j).eq.0) then

                 lj = lind(j,j)
                 ij = lind(i,j)
                 t  = 0.0d0
                 tx = 0.0d0
   
                 do kk=1,3
                     t = t + 4.0d0*rij(ij,kk)**2 - rij(li,kk)**2
     &                  - rij(lj,kk)**2
     &                  + 2.0d0*rij(li,kk)*rij(lj,kk)
                     tx = tx + rij(ij,kk)*(rij(lj,kk)-rij(li,kk))
                 end do

                 if (dabs(t).gt.1.d-10.or.dabs(tx).gt.1.d-10) then

                    tx = 4.0d0*tx
                    t = datan2(tx,t)/4.0d0
                    sign = 1.0d0
                    if (t.gt.0.0d0) sign = -1.0d0
                    t = t + sign*shift
                    itim = 0

                    s = dsin(t)
                    c = dcos(t)
                    rin = 0.0d0
                    itim = itim + 1

                    do kk=1,3
                        qpi = c*c*rij(li,kk) + s*s*rij(lj,kk)
     &                        + 2.0d0*c*s*rij(ij,kk)
                        qpj = c*c*rij(lj,kk)+s*s*rij(li,kk)
     &                        - 2.0d0*c*s*rij(ij,kk)
                        rin = rin + qpi*qpi+qpj*qpj - rij(li,kk)**2
     &                            - rij(lj,kk)**2
                    end do

                    ttest = dabs(t) - shift

                    if (dabs(t).le.1.d-8.or.dabs(ttest).lt.1.d-8.or.
     &                  rin.ge.(-1.d-8)) then
                        if (rin.gt.rm) then
                            rm = rin
                            jm = j
                            sm = s
                            cm = c
                            tm = t
                        endif
                    else

                       if (debug) then
                         write(6,'(a,i4,a,i4,a,f15.8,a,f15.8)') 
     &               ' BoyLoc:  No rotation increases integrals for I=',
     &                    i,' J=',j,' Theta=',t,' Change=',rin
                       endif

                    endif

                 endif

               endif

             end do

c end  loop over j

             rin = rm
             j = jm
             s = sm
             c = cm
             t = tm

             ij = lind(i,j)
             lj = lind(j,j)

             if (nor(i,j).eq.0) then

                change = change + t*t

                do kk=1,3
                    qpi = c*c*rij(li,kk) + s*s*rij(lj,kk)
     &                                   + 2.0d0*c*s*rij(ij,kk)
                    qpj = c*c*rij(lj,kk) + s*s*rij(li,kk)
     &                                   - 2.0d0*c*s*rij(ij,kk)
                    qpij = (c*c-s*s)*rij(ij,kk) + 
     &                      c*s*(rij(lj,kk)-rij(li,kk))

                    do k=1,norb
                        if (i.ne.k.and.j.ne.k) then
                           ik = lind(i,k)
                           jk = lind(j,k)
                           qpix(k) = c*rij(ik,kk) + s*rij(jk,kk)
                           qpjx(k) = c*rij(jk,kk) - s*rij(ik,kk)
                        endif
                    end do

                    do k=1,norb
                        if (i.ne.k.and.j.ne.k) then
                           ik = lind(i,k)
                           jk = lind(j,k)
                           rij(ik,kk) = qpix(k)
                           rij(jk,kk) = qpjx(k)
                        endif
                    end do

                    rin = rin + qpi + qpj - rij(li,kk) - rij(lj,kk)
                    rij(li,kk) = qpi
                    rij(lj,kk) = qpj
                    rij(ij,kk) = qpij

                end do

                do k=1,norb
                    c1 =  c*cl(k,i) + s*cl(k,j)
                    c2 = -s*cl(k,i) + c*cl(k,j)
                    cl(k,i) = c1
                    cl(k,j) = c2
                end do

             endif

           endif
         end do

c end loop over iii
 
         change = dsqrt(2.0d0*change/dble(norb*(norb-1)))

   
         if (debug) then

            write(6,'(a,i4,a,f15.6)') 
     &              ' BoyLoc:  Iteration',iter,' Change=',change


         endif

       end do

      end do

      if (iter.gt.150) then
            write(6,'(a)') 
     &        ' Localization failed. after 3 tries of 150 iterations.'
      else
         if (change.le.1.d-6) then
           write(6,'(a,i4,a)')
     &          ' Localization complete after',iter,' iterations.'
         endif
      endif

      do i=1,norb

          do k=1,nbasis
             qpix(k) = 0.0d0
          end do

          do j=1,norb
              jj = j

              do k=1,nbasis
                  qpix(k) = qpix(k) + cc(kind(k,jj))*cl(j,i)
              end do

          end do

          do k=1,nbasis
             cl(k,i) = qpix(k)
          end do

      end do

c     put rotated orbitals back 

      do i=1,norb
          ii = i

          do k=1,nbasis
             cc(kind(k,ii)) = cl(k,i)
          end do
         
      end do

      if (debug) then

         print*,"rotated orbitals"
         print*," "
         call prev(cc,nbasis,nbasis,mxorb)

         call caltra(norb,ntt,rij,tracfa)
         write(6,'(a,f15.8)') ' Final TraceA=',tracfa

      endif

      call caldrv(norb,ntt,rij,dmax,imax,jmax)

      if (debug) then
          write(6,'(a,i3,a,i3,a,f15.8)') 
     &    ' Largest derivative:  DDelta(',imax,',',jmax,
     &    ')=', dmax
      endif

      if (dmax.gt.1.d-8) then
          write(6,'(a,i3,a,i3,a,f15.8)') 
     &    ' Internal error in BoyLoc, DDelta(',imax,',',jmax,
     &    ')=',dmax
      endif

      return
      end

      integer function irnumb(max)
      implicit double precision (a-h,o-z), integer (i-n)

      irnumb = int(random()*dble(max)+1.0d0)

      return
      end

      subroutine caltra(norb,ntt,rij,traca)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension rij(ntt,3)

      lind(i,j) = ((max0(i,j)*(max0(i,j)-1)/2)+min0(i,j))

      traca = 0.0d0

      ii = 0
      do i=1,norb
         ii = ii + i
         if (i.ne.1) then
         traca = traca + rij(ii,1)**2 + rij(ii,2)**2 + rij(ii,3)**2
         endif
      end do

      return
      end

      subroutine caldrv(norb,ntt,rij,dmax,imax,jmax)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension rij(ntt,3)

      lind(i,j) = ((max0(i,j)*(max0(i,j)-1)/2)+min0(i,j))

      dmax = -1.0d0

      do i=2,norb

          li = lind(i,i)

          do j=1,i-1

              lj = lind(j,j)
              ij = lind(i,j)
              dt = 0.0d0

              do kk=1,3
                  dt = dt + rij(ij,kk)*(rij(li,kk)-rij(lj,kk))
              end do

              if (dabs(dt).gt.dmax) then
                  dmax = dabs(dt)
                  imax = i
                  jmax = j
              endif

           end do
      end do

      dmax = 2.0d0*dmax

      return
      end

      double precision function random()
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (mplier=16807,modlus=2147483647,mobymp=127773,
     &           momdmp=2836)
      common  /seed/jseed,ifrst,nextn
c     mseed comes from alloc.c
      integer hvlue, lvlue, testv, nextn, mseed
 
      if (ifrst.eq.0) then
        jseed = mseed()
        nextn = jseed
        ifrst = 1
      endif
 
      hvlue = nextn / mobymp
      lvlue = mod(nextn, mobymp)
      testv = mplier*lvlue - momdmp*hvlue

      if (testv.gt.0) then
        nextn = testv
      else
        nextn = testv + modlus
      endif

      random = dble(nextn)/dble(modlus)
 
      return
      end

      integer function ncoree()
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      dimension ncore(103)
      data ncore /2*0,2*1,6*1,2*5,6*5,2*9,10*9,9,5*14,
     &            2*18,10*18,18,5*23,2*27,14*27,10*34,
     &            34,5*39,2*43,14*43,50/
 
      ncoree = 0

      do i=1,natoms
          ncoree = ncoree + ncore(nat(i))
      end do

      return
      end

      subroutine setnor(nocc,nbasis,nor,focc)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension focc(*),nor(nbasis,nbasis)

      nocc = 0
      do i=1,nbasis
         if (focc(i).ne.0.d0) then 
            nocc = nocc + 1
         endif
      end do

      do i=1,nbasis
         do j=1,nbasis
            nor(i,j) = 1
         end do
      end do

      mcore = ncoree()
      do i=mcore+1,nocc
         do j=mcore+1,nocc
            nor(i,j) = 0
         end do
      end do

      return
      end

      subroutine pregrd(npts1,npts2,npts3,xden,yden,zden)
      implicit double precision(a-h,o-z)
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat
      common /grdhlp/ mx3d,mx3d2
      dimension v3(3),xden(*),yden(*),zden(*)

      v3(1) = cx
      v3(2) = cy
      v3(3) = cz

      call vnrm(v3)

      radius1 = 0.5d0 * r(1)
      radius2 = 0.5d0 * r(2)
      radius3 = 0.5d0 * r(3)
      rf1     = r(1) / dble(npts1-1)
      rf2     = r(2) / dble(npts2-1)
      rf3     = r(3) / dble(npts3-1)

      do kc=1,npts3
          z2 = -radius3+dble(kc-1)*rf3
          do ic=1,npts1
             x2 = -radius1+dble(ic-1)*rf1
             do jc=1,npts2
                y2 = -radius2+dble(jc-1)*rf2

                iadrs = (kc-1)*mx3d2 + ij

                xden(iadrs)  = v1(1)*x2 + v2(1)*y2 + px - v3(1)*z2
                yden(iadrs)  = v1(2)*x2 + v2(2)*y2 + py - v3(2)*z2
                zden(iadrs)  = v1(3)*x2 + v2(3)*y2 + pz - v3(3)*z2

             end do
          end do
      end do

      return
      end

