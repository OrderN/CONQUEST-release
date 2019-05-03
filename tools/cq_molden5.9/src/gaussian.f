      subroutine gaussian(x,y,z,psi,norbs,ipreca,ic,jc)
      implicit double precision (a-h,o-z)
      parameter (numpre=500)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common/prefa/ prexpx(numpre,numcex),prexpy(numpre,numcex),
     &              prexpz(numcex)
      common /congau/ d3,d5,d7,g1,g2,g3,g4,g5,g6,r1,r2,r3,r4,
     &                r3ov2,z1,z2,z3,igcon
      dimension psi(*)

      call gaucon

      do i=1,norbs
        psi(i) = 0.0d0
      end do

      ipnt = 0
      ico = 1

      do 420 ishell=1,nshell

          xk = x - gx(ishell)
          yk = y - gy(ishell)
          zk = z - gz(ishell)
          rr = xk*xk + yk*yk + zk*zk

          do 471 igauss=1,shelln(ishell)

          ipnth = ipnt + igauss

          if (ipreca.eq.1) then
             ex = prexpx(ic,ipnth) *
     &            prexpy(jc,ipnth) *
     &            prexpz(ipnth)

          else
             ex = dexp(-exx(shella(ishell)+igauss-1)*rr)
          endif

          itel = 0
          i = shellt(ishell)
          if (i.eq.4) goto 10
          if (i.eq.3) goto 9
          if ((i.eq.2).and.(shellc(ishell).eq.2)) goto 8
          if ((i.eq.1).and.(shellc(ishell).eq.1)) goto 6

c************* S-FUNCTIONS ******************************

          cons = ex * c1(shella(ishell)+igauss-1)
          psi(ico+itel) = psi(ico+itel) + cons
          itel = itel + 1
          if (i.eq.0) goto 471

C************* P-FUNCTIONS ******************************

   6      conp = ex * c2(shella(ishell)+igauss-1)
          psi(ico+itel) = psi(ico+itel) + xk*conp
          itel = itel + 1
          psi(ico+itel) = psi(ico+itel) + yk*conp
          itel = itel + 1
          psi(ico+itel) = psi(ico+itel) + zk*conp
          itel = itel + 1
          if (i.eq.1) goto 471

C****** D-FUNCTIONS **********************************  

    8     cond = ex * c3(shladf(ishell)+igauss-1)
          if (ido5d.eq.1) then
c D 0
             psi(ico+itel) = psi(ico+itel) +
     &          (-0.5d0*(xk*xk+yk*yk)+zk*zk)*cond
             itel = itel + 1
c D+1
             psi(ico+itel) = psi(ico+itel) +xk*zk*cond*d3
             itel = itel + 1
c D-1
             psi(ico+itel) = psi(ico+itel) +yk*zk*cond*d3
             itel = itel + 1
c D+2
             psi(ico+itel) = psi(ico+itel) +
     &          r3ov2*(xk*xk-yk*yk)*cond
             itel = itel + 1
c D-2
             psi(ico+itel) = psi(ico+itel) +xk*yk*cond*d3
             itel = itel + 1
          else
             psi(ico+itel) = psi(ico+itel) +xk*xk*cond
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) +yk*yk*cond
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) +zk*zk*cond
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) +xk*yk*cond*d3
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) +xk*zk*cond*d3
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) +yk*zk*cond*d3
             itel = itel + 1
          endif
          if(i.eq.2) goto 471

C******* F-FUNCTIONS *********************************    

    9     conf = ex * c4(shladf(ishell)+igauss-1)

          if (ido7f.eq.1) then

             xxx = xk*xk*xk
             yyy = yk*yk*yk
             zzz = zk*zk*zk
             xyy = xk*yk*yk*d5
             xxy = xk*xk*yk*d5
             xxz = xk*xk*zk*d5
             xzz = xk*zk*zk*d5
             yzz = yk*zk*zk*d5
             yyz = yk*yk*zk*d5
             xyz = xk*yk*zk*d5*d3
c F 0
             psi(ico+itel) = psi(ico+itel) + 
     &         (zzz-r2*(xxz+yyz))*conf
             itel = itel + 1
c F+1
             psi(ico+itel) = psi(ico+itel) + 
     &         (z1*xzz-xxx-z2*xyy)*r4*conf
             itel = itel + 1
c F-1
             psi(ico+itel) = psi(ico+itel) + 
     &         (z1*yzz-yyy-z2*xxy)*r4*conf
             itel = itel + 1
c F+2
             psi(ico+itel) = psi(ico+itel) + 
     &         (xxz-yyz)*r3*conf
             itel = itel + 1
c F-2
             psi(ico+itel) = psi(ico+itel) + 
     &         xyz*conf
             itel = itel + 1
c F+3
             psi(ico+itel) = psi(ico+itel) + 
     &         (xxx-z3*xyy)*r1*conf
             itel = itel + 1
c F-3
             psi(ico+itel) = psi(ico+itel) + 
     &         (z3*xxy-yyy)*r1*conf
             itel = itel + 1

          else

             psi(ico+itel) = psi(ico+itel) + xk*xk*xk*conf
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + yk*yk*yk*conf
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + zk*zk*zk*conf
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + xk*yk*yk*conf*d5
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + xk*xk*yk*conf*d5
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + xk*xk*zk*conf*d5
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + xk*zk*zk*conf*d5
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + yk*zk*zk*conf*d5
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + yk*yk*zk*conf*d5
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + xk*yk*zk*conf*d3*d5
             itel = itel + 1

          endif

          goto 471

C******* G-FUNCTIONS *********************************    

   10     cong = ex * c5(shladf(ishell)+igauss-1)

          if (ido9g.eq.1) then

             xxxx = xk*xk*xk*xk
             yyyy = yk*yk*yk*yk
             zzzz = zk*zk*zk*zk
             xxxy = xk*xk*xk*yk
             xxxz = xk*xk*xk*zk
             yyyx = yk*yk*yk*xk
             yyyz = yk*yk*yk*zk
             zzzx = zk*zk*zk*xk
             zzzy = zk*zk*zk*yk
             xxyy = xk*xk*yk*yk
             xxzz = xk*xk*zk*zk
             yyzz = yk*yk*zk*zk
             xxyz = xk*xk*yk*zk
             yyxz = yk*yk*xk*zk
             zzxy = zk*zk*xk*yk
c G 0
             psi(ico+itel) = psi(ico+itel) + 
     &         (0.375d0*xxxx + 0.375d0*yyyy + zzzz
     &          - 3.0d0*xxzz - 3.0d0*yyzz  + 0.75d0*xxyy)*cong
             itel = itel + 1
c G+1
             psi(ico+itel) = psi(ico+itel) + 
     &         g1*(4.0d0*zzzx - 3.0d0*yyxz - 3.0d0*xxxz)*cong
             itel = itel + 1
c G-1
             psi(ico+itel) = psi(ico+itel) + 
     &         g1*(4.0d0*zzzy - 3.0d0*xxyz - 3.0d0*yyyz)*cong
             itel = itel + 1
c G+2
             psi(ico+itel) = psi(ico+itel) + 
     &         g2*(6.0d0*xxzz - 6.0d0*yyzz - xxxx + yyyy)*cong
             itel = itel + 1
c G-2
             psi(ico+itel) = psi(ico+itel) + 
     &         g3*(6.0d0*zzxy - xxxy - yyyx)*cong
             itel = itel + 1
c G+3
             psi(ico+itel) = psi(ico+itel) + 
     &         g4*(xxxz - 3.0d0*yyxz)*cong
             itel = itel + 1
c G-3
             psi(ico+itel) = psi(ico+itel) + 
     &         g4*(3.0d0*xxyz - yyyz)*cong
             itel = itel + 1
c G+4
             psi(ico+itel) = psi(ico+itel) + 
     &         g5*(xxxx - 6.0d0*xxyy + yyyy)*cong
             itel = itel + 1
c G-4
             psi(ico+itel) = psi(ico+itel) + 
     &         g6*(xxxy - yyyx)*cong
             itel = itel + 1

          else

c xxxx
             psi(ico+itel) = psi(ico+itel) + xk*xk*xk*xk*cong
             itel = itel + 1
c yyyy
             psi(ico+itel) = psi(ico+itel) + yk*yk*yk*yk*cong
             itel = itel + 1
c zzzz
             psi(ico+itel) = psi(ico+itel) + zk*zk*zk*zk*cong
             itel = itel + 1
c xxxy
             psi(ico+itel) = psi(ico+itel) + xk*xk*xk*yk*cong*d7
             itel = itel + 1
c xxxz
             psi(ico+itel) = psi(ico+itel) + xk*xk*xk*zk*cong*d7
             itel = itel + 1
c yyyx
             psi(ico+itel) = psi(ico+itel) + yk*yk*yk*xk*cong*d7
             itel = itel + 1
c yyyz
             psi(ico+itel) = psi(ico+itel) + yk*yk*yk*zk*cong*d7
             itel = itel + 1
c zzzx
             psi(ico+itel) = psi(ico+itel) + zk*zk*zk*xk*cong*d7
             itel = itel + 1
c zzzy
             psi(ico+itel) = psi(ico+itel) + zk*zk*zk*yk*cong*d7
             itel = itel + 1
c xxyy
             psi(ico+itel) = psi(ico+itel) + xk*xk*yk*yk*cong*d5*d7/d3
             itel = itel + 1
c xxzz
             psi(ico+itel) = psi(ico+itel) + xk*xk*zk*zk*cong*d5*d7/d3
             itel = itel + 1
c yyzz
             psi(ico+itel) = psi(ico+itel) + yk*yk*zk*zk*cong*d5*d7/d3
             itel = itel + 1
c xxyz
             psi(ico+itel) = psi(ico+itel) + xk*xk*yk*zk*cong*d5*d7
             itel = itel + 1
c yyxz
             psi(ico+itel) = psi(ico+itel) + yk*yk*xk*zk*cong*d5*d7
             itel = itel + 1
c zzxy
             psi(ico+itel) = psi(ico+itel) + zk*zk*xk*yk*cong*d5*d7
             itel = itel + 1

          endif

          goto 471

471       continue

          ipnt = ipnt + shelln(ishell)
          ico = ico + itel

420   continue

      if (ico-1.ne.norbs) 
     & call inferr('GAUSSIAN: number of orbitals incorrect',1)

      return
      end

      subroutine gaucon
      implicit double precision (a-h,o-z)
      common /congau/ d3,d5,d7,g1,g2,g3,g4,g5,g6,r1,r2,r3,r4,
     &                r3ov2,z1,z2,z3,igcon

      if (igcon.eq.1) return

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

      igcon = 1

      return
      end

      subroutine gcres
      implicit double precision (a-h,o-z)
      common /congau/ d3,d5,d7,g1,g2,g3,g4,g5,g6,r1,r2,r3,r4,
     &                r3ov2,z1,z2,z3,igcon

      igcon = 0

      return
      end

      subroutine denhes(x,y,z,psi,grd,hess,norbs,dolap)
      implicit double precision (a-h,o-z)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      logical dolap
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /gauhlp/ d6to5(6,5),f10to7(10,7),g15to9(15,9),
     &                plmn(35),lmn(3,35)
      dimension psi(*),grd(3,*),hess(6,*),d(3),d2(3),
     &          tpsi(10),tgrd(3,10),thess(6,10)


      call denini

      dsmall = 1.0d-10
      do i=1,norbs
        psi(i) = 0.0d0
      end do

      if (dolap) then
         do i=1,norbs
           do j=1,3
              grd(j,i) = 0.0d0
              hess(j,i) = 0.0d0
              hess(j+3,i) = 0.0d0
           end do
         end do
      endif

      ico = 1

      do 420 ishell=1,nshell

          d(1) = x - gx(ishell)
          d(2) = y - gy(ishell)
          d(3) = z - gz(ishell)

          do j=1,3
             if (dabs(d(j)).lt.dsmall) d(j) = dsmall
          end do

          r2 = 0.0d0

          do j=1,3
             d2(j) = d(j)*d(j)
             r2 = r2 + d2(j)
          end do

          i = shellt(ishell)

          do 471 igauss=1,shelln(ishell)

          expon = exx(shella(ishell)+igauss-1)
          ex = dexp(-expon*r2)
          ex2 = 2.0d0*expon
          itel = 0

          if (i.eq.4) goto 10
          if (i.eq.3) goto 9
          if ((i.eq.2).and.(shellc(ishell).eq.2)) goto 8
          if ((i.eq.1).and.(shellc(ishell).eq.1)) goto 6

c************* S-FUNCTIONS ******************************

          cons = ex * c1(shella(ishell)+igauss-1)
          indh = ico + itel
          call calcl(d,d2,psi(indh),grd(1,indh),hess(1,indh),
     &               cons,ex2,1,1,dolap)
          itel = itel + 1
          if (i.eq.0) goto 471

C************* P-FUNCTIONS ******************************

   6      conp = ex * c2(shella(ishell)+igauss-1)
          indh = ico + itel
          call calcl(d,d2,psi(indh),grd(1,indh),hess(1,indh),
     &               conp,ex2,2,4,dolap)
          itel = itel + 3
          if (i.eq.1) goto 471

C****** D-FUNCTIONS **********************************  

    8     cond = ex * c3(shladf(ishell)+igauss-1)

          if (ido5d.eq.1) then

             call cleart(tpsi,tgrd,thess,dolap)
             call calcl(d,d2,tpsi,tgrd,thess,cond,ex2,5,10,dolap)

c D 0, D+1, D-1, D+2, D-2

             do j=1,5
                ipos = ico + itel
                do k=1,6
                   psi(ipos) = psi(ipos) + d6to5(k,j)*tpsi(k)
                   if (dolap) then
                      do l=1,3
                         grd(l,ipos) = grd(l,ipos) + 
     &                                 d6to5(k,j)*tgrd(l,k)
                      end do
                      do l=1,6
                         hess(l,ipos) = hess(l,ipos) + 
     &                                  d6to5(k,j)*thess(l,k)
                      end do
                   endif
                end do
                itel = itel + 1
             end do

          else

c xx,yy,zz,xy,xz,yz

             ipos = ico + itel
             call calcl(d,d2,psi(ipos),grd(1,ipos),hess(1,ipos),
     &                  cond,ex2,5,10,dolap)
             itel = itel + 6

          endif

          if (i.eq.2) goto 471

C******* F-FUNCTIONS *********************************    

    9     conf = ex * c4(shladf(ishell)+igauss-1)

          if (ido7f.eq.1) then

             call cleart(tpsi,tgrd,thess,dolap)
             call calcl(d,d2,tpsi,tgrd,thess,conf,ex2,11,20,dolap)

c F 0, F+1, F-1, F+2, F-2, F+3, F-3

             do j=1,7
                ipos = ico + itel
                do k=1,10
                   psi(ipos) = psi(ipos) + f10to7(k,j)*tpsi(k)
                   if (dolap) then
                      do l=1,3
                         grd(l,ipos) = grd(l,ipos) + 
     &                                 f10to7(k,j)*tgrd(l,k)
                      end do
                      do l=1,6
                         hess(l,ipos) = hess(l,ipos) + 
     &                                  f10to7(k,j)*thess(l,k)
                      end do
                   endif
                end do
                itel = itel + 1
             end do

          else

c xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz

             ipos = ico + itel
             call calcl(d,d2,psi(ipos),grd(1,ipos),hess(1,ipos),
     &                  conf,ex2,11,20,dolap)
             itel = itel + 10

          endif

          goto 471

C******* G-FUNCTIONS *********************************    

   10     cong = ex * c4(shladf(ishell)+igauss-1)

          if (ido9g.eq.1) then

             call cleart(tpsi,tgrd,thess,dolap)
             call calcl(d,d2,tpsi,tgrd,thess,cong,ex2,21,35,dolap)

c G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4

             do j=1,9
                ipos = ico + itel
                do k=1,15
                   psi(ipos) = psi(ipos) + g15to9(k,j)*tpsi(k)
                   if (dolap) then
                      do l=1,3
                         grd(l,ipos) = grd(l,ipos) + 
     &                                 g15to9(k,j)*tgrd(l,k)
                      end do
                      do l=1,6
                         hess(l,ipos) = hess(l,ipos) + 
     &                                  g15to9(k,j)*thess(l,k)
                      end do
                   endif
                end do
                itel = itel + 1
             end do

          else

c xxxx, yyyy, zzzz, xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz
c xxyz yyxz zzxy

             ipos = ico + itel
             call calcl(d,d2,psi(ipos),grd(1,ipos),hess(1,ipos),
     &                  cong,ex2,21,35,dolap)
          endif
c
471       continue
c
c
          ico = ico + itel
  420     continue
      if (ico-1.ne.norbs) 
     & call inferr('GAUSSIAN: number of orbitals incorrect',1)

      return
      end

      subroutine denini
      implicit double precision (a-h,o-z)
      common /gauhlp/ d6to5(6,5),f10to7(10,7),g15to9(15,9),
     &                plmn(35),lmn(3,35)
      data plmn /35*1.0d0/
      data ((lmn(j,i),j=1,3),i=1,35) /
     & 0,0,0, 1,0,0, 0,1,0, 0,0,1, 2,0,0, 0,2,0, 0,0,2, 1,1,0,
     & 1,0,1, 0,1,1, 3,0,0, 0,3,0, 0,0,3, 1,2,0, 2,1,0, 2,0,1,
     & 1,0,2, 0,1,2, 0,2,1, 1,1,1, 4,0,0, 0,4,0, 0,0,4, 3,1,0,
     & 3,0,1, 1,3,0, 0,3,1, 1,0,3, 0,1,3, 2,2,0, 2,0,2, 0,2,2,
     & 2,1,1, 1,2,1, 1,1,2/

      if (d6to5(3,1).eq.1.0d0) return

      d3 = dsqrt(3.0d0)
      d5 = dsqrt(5.0d0)
      d7 = dsqrt(7.0d0)

      plmn(8) = d3
      plmn(9) = d3
      plmn(10) = d3

      plmn(14) = d5
      plmn(15) = d5
      plmn(16) = d5
      plmn(17) = d5
      plmn(18) = d5
      plmn(19) = d5
      plmn(20) = d5*d3
      plmn(24) = d7
      plmn(25) = d7
      plmn(26) = d7
      plmn(27) = d7
      plmn(28) = d7
      plmn(29) = d7
      plmn(30) = d5*d7/d3
      plmn(31) = d5*d7/d3
      plmn(32) = d5*d7/d3
      plmn(33) = d5*d7
      plmn(34) = d5*d7
      plmn(35) = d5*d7

      z1 = 4.0d0/d5
      z2 = 1.0d0/d5
      z3 = 3.0d0/d5
      r1 = 0.5d0*dsqrt(5.0d0/2.0d0)
      r2 = 1.5d0/d5
      r3ov2 = 0.5d0*d3
      r3 = r3ov2
      r4 = 0.5d0*dsqrt(3.0d0/2.0d0)

      do i=1,5
        do j=1,6
           d6to5(j,i) = 0.0d0
        end do
      end do

      d6to5(1,1) = -0.5d0
      d6to5(2,1) = -0.5d0
      d6to5(3,1) = 1.0d0
      d6to5(5,2) = 1.0d0
      d6to5(6,3) = 1.0d0
      d6to5(1,4) = r3ov2
      d6to5(2,4) = -r3ov2
      d6to5(4,5) = 1.0d0

      do i=1,7
        do j=1,10
           f10to7(j,i) = 0.0d0
        end do
      end do

      f10to7(3,1) = 1.0d0
      f10to7(6,1) = -r2
      f10to7(9,1) = -r2
      f10to7(1,2) = -r4
      f10to7(4,2) = -z2*r4
      f10to7(7,2) = z1*r4
      f10to7(2,3) = -r4
      f10to7(5,3) = -z2*r4
      f10to7(8,3) = z1*r4
      f10to7(6,4) = r3
      f10to7(9,4) = -r3
      f10to7(10,5) = 1.0d0
      f10to7(1,6) = r1
      f10to7(4,6) = -z3*r1
      f10to7(2,7) = -r1
      f10to7(5,7) = z3*r1

      do i=1,9
        do j=1,15
           g15to9(j,i) = 0.0d0
        end do
      end do

c G0, g+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4

      g1 = dsqrt(5.0d0/8.0d0)
      g2 = dsqrt(5.0d0/16.0d0)
      g3 = dsqrt(5.0d0/4.0d0)
      g4 = dsqrt(35.0d0/8.0d0)
      g5 = dsqrt(35.0d0/64.0d0)
      g6 = dsqrt(35.0d0/4.0d0)
c G 0
      g15to9(1,1)  = 0.375d0
      g15to9(2,1)  = 0.375d0
      g15to9(3,1)  = 1.0d0
      g15to9(10,1) = 0.75d0
      g15to9(11,1) = -3.0d0
      g15to9(12,1) = -3.0d0
c G+1
      g15to9(5,2)  = -3.0d0*g1
      g15to9(8,2)  =  4.0d0*g1
      g15to9(14,2) = -3.0d0*g1
c G-1
      g15to9(7,3)  = -3.0d0*g1
      g15to9(9,3)  =  4.0d0*g1
      g15to9(13,3) = -3.0d0*g1
c G+2
      g15to9(1,4)  = -1.0d0*g2
      g15to9(2,4)  =  1.0d0*g2
      g15to9(11,4) =  6.0d0*g2
      g15to9(12,4) = -6.0d0*g2
c G-2
      g15to9(4,5)  = -1.0d0*g3
      g15to9(6,5)  = -1.0d0*g3
      g15to9(15,5) =  6.0d0*g3
c G+3
      g15to9(5,6)  = -3.0d0*g4
      g15to9(14,6) =  1.0d0*g4
c G-3
      g15to9(7,7)  = -1.0d0*g4
      g15to9(13,7) =  3.0d0*g4
c G+4
      g15to9(1,8)  =  1.0d0*g5
      g15to9(2,8)  =  1.0d0*g5
      g15to9(10,8) = -6.0d0*g5
c G-4
      g15to9(4,9)  =  1.0d0*g6
      g15to9(6,9)  = -1.0d0*g6

      return
      end

      subroutine calcl(d,d2,tpsi,tgrd,thess,con,ex2,ibeg,iend,dolap)
      implicit double precision (a-h,o-z)
      logical dolap
      common /gauhlp/ d6to5(6,5),f10to7(10,7),g15to9(15,9),
     &                plmn(35),lmn(3,35)
      dimension tpsi(*),tgrd(3,*),thess(6,*),d(3),d2(3),f(3),f2(3)

      do i=1,iend-ibeg+1

          indh = ibeg+i-1
          gau = 1.0d0
          do j=1,3
             l = lmn(j,indh)
c             if (l.ne.0) gau = gau * d(j)**l
             gau = gau * d(j)**l
          end do
          gau = gau * con*plmn(indh)
          tpsi(i) = tpsi(i) + gau
          
          if (.not.dolap) goto 100

          do j=1,3
              l = lmn(j,indh)
              f(j) = dble(l)/d(j) - ex2*d(j)
              f2(j) = -dble(l)/d2(j) - ex2
              tgrd(j,i) = tgrd(j,i) + f(j)*gau
          end do

          do j=1,6
              jh = j+4
              fact = 1.0d0
              do k=1,3
                 l = lmn(k,jh)
c                 if (l.ne.0) fact = fact * f(k)**l
                 fact = fact * f(k)**l
              end do
              if (j.le.3) fact = fact + f2(j)
              thess(j,i) = thess(j,i) + fact*gau
          end do

100       continue
      end do

      return
      end

      subroutine cleart(tpsi,tgrd,thess,dolap)
      implicit double precision (a-h,o-z)
      logical dolap
      dimension tpsi(*),tgrd(3,*),thess(6,*)

      do i=1,10
         tpsi(i) = 0.0d0
      end do

      if (.not.dolap) return

      do i=1,10

         do j=1,3
            tgrd(j,i) = 0.0d0
         end do

         do j=1,6
            thess(j,i) = 0.0d0
         end do

      end do

      return
      end

      subroutine calhed(psi,grd,hess,den,g,hes,p)
c THIS IS REALLY calhes
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /orbhlp/ mxorb,iuhf,ispd

      dimension p(*),psi(*),grd(3,*),hess(6,*), g(3),hes(6)

      den = 0.0d0

      do k=1,3
         g(k) = 0.0d0
      end do

      do k=1,6
         hes(k) = 0.0d0
      end do

      do i=1,norbs
          do j=1,i
             pij = p((j-1)*mxorb+i)
             if (i.ne.j) pij = pij + pij
             den = den + psi(i)*psi(j)*pij
             do k=1,3
               g(k) = g(k) + (psi(i)*grd(k,j) + psi(j)*grd(k,i))*pij
             end do
             hes(1) = hes(1) + (psi(i)*hess(1,j) + psi(j)*hess(1,i)
     &                       + 2.0d0*grd(1,i)*grd(1,j))*pij
             hes(2) = hes(2) + (psi(i)*hess(2,j) + psi(j)*hess(2,i)
     &                       + 2.0d0*grd(2,i)*grd(2,j))*pij
             hes(3) = hes(3) + (psi(i)*hess(3,j) + psi(j)*hess(3,i)
     &                       + 2.0d0*grd(3,i)*grd(3,j))*pij
             hes(4) = hes(4) + pij*(
     &                grd(1,i)*grd(2,j) + psi(j)*hess(4,i) +
     &                grd(2,i)*grd(1,j) + psi(i)*hess(4,j))
             hes(5) = hes(5) + pij*(
     &                grd(1,i)*grd(3,j) + psi(j)*hess(5,i) +
     &                grd(3,i)*grd(1,j) + psi(i)*hess(5,j))
             hes(6) = hes(6) + pij*(
     &                grd(2,i)*grd(3,j) + psi(j)*hess(6,i) +
     &                grd(3,i)*grd(2,j) + psi(i)*hess(6,j))
          end do
      end do

      return
      end

      subroutine gaudxyz(x,y,z,psi,dxpsi,dypsi,dzpsi,
     &                    norbs,ipreca,ic,jc)
      implicit double precision (a-h,o-z)
      parameter (numpre=500)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common/prefa/ prexpx(numpre,numcex),prexpy(numpre,numcex),
     &              prexpz(numcex)
      common /congau/ d3,d5,d7,g1,g2,g3,g4,g5,g6,r1,r2,r3,r4,
     &                r3ov2,z1,z2,z3,igcon
      dimension psi(*),dxpsi(*),dypsi(*),dzpsi(*)

      call gaucon

      do i=1,norbs
        psi(i) = 0.0d0
        dxpsi(i) = 0.0d0
        dypsi(i) = 0.0d0
        dzpsi(i) = 0.0d0
      end do

      ipnt = 0
      ico = 1

      do 420 ishell=1,nshell

          xk = x - gx(ishell)
          yk = y - gy(ishell)
          zk = z - gz(ishell)
          rr = xk*xk + yk*yk + zk*zk
          xk2 = xk*xk
          yk2 = yk*yk
          zk2 = zk*zk
          xyk = xk*yk
          xzk = xk*zk
          yzk = yk*zk
          xyzk = xk*yk*zk
          xk3 = xk*xk*xk
          yk3 = yk*yk*yk
          zk3 = zk*zk*zk

          do 471 igauss=1,shelln(ishell)

          ipnth = ipnt + igauss

          ex = dexp(-exx(shella(ishell)+igauss-1)*rr)
          exh = exx(shella(ishell)+igauss-1)

          itel = 0
          i = shellt(ishell)
          if (i.eq.4) goto 10
          if (i.eq.3) goto 9
          if ((i.eq.2).and.(shellc(ishell).eq.2)) goto 8
          if ((i.eq.1).and.(shellc(ishell).eq.1)) goto 6

c************* S-FUNCTIONS ******************************

          cons = ex * c1(shella(ishell)+igauss-1)
          psi(ico+itel) = psi(ico+itel) + cons
          dxpsi(ico+itel) = dxpsi(ico+itel) 
     &                      - 2.0d0*xk*exh*cons
          dypsi(ico+itel) = dypsi(ico+itel) 
     &                      - 2.0d0*yk*exh*cons
          dzpsi(ico+itel) = dzpsi(ico+itel) 
     &                      - 2.0d0*zk*exh*cons
          itel = itel + 1
          if (i.eq.0) goto 471

C************* P-FUNCTIONS ******************************

   6      conp = ex * c2(shella(ishell)+igauss-1)
          psi(ico+itel) = psi(ico+itel) + xk*conp
          dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*xk2*exh)*conp
          dypsi(ico+itel) = dypsi(ico+itel) 
     &                             - 2.0d0*xyk*exh*conp
          dzpsi(ico+itel) = dzpsi(ico+itel) 
     &                             - 2.0d0*xzk*exh*conp
          itel = itel + 1

          psi(ico+itel) = psi(ico+itel) + yk*conp
          dxpsi(ico+itel) = dxpsi(ico+itel) 
     &                             - 2.0d0*xyk*exh*conp
          dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*yk2*exh)*conp
          dzpsi(ico+itel) = dzpsi(ico+itel) 
     &                             - 2.0d0*yzk*exh*conp
          itel = itel + 1

          psi(ico+itel) = psi(ico+itel) + zk*conp
          dxpsi(ico+itel) = dxpsi(ico+itel) 
     &                             - 2.0d0*xzk*exh*conp
          dypsi(ico+itel) = dypsi(ico+itel) 
     &                             - 2.0d0*yzk*exh*conp
          dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*zk2*exh)*conp
          itel = itel + 1
          if (i.eq.1) goto 471

C****** D-FUNCTIONS **********************************  

    8     cond = ex * c3(shladf(ishell)+igauss-1)
          if (ido5d.eq.1) then
c D 0
             psi(ico+itel) = psi(ico+itel) +
     &          (-0.5d0*(xk*xk+yk*yk)+zk*zk)*cond
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &          (-1.0d0 + (xk2+yk2-2.0d0*zk2)*exh)*xk*cond
             dypsi(ico+itel) = dypsi(ico+itel) +
     &          (-1.0d0 + (xk2+yk2-2.0d0*zk2)*exh)*yk*cond
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &          ( 2.0d0 + (xk2+yk2-2.0d0*zk2)*exh)*zk*cond
             itel = itel + 1
c D+1
             psi(ico+itel) = psi(ico+itel) +xk*zk*cond*d3
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*xk2*exh)*zk*cond*d3
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                             - 2.0d0*xyzk*exh*cond*d3
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*zk2*exh)*xk*cond*d3
             itel = itel + 1
c D-1
             psi(ico+itel) = psi(ico+itel) +yk*zk*cond*d3
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                             - 2.0d0*xyzk*exh*cond*d3
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*yk2*exh)*zk*cond*d3
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*zk2*exh)*yk*cond*d3
             itel = itel + 1
c D+2
             psi(ico+itel) = psi(ico+itel) +
     &          r3ov2*(xk*xk-yk*yk)*cond
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &          r3ov2*(1.0d0 - (xk2-yk2)*exh)*2.0d0*xk*cond
             dypsi(ico+itel) = dypsi(ico+itel) -
     &          r3ov2*(1.0d0 + (xk2-yk2)*exh)*2.0d0*yk*cond
             dzpsi(ico+itel) = dzpsi(ico+itel) -
     &          r3ov2*(xk2-yk2)*exh*2.0d0*zk*cond
             itel = itel + 1
c D-2
             psi(ico+itel) = psi(ico+itel) + xk*yk*cond*d3
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*xk2*exh)*yk*cond*d3
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*yk2*exh)*xk*cond*d3
             dzpsi(ico+itel) = dzpsi(ico+itel)  
     &                             - 2.0d0*xyzk*exh*cond*d3
             itel = itel + 1
          else
             psi(ico+itel) = psi(ico+itel) + xk2*cond
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (2.0d0 - 2.0d0*xk2*exh)*xk*cond
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                              -2.0d0*xk2*exh*yk*cond
             dzpsi(ico+itel) = dzpsi(ico+itel)  
     &                             - 2.0d0*xk2*exh*zk*cond
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) +yk2*cond
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                             - 2.0d0*yk2*exh*xk*cond
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (2.0d0  -2.0d0*yk2*exh)*yk*cond
             dzpsi(ico+itel) = dzpsi(ico+itel)  
     &                             - 2.0d0*yk2*exh*zk*cond
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) +zk2*cond
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                             - 2.0d0*zk2*exh*xk*cond
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                              -2.0d0*zk2*exh*yk*cond
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (2.0d0 - 2.0d0*zk2*exh)*zk*cond
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) +xyk*cond*d3
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*xk2*exh)*yk*cond*d3
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*yk2*exh)*xk*cond*d3
             dzpsi(ico+itel) = dzpsi(ico+itel) 
     &                             - 2.0d0*xyzk*exh*cond*d3
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) +xzk*cond*d3
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*xk2*exh)*zk*cond*d3
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                             - 2.0d0*xyzk*exh*cond*d3
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*zk2*exh)*xk*cond*d3
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) +yzk*cond*d3
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                             - 2.0d0*xyzk*exh*cond*d3
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*yk2*exh)*zk*cond*d3
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*zk2*exh)*yk*cond*d3
             itel = itel + 1
          endif
          if(i.eq.2) goto 471

C******* F-FUNCTIONS *********************************    

    9     conf = ex * c4(shladf(ishell)+igauss-1)

          if (ido7f.eq.1) then

             xxx = xk*xk*xk
             yyy = yk*yk*yk
             zzz = zk*zk*zk
             xyy = xk*yk*yk*d5
             xxy = xk*xk*yk*d5
             xxz = xk*xk*zk*d5
             xzz = xk*zk*zk*d5
             yzz = yk*zk*zk*d5
             yyz = yk*yk*zk*d5
             xyz = xk*yk*zk*d5*d3
c F 0
             psi(ico+itel) = psi(ico+itel) + 
     &         (zzz-r2*(xxz+yyz))*conf
             itel = itel + 1
c F+1
             psi(ico+itel) = psi(ico+itel) + 
     &         (z1*xzz-xxx-z2*xyy)*r4*conf
             itel = itel + 1
c F-1
             psi(ico+itel) = psi(ico+itel) + 
     &         (z1*yzz-yyy-z2*xxy)*r4*conf
             itel = itel + 1
c F+2
             psi(ico+itel) = psi(ico+itel) + 
     &         (xxz-yyz)*r3*conf
             itel = itel + 1
c F-2
             psi(ico+itel) = psi(ico+itel) + 
     &         xyz*conf
             itel = itel + 1
c F+3
             psi(ico+itel) = psi(ico+itel) + 
     &         (xxx-z3*xyy)*r1*conf
             itel = itel + 1
c F-3
             psi(ico+itel) = psi(ico+itel) + 
     &         (z3*xxy-yyy)*r1*conf
             itel = itel + 1

          else

             psi(ico+itel) = psi(ico+itel) + xk3*conf
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (3.0d0 - 2.0d0*xk2*exh)*xk2*conf
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                              -2.0d0*xk3*exh*yk*conf
             dzpsi(ico+itel) = dzpsi(ico+itel)  
     &                             - 2.0d0*xk3*exh*zk*conf
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + yk3*conf
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                              -2.0d0*yk3*exh*xk*conf
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (3.0d0 - 2.0d0*yk2*exh)*yk2*conf
             dzpsi(ico+itel) = dzpsi(ico+itel)  
     &                             - 2.0d0*yk3*exh*zk*conf
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + zk3*conf
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                             - 2.0d0*zk3*exh*xk*conf
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                             - 2.0d0*zk3*exh*zk*conf
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (3.0d0 - 2.0d0*zk2*exh)*zk2*conf
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + xk*yk2*conf*d5
c this one is hard
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*xk2*exh)*yk2*conf*d5
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 - yk2*exh)*2.0d0*xyk*conf*d5
             dzpsi(ico+itel) = dzpsi(ico+itel)  
     &                             - 2.0d0*zxk*exh*yk2*conf*d5
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + xk2*yk*conf*d5
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 - xk2*exh)*2.0d0*xyk*conf*d5
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*yk2*exh)*2.0d0*xk2*conf*d5
             dzpsi(ico+itel) = dzpsi(ico+itel)  
     &                             - 2.0d0*xyk*exh*xk2*conf*d5
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + xk2*zk*conf*d5
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 - xk2*exh)*2.0d0*xzk*conf*d5
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                             - 2.0d0*yzk*exh*xk2*conf*d5
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*zk2*exh)*xk2*conf*d5
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + xk*zk2*conf*d5
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*xk2*exh)*zk2*conf*d5
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                             - 2.0d0*xyk*exh*zk2*conf*d5
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 - zk2*exh)*2.0d0*xzk*conf*d5
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + yk*zk2*conf*d5
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                              -2.0d0*xyk*exh*zk2*conf*d5
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 - yk2*exh)*2.0d0*zk2*conf*d5
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 - zk2*exh)*2.0d0*yzk*conf*d5
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + yk2*zk*conf*d5
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                              -2.0d0*xzk*exh*yk2*conf*d5
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 - yk2*exh)*2.0d0*yzk*conf*d5
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*zk2*exh)*yk2*conf*d5
             itel = itel + 1
             psi(ico+itel) = psi(ico+itel) + xyzk*conf*d3*d5
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*xk2*exh)*yzk*conf*d3*d5
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*yk2*exh)*xzk*conf*d3*d5
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 - 2.0d0*zk2*exh)*xyk*conf*d3*d5
             itel = itel + 1

          endif

          goto 471

C******* G-FUNCTIONS *********************************    

   10     cong = ex * c5(shladf(ishell)+igauss-1)

          if (ido9g.eq.1) then

             xxxx = xk*xk*xk*xk
             yyyy = yk*yk*yk*yk
             zzzz = zk*zk*zk*zk
             xxxy = xk*xk*xk*yk
             xxxz = xk*xk*xk*zk
             yyyx = yk*yk*yk*xk
             yyyz = yk*yk*yk*zk
             zzzx = zk*zk*zk*xk
             zzzy = zk*zk*zk*yk
             xxyy = xk*xk*yk*yk
             xxzz = xk*xk*zk*zk
             yyzz = yk*yk*zk*zk
             xxyz = xk*xk*yk*zk
             yyxz = yk*yk*xk*zk
             zzxy = zk*zk*xk*yk
c G 0
             psi(ico+itel) = psi(ico+itel) + 
     &         (0.375d0*xxxx + 0.375d0*yyyy + zzzz
     &          - 3.0d0*xxzz - 3.0d0*yyzz  + 0.75d0*xxyy)*cong
             itel = itel + 1
c G+1
             psi(ico+itel) = psi(ico+itel) + 
     &         g1*(4.0d0*zzzx - 3.0d0*yyxz - 3.0d0*xxxz)*cong
             itel = itel + 1
c G-1
             psi(ico+itel) = psi(ico+itel) + 
     &         g1*(4.0d0*zzzy - 3.0d0*xxyz - 3.0d0*yyyz)*cong
             itel = itel + 1
c G+2
             psi(ico+itel) = psi(ico+itel) + 
     &         g2*(6.0d0*xxzz - 6.0d0*yyzz - xxxx + yyyy)*cong
             itel = itel + 1
c G-2
             psi(ico+itel) = psi(ico+itel) + 
     &         g3*(6.0d0*zzxy - xxxy - yyyx)*cong
             itel = itel + 1
c G+3
             psi(ico+itel) = psi(ico+itel) + 
     &         g4*(xxxz - 3.0d0*yyxz)*cong
             itel = itel + 1
c G-3
             psi(ico+itel) = psi(ico+itel) + 
     &         g4*(3.0d0*xxyz - yyyz)*cong
             itel = itel + 1
c G+4
             psi(ico+itel) = psi(ico+itel) + 
     &         g5*(xxxx - 6.0d0*xxyy + yyyy)*cong
             itel = itel + 1
c G-4
             psi(ico+itel) = psi(ico+itel) + 
     &         g6*(xxxy - yyyx)*cong
             itel = itel + 1

          else

c xxxx
             psi(ico+itel) = psi(ico+itel) + xk*xk*xk*xk*cong
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (2.0d0 - xk2*exh)*2.0d0*xk3*cong
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                             - 2.0d0*xyk*exh*xk3*cong
             dzpsi(ico+itel) = dzpsi(ico+itel)  
     &                             - 2.0d0*xzk*exh*xk3*cong
             itel = itel + 1
c yyyy
             psi(ico+itel) = psi(ico+itel) + yk*yk*yk*yk*cong
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                             - 2.0d0*yxk*exh*yk3*cong
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (2.0d0 - yk2*exh)*2.0d0*yk3*cong
             dzpsi(ico+itel) = dzpsi(ico+itel)  
     &                             - 2.0d0*yzk*exh*yk3*cong
             itel = itel + 1
c zzzz
             psi(ico+itel) = psi(ico+itel) + zk*zk*zk*zk*cong
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                             - 2.0d0*xzk*exh*zk3*cong
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                             - 2.0d0*yzk*exh*zk3*cong
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (2.0d0 - zk2*exh)*2.0d0*zk3*cong
             itel = itel + 1
c xxxy
             psi(ico+itel) = psi(ico+itel) + xk*xk*xk*yk*cong*d7
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (3.0d0 -2.0d0*xk2*exh)*xk2*yk*cong*d7
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 -2.0d0*yk2*exh)*xk3*cong*d7
             dzpsi(ico+itel) = dzpsi(ico+itel)  
     &                            - 2.0d0*yzk*exh*xk3*cong*d7
             itel = itel + 1
c xxxz
             psi(ico+itel) = psi(ico+itel) + xk*xk*xk*zk*cong*d7
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (3.0d0 -2.0d0*xk2*exh)*xk2*zk*cong*d7
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                            - 2.0d0*yzk*exh*xk3*cong*d7
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 -2.0d0*zk2*exh)*xk3*cong*d7
             itel = itel + 1
c yyyx
             psi(ico+itel) = psi(ico+itel) + yk*yk*yk*xk*cong*d7
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 -2.0d0*xk2*exh)*yk3*cong*d7
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (3.0d0 -2.0d0*yk2*exh)*yk2*xk*cong*d7
             dzpsi(ico+itel) = dzpsi(ico+itel)  
     &                            - 2.0d0*xzk*exh*yk3*cong*d7
             itel = itel + 1
c yyyz
             psi(ico+itel) = psi(ico+itel) + yk*yk*yk*zk*cong*d7
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                            - 2.0d0*xzk*exh*yk3*cong*d7
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (3.0d0 -2.0d0*yk2*exh)*yk2*zk*cong*d7
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 -2.0d0*zk2*exh)*yk3*cong*d7
             itel = itel + 1
c zzzx
             psi(ico+itel) = psi(ico+itel) + zk*zk*zk*xk*cong*d7
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 -2.0d0*xk2*exh)*zk3*cong*d7
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                            - 2.0d0*xyk*exh*zk3*cong*d7
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (3.0d0 -2.0d0*zk2*exh)*zk2*xk*cong*d7
             itel = itel + 1
c zzzy
             psi(ico+itel) = psi(ico+itel) + zk*zk*zk*yk*cong*d7
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                            - 2.0d0*xyk*exh*zk3*cong*d7
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 -2.0d0*yk2*exh)*zk3*cong*d7
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (3.0d0 -2.0d0*zk2*exh)*zk2*yk*cong*d7
             itel = itel + 1
c xxyy
             psi(ico+itel) = psi(ico+itel) + xk*xk*yk*yk*cong*d5*d7/d3
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 -xk2*exh)*2.0d0*yk2*xk*cong*d5*d7/d3
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 -yk2*exh)*2.0d0*xk2*yk*cong*d5*d7/d3
             dzpsi(ico+itel) = dzpsi(ico+itel)  
     &                            - 2.0d0*zk*exh*xk2*yk2*cong*d5*d7/d3
             itel = itel + 1
c xxzz
             psi(ico+itel) = psi(ico+itel) + xk*xk*zk*zk*cong*d5*d7/d3
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 -xk2*exh)*2.0d0*zk2*xk*cong*d5*d7/d3
             dypsi(ico+itel) = dypsi(ico+itel)  
     &                            - 2.0d0*yk*exh*xk2*zk2*cong*d5*d7/d3
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 -zk2*exh)*2.0d0*xk2*zk*cong*d5*d7/d3
             itel = itel + 1
c yyzz
             psi(ico+itel) = psi(ico+itel) + yk*yk*zk*zk*cong*d5*d7/d3
             dxpsi(ico+itel) = dxpsi(ico+itel)  
     &                            - 2.0d0*xk*exh*yk2*zk2*cong*d5*d7/d3
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 -yk2*exh)*2.0d0*zk2*yk*cong*d5*d7/d3
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 -zk2*exh)*2.0d0*yk2*zk*cong*d5*d7/d3
             itel = itel + 1
c xxyz
             psi(ico+itel) = psi(ico+itel) + xk*xk*yk*zk*cong*d5*d7
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 -xk2*exh)*2.0d0*xyzk*cong*d5*d7
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 -2.0d0*yk2*exh)*xk2*zk*cong*d5*d7
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 -2.0d0*zk2*exh)*xk2*yk*cong*d5*d7
             itel = itel + 1
c yyxz
             psi(ico+itel) = psi(ico+itel) + yk*yk*xk*zk*cong*d5*d7
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 -2.0d0*xk2*exh)*yk2*zk*cong*d5*d7
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 -yk2*exh)*2.0d0*xyzk*cong*d5*d7
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 -2.0d0*zk2*exh)*yk2*xk*cong*d5*d7
             itel = itel + 1
c zzxy
             psi(ico+itel) = psi(ico+itel) + zk*zk*xk*yk*cong*d5*d7
             dxpsi(ico+itel) = dxpsi(ico+itel) +
     &                      (1.0d0 -2.0d0*xk2*exh)*zk2*yk*cong*d5*d7
             dypsi(ico+itel) = dypsi(ico+itel) +
     &                      (1.0d0 -2.0d0*yk2*exh)*zk2*xk*cong*d5*d7
             dzpsi(ico+itel) = dzpsi(ico+itel) +
     &                      (1.0d0 -zk2*exh)*2.0d0*xyzk*cong*d5*d7
             itel = itel + 1

          endif

          goto 471

471       continue

          ipnt = ipnt + shelln(ishell)
          ico = ico + itel

420   continue

      if (ico-1.ne.norbs) 
     & call inferr('GAUSSIAN: number of orbitals incorrect',1)

      return
      end

      subroutine exxtod(npr,norbs,vectrs,vectrb,focc,focb,
     &                  eiga,eigb)

c Write wfn file form gaussian shell structure
c Only works for cartesian gaussians (not the pure D,F, or G)
c G not implemented (iupp and ilow have, pnorm)

      implicit double precision (a-h,o-z)
      parameter (numpre=500)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      parameter (numatm=2000)
      parameter (mxel=100)

      
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp

      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp

      common /congau/ d3,d5,d7,g1,g2,g3,g4,g5,g6,r1,r2,r3,r4,
     &                r3ov2,z1,z2,z3,igcon
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /orbhlp/ mxorb,iuhf,ispd
      common /moldat/ natoms, nummo, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)
      character*2 elemnt,swel
      common /elem/  elemnt(mxel)
      logical opfil
      real eiga,eigb
      integer genaos

      dimension ilow(5,5),iupp(5)
      dimension itdat(35),pnorm(35)
      dimension vectrs(*),vectrb(*),focc(*),focb(*),eiga(*),eigb(*)
      dimension cmop(npr,norbs), cmopb(npr,norbs),
     &          itype(npr),icent(npr),expon(npr)
      dimension swel(mxel)

      data itdat/1,2,3,4,5,6,7,8,9,10,11,12,13,17,14,15,18,19,16,20,
     &           21,22,23,24,25,26,27,28,29,30,31,32,33,34,35/
      data iupp/1,4,10,20,35/
      data ilow/5*1,1,2,5,11,21,1,1,5,11,21,1,1,5,11,21,
     &          1,1,5,11,21/
      data pnorm/35*1.0d0/
 
      if (ido5d.eq.1.or.ido7f.eq.1.or.ido9g.eq.1) then
         print*,"Write of wfn not supported for use with pure D,F,G"
         return
      endif

      idum = genaos(.false.)

      do i=1,mxel
         swel(i) = '  '
         if (elemnt(i)(1:1).eq.' ') then
            swel(i)(1:1) = elemnt(i)(2:2)
         else
            swel(i) = elemnt(i)
         endif
      end do

      pnorm(8)  = d3
      pnorm(9)  = d3
      pnorm(10) = d3
      pnorm(14) = d5
      pnorm(15) = d5
      pnorm(16) = d5
      pnorm(17) = d5
      pnorm(18) = d5
      pnorm(19) = d5
      pnorm(20) = sqrt(15.0d0)
      pnorm(24) = d7
      pnorm(25) = d7
      pnorm(26) = d7
      pnorm(27) = d7
      pnorm(28) = d7
      pnorm(29) = d7
      pnorm(30) = d5*d7/d3
      pnorm(31) = d5*d7/d3
      pnorm(32) = d5*d7/d3
      pnorm(33) = d5*d7
      pnorm(34) = d5*d7
      pnorm(35) = d5*d7

      ipr = 0
      jpr = 0

      do is=1,nshell
        istart = ilow(shellt(is)+1,shellc(is)+1)
        iend = iupp(shellt(is)+1)

        do i=istart,iend

          mu = aos(is) - 1 + i

          do j=1,shelln(is)

            ipr = ipr + 1
            itype(ipr) = itdat(i)
            icent(ipr) = jan(is)
            expon(ipr) = exx(shella(is) + j-1)

            if (i.eq.1) then
              fact = pnorm(i) * c1(shella(is)+j-1)
            elseif (i.le.4) then
              fact = pnorm(i) * c2(shella(is)+j-1)
            elseif (i.le.10) then
              fact = pnorm(i) * c3(shladf(is)+j-1)
            elseif (i.le.20) then
              fact = pnorm(i) * c4(shladf(is)+j-1)
            elseif (i.le.35) then
              fact = pnorm(i) * c5(shladf(is)+j-1)
            else
              print*,'Illegal I in exxtop.'
            endIf

            do mo=1,norbs
              cmop(ipr,mo) = fact*vectrs((mo-1)*mxorb+mu)

              if (iuhf.eq.1) then
                 cmopb(ipr,mo) = fact*vectrb((mo-1)*mxorb+mu)
              endif
            end do

          end do
        end do
      end do
      
      escf = 0.0d0
      vt = 1.0d0

      iunit = 56
      if (.not.opfil(56,'molden.wfn',10,1,0,0)) return

      nmo = 0
      do i=1,norbs
         if (abs(focc(i)).gt.1.d-9) nmo = nmo  + 1
         if (abs(focb(i)).gt.1.d-9) nmo = nmo  + 1
      end do

      write(iunit,'(a)') 'Molden generated WFN file'

      write(iunit,
     & "('GAUSSIAN',i15,' MOL ORBITALS',i7,' PRIMITIVES',i9,' NUCLEI')") 
     &           nmo, npr, natoms

      write(iunit,
     &  "(2x,a2,i4,4x,'(CENTRE',I3,')',1x,3f12.8,'  CHARGE =',f5.1)") 
     &  (swel(nat(i)),i,i,(xyz(j,i),j=1,3),dble(nat(i)),i=1,natoms)

      write(iunit,"('CENTRE ASSIGNMENTS  ',20i3)") icent
      write(iunit,"('TYPE ASSIGNMENTS    ',20i3)") itype
      write(iunit,"('EXPONENTS ',5d14.7)") expon

      do i=1,norbs
        if (abs(focc(i)).gt.1.d-9) then
         write(iunit,"('MO',i5,'     MO 0.0        OCC NO = ',f12.7,
     &                '  ORB. ENERGY =', f12.7)") 
     &               i, focc(i), eiga(i)
         write(iunit,'(5d16.8)') (cmop(ipr,i),ipr=1,npr)
        endif
      end do

      if (iuhf.eq.1) then
        do i=1,norbs
         if (abs(focb(i)).gt.1.d-9) then
         write(iunit,"('MO',i5,'     MO 0.0        OCC NO = ',f12.7,
     &                '  ORB. ENERGY =', f12.7)") 
     &               (i+norbs), focb(i), eigb(i)
         write(iunit,'(5d16.8)') (cmopb(ipr,i),ipr=1,npr)
         endif
        end do
      endif

      write(iunit,"('END DATA',/,
     & ' THE  HF ENERGY = ',1f19.12,' THE VIRIAL(-V/T)= ',1f12.8)")
     &       escf,vt

      return
      end

      subroutine calnpr(npr)
      implicit double precision (a-h,o-z)
      parameter (numpre=500)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp

      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp

      npr = 0
      do is = 1,nshell
        isc = shellc(is)
        ist = shellt(is)
        if (ist.eq.1.and.isc.ne.1) then
           itmp = 4
        else
           itmp = ((ist+2)*(ist+1))/2
        endif
        npr = npr + itmp*shelln(is)
      end do

      return
      end

      subroutine wrtwfn
      parameter (numatm=2000)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)

      call calnpr(npr)
      call exxtop(npr,norbs)

      return
      end

