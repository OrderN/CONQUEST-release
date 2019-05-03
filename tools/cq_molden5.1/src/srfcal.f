      subroutine srfcal(idomap,denn,fmap,coo,lwrit,ityp,icont,ncont)
c
c     routine calculates the z-values of the plot plane
c     output array dens and step
c
      implicit double precision (a-h,o-z)
      parameter (numat1=20000)
      parameter (mxmmul=100)
      common /athlp/ iatoms, mxnat

      common /comsrf/  vo(3), r(3),v1(3),v2(3),v3(3), wo(3),
     &                 sl(3),isl
      common /vropt/ ivtwo,ihand,ivadd
      common /grdhlp/ mx3d,mx3d2
      common /types/ iff
      integer srfmap,srfloc
      integer*2 ityp
      common /hlpsrf/ npts1,npts2,npts3,srfmap,srfloc,ifogl,itsrf
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      logical dozme
      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme
      common /multim/ imulm, nmulm,ihasqm(mxmmul)

      dimension denn(*),fmap(*),coo(3,*),icont(*),lwrit(*),ityp(*)
      dimension rt(3)

      call curs(1)

      npts1t = npts1
      npts2t = npts2
      npts3t = npts3

      if (ipdbon.eq.1) then
         npts1 = 45
         npts2 = 45
         npts3 = 45
      endif

      itsrf = 1
      if (idomap.ne.0) then
         if (ipdbon.eq.1) then
            if (ihasq.eq.0) then
               ifftmp = iff
               do i=1,iatoms
                  lwrit(i) = ityp(i)
               end do
               iff = 3
               call dotyp(0)
               call chkbck(0)
               iff = ifftmp
               do i=1,iatoms
                  ityp(i) = lwrit(i)
               end do
            endif
         else
            if (ihasq.eq.0) call eem(3,1,istat)
         endif
      endif

      call defsrf

      radius1 = 0.5d0 * r(1)
      radius2 = 0.5d0 * r(2)
      radius3 = 0.5d0 * r(3)
      rf1     = radius1 / (npts1-1.d0) * 2.d0
      rf2     = radius2 / (npts2-1.d0) * 2.d0
      rf3     = radius3 / (npts3-1.d0) * 2.d0

      radmx = radius1
      if (radius2.gt.radmx) radmx = radius2
      if (radius2.gt.radmx) radmx = radius3

      ncont = 0
      do i=1,iatoms
         xk = (vo(1) - coo(1,i))
         yk = (vo(2) - coo(2,i))
         zk = (vo(3) - coo(3,i))
         r2 = xk*xk+yk*yk+zk*zk
         if (r2.lt.2*radmx*radmx) then
            if (ncont.le.numat1) then
               ncont = ncont + 1
               icont(ncont) = i
            endif
         endif
      end do

c v1, v2 , v3 zijn genormeerde (norm is 1) richtings vectoren

      do kc=1,npts3
         z2=-radius3+(kc-1)*rf3
         ij     = 0
         do ic=1,npts2
             y2=-radius2+(ic-1)*rf2
             do jc=1,npts1
                 ij =  ij + 1
                 x2 = -radius1+(jc-1)*rf1
                 rt(1) = v1(1)*x2 + v2(1)*y2 + v3(1)*z2 + vo(1)
                 rt(2) = v1(2)*x2 + v2(2)*y2 + v3(2)*z2 + vo(2)
                 rt(3) = v1(3)*x2 + v2(3)*y2 + v3(3)*z2 + vo(3)
                 call srfden(rt(1),rt(2),rt(3),psi)
                 if (psi.gt.rmax) rmax = psi
                 denn((kc-1)*mx3d2 + ij) = psi
                 if (idomap.ne.0) then
                    if (idomap.eq.1) then
                       call espot(rt(1),rt(2),rt(3),f,0)
                    elseif (idomap.eq.2.and.ihasq.ne.0) then
                       call calc(rt(1),rt(2),rt(3),f)
                    elseif (idomap.eq.3) then
                       call clmons(rt,f,srfloc)
                    endif
                    fmap((kc-1)*mx3d2 + ij) = f
                 endif
             end do
         end do
      end do
c
c------------- end grid points loop -----------------------
c
c maak origin oorsprong zodat we makkelijk kunnen roteren en
c dan pas de translatie naar de oorsprong kunnen doen
c
      call setcod(0.05d0)

      mapit = 0
      if (idomap.ne.0) mapit = 1
      if (ihasq.eq.0.and.idomap.eq.3) mapit = 0

      call mapsrf(denn,fmap,mapit)

      npts1 = npts1t
      npts2 = npts2t
      npts3 = npts3t

      return
      end

      subroutine srfrd(denn,fmap,coo)
      implicit double precision (a-h,o-z)
      dimension denn(*),fmap(*),coo(3,*)

c      call setcod(0.05d0)

      mapit = 0
      call mapsrf(denn,fmap,mapit)

      return
      end

      subroutine mapsrf(denn,fmap,mapit)

c     routine calculates the z-values of the plot plane
c     output array dens and step

      implicit double precision (a-h,o-z)
      logical dolap, lapdbl
      common /vropt/ ivtwo,ihand,ivadd
      integer srfmap,srfloc
      common /hlpsrf/ npts1,npts2,npts3,srfmap,srfloc,ifogl,itsrf
      common /psa/ psa, tsa, exs, pol, pol2, epmin, epmax, ipsa,
     &             icpsa, idtpsa
      dimension denn(*),fmap(*)
      dimension origin(3),rt(3)
c------------- end grid points loop -----------------------
c
c maak origin oorsprong zodat we makkelijk kunnen roteren en
c dan pas de translatie naar de oorsprong kunnen doen
c
      call curs(1)

      toang = 0.52917706d0
      toang3 = toang*toang*toang

      valcnt = 0.05d0
      if (ipsa.eq.1) valcnt = 0.01 / toang

      do i=1,3
         origin(i) = 0.0d0
         rt(i) = 1.0d0
      end do

      ipsi = 0
      dolap =.false.
      lapdbl = .false.
      itrans = 0
      iun = 25

      ivtwot = ivtwo
      ivtwo = 4

c      call sribcol(7)

      call mcubes(npts1,npts2,npts3,denn,fmap,mapit,origin,valcnt,
     &                      -valcnt,rt,ipsi,dolap,lapdbl,itrans,iun)

      ivtwo = ivtwot
      call curs(0)

      return
      end

      subroutine srfded(x,y,z,f,ianz,isurf,coo,icont,ncont)
      implicit double precision (a-h,o-z)
      parameter (nsupp=15)
      parameter (mxel=100)
      common /srfcom/ nsupa(mxel),ngel(nsupp),npq(nsupp),nels(2,nsupp),
     &                sexp(2,nsupp),snrm(2,nsupp)
      dimension ianz(*),isurf(*),coo(3,*),icont(*)

      data ngel  /1,1,1,1,1,1,1,1,1,2,2,2,2,2,2/
      data npq   /0,1,1,1,1,1,1,1,1,1,1,1,1,3,4/
      data nels  /1,0,
     &            1,0,
     &            2,0,
     &            3,0,
     &            4,0,
     &            5,0,
     &            6,0,
     &            7,0,
     &            3,0,
     &            2,2,
     &            2,3,
     &            2,4,
     &            2,5,
     &            2,5,
     &            2,5/
      data sexp  /2.51710141d0,0.0d0,
     &            1.32733145d0,0.0d0,
     &            1.89771849d0,0.0d0,
     &            2.84749616d0,0.0d0,
     &            3.37802055d0,0.0d0,
     &            4.26257495d0,0.0d0,
     &            5.10217946d0,0.0d0,
     &            5.38296417d0,0.0d0,
     &            2.72912143d0,0.0d0,
     &            2.48690111d0,3.23138631d0,
     &            3.98498017d0,3.37475693d0,
     &            4.37094906d0,3.79680895d0,
     &            7.15208053d0,3.84805364d0,
     &            7.28371554d0,4.15598308d0,
     &            4.29535666d0,4.09983616d0/
      data snrm  /2.25307905d0,0.0d0,
     &            0.661170083d0,0.0d0,
     &            1.61601061d0,0.0d0,
     &            4.45679442d0,0.0d0,
     &            6.83156889d0,0.0d0,
     &            12.2192572d0,0.0d0,
     &            19.1537608d0,0.0d0,
     &            21.8987062d0,0.0d0,
     &            4.0079462d0,0.0d0,
     &            3.17695472d0,10.5899771d0,
     &            10.325953d0,11.8039962d0,
     &            13.0108029d0,15.8478254d0,
     &            44.5600556d0,16.3879870d0,
     &            241.467384d0,33.4857583d0,
     &            20.3083480d0,27.2240458d0/
      data nsupa
     &  /1,0,2,3,4,5,6,7,8,0,0,0,9,10,11,12,13,17*0,14,17*0,15,47*0/
c H  1
c 
c Li 2
c Be 3
c B  4
c C  5
c N  6
c O  7
c F  8
c Ne?
c Na?
c Mg?
c Al 9
c Si 10
c P  11
c S  12
c CL 13
c ....
c Br 14
c I  15

c map elemnts to supported atoms array 

c loop over atoms

      toang = 0.52917706d0
      cutoff = 6.0d0
      f = 0.0d0

      do l=1,ncont

         i = icont(l)
         isupe = nsupa(ianz(i))
         if (isupe.lt.0.or.isupe.gt.nsupp) isupe = 0
         if (isurf(i).ne.0.and.isupe.ne.0) then

            xk = (x - coo(1,i))*toang
            yk = (y - coo(2,i))*toang
            zk = (z - coo(3,i))*toang

            r2 = xk*xk+yk*yk+zk*zk

c only calculate contribution if distance is less than a cutoff

            if (r2.lt.cutoff) then
               r2 = r2+1.d-10
               r1 = dsqrt(r2)
               r22 = 1.0d0
               if (npq(isupe).eq.1) r22 = r2
               s1 = snrm(1,isupe)
               s2 = snrm(2,isupe)
               g = nels(1,isupe)*s1*s1*r22*fexp(-2.0d0*sexp(1,isupe)*r1)
               if (ngel(isupe).eq.2)
     &           g = g + 
     &             nels(2,isupe)*s2*s2*r22*fexp(-2.0d0*sexp(2,isupe)*r1)
               f = f + g
            endif

         endif

      end do
      
      return
      end

      subroutine defsrd(isurf,ianz,coo)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      common /comsrf/  vo(3), r(3),v1(3),v2(3),v3(3), wo(3),
     &                 sl(3),isl
      dimension isurf(*),coo(3,*),ianz(*)

c vectors v1,v2,v3 need to be normalised 

      do i=1,3
        vo(i) = 0.0d0
        v1(i) = 0.0d0
        v2(i) = 0.0d0
        v3(i) = 0.0d0
      end do

      v1(1) = 1.0d0
      v2(2) = 1.0d0
      v3(3) = 1.0d0

      sl(1) = 1.0d0
      sl(2) = 1.0d0
      sl(3) = 1.0d0

      wo(1) = -0.5d0
      wo(2) = -0.5d0
      wo(3) = -0.5d0

      call cntsrf(vo,r,coo,ianz,isurf,iatoms)

      return
      end

      subroutine cntsrf(vec,r,coo,ianz,isurf,iatoms)
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      dimension vec(3), r(3), coo(3,*), isurf(*),ianz(*)

      do i=1,3
         vec(i) = 0.0d0
         r(i) = 0.0d0
      end do

      if (iatoms.le.0) return
      
      itel = 0
      do i=1,iatoms
         if (isurf(i).ne.0.and.ianz(i).ne.100) then
            itel = itel + 1
            do j=1,3
               vec(j) = vec(j) + coo(j,i)
            end do
         endif
      end do
 
      do j=1,3
         vec(j) = vec(j) / dble(itel)
      end do

      do i=1,iatoms
         if (isurf(i).ne.0.and.ianz(i).ne.100) then
            do j=1,3
               rmax = dabs(vec(j) - coo(j,i))
               if (rmax.gt.r(j)) r(j) = rmax
            end do
         endif
      end do

      do i=1,3
         r(i) = r(i)*2.0d0 + 6.0d0
      end do

      return
      end

      subroutine rotbck(w1,w2,w3,wn)
      implicit double precision (a-h,o-z)
      common /comsrf/  vo(3), r(3),v1(3),v2(3),v3(3), wo(3),
     &                 sl(3),isl
      dimension wn(3)

      do i=1,3
         wn(i) = v1(i)*(wo(1)+w1)*r(1) + 
     &           v2(i)*(wo(2)+w2)*r(2) +
     &           v3(i)*(wo(3)+w3)*r(3) + vo(i)
      end do

      return
      end

      subroutine rttbck(w1,w2,w3,wn)
      implicit double precision (a-h,o-z)
      common /comsrf/  vo(3), r(3),v1(3),v2(3),v3(3), wo(3),
     &                 sl(3),isl
      dimension wn(3)

      do i=1,3
         wn(i) = v1(i)*(-0.5d0+w1)*r(1) + 
     &           v2(i)*(-0.5d0+w2)*r(2) +
     &           v3(i)*(w3)*r(3) + vo(i)
      end do

      return
      end

      subroutine rtgbck(w1,w2,w3,wn)
      implicit double precision (a-h,o-z)
      common /comsrf/  vo(3), r(3),v1(3),v2(3),v3(3), wo(3),
     &                 sl(3),isl
      dimension wn(3)

      do i=1,3
         wn(i) = v1(i)*w1*sl(i) + 
     &           v2(i)*w2*sl(i) +
     &           v3(i)*w3*sl(i)
      end do

      return
      end

      subroutine cvtcom
      implicit double precision (a-h,o-z)
      common /comsrf/ vo(3),rt(3),v1t(3),v2t(3),v3t(3),wo(3),
     &                sl(3),isl
      common /plane/  px, py, pz, cx, cy, cz, r(3),v1(3),v2(3),iplat

      vo(1) = px
      vo(2) = py
      vo(3) = pz
      v3t(1) = -cx
      v3t(2) = -cy
      v3t(3) = -cz
      wo(1) = -0.5d0
      wo(2) = -0.5d0
      wo(3) = -0.5d0


      do i=1,3
         v1t(i) = v2(i)
         v2t(i) = v1(i)
      end do

      sl(1) = 1.0d0
      sl(2) = 1.0d0
      sl(3) = 1.0d0

      rt(1) = r(2)
      rt(2) = r(1)
      rt(3) = r(3)

      call vsc1(v1t,1.0d0,1.0d-4)
      call vsc1(v2t,1.0d0,1.0d-4)
      call vsc1(v3t,1.0d0,1.0d-4)

      return
      end
