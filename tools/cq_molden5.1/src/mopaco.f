      subroutine mopacd(istats,mopopt,irtype,ianz,iatclr,iconn,coo,qat,
     &                  nnat,icent,nspg,ichx,nopr,ir,it,
     &                  xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxcon=10)
      common /athlp/  iatoms, mxnat
      common /coord/  xyz(3,numatm)
      common /charge/ dipo(3),ihasq,ihsdp,iqon,idipon
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      parameter (mxel=100)
      character*2 elemnt
      common /elem/elemnt(mxel)
      common /curlin/ line
      character*137 string,line,str
      character*2 eltmp,tstr
      integer getlin
      common /xyzopt/ ixyz,ipdbwh,iambch,nwramb
      integer*2 ir,it
      logical chkesp
      dimension rr(3,3),coo(3,*),qat(*),ianz(*),iatclr(*),
     &          iconn(mxcon+1,*),ir(3,3,192),it(3,192)

      istats = 1
      ihasq = 0
      todeg = 45.0d0 / datan(1.0d0)

      rewind iun2
      isol = 0
      call search(line,'THE SYSTEM IS A  SOLID',istat)
      if (istat.eq.1) isol = 1
      rewind iun2

      irtype = 0
      if (mopopt.eq.1) then
         call searchv(line,'MOLECULE IN FORCE CALCULATION',
     &         'INTRINSIC REACTION COORDINATE',
     &         'POINTS ON REACTION COORDINATE',
     &         'FIRST VARIABLE   SECOND',
     &         'NUMBER OF POINTS IN PATH',istat)
         if (istat.eq.1) then
             if (index(line,'REACTION').ne.0.or.
     &           index(line,'POINTS IN PATH').ne.0.or.
     &           index(line,'VARIABLE').ne.0) goto 200
             irtype = 4
             natoms = 0
             call readel(line,3)
             do while (getlin(0).eq.1)
                if (linlen(line).le.1) return
                natoms = natoms + 1
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.2) goto 200
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.1.and.ktype.ne.2) goto 200
                if (ktype.eq.1) then
                   nat(natoms) = iatnum(str,nstr)
                else
                   nat(natoms) = itype
                endif
                do i=1,3
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.ne.3) goto 200
                   xyz(i,natoms) = rtype
                end do
             end do
             return
         endif
      endif

      rewind iun2
      ifrst = 0


c       skip the initial coordinates

      call search(string,'          Cartesian Coordinates',istat)
      if (istat.ne.1) goto 200

c       read the optimised coordinates

      call search(string,'          Cartesian Coordinates',istat)
      if (istat.ne.1) then
        ifrst = 1
        rewind iun2
        call search(string,'          Cartesian Coordinates',istat)
      endif

      read(iun2,'(a)',end=200) string
      if (mopopt.lt.3) then
         read(iun2,'(a)',end=200) string
         read(iun2,'(a)',end=200) string
      endif


      i = 1
      if (mopopt.eq.1.or.mopopt.eq.3) then
         do while (.true.)
          read(iun2,'(a)',err=200) string
          ichs = ifblen(string)
          if (ifrst.eq.1) then
           if (ichs.ne.50) goto 100
           if (mopopt.eq.3) then
              read(string,'(14x,i2,4x,3(f10.4))',err=200)
     &        nat(i),(xyz(j,i),j=1,3)
           else
              read(string,'(14x,a2,4x,3(f10.4))',err=200)
     &        eltmp,(xyz(j,i),j=1,3)
              do j=1,99
                if (eltmp.eq.elemnt(j)) nat(i) = j
              end do
           endif
          else
           if (ichs.ne.60) goto 100
           read(string,'(14x,a2,14x,3(f10.4))',err=200)
     &     eltmp,(xyz(j,i),j=1,3)
           if (mopopt.eq.3) then
              tstr = eltmp(1:2)
              ii = ichar(tstr(2:2))
              if (ii.lt.65.or.ii.gt.122.or.(ii.lt.97.and.ii.gt.90))
     &        then
                 tstr(2:2) = tstr(1:1)
                 tstr(1:1) = ' '
              endif
              eltmp = tstr(1:2)
           endif
           do j=1,99
             if (eltmp.eq.elemnt(j)) nat(i) = j
           end do
          endif
           i = i + 1
         end do
      elseif (mopopt.eq.2) then
         do while (.true.)
           read(iun2,'(a)',err=200) string
           ichs = ifblen(string)
           if (ichs.lt.2) goto 100
            if (ichs.lt.36) then
              read(string,'(4x,i5,3(f8.4,1x))',err=200)
     &        nat(i),(xyz(j,i),j=1,3)
              i = i + 1
            else
              read(string,'(4x,i5,3(f8.4,1x),5x,i5,2x,3(f8.4,1x))',
     &        err=200)nat(i),(xyz(j,i),j=1,3),nat(i+1),
     &        (xyz(j,i+1),j=1,3)
              i = i + 2
            endif
         end do
      endif

100   continue
      natoms = i - 1

      if (mopopt.eq.1.and.isol.eq.1) then
         call scback(line,' Tv ',istat)
         if (istat.eq.1) then
             backspace iun2
             backspace iun2
             do ii=1,3
               if (getlin(0).eq.1) then
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   kk = 0
                   do jj=1,6
                      ktype = nxtwrd(str,nstr,itype,rtype)
                      if (ktype.eq.3.and.kk.lt.3) then
                         kk = kk + 1
                         rr(kk,ii) = rtype
                      endif
                   end do
               endif
             end do
             nspg = 1
             a = vlen(rr(1,1))
             b = vlen(rr(1,2))
             c = vlen(rr(1,3))
             call impsc(rr(1,2),rr(1,3),csa)
             call impsc(rr(1,1),rr(1,3),csb)
             call impsc(rr(1,1),rr(1,2),csc)
             alpha = dacos(csa)*todeg
             beta = dacos(csb)*todeg
             gamma = dacos(csc)*todeg
             call prcell(nspg,a,b,c,alpha,beta,gamma)
             call setop(xa,ya,yb,za,zb,zc,a,b,c,alpha,beta,gamma,1)
             call cprot(nspg,nopr,icent,ir,it,.false.)
             call xyzcoo(1,1,0)
             call doconn

             iatoms = natoms
             nnat = iatoms
             nstor = mxnat-iatoms
             do i=1,iatoms
                call crt2fr(xyz(1,i),coo(1,nstor+i),xa,ya,yb,za,zb,zc)
                ianz(nstor+i) = ianz(i)
                do j=1,iconn(1,i)+1
                   iconn(j,nstor+i) = iconn(j,i)
                end do
                iatclr(nstor+i) = iatclr(i)

             end do

c             call prop(nopr,ir,it)

             iftyp = 6
             ichx = 1
             call fdat(1,1,0,0,0,0)

         endif
      endif

      if (mopopt.eq.1) then
         chkesp = .false.
120      rewind iun2
         if (chkesp) then
            call search(string,'ELECTROSTATIC POTENTIAL CHARGES',
     &                  istat)
         else
            call search(string,'NET ATOMIC CHARGES AND',istat)
         endif
         if (istat.eq.1) then
             ihasq = 1
             qtot = 0.0d0
             call readel(line,2)
             do i=1,natoms
                if (getlin(0).ne.1) goto 150
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.2) goto 150
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.1) goto 150
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.ne.3) goto 150
                qat(i) = rtype
                qtot = qtot + qat(i)
             end do
c             write(iun3,*) 'total charge ',qtot
         endif
         if (.not.chkesp) then
            chkesp = .true.
            goto 120
         endif
      endif

      return
150   ihasq = 0
      return

200   istats = 0
      return
      end

      integer function iatnum(str,nstr)
      parameter (mxel=100)
      character*(*) str
      character*2 catom,tolowf,elemnt
      common /elem/elemnt(mxel)

      iatnum = 100
      if (nstr.eq.1) then
          catom(1:1) = ' '
          catom(2:2) = str(1:1)
      else
          catom = str(1:2)
      endif
      do j=1,mxel
          if (tolowf(catom) .eq. tolowf(elemnt(j))) iatnum = j
      end do

      return
      end
