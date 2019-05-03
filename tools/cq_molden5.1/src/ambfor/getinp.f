      subroutine getind(igtinp,
     &                         iconn,iopt,iresid,ityp,coo,q)
      implicit real (a-h,o-z), integer ( i-n)
      integer getlin
      character*137 str
      character*2 catom, catomt, tolowf,iel
      parameter (mxcon=10)
      parameter (numres=50000)
      parameter (mxion=2000)
      common /athlp/ iatoms, mxnat
      integer*2 ityp
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      common /opts/ idebug,imon,iarc,ilog,iout,osingl,dolbfgs,oqscal
      character*137 line
      common /curlin/ line
      logical box,cell,fast,adjatm
      common /pbc/ abc(3),abc2(3),angles(3),box,cell,fast
      common /residu/  qtot,ihsres,
     &                 ires(numres),ibeg(numres),iend(numres)
      common /h2oer/ numwat,natnow,nion,iontyp,ionpl,iw(mxion),niwat
      logical dall,osingl,gnreal,dolbfgs,oqscal
      dimension coo(3,*),iconn(mxcon+1,*),q(*),ityp(*),iopt(*)
      dimension iresid(*)
      dimension dmn(3),dmx(3),ddif(3)

      big = 1000000.0e0
      small = -1000000.0e0

      do i=1,3
         dmn(i) = big
         dmx(i) = small
      end do

c      idebug = 0
      iform = 0
      igtinp = 1
      dall = .false.
      ifndhn = 0

      if (idebug.eq.1) write(iun5,*) 'subroutine getinp'

c get number of atoms

      line = ''
      if (getlin(0).eq.1) then

         if (icdex(line,'[AMBFOR]').ne.0) then
            iform = 1
            ic = icdex(line,'box')
            if (ic.ne.0) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)
               box = .true.
               if (gnreal(abc,3,.false.)) then
                  do i=1,3
                     abc2(i) = 0.5e0*abc(i)
                  end do
               else
                  box = .false.
               endif
            endif
            
            ic = icdex(line,'cell')
            if (ic.ne.0) then
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)
               cell = .true.
               if (gnreal(abc,3,.false.)) then
                  do i=1,3
                     abc2(i) = 0.5e0*abc(i)
                  end do
                  if (gnreal(angles,3,.false.)) then
                     call setop(abc(1),abc(2),abc(3),
     &                          angles(1),angles(2),angles(3))
                  else
                     cell = .false.
                  endif
               else
                  cell = .false.
               endif
            endif
            
            if (getlin(0).eq.1) idum = 1
         endif

         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            iatoms = itype
            if (iatoms.gt.mxnat) then
               dall = .true.
               goto 100
            endif
         else
            goto 100
         endif
      else
         goto 100
      endif

110   continue

      do iat=1,iatoms
         if (getlin(0).eq.1) then
           ktype = nxtwrd(str,nstr,itype,rtype)
           if (ktype.ne.2) goto 100
           if (iform.eq.1) then
               iopt(iat) = itype
           else
               iopt(iat) = 1
           endif
           ktype = nxtwrd(str,nstr,itype,rtype)
           if (ktype.eq.1) then
                 if (nstr.eq.1) then
                    catomt(1:1) = str(1:1)
                    catomt(2:2) = ' '
                 else
                    catomt = str(1:2)
                    if (catomt(2:2).eq.'+'.or.catomt(2:2).eq.'-'
     &                 .or.catomt(2:2).eq.'*') catomt(2:2) = ' '
                 endif
                 catom = tolowf(catomt)
           else
              goto 100
           endif

           do i=1,3
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.3) then
                 coo(i,iat) = rtype
                 if (box) then
                    if (rtype.lt.dmn(i)) dmn(i) = rtype
                    if (rtype.gt.dmx(i)) dmx(i) = rtype
                 endif
              else
                 goto 100
              endif
           end do

           ktype = nxtwrd(str,nstr,itype,rtype)
           if (ktype.eq.2) then

c for atom types 2001-2011 (tinker convention) -> 649-659
              if (itype.gt.2000) itype = itype - 1352

              ityp(iat) = itype

              if (itype.lt.0) then

c this is a gaff atom type, so read in an extra charge field                 

                 ktype = nxtwrd(str,nstr,itype,rtype)
                 if (ktype.eq.3) then
                    q(iat) = rtype
                 else
                    goto 100
                 endif
              endif

           else
              goto 100
           endif
           nc = 0
           do i=1,mxcon
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (ktype.eq.2) then
                 iconn(i+1,iat) = itype
                 nc = nc + 1
              endif
           end do
           iconn(1,iat) = nc
        else
          goto 100
        endif
      end do

      natnow = 0
      numwat = 0
      if (box) then

        iat = 0
        do while (iat.lt.iatoms)
          iat = iat + 1
          it1 = ityp(iat)
          it2 = ityp(iat+1)
          it3 = ityp(iat+2)
          if (it1.eq.649.and.it2.eq.650.and.it3.eq.650) then
              if (natnow.eq.0) natnow = iat
              iat = iat + 2
          else
              natnow = 0
          endif
        end do

        if (natnow.ne.0) natnow = natnow - 1

        numwat = (iatoms - natnow)/3

      endif
      
      
      if (getlin(0).eq.1) then
         if (icdex(line,'[RESIDUES]').ne.0) then
             ihsres = 0
             do while (getlin(0).eq.1)
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.2) then
                   ihsres = ihsres + 1
                   irs = itype
                   ires(ihsres) = irs
                endif
                ktype = nxtwrd(str,nstr,itype,rtype)
                if (ktype.eq.2) then
                   istrt = itype
                   ibeg(ihsres) = istrt
                   if (ihsres.gt.1) iend(ihsres-1) = istrt - 1
                endif
                do i=istrt,iatoms
                   iresid(i) = irs
                end do
             end do
             iend(ihsres) = iatoms
         endif
      endif

      if (box) then
         adjatm = .false.
         do i=1,3
            ddif(i) = dmx(i) - dmn(i)
            if (ddif(i).gt.abc(i)) adjatm = .true.
         end do
         if (adjatm) then
            if (idebug.eq.1) then
               print*,'adjusted atom coordinates to fit in box:'
               print*,''
               print*,'min x,y,z ',(dmn(i),i=1,3)
               print*,'max x,y,z ',(dmx(i),i=1,3)
               print*,'max-min x,y,z ',(dmx(i)-dmn(i),i=1,3)
               print*,'abc ',(abc(i),i=1,3)
               print*,''
            endif
            call appbnd(coo,ityp)
         endif
      endif

      return

100   if (dall) then
         rewind iun2
         igtinp = -1
         return
      endif

      igtinp = 0

      return
      end

      integer function getlin(ieq)
      character*137 line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      common /curlin/ line
      common /mflin/  linmf

      getlin = 1
      read(iun2,'(a)',err=100,end=100) line
      linmf = linmf + 1
      do i=1,137
         if (ichar(line(i:i)).eq.9) line(i:i) = ' '
      end do

      if (ieq.eq.1.or.ieq.eq.2) then
         do i=1,137
            if (ichar(line(i:i)).eq.61) line(i:i) = ' '
         end do
      endif

      if (ieq.eq.2) then
         do i=1,137
            ii = ichar(line(i:i))
            if (ii.eq.40.or.ii.eq.41.or.ii.eq.34) line(i:i) = ' '
         end do
      endif

      return
100   getlin = 0
      return
      end

      integer function nxtwrd(string,strlen,itype,rtype)
c
c      string           nxtwrd = 1 
c      integer          nxtwrd = 2 
c      real             nxtwrd = 3 
c      no word          nxtwrd = 0
c
      character*(*) string
      integer itype,strlen
      real rtype
      real reada
      logical chkstr
      character*137 line
      common /curlin/ line

      nxtwrd = 0
      
      llen = linlen(line)
      if (llen.eq.0) return

      do while (line(1:1).eq.' ')
          line = line(2:)
      end do

      iend = index(line,' ')
      if (iend.eq.0) then
         iend = llen
      else
         iend = iend - 1
      endif
      if (chkstr(line,iend)) then
           nxtwrd = 1
           string = line(1:iend)
           strlen = iend
      elseif (index(line(1:iend),'.').ne.0) then
           nxtwrd = 3
           rtype = reada(line,1,iend)
      else
           nxtwrd = 2
           itype = reada(line,1,iend)
      endif

      line = line(iend+1:)

      return
      end

      logical function chkstr(line,iend)
      character*(*) line
      chkstr = .false.

      ie    = ichar('e')
      iee   = ichar('E')
      id    = ichar('d')
      idd   = ichar('D')
      nine  = ichar('9')
      izero = ichar('0')
      minus = ichar('-')
      iplus = ichar('+')
      idot  = ichar('.')
      icomma = ichar(',')
      islash  = ichar('/')

      ihase = 0
      idig = 0
      do i=1,iend
         n = ichar(line(i:i))
         if ((n.eq.ie.or.n.eq.iee.or.n.eq.id.or.n.eq.idd)
     &      .and.ihase.eq.0.and.idig.eq.1) then
             n = izero
             ihase = 1
         endif
         if (n.lt.iplus.or.n.gt.nine.or.n.eq.islash
     &       .or.n.eq.icomma) goto 100
         idig = 1
      end do

      n = ichar(line(1:1))
      n2 = ichar(line(2:2))
      if (iend.eq.1) then
         if (n.eq.minus) goto 100
         if (n.eq.iplus) goto 100
         if (n.eq.ie.or.n.eq.iee) goto 100
         if (n.eq.id.or.n.eq.idd) goto 100
      elseif (iend.gt.1) then
         if (n.eq.minus.and.n2.eq.minus) goto 100
      endif

      return
100   chkstr = .true.
      return
      end

      integer function linlen(line)
      character*(*) line

      linlen = 0

      do i=len(line),1,-1
         n = ichar(line(i:i))
         if (n.gt.32.and.n.le.126) goto 100
      end do

      return
100   linlen = i
      return
      end

      real function reada(a,istart,ll)
      implicit real (a-h,o-z)
      character*1 a(ll)
      logical multi,scient

      nine  = ichar('9')
      izero = ichar('0')
      minus = ichar('-')
      iplus = ichar('+')
      idot  = ichar('.')
      ie    = ichar('E')
      iee   = ichar('e')
      id    = ichar('D')
      idd   = ichar('d')
      idig  = 0
      k1    = 0
      k2    = 0
      k3    = 0
      one   = 1.d0
      x     = 1.d0
      do j=istart,ll
         n = ichar(a(j))
         if((n.le.nine.and.n.ge.izero).or. n.eq.minus.or.n.eq.idot)
     &   goto 7
      end do
      reada = 0.d0
      return

7     continue
c     before the dot
      do i=j,ll
         n = ichar(a(i))
         if (n.le.nine.and.n.ge.izero) then
            idig = idig + 1
            if (idig.gt.10) then
               ii = i
               goto 90
            endif
            k1 = k1*10 + n - izero
         elseif (n.eq.minus.and.i.eq.j) then
             one = -1.d0
         elseif (n.eq.idot) then
             goto 30
         else
             ii = i
             goto 90
         endif
      end do
30    continue
c     after the dot
      idig = 0
      do ii=i+1,ll
         n = ichar(a(ii))
         if (n.le.nine.and.n.ge.izero) then
            idig = idig + 1
            if (idig.gt.9) goto 90
            k2 = k2*10 + n - izero
            x = x /10
         elseif (n.eq.minus.and.ii.eq.i) then
            x = -x
         else
            goto 90
         endif
      end do
90    continue
      scient = .false.
      do jj=ii,ll
         n = ichar(a(jj))
         if (n.eq.ie.or.n.eq.iee.or.n.eq.id.or.n.eq.idd) goto 95
      end do
      goto 100
95    continue
         if (jj.lt.159) then
            scient = .true.
            n = ichar(a(jj+1))
            if (n.eq.minus.or.n.eq.iplus) then
                multi = .true.
                if (n.eq.minus) multi = .false.
                idig = 0
                do j=jj+2,ll
                   n = ichar(a(j))
                   if (n.le.nine.and.n.ge.izero) then
                       idig = idig + 1
                       if (idig.gt.9) goto 90
                       k3 = k3*10 + n - izero
                   else
                      goto 100
                   endif
                end do
            endif
         endif

100   continue
      reada = one * ( k1 + k2 * x)
      if (scient) then
          if (multi) then
              reada = reada * 10**k3
          else
              if (k3.gt.25) then
                 reada = 0.0e0
              else
                 reada = reada / 10**k3
              endif
          endif
      endif

      return
      end

      logical function gnreal(r,n,doget)
      implicit real (a-h,o-z)
      character*137 line,str
      common /curlin/ line
      integer getlin
      logical doget
      dimension r(*)

      gnreal = .true.

      if (doget) then
         if (getlin(0).ne.1) gnreal = .false.
      endif

      if (gnreal) then
          do i=1,n
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.3) then
                r(i) = rtype
             elseif (ktype.eq.2) then
                r(i) = dble(itype)
             else
                gnreal = .false.
             endif
          end do
      endif

      return
      end
    
      subroutine getpar
      implicit real (a-h,o-z), integer ( i-n)
      integer getlin
      character*137 str
      character*137 line
      common /curlin/ line
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      parameter (mxgff=72)
      parameter (mxamb=1590)
      integer amb2gf
      common /gftyp/  gfvdw(2,mxgff),amb2gf(mxamb)
      logical opfil,opar

      iun2t = iun2
      opar = .false.
      if (opfil(62,'param',5,1,1,1)) then
         iun2 = 62
         if (getlin(0).eq.1) then
            if (icdex(line,'[PARAM]').ne.0) opar = .true.
         endif
         if (opar) then
            do while (.true.)
                if (getlin(0).eq.1) then
                   it = 0
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.eq.2) then
                      it = itype
                   else
                      goto 100
                   endif
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.eq.3) then
                      par1 = rtype
                   else
                      goto 100
                   endif
                   ktype = nxtwrd(str,nstr,itype,rtype)
                   if (ktype.eq.3) then
                      par2 = rtype
                   else
                      goto 100
                   endif
                   if (it.gt.0.and.it.le.mxgff) then
                      print*,'old ityp ',it,gfvdw(1,it),gfvdw(2,it)
                      gfvdw(1,it) = par1
                      gfvdw(2,it) = par2
                      print*,'new ityp ',it,gfvdw(1,it),gfvdw(2,it)
                   endif
                else
                   goto 100
                endif
            end do
         endif
      else
         return
      endif

      iun2 = iun2t
      return

100   close(62)
      iun2 = iun2t

      return
      end
    
