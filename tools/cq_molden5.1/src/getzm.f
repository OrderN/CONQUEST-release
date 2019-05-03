      subroutine getzd(natoms,izo,ido,istatz,cc,ianc,
     &              bl,alph,bet,ibl,ialph,ibet,imap,ianz,iz,
     &              c,cz,alpha,beta,ian)

c this is really getzm

      implicit double precision (a-h,o-z), integer ( i-n)
      integer getlin
      character*137 str
      character*2 catom,catomt,tolowf,iel
      parameter (maxat=1000)
      parameter (maxat3=maxat*3)
      logical oparse,ottest,oerror,ovar
      common /zmfrst/ ihaszm, nz, mxzat
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      dimension ian(*),c(3,*),cz(3,*),alpha(*),beta(*)
      dimension ianc(*),cc(3,*)
      dimension bl(*),alph(*),bet(*),ibl(*),ialph(*),ibet(*),
     &          imap(*),ianz(*),iz(4,*)

      dimension iel(100)
      character*10 varnam,vartmp
      character*1 stmp
      dimension varnam(maxat3),ivar(maxat3),rvar(maxat3),ilnk(maxat3)
      character*137 line
      common /curlin/ line
      data iel/'bq',
     &         'h ', 'he',
     &         'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne',
     &         'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar',
     &         'k ', 'ca',
     &                     'sc', 'ti', 'v ', 'cr', 'mn',
     &                     'fe', 'co', 'ni', 'cu', 'zn',
     &                     'ga', 'ge', 'as', 'se', 'br', 'kr',
     & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     & 'in','sn','sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     & 'pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     & 'ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po',
     & 'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm',
     & 'bk','cf','x '/

      istatz = 0
      oparse = .false.
      maxnz = mxzat
      ottest=.true.
      if (izo.eq.0) then
         do i=1,3
            do j=1,4
             iz(j,i) = 0
            end do
         end do
         nz = 0
      endif
      ikeyzm = 0

      nzz = nz

100   continue
      ielins = 0
      if (izo.eq.1) then
         nz = nzz
      else
         nz = 0
      endif
      nzzt = 0
      izmat = 0
      do while (getlin(0).eq.1)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.1) then
           if ((str(1:4).eq.'vari'.or.str(1:4).eq.'VARI'
     &      .or.str(1:4).eq.'cons'.or.str(1:4).eq.'CONS'
     &      .or.str(1:4).eq.'end'.or.str(1:4).eq.'END').and.izmat.eq.1)
     &     then
             ielins = 1
             goto 200
           endif
           if (nstr.le.10.and.
     &     (str(1:4).eq.'zmat'.or.str(1:4).eq.'ZMAT')) then
             nz = nzz
             nzzt = 0
             ikeyzm = 1
             izmat = 1
             goto 150
           endif
           if (nstr.le.10) then
              nz = nz + 1
              nzzt = nzzt + 1
              numv = 3
              if (nzzt.le.3) numv = nzzt - 1
c
c Atom String
c
              if (oparse) then
                 if (nstr.eq.1) then
                    catomt(1:1) = str(1:1)
                    catomt(2:2) = ' '
                 else
                    catomt = str(1:2)
                    if (catomt(2:2).eq.'-') catomt(2:2) = ' '
                 endif
                 catom = tolowf(catomt)
                 do j=1,100
                    if (catom .eq. iel(j)) ianz(nz) = j - 1
                 end do
                 if (catom.eq.'xx') ianz(nz) = 99
              endif
              do j=1,numv
c
c Connectivity
c
                 if (nxtwrd(str,nstr,itype,rtype).ne.2) then
                    if (nzzt.eq.2) backspace iun2
                    goto 100
                 else
                    if (oparse) iz(j,nz) = itype
                 endif
c
c Variable
c
                 ktype = nxtwrd(str,nstr,itype,rtype)
                 if (ktype.eq.0) then
                    if (nzzt.eq.2) backspace iun2
                    goto 100
                 elseif (ktype.eq.1.and.nstr.gt.10) then
                    if (nzzt.eq.2) backspace iun2
                    goto 100
                 elseif (ktype.eq.1.and.oparse) then
                    ovar = .false.
                    do k=1,nvars
                       if (str(1:nstr).eq.varnam(k)(1:ivar(k))) then
                          tmpvar = rvar(k)
                          if (ilnk(k).le.0) then
                             ivart = ilnk(k)+1
                             ilnk(k) = nz
                          else
                             ivart = ilnk(k)
                          endif
                          ovar = .true.
                       endif
                       if (nstr.gt.1) then
                          if (str(2:nstr).eq.varnam(k)(1:ivar(k)).and.
     &                        str(1:1).eq.'-') then
                             tmpvar = -rvar(k)
                             if (ilnk(k).le.0) then
                                ivart = ilnk(k)+1
                                ilnk(k) = nz
                             else
                                ivart = -ilnk(k)
                             endif
                             ovar = .true.
                          endif
                          if (str(2:nstr).eq.varnam(k)(1:ivar(k)).and.
     &                        str(1:1).eq.'+') then
                             tmpvar = rvar(k)
                             if (ilnk(k).le.0) then
                                ivart = ilnk(k)+1
                                ilnk(k) = nz
                             else
                                ivart = ilnk(k)
                             endif
                             ovar = .true.
                          endif
                       endif
                    end do
                    if (.not.ovar.and.nvars.ne.0) then
c                       if (nz.le.2) return
                       if (nzzt.le.2) then
                           if (nzzt.eq.2) backspace iun2
                           goto 100
                       endif
                       call inferr(
     &                   'Gauss/Gamess: Missing Variable Name !',1)
                       call inferr(
     &                   'Variable: '//str(1:nstr),1)
                       return
                    endif
                 elseif (ktype.eq.2.and.oparse) then
                     tmpvar = 1.0d0*itype 
                     ivart = 0
                 elseif (ktype.eq.3.and.oparse) then
                     tmpvar = rtype 
                     ivart = 0
                 endif
                 if (oparse) then
                    if (j.eq.1) then
                       bl(nz) = tmpvar
                       if (izo.eq.1) then
                          ibl(nz) = 1
                       else
                          ibl(nz) = ivart
                       endif
                    elseif (j.eq.2) then
                       alph(nz) = tmpvar 
                       if (izo.eq.1) then
                          ialph(nz) = 1
                       else
                          ialph(nz) = ivart
                       endif
                    elseif (j.eq.3) then
                       bet(nz) = tmpvar
                       if (izo.eq.1) then
                          ibet(nz) = 1
                       else
                          ibet(nz) = ivart
                       endif
                    endif
                 endif
              end do
c
c Check for Gamess ITYPE
c
              ktype = nxtwrd(str,nstr,itype,rtype)
              if (oparse) iz(4,nz) = 0
              if (ktype.ne.0) then
                 if (ktype.eq.2) then
                     if (oparse.and.(abs(itype).eq.1.or.itype.eq.0)) 
     &                   iz(4,nz) = itype
                 elseif (ktype.eq.3) then
                     if (oparse.and.abs(rtype).eq.1.0d0)
     &                   iz(4,nz) = int(rtype)
                 elseif (ktype.eq.1.and.nstr.eq.1) then
                     stmp = str(1:1)
                     if (stmp.eq.'L'.or.stmp.eq.'l'.or.
     &                   stmp.eq.'H'.or.stmp.eq.'h'.or.
     &                   stmp.eq.'M'.or.stmp.eq.'m') then
                     else
                        if (nzzt.eq.2) backspace iun2
                        goto 100
                     endif
                 else
                     if (nzzt.eq.2) backspace iun2
                     goto 100
                 endif
              endif
           endif
c
c Empty Line
c
         elseif (ktype.eq.0.and.nzzt.gt.1) then
           ielins = 1
           goto 200
         else
c
c First word is not a string
c
           goto 100
         endif
150      continue
      end do
c
c Out of Lines, Didnt Find Zmat
c
      
      return

200   continue
c
c Found Zmat
c
      if (oparse) goto 500

c
c Get Variables Constants if any
c
300   if (ielins.eq.1) then
        nvars = 0
        do while (getlin(1).eq.1)
          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.1) then
            if (str(1:4).eq.'cons'.or.str(1:4).eq.'CONS') then
               ielins = 2
            elseif (str(1:3).eq.'end'.or.str(1:3).eq.'END') then
               goto 400
            else
               if (nstr.gt.10) goto 400
               vartmp = str
               ntmp = nstr
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.ne.2.and.ktype.ne.3) goto 400
               nvars = nvars + 1
               varnam(nvars) = vartmp
               ivar(nvars) = ntmp
               ilnk(nvars) = 0
               if (ielins.eq.2) ilnk(nvars) = -1
               if (ktype.eq.2) rvar(nvars) = 1.0d0*itype
               if (ktype.eq.3) rvar(nvars) = rtype
            endif
          elseif (ktype.eq.0) then
            ielins = ielins + 1
            if (ielins.eq.3) goto 400
          else
             goto 400
          endif
        end do
      endif

400   if (ido.eq.1) then
          icnt = 0
          istmp = 0
          do while (getlin(1).eq.1)
             ktype = nxtwrd(str,nstr,itype,rtype)
             if (ktype.eq.1) then
                if (nstr.ge.3) then
                   if (icdex(str,'map').ne.0) istmp = 1
                endif
             elseif (ktype.eq.2.and.istmp.eq.1) then
                do i=1,icnt
                   if (itype.eq.imap(i)) then
                      call inferr('Two map numbers are equal!!',0)
                      return
                   endif
                end do
                icnt = icnt + 1
                imap(icnt) = itype
             endif
          end do
          if (icnt.ne.nz) then
              call inferr('Not enough map numbers !!',0)
              return
          endif
      endif

      if (ikeyzm.eq.1) then
         call scback(line,'zmat',istat)
      else
         rewind iun2
      endif
      
      oparse = .true.
      goto 100

500   if (nz.le.2) return
      istatz = 1
c      do i=1,nz
c         print*,ianz(i),bl(i),alph(i),bet(i),(iz(j,i),j=1,4)
c      end do
      if (izo.eq.1.or.ido.eq.1) return
      call stoc(maxnz,nz,0,0,0,ianz,iz,bl,alph,bet,
     &       ottest,natoms,ian,c,cz,imap,alpha,beta,
     &       oerror,.true.,.true.)
      if (oerror) then
         istatz = 0
      else
         ihaszm = 1
         do i=1,natoms
            do j=1,3
               cc(j,i) = c(j,i)
            end do
            ianc(i) = ian(i)
         end do
      endif

      return
      end

