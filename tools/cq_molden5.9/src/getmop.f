      subroutine getmdp(natoms,heat,izo,ido,istatz,cc,ianc,
     &              bl,alph,bet,ibl,ialph,ibet,imap,ianz,iz,
     &              c,cz,alpha,beta,ian)

c this is really getmop

      implicit double precision (a-h,o-z), integer ( i-n)
      integer getlin
      character*137 string
      character*2 catom,catomt,tolowf,iel
      logical first,ottest, oerror,zmstar,getnxt,gnreal

      common /zmfrst/ ihaszm, nz, mxzat

      dimension ian(*),c(3,*),cz(3,*),alpha(*),beta(*)
      dimension ianc(*),cc(3,*)
      dimension bl(*),alph(*),bet(*),ibl(*),ialph(*),ibet(*),
     &          imap(*),ianz(*),iz(4,*)
      dimension izt(4),iztt(3,3),zzz(3),idm(3),tmp(3)
      dimension iel(100)
      character*137 coplin
      character*137 line
      common /curlin/ line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
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
      ifnd = 0
      irwd = 0
      toang = 0.52917706d0

5     continue

      first = .true.
      if (izo.eq.1) then
          nzz = nz
      else
          nz = 0
      endif
      maxnz = mxzat
      ottest=.true.
      heat = 0.0d0

10    continue
      zmstar = .false.
      nztt = 0
      do while (getlin(0).eq.1)
          if (icdex(line,'[AMB').ne.0) return
          i1 = index(line,'(')
          i2 = index(line,')')
          if (i1.ne.0.and.i2.ne.0.and.i2.gt.i1) then
               line = line(1:i1-1)//line(i2+1:)
          endif
          if (index(line,'VARIABLE        FUNCTION').ne.0) then
               nword = 3
               if (index(line,'SECOND').ne.0) nword = 4
               if (getlin(0).eq.0) goto 10
               do ii=1,nword
                   ktype = nxtwrd(string,nstr,itype,rtype)
               end do
               if (ktype.eq.3) then
                   heat = rtype
               endif
               goto 10
          endif
          if (index(line,'HEAT OF FORMATION').ne.0) then
             do while(.true.)
                ktype = nxtwrd(string,nstr,itype,rtype)
                if (ktype.eq.3) heat = rtype
                if (ktype.eq.0.or.ktype.eq.3) goto 10
             end do
          endif
          if (index(line,'POINT   POTENTIAL').ne.0.or.
     &        index(line,'POINT  POTENTIAL').ne.0) then
               in = 2
               if (index(line,'FEMTOSECONDS').ne.0) in = 3
               if (getlin(0).eq.0) goto 10
c               do ii=1,4
c Armando Juan Navarro Vazquez says total energy not interresting
c Use potential energy
c
               do ii=1,in
                   ktype = nxtwrd(string,nstr,itype,rtype)
               end do
               if (ktype.eq.3) heat = rtype
               goto 10
          endif
          if (index(line,'FINAL GEOMETRY OBTAINED').ne.0) then
             call redel(line,4)
             coplin = line
             iws = 0
             ktype = -1
             do while (ktype.ne.0)
                ktype = nxtwrd(string,nstr,itype,rtype)
                if (ktype.ne.0) iws = iws + 1
             end do
             line = coplin
             if (iws.eq.7) then
               natoms = 0
               do while(.true.)
                 ktype = nxtwrd(string,nstr,itype,rtype)
                 if (ktype.eq.1.and.(nstr.eq.1.or.nstr.eq.2)) then
                      if (nstr.eq.1) then
                         catomt(1:1) = string(1:1)
                         catomt(2:2) = ' '
                      else
                         catomt = string(1:2)
                      endif
                 else
                      if (natoms.ne.0) then
                         istatz = 2
                         return
                      endif
                      goto 100
                 endif

                 iatmp = 0
                 catom = tolowf(catomt)
                 do i=1,100
                    if ( catom .eq. iel(i) ) iatmp = i - 1
                 end do
                 if (catom.eq.'xx') iatmp = 99
                 if (catom.eq.'tv') iatmp = 99
                 if (iatmp.eq.0.or.iatmp.gt.100) goto 100

                 natoms = natoms + 1
                 ianc(natoms) = iatmp

                 do ii=1,3
                    ktype = nxtwrd(string,nstr,itype,rtype)
                    if (ktype.eq.3) then
                       cc(ii,natoms) = rtype / toang
                    else
                       goto 100
                    endif
                    ktype = nxtwrd(string,nstr,itype,rtype)
                    if (ktype.ne.2)  goto 100
                 end do
                 if (getlin(0).ne.1) return
               end do
             endif
          endif
          if (index(line,'scf done: ').ne.0) then
             do while(.true.)
                ktype = nxtwrd(string,nstr,itype,rtype)
                if (ktype.eq.3) heat = rtype
                if (ktype.eq.0.or.ktype.eq.3) goto 10
             end do
          endif
          do i=1,3
             idm(i) = 0
             zzz(i) = 0.0d0
             izt(i) = 0
          end do
          ktype = nxtwrd(string,nstr,itype,rtype)
          if (ktype.eq.2) then
             if (itype.eq.1.and.first) zmstar = .true.
             if (zmstar) then
                 ktype = nxtwrd(string,nstr,itype,rtype)
             else
                 goto 100
             endif
          endif
c element label
          if (ktype.eq.1.and.(nstr.eq.1.or.nstr.eq.2)) then
               if (nstr.eq.1) then
                  catomt(1:1) = string(1:1)
                  catomt(2:2) = ' '
               else
                  catomt = string(1:2)
               endif
          else
               goto 100
          endif

c check which zmstar type

          if (zmstar.and.nztt.eq.0) then
               izmtyp = 0
               if (gnreal(tmp,3,.false.)) izmtyp = 1
          endif

          num = nztt
          if (num.gt.3) num = 3
          if (.not.zmstar.or.(zmstar.and.izmtyp.eq.1)) num = 3

          getnxt = .true.
          do i=1,num
             if (getnxt) ktype = nxtwrd(string,nstr,itype,rtype)

c bl,aplh,bet

             getnxt = .true.
             if (ktype.eq.2.or.ktype.eq.3) then
                if (ktype.eq.2) then
                   zzz(i) = 1.0d0*itype
                else
                   zzz(i) = rtype
                endif
             else
                if (zmstar) then
                   if (ktype.eq.0) then
c                     read in connectivity
                      goto 200
                   else
                      goto 100
                   endif
                else
                   goto 100
                endif
             endif

c optimize flag 1-3

             ktype = nxtwrd(string,nstr,itype,rtype)
             if (zmstar) then
                if (ktype.eq.1) then
                   if (nstr.eq.1) then
                      if (string(1:1).eq.'*'.or.
     &                    string(1:1).eq.'+') then
                          idm(i) = 1
                      else
                          goto 100
                      endif
                   else
                      goto 100
                   endif
                else
                   getnxt = .false.
                   idm(i) = 0
                endif
             else
                if (ktype.eq.2) then
                   idm(i) = abs(itype)
                else
                   goto 100
                endif
             endif
          end do
200       continue
c
c connectivity
c
          num = nztt
          if (num.gt.3) num = 3
          if (.not.zmstar) num = 3

          do i=1,num
             if (getnxt) ktype = nxtwrd(string,nstr,itype,rtype)
             getnxt = .true.
             if (ktype.eq.2) then
                izt(i) = itype
             endif
          end do

c check first three atoms

          if (izo.ne.1) then
             ioke = 1
             if (nztt.le.2) then
                do i=1,3
                   iztt(i,nztt+1) = izt(i)
                end do
             endif
             if (nztt.le.2) then
                if (iztt(1,1).ne.0) ioke = 0
                if (iztt(2,1).ne.0) ioke = 0
                if (iztt(3,1).ne.0) ioke = 0
             endif
             if (nztt.ge.1.and.nztt.le.2) then
                if (iztt(1,2).le.0) ioke = 0
                if (iztt(2,2).ne.0) ioke = 0
                if (iztt(3,2).ne.0) ioke = 0
             endif
             if (nztt.eq.2) then
                if (iztt(1,3).le.0) ioke = 0
                if (iztt(2,3).le.0) ioke = 0
                if (iztt(3,3).ne.0) ioke = 0
             endif
             if (ioke.ne.1) then
                istatz = 0
                first = .true.
                goto 5
             endif
          endif

c check atom4 and upwards

          iztot = 0
          do i=1,3
             iztot = iztot + izt(i)
          end do
          if (nztt.gt.4.and.iztot.eq.0) then
             istatz = 0
             first = .true.
             goto 5
          endif

300       continue
          if (zmstar.and.heat.eq.0.0d0.and.ifnd.eq.0) then
             ifnd = 1
             goto 100
          endif
c
c copy temporary variables into z-mat
c
          iatmp = 0
          catom = tolowf(catomt)
          do i=1,100
             if ( catom .eq. iel(i) ) iatmp = i - 1
          end do
          if (catom.eq.'xx') iatmp = 99
          if (catom.eq.'tv') iatmp = 99
          if (iatmp.eq.0.or.iatmp.gt.100) goto 100
          first = .false.
          if (nztt.ge.1) istatz = 1
          nz = nz + 1
          nztt = nztt + 1
          ianz(nz) = iatmp
          do i=1,3
             iz(i,nz) = izt(i)
          end do

          iz(4,nz) = 0
          bl(nz) = zzz(1)
          ibl(nz) = idm(1)
          alph(nz) = zzz(2)
          ialph(nz) = idm(2)
          bet(nz) = zzz(3)
          ibet(nz) = idm(3)
      end do
      first = .false.

100   if (first) goto 10
      if (izo.eq.1) then
          if (istatz.eq.0) nz = nzz
          return
      endif
      if (istatz.eq.1) then
c          call prtzm
          ihaszm = 1
          call stocc(maxnz,nz,iz,bl,alph,bet,alpha,beta,ierror)
          if (ierror.lt.5) then
             call stoc(maxnz,nz,0,0,0,ianz,iz,bl,alph,bet,
     &       ottest,natoms,ian,c,cz,imap,alpha,beta,
     &       oerror,.true.,.true.)
             if (oerror) then
                istatz = 0
                return
             endif
          else
             istatz = 0
             return
          endif
          do i=1,natoms
             do j=1,3
                cc(j,i) = c(j,i)
             end do
             ianc(i) = ian(i)
          end do
      else
          if (ifnd.eq.1.and.irwd.eq.0.and.ido.eq.1) then
             irwd = 1
             call rewfil
             goto 5
          endif
      endif

      return
      end
