      subroutine rdbad(idebug,dfree,istats,ityp)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      parameter (mxonh=100)
      character*3 cnam
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp,dfree
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)
      common /oniomh/ ioni,nion,ionih(mxonh),natonh(mxonh),
     &                iomap(numatm),icntat(mxonh),fct(mxonh),
     &                xyzi(3,mxonh)
      integer*2 ityp
      equivalence (nshell,i)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character*137 line,str
      common /curlin/ line
      integer getlin
      logical gnreal
      dimension r(4),ityp(*)


      istats = 1
      ig03 = 0

      if (idebug.eq.1)
     &    call inferr('looking for gaussian basis-set',0)

1030  format(4d18.10)
1040  format(4x,4d18.10)
      if (dfree.ne.0) then
         if (dfree.eq.1) then
            call searchu(line,'[GTO]',istat)
         elseif (dfree.eq.2) then
            call srchmf(line,'$BASIS',istat)
         endif
      else
         call seardu(line,'Basis set in the form',
     &                    'AO basis set:',istat)
         if (istat.eq.0) then
             call inferr(
     &       'Use keyword gfinput in gaussian inputfile',1)
             goto 1000
         endif
         if (index(line,'AO basis set').ne.0) ig03 = 1
         if (index(line,'AO basis set:').ne.0) ig03 = 2
      endif

      i = 0
      k = 0
      if (ioni.eq.1) then
         do m=1,numatm
            iomap(m) = 0
            ityp(m) = 20000
         end do
      endif
      mm = 1
      mmdf = 1
      jann = 0
      jmode = 1

10    if (dfree.ne.0) then
         if (getlin(0).eq.1) idum = 0
      else
         read(iun2,'(a)',err=1000) line
      endif

      if (dfree.ne.0) then

         ktype = nxtwrd(str,nstr,itype,rtype)

         if (ktype.eq.1) then
            if (nstr.ge.1) then
               if (index(str,'[').ne.0) return
               if (dfree.eq.2) then
                  if (icdex(str,'$END').ne.0) return
               endif
               if (str(1:1).eq.'#'.or.str(1:1).eq.'*') then
                  if (dfree.eq.2) jmode = 1
                  goto 10
               endif
               if (jmode.eq.1.and.dfree.eq.2) then
                  jann = jann + 1
                  jmode = 0
                  goto 10
               endif
            endif
         elseif (ktype.eq.2) then
            jann = itype
            goto 10
         elseif (ktype.eq.0) then
            goto 10
         endif

      else
         if (ig03.eq.2.and.line(2:5).ne.'Atom') return
         if (i.eq.0.or.line(2:5).eq.'****'.or.line(2:5).eq.'Atom') then
            if (line(2:5).eq.'Atom') then
               ig03 = 2
c example: Atom B3       Shell     8 SP   2    bf   20 -    23
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)
c            break up string and extract jann
               do kk=1,nstr
                  ik = ichar(str(kk:kk))
                  if (ik.lt.48.or.ik.gt.57) str(kk:kk) = ' '
               end do
               do kk=1,nstr
                  if (str(kk:kk).ne.' ') istrt = kk
               end do
               jann = reada(str,1,nstr)
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.1) then
                  cnam = '   '
                  if (nstr.gt.3) nstr = 3
                  cnam(1:nstr) = str(1:nstr)
                  call tocap(cnam,3)
               else
                  goto 1000
               endif
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.2) then
                  nsh = itype
               else
                  goto 1000
               endif
            else
               if (line(2:5).eq.'****') then
                   read(iun2,'(a)',err=1000) line
                   if (linlen(line).eq.0) return
               endif
               if (line(7:7).ne.' ') then
                  read(line,'(i7)',err=1000) jann
               elseif (line(4:4).ne.' ') then
                  read(line,'(i4)',err=1000) jann
               else
                  read(line,'(i3)',err=1000) jann
               endif
               if (ioni.eq.1) then
                  k = k + 1
                  iomap(k) = jann
                  nat(k) = nat(jann)
                  ityp(jann) = ityp(jann) - 20000
                  do ll=1,3
                     xyz(ll,k) = xyz(ll,jann)
                  end do
                  do jj=1,nion
                      if (jann.eq.ionih(jj)) then
                         nat(k) = natonh(jj)
                         do ll=1,3
                            xyz(ll,k) = xyzi(ll,jj)
                         end do
                      endif
                  end do
               endif
               read(iun2,'(a)',err=1000) line
            endif
         endif
      endif

      if (dfree.eq.0.and.ig03.eq.1) then
         ktype = nxtwrd(str,nstr,itype,rtype)
      endif

      if (dfree.ne.0.or.ig03.eq.1) then
      
         cnam = '   '
         if (nstr.gt.3) nstr = 3
         cnam(1:nstr) = str(1:nstr)
         call tocap(cnam,3)
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.2) then
            nsh = itype
         endif
      else
c         if (ig03.gt.0) then
c            read(line,'(1x,a3,1x,i3)',err=1000) cnam, nsh
c         else
c            read(line,'(1x,a3,1x,i2)',err=1000) cnam, nsh
c         endif
         if (ig03.eq.0) then
            read(line,'(1x,a3,1x,i2)',err=1000) cnam, nsh
         endif
      endif

      i = i + 1
      shelln(i)  = nsh
      if (ioni.eq.1) then
         jan(i)     = k
         gx(i)      = xyz(1,k)
         gy(i)      = xyz(2,k)
         gz(i)      = xyz(3,k)
      else
         jan(i)     = jann
         gx(i)      = xyz(1,jann)
         gy(i)      = xyz(2,jann)
         gz(i)      = xyz(3,jann)
      endif
      call settc(cnam,shellt(i),shellc(i))
      shella(i)  = mm
      shladf(i)  = mmdf
         do 20 j = 1, shelln(i)
           if (shellt(i).eq.0) then
               if (dfree.ne.0) then
                  if (gnreal(r,2,.true.)) then
                     exx(mm) = r(1)
                     c1(mm)  = r(2)
                  endif
               else
                  if (ig03.gt.0) then
                     read(iun2,1040,err=1000) exx(mm), c1(mm)
                  else
                     read(iun2,1030,err=1000) exx(mm), c1(mm)
                  endif
               endif
           endif
           if (shellt(i).eq.1.and.shellc(i).eq.1) then
               if (dfree.ne.0) then
                  if (gnreal(r,2,.true.)) then
                     exx(mm) = r(1)
                     c2(mm)  = r(2)
                  endif
               else
                  if (ig03.gt.0) then
                     read(iun2,1040,err=1000) exx(mm), c2(mm)
                  else
                     read(iun2,1030,err=1000) exx(mm), c2(mm)
                  endif
               endif
           endif
           if (shellt(i).eq.1.and.shellc(i).ne.1) then
               if (dfree.ne.0) then
                  if (gnreal(r,3,.true.)) then
                     exx(mm) = r(1)
                     c1(mm)  = r(2)
                     c2(mm)  = r(3)
                  endif
               else
                  if (ig03.gt.0) then
                     read(iun2,1040,err=1000) exx(mm), c1(mm), c2(mm)
                  else
                     read(iun2,1030,err=1000) exx(mm), c1(mm), c2(mm)
                  endif
               endif
           endif
           if (shellt(i).eq.2.and.shellc(i).eq.0) then
               if (dfree.ne.0) then
                  if (gnreal(r,4,.true.)) then
                     exx(mm) = r(1)
                     c1(mm)  = r(2)
                     c2(mm)  = r(3)
                     c3(mmdf)= r(4)
                  endif
               else
                  if (ig03.gt.0) then
                     read(iun2,1040,err=1000) 
     &                        exx(mm), c1(mm), c2(mm), c3(mmdf)
                  else
                     read(iun2,1030,err=1000) 
     &                        exx(mm), c1(mm), c2(mm), c3(mmdf)
                  endif
               endif
               mmdf = mmdf + 1
           endif
           if (shellt(i).eq.2.and.shellc(i).ne.0) then
               if (dfree.ne.0) then
                  if (gnreal(r,2,.true.)) then
                     exx(mm) = r(1)
                     c3(mmdf)= r(2)
                  endif
               else
                  if (ig03.gt.0) then
                     read(iun2,1040,err=1000) exx(mm), c3(mmdf)
                  else
                     read(iun2,1030,err=1000) exx(mm), c3(mmdf)
                  endif
               endif
               mmdf = mmdf + 1
           endif
           if (shellt(i).eq.3) then 
               if (dfree.ne.0) then
                  if (gnreal(r,2,.true.)) then
                     exx(mm) = r(1)
                     c4(mmdf)= r(2)
                  endif
               else
                  if (ig03.gt.0) then
                     read(iun2,1040,err=1000) exx(mm), c4(mmdf)
                  else
                     read(iun2,1030,err=1000) exx(mm), c4(mmdf)
                  endif
               endif
               mmdf = mmdf + 1
           endif
           if (shellt(i).eq.4) then 
               if (dfree.ne.0) then
                  if (gnreal(r,2,.true.)) then
                     exx(mm) = r(1)
                     c5(mmdf)= r(2)
                  endif
               else
                  if (ig03.gt.0) then
                     read(iun2,1040,err=1000) exx(mm), c5(mmdf)
                  else
                     read(iun2,1030,err=1000) exx(mm), c5(mmdf)
                  endif
               endif
               mmdf = mmdf + 1
           endif
   20      mm = mm + 1
      goto 10

1000  call inferr('error while reading Basis set!',1)
      istats = 0
      return
      end
