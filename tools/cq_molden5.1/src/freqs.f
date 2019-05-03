      subroutine getfd(istat,coo)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      character*137 lstr,str
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /athlp/ iatoms, mxnat
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      character*5 frsym
      common /frstr/ frsym(maxfrq)
      integer getlin
      logical g94,dosym
      dimension coo(3,*)

      istat = 1

      g94 = .false.
      dosym = .true.
      ihasi = 2
      ivibs = 0
      idirct = 1
      nframe = 5
      iframe = 0

      rewind iun2

      call iatnox(natoms)

      do i=1,natoms
         do j=1,3
             fcoo(j,i) = coo(j,i)
         end do
      end do

c     Read in Gaussian Frequencies

      nvibs = 3*natoms - 6
      if (natoms.eq.1) nvibs = 0
      if (natoms.eq.2) nvibs = 1

13    call search(lstr,'Frequencies --',istat)
      if (istat.eq.0) goto 100
      if (ivibs.gt.nvibs) goto 100
      g94 = (index(lstr,'Frequencies --').eq.2)

      if (dosym) then
         backspace iun2
         backspace iun2
         idum = getlin(0)
         if (g94) then
            nv = 3
         else
            nv = 5
         endif
         do i=1,nv
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.1) then
               ns = nstr
               if (ns.gt.3) ns = 3
               frsym(ivibs+i) = '   '
               frsym(ivibs+i) = str(1:ns)
            else
            endif
         end do
         call readel(lstr,1)
      endif

      if (g94) then
         indfr = index(lstr,'Frequencies --')
         lstr = lstr(indfr+14:)
         read(lstr,'(3(f10.4,13x))',err=100,end=100) 
     &          (freq(i),i=ivibs+1,ivibs+3)
         call readel(lstr,3)
         lstr = lstr(indfr+14:)
         read(lstr,'(3(f10.4,13x))',err=100,end=100) 
     &          (frint(i),i=ivibs+1,ivibs+3)
         call readel(lstr,1)
         if (index(lstr,'Raman').eq.0) then
            ihasi = 1
         else
            lstr = lstr(indfr+14:)
            read(lstr,'(3(f10.4,13x))',err=100,end=100) 
     &          (ramint(i),i=ivibs+1,ivibs+3)
         endif
         ivibs = ivibs + 3
      else
         indfr = index(lstr,'Frequencies --')
         lstr = lstr(indfr+16:)
         read(lstr,'(5f10.4)',err=100,end=100) 
     &          (freq(i),i=ivibs+1,ivibs+5)
         call readel(lstr,3)
         lstr = lstr(indfr+16:)
         read(lstr,'(5f10.4)',err=100,end=100) 
     &          (frint(i),i=ivibs+1,ivibs+5)
         call readel(lstr,1)
         if (index(lstr,'Raman').eq.0) then
            ihasi = 1
         else
            lstr = lstr(indfr+16:)
            read(lstr,'(5f10.4)',err=100,end=100) 
     &          (ramint(i),i=ivibs+1,ivibs+5)
         endif
         ivibs = ivibs + 5
      endif
      goto 13

100   if (ivibs.eq.0) istat = 0
200   if (freq(ivibs).eq.0.0d0.and.ivibs.ne.0) then
         ivibs = ivibs - 1
         goto 200
      endif
      nfreq = ivibs
      ihasi = 2
      
      if (dosym) then
         ihasi = -ihasi
         do i=1,ivibs
            call parsfn(frsym(i),3,12)
         end do
      endif

      return
      end 

      subroutine ncoord(idebug,ifreq,istat)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      character*137 lstr
c      character*2 gstr
      character*5 tstr
      real freqt
      dimension freqt(5),frtt(3,3)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      logical g94,g09

c     Get ifreq 'th Norm. Cordinates from Gaussian Output

      istat = 1
      ifroz = 0
      g94 = .false.
      g09 = .false.

      rewind iun2
      ilegnd = 3
      call iatnox(iatoms)
      nvibs = nfreq

      call search(lstr,'atoms frozen in the vibrational analysis',istat)
      if (istat.ne.0) then
         read(lstr,'(i5)') nfroz
         ifroz = 1
      else
         rewind iun2
      endif

      call search(lstr,'Frequencies -- ',istat)
      if (istat.eq.0) goto 100
      g94 = (index(lstr,'Frequencies --').eq.2)

c Find out length of legend

      if (g94) then
         do i=1,8
             read(iun2,'(a)',end=100,err=100) lstr
             g09 = (index(lstr,'Atom  AN      X').eq.3)
             if (index(lstr,'Atom AN      X').ne.0.or.
     &           index(lstr,'Atom  AN      X').ne.0) goto 15
             ilegnd = ilegnd + 1
         end do
15       continue
      else
          ilegnd = 8
      endif

      jatoms = iatoms
      if (ifroz.eq.1.and.g09) jatoms = iatoms - nfroz
      
      call scback(lstr,'and normal coordinates:',istat)
      call readel(lstr,1)
      if (istat.eq.0) goto 100
      if (g94) then
         ioff = ifreq - ((ifreq-1)/3)*3
      else
         ioff = ifreq - ((ifreq-1)/5)*5
      endif


      ivibs = 0
10    read(iun2,'(a)',end=100,err=100) lstr
      if (ivibs.ge.nvibs) goto 20
      if (index(lstr,'and normal coordinates:').ne.0) goto 20
      tstr = '     '
      write(tstr,'(1a,i3,1a)') ' ',ifreq,' '
c      tstr = ' '//gstr(ifreq)//' '

      if (index(lstr,tstr).ne.0) then
         call readel(lstr,ilegnd)
         if (g94) then
            do i=1,jatoms
               if (g09) then
c                  read(iun2,'(10x,3(2x,3f7.2))',err=100,end=100) 
c     &             ((frtt(j,k),j=1,3),k=1,3)
                  read(iun2,'(i6,4x,3(2x,3f7.2))',err=100,end=100) 
     &             iattmp,((frtt(j,k),j=1,3),k=1,3)
                  do j=1,3
                      a(j,iattmp) = frtt(j,ioff)
                  end do
               else
                  read(iun2,'(8x,3(2x,3f7.2))',err=100,end=100) 
     &             ((frtt(j,k),j=1,3),k=1,3)
                  do j=1,3
                      a(j,i) = frtt(j,ioff)
                  end do
               endif
            end do
         else
            do i=1,iatoms
               do j=1,3
                   read(iun2,'(23x,5f10.5)',err=100,end=100) 
     &                  (freqt(k),k=1,5)
                   a(j,i) = freqt(ioff)
               end do
            end do
         endif
      else
         if (g94) then
            call readel(lstr,jatoms+ilegnd)
         else
            call readel(lstr,iatoms*3+ilegnd)
         endif
      endif
      if (g94) then
         ivibs = ivibs + 3
      else
         ivibs = ivibs + 5
      endif
      goto 10

20    if (idebug.eq.1) call prtfr(ifreq)
      return

100   istat = 0
      call inferr('Error reading Norm. Coords. !',0)
      return
      end

      subroutine ggetfd(istat,coo)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      character*137 line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /athlp/ iatoms, mxnat
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      dimension coo(3,*)

cd     write(iun3,'(a)')'enter subroutine ggetfr'
      istat = 1

      ivibs = 0
      idirct = 1
      nframe = 5
      iframe = 0
      rewind iun2

      call iatnox(natoms)

      do i=1,natoms
        do j=1,3
           fcoo(j,i) = coo(j,i)
        end do
      end do

c     Read in Gamess Frequencies

      nvibs = 3*natoms - 6
      if (natoms.eq.1) nvibs = 0
      if (natoms.eq.2) nvibs = 1
      call search(line,'cartesians to normal mode',istat)
      if (istat.eq.0) goto 10

c...  This is a Hessian run
      call readel(line,6)
50    call readel(line,1)
      if (ivibs.lt.nvibs+6) then
         read(line,'(15x,8f11.2)',end=100,err=100)
     &        (freq(i),i=ivibs+1,ivibs+8)
         ivibs = ivibs + 8
         call readel(line,3*iatoms+7)
         goto 50
      endif
      n = ivibs
      ivibs = 0
      do i=1,n
         if (freq(i).ne.0.0d0) then
            ivibs=ivibs+1
            freq(ivibs)=freq(i)
         endif
      enddo
      goto 200

10    continue
c...  This is a Force run
      call search(line,'frequencies ----',istat)
      if (istat.eq.0) goto 100
      if (ivibs.gt.nvibs) goto 100
      line = line(index(line,'frequencies ----')+22:)
      read(line,'(9f10.4)',err=100,end=100) 
     &     (freq(i),i=ivibs+1,ivibs+9)
      ivibs = ivibs + 9
      goto 10

100   if (ivibs.eq.0) istat = 0
200   if (freq(ivibs).eq.0.0d0.and.ivibs.ne.0) then
         ivibs = ivibs - 1
         goto 200
      endif
      nfreq = ivibs
      ihasi = 0
cd     write(iun3,'(a)')'leave subroutine ggetfr'
      return
      end 

      subroutine ncoorg(idebug,ifreq,istat)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      character*137 lstr
      character*2 gstr
      character*4 tstr
      real freqt
      dimension freqt(9)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi

c     Get ifreq 'th Norm. Cordinates from Gamess Output
cd     write(iun3,'(a)')'enter subroutine ncoorg'

      istat = 1

      rewind iun2
      call iatnox(iatoms)
      nvibs = nfreq

      call search(lstr,'cartesians to normal mode',istat)
      if (istat.eq.0) goto 30
c
c...  This is a Hessian run
c
      if (freq(ifreq).lt.0.0d0) then
         ifreq1 = ifreq
      else
         ifreq1 = ifreq+6
      endif
      ioff = ifreq1 - ((ifreq1-1)/8)*8
      ivibs = 0
      tstr = ' '//gstr(ifreq1)//' '
      call readel(lstr,3)
40    read(iun2,'(a)',end=100,err=100) lstr
      if (ivibs.lt.nvibs+6) then
         if (index(lstr,tstr).ne.0) then
            call readel(lstr,5)
            do i=1,iatoms
               do j=1,3
                  read(iun2,'(15x,8f11.6)',err=100,end=100)
     &                 (freqt(k),k=1,8)
                  a(j,i) = freqt(ioff)
               end do
            end do
            call readel(lstr,2)
         else
            call readel(lstr,iatoms*3+7)
         endif
         ivibs = ivibs+8
         goto 40
      endif
      goto 20

30    continue
c
c...  This is a Force run
c
      ioff = ifreq - ((ifreq-1)/9)*9
      call search(lstr,'and normalised normal coordinates',istat)
      if (istat.eq.0) goto 100
      call readel(lstr,3)

      ivibs = 0
10    read(iun2,'(a)',end=100,err=100) lstr
      if (ivibs.gt.nvibs) goto 20
      if (index(lstr,'===========================').ne.0) goto 20
      tstr = ' '//gstr(ifreq)//' '

      if (index(lstr,tstr).ne.0) then
          call readel(lstr,7)
          do i=1,iatoms
             do j=1,3
                 read(iun2,'(23x,9f10.5)',err=100,end=100) 
     &                (freqt(k),k=1,9)
                 a(j,i) = freqt(ioff)
             end do
          end do
          call readel(lstr,2)
      else
          call readel(lstr,iatoms*3+9)
      endif
      ivibs = ivibs + 9
      goto 10

20    if (idebug.eq.1) call prtfr(ifreq)
cd     write(iun3,'(a)')'leave subroutine ncoorg'
      return

100   istat = 0
      call inferr('Error reading Norm. Coords. !',0)
cd     write(iun3,'(a)')'leave subroutine ncoorg'
      return
      end

      subroutine getfrd(istat,coo,ianz)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      character*137 lstr, str
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /athlp/ iatoms, mxnat
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      common /frmop/ imop6
      character*5 frsym
      common /frstr/ frsym(maxfrq)
      integer getlin
      logical dosym
      character*2 iel,catomt,catom,tolowf
      dimension iel(100),xyz(3)
      dimension coo(3,*),ianz(*)
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

      istat = 1

      imop6 = 0
      ivibs = 0
      idirct = 1
      nframe = 5
      iframe = 0
      ihasi = 0
      toau = 0.52917706d0
      dosym = .true.

      rewind iun2

      call search(lstr,'ORIENTATION OF MOLECULE IN FORCE',istat)
      if (istat.eq.0) goto 1000

      icoord = 0
      iatoms = 0
      do while (getlin(0).eq.1)
         if (nxtwrd(str,nstr,itype,rtype).ne.2) goto 100
         ktype = nxtwrd(str,nstr,itype,rtype)
         if (ktype.eq.1.and.(nstr.eq.1.or.nstr.eq.2)) then
             if (nstr.eq.1) then
                 catomt(1:1) = str(1:1)
                 catomt(2:2) = ' '
             else
                 catomt = str(1:2)
             endif
             catom = tolowf(catomt)
         else
             if (ktype.eq.2) then
                 imop6 = 1
                 ianzt = itype
             else
                 goto 100
             endif
         endif
         do j=1,3
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.3) then
                xyz(j) = rtype
            else
                goto 100
            endif
         end do
         if (icoord.eq.0) icoord = 1
         iatoms = iatoms + 1
         if (imop6.eq.1) then
            ianz(iatoms) = ianzt
         else
            do i=1,100
               if ( catom .eq. iel(i) ) ianz(iatoms) = i - 1
            end do
         endif
         do j=1,3
            coo(j,iatoms) = xyz(j) / toau
         end do
         goto 150
100      if (icoord.eq.1) goto 200
150      continue
      end do

200   continue


      do i=1,iatoms
         do j=1,3
             fcoo(j,i) = coo(j,i)
         end do
      end do


c     Read in Mopac Frequencies

      nvibs = 3*iatoms - 6
      if (iatoms.eq.1) nvibs = 0
      if (iatoms.eq.2) nvibs = 1

c      call search(lstr,'NORMAL COORDINATE ANALYSIS',istat)
      call search(lstr,'MASS-WEIGHTED COORDINATE ANALYSIS',istat)
13    call search(lstr,'Root No.',istat)
      if (istat.eq.0) goto 1000
         if (ivibs.eq.nvibs) goto 1000
         if (imop6.eq.1) then
            call readel(lstr,1)
         else
            if (dosym) then
               call readel(lstr,2)
               read(lstr,'(12x,10(4x,a4))',err=1000,end=1000)
     &               (frsym(ivibs+i),i=1,8)
               call readel(lstr,1)
            else
               call readel(lstr,3)
            endif
         endif
         idum = getlin(0)
         do while (nxtwrd(str,nstr,itype,rtype).eq.3) 
            ivibs = ivibs + 1
            freq(ivibs) = rtype
         end do
      goto 13

1000  if (ivibs.eq.0) istat = 0
      nfreq = ivibs

c find transition dipole moments (T-DIPOLE)

      rewind iun2
      do i=1,nfreq
         frint(i) = -1.0d0
      end do
      call search(lstr,'DESCRIPTION OF VIBRATIONS',istat)
      if (istat.eq.1) then
         do i=1,nfreq
            do while (getlin(0).eq.1)
               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.1.and.nstr.eq.9) then
                  if (str.eq.'VIBRATION') then
                     ktype = nxtwrd(str,nstr,itype,rtype)
                     if (ktype.eq.2) then
                        ivib = itype
                        call readel(lstr,1)
                        idum = getlin(0)
                        if (idum.eq.1) then
                           ktype = nxtwrd(str,nstr,itype,rtype)
                           if (ktype.eq.1.and.nstr.eq.8) then
                              if (str.eq.'T-DIPOLE') then
                                 ktype = nxtwrd(str,nstr,itype,rtype)
                                 if (ktype.eq.3) then
                                  if (ivib.gt.0.and.ivib.le.maxfrq) then
                                    frint(ivib) = rtype
                                  endif
                                 else
                                    goto 3000
                                 endif
                              else
                                 goto 3000
                              endif
                           else
                              goto 3000
                           endif
                        else
                           goto 3000
                        endif
                     endif
                     goto 2000
                  else
                     goto 3000
                  endif
               endif
            end do

c found frequency T-DIPOLE

2000        continue
         end do
      else
         goto 3000
      endif
      ihasi = 1


      if (dosym) then
         ihasi = -ihasi
         do i=1,ivibs
            call parsfn(frsym(i),3,12)
         end do
      endif

      return

3000  ihasi = 0
      ihasi = 1
      return
      end 

      subroutine mcoord(idebug,ifreq,istat)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      character*137 lstr,line
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /frmop/ imop6
      common /curlin/ line
      integer getlin

c     Get ifreq 'th Norm. Cordinates from Mopac Output

      istat = 1

      rewind iun2
      toau = 0.52917706d0
      if (imop6.eq.1) then
          nvlin = 6
      else
          nvlin = 8
      endif
      call iatnox(iatoms)
      nvibs = nfreq
      ioff = ifreq - ((ifreq-1)/nvlin)*nvlin
c      call search(lstr,'NORMAL COORDINATE ANALYSIS',istat)
c
c      if you choose for the above the algorithm still comes up
c      with the MASS-WEIGHTED set, have to look at this some time
c
      call search(lstr,'MASS-WEIGHTED COORDINATE ANALYSIS',istat)
      if (istat.eq.0) goto 100

      ivibs = 0
10    call search(lstr,'Root No.',istat)
      if (istat.eq.0.or.ivibs.ge.nvibs) goto 20
      backspace iun2
      n = getlin(0)
      do while (.true.)
         n = nxtwrd(lstr,nstr,itype,rtype)
         if (n.eq.0) then
            goto 10
         elseif (n.eq.2.and.itype.eq.ifreq) then
            goto 15
         endif
      end do
15    continue
      if (imop6.eq.1) then
         call readel(lstr,3)
      else
         call readel(lstr,5)
      endif

      do i=1,iatoms
         do j=1,3
             idum = getlin(0)
             if (imop6.eq.1.and.line(1:1).eq.'1') idum = getlin(0)
             do n=1,ioff+1
                ktype = nxtwrd(lstr,nstr,itype,rtype)
             end do
             if (ktype.eq.3) then
                a(j,i) = rtype / toau
             endif
         end do
      end do
      ivibs = ivibs + nvlin
      goto 10

20    if (idebug.eq.1) call prtfr(ifreq)
      return

100   istat = 0
      call inferr('Error reading Norm. Coords. !',0)
      return
      end

      subroutine getfrad(istat,ianz)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      parameter (mxel=100)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /athlp/ iatoms, mxnat
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      character*2 elemnt
      common /elem/elemnt(mxel)
      character*137 line,str
      common /curlin/ line
      integer getlin
      logical gnreal
      character*2 tocapf
      dimension r(3)
      dimension ianz(*)


      istat =  1

      ivibs = 0
      idirct = 1
      nframe = 5
      iframe = 0

      call rewmf

      call srchmf(line,'[FREQ]',istats)
      if (istats.eq.0) goto 100
 
      do while (getlin(0).eq.1)

          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.3) then
             ivibs = ivibs + 1
             freq(ivibs) = rtype
          elseif (ktype.eq.1) then
             if (str(1:1).eq.'[') goto 10
          endif

      end do

10    nfreq = ivibs

      call rewmf
      call srchmf(line,'[FR-COORD]',istats)
      if (istats.eq.0) goto 100

      iat = 0
      do while (getlin(0).eq.1)

          ktype = nxtwrd(str,nstr,itype,rtype)

          if (ktype.eq.1) then
              if (str(1:1).eq.'[') goto 20
              if (nstr.eq.1) then
                  str(2:2) = str(1:1)
                  str(1:1) = ' '
                  nstr = 2
              endif

              ifound = 0
              if (nstr.eq.2) then
                 do i=1,mxel
                    if (tocapf(str(1:2)).eq.tocapf(elemnt(i))) 
     &                  ifound = i
                 end do
                 if (ifound.ne.0) then
                    iat = iat + 1
                    ianz(iat) = ifound
                    if (gnreal(r,3,.false.)) then
                       do i=1,3
                          fcoo(i,iat) = r(i)
                       end do
                    endif
                 endif
              endif

          endif

      end do

20    iatoms = iat

      call resfr

      return

100   istat = 0
      call inferr('Error reading Norm. Coords. !',0)
      return
      end

      subroutine getint(istat)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      character*137 line,str
      common /curlin/ line
      integer getlin

      istat =  1

      i = 0
      ihsv = 0
      ihsram = 0
      ihasi = 0

      call rewmf

      call srchmf(line,'[INT]',istats)
      if (istats.eq.0) then
          istat = 0
          return
      endif
 
      do while (getlin(0).eq.1)

          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.3) then
             i = i + 1
             frint(i) = rtype
             ihsv = 1
          elseif (ktype.eq.1) then
             if (str(1:1).eq.'[') goto 10
          endif

          ktype = nxtwrd(str,nstr,itype,rtype)
          if (ktype.eq.3) then
             ramint(i) = rtype
             ihsram = 1
          elseif (ktype.eq.1) then
             if (str(1:1).eq.'[') goto 10
          elseif (ktype.eq.0) then
          endif

      end do

10    ihasi = ihsv + ihsram

      return
      end

      subroutine acoord(idebug,ifreq,istat)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      character*137 line,str
      common /curlin/ line
      integer getlin
      logical gnreal
      dimension r(3)

      istat = 1
      call rewmf
      call iatnox(iatoms)
      call srchmf(line,'[FR-NORM',istats)
      if (istats.eq.0) goto 100

      do while (getlin(0).eq.1)

          ktype = nxtwrd(str,nstr,itype,rtype)

          if (ktype.eq.1) then
             if (nstr.ge.3) then
                 if (icdex(str,'vib').ne.0) then
                    ktype = nxtwrd(str,nstr,itype,rtype)
                    if (ktype.eq.2) then
                       if (itype.eq.ifreq) goto 50
                    else
                       goto 100
                    endif
                 endif
             endif
          endif

      end do

      goto 100

50    do i=1,iatoms
         
          if (gnreal(r,3,.true.)) then
             do j=1,3
                a(j,i) = r(j)
             end do
          else
             goto 60
          endif

      end do

60    if (idebug.eq.1) call prtfr(ifreq)

      return

100   istat = 0
      call inferr('Error reading Norm. Coords. !',0)
      return
      end

      subroutine frqstr(ifrq,movfil,movlen)
      implicit double precision (a-h,o-z)
      character*11 movfil
      character*5 hstr,tstr

      tstr = hstr(ifrq)
      i = index(tstr,')')
      movlen = 4+i-2+4
      movfil = 'freq'//tstr(2:i-1)//'.xyz'

      return
      end

      subroutine nxtpnd(iopt,fancy,atcol,dolabs,backb,wpnt,
     &                  coo)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      common /athlp/ iatoms, mxnat
      common /pnthlp/ ipoints,ipnt
      logical wpnt
      integer dolabs,fancy,atcol,backb
      character*11 movfil
      character*5 hstr,tstr
      dimension coo(3,*)

      call iatnox(natoms)

      if (wpnt) then
         call frqstr(ifrq,movfil,movlen)
         if (iopt.eq.5) movfil(movlen-2:movlen) = 'wrl'
         if (iopt.eq.6) movfil = 'molden.ogl'
      endif

      iframe = iframe + idirct
      if (iframe.ge.nframe.or.iframe.le.-nframe) then
         idirct = -1*idirct
         iloop = iloop + 1
      endif

      do i=1,natoms
        do j=1,3
           coo(j,i) = fcoo(j,i) + iframe*frmul*a(j,i)
        end do
      end do

      if (wpnt.and.iloop.lt.3) then

        if (iloop.eq.1.and.iframe.eq.nframe) then
            call wrpnt(movfil,movlen,iopt,1,nframe*2,
     &                 fancy,atcol,dolabs,backb)
        elseif (iloop.eq.2.and.iframe.eq.nframe-1) then
            ipts = ipoints
            ipoints = nframe*4
            call wrpnt(movfil,movlen,iopt,4,0,
     &                 fancy,atcol,dolabs,backb)
            ipoints = ipts
            call inferr('wrote file: '//movfil,0)
            if (iopt.eq.6) call viewer
        else
            call wrpnt(movfil,movlen,iopt,3,0,
     &                 fancy,atcol,dolabs,backb)
        endif

      endif

      return
      end

      subroutine scalfd(iopt,fancy,atcol,dolabs,backb,wpnt,
     &                  t,coo,scal,scali,smag)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      logical wpnt
      integer dolabs,fancy,atcol,backb
      dimension t(3),coo(3,*)

      idirct = 1
      iframe = nframe - 1
      call nxtpnt(iopt,fancy,atcol,dolabs,backb,wpnt)

      call iatnox(natoms)

      scali = 0.0d0
      do i=1,natoms
         dist = 0.0d0
         do j=1,3
             dt = coo(j,i)-t(j)
             dist = dist + dt*dt
         end do
         if (dist.gt.scali) scali = dist
      end do
      scali = dsqrt(scali)
      scal = scali * 2.4d0 * smag

      return
      end

      subroutine resfd(coo)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      dimension coo(3,*)

      call iatnox(natoms)

      do i=1,natoms
        do j=1,3
           coo(j,i) = fcoo(j,i)
        end do
      end do

      return
      end

      subroutine tofcod(coo)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      dimension coo(3,*)

      call iatnox(natoms)

      do i=1,natoms
        do j=1,3
           fcoo(j,i) = coo(j,i)
        end do
      end do

      return
      end

      subroutine prtfr(ifreq)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      call iatnox(iatoms)

      write(iun3,*) ' '
      write(iun3,*) 'Coordinates used in frequency calculation (au)'
      write(iun3,*) '============='
      do i=1,iatoms
         write(iun3,*) (fcoo(j,i),j=1,3)
      end do

      write(iun3,*) ' '
      write(iun3,*) 'Norm. Mode ',ifreq,' (au)'
      write(iun3,*) '============='

      do i=1,iatoms
         write(iun3,*) (a(j,i),j=1,3)
      end do

      return
      end

      logical function negfrq(ifreq)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi

      negfrq = .false.

      if (freq(ifreq).lt.0.0e0) negfrq = .true.

      return
      end

      subroutine iatnod(iatnox,ianz)
      implicit double precision (a-h,o-z)
      common /athlp/ iatoms, mxnat
      dimension ianz(*)

      iatnox = 0
      do i=1,iatoms
         if (ianz(i).ne.99) iatnox = iatnox + 1
      end do

      return
      end

      subroutine calspc
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      real freq,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      logical opfil

      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      if (opfil(46,'fr.out',6,1,0,0)) then
         do i=1,nfreq
            write(46,'(f10.4,2x,f10.4)') freq(i),frint(i)
         end do
         close(46)
      endif

      return
      end

      subroutine gttrns(istats)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /athlp/ iatoms, mxnat
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifrq,normc
      character*5 frsym
      common /frstr/ frsym(maxfrq)
      character*137 line,str
      common /curlin/ line
      integer getlin
      double precision reada
      logical g94,dosym

      istats = 1

      ihasi = -2
      nfreq = 0
      do while(.true.)
         call search(line,"Excited State ",istat)
         if (istat.ne.0) then
            if (index(line,":").ne.0) then
                nfreq = nfreq + 1
                do i=1,4
                   ktype = nxtwrd(str,nstr,itype,rtype)
                end do
                if (ktype.eq.1) then
                    ind = index(str,"Singlet")
                    if (ind.eq.1) then
                       str = 'S'//str(8:)
                    endif
                    ind = index(str,"Triplet")
                    if (ind.eq.1) then
                       str = 'T'//str(8:)
                    endif
                    n = len(str)
                    frsym(nfreq) = str(1:n)
                endif
                do i=1,3
                   ktype = nxtwrd(str,nstr,itype,rtype)
                end do
                if (ktype.eq.3) then
                    freq(nfreq) = rtype
                endif
                do i=1,2
                   ktype = nxtwrd(str,nstr,itype,rtype)
                end do
                if (ktype.eq.1) then
                   frint(nfreq) = reada(str,3,nstr)
                endif
            endif
         else
            goto 10
         endif
      end do

10    continue

      do i=1,nfreq
         n = len(frsym(i))
         call parsfn(frsym(i),n,12)
      end do

      return
      end

      subroutine nrmi(ifrq)
      implicit double precision (a-h,o-z)
      logical wmovie,dconpt

      wmovie = .false.
      dconpt = .false.

      call setnrm(ifrq,0,wmovie,dconpt)

      return
      end

      subroutine setnrm(ifrq,iwopt,wmovie,dconpt)
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      character*137 line,str
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi
      common /athlp/ iatoms, mxnat
      common /frwrk/ frthr,frmul,idirct,iframe,nframe,iloop,ifreq,normc
      character*5 frsym
      common /frstr/ frsym(maxfrq)
      logical ctoz,molpot,elpot,chpot
      common /choic/  iftyp,isbin,ctoz,molpot,elpot,chpot
      integer dolabs,fancy,persp,shade,atcol,fyesno,backb
      common /displ/ fancy,shade,atcol,dolabs,persp,irtcel,
     &               ifd,fyesno,backb,logo
      common /getpnt/ irtype,ipdbon,ipdbgro,ifav,ioxyz,
     &                iconv,ircus,dozme
      logical dconpt,negfrq,wmovie,dozme

      idebug = 0
      normc = 1
      ifreq = ifrq
      dconpt = negfrq(ifreq)
      if (irtype.eq.4) then
           iloop = 0
           call resfr
           call doconn

           if (iftyp.eq.1) then
              call mcoord(idebug,ifreq,istat)
           elseif (iftyp.eq.2.or.iftyp.eq.3) then
              if (iftyp.eq.3) then
                  call ucoorg(idebug,ifreq,istat)
              else
                  call ncoorg(idebug,ifreq,istat)
              endif
           elseif (iftyp.eq.4) then
              call curs(1)
              rewind iun2
              call search(line,'Frequencies --',istat)
              call scback(line,'orientation:',istat)
              if (istat.eq.1) then
                   backspace iun2
                   call rdcor(idebug,istat)
                   if (istat.eq.1) then
                       call xyzcoo(1,0,0)
                       call tofcoo
                       call doconn
                   endif
              endif
              call ncoord(idebug,ifreq,istat)
              call curs(0)
           elseif (iftyp.eq.5) then
              call acoord(idebug,ifreq,istat)
           elseif (iftyp.eq.7) then
              call cpmdcoorg(idebug,ifreq,istat)
           elseif (iftyp.eq.8) then
              call qcoord(idebug,ifreq,istat)
           elseif (iftyp.eq.9) then
              call ocoord(idebug,ifreq,istat)
           elseif (iftyp.eq.15) then
              call nwcord(idebug,ifreq,istat)
           endif
           call scalfr(iwopt,fancy,atcol,dolabs,backb,wmovie)
      endif

      return
      end

      integer function nfrqs()
      implicit double precision (a-h,o-z)
      parameter (maxfat=1000)
      parameter (maxfrq=maxfat*3)
      real freq,frmul,a
      common /freq/ freq(maxfrq),a(3,maxfat),fcoo(3,maxfat),
     &              frint(maxfrq),ramint(maxfrq),nfreq,ihasi

      nfrqs = nfreq

      return
      end
