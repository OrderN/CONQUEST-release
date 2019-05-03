      subroutine files(oempty,idocub)
      implicit double precision (a-h,o-z)
      character keywrd*320, keyori*320
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      logical oempty
      character fniun*256
c      character temps*80
      common /fnunit/ fniun
      common /keywrd/ keywrd,keyori

c
c----- process FILE or CUBE directive --------------------
c
      idocub = 0
      i = index(keywrd,'FILE')
      if (i.eq.0) then
         i = index(keywrd,'GAUCUB')
         idocub = 1
      endif
      if (i.ne.0) then
         ind = i + 4
         i   = index(keywrd(ind:320),'=')
         do j=i+ind,320
            if (ichar(keywrd(j:j)).ne.32) goto 2170
         end do
2170     ind1 = j
         do j=ind1,320
            if (ichar(keywrd(j:j)).eq.32) goto 2190
         end do
2190     ind2 = j-1
         open(unit=iun2,form='formatted',file=keyori(ind1:ind2),
     &        status='old',err=2200)
         fniun = keyori(ind1:ind2)
         do j=ind1,ind2
            keyori(j:j) = ' '
            keywrd(j:j) = ' '
         end do
      else
         call inferr('Keyword FILE not supplied !!',1)
         oempty = .true.
      endif

      return

2200  write(iun3,*) 
     &   'File '//keyori(ind1:ind2)//' does not exist !!!!!!!'
c      call inferr(temps,1)
      oempty = .true.
      return
      end

      subroutine filmap(mapit)
      implicit double precision (a-h,o-z)
      character keywrd*320, keyori*320
      common /keywrd/ keywrd,keyori
      character*80 mapfil
      character*80 grdfil
      common /maphlp/ mapfil,grdfil

c
c----- process MAPFIL directive --------------------
c
      mapit = 0
      i = index(keywrd,'MAPFIL')
      if (i.ne.0) then
         mapit = 1
         ind = i + 4
         i   = index(keywrd(ind:320),'=')
         do j=i+ind,320
            if (ichar(keywrd(j:j)).ne.32) goto 2170
         end do
2170     ind1 = j
         do j=ind1,320
            if (ichar(keywrd(j:j)).eq.32) goto 2190
         end do
2190     ind2 = j-1
         mapfil = keyori(ind1:ind2)
      endif

      return
      end

