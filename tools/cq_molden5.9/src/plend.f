      subroutine plend(title,notitle)
      implicit double precision (a-h,o-z)
      character title*80, keywrd*320, keyhlp*80, keyori*320
      character esc, eot
      logical notitle
      common /keywrd/ keywrd,keyori
      common /rdwr/   iun1,iun2,iun3,iun4,iun5


      esc = char(27)
      eot = char(4)

      keyhlp(1:80)  = title(1:80)
      keyori(1:320) = keywrd(1:320)
      n1 = 1
      n2 = 1
      do i=1,320
         if (i.lt.80) then
           if (title(i:i).eq.'('.or.title(i:i).eq.')') then
             keyhlp(i+n1:80)       = title(i:80-n1)
             keyhlp(i+n1-1:i+n1-1) = char(92)
             n1 = n1 + 1
           endif
         endif
         if (keywrd(i:i).eq.'('.or.keywrd(i:i).eq.')') then
           keyori(i+n2:320)      = keywrd(i:320-n2)
           keyori(i+n2-1:i+n2-1) = char(92)
           n2 = n2 + 1
         endif
      end do

      if (.not.notitle) then
         write(iun4,'(''titleandlogo {'')')
         write(iun4,'(''labelcol'')')
         write(iun4,'(''2100 1700 m'')')
         write(iun4,'(''('',a,'') show'')')keyhlp(1:40)
         write(iun4,'(''2100 1640 m'')')
         write(iun4,'(''('',a,'') show'')')keyhlp(41:80)
         write(iun4,'(''2100 1580 m'')')
         write(iun4,'(''('',a,'') show'')')keyori(1:40)
         write(iun4,'(''2100 1520 m'')')
         write(iun4,'(''('',a,'') show'')')keyori(41:80)
         write(iun4,'(''2100 1460 m'')')
         write(iun4,'(''('',a,'') show'')')keywrd(81:120)
         write(iun4,'(''2100 1400 m'')')
         write(iun4,'(''('',a,'') show'')')keywrd(121:160)
         write(iun4,'(''/print { 0 0 moveto '',
     &       ''4 0 (MOLDEN) ashow } def'')')
         write(iun4,'(''2200 1850 translate'')')
         write(iun4,'(''4.0 4.0 scale'')')
         write(iun4,'(''.95 -.05 0'')')
         write(iun4,
     &       '(''{setgray print -1.5 +1.5 translate } for'')')
         write(iun4,'(''1 setgray print'')')
         write(iun4,'(''} if'')')
      endif
      call plpend

      return
      end

      subroutine plpend
      implicit double precision (a-h,o-z)
      common /rdwr/   iun1,iun2,iun3,iun4,iun5

      write(iun4,'(''/#copies    1 def'')')
      write(iun4,'(''showpage'')')
      write(iun4,'(''39.00 753.00 translate'')')
      write(iun4,'(''-90.00 rotate'')')
      write(iun4,'(''1 setlinewidth'')')
      write(iun4,'(''2 setlinecap'')')
      write(iun4,'(''0.240000 0.240000 scale'')')
      write(iun4,'(''n'')')

      return
      end
